#!/usr/bin/env python3
"""
UV-Index Schweiz — HSURF-Fetcher (Stufe 1).

Holt die Geländehöhe (HSURF) aus dem STATISCHEN horizontal_constants-File der
ICON-CH1-Collection und rastert sie auf ein regelmässiges Lon/Lat-Gitter über
der Schweiz. Ergebnis sind zwei Dateien, die ins UV-Projekt committet werden:

    hsurf_ch.bin    int16 Little-Endian, Höhe in Metern, zeilenweise
    hsurf_ch.json   Metadaten (BBox, Dimensionen, Zeilenreihenfolge, nodata)

Warum ein externer Worker: Infomaniak Shared Hosting hat kein eccodes und
exec() ist deaktiviert → GRIB2 lässt sich dort nicht dekodieren. Muster wie
weitsicht-forecast-fetcher / wetteralarm-hail-fetcher.

Warum nur EINMAL: Gelände ist statisch. Dieser Worker läuft ausschliesslich
per workflow_dispatch, bewusst ohne `schedule:`. Nach dem Commit der Ausgabe
ins UV-Projekt besteht KEINE Laufzeitabhängigkeit zu diesem Repo mehr.

Warum ICON-CH1 statt CH2: CH1 hat 1 km statt 2.1 km Maschenweite. Der kürzere
Vorhersagehorizont (33 h vs. 120 h) spielt hier keine Rolle — Gelände hat
keinen Vorhersagehorizont.

Env (alle optional):
    GRID_KM         Rasterweite in km, default 1.0
    BBOX            "west,south,east,north", default 5.90,45.75,10.55,47.85
    OUT_PREFIX      Dateiname ohne Endung, default hsurf_ch

Exit-Codes: 0 ok · 1 Fehler
"""
from __future__ import annotations

import json
import logging
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("fetch-hsurf")

HERE = Path(__file__).resolve().parent

STAC_BASE = (
    "https://data.geo.admin.ch/api/stac/v1/collections/"
    "ch.meteoschweiz.ogd-forecasting-icon-ch1"
)
HORIZONTAL_CONSTANTS = "horizontal_constants_icon-ch1-eps.grib2"

# Schweiz mit Rand — deckungsgleich mit den Bounds im UV-Projekt.
DEFAULT_BBOX = (5.90, 45.75, 10.55, 47.85)   # west, south, east, north
NODATA = -32768

# Plausibilitätsgrenzen für die Schweiz (Lago Maggiore 193 m … Dufourspitze
# 4634 m). Auf 1 km Maschenweite ist das Gelände geglättet, die Extremwerte
# werden also nicht ganz erreicht.
SANITY_MIN_M = -100
SANITY_MAX_M = 4700

# Referenzpunkte für die Ausgabe-Kontrolle im CI-Log: (Name, lat, lon, ca. Höhe)
CHECKPOINTS = [
    ("Bern",            46.948,  7.447,  540),
    ("Zuerich",         47.377,  8.540,  410),
    ("Locarno",         46.170,  8.795,  200),
    ("Zermatt",         46.021,  7.749, 1610),
    ("Jungfraujoch",    46.547,  7.985, 3460),
    ("Grosser Aletsch", 46.500,  8.050, 2500),
]


# ─────────────────────────────────────────────────────────────────────────────
# STAC → statisches constants-File
# ─────────────────────────────────────────────────────────────────────────────
def fetch_constants_url() -> str:
    """Href des horizontal_constants-Assets aus der STAC-Collection."""
    r = requests.get(f"{STAC_BASE}/assets", timeout=30)
    r.raise_for_status()
    assets = r.json().get("assets", [])
    for a in assets:
        if a.get("id") == HORIZONTAL_CONSTANTS:
            return a["href"]
    known = ", ".join(sorted(a.get("id", "?") for a in assets))
    raise RuntimeError(
        f"{HORIZONTAL_CONSTANTS} nicht in STAC gefunden. Vorhandene Assets: {known}"
    )


def download(url: str, dest: Path) -> None:
    log.info("Lade %s …", dest.name)
    with requests.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(1 << 16):
                if chunk:
                    f.write(chunk)
    log.info("  %.1f MB", dest.stat().st_size / 1e6)


def load_mesh(path: Path):
    """
    Liest CLAT/CLON/HSURF aus dem GRIB2 über die eccodes-Low-Level-API.

    ICON rechnet auf einem unstrukturierten Dreiecksgitter: pro Zelle je ein
    Lat-, Lon- und Höhenwert, kein reguläres Raster. Rückgabe sind drei
    1D-Arrays gleicher Länge.
    """
    import eccodes

    lat_names = {"clat", "tlat", "rlat", "latitude", "lat"}
    lon_names = {"clon", "tlon", "rlon", "longitude", "lon"}
    hsurf_names = {"hsurf", "h", "orog", "fis"}   # FIS = Geopotential → /g

    clat = clon = hsurf = None
    hsurf_is_geopot = False
    seen: list[str] = []

    with path.open("rb") as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                short = eccodes.codes_get(gid, "shortName")
                seen.append(short)
                lc = short.lower()
                if clat is None and lc in lat_names:
                    clat = np.asarray(eccodes.codes_get_array(gid, "values"), dtype=np.float64)
                elif clon is None and lc in lon_names:
                    clon = np.asarray(eccodes.codes_get_array(gid, "values"), dtype=np.float64)
                elif hsurf is None and lc in hsurf_names:
                    hsurf = np.asarray(eccodes.codes_get_array(gid, "values"), dtype=np.float64)
                    hsurf_is_geopot = (lc == "fis")
            finally:
                eccodes.codes_release(gid)

    log.info("shortNames im constants-File: %s", ", ".join(seen))

    if clat is None or clon is None:
        raise RuntimeError(f"CLAT/CLON nicht gefunden (shortNames: {seen})")
    if hsurf is None:
        raise RuntimeError(f"HSURF nicht gefunden (shortNames: {seen})")

    # ICON liefert die Koordinaten teils in Radiant.
    if np.nanmax(np.abs(clat)) <= math.pi + 0.01:
        log.info("Koordinaten liegen in Radiant vor — rechne nach Grad um.")
        clat = np.degrees(clat)
        clon = np.degrees(clon)

    if hsurf_is_geopot:
        log.info("Feld ist Geopotential (FIS) — teile durch 9.80665.")
        hsurf = hsurf / 9.80665

    log.info(
        "Mesh: %d Zellen | lat %.3f..%.3f | lon %.3f..%.3f | HSURF %.0f..%.0f m",
        clat.size, np.nanmin(clat), np.nanmax(clat),
        np.nanmin(clon), np.nanmax(clon),
        np.nanmin(hsurf), np.nanmax(hsurf),
    )
    return clon, clat, hsurf


# ─────────────────────────────────────────────────────────────────────────────
# Unstrukturiertes Mesh → regelmässiges Raster
# ─────────────────────────────────────────────────────────────────────────────
def build_grid(clon, clat, hsurf, bbox, grid_km):
    """
    Rastert das Mesh per Nearest-Cell auf ein regelmässiges Lon/Lat-Gitter.

    ZEILENREIHENFOLGE: Zeile 0 ist der NÖRDLICHSTE Streifen (Bildkonvention),
    danach absteigend nach Süden. Das steht auch in der JSON-Metadatei und muss
    beim Auslesen beachtet werden — eine verwechselte Y-Achse führt zu einem
    gespiegelten Geländemodell, das im Kartenbild plausibel aussieht und
    trotzdem falsch ist.
    """
    from scipy.spatial import cKDTree

    west, south, east, north = bbox
    lat0 = (south + north) / 2.0

    dlat = grid_km / 111.32
    dlon = grid_km / (111.32 * math.cos(math.radians(lat0)))

    ny = int(round((north - south) / dlat))
    nx = int(round((east - west) / dlon))

    # Zellmittelpunkte
    lats = north - (np.arange(ny) + 0.5) * dlat      # absteigend: Nord → Süd
    lons = west + (np.arange(nx) + 0.5) * dlon

    log.info("Zielraster: %d x %d Zellen (%.2f km) | dlon %.5f° dlat %.5f°",
             nx, ny, grid_km, dlon, dlat)

    # Mesh auf die BBox vorfiltern (mit Rand, damit Randzellen Nachbarn finden)
    pad = 0.2
    m = (
        (clon >= west - pad) & (clon <= east + pad)
        & (clat >= south - pad) & (clat <= north + pad)
    )
    if not np.any(m):
        raise RuntimeError("Kein Mesh-Punkt in der BBox — falsche Collection?")
    log.info("Mesh-Zellen in der BBox: %d von %d", int(m.sum()), clon.size)

    # Lokale äquidistante Projektion, damit Nearest metrisch stimmt
    kx = math.cos(math.radians(lat0))
    tree = cKDTree(np.column_stack([clon[m] * kx, clat[m]]))
    src_h = hsurf[m]

    grid_lon, grid_lat = np.meshgrid(lons, lats)
    query = np.column_stack([grid_lon.ravel() * kx, grid_lat.ravel()])

    # Maximale Suchdistanz: 3 Maschenweiten, sonst nodata
    max_dist = 3.0 * max(dlon * kx, dlat)
    dist, idx = tree.query(query, distance_upper_bound=max_dist)

    out = np.full(query.shape[0], np.nan, dtype=np.float64)
    hit = np.isfinite(dist)
    out[hit] = src_h[idx[hit]]

    misses = int((~hit).sum())
    if misses:
        log.warning("%d Zellen ohne Mesh-Nachbarn → nodata", misses)

    return lons, lats, out.reshape(ny, nx), dlon, dlat


def sanity_check(lons, lats, grid) -> None:
    valid = grid[np.isfinite(grid)]
    if valid.size == 0:
        raise RuntimeError("Raster enthält ausschliesslich nodata.")

    lo, hi = float(valid.min()), float(valid.max())
    log.info("Raster-Höhen: %.0f .. %.0f m (Median %.0f)", lo, hi, float(np.median(valid)))

    if lo < SANITY_MIN_M or hi > SANITY_MAX_M:
        raise RuntimeError(
            f"Höhen ausserhalb des Plausibelbereichs ({lo:.0f}..{hi:.0f} m) — "
            "vermutlich falsche Einheit oder falsches Feld."
        )

    log.info("Kontrollpunkte (erwartet ist die Grössenordnung, nicht der exakte Wert):")
    for name, lat, lon, expect in CHECKPOINTS:
        iy = int(np.argmin(np.abs(lats - lat)))
        ix = int(np.argmin(np.abs(lons - lon)))
        got = grid[iy, ix]
        got_s = "nodata" if not np.isfinite(got) else f"{got:6.0f} m"
        log.info("  %-16s erwartet ~%4d m   gelesen %s", name, expect, got_s)


def write_output(lons, lats, grid, dlon, dlat, bbox, grid_km, prefix: str) -> None:
    west, south, east, north = bbox
    ny, nx = grid.shape

    ints = np.where(np.isfinite(grid), np.rint(grid), NODATA).astype("<i2")

    bin_path = HERE / f"{prefix}.bin"
    bin_path.write_bytes(ints.tobytes(order="C"))

    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "MeteoSwiss ICON-CH1-EPS horizontal_constants (HSURF)",
        "source_collection": "ch.meteoschweiz.ogd-forecasting-icon-ch1",
        "licence": "MeteoSwiss Open Government Data",
        "variable": "surface_altitude",
        "unit": "m",
        "dtype": "int16",
        "byte_order": "little",
        "nodata": NODATA,
        "row_order": "north_to_south",
        "col_order": "west_to_east",
        "note": (
            "Zeile 0 ist der noerdlichste Streifen. Werte sind Zellmittelpunkte: "
            "lat = north - (row + 0.5) * dlat, lon = west + (col + 0.5) * dlon."
        ),
        "bbox": {"west": west, "south": south, "east": east, "north": north},
        "nx": nx,
        "ny": ny,
        "dlon": round(dlon, 8),
        "dlat": round(dlat, 8),
        "grid_km": grid_km,
        "bytes": bin_path.stat().st_size,
    }

    json_path = HERE / f"{prefix}.json"
    json_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    log.info("Geschrieben: %s (%.0f KB) + %s",
             bin_path.name, bin_path.stat().st_size / 1024, json_path.name)


def main() -> int:
    grid_km = float(os.environ.get("GRID_KM", "1.0"))
    prefix = os.environ.get("OUT_PREFIX", "hsurf_ch")

    bbox_raw = os.environ.get("BBOX", "").strip()
    if bbox_raw:
        parts = [float(p) for p in bbox_raw.split(",")]
        if len(parts) != 4:
            log.error("BBOX braucht vier Werte: west,south,east,north")
            return 1
        bbox = tuple(parts)
    else:
        bbox = DEFAULT_BBOX

    log.info("BBox %s | Raster %.2f km | Prefix %s", bbox, grid_km, prefix)

    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            url = fetch_constants_url()
            path = tmp / HORIZONTAL_CONSTANTS
            download(url, path)
            clon, clat, hsurf = load_mesh(path)

        lons, lats, grid, dlon, dlat = build_grid(clon, clat, hsurf, bbox, grid_km)
        sanity_check(lons, lats, grid)
        write_output(lons, lats, grid, dlon, dlat, bbox, grid_km, prefix)

    except Exception as e:                                   # noqa: BLE001
        log.error("Fehlgeschlagen: %s", e)
        return 1

    log.info("Fertig.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
