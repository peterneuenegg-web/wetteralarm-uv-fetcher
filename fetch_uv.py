#!/usr/bin/env python3
"""
UV-Raster (Stufe 2) — CAMS + ICON → 1-km-UV-Index-Feld für die Schweiz.

Rechnet den UV-Index flächendeckend, statt ihn aus Ortschaften zu
interpolieren. Modell und Begründung: UV/Doku/SPEC_uv-raster.md.

    UVI = UVI_klar(CAMS) · Höhenterm(ICON HSURF) · Bewölkung(ICON DURSUN) · Schnee

Ausgabe je Prognosetag eine uint8-Rasterdatei auf demselben Gitter wie
hsurf_ch.bin, danach POST an den Ingest-Endpoint des UV-Projekts.

Warum externer Worker: Infomaniak hat kein eccodes, exec() ist gesperrt.
Defensive xarray-Helfer übernommen aus weitsicht-forecast-fetcher — die
exakten Dim-/Koordinatennamen von meteodata-lab sind dort verifiziert.

Env (GitHub Secrets):
  ADS_API_KEY         Copernicus-ADS-Token (Konto geteilt mit Weitsicht)
  ADS_URL             optional, default https://ads.atmosphere.copernicus.eu/api
  UV_INGEST_URL       z.B. https://tool.wetteralarm.ch/uv/stage/api/raster-ingest.php
  UV_INGEST_TOKEN     muss RASTER_INGEST_TOKEN in der Infomaniak .env matchen
  FORECAST_DAYS       optional, default 5 (CAMS und ICON-CH2 reichen 120 h)
  USE_SNOW            optional, "false" schaltet den Albedo-Term ab
  DRY_RUN             optional, "true" = rechnen und schreiben, nicht senden

Exit-Codes: 0 ok · 1 Fehler · 2 Config fehlt
"""
from __future__ import annotations

import json
import logging
import math
import os
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import requests

from uv_model import encode_byte, solar_noon_utc_hour, uv_index

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("fetch-uv")

HERE = Path(__file__).resolve().parent

CH_AREA = [48.0, 5.5, 45.5, 11.0]          # N, W, S, E für CAMS
CAMS_DATASET = "cams-global-atmospheric-composition-forecasts"
VAR_UV_CS = "uv_biologically_effective_dose_clear_sky"
VAR_UV_ALL = "uv_biologically_effective_dose"

# Fenster um den Sonnenhöchststand, über das die Sonnenscheindauer
# ausgewertet wird. Ein Punktwert wäre anfällig für einzelne Wolken.
SUN_WINDOW_HOURS = 1

# Referenzhöhe: HSURF auf CAMS-Skala geglättet. 45 km bei ~1 km Zellen.
SMOOTH_CELLS = 45


# ─────────────────────────────────────────────────────────────────────────────
# Zielgitter
# ─────────────────────────────────────────────────────────────────────────────
def load_grid():
    """Gitter und Geländehöhe aus der Ausgabe von fetch_hsurf.py."""
    meta = json.loads((HERE / "hsurf_ch.json").read_text(encoding="utf-8"))
    raw = (HERE / "hsurf_ch.bin").read_bytes()

    expected = meta["nx"] * meta["ny"] * 2
    if len(raw) != expected:
        raise RuntimeError(f"hsurf_ch.bin: {len(raw)} Bytes, erwartet {expected}")
    if meta.get("row_order") != "north_to_south":
        raise RuntimeError(f"unerwartete row_order: {meta.get('row_order')}")

    hsurf = np.frombuffer(raw, dtype="<i2").astype(np.float64).reshape(meta["ny"], meta["nx"])
    hsurf = np.where(hsurf == meta["nodata"], np.nan, hsurf)

    b = meta["bbox"]
    lats = b["north"] - (np.arange(meta["ny"]) + 0.5) * meta["dlat"]   # Nord → Süd
    lons = b["west"] + (np.arange(meta["nx"]) + 0.5) * meta["dlon"]

    log.info("Zielgitter %d x %d | HSURF %.0f..%.0f m",
             meta["nx"], meta["ny"], np.nanmin(hsurf), np.nanmax(hsurf))
    return meta, lats, lons, hsurf


def reference_heights(hsurf: np.ndarray) -> np.ndarray:
    """
    HSURF auf CAMS-Skala geglättet — die Höhe, auf der uvbedcs bereits gilt.

    Der Höhenterm bezieht sich auf die DIFFERENZ zu dieser Fläche. Gegen
    Meeresniveau gerechnet käme der Höheneffekt doppelt hinein, weil CAMS ihn
    auf seiner eigenen groben Orografie schon zum Teil enthält.

    Gleitendes Mittel per Summed-Area-Table, NaN-tolerant.
    """
    valid = np.isfinite(hsurf)
    filled = np.where(valid, hsurf, 0.0)

    def boxsum(a):
        c = np.cumsum(np.cumsum(a, axis=0), axis=1)
        c = np.pad(c, ((1, 0), (1, 0)), mode="constant")
        ny, nx = a.shape
        r = SMOOTH_CELLS // 2
        y0 = np.clip(np.arange(ny) - r, 0, ny)
        y1 = np.clip(np.arange(ny) + r + 1, 0, ny)
        x0 = np.clip(np.arange(nx) - r, 0, nx)
        x1 = np.clip(np.arange(nx) + r + 1, 0, nx)
        return (c[np.ix_(y1, x1)] - c[np.ix_(y0, x1)]
                - c[np.ix_(y1, x0)] + c[np.ix_(y0, x0)])

    total = boxsum(filled)
    count = boxsum(valid.astype(np.float64))

    with np.errstate(invalid="ignore", divide="ignore"):
        ref = np.where(count > 0, total / count, np.nan)

    log.info("Referenzhöhen (CAMS-Skala): %.0f..%.0f m", np.nanmin(ref), np.nanmax(ref))
    return ref


# ─────────────────────────────────────────────────────────────────────────────
# CAMS
# ─────────────────────────────────────────────────────────────────────────────
def latest_cams_run() -> tuple[str, str]:
    """Neuester 00/12-UTC-Lauf, mit ~8 h Publikationsversatz (wie Weitsicht)."""
    override = os.environ.get("REFERENCE_DATETIME", "").strip()
    if override:
        dt = datetime.strptime(override, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    else:
        now = datetime.now(timezone.utc) - timedelta(hours=8)
        dt = now.replace(hour=12 if now.hour >= 12 else 0, minute=0, second=0, microsecond=0)
    return dt.strftime("%Y-%m-%d"), f"{dt.hour:02d}:00"


def download_cams(date_str: str, run_time: str, leadtimes: list[int], dest: Path) -> None:
    import cdsapi

    # Leeres Secret ("" statt fehlend) auf den Default fallen lassen — sonst
    # baut cdsapi eine URL ohne Schema.
    ads_url = os.environ.get("ADS_URL", "").strip() or "https://ads.atmosphere.copernicus.eu/api"
    key = os.environ.get("ADS_API_KEY", "").strip()
    if not key:
        raise SystemExit(2)

    client = cdsapi.Client(url=ads_url, key=key)
    request = {
        "variable": [VAR_UV_CS, VAR_UV_ALL],
        "date": f"{date_str}/{date_str}",
        "time": run_time,
        "leadtime_hour": [str(h) for h in leadtimes],
        "type": "forecast",
        "area": CH_AREA,
        # Neue CADS/ADS-API erwartet data_format; "format" ist veraltet.
        "data_format": "netcdf",
    }
    log.info("CAMS: %s %s leadtimes=%s", date_str, run_time, leadtimes)
    client.retrieve(CAMS_DATASET, request, str(dest))

    if not dest.exists() or dest.stat().st_size == 0:
        raise RuntimeError("CAMS-Download leer")


def pick_var(ds, *needles):
    """Variablenname im NetCDF finden — CAMS benennt teils um."""
    for name in ds.data_vars:
        low = name.lower()
        if all(n in low for n in needles):
            return name
    return None


def cams_field(ds, varname: str, lead_index: int, lats, lons) -> np.ndarray:
    """Ein CAMS-Feld bilinear auf das Zielgitter bringen."""
    da = ds[varname]

    # Lead-Dimension herausgreifen, alle übrigen Nicht-Raum-Dims auf 0.
    for dim in list(da.dims):
        if dim in ("latitude", "lat", "longitude", "lon"):
            continue
        da = da.isel({dim: min(lead_index, da.sizes[dim] - 1)}) if dim in ("forecast_period", "step", "time", "valid_time") else da.isel({dim: 0})
        lead_index = 0

    src_lat = np.asarray(da["latitude" if "latitude" in da.coords else "lat"].values, dtype=np.float64)
    src_lon = np.asarray(da["longitude" if "longitude" in da.coords else "lon"].values, dtype=np.float64)
    values = np.asarray(da.values, dtype=np.float64)

    # Für np.interp müssen die Achsen aufsteigend sein.
    if src_lat[0] > src_lat[-1]:
        src_lat = src_lat[::-1]
        values = values[::-1, :]
    if src_lon[0] > src_lon[-1]:
        src_lon = src_lon[::-1]
        values = values[:, ::-1]

    # Bilinear: erst entlang lon je Quellzeile, dann entlang lat.
    tmp = np.empty((values.shape[0], lons.size), dtype=np.float64)
    for i in range(values.shape[0]):
        tmp[i, :] = np.interp(lons, src_lon, values[i, :])

    out = np.empty((lats.size, lons.size), dtype=np.float64)
    for j in range(lons.size):
        out[:, j] = np.interp(lats, src_lat, tmp[:, j])

    return out


# ─────────────────────────────────────────────────────────────────────────────
# ICON
# ─────────────────────────────────────────────────────────────────────────────
def _ogd():
    from meteodatalab.ogd_api import Collection, Request, get_from_ogd
    return Request, Collection, get_from_ogd


def icon_latest_ref(collection_id: str) -> str:
    """Neuesten verfügbaren Lauf aus dem STAC-Katalog lesen."""
    url = (f"https://data.geo.admin.ch/api/stac/v1/collections/"
           f"{collection_id}/items?limit=1")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    feats = r.json().get("features", [])
    if not feats:
        raise RuntimeError(f"kein Lauf in {collection_id}")

    props = feats[0].get("properties", {})
    ref = props.get("forecast:reference_datetime") or props.get("datetime")
    if not ref:
        raise RuntimeError(f"keine reference_datetime in {collection_id}")
    return str(ref).split("/")[0]


def icon_field(variable: str, collection, ref_iso: str, horizon: timedelta,
               lats, lons, tree_cache: dict) -> np.ndarray:
    """
    Ein ICON-Feld per Nearest-Cell auf das Zielgitter bringen.

    ICON rechnet auf einem unstrukturierten Dreiecksgitter — deshalb KDTree
    statt Interpolation auf regulären Achsen.
    """
    from scipy.spatial import cKDTree

    Request, _, get_from_ogd = _ogd()
    da = get_from_ogd(Request(
        collection=collection,
        variable=variable,
        reference_datetime=ref_iso,
        perturbed=False,
        horizon=horizon,
    ))

    lon_src, lat_src = find_lonlat(da)
    cell_dim = find_cell_dim(da, lon_src.size)

    sel = {d: 0 for d in da.dims if d != cell_dim}
    arr = np.asarray((da.isel(**sel) if sel else da).values, dtype=np.float64).ravel()

    key = (lon_src.size, round(float(lon_src[0]), 6))
    if key not in tree_cache:
        lat0 = float(np.mean(lat_src))
        kx = math.cos(math.radians(lat0))
        tree = cKDTree(np.column_stack([lon_src * kx, lat_src]))
        gl, gt = np.meshgrid(lons, lats)
        _, idx = tree.query(np.column_stack([gl.ravel() * kx, gt.ravel()]))
        tree_cache[key] = idx
    idx = tree_cache[key]

    return arr[idx].reshape(lats.size, lons.size)


def find_lonlat(da):
    """(lon, lat) der Mesh-Zellen als 1D float64. Aus weitsicht-forecast-fetcher."""
    lon = lat = None
    for cand in ("longitude", "lon", "clon"):
        if cand in da.coords:
            lon = np.asarray(da.coords[cand].values, dtype=np.float64).ravel()
            break
    for cand in ("latitude", "lat", "clat"):
        if cand in da.coords:
            lat = np.asarray(da.coords[cand].values, dtype=np.float64).ravel()
            break
    if lon is None or lat is None:
        raise RuntimeError(f"lon/lat nicht gefunden; coords={list(da.coords.keys())}")
    if np.nanmax(np.abs(lat)) <= math.pi + 0.01:
        lon = np.degrees(lon)
        lat = np.degrees(lat)
    return lon, lat


def find_cell_dim(da, ncells: int) -> str:
    for d in da.dims:
        if da.sizes[d] == ncells:
            return d
    raise RuntimeError(f"Zell-Dimension ({ncells}) nicht in dims={da.dims}")


# ─────────────────────────────────────────────────────────────────────────────
# Hauptlauf
# ─────────────────────────────────────────────────────────────────────────────
def main() -> int:
    days = int(os.environ.get("FORECAST_DAYS", "5"))
    use_snow = os.environ.get("USE_SNOW", "true").strip().lower() != "false"
    dry_run = os.environ.get("DRY_RUN", "").strip().lower() == "true"

    meta, lats, lons, hsurf = load_grid()
    href = reference_heights(hsurf)
    lon_centre = float(np.mean(lons))

    cams_date, cams_time = latest_cams_run()
    run_dt = datetime.strptime(f"{cams_date} {cams_time}", "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)

    # Zielzeitpunkte: Sonnenhöchststand jedes Prognosetags.
    targets = []
    for d in range(days):
        day = (run_dt + timedelta(days=d)).date()
        noon_h = solar_noon_utc_hour(day.timetuple().tm_yday, lon_centre)
        target = datetime.combine(day, datetime.min.time(), tzinfo=timezone.utc) \
                 + timedelta(hours=round(noon_h))
        lead = int(round((target - run_dt).total_seconds() / 3600))
        if lead < 0:
            log.info("Tag %s liegt vor dem Lauf — übersprungen.", day)
            continue
        targets.append({"date": day.isoformat(), "lead": lead, "target": target})

    if not targets:
        log.error("Keine gültigen Zielzeitpunkte.")
        return 1

    log.info("Zielzeitpunkte: %s",
             ", ".join(f"{t['date']}@+{t['lead']}h" for t in targets))

    # ---------- CAMS ----------
    leads = sorted({t["lead"] for t in targets})
    with tempfile.TemporaryDirectory() as td:
        nc = Path(td) / "cams.nc"
        download_cams(cams_date, cams_time, leads, nc)

        import xarray as xr
        ds = xr.open_dataset(nc)
        log.info("CAMS NetCDF: vars=%s dims=%s", list(ds.data_vars), dict(ds.sizes))

        name_cs = pick_var(ds, "uvbedcs") or pick_var(ds, "clear")
        name_all = pick_var(ds, "uvbed")
        if name_cs is None:
            raise RuntimeError(f"uvbedcs nicht im NetCDF: {list(ds.data_vars)}")
        log.info("CAMS-Variablen: klar=%s bewölkt=%s", name_cs, name_all)

        for t in targets:
            t["uv_cs"] = cams_field(ds, name_cs, leads.index(t["lead"]), lats, lons)
            if name_all:
                t["uv_all"] = cams_field(ds, name_all, leads.index(t["lead"]), lats, lons)
        ds.close()

    # ---------- ICON ----------
    _, Collection, _ = _ogd()
    tree_cache: dict = {}

    ref_ch1 = icon_latest_ref("ch.meteoschweiz.ogd-forecasting-icon-ch1")
    ref_ch2 = icon_latest_ref("ch.meteoschweiz.ogd-forecasting-icon-ch2")
    ref_ch1_dt = datetime.fromisoformat(ref_ch1.replace("Z", "+00:00"))
    ref_ch2_dt = datetime.fromisoformat(ref_ch2.replace("Z", "+00:00"))
    log.info("ICON-Läufe: CH1 %s | CH2 %s", ref_ch1, ref_ch2)

    for t in targets:
        # CH1 hat 1 km, reicht aber nur 33 h — danach CH2 mit 2.1 km.
        lead_ch1 = (t["target"] - ref_ch1_dt).total_seconds() / 3600
        if 0 <= lead_ch1 <= 33:
            coll, ref_dt, label = Collection.ICON_CH1, ref_ch1_dt, "ICON-CH1"
        else:
            coll, ref_dt, label = Collection.ICON_CH2, ref_ch2_dt, "ICON-CH2"

        h_mid = t["target"] - ref_dt
        h_lo = h_mid - timedelta(hours=SUN_WINDOW_HOURS)
        h_hi = h_mid + timedelta(hours=SUN_WINDOW_HOURS)

        if h_lo.total_seconds() < 0:
            h_lo = timedelta(0)

        try:
            d1 = icon_field("DURSUN", coll, ref_dt.isoformat(), h_lo, lats, lons, tree_cache)
            d2 = icon_field("DURSUN", coll, ref_dt.isoformat(), h_hi, lats, lons, tree_cache)
            m1 = icon_field("DURSUN_M", coll, ref_dt.isoformat(), h_lo, lats, lons, tree_cache)
            m2 = icon_field("DURSUN_M", coll, ref_dt.isoformat(), h_hi, lats, lons, tree_cache)

            possible = np.maximum(m2 - m1, 0.0)
            actual = np.maximum(d2 - d1, 0.0)
            with np.errstate(invalid="ignore", divide="ignore"):
                rel = np.where(possible > 0, np.clip(actual / possible, 0.0, 1.0), 0.0)

            snow = icon_field("SNOWC", coll, ref_dt.isoformat(), h_mid, lats, lons, tree_cache) \
                   if use_snow else np.zeros_like(rel)

            t["model"] = label
            t["rel_sun"] = rel
            t["snow"] = snow
            log.info("%s %s: rel_sun %.2f..%.2f (Mittel %.2f)",
                     t["date"], label, float(np.nanmin(rel)), float(np.nanmax(rel)),
                     float(np.nanmean(rel)))
        except Exception as e:                                   # noqa: BLE001
            log.error("ICON für %s fehlgeschlagen: %s", t["date"], e)
            return 1

    # ---------- Modell ----------
    payload_days = []
    for t in targets:
        uvi = np.zeros_like(t["uv_cs"])
        it = np.nditer(uvi, flags=["multi_index"], op_flags=["writeonly"])
        while not it.finished:
            y, x = it.multi_index
            h = hsurf[y, x]
            if np.isnan(h):
                it[0] = 0.0
            else:
                it[0] = uv_index(
                    float(t["uv_cs"][y, x]),
                    float(h),
                    float(href[y, x]) if np.isfinite(href[y, x]) else float(h),
                    float(t["rel_sun"][y, x]),
                    float(t["snow"][y, x]),
                    use_snow=use_snow,
                )
            it.iternext()

        raster = np.vectorize(encode_byte)(uvi).astype(np.uint8)

        fname = f"uvraster_{t['date']}.bin"
        (HERE / fname).write_bytes(raster.tobytes(order="C"))

        stats = {
            "min": round(float(uvi.min()), 2),
            "max": round(float(uvi.max()), 2),
            "mean": round(float(uvi.mean()), 2),
        }
        log.info("%s: UVI %s  -> %s", t["date"], stats, fname)

        # Kontrolle gegen CAMS bewölkt (45 km): Flächenmittel sollte grob passen.
        if "uv_all" in t:
            cams_mean = float(np.mean(t["uv_all"])) * 40.0
            log.info("  Flächenmittel eigen %.2f vs CAMS bewölkt %.2f (Faktor %.2f)",
                     stats["mean"], cams_mean,
                     stats["mean"] / cams_mean if cams_mean else float("nan"))

        payload_days.append({"date": t["date"], "model": t["model"],
                             "file": fname, **stats})

    out_meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cams_run": f"{cams_date}T{cams_time}Z",
        "icon_ch1_run": ref_ch1,
        "icon_ch2_run": ref_ch2,
        "use_snow": use_snow,
        "bbox": meta["bbox"], "nx": meta["nx"], "ny": meta["ny"],
        "dlon": meta["dlon"], "dlat": meta["dlat"],
        "row_order": "north_to_south", "col_order": "west_to_east",
        "dtype": "uint8", "scale": 0.1,
        "note": "Wert = Byte / 10 = UV-Index. Gitter identisch mit hsurf_ch.bin.",
        "days": payload_days,
    }
    (HERE / "uvraster_meta.json").write_text(
        json.dumps(out_meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    log.info("Meta geschrieben: uvraster_meta.json")

    if dry_run:
        log.info("DRY_RUN — kein Versand.")
        return 0

    return send(out_meta, payload_days)


def send(out_meta: dict, days: list[dict]) -> int:
    url = os.environ.get("UV_INGEST_URL", "").strip()
    token = os.environ.get("UV_INGEST_TOKEN", "").strip()
    if not url or not token:
        log.error("UV_INGEST_URL / UV_INGEST_TOKEN fehlen.")
        return 2

    files = {"meta": ("uvraster_meta.json",
                      json.dumps(out_meta).encode("utf-8"), "application/json")}
    for d in days:
        files[d["file"]] = (d["file"], (HERE / d["file"]).read_bytes(),
                            "application/octet-stream")

    r = requests.post(url, files=files, headers={"X-Ingest-Token": token}, timeout=180)
    log.info("Ingest: HTTP %s — %s", r.status_code, r.text[:300])
    return 0 if r.ok else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as e:                                       # noqa: BLE001
        log.error("Fehlgeschlagen: %s", e)
        sys.exit(1)
