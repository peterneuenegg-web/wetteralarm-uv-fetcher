# wetteralarm-uv-fetcher

Externer Worker für das Projekt **UV-Index Schweiz**. Holt Daten aus Quellen,
die auf Infomaniak Shared Hosting nicht verarbeitbar sind, und stellt sie dem
UV-Backend bereit.

**Stand:** Stufe 1 umgesetzt (Geländehöhen-Raster). Stufe 2 geplant.

---

## Warum es dieses Repo gibt

Infomaniak Shared Hosting hat **kein eccodes**, und `exec()` ist per
`disable_functions` gesperrt. GRIB2 lässt sich dort nicht dekodieren. Muster
wie [`weitsicht-forecast-fetcher`](https://github.com/peterneuenegg-web/weitsicht-forecast-fetcher)
und `wetteralarm-hail-fetcher`.

**Eigenes Repo statt Mitbenutzung von Weitsicht:** Die beiden Projekte haben
unterschiedliche Lebenszyklen. Würde Weitsicht eingestellt oder auf andere
Referenzpunkte umgebaut, müsste jemand daran denken, dass UV mit daran hängt.
Der Preis dieser Entkopplung: Fehlerkorrekturen an den gemeinsamen Codemustern
müssen in beiden Repos nachgezogen werden.

## Stufe 1 — Geländehöhen-Raster (umgesetzt)

`fetch_hsurf.py` liest die Geländehöhe `HSURF` aus dem **statischen**
`horizontal_constants`-File der ICON-CH1-Collection und rastert sie auf ein
regelmässiges Lon/Lat-Gitter über der Schweiz.

**ICON-CH1 statt CH2:** 1 km statt 2,1 km Maschenweite. Der kürzere
Vorhersagehorizont von CH1 (33 h gegenüber 120 h) ist hier bedeutungslos —
Gelände hat keinen Vorhersagehorizont.

### Keine Laufzeitabhängigkeit

Gelände ist statisch. Der Workflow läuft **nur manuell** und muss nach der
Erstellung praktisch nie wieder laufen. Sobald die Ausgabe im UV-Projekt liegt,
ist dieses Repo für Stufe 1 ein reines Build-Werkzeug — das UV-Backend
funktioniert auch dann weiter, wenn es dieses Repo nicht mehr gäbe.

### Ausgabe

| Datei | Inhalt |
|---|---|
| `hsurf_ch.bin` | int16 Little-Endian, Höhe in Metern, zeilenweise, ~160 KB |
| `hsurf_ch.json` | BBox, Dimensionen, `dlon`/`dlat`, nodata, Zeilenreihenfolge |

**Zeile 0 ist der nördlichste Streifen** (Bildkonvention), Spalte 0 die
westlichste. Zellmittelpunkte:

```
lat = north - (row + 0.5) * dlat
lon = west  + (col + 0.5) * dlon
```

> Eine verwechselte Y-Achse liefert ein gespiegeltes Geländemodell, das im
> Kartenbild plausibel aussieht und trotzdem falsch ist. Genau dieser Fehler
> hat im Hagel-Worker einmal einen 60-km-Versatz erzeugt. Die Reihenfolge steht
> deshalb zusätzlich in `hsurf_ch.json` unter `row_order`.

Standard-BBox `5.90, 45.75, 10.55, 47.85` (west, south, east, north),
1 km Raster → 354 × 234 Zellen. Fehlende Werte sind `-32768`.

### Auslesen in PHP

```php
$meta = json_decode(file_get_contents('hsurf_ch.json'), true, 512, JSON_THROW_ON_ERROR);
$raw  = file_get_contents('hsurf_ch.bin');
$h    = unpack('s*', $raw);   // 1-basiert!

function heightAt(array $meta, array $h, float $lat, float $lon): ?int {
    $b = $meta['bbox'];
    $row = (int) floor(($b['north'] - $lat) / $meta['dlat']);
    $col = (int) floor(($lon - $b['west']) / $meta['dlon']);
    if ($row < 0 || $row >= $meta['ny'] || $col < 0 || $col >= $meta['nx']) return null;
    $v = $h[$row * $meta['nx'] + $col + 1];
    return $v === $meta['nodata'] ? null : $v;
}
```

`unpack('s*')` liest Host-Byte-Order. Auf x86 und ARM ist das Little-Endian und
passt; auf einer Big-Endian-Maschine müsste stattdessen `v*` mit
Vorzeichenkorrektur verwendet werden.

### Ausführen

Actions-Tab → **fetch-hsurf** → *Run workflow*. Optional `grid_km` und `bbox`
überschreiben. Danach das Artifact `hsurf-ch` herunterladen und beide Dateien

- in dieses Repo committen (Nachvollziehbarkeit) und
- nach `UV/Applikation/data/` im UV-Projekt kopieren.

Das CI-Log gibt Kontrollpunkte aus (Bern ~540 m, Zermatt ~1610 m,
Jungfraujoch ~3460 m). Weichen die stark ab, stimmt etwas mit Einheit,
Feldauswahl oder Achsenreihenfolge nicht — dann **nicht** deployen.

### Lokal testen

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python fetch_hsurf.py
```

Braucht Python 3.10–3.12 (**nicht** 3.13, `meteodata-lab 0.7.2` unterstützt es
nicht). Kein Token nötig — die OGD-Daten sind offen.

## Stufe 2 — eigenes UV-Feld (Code steht, ungetestet)

Rechnet den UV-Index flächig statt ihn aus Talstationen zu interpolieren.
Vollständige Spezifikation im UV-Projekt: `Doku/SPEC_uv-raster.md`.

```
UVI = UVI_klar(CAMS) · Höhenterm(HSURF) · Bewölkung(DURSUN) · Schnee(SNOWC)
```

| Baustein | Quelle | Auflösung |
|---|---|---|
| UV klar (Ozon, Sonnenstand, Aerosol) | CAMS `uv_biologically_effective_dose_clear_sky` | ~45 km |
| Bewölkung | ICON `DURSUN` / `DURSUN_M` | 1 km / 2,1 km |
| Höhe | `HSURF` aus Stufe 1 | 1 km |
| Schnee | ICON `SNOWC` | 1 km / 2,1 km |

Modellwahl je Tag: **ICON-CH1 für Tag 1** (33 h Horizont, 1 km),
**ICON-CH2 für Tag 2–5** (120 h, 2,1 km). Weiter als 5 Tage reicht keine der
beiden Quellen.

### Dateien

| Datei | Zweck |
|---|---|
| `uv_model.py` | Reine Rechenlogik, **kein Netzzugriff** |
| `test_uv_model.py` | 27 Zusicherungen dazu, laufen ohne ADS und ohne eccodes |
| `fetch_uv.py` | Abruf, Regridding, Modell, Versand |

Die Physik ist bewusst von der Netzlogik getrennt: Sie ist der Teil, der
stillschweigend falsch sein kann, und gehört deshalb prüfbar. `test_uv_model.py`
läuft im Workflow **vor** jedem Lauf.

### Voraussetzung

`fetch_uv.py` liest `hsurf_ch.bin` und `hsurf_ch.json` aus diesem Repo — sie
definieren das Zielgitter. Sie müssen also committet sein (Ergebnis von
`fetch-hsurf`), sonst bricht der Lauf sofort ab.

### Secrets

| Secret | Zweck |
|---|---|
| `ADS_API_KEY` | Copernicus-ADS-Token. **Konto geteilt mit Weitsicht** — bei Rotation beide Repos nachziehen |
| `UV_INGEST_URL_STAGE` / `_PROD` | Ziel-Endpoint `…/api/raster-ingest.php` |
| `UV_INGEST_TOKEN_STAGE` / `_PROD` | Muss `RASTER_INGEST_TOKEN` in der Infomaniak-`.env` matchen |

### Erster Lauf

Mit `dry_run: true` starten. Dann wird gerechnet und als Artifact
hochgeladen, aber nichts gesendet. Im Log prüfen:

- `CAMS NetCDF: vars=…` — heissen die Variablen wie erwartet?
- `rel_sun … (Mittel …)` — plausibel zwischen 0 und 1?
- `UVI {min, max, mean}` — Mittelland im Hochsommer 6–8, Gipfel 9–11?
- `Flächenmittel eigen … vs CAMS bewölkt …` — **Faktor sollte nahe 1 liegen.**
  Weicht er stark ab, stimmt der Bewölkungsfaktor oder die Einheitenumrechnung
  nicht. Die ADS-Doku nennt die Einheit von `uvbedcs` „dimensionless", der
  Code rechnet mit W/m² und Faktor 40 — genau das ist hier zu verifizieren.

Erst wenn diese vier Punkte stimmen, ohne `dry_run` laufen lassen.

## Betrieb

**Trigger.** Sollte Stufe 2 kommen, wird der Workflow **nicht** per
`schedule:` laufen, sondern vom Infomaniak-Cron des UV-Projekts per
`workflow_dispatch` angestossen — Vorlage
`Weitsicht/cron/trigger_forecast_fetcher.php`. Zwei Gründe:

- GitHub-Cron ist Best-Effort; im Hagel-Worker wurden bei einem
  5-Minuten-Schedule real nur alle 1–3 Stunden Läufe beobachtet.
- GitHub deaktiviert in öffentlichen Repos **geplante** Workflows nach 60 Tagen
  ohne Repository-Aktivität. Ohne `schedule:` greift die Regel nicht.

**Kosten.** Öffentliches Repo auf Standard-Runnern: keine. Copernicus ADS und
MeteoSchweiz OGD sind ebenfalls kostenlos.

**Sicherheit.** Stufe 1 braucht keine Secrets. Für Stufe 2 gilt: Secrets
ausschliesslich als Repository-Secrets, nie in Logs (`::add-mask::`), Actions
auf vollen Commit-SHA gepinnt, `permissions: contents: read` als Default.

**Geteilte Abhängigkeit trotz getrennter Repos.** Nutzt Stufe 2 denselben
Copernicus-ADS-Account wie Weitsicht, hängen beide Projekte an diesem Konto —
wird der Token rotiert oder das Konto gesperrt, stehen beide still. Entweder
ein eigenes ADS-Konto für UV anlegen oder die Abhängigkeit in beiden READMEs
festhalten.

## Verwandte Repos

| Repo | Zweck |
|---|---|
| `UV` | Das Backend, das die Ausgabe konsumiert |
| `weitsicht-forecast-fetcher` | Vorlage für Codemuster und Trigger-Anbindung |
| `wetteralarm-hail-fetcher` | Ursprung des Worker-Musters |
