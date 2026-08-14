#!/usr/bin/env python3
"""
Tests für uv_model.py — laufen ohne ADS-Zugang, ohne eccodes, ohne Netz.

Aufruf:  python test_uv_model.py     (Exit 0 = alles grün)
Im CI vor jedem Lauf, damit ein Rechenfehler nicht erst im Kartenbild auffällt.
"""
from __future__ import annotations

import sys

from uv_model import (
    UVBED_TO_UVI,
    altitude_factor,
    cloud_modification_factor,
    deaccumulate,
    encode_byte,
    equation_of_time_minutes,
    relative_sunshine,
    snow_factor,
    solar_noon_utc_hour,
    uv_index,
)

failures: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        failures.append(name)


def close(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol


print("Zeitgleichung und Sonnenmittag")
# Die Zeitgleichung schwankt übers Jahr etwa zwischen -14 und +16 Minuten.
eots = [equation_of_time_minutes(d) for d in range(1, 366)]
check("EoT im plausiblen Bereich", -17 < min(eots) and max(eots) < 17,
      f"min {min(eots):.1f} max {max(eots):.1f}")

# Schweiz, Mitte August: Sonnenmittag gegen 11:30 UTC (Sommerzeit 13:30 lokal).
noon = solar_noon_utc_hour(226, 8.2275)
check("Sonnenmittag CH im August ~11.5 UTC", 11.2 < noon < 11.8, f"{noon:.2f}")

# Weiter östlich ist der Sonnenmittag früher.
check("Sonnenmittag im Osten frueher",
      solar_noon_utc_hour(226, 10.5) < solar_noon_utc_hour(226, 6.0))

print("\nEntmittelung kumulierter Felder")
check("Differenz zweier Kumulativwerte", close(deaccumulate(3600, 1800), 1800))
check("Negative Differenz wird geklemmt", close(deaccumulate(100, 500), 0.0))

print("\nRelative Sonnenscheindauer")
# Zwei Stunden Fenster, davon eine Stunde Sonne.
check("halbes Fenster sonnig",
      close(relative_sunshine(0, 3600, 0, 7200), 0.5))
check("voll sonnig", close(relative_sunshine(0, 7200, 0, 7200), 1.0))
check("bedeckt", close(relative_sunshine(1000, 1000, 0, 7200), 0.0))
check("Ueberschuss wird auf 1 geklemmt",
      close(relative_sunshine(0, 9000, 0, 7200), 1.0))
check("Nacht (Nenner 0) ergibt 0",
      close(relative_sunshine(0, 0, 500, 500), 0.0))

print("\nBewoelkungsfaktor")
check("klar -> 1.0", close(cloud_modification_factor(1.0), 1.0))
check("bedeckt -> Sockel 0.25", close(cloud_modification_factor(0.0), 0.25))
check("halb -> 0.625", close(cloud_modification_factor(0.5), 0.625))
check("monoton steigend",
      all(cloud_modification_factor(i / 10) < cloud_modification_factor((i + 1) / 10)
          for i in range(10)))

print("\nHoehenfaktor")
check("auf Referenzhoehe neutral", close(altitude_factor(1500, 1500), 1.0))
check("1000 m ueber Referenz -> +7 %", close(altitude_factor(2500, 1500), 1.07))
check("unter Referenz -> kleiner 1", altitude_factor(500, 1500) < 1.0)
# Gegen Meeresniveau statt gegen die CAMS-Orografie waere der Effekt am
# Gipfel mehr als doppelt so gross — das ist der Doppelzaehl-Fehler.
check("Bezugsfehler waere erheblich",
      altitude_factor(3000, 0) - altitude_factor(3000, 1300) > 0.08)

print("\nSchneefaktor")
check("kein Schnee -> 1.0", close(snow_factor(0), 1.0))
check("volle Decke -> 1.25", close(snow_factor(100), 1.25))
check("halbe Decke -> 1.125", close(snow_factor(50), 1.125))
check("Werte ueber 100 % geklemmt", close(snow_factor(150), 1.25))

print("\nUV-Index gesamt")
# Klarhimmel-Dosisleistung, die einem UVI 8 auf Referenzhoehe entspricht.
uvbed_8 = 8.0 / UVBED_TO_UVI

u = uv_index(uvbed_8, 1300, 1300, 1.0, 0.0)
check("klar, auf Referenzhoehe, kein Schnee -> UVI 8", close(u, 8.0, 1e-9), f"{u:.3f}")

u_high = uv_index(uvbed_8, 3300, 1300, 1.0, 0.0)
check("2000 m hoeher -> +14 %", close(u_high, 8.0 * 1.14, 1e-9), f"{u_high:.3f}")

u_cloud = uv_index(uvbed_8, 1300, 1300, 0.0, 0.0)
check("bedeckt -> auf Sockel", close(u_cloud, 8.0 * 0.25, 1e-9), f"{u_cloud:.3f}")

u_snow = uv_index(uvbed_8, 1300, 1300, 1.0, 100.0)
check("volle Schneedecke -> +25 %", close(u_snow, 10.0, 1e-9), f"{u_snow:.3f}")

check("Schnee abschaltbar",
      close(uv_index(uvbed_8, 1300, 1300, 1.0, 100.0, use_snow=False), 8.0, 1e-9))
check("nie negativ", uv_index(0.0, 0, 4000, 0.0, 0.0) >= 0.0)

# Groessenordnung: Gipfel im Hochsommer, klar, mit Schnee.
peak = uv_index(9.0 / UVBED_TO_UVI, 3500, 1300, 1.0, 80.0)
check("Gipfelwert im plausiblen Bereich 10..14", 10.0 <= peak <= 14.0, f"{peak:.2f}")

print("\nByte-Kodierung")
check("UVI 8.0 -> 80", encode_byte(8.0) == 80)
check("UVI 0 -> 0", encode_byte(0.0) == 0)
check("Rundung auf 0.1", encode_byte(7.26) == 73)
check("Saettigung bei 25.5", encode_byte(99.0) == 255)
check("negativ -> 0", encode_byte(-1.0) == 0)
# Verlustfreiheit auf dem relevanten Bereich
check("Rundtrip 0..15 auf 0.1 genau",
      all(abs(encode_byte(v / 10) / 10 - v / 10) < 0.051 for v in range(0, 151)))

print()
if failures:
    print(f"{len(failures)} Test(s) fehlgeschlagen: {', '.join(failures)}")
    sys.exit(1)

print("Alle Tests bestanden.")
sys.exit(0)
