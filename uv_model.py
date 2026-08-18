#!/usr/bin/env python3
"""
UV-Modell — reine Rechenlogik, ohne Netzzugriff.

Bewusst von fetch_uv.py getrennt: Diese Funktionen sind ohne ADS-Zugang und
ohne GRIB-Bibliotheken testbar (siehe test_uv_model.py). Die Physik ist der
Teil, der stillschweigend falsch sein kann — der gehört prüfbar.

Modell (Details und Begründung: UV/Doku/SPEC_uv-raster.md):

    UVI = UVI_klar · A(Höhe) · CMF(Bewölkung) · S(Schnee)
"""
from __future__ import annotations

import math

# UV-Zunahme pro 1000 Höhenmeter (Literatur 6–10 %).
ALTITUDE_GRADIENT = 0.07

# Bewölkungsfaktor in Ångström-Form: CMF = CMF_MIN + (1-CMF_MIN) · rel_sun.
# Der Sockel bildet ab, dass UV auch bei geschlossener Decke durch Streuung
# den Boden erreicht — UV wird deutlich weniger gedämpft als Gesamtstrahlung.
#
# KALIBRIERT am 2026-08-18 gegen CAMS: Aus dem Verhältnis uvbed/uvbedcs lässt
# sich der Bewölkungsfaktor ablesen, den CAMS selbst ansetzt. Über fünf Tage
# mit rel_sun von 0.06 bis 0.99 ergibt die Kleinste-Quadrate-Anpassung einen
# Sockel von 0.47 (RMS 0.062). Der ursprüngliche Schätzwert 0.25 dämpfte bei
# geschlossener Decke etwa doppelt so stark wie CAMS und lag auch unter dem
# Literaturbereich (UV unter Bewölkung typisch 30–50 % des Klarhimmelwerts).
#
# Nachkalibrieren, wenn mehr Tage vorliegen: Die Zeile "Bewölkungsfaktor:
# CAMS x vs eigen y" im Worker-Log ist die Datengrundlage.
CMF_MIN = 0.47

# Zuschlag bei vollständiger Schneedecke (Albedo-Rückstreuung).
SNOW_ALBEDO_BONUS = 0.25

# CAMS liefert die UV-Dosisleistung in W/m². Der UV-Index ist per Definition
# das 40-fache davon.
UVBED_TO_UVI = 40.0


def equation_of_time_minutes(day_of_year: int) -> float:
    """
    Zeitgleichung in Minuten (Näherung nach Spencer/Cooper).

    Positiv = wahre Sonnenzeit läuft der mittleren voraus.
    """
    b = 2.0 * math.pi * (day_of_year - 81) / 364.0
    return 9.87 * math.sin(2 * b) - 7.53 * math.cos(b) - 1.5 * math.sin(b)


def solar_noon_utc_hour(day_of_year: int, longitude_deg: float) -> float:
    """
    Wahrer Sonnenhöchststand als UTC-Stunde (Dezimal).

    Der UV-Index ist konventionell das Tagesmaximum, das rund um diesen
    Zeitpunkt auftritt. Gerechnet wird deshalb pro Tag genau ein Zeitschritt,
    nicht der Tagesmittelwert.
    """
    eot_h = equation_of_time_minutes(day_of_year) / 60.0
    return 12.0 - longitude_deg / 15.0 - eot_h


def deaccumulate(value_t2: float, value_t1: float) -> float:
    """
    Differenz zweier kumulierter Werte (Accumulation/Summation seit
    Referenzzeit).

    ICON liefert DURSUN und DURSUN_M kumulativ ab Modellstart. Wer die
    Rohwerte direkt verwendet, mittelt unbemerkt über den ganzen Vorlauf und
    bekommt dadurch zu tiefe Mittagswerte.

    Negative Differenzen (Rundung, Modellartefakte) werden auf 0 geklemmt.
    """
    return max(0.0, value_t2 - value_t1)


def relative_sunshine(dursun_t1: float, dursun_t2: float,
                      dursun_max_t1: float, dursun_max_t2: float) -> float:
    """
    Relative Sonnenscheindauer im Fenster [t1, t2], geklemmt auf 0..1.

    Nenner ist das astronomisch mögliche Maximum im selben Fenster. Ist es
    null (Nacht, Polarwinter), gibt es keine Information — dann 0.
    """
    actual = deaccumulate(dursun_t2, dursun_t1)
    possible = deaccumulate(dursun_max_t2, dursun_max_t1)

    if possible <= 0.0:
        return 0.0

    return min(1.0, actual / possible)


def cloud_modification_factor(rel_sunshine: float) -> float:
    """Bewölkungsfaktor aus relativer Sonnenscheindauer."""
    r = min(1.0, max(0.0, rel_sunshine))
    return CMF_MIN + (1.0 - CMF_MIN) * r


def altitude_factor(height_m: float, reference_height_m: float,
                    gradient: float = ALTITUDE_GRADIENT) -> float:
    """
    Höhenfaktor gegenüber der Referenzhöhe.

    WICHTIG: Bezug ist die Höhe, auf der das Klarhimmel-Feld bereits gilt —
    also die geglättete Orografie des CAMS-Gitters, NICHT Meeresniveau. Sonst
    kommt der Höheneffekt doppelt hinein, weil CAMS ihn zum Teil schon enthält.
    """
    return 1.0 + gradient * (height_m - reference_height_m) / 1000.0


def snow_factor(snow_cover_percent: float,
                bonus: float = SNOW_ALBEDO_BONUS) -> float:
    """
    Albedo-Zuschlag durch Schneedecke.

    ICON liefert SNOWC in Prozent. Frischer Schnee reflektiert UV stark; im
    Hochgebirge sind Zuschläge in dieser Grössenordnung real.
    """
    c = min(100.0, max(0.0, snow_cover_percent)) / 100.0
    return 1.0 + bonus * c


def uv_index(uvbed_clearsky: float, height_m: float, reference_height_m: float,
             rel_sunshine: float, snow_cover_percent: float,
             use_snow: bool = True) -> float:
    """
    UV-Index für eine Rasterzelle.

    :param uvbed_clearsky: CAMS uvbedcs in W/m²
    :param height_m: Geländehöhe aus ICON HSURF
    :param reference_height_m: geglättete Höhe auf CAMS-Skala
    :param rel_sunshine: relative Sonnenscheindauer 0..1
    :param snow_cover_percent: ICON SNOWC in Prozent
    """
    uvi = UVBED_TO_UVI * uvbed_clearsky
    uvi *= altitude_factor(height_m, reference_height_m)
    uvi *= cloud_modification_factor(rel_sunshine)

    if use_snow:
        uvi *= snow_factor(snow_cover_percent)

    return max(0.0, uvi)


def encode_byte(uvi: float) -> int:
    """
    UV-Index als ein Byte: UVI · 10, also 0.0 bis 25.5.

    Werte darüber gibt es auf der Erde nicht (Rekord ~43 in den Anden liegt
    ausserhalb jeder Schweizer Relevanz), die Sättigung ist damit unkritisch.
    """
    return min(255, max(0, int(round(uvi * 10.0))))
