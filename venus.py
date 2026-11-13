"""Nilakantha's revised model for Venus (interior planet).
Sources: [S98] vv. 1.15-16 (revolutions), 1.35-38 (dhruva: +1 rasi 6 deg 13'), 1.40 (apogee 80 deg), 2.69-72;
         [RSS11] Appendix F, eq. F.28, Table F.2 (manda ratio), Table F.3 (sighra ratio), eqs. F.14a, F.15 (Madhava's karna).
"""
import numpy as np, warnings
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricTrueEcliptic, solar_system_ephemeris
from model import JD_KALI_EPOCH, jd_julian_calendar, CIVIL_DAYS, sind, cosd, asind, mean_longitude
from modern_elements import elements
warnings.filterwarnings("ignore")
wrap = lambda x: (x + 180) % 360 - 180

REV_VENUS, DHRUVA_VENUS, APOGEE_VENUS, R = 7_022_268, 36 + 13 / 60, 80.0, 3438.0
E_EARTH, E_VENUS = elements("earth")["e"], elements("venus")["e"]
RATIO_VENUS = elements("venus")["a"] / elements("earth")["a"]      # a_Venus / a_Earth

def mean_venus(ahar):
    return (360.0 * ((ahar * REV_VENUS / CIVIL_DAYS) % 1.0) + DHRUVA_VENUS) % 360.0

def sun(ahar, modern_e=False):
    m = mean_longitude("sun", ahar)
    return (m - asind((2 * E_EARTH if modern_e else 3 / 80) * sind(m - 78.0))) % 360

def venus(ahar, apogee=APOGEE_VENUS, manda="text", karna="madhava", sighra_r="text",
          true_sun_sighrocca=False, earth_dist=False, sun_modern_e=False):
    t0, ms = mean_venus(ahar), mean_longitude("sun", ahar)
    ts = sun(ahar, sun_modern_e) if true_sun_sighrocca else ms
    a = t0 - apogee                                                  # manda anomaly
    q = 1 / (14 + R * np.abs(sind(a)) / 240) if manda == "text" else 2 * E_VENUS + 0 * a
    tms = t0 + asind(-q * sind(a))                                   # manda-sphuta (true heliocentric Venus)
    if karna == "madhava":                                           # iterated manda-karna, Madhava (F.14a, F.15)
        K = 1 / (np.sqrt(1 - (q * sind(a)) ** 2) - q * cosd(a))
    elif karna == "modern":                                          # physical Sun-Venus distance, first order
        K = 1 + E_VENUS * cosd(tms - apogee)
    else:
        K = 1.0
    rho = (59 - 2 * np.abs(sind(tms - ts))) / 80 if sighra_r == "text" else RATIO_VENUS + 0 * tms
    r_earth = 1 + E_EARTH * cosd(ms - 78.0) if earth_dist else 1.0
    x = rho * K / r_earth                                            # Sun-Venus distance in units of Earth-Sun distance
    return (ts + asind(x * sind(tms - ts) / np.sqrt(sind(tms - ts) ** 2 + (cosd(tms - ts) + x) ** 2))) % 360

# Dates and modern positions (same grid as Jupiter and Saturn)
jd0, jd1 = jd_julian_calendar(1400, 1, 1), jd_julian_calendar(1700, 1, 1)
ahar = np.arange(np.ceil(jd0 - JD_KALI_EPOCH), jd1 - JD_KALI_EPOCH, 5.0)
jd = JD_KALI_EPOCH + ahar
year = 1400 + (jd - jd0) / 365.25
t = Time(jd, format="jd", scale="tt")
with solar_system_ephemeris.set("builtin"):
    lon = {b: get_body(b, t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg for b in ("sun", "venus")}
elong_modern = wrap(lon["venus"] - lon["sun"])

P_V, P_SYN = 224.70080 / 365.25636, 1 / (365.25636 / 224.70080 - 1)   # Venus orbit; synodic period (years)
PERIODS = {"Earth 1.0": 1.0, f"Venus orb {P_V:.3f}": P_V, f"synodic {P_SYN:.3f}": P_SYN,
           f"{1 / (1 / P_V + 1):.3f} yr": 1 / (1 / P_V + 1), "synodic/2": P_SYN / 2}
def harmonic_amplitudes(y, err):
    cols = [np.ones_like(y), y - y.mean()]
    for P in PERIODS.values():
        cols += [np.sin(2 * np.pi * y / P), np.cos(2 * np.pi * y / P)]
    coef, *_ = np.linalg.lstsq(np.stack(cols, 1), err, rcond=None)
    return [np.hypot(coef[2 + 2 * i], coef[3 + 2 * i]) for i in range(len(PERIODS))]

FRAME_OFFSET = 22.56          # J2000 ecliptic minus Nilakantha frame, c.1500 (measured in saturn_apogee.py)
APOGEE_MODERN = (elements("venus")["perihelion_J2000"] + 180 - FRAME_OFFSET) % 360

if __name__ == "__main__":
    FULL = {"true_sun_sighrocca": True, "earth_dist": True, "manda": "modern", "karna": "modern", "apogee": APOGEE_MODERN}
    VARIANTS = [
        ("baseline (as in text)",            {}, {}),
        ("his Sun -> modern Sun e",          {}, {"modern_e": True}),
        ("Earth ellipse (true Sun + dist)",  {"true_sun_sighrocca": True, "earth_dist": True}, {}),
        ("constant modern sighra radius",    {"sighra_r": "constant"}, {}),
        ("modern manda (2e = 0.014)",        {"manda": "modern"}, {}),
        ("no manda-karna (K = 1)",           {"karna": "none"}, {}),
        ("modern manda + physical distance", {"manda": "modern", "karna": "modern"}, {}),
        ("Earth ellipse + modern Venus orbit", FULL, {}),
        ("  ... + modern Sun e everywhere", {**FULL, "sun_modern_e": True}, {"modern_e": True}),
    ]
    print(f"Modern Venus apogee in his frame: {APOGEE_MODERN:.1f} deg (his value: {APOGEE_VENUS:.0f})\n")
    print(f"{'Venus variant':34s} {'RMS':>6s} {'max':>6s} " + " ".join(f"{k:>16s}" for k in PERIODS))
    for name, vkw, skw in VARIANTS:
        err = wrap(wrap(venus(ahar, **vkw) - sun(ahar, **skw)) - elong_modern)
        print(f"{name:34s} {np.sqrt((err**2).mean()):6.3f} {np.abs(err).max():6.2f} "
              + " ".join(f"{a:16.3f}" for a in harmonic_amplitudes(year, err)))