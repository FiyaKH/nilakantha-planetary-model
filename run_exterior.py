"""Ablation study for Jupiter and Saturn with all modern values from modern_elements.py (JPL Table 2a)."""
import numpy as np, warnings
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricTrueEcliptic, solar_system_ephemeris
from model import JD_KALI_EPOCH, jd_julian_calendar
from exterior import planet, sun_variant, PLANETS
from modern_elements import elements
warnings.filterwarnings("ignore")
wrap = lambda x: (x + 180) % 360 - 180
FRAME_OFFSET = 22.56   # J2000 ecliptic minus Nilakantha frame, c.1500 (measured in saturn_apogee.py)

jd0, jd1 = jd_julian_calendar(1400, 1, 1), jd_julian_calendar(1700, 1, 1)
ahar = np.arange(np.ceil(jd0 - JD_KALI_EPOCH), jd1 - JD_KALI_EPOCH, 5.0)
jd = JD_KALI_EPOCH + ahar
year = 1400 + (jd - jd0) / 365.25
t = Time(jd, format="jd", scale="tt")
with solar_system_ephemeris.set("builtin"):
    lon = {b: get_body(b, t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg
           for b in ("sun", "jupiter", "saturn")}

def periods(name):
    Po = PLANETS[name]["period"]; Ps = 1 / (1 - 1 / Po)
    return {"1.0 yr": 1.0, "synodic": Ps, "sideband": 1 / (1 / Ps - 1 / Po), "orbit": Po, "orbit/2": Po / 2, "orbit/3": Po / 3}

def harmonic_amplitudes(y, err, P):
    cols = [np.ones_like(y), y - y.mean()]
    for p in P.values():
        cols += [np.sin(2 * np.pi * y / p), np.cos(2 * np.pi * y / p)]
    coef, *_ = np.linalg.lstsq(np.stack(cols, 1), err, rcond=None)
    return [np.hypot(coef[2 + 2 * i], coef[3 + 2 * i]) for i in range(len(P))]

for name in ("jupiter", "saturn"):
    apo_mod = (elements(name)["perihelion_J2000"] + 180 - FRAME_OFFSET) % 360
    elong_modern = wrap(lon[name] - lon["sun"])
    ALL = dict(true_sun_sighrocca=True, earth_dist=True, modern_manda=True, constant_sighra=True,
               apogee=apo_mod, sun_modern_e=True)
    VARIANTS = [
        ("baseline (as in text)",           {}, {}),
        ("his Sun -> modern Sun e",         {}, {"modern_e": True}),
        ("Earth ellipse (true Sun + dist)", {"true_sun_sighrocca": True, "earth_dist": True}, {}),
        ("constant modern sighra radius",   {"constant_sighra": True}, {}),
        ("modern manda (2e)",               {"modern_manda": True}, {}),
        (f"modern apogee ({apo_mod:.1f})",  {"apogee": apo_mod}, {}),
        ("two-step",                        {"process": "two"}, {}),
        ("two-step + planet distance",      {"process": "two", "planet_dist": True}, {}),
        ("all of the above modern",         ALL, {"modern_e": True}),
    ]
    P = periods(name)
    print(f"\n{name.upper()}  (his apogee {PLANETS[name]['apogee']:.0f}, modern {apo_mod:.1f}; "
          f"periods: synodic {P['synodic']:.3f}, sideband {P['sideband']:.3f}, orbit {P['orbit']:.2f} yr)")
    print(f"{'variant':34s} {'RMS':>6s} " + " ".join(f"{k:>9s}" for k in P))
    for label, pkw, skw in VARIANTS:
        err = wrap(wrap(planet(name, ahar, **pkw) - sun_variant(ahar, **skw)) - elong_modern)
        print(f"{label:34s} {np.sqrt((err**2).mean()):6.3f} " + " ".join(f"{a:9.3f}" for a in harmonic_amplitudes(year, err, P)))