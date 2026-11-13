import numpy as np, warnings
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricTrueEcliptic, solar_system_ephemeris
from model import JD_KALI_EPOCH, jd_julian_calendar
from exterior import planet, sun_variant, PLANETS
warnings.filterwarnings("ignore")
wrap = lambda x: (x + 180) % 360 - 180

# Same dates as Jupiter: mean sunrise at Lanka, every 5 days, 1400-1700 CE
jd0, jd1 = jd_julian_calendar(1400, 1, 1), jd_julian_calendar(1700, 1, 1)
ahar = np.arange(np.ceil(jd0 - JD_KALI_EPOCH), jd1 - JD_KALI_EPOCH, 5.0)
jd = JD_KALI_EPOCH + ahar
year = 1400 + (jd - jd0) / 365.25
t = Time(jd, format="jd", scale="tt")
with solar_system_ephemeris.set("builtin"):
    lon = {b: get_body(b, t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg for b in ("sun", "saturn")}
elong_modern = wrap(lon["saturn"] - lon["sun"])

Ps, Psyn = PLANETS["saturn"]["period"], 1 / (1 - 1 / PLANETS["saturn"]["period"])
PERIODS = {"1.0 yr": 1.0, f"syn {Psyn:.3f}": Psyn, f"{1/(1/Psyn - 1/Ps):.3f} yr": 1 / (1 / Psyn - 1 / Ps),
           f"orbit {Ps:.2f}": Ps, f"orbit/2": Ps / 2, f"orbit/3": Ps / 3}
def harmonic_amplitudes(y, err):
    cols = [np.ones_like(y), y - y.mean()]
    for P in PERIODS.values():
        cols += [np.sin(2 * np.pi * y / P), np.cos(2 * np.pi * y / P)]
    coef, *_ = np.linalg.lstsq(np.stack(cols, 1), err, rcond=None)
    return [np.hypot(coef[2 + 2 * i], coef[3 + 2 * i]) for i in range(len(PERIODS))]

VARIANTS = [
    ("baseline (as in text)",          {}, {}),
    ("his Sun -> modern Sun e",         {}, {"modern_e": True}),
    ("Earth ellipse (true Sun + dist)", {"true_sun_sighrocca": True, "earth_dist": True}, {}),
    ("constant modern sighra radius",   {"constant_sighra": True}, {}),
    ("modern manda (2e)",               {"modern_manda": True}, {}),
    ("two-step",                        {"process": "two"}, {}),
    ("two-step + Saturn distance",      {"process": "two", "planet_dist": True}, {}),
]
print(f"{'Saturn variant':34s} {'RMS':>6s} " + " ".join(f"{k:>12s}" for k in PERIODS))
for name, pkw, skw in VARIANTS:
    err = wrap(wrap(planet("saturn", ahar, **pkw) - sun_variant(ahar, **skw)) - elong_modern)
    print(f"{name:34s} {np.sqrt((err**2).mean()):6.3f} " + " ".join(f"{a:12.3f}" for a in harmonic_amplitudes(year, err)))