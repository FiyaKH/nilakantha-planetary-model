import numpy as np, warnings, io, contextlib
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricMeanEcliptic, solar_system_ephemeris
from model import JD_KALI_EPOCH, jd_julian_calendar, true_sun
from exterior import PLANETS, planet, sun_variant
warnings.filterwarnings("ignore")
with contextlib.redirect_stdout(io.StringIO()):          # reuse saturn.py's data without its printout
    import saturn as S
wrap = lambda x: (x + 180) % 360 - 180

# 1. Frame offset: Sun's J2000-ecliptic longitude minus Nilakantha's, averaged over 1499-1500
jd0 = jd_julian_calendar(1499, 1, 1)
ahar = np.arange(np.ceil(jd0 - JD_KALI_EPOCH), jd0 + 730 - JD_KALI_EPOCH, 2.0)
t = Time(JD_KALI_EPOCH + ahar, format="jd", scale="tt")
with solar_system_ephemeris.set("builtin"):
    lj = get_body("sun", t).transform_to(GeocentricMeanEcliptic(equinox="J2000")).lon.deg
c = wrap(lj - true_sun(ahar)).mean()
print(f"Frame offset (J2000 ecliptic minus Nilakantha), c.1500: {c:.2f} deg\n")

# 2. Modern apogees from JPL Table 2a (Standish & Williams), converted to his frame
T = -5.0                                                  # centuries from J2000 back to 1500
JPL = {"sun": (102.93005885, 0.31795260, 0), "jupiter": (14.27495244, 0.18199196, 180),
       "saturn": (92.86136063, 0.54179478, 180)}          # perihelion longitude, rate (deg/century), +180 for aphelion
his = {"sun": 78.0, "jupiter": 172.0, "saturn": 240.0}
apo = {}
for b, (p0, rate, add) in JPL.items():
    apo[b] = (p0 + rate * T + add - c) % 360
    print(f"{b:8s} modern apogee in his frame {apo[b]:7.2f} | his {his[b]:6.1f} | difference {wrap(his[b] - apo[b]):+.2f}")

# 3. Saturn ablation with the modern apogee
print()
for label, a, mod_e in [("his apogee 240.0", 240.0, False),
                        (f"modern apogee {apo['saturn']:.1f}", apo["saturn"], False),
                        ("modern apogee + modern e", apo["saturn"], True)]:
    PLANETS["saturn"]["apogee"] = a
    err = wrap(wrap(planet("saturn", S.ahar, modern_manda=mod_e) - sun_variant(S.ahar)) - S.elong_modern)
    amps = S.harmonic_amplitudes(S.year, err)
    print(f"{label:28s} RMS {np.sqrt((err**2).mean()):.3f} | 29.46-yr amplitude {amps[3]:.3f}")
PLANETS["saturn"]["apogee"] = 240.0