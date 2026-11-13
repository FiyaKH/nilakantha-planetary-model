import numpy as np, warnings
from astropy.time import Time
from astropy.coordinates import get_body, GeocentricTrueEcliptic, solar_system_ephemeris
import astropy.units as u
from model import JD_KALI_EPOCH, true_sun, true_jupiter, jd_julian_calendar
warnings.filterwarnings("ignore")

# Dates: mean sunrise at Lanka, every 5 days, 1400-1700 CE (Julian calendar years)
jd0 = jd_julian_calendar(1400, 1, 1); jd1 = jd_julian_calendar(1700, 1, 1)
ahar = np.arange(np.ceil(jd0 - JD_KALI_EPOCH), jd1 - JD_KALI_EPOCH, 5.0)   # integer ahargana
jd = JD_KALI_EPOCH + ahar
year = 1400 + (jd - jd0) / 365.25

# Modern positions (built-in ERFA ephemeris), apparent geocentric ecliptic longitude of date
t = Time(jd, format="jd", scale="tt")          # Delta-T (~few min) ignored: < 0.001 deg for Jupiter
with solar_system_ephemeris.set("builtin"):
    lon = {b: get_body(b, t).transform_to(GeocentricTrueEcliptic(equinox=t)).lon.deg for b in ("sun", "jupiter")}

wrap = lambda x: (x + 180) % 360 - 180
elong_model  = wrap(true_jupiter(ahar) - true_sun(ahar))
elong_modern = wrap(lon["jupiter"] - lon["sun"])
err = wrap(elong_model - elong_modern)

np.savez("jupiter_1400_1700.npz", year=year, ahar=ahar, jd=jd, err=err,
         elong_model=elong_model, elong_modern=elong_modern,
         jup_model=true_jupiter(ahar), sun_model=true_sun(ahar),
         jup_modern=lon["jupiter"], sun_modern=lon["sun"])
print(f"N = {len(err)} dates")
print(f"Elongation error (model - modern), deg: mean {err.mean():+.3f}, std {err.std():.3f}, "
      f"max |err| {np.abs(err).max():.3f}, RMS {np.sqrt((err**2).mean()):.3f}")
for y0 in (1400, 1500, 1600, 1699):
    m = (year >= y0) & (year < y0 + 1)
    print(f"  {y0}: mean err {err[m].mean():+.3f} deg")
