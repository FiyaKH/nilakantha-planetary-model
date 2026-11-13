"""Nilakantha's four-step model for the exterior planets (Jupiter, Saturn), generalised.
Sources: [S98] Sarma & Narasimhan 1998, vv. 1.15-16, 1.35-40, 2.60-68; [RSS11] Appendix F, eq. F.27, Tables F.1-F.3.
"""
import numpy as np
from model import CIVIL_DAYS, sind, cosd, asind, true_sun, mean_longitude as _ml
from modern_elements import elements

PLANETS = {
    "jupiter": dict(rev=364_180, dhruva=12 + 10 / 60, apogee=172.0,
                    manda=lambda a: (7 + np.abs(sind(a))) / 82,
                    sighra=lambda a: (16 - np.abs(sind(a))) / 80,
                    e_modern=elements("jupiter")["e"], a_modern=elements("jupiter")["a"], period=11.862),
    "saturn":  dict(rev=146_612, dhruva=347 + 20 / 60, apogee=240.0,     # 11 rasi 17 deg 20'
                    manda=lambda a: 39 / 320 + 0 * a,                    # (1/8)(1 - 1/40), constant
                    sighra=lambda a: (9 - np.abs(sind(a))) / 80,
                    e_modern=elements("saturn")["e"], a_modern=elements("saturn")["a"], period=29.457),
}
E_EARTH = elements("earth")["e"]

def mean_planet(name, ahar):
    revs = ahar * PLANETS[name]["rev"] / CIVIL_DAYS
    return (360.0 * (revs % 1.0) + PLANETS[name]["dhruva"]) % 360.0

def sun_variant(ahar, modern_e=False, ms=None):
    m = _ml("sun", ahar) if ms is None else ms
    return (m - asind((2 * E_EARTH if modern_e else 3 / 80) * sind(m - 78.0))) % 360

def planet(name, ahar, true_sun_sighrocca=False, earth_dist=False, planet_dist=False,
           modern_manda=False, constant_sighra=False, process="four", apogee=None, sun_modern_e=False, t0=None, ms=None):
    P = PLANETS[name]
    t0 = mean_planet(name, ahar) if t0 is None else t0
    ms = _ml("sun", ahar) if ms is None else ms
    tm = P["apogee"] if apogee is None else apogee
    ts = sun_variant(ahar, sun_modern_e, ms) if true_sun_sighrocca else ms
    mr = (lambda a: 2 * P["e_modern"] + 0 * a) if modern_manda else P["manda"]
    sr = (lambda a: 1 / P["a_modern"] + 0 * a) if constant_sighra else P["sighra"]
    r_earth = 1 + E_EARTH * cosd(ms - 78.0) if earth_dist else 1.0

    def sighra_corr(theta, K=1.0):
        r = sr(theta - ts) * r_earth / K
        k = np.sqrt((r * sind(theta - ts)) ** 2 + (1 + r * cosd(theta - ts)) ** 2)
        return asind(r * sind(ts - theta) / k)

    if process == "four":
        t1 = t0 + 0.5 * asind(-mr(t0 - tm) * sind(t0 - tm))
        t2 = t1 + 0.5 * sighra_corr(t1)
        tms = t0 + asind(-mr(t2 - tm) * sind(t2 - tm))
    else:
        tms = t0 + asind(-mr(t0 - tm) * sind(t0 - tm))
    K = 1 + P["e_modern"] * cosd(tms - tm) if planet_dist else 1.0
    return (tms + sighra_corr(tms, K)) % 360