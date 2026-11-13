"""Nilakantha's planetary model (Tantrasangraha, c. 1500 CE) for the Sun and Jupiter.

Sources:
  [S98] K.V. Sarma & V.S. Narasimhan, Tantrasamgraha, Suppl. Indian J. Hist. Sci. 33(1), 1998, pp. S27-S31, S46-S48
  [RSS11] Ramasubramanian, Srinivas & Sriram, Appendix F of Tantrasangraha (Springer 2011), eqs. F.16-F.27
"""
import math
import numpy as np

# ---- Parameters from the text --------------------------------------------
CIVIL_DAYS = 1_577_917_200          # civil days in a mahayuga   [S98, v.1.18-22]
REV = {"sun": 4_320_000,            # revolutions per mahayuga  [S98, v.1.15-16]
       "jupiter": 364_180}
DHRUVA = {"sun": 0.0,               # position at Kali epoch (deg); Jupiter: +12 deg 10'  [S98, v.1.35-38]
          "jupiter": 12 + 10 / 60}
APOGEE = {"sun": 78.0,              # mandocca (deg)            [S98, v.1.40]
          "jupiter": 172.0}

# Kali epoch: mean sunrise (6h local mean time) at Lanka, Ujjain meridian (75.77 E), 18 Feb 3102 BCE (Julian)
UJJAIN_LON = 75.77

def jd_julian_calendar(year, month, day):
    """Julian day at 0h UT of a Julian-calendar date (astronomical year numbering; Meeus, ch. 7)."""
    if month <= 2:
        year, month = year - 1, month + 12
    return math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day - 1524.5

JD_KALI_EPOCH = jd_julian_calendar(-3101, 2, 18) + (6.0 - UJJAIN_LON / 15.0) / 24.0   # UT of sunrise at Lanka

def sind(x): return np.sin(np.radians(x))
def cosd(x): return np.cos(np.radians(x))
def asind(x): return np.degrees(np.arcsin(x))

def mean_longitude(body, ahargana):
    """Mean sidereal longitude (deg) after `ahargana` civil days since the Kali epoch [S98, v.1.26-28]."""
    revs = ahargana * REV[body] / CIVIL_DAYS
    return (360.0 * (revs % 1.0) + DHRUVA[body]) % 360.0

# ---- Sun: manda correction only, r0/R = 3/80  [RSS11, Table F.2] ------------
def true_sun(ahargana):
    m = mean_longitude("sun", ahargana)
    return (m - asind((3 / 80) * sind(m - APOGEE["sun"]))) % 360.0

# ---- Jupiter: four-step process  [S98, v.2.60-68; RSS11, eq. F.27] ------------
def manda_ratio_jupiter(arg):  return (7 + np.abs(sind(arg))) / 82
def sighra_ratio_jupiter(arg): return (16 - np.abs(sind(arg))) / 80

def sighra_karna(rho, d):
    """Normalised sighra hypotenuse (R = 1); d = planet - sighrocca."""
    return np.sqrt((rho * sind(d)) ** 2 + (1 + rho * cosd(d)) ** 2)

def true_jupiter(ahargana, return_steps=False):
    t0 = mean_longitude("jupiter", ahargana)
    ts = mean_longitude("sun", ahargana)                 # sighrocca = mean Sun
    tm = APOGEE["jupiter"]

    # Step 1: half manda correction
    t1 = t0 + 0.5 * asind(-manda_ratio_jupiter(t0 - tm) * sind(t0 - tm))
    # Step 2: half sighra correction
    rho = sighra_ratio_jupiter(t1 - ts)
    t2 = t1 + 0.5 * asind(rho * sind(ts - t1) / sighra_karna(rho, t1 - ts))
    # Step 3: full manda correction (anomaly from t2), applied to the mean planet
    tms = t0 + asind(-manda_ratio_jupiter(t2 - tm) * sind(t2 - tm))
    # Step 4: full sighra correction applied to the manda-sphuta
    rho = sighra_ratio_jupiter(tms - ts)
    true = tms + asind(rho * sind(ts - tms) / sighra_karna(rho, tms - ts))

    true = true % 360.0
    if return_steps:
        return dict(mean=t0, mean_sun=ts, t1=t1 % 360, t2=t2 % 360, manda_sphuta=tms % 360, true=true)
    return true
