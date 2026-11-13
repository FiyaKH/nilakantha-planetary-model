"""Step 2: sparse regression on the final residuals (after every physical fix).
A large library of sine/cosine terms in the model's own angles is fitted with an L1 (Lasso) penalty,
so only the few terms that really matter survive. Selected terms are then refitted without penalty
to report unbiased amplitudes.
"""
import numpy as np, io, contextlib, warnings
from sklearn.linear_model import LassoCV, LinearRegression
warnings.filterwarnings("ignore")
with contextlib.redirect_stdout(io.StringIO()):
    import run_exterior as X
    import venus as V
from exterior import planet, sun_variant, mean_planet
from modern_elements import elements, mean_longitude_J2000
from model import JD_KALI_EPOCH, sind, cosd, asind, mean_longitude
wrap = X.wrap

T = (JD_KALI_EPOCH + X.ahar - 2451545.0) / 36525.0
OFF = X.FRAME_OFFSET
ms_mod = (mean_longitude_J2000("earth", T) + 180 - OFF) % 360
sun_apo = (elements("earth")["perihelion_J2000"] - OFF) % 360            # Sun's apogee = Earth's perihelion direction

def final_residual(name):
    """Residual after every physical fix (all modern ingredients, incl. mean longitudes and long-period terms)."""
    if name == "venus":
        t0 = (mean_longitude_J2000("venus", T) - OFF) % 360
        apo = V.APOGEE_MODERN
        ts = ms_mod - asind(2 * V.E_EARTH * sind(ms_mod - sun_apo))
        tms = t0 + asind(-2 * V.E_VENUS * sind(t0 - apo))
        x = V.RATIO_VENUS * (1 + V.E_VENUS * cosd(tms - apo)) / (1 + V.E_EARTH * cosd(ms_mod - sun_apo))
        ven = ts + asind(x * sind(tms - ts) / np.sqrt(sind(tms - ts) ** 2 + (cosd(tms - ts) + x) ** 2))
        return wrap(wrap(ven - ts) - V.elong_modern), t0, apo
    apo = (elements(name)["perihelion_J2000"] + 180 - OFF) % 360
    t0 = (mean_longitude_J2000(name, T) - OFF) % 360
    kw = dict(true_sun_sighrocca=True, earth_dist=True, modern_manda=True, constant_sighra=True,
              apogee=apo, sun_modern_e=True, t0=t0, ms=ms_mod)
    e = wrap(wrap(planet(name, X.ahar, **kw) - sun_variant(X.ahar, True, ms_mod)) - wrap(X.lon[name] - X.lon["sun"]))
    return e, t0, apo

def library(t0, apo):
    """Candidate terms: harmonics of the planet's anomaly M, Earth's anomaly E, the synodic angle D, and their combinations."""
    M = np.radians(t0 - (apo - 180))          # planet's mean anomaly from perihelion
    E = np.radians(ms_mod - sun_apo)          # Sun's (Earth's) mean anomaly
    D = np.radians(t0 - ms_mod)               # synodic angle
    args = {}
    for k in (1, 2, 3):
        args[f"{k}M"], args[f"{k}E"], args[f"{k}D"] = k * M, k * E, k * D
    for a, b in [("D", "M"), ("D", "E"), ("M", "E"), ("2D", "M"), ("2D", "E")]:
        A = {"D": D, "M": M, "E": E, "2D": 2 * D}
        args[f"{a}+{b}"], args[f"{a}-{b}"] = A[a] + A[b], A[a] - A[b]
    cols, names = [], []
    for n, ang in args.items():
        cols += [np.sin(ang), np.cos(ang)]; names += [f"sin({n})", f"cos({n})"]
    return np.stack(cols, 1), names, list(args)

for name in ("jupiter", "saturn", "venus"):
    err, t0, apo = final_residual(name)
    Xlib, names, groups = library(t0, apo)
    Xfull = np.column_stack([Xlib, T - T.mean()])                       # include a linear drift term
    lasso = LassoCV(cv=5, random_state=0).fit(Xfull, err)
    keep = np.flatnonzero(np.abs(lasso.coef_[:-1]) > 1e-4)
    ols = LinearRegression().fit(np.column_stack([Xlib[:, keep], T - T.mean()]), err)
    resid = err - ols.predict(np.column_stack([Xlib[:, keep], T - T.mean()]))
    coef = dict(zip([names[i] for i in keep], ols.coef_[:-1]))
    amps = {}
    for g in groups:                                                       # combine sin/cos pairs into amplitudes
        s, c = coef.get(f"sin({g})", 0.0), coef.get(f"cos({g})", 0.0)
        if s or c: amps[g] = np.hypot(s, c)
    top = sorted(amps.items(), key=lambda kv: -kv[1])
    print(f"\n{name.upper()}: residual RMS {np.sqrt((err**2).mean()):.3f} deg -> {np.sqrt((resid**2).mean()):.3f} deg "
          f"after {len(keep)} selected terms (of {Xlib.shape[1]} candidates)")
    for g, a in top[:6]:
        print(f"   {g:6s} amplitude {a:.3f} deg")