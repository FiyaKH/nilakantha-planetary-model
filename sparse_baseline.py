"""Step 2b: can sparse regression find the error structure automatically, without our physical hypotheses?
Fits the ORIGINAL errors (model exactly as in the text) with a Fourier library in the model's own angles,
and compares the selected terms with the ablation results. For Venus, a physics-informed variant adds the
same terms divided by the Earth-Venus distance (errors are magnified when Venus is close to Earth).
"""
import numpy as np, io, contextlib, warnings
from sklearn.linear_model import LassoCV, LinearRegression
warnings.filterwarnings("ignore")
with contextlib.redirect_stdout(io.StringIO()):
    import run_exterior as X
    import venus as V
from exterior import planet, sun_variant, mean_planet, PLANETS
from model import mean_longitude, JD_KALI_EPOCH

Tc = (JD_KALI_EPOCH + X.ahar - 2451545.0) / 36525.0
Tc = Tc - Tc.mean()
ms = mean_longitude("sun", X.ahar)

def library(t0, apo):
    M = np.radians(t0 - (apo - 180)); E = np.radians(ms - 78.0); D = np.radians(t0 - ms)
    A = {"D": D, "M": M, "E": E, "2D": 2 * D}; args = {}
    for k in (1, 2, 3):
        args[f"{k}M"], args[f"{k}E"], args[f"{k}D"] = k * M, k * E, k * D
    for a, b in [("D", "M"), ("D", "E"), ("M", "E"), ("2D", "M"), ("2D", "E")]:
        args[f"{a}+{b}"], args[f"{a}-{b}"] = A[a] + A[b], A[a] - A[b]
    L = np.stack([f(g) for g in args.values() for f in (np.sin, np.cos)], 1)
    names = [f"{fn}({g})" for g in args for fn in ("sin", "cos")]
    return L, names, list(args), D

def sparse_fit(err, L, names, groups, label, top=6):
    las = LassoCV(cv=5, random_state=0).fit(np.column_stack([L, Tc]), err)
    keep = np.flatnonzero(np.abs(las.coef_[:-1]) > 1e-4)
    ols = LinearRegression().fit(np.column_stack([L[:, keep], Tc]), err)
    resid = err - ols.predict(np.column_stack([L[:, keep], Tc]))
    coef = dict(zip([names[i] for i in keep], ols.coef_[:-1]))
    amps = sorted(((g, np.hypot(coef.get(f"sin({g})", 0), coef.get(f"cos({g})", 0))) for g in groups), key=lambda kv: -kv[1])
    print(f"{label}: RMS {np.sqrt((err**2).mean()):.3f} -> {np.sqrt((resid**2).mean()):.3f} deg with {len(keep)} terms")
    print("   largest: " + ", ".join(f"{g} {a:.3f}" for g, a in amps[:top] if a > 0))

for name in ("jupiter", "saturn"):
    err = X.wrap(X.wrap(planet(name, X.ahar) - sun_variant(X.ahar)) - X.wrap(X.lon[name] - X.lon["sun"]))
    L, names, groups, _ = library(mean_planet(name, X.ahar), PLANETS[name]["apogee"])
    sparse_fit(err, L, names, groups, f"{name} (as in text)")

err = V.wrap(V.wrap(V.venus(V.ahar) - V.sun(V.ahar)) - V.elong_modern)
L, names, groups, D = library(V.mean_venus(V.ahar), V.APOGEE_VENUS)
sparse_fit(err, L, names, groups, "venus (as in text), plain library")
dist = np.sqrt(1 + V.RATIO_VENUS ** 2 + 2 * V.RATIO_VENUS * np.cos(D))          # Earth-Venus distance, mean orbits
# Same terms, each divided by the Earth-Venus distance (amplitudes are then in degrees per unit 1/distance)
sparse_fit(err, L / dist[:, None], names, groups, "venus (as in text), terms / Earth-Venus distance")