"""Figures for the README: error time series (as in text vs after all fixes) and RMS reduction by stage."""
import os, io, contextlib, warnings
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
with contextlib.redirect_stdout(io.StringIO()):
    import sparse_residual as S
    import venus as V
from exterior import planet, sun_variant
from modern_elements import elements, mean_longitude_J2000
from model import sind, cosd, asind
X, wrap, T, OFF = S.X, S.wrap, S.T, S.X.FRAME_OFFSET
os.makedirs("figures", exist_ok=True)
rms = lambda e: float(np.sqrt((e ** 2).mean()))

def exterior_stages(name):
    em = wrap(X.lon[name] - X.lon["sun"])
    apo = (elements(name)["perihelion_J2000"] + 180 - OFF) % 360
    ALL = dict(true_sun_sighrocca=True, earth_dist=True, modern_manda=True, constant_sighra=True, apogee=apo, sun_modern_e=True)
    err = lambda kw, ms=None: wrap(wrap(planet(name, X.ahar, **kw) - sun_variant(X.ahar, bool(kw), ms)) - em)
    s0 = err({})
    s1 = err(ALL)
    s2 = err(dict(ALL, t0=(mean_longitude_J2000(name, T, False) - OFF) % 360, ms=S.ms_mod), S.ms_mod)
    s3 = err(dict(ALL, t0=(mean_longitude_J2000(name, T, True) - OFF) % 360, ms=S.ms_mod), S.ms_mod)
    return [s0, s1, s2, s3]

def venus_stages():
    s0 = wrap(wrap(V.venus(V.ahar) - V.sun(V.ahar)) - V.elong_modern)
    FULL = dict(true_sun_sighrocca=True, earth_dist=True, manda="modern", karna="modern",
                apogee=V.APOGEE_MODERN, sun_modern_e=True, sighra_r="constant")
    s1 = wrap(wrap(V.venus(V.ahar, **FULL) - V.sun(V.ahar, modern_e=True)) - V.elong_modern)
    s3 = S.final_residual("venus")[0]
    return [s0, s1, s3, s3]                      # no long-period terms for Venus: last two stages coincide

stages = {"Jupiter": exterior_stages("jupiter"), "Saturn": exterior_stages("saturn"), "Venus": venus_stages()}
labels = ["as in text", "+ orbit fixes\n(Sun, apogee, e, sighra,\nEarth ellipse)", "+ modern mean\nlongitudes", "+ long-period terms\n(great inequality)"]

# Figure 1: error time series, as in text vs after all fixes
fig, ax = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
for a, (p, s) in zip(ax, stages.items()):
    a.plot(X.year, s[0], lw=0.4, color="#9db4c8", label=f"as in text (RMS {rms(s[0]):.2f} deg)")
    a.plot(X.year, s[-1], lw=0.5, color="#1f4e79", label=f"after all fixes (RMS {rms(s[-1]):.2f} deg)")
    a.axhline(0, color="k", lw=0.5); a.set_ylabel(f"{p}\nerror (deg)"); a.legend(loc="upper right", fontsize=8)
ax[-1].set_xlabel("Year (CE)")
fig.suptitle("Angular distance from the Sun: Nilakantha's model minus modern ephemeris, 1400-1700 CE", fontsize=11)
plt.tight_layout(); plt.savefig("figures/error_timeseries.png", dpi=130)

# Figure 2: RMS by stage
fig, ax = plt.subplots(figsize=(9, 4.5))
w, xs = 0.25, np.arange(len(labels))
for i, (p, s) in enumerate(stages.items()):
    vals = [rms(e) for e in s]
    ax.bar(xs + (i - 1) * w, vals, w, label=p)
    for xx, v in zip(xs + (i - 1) * w, vals):
        ax.text(xx, v + 0.01, f"{v:.2f}", ha="center", fontsize=7)
ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("RMS error (deg)"); ax.set_title("Error remaining as each of Nilakantha's ingredients is replaced by its modern value")
ax.legend(); plt.tight_layout(); plt.savefig("figures/error_reduction.png", dpi=130)

for p, s in stages.items():
    print(p, " -> ".join(f"{rms(e):.3f}" for e in s))
print("Saved figures/error_timeseries.png and figures/error_reduction.png")