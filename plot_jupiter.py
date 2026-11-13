import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = np.load("jupiter_1400_1700.npz"); y, err = d["year"], d["err"]

fig, ax = plt.subplots(2, 1, figsize=(10, 7.5))
ax[0].plot(y, err, lw=0.4, color="#1f4e79")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xlabel("Year (CE)"); ax[0].set_ylabel("Error (deg)")
ax[0].set_title("Jupiter's elongation from the Sun: Nilakantha's model minus modern ephemeris")
ax[0].text(0.01, 0.04, f"RMS {np.sqrt((err**2).mean()):.2f} deg, max |error| {np.abs(err).max():.2f} deg",
           transform=ax[0].transAxes)

r = err - np.polyval(np.polyfit(y, err, 1), y)
dt = 5 / 365.25
A = np.abs(np.fft.rfft(r * np.hanning(len(r)))); f = np.fft.rfftfreq(len(r), dt)
A = A / A.max()
ax[1].semilogx(1 / f[1:], A[1:], color="#1f4e79", lw=0.8)
for P, lab, h, ha in [(11.86, "Jupiter's orbit\n11.86 yr", 0.85, "left"), (1.0, "Earth-Sun\n1 yr", 0.85, "right"),
               (1.092, "Jupiter synodic\n1.09 yr", 0.62, "left"), (11.86 / 2, "orbit/2\n5.93 yr", 0.62, "left"),
               (11.86 / 3, "orbit/3\n3.95 yr", 0.85, "left")]:
    ax[1].axvline(P, color="#c00000", ls="--", lw=0.7)
    ax[1].text(P * (1.03 if ha == "left" else 0.97), h, lab, color="#c00000", fontsize=8, va="top", ha=ha)
ax[1].set_xlim(0.25, 60); ax[1].set_xlabel("Period (years, log scale)")
ax[1].set_ylabel("Relative amplitude"); ax[1].set_title("Which cycles the error follows (FFT of detrended error)")
plt.tight_layout(); plt.savefig("jupiter_error.png", dpi=130)
print("saved")
