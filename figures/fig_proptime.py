# -*- coding: utf-8 -*-
"""fig_proptime.py -- Proper-time / affine functionals vs rapidity.

Reproduces the proper-time figure of Sec. III.  For a fixed accelerator
parameter a = 10 and eta in [0, 3] we plot, using the EXACT CLOSED FORMS
(Proposition "Closed forms", verified in verify_sec3_closedforms.py):

    t(eta)        = sinh(eta)                                  (coordinate time)
    tau1(eta,a)   = asinh(sqrt(1+a^2) sinh eta) / sqrt(1+a^2)  (proper-time-like)
    tau2(eta,a)   = (a/sqrt(1+a^2)) [ acosh( sqrt(1+a^2)/a cosh eta )
                                      - acosh( sqrt(1+a^2)/a ) ] (affine-like)
    tau_ph(eta)   = eta                                        (photon, dashed)

Any interior crossing between curves is located numerically (sign change +
mpmath.findroot refinement) and annotated.  No crossing is asserted that is
not actually present.

Outputs (in this directory): fig_proptime.pdf, fig_proptime.png (300 dpi).

Tested with Python 3.13, numpy 2.4, matplotlib 3.10, mpmath 1.3.
"""
import os
import numpy as np
import mpmath as mp
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# --- parameters --------------------------------------------------------------
A = 10.0
ETA_MIN, ETA_MAX = 0.0, 3.0
S1 = np.sqrt(1.0 + A * A)

# --- closed-form curves (vectorized) -----------------------------------------
def t_coord(eta):
    return np.sinh(eta)

def tau1(eta):
    return np.arcsinh(S1 * np.sinh(eta)) / S1

def tau2(eta):
    return (A / S1) * (np.arccosh(S1 / A * np.cosh(eta)) - np.arccosh(S1 / A))

def tau_ph(eta):
    return eta

# --- numeric crossing finder (honest: only real interior crossings) ----------
mp.mp.dps = 30
_S1 = mp.sqrt(1 + mp.mpf(A) ** 2)
_mp = {
    "t": lambda e: mp.sinh(e),
    "tau1": lambda e: mp.asinh(_S1 * mp.sinh(e)) / _S1,
    "tau2": lambda e: (mp.mpf(A) / _S1)
    * (mp.acosh(_S1 / mp.mpf(A) * mp.cosh(e)) - mp.acosh(_S1 / mp.mpf(A))),
    "tau_ph": lambda e: e,
}

def find_crossings():
    """Return list of (name_i, name_j, eta*, value) for interior crossings."""
    grid = np.linspace(ETA_MIN + 1e-9, ETA_MAX, 40000)
    fns = {"t": t_coord, "tau1": tau1, "tau2": tau2, "tau_ph": tau_ph}
    names = list(fns)
    out = []
    for ii in range(len(names)):
        for jj in range(ii + 1, len(names)):
            ni, nj = names[ii], names[jj]
            d = fns[ni](grid) - fns[nj](grid)
            idx = np.where(np.diff(np.sign(d)) != 0)[0]
            for k in idx:
                if grid[k] <= 1e-3:  # skip the trivial common origin
                    continue
                g = lambda e: _mp[ni](e) - _mp[nj](e)
                try:
                    root = mp.findroot(g, mp.mpf(grid[k]))
                except Exception:
                    continue
                if ETA_MIN < float(root) < ETA_MAX:
                    out.append((ni, nj, float(root), float(_mp[ni](root))))
    return out

# --- figure ------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.6,
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "font.family": "serif",
})

fig, ax = plt.subplots(figsize=(3.4, 2.9))
eta = np.linspace(ETA_MIN, ETA_MAX, 800)

ax.plot(eta, t_coord(eta), color="#1b1b1b", label=r"$t=\sinh\eta$ (coord. time)")
ax.plot(eta, tau1(eta), color="#0072B2",
        label=r"$\tau_1(\eta,a)$ (proper-time-like)")
ax.plot(eta, tau2(eta), color="#D55E00",
        label=r"$\tau_2(\eta,a)$ (affine-like)")
ax.plot(eta, tau_ph(eta), color="#555555", ls="--", lw=1.3,
        label=r"$\tau_{\mathrm{ph}}=\eta$ (photon)")

# annotate real crossings only
crossings = find_crossings()
for (ni, nj, ex, ey) in crossings:
    ax.plot([ex], [ey], marker="o", ms=5, mfc="none", mec="k", mew=1.2, zorder=5)
    ax.annotate(
        r"$\tau_1=\tau_2$" + "\n" + r"$\eta^\ast\!=\!{:.3f}$".format(ex),
        xy=(ex, ey), xytext=(ex + 0.35, ey + 1.15),
        fontsize=8, ha="left",
        arrowprops=dict(arrowstyle="->", lw=0.8, color="k"))

ax.set_xlabel(r"Rapidity $\eta$")
ax.set_ylabel(r"proper-time / affine functional")
ax.set_xlim(ETA_MIN, ETA_MAX)
ax.set_ylim(0, None)
ax.text(0.03, 0.97, r"$a={:.0f}$".format(A), transform=ax.transAxes,
        ha="left", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", lw=0.6))
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.90), frameon=False,
          handlelength=1.6, labelspacing=0.3)
ax.grid(True, which="major", ls=":", lw=0.5, color="0.85")
fig.tight_layout(pad=0.4)

pdf_path = os.path.join(HERE, "fig_proptime.pdf")
png_path = os.path.join(HERE, "fig_proptime.png")
fig.savefig(pdf_path)
fig.savefig(png_path, dpi=300)
plt.close(fig)

print("fig_proptime: interior crossings found = {}".format(len(crossings)))
for (ni, nj, ex, ey) in crossings:
    print("  {} x {} at eta = {:.6f}, value = {:.6f}".format(ni, nj, ex, ey))
print("wrote {}".format(pdf_path))
print("wrote {}".format(png_path))
