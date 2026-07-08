# -*- coding: utf-8 -*-
"""fig_finterp.py -- The interpolating velocity functional f(eta,a).

Plots the eigenfunction-derived velocity functional of Sec. III

    f(eta,a) = a*tanh(eta) / sqrt(1 + a^2 tanh^2 eta),   a > 0,

for a in {0.5, 1, 2, 5, 10, 50} over eta in [-3, 3], together with the
massive reference v = tanh(eta) and the photon plateau v = +/- 1.  As a grows,
f approaches the null plateau |v| = 1 monotonically (ceiling velocity
v_inf(a) = a/sqrt(1+a^2), verified in verify_sec3_closedforms.py).

Outputs (in this directory): fig_finterp.pdf, fig_finterp.png (300 dpi).

Tested with Python 3.13, numpy 2.4, matplotlib 3.10.
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable

HERE = os.path.dirname(os.path.abspath(__file__))

A_VALUES = [0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
ETA_MIN, ETA_MAX = -3.0, 3.0

def f(eta, a):
    t = np.tanh(eta)
    return a * t / np.sqrt(1.0 + a * a * t * t)

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.6,
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "font.family": "serif",
})

fig, ax = plt.subplots(figsize=(3.4, 3.0))
eta = np.linspace(ETA_MIN, ETA_MAX, 800)

# colour-code by a on a log scale
norm = LogNorm(vmin=min(A_VALUES), vmax=max(A_VALUES))
cmap = mpl.colormaps["viridis"]

for a in A_VALUES:
    ax.plot(eta, f(eta, a), color=cmap(norm(a)), lw=1.7,
            label=r"$a={:g}$".format(a))

# massive reference v = tanh(eta)
ax.plot(eta, np.tanh(eta), color="k", ls="-.", lw=1.2,
        label=r"$v=\tanh\eta$")

# photon plateaus v = +/- 1
ax.axhline(1.0, color="0.35", ls="--", lw=1.1)
ax.axhline(-1.0, color="0.35", ls="--", lw=1.1)
ax.text(ETA_MAX - 0.05, 1.0, r"$v=+1$ (photon)", ha="right", va="bottom",
        fontsize=7.5, color="0.35")
ax.text(ETA_MIN + 0.05, -1.0, r"$v=-1$ (photon)", ha="left", va="top",
        fontsize=7.5, color="0.35")

ax.set_xlabel(r"Rapidity $\eta$")
ax.set_ylabel(r"velocity functional $f(\eta,a)$")
ax.set_xlim(ETA_MIN, ETA_MAX)
ax.set_ylim(-1.18, 1.18)
ax.grid(True, ls=":", lw=0.5, color="0.85")

# colourbar for a, plus a compact legend for the reference curves
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.05, ticks=A_VALUES)
cbar.ax.set_yticklabels([r"${:g}$".format(a) for a in A_VALUES])
cbar.set_label(r"$a$ (log scale)", fontsize=8)
cbar.ax.tick_params(labelsize=7.5)

ax.legend(loc="lower right", frameon=False, ncol=1, handlelength=1.5,
          labelspacing=0.25, fontsize=7.2)
fig.tight_layout(pad=0.4)

pdf_path = os.path.join(HERE, "fig_finterp.pdf")
png_path = os.path.join(HERE, "fig_finterp.png")
fig.savefig(pdf_path)
fig.savefig(png_path, dpi=300)
plt.close(fig)

# quick monotone-approach report
print("fig_finterp: ceiling velocity v_inf(a) = a/sqrt(1+a^2):")
for a in A_VALUES:
    print("  a={:>5g}  v_inf={:.6f}".format(a, a / np.sqrt(1 + a * a)))
print("wrote {}".format(pdf_path))
print("wrote {}".format(png_path))
