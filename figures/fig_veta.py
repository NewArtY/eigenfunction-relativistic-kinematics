#!/usr/bin/env python3
"""
FIG 1 (key figure): Longitudinal velocity vs rapidity.

Massive particle:  v = tanh(eta)   (solid)
Photon:            v = 1            (solid horizontal)

Honest statement: the curves do NOT coincide on any finite interval
(tanh never equals 1). They agree ASYMPTOTICALLY: tanh(eta) -> +/-1 as
|eta| -> infinity. We shade where |1 - tanh(eta)| < epsilon (epsilon=0.01)
and mark eta_* = artanh(1 - epsilon).

Publication target: Physical Review D, single-column (~3.4 in wide).
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- rcParams
mpl.rcParams.update({
    "text.usetex": False,          # do NOT depend on a LaTeX install
    "mathtext.fontset": "cm",      # Computer-Modern-like math
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "pdf.fonttype": 42,            # embed TrueType (editable / PRD-safe)
    "ps.fonttype": 42,
})

# colorblind-friendly (Wong palette)
C_MASSIVE = "#0072B2"   # blue
C_PHOTON  = "#D55E00"   # vermillion
C_ETASTAR = "#333333"   # near-black
C_SHADE   = "#009E73"   # bluish green

# ---------------------------------------------------------------- data
epsilon = 0.01
eta_star = np.arctanh(1.0 - epsilon)     # ~ 2.6467
eta = np.linspace(-4.0, 4.0, 2001)
v = np.tanh(eta)

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(3.4, 2.7))

# shaded asymptotic-agreement bands (|1 - |tanh|| < eps  -> |eta| > eta_star)
ax.axvspan(eta_star, 4.0, color=C_SHADE, alpha=0.14, lw=0, zorder=0)
ax.axvspan(-4.0, -eta_star, color=C_SHADE, alpha=0.14, lw=0, zorder=0)

# photon line v = 1  (and the mirror v = -1 asymptote, light)
ax.axhline(1.0, color=C_PHOTON, lw=1.4, zorder=3)
ax.axhline(-1.0, color=C_PHOTON, lw=0.8, ls=(0, (4, 3)), alpha=0.55, zorder=2)

# massive curve v = tanh(eta)
ax.plot(eta, v, color=C_MASSIVE, lw=1.6, zorder=4)

# eta_* vertical dashed markers
for xs in (eta_star, -eta_star):
    ax.axvline(xs, color=C_ETASTAR, lw=0.8, ls=(0, (5, 3)), alpha=0.8, zorder=3)

# zero axes (faint)
ax.axhline(0.0, color="0.75", lw=0.5, zorder=1)
ax.axvline(0.0, color="0.75", lw=0.5, zorder=1)

# ---------------------------------------------------------------- annotations
# photon line label
ax.annotate(r"photon $v=1$", xy=(-3.5, 1.0), xytext=(-3.85, 1.075),
            color=C_PHOTON, fontsize=8.5, ha="left", va="bottom")

# massive curve label
ax.annotate(r"$v=\tanh\eta$", xy=(1.05, np.tanh(1.05)),
            xytext=(0.35, 0.30), color=C_MASSIVE, fontsize=9,
            ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=C_MASSIVE,
                            lw=0.7, shrinkA=0, shrinkB=2))

# asymptote v -> 1 label with arrow to the curve near eta_star
ax.annotate(r"$v\to 1$", xy=(3.15, np.tanh(3.15)),
            xytext=(1.7, 0.55), fontsize=8.5, color="0.25",
            ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color="0.35", lw=0.7,
                            connectionstyle="arc3,rad=-0.25",
                            shrinkA=1, shrinkB=1))

# eta_* label
ax.annotate(r"$\eta_\ast=\mathrm{artanh}(1-\epsilon)\approx 2.65$",
            xy=(eta_star, -0.72), xytext=(eta_star - 0.2, -0.72),
            fontsize=7.8, color=C_ETASTAR, ha="right", va="center",
            arrowprops=dict(arrowstyle="-", color=C_ETASTAR, lw=0.6,
                            shrinkA=0, shrinkB=0))

# agreement-band label
ax.text(3.3, 0.02, "agreement\nto within 1%", fontsize=7.6,
        color=C_SHADE, ha="center", va="center", linespacing=1.05,
        fontweight="bold")

# ---------------------------------------------------------------- axes cosmetics
ax.set_xlim(-4.0, 4.0)
ax.set_ylim(-1.15, 1.15)
ax.set_xlabel(r"Rapidity $\eta$")
ax.set_ylabel(r"Velocity $v$")

ax.set_xticks(np.arange(-4, 5, 1))
ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

# clean spines
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_position(("outward", 2))
ax.spines["bottom"].set_position(("outward", 2))

fig.tight_layout(pad=0.4)

# ---------------------------------------------------------------- save
out = ("C:/Users/artio/Documents/Scientific work/"
       "Физика/Статьи/2026.07/PRD/code/figures/fig_veta")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=300)
plt.close(fig)

print("fig_veta: eta_star =", eta_star)
print("fig_veta: 1 - tanh(eta_star) =", 1.0 - np.tanh(eta_star))
print("fig_veta: saved PDF + PNG")
