#!/usr/bin/env python3
"""
FIG 4 (geometry schematic): the mass shell as a hyperbola in the
(p/m, E/m) plane.

Coordinates:  Phi2 = sinh(eta) = p/m   (horizontal)
              Phi1 = cosh(eta) = E/m   (vertical)

Upper branch of the unit hyperbola  Phi1^2 - Phi2^2 = 1  is the mass shell.
Null asymptotes  E = +/- p  form the light cone (45 deg, equal aspect).
A point P at rapidity eta0 slides along the curve under a boost eta0 -> eta0+delta.
As eta -> infinity the point approaches the upper-right asymptote: the massless limit.

Publication target: Physical Review D, single-column (~3.4 in wide).
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- rcParams
mpl.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 10,
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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# colorblind-friendly (Wong palette)
C_HYP    = "#0072B2"   # blue   -- mass shell
C_CONE   = "#D55E00"   # vermillion -- light cone / null
C_POINT  = "#000000"   # marker P
C_BOOST  = "#009E73"   # boost arrow
C_FAINT  = "#9970AB"   # massless-limit faded point

# ---------------------------------------------------------------- data
xmin, xmax = -2.5, 2.5
ymin, ymax = 0.0, 3.0

# hyperbola upper branch
eta = np.linspace(-2.4, 2.4, 800)
px = np.sinh(eta)
Ey = np.cosh(eta)

# marked point P and boosted point
eta0 = 0.9
delta = 0.55
P  = (np.sinh(eta0),        np.cosh(eta0))
Pb = (np.sinh(eta0+delta),  np.cosh(eta0+delta))

# massless-limit faded point (far up, near asymptote)
eta_far = 1.60
Pf = (np.sinh(eta_far), np.cosh(eta_far))

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(figsize=(3.4, 3.05))
ax.set_aspect("equal")

# light cone / null asymptotes  E = +/- p
xline = np.array([xmin, xmax])
ax.plot(xline,  xline, color=C_CONE, lw=1.0, ls=(0, (5, 3)), zorder=2)
ax.plot(xline, -xline, color=C_CONE, lw=1.0, ls=(0, (5, 3)), zorder=2)
# shade the causal (timelike, E>|p|) interior very lightly
xf = np.linspace(xmin, xmax, 200)
ax.fill_between(xf, np.abs(xf), ymax, color=C_CONE, alpha=0.05, lw=0, zorder=0)

# mass-shell hyperbola (bold)
ax.plot(px, Ey, color=C_HYP, lw=2.0, zorder=4)

# vertex label
ax.plot(0, 1, marker="o", ms=3.2, color=C_HYP, zorder=5)

# ---- boost: small arrow along the curve from P toward Pb
eta_arc = np.linspace(eta0, eta0 + delta, 40)
ax.annotate("", xy=(np.sinh(eta0+delta), np.cosh(eta0+delta)),
            xytext=(np.sinh(eta0+delta-0.14), np.cosh(eta0+delta-0.14)),
            arrowprops=dict(arrowstyle="-|>", color=C_BOOST, lw=1.4,
                            mutation_scale=11), zorder=6)
ax.plot(np.sinh(eta_arc), np.cosh(eta_arc), color=C_BOOST, lw=2.2,
        alpha=0.9, solid_capstyle="round", zorder=5)

# ---- markers
ax.plot(*P, marker="o", ms=5.0, color=C_POINT, zorder=7)
ax.plot(*Pb, marker="o", ms=4.0, mfc="white", mec=C_BOOST,
        mew=1.2, zorder=7)
# faded massless-limit point
ax.plot(*Pf, marker="o", ms=4.2, color=C_FAINT, alpha=0.55, zorder=6)

# ---------------------------------------------------------------- annotations
# mass shell
ax.annotate(r"$\Phi_1^{2}-\Phi_2^{2}=1$", xy=(-1.55, np.cosh(np.arcsinh(-1.55))),
            xytext=(-2.35, 2.55), color=C_HYP, fontsize=8.8,
            ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color=C_HYP, lw=0.7,
                            connectionstyle="arc3,rad=0.2",
                            shrinkA=1, shrinkB=2))
ax.text(-2.35, 2.30, "mass shell", color=C_HYP, fontsize=7.8,
        ha="left", va="center")

# point P
ax.annotate(r"$P=(\sinh\eta_0,\ \cosh\eta_0)$",
            xy=P, xytext=(P[0] + 0.15, P[1] - 0.42),
            fontsize=7.8, color=C_POINT, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=C_POINT, lw=0.6,
                            shrinkA=0, shrinkB=3))
ax.text(P[0] + 0.15, P[1] - 0.66, r"$\eta_0=0.9$", fontsize=7.4,
        color=C_POINT, ha="left", va="center")

# boost label
ax.annotate("boost", xy=(np.sinh(eta0+0.30), np.cosh(eta0+0.30)),
            xytext=(np.sinh(eta0+0.30) + 0.30, np.cosh(eta0+0.30) - 0.05),
            fontsize=7.8, color=C_BOOST, ha="left", va="center",
            arrowprops=dict(arrowstyle="->", color=C_BOOST, lw=0.7,
                            shrinkA=1, shrinkB=2))
ax.text(np.sinh(eta0+0.30) + 0.30, np.cosh(eta0+0.30) - 0.27,
        r"$\eta_0\!\to\!\eta_0+\delta$", fontsize=7.0, color=C_BOOST,
        ha="left", va="center")

# null directions phi_+  and  phi_- (placed low on each asymptote, in the clear)
ax.text(0.42, 0.54, r"$\phi_{+}\!:\ E=+p$", color=C_CONE, fontsize=7.8,
        ha="left", va="bottom", rotation=45, rotation_mode="anchor")
ax.text(-0.70, 0.82, r"$\phi_{-}\!:\ E=-p$", color=C_CONE, fontsize=7.8,
        ha="right", va="bottom", rotation=-45, rotation_mode="anchor")

# massless limit
ax.annotate(r"massless limit ($\eta\to\infty$)",
            xy=Pf, xytext=(Pf[0] - 0.15, Pf[1] + 0.28),
            fontsize=7.4, color=C_FAINT, ha="right", va="center",
            arrowprops=dict(arrowstyle="->", color=C_FAINT, lw=0.7,
                            shrinkA=1, shrinkB=2))

# ---------------------------------------------------------------- axes cosmetics
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel(r"$\Phi_2=\sinh\eta=p/m$")
ax.set_ylabel(r"$\Phi_1=\cosh\eta=E/m$")

ax.set_xticks(np.arange(-2, 3, 1))
ax.set_yticks(np.arange(0, 4, 1))

for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout(pad=0.4)

# ---------------------------------------------------------------- save
out = ("C:/Users/artio/Documents/Scientific work/"
       "Физика/Статьи/2026.07/PRD/code/figures/fig_hyperbola")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=300)
plt.close(fig)

print("fig_hyperbola: P =", P)
print("fig_hyperbola: boosted point =", Pb)
print("fig_hyperbola: massless-limit point =", Pf)
print("fig_hyperbola: saved PDF + PNG")
