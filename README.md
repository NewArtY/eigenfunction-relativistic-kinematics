# A Unified Eigenfunction Approach to Massive and Massless Particle Kinematics in Special Relativity — verification & figure code

This repository contains the reproducible symbolic-verification and
figure-generation code accompanying the Physical Review D theory paper

> **A Unified Eigenfunction Approach to Massive and Massless Particle
> Kinematics in Special Relativity**
> N. S. Akintsov, A. P. Nevecheria, S. N. Andreev, and Qing-Hua Qin.

## Summary

The paper recasts relativistic single-particle kinematics as an
eigenfunction problem: the boost eigenfunctions `Phi_1(eta)=cosh eta`,
`Phi_2(eta)=sinh eta` (eigenvalue `+1` of `d^2/deta^2`) and the null light-cone
basis `phi_pm=e^{+/-eta}` generate the energy-momentum, the Einstein velocity
addition law, and the `SO(1,1)` boost group. A one-parameter velocity
functional

```
f(eta, a) = a*tanh(eta) / sqrt(1 + a^2 * tanh^2(eta)),   a > 0,
```

interpolates between the massive reference velocity `v = tanh(eta)` (as
`a -> 0` the leading slope is `a`) and the null plateau `|v| = 1` (as
`a -> infinity`, with ceiling velocity `v_inf(a) = a/sqrt(1+a^2)`). Its
proper-time-like and affine primitives `tau_1`, `tau_2` and the approximate
trajectory `x_approx` admit exact closed forms whose `a -> infinity` limits
reproduce the massless/null kinematics. This code verifies every one of those
identities symbolically and at high numerical precision, recomputes Table I,
and regenerates the figures.

The repository is archived on **Zenodo** and cited from the article under DOI
[10.5281/zenodo.21269570](https://doi.org/10.5281/zenodo.21269570), which is
also recorded in `CITATION.cff` and `.zenodo.json`.

## Repository layout

| Path | What it does |
|------|--------------|
| `verification/verify_sec2.py` | SymPy proof of the Sec. II / Appendix A identities: eigenvalue equation, Wronskian `= 1`, null basis `phi_pm`, invariant `I = 1`, `v = tanh eta`, energy-momentum `E^2-p^2=m^2` and `p_pm = m e^{+/-eta}`, Einstein velocity addition, and the `SO(1,1)` boost representation (det `= 1`, `B^T g B = g`, `exp(eta k)=B`, group law, null eigenvectors). **40 checks.** |
| `verification/verify_sec3.py` | Symbolic + `mpmath` (dps=50) verification of the Sec. III functional `f(eta,a)`: `f(0,a)=0`, oddness, `1-f^2 = 1/(1+a^2 tanh^2 eta)` (so `|f|<1`), massless limit `lim_{a->oo} f = 1`, small-`eta` expansions of `f`, `tau_2`, `x_approx`, and exact-motion identities. **14 checks.** |
| `verification/verify_sec3_closedforms.py` | Proves the closed-form primitives by differentiation (square-and-sign symbolic proof + 40-digit numeric grid), checks boundary values at `eta=0`, the `a -> infinity` massless limits of `tau_1`, `tau_2`, `x_approx`, and the ceiling velocity `v_inf(a)=a/sqrt(1+a^2)`. |
| `compute_table1.py` | Independent high-precision (`mpmath`, dps=40) recomputation of **Table I** for `a in {5,10,50}`; compares against the draft values and writes `figures/table1.csv`. |
| `figures/fig_proptime.py` | Figure: for `a=10`, `eta in [0,3]`, plots `t=sinh eta`, the closed-form `tau_1(eta,a)`, `tau_2(eta,a)`, and the photon line `tau_ph=eta`; numerically locates and annotates the real interior crossing. Writes `fig_proptime.pdf` / `.png`. |
| `figures/fig_finterp.py` | Figure: `f(eta,a)` for `a in {0.5,1,2,5,10,50}` over `eta in [-3,3]`, with the massive reference `v=tanh eta` and photon plateaus `v=+/-1`, showing the monotone approach to the null plateau. Writes `fig_finterp.pdf` / `.png`. |
| `figures/table1.csv` | Machine-readable Table I output. |
| `run_all.py` | Driver: runs every verification and figure script in order, checks that all artifacts were produced, prints a summary, and exits nonzero on any failure. |

## Environment

- Python **3.13**
- sympy **1.14**
- numpy **2.4** (≥ 2.0)
- matplotlib **3.10** (≥ 3.9)
- mpmath **1.3**

Install the dependencies with:

```
pip install -r requirements.txt
```

## Reproducing all results

From the repository root (`.../PRD/code`), run the single driver:

```
python run_all.py
```

This reproduces everything and prints a PASS/FAIL summary. To run the steps
individually:

```
python verification/verify_sec2.py
python verification/verify_sec3.py
python verification/verify_sec3_closedforms.py
python compute_table1.py
python figures/fig_proptime.py
python figures/fig_finterp.py
```

Each verification script prints `PASS`/`FAIL` per identity and exits `0` only
if all checks pass. The figure scripts write vector **PDF** (for the
manuscript) and **PNG** (300 dpi) into `figures/`.

## Verification

All verification scripts currently pass:

- `verify_sec2.py` — **40 checks**, all PASS, exit 0.
- `verify_sec3.py` — **14 checks**, all PASS, exit 0.
- `verify_sec3_closedforms.py` — all closed-form / limit checks PASS, exit 0.

## Paper equation / theorem → verifying script

| Paper item | Verified by |
|------------|-------------|
| Eigenvalue equation `d^2 Phi/deta^2 = Phi` for `Phi_1=cosh`, `Phi_2=sinh` | `verify_sec2.py` (a) |
| Wronskian `W[Phi_1,Phi_2] = 1` | `verify_sec2.py` (b) |
| Null light-cone basis `phi_pm = e^{+/-eta} = cosh +/- sinh`, `K phi_pm = +/- phi_pm` | `verify_sec2.py` (c) |
| Invariant `I = Phi_1^2 - Phi_2^2 = 1`, `v = Phi_2/Phi_1 = tanh eta` | `verify_sec2.py` (d) |
| Energy-momentum `E^2 - p^2 = m^2`, `p_pm = m e^{+/-eta} = E +/- p`, boost covariance, massless limit `p_- = m^2/P -> 0` | `verify_sec2.py` (d') |
| Einstein velocity addition `tanh(a+b) = (tanh a + tanh b)/(1+tanh a tanh b)` | `verify_sec2.py` (e) |
| `SO(1,1)` boost: `B(eta)(cosh,sinh)^T=(cosh(+),sinh(+))`, `det B=1`, `B^T g B=g`, `exp(eta k)=B`, group law `B(a)B(b)=B(a+b)`, null eigenvectors of `k` | `verify_sec2.py` (f) |
| Functional `f(eta,a)`: `f(0,a)=0`, oddness, `1-f^2 = 1/(1+a^2 tanh^2 eta)`, `|f|<1` | `verify_sec3.py` |
| Massless limit `lim_{a->oo} f = 1` | `verify_sec3.py`, `verify_sec3_closedforms.py` |
| Small-`eta` expansions of `f`, `tau_2`, `x_approx` (`~ a eta`, `~ (a/2)eta^2`) | `verify_sec3.py` |
| Closed forms `tau_1 = asinh(sqrt(1+a^2) sinh eta)/sqrt(1+a^2)`, `tau_2 = (a/sqrt(1+a^2))[acosh(sqrt(1+a^2)/a cosh eta) - acosh(sqrt(1+a^2)/a)]`, `x_approx` (verified by differentiation) | `verify_sec3_closedforms.py` |
| Boundary values `tau_1(0)=tau_2(0)=x_approx(0)=0` | `verify_sec3_closedforms.py` |
| `a -> infinity` limits: `tau_1 -> 0`, `tau_2 -> eta`, `x_approx -> sinh eta` (null `x = t`) | `verify_sec3_closedforms.py` |
| Ceiling velocity `v_inf(a) = a/sqrt(1+a^2)`, `1 - v_inf^2 = 1/(1+a^2)` | `verify_sec3_closedforms.py` |
| Table I (`x_exact`, `x_approx`, `delta_x`, `tau_exact`, `tau_2`, `delta_tau`, `tau_ph`) | `compute_table1.py` → `figures/table1.csv` |
| Proper-time / affine figure (`t`, `tau_1`, `tau_2`, `tau_ph`) | `figures/fig_proptime.py` |
| Velocity-interpolation figure (`f(eta,a)` → plateau) | `figures/fig_finterp.py` |

## Citation

See `CITATION.cff`. Please cite both the software (via its Zenodo DOI) and the
Physical Review D article.

## License

MIT — see `LICENSE`. Copyright (c) 2026 N. S. Akintsov, A. P. Nevecheria,
S. N. Andreev, and Qing-Hua Qin.
