#!/usr/bin/env python3
"""verify_sec3_closedforms.py

Symbolic verification of the closed-form primitives of Sec. III
(Proposition "Closed forms"): the antiderivatives of the three functionals
tau1, tau2, x_approx are checked by differentiation (must return the integrand)
and their a -> infinity limits are checked against the massless-limit theorems.

Tested with Python 3.13 and sympy 1.14.
"""

import sys
import sympy as sp

eta, a = sp.symbols('eta a', positive=True)

_failures = []


def check(name, expr):
    residual = sp.simplify(expr)
    ok = (residual == 0)
    print("[{}] {}".format("PASS" if ok else "FAIL", name))
    if not ok:
        print("     residual = {}".format(residual))
        _failures.append(name)


# integrands
tanh = sp.tanh(eta)
f = a * tanh / sp.sqrt(1 + a**2 * tanh**2)
integrand_tau1 = sp.sqrt(1 - f**2)          # = cosh/sqrt(cosh^2 + a^2 sinh^2)
integrand_tau2 = f
integrand_x = f * sp.cosh(eta)

s1 = sp.sqrt(1 + a**2)

# closed forms (eta > 0 branch)
tau1_cf = sp.asinh(s1 * sp.sinh(eta)) / s1
tau2_cf = (a / s1) * (sp.acosh(s1 / a * sp.cosh(eta)) - sp.acosh(s1 / a))
x_cf = a / (1 + a**2) * (sp.sqrt(sp.cosh(eta)**2 + a**2 * sp.sinh(eta)**2) - 1)

# --- differentiation checks --------------------------------------------------
# For eta>0 both the derivative of each (increasing) closed form and the
# corresponding integrand are strictly positive, so equality is proved by
# F'^2 - g^2 == 0 (symbolic) together with the shared positive sign. sympy's
# simplify cannot collapse the raw nested-sqrt difference, but it collapses the
# squared difference. We additionally confirm on a high-precision numeric grid.

def check_deriv(name, F, g):
    d = sp.diff(F, eta)
    # square-and-sign proof (valid since both sides > 0 for eta>0, a>0)
    check(name + " [F'^2 = g^2, both > 0]", sp.simplify(d**2 - g**2))
    # numeric grid confirmation to 30 digits
    import mpmath as mp
    mp.mp.dps = 40
    dl = sp.lambdify((eta, a), d, 'mpmath')
    gl = sp.lambdify((eta, a), g, 'mpmath')
    worst = mp.mpf(0)
    for ev in [0.2, 0.5, 1.0, 2.0, 4.0]:
        for av in [0.5, 1.0, 5.0, 10.0, 50.0]:
            worst = max(worst, abs(dl(ev, av) - gl(ev, av)))
    ok = worst < mp.mpf(10) ** (-25)
    print("[{}] {} [numeric grid, max|F'-g| = {:.2e}]".format(
        "PASS" if ok else "FAIL", name, float(worst)))
    if not ok:
        _failures.append(name + " numeric")

check_deriv("d/deta tau1_closed = sqrt(1 - f^2)", tau1_cf, integrand_tau1)
check_deriv("d/deta tau2_closed = f", tau2_cf, integrand_tau2)
check_deriv("d/deta x_approx_closed = f cosh eta", x_cf, integrand_x)

# --- boundary values at eta = 0 (each functional vanishes) -------------------
check("tau1_closed(0) = 0", tau1_cf.subs(eta, 0))
check("tau2_closed(0) = 0", tau2_cf.subs(eta, 0))
check("x_approx_closed(0) = 0", x_cf.subs(eta, 0))

# --- a -> infinity limits (massless / null limit), eta > 0 fixed -------------
# Substitute a fixed positive numeric eta to evaluate the limits cleanly.
eta0 = sp.Rational(1)  # eta = 1 > 0
check("lim_{a->oo} tau1_closed = 0 (proper time collapses)",
      sp.limit(tau1_cf.subs(eta, eta0), a, sp.oo))
check("lim_{a->oo} tau2_closed = eta (affine -> |eta|)",
      sp.limit(tau2_cf.subs(eta, eta0), a, sp.oo) - eta0)
check("lim_{a->oo} x_approx_closed = sinh(eta)  (x -> t, null x=t)",
      sp.limit(x_cf.subs(eta, eta0), a, sp.oo) - sp.sinh(eta0))

# --- ceiling velocity and effective-mass ceiling -----------------------------
check("v_inf(a) = a/sqrt(1+a^2)", sp.limit(f, eta, sp.oo) - a / s1)
check("1 - v_inf^2 = 1/(1+a^2)", 1 - (a / s1)**2 - 1 / (1 + a**2))

# --- summary -----------------------------------------------------------------
print()
if _failures:
    print("{} check(s) FAILED".format(len(_failures)))
    sys.exit(1)
print("All closed-form checks passed.")
sys.exit(0)
