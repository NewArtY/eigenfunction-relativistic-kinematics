#!/usr/bin/env python3
"""verify_sec2.py

Symbolic verification of the identities in Sec. II / Appendix A of the paper
"A Unified Eigenfunction Approach to Massive and Massless Particle Kinematics
in Special Relativity."

Every identity is reduced to a normal form and compared to zero with
sympy.simplify.  The script prints PASS/FAIL for each check and exits with a
nonzero status if any check fails.

Tested with Python 3.13 and sympy 1.14.
"""

import sys
import sympy as sp

# --- symbols -----------------------------------------------------------------
eta, eta0, a, b, u, v, m, P = sp.symbols('eta eta0 a b u v m P', real=True)

_failures = []


def check(name, expr):
    """Assert `expr` is identically zero.

    Normalize hyperbolic/exponential forms to a common basis (rewrite to exp)
    before simplifying, so that identities mixing cosh/sinh with exp collapse
    reliably regardless of scalar prefactors.
    """
    residual = sp.simplify(sp.expand(sp.expand_trig(expr)).rewrite(sp.exp))
    if residual != 0:
        residual = sp.simplify(residual)
    ok = (residual == 0)
    print("[{}] {}".format("PASS" if ok else "FAIL", name))
    if not ok:
        print("     residual = {}".format(residual))
        _failures.append(name)


# --- definitions -------------------------------------------------------------
Phi1 = sp.cosh(eta)          # Phi_1(eta) = cosh eta
Phi2 = sp.sinh(eta)          # Phi_2(eta) = sinh eta


def X(f):                    # X = d^2/d eta^2
    return sp.diff(f, eta, 2)


def K(f):                    # K = d/d eta
    return sp.diff(f, eta, 1)


# --- (a) eigenvalue equation  X Phi = Phi  (lambda = 1) ----------------------
check("(a) X Phi1 = Phi1", X(Phi1) - Phi1)
check("(a) X Phi2 = Phi2", X(Phi2) - Phi2)

# --- (b) Wronskian W[Phi1, Phi2] = 1 -----------------------------------------
W = Phi1 * K(Phi2) - Phi2 * K(Phi1)
check("(b) Wronskian W[Phi1,Phi2] = 1", W - 1)

# --- (c) null / light-cone basis phi_pm = e^{+- eta} -------------------------
phi_plus = sp.exp(eta)
phi_minus = sp.exp(-eta)
check("(c) phi_+ = Phi1 + Phi2", phi_plus - (Phi1 + Phi2))
check("(c) phi_- = Phi1 - Phi2", phi_minus - (Phi1 - Phi2))
check("(c) K phi_+ = +phi_+", K(phi_plus) - phi_plus)
check("(c) K phi_- = -phi_-", K(phi_minus) + phi_minus)
check("(c) Phi1 = (phi_+ + phi_-)/2", Phi1 - (phi_plus + phi_minus) / 2)
check("(c) Phi2 = (phi_+ - phi_-)/2", Phi2 - (phi_plus - phi_minus) / 2)

# --- (d) invariant I = 1 and velocity v = tanh eta ---------------------------
I = Phi1**2 - Phi2**2
check("(d) I = Phi1^2 - Phi2^2 = 1", I - 1)
check("(d) v = Phi2/Phi1 = tanh eta", Phi2 / Phi1 - sp.tanh(eta))

# --- physical momentum: E = m cosh, p = m sinh, p_pm = m e^{+- eta} -----------
E = m * sp.cosh(eta)
p = m * sp.sinh(eta)
pplus = m * sp.exp(eta)
pminus = m * sp.exp(-eta)
check("(d') E^2 - p^2 = m^2", E**2 - p**2 - m**2)
check("(d') p_+ p_- = m^2", pplus * pminus - m**2)
check("(d') p_+ = m e^{eta} = E + p", pplus - (E + p))
check("(d') p_- = m e^{-eta} = E - p", pminus - (E - p))

# relative invariance (covariance) of p_pm under boost eta -> eta + eta0
check("(d') p_+ -> e^{eta0} p_+",
      m * sp.exp(eta + eta0) - sp.exp(eta0) * (m * sp.exp(eta)))
check("(d') p_- -> e^{-eta0} p_-",
      m * sp.exp(-(eta + eta0)) - sp.exp(-eta0) * (m * sp.exp(-eta)))
check("(d') product invariant p_+ p_- unchanged",
      (m * sp.exp(eta + eta0)) * (m * sp.exp(-(eta + eta0))) - m**2)

# massless limit: with p_+ = P held fixed, p_- = m^2/P -> 0 as m -> 0
check("(d') massless limit p_- = m^2/P -> 0", sp.limit(m**2 / P, m, 0))

# --- (e) Einstein velocity addition ------------------------------------------
check("(e) tanh(a+b) = (tanh a + tanh b)/(1 + tanh a tanh b)",
      sp.tanh(a + b) - (sp.tanh(a) + sp.tanh(b)) / (1 + sp.tanh(a) * sp.tanh(b)))
check("(e) tanh(atanh u + atanh v) = (u+v)/(1+uv)",
      sp.expand_trig(sp.tanh(sp.atanh(u) + sp.atanh(v))) - (u + v) / (1 + u * v))

# --- (f) boost additivity via the 2x2 vector representation ------------------
B = sp.Matrix([[sp.cosh(eta0), sp.sinh(eta0)],
               [sp.sinh(eta0), sp.cosh(eta0)]])
vec = sp.Matrix([sp.cosh(eta), sp.sinh(eta)])
out = B * vec
check("(f) B(eta0).(cosh,sinh)^T gives cosh(eta+eta0)",
      out[0] - sp.cosh(eta + eta0))
check("(f) B(eta0).(cosh,sinh)^T gives sinh(eta+eta0)",
      out[1] - sp.sinh(eta + eta0))

# SO(1,1): det B = 1 and B^T g B = g with g = diag(1,-1)
g = sp.diag(1, -1)
check("(f) det B(eta0) = 1", B.det() - 1)
gform = (B.T * g * B) - g
for i in range(2):
    for j in range(2):
        check("(f) B^T g B = g entry [{},{}]".format(i, j), gform[i, j])

# generator k = [[0,1],[1,0]],  exp(eta0 k) = B(eta0)
k = sp.Matrix([[0, 1], [1, 0]])
expm = (k * eta0).exp()
check("(f) exp(eta0 k)[0,0] = cosh eta0", expm[0, 0] - sp.cosh(eta0))
check("(f) exp(eta0 k)[0,1] = sinh eta0", expm[0, 1] - sp.sinh(eta0))
check("(f) exp(eta0 k)[1,0] = sinh eta0", expm[1, 0] - sp.sinh(eta0))
check("(f) exp(eta0 k)[1,1] = cosh eta0", expm[1, 1] - sp.cosh(eta0))

# group law B(a) B(b) = B(a+b)
Ba = sp.Matrix([[sp.cosh(a), sp.sinh(a)], [sp.sinh(a), sp.cosh(a)]])
Bb = sp.Matrix([[sp.cosh(b), sp.sinh(b)], [sp.sinh(b), sp.cosh(b)]])
Bab = sp.Matrix([[sp.cosh(a + b), sp.sinh(a + b)],
                 [sp.sinh(a + b), sp.cosh(a + b)]])
prod = Ba * Bb - Bab
for i in range(2):
    for j in range(2):
        check("(f) group law B(a)B(b)=B(a+b) entry [{},{}]".format(i, j),
              prod[i, j])

# null eigenvectors of the generator k: k e_pm = +- e_pm
e_plus = sp.Matrix([1, 1])
e_minus = sp.Matrix([1, -1])
check("(f) k e_+ = +e_+ entry 0", (k * e_plus - e_plus)[0])
check("(f) k e_+ = +e_+ entry 1", (k * e_plus - e_plus)[1])
check("(f) k e_- = -e_- entry 0", (k * e_minus + e_minus)[0])
check("(f) k e_- = -e_- entry 1", (k * e_minus + e_minus)[1])

# --- summary -----------------------------------------------------------------
print()
if _failures:
    print("{} check(s) FAILED:".format(len(_failures)))
    for name in _failures:
        print("  - {}".format(name))
    sys.exit(1)

print("All checks passed.")
sys.exit(0)
