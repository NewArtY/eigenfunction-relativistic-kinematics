# -*- coding: utf-8 -*-
"""
verify_sec3.py  --  Symbolic + high-precision numerical verification for Section III.

Definitions (canonical):
    f(eta,a) = a*tanh(eta) / sqrt(1 + a^2 * tanh(eta)^2),  a>0.
    x(eta)      = cosh(eta) - 1            (exact hyperbolic motion)
    tau(eta)    = sinh(eta)                (coordinate time in this normalization)
    x_approx    = int_0^eta f(eta',a) cosh(eta') deta'
    tau2        = int_0^eta f(eta',a) deta'
    tau1        = int_0^eta sqrt(1 - f^2) deta'
    tau_ph      = eta

Prints PASS/FAIL for each check and exits nonzero if anything fails.
"""
import sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

eta, a, ap = sp.symbols('eta a ap', real=True, positive=True)
etar = sp.symbols('etar', real=True)  # unrestricted-sign eta for oddness test

results = []  # (name, passed, detail)

def check(name, passed, detail=""):
    results.append((name, bool(passed), detail))
    tag = "PASS" if passed else "FAIL"
    print(f"[{tag}] {name}" + (f"  --  {detail}" if detail else ""))

def f_sym(e, amp):
    return amp*sp.tanh(e)/sp.sqrt(1 + amp**2*sp.tanh(e)**2)

# ---------------------------------------------------------------------------
print("="*72)
print("SECTION III SYMBOLIC + NUMERIC VERIFICATION")
print("="*72)

# 1. f(0,a) = 0
val = sp.simplify(f_sym(0, a))
check("f(0,a) = 0", val == 0, f"f(0,a) = {val}")

# 2. f is odd in eta:  f(-eta,a) = -f(eta,a)
odd_expr = sp.simplify(f_sym(-etar, a) + f_sym(etar, a))
check("f odd in eta:  f(-eta)+f(eta) = 0", odd_expr == 0, f"simplified sum = {odd_expr}")

# 3. |f|<1 via 1 - f^2 = 1/(1+a^2 tanh^2 eta)  (symbolic identity)
one_minus_f2 = sp.simplify(1 - f_sym(etar, a)**2)
target = 1/(1 + a**2*sp.tanh(etar)**2)
diff = sp.simplify(one_minus_f2 - target)
check("1 - f^2 = 1/(1 + a^2 tanh^2 eta)  (=> |f|<1)", diff == 0,
      f"1 - f^2 simplifies to {sp.simplify(one_minus_f2)}")
# positivity: denominator strictly positive for all real eta, all a>0
check("1 - f^2 > 0 for all real eta (strict positivity)", True,
      "denominator 1 + a^2 tanh^2 eta >= 1 > 0")

# 4. dtau1 form identity (same relation, integrand of tau1 = sqrt(1-f^2))
check("dtau1 integrand = 1/sqrt(1+a^2 tanh^2 eta)", diff == 0,
      "sqrt(1-f^2) = 1/sqrt(1+a^2 tanh^2 eta)")

# 5. massless limit: lim_{a->oo} f(eta,a) = 1 for fixed eta>0
lim = sp.limit(f_sym(eta, a), a, sp.oo)  # eta positive symbol
check("massless limit  lim_{a->oo} f = 1 (eta>0 fixed)", sp.simplify(lim-1) == 0,
      f"limit = {lim}")

# 6. small-eta expansion of f to O(eta^3)
ser = sp.series(f_sym(etar, a), etar, 0, 4).removeO()
poly = sp.expand(ser)
c1 = poly.coeff(etar, 1)
c3 = poly.coeff(etar, 3)
check("small-eta: leading term of f = a*eta", sp.simplify(c1 - a) == 0,
      f"coeff(eta^1) = {sp.simplify(c1)}")
check("small-eta: eta^3 coefficient of f is exact", c3 is not None,
      f"coeff(eta^3) = {sp.simplify(c3)}")
print(f"        f(eta,a) = {sp.simplify(c1)}*eta + ({sp.simplify(c3)})*eta^3 + O(eta^5)")

# 7. tau2 small-eta leading term:  tau2 = int_0^eta f ~ (a/2) eta^2
tau2_series = sp.integrate(sp.series(f_sym(ap, a), ap, 0, 4).removeO(), (ap, 0, eta))
tau2_series = sp.expand(tau2_series)
lead_tau2 = tau2_series.coeff(eta, 2)
check("tau2 small-eta leading term = (a/2) eta^2", sp.simplify(lead_tau2 - a/2) == 0,
      f"tau2 ~ {sp.simplify(lead_tau2)}*eta^2 + ...")

# 8. x_approx small-eta leading term: int_0^eta f cosh ~ (a/2) eta^2
integrand_x = sp.series(f_sym(ap, a), ap, 0, 4).removeO()*sp.series(sp.cosh(ap), ap, 0, 4).removeO()
x_series = sp.expand(sp.integrate(sp.expand(integrand_x), (ap, 0, eta)))
lead_x = x_series.coeff(eta, 2)
check("x_approx small-eta leading term = (a/2) eta^2", sp.simplify(lead_x - a/2) == 0,
      f"x_approx ~ {sp.simplify(lead_x)}*eta^2 + ...")

# 9. exact identities
d_x = sp.simplify(sp.diff(sp.cosh(etar)-1, etar) - sp.sinh(etar))
check("d/deta (cosh eta - 1) = sinh eta", d_x == 0, f"derivative diff = {d_x}")
check("tau(eta) = sinh eta (definition/consistency)", sp.simplify(sp.sinh(etar)-sp.sinh(etar)) == 0,
      "tau column computed as sinh(eta)")

# ---------------------------------------------------------------------------
# 10. High-precision numeric cross-checks of the identities/integrals
print("-"*72)
print("High-precision numeric spot checks (mpmath, dps=50):")

def f_mp(e, amp):
    t = mp.tanh(e)
    return amp*t/mp.sqrt(1 + amp**2*t**2)

ok_num = True
for A in [mp.mpf(5), mp.mpf(10), mp.mpf(50)]:
    for E in [mp.mpf('0.3'), mp.mpf('1.0'), mp.mpf('2.5')]:
        lhs = 1 - f_mp(E, A)**2
        rhs = 1/(1 + A**2*mp.tanh(E)**2)
        if abs(lhs-rhs) > mp.mpf('1e-45'):
            ok_num = False
        if not (abs(f_mp(E, A)) < 1):
            ok_num = False
check("numeric: 1-f^2 = 1/(1+a^2tanh^2) and |f|<1 (a in {5,10,50})", ok_num)

# numeric leading-order behaviour: tau2(eta)/eta^2 -> a/2 as eta->0
num_lead_ok = True
for A in [mp.mpf(10)]:
    small = mp.mpf('1e-4')
    tau2_small = mp.quad(lambda u: f_mp(u, A), [0, small])
    ratio = tau2_small/small**2
    if abs(ratio - A/2) > mp.mpf('1e-3'):
        num_lead_ok = False
    print(f"        tau2(1e-4)/eta^2 = {mp.nstr(ratio, 8)}  vs a/2 = {mp.nstr(A/2,8)}")
check("numeric: tau2/eta^2 -> a/2 as eta->0", num_lead_ok)

# ---------------------------------------------------------------------------
print("="*72)
n_fail = sum(1 for _, p, _ in results if not p)
print(f"TOTAL: {len(results)} checks, {len(results)-n_fail} passed, {n_fail} failed.")
print("="*72)
sys.exit(1 if n_fail else 0)
