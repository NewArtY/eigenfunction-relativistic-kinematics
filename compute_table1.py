# -*- coding: utf-8 -*-
"""
compute_table1.py -- Independent high-precision recomputation of Table I.

Definitions:
    f(eta,a)  = a*tanh(eta) / sqrt(1 + a^2 tanh^2 eta)
    x_exact   = cosh(eta) - 1
    x_approx  = int_0^eta f(eta',a) cosh(eta') deta'
    tau_exact = sinh(eta)
    tau2      = int_0^eta f(eta',a) deta'
    tau_ph    = eta
    delta_x   = |x_approx - x_exact| / |x_exact| * 100 %
    delta_tau = |tau2 - tau_exact| / |tau_exact| * 100 %

High-precision quadrature via mpmath (dps=40).
"""
import os
import csv
import mpmath as mp

mp.mp.dps = 40

def f_mp(e, a):
    t = mp.tanh(e)
    return a*t/mp.sqrt(1 + a**2*t**2)

def x_approx(eta, a):
    return mp.quad(lambda u: f_mp(u, a)*mp.cosh(u), [0, eta])

def tau2(eta, a):
    return mp.quad(lambda u: f_mp(u, a), [0, eta])

def tau1(eta, a):
    return mp.quad(lambda u: mp.sqrt(1 - f_mp(u, a)**2), [0, eta])

ETAS = [mp.mpf(s) for s in ['0.5', '1.0', '1.5', '2.0', '3.0', '4.0']]
AS = [mp.mpf(10), mp.mpf(5), mp.mpf(50)]

# Draft values (a=10) to check against
DRAFT = {
    '0.5': dict(x_approx=0.42887, dx=236.0, tau2=0.40880, dtau=21.6),
    '1.0': dict(x_approx=1.07454, dx=97.9,  tau2=0.90219, dtau=23.2),
    '1.5': dict(x_approx=2.02201, dx=49.5,  tau2=1.39868),
    '2.0': dict(x_approx=3.51121, dx=27.1,  tau2=1.89587),
    '3.0': dict(x_approx=9.86964, dx=8.8,   tau2=2.89075),
    '4.0': dict(x_approx=27.05565, dx=2.8,  tau2=3.88576),
}

def compute_rows(a):
    rows = []
    for eta in ETAS:
        xe = mp.cosh(eta) - 1
        xa = x_approx(eta, a)
        dx = abs(xa - xe)/abs(xe)*100
        te = mp.sinh(eta)
        t2 = tau2(eta, a)
        dtau = abs(t2 - te)/abs(te)*100
        rows.append(dict(eta=eta, x_exact=xe, x_approx=xa, dx=dx,
                         tau_exact=te, tau2=t2, dtau=dtau, tau_ph=eta))
    return rows

def print_table(a, rows):
    print("\n" + "="*100)
    print(f"TABLE  (a = {mp.nstr(a,4)})")
    print("="*100)
    hdr = f"{'eta':>5} | {'x_exact':>12} | {'x_approx':>12} | {'d_x %':>9} | {'tau_exact':>12} | {'tau2':>12} | {'d_tau %':>9} | {'tau_ph':>6}"
    print(hdr)
    print("-"*len(hdr))
    for r in rows:
        print(f"{mp.nstr(r['eta'],3):>5} | {mp.nstr(r['x_exact'],8):>12} | "
              f"{mp.nstr(r['x_approx'],8):>12} | {mp.nstr(r['dx'],5):>9} | "
              f"{mp.nstr(r['tau_exact'],8):>12} | {mp.nstr(r['tau2'],8):>12} | "
              f"{mp.nstr(r['dtau'],5):>9} | {mp.nstr(r['tau_ph'],3):>6}")

def compare_draft(rows):
    print("\n" + "="*100)
    print("COMPARISON vs AI DRAFT  (a = 10)   [CORRECT = agrees to ~4 sig figs]")
    print("="*100)
    print(f"{'eta':>5} | {'quantity':>9} | {'draft':>12} | {'computed':>14} | {'verdict':>8}")
    print("-"*66)
    def verdict(comp, draft, rtol=5e-4):
        if draft == 0:
            return "n/a"
        return "CORRECT" if abs(comp-mp.mpf(draft))/abs(mp.mpf(draft)) <= rtol else "WRONG"
    for r in rows:
        key = mp.nstr(r['eta'], 2)
        # normalize key like '0.5','1.0'
        keymap = {'0.5':'0.5','1.0':'1.0','1.5':'1.5','2.0':'2.0','3.0':'3.0','4.0':'4.0'}
        # find matching draft key by float
        dk = None
        for k in DRAFT:
            if abs(float(r['eta']) - float(k)) < 1e-9:
                dk = k; break
        d = DRAFT[dk]
        print(f"{dk:>5} | {'x_approx':>9} | {d['x_approx']:>12.5f} | {mp.nstr(r['x_approx'],9):>14} | {verdict(r['x_approx'], d['x_approx']):>8}")
        print(f"{'':>5} | {'d_x %':>9} | {d['dx']:>12.4g} | {mp.nstr(r['dx'],6):>14} | {verdict(r['dx'], d['dx'], rtol=5e-3):>8}")
        print(f"{'':>5} | {'tau2':>9} | {d['tau2']:>12.5f} | {mp.nstr(r['tau2'],9):>14} | {verdict(r['tau2'], d['tau2']):>8}")
        if 'dtau' in d:
            print(f"{'':>5} | {'d_tau %':>9} | {d['dtau']:>12.4g} | {mp.nstr(r['dtau'],6):>14} | {verdict(r['dtau'], d['dtau'], rtol=5e-3):>8}")
        print("-"*66)

def write_csv(all_rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['a', 'eta', 'x_exact', 'x_approx', 'delta_x_percent',
                    'tau_exact', 'tau2', 'delta_tau_percent', 'tau_ph'])
        for a, rows in all_rows:
            for r in rows:
                w.writerow([mp.nstr(a,4), mp.nstr(r['eta'],4),
                            mp.nstr(r['x_exact'],16), mp.nstr(r['x_approx'],16),
                            mp.nstr(r['dx'],10), mp.nstr(r['tau_exact'],16),
                            mp.nstr(r['tau2'],16), mp.nstr(r['dtau'],10),
                            mp.nstr(r['tau_ph'],16)])

if __name__ == '__main__':
    all_rows = []
    rows10 = None
    for a in AS:
        rows = compute_rows(a)
        all_rows.append((a, rows))
        print_table(a, rows)
        if a == 10:
            rows10 = rows
    compare_draft(rows10)

    csv_path = os.path.join(os.path.dirname(__file__), 'figures', 'table1.csv')
    write_csv(all_rows, csv_path)
    print(f"\nWrote CSV -> {csv_path}")

    # Paper-ready Table I (a=10)
    print("\n" + "="*100)
    print("PAPER-READY TABLE I  (a = 10)")
    print("="*100)
    print(f"{'eta':>5} | {'x_exact':>10} | {'x_approx':>10} | {'d_x (%)':>8} | {'tau_exact':>10} | {'tau2':>10} | {'d_tau (%)':>9} | {'tau_ph':>6}")
    print("-"*84)
    for r in rows10:
        print(f"{float(r['eta']):>5.1f} | {float(r['x_exact']):>10.4f} | "
              f"{float(r['x_approx']):>10.4f} | {float(r['dx']):>8.1f} | "
              f"{float(r['tau_exact']):>10.4f} | {float(r['tau2']):>10.4f} | "
              f"{float(r['dtau']):>9.1f} | {float(r['tau_ph']):>6.1f}")
