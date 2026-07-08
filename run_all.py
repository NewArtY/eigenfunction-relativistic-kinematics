# -*- coding: utf-8 -*-
"""run_all.py -- Reproduce every result and figure of the paper.

Driver that runs, in order, all symbolic/numeric verification scripts, the
Table I recomputation, and both figure scripts.  Each step is executed as a
subprocess; PASS/FAIL is reported per step and a final summary is printed.
The driver exits with a nonzero status if any step fails, so it doubles as a
continuous-integration entry point.

Usage:
    python run_all.py

Tested with Python 3.13, sympy 1.14, numpy 2.4, matplotlib 3.10, mpmath 1.3.
"""
import os
import sys
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, path relative to HERE)
STEPS = [
    ("verify_sec2 (40 checks)", os.path.join("verification", "verify_sec2.py")),
    ("verify_sec3 (14 checks)", os.path.join("verification", "verify_sec3.py")),
    ("verify_sec3_closedforms", os.path.join("verification",
                                             "verify_sec3_closedforms.py")),
    ("compute_table1 (Table I)", "compute_table1.py"),
    ("fig_veta.pdf", os.path.join("figures", "fig_veta.py")),
    ("fig_proptime.pdf", os.path.join("figures", "fig_proptime.py")),
    ("fig_finterp.pdf", os.path.join("figures", "fig_finterp.py")),
    ("fig_hyperbola.pdf", os.path.join("figures", "fig_hyperbola.py")),
]

# figures that must exist afterwards
EXPECTED_ARTIFACTS = [
    os.path.join("figures", "table1.csv"),
    os.path.join("figures", "fig_veta.pdf"),
    os.path.join("figures", "fig_veta.png"),
    os.path.join("figures", "fig_proptime.pdf"),
    os.path.join("figures", "fig_proptime.png"),
    os.path.join("figures", "fig_finterp.pdf"),
    os.path.join("figures", "fig_finterp.png"),
    os.path.join("figures", "fig_hyperbola.pdf"),
    os.path.join("figures", "fig_hyperbola.png"),
]


def run_step(label, relpath):
    path = os.path.join(HERE, relpath)
    print("=" * 72)
    print(">>> {}".format(label))
    print("    $ python {}".format(relpath.replace(os.sep, "/")))
    print("-" * 72)
    t0 = time.time()
    proc = subprocess.run([sys.executable, path], cwd=HERE)
    dt = time.time() - t0
    ok = (proc.returncode == 0)
    print("-" * 72)
    print("<<< {} : {} (rc={}, {:.1f}s)".format(
        label, "PASS" if ok else "FAIL", proc.returncode, dt))
    print()
    return ok


def main():
    print("#" * 72)
    print("# run_all.py  --  reproducing all verification and figure outputs")
    print("# python: {}".format(sys.version.split()[0]))
    print("#" * 72)
    print()

    results = []
    for label, relpath in STEPS:
        results.append((label, run_step(label, relpath)))

    # artifact existence check
    print("=" * 72)
    print("ARTIFACT CHECK")
    print("-" * 72)
    artifacts_ok = True
    for rel in EXPECTED_ARTIFACTS:
        p = os.path.join(HERE, rel)
        exists = os.path.isfile(p)
        size = os.path.getsize(p) if exists else 0
        artifacts_ok = artifacts_ok and exists
        print("  [{}] {} ({} bytes)".format(
            "OK" if exists else "MISSING", rel.replace(os.sep, "/"), size))
    print()

    # summary
    print("=" * 72)
    print("SUMMARY")
    print("-" * 72)
    n_fail = 0
    for label, ok in results:
        print("  [{}] {}".format("PASS" if ok else "FAIL", label))
        if not ok:
            n_fail += 1
    print("  [{}] all expected artifacts present".format(
        "PASS" if artifacts_ok else "FAIL"))
    print("=" * 72)

    total_fail = n_fail + (0 if artifacts_ok else 1)
    if total_fail:
        print("RESULT: {} step(s) FAILED.".format(total_fail))
        return 1
    print("RESULT: all {} steps passed; all artifacts produced.".format(
        len(results)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
