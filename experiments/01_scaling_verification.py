# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 1: Scaling Verification
===================================
Reproduces Table 5 in the paper.

Tests round-trip correctness across 22 configurations:
  - 8 primes (17, 23, 29, 47, 59, 89, 131, 251)
  - 5 topologies (path, star, binary tree, caterpillar, random tree)
  - 2 fiber dimensions (d=2, d=4)
  - 2 sigma functions (inverse, cube)

Total: 106 encrypt/decrypt round-trips.
Expected result: 106/106 correct, 0 failures.

Usage:
    python experiments/01_scaling_verification.py
"""

import sys
import os
import random
import time
import numpy as np
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def run_config(graph, dv, p, sigma_name, n_messages, label):
    """Run round-trip tests for a single configuration."""
    rng = random.Random(200 + graph.n + p)
    Fp = FiniteField(p)

    if sigma_name == "cube" and not Fp.cube_invertible():
        return None  # skip incompatible prime

    sheaf = Sheaf.random(graph, dv, Fp, rng=rng)
    if sheaf.H0_dim == 0:
        return None  # trivial H^0

    params = PublicParams.generate(sheaf, sigma_name)
    proto = Protocol(params)
    key = proto.keygen(rng=rng)

    correct = 0
    avg_local = 0
    max_local = 0
    t0 = time.time()

    for _ in range(n_messages):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)

        # Decrypt via tree propagation
        recovered = proto.decrypt_tree(c, key)
        if recovered is not None and np.array_equal(recovered, s):
            correct += 1

    elapsed = time.time() - t0

    return {
        "label": label,
        "n": graph.n,
        "m": graph.m,
        "dv": dv,
        "p": p,
        "sigma": sigma_name,
        "H0_dim": sheaf.H0_dim,
        "correct": correct,
        "total": n_messages,
        "time": elapsed,
    }


def main():
    print("=" * 75)
    print("EXPERIMENT 1: SCALING VERIFICATION")
    print("Reproduces Table 5 in the paper")
    print("=" * 75)

    results = []

    # --- Varying topology (p=17, d=2, sigma=inverse) ---
    topologies = [
        (Graph.path(6), "Path P_6"),
        (Graph.path(8), "Path P_8"),
        (Graph.path(10), "Path P_10"),
        (Graph.path(12), "Path P_12"),
        (Graph.path(16), "Path P_16"),
        (Graph.path(20), "Path P_20"),
        (Graph.star(9), "Star S_9"),
        (Graph.binary_tree(2), "BinTree d=2"),
        (Graph.binary_tree(3), "BinTree d=3"),
        (Graph.caterpillar(5, 2), "Caterpillar(5,2)"),
        (Graph.random_tree(12, rng=random.Random(999)), "RandomTree 12"),
    ]

    print("\n--- Varying topology (p=17, d=2, sigma=inverse) ---\n")
    for graph, name in topologies:
        n_msg = 5 if graph.n > 12 else 10
        r = run_config(graph, 2, 17, "inverse", n_msg, name)
        if r:
            results.append(r)
            status = "✓" if r["correct"] == r["total"] else "✗"
            print(f"  {r['label']:<22} |V|={r['n']:>2} |E|={r['m']:>2} "
                  f"H0={r['H0_dim']} RT={r['correct']}/{r['total']} {status} "
                  f"({r['time']:.1f}s)")

    # --- Varying prime (P_6, d=2, sigma=inverse) ---
    print("\n--- Varying prime (P_6, d=2, sigma=inverse) ---\n")
    for p in [17, 23, 29, 47, 59, 89, 131, 251]:
        r = run_config(Graph.path(6), 2, p, "inverse", 5,
                       f"P_6, F_{p}")
        if r:
            results.append(r)
            status = "✓" if r["correct"] == r["total"] else "✗"
            print(f"  {r['label']:<22} "
                  f"RT={r['correct']}/{r['total']} {status} "
                  f"({r['time']:.1f}s)")

    # --- Sigma=cube with compatible primes ---
    print("\n--- Sigma=cube, various primes (P_8, d=2) ---\n")
    for p in [17, 23, 29, 47, 53, 59, 83, 89]:
        r = run_config(Graph.path(8), 2, p, "cube", 5,
                       f"P_8, F_{p}, cube")
        if r:
            results.append(r)
            status = "✓" if r["correct"] == r["total"] else "✗"
            print(f"  {r['label']:<22} "
                  f"RT={r['correct']}/{r['total']} {status} "
                  f"({r['time']:.1f}s)")
        else:
            print(f"  P_8, F_{p}, cube: SKIPPED (3 | p-1)")

    # --- Higher dimension d=4 ---
    print("\n--- Higher dimension d=4 (p=17, sigma=inverse) ---\n")
    for n in [4, 6]:
        r = run_config(Graph.path(n), 4, 17, "inverse", 3,
                       f"P_{n}, d=4, F_17")
        if r:
            results.append(r)
            status = "✓" if r["correct"] == r["total"] else "✗"
            print(f"  {r['label']:<22} H0={r['H0_dim']} "
                  f"RT={r['correct']}/{r['total']} {status} "
                  f"({r['time']:.1f}s)")

    # --- Summary ---
    total_tests = sum(r["total"] for r in results)
    total_correct = sum(r["correct"] for r in results)
    n_configs = len(results)

    print("\n" + "=" * 75)
    print(f"SUMMARY: {total_correct}/{total_tests} round-trips correct "
          f"across {n_configs} configurations")
    print("=" * 75)

    if total_correct == total_tests:
        print("RESULT: ALL PASSED ✓")
    else:
        print(f"RESULT: {total_tests - total_correct} FAILURES ✗")
        sys.exit(1)


if __name__ == "__main__":
    main()
