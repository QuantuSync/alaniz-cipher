# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 4: Sigma Selection Analysis
========================================
Reproduces Table 6 in the paper.

Evaluates five candidate σ functions on:
  1. Bijectivity over F_p^d
  2. Algebraic degree
  3. Differential uniformity
  4. Round-trip correctness

Usage:
    python experiments/04_sigma_selection.py
"""

import sys
import os
import random
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


P = 17


def sigma_inv(x):
    return np.array([pow(int(xi), P - 2, P) if int(xi) % P != 0 else 0
                     for xi in x])


def sigma_cube(x):
    return np.array([pow(int(xi), 3, P) for xi in x])


def sigma_chebyshev(x):
    """T_3(x) = 4x^3 - 3x"""
    return np.array([(4 * int(xi) ** 3 - 3 * int(xi)) % P for xi in x])


def sigma_mixed_quad(x):
    """(x0^2 + x1, x0*x1 + x1^2)"""
    return np.array([(int(x[0]) ** 2 + int(x[1])) % P,
                     (int(x[0]) * int(x[1]) + int(x[1]) ** 2) % P])


def sigma_lai_massey(x):
    """Lai-Massey with f(t) = t^3"""
    t = (int(x[0]) - int(x[1])) % P
    ft = pow(t, 3, P)
    return np.array([(int(x[0]) + ft) % P, (int(x[1]) + ft) % P])


def main():
    print("=" * 70)
    print("EXPERIMENT 4: SIGMA SELECTION ANALYSIS")
    print("Reproduces Table 6 in the paper")
    print("=" * 70)

    sigmas = {
        "inverse (x^{p-2})": sigma_inv,
        "cube (x^3)": sigma_cube,
        "Chebyshev T_3": sigma_chebyshev,
        "Mixed quadratic": sigma_mixed_quad,
        "Lai-Massey": sigma_lai_massey,
    }

    d = 2

    for name, sigma in sigmas.items():
        print(f"\n--- σ = {name} ---")

        # 1. Bijectivity
        images = set()
        for x0 in range(P):
            for x1 in range(P):
                y = sigma(np.array([x0, x1]))
                images.add((int(y[0]), int(y[1])))
        bijective = len(images) == P * P
        print(f"  Bijective: {bijective} ({len(images)}/{P * P} images)")

        # 2. Differential uniformity
        max_diff = 0
        for dx0 in range(P):
            for dx1 in range(P):
                if dx0 == 0 and dx1 == 0:
                    continue
                output_diffs = defaultdict(int)
                for x0 in range(P):
                    for x1 in range(P):
                        x = np.array([x0, x1])
                        xp = np.array([(x0 + dx0) % P, (x1 + dx1) % P])
                        dy = tuple((sigma(xp) - sigma(x)) % P)
                        output_diffs[dy] += 1
                local_max = max(output_diffs.values())
                max_diff = max(max_diff, local_max)
        print(f"  Differential uniformity: {max_diff} (lower = better)")

        # 3. Round-trip test with protocol
        if bijective or name == "Mixed quadratic":
            # Test with protocol using brute-force sigma
            Fp = FiniteField(P)
            G = Graph.path(6)
            rng = random.Random(42)
            sheaf = Sheaf.random(G, d, Fp, rng=rng)

            A_key = {v: Fp.random_gl(d, rng=rng) for v in G.nodes}
            B_key = {v: Fp.random_gl(d, rng=rng) for v in G.nodes}

            correct = 0
            n_tests = 11

            for _ in range(n_tests):
                s = sheaf.random_section(rng=rng)
                # Encrypt manually with this sigma
                c = np.zeros(sheaf.C0_dim, dtype=int)
                for v in G.nodes:
                    sv = sheaf.get_node_value(s, v)
                    y = Fp.mat_vec(A_key[v], sv)
                    z = sigma(y)
                    t = Fp.mat_vec(B_key[v], z)
                    cv = np.array([(int(y[i]) + int(t[i])) % P
                                   for i in range(d)])
                    sheaf.set_node_value(c, v, cv)

                # Decrypt via brute-force at root + propagation
                R_to = sheaf.tree_propagation_maps(root=0)
                root_sols = []
                c_root = sheaf.get_node_value(c, 0)
                for y0 in range(P):
                    for y1 in range(P):
                        y = np.array([y0, y1])
                        z = sigma(y)
                        t = Fp.mat_vec(B_key[0], z)
                        res = np.array([(int(y[i]) + int(t[i])) % P
                                        for i in range(d)])
                        if np.array_equal(res, c_root):
                            root_sols.append(y)

                A0_inv = Fp.mat_inv(A_key[0])
                found = False
                for y_r in root_sols:
                    x_r = Fp.mat_vec(A0_inv, y_r)
                    cand = np.zeros(sheaf.C0_dim, dtype=int)
                    for v in G.nodes:
                        sheaf.set_node_value(cand, v,
                                             Fp.mat_vec(R_to[v], x_r))
                    # Verify
                    ok = True
                    for v in G.nodes:
                        sv = sheaf.get_node_value(cand, v)
                        y = Fp.mat_vec(A_key[v], sv)
                        z = sigma(y)
                        t = Fp.mat_vec(B_key[v], z)
                        cv_check = np.array([(int(y[i]) + int(t[i])) % P
                                             for i in range(d)])
                        if not np.array_equal(cv_check,
                                              sheaf.get_node_value(c, v)):
                            ok = False
                            break
                    if ok and np.array_equal(cand, s):
                        found = True
                        break
                if found:
                    correct += 1

            print(f"  Round-trip: {correct}/{n_tests}")
        else:
            print(f"  Round-trip: SKIPPED (not bijective)")

    print(f"""
{'=' * 70}
SUMMARY
{'=' * 70}
  RECOMMENDED:
    σ_inv (x^{{p-2}}): Maximum algebraic degree, proven in AES,
                       bijective for all primes.
    σ_cube (x^3):     Clean degree-3, field-independent bijectivity
                       (when gcd(3,p-1)=1), used in MiMC/Poseidon.

  ELIMINATED:
    Chebyshev T_3:    NOT bijective (121/{P*P} images).
    Mixed quadratic:  NOT bijective (189/{P*P} images), degree only 2.
    Lai-Massey:       Maximal diff. uniformity ({P*P}), near-linear.
""")


if __name__ == "__main__":
    main()
