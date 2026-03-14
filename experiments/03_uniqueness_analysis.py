# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 3: Decryption Uniqueness Analysis
==============================================
Reproduces the uniqueness results in Section 7 of the paper.

Measures:
  1. Local preimage distribution of g_v(y) = y + B*sigma(y)
  2. Collision probability gamma (Lemma 7.2)
  3. Tree propagation: impostor elimination rate
  4. Brute-force uniqueness confirmation

Usage:
    python experiments/03_uniqueness_analysis.py
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


def sigma_inv(y, p):
    return np.array([pow(int(yi), p - 2, p) if int(yi) % p != 0 else 0
                     for yi in y])


def sigma_cube(y, p):
    return np.array([pow(int(yi), 3, p) for yi in y])


def main():
    print("=" * 70)
    print("EXPERIMENT 3: DECRYPTION UNIQUENESS ANALYSIS")
    print("Reproduces Section 7 results")
    print("=" * 70)

    # ================================================================
    # PART 1: Local preimage distribution
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 1: PREIMAGE DISTRIBUTION OF g_v(y) = y + B·σ(y)")
    print("=" * 70)

    for p in [17, 23, 29, 47]:
        Fp = FiniteField(p)
        d = 2
        n_samples = 200

        # Distribution over random (B, c) pairs
        sol_counts = defaultdict(int)
        for trial in range(n_samples):
            rng = random.Random(trial + p * 1000)
            B = Fp.random_gl(d, rng=rng)
            c = np.array([rng.randint(0, p - 1) for _ in range(d)])

            # Count solutions to y + B*sigma(y) = c
            count = 0
            for y0 in range(p):
                for y1 in range(p):
                    y = np.array([y0, y1])
                    z = sigma_inv(y, p)
                    t = Fp.mat_vec(B, z)
                    res = np.array([(int(y[i]) + int(t[i])) % p
                                    for i in range(d)])
                    if np.array_equal(res, c):
                        count += 1
            sol_counts[count] += 1

        expected = sum(k * v for k, v in sol_counts.items()) / n_samples

        print(f"\n  F_{p}, σ=inverse, d=2 ({n_samples} random pairs):")
        for k in sorted(sol_counts.keys()):
            pct = sol_counts[k] / n_samples * 100
            print(f"    {k} solutions: {sol_counts[k]:4d} ({pct:5.1f}%)")
        print(f"    E[preimages] = {expected:.2f}")

    # ================================================================
    # PART 2: Collision probability gamma
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 2: COLLISION PROBABILITY γ (Lemma 7.2)")
    print("=" * 70)

    for p in [17, 23, 29]:
        Fp = FiniteField(p)
        d = 2

        total_pairs = 0
        collisions = 0

        for trial in range(100):
            rng = random.Random(trial + p * 9000)
            B = Fp.random_gl(d, rng=rng)

            # Build image table
            image_of = {}
            for y0 in range(p):
                for y1 in range(p):
                    y = np.array([y0, y1])
                    z = sigma_inv(y, p)
                    t = Fp.mat_vec(B, z)
                    result = tuple((int(y[i]) + int(t[i])) % p
                                   for i in range(d))
                    image_of[(y0, y1)] = result

            # Sample pairs and check collisions
            for _ in range(50):
                y1 = np.array([rng.randint(0, p - 1), rng.randint(0, p - 1)])
                y2 = np.array([rng.randint(0, p - 1), rng.randint(0, p - 1)])
                if np.array_equal(y1, y2):
                    continue
                total_pairs += 1
                if image_of[tuple(y1)] == image_of[tuple(y2)]:
                    collisions += 1

        gamma_emp = collisions / total_pairs if total_pairs > 0 else 0
        gamma_theory = 1 / p ** d

        print(f"\n  F_{p}, d={d}:")
        print(f"    Empirical γ  = {collisions}/{total_pairs} = {gamma_emp:.6f}")
        print(f"    Theoretical  = 1/p^d = {gamma_theory:.6f}")
        print(f"    Bezout bound = (3^d - 1)/(p^d - 1) = "
              f"{(3**d - 1) / (p**d - 1):.6f}")

    # ================================================================
    # PART 3: Tree propagation — impostor elimination
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 3: TREE PROPAGATION — IMPOSTOR ELIMINATION")
    print("=" * 70)

    for p in [17, 23, 29, 47]:
        Fp = FiniteField(p)
        d = 2
        n_nodes = 8
        G = Graph.path(n_nodes)
        rng = random.Random(p * 100)

        sheaf = Sheaf.random(G, d, Fp, rng=rng)
        R_to = sheaf.tree_propagation_maps(root=0)

        A_key = {v: Fp.random_gl(d, rng=rng) for v in G.nodes}
        B_key = {v: Fp.random_gl(d, rng=rng) for v in G.nodes}

        def encrypt_node(sv, Av, Bv):
            y = Fp.mat_vec(Av, sv)
            z = sigma_inv(y, p)
            t = Fp.mat_vec(Bv, z)
            return np.array([(int(y[i]) + int(t[i])) % p for i in range(d)])

        n_tests = 20
        all_unique = True
        impostor_fail_nodes = []

        for _ in range(n_tests):
            s_root = np.array([rng.randint(0, p - 1) for _ in range(d)])
            s_true = {v: Fp.mat_vec(R_to[v], s_root) for v in G.nodes}
            c = {v: encrypt_node(s_true[v], A_key[v], B_key[v])
                 for v in G.nodes}

            # Solve at root
            root_sols = []
            for y0 in range(p):
                for y1 in range(p):
                    y = np.array([y0, y1])
                    z = sigma_inv(y, p)
                    t = Fp.mat_vec(B_key[0], z)
                    res = np.array([(int(y[i]) + int(t[i])) % p
                                    for i in range(d)])
                    if np.array_equal(res, c[0]):
                        root_sols.append(y)

            A0_inv = Fp.mat_inv(A_key[0])
            n_pass = 0

            for y_r in root_sols:
                x_r = Fp.mat_vec(A0_inv, y_r)
                is_true = np.array_equal(x_r, s_true[0])

                # Propagate and verify
                passed_all = True
                fail_node = None
                for v in G.nodes:
                    sv = Fp.mat_vec(R_to[v], x_r)
                    cv = encrypt_node(sv, A_key[v], B_key[v])
                    if not np.array_equal(cv, c[v]):
                        passed_all = False
                        fail_node = v
                        break

                if passed_all:
                    n_pass += 1
                elif not is_true:
                    impostor_fail_nodes.append(fail_node)

            if n_pass != 1:
                all_unique = False

        print(f"\n  F_{p}, P_8 ({n_tests} tests):")
        print(f"    Root solutions (avg): ~{len(root_sols)}")
        print(f"    All unique: {all_unique}")
        if impostor_fail_nodes:
            from collections import Counter
            dist = Counter(impostor_fail_nodes)
            print(f"    Impostor first-failure node: {dict(dist)}")
        else:
            print(f"    No impostors survived (all failed immediately)")

    # ================================================================
    # PART 4: Brute-force uniqueness confirmation
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 4: BRUTE-FORCE UNIQUENESS CONFIRMATION")
    print("=" * 70)

    Fp = FiniteField(17)
    G = Graph.path(6)
    rng = random.Random(12345)
    sheaf = Sheaf.random(G, 2, Fp, rng=rng)
    params = PublicParams.generate(sheaf, "inverse")
    proto = Protocol(params)
    key = proto.keygen(rng=rng)

    n_tests = 10
    all_unique = True

    for t in range(n_tests):
        s = sheaf.random_section(rng=rng)
        c = proto.encrypt(s, key)
        valid = proto.decrypt_bruteforce(c, key)

        unique = len(valid) == 1 and np.array_equal(valid[0], s)
        if not unique:
            all_unique = False
        print(f"  Test {t}: brute-force solutions = {len(valid)}, "
              f"unique & correct = {unique}")

    print(f"\n  Result: {'ALL UNIQUE ✓' if all_unique else 'AMBIGUITY DETECTED ✗'}")


if __name__ == "__main__":
    main()
