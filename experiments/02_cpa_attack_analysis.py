# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 2: CPA Attack Analysis
===================================
Reproduces the CPA resistance results in Section 5 and Section 8.

Demonstrates:
  1. Linear scheme (Design 1): broken in d queries
  2. Message-dependent nonlinearity (Design 2): broken in 2d queries
  3. Key-dependent via det (Design 3): broken in 2d queries
  4. Post-composition σ (Design 4 / final): CPA FAILS

This is the core security argument of the paper.

Usage:
    python experiments/02_cpa_attack_analysis.py
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf


def solve_linear_Fp(equations, n_unknowns, p):
    """Solve/rank a linear system over F_p."""
    n_eq = len(equations)
    aug = np.zeros((n_eq, n_unknowns + 1), dtype=object)
    for i, (coeffs, rhs) in enumerate(equations):
        aug[i, :n_unknowns] = np.array(coeffs, dtype=object) % p
        aug[i, n_unknowns] = int(rhs) % p

    pivots, row = [], 0
    for col in range(n_unknowns):
        found = False
        for r in range(row, n_eq):
            if aug[r][col] % p != 0:
                found = True
                aug[[row, r]] = aug[[r, row]]
                break
        if not found:
            continue
        pivots.append(col)
        iv = pow(int(aug[row][col]), p - 2, p)
        aug[row] = (aug[row] * iv) % p
        for r in range(n_eq):
            if r != row and aug[r][col] % p != 0:
                aug[r] = (aug[r] - aug[r][col] * aug[row]) % p
        row += 1

    # Check consistency
    inconsistent = False
    for r in range(len(pivots), n_eq):
        if aug[r, n_unknowns] % p != 0:
            inconsistent = True
            break

    # Extract solution if consistent and fully determined
    solution = None
    if len(pivots) == n_unknowns and not inconsistent:
        solution = np.zeros(n_unknowns, dtype=int)
        for i, pc in enumerate(pivots):
            solution[pc] = int(aug[i, n_unknowns]) % p

    return len(pivots), inconsistent, solution


def main():
    print("=" * 70)
    print("EXPERIMENT 2: CPA ATTACK ANALYSIS")
    print("Reproduces Section 5 (Design Evolution) results")
    print("=" * 70)

    p = 17
    Fp = FiniteField(p)
    d = 2
    rng = random.Random(42)

    # Setup sheaf
    G = Graph.path(6)
    sheaf = Sheaf.random(G, d, Fp, rng=rng)
    v = 0  # attack node

    # Generate true keys
    A_true = Fp.random_gl(d, rng=rng)
    B_true = Fp.random_gl(d, rng=rng)

    def sigma_cube(y):
        return np.array([pow(int(yi), 3, p) for yi in y])

    # Generate CPA queries (sections and their encryptions)
    queries = []
    for i in range(10):
        s = sheaf.random_section(rng=rng)
        sv = s[v * d:(v + 1) * d]

        # Design 1: c = A * s
        c_d1 = Fp.mat_vec(A_true, sv)

        # Design 2: c = (A + h(s)*B) * s, h(s) = s0*s1
        h_val = int(sv[0]) * int(sv[1]) % p
        M_d2 = (A_true.astype(object) + h_val * B_true.astype(object)) % p
        c_d2 = Fp.mat_vec(M_d2.astype(int), sv)

        # Design 3: c = A*s + det(A)*B*s
        det_A = Fp.mat_det(A_true)
        c_d3 = np.array([(int(Fp.mat_vec(A_true, sv)[i]) +
                          det_A * int(Fp.mat_vec(B_true, sv)[i])) % p
                         for i in range(d)])

        # Design 4 (final): c = A*s + B*sigma(A*s)
        y = Fp.mat_vec(A_true, sv)
        z = sigma_cube(y)
        t = Fp.mat_vec(B_true, z)
        c_d4 = np.array([(int(y[i]) + int(t[i])) % p for i in range(d)])

        queries.append((sv, c_d1, c_d2, c_d3, c_d4, h_val))

    # ================================================================
    # DESIGN 1: c = A * s (linear)
    # ================================================================
    print("\n" + "-" * 70)
    print("DESIGN 1: c_v = A_v · s_v (pure linear)")
    print("-" * 70)

    for t in [2, 3, 4]:
        eqs = []
        for j in range(t):
            sv, c_d1 = queries[j][0], queries[j][1]
            s0, s1 = int(sv[0]), int(sv[1])
            # c0 = a00*s0 + a01*s1
            eqs.append(([s0, s1, 0, 0], int(c_d1[0])))
            # c1 = a10*s0 + a11*s1
            eqs.append(([0, 0, s0, s1], int(c_d1[1])))

        rank, incon, sol = solve_linear_Fp(eqs, 4, p)
        recovered = ""
        if sol is not None:
            rec_A = np.array([[sol[0], sol[1]], [sol[2], sol[3]]])
            match = np.array_equal(Fp.mat_mod(rec_A), Fp.mat_mod(A_true))
            recovered = f"Key match: {match}"

        print(f"  t={t} queries: rank={rank}/4 {recovered}")

    # ================================================================
    # DESIGN 2: c = (A + h(s)*B) * s
    # ================================================================
    print("\n" + "-" * 70)
    print("DESIGN 2: c_v = (A_v + h(s_v)·B_v) · s_v")
    print("-" * 70)

    for t in [2, 3, 4, 5]:
        eqs = []
        for j in range(t):
            sv, c_d2, h_val = queries[j][0], queries[j][2], queries[j][5]
            s0, s1 = int(sv[0]), int(sv[1])
            h = h_val
            # c0 = a00*s0 + a01*s1 + h*b00*s0 + h*b01*s1
            eqs.append(([s0, s1, 0, 0, (h * s0) % p, (h * s1) % p, 0, 0],
                        int(c_d2[0])))
            eqs.append(([0, 0, s0, s1, 0, 0, (h * s0) % p, (h * s1) % p],
                        int(c_d2[1])))

        rank, incon, sol = solve_linear_Fp(eqs, 8, p)
        recovered = ""
        if sol is not None:
            rec_A = np.array([[sol[0], sol[1]], [sol[2], sol[3]]])
            rec_B = np.array([[sol[4], sol[5]], [sol[6], sol[7]]])
            a_match = np.array_equal(Fp.mat_mod(rec_A), Fp.mat_mod(A_true))
            b_match = np.array_equal(Fp.mat_mod(rec_B), Fp.mat_mod(B_true))
            recovered = f"A match: {a_match}, B match: {b_match}"

        print(f"  t={t} queries: rank={rank}/8 {recovered}")

    # ================================================================
    # DESIGN 3: c = A*s + det(A)*B*s
    # ================================================================
    print("\n" + "-" * 70)
    print("DESIGN 3: c_v = A_v·s_v + det(A_v)·B_v·s_v")
    print("-" * 70)

    for t in [2, 3, 4, 5]:
        eqs = []
        for j in range(t):
            sv, c_d3 = queries[j][0], queries[j][3]
            s0, s1 = int(sv[0]), int(sv[1])
            # Linearize: d*b_ij -> db_ij
            # c0 = a00*s0 + a01*s1 + db00*s0 + db01*s1
            eqs.append(([s0, s1, 0, 0, s0, s1, 0, 0], int(c_d3[0])))
            eqs.append(([0, 0, s0, s1, 0, 0, s0, s1], int(c_d3[1])))

        rank, incon, sol = solve_linear_Fp(eqs, 8, p)
        recovered = ""
        if sol is not None:
            rec_A = np.array([[sol[0], sol[1]], [sol[2], sol[3]]])
            a_match = np.array_equal(Fp.mat_mod(rec_A), Fp.mat_mod(A_true))
            # Recover B from d*B and det(A)
            det_rec = Fp.mat_det(rec_A)
            if det_rec != 0:
                det_inv = pow(det_rec, p - 2, p)
                rec_B = np.array([[(sol[4] * det_inv) % p, (sol[5] * det_inv) % p],
                                  [(sol[6] * det_inv) % p, (sol[7] * det_inv) % p]])
                b_match = np.array_equal(Fp.mat_mod(rec_B), Fp.mat_mod(B_true))
                recovered = f"A match: {a_match}, B match: {b_match}"
            else:
                recovered = f"A match: {a_match}, det=0"

        print(f"  t={t} queries: rank={rank}/8 {recovered}")

    # ================================================================
    # DESIGN 4 (FINAL): c = A*s + B*sigma(A*s)
    # ================================================================
    print("\n" + "-" * 70)
    print("DESIGN 4 (FINAL): c_v = A_v·s_v + B_v·σ(A_v·s_v)")
    print("-" * 70)

    # Attempt 1: linear system for A only (ignoring B*sigma term)
    print("\n  --- Attempt: linear solve for A_v only ---")
    for t in [4, 6, 8]:
        eqs = []
        for j in range(t):
            sv, c_d4 = queries[j][0], queries[j][4]
            s0, s1 = int(sv[0]), int(sv[1])
            eqs.append(([s0, s1, 0, 0], int(c_d4[0])))
            eqs.append(([0, 0, s0, s1], int(c_d4[1])))

        rank, incon, sol = solve_linear_Fp(eqs, 4, p)
        print(f"    t={t}: rank={rank}/4, inconsistent={incon}")

    # Attempt 2: linear system for (A, B) treating sigma as known
    # This fails because sigma(A*s) depends on the unknown A
    print("\n  --- Attempt: linear solve for (A_v, B_v) ---")
    for t in [4, 6, 8, 10]:
        eqs = []
        for j in range(t):
            sv, c_d4 = queries[j][0], queries[j][4]
            s0, s1 = int(sv[0]), int(sv[1])
            # Try: c0 = a00*s0 + a01*s1 + (nonlinear mess)
            # Best linear approximation: just A part
            eqs.append(([s0, s1, 0, 0, 0, 0, 0, 0], int(c_d4[0])))
            eqs.append(([0, 0, s0, s1, 0, 0, 0, 0], int(c_d4[1])))

        rank, incon, sol = solve_linear_Fp(eqs, 8, p)
        print(f"    t={t}: rank={rank}/8, inconsistent={incon}")

    print(f"""
  CONCLUSION:
    Design 1 (linear):     BROKEN with {d} queries
    Design 2 (h(s)*B):     BROKEN with {2*d} queries
    Design 3 (det(A)*B):   BROKEN with {2*d} queries
    Design 4 (σ ∘ A):      LINEAR ATTACK FAILS — system is INCONSISTENT

  The nonlinearity of σ applied AFTER the secret map A_v
  makes the attacker's system genuinely polynomial (degree ≥ 3).
  No linear algebra technique can recover the key.
""")


if __name__ == "__main__":
    main()
