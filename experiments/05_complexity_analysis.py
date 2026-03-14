# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 5: Algebraic Attack Complexity
==========================================
Reproduces Table 4 (XL complexity) in the paper.

Analyzes:
  1. Attacker's polynomial system structure (degree, variables)
  2. XL attack: minimum degree for success, matrix size, cost
  3. Bezout bounds per configuration
  4. z-substitution trade-off for sigma=inverse

Usage:
    python experiments/05_complexity_analysis.py
"""

import sys
from math import comb, log2


def main():
    print("=" * 70)
    print("EXPERIMENT 5: ALGEBRAIC ATTACK COMPLEXITY")
    print("Reproduces Table 4 in the paper")
    print("=" * 70)

    # ================================================================
    # XL ANALYSIS
    # ================================================================
    print("\n" + "=" * 70)
    print("XL ATTACK ANALYSIS PER NODE")
    print("=" * 70)

    configs = [
        # (label, dv, sigma, eq_degree, n_vars_fn, n_eq_fn)
        ("σ=cube, d=2",
         2, "cube", 4,
         lambda d, t: 2 * d ** 2,  # 8 vars
         lambda d, t: d * t),  # 2t equations

        ("σ=cube, d=4",
         4, "cube", 4,
         lambda d, t: 2 * d ** 2,  # 32 vars
         lambda d, t: d * t),

        ("σ=inv(z-sub), d=2",
         2, "inv", 2,
         lambda d, t: 2 * d ** 2 + 2 * d * t,  # 8 + 2*2*t
         lambda d, t: d * t + d * t),  # enc + z constraints

        ("σ=inv(z-sub), d=4",
         4, "inv", 2,
         lambda d, t: 2 * d ** 2 + 2 * d * t,  # 32 + 2*4*t
         lambda d, t: d * t + d * t),

        ("σ=inv(z-sub), d=8",
         8, "inv", 2,
         lambda d, t: 2 * d ** 2 + 2 * d * t,
         lambda d, t: d * t + d * t),
    ]

    for label, dv, sigma, eq_deg, n_vars_fn, n_eq_fn in configs:
        print(f"\n--- {label} ---")

        # Find minimum t where system is determined
        t_min = 1
        while True:
            nv = n_vars_fn(dv, t_min)
            ne = n_eq_fn(dv, t_min)
            if ne >= nv - dv:  # approximately determined
                break
            t_min += 1
            if t_min > 50:
                break

        nv = n_vars_fn(dv, t_min)
        ne = n_eq_fn(dv, t_min)

        print(f"  Min queries t={t_min}: {nv} variables, "
              f"{ne} equations of degree {eq_deg}")

        # Bezout bound
        bezout = eq_deg ** nv
        bezout_log2 = nv * log2(eq_deg)
        print(f"  Bezout bound: {eq_deg}^{nv} ≈ 2^{bezout_log2:.0f}")

        # XL degree sweep
        print(f"  {'D':>4} {'Monomials':>14} {'XL equations':>16} "
              f"{'Sufficient?':>12} {'Cost (log2)':>12}")

        xl_succeeded = False
        for D in range(eq_deg, 20):
            monomials = comb(nv + D, D)
            if D >= eq_deg:
                xl_eq = ne * comb(nv + D - eq_deg, D - eq_deg)
            else:
                xl_eq = 0

            sufficient = xl_eq >= monomials
            mat_log2 = log2(monomials) if monomials > 1 else 0
            cost_log2 = 2.37 * mat_log2

            marker = " <<<" if sufficient and not xl_succeeded else ""
            print(f"  {D:>4} {monomials:>14} {xl_eq:>16} "
                  f"{'YES' if sufficient else 'no':>12} "
                  f"{'~2^' + str(int(cost_log2)):>12}{marker}")

            if sufficient and not xl_succeeded:
                xl_succeeded = True
                print(f"\n  → XL succeeds at degree {D}")
                print(f"  → Matrix dimension: {monomials}")
                print(f"  → Linear algebra cost: ~2^{int(cost_log2)}")
                break

        if not xl_succeeded:
            print(f"\n  → XL does NOT converge up to degree 19")
            print(f"  → Attack cost: > 2^{int(2.37 * log2(comb(nv + 19, 19)))}")

    # ================================================================
    # PARAMETER RECOMMENDATIONS
    # ================================================================
    print("\n\n" + "=" * 70)
    print("RECOMMENDED PARAMETERS (Table 3 in paper)")
    print("=" * 70)

    params = [
        ("Demo", 2, 5, 8, "cube", 38, 19),
        ("Standard", 4, 31, 16, "cube", 93, 47),
        ("PQ-128", 8, 61, 32, "inv", 200, 100),
        ("PQ-256", 8, 127, 64, "inv", 400, 200),
    ]

    print(f"\n{'Level':<12} {'d':>3} {'log2(p)':>8} {'n':>4} {'σ':>6} "
          f"{'Classical':>12} {'Quantum':>10} {'SK (bytes)':>12}")
    print("-" * 75)

    for name, dv, p_bits, n, sigma, cl, qu in params:
        sk_bytes = 2 * n * dv ** 2 * p_bits // 8
        print(f"{name:<12} {dv:>3} {p_bits:>8} {n:>4} {sigma:>6} "
              f"{'~2^' + str(cl):>12} {'~2^' + str(qu):>10} "
              f"{sk_bytes:>12,}")

    # ================================================================
    # z-SUBSTITUTION TRADE-OFF
    # ================================================================
    print(f"\n\n" + "=" * 70)
    print("z-SUBSTITUTION TRADE-OFF (σ_inv vs σ_cube)")
    print("=" * 70)

    for dv in [4, 8]:
        t_min_cube = 2 * dv ** 2 // dv
        t_min_inv = t_min_cube

        nv_cube = 2 * dv ** 2
        nv_inv = 2 * dv ** 2 + 2 * dv * t_min_inv

        bezout_cube = 4 ** nv_cube
        bezout_inv = 2 ** nv_inv

        print(f"\n  d={dv}:")
        print(f"    σ_cube: {nv_cube} vars, degree 4, "
              f"Bezout 4^{nv_cube} = 2^{nv_cube * 2}")
        print(f"    σ_inv:  {nv_inv} vars, degree 2, "
              f"Bezout 2^{nv_inv}")
        print(f"    → σ_inv is {'HARDER' if nv_inv > nv_cube * 2 else 'comparable'} "
              f"despite lower degree")


if __name__ == "__main__":
    main()
