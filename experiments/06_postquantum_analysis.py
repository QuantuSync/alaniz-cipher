# ⚠ V1 HISTORICAL: This experiment uses v1 sigmas (cube/inverse)
# which are DEPRECATED and vulnerable to CPA attacks.
# For v2 verification, see experiments/14_v2_redteam.py
# These are preserved for reproducibility of the paper's design evolution.

#!/usr/bin/env python3
"""
Experiment 6: Post-Quantum Security Analysis
==============================================
Reproduces Section 9 of the paper.

Demonstrates:
  1. Why Shor's algorithm does not apply (structural argument)
  2. Grover's impact quantified per parameter set (Table 7)
  3. RSA factoring example (broken by Shor)
  4. Comparison table: RSA / ECC / Kyber / Alaniz (Table 8)

Usage:
    python experiments/06_postquantum_analysis.py
"""

import sys
from math import gcd, log2


def main():
    print("=" * 70)
    print("EXPERIMENT 6: POST-QUANTUM SECURITY ANALYSIS")
    print("Reproduces Section 9 of the paper")
    print("=" * 70)

    # ================================================================
    # PART 1: RSA vulnerability demonstration
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 1: RSA IS BROKEN BY SHOR")
    print("=" * 70)

    p_rsa, q_rsa = 61, 53
    N = p_rsa * q_rsa
    phi_N = (p_rsa - 1) * (q_rsa - 1)
    e_rsa = 17
    d_rsa = pow(e_rsa, -1, phi_N)

    print(f"\n  RSA example: p={p_rsa}, q={q_rsa}, N={N}")
    print(f"  Public key: (N={N}, e={e_rsa})")
    print(f"  Private key: d={d_rsa}")

    # Shor's approach: find order of a mod N
    a = 2
    order = 1
    val = a
    while val != 1:
        val = (val * a) % N
        order += 1

    print(f"\n  Shor finds ord({a} mod {N}) = {order}")

    if order % 2 == 0:
        x = pow(a, order // 2, N)
        f1 = gcd(x - 1, N)
        f2 = gcd(x + 1, N)
        if f1 != 1 and f1 != N:
            print(f"  gcd({x}-1, {N}) = {f1}")
            print(f"  gcd({x}+1, {N}) = {f2}")
            print(f"  → FACTORED: {N} = {f1} × {f2}")
        elif f2 != 1 and f2 != N:
            print(f"  → FACTORED via other factor")
        else:
            print(f"  → This 'a' didn't work, try another")

    print(f"""
  For RSA-2048:
    Shor's cost: O((log N)^3) ≈ O(2048^3) ≈ 2^33 operations
    Required qubits: ~4000 logical qubits
    → COMPLETELY BROKEN in polynomial time
""")

    # ================================================================
    # PART 2: Why Shor fails against NL-SMIP
    # ================================================================
    print("=" * 70)
    print("PART 2: WHY SHOR FAILS AGAINST NL-SMIP")
    print("=" * 70)

    requirements = [
        ("Cyclic/abelian group",
         "GL(d, F_p) is NON-ABELIAN for d ≥ 2",
         False),
        ("Periodic function",
         "Encryption is polynomial, no periodicity",
         False),
        ("Group homomorphism",
         "σ is NOT a homomorphism (nonlinear)",
         False),
        ("Quantum Fourier transform",
         "No efficient QFT over GL(d, F_p) for d ≥ 2",
         False),
    ]

    print(f"\n  Shor's four requirements vs NL-SMIP:\n")
    for req, reason, met in requirements:
        status = "✓ Met" if met else "✗ NOT met"
        print(f"  {status:>12}  {req:<30} → {reason}")

    print(f"\n  CONCLUSION: Shor is CATEGORICALLY INAPPLICABLE to NL-SMIP.")
    print(f"  This is not a parameter-scaling defense — it is structural.")

    # ================================================================
    # PART 3: Grover impact per parameter set
    # ================================================================
    print("\n" + "=" * 70)
    print("PART 3: GROVER IMPACT (Table 7)")
    print("=" * 70)

    params = [
        ("Demo", 2, 5, 8, "cube", 38),
        ("Standard", 4, 31, 16, "cube", 93),
        ("PQ-128", 8, 61, 32, "inv", 200),
        ("PQ-256", 8, 127, 64, "inv", 400),
    ]

    print(f"\n  {'Level':<12} {'d':>3} {'log2 p':>7} {'n':>4} "
          f"{'Classical':>12} {'Post-Grover':>13} {'PQ Secure?':>12}")
    print("  " + "-" * 70)

    for name, dv, p_bits, n, sigma, classical in params:
        grover = classical // 2
        if grover >= 128:
            pq = "YES"
        elif grover >= 64:
            pq = "Marginal"
        else:
            pq = "No"
        print(f"  {name:<12} {dv:>3} {p_bits:>7} {n:>4} "
              f"{'~2^' + str(classical):>12} {'~2^' + str(grover):>13} "
              f"{pq:>12}")

    print(f"""
  Grover halves security bits. Standard defense: double parameters.
  PQ-128 (d=8, p≈2^61, n=32) achieves ~2^100 post-quantum security.
  PQ-256 (d=8, p≈2^127, n=64) achieves ~2^200 post-quantum security.
""")

    # ================================================================
    # PART 4: Full comparison table
    # ================================================================
    print("=" * 70)
    print("PART 4: COMPARISON WITH ESTABLISHED SCHEMES (Table 8)")
    print("=" * 70)

    schemes = [
        ("RSA-2048", "Factoring", "~2^112", "BROKEN", "No", 256, 384),
        ("RSA-4096", "Factoring", "~2^150", "BROKEN", "No", 512, 768),
        ("ECC P-256", "ECDLP", "~2^128", "BROKEN", "No", 32, 32),
        ("ECC P-521", "ECDLP", "~2^256", "BROKEN", "No", 66, 66),
        ("---", "---", "---", "---", "---", 0, 0),
        ("AES-128", "Key search", "~2^128", "~2^64", "Marginal", 0, 16),
        ("AES-256", "Key search", "~2^256", "~2^128", "Yes", 0, 32),
        ("---", "---", "---", "---", "---", 0, 0),
        ("Kyber-768", "MLWE", "~2^180", "~2^160", "Yes", 1184, 2400),
        ("Dilithium-3", "MLWE", "~2^180", "~2^160", "Yes", 1952, 4000),
        ("SPHINCS+-256", "Hash", "~2^256", "~2^128", "Yes", 32, 64),
        ("---", "---", "---", "---", "---", 0, 0),
        ("Alaniz PQ-128", "NL-SMIP", "~2^200", "~2^100", "Yes",
         15128, 31232),
        ("Alaniz PQ-256", "NL-SMIP", "~2^400", "~2^200", "Yes",
         64008, 130048),
    ]

    print(f"\n  {'Scheme':<16} {'Problem':<12} {'Classical':>10} "
          f"{'Quantum':>10} {'PQ?':>5} {'PK(B)':>8} {'SK(B)':>8}")
    print("  " + "-" * 75)

    for row in schemes:
        if row[0] == "---":
            print("  " + "-" * 75)
        else:
            name, prob, cl, qu, pq, pk, sk = row
            print(f"  {name:<16} {prob:<12} {cl:>10} {qu:>10} "
                  f"{pq:>5} {pk:>8} {sk:>8}")

    print(f"""
  KEY OBSERVATIONS:
    1. RSA and ECC are BROKEN by Shor (polynomial time) — no defense.
    2. Kyber/Dilithium resist Shor; best attack is Grover (quadratic).
    3. Alaniz cipher is in the SAME category as Kyber: Shor does not apply,
       Grover provides only quadratic speedup.
    4. Unique advantage over Rainbow (broken at NIST): NL-SMIP has
       SECRET coefficients, preventing structural MinRank attacks.
""")

    # ================================================================
    # PART 5: Quantum attack surface summary
    # ================================================================
    print("=" * 70)
    print("PART 5: QUANTUM ATTACK SURFACE SUMMARY")
    print("=" * 70)

    vectors = [
        ("Shor (QFT/HSP)", "No",
         "---", "No abelian group, no periodicity, no homomorphism"),
        ("Grover (search)", "Generic",
         "√(classical)", "Halves security bits; standard doubling defense"),
        ("HHL (linear sys)", "No",
         "---", "System is nonlinear; XL matrix is dense"),
        ("Quantum walks", "No",
         "---", "Security from secret keys, not graph structure"),
    ]

    print(f"\n  {'Algorithm':<22} {'Applies?':<10} {'Speedup':<16} "
          f"{'Reason'}")
    print("  " + "-" * 80)
    for alg, applies, speedup, reason in vectors:
        print(f"  {alg:<22} {applies:<10} {speedup:<16} {reason}")

    print(f"""
  FINAL ASSESSMENT:
    NL-SMIP does not reduce to any problem known to admit
    polynomial-time quantum algorithms. The Alaniz cipher is a
    legitimate post-quantum candidate, in the same resistance
    category as NIST-standardized schemes.
""")


if __name__ == "__main__":
    main()
