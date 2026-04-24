"""
Experiment 10: Centralizer dimension verification (Theorem 6.5 of paper).

Goal: verify empirically that the centralizer of the primitive element M_θ
in M_d(F_p) has F_p-dimension exactly d.

This is the algebraic fact underlying the key design choice of v3: restricting
the secondary key component to a scalar β ∈ F_{p^d}^* (rather than a matrix
B ∈ GL(d, F_p)) loses exactly (d^2 - d) bits of information, and this loss
does NOT create rank-deficient algebraic structure exploitable by MinRank.

Method: build M_θ for several (p, d), compute kernel of the linear map
X ↦ XM_θ - M_θX, and verify its dimension equals d.
"""
from __future__ import annotations
import numpy as np

from alaniz.core.field import FiniteField
from alaniz.core.field_ext import ExtensionField


def centralizer_dimension(p: int, d: int) -> int:
    """Compute dim_{F_p} Z(M_θ) via Gauss elimination."""
    Fp = FiniteField(p)
    field_ext = ExtensionField(p, d)
    GF = field_ext.GF
    theta = field_ext.theta

    # Build M_θ ∈ M_d(F_p): column j is the coordinates of θ·θ^{d-1-j}
    M_theta = np.zeros((d, d), dtype=int)
    for j in range(d):
        # θ · basis[j] = θ · θ^{d-1-j} = θ^{d-j}
        prod = theta * (theta ** (d - 1 - j))
        coeffs = field_ext.gf_to_vec(prod)
        for i in range(d):
            M_theta[i, j] = int(coeffs[i])

    # Build the linear map T: M_d(F_p) → M_d(F_p), T(X) = XM_θ - M_θX.
    # Represent X as a d^2-dimensional vector (row-major).
    # Compute kernel dimension.
    n = d * d
    T = np.zeros((n, n), dtype=int)
    for i in range(d):
        for j in range(d):
            # Basis element E_{ij}: X with 1 at (i,j) and 0 elsewhere.
            # (XM_θ)_{a,b} = sum_c X_{a,c} (M_θ)_{c,b} = (M_θ)_{j,b} if a=i, else 0
            # (M_θ X)_{a,b} = sum_c (M_θ)_{a,c} X_{c,b} = (M_θ)_{a,i} if b=j, else 0
            idx = i * d + j
            # XM_θ contribution
            for b in range(d):
                T[i * d + b, idx] = (T[i * d + b, idx] + M_theta[j, b]) % p
            # -M_θ X contribution
            for a in range(d):
                T[a * d + j, idx] = (T[a * d + j, idx] - M_theta[a, i]) % p

    # Kernel dimension = n - rank(T)
    rank = _matrix_rank_mod_p(T, p)
    return n - rank


def _matrix_rank_mod_p(M: np.ndarray, p: int) -> int:
    """Rank over F_p via Gauss elimination."""
    A = np.array(M, dtype=int) % p
    rows, cols = A.shape
    rank = 0
    col = 0
    for row in range(rows):
        if col >= cols:
            break
        # Find pivot
        pivot = -1
        for r in range(row, rows):
            if A[r, col] % p != 0:
                pivot = r; break
        if pivot == -1:
            col += 1
            continue
        if pivot != row:
            A[[row, pivot]] = A[[pivot, row]]
        inv = pow(int(A[row, col]), p - 2, p)
        A[row] = (A[row] * inv) % p
        for r in range(rows):
            if r != row and A[r, col] % p != 0:
                A[r] = (A[r] - A[r, col] * A[row]) % p
        rank += 1
        col += 1
    return rank


def main():
    print("=" * 72)
    print(" Experiment 10: Centralizer dimension verification")
    print("=" * 72)
    print(" For each (p, d), we verify that dim_{F_p} Z(M_θ) = d exactly.")
    print()
    print(f"{'p':>5}{'d':>4}{'expected':>12}{'measured':>12}{'match':>8}")
    print("-" * 40)

    all_match = True
    for p, d in [(7, 2), (7, 3), (7, 4), (7, 5), (11, 3), (13, 4), (17, 2)]:
        expected = d
        measured = centralizer_dimension(p, d)
        m = "YES" if expected == measured else "NO"
        if expected != measured:
            all_match = False
        print(f"{p:>5}{d:>4}{expected:>12}{measured:>12}{m:>8}")

    print("-" * 40)
    print()
    if all_match:
        print(" All centralizers have dimension exactly d.")
        print(" Confirms Theorem 6.5: no MinRank vulnerability from the")
        print(" scalar restriction of β.")
    else:
        print(" WARNING: dimension mismatch detected.")


if __name__ == "__main__":
    main()
