"""
Finite field arithmetic over F_p.

Provides modular arithmetic, matrix operations, and random
sampling over prime fields. All operations are exact (no floats).
"""

import numpy as np
from typing import Optional


class FiniteField:
    """Arithmetic over the prime field F_p."""

    def __init__(self, p: int):
        if p < 2:
            raise ValueError(f"Prime must be >= 2, got {p}")
        self.p = p

    def mod(self, x: int) -> int:
        return int(x) % self.p

    def add(self, a: int, b: int) -> int:
        return (int(a) + int(b)) % self.p

    def sub(self, a: int, b: int) -> int:
        return (int(a) - int(b)) % self.p

    def mul(self, a: int, b: int) -> int:
        return (int(a) * int(b)) % self.p

    def inv(self, x: int) -> int:
        """Multiplicative inverse via Fermat's little theorem."""
        x = int(x) % self.p
        if x == 0:
            raise ZeroDivisionError("Cannot invert 0 in F_p")
        return pow(x, self.p - 2, self.p)

    def inv_or_zero(self, x: int) -> int:
        """Multiplicative inverse, mapping 0 → 0."""
        x = int(x) % self.p
        return pow(x, self.p - 2, self.p) if x != 0 else 0

    def pow(self, base: int, exp: int) -> int:
        return pow(int(base), int(exp), self.p)

    # --- Matrix operations (object dtype for exact arithmetic) ---

    def mat_mod(self, M: np.ndarray) -> np.ndarray:
        return np.vectorize(lambda x: int(x) % self.p)(M.astype(object))

    def mat_mul(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return self.mat_mod(A.astype(object) @ B.astype(object))

    def mat_vec(self, A: np.ndarray, v: np.ndarray) -> np.ndarray:
        return self.mat_mod(
            A.astype(object) @ v.astype(object).reshape(-1, 1)
        ).flatten()

    def mat_inv(self, M: np.ndarray) -> np.ndarray:
        """Invert a d×d matrix over F_p via Gauss-Jordan."""
        d = M.shape[0]
        aug = np.zeros((d, 2 * d), dtype=object)
        for i in range(d):
            for j in range(d):
                aug[i, j] = int(M[i, j]) % self.p
            aug[i, d + i] = 1

        for col in range(d):
            pivot = -1
            for r in range(col, d):
                if aug[r, col] % self.p != 0:
                    pivot = r
                    break
            if pivot == -1:
                raise ValueError("Singular matrix — not invertible over F_p")
            aug[[col, pivot]] = aug[[pivot, col]]
            iv = pow(int(aug[col, col]), self.p - 2, self.p)
            aug[col] = (aug[col] * iv) % self.p
            for r in range(d):
                if r != col and aug[r, col] % self.p != 0:
                    aug[r] = (aug[r] - aug[r, col] * aug[col]) % self.p

        return aug[:, d:].astype(int) % self.p

    def mat_det(self, M: np.ndarray) -> int:
        """Determinant over F_p via row reduction."""
        d = M.shape[0]
        A = np.array(M, dtype=object) % self.p
        sign = 1
        for col in range(d):
            pivot = -1
            for r in range(col, d):
                if A[r, col] % self.p != 0:
                    pivot = r
                    break
            if pivot == -1:
                return 0
            if pivot != col:
                A[[col, pivot]] = A[[pivot, col]]
                sign *= -1
            iv = pow(int(A[col, col]), self.p - 2, self.p)
            for r in range(col + 1, d):
                if A[r, col] % self.p != 0:
                    factor = (A[r, col] * iv) % self.p
                    A[r] = (A[r] - factor * A[col]) % self.p
        result = sign
        for i in range(d):
            result = (result * int(A[i, i])) % self.p
        return int(result) % self.p

    def random_gl(self, d: int, rng=None) -> np.ndarray:
        """Sample a uniformly random invertible d×d matrix over F_p."""
        import random as _random
        _rng = rng or _random
        while True:
            M = np.array(
                [[_rng.randint(0, self.p - 1) for _ in range(d)] for _ in range(d)],
                dtype=object,
            )
            try:
                self.mat_inv(M)
                return M.astype(int) % self.p
            except ValueError:
                continue

    def random_vec(self, d: int, rng=None) -> np.ndarray:
        import random as _random
        _rng = rng or _random
        return np.array([_rng.randint(0, self.p - 1) for _ in range(d)])

    def kernel(self, M: np.ndarray) -> list[np.ndarray]:
        """Compute kernel of M over F_p via row reduction."""
        rows, cols = M.shape
        A = np.array(M, dtype=object) % self.p
        pivot_cols, row = [], 0
        for col in range(cols):
            found = False
            for r in range(row, rows):
                if A[r][col] % self.p != 0:
                    found = True
                    A[[row, r]] = A[[r, row]]
                    break
            if not found:
                continue
            pivot_cols.append(col)
            iv = pow(int(A[row][col]), self.p - 2, self.p)
            A[row] = (A[row] * iv) % self.p
            for r in range(rows):
                if r != row and A[r][col] % self.p != 0:
                    A[r] = (A[r] - A[r][col] * A[row]) % self.p
            row += 1
        free_cols = [c for c in range(cols) if c not in pivot_cols]
        vecs = []
        for fc in free_cols:
            vec = np.zeros(cols, dtype=object)
            vec[fc] = 1
            for i, pc in enumerate(pivot_cols):
                vec[pc] = (-A[i][fc]) % self.p
            vecs.append((vec % self.p).astype(int))
        return vecs

    def cube_invertible(self) -> bool:
        """Check if x^3 is a bijection on F_p (requires gcd(3, p-1) = 1)."""
        from math import gcd
        return gcd(3, self.p - 1) == 1

    def __repr__(self):
        return f"F_{self.p}"
