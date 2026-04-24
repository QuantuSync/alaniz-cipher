"""
Extension field arithmetic over F_{p^d}.

Wraps the `galois` library to provide:
  - Field construction F_{p^d} = F_p[x]/(f(x)) for irreducible f
  - Primitive element θ and the ordered basis (θ^{d-1}, ..., θ^0)
  - Isomorphism ι: F_p^d ↔ F_{p^d} as F_p-vector spaces
  - Polynomial factorization via Cantor-Zassenhaus (for decryption)

Used by the v3 protocol, which performs the nonlinear operation π_e(x) = x^e
inside F_{p^d} rather than component-wise in F_p^d.
"""

from __future__ import annotations
import numpy as np
from math import gcd

try:
    import galois
except ImportError as e:
    raise ImportError(
        "The `galois` package is required for v3 (F_{p^d} arithmetic). "
        "Install with: pip install galois>=0.3"
    ) from e


class ExtensionField:
    """Arithmetic over F_{p^d} with the canonical F_p-basis."""

    def __init__(self, p: int, d: int):
        if d < 1:
            raise ValueError(f"Extension degree must be >= 1, got {d}")
        self.p = p
        self.d = d
        try:
            self.GF = galois.GF(p ** d)
        except LookupError:
            # Fall back to explicit irreducible polynomial
            irr = galois.irreducible_poly(p, d)
            self.GF = galois.GF(p ** d, irreducible_poly=irr)
        self.theta = self.GF.primitive_element

    # ---------------- Isomorphism F_p^d <-> F_{p^d} ----------------

    def vec_to_gf(self, v):
        """ι: F_p^d → F_{p^d} using basis (θ^{d-1}, ..., θ^0)."""
        g = self.GF(0)
        for i in range(self.d):
            g += self.GF(int(v[i]) % self.p) * (self.theta ** (self.d - 1 - i))
        return g

    def gf_to_vec(self, g) -> np.ndarray:
        """ι^{-1}: F_{p^d} → F_p^d."""
        coeffs = g.vector()
        arr = np.array([int(x) for x in coeffs], dtype=np.int64)
        if len(arr) < self.d:
            arr = np.concatenate([np.zeros(self.d - len(arr), dtype=np.int64), arr])
        return arr

    # ---------------- Exponent selection ----------------

    @staticmethod
    def find_secure_exponent(p: int, d: int, min_e: int = 17) -> int:
        """
        Find the least prime e >= min_e coprime with p^d - 1.

        For v3, e = 17 is the default; for (p, d) where gcd(17, p^d-1) != 1,
        the next prime (19, 23, ...) is used.
        """
        N = p ** d - 1
        candidate = min_e
        while candidate < min_e + 200:
            if _is_prime(candidate) and gcd(candidate, N) == 1:
                return candidate
            candidate += 1
        raise ValueError(
            f"No secure exponent found in [{min_e}, {min_e + 200}) "
            f"for p={p}, d={d}"
        )

    # ---------------- Polynomial factorization (for decryption) ----------------

    def find_roots(self, coeffs_high_to_low: list) -> list:
        """
        Find all roots of a polynomial in F_{p^d}[t] via Cantor-Zassenhaus.

        Args:
            coeffs_high_to_low: list of GF elements, coefficient of t^deg first.

        Returns:
            List of roots (GF elements) in F_{p^d}.
        """
        poly = galois.Poly(coeffs_high_to_low, field=self.GF)
        return list(poly.roots())

    # ---------------- Uniform sampling ----------------

    def random_nonzero(self, rng=None) -> "galois.FieldArray":
        """Uniform sample from F_{p^d}^* = F_{p^d} \\ {0}."""
        import random as _random
        r = rng or _random
        while True:
            val = r.randint(0, self.p ** self.d - 1)
            g = self.GF(val)
            if g != self.GF(0):
                return g

    def random_not_in(self, exclude: list, rng=None) -> "galois.FieldArray":
        """Uniform sample from F_{p^d} \\ exclude."""
        import random as _random
        r = rng or _random
        exclude_set = {int(x) for x in exclude}
        while True:
            val = r.randint(0, self.p ** self.d - 1)
            if val not in exclude_set:
                return self.GF(val)

    def __repr__(self):
        return f"F_{{{self.p}^{self.d}}}"


def _is_prime(n: int) -> bool:
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    i = 3
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True
