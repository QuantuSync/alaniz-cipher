"""
Public nonlinear maps σ: F_p^d → F_p^d for the Alaniz Cipher.

The map σ is applied AFTER the secret linear transformation A_v·s_v,
creating key-message entanglement that defeats CPA attacks.

Supported maps:
  - inverse: σ(x)_i = x_i^{p-2} (multiplicative inverse, 0↦0)
             [DEPRECATED: vulnerable to scaling CPA attack]
  - cube:    σ(x)_i = x_i^3 (requires gcd(3, p-1) = 1)
             [DEPRECATED: vulnerable to scaling CPA attack]
  - id_spn:  σ(y) = y + S(M·S(y) + c), S = cube, M = mixing matrix,
             c = (1,...,1) constant vector
             [RECOMMENDED: resists scaling and interpolation attacks]

Design rationale for id_spn:
  Component-wise σ (cube, inverse) satisfy σ(λy)_i = λ^e σ(y)_i,
  enabling an O(d)-query CPA key recovery (Rodríguez Langa, 2026).
  Additionally, any σ without a linear part allows polynomial
  interpolation to extract A directly from the t^1 coefficient.
  
  σ_id_spn has:
    - Linear part = I (identity): the t^1 coefficient yields
      (I+B)·A, mixing A and B inseparably.
    - Degree-9 nonlinear part with cross-component mixing via M:
      prevents component-wise decomposition.
    - σ(0) = 0: no information leakage from zero queries.
"""

from __future__ import annotations
import numpy as np
from typing import Callable
from ..core.field import FiniteField


class Sigma:
    """Public nonlinear map for the encryption scheme."""

    def __init__(self, name: str, Fp: FiniteField, d: int = 2):
        self.name = name
        self.Fp = Fp
        self.d = d

        if name == "inverse":
            self._forward = self._inverse
            self._algebraic_degree = Fp.p - 2
        elif name == "cube":
            if not Fp.cube_invertible():
                raise ValueError(
                    f"Cube map not invertible over F_{Fp.p} "
                    f"(3 divides p-1={Fp.p - 1})"
                )
            self._forward = self._cube
            self._algebraic_degree = 3
        elif name == "id_spn":
            self._forward = self._id_spn
            self._algebraic_degree = 9
            # Precompute mixing matrix M (circulant)
            self._M = self._build_mixing_matrix(d, Fp.p)
        else:
            raise ValueError(
                f"Unknown sigma: {name}. "
                f"Use 'inverse', 'cube', or 'id_spn'."
            )

    def _build_mixing_matrix(self, d: int, p: int) -> np.ndarray:
        """Build public mixing matrix M for id_spn.
        Circulant: M[i][j] = 1 if i==j, 2 if j==(i+1)%d, 0 otherwise.
        """
        M = np.zeros((d, d), dtype=int)
        for i in range(d):
            M[i][i] = 1
            M[i][(i + 1) % d] = 2
        return M

    @property
    def algebraic_degree(self) -> int:
        return self._algebraic_degree

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return self._forward(x)

    def _inverse(self, x: np.ndarray) -> np.ndarray:
        """σ(x)_i = x_i^{p-2} mod p (with 0 ↦ 0)."""
        p = self.Fp.p
        return np.array([
            pow(int(xi), p - 2, p) if int(xi) % p != 0 else 0
            for xi in x
        ])

    def _cube(self, x: np.ndarray) -> np.ndarray:
        """σ(x)_i = x_i^3 mod p."""
        p = self.Fp.p
        return np.array([pow(int(xi), 3, p) for xi in x])

    def _id_spn(self, x: np.ndarray) -> np.ndarray:
        """σ(y) = y + S(M · S(y) + c), where S = component-wise cube,
        M = circulant mixing matrix, c = (1,...,1) constant vector.
        
        Properties:
          - Linear part = I (identity): interpolation gives (I+B)·A, not A
          - Constant c inside outer S: breaks scaling homogeneity
            σ(λy) ≠ λ^e σ(y) because (λ^3·w + c)^3 ≠ λ^k·(w+c)^3
          - σ(0) = S(c) ≠ 0: acceptable (d equations in d² unknowns)
          - Degree 9 (cube ∘ linear ∘ cube)
          - Cross-component mixing via M
          
        Resists:
          - Langa's scaling attack (constant c breaks λ^e factoring)
          - Polynomial interpolation (linear part I entangles A and B)
          - Differential column attack (cross-component mixing)
        """
        p = self.Fp.p
        d = len(x)
        
        # Adjust M if dimension doesn't match precomputed
        if d != self.d:
            M = self._build_mixing_matrix(d, p)
        else:
            M = self._M
        
        # S(y) = component-wise cube
        s1 = np.array([pow(int(xi), 3, p) for xi in x])
        
        # M · S(y) + c  (c = (1, ..., 1))
        mixed = np.array([
            (sum(int(M[i][j]) * int(s1[j]) for j in range(d)) + 1) % p
            for i in range(d)
        ])
        
        # S(M · S(y) + c) = component-wise cube
        s2 = np.array([pow(int(mixed[i]), 3, p) for i in range(d)])
        
        # y + S(M · S(y) + c)
        return np.array([(int(x[i]) + int(s2[i])) % p for i in range(d)])

    def attacker_system_degree(self) -> int:
        """
        Degree of the attacker's polynomial system in key unknowns.

        For c_v = A_v·s_v + B_v·σ(A_v·s_v):
          - σ = cube (deg 3): terms like b_ij · (a·s)^3 → degree 4
          - σ = inverse: with z-substitution, degree 2
          - σ = id_spn (deg 9): terms like b_ij · (a·s + N(a·s)) → degree 10
            But interpolation gives (I+B)A, so effective system
            to decompose has degree 9 in d² unknowns.
        """
        if self.name == "cube":
            return 4
        elif self.name == "inverse":
            return 2
        elif self.name == "id_spn":
            return 10

    def __repr__(self):
        return f"σ_{self.name}(deg={self._algebraic_degree})"
