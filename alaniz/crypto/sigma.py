"""
Public nonlinear maps σ: F_p^d → F_p^d for the Protocolo Alaniz.

The map σ is applied AFTER the secret linear transformation A_v·s_v,
creating key-message entanglement that defeats CPA linearization.

Supported maps:
  - inverse: σ(x)_i = x_i^{p-2} (multiplicative inverse, 0↦0)
  - cube:    σ(x)_i = x_i^3 (requires gcd(3, p-1) = 1)
"""

from __future__ import annotations
import numpy as np
from typing import Callable
from ..core.field import FiniteField


class Sigma:
    """Public nonlinear map for the encryption scheme."""

    def __init__(self, name: str, Fp: FiniteField):
        self.name = name
        self.Fp = Fp

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
        else:
            raise ValueError(f"Unknown sigma: {name}. Use 'inverse' or 'cube'.")

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

    def attacker_system_degree(self) -> int:
        """
        Degree of the attacker's polynomial system in key unknowns.

        For c_v = A_v·s_v + B_v·σ(A_v·s_v):
          - σ = cube (deg 3): terms like b_ij · (a·s)^3 → degree 4
          - σ = inverse: with z-substitution, degree 2
        """
        if self.name == "cube":
            return 4  # degree 3 of σ × degree 1 of B_v
        elif self.name == "inverse":
            return 2  # after z = y^{-1} substitution

    def __repr__(self):
        return f"σ_{self.name}(deg={self._algebraic_degree})"
