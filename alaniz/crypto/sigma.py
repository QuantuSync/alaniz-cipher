"""
Public nonlinear maps σ: F_p^d → F_p^d for the Alaniz Cipher.

Supported maps:
  - inverse: σ(x)_i = x_i^{p-2}
             [DEPRECATED v1: vulnerable to scaling CPA attack]
  - cube:    σ(x)_i = x_i^3
             [DEPRECATED v1: vulnerable to scaling CPA attack]
  - id_spn:  σ(y) = y + S(M·S(y) + c), S = cube, c = (1,...,1)
             [DEPRECATED v2: broken by polynomial interpolation + Gröbner.
             Image coverage only ~60% (not a permutation of F_p^d).]
  - monomial_power: σ(y) = y + ι^{-1}(π_e(L·ι(y) + 1)), π_e(x) = x^e in F_{p^d}
             [RECOMMENDED v3: vector-valued, defeats interpolation
             under probabilistic encryption, APN or nearly-APN properties.]

The v3 map `monomial_power` operates in the multiplicative structure of
F_{p^d}, mixing all d coordinates through field multiplication. This
defeats the scaling and interpolation attacks that broke v1 and v2.
"""

from __future__ import annotations
import numpy as np
from ..core.field import FiniteField


class Sigma:
    """Public nonlinear map for the encryption scheme."""

    def __init__(self, name: str, Fp: FiniteField, d: int = 2,
                 field_ext=None, L_gf=None, exponent: int = None):
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
            self._M = self._build_mixing_matrix(d, Fp.p)
        elif name == "monomial_power":
            # v3 sigma: vectorial in F_{p^d}
            if field_ext is None:
                from .field_ext import ExtensionField
                field_ext = ExtensionField(Fp.p, d)
            self.field_ext = field_ext
            if exponent is None:
                from .field_ext import ExtensionField
                exponent = ExtensionField.find_secure_exponent(Fp.p, d)
            self.exponent = exponent
            if L_gf is None:
                L_gf = field_ext.random_nonzero()
            self.L_gf = L_gf
            self.c_gf = field_ext.GF(1)
            self._forward = self._monomial_power
            self._algebraic_degree = exponent
        else:
            raise ValueError(
                f"Unknown sigma: {name}. "
                f"Use 'inverse', 'cube', 'id_spn', or 'monomial_power'."
            )

    def _build_mixing_matrix(self, d: int, p: int) -> np.ndarray:
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

    # ---------------- v1 sigmas (deprecated) ----------------

    def _inverse(self, x: np.ndarray) -> np.ndarray:
        p = self.Fp.p
        return np.array([
            pow(int(xi), p - 2, p) if int(xi) % p != 0 else 0
            for xi in x
        ])

    def _cube(self, x: np.ndarray) -> np.ndarray:
        p = self.Fp.p
        return np.array([pow(int(xi), 3, p) for xi in x])

    # ---------------- v2 sigma (deprecated) ----------------

    def _id_spn(self, x: np.ndarray) -> np.ndarray:
        """σ_{SPN}(y) = y + S(M·S(y) + c), component-wise cube + mixing."""
        p = self.Fp.p
        d = len(x)
        M = self._M if d == self.d else self._build_mixing_matrix(d, p)

        s1 = np.array([pow(int(xi), 3, p) for xi in x])
        mixed = np.array([
            (sum(int(M[i][j]) * int(s1[j]) for j in range(d)) + 1) % p
            for i in range(d)
        ])
        s2 = np.array([pow(int(mixed[i]), 3, p) for i in range(d)])
        return np.array([(int(x[i]) + int(s2[i])) % p for i in range(d)])

    # ---------------- v3 sigma (RECOMMENDED) ----------------

    def _monomial_power(self, y: np.ndarray) -> np.ndarray:
        """
        v3 σ: σ(y) = y + ι^{-1}(π_e(L·ι(y) + 1))
        where π_e(x) = x^e in F_{p^d}.

        Vector-valued: mixes all d coordinates through field multiplication.
        Defeats scaling (Langa) and interpolation (v2 CPA) attacks.
        """
        y_gf = self.field_ext.vec_to_gf(y)
        inner = self.L_gf * y_gf + self.c_gf
        powered = inner ** self.exponent
        out_gf = y_gf + powered  # in F_{p^d}
        # Extract the "nonlinear tail" ι^{-1}(π_e(L·ι(y)+1)) and add to y
        tail = self.field_ext.gf_to_vec(powered)
        return np.array([(int(y[i]) + int(tail[i])) % self.Fp.p
                         for i in range(self.d)])

    # ---------------- metadata ----------------

    def attacker_system_degree(self) -> int:
        """Degree of the attacker's polynomial system in key unknowns."""
        if self.name == "cube":
            return 4
        elif self.name == "inverse":
            return 2
        elif self.name == "id_spn":
            return 10
        elif self.name == "monomial_power":
            return self.exponent + 1
        return -1

    def __repr__(self):
        if self.name == "monomial_power":
            return f"σ_monomial(e={self.exponent}, F_{{{self.Fp.p}^{self.d}}})"
        return f"σ_{self.name}(deg={self._algebraic_degree})"
