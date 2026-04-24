"""
Alaniz Cipher protocols.

This module provides two protocol classes:
  - Protocol (v2, DEPRECATED): deterministic scheme with σ_SPN, matrices
    A_v, B_v. Broken by interpolation + Gröbner attack (see paper v3,
    Section 3). Retained for historical reproducibility.
  - ProtocolV3 (RECOMMENDED): probabilistic scheme with vector-valued
    σ in F_{p^d}, scalar β_v ∈ F_{p^d}^*, polynomial-time decryption
    via Cantor-Zassenhaus factorization.

For new code, use ProtocolV3. For reproducing v2 results, use Protocol.
"""

from __future__ import annotations
import numpy as np
import random as _random
import warnings
from dataclasses import dataclass, field as _field
from typing import Optional
from itertools import product

from ..core.field import FiniteField
from ..core.sheaf import Sheaf
from .sigma import Sigma


# ============================================================================
# v2 Protocol (DEPRECATED, retained for historical reproducibility)
# ============================================================================

@dataclass
class KeyPair:
    """v2 secret key: invertible matrices (A_v, B_v) per node."""
    A: dict
    B: dict


@dataclass
class PublicParams:
    """v2 public parameters."""
    sheaf: Sheaf
    sigma: Sigma

    @classmethod
    def generate(cls, sheaf: Sheaf, sigma_name: str = "id_spn") -> "PublicParams":
        return cls(sheaf=sheaf,
                   sigma=Sigma(sigma_name, sheaf.Fp, d=sheaf.dv))


class Protocol:
    """
    Alaniz v2 Protocol (DEPRECATED).

    Encryption: c_v = A_v·s_v + B_v·σ(A_v·s_v), deterministic.
    Known broken by polynomial interpolation + Gröbner basis attack.
    See paper v3 Section 3 and experiments/14_v2_redteam.py for details.
    """

    def __init__(self, params: PublicParams):
        warnings.warn(
            "Protocol (v2) is deprecated: it is broken by the interpolation "
            "+ Gröbner attack of paper v3. Use ProtocolV3 for new code.",
            DeprecationWarning, stacklevel=2,
        )
        self.params = params
        self.sheaf = params.sheaf
        self.Fp = params.sheaf.Fp
        self.sigma = params.sigma
        self.dv = params.sheaf.dv

    def keygen(self, rng=None) -> KeyPair:
        r = rng or _random
        A = {v: self.Fp.random_gl(self.dv, rng=r)
             for v in self.sheaf.graph.nodes}
        B = {v: self.Fp.random_gl(self.dv, rng=r)
             for v in self.sheaf.graph.nodes}
        return KeyPair(A=A, B=B)

    def encrypt_node(self, sv: np.ndarray, Av: np.ndarray,
                     Bv: np.ndarray) -> np.ndarray:
        y = self.Fp.mat_vec(Av, sv)
        z = self.sigma(y)
        t = self.Fp.mat_vec(Bv, z)
        return np.array([
            (int(y[i]) + int(t[i])) % self.Fp.p
            for i in range(self.dv)
        ])

    def encrypt(self, s: np.ndarray, key: KeyPair) -> np.ndarray:
        if not self.sheaf.is_global_section(s):
            raise ValueError("Plaintext must be a global section (δ^0(s) = 0)")
        c = np.zeros(self.sheaf.C0_dim, dtype=int)
        for v in self.sheaf.graph.nodes:
            sv = self.sheaf.get_node_value(s, v)
            cv = self.encrypt_node(sv, key.A[v], key.B[v])
            self.sheaf.set_node_value(c, v, cv)
        return c

    def _solve_y_at_node(self, cv: np.ndarray, Bv: np.ndarray) -> list:
        p, d = self.Fp.p, self.dv
        solutions = []
        if d == 2:
            for y0 in range(p):
                for y1 in range(p):
                    y = np.array([y0, y1])
                    z = self.sigma(y)
                    t = self.Fp.mat_vec(Bv, z)
                    res = np.array([(int(y[i]) + int(t[i])) % p for i in range(d)])
                    if np.array_equal(res, cv):
                        solutions.append(y)
        elif d <= 4 and p <= 31:
            for combo in product(range(p), repeat=d):
                y = np.array(combo)
                z = self.sigma(y)
                t = self.Fp.mat_vec(Bv, z)
                res = np.array([(int(y[i]) + int(t[i])) % p for i in range(d)])
                if np.array_equal(res, cv):
                    solutions.append(y)
        else:
            raise NotImplementedError(
                f"v2 brute-force solver not implemented for d={d}, p={p}."
            )
        return solutions

    def decrypt_tree(self, c: np.ndarray, key: KeyPair,
                     root: int = 0) -> Optional[np.ndarray]:
        if not self.sheaf.graph.is_tree:
            raise ValueError("Tree decryption requires a tree graph")
        R_to = self.sheaf.tree_propagation_maps(root)
        c_root = self.sheaf.get_node_value(c, root)
        y_solutions = self._solve_y_at_node(c_root, key.B[root])
        A_root_inv = self.Fp.mat_inv(key.A[root])

        for y_r in y_solutions:
            x_r = self.Fp.mat_vec(A_root_inv, y_r)
            s_candidate = np.zeros(self.sheaf.C0_dim, dtype=int)
            for v in self.sheaf.graph.nodes:
                s_v = self.Fp.mat_vec(R_to[v], x_r)
                self.sheaf.set_node_value(s_candidate, v, s_v)
            all_ok = True
            for v in self.sheaf.graph.nodes:
                sv = self.sheaf.get_node_value(s_candidate, v)
                cv_check = self.encrypt_node(sv, key.A[v], key.B[v])
                cv_actual = self.sheaf.get_node_value(c, v)
                if not np.array_equal(cv_check, cv_actual):
                    all_ok = False
                    break
            if all_ok:
                return s_candidate
        return None

    def decrypt_bruteforce(self, c: np.ndarray, key: KeyPair) -> list:
        local_sols = {}
        for v in self.sheaf.graph.nodes:
            cv = self.sheaf.get_node_value(c, v)
            y_sols = self._solve_y_at_node(cv, key.B[v])
            A_inv = self.Fp.mat_inv(key.A[v])
            local_sols[v] = [self.Fp.mat_vec(A_inv, y) for y in y_sols]
        valid = []
        nodes = list(self.sheaf.graph.nodes)
        for combo in product(*[local_sols[v] for v in nodes]):
            cand = np.zeros(self.sheaf.C0_dim, dtype=int)
            for v, sol in zip(nodes, combo):
                self.sheaf.set_node_value(cand, v, sol)
            if self.sheaf.is_global_section(cand):
                valid.append(cand)
        return valid

    def encode(self, message: bytes) -> np.ndarray:
        k = self.sheaf.H0_dim
        p = self.Fp.p
        bits_per_coeff = (p - 1).bit_length()
        capacity_bits = k * bits_per_coeff
        if len(message) * 8 > capacity_bits:
            raise ValueError(f"Message too long: {len(message)*8} bits > {capacity_bits}")
        msg_int = int.from_bytes(message, "big")
        coeffs = []
        for _ in range(k):
            coeffs.append(msg_int % p)
            msg_int //= p
        return self.sheaf.section_from_coeffs(coeffs)

    def decode(self, s: np.ndarray) -> bytes:
        k = self.sheaf.H0_dim
        p = self.Fp.p
        basis_mat = np.zeros((self.sheaf.C0_dim, k), dtype=object)
        for i, bv in enumerate(self.sheaf.H0_basis):
            basis_mat[:, i] = bv.astype(object)
        coeffs = self._solve_basis_coeffs(basis_mat, s)
        msg_int = 0
        for i in range(k - 1, -1, -1):
            msg_int = msg_int * p + int(coeffs[i])
        n_bytes = (msg_int.bit_length() + 7) // 8
        return msg_int.to_bytes(max(n_bytes, 1), "big")

    def _solve_basis_coeffs(self, basis_mat: np.ndarray, s: np.ndarray) -> list:
        k = basis_mat.shape[1]
        p = self.Fp.p
        M = basis_mat.astype(object) % p
        s_obj = s.astype(object) % p
        sub_M = M[:k, :]
        sub_s = s_obj[:k]
        aug = np.zeros((k, k + 1), dtype=object)
        aug[:, :k] = sub_M % p
        aug[:, k] = sub_s % p
        for col in range(k):
            pivot = -1
            for r in range(col, k):
                if aug[r, col] % p != 0:
                    pivot = r; break
            if pivot == -1:
                raise ValueError("Basis matrix degenerate at selected rows")
            aug[[col, pivot]] = aug[[pivot, col]]
            iv = pow(int(aug[col, col]), p - 2, p)
            aug[col] = (aug[col] * iv) % p
            for r in range(k):
                if r != col and aug[r, col] % p != 0:
                    aug[r] = (aug[r] - aug[r, col] * aug[col]) % p
        return [int(aug[i, k]) % p for i in range(k)]


# ============================================================================
# v3 Protocol (RECOMMENDED)
# ============================================================================

@dataclass
class KeyPairV3:
    """v3 secret key: matrices A_v ∈ GL(d, F_p) and scalars β_v ∈ F_{p^d}^*."""
    A: dict           # v -> np.ndarray (d×d matrix in F_p)
    beta: dict        # v -> field_ext.GF element
    s_reject: bytes   # for FO implicit rejection


@dataclass
class PublicParamsV3:
    """v3 public parameters."""
    sheaf: Sheaf
    sigma: Sigma
    field_ext: object
    L_gf: object
    exponent: int

    @classmethod
    def generate(cls, sheaf: Sheaf, exponent: int = None,
                 L_gf=None) -> "PublicParamsV3":
        from ..core.field_ext import ExtensionField
        p, d = sheaf.Fp.p, sheaf.dv
        field_ext = ExtensionField(p, d)
        if exponent is None:
            exponent = ExtensionField.find_secure_exponent(p, d)
        if L_gf is None:
            L_gf = field_ext.random_nonzero()
        sigma = Sigma("monomial_power", sheaf.Fp, d=d,
                      field_ext=field_ext, L_gf=L_gf, exponent=exponent)
        return cls(sheaf=sheaf, sigma=sigma, field_ext=field_ext,
                   L_gf=L_gf, exponent=exponent)


class ProtocolV3:
    """
    Alaniz v3 Protocol (RECOMMENDED).

    Encryption (per node, probabilistic):
      r_v = PRG(nonce, v)
      u_v = ι(A_v·s_v + r_v)   ∈ F_{p^d}
      w_v = β_v·u_v + (β_v-1)·(L·u_v + 1)^e   ∈ F_{p^d}
      c_v = ι^{-1}(w_v) - r_v

    Decryption: solves γ·τ^e + α·τ - (c' + α) = 0 in F_{p^d}[τ] via
    Cantor-Zassenhaus factorization, then verifies through the cohomology
    constraint at all n-1 non-root vertices.

    Target: IND-CPA as PKE; IND-CCA after the FO transform (see kem.py).
    """

    def __init__(self, params: PublicParamsV3):
        self.params = params
        self.sheaf = params.sheaf
        self.Fp = params.sheaf.Fp
        self.sigma = params.sigma
        self.field_ext = params.field_ext
        self.L_gf = params.L_gf
        self.exponent = params.exponent
        self.dv = params.sheaf.dv

    def keygen(self, rng=None, root: int = 0) -> KeyPairV3:
        """
        Generate a v3 secret key.

        Samples:
          - A_v ∈ GL(d, F_p) uniform per node
          - β_v ∈ F_{p^d}^* \\ {1} uniform per node
          - s_reject: 32 random bytes for FO implicit rejection
        """
        import os
        r = rng or _random
        A = {v: self.Fp.random_gl(self.dv, rng=r)
             for v in self.sheaf.graph.nodes}
        GF = self.field_ext.GF
        beta = {
            v: self.field_ext.random_not_in([0, 1], rng=r)
            for v in self.sheaf.graph.nodes
        }
        s_reject = os.urandom(32)
        return KeyPairV3(A=A, beta=beta, s_reject=s_reject)

    # ---- Encryption ----

    def encrypt_node(self, sv: np.ndarray, Av: np.ndarray,
                     beta_v, r_v: np.ndarray) -> np.ndarray:
        """Single-node encryption. Returns c_v ∈ F_p^d."""
        y = self.Fp.mat_vec(Av, sv)
        a = np.array([(int(y[i]) + int(r_v[i])) % self.Fp.p
                      for i in range(self.dv)])
        u_gf = self.field_ext.vec_to_gf(a)
        t_gf = self.L_gf * u_gf + self.field_ext.GF(1)
        powered = t_gf ** self.exponent
        w_gf = beta_v * u_gf + (beta_v - self.field_ext.GF(1)) * powered
        w_vec = self.field_ext.gf_to_vec(w_gf)
        return np.array([(int(w_vec[i]) - int(r_v[i])) % self.Fp.p
                         for i in range(self.dv)])

    def encrypt(self, s: np.ndarray, key: KeyPairV3,
                nonce: bytes = None) -> tuple:
        """
        Encrypt a global section s ∈ H^0.

        Returns (nonce, c) where c ∈ ⊕_v F_p^{d}.
        """
        import os
        from .prg import prg_derive
        if nonce is None:
            nonce = os.urandom(16)
        if not self.sheaf.is_global_section(s):
            raise ValueError("Plaintext must be a global section (δ^0(s) = 0)")
        c = np.zeros(self.sheaf.C0_dim, dtype=int)
        for v in self.sheaf.graph.nodes:
            sv = self.sheaf.get_node_value(s, v)
            r_v = prg_derive(nonce, v, self.dv, self.Fp.p)
            cv = self.encrypt_node(sv, key.A[v], key.beta[v], r_v)
            self.sheaf.set_node_value(c, v, cv)
        return nonce, c

    # ---- Decryption ----

    def decrypt(self, nonce: bytes, c: np.ndarray, key: KeyPairV3,
                root: int = 0) -> Optional[np.ndarray]:
        """
        Decrypt (nonce, c) via univariate factorization in F_{p^d}[τ].

        1. At root: form F(τ) = γ·τ^e + α·τ - (c' + α), find roots.
        2. For each root: recover candidate s, propagate through the sheaf.
        3. Verify encryption consistency at all vertices; return the
           unique surviving candidate.
        """
        from .prg import prg_derive
        GF = self.field_ext.GF

        r_root = prg_derive(nonce, root, self.dv, self.Fp.p)
        c_root = self.sheaf.get_node_value(c, root)
        c_root_with_r = np.array([
            (int(c_root[i]) + int(r_root[i])) % self.Fp.p
            for i in range(self.dv)
        ])
        c_prime = self.field_ext.vec_to_gf(c_root_with_r)

        beta = key.beta[root]
        alpha = beta * (self.L_gf ** -1)
        gamma = beta - GF(1)

        # Build F(τ) = γ·τ^e + α·τ - (c' + α)
        coeffs_high_to_low = [GF(0)] * (self.exponent + 1)
        coeffs_high_to_low[0] = gamma
        coeffs_high_to_low[self.exponent - 1] = alpha
        coeffs_high_to_low[self.exponent] = -(c_prime + alpha)

        roots = self.field_ext.find_roots(coeffs_high_to_low)

        R_to = self.sheaf.tree_propagation_maps(root)
        A_root_inv = self.Fp.mat_inv(key.A[root])

        for tau in roots:
            u_gf = (tau - GF(1)) * (self.L_gf ** -1)
            a = self.field_ext.gf_to_vec(u_gf)
            y_r = np.array([(int(a[i]) - int(r_root[i])) % self.Fp.p
                            for i in range(self.dv)])
            s_cand_root = self.Fp.mat_vec(A_root_inv, y_r)

            # Propagate and verify
            s_candidate = np.zeros(self.sheaf.C0_dim, dtype=int)
            for v in self.sheaf.graph.nodes:
                s_v = self.Fp.mat_vec(R_to[v], s_cand_root)
                self.sheaf.set_node_value(s_candidate, v, s_v)

            all_ok = True
            for v in self.sheaf.graph.nodes:
                sv = self.sheaf.get_node_value(s_candidate, v)
                r_v = prg_derive(nonce, v, self.dv, self.Fp.p)
                cv_check = self.encrypt_node(sv, key.A[v], key.beta[v], r_v)
                cv_actual = self.sheaf.get_node_value(c, v)
                if not np.array_equal(cv_check, cv_actual):
                    all_ok = False
                    break

            if all_ok:
                return s_candidate
        return None

    # ---- Encoding (shared with v2) ----

    def encode(self, message: bytes) -> np.ndarray:
        k = self.sheaf.H0_dim
        p = self.Fp.p
        bits_per_coeff = (p - 1).bit_length()
        capacity_bits = k * bits_per_coeff
        if len(message) * 8 > capacity_bits:
            raise ValueError(
                f"Message too long: {len(message)*8} bits > {capacity_bits}"
            )
        msg_int = int.from_bytes(message, "big")
        coeffs = []
        for _ in range(k):
            coeffs.append(msg_int % p)
            msg_int //= p
        return self.sheaf.section_from_coeffs(coeffs)

    def decode(self, s: np.ndarray) -> bytes:
        k = self.sheaf.H0_dim
        p = self.Fp.p
        basis_mat = np.zeros((self.sheaf.C0_dim, k), dtype=object)
        for i, bv in enumerate(self.sheaf.H0_basis):
            basis_mat[:, i] = bv.astype(object)
        # Use the v2 solver
        proto_v2 = Protocol.__new__(Protocol)
        proto_v2.sheaf = self.sheaf
        proto_v2.Fp = self.Fp
        coeffs = proto_v2._solve_basis_coeffs(basis_mat, s)
        msg_int = 0
        for i in range(k - 1, -1, -1):
            msg_int = msg_int * p + int(coeffs[i])
        n_bytes = (msg_int.bit_length() + 7) // 8
        return msg_int.to_bytes(max(n_bytes, 1), "big")
