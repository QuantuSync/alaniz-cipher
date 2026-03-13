"""
Protocolo Alaniz v1.0 — SHEC (Sheaf Homology Encryption Cipher)

Encryption:
  c_v = A_v · s_v + B_v · σ(A_v · s_v)

where s ∈ H^0(G, F_0) is a global section, (A_v, B_v) are secret
GL(d, F_p) matrices per node, and σ is a public nonlinear map.

Decryption:
  1. Solve at root: find y_r s.t. y + B_r·σ(y) = c_r
  2. For each solution: x_r = A_r^{-1}·y_r, propagate to all nodes
  3. Verify: the true plaintext passes encryption check at all nodes

Correctness: Dec(Enc(s)) = s with probability 1 - L·γ^{n-1}
where L ≈ 2 is the root solution count and γ ≈ 1/p^d.
"""

from __future__ import annotations
import numpy as np
import random as _random
from dataclasses import dataclass, field
from typing import Optional
from itertools import product

from ..core.field import FiniteField
from ..core.sheaf import Sheaf
from .sigma import Sigma


@dataclass
class KeyPair:
    """Secret key: invertible matrices (A_v, B_v) per node."""
    A: dict[int, np.ndarray]
    B: dict[int, np.ndarray]


@dataclass
class PublicParams:
    """Public parameters of the protocol."""
    sheaf: Sheaf
    sigma: Sigma

    @classmethod
    def generate(cls, sheaf: Sheaf, sigma_name: str = "id_spn") -> PublicParams:
        return cls(sheaf=sheaf,
                   sigma=Sigma(sigma_name, sheaf.Fp, d=sheaf.dv))


class Protocol:
    """Protocolo Alaniz v1.0 — Full encryption scheme."""

    def __init__(self, params: PublicParams):
        self.params = params
        self.sheaf = params.sheaf
        self.Fp = params.sheaf.Fp
        self.sigma = params.sigma
        self.dv = params.sheaf.dv

    # --- Key Generation ---

    def keygen(self, rng=None) -> KeyPair:
        """Generate secret key: (A_v, B_v) ∈ GL(d, F_p)² per node."""
        r = rng or _random
        A = {v: self.Fp.random_gl(self.dv, rng=r)
             for v in self.sheaf.graph.nodes}
        B = {v: self.Fp.random_gl(self.dv, rng=r)
             for v in self.sheaf.graph.nodes}
        return KeyPair(A=A, B=B)

    # --- Encryption ---

    def encrypt_node(self, sv: np.ndarray, Av: np.ndarray,
                     Bv: np.ndarray) -> np.ndarray:
        """
        Encrypt at a single node:
          c_v = A_v · s_v + B_v · σ(A_v · s_v)
        """
        y = self.Fp.mat_vec(Av, sv)
        z = self.sigma(y)
        t = self.Fp.mat_vec(Bv, z)
        return np.array([
            (int(y[i]) + int(t[i])) % self.Fp.p
            for i in range(self.dv)
        ])

    def encrypt(self, s: np.ndarray, key: KeyPair) -> np.ndarray:
        """
        Encrypt a global section s ∈ H^0(G, F_0).

        Returns ciphertext c ∈ ⊕_v F_p^{d_v}.
        """
        if not self.sheaf.is_global_section(s):
            raise ValueError("Plaintext must be a global section (δ^0(s) = 0)")

        c = np.zeros(self.sheaf.C0_dim, dtype=int)
        for v in self.sheaf.graph.nodes:
            sv = self.sheaf.get_node_value(s, v)
            cv = self.encrypt_node(sv, key.A[v], key.B[v])
            self.sheaf.set_node_value(c, v, cv)
        return c

    # --- Decryption ---

    def _solve_y_at_node(self, cv: np.ndarray,
                         Bv: np.ndarray) -> list[np.ndarray]:
        """
        Find all y ∈ F_p^d such that y + B_v · σ(y) = c_v.

        Brute force over F_p^d (feasible for small p, d).
        """
        p, d = self.Fp.p, self.dv
        solutions = []

        if d == 2:
            for y0 in range(p):
                for y1 in range(p):
                    y = np.array([y0, y1])
                    z = self.sigma(y)
                    t = self.Fp.mat_vec(Bv, z)
                    res = np.array([
                        (int(y[i]) + int(t[i])) % p for i in range(d)
                    ])
                    if np.array_equal(res, cv):
                        solutions.append(y)
        elif d <= 4 and p <= 31:
            # Brute force in y-space for small parameters
            for combo in product(range(p), repeat=d):
                y = np.array(combo)
                z = self.sigma(y)
                t = self.Fp.mat_vec(Bv, z)
                res = np.array([
                    (int(y[i]) + int(t[i])) % p for i in range(d)
                ])
                if np.array_equal(res, cv):
                    solutions.append(y)
        else:
            raise NotImplementedError(
                f"Decryption solver not implemented for d={d}, p={p}. "
                f"Needs Gröbner basis or Newton-Hensel lifting."
            )
        return solutions

    def decrypt_tree(self, c: np.ndarray, key: KeyPair,
                     root: int = 0) -> Optional[np.ndarray]:
        """
        Decrypt using tree propagation (efficient algorithm).

        1. Solve y + B_r·σ(y) = c_r at the root
        2. For each y_r: compute x_r = A_r^{-1}·y_r
        3. Propagate x_r to all nodes via restriction maps
        4. Verify encryption at all non-root nodes

        Complexity: O(L · n · poly(d, log p))
        where L ≈ 2 is the root solution count.

        Returns the unique valid plaintext, or None if decryption fails.
        """
        if not self.sheaf.graph.is_tree:
            raise ValueError("Tree decryption requires a tree graph")

        R_to = self.sheaf.tree_propagation_maps(root)
        c_root = self.sheaf.get_node_value(c, root)

        # Step 1: solve at root in y-space
        y_solutions = self._solve_y_at_node(c_root, key.B[root])

        # Step 2-4: for each root solution, propagate and verify
        A_root_inv = self.Fp.mat_inv(key.A[root])

        for y_r in y_solutions:
            x_r = self.Fp.mat_vec(A_root_inv, y_r)

            # Propagate to all nodes
            s_candidate = np.zeros(self.sheaf.C0_dim, dtype=int)
            for v in self.sheaf.graph.nodes:
                s_v = self.Fp.mat_vec(R_to[v], x_r)
                self.sheaf.set_node_value(s_candidate, v, s_v)

            # Verify at all nodes
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

        return None  # No valid solution found (should not happen)

    def decrypt_bruteforce(self, c: np.ndarray,
                           key: KeyPair) -> list[np.ndarray]:
        """
        Decrypt via brute-force local solve + cohomological filter.

        Slower than tree propagation but works for any graph topology.
        Returns list of valid global sections mapping to c.
        """
        local_sols = {}
        for v in self.sheaf.graph.nodes:
            cv = self.sheaf.get_node_value(c, v)
            y_sols = self._solve_y_at_node(cv, key.B[v])
            A_inv = self.Fp.mat_inv(key.A[v])
            local_sols[v] = [self.Fp.mat_vec(A_inv, y) for y in y_sols]

        # Cohomological filter
        valid = []
        nodes = list(self.sheaf.graph.nodes)
        for combo in product(*[local_sols[v] for v in nodes]):
            cand = np.zeros(self.sheaf.C0_dim, dtype=int)
            for v, sol in zip(nodes, combo):
                self.sheaf.set_node_value(cand, v, sol)
            if self.sheaf.is_global_section(cand):
                valid.append(cand)
        return valid

    # --- Encoding / Decoding ---

    def encode(self, message: bytes) -> np.ndarray:
        """
        Encode a byte message as a global section in H^0.

        Maps message bytes to F_p coefficients in the H^0 basis.
        """
        k = self.sheaf.H0_dim
        p = self.Fp.p

        # Compute capacity
        bits_per_coeff = (p - 1).bit_length()  # conservative
        capacity_bits = k * bits_per_coeff
        if len(message) * 8 > capacity_bits:
            raise ValueError(
                f"Message too long: {len(message)*8} bits > "
                f"capacity {capacity_bits} bits "
                f"(H^0 dim={k}, p={p})"
            )

        # Convert bytes → integer → F_p coefficients
        msg_int = int.from_bytes(message, "big")
        coeffs = []
        for _ in range(k):
            coeffs.append(msg_int % p)
            msg_int //= p

        return self.sheaf.section_from_coeffs(coeffs)

    def decode(self, s: np.ndarray) -> bytes:
        """
        Decode a global section back to bytes.

        Extracts F_p coefficients from the H^0 basis representation.
        """
        k = self.sheaf.H0_dim
        p = self.Fp.p

        # Project onto H^0 basis to get coefficients
        # For tree graphs with our convention, s_root determines everything
        # and the basis vectors are linearly independent at the root
        # Use least-squares-like approach over F_p

        # Build basis matrix at root (or full)
        basis_mat = np.zeros((self.sheaf.C0_dim, k), dtype=object)
        for i, bv in enumerate(self.sheaf.H0_basis):
            basis_mat[:, i] = bv.astype(object)

        # Solve basis_mat · coeffs = s over F_p
        # Use the first k linearly independent rows
        coeffs = self._solve_basis_coeffs(basis_mat, s)

        # Convert coefficients → integer → bytes
        msg_int = 0
        for i in range(k - 1, -1, -1):
            msg_int = msg_int * p + int(coeffs[i])

        # Convert to bytes
        n_bytes = (msg_int.bit_length() + 7) // 8
        return msg_int.to_bytes(max(n_bytes, 1), "big")

    def _solve_basis_coeffs(self, basis_mat: np.ndarray,
                            s: np.ndarray) -> list[int]:
        """Solve basis_mat · c = s over F_p for coefficients c."""
        k = basis_mat.shape[1]
        p = self.Fp.p

        # Extract a k×k invertible submatrix
        M = basis_mat.astype(object) % p
        s_obj = s.astype(object) % p

        # Try first k rows
        sub_M = M[:k, :]
        sub_s = s_obj[:k]

        # Solve via row reduction
        aug = np.zeros((k, k + 1), dtype=object)
        aug[:, :k] = sub_M % p
        aug[:, k] = sub_s % p

        for col in range(k):
            pivot = -1
            for r in range(col, k):
                if aug[r, col] % p != 0:
                    pivot = r
                    break
            if pivot == -1:
                # Need different rows — fall back to scanning
                raise ValueError("Basis matrix degenerate at selected rows")
            aug[[col, pivot]] = aug[[pivot, col]]
            iv = pow(int(aug[col, col]), p - 2, p)
            aug[col] = (aug[col] * iv) % p
            for r in range(k):
                if r != col and aug[r, col] % p != 0:
                    aug[r] = (aug[r] - aug[r, col] * aug[col]) % p

        return [int(aug[i, k]) % p for i in range(k)]
