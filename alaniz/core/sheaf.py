"""
Cellular sheaves over graphs with coefficients in F_p.

A sheaf F on graph G assigns:
  - To each node v: a vector space F(v) ≅ F_p^{d_v}
  - To each edge e=(u,v): restriction maps ρ_{v←e}: F(e) → F(v)

The 0-th cohomology H^0(G, F) = ker(δ^0) captures global sections:
assignments of vectors to nodes consistent with all restriction maps.
"""

from __future__ import annotations
import numpy as np
import random as _random
from dataclasses import dataclass, field
from typing import Optional

from .field import FiniteField
from .graph import Graph


@dataclass
class Sheaf:
    """Cellular sheaf over a graph with uniform fiber dimension."""

    graph: Graph
    dv: int
    Fp: FiniteField
    rho: dict = field(default_factory=dict, repr=False)

    # Computed attributes
    delta0: np.ndarray = field(default=None, repr=False)
    H0_basis: list = field(default_factory=list, repr=False)

    @property
    def H0_dim(self) -> int:
        return len(self.H0_basis)

    @property
    def C0_dim(self) -> int:
        """Dimension of C^0 = ⊕_v F(v)."""
        return self.graph.n * self.dv

    @property
    def C1_dim(self) -> int:
        """Dimension of C^1 = ⊕_e F(e)."""
        return self.graph.m * self.dv

    @classmethod
    def random(cls, graph: Graph, dv: int, Fp: FiniteField,
               rng=None) -> Sheaf:
        """
        Construct a sheaf with random invertible restriction maps.

        Convention: for edge e=(u,v),
          ρ_{u←e} = I  (identity)
          ρ_{v←e} = R_e  (random GL matrix)

        This ensures dim H^0 = dv for tree graphs.
        """
        r = rng or _random
        sheaf = cls(graph=graph, dv=dv, Fp=Fp)

        for e in graph.edges:
            u, v = e
            sheaf.rho[(e, u)] = np.eye(dv, dtype=int)
            sheaf.rho[(e, v)] = Fp.random_gl(dv, rng=r)

        sheaf._build_coboundary()
        sheaf._compute_H0()
        return sheaf

    def _build_coboundary(self):
        """Build the coboundary operator δ^0: C^0 → C^1."""
        n, m, d = self.graph.n, self.graph.m, self.dv
        self.delta0 = np.zeros((m * d, n * d), dtype=int)

        for idx, e in enumerate(self.graph.edges):
            u, v = e
            for i in range(d):
                for j in range(d):
                    self.delta0[idx * d + i][v * d + j] = self.Fp.mod(
                        self.rho[(e, v)][i][j]
                    )
                    self.delta0[idx * d + i][u * d + j] = self.Fp.mod(
                        -self.rho[(e, u)][i][j]
                    )
        self.delta0 = self.Fp.mat_mod(self.delta0)

    def _compute_H0(self):
        """Compute H^0(G, F) = ker(δ^0) over F_p."""
        self.H0_basis = self.Fp.kernel(self.delta0)

    def random_section(self, rng=None) -> np.ndarray:
        """Sample a random element of H^0(G, F)."""
        r = rng or _random
        s = np.zeros(self.C0_dim, dtype=object)
        for bv in self.H0_basis:
            coeff = r.randint(0, self.Fp.p - 1)
            s = (s + coeff * bv.astype(object)) % self.Fp.p
        return s.astype(int)

    def section_from_coeffs(self, coeffs: list[int]) -> np.ndarray:
        """Build section from coefficients in the H^0 basis."""
        if len(coeffs) != self.H0_dim:
            raise ValueError(
                f"Need {self.H0_dim} coefficients, got {len(coeffs)}"
            )
        s = np.zeros(self.C0_dim, dtype=object)
        for c, bv in zip(coeffs, self.H0_basis):
            s = (s + c * bv.astype(object)) % self.Fp.p
        return s.astype(int)

    def is_global_section(self, s: np.ndarray) -> bool:
        """Check s ∈ H^0(G, F), i.e., δ^0(s) = 0."""
        check = self.Fp.mat_mod(
            self.delta0.astype(object) @ s.astype(object).reshape(-1, 1)
        )
        return bool(np.all(check == 0))

    def get_node_value(self, s: np.ndarray, v: int) -> np.ndarray:
        """Extract the fiber value s_v from a global section."""
        return s[v * self.dv : (v + 1) * self.dv].copy()

    def set_node_value(self, s: np.ndarray, v: int, val: np.ndarray):
        """Set the fiber value s_v in a section vector."""
        s[v * self.dv : (v + 1) * self.dv] = val

    def compose_restriction(self, path: list[int]) -> np.ndarray:
        """
        Compose restriction maps along a path in the graph.

        For tree graphs, this gives R_{root→v} such that
        s_v = R_{root→v} · s_root for any global section s.
        """
        R = np.eye(self.dv, dtype=int)
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            # Find the edge
            e = (u, v) if (u, v) in self.graph.edges else (v, u)
            # δ^0(s)(e) = ρ_{v←e}·s_v - ρ_{u←e}·s_u = 0
            # With our convention ρ_{u←e}=I: ρ_{v←e}·s_v = s_u
            # So s_v = ρ_{v←e}^{-1} · s_u
            rho_v = self.rho[(e, v)]
            rho_v_inv = self.Fp.mat_inv(rho_v)
            R = self.Fp.mat_mul(rho_v_inv, R)
        return R

    def tree_propagation_maps(self, root: int = 0) -> dict[int, np.ndarray]:
        """
        For tree graphs, compute R_{root→v} for all nodes.

        Returns dict mapping node → composed restriction matrix.
        """
        if not self.graph.is_tree:
            raise ValueError("Tree propagation requires a tree graph")

        adj = self.graph.adjacency_list()
        R_to = {root: np.eye(self.dv, dtype=int)}
        visited = {root}
        queue = [root]

        while queue:
            u = queue.pop(0)
            for v in adj[u]:
                if v in visited:
                    continue
                visited.add(v)
                queue.append(v)
                e = (u, v) if (u, v) in self.graph.edges else (v, u)
                rho_v = self.rho[(e, v)]
                rho_u = self.rho[(e, u)]
                # s_v from s_u: ρ_{v←e}·s_v = ρ_{u←e}·s_u
                # s_v = ρ_{v←e}^{-1} · ρ_{u←e} · s_u
                rho_v_inv = self.Fp.mat_inv(rho_v)
                transfer = self.Fp.mat_mul(rho_v_inv, rho_u)
                R_to[v] = self.Fp.mat_mul(transfer, R_to[u])

        return R_to

    def __repr__(self):
        return (
            f"Sheaf(graph={self.graph}, dv={self.dv}, "
            f"F={self.Fp}, dim H^0={self.H0_dim})"
        )
