"""
Cellular sheaves over graphs with coefficients in F_p.

A sheaf F on graph G assigns:
  - To each node v: a vector space F(v) ≅ F_p^{d_v}
  - To each edge e=(u,v): restriction maps ρ_{v←e}: F(e) → F(v)

The 0-th cohomology H^0(G, F) = ker(δ^0) captures global sections:
assignments of vectors to nodes consistent with all restriction maps.

For v3, constructor `random_with_cohomology` samples ρ_e freely on a
spanning tree and computes the remaining edges from the cycle constraints
∏_{e ∈ γ} ρ_e^{±1} = I, ensuring dim H^0 = d_v even when β_1 >= 1.
"""

from __future__ import annotations
import numpy as np
import random as _random
from dataclasses import dataclass, field

from .field import FiniteField
from .graph import Graph


@dataclass
class Sheaf:
    """Cellular sheaf over a graph with uniform fiber dimension."""

    graph: Graph
    dv: int
    Fp: FiniteField
    rho: dict = field(default_factory=dict, repr=False)

    delta0: np.ndarray = field(default=None, repr=False)
    H0_basis: list = field(default_factory=list, repr=False)

    @property
    def H0_dim(self) -> int:
        return len(self.H0_basis)

    @property
    def C0_dim(self) -> int:
        return self.graph.n * self.dv

    @property
    def C1_dim(self) -> int:
        return self.graph.m * self.dv

    @classmethod
    def random(cls, graph: Graph, dv: int, Fp: FiniteField,
               rng=None) -> "Sheaf":
        """
        Construct a sheaf with random invertible restriction maps.

        For tree graphs: ensures dim H^0 = dv.
        For graphs with β_1 >= 1: may yield dim H^0 < dv (generic restrictions
        do not satisfy cycle constraints). Use `random_with_cohomology`
        for graphs with cycles when dim H^0 = dv is required.
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

    @classmethod
    def random_with_cohomology(cls, graph: Graph, dv: int, Fp: FiniteField,
                                root: int = 0, rng=None) -> "Sheaf":
        """
        Sample random restrictions satisfying cycle constraints.

        Procedure:
          1. Compute spanning tree T from `root`.
          2. Sample ρ_e ∈ GL(d, F_p) uniformly for each tree edge.
          3. For each non-tree edge e = (u, v): compute the propagated
             restriction from u to v through T and set ρ_e to match,
             ensuring ∏ ρ = I on the fundamental cycle through e.

        Guarantees dim H^0 = dv regardless of β_1.
        """
        r = rng or _random
        sheaf = cls(graph=graph, dv=dv, Fp=Fp)

        tree_edges, non_tree_edges = graph.spanning_tree(root)
        tree_set = set(tree_edges)

        # Step 1: sample tree edges freely
        for e in tree_edges:
            u, v = e
            sheaf.rho[(e, u)] = np.eye(dv, dtype=int)
            sheaf.rho[(e, v)] = Fp.random_gl(dv, rng=r)

        # Step 2: compute propagation matrices R[root -> v] on the tree
        R_to = sheaf._compute_tree_propagation(tree_edges, root)

        # Step 3: constrain non-tree edges
        for e in non_tree_edges:
            u, v = e
            # A section s with s_root = x has s_u = R_to[u] x and s_v = R_to[v] x.
            # Edge constraint: ρ_{v←e} s_v = ρ_{u←e} s_u for all x.
            # With ρ_{u←e} = I: ρ_{v←e} R_to[v] = R_to[u]
            #   => ρ_{v←e} = R_to[u] R_to[v]^{-1}
            sheaf.rho[(e, u)] = np.eye(dv, dtype=int)
            R_v_inv = Fp.mat_inv(R_to[v])
            sheaf.rho[(e, v)] = Fp.mat_mul(R_to[u], R_v_inv)

        sheaf._build_coboundary()
        sheaf._compute_H0()
        return sheaf

    def _compute_tree_propagation(self, tree_edges: tuple,
                                   root: int) -> dict:
        """BFS on the spanning tree to compute R[root -> v] for all v."""
        Fp = self.Fp
        dv = self.dv
        adj = {v: [] for v in self.graph.nodes}
        for u, v in tree_edges:
            adj[u].append(v)
            adj[v].append(u)
        R_to = {root: np.eye(dv, dtype=int)}
        visited = {root}
        queue = [root]
        while queue:
            u = queue.pop(0)
            for v in adj[u]:
                if v in visited:
                    continue
                visited.add(v)
                queue.append(v)
                e = (min(u, v), max(u, v))
                # s_v = ρ_{v←e}^{-1} · ρ_{u←e} · s_u (both sampled above)
                rho_v = self.rho[(e, v)]
                rho_u = self.rho[(e, u)]
                transfer = Fp.mat_mul(Fp.mat_inv(rho_v), rho_u)
                R_to[v] = Fp.mat_mul(transfer, R_to[u])
        return R_to

    def _build_coboundary(self):
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
        self.H0_basis = self.Fp.kernel(self.delta0)

    def random_section(self, rng=None) -> np.ndarray:
        r = rng or _random
        s = np.zeros(self.C0_dim, dtype=object)
        for bv in self.H0_basis:
            coeff = r.randint(0, self.Fp.p - 1)
            s = (s + coeff * bv.astype(object)) % self.Fp.p
        return s.astype(int)

    def section_from_coeffs(self, coeffs) -> np.ndarray:
        if len(coeffs) != self.H0_dim:
            raise ValueError(
                f"Need {self.H0_dim} coefficients, got {len(coeffs)}"
            )
        s = np.zeros(self.C0_dim, dtype=object)
        for c, bv in zip(coeffs, self.H0_basis):
            s = (s + c * bv.astype(object)) % self.Fp.p
        return s.astype(int)

    def is_global_section(self, s: np.ndarray) -> bool:
        check = self.Fp.mat_mod(
            self.delta0.astype(object) @ s.astype(object).reshape(-1, 1)
        )
        return bool(np.all(check == 0))

    def get_node_value(self, s: np.ndarray, v: int) -> np.ndarray:
        return s[v * self.dv : (v + 1) * self.dv].copy()

    def set_node_value(self, s: np.ndarray, v: int, val: np.ndarray):
        s[v * self.dv : (v + 1) * self.dv] = val

    def compose_restriction(self, path) -> np.ndarray:
        R = np.eye(self.dv, dtype=int)
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            e = (u, v) if (u, v) in self.graph.edges else (v, u)
            rho_v = self.rho[(e, v)]
            rho_v_inv = self.Fp.mat_inv(rho_v)
            R = self.Fp.mat_mul(rho_v_inv, R)
        return R

    def tree_propagation_maps(self, root: int = 0) -> dict:
        """Compute R_{root→v} for all nodes on a tree (or spanning tree)."""
        if self.graph.is_tree:
            tree_edges = self.graph.edges
        else:
            tree_edges, _ = self.graph.spanning_tree(root)

        adj = {v: [] for v in self.graph.nodes}
        for u, v in tree_edges:
            adj[u].append(v)
            adj[v].append(u)

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
                e = (min(u, v), max(u, v))
                rho_v = self.rho[(e, v)]
                rho_u = self.rho[(e, u)]
                rho_v_inv = self.Fp.mat_inv(rho_v)
                transfer = self.Fp.mat_mul(rho_v_inv, rho_u)
                R_to[v] = self.Fp.mat_mul(transfer, R_to[u])
        return R_to

    def __repr__(self):
        return (
            f"Sheaf(graph={self.graph}, dv={self.dv}, "
            f"F={self.Fp}, dim H^0={self.H0_dim})"
        )
