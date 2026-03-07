"""
Graph structures for sheaf-based cryptography.

Supports paths, trees, stars, and custom graphs. The protocol
requires tree-like graphs (n > m) for non-trivial H^0.
"""

from __future__ import annotations
import random as _random
from dataclasses import dataclass


@dataclass(frozen=True)
class Graph:
    """Undirected graph represented by node and edge lists."""

    nodes: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]

    @property
    def n(self) -> int:
        return len(self.nodes)

    @property
    def m(self) -> int:
        return len(self.edges)

    @property
    def is_tree(self) -> bool:
        return self.m == self.n - 1 and self.n > 0

    @property
    def expected_H0_dim(self) -> int:
        """For generic restriction maps, dim H^0 ≈ max(0, n - m) * d."""
        return max(0, self.n - self.m)

    # --- Constructors ---

    @staticmethod
    def path(n: int) -> Graph:
        """Path graph P_n: n nodes, n-1 edges."""
        return Graph(
            nodes=tuple(range(n)),
            edges=tuple((i, i + 1) for i in range(n - 1)),
        )

    @staticmethod
    def star(n: int) -> Graph:
        """Star graph S_n: central node 0, n-1 leaves."""
        return Graph(
            nodes=tuple(range(n)),
            edges=tuple((0, i) for i in range(1, n)),
        )

    @staticmethod
    def binary_tree(depth: int) -> Graph:
        """Complete binary tree of given depth."""
        n = 2 ** (depth + 1) - 1
        nodes = tuple(range(n))
        edges = []
        for i in range(n):
            left, right = 2 * i + 1, 2 * i + 2
            if left < n:
                edges.append((i, left))
            if right < n:
                edges.append((i, right))
        return Graph(nodes=nodes, edges=tuple(edges))

    @staticmethod
    def caterpillar(spine: int, pendants: int = 1) -> Graph:
        """Caterpillar graph: path spine with pendant vertices."""
        nodes_list = list(range(spine))
        edges_list = [(i, i + 1) for i in range(spine - 1)]
        idx = spine
        for i in range(spine):
            for _ in range(pendants):
                nodes_list.append(idx)
                edges_list.append((i, idx))
                idx += 1
        return Graph(nodes=tuple(nodes_list), edges=tuple(edges_list))

    @staticmethod
    def random_tree(n: int, rng=None) -> Graph:
        """Random tree on n nodes via random Prüfer sequence."""
        r = rng or _random
        if n <= 2:
            edges = [(0, 1)] if n == 2 else []
            return Graph(nodes=tuple(range(n)), edges=tuple(edges))
        # Prüfer → tree
        seq = [r.randint(0, n - 1) for _ in range(n - 2)]
        degree = [1] * n
        for v in seq:
            degree[v] += 1
        edges = []
        for v in seq:
            for u in range(n):
                if degree[u] == 1:
                    edges.append((min(u, v), max(u, v)))
                    degree[u] -= 1
                    degree[v] -= 1
                    break
        last = [i for i in range(n) if degree[i] == 1]
        if len(last) == 2:
            edges.append((min(last), max(last)))
        return Graph(nodes=tuple(range(n)), edges=tuple(sorted(edges)))

    def adjacency_list(self) -> dict[int, list[int]]:
        adj = {v: [] for v in self.nodes}
        for u, v in self.edges:
            adj[u].append(v)
            adj[v].append(u)
        return adj

    def __repr__(self):
        return f"Graph(|V|={self.n}, |E|={self.m})"
