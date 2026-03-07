"""
Protocolo Alaniz — Test Suite

Tests:
  1. Round-trip correctness: Dec(Enc(m)) = m
  2. Decryption uniqueness: single valid global section
  3. CPA resistance: attacker's linear system is inconsistent
  4. Scaling: multiple primes, dimensions, topologies
  5. Byte encoding/decoding round-trip
"""

import sys
import os
import numpy as np
import random
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.sigma import Sigma
from alaniz.crypto.protocol import Protocol, PublicParams, KeyPair


# ============================================================
# FIXTURES
# ============================================================

def make_protocol(n_nodes=6, dv=2, p=17, sigma_name="inverse", seed=42):
    """Helper to build a protocol instance."""
    rng = random.Random(seed)
    Fp = FiniteField(p)
    G = Graph.path(n_nodes)
    sheaf = Sheaf.random(G, dv, Fp, rng=rng)
    params = PublicParams.generate(sheaf, sigma_name)
    proto = Protocol(params)
    key = proto.keygen(rng=rng)
    return proto, key, rng


# ============================================================
# TEST 1: ROUND-TRIP CORRECTNESS
# ============================================================

class TestRoundTrip:
    """Dec(Enc(s)) = s for all valid sections."""

    def test_basic_roundtrip(self):
        proto, key, rng = make_protocol()
        for _ in range(10):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            recovered = proto.decrypt_tree(c, key)
            assert recovered is not None, "Decryption returned None"
            assert np.array_equal(recovered, s), "Round-trip failed"

    def test_zero_section(self):
        proto, key, rng = make_protocol()
        s = proto.sheaf.section_from_coeffs([0] * proto.sheaf.H0_dim)
        c = proto.encrypt(s, key)
        recovered = proto.decrypt_tree(c, key)
        assert np.array_equal(recovered, s)

    def test_basis_sections(self):
        proto, key, rng = make_protocol()
        for i in range(proto.sheaf.H0_dim):
            coeffs = [0] * proto.sheaf.H0_dim
            coeffs[i] = 1
            s = proto.sheaf.section_from_coeffs(coeffs)
            c = proto.encrypt(s, key)
            recovered = proto.decrypt_tree(c, key)
            assert np.array_equal(recovered, s), f"Basis section {i} failed"

    @pytest.mark.parametrize("p", [17, 23, 29, 47])
    def test_various_primes(self, p):
        sigma = "inverse"
        proto, key, rng = make_protocol(p=p, sigma_name=sigma)
        for _ in range(5):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            recovered = proto.decrypt_tree(c, key)
            assert np.array_equal(recovered, s)

    @pytest.mark.parametrize("p", [17, 23, 29])
    def test_cube_sigma(self, p):
        Fp = FiniteField(p)
        if not Fp.cube_invertible():
            pytest.skip(f"Cube not invertible over F_{p}")
        proto, key, rng = make_protocol(p=p, sigma_name="cube")
        for _ in range(5):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            recovered = proto.decrypt_tree(c, key)
            assert np.array_equal(recovered, s)


# ============================================================
# TEST 2: TOPOLOGIES
# ============================================================

class TestTopologies:
    """Round-trip works on various tree topologies."""

    def _test_topology(self, graph, p=17, sigma="inverse", n_tests=5):
        rng = random.Random(123)
        Fp = FiniteField(p)
        sheaf = Sheaf.random(graph, 2, Fp, rng=rng)
        if sheaf.H0_dim == 0:
            pytest.skip("H^0 is trivial for this graph")
        params = PublicParams.generate(sheaf, sigma)
        proto = Protocol(params)
        key = proto.keygen(rng=rng)
        for _ in range(n_tests):
            s = sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            recovered = proto.decrypt_tree(c, key)
            assert np.array_equal(recovered, s)

    def test_path_8(self):
        self._test_topology(Graph.path(8))

    def test_path_12(self):
        self._test_topology(Graph.path(12))

    def test_path_20(self):
        self._test_topology(Graph.path(20))

    def test_star_9(self):
        self._test_topology(Graph.star(9))

    def test_binary_tree_depth2(self):
        self._test_topology(Graph.binary_tree(2))

    def test_binary_tree_depth3(self):
        self._test_topology(Graph.binary_tree(3))

    def test_caterpillar(self):
        self._test_topology(Graph.caterpillar(5, 2))

    def test_random_tree(self):
        rng = random.Random(999)
        G = Graph.random_tree(12, rng=rng)
        self._test_topology(G)


# ============================================================
# TEST 3: DECRYPTION UNIQUENESS
# ============================================================

class TestUniqueness:
    """Brute-force decryption returns exactly one valid section."""

    def test_unique_bruteforce(self):
        proto, key, rng = make_protocol(n_nodes=6, p=17)
        for _ in range(5):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            valid = proto.decrypt_bruteforce(c, key)
            assert len(valid) == 1, f"Expected 1 solution, got {len(valid)}"
            assert np.array_equal(valid[0], s)

    def test_tree_matches_bruteforce(self):
        """Tree propagation and brute force give the same answer."""
        proto, key, rng = make_protocol(n_nodes=6, p=17)
        for _ in range(5):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            tree_result = proto.decrypt_tree(c, key)
            bf_results = proto.decrypt_bruteforce(c, key)
            assert len(bf_results) == 1
            assert np.array_equal(tree_result, bf_results[0])


# ============================================================
# TEST 4: CPA RESISTANCE
# ============================================================

class TestCPAResistance:
    """
    Verify that the attacker's system is nonlinear and
    linear CPA recovery fails.
    """

    def test_linear_attack_fails(self):
        """
        CPA with known plaintexts: the linear system for recovering
        the key should be INCONSISTENT (because the true system
        has nonlinear terms from σ).
        """
        proto, key, rng = make_protocol(n_nodes=6, p=17, sigma_name="inverse")
        Fp = proto.Fp
        dv = proto.dv
        v = 0  # attack node 0

        # Collect 8 CPA queries (enough for linear system if it were linear)
        queries = []
        for i in range(8):
            s = proto.sheaf.random_section(rng=rng)
            c = proto.encrypt(s, key)
            sv = proto.sheaf.get_node_value(s, v)
            cv = proto.sheaf.get_node_value(c, v)
            queries.append((sv, cv))

        # Build linear system: assuming c_v = M_v · s_v
        # (which would be the case without σ)
        # 2 equations per query, 4 unknowns (entries of M_v)
        n_eq = len(queries) * dv
        n_unk = dv * dv
        aug = np.zeros((n_eq, n_unk + 1), dtype=object)

        for j, (sv, cv) in enumerate(queries):
            for i in range(dv):
                for k in range(dv):
                    aug[j * dv + i, i * dv + k] = int(sv[k]) % Fp.p
                aug[j * dv + i, n_unk] = int(cv[i]) % Fp.p

        # Row reduce
        p = Fp.p
        pivots, row = [], 0
        A = aug.copy()
        for col in range(n_unk):
            found = False
            for r in range(row, n_eq):
                if A[r, col] % p != 0:
                    found = True
                    A[[row, r]] = A[[r, row]]
                    break
            if not found:
                continue
            pivots.append(col)
            iv = pow(int(A[row, col]), p - 2, p)
            A[row] = (A[row] * iv) % p
            for r in range(n_eq):
                if r != row and A[r, col] % p != 0:
                    A[r] = (A[r] - A[r, col] * A[row]) % p
            row += 1

        # Check consistency: if linear model were correct, system would
        # be consistent. With σ, it should be INCONSISTENT.
        inconsistent = False
        for r in range(len(pivots), n_eq):
            if A[r, n_unk] % p != 0:
                inconsistent = True
                break

        assert inconsistent, (
            "Linear CPA system is consistent — σ may not be providing "
            "sufficient nonlinearity!"
        )


# ============================================================
# TEST 5: BYTE ENCODING
# ============================================================

class TestEncoding:
    """Encode and decode byte messages."""

    def test_short_message(self):
        proto, key, rng = make_protocol(p=17)
        # p=17, H0_dim=2 → capacity ~10 bits → 1 byte max
        msg = b"\x05"
        s = proto.encode(msg)
        assert proto.sheaf.is_global_section(s)
        c = proto.encrypt(s, key)
        s_dec = proto.decrypt_tree(c, key)
        recovered = proto.decode(s_dec)
        assert recovered == msg

    def test_single_byte(self):
        proto, key, rng = make_protocol(p=251)
        msg = b"A"
        s = proto.encode(msg)
        c = proto.encrypt(s, key)
        s_dec = proto.decrypt_tree(c, key)
        recovered = proto.decode(s_dec)
        assert recovered == msg


# ============================================================
# TEST 6: STRUCTURAL PROPERTIES
# ============================================================

class TestStructure:
    """Verify structural properties of the mathematical objects."""

    def test_H0_dimension_path(self):
        """Path graph P_n with d_v = d gives dim H^0 = d."""
        for n in [4, 6, 8, 12]:
            Fp = FiniteField(17)
            G = Graph.path(n)
            sheaf = Sheaf.random(G, 2, Fp, rng=random.Random(42))
            assert sheaf.H0_dim == 2, f"P_{n}: expected dim H^0 = 2, got {sheaf.H0_dim}"

    def test_H0_dimension_star(self):
        Fp = FiniteField(17)
        G = Graph.star(9)
        sheaf = Sheaf.random(G, 2, Fp, rng=random.Random(42))
        assert sheaf.H0_dim == 2

    def test_sections_are_global(self):
        Fp = FiniteField(17)
        G = Graph.path(8)
        sheaf = Sheaf.random(G, 2, Fp, rng=random.Random(42))
        rng = random.Random(100)
        for _ in range(20):
            s = sheaf.random_section(rng=rng)
            assert sheaf.is_global_section(s)

    def test_tree_propagation_consistent(self):
        """Sections built by tree propagation are valid global sections."""
        Fp = FiniteField(17)
        G = Graph.path(8)
        sheaf = Sheaf.random(G, 2, Fp, rng=random.Random(42))
        R_to = sheaf.tree_propagation_maps(root=0)

        rng = random.Random(77)
        for _ in range(10):
            x_root = Fp.random_vec(2, rng=rng)
            s = np.zeros(sheaf.C0_dim, dtype=int)
            for v in G.nodes:
                sheaf.set_node_value(s, v, Fp.mat_vec(R_to[v], x_root))
            assert sheaf.is_global_section(s)

    def test_plaintext_must_be_section(self):
        """Encrypting a non-section raises an error."""
        proto, key, rng = make_protocol()
        bad_s = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        with pytest.raises(ValueError, match="global section"):
            proto.encrypt(bad_s, key)

    def test_cube_requires_compatible_prime(self):
        """Cube sigma rejects primes where 3 | p-1."""
        Fp = FiniteField(31)  # 31 - 1 = 30, 3 | 30
        with pytest.raises(ValueError, match="not invertible"):
            Sigma("cube", Fp)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
