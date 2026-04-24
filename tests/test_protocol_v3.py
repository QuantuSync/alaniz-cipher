"""Tests for ProtocolV3 and the FO-KEM."""
from __future__ import annotations
import os
import numpy as np
import pytest
import random as _random

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import ProtocolV3, PublicParamsV3


# ----------------------------------------------------------------- Fixtures

@pytest.fixture
def small_proto():
    """d=2, p=17, 4-node cycle."""
    Fp = FiniteField(17)
    graph = Graph.cycle(4)
    rng = _random.Random(42)
    sheaf = Sheaf.random_with_cohomology(graph, 2, Fp, rng=rng)
    params = PublicParamsV3.generate(sheaf)
    proto = ProtocolV3(params)
    return proto, rng


@pytest.fixture
def small_proto_d3():
    """d=3, p=17, 4-node cycle."""
    Fp = FiniteField(17)
    graph = Graph.cycle(4)
    rng = _random.Random(42)
    sheaf = Sheaf.random_with_cohomology(graph, 3, Fp, rng=rng)
    params = PublicParamsV3.generate(sheaf)
    proto = ProtocolV3(params)
    return proto, rng


# -------------------------------------------------------------- KeyGen tests

def test_keygen_produces_required_fields(small_proto):
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    assert set(key.A.keys()) == set(proto.sheaf.graph.nodes)
    assert set(key.beta.keys()) == set(proto.sheaf.graph.nodes)
    assert len(key.s_reject) == 32


def test_keygen_beta_never_one_or_zero(small_proto):
    proto, rng = small_proto
    GF = proto.field_ext.GF
    for seed in range(20):
        key = proto.keygen(rng=_random.Random(seed))
        for v in proto.sheaf.graph.nodes:
            assert key.beta[v] != GF(0)
            assert key.beta[v] != GF(1)


def test_keygen_A_invertible(small_proto):
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    for v in proto.sheaf.graph.nodes:
        # Must not raise
        _ = proto.Fp.mat_inv(key.A[v])


# -------------------------------------------------------- Round-trip tests

def test_roundtrip_basic(small_proto):
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    for _ in range(10):
        s = proto.sheaf.random_section(rng=rng)
        nonce, c = proto.encrypt(s, key)
        s_rec = proto.decrypt(nonce, c, key)
        assert s_rec is not None, "Decryption returned None"
        assert np.array_equal(s_rec, s), "Decrypted section differs from plaintext"


def test_roundtrip_d3(small_proto_d3):
    proto, rng = small_proto_d3
    key = proto.keygen(rng=rng)
    for _ in range(5):
        s = proto.sheaf.random_section(rng=rng)
        nonce, c = proto.encrypt(s, key)
        s_rec = proto.decrypt(nonce, c, key)
        assert s_rec is not None
        assert np.array_equal(s_rec, s)


def test_encrypt_requires_global_section(small_proto):
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    # Random vector not in H^0 (likely)
    bogus = np.random.randint(0, proto.Fp.p, size=proto.sheaf.C0_dim)
    if proto.sheaf.is_global_section(bogus):
        # Unlikely but possible; just skip
        pytest.skip("Random vector happened to be a global section")
    with pytest.raises(ValueError):
        proto.encrypt(bogus, key)


def test_encrypt_deterministic_given_nonce(small_proto):
    """Same (s, nonce) → same ciphertext (required for the FO transform)."""
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    s = proto.sheaf.random_section(rng=rng)
    nonce = b"\x42" * 16
    _, c1 = proto.encrypt(s, key, nonce=nonce)
    _, c2 = proto.encrypt(s, key, nonce=nonce)
    assert np.array_equal(c1, c2)


def test_different_nonce_different_ciphertext(small_proto):
    proto, rng = small_proto
    key = proto.keygen(rng=rng)
    s = proto.sheaf.random_section(rng=rng)
    _, c1 = proto.encrypt(s, key, nonce=b"\x00" * 16)
    _, c2 = proto.encrypt(s, key, nonce=b"\xff" * 16)
    assert not np.array_equal(c1, c2)


# ------------------------------------------------------------------ KEM tests

def test_kem_roundtrip(small_proto):
    from alaniz.crypto.kem import KEM
    proto, rng = small_proto
    kem = KEM(proto)
    sk = kem.keygen(rng=rng)
    for _ in range(5):
        K_send, ct = kem.encaps(sk)
        K_recv = kem.decaps(sk, ct)
        assert K_send == K_recv


def test_kem_implicit_rejection(small_proto):
    from alaniz.crypto.kem import KEM
    proto, rng = small_proto
    kem = KEM(proto)
    sk = kem.keygen(rng=rng)
    K_orig, (nonce, c) = kem.encaps(sk)
    # Corrupt c
    c_bad = c.copy()
    c_bad[0] = (int(c_bad[0]) + 1) % proto.Fp.p
    K_bad = kem.decaps(sk, (nonce, c_bad))
    assert K_bad != K_orig
    # Deterministic rejection
    K_bad2 = kem.decaps(sk, (nonce, c_bad))
    assert K_bad == K_bad2
