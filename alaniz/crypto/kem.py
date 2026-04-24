"""
Fujisaki-Okamoto transform on Alaniz v3 → IND-CCA KEM.

Implements FO^⊥ with implicit rejection (Hofheinz-Hövelmanns-Kiltz, TCC 2017).

KEM.Encaps(pk):
  m ← Unif(F_p^d)
  (K, N) = SHAKE-256(m || pk)
  c = V3.Encrypt(pk, m; N)   # deterministic with N
  return (K, c)

KEM.Decaps(sk, (nonce, c)):
  m' = V3.Decrypt(sk, (nonce, c))
  (K', N') = SHAKE-256(m' || pk)
  if m' != ⊥ and re-encrypt matches (nonce, c):
      return K'
  else:
      return SHAKE-256(s_reject || nonce || c)   # implicit rejection
"""

from __future__ import annotations
import os
import hashlib
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple

from .protocol import ProtocolV3, KeyPairV3, PublicParamsV3


def _hash_split(m_bytes: bytes, pk_bytes: bytes,
                K_len: int = 32, nonce_len: int = 16) -> Tuple[bytes, bytes]:
    """(K, N) ← SHAKE-256(tag || m || pk)."""
    shake = hashlib.shake_256()
    shake.update(b"alaniz-v3-kem-H")
    shake.update(m_bytes)
    shake.update(pk_bytes)
    out = shake.digest(K_len + nonce_len)
    return out[:K_len], out[K_len:K_len + nonce_len]


def _hash_reject(s_reject: bytes, nonce: bytes, c_bytes: bytes,
                 K_len: int = 32) -> bytes:
    """Implicit rejection key: SHAKE-256(tag || s_reject || nonce || c)."""
    shake = hashlib.shake_256()
    shake.update(b"alaniz-v3-kem-Hrej")
    shake.update(s_reject)
    shake.update(nonce)
    shake.update(c_bytes)
    return shake.digest(K_len)


def _section_to_bytes(s: np.ndarray, p: int) -> bytes:
    """Canonical byte encoding of a section (or message vector)."""
    bits = (p - 1).bit_length()
    bytes_per = max(1, (bits + 7) // 8)
    out = bytearray()
    for val in s:
        out.extend(int(val).to_bytes(bytes_per, "big"))
    return bytes(out)


def _ciphertext_to_bytes(nonce: bytes, c: np.ndarray, p: int) -> bytes:
    return nonce + _section_to_bytes(c, p)


def _public_params_bytes(proto: ProtocolV3) -> bytes:
    """Canonical encoding of the public params."""
    return (
        f"alaniz-v3|p={proto.Fp.p}|d={proto.dv}|n={proto.sheaf.graph.n}|"
        f"L={int(proto.L_gf)}|e={proto.exponent}"
    ).encode()


class KEM:
    """
    IND-CCA KEM via Fujisaki-Okamoto transform on ProtocolV3.

    Typical use:
        params = PublicParamsV3.generate(sheaf)
        proto  = ProtocolV3(params)
        kem    = KEM(proto)
        sk     = kem.keygen()
        K, ct  = kem.encaps(sk)   # sender
        K'     = kem.decaps(sk, ct)   # receiver (same key if valid)
    """

    def __init__(self, proto: ProtocolV3):
        self.proto = proto
        self._pk_bytes = _public_params_bytes(proto)

    def keygen(self, rng=None) -> KeyPairV3:
        return self.proto.keygen(rng=rng)

    def encaps(self, sk: KeyPairV3, rng=None) -> Tuple[bytes, Tuple[bytes, np.ndarray]]:
        """Encapsulation: returns (shared_key, ciphertext)."""
        # Sample random message m in H^0
        s = self.proto.sheaf.random_section(rng=rng)
        m_bytes = _section_to_bytes(s, self.proto.Fp.p)
        # Derive (K, N) = H(m, pk)
        K, nonce = _hash_split(m_bytes, self._pk_bytes)
        # Deterministic encryption with N
        _, c = self.proto.encrypt(s, sk, nonce=nonce)
        return K, (nonce, c)

    def decaps(self, sk: KeyPairV3, ct: Tuple[bytes, np.ndarray]) -> bytes:
        """Decapsulation: returns shared_key (never ⊥; uses implicit rejection)."""
        nonce, c = ct
        c_bytes = _ciphertext_to_bytes(nonce, c, self.proto.Fp.p)
        # Decrypt
        s_rec = self.proto.decrypt(nonce, c, sk)
        if s_rec is None:
            return _hash_reject(sk.s_reject, nonce, c_bytes)
        # Re-encrypt and compare
        m_bytes = _section_to_bytes(s_rec, self.proto.Fp.p)
        K_rec, nonce_rec = _hash_split(m_bytes, self._pk_bytes)
        _, c_rec = self.proto.encrypt(s_rec, sk, nonce=nonce_rec)
        if nonce_rec == nonce and np.array_equal(c_rec, c):
            return K_rec
        else:
            return _hash_reject(sk.s_reject, nonce, c_bytes)
