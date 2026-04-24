"""
Pseudorandom generator for nonce-derived randomness in the v3 protocol.

Instantiates PRG: {0,1}^λ × [n] → F_p^d via SHAKE-256, consuming bytes
from an extendable-output function keyed by (nonce || vertex_index).

This is the randomness source r_v = PRG(N, v) used inside the v3
encryption equation, which makes encryption probabilistic and defeats
the polynomial interpolation attack that broke v1/v2.
"""

from __future__ import annotations
import hashlib
import numpy as np
from math import ceil, log2


def prg_derive(nonce: bytes, vertex: int, d: int, p: int) -> np.ndarray:
    """
    Derive a d-dimensional vector in F_p^d from (nonce, vertex).

    Uses SHAKE-256 with input nonce || little-endian(vertex, 4 bytes) and
    rejection sampling per coordinate to avoid modulo bias.

    Args:
        nonce: arbitrary bytes (typically 16 or 32 bytes).
        vertex: integer index 0 <= vertex < 2^32.
        d: dimension.
        p: prime modulus.

    Returns:
        numpy int64 array of length d with entries in [0, p).
    """
    if vertex < 0 or vertex >= (1 << 32):
        raise ValueError(f"vertex out of range: {vertex}")

    bits = ceil(log2(p)) if p > 1 else 1
    bytes_per = max(1, (bits + 7) // 8 + 1)  # extra byte to reduce bias
    mask = (1 << (bytes_per * 8)) - 1

    shake = hashlib.shake_256()
    shake.update(b"alaniz-v3-prg")
    shake.update(nonce)
    shake.update(int(vertex).to_bytes(4, "little"))

    # Stream enough bytes upfront, extend on rejection
    out = shake.digest(d * bytes_per * 2)
    idx = 0

    result = np.zeros(d, dtype=np.int64)
    for i in range(d):
        while True:
            if idx + bytes_per > len(out):
                # Extend output if rejection ate through the buffer
                shake2 = hashlib.shake_256()
                shake2.update(b"alaniz-v3-prg-extend")
                shake2.update(nonce)
                shake2.update(int(vertex).to_bytes(4, "little"))
                shake2.update(int(i).to_bytes(4, "little"))
                out = out + shake2.digest(d * bytes_per * 4)
            chunk = out[idx:idx + bytes_per]
            idx += bytes_per
            val = int.from_bytes(chunk, "little") & mask
            # Rejection sampling: accept if val < floor(mask+1 / p) * p
            threshold = ((mask + 1) // p) * p
            if val < threshold:
                result[i] = val % p
                break
    return result


def prg_bytes(nonce: bytes, tag: bytes, length: int) -> bytes:
    """Generic SHAKE-256 output for the KEM construction."""
    shake = hashlib.shake_256()
    shake.update(tag)
    shake.update(nonce)
    return shake.digest(length)
