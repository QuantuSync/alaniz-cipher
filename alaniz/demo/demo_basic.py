#!/usr/bin/env python3
"""
Protocolo Alaniz v1.0 — Interactive Demo

Demonstrates the full encryption/decryption pipeline:
  1. Setup: graph, sheaf, public parameters
  2. Key generation
  3. Message encoding → global section
  4. Encryption
  5. Decryption via tree propagation
  6. Message recovery
"""

import sys
import os
import random
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import Protocol, PublicParams


def main():
    print("=" * 60)
    print("  PROTOCOLO ALANIZ v1.0")
    print("  Sheaf Homology Encryption Cipher (SHEC)")
    print("=" * 60)

    # --- Parameters ---
    p = 17
    dv = 2
    n_nodes = 8
    sigma_name = "inverse"
    seed = 2024

    rng = random.Random(seed)
    Fp = FiniteField(p)

    print(f"\n  Field:     F_{p}")
    print(f"  Fiber dim: d_v = {dv}")
    print(f"  Sigma:     {sigma_name}")

    # --- Graph ---
    G = Graph.path(n_nodes)
    print(f"\n  Graph:     P_{n_nodes} (path, {G.n} nodes, {G.m} edges)")

    # --- Sheaf ---
    sheaf = Sheaf.random(G, dv, Fp, rng=rng)
    print(f"  Sheaf:     dim H^0 = {sheaf.H0_dim}")
    print(f"  Capacity:  {sheaf.H0_dim} coefficients in F_{p}")

    # --- Protocol ---
    params = PublicParams.generate(sheaf, sigma_name)
    proto = Protocol(params)

    # --- Key Generation ---
    print("\n--- KEY GENERATION ---")
    t0 = time.time()
    key = proto.keygen(rng=rng)
    t_keygen = time.time() - t0
    print(f"  Generated {G.n} pairs (A_v, B_v) in GL({dv}, F_{p})")
    print(f"  Time: {t_keygen*1000:.1f} ms")

    # --- Message ---
    message = b"\xAB"  # 1 byte for p=17, H0_dim=2
    print(f"\n--- MESSAGE ---")
    print(f"  Original: 0x{message.hex()}")
    print(f"  Bytes:    {message.hex()}")

    # --- Encode ---
    print(f"\n--- ENCODE ---")
    s = proto.encode(message)
    print(f"  Section s ∈ H^0(G, F_0):")
    for v in G.nodes:
        sv = sheaf.get_node_value(s, v)
        print(f"    s_{v} = {sv}")
    print(f"  Global section check: δ^0(s) = 0 ✓" if sheaf.is_global_section(s)
          else "  ERROR: not a global section!")

    # --- Encrypt ---
    print(f"\n--- ENCRYPT ---")
    print(f"  c_v = A_v · s_v + B_v · σ(A_v · s_v)")
    t0 = time.time()
    c = proto.encrypt(s, key)
    t_enc = time.time() - t0
    for v in G.nodes:
        cv = sheaf.get_node_value(c, v)
        print(f"    c_{v} = {cv}")
    print(f"  Time: {t_enc*1000:.1f} ms")

    # --- Decrypt (tree propagation) ---
    print(f"\n--- DECRYPT (tree propagation) ---")
    t0 = time.time()
    s_recovered = proto.decrypt_tree(c, key)
    t_dec = time.time() - t0
    print(f"  Recovered section:")
    for v in G.nodes:
        sv = sheaf.get_node_value(s_recovered, v)
        print(f"    s_{v} = {sv}")
    print(f"  Match: {'✓' if s_recovered is not None and all(s_recovered == s) else '✗'}")
    print(f"  Time: {t_dec*1000:.1f} ms")

    # --- Decode ---
    print(f"\n--- DECODE ---")
    recovered_msg = proto.decode(s_recovered)
    print(f'  Recovered: 0x{recovered_msg.hex()}')
    print(f"  Match: {'✓' if recovered_msg == message else '✗'}")

    # --- Summary ---
    print(f"\n{'=' * 60}")
    print(f"  SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Round-trip: {'CORRECT ✓' if recovered_msg == message else 'FAILED ✗'}")
    print(f"  KeyGen:     {t_keygen*1000:.1f} ms")
    print(f"  Encrypt:    {t_enc*1000:.1f} ms")
    print(f"  Decrypt:    {t_dec*1000:.1f} ms")
    print(f"  Plaintext:  {len(message)} bytes")
    print(f"  Ciphertext: {len(c) * ((p-1).bit_length() + 7) // 8} bytes")
    print(f"  Expansion:  {len(c) / sheaf.H0_dim:.1f}x (field elements)")
    print()


if __name__ == "__main__":
    main()
