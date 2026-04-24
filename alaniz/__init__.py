"""Alaniz Cipher — Sheaf-based post-quantum encryption.

v3 (current): probabilistic scheme with vector-valued σ in F_{p^d}
              and polynomial-time decryption via univariate factorization.
              IND-CPA as PKE, IND-CCA after Fujisaki-Okamoto transform.

v2 (deprecated): broken by interpolation + Gröbner attack.
v1 (deprecated): broken by Langa's scaling attack.
"""
__version__ = "3.0.0"
