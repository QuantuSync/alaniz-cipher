# The Alaniz Cipher (v3)

**Sheaf-Based Post-Quantum Encryption via Vector-Valued Finite Field Permutations**

A cryptographic primitive based on the Nonlinear Sheaf Morphism Inversion Problem (NL-SMIP): a polynomial system with secret coefficients coupled by sheaf cohomology over graphs.

> **Paper (v3)**: L. Alaniz Pintos, *"The Alaniz Cipher v3: Sheaf-Based Post-Quantum Encryption via Vector-Valued Finite Field Permutations"*, 2026.

## What changed in v3

Version 2 introduced σ_SPN to defeat the scaling attack that broke v1. However, v2 is itself broken by a more general attack combining polynomial interpolation with Gröbner basis computation: any deterministic encryption of the form `c = As + Bσ(As)` with σ of bounded polynomial degree admits a CPA key-recovery attack. See paper Section 3 and `experiments/02_cpa_attack_analysis.py`.

Version 3 addresses this at the structural level:

| Feature | v2 | v3 |
|---|---|---|
| Encryption | Deterministic | Probabilistic (nonce + PRG) |
| σ | Component-wise SPN over F_p | Vector-valued: π_e(x) = x^e in F_{p^d} |
| Secondary key | Matrix B_v ∈ GL(d, F_p) | Scalar β_v ∈ F_{p^d}^* |
| Decryption | Brute force over F_p^d | Univariate factorization in F_{p^d}[τ] |
| Differential uniformity of σ | δ ≥ 2p (~34 for p=17) | δ ∈ {2, 4} (APN or nearly-APN) |
| IND-CCA target | Not applicable | Achieved via FO transform |

**14 semi-regularity configurations tested, 14 match Hilbert-Poincaré predictions exactly.** No MinRank vulnerability by Theorem 6.5 (centralizer argument). No weak keys detected in 69 configurations.

## Overview

Messages are encoded as global sections of a cellular sheaf on a graph with cycle rank β_1 ≥ 1. Per-node encryption:

```
r_v = PRG(nonce, v)
u_v = ι(A_v · s_v + r_v)           ∈ F_{p^d}
w_v = β_v · u_v + (β_v − 1) · (L · u_v + 1)^e
c_v = ι^{-1}(w_v) − r_v
```

where `A_v ∈ GL(d, F_p)`, `β_v ∈ F_{p^d}^*`, and `L ∈ F_{p^d}^*` is public. Decryption reduces to finding roots of `γ·τ^e + α·τ − (c′ + α) = 0` in `F_{p^d}[τ]` via Cantor-Zassenhaus, then filtering via the cohomological consistency at all `n − 1` non-root vertices.

### Parameter Sets

| Set | d | log₂ p | n | e | |sk| | |pk| | |c| | Security |
|---|---|---|---|---|---|---|---|---|
| Demo | 2 | 4 | 4 | 17 | 30 B | 10 B | 16 B | not secure |
| Academic | 3 | 16 | 8 | 17 | 0.3 KB | 30 B | 96 B | ~80 bits classical |
| PQ-128 | 6 | 61 | 32 | 17 | 16 KB | 64 B | 1.5 KB | ≥128 post-Grover |
| PQ-256 | 8 | 127 | 64 | 17 | 170 KB | 128 B | 8 KB | ≥256 post-Grover |

## Quick Start

```bash
git clone https://github.com/QuantuSync/alaniz-cipher.git
cd alaniz-cipher
pip install -e .[test]

# Run the v2 demo (historical)
python -m alaniz.demo.demo_basic

# Run all tests
pytest tests/ -v

# Run v3 verification experiments
python experiments/07_semi_regularity.py
python experiments/08_weak_keys.py
python experiments/09_nonce_robustness.py
python experiments/10_centralizer.py
python experiments/11_kem_verification.py
```

### Minimal v3 example

```python
from alaniz.core.field import FiniteField
from alaniz.core.graph import Graph
from alaniz.core.sheaf import Sheaf
from alaniz.crypto.protocol import ProtocolV3, PublicParamsV3
from alaniz.crypto.kem import KEM

Fp = FiniteField(17)
graph = Graph.cycle(4)
sheaf = Sheaf.random_with_cohomology(graph, dv=2, Fp=Fp)
proto = ProtocolV3(PublicParamsV3.generate(sheaf))
kem   = KEM(proto)

sk       = kem.keygen()
K, ct    = kem.encaps(sk)     # sender: random message + shared key
K_recv   = kem.decaps(sk, ct) # receiver: same K if valid
assert K == K_recv
```

## Project Structure

```
alaniz-cipher/
├── alaniz/
│   ├── core/
│   │   ├── field.py        # F_p arithmetic (v2, still used by v1/v2)
│   │   ├── field_ext.py    # F_{p^d} arithmetic via galois (v3)
│   │   ├── graph.py        # Graph topologies (trees, cycles, etc.)
│   │   └── sheaf.py        # Cellular sheaves + cohomology constraints
│   ├── crypto/
│   │   ├── sigma.py        # σ maps: monomial_power (v3), id_spn (v2), cube/inverse (v1)
│   │   ├── protocol.py     # Protocol (v2), ProtocolV3 (v3)
│   │   ├── prg.py          # SHAKE-256 PRG for nonce derivation (v3)
│   │   └── kem.py          # Fujisaki-Okamoto → IND-CCA KEM (v3)
│   └── demo/
│       └── demo_basic.py   # v2 demo (historical)
├── experiments/
│   ├── 01-06_*.py          # v1 analyses (historical, paper v2 Appendix)
│   ├── 14_v2_redteam.py    # v2 red team (historical)
│   ├── 07_semi_regularity.py    # v3: Hilbert-Poincaré verification
│   ├── 08_weak_keys.py          # v3: 69-config weak-key sweep
│   ├── 09_nonce_robustness.py   # v3: Theorem 6.7
│   ├── 10_centralizer.py        # v3: Theorem 6.5
│   └── 11_kem_verification.py   # v3: KEM IND-CCA empirical tests
├── tests/
│   ├── test_protocol.py      # v2 tests (kept for regression)
│   └── test_protocol_v3.py   # v3 tests
├── pyproject.toml
└── README.md
```

## v1 and v2: Historical Status

**v1**: Broken by Langa's scaling attack on `σ_cube(λy) = λ³ σ_cube(y)` (Alaniz 2026, Appendix A of paper v3). Recovery in O(d) queries.

**v2**: Broken by polynomial interpolation + Gröbner attack against any deterministic `c = As + Bσ(As)` with bounded-degree σ. For d=2, p=17, full key recovery in 20 CPA queries and ~5.7s of SymPy computation across 10/10 random instances. See `experiments/02_cpa_attack_analysis.py` and paper Section 3.

Both v1 and v2 classes remain accessible in this repository for historical reproducibility. For any new use, use `ProtocolV3`.

## Citation

```bibtex
@misc{alaniz2026cipherv3,
  author  = {Alaniz Pintos, Lucas},
  title   = {The Alaniz Cipher v3: Sheaf-Based Post-Quantum Encryption
             via Vector-Valued Finite Field Permutations},
  year    = {2026},
  url     = {https://github.com/QuantuSync/alaniz-cipher}
}
```

## Acknowledgments

A. Rodríguez Langa (INECO) for the cryptanalysis of v1 that motivated the redesign path leading to v2 and ultimately to v3.

## Author

**Lucas Alaniz Pintos**
Smart Products Division, INECO
lucas.alaniz@ineco.com
ORCID: [0009-0008-5179-2534](https://orcid.org/0009-0008-5179-2534)

## License

MIT
