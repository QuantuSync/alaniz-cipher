# The Alaniz Cipher (v2)

**Encryption from Nonlinear Sheaf Morphisms over Graphs**

A cryptographic primitive based on the Nonlinear Sheaf Morphism Inversion Problem (NL-SMIP): a polynomial system with secret coefficients coupled by sheaf cohomology over graphs.

> **Paper**: L. Alaniz Pintos, *"The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs"*, 2026.  
> **Preprint (DOI)**: [10.5281/zenodo.19020097](https://doi.org/10.5281/zenodo.19020097)  
> **ResearchGate**: [Publication 401930354](https://www.researchgate.net/publication/401930354_The_Alaniz_Cipher_Encryption_from_Nonlinear_Sheaf_Morphisms_over_Graphs_SPN_Fix)

## What changed in v2

Version 1 used component-wise σ functions (cube, multiplicative inverse) that were vulnerable to an O(d)-query CPA key-recovery attack discovered by A. Rodríguez Langa (INECO, 2026). A deeper analysis revealed that **any σ without a linear part** is broken by polynomial interpolation.

Version 2 introduces σ_SPN(y) = y + S(M·S(y) + c), a two-round SPN with three simultaneous defenses:

| Component | Defeats | Mechanism |
|---|---|---|
| Linear part (y) | Interpolation attack | coeff[t¹] = (I+B)A, not A |
| Mixing matrix (M) | Differential column attack | σ(y)ᵢ depends on all components |
| Constant (c = (1,...,1)) | Scaling attack (Langa) | (λ³w+1)³ ≠ λᵉ(w+1)³ |

**25 attack vectors tested. 0 succeeded.** See Experiment 14 for the full red team.

## Overview

Messages are encoded as global sections of a cellular sheaf — elements of H⁰(G, F₀) — and encrypted via:

$$c_v = A_v \cdot s_v + B_v \cdot \sigma(A_v \cdot s_v)$$

where (Aᵥ, Bᵥ) ∈ GL(d, 𝔽ₚ)² are secret matrices per node and σ = σ_SPN is a public nonlinear map. Decryption uses tree propagation with cohomological filtering in O(n) time.

### Security (NL-SMIP)

The attacker's best strategy is polynomial interpolation to obtain (I+B)A, then Bézout-bounded decomposition to recover A and B separately. The binding security estimate is 9^(d²):

| Parameter Set | d | log₂ p | n | Classical | Post-Grover | PQ? |
|---|---|---|---|---|---|---|
| Demo | 2 | 5 | 8 | ~2¹³ | ~2⁷ | No |
| Standard | 4 | 31 | 16 | ~2⁵¹ | ~2²⁶ | No |
| PQ-128 | 8 | 61 | 32 | ~2²⁰³ | ~2¹⁰¹ | Yes |
| PQ-256 | 12 | 127 | 64 | ~2⁴⁵⁷ | ~2²²⁸ | Yes |

## Quick Start

```bash
git clone https://github.com/QuantuSync/alaniz-cipher.git
cd alaniz-cipher

# Run the demo (uses σ_SPN by default)
python -m alaniz.demo.demo_basic

# Run all tests (29 tests)
pip install pytest numpy
pytest tests/ -v

# Run v2 red team (25 attacks, 65 round-trips)
python experiments/14_v2_redteam.py

# Run all experiments
for f in experiments/0*.py experiments/14*.py; do python "$f"; done
```

## Project Structure

```
alaniz-cipher/
├── alaniz/
│   ├── core/
│   │   ├── field.py        # 𝔽ₚ arithmetic, matrix operations
│   │   ├── graph.py        # Graph topologies (paths, trees, stars)
│   │   └── sheaf.py        # Cellular sheaves, H⁰, coboundary, tree propagation
│   ├── crypto/
│   │   ├── sigma.py        # Nonlinear maps: id_spn (v2), cube, inverse (v1, deprecated)
│   │   └── protocol.py     # KeyGen, Encrypt, Decrypt, Encode/Decode
│   └── demo/
│       └── demo_basic.py   # Interactive walkthrough
├── experiments/
│   ├── 01_scaling_verification.py   # 176/176 round-trips (v1 sigmas)
│   ├── 02_cpa_attack_analysis.py    # 4 broken designs (v1)
│   ├── 03_uniqueness_analysis.py    # Preimage distribution, γ bound
│   ├── 04_sigma_selection.py        # σ candidate evaluation
│   ├── 05_complexity_analysis.py    # XL attack costs (v1 estimates)
│   ├── 06_postquantum_analysis.py   # Quantum attack surface
│   └── 14_v2_redteam.py            # v2: 25 attacks, 65 round-trips, σ_SPN
├── tests/
│   └── test_protocol.py    # Full test suite (29 tests, σ_SPN default)
├── pyproject.toml
└── README.md
```

**Note:** Experiments 01–06 use v1 sigmas (cube, inverse) for historical reproducibility of the paper's design evolution section. Experiment 14 is the authoritative v2 verification.

## Citation

```bibtex
@misc{alaniz2026cipher,
  author       = {Alaniz Pintos, Lucas},
  title        = {The Alaniz Cipher: Encryption from Nonlinear
                  Sheaf Morphisms over Graphs},
  year         = {2026},
  doi          = {10.5281/zenodo.19020097},
  url          = {https://doi.org/10.5281/zenodo.19020097},
  publisher    = {Zenodo}
}
```

## Acknowledgments

A. Rodríguez Langa (INECO) for the external security analysis that identified the scaling homogeneity vulnerability in v1, directly motivating the σ_SPN construction.

## Author

**Dr. Lucas Alaniz Pintos**  
Smart Products Division, INECO  
lucas.alaniz@ineco.com  
ORCID: [0009-0008-5179-2534](https://orcid.org/0009-0008-5179-2534)

## License

MIT
