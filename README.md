# The Alaniz Cipher

**Encryption from Nonlinear Sheaf Morphisms over Graphs**

A cryptographic primitive based on the Nonlinear Sheaf Morphism Inversion Problem (NL-SMIP): a polynomial system with secret coefficients coupled by sheaf cohomology over graphs.

> **Paper**: L. Alaniz Pintos, *"The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs"*, IACR ePrint, 2026.  
> <!-- Update the line below when you have the ePrint number -->
> <!-- https://eprint.iacr.org/2026/XXXX -->

## Overview

Messages are encoded as global sections of a cellular sheaf — elements of H⁰(G, F₀) — and encrypted via:

$$c_v = A_v \cdot s_v + B_v \cdot \sigma(A_v \cdot s_v)$$

where (Aᵥ, Bᵥ) ∈ GL(d, 𝔽ₚ)² are secret matrices per node and σ is a public nonlinear map (multiplicative inverse or cube map). Decryption uses tree propagation with cohomological filtering in O(n) time.

### Key Properties

| Property | Status | Evidence |
|---|---|---|
| Round-trip correctness | ✓ | 106/106 configurations, 0 failures |
| Decryption uniqueness | ✓ | Pr[ambiguity] ≤ (L−1) · γⁿ⁻¹ (Theorem 7.3) |
| CPA resistance | ✓ | Degree-3+ polynomial system, linearization fails |
| Post-quantum resistance | ✓ | No reduction to HSP; Grover is only quantum threat |

### Security (NL-SMIP)

The attacker faces a polynomial system whose **coefficients are secret** (unlike MQ/Rainbow where the system is public). NL-SMIP is proven at least NP-hard via reduction from MQ (Lemma 6.3 in the paper).

| Parameter Set | dᵥ | log₂ p | n | σ | Classical | Quantum |
|---|---|---|---|---|---|---|
| Demo | 2 | 5 | 8 | cube | ~2³⁸ | ~2¹⁹ |
| Standard | 4 | 31 | 16 | cube | ~2⁹³ | ~2⁴⁷ |
| PQ-128 | 8 | 61 | 32 | inverse | ~2²⁰⁰ | ~2¹⁰⁰ |
| PQ-256 | 8 | 127 | 64 | inverse | ~2⁴⁰⁰ | ~2²⁰⁰ |

## Quick Start

```bash
# Clone
git clone https://github.com/YOUR_USERNAME/alaniz-cipher.git
cd alaniz-cipher

# Run the demo
python -m alaniz.demo.demo_basic

# Run all tests (29 tests, ~5 seconds)
pip install pytest
pytest tests/ -v
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
│   │   ├── sigma.py        # Nonlinear maps σ (inverse, cube)
│   │   └── protocol.py     # KeyGen, Encrypt, Decrypt, Encode/Decode
│   └── demo/
│       └── demo_basic.py   # Interactive walkthrough
├── tests/
│   └── test_protocol.py    # Full test suite
├── pyproject.toml
└── README.md
```

## Test Suite

The 29 automated tests cover:

- **Round-trip correctness**: 10 random messages + zero section + basis sections
- **Multiple primes**: 𝔽₁₇, 𝔽₂₃, 𝔽₂₉, 𝔽₄₇
- **Both σ functions**: inverse (x^{p−2}) and cube (x³)
- **8 graph topologies**: paths, stars, binary trees, caterpillars, random trees
- **Decryption uniqueness**: brute-force confirms single valid section
- **CPA resistance**: linear system inconsistency verified
- **Structural properties**: H⁰ dimension, section validity, tree propagation

## Requirements

- Python ≥ 3.11
- NumPy ≥ 1.24
- pytest ≥ 7.0 (for tests)

## Citation

```bibtex
@misc{alaniz2026cipher,
  author = {Alaniz Pintos, Lucas},
  title = {The Alaniz Cipher: Encryption from Nonlinear Sheaf Morphisms over Graphs},
  year = {2026},
  howpublished = {IACR ePrint},
  note = {https://eprint.iacr.org/2026/XXXX}
}
```

## Author

**Dr. Lucas Alaniz Pintos**  
Gerencia de Smart Products, INECO  
lucas.alaniz@ineco.com  
ORCID: [0009-0008-5179-2534](https://orcid.org/0009-0008-5179-2534)

## License

MIT
