# specialist-maths — Mathematics Specialist

**Author:** Michael Fuhrmann
**Build:** 199 | **Last updated:** 2026-03-13 | **Tests:** 10,403+

> **[Deutsche Version → README.de.md](README.de.md)**

---

## Disclaimer

> **This project is exclusively an exploration of the current capabilities of
> [Claude Code](https://claude.ai/claude-code) (Anthropic's AI coding assistant).**
>
> All mathematical papers, proofs, source code, Python modules, and test suites
> in this repository were **generated entirely by Claude Code** (models:
> claude-opus-4-6, claude-sonnet-4-6) under the direction of Michael Fuhrmann.
>
> The content is intended to demonstrate what an AI coding assistant can produce
> when applied to a complex, long-running scientific project. **No human wrote
> any of the mathematical papers or source code directly.**
>
> Mathematical claims should be treated with appropriate skepticism — while
> significant effort was invested in auditing correctness, AI-generated
> mathematical content may contain errors.

---

## Overview

This project systematically explores mathematical expertise across all areas —
from fundamentals to advanced research topics.

**Ultimate goal:** Investigation of the Millennium Prize Problems (RH, Goldbach,
BSD, P≠NP).

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

Free to use, modify, and distribute, provided the original author
(**Michael Fuhrmann**) is credited.

---

## Requirements

- Python 3.13+
- sympy >= 1.14, numpy >= 2.2, scipy >= 1.17, matplotlib >= 3.10
- mpmath, numba, networkx, gmpy2, pysat, z3-solver

## Installation

```bash
pip install sympy numpy scipy matplotlib mpmath numba networkx gmpy2 python-sat z3-solver pytest pytest-xdist --break-system-packages
```

## Running Tests

```bash
# Fast (parallel):
python3 -m pytest tests/ -n auto

# Verbose:
python3 -m pytest tests/ -v --timeout=60
```

---

## Project Structure

```
specialist-maths/
├── src/                               # 120+ Python modules
├── tests/                             # 10,403+ unit tests (TDD, pytest-xdist)
├── papers/                            # 97 LaTeX papers (Batches 1–25, EN + DE)
│   └── reviewed/                      # Batches 1–9 fully audited
├── reviews/                           # Audit reports
├── research/                          # Mathematical research & notes
├── dev-log/                           # Development diary
├── debugging/                         # Debug scripts
├── build/                             # Executables
├── claude-generated/                  # Other generated files
├── FEATURES.md                        # Feature list
├── BUGS.md                            # Bug tracker
├── OPTIMIZE.md                        # Optimization ideas
└── LICENSE                            # MIT License
```

---

## Mathematical Topics

### Foundations
1. **Algebra** — Polynomials, equations, number theory
2. **Analysis** — Differentiation, integration, Taylor series
3. **Linear Algebra** — Matrices, vectors, eigenvalues, LU/QR/SVD
4. **Statistics** — Distributions, hypothesis tests, Monte Carlo
5. **Differential Equations** — ODE/PDE, Euler, RK4/RK45, Laplace

### Advanced
6. **Complex Analysis** — Riemann zeta, Gamma function, zeros
7. **Analytic Number Theory** — π(x), Li(x), PNT, Dirichlet L-series
8. **Algebraic Number Theory** — Number fields, rings of integers, ideal classes
9. **Modular Forms** — SL(2,Z), Eisenstein series, Delta function, Hecke operators
10. **p-adic Numbers** — p-adic valuation, PAdicNumber class, Hensel's lemma
11. **L-Functions** — Dirichlet L-series, elliptic curves, BSD conjecture
12. **Fourier Analysis** — DFT/FFT, Fourier series, STFT, window functions

### Millennium Problems & Open Conjectures
13. **Riemann Hypothesis** — Zeros of ζ, critical line
14. **BSD Conjecture** — Elliptic curves, L-functions, ranks
15. **Goldbach** — Singular series, Hardy-Littlewood, circle method
16. **P vs NP** — Complexity classes, reductions, heuristics
17. **Yang-Mills** — Existence + mass gap
18. **Hodge Conjecture** — Algebraic geometry, cohomology
19. **Navier-Stokes** — Regularity, weak solutions, turbulence

### Group B — Further Open Problems (Papers 72–97)
- Freiman-PFR (proven: Gowers-Green-Manners-Tao 2023)
- Erdős-Ko-Rado (proven 1961)
- Lonely Runner, Graceful Tree
- Mandelbrot MLC, Conway knot, Whitehead asphericity
- Jacobian conjecture, Hadwiger n≥7, Andrews-Curtis, Beal
- Fontaine-Mazur, Lehmer-Mahler, Elliptic curve rank
- Odd perfect numbers, Mersenne primes, Bunyakovsky
- Donaldson, Gromov, Hartshorne, Uniformization
- Grothendieck standard conjectures, Kontsevich integral

---

## Papers Status (Build 199)

**97 papers** in 25 batches (EN + DE each), all fully audited and **print-ready**.

| Batch | Papers | Topics |
|-------|--------|--------|
| 1–4 | 1–20 | Giuga, Lehmer, Wilson, Sieve, Goldbach, Vinogradov, Circle method, Waring |
| 5–6 | 21–28 | Riemann Hypothesis, Elliptic curves + BSD |
| 7–8 | 29–36 | Collatz-Tao, Modular forms, abc, Navier-Stokes |
| 9–10 | 37–39 | Algebraic number theory, Iwasawa, Langlands |
| 11–13 | 40–51 | Giuga 4-prime, Group theory, Additive NT, Algorithmic, Twin primes |
| 14–15 | 52–59 | Algebraic topology, Differential geometry, Special functions, Transcendence |
| 16–18 | 60–71 | Hecke, Automorphic forms, Galois, Combinatorics, Yang-Mills, Hodge, P vs NP |
| 19–25 | 72–97 | Group B — 26 further open problems |

**Audit result (Build 198):**
- 30+ mathematical errors corrected (wrong coefficients, theorems vs. conjectures)
- 15+ LaTeX structural errors fixed
- 10+ bibliographic errors corrected
- All papers updated to current research state (2023–2024)

---

## Usage Examples

### Algebra
```python
from algebra import Polynomial, is_prime, euler_phi

p = Polynomial([1, -5, 6])  # x² - 5x + 6
print(p.evaluate(2))        # 0 (root)
print(is_prime(97))         # True
print(euler_phi(12))        # 4
```

### Number Theory
```python
from analytic_number_theory import prime_pi

print(prime_pi(1000))  # 168 (exact)
```

### Modular Forms
```python
from modular_forms import delta_function

tau_vals = delta_function(coefficients=10)
print(tau_vals[1])   # 1  (τ(1) = 1)
print(tau_vals[2])   # -24 (τ(2) = -24)
```

### Proof Attempts
```python
from beweisversuche import verify_no_3prime_giuga

result = verify_no_3prime_giuga()
print(result)  # Proven: True
```

---

## Heavy Compute Results

| Conjecture | Range | Result | Duration |
|-----------|-------|--------|----------|
| Erdős-Straus | p ≤ 10⁹ | 0 counterexamples (50,847,534 primes) | 401s |
| Giuga 4-prime | product ≤ 10¹⁵ | 0 solutions | 9.12h (20 cores) |
| Kurepa | p ≤ 10⁷ | 664,578 primes, 0 violations | — |
| Lehmer τ | n ≤ 50,000 | τ(n) ≠ 0 for all n | — |
| Brocard-Ramanujan | n ≤ 100,000 | Only {4, 5, 7} | — |

---

*Author: Michael Fuhrmann*
*Build: 199 | Date: 2026-03-13 | Tests: 10,403+ | Papers: 97 (all print-ready)*
*Generated by [Claude Code](https://claude.ai/claude-code) — see Disclaimer above*
