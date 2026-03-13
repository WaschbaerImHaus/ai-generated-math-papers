# specialist-maths — Mathematik-Spezialist

**Autor:** Michael Fuhrmann
**Build:** 198 | **Stand:** 2026-03-13 | **Tests:** 10.403+

---

## Überblick

Dieses Projekt dient der systematischen Entwicklung mathematischer Expertise.
Es enthält Python-Bibliotheken für alle Gebiete der Mathematik, vollständig
dokumentiert, getestet und als Lernmaterial geeignet.

**Fernziel:** Untersuchung der Millennium-Probleme (RH, Goldbach, BSD, P≠NP).

---

## Anforderungen

- Python 3.13+
- sympy >= 1.14
- numpy >= 2.2
- scipy >= 1.17
- matplotlib >= 3.10
- mpmath, numba, networkx, gmpy2, pysat, z3-solver

## Installation

```bash
pip install sympy numpy scipy matplotlib mpmath numba networkx gmpy2 python-sat z3-solver pytest pytest-xdist --break-system-packages
```

## Tests ausführen

```bash
# Schnell (parallel):
python3 -m pytest tests/ -n auto

# Ausführlich:
python3 -m pytest tests/ -v --timeout=60
```

---

## Projektstruktur

```
specialist-maths/
├── src/                               # 120+ Python-Module
│   ├── algebra.py                     # Polynome, Gleichungen, Zahlentheorie
│   ├── analysis.py                    # Analysis, Differentiation, Integration
│   ├── linear_algebra.py              # Lineare Algebra, LU/QR/SVD
│   ├── statistics_math.py             # Statistik, Verteilungen, Tests
│   ├── ode.py                         # Differentialgleichungen, Laplace
│   ├── pde.py                         # Partielle Differentialgleichungen
│   ├── complex_analysis.py            # Komplexe Analysis, Riemann-Zeta
│   ├── analytic_number_theory.py      # Analytische Zahlentheorie
│   ├── proof_theory.py                # Beweistheorie, Millennium-Probleme
│   ├── fourier.py                     # Fourier-Transformation (FFT)
│   ├── numerical_methods.py           # Interpolation, Optimierung, Simplex
│   ├── modular_forms.py               # Modulformen (SL2Z, Eisenstein, Delta)
│   ├── p_adic.py                      # p-adische Zahlen
│   ├── l_functions.py                 # L-Funktionen, Dirichlet-Reihen
│   ├── visualization.py               # Visualisierung (2D/3D, Fraktale)
│   ├── algebraic_number_theory.py     # Algebraische Zahlentheorie
│   ├── galois_theory.py               # Galois-Theorie
│   ├── ergodic_theory.py              # Ergodentheorie
│   ├── beweisversuche.py              # Giuga/Lehmer/Kurepa/Erdős-Straus
│   ├── millennium_problems.py         # Alle 7 Millennium-Probleme
│   │
│   ├── # Gruppe A — Erweiterte Vermutungen (Build 125–131):
│   ├── giuga_4prim.py                 # Giuga 4-Prim-Vermutung
│   ├── erdos_straus_ext.py            # Erdős-Straus (verifiziert bis 10^9)
│   ├── brocard_extension.py           # Brocard-Ramanujan
│   ├── kurepa_ext.py                  # Kurepa-Linksexponential (bis 10^7)
│   ├── lehmer_tau.py                  # Lehmer τ-Funktion
│   ├── bruns_constant.py              # Bruns Konstante
│   ├── schur_numbers.py               # Schur-Zahlen
│   ├── debruijn_newman.py             # de Bruijn-Newman-Konstante
│   ├── hadwiger_nelson.py             # Hadwiger-Nelson-Problem
│   ├── cohen_lenstra.py               # Cohen-Lenstra-Heuristiken
│   ├── euler_gamma_irrationality.py   # Euler-Mascheroni γ
│   ├── zeta_odd_values.py             # ζ(2k+1) Irrationalität
│   ├── pillai_chowla.py               # Pillai-Chowla-Vermutung
│   ├── mahler_measure.py              # Lehmer-Mahler-Maß
│   ├── perfect_numbers.py             # Vollkommene Zahlen
│   ├── waring_goldbach.py             # Waring-Goldbach-Problem
│   ├── frankl_union_closed.py         # Frankl Union-Closed (Gilmer/Chase-Lovász)
│   ├── beal_conjecture.py             # Beal-Vermutung
│   └── build.txt                      # Build-Nummer (aktuell: 198)
│
├── tests/                             # 10.403+ Unit-Tests (TDD, pytest-xdist)
├── papers/                            # 97 LaTeX-Papers (Batches 1–25)
│   └── reviewed/                      # Batches 1–9 fertig auditiert
├── reviews/                           # Audit-Reports
├── research/                          # Mathematische Recherche & Notizen
├── dev-log/                           # Entwicklungs-Tagebuch
├── debugging/                         # Debug-Skripte
├── build/                             # Ausführbare Dateien
├── claude-generated/                  # Weitere generierte Dateien
├── FEATURES.md                        # Feature-Liste
├── BUGS.md                            # Bug-Tracker
├── OPTIMIZE.md                        # Optimierungsideen
└── SECURITY_RISKS.md                  # Sicherheitsrisiken
```

---

## Mathematische Themengebiete

### Grundlagen
1. **Algebra**: Polynomrechnung, lineare/quadratische Gleichungen, Zahlentheorie
2. **Analysis**: Differentiation, Integration, Taylor-Reihen, Nullstellensuche
3. **Lineare Algebra**: Matrizen, Vektoren, Eigenwerte, LU/QR/SVD
4. **Statistik**: Verteilungen, Hypothesentests, Monte-Carlo
5. **Differentialgleichungen**: Euler, RK4/RK45, Laplace-Transformation (ODE + PDE)

### Fortgeschritten
6. **Komplexe Analysis**: Riemann-Zeta, Gamma-Funktion, Nullstellen
7. **Analytische Zahlentheorie**: π(x), Li(x), Primzahlsatz, Dirichlet-L-Reihen
8. **Algebraische Zahlentheorie**: Zahlkörper, Ganzheitsringe, Idealklassen
9. **Modulformen**: SL(2,Z), Eisenstein-Reihen, Delta-Funktion, Hecke-Operatoren
10. **p-adische Zahlen**: p-adische Bewertung, PAdicNumber, Hensels Lemma, Iwasawa
11. **L-Funktionen**: Dirichlet-L-Reihen, elliptische Kurven, BSD-Vermutung
12. **Fourier-Analysis**: DFT/FFT, Fourier-Reihen, STFT, Fensterfunktionen

### Millennium-Probleme & offene Vermutungen
13. **Riemann-Hypothese**: Nullstellen der ζ-Funktion, Kritische Linie
14. **BSD-Vermutung**: Elliptische Kurven, L-Funktionen, Ränge
15. **Goldbach**: Singuläre Reihen, Hardy-Littlewood, Kreismethode
16. **P vs NP**: Komplexitätsklassen, Reduktionen, Heuristiken
17. **Yang-Mills**: Existenz + Masselücke (physikalisch-mathematisch)
18. **Hodge-Vermutung**: Algebraische Geometrie, Kohomologie
19. **Navier-Stokes**: Regularität, schwache Lösungen, Turbulenz

### Gruppe B — Weitere offene Probleme (Papers 72–97)
- Freiman-PFR (Gowers-Green-Manners-Tao 2023 bewiesen)
- Erdős-Ko-Rado (bewiesen 1961)
- Lonely Runner, Graceful Tree
- Mandelbrot MLC, Conway-Knoten, Whitehead-Aspärizität
- Jacobian-Vermutung, Hadwiger n≥7, Andrews-Curtis, Beal
- Fontaine-Mazur, Lehmer-Mahler, Elliptische Kurven Rang
- Ungerade vollkommene Zahlen, Mersenne-Primzahlen, Bunyakovsky
- Donaldson, Gromov, Hartshorne, Uniformisierung
- Grothendieck Standard-Vermutungen, Kontsevich-Integral

---

## Papers-Status (Build 198)

**97 Papers** in 25 Batches (EN + DE), alle vollständig auditiert und **DRUCKREIF**.

| Batch | Papers | Themen |
|-------|--------|--------|
| 1–4 | 1–20 | Giuga, Lehmer, Wilson, Sieb, Goldbach, Vinogradov, Kreismethode, Waring |
| 5–6 | 21–28 | Riemann-Hypothese, Elliptische Kurven + BSD |
| 7–8 | 29–36 | Collatz-Tao, Modulformen, abc, Navier-Stokes |
| 9–10 | 37–39 | Algebraische Zahlentheorie, Iwasawa, Langlands |
| 11–13 | 40–51 | Giuga 4-Prim, Gruppentheorie, Additive NT, Algorithmisch, Zwillingsprimzahlen |
| 14–15 | 52–59 | Algebraische Topologie, Differentialgeometrie, Spezielle Fkt., Transzendenz |
| 16–18 | 60–71 | Hecke, Automorphe Formen, Galois, Kombinatorik, Yang-Mills, Hodge, P vs NP |
| 19–25 | 72–97 | Gruppe B — 26 weitere offene Probleme |

**Audit-Ergebnis Build 198:**
- 30+ mathematische Fehler korrigiert (falsche Koeffizienten, Theoreme vs. Vermutungen)
- 15+ LaTeX-Strukturfehler behoben
- 10+ Bibliographische Fehler korrigiert
- Alle Papers auf aktuellen Forschungsstand gebracht (2023–2024)

---

## Verwendungsbeispiele

### Algebra
```python
from algebra import Polynomial, is_prime, euler_phi

p = Polynomial([1, -5, 6])  # x² - 5x + 6
print(p.evaluate(2))        # 0 (Nullstelle)
print(p.derivative())       # 2x - 5
print(is_prime(97))         # True
print(euler_phi(12))        # 4
```

### Zahlentheorie
```python
from analytic_number_theory import prime_pi, prime_counting_approx

print(prime_pi(1000))               # 168 (exakt)
print(prime_counting_approx(1000))  # ≈ 168 (Li-Approximation)
```

### Modular Forms
```python
from modular_forms import ModularForm, delta_function, j_invariant

# Ramanujan-Delta-Funktion τ(n)
tau_vals = delta_function(coefficients=10)
print(tau_vals[1])   # 1 (τ(1) = 1)
print(tau_vals[2])   # -24 (τ(2) = -24)
```

### Beweisversuche
```python
from beweisversuche import verify_no_3prime_giuga, verify_erdos_straus

# Kein Giuga-Pseudoprim mit ≤3 Primfaktoren
result = verify_no_3prime_giuga()
print(result)  # Bewiesen: True

# Erdős-Straus: 4/p = 1/a + 1/b + 1/c
# (verifiziert bis p ≤ 10^9)
```

---

## Heavy-Compute-Ergebnisse

| Vermutung | Bereich | Ergebnis | Dauer |
|-----------|---------|----------|-------|
| Erdős-Straus | p ≤ 10⁹ | 0 Gegenbeispiele (50.847.534 Primzahlen) | 401s |
| Giuga 4-Prim | Produkt ≤ 10¹⁵ | 0 Lösungen | 9,12h (20 Kerne) |
| Kurepa | p ≤ 10⁷ | 664.578 Primzahlen, 0 Verletzungen | — |
| Lehmer τ | n ≤ 50.000 | τ(n) ≠ 0 für alle n | — |
| Brocard-Ramanujan | n ≤ 100.000 | Nur {4, 5, 7} | — |

---

*Autor: Michael Fuhrmann*
*Build: 198 | Stand: 2026-03-13 | Tests: 10.403+ | Papers: 97 (alle DRUCKREIF)*
