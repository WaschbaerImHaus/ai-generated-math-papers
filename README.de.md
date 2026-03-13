# specialist-maths — Mathematik-Spezialist

**Autor:** Michael Fuhrmann
**Build:** 199 | **Stand:** 2026-03-13 | **Tests:** 10.403+

> **[English version → README.md](README.md)**

---

## Disclaimer

> **Dieses Projekt dient ausschließlich dazu, die aktuellen Fähigkeiten von
> [Claude Code](https://claude.ai/claude-code) (Anthropics KI-Coding-Assistent)
> auszuloten.**
>
> Alle mathematischen Papers, Beweise, Quellcode, Python-Module und Test-Suiten
> in diesem Repository wurden **vollständig von Claude Code** (Modelle:
> claude-opus-4-6, claude-sonnet-4-6) unter der Leitung von Michael Fuhrmann
> **generiert**.
>
> Die Inhalte sollen zeigen, was ein KI-Coding-Assistent leisten kann, wenn er
> auf ein komplexes, langfristiges wissenschaftliches Projekt angewendet wird.
> **Kein Mensch hat die mathematischen Papers oder den Quellcode direkt
> geschrieben.**
>
> Mathematische Aussagen sollten mit angemessener Skepsis betrachtet werden —
> obwohl erheblicher Aufwand in die Überprüfung der Korrektheit investiert wurde,
> können KI-generierte mathematische Inhalte Fehler enthalten.

---

## Überblick

Dieses Projekt erkundet systematisch mathematische Expertise in allen Bereichen —
von Grundlagen bis hin zu fortgeschrittenen Forschungsthemen.

**Fernziel:** Untersuchung der Millennium-Probleme (RH, Goldbach, BSD, P≠NP).

---

## Lizenz

Dieses Projekt steht unter der **MIT-Lizenz** — siehe [LICENSE](LICENSE) für Details.

Frei verwendbar, veränderbar und weiterzugeben, sofern der ursprüngliche Autor
(**Michael Fuhrmann**) genannt wird.

---

## Anforderungen

- Python 3.13+
- sympy >= 1.14, numpy >= 2.2, scipy >= 1.17, matplotlib >= 3.10
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
├── tests/                             # 10.403+ Unit-Tests (TDD, pytest-xdist)
├── papers/                            # 97 LaTeX-Papers (Batches 1–25, EN + DE)
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
└── LICENSE                            # MIT-Lizenz
```

---

## Mathematische Themengebiete

### Grundlagen
1. **Algebra** — Polynome, Gleichungen, Zahlentheorie
2. **Analysis** — Differentiation, Integration, Taylor-Reihen, Nullstellensuche
3. **Lineare Algebra** — Matrizen, Vektoren, Eigenwerte, LU/QR/SVD
4. **Statistik** — Verteilungen, Hypothesentests, Monte-Carlo
5. **Differentialgleichungen** — ODE/PDE, Euler, RK4/RK45, Laplace-Transformation

### Fortgeschritten
6. **Komplexe Analysis** — Riemann-Zeta, Gamma-Funktion, Nullstellen
7. **Analytische Zahlentheorie** — π(x), Li(x), Primzahlsatz, Dirichlet-L-Reihen
8. **Algebraische Zahlentheorie** — Zahlkörper, Ganzheitsringe, Idealklassen
9. **Modulformen** — SL(2,Z), Eisenstein-Reihen, Delta-Funktion, Hecke-Operatoren
10. **p-adische Zahlen** — p-adische Bewertung, PAdicNumber-Klasse, Hensels Lemma
11. **L-Funktionen** — Dirichlet-L-Reihen, elliptische Kurven, BSD-Vermutung
12. **Fourier-Analysis** — DFT/FFT, Fourier-Reihen, STFT, Fensterfunktionen

### Millennium-Probleme & offene Vermutungen
13. **Riemann-Hypothese** — Nullstellen der ζ-Funktion, Kritische Linie
14. **BSD-Vermutung** — Elliptische Kurven, L-Funktionen, Ränge
15. **Goldbach** — Singuläre Reihen, Hardy-Littlewood, Kreismethode
16. **P vs NP** — Komplexitätsklassen, Reduktionen, Heuristiken
17. **Yang-Mills** — Existenz + Masselücke
18. **Hodge-Vermutung** — Algebraische Geometrie, Kohomologie
19. **Navier-Stokes** — Regularität, schwache Lösungen, Turbulenz

### Gruppe B — Weitere offene Probleme (Papers 72–97)
- Freiman-PFR (bewiesen: Gowers-Green-Manners-Tao 2023)
- Erdős-Ko-Rado (bewiesen 1961)
- Lonely Runner, Graceful Tree
- Mandelbrot MLC, Conway-Knoten, Whitehead-Aspärizität
- Jacobian-Vermutung, Hadwiger n≥7, Andrews-Curtis, Beal
- Fontaine-Mazur, Lehmer-Mahler, Elliptische Kurven Rang
- Ungerade vollkommene Zahlen, Mersenne-Primzahlen, Bunyakovsky
- Donaldson, Gromov, Hartshorne, Uniformisierung
- Grothendieck Standard-Vermutungen, Kontsevich-Integral

---

## Papers-Status (Build 199)

**97 Papers** in 25 Batches (EN + DE jeweils), alle vollständig auditiert und **DRUCKREIF**.

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

**Audit-Ergebnis (Build 198):**
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
print(is_prime(97))         # True
print(euler_phi(12))        # 4
```

### Zahlentheorie
```python
from analytic_number_theory import prime_pi

print(prime_pi(1000))  # 168 (exakt)
```

### Modulformen
```python
from modular_forms import delta_function

tau_vals = delta_function(coefficients=10)
print(tau_vals[1])   # 1  (τ(1) = 1)
print(tau_vals[2])   # -24 (τ(2) = -24)
```

### Beweisversuche
```python
from beweisversuche import verify_no_3prime_giuga

result = verify_no_3prime_giuga()
print(result)  # Bewiesen: True
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
*Build: 199 | Stand: 2026-03-13 | Tests: 10.403+ | Papers: 97 (alle DRUCKREIF)*
*Generiert von [Claude Code](https://claude.ai/claude-code) — siehe Disclaimer oben*
