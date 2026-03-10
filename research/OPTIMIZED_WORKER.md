# OPTIMIZED_WORKER.md - Mathematik-Spezialist
## Aktualisiert: 2026-03-10 (Build 12)

## Was ich über dieses Projekt weiß

### Ziel
Ich werde zum Mathematik-Spezialisten ausgebildet. Das Projekt dient als
Wissensbasis, Übungsumgebung und Dokumentationssammlung für mathematische Konzepte.
Jedes Modul enthält Algorithmen mit ausführlicher mathematischer Dokumentation,
sodass der Code auch als Lernmaterial dient. Das Fernziel ist die Untersuchung
der Millennium-Probleme (insbesondere Riemann-Hypothese und Goldbach-Vermutung).

### Implementierter Stack (Build 12)
- **Python 3.13** mit sympy 1.14, numpy 2.2, scipy 1.17, matplotlib 3.10, pytest 9.0, mpmath
- **14 Python-Module**, 641+ Tests

| Modul | Kernfunktionen |
|-------|---------------|
| algebra.py | Polynome, Gleichungen, Zahlentheorie, Diophant, Reziprozität, RSA |
| analysis.py | Differentiation, Integration, Grenzwerte, Partialbrüche |
| linear_algebra.py | Matrizen, LU/QR/SVD/Givens, Eigenwerte/-vektoren |
| statistics_math.py | Verteilungen, Hypothesentests, Bayes, Monte-Carlo |
| ode.py | Euler, RK4/RK45, Laplace/inverse Laplace |
| complex_analysis.py | ζ(s), Γ(z), ξ, Z-Funktion, Nullstellen |
| analytic_number_theory.py | π(x), Li(x), Λ(n), θ(x), Dirichlet |
| proof_theory.py | Collatz, Goldbach+Kreismethode, CRT, Miller-Rabin |
| fourier.py | DFT/FFT (Cooley-Tukey), STFT, Fensterfunktionen |
| numerical_methods.py | Interpolation, BFGS, Simplex, Gradientenverfahren |
| modular_forms.py | SL(2,Z), Eisenstein, Δ, j-Invariante, Hecke |
| p_adic.py | p-adische Zahlen, Hensel-Lift, Ostrowski |
| visualization.py | 2D/3D-Plotter, Vektorfeld, Fraktale |
| **millennium_problems.py** | **7 Millennium-Probleme: RH, Goldbach, P≠NP, NS, YM, BSD** |

### Mathematische Erkenntnisse

#### Numerische Stabilität
- 1. numerische Ableitung: optimales h ≈ ε^(1/3) ≈ 6×10⁻⁶
- 2. Ableitung: h ≈ ε^(1/4) ≈ 1.2×10⁻⁴ (wegen h² im Nenner)
- Finite Differenzen k-ter Ordnung für k > 8 INSTABIL → SymPy verwenden
- QR via Householder stabiler als Gram-Schmidt; SVD robustester Weg für Eigenvektoren

#### Komplexe Analysis / Zeta-Funktion
- ζ(s) für Re(s) > 1: Euler-Maclaurin-Formel
- ζ(s) für 0 < Re(s) ≤ 1: η(s)/(1-2^{1-s}) mit Euler-Knopp (60 Terme → Maschinengenauigkeit)
- ζ(s) für Re(s) ≤ 0: Funktionalgleichung (Spiegelung an Re=1/2)
- ξ-Symmetrie ξ(s) = ξ(1-s) auf ~10⁻¹⁵ verifiziert

#### Modulformen
- Fundamentalbereich: |τ| > 1, |Re(τ)| ≤ 1/2, Im(τ) > 0
- SL(2,Z) erzeugt durch S: τ→-1/τ und T: τ→τ+1
- Δ(τ) = q∏(1-qⁿ)²⁴, q = e^{2πiτ}, konvergiert nur für Im(τ) > 0
- Shimura-Taniyama-Wiles: jede elliptische Kurve /ℚ ist modular

#### p-adische Zahlen
- Ostrowski: Alle nicht-trivialen Absolutwerte auf ℚ sind |·|∞ oder |·|p
- Produktformel: |n|∞ · ∏_p |n|p = 1 für alle n ≠ 0 ∈ ℚ
- Hensel-Lifting: Wurzeln mod p → mod p^k (Newton in ℤ_p)

#### Diophantische Gleichungen
- Lineares Diophant: lösbar ⟺ gcd(a,b) | c
- Pell x²-Dy²=1: Fundamentallösung via Kettenbruch √D
- Zwei-Quadrate-Satz: n=a²+b² ⟺ alle Primteiler ≡ 3 (mod 4) in gerader Potenz

### Strategische Ausrichtung: Millennium-Probleme

#### Riemann-Hypothese
- Werkzeuge: ζ(s), ξ(s), Z(t), N(T)-Formel, Nullstellensuche
- Status: >10^13 Nullstellen auf Re=1/2 verifiziert (empirisch)

#### Goldbach-Vermutung
- Werkzeuge: Goldbach-Zerlegung, Kreismethode, Hardy-Littlewood-Schätzung
- Singuläre Reihe S(n) quantifiziert Erwartungswert der Zerlegungen

### Brainstorming: Verbindungen

**Hardy-Littlewood ↔ Goldbach**: Kreismethode approximiert r₂(n) via Integral
über den Einheitskreis (Hauptbogen = Farey-Folge + Nebenbogen).

**Modulformen ↔ Zahlentheorie**: Ramanujan-Tau τ(n) aus Δ-Koeffizienten;
BSD-Vermutung verbindet L-Funktionen elliptischer Kurven mit Modulformen.

**p-adische ↔ Primzahlen**: p-adische L-Funktionen verallgemeinern Dirichlet-L-Reihen.

**Fourier ↔ Zahlentheorie**: DFT auf ℤ/nℤ; Dirichlet-Reihen auf ℤ; Modulformen auf ℝ.

### Build 11-12 Erkenntnisse (2026-03-10)

#### millennium_problems.py – Implementierte Werkzeuge

**Riemann-Hypothese (Build 12):**
- `riemann_zeros_mpmath(n, dps)`: mpmath.zetazero(k) liefert Nullstellen auf >50 Stellen genau
- Erste bekannte Nullstelle: 1/2 + 14.134725...i
- >10^13 Nullstellen verifiziert, Montgomery-Odlyzko: Nullstellenabstände ~ GUE-Statistik
- Gram-Gesetz gilt nur "oft" (nicht immer): Gram-Versagen bekannt ab g_126

**Goldbach (Build 12):**
- Schwache Goldbach (Helfgott 2013) BEWIESEN für alle ungeraden n > 5
- Starke Goldbach: Bis 4×10^18 verifiziert (Oliveira e Silva 2014)
- Hardy-Littlewood r₂(n)-Schätzung: C₂ ≈ 0.6601618 (Produktformel konvergiert langsam)

**P vs NP (Build 12):**
- 3-SAT Backtracking: NP-vollständig (Cook-Levin 1971)
- TSP Brute-Force: (n-1)! Permutationen, nur für n ≤ 12 praktikabel
- Greedy Nearest-Neighbor: O(n²), Approximationsverhältnis O(log n)
- Komplexitätsklassen-Vergleich: n=100 → n²=10⁴, 2^100≈10^30, 100!≈10^157

**Navier-Stokes (Build 12):**
- Chorin-Projektion (Druckpoisson + Geschwindigkeitskorrektur)
- 2D-NS: Glattheit bewiesen (Ladyzhenskaya 1969); 3D: offen
- CFL-Bedingung: dt ≤ min(dx,dy)² / (4ν) für Stabilität
- Enstrophie ε = 1/2 ∫|∇×u|² als Wirbelmaß

**Yang-Mills (Build 12):**
- Gitter-Eichtheorie (Wilson 1974): Links = U(1)-Phasen
- Metropolis-Monte-Carlo für Thermalisierung der Eichfelder
- Area-Law W(r) ~ exp(-σ·r²) als numerische Evidenz für Massenlücke

**BSD-Vermutung (Build 12):**
- Elliptische Kurven über F_p: Hasse-Weil-Schranke |#E(F_p)-(p+1)| ≤ 2√p
- BSD-Produkt-Schätzung: Π_{p≤X} (#E(F_p)/p) ≈ C·(log X)^r
- Beweis: Nur für Rang 0,1 (Coates-Wiles 1977, Kolyvagin 1990)

#### Numerische Erkenntnisse (Build 11-12)
- mpmath.zetazero(k): Deutlich schneller als selbst implementierte Nullstellensuche
- Upwind-Schema für Advektionsterm in NS: stabiler als zentrale Differenzen
- Gram-Gesetz Erfolgsquote bei n=20: typisch ~75-85% (nicht 100%)
- Wilson-Loop bei g=1 für 4×4-Gitter: W(r=1) ≈ 0.5-0.8 je nach Kopplung

### Brainstorming: Neue Verbindungen (Build 12)

**Montgomery-Odlyzko ↔ Zufallsmatrizen**: Nullstellenabstände der ζ-Funktion haben
exakt dieselbe statistische Verteilung wie Eigenwerte zufälliger hermitescher Matrizen
(GUE). Dies ist völlig unerwartet und unverstanden – eine tiefe Verbindung zwischen
Zahlentheorie und mathematischer Physik.

**Yang-Mills ↔ Navier-Stokes**: Beide sind PDE-Existenzprobleme (Millennium-Probleme),
beide involvieren nichtlineare Terme, beide haben in niedrigeren Dimensionen Lösungen.
Strukturelle Analogie: Yang-Mills Wirkung ~ kinetische Energie in NS.

**BSD ↔ Modulformen**: Wiles' Beweis von Fermats Letztem Satz nutzte Shimura-Taniyama
(elliptische Kurven sind modular). BSD-L-Funktion = spezielle Modular-L-Funktion.

**3-SAT ↔ Riemann-Hypothese**: Könnten NP-harte Probleme einen Zusammenhang mit
analytischen Nullstellenproblemen haben? Freie Spekulation, aber interessant für Forschung.

### Nächste Entwicklungsschritte (Build 13+)
1. Jupyter-Integration: Interaktiver REPL-Modus für Millennium-Explorationen
2. p-adische L-Funktionen: Kubota-Leopoldt, Iwasawa-Algebra
3. Explizite Formel ψ(x) = x - Σ_{ρ} x^ρ/ρ - log(2π) vertiefen
4. BSD: Elliptische-Kurven-Arithmetik (Gruppenstruktur, Punktaddition)
5. Yang-Mills: SU(2) statt U(1) auf Gitter (2×2-Matrizen als Links)
6. NumPy-Vektorisierung der Kern-Schleifen (10-100x Speed)
7. Modulformen Vertiefung: Cusp-Formen (Γ₀(N)), Theta-Reihen θ(τ) = Σqⁿ²
