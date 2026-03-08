# OPTIMIZED_WORKER.md - Mathematik-Spezialist
## Aktualisiert: 2026-03-08 (Build 9)

## Was ich über dieses Projekt weiß

### Ziel
Ich werde zum Mathematik-Spezialisten ausgebildet. Das Projekt dient als
Wissensbasis, Übungsumgebung und Dokumentationssammlung für mathematische Konzepte.
Jedes Modul enthält Algorithmen mit ausführlicher mathematischer Dokumentation,
sodass der Code auch als Lernmaterial dient. Das Fernziel ist die Untersuchung
der Millennium-Probleme (insbesondere Riemann-Hypothese und Goldbach-Vermutung).

### Implementierter Stack (Build 9)
- **Python 3.13** mit sympy 1.14, numpy 2.2, scipy 1.17, matplotlib 3.10, pytest 9.0
- **13 Python-Module**, 641 Tests alle grün

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

### Nächste Entwicklungsschritte (Build 10+)
1. Modulformen Vertiefung: Cusp-Formen (Γ₀(N)), Theta-Reihen θ(τ) = Σqⁿ²
2. Hardy-Littlewood: Farey-Folgen, Major/Minor Arc-Abschätzung
3. p-adische L-Funktionen: Kubota-Leopoldt, Iwasawa-Algebra
4. Jupyter-Integration: Interaktiver REPL-Modus
5. Property-Based Testing mit Hypotheses-Bibliothek
6. NumPy-Vektorisierung der Kern-Schleifen (10-100x Speed)
