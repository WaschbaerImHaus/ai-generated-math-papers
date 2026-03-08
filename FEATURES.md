# FEATURES.md - specialist-maths

## Implementierte Features

### Algebra (src/algebra.py) - Build 2
- [x] Polynomklasse mit Horner-Schema-Auswertung
- [x] Polynom-Addition, Subtraktion, Multiplikation
- [x] Polynom-Ableitung (Potenzregel)
- [x] Polynom String-Darstellung
- [x] Linearer Gleichungslöser (mit Sonderfällen: keine/unendlich viele Lösungen)
- [x] Quadratischer Gleichungslöser (inkl. komplexe Wurzeln)
- [x] Euklidischer Algorithmus (ggT, kgV)
- [x] Erweiterter euklidischer Algorithmus (Bezout-Koeffizienten)
- [x] Modulares Inverses
- [x] Primzahltest (optimiert: 6k±1-Methode)
- [x] Primfaktorzerlegung
- [x] Eulersche Phi-Funktion

### Analysis (src/analysis.py) - Build 2
- [x] Numerische Ableitung 1. Ordnung (zentraler Differenzenquotient, h≈ε^1/3)
- [x] Numerische Ableitung 2. Ordnung (mit optimaler Schrittweite h≈ε^1/4)
- [x] Numerische Integration (Simpson-Regel, O(h⁴) Genauigkeit)
- [x] Newton-Raphson-Verfahren (Nullstellensuche, quadratische Konvergenz)
- [x] Bisektionsverfahren (Nullstellensuche, garantierte Konvergenz)
- [x] Taylor-Reihen-Entwicklung (symbolisch via SymPy)

### Lineare Algebra (src/linear_algebra.py) - Build 2
- [x] Vektor: Skalarprodukt, Kreuzprodukt, Norm, Normierung
- [x] Vektor: Addition, Subtraktion, Skalarmultiplikation
- [x] Matrix: Erstellen, Identity, Zeros
- [x] Matrix: Multiplikation (matmul)
- [x] Matrix: Transponierung, Spur
- [x] Matrix: Determinante (Gauss mit Pivotsuche)
- [x] Matrix: Inverse (Gauss-Jordan)
- [x] Matrix: LGS lösen (Gauss mit Rücksubstitution)
- [x] Matrix: Eigenwerte (analytisch 2×2, QR-Iteration n×n)
- [x] Gram-Schmidt-Orthogonalisierung (mit ONB-Option)

### Statistik & Wahrscheinlichkeit (src/statistics_math.py) - Build 2
- [x] Deskriptive Statistik: mean, median, mode, variance, std_dev
- [x] Lage- und Streuungsmaße: quartiles, iqr, skewness, kurtosis
- [x] Normalverteilung: PDF, CDF, PPF (Quantilfunktion)
- [x] Binomialverteilung: PMF, CDF
- [x] Poisson-Verteilung: PMF, CDF
- [x] Hypothesentests: Ein-/Zwei-Stichproben-t-Test, Chi-Quadrat-Test
- [x] Bayes-Theorem
- [x] Monte-Carlo-Simulation (π-Schätzung)

### Differentialgleichungen (src/ode.py) - Build 2
- [x] Euler-Verfahren (explizit, 1. Ordnung)
- [x] Runge-Kutta 4. Ordnung (klassisch, O(h⁴))
- [x] Runge-Kutta-Fehlberg RK45 (adaptiv, Schrittweitensteuerung)
- [x] Lineare ODE mit konstanten Koeffizienten (Zustandsraum + RK4)
- [x] Numerische Laplace-Transformation (Simpson-Regel)
- [x] Numerische inverse Laplace-Transformation (Stehfest-Algorithmus)

## Implementierte Features (Fortsetzung)

### Beweistheorie & Offene Vermutungen (src/proof_theory.py) - Build 4
- [x] Collatz-Vermutung: Folge, Stoppzeit, Maximalwert, Bereichs-Verifikation
- [x] Goldbach-Vermutung: Zerlegung, alle Zerlegungen, Bereichs-Verifikation
- [x] Zwillingsprimzahlen: Suche, Zählung, Hardy-Littlewood-Konstante
- [x] Riemann-Zeta-Funktion: Partialsumme, Euler-Maclaurin, numerische Nullstellensuche
- [x] Sieb des Eratosthenes (O(n log log n))
- [x] Miller-Rabin-Primzahltest (deterministisch bis 3.3×10^24)
- [x] Legendre-Symbol (quadratische Reste)
- [x] Jacobi-Symbol (verallgemeinertes Legendre)
- [x] Chinesischer Restsatz (CRT)
- [x] ProofByInduction-Klasse (Basisfall, empirische Verifikation)
- [x] Vermutungsstatus-Report (alle 7 Millennium-Probleme + weitere)

### Komplexe Analysis (src/complex_analysis.py) - Build 5
- [x] Gamma-Funktion Γ(z) via Lanczos-Approximation (15 Stellen Genauigkeit)
- [x] log_gamma(z) via Stirling-Reihe (numerisch stabil für großes |z|)
- [x] Vollständige Riemann-Zeta ζ(s) über analytische Fortsetzung (3-Bereich-Strategie)
- [x] ξ-Funktion (symmetrisch: ξ(s)=ξ(1-s) auf ~10^-15 verifiziert)
- [x] Riemann-Siegel-Z-Funktion (reellwertig auf kritischer Geraden)
- [x] Nullstellensuche via find_zeta_zeros() mit Bisektionsverfeinerung
- [x] N(T)-Formel (Riemann-von-Mangoldt: Anzahl Nullstellen bis |Im| ≤ T)
- [x] Cauchy-Integralsatz (numerische Mittelwerteigenschaft)
- [x] Residuensatz (numerisch)

### Fourier-Analysis (src/fourier.py) - Build 6
- [x] DFT/IDFT (diskrete Fourier-Transformation, O(N²))
- [x] FFT/IFFT (Cooley-Tukey radix-2, O(N log N))
- [x] fft_padded (Zero-Padding auf nächste Zweierpotenz)
- [x] Fourier-Reihen: fourier_coefficients, fourier_series_eval
- [x] Fensterfunktionen: Hanning, Hamming, Blackman, apply_window
- [x] Leistungsdichtespektrum (PSD), dominant_frequency
- [x] Kurzzeit-Fourier-Transformation (STFT)

### Numerische Methoden (src/numerical_methods.py) - Build 6
- [x] Lagrange-Interpolation
- [x] Newton-Interpolation (dividierte Differenzen, inkrementell)
- [x] Kubische Splines (natürlich, Tridiagonal-LGS, Thomas-Algorithmus)
- [x] Gradient Descent (mit numerischem Gradient-Fallback)
- [x] Goldener-Schnitt-Suche (unimodale Funktionen)
- [x] Numerischer Gradient (zentrale Differenzen)
- [x] Simplex-Algorithmus (lineare Programmierung, Minimierungsform)

### Lineare Algebra Erweiterungen (src/linear_algebra.py) - Build 6
- [x] LU-Zerlegung (Doolittle mit Teilpivotisierung, PA=LU)
- [x] QR-Zerlegung (Householder-Reflexionen, numerisch stabil)
- [x] Singulärwertzerlegung SVD (A = U·Σ·Vᵀ)
- [x] Matrixrang (via SVD, robust gegen Fast-Singularität)
- [x] Konditionszahl κ = σ_max/σ_min

### Build 8 Erweiterungen (2026-03-08)
- [x] Eigenvektoren via SVD-Kernel (linear_algebra.py) – 5 Tests
- [x] BFGS quasi-Newton Optimierer (numerical_methods.py) – 5 Tests
- [x] RSA-Kryptosystem: keygen/encrypt/decrypt (algebra.py) – 5 Tests
- [x] Umfassende .md-Dokumentation für alle 10 Module (Build 7)
- **Gesamttests: 333/333 grün**

### Build 9 Erweiterungen (2026-03-08)
- [x] Kreismethode (Hardy-Littlewood) für Goldbach (proof_theory.py) – 20 Tests
- [x] Modulformen Grundlagen → neue modular_forms.py (SL2Z, Eisenstein, Delta, j-Invariante, Hecke, Shimura-Taniyama) – 40 Tests
- [x] p-adische Zahlen → neue p_adic.py (Bewertung, Norm, PAdicNumber, Hensel, Exp/Log, Ostrowski) – 47 Tests
- [x] Diophantische Gleichungen (algebra.py): linear, Pell, Pythagoras, Zwei-Quadrate, Markov – 66 Tests
- [x] Quadratisches Reziprozitätsgesetz vollständig (is_quadratic_residue, Tonelli-Shanks, Cipolla) – Teil der 66 Tests
- [x] Analysis: Grenzwerte symbolisch, L'Hôpital, Partialbrüche, Uneigentliche Integrale, Cauchy-PV – 40 Tests
- [x] Givens-Rotationen für QR-Zerlegung (linear_algebra.py) – 27 Tests
- [x] Vollständiges Visualisierungsmodul → neue visualization.py (2D/3D, Vektorfeld, Phasenraum, Fraktale) – 48 Tests
- **Gesamttests: 641/641 grün**

### Analytische Zahlentheorie (src/analytic_number_theory.py) - Build 5
- [x] Primzahlzählfunktion π(x) (exakt via Sieb)
- [x] Logarithmisches Integral Li(x) (Gauß/Riemann-Approximation)
- [x] Vergleich π(x) vs Li(x) vs x/ln(x) (Primzahlsatz mit Fehlerterm)
- [x] Von-Mangoldt-Funktion Λ(n)
- [x] Tschebyschow-Funktionen θ(x), ψ(x)
- [x] Dirichlet-Charaktere und L-Funktionen
- [x] Primzahlen in arithmetischen Restklassen (Dirichlet-Theorem)
- [x] Semiprime und Chen-Zerlegung (fast-Primzahlen)
- [x] Primzahl-Lücken und Cramér-Vermutung
- [x] Riemanns Explizite Formel π(x) ≈ Li(x) - Σ Li(x^ρ) (Nullstellen-Approximation)
- [x] Sieb-Methoden: Brun-Sieb, Selberg-Sieb, Chen-Zerlegung

## Geplante Features

### Beweistheorie Erweiterungen (Priorität: HÖCHSTE – Fernziel: Millennium-Beweise)
- [x] Kreismethode (Hardy-Littlewood) für Goldbach – Build 9
- [x] Modulformen (Grundlagen) – Build 9
- [x] p-adische Zahlen – Build 9

### Zahlentheorie (Priorität: hoch)
- [x] RSA-Verschlüsselung (keygen, encrypt, decrypt) – Build 8
- [x] Diophantische Gleichungen (linear, Pell, Pythagoras, Zwei-Quadrate, Markov) – Build 9
- [x] Quadratisches Reziprozitätsgesetz (vollständig) – Build 9

### Analysis Erweiterungen (Priorität: mittel)
- [x] Grenzwertberechnung (symbolisch via SymPy) – Build 9
- [x] Fourier-Reihen und -Transformation (DFT/FFT) – Build 6
- [x] Partialbruchzerlegung – Build 9
- [x] Uneigentliche Integrale – Build 9

### Lineare Algebra Erweiterungen (Priorität: mittel)
- [x] LU-Zerlegung (Doolittle mit Teilpivotisierung) – Build 6
- [x] QR-Zerlegung (Householder) – Build 6
- [x] Singulärwertzerlegung (SVD) – Build 6
- [x] Eigenvektoren via SVD-Kernel (zu vorhandenen Eigenwerten) – Build 8
- [x] Givens-Rotationen für QR – Build 9

### Numerische Methoden (Priorität: mittel)
- [x] Interpolation (Lagrange, Newton, kubische Splines) – Build 6
- [x] Optimierung (Gradient Descent, Goldener Schnitt) – Build 6
- [x] Lineare Programmierung (Simplex-Algorithmus) – Build 6
- [x] BFGS quasi-Newton Optimierung (mit Armijo-Line-Search) – Build 8

### Visualisierung (Priorität: niedrig)
- [x] Funktionsplotter 2D (matplotlib) – Build 9
- [x] Funktionsplotter 3D (matplotlib 3D) – Build 9
- [x] Vektorfelddarstellung – Build 9
- [x] Phasenraum-Visualisierung (für ODE) – Build 9
- [x] Fraktal-Generator (Mandelbrot, Julia, Sierpinski, Newton) – Build 9

## Offene Features (nach Build 9)
- [ ] Modulformen Vertiefung: Cusp-Formen, Theta-Reihen
- [ ] Hardy-Littlewood-Kreismethode Vertiefung: Farey-Folgen
- [ ] Jupyter-Notebook-Integration / REPL-Modus
- [ ] LaTeX-Export für Ergebnisse
- [ ] Property-Based Testing (Hypotheses-Bibliothek)
- [ ] NumPy-Vektorisierung der Kernalgorithmen
