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

## Geplante Features

### Beweistheorie Erweiterungen (Priorität: HÖCHSTE – Fernziel: Millennium-Beweise)
- [ ] Analytische Zahlentheorie (Primzahlsatz mit Fehlerterm, Dirichlet-L-Reihen)
- [ ] Riemann-Zeta via Euler-Produkt und analytische Fortsetzung
- [ ] Funktionalgleichung ζ(s) = χ(s)ζ(1-s) implementieren
- [ ] Explizite Formel Riemann: π(x) = Li(x) - Σ Li(x^ρ) + ...
- [ ] N(T)-Formel: Anzahl der Nullstellen mit |Im(ρ)| ≤ T
- [ ] Sieb-Methoden: Brun-Sieb, Selberg-Sieb, Chen-Theorem
- [ ] Kreismethode (Hardy-Littlewood) für Goldbach
- [ ] Modulformen (Grundlagen)
- [ ] p-adische Zahlen

### Zahlentheorie (Priorität: hoch)
- [ ] RSA-Verschlüsselung (als Anwendungsbeispiel)
- [ ] Diophantische Gleichungen (lineare, quadratische)
- [ ] Quadratisches Reziprozitätsgesetz (vollständig)

### Analysis Erweiterungen (Priorität: mittel)
- [ ] Grenzwertberechnung (symbolisch via SymPy)
- [ ] Fourier-Reihen und -Transformation (DFT/FFT)
- [ ] Partialbruchzerlegung
- [ ] Uneigentliche Integrale

### Lineare Algebra Erweiterungen (Priorität: mittel)
- [ ] LU-Zerlegung (direkte Methode)
- [ ] QR-Zerlegung (Householder/Givens)
- [ ] Singulärwertzerlegung (SVD)
- [ ] Eigenvektoren (zu vorhandenen Eigenwerten)

### Numerische Methoden (Priorität: mittel)
- [ ] Interpolation (Lagrange, Newton, kubische Splines)
- [ ] Optimierung (Gradient Descent, BFGS)
- [ ] Lineare Programmierung (Simplex-Algorithmus)

### Visualisierung (Priorität: niedrig)
- [ ] Funktionsplotter 2D (matplotlib)
- [ ] Funktionsplotter 3D (matplotlib 3D)
- [ ] Vektorfelddarstellung
- [ ] Phasenraum-Visualisierung (für ODE)
- [ ] Fraktal-Generator (Mandelbrot, Julia-Menge)
