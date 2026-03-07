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

## Geplante Features

### Zahlentheorie (Priorität: hoch)
- [ ] Primzahl-Sieb des Eratosthenes (für große Bereiche)
- [ ] Miller-Rabin-Primzahltest (für sehr große Zahlen)
- [ ] Chinesischer Restsatz (CRT)
- [ ] RSA-Verschlüsselung (als Anwendungsbeispiel der Zahlentheorie)
- [ ] Diophantische Gleichungen (lineare, quadratische)

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
