# FEATURES.md - specialist-maths

## Implementierte Features

### Statistik (src/statistics_math.py) - Build 2
- [x] Deskriptive Statistik: mean, median, mode, variance, std_dev
- [x] Lage- und Streuungsmaße: quartiles, iqr, skewness, kurtosis
- [x] Normalverteilung: PDF, CDF, PPF (Quantilfunktion)
- [x] Binomialverteilung: PMF, CDF
- [x] Poisson-Verteilung: PMF, CDF
- [x] Hypothesentests: Ein-/Zwei-Stichproben-t-Test, Chi-Quadrat-Test
- [x] Bayes-Theorem
- [x] Monte-Carlo-Simulation (pi-Schätzung)

### Differentialgleichungen (src/ode.py) - Build 2
- [x] Euler-Verfahren (explizit, 1. Ordnung)
- [x] Runge-Kutta 4. Ordnung (klassisch)
- [x] Runge-Kutta-Fehlberg RK45 (adaptiv)
- [x] Lineare ODE mit konstanten Koeffizienten (Zustandsraum + RK4)
- [x] Numerische Laplace-Transformation (Simpson-Regel)
- [x] Numerische inverse Laplace-Transformation (Stehfest-Algorithmus)

### Algebra (src/algebra.py) - Build 1
- [x] Polynomklasse mit Horner-Schema-Auswertung
- [x] Polynom-Addition, Subtraktion, Multiplikation
- [x] Polynom-Ableitung (Potenzregel)
- [x] Polynom String-Darstellung
- [x] Linearer Gleichungslöser (mit Sonderfällen)
- [x] Quadratischer Gleichungslöser (inkl. komplexe Wurzeln)
- [x] Euklidischer Algorithmus (ggT, kgV)
- [x] Erweiterter euklidischer Algorithmus (Bezout)
- [x] Modulares Inverses
- [x] Primzahltest (optimiert: 6k±1-Methode)
- [x] Primfaktorzerlegung
- [x] Eulersche Phi-Funktion

### Analysis (src/analysis.py) - Build 1
- [x] Numerische Ableitung 1. Ordnung (zentraler Differenzenquotient)
- [x] Numerische Ableitung 2. Ordnung (mit optimaler Schrittweite)
- [x] Numerische Integration (Simpson-Regel, O(h^4) Genauigkeit)
- [x] Newton-Raphson-Verfahren (Nullstellensuche)
- [x] Bisektionsverfahren (Nullstellensuche)
- [x] Taylor-Reihen-Entwicklung (symbolisch via SymPy)

### Lineare Algebra (src/linear_algebra.py) - Build 1
- [x] Vektor: Skalarprodukt, Kreuzprodukt, Norm, Normierung
- [x] Vektor: Addition, Subtraktion, Skalarmultiplikation
- [x] Matrix: Erstellen, Identity, Zeros
- [x] Matrix: Multiplikation (matmul)
- [x] Matrix: Transponierung, Spur
- [x] Matrix: Determinante (Gauss mit Pivotsuche)
- [x] Matrix: Inverse (Gauss-Jordan)
- [x] Matrix: LGS lösen (Gauss mit Rücksubstitution)
- [x] Matrix: Eigenwerte (analytisch 2x2, QR-Iteration n×n)
- [x] Gram-Schmidt-Orthogonalisierung (mit ONB-Option)

## Geplante Features

### Algebra
- [ ] Polynomoperationen (Addition, Multiplikation, Division mit Rest)
- [ ] Gleichungslöser (linear, quadratisch, kubisch, quartic)
- [ ] Polynom-Faktorisierung
- [ ] Gruppentheorie-Grundlagen (Ordnung, Untergruppen, Homomorphismen)
- [ ] Matrizenrechnung über endlichen Körpern

### Analysis
- [ ] Grenzwertberechnung (symbolisch + numerisch)
- [ ] Ableitungsrechner (Kettenregel, Produktregel, etc.)
- [ ] Integralrechner (unbestimmt, bestimmt, uneigentlich)
- [ ] Taylor- und Maclaurin-Reihen
- [ ] Fourier-Reihen und -Transformation
- [ ] Partialbruchzerlegung

### Lineare Algebra
- [ ] Gauss-Jordan-Elimination
- [ ] LU-, QR-, SVD-Zerlegung
- [ ] Eigenwert- und Eigenvektorberechnung
- [ ] Determinantenberechnung (verschiedene Methoden)
- [ ] Vektorraum-Operationen (Orthogonalisierung, Gram-Schmidt)

### Zahlentheorie
- [ ] Primzahl-Sieb (Eratosthenes, Atkin)
- [ ] Erweiterter Euklidischer Algorithmus
- [ ] Modulare Arithmetik (Chinesischer Restsatz)
- [ ] RSA-Verschlüsselung (als Anwendungsbeispiel)
- [ ] Diophantische Gleichungen

### Statistik & Wahrscheinlichkeit
- [ ] Deskriptive Statistik
- [ ] Wahrscheinlichkeitsverteilungen
- [ ] Hypothesentests (t-Test, Chi-Quadrat, etc.)
- [ ] Bayes-Theorem-Rechner
- [ ] Monte-Carlo-Simulation

### Differentialgleichungen
- [ ] ODE-Löser (Runge-Kutta, Adams-Bashforth)
- [ ] Lineare ODE mit konstanten Koeffizienten
- [ ] Laplace-Transformation
- [ ] Phasenportraits und Stabilität

### Numerische Methoden
- [ ] Bisektionsverfahren
- [ ] Newton-Raphson-Verfahren
- [ ] Numerische Integration (Trapez, Simpson, Gauss)
- [ ] Interpolation (Lagrange, Newton, Splines)
- [ ] Optimierungsverfahren (Gradient Descent)

### Visualisierung
- [ ] Funktionsplotter (2D und 3D)
- [ ] Vektorfelddarstellung
- [ ] Phasenraum-Visualisierung
- [ ] Fraktal-Generator (Mandelbrot, Julia)
