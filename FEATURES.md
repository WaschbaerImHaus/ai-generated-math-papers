# FEATURES.md - specialist-maths

## Implementierte Features

### Algebra (src/algebra_core.py, algebra_numbertheory.py, algebra_diophantine.py) - Build 2-11
- [x] Polynomklasse mit Horner-Schema-Auswertung
- [x] Polynom-Addition, Subtraktion, Multiplikation
- [x] Polynom-Ableitung (Potenzregel)
- [x] Linearer Gleichungslöser (mit Sonderfällen: keine/unendlich viele Lösungen)
- [x] Quadratischer Gleichungslöser (inkl. komplexe Wurzeln)
- [x] Euklidischer Algorithmus (ggT, kgV)
- [x] Erweiterter euklidischer Algorithmus (Bezout-Koeffizienten)
- [x] Modulares Inverses
- [x] Primzahltest (optimiert: 6k±1-Methode)
- [x] Primfaktorzerlegung
- [x] Eulersche Phi-Funktion
- [x] RSA-Kryptographie (keygen, encrypt, decrypt)
- [x] Diophantische Gleichungen: linear, Pell, Pythagoreer, Zwei-Quadrate, Markov
- [x] Quadratische Reste: is_quadratic_residue, quadratic_residues, quadratic_reciprocity
- [x] Wurzelberechnung mod p: Tonelli-Shanks, Cipolla

### Analysis (src/analysis.py) - Build 2-9
- [x] Numerische Ableitung 1./2. Ordnung (zentraler Differenzenquotient)
- [x] Numerische Integration (Simpson-Regel, O(h⁴))
- [x] Newton-Raphson, Bisektionsverfahren
- [x] Taylor-Reihen-Entwicklung (symbolisch via SymPy)
- [x] Symbolische Grenzwerte (symbolic_limit, lhopital_applicable, limit_comparison)
- [x] Partialbruchzerlegung (symbolisch und numerisch)
- [x] Uneigentliche Integrale (numerisch und symbolisch)
- [x] Cauchy-Hauptwert

### Lineare Algebra (src/vectors.py, matrix_ops.py, matrix_decomp.py) - Build 2-8
- [x] Vektor: Skalarprodukt, Kreuzprodukt, Norm, Normierung
- [x] Matrix: Determinante, Inverse, LGS lösen
- [x] Matrix: Eigenwerte und Eigenvektoren (QR-Iteration + SVD-Kern)
- [x] Gram-Schmidt-Orthogonalisierung
- [x] LU-Zerlegung (Doolittle + Teilpivot), QR-Zerlegung (Householder)
- [x] SVD, Rang, Konditionszahl
- [x] Givens-Rotation, Givens-QR, Kleinste-Quadrate

### Statistik & Wahrscheinlichkeit (src/statistics_math.py) - Build 2
- [x] Deskriptive Statistik, Verteilungsfunktionen (Normal, Binomial, Poisson)
- [x] Hypothesentests, Bayes-Theorem, Monte-Carlo-Simulation

### Differentialgleichungen (src/ode.py) - Build 2
- [x] Euler, Runge-Kutta 4, RK45 (adaptiv), Zustandsraum, Laplace/ILT

### Beweistheorie (src/proof_theory.py) - Build 4
- [x] Collatz, Goldbach, Zwillingsprimzahlen, Riemann-Zeta, Miller-Rabin
- [x] Legendre/Jacobi-Symbol, CRT, Hardy-Littlewood-Kreismethode

### Komplexe Analysis (src/complex_analysis.py) - Build 5
- [x] Gamma (Lanczos), vollständige Riemann-Zeta ζ(s), ξ-Funktion
- [x] Riemann-Siegel-Z, Nullstellensuche, N(T)-Formel, Cauchy-Integral
- [x] **NEU Build 11**: Arbitrary Precision via mpmath: riemann_zeta_mpmath, gamma_mpmath
- [x] **NEU Build 11**: riemann_siegel_z_mpmath, find_zeta_zeros_mpmath, verify_riemann_hypothesis_range_mpmath, li_function_mpmath

### Analytische Zahlentheorie (src/analytic_number_theory.py) - Build 5
- [x] π(x), Li(x), von-Mangoldt-Funktion Λ(n), θ(x)/ψ(x)
- [x] Dirichlet-L-Reihen, Chen-Primzahlsatz, Explizite Formel

### Fourier-Analysis (src/fourier.py) - Build 6
- [x] DFT, FFT (Cooley-Tukey), IFFT, Fourier-Reihen
- [x] Fensterfunktionen (Hann, Hamming, Blackman), STFT

### Numerische Methoden (src/numerical_methods.py) - Build 6
- [x] Lagrange/Newton/CubicSpline-Interpolation
- [x] Gradientenverfahren, Goldener Schnitt, Simplex-Verfahren, BFGS

### Modulformen (src/modular_forms.py) - Build 9-11
- [x] ModularGroup SL(2,Z), Eisenstein-Reihen, E4/E6, Delta, j-Invariante
- [x] Ramanujan-Tau, Hecke-Operatoren, Shimura-Taniyama, Cusp-Formen, Theta-Reihen
- [x] **NEU Build 11**: Hecke-Algebra Vertiefung: hecke_algebra_structure, hecke_eigenform, petersson_inner_product_estimate
- [x] **NEU Build 11**: L-Funktionen elliptischer Kurven: l_function_elliptic_curve, trace_of_frobenius, birch_swinnerton_dyer_bsd_estimate
- [x] **NEU Build 11**: L-Funktion Modulformen: l_function_modular_form, functional_equation_l_function

### p-adische Zahlen (src/p_adic.py) - Build 9-11
- [x] p_adic_valuation, p_adic_norm, PAdicNumber, hensel_lift, p_adic_exp/log, ostrowski_theorem_demo
- [x] **NEU Build 11**: bernoulli_number, generalized_bernoulli_numbers
- [x] **NEU Build 11**: kubota_leopoldt_l_function, p_adic_zeta_function, kummer_congruence_check, iwasawa_mu_lambda_invariants

### Visualisierung (src/visualization.py) - Build 9-11
- [x] 2D/3D Plotter, Vektorfeld, Phasenporträt
- [x] Fraktale: Mandelbrot/Julia/Sierpinski/Newton
- [x] **NEU Build 11**: animate_ode_trajectory, animate_phase_portrait, animate_wave_equation (GIF-Export)
- [x] **NEU Build 11**: save_figure_svg, save_figure_pdf, plot_and_export, export_all_formats
- [x] **NEU Build 11**: mandelbrot_smooth (Smooth Iteration Count, vollständig NumPy-vektorisiert)

### REPL-Modus (src/repl.py) - Build 10
- [x] Interaktive Kommandozeile, LaTeX-Export

## Optimierungen - Build 11
- [x] config.py: Globale Konstanten (kein Magic Numbers)
- [x] __init__.py: src/ als Python-Package
- [x] sympify()-Härtung: Alle unsicheren Aufrufe in analysis.py gesichert
- [x] Modulaufteilung linear_algebra.py → vectors.py, matrix_ops.py, matrix_decomp.py
- [x] Modulaufteilung algebra.py → algebra_core.py, algebra_numbertheory.py, algebra_diophantine.py
- [x] Mandelbrot/Julia: Vollständige NumPy-Vektorisierung + Smooth Coloring

## Offene Aufgaben
_(keine bekannten offenen Features mehr)_
