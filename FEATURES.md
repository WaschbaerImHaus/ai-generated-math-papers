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

### Topologie & Geometrie (src/topology.py) - Build 12
- [x] Standardmetriken: Euklidisch, Manhattan, Chebyshev, Diskret, Lp
- [x] MetricSpace-Klasse mit Axiomprüfung (alle 4 Axiome)
- [x] Offene Kugeln, offene Mengen (diskret)
- [x] Cauchy-Folgen Erkennung
- [x] Epsilon-Zusammenhang (Union-Find)
- [x] Hausdorff-Abstand H(A,B)
- [x] Durchmesser einer Punktmenge
- [x] ParametricCurve: Bogenlänge, Krümmung, Umlaufzahl, Abgeschlossenheit
- [x] Vorgefertigte Kurven: Kreis, Lissajous, Helix
- [x] Euler-Charakteristik, Geschlecht
- [x] Betti-Zahlen für Graphen (β₀, β₁)
- [x] Box-Counting-Dimension
- [x] Hausdorff-Dimension: Cantor-Menge, Sierpinski-Dreieck

### Graphentheorie & Kombinatorik (src/graph_theory.py) - Build 12
- [x] Graph-Klasse (gerichtet/ungerichtet, gewichtet, Adjazenzliste)
- [x] Grundoperationen: add/remove vertex/edge, degree, Adjazenzmatrix
- [x] Grapheigenschaften: Zusammenhang, Zyklen, Baumtest, Bipartitheit
- [x] Komplementgraph, induzierter Teilgraph
- [x] BFS (Breitensuche) mit Abständen und Vorgänger
- [x] DFS (Tiefensuche) mit Entdeckungs-/Abschlusszeiten
- [x] Dijkstra-Algorithmus (MinHeap, O((V+E)log V))
- [x] Bellman-Ford (negative Gewichte, Zykluserkennung)
- [x] Floyd-Warshall (alle kürzesten Wege, O(V³))
- [x] Kruskal-MST (Union-Find, O(E log E))
- [x] Prim-MST (MinHeap, O((V+E) log V))
- [x] Topologische Sortierung (Kahn-Algorithmus)
- [x] Greedy-Graphfärbung, chromatische Zahl (Obergrenze)
- [x] Euler-Kreis und Euler-Pfad Erkennung
- [x] Hamiltonpfad-Suche (Backtracking)
- [x] Kliquenzahl ω(G) (Bron-Kerbosch)
- [x] Unabhängigkeitszahl α(G) = ω(Komplement)
- [x] Graphdichte
- [x] Graphkonstruktoren: K_n, C_n, P_n, K_{m,n}, Petersen, Gitter
- [x] Binomialkoeffizient, Stirling-Zahlen (2. Art), Bell-Zahlen
- [x] Catalan-Zahlen, Derangements, Ganzzahlpartitionen, Multinomialkoeffizient

### Elliptische Kurven (src/elliptic_curves.py) - Build 16
- [x] ECPoint: Punkt auf elliptischer Kurve (Weierstraß-Normalform)
- [x] ECPoint.infinity(): Punkt im Unendlichen (Neutralelement)
- [x] Chord-and-Tangent-Additionsgesetz (alle 5 Fälle)
- [x] Double-and-Add Skalarmultiplikation O(log n)
- [x] EllipticCurve: y² = x³ + ax + b über ℝ
- [x] Diskriminante Δ = -16(4a³ + 27b²)
- [x] j-Invariante: j = -1728·(4a)³/Δ
- [x] Punktberechnung über E(F_p): Brute-Force + Tonelli-Shanks
- [x] Gruppenordnung über F_p, Frobenius-Spur a_p
- [x] Hasse-Schranke: |#E(F_p) - (p+1)| ≤ 2√p
- [x] EllipticCurveModP: Arithmetik über endlichem Körper F_p
- [x] Baby-Step Giant-Step für diskreten Logarithmus O(√n)
- [x] Generator-Test für E(F_p)
- [x] ECCKeyExchange: ECDH-Schlüsselaustausch Protokoll
- [x] lenstra_ecm_factorization: ECM via projektive Koordinaten
- [x] nagell_lutz_theorem: Ganzzahlige Torsionspunkte
- [x] mordell_weil_rank_estimate: Schaetzung des Mordell-Weil-Rangs
- [x] l_function_rank_order: BSD L-Funktions Nullstellenordnung
- [x] is_supersingular / endomorphism_ring_type
- [x] secp256k1(): Bitcoin-Kurve
- [x] curve25519(): Bernstein-Kurve (Signal, WireGuard, TLS 1.3)
- [x] example_bsd_curve(): y² = x³ - x (Rang 0)
- [x] congruent_number_curve(n): Verbindung Geometrie-Zahlentheorie
- [x] Dokumentation elliptic_curves.md mit KaTeX-Formeln, BSD, Shimura-Taniyama-Wiles

### Logging-System (src/math_logger.py) - Build 14
- [x] MathLogger-Klasse mit konfigurierbaren Leveln (DEBUG/INFO/WARNING/ERROR)
- [x] Ausgabe auf Konsole und/oder Datei (logs/-Verzeichnis)
- [x] step(), result(), convergence_warning(), matrix_step(), timing(), section()
- [x] Globale Hilfsfunktionen: get_logger(), enable_debug_logging(), disable_logging()
- [x] newton_raphson() in analysis.py um verbose=True erweitert

### Schritt-für-Schritt-Modus (src/step_by_step.py) - Build 14
- [x] newton_raphson_steps(), bisection_steps(), gauss_elimination_steps()
- [x] euclidean_algorithm_steps(), rsa_steps(), prime_factorization_steps(), lu_decomposition_steps()
- [x] format_steps_text() und format_steps_html() mit KaTeX-Formeln

### Web-Interface: Schritt-für-Schritt (Build 14)
- [x] GET /steps: Neue interaktive Seite mit Fade-in-Animationen
- [x] POST /api/steps/{newton,bisection,gauss,gcd,rsa,prime_factorization}
- [x] Navigationseintrag "Schritt-für-Schritt" in Sidebar

### L-Funktionen-Infrastruktur (src/l_functions.py) - Build 22
- [x] dirichlet_characters(q): Alle φ(q) Dirichlet-Charaktere mod q (Primitivwurzel-Methode)
- [x] chi_value(n, q, idx): Einzelner Charakterwert χ_idx(n)
- [x] is_primitive_character(): Gauß-Summen-Kriterium |τ(χ)|² = q
- [x] dirichlet_l_function(): L(s,χ) = Σ χ(n)/n^s mit Euler-Knopp für Re(s)≤1
- [x] dirichlet_l_function_zeros(): Nullstellensuche auf Re(s)=1/2 (GRH-Streifen)
- [x] gauss_sum(q, idx): τ(χ) = Σ_{n=1}^q χ(n)·e^{2πin/q}
- [x] l_function_conductor(): Führer f(χ) als kleinster induzierender Teiler
- [x] l_function_special_values(): L(1,χ), L(0,χ), gerade/ungerade Klassifikation
- [x] functional_equation_check(): Λ(s,χ) = ε(χ)·Λ(1-s,χ̄) mit Root Number ε
- [x] hecke_l_function(): L(s,f) = Σ a_n/n^s für Modulformen (Ramanujan-Delta)
- [x] completed_l_function(): Λ(s,f) = (2π)^{-s}Γ(s)L(s,f) (Funktionalgleichung)
- [x] elliptic_l_function_approx(): L(E,s) via Euler-Produkt, Frobenius-Spur a_p
- [x] bsd_rank_from_zeros(): BSD-Rang-Schätzung via L(E,1)≈0 Kriterium
- [x] 72 Tests in tests/test_l_functions_dirichlet.py, alle grün
- [x] Dokumentiert: Leibniz-Reihe L(1,χ₋₄)=π/4, |τ|²=q, BSD-Vermutung

### Moduln über Ringen (src/modules_algebra.py) - Build 27
- [x] smith_normal_form(): Smith-Normalform A = U·D·V mit Teilbarkeitskette d₁|d₂|...|dᵣ
- [x] module_from_matrix(): ℤ-Modul M = ℤⁿ/Im(A) aus Relationsmatrix
- [x] structure_theorem_abelian_groups(): Klassifikation aller abelschen Gruppen der Ordnung n
- [x] tensor_product_modules(): ℤ_m ⊗ ℤ_n ≅ ℤ_{gcd(m,n)}
- [x] hom_module(): |Hom(ℤ_m, ℤ_n)| = gcd(m,n)
- [x] exact_sequence_check(): Exaktheit 0 → A → B → C → 0 prüfen
- [x] free_resolution(): Freie Auflösung via Smith-Normalform
- [x] Module-Klasse: rank(), torsion_submodule(), free_part(), is_free(), is_finitely_generated()
- [x] 42 Tests in tests/test_modules_algebra.py, alle grün

### Kommutative Algebra (src/commutative_algebra.py) - Build 27
- [x] localization(): S^{-1}R Lokalisierung an multiplikativer Menge
- [x] prime_spectrum(): Spec(ℤ/nℤ) = {(p) : p Primteiler von n}
- [x] nakayama_lemma_check(): Nakayama M = I·M ⟹ M = 0 Verifikation
- [x] integral_closure(): Ganzer Abschluss von ℤ in ℚ(√d), Diskriminante, PID-Test
- [x] class_group_estimate(): Klassenzahl h(d), Minkowski-Schranke, PID-Nachweis
- [x] noether_normalization(): Krull-Dimension, algebraisch unabhängige Elemente
- [x] hilbert_basis_theorem_verify(): Endliche Erzeugbarkeit im Polynomring über ℤ_p
- [x] 32 Tests in tests/test_commutative_algebra.py, alle grün

## Offene Aufgaben
_(keine bekannten offenen Features mehr)_

### Algebra Vertiefung (Build 33-34)
- [x] Darstellungstheorie (src/representation_theory.py): Charaktertafeln, Schur-Orthogonalität, irreduzible Darstellungen
- [x] Verbandstheorie (src/lattice_theory.py): Halbverbände, Verbände, distributive/modulare Verbände, Boolesche Algebren
- [x] Multilineare Algebra (src/multilinear_algebra.py): Tensorprodukte, äußere Algebra, symmetrische Algebra, Tensorkomponenten
- [x] Homologische Algebra (src/homological_algebra.py): Kettenkomplexe, Homologiegruppen, exakte Sequenzen, Ext/Tor
- [x] Algebraische Strukturen (src/algebraic_structures.py): Magma/Halbgruppe/Monoid/Gruppe-Hierarchie, freie Strukturen
- [x] Invariantentheorie (src/invariant_theory.py): Reynolds-Operator, Molien-Reihe, Hilbert-Basis, Newton-Identitäten
- [x] Universelle Algebra (src/universal_algebra.py): Signaturen, Varietäten, freie Algebren, Birkhoff-Theorem

### Logische Grundlagen (Build 33-34)
- [x] Modelltheorie (src/model_theory.py): Strukturen, Interpretationen, Elementar-Äquivalenz, Compactness-Theorem
- [x] Formale Beweistheorie (src/proof_theory_formal.py): Sequenzenkalkül LK, Natürliches Schließen, Hilbert-Kalkül, Schnittelimination
- [x] Rekursionstheorie (src/recursion_theory.py): Turing-Maschinen, μ-rekursive Funktionen, Halteproblem, Arithmetische Hierarchie

### Analysis Vertiefung (Build 34)
- [x] Maßtheorie (src/measure_theory.py): σ-Algebren, Lebesgue-Maß, Lebesgue-Integral, Konvergenzsätze, Cantor-Funktion
- [x] Spezielle Funktionen (src/special_functions.py): Bessel, Legendre, Airy, Hypergeometrisch, Orthogonalpolynome
- [x] Funktionalanalysis (src/functional_analysis.py): Banach/Hilbert-Räume, lineare Operatoren, Spektraltheorie, Fixpunktsätze
- [x] Partielle DGL (src/pde.py): Wärmegleichung, Wellengleichung, Laplace/Poisson, Schrödinger, FEM, Charakteristiken
- [x] Operatoralgebren (src/operator_algebras.py): C*-Algebren, Von-Neumann-Algebren, Gelfand-Darstellung, K-Theorie
