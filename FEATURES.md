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

### Builds 57-62 (2026-03-11)

#### Sicherheit (Build 57)
- [x] repl.py: Sicherheitsaudit-Kommentar erweitert (vollständige Dokumentation der Sicherheitseigenschaften)
- [x] analysis.py: _safe_parse() loggt WARNING bei Fallback auf sympify() (verifiziert)

#### Performance (Build 58)
- [x] fourier.py: polynomial_multiply_fft() – FFT-basierte Polynommultiplikation O(n log n) via numpy.fft
- [x] fourier.py: polynomial_multiply_naive() – O(n²) Referenzimplementierung für Vergleiche
- [x] algebra_numbertheory.py: euler_phi() nutzt prime_factorization() mit @lru_cache (verifiziert)
- [x] tensor_geometry.py: Christoffel-Symbol-Cache bereits implementiert (verifiziert)

#### Architektur (Build 59)
- [x] category_theory.py: Adjunction-Klasse (F ⊣ G, unit, counit, Dreieck-Identitäten)
- [x] category_theory.py: DiscreteCategory-Klasse (nur Identitätsmorphismen)
- [x] category_theory.py: FinSetCategory-Klasse (endliche Mengen, enumerate_functions)
- [x] category_theory.py: discrete_category() Fabrikfunktion
- [x] category_theory.py: yoneda_lemma_demo() Demonstration
- [x] plugin_registry.py: DEFAULT_REGISTRY mit 25 weiteren Modulen erweitert

#### Visualisierung (Build 60)
- [x] visualization.py: plot_adaptive_grid() – Adaptiver Gitterplot via Bisektions-Strategie
- [x] visualization.py: create_interactive_plot() – Interaktiver Plot mit matplotlib.widgets.Slider

#### Mathematik (Build 61)
- [x] arbitrary_precision.py: zeta_zeros_mpmath(n, prec) – Erste n Riemann-Nullstellen
- [x] arbitrary_precision.py: verify_riemann_hypothesis_mpmath(N, prec) – Numerische RH-Verifikation
- [x] arbitrary_precision.py: pi_mpmath(prec) – π auf prec Stellen
- [x] notebook_utils.py: display_matrix_html(M) – IPython-kompatible Matrixanzeige
- [x] notebook_utils.py: display_polynomial_latex(coeffs) – IPython-kompatible Polynomdarstellung

#### Neue Papers – Satz von Wilson (Build 62)
- [x] papers/batch2/paper8_wilson_theorem.tex (EN): 5 Sätze (Wilson, Fakultätsrest, Primzahlpotenzen, Quotient, allgemein)
- [x] papers/batch2/paper8_wilson_satz_de.tex (DE): Vollständige deutsche Übersetzung

#### Neue Papers – Wilson für Primzahlpotenzen (Build 63)
- [x] papers/batch2/paper9_wilson_prime_powers.tex (EN): Vollst. Beweis via Paarung + Struktur
- [x] papers/batch2/paper9_wilson_primzahlpotenzen_de.tex (DE): Deutsche Version

#### Neue Papers – Wilson-Quotient und Wilson-Primzahlen (Build 64)
- [x] papers/batch2/paper10_wilson_quotient.tex (EN): W_p, Wilson-Primzahlen, Bernoulli, Heuristik
- [x] papers/batch2/paper10_wilson_quotient_de.tex (DE): Deutsche Version

#### Neue Papers – Allgemeiner Wilson für abelsche Gruppen (Build 65)
- [x] papers/batch2/paper11_wilson_abelian.tex (EN): Gruppentheoretische Verallgemeinerung
- [x] papers/batch2/paper11_wilson_abelsch_de.tex (DE): Deutsche Version

#### Neue Papers – Anwendungen des Wilson-Satzes (Build 66)
- [x] papers/batch2/paper12_wilson_applications.tex (EN): 4 Anwendungen inkl. QR, Fermat, Wolstenholme
- [x] papers/batch2/paper12_wilson_anwendungen_de.tex (DE): Deutsche Version
- [x] tests/test_wilson_theorems.py: 78 Tests, alle grün (Wilson-Themen)

### Offene Vermutungen — Gruppe A (Build 125)

#### Giuga 4-Prim-Fall (src/giuga_4prim.py)
- [x] `berechne_giuga_bedingung(n)` — schwache + starke Giuga-Bedingung, erkennt Giuga-Zahl vs. Pseudoprime
- [x] `numerische_suche_4prim(grenze)` — Suche über 4-Prim-Produkte bis Grenze; bis 10^6: kein Pseudoprime
- [x] `Giuga4PrimBeweis.fall1_p_gleich_2_analyse()` — Schrankenanalyse für p=2 (n=2·q·r·s)
- [x] `Giuga4PrimBeweis.fall2_p_gleich_3_analyse()` — Schrankenanalyse für p=3
- [x] `Giuga4PrimBeweis.fall3_p_ungerade_analyse()` — allgemeiner ungerade Fall
- [x] `Giuga4PrimBeweis.schranken_analyse()` — untere Schranken für alle Primteiler
- [x] `Giuga4PrimBeweis.korollar_kein_4prim_pseudoprime()` — Status: numerisch bis 10^6, kein elem. Beweis
- [x] 50 pytest-Tests (tests/test_giuga_4prim.py)
- **Befund**: 3-Prim-Widerspruch nicht auf k=4 übertragbar; 4. Faktor absorbiert Widerspruch; Literatur schließt k=4 implizit aus (Borwein 1996: min. 13635 Primteiler)

#### Erdős-Straus Erweiterung (src/erdos_straus_ext.py)
- [x] `ErdosStrausExt.loese_unit_fraction(n)` — 3-stufiger Algorithmus (Formel → Teiler → Brute-Force)
- [x] `beweis_klasse_mod4_3(p)` — **Formaler Beweis**: 4/p = 1/((p+1)/4) + 1/((p+1)/4·p) für p≡3(mod 4) ✓
- [x] `beweis_klasse_mod4_1(p)` — Konstruktive Lösung für p≡1(mod 4)
- [x] `beweis_klasse_mod3`, `beweis_klasse_mod_12` — kombinierte Restklassenbeweise
- [x] `vollstaendige_analyse`, `suche_ausnahmen`, `parametrische_loesungen`
- [x] 89 pytest-Tests (tests/test_erdos_brocard.py)
- **Befund**: Für p≡3(mod 4) vollständig konstruktiv bewiesen (2-Term-Zerlegung, exakt algebraisch)

#### Brocard-Ramanujan Numerische Erweiterung (src/brocard_extension.py)
- [x] `BrocardExtension.numerische_suche(n_max)` — exakte Quadratwurzelprüfung via `math.isqrt`; n≤1000 in 11ms; nur n=4,5,7 gefunden
- [x] `modular_ausschluss`, `finde_ausschlusskandidaten`, `analysiere_restklassen`
- [x] `schranken_analyse` — Stirling-Schranken für m(n)
- [x] `p_adische_analyse` — Legendre-Formel, v_p(n!+1)-Analyse

#### Kurepa-Vermutung Erweiterung (src/kurepa_ext.py)
- [x] `KurepaExt.numerische_verifikation(p_max)` — bis p=50000 (5133 Primzahlen) verifiziert, ~7s
- [x] `p_adische_bewertung_analyse`, `restklassen_struktur`, `wilsons_verbindung`
- [x] Rekursion !p ≡ !(p-1) - 1 (mod p) formal dokumentiert
- [x] 93 pytest-Tests (tests/test_kurepa_tau_bruns.py)

#### Lehmer τ(n) ≠ 0 Analyse (src/lehmer_tau.py)
- [x] `LehmerTauAnalyse.berechne_tau(n)` — exakte Koeffizientenextraktion aus Δ(q) = q·∏(1-qⁿ)²⁴
- [x] `verifiziere_lehmer(n_max)` — τ(n) ≠ 0 für n ≤ 1000 bestätigt
- [x] `kongruenz_analyse` — Kongruenzen mod 2,3,5,7,11,13,23,691
- [x] `deligne_schranke` — |τ(p)| ≤ 2p^{11/2} geprüft (Deligne 1974)
- [x] `ausschluss_durch_kongruenzen` — systematische Ausschlussanalyse

#### Bruns Konstante (src/bruns_constant.py)
- [x] `BrunsKonstante.bruns_summe(grenze)` — partielle Summe B₂(x); B₂(10^6) ≈ 1.7108
- [x] `hochpraezise_berechnung` — mpmath-basierte Hochpräzision
- [x] `hardylittlewood_vorhersage` — C₂ = ∏_{p>2} p(p-2)/(p-1)² ≈ 0.6601618 berechnet
- [x] Asymptotische Extrapolation liefert ≈ 1.9022 (bekannt: 1.9021605831...)

#### Schur-Zahlen S(6) (src/schur_numbers.py)
- [x] `SchurZahlen.partition_s5()` — explizite Partition {1,...,160} in 5 sum-freie Klassen
- [x] `untere_schranke_s6()` — konstruktive Partition bis N=536 verifiziert
- [x] `schranken_analyse()` — 536 ≤ S(6) ≤ 1836 dokumentiert
- [x] Greedy + Monte-Carlo-Suche nach verbesserten unteren Schranken
- [x] SAT-Encoding-Beschreibung implementiert
- [x] 115 pytest-Tests (tests/test_schur_debruijn.py)

#### de Bruijn-Newman-Konstante (src/debruijn_newman.py)
- [x] `DeBruijnNewman.H_t(x, t)` — numerische Berechnung von H_t via mpmath
- [x] `suche_nullstellen_H_t` — Bisektion-basierte Nullstellensuche
- [x] `newman_vermutung_status()` — Λ ≥ 0 bewiesen (Rodgers-Tao 2018) dokumentiert
- [x] Schranke 0 ≤ Λ ≤ 0.2 (Platt-Trudgian 2021) implementiert
- [x] Verbindung H_0 ↔ ½ξ(½+ix) formal dokumentiert
- [x] Historische Schranken-Tabelle: Λ-Bounds von 1990 bis 2021

#### Gesamt neue Tests (Build 125): 347 Tests (50 + 89 + 93 + 115)

#### Erdős-Moser-Vermutung (src/erdos_moser.py) - Build 122
- [x] `ErdosMoser.is_solution(k, m)` — Prüft S_k(m-1) = m^k
- [x] `find_solutions(k_max, m_max)` — Exhaustive Suche (Ergebnis: nur (1,3))
- [x] `verify_k1_solution()` — Analytischer Beweis für k=1, m=3
- [x] `modular_residues(k, q)` — Modulare Ausschluss-Analyse
- [x] `is_excluded_by_modulus(k, m, q)` — Einzelprüfung per Modulus
- [x] `bernoulli_number(n)` — Exakte Bernoulli-Zahlen B_n via sympy
- [x] `faulhaber_formula(N, k)` — S_k(N) via Faulhaberscher Formel
- [x] `mosers_lower_bound()` — Moser 1953: m > 10^(10^6), Gallot-Moree-Zudilin 2011
- [x] CONJECTURE-Status korrekt dokumentiert

#### Artin-Vermutung (src/artin_primitive_roots.py) - Build 122
- [x] `ArtinPrimitiveRoots.is_primitive_root(a, p)` — ord_p(a) = p-1 Prüfung
- [x] `primitive_root_check_detailed(a, p)` — Detailanalyse mit Ordnung
- [x] `empirical_density(a, limit)` — Empirische Dichte vs Artin-Konstante
- [x] `artin_constant(num_primes)` — C_A = ∏(1-1/(p(p-1))) ≈ 0.3739558...
- [x] `hooley_density(a)` — GRH-basierte Dichte (Hooley 1967)
- [x] `heath_brown_result()` — Unbedingtes Theorem (Heath-Brown 1986)
- [x] Basen a=2,3,5,7 vollständig analysiert
- [x] CONJECTURE-Status (unter GRH von Hooley bewiesen)

#### Erdős-Selfridge-Vermutung (src/erdos_selfridge.py) - Build 122
- [x] `ErdosSelfridge.is_binomial_prime_power(n, k)` — C(n,k) = p^a Prüfung
- [x] `verify_up_to(n_max)` — Verifikation bis n=200 (keine Gegenbeispiele)
- [x] `sylvester_theorem_check(n, k)` — Sylvester 1892: Primteiler > k wenn n≥2k
- [x] `verify_sylvester(n_max)` — Massenverifikation des Theorems
- [x] `analyze_k1(n_max)`, `analyze_k2(n_max)` — Spezialfälle
- [x] `granville_special_cases()` — Granvilles Teilbeweise dokumentiert
- [x] CONJECTURE-Status korrekt (k≥2, n≥k+2)

---

### Build 127: Mahler-Maß, Vollkommene Zahlen, Waring-Goldbach

#### Mahler-Maß (src/mahler_measure.py)
- [x] `MahlerMeasure.compute_product_formula()` — M(P) = |aₙ|·∏_{|αᵢ|>1}|αᵢ| (Produktformel)
- [x] `MahlerMeasure.compute_jensen_integral(n)` — M(P) via exp(∫₀¹ log|P(e^{2πit})|dt)
- [x] `MahlerMeasure.logarithmic_mahler_measure()` — m(P) = log M(P)
- [x] `MahlerMeasure.is_kronecker()` — M(P)=1 ↔ Kronecker-Theorem
- [x] `MahlerMeasure.is_reciprocal()` — palindromisch/antipalindromisch
- [x] `MahlerMeasure.lehmer_polynomial()` — Lehmer-Dezim-Polynom, M≈1.17628
- [x] `MahlerMeasure.smyth_polynomial()` — x³-x-1, Smyth-Konstante 1.3247...
- [x] `MahlerMeasure.cyclotomic(n)` — Kreisteilungspolynome (M=1)
- [x] `SchurSiegelSmythTrace.trace_ratio()` — Spurverhältnis total positiver alg. Zahlen
- [x] Lehrers Problem korrekt als CONJECTURE markiert (offen seit 1933)

#### Vollkommene Zahlen (src/perfect_numbers.py)
- [x] `SigmaFunction.sigma(n)` — Teiler-Summenfunktion σ(n)
- [x] `SigmaFunction.multiplicativity_check()` — Multiplikativität verifiziert
- [x] `SigmaFunction.is_perfect/abundant/deficient()` — Klassifikation
- [x] `SigmaFunction.abundancy_index()` — σ(n)/n
- [x] `EvenPerfectNumbers.first_n_even_perfect(n)` — erste n gerade vollkommene Zahlen
- [x] `EvenPerfectNumbers.is_mersenne_prime(p)` — Mersenne-Primzahltest
- [x] `EvenPerfectNumbers.verify_euclid_euler(p)` — vollständige Theorem-Verifikation
- [x] `OddPerfectNumberBounds.touchard_congruence(n)` — n ≡ 1(mod 12) oder n ≡ 9(mod 36)
- [x] `OddPerfectNumberBounds.nielsen_bound_check(n)` — ≥9 verschiedene Primfaktoren
- [x] `OddPerfectNumberBounds.euler_form_check(n)` — n = p^{4a+1}·m² Prüfung
- [x] Ungerade vollkommene Zahlen korrekt als CONJECTURE markiert

#### Waring-Goldbach (src/waring_goldbach.py)
- [x] `WaringGoldbachK2.represent_as_prime_squares()` — Darstellung als ≤5 Primzahlquadrate
- [x] `WaringGoldbachK2.count_representations()` — Anzahl der Darstellungen
- [x] `WaringGoldbachK2.minimal_representation()` — minimale Termsanzahl
- [x] `WaringGoldbachGeneral.g_bound_hua()` — Hua (1938) Schranken G(k)
- [x] `WaringGoldbachGeneral.g_bound_vinogradov_wooley()` — Vinogradov-Wooley
- [x] `WaringGoldbachGeneral.goldbach_conjecture_check()` — Goldbach (CONJECTURE)
- [x] `WaringGoldbachGeneral.vinogradov_three_primes_check()` — Vinogradov/Helfgott
- [x] Goldbach-Vermutung korrekt als CONJECTURE markiert
- [x] 111 Tests in tests/test_mahler_perfect_waring.py

### Graceful Trees (src/graceful_trees.py) — Build 129
- [x] `GracefulTreeConjecture.is_graceful(tree, labeling)` — Labeling-Prüfung
- [x] `GracefulTreeConjecture.find_graceful_labeling(tree)` — Backtracking (n≤10)
- [x] `GracefulTreeConjecture.generate_all_trees(n)` — Nicht-isomorphe Bäume
- [x] `GracefulTreeConjecture.graceful_path(n)` — Explizit: Pfade P_n
- [x] `GracefulTreeConjecture.graceful_star(k)` — Explizit: Sterne K_{1,k}
- [x] `GracefulTreeConjecture.graceful_caterpillar(...)` — Caterpillar-Bäume
- [x] `GracefulTreeConjecture.graceful_wheel(n)` — Wheel-Graphen
- [x] `GracefulTreeConjecture.is_caterpillar(tree)` — Caterpillar-Prüfung
- [x] `GracefulTreeConjecture.count_graceful_trees(max_n)` — Statistik
- [x] `GracefulTreeConjecture.verify_known_classes()` — Verifikation bekannter Klassen
- [x] 106 Tests in tests/test_graceful_ramsey.py

### Ramsey-Zahlen (src/ramsey_numbers.py) — Build 129
- [x] `RamseyNumbers.get(s,t)` — Lookup bekannter exakter Werte
- [x] `RamseyNumbers.get_bounds(s,t)` — Schranken (lower, upper)
- [x] `RamseyNumbers.has_clique(graph, size)` — Clique-Prüfung
- [x] `RamseyNumbers.has_independent_set(graph, size)` — I-Menge-Prüfung
- [x] `RamseyNumbers.is_ramsey_witness(graph, s, t)` — Untergrenzen-Zeugnis
- [x] `RamseyNumbers.build_paley_graph(q)` — Paley-Graph P(q) für q≡1(mod 4)
- [x] `RamseyNumbers.build_r55_lower_bound_witness()` — P(29): Zeuge R(5,5)>29
- [x] `RamseyNumbers.verify_r33_equals_6()` — R(3,3)=6 vollständige Verifikation
- [x] `RamseyNumbers.verify_via_sat(s,t,n)` — SAT-basierte Verifikation
- [x] `RamseyBounds.upper_binomial(s,t)` — C(s+t-2,s-1)
- [x] `RamseyBounds.lower_erdos(k)` — floor(2^{k/2})
- [x] `RamseyBounds.upper_sah_2023(k)` — (4-ε)^k
- [x] `RamseyBounds.upper_spencer_1975(k)` — c·k·2^{k/2}/sqrt(ln k)
- [x] `RamseyBounds.diagonal_bounds(k_max)` — Übersichtstabelle
- [x] `RamseyBounds.multiplicity(s,t,n)` — Goodman-Formel

### Papers — Gruppe B Vermutungen (Batches 19–25) — Build 164–171
**52 neue LaTeX-Papers (26 Themen × EN + DE)**

#### Batch 19 — Additive Kombinatorik & Kombinatorik (Build 164)
- [x] paper72: Freiman-Struktursatz — Plünnecke-Ruzsa, Sanders, PFR über F₂ⁿ (Gowers-Green-Manners-Tao 2023)
- [x] paper73: Erdős-Ko-Rado-Verallgemeinerungen — Frankl t-intersecting, Ahlswede-Khachatrian, q-Analog
- [x] paper74: Lonely Runner Vermutung — n≤7 bewiesen, Diophantische Äquivalenz
- [x] paper75: Graceful Tree Vermutung (Ringel-Kotzig) — α-Labelings, Rosa-Satz

#### Batch 20 — Topologie & Dynamik (Build 165)
- [x] paper76: Mandelbrot-MLC-Vermutung — Douady-Hubbard, Yoccoz-Puzzle, Shishikura
- [x] paper77: Conway-Knoten Scheibigkeit — Piccirillo 2020 (glatt NICHT scheibig), topol. offen
- [x] paper78: Whitehead-Asphärizitätsvermutung — Lyndon-Satz, D(2)-Problem, seit 1941 offen
- [x] paper79: Temperley-Lieb / Jones-Polynom — Khovanov, Jones-Unknoten-Vermutung, Volumen-Vermutung

#### Batch 21 — Algebraische Geometrie & Algebra (Build 166)
- [x] paper80: Jacobian-Vermutung — Bass-Connell-Wright, Pinchuk (reell falsch), Dixmier-Äquivalenz
- [x] paper81: Hadwiger-Vermutung n≥7 — k≤6 bewiesen, Kostochka/Thomason, Kühn-Osthus
- [x] paper82: Andrews-Curtis-Vermutung — AC-Züge, Akbulut-Kirby, 4-Mannigfaltigkeitstopologie
- [x] paper83: Beal-Vermutung — FLT als Spezialfall, Darmon-Granville, abc-Verbindung

#### Batch 22 — Zahlentheorie (Galois & L-Funktionen) (Build 167)
- [x] paper84: Fontaine-Mazur-Vermutung — p-adische Galois-Darstellungen, Kisin/Emerton, Taylor
- [x] paper85: Lehmer-Problem / Mahler-Maß — Dobrowolski-Schranke, Deninger-Formel
- [x] paper86: Rang elliptischer Kurven (unbeschränkt) — Bhargava-Shankar, Goldfeld-Vermutung
- [x] paper87: Bateman-Horn-Vermutung — Singuläre Reihe, Green-Tao & Maynard-Tao (bewiesen)

#### Batch 23 — Zahlentheorie II (Build 168)
- [x] paper88: Ungerade vollkommene Zahlen — Euler-Struktur, Ochem-Rao n>10^1500
- [x] paper89: Mersenne-Primzahlen / gerade vollkommene Zahlen — Lucas-Lehmer, Lenstra-Wagstaff
- [x] paper90: Bunyakovsky-Vermutung — Iwaniec (n²+1 ist P₂), lokale Obstruktionen
- [x] paper91: Erdős-Gallai-Satz (bewiesen) + offene Verallgemeinerungen

#### Batch 24 — Geometrie & Differentialgeometrie (Build 169)
- [x] paper92: Donaldson 4-Mannigfaltigkeiten — Seiberg-Witten, exotisches ℝ⁴, glattes Poincaré offen
- [x] paper93: Gromov-Füllungsradius / Systolische Geometrie — Pu, Loewner, Katz
- [x] paper94: Hartshorne-Vermutungen — Quillen-Suslin (bewiesen), Horrocks-Mumford, Kodim-2 offen
- [x] paper95: Uniformisierung höherer Dimensionen — Yau-Vermutung, Mori-Siu-Yau (kompakt bewiesen)

#### Batch 25 — Algebraische Geometrie Final (Build 170)
- [x] paper96: Grothendieck-Standard-Vermutungen A/B/C/D — vollständiges Implikationsdiagramm, Motive
- [x] paper97: Kontsevich-Integral / Vassiliev-Invarianten — Drinfeld-Assoziator, Wheeling-Satz
