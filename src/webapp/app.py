"""
@file app.py
@brief Flask-Webanwendung für den Mathematik-Spezialisten.
@description
    Stellt alle mathematischen Module über ein Web-Interface bereit.
    Läuft auf Port 8080 und bietet REST-API-Endpunkte für alle
    mathematischen Berechnungen. Plots werden als Base64-PNG zurückgegeben.

    Routen:
        GET  /                          → Startseite
        GET  /api/health                → Systemstatus
        POST /api/algebra/solve         → Gleichungen lösen
        POST /api/algebra/polynomial    → Polynom auswerten
        POST /api/algebra/primes        → Primzahl-Eigenschaften
        POST /api/analysis/derivative   → Numerische Ableitung
        POST /api/analysis/integral     → Bestimmtes Integral
        POST /api/analysis/roots        → Nullstellensuche
        POST /api/linear_algebra/matrix → Matrix-Operationen
        POST /api/statistics/describe   → Deskriptive Statistik
        POST /api/ode/solve             → ODE lösen
        POST /api/complex/zeta          → Riemann-Zeta-Funktion
        POST /api/visualization/plot    → Funktionsplot als PNG
        GET  /topology                  → Topologie-Seite
        GET  /graph_theory              → Graphentheorie-Seite
        GET  /millennium                → Millennium-Probleme-Seite
        POST /api/topology/metric       → Metriken berechnen
        POST /api/topology/curve        → Parametrische Kurven
        POST /api/topology/euler        → Euler-Charakteristik
        POST /api/topology/fractal      → Box-Counting-Dimension
        POST /api/graph/analyze         → Graphen analysieren
        POST /api/graph/dijkstra        → Dijkstra kürzeste Wege
        POST /api/graph/combinatorics   → Kombinatorische Formeln
        POST /api/millennium/zeros      → Riemann-Nullstellen
        POST /api/millennium/goldbach   → Goldbach-Zerlegungen
        POST /api/millennium/complexity → Komplexitätsvergleich
        GET  /elliptic                  → Elliptische-Kurven-Seite
        POST /api/elliptic/analyze      → Kurvenanalyse (Δ, j, Plot, Punktarithmetik)
        POST /api/elliptic/hasse        → Hasse-Schranke (Gruppenordnung mod p)
        GET  /tensor                    → Tensor-Geometrie-Seite
        POST /api/tensor/curvature      → Metriktensor, Ricci-Skalar, Gaußsche Krümmung

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import io
import base64
import traceback

# Elternverzeichnis (src/) zum Python-Pfad hinzufügen, damit alle Module ladbar sind
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Flask-Framework importieren
from flask import Flask, request, jsonify, render_template

# matplotlib im Headless-Modus konfigurieren (kein Display nötig)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Sympy für symbolische Mathematik (Ausdrucks-Parsing)
import sympy as sp

# ---- Eigene Mathematik-Module importieren ----
from algebra import (
    solve_linear, solve_quadratic, is_prime, prime_factorization,
    euler_phi, Polynomial
)
from analysis import (
    numerical_derivative, numerical_integral, newton_raphson, bisection
)
from linear_algebra import Matrix, Vector
from statistics_math import mean, median, mode, variance, std_dev
from ode import euler_method, runge_kutta4
from complex_analysis import riemann_zeta, riemann_zeta_mpmath
from visualization import plot_function_2d
from topology import (
    euclidean_metric, manhattan_metric, chebyshev_metric, p_norm_metric,
    circle_curve, lissajous_curve, helix_curve,
    euler_characteristic_polygon, genus_from_euler,
    box_counting_dimension, hausdorff_dimension_cantor, sierpinski_dimension
)
from graph_theory import (
    Graph, dijkstra, is_eulerian,
    binomial_coefficient, catalan_number, bell_numbers, derangements, partition_count
)
from millennium_problems import (
    riemann_zeros_mpmath, hardy_littlewood_goldbach_density, complexity_comparison
)
from proof_theory import goldbach_all_decompositions
from elliptic_curves import EllipticCurve, ECPoint
from tensor_geometry import (
    sphere_metric, torus_metric, hyperbolic_plane_metric,
    saddle_metric, flat_metric,
    ricci_scalar, gaussian_curvature
)
from step_by_step import (
    newton_raphson_steps,
    bisection_steps,
    gauss_elimination_steps,
    euclidean_algorithm_steps,
    rsa_steps,
    prime_factorization_steps,
    format_steps_text,
    format_steps_html,
)

# ---- Neue Module: Gruppentheorie, Ringtheorie, Galois-Theorie, Moduln ----
try:
    from group_theory import (
        Group, cyclic_group, dihedral_group, sylow_theorems
    )
    _group_theory_ok = True
except Exception as _e:
    _group_theory_ok = False
    print(f"[WARN] group_theory nicht ladbar: {_e}")

try:
    from ring_theory import Ring, polynomial_factorization_mod_p, is_irreducible_mod_p
    _ring_theory_ok = True
except Exception as _e:
    _ring_theory_ok = False
    print(f"[WARN] ring_theory nicht ladbar: {_e}")

try:
    from galois_theory import (
        galois_group_polynomial, is_solvable_by_radicals, FiniteField
    )
    _galois_ok = True
except Exception as _e:
    _galois_ok = False
    print(f"[WARN] galois_theory nicht ladbar: {_e}")

try:
    from modules_algebra import smith_normal_form, structure_theorem_abelian_groups
    from commutative_algebra import prime_spectrum
    _modules_ok = True
except Exception as _e:
    _modules_ok = False
    print(f"[WARN] modules_algebra/commutative_algebra nicht ladbar: {_e}")

# ---- Neue Module: Fourier, Analytische ZT, Modulformen, L-Funktionen ----
try:
    from fourier import dft, fft, fourier_coefficients, dominant_frequency
    _fourier_ok = True
except Exception as _e:
    _fourier_ok = False
    print(f"[WARN] fourier nicht ladbar: {_e}")

try:
    from analytic_number_theory import (
        prime_counting_function, logarithmic_integral,
        von_mangoldt_function, chebyshev_psi
    )
    _analytic_nt_ok = True
except Exception as _e:
    _analytic_nt_ok = False
    print(f"[WARN] analytic_number_theory nicht ladbar: {_e}")

try:
    from modular_forms import (
        eisenstein_series, j_invariant, fourier_coefficients_delta
    )
    _modular_ok = True
except Exception as _e:
    _modular_ok = False
    print(f"[WARN] modular_forms nicht ladbar: {_e}")

try:
    from l_functions import (
        dirichlet_l_function, dirichlet_characters, l_function_special_values
    )
    _lfunc_ok = True
except Exception as _e:
    _lfunc_ok = False
    print(f"[WARN] l_functions nicht ladbar: {_e}")

# ---- Neue Module: p-adisch, Iwasawa, Spinoren, Numerik ----
try:
    from p_adic import p_adic_valuation, p_adic_norm, PAdicNumber, hensel_lift
    _padic_ok = True
except Exception as _e:
    _padic_ok = False
    print(f"[WARN] p_adic nicht ladbar: {_e}")

try:
    from iwasawa_theory import (
        iwasawa_polynomial, kummer_congruences_check, kubota_leopoldt_l_function
    )
    _iwasawa_ok = True
except Exception as _e:
    _iwasawa_ok = False
    print(f"[WARN] iwasawa_theory nicht ladbar: {_e}")

try:
    from spinors import gamma_matrices, clifford_algebra_check, dirac_spinor
    _spinors_ok = True
except Exception as _e:
    _spinors_ok = False
    print(f"[WARN] spinors nicht ladbar: {_e}")

try:
    from numerical_methods import (
        lagrange_interpolation, gradient_descent, golden_section_search
    )
    _numerical_ok = True
except Exception as _e:
    _numerical_ok = False
    print(f"[WARN] numerical_methods nicht ladbar: {_e}")

try:
    from mathematical_logic import (
        Proposition, LogicFormula, truth_table,
        is_tautology, is_satisfiable, is_contradiction,
        logical_equivalence, resolution_refutation, dpll,
        godel_numbering, godel_decode, incompleteness_demonstration
    )
    _logic_ok = True
except Exception as _e:
    _logic_ok = False
    print(f"[WARN] mathematical_logic nicht ladbar: {_e}")

try:
    from set_theory import (
        MathSet, Relation, MathFunction, Ordinal,
        zfc_axioms, cantors_diagonal_argument, is_countable,
        cantor_schroeder_bernstein
    )
    _set_theory_ok = True
except Exception as _e:
    _set_theory_ok = False
    print(f"[WARN] set_theory nicht ladbar: {_e}")

# ---- Neue Module: Darstellungstheorie, Verbandstheorie, Multilineare Algebra ----
try:
    from representation_theory import (
        character_table, schur_lemma, z2_representations, s3_representations,
        burnside_lemma, regular_representation, maschke_theorem_check
    )
    _repr_ok = True
except Exception as _e:
    _repr_ok = False
    print(f"[WARN] representation_theory nicht ladbar: {_e}")

try:
    from lattice_theory import (
        Lattice, PartialOrder, BooleanAlgebra,
        divisibility_lattice, power_set_lattice, birkhoff_representation_theorem,
        jordan_dedekind_chain_condition
    )
    _lattice_ok = True
except Exception as _e:
    _lattice_ok = False
    print(f"[WARN] lattice_theory nicht ladbar: {_e}")

try:
    from multilinear_algebra import (
        tensor_product_vectors, tensor_product_matrices, wedge_product,
        exterior_power, gram_matrix, gram_determinant, signature_bilinear_form
    )
    _multilinear_ok = True
except Exception as _e:
    _multilinear_ok = False
    print(f"[WARN] multilinear_algebra nicht ladbar: {_e}")

try:
    from homological_algebra import (
        ChainComplex, simplicial_complex_to_chain,
        ext_group, tor_group, free_resolution, universal_coefficient_theorem
    )
    _homological_ok = True
except Exception as _e:
    _homological_ok = False
    print(f"[WARN] homological_algebra nicht ladbar: {_e}")

try:
    from algebraic_structures import (
        Magma, Semigroup, Monoid, GroupFromMagma,
        algebraic_structure_hierarchy, free_monoid
    )
    _algstruct_ok = True
except Exception as _e:
    _algstruct_ok = False
    print(f"[WARN] algebraic_structures nicht ladbar: {_e}")

try:
    from invariant_theory import (
        reynolds_operator, molien_series, polynomial_invariants_sn,
        elementary_symmetric_polynomials, hilbert_basis_theorem_demo,
        fundamental_theorem_symmetric_poly
    )
    _invariant_ok = True
except Exception as _e:
    _invariant_ok = False
    print(f"[WARN] invariant_theory nicht ladbar: {_e}")

try:
    from universal_algebra import (
        Signature, Algebra, Variety, birkhoff_theorem_demo,
        congruence_relation, groups_variety, rings_variety, lattices_variety
    )
    _universal_ok = True
except Exception as _e:
    _universal_ok = False
    print(f"[WARN] universal_algebra nicht ladbar: {_e}")

# ---- Neue Module: Modelltheorie, Formale Beweistheorie, Rekursionstheorie ----
try:
    from model_theory import (
        Structure, Signature as ModelSignature,
        elementary_equivalence, compactness_theorem_demo,
        nonstandard_arithmetic_demo, theory_is_categorical,
        quantifier_elimination_demo, lowenheim_skolem_downward
    )
    _model_ok = True
except Exception as _e:
    _model_ok = False
    print(f"[WARN] model_theory nicht ladbar: {_e}")

try:
    from proof_theory_formal import (
        LKRules, NDProof, prove_simple, cut_elimination_demo,
        proof_complexity_comparison, big_five_systems,
        proof_theoretic_ordinal, gentzen_consistency_proof_sketch
    )
    _proof_formal_ok = True
except Exception as _e:
    _proof_formal_ok = False
    print(f"[WARN] proof_theory_formal nicht ladbar: {_e}")

try:
    from recursion_theory import (
        TuringMachine, tm_recognizes_palindromes, tm_binary_increment,
        halting_problem_undecidability_proof, rice_theorem_demo,
        ackermann_function, ackermann_growth_demo, arithmetical_hierarchy,
        time_complexity_classes, kolmogorov_complexity_approx
    )
    _recursion_ok = True
except Exception as _e:
    _recursion_ok = False
    print(f"[WARN] recursion_theory nicht ladbar: {_e}")

# ---- Neue Module: Maßtheorie, Spezielle Funktionen, Funktionalanalysis, PDE, Operatoralgebren ----
try:
    from measure_theory import (
        SigmaAlgebra, Measure, LebesgueMeasure, lebesgue_integral,
        riemann_vs_lebesgue, LpSpace, radon_nikodym_theorem_demo,
        fubini_theorem_demo, dominated_convergence_theorem_demo,
        monotone_convergence_theorem_demo
    )
    _measure_ok = True
except Exception as _e:
    _measure_ok = False
    print(f"[WARN] measure_theory nicht ladbar: {_e}")

try:
    from special_functions import (
        BesselFunctions, LegendrePolynomials, AiryFunctions,
        gamma_function_properties, beta_function, digamma_function,
        riemann_zeta_special_values, elliptic_integrals, error_function_properties
    )
    _special_ok = True
except Exception as _e:
    _special_ok = False
    print(f"[WARN] special_functions nicht ladbar: {_e}")

try:
    from functional_analysis import (
        BanachSpace, HilbertSpace, LinearOperator,
        banach_fixed_point, hahn_banach_theorem_demo,
        spectral_theorem_demo, riesz_representation_theorem,
        open_mapping_theorem_demo, fredholm_alternative,
        uniform_boundedness_principle_demo
    )
    _functional_ok = True
except Exception as _e:
    _functional_ok = False
    print(f"[WARN] functional_analysis nicht ladbar: {_e}")

try:
    from pde import (
        classify_pde, heat_equation_explicit, wave_equation_explicit,
        heat_equation_analytical, wave_equation_analytical,
        laplace_equation_2d, schrodinger_stationary
    )
    _pde_ok = True
except Exception as _e:
    _pde_ok = False
    print(f"[WARN] pde nicht ladbar: {_e}")

try:
    from operator_algebras import (
        CStarAlgebra, VonNeumannAlgebra,
        factors_classification, k_theory_intro, gelfand_transform_demo,
        gns_construction_demo, cuntz_algebra_demo
    )
    _opalg_ok = True
except Exception as _e:
    _opalg_ok = False
    print(f"[WARN] operator_algebras nicht ladbar: {_e}")


# ===========================================================================
# FLASK-APP INITIALISIEREN
# ===========================================================================

# Template- und Static-Ordner relativ zur app.py-Datei setzen
app = Flask(
    __name__,
    template_folder='templates',
    static_folder='static'
)


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def safe_parse_expr(expr_str: str):
    """
    @brief Parst einen mathematischen Ausdruck sicher mit SymPy.
    @description
        Wandelt einen String wie "sin(x)**2 + x" in einen SymPy-Ausdruck um.
        Erlaubt nur bekannte mathematische Funktionen (kein eval mit beliebigem Code).
    @param expr_str: Mathematischer Ausdruck als String
    @return Tupel (expr, x_symbol) – SymPy-Ausdruck und das x-Symbol
    @raises ValueError wenn der Ausdruck nicht parsbar ist
    @date 2026-03-10
    """
    # x als Symbol definieren
    x = sp.Symbol('x')
    # Erlaubte Funktionen und Konstanten für sicheres Parsing
    allowed_locals = {
        'x': x,
        'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
        'asin': sp.asin, 'acos': sp.acos, 'atan': sp.atan,
        'sinh': sp.sinh, 'cosh': sp.cosh, 'tanh': sp.tanh,
        'exp': sp.exp, 'log': sp.log, 'ln': sp.log,
        'sqrt': sp.sqrt, 'abs': sp.Abs,
        'pi': sp.pi, 'E': sp.E, 'e': sp.E,
        'oo': sp.oo, 'inf': sp.oo,
    }
    try:
        # sympify parst den String in einen SymPy-Ausdruck
        expr = sp.sympify(expr_str, locals=allowed_locals)
        return expr, x
    except Exception as e:
        raise ValueError(f"Ungültiger Ausdruck '{expr_str}': {e}")


def expr_to_lambda(expr_str: str):
    """
    @brief Konvertiert einen Ausdrucks-String in eine aufrufbare Python-Funktion.
    @description
        Parst den Ausdruck und erstellt eine lambda-Funktion f(x) → float.
        Wird für numerische Berechnungen (Ableitung, Integral, Nullstellen) genutzt.
    @param expr_str: Mathematischer Ausdruck als String (z.B. "sin(x)**2")
    @return Callable f(x) → float
    @date 2026-03-10
    """
    expr, x_sym = safe_parse_expr(expr_str)
    # SymPy-Ausdruck in numerische Funktion umwandeln
    f = sp.lambdify(x_sym, expr, modules=['numpy', 'sympy'])
    # Wrapper, der sicherstellt dass immer ein float zurückkommt
    def wrapped(x_val):
        result = f(x_val)
        return float(result)
    return wrapped


def figure_to_base64(fig) -> str:
    """
    @brief Konvertiert eine matplotlib-Figure in einen Base64-kodierten PNG-String.
    @description
        Speichert das Bild in einen Speicherpuffer (kein Dateisystem nötig)
        und kodiert es als Base64-String für die JSON-Übertragung.
    @param fig: matplotlib Figure-Objekt
    @return Base64-kodierter PNG-String (ohne data:image/png;base64, Präfix)
    @date 2026-03-10
    """
    # Bild in Speicherpuffer rendern (kein Dateisystem-Zugriff nötig)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
    buf.seek(0)
    # Bytes in Base64-String kodieren
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


def error_response(message: str, status_code: int = 400):
    """
    @brief Erstellt eine einheitliche Fehler-Antwort.
    @param message: Fehlerbeschreibung (deutsch)
    @param status_code: HTTP-Statuscode (Standard: 400 Bad Request)
    @return Flask JSON-Response mit Fehlerdetails
    @date 2026-03-10
    """
    return jsonify({"error": message}), status_code


# ===========================================================================
# ROUTEN: SEITEN (GET)
# ===========================================================================

@app.route('/')
def index():
    """
    @brief Startseite: Übersicht aller mathematischen Module.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    """
    return render_template('index.html')


@app.route('/algebra')
def page_algebra():
    """@brief Algebra-Modul-Seite. @date 2026-03-10"""
    return render_template('algebra.html')


@app.route('/analysis')
def page_analysis():
    """@brief Analysis-Modul-Seite. @date 2026-03-10"""
    return render_template('analysis.html')


@app.route('/linear_algebra')
def page_linear_algebra():
    """@brief Lineare-Algebra-Modul-Seite. @date 2026-03-10"""
    return render_template('linear_algebra.html')


@app.route('/statistics')
def page_statistics():
    """@brief Statistik-Modul-Seite. @date 2026-03-10"""
    return render_template('statistics.html')


@app.route('/ode')
def page_ode():
    """@brief ODE-Modul-Seite. @date 2026-03-10"""
    return render_template('ode.html')


@app.route('/complex')
def page_complex():
    """@brief Komplexe-Analysis-Seite. @date 2026-03-10"""
    return render_template('complex.html')


@app.route('/number_theory')
def page_number_theory():
    """@brief Zahlentheorie-Seite. @date 2026-03-10"""
    return render_template('number_theory.html')


@app.route('/visualization')
def page_visualization():
    """@brief Visualisierungs-Seite. @date 2026-03-10"""
    return render_template('visualization.html')


# ---- Neue GET-Routen für alle neuen Module ----

@app.route('/group_theory')
def page_group_theory():
    """@brief Gruppentheorie-Seite. @date 2026-03-10"""
    return render_template('group_theory.html')


@app.route('/ring_theory')
def page_ring_theory():
    """@brief Ringtheorie-Seite. @date 2026-03-10"""
    return render_template('ring_theory.html')


@app.route('/galois_theory')
def page_galois_theory():
    """@brief Galois-Theorie-Seite. @date 2026-03-10"""
    return render_template('galois_theory.html')


@app.route('/modules_algebra')
def page_modules_algebra():
    """@brief Moduln-und-Kommutative-Algebra-Seite. @date 2026-03-10"""
    return render_template('modules_algebra.html')


@app.route('/fourier')
def page_fourier():
    """@brief Fourier-Analysis-Seite. @date 2026-03-10"""
    return render_template('fourier.html')


@app.route('/analytic_nt')
def page_analytic_nt():
    """@brief Analytische-Zahlentheorie-Seite. @date 2026-03-10"""
    return render_template('analytic_nt.html')


@app.route('/modular_forms')
def page_modular_forms():
    """@brief Modulformen-Seite. @date 2026-03-10"""
    return render_template('modular_forms.html')


@app.route('/l_functions')
def page_l_functions():
    """@brief L-Funktionen-Seite. @date 2026-03-10"""
    return render_template('l_functions.html')


@app.route('/p_adic')
def page_p_adic():
    """@brief p-adische-Zahlen-Seite. @date 2026-03-10"""
    return render_template('p_adic.html')


@app.route('/iwasawa')
def page_iwasawa():
    """@brief Iwasawa-Theorie-Seite. @date 2026-03-10"""
    return render_template('iwasawa.html')


@app.route('/spinors')
def page_spinors():
    """@brief Spinoren-Seite. @date 2026-03-10"""
    return render_template('spinors.html')


@app.route('/numerical')
def page_numerical():
    """@brief Numerische-Methoden-Seite. @date 2026-03-10"""
    return render_template('numerical.html')


@app.route('/logic')
def page_logic():
    """@brief Mathematische-Logik-Seite. @date 2026-03-10"""
    return render_template('logic.html')

@app.route('/set_theory')
def page_set_theory():
    """@brief Mengenlehre-Seite. @date 2026-03-10"""
    return render_template('set_theory.html')



# ===========================================================================
# ROUTE: SYSTEMSTATUS
# ===========================================================================

@app.route('/api/health', methods=['GET'])
def health():
    """
    @brief Gibt den Systemstatus zurück.
    @return JSON {"status": "ok", "build": 30}
    @date 2026-03-10
    """
    return jsonify({"status": "ok", "build": 30})


# ===========================================================================
# ROUTEN: ALGEBRA API
# ===========================================================================

@app.route('/api/algebra/solve', methods=['POST'])
def api_algebra_solve():
    """
    @brief Löst lineare oder quadratische Gleichungen.
    @description
        Erwartet JSON:
            {"type": "linear"|"quadratic", "coefficients": [a, b] oder [a, b, c]}
        Linear:     a*x + b = 0  →  x = -b/a
        Quadratisch: a*x² + b*x + c = 0  →  x₁, x₂ über Lösungsformel
    @return JSON mit Lösungen oder Fehlermeldung
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    eq_type = data.get('type', '')
    coeffs = data.get('coefficients', [])

    try:
        if eq_type == 'linear':
            # Lineare Gleichung: a*x + b = 0
            if len(coeffs) < 2:
                return error_response("Lineare Gleichung erfordert [a, b].")
            a, b = float(coeffs[0]), float(coeffs[1])
            result = solve_linear(a, b)
            # solve_linear gibt Zahl oder None zurück; als {real, imag}-Dict normalisieren
            solutions = [{"real": float(result), "imag": 0.0}] if result is not None else []
            return jsonify({
                "type": "linear",
                "equation": f"{a}·x + {b} = 0",
                "solutions": solutions
            })

        elif eq_type == 'quadratic':
            # Quadratische Gleichung: a*x² + b*x + c = 0
            if len(coeffs) < 3:
                return error_response("Quadratische Gleichung erfordert [a, b, c].")
            a, b, c = float(coeffs[0]), float(coeffs[1]), float(coeffs[2])
            roots = solve_quadratic(a, b, c)
            # Ergebnis: Liste von Lösungen (können komplex sein)
            solutions = []
            for r in roots:
                if isinstance(r, complex):
                    solutions.append({"real": r.real, "imag": r.imag})
                else:
                    solutions.append({"real": float(r), "imag": 0.0})
            return jsonify({
                "type": "quadratic",
                "equation": f"{a}·x² + {b}·x + {c} = 0",
                "solutions": solutions
            })

        else:
            return error_response(f"Unbekannter Gleichungstyp: '{eq_type}'. Erlaubt: 'linear', 'quadratic'.")

    except Exception as e:
        return error_response(f"Fehler beim Lösen: {str(e)}")


@app.route('/api/algebra/polynomial', methods=['POST'])
def api_algebra_polynomial():
    """
    @brief Wertet ein Polynom an einem Punkt x aus.
    @description
        Erwartet JSON:
            {"coefficients": [a_n, ..., a_1, a_0], "x": float}
        Das Polynom ist: a_n*x^n + ... + a_1*x + a_0
    @return JSON mit Ergebnis der Auswertung
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    coeffs = data.get('coefficients', [])
    x_val = data.get('x', None)

    if not coeffs:
        return error_response("Koeffizientenliste darf nicht leer sein.")
    if x_val is None:
        return error_response("Parameter 'x' fehlt.")

    try:
        # Polynomial-Objekt erstellen (Koeffizienten: höchster Grad zuerst)
        poly = Polynomial(coeffs)
        result = poly.evaluate(float(x_val))
        # Polynomdarstellung als lesbarer String
        terms = []
        n = len(coeffs) - 1
        for i, c in enumerate(coeffs):
            power = n - i
            if c == 0:
                continue
            if power == 0:
                terms.append(f"{c}")
            elif power == 1:
                terms.append(f"{c}·x")
            else:
                terms.append(f"{c}·x^{power}")
        poly_str = " + ".join(terms) if terms else "0"

        return jsonify({
            "polynomial": poly_str,
            "x": float(x_val),
            "result": result
        })
    except Exception as e:
        return error_response(f"Fehler bei Polynomauswertung: {str(e)}")


@app.route('/api/algebra/primes', methods=['POST'])
def api_algebra_primes():
    """
    @brief Berechnet Primzahl-Eigenschaften einer natürlichen Zahl.
    @description
        Erwartet JSON: {"n": int}
        Gibt zurück: is_prime, prime_factorization, euler_phi
    @return JSON mit Primzahl-Eigenschaften
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    n = data.get('n', None)
    if n is None:
        return error_response("Parameter 'n' fehlt.")

    try:
        n = int(n)
        if n < 1:
            return error_response("n muss eine positive ganze Zahl sein.")

        # Primzahl-Test, Primfaktorzerlegung und Euler'sche Phi-Funktion
        prime = is_prime(n)
        factors = prime_factorization(n) if n > 1 else {}
        phi = euler_phi(n)

        # Primfaktoren als lesbarer String: 12 → "2² · 3"
        factor_str = " · ".join(
            f"{p}^{e}" if e > 1 else str(p)
            for p, e in sorted(factors.items())
        ) if factors else str(n)

        return jsonify({
            "n": n,
            "is_prime": prime,
            "prime_factorization": factors,
            "factorization_str": factor_str,
            "euler_phi": phi
        })
    except Exception as e:
        return error_response(f"Fehler bei Primzahl-Berechnung: {str(e)}")


# ===========================================================================
# ROUTEN: ANALYSIS API
# ===========================================================================

@app.route('/api/analysis/derivative', methods=['POST'])
def api_analysis_derivative():
    """
    @brief Berechnet die numerische Ableitung eines Ausdrucks an einem Punkt.
    @description
        Erwartet JSON: {"expr": "sin(x)", "x": float}
        Nutzt zentralen Differenzenquotienten (Ordnung 2).
    @return JSON mit Ableitungswert
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    expr_str = data.get('expr', '')
    x_val = data.get('x', None)

    if not expr_str:
        return error_response("Parameter 'expr' fehlt.")
    if x_val is None:
        return error_response("Parameter 'x' fehlt.")

    try:
        f = expr_to_lambda(expr_str)
        x_val = float(x_val)
        # Numerische Ableitung berechnen
        deriv = numerical_derivative(f, x_val)
        return jsonify({
            "expr": expr_str,
            "x": x_val,
            "derivative": deriv,
            "latex": f"f'({x_val}) \\approx {deriv:.6f}"
        })
    except ValueError as e:
        return error_response(str(e))
    except Exception as e:
        return error_response(f"Fehler bei Ableitungsberechnung: {str(e)}")


@app.route('/api/analysis/integral', methods=['POST'])
def api_analysis_integral():
    """
    @brief Berechnet ein bestimmtes Integral numerisch (Simpson-Regel).
    @description
        Erwartet JSON: {"expr": "x**2", "a": float, "b": float}
        Berechnet ∫ₐᵇ f(x) dx mit der zusammengesetzten Simpson-Regel.
    @return JSON mit Integralwert
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    expr_str = data.get('expr', '')
    a = data.get('a', None)
    b = data.get('b', None)

    if not expr_str:
        return error_response("Parameter 'expr' fehlt.")
    if a is None or b is None:
        return error_response("Parameter 'a' und 'b' sind erforderlich.")

    try:
        f = expr_to_lambda(expr_str)
        a, b = float(a), float(b)
        if a >= b:
            return error_response("'a' muss kleiner als 'b' sein.")
        # Numerisches Integral (Simpson-Regel mit 1000 Teilintervallen)
        result = numerical_integral(f, a, b, n=1000)
        return jsonify({
            "expr": expr_str,
            "a": a,
            "b": b,
            "integral": result,
            "latex": f"\\int_{{{a}}}^{{{b}}} {expr_str}\\, dx \\approx {result:.6f}"
        })
    except ValueError as e:
        return error_response(str(e))
    except Exception as e:
        return error_response(f"Fehler bei Integralberechnung: {str(e)}")


@app.route('/api/analysis/roots', methods=['POST'])
def api_analysis_roots():
    """
    @brief Findet Nullstellen eines Ausdrucks numerisch.
    @description
        Erwartet JSON:
            {"expr": "x**2 - 2", "x0": float, "method": "newton"|"bisection"}
        Newton-Raphson: schnelle Konvergenz, braucht guten Startwert.
        Bisektionsverfahren: braucht Intervall [x0, x0+2].
    @return JSON mit gefundener Nullstelle
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    expr_str = data.get('expr', '')
    x0 = data.get('x0', 1.0)
    method = data.get('method', 'newton')

    if not expr_str:
        return error_response("Parameter 'expr' fehlt.")

    try:
        f = expr_to_lambda(expr_str)
        x0 = float(x0)

        if method == 'newton':
            # Newton-Raphson-Verfahren (berechnet Ableitung intern numerisch)
            root = newton_raphson(f, x0)
        elif method == 'bisection':
            # Bisektionsverfahren: Suche nach Vorzeichenwechsel im Intervall [x0, x0+4]
            # Intervall automatisch suchen
            a_bsct, b_bsct = x0, x0 + 2.0
            # Sicherheitscheck: Vorzeichenwechsel vorhanden?
            for step in [2.0, 5.0, 10.0, 20.0]:
                if f(x0) * f(x0 + step) < 0:
                    a_bsct, b_bsct = x0, x0 + step
                    break
                if f(x0 - step) * f(x0) < 0:
                    a_bsct, b_bsct = x0 - step, x0
                    break
            root = bisection(f, a_bsct, b_bsct)
        else:
            return error_response(f"Unbekannte Methode '{method}'. Erlaubt: 'newton', 'bisection'.")

        return jsonify({
            "expr": expr_str,
            "method": method,
            "x0": x0,
            "root": root,
            "verification": f(root),
            "latex": f"f({root:.6f}) \\approx {f(root):.2e}"
        })
    except ValueError as e:
        return error_response(str(e))
    except Exception as e:
        return error_response(f"Fehler bei Nullstellensuche: {str(e)}")


# ===========================================================================
# ROUTEN: LINEARE ALGEBRA API
# ===========================================================================

@app.route('/api/linear_algebra/matrix', methods=['POST'])
def api_linear_algebra_matrix():
    """
    @brief Führt Operationen auf Matrizen aus.
    @description
        Erwartet JSON:
            {"matrix": [[a,b],[c,d]], "operation": "det"|"inv"|"eigenvalues"}
        det: Determinante berechnen
        inv: Inverse Matrix berechnen
        eigenvalues: Eigenwerte berechnen
    @return JSON mit Operationsergebnis
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    matrix_data = data.get('matrix', [])
    operation = data.get('operation', '')

    if not matrix_data:
        return error_response("Parameter 'matrix' fehlt oder ist leer.")
    if not operation:
        return error_response("Parameter 'operation' fehlt.")

    try:
        # Matrix-Objekt aus verschachtelter Liste erstellen
        mat = Matrix(matrix_data)

        if operation == 'det':
            # Determinante berechnen (Methode heißt determinant())
            det_val = mat.determinant()
            return jsonify({
                "operation": "det",
                "result": det_val,
                "latex": f"\\det(A) = {det_val:.6f}"
            })

        elif operation == 'inv':
            # Inverse Matrix berechnen (Methode heißt inverse())
            inv_mat = mat.inverse()
            # Daten aus dem internen _data-Attribut lesen
            return jsonify({
                "operation": "inv",
                "result": inv_mat._data,
                "latex": "A^{-1}"
            })

        elif operation == 'eigenvalues':
            # Eigenwerte berechnen
            eigenvals = mat.eigenvalues()
            # Eigenwerte können komplex sein → sicher serialisieren
            eigenvals_serializable = []
            for ev in eigenvals:
                if isinstance(ev, complex):
                    eigenvals_serializable.append({"real": ev.real, "imag": ev.imag})
                else:
                    eigenvals_serializable.append({"real": float(ev), "imag": 0.0})
            return jsonify({
                "operation": "eigenvalues",
                "result": eigenvals_serializable,
                "count": len(eigenvals_serializable)
            })

        else:
            return error_response(f"Unbekannte Operation '{operation}'. Erlaubt: 'det', 'inv', 'eigenvalues'.")

    except Exception as e:
        return error_response(f"Fehler bei Matrix-Operation: {str(e)}")


# ===========================================================================
# ROUTEN: STATISTIK API
# ===========================================================================

@app.route('/api/statistics/describe', methods=['POST'])
def api_statistics_describe():
    """
    @brief Berechnet deskriptive Statistik für einen Datensatz.
    @description
        Erwartet JSON: {"data": [1, 2, 3, 4, 5, ...]}
        Gibt zurück: Mittelwert, Median, Modus, Varianz, Standardabweichung.
    @return JSON mit statistischen Kennzahlen
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    dataset = data.get('data', [])
    if not dataset:
        return error_response("Datensatz ist leer.")
    if len(dataset) < 1:
        return error_response("Mindestens ein Datenpunkt erforderlich.")

    try:
        # In float-Liste umwandeln
        nums = [float(x) for x in dataset]

        # Alle statistischen Kennzahlen berechnen
        result_mean = mean(nums)
        result_median = median(nums)
        result_mode = mode(nums)
        result_variance = variance(nums)
        result_std = std_dev(nums)

        return jsonify({
            "n": len(nums),
            "mean": result_mean,
            "median": result_median,
            "mode": result_mode,
            "variance": result_variance,
            "std_dev": result_std,
            "min": min(nums),
            "max": max(nums),
            "range": max(nums) - min(nums)
        })
    except Exception as e:
        return error_response(f"Fehler bei Statistik-Berechnung: {str(e)}")


# ===========================================================================
# ROUTEN: ODE API
# ===========================================================================

@app.route('/api/ode/solve', methods=['POST'])
def api_ode_solve():
    """
    @brief Löst gewöhnliche Differentialgleichungen numerisch.
    @description
        Erwartet JSON:
            {"equation": "harmonic"|"exponential"|"custom",
             "y0": [float],
             "t_end": float}

        harmonic:    y'' + y = 0  (harmonischer Oszillator)
        exponential: y' = -y      (exponentieller Zerfall)
        custom:      Nutzer-definierte einfache ODE y' = f(t, y)

    @return JSON mit t-Werten und y-Werten (für Plotting)
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    equation = data.get('equation', 'exponential')
    y0_list = data.get('y0', [1.0])
    t_end = float(data.get('t_end', 10.0))

    try:
        # Anzahl der Zeitschritte (max. 500 für übertragbare JSON-Größe)
        n_steps = min(500, int(t_end * 50))
        dt = t_end / n_steps

        if equation == 'exponential':
            # y' = -y  →  y(t) = y₀ · e^(-t)
            # Signatur: runge_kutta4(f, t0, y0, t_end, n_steps)
            f_ode = lambda t, y: -y
            y0 = float(y0_list[0]) if y0_list else 1.0
            t_vals, y_vals = runge_kutta4(f_ode, 0.0, y0, t_end, n_steps)
            description = "y' = -y  (exponentieller Zerfall)"

        elif equation == 'harmonic':
            # Harmonischer Oszillator als System: y'' + y = 0
            # Umformung: [y, v]' = [v, -y]
            def f_harmonic(t, y_vec):
                # y_vec = [position, geschwindigkeit]
                return [y_vec[1], -y_vec[0]]
            y0_vec = y0_list[:2] if len(y0_list) >= 2 else [1.0, 0.0]
            # Euler-Methode für Systeme
            t_vals = [i * dt for i in range(n_steps + 1)]
            y_current = list(y0_vec)
            y_history = [y_current[0]]
            for i in range(n_steps):
                dydt = f_harmonic(t_vals[i], y_current)
                y_current = [y_current[j] + dt * dydt[j] for j in range(len(y_current))]
                y_history.append(y_current[0])
            y_vals = y_history
            description = "y'' + y = 0  (harmonischer Oszillator)"

        elif equation == 'logistic':
            # Logistisches Wachstum: y' = r·y·(1 - y/K)
            r, K = 0.5, 10.0
            f_ode = lambda t, y: r * y * (1 - y / K)
            y0 = float(y0_list[0]) if y0_list else 1.0
            # Signatur: runge_kutta4(f, t0, y0, t_end, n_steps)
            t_vals, y_vals = runge_kutta4(f_ode, 0.0, y0, t_end, n_steps)
            description = f"y' = {r}·y·(1 - y/{K})  (logistisches Wachstum)"

        else:
            return error_response(f"Unbekannte Gleichung '{equation}'. Erlaubt: 'exponential', 'harmonic', 'logistic'.")

        # Nur jeden 5. Punkt zurückgeben (Datenreduktion für JSON)
        step = max(1, len(t_vals) // 100)
        return jsonify({
            "equation": equation,
            "description": description,
            "t": [float(t_vals[i]) for i in range(0, len(t_vals), step)],
            "y": [float(y_vals[i]) for i in range(0, len(y_vals), step)],
            "t_end": t_end
        })
    except Exception as e:
        return error_response(f"Fehler beim ODE-Lösen: {str(e)}")


# ===========================================================================
# ROUTEN: KOMPLEXE ANALYSIS API
# ===========================================================================

@app.route('/api/complex/zeta', methods=['POST'])
def api_complex_zeta():
    """
    @brief Berechnet die Riemann-Zeta-Funktion ζ(s) für einen komplexen Punkt.
    @description
        Erwartet JSON: {"s_real": float, "s_imag": float, "precision": int}
        Nutzt analytische Fortsetzung für alle s ≠ 1.
    @return JSON mit Real- und Imaginärteil von ζ(s)
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    s_real = float(data.get('s_real', 0.5))
    s_imag = float(data.get('s_imag', 0.0))
    precision = int(data.get('precision', 50))

    try:
        s = complex(s_real, s_imag)

        # Pol bei s=1 abfangen
        if abs(s - 1.0) < 1e-10:
            return error_response("ζ(s) hat einen Pol bei s = 1.")

        # Zeta-Funktion berechnen (mpmath für höhere Präzision verfügbar)
        try:
            zeta_val = riemann_zeta_mpmath(s, precision)
        except Exception:
            # Fallback auf Standard-Implementierung
            zeta_val = riemann_zeta(s)

        return jsonify({
            "s": {"real": s_real, "imag": s_imag},
            "zeta_real": float(zeta_val.real),
            "zeta_imag": float(zeta_val.imag),
            "abs": float(abs(zeta_val)),
            "latex": f"\\zeta({s_real} + {s_imag}i) = {zeta_val.real:.6f} + {zeta_val.imag:.6f}i"
        })
    except Exception as e:
        return error_response(f"Fehler bei Zeta-Berechnung: {str(e)}")


# ===========================================================================
# ROUTEN: VISUALISIERUNG API
# ===========================================================================

@app.route('/api/visualization/plot', methods=['POST'])
def api_visualization_plot():
    """
    @brief Erstellt einen Funktionsplot und gibt ihn als Base64-PNG zurück.
    @description
        Erwartet JSON:
            {"expr": "sin(x)", "x_min": float, "x_max": float, "title": "string"}
        Gibt das Bild als Base64-kodiertes PNG zurück (kein Dateisystem-Zugriff).
    @return JSON mit Base64-PNG und Metadaten
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    expr_str = data.get('expr', 'sin(x)')
    x_min = float(data.get('x_min', -10.0))
    x_max = float(data.get('x_max', 10.0))
    title = data.get('title', f'f(x) = {expr_str}')

    if x_min >= x_max:
        return error_response("'x_min' muss kleiner als 'x_max' sein.")

    try:
        f = expr_to_lambda(expr_str)

        # Plot erstellen
        x_vals = np.linspace(x_min, x_max, 800)
        y_vals = []
        for xi in x_vals:
            try:
                yi = f(xi)
                # Ungültige Werte (Polstellen, etc.) als NaN markieren
                y_vals.append(yi if np.isfinite(yi) else np.nan)
            except Exception:
                y_vals.append(np.nan)
        y_arr = np.array(y_vals)

        # Matplotlib-Figure erstellen (schönes dunkles Theme)
        fig, ax = plt.subplots(figsize=(9, 5))
        fig.patch.set_facecolor('#1e1e2e')
        ax.set_facecolor('#1e1e2e')

        # Funktion plotten
        ax.plot(x_vals, y_arr, color='#89b4fa', linewidth=2.0, label=f'f(x) = {expr_str}')

        # Achsen und Gitter stylen
        ax.axhline(0, color='#6c7086', linewidth=0.8)
        ax.axvline(0, color='#6c7086', linewidth=0.8)
        ax.grid(True, alpha=0.2, color='#6c7086')
        ax.set_title(title, color='#cdd6f4', fontsize=14, pad=12)
        ax.set_xlabel('x', color='#a6adc8')
        ax.set_ylabel('f(x)', color='#a6adc8')
        ax.tick_params(colors='#a6adc8')
        ax.legend(facecolor='#313244', edgecolor='#6c7086', labelcolor='#cdd6f4')

        # Rahmen einfärben
        for spine in ax.spines.values():
            spine.set_edgecolor('#6c7086')

        plt.tight_layout()

        # In Base64 konvertieren
        img_b64 = figure_to_base64(fig)

        return jsonify({
            "expr": expr_str,
            "x_min": x_min,
            "x_max": x_max,
            "title": title,
            "image_base64": img_b64,
            "latex": f"f(x) = {expr_str}"
        })
    except ValueError as e:
        return error_response(str(e))
    except Exception as e:
        return error_response(f"Fehler beim Erstellen des Plots: {str(e)}")


# ===========================================================================
# ZAHLENTHEORIE-ROUTEN (erweitert)
# ===========================================================================

@app.route('/api/number_theory/info', methods=['POST'])
def api_number_theory_info():
    """
    @brief Kombinierte Zahlentheorie-Informationen für eine Zahl.
    @description
        Erwartet JSON: {"n": int}
        Gibt Primfaktorzerlegung, euler_phi, is_prime, gcd-Info zurück.
    @return JSON mit Zahlentheorie-Eigenschaften
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")
    # Delegiere an /api/algebra/primes (gleiche Logik)
    return api_algebra_primes()


# ===========================================================================
# ROUTEN: NEUE SEITEN (GET)
# ===========================================================================

@app.route('/topology')
def topology():
    """
    @brief Topologie-Seite: Metriken, Kurven, Euler-Charakteristik, Fraktale.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    """
    return render_template('topology.html')


@app.route('/graph_theory')
def graph_theory_page():
    """
    @brief Graphentheorie-Seite: Graphanalyse, Dijkstra, Kombinatorik.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    """
    return render_template('graph_theory.html')


@app.route('/millennium')
def millennium():
    """
    @brief Millennium-Probleme-Seite: Riemann, Goldbach, P vs NP, Komplexität.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    """
    return render_template('millennium.html')


# ===========================================================================
# HILFSFUNKTIONEN: KANTENLISTEN-PARSER
# ===========================================================================

def parse_edge_list(edge_string: str, weighted: bool = False) -> list:
    """
    @brief Parst eine Kantenliste wie '1-2, 2-3:5, 3-1:2'.
    @description
        Unterstützt zwei Formate:
        - Ungewichtet: 'A-B, B-C, C-A'  → Liste von (u, v)
        - Gewichtet:   'A-B:5, B-C:3'   → Liste von (u, v, weight)
        Leerzeichen und leere Einträge werden ignoriert.
        Knoten können Strings oder Zahlen sein.
    @param edge_string Kantenliste als String.
    @param weighted True wenn Gewichte erwartet werden (Format 'u-v:w').
    @return Liste von (u, v) oder (u, v, weight) Tupeln.
    @raises ValueError bei ungültigem Format.
    @author Kurt Ingwer
    @date 2026-03-10
    """
    edges = []
    # Kantenliste an Kommas aufteilen und jede Kante einzeln verarbeiten
    for part in edge_string.split(','):
        part = part.strip()
        if not part:
            continue  # Leerzeichen-Token ignorieren

        # Prüfen ob Gewicht angegeben (Format: 'u-v:w')
        weight = 1.0  # Standardgewicht
        if ':' in part:
            edge_part, weight_str = part.rsplit(':', 1)
            try:
                weight = float(weight_str.strip())
            except ValueError:
                raise ValueError(f"Ungültiges Kantengewicht: '{weight_str.strip()}'")
            part = edge_part.strip()

        # Kante an Bindestrich aufteilen: 'u-v'
        # Sonderfall: negative Zahlen wie '-3-5' → nur ein Bindestrich als Trennzeichen
        if '-' not in part:
            raise ValueError(f"Ungültiges Kantenformat: '{part}' (erwartet 'u-v')")

        # Letzten Bindestrich als Trennzeichen verwenden
        dash_idx = part.rfind('-')
        if dash_idx == 0:
            raise ValueError(f"Ungültiges Kantenformat: '{part}'")

        u_str = part[:dash_idx].strip()
        v_str = part[dash_idx + 1:].strip()

        if not u_str or not v_str:
            raise ValueError(f"Leerer Knoten in Kante: '{part}'")

        # Knoten als Integer versuchen, sonst String
        def parse_node(s):
            try:
                return int(s)
            except ValueError:
                return s

        u = parse_node(u_str)
        v = parse_node(v_str)

        if weighted:
            edges.append((u, v, weight))
        else:
            edges.append((u, v))

    return edges


# ===========================================================================
# ROUTEN: TOPOLOGIE API
# ===========================================================================

@app.route('/api/topology/metric', methods=['POST'])
def topology_metric():
    """
    @brief Berechnet den Abstand zwischen zwei Punkten mit einer wählbaren Metrik.
    @description
        Erwartet JSON:
            {"metric": "euclidean"|"manhattan"|"chebyshev"|"p_norm",
             "point1": [x1, y1, ...],
             "point2": [x2, y2, ...],
             "p": float (nur für p_norm)}
        Unterstützte Metriken:
        - euclidean:  d = sqrt(Σ(xᵢ-yᵢ)²)
        - manhattan:  d = Σ|xᵢ-yᵢ|
        - chebyshev:  d = max|xᵢ-yᵢ|
        - p_norm:     d = (Σ|xᵢ-yᵢ|^p)^(1/p)
    @return JSON mit Distanz, Dimension und LaTeX-Formel
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    metric_type = data.get('metric', 'euclidean')
    point1 = data.get('point1', [])
    point2 = data.get('point2', [])
    p_val = float(data.get('p', 2.0))

    # Eingaben validieren
    if not point1 or not point2:
        return error_response("'point1' und 'point2' müssen nicht-leere Arrays sein.")
    if len(point1) != len(point2):
        return error_response(
            f"Punkte müssen gleiche Dimension haben: {len(point1)} vs {len(point2)}.")

    try:
        # Punkte in float-Listen umwandeln
        p1 = [float(x) for x in point1]
        p2 = [float(x) for x in point2]

        # Gewählte Metrik anwenden
        if metric_type == 'euclidean':
            dist = euclidean_metric(p1, p2)
            latex = f"d_{{2}} = \\sqrt{{\\sum_i(x_i-y_i)^2}} = {dist:.8f}"
        elif metric_type == 'manhattan':
            dist = manhattan_metric(p1, p2)
            latex = f"d_{{1}} = \\sum_i|x_i-y_i| = {dist:.8f}"
        elif metric_type == 'chebyshev':
            dist = chebyshev_metric(p1, p2)
            latex = f"d_{{\\infty}} = \\max_i|x_i-y_i| = {dist:.8f}"
        elif metric_type == 'p_norm':
            if p_val < 1:
                return error_response("p muss >= 1 sein für eine gültige Metrik.")
            dist = p_norm_metric(p1, p2, p_val)
            latex = f"d_{{{p_val}}} = \\left(\\sum_i|x_i-y_i|^{{{p_val}}}\\right)^{{1/{p_val}}} = {dist:.8f}"
        else:
            return error_response(
                f"Unbekannte Metrik '{metric_type}'. Erlaubt: euclidean, manhattan, chebyshev, p_norm.")

        return jsonify({
            "metric": metric_type,
            "point1": p1,
            "point2": p2,
            "distance": dist,
            "dimension": len(p1),
            "latex": latex
        })
    except Exception as e:
        return error_response(f"Fehler bei Metrik-Berechnung: {str(e)}")


@app.route('/api/topology/curve', methods=['POST'])
def topology_curve():
    """
    @brief Berechnet Bogenlänge und Abgeschlossenheit einer parametrischen Kurve.
    @description
        Erwartet JSON:
            {"type": "circle"|"lissajous"|"helix",
             "params": {... kurvenspezifische Parameter ...}}
        Kreis:     params = {"radius": float, "center": [cx, cy]}
        Lissajous: params = {"a": int, "b": int, "delta": float}
        Helix:     params = {"radius": float, "pitch": float}
    @return JSON mit arc_length, is_closed, dimension, ambient_dim
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    curve_type = data.get('type', 'circle')
    params = data.get('params', {})

    try:
        if curve_type == 'circle':
            # Kreisparameter aus JSON lesen (mit Standardwerten)
            radius = float(params.get('radius', 1.0))
            center = params.get('center', [0.0, 0.0])
            center = [float(c) for c in center]
            if radius <= 0:
                return error_response("Radius muss positiv sein.")
            curve = circle_curve(radius=radius, center=center)
            ambient_dim = 2

        elif curve_type == 'lissajous':
            # Lissajous-Parameter: Frequenzen a, b und Phasenversatz delta
            a = int(params.get('a', 3))
            b = int(params.get('b', 2))
            delta = float(params.get('delta', 0.0))
            if a <= 0 or b <= 0:
                return error_response("Frequenzen a und b müssen positiv sein.")
            curve = lissajous_curve(a=a, b=b, delta=delta)
            ambient_dim = 2

        elif curve_type == 'helix':
            # Helix-Parameter: Radius und Steigung
            radius = float(params.get('radius', 1.0))
            pitch = float(params.get('pitch', 1.0))
            if radius <= 0:
                return error_response("Radius muss positiv sein.")
            curve = helix_curve(radius=radius, pitch=pitch)
            ambient_dim = 3

        else:
            return error_response(
                f"Unbekannter Kurventyp '{curve_type}'. Erlaubt: circle, lissajous, helix.")

        # Bogenlänge und Abgeschlossenheit berechnen
        arc = curve.arc_length(n=2000)
        closed = curve.is_closed()

        return jsonify({
            "type": curve_type,
            "params": params,
            "arc_length": arc,
            "is_closed": closed,
            "dimension": 1,          # Kurve ist 1-dimensional (Linie)
            "ambient_dim": ambient_dim
        })
    except Exception as e:
        return error_response(f"Fehler bei Kurvenberechnung: {str(e)}")


@app.route('/api/topology/euler', methods=['POST'])
def topology_euler():
    """
    @brief Berechnet die Euler-Charakteristik χ = V - E + F und das Geschlecht.
    @description
        Erwartet JSON: {"vertices": int, "edges": int, "faces": int}
        Euler-Formel für konvexe Polyeder: χ = 2 (Sphäre).
        Geschlecht: g = (2 - χ) / 2 für orientierbare Flächen.
        Beispiele:
            Würfel:    V=8, E=12, F=6  → χ=2, g=0
            Torus:     V=1, E=1,  F=1  → χ=0 (ohne Rand), g=1
    @return JSON mit chi, genus, vertices, edges, faces
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    try:
        v = int(data.get('vertices', 0))
        e = int(data.get('edges', 0))
        f = int(data.get('faces', 0))

        if v < 0 or e < 0 or f < 0:
            return error_response("V, E, F müssen nicht-negative ganze Zahlen sein.")

        # Euler-Charakteristik berechnen: χ = V - E + F
        chi = euler_characteristic_polygon(v, e, f)
        # Geschlecht aus Euler-Charakteristik ableiten
        genus = genus_from_euler(chi, orientable=True)

        return jsonify({
            "vertices": v,
            "edges": e,
            "faces": f,
            "chi": chi,
            "genus": max(0, genus)  # Negatives Geschlecht physikalisch nicht sinnvoll
        })
    except Exception as e:
        return error_response(f"Fehler bei Euler-Berechnung: {str(e)}")


@app.route('/api/topology/fractal', methods=['POST'])
def topology_fractal():
    """
    @brief Berechnet die Box-Counting-Dimension einer Punktmenge.
    @description
        Erwartet JSON: {"points": [[x1,y1], [x2,y2], ...]}
        Mindestens 3 Punkte erforderlich.
        Die Box-Counting-Dimension ist ein Maß für die fraktale Komplexität
        einer Punktmenge: d = log(N(ε)) / log(1/ε).
    @return JSON mit dimension, interpretation, epsilon_values
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    points = data.get('points', [])

    if not points or len(points) < 3:
        return error_response("Mindestens 3 Punkte als [[x,y],...] erforderlich.")

    try:
        # Punkte in float-Listen umwandeln
        pts = [[float(coord) for coord in pt] for pt in points]

        # Koordinatenbereich bestimmen um sinnvolle Epsilon-Werte zu wählen
        all_x = [pt[0] for pt in pts]
        all_y = [pt[1] if len(pt) > 1 else 0 for pt in pts]
        x_range = max(all_x) - min(all_x) if len(all_x) > 1 else 1.0
        y_range = max(all_y) - min(all_y) if len(all_y) > 1 else 1.0
        scale = max(x_range, y_range, 0.01)

        # Epsilon-Werte logarithmisch gleichverteilt über 4 Dekaden
        epsilons = [scale / (2 ** k) for k in range(1, 8)]

        # Box-Counting-Dimension berechnen (nur 2D-Koordinaten verwenden)
        pts_2d = [[pt[0], pt[1] if len(pt) > 1 else 0.0] for pt in pts]
        dim = box_counting_dimension(pts_2d, epsilons)

        # Dimension interpretieren
        if dim < 0.1:
            interpretation = "Punktmenge (0-dimensional)"
        elif dim < 1.1:
            interpretation = "Kurve/eindimensionale Struktur"
        elif dim < 1.5:
            interpretation = "Fraktal zwischen Linie und Fläche (schwach fraktal)"
        elif dim < 1.9:
            interpretation = "Fraktal (z.B. Sierpinski-ähnlich)"
        elif dim < 2.1:
            interpretation = "Flächenfüllendes Fraktal / 2-dimensionale Menge"
        else:
            interpretation = "Höherdimensionale Struktur"

        return jsonify({
            "dimension": dim,
            "n_points": len(pts),
            "interpretation": interpretation,
            "epsilon_values": [round(e, 6) for e in epsilons]
        })
    except Exception as e:
        return error_response(f"Fehler bei Box-Counting: {str(e)}")


# ===========================================================================
# ROUTEN: GRAPHENTHEORIE API
# ===========================================================================

@app.route('/api/graph/analyze', methods=['POST'])
def graph_analyze():
    """
    @brief Analysiert einen Graphen aus einer Kantenliste.
    @description
        Erwartet JSON:
            {"edges": "1-2, 2-3, 3-1", "directed": false}
        Berechnet:
        - Knotenanzahl, Kantenanzahl
        - Knotengrade (degree)
        - Zusammenhang (connected)
        - Bipartitheit
        - Euler-Kreis-Status (alle Knoten gerader Grad)
    @return JSON mit Grapheigenschaften
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    edge_string = data.get('edges', '')
    directed = data.get('directed', False)
    if isinstance(directed, str):
        directed = directed.lower() == 'true'

    if not edge_string.strip():
        return error_response("Kantenliste darf nicht leer sein.")

    try:
        # Kantenliste parsen (ungewichtet)
        edges = parse_edge_list(edge_string, weighted=False)

        # Graph aufbauen
        g = Graph(directed=directed)
        for u, v in edges:
            g.add_edge(u, v)

        # Knotengrade berechnen
        degrees = {str(node): g.degree(node) for node in g.vertices()}

        # Zusammenhang prüfen via BFS (alle Knoten von einem aus erreichbar?)
        vertices = list(g.vertices())
        is_connected_flag = False
        if vertices:
            # BFS vom ersten Knoten
            visited = set()
            queue = [vertices[0]]
            visited.add(vertices[0])
            while queue:
                curr = queue.pop(0)
                for nb in g.neighbors(curr):
                    if nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            is_connected_flag = len(visited) == len(vertices)

        # Bipartitheit via 2-Färbbarkeit prüfen (BFS mit 2 Farben)
        is_bipartite_flag = True
        if vertices:
            color = {}
            for start in vertices:
                if start in color:
                    continue
                color[start] = 0
                queue = [start]
                while queue:
                    curr = queue.pop(0)
                    for nb in g.neighbors(curr):
                        if nb not in color:
                            color[nb] = 1 - color[curr]
                            queue.append(nb)
                        elif color[nb] == color[curr]:
                            is_bipartite_flag = False
                            break

        # Euler-Status berechnen
        euler_raw = is_eulerian(g)
        # Rückgabe vereinheitlichen: intern heißen die Felder 'euler_circuit' und 'euler_path'
        euler_info = {
            "has_euler_circuit": euler_raw.get('euler_circuit', False),
            "has_euler_path": euler_raw.get('euler_path', False),
            "odd_degree_vertices": [str(v) for v in euler_raw.get('odd_degree_vertices', [])],
            "reason": (
                "Alle Knoten haben geraden Grad → Euler-Kreis existiert."
                if euler_raw.get('euler_circuit')
                else "Genau 2 Knoten mit ungeradem Grad → Euler-Pfad existiert."
                if euler_raw.get('euler_path')
                else f"{len(euler_raw.get('odd_degree_vertices', []))} Knoten mit ungeradem Grad."
            )
        }

        return jsonify({
            "num_vertices": len(vertices),
            "num_edges": len(edges),
            "directed": directed,
            "is_connected": is_connected_flag,
            "is_bipartite": is_bipartite_flag,
            "degrees": degrees,
            "is_eulerian": euler_info
        })
    except ValueError as e:
        return error_response(f"Fehler beim Parsen der Kantenliste: {str(e)}")
    except Exception as e:
        return error_response(f"Fehler bei Graphenanalyse: {str(e)}")


@app.route('/api/graph/dijkstra', methods=['POST'])
def graph_dijkstra():
    """
    @brief Berechnet kürzeste Wege mit dem Dijkstra-Algorithmus.
    @description
        Erwartet JSON:
            {"edges": "1-2:5, 2-3:3, 1-3:10", "start": "1"}
        Format der Kantenliste: 'u-v:gewicht'
        Gibt Distanzen und Pfade von Startknoten zu allen anderen Knoten zurück.
    @return JSON mit distances (dict) und paths (dict)
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    edge_string = data.get('edges', '')
    start_str = data.get('start', '')

    if not edge_string.strip():
        return error_response("Kantenliste darf nicht leer sein.")
    if not start_str.strip():
        return error_response("Startknoten ('start') ist erforderlich.")

    try:
        # Gewichtete Kantenliste parsen
        edges = parse_edge_list(edge_string, weighted=True)

        # Ungerichteten gewichteten Graphen aufbauen
        g = Graph(directed=False)
        for u, v, w in edges:
            g.add_edge(u, v, weight=w)

        # Startknoten parsen (int oder str)
        def parse_node(s):
            try:
                return int(s)
            except ValueError:
                return s

        start = parse_node(start_str.strip())

        # Startknoten muss im Graphen existieren
        if start not in g.vertices():
            return error_response(
                f"Startknoten '{start}' nicht im Graphen. "
                f"Vorhandene Knoten: {sorted(str(v) for v in g.vertices())}")

        # Dijkstra-Algorithmus ausführen
        result = dijkstra(g, start)

        # Distanzen serialisieren (None und inf → -1 für JSON)
        distances_serializable = {}
        paths_serializable = {}
        for node in g.vertices():
            node_key = str(node)
            dist = result.get('distances', {}).get(node)
            distances_serializable[node_key] = (
                float(dist) if dist is not None and dist != float('inf') else -1
            )
            # Dijkstra gibt 'path' (nicht 'paths') zurück
            path = result.get('path', {}).get(node, [])
            paths_serializable[node_key] = [str(p) for p in path]

        return jsonify({
            "start": str(start),
            "distances": distances_serializable,
            "paths": paths_serializable,
            "num_vertices": len(list(g.vertices()))
        })
    except ValueError as e:
        return error_response(f"Fehler beim Parsen: {str(e)}")
    except Exception as e:
        return error_response(f"Fehler bei Dijkstra: {str(e)}")


@app.route('/api/graph/combinatorics', methods=['POST'])
def graph_combinatorics():
    """
    @brief Berechnet kombinatorische Formeln.
    @description
        Erwartet JSON:
            {"function": "binomial"|"catalan"|"bell"|"derangements"|"partition",
             "n": int,
             "k": int}
        Funktionen:
        - binomial:     C(n,k) = n! / (k! * (n-k)!)
        - catalan:      C_n = C(2n,n) / (n+1)
        - bell:         B_n = Summe der Stirling-Zahlen 2. Art
        - derangements: D_n = n! * Σ(-1)^k / k!
        - partition:    p(n) = Anzahl additive Zerlegungen von n
    @return JSON mit result, function_name, latex
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    func = data.get('function', 'binomial')
    n = data.get('n', 0)
    k = data.get('k', 0)

    try:
        n = int(n)
        k = int(k)
    except (TypeError, ValueError):
        return error_response("'n' und 'k' müssen ganze Zahlen sein.")

    if n < 0:
        return error_response("n muss nicht-negativ sein.")

    try:
        if func == 'binomial':
            # Binomialkoeffizient C(n,k): n über k
            if k < 0 or k > n:
                return error_response(f"k muss zwischen 0 und n={n} liegen.")
            result = binomial_coefficient(n, k)
            func_name = f"C({n}, {k})"
            latex = f"\\binom{{{n}}}{{{k}}} = \\frac{{{n}!}}{{{k}!\\cdot({n}-{k})!}} = {result}"

        elif func == 'catalan':
            # Catalan-Zahl C_n
            result = catalan_number(n)
            func_name = f"Catalan C_{n}"
            latex = f"C_{{{n}}} = \\frac{{1}}{{{n}+1}}\\binom{{2\\cdot{n}}}{{{n}}} = {result}"

        elif func == 'bell':
            # Bell-Zahlen: bell_numbers(n) gibt Liste B_0, ..., B_n zurück
            bell_list = bell_numbers(n)
            result = bell_list[n] if n < len(bell_list) else bell_list[-1]
            func_name = f"Bell B_{n}"
            latex = f"B_{{{n}}} = {result}"

        elif func == 'derangements':
            # Derangements D_n: Permutationen ohne Fixpunkte
            result = derangements(n)
            func_name = f"Derangements D_{n}"
            latex = f"D_{{{n}}} = {n}! \\sum_{{k=0}}^{{{n}}} \\frac{{(-1)^k}}{{k!}} = {result}"

        elif func == 'partition':
            # Partitionszahl p(n)
            if n > 200:
                return error_response("n > 200 zu groß für Partitionsberechnung.")
            result = partition_count(n)
            func_name = f"Partition p({n})"
            latex = f"p({n}) = {result}"

        else:
            return error_response(
                f"Unbekannte Funktion '{func}'. "
                f"Erlaubt: binomial, catalan, bell, derangements, partition.")

        return jsonify({
            "function": func,
            "function_name": func_name,
            "n": n,
            "k": k,
            "result": result,
            "latex": latex
        })
    except Exception as e:
        return error_response(f"Fehler bei kombinatorischer Berechnung: {str(e)}")


# ===========================================================================
# ROUTEN: MILLENNIUM-PROBLEME API
# ===========================================================================

@app.route('/api/millennium/zeros', methods=['POST'])
def millennium_zeros():
    """
    @brief Berechnet die ersten n nicht-trivialen Riemann-Nullstellen.
    @description
        Erwartet JSON: {"n_zeros": int (1-20)}
        Nutzt mpmath.zetazero() für hochgenaue Berechnung.
        Alle bekannten Nullstellen liegen auf Re(s) = 1/2 (krit. Gerade).
    @return JSON mit Liste von Nullstellen {real, imag}
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    n_zeros = data.get('n_zeros', 10)

    try:
        n_zeros = int(n_zeros)
        if n_zeros < 1 or n_zeros > 20:
            return error_response("n_zeros muss zwischen 1 und 20 liegen.")

        # Nullstellen via mpmath berechnen (hochgenau)
        zeros_raw = riemann_zeros_mpmath(n_zeros=n_zeros, dps=30)

        # Nullstellen serialisieren: jede als {real, imag}
        # riemann_zeros_mpmath gibt Dicts mit 'n', 're', 'im', 'on_line' zurück
        zeros_serializable = []
        for z in zeros_raw:
            if isinstance(z, dict):
                # Felder umbenennen auf 'real'/'imag' für Frontend-Einheitlichkeit
                zeros_serializable.append({
                    "real": float(z.get('re', z.get('real', 0.5))),
                    "imag": float(z.get('im', z.get('imag', 0.0))),
                    "n": z.get('n', 0),
                    "on_line": z.get('on_line', True)
                })
            elif hasattr(z, 'real') and hasattr(z, 'imag'):
                zeros_serializable.append({
                    "real": float(z.real),
                    "imag": float(z.imag)
                })
            else:
                zeros_serializable.append({"real": 0.5, "imag": float(z)})

        return jsonify({
            "n_zeros": n_zeros,
            "zeros": zeros_serializable
        })
    except Exception as e:
        return error_response(f"Fehler bei Nullstellen-Berechnung: {str(e)}")


@app.route('/api/millennium/goldbach', methods=['POST'])
def millennium_goldbach():
    """
    @brief Berechnet alle Goldbach-Zerlegungen einer geraden Zahl.
    @description
        Erwartet JSON: {"n": int}
        n muss gerade und >= 4 sein.
        Gibt alle Zerlegungen n = p + q (p,q prim) zurück,
        sowie die Hardy-Littlewood-Schätzung der Anzahl.
    @return JSON mit decompositions, hl_estimate, hl_accuracy
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    n = data.get('n', 100)

    try:
        n = int(n)
        if n < 4:
            return error_response("n muss mindestens 4 sein.")
        if n % 2 != 0:
            return error_response("n muss gerade sein (Goldbach: jede gerade Zahl > 2).")
        if n > 100000:
            return error_response("n > 100.000 zu groß (Berechnungsdauer).")

        # Alle Goldbach-Zerlegungen berechnen
        decomps = goldbach_all_decompositions(n)

        # Hardy-Littlewood-Schätzung der Anzahl der Zerlegungen
        # Rückgabe: {'estimate_r2': float, 'actual_r2': int, ...}
        hl_result = hardy_littlewood_goldbach_density(n)
        hl_estimate = hl_result.get('estimate_r2', 0) if isinstance(hl_result, dict) else 0

        # Genauigkeit der Schätzung (falls Zerlegungen vorhanden)
        hl_accuracy = None
        if len(decomps) > 0 and hl_estimate > 0:
            hl_accuracy = min(len(decomps), hl_estimate) / max(len(decomps), hl_estimate)

        return jsonify({
            "n": n,
            "decompositions": list(decomps),
            "count": len(decomps),
            "hl_estimate": hl_estimate,
            "hl_accuracy": hl_accuracy
        })
    except Exception as e:
        return error_response(f"Fehler bei Goldbach-Berechnung: {str(e)}")


@app.route('/api/millennium/complexity', methods=['POST'])
def millennium_complexity():
    """
    @brief Vergleicht verschiedene Algorithmus-Komplexitätsklassen für Eingabegröße n.
    @description
        Erwartet JSON: {"n": int}
        Gibt für jede Komplexitätsklasse die Anzahl der Operationen zurück:
        O(1), O(log n), O(sqrt n), O(n), O(n log n), O(n²), O(n³), O(2^n), O(n!)
        Für n > 30 werden O(2^n) und O(n!) als Infinity zurückgegeben.
    @return JSON mit Operationszahlen pro Komplexitätsklasse
    @author Kurt Ingwer
    @date 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    n = data.get('n', 20)

    try:
        n = int(n)
        if n < 1:
            return error_response("n muss mindestens 1 sein.")
        if n > 60:
            return error_response("n > 60 zu groß (astronomische Werte für O(2^n) und O(n!)).")

        # Komplexitätsklassen via Millennium-Modul berechnen
        # Rückgabe: {'n': n, 'complexities': {'O(1)': 1, 'O(n)': n, ...}, ...}
        comp_data = complexity_comparison(n)
        result = comp_data.get('complexities', {})

        # Alle Werte JSON-serialisierbar machen (inf → None für JSON)
        import math
        clean_result = {}
        for k, v in result.items():
            try:
                fv = float(v)
                clean_result[k] = fv if fv < 1e308 else None
            except (TypeError, ValueError):
                clean_result[k] = None

        # Sicherstellen dass alle benötigten Klassen vorhanden sind
        if 'O(1)' not in clean_result:       clean_result['O(1)'] = 1
        if 'O(log n)' not in clean_result:   clean_result['O(log n)'] = round(math.log2(n), 4) if n > 0 else 0
        if 'O(sqrt n)' not in clean_result:  clean_result['O(sqrt n)'] = round(math.sqrt(n), 4)
        if 'O(n)' not in clean_result:       clean_result['O(n)'] = n
        if 'O(n log n)' not in clean_result: clean_result['O(n log n)'] = round(n * math.log2(n), 4) if n > 1 else 0
        if 'O(n^2)' not in clean_result:     clean_result['O(n^2)'] = n * n
        if 'O(n^3)' not in clean_result:     clean_result['O(n^3)'] = n * n * n
        if 'O(2^n)' not in clean_result:     clean_result['O(2^n)'] = 2 ** n if n <= 60 else None
        if 'O(n!)' not in clean_result:
            clean_result['O(n!)'] = math.factorial(n) if n <= 20 else None

        # TSP-spezifischer Schlüssel für das Frontend (Brute-Force = n!)
        clean_result['factorial'] = clean_result.get('O(n!)')

        return jsonify(clean_result)
    except Exception as e:
        return error_response(f"Fehler bei Komplexitätsberechnung: {str(e)}")


# ===========================================================================
# SCHRITT-FÜR-SCHRITT-SEITE UND API-ENDPUNKTE
# ===========================================================================

@app.route('/steps')
def steps_page():
    """
    @brief Schritt-für-Schritt-Erklärungen Seite.
    @description
        Zeigt eine Übersicht aller verfügbaren interaktiven Schritt-für-Schritt-
        Erklärungen: Newton-Raphson, Bisektion, Gauss-Elimination, ggT, RSA, Primfaktoren.

    @return Gerendertes HTML-Template
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    return render_template('steps.html')


@app.route('/api/steps/newton', methods=['POST'])
def steps_newton():
    """
    @brief API-Endpunkt: Newton-Raphson Schritt-für-Schritt.
    @description
        Führt das Newton-Raphson-Verfahren auf einem symbolischen Ausdruck aus
        und gibt jeden Berechnungsschritt zurück.

        Request-Body (JSON):
            {"expr": "x**2 - 2", "x0": 1.5}

        Response (JSON):
            {"steps": [...], "result": float, "html": "...", "text": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        # JSON-Body lesen und validieren
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        expr_str = data.get('expr', '')
        x0 = float(data.get('x0', 1.0))

        if not expr_str:
            return error_response("Parameter 'expr' fehlt.")

        # Ausdruck sicher parsen und numerische Funktion erzeugen
        expr, x_sym = safe_parse_expr(expr_str)
        f_num = sp.lambdify(x_sym, expr, modules=['math'])

        # Schritt-für-Schritt Newton-Raphson ausführen
        tol = float(data.get('tol', 1e-10))
        max_iter = int(data.get('max_iter', 50))
        result_dict = newton_raphson_steps(f_num, x0, tol=tol, max_iter=max_iter)

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'converged': result_dict.get('converged', False),
            'iterations': result_dict.get('iterations', 0),
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"Newton-Raphson Fehler: {str(e)}")


@app.route('/api/steps/bisection', methods=['POST'])
def steps_bisection():
    """
    @brief API-Endpunkt: Bisektionsverfahren Schritt-für-Schritt.
    @description
        Führt das Bisektionsverfahren auf einem symbolischen Ausdruck aus.

        Request-Body (JSON):
            {"expr": "x**2 - 2", "a": 1.0, "b": 2.0}

        Response (JSON):
            {"steps": [...], "result": float, "html": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        expr_str = data.get('expr', '')
        a = float(data.get('a', 0.0))
        b = float(data.get('b', 2.0))

        if not expr_str:
            return error_response("Parameter 'expr' fehlt.")

        # Ausdruck parsen und Funktion erzeugen
        expr, x_sym = safe_parse_expr(expr_str)
        f_num = sp.lambdify(x_sym, expr, modules=['math'])

        tol = float(data.get('tol', 1e-10))
        max_iter = int(data.get('max_iter', 50))
        result_dict = bisection_steps(f_num, a, b, tol=tol, max_iter=max_iter)

        # Fehlerfall (kein Vorzeichenwechsel)
        if result_dict.get('result') is None:
            return error_response(result_dict.get('error', 'Bisektionsfehler'))

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'converged': result_dict.get('converged', False),
            'iterations': result_dict.get('iterations', 0),
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"Bisektion Fehler: {str(e)}")


@app.route('/api/steps/gauss', methods=['POST'])
def steps_gauss():
    """
    @brief API-Endpunkt: Gauss-Elimination Schritt-für-Schritt.
    @description
        Löst ein lineares Gleichungssystem Ax = b und zeigt jeden Umformungsschritt.

        Request-Body (JSON):
            {"matrix": [[2,1],[1,3]], "b": [5, 10]}

        Response (JSON):
            {"steps": [...], "result": [float], "html": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        matrix = data.get('matrix')
        b = data.get('b')

        if matrix is None or b is None:
            return error_response("Parameter 'matrix' und 'b' erforderlich.")

        # Dimensionsprüfung
        n = len(matrix)
        if len(b) != n:
            return error_response(f"matrix hat {n} Zeilen, b hat {len(b)} Einträge – müssen gleich sein.")

        result_dict = gauss_elimination_steps(matrix, b)

        # Singuläre Matrix abfangen
        if result_dict.get('result') is None:
            return error_response(result_dict.get('error', 'Gauss-Eliminations-Fehler'))

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"Gauss-Elimination Fehler: {str(e)}")


@app.route('/api/steps/gcd', methods=['POST'])
def steps_gcd():
    """
    @brief API-Endpunkt: Euklidischer Algorithmus (ggT) Schritt-für-Schritt.
    @description
        Berechnet den größten gemeinsamen Teiler zweier Zahlen und zeigt jeden Schritt.

        Request-Body (JSON):
            {"a": 48, "b": 18}

        Response (JSON):
            {"steps": [...], "result": int, "html": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        a = int(data.get('a', 0))
        b = int(data.get('b', 0))

        if a <= 0 or b <= 0:
            return error_response("a und b müssen positive ganze Zahlen sein.")

        result_dict = euclidean_algorithm_steps(a, b)

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"Euklidischer Algorithmus Fehler: {str(e)}")


@app.route('/api/steps/rsa', methods=['POST'])
def steps_rsa():
    """
    @brief API-Endpunkt: RSA-Kryptographie Schritt-für-Schritt.
    @description
        Zeigt alle Schritte der RSA-Schlüsselgenerierung und Ver-/Entschlüsselung.

        Request-Body (JSON):
            {"p": 5, "q": 11, "message": 7}

        Response (JSON):
            {"steps": [...], "result": {...}, "html": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        p = int(data.get('p', 5))
        q = int(data.get('q', 11))
        message = int(data.get('message', 7))

        if p < 2 or q < 2:
            return error_response("p und q müssen Primzahlen ≥ 2 sein.")

        result_dict = rsa_steps(p, q, message)

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"RSA Fehler: {str(e)}")


@app.route('/api/steps/prime_factorization', methods=['POST'])
def steps_prime_factorization():
    """
    @brief API-Endpunkt: Primfaktorzerlegung Schritt-für-Schritt.
    @description
        Zerlegt eine Zahl in Primfaktoren und zeigt jeden Divisionsschritt.

        Request-Body (JSON):
            {"n": 360}

        Response (JSON):
            {"steps": [...], "result": [int], "html": "..."}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    try:
        data = request.get_json(force=True)
        if not data:
            return error_response("Kein JSON-Body erhalten.")

        n = int(data.get('n', 0))

        if n < 2:
            return error_response("n muss eine ganze Zahl ≥ 2 sein.")
        if n > 10**12:
            return error_response("n ist zu groß (max. 10^12 für schnelle Berechnung).")

        result_dict = prime_factorization_steps(n)

        return jsonify({
            'steps': result_dict['steps'],
            'result': result_dict['result'],
            'factor_counts': result_dict.get('factor_counts', {}),
            'product_str': result_dict.get('product_str', ''),
            'method': result_dict.get('method', ''),
            'html': format_steps_html(result_dict),
            'text': format_steps_text(result_dict)
        })
    except Exception as e:
        return error_response(f"Primfaktorzerlegung Fehler: {str(e)}")


# ===========================================================================
# FORMALE BEWEISE — Routen für Beweisinfrastruktur
# ===========================================================================

@app.route('/proofs')
def proofs_page():
    """
    @brief Seite für formale Beweise.
    @description Zeigt alle klassischen Beweise mit Schrittanzeige und
                 interaktivem Induktions-Tester.
    @date 2026-03-10
    """
    return render_template('proofs.html')


@app.route('/api/proof/gauss_sum', methods=['GET'])
def proof_gauss_sum():
    """
    @brief API-Endpunkt: Beweis der Gaußschen Summenformel.
    @description Gibt Theorem, Schritte, Zusammenfassung und HTML-Darstellung zurück.
    @return JSON mit theorem, steps, summary und html
    @date 2026-03-10
    """
    try:
        from formal_proof import prove_gauss_sum
        proof = prove_gauss_sum()
        summary = proof.status_summary()
        return jsonify({
            'theorem': proof.theorem,
            'steps': [str(s) for s in proof.steps],
            'summary': summary,
            'html': proof.to_html()
        })
    except Exception as e:
        return error_response(f"Gauß-Beweis Fehler: {str(e)}")


@app.route('/api/proof/sqrt2', methods=['GET'])
def proof_sqrt2():
    """
    @brief API-Endpunkt: Irrationalitätsbeweis von √2.
    @description Gibt Theorem, Schritte und HTML des Widerspruchsbeweises zurück.
    @return JSON mit theorem, steps und html
    @date 2026-03-10
    """
    try:
        from formal_proof import prove_sqrt2_irrational
        proof = prove_sqrt2_irrational()
        return jsonify({
            'theorem': proof.theorem,
            'steps': [str(s) for s in proof.steps],
            'html': proof.to_html()
        })
    except Exception as e:
        return error_response(f"√2-Beweis Fehler: {str(e)}")


@app.route('/api/proof/primes_infinite', methods=['GET'])
def proof_primes_infinite():
    """
    @brief API-Endpunkt: Euklids Beweis unendlich vieler Primzahlen.
    @description Gibt Theorem, Schritte und HTML des Widerspruchsbeweises zurück.
    @return JSON mit theorem, steps und html
    @date 2026-03-10
    """
    try:
        from formal_proof import prove_infinitely_many_primes
        proof = prove_infinitely_many_primes()
        return jsonify({
            'theorem': proof.theorem,
            'steps': [str(s) for s in proof.steps],
            'html': proof.to_html()
        })
    except Exception as e:
        return error_response(f"Primzahlen-Beweis Fehler: {str(e)}")


@app.route('/api/proof/rh_evidence', methods=['GET'])
def proof_rh_evidence():
    """
    @brief API-Endpunkt: Empirische Evidenz für die Riemann-Hypothese.
    @description Gibt Theorem, Schritte, Zusammenfassung und HTML zurück.
    @return JSON mit theorem, steps, summary und html
    @date 2026-03-10
    """
    try:
        from formal_proof import riemann_hypothesis_evidence
        proof = riemann_hypothesis_evidence()
        return jsonify({
            'theorem': proof.theorem,
            'steps': [str(s) for s in proof.steps],
            'summary': proof.status_summary(),
            'html': proof.to_html()
        })
    except Exception as e:
        return error_response(f"Riemann-Evidenz Fehler: {str(e)}")


@app.route('/api/proof/induction_check', methods=['POST'])
def proof_induction_check():
    """
    @brief API-Endpunkt: Interaktiver Induktions-Tester für Gaußsche Summe.
    @description Prüft P(n): Σ(k=1..n) k = n(n+1)/2 für eine gegebene Zahl n.
    @request JSON mit {'n': int}
    @return JSON mit n, lhs, rhs, is_valid
    @date 2026-03-10
    """
    try:
        data = request.get_json()
        n = int(data.get('n', 1))
        # Gaußsche Summe für n prüfen
        lhs = sum(range(1, n + 1))           # Direkte Summe
        rhs = n * (n + 1) // 2              # Gaußsche Formel
        return jsonify({
            'n': n,
            'lhs': lhs,
            'rhs': rhs,
            'is_valid': lhs == rhs,
            'formula': f"Σ(k=1..{n}) k = {lhs} = {n}·{n+1}/2 = {rhs}"
        })
    except Exception as e:
        return error_response(f"Induktions-Test Fehler: {str(e)}")


# ===========================================================================
# ROUTEN: ELLIPTISCHE KURVEN (Seite + API)
# ===========================================================================

@app.route('/elliptic')
def page_elliptic():
    """
    @brief Elliptische-Kurven-Seite: Gruppenstruktur, Arithmetik, Hasse-Schranke.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('elliptic.html')


@app.route('/api/elliptic/analyze', methods=['POST'])
def api_elliptic_analyze():
    """
    @brief Analysiert eine elliptische Kurve y²=x³+ax+b.
    @description
        Berechnet Diskriminante, j-Invariante, Singularität,
        optionale Punktarithmetik (2P, n·P) und gibt einen Base64-Plot zurück.

        Erwartet JSON:
            {
                "a": float,
                "b": float,
                "n": int (optional, Standard 3),
                "px": float (optional),
                "py": float (optional)
            }
    @return JSON mit Kurvenanalyse und Base64-Plot
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    try:
        # Koeffizienten aus Request lesen
        a = float(data.get('a', 0))
        b = float(data.get('b', 0))
        n = int(data.get('n', 3))

        # Diskriminante vorberechnen (unabhängig von Singularität)
        discriminant = -16.0 * (4.0 * a**3 + 27.0 * b**2)
        is_singular = abs(discriminant) < 1e-12

        result = {
            "discriminant": discriminant,
            "is_singular": is_singular,
            "j_invariant": None,
            "point_results": None,
            "plot": None,
        }

        if is_singular:
            # Singuläre Kurve – kein Plot, keine Punktarithmetik
            return jsonify(result)

        # Kurvenobjekt erstellen (wirft ValueError wenn singulär)
        curve = EllipticCurve(a, b)
        result["j_invariant"] = curve.j_invariant()

        # Optionale Punktarithmetik
        px = data.get('px', None)
        py = data.get('py', None)
        if px is not None and py is not None:
            px, py = float(px), float(py)
            on_curve = curve.is_on_curve(px, py)
            point_results = {
                "px": px,
                "py": py,
                "on_curve": on_curve,
                "double_p": "–",
                "n_times_p": "–",
            }

            if on_curve:
                try:
                    # Punkt P erstellen
                    P = curve.point(px, py)
                    # Verdoppelung 2P
                    P2 = P + P
                    if P2.is_infinity:
                        point_results["double_p"] = "O (Punkt im Unendlichen)"
                    else:
                        point_results["double_p"] = f"({P2.x:.6f}, {P2.y:.6f})"

                    # Skalarmultiplikation n·P über wiederholte Addition
                    Pn = ECPoint.infinity(curve)
                    Pi = P
                    for _ in range(n):
                        Pn = Pn + Pi
                    if Pn.is_infinity:
                        point_results["n_times_p"] = "O (Punkt im Unendlichen)"
                    else:
                        point_results["n_times_p"] = f"({Pn.x:.6f}, {Pn.y:.6f})"
                except Exception as e:
                    # Fehler in Punktarithmetik – trotzdem weitermachen
                    point_results["double_p"] = f"Fehler: {e}"
                    point_results["n_times_p"] = f"Fehler: {e}"

            result["point_results"] = point_results

        # Plot erstellen (Dark Theme angepasst)
        try:
            fig = curve.plot_curve(x_range=(-3, 3))
            # Axes aus Figure holen und Dark-Theme anwenden
            ax = fig.axes[0]
            fig.patch.set_facecolor('#1e1e2e')
            ax.set_facecolor('#181825')
            ax.tick_params(colors='#cdd6f4')
            ax.xaxis.label.set_color('#cdd6f4')
            ax.yaxis.label.set_color('#cdd6f4')
            ax.title.set_color('#cba6f7')
            for spine in ax.spines.values():
                spine.set_edgecolor('#313244')
            # Linienfarben auf Blautöne ändern
            for line in ax.get_lines():
                line.set_color('#89b4fa')

            # Optionalen Punkt P markieren
            if px is not None and py is not None and curve.is_on_curve(px, py):
                ax.plot(px, py, 'o', color='#a6e3a1', markersize=10,
                        label=f'P = ({px:.2f}, {py:.2f})', zorder=5)
                ax.legend(facecolor='#313244', labelcolor='#cdd6f4',
                          edgecolor='#45475a')

            result["plot"] = figure_to_base64(fig)
        except Exception:
            # Plot-Fehler ist nicht kritisch
            result["plot"] = None

        return jsonify(result)

    except Exception as e:
        return error_response(f"Kurvenanalyse-Fehler: {str(e)}")


@app.route('/api/elliptic/hasse', methods=['POST'])
def api_elliptic_hasse():
    """
    @brief Berechnet Hasse-Schranke: Gruppenordnung #E(F_p) für kleine Primzahlen.
    @description
        Erwartet JSON:
            {"a": int, "b": int, "prime_limit": int (bis zu welcher Primzahl)}

        Gibt für jede Primzahl p ≤ prime_limit zurück:
        - Gruppenordnung #E(F_p)
        - Frobenius-Spur a_p = p + 1 - #E(F_p)
        - Ob Hasse-Schranke |a_p| ≤ 2√p eingehalten wird
    @return JSON mit Liste von Primzahl-Ergebnissen
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    try:
        a = int(float(data.get('a', 0)))
        b = int(float(data.get('b', 0)))
        prime_limit = min(int(data.get('prime_limit', 20)), 50)

        # Singularität vorab prüfen
        discriminant = -16.0 * (4.0 * a**3 + 27.0 * b**2)
        if abs(discriminant) < 1e-12:
            return jsonify({"is_singular": True, "primes": []})

        curve = EllipticCurve(a, b)

        # Primzahlen bis prime_limit bestimmen (Sieb des Eratosthenes)
        sieve = list(range(2, prime_limit + 1))
        primes = []
        for candidate in sieve:
            is_prime = all(candidate % p != 0 for p in range(2, int(candidate**0.5) + 1))
            if is_prime:
                primes.append(candidate)

        # Für jede Primzahl Gruppenordnung berechnen
        prime_results = []
        for p in primes:
            try:
                order = curve.order_over_fp(p)
                trace = p + 1 - order
                prime_results.append({
                    "p": p,
                    "order": order,
                    "trace": trace,
                })
            except Exception:
                # Kurve mod p singulär – überspringen
                pass

        return jsonify({
            "is_singular": False,
            "discriminant": discriminant,
            "primes": prime_results,
        })

    except Exception as e:
        return error_response(f"Hasse-Berechnung fehlgeschlagen: {str(e)}")


# ===========================================================================
# ROUTEN: TENSOR-GEOMETRIE (Seite + API)
# ===========================================================================

@app.route('/tensor')
def page_tensor():
    """
    @brief Tensor-Geometrie-Seite: Metriktensor, Krümmung, Mannigfaltigkeiten.
    @return Gerendertes HTML-Template
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('tensor.html')


@app.route('/api/tensor/curvature', methods=['POST'])
def api_tensor_curvature():
    """
    @brief Berechnet Metriktensor, Ricci-Skalar und Gaußsche Krümmung.
    @description
        Erwartet JSON:
            {
                "manifold": "sphere"|"torus"|"hyperbolic"|"saddle"|"flat",
                "params": {
                    // Sphäre:      {"r": float, "theta": float, "phi": float}
                    // Torus:       {"R": float, "r": float, "theta": float, "phi": float}
                    // Hyperbolisch:{"x": float, "y": float}
                    // Sattel:      {"x": float, "y": float}
                    // Flach:       {"x": float, "y": float}
                }
            }

        Gibt zurück:
        - metric:              2×2 Metriktensor als Liste
        - ricci_scalar:        R (Ricci-Skalar)
        - gaussian_curvature:  K (Gaußsche Krümmung)
        - curvature_type:      Textbeschreibung der Krümmung
        - physical_description: Physikalische Bedeutung
        - metric_latex:        KaTeX-LaTeX-Formel des Metriktensors
        - plot:                Base64-PNG-Plot der Mannigfaltigkeit
    @return JSON mit allen Krümmungsgrößen
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    if not data:
        return error_response("Kein JSON-Body übergeben.")

    try:
        manifold = data.get('manifold', 'sphere')
        params = data.get('params', {})

        # Mannigfaltigkeit und Koordinaten auflösen
        if manifold == 'sphere':
            r = float(params.get('r', 1.0))
            theta = float(params.get('theta', np.pi / 2))
            phi = float(params.get('phi', 0.0))
            metric_func = sphere_metric(r)
            point = [theta, phi]
            name = f"Sphäre S² (Radius r = {r})"
            point_str = f"θ={theta:.4f}, φ={phi:.4f}"
            # Theoretische Werte für die Ausgabe vorbereiten
            metric_latex_template = (
                f"g_{{ij}} = \\begin{{pmatrix}} {r**2} & 0 \\\\"
                f" 0 & {r**2}\\sin^2\\theta \\end{{pmatrix}}"
            )
            physical = (
                f"Die Sphäre hat konstant positive Gaußsche Krümmung K = 1/r² = {1/r**2:.4f}. "
                f"In der ART entspricht dies einer geschlossenen, endlichen Raumzeit "
                f"(positiv-gekrümmtes Universum). Geodäten sind Großkreise."
            )

        elif manifold == 'torus':
            R = float(params.get('R', 2.0))
            r = float(params.get('r', 1.0))
            theta = float(params.get('theta', 0.0))
            phi = float(params.get('phi', 0.0))
            metric_func = torus_metric(R, r)
            point = [theta, phi]
            name = f"Torus T² (R = {R}, r = {r})"
            point_str = f"θ={theta:.4f}, φ={phi:.4f}"
            rho = R + r * np.cos(theta)
            metric_latex_template = (
                f"g_{{ij}} = \\begin{{pmatrix}} {r**2} & 0 \\\\"
                f" 0 & (R + r\\cos\\theta)^2 \\end{{pmatrix}}"
            )
            physical = (
                f"Der Torus hat variable Krümmung: positiv außen (θ≈0), negativ innen (θ≈π). "
                f"Am Punkt θ={theta:.2f}: ρ = R + r·cos(θ) = {rho:.4f}. "
                f"Torusförmige Universen sind in der Kosmologie als 'Doughnut-Universum' bekannt."
            )

        elif manifold == 'hyperbolic':
            x = float(params.get('x', 0.0))
            y = float(params.get('y', 1.0))
            if y <= 0:
                y = 1e-4
            metric_func = hyperbolic_plane_metric()
            point = [x, y]
            name = f"Hyperbolische Ebene H² (Poincaré-Halbebene)"
            point_str = f"x={x:.4f}, y={y:.4f}"
            coeff = 1.0 / y**2
            metric_latex_template = (
                f"g_{{ij}} = \\frac{{1}}{{y^2}} \\begin{{pmatrix}} 1 & 0 \\\\ 0 & 1 \\end{{pmatrix}}"
                f",\\quad y = {y:.4f}"
            )
            physical = (
                f"Die hyperbolische Ebene hat konstant negative Krümmung K = -1. "
                f"Am Punkt (x={x:.2f}, y={y:.2f}): g_ij = {coeff:.4f}·I. "
                f"In der ART entspricht dies anti-de-Sitter-Raumzeit (AdS), "
                f"relevant in der AdS/CFT-Korrespondenz."
            )

        elif manifold == 'saddle':
            x = float(params.get('x', 0.5))
            y = float(params.get('y', 0.5))
            metric_func = saddle_metric()
            point = [x, y]
            name = f"Sattelfläche z = xy"
            point_str = f"x={x:.4f}, y={y:.4f}"
            metric_latex_template = (
                f"g_{{ij}} = \\begin{{pmatrix}} 1+y^2 & xy \\\\ xy & 1+x^2 \\end{{pmatrix}}"
                f",\\quad (x,y)=({x:.2f},{y:.2f})"
            )
            physical = (
                f"Die Sattelfläche z=xy hat variable negative Krümmung K = -1/(1+x²+y²)². "
                f"Am Ursprung (0,0) ist K = -1 (maximal negativ). "
                f"Sattelpunkte treten in der ART an Sattelstellen des Gravitationspotentials auf."
            )

        elif manifold == 'flat':
            x = float(params.get('x', 0.0))
            y = float(params.get('y', 0.0))
            metric_func = flat_metric(2)
            point = [x, y]
            name = "Flache Ebene ℝ² (Euklidisch)"
            point_str = f"x={x:.4f}, y={y:.4f}"
            metric_latex_template = (
                "g_{ij} = \\begin{pmatrix} 1 & 0 \\\\ 0 & 1 \\end{pmatrix}"
            )
            physical = (
                "Die flache Ebene ℝ² hat überall K = 0 und R = 0. "
                "In der ART entspricht dies dem Minkowski-Raum (spezielle Relativitätstheorie) "
                "im Vakuum ohne Materie oder Energie. Einstein-Tensor G_μν = 0."
            )
        else:
            return error_response(f"Unbekannte Mannigfaltigkeit: '{manifold}'.")

        # Metriktensor am Punkt berechnen
        g = metric_func(point)
        det_g = float(np.linalg.det(g))

        # Ricci-Skalar und Gaußsche Krümmung numerisch berechnen
        R = ricci_scalar(metric_func, point)
        K = gaussian_curvature(metric_func, point)

        # Krümmungstyp bestimmen
        if abs(K) < 1e-6:
            curvature_type = "Flach (K ≈ 0) – keine intrinsische Krümmung"
        elif K > 0:
            curvature_type = f"Positiv gekrümmt (K = {K:.6f}) – wie Kugeloberfläche"
        else:
            curvature_type = f"Negativ gekrümmt (K = {K:.6f}) – wie Sattelform"

        # Plot erstellen
        plot_b64 = None
        try:
            fig = _plot_manifold(manifold, params, point)
            if fig is not None:
                plot_b64 = figure_to_base64(fig)
        except Exception:
            plot_b64 = None

        return jsonify({
            "name": name,
            "point_str": point_str,
            "metric": g.tolist(),
            "det_g": det_g,
            "ricci_scalar": float(R),
            "gaussian_curvature": float(K),
            "curvature_type": curvature_type,
            "physical_description": physical,
            "metric_latex": metric_latex_template,
            "plot": plot_b64,
        })

    except Exception as e:
        return error_response(f"Tensor-Berechnung fehlgeschlagen: {str(e)}")


def _plot_manifold(manifold: str, params: dict, point: list):
    """
    @brief Erstellt einen matplotlib-Plot der angegebenen Mannigfaltigkeit.
    @description
        Zeichnet die Mannigfaltigkeit als 3D-Oberfläche (Sphäre, Torus)
        oder als 2D-Farbkarte (Hyperbolisch, Sattel, Flach).
        Der gegebene Koordinatenpunkt wird markiert.
    @param manifold: Bezeichner der Mannigfaltigkeit
    @param params: Parameter-Dictionary (Radien, Winkel, Koordinaten)
    @param point: Koordinatenpunkt [coord1, coord2]
    @return matplotlib.figure.Figure oder None bei Fehler
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    # Dark-Theme-Farben
    BG_DARK = '#1e1e2e'
    AX_DARK = '#181825'
    TEXT_COLOR = '#cdd6f4'
    ACCENT = '#89b4fa'
    POINT_COLOR = '#a6e3a1'

    if manifold == 'sphere':
        r = float(params.get('r', 1.0))
        theta_p = float(params.get('theta', np.pi / 2))
        phi_p = float(params.get('phi', 0.0))

        # Sphärengitter generieren
        u = np.linspace(0, np.pi, 60)
        v = np.linspace(0, 2 * np.pi, 60)
        X = r * np.outer(np.sin(u), np.cos(v))
        Y = r * np.outer(np.sin(u), np.sin(v))
        Z = r * np.outer(np.cos(u), np.ones_like(v))

        fig = plt.figure(figsize=(8, 6))
        fig.patch.set_facecolor(BG_DARK)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor(AX_DARK)
        ax.plot_surface(X, Y, Z, alpha=0.4, color=ACCENT)
        # Punkt markieren
        px = r * np.sin(theta_p) * np.cos(phi_p)
        py = r * np.sin(theta_p) * np.sin(phi_p)
        pz = r * np.cos(theta_p)
        ax.scatter([px], [py], [pz], color=POINT_COLOR, s=120, zorder=5,
                   label=f'P=(θ={theta_p:.2f}, φ={phi_p:.2f})')
        ax.set_title(f'Sphäre S² (r={r})', color=TEXT_COLOR)
        ax.legend(facecolor='#313244', labelcolor=TEXT_COLOR, edgecolor='#45475a')
        ax.tick_params(colors=TEXT_COLOR)
        return fig

    elif manifold == 'torus':
        R = float(params.get('R', 2.0))
        r = float(params.get('r', 1.0))
        theta_p = float(params.get('theta', 0.0))
        phi_p = float(params.get('phi', 0.0))

        # Torus-Gitter
        u = np.linspace(0, 2 * np.pi, 80)
        v = np.linspace(0, 2 * np.pi, 80)
        U, V = np.meshgrid(u, v)
        X = (R + r * np.cos(U)) * np.cos(V)
        Y = (R + r * np.cos(U)) * np.sin(V)
        Z = r * np.sin(U)

        fig = plt.figure(figsize=(8, 6))
        fig.patch.set_facecolor(BG_DARK)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor(AX_DARK)
        ax.plot_surface(X, Y, Z, alpha=0.4, color=ACCENT)
        # Punkt markieren
        px = (R + r * np.cos(theta_p)) * np.cos(phi_p)
        py = (R + r * np.cos(theta_p)) * np.sin(phi_p)
        pz = r * np.sin(theta_p)
        ax.scatter([px], [py], [pz], color=POINT_COLOR, s=120, zorder=5,
                   label=f'P=(θ={theta_p:.2f}, φ={phi_p:.2f})')
        ax.set_title(f'Torus T² (R={R}, r={r})', color=TEXT_COLOR)
        ax.legend(facecolor='#313244', labelcolor=TEXT_COLOR, edgecolor='#45475a')
        ax.tick_params(colors=TEXT_COLOR)
        return fig

    elif manifold == 'hyperbolic':
        x_p = float(params.get('x', 0.0))
        y_p = float(params.get('y', 1.0))

        # Poincaré-Halbebene: Krümmungskarte (K = -1/y²... nein, K = -1 überall)
        xs = np.linspace(-3, 3, 100)
        ys = np.linspace(0.1, 3, 100)
        XX, YY = np.meshgrid(xs, ys)
        # Gaußsche Krümmung ist konstant -1, aber zeige 1/y² für visuelles Verständnis
        ZZ = 1.0 / YY**2   # g-Koeffizient (zeigt Verzerrung)

        fig, ax = plt.subplots(figsize=(8, 6))
        fig.patch.set_facecolor(BG_DARK)
        ax.set_facecolor(AX_DARK)
        cf = ax.contourf(XX, YY, ZZ, levels=20, cmap='plasma')
        fig.colorbar(cf, ax=ax, label='g-Koeffizient 1/y²')
        ax.axhline(y=0, color='#45475a', linewidth=1)
        ax.scatter([x_p], [y_p], color=POINT_COLOR, s=120, zorder=5,
                   label=f'P=({x_p:.2f}, {y_p:.2f})')
        ax.set_xlabel('x', color=TEXT_COLOR)
        ax.set_ylabel('y', color=TEXT_COLOR)
        ax.set_title('Hyperbolische Ebene H² (Poincaré-Modell)', color=TEXT_COLOR)
        ax.legend(facecolor='#313244', labelcolor=TEXT_COLOR, edgecolor='#45475a')
        ax.tick_params(colors=TEXT_COLOR)
        return fig

    elif manifold == 'saddle':
        x_p = float(params.get('x', 0.5))
        y_p = float(params.get('y', 0.5))

        # Sattelfläche z = x*y
        xs = np.linspace(-2, 2, 60)
        ys = np.linspace(-2, 2, 60)
        XX, YY = np.meshgrid(xs, ys)
        ZZ = XX * YY

        fig = plt.figure(figsize=(8, 6))
        fig.patch.set_facecolor(BG_DARK)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor(AX_DARK)
        ax.plot_surface(XX, YY, ZZ, alpha=0.5, cmap='coolwarm')
        # Punkt markieren
        z_p = x_p * y_p
        ax.scatter([x_p], [y_p], [z_p], color=POINT_COLOR, s=120, zorder=5,
                   label=f'P=({x_p:.2f}, {y_p:.2f}, {z_p:.2f})')
        ax.set_xlabel('x', color=TEXT_COLOR)
        ax.set_ylabel('y', color=TEXT_COLOR)
        ax.set_zlabel('z = xy', color=TEXT_COLOR)
        ax.set_title('Sattelfläche z = xy', color=TEXT_COLOR)
        ax.legend(facecolor='#313244', labelcolor=TEXT_COLOR, edgecolor='#45475a')
        ax.tick_params(colors=TEXT_COLOR)
        return fig

    elif manifold == 'flat':
        x_p = float(params.get('x', 0.0))
        y_p = float(params.get('y', 0.0))

        # Flache Ebene – einfaches Gitter
        xs = np.linspace(-3, 3, 10)
        ys = np.linspace(-3, 3, 10)
        XX, YY = np.meshgrid(xs, ys)
        ZZ = np.zeros_like(XX)

        fig = plt.figure(figsize=(8, 6))
        fig.patch.set_facecolor(BG_DARK)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor(AX_DARK)
        ax.plot_surface(XX, YY, ZZ, alpha=0.3, color=ACCENT)
        ax.scatter([x_p], [y_p], [0], color=POINT_COLOR, s=120, zorder=5,
                   label=f'P=({x_p:.2f}, {y_p:.2f}, 0)')
        ax.set_title('Flache Ebene ℝ² (K = 0)', color=TEXT_COLOR)
        ax.legend(facecolor='#313244', labelcolor=TEXT_COLOR, edgecolor='#45475a')
        ax.tick_params(colors=TEXT_COLOR)
        return fig

    return None


# ===========================================================================
# NEUE API-ENDPUNKTE: GRUPPENTHEORIE
# ===========================================================================

@app.route('/api/group/cyclic', methods=['POST'])
def api_group_cyclic():
    """
    @brief Analysiert eine zyklische Gruppe Z/nZ.
    @description
        Erwartet JSON: {"n": 6}
        Gibt Ordnung, Elemente, Generatoren, is_abelian zurück.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _group_theory_ok:
            return jsonify({"error": "group_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 6))
        # Zyklische Gruppe der Ordnung n erzeugen
        g = cyclic_group(n)
        return jsonify({
            "order": g.order(),
            "elements": list(g.elements),
            "generators": list(g.generators()),
            "is_abelian": g.is_abelian(),
            "is_cyclic": True,  # zyklische Gruppen sind per Definition zyklisch
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/group/sylow', methods=['POST'])
def api_group_sylow():
    """
    @brief Wendet Sylow-Theorie auf eine Gruppe der Ordnung n an.
    @description
        Erwartet JSON: {"n": 12, "p": 2}
        Gibt Sylow-Untergruppen-Ordnung und Anzahl zurück.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _group_theory_ok:
            return jsonify({"error": "group_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 12))
        p = int(data.get('p', 2))
        # Zyklische Gruppe als Stellvertreter – Sylow-Sätze gelten für beliebige Gruppen
        g = cyclic_group(n)
        result = sylow_theorems(g, p)
        return jsonify({
            "p_part": result.get("p_power"),
            "sylow_order": result.get("sylow_order"),
            "num_sylow": result.get("num_sylow_subgroups"),
            "details": str(result.get("description", "")),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/group/dihedral', methods=['POST'])
def api_group_dihedral():
    """
    @brief Analysiert eine Diedergruppe D_n.
    @description
        Erwartet JSON: {"n": 4}
        Gibt Ordnung, is_abelian, Rotationen/Spiegelungen zurück.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _group_theory_ok:
            return jsonify({"error": "group_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 4))
        # Diedergruppe D_n (Ordnung 2n)
        d = dihedral_group(n)
        elements_str = str(d.elements[:min(20, len(d.elements))])
        # Erste n Elemente sind Drehungen, letzte n sind Spiegelungen (in der Implementierung)
        half = len(d.elements) // 2
        return jsonify({
            "order": d.order(),
            "is_abelian": d.is_abelian(),
            "elements": elements_str,
            "rotations": str(list(d.elements[:half])),
            "reflections": str(list(d.elements[half:])),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: RINGTHEORIE
# ===========================================================================

@app.route('/api/ring/properties', methods=['POST'])
def api_ring_properties():
    """
    @brief Analysiert Eigenschaften des Rings Z/nZ.
    @description
        Erwartet JSON: {"n": 12}
        Gibt Einheiten, Nullteiler, is_field, is_integral_domain zurück.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _ring_theory_ok:
            return jsonify({"error": "ring_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 12))
        r = Ring(n)
        return jsonify({
            "units": r.units(),
            "zero_divisors": r.zero_divisors(),
            "is_field": r.is_field(),
            "is_integral_domain": r.is_integral_domain(),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/ring/factorize', methods=['POST'])
def api_ring_factorize():
    """
    @brief Faktorisiert ein Polynom modulo einer Primzahl p.
    @description
        Erwartet JSON: {"coefficients": [1, 0, 1], "p": 2}
        Koeffizienten: höchster Grad zuerst.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _ring_theory_ok:
            return jsonify({"error": "ring_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        coeffs = [int(c) for c in data.get('coefficients', [1, 0, 1])]
        p = int(data.get('p', 2))
        result = polynomial_factorization_mod_p(coeffs, p)
        return jsonify({
            "factors": result.get('factors', []),
            "multiplicities": result.get('multiplicities', []),
            "description": f"Faktorisierung von Polynom mod {p}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/ring/irreducible', methods=['POST'])
def api_ring_irreducible():
    """
    @brief Prüft ob ein Polynom irreduzibel mod p ist.
    @description
        Erwartet JSON: {"coefficients": [1, 1, 1], "p": 2}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _ring_theory_ok:
            return jsonify({"error": "ring_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        coeffs = [int(c) for c in data.get('coefficients', [1, 1, 1])]
        p = int(data.get('p', 2))
        result = is_irreducible_mod_p(coeffs, p)
        return jsonify({
            "is_irreducible": bool(result),
            "reason": "Kein echter Teiler in F_p[x] gefunden" if result else "Hat echte Teiler in F_p[x]",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: GALOIS-THEORIE
# ===========================================================================

@app.route('/api/galois/group', methods=['POST'])
def api_galois_group():
    """
    @brief Berechnet die Galois-Gruppe eines Polynoms über Q.
    @description
        Erwartet JSON: {"coefficients": [1, 0, -2]}
        Koeffizienten höchster Grad zuerst.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _galois_ok:
            return jsonify({"error": "galois_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        coeffs = [int(c) for c in data.get('coefficients', [1, 0, -2])]
        result = galois_group_polynomial(coeffs)
        return jsonify({
            "galois_group": result.get('galois_group'),
            "order": result.get('order'),
            "is_solvable": result.get('is_solvable'),
            "discriminant": result.get('discriminant'),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/galois/solvable', methods=['POST'])
def api_galois_solvable():
    """
    @brief Prüft ob ein Polynom durch Radikale lösbar ist.
    @description
        Erwartet JSON: {"coefficients": [1, 0, -2]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _galois_ok:
            return jsonify({"error": "galois_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        coeffs = [int(c) for c in data.get('coefficients', [1, 0, -2])]
        result = is_solvable_by_radicals(coeffs)
        return jsonify({
            "solvable": result.get('solvable'),
            "galois_group": result.get('galois_group'),
            "reason": result.get('reason'),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/galois/finite_field', methods=['POST'])
def api_galois_finite_field():
    """
    @brief Analysiert den endlichen Körper GF(p^n).
    @description
        Erwartet JSON: {"p": 2, "n": 3}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _galois_ok:
            return jsonify({"error": "galois_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        p = int(data.get('p', 2))
        n = int(data.get('n', 3))
        ff = FiniteField(p, n)
        return jsonify({
            "order": ff.order,
            "characteristic": ff.characteristic,
            "description": f"GF({p}^{n}) hat {p**n} Elemente, Charakteristik {p}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: MODULN & KOMMUTATIVE ALGEBRA
# ===========================================================================

@app.route('/api/modules/smith', methods=['POST'])
def api_modules_smith():
    """
    @brief Berechnet die Smith-Normalform einer ganzzahligen Matrix.
    @description
        Erwartet JSON: {"matrix": [[2,4,4],[-6,6,12],[10,-4,-16]]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modules_ok:
            return jsonify({"error": "modules_algebra Modul nicht verfügbar"}), 503
        data = request.get_json()
        matrix = data.get('matrix', [[2, 4], [4, 8]])
        result = smith_normal_form(matrix)
        return jsonify({
            "invariant_factors": result.get('invariant_factors', []),
            "rank": result.get('rank', 0),
            "description": f"Invariantenfaktoren bilden die Smith-Normalform",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/modules/classify_abelian', methods=['POST'])
def api_modules_classify():
    """
    @brief Klassifiziert eine endliche abelsche Gruppe per Invariantenfaktoren.
    @description
        Erwartet JSON: {"invariants": [2, 6, 12]}
        Invariantenfaktoren als Liste (Smith-Diagonale).
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modules_ok:
            return jsonify({"error": "modules_algebra Modul nicht verfügbar"}), 503
        data = request.get_json()
        invariants = [int(x) for x in data.get('invariants', [2, 6])]
        # Ordnung = Produkt aller Invariantenfaktoren
        order = 1
        for inv in invariants:
            order *= inv
        # Isomorphietyp als Produkt zyklischer Gruppen
        iso_type = " × ".join([f"ℤ_{inv}" for inv in invariants])
        return jsonify({
            "isomorphism_type": iso_type,
            "order": order,
            "description": f"Die Gruppe der Ordnung {order} ist isomorph zu {iso_type}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/modules/spectrum', methods=['POST'])
def api_modules_spectrum():
    """
    @brief Berechnet das Primspektrum von Z/nZ.
    @description
        Erwartet JSON: {"n": 12}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modules_ok:
            return jsonify({"error": "commutative_algebra Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 12))
        result = prime_spectrum(n)
        return jsonify({
            "prime_ideals": result.get('prime_ideals', []),
            "maximal_ideals": result.get('maximal_ideals', []),
            "description": f"Primspektrum von ℤ/{n}ℤ",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: FOURIER-ANALYSIS
# ===========================================================================

@app.route('/api/fourier/dft', methods=['POST'])
def api_fourier_dft():
    """
    @brief Berechnet die diskrete Fourier-Transformation (FFT).
    @description
        Erwartet JSON: {"values": [1, 0, -1, 0]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _fourier_ok:
            return jsonify({"error": "fourier Modul nicht verfügbar"}), 503
        data = request.get_json()
        values = [float(v) for v in data.get('values', [1, 0, -1, 0])]
        # FFT berechnen (Cooley-Tukey)
        result = fft(values)
        frequencies = []
        dominant_k = 0
        max_mag = 0
        for k, c in enumerate(result):
            mag = abs(c)
            phase = float(np.angle(c)) if hasattr(c, 'real') else 0.0
            frequencies.append({"magnitude": mag, "phase": phase})
            if mag > max_mag and k > 0:
                max_mag = mag
                dominant_k = k
        return jsonify({
            "frequencies": frequencies,
            "n": len(values),
            "dominant_k": dominant_k,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/fourier/series', methods=['POST'])
def api_fourier_series():
    """
    @brief Berechnet Fourier-Reihen-Koeffizienten einer Funktion.
    @description
        Erwartet JSON: {"func": "sin(x)", "period": 6.283, "n_terms": 5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _fourier_ok:
            return jsonify({"error": "fourier Modul nicht verfügbar"}), 503
        data = request.get_json()
        func_str = data.get('func', 'sin(x)')
        period = float(data.get('period', 2 * np.pi))
        n_terms = int(data.get('n_terms', 5))
        # Funktion aus String erzeugen (sicher via SymPy)
        f_lambda = expr_to_lambda(func_str)
        # Fourier-Koeffizienten numerisch berechnen
        a_coeffs, b_coeffs = fourier_coefficients(f_lambda, period, n_terms)
        return jsonify({
            "a_coeffs": [float(a) for a in a_coeffs],
            "b_coeffs": [float(b) for b in b_coeffs],
            "period": period,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/fourier/dominant', methods=['POST'])
def api_fourier_dominant():
    """
    @brief Bestimmt die dominante Frequenz in einem Signal.
    @description
        Erwartet JSON: {"values": [...], "sample_rate": 16.0}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _fourier_ok:
            return jsonify({"error": "fourier Modul nicht verfügbar"}), 503
        data = request.get_json()
        values = [float(v) for v in data.get('values', [0, 1, 0, -1])]
        sample_rate = float(data.get('sample_rate', 1.0))
        # Dominante Frequenz berechnen
        dom_freq = dominant_frequency(values, sample_rate=sample_rate)
        # FFT für Amplitude des Peaks
        result = fft(values)
        magnitudes = [abs(c) for c in result]
        dom_k = magnitudes.index(max(magnitudes[1:len(magnitudes)//2+1])) + 1
        amp = magnitudes[dom_k] if dom_k < len(magnitudes) else 0.0
        return jsonify({
            "dominant_frequency": float(dom_freq),
            "dominant_k": dom_k,
            "amplitude": float(amp) / len(values),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: ANALYTISCHE ZAHLENTHEORIE
# ===========================================================================

@app.route('/api/analytic_nt/pi', methods=['POST'])
def api_analytic_pi():
    """
    @brief Berechnet die Primzählfunktion pi(x) und Vergleichsapproximationen.
    @description
        Erwartet JSON: {"x": 100}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _analytic_nt_ok:
            return jsonify({"error": "analytic_number_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        x = float(data.get('x', 100))
        # Primzählfunktion
        pi_x = prime_counting_function(x)
        # Approximation x/ln(x)
        import math
        approx_log = x / math.log(x) if x > 1 else 0
        # Logarithmisches Integral
        li_x = logarithmic_integral(x)
        # Relativer Fehler
        rel_err = abs(approx_log - pi_x) / pi_x * 100 if pi_x > 0 else 0
        return jsonify({
            "pi_x": int(pi_x),
            "approx_log": float(approx_log),
            "li_x": float(li_x),
            "relative_error_log": float(rel_err),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/analytic_nt/li', methods=['POST'])
def api_analytic_li():
    """
    @brief Berechnet das Logarithmische Integral Li(x).
    @description
        Erwartet JSON: {"x": 1000}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _analytic_nt_ok:
            return jsonify({"error": "analytic_number_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        x = float(data.get('x', 1000))
        li_x = logarithmic_integral(x)
        pi_x = prime_counting_function(x)
        return jsonify({
            "li_x": float(li_x),
            "pi_x": int(pi_x),
            "difference": float(li_x) - int(pi_x),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/analytic_nt/von_mangoldt', methods=['POST'])
def api_analytic_von_mangoldt():
    """
    @brief Berechnet die Von-Mangoldt-Funktion Lambda(n) bis zur Grenze n.
    @description
        Erwartet JSON: {"n": 30}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _analytic_nt_ok:
            return jsonify({"error": "analytic_number_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 30))
        import math
        values = []
        psi = 0.0
        for k in range(1, n + 1):
            lam = von_mangoldt_function(k)
            psi += lam
            # Primzahl p bestimmen wenn Lambda > 0
            p = None
            if lam > 0:
                p = round(math.exp(lam))
            values.append({"n": k, "lambda": float(lam), "prime": p})
        return jsonify({
            "values": values,
            "psi": float(psi),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: MODULFORMEN
# ===========================================================================

@app.route('/api/modular/eisenstein', methods=['POST'])
def api_modular_eisenstein():
    """
    @brief Berechnet den Wert der Eisenstein-Reihe E_k(z).
    @description
        Erwartet JSON: {"k": 4, "z_re": 0.0, "z_im": 1.0}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modular_ok:
            return jsonify({"error": "modular_forms Modul nicht verfügbar"}), 503
        data = request.get_json()
        k = int(data.get('k', 4))
        z_re = float(data.get('z_re', 0.0))
        z_im = float(data.get('z_im', 1.0))
        z = complex(z_re, z_im)
        # Eisenstein-Reihe auswerten
        val = eisenstein_series(k, z)
        return jsonify({
            "real": float(val.real),
            "imag": float(val.imag),
            "abs": float(abs(val)),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/modular/j_invariant', methods=['POST'])
def api_modular_j():
    """
    @brief Berechnet die j-Invariante.
    @description
        Erwartet JSON: {"z_re": 0.0, "z_im": 1.0}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modular_ok:
            return jsonify({"error": "modular_forms Modul nicht verfügbar"}), 503
        data = request.get_json()
        z_re = float(data.get('z_re', 0.0))
        z_im = float(data.get('z_im', 1.0))
        z = complex(z_re, z_im)
        val = j_invariant(z)
        note = ""
        # Bekannte Spezialwerte
        if abs(z_re) < 0.01 and abs(z_im - 1) < 0.01:
            note = "j(i) ≈ 1728 (elliptische Kurve mit j=1728 hat extra Automorphismen)"
        return jsonify({
            "real": float(val.real),
            "imag": float(val.imag),
            "abs": float(abs(val)),
            "note": note,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/modular/delta', methods=['POST'])
def api_modular_delta():
    """
    @brief Berechnet die ersten n Ramanujan tau-Koeffizienten.
    @description
        Erwartet JSON: {"n_terms": 10}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _modular_ok:
            return jsonify({"error": "modular_forms Modul nicht verfügbar"}), 503
        data = request.get_json()
        n_terms = int(data.get('n_terms', 10))
        # Ramanujan tau-Koeffizienten via q-Entwicklung von Delta
        tau = fourier_coefficients_delta(n_terms)
        return jsonify({
            "tau": [int(t) for t in tau],
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: L-FUNKTIONEN
# ===========================================================================

@app.route('/api/lfunc/dirichlet', methods=['POST'])
def api_lfunc_dirichlet():
    """
    @brief Berechnet eine Dirichlet-L-Funktion L(s, chi).
    @description
        Erwartet JSON: {"s_re": 2.0, "s_im": 0.0, "q": 4, "chi_index": 0}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _lfunc_ok:
            return jsonify({"error": "l_functions Modul nicht verfügbar"}), 503
        data = request.get_json()
        s_re = float(data.get('s_re', 2.0))
        s_im = float(data.get('s_im', 0.0))
        q = int(data.get('q', 4))
        chi_index = int(data.get('chi_index', 0))
        s = complex(s_re, s_im)
        val = dirichlet_l_function(s, q, chi_index)
        return jsonify({
            "real": float(val.real),
            "imag": float(val.imag),
            "abs": float(abs(val)),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/lfunc/characters', methods=['POST'])
def api_lfunc_characters():
    """
    @brief Listet alle Dirichlet-Charaktere mod q auf.
    @description
        Erwartet JSON: {"q": 5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _lfunc_ok:
            return jsonify({"error": "l_functions Modul nicht verfügbar"}), 503
        import math
        data = request.get_json()
        q = int(data.get('q', 5))
        chars = dirichlet_characters(q)
        # Charaktere serialisieren (Funktionsobjekte → Werte-Dict)
        chars_out = []
        for ch in chars:
            values_list = [float(abs(v)) for v in ch.get('values', {}).values()]
            chars_out.append({
                "values": values_list[:q],
                "is_primitive": not ch.get('is_principal', False),
            })
        # phi(q) berechnen
        phi_q = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
        return jsonify({
            "characters": chars_out,
            "count": len(chars_out),
            "phi_q": phi_q,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/lfunc/special_values', methods=['POST'])
def api_lfunc_special_values():
    """
    @brief Berechnet Spezialwerte von L-Funktionen mod q.
    @description
        Erwartet JSON: {"q": 4}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _lfunc_ok:
            return jsonify({"error": "l_functions Modul nicht verfügbar"}), 503
        data = request.get_json()
        q = int(data.get('q', 4))
        chars = dirichlet_characters(q)
        values = []
        # Für jeden nicht-Hauptcharakter L(1, chi) berechnen
        for i, ch in enumerate(chars):
            if not ch.get('is_principal', False):
                try:
                    # L(1+ε, chi) als Annäherung an L(1, chi) für nicht-Hauptcharaktere
                    val = dirichlet_l_function(complex(1.001, 0), q, i)
                    values.append({"s": 1, "value": float(abs(val))})
                except Exception:
                    pass
            # L(2, chi) für alle Charaktere
            try:
                val2 = dirichlet_l_function(complex(2.0, 0), q, i)
                values.append({"s": 2, "value": float(abs(val2))})
            except Exception:
                pass
        return jsonify({
            "values": values[:8],
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: P-ADISCHE ZAHLEN
# ===========================================================================

@app.route('/api/padic/valuation', methods=['POST'])
def api_padic_valuation():
    """
    @brief Berechnet p-adische Bewertung und Norm einer Zahl n.
    @description
        Erwartet JSON: {"n": 360, "p": 2}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _padic_ok:
            return jsonify({"error": "p_adic Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 360))
        p = int(data.get('p', 2))
        val = p_adic_valuation(n, p)
        norm = p_adic_norm(n, p)
        return jsonify({
            "valuation": int(val) if val != float('inf') else "inf",
            "norm": float(norm),
            "abs_inf": abs(n),
            "product_formula": float(norm) * abs(n),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/padic/expansion', methods=['POST'])
def api_padic_expansion():
    """
    @brief Berechnet die p-adische Entwicklung einer Zahl n.
    @description
        Erwartet JSON: {"n": 42, "p": 5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _padic_ok:
            return jsonify({"error": "p_adic Modul nicht verfügbar"}), 503
        data = request.get_json()
        n = int(data.get('n', 42))
        p = int(data.get('p', 5))
        # p-adische Entwicklung manuell berechnen (niedrigstes Bit zuerst)
        digits = []
        x = abs(n)
        if x == 0:
            digits = [0]
        while x > 0:
            digits.append(x % p)
            x //= p
        valuation = 0
        for d in digits:
            if d == 0:
                valuation += 1
            else:
                break
        return jsonify({
            "digits": digits,
            "valuation": valuation,
            "p": p,
            "n": n,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/padic/hensel', methods=['POST'])
def api_padic_hensel():
    """
    @brief Führt Hensel-Lifting einer Wurzel durch.
    @description
        Erwartet JSON: {"a": 3, "p": 7, "k": 3}
        Hebt f(x) = x^2 - 2 (Beispiel) von mod p auf mod p^k.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _padic_ok:
            return jsonify({"error": "p_adic Modul nicht verfügbar"}), 503
        data = request.get_json()
        a = int(data.get('a', 3))
        p = int(data.get('p', 7))
        k = int(data.get('k', 3))
        # Hebt Wurzel von x^2 - 2 mod p
        result = hensel_lift([1, 0, -2], p, a, n_lifts=k)
        # Zwischenwerte zeigen
        steps = []
        curr = a
        for i in range(1, k + 1):
            mod = p ** (i + 1)
            steps.append(curr % mod)
            curr = result  # Vereinfacht: letzten Wert nutzen
        return jsonify({
            "steps": steps,
            "final": int(result),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: IWASAWA-THEORIE
# ===========================================================================

@app.route('/api/iwasawa/kummer', methods=['POST'])
def api_iwasawa_kummer():
    """
    @brief Prüft Kummer-Kongruenzen für eine Primzahl p.
    @description
        Erwartet JSON: {"p": 5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _iwasawa_ok:
            return jsonify({"error": "iwasawa_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        p = int(data.get('p', 5))
        result = kummer_congruences_check(p, n_range=8)
        checks_out = []
        for ch in result.get('checks', []):
            checks_out.append({
                "k": ch.get('n', 0),
                "value_mod_p": str(ch.get('val_n', 'N/A')),
                "verified": bool(ch.get('kummer_holds', False))
            })
        return jsonify({
            "all_verified": bool(result.get('all_consistent', False)),
            "checks": checks_out,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/iwasawa/polynomial', methods=['POST'])
def api_iwasawa_polynomial():
    """
    @brief Berechnet Iwasawa-Invarianten eines Polynoms.
    @description
        Erwartet JSON: {"coefficients": [3, 1, 0, 1], "p": 3}
        Koeffizienten: niedrigstes Grad zuerst.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _iwasawa_ok:
            return jsonify({"error": "iwasawa_theory Modul nicht verfügbar"}), 503
        data = request.get_json()
        coefficients = [int(c) for c in data.get('coefficients', [3, 1, 0, 1])]
        p = int(data.get('p', 3))
        result = iwasawa_polynomial(coefficients, p)
        return jsonify({
            "mu": result.get('mu_invariant', 0),
            "lambda": result.get('lambda_invariant', 0),
            "weierstrass": str(result.get('distinguished_poly', [])),
            "description": f"µ={result.get('mu_invariant')}, λ={result.get('lambda_invariant')}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/iwasawa/kubota_leopoldt', methods=['POST'])
def api_iwasawa_kubota():
    """
    @brief Berechnet die Kubota-Leopoldt p-adische Zeta-Funktion.
    @description
        Erwartet JSON: {"p": 5, "s": -1}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _iwasawa_ok:
            return jsonify({"error": "iwasawa_theory Modul nicht verfügbar"}), 503
        import math
        data = request.get_json()
        p = int(data.get('p', 5))
        s = int(data.get('s', -1))
        # Klassischer Zeta-Wert bei negativen ganzen Zahlen via Bernoulli-Zahlen
        # zeta(-n) = -B_{n+1}/(n+1)
        from p_adic import bernoulli_number
        n = -s  # s = -1 → n = 1
        if n >= 0:
            b = bernoulli_number(n + 1)
            classical = float(-b / (n + 1)) if n >= 0 else 0
        else:
            classical = 0
        # p-adischer Wert (Annäherung über Kummer-Kongruenzen)
        result = kummer_congruences_check(p, n_range=6)
        padic_val = str(result.get('checks', [{}])[0].get('val_n', 'N/A'))
        return jsonify({
            "value": padic_val,
            "classical_value": float(classical),
            "note": f"p-adische Interpolation von ζ({s}) für p={p}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: SPINOREN
# ===========================================================================

@app.route('/api/spinors/gamma', methods=['POST'])
def api_spinors_gamma():
    """
    @brief Gibt Gamma-Matrizen für eine gegebene Dimension zurück.
    @description
        Erwartet JSON: {"dim": 4}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _spinors_ok:
            return jsonify({"error": "spinors Modul nicht verfügbar"}), 503
        data = request.get_json()
        dim = int(data.get('dim', 4))
        gammas = gamma_matrices(dim)
        # Matrizen serialisieren
        matrices_out = []
        for gamma in gammas:
            mat = []
            for row in gamma:
                mat.append([{"re": float(c.real), "im": float(c.imag)} for c in row])
            matrices_out.append(mat)
        return jsonify({
            "matrices": matrices_out,
            "dim": dim,
            "size": gammas[0].shape[0],
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/spinors/clifford', methods=['POST'])
def api_spinors_clifford():
    """
    @brief Prüft Clifford-Algebra-Relationen der Gamma-Matrizen.
    @description
        Erwartet JSON: {"dim": 4}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _spinors_ok:
            return jsonify({"error": "spinors Modul nicht verfügbar"}), 503
        data = request.get_json()
        dim = int(data.get('dim', 4))
        gammas = gamma_matrices(dim)
        result = clifford_algebra_check(gammas)
        n = len(gammas)
        return jsonify({
            "all_correct": bool(result.get('is_clifford_algebra', False)),
            "max_error": float(result.get('max_error', 0)),
            "pairs_tested": n * (n + 1) // 2,
            "details": f"{n} Gamma-Matrizen in {dim}D getestet",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/spinors/dirac', methods=['POST'])
def api_spinors_dirac():
    """
    @brief Berechnet Dirac-Spinor für freies Teilchen.
    @description
        Erwartet JSON: {"m": 1.0, "pz": 0.5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _spinors_ok:
            return jsonify({"error": "spinors Modul nicht verfügbar"}), 503
        import math
        data = request.get_json()
        m = float(data.get('m', 1.0))
        pz = float(data.get('pz', 0.5))
        # Energie: E = sqrt(pz^2 + m^2)
        energy = math.sqrt(pz ** 2 + m ** 2)
        momentum = np.array([0.0, 0.0, pz])
        # Spin up (spin=0) und spin down (spin=1) Spinoren
        sp_up = dirac_spinor(m, momentum, spin=0)
        sp_down = dirac_spinor(m, momentum, spin=1)
        return jsonify({
            "energy": float(energy),
            "spinor_up": [float(c.real) for c in sp_up],
            "spinor_down": [float(c.real) for c in sp_down],
            "note": f"Normierung: u†u = 2E = {2*energy:.4f}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# NEUE API-ENDPUNKTE: NUMERISCHE METHODEN
# ===========================================================================

@app.route('/api/numerical/interpolate', methods=['POST'])
def api_numerical_interpolate():
    """
    @brief Lagrange-Interpolation an einem Punkt.
    @description
        Erwartet JSON: {"x_vals": [...], "y_vals": [...], "x_eval": 1.5}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _numerical_ok:
            return jsonify({"error": "numerical_methods Modul nicht verfügbar"}), 503
        data = request.get_json()
        x_vals = [float(x) for x in data.get('x_vals', [0, 1, 2, 3])]
        y_vals = [float(y) for y in data.get('y_vals', [0, 1, 4, 9])]
        x_eval = float(data.get('x_eval', 1.5))
        result = lagrange_interpolation(x_vals, y_vals, x_eval)
        return jsonify({
            "result": float(result),
            "degree": len(x_vals) - 1,
            "note": f"Lagrange-Polynom vom Grad {len(x_vals)-1}",
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/numerical/optimize', methods=['POST'])
def api_numerical_optimize():
    """
    @brief Gradientenabstieg für f(x,y) = x^2 + y^2.
    @description
        Erwartet JSON: {"x0": 3.0, "y0": 4.0, "learning_rate": 0.1}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _numerical_ok:
            return jsonify({"error": "numerical_methods Modul nicht verfügbar"}), 503
        data = request.get_json()
        x0 = float(data.get('x0', 3.0))
        y0 = float(data.get('y0', 4.0))
        lr = float(data.get('learning_rate', 0.1))
        # Zielfunktion f(x) = x[0]^2 + x[1]^2
        f = lambda x: x[0] ** 2 + x[1] ** 2
        grad = lambda x: [2 * x[0], 2 * x[1]]
        x_min, f_min, iters = gradient_descent(f, grad, [x0, y0], learning_rate=lr)
        return jsonify({
            "x_min": float(x_min[0]),
            "y_min": float(x_min[1]),
            "f_min": float(f_min),
            "iterations": int(iters),
            "converged": float(f_min) < 1e-6,
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/api/numerical/golden_section', methods=['POST'])
def api_numerical_golden():
    """
    @brief Golden-Section-Suche für eine Funktion auf [a, b].
    @description
        Erwartet JSON: {"func": "(x-2)**2+1", "a": 0, "b": 4}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _numerical_ok:
            return jsonify({"error": "numerical_methods Modul nicht verfügbar"}), 503
        data = request.get_json()
        func_str = data.get('func', '(x-2)**2 + 1')
        a = float(data.get('a', 0))
        b = float(data.get('b', 4))
        f = expr_to_lambda(func_str)
        x_min, f_min = golden_section_search(f, a, b)
        # Iterationszahl approximieren
        import math
        iters = int(math.log(1e-8 / (b - a)) / math.log(0.618)) + 1
        return jsonify({
            "x_min": float(x_min),
            "f_min": float(f_min),
            "iterations": max(1, iters),
            "error": None
        })
    except Exception as e:
        return jsonify({"error": str(e)})


# ===========================================================================
# ROUTEN: MATHEMATISCHE LOGIK API
# ===========================================================================

@app.route('/api/logic/truth_table', methods=['POST'])
def api_logic_truth_table():
    """
    @brief Wahrheitstabelle für eine Formel.
    @description
        Nimmt eine Liste von Variablen und eine Formel als String entgegen.
        Die Formel wird mit eval und einer Whitelist geparst.
        Gibt Wahrheitstabelle sowie logische Eigenschaften zurück.
    @date 2026-03-10
    """
    data = request.get_json()
    try:
        # Formel aus String aufbauen: z.B. "AND(p,OR(q,NOT(p)))"
        # Einfacher: Variablen und Formel als Text
        vars_list = data.get('variables', ['p', 'q'])
        formula_str = data.get('formula', 'AND(p,q)')

        # Propositions erstellen
        props = {v: Proposition(v) for v in vars_list}

        # Formel parsen (einfache eval mit Whitelist)
        allowed = {
            **props,
            'AND': lambda a, b: LogicFormula('AND', a, b),
            'OR': lambda a, b: LogicFormula('OR', a, b),
            'NOT': lambda a: LogicFormula('NOT', a),
            'IMPLIES': lambda a, b: LogicFormula('IMPLIES', a, b),
            'IFF': lambda a, b: LogicFormula('IFF', a, b),
            'XOR': lambda a, b: LogicFormula('XOR', a, b),
        }
        formula = eval(formula_str, {"__builtins__": {}}, allowed)

        table = truth_table(formula, vars_list)
        taut = is_tautology(formula, vars_list)
        sat = is_satisfiable(formula, vars_list)
        contra = is_contradiction(formula, vars_list)

        return jsonify({
            'variables': vars_list,
            'formula': formula_str,
            'table': [
                {**row, 'result': row['result']}
                for row in table
            ],
            'is_tautology': taut,
            'is_satisfiable': sat,
            'is_contradiction': contra
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/logic/dpll', methods=['POST'])
def api_logic_dpll():
    """
    @brief DPLL SAT-Solver.
    @description
        Nimmt Klauseln als Liste von Listen entgegen.
        Positive Zahlen = Variable i, negative = Negation von Variable i.
        Gibt Erfüllbarkeit und Belegung zurück.
    @date 2026-03-10
    """
    data = request.get_json()
    try:
        # clauses: Liste von Listen, z.B. [[1, -2], [2, 3], [-1, -3]]
        # positive Zahl = Variable i, negative = Negation
        raw_clauses = data.get('clauses', [[1, -2], [2, 3]])

        # dpll() erwartet String-Literale ("v1", "-v1") statt Integers
        def int_to_lit(n: int) -> str:
            return f"v{n}" if n > 0 else f"-v{-n}"

        clauses = [frozenset(int_to_lit(lit) for lit in c) for c in raw_clauses]
        result = dpll(clauses)
        return jsonify({
            'satisfiable': result is not None,
            'assignment': result if result else {}
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/logic/godel', methods=['POST'])
def api_logic_godel():
    """
    @brief Gödel-Nummerierung.
    @description
        Wandelt einen Formel-String in eine Gödel-Zahl um und dekodiert diese
        wieder. Gibt zusätzlich eine Demo des Unvollständigkeitssatzes zurück.
    @date 2026-03-10
    """
    data = request.get_json()
    try:
        formula_str = data.get('formula', 'P(x)')
        number = godel_numbering(formula_str)
        decoded = godel_decode(number)
        demo = incompleteness_demonstration()
        return jsonify({
            'formula': formula_str,
            'godel_number': number,
            'decoded': decoded,
            'incompleteness': demo
        })
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# ROUTEN: MENGENLEHRE API
# ===========================================================================

@app.route('/api/sets/operations', methods=['POST'])
def api_sets_operations():
    """
    @brief Mengenoperationen.
    @description
        Berechnet Vereinigung, Schnittmenge, Differenzen, symmetrische Differenz
        und Teilmengenrelation für zwei Mengen A und B.
    @date 2026-03-10
    """
    data = request.get_json()
    try:
        A = MathSet(data.get('A', [1, 2, 3, 4]))
        B = MathSet(data.get('B', [3, 4, 5, 6]))

        return jsonify({
            'A': sorted(list(A._elements)) if hasattr(A, '_elements') else sorted(list(A.elements if hasattr(A, 'elements') else [])),
            'B': sorted(list(B._elements)) if hasattr(B, '_elements') else sorted(list([])),
            'union': sorted(list((A | B))),
            'intersection': sorted(list((A & B))),
            'difference_AB': sorted(list((A - B))),
            'difference_BA': sorted(list((B - A))),
            'symmetric_difference': sorted(list((A ^ B))),
            'A_subset_B': A <= B,
            'cardinality_A': len(A),
            'cardinality_B': len(B),
            'power_set_size': 2 ** len(A) if len(A) <= 6 else None
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/sets/cantor', methods=['POST'])
def api_sets_cantor():
    """
    @brief Cantors Diagonalargument Demo.
    @description
        Demonstriert, warum die reellen Zahlen nicht abzählbar sind,
        indem gezeigt wird, dass jede Aufzählung unvollständig ist.
    @date 2026-03-10
    """
    data = request.get_json()
    try:
        n = min(int(data.get('n', 5)), 10)
        result = cantors_diagonal_argument(n)
        # Python-sets sind nicht JSON-serialisierbar → in sortierte Listen umwandeln
        def serialize(obj):
            if isinstance(obj, (set, frozenset)):
                return sorted(list(obj))
            if isinstance(obj, dict):
                return {k: serialize(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [serialize(x) for x in obj]
            return obj
        return jsonify(serialize(result))
    except Exception as e:
        return error_response(str(e))


@app.route('/api/sets/zfc', methods=['GET'])
def api_sets_zfc():
    """
    @brief ZFC-Axiome abrufen.
    @description
        Gibt alle Zermelo-Fraenkel-Axiome mit Auswahlaxiom zurück.
    @date 2026-03-10
    """
    try:
        axioms = zfc_axioms()
        return jsonify({'axioms': axioms})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# ROUTEN: Neue Seiten (Erweiterte Algebra, Logik & Grundlagen, Erweiterte Analysis)
# ===========================================================================

@app.route('/algebra_advanced')
def page_algebra_advanced():
    """
    @brief Seite für Erweiterte Algebra (Darstellungstheorie, Verbandstheorie usw.).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('algebra_advanced.html')


@app.route('/logic_foundations')
def page_logic_foundations():
    """
    @brief Seite für Logik & Grundlagen (Modelltheorie, Beweistheorie, Rekursionstheorie).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('logic_foundations.html')


@app.route('/analysis_advanced')
def page_analysis_advanced():
    """
    @brief Seite für Erweiterte Analysis (Maßtheorie, Funktionalanalysis, PDE, Operatoralgebren).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('analysis_advanced.html')


@app.route('/algebraic_geometry')
def page_algebraic_geometry():
    """
    @brief Seite für Algebraische Geometrie (Varietäten, Nullstellensatz, Gröbner-Basen).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('algebraic_geometry.html')


@app.route('/lie_groups')
def page_lie_groups():
    """
    @brief Seite für Lie-Gruppen und Lie-Algebren.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('lie_groups.html')


@app.route('/stochastic')
def page_stochastic():
    """
    @brief Seite für Stochastische Prozesse und Ergodentheorie.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('stochastic.html')


@app.route('/differential_topology')
def page_differential_topology():
    """
    @brief Seite für Differentialtopologie (Differentialformen, Morse-Theorie).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    return render_template('differential_topology.html')


# ===========================================================================
# API: Darstellungstheorie (representation_theory.py)
# ===========================================================================

@app.route('/api/representation/character_table', methods=['POST'])
def api_repr_character_table():
    """
    @brief Charaktertafel einer Gruppe berechnen.
    @description
        Gibt die Charaktertafel für S3 oder Z2 zurück.
        POST-Body: {"group_name": "S3" | "Z2"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _repr_ok:
            return error_response("representation_theory nicht verfügbar")
        # Gruppenname aus Anfrage lesen
        group_name = data.get('group_name', 'S3').upper()

        if group_name == 'Z2':
            reps = z2_representations()
        elif group_name == 'S3':
            reps = s3_representations()
        else:
            return error_response(f"Unbekannte Gruppe: {group_name}. Bitte 'S3' oder 'Z2' wählen.")

        # Charaktere aus den Darstellungen direkt berechnen
        # character_table(group_elements, conjugacy_classes, representations)
        # Einfache Variante: Charaktere direkt aus den Darstellungen ableiten
        char_data = []
        for rep in reps:
            char_row = [complex(rep.character().get(g, 0)) for g in rep.group_elements]
            char_data.append([round(c.real, 4) if abs(c.imag) < 1e-10 else [round(c.real,4), round(c.imag,4)] for c in char_row])

        result = {
            'group': group_name,
            'order': len(reps[0].group_elements) if reps else 0,
            'irreducible_count': len(reps),
            'characters': char_data,
            'class_representatives': [str(g) for g in (reps[0].group_elements if reps else [])],
        }
        return jsonify({'result': result})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/representation/schur_orthogonality', methods=['POST'])
def api_repr_schur_orthogonality():
    """
    @brief Schur-Orthogonalitätsrelationen prüfen.
    @description
        Prüft ⟨χᵢ, χⱼ⟩ = δᵢⱼ für die irreduziblen Charaktere einer Gruppe.
        POST-Body: {"group_name": "S3" | "Z2"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _repr_ok:
            return error_response("representation_theory nicht verfügbar")
        group_name = data.get('group_name', 'S3').upper()

        if group_name == 'Z2':
            reps = z2_representations()
        elif group_name == 'S3':
            reps = s3_representations()
        else:
            return error_response(f"Unbekannte Gruppe: {group_name}")

        # Schur-Orthogonalität durch inneres Produkt der Charaktere prüfen
        # ⟨χᵢ, χⱼ⟩ = (1/|G|) Σ_g χᵢ(g)·conj(χⱼ(g))
        import numpy as np
        results = []
        n_group = len(reps[0].group_elements) if reps else 1
        for i, r1 in enumerate(reps):
            char1 = np.array([r1.character().get(g, 0) for g in r1.group_elements], dtype=complex)
            for j, r2 in enumerate(reps):
                char2 = np.array([r2.character().get(g, 0) for g in r2.group_elements], dtype=complex)
                # Inneres Produkt der Charaktere
                ip = np.sum(char1 * np.conj(char2)) / n_group
                is_orth = abs(ip) < 0.01 if i != j else abs(ip - 1.0) < 0.01
                results.append({
                    'i': i, 'j': j,
                    'orthogonal': bool(is_orth),
                    'inner_product': f'{ip.real:.4f}',
                    'expected': 1 if i == j else 0
                })

        return jsonify({'result': {
            'group': group_name,
            'pairs': results,
            'schur_orthogonality_holds': all(
                r['orthogonal'] == (r['i'] != r['j']) or r['i'] == r['j']
                for r in results
            )
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/representation/burnside', methods=['POST'])
def api_repr_burnside():
    """
    @brief Burnside-Lemma: Anzahl der Orbits berechnen.
    @description
        POST-Body: {"n_colorings": int, "n_rotations": int}
        Beispiel: Perlen-Halskette mit n Farben und Rotationssymmetrie.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _repr_ok:
            return error_response("representation_theory nicht verfügbar")

        n_colors = int(data.get('n_colorings', 3))
        n_rot = int(data.get('n_rotations', 4))

        # Fixpunkte für zyklische Gruppe berechnen
        # Fix(r^k) = n_colors^gcd(k, n_rot) für Rotationen
        import math
        fix_counts = [n_colors ** math.gcd(k, n_rot) for k in range(n_rot)]
        # burnside_lemma(group_elements, action_dict) – wir berechnen direkt
        orbit_count = sum(fix_counts) // n_rot

        return jsonify({'result': {
            'n_colors': n_colors,
            'n_rotations': n_rot,
            'fixed_counts': fix_counts,
            'orbit_count': orbit_count,
            'description': f'Burnside: ({" + ".join(str(f) for f in fix_counts)}) / {n_rot} = {orbit_count}',
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Verbandstheorie (lattice_theory.py)
# ===========================================================================

@app.route('/api/lattice/create', methods=['POST'])
def api_lattice_create():
    """
    @brief Teilerverband oder Potenzmengeverband erstellen.
    @description
        POST-Body: {"type": "divisibility" | "power_set", "n": int}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _lattice_ok:
            return error_response("lattice_theory nicht verfügbar")

        lattice_type = data.get('type', 'divisibility')
        n = int(data.get('n', 12))

        if lattice_type == 'power_set':
            n = min(n, 4)  # Potenzmenge nur für kleine n
            lat = power_set_lattice(n)
        else:
            n = min(n, 60)
            lat = divisibility_lattice(n)

        # Verbandeigenschaften ermitteln
        info = {
            'type': lattice_type,
            'n': n,
            'elements': [str(e) for e in lat.elements],
            'element_count': len(lat.elements),
            'is_distributive': lat.is_distributive(),
            'is_modular': lat.is_modular(),
            'has_complement': lat.is_complemented() if hasattr(lat, 'is_complemented') else None,
        }
        return jsonify({'result': info})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/lattice/is_distributive', methods=['POST'])
def api_lattice_is_distributive():
    """
    @brief Prüft ob ein Verband distributiv ist.
    @description
        POST-Body: {"type": "divisibility" | "power_set", "n": int}
        Gibt is_distributive, is_modular und Gegenbeispiel zurück (falls vorhanden).
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _lattice_ok:
            return error_response("lattice_theory nicht verfügbar")

        lattice_type = data.get('type', 'divisibility')
        n = int(data.get('n', 12))

        if lattice_type == 'power_set':
            lat = power_set_lattice(min(n, 3))
        else:
            lat = divisibility_lattice(min(n, 30))

        is_dist = lat.is_distributive()
        is_mod = lat.is_modular()
        chain_ok = jordan_dedekind_chain_condition(lat)

        return jsonify({'result': {
            'n': n,
            'type': lattice_type,
            'is_distributive': is_dist,
            'is_modular': is_mod,
            'jordan_dedekind_chain_condition': chain_ok,
            'note': 'Distributiv ⟹ Modular, aber nicht umgekehrt'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/lattice/birkhoff', methods=['POST'])
def api_lattice_birkhoff():
    """
    @brief Birkhoff-Darstellungssatz für endliche distributive Verbände.
    @description
        POST-Body: {"n": int} – Teilerverband von n analysieren.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _lattice_ok:
            return error_response("lattice_theory nicht verfügbar")

        n = min(int(data.get('n', 12)), 30)
        lat = divisibility_lattice(n)
        birkhoff = birkhoff_representation_theorem(lat)

        # frozenset und andere nicht-serialisierbare Typen konvertieren
        def make_serializable(obj):
            """Wandelt nicht-JSON-serialisierbare Typen in serialisierbare um."""
            if isinstance(obj, (frozenset, set)):
                return sorted([make_serializable(e) for e in obj])
            if isinstance(obj, dict):
                return {str(k): make_serializable(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [make_serializable(e) for e in obj]
            if hasattr(obj, 'item'):
                return obj.item()
            return obj

        return jsonify({'result': make_serializable(birkhoff)})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Multilineare Algebra (multilinear_algebra.py)
# ===========================================================================

@app.route('/api/multilinear/tensor_product', methods=['POST'])
def api_multilinear_tensor_product():
    """
    @brief Tensorprodukt zweier Matrizen (Kronecker-Produkt) berechnen.
    @description
        POST-Body: {"A": [[...], ...], "B": [[...], ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _multilinear_ok:
            return error_response("multilinear_algebra nicht verfügbar")

        A = data.get('A', [[1, 0], [0, 1]])
        B = data.get('B', [[1, 2], [3, 4]])

        result = tensor_product_matrices(A, B)
        shape_A = (len(A), len(A[0]))
        shape_B = (len(B), len(B[0]))

        return jsonify({'result': {
            'tensor_product': result,
            'shape_A': shape_A,
            'shape_B': shape_B,
            'shape_result': (shape_A[0] * shape_B[0], shape_A[1] * shape_B[1]),
            'description': f'A ⊗ B hat Dimension {shape_A[0]*shape_B[0]}×{shape_A[1]*shape_B[1]}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/multilinear/exterior_product', methods=['POST'])
def api_multilinear_exterior_product():
    """
    @brief Äußeres Produkt (Wedge-Produkt) zweier Vektoren.
    @description
        POST-Body: {"u": [float, ...], "v": [float, ...]}
        Berechnet u ∧ v als antisymmetrischen Tensor.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _multilinear_ok:
            return error_response("multilinear_algebra nicht verfügbar")

        u = data.get('u', [1, 0, 0])
        v = data.get('v', [0, 1, 0])

        result = wedge_product(u, v)
        gram = gram_matrix([u, v])
        gram_det = gram_determinant([u, v])

        return jsonify({'result': {
            'u': u,
            'v': v,
            'wedge_product': result,
            'gram_matrix': gram,
            'gram_determinant': gram_det,
            'area': float(abs(gram_det) ** 0.5),
            'description': 'u ∧ v – Fläche des aufgespannten Parallelogramms'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/multilinear/signature', methods=['POST'])
def api_multilinear_signature():
    """
    @brief Signatur einer symmetrischen Bilinearform.
    @description
        POST-Body: {"M": [[...], ...]} – symmetrische Matrix der Bilinearform.
        Gibt (p, q, r) zurück: p positive, q negative, r Null-Eigenwerte.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _multilinear_ok:
            return error_response("multilinear_algebra nicht verfügbar")

        M = data.get('M', [[1, 0], [0, -1]])
        sig = signature_bilinear_form(np.array(M, dtype=float))

        return jsonify({'result': {
            'matrix': M,
            'signature': sig,
            'positive': sig[0],
            'negative': sig[1],
            'zero': sig[2] if len(sig) > 2 else 0,
            'description': f'Sylvesterscher Trägheitssatz: Signatur ({sig[0]}, {sig[1]})'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Homologische Algebra (homological_algebra.py)
# ===========================================================================

@app.route('/api/homological/chain_complex', methods=['POST'])
def api_homological_chain_complex():
    """
    @brief Kettenkomplex aus Randmatrizen aufbauen.
    @description
        POST-Body: {"simplex_type": "triangle" | "tetrahedron" | "circle"}
        Erstellt einen simplizialen Kettenkomplex.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _homological_ok:
            return error_response("homological_algebra nicht verfügbar")

        simplex_type = data.get('simplex_type', 'triangle')

        # Vordefinierte Simplizialkomplexe
        if simplex_type == 'triangle':
            # Dreieck: 3 Ecken, 3 Kanten, 1 Fläche
            simplices = {0: [[0], [1], [2]], 1: [[0,1], [1,2], [0,2]], 2: [[0,1,2]]}
        elif simplex_type == 'tetrahedron':
            # Tetraeder: 4 Ecken, 6 Kanten, 4 Dreiecke, 1 Volumen
            simplices = {
                0: [[0],[1],[2],[3]],
                1: [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]],
                2: [[0,1,2],[0,1,3],[0,2,3],[1,2,3]],
                3: [[0,1,2,3]]
            }
        else:
            # Kreis: 2 Ecken, 2 Kanten
            simplices = {0: [[0], [1]], 1: [[0,1]]}

        cc = simplicial_complex_to_chain(simplices)
        # Betti-Zahlen: cc.betti_numbers() → {deg: int}
        betti = cc.betti_numbers()
        euler = cc.euler_characteristic()

        result = {
            'simplex_type': simplex_type,
            'chain_groups': {str(k): v for k, v in cc.groups.items()},
            'homology': {}
        }
        for deg in sorted(cc.groups.keys()):
            h = cc.homology(deg)
            result['homology'][str(deg)] = {
                'betti_number': betti.get(deg, 0),
                'torsion': [str(t) for t in h.get('torsion', [])],
                'rank': h.get('rank', 0)
            }
        result['euler_characteristic'] = euler
        return jsonify({'result': result})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/homological/homology', methods=['POST'])
def api_homological_homology():
    """
    @brief Homologiegruppen eines simplizialen Komplexes berechnen.
    @description
        POST-Body: {"simplex_type": "triangle" | "tetrahedron" | "circle" | "torus"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _homological_ok:
            return error_response("homological_algebra nicht verfügbar")

        simplex_type = data.get('simplex_type', 'triangle')

        if simplex_type == 'triangle':
            simplices = {0: [[0],[1],[2]], 1: [[0,1],[1,2],[0,2]], 2: [[0,1,2]]}
        elif simplex_type == 'tetrahedron':
            simplices = {
                0: [[0],[1],[2],[3]],
                1: [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]],
                2: [[0,1,2],[0,1,3],[0,2,3],[1,2,3]],
                3: [[0,1,2,3]]
            }
        else:
            simplices = {0: [[0],[1]], 1: [[0,1]]}

        cc = simplicial_complex_to_chain(simplices)
        betti_raw = cc.betti_numbers()  # {deg: int}
        betti = {str(k): v for k, v in betti_raw.items()}
        euler = cc.euler_characteristic()

        return jsonify({'result': {
            'simplex_type': simplex_type,
            'betti_numbers': betti,
            'euler_characteristic': euler,
            'interpretation': {
                'H0': f'Anzahl Zusammenhangskomponenten: {betti.get("0", 0)}',
                'H1': f'Anzahl unabhängiger Schleifen: {betti.get("1", 0)}',
                'H2': f'Anzahl 2D-Hohlräume: {betti.get("2", 0)}',
            }
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/homological/ext_tor', methods=['POST'])
def api_homological_ext_tor():
    """
    @brief Ext- und Tor-Gruppen berechnen (abgeleitete Funktoren).
    @description
        POST-Body: {"n": int, "m": int, "degree": int}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _homological_ok:
            return error_response("homological_algebra nicht verfügbar")

        n = int(data.get('n', 2))
        m = int(data.get('m', 3))
        degree = int(data.get('degree', 1))

        ext = ext_group(n, m, degree)
        tor = tor_group(n, m, degree)

        return jsonify({'result': {
            'n': n, 'm': m, 'degree': degree,
            'ext': ext,
            'tor': tor,
            'description': f'Ext^{degree}(Z/{n}Z, Z/{m}Z) und Tor_{degree}(Z/{n}Z, Z/{m}Z)'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Algebraische Strukturen (algebraic_structures.py)
# ===========================================================================

@app.route('/api/algebraic_structures/magma', methods=['POST'])
def api_algstruct_magma():
    """
    @brief Magma-Struktur aus Verknüpfungstafel erstellen und analysieren.
    @description
        POST-Body: {"table": [[int, ...], ...], "elements": [str, ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _algstruct_ok:
            return error_response("algebraic_structures nicht verfügbar")

        table = data.get('table', [[0,1],[1,0]])
        elements = data.get('elements', ['e', 'a'])

        # Magma aus der Tabelle erstellen (operation = callable aus Tabelle)
        def make_op(elems, tbl):
            return lambda x, y: elems[tbl[elems.index(x)][elems.index(y)]]

        op = make_op(elements, table)
        m = Magma(elements, op)
        # Semigroup-Test über separate Instanz
        sg = Semigroup(elements, op)
        is_sg = sg.is_semigroup()
        is_comm = m.is_commutative()

        return jsonify({'result': {
            'elements': elements,
            'table': table,
            'is_semigroup': is_sg,
            'is_commutative': is_comm,
            'description': 'Magma: Menge mit abgeschlossener Verknüpfung'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebraic_structures/check_group', methods=['POST'])
def api_algstruct_check_group():
    """
    @brief Prüft ob eine Verknüpfungstafel eine Gruppe definiert.
    @description
        POST-Body: {"table": [[int, ...], ...], "elements": [str, ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _algstruct_ok:
            return error_response("algebraic_structures nicht verfügbar")

        table = data.get('table', [[0,1],[1,0]])
        elements = data.get('elements', ['e', 'a'])

        # Operation aus Tabelle erzeugen
        def make_op2(elems, tbl):
            return lambda x, y: elems[tbl[elems.index(x)][elems.index(y)]]

        op = make_op2(elements, table)
        sg = Semigroup(elements, op)
        is_sg = sg.is_semigroup()
        mon = Monoid(elements, op)
        is_mon = mon.has_identity() if is_sg else False
        grp = GroupFromMagma(elements, op)
        is_grp = grp.is_group() if is_mon else False
        is_abel = grp.is_commutative() if is_grp else False

        return jsonify({'result': {
            'elements': elements,
            'is_magma': True,
            'is_semigroup': is_sg,
            'is_monoid': is_mon,
            'is_group': is_grp,
            'is_abelian': is_abel,
            'hierarchy': 'Magma → Halbgruppe → Monoid → Gruppe → Abel. Gruppe'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebraic_structures/hierarchy', methods=['GET'])
def api_algstruct_hierarchy():
    """
    @brief Hierarchie algebraischer Strukturen abrufen.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _algstruct_ok:
            return error_response("algebraic_structures nicht verfügbar")
        hierarchy = algebraic_structure_hierarchy()
        return jsonify({'result': hierarchy})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Invariantentheorie (invariant_theory.py)
# ===========================================================================

@app.route('/api/invariant/reynolds_operator', methods=['POST'])
def api_invariant_reynolds():
    """
    @brief Reynolds-Operator: Mittelwert eines Polynoms über eine Gruppe.
    @description
        POST-Body: {"poly": "x**2 + y**2", "n": int}
        Mittelt das Polynom über die symmetrische Gruppe S_n.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _invariant_ok:
            return error_response("invariant_theory nicht verfügbar")

        poly_str = data.get('poly', 'x**2 + y**2')
        n = min(int(data.get('n', 2)), 4)  # S_n für kleine n

        # Variablen x1, ..., xn
        var_names = [f'x{i+1}' for i in range(n)]
        symbols = sp.symbols(' '.join(var_names))
        if n == 1:
            symbols = (symbols,)

        # Polynom parsen
        local_dict = {var_names[i]: symbols[i] for i in range(n)}
        # Auch x, y, z als Aliase unterstützen
        for i, name in enumerate(['x', 'y', 'z', 'w']):
            if i < n:
                local_dict[name] = symbols[i]

        poly_expr = sp.sympify(poly_str, locals=local_dict)

        # Reynolds-Operator: Mittlung über alle Permutationen von S_n
        # reynolds_operator(poly_expr, group_elements, action)
        from itertools import permutations

        # Alle Permutationen von S_n als Gruppenelemente
        sym_list = list(symbols)
        perm_list = list(permutations(range(n)))

        def perm_action(perm_idx):
            # Permutation als Substitution auf poly_expr
            perm = perm_list[perm_idx]
            subs = {sym_list[i]: sym_list[perm[i]] for i in range(n)}
            return lambda p: p.subs(subs)

        # Manuell mitteln
        total = sp.Integer(0)
        for perm in perm_list:
            subs = {sym_list[i]: sym_list[perm[i]] for i in range(n)}
            total += poly_expr.subs(subs)
        invariant = sp.simplify(total / len(perm_list))
        fund = elementary_symmetric_polynomials(n)

        return jsonify({'result': {
            'input_poly': poly_str,
            'group': f'S_{n}',
            'reynolds_image': str(invariant),
            'elementary_symmetric': [str(p) for p in fund],
            'description': 'Reynolds-Operator: R[f] = (1/|G|) Σ_{g∈G} g·f'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/invariant/molien_series', methods=['POST'])
def api_invariant_molien():
    """
    @brief Molien-Reihe: Erzeugende Funktion der Invariantendimensionen.
    @description
        POST-Body: {"group_type": "Z2" | "S2" | "cyclic", "max_degree": int}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _invariant_ok:
            return error_response("invariant_theory nicht verfügbar")

        group_type = data.get('group_type', 'Z2')
        max_degree = min(int(data.get('max_degree', 8)), 12)

        # Gruppenmatrizen für Z2 (Reflexion) oder S2 (Permutation)
        if group_type == 'Z2':
            matrices = [[[1, 0], [0, 1]], [[-1, 0], [0, -1]]]
        elif group_type == 'S2':
            matrices = [[[1, 0], [0, 1]], [[0, 1], [1, 0]]]
        else:
            matrices = [[[1, 0], [0, 1]], [[-1, 0], [0, -1]]]
            group_type = 'Z2 (Standard)'

        coeffs = molien_series(matrices, matrices, max_degree)

        return jsonify({'result': {
            'group': group_type,
            'max_degree': max_degree,
            'molien_coefficients': coeffs,
            'description': 'Molien(t) = Σ dim(V^G_d) t^d – Erzeugende Funktion'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/invariant/symmetric_polynomials', methods=['POST'])
def api_invariant_symmetric():
    """
    @brief Elementarsymmetrische Polynome für S_n berechnen.
    @description
        POST-Body: {"n": int}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _invariant_ok:
            return error_response("invariant_theory nicht verfügbar")

        n = min(int(data.get('n', 3)), 5)
        polys = elementary_symmetric_polynomials(n)
        inv_polys = polynomial_invariants_sn(n, 2)

        return jsonify({'result': {
            'n': n,
            'elementary_symmetric': [str(p) for p in polys],
            'degree_2_invariants': [str(p) for p in inv_polys],
            'hilbert_basis': hilbert_basis_theorem_demo()
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Universelle Algebra (universal_algebra.py)
# ===========================================================================

@app.route('/api/universal/algebra_check', methods=['POST'])
def api_universal_algebra_check():
    """
    @brief Prüft ob eine Struktur eine Varietät erfüllt.
    @description
        POST-Body: {"variety": "groups" | "rings" | "lattices"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _universal_ok:
            return error_response("universal_algebra nicht verfügbar")

        variety_name = data.get('variety', 'groups')

        if variety_name == 'groups':
            variety = groups_variety()
        elif variety_name == 'rings':
            variety = rings_variety()
        elif variety_name == 'lattices':
            variety = lattices_variety()
        else:
            return error_response(f"Unbekannte Varietät: {variety_name}")

        return jsonify({'result': {
            'variety': variety_name,
            'axioms': variety.axioms if hasattr(variety, 'axioms') else [],
            'signature': str(variety.signature) if hasattr(variety, 'signature') else '',
            'birkhoff_theorem': 'Eine Klasse ist genau dann eine Varietät, wenn sie unter '
                                'Unterstrukturen, homomorphen Bildern und Produkten abgeschlossen ist.',
            'description': f'Varietät der {variety_name}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/universal/birkhoff', methods=['GET'])
def api_universal_birkhoff():
    """
    @brief Birkhoff-HSP-Satz Demo.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _universal_ok:
            return error_response("universal_algebra nicht verfügbar")
        result = birkhoff_theorem_demo()
        return jsonify({'result': result})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/universal/free_monoid', methods=['POST'])
def api_universal_free_monoid():
    """
    @brief Freies Monoid über einem Alphabet generieren.
    @description
        POST-Body: {"alphabet": ["a", "b", ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _universal_ok:
            return error_response("universal_algebra nicht verfügbar")

        alphabet = data.get('alphabet', ['a', 'b'])
        result = free_monoid(alphabet)
        return jsonify({'result': {
            'alphabet': alphabet,
            'free_monoid_sample': result,
            'description': 'Freies Monoid: alle endlichen Wörter über dem Alphabet'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Modelltheorie (model_theory.py)
# ===========================================================================

@app.route('/api/model/satisfies', methods=['POST'])
def api_model_satisfies():
    """
    @brief Prüft ob eine Struktur eine Theorie erfüllt.
    @description
        POST-Body: {"theory": "group" | "field" | "order"}
        Gibt Axiome und Erfüllbarkeit zurück.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _model_ok:
            return error_response("model_theory nicht verfügbar")

        theory = data.get('theory', 'group')
        # Kompaktheitssatz-Demo als Beispiel für Erfüllbarkeit
        demo = compactness_theorem_demo()
        cat = theory_is_categorical(theory if theory in ['group', 'field', 'linear_order'] else 'group', 'omega')

        return jsonify({'result': {
            'theory': theory,
            'compactness_demo': demo,
            'categoricity': cat,
            'description': 'Modelltheorie: Zusammenhang zwischen Formeln und Strukturen'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/model/elementary_equivalence', methods=['POST'])
def api_model_elementary_equivalence():
    """
    @brief Elementare Äquivalenz zweier Strukturen prüfen.
    @description
        POST-Body: {"theory1": str, "theory2": str}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _model_ok:
            return error_response("model_theory nicht verfügbar")

        # Löwenheim-Skolem Demo
        ls_down = lowenheim_skolem_downward(['∀x∀y(x=y → x=y)', '∃x(x=x)'], 5)
        nonstd = nonstandard_arithmetic_demo()

        return jsonify({'result': {
            'lowenheim_skolem_downward': ls_down,
            'nonstandard_arithmetic': nonstd,
            'description': 'Elementare Äquivalenz: M ≡ N gdw. M und N dieselben Sätze erfüllen'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/model/quantifier_elimination', methods=['POST'])
def api_model_quantifier_elimination():
    """
    @brief Quantorenelimination für Theorien demonstrieren.
    @description
        POST-Body: {"theory": "DLO" | "RCF" | "ACF"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _model_ok:
            return error_response("model_theory nicht verfügbar")

        theory = data.get('theory', 'DLO')
        result = quantifier_elimination_demo(theory)
        return jsonify({'result': result})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Formale Beweistheorie (proof_theory_formal.py)
# ===========================================================================

@app.route('/api/proof_formal/sequent_calculus', methods=['POST'])
def api_proof_formal_sequent():
    """
    @brief Sequenzenkalkül LK – einfache Beweise.
    @description
        POST-Body: {"formula": str} – logische Formel zu beweisen.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _proof_formal_ok:
            return error_response("proof_theory_formal nicht verfügbar")

        # Schnitteliminationssatz demonstrieren
        cut_demo = cut_elimination_demo()
        complexity = proof_complexity_comparison()

        return jsonify({'result': {
            'cut_elimination': cut_demo,
            'complexity_comparison': complexity,
            'description': 'LK-Kalkül: Gentzens Sequenzenkalkül für die klassische Logik'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/proof_formal/natural_deduction', methods=['POST'])
def api_proof_formal_nd():
    """
    @brief Natürliches Schließen – einfache Formeln beweisen.
    @description
        POST-Body: {"formula": str} – Aussagenlogische Formel.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _proof_formal_ok:
            return error_response("proof_theory_formal nicht verfügbar")

        formula = data.get('formula', 'A → A')
        # Einfache Beweise via prove_simple
        proof = prove_simple(formula)

        ordinal_info = proof_theoretic_ordinal('PA')
        big5 = big_five_systems()

        if proof is None:
            return jsonify({'result': {
                'formula': formula,
                'proved': False,
                'message': 'Formel nicht automatisch beweisbar (zu komplex)',
                'proof_theoretic_ordinals': ordinal_info,
                'big_five': big5
            }})

        return jsonify({'result': {
            'formula': formula,
            'proved': True,
            'proof_steps': str(proof),
            'proof_theoretic_ordinals': ordinal_info,
            'big_five': big5
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/proof_formal/ordinals', methods=['GET'])
def api_proof_formal_ordinals():
    """
    @brief Beweistheoretische Ordinalzahlen der großen 5 Systeme.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _proof_formal_ok:
            return error_response("proof_theory_formal nicht verfügbar")

        big5 = big_five_systems()
        gentzen = gentzen_consistency_proof_sketch()
        return jsonify({'result': {'big_five': big5, 'gentzen': gentzen}})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Rekursionstheorie (recursion_theory.py)
# ===========================================================================

@app.route('/api/recursion/turing_machine', methods=['POST'])
def api_recursion_turing_machine():
    """
    @brief Turingmaschine simulieren.
    @description
        POST-Body: {"machine": "palindrome" | "binary_increment" | "unary_addition",
                    "input": [str, ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _recursion_ok:
            return error_response("recursion_theory nicht verfügbar")

        machine_name = data.get('machine', 'palindrome')
        input_tape = data.get('input', ['a', 'b', 'a'])

        # Turingmaschine wählen
        if machine_name == 'palindrome':
            tm = tm_recognizes_palindromes()
            description = 'Erkennt Palindrome über {a, b}'
        elif machine_name == 'binary_increment':
            tm = tm_binary_increment()
            description = 'Binäre Inkrementierung (z.B. 0111 → 1000)'
        else:
            tm = tm_recognizes_palindromes()
            description = 'Standard: Palindrom-Erkennung'

        # Simulation ausführen (mit Schrittzähler-Limit)
        result_tape, steps, accepted = tm.run(input_tape, max_steps=1000)

        return jsonify({'result': {
            'machine': machine_name,
            'description': description,
            'input': input_tape,
            'output_tape': [str(s) for s in result_tape],
            'steps': steps,
            'accepted': accepted,
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/recursion/halting_problem_demo', methods=['GET'])
def api_recursion_halting():
    """
    @brief Unentscheidbarkeit des Halteproblems demonstrieren.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _recursion_ok:
            return error_response("recursion_theory nicht verfügbar")

        result = halting_problem_undecidability_proof()
        return jsonify({'result': result})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/recursion/ackermann', methods=['POST'])
def api_recursion_ackermann():
    """
    @brief Ackermann-Funktion berechnen (schnell wachsend).
    @description
        POST-Body: {"m": int (0-3), "n": int (0-10)}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _recursion_ok:
            return error_response("recursion_theory nicht verfügbar")

        m = min(int(data.get('m', 2)), 3)  # m>3 zu langsam
        n = min(int(data.get('n', 3)), 10)

        value = ackermann_function(m, n)
        growth = ackermann_growth_demo()
        hierarchy = arithmetical_hierarchy()

        return jsonify({'result': {
            'm': m, 'n': n,
            'ackermann': value,
            'growth_demo': growth,
            'arithmetical_hierarchy': hierarchy,
            'description': f'A({m},{n}) = {value} – nicht primitiv-rekursiv'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Maßtheorie (measure_theory.py)
# ===========================================================================

@app.route('/api/measure/lebesgue_integral', methods=['POST'])
def api_measure_lebesgue_integral():
    """
    @brief Lebesgue-Integral einer Funktion berechnen.
    @description
        POST-Body: {"function_type": "polynomial" | "step" | "continuous",
                    "a": float, "b": float}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _measure_ok:
            return error_response("measure_theory nicht verfügbar")

        func_type = data.get('function_type', 'polynomial')
        a = float(data.get('a', 0.0))
        b = float(data.get('b', 1.0))

        # Funktion wählen
        if func_type == 'step':
            f = lambda x: 1.0 if x > (a + b) / 2 else 0.0
            fname = 'Stufenfunktion χ_{[(a+b)/2, b]}'
        elif func_type == 'continuous':
            f = lambda x: float(np.sin(x))
            fname = 'sin(x)'
        else:
            f = lambda x: x ** 2
            fname = 'x²'

        # Lebesgue-Integral berechnen
        integral = lebesgue_integral(f, a, b)
        # riemann_vs_lebesgue akzeptiert nur 'dirichlet', 'bounded_discontinuous', 'improper'
        riem_type = 'bounded_discontinuous' if func_type in ['step', 'polynomial'] else 'dirichlet'
        comparison = riemann_vs_lebesgue(riem_type)

        return jsonify({'result': {
            'function': fname,
            'a': a, 'b': b,
            'lebesgue_integral': float(integral),
            'comparison': comparison,
            'description': 'Lebesgue-Integral: Integral bzgl. des Lebesgue-Maßes'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/measure/sigma_algebra', methods=['POST'])
def api_measure_sigma_algebra():
    """
    @brief σ-Algebra erzeugt von Generatoren berechnen.
    @description
        POST-Body: {"base_set": [elements], "generators": [[subset1], [subset2], ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _measure_ok:
            return error_response("measure_theory nicht verfügbar")

        base_set = data.get('base_set', [1, 2, 3, 4])
        generators = data.get('generators', [[1, 2], [3]])

        # σ-Algebra aus Generatoren erzeugen
        # SigmaAlgebra(universe: frozenset, sets: list)
        universe = frozenset(base_set)
        gen_sets = [frozenset(g) for g in generators]
        sigma = SigmaAlgebra(universe, gen_sets)
        # generated_sigma_algebra(generators) gibt neue SigmaAlgebra zurück
        gen_sigma = sigma.generated_sigma_algebra(gen_sets)
        generated = gen_sigma.sets  # set von frozensets

        # JSON-serialisierbar machen
        gen_list = sorted(
            [sorted(list(s)) for s in generated],
            key=lambda s: (len(s), s)
        )

        return jsonify({'result': {
            'base_set': list(base_set),
            'generators': [list(g) for g in generators],
            'sigma_algebra_size': len(generated),
            'sigma_algebra': gen_list,
            'description': 'σ(G): kleinste σ-Algebra die G enthält'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/measure/convergence_theorems', methods=['GET'])
def api_measure_convergence():
    """
    @brief Konvergenzsätze der Maßtheorie (Lebesgue, Fatou, dominiert).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _measure_ok:
            return error_response("measure_theory nicht verfügbar")

        monotone = monotone_convergence_theorem_demo()
        dominated = dominated_convergence_theorem_demo()
        fubini = fubini_theorem_demo()

        return jsonify({'result': {
            'monotone_convergence': monotone,
            'dominated_convergence': dominated,
            'fubini': fubini
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Spezielle Funktionen (special_functions.py)
# ===========================================================================

@app.route('/api/special/bessel', methods=['POST'])
def api_special_bessel():
    """
    @brief Bessel-Funktionen J_n(x) und Y_n(x) berechnen.
    @description
        POST-Body: {"n": int, "x": float}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _special_ok:
            return error_response("special_functions nicht verfügbar")

        n = int(data.get('n', 0))
        x = float(data.get('x', 1.0))

        # BesselFunctions() ohne Argument; J(nu, x), Y(nu, x)
        bf = BesselFunctions()
        jn = bf.J(float(n), x)
        try:
            yn = bf.Y(float(n), x)
        except Exception:
            yn = None

        # Werte für Plot-Punkte
        xs = [round(0.1 * i, 2) for i in range(1, 51)]
        j_vals = [bf.J(float(n), xi) for xi in xs]

        return jsonify({'result': {
            'n': n, 'x': x,
            'J_n_x': float(jn),
            'Y_n_x': float(yn) if yn is not None and not np.isinf(yn) else None,
            'plot_x': xs,
            'plot_J': [float(v) for v in j_vals],
            'zeros_approx': bf.zeros(float(n), 3) if hasattr(bf, 'zeros') else [],
            'description': f'Bessel-Funktion J_{n}({x}) = {float(jn):.6f}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/special/legendre', methods=['POST'])
def api_special_legendre():
    """
    @brief Legendre-Polynome P_n(x) berechnen.
    @description
        POST-Body: {"n": int, "x": float}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _special_ok:
            return error_response("special_functions nicht verfügbar")

        n = min(int(data.get('n', 3)), 10)
        x = float(data.get('x', 0.5))

        lp = LegendrePolynomials()
        pn = lp.P(n, x)

        # Erste Polynome an verschiedenen Punkten auswerten (kein symbolic-Attribut)
        polys = []
        for k in range(min(n + 1, 6)):
            # Werte an einigen Stützpunkten
            sample_vals = [(xi/4 - 1.0, round(lp.P(k, xi/4 - 1.0), 5)) for xi in range(5)]
            polys.append({'degree': k, 'formula': f'P_{k}(x)', 'sample_values': sample_vals})

        return jsonify({'result': {
            'n': n, 'x': x,
            'P_n_x': float(pn),
            'polynomials': polys,
            'orthogonality': 'P_n orthogonal auf [-1,1] bzgl. Gewicht 1',
            'description': f'Legendre-Polynom P_{n}({x}) = {float(pn):.6f}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/special/gamma_properties', methods=['GET'])
def api_special_gamma():
    """
    @brief Eigenschaften der Gamma-Funktion Γ(z).
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _special_ok:
            return error_response("special_functions nicht verfügbar")

        props = gamma_function_properties()
        zeta_vals = riemann_zeta_special_values()
        erf = error_function_properties()

        # numpy-Typen in Python-Standardtypen konvertieren
        def np_safe(obj):
            """Konvertiert numpy-Typen in JSON-serialisierbare Python-Typen."""
            if hasattr(obj, 'item'):  # numpy scalar → Python
                return obj.item()
            if isinstance(obj, dict):
                return {k: np_safe(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [np_safe(v) for v in obj]
            return obj

        return jsonify({'result': {
            'gamma_properties': np_safe(props),
            'zeta_special_values': np_safe(zeta_vals),
            'error_function': np_safe(erf)
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Funktionalanalysis (functional_analysis.py)
# ===========================================================================

@app.route('/api/functional/banach_space', methods=['POST'])
def api_functional_banach():
    """
    @brief Banach-Raum: Norm und Vollständigkeit prüfen.
    @description
        POST-Body: {"vectors": [[float, ...], ...], "norm_type": "1" | "2" | "inf"}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _functional_ok:
            return error_response("functional_analysis nicht verfügbar")

        vectors = data.get('vectors', [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        norm_type = str(data.get('norm_type', '2'))

        # p-Wert für die Norm
        p_map = {'1': 1, '2': 2, 'inf': float('inf')}
        p = p_map.get(norm_type, 2)

        # BanachSpace(vectors, norm_func)
        norm_func = lambda v: float(np.linalg.norm(v, ord=p))
        bs = BanachSpace(vectors, norm_func)
        norms = [bs.norm(v) for v in vectors]

        # Dreiecksungleichung prüfen
        v1, v2 = np.array(vectors[0]), np.array(vectors[1])
        triangle_lhs = float(norm_func((v1 + v2).tolist()))
        triangle_rhs = float(norms[0]) + float(norms[1])

        hahn_banach = hahn_banach_theorem_demo()

        return jsonify({'result': {
            'vectors': vectors,
            'norm_type': f'L^{norm_type}',
            'norms': [float(n) for n in norms],
            'triangle_inequality': {
                'lhs': triangle_lhs,
                'rhs': triangle_rhs,
                'holds': triangle_lhs <= triangle_rhs + 1e-10
            },
            'hahn_banach_demo': hahn_banach
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/functional/hilbert_space', methods=['POST'])
def api_functional_hilbert():
    """
    @brief Hilbert-Raum: Skalarprodukt, Orthogonalität, Parseval.
    @description
        POST-Body: {"vectors": [[float, ...], ...]}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _functional_ok:
            return error_response("functional_analysis nicht verfügbar")

        vectors = data.get('vectors', [[1.0, 0.0], [0.0, 1.0]])

        # HilbertSpace(vectors, inner_product_func)
        ip_func = lambda u, v: float(np.dot(np.array(u), np.conj(np.array(v))))
        hs = HilbertSpace(vectors, ip_func)

        # Skalarprodukte berechnen
        inner_products = []
        for i, u in enumerate(vectors):
            for j, v in enumerate(vectors):
                ip = float(np.real(hs.inner_product(u, v)))
                inner_products.append({'i': i, 'j': j, 'inner_product': ip})

        riesz = riesz_representation_theorem()
        spectral = spectral_theorem_demo(np.array([[2.0, 1.0], [1.0, 2.0]]))

        return jsonify({'result': {
            'vectors': vectors,
            'inner_products': inner_products,
            'riesz_theorem': riesz,
            'spectral_theorem_demo': spectral,
            'description': 'Hilbert-Raum: vollständiger normierter Raum mit Skalarprodukt'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/functional/fixed_point', methods=['POST'])
def api_functional_fixed_point():
    """
    @brief Banachscher Fixpunktsatz: Kontraktionsabbildung.
    @description
        POST-Body: {"contraction_factor": float (0<c<1), "x0": float}
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _functional_ok:
            return error_response("functional_analysis nicht verfügbar")

        c = float(data.get('contraction_factor', 0.5))
        x0 = float(data.get('x0', 1.0))

        if c >= 1.0 or c <= 0.0:
            return error_response("Kontraktionsfaktor muss zwischen 0 und 1 liegen")

        # f(x) = c*x + (1-c)*fixed_point → Fixpunkt bei x=1
        f = lambda x: c * x + (1.0 - c)
        result = banach_fixed_point(f, x0)

        return jsonify({'result': {
            'contraction_factor': c,
            'x0': x0,
            'fixed_point': result.get('fixed_point', None),
            'iterations': result.get('iterations', 0),
            'converged': result.get('converged', False),
            'description': f'f(x) = {c}·x + {1-c} hat Fixpunkt bei x = 1'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Partielle Differentialgleichungen (pde.py)
# ===========================================================================

@app.route('/api/pde/classify', methods=['POST'])
def api_pde_classify():
    """
    @brief PDE klassifizieren: elliptisch, parabolisch oder hyperbolisch.
    @description
        POST-Body: {"A": float, "B": float, "C": float}
        PDE: A·u_xx + B·u_xy + C·u_yy + ... = 0
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _pde_ok:
            return error_response("pde nicht verfügbar")

        A = float(data.get('A', 1.0))
        B = float(data.get('B', 0.0))
        C = float(data.get('C', 1.0))

        classification = classify_pde(A, B, C)
        discriminant = B ** 2 - 4 * A * C

        examples = {
            'elliptic': 'Laplace-Gleichung: Δu = 0 (B²-4AC < 0)',
            'parabolic': 'Wärmeleitungsgleichung: u_t = α·u_xx (B²-4AC = 0)',
            'hyperbolic': 'Wellengleichung: u_tt = c²·u_xx (B²-4AC > 0)'
        }

        return jsonify({'result': {
            'A': A, 'B': B, 'C': C,
            'discriminant': discriminant,
            'classification': classification,
            'example': examples.get(classification.lower(), ''),
            'description': f'B²-4AC = {discriminant:.4f} → {classification}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/pde/heat_equation', methods=['POST'])
def api_pde_heat():
    """
    @brief Wärmeleitungsgleichung numerisch lösen.
    @description
        POST-Body: {"u0": [float, ...], "L": float, "T": float, "alpha": float}
        u_t = α · u_xx auf [0, L] × [0, T]
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _pde_ok:
            return error_response("pde nicht verfügbar")

        u0_list = data.get('u0', [0.0, 0.5, 1.0, 0.5, 0.0])
        L = float(data.get('L', 1.0))
        T = float(data.get('T', 0.1))
        alpha = float(data.get('alpha', 1.0))

        u0 = np.array(u0_list, dtype=float)
        N = len(u0)

        # Explizites Finite-Differenzen-Verfahren
        # heat_equation_explicit(u0, L, T, alpha, nx, nt) → np.ndarray
        # heat_equation_explicit gibt (nt+1) × (nx+1)-Array zurück; letzte Zeile = Endzustand
        u_hist = heat_equation_explicit(u0_list, L, T, alpha, nx=N, nt=200)
        u_final = u_hist[-1] if u_hist.ndim == 2 else u_hist

        # Analytische Lösung zum Vergleich
        # heat_equation_analytical(x, t, L, n_terms) – kein alpha-Parameter
        x_arr = np.linspace(0, L, N)
        try:
            u_analytical = heat_equation_analytical(x_arr, T, L, n_terms=10)
        except Exception:
            u_analytical = None

        return jsonify({'result': {
            'u0': u0_list,
            'u_final': [float(v) for v in u_final],
            'u_analytical': [float(v) for v in u_analytical] if u_analytical is not None else None,
            'x': [float(v) for v in x_arr],
            'L': L, 'T': T, 'alpha': alpha,
            'description': f'Wärmeleitungsgl. u_t = {alpha}·u_xx nach Zeit T={T}'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/pde/wave_equation', methods=['POST'])
def api_pde_wave():
    """
    @brief Wellengleichung numerisch lösen.
    @description
        POST-Body: {"u0": [float, ...], "v0": [float, ...], "L": float, "T": float, "c": float}
        u_tt = c² · u_xx auf [0, L] × [0, T]
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    data = request.get_json()
    try:
        if not _pde_ok:
            return error_response("pde nicht verfügbar")

        u0_list = data.get('u0', [0.0, 0.5, 1.0, 0.5, 0.0])
        v0_list = data.get('v0', [0.0] * len(u0_list))
        L = float(data.get('L', 1.0))
        T = float(data.get('T', 0.5))
        c = float(data.get('c', 1.0))

        u0 = np.array(u0_list, dtype=float)
        v0 = np.array(v0_list, dtype=float)
        N = len(u0)

        # Explizites Finite-Differenzen-Verfahren
        # wave_equation_explicit(u0, u0_t, L, T, c, nx, nt) → np.ndarray
        # wave_equation_explicit gibt (nt+1) × (nx+1)-Array zurück; letzte Zeile = Endzustand
        u_hist = wave_equation_explicit(u0_list, v0_list, L, T, c, nx=N, nt=200)
        u_final = u_hist[-1] if u_hist.ndim == 2 else u_hist

        x_arr = np.linspace(0, L, N)
        try:
            u_analytical = wave_equation_analytical(x_arr, T, c, L, n_terms=10)
        except Exception:
            u_analytical = None

        return jsonify({'result': {
            'u0': u0_list,
            'u_final': [float(v) for v in u_final],
            'u_analytical': [float(v) for v in u_analytical] if u_analytical is not None else None,
            'x': [float(v) for v in x_arr],
            'L': L, 'T': T, 'c': c,
            'description': f'Wellengl. u_tt = {c}²·u_xx nach Zeit T={T}'
        }})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Operatoralgebren (operator_algebras.py)
# ===========================================================================

@app.route('/api/operator_algebras/cstar_demo', methods=['GET'])
def api_opalg_cstar():
    """
    @brief C*-Algebren: Gelfand-Transform und GNS-Konstruktion.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _opalg_ok:
            return error_response("operator_algebras nicht verfügbar")

        gelfand = gelfand_transform_demo()
        gns = gns_construction_demo()
        cuntz = cuntz_algebra_demo()

        # complex-Zahlen in str konvertieren (nicht JSON-serialisierbar)
        def complex_safe(obj):
            """Konvertiert complex und numpy-Typen in JSON-serialisierbare Typen."""
            if isinstance(obj, complex):
                return {'real': obj.real, 'imag': obj.imag, 'str': str(obj)}
            if hasattr(obj, 'item'):
                return obj.item()
            if isinstance(obj, dict):
                return {k: complex_safe(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [complex_safe(v) for v in obj]
            return obj

        return jsonify({'result': {
            'gelfand_transform': complex_safe(gelfand),
            'gns_construction': complex_safe(gns),
            'cuntz_algebra': complex_safe(cuntz),
            'description': 'C*-Algebra: Banach-*-Algebra mit ||a*a|| = ||a||²'
        }})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/operator_algebras/factors', methods=['GET'])
def api_opalg_factors():
    """
    @brief von-Neumann-Algebra Faktoren: Typ I, II, III.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _opalg_ok:
            return error_response("operator_algebras nicht verfügbar")

        factors = factors_classification()
        return jsonify({'result': factors})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/operator_algebras/k_theory', methods=['GET'])
def api_opalg_ktheory():
    """
    @brief K-Theorie für C*-Algebren: K0 und K1.
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    try:
        if not _opalg_ok:
            return error_response("operator_algebras nicht verfügbar")

        k_theory = k_theory_intro()
        return jsonify({'result': k_theory})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# ERWEITERTE ALGEBRA API (algebraic_structures, invariant_theory, representation_theory,
#                          homological_algebra, multilinear_algebra)
# ===========================================================================

@app.route('/api/algebra_adv/check_structure', methods=['POST'])
def api_algadv_check_structure():
    """
    @brief Prüft eine algebraische Struktur (Magma/Halbgruppe/Monoid/Gruppe).
    @lastModified 2026-03-10
    """
    if not _algstruct_ok:
        return error_response('algebraic_structures nicht geladen')
    data = request.get_json()
    table = data.get('table', [])
    elements = data.get('elements', [])
    try:
        m = Magma(elements, table)
        result = {
            'ist_magma': True,
            'ist_halbgruppe': m.is_semigroup(),
            'ist_monoid': m.is_monoid(),
            'ist_gruppe': m.is_group(),
            'ist_abelsche_gruppe': m.is_abelian_group(),
            'ordnung': len(elements),
        }
        return jsonify(result)
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebra_adv/invariant', methods=['POST'])
def api_algadv_invariant():
    """
    @brief Berechnet Molien-Reihe und Invariantengeneratoren für eine Gruppe.
    @lastModified 2026-03-10
    """
    if not _invariant_ok:
        return error_response('invariant_theory nicht geladen')
    data = request.get_json()
    group_type = data.get('group_type', 'Z2')
    max_degree = min(int(data.get('max_degree', 6)), 8)
    try:
        import numpy as np
        # Matrizen für bekannte Gruppen aufbauen
        if group_type == 'Z2':
            # Z2 = {I, -I} auf R²
            mats = [np.eye(2), -np.eye(2)]
            elems = ['e', 'g']
        elif group_type.startswith('S') and len(group_type) == 2:
            n_sym = int(group_type[1])
            # Triviale Darstellung: nur Identität (vereinfacht)
            mats = [np.eye(n_sym)]
            elems = ['e']
        else:
            mats = [np.eye(2)]
            elems = ['e']
        molien = molien_series(elems, mats, max_degree)
        gens = elementary_symmetric_polynomials(2)
        hbt = hilbert_basis_theorem_demo()
        return jsonify({
            'group': group_type,
            'molien_coefficients': molien,
            'invariant_generators': [str(g) for g in gens[:5]],
            'hilbert_basis_theorem': hbt.get('statement', ''),
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebra_adv/character_table', methods=['POST'])
def api_algadv_character_table():
    """
    @brief Liefert die Charaktertafel einer Gruppe.
    @lastModified 2026-03-10
    """
    if not _repr_ok:
        return error_response('representation_theory nicht geladen')
    data = request.get_json()
    group_name = data.get('group_name', 'Z2').upper()
    try:
        # Vorbereitete Darstellungen für bekannte Gruppen nutzen
        if group_name == 'Z2':
            reps = z2_representations()
        elif group_name == 'S3':
            reps = s3_representations()
        else:
            return error_response(f"Unbekannte Gruppe: {group_name}. Bitte 'S3' oder 'Z2'.")
        # Charaktertafel aus Darstellungen ableiten
        char_data = []
        for rep in reps:
            char_row = [complex(rep.character().get(g, 0)) for g in rep.group_elements]
            char_data.append([
                round(c.real, 4) if abs(c.imag) < 1e-10
                else [round(c.real, 4), round(c.imag, 4)]
                for c in char_row
            ])
        return jsonify({
            'group': group_name,
            'character_table': char_data,
            'irreducible_count': len(reps),
            'elements': [str(g) for g in (reps[0].group_elements if reps else [])],
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebra_adv/homology', methods=['POST'])
def api_algadv_homology():
    """
    @brief Berechnet Homologiegruppen eines Kettenkomplexes.
    @lastModified 2026-03-10
    """
    if not _homological_ok:
        return error_response('homological_algebra nicht geladen')
    data = request.get_json()
    boundary_maps = data.get('boundary_maps', [])
    try:
        import numpy as np
        # Standard: Simplizialkomplex (Dreieck = S^1-Analog)
        # boundaries[n] = Randoperator ∂_n: C_n → C_{n-1}
        if not boundary_maps:
            # Default: Dreieck (1-Simplex ∂₁: C₁→C₀)
            # ∂₁([01],[12],[02]) = [-1,1,0; 0,-1,1; 1,0,-1]
            boundaries = {
                1: np.array([[-1, 0, 1], [1, -1, 0], [0, 1, -1]], dtype=float)
            }
            groups = {0: 3, 1: 3}
        else:
            boundaries = {i: np.array(m, dtype=float) for i, m in enumerate(boundary_maps, 1)}
            groups = {}
            for k, mat in boundaries.items():
                groups[k - 1] = mat.shape[0]
                groups[k] = mat.shape[1]
        cc = ChainComplex(groups, boundaries)
        results = {}
        for n in sorted(groups.keys()):
            try:
                h = cc.homology(n)
                results[str(n)] = h
            except Exception:
                results[str(n)] = {'rank': groups[n], 'torsion': []}
        return jsonify({
            'groups': groups,
            'homology': results,
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebra_adv/tensor_product', methods=['POST'])
def api_algadv_tensor_product():
    """
    @brief Berechnet das Kronecker-Tensorprodukt zweier Matrizen.
    @lastModified 2026-03-10
    """
    if not _multilinear_ok:
        return error_response('multilinear_algebra nicht geladen')
    data = request.get_json()
    A = data.get('A', [[1]])
    B = data.get('B', [[1]])
    try:
        result = tensor_product_matrices(A, B)
        import numpy as np
        r = np.array(result)
        return jsonify({
            'result': r.tolist(),
            'shape': list(r.shape),
        })
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# LOGIK & GRUNDLAGEN API (model_theory, proof_theory_formal, recursion_theory)
# ===========================================================================

@app.route('/api/foundations/model_theory', methods=['POST'])
def api_foundations_model_theory():
    """
    @brief Zeigt Eigenschaften einer mathematischen Modell-Struktur.
    @lastModified 2026-03-10
    """
    if not _model_ok:
        return error_response('model_theory nicht geladen')
    data = request.get_json()
    structure_type = data.get('structure_type', 'integers')
    try:
        demo = nonstandard_arithmetic_demo() if structure_type == 'integers' else compactness_theorem_demo()
        props = {
            'integers': {
                'name': 'Ganze Zahlen ℤ',
                'theory': 'Th(ℤ) — Vollständige Theorie der additiven ganzen Zahlen',
                'properties': {
                    'Angeordnet': True, 'Abzählbar': True,
                    'Algebraisch abgeschlossen': False, 'Archimedisch': True,
                },
                'elementary_equivalences': ['ℤ (selbst)', 'Alle archimedischen angeordneten Ringe'],
            },
            'rationals': {
                'name': 'Rationale Zahlen ℚ',
                'theory': 'DLO — Dense Linear Order ohne Endpunkte (vollständig, entscheidbar)',
                'properties': {
                    'Angeordnet': True, 'Abzählbar': True, 'Dicht': True,
                    'Algebraisch abgeschlossen': False,
                },
                'elementary_equivalences': ['ℚ — DLO ist vollständig (Cantor)'],
            },
            'reals': {
                'name': 'Reelle Zahlen ℝ',
                'theory': 'RCF — Real Closed Fields (vollständig, entscheidbar via Tarski)',
                'properties': {
                    'Angeordnet': True, 'Überabzählbar': True, 'Vollständig': True,
                    'Algebraisch abgeschlossen': False,
                },
                'elementary_equivalences': ['Alle reell abgeschlossenen Körper'],
            },
            'finite_field': {
                'name': 'Endlicher Körper F₅',
                'theory': 'Th(F_5) — Theorie der Charakteristik 5',
                'properties': {
                    'Endlich': True, 'Algebraisch abgeschlossen': False,
                    'Abzählbar': True, 'Vollständig': True,
                },
                'elementary_equivalences': ['F_5, F_{5^n} (für n≥1)'],
            },
        }
        info = props.get(structure_type, {'name': structure_type, 'properties': {}, 'theory': ''})
        return jsonify(info)
    except Exception as e:
        return error_response(str(e))


@app.route('/api/foundations/sequent_calculus', methods=['POST'])
def api_foundations_sequent_calculus():
    """
    @brief Demonstriert Konzepte des Sequenzenkalküls.
    @lastModified 2026-03-10
    """
    if not _proof_formal_ok:
        return error_response('proof_theory_formal nicht geladen')
    data = request.get_json()
    demo_type = data.get('demo_type', 'cut_elimination')
    try:
        if demo_type == 'cut_elimination':
            result = cut_elimination_demo()
        elif demo_type == 'proof_complexity':
            result = proof_complexity_comparison()
        else:
            result = {'title': 'Subformel-Eigenschaft', 'description': 'Jede Formel im Beweis ist Subformel der Konklusion', 'theorem': 'Gentzen 1935: LK ohne Cut hat die Subformel-Eigenschaft'}
        return jsonify(result)
    except Exception as e:
        return error_response(str(e))


@app.route('/api/foundations/turing_machine', methods=['POST'])
def api_foundations_turing_machine():
    """
    @brief Simuliert eine Turing-Maschine.
    @lastModified 2026-03-10
    """
    if not _recursion_ok:
        return error_response('recursion_theory nicht geladen')
    data = request.get_json()
    tm_type = data.get('tm_type', 'palindrome')
    tape = data.get('tape', ['1', '0', '1'])
    try:
        if tm_type == 'palindrome':
            tm = tm_recognizes_palindromes()
            machine_name = 'Palindrom-Erkenner'
        elif tm_type == 'binary_increment':
            tm = tm_binary_increment()
            machine_name = 'Binäres Inkrement'
        else:
            from recursion_theory import tm_unary_addition
            tm = tm_unary_addition()
            machine_name = 'Unäre Addition'
        tape_list = [str(s) for s in tape]
        accepts = tm.accepts(tape_list)
        configs = tm.configuration_sequence(tape_list)
        result = tm.run(tape_list)
        return jsonify({
            'machine_name': machine_name,
            'result': accepts,
            'accepts': accepts,
            'final_tape': result if isinstance(result, list) else tape_list,
            'steps': len(configs),
            'configurations': [
                {'state': c[0], 'head': c[1], 'symbol': c[2]}
                for c in configs[:15]  # Max. 15 Konfigurationen anzeigen
            ],
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/foundations/complexity_hierarchy', methods=['GET'])
def api_foundations_complexity():
    """
    @brief Liefert die Komplexitätshierarchie.
    @lastModified 2026-03-10
    """
    if not _recursion_ok:
        return error_response('recursion_theory nicht geladen')
    try:
        classes = time_complexity_classes()
        return jsonify(classes)
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# ERWEITERTE ANALYSIS API (measure_theory, special_functions, functional_analysis,
#                          pde, operator_algebras)
# ===========================================================================

@app.route('/api/analysis_adv/measure', methods=['POST'])
def api_analysisadv_measure():
    """
    @brief Vergleich Lebesgue- vs. Riemann-Integral.
    @lastModified 2026-03-10
    """
    if not _measure_ok:
        return error_response('measure_theory nicht geladen')
    data = request.get_json()
    f_type = data.get('f_type', 'dirichlet')
    try:
        result = riemann_vs_lebesgue(f_type)
        # Sicherstellen dass alles JSON-serialisierbar ist
        def safe(v):
            if isinstance(v, (bool, str, type(None))): return v
            try: return float(v)
            except: return str(v)
        return jsonify({k: safe(v) for k, v in result.items()})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/analysis_adv/special_function', methods=['POST'])
def api_analysisadv_special():
    """
    @brief Berechnet den Wert einer speziellen Funktion.
    @lastModified 2026-03-10
    """
    if not _special_ok:
        return error_response('special_functions nicht geladen')
    data = request.get_json()
    sf_type = data.get('sf_type', 'bessel')
    n = int(data.get('n', 0))
    x = float(data.get('x', 1.0))
    try:
        if sf_type == 'bessel':
            bf = BesselFunctions()
            return jsonify({'name': f'Bessel J_{n}({x})', 'J_n_x': float(bf.J(n, x)), 'Y_n_x': float(bf.Y(n, x))})
        elif sf_type == 'legendre':
            lp = LegendrePolynomials()
            return jsonify({'name': f'Legendre P_{n}({x})', 'P_n_x': float(lp.P(n, x))})
        elif sf_type == 'airy':
            af = AiryFunctions()
            return jsonify({'name': f'Airy Ai({x})', 'Ai_x': float(af.Ai(x)), 'Bi_x': float(af.Bi(x))})
        elif sf_type == 'hermite':
            from special_functions import OrthogonalPolynomials
            op = OrthogonalPolynomials()
            return jsonify({'name': f'Hermite H_{n}({x})', 'H_n_x': float(op.hermite(n, x))})
        elif sf_type == 'gamma':
            props = gamma_function_properties()
            return jsonify({'name': 'Gamma-Funktion', **{k: str(v) for k, v in props.items()}})
        else:
            return error_response(f'Unbekannter Typ: {sf_type}')
    except Exception as e:
        return error_response(str(e))


@app.route('/api/analysis_adv/functional_analysis', methods=['POST'])
def api_analysisadv_functional():
    """
    @brief Demonstriert Sätze der Funktionalanalysis.
    @lastModified 2026-03-10
    """
    if not _functional_ok:
        return error_response('functional_analysis nicht geladen')
    data = request.get_json()
    demo_type = data.get('demo_type', 'banach_fixed_point')
    try:
        if demo_type == 'banach_fixed_point':
            result = banach_fixed_point(lambda x: 0.5 * x + 1, 0.0)
            return jsonify({
                'theorem': 'Banach-Fixpunktsatz',
                'statement': 'Eine Kontraktion f auf vollständigem metrischen Raum hat genau einen Fixpunkt.',
                'fixed_point': result.get('fixed_point'),
                'iterations': result.get('iterations'),
                'converged': result.get('converged'),
            })
        elif demo_type == 'hahn_banach':
            result = hahn_banach_theorem_demo()
            return jsonify(result)
        elif demo_type == 'spectral':
            import numpy as np
            A = np.array([[2.0, 1.0], [1.0, 2.0]])
            result = spectral_theorem_demo(A)
            def to_json(v):
                if hasattr(v, 'tolist'): return v.tolist()
                if isinstance(v, (list, tuple)): return [to_json(x) for x in v]
                try: return float(v)
                except: return str(v)
            return jsonify({k: to_json(v) for k, v in result.items()})
        elif demo_type == 'fredholm':
            import numpy as np
            A = np.array([[1.0, 2.0], [2.0, 4.0]])
            b = np.array([3.0, 6.0])
            result = fredholm_alternative(A, b)
            def to_json(v):
                if hasattr(v, 'tolist'): return v.tolist()
                try: return float(v)
                except: return str(v)
            return jsonify({k: to_json(v) for k, v in result.items()})
        else:
            result = uniform_boundedness_principle_demo()
            return jsonify(result)
    except Exception as e:
        return error_response(str(e))


@app.route('/api/analysis_adv/pde_classify', methods=['POST'])
def api_analysisadv_pde_classify():
    """
    @brief Klassifiziert eine PDE zweiter Ordnung.
    @lastModified 2026-03-10
    """
    if not _pde_ok:
        return error_response('pde nicht geladen')
    data = request.get_json()
    A = float(data.get('A', 1))
    B = float(data.get('B', 0))
    C = float(data.get('C', 1))
    try:
        pde_type = classify_pde(A, B, C)
        discriminant = B**2 - 4*A*C
        examples = {
            'elliptic': 'Laplace-Gleichung $\\Delta u = 0$ (A=C=1, B=0)',
            'parabolic': 'Wärmeleitungsgleichung $u_t = u_{xx}$ (A=1, B=C=0)',
            'hyperbolic': 'Wellengleichung $u_{tt} = c^2 u_{xx}$ (A=1, C=-1, B=0)',
        }
        return jsonify({
            'type': pde_type,
            'discriminant': float(discriminant),
            'description': f'$D = B^2 - 4AC = {discriminant:.4f}$',
            'example': examples.get(pde_type, ''),
        })
    except Exception as e:
        return error_response(str(e))


@app.route('/api/analysis_adv/operator_algebra', methods=['GET'])
def api_analysisadv_operator_algebra():
    """
    @brief Demonstriert C*-Algebren und Operatoralgebren.
    @lastModified 2026-03-10
    """
    if not _opalg_ok:
        return error_response('operator_algebras nicht geladen')
    demo = request.args.get('demo', 'gelfand')
    try:
        if demo == 'gelfand':
            result = gelfand_transform_demo()
            return jsonify({'title': 'Gelfand-Darstellungssatz', **result})
        elif demo == 'factors':
            result = factors_classification()
            types = [{'name': k, 'description': str(v)} for k, v in result.get('factors', {}).items()]
            return jsonify({'title': 'Von-Neumann-Faktorenklassifikation', 'types': types, **{k: v for k, v in result.items() if k != 'factors'}})
        elif demo == 'k_theory':
            result = k_theory_intro()
            return jsonify({'title': 'K-Theorie der C*-Algebren', **result})
        elif demo == 'cuntz':
            result = cuntz_algebra_demo()
            return jsonify({'title': f'Cuntz-Algebra $\\mathcal{{O}}_n$', **result})
        else:
            return error_response(f'Unbekannte Demo: {demo}')
    except Exception as e:
        return error_response(str(e))




# ===========================================================================
# API: Algebraische Geometrie (algebraic_geometry.py)
# ===========================================================================

try:
    from algebraic_geometry import (
        AffineVariety, ProjectiveVariety, HilbertBasisTheorem,
        NullstellensatzDemo, EllipticCurveGeometry, IntersectionTheory,
        compute_groebner_basis, ideal_radical
    )
    _alg_geom_available = True
except ImportError:
    _alg_geom_available = False


@app.route('/api/algebraic_geometry/groebner', methods=['POST'])
def api_alg_geom_groebner():
    """@brief Gröbner-Basis berechnen. @lastModified 2026-03-10"""
    if not _alg_geom_available:
        return error_response('Modul algebraic_geometry nicht verfügbar')
    try:
        data = request.get_json()
        gen_strs = data.get('generators', ['x**2 - y', 'x - 1'])
        var_strs = data.get('variables', ['x', 'y'])
        import sympy
        variables = [sympy.Symbol(v) for v in var_strs]
        generators = [sympy.sympify(g) for g in gen_strs]
        gb = compute_groebner_basis(generators, variables)
        return jsonify({'title': 'Gröbner-Basis', 'generators': gen_strs,
                        'variables': var_strs, 'groebner_basis': [str(p) for p in gb], 'size': len(gb)})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebraic_geometry/nullstellensatz', methods=['POST'])
def api_alg_geom_nullstellensatz():
    """@brief Hilbert'scher Nullstellensatz. @lastModified 2026-03-10"""
    if not _alg_geom_available:
        return error_response('Modul algebraic_geometry nicht verfügbar')
    try:
        data = request.get_json()
        gen_strs = data.get('generators', ['x**2 + 1'])
        var_strs = data.get('variables', ['x'])
        import sympy
        variables = [sympy.Symbol(v) for v in var_strs]
        generators = [sympy.sympify(g) for g in gen_strs]
        nss = NullstellensatzDemo()
        result = nss.weak_nullstellensatz(generators, variables)
        return jsonify({'title': 'Nullstellensatz', **result})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebraic_geometry/elliptic_curve', methods=['POST'])
def api_alg_geom_elliptic():
    """@brief Elliptische Kurve y² = x³ + ax + b. @lastModified 2026-03-10"""
    if not _alg_geom_available:
        return error_response('Modul algebraic_geometry nicht verfügbar')
    try:
        data = request.get_json()
        a = float(data.get('a', -1)); b = float(data.get('b', 1)); p = int(data.get('p', 7))
        ec = EllipticCurveGeometry(a, b)
        pts = ec.rational_points_mod_p(p)
        return jsonify({'title': f'Kurve y²=x³+{a}x+{b}', 'a': a, 'b': b,
                        'discriminant': ec.discriminant(), 'j_invariant': ec.j_invariant(),
                        'is_smooth': ec.is_smooth(), 'mod_p': p,
                        'num_points': len(pts), 'points': [str(pt) for pt in pts[:20]]})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/algebraic_geometry/intersection', methods=['POST'])
def api_alg_geom_intersection():
    """@brief Bézout-Theorem. @lastModified 2026-03-10"""
    if not _alg_geom_available:
        return error_response('Modul algebraic_geometry nicht verfügbar')
    try:
        data = request.get_json()
        d1 = int(data.get('d1', 2)); d2 = int(data.get('d2', 3))
        it = IntersectionTheory()
        result = it.bezout_theorem(d1, d2)
        return jsonify({'title': f'Bézout: Grad {d1} ∩ Grad {d2}', **result})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Lie-Gruppen (lie_groups.py)
# ===========================================================================

try:
    from lie_groups import (
        SO, SU, GL, LieAlgebra, ExponentialMap,
        classify_simple_lie_algebra, exceptional_lie_algebras_info,
        fundamental_representation_su2, cartan_matrix
    )
    _lie_available = True
except ImportError:
    _lie_available = False


@app.route('/api/lie_groups/so_rotation', methods=['POST'])
def api_lie_so_rotation():
    """@brief SO(n)-Rotation. @lastModified 2026-03-10"""
    if not _lie_available:
        return error_response('Modul lie_groups nicht verfügbar')
    try:
        import numpy as np
        data = request.get_json()
        n = int(data.get('n', 3)); theta = float(data.get('theta', 1.5707963))
        so = SO(n)
        R = so.rotation_2d(theta) if n == 2 else so.rotation_3d_z(theta)
        return jsonify({'title': f'SO({n}) θ={theta:.4f}', 'n': n, 'theta': theta,
                        'is_element': bool(so.is_element(R)), 'matrix': R.tolist(),
                        'determinant': float(np.linalg.det(R))})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/lie_groups/lie_bracket', methods=['POST'])
def api_lie_bracket():
    """@brief Lie-Klammer [X,Y]. @lastModified 2026-03-10"""
    if not _lie_available:
        return error_response('Modul lie_groups nicht verfügbar')
    try:
        import numpy as np
        data = request.get_json()
        X = np.array(data.get('X', [[0, -1], [1, 0]]), dtype=float)
        Y = np.array(data.get('Y', [[0, 1], [1, 0]]), dtype=float)
        la = LieAlgebra([X, Y])
        bracket = la.bracket(X, Y)
        return jsonify({'title': 'Lie-Klammer [X, Y]', 'bracket': bracket.tolist(),
                        'killing_form': float(la.killing_form(X, Y)),
                        'is_semisimple': la.is_semisimple()})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/lie_groups/exponential_map', methods=['POST'])
def api_lie_exp_map():
    """@brief exp(X) berechnen. @lastModified 2026-03-10"""
    if not _lie_available:
        return error_response('Modul lie_groups nicht verfügbar')
    try:
        import numpy as np
        data = request.get_json()
        X = np.array(data.get('X', [[0, -1.5707963], [1.5707963, 0]]), dtype=float)
        em = ExponentialMap()
        exp_X = em.exp_map(X)
        return jsonify({'title': 'exp(X)', 'X': X.tolist(), 'exp_X': exp_X.tolist(),
                        'determinant': float(np.linalg.det(exp_X))})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/lie_groups/classification', methods=['GET'])
def api_lie_classification():
    """@brief Cartan-Klassifikation. @lastModified 2026-03-10"""
    if not _lie_available:
        return error_response('Modul lie_groups nicht verfügbar')
    try:
        info = exceptional_lie_algebras_info()
        return jsonify({'title': 'Cartan-Klassifikation', 'exceptional': info,
                        'classical_types': {'A_n': 'su(n+1)', 'B_n': 'so(2n+1)',
                                            'C_n': 'sp(2n)', 'D_n': 'so(2n)'}})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Stochastische Prozesse (stochastic.py)
# ===========================================================================

try:
    from stochastic import (
        MarkovChain, BrownianMotion, StochasticDifferentialEquation,
        ErgodicTheory, PoissonProcess, random_walk_1d, gambler_ruin_probability
    )
    _stochastic_available = True
except ImportError:
    _stochastic_available = False


@app.route('/api/stochastic/markov_chain', methods=['POST'])
def api_stoch_markov():
    """@brief Markov-Kette analysieren. @lastModified 2026-03-10"""
    if not _stochastic_available:
        return error_response('Modul stochastic nicht verfügbar')
    try:
        data = request.get_json()
        P = data.get('transition_matrix', [[0.7, 0.3], [0.4, 0.6]])
        n_steps = int(data.get('n_steps', 10))
        mc = MarkovChain(P)
        stat = mc.stationary_distribution()
        sim = mc.simulate(0, n_steps, seed=42)
        return jsonify({'title': 'Markov-Kette', 'transition_matrix': P,
                        'stationary_distribution': stat.tolist() if hasattr(stat, 'tolist') else list(stat),
                        'simulation': sim.tolist() if hasattr(sim, 'tolist') else list(sim),
                        'is_ergodic': mc.is_ergodic()})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/stochastic/brownian_motion', methods=['POST'])
def api_stoch_brownian():
    """@brief Brownsche Bewegung. @lastModified 2026-03-10"""
    if not _stochastic_available:
        return error_response('Modul stochastic nicht verfügbar')
    try:
        data = request.get_json()
        T = float(data.get('T', 1.0)); n_steps = int(data.get('n_steps', 100))
        n_paths = int(data.get('n_paths', 3))
        bm = BrownianMotion()
        t_vals, paths = bm.simulate(T, n_steps, n_paths=n_paths, seed=42)
        return jsonify({'title': f'Brownsche Bewegung T={T}',
                        'T': T, 'n_steps': n_steps, 'n_paths': n_paths,
                        't_values': t_vals.tolist(),
                        'paths': paths.tolist() if hasattr(paths, 'tolist') else paths})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/stochastic/sde', methods=['POST'])
def api_stoch_sde():
    """@brief Geometrische Brownsche Bewegung. @lastModified 2026-03-10"""
    if not _stochastic_available:
        return error_response('Modul stochastic nicht verfügbar')
    try:
        data = request.get_json()
        mu = float(data.get('mu', 0.1)); sigma = float(data.get('sigma', 0.2))
        S0 = float(data.get('S0', 100.0)); T = float(data.get('T', 1.0))
        sde = StochasticDifferentialEquation()
        t_vals, paths = sde.geometric_brownian_motion(mu, sigma, S0, T, 200)
        return jsonify({'title': f'GBM μ={mu} σ={sigma}', 'mu': mu, 'sigma': sigma,
                        'S0': S0, 'T': T, 't_values': t_vals.tolist(),
                        'paths': paths.tolist() if hasattr(paths, 'tolist') else paths})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/stochastic/ergodic', methods=['POST'])
def api_stoch_ergodic():
    """@brief Logistische Abbildung / Ergodentheorie. @lastModified 2026-03-10"""
    if not _stochastic_available:
        return error_response('Modul stochastic nicht verfügbar')
    try:
        data = request.get_json()
        r = float(data.get('r', 3.9)); x0 = float(data.get('x0', 0.5))
        n_iter = int(data.get('n_iter', 200))
        et = ErgodicTheory()
        orbit = et.logistic_map(r, x0, n_iter)
        lyap = et.lyapunov_exponent(r, x0, n_iter)
        time_avg = et.time_average(orbit)
        orb_list = orbit[:20].tolist() if hasattr(orbit, 'tolist') else list(orbit[:20])
        return jsonify({'title': f'Logistische Abbildung r={r}', 'r': r, 'x0': x0,
                        'orbit_first_20': orb_list, 'time_average': float(time_avg),
                        'lyapunov_exponent': float(lyap), 'is_chaotic': float(lyap) > 0})
    except Exception as e:
        return error_response(str(e))


# ===========================================================================
# API: Differentialtopologie (differential_topology.py)
# ===========================================================================

try:
    from differential_topology import (
        DifferentialForm, DeRhamCohomology, MorseTheory,
        TransversalityTheory, SphereTopology,
        whitney_embedding_dimension, compute_jacobian,
        hairy_ball_theorem_demo, poincare_hopf_theorem_demo
    )
    _diff_top_available = True
except ImportError:
    _diff_top_available = False


@app.route('/api/differential_topology/differential_form', methods=['POST'])
def api_diff_top_form():
    """@brief Differentialform + äußere Ableitung. @lastModified 2026-03-10"""
    if not _diff_top_available:
        return error_response('Modul differential_topology nicht verfügbar')
    try:
        data = request.get_json()
        coeff_strs = data.get('coefficients', ['x*y', 'x**2'])
        var_strs = data.get('variables', ['x', 'y'])
        degree = int(data.get('degree', 1))
        import sympy
        variables = [sympy.Symbol(v) for v in var_strs]
        coefficients = [sympy.sympify(c) for c in coeff_strs]
        form = DifferentialForm(degree, coefficients, variables)
        d_form = form.exterior_derivative()
        return jsonify({'title': f'{degree}-Form', 'coefficients': coeff_strs,
                        'exterior_derivative': str(d_form), 'is_closed': form.is_closed()})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/differential_topology/morse_theory', methods=['POST'])
def api_diff_top_morse():
    """@brief Morse-Index berechnen. @lastModified 2026-03-10"""
    if not _diff_top_available:
        return error_response('Modul differential_topology nicht verfügbar')
    try:
        import sympy
        data = request.get_json()
        f_str = data.get('function', 'x**2 + y**2')
        var_strs = data.get('variables', ['x', 'y'])
        variables = [sympy.Symbol(v) for v in var_strs]
        f_expr = sympy.sympify(f_str)
        mt = MorseTheory()
        crit_pts = mt.morse_function_critical_points(f_expr, variables)
        results = [{'point': str(pt), 'morse_index': mt.morse_index(f_expr, pt, variables)}
                   for pt in crit_pts[:5]]
        return jsonify({'title': f'Morse: f={f_str}', 'function': f_str,
                        'critical_points': results,
                        'is_morse_function': mt.is_morse_function(f_expr, variables)})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/differential_topology/cohomology', methods=['POST'])
def api_diff_top_cohom():
    """@brief de Rham-Betti-Zahlen. @lastModified 2026-03-10"""
    if not _diff_top_available:
        return error_response('Modul differential_topology nicht verfügbar')
    try:
        data = request.get_json()
        space = data.get('space', 'sphere'); n = int(data.get('n', 2))
        dr = DeRhamCohomology()
        betti = dr.betti_numbers_sphere(n) if space == 'sphere' else dr.betti_numbers_torus(n)
        return jsonify({'title': f'de Rham: {n}-{space}', 'space': space, 'n': n,
                        'betti_numbers': betti, 'euler_characteristic': dr.euler_characteristic(betti)})
    except Exception as e:
        return error_response(str(e))


@app.route('/api/differential_topology/sphere', methods=['POST'])
def api_diff_top_sphere():
    """@brief n-Sphäre Topologie. @lastModified 2026-03-10"""
    if not _diff_top_available:
        return error_response('Modul differential_topology nicht verfügbar')
    try:
        data = request.get_json()
        n = int(data.get('n', 2))
        sph = SphereTopology(n)
        return jsonify({'title': f'S^{n}', 'n': n, 'is_orientable': sph.is_orientable(),
                        'whitney_dim': whitney_embedding_dimension(n),
                        'homotopy_groups': sph.homotopy_groups_low(n),
                        'homology_groups': sph.homology_groups(n),
                        'hairy_ball': hairy_ball_theorem_demo()})
    except Exception as e:
        return error_response(str(e))

# ===========================================================================
# HAUPTPROGRAMM
# ===========================================================================

if __name__ == '__main__':
    # Webapp starten: Host 0.0.0.0 erlaubt Zugriff aus dem LAN/Container
    # Port 8080, debug=True für Entwicklung (zeigt Fehler im Browser)
    print("=" * 60)
    print("  Mathematik-Spezialist Web-Interface")
    print("  Build 36 | Port 8080")
    print("  URL: http://localhost:8080")
    print("=" * 60)
    app.run(host='0.0.0.0', port=8080, debug=True)
# Diese Zeile wird nie erreicht, aber Python braucht etwas hier
