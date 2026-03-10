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


# ===========================================================================
# ROUTE: SYSTEMSTATUS
# ===========================================================================

@app.route('/api/health', methods=['GET'])
def health():
    """
    @brief Gibt den Systemstatus zurück.
    @return JSON {"status": "ok", "build": 11}
    @date 2026-03-10
    """
    return jsonify({"status": "ok", "build": 13})


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
            # solve_linear gibt Zahl oder None zurück
            solutions = [result] if result is not None else []
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
# HAUPTPROGRAMM
# ===========================================================================

if __name__ == '__main__':
    # Webapp starten: Host 0.0.0.0 erlaubt Zugriff aus dem LAN/Container
    # Port 8080, debug=True für Entwicklung (zeigt Fehler im Browser)
    print("=" * 60)
    print("  Mathematik-Spezialist Web-Interface")
    print("  Build 15 | Port 8080")
    print("  URL: http://localhost:8080")
    print("=" * 60)
    app.run(host='0.0.0.0', port=8080, debug=True)
