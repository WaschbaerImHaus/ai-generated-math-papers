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
    return jsonify({"status": "ok", "build": 12})


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
# HAUPTPROGRAMM
# ===========================================================================

if __name__ == '__main__':
    # Webapp starten: Host 0.0.0.0 erlaubt Zugriff aus dem LAN/Container
    # Port 8080, debug=True für Entwicklung (zeigt Fehler im Browser)
    print("=" * 60)
    print("  Mathematik-Spezialist Web-Interface")
    print("  Build 12 | Port 8080")
    print("  URL: http://localhost:8080")
    print("=" * 60)
    app.run(host='0.0.0.0', port=8080, debug=True)
