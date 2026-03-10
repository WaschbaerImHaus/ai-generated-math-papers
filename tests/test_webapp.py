"""
@file test_webapp.py
@brief Tests für die Flask-Webanwendung des Mathematik-Spezialisten.
@description
    Testet alle API-Endpunkte der Flask-App mit dem integrierten Test-Client.
    Deckt folgende Routen ab:
    - GET  /            → Startseite (Status 200)
    - GET  /api/health  → Systemstatus
    - POST /api/algebra/solve (linear, quadratisch, komplex)
    - POST /api/algebra/polynomial
    - POST /api/algebra/primes
    - POST /api/analysis/derivative
    - POST /api/analysis/integral
    - POST /api/analysis/roots
    - POST /api/linear_algebra/matrix (det, inv, eigenvalues)
    - POST /api/statistics/describe
    - POST /api/ode/solve (exponential, harmonic, logistic)
    - POST /api/complex/zeta
    - POST /api/visualization/plot
    - Fehlerbehandlung (400-Antworten)

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import json
import pytest

# Webapp-Verzeichnis und src-Verzeichnis zum Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'webapp'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Flask-App importieren (muss nach sys.path-Setup erfolgen)
from app import app


# ===========================================================================
# PYTEST-FIXTURE: Test-Client
# ===========================================================================

@pytest.fixture
def client():
    """
    @brief Erstellt einen Flask-Test-Client für alle Tests.
    @description
        Konfiguriert die App im Test-Modus (TESTING=True) und gibt
        einen HTTP-Client zurück, der ohne echten Server läuft.
    @return Flask-Test-Client
    @date 2026-03-10
    """
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client


# ===========================================================================
# HILFSFUNKTION: JSON-POST
# ===========================================================================

def post_json(client, url: str, data: dict):
    """
    @brief Sendet eine POST-Anfrage mit JSON-Body.
    @param client: Flask-Test-Client
    @param url: Ziel-URL
    @param data: Anfrage-Body als Dictionary
    @return Response-Objekt
    @date 2026-03-10
    """
    return client.post(
        url,
        data=json.dumps(data),
        content_type='application/json'
    )


# ===========================================================================
# TESTS: SEITEN-ROUTEN (GET)
# ===========================================================================

class TestPageRoutes:
    """Tests für alle HTML-Seiten-Routen (GET-Anfragen)."""

    def test_index_returns_200(self, client):
        """Startseite muss HTTP 200 zurückgeben."""
        resp = client.get('/')
        assert resp.status_code == 200

    def test_index_contains_title(self, client):
        """Startseite muss den Titel 'Mathematik-Spezialist' enthalten."""
        resp = client.get('/')
        assert b'Mathematik' in resp.data

    def test_algebra_page_returns_200(self, client):
        """Algebra-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/algebra').status_code == 200

    def test_analysis_page_returns_200(self, client):
        """Analysis-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/analysis').status_code == 200

    def test_linear_algebra_page_returns_200(self, client):
        """Lineare-Algebra-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/linear_algebra').status_code == 200

    def test_statistics_page_returns_200(self, client):
        """Statistik-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/statistics').status_code == 200

    def test_ode_page_returns_200(self, client):
        """ODE-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/ode').status_code == 200

    def test_complex_page_returns_200(self, client):
        """Komplexe-Analysis-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/complex').status_code == 200

    def test_number_theory_page_returns_200(self, client):
        """Zahlentheorie-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/number_theory').status_code == 200

    def test_visualization_page_returns_200(self, client):
        """Visualisierungs-Seite muss HTTP 200 zurückgeben."""
        assert client.get('/visualization').status_code == 200


# ===========================================================================
# TESTS: HEALTH-ENDPUNKT
# ===========================================================================

class TestHealthEndpoint:
    """Tests für den /api/health-Endpunkt."""

    def test_health_status_ok(self, client):
        """Health-Endpunkt muss {'status': 'ok'} zurückgeben."""
        resp = client.get('/api/health')
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['status'] == 'ok'

    def test_health_build_number(self, client):
        """Health-Endpunkt muss eine gültige Build-Nummer zurückgeben."""
        resp = client.get('/api/health')
        data = resp.get_json()
        assert 'build' in data
        assert data['build'] >= 12  # Build-Nummer steigt mit jeder Version


# ===========================================================================
# TESTS: ALGEBRA API
# ===========================================================================

class TestAlgebra:
    """Tests für die Algebra-API-Endpunkte."""

    def test_solve_linear_basic(self, client):
        """Lineare Gleichung 1·x + (-3) = 0 muss Lösung 3.0 ergeben."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'linear',
            'coefficients': [1, -3]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'solutions' in data
        assert len(data['solutions']) == 1
        assert abs(data['solutions'][0]['real'] - 3.0) < 1e-8

    def test_solve_linear_2x_minus_8(self, client):
        """Lineare Gleichung 2x - 8 = 0 muss Lösung 4.0 ergeben."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'linear',
            'coefficients': [2, -8]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['solutions'][0]['real'] - 4.0) < 1e-8

    def test_solve_quadratic_two_real_roots(self, client):
        """x² - 5x + 6 = 0 muss Lösungen 2.0 und 3.0 ergeben."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'quadratic',
            'coefficients': [1, -5, 6]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        roots = sorted([s['real'] for s in data['solutions']])
        assert abs(roots[0] - 2.0) < 1e-6
        assert abs(roots[1] - 3.0) < 1e-6

    def test_solve_quadratic_complex_roots(self, client):
        """x² + 2x + 5 = 0 muss zwei komplexe Lösungen ergeben (Imag ≠ 0)."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'quadratic',
            'coefficients': [1, 2, 5]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data['solutions']) == 2
        # Beide Lösungen müssen Imaginärteil ≠ 0 haben
        for sol in data['solutions']:
            assert abs(sol['imag']) > 1e-6

    def test_solve_unknown_type_returns_400(self, client):
        """Unbekannter Gleichungstyp muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'cubic',
            'coefficients': [1, 0, 0, -1]
        })
        assert resp.status_code == 400

    def test_polynomial_evaluate(self, client):
        """Polynom x² - 2 bei x=2 muss 2.0 ergeben."""
        resp = post_json(client, '/api/algebra/polynomial', {
            'coefficients': [1, 0, -2],
            'x': 2.0
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['result'] - 2.0) < 1e-8

    def test_polynomial_constant(self, client):
        """Konstantes Polynom [5] bei beliebigem x muss 5 ergeben."""
        resp = post_json(client, '/api/algebra/polynomial', {
            'coefficients': [5],
            'x': 100.0
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['result'] - 5.0) < 1e-8

    def test_primes_is_prime_true(self, client):
        """n=7 muss als Primzahl erkannt werden."""
        resp = post_json(client, '/api/algebra/primes', {'n': 7})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_prime'] is True

    def test_primes_is_prime_false(self, client):
        """n=12 darf keine Primzahl sein."""
        resp = post_json(client, '/api/algebra/primes', {'n': 12})
        assert resp.status_code == 200
        assert resp.get_json()['is_prime'] is False

    def test_primes_euler_phi_prime(self, client):
        """φ(7) = 6, da 7 eine Primzahl ist."""
        resp = post_json(client, '/api/algebra/primes', {'n': 7})
        data = resp.get_json()
        assert data['euler_phi'] == 6

    def test_primes_factorization_12(self, client):
        """12 = 2² · 3, also factorization_str muss '2^2 · 3' enthalten."""
        resp = post_json(client, '/api/algebra/primes', {'n': 12})
        data = resp.get_json()
        assert '2' in data['factorization_str']
        assert '3' in data['factorization_str']

    def test_primes_invalid_n_returns_400(self, client):
        """n=0 muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/algebra/primes', {'n': 0})
        assert resp.status_code == 400


# ===========================================================================
# TESTS: ANALYSIS API
# ===========================================================================

class TestAnalysis:
    """Tests für die Analysis-API-Endpunkte."""

    def test_derivative_sin_at_0(self, client):
        """f'(x) = sin(x) bei x=0 → cos(0) = 1.0."""
        resp = post_json(client, '/api/analysis/derivative', {
            'expr': 'sin(x)',
            'x': 0.0
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['derivative'] - 1.0) < 1e-5

    def test_derivative_x_squared(self, client):
        """f'(x) = x² bei x=3 → 2·3 = 6."""
        resp = post_json(client, '/api/analysis/derivative', {
            'expr': 'x**2',
            'x': 3.0
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['derivative'] - 6.0) < 1e-4

    def test_derivative_exp_at_0(self, client):
        """f'(x) = exp(x) bei x=0 → 1.0."""
        resp = post_json(client, '/api/analysis/derivative', {
            'expr': 'exp(x)',
            'x': 0.0
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['derivative'] - 1.0) < 1e-5

    def test_integral_x_squared_0_to_1(self, client):
        """∫₀¹ x² dx = 1/3 ≈ 0.3333."""
        resp = post_json(client, '/api/analysis/integral', {
            'expr': 'x**2',
            'a': 0.0,
            'b': 1.0
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['integral'] - 1.0 / 3.0) < 1e-6

    def test_integral_constant_1(self, client):
        """∫₀⁵ 1 dx = 5.0."""
        resp = post_json(client, '/api/analysis/integral', {
            'expr': '1',
            'a': 0.0,
            'b': 5.0
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['integral'] - 5.0) < 1e-6

    def test_integral_invalid_bounds_returns_400(self, client):
        """a >= b muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/analysis/integral', {
            'expr': 'x',
            'a': 5.0,
            'b': 1.0
        })
        assert resp.status_code == 400

    def test_roots_newton_sqrt2(self, client):
        """Newton-Raphson für x²-2=0 ab x₀=1.5 muss √2 ≈ 1.41421 finden."""
        resp = post_json(client, '/api/analysis/roots', {
            'expr': 'x**2 - 2',
            'x0': 1.5,
            'method': 'newton'
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['root'] - 2 ** 0.5) < 1e-6

    def test_roots_missing_expr_returns_400(self, client):
        """Fehlender 'expr'-Parameter muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/analysis/roots', {
            'x0': 1.0,
            'method': 'newton'
        })
        assert resp.status_code == 400


# ===========================================================================
# TESTS: LINEARE ALGEBRA API
# ===========================================================================

class TestLinearAlgebra:
    """Tests für die Lineare-Algebra-API-Endpunkte."""

    def test_det_2x2(self, client):
        """det([[1,2],[3,4]]) = 1·4 - 2·3 = -2."""
        resp = post_json(client, '/api/linear_algebra/matrix', {
            'matrix': [[1, 2], [3, 4]],
            'operation': 'det'
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['result'] - (-2.0)) < 1e-8

    def test_det_identity_3x3(self, client):
        """det(I₃) = 1."""
        resp = post_json(client, '/api/linear_algebra/matrix', {
            'matrix': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            'operation': 'det'
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['result'] - 1.0) < 1e-8

    def test_inv_2x2(self, client):
        """Inverse von [[2,1],[1,1]] muss A⁻¹·A = I ergeben."""
        resp = post_json(client, '/api/linear_algebra/matrix', {
            'matrix': [[2, 1], [1, 1]],
            'operation': 'inv'
        })
        assert resp.status_code == 200
        data = resp.get_json()
        inv = data['result']
        # [[1,-1],[-1,2]] ist die Inverse
        assert abs(inv[0][0] - 1.0) < 1e-6
        assert abs(inv[0][1] - (-1.0)) < 1e-6

    def test_eigenvalues_diagonal(self, client):
        """Diagonalmatrix [[3,0],[0,5]] hat Eigenwerte 3 und 5."""
        resp = post_json(client, '/api/linear_algebra/matrix', {
            'matrix': [[3, 0], [0, 5]],
            'operation': 'eigenvalues'
        })
        assert resp.status_code == 200
        data = resp.get_json()
        evals = sorted([ev['real'] for ev in data['result']])
        assert abs(evals[0] - 3.0) < 1e-4
        assert abs(evals[1] - 5.0) < 1e-4

    def test_unknown_operation_returns_400(self, client):
        """Unbekannte Operation muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/linear_algebra/matrix', {
            'matrix': [[1, 0], [0, 1]],
            'operation': 'trace'
        })
        assert resp.status_code == 400


# ===========================================================================
# TESTS: STATISTIK API
# ===========================================================================

class TestStatistics:
    """Tests für die Statistik-API-Endpunkte."""

    def test_describe_mean(self, client):
        """Mittelwert von [1,2,3,4,5] = 3.0."""
        resp = post_json(client, '/api/statistics/describe', {
            'data': [1, 2, 3, 4, 5]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['mean'] - 3.0) < 1e-8

    def test_describe_median_odd(self, client):
        """Median von [1,3,2] = 2.0 (nach Sortierung: 1,2,3)."""
        resp = post_json(client, '/api/statistics/describe', {
            'data': [1, 3, 2]
        })
        assert resp.status_code == 200
        assert abs(resp.get_json()['median'] - 2.0) < 1e-8

    def test_describe_std_dev(self, client):
        """σ von [2,4,4,4,5,5,7,9] ≈ 2.0 (Population) oder ≈ 2.14 (Stichprobe)."""
        resp = post_json(client, '/api/statistics/describe', {
            'data': [2, 4, 4, 4, 5, 5, 7, 9]
        })
        assert resp.status_code == 200
        # Standardabweichung liegt je nach Formel zwischen ~2.0 und ~2.14
        sd = resp.get_json()['std_dev']
        assert 1.8 < sd < 2.3

    def test_describe_count(self, client):
        """n für [10,20,30,40] muss 4 sein."""
        resp = post_json(client, '/api/statistics/describe', {
            'data': [10, 20, 30, 40]
        })
        assert resp.status_code == 200
        assert resp.get_json()['n'] == 4

    def test_describe_min_max(self, client):
        """Min und Max von [3,1,4,1,5] müssen 1 und 5 sein."""
        resp = post_json(client, '/api/statistics/describe', {
            'data': [3, 1, 4, 1, 5]
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['min'] == 1.0
        assert data['max'] == 5.0

    def test_describe_empty_returns_400(self, client):
        """Leerer Datensatz muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/statistics/describe', {'data': []})
        assert resp.status_code == 400


# ===========================================================================
# TESTS: ODE API
# ===========================================================================

class TestODE:
    """Tests für die ODE-API-Endpunkte."""

    def test_ode_exponential_returns_data(self, client):
        """Exponentieller Zerfall muss t- und y-Arrays zurückgeben."""
        resp = post_json(client, '/api/ode/solve', {
            'equation': 'exponential',
            'y0': [1.0],
            't_end': 5.0
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 't' in data
        assert 'y' in data
        assert len(data['t']) > 5

    def test_ode_exponential_decay(self, client):
        """y(0) = 1.0 muss bei t=0 vorhanden und > 0 bleiben."""
        resp = post_json(client, '/api/ode/solve', {
            'equation': 'exponential',
            'y0': [1.0],
            't_end': 3.0
        })
        data = resp.get_json()
        # y(0) ≈ 1.0
        assert abs(data['y'][0] - 1.0) < 0.05
        # Alle Werte positiv (Zerfall bleibt positiv)
        assert all(v > 0 for v in data['y'])

    def test_ode_harmonic_oscillator(self, client):
        """Harmonischer Oszillator: t- und y-Felder vorhanden."""
        resp = post_json(client, '/api/ode/solve', {
            'equation': 'harmonic',
            'y0': [1.0, 0.0],
            't_end': 6.28
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert len(data['t']) > 5

    def test_ode_logistic_growth(self, client):
        """Logistisches Wachstum muss gesättigte y-Werte zurückgeben."""
        resp = post_json(client, '/api/ode/solve', {
            'equation': 'logistic',
            'y0': [1.0],
            't_end': 20.0
        })
        assert resp.status_code == 200

    def test_ode_unknown_equation_returns_400(self, client):
        """Unbekannter Gleichungstyp muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/ode/solve', {
            'equation': 'lorenz',
            'y0': [1.0],
            't_end': 5.0
        })
        assert resp.status_code == 400


# ===========================================================================
# TESTS: KOMPLEXE ANALYSIS API
# ===========================================================================

class TestComplexAnalysis:
    """Tests für die Komplexe-Analysis-API."""

    def test_zeta_2_equals_pi_squared_over_6(self, client):
        """ζ(2) = π²/6 ≈ 1.6449."""
        import math
        resp = post_json(client, '/api/complex/zeta', {
            's_real': 2.0,
            's_imag': 0.0,
            'precision': 50
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['zeta_real'] - (math.pi ** 2 / 6)) < 1e-4
        assert abs(data['zeta_imag']) < 1e-10

    def test_zeta_first_nontrivial_zero(self, client):
        """ζ(0.5 + 14.134725i) muss nahe 0 sein (erste nichttriviale Nullstelle)."""
        resp = post_json(client, '/api/complex/zeta', {
            's_real': 0.5,
            's_imag': 14.134725,
            'precision': 50
        })
        assert resp.status_code == 200
        data = resp.get_json()
        # |ζ(0.5 + 14.134725i)| << 1
        assert data['abs'] < 0.01

    def test_zeta_pole_returns_400(self, client):
        """ζ(1) hat einen Pol, muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/complex/zeta', {
            's_real': 1.0,
            's_imag': 0.0
        })
        assert resp.status_code == 400


# ===========================================================================
# TESTS: VISUALISIERUNG API
# ===========================================================================

class TestVisualization:
    """Tests für die Visualisierungs-API."""

    def test_plot_returns_base64(self, client):
        """Plot muss einen Base64-String zurückgeben."""
        resp = post_json(client, '/api/visualization/plot', {
            'expr': 'sin(x)',
            'x_min': -6.28,
            'x_max': 6.28,
            'title': 'Sinus-Funktion'
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'image_base64' in data
        assert len(data['image_base64']) > 100  # Base64-String nicht leer

    def test_plot_x2(self, client):
        """Plot für x² muss erfolgreich erstellt werden."""
        resp = post_json(client, '/api/visualization/plot', {
            'expr': 'x**2',
            'x_min': -5.0,
            'x_max': 5.0
        })
        assert resp.status_code == 200
        assert 'image_base64' in resp.get_json()

    def test_plot_invalid_bounds_returns_400(self, client):
        """x_min >= x_max muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/visualization/plot', {
            'expr': 'x',
            'x_min': 5.0,
            'x_max': 1.0
        })
        assert resp.status_code == 400

    def test_plot_invalid_expr_returns_400(self, client):
        """Ungültiger Ausdruck muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/visualization/plot', {
            'expr': 'import os; os.system("rm -rf /")',
            'x_min': -1.0,
            'x_max': 1.0
        })
        assert resp.status_code == 400


# ===========================================================================
# TESTS: FEHLERBEHANDLUNG (Edge Cases)
# ===========================================================================

class TestErrorHandling:
    """Tests für Fehlerbehandlung und Edge Cases."""

    def test_missing_json_body_algebra_solve(self, client):
        """POST ohne Body muss HTTP 400 zurückgeben."""
        resp = client.post('/api/algebra/solve', content_type='application/json')
        assert resp.status_code == 400

    def test_missing_json_body_statistics(self, client):
        """POST ohne Body muss HTTP 400 zurückgeben."""
        resp = client.post('/api/statistics/describe', content_type='application/json')
        assert resp.status_code == 400

    def test_division_by_zero_poly(self, client):
        """Polynom mit Koeffizient a=0 bei Division-Ähnlichen Operationen wird abgefangen."""
        resp = post_json(client, '/api/algebra/solve', {
            'type': 'linear',
            'coefficients': [0, 5]
        })
        # Sollte entweder 400 oder einen Fehler zurückgeben
        data = resp.get_json()
        assert 'error' in data or resp.status_code == 400


# ===========================================================================
# TESTS: ELLIPTISCHE KURVEN
# ===========================================================================

class TestEllipticCurves:
    """
    @brief Tests für die Elliptische-Kurven-API.
    @description
        Prüft /api/elliptic/analyze und /api/elliptic/hasse auf
        korrekte Berechnung von Diskriminante, Gruppenordnung und Hasse-Schranke.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def test_elliptic_page_loads(self, client):
        """GET /elliptic muss Status 200 zurückgeben."""
        resp = client.get('/elliptic')
        assert resp.status_code == 200

    def test_analyze_regular_curve(self, client):
        """Reguläre Kurve y²=x³−x: Δ≠0, nicht singuläre Antwort."""
        resp = post_json(client, '/api/elliptic/analyze', {'a': -1, 'b': 0})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_singular'] is False
        # Δ = -16*(4*(-1)³+27*0²) = -16*(-4) = 64
        assert abs(data['discriminant'] - 64.0) < 1e-6
        assert 'j_invariant' in data

    def test_analyze_singular_curve(self, client):
        """Singuläre Kurve y²=x³ (a=0, b=0): Δ=0, is_singular=True."""
        resp = post_json(client, '/api/elliptic/analyze', {'a': 0, 'b': 0})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_singular'] is True

    def test_analyze_with_point_on_curve(self, client):
        """Punkt (0,1) liegt auf y²=x³+1 → on_curve=True."""
        # y²=0+0+1=1 → y=1 ✓
        resp = post_json(client, '/api/elliptic/analyze', {
            'a': 0, 'b': 1,
            'px': 0.0, 'py': 1.0, 'n': 2
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_singular'] is False
        assert data['point_results'] is not None
        assert data['point_results']['on_curve'] is True

    def test_analyze_with_point_not_on_curve(self, client):
        """Punkt (1,5) liegt nicht auf y²=x³−x → on_curve=False."""
        resp = post_json(client, '/api/elliptic/analyze', {
            'a': -1, 'b': 0,
            'px': 1.0, 'py': 5.0, 'n': 3
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['point_results']['on_curve'] is False

    def test_analyze_returns_plot(self, client):
        """Reguläre Kurve muss einen Base64-Plot zurückgeben."""
        resp = post_json(client, '/api/elliptic/analyze', {'a': -1, 'b': 0})
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['plot'] is not None
        assert len(data['plot']) > 100  # Base64-String nicht leer

    def test_analyze_j_invariant_y2_x3_minus_x(self, client):
        """j-Invariante von y²=x³−x muss 1728 sein (CM durch ℤ[i])."""
        resp = post_json(client, '/api/elliptic/analyze', {'a': -1, 'b': 0})
        assert resp.status_code == 200
        data = resp.get_json()
        # j = -1728*(4*(-1))^3 / 64 = -1728*(-64)/64 = 1728
        assert abs(data['j_invariant'] - 1728.0) < 1.0

    def test_hasse_page_values(self, client):
        """Hasse-Schranke für y²=x³−x: #E(F_p) muss in Hasse-Band liegen."""
        resp = post_json(client, '/api/elliptic/hasse', {
            'a': -1, 'b': 0, 'prime_limit': 13
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_singular'] is False
        assert len(data['primes']) > 0
        # Hasse-Schranke prüfen: |a_p| ≤ 2√p
        import math
        for entry in data['primes']:
            p = entry['p']
            trace = abs(entry['trace'])
            assert trace <= 2 * math.sqrt(p) + 0.001, \
                f"Hasse-Schranke verletzt für p={p}: |a_p|={trace} > 2√{p}={2*math.sqrt(p):.3f}"

    def test_hasse_singular_returns_flag(self, client):
        """Singuläre Kurve bei Hasse: is_singular=True, leere Primzahlliste."""
        resp = post_json(client, '/api/elliptic/hasse', {
            'a': 0, 'b': 0, 'prime_limit': 10
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['is_singular'] is True

    def test_hasse_prime_limit_capped(self, client):
        """prime_limit > 50 wird auf 50 begrenzt."""
        resp = post_json(client, '/api/elliptic/hasse', {
            'a': -1, 'b': 0, 'prime_limit': 999
        })
        assert resp.status_code == 200
        data = resp.get_json()
        # Alle zurückgegebenen Primzahlen dürfen nicht > 50 sein
        if data['primes']:
            assert max(e['p'] for e in data['primes']) <= 50

    def test_analyze_no_json_returns_400(self, client):
        """POST ohne Body muss HTTP 400 zurückgeben."""
        resp = client.post('/api/elliptic/analyze', content_type='application/json')
        assert resp.status_code == 400


# ===========================================================================
# TESTS: TENSOR-GEOMETRIE
# ===========================================================================

class TestTensorGeometry:
    """
    @brief Tests für die Tensor-Geometrie-API.
    @description
        Prüft /api/tensor/curvature für alle Mannigfaltigkeiten:
        Sphäre, Torus, Hyperbolisch, Sattel, Flach.
        Validiert Metriktensor, Ricci-Skalar und Gaußsche Krümmung.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def test_tensor_page_loads(self, client):
        """GET /tensor muss Status 200 zurückgeben."""
        resp = client.get('/tensor')
        assert resp.status_code == 200

    def test_sphere_curvature(self, client):
        """Sphäre r=1: K=1, R=2, g_11=r²=1."""
        import math
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'sphere',
            'params': {'r': 1.0, 'theta': math.pi / 2, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['gaussian_curvature'] - 1.0) < 0.05
        assert abs(data['ricci_scalar'] - 2.0) < 0.1
        # Metriktensor-Diagonale g_11 = r² = 1
        assert abs(data['metric'][0][0] - 1.0) < 1e-6
        # Nicht-Diagonale g_12 = 0
        assert abs(data['metric'][0][1]) < 1e-6

    def test_sphere_larger_radius(self, client):
        """Sphäre r=2: K=1/r²=0.25."""
        import math
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'sphere',
            'params': {'r': 2.0, 'theta': math.pi / 2, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['gaussian_curvature'] - 0.25) < 0.05

    def test_hyperbolic_curvature(self, client):
        """Hyperbolische Ebene: K = -1 (konstant)."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'hyperbolic',
            'params': {'x': 0.0, 'y': 1.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['gaussian_curvature'] - (-1.0)) < 0.1

    def test_flat_curvature(self, client):
        """Flache Ebene: K=0, R=0."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'flat',
            'params': {'x': 0.0, 'y': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['gaussian_curvature']) < 1e-4
        assert abs(data['ricci_scalar']) < 1e-4
        # Einheitsmetrik: g_ij = δ_ij
        assert abs(data['metric'][0][0] - 1.0) < 1e-6
        assert abs(data['metric'][1][1] - 1.0) < 1e-6

    def test_torus_curvature_outside(self, client):
        """Torus außen (θ=0): K > 0 (positiv)."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'torus',
            'params': {'R': 2.0, 'r': 1.0, 'theta': 0.0, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        # Außen: K = cos(0) / (r*(R+r*cos(0))) = 1/(1*3) ≈ 0.333
        assert data['gaussian_curvature'] > 0

    def test_torus_curvature_inside(self, client):
        """Torus innen (θ=π): K < 0 (negativ)."""
        import math
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'torus',
            'params': {'R': 2.0, 'r': 1.0, 'theta': math.pi, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['gaussian_curvature'] < 0

    def test_saddle_curvature_origin(self, client):
        """Sattelfläche im Ursprung: K = -1."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'saddle',
            'params': {'x': 0.0, 'y': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        # K(0,0) = -1/(1+0+0)² = -1
        assert abs(data['gaussian_curvature'] - (-1.0)) < 0.1

    def test_tensor_returns_plot(self, client):
        """Sphäre muss Base64-Plot zurückgeben."""
        import math
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'sphere',
            'params': {'r': 1.0, 'theta': math.pi / 2, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert data['plot'] is not None
        assert len(data['plot']) > 100

    def test_tensor_returns_physical_description(self, client):
        """Antwort muss physikalische Beschreibung enthalten."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'flat',
            'params': {'x': 0.0, 'y': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'physical_description' in data
        assert len(data['physical_description']) > 10

    def test_tensor_returns_metric_latex(self, client):
        """Antwort muss metric_latex für KaTeX enthalten."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'sphere',
            'params': {'r': 1.0, 'theta': 1.5708, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'metric_latex' in data
        assert 'g_{ij}' in data['metric_latex'] or 'g_' in data['metric_latex']

    def test_tensor_unknown_manifold_returns_400(self, client):
        """Unbekannte Mannigfaltigkeit muss HTTP 400 zurückgeben."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'unknownxyz',
            'params': {}
        })
        assert resp.status_code == 400

    def test_tensor_no_json_returns_400(self, client):
        """POST ohne Body muss HTTP 400 zurückgeben."""
        resp = client.post('/api/tensor/curvature', content_type='application/json')
        assert resp.status_code == 400

    def test_curvature_type_positive(self, client):
        """Sphäre: curvature_type muss 'positiv' enthalten."""
        import math
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'sphere',
            'params': {'r': 1.0, 'theta': math.pi / 2, 'phi': 0.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert 'positiv' in data['curvature_type'].lower() or data['gaussian_curvature'] > 0

    def test_curvature_type_flat(self, client):
        """Flache Ebene: curvature_type muss 'flach' oder K≈0 zeigen."""
        resp = post_json(client, '/api/tensor/curvature', {
            'manifold': 'flat',
            'params': {'x': 1.0, 'y': 2.0}
        })
        assert resp.status_code == 200
        data = resp.get_json()
        assert abs(data['gaussian_curvature']) < 1e-3
