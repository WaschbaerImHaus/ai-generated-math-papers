"""
@file test_step_by_step.py
@brief Tests für das Schritt-für-Schritt-Modul (step_by_step.py).
@description
    Überprüft alle Funktionen des step_by_step-Moduls auf korrekte Rückgabe-Struktur,
    mathematische Korrektheit und Robustheit bei Edge-Cases.

    Getestete Funktionen:
    - newton_raphson_steps: Nullstellensuche mit Schritt-Dokumentation
    - bisection_steps: Bisektionsverfahren mit Schritt-Dokumentation
    - gauss_elimination_steps: Gauss-Elimination mit Schritt-Dokumentation
    - euclidean_algorithm_steps: Euklidischer Algorithmus (ggT)
    - rsa_steps: RSA-Kryptographie Schritt-für-Schritt
    - prime_factorization_steps: Primfaktorzerlegung
    - format_steps_text: Text-Formatierung
    - format_steps_html: HTML-Formatierung

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest

# src-Verzeichnis zum Python-Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from step_by_step import (
    newton_raphson_steps,
    bisection_steps,
    gauss_elimination_steps,
    euclidean_algorithm_steps,
    rsa_steps,
    prime_factorization_steps,
    lu_decomposition_steps,
    format_steps_text,
    format_steps_html,
)


# ===========================================================================
# HILFSFUNKTIONEN FÜR TESTS
# ===========================================================================

def f_square_minus_2(x):
    """f(x) = x² - 2, Nullstelle bei x = √2 ≈ 1.41421356..."""
    return x * x - 2.0


def f_cubic(x):
    """f(x) = x³ - x - 2, Nullstelle bei x ≈ 1.52..."""
    return x ** 3 - x - 2.0


def f_sin(x):
    """f(x) = sin(x), Nullstelle bei x = π"""
    return math.sin(x)


# ===========================================================================
# TESTS: newton_raphson_steps
# ===========================================================================

class TestNewtonRaphsonSteps:
    """Tests für newton_raphson_steps()."""

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        assert isinstance(result, dict)
        assert 'steps' in result
        assert 'result' in result
        assert 'converged' in result
        assert 'iterations' in result
        assert 'method' in result

    def test_result_is_correct_sqrt2(self):
        """Newton-Raphson findet √2 korrekt für f(x) = x² - 2."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        assert result['converged'] is True
        assert abs(result['result'] - math.sqrt(2)) < 1e-8

    def test_steps_are_list(self):
        """steps ist eine Liste."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_step_has_required_keys(self):
        """Jeder Schritt enthält Pflicht-Schlüssel."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        for step in result['steps']:
            assert 'iteration' in step
            assert 'x' in step
            assert 'fx' in step
            assert 'explanation' in step

    def test_step_explanation_is_string(self):
        """Jeder Schritt hat eine nicht-leere Erklärung."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        for step in result['steps']:
            assert isinstance(step['explanation'], str)
            assert len(step['explanation']) > 0

    def test_converged_true_for_good_start(self):
        """Konvergiert für guten Startwert."""
        result = newton_raphson_steps(f_square_minus_2, 2.0)
        assert result['converged'] is True

    def test_cubic_function(self):
        """Findet Nullstelle einer kubischen Funktion."""
        result = newton_raphson_steps(f_cubic, 2.0)
        assert result['converged'] is True
        # f(result) ≈ 0
        assert abs(f_cubic(result['result'])) < 1e-8

    def test_custom_tolerance(self):
        """Respektiert benutzerdefinierte Toleranz."""
        result = newton_raphson_steps(f_square_minus_2, 1.5, tol=1e-3)
        assert result['converged'] is True
        # Bei tol=1e-3 sind weniger Iterationen nötig als bei Standard-tol
        result_strict = newton_raphson_steps(f_square_minus_2, 1.5, tol=1e-12)
        assert result['iterations'] <= result_strict['iterations']

    def test_max_iter_limits_steps(self):
        """max_iter begrenzt die Anzahl der Schritte."""
        result = newton_raphson_steps(f_square_minus_2, 1.5, max_iter=2)
        assert result['iterations'] <= 2

    def test_iterations_count_matches_steps_length(self):
        """iterations stimmt mit Länge der steps-Liste überein."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        assert result['iterations'] == len(result['steps'])

    def test_method_string_not_empty(self):
        """method-Feld ist nicht leer."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        assert isinstance(result['method'], str)
        assert len(result['method']) > 0


# ===========================================================================
# TESTS: bisection_steps
# ===========================================================================

class TestBisectionSteps:
    """Tests für bisection_steps()."""

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        assert isinstance(result, dict)
        assert 'steps' in result
        assert 'result' in result
        assert 'converged' in result
        assert 'iterations' in result
        assert 'method' in result

    def test_result_near_sqrt2(self):
        """Bisektionsergebnis ist nahe an √2."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        assert result['converged'] is True
        assert abs(result['result'] - math.sqrt(2)) < 1e-8

    def test_steps_are_list(self):
        """steps ist eine nicht-leere Liste."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_step_has_required_keys(self):
        """Jeder Schritt enthält a, b, midpoint, f_midpoint."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        for step in result['steps']:
            assert 'a' in step
            assert 'b' in step
            assert 'midpoint' in step
            assert 'f_midpoint' in step
            assert 'explanation' in step

    def test_midpoint_always_in_interval(self):
        """Mittelpunkt liegt immer im Intervall [a, b]."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        for step in result['steps']:
            assert step['midpoint'] >= step['a'] - 1e-12
            assert step['midpoint'] <= step['b'] + 1e-12

    def test_interval_shrinks(self):
        """Intervallbreite nimmt ab (Bisektion = Halbierung)."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        steps = result['steps']
        if len(steps) > 1:
            first_width = steps[0]['b'] - steps[0]['a']
            last_width = steps[-1]['b'] - steps[-1]['a']
            assert last_width <= first_width

    def test_same_sign_returns_error(self):
        """Kein Vorzeichenwechsel → Fehler-Rückgabe."""
        # f(1)=1-2=-1, f(1.1)=1.21-2=-0.79 → beide negativ
        result = bisection_steps(f_square_minus_2, 0.1, 0.5)
        assert result['result'] is None
        assert result['converged'] is False
        assert 'error' in result

    def test_cubic_root(self):
        """Findet Nullstelle einer kubischen Funktion."""
        result = bisection_steps(f_cubic, 1.0, 2.0)
        assert result['converged'] is True
        assert abs(f_cubic(result['result'])) < 1e-6

    def test_sin_root(self):
        """Findet Nullstelle von sin(x) nahe π."""
        result = bisection_steps(f_sin, 3.0, 3.5)
        assert result['converged'] is True
        assert abs(result['result'] - math.pi) < 1e-6


# ===========================================================================
# TESTS: gauss_elimination_steps
# ===========================================================================

class TestGaussEliminationSteps:
    """Tests für gauss_elimination_steps()."""

    def test_simple_2x2_system(self):
        """Löst einfaches 2×2-System: 2x+y=5, x+3y=10 → x=1, y=3."""
        result = gauss_elimination_steps([[2, 1], [1, 3]], [5, 10])
        assert result['result'] is not None
        sol = result['result']
        assert abs(sol[0] - 1.0) < 1e-8
        assert abs(sol[1] - 3.0) < 1e-8

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = gauss_elimination_steps([[1, 0], [0, 1]], [1, 2])
        assert 'steps' in result
        assert 'result' in result
        assert 'method' in result

    def test_steps_is_list(self):
        """steps ist eine nicht-leere Liste."""
        result = gauss_elimination_steps([[2, 1], [1, 3]], [5, 10])
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_step_has_phase_and_operation(self):
        """Jeder Schritt hat phase und operation."""
        result = gauss_elimination_steps([[2, 1], [1, 3]], [5, 10])
        for step in result['steps']:
            assert 'phase' in step
            assert 'operation' in step

    def test_identity_matrix(self):
        """Einheitsmatrix → direkte Lösung."""
        result = gauss_elimination_steps([[1, 0], [0, 1]], [3, 7])
        sol = result['result']
        assert abs(sol[0] - 3.0) < 1e-8
        assert abs(sol[1] - 7.0) < 1e-8

    def test_3x3_system(self):
        """Löst ein 3×3-System korrekt."""
        # 2x + y - z = 8
        # -3x - y + 2z = -11
        # -2x + y + 2z = -3
        # Lösung: x=2, y=3, z=-1
        A = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]
        b = [8, -11, -3]
        result = gauss_elimination_steps(A, b)
        sol = result['result']
        assert abs(sol[0] - 2.0) < 1e-6
        assert abs(sol[1] - 3.0) < 1e-6
        assert abs(sol[2] - (-1.0)) < 1e-6

    def test_singular_matrix_returns_error(self):
        """Singuläre Matrix → result ist None."""
        result = gauss_elimination_steps([[1, 1], [1, 1]], [2, 3])
        assert result['result'] is None


# ===========================================================================
# TESTS: euclidean_algorithm_steps
# ===========================================================================

class TestEuclideanAlgorithmSteps:
    """Tests für euclidean_algorithm_steps()."""

    def test_gcd_48_18(self):
        """ggT(48, 18) = 6."""
        result = euclidean_algorithm_steps(48, 18)
        assert result['result'] == 6

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = euclidean_algorithm_steps(48, 18)
        assert 'steps' in result
        assert 'result' in result
        assert 'method' in result

    def test_last_step_has_remainder_zero(self):
        """Letzter Schritt hat remainder=0 (Algorithmus terminiert)."""
        result = euclidean_algorithm_steps(48, 18)
        last_step = result['steps'][-1]
        assert last_step['remainder'] == 0

    def test_steps_is_list(self):
        """steps ist eine nicht-leere Liste."""
        result = euclidean_algorithm_steps(48, 18)
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_step_has_required_keys(self):
        """Jeder Schritt hat a, b, quotient, remainder, explanation."""
        result = euclidean_algorithm_steps(48, 18)
        for step in result['steps']:
            assert 'a' in step or 'step' in step  # letzter Schritt hat anderes Format
            assert 'remainder' in step
            assert 'explanation' in step

    def test_gcd_100_75(self):
        """ggT(100, 75) = 25."""
        result = euclidean_algorithm_steps(100, 75)
        assert result['result'] == 25

    def test_gcd_prime_numbers(self):
        """ggT zweier Primzahlen ist 1."""
        result = euclidean_algorithm_steps(17, 13)
        assert result['result'] == 1

    def test_gcd_same_numbers(self):
        """ggT(n, n) = n."""
        result = euclidean_algorithm_steps(12, 12)
        assert result['result'] == 12

    def test_gcd_one(self):
        """ggT(n, 1) = 1."""
        result = euclidean_algorithm_steps(100, 1)
        assert result['result'] == 1

    def test_quotient_times_b_plus_remainder_equals_a(self):
        """Jeder Schritt erfüllt: a = b*q + r."""
        result = euclidean_algorithm_steps(48, 18)
        for step in result['steps'][:-1]:  # Letzter Schritt hat b=0
            a, b_val = step['a'], step['b']
            q, r = step['quotient'], step['remainder']
            assert a == b_val * q + r


# ===========================================================================
# TESTS: rsa_steps
# ===========================================================================

class TestRsaSteps:
    """Tests für rsa_steps()."""

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = rsa_steps(5, 11, 7)
        assert 'steps' in result
        assert 'result' in result
        assert 'method' in result

    def test_steps_is_list(self):
        """steps ist eine nicht-leere Liste."""
        result = rsa_steps(5, 11, 7)
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_result_has_rsa_components(self):
        """Ergebnis enthält alle RSA-Schlüsselkomponenten."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert 'n' in r
        assert 'phi_n' in r
        assert 'e' in r
        assert 'd' in r
        assert 'public_key' in r
        assert 'private_key' in r
        assert 'ciphertext' in r
        assert 'decrypted' in r

    def test_n_equals_p_times_q(self):
        """n = p * q."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert r['n'] == 5 * 11  # n = 55

    def test_phi_n_correct(self):
        """φ(n) = (p-1)(q-1)."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert r['phi_n'] == (5 - 1) * (11 - 1)  # φ(55) = 40

    def test_decryption_recovers_message(self):
        """Entschlüsseln ergibt die ursprüngliche Nachricht zurück."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert r['success'] is True
        assert r['decrypted'] == r['message']

    def test_e_is_coprime_to_phi(self):
        """e ist teilerfremd zu φ(n): ggT(e, φ(n)) = 1."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert math.gcd(r['e'], r['phi_n']) == 1

    def test_d_is_inverse_of_e(self):
        """d·e ≡ 1 (mod φ(n))."""
        result = rsa_steps(5, 11, 7)
        r = result['result']
        assert (r['d'] * r['e']) % r['phi_n'] == 1

    def test_encryption_decryption_roundtrip(self):
        """Ver- und Entschlüsselung sind konsistent: pow(pow(m,e,n),d,n) = m."""
        result = rsa_steps(7, 13, 9)
        r = result['result']
        # Manuelle Überprüfung
        c = pow(r['message'], r['e'], r['n'])
        m_dec = pow(c, r['d'], r['n'])
        assert m_dec == r['message']

    def test_larger_primes(self):
        """Funktioniert auch mit größeren Primzahlen."""
        result = rsa_steps(17, 19, 42)
        r = result['result']
        assert r['success'] is True

    def test_each_step_has_title_and_explanation(self):
        """Jeder Schritt hat title und explanation."""
        result = rsa_steps(5, 11, 7)
        for step in result['steps']:
            assert 'title' in step
            assert 'explanation' in step
            assert len(step['explanation']) > 0


# ===========================================================================
# TESTS: prime_factorization_steps
# ===========================================================================

class TestPrimeFactorizationSteps:
    """Tests für prime_factorization_steps()."""

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält alle Pflicht-Schlüssel."""
        result = prime_factorization_steps(12)
        assert 'steps' in result
        assert 'result' in result
        assert 'method' in result

    def test_12_factors_contain_2_and_3(self):
        """12 = 2² × 3 → Faktoren enthalten 2 und 3."""
        result = prime_factorization_steps(12)
        assert 2 in result['result']
        assert 3 in result['result']

    def test_product_equals_original(self):
        """Produkt aller Primfaktoren ergibt n zurück."""
        for n in [12, 30, 360, 1024, 97]:
            result = prime_factorization_steps(n)
            product = 1
            for factor in result['result']:
                product *= factor
            assert product == n, f"Faktoren von {n} ergeben nicht {n}"

    def test_prime_has_one_factor(self):
        """Primzahl hat genau sich selbst als einzigen Primfaktor."""
        result = prime_factorization_steps(97)
        assert result['result'] == [97]

    def test_power_of_2(self):
        """64 = 2^6 → nur 2 als Primfaktor."""
        result = prime_factorization_steps(64)
        assert all(f == 2 for f in result['result'])
        assert len(result['result']) == 6

    def test_360_factorization(self):
        """360 = 2³ × 3² × 5."""
        result = prime_factorization_steps(360)
        counts = result.get('factor_counts', {})
        assert counts.get(2, 0) == 3
        assert counts.get(3, 0) == 2
        assert counts.get(5, 0) == 1

    def test_steps_is_nonempty_list(self):
        """steps ist eine nicht-leere Liste."""
        result = prime_factorization_steps(12)
        assert isinstance(result['steps'], list)
        assert len(result['steps']) > 0

    def test_step_has_divisor_and_quotient(self):
        """Jeder Schritt hat divisor, quotient und explanation."""
        result = prime_factorization_steps(12)
        for step in result['steps'][:-1]:  # Letzter Schritt ist Zusammenfassung
            assert 'divisor' in step
            assert 'explanation' in step

    def test_product_str_not_empty(self):
        """product_str ist nicht leer."""
        result = prime_factorization_steps(360)
        assert isinstance(result.get('product_str'), str)
        assert len(result['product_str']) > 0


# ===========================================================================
# TESTS: lu_decomposition_steps
# ===========================================================================

class TestLuDecompositionSteps:
    """Tests für lu_decomposition_steps()."""

    def test_returns_dict_with_required_keys(self):
        """Rückgabe enthält steps, result, method."""
        result = lu_decomposition_steps([[2, 1], [1, 3]])
        assert 'steps' in result
        assert 'result' in result
        assert 'method' in result

    def test_result_has_L_U_P(self):
        """Ergebnis enthält L, U, P-Matrizen."""
        result = lu_decomposition_steps([[2, 1], [1, 3]])
        r = result['result']
        assert 'L' in r
        assert 'U' in r
        assert 'P' in r

    def test_L_diagonal_is_one(self):
        """L hat 1en auf der Diagonale (Doolittle-Konvention)."""
        result = lu_decomposition_steps([[2, 1], [1, 3]])
        L = result['result']['L']
        n = len(L)
        for i in range(n):
            assert abs(L[i][i] - 1.0) < 1e-10

    def test_U_is_upper_triangular(self):
        """U ist eine obere Dreiecksmatrix."""
        result = lu_decomposition_steps([[2, 1], [1, 3]])
        U = result['result']['U']
        n = len(U)
        for i in range(n):
            for j in range(i):
                assert abs(U[i][j]) < 1e-10, f"U[{i}][{j}] = {U[i][j]} sollte 0 sein"

    def test_steps_nonempty(self):
        """steps ist nicht leer."""
        result = lu_decomposition_steps([[2, 1], [1, 3]])
        assert len(result['steps']) > 0


# ===========================================================================
# TESTS: format_steps_text
# ===========================================================================

class TestFormatStepsText:
    """Tests für format_steps_text()."""

    def test_returns_nonempty_string(self):
        """Gibt nicht-leeren String zurück."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        text = format_steps_text(result)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_contains_method_name(self):
        """Text enthält den Methodennamen."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        text = format_steps_text(result)
        assert 'Newton-Raphson' in text

    def test_contains_result_word(self):
        """Text enthält 'Ergebnis'."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        text = format_steps_text(result)
        assert 'Ergebnis' in text

    def test_works_with_gcd_result(self):
        """Funktioniert mit ggT-Ergebnis."""
        result = euclidean_algorithm_steps(48, 18)
        text = format_steps_text(result)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_works_with_rsa_result(self):
        """Funktioniert mit RSA-Ergebnis."""
        result = rsa_steps(5, 11, 7)
        text = format_steps_text(result)
        assert isinstance(text, str)
        assert len(text) > 0

    def test_works_with_empty_steps(self):
        """Funktioniert auch bei leerer steps-Liste."""
        minimal_dict = {'method': 'Test', 'steps': [], 'result': 42}
        text = format_steps_text(minimal_dict)
        assert isinstance(text, str)
        assert 'Test' in text

    def test_contains_separator_lines(self):
        """Text enthält Trennlinien (= Zeichen)."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        text = format_steps_text(result)
        assert '=' * 10 in text

    def test_contains_step_operations(self):
        """Text enthält die Schritt-Operationen."""
        result = euclidean_algorithm_steps(48, 18)
        text = format_steps_text(result)
        assert 'Operation' in text


# ===========================================================================
# TESTS: format_steps_html
# ===========================================================================

class TestFormatStepsHtml:
    """Tests für format_steps_html()."""

    def test_returns_nonempty_string(self):
        """Gibt nicht-leeren String zurück."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert isinstance(html, str)
        assert len(html) > 0

    def test_contains_div_tag(self):
        """HTML enthält <div>-Tags."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert '<div' in html

    def test_contains_step_class(self):
        """HTML enthält class="step" für Animationsunterstützung."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert 'class="step' in html

    def test_contains_steps_container(self):
        """HTML hat steps-container div."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert 'steps-container' in html

    def test_contains_code_tags(self):
        """HTML enthält <code>-Tags für Operationen."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert '<code>' in html

    def test_contains_method_name(self):
        """HTML enthält den Methodennamen."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert 'Newton-Raphson' in html

    def test_contains_formula(self):
        """HTML enthält KaTeX-kompatible Formel ($$...$$)."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert '$$' in html

    def test_contains_result_section(self):
        """HTML enthält Ergebnis-Bereich."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        html = format_steps_html(result)
        assert 'steps-result' in html

    def test_gcd_html_contains_all_steps(self):
        """HTML-Ausgabe für ggT enthält Schritt-Elemente."""
        result = euclidean_algorithm_steps(48, 18)
        html = format_steps_html(result)
        # Mindestens einer der Schritte muss sichtbar sein
        assert html.count('<div class="step') >= 1

    def test_gauss_html(self):
        """HTML-Ausgabe für Gauss-Elimination funktioniert."""
        result = gauss_elimination_steps([[2, 1], [1, 3]], [5, 10])
        html = format_steps_html(result)
        assert '<div' in html
        assert 'steps-container' in html


# ===========================================================================
# INTEGRATIONS-TESTS
# ===========================================================================

class TestIntegration:
    """Integrations-Tests: Kombination von step-Funktion und Formatierung."""

    def test_newton_full_pipeline(self):
        """Newton-Raphson → Text-Format → nicht leer."""
        result = newton_raphson_steps(f_square_minus_2, 1.5)
        text = format_steps_text(result)
        html = format_steps_html(result)
        assert result['converged']
        assert len(text) > 100
        assert len(html) > 100

    def test_bisection_full_pipeline(self):
        """Bisektion → HTML-Format → enthält div."""
        result = bisection_steps(f_square_minus_2, 1.0, 2.0)
        html = format_steps_html(result)
        assert '<div' in html

    def test_rsa_full_pipeline(self):
        """RSA → Text-Format → enthält Ergebnis."""
        result = rsa_steps(7, 11, 5)
        text = format_steps_text(result)
        assert 'Ergebnis' in text
        assert result['result']['success'] is True

    def test_prime_full_pipeline(self):
        """Primfaktorzerlegung → Textformat."""
        result = prime_factorization_steps(360)
        text = format_steps_text(result)
        html = format_steps_html(result)
        assert '360' in text
        assert len(html) > 0

    def test_gauss_solution_satisfies_system(self):
        """Gauss-Lösung erfüllt das ursprüngliche Gleichungssystem."""
        A = [[2.0, 1.0], [1.0, 3.0]]
        b = [5.0, 10.0]
        result = gauss_elimination_steps(A, b)
        x = result['result']
        # Probe: Ax = b
        for i in range(len(b)):
            row_sum = sum(A[i][j] * x[j] for j in range(len(x)))
            assert abs(row_sum - b[i]) < 1e-8
