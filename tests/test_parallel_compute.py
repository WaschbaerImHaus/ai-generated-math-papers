"""
@file test_parallel_compute.py
@brief Tests für das parallel_compute-Modul.
@description
    Überprüft die parallele symbolische und numerische Berechnung.
    Testet Korrektheit der Ergebnisse, Fehlerbehandlung und Edge-Cases.

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import pytest
import sympy as sp

# Pfad zum src/-Verzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from parallel_compute import (
    parallel_symbolic,
    parallel_integrate,
    parallel_limit,
    parallel_derivative,
    _integrate_single,
    _limit_single,
    _derivative_single,
)


# ===========================================================================
# TESTS: parallel_integrate
# ===========================================================================

class TestParallelIntegrate:
    """Tests für die parallele symbolische Integration."""

    def test_simple_polynomial_integration(self):
        """Überprüft die Integration einfacher Polynome: ∫x² dx = x³/3."""
        results = parallel_integrate(['x**2'], x='x', max_workers=1)
        assert 'x**2' in results
        # Ergebnis als SymPy-Ausdruck parsen und prüfen
        x = sp.Symbol('x')
        expected = sp.integrate(x**2, x)
        actual = sp.sympify(results['x**2'])
        # Differenz muss 0 sein (symbolisch)
        diff = sp.simplify(actual - expected)
        assert diff == 0, f"Integral falsch: {results['x**2']} ≠ {expected}"

    def test_trigonometric_integration(self):
        """Überprüft: ∫sin(x) dx = -cos(x)."""
        results = parallel_integrate(['sin(x)'], x='x', max_workers=1)
        x = sp.Symbol('x')
        expected = sp.integrate(sp.sin(x), x)
        actual = sp.sympify(results['sin(x)'])
        diff = sp.simplify(actual - expected)
        assert diff == 0, f"Integral falsch: {actual} ≠ {expected}"

    def test_exponential_integration(self):
        """Überprüft: ∫exp(x) dx = exp(x)."""
        results = parallel_integrate(['exp(x)'], x='x', max_workers=1)
        x = sp.Symbol('x')
        expected = sp.integrate(sp.exp(x), x)
        actual = sp.sympify(results['exp(x)'])
        diff = sp.simplify(actual - expected)
        assert diff == 0

    def test_multiple_expressions_parallel(self):
        """Überprüft die parallele Integration mehrerer Ausdrücke gleichzeitig."""
        exprs = ['x**2', 'sin(x)', 'exp(x)', 'x**3 + 2*x']
        results = parallel_integrate(exprs, x='x', max_workers=2)

        # Alle Ausdrücke müssen im Ergebnis sein
        assert set(results.keys()) == set(exprs)

        x = sp.Symbol('x')
        for expr_str in exprs:
            expected = sp.integrate(sp.sympify(expr_str), x)
            actual = sp.sympify(results[expr_str])
            diff = sp.simplify(actual - expected)
            assert diff == 0, f"Fehlerhaftes Integral für {expr_str}"

    def test_empty_list(self):
        """Leere Liste ergibt leeres Dictionary."""
        results = parallel_integrate([], x='x', max_workers=1)
        assert results == {}


# ===========================================================================
# TESTS: parallel_limit
# ===========================================================================

class TestParallelLimit:
    """Tests für die parallele Grenzwertberechnung."""

    def test_sinc_limit(self):
        """Überprüft: lim_{x→0} sin(x)/x = 1."""
        results = parallel_limit(['sin(x)/x'], x='x', point='0', max_workers=1)
        assert 'sin(x)/x' in results
        actual = sp.sympify(results['sin(x)/x'])
        assert sp.simplify(actual - 1) == 0, f"Erwartet 1, erhalten: {actual}"

    def test_polynomial_limit(self):
        """Überprüft: lim_{x→2} x² = 4."""
        results = parallel_limit(['x**2'], x='x', point='2', max_workers=1)
        actual = sp.sympify(results['x**2'])
        assert sp.simplify(actual - 4) == 0

    def test_limit_at_infinity(self):
        """Überprüft: lim_{x→∞} 1/x = 0."""
        results = parallel_limit(['1/x'], x='x', point='oo', max_workers=1)
        actual = sp.sympify(results['1/x'])
        assert sp.simplify(actual) == 0

    def test_multiple_limits_parallel(self):
        """Überprüft parallele Berechnung mehrerer Grenzwerte."""
        exprs = ['sin(x)/x', 'x**2', '(1-cos(x))/x**2']
        results = parallel_limit(exprs, x='x', point='0', max_workers=2)

        assert set(results.keys()) == set(exprs)

        # Erwartete Ergebnisse
        expected = {
            'sin(x)/x': 1,
            'x**2': 0,
            '(1-cos(x))/x**2': sp.Rational(1, 2),
        }
        x = sp.Symbol('x')
        for expr_str, exp_val in expected.items():
            actual = sp.sympify(results[expr_str])
            diff = sp.simplify(actual - exp_val)
            assert diff == 0, f"Grenzwert falsch für {expr_str}: {actual} ≠ {exp_val}"

    def test_empty_limit_list(self):
        """Leere Liste ergibt leeres Dictionary."""
        results = parallel_limit([], max_workers=1)
        assert results == {}


# ===========================================================================
# TESTS: parallel_derivative
# ===========================================================================

class TestParallelDerivative:
    """Tests für die parallele Ableitungsberechnung."""

    def test_polynomial_derivative(self):
        """Überprüft: d/dx x³ = 3x²."""
        results = parallel_derivative(['x**3'], x='x', order=1, max_workers=1)
        x = sp.Symbol('x')
        expected = 3 * x**2
        actual = sp.sympify(results['x**3'])
        diff = sp.simplify(actual - expected)
        assert diff == 0, f"Ableitung falsch: {actual} ≠ {expected}"

    def test_second_derivative(self):
        """Überprüft: d²/dx² x⁴ = 12x²."""
        results = parallel_derivative(['x**4'], x='x', order=2, max_workers=1)
        x = sp.Symbol('x')
        expected = 12 * x**2
        actual = sp.sympify(results['x**4'])
        diff = sp.simplify(actual - expected)
        assert diff == 0

    def test_trig_derivative(self):
        """Überprüft: d/dx sin(x) = cos(x)."""
        results = parallel_derivative(['sin(x)'], x='x', order=1, max_workers=1)
        x = sp.Symbol('x')
        expected = sp.cos(x)
        actual = sp.sympify(results['sin(x)'])
        diff = sp.simplify(actual - expected)
        assert diff == 0

    def test_multiple_derivatives_parallel(self):
        """Überprüft parallele Ableitungsberechnung mehrerer Ausdrücke."""
        exprs = ['x**3', 'sin(x)', 'exp(2*x)']
        results = parallel_derivative(exprs, x='x', order=2, max_workers=2)

        assert set(results.keys()) == set(exprs)

        x = sp.Symbol('x')
        expected_map = {
            'x**3': sp.diff(x**3, x, 2),
            'sin(x)': sp.diff(sp.sin(x), x, 2),
            'exp(2*x)': sp.diff(sp.exp(2*x), x, 2),
        }
        for expr_str, exp_val in expected_map.items():
            actual = sp.sympify(results[expr_str])
            diff = sp.simplify(actual - exp_val)
            assert diff == 0, f"2. Ableitung falsch für {expr_str}"

    def test_empty_derivative_list(self):
        """Leere Liste ergibt leeres Dictionary."""
        results = parallel_derivative([], max_workers=1)
        assert results == {}


# ===========================================================================
# TESTS: parallel_symbolic (allgemeine Parallelisierung)
# ===========================================================================

class TestParallelSymbolic:
    """Tests für die allgemeine parallele Ausführung."""

    def test_custom_function_parallel(self):
        """Überprüft parallele Ausführung einer benutzerdefinierten Funktion."""
        # Verwende _integrate_single als Beispiel-Funktion
        args_list = [('x**2', 'x'), ('sin(x)', 'x')]
        results = parallel_symbolic(_integrate_single, args_list, max_workers=2)

        # Alle Argumente müssen im Ergebnis sein
        assert len(results) == 2

        # Ergebnisse prüfen
        result_dict = {args[0]: val[1] for args, val in results.items()
                      if not isinstance(val, Exception)}

        x = sp.Symbol('x')
        assert sp.simplify(sp.sympify(result_dict['x**2']) - x**3/3) == 0
        assert sp.simplify(sp.sympify(result_dict['sin(x)']) + sp.cos(x)) == 0

    def test_single_worker(self):
        """max_workers=1 entspricht sequenzieller Ausführung."""
        args_list = [('x', 'x'), ('x**2', 'x')]
        results = parallel_symbolic(_integrate_single, args_list, max_workers=1)
        assert len(results) == 2

    def test_more_workers_than_tasks(self):
        """Mehr Worker als Aufgaben → kein Fehler."""
        args_list = [('x**2', 'x')]
        results = parallel_symbolic(_integrate_single, args_list, max_workers=8)
        assert len(results) == 1
