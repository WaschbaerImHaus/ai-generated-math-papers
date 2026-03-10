"""
@file test_analysis.py
@brief Tests für das Analysis-Modul (TDD).
@description
    Testet numerische und symbolische Analysis-Funktionen:
    - Grenzwerte (numerisch)
    - Ableitungen (numerische Differentiation)
    - Numerische Integration
    - Taylor-Entwicklung
    - Nullstellensuche
@author Reisen macht Spass... mit Pia und Dirk e.Kfm.
@date 2026-03-05
"""

import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from analysis import (
    numerical_derivative, numerical_integral,
    newton_raphson, bisection, brent_method,
    taylor_series, is_continuous
)


class TestNumericalDerivative:
    """Tests für numerische Ableitung."""

    def test_derivative_linear(self):
        """Ableitung von f(x) = 3x ist 3 (überall)."""
        f = lambda x: 3 * x
        result = numerical_derivative(f, 5.0)
        assert abs(result - 3.0) < 1e-6

    def test_derivative_quadratic(self):
        """Ableitung von f(x) = x^2 ist 2x. Bei x=3: 6."""
        f = lambda x: x**2
        result = numerical_derivative(f, 3.0)
        assert abs(result - 6.0) < 1e-6

    def test_derivative_sin(self):
        """Ableitung von sin(x) ist cos(x). Bei x=0: cos(0) = 1."""
        result = numerical_derivative(math.sin, 0.0)
        assert abs(result - 1.0) < 1e-6

    def test_derivative_exp(self):
        """Ableitung von e^x ist e^x. Bei x=1: e."""
        result = numerical_derivative(math.exp, 1.0)
        assert abs(result - math.e) < 1e-6

    def test_second_derivative(self):
        """Zweite Ableitung von x^3 ist 6x. Bei x=2: 12."""
        f = lambda x: x**3
        result = numerical_derivative(f, 2.0, order=2)
        assert abs(result - 12.0) < 1e-4


class TestNumericalIntegral:
    """Tests für numerische Integration (Gauss-Quadratur / Simpson)."""

    def test_constant_function(self):
        """Integral von 5 von 0 bis 3 = 15."""
        f = lambda x: 5
        result = numerical_integral(f, 0, 3)
        assert abs(result - 15.0) < 1e-6

    def test_linear_function(self):
        """Integral von x von 0 bis 4 = 8 (Dreiecksfläche)."""
        f = lambda x: x
        result = numerical_integral(f, 0, 4)
        assert abs(result - 8.0) < 1e-6

    def test_quadratic_function(self):
        """Integral von x^2 von 0 bis 3 = [x^3/3]_0^3 = 9."""
        f = lambda x: x**2
        result = numerical_integral(f, 0, 3)
        assert abs(result - 9.0) < 1e-6

    def test_sine_full_period(self):
        """Integral von sin(x) von 0 bis 2pi = 0."""
        result = numerical_integral(math.sin, 0, 2 * math.pi)
        assert abs(result) < 1e-6

    def test_sine_half_period(self):
        """Integral von sin(x) von 0 bis pi = 2."""
        result = numerical_integral(math.sin, 0, math.pi)
        assert abs(result - 2.0) < 1e-6

    def test_exponential(self):
        """Integral von e^x von 0 bis 1 = e - 1."""
        result = numerical_integral(math.exp, 0, 1)
        expected = math.e - 1
        assert abs(result - expected) < 1e-6


class TestRootFinding:
    """Tests für Nullstellensuche."""

    def test_newton_quadratic(self):
        """Newton-Raphson für x^2 - 4 = 0 -> x = 2."""
        f = lambda x: x**2 - 4
        root = newton_raphson(f, x0=3.0)
        assert abs(root - 2.0) < 1e-8

    def test_newton_sin(self):
        """Newton-Raphson für sin(x) nahe pi -> Nullstelle bei pi."""
        root = newton_raphson(math.sin, x0=3.0)
        assert abs(root - math.pi) < 1e-8

    def test_bisection_simple(self):
        """Bisektionsverfahren für x^2 - 9 = 0 in [0, 5] -> x = 3."""
        f = lambda x: x**2 - 9
        root = bisection(f, 0, 5)
        assert abs(root - 3.0) < 1e-8

    def test_bisection_requires_sign_change(self):
        """Bisektionsverfahren erfordert Vorzeichenwechsel."""
        f = lambda x: x**2 + 1  # Keine reelle Nullstelle
        with pytest.raises(ValueError):
            bisection(f, -2, 2)


class TestTaylorSeries:
    """Tests für Taylor-Reihenentwicklung."""

    def test_exp_taylor(self):
        """Taylor-Reihe von e^x an x=0 bis Grad 10, ausgewertet bei x=1."""
        # e^x = 1 + x + x^2/2! + x^3/3! + ...
        result = taylor_series(math.exp, center=0, degree=15, evaluate_at=1.0)
        assert abs(result - math.e) < 1e-6

    def test_sin_taylor(self):
        """Taylor-Reihe von sin(x) an x=0, ausgewertet bei x=pi/4."""
        result = taylor_series(math.sin, center=0, degree=15, evaluate_at=math.pi/4)
        expected = math.sin(math.pi/4)
        assert abs(result - expected) < 1e-8

    def test_cos_taylor(self):
        """Taylor-Reihe von cos(x) an x=0, ausgewertet bei x=pi/3."""
        result = taylor_series(math.cos, center=0, degree=15, evaluate_at=math.pi/3)
        expected = math.cos(math.pi/3)
        assert abs(result - expected) < 1e-8


class TestBrentMethod:
    """
    Tests für die Brent-Methode (Hybridverfahren aus Bisection + Sekante + IQI).
    Die Brent-Methode kombiniert die Konvergenzsicherheit der Bisection
    mit der Geschwindigkeit des Sekantenverfahrens und der IQI.
    """

    def test_brent_simple_root(self):
        """x^2 - 2 hat die Wurzel sqrt(2) ≈ 1.41421356 im Intervall [1, 2]."""
        f = lambda x: x**2 - 2
        root = brent_method(f, 1.0, 2.0)
        assert abs(root - math.sqrt(2)) < 1e-10, (
            f"Erwartete Wurzel sqrt(2) ≈ {math.sqrt(2):.10f}, erhielt {root:.10f}"
        )

    def test_brent_known_root(self):
        """sin(x) hat eine Nullstelle bei x = π ≈ 3.14159 im Intervall [3, 4]."""
        root = brent_method(math.sin, 3.0, 4.0)
        assert abs(root - math.pi) < 1e-10, (
            f"Erwartete Nullstelle bei π ≈ {math.pi:.10f}, erhielt {root:.10f}"
        )

    def test_brent_negative_root(self):
        """x^2 - 4 hat negative Wurzel x = -2 im Intervall [-3, -1]."""
        f = lambda x: x**2 - 4
        root = brent_method(f, -3.0, -1.0)
        assert abs(root - (-2.0)) < 1e-10, (
            f"Erwartete negative Wurzel -2.0, erhielt {root:.10f}"
        )

    def test_brent_close_bracket(self):
        """Sehr enges Startintervall [1.4, 1.5] für Wurzel von x^2 - 2."""
        f = lambda x: x**2 - 2
        root = brent_method(f, 1.4, 1.5, tol=1e-12)
        assert abs(root - math.sqrt(2)) < 1e-10, (
            f"Enges Intervall: Erwartete sqrt(2) ≈ {math.sqrt(2):.12f}, erhielt {root:.12f}"
        )

    def test_brent_no_sign_change_raises(self):
        """Wirft ValueError wenn f(a) und f(b) dasselbe Vorzeichen haben."""
        f = lambda x: x**2 + 1  # Immer positiv – keine reelle Nullstelle
        with pytest.raises(ValueError, match="Kein Vorzeichenwechsel"):
            brent_method(f, 1.0, 3.0)

    def test_brent_tolerance(self):
        """Brent-Methode hält die angegebene Toleranz ein."""
        f = lambda x: x**3 - x - 2  # Nullstelle bei ≈ 1.5214
        tol = 1e-9
        root = brent_method(f, 1.0, 2.0, tol=tol)
        # Prüfen: f(root) ist nahe 0 und Nullstelle ist korrekt
        assert abs(f(root)) < 1e-6, (
            f"f(root) = {f(root):.2e} sollte kleiner als 1e-6 sein"
        )

    def test_brent_vs_bisection_same_root(self):
        """Brent-Methode und Bisection finden dieselbe Nullstelle (cos(x) in [1, 2])."""
        f = lambda x: math.cos(x)  # Nullstelle bei π/2 ≈ 1.5708 in [1, 2]
        root_brent = brent_method(f, 1.0, 2.0, tol=1e-12)
        root_bisection = bisection(f, 1.0, 2.0, tol=1e-12)
        # Beide Ergebnisse sollten sehr nah beieinander liegen
        assert abs(root_brent - root_bisection) < 1e-9, (
            f"Brent: {root_brent:.12f}, Bisection: {root_bisection:.12f} "
            f"– Differenz: {abs(root_brent - root_bisection):.2e}"
        )
        # Und beide nahe dem echten Wert π/2
        assert abs(root_brent - math.pi / 2) < 1e-10

    def test_brent_cubic(self):
        """Kubisches Polynom x^3 - 6x^2 + 11x - 6 hat Wurzeln bei x=1, 2, 3."""
        # Suche Nullstelle bei x = 3 im Intervall [2.5, 3.5]
        f = lambda x: x**3 - 6 * x**2 + 11 * x - 6
        root = brent_method(f, 2.5, 3.5, tol=1e-12)
        assert abs(root - 3.0) < 1e-10, (
            f"Erwartete Wurzel 3.0, erhielt {root:.12f}"
        )
        # Suche auch Wurzel bei x = 1 im Intervall [0.5, 1.5]
        root2 = brent_method(f, 0.5, 1.5, tol=1e-12)
        assert abs(root2 - 1.0) < 1e-10, (
            f"Erwartete Wurzel 1.0, erhielt {root2:.12f}"
        )


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
