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
    newton_raphson, bisection,
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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
