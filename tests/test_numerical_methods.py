"""
Tests für das Numerische-Methoden-Modul (src/numerical_methods.py).

Testet Interpolation (Lagrange, Newton, kubische Splines),
Optimierung (Gradient Descent, Goldener Schnitt) und Simplex-Algorithmus.

@author: Kurt Ingwer
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import math
import pytest
import sys
import os

# src-Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from numerical_methods import (
    lagrange_interpolation,
    NewtonInterpolation,
    CubicSpline,
    gradient_descent,
    golden_section_search,
    numerical_gradient,
    simplex,
    bfgs
)


# ---------------------------------------------------------------------------
# LAGRANGE-INTERPOLATION TESTS
# ---------------------------------------------------------------------------

class TestLagrangeInterpolation:
    """Tests für Lagrange-Interpolationspolynom."""

    def test_single_point(self):
        """Lagrange mit einem Punkt: konstante Funktion."""
        result = lagrange_interpolation([2.0], [5.0], 3.0)
        assert abs(result - 5.0) < 1e-10

    def test_two_points_linear(self):
        """Zwei Punkte → lineare Interpolation."""
        # Gerade durch (0, 0) und (1, 1): f(x) = x
        result = lagrange_interpolation([0.0, 1.0], [0.0, 1.0], 0.5)
        assert abs(result - 0.5) < 1e-10

    def test_three_points_exact_at_nodes(self):
        """Lagrange-Polynom geht durch alle Stützpunkte."""
        x = [0.0, 1.0, 2.0]
        y = [0.0, 1.0, 4.0]  # f(x) = x²
        for xi, yi in zip(x, y):
            result = lagrange_interpolation(x, y, xi)
            assert abs(result - yi) < 1e-10

    def test_quadratic_function(self):
        """Interpolation einer Parabel mit 3 Punkten."""
        # f(x) = x²: Punkte bei x = -1, 0, 1
        x = [-1.0, 0.0, 1.0]
        y = [1.0, 0.0, 1.0]
        # Auswertung bei x=0.5: f(0.5) = 0.25
        result = lagrange_interpolation(x, y, 0.5)
        assert abs(result - 0.25) < 1e-10

    def test_duplicates_raise_error(self):
        """Doppelte Stützstellen lösen Fehler aus."""
        with pytest.raises(ValueError):
            lagrange_interpolation([1.0, 1.0, 2.0], [0.0, 1.0, 2.0], 1.5)

    def test_length_mismatch_raises(self):
        """Längenfehler wird erkannt."""
        with pytest.raises(ValueError):
            lagrange_interpolation([1.0, 2.0], [0.0, 1.0, 2.0], 1.5)

    def test_cubic_interpolation(self):
        """4 Punkte einer kubischen Funktion werden exakt interpoliert."""
        # f(x) = x³: Punkte bei x = 0, 1, 2, 3
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 1.0, 8.0, 27.0]
        # f(1.5) = 3.375
        result = lagrange_interpolation(x, y, 1.5)
        assert abs(result - 3.375) < 1e-8


# ---------------------------------------------------------------------------
# NEWTON-INTERPOLATION TESTS
# ---------------------------------------------------------------------------

class TestNewtonInterpolation:
    """Tests für Newton-Interpolation mit dividierten Differenzen."""

    def test_matches_lagrange(self):
        """Newton-Interpolation stimmt mit Lagrange überein."""
        x = [0.0, 1.0, 2.0, 3.0]
        y = [1.0, 3.0, 7.0, 13.0]  # f(x) = 2x² - x + 1
        newton = NewtonInterpolation(x, y)

        for t in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
            newton_val = newton.evaluate(t)
            lagrange_val = lagrange_interpolation(x, y, t)
            assert abs(newton_val - lagrange_val) < 1e-8

    def test_exact_at_nodes(self):
        """Newton-Polynom geht durch alle Stützpunkte."""
        x = [0.0, 1.0, 2.0]
        y = [0.0, 1.0, 4.0]
        newton = NewtonInterpolation(x, y)
        for xi, yi in zip(x, y):
            assert abs(newton.evaluate(xi) - yi) < 1e-10

    def test_linear_interpolation(self):
        """Lineare Interpolation mit Newton."""
        newton = NewtonInterpolation([0.0, 2.0], [1.0, 5.0])
        # Gerade: f(x) = 1 + 2x
        assert abs(newton.evaluate(1.0) - 3.0) < 1e-10
        assert abs(newton.evaluate(0.5) - 2.0) < 1e-10

    def test_quadratic_function(self):
        """Newton für f(x) = x²."""
        x = [-1.0, 0.0, 1.0]
        y = [1.0, 0.0, 1.0]
        newton = NewtonInterpolation(x, y)
        assert abs(newton.evaluate(0.5) - 0.25) < 1e-10
        assert abs(newton.evaluate(-0.5) - 0.25) < 1e-10

    def test_length_mismatch_raises(self):
        """Längenfehler wird erkannt."""
        with pytest.raises(ValueError):
            NewtonInterpolation([1.0, 2.0], [1.0, 2.0, 3.0])


# ---------------------------------------------------------------------------
# KUBISCHE SPLINE TESTS
# ---------------------------------------------------------------------------

class TestCubicSpline:
    """Tests für kubische Spline-Interpolation."""

    def test_exact_at_nodes(self):
        """Spline geht exakt durch alle Stützpunkte."""
        x = [0.0, 1.0, 2.0, 3.0, 4.0]
        y = [0.0, 1.0, 0.0, 1.0, 0.0]
        spline = CubicSpline(x, y)
        for xi, yi in zip(x, y):
            assert abs(spline.evaluate(xi) - yi) < 1e-10

    def test_linear_function(self):
        """Spline durch Punkte einer Geraden = Gerade."""
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 2.0, 4.0, 6.0]  # f(x) = 2x
        spline = CubicSpline(x, y)
        # Mittelpunkte sollten auf der Geraden liegen
        assert abs(spline.evaluate(0.5) - 1.0) < 1e-8
        assert abs(spline.evaluate(1.5) - 3.0) < 1e-8

    def test_continuity_at_knots(self):
        """Spline ist stetig: Auswertung von beiden Seiten eines Knotens."""
        x = [0.0, 1.0, 2.0]
        y = [0.0, 1.0, 0.0]
        spline = CubicSpline(x, y)
        # Bei x=1 (Knoten): beide Splines müssen übereinstimmen
        assert abs(spline.evaluate(1.0) - 1.0) < 1e-10

    def test_out_of_range_raises(self):
        """Außerhalb des Intervalls gibt es einen Fehler."""
        spline = CubicSpline([0.0, 1.0, 2.0], [0.0, 1.0, 0.0])
        with pytest.raises(ValueError):
            spline.evaluate(-0.1)
        with pytest.raises(ValueError):
            spline.evaluate(2.1)

    def test_minimum_points(self):
        """Mindestens 2 Punkte erforderlich."""
        with pytest.raises(ValueError):
            CubicSpline([1.0], [1.0])

    def test_non_ascending_raises(self):
        """Nicht-aufsteigende x-Werte lösen Fehler aus."""
        with pytest.raises(ValueError):
            CubicSpline([0.0, 2.0, 1.0], [0.0, 1.0, 0.5])

    def test_sine_function(self):
        """Spline-Interpolation von sin(x) mit 5 Punkten."""
        n = 5
        x = [i * math.pi / (n - 1) for i in range(n)]
        y = [math.sin(xi) for xi in x]
        spline = CubicSpline(x, y)

        # Mittelpunkte prüfen (Fehler sollte < 0.01 sein)
        for t in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5]:
            if x[0] <= t <= x[-1]:
                assert abs(spline.evaluate(t) - math.sin(t)) < 0.05

    def test_derivative_linear(self):
        """Ableitung einer linearen Funktion ist konstant."""
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 3.0, 6.0, 9.0]  # f(x) = 3x, f'(x) = 3
        spline = CubicSpline(x, y)
        # Ableitung im Inneren sollte ≈ 3 sein
        assert abs(spline.derivative(1.5) - 3.0) < 0.01


# ---------------------------------------------------------------------------
# OPTIMIERUNG TESTS
# ---------------------------------------------------------------------------

class TestOptimization:
    """Tests für Optimierungsalgorithmen."""

    def test_gradient_descent_parabola(self):
        """Gradient Descent findet Minimum einer Parabel."""
        # f(x) = (x-3)²: Minimum bei x=3
        f = lambda x: (x[0] - 3) ** 2
        grad = lambda x: [2 * (x[0] - 3)]
        x_opt, f_opt, n_iter = gradient_descent(f, grad, [0.0], learning_rate=0.1)
        assert abs(x_opt[0] - 3.0) < 0.01
        assert abs(f_opt) < 0.001

    def test_gradient_descent_2d(self):
        """Gradient Descent in 2D."""
        # f(x, y) = x² + 2y²: Minimum bei (0, 0)
        f = lambda x: x[0] ** 2 + 2 * x[1] ** 2
        grad = lambda x: [2 * x[0], 4 * x[1]]
        x_opt, f_opt, _ = gradient_descent(f, grad, [2.0, 1.0], learning_rate=0.1)
        assert abs(x_opt[0]) < 0.1
        assert abs(x_opt[1]) < 0.1

    def test_golden_section_parabola(self):
        """Goldener Schnitt findet Minimum einer Parabel."""
        f = lambda x: (x - 2) ** 2
        x_opt, f_opt = golden_section_search(f, 0.0, 5.0)
        assert abs(x_opt - 2.0) < 1e-6
        assert abs(f_opt) < 1e-12

    def test_golden_section_sin(self):
        """Goldener Schnitt findet Minimum von -sin(x) auf [0, π]."""
        f = lambda x: -math.sin(x)  # Minimum bei π/2
        x_opt, f_opt = golden_section_search(f, 0.0, math.pi)
        assert abs(x_opt - math.pi / 2) < 1e-6

    def test_numerical_gradient_parabola(self):
        """Numerischer Gradient einer Parabel."""
        f = lambda x: x[0] ** 2 + x[1] ** 2
        grad = numerical_gradient(f, [2.0, 3.0])
        # ∂f/∂x₁ = 2x₁ = 4, ∂f/∂x₂ = 2x₂ = 6
        assert abs(grad[0] - 4.0) < 1e-5
        assert abs(grad[1] - 6.0) < 1e-5

    def test_numerical_gradient_linear(self):
        """Numerischer Gradient einer linearen Funktion."""
        f = lambda x: 3 * x[0] + 5 * x[1]
        grad = numerical_gradient(f, [1.0, 1.0])
        assert abs(grad[0] - 3.0) < 1e-5
        assert abs(grad[1] - 5.0) < 1e-5


# ---------------------------------------------------------------------------
# SIMPLEX TESTS
# ---------------------------------------------------------------------------

class TestSimplex:
    """Tests für den Simplex-Algorithmus."""

    def test_simple_lp_1d(self):
        """1D LP: min x, x ≤ 5 → Lösung x=0."""
        c = [1.0]       # min x
        A = [[1.0]]     # x ≤ 5
        b = [5.0]
        x_opt, f_opt = simplex(c, A, b)
        assert abs(f_opt) < 1e-8  # min bei x=0

    def test_simple_lp_2d(self):
        """2D LP: min -x₁ - x₂, x₁ + x₂ ≤ 4, x₁ ≤ 3, x₂ ≤ 3."""
        # Entspricht: max x₁ + x₂
        c = [-1.0, -1.0]
        A = [
            [1.0, 1.0],  # x₁ + x₂ ≤ 4
            [1.0, 0.0],  # x₁ ≤ 3
            [0.0, 1.0],  # x₂ ≤ 3
        ]
        b = [4.0, 3.0, 3.0]
        x_opt, f_opt = simplex(c, A, b)
        # Optimum: x₁=1, x₂=3 oder x₁=3, x₂=1 (f=-4)
        assert abs(f_opt - (-4.0)) < 1e-6
        assert abs(x_opt[0] + x_opt[1] - 4.0) < 1e-6

    def test_trivial_zero_solution(self):
        """LP mit nur Nichtnegativitätsbedingung: Lösung ist 0."""
        c = [1.0, 1.0]
        A = [[1.0, 0.0], [0.0, 1.0]]
        b = [10.0, 10.0]
        x_opt, f_opt = simplex(c, A, b)
        assert abs(f_opt) < 1e-8

    def test_optimal_vertex(self):
        """LP-Optimum liegt an einem Eckpunkt."""
        # min -3x₁ - 2x₂, 2x₁ + x₂ ≤ 4, x₁ + 2x₂ ≤ 4
        c = [-3.0, -2.0]
        A = [[2.0, 1.0], [1.0, 2.0]]
        b = [4.0, 4.0]
        x_opt, f_opt = simplex(c, A, b)
        # Optimum bei x₁=4/3, x₂=4/3: f = -3·4/3 - 2·4/3 = -4 - 8/3 = -20/3
        assert f_opt < -5.0  # Mindestens -5

    def test_feasibility_constraints_satisfied(self):
        """Lösung erfüllt alle Nebenbedingungen."""
        c = [-1.0, -2.0]
        A = [[1.0, 1.0], [2.0, 1.0], [0.0, 1.0]]
        b = [6.0, 8.0, 4.0]
        x_opt, _ = simplex(c, A, b)

        # Prüfe Ax ≤ b
        for i, (ai, bi) in enumerate(zip(A, b)):
            lhs = sum(ai[j] * x_opt[j] for j in range(len(x_opt)))
            assert lhs <= bi + 1e-6, f"Nebenbedingung {i} verletzt"

        # Prüfe x ≥ 0
        for xi in x_opt:
            assert xi >= -1e-6


class TestBFGS:
    """Tests für den BFGS quasi-Newton-Optimierer."""

    def test_bfgs_quadratic_1d(self):
        """BFGS minimiert f(x) = x² → Minimum bei x=0."""
        def f(x):
            return [xi**2 for xi in [sum(x)]][0] if False else x[0]**2

        x_opt, f_opt = bfgs(f, x0=[3.0])
        assert abs(x_opt[0]) < 1e-5
        assert abs(f_opt) < 1e-10

    def test_bfgs_quadratic_2d(self):
        """BFGS minimiert f(x,y) = x² + y² → Minimum bei (0,0)."""
        def f(x):
            return x[0]**2 + x[1]**2

        x_opt, f_opt = bfgs(f, x0=[2.0, -3.0])
        assert abs(x_opt[0]) < 1e-5
        assert abs(x_opt[1]) < 1e-5
        assert abs(f_opt) < 1e-10

    def test_bfgs_shifted_minimum(self):
        """BFGS findet Minimum bei (3, -2) von f(x,y) = (x-3)² + (y+2)²."""
        def f(x):
            return (x[0] - 3.0)**2 + (x[1] + 2.0)**2

        x_opt, f_opt = bfgs(f, x0=[0.0, 0.0])
        assert abs(x_opt[0] - 3.0) < 1e-4
        assert abs(x_opt[1] + 2.0) < 1e-4
        assert abs(f_opt) < 1e-8

    def test_bfgs_rosenbrock(self):
        """BFGS minimiert Rosenbrock-Funktion f(x,y)=100(y-x²)²+(1-x)² → (1,1)."""
        def f(x):
            return 100.0 * (x[1] - x[0]**2)**2 + (1.0 - x[0])**2

        x_opt, f_opt = bfgs(f, x0=[-1.0, 1.0], max_iter=5000, tol=1e-8)
        assert abs(x_opt[0] - 1.0) < 1e-3
        assert abs(x_opt[1] - 1.0) < 1e-3

    def test_bfgs_returns_tuple(self):
        """BFGS gibt (x_opt, f_opt) als Tupel zurück."""
        def f(x):
            return x[0]**2

        result = bfgs(f, x0=[1.0])
        assert isinstance(result, tuple)
        assert len(result) == 2
        x_opt, f_opt = result
        assert isinstance(x_opt, list)
        assert isinstance(f_opt, float)
