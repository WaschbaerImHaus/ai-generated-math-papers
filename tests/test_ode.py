"""
Tests für das ODE-Modul (Gewöhnliche Differentialgleichungen).

Autor: Reisen macht Spass... mit Pia und Dirk e.Kfm.
Erstellt: 2026-03-05
Letzte Änderung: 2026-03-05
"""

import pytest
import math
from src.ode import (
    euler_method,
    runge_kutta4,
    runge_kutta45,
    solve_linear_ode,
    laplace_transform,
    inverse_laplace,
)


# ─────────────────────────────────────────────────────────
# Euler-Verfahren
# ─────────────────────────────────────────────────────────

class TestEulerMethod:
    """Tests für das explizite Euler-Verfahren."""

    def test_euler_exponential_decay(self):
        # dy/dt = -y, y(0) = 1 -> Lösung: y(t) = e^(-t)
        f = lambda t, y: -y
        t, y = euler_method(f, t0=0.0, y0=1.0, t_end=1.0, n_steps=1000)
        assert y[-1] == pytest.approx(math.exp(-1), rel=0.01)  # 1% Fehler OK für Euler

    def test_euler_linear_growth(self):
        # dy/dt = 1, y(0) = 0 -> Lösung: y(t) = t
        f = lambda t, y: 1.0
        t, y = euler_method(f, t0=0.0, y0=0.0, t_end=2.0, n_steps=100)
        assert y[-1] == pytest.approx(2.0, rel=1e-4)

    def test_euler_returns_correct_length(self):
        f = lambda t, y: y
        t, y = euler_method(f, t0=0.0, y0=1.0, t_end=1.0, n_steps=50)
        assert len(t) == 51  # n_steps + 1 Punkte
        assert len(y) == 51

    def test_euler_initial_condition(self):
        f = lambda t, y: t
        t, y = euler_method(f, t0=0.0, y0=5.0, t_end=1.0, n_steps=10)
        assert y[0] == pytest.approx(5.0)
        assert t[0] == pytest.approx(0.0)

    def test_euler_constant_function(self):
        # dy/dt = 0 -> y = const
        f = lambda t, y: 0.0
        t, y = euler_method(f, t0=0.0, y0=3.14, t_end=5.0, n_steps=100)
        assert all(yi == pytest.approx(3.14) for yi in y)


# ─────────────────────────────────────────────────────────
# Runge-Kutta 4. Ordnung
# ─────────────────────────────────────────────────────────

class TestRungeKutta4:
    """Tests für das klassische Runge-Kutta-Verfahren 4. Ordnung."""

    def test_rk4_exponential_decay(self):
        # dy/dt = -y, y(0) = 1 -> Lösung: y(t) = e^(-t)
        f = lambda t, y: -y
        t, y = runge_kutta4(f, t0=0.0, y0=1.0, t_end=2.0, n_steps=100)
        assert y[-1] == pytest.approx(math.exp(-2), rel=1e-6)

    def test_rk4_exponential_growth(self):
        # dy/dt = y, y(0) = 1 -> Lösung: y(t) = e^t
        f = lambda t, y: y
        t, y = runge_kutta4(f, t0=0.0, y0=1.0, t_end=1.0, n_steps=100)
        assert y[-1] == pytest.approx(math.e, rel=1e-6)

    def test_rk4_harmonic_oscillator(self):
        # d^2x/dt^2 = -x als System: dy/dt = [y[1], -y[0]]
        # Anfangsbedingungen: x(0)=1, x'(0)=0 -> Lösung: x(t)=cos(t)
        f = lambda t, y: [y[1], -y[0]]
        t, y = runge_kutta4(f, t0=0.0, y0=[1.0, 0.0], t_end=2*math.pi, n_steps=1000)
        # x(2pi) = cos(2pi) = 1
        assert y[-1][0] == pytest.approx(1.0, rel=1e-5)

    def test_rk4_logistic_equation(self):
        # dy/dt = r*y*(1 - y/K), y(0) = y0
        # Lösung: y(t) = K / (1 + (K/y0 - 1)*e^(-r*t))
        r, K, y0 = 1.0, 10.0, 1.0
        f = lambda t, y: r * y * (1 - y / K)
        t, y = runge_kutta4(f, t0=0.0, y0=y0, t_end=5.0, n_steps=500)
        t_end = 5.0
        expected = K / (1 + (K/y0 - 1) * math.exp(-r * t_end))
        assert y[-1] == pytest.approx(expected, rel=1e-4)

    def test_rk4_vs_euler_accuracy(self):
        # RK4 muss genauer sein als Euler für gleiche Schrittzahl
        f = lambda t, y: -y
        n = 50
        _, y_euler = euler_method(f, 0.0, 1.0, 1.0, n)
        _, y_rk4   = runge_kutta4(f, 0.0, 1.0, 1.0, n)
        exact = math.exp(-1)
        assert abs(y_rk4[-1] - exact) < abs(y_euler[-1] - exact)

    def test_rk4_initial_condition_preserved(self):
        f = lambda t, y: y
        t, y = runge_kutta4(f, t0=0.0, y0=7.5, t_end=1.0, n_steps=10)
        assert y[0] == pytest.approx(7.5)


# ─────────────────────────────────────────────────────────
# Runge-Kutta-Fehlberg (adaptiv)
# ─────────────────────────────────────────────────────────

class TestRungeKutta45:
    """Tests für das adaptive Runge-Kutta-Fehlberg-Verfahren (RK45)."""

    def test_rk45_exponential(self):
        # dy/dt = -y, y(0) = 1 -> y(1) = e^(-1)
        f = lambda t, y: -y
        t, y = runge_kutta45(f, t0=0.0, y0=1.0, t_end=1.0, tol=1e-8)
        assert y[-1] == pytest.approx(math.exp(-1), rel=1e-6)

    def test_rk45_reaches_endpoint(self):
        f = lambda t, y: y
        t, y = runge_kutta45(f, t0=0.0, y0=1.0, t_end=2.0, tol=1e-6)
        assert t[-1] == pytest.approx(2.0, abs=1e-10)

    def test_rk45_monotone_time(self):
        # Zeitschritte müssen monoton steigen
        f = lambda t, y: -2 * t * y
        t, y = runge_kutta45(f, t0=0.0, y0=1.0, t_end=2.0, tol=1e-6)
        assert all(t[i] < t[i+1] for i in range(len(t)-1))


# ─────────────────────────────────────────────────────────
# Lineare ODE mit konstanten Koeffizienten
# ─────────────────────────────────────────────────────────

class TestLinearODE:
    """Tests für lineare ODE mit konstanten Koeffizienten."""

    def test_first_order_homogeneous(self):
        # y' + ay = 0 -> y(t) = C*e^(-at)
        # y'(t) + 2*y(t) = 0, y(0) = 3 -> y(t) = 3*e^(-2t)
        coeffs = [1.0, 2.0]  # y' + 2y = 0
        initial = [3.0]      # y(0) = 3
        t_vals = [0.0, 0.5, 1.0]
        y_vals = solve_linear_ode(coeffs, initial, t_vals, forcing=None)
        assert y_vals[0] == pytest.approx(3.0)
        assert y_vals[1] == pytest.approx(3 * math.exp(-1.0), rel=1e-4)
        assert y_vals[2] == pytest.approx(3 * math.exp(-2.0), rel=1e-4)

    def test_second_order_oscillator(self):
        # y'' + y = 0, y(0)=1, y'(0)=0 -> y(t) = cos(t)
        coeffs = [1.0, 0.0, 1.0]   # y'' + 0*y' + y = 0
        initial = [1.0, 0.0]       # y(0)=1, y'(0)=0
        t_vals = [0.0, math.pi/2, math.pi]
        y_vals = solve_linear_ode(coeffs, initial, t_vals, forcing=None)
        assert y_vals[0] == pytest.approx(1.0, abs=1e-6)
        assert y_vals[1] == pytest.approx(0.0, abs=1e-4)  # cos(pi/2) = 0
        assert y_vals[2] == pytest.approx(-1.0, abs=1e-4) # cos(pi) = -1

    def test_first_order_with_forcing(self):
        # y' + y = 1, y(0) = 0 -> y(t) = 1 - e^(-t)
        coeffs = [1.0, 1.0]    # y' + y
        initial = [0.0]
        t_vals = [0.0, 1.0, 2.0]
        forcing = lambda t: 1.0
        y_vals = solve_linear_ode(coeffs, initial, t_vals, forcing=forcing)
        assert y_vals[0] == pytest.approx(0.0, abs=1e-6)
        assert y_vals[1] == pytest.approx(1 - math.exp(-1), rel=1e-4)
        assert y_vals[2] == pytest.approx(1 - math.exp(-2), rel=1e-4)


# ─────────────────────────────────────────────────────────
# Laplace-Transformation
# ─────────────────────────────────────────────────────────

class TestLaplaceTransform:
    """Tests für die (numerische) Laplace-Transformation."""

    def test_laplace_constant(self):
        # L{1}(s) = 1/s
        f = lambda t: 1.0
        s = 2.0
        result = laplace_transform(f, s)
        assert result == pytest.approx(1.0 / s, rel=1e-3)

    def test_laplace_exponential(self):
        # L{e^(at)}(s) = 1/(s-a) für s > a
        a = 1.0
        f = lambda t: math.exp(a * t)
        s = 3.0
        result = laplace_transform(f, s)
        assert result == pytest.approx(1.0 / (s - a), rel=1e-3)

    def test_laplace_sine(self):
        # L{sin(omega*t)}(s) = omega / (s^2 + omega^2)
        omega = 2.0
        f = lambda t: math.sin(omega * t)
        s = 3.0
        result = laplace_transform(f, s)
        expected = omega / (s**2 + omega**2)
        assert result == pytest.approx(expected, rel=1e-3)

    def test_laplace_cosine(self):
        # L{cos(omega*t)}(s) = s / (s^2 + omega^2)
        omega = 1.0
        f = lambda t: math.cos(omega * t)
        s = 2.0
        result = laplace_transform(f, s)
        expected = s / (s**2 + omega**2)
        assert result == pytest.approx(expected, rel=1e-3)

    def test_inverse_laplace_exponential(self):
        # L^{-1}{1/(s+1)}(t) = e^(-t)
        F = lambda s: 1.0 / (s + 1)
        t = 1.0
        result = inverse_laplace(F, t)
        assert result == pytest.approx(math.exp(-1), rel=1e-2)

    def test_inverse_laplace_constant(self):
        # L^{-1}{1/s}(t) = 1 (Heaviside-Funktion)
        F = lambda s: 1.0 / s
        t = 0.5
        result = inverse_laplace(F, t)
        assert result == pytest.approx(1.0, rel=1e-2)
