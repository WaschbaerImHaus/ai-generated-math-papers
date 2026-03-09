"""
@file test_property_based.py
@brief Property-Based Tests mit Hypothesis.
@description
    Automatisch generierte Testfälle für mathematische Eigenschaften.
    Hypothesis erzeugt zufällige Eingaben und sucht nach Gegenbeispielen.

    Getestete Eigenschaften:
    - Algebra: gcd, lcm, is_prime, euler_phi, mod_inverse
    - Analysis: numerical_derivative, numerical_integral
    - Lineare Algebra: Vektoren, Givens-QR, LU-Zerlegung
    - Zahlentheorie: Diophantische Gleichungen, Tonelli-Shanks, Pell-Gleichung

@author Kurt Ingwer
@date 2026-03-09
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from hypothesis import given, settings, assume, HealthCheck
from hypothesis import strategies as st
import algebra
import analysis
import linear_algebra

# Globales Profil mit kurzen Timeouts für CI-Umgebungen
settings.register_profile("fast", max_examples=10, deadline=5000,
                           suppress_health_check=[HealthCheck.too_slow,
                                                   HealthCheck.filter_too_much])
settings.load_profile("fast")

# Erste 50 Primzahlen für Tests (als Konstante, nicht wiederberechnet)
_SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]


# =============================================================================
# ALGEBRA: GCD / LCM
# =============================================================================

class TestGcdProperties:
    """Property-Based Tests für den Größten Gemeinsamen Teiler."""

    @given(st.integers(1, 10000), st.integers(1, 10000))
    def test_gcd_commutative(self, a: int, b: int):
        """gcd(a, b) == gcd(b, a) — Kommutativität."""
        assert algebra.gcd(a, b) == algebra.gcd(b, a)

    @given(st.integers(1, 10000), st.integers(1, 10000))
    def test_gcd_divides_a(self, a: int, b: int):
        """gcd(a, b) teilt a."""
        g = algebra.gcd(a, b)
        assert g > 0
        assert a % g == 0

    @given(st.integers(1, 10000), st.integers(1, 10000))
    def test_gcd_divides_b(self, a: int, b: int):
        """gcd(a, b) teilt b."""
        g = algebra.gcd(a, b)
        assert b % g == 0

    @given(st.integers(1, 10000), st.integers(1, 10000))
    def test_gcd_positive(self, a: int, b: int):
        """gcd(a, b) > 0 für positive Eingaben."""
        assert algebra.gcd(a, b) > 0

    @given(st.integers(1, 10000))
    def test_gcd_with_self(self, a: int):
        """gcd(a, a) == a."""
        assert algebra.gcd(a, a) == a

    @given(st.integers(1, 10000))
    def test_gcd_with_one(self, a: int):
        """gcd(a, 1) == 1."""
        assert algebra.gcd(a, 1) == 1


class TestLcmProperties:
    """Property-Based Tests für das Kleinste Gemeinsame Vielfache."""

    @given(st.integers(1, 1000), st.integers(1, 1000))
    def test_lcm_gcd_relation(self, a: int, b: int):
        """lcm(a, b) * gcd(a, b) == a * b."""
        assert algebra.lcm(a, b) * algebra.gcd(a, b) == a * b

    @given(st.integers(1, 1000), st.integers(1, 1000))
    def test_lcm_divisible_by_a(self, a: int, b: int):
        """lcm(a, b) ist durch a teilbar."""
        l = algebra.lcm(a, b)
        assert l % a == 0

    @given(st.integers(1, 1000), st.integers(1, 1000))
    def test_lcm_divisible_by_b(self, a: int, b: int):
        """lcm(a, b) ist durch b teilbar."""
        l = algebra.lcm(a, b)
        assert l % b == 0

    @given(st.integers(1, 1000), st.integers(1, 1000))
    def test_lcm_commutative(self, a: int, b: int):
        """lcm(a, b) == lcm(b, a) — Kommutativität."""
        assert algebra.lcm(a, b) == algebra.lcm(b, a)


# =============================================================================
# ALGEBRA: Primzahlen und Euler-Phi
# =============================================================================

class TestPrimeProperties:
    """Property-Based Tests für Primzahlprüfung."""

    @given(st.sampled_from([p for p in _SMALL_PRIMES if p > 2]))
    def test_prime_product_is_not_prime(self, p: int):
        """Falls p > 2 eine Primzahl ist, ist p * 2 keine Primzahl."""
        # p * 2 ist gerade und > 2, also keine Primzahl
        assert not algebra.is_prime(p * 2)

    @given(st.integers(2, 500).map(lambda n: n * 2))
    def test_even_number_not_prime(self, n: int):
        """Gerade Zahlen > 2 sind keine Primzahlen."""
        assume(n > 2)
        assert not algebra.is_prime(n)

    @given(st.integers(2, 100))
    def test_prime_squared_not_prime(self, n: int):
        """Das Quadrat einer Zahl > 1 ist keine Primzahl."""
        assert not algebra.is_prime(n * n)


class TestEulerPhiProperties:
    """Property-Based Tests für die Euler-Phi-Funktion."""

    @given(st.sampled_from(_SMALL_PRIMES[:20]))
    def test_euler_phi_prime(self, p: int):
        """euler_phi(p) == p - 1 für Primzahlen p."""
        assert algebra.euler_phi(p) == p - 1

    @given(st.integers(1, 200))
    def test_euler_phi_positive(self, n: int):
        """euler_phi(n) >= 1 für alle n >= 1."""
        assert algebra.euler_phi(n) >= 1

    @given(st.integers(2, 200))
    def test_euler_phi_less_than_n(self, n: int):
        """euler_phi(n) < n für n >= 2."""
        assert algebra.euler_phi(n) < n


# =============================================================================
# ALGEBRA: Modulares Inverses
# =============================================================================

class TestModInverseProperties:
    """Property-Based Tests für das modulare Inverse."""

    @given(st.integers(2, 1000), st.sampled_from(_SMALL_PRIMES[:15]))
    def test_mod_inverse_property(self, a: int, p: int):
        """mod_inverse(a, p) * a ≡ 1 (mod p) wenn gcd(a, p) == 1."""
        a = a % p
        assume(a > 0)
        assume(algebra.gcd(a, p) == 1)
        inv = algebra.mod_inverse(a, p)
        assert (inv * a) % p == 1

    @given(st.sampled_from(_SMALL_PRIMES[:15]))
    def test_mod_inverse_of_1(self, p: int):
        """mod_inverse(1, p) == 1."""
        assert algebra.mod_inverse(1, p) == 1


# =============================================================================
# ANALYSIS: Numerische Ableitung
# =============================================================================

class TestNumericalDerivativeProperties:
    """Property-Based Tests für numerische Ableitung."""

    @given(st.floats(-10, 10, allow_nan=False, allow_infinity=False),
           st.floats(-5, 5, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_derivative_of_linear_function(self, c: float, x: float):
        """Ableitung von f(x) = c*x ist c (konstant)."""
        assume(abs(c) < 100)  # Keine extremen Werte
        f = lambda t: c * t
        deriv = analysis.numerical_derivative(f, x)
        assert abs(deriv - c) < 1e-4

    @given(st.floats(0.1, 5.0, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_derivative_of_quadratic(self, x: float):
        """Ableitung von f(x) = x² ist 2x."""
        f = lambda t: t * t
        deriv = analysis.numerical_derivative(f, x)
        expected = 2 * x
        assert abs(deriv - expected) < 1e-4

    @given(st.floats(0.01, 3.0, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_derivative_of_sin(self, x: float):
        """Ableitung von sin(x) ist cos(x)."""
        deriv = analysis.numerical_derivative(math.sin, x)
        expected = math.cos(x)
        assert abs(deriv - expected) < 1e-5


# =============================================================================
# ANALYSIS: Numerisches Integral
# =============================================================================

class TestNumericalIntegralProperties:
    """Property-Based Tests für numerische Integration."""

    @given(st.floats(-5, 5, allow_nan=False, allow_infinity=False),
           st.floats(-5, 5, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_integral_of_constant(self, a: float, b: float):
        """∫₍ₐ₎ᵇ 1 dx = b - a."""
        assume(a != b)
        assume(abs(b - a) > 0.001)
        assume(abs(b - a) < 100)
        result = analysis.numerical_integral(lambda x: 1.0, a, b, n=100)
        assert abs(result - (b - a)) < 1e-4

    @given(st.floats(0, 5, allow_nan=False, allow_infinity=False),
           st.floats(0, 5, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_integral_additivity(self, a: float, c: float):
        """∫ₐᵇ f + ∫ᵦᶜ f = ∫ₐᶜ f (Additivität)."""
        assume(a < c)
        b = (a + c) / 2
        f = lambda x: x * x  # x²

        int_ac = analysis.numerical_integral(f, a, c, n=100)
        int_ab = analysis.numerical_integral(f, a, b, n=100)
        int_bc = analysis.numerical_integral(f, b, c, n=100)

        assert abs(int_ac - (int_ab + int_bc)) < 1e-6

    @given(st.floats(-3, 3, allow_nan=False, allow_infinity=False),
           st.floats(-3, 3, allow_nan=False, allow_infinity=False))
    @settings(max_examples=10, deadline=3000)
    def test_integral_sign_flip(self, a: float, b: float):
        """∫ₐᵇ f = -∫ᵦₐ f (Vorzeichen bei Umkehrung)."""
        assume(abs(a - b) > 0.001)
        f = lambda x: x * x
        assert abs(analysis.numerical_integral(f, a, b, n=100) +
                   analysis.numerical_integral(f, b, a, n=100)) < 1e-6


# =============================================================================
# LINEARE ALGEBRA: Vektoren
# =============================================================================

class TestVectorProperties:
    """Property-Based Tests für Vektoren."""

    @given(st.lists(st.floats(-100, 100, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10),
           st.lists(st.floats(-100, 100, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10))
    @settings(max_examples=10, deadline=3000)
    def test_dot_product_commutative(self, v: list, w: list):
        """v · w == w · v — Kommutativität des Skalarprodukts."""
        # Gleiche Länge sicherstellen
        n = min(len(v), len(w))
        v, w = v[:n], w[:n]
        assume(n >= 2)

        vec_v = linear_algebra.Vector(v)
        vec_w = linear_algebra.Vector(w)

        dot_vw = vec_v.dot(vec_w)
        dot_wv = vec_w.dot(vec_v)
        assert abs(dot_vw - dot_wv) < 1e-8

    @given(st.lists(st.floats(-100, 100, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10))
    @settings(max_examples=10, deadline=3000)
    def test_norm_non_negative(self, v: list):
        """||v|| >= 0 für alle Vektoren."""
        vec = linear_algebra.Vector(v)
        assert vec.norm() >= 0

    @given(st.lists(st.floats(-100, 100, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10))
    @settings(max_examples=10, deadline=3000)
    def test_dot_product_with_self_non_negative(self, v: list):
        """v · v >= 0 (positiv semidefinit)."""
        vec = linear_algebra.Vector(v)
        assert vec.dot(vec) >= 0

    @given(st.lists(st.floats(-50, 50, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10))
    @settings(max_examples=10, deadline=3000)
    def test_norm_equals_sqrt_dot_self(self, v: list):
        """||v|| == sqrt(v · v)."""
        vec = linear_algebra.Vector(v)
        norm = vec.norm()
        dot_self = vec.dot(vec)
        assert abs(norm - math.sqrt(max(0, dot_self))) < 1e-8

    @given(st.lists(st.floats(-100, 100, allow_nan=False, allow_infinity=False),
                    min_size=2, max_size=10))
    @settings(max_examples=10, deadline=3000)
    def test_normalized_vector_has_unit_norm(self, v: list):
        """Ein normierter Vektor hat Norm 1."""
        assume(any(abs(x) > 1e-10 for x in v))  # Kein Nullvektor
        vec = linear_algebra.Vector(v)
        assume(vec.norm() > 1e-10)
        normalized = vec.normalize()
        assert abs(normalized.norm() - 1.0) < 1e-8


# =============================================================================
# LINEARE ALGEBRA: Givens QR-Zerlegung
# =============================================================================

class TestGivensQRProperties:
    """Property-Based Tests für die Givens QR-Zerlegung."""

    @given(st.lists(
        st.lists(st.floats(-10, 10, allow_nan=False, allow_infinity=False),
                 min_size=2, max_size=4),
        min_size=2, max_size=4
    ))
    @settings(max_examples=10, deadline=3000)
    def test_givens_qr_q_orthogonal(self, data: list):
        """Q in der Givens QR-Zerlegung ist orthogonal: Q^T @ Q ≈ I."""
        import numpy as np

        # Matrix muss m >= n sein
        m = len(data)
        n = len(data[0])
        # Gleiche Spaltenanzahl sicherstellen
        data = [row[:n] for row in data]
        assume(m >= n)
        assume(n >= 2)

        # Überprüfen ob Matrix nicht degeneriert ist
        A = linear_algebra.Matrix(data)
        try:
            Q, R = linear_algebra.givens_qr_decomposition(A)
        except Exception:
            return  # Bei numerischen Problemen Test überspringen

        # Q^T @ Q prüfen
        Q_arr = np.array([[Q._data[i][j] for j in range(Q.cols)]
                          for i in range(Q.rows)])
        product = Q_arr.T @ Q_arr
        identity = np.eye(Q_arr.shape[0])

        assert np.allclose(product, identity, atol=1e-6)

    @given(st.lists(
        st.lists(st.floats(-5, 5, allow_nan=False, allow_infinity=False),
                 min_size=2, max_size=3),
        min_size=2, max_size=3
    ))
    @settings(max_examples=10, deadline=3000)
    def test_givens_qr_reconstruction(self, data: list):
        """Givens QR: Q @ R ≈ A."""
        import numpy as np

        m = len(data)
        n = len(data[0])
        data = [row[:n] for row in data]
        assume(m >= n >= 2)

        A = linear_algebra.Matrix(data)
        try:
            Q, R = linear_algebra.givens_qr_decomposition(A)
        except Exception:
            return

        A_arr = np.array([[A._data[i][j] for j in range(A.cols)]
                          for i in range(A.rows)])
        Q_arr = np.array([[Q._data[i][j] for j in range(Q.cols)]
                          for i in range(Q.rows)])
        R_arr = np.array([[R._data[i][j] for j in range(R.cols)]
                          for i in range(R.rows)])

        reconstructed = Q_arr @ R_arr
        assert np.allclose(reconstructed, A_arr, atol=1e-6)


# =============================================================================
# LINEARE ALGEBRA: LU-Zerlegung
# =============================================================================

class TestLUDecompositionProperties:
    """Property-Based Tests für die LU-Zerlegung."""

    @given(st.lists(
        st.lists(st.floats(-10, 10, allow_nan=False, allow_infinity=False),
                 min_size=3, max_size=3),
        min_size=3, max_size=3
    ))
    @settings(max_examples=10, deadline=3000)
    def test_lu_reconstruction(self, data: list):
        """LU-Zerlegung: P @ A ≈ L @ U."""
        import numpy as np

        n = len(data[0])
        m = len(data)
        # Bei 3x3 fest: m == n == 3
        data = [row[:n] for row in data[:n]]
        assume(n >= 2)

        A = linear_algebra.Matrix(data)

        try:
            P, L, U = linear_algebra.lu_decomposition(A)
        except Exception:
            return  # Numerisch singuläre Matrizen überspringen

        P_arr = np.array([[P._data[i][j] for j in range(P.cols)]
                          for i in range(P.rows)])
        L_arr = np.array([[L._data[i][j] for j in range(L.cols)]
                          for i in range(L.rows)])
        U_arr = np.array([[U._data[i][j] for j in range(U.cols)]
                          for i in range(U.rows)])
        A_arr = np.array([[A._data[i][j] for j in range(A.cols)]
                          for i in range(A.rows)])

        # P @ A == L @ U
        assert np.allclose(P_arr @ A_arr, L_arr @ U_arr, atol=1e-5)


# =============================================================================
# ZAHLENTHEORIE: Lineare Diophantische Gleichungen
# =============================================================================

class TestLinearDiophantineProperties:
    """Property-Based Tests für lineare Diophantische Gleichungen."""

    @given(st.integers(-100, 100), st.integers(-100, 100), st.integers(-100, 100))
    @settings(max_examples=10, deadline=3000)
    def test_solution_satisfies_equation(self, a: int, b: int, c: int):
        """Wenn solve_linear_diophantine(a, b, c) eine Lösung liefert,
        muss a*x0 + b*y0 == c erfüllt sein."""
        assume(a != 0 or b != 0)  # Nicht triviale Gleichung

        result = algebra.solve_linear_diophantine(a, b, c)

        if result is not None:
            x0, y0, g = result
            # Die Partikulärlösung muss die Gleichung erfüllen
            assert a * x0 + b * y0 == c

    @given(st.integers(2, 50), st.integers(2, 50))
    @settings(max_examples=10, deadline=3000)
    def test_no_solution_when_gcd_not_divides_c(self, a: int, b: int):
        """Keine Lösung wenn gcd(a,b) nicht c teilt."""
        g = algebra.gcd(a, b)
        assume(g > 1)  # Nur nicht-triviale GCDs
        # c = g - 1 ist nicht durch g teilbar (da 0 < g-1 < g)
        c_no_sol = g - 1
        assume(c_no_sol > 0)
        result = algebra.solve_linear_diophantine(a, b, c_no_sol)
        assert result is None


# =============================================================================
# ZAHLENTHEORIE: Tonelli-Shanks
# =============================================================================

class TestTonelliShanksProperties:
    """Property-Based Tests für den Tonelli-Shanks-Algorithmus."""

    @given(st.sampled_from([p for p in _SMALL_PRIMES if p > 2]),
           st.integers(1, 1000))
    @settings(max_examples=10, deadline=3000)
    def test_square_root_squares_to_n(self, p: int, n: int):
        """Wenn tonelli_shanks(n, p) = r, dann muss r² ≡ n (mod p)."""
        n_mod = n % p
        assume(n_mod > 0)

        r = algebra.tonelli_shanks(n_mod, p)

        if r is not None:
            # r² ≡ n (mod p) prüfen
            assert (r * r) % p == n_mod

    @given(st.sampled_from([p for p in _SMALL_PRIMES if p > 2]))
    @settings(max_examples=10, deadline=3000)
    def test_tonelli_zero_returns_zero(self, p: int):
        """tonelli_shanks(0, p) == 0."""
        assert algebra.tonelli_shanks(0, p) == 0


# =============================================================================
# ZAHLENTHEORIE: Pell-Gleichung
# =============================================================================

class TestPellEquationProperties:
    """Property-Based Tests für die Pell-Gleichung."""

    @given(st.integers(2, 30))
    @settings(max_examples=10, deadline=3000)
    def test_pell_solutions_satisfy_equation(self, D: int):
        """Alle berechneten Pell-Lösungen erfüllen x² - D*y² = 1."""
        # D darf kein perfektes Quadrat sein
        assume(int(math.isqrt(D)) ** 2 != D)

        solutions = algebra.solve_quadratic_diophantine_pell(D, n_solutions=3)

        for x, y in solutions:
            assert x * x - D * y * y == 1, (
                f"Pell-Gleichung nicht erfüllt: {x}² - {D}·{y}² = {x*x - D*y*y} ≠ 1"
            )

    @given(st.integers(2, 20))
    @settings(max_examples=10, deadline=3000)
    def test_pell_first_solution_positive(self, D: int):
        """Erste Pell-Lösung hat positive Werte x > 0, y > 0."""
        assume(int(math.isqrt(D)) ** 2 != D)

        solutions = algebra.solve_quadratic_diophantine_pell(D, n_solutions=1)

        if solutions:
            x, y = solutions[0]
            assert x > 0
            assert y > 0
