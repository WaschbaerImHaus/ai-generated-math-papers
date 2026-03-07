"""
@file test_algebra.py
@brief Tests für das Algebra-Modul (Test-Driven Development).
@description
    Testet alle Funktionen des Algebra-Moduls:
    - Polynomoperationen
    - Gleichungslöser
    - Faktorisierung
    - GCD/LCM
    - Modulare Arithmetik
@author Reisen macht Spass... mit Pia und Dirk e.Kfm.
@date 2026-03-05
"""

import sys
import os
import math

# Suchpfad für src-Modul setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from algebra import (
    Polynomial, solve_linear, solve_quadratic,
    gcd, lcm, extended_gcd, mod_inverse,
    is_prime, prime_factorization, euler_phi
)


class TestPolynomial:
    """Tests für die Polynomial-Klasse."""

    def test_creation(self):
        """Polynome werden korrekt aus Koeffizientenliste erstellt."""
        p = Polynomial([1, -3, 2])  # x^2 - 3x + 2
        assert p.degree == 2
        assert p.coefficients == [1, -3, 2]

    def test_zero_polynomial(self):
        """Das Nullpolynom hat Grad -1 oder 0."""
        p = Polynomial([0])
        assert p.evaluate(5) == 0

    def test_evaluate(self):
        """Polynomauswertung an bestimmten Punkten."""
        p = Polynomial([1, -3, 2])  # x^2 - 3x + 2
        assert p.evaluate(0) == 2    # f(0) = 2
        assert p.evaluate(1) == 0    # f(1) = 0 (Nullstelle)
        assert p.evaluate(2) == 0    # f(2) = 0 (Nullstelle)
        assert p.evaluate(3) == 2    # f(3) = 9 - 9 + 2 = 2

    def test_addition(self):
        """Addition zweier Polynome."""
        p1 = Polynomial([1, 2, 3])   # x^2 + 2x + 3
        p2 = Polynomial([2, -1, 0])  # 2x^2 - x
        result = p1 + p2             # 3x^2 + x + 3
        assert result.evaluate(1) == 7   # 3*1 + 1 + 3 = 7
        assert result.evaluate(0) == 3

    def test_subtraction(self):
        """Subtraktion zweier Polynome."""
        p1 = Polynomial([1, 2, 3])
        p2 = Polynomial([1, 2, 3])
        result = p1 - p2
        assert result.evaluate(100) == 0

    def test_multiplication(self):
        """Multiplikation zweier Polynome: (x-1)(x-2) = x^2 - 3x + 2."""
        p1 = Polynomial([1, -1])   # x - 1
        p2 = Polynomial([1, -2])   # x - 2
        result = p1 * p2           # x^2 - 3x + 2
        assert result.evaluate(1) == 0   # Nullstelle bei x=1
        assert result.evaluate(2) == 0   # Nullstelle bei x=2
        assert result.evaluate(0) == 2   # f(0) = 2

    def test_derivative(self):
        """Ableitung: (x^3)' = 3x^2."""
        p = Polynomial([1, 0, 0, 0])  # x^3
        dp = p.derivative()
        assert dp.evaluate(1) == 3   # 3*1^2 = 3
        assert dp.evaluate(2) == 12  # 3*4 = 12

    def test_str_representation(self):
        """Polynome werden korrekt als String dargestellt."""
        p = Polynomial([1, -3, 2])
        s = str(p)
        assert 'x' in s


class TestLinearEquations:
    """Tests für lineare Gleichungslöser."""

    def test_simple_linear(self):
        """ax + b = 0 lösen: 2x + 4 = 0 -> x = -2."""
        result = solve_linear(2, 4)
        assert result == -2.0

    def test_linear_no_solution(self):
        """0x + 5 = 0 hat keine Lösung."""
        with pytest.raises(ValueError, match="keine Lösung"):
            solve_linear(0, 5)

    def test_linear_infinite_solutions(self):
        """0x + 0 = 0 hat unendlich viele Lösungen."""
        with pytest.raises(ValueError, match="unendlich"):
            solve_linear(0, 0)


class TestQuadraticEquations:
    """Tests für quadratische Gleichungslöser."""

    def test_two_real_roots(self):
        """x^2 - 5x + 6 = 0 -> x=2 oder x=3."""
        roots = solve_quadratic(1, -5, 6)
        assert sorted(roots) == pytest.approx([2.0, 3.0])

    def test_one_real_root(self):
        """x^2 - 2x + 1 = 0 -> x=1 (doppelte Wurzel)."""
        roots = solve_quadratic(1, -2, 1)
        assert roots == pytest.approx([1.0, 1.0])

    def test_complex_roots(self):
        """x^2 + 1 = 0 -> komplexe Wurzeln."""
        roots = solve_quadratic(1, 0, 1)
        assert len(roots) == 2
        # Komplexe Wurzeln: +i und -i
        for root in roots:
            assert abs(root**2 + 1) < 1e-10


class TestNumberTheory:
    """Tests für zahlentheoretische Funktionen."""

    def test_gcd_basic(self):
        """ggT(12, 8) = 4."""
        assert gcd(12, 8) == 4

    def test_gcd_prime(self):
        """ggT zweier verschiedener Primzahlen = 1."""
        assert gcd(7, 13) == 1

    def test_gcd_with_zero(self):
        """ggT(0, n) = n."""
        assert gcd(0, 5) == 5

    def test_lcm_basic(self):
        """kgV(4, 6) = 12."""
        assert lcm(4, 6) == 12

    def test_extended_gcd(self):
        """Erweiterter euklidischer Algorithmus: ax + by = ggT(a,b)."""
        g, x, y = extended_gcd(35, 15)
        assert g == 5
        assert 35 * x + 15 * y == g

    def test_mod_inverse(self):
        """Modulares Inverses: 3 * mod_inverse(3, 7) ≡ 1 (mod 7)."""
        inv = mod_inverse(3, 7)
        assert (3 * inv) % 7 == 1

    def test_mod_inverse_no_inverse(self):
        """Kein modulares Inverses, wenn ggT(a,m) != 1."""
        with pytest.raises(ValueError):
            mod_inverse(4, 8)  # ggT(4,8)=4 != 1

    def test_is_prime(self):
        """Primzahltest für bekannte Werte."""
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 97]
        composites = [1, 4, 6, 8, 9, 15, 25, 100]
        for p in primes:
            assert is_prime(p), f"{p} sollte eine Primzahl sein"
        for c in composites:
            assert not is_prime(c), f"{c} sollte keine Primzahl sein"

    def test_prime_factorization(self):
        """Primfaktorzerlegung: 360 = 2^3 * 3^2 * 5."""
        factors = prime_factorization(360)
        assert factors == {2: 3, 3: 2, 5: 1}
        # Probe: Produkt ergibt 360
        product = 1
        for prime, exp in factors.items():
            product *= prime ** exp
        assert product == 360

    def test_euler_phi(self):
        """Eulersche Phi-Funktion: phi(12) = 4 (1,5,7,11 sind teilerfremd zu 12)."""
        assert euler_phi(1) == 1
        assert euler_phi(12) == 4
        assert euler_phi(7) == 6   # Primzahl: phi(p) = p-1
        assert euler_phi(36) == 12


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
