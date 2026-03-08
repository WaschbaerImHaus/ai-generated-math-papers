"""
@file test_p_adic.py
@brief Tests für das p-adische Zahlen-Modul (p_adic.py).
@description
    Testet:
    - p_adic_valuation
    - p_adic_norm
    - p_adic_distance
    - PAdicNumber (from_integer, from_fraction, __add__, __mul__, to_integer_approx, norm)
    - hensel_lift
    - p_adic_exp, p_adic_log
    - ostrowski_theorem_demo

@author Kurt Ingwer
@version 1.0
@since 2026-03-08
@lastModified 2026-03-08
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from p_adic import (
    p_adic_valuation,
    p_adic_norm,
    p_adic_distance,
    PAdicNumber,
    hensel_lift,
    p_adic_exp,
    p_adic_log,
    ostrowski_theorem_demo,
)


# ===========================================================================
# TESTS FÜR p_adic_valuation
# ===========================================================================

class TestPAdicValuation:
    """Tests für die p-adische Bewertung v_p(n)."""

    def test_v2_of_12_is_2(self):
        """v_2(12) = 2, da 12 = 4·3 = 2²·3."""
        assert p_adic_valuation(12, 2) == 2

    def test_v3_of_12_is_1(self):
        """v_3(12) = 1, da 12 = 4·3 = 2²·3¹."""
        assert p_adic_valuation(12, 3) == 1

    def test_v5_of_12_is_0(self):
        """v_5(12) = 0, da 5 ∤ 12."""
        assert p_adic_valuation(12, 5) == 0

    def test_v2_of_0_is_infinity(self):
        """v_p(0) = +∞ per Konvention."""
        assert p_adic_valuation(0, 2) == float('inf')

    def test_v2_of_1_is_0(self):
        """v_2(1) = 0, da 2 ∤ 1."""
        assert p_adic_valuation(1, 2) == 0

    def test_v2_of_8_is_3(self):
        """v_2(8) = 3, da 8 = 2³."""
        assert p_adic_valuation(8, 2) == 3

    def test_v7_of_343_is_3(self):
        """v_7(343) = 3, da 343 = 7³."""
        assert p_adic_valuation(343, 7) == 3

    def test_negative_n_same_as_positive(self):
        """v_p(-n) = v_p(n)."""
        assert p_adic_valuation(-12, 2) == p_adic_valuation(12, 2)

    def test_v2_of_prime_is_0(self):
        """v_2(p) = 0 für ungerade Primzahl p."""
        for p_prime in [3, 5, 7, 11, 13]:
            assert p_adic_valuation(p_prime, 2) == 0

    def test_raises_for_p_less_than_2(self):
        """Fehler wenn p < 2."""
        with pytest.raises(ValueError):
            p_adic_valuation(12, 1)

    def test_v_of_power_of_p(self):
        """v_p(p^k) = k."""
        for p in [2, 3, 5, 7]:
            for k in range(1, 6):
                assert p_adic_valuation(p**k, p) == k

    def test_additivity(self):
        """v_p(mn) = v_p(m) + v_p(n) für m, n ≠ 0."""
        for m, n in [(12, 18), (8, 9), (25, 7), (4, 6)]:
            for p in [2, 3]:
                v_mn = p_adic_valuation(m * n, p)
                v_m = p_adic_valuation(m, p)
                v_n = p_adic_valuation(n, p)
                assert v_mn == v_m + v_n, f"Additivität verletzt: v_{p}({m}·{n})"


# ===========================================================================
# TESTS FÜR p_adic_norm
# ===========================================================================

class TestPAdicNorm:
    """Tests für den p-adischen Betrag |n|_p."""

    def test_norm_12_base2_is_quarter(self):
        """|12|_2 = 2^{-2} = 0.25."""
        assert p_adic_norm(12, 2) == pytest.approx(0.25)

    def test_norm_12_base3(self):
        """|12|_3 = 3^{-1} = 1/3."""
        assert p_adic_norm(12, 3) == pytest.approx(1.0/3.0)

    def test_norm_0_is_0(self):
        """|0|_p = 0 per Definition."""
        assert p_adic_norm(0, 2) == 0.0
        assert p_adic_norm(0, 5) == 0.0

    def test_norm_1_is_1(self):
        """|1|_p = 1 für alle p."""
        for p in [2, 3, 5, 7]:
            assert p_adic_norm(1, p) == pytest.approx(1.0)

    def test_norm_multiplicative(self):
        """|mn|_p = |m|_p · |n|_p."""
        for m, n in [(12, 18), (8, 9)]:
            for p in [2, 3]:
                lhs = p_adic_norm(m * n, p)
                rhs = p_adic_norm(m, p) * p_adic_norm(n, p)
                assert lhs == pytest.approx(rhs, rel=1e-10)

    def test_norm_in_range_0_to_1_for_integers(self):
        """Für ganze Zahlen n ≥ 1 gilt: |n|_p ≤ 1."""
        for n in range(1, 50):
            for p in [2, 3, 5]:
                norm = p_adic_norm(n, p)
                assert 0 <= norm <= 1, f"|{n}|_{p} = {norm} nicht in [0,1]"

    def test_ultrametric_inequality(self):
        """|m+n|_p ≤ max(|m|_p, |n|_p)."""
        test_cases = [(3, 5, 2), (12, 8, 2), (9, 6, 3)]
        for m, n, p in test_cases:
            lhs = p_adic_norm(m + n, p)
            rhs = max(p_adic_norm(m, p), p_adic_norm(n, p))
            assert lhs <= rhs + 1e-12, f"Ultrametrik verletzt: |{m}+{n}|_{p}"


# ===========================================================================
# TESTS FÜR p_adic_distance
# ===========================================================================

class TestPAdicDistance:
    """Tests für den p-adischen Abstand d_p(x, y)."""

    def test_distance_symmetric(self):
        """d_p(x, y) = d_p(y, x)."""
        for x, y, p in [(12, 20, 2), (9, 18, 3), (25, 100, 5)]:
            assert p_adic_distance(x, y, p) == pytest.approx(p_adic_distance(y, x, p))

    def test_distance_self_is_zero(self):
        """d_p(x, x) = 0."""
        for x, p in [(12, 2), (100, 3), (0, 5)]:
            assert p_adic_distance(x, x, p) == 0.0

    def test_distance_0_and_p_power(self):
        """d_p(0, p^k) = p^{-k} (p^k ist nah bei 0 in p-adischer Metrik)."""
        for p in [2, 3, 5]:
            for k in [1, 2, 3]:
                dist = p_adic_distance(0, p**k, p)
                expected = float(p) ** (-k)
                assert dist == pytest.approx(expected)

    def test_ultrametric_triangle_inequality(self):
        """d_p(x, z) ≤ max(d_p(x, y), d_p(y, z))."""
        x, y, z, p = 0, 8, 12, 2
        d_xz = p_adic_distance(x, z, p)
        d_xy = p_adic_distance(x, y, p)
        d_yz = p_adic_distance(y, z, p)
        assert d_xz <= max(d_xy, d_yz) + 1e-12


# ===========================================================================
# TESTS FÜR PAdicNumber
# ===========================================================================

class TestPAdicNumber:
    """Tests für die PAdicNumber-Klasse."""

    def test_from_integer_zero(self):
        """PAdicNumber.from_integer(0, p) ergibt alle Nullziffern."""
        x = PAdicNumber.from_integer(0, 2)
        assert all(d == 0 for d in x.digits)

    def test_from_integer_one(self):
        """PAdicNumber.from_integer(1, p): Erste Ziffer ist 1, Rest 0."""
        x = PAdicNumber.from_integer(1, 2, precision=8)
        assert x.digits[0] == 1
        assert all(d == 0 for d in x.digits[1:])

    def test_from_integer_four_in_base2(self):
        """4 in Basis 2: Ziffern sind [0, 0, 1, 0, 0, ...]."""
        x = PAdicNumber.from_integer(4, 2, precision=8)
        assert x.digits[0] == 0  # 4 mod 2 = 0
        assert x.digits[1] == 0  # (4//2) mod 2 = 0
        assert x.digits[2] == 1  # (4//4) mod 2 = 1

    def test_roundtrip_positive_integer(self):
        """from_integer(n) → to_integer_approx() ≈ n."""
        for n in [1, 5, 13, 42, 127]:
            for p in [2, 3, 5]:
                x = PAdicNumber.from_integer(n, p, precision=20)
                recovered = x.to_integer_approx(terms=10)
                # Roundtrip: recovered ≡ n (mod p^10)
                assert recovered % (p**10) == n % (p**10), (
                    f"Roundtrip fehlgeschlagen für n={n}, p={p}: {recovered}"
                )

    def test_from_fraction_half_in_3adic(self):
        """
        1/2 in 3-adischen Zahlen: 2 · (1/2) ≡ 1 (mod 3^k).
        Prüfe, dass 2 · from_fraction(1,2,3) ≡ 1 (mod 3^5).
        """
        half = PAdicNumber.from_fraction(1, 2, 3, precision=10)
        # Multipliziere mit 2 (als from_integer(2, 3))
        two = PAdicNumber.from_integer(2, 3, precision=10)
        product = half * two
        # Das Produkt sollte ≡ 1 (mod 3^5) sein
        assert product.to_integer_approx(5) % (3**5) == 1 % (3**5)

    def test_add_integers_base2(self):
        """3 + 5 = 8 in 2-adischen Zahlen."""
        three = PAdicNumber.from_integer(3, 2, precision=10)
        five = PAdicNumber.from_integer(5, 2, precision=10)
        eight = PAdicNumber.from_integer(8, 2, precision=10)
        result = three + five
        # Vergleiche mod 2^8
        assert result.to_integer_approx(8) % (2**8) == 8 % (2**8)

    def test_add_different_p_raises(self):
        """Addition mit unterschiedlichen Primbasen wirft Fehler."""
        x = PAdicNumber([1, 0, 1], p=2)
        y = PAdicNumber([1, 0, 1], p=3)
        with pytest.raises(ValueError):
            x + y

    def test_mul_integers_base3(self):
        """2 · 4 = 8 in 3-adischen Zahlen."""
        two = PAdicNumber.from_integer(2, 3, precision=10)
        four = PAdicNumber.from_integer(4, 3, precision=10)
        result = two * four
        assert result.to_integer_approx(5) % (3**5) == 8 % (3**5)

    def test_norm_of_power_of_p(self):
        """Norm von p^k ist p^{-k}."""
        for p in [2, 3]:
            for k in [1, 2, 3]:
                x = PAdicNumber.from_integer(p**k, p, precision=15)
                norm = x.norm()
                expected = float(p) ** (-k)
                assert norm == pytest.approx(expected, rel=0.01), (
                    f"Norm von {p}^{k} in Basis {p}: {norm} ≠ {expected}"
                )

    def test_invalid_digit_raises(self):
        """Ziffer außerhalb [0, p-1] wirft ValueError."""
        with pytest.raises(ValueError):
            PAdicNumber([0, 2, 1], p=2)  # 2 ist keine gültige 2-adische Ziffer

    def test_repr_and_str(self):
        """__repr__ und __str__ liefern nicht-leere Strings."""
        x = PAdicNumber.from_integer(5, 2)
        assert len(repr(x)) > 0
        assert len(str(x)) > 0

    def test_to_integer_approx_zero(self):
        """to_integer_approx(0) = 0."""
        x = PAdicNumber.from_integer(0, 2)
        assert x.to_integer_approx(0) == 0


# ===========================================================================
# TESTS FÜR hensel_lift
# ===========================================================================

class TestHenselLift:
    """Tests für Hensels Lemma."""

    def test_sqrt_minus1_mod5(self):
        """
        x² ≡ -1 (mod 5): Startwurzel x₀ = 2 (da 2² = 4 ≡ -1 mod 5).
        Hebe auf zu x² ≡ -1 (mod 5^3).
        """
        # f(x) = x² + 1, Koeffizienten [1, 0, 1] (höchster Grad zuerst)
        root = hensel_lift([1, 0, 1], p=5, initial_root=2, n_lifts=2)
        # Prüfe: root² ≡ -1 (mod 5^3 = 125)
        assert (root ** 2 + 1) % 125 == 0, (
            f"Hensel-Lift fehlgeschlagen: {root}² + 1 ≡ {(root**2+1)%125} (mod 125)"
        )

    def test_sqrt_minus1_mod5_other_root(self):
        """
        x² ≡ -1 (mod 5): Andere Startwurzel x₀ = 3 (3² = 9 ≡ -1 mod 5).
        """
        root = hensel_lift([1, 0, 1], p=5, initial_root=3, n_lifts=2)
        assert (root ** 2 + 1) % 125 == 0

    def test_simple_linear_equation(self):
        """
        f(x) = x - 1: Wurzel x₀ = 1 mod p.
        Lift liefert x = 1 + k·p für passendes k.
        f = [1, -1], also x - 1 = 0 → x = 1.
        """
        root = hensel_lift([1, -1], p=3, initial_root=1, n_lifts=3)
        # f(1) = 1 - 1 = 0 mod 3^k für alle k
        assert (root - 1) % (3**4) == 0

    def test_raises_for_invalid_start(self):
        """Fehler wenn f(initial_root) ≢ 0 (mod p)."""
        # f(x) = x² + 1, start = 1: f(1) = 2 ≢ 0 (mod 5)
        with pytest.raises(ValueError):
            hensel_lift([1, 0, 1], p=5, initial_root=1, n_lifts=2)

    def test_cubic_equation(self):
        """
        f(x) = x³ - 2 modulo 7.
        3³ = 27 ≡ 6 ≡ -1 (mod 7) – kein offensichtlicher Kandidat.
        Suche eine Kubikwurzel von 2 mod 7: 4³ = 64 ≡ 1 (mod 7) – nein.
        Verwende f(x) = x² - 2: x=3, 3²=9≡2 (mod 7). Lift auf mod 49.
        """
        root = hensel_lift([1, 0, -2], p=7, initial_root=3, n_lifts=1)
        # Prüfe: root² ≡ 2 (mod 7²=49)
        assert (root**2 - 2) % 49 == 0

    def test_singular_root_raises(self):
        """Fehler wenn f'(root) ≡ 0 (mod p) (singuläre Wurzel)."""
        # f(x) = x², f'(x) = 2x. Bei x₀ = 0: f(0) = 0 ✓, f'(0) = 0 → Fehler
        with pytest.raises(ValueError):
            hensel_lift([1, 0, 0], p=3, initial_root=0, n_lifts=2)


# ===========================================================================
# TESTS FÜR ostrowski_theorem_demo
# ===========================================================================

class TestOstrowskiTheoremDemo:
    """Tests für die Demonstration des Satzes von Ostrowski."""

    def test_product_formula_n12(self):
        """Produktformel für n=12: |12|_∞ · Π_p |12|_p = 1."""
        result = ostrowski_theorem_demo(12)
        assert result['product_equals_one'] == True

    def test_product_formula_n1(self):
        """Produktformel für n=1: |1|_∞ = 1, alle |1|_p = 1, Produkt = 1."""
        result = ostrowski_theorem_demo(1)
        assert result['product_equals_one'] == True

    def test_product_formula_n_prime(self):
        """Produktformel für Primzahlen: |p|_∞ · |p|_p = p · (1/p) = 1."""
        for p in [2, 3, 5, 7, 11]:
            result = ostrowski_theorem_demo(p)
            assert result['product_equals_one'] == True, (
                f"Produktformel verletzt für n={p}: {result['product_formula_exact']}"
            )

    def test_abs_infinity_correct(self):
        """|n|_∞ = |n| (gewöhnlicher absoluter Wert)."""
        for n in [1, -5, 12, -100, 36]:
            result = ostrowski_theorem_demo(n)
            assert result['abs_infinity'] == abs(n)

    def test_raises_for_zero(self):
        """Fehler für n = 0."""
        with pytest.raises(ValueError):
            ostrowski_theorem_demo(0)

    def test_p_adic_norms_correct_for_n12(self):
        """Für n=12: |12|_2 = 0.25, |12|_3 = 1/3."""
        result = ostrowski_theorem_demo(12)
        assert result['p_adic_norms'][2] == pytest.approx(0.25)
        assert result['p_adic_norms'][3] == pytest.approx(1.0/3.0)
        assert result['p_adic_norms'][5] == pytest.approx(1.0)

    def test_returns_prime_factorization(self):
        """Rückgabe enthält korrekte Primfaktorzerlegung."""
        result = ostrowski_theorem_demo(12)
        assert result['prime_factorization'] == {2: 2, 3: 1}

    def test_n_is_stored(self):
        """Das n wird korrekt gespeichert."""
        result = ostrowski_theorem_demo(42)
        assert result['n'] == 42

    def test_product_formula_large_n(self):
        """Produktformel gilt auch für größere Zahlen."""
        for n in [60, 100, 360]:
            result = ostrowski_theorem_demo(n)
            assert result['product_equals_one'] == True, (
                f"Produktformel verletzt für n={n}"
            )

    def test_negative_n(self):
        """Produktformel gilt auch für negative n (|n|_∞ = |-n|_∞ = |n|)."""
        result = ostrowski_theorem_demo(-12)
        assert result['product_equals_one'] == True
        assert result['abs_infinity'] == 12


# ===========================================================================
# TESTS FÜR p_adic_exp UND p_adic_log
# ===========================================================================

class TestPAdicExpLog:
    """Tests für p-adische Exponential- und Logarithmusfunktion."""

    def test_exp_of_zero(self):
        """exp_p(0) = 1."""
        # 0 hat alle Ziffern = 0
        zero_digits = [0] * 8
        result = p_adic_exp(zero_digits, p=5, precision=5)
        # Ergebnis sollte 1 entsprechen: erste Ziffer 1, Rest 0
        assert result[0] == 1
        assert all(d == 0 for d in result[1:])

    def test_exp_returns_list_of_correct_length(self):
        """exp_p gibt Liste der Länge precision zurück."""
        digits = [0, 1] + [0] * 6  # x = p (Bewertung 1)
        result = p_adic_exp(digits, p=5, precision=6)
        assert len(result) == 6

    def test_log_of_one(self):
        """log_p(1) = 0."""
        # 1 hat Ziffern [1, 0, 0, ...]
        one_digits = [1] + [0] * 7
        result = p_adic_log(one_digits, p=5, precision=5)
        # Ergebnis sollte 0 sein
        assert result[0] == 0

    def test_log_returns_list_of_correct_length(self):
        """log_p gibt Liste der Länge precision zurück."""
        one_digits = [1] + [0] * 7
        result = p_adic_log(one_digits, p=3, precision=5)
        assert len(result) == 5
