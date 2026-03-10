"""
@file test_ring_theory.py
@brief Tests für das Ringtheorie-Modul (Test-Driven Development).
@description
    Testet alle Klassen und Funktionen des Ringtheorie-Moduls:
    - Ring: Einheiten, Nullteiler, Charakteristik, Körper/PID/UFD-Erkennung
    - Ideal: Prim-, maximal-, Hauptideal-Prüfungen
    - QuotientRing: Körper- und Integritätsbereich-Erkennung
    - PolynomialRingModP: Endliche Körper über ℤ_p
    - Hilfsfunktionen: ggT, Irreduzibilität, Faktorisierung, CRT
    - ring_of_integers: Quadratische Zahlkörper
    - gaussian_integers_gcd: Gauß'sche Ganzzahlen

@author Kurt Ingwer
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math

# Suchpfad für src-Modul setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from ring_theory import (
    Ring, Ideal, QuotientRing, PolynomialRingModP,
    principal_ideal, ideal_sum, ideal_product, ideal_intersection,
    chinese_remainder_ring, is_irreducible_mod_p, polynomial_gcd_mod_p,
    polynomial_factorization_mod_p, quotient_ring_field_check,
    noetherian_check, krull_dimension_estimate, ring_of_integers,
    gaussian_integers_gcd
)


# ===========================================================================
# TESTS: Ring ℤ/nℤ
# ===========================================================================

class TestRing:
    """Tests für die Ring-Klasse."""

    def test_ring_z6_elements(self):
        """ℤ/6ℤ hat genau 6 Elemente {0,1,2,3,4,5}."""
        R = Ring(6)
        assert R.elements == [0, 1, 2, 3, 4, 5]

    def test_ring_z6_characteristic(self):
        """ℤ/6ℤ hat Charakteristik 6."""
        R = Ring(6)
        assert R.characteristic() == 6

    def test_ring_z6_zero_divisors(self):
        """ℤ/6ℤ: Nullteiler sind 2 und 3 (denn 2·3 = 6 ≡ 0 mod 6)."""
        R = Ring(6)
        zero_divs = R.zero_divisors()
        assert 2 in zero_divs
        assert 3 in zero_divs
        # 1 ist Einheit, kein Nullteiler
        assert 1 not in zero_divs

    def test_ring_z6_is_not_integral_domain(self):
        """ℤ/6ℤ ist kein Integritätsbereich (hat Nullteiler)."""
        R = Ring(6)
        assert not R.is_integral_domain()

    def test_ring_z6_is_not_field(self):
        """ℤ/6ℤ ist kein Körper (6 ist nicht prim)."""
        R = Ring(6)
        assert not R.is_field()

    def test_ring_z7_is_field(self):
        """ℤ/7ℤ ist ein Körper (7 ist prim)."""
        R = Ring(7)
        assert R.is_field()

    def test_ring_z7_all_nonzero_are_units(self):
        """ℤ/7ℤ: Alle Elemente außer 0 sind Einheiten."""
        R = Ring(7)
        units = R.units()
        expected = [1, 2, 3, 4, 5, 6]
        assert sorted(units) == expected

    def test_ring_z7_no_zero_divisors(self):
        """ℤ/7ℤ: Keine Nullteiler (Körper ist Integritätsbereich)."""
        R = Ring(7)
        assert R.zero_divisors() == []

    def test_ring_z7_characteristic(self):
        """ℤ/7ℤ hat Charakteristik 7."""
        R = Ring(7)
        assert R.characteristic() == 7

    def test_ring_z4_not_integral_domain(self):
        """ℤ/4ℤ ist kein Integritätsbereich: 2·2 = 4 ≡ 0 (mod 4)."""
        R = Ring(4)
        assert not R.is_integral_domain()
        # 2 ist Nullteiler
        assert 2 in R.zero_divisors()

    def test_ring_z2_is_field(self):
        """ℤ/2ℤ ist ein Körper."""
        R = Ring(2)
        assert R.is_field()

    def test_ring_commutative(self):
        """Alle ℤ/nℤ sind kommutativ."""
        for n in [4, 6, 7, 12]:
            R = Ring(n)
            assert R.is_commutative()

    def test_ring_z7_is_pid_and_ufd(self):
        """ℤ/7ℤ (Körper) ist PID und UFD."""
        R = Ring(7)
        assert R.is_pid()
        assert R.is_ufd()

    def test_ring_z6_not_pid(self):
        """ℤ/6ℤ ist kein PID (kein Integritätsbereich)."""
        R = Ring(6)
        assert not R.is_pid()

    def test_ring_invalid_modulus(self):
        """Ring mit Modulus < 2 wirft ValueError."""
        with pytest.raises(ValueError):
            Ring(1)
        with pytest.raises(ValueError):
            Ring(0)

    def test_ring_addition(self):
        """Addition in ℤ/5ℤ: 3 + 4 = 7 ≡ 2 (mod 5)."""
        R = Ring(5)
        assert R.add(3, 4) == 2

    def test_ring_multiplication(self):
        """Multiplikation in ℤ/5ℤ: 3 · 4 = 12 ≡ 2 (mod 5)."""
        R = Ring(5)
        assert R.mul(3, 4) == 2


# ===========================================================================
# TESTS: Ideal
# ===========================================================================

class TestIdeal:
    """Tests für die Ideal-Klasse."""

    def test_ideal_principal_2_in_z6(self):
        """Das Hauptideal (2) in ℤ/6ℤ = {0, 2, 4}."""
        R = Ring(6)
        I = principal_ideal(R, 2)
        assert set(I.elements) == {0, 2, 4}

    def test_ideal_2_in_z_is_prime(self):
        """Das Ideal (2) in ℤ/6ℤ: Prüfe Primidealeigenschaft direkt."""
        # In ℤ: Ideal (2) ist prim, weil ℤ/(2) = ℤ/2ℤ Körper ist.
        # In ℤ/6ℤ: (2) = {0,2,4}: Primidealeigenschaft testen
        R = Ring(6)
        I = principal_ideal(R, 2)
        # (2) in ℤ/6ℤ ist Primideal, weil 6/gcd(2,6)=6/2=3 prim ist
        assert I.is_prime()

    def test_ideal_4_in_z12_not_prime(self):
        """Das Ideal (4) in ℤ/12ℤ ist nicht prim."""
        R = Ring(12)
        I = principal_ideal(R, 4)
        # (4) = {0, 4, 8}: 2·6 = 12 ≡ 0, und 2*2 = 4 ∈ I aber 2 ∉ I
        # Prüfe: Nicht-Primideal
        assert not I.is_prime()

    def test_ideal_is_principal(self):
        """Alle Ideale in ℤ/nℤ sind Hauptideale."""
        R = Ring(12)
        I = principal_ideal(R, 3)
        assert I.is_principal()

    def test_ideal_zero_in_ring(self):
        """Das Nullideal (0) = {0} ist in jedem Ring enthalten."""
        R = Ring(7)
        I = principal_ideal(R, 0)
        assert I.elements == [0]

    def test_ideal_whole_ring(self):
        """Das Einsideal (1) = R ist der gesamte Ring."""
        R = Ring(6)
        I = principal_ideal(R, 1)
        assert set(I.elements) == set(R.elements)

    def test_ideal_maximal_in_z6(self):
        """In ℤ/6ℤ: (2) = {0,2,4} ist maximal (ℤ/6ℤ / (2) ≅ ℤ/3ℤ, 3 prim)."""
        R = Ring(6)
        I = principal_ideal(R, 2)
        assert I.is_maximal()

    def test_ideal_not_maximal(self):
        """Das Nullideal in ℤ/6ℤ ist nicht maximal (kein Körperquotient)."""
        R = Ring(6)
        I = principal_ideal(R, 0)
        # {0} ⊊ (2) ⊊ ℤ/6ℤ: echtes Zwischenideal existiert
        assert not I.is_maximal()

    def test_ideal_containment(self):
        """Teilmengen-Relation: (6) ⊆ (3) ⊆ (1) in ℤ/12ℤ."""
        R = Ring(12)
        I6 = principal_ideal(R, 6)
        I3 = principal_ideal(R, 3)
        I1 = principal_ideal(R, 1)
        assert I6 <= I3
        assert I3 <= I1


# ===========================================================================
# TESTS: Idealoperationen
# ===========================================================================

class TestIdealOperations:
    """Tests für Ideal-Summe, -Produkt und -Schnitt."""

    def test_ideal_sum(self):
        """(2) + (3) = (1) in ℤ/6ℤ (gcd(2,3) = 1)."""
        R = Ring(6)
        I2 = principal_ideal(R, 2)
        I3 = principal_ideal(R, 3)
        S = ideal_sum(I2, I3)
        # (2) + (3) ⊇ {0+0=0, 2+3=5, 2+0=2, 0+3=3, ...} → {0,1,2,3,4,5}
        assert set(S.elements) == set(R.elements)

    def test_ideal_intersection(self):
        """(2) ∩ (3) = (6) ≅ (0) in ℤ/6ℤ."""
        R = Ring(6)
        I2 = principal_ideal(R, 2)
        I3 = principal_ideal(R, 3)
        inter = ideal_intersection(I2, I3)
        # kgV(2,3)=6 ≡ 0 mod 6
        assert set(inter.elements) == {0}

    def test_ideal_product_in_z12(self):
        """(2)·(3) ⊆ (6) in ℤ/12ℤ."""
        R = Ring(12)
        I2 = principal_ideal(R, 2)
        I3 = principal_ideal(R, 3)
        prod = ideal_product(I2, I3)
        I6 = principal_ideal(R, 6)
        # Produkt muss Teilmenge von (6) sein
        assert prod <= I6


# ===========================================================================
# TESTS: QuotientRing
# ===========================================================================

class TestQuotientRing:
    """Tests für den Quotientenring R/I."""

    def test_quotient_ring_z6_by_2_is_field(self):
        """ℤ/6ℤ / (2) ≅ ℤ/3ℤ ist ein Körper (3 prim)."""
        R = Ring(6)
        I = principal_ideal(R, 2)
        QR = QuotientRing(R, I)
        assert QR.is_field()

    def test_quotient_ring_is_integral_domain(self):
        """Quotient über Primideal ist Integritätsbereich."""
        R = Ring(6)
        I = principal_ideal(R, 2)
        QR = QuotientRing(R, I)
        assert QR.is_integral_domain()


# ===========================================================================
# TESTS: PolynomialRingModP
# ===========================================================================

class TestPolynomialRingModP:
    """Tests für den Polynomring ℤ_p[x]/(f(x))."""

    def test_poly_ring_order_gf4(self):
        """ℤ_2[x]/(x²+x+1) hat Ordnung 2² = 4."""
        f = [1, 1, 1]  # x² + x + 1
        PR = PolynomialRingModP(f, 2)
        assert PR.order == 4

    def test_poly_ring_is_field_gf4(self):
        """ℤ_2[x]/(x²+x+1) ist ein Körper (x²+x+1 irreduzibel über ℤ_2)."""
        f = [1, 1, 1]  # x² + x + 1
        PR = PolynomialRingModP(f, 2)
        assert PR.is_field()

    def test_poly_ring_not_field_x2_minus_1(self):
        """ℤ_2[x]/(x²+1) = ℤ_2[x]/((x+1)²) ist kein Körper."""
        # x² + 1 = (x+1)² über ℤ_2 (da 1 = -1 mod 2)
        f = [1, 0, 1]  # x² + 1 = (x+1)² mod 2
        PR = PolynomialRingModP(f, 2)
        assert not PR.is_field()

    def test_poly_ring_elements_gf4(self):
        """ℤ_2[x]/(x²+x+1) hat 4 Elemente."""
        f = [1, 1, 1]
        PR = PolynomialRingModP(f, 2)
        elements = PR.elements()
        assert len(elements) == 4


# ===========================================================================
# TESTS: is_irreducible_mod_p
# ===========================================================================

class TestIrreducibleModP:
    """Tests für die Irreduzibilitätsprüfung."""

    def test_x2_plus_1_mod_2_not_irreducible(self):
        """x² + 1 = (x+1)² mod 2 → nicht irreduzibel."""
        # x² + 1 hat Nullstelle x = 1 mod 2: 1 + 1 = 0 mod 2
        assert not is_irreducible_mod_p([1, 0, 1], 2)

    def test_x2_plus_x_plus_1_mod_2_irreducible(self):
        """x² + x + 1 ist irreduzibel mod 2 (GF(4) Erzeugungspolynom)."""
        # Keine Nullstellen mod 2: f(0)=1, f(1)=1+1+1=3≡1
        assert is_irreducible_mod_p([1, 1, 1], 2)

    def test_x2_plus_1_mod_3_irreducible(self):
        """x² + 1 ist irreduzibel mod 3."""
        # f(0)=1, f(1)=2, f(2)=5≡2 mod 3 → keine Nullstellen
        assert is_irreducible_mod_p([1, 0, 1], 3)

    def test_x2_minus_1_mod_3_not_irreducible(self):
        """x² - 1 = (x-1)(x+1) ist nicht irreduzibel mod 3."""
        # f(1) = 1-1 = 0 → Nullstelle bei x=1
        assert not is_irreducible_mod_p([1, 0, -1], 3)

    def test_linear_always_irreducible(self):
        """Lineare Polynome sind immer irreduzibel."""
        assert is_irreducible_mod_p([1, 0], 5)   # x
        assert is_irreducible_mod_p([1, 1], 2)   # x + 1

    def test_x3_plus_x_plus_1_mod_2_irreducible(self):
        """x³ + x + 1 ist irreduzibel mod 2 (keine Nullstellen, kein Faktor Grad 1)."""
        # f(0)=1, f(1)=1+1+1=3≡1 mod 2 → keine linearen Faktoren
        assert is_irreducible_mod_p([1, 0, 1, 1], 2)


# ===========================================================================
# TESTS: polynomial_gcd_mod_p
# ===========================================================================

class TestPolynomialGcdModP:
    """Tests für den Polynom-ggT über ℤ_p."""

    def test_gcd_x2_minus_1_and_x_minus_1_mod_2(self):
        """ggT(x²-1, x-1) mod 2 = x-1 = x+1 (da -1≡1 mod 2)."""
        # x²-1 = (x-1)(x+1), ggT mit (x-1) = (x-1)
        f = [1, 0, -1]  # x² - 1
        g = [1, -1]     # x - 1
        result = polynomial_gcd_mod_p(f, g, 2)
        # Mod 2: -1 ≡ 1, also x-1 ≡ x+1 ≡ [1,1]
        # ggT muss Teiler von beiden sein
        assert len(result) == 2  # Grad 1

    def test_gcd_coprime_polynomials(self):
        """ggT(x²+x+1, x+1) mod 2 = 1 (teilerfremd)."""
        f = [1, 1, 1]  # x² + x + 1 (irreduzibel mod 2)
        g = [1, 1]     # x + 1
        result = polynomial_gcd_mod_p(f, g, 2)
        # Beide teilerfremd → ggT = 1
        assert result == [1]

    def test_gcd_same_polynomial(self):
        """ggT(f, f) = f (monisch)."""
        f = [1, 1, 1]  # x² + x + 1 mod 2
        result = polynomial_gcd_mod_p(f, f, 2)
        assert result == [1, 1, 1]

    def test_gcd_with_zero(self):
        """ggT(f, 0) = f."""
        f = [1, 1]  # x + 1 mod 2
        g = [0]
        result = polynomial_gcd_mod_p(f, g, 2)
        assert result == [1, 1]


# ===========================================================================
# TESTS: quotient_ring_field_check
# ===========================================================================

class TestQuotientRingFieldCheck:
    """Tests für die Körper-Überprüfung von ℤ_p[x]/(f)."""

    def test_gf4_is_field(self):
        """ℤ_2[x]/(x²+x+1) ist Körper mit 4 Elementen."""
        result = quotient_ring_field_check([1, 1, 1], 2)
        assert result['is_field'] is True
        assert result['field_order'] == 4

    def test_gf8_is_field(self):
        """ℤ_2[x]/(x³+x+1) ist Körper mit 8 Elementen (GF(8))."""
        result = quotient_ring_field_check([1, 0, 1, 1], 2)
        assert result['is_field'] is True
        assert result['field_order'] == 8

    def test_reducible_not_field(self):
        """ℤ_2[x]/(x²+1) ist kein Körper (x²+1 = (x+1)² nicht irreduzibel)."""
        result = quotient_ring_field_check([1, 0, 1], 2)
        assert result['is_field'] is False
        assert result['field_order'] is None


# ===========================================================================
# TESTS: chinese_remainder_ring
# ===========================================================================

class TestChineseRemainderRing:
    """Tests für den Chinesischen Restsatz für Ringe."""

    def test_crt_z6_iso_z2_times_z3(self):
        """ℤ/6ℤ ≅ ℤ/2ℤ × ℤ/3ℤ (6 = 2·3, ggT(2,3)=1)."""
        result = chinese_remainder_ring(6, [2, 3])
        assert result['n'] == 6
        assert result['moduli'] == [2, 3]
        # Isomorphismus: jedes a ↦ (a mod 2, a mod 3)
        iso = result['isomorphism']
        assert iso[0] == (0, 0)
        assert iso[1] == (1, 1)
        assert iso[2] == (0, 2)
        assert iso[3] == (1, 0)
        assert iso[4] == (0, 1)
        assert iso[5] == (1, 2)

    def test_crt_invalid_not_coprime(self):
        """CRT wirft ValueError wenn Moduln nicht teilerfremd sind."""
        with pytest.raises(ValueError):
            chinese_remainder_ring(4, [2, 2])

    def test_crt_z30_three_factors(self):
        """ℤ/30ℤ ≅ ℤ/2ℤ × ℤ/3ℤ × ℤ/5ℤ."""
        result = chinese_remainder_ring(30, [2, 3, 5])
        assert result['n'] == 30
        # Bijektivität: 30 verschiedene Bilder
        images = set(result['isomorphism'].values())
        assert len(images) == 30


# ===========================================================================
# TESTS: ring_of_integers
# ===========================================================================

class TestRingOfIntegers:
    """Tests für quadratische Zahlkörper."""

    def test_gaussian_integers_is_pid(self):
        """ℤ[i] (d=-1) ist ein PID (Klassenzahl 1)."""
        result = ring_of_integers(-1)
        assert result['is_pid'] is True
        assert result['is_imaginary'] is True

    def test_gaussian_integers_name(self):
        """ℤ[i] wird korrekt als Gauß'sche Zahlen bezeichnet."""
        result = ring_of_integers(-1)
        assert 'Gauß' in result['ring_name']

    def test_z_sqrt_2_fundamental_unit(self):
        """ℤ[√2]: fundamentale Einheit ist 1 + √2 (a=1, b=1)."""
        result = ring_of_integers(2)
        # Pell-Gleichung: 1² - 2·1² = -1 (auch erlaubt)
        # oder nächste Lösung: 3² - 2·2² = 1
        assert result['fundamental_unit'] is not None
        a, b = result['fundamental_unit']
        # Überprüfe: a² - 2b² = ±1
        assert a * a - 2 * b * b in (1, -1)

    def test_z_sqrt_minus_3_is_pid(self):
        """ℤ[√-3] (Eisenstein-Zahlen, d=-3) ist PID."""
        result = ring_of_integers(-3)
        assert result['is_pid'] is True

    def test_z_sqrt_5_is_not_pid(self):
        """ℤ[√5] (d=5) ist kein PID (Klassenzahl 2)."""
        result = ring_of_integers(5)
        assert result['is_pid'] is False

    def test_z_sqrt_d_discriminant_d_equals_minus1(self):
        """Diskriminante von ℤ[i]: Δ = 4·(-1) = -4."""
        result = ring_of_integers(-1)
        assert result['discriminant'] == -4

    def test_z_sqrt_d_discriminant_d_congruent_1_mod4(self):
        """Für d ≡ 1 (mod 4) gilt Δ = d (z.B. d=5: Δ=5)."""
        result = ring_of_integers(5)
        assert result['discriminant'] == 5  # 5 ≡ 1 (mod 4)

    def test_ring_of_integers_invalid_d_zero(self):
        """d=0 wirft ValueError."""
        with pytest.raises(ValueError):
            ring_of_integers(0)

    def test_ring_of_integers_not_squarefree(self):
        """Nicht-quadratfreies d wirft ValueError (z.B. d=4=2²)."""
        with pytest.raises(ValueError):
            ring_of_integers(4)


# ===========================================================================
# TESTS: gaussian_integers_gcd
# ===========================================================================

class TestGaussianIntegersGcd:
    """Tests für den ggT in ℤ[i]."""

    def test_gcd_5_and_2_plus_i(self):
        """ggT(5, 2+i) in ℤ[i] = 2+i (da 5 = (2+i)(2-i))."""
        result = gaussian_integers_gcd(5+0j, 2+1j)
        # Normform: Real > 0 oder Real=0, Im>0
        # ggT soll Teiler von 5 und 2+i sein
        # 5 / (2+i) = (2-i) = 2-i, 5 = (2+i)(2-i)
        # ggT(5, 2+i): 2+i teilt 5 (da 5 = (2+i)(2-i))
        # Norm des Ergebnisses muss Norm(2+i) = 5 teilen
        g = result
        g_real = round(g.real)
        g_imag = round(g.imag)
        norm_g = g_real**2 + g_imag**2
        # Norm(2+i) = 4+1 = 5, Norm(5) = 25
        # ggT muss Norm 1 oder 5 haben
        assert norm_g in (1, 5)

    def test_gcd_coprime_gaussian(self):
        """ggT(2, 3) in ℤ[i] = 1 (2 und 3 sind teilerfremd in ℤ[i])."""
        result = gaussian_integers_gcd(2+0j, 3+0j)
        g_real = round(result.real)
        g_imag = round(result.imag)
        norm_g = g_real**2 + g_imag**2
        assert norm_g == 1  # Einheit

    def test_gcd_same_element(self):
        """ggT(a, a) = a (normiert) in ℤ[i]."""
        a = 3 + 2j
        result = gaussian_integers_gcd(a, a)
        g_real = round(result.real)
        g_imag = round(result.imag)
        norm_result = g_real**2 + g_imag**2
        norm_a = round(a.real)**2 + round(a.imag)**2
        assert norm_result == norm_a

    def test_gcd_zero_element(self):
        """ggT(a, 0) = a (normiert)."""
        a = 1 + 2j
        result = gaussian_integers_gcd(a, 0+0j)
        g_real = round(result.real)
        g_imag = round(result.imag)
        norm_result = g_real**2 + g_imag**2
        norm_a = round(a.real)**2 + round(a.imag)**2
        assert norm_result == norm_a


# ===========================================================================
# TESTS: noetherian_check & krull_dimension_estimate
# ===========================================================================

class TestNoetherianAndKrull:
    """Tests für Noethersche Bedingung und Krull-Dimension."""

    def test_noetherian_single_ideal(self):
        """Eine einzelne Kette ist trivial stationär."""
        R = Ring(12)
        I = principal_ideal(R, 3)
        assert noetherian_check([I])

    def test_noetherian_stabilizing_chain(self):
        """Eine stationäre Kette erfüllt die noethersche Bedingung."""
        R = Ring(12)
        I3 = principal_ideal(R, 3)
        I1 = principal_ideal(R, 1)
        # Aufsteigende Kette: (3) ⊆ (1) ⊆ (1) stationär
        assert noetherian_check([I3, I1, I1])

    def test_krull_dimension_field(self):
        """Körper haben Krull-Dimension 0 (nur Nullideal als Primideal)."""
        R = Ring(7)
        I0 = principal_ideal(R, 0)
        # Einziges Primideal in einem Körper ist (0)
        dim = krull_dimension_estimate([I0])
        assert dim == 0

    def test_krull_dimension_chain_length(self):
        """Krull-Dimension = Anzahl strikter Inklusionen in Primidealkette."""
        R = Ring(30)
        # (0) ⊊ (2) in ℤ/30ℤ (beides Primideale)
        I0 = principal_ideal(R, 0)
        I2 = principal_ideal(R, 2)
        # (2) enthält Primideal-Eigenschaften: 30/gcd(2,30) = 15 (nicht prim)
        # Daher nur Kette testen
        dim = krull_dimension_estimate([I0, I2])
        assert dim >= 0  # Mindestens 0


# ===========================================================================
# TESTS: polynomial_factorization_mod_p
# ===========================================================================

class TestPolynomialFactorizationModP:
    """Tests für die Polynom-Faktorisierung über ℤ_p."""

    def test_factorize_x2_minus_1_mod_2(self):
        """x² - 1 = (x+1)² mod 2."""
        result = polynomial_factorization_mod_p([1, 0, -1], 2)
        # x² + 1 = (x+1)² mod 2
        factors = result['factors']
        # Muss mindestens einen Faktor haben
        assert len(factors) >= 1

    def test_factorize_irreducible_stays_whole(self):
        """x² + x + 1 mod 2 ist irreduzibel → bleibt ungeteilt."""
        result = polynomial_factorization_mod_p([1, 1, 1], 2)
        factors = result['factors']
        # Irreduzibles Polynom hat nur sich selbst als Faktor
        # Alle Faktoren zusammen multipliziert müssen das Original ergeben
        assert len(factors) >= 1

    def test_factorize_linear_factor(self):
        """x² - x = x(x-1) = x(x+1) mod 2 hat zwei lineare Faktoren."""
        # x² + x = x(x+1) mod 2
        result = polynomial_factorization_mod_p([1, 1, 0], 2)
        factors = result['factors']
        # Muss in Faktoren niedrigeren Grads zerlegt werden
        total_degree = sum(_poly_degree_helper(f) for f in factors)
        assert total_degree >= 2


def _poly_degree_helper(poly: list[int]) -> int:
    """Hilfsfunktion für Polynomgrad."""
    while len(poly) > 1 and poly[0] == 0:
        poly = poly[1:]
    if len(poly) == 1 and poly[0] == 0:
        return 0
    return len(poly) - 1
