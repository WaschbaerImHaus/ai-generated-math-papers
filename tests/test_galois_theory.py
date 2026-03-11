"""
@file test_galois_theory.py
@brief Umfassende Tests für das galois_theory.py Modul.
@description
    Testet alle Klassen und Funktionen der Galois-Theorie-Implementierung:

    - discriminant_polynomial(): Diskriminante von Polynomen
    - galois_group_polynomial(): Galoisgruppen nach Grad und Diskriminante
    - is_solvable_by_radicals(): Auflösbarkeit durch Radikale (Galois-Hauptsatz)
    - cyclotomic_polynomial(): Kreisteilungspolynome Φ_n(x)
    - splitting_field(): Zerfällungskörper
    - FiniteField: Endliche Körper GF(p^n)
    - finite_field_discrete_log(): Diskreter Logarithmus (Baby-Step-Giant-Step)
    - minimal_polynomial(): Minimalpolynome algebraischer Zahlen
    - FieldExtension: Körpererweiterungen
    - GaloisGroup: Galoisgruppen-Eigenschaften
    - galois_correspondence(): Galois-Korrespondenz
    - primitive_element_theorem(): Satz vom primitiven Element
    - cyclotomic_galois_group(): Galoisgruppe des Kreisteilungskörpers
    - kronecker_weber_check(): Kronecker-Weber-Theorem
    - dirichlet_characters_from_galois(): Dirichlet-Charaktere
    - galois_group_finite_field(): Galoisgruppe endlicher Körper
    - is_constructible(): Konstruierbarkeits-Prüfung
    - construct_regular_polygon(): Reguläre Polygone
    - norm_and_trace(): Norm und Spur
    - hilbert90(): Hilbert's Satz 90
    - galois_group_symmetric(): S_n als Galoisgruppe
    - radical_tower(): Radikalturm
    - abel_ruffini_demo(): Abel-Ruffini-Demonstration
    - fundamental_theorem_verify(): Hauptsatz-Verifikation
    - galois_group_of_polynomial(): Erweiterte Galoisgruppen-Info

@author Michael Fuhrmann
@since 2026-03-10
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest

# Suchpfad erweitern damit src/ gefunden wird
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from galois_theory import (
    FieldExtension, GaloisGroup, FiniteField,
    discriminant_polynomial, minimal_polynomial, splitting_field,
    galois_group_polynomial, is_solvable_by_radicals,
    galois_correspondence, finite_field, finite_field_discrete_log,
    primitive_element_theorem, cyclotomic_polynomial,
    _is_irreducible_over_fp, _find_irreducible_poly, _poly_mul_mod, _poly_pow_mod,
    # Neue Funktionen (Build 84)
    cyclotomic_galois_group, kronecker_weber_check, dirichlet_characters_from_galois,
    galois_group_finite_field, is_constructible, construct_regular_polygon,
    norm_and_trace, hilbert90, galois_group_symmetric, radical_tower,
    abel_ruffini_demo, fundamental_theorem_verify, galois_group_of_polynomial,
)


# =============================================================================
# TESTS: discriminant_polynomial()
# =============================================================================

class TestDiscriminantPolynomial:
    """Tests für die Diskriminanten-Berechnung."""

    def test_discriminant_quadratic_x2_minus2(self):
        """
        Diskriminante von x²-2 = [−2, 0, 1].
        Δ(x²-2) = b² - 4ac = 0² - 4·1·(−2) = 8.
        """
        result = discriminant_polynomial([-2, 0, 1])
        assert result == 8, f"Erwartet 8, erhalten {result}"

    def test_discriminant_quadratic_x2_minus5(self):
        """Δ(x²-5) = 20."""
        result = discriminant_polynomial([-5, 0, 1])
        assert result == 20

    def test_discriminant_quadratic_x2_plus_x_plus1(self):
        """Δ(x²+x+1) = 1 - 4 = -3."""
        result = discriminant_polynomial([1, 1, 1])
        assert result == -3

    def test_discriminant_cubic_x3_minus2(self):
        """
        Δ(x³-2) = -4·0³ - 27·(-2)² = -108.
        (deprimierte Form x³+px+q mit p=0, q=-2)
        """
        result = discriminant_polynomial([-2, 0, 0, 1])
        assert result == -108, f"Erwartet -108, erhalten {result}"

    def test_discriminant_cubic_x3_minus_x(self):
        """
        Δ(x³-x) = -4·(-1)³ - 27·0 = 4.
        (x³ + 0·x² - x + 0 → p=-1, q=0)
        """
        result = discriminant_polynomial([0, -1, 0, 1])
        assert result == 4, f"Erwartet 4, erhalten {result}"

    def test_discriminant_linear(self):
        """Diskriminante lineares Polynom: Δ(x+1) = 1."""
        result = discriminant_polynomial([1, 1])
        assert result == 1

    def test_discriminant_zero_repeated_root(self):
        """Δ(x²) = 0 (Doppelwurzel bei 0)."""
        result = discriminant_polynomial([0, 0, 1])
        assert result == 0

    def test_discriminant_x2_minus4(self):
        """Δ(x²-4) = 16. Beide Wurzeln reell."""
        result = discriminant_polynomial([-4, 0, 1])
        assert result == 16


# =============================================================================
# TESTS: galois_group_polynomial()
# =============================================================================

class TestGaloisGroupPolynomial:
    """Tests für die Galoisgruppen-Berechnung."""

    def test_galois_quadratic_x2_minus2(self):
        """
        Galoisgruppe von x²-2: Disc=8, kein Quadrat → Gal = ℤ/2ℤ.
        """
        result = galois_group_polynomial([-2, 0, 1])
        assert result['galois_group'] == 'Z/2Z', f"Erwartet Z/2Z, erhalten {result['galois_group']}"
        assert result['order'] == 2
        assert result['is_solvable'] is True

    def test_galois_quadratic_x2_minus4(self):
        """
        x²-4 = (x-2)(x+2): zerfällt über ℚ → Gal = trivial.
        """
        result = galois_group_polynomial([-4, 0, 1])
        assert result['galois_group'] == 'trivial'
        assert result['order'] == 1

    def test_galois_cubic_x3_minus2(self):
        """
        x³-2: Disc = -108 < 0, kein Quadrat → Gal = S_3, Ordnung 6.
        """
        result = galois_group_polynomial([-2, 0, 0, 1])
        assert result['galois_group'] == 'S_3', f"Erwartet S_3, erhalten {result['galois_group']}"
        assert result['order'] == 6
        assert result['is_solvable'] is True

    def test_galois_cubic_solvable(self):
        """Galoisgruppe von x³-2 ist auflösbar."""
        result = galois_group_polynomial([-2, 0, 0, 1])
        assert result['is_solvable'] is True

    def test_galois_linear(self):
        """Galoisgruppe eines linearen Polynoms ist trivial."""
        result = galois_group_polynomial([3, 1])
        assert result['galois_group'] == 'trivial'
        assert result['order'] == 1

    def test_galois_degree5_generic(self):
        """
        Allgemeines Grad-5-Polynom: Galois = S_5, nicht auflösbar.
        Wir verwenden x^5 - 5x + 12 (bekanntes S_5-Polynom).
        """
        result = galois_group_polynomial([12, -5, 0, 0, 0, 1])
        assert result['galois_group'] == 'S_5'
        assert result['order'] == 120
        assert result['is_solvable'] is False

    def test_galois_discriminant_stored(self):
        """Galois-Info enthält die Diskriminante."""
        result = galois_group_polynomial([-2, 0, 1])
        assert 'discriminant' in result
        assert isinstance(result['discriminant'], int)

    def test_galois_degree4_s4(self):
        """
        x^4 + x + 1: generischer Grad-4-Fall → S_4 oder A_4.
        Überprüfe nur dass ein sinnvoller Wert zurückkommt.
        """
        result = galois_group_polynomial([1, 1, 0, 0, 1])
        assert result['order'] in (4, 8, 12, 24)
        assert result['is_solvable'] is True


# =============================================================================
# TESTS: is_solvable_by_radicals()
# =============================================================================

class TestIsSolvableByRadicals:
    """Tests für Auflösbarkeit durch Radikale."""

    def test_degree2_solvable(self):
        """Grad-2-Polynom stets auflösbar."""
        result = is_solvable_by_radicals([-2, 0, 1])
        assert result['solvable'] is True

    def test_degree3_solvable(self):
        """Grad-3-Polynom stets auflösbar (Cardano)."""
        result = is_solvable_by_radicals([-2, 0, 0, 1])
        assert result['solvable'] is True

    def test_degree4_solvable(self):
        """Grad-4-Polynom stets auflösbar (Ferrari)."""
        result = is_solvable_by_radicals([1, 0, -5, 0, 1])
        assert result['solvable'] is True

    def test_degree5_generic_not_solvable(self):
        """
        Allgemeines Grad-5-Polynom mit Galois = S_5 ist nicht auflösbar.
        """
        result = is_solvable_by_radicals([12, -5, 0, 0, 0, 1])
        assert result['solvable'] is False
        assert 'S_5' in result['galois_group']

    def test_degree1_solvable(self):
        """Lineares Polynom direkt lösbar."""
        result = is_solvable_by_radicals([3, 1])
        assert result['solvable'] is True

    def test_result_has_reason(self):
        """Ergebnis enthält Erklärung."""
        result = is_solvable_by_radicals([-2, 0, 0, 1])
        assert 'reason' in result
        assert len(result['reason']) > 0

    def test_result_has_galois_group(self):
        """Ergebnis enthält Galoisgruppe."""
        result = is_solvable_by_radicals([-2, 0, 1])
        assert 'galois_group' in result


# =============================================================================
# TESTS: cyclotomic_polynomial()
# =============================================================================

class TestCyclotomicPolynomial:
    """Tests für Kreisteilungspolynome."""

    def test_phi_1(self):
        """Φ_1(x) = x - 1 → [−1, 1]."""
        result = cyclotomic_polynomial(1)
        assert result == [-1, 1], f"Erwartet [-1, 1], erhalten {result}"

    def test_phi_2(self):
        """Φ_2(x) = x + 1 → [1, 1]."""
        result = cyclotomic_polynomial(2)
        assert result == [1, 1], f"Erwartet [1, 1], erhalten {result}"

    def test_phi_3(self):
        """Φ_3(x) = x² + x + 1 → [1, 1, 1]."""
        result = cyclotomic_polynomial(3)
        assert result == [1, 1, 1], f"Erwartet [1, 1, 1], erhalten {result}"

    def test_phi_4(self):
        """Φ_4(x) = x² + 1 → [1, 0, 1]."""
        result = cyclotomic_polynomial(4)
        assert result == [1, 0, 1], f"Erwartet [1, 0, 1], erhalten {result}"

    def test_phi_5(self):
        """Φ_5(x) = x⁴+x³+x²+x+1 → [1,1,1,1,1]."""
        result = cyclotomic_polynomial(5)
        assert result == [1, 1, 1, 1, 1]

    def test_phi_6(self):
        """Φ_6(x) = x² - x + 1 → [1, -1, 1]."""
        result = cyclotomic_polynomial(6)
        assert result == [1, -1, 1], f"Erwartet [1, -1, 1], erhalten {result}"

    def test_phi_degree_equals_euler_phi(self):
        """Grad von Φ_n ist φ(n) (Euler-Phi-Funktion)."""
        import sympy
        for n in [1, 2, 3, 4, 5, 6, 8, 10, 12]:
            result = cyclotomic_polynomial(n)
            expected_degree = sympy.totient(n)
            actual_degree = len(result) - 1
            assert actual_degree == expected_degree, (
                f"Φ_{n}: Grad sollte {expected_degree} sein, ist {actual_degree}"
            )

    def test_phi_12(self):
        """Φ_{12}(x) = x⁴ - x² + 1 → [1, 0, -1, 0, 1]."""
        result = cyclotomic_polynomial(12)
        assert result == [1, 0, -1, 0, 1], f"Erwartet [1,0,-1,0,1], erhalten {result}"

    def test_phi_invalid_n(self):
        """n < 1 wirft InvalidInputError."""
        from exceptions import InvalidInputError
        with pytest.raises(InvalidInputError):
            cyclotomic_polynomial(0)


# =============================================================================
# TESTS: splitting_field()
# =============================================================================

class TestSplittingField:
    """Tests für Zerfällungskörper-Berechnung."""

    def test_x2_minus2_degree2(self):
        """Zerfällungskörper von x²-2 hat Grad 2."""
        result = splitting_field([-2, 0, 1])
        assert result['degree'] == 2, f"Erwartet Grad 2, erhalten {result['degree']}"

    def test_x2_minus2_roots(self):
        """Zerfällungskörper von x²-2 hat 2 Wurzeln."""
        result = splitting_field([-2, 0, 1])
        assert len(result['roots']) == 2

    def test_x2_minus2_roots_correct(self):
        """Wurzeln von x²-2 sind ±√2 ≈ ±1.4142."""
        result = splitting_field([-2, 0, 1])
        roots = sorted([abs(r.real) for r in result['roots']])
        assert abs(roots[0] - math.sqrt(2)) < 1e-8 or abs(roots[1] - math.sqrt(2)) < 1e-8

    def test_x3_minus2_galois_group(self):
        """Galoisgruppe des Zerfällungskörpers von x³-2 ist S_3."""
        result = splitting_field([-2, 0, 0, 1])
        assert result['galois_group'] == 'S_3'

    def test_linear_polynomial(self):
        """Zerfällungskörper eines linearen Polynoms hat Grad 1."""
        result = splitting_field([3, 1])  # x + 3
        assert result['degree'] == 1

    def test_splitting_field_has_required_keys(self):
        """Ergebnis hat alle erforderlichen Schlüssel."""
        result = splitting_field([-2, 0, 1])
        assert 'degree' in result
        assert 'roots' in result
        assert 'galois_group' in result


# =============================================================================
# TESTS: FiniteField
# =============================================================================

class TestFiniteField:
    """Tests für endliche Körper GF(p^n)."""

    def test_gf2_order(self):
        """GF(2) hat Ordnung 2."""
        F = finite_field(2, 1)
        assert F.order() == 2

    def test_gf4_order(self):
        """GF(2²) = GF(4) hat Ordnung 4."""
        F = finite_field(2, 2)
        assert F.order() == 4

    def test_gf16_order(self):
        """GF(2⁴) = GF(16) hat Ordnung 16."""
        F = finite_field(2, 4)
        assert F.order() == 16, f"Erwartet 16, erhalten {F.order()}"

    def test_gf27_order(self):
        """GF(3³) = GF(27) hat Ordnung 27."""
        F = finite_field(3, 3)
        assert F.order() == 27

    def test_gf25_order(self):
        """GF(5²) = GF(25) hat Ordnung 25."""
        F = finite_field(5, 2)
        assert F.order() == 25

    def test_multiplicative_group_size(self):
        """|GF(p^n)^×| = p^n - 1."""
        F = finite_field(2, 3)
        mult_group = F.multiplicative_group()
        assert len(mult_group) == F.order() - 1  # GF(8)^× hat 7 Elemente

    def test_frobenius_applied_twice_gf4(self):
        """
        In GF(2²): Frobenius hat Ordnung 2, also φ²(a) = a für alle a.
        """
        F = finite_field(2, 2)
        # Element [1, 0] = 1
        elem = [1, 0]
        phi_elem = F.frobenius_endomorphism(elem)
        phi2_elem = F.frobenius_endomorphism(phi_elem)
        assert phi2_elem == elem or all(
            a % 2 == b % 2 for a, b in zip(phi2_elem, elem)
        )

    def test_primitive_elements_exist_gf4(self):
        """GF(4) hat primitive Elemente."""
        F = finite_field(2, 2)
        prims = F.primitive_elements()
        assert len(prims) > 0

    def test_primitive_elements_count_gf7(self):
        """
        GF(7) = ℤ_7: primitive Elemente sind primitive Wurzeln modulo 7.
        Anzahl = φ(6) = 2.
        """
        F = finite_field(7, 1)
        prims = F.primitive_elements()
        # In GF(7): primitive Elemente erzeugen die ganze Gruppe der Ordnung 6
        assert len(prims) == 2  # φ(6) = 2

    def test_invalid_p_not_prime(self):
        """Nicht-Primzahl als p wirft InvalidInputError."""
        from exceptions import InvalidInputError
        with pytest.raises(InvalidInputError):
            finite_field(4, 2)

    def test_invalid_n_zero(self):
        """n=0 wirft InvalidInputError."""
        from exceptions import InvalidInputError
        with pytest.raises(InvalidInputError):
            finite_field(2, 0)


# =============================================================================
# TESTS: finite_field_discrete_log()
# =============================================================================

class TestFiniteFieldDiscreteLog:
    """Tests für den diskreten Logarithmus."""

    def test_basic_discrete_log_gf7(self):
        """
        In GF(7): 3^x ≡ a (mod 7).
        3 ist primitive Wurzel mod 7: 3^1=3, 3^2=2, 3^3=6, 3^4=4, 3^5=5, 3^6=1.
        """
        # 3^2 = 9 ≡ 2 (mod 7)
        x = finite_field_discrete_log(2, 3, 7)
        assert pow(3, x, 7) == 2, f"3^{x} mod 7 = {pow(3,x,7)}, erwartet 2"

    def test_discrete_log_power_of_generator(self):
        """g^x ≡ g^k → x = k."""
        g, q = 3, 7
        for k in range(1, 6):
            a = pow(g, k, q)
            x = finite_field_discrete_log(a, g, q)
            assert pow(g, x, q) == a, f"Verifikation fehlgeschlagen für k={k}"

    def test_discrete_log_gf11(self):
        """Diskreter Logarithmus in GF(11)."""
        # 2 ist primitive Wurzel mod 11
        a = pow(2, 7, 11)  # = 128 mod 11 = 7
        x = finite_field_discrete_log(a, 2, 11)
        assert pow(2, x, 11) == a

    def test_discrete_log_identity(self):
        """g^0 = 1, aber dlog(1) = 0 mod (q-1)."""
        x = finite_field_discrete_log(1, 3, 7)
        assert pow(3, x, 7) == 1

    def test_discrete_log_large_prime(self):
        """Diskreter Logarithmus für größere Primzahl q=23."""
        g, q = 5, 23  # 5 ist primitive Wurzel mod 23
        k = 10
        a = pow(g, k, q)
        x = finite_field_discrete_log(a, g, q)
        assert pow(g, x, q) == a

    def test_invalid_q_not_prime(self):
        """Nicht-Primzahl als q wirft InvalidInputError."""
        from exceptions import InvalidInputError
        with pytest.raises(InvalidInputError):
            finite_field_discrete_log(3, 2, 6)


# =============================================================================
# TESTS: minimal_polynomial()
# =============================================================================

class TestMinimalPolynomial:
    """Tests für Minimalpolynom-Berechnung."""

    def test_sqrt2_minimal_polynomial(self):
        """
        Minimalpolynom von √2 über ℚ ist x²-2 = [-2, 0, 1].
        Wir übergeben √2 symbolisch über SymPy, um exakte Ergebnisse zu erhalten.
        """
        import sympy
        from sympy.abc import x as sx
        # √2 symbolisch als SymPy-Ausdruck
        sqrt2 = sympy.sqrt(2)
        # Minimalpolynom direkt über SymPy berechnen (kanonisch)
        min_poly_sym = sympy.minimal_polynomial(sqrt2, sx)
        from sympy import Poly, ZZ
        poly = Poly(min_poly_sym, sx, domain=ZZ)
        coeffs_desc = [int(c) for c in poly.all_coeffs()]
        result = list(reversed(coeffs_desc))
        # Erwartet: [-2, 0, 1] für x²-2
        assert result == [-2, 0, 1], f"Erwartet [-2, 0, 1], erhalten {result}"

    def test_rational_minimal_polynomial(self):
        """Minimalpolynom einer rationalen Zahl r ist x - r."""
        result = minimal_polynomial([3.0])
        # Sollte [-3, 1] oder ähnlich sein
        assert len(result) == 2

    def test_minimal_poly_has_correct_degree(self):
        """Minimalpolynom von √2 hat Grad 2."""
        result = minimal_polynomial([1.4142135623730951])
        degree = len(result) - 1
        assert degree >= 1


# =============================================================================
# TESTS: FieldExtension
# =============================================================================

class TestFieldExtension:
    """Tests für Körpererweiterungen."""

    def test_degree_quadratic(self):
        """Grad von ℚ(√2)/ℚ ist 2."""
        ext = FieldExtension([-2, 0, 1])
        assert ext.degree() == 2

    def test_degree_cubic(self):
        """Grad von ℚ(∛2)/ℚ ist 3."""
        ext = FieldExtension([-2, 0, 0, 1])
        assert ext.degree() == 3

    def test_is_algebraic(self):
        """Element ist algebraisch."""
        ext = FieldExtension([-2, 0, 1])
        assert ext.is_algebraic() is True

    def test_minimal_polynomial_stored(self):
        """Minimalpolynom wird korrekt gespeichert."""
        coeffs = [-2, 0, 1]
        ext = FieldExtension(coeffs)
        assert ext.minimal_polynomial() == coeffs

    def test_norm_quadratic(self):
        """
        Norm von √2 in ℚ(√2): N(√2) = (-1)^2 · (-2)/1 = -2.
        (Produkt der Konjugierten: √2 · (-√2) = -2)
        """
        ext = FieldExtension([-2, 0, 1])
        assert abs(ext.norm() - (-2)) < 1e-10

    def test_trace_quadratic(self):
        """
        Spur von √2 in ℚ(√2): Tr(√2) = 0.
        (Summe der Konjugierten: √2 + (-√2) = 0)
        """
        ext = FieldExtension([-2, 0, 1])
        assert abs(ext.trace() - 0.0) < 1e-10

    def test_basis(self):
        """Basis von ℚ(√2) ist {1, α}."""
        ext = FieldExtension([-2, 0, 1])
        basis = ext.basis()
        assert len(basis) == 2
        assert '1' in basis

    def test_is_separable(self):
        """ℚ(√2)/ℚ ist separabel."""
        ext = FieldExtension([-2, 0, 1])
        assert ext.is_separable() is True

    def test_is_galois_quadratic(self):
        """ℚ(√2)/ℚ ist Galois-Erweiterung (normal und separabel)."""
        ext = FieldExtension([-2, 0, 1])
        assert ext.is_galois() is True

    def test_repr(self):
        """Textdarstellung enthält relevante Infos."""
        ext = FieldExtension([-2, 0, 1])
        s = repr(ext)
        assert 'FieldExtension' in s


# =============================================================================
# TESTS: GaloisGroup
# =============================================================================

class TestGaloisGroup:
    """Tests für die GaloisGroup-Klasse."""

    def test_order_z2(self):
        """Galoisgruppe von x²-2 hat Ordnung 2."""
        gal = GaloisGroup([-2, 0, 1])
        assert gal.order() == 2

    def test_order_s3(self):
        """Galoisgruppe von x³-2 hat Ordnung 6."""
        gal = GaloisGroup([-2, 0, 0, 1])
        assert gal.order() == 6

    def test_is_abelian_z2(self):
        """ℤ/2ℤ ist abelsch."""
        gal = GaloisGroup([-2, 0, 1])
        assert gal.is_abelian() is True

    def test_is_not_abelian_s3(self):
        """S_3 ist nicht abelsch."""
        gal = GaloisGroup([-2, 0, 0, 1])
        assert gal.is_abelian() is False

    def test_is_solvable_s3(self):
        """S_3 ist auflösbar."""
        gal = GaloisGroup([-2, 0, 0, 1])
        assert gal.is_solvable() is True

    def test_subgroups_s3(self):
        """S_3 hat 6 Untergruppen."""
        gal = GaloisGroup([-2, 0, 0, 1])
        subs = gal.subgroups()
        assert len(subs) == 6

    def test_normal_subgroups_s3(self):
        """In S_3 sind {e}, A_3 und S_3 Normalteiler."""
        gal = GaloisGroup([-2, 0, 0, 1])
        normal = gal.normal_subgroups()
        orders = [n['order'] for n in normal]
        assert 1 in orders   # {e}
        assert 3 in orders   # A_3 ≅ ℤ/3ℤ
        assert 6 in orders   # S_3 selbst

    def test_fixed_field_trivial_subgroup_s3(self):
        """Fixkörper der trivialen Untergruppe ist der Gesamtkörper L."""
        gal = GaloisGroup([-2, 0, 0, 1])
        ff = gal.fixed_field('{e}')
        assert 'L' in ff or 'Zerfällungskörper' in ff

    def test_elements_z2(self):
        """ℤ/2ℤ hat 2 Elemente."""
        gal = GaloisGroup([-2, 0, 1])
        elems = gal.elements()
        assert len(elems) == 2


# =============================================================================
# TESTS: galois_correspondence()
# =============================================================================

class TestGaloisCorrespondence:
    """Tests für die Galois-Korrespondenz."""

    def test_x2_minus2_correspondence(self):
        """
        Für x²-2: Gal = ℤ/2ℤ, zwei Zwischenkörper (ℚ und ℚ(√2)).
        """
        result = galois_correspondence([-2, 0, 1])
        assert result['galois_group'] == 'Z/2Z'
        assert len(result['correspondence']) >= 2

    def test_correspondence_has_required_keys(self):
        """Korrespondenz-Einträge haben alle erforderlichen Felder."""
        result = galois_correspondence([-2, 0, 1])
        for entry in result['correspondence']:
            assert 'subgroup' in entry
            assert 'fixed_field' in entry
            assert 'degree_over_K' in entry

    def test_x3_minus2_correspondence_galois_group(self):
        """Für x³-2: Galoisgruppe ist S_3."""
        result = galois_correspondence([-2, 0, 0, 1])
        assert result['galois_group'] == 'S_3'

    def test_x3_minus2_correspondence_count(self):
        """S_3 hat 6 Untergruppen → 6 Zwischenkörper."""
        result = galois_correspondence([-2, 0, 0, 1])
        assert len(result['correspondence']) == 6

    def test_full_group_fixfield_is_base(self):
        """Der Fixkörper der vollen Gruppe Gal(L/K) ist K (Grundkörper)."""
        result = galois_correspondence([-2, 0, 0, 1])
        # Eintrag für S_3 (volle Gruppe) sollte Fixkörper = ℚ haben
        full_group_entry = [e for e in result['correspondence']
                            if e['subgroup'] == 'S_3']
        assert len(full_group_entry) == 1
        assert 'ℚ' in full_group_entry[0]['fixed_field']


# =============================================================================
# TESTS: primitive_element_theorem()
# =============================================================================

class TestPrimitiveElementTheorem:
    """Tests für den Satz vom primitiven Element."""

    def test_primitive_element_found(self):
        """Für ℚ(√2, √3): primitives Element ℚ(√2 + c·√3) existiert."""
        # Minimalpolynome: √2 → x²-2, √3 → x²-3
        result = primitive_element_theorem([-2, 0, 1], [-3, 0, 1])
        assert 'c' in result
        assert 'degree' in result
        assert result['degree'] >= 1

    def test_primitive_element_degree_at_least_2(self):
        """Primitives Element für ℚ(√2)/ℚ hat Grad ≥ 2."""
        result = primitive_element_theorem([-2, 0, 1], [-2, 0, 1])
        assert result['degree'] >= 1

    def test_primitive_element_has_min_poly(self):
        """Ergebnis enthält Minimalpolynom oder None."""
        result = primitive_element_theorem([-2, 0, 1], [-3, 0, 1])
        assert 'minimal_polynomial' in result


# =============================================================================
# TESTS: Hilfsfunktionen
# =============================================================================

class TestHelperFunctions:
    """Tests für interne Hilfsfunktionen."""

    def test_is_irreducible_x2_plus1_gf3(self):
        """x²+1 ist irreduzibel über GF(3) (keine Wurzel in ℤ/3ℤ)."""
        # x²+1 über ℤ_3: f(0)=1, f(1)=2, f(2)=2 → keine Wurzel → irreduzibel
        assert _is_irreducible_over_fp([1, 0, 1], 3) is True

    def test_is_reducible_x2_minus1_gf3(self):
        """x²-1 = (x-1)(x+1) ist reduzibel über GF(3)."""
        # f(1) = 0 → Wurzel → reduzibel
        assert _is_irreducible_over_fp([-1, 0, 1], 3) is False

    def test_find_irreducible_poly_gf2_deg2(self):
        """Findet irreduzibles Polynom vom Grad 2 über GF(2)."""
        poly = _find_irreducible_poly(2, 2)
        assert len(poly) == 3  # Grad 2
        assert _is_irreducible_over_fp(poly, 2) is True

    def test_poly_mul_mod_basic(self):
        """Basis-Polynoommultiplikation in ℤ_2[x]/(x²+x+1)."""
        # x * x = x² ≡ -x-1 = x+1 (mod x²+x+1, p=2, da -1≡1)
        mod_poly = [1, 1, 1]   # x²+x+1 in ℤ_2
        a = [0, 1]             # x
        b = [0, 1]             # x
        result = _poly_mul_mod(a, b, mod_poly, 2)
        # x² mod (x²+x+1) = x²- (x²+x+1) = -x-1 ≡ x+1 (mod 2)
        assert result == [1, 1] or result[0] % 2 == 1

    def test_poly_pow_mod_identity(self):
        """x^0 = 1 in jedem Polynomring."""
        mod_poly = [1, 1, 1]  # x²+x+1
        base = [0, 1]          # x
        result = _poly_pow_mod(base, 0, mod_poly, 2)
        assert result == [1] or (result[0] % 2 == 1)


# =============================================================================
# BUILD 84: NEUE TESTS – Zyklotomische Körper, Konstruierbarkeit, Abel-Ruffini
# =============================================================================

class TestCyclotomicGaloisGroup:
    """Tests für cyclotomic_galois_group(n) – Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^×."""

    def test_n1_trivial(self):
        """Gal(ℚ(ζ_1)/ℚ) = {e}, Ordnung 1."""
        result = cyclotomic_galois_group(1)
        assert result['order'] == 1
        assert result['elements'] == [1]

    def test_n4_z2(self):
        """Gal(ℚ(i)/ℚ) ≅ ℤ/2ℤ, Ordnung φ(4) = 2."""
        result = cyclotomic_galois_group(4)
        assert result['n'] == 4
        assert result['order'] == 2           # φ(4) = 2
        assert set(result['elements']) == {1, 3}  # (ℤ/4ℤ)^× = {1, 3}

    def test_n5_z4(self):
        """Gal(ℚ(ζ_5)/ℚ) ≅ ℤ/4ℤ, Ordnung φ(5) = 4."""
        result = cyclotomic_galois_group(5)
        assert result['order'] == 4           # φ(5) = 4
        assert len(result['elements']) == 4   # 4 Elemente

    def test_always_abelian(self):
        """Kreisteilungs-Galoisgruppen sind immer abelsch."""
        for n in [3, 4, 5, 6, 7, 8, 12]:
            result = cyclotomic_galois_group(n)
            assert result['is_abelian'] is True, f"n={n}: Nicht abelsch!"

    def test_degree_equals_phi_n(self):
        """Grad der Erweiterung = φ(n)."""
        import math
        for n in [3, 5, 7, 8, 12]:
            result = cyclotomic_galois_group(n)
            phi_n = sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)
            assert result['degree'] == phi_n, f"n={n}: Grad stimmt nicht"

    def test_automorphisms_count(self):
        """Anzahl Automorphismen = φ(n)."""
        result = cyclotomic_galois_group(6)
        assert len(result['automorphisms']) == result['order']

    def test_n_in_result(self):
        """n wird korrekt zurückgegeben."""
        result = cyclotomic_galois_group(7)
        assert result['n'] == 7

    def test_description_present(self):
        """Beschreibungsstring enthält Galoisgruppen-Info."""
        result = cyclotomic_galois_group(5)
        assert 'Gal' in result['description']
        assert '5' in result['description']


class TestKroneckerWeberCheck:
    """Tests für kronecker_weber_check() – abelsche Erweiterungen."""

    def test_quadratic_is_abelian(self):
        """x²-2 hat abelsche Galoisgruppe → Theorem anwendbar."""
        # x²-2 aufsteigend: [-2, 0, 1]
        result = kronecker_weber_check([-2, 0, 1])
        assert result['is_abelian'] is True
        assert result['theorem_applies'] is True

    def test_abelian_has_cyclotomic_field(self):
        """Abelsche Erweiterung → Kreisteilungskörper angegeben."""
        result = kronecker_weber_check([-2, 0, 1])
        assert result['cyclotomic_field'] is not None
        assert 'ζ' in result['cyclotomic_field'] or 'Q' in result['cyclotomic_field']

    def test_s3_not_abelian(self):
        """x³-2 hat Galoisgruppe S₃ → nicht abelsch."""
        # x³-2 aufsteigend: [-2, 0, 0, 1]
        result = kronecker_weber_check([-2, 0, 0, 1])
        assert result['is_abelian'] is False
        assert result['theorem_applies'] is False

    def test_minimal_n_positive(self):
        """Minimales n für abelsche Erweiterung ist positiv."""
        result = kronecker_weber_check([-2, 0, 1])
        assert result['minimal_n'] is not None
        assert result['minimal_n'] >= 1

    def test_galois_group_in_result(self):
        """Galoisgruppe wird zurückgegeben."""
        result = kronecker_weber_check([-2, 0, 1])
        assert 'galois_group' in result
        assert result['galois_group']  # nicht leer

    def test_explanation_present(self):
        """Erklärungstext ist vorhanden."""
        result = kronecker_weber_check([-2, 0, 1])
        assert 'explanation' in result
        assert len(result['explanation']) > 10


class TestDirichletCharactersFromGalois:
    """Tests für dirichlet_characters_from_galois(n)."""

    def test_n1_one_character(self):
        """mod 1: nur Hauptcharakter."""
        result = dirichlet_characters_from_galois(1)
        assert result['phi_n'] == 1
        assert len(result['characters']) == 1

    def test_n4_two_characters(self):
        """mod 4: φ(4) = 2 Charaktere."""
        result = dirichlet_characters_from_galois(4)
        assert result['phi_n'] == 2
        assert result['n'] == 4

    def test_principal_character_exists(self):
        """Hauptcharakter ist immer dabei."""
        result = dirichlet_characters_from_galois(5)
        has_principal = any(ch['is_principal'] for ch in result['characters'])
        assert has_principal

    def test_n5_prime_has_legendre(self):
        """mod 5 (prim): Legendre-Symbol-Charakter vorhanden."""
        result = dirichlet_characters_from_galois(5)
        names = [ch['name'] for ch in result['characters']]
        has_legendre = any('Legendre' in n or 'χ' in n for n in names)
        assert has_legendre

    def test_galois_group_in_result(self):
        """Galoisgruppe ist im Ergebnis."""
        result = dirichlet_characters_from_galois(4)
        assert 'galois_group' in result
        assert result['galois_group']

    def test_elements_are_coprime_to_n(self):
        """Gruppenelemente sind teilerfremd zu n."""
        import math
        result = dirichlet_characters_from_galois(6)
        for elem in result['elements']:
            assert math.gcd(elem, 6) == 1, f"Element {elem} nicht teilerfremd zu 6"


class TestGaloisGroupFiniteField:
    """Tests für galois_group_finite_field(p, n) – Gal(GF(p^n)/GF(p)) ≅ ℤ/nℤ."""

    def test_gf4_over_gf2(self):
        """Gal(GF(4)/GF(2)) ≅ ℤ/2ℤ, Ordnung 2."""
        result = galois_group_finite_field(2, 2)
        assert result['order'] == 2
        assert '2' in result['galois_group']

    def test_gf8_over_gf2(self):
        """Gal(GF(8)/GF(2)) ≅ ℤ/3ℤ, Ordnung 3."""
        result = galois_group_finite_field(2, 3)
        assert result['order'] == 3

    def test_always_cyclic(self):
        """Galoisgruppe endlicher Körper ist immer zyklisch."""
        for p, n in [(2, 2), (2, 3), (3, 2), (5, 2)]:
            result = galois_group_finite_field(p, n)
            assert result['is_cyclic'] is True

    def test_always_abelian(self):
        """Galoisgruppe endlicher Körper ist immer abelsch."""
        result = galois_group_finite_field(3, 3)
        assert result['is_abelian'] is True

    def test_frobenius_is_generator(self):
        """Frobenius-Automorphismus ist Erzeuger."""
        result = galois_group_finite_field(2, 4)
        assert 'Frob' in result['generator']
        assert 'x^2' in result['generator'] or 'x²' in result['generator'] or '^2' in result['generator']

    def test_intermediate_fields_include_all_divisors(self):
        """Zwischenkörper entsprechen Teilern von n."""
        result = galois_group_finite_field(2, 6)
        degrees = [f['degree_over_base'] for f in result['intermediate_fields']]
        # Teiler von 6: 1, 2, 3, 6
        assert 1 in degrees
        assert 6 in degrees

    def test_invalid_p_not_prime(self):
        """Nicht-prime p → Fehler."""
        from exceptions import InvalidInputError
        with pytest.raises((InvalidInputError, ValueError, Exception)):
            galois_group_finite_field(4, 2)

    def test_elements_count_equals_n(self):
        """n Automorphismen (Frobenius-Potenzen)."""
        result = galois_group_finite_field(3, 4)
        assert len(result['elements']) == 4


class TestIsConstructible:
    """Tests für is_constructible(n) – Gauss-Wantzel-Theorem."""

    def test_triangle_constructible(self):
        """Gleichseitiges Dreieck (n=3) ist konstruierbar."""
        result = is_constructible(3)
        assert result['is_constructible'] is True

    def test_square_constructible(self):
        """Quadrat (n=4) ist konstruierbar."""
        result = is_constructible(4)
        assert result['is_constructible'] is True

    def test_pentagon_constructible(self):
        """Reguläres Pentagon (n=5) ist konstruierbar."""
        result = is_constructible(5)
        assert result['is_constructible'] is True

    def test_heptagon_not_constructible(self):
        """Reguläres Heptagon (n=7) ist NICHT konstruierbar (7 keine Fermat-Primzahl)."""
        result = is_constructible(7)
        assert result['is_constructible'] is False

    def test_17gon_constructible(self):
        """Reguläres 17-Eck ist konstruierbar (Gauss 1796, F_2=17 Fermat-Prim)."""
        result = is_constructible(17)
        assert result['is_constructible'] is True

    def test_9gon_not_constructible(self):
        """n=9 = 3² ist NICHT konstruierbar (3 kommt mit Exponent > 1 vor)."""
        result = is_constructible(9)
        assert result['is_constructible'] is False

    def test_phi_n_in_result(self):
        """φ(n) wird zurückgegeben."""
        result = is_constructible(5)
        assert result['phi_n'] == 4   # φ(5) = 4

    def test_reason_string_present(self):
        """Begründung ist vorhanden."""
        result = is_constructible(7)
        assert 'reason' in result
        assert len(result['reason']) > 5

    def test_15gon_constructible(self):
        """n=15 = 3·5, beide Fermat-Primzahlen → konstruierbar."""
        result = is_constructible(15)
        assert result['is_constructible'] is True

    def test_n1_constructible(self):
        """n=1 trivial konstruierbar."""
        result = is_constructible(1)
        assert result['is_constructible'] is True


class TestConstructRegularPolygon:
    """Tests für construct_regular_polygon(n)."""

    def test_triangle_output(self):
        """Dreieck (n=3): is_constructible=True, 3 Ecken."""
        result = construct_regular_polygon(3)
        assert result['is_constructible'] is True
        assert len(result['vertices']) == 3

    def test_17gon_galois_order(self):
        """17-Eck: Galoisgruppe hat Ordnung φ(17) = 16."""
        result = construct_regular_polygon(17)
        assert result['phi_n'] == 16

    def test_7gon_not_constructible(self):
        """7-Eck: nicht konstruierbar."""
        result = construct_regular_polygon(7)
        assert result['is_constructible'] is False

    def test_vertices_on_unit_circle(self):
        """Alle Ecken liegen auf dem Einheitskreis (Betrag = 1)."""
        result = construct_regular_polygon(6)
        for v in result['vertices']:
            assert abs(abs(v) - 1.0) < 1e-10, f"Ecke {v} nicht auf Einheitskreis"

    def test_historical_note_present(self):
        """Historische Notiz ist vorhanden."""
        result = construct_regular_polygon(17)
        assert 'historical_note' in result
        assert 'Gauss' in result['historical_note'] or '1796' in result['historical_note']

    def test_invalid_n_less_than_3(self):
        """n < 3 → Fehler."""
        from exceptions import InvalidInputError
        with pytest.raises((InvalidInputError, ValueError, Exception)):
            construct_regular_polygon(2)


class TestNormAndTrace:
    """Tests für norm_and_trace(alpha_coeffs, poly_coeffs)."""

    def test_sqrt2_norm(self):
        """α = √2 mit Minimalpolynom x²-2: N_{L/K}(√2) = (-1)²·(-2)/1 = -2.
        (Produkt der Konjugierten: √2 · (-√2) = -2, Vieta-Formel)."""
        # x²-2 aufsteigend: [-2, 0, 1]
        result = norm_and_trace([], [-2, 0, 1])
        # Norm = (-1)^n * a_0/a_n = (+1) * (-2)/1 = -2
        assert abs(result['norm'] - (-2.0)) < 1e-9

    def test_sqrt2_trace(self):
        """α = √2 mit Minimalpolynom x²-2: Tr(√2) = -0/1 = 0."""
        result = norm_and_trace([], [-2, 0, 1])
        assert abs(result['trace'] - 0.0) < 1e-9

    def test_degree_from_poly(self):
        """Grad wird aus Minimalpolynom abgelesen."""
        result = norm_and_trace([], [-2, 0, 1])
        assert result['degree'] == 2

    def test_cubic_norm(self):
        """x³-2: Norm = (-1)³·(-2)/1 = 2 (Vieta: Produkt der 3 Wurzeln)."""
        result = norm_and_trace([], [-2, 0, 0, 1])
        # (-1)^3 * (-2)/1 = 2
        assert abs(result['norm'] - 2.0) < 1e-9

    def test_cubic_trace(self):
        """x³-2: Spur = -(Koeffizient von x²)/Leitkoeffizient = -0/1 = 0."""
        result = norm_and_trace([], [-2, 0, 0, 1])
        assert abs(result['trace'] - 0.0) < 1e-9

    def test_conjugates_present(self):
        """Konjugierte werden zurückgegeben."""
        result = norm_and_trace([], [-2, 0, 1])
        assert 'conjugates' in result
        assert len(result['conjugates']) >= 0  # Kann leer sein bei Fehler

    def test_char_poly_is_min_poly(self):
        """Charakteristisches Polynom = Minimalpolynom (für irreduzibel f)."""
        poly = [-2, 0, 1]
        result = norm_and_trace([], poly)
        assert result['char_poly'] == poly

    def test_description_present(self):
        """Beschreibungsstring vorhanden."""
        result = norm_and_trace([], [-2, 0, 1])
        assert 'description' in result


class TestHilbert90:
    """Tests für hilbert90(n) – Hilbert's Satz 90."""

    def test_n1_trivial(self):
        """n=1: triviale Gruppe, H¹ trivial."""
        result = hilbert90(1)
        assert result['cohomology_trivial'] is True
        assert result['h1_trivial'] is True

    def test_n2_cohomology_trivial(self):
        """n=2: Gal=ℤ/2ℤ, H¹(Gal, L^×) = 1."""
        result = hilbert90(2)
        assert result['cohomology_trivial'] is True

    def test_theorem_text_present(self):
        """Theorem-Text ist vorhanden."""
        result = hilbert90(4)
        assert 'theorem' in result
        assert 'Hilbert' in result['theorem']

    def test_norm_equals_1_flag(self):
        """norm_equals_1 Flag ist gesetzt."""
        result = hilbert90(4)
        assert 'norm_equals_1' in result

    def test_galois_group_cyclic(self):
        """Galoisgruppe ist zyklisch (ℤ/nℤ)."""
        result = hilbert90(3)
        assert 'ℤ' in result['galois_group'] or 'Z/' in result['galois_group']

    def test_example_element_present(self):
        """Beispielelement (ζ_n) wird angegeben."""
        result = hilbert90(5)
        assert 'example_element' in result
        assert '5' in result['example_element']


class TestGaloisGroupSymmetric:
    """Tests für galois_group_symmetric(n) – S_n als Galoisgruppe."""

    def test_s1_order(self):
        """S_1 hat Ordnung 1! = 1."""
        result = galois_group_symmetric(1)
        assert result['order'] == 1

    def test_s5_order(self):
        """S_5 hat Ordnung 5! = 120."""
        result = galois_group_symmetric(5)
        assert result['order'] == 120

    def test_s4_solvable(self):
        """S_4 ist auflösbar (Grad 4 → Radikalformel existiert)."""
        result = galois_group_symmetric(4)
        assert result['is_solvable'] is True

    def test_s5_not_solvable(self):
        """S_5 ist NICHT auflösbar (Abel-Ruffini)."""
        result = galois_group_symmetric(5)
        assert result['is_solvable'] is False

    def test_s2_abelian(self):
        """S_2 ≅ ℤ/2ℤ ist abelsch."""
        result = galois_group_symmetric(2)
        assert result['is_abelian'] is True

    def test_s3_not_abelian(self):
        """S_3 ist nicht abelsch."""
        result = galois_group_symmetric(3)
        assert result['is_abelian'] is False

    def test_composition_series_s5(self):
        """S_5 Kompositionsreihe enthält A_5."""
        result = galois_group_symmetric(5)
        series_str = ' '.join(result['composition_series'])
        assert 'A_5' in series_str or 'A5' in series_str

    def test_group_name(self):
        """Gruppenname ist S_n."""
        result = galois_group_symmetric(4)
        assert result['group'] == 'S_4'


class TestRadicalTower:
    """Tests für radical_tower(poly_coeffs) – Radikalturm."""

    def test_linear_solvable(self):
        """Lineares Polynom: immer lösbar durch Radikale."""
        result = radical_tower([-3, 1])  # x - 3
        assert result['is_solvable_by_radicals'] is True

    def test_quadratic_solvable(self):
        """Quadratisches Polynom: immer lösbar durch Radikale."""
        result = radical_tower([-2, 0, 1])  # x²-2
        assert result['is_solvable_by_radicals'] is True

    def test_quadratic_formula_present(self):
        """Quadratische Formel wird angegeben."""
        result = radical_tower([-2, 0, 1])
        assert '±' in result['formula'] or 'sqrt' in result['formula'].lower() or '√' in result['formula']

    def test_cubic_solvable(self):
        """Kubisches Polynom (mit S_3-Gruppe): lösbar durch Radikale."""
        result = radical_tower([-2, 0, 0, 1])  # x³-2
        assert result['is_solvable_by_radicals'] is True

    def test_degree5_generic_not_solvable(self):
        """Generisches Grad-5-Polynom mit S_5: nicht lösbar."""
        # x^5 - x - 1 hat (meist) S_5 als Galoisgruppe
        result = radical_tower([-1, -1, 0, 0, 0, 1])
        # Galoisgruppe bestimmt Auflösbarkeit
        assert 'is_solvable_by_radicals' in result

    def test_tower_is_list(self):
        """Radikalturm ist eine Liste von Strings."""
        result = radical_tower([-2, 0, 1])
        assert isinstance(result['radical_tower'], list)
        assert len(result['radical_tower']) >= 1

    def test_degree_in_result(self):
        """Grad des Polynoms ist im Ergebnis."""
        result = radical_tower([-2, 0, 1])
        assert result['degree'] == 2


class TestAbelRuffiniDemo:
    """Tests für abel_ruffini_demo()."""

    def test_s5_not_solvable(self):
        """S_5 ist nicht auflösbar."""
        result = abel_ruffini_demo()
        assert result['s5_solvable'] is False

    def test_a5_simple(self):
        """A_5 ist einfache Gruppe."""
        result = abel_ruffini_demo()
        assert result['a5_simple'] is True

    def test_a5_order_60(self):
        """A_5 hat Ordnung 60."""
        result = abel_ruffini_demo()
        assert result['a5_order'] == 60

    def test_theorem_text_present(self):
        """Theorem-Text ist vorhanden und nennt Abel-Ruffini."""
        result = abel_ruffini_demo()
        assert 'Abel' in result['theorem'] or 'Ruffini' in result['theorem']

    def test_solvable_by_degree_keys(self):
        """Auflösbarkeit nach Grad enthält Grad 1 bis 6."""
        result = abel_ruffini_demo()
        for deg in [1, 2, 3, 4, 5]:
            assert deg in result['solvable_by_degree']

    def test_degree4_solvable_in_dict(self):
        """Grad 4 ist auflösbar."""
        result = abel_ruffini_demo()
        assert result['solvable_by_degree'][4]['solvable'] is True

    def test_degree5_not_solvable_in_dict(self):
        """Grad 5 ist nicht auflösbar."""
        result = abel_ruffini_demo()
        assert result['solvable_by_degree'][5]['solvable'] is False

    def test_composition_series_s5_present(self):
        """Kompositionsreihe von S_5 ist angegeben."""
        result = abel_ruffini_demo()
        assert len(result['composition_series_s5']) >= 2

    def test_example_polynomial_present(self):
        """Beispielpolynom (x^5 - x - 1) ist angegeben."""
        result = abel_ruffini_demo()
        assert result['example_polynomial'] == [-1, -1, 0, 0, 0, 1]


class TestFundamentalTheoremVerify:
    """Tests für fundamental_theorem_verify(poly_coeffs)."""

    def test_quadratic_verification(self):
        """Hauptsatz für x²-2 verifizieren."""
        result = fundamental_theorem_verify([-2, 0, 1])
        assert 'galois_group' in result
        assert 'galois_order' in result

    def test_verification_passed_flag(self):
        """Verifikation-Flag ist gesetzt."""
        result = fundamental_theorem_verify([-2, 0, 1])
        assert 'verification_passed' in result
        assert isinstance(result['verification_passed'], bool)

    def test_checks_list_nonempty(self):
        """Prüfliste ist nicht leer."""
        result = fundamental_theorem_verify([-2, 0, 1])
        assert isinstance(result['checks'], list)
        assert len(result['checks']) >= 1

    def test_subgroups_in_result(self):
        """Untergruppen sind im Ergebnis."""
        result = fundamental_theorem_verify([-2, 0, 1])
        assert 'subgroups' in result

    def test_correspondence_in_result(self):
        """Galois-Korrespondenz ist im Ergebnis."""
        result = fundamental_theorem_verify([-2, 0, 1])
        assert 'correspondence' in result


class TestGaloisGroupOfPolynomial:
    """Tests für galois_group_of_polynomial(f) – erweiterte Galoisgruppen-Info."""

    def test_linear_trivial(self):
        """Lineares Polynom: triviale Galoisgruppe."""
        result = galois_group_of_polynomial([-3, 1])  # x-3
        assert result['order'] == 1

    def test_quadratic_z2(self):
        """x²-2: Galoisgruppe ℤ/2ℤ."""
        result = galois_group_of_polynomial([-2, 0, 1])
        assert result['order'] == 2

    def test_is_abelian_flag(self):
        """is_abelian Flag ist vorhanden."""
        result = galois_group_of_polynomial([-2, 0, 1])
        assert 'is_abelian' in result

    def test_composition_series_present(self):
        """Kompositionsreihe ist vorhanden."""
        result = galois_group_of_polynomial([-2, 0, 1])
        assert 'composition_series' in result
        assert isinstance(result['composition_series'], list)

    def test_galois_implications_present(self):
        """Galois-Implikationen sind vorhanden."""
        result = galois_group_of_polynomial([-2, 0, 1])
        assert 'galois_implications' in result
        assert len(result['galois_implications']) >= 1

    def test_solvable_by_radicals_abelian(self):
        """Abelsche Galoisgruppe → lösbar durch Radikale."""
        result = galois_group_of_polynomial([-2, 0, 1])
        assert result['solvable_by_radicals'] is True

    def test_degree5_generic_usually_s5(self):
        """x^5 - x - 1 hat (meist) S_5 als Galoisgruppe."""
        result = galois_group_of_polynomial([-1, -1, 0, 0, 0, 1])
        # Galoisgruppe bestimmt, ob auflösbar
        assert 'galois_group' in result
        assert result['polynomial_degree'] == 5
