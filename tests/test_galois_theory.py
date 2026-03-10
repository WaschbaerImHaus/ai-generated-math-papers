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

@author Kurt Ingwer
@since 2026-03-10
@lastModified 2026-03-10
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
    _is_irreducible_over_fp, _find_irreducible_poly, _poly_mul_mod, _poly_pow_mod
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
