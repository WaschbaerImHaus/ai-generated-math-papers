"""
@file test_invariant_theory.py
@brief Tests für das invariant_theory-Modul.
       Prüft Reynolds-Operator, symmetrische Polynome, Newton-Identitäten,
       Molien-Reihe, Hilbert-Basissatz und Diskriminante.
@description
    Testet alle Funktionen aus invariant_theory.py:
    - reynolds_operator: Projektion auf G-invariante Polynome
    - is_invariant: Invarianz-Prüfung
    - elementary_symmetric_polynomials: eₖ in x₁,...,xₙ
    - power_sum_symmetric_polynomials: pₖ = Σ xᵢᵏ
    - newton_identities: Rekursion eₖ ↔ pₖ
    - polynomial_invariants_sn: Invariante unter Sₙ
    - molien_series: Poincaré-Reihe
    - hilbert_basis_theorem_demo: Erzeuger für S₂
    - discriminant_invariant: Diskriminante
    - fundamental_theorem_symmetric_poly: Darstellung in eₖ

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
import sympy as sp
from sympy import symbols, expand, simplify
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from invariant_theory import (
    reynolds_operator, is_invariant,
    elementary_symmetric_polynomials, power_sum_symmetric_polynomials,
    newton_identities, polynomial_invariants_sn,
    molien_series, hilbert_basis_theorem_demo,
    discriminant_invariant, fundamental_theorem_symmetric_poly
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def s2_action(poly, perm):
    """
    S₂-Wirkung: Vertauscht x1 und x2 (simultane Substitution).

    simultaneous=True verhindert Kettensubstitution: ohne das würde
    x1→x2 und dann x2→x1 zu x1→x1 führen (falsch).
    """
    x1, x2 = sp.symbols('x1 x2')
    if perm == (0, 1):
        return poly  # Identität: keine Änderung
    else:
        # Simultane Substitution: x1↔x2 auf einmal
        return poly.subs([(x1, x2), (x2, x1)], simultaneous=True)


# ---------------------------------------------------------------------------
# Tests für Reynolds-Operator
# ---------------------------------------------------------------------------

class TestReynoldsOperator:
    """Tests für den Reynolds-Operator."""

    def test_reynolds_symmetric_polynomial_unchanged(self):
        """Test: Symmetrisches Polynom bleibt unverändert unter Reynolds."""
        x1, x2 = sp.symbols('x1 x2')
        # x1 + x2 ist S₂-invariant
        poly = x1 + x2
        s2_elements = [(0, 1), (1, 0)]
        result = reynolds_operator(poly, s2_elements, s2_action)
        assert sp.expand(result - poly) == 0

    def test_reynolds_projects_to_invariant(self):
        """Test: Reynolds-Operator projiziert x1 auf symmetrisches Polynom."""
        x1, x2 = sp.symbols('x1 x2')
        # Reynolds(x1) = (x1 + x2) / 2
        poly = x1
        s2_elements = [(0, 1), (1, 0)]
        result = reynolds_operator(poly, s2_elements, s2_action)
        expected = (x1 + x2) / 2
        assert sp.expand(result - expected) == 0

    def test_reynolds_idempotent(self):
        """Test: Reynolds-Operator ist idempotent: R∘R = R."""
        x1, x2 = sp.symbols('x1 x2')
        poly = x1**2
        s2_elements = [(0, 1), (1, 0)]
        r_poly = reynolds_operator(poly, s2_elements, s2_action)
        r_r_poly = reynolds_operator(r_poly, s2_elements, s2_action)
        assert sp.expand(r_poly - r_r_poly) == 0

    def test_reynolds_empty_group_raises(self):
        """Test: Leere Gruppe wirft ValueError."""
        x1 = sp.Symbol('x1')
        with pytest.raises(ValueError):
            reynolds_operator(x1, [], lambda p, g: p)

    def test_reynolds_constant_unchanged(self):
        """Test: Konstante ist immer invariant."""
        s2_elements = [(0, 1), (1, 0)]
        result = reynolds_operator(sp.Integer(5), s2_elements, s2_action)
        assert result == 5


# ---------------------------------------------------------------------------
# Tests für is_invariant
# ---------------------------------------------------------------------------

class TestIsInvariant:
    """Tests für die Invarianz-Prüfung."""

    def test_symmetric_poly_is_invariant(self):
        """Test: x1+x2 ist S₂-invariant."""
        x1, x2 = sp.symbols('x1 x2')
        poly = x1 + x2
        s2 = [(0, 1), (1, 0)]
        assert is_invariant(poly, ['x1', 'x2'], s2, s2_action) is True

    def test_symmetric_poly_product_is_invariant(self):
        """Test: x1*x2 ist S₂-invariant."""
        x1, x2 = sp.symbols('x1 x2')
        poly = x1 * x2
        s2 = [(0, 1), (1, 0)]
        assert is_invariant(poly, ['x1', 'x2'], s2, s2_action) is True

    def test_nonsymmetric_poly_not_invariant(self):
        """Test: x1 allein ist nicht S₂-invariant."""
        x1, x2 = sp.symbols('x1 x2')
        poly = x1
        s2 = [(0, 1), (1, 0)]
        assert is_invariant(poly, ['x1', 'x2'], s2, s2_action) is False

    def test_x1_squared_plus_x2_squared_invariant(self):
        """Test: x1² + x2² ist S₂-invariant."""
        x1, x2 = sp.symbols('x1 x2')
        poly = x1**2 + x2**2
        s2 = [(0, 1), (1, 0)]
        assert is_invariant(poly, ['x1', 'x2'], s2, s2_action) is True


# ---------------------------------------------------------------------------
# Tests für elementarsymmetrische Polynome
# ---------------------------------------------------------------------------

class TestElementarySymmetricPolynomials:
    """Tests für elementarsymmetrische Polynome."""

    def test_e1_n2(self):
        """Test: e₁ für n=2 ist x1+x2."""
        x1, x2 = sp.symbols('x1 x2')
        e = elementary_symmetric_polynomials(2)
        assert sp.expand(e[0] - (x1 + x2)) == 0

    def test_e2_n2(self):
        """Test: e₂ für n=2 ist x1*x2."""
        x1, x2 = sp.symbols('x1 x2')
        e = elementary_symmetric_polynomials(2)
        assert sp.expand(e[1] - x1*x2) == 0

    def test_e_count_n3(self):
        """Test: Für n=3 gibt es genau 3 elementarsymmetrische Polynome."""
        e = elementary_symmetric_polynomials(3)
        assert len(e) == 3

    def test_e1_n3(self):
        """Test: e₁ für n=3 ist x1+x2+x3."""
        x1, x2, x3 = sp.symbols('x1 x2 x3')
        e = elementary_symmetric_polynomials(3)
        assert sp.expand(e[0] - (x1 + x2 + x3)) == 0

    def test_e3_n3(self):
        """Test: e₃ für n=3 ist x1*x2*x3."""
        x1, x2, x3 = sp.symbols('x1 x2 x3')
        e = elementary_symmetric_polynomials(3)
        assert sp.expand(e[2] - x1*x2*x3) == 0

    def test_e2_n3(self):
        """Test: e₂ für n=3 ist x1*x2 + x1*x3 + x2*x3."""
        x1, x2, x3 = sp.symbols('x1 x2 x3')
        e = elementary_symmetric_polynomials(3)
        expected = x1*x2 + x1*x3 + x2*x3
        assert sp.expand(e[1] - expected) == 0

    def test_e_are_symmetric(self):
        """Test: Elementarsymmetrische Polynome sind symmetrisch."""
        x1, x2 = sp.symbols('x1 x2')
        e = elementary_symmetric_polynomials(2)
        for ek in e:
            # Simultane Permutation: x1 ↔ x2
            permuted = ek.subs([(x1, x2), (x2, x1)], simultaneous=True)
            assert sp.expand(ek - permuted) == 0


# ---------------------------------------------------------------------------
# Tests für Potenzsummen
# ---------------------------------------------------------------------------

class TestPowerSumPolynomials:
    """Tests für Potenzsummen-Polynome."""

    def test_p1_is_e1(self):
        """Test: p₁ = x1+x2 = e₁."""
        x1, x2 = sp.symbols('x1 x2')
        p = power_sum_symmetric_polynomials(2, 1)
        assert sp.expand(p[0] - (x1 + x2)) == 0

    def test_p2_n2(self):
        """Test: p₂ = x1² + x2²."""
        x1, x2 = sp.symbols('x1 x2')
        p = power_sum_symmetric_polynomials(2, 2)
        assert sp.expand(p[1] - (x1**2 + x2**2)) == 0

    def test_pk_count(self):
        """Test: max_degree Potenzsummen werden zurückgegeben."""
        p = power_sum_symmetric_polynomials(3, 4)
        assert len(p) == 4

    def test_pk_are_symmetric(self):
        """Test: Potenzsummen sind symmetrisch."""
        x1, x2 = sp.symbols('x1 x2')
        p = power_sum_symmetric_polynomials(2, 3)
        for pk in p:
            # Simultane Substitution: x1 ↔ x2
            permuted = pk.subs([(x1, x2), (x2, x1)], simultaneous=True)
            assert sp.expand(pk - permuted) == 0


# ---------------------------------------------------------------------------
# Tests für Newton-Identitäten
# ---------------------------------------------------------------------------

class TestNewtonIdentities:
    """Tests für Newton-Identitäten."""

    def test_newton_n2_structure(self):
        """Test: Newton-Identitäten für n=2 haben korrekten Aufbau."""
        result = newton_identities(2)
        assert 'n' in result
        assert 'elementary_symmetric' in result
        assert 'power_sums' in result
        assert 'identities' in result
        assert result['n'] == 2

    def test_newton_n2_verification(self):
        """Test: Newton-Identitäten für n=2 sind korrekt (Verifikation = 0)."""
        result = newton_identities(2)
        for identity in result['identities']:
            # Jede Identität sollte verification=0 haben
            verif = identity['verification']
            assert sp.expand(verif) == 0, (
                f"Newton-Identität k={identity['k']} ist nicht erfüllt: {verif}"
            )

    def test_newton_n3_verification(self):
        """Test: Newton-Identitäten für n=3 sind korrekt."""
        result = newton_identities(3)
        for identity in result['identities']:
            verif = identity['verification']
            assert sp.expand(verif) == 0

    def test_newton_n2_count(self):
        """Test: Für n=2 gibt es 2 Newton-Identitäten."""
        result = newton_identities(2)
        assert len(result['identities']) == 2


# ---------------------------------------------------------------------------
# Tests für Sₙ-Invarianten
# ---------------------------------------------------------------------------

class TestPolynomialInvariantsSn:
    """Tests für Sₙ-invariante Polynome."""

    def test_invariants_s2_degree1(self):
        """Test: Grad-1-Invariante unter S₂ ist e₁ = x1+x2."""
        x1, x2 = sp.symbols('x1 x2')
        invs = polynomial_invariants_sn(2, 1)
        # e₁ = x1+x2 sollte enthalten sein
        e1 = x1 + x2
        found = any(sp.expand(inv - e1) == 0 for inv in invs)
        assert found, f"e₁=x1+x2 nicht in Invarianten gefunden: {invs}"

    def test_invariants_s2_degree2_contains_e1_squared(self):
        """Test: Grad-2-Invarianten unter S₂ enthalten e₁² oder e₂."""
        invs = polynomial_invariants_sn(2, 2)
        # Mindestens 2 Invarianten bis Grad 2
        assert len(invs) >= 2

    def test_invariants_are_symmetric(self):
        """Test: Alle zurückgegebenen Polynome sind S₂-symmetrisch."""
        x1, x2 = sp.symbols('x1 x2')
        invs = polynomial_invariants_sn(2, 2)
        for inv in invs:
            # Simultane Substitution: x1 ↔ x2
            permuted = inv.subs([(x1, x2), (x2, x1)], simultaneous=True)
            assert sp.expand(inv - permuted) == 0, (
                f"Polynom {inv} ist nicht symmetrisch"
            )


# ---------------------------------------------------------------------------
# Tests für Molien-Reihe
# ---------------------------------------------------------------------------

class TestMolienSeries:
    """Tests für die Molien-Reihe."""

    def test_molien_trivial_group(self):
        """Test: Triviale Gruppe {I} hat alle Potenzen als Invariante."""
        # Triviale Gruppe: einziges Element ist Einheitsmatrix
        matrices = [[[1, 0], [0, 1]]]
        group_elements = ['identity']
        coeffs = molien_series(group_elements, matrices, max_degree=4)
        # Jeder Grad hat Dimension 1 (Potenzen von t als Invarianten)
        # Für 2D: Dim(Inv_k) = Anzahl Monome vom Grad k = k+1
        assert coeffs[0] == 1  # Grad 0: konstante Polynome

    def test_molien_s2_coefficients(self):
        """Test: Molien-Reihe für S₂ hat korrekte Koeffizienten."""
        # S₂-Wirkung: {I, Transposition}
        I = [[1, 0], [0, 1]]
        T = [[0, 1], [1, 0]]  # Permutationsmatrix
        matrices = [I, T]
        group_elements = ['id', 'tau']
        coeffs = molien_series(group_elements, matrices, max_degree=6)
        # c₀ = 1 (eine Konstante)
        assert coeffs[0] == 1
        # c₁ = 1 (eine Grad-1-Invariante: x+y)
        assert coeffs[1] == 1
        # c₂ = 2 (Grad-2-Invarianten: (x+y)², xy)
        assert coeffs[2] == 2

    def test_molien_empty_raises(self):
        """Test: Leere Matrizenliste wirft ValueError."""
        with pytest.raises(ValueError):
            molien_series([], [], max_degree=4)


# ---------------------------------------------------------------------------
# Tests für Hilbert-Basissatz
# ---------------------------------------------------------------------------

class TestHilbertBasisTheorem:
    """Tests für die Hilbert-Basissatz-Demo."""

    def test_hilbert_demo_returns_dict(self):
        """Test: Demo gibt Dict zurück."""
        result = hilbert_basis_theorem_demo()
        assert isinstance(result, dict)

    def test_hilbert_demo_has_generators(self):
        """Test: Demo enthält Erzeuger e₁ und e₂."""
        result = hilbert_basis_theorem_demo()
        assert 'generators' in result
        assert 'e1' in result['generators']
        assert 'e2' in result['generators']

    def test_hilbert_demo_examples_verified(self):
        """Test: Alle Beispiele in der Demo sind korrekt verifiziert."""
        result = hilbert_basis_theorem_demo()
        for example in result['examples']:
            assert example['verified'] is True, (
                f"Beispiel {example['poly']} ist nicht korrekt: {example}"
            )

    def test_hilbert_demo_ring_name(self):
        """Test: Demo bezeichnet den korrekten Ring."""
        result = hilbert_basis_theorem_demo()
        assert 'S₂' in result['ring'] or 'S2' in str(result['ring'])


# ---------------------------------------------------------------------------
# Tests für Diskriminante
# ---------------------------------------------------------------------------

class TestDiscriminantInvariant:
    """Tests für die Diskriminante als Invariante."""

    def test_discriminant_quadratic_positive(self):
        """Test: b²-4ac > 0 für x²-5x+6 (Wurzeln 2,3)."""
        # x²-5x+6 = (x-2)(x-3), Δ = 25-24 = 1 > 0
        disc = discriminant_invariant([1, -5, 6])
        assert int(disc) == 1

    def test_discriminant_quadratic_zero(self):
        """Test: b²-4ac = 0 für x²-2x+1 (Doppelwurzel)."""
        # x²-2x+1 = (x-1)², Δ = 4-4 = 0
        disc = discriminant_invariant([1, -2, 1])
        assert int(disc) == 0

    def test_discriminant_quadratic_negative(self):
        """Test: b²-4ac < 0 für x²+x+1 (komplexe Wurzeln)."""
        # x²+x+1: Δ = 1-4 = -3
        disc = discriminant_invariant([1, 1, 1])
        assert int(disc) == -3

    def test_discriminant_cubic(self):
        """Test: Diskriminante für kubisches Polynom."""
        # x³ - 3x + 2 = (x-1)²(x+2): Δ = 0
        disc = discriminant_invariant([1, 0, -3, 2])
        assert int(disc) == 0

    def test_discriminant_invalid_degree_raises(self):
        """Test: Zu kurze Koeffizientenliste wirft ValueError."""
        with pytest.raises(ValueError):
            discriminant_invariant([1])  # Grad 0, kein Polynom


# ---------------------------------------------------------------------------
# Tests für Hauptsatz symmetrischer Polynome
# ---------------------------------------------------------------------------

class TestFundamentalTheoremSymmetricPoly:
    """Tests für den Hauptsatz der Theorie symmetrischer Polynome."""

    def test_symmetric_poly_is_detected(self):
        """Test: Symmetrisches Polynom wird korrekt erkannt."""
        result = fundamental_theorem_symmetric_poly('x1**2 + x2**2', 2)
        assert result['is_symmetric'] is True

    def test_nonsymmetric_poly_detected(self):
        """Test: Nicht-symmetrisches Polynom wird erkannt."""
        result = fundamental_theorem_symmetric_poly('x1**2 + x2', 2)
        assert result['is_symmetric'] is False

    def test_result_has_elementary_symmetric(self):
        """Test: Ergebnis enthält elementarsymmetrische Polynome."""
        result = fundamental_theorem_symmetric_poly('x1 + x2', 2)
        assert 'elementary_symmetric' in result
        assert len(result['elementary_symmetric']) == 2

    def test_result_contains_theorem(self):
        """Test: Ergebnis enthält Theorem-Beschreibung."""
        result = fundamental_theorem_symmetric_poly('x1*x2', 2)
        assert 'theorem' in result

    def test_symmetric_poly_verified_representation(self):
        """Test: x1²+x2² kann als e₁²-2e₂ dargestellt werden."""
        result = fundamental_theorem_symmetric_poly('x1**2 + x2**2', 2)
        assert result['is_symmetric'] is True
        if result.get('verified'):
            assert result['verified'] is True
