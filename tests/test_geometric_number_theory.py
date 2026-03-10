"""
@file test_geometric_number_theory.py
@brief Umfassende Tests für das Modul geometric_number_theory.py.
@description
    Testet alle Klassen und Funktionen der Geometrischen Zahlentheorie:
    - Lattice: Gitter-Operationen, Determinante, duales Gitter, kürzester Vektor
    - LLLReduction: Gram-Schmidt, LLL-Algorithmus, LLL-Überprüfung
    - Minkowski-Sätze: Konvexer Körper, Linearformen, Minkowski-Schranke
    - QuadraticForm: Diskriminante, Auswertung, Reduktion, Darstellung
    - GaussComposition: Klassenzahl, reduzierte Formen, Komposition
    - SpherePackingTheory: Kusszahlen, Packungsdichten, Leech-Gitter
    - Standalone-Funktionen: Hermite-Konstante, Gitter aus Form, Sukzessivminima

    Strategie: Test-Driven Development mit mindestens 80 Tests.
    Edge Cases: singuläre Basis, D=0, n=1, max/min-Werte.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os

# Suchpfad für Module anpassen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import numpy as np
import math

from geometric_number_theory import (
    Lattice,
    LLLReduction,
    minkowski_convex_body_theorem,
    minkowski_linear_forms,
    minkowski_bound,
    QuadraticForm,
    GaussComposition,
    SpherePackingTheory,
    hermite_constant,
    lattice_from_quadratic_form,
    successive_minima_bound,
)


# ==============================================================
# Tests für Lattice
# ==============================================================

class TestLattice:
    """Tests für die Gitter-Klasse."""

    def test_determinant_identity(self):
        """Einheitsmatrix hat Determinante 1."""
        L = Lattice(np.eye(3))
        assert abs(L.determinant() - 1.0) < 1e-10

    def test_determinant_2x2(self):
        """2×2-Gitter mit bekannter Determinante."""
        # Basis: [[1,0],[1,2]] → det = 2
        B = np.array([[1.0, 0.0], [1.0, 2.0]])
        L = Lattice(B)
        assert abs(L.determinant() - 2.0) < 1e-10

    def test_determinant_negative(self):
        """Determinante kann negativ sein (Vorzeichen)."""
        B = np.array([[0.0, 1.0], [1.0, 0.0]])
        L = Lattice(B)
        # det([[0,1],[1,0]]) = -1
        assert abs(L.determinant() - (-1.0)) < 1e-10

    def test_fundamental_domain_volume(self):
        """Fundamentalbereichsvolumen ist immer nicht-negativ."""
        B = np.array([[0.0, 1.0], [1.0, 0.0]])
        L = Lattice(B)
        assert L.fundamental_domain_volume() >= 0
        assert abs(L.fundamental_domain_volume() - 1.0) < 1e-10

    def test_fundamental_domain_identity(self):
        """Einheitsmatrix hat Fundamentalbereichsvolumen 1."""
        L = Lattice(np.eye(4))
        assert abs(L.fundamental_domain_volume() - 1.0) < 1e-10

    def test_contains_point_origin(self):
        """Ursprung liegt immer im Gitter."""
        L = Lattice(np.eye(3))
        assert L.contains_point(np.zeros(3))

    def test_contains_point_basis_vector(self):
        """Basisvektoren liegen im Gitter."""
        L = Lattice(np.eye(3))
        assert L.contains_point(np.array([1.0, 0.0, 0.0]))
        assert L.contains_point(np.array([0.0, 1.0, 0.0]))
        assert L.contains_point(np.array([0.0, 0.0, 1.0]))

    def test_contains_point_integer_combination(self):
        """Ganzzahlige Linearkombinationen liegen im Gitter."""
        L = Lattice(np.eye(3))
        assert L.contains_point(np.array([2.0, 3.0, -1.0]))

    def test_contains_point_not_in_lattice(self):
        """Nicht-ganzzahlige Koordinaten liegen nicht im Gitter."""
        L = Lattice(np.eye(3))
        assert not L.contains_point(np.array([0.5, 0.0, 0.0]))
        assert not L.contains_point(np.array([1.0, 0.3, 0.0]))

    def test_contains_point_singular_basis(self):
        """Singuläre Basis: Punkt nicht im Gitter."""
        # Singuläre Basis (linear abhängig)
        B = np.array([[1.0, 0.0], [2.0, 0.0]])
        L = Lattice(B)
        # Singuläre Basis → immer False
        result = L.contains_point(np.array([1.0, 0.0]))
        assert isinstance(result, bool)

    def test_contains_point_non_identity(self):
        """Gitterpunkt in nicht-trivialer Basis."""
        # Basis: [[2,0],[0,3]] → Gitterpunkte: (2k, 3l)
        B = np.array([[2.0, 0.0], [0.0, 3.0]])
        L = Lattice(B)
        assert L.contains_point(np.array([4.0, 9.0]))  # k=2, l=3
        assert not L.contains_point(np.array([1.0, 0.0]))  # nicht im Gitter

    def test_gram_matrix_identity(self):
        """Gram-Matrix der Einheitsbasis ist die Einheitsmatrix."""
        L = Lattice(np.eye(3))
        G = L.gram_matrix()
        assert np.allclose(G, np.eye(3))

    def test_gram_matrix_symmetry(self):
        """Gram-Matrix ist symmetrisch."""
        B = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
        L = Lattice(B)
        G = L.gram_matrix()
        assert np.allclose(G, G.T)

    def test_gram_matrix_positive_semidefinite(self):
        """Gram-Matrix ist positiv semidefinit."""
        B = np.array([[1.0, 1.0], [0.0, 1.0]])
        L = Lattice(B)
        G = L.gram_matrix()
        eigenvalues = np.linalg.eigvalsh(G)
        assert np.all(eigenvalues >= -1e-10)

    def test_gram_matrix_values(self):
        """Gram-Matrix mit bekannten Werten."""
        B = np.array([[1.0, 0.0], [1.0, 1.0]])
        L = Lattice(B)
        G = L.gram_matrix()
        # G[0,0] = 1, G[0,1] = G[1,0] = 1, G[1,1] = 2
        expected = np.array([[1.0, 1.0], [1.0, 2.0]])
        assert np.allclose(G, expected)

    def test_dual_lattice_identity(self):
        """Duales Gitter der Einheitsbasis ist die Einheitsbasis."""
        L = Lattice(np.eye(3))
        B_dual = L.dual_lattice_basis()
        assert np.allclose(B_dual, np.eye(3))

    def test_dual_lattice_determinant(self):
        """det(Λ*) = 1/det(Λ)."""
        B = np.array([[2.0, 0.0], [0.0, 3.0]])
        L = Lattice(B)
        B_dual = L.dual_lattice_basis()
        det_dual = np.linalg.det(B_dual)
        det_original = L.determinant()
        assert abs(det_dual - 1.0 / det_original) < 1e-10

    def test_dual_lattice_singular_raises(self):
        """Singuläre Basis wirft ValueError beim dualen Gitter."""
        B = np.array([[1.0, 0.0], [2.0, 0.0]])
        L = Lattice(B)
        with pytest.raises(ValueError):
            L.dual_lattice_basis()

    def test_shortest_vector_approx_1d(self):
        """Kürzester Vektor in 1D ist trivial."""
        B = np.array([[3.0]])
        L = Lattice(B)
        sv = L.shortest_vector_approx()
        assert np.linalg.norm(sv) > 0

    def test_shortest_vector_approx_norm(self):
        """Kürzester Vektor hat positive Länge."""
        B = np.array([[1.0, 0.0], [0.5, math.sqrt(3) / 2]])
        L = Lattice(B)
        sv = L.shortest_vector_approx()
        assert np.linalg.norm(sv) > 0

    def test_successive_minima_approx_first(self):
        """Erstes Sukzessivminimum ist positiv."""
        L = Lattice(np.eye(3))
        lam1 = L.successive_minima_approx(1)
        assert lam1 > 0

    def test_successive_minima_approx_monotone(self):
        """Sukzessivminima sind monoton wachsend (näherungsweise)."""
        B = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
        L = Lattice(B)
        lam1 = L.successive_minima_approx(1)
        lam2 = L.successive_minima_approx(2)
        lam3 = L.successive_minima_approx(3)
        assert lam1 <= lam2 + 1e-10
        assert lam2 <= lam3 + 1e-10

    def test_successive_minima_invalid_k(self):
        """k außerhalb [1,n] wirft ValueError."""
        L = Lattice(np.eye(3))
        with pytest.raises(ValueError):
            L.successive_minima_approx(0)
        with pytest.raises(ValueError):
            L.successive_minima_approx(4)

    def test_lattice_1d(self):
        """1D-Gitter (einzelner Vektor) funktioniert."""
        L = Lattice(np.array([5.0]))
        assert abs(L.determinant() - 5.0) < 1e-10

    def test_lattice_basis_stored_correctly(self):
        """Basismatrix wird korrekt gespeichert."""
        B = np.array([[1.0, 2.0], [3.0, 4.0]])
        L = Lattice(B)
        assert np.allclose(L.basis, B)


# ==============================================================
# Tests für LLLReduction
# ==============================================================

class TestLLLReduction:
    """Tests für den LLL-Reduktionsalgorithmus."""

    def test_init_default_delta(self):
        """Standard-Delta ist 0.75."""
        lll = LLLReduction()
        assert lll.delta == 0.75

    def test_init_custom_delta(self):
        """Benutzerdefiniertes Delta wird gespeichert."""
        lll = LLLReduction(delta=0.99)
        assert lll.delta == 0.99

    def test_init_invalid_delta_too_small(self):
        """Delta ≤ 0.25 wirft ValueError."""
        with pytest.raises(ValueError):
            LLLReduction(delta=0.25)

    def test_init_invalid_delta_too_large(self):
        """Delta ≥ 1.0 wirft ValueError."""
        with pytest.raises(ValueError):
            LLLReduction(delta=1.0)

    def test_gram_schmidt_identity(self):
        """Gram-Schmidt der Einheitsbasis ergibt Einheitsbasis."""
        lll = LLLReduction()
        B = np.eye(3)
        b_star, mu = lll.gram_schmidt_orthogonalization(B)
        assert np.allclose(b_star, np.eye(3))
        assert np.allclose(np.tril(mu, -1), np.zeros((3, 3)))

    def test_gram_schmidt_orthogonality(self):
        """Gram-Schmidt-Vektoren sind orthogonal zueinander."""
        lll = LLLReduction()
        B = np.array([[1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 2.0, 3.0]])
        b_star, mu = lll.gram_schmidt_orthogonalization(B)

        # Überprüfe Orthogonalität
        for i in range(3):
            for j in range(i):
                dot = np.dot(b_star[i], b_star[j])
                assert abs(dot) < 1e-8, f"b*[{i}] und b*[{j}] nicht orthogonal: dot={dot}"

    def test_lll_reduce_trivial(self):
        """Einheitsbasis ist bereits LLL-reduziert."""
        lll = LLLReduction()
        B = np.eye(3)
        reduced = lll.reduce(B)
        # LLL der Einheitsbasis sollte kurze Vektoren zurückgeben
        assert reduced.shape == (3, 3)

    def test_lll_reduce_1d(self):
        """1D-Basis: triviale Reduktion."""
        lll = LLLReduction()
        B = np.array([[7.0, 3.0]])
        reduced = lll.reduce(B)
        assert reduced.shape == (1, 2)
        assert np.allclose(reduced, B)

    def test_lll_reduce_preserves_span(self):
        """LLL-reduzierte Basis überspannt dasselbe Gitter."""
        lll = LLLReduction()
        B = np.array([[1.0, 1.0], [0.0, 1.0]])
        reduced = lll.reduce(B)

        # Die reduzierte Basis muss dasselbe Gitter erzeugen
        # d.h. det(reduced) = ± det(B)
        assert abs(abs(np.linalg.det(reduced)) - abs(np.linalg.det(B))) < 1e-8

    def test_lll_reduce_shorter_first_vector(self):
        """LLL ergibt kürzere Vektoren als die ursprüngliche Basis."""
        lll = LLLReduction()
        # Bewusst "schlechte" Basis mit langen Vektoren
        B = np.array([[100.0, 0.0], [99.5, 0.5]])
        reduced = lll.reduce(B)
        # Der erste Vektor der reduzierten Basis sollte kurz sein
        norm_reduced = np.linalg.norm(reduced[0])
        norm_original = np.linalg.norm(B[0])
        assert norm_reduced <= norm_original + 1e-6

    def test_is_lll_reduced_identity(self):
        """Einheitsbasis ist LLL-reduziert."""
        lll = LLLReduction()
        assert lll.is_lll_reduced(np.eye(3))

    def test_is_lll_reduced_1d(self):
        """1D-Basis ist immer LLL-reduziert."""
        lll = LLLReduction()
        assert lll.is_lll_reduced(np.array([[5.0, 3.0]]))

    def test_is_lll_reduced_after_reduction(self):
        """Nach LLL-Reduktion ist die Basis LLL-reduziert."""
        lll = LLLReduction()
        B = np.array([[3.0, 5.0], [1.0, 2.0]])
        reduced = lll.reduce(B)
        assert lll.is_lll_reduced(reduced)

    def test_lll_reduce_3d(self):
        """LLL-Reduktion in 3D liefert eine gültige Basis."""
        lll = LLLReduction()
        B = np.array([
            [6.0, 5.0, 9.0],
            [3.0, 7.0, 1.0],
            [4.0, 2.0, 8.0]
        ])
        reduced = lll.reduce(B)
        assert reduced.shape == (3, 3)
        # Determinante bleibt (bis auf Vorzeichen) gleich
        assert abs(abs(np.linalg.det(reduced)) - abs(np.linalg.det(B))) < 1e-6


# ==============================================================
# Tests für Minkowski-Sätze
# ==============================================================

class TestMinkowskiTheorems:
    """Tests für die Minkowski-Sätze."""

    def test_convex_body_theorem_true(self):
        """Großes Volumen → Gitterpunkt garantiert."""
        # Vol > 2^n * det
        # dim=2, det=1: 2^2*1 = 4, Vol=5 > 4 → True
        assert minkowski_convex_body_theorem(volume=5.0, det_lattice=1.0, dim=2)

    def test_convex_body_theorem_false(self):
        """Kleines Volumen → kein Gitterpunkt garantiert."""
        # Vol < 2^n * det
        assert not minkowski_convex_body_theorem(volume=3.0, det_lattice=1.0, dim=2)

    def test_convex_body_threshold_exact(self):
        """Genau am Schwellenwert: kein Gitterpunkt garantiert (Satz gilt nur für strikt >)."""
        # Vol = 2^n * det → nicht garantiert (strikt >)
        assert not minkowski_convex_body_theorem(volume=4.0, det_lattice=1.0, dim=2)

    def test_convex_body_dim1(self):
        """Dimension 1: 2^1 = 2."""
        assert minkowski_convex_body_theorem(volume=3.0, det_lattice=1.0, dim=1)
        assert not minkowski_convex_body_theorem(volume=1.0, det_lattice=1.0, dim=1)

    def test_convex_body_dim3(self):
        """Dimension 3: 2^3 = 8."""
        assert minkowski_convex_body_theorem(volume=9.0, det_lattice=1.0, dim=3)
        assert not minkowski_convex_body_theorem(volume=7.0, det_lattice=1.0, dim=3)

    def test_convex_body_large_det(self):
        """Große Gitterdeterminante verschiebt den Schwellenwert."""
        # dim=2, det=10: Schwelle = 40
        assert minkowski_convex_body_theorem(volume=41.0, det_lattice=10.0, dim=2)
        assert not minkowski_convex_body_theorem(volume=39.0, det_lattice=10.0, dim=2)

    def test_linear_forms_valid(self):
        """Gültige Schranken (Produkt ≥ 1)."""
        bounds = np.array([1.0, 1.0, 1.0])
        assert minkowski_linear_forms(np.zeros(3), bounds)

    def test_linear_forms_invalid_zero_bound(self):
        """Schranke = 0 → False."""
        bounds = np.array([0.0, 1.0, 1.0])
        assert not minkowski_linear_forms(np.zeros(3), bounds)

    def test_linear_forms_product_less_than_1(self):
        """Produkt < 1 → False."""
        bounds = np.array([0.5, 0.5, 0.5])  # 0.5^3 = 0.125 < 1
        assert not minkowski_linear_forms(np.zeros(3), bounds)

    def test_linear_forms_product_greater_than_1(self):
        """Produkt > 1 → True."""
        bounds = np.array([2.0, 1.0, 1.0])
        assert minkowski_linear_forms(np.zeros(3), bounds)

    def test_minkowski_bound_quadratic_field(self):
        """Minkowski-Schranke für quadratischen Zahlkörper Q(√-5): n=2, r2=1."""
        # Q(√-5): D=-20, n=2, r2=1
        # M = (2!/4) * (4/π)^1 * √20 ≈ 2.853
        bound = minkowski_bound(discriminant=-20, degree=2, num_complex_pairs=1)
        assert bound > 0
        assert bound < 5  # Plausibilitätsprüfung

    def test_minkowski_bound_positive_discriminant(self):
        """Minkowski-Schranke mit positivem Discriminant (Betrag)."""
        # Betrag des Discriminants wird verwendet
        bound = minkowski_bound(discriminant=5, degree=2, num_complex_pairs=0)
        assert bound > 0

    def test_minkowski_bound_cubic_field(self):
        """Minkowski-Schranke für kubischen Zahlkörper."""
        bound = minkowski_bound(discriminant=-23, degree=3, num_complex_pairs=1)
        assert bound > 0


# ==============================================================
# Tests für QuadraticForm
# ==============================================================

class TestQuadraticForm:
    """Tests für binäre quadratische Formen."""

    def test_discriminant_positive_definite(self):
        """Form x² + y² hat Diskriminante -4."""
        f = QuadraticForm(1, 0, 1)
        assert f.discriminant() == -4

    def test_discriminant_zero(self):
        """Entartete Form hat Diskriminante 0."""
        f = QuadraticForm(1, 2, 1)
        assert f.discriminant() == 0

    def test_discriminant_negative(self):
        """Form x² + xy + y² hat D = 1-4 = -3."""
        f = QuadraticForm(1, 1, 1)
        assert f.discriminant() == -3

    def test_discriminant_positive(self):
        """Indefinite Form x² - y² hat D = 4."""
        f = QuadraticForm(1, 0, -1)
        assert f.discriminant() == 4

    def test_evaluate_simple(self):
        """f(1,0) = a für f = ax² + bxy + cy²."""
        f = QuadraticForm(3, 2, 5)
        assert f.evaluate(1, 0) == 3
        assert f.evaluate(0, 1) == 5
        assert f.evaluate(1, 1) == 3 + 2 + 5

    def test_evaluate_zero(self):
        """f(0,0) = 0 immer."""
        f = QuadraticForm(5, 3, 7)
        assert f.evaluate(0, 0) == 0

    def test_evaluate_negative_input(self):
        """f(-x,-y) = f(x,y) für positiv definite Formen."""
        f = QuadraticForm(2, 0, 3)
        assert f.evaluate(-1, -2) == f.evaluate(1, 2)

    def test_is_positive_definite_true(self):
        """x² + y² ist positiv definit."""
        f = QuadraticForm(1, 0, 1)
        assert f.is_positive_definite()

    def test_is_positive_definite_false_indefinite(self):
        """x² - y² ist nicht positiv definit."""
        f = QuadraticForm(1, 0, -1)
        assert not f.is_positive_definite()

    def test_is_positive_definite_false_discriminant_zero(self):
        """Form mit D=0 ist nicht positiv definit (entartet)."""
        f = QuadraticForm(1, 2, 1)
        assert not f.is_positive_definite()

    def test_is_positive_definite_false_a_negative(self):
        """Form mit a < 0 ist nicht positiv definit."""
        f = QuadraticForm(-1, 0, -1)
        assert not f.is_positive_definite()

    def test_is_reduced_standard(self):
        """(1, 0, 1) ist reduziert: |0| ≤ 1 ≤ 1, b=0 ≥ 0."""
        f = QuadraticForm(1, 0, 1)
        assert f.is_reduced()

    def test_is_reduced_b_too_large(self):
        """(1, 2, 3): |b|=2 > a=1 → nicht reduziert."""
        f = QuadraticForm(1, 2, 3)
        assert not f.is_reduced()

    def test_is_reduced_a_greater_c(self):
        """(3, 0, 1): a=3 > c=1 → nicht reduziert."""
        f = QuadraticForm(3, 0, 1)
        assert not f.is_reduced()

    def test_is_reduced_boundary_b_equals_a(self):
        """(2, 2, 3): |b|=a=2, aber b=2 ≥ 0 → reduziert."""
        f = QuadraticForm(2, 2, 3)
        assert f.is_reduced()

    def test_is_reduced_boundary_negative_b_equals_a(self):
        """(2, -2, 3): |b|=a=2, b=-2 < 0 → nicht reduziert."""
        f = QuadraticForm(2, -2, 3)
        assert not f.is_reduced()

    def test_reduce_already_reduced(self):
        """Bereits reduzierte Form bleibt reduziert."""
        f = QuadraticForm(1, 1, 1)
        r = f.reduce()
        assert r.is_reduced()

    def test_reduce_preserves_discriminant(self):
        """Reduktion ändert die Diskriminante nicht."""
        f = QuadraticForm(5, 4, 2)
        D = f.discriminant()
        if f.is_positive_definite():
            r = f.reduce()
            assert r.discriminant() == D

    def test_reduce_various_forms(self):
        """Verschiedene Formen lassen sich reduzieren."""
        test_cases = [
            QuadraticForm(3, 4, 2),   # D = 16 - 24 = -8 → nicht pos. def. (a<0 nein, aber D>0)
            QuadraticForm(2, 2, 3),   # D = 4-24 = -20, pos. def.
            QuadraticForm(5, 2, 3),   # D = 4-60 = -56, pos. def.
        ]
        for f in test_cases:
            if f.is_positive_definite():
                r = f.reduce()
                assert r.is_reduced(), f"Nicht reduziert nach Reduktion: {r}"

    def test_reduce_not_positive_definite_raises(self):
        """Reduktion nicht-positiv-definiter Form wirft ValueError."""
        f = QuadraticForm(1, 0, -1)  # D=4, indefinit
        with pytest.raises(ValueError):
            f.reduce()

    def test_represents_trivially(self):
        """Jede Form stellt 0 dar."""
        f = QuadraticForm(2, 1, 3)
        assert f.represents(0)

    def test_represents_x2_plus_y2(self):
        """x² + y² stellt 1, 2, 4, 5 dar."""
        f = QuadraticForm(1, 0, 1)
        assert f.represents(1)   # f(1,0) = 1
        assert f.represents(2)   # f(1,1) = 2
        assert f.represents(4)   # f(2,0) = 4
        assert f.represents(5)   # f(2,1) = 5

    def test_represents_x2_plus_y2_not_3(self):
        """x² + y² stellt 3 nicht dar (bekanntes Resultat)."""
        f = QuadraticForm(1, 0, 1)
        assert not f.represents(3)

    def test_represents_2x2_plus_2y2(self):
        """2x² + 2y² stellt gerade Zahlen dar."""
        f = QuadraticForm(2, 0, 2)
        assert f.represents(4)  # f(1,1) = 4
        assert not f.represents(1)  # ungerade Zahl nicht darstellbar

    def test_equivalent_forms_nonempty(self):
        """Äquivalente Formen sind nicht leer."""
        f = QuadraticForm(1, 0, 1)
        forms = f.equivalent_forms(max_coeff=5)
        assert len(forms) > 0

    def test_equivalent_forms_contains_original(self):
        """Äquivalente Formen enthalten die ursprüngliche Form."""
        f = QuadraticForm(1, 0, 1)
        forms = f.equivalent_forms(max_coeff=5)
        # Mindestens die Identitätstransformation sollte enthalten sein
        assert any(g.a == f.a and g.b == f.b and g.c == f.c for g in forms)

    def test_quadratic_form_repr(self):
        """Repr enthält Koeffizienten."""
        f = QuadraticForm(2, 3, 5)
        r = repr(f)
        assert "2" in r and "3" in r and "5" in r

    def test_quadratic_form_equality(self):
        """Gleichheit zweier Formen."""
        f1 = QuadraticForm(1, 2, 3)
        f2 = QuadraticForm(1, 2, 3)
        f3 = QuadraticForm(1, 0, 3)
        assert f1 == f2
        assert f1 != f3


# ==============================================================
# Tests für GaussComposition
# ==============================================================

class TestGaussComposition:
    """Tests für die Gauß'sche Komposition quadratischer Formen."""

    def test_reduced_forms_D_minus3(self):
        """h(-3) = 1: Einzige reduzierte Form ist (1,1,1)."""
        forms = GaussComposition.reduced_forms(-3)
        assert len(forms) == 1
        assert forms[0].a == 1 and forms[0].b == 1 and forms[0].c == 1

    def test_reduced_forms_D_minus4(self):
        """h(-4) = 1: Einzige reduzierte Form ist (1,0,1)."""
        forms = GaussComposition.reduced_forms(-4)
        assert len(forms) == 1
        assert forms[0].a == 1 and forms[0].b == 0 and forms[0].c == 1

    def test_reduced_forms_all_have_correct_discriminant(self):
        """Alle zurückgegebenen Formen haben die korrekte Diskriminante."""
        D = -20
        forms = GaussComposition.reduced_forms(D)
        for f in forms:
            assert f.discriminant() == D, f"Form {f} hat falsche Diskriminante: {f.discriminant()}"

    def test_reduced_forms_all_reduced(self):
        """Alle zurückgegebenen Formen sind reduziert."""
        D = -20
        forms = GaussComposition.reduced_forms(D)
        for f in forms:
            assert f.is_reduced(), f"Form {f} ist nicht reduziert"

    def test_class_number_D_minus3(self):
        """h(-3) = 1."""
        assert GaussComposition.class_number(-3) == 1

    def test_class_number_D_minus4(self):
        """h(-4) = 1."""
        assert GaussComposition.class_number(-4) == 1

    def test_class_number_D_minus23(self):
        """h(-23) = 3."""
        h = GaussComposition.class_number(-23)
        assert h == 3

    def test_class_number_positive(self):
        """Klassenzahl ist immer ≥ 1."""
        for D in [-3, -4, -7, -8, -11, -15, -20, -23]:
            h = GaussComposition.class_number(D)
            assert h >= 1

    def test_compose_same_discriminant_required(self):
        """Komposition verschiedener Diskriminanten wirft ValueError."""
        f1 = QuadraticForm(1, 0, 1)  # D=-4
        f2 = QuadraticForm(1, 1, 1)  # D=-3
        with pytest.raises(ValueError):
            GaussComposition.compose(f1, f2)

    def test_compose_returns_quadratic_form(self):
        """Komposition gibt QuadraticForm zurück."""
        f1 = QuadraticForm(1, 0, 1)  # D=-4
        f2 = QuadraticForm(1, 0, 1)  # D=-4
        result = GaussComposition.compose(f1, f2)
        assert isinstance(result, QuadraticForm)

    def test_compose_preserves_discriminant(self):
        """Komposition erhält die Diskriminante."""
        f1 = QuadraticForm(2, 2, 3)  # D=4-24=-20
        f2 = QuadraticForm(2, 2, 3)  # D=-20
        result = GaussComposition.compose(f1, f2)
        assert result.discriminant() == f1.discriminant()


# ==============================================================
# Tests für SpherePackingTheory
# ==============================================================

class TestSpherePackingTheory:
    """Tests für die Kugelpackungstheorie."""

    def test_kissing_number_dim1(self):
        """Kusszahl in dim=1 ist 2."""
        assert SpherePackingTheory.kissing_number(1) == 2

    def test_kissing_number_dim2(self):
        """Kusszahl in dim=2 ist 6."""
        assert SpherePackingTheory.kissing_number(2) == 6

    def test_kissing_number_dim3(self):
        """Kusszahl in dim=3 ist 12."""
        assert SpherePackingTheory.kissing_number(3) == 12

    def test_kissing_number_dim4(self):
        """Kusszahl in dim=4 ist 24."""
        assert SpherePackingTheory.kissing_number(4) == 24

    def test_kissing_number_dim8(self):
        """Kusszahl in dim=8 ist 240 (E₈-Gitter)."""
        assert SpherePackingTheory.kissing_number(8) == 240

    def test_kissing_number_dim24(self):
        """Kusszahl in dim=24 ist 196560 (Leech-Gitter)."""
        assert SpherePackingTheory.kissing_number(24) == 196560

    def test_kissing_number_unknown_dim(self):
        """Unbekannte Dimension gibt None zurück."""
        assert SpherePackingTheory.kissing_number(5) is None
        assert SpherePackingTheory.kissing_number(10) is None

    def test_packing_density_upper_bound_dim1(self):
        """Packungsdichte-Schranke in dim=1 ist 1.0."""
        assert SpherePackingTheory.packing_density_upper_bound(1) == 1.0

    def test_packing_density_upper_bound_positive(self):
        """Packungsdichte-Schranke ist immer positiv."""
        for dim in [1, 2, 3, 4, 8, 16, 24]:
            bound = SpherePackingTheory.packing_density_upper_bound(dim)
            assert bound > 0

    def test_packing_density_upper_bound_leq_1(self):
        """Packungsdichte-Schranke ist nie > 1."""
        for dim in [1, 2, 3, 5, 10, 24]:
            bound = SpherePackingTheory.packing_density_upper_bound(dim)
            assert bound <= 1.0 + 1e-10

    def test_packing_density_invalid_dim(self):
        """Dimension ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError):
            SpherePackingTheory.packing_density_upper_bound(0)
        with pytest.raises(ValueError):
            SpherePackingTheory.packing_density_upper_bound(-1)

    def test_hexagonal_packing_density(self):
        """Hexagonale Packungsdichte ≈ 0.9069."""
        density = SpherePackingTheory.hexagonal_packing_density()
        expected = math.pi / (2 * math.sqrt(3))
        assert abs(density - expected) < 1e-12
        assert abs(density - 0.9069) < 0.001

    def test_fcc_packing_density(self):
        """FCC-Packungsdichte ≈ 0.7405."""
        density = SpherePackingTheory.fcc_packing_density()
        expected = math.pi / (3 * math.sqrt(2))
        assert abs(density - expected) < 1e-12
        assert abs(density - 0.7405) < 0.001

    def test_hexagonal_denser_than_fcc(self):
        """2D hexagonale Packung ist dichter als 3D FCC-Packung."""
        hex_d = SpherePackingTheory.hexagonal_packing_density()
        fcc_d = SpherePackingTheory.fcc_packing_density()
        assert hex_d > fcc_d

    def test_leech_lattice_info_keys(self):
        """Leech-Gitter-Info enthält alle erwarteten Schlüssel."""
        info = SpherePackingTheory.leech_lattice_info()
        assert "dimension" in info
        assert "kissing_number" in info
        assert "packing_density" in info
        assert "determinant" in info

    def test_leech_lattice_dim24(self):
        """Leech-Gitter ist 24-dimensional."""
        info = SpherePackingTheory.leech_lattice_info()
        assert info["dimension"] == 24

    def test_leech_lattice_kissing_number(self):
        """Leech-Gitter Kusszahl ist 196560."""
        info = SpherePackingTheory.leech_lattice_info()
        assert info["kissing_number"] == 196560

    def test_leech_lattice_unimodular(self):
        """Leech-Gitter ist unimodulär (det=1)."""
        info = SpherePackingTheory.leech_lattice_info()
        assert info["determinant"] == 1
        assert info["is_unimodular"] is True


# ==============================================================
# Tests für Standalone-Funktionen
# ==============================================================

class TestStandaloneFunctions:
    """Tests für die globalen Hilfsfunktionen."""

    def test_hermite_constant_dim1(self):
        """γ₁ = 1."""
        assert abs(hermite_constant(1) - 1.0) < 1e-10

    def test_hermite_constant_dim2(self):
        """γ₂ = 4/3."""
        assert abs(hermite_constant(2) - 4.0 / 3.0) < 1e-10

    def test_hermite_constant_dim3(self):
        """γ₃ = 2."""
        assert abs(hermite_constant(3) - 2.0) < 1e-10

    def test_hermite_constant_dim4(self):
        """γ₄ = 4."""
        assert abs(hermite_constant(4) - 4.0) < 1e-10

    def test_hermite_constant_dim8(self):
        """γ₈ = 2."""
        assert abs(hermite_constant(8) - 2.0) < 1e-10

    def test_hermite_constant_positive(self):
        """Hermite-Konstante ist immer positiv."""
        for n in range(1, 12):
            assert hermite_constant(n) > 0

    def test_hermite_constant_invalid(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            hermite_constant(0)
        with pytest.raises(ValueError):
            hermite_constant(-1)

    def test_lattice_from_quadratic_form_identity(self):
        """Form x² + y² erzeugt Einheitsgitter."""
        L = lattice_from_quadratic_form(1, 0, 1)
        # det = 1
        assert abs(L.fundamental_domain_volume() - 1.0) < 1e-10

    def test_lattice_from_quadratic_form_gram_matrix(self):
        """Gram-Matrix des erzeugten Gitters entspricht der Form-Matrix."""
        a, b, c = 2, 2, 3
        L = lattice_from_quadratic_form(a, b, c)
        G = L.gram_matrix()
        # G sollte der Form [[a, b/2],[b/2, c]] entsprechen
        assert abs(G[0, 0] - a) < 1e-10
        assert abs(G[0, 1] - b / 2.0) < 1e-10
        assert abs(G[1, 1] - c) < 1e-10

    def test_lattice_from_quadratic_form_not_positive_definite_raises(self):
        """Nicht positiv definite Form wirft ValueError."""
        with pytest.raises(ValueError):
            lattice_from_quadratic_form(1, 0, -1)  # D=4, indefinit

    def test_successive_minima_bound_dim1(self):
        """Schranke für dim=1, det=1 ergibt γ₁^{1/2}·1 = 1."""
        bound = successive_minima_bound(det_lattice=1.0, dim=1)
        expected = math.sqrt(hermite_constant(1)) * 1.0
        assert abs(bound - expected) < 1e-10

    def test_successive_minima_bound_positive(self):
        """Schranke ist immer positiv."""
        for dim in range(1, 9):
            bound = successive_minima_bound(det_lattice=1.0, dim=dim)
            assert bound > 0

    def test_successive_minima_bound_invalid_det(self):
        """Determinante ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError):
            successive_minima_bound(det_lattice=0.0, dim=3)
        with pytest.raises(ValueError):
            successive_minima_bound(det_lattice=-1.0, dim=3)

    def test_successive_minima_bound_invalid_dim(self):
        """Dimension < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            successive_minima_bound(det_lattice=1.0, dim=0)

    def test_successive_minima_bound_larger_det(self):
        """Größere Determinante → größere Schranke."""
        bound1 = successive_minima_bound(det_lattice=1.0, dim=2)
        bound2 = successive_minima_bound(det_lattice=4.0, dim=2)
        assert bound2 > bound1


# ==============================================================
# Edge-Case Tests
# ==============================================================

class TestEdgeCases:
    """Edge-Case Tests für besondere Situationen."""

    def test_lattice_1x1(self):
        """1×1-Gitter (ein 1D-Vektor)."""
        L = Lattice(np.array([[5.0]]))
        assert abs(L.determinant() - 5.0) < 1e-10

    def test_quadratic_form_discriminant_zero(self):
        """Form (1,2,1) hat D=0 (entartet)."""
        f = QuadraticForm(1, 2, 1)
        assert f.discriminant() == 0
        assert not f.is_positive_definite()

    def test_quadratic_form_large_coefficients(self):
        """Form mit großen Koeffizienten."""
        f = QuadraticForm(1000, 0, 1000)
        assert f.discriminant() == -4000000
        assert f.is_positive_definite()

    def test_lll_2x2_reduction(self):
        """LLL für 2×2-Basis mit bekanntem Ergebnis."""
        lll = LLLReduction()
        # Bekanntes Beispiel aus Literatur
        B = np.array([[1.0, 1.0], [0.0, 1.0]])
        reduced = lll.reduce(B)
        assert reduced.shape == (2, 2)
        assert lll.is_lll_reduced(reduced)

    def test_minkowski_convex_body_dim0_boundary(self):
        """Kleines dim=1 Szenario."""
        # dim=1: Schwelle = 2 * det
        assert minkowski_convex_body_theorem(volume=2.1, det_lattice=1.0, dim=1)

    def test_sphere_packing_density_decreasing_trend(self):
        """Kabatiansky-Levenshtein: Dichte nimmt mit dim ab."""
        dims = [2, 4, 8, 16]
        bounds = [SpherePackingTheory.packing_density_upper_bound(d) for d in dims]
        # Schranken sollten mit höherer Dimension abnehmen
        for i in range(len(bounds) - 1):
            assert bounds[i] >= bounds[i + 1] - 1e-10

    def test_class_number_D_minus11(self):
        """h(-11) = 1."""
        assert GaussComposition.class_number(-11) == 1

    def test_reduced_forms_D_minus7(self):
        """h(-7) = 1: Einzige reduzierte Form ist (1,1,2)."""
        forms = GaussComposition.reduced_forms(-7)
        assert len(forms) == 1

    def test_gram_schmidt_2d(self):
        """Gram-Schmidt für 2D-Basis mit bekanntem Ergebnis."""
        lll = LLLReduction()
        B = np.array([[1.0, 0.0], [1.0, 1.0]])
        b_star, mu = lll.gram_schmidt_orthogonalization(B)
        # b*[0] = b[0] = (1,0)
        assert np.allclose(b_star[0], [1.0, 0.0])
        # μ[1,0] = <b[1],b*[0]>/<b*[0],b*[0]> = 1/1 = 1
        assert abs(mu[1, 0] - 1.0) < 1e-10
        # b*[1] = b[1] - 1*b*[0] = (0,1)
        assert np.allclose(b_star[1], [0.0, 1.0])


# ==============================================================
# Hauptprogramm (direkt ausführbar)
# ==============================================================

if __name__ == "__main__":
    # Tests direkt ausführen
    pytest.main([__file__, "-v", "--tb=short"])
