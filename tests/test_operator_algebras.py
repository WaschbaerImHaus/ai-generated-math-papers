"""
@file test_operator_algebras.py
@brief Tests für das Operatoralgebren-Modul.
@description
    Umfassende Tests für alle Klassen und Funktionen in src/operator_algebras.py:
    - CStarAlgebra (C*-Identität, Kommutativität, Spektrum)
    - VonNeumannAlgebra (Kommutant, Bikommutant, Zentrum, Faktor)
    - Projection (Idempotenz, Selbstadjungiertheit, partiell-isometrisch)
    - Gelfand-Isomorphismus
    - GNS-Konstruktion
    - Faktoren-Klassifikation
    - K-Theorie
    - Cuntz-Algebren
    - Spurklasse-Operatoren

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
import pytest
import sys
import os

# Sicherstellen, dass src/ im Pfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from operator_algebras import (
    CStarAlgebra,
    VonNeumannAlgebra,
    Projection,
    gelfand_transform_demo,
    gns_construction_demo,
    factors_classification,
    k_theory_intro,
    cuntz_algebra_demo,
    trace_class_operators,
)


# =============================================================================
# CStarAlgebra
# =============================================================================

class TestCStarAlgebra:
    """Tests für die CStarAlgebra-Klasse."""

    @pytest.fixture
    def m2_algebra(self):
        """M_2(ℂ): alle 2×2-Matrizen als Erzeuger."""
        generators = [
            np.array([[1, 0], [0, 0]], dtype=complex),
            np.array([[0, 1], [0, 0]], dtype=complex),
            np.array([[0, 0], [1, 0]], dtype=complex),
            np.array([[0, 0], [0, 1]], dtype=complex),
        ]
        return CStarAlgebra(generators)

    @pytest.fixture
    def diagonal_algebra(self):
        """Diagonale Algebra D_2: kommutativ."""
        generators = [
            np.array([[1, 0], [0, 0]], dtype=complex),
            np.array([[0, 0], [0, 1]], dtype=complex),
        ]
        return CStarAlgebra(generators)

    def test_init_valid(self, m2_algebra):
        """Initialisierung mit gültigen Erzeugern."""
        assert m2_algebra.dim == 2

    def test_init_invalid_empty(self):
        """Leere Erzeuger-Liste soll ValueError werfen."""
        with pytest.raises(ValueError):
            CStarAlgebra([])

    def test_init_inconsistent_dims(self):
        """Unterschiedliche Dimensionen sollen ValueError werfen."""
        with pytest.raises(ValueError):
            CStarAlgebra([
                np.eye(2, dtype=complex),
                np.eye(3, dtype=complex),
            ])

    def test_multiply(self, m2_algebra):
        """Matrizenmultiplikation a·b."""
        a = np.array([[1, 2], [3, 4]], dtype=complex)
        b = np.array([[5, 6], [7, 8]], dtype=complex)
        expected = a @ b
        result = m2_algebra.multiply(a, b)
        assert np.allclose(result, expected)

    def test_star(self, m2_algebra):
        """Involution a* = konjugierte Transposition."""
        a = np.array([[1 + 2j, 3], [4j, 5 - 1j]], dtype=complex)
        result = m2_algebra.star(a)
        expected = a.conj().T
        assert np.allclose(result, expected)

    def test_star_involutive(self, m2_algebra):
        """(a*)* = a."""
        a = np.array([[1 + 2j, 3], [4j, 5 - 1j]], dtype=complex)
        assert np.allclose(m2_algebra.star(m2_algebra.star(a)), a)

    def test_norm_positive_definite(self, m2_algebra):
        """Norm muss positiv definit sein: ‖a‖ ≥ 0, ‖a‖ = 0 ⟺ a = 0."""
        a = np.array([[1, 2], [3, 4]], dtype=complex)
        assert m2_algebra.norm(a) > 0
        assert m2_algebra.norm(np.zeros((2, 2), dtype=complex)) == pytest.approx(0.0)

    def test_cstar_identity(self, m2_algebra):
        """C*-Identität: ‖a*a‖ = ‖a‖² für alle a."""
        sample = [
            np.array([[1, 2], [3, 4]], dtype=complex),
            np.array([[1j, 0], [0, -1j]], dtype=complex),
            np.eye(2, dtype=complex),
        ]
        assert m2_algebra.is_cstar_algebra(sample)

    def test_cstar_identity_identity_matrix(self, m2_algebra):
        """C*-Identität für Einheitsmatrix: ‖I*I‖ = 1 = ‖I‖²."""
        I = np.eye(2, dtype=complex)
        lhs = m2_algebra.norm(m2_algebra.multiply(m2_algebra.star(I), I))
        rhs = m2_algebra.norm(I) ** 2
        assert abs(lhs - rhs) < 1e-10

    def test_is_commutative_diagonal(self, diagonal_algebra):
        """Diagonalmatrizen kommutieren."""
        assert diagonal_algebra.is_commutative()

    def test_is_not_commutative_m2(self, m2_algebra):
        """M_2(ℂ) ist nicht kommutativ."""
        assert not m2_algebra.is_commutative()

    def test_spectrum_identity(self, m2_algebra):
        """Spektrum der Einheitsmatrix = {1}."""
        I = np.eye(2, dtype=complex)
        spec = m2_algebra.spectrum(I)
        assert len(spec) == 1
        assert abs(spec[0] - 1.0) < 1e-10

    def test_spectrum_diagonal(self, m2_algebra):
        """Spektrum einer Diagonalmatrix = Diagonalelemente."""
        D = np.diag([2.0 + 0j, 3.0 + 0j])
        spec = m2_algebra.spectrum(D)
        spec_vals = sorted([abs(s) for s in spec])
        assert len(spec) == 2
        assert abs(spec_vals[0] - 2.0) < 1e-10
        assert abs(spec_vals[1] - 3.0) < 1e-10

    def test_spectral_radius(self, m2_algebra):
        """Spektralradius = max |Eigenwert|."""
        D = np.diag([2.0 + 0j, -3.0 + 0j])
        r = m2_algebra.spectral_radius(D)
        assert abs(r - 3.0) < 1e-10


# =============================================================================
# VonNeumannAlgebra
# =============================================================================

class TestVonNeumannAlgebra:
    """Tests für die VonNeumannAlgebra-Klasse."""

    @pytest.fixture
    def m2_vna(self):
        """Vollständige M_2(ℂ) als von-Neumann-Algebra."""
        generators = [
            np.array([[1, 0], [0, 0]], dtype=complex),
            np.array([[0, 1], [0, 0]], dtype=complex),
            np.array([[0, 0], [1, 0]], dtype=complex),
            np.array([[0, 0], [0, 1]], dtype=complex),
        ]
        return VonNeumannAlgebra(generators)

    @pytest.fixture
    def scalar_vna(self):
        """Skalare (Vielfache der Einheit): Zentrum von M_2."""
        generators = [np.eye(2, dtype=complex)]
        return VonNeumannAlgebra(generators)

    def test_commutant_of_all_matrices(self, m2_vna):
        """Kommutant von M_2(ℂ) = {λI}: nur Vielfache der Einheitsmatrix."""
        commutant = m2_vna.commutant(m2_vna.elements)
        # Alle Kommutant-Elemente müssen Vielfache von I sein
        n = 2
        I = np.eye(n, dtype=complex)
        for c in commutant:
            tr = np.trace(c)
            if abs(tr) > 1e-8:
                # Prüfe: c / (tr/n) ≈ I
                diff = c - (tr / n) * I
                assert np.linalg.norm(diff) < 1e-6 * np.linalg.norm(c)

    def test_center_is_scalars(self, m2_vna):
        """Zentrum von M_2(ℂ) = {λI}."""
        center = m2_vna.center()
        # Das Zentrum sollte nur Vielfache der Einheitsmatrix enthalten
        n = m2_vna.dim
        I = np.eye(n, dtype=complex)
        for z in center:
            tr = np.trace(z)
            if np.linalg.norm(z) > 1e-8:
                diff = z - (tr / n) * I
                assert np.linalg.norm(diff) < 1e-5 * np.linalg.norm(z)

    def test_is_factor_m2(self, m2_vna):
        """M_2(ℂ) ist ein Faktor (triviales Zentrum)."""
        assert m2_vna.is_factor()

    def test_bicommutant_property(self, m2_vna):
        """Bikommutant ist mindestens so groß wie der Span von M."""
        bicomm = m2_vna.bicommutant(m2_vna.elements)
        # Bikommutant sollte nicht leer sein
        assert len(bicomm) > 0

    def test_commutant_commutes(self, m2_vna):
        """Alle Elemente des Kommutanten müssen mit allen Erzeugern kommutieren."""
        commutant = m2_vna.commutant(m2_vna.elements[:2])
        for c in commutant:
            for gen in m2_vna.elements[:2]:
                comm = gen @ c - c @ gen
                assert np.linalg.norm(comm) < 1e-6


# =============================================================================
# Projection
# =============================================================================

class TestProjection:
    """Tests für die Projection-Klasse."""

    @pytest.fixture
    def proj_algebra(self):
        """Projektionsalgebra in M_2(ℂ)."""
        generators = [np.eye(2, dtype=complex)]
        return Projection(generators)

    def test_rank1_projection(self, proj_algebra):
        """Rang-1-Projektion P = diag(1,0) ist eine Projektion."""
        P = np.array([[1, 0], [0, 0]], dtype=complex)
        assert proj_algebra.is_projection(P)

    def test_identity_is_projection(self, proj_algebra):
        """Einheitsmatrix ist eine Projektion."""
        I = np.eye(2, dtype=complex)
        assert proj_algebra.is_projection(I)

    def test_zero_is_projection(self, proj_algebra):
        """Nullmatrix ist eine Projektion."""
        Z = np.zeros((2, 2), dtype=complex)
        assert proj_algebra.is_projection(Z)

    def test_non_projection(self, proj_algebra):
        """Allgemeine Matrix ist keine Projektion."""
        A = np.array([[1, 2], [0, 1]], dtype=complex)
        assert not proj_algebra.is_projection(A)

    def test_non_selfadjoint_not_projection(self, proj_algebra):
        """Nicht-selbstadjungierte Idempotente sind keine Projektionen."""
        # Schräge Projektion: A² = A aber A* ≠ A
        A = np.array([[1, 1], [0, 0]], dtype=complex)
        # A² = A prüfen
        assert np.allclose(A @ A, A)
        # Aber A* ≠ A
        assert not np.allclose(A.conj().T, A)
        assert not proj_algebra.is_projection(A)

    def test_range_space_rank1(self, proj_algebra):
        """Bildraum einer Rang-1-Projektion hat Dimension 1."""
        P = np.array([[1, 0], [0, 0]], dtype=complex)
        basis = proj_algebra.range_space(P)
        assert basis.shape[1] == 1  # Rang 1

    def test_range_space_identity(self, proj_algebra):
        """Bildraum der Einheitsmatrix hat Dimension n."""
        I = np.eye(2, dtype=complex)
        basis = proj_algebra.range_space(I)
        assert basis.shape[1] == 2

    def test_orthogonal_complement(self, proj_algebra):
        """I - P ist orthogonales Komplement."""
        P = np.array([[1, 0], [0, 0]], dtype=complex)
        Q = proj_algebra.orthogonal_complement(P)
        # Q = I - P
        expected = np.array([[0, 0], [0, 1]], dtype=complex)
        assert np.allclose(Q, expected)
        # Q ist auch eine Projektion
        assert proj_algebra.is_projection(Q)

    def test_partial_isometry(self, proj_algebra):
        """Isometrie ist partielle Isometrie."""
        # Isometrie: V = [[1,0],[0,1],[0,0]] würde 3D brauchen; hier nutzen wir unitäre Matrix
        # Für 2×2: jede unitäre Matrix ist partielle Isometrie mit V*V = I
        theta = np.pi / 4
        U = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ], dtype=complex)
        assert proj_algebra.partial_isometry(U)

    def test_not_partial_isometry(self, proj_algebra):
        """Allgemeine Matrix ist keine partielle Isometrie."""
        A = np.array([[1, 1], [0, 0]], dtype=complex)
        # A*A = [[1,0],[1,0]] ist keine Projektion
        assert not proj_algebra.partial_isometry(A)


# =============================================================================
# Gelfand-Isomorphismus
# =============================================================================

class TestGelfandTransformDemo:
    """Tests für gelfand_transform_demo()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = gelfand_transform_demo()
        assert isinstance(result, dict)

    def test_is_commutative(self):
        """Diagonalmatrizen müssen kommutativ sein."""
        result = gelfand_transform_demo()
        assert result['is_commutative'] is True

    def test_norms_equal(self):
        """Supremumsnorm = Operatornorm für Diagonalmatrizen."""
        result = gelfand_transform_demo()
        assert result['norms_equal'] is True

    def test_gelfand_values_count(self):
        """Gelfand-Werte haben die richtige Anzahl."""
        result = gelfand_transform_demo()
        assert len(result['gelfand_values']) == result['dimension']

    def test_description_present(self):
        """Beschreibung muss vorhanden sein."""
        result = gelfand_transform_demo()
        assert 'description' in result
        assert len(result['description']) > 0


# =============================================================================
# GNS-Konstruktion
# =============================================================================

class TestGNSConstructionDemo:
    """Tests für gns_construction_demo()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = gns_construction_demo()
        assert isinstance(result, dict)

    def test_basis_size(self):
        """Für M_2(ℂ): Basis hat 4 Elemente (n²)."""
        result = gns_construction_demo()
        assert result['basis_size'] == 4

    def test_gram_matrix_rank(self):
        """Gram-Matrix muss vollen Rang haben (linear unabhängige Basis)."""
        result = gns_construction_demo()
        assert result['gram_matrix_rank'] == 4

    def test_state_of_identity(self):
        """Zustand φ(I) = 1 (normierter Zustand)."""
        result = gns_construction_demo()
        assert abs(result['state_of_identity'] - 1.0) < 1e-10

    def test_norms_positive(self):
        """Alle Normen müssen positiv sein."""
        result = gns_construction_demo()
        assert result['hs_norm_of_test_op'] > 0
        assert result['operator_norm_of_test_op'] > 0


# =============================================================================
# Faktoren-Klassifikation
# =============================================================================

class TestFactorsClassification:
    """Tests für factors_classification()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = factors_classification()
        assert isinstance(result, dict)

    def test_type_i2_is_factor(self):
        """M_2(ℂ) ist ein Faktor (Type I_2)."""
        result = factors_classification()
        assert result['Type_I_2']['is_factor'] is True

    def test_type_i3_is_factor(self):
        """M_3(ℂ) ist ein Faktor (Type I_3)."""
        result = factors_classification()
        assert result['Type_I_3']['is_factor'] is True

    def test_normalized_trace_i2(self):
        """Normierte Spur Tr(I)/n = 1 für M_n(ℂ)."""
        result = factors_classification()
        assert abs(result['Type_I_2']['normalized_trace_of_I'] - 1.0) < 1e-10
        assert abs(result['Type_I_3']['normalized_trace_of_I'] - 1.0) < 1e-10

    def test_minimal_projection_trace_i2(self):
        """Minimale Projektion in M_2 hat normierte Spur 1/2."""
        result = factors_classification()
        assert abs(result['Type_I_2']['trace_of_minimal_projection'] - 0.5) < 1e-10

    def test_minimal_projection_trace_i3(self):
        """Minimale Projektion in M_3 hat normierte Spur 1/3."""
        result = factors_classification()
        assert abs(result['Type_I_3']['trace_of_minimal_projection'] - 1.0/3.0) < 1e-10

    def test_type_iii_info_present(self):
        """Typ III-Information muss vorhanden sein."""
        result = factors_classification()
        assert 'Type_III' in result
        assert 'description' in result['Type_III']


# =============================================================================
# K-Theorie
# =============================================================================

class TestKTheoryIntro:
    """Tests für k_theory_intro()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = k_theory_intro()
        assert isinstance(result, dict)

    def test_unitary_equivalence(self):
        """P1 und P2 müssen unitär äquivalent sein (gleicher K₀-Klasse)."""
        result = k_theory_intro()
        assert result['P1_P2_unitarily_equivalent'] is True

    def test_k0_classes_equal(self):
        """P1 und P2 haben die gleiche K₀-Klasse (beide Rang 1)."""
        result = k_theory_intro()
        assert result['P1_K0_class'] == result['P2_K0_class']
        assert result['P1_K0_class'] == 1

    def test_description_present(self):
        """Beschreibungen müssen vorhanden sein."""
        result = k_theory_intro()
        assert 'K0_description' in result
        assert 'K1_description' in result
        assert 'bott_periodicity' in result

    def test_examples_present(self):
        """Beispiele müssen vorhanden sein."""
        result = k_theory_intro()
        assert 'examples' in result
        assert len(result['examples']) > 0


# =============================================================================
# Cuntz-Algebren
# =============================================================================

class TestCuntzAlgebraDemo:
    """Tests für cuntz_algebra_demo()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = cuntz_algebra_demo()
        assert isinstance(result, dict)

    def test_cuntz_relations_satisfied(self):
        """Cuntz-Relationen müssen für das Toy-Modell erfüllt sein."""
        result = cuntz_algebra_demo()
        assert result['cuntz_relations_satisfied'] is True

    def test_isometry_conditions(self):
        """Sᵢ*Sᵢ = I (Isometrie-Bedingung)."""
        result = cuntz_algebra_demo()
        assert result['S1_star_S1_error'] < 1e-10
        assert result['S2_star_S2_error'] < 1e-10

    def test_completeness_relation(self):
        """S₁S₁* + S₂S₂* = I (Vollständigkeitsrelation)."""
        result = cuntz_algebra_demo()
        assert result['completeness_error'] < 1e-10

    def test_orthogonality(self):
        """S₁*S₂ = 0 (orthogonale Bilder)."""
        result = cuntz_algebra_demo()
        assert result['orthogonality_error'] < 1e-10

    def test_n_equals_2(self):
        """Demonstration für O₂ (n=2)."""
        result = cuntz_algebra_demo()
        assert result['n'] == 2

    def test_shape_info(self):
        """S₁, S₂ müssen korrekte Formen haben."""
        result = cuntz_algebra_demo()
        assert result['S1_shape'] == [4, 2]
        assert result['S2_shape'] == [4, 2]


# =============================================================================
# Spurklasse-Operatoren
# =============================================================================

class TestTraceClassOperators:
    """Tests für trace_class_operators()."""

    def test_returns_dict(self):
        """Rückgabe muss ein dict sein."""
        result = trace_class_operators(H_dim=4)
        assert isinstance(result, dict)

    def test_cyclicity(self):
        """Tr(AB) = Tr(BA) (Zyklizität der Spur)."""
        result = trace_class_operators(H_dim=4)
        assert result['cyclicity_holds'] is True
        assert result['cyclicity_error'] < 1e-8

    def test_norm_inequalities(self):
        """‖T‖ ≤ ‖T‖₂ ≤ ‖T‖₁."""
        result = trace_class_operators(H_dim=4)
        assert result['norm_inequalities_hold'] is True

    def test_positive_norms(self):
        """Alle Normen müssen nicht-negativ sein."""
        result = trace_class_operators(H_dim=4)
        assert result['trace_norm'] >= 0
        assert result['hilbert_schmidt_norm'] >= 0
        assert result['operator_norm'] >= 0

    def test_positive_operator_trace_equals_trace_norm(self):
        """Für positive Operatoren: Tr(T) = ‖T‖₁."""
        result = trace_class_operators(H_dim=4)
        assert result['trace_equals_trace_norm_for_positive_T'] is True

    def test_singular_values_count(self):
        """Anzahl der Singulärwerte = H_dim."""
        H_dim = 5
        result = trace_class_operators(H_dim=H_dim)
        assert len(result['singular_values']) == H_dim

    def test_schatten_classes_present(self):
        """Schatten-Klassen-Informationen müssen vorhanden sein."""
        result = trace_class_operators(H_dim=3)
        assert 'schatten_classes' in result
        assert 'S1_trace_class' in result['schatten_classes']
        assert 'S2_hilbert_schmidt' in result['schatten_classes']

    def test_different_dimensions(self):
        """Muss für verschiedene Dimensionen funktionieren."""
        for dim in [2, 3, 6, 8]:
            result = trace_class_operators(H_dim=dim)
            assert len(result['singular_values']) == dim
            assert result['norm_inequalities_hold'] is True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
