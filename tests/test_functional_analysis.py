"""
@file test_functional_analysis.py
@brief Umfassende Tests für das Modul functional_analysis.py.
@description
    Testet alle Klassen und Funktionen des Funktionalanalysis-Moduls:
    - NormedSpace: Norm-Axiome, Abstände, Einheitskugel
    - BanachSpace: Cauchy-Folgen, lp-Normen, Separabilität
    - HilbertSpace: Skalarprodukt, Gram-Schmidt, Projektion, Parseval, Bessel
    - LinearOperator: Anwendung, Normen, Adjungierte, Spektraleigenschaften
    - CompactOperator: Spektralzerlegung
    - Spektraltheorie: spectrum, spectral_radius, spectral_theorem_demo
    - Fredholm-Alternative
    - Funktionenräume: ContinuousFunctions, SobolevSpace, L2Space
    - Sätze: Hahn-Banach, offene Abbildung, abgeschlossener Graph,
             gleichmäßige Beschränktheit, schwache Konvergenz,
             Riesz-Darstellung, kompakte Operatoren
    - Fixpunktsätze: Banach, Schauder

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math

import numpy as np
import pytest

# Suchpfad so setzen, dass das src-Verzeichnis gefunden wird
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from functional_analysis import (
    NormedSpace,
    BanachSpace,
    HilbertSpace,
    LinearOperator,
    CompactOperator,
    spectrum,
    spectral_radius,
    spectral_theorem_demo,
    fredholm_alternative,
    ContinuousFunctions,
    SobolevSpace,
    L2Space,
    hahn_banach_theorem_demo,
    open_mapping_theorem_demo,
    closed_graph_theorem_demo,
    uniform_boundedness_principle_demo,
    weak_convergence_demo,
    riesz_representation_theorem,
    compact_operator_theory,
    banach_fixed_point,
    schauder_fixed_point_demo,
)


# ============================================================
# Hilfsfunktionen für die Tests
# ============================================================

def _euclidean_norm(v: np.ndarray) -> float:
    """Euklidische L2-Norm als Hilfsfunktion."""
    return float(np.linalg.norm(v))


def _standard_inner_product(u: np.ndarray, v: np.ndarray) -> float:
    """Standard-Skalarprodukt als Hilfsfunktion."""
    return float(np.dot(u, v))


def _make_hilbert_space(vectors: list) -> HilbertSpace:
    """Erstellt einen Standard-Hilbert-Raum (ℝⁿ) mit den gegebenen Vektoren."""
    return HilbertSpace(vectors, _standard_inner_product)


# ============================================================
# Tests: NormedSpace
# ============================================================

class TestNormedSpace:
    """
    @brief Tests für die Klasse NormedSpace.
    @description Prüft Norm-Berechnung, Norm-Axiome und Abstandsfunktion.
    """

    def test_norm_nonnegative(self):
        """Norm eines beliebigen Vektors muss nicht-negativ sein."""
        space = NormedSpace([[3.0, -4.0]], _euclidean_norm)
        assert space.norm([3.0, -4.0]) >= 0.0

    def test_norm_zero_vector(self):
        """Norm des Nullvektors muss 0 sein."""
        space = NormedSpace([[0.0, 0.0]], _euclidean_norm)
        assert space.norm([0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)

    def test_norm_euclidean_value(self):
        """Euklidische Norm von (3, 4) muss 5 sein (Pythagoras)."""
        space = NormedSpace([[3.0, 4.0]], _euclidean_norm)
        assert space.norm([3.0, 4.0]) == pytest.approx(5.0, rel=1e-10)

    def test_norm_homogeneity(self):
        """‖α·v‖ = |α|·‖v‖ für alle Skalare α."""
        space = NormedSpace([[1.0, 2.0, 3.0]], _euclidean_norm)
        v = [1.0, 2.0, 3.0]
        alpha = 3.0
        assert space.norm([alpha * x for x in v]) == pytest.approx(
            alpha * space.norm(v), rel=1e-10
        )

    def test_triangle_inequality(self):
        """Dreiecksungleichung: ‖u + v‖ ≤ ‖u‖ + ‖v‖."""
        space = NormedSpace([[1.0, 0.0], [0.0, 1.0]], _euclidean_norm)
        u = np.array([1.0, 0.0])
        v = np.array([0.0, 1.0])
        assert space.norm(u + v) <= space.norm(u) + space.norm(v) + 1e-12

    def test_distance_symmetry(self):
        """Abstand ist symmetrisch: d(u,v) = d(v,u)."""
        space = NormedSpace([[1.0, 2.0], [3.0, 4.0]], _euclidean_norm)
        u = [1.0, 2.0]
        v = [3.0, 4.0]
        assert space.distance(u, v) == pytest.approx(space.distance(v, u), rel=1e-10)

    def test_distance_self_is_zero(self):
        """Abstand eines Vektors zu sich selbst ist 0."""
        space = NormedSpace([[5.0, -2.0]], _euclidean_norm)
        v = [5.0, -2.0]
        assert space.distance(v, v) == pytest.approx(0.0, abs=1e-12)

    def test_distance_concrete_value(self):
        """Konkreter Abstand: d([0,0], [3,4]) = 5."""
        space = NormedSpace([[0.0, 0.0], [3.0, 4.0]], _euclidean_norm)
        assert space.distance([0.0, 0.0], [3.0, 4.0]) == pytest.approx(5.0, rel=1e-10)

    def test_is_normed_space_returns_true(self):
        """is_normed_space() muss True für eine gültige Norm zurückgeben."""
        space = NormedSpace([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], _euclidean_norm)
        assert space.is_normed_space() is True

    def test_unit_ball_points_have_norm_one(self):
        """Alle Punkte auf der Einheitskugel haben Norm ≈ 1."""
        space = NormedSpace([[1.0, 0.0]], _euclidean_norm)
        points = space.unit_ball(20)
        assert len(points) == 20
        for p in points:
            assert space.norm(p) == pytest.approx(1.0, rel=1e-6)

    def test_l1_norm(self):
        """L1-Norm von [1, -2, 3] = 6."""
        l1_norm = lambda v: float(np.sum(np.abs(v)))
        space = NormedSpace([[1.0, -2.0, 3.0]], l1_norm)
        assert space.norm([1.0, -2.0, 3.0]) == pytest.approx(6.0, rel=1e-10)


# ============================================================
# Tests: BanachSpace
# ============================================================

class TestBanachSpace:
    """
    @brief Tests für die Klasse BanachSpace.
    @description Prüft Cauchy-Folgen, lp-Normen und Separabilität.
    """

    def test_cauchy_sequence_converging(self):
        """Eine konvergierende Folge ist eine Cauchy-Folge."""
        space = BanachSpace([[1.0]], _euclidean_norm)
        # xₙ = [1/n, 1/n] → [0, 0]
        sequence = [[1.0 / (n + 1), 1.0 / (n + 1)] for n in range(100)]
        assert space.cauchy_sequence_test(sequence, tol=0.05) is True

    def test_cauchy_sequence_not_cauchy(self):
        """Eine oszillierende Folge ist keine Cauchy-Folge."""
        space = BanachSpace([[1.0]], _euclidean_norm)
        # Abwechselnd [0,0] und [1,0] — keine Cauchy-Folge
        sequence = [[float(i % 2), 0.0] for i in range(20)]
        assert space.cauchy_sequence_test(sequence, tol=0.1) is False

    def test_single_element_sequence_is_cauchy(self):
        """Eine einelementige Folge ist trivialerweise eine Cauchy-Folge."""
        space = BanachSpace([[1.0, 2.0]], _euclidean_norm)
        assert space.cauchy_sequence_test([[1.0, 2.0]]) is True

    def test_is_separable(self):
        """Endlich-dimensionale Räume sind separabel."""
        space = BanachSpace([[1.0, 0.0]], _euclidean_norm)
        assert space.is_separable() is True

    def test_lp_norm_l1(self):
        """ℓ¹-Norm von [1, 2, 3] = 6."""
        space = BanachSpace([[1.0, 2.0, 3.0]], _euclidean_norm)
        assert space.lp_norm([1.0, 2.0, 3.0], 1.0) == pytest.approx(6.0, rel=1e-10)

    def test_lp_norm_l2(self):
        """ℓ²-Norm von [3, 4] = 5."""
        space = BanachSpace([[3.0, 4.0]], _euclidean_norm)
        assert space.lp_norm([3.0, 4.0], 2.0) == pytest.approx(5.0, rel=1e-10)

    def test_lp_norm_linf(self):
        """ℓ∞-Norm von [1, -5, 3] = 5."""
        space = BanachSpace([[1.0, -5.0, 3.0]], _euclidean_norm)
        assert space.lp_norm([1.0, -5.0, 3.0], float('inf')) == pytest.approx(5.0, rel=1e-10)

    def test_lp_norm_invalid_p(self):
        """Negativer p-Wert muss einen ValueError auslösen."""
        space = BanachSpace([[1.0]], _euclidean_norm)
        with pytest.raises(ValueError):
            space.lp_norm([1.0, 2.0], -1.0)

    def test_lp_norm_zero_vector(self):
        """ℓᵖ-Norm des Nullvektors ist 0."""
        space = BanachSpace([[0.0, 0.0, 0.0]], _euclidean_norm)
        assert space.lp_norm([0.0, 0.0, 0.0], 2.0) == pytest.approx(0.0, abs=1e-12)


# ============================================================
# Tests: HilbertSpace
# ============================================================

class TestHilbertSpace:
    """
    @brief Tests für die Klasse HilbertSpace.
    @description Prüft Skalarprodukt, Orthogonalität, Gram-Schmidt,
                 Projektion, Parseval und Bessel-Ungleichung.
    """

    def test_inner_product_standard(self):
        """Standard-Skalarprodukt ⟨[1,2],[3,4]⟩ = 1·3 + 2·4 = 11."""
        H = _make_hilbert_space([[1.0, 2.0], [3.0, 4.0]])
        assert H.inner_product([1.0, 2.0], [3.0, 4.0]) == pytest.approx(11.0, rel=1e-10)

    def test_inner_product_symmetry(self):
        """⟨u,v⟩ = ⟨v,u⟩ (Symmetrie)."""
        H = _make_hilbert_space([[1.0, 3.0], [2.0, -1.0]])
        u, v = [1.0, 3.0], [2.0, -1.0]
        assert H.inner_product(u, v) == pytest.approx(H.inner_product(v, u), rel=1e-10)

    def test_inner_product_self_nonnegative(self):
        """⟨v,v⟩ ≥ 0 für alle v."""
        H = _make_hilbert_space([[2.0, -3.0, 1.0]])
        v = [2.0, -3.0, 1.0]
        assert H.inner_product(v, v) >= 0.0

    def test_inner_product_zero_vector(self):
        """⟨0,v⟩ = 0."""
        H = _make_hilbert_space([[0.0, 0.0], [1.0, 1.0]])
        assert H.inner_product([0.0, 0.0], [1.0, 1.0]) == pytest.approx(0.0, abs=1e-12)

    def test_is_orthogonal_perpendicular(self):
        """[1,0] und [0,1] sind orthogonal."""
        H = _make_hilbert_space([[1.0, 0.0], [0.0, 1.0]])
        assert H.is_orthogonal([1.0, 0.0], [0.0, 1.0]) is True

    def test_is_orthogonal_not_perpendicular(self):
        """[1,1] und [1,1] sind nicht orthogonal."""
        H = _make_hilbert_space([[1.0, 1.0]])
        assert H.is_orthogonal([1.0, 1.0], [1.0, 1.0]) is False

    def test_norm_from_inner_product(self):
        """Norm aus Skalarprodukt: ‖[3,4]‖ = 5."""
        H = _make_hilbert_space([[3.0, 4.0]])
        assert H.norm([3.0, 4.0]) == pytest.approx(5.0, rel=1e-10)

    def test_gram_schmidt_orthonormality(self):
        """Gram-Schmidt liefert orthonormale Vektoren (⟨eᵢ,eⱼ⟩ = δᵢⱼ)."""
        H = _make_hilbert_space([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
        onb = H.gram_schmidt_orthonormalize([[1.0, 1.0, 0.0],
                                              [0.0, 1.0, 1.0],
                                              [1.0, 0.0, 1.0]])
        # Jeder Vektor muss Norm 1 haben
        for e in onb:
            assert H.norm(e) == pytest.approx(1.0, rel=1e-8)
        # Verschiedene Vektoren müssen orthogonal sein
        for i in range(len(onb)):
            for j in range(i + 1, len(onb)):
                assert H.inner_product(onb[i], onb[j]) == pytest.approx(0.0, abs=1e-8)

    def test_gram_schmidt_single_vector(self):
        """Gram-Schmidt eines einzelnen Vektors ergibt normierten Vektor."""
        H = _make_hilbert_space([[3.0, 4.0]])
        onb = H.gram_schmidt_orthonormalize([[3.0, 4.0]])
        assert len(onb) == 1
        assert H.norm(onb[0]) == pytest.approx(1.0, rel=1e-10)

    def test_projection_onto_subspace(self):
        """Projektion von [1,1,0] auf span{[1,0,0]} = [1,0,0]."""
        H = _make_hilbert_space([[1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        proj = H.projection([1.0, 1.0, 0.0], [[1.0, 0.0, 0.0]])
        assert proj == pytest.approx([1.0, 0.0, 0.0], abs=1e-8)

    def test_projection_residual_orthogonal(self):
        """v - P(v) muss orthogonal zur Basis sein."""
        H = _make_hilbert_space([[2.0, 3.0, 1.0], [1.0, 0.0, 0.0]])
        v = [2.0, 3.0, 1.0]
        basis = [[1.0, 0.0, 0.0]]
        proj = H.projection(v, basis)
        residual = [v[i] - proj[i] for i in range(len(v))]
        assert H.inner_product(residual, basis[0]) == pytest.approx(0.0, abs=1e-8)

    def test_pythagoras_orthogonal_vectors(self):
        """Pythagoras gilt für orthogonale Vektoren: ‖u+v‖² = ‖u‖² + ‖v‖²."""
        H = _make_hilbert_space([[1.0, 0.0], [0.0, 1.0]])
        assert H.pythagoras([1.0, 0.0], [0.0, 1.0]) is True

    def test_pythagoras_nonorthogonal_vectors(self):
        """Pythagoras gilt NICHT für nicht-orthogonale Vektoren."""
        H = _make_hilbert_space([[1.0, 1.0], [1.0, 0.0]])
        assert H.pythagoras([1.0, 1.0], [1.0, 0.0]) is False

    def test_parseval_identity_onb(self):
        """Parseval-Identität gilt für eine vollständige ONB."""
        H = _make_hilbert_space([[1.0, 0.0, 0.0],
                                   [0.0, 1.0, 0.0],
                                   [0.0, 0.0, 1.0]])
        # Standardbasis ist eine ONB für ℝ³
        onb = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        v = [1.0, 2.0, 3.0]
        assert H.parseval_identity(v, onb) is True

    def test_bessel_inequality(self):
        """Bessel-Ungleichung: Σ|⟨v,eᵢ⟩|² ≤ ‖v‖²."""
        H = _make_hilbert_space([[1.0, 0.0, 0.0]])
        # Partielle ONB (nur ein Vektor)
        partial_onb = [[1.0, 0.0, 0.0]]
        v = [1.0, 2.0, 3.0]
        assert H.bessel_inequality(v, partial_onb) is True

    def test_riesz_representation(self):
        """Riesz-Repräsentant von f(x) = 2x₁ auf ℝ² ist [2, 0]."""
        H = _make_hilbert_space([[1.0, 0.0], [0.0, 1.0]])
        f = lambda x: 2.0 * x[0]
        basis = [[1.0, 0.0], [0.0, 1.0]]
        y = H.riesz_representation(f, basis)
        # f(x) = ⟨x, y⟩ → y = [2, 0]
        assert y == pytest.approx([2.0, 0.0], abs=1e-8)


# ============================================================
# Tests: LinearOperator
# ============================================================

class TestLinearOperator:
    """
    @brief Tests für die Klasse LinearOperator.
    @description Prüft Anwendung, Operatornorm, Adjungierte und
                 strukturelle Eigenschaften.
    """

    def test_apply_identity(self):
        """Identitätsoperator ändert keinen Vektor."""
        I = LinearOperator(np.eye(3))
        v = [1.0, 2.0, 3.0]
        assert I.apply(v) == pytest.approx(v, rel=1e-10)

    def test_apply_matrix_multiply(self):
        """Matrixmultiplikation: [[1,2],[3,4]] · [1,1] = [3,7]."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        op = LinearOperator(A)
        result = op.apply([1.0, 1.0])
        assert result == pytest.approx([3.0, 7.0], rel=1e-10)

    def test_apply_zero_vector(self):
        """T · 0 = 0 für jeden Operator."""
        A = np.array([[5.0, -1.0], [2.0, 3.0]])
        op = LinearOperator(A)
        result = op.apply([0.0, 0.0])
        assert result == pytest.approx([0.0, 0.0], abs=1e-12)

    def test_is_bounded_finite_matrix(self):
        """Endlich-dimensionale Operatoren sind immer beschränkt."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        op = LinearOperator(A)
        vectors = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        assert op.is_bounded(vectors) is True

    def test_is_compact(self):
        """Endlich-dimensionale Operatoren sind immer kompakt."""
        op = LinearOperator(np.array([[2.0, 1.0], [0.0, 3.0]]))
        assert op.is_compact() is True

    def test_operator_norm_l2(self):
        """Spektrale Operatornorm = größter Singulärwert."""
        A = np.diag([3.0, 1.0, 0.5])  # Eigenwerte = Singulärwerte (diagonal, positiv)
        op = LinearOperator(A)
        assert op.operator_norm() == pytest.approx(3.0, rel=1e-8)

    def test_adjoint_is_transpose(self):
        """Adjungierter Operator = Transponierte (reelle Matrizen)."""
        A = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        op = LinearOperator(A)
        adj = op.adjoint()
        assert np.allclose(adj.matrix, A.T, atol=1e-10)

    def test_is_self_adjoint_symmetric(self):
        """Symmetrische Matrizen sind selbstadjungiert."""
        A = np.array([[2.0, 1.0], [1.0, 3.0]])
        op = LinearOperator(A)
        assert op.is_self_adjoint() is True

    def test_is_self_adjoint_asymmetric(self):
        """Unsymmetrische Matrizen sind nicht selbstadjungiert."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        op = LinearOperator(A)
        assert op.is_self_adjoint() is False

    def test_is_unitary_rotation(self):
        """Rotationsmatrix ist unitär."""
        theta = math.pi / 4
        R = np.array([[math.cos(theta), -math.sin(theta)],
                      [math.sin(theta),  math.cos(theta)]])
        op = LinearOperator(R)
        assert op.is_unitary() is True

    def test_is_normal_self_adjoint(self):
        """Selbstadjungierte Operatoren sind normal."""
        A = np.array([[1.0, 2.0], [2.0, 3.0]])
        op = LinearOperator(A)
        assert op.is_normal() is True

    def test_is_projection_identity(self):
        """Identitätsmatrix ist eine Projektion (I² = I)."""
        op = LinearOperator(np.eye(3))
        assert op.is_projection() is True

    def test_is_projection_ortho_proj(self):
        """Orthogonalprojektion P = vvᵀ/‖v‖² ist idempotent."""
        # P = e₁·e₁ᵀ = [[1,0],[0,0]]
        P = np.array([[1.0, 0.0], [0.0, 0.0]])
        op = LinearOperator(P)
        assert op.is_projection() is True

    def test_kernel_invertible_matrix(self):
        """Invertierbarer Operator hat trivialen Kern ({0})."""
        A = np.array([[1.0, 0.0], [0.0, 1.0]])
        op = LinearOperator(A)
        kern = op.kernel()
        # Kern ist leer (numerisch keine Vektoren mit Singulärwert ≈ 0)
        assert len(kern) == 0

    def test_kernel_singular_matrix(self):
        """Singulärer Operator hat nicht-trivialen Kern."""
        # Rang-1-Matrix: alle Zeilen sind vielfache von [1,1]
        A = np.array([[1.0, 1.0], [2.0, 2.0]])
        op = LinearOperator(A)
        kern = op.kernel()
        assert len(kern) >= 1
        # Kern-Vektor muss T·v ≈ 0 ergeben
        for kv in kern:
            result = op.apply(kv)
            assert np.linalg.norm(result) == pytest.approx(0.0, abs=1e-8)

    def test_range_space_full_rank(self):
        """Bildraum einer vollen Rangmatrix hat Dimension = n."""
        A = np.array([[1.0, 0.0], [0.0, 1.0]])
        op = LinearOperator(A)
        rng = op.range_space()
        assert len(rng) == 2


# ============================================================
# Tests: CompactOperator
# ============================================================

class TestCompactOperator:
    """
    @brief Tests für die Klasse CompactOperator.
    @description Prüft Spektralzerlegung und Kompaktheitseigenschaften.
    """

    def test_spectral_decomposition_symmetric(self):
        """Spektralzerlegung einer symmetrischen Matrix ist gültig."""
        A = np.array([[2.0, 1.0], [1.0, 2.0]])
        op = CompactOperator(A)
        result = op.spectral_decomposition()
        assert result['decomposition_valid'] is True
        assert len(result['eigenvalues']) == 2
        assert len(result['eigenvectors']) == 2

    def test_spectral_decomposition_eigenvalues_real(self):
        """Eigenwerte einer symmetrischen Matrix sind reell."""
        A = np.array([[3.0, 1.0], [1.0, 3.0]])
        op = CompactOperator(A)
        result = op.spectral_decomposition()
        for ev in result['eigenvalues']:
            assert isinstance(ev, float)

    def test_spectral_decomposition_sorted_by_magnitude(self):
        """Eigenwerte sind nach fallendem Absolutwert sortiert."""
        A = np.diag([5.0, 1.0, 3.0])
        op = CompactOperator(A)
        result = op.spectral_decomposition()
        evs = [abs(x) for x in result['eigenvalues']]
        assert evs == sorted(evs, reverse=True)

    def test_compact_operator_inherits_is_compact(self):
        """CompactOperator erbt is_compact() → immer True."""
        op = CompactOperator(np.eye(2))
        assert op.is_compact() is True


# ============================================================
# Tests: Spektraltheorie (Funktionen)
# ============================================================

class TestSpectralTheory:
    """
    @brief Tests für die Spektraltheorie-Funktionen.
    @description Prüft spectrum(), spectral_radius(), spectral_theorem_demo().
    """

    def test_spectrum_returns_dict_keys(self):
        """spectrum() gibt dict mit korrekten Schlüsseln zurück."""
        op = LinearOperator(np.eye(3))
        result = spectrum(op)
        for key in ('point_spectrum', 'eigenvalues', 'spectral_radius',
                    'continuous_spectrum', 'residual_spectrum'):
            assert key in result

    def test_spectrum_identity_eigenvalues(self):
        """Spektrum der Identität = {1}."""
        op = LinearOperator(np.eye(3))
        result = spectrum(op)
        # Alle Eigenwerte ≈ 1 (komplex)
        for ev in result['eigenvalues']:
            assert abs(ev - 1.0) == pytest.approx(0.0, abs=1e-8)

    def test_spectrum_diagonal_matrix(self):
        """Spektrum einer Diagonalmatrix = Diagonaleinträge."""
        A = np.diag([2.0, -3.0, 5.0])
        op = LinearOperator(A)
        result = spectrum(op)
        evs = sorted([ev.real for ev in result['eigenvalues']])
        assert evs == pytest.approx(sorted([-3.0, 2.0, 5.0]), abs=1e-8)

    def test_spectral_radius_value(self):
        """Spektralradius der Diagonalmatrix = größter Absolutwert der Eigenwerte."""
        A = np.diag([2.0, -5.0, 3.0])
        op = LinearOperator(A)
        assert spectral_radius(op) == pytest.approx(5.0, rel=1e-8)

    def test_spectral_radius_identity(self):
        """Spektralradius der Identität = 1."""
        op = LinearOperator(np.eye(4))
        assert spectral_radius(op) == pytest.approx(1.0, rel=1e-8)

    def test_spectral_theorem_demo_returns_dict(self):
        """spectral_theorem_demo() gibt dict mit korrekten Schlüsseln zurück."""
        A = np.array([[2.0, 1.0], [1.0, 2.0]])
        result = spectral_theorem_demo(A)
        for key in ('eigenvalues', 'eigenvectors', 'reconstruction_error', 'is_valid'):
            assert key in result

    def test_spectral_theorem_demo_reconstruction(self):
        """Spektralsatz: Rekonstruktionsfehler muss < 1e-10 sein."""
        A = np.array([[3.0, 1.0], [1.0, 3.0]])
        result = spectral_theorem_demo(A)
        assert result['is_valid'] is True
        assert result['reconstruction_error'] == pytest.approx(0.0, abs=1e-8)

    def test_spectral_theorem_demo_eigenvalues_real(self):
        """Alle Eigenwerte selbstadjungierter Matrix sind reell."""
        A = np.array([[4.0, 2.0], [2.0, 1.0]])
        result = spectral_theorem_demo(A)
        for ev in result['eigenvalues']:
            assert isinstance(ev, float)


# ============================================================
# Tests: Fredholm-Alternative
# ============================================================

class TestFredholmAlternative:
    """
    @brief Tests für die Fredholm-Alternative fredholm_alternative().
    @description Prüft Lösbarkeit, Kern und Fredholm-Index.
    """

    def test_fredholm_solvable_full_rank(self):
        """Ax = b ist lösbar wenn A vollen Rang hat und b im Bildraum liegt."""
        A = np.eye(3)
        b = np.array([1.0, 2.0, 3.0])
        result = fredholm_alternative(A, b)
        assert result['solvable'] is True
        assert result['solution'] is not None

    def test_fredholm_solution_correct(self):
        """Lösung von Ix = b muss b selbst sein."""
        A = np.eye(3)
        b = np.array([1.0, 2.0, 3.0])
        result = fredholm_alternative(A, b)
        assert result['solution'] == pytest.approx([1.0, 2.0, 3.0], abs=1e-8)

    def test_fredholm_returns_dict_keys(self):
        """Rückgabe enthält alle erwarteten Schlüssel."""
        A = np.eye(2)
        b = np.array([1.0, 0.0])
        result = fredholm_alternative(A, b)
        for key in ('solvable', 'solution', 'ker_A', 'ker_A_star', 'fredholm_index'):
            assert key in result

    def test_fredholm_index_identity(self):
        """Fredholm-Index der Identität ist 0."""
        A = np.eye(3)
        b = np.zeros(3)
        result = fredholm_alternative(A, b)
        assert result['fredholm_index'] == 0

    def test_fredholm_kernel_singular_matrix(self):
        """Singuläre Matrix hat nicht-trivialen Kern."""
        A = np.array([[1.0, 1.0], [1.0, 1.0]])  # Rang 1
        b = np.array([0.0, 0.0])
        result = fredholm_alternative(A, b)
        assert len(result['ker_A']) >= 1


# ============================================================
# Tests: ContinuousFunctions
# ============================================================

class TestContinuousFunctions:
    """
    @brief Tests für ContinuousFunctions C[a,b].
    @description Prüft Supremumsnorm und Dichtheitsaussagen.
    """

    def test_sup_norm_positive(self):
        """Supremumsnorm ist nicht-negativ."""
        cf = ContinuousFunctions(0.0, 1.0)
        f_values = [0.5, -1.0, 0.3]
        assert cf.sup_norm(f_values) >= 0.0

    def test_sup_norm_value(self):
        """Supremumsnorm = Maximum der Beträge."""
        cf = ContinuousFunctions(0.0, 1.0)
        f_values = [1.0, -3.0, 2.0]
        assert cf.sup_norm(f_values) == pytest.approx(3.0, rel=1e-10)

    def test_sup_norm_zero_function(self):
        """Supremumsnorm der Nullfunktion = 0."""
        cf = ContinuousFunctions(0.0, 1.0)
        assert cf.sup_norm([0.0, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)

    def test_is_dense_in_L2(self):
        """C[a,b] liegt dicht in L²[a,b]."""
        cf = ContinuousFunctions(0.0, 1.0)
        assert cf.is_dense_in_L2() is True

    def test_stone_weierstrass_demo_returns_dict(self):
        """stone_weierstrass_demo() gibt dict mit korrekten Schlüsseln zurück."""
        cf = ContinuousFunctions(0.0, 1.0)
        result = cf.stone_weierstrass_demo()
        for key in ('function', 'degree', 'approximation_error', 'dense', 'statement'):
            assert key in result

    def test_stone_weierstrass_demo_dense_true(self):
        """Stone-Weierstraß: dense == True."""
        cf = ContinuousFunctions(0.0, 1.0)
        result = cf.stone_weierstrass_demo()
        assert result['dense'] is True


# ============================================================
# Tests: SobolevSpace
# ============================================================

class TestSobolevSpace:
    """
    @brief Tests für die Klasse SobolevSpace W^{k,p}(Ω).
    @description Prüft Sobolev-Norm und Einbettungssatz.
    """

    def test_sobolev_norm_nonnegative(self):
        """Sobolev-Norm ist nicht-negativ."""
        sob = SobolevSpace(1, 2, (0.0, 1.0))
        x = np.linspace(0, 1, 50).tolist()
        f = np.sin(np.linspace(0, 1, 50)).tolist()
        norm_val = sob.sobolev_norm(f, x, 1)
        assert norm_val >= 0.0

    def test_sobolev_norm_greater_than_l2(self):
        """W^{1,2}-Norm ≥ L²-Norm (wegen Ableitungsterm)."""
        sob = SobolevSpace(1, 2, (0.0, 1.0))
        x = np.linspace(0, 1, 100).tolist()
        f = np.sin(np.pi * np.linspace(0, 1, 100)).tolist()
        w11_norm = sob.sobolev_norm(f, x, 1)
        w00_norm = sob.sobolev_norm(f, x, 0)
        assert w11_norm >= w00_norm - 1e-10

    def test_embedding_theorem_returns_dict(self):
        """embedding_theorem() gibt dict mit korrekten Schlüsseln zurück."""
        sob = SobolevSpace(2, 2, (0.0, 1.0))
        result = sob.embedding_theorem()
        for key in ('space', 'dimension', 'embeddings', 'continuous_embedding', 'is_hilbert'):
            assert key in result

    def test_embedding_theorem_hilbert_for_p2(self):
        """W^{k,2} ist ein Hilbert-Raum."""
        sob = SobolevSpace(1, 2, (0.0, 1.0))
        result = sob.embedding_theorem()
        assert result['is_hilbert'] is True

    def test_embedding_theorem_not_hilbert_for_p1(self):
        """W^{k,1} ist kein Hilbert-Raum."""
        sob = SobolevSpace(1, 1, (0.0, 1.0))
        result = sob.embedding_theorem()
        assert result['is_hilbert'] is False


# ============================================================
# Tests: L2Space
# ============================================================

class TestL2Space:
    """
    @brief Tests für den L²[a,b]-Raum.
    @description Prüft L²-Skalarprodukt, L²-Norm und Orthogonalität.
    """

    def test_inner_product_nonnegative_same_function(self):
        """⟨f,f⟩ ≥ 0 für alle f."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 200).tolist()
        f = np.sin(np.pi * np.linspace(0, 1, 200)).tolist()
        assert l2.inner_product(f, f, x) >= 0.0

    def test_inner_product_symmetry(self):
        """⟨f,g⟩ = ⟨g,f⟩."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 200).tolist()
        f = np.sin(2 * np.pi * np.linspace(0, 1, 200)).tolist()
        g = np.cos(2 * np.pi * np.linspace(0, 1, 200)).tolist()
        assert l2.inner_product(f, g, x) == pytest.approx(
            l2.inner_product(g, f, x), abs=1e-8
        )

    def test_norm_l2_nonneg(self):
        """L²-Norm ist nicht-negativ."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 100).tolist()
        f = np.ones(100).tolist()
        assert l2.norm(f, x) >= 0.0

    def test_norm_l2_constant_one(self):
        """‖1‖_L²[0,1] = 1."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 1000).tolist()
        f = np.ones(1000).tolist()
        assert l2.norm(f, x) == pytest.approx(1.0, rel=1e-3)

    def test_norm_l2_zero_function(self):
        """‖0‖_L² = 0."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 100).tolist()
        f = np.zeros(100).tolist()
        assert l2.norm(f, x) == pytest.approx(0.0, abs=1e-12)

    def test_fourier_basis_orthonormality(self):
        """Fourier-ONB erfüllt ⟨eᵢ,eⱼ⟩ ≈ δᵢⱼ."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 1000)
        basis_funcs = l2.orthonormal_basis_fourier(3, 0.0, 1.0)
        x_list = x.tolist()
        # Überprüfe Orthonormalität der ersten 3 Basisvektoren
        for i in range(3):
            ei = basis_funcs[i](x).tolist()
            # Norm ≈ 1
            norm_i = l2.norm(ei, x_list)
            assert norm_i == pytest.approx(1.0, rel=1e-2)
            for j in range(i + 1, 3):
                ej = basis_funcs[j](x).tolist()
                ip_ij = l2.inner_product(ei, ej, x_list)
                assert ip_ij == pytest.approx(0.0, abs=1e-2)

    def test_projection_onto_polynomials(self):
        """Projektion auf Polynome liefert Liste gleicher Länge wie x."""
        l2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 50).tolist()
        f = np.sin(np.pi * np.linspace(0, 1, 50)).tolist()
        proj = l2.projection_onto_polynomials(f, x, degree=3)
        assert len(proj) == 50


# ============================================================
# Tests: Klassische Sätze
# ============================================================

class TestClassicalTheorems:
    """
    @brief Tests für die klassischen Sätze der Funktionalanalysis.
    @description Hahn-Banach, offene Abbildung, abgeschlossener Graph,
                 gleichmäßige Beschränktheit, schwache Konvergenz,
                 Riesz-Darstellung, kompakte Operatoren.
    """

    def test_hahn_banach_demo_returns_dict(self):
        """hahn_banach_theorem_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = hahn_banach_theorem_demo()
        for key in ('theorem', 'norm_preserved', 'original_norm',
                    'extension_norm', 'consequences'):
            assert key in result

    def test_hahn_banach_norm_preserved(self):
        """Hahn-Banach: Fortsetzung bewahrt die Norm."""
        result = hahn_banach_theorem_demo()
        assert result['norm_preserved'] is True

    def test_hahn_banach_norm_values(self):
        """Hahn-Banach: Original- und Fortsetzungsnorm = 2.0."""
        result = hahn_banach_theorem_demo()
        assert result['original_norm'] == pytest.approx(2.0, rel=1e-10)
        assert result['extension_norm'] == pytest.approx(2.0, rel=1e-10)

    def test_open_mapping_demo_returns_dict(self):
        """open_mapping_theorem_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = open_mapping_theorem_demo()
        for key in ('theorem', 'is_surjective', 'det', 'inverse_exists'):
            assert key in result

    def test_open_mapping_surjective(self):
        """Operator im Demo-Beispiel ist surjektiv (det ≠ 0)."""
        result = open_mapping_theorem_demo()
        assert result['is_surjective'] is True

    def test_open_mapping_inverse_exists(self):
        """Bijektiver Operator hat Inverse."""
        result = open_mapping_theorem_demo()
        assert result['inverse_exists'] is True

    def test_closed_graph_demo_returns_dict(self):
        """closed_graph_theorem_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = closed_graph_theorem_demo()
        for key in ('theorem', 'graph_closed', 'operator_bounded', 'convergence_errors'):
            assert key in result

    def test_closed_graph_closed(self):
        """Im Demo-Beispiel ist der Graph abgeschlossen."""
        result = closed_graph_theorem_demo()
        assert result['graph_closed'] is True

    def test_closed_graph_operator_bounded(self):
        """Endlich-dimensionaler Operator ist immer beschränkt."""
        result = closed_graph_theorem_demo()
        assert result['operator_bounded'] is True

    def test_uniform_boundedness_demo_returns_dict(self):
        """uniform_boundedness_principle_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = uniform_boundedness_principle_demo()
        for key in ('theorem', 'pointwise_bounded', 'uniformly_bounded',
                    'max_op_norm'):
            assert key in result

    def test_uniform_boundedness_both_bounded(self):
        """Banach-Steinhaus: punktweise beschränkt → gleichmäßig beschränkt."""
        result = uniform_boundedness_principle_demo()
        assert result['pointwise_bounded'] is True
        assert result['uniformly_bounded'] is True

    def test_weak_convergence_demo_returns_dict(self):
        """weak_convergence_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = weak_convergence_demo()
        for key in ('weak_values', 'strong_norms', 'weakly_converges_to_zero',
                    'strongly_converges_to_zero', 'conclusion'):
            assert key in result

    def test_weak_convergence_weak_decreasing(self):
        """Schwache Werte ⟨eₙ,y⟩ nehmen ab (yₙ → 0)."""
        result = weak_convergence_demo()
        assert result['weakly_converges_to_zero'] is True

    def test_weak_convergence_strong_norms_constant(self):
        """Starke Normen der Einheitsvektoren sind konstant 1."""
        result = weak_convergence_demo()
        for n in result['strong_norms']:
            assert n == pytest.approx(1.0, rel=1e-10)

    def test_riesz_representation_returns_dict(self):
        """riesz_representation_theorem() gibt dict mit korrekten Schlüsseln zurück."""
        result = riesz_representation_theorem()
        for key in ('theorem', 'riesz_representative', 'functional',
                    'norm_functional', 'verifications'):
            assert key in result

    def test_riesz_representation_norm(self):
        """‖f‖ = ‖y‖ = ‖[2,3,-1]‖ = √14."""
        result = riesz_representation_theorem()
        assert result['norm_functional'] == pytest.approx(math.sqrt(14.0), rel=1e-8)

    def test_riesz_representation_verifications(self):
        """Alle Verifikationen müssen f(x) = ⟨x,y⟩ bestätigen."""
        result = riesz_representation_theorem()
        for v in result['verifications']:
            assert v['equal'] is True

    def test_compact_operator_theory_returns_dict(self):
        """compact_operator_theory() gibt dict mit korrekten Schlüsseln zurück."""
        result = compact_operator_theory()
        for key in ('hilbert_schmidt_norm', 'is_compact', 'eigenvalues',
                    'spectral_radius', 'accumulation_point'):
            assert key in result

    def test_compact_operator_theory_is_compact(self):
        """Kompakter Operator: is_compact == True."""
        result = compact_operator_theory()
        assert result['is_compact'] is True

    def test_compact_operator_accumulation_point(self):
        """Einziger Häufungspunkt des Spektrums kompakter Operatoren ist 0."""
        result = compact_operator_theory()
        assert result['accumulation_point'] == pytest.approx(0.0, abs=1e-10)


# ============================================================
# Tests: Fixpunktsätze
# ============================================================

class TestFixedPointTheorems:
    """
    @brief Tests für die Fixpunktsätze (Banach, Schauder).
    @description Prüft Konvergenz, Fixpunkte und Rückgabeformat.
    """

    def test_banach_fixed_point_converges(self):
        """Banach-Fixpunktsatz: Kontraktion T(x) = x/2 konvergiert zu 0."""
        result = banach_fixed_point(lambda x: x / 2.0, x0=1.0, tol=1e-10, max_iter=1000)
        assert result['converged'] is True
        assert result['fixed_point'] == pytest.approx(0.0, abs=1e-8)

    def test_banach_fixed_point_cosine(self):
        """T(x) = cos(x) hat Fixpunkt bei x ≈ 0.7391 (Dottie-Zahl)."""
        result = banach_fixed_point(math.cos, x0=0.5, tol=1e-10, max_iter=1000)
        assert result['converged'] is True
        assert result['fixed_point'] == pytest.approx(0.7390851332, rel=1e-6)

    def test_banach_fixed_point_residual_small(self):
        """Nach Konvergenz: |T(x*) - x*| < Toleranz."""
        result = banach_fixed_point(lambda x: x / 3.0, x0=10.0, tol=1e-10)
        assert result['final_residual'] == pytest.approx(0.0, abs=1e-8)

    def test_banach_fixed_point_returns_history(self):
        """Rückgabe enthält Iterationshistorie."""
        result = banach_fixed_point(lambda x: x / 2.0, x0=1.0)
        assert 'history' in result
        assert len(result['history']) >= 1

    def test_banach_fixed_point_contraction_constant(self):
        """Kontraktionskonstante für T(x) = x/2 liegt in [0,1)."""
        result = banach_fixed_point(lambda x: x / 2.0, x0=1.0, tol=1e-12)
        q = result['contraction_constant_estimate']
        if q is not None:
            assert 0.0 <= q < 1.0 + 1e-6  # Numerische Toleranz

    def test_banach_fixed_point_returns_dict_keys(self):
        """Rückgabe enthält alle erwarteten Schlüssel."""
        result = banach_fixed_point(lambda x: x / 4.0, x0=8.0)
        for key in ('fixed_point', 'converged', 'iterations',
                    'history', 'final_residual', 'contraction_constant_estimate'):
            assert key in result

    def test_schauder_fixed_point_demo_returns_dict(self):
        """schauder_fixed_point_demo() gibt dict mit korrekten Schlüsseln zurück."""
        result = schauder_fixed_point_demo()
        for key in ('theorem', 'K_convex', 'K_compact',
                    'T_continuous', 'fixed_points_found'):
            assert key in result

    def test_schauder_fixed_point_demo_k_properties(self):
        """K = [0,1] ist konvex und kompakt."""
        result = schauder_fixed_point_demo()
        assert result['K_convex'] is True
        assert result['K_compact'] is True

    def test_schauder_fixed_point_demo_finds_fixed_points(self):
        """T(x) = x² hat Fixpunkte bei 0 und 1."""
        result = schauder_fixed_point_demo()
        found = sorted(result['fixed_points_found'])
        # Mindestens 0 oder 1 muss gefunden worden sein
        assert len(found) >= 1
        # Prüfe ob 0.0 oder 1.0 in den gefundenen Fixpunkten vorkommt
        has_zero = any(abs(fp) < 1e-6 for fp in found)
        has_one = any(abs(fp - 1.0) < 1e-6 for fp in found)
        assert has_zero or has_one

    def test_banach_no_convergence_non_contraction(self):
        """T(x) = 2x ist keine Kontraktion — konvergiert nicht für x₀ ≠ 0."""
        result = banach_fixed_point(lambda x: 2.0 * x, x0=1.0,
                                    tol=1e-10, max_iter=50)
        # Darf nicht konvergiert sein (oder Fixpunkt nur wenn x₀=0)
        assert result['converged'] is False
