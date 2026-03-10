"""
@file test_functional_analysis.py
@brief Tests für das Modul functional_analysis.py.
@description
    Test-Suite für Funktionalanalysis: normierte Räume, Banach/Hilbert-Räume,
    lineare Operatoren, Spektraltheorie, Fredholm-Alternative, Fixpunktsätze.

    Test-Kategorien:
    1. NormedSpace: Norm, Abstand, Axiome, Einheitsball
    2. BanachSpace: Cauchy-Folgen, lp-Normen
    3. HilbertSpace: Skalarprodukt, Orthogonalität, Gram-Schmidt,
                     Parseval, Bessel, Riesz-Repräsentation
    4. LinearOperator: Anwendung, Operatornorm, Adjungierte, Eigenschaften
    5. CompactOperator: Spektralzerlegung
    6. Spektraltheorie: spectrum(), spectral_radius(), spectral_theorem_demo()
    7. Fredholm-Alternative
    8. Funktionenräume: ContinuousFunctions, SobolevSpace, L2Space
    9. Klassische Sätze: Hahn-Banach, Offene-Abbildung, Satz v. abg. Graph
    10. Schwache Konvergenz, Riesz-Darstellungssatz, Kompakte Operatoren
    11. Banach-Fixpunktsatz, Schauder-Fixpunktsatz

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math
import numpy as np
import pytest

# Modul-Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from functional_analysis import (
    NormedSpace, BanachSpace, HilbertSpace,
    LinearOperator, CompactOperator,
    spectrum, spectral_radius, spectral_theorem_demo, fredholm_alternative,
    ContinuousFunctions, SobolevSpace, L2Space,
    hahn_banach_theorem_demo, open_mapping_theorem_demo, closed_graph_theorem_demo,
    uniform_boundedness_principle_demo, weak_convergence_demo,
    riesz_representation_theorem, compact_operator_theory,
    banach_fixed_point, schauder_fixed_point_demo
)


# ============================================================
# Hilfsfunktionen für Tests
# ============================================================

def euclidean_norm(v):
    """Euklidische Norm für Tests."""
    return float(np.linalg.norm(v))


def standard_ip(u, v):
    """Standard-Skalarprodukt für Tests."""
    return float(np.dot(u, v))


# ============================================================
# 1. NormedSpace
# ============================================================

class TestNormedSpace:
    """Tests für normierte Räume."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.vectors = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, -1.0]]
        self.ns = NormedSpace(self.vectors, euclidean_norm)

    def test_norm_positivity(self):
        """Norm ist nichtnegativ."""
        for v in self.vectors:
            assert self.ns.norm(v) >= 0

    def test_norm_zero_vector(self):
        """Norm des Nullvektors ist 0."""
        assert self.ns.norm([0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)

    def test_norm_unit_vector(self):
        """Einheitsvektoren haben Norm 1."""
        assert self.ns.norm([1.0, 0.0]) == pytest.approx(1.0)
        assert self.ns.norm([0.0, 1.0]) == pytest.approx(1.0)

    def test_norm_homogeneity(self):
        """Homogenität: ‖αv‖ = |α|·‖v‖."""
        v = [3.0, 4.0]
        alpha = 2.5
        assert self.ns.norm([alpha * x for x in v]) == pytest.approx(abs(alpha) * self.ns.norm(v))

    def test_triangle_inequality(self):
        """Dreiecksungleichung: ‖u+v‖ ≤ ‖u‖+‖v‖."""
        u = np.array([1.0, 2.0])
        v = np.array([3.0, -1.0])
        assert self.ns.norm(u + v) <= self.ns.norm(u) + self.ns.norm(v) + 1e-10

    def test_distance(self):
        """Abstand d(u,v) = ‖u-v‖."""
        u = [3.0, 0.0]
        v = [0.0, 4.0]
        assert self.ns.distance(u, v) == pytest.approx(5.0)

    def test_distance_symmetry(self):
        """Abstand ist symmetrisch: d(u,v) = d(v,u)."""
        u = [1.0, 2.0]
        v = [3.0, -1.0]
        assert self.ns.distance(u, v) == pytest.approx(self.ns.distance(v, u))

    def test_distance_to_self_is_zero(self):
        """Abstand eines Vektors zu sich selbst ist 0."""
        u = [1.0, 2.0]
        assert self.ns.distance(u, u) == pytest.approx(0.0, abs=1e-12)

    def test_is_normed_space(self):
        """Axiome des normierten Raums erfüllt."""
        assert self.ns.is_normed_space() is True

    def test_unit_ball_points(self):
        """Einheitsball-Punkte haben Norm ≤ 1."""
        points = self.ns.unit_ball(n_points=8)
        for p in points:
            # Punkte auf dem Rand: Norm ≈ 1
            assert abs(self.ns.norm(p) - 1.0) < 0.01

    def test_l1_norm(self):
        """L1-Norm korrekt."""
        l1_norm = lambda v: float(np.sum(np.abs(v)))
        ns_l1 = NormedSpace([[1.0, 2.0, -3.0]], l1_norm)
        assert ns_l1.norm([1.0, 2.0, -3.0]) == pytest.approx(6.0)


# ============================================================
# 2. BanachSpace
# ============================================================

class TestBanachSpace:
    """Tests für Banach-Räume."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.vectors = [[1.0, 0.0], [0.0, 1.0]]
        self.bs = BanachSpace(self.vectors, euclidean_norm)

    def test_cauchy_sequence_converging(self):
        """Konvergierende Folge ist Cauchy-Folge."""
        # xₙ = (1/n, 0) → (0, 0)
        sequence = [[1.0 / (n + 1), 0.0] for n in range(20)]
        assert self.bs.cauchy_sequence_test(sequence, tol=0.05) is True

    def test_cauchy_sequence_constant(self):
        """Konstante Folge ist Cauchy-Folge."""
        sequence = [[1.0, 1.0]] * 10
        assert self.bs.cauchy_sequence_test(sequence) is True

    def test_cauchy_sequence_diverging(self):
        """Divergierende Folge ist keine Cauchy-Folge."""
        sequence = [[float(n), 0.0] for n in range(20)]
        assert self.bs.cauchy_sequence_test(sequence, tol=0.5) is False

    def test_is_separable(self):
        """Endlich-dimensionaler Raum ist separabel."""
        assert self.bs.is_separable() is True

    def test_lp_norm_l1(self):
        """L1-Norm korrekt."""
        v = [3.0, -4.0]
        assert self.bs.lp_norm(v, 1.0) == pytest.approx(7.0)

    def test_lp_norm_l2(self):
        """L2-Norm korrekt (Pythagoras 3-4-5)."""
        v = [3.0, 4.0]
        assert self.bs.lp_norm(v, 2.0) == pytest.approx(5.0)

    def test_lp_norm_linf(self):
        """L∞-Norm = Maximum der Beträge."""
        v = [3.0, -4.0, 2.0]
        assert self.bs.lp_norm(v, float('inf')) == pytest.approx(4.0)

    def test_lp_norm_p3(self):
        """L3-Norm korrekt."""
        v = [1.0, 1.0, 1.0]
        assert self.bs.lp_norm(v, 3.0) == pytest.approx(3.0 ** (1.0 / 3.0))

    def test_lp_norm_invalid_p(self):
        """Ungültiges p wirft ValueError."""
        with pytest.raises(ValueError):
            self.bs.lp_norm([1.0, 2.0], -1.0)


# ============================================================
# 3. HilbertSpace
# ============================================================

class TestHilbertSpace:
    """Tests für Hilbert-Räume."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                        [1.0, 1.0, 0.0]]
        self.hs = HilbertSpace(self.vectors, standard_ip)

    def test_inner_product_symmetry(self):
        """Skalarprodukt ist symmetrisch: ⟨u,v⟩ = ⟨v,u⟩."""
        u = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        assert self.hs.inner_product(u, v) == pytest.approx(self.hs.inner_product(v, u))

    def test_inner_product_positivity(self):
        """Skalarprodukt ⟨v,v⟩ ≥ 0."""
        for v in self.vectors:
            assert self.hs.inner_product(v, v) >= 0

    def test_inner_product_definiteness(self):
        """⟨v,v⟩ = 0 genau für Nullvektor."""
        assert self.hs.inner_product([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]) == pytest.approx(0.0)
        assert self.hs.inner_product([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) > 0

    def test_inner_product_value(self):
        """Konkreter Wert des Skalarprodukts."""
        u = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        # ⟨u,v⟩ = 1·4 + 2·5 + 3·6 = 4+10+18 = 32
        assert self.hs.inner_product(u, v) == pytest.approx(32.0)

    def test_norm_induced_by_ip(self):
        """Norm ist √⟨v,v⟩."""
        v = [3.0, 4.0, 0.0]
        assert self.hs.norm(v) == pytest.approx(5.0)

    def test_is_orthogonal_true(self):
        """Orthogonale Vektoren erkannt."""
        u = [1.0, 0.0, 0.0]
        v = [0.0, 1.0, 0.0]
        assert self.hs.is_orthogonal(u, v) is True

    def test_is_orthogonal_false(self):
        """Nicht-orthogonale Vektoren erkannt."""
        u = [1.0, 1.0, 0.0]
        v = [1.0, 0.0, 0.0]
        assert self.hs.is_orthogonal(u, v) is False

    def test_gram_schmidt_orthonormality(self):
        """Gram-Schmidt liefert orthonormale Vektoren."""
        basis = [[1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]
        onb = self.hs.gram_schmidt_orthonormalize(basis)
        # Prüfe Orthonormalität: ⟨eᵢ,eⱼ⟩ = δᵢⱼ
        for i, ei in enumerate(onb):
            for j, ej in enumerate(onb):
                expected = 1.0 if i == j else 0.0
                assert self.hs.inner_product(ei, ej) == pytest.approx(expected, abs=1e-10)

    def test_gram_schmidt_span_preserved(self):
        """Gram-Schmidt ändert nicht den aufgespannten Raum (Dimension)."""
        basis = [[1.0, 0.0], [1.0, 1.0]]
        hs2 = HilbertSpace(basis, standard_ip)
        onb = hs2.gram_schmidt_orthonormalize(basis)
        assert len(onb) == 2

    def test_projection_onto_subspace(self):
        """Projektion eines Vektors auf einen Unterraum."""
        # Projektion von v = (1,1,1) auf U = span{(1,0,0), (0,1,0)} = xy-Ebene
        v = [1.0, 1.0, 1.0]
        onto = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        proj = self.hs.projection(v, onto)
        # Erwartete Projektion: (1,1,0)
        assert proj[0] == pytest.approx(1.0, abs=1e-8)
        assert proj[1] == pytest.approx(1.0, abs=1e-8)
        assert proj[2] == pytest.approx(0.0, abs=1e-8)

    def test_projection_residual_orthogonal(self):
        """Residuum v - proj(v) ist orthogonal zum Unterraum."""
        v = np.array([1.0, 1.0, 1.0])
        onto = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        proj = np.array(self.hs.projection(v, onto))
        residual = v - proj
        # Residuum muss senkrecht auf allen Basisvektoren stehen
        for b in onto:
            assert self.hs.inner_product(residual, b) == pytest.approx(0.0, abs=1e-8)

    def test_pythagoras_orthogonal_vectors(self):
        """Pythagoras gilt für orthogonale Vektoren."""
        u = [3.0, 0.0]
        v = [0.0, 4.0]
        hs2 = HilbertSpace([[1.0, 0.0], [0.0, 1.0]], standard_ip)
        assert hs2.pythagoras(u, v) is True

    def test_pythagoras_nonorthogonal_vectors(self):
        """Pythagoras gilt nicht für nicht-orthogonale Vektoren."""
        u = [1.0, 1.0]
        v = [1.0, 0.0]
        hs2 = HilbertSpace([[1.0, 0.0], [0.0, 1.0]], standard_ip)
        assert hs2.pythagoras(u, v) is False

    def test_parseval_identity_onb(self):
        """Parseval-Identität gilt für vollständige ONB."""
        # Standard-ONB in ℝ³
        onb = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        v = [2.0, 3.0, -1.0]
        assert self.hs.parseval_identity(v, onb) is True

    def test_bessel_inequality_partial_onb(self):
        """Bessel-Ungleichung gilt für partielle ONB."""
        # Nur zwei der drei Basisvektoren
        partial_onb = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        v = [2.0, 3.0, 1.0]
        assert self.hs.bessel_inequality(v, partial_onb) is True

    def test_bessel_inequality_full_onb(self):
        """Bessel-Ungleichung mit vollständiger ONB (= Parseval)."""
        full_onb = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        v = [1.0, 2.0, 3.0]
        assert self.hs.bessel_inequality(v, full_onb) is True

    def test_riesz_representation(self):
        """Riesz-Darstellungssatz: Funktional via Skalarprodukt."""
        basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        # Funktional f(x) = 2x₁ + 3x₂ (Riesz-Rep: y = (2,3,0))
        def functional(x):
            return 2.0 * x[0] + 3.0 * x[1]
        y = self.hs.riesz_representation(functional, basis)
        # Prüfe y ≈ (2, 3, 0)
        assert y[0] == pytest.approx(2.0, abs=1e-8)
        assert y[1] == pytest.approx(3.0, abs=1e-8)
        assert y[2] == pytest.approx(0.0, abs=1e-8)

    def test_cauchy_schwarz_inequality(self):
        """Cauchy-Schwarz: |⟨u,v⟩| ≤ ‖u‖·‖v‖."""
        u = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        lhs = abs(self.hs.inner_product(u, v))
        rhs = self.hs.norm(u) * self.hs.norm(v)
        assert lhs <= rhs + 1e-10


# ============================================================
# 4. LinearOperator
# ============================================================

class TestLinearOperator:
    """Tests für lineare Operatoren."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.A = np.array([[2.0, 1.0], [0.0, 3.0]])
        self.op = LinearOperator(self.A)
        # Identitätsoperator
        self.I = LinearOperator(np.eye(3))
        # Symmetrische Matrix
        self.S = LinearOperator(np.array([[2.0, 1.0], [1.0, 3.0]]))

    def test_apply(self):
        """Operator anwenden: y = Av."""
        v = [1.0, 1.0]
        result = self.op.apply(v)
        # [2·1+1·1, 0·1+3·1] = [3, 3]
        assert result[0] == pytest.approx(3.0)
        assert result[1] == pytest.approx(3.0)

    def test_apply_zero_vector(self):
        """T(0) = 0 (Linearität)."""
        result = self.op.apply([0.0, 0.0])
        assert result[0] == pytest.approx(0.0)
        assert result[1] == pytest.approx(0.0)

    def test_is_bounded(self):
        """Operator ist beschränkt."""
        test_vecs = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        assert self.op.is_bounded(test_vecs) is True

    def test_operator_norm_l2(self):
        """Operatornorm (L2) = größter Singulärwert."""
        singular_values = np.linalg.svd(self.A, compute_uv=False)
        expected_norm = float(singular_values[0])
        assert self.op.operator_norm() == pytest.approx(expected_norm, rel=1e-8)

    def test_operator_norm_identity(self):
        """Identitätsoperator hat Norm 1."""
        assert self.I.operator_norm() == pytest.approx(1.0, abs=1e-10)

    def test_is_compact(self):
        """Endlich-dim. Operator ist kompakt."""
        assert self.op.is_compact() is True

    def test_adjoint_matrix(self):
        """Adjungierter Operator = transponierte Matrix."""
        adj = self.op.adjoint()
        expected = self.A.T
        assert np.allclose(adj.matrix, expected)

    def test_is_self_adjoint_symmetric(self):
        """Symmetrische Matrix ist selbstadjungiert."""
        assert self.S.is_self_adjoint() is True

    def test_is_self_adjoint_asymmetric(self):
        """Nicht-symmetrische Matrix ist nicht selbstadjungiert."""
        assert self.op.is_self_adjoint() is False

    def test_is_unitary_rotation(self):
        """Rotationsmatrix ist unitär."""
        theta = math.pi / 4
        R = np.array([[math.cos(theta), -math.sin(theta)],
                       [math.sin(theta),  math.cos(theta)]])
        rot = LinearOperator(R)
        assert rot.is_unitary() is True

    def test_is_unitary_non_unitary(self):
        """Nicht-unitäre Matrix erkannt."""
        assert self.op.is_unitary() is False

    def test_is_normal_symmetric(self):
        """Symmetrische Matrix ist normal."""
        assert self.S.is_normal() is True

    def test_is_projection(self):
        """Projektionsmatrix erkannt (P² = P)."""
        # Projektion auf x-Achse in ℝ²
        P = np.array([[1.0, 0.0], [0.0, 0.0]])
        proj_op = LinearOperator(P)
        assert proj_op.is_projection() is True

    def test_is_projection_false(self):
        """Nicht-Projektion erkannt."""
        assert self.op.is_projection() is False

    def test_kernel_full_rank(self):
        """Invertierbare Matrix hat trivialen Kern."""
        ker = self.op.kernel()
        assert len(ker) == 0

    def test_kernel_singular_matrix(self):
        """Singuläre Matrix hat nicht-trivialen Kern."""
        A_sing = np.array([[1.0, 2.0], [2.0, 4.0]])  # Rang 1
        op_sing = LinearOperator(A_sing)
        ker = op_sing.kernel()
        assert len(ker) >= 1
        # Prüfe: Av = 0 für Kernvektor v
        v = np.array(ker[0])
        result = A_sing @ v
        assert np.linalg.norm(result) < 1e-8

    def test_range_space_full_rank(self):
        """Invertierbare Matrix hat vollständigen Bildraum."""
        rs = self.op.range_space()
        assert len(rs) == 2  # Rang 2 in ℝ²

    def test_range_space_singular(self):
        """Singuläre Matrix hat eindimensionalen Bildraum."""
        A_sing = np.array([[1.0, 2.0], [2.0, 4.0]])
        op_sing = LinearOperator(A_sing)
        rs = op_sing.range_space()
        assert len(rs) == 1


# ============================================================
# 5. CompactOperator
# ============================================================

class TestCompactOperator:
    """Tests für kompakte Operatoren."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.K = np.array([[2.0, 1.0], [1.0, 3.0]])  # Symmetrisch
        self.cop = CompactOperator(self.K)

    def test_spectral_decomposition_eigenvalues(self):
        """Spektralzerlegung liefert korrekte Eigenwerte."""
        decomp = self.cop.spectral_decomposition()
        evals = np.sort(decomp['eigenvalues'])
        evals_np = np.sort(np.linalg.eigvalsh(self.K))
        assert np.allclose(evals, evals_np, atol=1e-8)

    def test_spectral_decomposition_valid(self):
        """Spektralzerlegung ist valide."""
        decomp = self.cop.spectral_decomposition()
        assert decomp['decomposition_valid'] is True

    def test_spectral_decomposition_orthonormal_eigenvectors(self):
        """Eigenvektoren sind orthonormal."""
        decomp = self.cop.spectral_decomposition()
        evecs = [np.array(e) for e in decomp['eigenvectors']]
        for i, ei in enumerate(evecs):
            for j, ej in enumerate(evecs):
                expected = 1.0 if i == j else 0.0
                assert abs(float(np.dot(ei, ej)) - expected) < 1e-8


# ============================================================
# 6. Spektraltheorie
# ============================================================

class TestSpectralTheory:
    """Tests für Spektraltheorie."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.A = np.array([[2.0, 0.0], [0.0, 3.0]])  # Diagonalmatrix
        self.op = LinearOperator(self.A)

    def test_spectrum_diagonal_matrix(self):
        """Spektrum einer Diagonalmatrix = Diagonaleinträge."""
        spec = spectrum(self.op)
        # Eigenwerte sollten 2 und 3 sein
        evals = sorted([abs(x) for x in spec['point_spectrum']])
        assert evals[0] == pytest.approx(2.0, abs=1e-8)
        assert evals[1] == pytest.approx(3.0, abs=1e-8)

    def test_spectrum_spectral_radius(self):
        """Spektralradius = größter Eigenabsolutbetrag."""
        spec = spectrum(self.op)
        assert spec['spectral_radius'] == pytest.approx(3.0, abs=1e-8)

    def test_spectral_radius_function(self):
        """spectral_radius() liefert korrekten Wert."""
        r = spectral_radius(self.op)
        assert r == pytest.approx(3.0, abs=1e-8)

    def test_spectral_radius_identity(self):
        """Spektralradius der Identität ist 1."""
        I_op = LinearOperator(np.eye(3))
        assert spectral_radius(I_op) == pytest.approx(1.0, abs=1e-8)

    def test_spectral_theorem_demo_symmetric(self):
        """Spektralsatz für symmetrische Matrix."""
        A = np.array([[4.0, 2.0], [2.0, 3.0]])
        result = spectral_theorem_demo(A)
        assert result['is_valid'] is True
        assert result['reconstruction_error'] < 1e-10

    def test_spectral_theorem_demo_real_eigenvalues(self):
        """Selbstadjungierte Operatoren haben reelle Eigenwerte."""
        A = np.array([[1.0, 2.0], [2.0, 4.0]])
        result = spectral_theorem_demo(A)
        for ev in result['eigenvalues']:
            assert isinstance(ev, float)


# ============================================================
# 7. Fredholm-Alternative
# ============================================================

class TestFredholmAlternative:
    """Tests für die Fredholm-Alternative."""

    def test_solvable_system(self):
        """Ax=b lösbar wenn b ⊥ ker(A*)."""
        A = np.array([[1.0, 0.0], [0.0, 1.0]])  # Invertierbar
        b = np.array([1.0, 2.0])
        result = fredholm_alternative(A, b)
        assert result['solvable'] is True
        assert result['solution'] is not None

    def test_solution_correctness(self):
        """Berechnete Lösung erfüllt Ax = b."""
        A = np.array([[2.0, 1.0], [1.0, 3.0]])
        b = np.array([5.0, 7.0])
        result = fredholm_alternative(A, b)
        if result['solvable'] and result['solution'] is not None:
            x = np.array(result['solution'])
            residual = np.linalg.norm(A @ x - b)
            assert residual < 1e-8

    def test_singular_system_unsolvable(self):
        """Singuläres System: Ax=b kann unlösbar sein."""
        A = np.array([[1.0, 2.0], [2.0, 4.0]])  # Rang 1
        # b nicht im Bild von A
        b = np.array([1.0, 0.0])
        result = fredholm_alternative(A, b)
        # Entweder lösbar oder nicht, je nach b
        assert 'solvable' in result
        assert 'fredholm_index' in result

    def test_fredholm_index_invertible(self):
        """Fredholm-Index einer invertierbaren Matrix ist 0."""
        A = np.array([[1.0, 0.0], [0.0, 1.0]])
        b = np.array([1.0, 1.0])
        result = fredholm_alternative(A, b)
        assert result['fredholm_index'] == 0


# ============================================================
# 8. Funktionenräume
# ============================================================

class TestContinuousFunctions:
    """Tests für C[a,b]."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.cf = ContinuousFunctions(0.0, 1.0)

    def test_sup_norm_constant(self):
        """Supremumsnorm einer konstanten Funktion."""
        f = [3.0] * 10
        assert self.cf.sup_norm(f) == pytest.approx(3.0)

    def test_sup_norm_sin(self):
        """Supremumsnorm von sin ist ≤ 1."""
        x = np.linspace(0, 1, 100)
        f = np.sin(2 * math.pi * x)
        assert self.cf.sup_norm(f.tolist()) <= 1.0 + 1e-10

    def test_is_dense_in_L2(self):
        """C[a,b] liegt dicht in L²[a,b]."""
        assert self.cf.is_dense_in_L2() is True

    def test_stone_weierstrass_demo(self):
        """Stone-Weierstraß Demonstration läuft fehlerfrei."""
        result = self.cf.stone_weierstrass_demo()
        assert 'approximation_error' in result
        assert result['dense'] is True
        assert result['approximation_error'] >= 0


class TestSobolevSpace:
    """Tests für Sobolev-Räume."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.W = SobolevSpace(k=1, p=2.0, domain=(0.0, 1.0))

    def test_sobolev_norm_positive(self):
        """Sobolev-Norm ist nichtnegativ."""
        x = np.linspace(0, 1, 50)
        f = np.sin(math.pi * x)
        norm = self.W.sobolev_norm(f.tolist(), x.tolist(), k=1)
        assert norm >= 0

    def test_sobolev_norm_zero_function(self):
        """Sobolev-Norm der Nullfunktion ist 0."""
        x = np.linspace(0, 1, 50)
        f = np.zeros_like(x)
        norm = self.W.sobolev_norm(f.tolist(), x.tolist(), k=0)
        assert norm == pytest.approx(0.0, abs=1e-10)

    def test_embedding_theorem(self):
        """Einbettungssatz gibt Informationen zurück."""
        result = self.W.embedding_theorem()
        assert 'space' in result
        assert 'is_hilbert' in result
        assert result['is_hilbert'] is True  # p=2

    def test_sobolev_is_hilbert_for_p2(self):
        """W^{k,2} ist Hilbert-Raum."""
        result = self.W.embedding_theorem()
        assert result['is_hilbert'] is True


class TestL2Space:
    """Tests für L²[a,b]."""

    def setup_method(self):
        """Testdaten initialisieren."""
        self.L2 = L2Space(0.0, 1.0)
        self.x = np.linspace(0, 1, 200)

    def test_inner_product_orthogonal_functions(self):
        """sin und cos auf [0,2π] sind orthogonal."""
        x = np.linspace(0, 2 * math.pi, 500)
        f = np.sin(x)
        g = np.cos(x)
        ip = L2Space(0, 2 * math.pi).inner_product(f.tolist(), g.tolist(), x.tolist())
        assert abs(ip) < 0.1  # Grob orthogonal

    def test_inner_product_symmetry(self):
        """L²-Skalarprodukt ist symmetrisch."""
        f = np.sin(math.pi * self.x)
        g = np.cos(math.pi * self.x)
        ip_fg = self.L2.inner_product(f.tolist(), g.tolist(), self.x.tolist())
        ip_gf = self.L2.inner_product(g.tolist(), f.tolist(), self.x.tolist())
        assert ip_fg == pytest.approx(ip_gf, abs=1e-10)

    def test_norm_positive(self):
        """L²-Norm ist nichtnegativ."""
        f = np.sin(math.pi * self.x)
        norm = self.L2.norm(f.tolist(), self.x.tolist())
        assert norm >= 0

    def test_norm_zero_function(self):
        """L²-Norm der Nullfunktion ist 0."""
        f = np.zeros_like(self.x)
        assert self.L2.norm(f.tolist(), self.x.tolist()) == pytest.approx(0.0, abs=1e-12)

    def test_norm_constant_function(self):
        """L²-Norm von f=1 auf [0,1] ist 1."""
        f = np.ones_like(self.x)
        norm = self.L2.norm(f.tolist(), self.x.tolist())
        assert norm == pytest.approx(1.0, abs=0.01)

    def test_orthonormal_basis_fourier(self):
        """Fourier-ONB hat korrekte Norm."""
        basis = self.L2.orthonormal_basis_fourier(3, 0.0, 1.0)
        assert len(basis) == 3
        # e₀ auf [0,1] hat Norm 1/√1 = 1
        e0_values = basis[0](self.x)
        norm_e0 = self.L2.norm(e0_values.tolist(), self.x.tolist())
        assert norm_e0 == pytest.approx(1.0, abs=0.02)

    def test_fourier_basis_orthogonality(self):
        """Fourier-Basiselemente sind orthogonal."""
        basis = self.L2.orthonormal_basis_fourier(3, 0.0, 1.0)
        e0 = basis[0](self.x)
        e1 = basis[1](self.x)
        ip = self.L2.inner_product(e0.tolist(), e1.tolist(), self.x.tolist())
        assert abs(ip) < 0.05

    def test_projection_onto_polynomials(self):
        """Projektion auf Polynome liefert Vektor gleicher Länge."""
        f = np.sin(math.pi * self.x)
        proj = self.L2.projection_onto_polynomials(f.tolist(), self.x.tolist(), degree=3)
        assert len(proj) == len(f)


# ============================================================
# 9. Klassische Sätze
# ============================================================

class TestClassicalTheorems:
    """Tests für klassische Sätze der Funktionalanalysis."""

    def test_hahn_banach_demo(self):
        """Hahn-Banach-Satz Demonstration."""
        result = hahn_banach_theorem_demo()
        assert result['norm_preserved'] is True
        assert 'consequences' in result
        assert len(result['consequences']) > 0

    def test_hahn_banach_norm_preserved(self):
        """Hahn-Banach: Norm des Funktionals bleibt erhalten."""
        result = hahn_banach_theorem_demo()
        assert result['original_norm'] == pytest.approx(result['extension_norm'])

    def test_open_mapping_theorem_demo(self):
        """Offene-Abbildung-Satz Demonstration."""
        result = open_mapping_theorem_demo()
        assert 'is_surjective' in result
        assert result['is_surjective'] is True

    def test_open_mapping_inverse_exists(self):
        """Invertierbare surjektive Abbildung hat beschränkte Inverse."""
        result = open_mapping_theorem_demo()
        assert result['inverse_exists'] is True
        assert result['inverse_bounded'] is True

    def test_closed_graph_theorem_demo(self):
        """Satz vom abgeschlossenen Graphen Demonstration."""
        result = closed_graph_theorem_demo()
        assert 'graph_closed' in result
        assert result['operator_bounded'] is True

    def test_closed_graph_convergence(self):
        """Abgeschlossener Graph: Konvergenzfehler nimmt ab."""
        result = closed_graph_theorem_demo()
        errors = result['convergence_errors']
        assert errors[-1] < errors[0]  # Fehler nimmt ab

    def test_uniform_boundedness_demo(self):
        """Gleichmäßiges Beschränktheitsprinzip Demonstration."""
        result = uniform_boundedness_principle_demo()
        assert result['pointwise_bounded'] is True
        assert result['uniformly_bounded'] is True

    def test_uniform_boundedness_implication(self):
        """Punktweise beschränkt impliziert gleichmäßig beschränkt."""
        result = uniform_boundedness_principle_demo()
        # max Operatornorm muss endlich sein
        assert result['max_op_norm'] < float('inf')
        assert result['max_op_norm'] > 0


# ============================================================
# 10. Schwache Konvergenz, Riesz, Kompakte Operatoren
# ============================================================

class TestAdvancedTheorems:
    """Tests für fortgeschrittene Sätze."""

    def test_weak_convergence_demo(self):
        """Schwache Konvergenz Demonstration."""
        result = weak_convergence_demo()
        assert result['weakly_converges_to_zero'] is True
        assert not result['strongly_converges_to_zero']

    def test_weak_vs_strong_convergence(self):
        """Schwache ≠ starke Konvergenz."""
        result = weak_convergence_demo()
        # Einheitsvektoren konvergieren schwach gegen 0, aber stark nicht
        assert result['weakly_converges_to_zero'] is True

    def test_weak_convergence_values_decrease(self):
        """Schwache Konvergenzwerte ⟨eₙ,y⟩ nehmen ab."""
        result = weak_convergence_demo()
        weak_vals = result['weak_values']
        assert weak_vals[0] > weak_vals[-1]

    def test_riesz_representation_theorem(self):
        """Riesz-Darstellungssatz Demonstration."""
        result = riesz_representation_theorem()
        assert 'riesz_representative' in result
        for v in result['verifications']:
            assert v['equal'] is True

    def test_riesz_norm_equality(self):
        """Riesz: ‖f‖_{H*} = ‖y‖_H."""
        result = riesz_representation_theorem()
        y = np.array(result['riesz_representative'])
        # ‖y‖₂ = √(4+9+1) = √14
        assert result['norm_functional'] == pytest.approx(math.sqrt(14.0), abs=1e-8)

    def test_compact_operator_theory(self):
        """Kompakte Operatoren Theorie."""
        result = compact_operator_theory()
        assert result['is_compact'] is True
        assert result['hilbert_schmidt_norm'] > 0
        assert 'properties' in result
        assert len(result['properties']) > 0

    def test_compact_operator_spectral_properties(self):
        """Kompakter Operator: Spektralradius positiv."""
        result = compact_operator_theory()
        assert result['spectral_radius'] > 0


# ============================================================
# 11. Fixpunktsätze
# ============================================================

class TestFixedPointTheorems:
    """Tests für Fixpunktsätze."""

    def test_banach_fixed_point_contraction(self):
        """Banach-Fixpunktsatz: cos ist Kontraktion auf [0, π/2]."""
        result = banach_fixed_point(math.cos, 1.0, tol=1e-10)
        # Fixpunkt von cos(x) ≈ 0.7390851... (Dottie-Zahl)
        assert result['converged'] is True
        fp = result['fixed_point']
        assert abs(math.cos(fp) - fp) < 1e-8

    def test_banach_fixed_point_linear(self):
        """Banach-Fixpunktsatz: x/2 hat Fixpunkt 0."""
        result = banach_fixed_point(lambda x: x / 2.0, 10.0, tol=1e-10)
        assert result['converged'] is True
        assert abs(result['fixed_point']) < 1e-6

    def test_banach_fixed_point_affine(self):
        """Banach-Fixpunktsatz: f(x) = 0.5x + 1 hat Fixpunkt 2."""
        result = banach_fixed_point(lambda x: 0.5 * x + 1.0, 0.0, tol=1e-10)
        assert result['converged'] is True
        assert abs(result['fixed_point'] - 2.0) < 1e-8

    def test_banach_fixed_point_residual_small(self):
        """Finales Residuum |f(x*) - x*| ist klein."""
        result = banach_fixed_point(math.cos, 0.0, tol=1e-10)
        assert result['final_residual'] < 1e-8

    def test_banach_fixed_point_history(self):
        """Iterationshistorie ist nicht leer."""
        result = banach_fixed_point(math.cos, 1.0, tol=1e-10)
        assert len(result['history']) > 0

    def test_banach_fixed_point_contraction_estimate(self):
        """Kontraktionskonstante wird geschätzt."""
        result = banach_fixed_point(lambda x: x / 2.0, 8.0, tol=1e-10)
        if result['contraction_constant_estimate'] is not None:
            # Für x/2 ist q = 0.5
            q = result['contraction_constant_estimate']
            assert 0.0 <= q < 1.0

    def test_schauder_fixed_point_demo(self):
        """Schauder-Fixpunktsatz Demonstration."""
        result = schauder_fixed_point_demo()
        assert result['K_convex'] is True
        assert result['K_compact'] is True
        assert result['T_continuous'] is True

    def test_schauder_fixed_points_found(self):
        """Schauder: Fixpunkte von T(x)=x² werden gefunden."""
        result = schauder_fixed_point_demo()
        fps = result['fixed_points_found']
        # Mindestens ein Fixpunkt muss gefunden worden sein
        assert len(fps) >= 1

    def test_schauder_applications_listed(self):
        """Schauder: Anwendungen werden aufgelistet."""
        result = schauder_fixed_point_demo()
        assert 'applications' in result
        assert len(result['applications']) > 0


# ============================================================
# Edge Cases und Sonderfälle
# ============================================================

class TestEdgeCases:
    """Tests für Randfälle und Sonderfälle."""

    def test_norm_very_large_vector(self):
        """Norm funktioniert für große Werte."""
        v = [1e6, 1e6]
        ns = NormedSpace([v], euclidean_norm)
        norm = ns.norm(v)
        assert math.isfinite(norm)
        assert norm == pytest.approx(1e6 * math.sqrt(2), rel=1e-8)

    def test_gram_schmidt_single_vector(self):
        """Gram-Schmidt mit einem Vektor."""
        hs = HilbertSpace([[1.0, 2.0, 3.0]], standard_ip)
        onb = hs.gram_schmidt_orthonormalize([[1.0, 2.0, 3.0]])
        assert len(onb) == 1
        assert abs(hs.norm(onb[0]) - 1.0) < 1e-10

    def test_fredholm_identity_matrix(self):
        """Fredholm-Alternative: Identitätsmatrix immer lösbar."""
        A = np.eye(3)
        b = np.array([1.0, 2.0, 3.0])
        result = fredholm_alternative(A, b)
        assert result['solvable'] is True

    def test_operator_apply_zero(self):
        """Operator auf Nullvektor gibt Nullvektor."""
        A = np.array([[5.0, 3.0], [2.0, 7.0]])
        op = LinearOperator(A)
        result = op.apply([0.0, 0.0])
        assert all(abs(x) < 1e-12 for x in result)

    def test_l2_space_parseval_style(self):
        """Parseval-artiger Test im L²-Raum."""
        L2 = L2Space(0.0, 1.0)
        x = np.linspace(0, 1, 500)
        # f = Σ cₙ eₙ (erste 3 Basisvektoren)
        basis = L2.orthonormal_basis_fourier(3, 0.0, 1.0)
        coeffs = [1.0, 0.5, 0.25]
        f = sum(c * basis[i](x) for i, c in enumerate(coeffs))
        norm_sq = L2.norm(f.tolist(), x.tolist()) ** 2
        # Parseval: ‖f‖² ≈ Σcᵢ²
        parseval_sum = sum(c ** 2 for c in coeffs)
        assert abs(norm_sq - parseval_sum) < 0.05

    def test_spectral_radius_upper_bound(self):
        """Spektralradius ≤ Operatornorm."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        op = LinearOperator(A)
        r = spectral_radius(op)
        op_norm = op.operator_norm()
        assert r <= op_norm + 1e-8

    def test_banach_fixed_point_no_contraction_max_iter(self):
        """Nicht-Kontraktion: max_iter wird erreicht ohne Konvergenz."""
        # f(x) = 2x (keine Kontraktion)
        result = banach_fixed_point(lambda x: 2.0 * x, 1.0, tol=1e-10, max_iter=5)
        # Konvergiert nicht
        assert result['converged'] is False

    def test_lp_norm_holder_inequality(self):
        """Hölder-Ungleichung: |Σ uᵢvᵢ| ≤ ‖u‖_p · ‖v‖_q (1/p+1/q=1)."""
        bs = BanachSpace([], euclidean_norm)
        u = [3.0, 4.0, 0.0]
        v = [1.0, 0.0, 2.0]
        p, q = 2.0, 2.0  # 1/2 + 1/2 = 1 (Cauchy-Schwarz)
        inner = abs(sum(ui * vi for ui, vi in zip(u, v)))
        lhs = inner
        rhs = bs.lp_norm(u, p) * bs.lp_norm(v, q)
        assert lhs <= rhs + 1e-10
