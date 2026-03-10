"""
@file test_multilinear_algebra.py
@brief Tests für das Modul multilinear_algebra.py.
@description
    Testet alle Funktionen und Klassen des Multilinearen Algebra-Moduls:
    - Tensorprodukt (Vektoren, Matrizen, Klasse Tensor)
    - Äußere Algebra (Wedge-Produkt, äußere Potenzen, ExteriorAlgebra)
    - Gram-Matrix und Gram-Determinante
    - Symmetrische Algebra
    - Multilineare Formen (Symmetrie, Antisymmetrie, Rang, Signatur)

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math
import numpy as np
import pytest

# Quellverzeichnis zum Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from multilinear_algebra import (
    tensor_product_vectors,
    tensor_product_matrices,
    Tensor,
    tensor_contraction,
    outer_product,
    wedge_product,
    exterior_power,
    ExteriorAlgebra,
    gram_matrix,
    gram_determinant,
    SymmetricAlgebra,
    is_alternating,
    is_symmetric_tensor,
    multilinear_form_rank,
    signature_bilinear_form,
    _sign_permutation,
)


# ---------------------------------------------------------------------------
# Tests: Hilfsfunktionen
# ---------------------------------------------------------------------------

class TestSignPermutation:
    """Tests für die Vorzeichen-Berechnung von Permutationen."""

    def test_identity_is_even(self):
        """Identitätspermutation hat Vorzeichen +1."""
        assert _sign_permutation([0, 1, 2]) == 1

    def test_single_transposition_is_odd(self):
        """Einfache Transposition hat Vorzeichen -1."""
        assert _sign_permutation([1, 0, 2]) == -1
        assert _sign_permutation([0, 2, 1]) == -1

    def test_3_cycle_is_even(self):
        """3-Zykel (2 Transpositionen) hat Vorzeichen +1."""
        assert _sign_permutation([1, 2, 0]) == 1

    def test_reverse_4_is_even(self):
        """Umkehrpermutation [3,2,1,0] hat Vorzeichen (-1)^(4·3/2) = (-1)^6 = +1."""
        # [3,2,1,0]: 6 Inversionen → +1
        assert _sign_permutation([3, 2, 1, 0]) == 1

    def test_single_element(self):
        """Einelementige Permutation ist immer gerade."""
        assert _sign_permutation([0]) == 1


# ---------------------------------------------------------------------------
# Tests: Tensorprodukt
# ---------------------------------------------------------------------------

class TestTensorProductVectors:
    """Tests für das Tensorprodukt zweier Vektoren."""

    def test_basic_2d(self):
        """Tensorprodukt zweier 2D-Vektoren ergibt 2×2-Matrix."""
        u = [1.0, 2.0]
        v = [3.0, 4.0]
        result = tensor_product_vectors(u, v)
        # (u⊗v)_{ij} = u_i * v_j
        assert result[0][0] == pytest.approx(3.0)
        assert result[0][1] == pytest.approx(4.0)
        assert result[1][0] == pytest.approx(6.0)
        assert result[1][1] == pytest.approx(8.0)

    def test_shape(self):
        """Richtiges Format des Ergebnisses."""
        u = [1, 2, 3]
        v = [4, 5]
        result = tensor_product_vectors(u, v)
        assert len(result) == 3
        assert len(result[0]) == 2

    def test_standard_basis(self):
        """e₁⊗e₂ hat nur einen Nicht-Null-Eintrag."""
        e1 = [1.0, 0.0]
        e2 = [0.0, 1.0]
        result = tensor_product_vectors(e1, e2)
        assert result[0][0] == pytest.approx(0.0)
        assert result[0][1] == pytest.approx(1.0)
        assert result[1][0] == pytest.approx(0.0)
        assert result[1][1] == pytest.approx(0.0)

    def test_outer_product_relation(self):
        """Tensorprodukt = np.outer für Vektoren."""
        u = [1, -2, 3]
        v = [4, 0, -1]
        result = tensor_product_vectors(u, v)
        expected = np.outer(u, v)
        assert np.allclose(result, expected)


class TestTensorProductMatrices:
    """Tests für das Kronecker-Produkt (Tensorprodukt von Matrizen)."""

    def test_identity_kronecker(self):
        """I₂ ⊗ I₂ = I₄."""
        I = [[1, 0], [0, 1]]
        result = tensor_product_matrices(I, I)
        expected = np.eye(4)
        assert np.allclose(result, expected)

    def test_shape_2x2_kron_3x3(self):
        """2×2 ⊗ 3×3 ergibt 6×6-Matrix."""
        A = [[1, 2], [3, 4]]
        B = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        result = tensor_product_matrices(A, B)
        assert len(result) == 6
        assert len(result[0]) == 6

    def test_scalar_kronecker(self):
        """a·[[1]] ⊗ B = a·B."""
        A = [[3.0]]
        B = [[1, 2], [3, 4]]
        result = tensor_product_matrices(A, B)
        assert np.allclose(result, [[3, 6], [9, 12]])

    def test_known_kronecker(self):
        """Bekanntes Beispiel: [[1,2],[3,4]] ⊗ [[0,5],[6,7]]."""
        A = [[1, 2], [3, 4]]
        B = [[0, 5], [6, 7]]
        expected = np.kron(A, B)
        result = tensor_product_matrices(A, B)
        assert np.allclose(result, expected)


class TestTensorClass:
    """Tests für die Tensor-Klasse."""

    def test_create_tensor(self):
        """Erstelle einen einfachen (1,1)-Tensor."""
        components = np.eye(3)
        T = Tensor(components, covariant=1, contravariant=1)
        assert T.rank == 2
        assert T.covariant == 1
        assert T.contravariant == 1

    def test_add_tensors(self):
        """Addition zweier gleichartiger Tensoren."""
        A = np.ones((3, 3))
        B = np.eye(3)
        T1 = Tensor(A, 2, 0)
        T2 = Tensor(B, 2, 0)
        T_sum = T1 + T2
        expected = A + B
        assert np.allclose(T_sum.components, expected)

    def test_scalar_mul(self):
        """Skalarmultiplikation eines Tensors."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        T = Tensor(A, 2, 0)
        T2 = T * 3
        assert np.allclose(T2.components, 3 * A)

    def test_symmetrize_matrix(self):
        """Symmetrisierung einer Matrix ergibt symmetrische Matrix."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        T = Tensor(A, 2, 0)
        T_sym = T.symmetrize()
        # Ergebnis sollte symmetrisch sein
        assert np.allclose(T_sym.components, T_sym.components.T)
        # Kontrollwert: (A + Aᵀ) / 2
        expected = (A + A.T) / 2
        assert np.allclose(T_sym.components, expected)

    def test_antisymmetrize_matrix(self):
        """Antisymmetrisierung einer Matrix ergibt antisymmetrische Matrix."""
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        T = Tensor(A, 2, 0)
        T_anti = T.antisymmetrize()
        # Ergebnis sollte antisymmetrisch sein: T = -Tᵀ
        assert np.allclose(T_anti.components, -T_anti.components.T)
        expected = (A - A.T) / 2
        assert np.allclose(T_anti.components, expected)

    def test_contraction_trace(self):
        """Kontraktion einer Matrix über beide Indizes liefert die Spur."""
        A = np.array([[1.0, 2.0, 3.0],
                      [4.0, 5.0, 6.0],
                      [7.0, 8.0, 9.0]])
        T = Tensor(A, 1, 1)  # (1,1)-Tensor
        T_contracted = T.contract(0, 1)
        # Spur = 1 + 5 + 9 = 15
        assert T_contracted.components == pytest.approx(15.0)

    def test_type_error_add(self):
        """Addition von Tensoren verschiedenen Typs wirft Fehler."""
        T1 = Tensor(np.eye(2), 2, 0)
        T2 = Tensor(np.eye(2), 1, 1)
        with pytest.raises(ValueError):
            T1 + T2


class TestTensorContraction:
    """Tests für die standalone Kontraktion-Funktion."""

    def test_matrix_trace(self):
        """Spur einer 3×3-Matrix via Kontraktion."""
        A = np.array([[2, 0, 0], [0, 3, 0], [0, 0, 5]], dtype=float)
        result = tensor_contraction(A, 0, 1)
        assert result == pytest.approx(10.0)

    def test_same_index_error(self):
        """Gleiche Indexachsen werfen ValueError."""
        A = np.eye(3)
        with pytest.raises(ValueError):
            tensor_contraction(A, 1, 1)


class TestOuterProduct:
    """Tests für das äußere Produkt."""

    def test_basic(self):
        """Grundlegendes äußeres Produkt."""
        a = np.array([1.0, 2.0])
        b = np.array([3.0, 4.0])
        result = outer_product(a, b)
        expected = np.array([[3.0, 4.0], [6.0, 8.0]])
        assert np.allclose(result, expected)


# ---------------------------------------------------------------------------
# Tests: Äußere Algebra
# ---------------------------------------------------------------------------

class TestWedgeProduct:
    """Tests für das Wedge-Produkt."""

    def test_standard_basis_2d(self):
        """e₁∧e₂ in ℝ² hat Koeffizient 1."""
        e1 = [1.0, 0.0]
        e2 = [0.0, 1.0]
        result = wedge_product(e1, e2)
        # Einziger Basisvector: e₁∧e₂, Koeffizient = 1
        assert result[0] == pytest.approx(1.0)

    def test_anticommutativity(self):
        """u∧v = -(v∧u) (Antisymmetrie)."""
        u = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        uv = wedge_product(u, v)
        vu = wedge_product(v, u)
        assert np.allclose(uv, [-x for x in vu])

    def test_parallel_vectors(self):
        """u∧u = 0 für parallele Vektoren."""
        u = [1.0, 2.0, 3.0]
        result = wedge_product(u, u)
        assert np.allclose(result, 0)

    def test_3d_cross_product_relation(self):
        """In ℝ³ entsprechen die Komponenten von u∧v dem Kreuzprodukt."""
        u = [1.0, 0.0, 0.0]
        v = [0.0, 1.0, 0.0]
        result = wedge_product(u, v)
        # (u∧v)_{01}=1, (u∧v)_{02}=0, (u∧v)_{12}=0
        assert result[0] == pytest.approx(1.0)  # e₁∧e₂ Komponente
        assert result[1] == pytest.approx(0.0)  # e₁∧e₃ Komponente
        assert result[2] == pytest.approx(0.0)  # e₂∧e₃ Komponente

    def test_dimension_error(self):
        """Vektoren verschiedener Dimension werfen ValueError."""
        with pytest.raises(ValueError):
            wedge_product([1, 2], [1, 2, 3])


class TestExteriorPower:
    """Tests für die äußere Potenz."""

    def test_2_vectors_in_3d(self):
        """2-Potenz von Standardbasis-Vektoren in ℝ³."""
        e1 = [1.0, 0.0, 0.0]
        e2 = [0.0, 1.0, 0.0]
        result = exterior_power([e1, e2], 2)
        # e₁∧e₂ hat Koeffizient 1 bei Index (0,1)
        # Koeffizienten: {(0,1): 1, (0,2): 0, (1,2): 0}
        assert result[0] == pytest.approx(1.0)  # Index (0,1)
        assert result[1] == pytest.approx(0.0)  # Index (0,2)
        assert result[2] == pytest.approx(0.0)  # Index (1,2)

    def test_determinant_as_exterior_power(self):
        """k-Potenz für k = n liefert die Determinante."""
        v1 = [1.0, 0.0, 0.0]
        v2 = [0.0, 2.0, 0.0]
        v3 = [0.0, 0.0, 3.0]
        result = exterior_power([v1, v2, v3], 3)
        # Determinante der Diagonalmatrix = 1·2·3 = 6
        assert result[0] == pytest.approx(6.0)

    def test_linearly_dependent_zero(self):
        """Linear abhängige Vektoren ergeben 0."""
        u = [1.0, 2.0, 3.0]
        v = [2.0, 4.0, 6.0]  # v = 2u
        result = exterior_power([u, v], 2)
        assert np.allclose(result, 0, atol=1e-10)

    def test_wrong_count_error(self):
        """Falsche Anzahl Vektoren wirft ValueError."""
        with pytest.raises(ValueError):
            exterior_power([[1, 2], [3, 4]], 3)


class TestExteriorAlgebra:
    """Tests für die ExteriorAlgebra-Klasse."""

    def test_basis_0_forms(self):
        """Λ⁰(V) hat genau eine Basis: leeres Tupel."""
        ea = ExteriorAlgebra(3)
        basis = ea.basis_k_forms(0)
        assert len(basis) == 1

    def test_basis_1_forms(self):
        """Λ¹(V) hat Dimension n."""
        ea = ExteriorAlgebra(4)
        basis = ea.basis_k_forms(1)
        assert len(basis) == 4

    def test_basis_2_forms_3d(self):
        """Λ²(ℝ³) hat Dimension C(3,2) = 3."""
        ea = ExteriorAlgebra(3)
        basis = ea.basis_k_forms(2)
        assert len(basis) == 3

    def test_basis_top_form(self):
        """Λⁿ(ℝⁿ) hat genau eine Basis (Volumform)."""
        ea = ExteriorAlgebra(5)
        basis = ea.basis_k_forms(5)
        assert len(basis) == 1
        assert basis[0] == (0, 1, 2, 3, 4)

    def test_wedge_basic(self):
        """Keilprodukt zweier 1-Formen."""
        ea = ExteriorAlgebra(3)
        # e₁ als 1-Form
        form1 = {(0,): 1.0}
        form2 = {(1,): 1.0}
        result = ea.wedge(form1, form2)
        # e₁∧e₂ → Koeffizient bei (0,1)
        assert (0, 1) in result
        assert result[(0, 1)] == pytest.approx(1.0)

    def test_wedge_anticommutativity(self):
        """e₁∧e₂ = -(e₂∧e₁)."""
        ea = ExteriorAlgebra(3)
        form1 = {(0,): 1.0}
        form2 = {(1,): 1.0}
        r1 = ea.wedge(form1, form2)
        r2 = ea.wedge(form2, form1)
        assert (0, 1) in r1
        assert (0, 1) in r2
        assert r1[(0, 1)] == pytest.approx(-r2[(0, 1)])

    def test_wedge_same_form_zero(self):
        """α∧α = 0 für ungerade Grade (Antisymmetrie)."""
        ea = ExteriorAlgebra(3)
        form = {(0,): 1.0}
        result = ea.wedge(form, form)
        assert len(result) == 0  # Ergebnis ist die Nullform

    def test_hodge_dual_volume_form(self):
        """Hodge-Dual der Einheitsvolumform ist 1 (Skalar)."""
        ea = ExteriorAlgebra(3)
        # Volumform e₁∧e₂∧e₃
        vol_form = {(0, 1, 2): 1.0}
        dual = ea.hodge_dual(vol_form)
        # ⋆(e₁∧e₂∧e₃) = 1 (leeres Tupel, 0-Form)
        assert () in dual
        assert dual[()] == pytest.approx(1.0)

    def test_hodge_dual_1_form(self):
        """Hodge-Dual einer 1-Form in ℝ³ ergibt 2-Form."""
        ea = ExteriorAlgebra(3)
        # ⋆(e₁) = e₂∧e₃ in ℝ³
        form = {(0,): 1.0}
        dual = ea.hodge_dual(form)
        # e₁∧e₂∧e₃ mit e₁ weggelassen: e₂∧e₃ = (1,2)
        assert (1, 2) in dual


# ---------------------------------------------------------------------------
# Tests: Gram-Matrix
# ---------------------------------------------------------------------------

class TestGramMatrix:
    """Tests für die Gram-Matrix."""

    def test_orthonormal_basis(self):
        """Gram-Matrix einer ONB ist die Einheitsmatrix."""
        vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        G = gram_matrix(vectors)
        assert np.allclose(G, np.eye(3))

    def test_symmetry(self):
        """Gram-Matrix ist symmetrisch."""
        vectors = [[1.0, 2.0], [3.0, 4.0]]
        G = gram_matrix(vectors)
        assert np.allclose(np.array(G), np.array(G).T)

    def test_positive_semidefinite(self):
        """Gram-Matrix ist positiv semidefinit."""
        vectors = [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]]
        G = np.array(gram_matrix(vectors))
        eigenvalues = np.linalg.eigvalsh(G)
        assert np.all(eigenvalues >= -1e-10)


class TestGramDeterminant:
    """Tests für die Gram-Determinante."""

    def test_orthonormal_vectors(self):
        """Gram-Det von ONB = 1² = 1."""
        vectors = [[1.0, 0.0], [0.0, 1.0]]
        assert gram_determinant(vectors) == pytest.approx(1.0)

    def test_parallel_vectors_zero(self):
        """Linear abhängige Vektoren → Gram-Det = 0."""
        vectors = [[1.0, 0.0], [2.0, 0.0]]  # Parallel
        assert abs(gram_determinant(vectors)) < 1e-10

    def test_volume_squared(self):
        """Gram-Det = Volumen² des Parallelotops."""
        # Einheitsquadrat: Fläche = 1, Gram-Det = 1
        vectors = [[1.0, 0.0], [0.0, 1.0]]
        assert gram_determinant(vectors) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Tests: Symmetrische Algebra
# ---------------------------------------------------------------------------

class TestSymmetricAlgebra:
    """Tests für die symmetrische Algebra."""

    def test_dim_s1(self):
        """dim(S¹(V)) = dim(V) = n."""
        sa = SymmetricAlgebra(4)
        assert sa.symmetric_power_dim(1) == 4

    def test_dim_s0(self):
        """dim(S⁰(V)) = 1 (Skalare)."""
        sa = SymmetricAlgebra(5)
        assert sa.symmetric_power_dim(0) == 1

    def test_dim_s2_formula(self):
        """dim(S²(ℝⁿ)) = C(n+1,2) = n(n+1)/2."""
        sa = SymmetricAlgebra(3)
        # C(3+2-1, 2) = C(4,2) = 6
        assert sa.symmetric_power_dim(2) == 6

    def test_dim_s3(self):
        """dim(S³(ℝ²)) = C(2+3-1,3) = C(4,3) = 4."""
        sa = SymmetricAlgebra(2)
        assert sa.symmetric_power_dim(3) == 4

    def test_polarization_basic(self):
        """Polarisierung eines einfachen Monomials."""
        sa = SymmetricAlgebra(2)
        # p(x,y) = x² mit Koeffizient 1
        poly = {(2, 0): 1.0}  # x²
        result = sa.polarization(poly, 2)
        # Polarisiert: Koeffizient / (2!/2!) = 1/1 = 1 → nein
        # Faktor: k! / (e1! * e2!) = 2! / (2! * 0!) = 1
        assert (2, 0) in result


# ---------------------------------------------------------------------------
# Tests: Multilineare Formen
# ---------------------------------------------------------------------------

class TestIsAlternating:
    """Tests für die Alternierend-Prüfung."""

    def test_zero_tensor(self):
        """Nulltensor ist trivial alternierend."""
        T = np.zeros((3, 3))
        assert is_alternating(T)

    def test_antisymmetric_matrix(self):
        """Antisymmetrische Matrix ist alternierend."""
        T = np.array([[0, 1, -2], [-1, 0, 3], [2, -3, 0]], dtype=float)
        assert is_alternating(T)

    def test_symmetric_matrix_not_alternating(self):
        """Symmetrische Matrix (außer 0) ist nicht alternierend."""
        T = np.eye(3)
        assert not is_alternating(T)

    def test_scalar(self):
        """Skalar (0D) ist trivial alternierend."""
        T = np.array(5.0)
        assert is_alternating(T)


class TestIsSymmetricTensor:
    """Tests für die Symmetrie-Prüfung."""

    def test_identity_is_symmetric(self):
        """Einheitsmatrix ist symmetrisch."""
        T = np.eye(4)
        assert is_symmetric_tensor(T)

    def test_diagonal_is_symmetric(self):
        """Diagonalmatrix ist symmetrisch."""
        T = np.diag([1.0, 2.0, 3.0])
        assert is_symmetric_tensor(T)

    def test_antisymmetric_not_symmetric(self):
        """Antisymmetrische Matrix ist nicht symmetrisch (außer Nullmatrix)."""
        T = np.array([[0, 1], [-1, 0]], dtype=float)
        assert not is_symmetric_tensor(T)

    def test_random_symmetric(self):
        """Explizit symmetrisierte Matrix ist symmetrisch."""
        A = np.array([[1.0, 2.0, 3.0], [2.0, 5.0, 6.0], [3.0, 6.0, 9.0]])
        assert is_symmetric_tensor(A)


class TestMultilinearFormRank:
    """Tests für den Rang einer Bilinearform."""

    def test_identity_rank(self):
        """Einheitsmatrix hat vollen Rang."""
        M = np.eye(4)
        assert multilinear_form_rank(M) == 4

    def test_zero_rank(self):
        """Nullmatrix hat Rang 0."""
        M = np.zeros((3, 3))
        assert multilinear_form_rank(M) == 0

    def test_rank_2_matrix(self):
        """Matrix mit Rang 2."""
        M = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]], dtype=float)
        assert multilinear_form_rank(M) == 2


class TestSignatureBilinearForm:
    """Tests für die Sylvester-Signatur."""

    def test_euclidean_metric(self):
        """Euklidische Metrik: Signatur (n, 0)."""
        M = np.eye(3)
        p, q = signature_bilinear_form(M)
        assert p == 3
        assert q == 0

    def test_minkowski_metric(self):
        """Minkowski-Metrik: Signatur (1, 3) oder (3, 1)."""
        M = np.diag([1.0, -1.0, -1.0, -1.0])
        p, q = signature_bilinear_form(M)
        assert (p == 1 and q == 3) or (p == 3 and q == 1)

    def test_indefinite_form(self):
        """Indefinite Form: p > 0 und q > 0."""
        M = np.diag([1.0, -1.0])
        p, q = signature_bilinear_form(M)
        assert p == 1
        assert q == 1

    def test_negative_definite(self):
        """Negativ definite Form: Signatur (0, n)."""
        M = -np.eye(4)
        p, q = signature_bilinear_form(M)
        assert p == 0
        assert q == 4

    def test_degenerate_form(self):
        """Entartete Form: p + q < n."""
        M = np.diag([1.0, 0.0, -1.0])
        p, q = signature_bilinear_form(M)
        assert p == 1
        assert q == 1
        assert p + q == 2  # Rang ist 2, nicht 3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
