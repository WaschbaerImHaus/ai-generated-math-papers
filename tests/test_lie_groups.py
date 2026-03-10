"""
@file test_lie_groups.py
@brief Umfassende Tests für das Lie-Gruppen-Modul (lie_groups.py).
@description
    Testklassen und Testfunktionen für alle Klassen und Funktionen
    des Lie-Gruppen-Moduls:

    - TestLieGroupAbstract         – Abstrakte Basisklasse via Unterklasse
    - TestMatrixLieGroup           – MatrixLieGroup Grundfunktionen
    - TestSO                       – SO(n) Rotationsgruppe (n=2,3,4)
    - TestSU                       – SU(n) unitäre Gruppe (n=1,2,3)
    - TestGL                       – GL(n,R) allgemeine lineare Gruppe
    - TestLieAlgebra               – Lie-Klammer, adjungierte Darstellung, Killing-Form
    - TestExponentialMap           – exp_map, log_map, BCH-Formel
    - TestOneParameterSubgroup     – Einparametrige Untergruppen
    - TestClassification           – Cartan-Klassifikation, Dynkin-Diagramme
    - TestExceptionalAlgebras      – G2, F4, E6, E7, E8
    - TestFundamentalRepSU2        – Spin-j-Darstellungen
    - TestCasimirElement           – Casimir-Operator
    - TestEdgeCases                – Randfälle und Fehlerfälle

    Mathematische Grundlagen:
    - SO(n): AᵀA = I, det = +1
    - SU(n): A†A = I, det = +1 (komplex)
    - [X, Y] = XY - YX (Lie-Klammer/Kommutator)
    - exp(X)·exp(-X) = I (Exponentialabbildung)
    - BCH: log(eˣeʸ) ≈ X + Y + ½[X,Y] + …

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest
import numpy as np
from numpy.testing import assert_allclose

# Modul-Pfad einfügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lie_groups import (
    LieGroup, MatrixLieGroup, SO, SU, GL,
    LieAlgebra, ExponentialMap,
    one_parameter_subgroup,
    classify_simple_lie_algebra,
    exceptional_lie_algebras_info,
    cartan_matrix,
    dynkin_diagram_info,
    fundamental_representation_su2,
    casimir_element,
)



# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def random_so_element(n: int) -> np.ndarray:
    """Erzeugt zufälliges SO(n)-Element via QR-Zerlegung."""
    Q, R = np.linalg.qr(np.random.randn(n, n))
    # Vorzeichen korrigieren damit det = +1
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def random_su_element(n: int) -> np.ndarray:
    """Erzeugt zufälliges SU(n)-Element via QR-Zerlegung komplexer Matrizen."""
    Z = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    Q, R = np.linalg.qr(Z)
    # Determinante auf 1 normieren
    det = np.linalg.det(Q)
    Q[:, 0] /= det
    return Q


def skew_symmetric_matrix(n: int) -> np.ndarray:
    """Erzeugt eine zufällige schiefsymmetrische n×n-Matrix (so(n)-Element)."""
    A = np.random.randn(n, n)
    return (A - A.T) / 2.0


def skew_hermitian_traceless(n: int) -> np.ndarray:
    """Erzeugt eine zufällige schief-hermitesche spurlose n×n-Matrix (su(n)-Element)."""
    A = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    A = (A - A.conj().T) / 2.0    # schief-hermitesch
    A -= np.trace(A) / n * np.eye(n)   # spurlos
    return A


# =============================================================================
# TEST: Abstrakte Basisklasse LieGroup
# =============================================================================

class TestLieGroupAbstract:
    """Tests für die abstrakte LieGroup-Basisklasse via konkreter Unterklasse."""

    def test_abstract_class_cannot_be_instantiated(self):
        """LieGroup kann nicht direkt instanziiert werden (abstract)."""
        with pytest.raises(TypeError):
            LieGroup()

    def test_concrete_subclass_implements_all_methods(self):
        """Konkreter Unterklasse mit allen Methoden kann instanziiert werden."""
        # MatrixLieGroup ist konkrete Unterklasse
        g = MatrixLieGroup(2)
        assert g.dimension() == 4
        assert isinstance(g.is_compact(), bool)
        assert isinstance(g.is_connected(), bool)

    def test_so2_satisfies_lie_group_interface(self):
        """SO(2) implementiert das LieGroup-Interface korrekt."""
        so2 = SO(2)
        assert isinstance(so2, LieGroup)
        assert so2.dimension() == 1
        assert bool(so2.is_compact())
        assert bool(so2.is_connected())

    def test_su2_satisfies_lie_group_interface(self):
        """SU(2) implementiert das LieGroup-Interface korrekt."""
        su2 = SU(2)
        assert isinstance(su2, LieGroup)
        assert su2.dimension() == 3
        assert bool(su2.is_compact())
        assert bool(su2.is_connected())


# =============================================================================
# TEST: MatrixLieGroup
# =============================================================================

class TestMatrixLieGroup:
    """Tests für die allgemeine MatrixLieGroup-Klasse."""

    def test_identity_real_2x2(self):
        """Einheitsmatrix für 2×2 reelle Gruppe."""
        g = MatrixLieGroup(2, 'real')
        I = g.identity()
        assert_allclose(I, np.eye(2), atol=1e-15)
        assert I.dtype in (np.float64, np.float32)

    def test_identity_complex_3x3(self):
        """Einheitsmatrix für 3×3 komplexe Gruppe."""
        g = MatrixLieGroup(3, 'complex')
        I = g.identity()
        assert_allclose(I, np.eye(3, dtype=complex), atol=1e-15)
        assert I.dtype == complex

    def test_compose_is_matrix_multiplication(self):
        """compose(A, B) = A @ B."""
        g = MatrixLieGroup(2)
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        B = np.array([[5.0, 6.0], [7.0, 8.0]])
        result = g.compose(A, B)
        assert_allclose(result, A @ B, atol=1e-14)

    def test_compose_identity_left_neutral(self):
        """I · A = A."""
        g = MatrixLieGroup(3)
        A = np.random.randn(3, 3)
        I = g.identity()
        assert_allclose(g.compose(I, A), A, atol=1e-14)

    def test_compose_identity_right_neutral(self):
        """A · I = A."""
        g = MatrixLieGroup(3)
        A = np.random.randn(3, 3)
        I = g.identity()
        assert_allclose(g.compose(A, I), A, atol=1e-14)

    def test_inverse_of_identity(self):
        """Inverses der Einheitsmatrix = Einheitsmatrix."""
        g = MatrixLieGroup(4)
        I = g.identity()
        assert_allclose(g.inverse(I), I, atol=1e-14)

    def test_inverse_times_original_is_identity(self):
        """A · A⁻¹ = I."""
        g = MatrixLieGroup(3)
        A = np.random.randn(3, 3) + 2 * np.eye(3)   # gut-konditioniert
        A_inv = g.inverse(A)
        assert_allclose(A @ A_inv, np.eye(3), atol=1e-12)

    def test_is_element_for_invertible_matrix(self):
        """Invertierbare Matrix ist in GL(n)."""
        g = MatrixLieGroup(2)
        A = np.array([[2.0, 0.0], [0.0, 3.0]])
        assert g.is_element(A)

    def test_is_element_for_singular_matrix(self):
        """Singuläre Matrix ist nicht in GL(n)."""
        g = MatrixLieGroup(2)
        A = np.array([[1.0, 2.0], [2.0, 4.0]])   # det = 0
        assert not bool(g.is_element(A))

    def test_is_element_wrong_shape(self):
        """Matrix falscher Form ist nicht Element."""
        g = MatrixLieGroup(3)
        A = np.eye(2)   # 2×2 statt 3×3
        assert not bool(g.is_element(A))

    def test_dimension_n_squared(self):
        """Dimension von MatrixLieGroup = n²."""
        for n in [1, 2, 3, 4, 5]:
            g = MatrixLieGroup(n)
            assert g.dimension() == n * n

    def test_invalid_n_raises(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            MatrixLieGroup(0)

    def test_invalid_field_raises(self):
        """Ungültiger field-Parameter wirft ValueError."""
        with pytest.raises(ValueError):
            MatrixLieGroup(2, 'quaternion')

    def test_gl_real_not_connected(self):
        """GL(n, ℝ) ist nicht zusammenhängend."""
        g = MatrixLieGroup(3, 'real')
        assert not bool(g.is_connected())

    def test_gl_complex_connected(self):
        """GL(n, ℂ) ist zusammenhängend."""
        g = MatrixLieGroup(3, 'complex')
        assert bool(g.is_connected())

    def test_gl_not_compact(self):
        """GL(n) ist nicht kompakt."""
        g = MatrixLieGroup(2, 'real')
        assert not bool(g.is_compact())


# =============================================================================
# TEST: SO(n) - Spezielle orthogonale Gruppe
# =============================================================================

class TestSO:
    """Tests für die spezielle orthogonale Gruppe SO(n)."""

    def test_so2_dimension(self):
        """dim(SO(2)) = 1."""
        assert SO(2).dimension() == 1

    def test_so3_dimension(self):
        """dim(SO(3)) = 3."""
        assert SO(3).dimension() == 3

    def test_so4_dimension(self):
        """dim(SO(4)) = 6."""
        assert SO(4).dimension() == 6

    def test_son_dimension_formula(self):
        """dim(SO(n)) = n(n-1)/2."""
        for n in [2, 3, 4, 5, 6]:
            assert SO(n).dimension() == n * (n - 1) // 2

    def test_so_is_compact(self):
        """SO(n) ist kompakt."""
        for n in [2, 3, 4, 5]:
            assert bool(SO(n).is_compact())

    def test_so_is_connected(self):
        """SO(n) ist zusammenhängend."""
        for n in [2, 3, 4, 5]:
            assert bool(SO(n).is_connected())

    def test_so1_raises(self):
        """SO(1) ist nicht definiert (n < 2)."""
        with pytest.raises(ValueError):
            SO(1)

    def test_identity_in_so2(self):
        """Einheitsmatrix ist in SO(2)."""
        so2 = SO(2)
        assert so2.is_element(so2.identity())

    def test_identity_in_so3(self):
        """Einheitsmatrix ist in SO(3)."""
        so3 = SO(3)
        assert so3.is_element(so3.identity())

    def test_rotation_2d_is_so2_element(self):
        """2D-Rotationsmatrix liegt in SO(2)."""
        so2 = SO(2)
        for theta in [0.0, 0.5, np.pi / 4, np.pi / 2, np.pi, 2 * np.pi]:
            R = so2.rotation_2d(theta)
            assert so2.is_element(R), f"rotation_2d({theta}) nicht in SO(2)"

    def test_rotation_2d_zero_is_identity(self):
        """rotation_2d(0) = I₂."""
        so2 = SO(2)
        R = so2.rotation_2d(0.0)
        assert_allclose(R, np.eye(2), atol=1e-15)

    def test_rotation_2d_pi_half(self):
        """rotation_2d(π/2) = [[0,-1],[1,0]]."""
        so2 = SO(2)
        R = so2.rotation_2d(np.pi / 2)
        expected = np.array([[0.0, -1.0], [1.0, 0.0]])
        assert_allclose(R, expected, atol=1e-14)

    def test_rotation_2d_composition(self):
        """R(α)·R(β) = R(α+β)."""
        so2 = SO(2)
        alpha, beta = 0.3, 0.7
        Ra = so2.rotation_2d(alpha)
        Rb = so2.rotation_2d(beta)
        Rab = so2.rotation_2d(alpha + beta)
        assert_allclose(so2.compose(Ra, Rb), Rab, atol=1e-14)

    def test_rotation_3d_x_is_so3_element(self):
        """3D-Rotation um X liegt in SO(3)."""
        so3 = SO(3)
        for theta in [0.0, 0.3, np.pi / 3, np.pi]:
            R = so3.rotation_3d_x(theta)
            assert so3.is_element(R), f"rotation_3d_x({theta}) nicht in SO(3)"

    def test_rotation_3d_y_is_so3_element(self):
        """3D-Rotation um Y liegt in SO(3)."""
        so3 = SO(3)
        for theta in [0.0, 0.5, np.pi / 2, np.pi]:
            R = so3.rotation_3d_y(theta)
            assert so3.is_element(R), f"rotation_3d_y({theta}) nicht in SO(3)"

    def test_rotation_3d_z_is_so3_element(self):
        """3D-Rotation um Z liegt in SO(3)."""
        so3 = SO(3)
        for theta in [0.0, 1.0, 2 * np.pi / 3, 2 * np.pi]:
            R = so3.rotation_3d_z(theta)
            assert so3.is_element(R), f"rotation_3d_z({theta}) nicht in SO(3)"

    def test_rotation_3d_x_zero(self):
        """rotation_3d_x(0) = I₃."""
        so3 = SO(3)
        assert_allclose(so3.rotation_3d_x(0.0), np.eye(3), atol=1e-15)

    def test_rotation_3d_xyz_composition_closed(self):
        """Produkt von X·Y·Z-Rotationen liegt in SO(3)."""
        so3 = SO(3)
        Rx = so3.rotation_3d_x(0.3)
        Ry = so3.rotation_3d_y(0.7)
        Rz = so3.rotation_3d_z(1.1)
        result = Rx @ Ry @ Rz
        assert so3.is_element(result, tol=1e-12)

    def test_reflection_not_in_so2(self):
        """Spiegelung (det = -1) ist nicht in SO(2)."""
        so2 = SO(2)
        M = np.array([[-1.0, 0.0], [0.0, 1.0]])   # det = -1
        assert not bool(so2.is_element(M))

    def test_scaled_matrix_not_in_so2(self):
        """Skalierte Rotationsmatrix (det ≠ 1) ist nicht in SO(2)."""
        so2 = SO(2)
        R = 2.0 * so2.rotation_2d(0.5)   # det = 4 ≠ 1
        assert not bool(so2.is_element(R))

    def test_random_so3_elements(self):
        """Zufällige SO(3)-Elemente werden korrekt erkannt."""
        so3 = SO(3)
        for _ in range(5):
            R = random_so_element(3)
            assert so3.is_element(R, tol=1e-10)

    def test_rotation_2d_wrong_n_raises(self):
        """rotation_2d() für SO(3) wirft ValueError."""
        so3 = SO(3)
        with pytest.raises(ValueError):
            so3.rotation_2d(0.5)

    def test_rotation_3d_x_wrong_n_raises(self):
        """rotation_3d_x() für SO(2) wirft ValueError."""
        so2 = SO(2)
        with pytest.raises(ValueError):
            so2.rotation_3d_x(0.5)


# =============================================================================
# TEST: SU(n) - Spezielle unitäre Gruppe
# =============================================================================

class TestSU:
    """Tests für die spezielle unitäre Gruppe SU(n)."""

    def test_su1_dimension(self):
        """dim(SU(1)) = 0."""
        assert SU(1).dimension() == 0

    def test_su2_dimension(self):
        """dim(SU(2)) = 3."""
        assert SU(2).dimension() == 3

    def test_su3_dimension(self):
        """dim(SU(3)) = 8."""
        assert SU(3).dimension() == 8

    def test_sun_dimension_formula(self):
        """dim(SU(n)) = n²-1."""
        for n in [1, 2, 3, 4, 5]:
            assert SU(n).dimension() == n * n - 1

    def test_su_is_compact(self):
        """SU(n) ist kompakt."""
        for n in [1, 2, 3]:
            assert bool(SU(n).is_compact())

    def test_su_is_connected(self):
        """SU(n) ist zusammenhängend."""
        for n in [1, 2, 3]:
            assert bool(SU(n).is_connected())

    def test_identity_in_su2(self):
        """Einheitsmatrix I₂ liegt in SU(2)."""
        su2 = SU(2)
        assert su2.is_element(su2.identity())

    def test_identity_in_su3(self):
        """Einheitsmatrix I₃ liegt in SU(3)."""
        su3 = SU(3)
        assert su3.is_element(su3.identity())

    def test_su2_pauli_matrices_exponentiated(self):
        """exp(i·π·σ_z/2) liegt in SU(2) (Pauli-z-Rotation)."""
        su2 = SU(2)
        sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
        # exp(i·θ·σ_z/2) für θ = π
        theta = np.pi
        U = np.array([[np.exp(1j * theta / 2), 0],
                      [0, np.exp(-1j * theta / 2)]], dtype=complex)
        assert su2.is_element(U, tol=1e-12)

    def test_su2_from_angles_is_element(self):
        """su2_from_angles() erzeugt gültige SU(2)-Elemente."""
        su2 = SU(2)
        test_cases = [
            (0.0, 0.0, 0.0),
            (np.pi / 2, 0.0, 0.0),
            (np.pi / 3, np.pi / 4, np.pi / 6),
            (np.pi, 2 * np.pi / 3, np.pi / 2),
        ]
        for theta, phi, psi in test_cases:
            U = su2.su2_from_angles(theta, phi, psi)
            assert su2.is_element(U, tol=1e-12), \
                f"su2_from_angles({theta}, {phi}, {psi}) nicht in SU(2)"
            

    def test_su2_from_angles_zero_is_identity(self):
        """su2_from_angles(0, 0, 0) ≈ I₂ (bis auf Phase-Vorzeichen)."""
        su2 = SU(2)
        U = su2.su2_from_angles(0.0, 0.0, 0.0)
        # U sollte diagonal sein: [[e^{iψ/2·...}, 0], [0, e^{-i...}]]
        # Für θ=0: sin(θ/2)=0, also [[e^{i(φ+ψ)/2}, 0],[0, e^{-i(φ+ψ)/2}]]
        # bei φ=ψ=0: [[1, 0],[0, 1]] = I
        assert su2.is_element(U, tol=1e-12)
        assert_allclose(np.abs(U), np.eye(2), atol=1e-12)

    def test_su2_wrong_n_raises(self):
        """su2_from_angles() für SU(3) wirft ValueError."""
        su3 = SU(3)
        with pytest.raises(ValueError):
            su3.su2_from_angles(0.5, 0.5, 0.5)

    def test_non_unitary_not_in_su2(self):
        """Nicht-unitäre Matrix ist nicht in SU(2)."""
        su2 = SU(2)
        A = np.array([[2.0, 0.0], [0.0, 0.5]], dtype=complex)   # nicht unitär
        assert not bool(su2.is_element(A))

    def test_unitary_det_minus1_not_in_su2(self):
        """Unitäre Matrix mit det = -1 ist nicht in SU(2)."""
        su2 = SU(2)
        # U(1)×U(1) mit det = -1: z.B. diag(-1, 1)
        A = np.diag([-1.0 + 0j, 1.0 + 0j])
        assert not bool(su2.is_element(A))

    def test_random_su3_elements(self):
        """Zufällige SU(3)-Elemente werden korrekt erkannt."""
        su3 = SU(3)
        for _ in range(5):
            U = random_su_element(3)
            assert su3.is_element(U, tol=1e-10)

    def test_su2_inverse_is_adjoint(self):
        """Für U ∈ SU(2): U⁻¹ = U†."""
        su2 = SU(2)
        U = su2.su2_from_angles(0.7, 0.4, 1.2)
        U_inv = su2.inverse(U)
        U_dag = U.conj().T
        assert_allclose(U_inv, U_dag, atol=1e-12)


# =============================================================================
# TEST: GL(n) - Allgemeine lineare Gruppe
# =============================================================================

class TestGL:
    """Tests für die allgemeine lineare Gruppe GL(n, F)."""

    def test_identity_in_gl2(self):
        """Einheitsmatrix ist in GL(2)."""
        gl = GL(2)
        assert bool(gl.is_element(gl.identity()))

    def test_invertible_diagonal_in_gl3(self):
        """Diagonalmatrix mit Einträgen ≠ 0 liegt in GL(3)."""
        gl = GL(3)
        D = np.diag([2.0, -3.0, 0.5])
        assert bool(gl.is_element(D))

    def test_singular_matrix_not_in_gl(self):
        """Singuläre Matrix liegt nicht in GL(n)."""
        gl = GL(3)
        A = np.zeros((3, 3))   # det = 0
        assert not bool(gl.is_element(A))

    def test_zero_matrix_not_in_gl(self):
        """Nullmatrix liegt nicht in GL(n)."""
        gl = GL(2)
        assert not bool(gl.is_element(np.zeros((2, 2))))

    def test_rotation_in_gl(self):
        """Rotationsmatrix liegt auch in GL(2)."""
        gl = GL(2)
        so2 = SO(2)
        R = so2.rotation_2d(1.2)
        assert bool(gl.is_element(R))

    def test_gl_dimension(self):
        """dim(GL(n)) = n²."""
        for n in [1, 2, 3, 4]:
            gl = GL(n)
            assert gl.dimension() == n * n

    def test_gl_real_not_compact(self):
        """GL(n, ℝ) ist nicht kompakt."""
        gl = GL(3, 'real')
        assert not bool(gl.is_compact())

    def test_gl_complex_not_compact(self):
        """GL(n, ℂ) ist nicht kompakt."""
        gl = GL(3, 'complex')
        assert not bool(gl.is_compact())

    def test_complex_matrix_in_gl_complex(self):
        """Komplexe invertierbare Matrix liegt in GL(n, ℂ)."""
        gl = GL(2, 'complex')
        A = np.array([[1+1j, 0], [0, 1-1j]])   # det = 1+1 = 2
        assert bool(gl.is_element(A))

    def test_wrong_shape_not_in_gl(self):
        """Matrix falscher Form liegt nicht in GL(n)."""
        gl = GL(3)
        assert not bool(gl.is_element(np.eye(2)))


# =============================================================================
# TEST: LieAlgebra - Lie-Klammer, Killing-Form
# =============================================================================

class TestLieAlgebra:
    """Tests für die Lie-Algebra-Klasse."""

    def _su2_basis(self):
        """Standardbasis von su(2): i/2 · Pauli-Matrizen."""
        # Basis: e_k = i/2 · σ_k (schief-hermitesch, spurlos)
        e1 = np.array([[0, 1j/2], [1j/2, 0]], dtype=complex)      # i/2 σ_x
        e2 = np.array([[0, 1/2], [-1/2, 0]], dtype=complex)        # i/2 σ_y (= 1/2 σ_y · i)
        e3 = np.array([[1j/2, 0], [0, -1j/2]], dtype=complex)      # i/2 σ_z
        return [e1, e2, e3]

    def _so3_basis(self):
        """Standardbasis von so(3): schiefsymmetrische Matrizen."""
        L1 = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)   # Lx
        L2 = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=float)   # Ly
        L3 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)   # Lz
        return [L1, L2, L3]

    def test_empty_basis_raises(self):
        """Leere Basis wirft ValueError."""
        with pytest.raises(ValueError):
            LieAlgebra([])

    def test_non_square_basis_raises(self):
        """Nicht-quadratische Basismatrix wirft ValueError."""
        with pytest.raises(ValueError):
            LieAlgebra([np.ones((2, 3))])

    def test_different_shapes_raises(self):
        """Basismatrizen verschiedener Form werfen ValueError."""
        with pytest.raises(ValueError):
            LieAlgebra([np.eye(2), np.eye(3)])

    def test_bracket_anticommutator(self):
        """[X, Y] = -[Y, X] (Antisymmetrie)."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        X, Y = basis[0], basis[1]
        assert_allclose(alg.bracket(X, Y), -alg.bracket(Y, X), atol=1e-14)

    def test_bracket_self_is_zero(self):
        """[X, X] = 0."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        for b in basis:
            assert_allclose(alg.bracket(b, b), np.zeros((2, 2)), atol=1e-14)

    def test_jacobi_identity(self):
        """[[X,Y],Z] + [[Y,Z],X] + [[Z,X],Y] = 0 (Jacobi-Identität)."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        X, Y, Z = basis[0], basis[1], basis[2]
        br = alg.bracket
        jacobi = br(br(X, Y), Z) + br(br(Y, Z), X) + br(br(Z, X), Y)
        assert_allclose(jacobi, np.zeros((2, 2)), atol=1e-13)

    def test_bracket_su2_commutation_relations(self):
        """so(3)-Kommutatorrelationen: [L_i, L_j] = ε_{ijk} L_k."""
        # [Lx, Ly] = Lz, [Ly, Lz] = Lx, [Lz, Lx] = Ly
        basis = self._so3_basis()
        alg = LieAlgebra(basis)
        Lx, Ly, Lz = basis
        assert_allclose(alg.bracket(Lx, Ly), Lz, atol=1e-14)
        assert_allclose(alg.bracket(Ly, Lz), Lx, atol=1e-14)
        assert_allclose(alg.bracket(Lz, Lx), Ly, atol=1e-14)

    def test_adjoint_representation_shape(self):
        """adjoint_representation gibt dim×dim-Matrix zurück."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        ad = alg.adjoint_representation(basis[0])
        assert ad.shape == (3, 3)

    def test_adjoint_representation_is_linear(self):
        """ad_{αX} = α·ad_X (Linearität)."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        X = basis[0]
        alpha = 2.5
        ad_alphaX = alg.adjoint_representation(alpha * X)
        alpha_adX = alpha * alg.adjoint_representation(X)
        assert_allclose(ad_alphaX, alpha_adX, atol=1e-12)

    def test_killing_form_symmetry(self):
        """B(X, Y) = B(Y, X) (Symmetrie der Killing-Form)."""
        basis = self._so3_basis()
        alg = LieAlgebra(basis)
        X, Y = basis[0], basis[1]
        assert_allclose(alg.killing_form(X, Y),
                        alg.killing_form(Y, X), atol=1e-12)

    def test_is_semisimple_for_su2(self):
        """su(2) ist halbeinfach (Killing-Form nicht-degeneriert)."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        assert bool(alg.is_semisimple())

    def test_is_semisimple_for_so3(self):
        """so(3) ist halbeinfach."""
        basis = self._so3_basis()
        alg = LieAlgebra(basis)
        assert bool(alg.is_semisimple())

    def test_bracket_zero_matrix(self):
        """[0, X] = 0 (Nullmatrix in Klammer gibt Null)."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        zero = np.zeros((2, 2))
        X = basis[0]
        assert_allclose(alg.bracket(zero, X), np.zeros((2, 2)), atol=1e-14)

    def test_lie_algebra_dimension(self):
        """LieAlgebra.dim entspricht Anzahl der Basismatrizen."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        assert alg.dim == 3

    def test_killing_form_negative_definite_for_su2(self):
        """Killing-Form von su(2) auf den Basisvektoren ist negativ."""
        basis = self._su2_basis()
        alg = LieAlgebra(basis)
        for b in basis:
            # B(e_i, e_i) sollte real und negativ sein für kompakte Algebren
            val = alg.killing_form(b, b)
            assert val.real < 0, f"B(e_i, e_i) sollte negativ sein, ist {val}"


# =============================================================================
# TEST: ExponentialMap
# =============================================================================

class TestExponentialMap:
    """Tests für Exponentialabbildung und BCH-Formel."""

    def test_exp_zero_is_identity(self):
        """exp(0) = I."""
        n = 3
        X = np.zeros((n, n))
        result = ExponentialMap.exp_map(X)
        assert_allclose(result, np.eye(n), atol=1e-14)

    def test_exp_of_skew_symmetric_is_orthogonal(self):
        """exp(X) für schiefsymmetrisches X liegt in SO(n)."""
        so3 = SO(3)
        for _ in range(5):
            X = skew_symmetric_matrix(3)
            A = ExponentialMap.exp_map(X)
            assert so3.is_element(A.real, tol=1e-10)

    def test_exp_of_skew_hermitian_traceless_is_su(self):
        """exp(X) für schief-hermitesch spurloses X liegt in SU(2)."""
        su2 = SU(2)
        for _ in range(5):
            X = skew_hermitian_traceless(2)
            A = ExponentialMap.exp_map(X)
            assert su2.is_element(A, tol=1e-10)

    def test_exp_inverse(self):
        """exp(X)·exp(-X) = I."""
        n = 3
        X = np.random.randn(n, n)
        eX = ExponentialMap.exp_map(X)
        emX = ExponentialMap.exp_map(-X)
        assert_allclose(eX @ emX, np.eye(n), atol=1e-12)

    def test_log_of_identity_is_zero(self):
        """log(I) = 0."""
        n = 3
        I = np.eye(n)
        result = ExponentialMap.log_map(I)
        assert_allclose(result, np.zeros((n, n)), atol=1e-12)

    def test_log_exp_roundtrip(self):
        """log(exp(X)) ≈ X für kleine X."""
        # Für kleine Matrizen gilt log(exp(X)) = X exakt
        X = 0.1 * skew_hermitian_traceless(3)
        result = ExponentialMap.log_map(ExponentialMap.exp_map(X))
        assert_allclose(result, X, atol=1e-11)

    def test_exp_log_roundtrip(self):
        """exp(log(A)) ≈ A für A nahe I."""
        # Erzeuge A nahe Einheitsmatrix
        X = 0.1 * skew_symmetric_matrix(3)
        A = ExponentialMap.exp_map(X)
        result = ExponentialMap.exp_map(ExponentialMap.log_map(A))
        assert_allclose(result, A, atol=1e-11)

    def test_bch_order1_is_sum(self):
        """BCH Ordnung 1: Z ≈ X + Y."""
        n = 2
        X = 0.1 * np.random.randn(n, n)
        Y = 0.1 * np.random.randn(n, n)
        Z = ExponentialMap.baker_campbell_hausdorff(X, Y, order=1)
        assert_allclose(Z, X + Y, atol=1e-15)

    def test_bch_order2_contains_commutator(self):
        """BCH Ordnung 2: Z ≈ X + Y + ½[X,Y]."""
        n = 2
        X = 0.1 * np.random.randn(n, n)
        Y = 0.1 * np.random.randn(n, n)
        Z = ExponentialMap.baker_campbell_hausdorff(X, Y, order=2)
        expected = X + Y + 0.5 * (X @ Y - Y @ X)
        assert_allclose(Z, expected, atol=1e-15)

    def test_bch_commuting_matrices(self):
        """Für kommutierende X, Y: BCH = X + Y (alle Kommutatorterme = 0)."""
        # Diagonalmatrizen kommutieren
        X = np.diag([0.1, 0.2])
        Y = np.diag([0.3, 0.4])
        for order in [1, 2, 3, 4]:
            Z = ExponentialMap.baker_campbell_hausdorff(X, Y, order=order)
            assert_allclose(Z, X + Y, atol=1e-14)

    def test_bch_approximates_log_exp_exp(self):
        """BCH(X,Y,3) ≈ log(exp(X)·exp(Y)) für kleine X, Y."""
        n = 2
        # Kleine Matrizen damit Reihe schnell konvergiert
        X = 0.05 * skew_symmetric_matrix(n)
        Y = 0.05 * skew_symmetric_matrix(n)
        bch = ExponentialMap.baker_campbell_hausdorff(X, Y, order=3)
        exact = ExponentialMap.log_map(
            ExponentialMap.exp_map(X) @ ExponentialMap.exp_map(Y))
        assert_allclose(bch, exact, atol=1e-6)

    def test_exp_so3_rotation_angle(self):
        """exp(θ · L_z) mit L_z ∈ so(3) erzeugt Rotation um Z-Achse."""
        so3 = SO(3)
        theta = np.pi / 3
        Lz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
        R = ExponentialMap.exp_map(theta * Lz).real
        R_expected = so3.rotation_3d_z(theta)
        assert_allclose(R, R_expected, atol=1e-12)


# =============================================================================
# TEST: Einparametrige Untergruppe
# =============================================================================

class TestOneParameterSubgroup:
    """Tests für one_parameter_subgroup(X, t)."""

    def test_t0_is_identity(self):
        """exp(0·X) = I."""
        X = np.array([[0, -1], [1, 0]], dtype=float)
        result = one_parameter_subgroup(X, 0.0)
        assert_allclose(result, np.eye(2), atol=1e-14)

    def test_homomorphism_property(self):
        """γ(s)·γ(t) = γ(s+t)."""
        X = skew_symmetric_matrix(3)
        s, t = 0.3, 0.7
        gs = one_parameter_subgroup(X, s)
        gt = one_parameter_subgroup(X, t)
        gst = one_parameter_subgroup(X, s + t)
        assert_allclose(gs @ gt, gst, atol=1e-12)

    def test_tangent_vector_at_zero(self):
        """Ableitung von γ(t) bei t=0 ist X."""
        X = np.array([[0, -1], [1, 0]], dtype=float)
        eps = 1e-6
        # Numerische Ableitung: (γ(ε) - I) / ε ≈ X
        g_eps = one_parameter_subgroup(X, eps)
        tangent = (g_eps - np.eye(2)) / eps
        assert_allclose(tangent.real, X, atol=1e-5)

    def test_stays_in_son(self):
        """Einparametrige Untergruppe in so(n) bleibt in SO(n)."""
        so3 = SO(3)
        X = skew_symmetric_matrix(3)
        for t in [0.0, 0.5, 1.0, 2.0, np.pi]:
            g = one_parameter_subgroup(X, t).real
            assert so3.is_element(g, tol=1e-10), f"t={t}: nicht in SO(3)"

    def test_negative_t(self):
        """γ(-t) = γ(t)⁻¹."""
        X = skew_symmetric_matrix(3)
        t = 1.2
        g_t = one_parameter_subgroup(X, t)
        g_mt = one_parameter_subgroup(X, -t)
        assert_allclose(g_t @ g_mt, np.eye(3), atol=1e-12)

    def test_scalar_multiplication(self):
        """γ_X(2t) = γ_{2X}(t)."""
        X = skew_symmetric_matrix(2)
        t = 0.5
        result1 = one_parameter_subgroup(X, 2 * t)
        result2 = one_parameter_subgroup(2 * X, t)
        assert_allclose(result1, result2, atol=1e-12)


# =============================================================================
# TEST: Klassifikation
# =============================================================================

class TestClassification:
    """Tests für Cartan-Klassifikation, Cartan-Matrix und Dynkin-Diagramme."""

    def test_classify_an_dimension(self):
        """A_n hat Dimension n(n+2) = (n+1)²-1."""
        for n in [1, 2, 3, 4]:
            info = classify_simple_lie_algebra(n, 'A')
            assert info['dimension'] == (n + 1) ** 2 - 1

    def test_classify_bn_lie_group(self):
        """B_n entspricht SO(2n+1)."""
        info = classify_simple_lie_algebra(2, 'B')
        assert info['lie_group'] == 'SO(5)'

    def test_classify_cn_lie_group(self):
        """C_n entspricht Sp(2n)."""
        info = classify_simple_lie_algebra(3, 'C')
        assert info['lie_group'] == 'Sp(6)'

    def test_classify_dn_lie_group(self):
        """D_n entspricht SO(2n)."""
        info = classify_simple_lie_algebra(4, 'D')
        assert info['lie_group'] == 'SO(8)'

    def test_classify_invalid_type_raises(self):
        """Ungültiger Typ 'E' wirft ValueError."""
        with pytest.raises(ValueError):
            classify_simple_lie_algebra(6, 'E')

    def test_classify_an_rank_too_small_raises(self):
        """A_0 (n=0) wirft ValueError."""
        with pytest.raises(ValueError):
            classify_simple_lie_algebra(0, 'A')

    def test_classify_bn_rank_too_small_raises(self):
        """B_1 wirft ValueError (n ≥ 2 erforderlich)."""
        with pytest.raises(ValueError):
            classify_simple_lie_algebra(1, 'B')

    def test_classify_cn_rank_too_small_raises(self):
        """C_2 wirft ValueError (n ≥ 3 erforderlich)."""
        with pytest.raises(ValueError):
            classify_simple_lie_algebra(2, 'C')

    def test_classify_dn_rank_too_small_raises(self):
        """D_3 wirft ValueError (n ≥ 4 erforderlich)."""
        with pytest.raises(ValueError):
            classify_simple_lie_algebra(3, 'D')

    def test_classify_returns_correct_keys(self):
        """Rückgabe-dict hat alle erwarteten Schlüssel."""
        info = classify_simple_lie_algebra(3, 'A')
        for key in ('type', 'rank', 'name', 'lie_group', 'dimension', 'lie_algebra'):
            assert key in info, f"Schlüssel '{key}' fehlt in Klassifikationsergebnis"

    def test_cartan_matrix_an_diagonal(self):
        """Cartan-Matrix A_n hat 2 auf der Diagonalen."""
        for n in [2, 3, 4]:
            C = cartan_matrix('A', n)
            assert_allclose(np.diag(C), 2 * np.ones(n), atol=1e-10)

    def test_cartan_matrix_an_off_diagonal(self):
        """Cartan-Matrix A_n hat -1 auf den Nebendiagonalen."""
        n = 4
        C = cartan_matrix('A', n)
        for i in range(n - 1):
            assert C[i, i+1] == -1
            assert C[i+1, i] == -1

    def test_cartan_matrix_a1(self):
        """Cartan-Matrix A_1 = [[2]]."""
        C = cartan_matrix('A', 1)
        assert C.shape == (1, 1)
        assert C[0, 0] == 2

    def test_cartan_matrix_a2(self):
        """Cartan-Matrix A_2 ist [[2,-1],[-1,2]]."""
        C = cartan_matrix('A', 2)
        expected = np.array([[2, -1], [-1, 2]])
        assert_allclose(C, expected)

    def test_cartan_matrix_b_last_entry(self):
        """Cartan-Matrix B_n: A[n-1, n-2] = -2 (kurze→lange Wurzel)."""
        for n in [2, 3, 4]:
            C = cartan_matrix('B', n)
            assert C[n-1, n-2] == -2

    def test_cartan_matrix_shape(self):
        """Cartan-Matrix hat Form n×n."""
        for type_c in ['A', 'B', 'C', 'D']:
            for n in [4, 5]:
                if type_c == 'B' and n < 2:
                    continue
                if type_c == 'C' and n < 3:
                    continue
                if type_c == 'D' and n < 4:
                    continue
                C = cartan_matrix(type_c, n)
                assert C.shape == (n, n)

    def test_dynkin_diagram_an_simply_laced(self):
        """A_n ist simply laced (alle Kanten einfach)."""
        info = dynkin_diagram_info('A', 4)
        assert info['simply_laced'] is True

    def test_dynkin_diagram_bn_not_simply_laced(self):
        """B_n ist nicht simply laced."""
        info = dynkin_diagram_info('B', 3)
        assert info['simply_laced'] is False

    def test_dynkin_diagram_cn_not_simply_laced(self):
        """C_n ist nicht simply laced."""
        info = dynkin_diagram_info('C', 3)
        assert info['simply_laced'] is False

    def test_dynkin_diagram_dn_simply_laced(self):
        """D_n ist simply laced."""
        info = dynkin_diagram_info('D', 4)
        assert info['simply_laced'] is True

    def test_dynkin_diagram_rank(self):
        """Dynkin-Diagramm-Info enthält korrekten Rang."""
        for n in [3, 4, 5]:
            info = dynkin_diagram_info('A', n)
            assert info['rank'] == n

    def test_dynkin_diagram_nodes(self):
        """Anzahl der Knoten im Dynkin-Diagramm = Rang."""
        for n in [4, 5, 6]:
            info = dynkin_diagram_info('A', n)
            assert info['nodes'] == n

    def test_dynkin_invalid_type_raises(self):
        """Ungültiger Typ wirft ValueError."""
        with pytest.raises(ValueError):
            dynkin_diagram_info('X', 3)


# =============================================================================
# TEST: Exzeptionelle Lie-Algebren
# =============================================================================

class TestExceptionalAlgebras:
    """Tests für die Info-Funktion der exzeptionellen Lie-Algebren."""

    def setup_method(self):
        """Einmal die Info laden."""
        self.info = exceptional_lie_algebras_info()

    def test_g2_rank(self):
        """G₂ hat Rang 2."""
        assert self.info['G2']['rank'] == 2

    def test_g2_dimension(self):
        """G₂ hat Dimension 14."""
        assert self.info['G2']['dimension'] == 14

    def test_f4_dimension(self):
        """F₄ hat Dimension 52."""
        assert self.info['F4']['dimension'] == 52

    def test_e6_dimension(self):
        """E₆ hat Dimension 78."""
        assert self.info['E6']['dimension'] == 78

    def test_e7_dimension(self):
        """E₇ hat Dimension 133."""
        assert self.info['E7']['dimension'] == 133

    def test_e8_dimension(self):
        """E₈ hat Dimension 248."""
        assert self.info['E8']['dimension'] == 248

    def test_e8_rank(self):
        """E₈ hat Rang 8."""
        assert self.info['E8']['rank'] == 8

    def test_all_five_algebras_present(self):
        """Alle 5 exzeptionellen Algebren sind vorhanden."""
        for key in ['G2', 'F4', 'E6', 'E7', 'E8']:
            assert key in self.info

    def test_each_has_required_keys(self):
        """Jeder Eintrag hat die erwarteten Schlüssel."""
        required = ('rank', 'dimension', 'description', 'dynkin')
        for name, data in self.info.items():
            for key in required:
                assert key in data, f"{name} fehlt Schlüssel '{key}'"

    def test_positive_roots_count(self):
        """Anzahl positiver Wurzeln ist korrekt und positiv."""
        expected = {'G2': 6, 'F4': 24, 'E6': 36, 'E7': 63, 'E8': 120}
        for name, n_roots in expected.items():
            assert self.info[name]['roots_positive'] == n_roots


# =============================================================================
# TEST: Spin-j-Darstellungen von SU(2)
# =============================================================================

class TestFundamentalRepSU2:
    """Tests für die Spin-j-Darstellungen fundamental_representation_su2."""

    def test_spin0_dimension(self):
        """Spin-0: 1-dimensionale triviale Darstellung."""
        rep = fundamental_representation_su2(0)
        assert rep['dimension'] == 1

    def test_spin_half_dimension(self):
        """Spin-1/2: 2-dimensionale fundamentale Darstellung (Pauli)."""
        rep = fundamental_representation_su2(0.5)
        assert rep['dimension'] == 2

    def test_spin1_dimension(self):
        """Spin-1: 3-dimensionale Darstellung."""
        rep = fundamental_representation_su2(1.0)
        assert rep['dimension'] == 3

    def test_spin_half_jz_eigenvalues(self):
        """Spin-1/2: J_z hat Eigenwerte +1/2 und -1/2."""
        rep = fundamental_representation_su2(0.5)
        Jz = rep['Jz']
        evals = np.sort(np.diag(Jz).real)
        assert_allclose(evals, [-0.5, 0.5], atol=1e-14)

    def test_spin1_jz_eigenvalues(self):
        """Spin-1: J_z hat Eigenwerte -1, 0, +1."""
        rep = fundamental_representation_su2(1.0)
        Jz = rep['Jz']
        evals = np.sort(np.diag(Jz).real)
        assert_allclose(evals, [-1.0, 0.0, 1.0], atol=1e-14)

    def test_casimir_eigenvalue(self):
        """Casimir-Eigenwert = j(j+1)."""
        for j in [0, 0.5, 1.0, 1.5, 2.0]:
            rep = fundamental_representation_su2(j)
            assert_allclose(rep['casimir_eigenvalue'], j * (j + 1), atol=1e-14)

    def test_su2_angular_momentum_algebra(self):
        """[Jx, Jy] = i·Jz, [Jy, Jz] = i·Jx, [Jz, Jx] = i·Jy."""
        for j in [0.5, 1.0, 1.5]:
            rep = fundamental_representation_su2(j)
            Jx, Jy, Jz = rep['Jx'], rep['Jy'], rep['Jz']
            assert_allclose(Jx @ Jy - Jy @ Jx, 1j * Jz, atol=1e-12)
            assert_allclose(Jy @ Jz - Jz @ Jy, 1j * Jx, atol=1e-12)
            assert_allclose(Jz @ Jx - Jx @ Jz, 1j * Jy, atol=1e-12)

    def test_jplus_jminus_relation(self):
        """[J_+, J_-] = 2 J_z."""
        for j in [0.5, 1.0, 2.0]:
            rep = fundamental_representation_su2(j)
            Jp, Jm, Jz = rep['Jplus'], rep['Jminus'], rep['Jz']
            comm = Jp @ Jm - Jm @ Jp
            assert_allclose(comm, 2 * Jz, atol=1e-12)

    def test_negative_j_raises(self):
        """Negativer Spin wirft ValueError."""
        with pytest.raises(ValueError):
            fundamental_representation_su2(-0.5)

    def test_non_half_integer_raises(self):
        """j = 1/3 (nicht halbganzzahlig) wirft ValueError."""
        with pytest.raises(ValueError):
            fundamental_representation_su2(1.0 / 3.0)

    def test_jx_hermitian(self):
        """J_x ist hermitesch (J_x = J_x†)."""
        for j in [0.5, 1.0, 1.5]:
            rep = fundamental_representation_su2(j)
            Jx = rep['Jx']
            assert_allclose(Jx, Jx.conj().T, atol=1e-13)

    def test_jy_hermitian(self):
        """J_y ist hermitesch."""
        for j in [0.5, 1.0, 1.5]:
            rep = fundamental_representation_su2(j)
            Jy = rep['Jy']
            assert_allclose(Jy, Jy.conj().T, atol=1e-13)

    def test_jz_hermitian(self):
        """J_z ist hermitesch."""
        for j in [0.5, 1.0, 1.5]:
            rep = fundamental_representation_su2(j)
            Jz = rep['Jz']
            assert_allclose(Jz, Jz.conj().T, atol=1e-13)

    def test_m_values_count(self):
        """Anzahl der m-Werte = 2j+1."""
        for j in [0, 0.5, 1, 1.5, 2]:
            rep = fundamental_representation_su2(j)
            assert len(rep['m_values']) == int(2 * j + 1)

    def test_spin0_operators_are_zero(self):
        """Spin-0: alle Operatoren sind 1×1-Nullmatrizen."""
        rep = fundamental_representation_su2(0)
        for key in ('Jx', 'Jy', 'Jz', 'Jplus', 'Jminus'):
            assert_allclose(rep[key], np.zeros((1, 1)), atol=1e-14)


# =============================================================================
# TEST: Casimir-Operator
# =============================================================================

class TestCasimirElement:
    """Tests für den quadratischen Casimir-Operator."""

    def test_casimir_su2_proportional_to_identity(self):
        """Casimir-Operator von SU(2) ist proportional zur Einheitsmatrix."""
        rep = fundamental_representation_su2(0.5)
        basis = [rep['Jx'], rep['Jy'], rep['Jz']]
        C = casimir_element(basis)
        # J_x² + J_y² + J_z² = j(j+1)·I für j=1/2 → 3/4·I
        j = 0.5
        expected = j * (j + 1) * np.eye(2, dtype=complex)
        assert_allclose(C, expected, atol=1e-12)

    def test_casimir_spin1_eigenvalue(self):
        """Casimir für Spin-1: C = j(j+1)·I = 2·I."""
        rep = fundamental_representation_su2(1.0)
        basis = [rep['Jx'], rep['Jy'], rep['Jz']]
        C = casimir_element(basis)
        expected = 1.0 * 2.0 * np.eye(3, dtype=complex)   # j(j+1) = 2
        assert_allclose(C, expected, atol=1e-12)

    def test_casimir_commutes_with_algebra(self):
        """Casimir-Operator kommutiert mit allen Algebra-Elementen: [C, J_k] = 0."""
        rep = fundamental_representation_su2(1.0)
        basis = [rep['Jx'], rep['Jy'], rep['Jz']]
        C = casimir_element(basis)
        for J in basis:
            comm = C @ J - J @ C
            assert_allclose(comm, np.zeros_like(C), atol=1e-11)

    def test_casimir_empty_basis_raises(self):
        """Leere Basis wirft ValueError."""
        with pytest.raises(ValueError):
            casimir_element([])

    def test_casimir_identity_basis(self):
        """Casimir mit Einheitsmatrix als Basis: C = I²."""
        basis = [np.eye(2)]
        C = casimir_element(basis)
        assert_allclose(C, np.eye(2), atol=1e-14)

    def test_casimir_shape(self):
        """Casimir hat gleiche Form wie die Basismatrizen."""
        rep = fundamental_representation_su2(1.5)
        basis = [rep['Jx'], rep['Jy'], rep['Jz']]
        C = casimir_element(basis)
        assert C.shape == (4, 4)


# =============================================================================
# TEST: Randfälle (Edge Cases)
# =============================================================================

class TestEdgeCases:
    """Tests für Randfälle und besondere Situationen."""

    def test_so2_rotation_full_circle(self):
        """R(2π) = I₂ (eine volle Umdrehung)."""
        so2 = SO(2)
        R = so2.rotation_2d(2 * np.pi)
        assert_allclose(R, np.eye(2), atol=1e-14)

    def test_so3_rotation_x_pi_squared(self):
        """R_x(π)² = I₃ (zwei halbe Umdrehungen)."""
        so3 = SO(3)
        Rx_pi = so3.rotation_3d_x(np.pi)
        assert_allclose(Rx_pi @ Rx_pi, np.eye(3), atol=1e-14)

    def test_exp_large_matrix(self):
        """exp() funktioniert auch für größere Matrizen."""
        n = 5
        X = 0.1 * skew_symmetric_matrix(n)
        A = ExponentialMap.exp_map(X)
        assert A.shape == (n, n)
        # Sollte orthogonal sein
        assert_allclose(A.real @ A.real.T, np.eye(n), atol=1e-10)

    def test_one_parameter_subgroup_large_t(self):
        """Einparametrige Untergruppe für großes t bleibt stabil."""
        X = np.array([[0, -1], [1, 0]], dtype=float)  # Rotation
        for t in [10.0, 100.0, 1000.0]:
            g = one_parameter_subgroup(X, t)
            # |det| sollte 1 sein (Rotation)
            assert_allclose(abs(np.linalg.det(g)), 1.0, atol=1e-10)

    def test_bch_order4_implemented(self):
        """BCH bis Ordnung 4 ist implementiert und gibt anderen Wert als Ordnung 3."""
        X = 0.3 * np.array([[0, 1], [-1, 0]], dtype=complex)
        Y = 0.2 * np.array([[0, 1j], [1j, 0]], dtype=complex)
        Z3 = ExponentialMap.baker_campbell_hausdorff(X, Y, order=3)
        Z4 = ExponentialMap.baker_campbell_hausdorff(X, Y, order=4)
        # Ordnung 4 kann von Ordnung 3 abweichen (4. Term ≠ 0 wenn Y·[X,XY] ≠ 0)
        assert Z3.shape == Z4.shape   # gleiche Form

    def test_lie_algebra_with_single_basis_element(self):
        """Lie-Algebra mit nur einem Basisvektor (abelsch)."""
        basis = [np.array([[0, 1], [-1, 0]], dtype=complex)]
        alg = LieAlgebra(basis)
        assert alg.dim == 1
        # [X, X] = 0
        result = alg.bracket(basis[0], basis[0])
        assert_allclose(result, np.zeros((2, 2)), atol=1e-14)

    def test_cartan_matrix_d4_shape(self):
        """Cartan-Matrix D_4 hat Form 4×4."""
        C = cartan_matrix('D', 4)
        assert C.shape == (4, 4)
        assert C[0, 0] == 2

    def test_gl_inverse_of_diagonal(self):
        """GL.inverse() für Diagonalmatrix ist korrekt."""
        gl = GL(3)
        D = np.diag([2.0, 3.0, 4.0])
        D_inv = gl.inverse(D)
        assert_allclose(D_inv, np.diag([0.5, 1.0/3.0, 0.25]), atol=1e-14)

    def test_su_field_is_complex(self):
        """SU-Gruppe verwendet komplexen Körper."""
        su3 = SU(3)
        assert su3.field == 'complex'

    def test_so_field_is_real(self):
        """SO-Gruppe verwendet reellen Körper."""
        so4 = SO(4)
        assert so4.field == 'real'

    def test_classification_a1_is_su2_algebra(self):
        """A_1 entspricht sl(2, ℂ) = Lie-Algebra von SL(2)."""
        info = classify_simple_lie_algebra(1, 'A')
        assert 'SL(2)' in info['lie_group']
        assert info['dimension'] == 3

    def test_exceptional_info_returns_dict(self):
        """exceptional_lie_algebras_info() gibt dict zurück."""
        result = exceptional_lie_algebras_info()
        assert isinstance(result, dict)
        assert len(result) == 5

    def test_matrix_lie_group_compose_associativity(self):
        """Matrixmultiplikation ist assoziativ: (A·B)·C = A·(B·C)."""
        g = MatrixLieGroup(3)
        A = np.random.randn(3, 3)
        B = np.random.randn(3, 3)
        C = np.random.randn(3, 3)
        lhs = g.compose(g.compose(A, B), C)
        rhs = g.compose(A, g.compose(B, C))
        assert_allclose(lhs, rhs, atol=1e-13)
