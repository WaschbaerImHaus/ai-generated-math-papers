"""
@file test_matrix_decomp.py
@brief Umfassende Tests für das Modul matrix_decomp.py.
@description
    Testet alle Funktionen des Moduls:
    - lu_decomposition(): LU-Zerlegung mit Teilpivotisierung
    - qr_decomposition(): QR-Zerlegung via Householder-Reflexionen
    - svd(): Singulärwertzerlegung
    - matrix_rank(): Rang via SVD
    - condition_number(): Konditionszahl
    - givens_rotation_matrix(): Givens-Rotationsmatrix
    - givens_qr_decomposition(): QR via Givens-Rotationen
    - givens_solve_least_squares(): Kleinste-Quadrate via Givens-QR

    Mathematische Verifikation:
    - LU: P·A = L·U, L·U·x = P·b (Konsistenz)
    - QR: Q·R = A, Q^T·Q = I (Orthogonalität)
    - SVD: U·diag(σ)·V^T = A
    - Givens: G^T·G = I (Orthogonalität)

    Edge-Cases:
    - Nicht-quadratische Matrix bei LU → ValueError
    - Singuläre Matrix bei LU → SingularMatrixError
    - m < n bei QR → ValueError
    - Singuläre Matrix bei condition_number → SingularMatrixError
    - Ungültige Indizes bei givens_rotation_matrix → ValueError

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from matrix_decomp import (
    lu_decomposition,
    qr_decomposition,
    svd,
    matrix_rank,
    condition_number,
    givens_rotation_matrix,
    givens_qr_decomposition,
    givens_solve_least_squares,
)
from matrix_ops import Matrix
from vectors import Vector
from exceptions import SingularMatrixError


# =============================================================================
# TESTS: lu_decomposition()
# =============================================================================

class TestLUDecomposition:
    """Tests für die LU-Zerlegung."""

    def test_2x2_product_pa_equals_lu(self):
        """P·A = L·U muss gelten."""
        A = Matrix([[2, 1], [4, 3]])
        L, U, P, _ = lu_decomposition(A)
        # P·A
        PA = P @ A
        # L·U
        LU = L @ U
        for i in range(2):
            for j in range(2):
                assert abs(PA.get(i, j) - LU.get(i, j)) < 1e-10

    def test_3x3_product(self):
        """3×3-LU: P·A = L·U."""
        A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        L, U, P, _ = lu_decomposition(A)
        PA = P @ A
        LU = L @ U
        for i in range(3):
            for j in range(3):
                assert abs(PA.get(i, j) - LU.get(i, j)) < 1e-8

    def test_l_lower_triangular(self):
        """L ist untere Dreiecksmatrix mit Diagonale 1."""
        A = Matrix([[2, 1], [4, 3]])
        L, _, _, _ = lu_decomposition(A)
        n = L.rows
        for i in range(n):
            assert abs(L.get(i, i) - 1.0) < 1e-10  # Diagonale = 1
            for j in range(i + 1, n):
                assert abs(L.get(i, j)) < 1e-10  # Oberer Teil = 0

    def test_u_upper_triangular(self):
        """U ist obere Dreiecksmatrix."""
        A = Matrix([[2, 1], [4, 3]])
        _, U, _, _ = lu_decomposition(A)
        n = U.rows
        for i in range(n):
            for j in range(i):
                assert abs(U.get(i, j)) < 1e-10  # Unterer Teil = 0

    def test_returns_4_tuple(self):
        """LU gibt (L, U, P, n_swaps) zurück."""
        A = Matrix([[1, 2], [3, 4]])
        result = lu_decomposition(A)
        assert len(result) == 4

    def test_non_square_raises(self):
        """Nicht-quadratische Matrix → ValueError."""
        A = Matrix([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError):
            lu_decomposition(A)

    def test_singular_raises(self):
        """Singuläre Matrix → SingularMatrixError."""
        A = Matrix([[1, 2], [2, 4]])
        with pytest.raises(SingularMatrixError):
            lu_decomposition(A)

    def test_identity_lu(self):
        """LU-Zerlegung der Einheitsmatrix: L = U = P = I."""
        I = Matrix.identity(3)
        L, U, P, n_swaps = lu_decomposition(I)
        # U muss Einheitsmatrix sein
        for i in range(3):
            assert abs(U.get(i, i) - 1.0) < 1e-10


# =============================================================================
# TESTS: qr_decomposition()
# =============================================================================

class TestQRDecomposition:
    """Tests für die QR-Zerlegung via Householder."""

    def test_product_qr_equals_a(self):
        """Q · R = A."""
        A = Matrix([[1, 2], [3, 4]])
        Q, R = qr_decomposition(A)
        QR = Q @ R
        for i in range(2):
            for j in range(2):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-8

    def test_q_is_orthogonal(self):
        """Q^T · Q = I (Orthogonalität)."""
        A = Matrix([[1, 2], [3, 4]])
        Q, _ = qr_decomposition(A)
        QT = Q.transpose()
        I = QT @ Q
        n = Q.rows
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                assert abs(I.get(i, j) - expected) < 1e-8

    def test_r_upper_triangular(self):
        """R ist obere Dreiecksmatrix."""
        A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        _, R = qr_decomposition(A)
        for i in range(1, R.rows):
            for j in range(min(i, R.cols)):
                assert abs(R.get(i, j)) < 1e-8, f"R[{i},{j}] = {R.get(i, j)}"

    def test_3x3_product(self):
        """3×3-QR: Q·R = A."""
        A = Matrix([[12, -51, 4], [6, 167, -68], [-4, 24, -41]])
        Q, R = qr_decomposition(A)
        QR = Q @ R
        for i in range(3):
            for j in range(3):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-6

    def test_rectangular_matrix(self):
        """Rechteckige Matrix (m > n): QR-Zerlegung."""
        A = Matrix([[1, 2], [3, 4], [5, 6]])  # 3×2
        Q, R = qr_decomposition(A)
        QR = Q @ R
        for i in range(3):
            for j in range(2):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-8

    def test_m_less_than_n_raises(self):
        """m < n → ValueError."""
        A = Matrix([[1, 2, 3]])  # 1×3
        with pytest.raises(ValueError):
            qr_decomposition(A)


# =============================================================================
# TESTS: svd()
# =============================================================================

class TestSVD:
    """Tests für die Singulärwertzerlegung."""

    def test_reconstruction(self):
        """U · diag(σ) · V^T = A."""
        A = Matrix([[1, 2], [3, 4]])
        U, sigma, Vt = svd(A)
        # Rekonstruktion
        m, n = A.rows, A.cols
        S = Matrix.zeros(m, n)
        for i in range(min(m, n)):
            S.set(i, i, sigma[i])
        A_rec = U @ S @ Vt
        for i in range(m):
            for j in range(n):
                assert abs(A_rec.get(i, j) - A.get(i, j)) < 1e-8

    def test_singular_values_nonnegative(self):
        """Singulärwerte sind immer ≥ 0."""
        A = Matrix([[1, -2], [3, 4]])
        _, sigma, _ = svd(A)
        for s in sigma:
            assert s >= 0

    def test_singular_values_sorted(self):
        """Singulärwerte sind absteigend sortiert: σ₁ ≥ σ₂ ≥ ..."""
        A = Matrix([[3, 1, 0], [0, 2, 1]])
        _, sigma, _ = svd(A)
        for i in range(len(sigma) - 1):
            assert sigma[i] >= sigma[i + 1] - 1e-10

    def test_identity_singular_values(self):
        """Singulärwerte der Einheitsmatrix sind alle 1."""
        I = Matrix.identity(3)
        _, sigma, _ = svd(I)
        for s in sigma:
            assert abs(s - 1.0) < 1e-10

    def test_u_is_orthogonal(self):
        """U ist orthogonal: U^T·U = I."""
        A = Matrix([[1, 2], [3, 4]])
        U, _, _ = svd(A)
        UT = U.transpose()
        I = UT @ U
        for i in range(2):
            for j in range(2):
                expected = 1.0 if i == j else 0.0
                assert abs(I.get(i, j) - expected) < 1e-8


# =============================================================================
# TESTS: matrix_rank()
# =============================================================================

class TestMatrixRank:
    """Tests für den Matrixrang."""

    def test_rank_identity(self):
        """Rang der n×n-Einheitsmatrix ist n."""
        for n in [2, 3, 4]:
            I = Matrix.identity(n)
            assert matrix_rank(I) == n

    def test_rank_singular(self):
        """Singuläre Matrix hat Rang < n."""
        M = Matrix([[1, 2], [2, 4]])  # Zeile 2 = 2 * Zeile 1
        r = matrix_rank(M)
        assert r == 1

    def test_rank_full_rank(self):
        """Matrix mit vollem Rang."""
        A = Matrix([[1, 0], [0, 2]])
        assert matrix_rank(A) == 2

    def test_rank_zero_matrix(self):
        """Nullmatrix hat Rang 0."""
        Z = Matrix.zeros(3, 3)
        assert matrix_rank(Z) == 0

    def test_rank_rectangular(self):
        """Rang einer 3×2-Matrix (voller Spaltenrang)."""
        A = Matrix([[1, 0], [0, 1], [0, 0]])
        r = matrix_rank(A)
        assert r == 2


# =============================================================================
# TESTS: condition_number()
# =============================================================================

class TestConditionNumber:
    """Tests für die Konditionszahl."""

    def test_identity_condition_is_one(self):
        """Konditionszahl der Einheitsmatrix ist 1."""
        I = Matrix.identity(3)
        kappa = condition_number(I)
        assert abs(kappa - 1.0) < 1e-10

    def test_condition_number_positive(self):
        """Konditionszahl ist immer ≥ 1."""
        A = Matrix([[2, 1], [1, 3]])
        kappa = condition_number(A)
        assert kappa >= 1.0 - 1e-10

    def test_singular_raises(self):
        """Singuläre Matrix → SingularMatrixError."""
        M = Matrix([[1, 2], [2, 4]])
        with pytest.raises(SingularMatrixError):
            condition_number(M)

    def test_diagonal_matrix_condition(self):
        """Konditionszahl einer Diagonalmatrix: σ_max / σ_min."""
        # diag(3, 1): κ = 3/1 = 3
        D = Matrix([[3, 0], [0, 1]])
        kappa = condition_number(D)
        assert abs(kappa - 3.0) < 1e-10


# =============================================================================
# TESTS: givens_rotation_matrix()
# =============================================================================

class TestGivensRotationMatrix:
    """Tests für die Givens-Rotationsmatrix."""

    def test_is_orthogonal(self):
        """Givens-Matrix ist orthogonal: G^T · G = I."""
        G = givens_rotation_matrix(3, 0, 1, math.pi / 4)
        GT = G.transpose()
        I = GT @ G
        for i in range(3):
            for j in range(3):
                expected = 1.0 if i == j else 0.0
                assert abs(I.get(i, j) - expected) < 1e-10

    def test_det_is_one(self):
        """Determinante der Givens-Matrix ist 1."""
        G = givens_rotation_matrix(2, 0, 1, math.pi / 6)
        det = G.determinant()
        assert abs(det - 1.0) < 1e-10

    def test_rotation_entries(self):
        """Korrekte Einträge: G[i,i]=cos, G[j,j]=cos, G[i,j]=-sin, G[j,i]=sin."""
        theta = math.pi / 3
        G = givens_rotation_matrix(3, 0, 1, theta)
        c = math.cos(theta)
        s = math.sin(theta)
        assert abs(G.get(0, 0) - c) < 1e-10
        assert abs(G.get(1, 1) - c) < 1e-10
        assert abs(G.get(0, 1) - (-s)) < 1e-10
        assert abs(G.get(1, 0) - s) < 1e-10

    def test_invalid_indices_raises(self):
        """Ungültige Indizes → ValueError."""
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, 0, 5, 0.0)  # j=5 ≥ n=3
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, -1, 1, 0.0)  # i=-1 < 0

    def test_same_index_raises(self):
        """i == j → ValueError."""
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, 1, 1, 0.0)

    def test_dimension(self):
        """Ergebnismatrix hat korrekte Dimension."""
        G = givens_rotation_matrix(5, 1, 3, 0.5)
        assert G.rows == 5
        assert G.cols == 5


# =============================================================================
# TESTS: givens_qr_decomposition()
# =============================================================================

class TestGivensQRDecomposition:
    """Tests für die Givens-QR-Zerlegung."""

    def test_product_qr_equals_a(self):
        """Q · R = A."""
        A = Matrix([[1, 2], [3, 4]])
        Q, R = givens_qr_decomposition(A)
        QR = Q @ R
        for i in range(2):
            for j in range(2):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-8

    def test_q_is_orthogonal(self):
        """Q^T · Q = I."""
        A = Matrix([[1, 2], [3, 4]])
        Q, _ = givens_qr_decomposition(A)
        QT = Q.transpose()
        I_approx = QT @ Q
        n = Q.rows
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                assert abs(I_approx.get(i, j) - expected) < 1e-8

    def test_r_upper_triangular(self):
        """R ist obere Dreiecksmatrix (Unterdiagonale ≈ 0)."""
        A = Matrix([[1, 2], [3, 4]])
        _, R = givens_qr_decomposition(A)
        assert abs(R.get(1, 0)) < 1e-8

    def test_3x3(self):
        """3×3-Givens-QR: Q·R = A."""
        A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        Q, R = givens_qr_decomposition(A)
        QR = Q @ R
        for i in range(3):
            for j in range(3):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-7

    def test_m_less_than_n_raises(self):
        """m < n → ValueError."""
        A = Matrix([[1, 2, 3]])  # 1×3
        with pytest.raises(ValueError):
            givens_qr_decomposition(A)


# =============================================================================
# TESTS: givens_solve_least_squares()
# =============================================================================

class TestGivensSolveLeastSquares:
    """Tests für die Kleinste-Quadrate-Lösung via Givens-QR."""

    def test_exact_solution_square(self):
        """Quadratisches System: exakte Lösung."""
        # [[2,0],[0,3]] x = [4, 9] → x = [2, 3]
        A = Matrix([[2, 0], [0, 3]])
        b = [4.0, 9.0]
        x = givens_solve_least_squares(A, b)
        assert abs(x[0] - 2.0) < 1e-8
        assert abs(x[1] - 3.0) < 1e-8

    def test_overdetermined_consistency(self):
        """Überbestimmtes System: A·x ≈ b (im Sinn der kleinsten Quadrate)."""
        A = Matrix([[1, 0], [0, 1], [0, 0]])  # 3×2
        b = [1.0, 2.0, 0.5]  # 0.5 ist Fehler (keine exakte Lösung)
        x = givens_solve_least_squares(A, b)
        assert len(x) == 2
        # Beste Lösung: x ≈ [1, 2]
        assert abs(x[0] - 1.0) < 1e-8
        assert abs(x[1] - 2.0) < 1e-8

    def test_dimension_mismatch_raises(self):
        """Falsches b → ValueError."""
        A = Matrix([[1, 2], [3, 4]])
        b = [1.0, 2.0, 3.0]  # Zu lang
        with pytest.raises(ValueError):
            givens_solve_least_squares(A, b)

    def test_underdetermined_raises(self):
        """m < n → ValueError."""
        A = Matrix([[1, 2, 3]])  # 1×3
        b = [1.0]
        with pytest.raises(ValueError):
            givens_solve_least_squares(A, b)

    def test_residual_minimized(self):
        """Residuum ||A·x - b||² ist minimal."""
        import numpy as np
        A_m = Matrix([[1, 1], [1, 2], [1, 3]])
        b = [6.0, 5.0, 7.0]
        x = givens_solve_least_squares(A_m, b)
        A_np = np.array([[1, 1], [1, 2], [1, 3]], dtype=float)
        b_np = np.array(b, dtype=float)
        x_np = np.array(x, dtype=float)
        residual = np.linalg.norm(A_np @ x_np - b_np)
        # Numpy-Referenzlösung
        x_ref, _, _, _ = np.linalg.lstsq(A_np, b_np, rcond=None)
        residual_ref = np.linalg.norm(A_np @ x_ref - b_np)
        # Residuen sollten vergleichbar sein
        assert abs(residual - residual_ref) < 1e-6
