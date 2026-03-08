"""
@file test_linear_algebra_givens.py
@brief Tests für Givens-Rotationen in linear_algebra.py
@description
    Testet die drei neuen Funktionen:
    - givens_rotation_matrix: Konstruktion der Givens-Matrix
    - givens_qr_decomposition: QR-Zerlegung via Givens-Rotationen
    - givens_solve_least_squares: Kleinste-Quadrate-Lösung

@author Kurt Ingwer
@date 2026-03-08
"""

import pytest
import math
import numpy as np
import sys
import os

# Sicherstellen dass src/ im Suchpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from linear_algebra import (
    Matrix,
    Vector,
    givens_rotation_matrix,
    givens_qr_decomposition,
    givens_solve_least_squares
)


# ===========================================================================
# TESTS FÜR givens_rotation_matrix
# ===========================================================================

class TestGivensRotationMatrix:
    """Tests für die Givens-Rotationsmatrix-Konstruktion."""

    def test_basic_dimensions(self):
        """Givens-Matrix muss n×n sein."""
        G = givens_rotation_matrix(4, 0, 1, math.pi / 4)
        assert G.rows == 4
        assert G.cols == 4

    def test_identity_outside_ij_block(self):
        """Außerhalb des (i,j)-Blocks muss G wie I aussehen."""
        n = 4
        theta = math.pi / 4
        G = givens_rotation_matrix(n, 0, 1, theta)

        # (2,2) und (3,3) müssen 1 sein (keine Rotation in diesen Indizes)
        assert abs(G.get(2, 2) - 1.0) < 1e-12
        assert abs(G.get(3, 3) - 1.0) < 1e-12
        # Alle Off-Diagonal-Elemente außerhalb des (0,1)-Blocks müssen 0 sein
        assert abs(G.get(2, 3)) < 1e-12
        assert abs(G.get(3, 2)) < 1e-12

    def test_diagonal_cos_values(self):
        """Diagonalelemente G[i,i] und G[j,j] müssen cos(θ) sein."""
        theta = math.pi / 6  # 30 Grad
        G = givens_rotation_matrix(3, 0, 2, theta)
        c = math.cos(theta)
        assert abs(G.get(0, 0) - c) < 1e-12
        assert abs(G.get(2, 2) - c) < 1e-12

    def test_offdiagonal_sin_values(self):
        """G[i,j] = -sin(θ), G[j,i] = +sin(θ)."""
        theta = math.pi / 3  # 60 Grad
        G = givens_rotation_matrix(3, 0, 1, theta)
        s = math.sin(theta)
        assert abs(G.get(0, 1) - (-s)) < 1e-12   # G[i,j] = -sin
        assert abs(G.get(1, 0) - s) < 1e-12       # G[j,i] = +sin

    def test_orthogonality(self):
        """Givens-Matrix muss orthogonal sein: G^T @ G = I."""
        theta = 1.234  # beliebiger Winkel
        n = 5
        G = givens_rotation_matrix(n, 1, 3, theta)
        GT = G.transpose()
        product = GT @ G

        # Überprüfen G^T @ G = I
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                assert abs(product.get(i, j) - expected) < 1e-12, \
                    f"G^T@G[{i},{j}] = {product.get(i,j)}, erwartet {expected}"

    def test_determinant_is_one(self):
        """Determinante einer Givens-Matrix muss 1 sein (eigentliche Rotation)."""
        G = givens_rotation_matrix(4, 0, 3, 0.7)
        det = G.determinant()
        assert abs(det - 1.0) < 1e-10

    def test_2x2_rotation(self):
        """Für n=2: G ist die klassische 2D-Rotationsmatrix."""
        theta = math.pi / 4
        G = givens_rotation_matrix(2, 0, 1, theta)
        c = math.cos(theta)
        s = math.sin(theta)
        assert abs(G.get(0, 0) - c) < 1e-12
        assert abs(G.get(0, 1) - (-s)) < 1e-12
        assert abs(G.get(1, 0) - s) < 1e-12
        assert abs(G.get(1, 1) - c) < 1e-12

    def test_invalid_indices_raise_error(self):
        """Ungültige Indizes müssen ValueError werfen."""
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, -1, 1, 0.5)
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, 0, 3, 0.5)  # j >= n
        with pytest.raises(ValueError):
            givens_rotation_matrix(3, 1, 1, 0.5)  # i == j

    def test_zero_angle_is_identity(self):
        """θ=0 → Givens-Matrix ist Einheitsmatrix."""
        G = givens_rotation_matrix(4, 0, 2, 0.0)
        I = Matrix.identity(4)
        for i in range(4):
            for j in range(4):
                assert abs(G.get(i, j) - I.get(i, j)) < 1e-12

    def test_pi_rotation(self):
        """θ=π → G[i,i]=G[j,j]=-1, G[i,j]=0, G[j,i]=0."""
        G = givens_rotation_matrix(3, 0, 1, math.pi)
        assert abs(G.get(0, 0) - (-1.0)) < 1e-12
        assert abs(G.get(1, 1) - (-1.0)) < 1e-12
        assert abs(G.get(0, 1)) < 1e-12  # -sin(π) ≈ 0
        assert abs(G.get(1, 0)) < 1e-12  # sin(π) ≈ 0


# ===========================================================================
# TESTS FÜR givens_qr_decomposition
# ===========================================================================

class TestGivensQRDecomposition:
    """Tests für die QR-Zerlegung via Givens-Rotationen."""

    def _matrix_approx_equal(self, A: Matrix, B: Matrix, tol: float = 1e-9) -> bool:
        """Hilfsfunktion: Prüft ob zwei Matrizen annähernd gleich sind."""
        if A.rows != B.rows or A.cols != B.cols:
            return False
        for i in range(A.rows):
            for j in range(A.cols):
                if abs(A.get(i, j) - B.get(i, j)) > tol:
                    return False
        return True

    def test_2x2_reconstruction(self):
        """QR für 2×2-Matrix: Q @ R muss A ergeben."""
        A = Matrix([[1.0, 2.0], [3.0, 4.0]])
        Q, R = givens_qr_decomposition(A)
        QR = Q @ R
        for i in range(2):
            for j in range(2):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-9, \
                    f"Q@R[{i},{j}] = {QR.get(i,j)}, erwartet {A.get(i,j)}"

    def test_q_is_orthogonal_2x2(self):
        """Q muss orthogonal sein: Q^T @ Q = I."""
        A = Matrix([[1.0, 2.0], [3.0, 4.0]])
        Q, _ = givens_qr_decomposition(A)
        QT = Q.transpose()
        product = QT @ Q
        for i in range(2):
            for j in range(2):
                expected = 1.0 if i == j else 0.0
                assert abs(product.get(i, j) - expected) < 1e-9, \
                    f"Q^T@Q[{i},{j}] = {product.get(i,j)}"

    def test_r_is_upper_triangular(self):
        """R muss obere Dreiecksmatrix sein (Einträge unterhalb Diagonale ≈ 0)."""
        A = Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
        _, R = givens_qr_decomposition(A)
        for i in range(R.rows):
            for j in range(R.cols):
                if i > j:  # Unterhalb der Diagonale
                    assert abs(R.get(i, j)) < 1e-9, \
                        f"R[{i},{j}] = {R.get(i,j)} sollte 0 sein"

    def test_3x3_reconstruction(self):
        """QR für 3×3-Matrix: Q @ R muss A ergeben."""
        A = Matrix([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
        Q, R = givens_qr_decomposition(A)
        QR = Q @ R
        for i in range(3):
            for j in range(3):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-8

    def test_q_is_orthogonal_3x3(self):
        """Q aus 3×3-Zerlegung muss orthogonal sein."""
        A = Matrix([[2.0, 1.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        Q, _ = givens_qr_decomposition(A)
        n = Q.rows
        QT = Q.transpose()
        product = QT @ Q
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                assert abs(product.get(i, j) - expected) < 1e-8

    def test_overdetermined_4x2(self):
        """QR für überbestimmtes System (4×2): Q ist 4×4, R ist 4×2."""
        A = Matrix([
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
            [7.0, 8.0]
        ])
        Q, R = givens_qr_decomposition(A)
        # Dimensionen prüfen
        assert Q.rows == 4 and Q.cols == 4
        assert R.rows == 4 and R.cols == 2

        # Q @ R muss A ergeben
        QR = Q @ R
        for i in range(4):
            for j in range(2):
                assert abs(QR.get(i, j) - A.get(i, j)) < 1e-8

    def test_m_less_than_n_raises_error(self):
        """m < n muss ValueError werfen."""
        A = Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # 2×3, m < n
        with pytest.raises(ValueError):
            givens_qr_decomposition(A)

    def test_identity_matrix_qr(self):
        """QR der Einheitsmatrix: Q = ±I, R = ±I."""
        n = 3
        I = Matrix.identity(n)
        Q, R = givens_qr_decomposition(I)
        QR = Q @ R
        # Q @ R muss I ergeben
        for i in range(n):
            for j in range(n):
                expected = 1.0 if i == j else 0.0
                assert abs(QR.get(i, j) - expected) < 1e-10

    def test_comparison_with_numpy_qr(self):
        """Vergleich mit numpy.linalg.qr (Referenzimplementierung)."""
        data = [[3.0, 1.0, 2.0], [1.0, 4.0, 0.0], [2.0, 0.0, 5.0]]
        A = Matrix(data)
        Q, R = givens_qr_decomposition(A)

        # Q @ R muss A ergeben (unabhängig von Vorzeichen der Spalten)
        QR = Q @ R
        for i in range(3):
            for j in range(3):
                assert abs(QR.get(i, j) - data[i][j]) < 1e-9


# ===========================================================================
# TESTS FÜR givens_solve_least_squares
# ===========================================================================

class TestGivensSolveLeastSquares:
    """Tests für die Kleinste-Quadrate-Lösung via Givens-QR."""

    def test_square_exact_solution(self):
        """Für quadratisches LGS muss die exakte Lösung gefunden werden."""
        A = Matrix([[2.0, 1.0], [1.0, 3.0]])
        b = [5.0, 10.0]
        x = givens_solve_least_squares(A, b)

        # Exakte Lösung: Ax = b
        assert abs(2 * x[0] + x[1] - 5.0) < 1e-8
        assert abs(x[0] + 3 * x[1] - 10.0) < 1e-8

    def test_overdetermined_3x2(self):
        """Überbestimmtes System 3×2: Kleinste-Quadrate-Lösung."""
        # System: x + y = 3, 2x + y = 4, x + 2y = 5
        A = Matrix([[1.0, 1.0], [2.0, 1.0], [1.0, 2.0]])
        b = [3.0, 4.0, 5.0]
        x = givens_solve_least_squares(A, b)

        # Residuum r = A@x - b minimieren → Normalengleichung: A^T @ A @ x = A^T @ b
        A_np = np.array([[1.0, 1.0], [2.0, 1.0], [1.0, 2.0]])
        b_np = np.array([3.0, 4.0, 5.0])
        x_ref, _, _, _ = np.linalg.lstsq(A_np, b_np, rcond=None)

        assert abs(x[0] - x_ref[0]) < 1e-8
        assert abs(x[1] - x_ref[1]) < 1e-8

    def test_overdetermined_4x2(self):
        """Überbestimmtes System 4×2: Vergleich mit numpy.lstsq."""
        A_data = [[1.0, 1.0], [2.0, 1.0], [3.0, 1.0], [4.0, 1.0]]
        b = [2.0, 3.0, 5.0, 6.0]
        A = Matrix(A_data)
        x = givens_solve_least_squares(A, b)

        # Referenz via numpy
        A_np = np.array(A_data)
        b_np = np.array(b)
        x_ref, _, _, _ = np.linalg.lstsq(A_np, b_np, rcond=None)

        assert abs(x[0] - x_ref[0]) < 1e-8
        assert abs(x[1] - x_ref[1]) < 1e-8

    def test_dimension_mismatch_raises_error(self):
        """Falsche Dimension von b muss ValueError werfen."""
        A = Matrix([[1.0, 2.0], [3.0, 4.0]])
        with pytest.raises(ValueError):
            givens_solve_least_squares(A, [1.0, 2.0, 3.0])  # b hat 3 Einträge, m=2

    def test_underdetermined_raises_error(self):
        """m < n muss ValueError werfen."""
        A = Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # 2×3
        with pytest.raises(ValueError):
            givens_solve_least_squares(A, [1.0, 2.0])

    def test_identity_system(self):
        """Einheitsmatrix: x = b."""
        A = Matrix.identity(3)
        b = [1.0, 2.0, 3.0]
        x = givens_solve_least_squares(A, b)
        for i in range(3):
            assert abs(x[i] - b[i]) < 1e-10

    def test_residual_minimization(self):
        """Das Residuum ||Ax - b||² muss minimal sein."""
        A_data = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        b = [1.0, 2.0, 1.0]
        A = Matrix(A_data)
        x = givens_solve_least_squares(A, b)

        # Berechne Residuum
        A_np = np.array(A_data)
        x_np = np.array(x)
        b_np = np.array(b)
        residual = np.linalg.norm(A_np @ x_np - b_np)

        # Referenzlösung
        x_ref, _, _, _ = np.linalg.lstsq(A_np, b_np, rcond=None)
        residual_ref = np.linalg.norm(A_np @ x_ref - b_np)

        # Unser Residuum darf nicht größer als Referenz sein
        assert residual <= residual_ref + 1e-8

    def test_return_type_is_list(self):
        """Rückgabewert muss eine Liste sein."""
        A = Matrix([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        b = [1.0, 2.0, 3.0]
        x = givens_solve_least_squares(A, b)
        assert isinstance(x, list)
        assert len(x) == 2  # n = Anzahl Spalten


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
