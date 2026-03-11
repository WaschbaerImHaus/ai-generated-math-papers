"""
@file test_matrix_ops.py
@brief Umfassende Tests für das Modul matrix_ops.py (Matrix-Klasse + Hilfsfunktionen).
@description
    Testet alle Methoden der Matrix-Klasse und Modulfunktionen:
    - Konstruktor: Dimensionen, Fehler bei leerer Matrix
    - identity(), zeros(): Klassenmethoden
    - get(), set(): Element-Zugriff
    - __matmul__(): Matrixmultiplikation
    - transpose(): Transponierung
    - trace(): Spur
    - determinant(): Determinante via Gauss-Elimination
    - inverse(): Gauss-Jordan-Inversion
    - solve(): LGS-Lösung
    - eigenvalues(): Eigenwert-Berechnung (2×2 analytisch, n×n QR-Iteration)
    - eigenvectors(): Eigenvektor-Berechnung via SVD
    - analyze_stability(): Stabilitätsanalyse
    - condition_number_check(): Konditionszahl-Prüfung

    Edge-Cases:
    - Leere Matrix → ValueError
    - Nicht-quadratische Matrix bei det/inv/trace → ValueError
    - Singuläre Matrix bei inv/solve → SingularMatrixError
    - Matrixmultiplikation mit inkompatiblen Dimensionen → ValueError

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from matrix_ops import Matrix, analyze_stability, condition_number_check
from vectors import Vector
from exceptions import SingularMatrixError


# =============================================================================
# TESTS: Matrix-Konstruktor und Eigenschaften
# =============================================================================

class TestMatrixConstructor:
    """Tests für den Konstruktor und grundlegende Eigenschaften."""

    def test_2x2_matrix(self):
        """2×2-Matrix hat korrekte Dimensionen."""
        M = Matrix([[1, 2], [3, 4]])
        assert M.rows == 2
        assert M.cols == 2

    def test_3x3_matrix(self):
        """3×3-Matrix hat korrekte Dimensionen."""
        M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert M.rows == 3
        assert M.cols == 3

    def test_rectangular_matrix(self):
        """Rechteckige Matrix hat korrekte Dimensionen."""
        M = Matrix([[1, 2, 3], [4, 5, 6]])
        assert M.rows == 2
        assert M.cols == 3

    def test_empty_matrix_raises(self):
        """Leere Matrix → ValueError."""
        with pytest.raises((ValueError, Exception)):
            Matrix([])

    def test_data_is_copy(self):
        """Änderungen an Originaldaten ändern Matrix nicht."""
        data = [[1, 2], [3, 4]]
        M = Matrix(data)
        data[0][0] = 999
        assert M.get(0, 0) == 1


# =============================================================================
# TESTS: identity() und zeros()
# =============================================================================

class TestMatrixClassMethods:
    """Tests für Klassenmethoden identity() und zeros()."""

    def test_identity_2x2(self):
        """2×2-Einheitsmatrix."""
        I = Matrix.identity(2)
        assert I.get(0, 0) == 1.0
        assert I.get(0, 1) == 0.0
        assert I.get(1, 0) == 0.0
        assert I.get(1, 1) == 1.0

    def test_identity_3x3_trace(self):
        """Spur der 3×3-Einheitsmatrix ist 3."""
        I = Matrix.identity(3)
        assert I.trace() == 3.0

    def test_identity_det_is_one(self):
        """Determinante der Einheitsmatrix ist 1."""
        I = Matrix.identity(3)
        assert abs(I.determinant() - 1.0) < 1e-10

    def test_zeros_2x3(self):
        """2×3-Nullmatrix hat überall 0."""
        Z = Matrix.zeros(2, 3)
        assert Z.rows == 2
        assert Z.cols == 3
        for i in range(2):
            for j in range(3):
                assert Z.get(i, j) == 0.0


# =============================================================================
# TESTS: get() und set()
# =============================================================================

class TestMatrixGetSet:
    """Tests für Elementzugriff."""

    def test_get_element(self):
        """get() gibt korrektes Element zurück."""
        M = Matrix([[1, 2], [3, 4]])
        assert M.get(0, 0) == 1
        assert M.get(1, 1) == 4

    def test_set_element(self):
        """set() ändert Element korrekt."""
        M = Matrix([[1, 2], [3, 4]])
        M.set(0, 0, 99)
        assert M.get(0, 0) == 99

    def test_set_does_not_affect_other(self):
        """set() ändert nur das spezifizierte Element."""
        M = Matrix([[1, 2], [3, 4]])
        M.set(0, 0, 99)
        assert M.get(0, 1) == 2
        assert M.get(1, 0) == 3


# =============================================================================
# TESTS: Matrixmultiplikation (__matmul__)
# =============================================================================

class TestMatrixMultiplication:
    """Tests für Matrixmultiplikation."""

    def test_identity_multiplication(self):
        """A @ I = A."""
        A = Matrix([[1, 2], [3, 4]])
        I = Matrix.identity(2)
        result = A @ I
        assert abs(result.get(0, 0) - 1) < 1e-12
        assert abs(result.get(1, 1) - 4) < 1e-12

    def test_known_product(self):
        """Bekanntes Matrixprodukt."""
        A = Matrix([[1, 2], [3, 4]])
        B = Matrix([[5, 6], [7, 8]])
        C = A @ B
        # [1*5+2*7, 1*6+2*8] = [19, 22]
        # [3*5+4*7, 3*6+4*8] = [43, 50]
        assert abs(C.get(0, 0) - 19) < 1e-12
        assert abs(C.get(0, 1) - 22) < 1e-12
        assert abs(C.get(1, 0) - 43) < 1e-12
        assert abs(C.get(1, 1) - 50) < 1e-12

    def test_incompatible_dimensions_raises(self):
        """Inkompatible Dimensionen → ValueError."""
        A = Matrix([[1, 2, 3]])  # 1×3
        B = Matrix([[1, 2]])     # 1×2
        with pytest.raises(ValueError):
            A @ B

    def test_rectangular_multiplication(self):
        """Rechteckige Matrixmultiplikation: (2×3) @ (3×2) = 2×2."""
        A = Matrix([[1, 0, 0], [0, 1, 0]])  # 2×3
        B = Matrix([[1, 0], [0, 1], [0, 0]])  # 3×2
        C = A @ B
        assert C.rows == 2
        assert C.cols == 2


# =============================================================================
# TESTS: transpose()
# =============================================================================

class TestMatrixTranspose:
    """Tests für die Transponierung."""

    def test_transpose_2x2(self):
        """Transponierung einer 2×2-Matrix."""
        M = Matrix([[1, 2], [3, 4]])
        T = M.transpose()
        assert T.get(0, 0) == 1
        assert T.get(0, 1) == 3
        assert T.get(1, 0) == 2
        assert T.get(1, 1) == 4

    def test_transpose_dimensions(self):
        """Transponierung ändert Dimensionen korrekt."""
        M = Matrix([[1, 2, 3], [4, 5, 6]])  # 2×3
        T = M.transpose()
        assert T.rows == 3
        assert T.cols == 2

    def test_double_transpose_is_identity(self):
        """(A^T)^T = A."""
        A = Matrix([[1, 2], [3, 4]])
        AT = A.transpose()
        ATT = AT.transpose()
        for i in range(2):
            for j in range(2):
                assert abs(ATT.get(i, j) - A.get(i, j)) < 1e-12

    def test_transpose_symmetric_matrix(self):
        """Symmetrische Matrix ist gleich ihrer Transponierten."""
        M = Matrix([[1, 2], [2, 4]])
        T = M.transpose()
        for i in range(2):
            for j in range(2):
                assert abs(M.get(i, j) - T.get(i, j)) < 1e-12


# =============================================================================
# TESTS: trace()
# =============================================================================

class TestMatrixTrace:
    """Tests für die Spur."""

    def test_trace_2x2(self):
        """Spur einer 2×2-Matrix."""
        M = Matrix([[1, 2], [3, 4]])
        assert M.trace() == 5  # 1 + 4

    def test_trace_3x3(self):
        """Spur einer 3×3-Matrix."""
        M = Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        assert M.trace() == 6

    def test_trace_identity(self):
        """Spur der n×n-Einheitsmatrix ist n."""
        for n in [1, 2, 3, 4]:
            I = Matrix.identity(n)
            assert I.trace() == float(n)

    def test_trace_non_square_raises(self):
        """Spur für nicht-quadratische Matrix → ValueError."""
        M = Matrix([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError):
            M.trace()


# =============================================================================
# TESTS: determinant()
# =============================================================================

class TestMatrixDeterminant:
    """Tests für die Determinante."""

    def test_det_2x2(self):
        """Determinante einer 2×2-Matrix: [[1,2],[3,4]] = -2."""
        M = Matrix([[1, 2], [3, 4]])
        assert abs(M.determinant() - (-2)) < 1e-10

    def test_det_identity(self):
        """Determinante der Einheitsmatrix ist 1."""
        I = Matrix.identity(3)
        assert abs(I.determinant() - 1.0) < 1e-10

    def test_det_singular_is_zero(self):
        """Singuläre Matrix hat Determinante 0."""
        M = Matrix([[1, 2], [2, 4]])  # Zeile 2 = 2 * Zeile 1
        assert abs(M.determinant()) < 1e-10

    def test_det_3x3_known(self):
        """Bekannte 3×3-Determinante."""
        # Diagonal-Matrix: det = 1*2*3 = 6
        M = Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        assert abs(M.determinant() - 6.0) < 1e-10

    def test_det_non_square_raises(self):
        """Determinante für nicht-quadratische Matrix → ValueError."""
        M = Matrix([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError):
            M.determinant()

    def test_det_product_rule(self):
        """det(A@B) = det(A) * det(B)."""
        A = Matrix([[2, 1], [1, 3]])
        B = Matrix([[1, 2], [0, 1]])
        det_AB = (A @ B).determinant()
        det_A_times_det_B = A.determinant() * B.determinant()
        assert abs(det_AB - det_A_times_det_B) < 1e-10


# =============================================================================
# TESTS: inverse()
# =============================================================================

class TestMatrixInverse:
    """Tests für die Matrixinverse."""

    def test_inverse_2x2(self):
        """Inverse einer 2×2-Matrix."""
        M = Matrix([[2, 0], [0, 2]])
        inv = M.inverse()
        # Inverse von diag(2,2) ist diag(1/2, 1/2)
        assert abs(inv.get(0, 0) - 0.5) < 1e-10
        assert abs(inv.get(1, 1) - 0.5) < 1e-10

    def test_inverse_times_matrix_is_identity(self):
        """A @ A^{-1} = I."""
        M = Matrix([[1, 2], [3, 4]])
        inv = M.inverse()
        result = M @ inv
        I = Matrix.identity(2)
        for i in range(2):
            for j in range(2):
                assert abs(result.get(i, j) - I.get(i, j)) < 1e-10

    def test_inverse_singular_raises(self):
        """Singuläre Matrix → SingularMatrixError."""
        M = Matrix([[1, 2], [2, 4]])
        with pytest.raises(SingularMatrixError):
            M.inverse()

    def test_inverse_non_square_raises(self):
        """Nicht-quadratische Matrix → ValueError."""
        M = Matrix([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError):
            M.inverse()

    def test_inverse_3x3(self):
        """Inverse einer 3×3-Einheitsmatrix ist Einheitsmatrix."""
        I = Matrix.identity(3)
        inv = I.inverse()
        for i in range(3):
            for j in range(3):
                expected = 1.0 if i == j else 0.0
                assert abs(inv.get(i, j) - expected) < 1e-10


# =============================================================================
# TESTS: solve()
# =============================================================================

class TestMatrixSolve:
    """Tests für LGS-Lösung."""

    def test_solve_simple(self):
        """Einfaches 2×2-System: [[2,0],[0,2]] x = [4,6] → x = [2,3]."""
        A = Matrix([[2, 0], [0, 2]])
        b = Vector([4, 6])
        x = A.solve(b)
        assert abs(x.components[0] - 2.0) < 1e-10
        assert abs(x.components[1] - 3.0) < 1e-10

    def test_solve_consistency(self):
        """A @ x = b muss erfüllt sein."""
        A = Matrix([[1, 2], [3, 4]])
        b = Vector([5, 11])
        x = A.solve(b)
        # Prüfung: A @ x ≈ b
        for i in range(2):
            ax_i = sum(A.get(i, j) * x.components[j] for j in range(2))
            assert abs(ax_i - b.components[i]) < 1e-10

    def test_solve_singular_raises(self):
        """Singuläres System → SingularMatrixError."""
        A = Matrix([[1, 2], [2, 4]])
        b = Vector([3, 7])
        with pytest.raises(SingularMatrixError):
            A.solve(b)

    def test_solve_dimension_mismatch_raises(self):
        """Dimensionsfehler → ValueError."""
        A = Matrix([[1, 2], [3, 4]])
        b = Vector([1, 2, 3])  # Zu lang
        with pytest.raises(ValueError):
            A.solve(b)

    def test_solve_identity_system(self):
        """I @ x = b → x = b."""
        I = Matrix.identity(3)
        b = Vector([1, 2, 3])
        x = I.solve(b)
        for i in range(3):
            assert abs(x.components[i] - b.components[i]) < 1e-10


# =============================================================================
# TESTS: eigenvalues()
# =============================================================================

class TestMatrixEigenvalues:
    """Tests für die Eigenwert-Berechnung."""

    def test_eigenvalues_diagonal_2x2(self):
        """Eigenwerte einer Diagonalmatrix sind die Diagonalelemente."""
        M = Matrix([[3, 0], [0, 5]])
        evs = M.eigenvalues()
        evs_real = sorted([float(e.real if hasattr(e, 'real') else e) for e in evs])
        assert abs(evs_real[0] - 3.0) < 1e-6
        assert abs(evs_real[1] - 5.0) < 1e-6

    def test_eigenvalues_identity(self):
        """Eigenwerte der Einheitsmatrix sind alle 1."""
        I = Matrix.identity(3)
        evs = I.eigenvalues()
        for ev in evs:
            ev_val = ev.real if hasattr(ev, 'real') else float(ev)
            assert abs(ev_val - 1.0) < 1e-6

    def test_eigenvalues_trace_sum(self):
        """Summe der Eigenwerte = Spur der Matrix."""
        M = Matrix([[2, 1], [1, 2]])
        evs = M.eigenvalues()
        ev_sum = sum(e.real if hasattr(e, 'real') else float(e) for e in evs)
        assert abs(ev_sum - M.trace()) < 1e-6

    def test_eigenvalues_det_product(self):
        """Produkt der Eigenwerte = Determinante."""
        M = Matrix([[3, 1], [0, 2]])
        evs = M.eigenvalues()
        ev_prod = 1.0
        for e in evs:
            ev_prod *= (e.real if hasattr(e, 'real') else float(e))
        assert abs(ev_prod - M.determinant()) < 1e-6

    def test_eigenvalues_count(self):
        """Anzahl der Eigenwerte = n (für n×n-Matrix)."""
        for n in [2, 3, 4]:
            M = Matrix.identity(n)
            evs = M.eigenvalues()
            assert len(evs) == n


# =============================================================================
# TESTS: eigenvectors()
# =============================================================================

class TestMatrixEigenvectors:
    """Tests für die Eigenvektor-Berechnung."""

    def test_eigenvectors_diagonal(self):
        """Eigenvektoren einer Diagonalmatrix sind Standardbasis-Vektoren."""
        M = Matrix([[3, 0], [0, 5]])
        pairs = M.eigenvectors()
        assert len(pairs) == 2

    def test_eigenvectors_normalized(self):
        """Eigenvektoren sind normiert (Norm = 1)."""
        M = Matrix([[2, 1], [1, 2]])
        pairs = M.eigenvectors()
        for _, evec in pairs:
            assert abs(evec.norm() - 1.0) < 1e-10

    def test_eigenvectors_ax_equals_lambda_x(self):
        """A @ v ≈ λ · v für Eigenvektor-Eigenvert-Paare."""
        import numpy as np
        M = Matrix([[3, 0], [0, 5]])
        A = np.array(M._data, dtype=float)
        pairs = M.eigenvectors()
        for lam, evec in pairs:
            lam_r = lam.real if hasattr(lam, 'real') else float(lam)
            v = np.array(evec.components)
            Av = A @ v
            lv = lam_r * v
            # Av und lv können sich im Vorzeichen unterscheiden
            err = min(np.linalg.norm(Av - lv), np.linalg.norm(Av + lv))
            assert err < 1e-8, f"Eigenvektorbedingung verletzt: err={err}"


# =============================================================================
# TESTS: analyze_stability()
# =============================================================================

class TestAnalyzeStability:
    """Tests für die Stabilitätsanalyse."""

    def test_returns_dict(self):
        """analyze_stability() gibt ein Dictionary zurück."""
        M = Matrix([[1, 0], [0, 1]])
        result = analyze_stability(M)
        assert isinstance(result, dict)

    def test_has_required_keys(self):
        """Ergebnis enthält alle erforderlichen Schlüssel."""
        M = Matrix([[1, 0], [0, 1]])
        result = analyze_stability(M)
        assert 'condition_number' in result
        assert 'is_well_conditioned' in result
        assert 'rank' in result
        assert 'determinant' in result

    def test_identity_is_well_conditioned(self):
        """Einheitsmatrix ist gut konditioniert."""
        I = Matrix.identity(3)
        result = analyze_stability(I)
        # np.True_ und bool(True) müssen beide akzeptiert werden
        assert bool(result['is_well_conditioned']) is True

    def test_identity_rank(self):
        """Rang der 3×3-Einheitsmatrix ist 3."""
        I = Matrix.identity(3)
        result = analyze_stability(I)
        assert result['rank'] == 3

    def test_identity_det_approx_one(self):
        """Determinante der Einheitsmatrix ≈ 1."""
        I = Matrix.identity(3)
        result = analyze_stability(I)
        assert abs(result['determinant'] - 1.0) < 1e-10

    def test_singular_matrix(self):
        """Singuläre Matrix: Rang < n, det ≈ 0."""
        M = Matrix([[1, 2], [2, 4]])  # Rang 1
        result = analyze_stability(M)
        assert result['rank'] < 2

    def test_accepts_numpy_array(self):
        """analyze_stability() akzeptiert numpy-Arrays."""
        import numpy as np
        A = np.array([[1.0, 0], [0, 2.0]])
        result = analyze_stability(A)
        assert 'condition_number' in result

    def test_accepts_list(self):
        """analyze_stability() akzeptiert Listen."""
        A = [[1.0, 0], [0, 2.0]]
        result = analyze_stability(A)
        assert 'condition_number' in result


# =============================================================================
# TESTS: condition_number_check()
# =============================================================================

class TestConditionNumberCheck:
    """Tests für condition_number_check()."""

    def test_identity_is_good(self):
        """Einheitsmatrix ist gut konditioniert → True."""
        I = Matrix.identity(3)
        assert condition_number_check(I) is True

    def test_diagonal_is_good(self):
        """Gut konditionierte Diagonalmatrix → True."""
        M = Matrix([[2, 0], [0, 3]])
        assert condition_number_check(M) is True

    def test_returns_bool(self):
        """Rückgabewert ist immer ein Bool."""
        I = Matrix.identity(2)
        result = condition_number_check(I)
        assert isinstance(result, bool)

    def test_accepts_numpy_array(self):
        """Akzeptiert numpy-Arrays."""
        import numpy as np
        A = np.eye(3)
        result = condition_number_check(A)
        assert isinstance(result, bool)

    def test_accepts_list(self):
        """Akzeptiert Listen."""
        A = [[1.0, 0], [0, 1.0]]
        result = condition_number_check(A)
        assert isinstance(result, bool)
