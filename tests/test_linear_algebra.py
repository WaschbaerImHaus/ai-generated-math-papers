"""
@file test_linear_algebra.py
@brief Tests für das Lineare-Algebra-Modul (TDD).
@description
    Testet Matrizen- und Vektoroperationen:
    - Matrixmultiplikation, Transponierung, Spur
    - Determinante (Laplace-Entwicklung und LU-Zerlegung)
    - Gauss-Jordan-Elimination
    - Eigenwerte und Eigenvektoren
    - Gram-Schmidt-Orthogonalisierung
    - Vektoroperationen (Kreuzprodukt, Skalarprodukt, Norm)
@author Reisen macht Spass... mit Pia und Dirk e.Kfm.
@date 2026-03-05
"""

import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from linear_algebra import (
    Matrix, Vector,
    gram_schmidt
)


class TestVector:
    """Tests für Vektoroperationen."""

    def test_dot_product(self):
        """Skalarprodukt: [1,2,3] · [4,5,6] = 4+10+18 = 32."""
        v1 = Vector([1, 2, 3])
        v2 = Vector([4, 5, 6])
        assert v1.dot(v2) == 32

    def test_norm(self):
        """Euklidische Norm: ||[3,4]|| = 5."""
        v = Vector([3, 4])
        assert abs(v.norm() - 5.0) < 1e-10

    def test_normalize(self):
        """Normierter Vektor hat Länge 1."""
        v = Vector([1, 2, 3])
        n = v.normalize()
        assert abs(n.norm() - 1.0) < 1e-10

    def test_cross_product(self):
        """Kreuzprodukt: [1,0,0] × [0,1,0] = [0,0,1]."""
        v1 = Vector([1, 0, 0])
        v2 = Vector([0, 1, 0])
        result = v1.cross(v2)
        assert result.components == [0, 0, 1]

    def test_cross_product_anticommutative(self):
        """Kreuzprodukt ist antikommutativ: v×w = -(w×v)."""
        v1 = Vector([1, 2, 3])
        v2 = Vector([4, 5, 6])
        c1 = v1.cross(v2)
        c2 = v2.cross(v1)
        for a, b in zip(c1.components, c2.components):
            assert abs(a + b) < 1e-10

    def test_addition(self):
        """Vektoraddition komponentenweise."""
        v1 = Vector([1, 2, 3])
        v2 = Vector([4, 5, 6])
        result = v1 + v2
        assert result.components == [5, 7, 9]

    def test_scalar_multiply(self):
        """Skalarmultiplikation."""
        v = Vector([1, 2, 3])
        result = v * 3
        assert result.components == [3, 6, 9]


class TestMatrix:
    """Tests für Matrixoperationen."""

    def test_creation(self):
        """Matrix aus Zeilenlisten erstellen."""
        m = Matrix([[1, 2], [3, 4]])
        assert m.rows == 2
        assert m.cols == 2
        assert m.get(0, 0) == 1
        assert m.get(1, 1) == 4

    def test_multiplication_identity(self):
        """Multiplikation mit Einheitsmatrix ergibt original."""
        m = Matrix([[1, 2], [3, 4]])
        eye = Matrix.identity(2)
        result = m @ eye
        assert result.get(0, 0) == 1
        assert result.get(0, 1) == 2
        assert result.get(1, 0) == 3
        assert result.get(1, 1) == 4

    def test_multiplication_basic(self):
        """Matrixmultiplikation: [[1,2],[3,4]] @ [[5,6],[7,8]] = [[19,22],[43,50]]."""
        a = Matrix([[1, 2], [3, 4]])
        b = Matrix([[5, 6], [7, 8]])
        result = a @ b
        assert result.get(0, 0) == 19
        assert result.get(0, 1) == 22
        assert result.get(1, 0) == 43
        assert result.get(1, 1) == 50

    def test_transpose(self):
        """Transponierung: [[1,2,3],[4,5,6]]^T = [[1,4],[2,5],[3,6]]."""
        m = Matrix([[1, 2, 3], [4, 5, 6]])
        t = m.transpose()
        assert t.rows == 3
        assert t.cols == 2
        assert t.get(0, 0) == 1
        assert t.get(1, 0) == 2
        assert t.get(2, 0) == 3
        assert t.get(0, 1) == 4

    def test_trace(self):
        """Spur (Summe der Diagonalelemente): Spur([[1,2],[3,4]]) = 5."""
        m = Matrix([[1, 2], [3, 4]])
        assert m.trace() == 5

    def test_determinant_2x2(self):
        """Determinante 2x2: det([[a,b],[c,d]]) = ad - bc."""
        m = Matrix([[3, 8], [4, 6]])
        # 3*6 - 8*4 = 18 - 32 = -14
        assert abs(m.determinant() - (-14)) < 1e-10

    def test_determinant_3x3(self):
        """Determinante einer 3x3-Matrix (bekannter Wert)."""
        m = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        # det = 1*(50-48) - 2*(40-42) + 3*(32-35) = 2 + 4 - 9 = -3
        assert abs(m.determinant() - (-3)) < 1e-10

    def test_determinant_singular(self):
        """Singuläre Matrix hat Determinante 0."""
        m = Matrix([[1, 2], [2, 4]])
        assert abs(m.determinant()) < 1e-10

    def test_inverse(self):
        """Inverse: A @ A^(-1) = I."""
        m = Matrix([[1, 2], [3, 4]])
        inv = m.inverse()
        result = m @ inv
        # Sollte Einheitsmatrix ergeben
        assert abs(result.get(0, 0) - 1) < 1e-10
        assert abs(result.get(0, 1)) < 1e-10
        assert abs(result.get(1, 0)) < 1e-10
        assert abs(result.get(1, 1) - 1) < 1e-10

    def test_inverse_singular_raises(self):
        """Singuläre Matrix hat keine Inverse."""
        m = Matrix([[1, 2], [2, 4]])
        with pytest.raises(ValueError):
            m.inverse()

    def test_solve_linear_system(self):
        """Lineares Gleichungssystem Ax = b lösen."""
        # 2x + y = 5
        # x + 3y = 10
        A = Matrix([[2, 1], [1, 3]])
        b = Vector([5, 10])
        x = A.solve(b)
        # Lösung: x=1, y=3
        assert abs(x.components[0] - 1.0) < 1e-10
        assert abs(x.components[1] - 3.0) < 1e-10

    def test_eigenvalues_symmetric(self):
        """Eigenwerte einer symmetrischen 2x2-Matrix."""
        # [[2, 1], [1, 2]] -> Eigenwerte: 1 und 3
        m = Matrix([[2, 1], [1, 2]])
        eigenvalues = m.eigenvalues()
        eigenvalues_sorted = sorted([abs(v) for v in eigenvalues])
        assert abs(eigenvalues_sorted[0] - 1.0) < 1e-6
        assert abs(eigenvalues_sorted[1] - 3.0) < 1e-6


class TestGramSchmidt:
    """Tests für das Gram-Schmidt-Orthogonalisierungsverfahren."""

    def test_orthogonality(self):
        """Ergebnis des Gram-Schmidt-Verfahrens ist orthogonal."""
        vectors = [
            Vector([1, 1, 0]),
            Vector([1, 0, 1]),
            Vector([0, 1, 1]),
        ]
        orthogonal = gram_schmidt(vectors)
        # Alle Paare sollten senkrecht stehen (Skalarprodukt ≈ 0)
        for i in range(len(orthogonal)):
            for j in range(i + 1, len(orthogonal)):
                dot = orthogonal[i].dot(orthogonal[j])
                assert abs(dot) < 1e-10, f"Vektor {i} und {j} nicht orthogonal: {dot}"

    def test_orthonormality(self):
        """Gram-Schmidt (normiert) liefert ONB: alle Vektoren haben Länge 1."""
        vectors = [Vector([1, 0, 0]), Vector([1, 1, 0])]
        orthonormal = gram_schmidt(vectors, normalize=True)
        for v in orthonormal:
            assert abs(v.norm() - 1.0) < 1e-10


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
