"""
@file test_vectors.py
@brief Umfassende Tests für das Modul vectors.py (Vector-Klasse + gram_schmidt).
@description
    Testet alle Methoden der Vector-Klasse:
    - Konstruktor, dim, components
    - dot(): Skalarprodukt
    - norm(): Euklidische Norm
    - normalize(): Normierung zum Einheitsvektor
    - cross(): Kreuzprodukt (3D)
    - Arithmetik: __add__, __sub__, __mul__, __rmul__
    - __len__, __getitem__, __repr__
    - gram_schmidt(): Orthogonalisierung / Orthonormalisierung

    Edge-Cases:
    - Nullvektor-Normierung → ValueError
    - Kreuzprodukt für nicht-3D-Vektoren → ValueError
    - Verschiedene Dimensionen bei add/sub/dot → ValueError
    - Linear abhängige Vektoren in gram_schmidt (werden übersprungen)

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from vectors import Vector, gram_schmidt


# =============================================================================
# TESTS: Vector-Konstruktor und Eigenschaften
# =============================================================================

class TestVectorConstructor:
    """Tests für Konstruktor und grundlegende Eigenschaften."""

    def test_dim_1d(self):
        """1D-Vektor hat Dimension 1."""
        v = Vector([5.0])
        assert v.dim == 1

    def test_dim_3d(self):
        """3D-Vektor hat Dimension 3."""
        v = Vector([1, 2, 3])
        assert v.dim == 3

    def test_components_stored_correctly(self):
        """Komponenten werden korrekt gespeichert."""
        v = Vector([1.0, 2.0, 3.0])
        assert v.components == [1.0, 2.0, 3.0]

    def test_components_are_copy(self):
        """Änderungen an der Eingabeliste ändern den Vektor nicht."""
        lst = [1.0, 2.0, 3.0]
        v = Vector(lst)
        lst[0] = 999.0
        assert v.components[0] == 1.0

    def test_len(self):
        """__len__ gibt die Dimension zurück."""
        v = Vector([1, 2, 3, 4])
        assert len(v) == 4

    def test_getitem(self):
        """__getitem__ gibt die korrekte Komponente zurück."""
        v = Vector([10, 20, 30])
        assert v[0] == 10
        assert v[1] == 20
        assert v[2] == 30

    def test_repr(self):
        """__repr__ gibt lesbaren String zurück."""
        v = Vector([1, 2])
        r = repr(v)
        assert 'Vector' in r


# =============================================================================
# TESTS: Vector.dot()
# =============================================================================

class TestVectorDot:
    """Tests für das Skalarprodukt."""

    def test_standard_dot(self):
        """Skalarprodukt bekannter Vektoren: [1,2,3]·[4,5,6] = 32."""
        v = Vector([1, 2, 3])
        w = Vector([4, 5, 6])
        assert v.dot(w) == 32

    def test_dot_commutative(self):
        """Skalarprodukt ist kommutativ: v·w = w·v."""
        v = Vector([1, 2, 3])
        w = Vector([4, 5, 6])
        assert v.dot(w) == w.dot(v)

    def test_dot_with_zero_vector(self):
        """Skalarprodukt mit Nullvektor ist 0."""
        v = Vector([1, 2, 3])
        z = Vector([0, 0, 0])
        assert v.dot(z) == 0

    def test_dot_orthogonal(self):
        """Senkrechte Vektoren haben Skalarprodukt 0."""
        v = Vector([1, 0, 0])
        w = Vector([0, 1, 0])
        assert v.dot(w) == 0

    def test_dot_self_is_norm_squared(self):
        """v·v = ||v||²."""
        v = Vector([3, 4])
        assert abs(v.dot(v) - v.norm()**2) < 1e-12

    def test_dot_dimension_mismatch_raises(self):
        """Unterschiedliche Dimensionen → ValueError."""
        v = Vector([1, 2])
        w = Vector([1, 2, 3])
        with pytest.raises(ValueError):
            v.dot(w)

    def test_dot_1d(self):
        """1D-Skalarprodukt: [3]·[4] = 12."""
        v = Vector([3])
        w = Vector([4])
        assert v.dot(w) == 12

    def test_dot_negative_components(self):
        """Negative Komponenten korrekt verarbeitet."""
        v = Vector([-1, -2, -3])
        w = Vector([1, 2, 3])
        assert v.dot(w) == -14


# =============================================================================
# TESTS: Vector.norm()
# =============================================================================

class TestVectorNorm:
    """Tests für die euklidische Norm."""

    def test_norm_pythagoras(self):
        """Pythagoreisches Tripel: ||[3,4]|| = 5."""
        v = Vector([3, 4])
        assert abs(v.norm() - 5.0) < 1e-12

    def test_norm_unit_vector(self):
        """Norm eines Einheitsvektors ist 1."""
        v = Vector([1, 0, 0])
        assert abs(v.norm() - 1.0) < 1e-12

    def test_norm_zero_vector(self):
        """Norm des Nullvektors ist 0."""
        v = Vector([0, 0, 0])
        assert v.norm() == 0.0

    def test_norm_3d(self):
        """||[1,1,1]|| = √3."""
        v = Vector([1, 1, 1])
        assert abs(v.norm() - math.sqrt(3)) < 1e-12

    def test_norm_negative_components(self):
        """Norm ist immer nicht-negativ."""
        v = Vector([-3, -4])
        assert abs(v.norm() - 5.0) < 1e-12


# =============================================================================
# TESTS: Vector.normalize()
# =============================================================================

class TestVectorNormalize:
    """Tests für die Normierung."""

    def test_normalized_has_unit_norm(self):
        """Normierter Vektor hat Norm 1."""
        v = Vector([3, 4])
        n = v.normalize()
        assert abs(n.norm() - 1.0) < 1e-12

    def test_normalized_direction(self):
        """Normierter Vektor zeigt in die gleiche Richtung."""
        v = Vector([3, 4])
        n = v.normalize()
        # n = (3/5, 4/5)
        assert abs(n.components[0] - 0.6) < 1e-12
        assert abs(n.components[1] - 0.8) < 1e-12

    def test_normalize_zero_vector_raises(self):
        """Normierung des Nullvektors → ValueError."""
        v = Vector([0, 0, 0])
        with pytest.raises(ValueError):
            v.normalize()

    def test_normalize_preserves_dim(self):
        """Normierung bewahrt die Dimension."""
        v = Vector([1, 2, 3, 4])
        n = v.normalize()
        assert n.dim == 4

    def test_normalize_already_unit(self):
        """Bereits normierter Vektor bleibt normiert."""
        v = Vector([1, 0, 0])
        n = v.normalize()
        assert abs(n.norm() - 1.0) < 1e-12


# =============================================================================
# TESTS: Vector.cross()
# =============================================================================

class TestVectorCross:
    """Tests für das Kreuzprodukt."""

    def test_cross_standard(self):
        """[1,0,0] × [0,1,0] = [0,0,1]."""
        v = Vector([1, 0, 0])
        w = Vector([0, 1, 0])
        result = v.cross(w)
        assert abs(result.components[0] - 0) < 1e-12
        assert abs(result.components[1] - 0) < 1e-12
        assert abs(result.components[2] - 1) < 1e-12

    def test_cross_anticommutative(self):
        """Kreuzprodukt ist antikommutativ: v×w = -(w×v)."""
        v = Vector([1, 2, 3])
        w = Vector([4, 5, 6])
        vxw = v.cross(w)
        wxv = w.cross(v)
        for i in range(3):
            assert abs(vxw.components[i] + wxv.components[i]) < 1e-12

    def test_cross_with_itself_is_zero(self):
        """v × v = 0."""
        v = Vector([1, 2, 3])
        result = v.cross(v)
        for c in result.components:
            assert abs(c) < 1e-12

    def test_cross_perpendicular_to_both(self):
        """Kreuzprodukt ist senkrecht zu beiden Eingangsvektoren."""
        v = Vector([1, 2, 3])
        w = Vector([4, 5, 6])
        c = v.cross(w)
        assert abs(c.dot(v)) < 1e-12
        assert abs(c.dot(w)) < 1e-12

    def test_cross_wrong_dimension_raises(self):
        """Kreuzprodukt für nicht-3D → ValueError."""
        v = Vector([1, 2])
        w = Vector([3, 4])
        with pytest.raises(ValueError):
            v.cross(w)

    def test_cross_result_is_3d(self):
        """Ergebnis des Kreuzprodukts ist 3-dimensional."""
        v = Vector([1, 0, 0])
        w = Vector([0, 0, 1])
        result = v.cross(w)
        assert result.dim == 3


# =============================================================================
# TESTS: Arithmetische Operationen
# =============================================================================

class TestVectorArithmetic:
    """Tests für +, -, * Operationen."""

    def test_addition(self):
        """Vektoraddition komponentenweise."""
        v = Vector([1, 2, 3])
        w = Vector([4, 5, 6])
        result = v + w
        assert result.components == [5, 7, 9]

    def test_subtraction(self):
        """Vektorsubtraktion komponentenweise."""
        v = Vector([4, 5, 6])
        w = Vector([1, 2, 3])
        result = v - w
        assert result.components == [3, 3, 3]

    def test_scalar_multiplication(self):
        """Skalarmultiplikation: v * 2."""
        v = Vector([1, 2, 3])
        result = v * 2
        assert result.components == [2, 4, 6]

    def test_right_scalar_multiplication(self):
        """Rechtsseitige Multiplikation: 2 * v."""
        v = Vector([1, 2, 3])
        result = 2 * v
        assert result.components == [2, 4, 6]

    def test_add_dimension_mismatch_raises(self):
        """Addition bei verschiedenen Dimensionen → ValueError."""
        v = Vector([1, 2])
        w = Vector([1, 2, 3])
        with pytest.raises(ValueError):
            v + w

    def test_sub_dimension_mismatch_raises(self):
        """Subtraktion bei verschiedenen Dimensionen → ValueError."""
        v = Vector([1, 2])
        w = Vector([1, 2, 3])
        with pytest.raises(ValueError):
            v - w

    def test_add_preserves_dim(self):
        """Addition bewahrt die Dimension."""
        v = Vector([1, 2, 3])
        w = Vector([0, 0, 0])
        result = v + w
        assert result.dim == 3

    def test_scalar_zero(self):
        """Multiplikation mit 0 ergibt Nullvektor."""
        v = Vector([1, 2, 3])
        result = v * 0
        assert all(c == 0 for c in result.components)


# =============================================================================
# TESTS: gram_schmidt()
# =============================================================================

class TestGramSchmidt:
    """Tests für die Gram-Schmidt-Orthogonalisierung."""

    def test_already_orthogonal(self):
        """Bereits orthogonale Vektoren bleiben orthogonal."""
        v1 = Vector([1, 0, 0])
        v2 = Vector([0, 1, 0])
        v3 = Vector([0, 0, 1])
        result = gram_schmidt([v1, v2, v3])
        assert len(result) == 3

    def test_orthogonality_after(self):
        """Ergebnis-Vektoren sind paarweise orthogonal."""
        v1 = Vector([1, 1, 0])
        v2 = Vector([0, 1, 1])
        v3 = Vector([1, 0, 1])
        result = gram_schmidt([v1, v2, v3])
        for i in range(len(result)):
            for j in range(i + 1, len(result)):
                dot = result[i].dot(result[j])
                assert abs(dot) < 1e-10, f"u_{i} und u_{j} sind nicht orthogonal: {dot}"

    def test_normalize_true(self):
        """Mit normalize=True: Alle Vektoren haben Norm 1 (Orthonormalbasis)."""
        v1 = Vector([1, 1, 0])
        v2 = Vector([0, 1, 1])
        result = gram_schmidt([v1, v2], normalize=True)
        for u in result:
            assert abs(u.norm() - 1.0) < 1e-10

    def test_normalize_false(self):
        """Mit normalize=False: Vektoren haben nicht zwingend Norm 1."""
        v1 = Vector([3, 4, 0])
        result = gram_schmidt([v1], normalize=False)
        assert len(result) == 1

    def test_linear_dependent_skipped(self):
        """Linear abhängige Vektoren werden übersprungen."""
        v1 = Vector([1, 0, 0])
        v2 = Vector([2, 0, 0])  # linear abhängig von v1
        result = gram_schmidt([v1, v2])
        # Nur 1 Vektor darf im Ergebnis sein
        assert len(result) == 1

    def test_single_vector(self):
        """Ein einzelner Vektor bleibt erhalten."""
        v = Vector([3, 4, 0])
        result = gram_schmidt([v])
        assert len(result) == 1
        assert abs(result[0].norm() - v.norm()) < 1e-10

    def test_empty_list(self):
        """Leere Eingabe → leere Ausgabe."""
        result = gram_schmidt([])
        assert result == []

    def test_2d_basis(self):
        """Orthogonalisierung einer 2D-Basis."""
        v1 = Vector([2, 1])
        v2 = Vector([1, 3])
        result = gram_schmidt([v1, v2])
        assert len(result) == 2
        # Orthogonalitätsprüfung
        assert abs(result[0].dot(result[1])) < 1e-10
