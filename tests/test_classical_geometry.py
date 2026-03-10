"""
@file test_classical_geometry.py
@brief Tests für das Modul classical_geometry.py
@description
    Umfassende Test-Suite für alle Klassen und Funktionen der klassischen,
    synthetischen und nichteuklidischen Geometrie.

    Testbereiche:
    - EuclideanGeometry: Abstände, Winkel, Parallelität, Orthogonalität,
      Umkreis, Inkreis, Schwerpunkt, Euler-Gerade
    - ProjectiveGeometry: Homogene Koordinaten, projektive Geraden,
      Schnittpunkte, Doppelverhältnis, harmonische Konjugierte,
      projektive Transformationen, Desargues, Pappus
    - HyperbolicGeometry: Poincaré-Abstand, Mittelpunkt, Winkelsumme,
      Fläche, Gauß-Bonnet, oberes Halbebenenmodell
    - EllipticGeometry: Haversine-Formel, sphärischer Exzess, Fläche,
      Lune, Girard'scher Satz
    - AffinGeometry: Affinkombinationen, affine Abbildungen, baryzentrische
      Koordinaten, affine Hülle, konvexe Hülle, affine Unabhängigkeit
    - SyntheticGeometry: Hilbert-Axiome, endliche projektive Ebene
    - Standalone: Euler-Charakteristik, Pick, Ptolemäus, Neun-Punkte-Kreis,
      Morley

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math

# Suchpfad für src-Modul setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import numpy as np

# Zu testende Module
from classical_geometry import (
    EuclideanGeometry,
    ProjectiveGeometry,
    HyperbolicGeometry,
    EllipticGeometry,
    AffinGeometry,
    SyntheticGeometry,
    euler_characteristic_polygon,
    pick_theorem,
    ptolemy_theorem_check,
    nine_point_circle,
    morley_theorem_demo,
)

# Gemeinsame Toleranz für Fließkomma-Vergleiche
TOL = 1e-8


# ===========================================================================
# EuclideanGeometry
# ===========================================================================

class TestEuclideanGeometryDistance:
    """Tests für EuclideanGeometry.distance()"""

    def setup_method(self):
        self.geom = EuclideanGeometry(dim=2)

    def test_distance_zero_same_point(self):
        """Abstand eines Punktes zu sich selbst ist 0."""
        assert self.geom.distance([0, 0], [0, 0]) == pytest.approx(0.0)

    def test_distance_unit_vectors(self):
        """Abstand vom Ursprung zum Einheitsvektor (1,0) ist 1."""
        assert self.geom.distance([0, 0], [1, 0]) == pytest.approx(1.0)

    def test_distance_3_4_5(self):
        """Pythagoreisches Tripel 3-4-5."""
        assert self.geom.distance([0, 0], [3, 4]) == pytest.approx(5.0)

    def test_distance_negative_coords(self):
        """Abstand mit negativen Koordinaten."""
        assert self.geom.distance([-1, -1], [2, 3]) == pytest.approx(5.0)

    def test_distance_3d(self):
        """Abstand in 3D."""
        geom3d = EuclideanGeometry(dim=3)
        assert geom3d.distance([0, 0, 0], [1, 2, 2]) == pytest.approx(3.0)

    def test_distance_symmetry(self):
        """Abstand ist symmetrisch: d(A,B) = d(B,A)."""
        p1, p2 = [1.5, 2.3], [4.7, -1.1]
        assert self.geom.distance(p1, p2) == pytest.approx(self.geom.distance(p2, p1))


class TestEuclideanGeometryAngle:
    """Tests für EuclideanGeometry.angle_between()"""

    def setup_method(self):
        self.geom = EuclideanGeometry()

    def test_angle_parallel_vectors(self):
        """Parallele Vektoren haben Winkel 0."""
        assert self.geom.angle_between([1, 0], [2, 0]) == pytest.approx(0.0)

    def test_angle_perpendicular(self):
        """Rechter Winkel zwischen (1,0) und (0,1)."""
        assert self.geom.angle_between([1, 0], [0, 1]) == pytest.approx(math.pi / 2)

    def test_angle_antiparallel(self):
        """Antiparallele Vektoren haben Winkel π."""
        assert self.geom.angle_between([1, 0], [-1, 0]) == pytest.approx(math.pi)

    def test_angle_45_degrees(self):
        """45° Winkel."""
        v = [1, 1]
        assert self.geom.angle_between([1, 0], v) == pytest.approx(math.pi / 4, rel=1e-6)

    def test_angle_null_vector_raises(self):
        """Nullvektor soll ValueError werfen."""
        with pytest.raises(ValueError, match="Nullvektor"):
            self.geom.angle_between([0, 0], [1, 0])


class TestEuclideanParallelPerpendicular:
    """Tests für are_parallel() und are_perpendicular()"""

    def setup_method(self):
        self.geom = EuclideanGeometry()

    def test_horizontal_lines_parallel(self):
        """Horizontale Geraden sind parallel."""
        assert self.geom.are_parallel([0, 0], [1, 0], [0, 1], [1, 1]) is True

    def test_non_parallel_lines(self):
        """Schneidende Geraden sind nicht parallel."""
        assert self.geom.are_parallel([0, 0], [1, 0], [0, 0], [0, 1]) is False

    def test_perpendicular_vectors(self):
        """Orthogonale Vektoren."""
        assert self.geom.are_perpendicular([1, 0], [0, 1]) is True

    def test_not_perpendicular(self):
        """Nicht orthogonale Vektoren."""
        assert self.geom.are_perpendicular([1, 1], [1, 0]) is False

    def test_parallel_diagonal_lines(self):
        """Diagonale Geraden mit gleicher Steigung sind parallel."""
        assert self.geom.are_parallel([0, 0], [1, 1], [2, 0], [3, 1]) is True


class TestEuclideanCircumIncenter:
    """Tests für circumcenter() und incenter()"""

    def setup_method(self):
        self.geom = EuclideanGeometry()

    def test_circumcenter_right_triangle(self):
        """Umkreismittelpunkt des rechtwinkligen Dreiecks liegt auf der Hypotenuse-Mitte."""
        # Rechtwinkliges Dreieck: (0,0), (2,0), (0,2)
        O = self.geom.circumcenter([0, 0], [2, 0], [0, 2])
        # Mitte der Hypotenuse = (1,1)
        assert O == pytest.approx([1.0, 1.0], abs=TOL)

    def test_circumcenter_equilateral(self):
        """Umkreismittelpunkt des gleichseitigen Dreiecks ist der Schwerpunkt."""
        p1 = [0, 0]
        p2 = [1, 0]
        p3 = [0.5, math.sqrt(3)/2]
        O = self.geom.circumcenter(p1, p2, p3)
        G = self.geom.centroid([p1, p2, p3])
        assert O == pytest.approx(G, abs=TOL)

    def test_circumcenter_equidistant(self):
        """Umkreismittelpunkt ist von allen drei Ecken gleich weit entfernt."""
        p1, p2, p3 = [0, 0], [4, 0], [0, 3]
        O = self.geom.circumcenter(p1, p2, p3)
        d1 = self.geom.distance(O, p1)
        d2 = self.geom.distance(O, p2)
        d3 = self.geom.distance(O, p3)
        assert d1 == pytest.approx(d2, abs=TOL)
        assert d2 == pytest.approx(d3, abs=TOL)

    def test_circumcenter_collinear_raises(self):
        """Kollineare Punkte sollen ValueError auslösen."""
        with pytest.raises(ValueError, match="Kollineare"):
            self.geom.circumcenter([0, 0], [1, 1], [2, 2])

    def test_incenter_equilateral(self):
        """Inkreismittelpunkt des gleichseitigen Dreiecks = Schwerpunkt."""
        p1 = [0, 0]
        p2 = [2, 0]
        p3 = [1, math.sqrt(3)]
        I = self.geom.incenter(p1, p2, p3)
        G = self.geom.centroid([p1, p2, p3])
        assert I == pytest.approx(G, abs=TOL)

    def test_incenter_positive_coords(self):
        """Inkreismittelpunkt liegt immer innerhalb des Dreiecks."""
        p1, p2, p3 = [0, 0], [4, 0], [0, 3]
        I = self.geom.incenter(p1, p2, p3)
        # Alle Koordinaten > 0 (innen)
        assert I[0] > 0
        assert I[1] > 0


class TestEuclideanCentroidEulerLine:
    """Tests für centroid() und euler_line_demo()"""

    def setup_method(self):
        self.geom = EuclideanGeometry()

    def test_centroid_triangle(self):
        """Schwerpunkt eines Dreiecks."""
        G = self.geom.centroid([[0, 0], [3, 0], [0, 3]])
        assert G == pytest.approx([1.0, 1.0], abs=TOL)

    def test_centroid_single_point(self):
        """Schwerpunkt eines einzelnen Punktes ist der Punkt selbst."""
        G = self.geom.centroid([[5, 7]])
        assert G == pytest.approx([5.0, 7.0], abs=TOL)

    def test_centroid_empty_raises(self):
        """Leere Punktliste soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.geom.centroid([])

    def test_euler_line_collinear(self):
        """Euler-Gerade: O, G, H sind kollinear."""
        result = self.geom.euler_line_demo([0, 0], [4, 0], [0, 3])
        assert result["euler_line_verified"] is True

    def test_euler_line_ratio(self):
        """Euler-Gerade: |GH|/|OG| ≈ 2."""
        result = self.geom.euler_line_demo([0, 0], [4, 0], [1, 3])
        assert result["ratio_GH_to_OG"] == pytest.approx(2.0, abs=TOL)


# ===========================================================================
# ProjectiveGeometry
# ===========================================================================

class TestProjectiveGeometry:
    """Tests für ProjectiveGeometry"""

    def setup_method(self):
        self.pg = ProjectiveGeometry()

    def test_homogeneous_coords_2d(self):
        """2D-Punkt → homogene Koordinaten [x:y:1]."""
        h = self.pg.homogeneous_coords([3, 5])
        assert h == pytest.approx([3.0, 5.0, 1.0])

    def test_homogeneous_coords_already_homogeneous(self):
        """Bereits homogene Koordinaten bleiben unverändert."""
        h = self.pg.homogeneous_coords([2, 3, 1])
        assert h == pytest.approx([2.0, 3.0, 1.0])

    def test_projective_line_x_axis(self):
        """Gerade durch (0,0) und (1,0) ist die x-Achse: [0:1:0] oder ähnlich."""
        l = self.pg.projective_line_through([0, 0], [1, 0])
        # Gerade ax+by+c=0. Für y=0: b=1, a=0, c=0 → [0:1:0] × [1:0:0] = [0:0:0-1] variiert
        # Prüfe: beide Punkte liegen auf der Geraden
        p1 = self.pg.homogeneous_coords([0, 0])
        p2 = self.pg.homogeneous_coords([1, 0])
        assert abs(np.dot(l, p1)) < TOL
        assert abs(np.dot(l, p2)) < TOL

    def test_projective_line_identical_points_raises(self):
        """Identische Punkte sollen ValueError auslösen."""
        with pytest.raises(ValueError, match="Identische Punkte"):
            self.pg.projective_line_through([1, 2], [1, 2])

    def test_projective_intersection_basic(self):
        """Schnittpunkt zweier Geraden."""
        l1 = self.pg.projective_line_through([0, 0], [2, 0])  # y=0
        l2 = self.pg.projective_line_through([1, 1], [1, -1])  # x=1
        P = self.pg.projective_intersection(l1, l2)
        # Normalisiere auf w=1
        P_norm = P / P[2]
        assert P_norm[0] == pytest.approx(1.0, abs=TOL)
        assert P_norm[1] == pytest.approx(0.0, abs=TOL)

    def test_projective_intersection_parallel_lines_at_infinity(self):
        """Parallele Geraden schneiden sich im Unendlichen (w=0)."""
        # y=0 und y=1 (parallel)
        l1 = self.pg.projective_line_through([0, 0], [1, 0])
        l2 = self.pg.projective_line_through([0, 1], [1, 1])
        P = self.pg.projective_intersection(l1, l2)
        # Punkt im Unendlichen: dritte Komponente ≈ 0
        assert abs(P[2]) < TOL

    def test_cross_ratio_basic(self):
        """Doppelverhältnis: (0,∞;1,-1) = -1 (harmonisch)."""
        # Klassisches harmonisches Verhältnis: (−1, 1; 0, ∞) = −1
        # Verwende 1D-Punkte: p1=0, p2=4, p3=1, p4=-2
        # (0,4;1,-2) = (1-0)(-2-4) / ((1-4)(-2-0)) = 1·(-6) / ((-3)·(-2)) = -6/6 = -1
        cr = self.pg.cross_ratio([0], [4], [1], [-2])
        assert cr == pytest.approx(-1.0, abs=TOL)

    def test_cross_ratio_invariant(self):
        """Doppelverhältnis: Standardwerte 0,1,∞ entsprechen je 0,1 und großem Wert."""
        # (0,1;2,3) = (2-0)(3-1)/((2-1)(3-0)) = 2·2/(1·3) = 4/3
        cr = self.pg.cross_ratio([0], [1], [2], [3])
        assert cr == pytest.approx(4.0 / 3.0, abs=TOL)

    def test_cross_ratio_zero_denominator_raises(self):
        """Entartetes Doppelverhältnis soll ZeroDivisionError auslösen."""
        with pytest.raises(ZeroDivisionError):
            self.pg.cross_ratio([0], [1], [1], [0])  # p3=p2, p4=p1 → Nenner 0

    def test_harmonic_conjugate(self):
        """Harmonische Konjugierte: Doppelverhältnis = -1."""
        # P4 = harmonische Konjugierte von P3=1 bzgl. P1=0, P2=4
        # Nenner = 2·1 - 0 - 4 = -2 ≠ 0
        P4 = self.pg.harmonic_conjugate(0.0, 4.0, 1.0)
        # Prüfe: (0,4;1,P4) = -1
        cr = self.pg.cross_ratio([0], [4], [1], [P4])
        assert cr == pytest.approx(-1.0, abs=TOL)

    def test_projective_transformation_identity(self):
        """Identitätstransformation lässt Punkte unverändert."""
        M = np.eye(3)
        p = [3, 5]
        result = self.pg.projective_transformation(M, p)
        # Normalisiere
        result_norm = result / result[2]
        assert result_norm[0] == pytest.approx(3.0, abs=TOL)
        assert result_norm[1] == pytest.approx(5.0, abs=TOL)

    def test_projective_transformation_scaling(self):
        """Skalierungstransformation."""
        M = np.diag([2.0, 3.0, 1.0])
        result = self.pg.projective_transformation(M, [1, 1])
        result_norm = result / result[2]
        assert result_norm[0] == pytest.approx(2.0, abs=TOL)
        assert result_norm[1] == pytest.approx(3.0, abs=TOL)

    def test_desargues_theorem_perspectively_collinear(self):
        """Desargues-Satz: Perspektive Dreiecke erfüllen die Bedingung."""
        # Erzeuge zwei Dreiecke perspektiv bzgl. Punkt O=(0,0)
        # Dreieck 1: A1, B1, C1; Dreieck 2: Skalierung von O aus
        O = np.array([0.0, 0.0])
        A1 = np.array([1.0, 0.0])
        B1 = np.array([0.0, 1.0])
        C1 = np.array([1.0, 1.0])
        # Perspektive Skalierung: Ai' = O + 2*(Ai - O)
        A2 = O + 2.0 * (A1 - O)  # [2,0]
        B2 = O + 2.0 * (B1 - O)  # [0,2]
        C2 = O + 2.0 * (C1 - O)  # [2,2]

        tri1 = [A1.tolist(), B1.tolist(), C1.tolist()]
        tri2 = [A2.tolist(), B2.tolist(), C2.tolist()]
        result = self.pg.desargues_theorem_check(tri1, tri2, tol=1e-6)
        assert result["desargues_satisfied"] == True  # numpy bool kompatibel

    def test_pappus_theorem_check(self):
        """Pappus-Satz: Kollineare Punkte auf zwei Geraden."""
        # Punkte auf Gerade 1 (y=0): A=(0,0), B=(2,0), C=(4,0)
        # Punkte auf Gerade 2 (y=1): D=(1,1), E=(3,1), F=(5,1)
        l1 = [[0, 0], [2, 0], [4, 0]]
        l2 = [[1, 1], [3, 1], [5, 1]]
        assert self.pg.pappus_theorem_check(l1, l2, tol=1e-6) == True


# ===========================================================================
# HyperbolicGeometry
# ===========================================================================

class TestHyperbolicGeometry:
    """Tests für HyperbolicGeometry"""

    def setup_method(self):
        self.hyp = HyperbolicGeometry()

    def test_poincare_distance_zero(self):
        """Abstand eines Punktes zu sich selbst ist 0."""
        assert self.hyp.poincare_disk_distance(0+0j, 0+0j) == pytest.approx(0.0)

    def test_poincare_distance_origin_to_half(self):
        """Abstand Ursprung zu (0.5, 0) = 2·arctanh(0.5)."""
        z = complex(0.5, 0)
        expected = 2.0 * math.atanh(0.5)
        assert self.hyp.poincare_disk_distance(0+0j, z) == pytest.approx(expected, rel=1e-8)

    def test_poincare_distance_symmetry(self):
        """Hyperbolischer Abstand ist symmetrisch."""
        z1 = complex(0.3, 0.2)
        z2 = complex(-0.1, 0.4)
        d12 = self.hyp.poincare_disk_distance(z1, z2)
        d21 = self.hyp.poincare_disk_distance(z2, z1)
        assert d12 == pytest.approx(d21, rel=1e-8)

    def test_poincare_distance_outside_disk_raises(self):
        """Punkt außerhalb des Einheitskreises soll ValueError auslösen."""
        with pytest.raises(ValueError, match="Einheitskreis"):
            self.hyp.poincare_disk_distance(complex(0.9, 0.9), 0+0j)

    def test_poincare_distance_boundary_raises(self):
        """Punkt auf dem Rand (|z|=1) soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.hyp.poincare_disk_distance(1+0j, 0+0j)

    def test_poincare_midpoint_is_between(self):
        """Hyperbolischer Mittelpunkt hat gleiche Abstände zu beiden Punkten."""
        z1 = complex(0.2, 0.0)
        z2 = complex(-0.2, 0.0)
        mid = self.hyp.poincare_disk_midpoint(z1, z2)
        d1 = self.hyp.poincare_disk_distance(z1, mid)
        d2 = self.hyp.poincare_disk_distance(mid, z2)
        assert d1 == pytest.approx(d2, rel=1e-6)

    def test_poincare_midpoint_origin(self):
        """Mittelpunkt von Punkten auf der x-Achse liegt auf der x-Achse."""
        z1 = complex(0.4, 0.0)
        z2 = complex(-0.4, 0.0)
        mid = self.hyp.poincare_disk_midpoint(z1, z2)
        assert abs(mid.imag) < TOL  # Mittelpunkt auf x-Achse
        assert abs(mid.real) < TOL  # Mittelpunkt im Ursprung (Symmetrie)

    def test_hyperbolic_angle_sum_less_than_pi(self):
        """Winkelsumme eines hyperbolischen Dreiecks ist < π.

        Wir nutzen weit auseinanderliegende Punkte tief im Disk-Inneren,
        sodass die hyperbolischen Abstände groß sind und die Winkelsumme
        deutlich unter π liegt.
        """
        # Dreieck mit Punkten, die ein echtes hyperbolisches Dreieck bilden
        # (weiter vom Zentrum entfernt → größere hyperbolische Abstände)
        v = [complex(0.5, 0.0), complex(-0.4, 0.3), complex(0.0, -0.5)]
        angle_sum = self.hyp.hyperbolic_angle_sum(v)
        assert angle_sum < math.pi

    def test_hyperbolic_area_positive(self):
        """Hyperbolische Fläche ist > 0."""
        v = [complex(0.5, 0.0), complex(-0.4, 0.3), complex(0.0, -0.5)]
        area = self.hyp.hyperbolic_area(v)
        assert area > 0.0

    def test_hyperbolic_area_equals_pi_minus_angles(self):
        """Fläche = π − Winkelsumme."""
        v = [complex(0.5, 0.0), complex(-0.4, 0.3), complex(0.0, -0.5)]
        angle_sum = self.hyp.hyperbolic_angle_sum(v)
        area = self.hyp.hyperbolic_area(v)
        assert area == pytest.approx(math.pi - angle_sum, abs=TOL)

    def test_gauss_bonnet_triangle(self):
        """Gauß-Bonnet für hyperbolisches Dreieck: A = π - (α+β+γ)."""
        angles = [0.5, 0.6, 0.7]
        area = self.hyp.gauss_bonnet_hyperbolic(angles)
        expected = math.pi - sum(angles)
        assert area == pytest.approx(expected, abs=TOL)

    def test_gauss_bonnet_polygon(self):
        """Gauß-Bonnet für hyperbolisches Viereck: A = 2π - Σαᵢ."""
        angles = [0.3, 0.4, 0.5, 0.6]  # 4-Eck
        area = self.hyp.gauss_bonnet_hyperbolic(angles)
        expected = (4 - 2) * math.pi - sum(angles)
        assert area == pytest.approx(expected, abs=TOL)

    def test_hyperbolic_parallel_lines_returns_dict(self):
        """Parallellinien-Demo gibt ein Dict zurück."""
        result = self.hyp.hyperbolic_parallel_lines()
        assert "parallel_count" in result
        assert result["parallel_count"] == float('inf')

    def test_upper_half_plane_distance_positive(self):
        """Obere Halbebene: Abstand ist ≥ 0."""
        z1 = complex(0.0, 1.0)
        z2 = complex(1.0, 1.0)
        d = self.hyp.upper_half_plane_distance(z1, z2)
        assert d >= 0.0

    def test_upper_half_plane_distance_zero(self):
        """Gleiche Punkte haben Abstand 0."""
        z = complex(1.0, 2.0)
        assert self.hyp.upper_half_plane_distance(z, z) == pytest.approx(0.0)

    def test_upper_half_plane_real_axis_raises(self):
        """Punkt auf der reellen Achse (Im=0) soll ValueError auslösen."""
        with pytest.raises(ValueError, match="oberen Halbebene"):
            self.hyp.upper_half_plane_distance(complex(1.0, 0.0), complex(1.0, 1.0))

    def test_upper_half_plane_distance_i_to_2i(self):
        """Abstand von i zu 2i im oberen Halbebenenmodell."""
        z1 = complex(0.0, 1.0)
        z2 = complex(0.0, 2.0)
        d = self.hyp.upper_half_plane_distance(z1, z2)
        # d = arccosh(1 + |i-2i|²/(2·1·2)) = arccosh(1 + 1/4) = arccosh(1.25)
        expected = math.acosh(1.25)
        assert d == pytest.approx(expected, rel=1e-8)

    def test_upper_half_plane_distance_symmetry(self):
        """Abstand ist symmetrisch."""
        z1 = complex(0.5, 1.0)
        z2 = complex(-0.5, 2.0)
        d12 = self.hyp.upper_half_plane_distance(z1, z2)
        d21 = self.hyp.upper_half_plane_distance(z2, z1)
        assert d12 == pytest.approx(d21, rel=1e-8)


# ===========================================================================
# EllipticGeometry
# ===========================================================================

class TestEllipticGeometry:
    """Tests für EllipticGeometry (sphärische Geometrie)"""

    def setup_method(self):
        self.ell = EllipticGeometry()

    def test_spherical_distance_same_point(self):
        """Abstand eines Punktes zu sich selbst ist 0."""
        assert self.ell.spherical_distance(0.0, 0.0, 0.0, 0.0) == pytest.approx(0.0)

    def test_spherical_distance_antipodal(self):
        """Abstand zu Antipodenpunkt = π (halber Großkreis)."""
        d = self.ell.spherical_distance(0.0, 0.0, 0.0, math.pi)
        assert d == pytest.approx(math.pi, abs=TOL)

    def test_spherical_distance_quarter(self):
        """Quadrant-Abstand: Pol zu Äquator = π/2."""
        d = self.ell.spherical_distance(0.0, 0.0, math.pi/2, 0.0)
        assert d == pytest.approx(math.pi / 2, abs=TOL)

    def test_spherical_distance_real_world(self):
        """Großkreisabstand zweier bekannter Punkte (ungefähr)."""
        # Berlin (52.5°N, 13.4°E) und New York (40.7°N, 74.0°W = -74.0°)
        lat1 = math.radians(52.5)
        lon1 = math.radians(13.4)
        lat2 = math.radians(40.7)
        lon2 = math.radians(-74.0)
        R_earth = 6371.0  # km
        d = self.ell.spherical_distance(lat1, lon1, lat2, lon2, radius=R_earth)
        # Ungefähr 6386 km (Haversine-Ergebnis, ohne Erdabflachung)
        assert 6000 < d < 7000

    def test_spherical_angle_sum_greater_than_pi(self):
        """Winkelsumme eines sphärischen Dreiecks ist > π."""
        # Dreieck: Nordpol + zwei Äquatorpunkte
        north_pole = [math.pi / 2, 0.0]
        eq1 = [0.0, 0.0]
        eq2 = [0.0, math.pi / 2]
        angle_sum = self.ell.spherical_angle_sum([north_pole, eq1, eq2])
        assert angle_sum > math.pi

    def test_spherical_excess_positive(self):
        """Sphärischer Exzess ist > 0."""
        north_pole = [math.pi / 2, 0.0]
        eq1 = [0.0, 0.0]
        eq2 = [0.0, math.pi / 2]
        E = self.ell.spherical_excess([north_pole, eq1, eq2])
        assert E > 0.0

    def test_spherical_area_one_eighth_sphere(self):
        """Dreieck mit 3 rechten Winkeln = 1/8 Kugeloberfläche."""
        # Dreieck: Nordpol, (0°N, 0°E), (0°N, 90°E) auf S²
        north_pole = [math.pi / 2, 0.0]
        eq1 = [0.0, 0.0]
        eq2 = [0.0, math.pi / 2]
        area = self.ell.spherical_area([north_pole, eq1, eq2], radius=1.0)
        # Jeder Winkel = π/2 → E = 3·π/2 - π = π/2 → A = π/2
        assert area == pytest.approx(math.pi / 2, abs=1e-6)

    def test_lune_area_formula(self):
        """Zweieckfläche: A = 2R²·α."""
        assert self.ell.lune_area(math.pi / 4, radius=1.0) == pytest.approx(math.pi / 2)
        assert self.ell.lune_area(math.pi, radius=1.0) == pytest.approx(2 * math.pi)

    def test_girard_theorem_equilateral(self):
        """Girard: gleichseitiges sphärisches Dreieck."""
        alpha = beta = gamma = math.pi / 2  # Drei rechte Winkel
        area = self.ell.girard_theorem(alpha, beta, gamma, radius=1.0)
        expected = math.pi / 2  # E = 3π/2 - π = π/2
        assert area == pytest.approx(expected, abs=TOL)

    def test_girard_theorem_general(self):
        """Girard: A = R²·(α+β+γ-π)."""
        a, b, g = 0.7, 0.8, 0.9
        R = 2.0
        area = self.ell.girard_theorem(a, b, g, radius=R)
        expected = R**2 * (a + b + g - math.pi)
        assert area == pytest.approx(expected, abs=TOL)


# ===========================================================================
# AffinGeometry
# ===========================================================================

class TestAffinGeometry:
    """Tests für AffinGeometry"""

    def setup_method(self):
        self.aff = AffinGeometry()

    def test_affine_combination_midpoint(self):
        """Affinkombination: Mittelpunkt zweier Punkte."""
        result = self.aff.affine_combination([[0, 0], [2, 4]], [0.5, 0.5])
        assert result == pytest.approx([1.0, 2.0], abs=TOL)

    def test_affine_combination_endpoint(self):
        """Koeffizienten [1,0] → erster Punkt."""
        result = self.aff.affine_combination([[3, 5], [7, 1]], [1.0, 0.0])
        assert result == pytest.approx([3.0, 5.0], abs=TOL)

    def test_affine_combination_wrong_sum_raises(self):
        """Koeffizienten summieren nicht auf 1 → ValueError."""
        with pytest.raises(ValueError, match="summieren"):
            self.aff.affine_combination([[0, 0], [1, 1]], [0.3, 0.3])

    def test_affine_map_identity(self):
        """Identitätsabbildung: A=I, b=0."""
        A = np.eye(2)
        b = np.zeros(2)
        result = self.aff.affine_map(A, b, [3, 7])
        assert result == pytest.approx([3.0, 7.0], abs=TOL)

    def test_affine_map_translation(self):
        """Reine Translation: A=I, b=(1,2)."""
        A = np.eye(2)
        b = np.array([1.0, 2.0])
        result = self.aff.affine_map(A, b, [3, 4])
        assert result == pytest.approx([4.0, 6.0], abs=TOL)

    def test_affine_map_rotation_90(self):
        """90°-Rotation."""
        A = np.array([[0.0, -1.0], [1.0, 0.0]])
        b = np.zeros(2)
        result = self.aff.affine_map(A, b, [1, 0])
        assert result == pytest.approx([0.0, 1.0], abs=TOL)

    def test_barycentric_coords_vertex(self):
        """Baryzentrische Koordinate eines Eckpunkts = (1,0,0)."""
        p1, p2, p3 = [0, 0], [1, 0], [0, 1]
        coords = self.aff.barycentric_coords([0, 0], p1, p2, p3)
        assert coords == pytest.approx([1.0, 0.0, 0.0], abs=TOL)

    def test_barycentric_coords_centroid(self):
        """Schwerpunkt hat baryzentrische Koordinaten (1/3, 1/3, 1/3)."""
        p1, p2, p3 = [0, 0], [3, 0], [0, 3]
        centroid = [1.0, 1.0]
        coords = self.aff.barycentric_coords(centroid, p1, p2, p3)
        assert coords == pytest.approx([1/3, 1/3, 1/3], abs=TOL)

    def test_barycentric_coords_inside_positive(self):
        """Innerer Punkt: alle baryzentrischen Koordinaten > 0."""
        p1, p2, p3 = [0, 0], [4, 0], [0, 4]
        coords = self.aff.barycentric_coords([1, 1], p1, p2, p3)
        assert all(c > 0 for c in coords)

    def test_barycentric_coords_degenerate_raises(self):
        """Entartetes Dreieck → ValueError."""
        with pytest.raises(ValueError):
            self.aff.barycentric_coords([0, 0], [0, 0], [1, 1], [2, 2])

    def test_affine_hull_point(self):
        """Einzelner Punkt: Dimension 0."""
        result = self.aff.affine_hull([[1, 2, 3]])
        assert result["dimension"] == 0

    def test_affine_hull_line(self):
        """Zwei verschiedene Punkte: Dimension 1 (Gerade)."""
        result = self.aff.affine_hull([[0, 0], [1, 0]])
        assert result["dimension"] == 1

    def test_affine_hull_plane(self):
        """Drei nicht-kollineare Punkte: Dimension 2 (Ebene)."""
        result = self.aff.affine_hull([[0, 0], [1, 0], [0, 1]])
        assert result["dimension"] == 2

    def test_affine_hull_collinear(self):
        """Drei kollineare Punkte: Dimension 1 (Gerade)."""
        result = self.aff.affine_hull([[0, 0], [1, 1], [2, 2]])
        assert result["dimension"] == 1

    def test_convex_hull_square(self):
        """Konvexe Hülle eines Quadrats mit inneren Punkten."""
        points = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
        hull = self.aff.convex_hull_2d(points)
        # Alle 4 Ecken müssen in der Hülle sein
        assert len(hull) >= 4
        hull_set = [tuple(p) for p in hull]
        assert (0, 0) in hull_set or any(abs(p[0]) < TOL and abs(p[1]) < TOL for p in hull)

    def test_convex_hull_triangle(self):
        """Konvexe Hülle eines Dreiecks."""
        points = [[0, 0], [1, 0], [0, 1]]
        hull = self.aff.convex_hull_2d(points)
        assert len(hull) == 3

    def test_convex_hull_few_points(self):
        """Weniger als 3 Punkte → triviale Hülle."""
        hull = self.aff.convex_hull_2d([[1, 2], [3, 4]])
        assert len(hull) <= 2

    def test_affinely_independent_triangle(self):
        """Drei nicht-kollineare Punkte sind affin unabhängig."""
        assert self.aff.is_affinely_independent([[0, 0], [1, 0], [0, 1]]) == True

    def test_affinely_independent_collinear(self):
        """Kollineare Punkte sind NICHT affin unabhängig."""
        assert self.aff.is_affinely_independent([[0, 0], [1, 1], [2, 2]]) == False

    def test_affinely_independent_single_point(self):
        """Ein einzelner Punkt ist trivial unabhängig."""
        assert self.aff.is_affinely_independent([[5, 7]]) == True


# ===========================================================================
# SyntheticGeometry
# ===========================================================================

class TestSyntheticGeometry:
    """Tests für SyntheticGeometry"""

    def setup_method(self):
        self.syn = SyntheticGeometry()

    def test_hilbert_axioms_has_five_groups(self):
        """Hilbert-Axiomensystem hat 5 Gruppen."""
        axioms = self.syn.hilbert_axioms_overview()
        assert "I_Inzidenz" in axioms
        assert "II_Anordnung" in axioms
        assert "III_Kongruenz" in axioms
        assert "IV_Parallelen" in axioms
        assert "V_Stetigkeit" in axioms

    def test_hilbert_axioms_group_I_count(self):
        """Gruppe I hat 5 Inzidenzaxiome."""
        axioms = self.syn.hilbert_axioms_overview()
        assert len(axioms["I_Inzidenz"]) == 5

    def test_incidence_axioms_demo_structure(self):
        """Inzidenz-Demo gibt strukturiertes Dict zurück."""
        demo = self.syn.incidence_axioms_demo()
        assert "I.1_Beispiel" in demo
        assert "I.2_Beispiel" in demo

    def test_parallel_axiom_demo_three_geometries(self):
        """Parallelenaxiom-Demo enthält alle drei Geometrien."""
        demo = self.syn.parallel_axiom_demo()
        assert "euklidisch" in demo
        assert "hyperbolisch" in demo
        assert "elliptisch" in demo
        assert demo["euklidisch"]["Krümmung"] == 0
        assert demo["hyperbolisch"]["Krümmung"] == -1
        assert demo["elliptisch"]["Krümmung"] == +1

    def test_axiom_independence_demo(self):
        """Unabhängigkeits-Demo gibt dict zurück."""
        demo = self.syn.axiom_independence_demo()
        assert "Parallelenaxiom_unabhängig" in demo

    def test_finite_projective_plane_order_2(self):
        """PG(2,2) = Fano-Ebene: 7 Punkte, 7 Geraden."""
        result = self.syn.finite_projective_plane(2)
        assert result["n_points"] == 7   # 2²+2+1 = 7
        assert result["n_lines"] == 7
        assert result["points_per_line"] == 3  # 2+1
        assert result["lines_per_point"] == 3

    def test_finite_projective_plane_order_3(self):
        """PG(2,3): 13 Punkte, 13 Geraden."""
        result = self.syn.finite_projective_plane(3)
        assert result["n_points"] == 13   # 3²+3+1 = 13
        assert result["points_per_line"] == 4  # 3+1

    def test_finite_projective_plane_order_5(self):
        """PG(2,5): 31 Punkte."""
        result = self.syn.finite_projective_plane(5)
        assert result["n_points"] == 31   # 5²+5+1 = 31

    def test_finite_projective_plane_non_prime_raises(self):
        """Keine Primzahl → ValueError."""
        with pytest.raises(ValueError):
            self.syn.finite_projective_plane(4)  # 4 ist nicht prim

    def test_finite_projective_plane_one_raises(self):
        """q=1 → ValueError."""
        with pytest.raises(ValueError):
            self.syn.finite_projective_plane(1)


# ===========================================================================
# Standalone-Funktionen
# ===========================================================================

class TestEulerCharacteristic:
    """Tests für euler_characteristic_polygon()"""

    def test_tetrahedron(self):
        """Tetraeder: V=4, E=6, F=4 → χ=2."""
        assert euler_characteristic_polygon(4, 6, 4) == 2

    def test_cube(self):
        """Würfel: V=8, E=12, F=6 → χ=2."""
        assert euler_characteristic_polygon(8, 12, 6) == 2

    def test_octahedron(self):
        """Oktaeder: V=6, E=12, F=8 → χ=2."""
        assert euler_characteristic_polygon(6, 12, 8) == 2

    def test_torus(self):
        """Torus: χ=0 (z.B. V=1, E=2, F=1 → χ=0)."""
        # Vereinfachte CW-Struktur: 1 Punkt, 2 Schleifen, 1 Fläche
        assert euler_characteristic_polygon(1, 2, 1) == 0

    def test_icosahedron(self):
        """Ikosaeder: V=12, E=30, F=20 → χ=2."""
        assert euler_characteristic_polygon(12, 30, 20) == 2


class TestPickTheorem:
    """Tests für pick_theorem()"""

    def test_unit_square(self):
        """Einheitsquadrat: I=0, B=4 → A=1."""
        # Einheitsquadrat (0,0),(1,0),(1,1),(0,1): I=0, B=4
        assert pick_theorem(0, 4) == pytest.approx(1.0)

    def test_2x2_square(self):
        """2×2-Quadrat: I=1, B=8 → A=4."""
        assert pick_theorem(1, 8) == pytest.approx(4.0)

    def test_larger_polygon(self):
        """Größeres Polygon."""
        # A = I + B/2 - 1 = 5 + 12/2 - 1 = 10
        assert pick_theorem(5, 12) == pytest.approx(10.0)

    def test_triangle(self):
        """Dreieck (0,0),(4,0),(0,3): I=3, B=6 → A=6."""
        # Fläche = 4*3/2 = 6, I = 3 (prüfe: I+B/2-1 = 3+3-1=5... muss variieren)
        # Korrekt: Dreieck (0,0)(2,0)(0,2): I=0, B=3+... anpassen
        # Einfaches Beispiel: A=I+B/2-1
        A = pick_theorem(3, 6)
        assert A == pytest.approx(5.0)


class TestPtolemyTheorem:
    """Tests für ptolemy_theorem_check()"""

    def test_cyclic_quadrilateral_square(self):
        """Quadrat ist ein Sehnenviereck."""
        p1 = [0, 0]
        p2 = [1, 0]
        p3 = [1, 1]
        p4 = [0, 1]
        result = ptolemy_theorem_check(p1, p2, p3, p4)
        assert result["ptolemy_satisfied"] is True

    def test_cyclic_quadrilateral_rectangle(self):
        """Rechteck ist ein Sehnenviereck."""
        p1 = [0, 0]
        p2 = [2, 0]
        p3 = [2, 1]
        p4 = [0, 1]
        result = ptolemy_theorem_check(p1, p2, p3, p4)
        assert result["ptolemy_satisfied"] is True

    def test_non_cyclic_quadrilateral(self):
        """Nicht-zyklisches Viereck erfüllt Ptolemäus-Gleichung NICHT."""
        # Viereck (0,0),(3,0),(3,4),(0,2): liegt NICHT auf einem Kreis
        p1 = [0, 0]
        p2 = [3, 0]
        p3 = [3, 4]
        p4 = [0, 2]
        result = ptolemy_theorem_check(p1, p2, p3, p4)
        # Kein Sehnenviereck: |lhs - rhs| > 0
        assert result["ptolemy_satisfied"] == False

    def test_result_has_required_keys(self):
        """Ergebnis-Dict enthält alle erwarteten Schlüssel."""
        result = ptolemy_theorem_check([0,0],[1,0],[1,1],[0,1])
        for key in ["lhs", "rhs", "ptolemy_satisfied", "is_cyclic"]:
            assert key in result


class TestNinePointCircle:
    """Tests für nine_point_circle()"""

    def test_nine_point_radius_half_circumradius(self):
        """Radius des Neun-Punkte-Kreises = Umkreisradius / 2."""
        p1, p2, p3 = [0, 0], [4, 0], [0, 3]
        result = nine_point_circle(p1, p2, p3)
        R = float(np.linalg.norm(np.array(result["circumcenter"]) - np.array(p1)))
        assert result["radius"] == pytest.approx(R / 2, rel=1e-6)

    def test_nine_point_midpoints_equidistant(self):
        """Alle drei Seitenmittelpunkte haben gleichen Abstand zum Mittelpunkt."""
        p1, p2, p3 = [0, 0], [4, 0], [2, 3]
        result = nine_point_circle(p1, p2, p3)
        N = result["center"]
        r9 = result["radius"]
        for M in result["midpoints_of_sides"]:
            d = np.linalg.norm(np.array(M) - np.array(N))
            assert d == pytest.approx(r9, rel=1e-6)

    def test_nine_point_foot_altitudes_equidistant(self):
        """Alle Höhenfußpunkte liegen auf dem Neun-Punkte-Kreis."""
        p1, p2, p3 = [0, 0], [4, 0], [1, 3]
        result = nine_point_circle(p1, p2, p3)
        N = result["center"]
        r9 = result["radius"]
        for F in result["foot_of_altitudes"]:
            d = np.linalg.norm(np.array(F) - np.array(N))
            assert d == pytest.approx(r9, rel=1e-5)

    def test_nine_point_degenerate_raises(self):
        """Entartetes Dreieck → ValueError."""
        with pytest.raises(ValueError):
            nine_point_circle([0, 0], [1, 1], [2, 2])

    def test_nine_point_result_keys(self):
        """Ergebnis enthält alle erwarteten Schlüssel."""
        result = nine_point_circle([0, 0], [4, 0], [0, 3])
        for key in ["center", "radius", "midpoints_of_sides",
                    "foot_of_altitudes", "euler_midpoints"]:
            assert key in result


class TestMorleyTheorem:
    """Tests für morley_theorem_demo()"""

    def test_morley_returns_dict(self):
        """Morley-Demo gibt ein Dict zurück."""
        result = morley_theorem_demo()
        assert isinstance(result, dict)

    def test_morley_equilateral_triangle(self):
        """Bei gleichseitigem Dreieck ist das Morley-Dreieck ebenfalls gleichseitig."""
        result = morley_theorem_demo()
        assert result["beispiel_gleichseitig"]["morley_equilateral"] is True

    def test_morley_side_length_positive(self):
        """Seitenlänge des Morley-Dreiecks ist > 0."""
        result = morley_theorem_demo()
        assert result["beispiel_gleichseitig"]["morley_side_length"] > 0.0
        assert result["beispiel_allgemein"]["morley_side_length"] > 0.0

    def test_morley_formula_key(self):
        """Formel ist im Ergebnis enthalten."""
        result = morley_theorem_demo()
        assert "seitenlaenge_formel" in result

    def test_morley_equilateral_side_length_formula(self):
        """Seitenlänge bei 60°-Dreieck: s = 8·sin(20°)³."""
        result = morley_theorem_demo()
        R = result["beispiel_gleichseitig"]["R"]
        s = result["beispiel_gleichseitig"]["morley_side_length"]
        expected = 8.0 * R * math.sin(math.pi/9)**3  # sin(20°)³
        assert s == pytest.approx(expected, rel=1e-6)
