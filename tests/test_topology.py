"""
@file test_topology.py
@brief Tests für das Topologie-Modul (topology.py).
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import pytest
from src.topology import (
    # Standardmetriken
    euclidean_metric,
    manhattan_metric,
    chebyshev_metric,
    discrete_metric,
    p_norm_metric,
    # Metrischer Raum
    MetricSpace,
    # Topologische Eigenschaften
    is_connected,
    hausdorff_distance,
    compute_diameter,
    is_compact_discrete,
    # Parametrische Kurven
    ParametricCurve,
    circle_curve,
    lissajous_curve,
    helix_curve,
    # Simpliziale Homologie
    euler_characteristic_polygon,
    genus_from_euler,
    betti_numbers_graph,
    # Fraktale Dimensionen
    box_counting_dimension,
    hausdorff_dimension_cantor,
    sierpinski_dimension,
)


# ============================================================
# Tests: Euklidische Metrik
# ============================================================

class TestEuclideanMetric:
    """Tests für die euklidische Metrik."""

    def test_basic_3_4_5(self):
        """3-4-5-Dreieck: Hypotenuse = 5."""
        assert abs(euclidean_metric([0, 0], [3, 4]) - 5.0) < 1e-10

    def test_same_point_is_zero(self):
        """Abstand eines Punktes zu sich selbst ist 0."""
        assert euclidean_metric([1, 2, 3], [1, 2, 3]) == 0.0

    def test_symmetry(self):
        """d(x,y) = d(y,x)."""
        x, y = [1, 2], [4, 6]
        assert abs(euclidean_metric(x, y) - euclidean_metric(y, x)) < 1e-15

    def test_1d(self):
        """Euklidisch in 1D = absoluter Betrag."""
        assert abs(euclidean_metric([0], [5]) - 5.0) < 1e-10

    def test_3d(self):
        """Euklidisch in 3D: sqrt(1+4+9) = sqrt(14)."""
        assert abs(euclidean_metric([0, 0, 0], [1, 2, 3]) - math.sqrt(14)) < 1e-10

    def test_non_negative(self):
        """Nichtnegativität."""
        assert euclidean_metric([3, -5], [-1, 2]) >= 0


# ============================================================
# Tests: Manhattan-Metrik
# ============================================================

class TestManhattanMetric:
    """Tests für die Manhattan-Metrik."""

    def test_basic(self):
        """Manhattan-Abstand von (0,0) nach (3,4) = 7."""
        assert manhattan_metric([0, 0], [3, 4]) == 7.0

    def test_same_point(self):
        """Abstand zu sich selbst = 0."""
        assert manhattan_metric([5, 5], [5, 5]) == 0.0

    def test_1d(self):
        """1D: |x - y|."""
        assert manhattan_metric([0], [10]) == 10.0

    def test_symmetry(self):
        """Symmetrie."""
        x, y = [2, 3], [5, 8]
        assert manhattan_metric(x, y) == manhattan_metric(y, x)

    def test_always_geq_euclidean(self):
        """Manhattan ≥ Euklidisch (allgemeine Ungleichung)."""
        x, y = [1, 2, 3], [4, 6, 8]
        assert manhattan_metric(x, y) >= euclidean_metric(x, y) - 1e-10


# ============================================================
# Tests: Chebyshev-Metrik
# ============================================================

class TestChebyshevMetric:
    """Tests für die Chebyshev-Metrik."""

    def test_basic(self):
        """Chebyshev: max(|3-0|, |4-0|) = 4."""
        assert chebyshev_metric([0, 0], [3, 4]) == 4.0

    def test_same_point(self):
        assert chebyshev_metric([1, 1], [1, 1]) == 0.0

    def test_symmetry(self):
        x, y = [2, 7], [5, 3]
        assert chebyshev_metric(x, y) == chebyshev_metric(y, x)

    def test_always_leq_manhattan(self):
        """Chebyshev ≤ Manhattan."""
        x, y = [1, 2], [4, 6]
        assert chebyshev_metric(x, y) <= manhattan_metric(x, y) + 1e-10


# ============================================================
# Tests: Diskrete Metrik
# ============================================================

class TestDiscreteMetric:
    """Tests für die diskrete Metrik."""

    def test_equal_elements(self):
        """Gleiche Elemente: Abstand 0."""
        assert discrete_metric(5, 5) == 0.0

    def test_unequal_elements(self):
        """Verschiedene Elemente: Abstand 1."""
        assert discrete_metric(5, 6) == 1.0

    def test_strings(self):
        """Funktioniert auch für Strings."""
        assert discrete_metric('a', 'a') == 0.0
        assert discrete_metric('a', 'b') == 1.0


# ============================================================
# Tests: Lp-Metrik
# ============================================================

class TestPNormMetric:
    """Tests für die Lp-Metrik."""

    def test_p2_equals_euclidean(self):
        """p=2 entspricht der euklidischen Metrik."""
        x, y = [0, 0], [3, 4]
        assert abs(p_norm_metric(x, y, 2) - euclidean_metric(x, y)) < 1e-10

    def test_p1_equals_manhattan(self):
        """p=1 entspricht der Manhattan-Metrik."""
        x, y = [0, 0], [3, 4]
        assert abs(p_norm_metric(x, y, 1) - manhattan_metric(x, y)) < 1e-10

    def test_p_less_than_1_raises(self):
        """p < 1 ist ungültig → ValueError."""
        with pytest.raises(ValueError):
            p_norm_metric([0], [1], 0.5)

    def test_p3(self):
        """p=3: (|1|^3 + |2|^3)^{1/3} = (9)^{1/3}."""
        result = p_norm_metric([0, 0], [1, 2], 3)
        expected = (1 + 8) ** (1 / 3)
        assert abs(result - expected) < 1e-10


# ============================================================
# Tests: MetricSpace
# ============================================================

class TestMetricSpace:
    """Tests für die MetricSpace-Klasse."""

    def test_verify_axioms_euclidean(self):
        """Euklidische Metrik erfüllt alle 4 Axiome."""
        space = MetricSpace(euclidean_metric)
        points = [[0, 0], [1, 0], [0, 1]]
        result = space.verify_metric_axioms(points)

        assert result['non_negativity'] is True
        assert result['identity'] is True
        assert result['symmetry'] is True
        assert result['triangle_inequality'] is True

    def test_verify_axioms_all_true(self):
        """Alle Axiome für Manhattan-Metrik mit mehreren Punkten."""
        space = MetricSpace(manhattan_metric)
        points = [[0, 0], [1, 0], [0, 1], [1, 1], [-1, 2]]
        result = space.verify_metric_axioms(points)
        assert all(result.values())

    def test_open_ball_euclidean(self):
        """Offene Kugel um Ursprung mit Radius 2."""
        space = MetricSpace(euclidean_metric)
        points = [[0, 0], [1, 0], [1, 1], [3, 0]]  # [3,0] liegt außerhalb
        ball = space.open_ball([0, 0], 2.0, points)
        assert [0, 0] in ball
        assert [1, 0] in ball
        assert [3, 0] not in ball

    def test_open_ball_boundary_excluded(self):
        """Randpunkt liegt NICHT in der offenen Kugel."""
        space = MetricSpace(euclidean_metric)
        points = [[0, 0], [1, 0]]
        # Genau auf dem Rand: d = 1.0, Radius = 1.0 → nicht drin
        ball = space.open_ball([0, 0], 1.0, points)
        assert [1, 0] not in ball

    def test_cauchy_sequence_convergent(self):
        """Konvergente Folge ist Cauchy."""
        space = MetricSpace(euclidean_metric)
        # x_n = 1/n → 0 (Cauchy)
        seq = [[1.0 / n, 0] for n in range(1, 101)]
        assert space.is_cauchy_sequence(seq, tol=0.1) is True

    def test_cauchy_sequence_divergent(self):
        """Divergente Folge ist nicht Cauchy."""
        space = MetricSpace(euclidean_metric)
        # x_n = n → ∞ (nicht Cauchy)
        seq = [[float(n), 0] for n in range(100)]
        assert space.is_cauchy_sequence(seq, tol=1e-5) is False


# ============================================================
# Tests: Topologische Eigenschaften
# ============================================================

class TestTopologicalProperties:
    """Tests für topologische Eigenschaften."""

    def test_is_connected_single_point(self):
        """Einzelner Punkt ist immer zusammenhängend."""
        assert is_connected([[0, 0]], euclidean_metric, 1.0) is True

    def test_is_connected_true(self):
        """Nahe beieinander liegende Punkte: zusammenhängend."""
        points = [[0, 0], [0.5, 0], [1.0, 0]]
        assert is_connected(points, euclidean_metric, 0.6) is True

    def test_is_connected_false(self):
        """Weit auseinanderliegende Punkte: nicht zusammenhängend."""
        points = [[0, 0], [10, 0]]
        assert is_connected(points, euclidean_metric, 1.0) is False

    def test_hausdorff_distance_same_sets(self):
        """Hausdorff-Abstand einer Menge zu sich selbst = 0."""
        pts = [[0], [1], [2]]
        assert hausdorff_distance(pts, pts, euclidean_metric) == 0.0

    def test_hausdorff_distance_shifted(self):
        """Hausdorff-Abstand zwischen {0} und {1} = 1."""
        a = [[0]]
        b = [[1]]
        assert abs(hausdorff_distance(a, b, euclidean_metric) - 1.0) < 1e-10

    def test_hausdorff_distance_symmetry(self):
        """Hausdorff-Abstand ist symmetrisch."""
        a = [[0], [1]]
        b = [[0.5], [1.5]]
        assert abs(hausdorff_distance(a, b, euclidean_metric) -
                   hausdorff_distance(b, a, euclidean_metric)) < 1e-10

    def test_compute_diameter_empty(self):
        """Leere Menge hat Durchmesser 0."""
        assert compute_diameter([], euclidean_metric) == 0.0

    def test_compute_diameter_square(self):
        """Einheitsquadrat: Durchmesser = sqrt(2)."""
        pts = [[0, 0], [1, 0], [1, 1], [0, 1]]
        assert abs(compute_diameter(pts, euclidean_metric) - math.sqrt(2)) < 1e-10

    def test_is_compact_discrete(self):
        """Endliche diskrete Menge ist immer kompakt."""
        pts = [[i, 0] for i in range(10)]
        assert is_compact_discrete(pts, euclidean_metric, 1.0) is True


# ============================================================
# Tests: ParametricCurve
# ============================================================

class TestParametricCurve:
    """Tests für parametrische Kurven."""

    def test_evaluate_circle(self):
        """Kreis: γ(0) = (1, 0)."""
        curve = circle_curve(radius=1.0)
        pt = curve.evaluate(0.0)
        assert abs(pt[0] - 1.0) < 1e-10
        assert abs(pt[1] - 0.0) < 1e-10

    def test_evaluate_circle_half(self):
        """Kreis: γ(π) = (-1, 0)."""
        curve = circle_curve(radius=1.0)
        pt = curve.evaluate(math.pi)
        assert abs(pt[0] - (-1.0)) < 1e-10
        assert abs(pt[1] - 0.0) < 1e-9

    def test_arc_length_circle(self):
        """Kreisumfang ≈ 2π."""
        curve = circle_curve(radius=1.0)
        length = curve.arc_length(n=10000)
        assert abs(length - 2 * math.pi) < 1e-3

    def test_arc_length_circle_radius2(self):
        """Kreisumfang mit r=2: ≈ 4π."""
        curve = circle_curve(radius=2.0)
        length = curve.arc_length(n=10000)
        assert abs(length - 4 * math.pi) < 1e-3

    def test_is_closed_circle(self):
        """Kreis ist eine geschlossene Kurve."""
        curve = circle_curve()
        assert curve.is_closed() is True

    def test_is_closed_helix(self):
        """Helix über [0, 2π]: Anfangs- und Endpunkt identisch in x,y,
        aber nicht in z → nicht geschlossen."""
        curve = helix_curve(radius=1.0, pitch=1.0)
        # z(0)=0, z(2π)=1 → nicht geschlossen
        assert curve.is_closed() is False

    def test_curvature_circle(self):
        """Krümmung des Einheitskreises = 1.0 überall."""
        curve = circle_curve(radius=1.0)
        kappa = curve.curvature(1.0)
        assert abs(kappa - 1.0) < 1e-4

    def test_curvature_circle_radius2(self):
        """Krümmung des Kreises mit r=2: κ = 1/2."""
        curve = circle_curve(radius=2.0)
        kappa = curve.curvature(1.0)
        assert abs(kappa - 0.5) < 1e-4

    def test_winding_number_inside(self):
        """Einheitskreis umläuft den Ursprung einmal."""
        curve = circle_curve(radius=1.0)
        wn = curve.winding_number([0.0, 0.0])
        assert wn == 1

    def test_winding_number_outside(self):
        """Einheitskreis umläuft Punkt außerhalb 0 mal."""
        curve = circle_curve(radius=1.0)
        wn = curve.winding_number([5.0, 0.0])
        assert wn == 0

    def test_lissajous_closed(self):
        """Lissajous-Figur über [0, 2π] ist (annähernd) geschlossen."""
        curve = lissajous_curve(3, 2)
        # sin(3*0+0)=0, sin(3*2π+0)=0; sin(2*0)=0, sin(2*2π)=0 → geschlossen
        assert curve.is_closed(tol=1e-8) is True

    def test_helix_arc_length(self):
        """Helix-Bogenlänge ≈ sqrt((2π)² + 1²) pro Windung."""
        # Bogenlänge der Helix: L = sqrt(r²·(2π)² + pitch²) für eine Windung
        r, pitch = 1.0, 1.0
        expected = math.sqrt((r * 2 * math.pi) ** 2 + pitch ** 2)
        curve = helix_curve(radius=r, pitch=pitch)
        length = curve.arc_length(n=10000)
        assert abs(length - expected) < 1e-2


# ============================================================
# Tests: Simpliziale Homologie
# ============================================================

class TestSimplicialHomology:
    """Tests für Euler-Charakteristik und Betti-Zahlen."""

    def test_euler_cube(self):
        """Würfel: V=8, E=12, F=6 → χ = 8-12+6 = 2."""
        assert euler_characteristic_polygon(8, 12, 6) == 2

    def test_euler_tetrahedron(self):
        """Tetraeder: V=4, E=6, F=4 → χ = 4-6+4 = 2."""
        assert euler_characteristic_polygon(4, 6, 4) == 2

    def test_euler_torus(self):
        """Torus: χ = 0 (z.B. V=1, E=2, F=1 für Zellzerlegung)."""
        # Standard Torus-Triangulierung: V-E+F = 0
        # Einfachste Beispiel: V=1, E=2, F=1 → χ=0
        assert euler_characteristic_polygon(1, 2, 1) == 0

    def test_genus_sphere(self):
        """Sphäre: χ=2 → g=0."""
        assert genus_from_euler(2, orientable=True) == 0

    def test_genus_torus(self):
        """Torus: χ=0 → g=1."""
        assert genus_from_euler(0, orientable=True) == 1

    def test_genus_double_torus(self):
        """Doppeltorus: χ=-2 → g=2."""
        assert genus_from_euler(-2, orientable=True) == 2

    def test_genus_non_orientable(self):
        """Projektive Ebene (nicht-orientierbar): χ=1 → g=1."""
        assert genus_from_euler(1, orientable=False) == 1

    def test_betti_path_graph(self):
        """Pfadgraph P_3 (3 Knoten, 2 Kanten): β0=1, β1=0."""
        # Adjazenzmatrix: 0-1-2
        adj = [
            [0, 1, 0],
            [1, 0, 1],
            [0, 1, 0],
        ]
        result = betti_numbers_graph(adj)
        assert result['beta_0'] == 1
        assert result['beta_1'] == 0

    def test_betti_cycle_graph(self):
        """Dreieck C_3: β0=1, β1=1 (ein unabhängiger Zyklus)."""
        adj = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        result = betti_numbers_graph(adj)
        assert result['beta_0'] == 1
        assert result['beta_1'] == 1

    def test_betti_disconnected(self):
        """Zwei isolierte Knoten: β0=2, β1=0."""
        adj = [
            [0, 0],
            [0, 0],
        ]
        result = betti_numbers_graph(adj)
        assert result['beta_0'] == 2
        assert result['beta_1'] == 0


# ============================================================
# Tests: Fraktale Dimensionen
# ============================================================

class TestFractalDimensions:
    """Tests für fraktale Dimensionen."""

    def test_hausdorff_cantor(self):
        """Hausdorff-Dimension der Cantor-Menge ≈ 0.6309."""
        d = hausdorff_dimension_cantor()
        assert abs(d - math.log(2) / math.log(3)) < 1e-15
        assert abs(d - 0.6309) < 1e-4

    def test_sierpinski_dimension(self):
        """Hausdorff-Dimension des Sierpinski-Dreiecks ≈ 1.585."""
        d = sierpinski_dimension()
        assert abs(d - math.log(3) / math.log(2)) < 1e-15
        assert abs(d - 1.585) < 1e-3

    def test_box_counting_returns_float(self):
        """box_counting_dimension gibt einen float zurück."""
        # Einfacher Test: Punkte auf einer Linie
        points = [[i * 0.1, 0] for i in range(100)]
        eps_values = [0.5, 0.25, 0.1, 0.05]
        result = box_counting_dimension(points, eps_values)
        assert isinstance(result, float)

    def test_box_counting_line(self):
        """Punkte auf einer Linie haben Dimension ≈ 1."""
        # Viele Punkte auf der x-Achse
        points = [[i * 0.01, 0] for i in range(1000)]
        eps_values = [0.2, 0.1, 0.05, 0.02, 0.01]
        d = box_counting_dimension(points, eps_values)
        # Dimension einer Linie ist ≈ 1
        assert 0.8 < d < 1.2

    def test_box_counting_plane(self):
        """Punkte auf einem 2D-Gitter haben Dimension ≈ 2."""
        # Punkte auf einem Quadratgitter
        points = [[i * 0.1, j * 0.1] for i in range(20) for j in range(20)]
        eps_values = [0.5, 0.3, 0.2, 0.1, 0.05]
        d = box_counting_dimension(points, eps_values)
        # Dimension eines 2D-Gitters ≈ 2
        assert 1.5 < d < 2.5

    def test_box_counting_empty(self):
        """Leere Punktmenge gibt 0.0 zurück."""
        result = box_counting_dimension([], [0.1, 0.05])
        assert result == 0.0

    def test_cantor_less_than_1(self):
        """Cantor-Dimension < 1 (gebrochene Dimension)."""
        assert hausdorff_dimension_cantor() < 1.0

    def test_sierpinski_between_1_and_2(self):
        """Sierpinski-Dimension liegt zwischen 1 und 2."""
        d = sierpinski_dimension()
        assert 1.0 < d < 2.0


# ============================================================
# Tests: Metrische Axiome (detailliert)
# ============================================================

class TestMetricAxioms:
    """Detaillierte Tests der vier Metrik-Axiome."""

    def test_non_negativity_all_metrics(self):
        """Alle Standardmetriken sind nicht-negativ."""
        x, y = [1, 2, 3], [4, 6, 8]
        assert euclidean_metric(x, y) >= 0
        assert manhattan_metric(x, y) >= 0
        assert chebyshev_metric(x, y) >= 0

    def test_identity_euclidean(self):
        """d(x,x) = 0 für euklidische Metrik."""
        pt = [1, 2, 3]
        assert euclidean_metric(pt, pt) == 0.0

    def test_triangle_inequality_euclidean(self):
        """Dreiecksungleichung: d(x,z) ≤ d(x,y) + d(y,z)."""
        x, y, z = [0, 0], [3, 0], [3, 4]
        d_xz = euclidean_metric(x, z)  # sqrt(9+16) = 5
        d_xy = euclidean_metric(x, y)  # 3
        d_yz = euclidean_metric(y, z)  # 4
        assert d_xz <= d_xy + d_yz + 1e-10
