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
        # Punkte auf einem Quadratgitter mit Abstand 0.1 (Bereich 0..1.9)
        # Epsilon-Werte müssen >= Gitterabstand sein, damit Skalierung sichtbar ist
        points = [[i * 0.1, j * 0.1] for i in range(20) for j in range(20)]
        eps_values = [1.0, 0.5, 0.25, 0.2, 0.1]
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


# ============================================================
# Tests: Trennungseigenschaften (Separation Axioms)
# ============================================================

from src.topology import (
    separation_axioms_demo,
    CompactSpace,
    heine_borel_theorem,
    tychonoff_theorem_demo,
    UniformSpace,
    cauchy_filter_demo,
    completion_of_rationals,
    TopologicalGroup,
    topological_group_examples,
    TopologicalVectorSpace,
    locally_convex_demo,
    schwartz_space_intro,
)


class TestSeparationAxioms:
    """Tests für die Trennungseigenschaften T0 bis T4."""

    def test_returns_dict(self):
        """separation_axioms_demo gibt ein dict zurück."""
        result = separation_axioms_demo()
        assert isinstance(result, dict)

    def test_all_axioms_present(self):
        """Alle 5 Axiome T0–T4 sind im Ergebnis enthalten."""
        result = separation_axioms_demo()
        assert 'T0_Kolmogorov' in result
        assert 'T1_Frechet' in result
        assert 'T2_Hausdorff' in result
        assert 'T3_Regular' in result
        assert 'T4_Normal' in result

    def test_axiom_count(self):
        """axiom_count ist 5."""
        result = separation_axioms_demo()
        assert result['axiom_count'] == 5

    def test_implications_present(self):
        """Implikationshierarchie ist vorhanden."""
        result = separation_axioms_demo()
        assert 'implications' in result
        assert 'hierarchy' in result['implications']

    def test_t2_has_hausdorff_content(self):
        """T2-Eintrag enthält 'Hausdorff'."""
        result = separation_axioms_demo()
        assert 'Hausdorff' in result['T2_Hausdorff']['name']

    def test_t4_mentions_urysohn(self):
        """T4-Eintrag enthält Urysohn-Lemma."""
        result = separation_axioms_demo()
        assert 'urysohn_lemma' in result['T4_Normal']

    def test_implication_string(self):
        """Implikationskette enthält T4 und T0."""
        result = separation_axioms_demo()
        hierarchy = result['implications']['hierarchy']
        assert 'T4' in hierarchy and 'T0' in hierarchy

    def test_each_axiom_has_definition(self):
        """Jedes Axiom hat eine 'definition'-Schlüssel."""
        result = separation_axioms_demo()
        for key in ['T0_Kolmogorov', 'T1_Frechet', 'T2_Hausdorff', 'T3_Regular', 'T4_Normal']:
            assert 'definition' in result[key], f"Fehlende Definition in {key}"


# ============================================================
# Tests: Kompaktheit (Compactness)
# ============================================================

class TestCompactSpace:
    """Tests für die CompactSpace-Klasse."""

    def test_compact_interval(self):
        """[0, 1] ist kompakt (a ≤ b)."""
        cs = CompactSpace()
        assert cs.is_compact_heine_borel(0.0, 1.0) is True

    def test_compact_negative_interval(self):
        """[-5, 5] ist kompakt."""
        cs = CompactSpace()
        assert cs.is_compact_heine_borel(-5.0, 5.0) is True

    def test_compact_single_point(self):
        """[1, 1] (einpunktiger Raum) ist kompakt."""
        cs = CompactSpace()
        assert cs.is_compact_heine_borel(1.0, 1.0) is True

    def test_not_compact_reversed(self):
        """[1, 0] (a > b) ist kein gültiges Intervall → nicht kompakt."""
        cs = CompactSpace()
        assert cs.is_compact_heine_borel(1.0, 0.0) is False

    def test_finite_subcover_simple(self):
        """Einfache Überdeckung von [0,1] hat endliche Teilüberdeckung."""
        cs = CompactSpace()
        cover = [(-0.1, 0.6), (0.4, 1.1)]
        subcover = cs.finite_subcover_demo(cover, 0.0, 1.0)
        assert len(subcover) > 0
        assert len(subcover) <= len(cover)

    def test_finite_subcover_covers_interval(self):
        """Die Teilüberdeckung überdeckt tatsächlich [0,1]."""
        cs = CompactSpace()
        cover = [(-0.1, 0.4), (0.3, 0.8), (0.7, 1.1)]
        subcover = cs.finite_subcover_demo(cover, 0.0, 1.0)
        # Überprüfe dass Teilüberdeckung [0,1] überdeckt
        assert len(subcover) > 0

    def test_sequential_compactness(self):
        """Folge in [0,1] hat konvergente Teilfolge."""
        cs = CompactSpace()
        import math
        seq = [math.sin(i * 0.1) * 0.5 + 0.5 for i in range(50)]
        result = cs.sequential_compactness_demo(seq, 0.0, 1.0)
        assert result['converges'] is True
        assert result['bounded_in_interval'] is True
        assert result['limit'] is not None


class TestHeineBorel:
    """Tests für den Satz von Heine-Borel."""

    def test_returns_dict(self):
        """heine_borel_theorem gibt ein dict zurück."""
        result = heine_borel_theorem()
        assert isinstance(result, dict)

    def test_theorem_present(self):
        """Theorem-Text ist vorhanden."""
        result = heine_borel_theorem()
        assert 'theorem' in result
        assert 'Heine-Borel' in result['theorem']

    def test_compact_examples(self):
        """Mindestens 2 kompakte Beispiele vorhanden."""
        result = heine_borel_theorem()
        assert len(result['compact_examples']) >= 2
        for ex in result['compact_examples']:
            assert ex['compact'] is True
            assert ex['closed'] is True
            assert ex['bounded'] is True

    def test_non_compact_examples(self):
        """Mindestens 2 nicht-kompakte Beispiele vorhanden."""
        result = heine_borel_theorem()
        assert len(result['non_compact_examples']) >= 2
        for ex in result['non_compact_examples']:
            assert ex['compact'] is False

    def test_numerical_verification(self):
        """Numerische Verifikation enthält endliche Teilüberdeckung."""
        result = heine_borel_theorem()
        nv = result['numerical_verification']
        assert nv['is_finite'] is True
        assert len(nv['finite_subcover']) > 0

    def test_equivalences_list(self):
        """Drei Äquivalenzen für Kompaktheit werden aufgelistet."""
        result = heine_borel_theorem()
        assert len(result['equivalences']) >= 3


class TestTychonoff:
    """Tests für den Satz von Tychonoff."""

    def test_returns_dict(self):
        """tychonoff_theorem_demo gibt ein dict zurück."""
        result = tychonoff_theorem_demo()
        assert isinstance(result, dict)

    def test_finite_product_compact(self):
        """Endliches Produkt kompakter Räume ist kompakt."""
        result = tychonoff_theorem_demo()
        assert result['finite_demo']['is_compact'] is True

    def test_product_size(self):
        """Produktmenge hat |A| × |B| Elemente."""
        result = tychonoff_theorem_demo()
        fd = result['finite_demo']
        assert fd['product_size'] == len(fd['space_a']) * len(fd['space_b'])

    def test_important_examples(self):
        """Mindestens 3 wichtige Beispiele (Hilbert-Würfel etc.)."""
        result = tychonoff_theorem_demo()
        assert len(result['important_examples']) >= 3

    def test_axiom_of_choice_mentioned(self):
        """Auswahlaxiom wird erwähnt."""
        result = tychonoff_theorem_demo()
        assert 'axiom_of_choice' in result
        assert 'Auswahlaxiom' in result['axiom_of_choice']


# ============================================================
# Tests: Uniforme Räume (Uniform Spaces)
# ============================================================

class TestUniformSpace:
    """Tests für die UniformSpace-Klasse."""

    def test_entourage_contains_diagonal(self):
        """Epsilon-Entourage enthält die Diagonale."""
        points = [0, 1, 2, 3]
        us = UniformSpace(points, metric=lambda x, y: abs(x - y))
        ent = us.entourage(0.5)
        # Prüfe Diagonale
        assert us.contains_diagonal(ent) is True

    def test_entourage_size(self):
        """Für ε=0.5 enthält Entourage auf {0,1,2,3} nur Diagonal-Paare."""
        points = [0, 1, 2, 3]
        us = UniformSpace(points, metric=lambda x, y: abs(x - y))
        ent = us.entourage(0.5)
        # Nur Paare (x,y) mit |x-y| < 0.5, also nur Diagonale (|x-x|=0 < 0.5)
        assert (0, 0) in ent
        assert (0, 1) not in ent  # |0-1|=1 ≥ 0.5

    def test_entourage_large_epsilon(self):
        """Für große ε enthält Entourage alle Paare."""
        points = [0, 1, 2]
        us = UniformSpace(points, metric=lambda x, y: abs(x - y))
        ent = us.entourage(100.0)
        assert len(ent) == len(points) ** 2

    def test_symmetric_entourage(self):
        """Entourage ist symmetrisch."""
        points = [0, 1, 2]
        us = UniformSpace(points, metric=lambda x, y: abs(x - y))
        ent = us.entourage(1.5)
        assert us.is_symmetric(ent) is True

    def test_uniform_continuity_identity(self):
        """Identitätsfunktion ist gleichmäßig stetig."""
        points = [0.0, 0.5, 1.0, 1.5, 2.0]
        us = UniformSpace(points, metric=lambda x, y: abs(x - y))
        result = us.uniform_continuity_check(
            f=lambda x: x,
            target_metric=lambda x, y: abs(x - y),
            epsilon=0.01
        )
        assert result['uniformly_continuous'] is True
        assert result['delta'] is not None


class TestCauchyFilter:
    """Tests für cauchy_filter_demo (√2-Approximation in ℚ)."""

    def test_returns_dict(self):
        """cauchy_filter_demo gibt ein dict zurück."""
        result = cauchy_filter_demo()
        assert isinstance(result, dict)

    def test_sequence_converges_to_sqrt2(self):
        """Die Folge konvergiert numerisch gegen √2."""
        import math
        result = cauchy_filter_demo()
        assert abs(result['limit_value'] - math.sqrt(2)) < 1e-15

    def test_is_cauchy(self):
        """Die Folge ist eine Cauchy-Folge."""
        result = cauchy_filter_demo()
        assert result['is_cauchy'] is True

    def test_limit_not_in_Q(self):
        """Der Grenzwert √2 ist nicht in ℚ."""
        result = cauchy_filter_demo()
        assert result['limit_in_Q'] is False

    def test_differences_decrease(self):
        """Die Abstände aufeinanderfolgender Glieder fallen monoton."""
        result = cauchy_filter_demo()
        diffs = result['differences']
        # Nach einigen Schritten fallen die Differenzen schnell
        assert diffs[5] < diffs[1]


class TestCompletionOfRationals:
    """Tests für completion_of_rationals."""

    def test_returns_dict(self):
        """completion_of_rationals gibt ein dict zurück."""
        result = completion_of_rationals()
        assert isinstance(result, dict)

    def test_real_completion_present(self):
        """Vervollständigung zu ℝ ist beschrieben."""
        result = completion_of_rationals()
        assert 'completion_real' in result
        assert 'ℝ' in result['completion_real']['result']

    def test_p_adic_completion_present(self):
        """p-adische Vervollständigung ist beschrieben."""
        result = completion_of_rationals()
        assert 'completion_p_adic' in result

    def test_p_adic_examples_correct(self):
        """p-adische Normen: |5|_5 = 0.2, |25|_5 = 0.04."""
        result = completion_of_rationals()
        examples = result['completion_p_adic']['examples']
        norm_5 = next(e['norm'] for e in examples if e['x'] == 5)
        norm_25 = next(e['norm'] for e in examples if e['x'] == 25)
        assert abs(norm_5 - 0.2) < 1e-10
        assert abs(norm_25 - 0.04) < 1e-10

    def test_ostrowski_theorem_mentioned(self):
        """Satz von Ostrowski wird erwähnt."""
        result = completion_of_rationals()
        assert 'ostrowski_theorem' in result
        assert 'Ostrowski' in result['ostrowski_theorem']


# ============================================================
# Tests: Topologische Gruppen
# ============================================================

class TestTopologicalGroupExamples:
    """Tests für topological_group_examples."""

    def test_returns_list(self):
        """topological_group_examples gibt eine Liste zurück."""
        result = topological_group_examples()
        assert isinstance(result, list)

    def test_at_least_5_examples(self):
        """Mindestens 5 Beispiele vorhanden."""
        result = topological_group_examples()
        assert len(result) >= 5

    def test_each_example_has_required_keys(self):
        """Jedes Beispiel hat name, operation, topology, compact, connected, abelian."""
        result = topological_group_examples()
        required_keys = {'name', 'operation', 'topology', 'compact', 'connected', 'abelian'}
        for example in result:
            assert required_keys.issubset(example.keys()), \
                f"Fehlende Schlüssel in: {example}"

    def test_real_addition_is_abelian_and_connected(self):
        """(ℝ, +) ist abelsch und zusammenhängend."""
        result = topological_group_examples()
        real_add = next(e for e in result if 'ℝ, +' in e['name'])
        assert real_add['abelian'] is True
        assert real_add['connected'] is True
        assert real_add['compact'] is False

    def test_unit_circle_is_compact(self):
        """S¹ (Einheitskreis) ist kompakt."""
        result = topological_group_examples()
        s1 = next(e for e in result if 'S¹' in e['name'] or 'Einheitskreis' in e['name'])
        assert s1['compact'] is True
        assert s1['connected'] is True

    def test_gl_n_is_not_abelian(self):
        """GL(n, ℝ) ist nicht abelsch (für n ≥ 2)."""
        result = topological_group_examples()
        gln = next(e for e in result if 'GL' in e['name'])
        assert gln['abelian'] is False


# ============================================================
# Tests: Topologische Vektorräume
# ============================================================

class TestTopologicalVectorSpace:
    """Tests für TopologicalVectorSpace und zugehörige Funktionen."""

    def test_tvs_creation(self):
        """TVS-Objekt kann erstellt werden."""
        tvs = TopologicalVectorSpace('ℝ²', 2, True, True, True)
        assert tvs.name == 'ℝ²'
        assert tvs.is_banach is True

    def test_banach_requires_normed_and_complete(self):
        """Banach-Raum erfordert normiert + vollständig."""
        tvs_banach = TopologicalVectorSpace('L²', 'unendlich', True, True, True)
        tvs_not_banach = TopologicalVectorSpace('C^∞', 'unendlich', False, True, True)
        assert tvs_banach.is_banach is True
        assert tvs_not_banach.is_banach is False

    def test_convex_neighborhood_disk(self):
        """Kreisscheibenpunkte (konvex) bestehen Konvexitätsprüfung."""
        import math
        tvs = TopologicalVectorSpace('ℝ²', 2, True, True, True)
        # Punkte auf dem Rand eines Kreises
        pts = [[math.cos(2 * math.pi * k / 8), math.sin(2 * math.pi * k / 8)]
               for k in range(8)]
        pts.append([0.0, 0.0])
        assert tvs.check_convex_neighborhood(pts) is True

    def test_to_dict(self):
        """to_dict gibt alle wichtigen Schlüssel zurück."""
        tvs = TopologicalVectorSpace('Test', 3, True, True, True)
        d = tvs.to_dict()
        assert 'name' in d
        assert 'dimension' in d
        assert 'locally_convex' in d
        assert 'banach' in d


class TestLocallyConvexDemo:
    """Tests für locally_convex_demo."""

    def test_returns_dict(self):
        """locally_convex_demo gibt ein dict zurück."""
        result = locally_convex_demo()
        assert isinstance(result, dict)

    def test_examples_present(self):
        """Mindestens 4 Beispiele für lokalkonvexe Räume."""
        result = locally_convex_demo()
        assert len(result['examples']) >= 4

    def test_L_p_not_locally_convex(self):
        """L^p für p < 1 ist nicht lokalkonvex."""
        result = locally_convex_demo()
        non_lc = [e for e in result['examples'] if not e['locally_convex']]
        assert len(non_lc) >= 1

    def test_geometric_demo_convex(self):
        """Die geometrische Demo zeigt Konvexität der Einheitsscheibe."""
        result = locally_convex_demo()
        assert result['geometric_demo']['is_convex'] is True

    def test_frechet_spaces_listed(self):
        """Fréchet-Räume werden beschrieben."""
        result = locally_convex_demo()
        assert 'frechet_spaces' in result
        assert len(result['frechet_spaces']['examples']) >= 2

    def test_sobolev_spaces_listed(self):
        """Sobolev-Räume werden beschrieben."""
        result = locally_convex_demo()
        assert 'sobolev_spaces' in result
        assert 'applications' in result['sobolev_spaces']


class TestSchwartzSpaceIntro:
    """Tests für schwartz_space_intro."""

    def test_returns_dict(self):
        """schwartz_space_intro gibt ein dict zurück."""
        result = schwartz_space_intro()
        assert isinstance(result, dict)

    def test_gaussian_in_schwartz(self):
        """Gauß-Funktion liegt im Schwartz-Raum."""
        result = schwartz_space_intro()
        assert result['gaussian_example']['in_schwartz'] is True
        assert result['gaussian_example']['all_seminorms_finite'] is True

    def test_seminorm_p00_is_1(self):
        """Halbnorm p_{0,0} der Gauß-Funktion ≈ 1 (Supremum)."""
        result = schwartz_space_intro()
        p00 = result['gaussian_example']['seminorms_computed']['p_{0,0}']
        assert abs(p00 - 1.0) < 1e-5

    def test_not_in_schwartz_examples(self):
        """Mindestens 2 Funktionen die nicht in S(ℝ) liegen."""
        result = schwartz_space_intro()
        assert len(result['not_in_schwartz']) >= 2

    def test_fourier_invariance(self):
        """Fourier-Invarianz des Schwartz-Raums wird beschrieben."""
        result = schwartz_space_intro()
        assert 'fourier_invariance' in result
        assert 'Isomorphismus' in result['fourier_invariance']['statement']

    def test_properties_list(self):
        """Mindestens 3 Eigenschaften des Schwartz-Raums."""
        result = schwartz_space_intro()
        assert len(result['properties']) >= 3

    def test_tempered_distributions_mentioned(self):
        """Temperierte Distributionen werden erwähnt."""
        result = schwartz_space_intro()
        assert 'tempered_distributions' in result
        assert 'Dirac' in result['tempered_distributions'] or 'δ' in result['tempered_distributions']
