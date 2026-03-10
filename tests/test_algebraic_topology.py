"""
@file test_algebraic_topology.py
@brief Tests für die algebraische Topologie-Erweiterungen in topology.py.
@author Kurt Ingwer
@lastModified 2026-03-10

Testet simpliziale Homologie, de Rham-Kohomologie, Stokes-Satz
und den Brouwerschen Fixpunktsatz numerisch/algebraisch.
"""

import math
import pytest
from src.topology import (
    simplicial_homology,
    betti_numbers_from_adjacency,
    fundamental_group_free_generators,
    de_rham_cohomology_circle,
    de_rham_cohomology_sphere,
    stokes_theorem_verify,
    brouwer_fixed_point_evidence,
)


# ============================================================
# Hilfsfunktionen: Standard-Simplizialkomplexe
# ============================================================

def _triangle_simplices() -> dict:
    """
    Dreieck (Rand eines 2-Simplex):
    Knoten: 0, 1, 2
    Kanten: (0,1), (1,2), (0,2)
    Keine 2-Flächen (nur Rand)
    Topologie: S¹  →  β₀=1, β₁=1
    """
    return {
        0: [(0,), (1,), (2,)],
        1: [(0, 1), (1, 2), (0, 2)],
    }


def _filled_triangle_simplices() -> dict:
    """
    Ausgefülltes Dreieck (2-Simplex inklusive Innenfläche):
    Knoten: 0, 1, 2
    Kanten: (0,1), (1,2), (0,2)
    Flächen: (0,1,2)
    Topologie: D²  →  β₀=1, β₁=0, β₂=0
    """
    return {
        0: [(0,), (1,), (2,)],
        1: [(0, 1), (1, 2), (0, 2)],
        2: [(0, 1, 2)],
    }


def _torus_simplices() -> dict:
    """
    Minimale Triangulierung des Torus mit 7 Knoten (Möbius-Heawood).
    Knoten: 0..6
    Topologie: Torus T²  →  β₀=1, β₁=2, β₂=1
    Euler-Charakteristik: χ = 1 - 2 + 1 = 0

    Standardmäßige minimale Triangulierung des Torus (7 Knoten, 21 Kanten, 14 Dreiecke).
    """
    vertices = [(i,) for i in range(7)]
    # Kanten der Heawood-Triangulierung des Torus
    edges = [
        (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6),
        (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
        (2, 3), (2, 4), (2, 5), (2, 6),
        (3, 4), (3, 5), (3, 6),
        (4, 5), (4, 6),
        (5, 6),
    ]
    # 14 Dreiecke der Standardtriangulierung
    triangles = [
        (0, 1, 2), (0, 1, 3), (0, 2, 6), (0, 3, 4),
        (0, 4, 5), (0, 5, 6), (1, 2, 4), (1, 3, 5),
        (1, 4, 6), (1, 5, 6), (2, 3, 5), (2, 4, 5),
        (2, 3, 6), (3, 4, 6),
    ]
    return {
        0: vertices,
        1: edges,
        2: triangles,
    }


def _sphere_s2_simplices() -> dict:
    """
    Triangulierte S² als Oberfläche eines Tetraeders (4 Dreiecke).
    Knoten: 0, 1, 2, 3
    Kanten: alle 6 Paare
    Dreiecke: 4 Seiten des Tetraeders
    Topologie: S²  →  β₀=1, β₁=0, β₂=1
    Euler-Charakteristik: χ = 4 - 6 + 4 = 2
    """
    return {
        0: [(0,), (1,), (2,), (3,)],
        1: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
        2: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
    }


# ============================================================
# Tests: Simpliziale Homologie – Dreieck (S¹)
# ============================================================

class TestSimplicialHomologyTriangle:
    """Tests für den Rand eines Dreiecks (= Kreis S¹)."""

    def test_triangle_beta0(self):
        """Dreieck (Rand): β₀ = 1 (zusammenhängend)."""
        result = simplicial_homology(_triangle_simplices())
        assert result['h0'] == 1, f"β₀ erwartet 1, erhalten {result['h0']}"

    def test_triangle_beta1(self):
        """Dreieck (Rand): β₁ = 1 (ein Zyklus: der Rand selbst)."""
        result = simplicial_homology(_triangle_simplices())
        assert result['h1'] == 1, f"β₁ erwartet 1, erhalten {result['h1']}"

    def test_triangle_euler(self):
        """Dreieck (Rand): χ = β₀ - β₁ = 1 - 1 = 0."""
        result = simplicial_homology(_triangle_simplices())
        assert result['euler_characteristic'] == 0

    def test_filled_triangle_beta0(self):
        """Ausgefülltes Dreieck: β₀ = 1 (zusammenhängend)."""
        result = simplicial_homology(_filled_triangle_simplices())
        assert result['h0'] == 1

    def test_filled_triangle_beta1(self):
        """Ausgefülltes Dreieck: β₁ = 0 (der Zyklus wird durch die Fläche gefüllt)."""
        result = simplicial_homology(_filled_triangle_simplices())
        assert result['h1'] == 0, f"β₁ erwartet 0, erhalten {result['h1']}"

    def test_filled_triangle_euler(self):
        """Ausgefülltes Dreieck (D²): χ = 1."""
        result = simplicial_homology(_filled_triangle_simplices())
        assert result['euler_characteristic'] == 1


# ============================================================
# Tests: Simpliziale Homologie – Sphäre S²
# ============================================================

class TestSimplicialHomologySphere:
    """Tests für die 2-Sphäre als Tetraeder-Oberfläche."""

    def test_sphere_beta0(self):
        """Sphäre S²: β₀ = 1 (zusammenhängend)."""
        result = simplicial_homology(_sphere_s2_simplices())
        assert result['h0'] == 1

    def test_sphere_beta1(self):
        """Sphäre S²: β₁ = 0 (keine nicht-kontrahierbaren Schleifen)."""
        result = simplicial_homology(_sphere_s2_simplices())
        assert result['h1'] == 0, f"β₁ erwartet 0, erhalten {result['h1']}"

    def test_sphere_beta2(self):
        """Sphäre S²: β₂ = 1 (eine 2D-Hohlraum)."""
        result = simplicial_homology(_sphere_s2_simplices())
        assert result['h2'] == 1, f"β₂ erwartet 1, erhalten {result['h2']}"

    def test_sphere_euler(self):
        """Sphäre S²: χ = β₀ - β₁ + β₂ = 1 - 0 + 1 = 2."""
        result = simplicial_homology(_sphere_s2_simplices())
        assert result['euler_characteristic'] == 2, (
            f"χ erwartet 2, erhalten {result['euler_characteristic']}"
        )

    def test_sphere_betti_list_length(self):
        """Sphäre S²: Betti-Zahlen-Liste hat mindestens 3 Einträge."""
        result = simplicial_homology(_sphere_s2_simplices())
        assert len(result['betti_numbers']) >= 3


# ============================================================
# Tests: Simpliziale Homologie – Torus
# ============================================================

class TestSimplicialHomologyTorus:
    """Tests für den Torus T² (minimale Triangulierung)."""

    def test_torus_beta0(self):
        """Torus: β₀ = 1 (zusammenhängend)."""
        result = simplicial_homology(_torus_simplices())
        assert result['h0'] == 1

    def test_torus_beta1(self):
        """Torus: β₁ = 2 (zwei unabhängige Zyklen: Längen- und Breitengrad)."""
        result = simplicial_homology(_torus_simplices())
        assert result['h1'] == 2, f"β₁ erwartet 2, erhalten {result['h1']}"

    def test_torus_beta2(self):
        """Torus: β₂ = 1 (die Fundamentalklasse)."""
        result = simplicial_homology(_torus_simplices())
        assert result['h2'] == 1, f"β₂ erwartet 1, erhalten {result['h2']}"

    def test_torus_euler(self):
        """Torus: χ = β₀ - β₁ + β₂ = 1 - 2 + 1 = 0."""
        result = simplicial_homology(_torus_simplices())
        assert result['euler_characteristic'] == 0, (
            f"χ erwartet 0, erhalten {result['euler_characteristic']}"
        )


# ============================================================
# Tests: Betti-Zahlen aus Adjazenzmatrix
# ============================================================

class TestBettiNumbersFromAdjacency:
    """Tests für betti_numbers_from_adjacency."""

    def test_circle_graph_beta0(self):
        """Kreis C₃ (Dreiecksgraph): β₀ = 1."""
        # Dreieck: 3 Knoten, 3 Kanten
        adj = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        result = betti_numbers_from_adjacency(adj)
        assert result['beta_0'] == 1

    def test_circle_graph_beta1(self):
        """Kreis C₃: β₁ = 1 (ein Zyklus)."""
        adj = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        result = betti_numbers_from_adjacency(adj)
        assert result['beta_1'] == 1

    def test_circle_euler(self):
        """Kreis C₃: χ = β₀ - β₁ = 1 - 1 = 0."""
        adj = [
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ]
        result = betti_numbers_from_adjacency(adj)
        assert result['euler_characteristic'] == 0

    def test_tree_beta1_zero(self):
        """Baum (Pfadgraph P₃): β₁ = 0 (azyklisch)."""
        adj = [
            [0, 1, 0],
            [1, 0, 1],
            [0, 1, 0],
        ]
        result = betti_numbers_from_adjacency(adj)
        assert result['beta_1'] == 0

    def test_disconnected_beta0(self):
        """Zwei isolierte Knoten: β₀ = 2."""
        adj = [[0, 0], [0, 0]]
        result = betti_numbers_from_adjacency(adj)
        assert result['beta_0'] == 2

    def test_empty_graph(self):
        """Leerer Graph: alle Werte 0."""
        result = betti_numbers_from_adjacency([])
        assert result['beta_0'] == 0
        assert result['beta_1'] == 0


# ============================================================
# Tests: Fundamentalgruppe
# ============================================================

class TestFundamentalGroup:
    """Tests für fundamental_group_free_generators."""

    def test_circle_free_generators(self):
        """Kreis S¹: 1 freier Generator der Fundamentalgruppe."""
        result = fundamental_group_free_generators(_triangle_simplices())
        assert result['free_generators'] == 1

    def test_circle_not_simply_connected(self):
        """Kreis S¹: nicht einfach zusammenhängend."""
        result = fundamental_group_free_generators(_triangle_simplices())
        assert result['is_simply_connected'] is False

    def test_sphere_simply_connected(self):
        """Sphäre S²: einfach zusammenhängend (π₁(S²) = 1)."""
        result = fundamental_group_free_generators(_sphere_s2_simplices())
        assert result['is_simply_connected'] is True

    def test_sphere_free_generators_zero(self):
        """Sphäre S²: 0 freie Generatoren."""
        result = fundamental_group_free_generators(_sphere_s2_simplices())
        assert result['free_generators'] == 0

    def test_filled_triangle_simply_connected(self):
        """Ausgefülltes Dreieck D²: einfach zusammenhängend."""
        result = fundamental_group_free_generators(_filled_triangle_simplices())
        assert result['is_simply_connected'] is True

    def test_comment_is_string(self):
        """Kommentar ist ein nicht-leerer String."""
        result = fundamental_group_free_generators(_sphere_s2_simplices())
        assert isinstance(result['comment'], str)
        assert len(result['comment']) > 0


# ============================================================
# Tests: de Rham-Kohomologie – S¹
# ============================================================

class TestDeRhamCircle:
    """Tests für de Rham-Kohomologie des Kreises."""

    def test_h0_is_r(self):
        """H^0_dR(S¹) = ℝ (konstante Funktionen)."""
        result = de_rham_cohomology_circle()
        assert result['h0'] == 'ℝ'

    def test_h1_is_r(self):
        """H^1_dR(S¹) = ℝ (erzeugt von dθ)."""
        result = de_rham_cohomology_circle()
        assert result['h1'] == 'ℝ'

    def test_h0_dim(self):
        """dim(H^0(S¹)) = 1."""
        result = de_rham_cohomology_circle()
        assert result['h0_dim'] == 1

    def test_h1_dim(self):
        """dim(H^1(S¹)) = 1."""
        result = de_rham_cohomology_circle()
        assert result['h1_dim'] == 1

    def test_dtheta_not_exact(self):
        """dθ ist nicht exakt auf S¹ (Integral ≠ 0)."""
        result = de_rham_cohomology_circle()
        assert result['is_dtheta_exact'] is False

    def test_integral_dtheta_approx_2pi(self):
        """∮_{S¹} dθ ≈ 2π."""
        result = de_rham_cohomology_circle()
        assert abs(result['integral_dtheta'] - 2 * math.pi) < 0.01

    def test_euler_char_circle(self):
        """χ(S¹) = 0."""
        result = de_rham_cohomology_circle()
        assert result['euler_characteristic'] == 0

    def test_betti_numbers_circle(self):
        """Betti-Zahlen: [1, 1]."""
        result = de_rham_cohomology_circle()
        assert result['betti_numbers'] == [1, 1]


# ============================================================
# Tests: de Rham-Kohomologie – Sⁿ
# ============================================================

class TestDeRhamSphere:
    """Tests für de Rham-Kohomologie der n-Sphäre."""

    def test_s2_h0_dim(self):
        """dim(H^0(S²)) = 1."""
        result = de_rham_cohomology_sphere(n=2)
        assert result['h0_dim'] == 1

    def test_s2_h1_is_zero(self):
        """dim(H^1(S²)) = 0 (S² ist einfach zusammenhängend)."""
        result = de_rham_cohomology_sphere(n=2)
        betti = result['betti_numbers']
        assert betti[1] == 0, f"H^1(S²) erwartet 0, erhalten {betti[1]}"

    def test_s2_h2_is_r(self):
        """dim(H^2(S²)) = 1 (Orientierungsklasse)."""
        result = de_rham_cohomology_sphere(n=2)
        assert result['hn_dim'] == 1

    def test_s2_euler_char(self):
        """χ(S²) = 1 + (-1)² = 2."""
        result = de_rham_cohomology_sphere(n=2)
        assert result['euler_characteristic'] == 2

    def test_s1_euler_char(self):
        """χ(S¹) = 1 + (-1)¹ = 0."""
        result = de_rham_cohomology_sphere(n=1)
        assert result['euler_characteristic'] == 0

    def test_s3_euler_char(self):
        """χ(S³) = 1 + (-1)³ = 0."""
        result = de_rham_cohomology_sphere(n=3)
        assert result['euler_characteristic'] == 0

    def test_s2_betti_numbers(self):
        """S²: Betti-Zahlen [1, 0, 1]."""
        result = de_rham_cohomology_sphere(n=2)
        betti = result['betti_numbers']
        assert betti[0] == 1
        assert betti[1] == 0
        assert betti[2] == 1

    def test_dimension_field(self):
        """Dimensionsfeld wird korrekt gesetzt."""
        result = de_rham_cohomology_sphere(n=3)
        assert result['dimension'] == 3

    def test_invalid_dimension_raises(self):
        """Dimension n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            de_rham_cohomology_sphere(n=0)


# ============================================================
# Tests: Stokes-Satz (Grünscher Satz)
# ============================================================

class TestStokesTheorem:
    """Tests für die numerische Verifikation des Stokes-Satzes."""

    def test_stokes_verified(self):
        """Linien- und Flächenintegral stimmen überein (relativer Fehler < 1%)."""
        result = stokes_theorem_verify()
        assert result['stokes_verified'] is True, (
            f"Stokes nicht verifiziert. Rel. Fehler: {result['relative_error']:.4f}"
        )

    def test_line_integral_approx_2pi(self):
        """Linienintegral ≈ 2π ≈ 6.2832 für ω = -y dx + x dy."""
        result = stokes_theorem_verify()
        expected = 2 * math.pi
        assert abs(result['line_integral'] - expected) < 0.1, (
            f"Linienintegral erwartet ≈ {expected:.4f}, erhalten {result['line_integral']:.4f}"
        )

    def test_area_integral_approx_2pi(self):
        """Flächenintegral ≈ 2π (∫∫ 2 dA = 2·π·1²)."""
        result = stokes_theorem_verify()
        expected = 2 * math.pi
        assert abs(result['area_integral'] - expected) < 0.01, (
            f"Flächenintegral erwartet ≈ {expected:.4f}, erhalten {result['area_integral']:.4f}"
        )

    def test_relative_error_small(self):
        """Relativer Fehler < 0.01 (1%)."""
        result = stokes_theorem_verify()
        assert result['relative_error'] < 0.01

    def test_green_curl_value(self):
        """Rotor von ω = -y dx + x dy ist ∂Q/∂x - ∂P/∂y = 1-(-1) = 2."""
        result = stokes_theorem_verify()
        assert abs(result['green_curl'] - 2.0) < 1e-10


# ============================================================
# Tests: Brouwerscher Fixpunktsatz
# ============================================================

class TestBrouwerFixedPoint:
    """Tests für die numerische Evidenz des Brouwerschen Fixpunktsatzes."""

    def test_fixpoint_ratio_above_95_percent(self):
        """Mehr als 95% der Abbildungen haben einen gefundenen Fixpunkt."""
        result = brouwer_fixed_point_evidence(n_trials=500)
        assert result['found_ratio'] > 0.95, (
            f"Fixpunktquote zu niedrig: {result['found_ratio']:.3f} < 0.95"
        )

    def test_brouwer_confirmed(self):
        """Flag 'brouwer_confirmed' ist True."""
        result = brouwer_fixed_point_evidence(n_trials=300)
        assert result['brouwer_confirmed'] is True

    def test_n_trials_matches(self):
        """Anzahl der Versuche stimmt mit n_trials überein."""
        result = brouwer_fixed_point_evidence(n_trials=100)
        assert result['n_trials'] == 100

    def test_n_found_leq_n_trials(self):
        """Gefundene Fixpunkte ≤ Gesamtversuche."""
        result = brouwer_fixed_point_evidence(n_trials=200)
        assert result['n_found'] <= result['n_trials']

    def test_method_string(self):
        """'method'-Feld ist ein nicht-leerer String."""
        result = brouwer_fixed_point_evidence(n_trials=10)
        assert isinstance(result['method'], str)
        assert len(result['method']) > 0
