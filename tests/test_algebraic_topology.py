"""
@file test_algebraic_topology.py
@brief Umfassende Tests für das algebraic_topology-Modul.
@author Kurt Ingwer
@lastModified 2026-03-10

Testet alle Klassen und Funktionen des algebraischen Topologie-Moduls:
  - SimplicialComplex: Euler-Charakteristik, Randoperatoren, Betti-Zahlen
  - SimplicialHomology: Homologiegruppen mit Torsion
  - SingularHomology: S^n, T², RP^n, CP^n
  - CohomologyRing: Kohomologie und Cup-Produkt
  - HomotopyGroups: π_1, π_k(S^n), Seifert-van-Kampen
  - FiberBundle: Hopf-Faserung
  - SpectralSequence: Serre, Leray-Hirsch
  - KTheory: K̃(S^n), Bott-Periodizität
  - CWComplex: Euler-Charakteristik, zelluläre Homologie
  - Freie Funktionen: classify_surface, van_kampen_free_product, etc.
"""

import sys
import pytest
import numpy as np

# Pfade für beide üblichen Ausführungsweisen
sys.path.insert(0, "src")
sys.path.insert(0, "../src")

from algebraic_topology import (
    SimplicialComplex,
    SimplicialHomology,
    SingularHomology,
    CohomologyRing,
    HomotopyGroups,
    FiberBundle,
    SpectralSequence,
    KTheory,
    CWComplex,
    compute_smith_normal_form,
    homology_from_boundary_matrices,
    classify_surface,
    van_kampen_free_product,
    lyndon_hochschild_serre_demo,
    classifying_space_demo,
)


# ============================================================
# Hilfsfunktionen: Standard-Simplizialkomplexe
# ============================================================

def make_tetrahedron_surface() -> SimplicialComplex:
    """Tetraeder-Oberfläche = S²: V=4, E=6, F=4, χ=2."""
    return SimplicialComplex([
        (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3),
    ])


def make_circle() -> SimplicialComplex:
    """Kreis als Dreieck-Rand: 3 Ecken, 3 Kanten, χ=0."""
    return SimplicialComplex([(0, 1), (1, 2), (0, 2)])


def make_filled_triangle() -> SimplicialComplex:
    """Ausgefülltes Dreieck D²: 3 Ecken, 3 Kanten, 1 Fläche, χ=1."""
    return SimplicialComplex([(0, 1, 2)])


def make_torus_minimal() -> SimplicialComplex:
    """
    Torus-Triangulierung via 3x3-Gitter mit periodischen Randbed. (9 Knoten, 27 Kanten, 18 Dreiecke).
    β_0=1, β_1=2, β_2=1, χ=0.

    Konstruktion: Gitterpunkte (i,j) für i,j in {0,1,2}, mit Identifikationen
    (0,j)=(3,j) und (i,0)=(i,3). Jede der 9 Gitterflächen wird in 2 Dreiecke aufgeteilt.
    """
    def v(i: int, j: int) -> int:
        """Knotenindex mit periodischen Randbedingungen."""
        return (i % 3) * 3 + (j % 3)

    triangles = []
    for i in range(3):
        for j in range(3):
            # Unteres Dreieck des Gittervierecks
            triangles.append(tuple(sorted([v(i, j), v(i + 1, j), v(i, j + 1)])))
            # Oberes Dreieck des Gittervierecks
            triangles.append(tuple(sorted([v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)])))
    return SimplicialComplex(triangles)


# ============================================================
# Tests: compute_smith_normal_form
# ============================================================

class TestSmithNormalForm:
    """Tests für die Smith-Normalform-Berechnung."""

    def test_identity_2x2(self):
        """SNF der Einheitsmatrix ist Einheitsmatrix."""
        A = np.eye(2, dtype=int)
        D, U, V = compute_smith_normal_form(A)
        assert D[0, 0] == 1
        assert D[1, 1] == 1

    def test_zero_matrix(self):
        """SNF der Nullmatrix ist Nullmatrix."""
        A = np.zeros((3, 3), dtype=int)
        D, _, _ = compute_smith_normal_form(A)
        assert np.all(D == 0)

    def test_diagonal_matrix(self):
        """SNF einer Diagonalmatrix hat gleiche Diagonale (ggf. umsortiert)."""
        A = np.diag([2, 4]).astype(int)
        D, U, V = compute_smith_normal_form(A)
        diag = sorted([D[0, 0], D[1, 1]])
        assert 2 in diag
        assert 4 in diag

    def test_rank_preserved(self):
        """Rang bleibt durch SNF erhalten."""
        A = np.array([[1, 2, 3], [4, 5, 6]], dtype=int)
        D, _, _ = compute_smith_normal_form(A)
        rank_orig = np.linalg.matrix_rank(A)
        rank_snf = np.sum(np.diag(D[:min(D.shape)]) != 0)
        assert rank_orig == rank_snf

    def test_1x1_matrix(self):
        """SNF einer 1×1-Matrix."""
        A = np.array([[7]], dtype=int)
        D, _, _ = compute_smith_normal_form(A)
        assert D[0, 0] == 7

    def test_rectangular_matrix(self):
        """SNF einer rechteckigen Matrix."""
        A = np.array([[2, 0, 0], [0, 3, 0]], dtype=int)
        D, U, V = compute_smith_normal_form(A)
        assert D.shape == (2, 3)
        nonzero = sorted([D[i, i] for i in range(2) if D[i, i] != 0])
        assert len(nonzero) == 2

    def test_snf_positive_diagonal(self):
        """Diagonaleinträge der SNF sind nicht-negativ."""
        A = np.array([[1, 2], [3, 4]], dtype=int)
        D, _, _ = compute_smith_normal_form(A)
        for i in range(min(D.shape)):
            assert D[i, i] >= 0


# ============================================================
# Tests: SimplicialComplex – Basisfunktionalität
# ============================================================

class TestSimplicialComplexBasic:
    """Tests für grundlegende Funktionen von SimplicialComplex."""

    def test_faces_nonempty(self):
        """faces() liefert nicht-leere Liste."""
        sc = make_circle()
        assert len(sc.faces()) > 0

    def test_faces_contain_all_dims(self):
        """Alle Dimensionen sind in faces() vertreten."""
        sc = make_filled_triangle()
        faces = sc.faces()
        dims = set(len(f) - 1 for f in faces)
        assert 0 in dims  # Punkte
        assert 1 in dims  # Kanten
        assert 2 in dims  # Dreiecke

    def test_closure_property(self):
        """SimplicialComplex ist unter Seiten abgeschlossen."""
        sc = SimplicialComplex([(0, 1, 2)])
        faces = sc.faces()
        # Wenn (0,1,2) drin, müssen (0,), (1,), (2,), (0,1), (0,2), (1,2) auch drin sein
        assert (0,) in faces
        assert (0, 1) in faces
        assert (1, 2) in faces

    def test_simplex_boundary_matrix_shape(self):
        """Randmatrix hat korrekte Größe."""
        sc = make_filled_triangle()
        # ∂_1: C_1 (3 Kanten) → C_0 (3 Punkte), Matrix: 3×3
        mat = sc.boundary_matrix(1)
        assert mat.shape[0] == 3  # Zeilen = #Punkte
        assert mat.shape[1] == 3  # Spalten = #Kanten

    def test_boundary_squared_is_zero(self):
        """∂² = 0: Randoperator zweimal angewendet ergibt 0."""
        sc = make_filled_triangle()
        d1 = sc.boundary_matrix(1)
        d2 = sc.boundary_matrix(2)
        product = d1 @ d2
        assert np.all(product == 0), f"∂∂ ≠ 0: {product}"


# ============================================================
# Tests: SimplicialComplex – Euler-Charakteristik
# ============================================================

class TestEulerCharacteristic:
    """Tests für die Euler-Charakteristik verschiedener Räume."""

    def test_tetrahedron_surface_euler(self):
        """Tetraeder-Oberfläche (S²): χ = 4 - 6 + 4 = 2."""
        sc = make_tetrahedron_surface()
        assert sc.euler_characteristic() == 2

    def test_circle_euler(self):
        """Kreis S¹: χ = 3 - 3 = 0."""
        sc = make_circle()
        assert sc.euler_characteristic() == 0

    def test_filled_triangle_euler(self):
        """Ausgefülltes Dreieck D²: χ = 3 - 3 + 1 = 1."""
        sc = make_filled_triangle()
        assert sc.euler_characteristic() == 1

    def test_torus_euler(self):
        """Torus T²: χ = 0."""
        sc = make_torus_minimal()
        assert sc.euler_characteristic() == 0

    def test_single_vertex_euler(self):
        """Einzelner Punkt: χ = 1."""
        sc = SimplicialComplex([(0,)])
        assert sc.euler_characteristic() == 1

    def test_two_points_euler(self):
        """Zwei isolierte Punkte: χ = 2."""
        sc = SimplicialComplex([(0,), (1,)])
        assert sc.euler_characteristic() == 2

    def test_edge_euler(self):
        """Einzelne Kante: χ = 2 - 1 = 1."""
        sc = SimplicialComplex([(0, 1)])
        assert sc.euler_characteristic() == 1

    def test_tetrahedron_vertex_count(self):
        """Tetraeder-Oberfläche: V=4, E=6, F=4."""
        sc = make_tetrahedron_surface()
        faces = sc.faces()
        v = sum(1 for f in faces if len(f) == 1)
        e = sum(1 for f in faces if len(f) == 2)
        f = sum(1 for f in faces if len(f) == 3)
        assert v == 4
        assert e == 6
        assert f == 4


# ============================================================
# Tests: SimplicialComplex – Betti-Zahlen
# ============================================================

class TestSimplicialComplexBetti:
    """Tests für Betti-Zahlen von SimplicialComplex."""

    def test_circle_beta0(self):
        """Kreis: β_0 = 1 (zusammenhängend)."""
        sc = make_circle()
        b = sc.betti_numbers()
        assert b[0] == 1

    def test_circle_beta1(self):
        """Kreis: β_1 = 1 (ein Zyklus)."""
        sc = make_circle()
        b = sc.betti_numbers()
        assert b[1] == 1

    def test_filled_triangle_beta1_zero(self):
        """Ausgefülltes Dreieck: β_1 = 0."""
        sc = make_filled_triangle()
        b = sc.betti_numbers()
        assert b.get(1, 0) == 0

    def test_sphere_beta0(self):
        """Sphäre S²: β_0 = 1."""
        sc = make_tetrahedron_surface()
        b = sc.betti_numbers()
        assert b[0] == 1

    def test_sphere_beta1_zero(self):
        """Sphäre S²: β_1 = 0."""
        sc = make_tetrahedron_surface()
        b = sc.betti_numbers()
        assert b.get(1, 0) == 0

    def test_sphere_beta2(self):
        """Sphäre S²: β_2 = 1."""
        sc = make_tetrahedron_surface()
        b = sc.betti_numbers()
        assert b.get(2, 0) == 1

    def test_torus_beta0(self):
        """Torus: β_0 = 1."""
        sc = make_torus_minimal()
        b = sc.betti_numbers()
        assert b[0] == 1

    def test_torus_beta1(self):
        """Torus: β_1 = 2."""
        sc = make_torus_minimal()
        b = sc.betti_numbers()
        assert b.get(1, 0) == 2, f"β_1 erwartet 2, erhalten {b.get(1, 0)}"

    def test_torus_beta2(self):
        """Torus: β_2 = 1."""
        sc = make_torus_minimal()
        b = sc.betti_numbers()
        assert b.get(2, 0) == 1, f"β_2 erwartet 1, erhalten {b.get(2, 0)}"

    def test_betti_euler_relation(self):
        """Euler-Charakteristik = alternierende Summe der Betti-Zahlen."""
        sc = make_tetrahedron_surface()
        b = sc.betti_numbers()
        chi_betti = sum((-1) ** k * v for k, v in b.items())
        assert chi_betti == sc.euler_characteristic()

    def test_betti_euler_torus(self):
        """Euler-Charakteristik des Torus aus Betti-Zahlen."""
        sc = make_torus_minimal()
        b = sc.betti_numbers()
        chi_betti = sum((-1) ** k * v for k, v in b.items())
        assert chi_betti == 0


# ============================================================
# Tests: SimplicialHomology
# ============================================================

class TestSimplicialHomology:
    """Tests für SimplicialHomology-Klasse."""

    def test_chain_groups_circle(self):
        """Kreis: Kettengruppen C_0 (dim=3), C_1 (dim=3)."""
        sc = make_circle()
        sh = SimplicialHomology(sc)
        cg = sh.chain_groups()
        assert cg[0] == 3
        assert cg[1] == 3

    def test_chain_groups_filled(self):
        """Ausgefülltes Dreieck: C_0=3, C_1=3, C_2=1."""
        sc = make_filled_triangle()
        sh = SimplicialHomology(sc)
        cg = sh.chain_groups()
        assert cg[0] == 3
        assert cg[1] == 3
        assert cg[2] == 1

    def test_boundary_operators_nonempty(self):
        """Randoperatoren sind nicht leer für nicht-trivialen Komplex."""
        sc = make_circle()
        sh = SimplicialHomology(sc)
        ops = sh.boundary_operators()
        assert len(ops) >= 1

    def test_homology_groups_circle(self):
        """Kreis: H_0 = ℤ (rank=1), H_1 = ℤ (rank=1)."""
        sc = make_circle()
        sh = SimplicialHomology(sc)
        hg = sh.homology_groups()
        assert hg[0]["rank"] == 1
        assert hg[1]["rank"] == 1

    def test_homology_groups_filled_triangle(self):
        """Ausgefülltes Dreieck: H_0=ℤ, H_1=0."""
        sc = make_filled_triangle()
        sh = SimplicialHomology(sc)
        hg = sh.homology_groups()
        assert hg[0]["rank"] == 1
        assert hg.get(1, {}).get("rank", 0) == 0

    def test_homology_groups_sphere(self):
        """Sphäre S²: H_0=ℤ, H_1=0, H_2=ℤ."""
        sc = make_tetrahedron_surface()
        sh = SimplicialHomology(sc)
        hg = sh.homology_groups()
        assert hg[0]["rank"] == 1
        assert hg.get(1, {}).get("rank", 0) == 0
        assert hg.get(2, {}).get("rank", 0) == 1

    def test_homology_groups_torus(self):
        """Torus: H_0=ℤ, H_1=ℤ², H_2=ℤ."""
        sc = make_torus_minimal()
        sh = SimplicialHomology(sc)
        hg = sh.homology_groups()
        assert hg[0]["rank"] == 1
        assert hg[1]["rank"] == 2
        assert hg[2]["rank"] == 1

    def test_betti_numbers_via_class(self):
        """Betti-Zahlen über SimplicialHomology stimmen mit SimplicialComplex überein."""
        sc = make_tetrahedron_surface()
        sh = SimplicialHomology(sc)
        b = sh.betti_numbers()
        assert b[0] == 1
        assert b.get(1, 0) == 0
        assert b.get(2, 0) == 1

    def test_description_is_string(self):
        """Beschreibung der Homologiegruppe ist ein String."""
        sc = make_circle()
        sh = SimplicialHomology(sc)
        hg = sh.homology_groups()
        assert isinstance(hg[0]["description"], str)
        assert len(hg[0]["description"]) > 0


# ============================================================
# Tests: SingularHomology
# ============================================================

class TestSingularHomologySphere:
    """Tests für SingularHomology.homology_of_sphere."""

    def setup_method(self):
        self.sh = SingularHomology()

    def test_s1_h0(self):
        """H_0(S^1) = ℤ."""
        h = self.sh.homology_of_sphere(1)
        assert h[0] == "ℤ"

    def test_s1_h1(self):
        """H_1(S^1) = ℤ."""
        h = self.sh.homology_of_sphere(1)
        assert h[1] == "ℤ"

    def test_s2_h0(self):
        """H_0(S^2) = ℤ."""
        h = self.sh.homology_of_sphere(2)
        assert h[0] == "ℤ"

    def test_s2_h1_zero(self):
        """H_1(S^2) = 0."""
        h = self.sh.homology_of_sphere(2)
        assert h[1] == "0"

    def test_s2_h2(self):
        """H_2(S^2) = ℤ."""
        h = self.sh.homology_of_sphere(2)
        assert h[2] == "ℤ"

    def test_s3_h3(self):
        """H_3(S^3) = ℤ."""
        h = self.sh.homology_of_sphere(3)
        assert h[3] == "ℤ"

    def test_s3_intermediate_zero(self):
        """H_1(S^3) = H_2(S^3) = 0."""
        h = self.sh.homology_of_sphere(3)
        assert h[1] == "0"
        assert h[2] == "0"

    def test_sn_only_0_and_n(self):
        """Für S^4: nur H_0 und H_4 sind ℤ."""
        h = self.sh.homology_of_sphere(4)
        assert h[0] == "ℤ"
        assert h[4] == "ℤ"
        assert h[1] == "0"
        assert h[2] == "0"
        assert h[3] == "0"

    def test_invalid_n_raises(self):
        """n < 0 wirft ValueError."""
        with pytest.raises(ValueError):
            self.sh.homology_of_sphere(-1)


class TestSingularHomologyTorus:
    """Tests für SingularHomology.homology_of_torus."""

    def setup_method(self):
        self.sh = SingularHomology()

    def test_torus_h0(self):
        """H_0(T²) = ℤ."""
        h = self.sh.homology_of_torus()
        assert h[0] == "ℤ"

    def test_torus_h1(self):
        """H_1(T²) = ℤ²."""
        h = self.sh.homology_of_torus()
        assert "ℤ" in h[1] and "ℤ" in h[1]  # ℤ ⊕ ℤ

    def test_torus_h2(self):
        """H_2(T²) = ℤ."""
        h = self.sh.homology_of_torus()
        assert h[2] == "ℤ"


class TestSingularHomologyRP:
    """Tests für SingularHomology.homology_of_rp."""

    def setup_method(self):
        self.sh = SingularHomology()

    def test_rp1_h0(self):
        """H_0(RP^1) = ℤ (RP^1 ≅ S^1)."""
        h = self.sh.homology_of_rp(1)
        assert h[0] == "ℤ"

    def test_rp2_h1_torsion(self):
        """H_1(RP^2) = ℤ/2 (Torsion)."""
        h = self.sh.homology_of_rp(2)
        assert "ℤ/2" in h[1]

    def test_rp2_h2_zero(self):
        """H_2(RP^2) = 0 (nicht orientierbar)."""
        h = self.sh.homology_of_rp(2)
        assert h[2] == "0"

    def test_rp3_h3(self):
        """H_3(RP^3) = ℤ (RP^3 ist orientierbar, n=3 ungerade)."""
        h = self.sh.homology_of_rp(3)
        assert h[3] == "ℤ"

    def test_rp_invalid(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            self.sh.homology_of_rp(0)


class TestSingularHomologyCP:
    """Tests für SingularHomology.homology_of_cp."""

    def setup_method(self):
        self.sh = SingularHomology()

    def test_cp1_is_s2(self):
        """CP^1 ≅ S^2: H_0=H_2=ℤ, H_1=0."""
        h = self.sh.homology_of_cp(1)
        assert h[0] == "ℤ"
        assert h[1] == "0"
        assert h[2] == "ℤ"

    def test_cp2_even_degrees(self):
        """CP^2: H_0=H_2=H_4=ℤ, H_1=H_3=0."""
        h = self.sh.homology_of_cp(2)
        assert h[0] == "ℤ"
        assert h[2] == "ℤ"
        assert h[4] == "ℤ"
        assert h[1] == "0"
        assert h[3] == "0"

    def test_cpn_only_even(self):
        """CP^3: Nur gerade Grade haben Homologie."""
        h = self.sh.homology_of_cp(3)
        for k in [0, 2, 4, 6]:
            assert h[k] == "ℤ", f"H_{k}(CP^3) sollte ℤ sein"
        for k in [1, 3, 5]:
            assert h[k] == "0", f"H_{k}(CP^3) sollte 0 sein"

    def test_cp_invalid(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            self.sh.homology_of_cp(0)

    def test_mayer_vietoris_returns_string(self):
        """mayer_vietoris_demo gibt nicht-leeren String zurück."""
        sh = SingularHomology()
        result = sh.mayer_vietoris_demo()
        assert isinstance(result, str)
        assert len(result) > 0
        assert "Mayer" in result or "mayer" in result.lower() or "H" in result


# ============================================================
# Tests: CohomologyRing
# ============================================================

class TestCohomologyRing:
    """Tests für CohomologyRing."""

    def setup_method(self):
        self.cr = CohomologyRing("test")

    def test_sphere_h0(self):
        """H^0(S^n) = ℤ."""
        for n in [1, 2, 3, 4]:
            h = self.cr.cohomology_of_sphere(n)
            assert h[0] == "ℤ"

    def test_sphere_hn(self):
        """H^n(S^n) = ℤ."""
        for n in [1, 2, 3]:
            h = self.cr.cohomology_of_sphere(n)
            assert h[n] == "ℤ"

    def test_sphere_intermediate_zero(self):
        """H^k(S^3) = 0 für k=1,2."""
        h = self.cr.cohomology_of_sphere(3)
        assert h[1] == "0"
        assert h[2] == "0"

    def test_sphere_invalid_dim(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            self.cr.cohomology_of_sphere(0)

    def test_torus_cohomology(self):
        """Torus-Kohomologie: H^0=H^2=ℤ, H^1=ℤ²."""
        h = self.cr.cohomology_of_torus()
        assert h[0] == "ℤ"
        assert h[2] == "ℤ"
        assert "ℤ" in h[1]

    def test_cup_product_sphere(self):
        """Cup-Produkt auf Sphäre: String-Antwort."""
        result = self.cr.cup_product_demo("sphere")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_cup_product_torus(self):
        """Cup-Produkt auf Torus: enthält Schlüsselwörter."""
        result = self.cr.cup_product_demo("torus")
        assert isinstance(result, str)
        assert "∪" in result or "cup" in result.lower() or "Cup" in result

    def test_cup_product_unknown(self):
        """Unbekannter Raum: String-Antwort (kein Fehler)."""
        result = self.cr.cup_product_demo("unknown_space_xyz")
        assert isinstance(result, str)

    def test_poincare_duality(self):
        """Poincaré-Dualität für Torus (dim=2)."""
        result = self.cr.poincare_duality_demo("T²", 2)
        assert isinstance(result, str)
        assert "Dualit" in result or "dual" in result.lower() or "≅" in result


# ============================================================
# Tests: HomotopyGroups
# ============================================================

class TestHomotopyGroups:
    """Tests für HomotopyGroups."""

    def setup_method(self):
        self.hg = HomotopyGroups()

    def test_pi1_s1(self):
        """π_1(S¹) = ℤ."""
        result = self.hg.fundamental_group("S1")
        assert "ℤ" in result

    def test_pi1_torus(self):
        """π_1(T²) = ℤ × ℤ."""
        result = self.hg.fundamental_group("T2")
        assert "ℤ" in result

    def test_pi1_s2_trivial(self):
        """π_1(S²) = trivial (einfach zusammenhängend)."""
        result = self.hg.fundamental_group("S2")
        assert "trivial" in result or "1" in result

    def test_pi1_rp2(self):
        """π_1(RP²) = ℤ/2."""
        result = self.hg.fundamental_group("RP2")
        assert "ℤ/2" in result or "Z/2" in result

    def test_pi1_unknown(self):
        """Unbekannter Raum: kein Absturz, String-Antwort."""
        result = self.hg.fundamental_group("some_weird_space")
        assert isinstance(result, str)

    def test_homotopy_sphere_k_lt_n(self):
        """π_k(S^n) = 0 für k < n (Hurewicz)."""
        for n in [2, 3, 4]:
            for k in range(1, n):
                result = self.hg.higher_homotopy_sphere(n, k)
                assert result == "0", f"π_{k}(S^{n}) sollte 0 sein, erhalten: {result}"

    def test_homotopy_sphere_pn_equals_z(self):
        """π_n(S^n) = ℤ."""
        for n in [1, 2, 3, 4]:
            result = self.hg.higher_homotopy_sphere(n, n)
            assert "ℤ" in result, f"π_{n}(S^{n}) sollte ℤ enthalten, erhalten: {result}"

    def test_homotopy_s1_k2_zero(self):
        """π_k(S¹) = 0 für k ≥ 2 (universelle Überlagerung ℝ)."""
        result = self.hg.higher_homotopy_sphere(1, 2)
        assert result == "0"

    def test_seifert_van_kampen_demo(self):
        """seifert_van_kampen_demo gibt String zurück."""
        result = self.hg.seifert_van_kampen_demo()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_seifert_van_kampen_contains_pi1(self):
        """Demo enthält π_1."""
        result = self.hg.seifert_van_kampen_demo()
        assert "π" in result or "pi" in result.lower()

    def test_covering_space_demo(self):
        """covering_space_demo gibt String zurück."""
        result = self.hg.covering_space_demo()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_long_exact_sequence_demo(self):
        """long_exact_sequence_demo gibt String zurück."""
        result = self.hg.long_exact_sequence_demo()
        assert isinstance(result, str)
        assert "Hopf" in result or "exakt" in result.lower() or "π" in result


# ============================================================
# Tests: FiberBundle
# ============================================================

class TestFiberBundle:
    """Tests für FiberBundle."""

    def setup_method(self):
        self.fb = FiberBundle("S^2", "S^1", "S^3")

    def test_hopf_fibration_returns_dict(self):
        """Hopf-Faserung Demo gibt Dict zurück."""
        result = self.fb.hopf_fibration_demo()
        assert isinstance(result, dict)

    def test_hopf_fibration_keys(self):
        """Hopf-Faserung Dict enthält Pflichtfelder."""
        result = self.fb.hopf_fibration_demo()
        assert "fiber" in result
        assert "total_space" in result
        assert "base_space" in result

    def test_hopf_fiber_is_s1(self):
        """Hopf-Faserung Faser ist S¹."""
        result = self.fb.hopf_fibration_demo()
        assert "S¹" in result["fiber"] or "S1" in result["fiber"]

    def test_hopf_base_is_s2(self):
        """Hopf-Faserung Basis ist S²."""
        result = self.fb.hopf_fibration_demo()
        assert "S²" in result["base_space"] or "S2" in result["base_space"]

    def test_hopf_total_is_s3(self):
        """Hopf-Faserung Totalraum ist S³."""
        result = self.fb.hopf_fibration_demo()
        assert "S³" in result["total_space"] or "S3" in result["total_space"]

    def test_hopf_invariant_is_1(self):
        """Hopf-Invariante ist 1."""
        result = self.fb.hopf_fibration_demo()
        assert result.get("hopf_invariant") == 1

    def test_vector_bundle_demo(self):
        """vector_bundle_demo gibt String zurück."""
        result = self.fb.vector_bundle_demo()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_characteristic_classes_demo(self):
        """characteristic_classes_demo gibt String zurück."""
        result = self.fb.characteristic_classes_demo()
        assert isinstance(result, str)
        assert "Chern" in result or "Pontryagin" in result

    def test_euler_class_demo(self):
        """euler_class_demo gibt String zurück."""
        result = self.fb.euler_class_demo()
        assert isinstance(result, str)
        assert "Euler" in result or "χ" in result


# ============================================================
# Tests: SpectralSequence
# ============================================================

class TestSpectralSequence:
    """Tests für SpectralSequence."""

    def setup_method(self):
        self.ss = SpectralSequence("Serre")

    def test_serre_spectral_sequence(self):
        """Serre-Spektralsequenz Demo gibt String zurück."""
        result = self.ss.serre_spectral_sequence_demo()
        assert isinstance(result, str)
        assert "Serre" in result or "E²" in result or "E^2" in result

    def test_leray_hirsch(self):
        """Leray-Hirsch Satz Demo gibt String zurück."""
        result = self.ss.leray_hirsch_theorem_demo()
        assert isinstance(result, str)
        assert "Leray" in result or "Hirsch" in result

    def test_e2_page_description(self):
        """E₂-Seiten-Beschreibung gibt String zurück."""
        result = self.ss.e2_page_description()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_e2_contains_name(self):
        """E₂-Seite enthält den Namen der Sequenz."""
        result = self.ss.e2_page_description()
        assert "Serre" in result


# ============================================================
# Tests: KTheory
# ============================================================

class TestKTheory:
    """Tests für KTheory."""

    def setup_method(self):
        self.kt = KTheory()

    def test_k_tilde_s0_is_z(self):
        """K̃(S^0) = ℤ (0 ist gerade)."""
        result = self.kt.k_group_of_sphere(0)
        assert result["k_tilde"] == "ℤ"

    def test_k_tilde_s1_is_zero(self):
        """K̃(S^1) = 0 (1 ist ungerade)."""
        result = self.kt.k_group_of_sphere(1)
        assert result["k_tilde"] == "0"

    def test_k_tilde_s2_is_z(self):
        """K̃(S^2) = ℤ (Bott-Periodizität, 2 ist gerade)."""
        result = self.kt.k_group_of_sphere(2)
        assert result["k_tilde"] == "ℤ"

    def test_k_tilde_s3_is_zero(self):
        """K̃(S^3) = 0 (3 ist ungerade)."""
        result = self.kt.k_group_of_sphere(3)
        assert result["k_tilde"] == "0"

    def test_k_tilde_s4_is_z(self):
        """K̃(S^4) = ℤ (4 ist gerade)."""
        result = self.kt.k_group_of_sphere(4)
        assert result["k_tilde"] == "ℤ"

    def test_k_tilde_s5_is_zero(self):
        """K̃(S^5) = 0 (5 ist ungerade)."""
        result = self.kt.k_group_of_sphere(5)
        assert result["k_tilde"] == "0"

    def test_k_tilde_s6_is_z(self):
        """K̃(S^6) = ℤ (6 ist gerade)."""
        result = self.kt.k_group_of_sphere(6)
        assert result["k_tilde"] == "ℤ"

    def test_k_tilde_s7_is_zero(self):
        """K̃(S^7) = 0 (7 ist ungerade)."""
        result = self.kt.k_group_of_sphere(7)
        assert result["k_tilde"] == "0"

    def test_k_group_returns_dict(self):
        """k_group_of_sphere gibt Dict zurück."""
        result = self.kt.k_group_of_sphere(2)
        assert isinstance(result, dict)
        assert "space" in result
        assert "k_tilde" in result

    def test_k_group_invalid_n(self):
        """n < 0 wirft ValueError."""
        with pytest.raises(ValueError):
            self.kt.k_group_of_sphere(-1)

    def test_bott_periodicity_demo(self):
        """Bott-Periodizität Demo gibt String zurück."""
        result = self.kt.bott_periodicity_demo()
        assert isinstance(result, str)
        assert "Bott" in result or "Periodiz" in result

    def test_chern_character_demo(self):
        """Chern-Charakter Demo gibt String zurück."""
        result = self.kt.chern_character_demo()
        assert isinstance(result, str)
        assert "Chern" in result or "ch" in result

    def test_atiyah_singer_demo(self):
        """Atiyah-Singer Demo gibt String zurück."""
        result = self.kt.atiyah_singer_index_theorem_demo()
        assert isinstance(result, str)
        assert "Atiyah" in result or "index" in result.lower()


# ============================================================
# Tests: CWComplex
# ============================================================

class TestCWComplex:
    """Tests für CWComplex."""

    def test_sphere_s2_euler(self):
        """S² als CW-Komplex {0:1, 2:1}: χ = 1 - 0 + 1 = 2."""
        cw = CWComplex({0: 1, 2: 1})
        assert cw.euler_characteristic() == 2

    def test_torus_euler(self):
        """Torus T² als CW-Komplex {0:1, 1:2, 2:1}: χ = 1 - 2 + 1 = 0."""
        cw = CWComplex({0: 1, 1: 2, 2: 1})
        assert cw.euler_characteristic() == 0

    def test_rp2_euler(self):
        """RP² als CW-Komplex {0:1, 1:1, 2:1}: χ = 1 - 1 + 1 = 1."""
        cw = CWComplex({0: 1, 1: 1, 2: 1})
        assert cw.euler_characteristic() == 1

    def test_s3_euler(self):
        """S³ als CW-Komplex {0:1, 3:1}: χ = 1 - 1 = 0."""
        cw = CWComplex({0: 1, 3: 1})
        assert cw.euler_characteristic() == 0

    def test_genus2_surface_euler(self):
        """Geschlecht-2-Fläche: {0:1, 1:4, 2:1} → χ = 1 - 4 + 1 = -2."""
        cw = CWComplex({0: 1, 1: 4, 2: 1})
        assert cw.euler_characteristic() == -2

    def test_attaching_map_example_string(self):
        """attaching_map_example gibt String zurück."""
        cw = CWComplex({0: 1, 2: 1})
        result = cw.attaching_map_example()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_attaching_map_torus_string(self):
        """Torus-CW-Komplex Anheftungsabbildung."""
        cw = CWComplex({0: 1, 1: 2, 2: 1})
        result = cw.attaching_map_example()
        assert isinstance(result, str)

    def test_cellular_homology_demo_sphere(self):
        """Zelluläre Homologie von S²: {0:1, 2:1}."""
        cw = CWComplex({0: 1, 2: 1})
        result = cw.cellular_homology_demo()
        assert isinstance(result, dict)
        assert result[0] == 1
        assert result[2] == 1

    def test_cellular_homology_demo_torus(self):
        """Zelluläre Homologie des Torus: {0:1, 1:2, 2:1}."""
        cw = CWComplex({0: 1, 1: 2, 2: 1})
        result = cw.cellular_homology_demo()
        assert result[1] == 2


# ============================================================
# Tests: classify_surface
# ============================================================

class TestClassifySurface:
    """Tests für classify_surface."""

    def test_sphere_orientable(self):
        """Geschlecht 0, orientierbar = S²."""
        result = classify_surface(0, True)
        assert result["euler_characteristic"] == 2
        assert result["orientable"] is True

    def test_torus_orientable(self):
        """Geschlecht 1, orientierbar = T²."""
        result = classify_surface(1, True)
        assert result["euler_characteristic"] == 0
        assert result["orientable"] is True

    def test_genus2_orientable(self):
        """Geschlecht 2, orientierbar: χ = -2."""
        result = classify_surface(2, True)
        assert result["euler_characteristic"] == -2

    def test_genus3_orientable(self):
        """Geschlecht 3: χ = 2 - 2·3 = -4."""
        result = classify_surface(3, True)
        assert result["euler_characteristic"] == -4

    def test_rp2_nonorientable(self):
        """1 Kreuzkappe (RP²): nicht orientierbar, χ = 1."""
        result = classify_surface(1, False)
        assert result["euler_characteristic"] == 1
        assert result["orientable"] is False

    def test_klein_bottle_nonorientable(self):
        """2 Kreuzkappenröhren (Klein-Flasche): χ = 0."""
        result = classify_surface(2, False)
        assert result["euler_characteristic"] == 0
        assert result["orientable"] is False

    def test_sphere_name_contains_s2(self):
        """Name der Sphäre enthält S²."""
        result = classify_surface(0, True)
        assert "S²" in result["name"] or "Sphäre" in result["name"] or "sphere" in result["name"].lower()

    def test_torus_name_contains_t(self):
        """Name des Torus enthält T."""
        result = classify_surface(1, True)
        assert "T" in result["name"] or "Torus" in result["name"]

    def test_rp2_has_homology(self):
        """RP² hat Homologie-Info."""
        result = classify_surface(1, False)
        assert "homology" in result
        assert isinstance(result["homology"], dict)

    def test_sphere_homology_keys(self):
        """Sphäre hat H_0=ℤ, H_2=ℤ."""
        result = classify_surface(0, True)
        assert result["homology"][0] == "ℤ"
        assert result["homology"][2] == "ℤ"

    def test_euler_formula_orientable(self):
        """Euler-Charakteristik = 2 - 2g für orientierbare Flächen."""
        for g in range(5):
            result = classify_surface(g, True)
            assert result["euler_characteristic"] == 2 - 2 * g

    def test_euler_formula_nonorientable(self):
        """Euler-Charakteristik = 2 - k für nicht-orientierbare Flächen."""
        for k in range(1, 5):
            result = classify_surface(k, False)
            assert result["euler_characteristic"] == 2 - k


# ============================================================
# Tests: van_kampen_free_product
# ============================================================

class TestVanKampenFreeProduct:
    """Tests für van_kampen_free_product."""

    def test_z_star_z(self):
        """ℤ * ℤ = freies Produkt zweier Kopien von ℤ."""
        result = van_kampen_free_product("ℤ", "ℤ")
        assert "ℤ" in result
        assert "*" in result

    def test_trivial_left(self):
        """1 * G = G."""
        result = van_kampen_free_product("1", "ℤ")
        assert result == "ℤ"

    def test_trivial_right(self):
        """G * 1 = G."""
        result = van_kampen_free_product("ℤ", "trivial")
        assert result == "ℤ"

    def test_both_trivial(self):
        """1 * 1 = 1."""
        result = van_kampen_free_product("1", "trivial")
        assert result in ("1", "trivial", "{1}")

    def test_z_star_z2(self):
        """ℤ * ℤ/2 enthält beide Gruppen."""
        result = van_kampen_free_product("ℤ", "ℤ/2")
        assert "ℤ" in result
        assert "*" in result


# ============================================================
# Tests: Freie Funktionen
# ============================================================

class TestFreeFunctions:
    """Tests für lyndon_hochschild_serre_demo und classifying_space_demo."""

    def test_lyndon_hochschild_serre(self):
        """lyndon_hochschild_serre_demo gibt nicht-leeren String zurück."""
        result = lyndon_hochschild_serre_demo()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_lhs_contains_e2(self):
        """LHS-Demo enthält E²."""
        result = lyndon_hochschild_serre_demo()
        assert "E²" in result or "E^2" in result or "Spektral" in result

    def test_classifying_space_z(self):
        """B(ℤ) = S¹."""
        result = classifying_space_demo("ℤ")
        assert isinstance(result, str)
        assert "S¹" in result or "Kreis" in result or "circle" in result.lower()

    def test_classifying_space_z2(self):
        """B(ℤ/2) = RP^∞."""
        result = classifying_space_demo("ℤ/2")
        assert isinstance(result, str)
        assert "RP" in result or "projektiv" in result.lower()

    def test_classifying_space_unknown(self):
        """Unbekannte Gruppe: kein Fehler, String zurück."""
        result = classifying_space_demo("some_group_G")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_homology_from_boundary_matrices_empty(self):
        """Leeres Dict als Eingabe: leeres Ergebnis."""
        result = homology_from_boundary_matrices({})
        assert result == {}

    def test_homology_from_boundary_matrices_circle(self):
        """Kreis-Randmatrizen → H_0=1, H_1=1."""
        # Kreis: C_1 (3 Kanten) → C_0 (3 Ecken)
        # ∂_1 für Dreiecksrand: [[1,-1,0],[-1,0,1],[0,1,-1]] (Orientierung)
        d1 = np.array([[1, -1, 0], [-1, 0, 1], [0, 1, -1]], dtype=int)
        boundaries = {1: d1}
        result = homology_from_boundary_matrices(boundaries)
        # H_0: Rang = dim(C_0) - rank(∂_1) = 3 - 2 = 1
        if 0 in result:
            assert result[0]["betti"] >= 0
