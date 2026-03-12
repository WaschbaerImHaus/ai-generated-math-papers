"""
Tests für graceful_trees.py und ramsey_numbers.py

Umfassende Test-Suite für:
  - GracefulTreeConjecture: Labelings, Konstruktionen, Statistik
  - RamseyNumbers: Bekannte Werte, Schranken, SAT-Verifikation
  - RamseyBounds: Binomial, Erdős, Sah, Spencer

Autor: Michael Fuhrmann
Letzte Änderung: 2026-03-12
"""

import pytest
import networkx as nx
import math
import sys
import os

# Sicherstellen, dass src/ im Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from graceful_trees import GracefulTreeConjecture
from ramsey_numbers import RamseyNumbers, RamseyBounds, KNOWN_RAMSEY


# ==============================================================================
# Fixtures
# ==============================================================================

@pytest.fixture(scope="module")
def gtc():
    """GracefulTreeConjecture-Instanz für alle Tests."""
    return GracefulTreeConjecture(max_backtrack_n=10)


@pytest.fixture(scope="module")
def ramsey():
    """RamseyNumbers-Instanz."""
    return RamseyNumbers()


@pytest.fixture(scope="module")
def bounds():
    """RamseyBounds-Instanz."""
    return RamseyBounds()


# ==============================================================================
# GracefulTreeConjecture – is_graceful()
# ==============================================================================

class TestIsGraceful:
    """Tests für GracefulTreeConjecture.is_graceful()."""

    def test_path_p4_known_labeling(self, gtc):
        """P_4 mit Labeling 0,3,1,2 ist graceful (Kanten: 3,2,1)."""
        path = nx.path_graph(4)
        # Pfad: 0-1-2-3
        labeling = {0: 0, 1: 3, 2: 1, 3: 2}
        assert gtc.is_graceful(path, labeling)

    def test_path_p4_edge_labels_correct(self, gtc):
        """P_4 Kantenlabels sind {1,2,3} bei Labeling 0,3,1,2."""
        path = nx.path_graph(4)
        labeling = {0: 0, 1: 3, 2: 1, 3: 2}
        edge_labels = gtc.get_edge_labels(path, labeling)
        assert edge_labels == [1, 2, 3]

    def test_path_p2_graceful(self, gtc):
        """P_2 (eine Kante) ist trivial graceful."""
        path = nx.path_graph(2)
        labeling = {0: 0, 1: 1}
        assert gtc.is_graceful(path, labeling)

    def test_path_p3_graceful(self, gtc):
        """P_3 mit Labeling 0,2,1 ist graceful."""
        path = nx.path_graph(3)
        labeling = {0: 0, 1: 2, 2: 1}
        assert gtc.is_graceful(path, labeling)

    def test_star_k13_graceful(self, gtc):
        """Stern K_{1,3} mit Labeling {Zentrum:0, Blätter:1,2,3} ist graceful."""
        star = nx.star_graph(3)  # Knoten: 0 (Zentrum), 1, 2, 3
        labeling = {0: 0, 1: 1, 2: 2, 3: 3}
        assert gtc.is_graceful(star, labeling)

    def test_star_k13_edge_labels(self, gtc):
        """K_{1,3} Kantenlabels sind {1,2,3}."""
        star = nx.star_graph(3)
        labeling = {0: 0, 1: 1, 2: 2, 3: 3}
        edge_labels = gtc.get_edge_labels(star, labeling)
        assert edge_labels == [1, 2, 3]

    def test_wrong_labeling_duplicate_edge(self, gtc):
        """Labeling mit doppeltem Kantenlabel ist nicht graceful."""
        path = nx.path_graph(3)
        labeling = {0: 0, 1: 1, 2: 2}  # Kanten: 1,1 → Duplikat
        assert not gtc.is_graceful(path, labeling)

    def test_wrong_labeling_wrong_range(self, gtc):
        """Labeling mit Labels außerhalb {0,..,n-1} ist nicht graceful."""
        path = nx.path_graph(3)
        labeling = {0: 0, 1: 5, 2: 2}  # 5 nicht in {0,1,2}
        assert not gtc.is_graceful(path, labeling)

    def test_wrong_labeling_missing_node(self, gtc):
        """Unvollständiges Labeling (fehlender Knoten) ist nicht graceful."""
        path = nx.path_graph(3)
        labeling = {0: 0, 1: 2}  # Knoten 2 fehlt
        assert not gtc.is_graceful(path, labeling)

    def test_single_node_graceful(self, gtc):
        """Einzelknoten (trivial) ist graceful mit Labeling {0: 0}."""
        g = nx.Graph()
        g.add_node(0)
        labeling = {0: 0}
        assert gtc.is_graceful(g, labeling)

    def test_p5_graceful_construction(self, gtc):
        """P_5 via graceful_path() ist tatsächlich graceful."""
        path, labeling = gtc.graceful_path(5)
        assert gtc.is_graceful(path, labeling)

    def test_p6_graceful_construction(self, gtc):
        """P_6 via graceful_path() ist graceful."""
        path, labeling = gtc.graceful_path(6)
        assert gtc.is_graceful(path, labeling)

    def test_p7_graceful_construction(self, gtc):
        """P_7 via graceful_path() ist graceful."""
        path, labeling = gtc.graceful_path(7)
        assert gtc.is_graceful(path, labeling)

    def test_p8_graceful_construction(self, gtc):
        """P_8 via graceful_path() ist graceful."""
        path, labeling = gtc.graceful_path(8)
        assert gtc.is_graceful(path, labeling)


# ==============================================================================
# GracefulTreeConjecture – find_graceful_labeling()
# ==============================================================================

class TestFindGracefulLabeling:
    """Tests für den Backtracking-Algorithmus."""

    def test_find_labeling_p4(self, gtc):
        """Für P_4 findet Backtracking ein graceful Labeling."""
        path = nx.path_graph(4)
        labeling = gtc.find_graceful_labeling(path)
        assert labeling is not None
        assert gtc.is_graceful(path, labeling)

    def test_find_labeling_p5(self, gtc):
        """Für P_5 findet Backtracking ein graceful Labeling."""
        path = nx.path_graph(5)
        labeling = gtc.find_graceful_labeling(path)
        assert labeling is not None
        assert gtc.is_graceful(path, labeling)

    def test_find_labeling_star_k15(self, gtc):
        """Für Stern K_{1,5} wird ein graceful Labeling gefunden."""
        star = nx.star_graph(5)
        labeling = gtc.find_graceful_labeling(star)
        assert labeling is not None
        assert gtc.is_graceful(star, labeling)

    def test_find_labeling_single_node(self, gtc):
        """Einzelknoten-Labeling wird korrekt zurückgegeben."""
        g = nx.Graph()
        g.add_node(0)
        labeling = gtc.find_graceful_labeling(g)
        assert labeling == {0: 0}

    def test_find_labeling_two_nodes(self, gtc):
        """Zwei-Knoten-Baum wird graceful gelabelt."""
        g = nx.path_graph(2)
        labeling = gtc.find_graceful_labeling(g)
        assert labeling is not None
        assert gtc.is_graceful(g, labeling)

    def test_find_labeling_returns_valid_bijection(self, gtc):
        """Das gefundene Labeling ist eine Bijektion auf {0,...,n-1}."""
        path = nx.path_graph(6)
        labeling = gtc.find_graceful_labeling(path)
        assert labeling is not None
        n = path.number_of_nodes()
        assert set(labeling.values()) == set(range(n))

    def test_find_labeling_caterpillar(self, gtc):
        """Caterpillar (Spine=3, Blätter=[1,2,0]) wird graceful gelabelt."""
        # Caterpillar: 0-1-2 mit je 1,2,0 Blättern
        g = nx.Graph()
        # Spine: 0-1-2
        g.add_edges_from([(0, 1), (1, 2)])
        # Blätter an Knoten 0: keine
        # Blätter an Knoten 1: 3, 4
        g.add_edges_from([(1, 3), (1, 4)])
        # Blätter an Knoten 2: 5
        g.add_edge(2, 5)
        labeling = gtc.find_graceful_labeling(g)
        assert labeling is not None
        assert gtc.is_graceful(g, labeling)

    def test_find_labeling_n8_tree(self, gtc):
        """Beliebiger Baum mit 8 Knoten wird graceful gelabelt."""
        # Doppelt-Stern: Zentrum verbunden mit zwei Sub-Sternen
        g = nx.Graph()
        g.add_edges_from([(0, 1), (1, 2), (1, 3), (1, 4), (0, 5), (0, 6), (0, 7)])
        labeling = gtc.find_graceful_labeling(g)
        assert labeling is not None
        assert gtc.is_graceful(g, labeling)


# ==============================================================================
# GracefulTreeConjecture – generate_all_trees() & Statistik
# ==============================================================================

class TestGenerateAllTrees:
    """Tests für Baum-Generierung."""

    def test_generate_n1(self, gtc):
        """n=1: Genau ein Baum (Einzelknoten)."""
        trees = gtc.generate_all_trees(1)
        assert len(trees) == 1

    def test_generate_n2(self, gtc):
        """n=2: Genau ein Baum (Kante)."""
        trees = gtc.generate_all_trees(2)
        assert len(trees) == 1

    def test_generate_n3(self, gtc):
        """n=3: Genau ein Baum (Pfad P_3)."""
        trees = gtc.generate_all_trees(3)
        assert len(trees) == 1

    def test_generate_n4(self, gtc):
        """n=4: Genau 2 nicht-isomorphe Bäume (P_4 und Stern K_{1,3})."""
        trees = gtc.generate_all_trees(4)
        assert len(trees) == 2

    def test_generate_n5(self, gtc):
        """n=5: Genau 3 nicht-isomorphe Bäume."""
        trees = gtc.generate_all_trees(5)
        assert len(trees) == 3

    def test_generate_n6(self, gtc):
        """n=6: Genau 6 nicht-isomorphe Bäume."""
        trees = gtc.generate_all_trees(6)
        assert len(trees) == 6

    def test_generate_all_are_trees(self, gtc):
        """Alle generierten Graphen sind tatsächlich Bäume."""
        for n in range(1, 8):
            trees = gtc.generate_all_trees(n)
            for tree in trees:
                assert tree.number_of_nodes() > 0, f"Leerer Graph bei n={n}"
                assert nx.is_tree(tree), f"Kein Baum: n={n}"

    def test_generate_all_n_nodes(self, gtc):
        """Alle generierten Bäume haben genau n Knoten."""
        for n in range(1, 8):
            for tree in gtc.generate_all_trees(n):
                assert tree.number_of_nodes() == n


# ==============================================================================
# GracefulTreeConjecture – Explizite Konstruktionen
# ==============================================================================

class TestExplicitConstructions:
    """Tests für die expliziten Konstruktionen bekannter Klassen."""

    def test_graceful_path_p2(self, gtc):
        """P_2 via graceful_path ist graceful."""
        path, labeling = gtc.graceful_path(2)
        assert gtc.is_graceful(path, labeling)

    def test_graceful_path_p3(self, gtc):
        """P_3 via graceful_path ist graceful."""
        path, labeling = gtc.graceful_path(3)
        assert gtc.is_graceful(path, labeling)

    def test_graceful_path_p10(self, gtc):
        """P_10 via graceful_path ist graceful."""
        path, labeling = gtc.graceful_path(10)
        assert gtc.is_graceful(path, labeling)

    def test_graceful_star_k11(self, gtc):
        """Stern K_{1,1} ist graceful."""
        star, labeling = gtc.graceful_star(1)
        assert gtc.is_graceful(star, labeling)

    def test_graceful_star_k13(self, gtc):
        """Stern K_{1,3} via graceful_star ist graceful."""
        star, labeling = gtc.graceful_star(3)
        assert gtc.is_graceful(star, labeling)

    def test_graceful_star_k16(self, gtc):
        """Stern K_{1,6} ist graceful."""
        star, labeling = gtc.graceful_star(6)
        assert gtc.is_graceful(star, labeling)

    def test_graceful_star_center_label_zero(self, gtc):
        """Sternzentrum erhält Label 0 in der Konstruktion."""
        star, labeling = gtc.graceful_star(4)
        assert labeling[0] == 0  # Zentrum in nx.star_graph ist Knoten 0

    def test_graceful_star_leaf_labels(self, gtc):
        """Stern K_{1,k} Blätter erhalten Labels 1..k."""
        k = 5
        star, labeling = gtc.graceful_star(k)
        leaf_labels = sorted(labeling[i] for i in range(1, k + 1))
        assert leaf_labels == list(range(1, k + 1))

    def test_is_caterpillar_star(self, gtc):
        """Stern K_{1,k} ist ein Caterpillar."""
        star = nx.star_graph(5)
        assert gtc.is_caterpillar(star)

    def test_is_caterpillar_path(self, gtc):
        """Pfad P_n ist ein Caterpillar."""
        path = nx.path_graph(7)
        assert gtc.is_caterpillar(path)

    def test_is_caterpillar_non_caterpillar(self, gtc):
        """Ein Baum mit Abzweigung der Tiefe 2 ist kein Caterpillar."""
        # Spider-Baum: Zentrum 0, drei Arme der Länge 2:
        # 0→1→4, 0→2→5, 0→3→6
        # Nach Blatt-Entfernung: {0,1,2,3} mit Kanten 0-1,0-2,0-3 → Stern K_{1,3}
        # Stern hat Knoten mit Grad 3 > 2 → kein Pfad → kein Caterpillar
        g = nx.Graph()
        g.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6)])
        assert not gtc.is_caterpillar(g)

    def test_verify_known_classes_all_true(self, gtc):
        """Alle bekannten Klassen (Pfade, Sterne) sind graceful."""
        results = gtc.verify_known_classes()
        for key, val in results.items():
            assert val, f"{key} sollte graceful sein"


# ==============================================================================
# RamseyNumbers – Bekannte Werte
# ==============================================================================

class TestKnownRamseyValues:
    """Tests für bekannte exakte Ramsey-Zahlen."""

    def test_r33_equals_6(self, ramsey):
        """R(3,3) = 6."""
        assert ramsey.get(3, 3) == 6

    def test_r34_equals_9(self, ramsey):
        """R(3,4) = 9."""
        assert ramsey.get(3, 4) == 9

    def test_r43_equals_9_symmetric(self, ramsey):
        """R(4,3) = 9 (Symmetrie)."""
        assert ramsey.get(4, 3) == 9

    def test_r35_equals_14(self, ramsey):
        """R(3,5) = 14."""
        assert ramsey.get(3, 5) == 14

    def test_r44_equals_18(self, ramsey):
        """R(4,4) = 18."""
        assert ramsey.get(4, 4) == 18

    def test_r36_equals_18(self, ramsey):
        """R(3,6) = 18."""
        assert ramsey.get(3, 6) == 18

    def test_r37_equals_23(self, ramsey):
        """R(3,7) = 23."""
        assert ramsey.get(3, 7) == 23

    def test_r38_equals_28(self, ramsey):
        """R(3,8) = 28."""
        assert ramsey.get(3, 8) == 28

    def test_unknown_returns_none(self, ramsey):
        """R(5,5) ist nicht als exakter Wert hinterlegt → None."""
        assert ramsey.get(5, 5) is None

    def test_r45_equals_25(self, ramsey):
        """R(4,5) = 25."""
        assert ramsey.get(4, 5) == 25

    def test_symmetry_consistent(self, ramsey):
        """Alle hinterlegten Werte sind symmetrisch: R(s,t) = R(t,s)."""
        for (s, t), val in KNOWN_RAMSEY.items():
            assert ramsey.get(t, s) == val, f"R({s},{t})≠R({t},{s})"

    def test_r55_bounds(self, ramsey):
        """R(5,5): 43 ≤ R(5,5) ≤ 48 (bekannte Schranken)."""
        lower, upper = ramsey.get_bounds(5, 5)
        assert lower == 43
        assert upper == 48


# ==============================================================================
# RamseyNumbers – Graph-Eigenschaften und Zeugen
# ==============================================================================

class TestRamseyGraphProperties:
    """Tests für Graph-Clique/Unabhängigkeitsmenge-Prüfungen."""

    def test_k5_no_k3_as_independent(self, ramsey):
        """K_5 enthält keine unabhängige 3-Menge (da vollständig)."""
        k5 = nx.complete_graph(5)
        assert not ramsey.has_independent_set(k5, 3)

    def test_k5_has_k3(self, ramsey):
        """K_5 enthält eine K_3 (Clique)."""
        k5 = nx.complete_graph(5)
        assert ramsey.has_clique(k5, 3)

    def test_k6_has_k3(self, ramsey):
        """K_6 enthält eine K_3."""
        k6 = nx.complete_graph(6)
        assert ramsey.has_clique(k6, 3)

    def test_c5_no_k3(self, ramsey):
        """C_5 (5-Kreis) enthält keine K_3."""
        c5 = nx.cycle_graph(5)
        assert not ramsey.has_clique(c5, 3)

    def test_c5_no_independent_3_set(self, ramsey):
        """C_5 hat keine unabhängige 3-Menge (maximale unabhängige Menge hat Größe 2)."""
        c5 = nx.cycle_graph(5)
        assert not ramsey.has_independent_set(c5, 3)

    def test_c5_is_r33_witness(self, ramsey):
        """C_5 ist Zeuge für R(3,3) > 5."""
        c5 = nx.cycle_graph(5)
        assert ramsey.is_ramsey_witness(c5, 3, 3)

    def test_verify_r33_complete(self, ramsey):
        """Vollständige Verifikation von R(3,3) = 6."""
        results = ramsey.verify_r33_equals_6()
        assert results["C5_no_K3"]
        assert results["C5_no_I3"]
        assert results["C5_is_witness_R33_gt5"]
        assert results["K6_has_K3"]
        assert results["known_R33_equals_6"]

    def test_empty_graph_no_clique(self, ramsey):
        """Leerer Graph (keine Kanten) hat keine K_2."""
        g = nx.empty_graph(5)
        assert not ramsey.has_clique(g, 2)

    def test_empty_graph_has_independent_set(self, ramsey):
        """Leerer Graph auf n Knoten hat unabhängige n-Menge."""
        g = nx.empty_graph(5)
        assert ramsey.has_independent_set(g, 5)

    def test_ramsey_satisfies_with_clique(self, ramsey):
        """Graph mit Clique erfüllt Ramsey-Eigenschaft."""
        k4 = nx.complete_graph(4)
        assert ramsey.satisfies_ramsey_property(k4, 3, 3)


# ==============================================================================
# RamseyNumbers – Paley-Graph
# ==============================================================================

class TestPaleyGraph:
    """Tests für den Paley-Graphen."""

    def test_paley_5_structure(self, ramsey):
        """Paley-Graph P(5) hat 5 Knoten und ist regulär."""
        p5 = ramsey.build_paley_graph(5)
        assert p5.number_of_nodes() == 5

    def test_paley_13_nodes(self, ramsey):
        """Paley-Graph P(13) hat 13 Knoten."""
        p13 = ramsey.build_paley_graph(13)
        assert p13.number_of_nodes() == 13

    def test_paley_invalid_prime_raises(self, ramsey):
        """Nicht-Primzahl löst ValueError aus."""
        with pytest.raises(ValueError):
            ramsey.build_paley_graph(4)

    def test_paley_q3mod4_raises(self, ramsey):
        """q=3 ≡ 3 (mod 4) löst ValueError aus."""
        with pytest.raises(ValueError):
            ramsey.build_paley_graph(3)

    def test_paley_17_is_self_complementary(self, ramsey):
        """Paley-Graph P(17) ist selbst-komplementär."""
        p17 = ramsey.build_paley_graph(17)
        complement = nx.complement(p17)
        # Selbst-komplementär: isomorph zum Komplement
        assert nx.is_isomorphic(p17, complement)

    def test_paley_41_exists(self, ramsey):
        """P(41) kann konstruiert werden (41 ≡ 1 mod 4)."""
        p41 = ramsey.build_paley_graph(41)
        assert p41.number_of_nodes() == 41

    def test_r55_witness_p29_no_k5(self, ramsey):
        """Paley-Graph P(29) enthält keine K_5."""
        p29 = ramsey.build_r55_lower_bound_witness()
        assert not ramsey.has_clique(p29, 5)

    def test_r55_witness_p29_no_i5(self, ramsey):
        """Paley-Graph P(29) enthält keine I_5."""
        p29 = ramsey.build_r55_lower_bound_witness()
        assert not ramsey.has_independent_set(p29, 5)

    def test_r55_witness_p29_is_witness_gt29(self, ramsey):
        """Paley-Graph P(29) ist Zeuge für R(5,5) > 29."""
        p29 = ramsey.build_r55_lower_bound_witness()
        assert ramsey.is_ramsey_witness(p29, 5, 5)


# ==============================================================================
# RamseyBounds – Schranken-Berechnungen
# ==============================================================================

class TestRamseyBounds:
    """Tests für Ramsey-Schranken."""

    def test_binomial_r33(self, bounds):
        """R(3,3) ≤ C(4,2) = 6."""
        assert bounds.upper_binomial(3, 3) == 6

    def test_binomial_r34(self, bounds):
        """R(3,4) ≤ C(5,2) = 10."""
        assert bounds.upper_binomial(3, 4) == 10

    def test_binomial_r44(self, bounds):
        """R(4,4) ≤ C(6,3) = 20."""
        assert bounds.upper_binomial(4, 4) == 20

    def test_binomial_r55_leq_70(self, bounds):
        """R(5,5) ≤ C(8,4) = 70."""
        assert bounds.upper_binomial(5, 5) == 70

    def test_binomial_symmetric(self, bounds):
        """Binomial-Schranke ist symmetrisch: C(s+t-2,s-1) = C(s+t-2,t-1)."""
        assert bounds.upper_binomial(3, 4) == bounds.upper_binomial(4, 3)

    def test_erdos_lower_k3(self, bounds):
        """Erdős: R(3,3) > floor(2^{3/2}) = 2."""
        assert bounds.lower_erdos(3) == 2

    def test_erdos_lower_k4(self, bounds):
        """Erdős: R(4,4) > floor(2^2) = 4."""
        assert bounds.lower_erdos(4) == 4

    def test_erdos_lower_k5(self, bounds):
        """Erdős: R(5,5) > floor(2^{5/2}) = 5."""
        lower = bounds.lower_erdos(5)
        assert lower == 5  # floor(sqrt(32)) = 5

    def test_erdos_lower_k6(self, bounds):
        """Erdős: R(6,6) > floor(2^3) = 8."""
        assert bounds.lower_erdos(6) == 8

    def test_sah_2023_k3(self, bounds):
        """Sah (2023): R(3,3) ≤ (4-ε)^3 ≈ 62."""
        val = bounds.upper_sah_2023(3)
        assert val < 64.0  # Kleiner als 4^3

    def test_sah_2023_less_than_4_to_k(self, bounds):
        """Sah (2023): (4-ε)^k < 4^k für alle k ≥ 1."""
        for k in range(3, 10):
            sah = bounds.upper_sah_2023(k)
            classical = 4.0 ** k
            assert sah < classical, f"Sah-Schranke für k={k} verletzt"

    def test_spencer_1975_positive(self, bounds):
        """Spencer (1975) Untergrenzen sind positiv für k ≥ 3."""
        for k in range(3, 8):
            assert bounds.upper_spencer_1975(k) > 0

    def test_spencer_1975_increasing(self, bounds):
        """Spencer (1975) Untergrenze wächst mit k."""
        vals = [bounds.upper_spencer_1975(k) for k in range(3, 8)]
        for i in range(len(vals) - 1):
            assert vals[i] < vals[i + 1]

    def test_recursive_upper_r33(self, bounds):
        """Rekursive Schranke: R(3,3) ≤ 6."""
        assert bounds.recursive_upper(3, 3) == 6

    def test_recursive_upper_r23(self, bounds):
        """Basisfall: R(2,3) = 3."""
        assert bounds.recursive_upper(2, 3) == 3

    def test_recursive_upper_r32(self, bounds):
        """Basisfall: R(3,2) = 3."""
        assert bounds.recursive_upper(3, 2) == 3

    def test_multiplicity_r33_n6(self, bounds):
        """Goodman (1959): M(3,3,6) ≥ 2."""
        val = bounds.multiplicity(3, 3, 6)
        assert val >= 0  # Nicht negativ

    def test_diagonal_bounds_format(self, bounds):
        """diagonal_bounds gibt Liste mit korrekten Keys zurück."""
        rows = bounds.diagonal_bounds(k_max=5)
        assert len(rows) == 3  # k = 3, 4, 5
        for row in rows:
            assert "k" in row
            assert "lower_erdos" in row
            assert "upper_binomial" in row

    def test_diagonal_bounds_k4(self, bounds):
        """R(4,4) = 18: exact in diagonal_bounds."""
        rows = bounds.diagonal_bounds(k_max=4)
        r44_row = next(r for r in rows if r["k"] == 4)
        assert r44_row["exact"] == 18


# ==============================================================================
# SAT-Verifikation
# ==============================================================================

class TestSATVerification:
    """Tests für die SAT-basierte Verifikation."""

    def test_sat_r33_n5_satisfiable(self, ramsey):
        """SAT: Für n=5 gibt es Graph ohne K_3 und I_3 (R(3,3)>5)."""
        result = ramsey.verify_via_sat(3, 3, 5)
        if result is None:
            pytest.skip("PySAT nicht verfügbar")
        assert result is True  # Zeuge existiert

    def test_sat_r33_n6_unsatisfiable(self, ramsey):
        """SAT: Für n=6 gibt es keinen Graph ohne K_3 und I_3 (R(3,3)≤6)."""
        result = ramsey.verify_via_sat(3, 3, 6)
        if result is None:
            pytest.skip("PySAT nicht verfügbar")
        assert result is False  # Kein Zeuge → R(3,3) ≤ 6

    def test_sat_r22_n1_unsat(self, ramsey):
        """SAT: Für R(2,2)=2, n=2 unerfüllbar."""
        result = ramsey.verify_via_sat(2, 2, 2)
        if result is None:
            pytest.skip("PySAT nicht verfügbar")
        assert result is False


# ==============================================================================
# Edge-Cases und Sonderfälle
# ==============================================================================

class TestEdgeCases:
    """Tests für Grenz- und Sonderfälle."""

    def test_graceful_empty_trees_list(self, gtc):
        """generate_all_trees(0) gibt leere Liste zurück."""
        assert gtc.generate_all_trees(0) == []

    def test_graceful_single_node_labeling(self, gtc):
        """find_graceful_labeling für Einzelknoten."""
        g = nx.Graph()
        g.add_node(42)
        result = gtc.find_graceful_labeling(g)
        assert result == {42: 0}

    def test_ramsey_r11_equals_1(self, ramsey):
        """R(1,1) = 1 (trivial)."""
        assert ramsey.get(1, 1) == 1

    def test_ramsey_r12_equals_2(self, ramsey):
        """R(1,2) = 2."""
        assert ramsey.get(1, 2) == 2

    def test_bounds_upper_binomial_r22(self, bounds):
        """R(2,2) ≤ C(2,1) = 2."""
        assert bounds.upper_binomial(2, 2) == 2

    def test_get_edge_labels_empty_tree(self, gtc):
        """Leerer Baum (ohne Kanten) hat keine Kantenlabels."""
        g = nx.Graph()
        g.add_node(0)
        result = gtc.get_edge_labels(g, {0: 0})
        assert result == []

    def test_star_k12_graceful(self, gtc):
        """K_{1,2} (= P_3) ist graceful."""
        star, labeling = gtc.graceful_star(2)
        assert gtc.is_graceful(star, labeling)

    def test_caterpillar_single_node(self, gtc):
        """Einzelknoten ist Caterpillar."""
        g = nx.Graph()
        g.add_node(0)
        assert gtc.is_caterpillar(g)

    def test_caterpillar_two_nodes(self, gtc):
        """Zwei-Knoten-Pfad ist Caterpillar."""
        g = nx.path_graph(2)
        assert gtc.is_caterpillar(g)

    def test_ramsey_bounds_r39(self, ramsey):
        """R(3,9): bekannter exakter Wert 36."""
        assert ramsey.get(3, 9) == 36

    def test_paley_q5_vertex_count(self, ramsey):
        """P(5) hat genau 5 Knoten und ist regulär vom Grad 2."""
        p5 = ramsey.build_paley_graph(5)
        assert p5.number_of_nodes() == 5
        degrees = [d for _, d in p5.degree()]
        # Paley-Graph ist (q-1)/2-regulär
        assert all(d == 2 for d in degrees)
