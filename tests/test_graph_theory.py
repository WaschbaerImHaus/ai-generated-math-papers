"""
@file test_graph_theory.py
@brief Tests für das Graphentheorie-Modul (graph_theory.py).
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
from src.graph_theory import (
    # Graph-Klasse
    Graph,
    # Traversierung
    bfs,
    dfs,
    dijkstra,
    bellman_ford,
    floyd_warshall,
    # Spanning Trees
    kruskal_mst,
    prim_mst,
    topological_sort,
    # Graphenanalyse
    chromatic_number_greedy,
    graph_coloring_greedy,
    is_eulerian,
    is_hamiltonian_path,
    clique_number,
    independence_number,
    graph_density,
    # Graphenkonstruktoren
    complete_graph,
    cycle_graph,
    path_graph,
    bipartite_complete_graph,
    petersen_graph,
    grid_graph,
    # Kombinatorik
    binomial_coefficient,
    stirling_numbers_second_kind,
    bell_numbers,
    catalan_number,
    derangements,
    partition_count,
    multinomial_coefficient,
)


# ============================================================
# Tests: Graph-Klasse – Grundoperationen
# ============================================================

class TestGraphBasicOperations:
    """Tests für Grundoperationen der Graph-Klasse."""

    def test_add_vertex(self):
        """Knoten hinzufügen."""
        g = Graph()
        g.add_vertex(1)
        g.add_vertex(2)
        assert 1 in g.vertices()
        assert 2 in g.vertices()

    def test_add_edge_undirected(self):
        """Kante in ungerichtetem Graphen (beide Richtungen)."""
        g = Graph(directed=False)
        g.add_edge(1, 2)
        assert 2 in g.adj[1]
        assert 1 in g.adj[2]  # Auch Rückkante

    def test_add_edge_directed(self):
        """Kante in gerichtetem Graphen (nur eine Richtung)."""
        g = Graph(directed=True)
        g.add_edge(1, 2)
        assert 2 in g.adj[1]
        assert 1 not in g.adj.get(2, set())

    def test_remove_vertex(self):
        """Knoten und zugehörige Kanten entfernen."""
        g = Graph()
        g.add_edge(1, 2)
        g.add_edge(1, 3)
        g.remove_vertex(1)
        assert 1 not in g.vertices()
        assert 1 not in g.adj.get(2, set())
        assert 1 not in g.adj.get(3, set())

    def test_remove_edge(self):
        """Kante entfernen."""
        g = Graph()
        g.add_edge(1, 2)
        g.remove_edge(1, 2)
        assert 2 not in g.adj[1]
        assert 1 not in g.adj[2]

    def test_vertices_returns_list(self):
        """vertices() gibt Liste zurück."""
        g = Graph()
        for i in [3, 1, 4, 1, 5]:
            g.add_vertex(i)
        verts = g.vertices()
        assert isinstance(verts, list)
        # Keine Duplikate
        assert len(verts) == len(set(verts))

    def test_edges_undirected(self):
        """edges() gibt jede Kante einmal zurück (ungerichtet)."""
        g = Graph(directed=False)
        g.add_edge(1, 2)
        g.add_edge(2, 3)
        g.add_edge(1, 3)
        edges = g.edges()
        assert len(edges) == 3

    def test_degree_calculation(self):
        """Grad eines Knotens."""
        g = Graph()
        g.add_edge(1, 2)
        g.add_edge(1, 3)
        g.add_edge(1, 4)
        assert g.degree(1) == 3
        assert g.degree(2) == 1

    def test_degree_sequence(self):
        """Gradfolge ist absteigend sortiert."""
        g = complete_graph(4)
        seq = g.degree_sequence()
        assert seq == sorted(seq, reverse=True)
        assert all(d == 3 for d in seq)  # K_4: alle Grad 3

    def test_adjacency_matrix_shape(self):
        """Adjazenzmatrix hat korrekte Größe."""
        g = Graph()
        g.add_edge(0, 1)
        g.add_edge(1, 2)
        matrix = g.adjacency_matrix()
        assert len(matrix) == 3
        assert all(len(row) == 3 for row in matrix)

    def test_adjacency_matrix_symmetric(self):
        """Adjazenzmatrix ist symmetrisch für ungerichtete Graphen."""
        g = path_graph(4)
        m = g.adjacency_matrix()
        n = len(m)
        for i in range(n):
            for j in range(n):
                assert m[i][j] == m[j][i]


# ============================================================
# Tests: Grapheneigenschaften
# ============================================================

class TestGraphProperties:
    """Tests für topologische Eigenschaften von Graphen."""

    def test_is_connected_path(self):
        """Pfadgraph ist zusammenhängend."""
        g = path_graph(5)
        assert g.is_connected() is True

    def test_is_connected_isolated(self):
        """Isolierte Knoten → nicht zusammenhängend."""
        g = Graph()
        g.add_vertex(1)
        g.add_vertex(2)
        # Keine Kanten → 2 Komponenten
        assert g.is_connected() is False

    def test_connected_components_count(self):
        """Anzahl der Zusammenhangskomponenten."""
        g = Graph()
        g.add_edge(1, 2)
        g.add_edge(3, 4)
        # Zwei Komponenten: {1,2} und {3,4}
        components = g.connected_components()
        assert len(components) == 2

    def test_has_cycle_true(self):
        """Kreisgraph C_3 hat einen Zyklus."""
        g = cycle_graph(3)
        assert g.has_cycle() is True

    def test_has_cycle_false(self):
        """Baum (Pfadgraph) hat keinen Zyklus."""
        g = path_graph(5)
        assert g.has_cycle() is False

    def test_is_tree_true(self):
        """Pfadgraph ist ein Baum."""
        g = path_graph(5)
        assert g.is_tree() is True

    def test_is_tree_false_cycle(self):
        """Kreisgraph ist kein Baum (hat Zyklus)."""
        g = cycle_graph(4)
        assert g.is_tree() is False

    def test_is_bipartite_bipartite(self):
        """Bipartiter vollständiger Graph ist bipartit."""
        g = bipartite_complete_graph(3, 3)
        assert g.is_bipartite() is True

    def test_is_bipartite_cycle_even(self):
        """Geradzahliger Kreis C_4 ist bipartit."""
        g = cycle_graph(4)
        assert g.is_bipartite() is True

    def test_is_bipartite_cycle_odd(self):
        """Ungeradzahliger Kreis C_3 ist nicht bipartit."""
        g = cycle_graph(3)
        assert g.is_bipartite() is False

    def test_complement(self):
        """Komplementgraph: K_n und sein Komplement haben disjunkte Kantenmengen."""
        g = complete_graph(4)
        comp = g.complement()
        # Komplement von K_4 hat keine Kanten
        assert len(comp.edges()) == 0

    def test_subgraph(self):
        """Teilgraph enthält nur die angegebenen Knoten."""
        g = complete_graph(5)
        sub = g.subgraph([0, 1, 2])
        assert set(sub.vertices()) == {0, 1, 2}


# ============================================================
# Tests: BFS
# ============================================================

class TestBFS:
    """Tests für Breitensuche."""

    def test_bfs_distances_path_graph(self):
        """BFS auf Pfadgraph: korrekte Abstände."""
        g = path_graph(5)  # 0-1-2-3-4
        result = bfs(g, 0)
        assert result['distances'][0] == 0
        assert result['distances'][1] == 1
        assert result['distances'][2] == 2
        assert result['distances'][3] == 3
        assert result['distances'][4] == 4

    def test_bfs_visited_all(self):
        """BFS besucht alle Knoten in zusammenhängendem Graphen."""
        g = complete_graph(5)
        result = bfs(g, 0)
        assert len(result['visited']) == 5

    def test_bfs_parent_tree(self):
        """BFS-Vorgänger bilden einen Baum."""
        g = path_graph(4)
        result = bfs(g, 0)
        assert result['parent'][0] is None
        assert result['parent'][1] == 0
        assert result['parent'][2] == 1
        assert result['parent'][3] == 2

    def test_bfs_start_distance_zero(self):
        """Abstand zum Startknoten ist 0."""
        g = cycle_graph(5)
        result = bfs(g, 2)
        assert result['distances'][2] == 0


# ============================================================
# Tests: DFS
# ============================================================

class TestDFS:
    """Tests für Tiefensuche."""

    def test_dfs_visits_all(self):
        """DFS besucht alle Knoten im zusammenhängenden Graphen."""
        g = path_graph(5)
        result = dfs(g, 0)
        assert len(result['visited']) == 5

    def test_dfs_discovery_times(self):
        """Entdeckungszeitpunkte sind eindeutig und positiv."""
        g = path_graph(4)
        result = dfs(g, 0)
        times = list(result['discovery'].values())
        assert len(times) == len(set(times))  # Eindeutig
        assert all(t > 0 for t in times)


# ============================================================
# Tests: Dijkstra
# ============================================================

class TestDijkstra:
    """Tests für den Dijkstra-Algorithmus."""

    def test_dijkstra_path_graph(self):
        """Dijkstra auf Pfadgraph findet kürzeste Wege."""
        g = path_graph(5)
        result = dijkstra(g, 0)
        assert result['distances'][0] == 0
        assert result['distances'][1] == 1
        assert result['distances'][4] == 4

    def test_dijkstra_weighted(self):
        """Dijkstra mit gewichteten Kanten."""
        g = Graph(directed=False)
        g.add_edge(0, 1, weight=1.0)
        g.add_edge(0, 2, weight=4.0)
        g.add_edge(1, 2, weight=2.0)
        g.add_edge(1, 3, weight=5.0)
        g.add_edge(2, 3, weight=1.0)
        # Kürzester Weg 0→3: 0→1→2→3, Länge = 1+2+1 = 4
        result = dijkstra(g, 0)
        assert abs(result['distances'][3] - 4.0) < 1e-10

    def test_dijkstra_path_reconstruction(self):
        """Pfadrekonstruktion liefert korrekten Pfad."""
        g = path_graph(4)  # 0-1-2-3
        result = dijkstra(g, 0)
        assert result['path'][3] == [0, 1, 2, 3]

    def test_dijkstra_unreachable(self):
        """Unerreichbarer Knoten hat Abstand inf."""
        g = Graph()
        g.add_vertex(0)
        g.add_vertex(1)
        # Keine Kante zwischen 0 und 1
        result = dijkstra(g, 0)
        assert result['distances'][1] == float('inf')


# ============================================================
# Tests: Bellman-Ford
# ============================================================

class TestBellmanFord:
    """Tests für den Bellman-Ford-Algorithmus."""

    def test_bellman_ford_simple(self):
        """Bellman-Ford auf einfachem Graphen."""
        g = path_graph(4)
        result = bellman_ford(g, 0)
        assert result['distances'][3] == 3
        assert result['negative_cycle'] is False

    def test_bellman_ford_negative_weights(self):
        """Bellman-Ford mit negativen Gewichten."""
        g = Graph(directed=True)
        g.add_edge(0, 1, weight=4.0)
        g.add_edge(0, 2, weight=5.0)
        g.add_edge(1, 3, weight=3.0)
        g.add_edge(2, 1, weight=-6.0)  # Negative Kante!
        # Kürzester Weg 0→1: 0→2→1 = 5 + (-6) = -1
        result = bellman_ford(g, 0)
        assert result['distances'][1] == -1.0
        assert result['negative_cycle'] is False


# ============================================================
# Tests: Floyd-Warshall
# ============================================================

class TestFloydWarshall:
    """Tests für den Floyd-Warshall-Algorithmus."""

    def test_floyd_warshall_path_graph(self):
        """Floyd-Warshall auf Pfadgraph."""
        g = path_graph(3)  # 0-1-2
        dist = floyd_warshall(g)
        # Reihenfolge der Knoten: 0, 1, 2
        # d(0,2) = 2
        assert dist[0][2] == 2.0
        assert dist[2][0] == 2.0

    def test_floyd_warshall_diagonal_zero(self):
        """Diagonale der Distanzmatrix ist 0."""
        g = complete_graph(4)
        dist = floyd_warshall(g)
        for i in range(len(dist)):
            assert dist[i][i] == 0.0


# ============================================================
# Tests: Kruskal und Prim
# ============================================================

class TestSpanningTrees:
    """Tests für minimale Spannbäume."""

    def test_kruskal_mst_edge_count(self):
        """MST hat genau V-1 Kanten."""
        g = complete_graph(5)
        mst = kruskal_mst(g)
        n_edges = len(mst.edges())
        assert n_edges == 4  # V - 1 = 5 - 1 = 4

    def test_kruskal_mst_is_tree(self):
        """Kruskal-Ergebnis ist ein Baum."""
        g = Graph()
        g.add_edge(0, 1, weight=2.0)
        g.add_edge(0, 2, weight=3.0)
        g.add_edge(1, 2, weight=1.0)
        g.add_edge(1, 3, weight=4.0)
        mst = kruskal_mst(g)
        assert mst.is_tree() is True

    def test_prim_mst_edge_count(self):
        """Prim ergibt ebenfalls V-1 Kanten."""
        g = complete_graph(5)
        mst = prim_mst(g)
        assert len(mst.edges()) == 4

    def test_prim_mst_is_connected(self):
        """Prim-MST ist zusammenhängend."""
        g = Graph()
        g.add_edge(0, 1, weight=1.0)
        g.add_edge(1, 2, weight=2.0)
        g.add_edge(0, 2, weight=5.0)
        mst = prim_mst(g)
        assert mst.is_connected() is True


# ============================================================
# Tests: Topologische Sortierung
# ============================================================

class TestTopologicalSort:
    """Tests für topologische Sortierung."""

    def test_topological_sort_dag(self):
        """Topologische Sortierung eines DAGs."""
        g = Graph(directed=True)
        g.add_edge(0, 1)
        g.add_edge(0, 2)
        g.add_edge(1, 3)
        g.add_edge(2, 3)
        order = topological_sort(g)
        assert len(order) == 4
        # 0 muss vor 1 und 2 kommen; 1 und 2 müssen vor 3 kommen
        assert order.index(0) < order.index(1)
        assert order.index(0) < order.index(2)
        assert order.index(1) < order.index(3)
        assert order.index(2) < order.index(3)

    def test_topological_sort_cycle_returns_empty(self):
        """Topologische Sortierung mit Zyklus gibt [] zurück."""
        g = Graph(directed=True)
        g.add_edge(0, 1)
        g.add_edge(1, 2)
        g.add_edge(2, 0)  # Zyklus!
        order = topological_sort(g)
        assert order == []

    def test_topological_sort_linear(self):
        """Linearer DAG: einzige Reihenfolge."""
        g = Graph(directed=True)
        g.add_edge(0, 1)
        g.add_edge(1, 2)
        g.add_edge(2, 3)
        order = topological_sort(g)
        assert order == [0, 1, 2, 3]


# ============================================================
# Tests: Graphenanalyse
# ============================================================

class TestGraphAnalysis:
    """Tests für Graphenanalyse-Algorithmen."""

    def test_chromatic_number_complete_k4(self):
        """K_4 benötigt 4 Farben: χ(K_4) = 4."""
        g = complete_graph(4)
        assert chromatic_number_greedy(g) == 4

    def test_chromatic_number_bipartite(self):
        """Bipartiter Graph kann mit 2 Farben gefärbt werden."""
        g = bipartite_complete_graph(3, 3)
        assert chromatic_number_greedy(g) == 2

    def test_chromatic_number_path(self):
        """Pfadgraph braucht ≤ 2 Farben (bipartit)."""
        g = path_graph(6)
        assert chromatic_number_greedy(g) <= 2

    def test_graph_coloring_valid(self):
        """Keine zwei benachbarten Knoten haben dieselbe Farbe."""
        g = complete_graph(5)
        coloring = graph_coloring_greedy(g)
        for u in g.adj:
            for v in g.adj[u]:
                assert coloring[u] != coloring[v], \
                    f"Knoten {u} und {v} haben beide Farbe {coloring[u]}!"

    def test_graph_coloring_cycle_3(self):
        """C_3 braucht genau 3 Farben."""
        g = cycle_graph(3)
        coloring = graph_coloring_greedy(g)
        assert len(set(coloring.values())) == 3

    def test_is_eulerian_complete_k4(self):
        """K_4: alle Knotengrade = 3 (ungerade) → kein Euler-Kreis, kein Euler-Pfad."""
        g = complete_graph(4)  # Grad 3 für alle Knoten
        result = is_eulerian(g)
        # K_4: 4 Knoten haben ungeraden Grad → kein Euler-Kreis, kein Euler-Pfad (4 ≠ 2)
        assert result['euler_circuit'] is False
        assert result['euler_path'] is False

    def test_is_eulerian_circuit_exists(self):
        """Graph mit allen geraden Graden hat Euler-Kreis."""
        # C_4 hat alle Knoten mit Grad 2 (gerade)
        g = cycle_graph(4)
        result = is_eulerian(g)
        assert result['euler_circuit'] is True
        assert result['euler_path'] is False

    def test_is_eulerian_path_exists(self):
        """Graph mit genau 2 Knoten ungeradem Grad hat Euler-Pfad."""
        # Pfadgraph P_4: Endknoten haben Grad 1 (ungerade), die anderen Grad 2
        g = path_graph(4)
        result = is_eulerian(g)
        assert result['euler_path'] is True
        assert result['euler_circuit'] is False
        assert len(result['odd_degree_vertices']) == 2

    def test_is_hamiltonian_path_complete(self):
        """K_4 hat einen Hamiltonpfad."""
        g = complete_graph(4)
        result = is_hamiltonian_path(g)
        assert result['has_path'] is True
        assert len(result['path']) == 4

    def test_is_hamiltonian_path_single(self):
        """Graph mit einem Knoten: trivialer Hamiltonpfad."""
        g = Graph()
        g.add_vertex(0)
        result = is_hamiltonian_path(g)
        assert result['has_path'] is True

    def test_clique_number_complete(self):
        """Kliquenzahl von K_n = n."""
        assert clique_number(complete_graph(4)) == 4

    def test_clique_number_path(self):
        """Pfadgraph hat Kliquenzahl 2 (nur Kanten)."""
        g = path_graph(5)
        assert clique_number(g) == 2

    def test_independence_number_complete(self):
        """α(K_n) = 1 (maximale unabhängige Menge: ein Knoten)."""
        assert independence_number(complete_graph(4)) == 1

    def test_independence_number_path(self):
        """Unabhängigkeitszahl des Pfadgraphen P_4."""
        g = path_graph(4)  # 0-1-2-3 → max. unabhängige Menge: {0,2} oder {1,3} → α=2
        assert independence_number(g) == 2

    def test_graph_density_complete(self):
        """Dichte von K_n = 1."""
        assert abs(graph_density(complete_graph(5)) - 1.0) < 1e-10

    def test_graph_density_empty(self):
        """Dichte eines Graphen ohne Kanten = 0."""
        g = Graph()
        for i in range(5):
            g.add_vertex(i)
        assert graph_density(g) == 0.0


# ============================================================
# Tests: Graphen-Konstruktoren
# ============================================================

class TestGraphConstructors:
    """Tests für vorgefertigte Graphen."""

    def test_complete_graph_k4_edges(self):
        """K_4 hat genau 6 Kanten."""
        g = complete_graph(4)
        assert len(g.edges()) == 6

    def test_complete_graph_k5_edges(self):
        """K_5 hat genau 10 Kanten: C(5,2) = 10."""
        g = complete_graph(5)
        assert len(g.edges()) == 10

    def test_complete_graph_degree(self):
        """Jeder Knoten in K_n hat Grad n-1."""
        g = complete_graph(5)
        for v in g.vertices():
            assert g.degree(v) == 4

    def test_cycle_graph_edges(self):
        """C_n hat genau n Kanten."""
        g = cycle_graph(6)
        assert len(g.edges()) == 6

    def test_cycle_graph_degree(self):
        """Jeder Knoten in C_n hat Grad 2."""
        g = cycle_graph(5)
        for v in g.vertices():
            assert g.degree(v) == 2

    def test_path_graph_edges(self):
        """P_n hat n-1 Kanten."""
        g = path_graph(5)
        assert len(g.edges()) == 4

    def test_path_graph_endpoints_degree_1(self):
        """Endknoten des Pfadgraphen haben Grad 1."""
        g = path_graph(4)
        assert g.degree(0) == 1
        assert g.degree(3) == 1

    def test_bipartite_complete_edges(self):
        """K_{m,n} hat m*n Kanten."""
        g = bipartite_complete_graph(3, 3)
        assert len(g.edges()) == 9

    def test_bipartite_complete_is_bipartite(self):
        """K_{m,n} ist bipartit."""
        g = bipartite_complete_graph(3, 4)
        assert g.is_bipartite() is True

    def test_petersen_graph_vertices_edges(self):
        """Petersen-Graph: 10 Knoten, 15 Kanten."""
        g = petersen_graph()
        assert len(g.vertices()) == 10
        assert len(g.edges()) == 15

    def test_petersen_graph_3_regular(self):
        """Petersen-Graph ist 3-regulär."""
        g = petersen_graph()
        for v in g.vertices():
            assert g.degree(v) == 3

    def test_grid_graph_nodes(self):
        """m×n Gittergraph hat m*n Knoten."""
        g = grid_graph(3, 4)
        assert len(g.vertices()) == 12

    def test_grid_graph_is_bipartite(self):
        """Gittergraph ist bipartit (Schachbrettfärbung)."""
        g = grid_graph(4, 4)
        assert g.is_bipartite() is True


# ============================================================
# Tests: Kombinatorik
# ============================================================

class TestCombinatorics:
    """Tests für kombinatorische Formeln."""

    def test_binomial_basic(self):
        """C(10, 3) = 120."""
        assert binomial_coefficient(10, 3) == 120

    def test_binomial_zero(self):
        """C(n, 0) = 1."""
        assert binomial_coefficient(10, 0) == 1
        assert binomial_coefficient(0, 0) == 1

    def test_binomial_n_choose_n(self):
        """C(n, n) = 1."""
        assert binomial_coefficient(7, 7) == 1

    def test_binomial_symmetry(self):
        """C(n, k) = C(n, n-k)."""
        assert binomial_coefficient(10, 3) == binomial_coefficient(10, 7)

    def test_binomial_pascal(self):
        """Pascalsches Dreieck: C(n,k) = C(n-1,k-1) + C(n-1,k)."""
        n, k = 8, 3
        assert binomial_coefficient(n, k) == (binomial_coefficient(n-1, k-1) +
                                               binomial_coefficient(n-1, k))

    def test_binomial_invalid_k(self):
        """C(n, k) = 0 für k > n."""
        assert binomial_coefficient(5, 7) == 0

    def test_stirling_s21(self):
        """S(2,1) = 1 (eine Partition von 2 Elementen in 1 Block)."""
        assert stirling_numbers_second_kind(2, 1) == 1

    def test_stirling_s22(self):
        """S(2,2) = 1."""
        assert stirling_numbers_second_kind(2, 2) == 1

    def test_stirling_s32(self):
        """S(3,2) = 3 (drei Wege 3 Elemente in 2 Blöcke zu teilen)."""
        assert stirling_numbers_second_kind(3, 2) == 3

    def test_stirling_s43(self):
        """S(4,3) = 6."""
        assert stirling_numbers_second_kind(4, 3) == 6

    def test_bell_numbers_first_five(self):
        """Bell-Zahlen B_0 bis B_4."""
        b = bell_numbers(4)
        assert b == [1, 1, 2, 5, 15]

    def test_bell_numbers_b5(self):
        """B_5 = 52."""
        b = bell_numbers(5)
        assert b[5] == 52

    def test_bell_length(self):
        """bell_numbers(n) gibt n+1 Werte zurück."""
        b = bell_numbers(5)
        assert len(b) == 6

    def test_catalan_c5(self):
        """C_5 = 42."""
        assert catalan_number(5) == 42

    def test_catalan_c0(self):
        """C_0 = 1."""
        assert catalan_number(0) == 1

    def test_catalan_first_values(self):
        """Erste Catalan-Zahlen: 1, 1, 2, 5, 14, 42."""
        expected = [1, 1, 2, 5, 14, 42]
        for i, val in enumerate(expected):
            assert catalan_number(i) == val

    def test_derangements_d0(self):
        """!0 = 1 (leere Permutation ist Derangement)."""
        assert derangements(0) == 1

    def test_derangements_d1(self):
        """!1 = 0 (die einzige Permutation ist der Fixpunkt)."""
        assert derangements(1) == 0

    def test_derangements_d4(self):
        """!4 = 9."""
        assert derangements(4) == 9

    def test_derangements_d5(self):
        """!5 = 44."""
        assert derangements(5) == 44

    def test_partition_count_p0(self):
        """p(0) = 1 (leere Partition)."""
        assert partition_count(0) == 1

    def test_partition_count_p4(self):
        """p(4) = 5 (vier Partitionen von 4)."""
        assert partition_count(4) == 5

    def test_partition_count_p10(self):
        """p(10) = 42."""
        assert partition_count(10) == 42

    def test_multinomial_basic(self):
        """Multinomialkoeffizient: 4! / (2! * 1! * 1!) = 12."""
        assert multinomial_coefficient(4, [2, 1, 1]) == 12

    def test_multinomial_binomial(self):
        """Spezialfall p=2: entspricht Binomialkoeffizient."""
        n, k = 10, 3
        assert multinomial_coefficient(n, [k, n - k]) == binomial_coefficient(n, k)

    def test_multinomial_wrong_sum_raises(self):
        """Falsche Gruppengrößen → ValueError."""
        with pytest.raises(ValueError):
            multinomial_coefficient(5, [2, 2])  # 2+2 = 4 ≠ 5


# ============================================================
# Tests: Erweiterte Edge Cases
# ============================================================

class TestEdgeCases:
    """Tests für Randfälle und besondere Situationen."""

    def test_empty_graph(self):
        """Leerer Graph verhält sich korrekt."""
        g = Graph()
        assert g.vertices() == []
        assert g.edges() == []
        assert g.is_connected() is True  # Leerer Graph: trivial zusammenhängend

    def test_single_vertex_graph(self):
        """Graph mit einem Knoten."""
        g = Graph()
        g.add_vertex(42)
        assert g.is_connected() is True
        assert g.has_cycle() is False
        assert g.is_tree() is True

    def test_self_loop_detection(self):
        """Graph mit Schleifen (Kante von v zu v)."""
        g = Graph(directed=True)
        g.add_edge(0, 0)  # Selbstschleife
        assert g.has_cycle() is True

    def test_large_complete_graph_density(self):
        """K_10: Dichte = 1."""
        g = complete_graph(10)
        assert abs(graph_density(g) - 1.0) < 1e-10

    def test_binomial_large(self):
        """Binomialkoeffizient für größere Werte."""
        # C(20, 10) = 184756
        assert binomial_coefficient(20, 10) == 184756

    def test_catalan_large(self):
        """Catalan-Zahl für n=10."""
        # C_10 = 16796
        assert catalan_number(10) == 16796
