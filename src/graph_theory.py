"""
@file graph_theory.py
@brief Graphentheorie – Graphen, Traversierungsalgorithmen, Spanning Trees,
       Graphenanalyse und kombinatorische Formeln.
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import heapq
import math
from typing import Any


# ============================================================
# Graph-Klasse
# ============================================================

class Graph:
    """
    @brief Ungerichteter oder gerichteter Graph G = (V, E).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Interne Darstellung: Adjazenzliste (dict: vertex → set of neighbors).
    Kantengewichte werden in einem separaten Dict gespeichert.

    Mathematisch: G = (V, E) mit V = Knotenmenge, E ⊆ V × V = Kantenmenge.

    Für ungerichtete Graphen gilt: (u,v) ∈ E ⟺ (v,u) ∈ E.
    Für gerichtete Graphen (Digraphen): (u,v) ≠ (v,u) möglich.
    """

    def __init__(self, directed: bool = False):
        """
        @brief Initialisiert einen leeren Graphen.
        @param directed True für gerichteten Graphen, False für ungerichteten.
        """
        self.directed = directed                        # Gerichtet oder ungerichtet
        self.adj: dict[Any, set] = {}                  # Adjazenzliste: {knoten: {nachbarn}}
        self.weights: dict[tuple, float] = {}          # Kantengewichte: {(u,v): weight}

    def add_vertex(self, v: int | str) -> None:
        """
        @brief Fügt einen Knoten zum Graphen hinzu.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param v Knoten (beliebiger hashbarer Typ: int, str, tuple, ...).
        """
        # Neuen Knoten mit leerer Nachbarschaftsmenge anlegen
        if v not in self.adj:
            self.adj[v] = set()

    def add_edge(self, u: int | str, v: int | str, weight: float = 1.0) -> None:
        """
        @brief Fügt eine Kante (u,v) mit optionalem Gewicht hinzu.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bei ungerichteten Graphen wird auch (v,u) hinzugefügt.
        Knoten werden automatisch angelegt falls noch nicht vorhanden.

        @param u Startknoten.
        @param v Endknoten.
        @param weight Kantengewicht (Standard: 1.0).
        """
        # Knoten anlegen falls noch nicht vorhanden
        self.add_vertex(u)
        self.add_vertex(v)

        # Kante eintragen
        self.adj[u].add(v)
        self.weights[(u, v)] = weight

        # Bei ungerichtetem Graph: Kante in beide Richtungen
        if not self.directed:
            self.adj[v].add(u)
            self.weights[(v, u)] = weight

    def remove_vertex(self, v: int | str) -> None:
        """
        @brief Entfernt einen Knoten und alle angrenzenden Kanten.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param v Zu entfernender Knoten.
        """
        if v not in self.adj:
            return

        # Alle Kanten zu v aus den Nachbarschaftslisten entfernen
        for neighbor in list(self.adj[v]):
            self.adj[neighbor].discard(v)
            self.weights.pop((neighbor, v), None)
            self.weights.pop((v, neighbor), None)

        # Knoten selbst entfernen
        del self.adj[v]

    def remove_edge(self, u: int | str, v: int | str) -> None:
        """
        @brief Entfernt die Kante (u,v).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bei ungerichteten Graphen wird auch (v,u) entfernt.

        @param u Startknoten.
        @param v Endknoten.
        """
        if u in self.adj:
            self.adj[u].discard(v)
        self.weights.pop((u, v), None)

        if not self.directed:
            if v in self.adj:
                self.adj[v].discard(u)
            self.weights.pop((v, u), None)

    def vertices(self) -> list[int | str]:
        """
        @brief Gibt alle Knoten des Graphen zurück.
        @return Sortierte Liste aller Knoten.
        """
        return sorted(self.adj.keys(), key=lambda x: (str(type(x)), x))

    def edges(self) -> list[tuple]:
        """
        @brief Gibt alle Kanten als Liste von Tupeln zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bei ungerichteten Graphen wird jede Kante nur einmal zurückgegeben.

        @return Liste von (u, v, weight)-Tupeln.
        """
        result = []
        seen = set()

        for u in self.adj:
            for v in self.adj[u]:
                # Bei ungerichteten Graphen: jede Kante nur einmal
                if not self.directed:
                    edge_key = (min(str(u), str(v)), max(str(u), str(v)))
                    if edge_key in seen:
                        continue
                    seen.add(edge_key)
                w = self.weights.get((u, v), 1.0)
                result.append((u, v, w))

        return result

    def neighbors(self, v: int | str) -> list[int | str]:
        """
        @brief Gibt die Nachbarn eines Knotens zurück.
        @param v Knoten.
        @return Sortierte Liste der Nachbarn.
        """
        if v not in self.adj:
            return []
        return sorted(self.adj[v], key=lambda x: (str(type(x)), x))

    def degree(self, v: int | str) -> int:
        """
        @brief Grad (Anzahl der Kanten) eines Knotens.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Bei gerichteten Graphen: Ausgangsgrad (out-degree).

        @param v Knoten.
        @return Grad des Knotens.
        """
        if v not in self.adj:
            return 0
        return len(self.adj[v])

    def degree_sequence(self) -> list[int]:
        """
        @brief Sortierte Liste aller Knotengrade (absteigend).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Gradfolge ist eine topologische Invariante des Graphen.
        Erdős-Gallai-Theorem: Charakterisiert realisierbare Gradfolgen.

        @return Absteigend sortierte Liste der Knotengrade.
        """
        return sorted([self.degree(v) for v in self.adj], reverse=True)

    def adjacency_matrix(self) -> list[list[float]]:
        """
        @brief Adjazenzmatrix als 2D-Liste.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        A[i][j] = Gewicht der Kante (v_i, v_j), oder 0 wenn keine Kante.

        @return 2D-Liste der Größe n×n.
        """
        verts = self.vertices()
        n = len(verts)
        idx = {v: i for i, v in enumerate(verts)}

        # Matrix mit Nullen initialisieren
        matrix = [[0.0] * n for _ in range(n)]

        for u in self.adj:
            for v in self.adj[u]:
                # Kantengewicht eintragen
                matrix[idx[u]][idx[v]] = self.weights.get((u, v), 1.0)

        return matrix

    def adjacency_list(self) -> dict[int | str, list[int | str]]:
        """
        @brief Adjazenzliste als Dict.
        @return Dict {knoten: sorted_liste_der_nachbarn}.
        """
        return {v: self.neighbors(v) for v in self.adj}

    def is_connected(self) -> bool:
        """
        @brief Prüft ob der Graph zusammenhängend ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Verwendet BFS von einem Startknoten. Ein Graph ist zusammenhängend,
        wenn von jedem Knoten aus jeder andere Knoten erreichbar ist.

        @return True wenn der Graph zusammenhängend ist.
        """
        if not self.adj:
            return True  # Leerer Graph gilt als zusammenhängend

        # BFS von einem beliebigen Startknoten
        start = next(iter(self.adj))
        visited = self._bfs_visited(start)

        # Alle Knoten müssen besucht worden sein
        return len(visited) == len(self.adj)

    def _bfs_visited(self, start: int | str) -> set:
        """
        @brief Interne BFS-Hilfsmethode – gibt Menge besuchter Knoten zurück.
        @param start Startknoten.
        @return Menge aller von start erreichbaren Knoten.
        """
        visited = {start}
        queue = [start]

        while queue:
            node = queue.pop(0)
            for neighbor in self.adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        return visited

    def connected_components(self) -> list[set]:
        """
        @brief Findet alle Zusammenhangskomponenten.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Verwendet wiederholte BFS für unbesuchte Knoten.

        @return Liste von Mengen, eine Menge pro Komponente.
        """
        visited_all = set()
        components = []

        for v in self.adj:
            if v not in visited_all:
                # Neue Komponente via BFS erkunden
                component = self._bfs_visited(v)
                components.append(component)
                visited_all |= component

        return components

    def has_cycle(self) -> bool:
        """
        @brief Prüft ob der Graph einen Zyklus enthält.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        DFS-basierte Zykluserkennung:
        - Ungerichtet: Rückkante zu einem bereits besuchten Knoten (nicht Vater)
        - Gerichtet: Rückkante im DFS-Baum

        @return True wenn ein Zyklus vorhanden ist.
        """
        visited = set()
        rec_stack = set()  # Nur für gerichtete Graphen: Rekursionskeller

        def dfs_has_cycle(v: Any, parent: Any) -> bool:
            """DFS Hilfsfunktion zur Zykluserkennung."""
            visited.add(v)
            if self.directed:
                rec_stack.add(v)

            for neighbor in self.adj[v]:
                if neighbor not in visited:
                    if dfs_has_cycle(neighbor, v):
                        return True
                elif self.directed and neighbor in rec_stack:
                    # Gerichtet: Rückkante im Rekursionskeller → Zyklus
                    return True
                elif not self.directed and neighbor != parent:
                    # Ungerichtet: Besuchter Nachbar der nicht Vater ist → Zyklus
                    return True

            if self.directed:
                rec_stack.discard(v)
            return False

        for v in self.adj:
            if v not in visited:
                if dfs_has_cycle(v, None):
                    return True

        return False

    def is_tree(self) -> bool:
        """
        @brief Prüft ob der Graph ein Baum ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Baum ist ein zusammenhängender, azyklischer Graph.
        Äquivalent: n Knoten, n-1 Kanten, zusammenhängend.

        @return True wenn der Graph ein Baum ist.
        """
        return self.is_connected() and not self.has_cycle()

    def is_bipartite(self) -> bool:
        """
        @brief Prüft ob der Graph bipartit ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Graph ist bipartit, wenn seine Knoten in zwei disjunkte Mengen A, B
        aufgeteilt werden können, sodass alle Kanten zwischen A und B verlaufen.

        Äquivalent: 2-Färbbarkeit (via BFS) oder: kein ungerader Zyklus.

        @return True wenn der Graph bipartit ist.
        """
        # 2-Färbung via BFS
        color = {}

        for start in self.adj:
            if start in color:
                continue

            color[start] = 0
            queue = [start]

            while queue:
                node = queue.pop(0)
                for neighbor in self.adj[node]:
                    if neighbor not in color:
                        # Nachbar mit anderer Farbe einfärben
                        color[neighbor] = 1 - color[node]
                        queue.append(neighbor)
                    elif color[neighbor] == color[node]:
                        # Gleiche Farbe → nicht bipartit
                        return False

        return True

    def complement(self) -> 'Graph':
        """
        @brief Erzeugt den Komplementgraphen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Im Komplementgraphen G̅ sind genau diejenigen Paare verbunden,
        die in G nicht verbunden sind:
        E(G̅) = {(u,v) : u≠v und (u,v) ∉ E(G)}

        @return Komplementgraph als neues Graph-Objekt.
        """
        comp = Graph(directed=self.directed)
        verts = self.vertices()

        # Alle Knoten hinzufügen
        for v in verts:
            comp.add_vertex(v)

        # Alle fehlenden Kanten hinzufügen
        for i, u in enumerate(verts):
            for j, v in enumerate(verts):
                if u != v:
                    if self.directed:
                        if v not in self.adj.get(u, set()):
                            comp.add_edge(u, v)
                    else:
                        if j > i and v not in self.adj.get(u, set()):
                            comp.add_edge(u, v)

        return comp

    def subgraph(self, vertices: list[int | str]) -> 'Graph':
        """
        @brief Erzeugt den durch vertices induzierten Teilgraphen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Der induzierte Teilgraph enthält alle Kanten aus G,
        deren beide Endpunkte in vertices liegen.

        @param vertices Liste der Knoten des Teilgraphen.
        @return Induzierter Teilgraph als neues Graph-Objekt.
        """
        sub = Graph(directed=self.directed)
        vertex_set = set(vertices)

        # Knoten hinzufügen
        for v in vertices:
            sub.add_vertex(v)

        # Nur Kanten innerhalb der Knotenmenge hinzufügen
        for u in vertices:
            for v in self.adj.get(u, set()):
                if v in vertex_set:
                    w = self.weights.get((u, v), 1.0)
                    if self.directed or u < v or str(u) < str(v):
                        sub.add_edge(u, v, w)

        return sub


# ============================================================
# Traversierung und kürzeste Wege
# ============================================================

def bfs(graph: Graph, start: int | str) -> dict[str, dict]:
    """
    @brief Breitensuche (Breadth-First Search).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    BFS erkundet alle Knoten schichtweise nach ihrer Entfernung vom Startknoten.
    Findet kürzeste Wege (in ungewichteten Graphen).

    Zeitkomplexität: O(V + E)
    Raumkomplexität: O(V)

    @param graph Der Graph.
    @param start Startknoten.
    @return Dict mit 'visited' (Liste), 'distances' (Dict), 'parent' (Dict).
    """
    visited = []         # Reihenfolge der Besuchung
    distances = {}       # Abstände vom Startknoten
    parent = {}          # Vorgänger im BFS-Baum

    # Startknoten initialisieren
    distances[start] = 0
    parent[start] = None
    queue = [start]
    visited_set = {start}

    while queue:
        # Vordersten Knoten aus der Warteschlange nehmen
        node = queue.pop(0)
        visited.append(node)

        # Alle Nachbarn besuchen
        for neighbor in sorted(graph.adj.get(node, set()), key=lambda x: (str(type(x)), x)):
            if neighbor not in visited_set:
                visited_set.add(neighbor)
                distances[neighbor] = distances[node] + 1
                parent[neighbor] = node
                queue.append(neighbor)

    return {
        'visited': visited,
        'distances': distances,
        'parent': parent
    }


def dfs(graph: Graph, start: int | str) -> dict[str, dict]:
    """
    @brief Tiefensuche (Depth-First Search).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    DFS erkundet Pfade so tief wie möglich, bevor es zurückspringt.
    Nützlich für Zykluserkennung, topologische Sortierung, etc.

    Zeitkomplexität: O(V + E)
    Raumkomplexität: O(V) (Rekursionstiefe)

    @param graph Der Graph.
    @param start Startknoten.
    @return Dict mit 'visited' (Liste), 'discovery' (Dict), 'finish' (Dict).
    """
    visited = []       # Reihenfolge der Entdeckung
    discovery = {}     # Entdeckungszeitpunkt
    finish = {}        # Abschluss-Zeitpunkt
    time_counter = [0]  # Mutable für Closure

    def dfs_visit(v: Any) -> None:
        """Rekursive DFS-Hilfsfunktion."""
        time_counter[0] += 1
        discovery[v] = time_counter[0]
        visited.append(v)

        # Nachbarn in sortierter Reihenfolge besuchen
        for neighbor in sorted(graph.adj.get(v, set()), key=lambda x: (str(type(x)), x)):
            if neighbor not in discovery:
                dfs_visit(neighbor)

        time_counter[0] += 1
        finish[v] = time_counter[0]

    # DFS starten (auch für nicht zusammenhängende Graphen)
    if start in graph.adj:
        dfs_visit(start)

    return {
        'visited': visited,
        'discovery': discovery,
        'finish': finish
    }


def dijkstra(graph: Graph, start: int | str) -> dict[str, dict]:
    """
    @brief Dijkstra-Algorithmus für kürzeste Wege (nicht-negative Gewichte).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Findet die kürzesten Wege von einem Startknoten zu allen anderen Knoten.
    Funktioniert NICHT mit negativen Kantengewichten (→ Bellman-Ford).

    Algorithmus: Greedy mit MinHeap (Prioritätswarteschlange).
    Zeitkomplexität: O((V+E) log V)

    @param graph Der Graph mit nicht-negativen Kantengewichten.
    @param start Startknoten.
    @return Dict mit 'distances' (Dict) und 'path' (Dict: v → Pfad als Liste).
    """
    # Initialisierung: alle Abstände = unendlich
    distances = {v: float('inf') for v in graph.adj}
    distances[start] = 0
    parent = {v: None for v in graph.adj}

    # MinHeap: (aktueller_abstand, knoten)
    heap = [(0, start)]

    while heap:
        # Knoten mit kleinstem bisher bekannten Abstand wählen
        dist_u, u = heapq.heappop(heap)

        # Veralteten Eintrag überspringen
        if dist_u > distances[u]:
            continue

        # Alle Nachbarn relaxieren (Abstand ggf. verkürzen)
        for v in graph.adj.get(u, set()):
            weight = graph.weights.get((u, v), 1.0)
            new_dist = distances[u] + weight

            if new_dist < distances[v]:
                distances[v] = new_dist
                parent[v] = u
                heapq.heappush(heap, (new_dist, v))

    # Pfade rekonstruieren
    def reconstruct_path(target: Any) -> list:
        """Pfad vom Start zum Zielknoten via Vorgängerarray rekonstruieren."""
        if distances[target] == float('inf'):
            return []  # Kein Pfad vorhanden
        path = []
        current = target
        while current is not None:
            path.append(current)
            current = parent[current]
        return list(reversed(path))

    paths = {v: reconstruct_path(v) for v in graph.adj}

    return {'distances': distances, 'path': paths}


def bellman_ford(graph: Graph, start: int | str) -> dict[str, dict | bool]:
    """
    @brief Bellman-Ford-Algorithmus für kürzeste Wege.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Funktioniert auch mit negativen Kantengewichten und erkennt
    negative Zyklen.

    Algorithmus: Relaxiere alle |V|-1 mal alle Kanten.
    Falls danach noch eine Relaxierung möglich ist → negativer Zyklus.

    Zeitkomplexität: O(V·E)

    @param graph Der Graph (auch mit negativen Gewichten erlaubt).
    @param start Startknoten.
    @return Dict mit 'distances' (Dict) und 'negative_cycle' (bool).
    """
    # Initialisierung
    vertices = graph.vertices()
    distances = {v: float('inf') for v in vertices}
    distances[start] = 0
    parent = {v: None for v in vertices}

    # Alle Kanten (u,v) als Liste (für Iteration)
    edge_list = [(u, v, graph.weights.get((u, v), 1.0)) for u in graph.adj for v in graph.adj[u]]

    # V-1 Relaxierungsrunden
    for _ in range(len(vertices) - 1):
        updated = False
        for u, v, w in edge_list:
            if distances[u] != float('inf') and distances[u] + w < distances[v]:
                distances[v] = distances[u] + w
                parent[v] = u
                updated = True
        if not updated:
            break  # Frühzeitig abbrechen wenn keine Änderungen

    # Negativen Zyklus prüfen: weitere Relaxierung möglich?
    negative_cycle = False
    for u, v, w in edge_list:
        if distances[u] != float('inf') and distances[u] + w < distances[v]:
            negative_cycle = True
            break

    return {
        'distances': distances,
        'negative_cycle': negative_cycle
    }


def floyd_warshall(graph: Graph) -> list[list[float]]:
    """
    @brief Floyd-Warshall-Algorithmus: Alle kürzesten Wege.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Berechnet die kürzesten Wege zwischen ALLEN Knotenpaaren.
    Auch mit negativen Gewichten, aber nicht negativen Zyklen.

    Algorithmus: Dynamische Programmierung.
    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]) für alle k.

    Zeitkomplexität: O(V³)
    Raumkomplexität: O(V²)

    @param graph Der Graph.
    @return Distanzmatrix als 2D-Liste (Indizes entsprechen graph.vertices()).
    """
    verts = graph.vertices()
    n = len(verts)
    idx = {v: i for i, v in enumerate(verts)}

    # Distanzmatrix initialisieren
    INF = float('inf')
    dist = [[INF] * n for _ in range(n)]

    # Diagonale auf 0 setzen
    for i in range(n):
        dist[i][i] = 0.0

    # Bekannte Kanten eintragen
    for u in graph.adj:
        for v in graph.adj[u]:
            w = graph.weights.get((u, v), 1.0)
            dist[idx[u]][idx[v]] = w

    # Floyd-Warshall-Kern: für jeden Zwischenknoten k relaxieren
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i][k] + dist[k][j] < dist[i][j]:
                    dist[i][j] = dist[i][k] + dist[k][j]

    return dist


# ============================================================
# Spanning Trees und Spezielle Algorithmen
# ============================================================

def kruskal_mst(graph: Graph) -> Graph:
    """
    @brief Kruskal-Algorithmus: Minimaler Spannbaum.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Wählt iterativ die günstigste Kante, die keinen Zyklus erzeugt.
    Verwendet Union-Find für effiziente Zykluserkennung.

    Greedy-Eigenschaft: Lokal optimale Wahl → global optimaler Spannbaum.
    Korrektheit: Schnitteigentschaft und Zyklus-Eigenschaft.

    Zeitkomplexität: O(E log E) durch Sortierung der Kanten.

    @param graph Der Graph (ungerichtet, gewichtet).
    @return Minimaler Spannbaum als neues Graph-Objekt.
    """
    mst = Graph(directed=False)

    # Alle Knoten in den MST übernehmen
    for v in graph.vertices():
        mst.add_vertex(v)

    # Kanten aufsteigend nach Gewicht sortieren
    edges_sorted = sorted(graph.edges(), key=lambda e: e[2])

    # Union-Find-Struktur initialisieren
    verts = graph.vertices()
    parent = {v: v for v in verts}
    rank = {v: 0 for v in verts}

    def find(x: Any) -> Any:
        """Findet Wurzel der Komponente mit Pfadkompression."""
        if parent[x] != x:
            parent[x] = find(parent[x])  # Pfadkompression
        return parent[x]

    def union(x: Any, y: Any) -> bool:
        """
        Vereinigt die Komponenten von x und y.
        Returns True wenn sie vorher verschieden waren.
        """
        rx, ry = find(x), find(y)
        if rx == ry:
            return False  # Selbe Komponente → Zyklus!
        # Union by Rank: kleineren Baum unter größeren hängen
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1
        return True

    # Kanten greedy hinzufügen wenn kein Zyklus entsteht
    for u, v, w in edges_sorted:
        if union(u, v):
            mst.add_edge(u, v, w)

    return mst


def prim_mst(graph: Graph, start: int | str | None = None) -> Graph:
    """
    @brief Prim-Algorithmus: Minimaler Spannbaum.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Startet von einem Knoten und wächst den MST durch sukzessives Hinzufügen
    der billigsten Kante zur bereits aufgebauten Menge.

    Algorithmus: MinHeap über "Kreuzungskanten" (Kanten zwischen MST und Rest).
    Zeitkomplexität: O((V+E) log V)

    @param graph Der Graph (ungerichtet, gewichtet).
    @param start Startknoten (falls None: erster Knoten).
    @return Minimaler Spannbaum als neues Graph-Objekt.
    """
    mst = Graph(directed=False)
    if not graph.adj:
        return mst

    # Startknoten festlegen
    if start is None:
        start = graph.vertices()[0]

    # Alle Knoten in den MST übernehmen
    for v in graph.vertices():
        mst.add_vertex(v)

    in_mst = {start}         # Bereits im MST enthaltene Knoten
    # MinHeap: (kantengewicht, startknoten, endknoten)
    heap = []
    for v in graph.adj.get(start, set()):
        w = graph.weights.get((start, v), 1.0)
        heapq.heappush(heap, (w, start, v))

    while heap and len(in_mst) < len(graph.adj):
        w, u, v = heapq.heappop(heap)

        # Wenn v bereits im MST: Kante überspringen
        if v in in_mst:
            continue

        # Kante zum MST hinzufügen
        in_mst.add(v)
        mst.add_edge(u, v, w)

        # Neue Kreuzungskanten von v aus in den Heap
        for neighbor in graph.adj.get(v, set()):
            if neighbor not in in_mst:
                nw = graph.weights.get((v, neighbor), 1.0)
                heapq.heappush(heap, (nw, v, neighbor))

    return mst


def topological_sort(graph: Graph) -> list[int | str]:
    """
    @brief Topologische Sortierung für gerichtete azyklische Graphen (DAGs).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Kahn-Algorithmus: Iterative Entfernung von Knoten mit Eingangsgrad 0.

    Eine topologische Sortierung ordnet Knoten so an, dass alle Kanten
    von links nach rechts zeigen.

    Zeitkomplexität: O(V + E)

    @param graph Der gerichtete Graph (muss DAG sein).
    @return Liste der Knoten in topologischer Reihenfolge,
            oder [] wenn ein Zyklus vorhanden ist.
    """
    # Eingangsgrad für jeden Knoten berechnen
    in_degree = {v: 0 for v in graph.adj}
    for u in graph.adj:
        for v in graph.adj[u]:
            in_degree[v] = in_degree.get(v, 0) + 1

    # Alle Knoten mit Eingangsgrad 0 in die Warteschlange
    queue = sorted([v for v in in_degree if in_degree[v] == 0],
                   key=lambda x: (str(type(x)), x))
    result = []

    while queue:
        # Lexikographisch kleinsten Knoten wählen (für Determinismus)
        u = queue.pop(0)
        result.append(u)

        # Eingangsgrade der Nachbarn verringern
        for v in sorted(graph.adj.get(u, set()), key=lambda x: (str(type(x)), x)):
            in_degree[v] -= 1
            if in_degree[v] == 0:
                # Neuer Knoten mit Eingangsgrad 0: zur Warteschlange hinzufügen
                queue.append(v)
                queue.sort(key=lambda x: (str(type(x)), x))

    # Wenn nicht alle Knoten besucht: Zyklus vorhanden
    if len(result) != len(graph.adj):
        return []  # Kein gültiger DAG

    return result


# ============================================================
# Graphenanalyse und Färbung
# ============================================================

def chromatic_number_greedy(graph: Graph) -> int:
    """
    @brief Greedy-Algorithmus zur Bestimmung der chromatischen Zahl χ(G).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die chromatische Zahl χ(G) ist die minimale Anzahl Farben,
    die benötigt wird um G zu färben (keine zwei benachbarten
    Knoten gleiche Farbe).

    Dieser Greedy-Algorithmus liefert eine OBERE SCHRANKE.
    Das optimale χ(G) zu finden ist NP-schwer.

    @param graph Der Graph.
    @return Anzahl der verwendeten Farben (obere Schranke für χ(G)).
    """
    coloring = graph_coloring_greedy(graph)
    if not coloring:
        return 0
    return max(coloring.values()) + 1  # Farben: 0, 1, ..., k-1 → k Farben


def graph_coloring_greedy(graph: Graph) -> dict[int | str, int]:
    """
    @brief Greedy-Graph-Färbung.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Weist jedem Knoten die kleinste verfügbare Farbe zu,
    die keiner seiner Nachbarn hat.

    @param graph Der Graph.
    @return Dict {knoten: farbnummer} (Farben: 0, 1, 2, ...).
    """
    coloring = {}

    # Knoten nach Grad absteigend sortieren (häufig bessere Färbung)
    verts_sorted = sorted(graph.vertices(),
                          key=lambda v: -graph.degree(v))

    for v in verts_sorted:
        # Farben der Nachbarn sammeln
        neighbor_colors = {coloring[n] for n in graph.adj.get(v, set()) if n in coloring}

        # Kleinste verfügbare Farbe finden
        color = 0
        while color in neighbor_colors:
            color += 1

        coloring[v] = color

    return coloring


def is_eulerian(graph: Graph) -> dict[str, bool | list]:
    """
    @brief Prüft ob der Graph einen Euler-Kreis oder Euler-Pfad besitzt.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Euler-Kreis: Jede Kante genau einmal, Anfangs = Endknoten.
    Euler-Pfad: Jede Kante genau einmal, Anfangs ≠ Endknoten.

    Kriterien (für ungerichtete Graphen, nach Euler 1736):
    - Euler-Kreis ⟺ zusammenhängend + ALLE Knotengrade gerade
    - Euler-Pfad ⟺ zusammenhängend + GENAU 2 Knoten mit ungeradem Grad

    @param graph Der Graph.
    @return Dict mit 'euler_circuit', 'euler_path' und 'odd_degree_vertices'.
    """
    # Knoten mit ungeradem Grad finden
    odd_degree_verts = [v for v in graph.adj if graph.degree(v) % 2 == 1]
    connected = graph.is_connected()

    # Euler-Kreis: zusammenhängend + alle Grade gerade
    has_euler_circuit = connected and len(odd_degree_verts) == 0

    # Euler-Pfad: zusammenhängend + genau 2 Knoten mit ungeradem Grad
    has_euler_path = connected and len(odd_degree_verts) == 2

    return {
        'euler_circuit': has_euler_circuit,
        'euler_path': has_euler_path,
        'odd_degree_vertices': odd_degree_verts
    }


def is_hamiltonian_path(graph: Graph) -> dict[str, bool | list]:
    """
    @brief Sucht einen Hamiltonpfad (besucht alle Knoten genau einmal).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Das Hamiltonpfad-Problem ist NP-vollständig.
    Dieser Backtracking-Algorithmus hat exponentielle Laufzeit.

    WARNUNG: Nur für kleine Graphen (|V| <= 20) sinnvoll!

    @param graph Der Graph.
    @return Dict mit 'has_path' (bool) und 'path' (Liste oder []).
    """
    verts = graph.vertices()
    n = len(verts)

    if n == 0:
        return {'has_path': False, 'path': []}

    def backtrack(path: list, visited: set) -> list | None:
        """Backtracking-Suche nach einem Hamiltonpfad."""
        if len(path) == n:
            return path  # Alle Knoten besucht → Pfad gefunden!

        current = path[-1]
        for neighbor in sorted(graph.adj.get(current, set()), key=lambda x: (str(type(x)), x)):
            if neighbor not in visited:
                visited.add(neighbor)
                path.append(neighbor)
                result = backtrack(path, visited)
                if result is not None:
                    return result
                # Backtrack: letzten Schritt rückgängig machen
                path.pop()
                visited.discard(neighbor)

        return None  # Kein Pfad von diesem Zustand aus

    # Jeden Startknoten versuchen
    for start in verts:
        visited = {start}
        result = backtrack([start], visited)
        if result is not None:
            return {'has_path': True, 'path': result}

    return {'has_path': False, 'path': []}


def clique_number(graph: Graph) -> int:
    """
    @brief Kliquenzahl ω(G): Größte vollständige Teilmenge.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine Clique ist eine vollständige Teilmenge (alle Paare verbunden).
    Die Kliquenzahl ω(G) ist die maximale Cliquengröße.

    Algorithmus: Bron-Kerbosch (Backtracking).
    Zeitkomplexität: O(3^{n/3}) – exponentiell!
    Nur für kleine Graphen geeignet.

    @param graph Der Graph.
    @return Kliquenzahl ω(G).
    """
    max_clique_size = [0]  # Mutable für Closure

    def bron_kerbosch(R: set, P: set, X: set) -> None:
        """
        Bron-Kerbosch-Algorithmus zum Finden aller maximalen Cliquen.
        R: aktuelle Clique, P: Kandidaten, X: bereits verarbeitete Knoten.
        """
        if not P and not X:
            # Maximale Clique gefunden
            if len(R) > max_clique_size[0]:
                max_clique_size[0] = len(R)
            return

        # Pivot-Knoten wählen (Heuristik: maximaler Grad in P ∪ X)
        pivot = max(P | X, key=lambda v: len(graph.adj.get(v, set()) & P)) if (P | X) else None

        # Kandidaten = P ohne Nachbarn des Pivots
        for v in list(P - (graph.adj.get(pivot, set()) if pivot else set())):
            neighbors_v = graph.adj.get(v, set())
            bron_kerbosch(
                R | {v},
                P & neighbors_v,
                X & neighbors_v
            )
            P.discard(v)
            X.add(v)

    verts_set = set(graph.vertices())
    bron_kerbosch(set(), verts_set, set())
    return max_clique_size[0]


def independence_number(graph: Graph) -> int:
    """
    @brief Unabhängigkeitszahl α(G) = ω(Komplement(G)).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine unabhängige Menge enthält keine zwei benachbarten Knoten.
    Die Unabhängigkeitszahl ist die Kliquenzahl des Komplementgraphen.

    Relation: α(G) = ω(G̅)

    @param graph Der Graph.
    @return Unabhängigkeitszahl α(G).
    """
    # Kliquenzahl des Komplementgraphen berechnen
    complement = graph.complement()
    return clique_number(complement)


def graph_density(graph: Graph) -> float:
    """
    @brief Dichte des Graphen: d = |E| / (|V|·(|V|-1)/2).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Graphendichte gibt an, wie "vollständig" ein Graph ist.
    d = 0: leerer Graph, d = 1: vollständiger Graph.

    Mathematisch: d = |E| / C(|V|, 2) = 2|E| / (|V|·(|V|-1))

    @param graph Der Graph.
    @return Dichte zwischen 0.0 und 1.0.
    """
    n = len(graph.adj)
    if n < 2:
        return 0.0

    # Kanten zählen
    n_edges = len(graph.edges())

    # Maximale Kantenanzahl für ungerichteten Graphen
    max_edges = n * (n - 1) / 2

    return n_edges / max_edges


# ============================================================
# Bekannte Graphen-Konstruktoren
# ============================================================

def complete_graph(n: int) -> Graph:
    """
    @brief K_n: Vollständiger Graph mit n Knoten.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    In K_n ist jeder Knoten mit jedem anderen verbunden.
    |E(K_n)| = C(n,2) = n·(n-1)/2

    @param n Anzahl der Knoten (0, 1, 2, ..., n-1).
    @return Vollständiger Graph K_n.
    """
    g = Graph(directed=False)

    # Alle Knoten hinzufügen
    for i in range(n):
        g.add_vertex(i)

    # Alle möglichen Kanten hinzufügen
    for i in range(n):
        for j in range(i + 1, n):
            g.add_edge(i, j)

    return g


def cycle_graph(n: int) -> Graph:
    """
    @brief C_n: Kreisgraph mit n Knoten.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die n Knoten liegen auf einem Kreis: 0-1-2-...(n-1)-0.
    |E(C_n)| = n

    @param n Anzahl der Knoten (n ≥ 3).
    @return Kreisgraph C_n.
    """
    g = Graph(directed=False)

    for i in range(n):
        g.add_vertex(i)

    # Kanten im Kreis: 0→1, 1→2, ..., (n-2)→(n-1), (n-1)→0
    for i in range(n):
        g.add_edge(i, (i + 1) % n)

    return g


def path_graph(n: int) -> Graph:
    """
    @brief P_n: Pfadgraph mit n Knoten.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die n Knoten liegen auf einem Pfad: 0-1-2-...(n-1).
    |E(P_n)| = n-1

    @param n Anzahl der Knoten (n ≥ 1).
    @return Pfadgraph P_n.
    """
    g = Graph(directed=False)

    for i in range(n):
        g.add_vertex(i)

    # Kanten auf dem Pfad: 0→1, 1→2, ..., (n-2)→(n-1)
    for i in range(n - 1):
        g.add_edge(i, i + 1)

    return g


def bipartite_complete_graph(m: int, n: int) -> Graph:
    """
    @brief K_{m,n}: Vollständiger bipartiter Graph.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Knoten: Gruppe A = {a0, ..., a(m-1)}, Gruppe B = {b0, ..., b(n-1)}.
    Kanten: Alle Paare (ai, bj).
    |E(K_{m,n})| = m·n

    @param m Anzahl Knoten in Gruppe A.
    @param n Anzahl Knoten in Gruppe B.
    @return Vollständiger bipartiter Graph K_{m,n}.
    """
    g = Graph(directed=False)

    # Gruppe A: Knoten a0, a1, ...
    group_a = [f'a{i}' for i in range(m)]
    # Gruppe B: Knoten b0, b1, ...
    group_b = [f'b{j}' for j in range(n)]

    for v in group_a + group_b:
        g.add_vertex(v)

    # Alle Verbindungen zwischen A und B
    for a in group_a:
        for b in group_b:
            g.add_edge(a, b)

    return g


def petersen_graph() -> Graph:
    """
    @brief Der Petersen-Graph.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Petersen-Graph ist ein klassisches Gegenbeispiel in der Graphentheorie.
    Er ist:
    - 3-regulär (jeder Knoten hat Grad 3)
    - Nicht hamiltonsch (kein Hamiltonkreis)
    - Nicht planar (Kuratowski: enthält K_{3,3} als Minor)
    - 3-kanten-zusammenhängend

    10 Knoten, 15 Kanten.

    @return Der Petersen-Graph.
    """
    g = Graph(directed=False)

    for i in range(10):
        g.add_vertex(i)

    # Äußerer Fünfeck: 0-1-2-3-4-0
    outer = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    # Innerer Pentagramm: 5-7-9-6-8-5 (jeder Schritt überspringt einen Knoten)
    inner = [(5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]
    # Speichen: äußerer ↔ innerer Ring
    spokes = [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]

    for u, v in outer + inner + spokes:
        g.add_edge(u, v)

    return g


def grid_graph(m: int, n: int) -> Graph:
    """
    @brief m×n Gittergraph.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Knoten: (i, j) mit 0 ≤ i < m und 0 ≤ j < n.
    Kanten: Horizontal (i,j)-(i,j+1) und vertikal (i,j)-(i+1,j).

    @param m Anzahl Zeilen.
    @param n Anzahl Spalten.
    @return Gittergraph.
    """
    g = Graph(directed=False)

    # Alle Gitterpunkte anlegen
    for i in range(m):
        for j in range(n):
            g.add_vertex((i, j))

    # Horizontale Kanten
    for i in range(m):
        for j in range(n - 1):
            g.add_edge((i, j), (i, j + 1))

    # Vertikale Kanten
    for i in range(m - 1):
        for j in range(n):
            g.add_edge((i, j), (i + 1, j))

    return g


# ============================================================
# Kombinatorik
# ============================================================

def binomial_coefficient(n: int, k: int) -> int:
    """
    @brief Binomialkoeffizient C(n,k) = n! / (k! * (n-k)!).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Berechnet via dynamischer Programmierung (Pascalsches Dreieck)
    ohne große Zwischenergebnisse.

    Mathematisch:
        C(n, k) = \\frac{n!}{k! \\cdot (n-k)!}

    Eigenschaften:
    - C(n,0) = C(n,n) = 1
    - C(n,k) = C(n-1,k-1) + C(n-1,k)  (Pascalsches Dreieck)

    @param n Gesamtzahl.
    @param k Auswahl.
    @return Binomialkoeffizient als int.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1

    # Symmetrie ausnutzen: C(n,k) = C(n, n-k)
    k = min(k, n - k)

    # Iterative Berechnung (verhindert Überlauf gegenüber Fakultät)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)

    return result


def stirling_numbers_second_kind(n: int, k: int) -> int:
    """
    @brief Stirling-Zahlen zweiter Art S(n,k).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    S(n,k) ist die Anzahl der Möglichkeiten, eine n-elementige Menge
    in genau k nichtleere Teilmengen zu partitionieren.

    Rekursion:
        S(n, k) = k·S(n-1, k) + S(n-1, k-1)
    Randbedingungen:
        S(0, 0) = 1, S(n, 0) = 0 (n>0), S(0, k) = 0 (k>0)

    @param n Größe der Menge.
    @param k Anzahl der Teilmengen.
    @return Stirling-Zahl S(n,k).
    """
    if n == 0 and k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    if k > n:
        return 0

    # Dynamische Programmierung: dp[i][j] = S(i,j)
    dp = [[0] * (k + 1) for _ in range(n + 1)]
    dp[0][0] = 1

    for i in range(1, n + 1):
        for j in range(1, min(i, k) + 1):
            # Rekursionsformel: k·S(n-1,k) + S(n-1,k-1)
            dp[i][j] = j * dp[i - 1][j] + dp[i - 1][j - 1]

    return dp[n][k]


def bell_numbers(n: int) -> list[int]:
    """
    @brief Bell-Zahlen B_0, B_1, ..., B_n.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    B_n ist die Anzahl der Partitionen einer n-elementigen Menge.
    (Summe aller Stirling-Zahlen zweiter Art für festes n)

    Berechnung via Bell-Dreieck:
        B_0 = 1, B_{n+1} = Σ_{k=0}^{n} C(n,k) · B_k

    Erste Werte: B_0=1, B_1=1, B_2=2, B_3=5, B_4=15, B_5=52

    @param n Maximaler Index.
    @return Liste [B_0, B_1, ..., B_n].
    """
    if n < 0:
        return []

    # Bell-Dreieck: B[i] = Zahl am Anfang der i-ten Zeile
    B = [0] * (n + 1)
    B[0] = 1

    for i in range(1, n + 1):
        # B[i] = Σ_{k=0}^{i-1} C(i-1, k) * B[k]
        B[i] = sum(binomial_coefficient(i - 1, k) * B[k] for k in range(i))

    return B


def catalan_number(n: int) -> int:
    """
    @brief Catalan-Zahl C_n = C(2n,n) / (n+1).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Catalan-Zahlen zählen viele kombinatorische Strukturen:
    - Korrekt geklammerte Ausdrücke mit n Klammerpaaren
    - Binäre Suchbäume mit n+1 Blättern
    - Triangulierungen eines (n+2)-Ecks
    - Monotone Pfade unter der Diagonale in einem n×n-Gitter

    Formel: C_n = C(2n, n) / (n+1) = (2n)! / (n! · (n+1)!)

    Erste Werte: C_0=1, C_1=1, C_2=2, C_3=5, C_4=14, C_5=42

    @param n Index (n ≥ 0).
    @return Catalan-Zahl C_n.
    """
    if n < 0:
        return 0
    # C_n = C(2n, n) / (n+1)
    return binomial_coefficient(2 * n, n) // (n + 1)


def derangements(n: int) -> int:
    """
    @brief Anzahl der Permutationen ohne Fixpunkt (Subfakultät !n).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine Derangement ist eine Permutation σ mit σ(i) ≠ i für alle i.

    Formel via Inklusion-Exklusion:
        !n = n! · Σ_{k=0}^{n} (-1)^k / k!

    Rekursion:
        !0 = 1, !1 = 0
        !n = (n-1) · (!(n-1) + !(n-2))

    Asymptotisch: !n ≈ n!/e

    @param n Anzahl der Elemente.
    @return Anzahl der Derangements.
    """
    if n == 0:
        return 1
    if n == 1:
        return 0

    # Rekursive Formel via dynamischer Programmierung
    d_prev2 = 1  # !0
    d_prev1 = 0  # !1
    for i in range(2, n + 1):
        d_curr = (i - 1) * (d_prev1 + d_prev2)
        d_prev2 = d_prev1
        d_prev1 = d_curr

    return d_prev1


def partition_count(n: int) -> int:
    """
    @brief Anzahl der Ganzzahlpartitionen p(n).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    p(n) ist die Anzahl der Möglichkeiten, n als Summe positiver
    ganzer Zahlen zu schreiben (Reihenfolge egal).

    Beispiel: p(4) = 5: 4, 3+1, 2+2, 2+1+1, 1+1+1+1

    Algorithmus: Dynamische Programmierung mit Euler-Produktformel-Idee.
    dp[j] += dp[j-i] für i = 1 bis j.

    @param n Zu partitionierende Zahl.
    @return Anzahl der Partitionen p(n).
    """
    if n < 0:
        return 0
    if n == 0:
        return 1

    # dp[j] = Anzahl Partitionen von j
    dp = [0] * (n + 1)
    dp[0] = 1

    # Für jede Summanden-Größe i die Partitionen aktualisieren
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            dp[j] += dp[j - i]

    return dp[n]


def multinomial_coefficient(n: int, groups: list[int]) -> int:
    """
    @brief Multinomialkoeffizient n! / (k1! * k2! * ... * km!).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Verallgemeinert den Binomialkoeffizienten auf mehrere Gruppen.

    Anzahl der Anordnungen von n Elementen in Gruppen der Größen k1, k2, ...:

        M(n; k1, ..., km) = n! / (k1! · k2! · ... · km!)

    Voraussetzung: k1 + k2 + ... + km = n

    @param n Gesamtzahl.
    @param groups Liste der Gruppengrößen [k1, k2, ..., km].
    @return Multinomialkoeffizient als int.
    @raises ValueError wenn Summe der Gruppen ≠ n.
    """
    if sum(groups) != n:
        raise ValueError(f"Summe der Gruppengrößen {sum(groups)} muss gleich n={n} sein.")

    # Schrittweise Binomialkoeffizienten multiplizieren
    result = 1
    remaining = n
    for k in groups:
        result *= binomial_coefficient(remaining, k)
        remaining -= k

    return result
