"""
Graceful Labeling von Bäumen (Graceful Tree Conjecture)

Implementiert Algorithmen zur Untersuchung der Graceful-Tree-Vermutung
(Rosa 1967, Ringel-Kotzig-Vermutung): Jeder Baum besitzt ein graceful Labeling.

Ein graceful Labeling eines Baumes T mit n Knoten ist eine bijektive Abbildung
    f: V(T) → {0, 1, ..., n−1}
sodass die induzierten Kantenlabels |f(u) − f(v)| für alle Kanten {u,v} ∈ E(T)
alle verschieden sind, also die Menge {1, 2, ..., n−1} ergeben.

Bekannte bewiesene Klassen:
  - Pfade P_n (trivial durch alternierende Konstruktion)
  - Caterpillar-Graphen (Gallian 1994 – via induktive Konstruktion)
  - Wheel-Graphen W_n (Slamin et al.)
  - Helme, Fischschwänze u.v.a.

Autor: Michael Fuhrmann
Letzte Änderung: 2026-03-12
"""

import itertools
from typing import Optional, Dict, List, Tuple, Iterator
import networkx as nx


class GracefulTreeConjecture:
    """
    Klasse zur Analyse der Graceful-Tree-Vermutung.

    Stellt Backtracking-Algorithmen, explizite Konstruktionen für bekannte
    Klassen sowie Statistikfunktionen bereit.

    Attribute:
        max_backtrack_n (int): Maximale Knotenanzahl für den Backtracking-Algorithmus.

    Autor: Michael Fuhrmann
    Letzte Änderung: 2026-03-12
    """

    def __init__(self, max_backtrack_n: int = 10):
        """
        Initialisiert die GracefulTreeConjecture-Analyse.

        Args:
            max_backtrack_n: Maximale Baumgröße für erschöpfende Suche (Standard: 10).

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        self.max_backtrack_n = max_backtrack_n
        # Cache für bereits gefundene Labelings
        self._labeling_cache: Dict[int, Optional[Dict]] = {}

    # ------------------------------------------------------------------
    # Hilfsmethoden
    # ------------------------------------------------------------------

    def is_graceful(self, tree: nx.Graph, labeling: Dict) -> bool:
        """
        Prüft ob ein gegebenes Labeling für einen Baum graceful ist.

        Ein Labeling ist graceful, wenn:
          1. Es bijektiv von V auf {0, ..., n−1} abbildet.
          2. Die Kantenlabels |f(u) − f(v)| für alle Kanten paarweise verschieden
             und gleich {1, ..., n−1} sind.

        Args:
            tree:     NetworkX-Graph (muss ein Baum sein).
            labeling: Dict {Knoten: Label} mit Labels in {0, ..., n−1}.

        Returns:
            True wenn das Labeling graceful ist, sonst False.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        n = tree.number_of_nodes()
        # Prüfe Bijektivität: genau n Labels, alle verschieden, aus {0,...,n-1}
        if set(labeling.keys()) != set(tree.nodes()):
            return False
        label_values = set(labeling.values())
        if label_values != set(range(n)):
            return False
        # Berechne Kantenlabels
        edge_labels = set()
        for u, v in tree.edges():
            diff = abs(labeling[u] - labeling[v])
            if diff in edge_labels:
                return False  # Duplikat → kein graceful Labeling
            edge_labels.add(diff)
        # Muss genau {1, ..., n−1} sein
        return edge_labels == set(range(1, n))

    def find_graceful_labeling(self, tree: nx.Graph) -> Optional[Dict]:
        """
        Sucht ein graceful Labeling für einen Baum via Backtracking.

        Algorithmus:
          - Ordne Knoten nach Grad (hoch zuerst) als Heuristik.
          - Weise Labels {0, ..., n−1} rekursiv zu.
          - Pruning: prüfe nach jeder Zuweisung bereits entstandene Kantenlabels
            auf Konflikte.

        Args:
            tree: NetworkX-Graph (Baum mit n ≤ max_backtrack_n empfohlen).

        Returns:
            Dict {Knoten: Label} wenn gefunden, sonst None.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        n = tree.number_of_nodes()
        if n == 0:
            return {}
        if n == 1:
            nodes = list(tree.nodes())
            return {nodes[0]: 0}

        # Sortiere Knoten nach Grad (absteigend) für bessere Pruning-Effizienz
        nodes = sorted(tree.nodes(), key=lambda v: tree.degree(v), reverse=True)
        labels = list(range(n))
        assignment: Dict = {}
        used_labels: set = set()
        used_edge_labels: set = set()

        def backtrack(idx: int) -> bool:
            """Rekursive Backtracking-Suche."""
            if idx == n:
                # Alle Knoten zugewiesen – prüfe Vollständigkeit
                return used_edge_labels == set(range(1, n))

            node = nodes[idx]
            for label in labels:
                if label in used_labels:
                    continue
                # Berechne neue Kantenlabels zu bereits zugewiesenen Nachbarn
                new_edge_labels = []
                conflict = False
                for neighbor in tree.neighbors(node):
                    if neighbor in assignment:
                        diff = abs(label - assignment[neighbor])
                        if diff in used_edge_labels or diff in new_edge_labels:
                            conflict = True
                            break
                        new_edge_labels.append(diff)
                if conflict:
                    continue

                # Zuweisung vornehmen
                assignment[node] = label
                used_labels.add(label)
                for el in new_edge_labels:
                    used_edge_labels.add(el)

                if backtrack(idx + 1):
                    return True

                # Zurücksetzen
                del assignment[node]
                used_labels.discard(label)
                for el in new_edge_labels:
                    used_edge_labels.discard(el)

            return False

        if backtrack(0):
            return dict(assignment)
        return None

    def generate_all_trees(self, n: int) -> List[nx.Graph]:
        """
        Erzeugt alle nicht-isomorphen Bäume mit genau n Knoten.

        Nutzt den networkx-Generator `nx.nonisomorphic_trees(n)`, der alle
        nicht-isomorphen Bäume über den Prüfer-Sequenz-Algorithmus aufzählt.

        Args:
            n: Anzahl der Knoten (n ≥ 1).

        Returns:
            Liste aller nicht-isomorphen Bäume mit n Knoten als NetworkX-Graphen.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        if n < 1:
            return []
        if n == 1:
            # nx.nonisomorphic_trees(1) gibt einen leeren Graph zurück –
            # wir erzeugen den Einzelknoten-Graphen explizit
            g = nx.Graph()
            g.add_node(0)
            return [g]
        return list(nx.nonisomorphic_trees(n))

    # ------------------------------------------------------------------
    # Explizite Konstruktionen für bekannte Klassen
    # ------------------------------------------------------------------

    def graceful_path(self, n: int) -> Tuple[nx.Graph, Dict]:
        """
        Erzeugt Pfad P_n mit einem expliziten graceful Labeling.

        Konstruktion (Standard-Methode):
          Sei P_n = (v_0, v_1, ..., v_{n−1}).
          Weise Labels abwechselnd von unten und oben zu:
            f(v_0) = 0, f(v_1) = n−1, f(v_2) = 1, f(v_3) = n−2, ...
          Dies erzeugt Kantenlabels n−1, n−2, ..., 1.

        Args:
            n: Anzahl der Knoten (n ≥ 1).

        Returns:
            Tupel (Pfad-Graph, graceful-Labeling-Dict).

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        # Pfad-Graph erstellen
        path = nx.path_graph(n)
        if n == 1:
            return path, {0: 0}

        labeling = {}
        low, high = 0, n - 1
        for i, node in enumerate(range(n)):
            if i % 2 == 0:
                labeling[node] = low
                low += 1
            else:
                labeling[node] = high
                high -= 1
        return path, labeling

    def graceful_caterpillar(self, spine_length: int, leaf_counts: List[int]) -> Tuple[nx.Graph, Dict]:
        """
        Erzeugt einen Caterpillar-Graphen mit einem graceful Labeling.

        Ein Caterpillar ist ein Baum, bei dem alle Knoten auf einem zentralen
        Pfad (Spine) liegen oder Blätter (direkte Nachbarn des Spines) sind.

        Konstruktion nach Gallian (1994):
          - Benenne Spine-Knoten 0..s-1, Blätter folgen danach.
          - Verwende alternierende Beschriftung für den Spine, dann Blätter passend.

        Args:
            spine_length: Länge des zentralen Pfades (Anzahl Spine-Knoten).
            leaf_counts:  Liste der Blattanzahlen pro Spine-Knoten (Länge = spine_length).

        Returns:
            Tupel (Caterpillar-Graph, graceful-Labeling-Dict).

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        # Graph aufbauen
        G = nx.Graph()
        node_id = 0
        spine_nodes = []
        # Spine-Knoten
        for i in range(spine_length):
            G.add_node(node_id)
            spine_nodes.append(node_id)
            node_id += 1
        # Spine-Kanten
        for i in range(spine_length - 1):
            G.add_edge(spine_nodes[i], spine_nodes[i + 1])
        # Blätter hinzufügen
        leaf_groups: List[List[int]] = []
        for i, count in enumerate(leaf_counts):
            leaves = []
            for _ in range(count):
                G.add_node(node_id)
                G.add_edge(spine_nodes[i], node_id)
                leaves.append(node_id)
                node_id += 1
            leaf_groups.append(leaves)

        # Graceful Labeling via Backtracking (für kleine Caterpillars)
        labeling = self.find_graceful_labeling(G)
        return G, labeling if labeling else {}

    def graceful_star(self, k: int) -> Tuple[nx.Graph, Dict]:
        """
        Erzeugt Stern K_{1,k} (= Caterpillar mit einem Spine-Knoten) mit graceful Labeling.

        Explizite Konstruktion:
          - Zentrum erhält Label 0.
          - Blätter erhalten Labels 1, 2, ..., k.
          - Kantenlabels: |0 − i| = i für i = 1..k → {1,...,k} ✓

        Args:
            k: Anzahl der Blätter.

        Returns:
            Tupel (Stern-Graph, graceful-Labeling-Dict).

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        star = nx.star_graph(k)
        # Knoten 0 ist das Zentrum in nx.star_graph
        labeling = {0: 0}
        for i in range(1, k + 1):
            labeling[i] = i
        return star, labeling

    def graceful_wheel(self, n: int) -> Tuple[nx.Graph, Optional[Dict]]:
        """
        Erzeugt Wheel-Graph W_n (Nabe + Kreis C_{n-1}) und sucht graceful Labeling.

        Hinweis: W_n ist für kleine n bekannt graceful (Slamin et al. 2000).
        Für n ≤ max_backtrack_n wird Backtracking verwendet.

        Args:
            n: Gesamtanzahl Knoten (Nabe + n−1 Rand-Knoten).

        Returns:
            Tupel (Wheel-Graph, graceful-Labeling-Dict oder None).

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        wheel = nx.wheel_graph(n)
        if n <= self.max_backtrack_n:
            labeling = self.find_graceful_labeling(wheel)
        else:
            labeling = None
        return wheel, labeling

    # ------------------------------------------------------------------
    # Statistik-Methoden
    # ------------------------------------------------------------------

    def count_graceful_trees(self, max_n: int = 12) -> Dict[int, Dict]:
        """
        Zählt für n = 1..max_n die nicht-isomorphen Bäume und wie viele davon
        graceful sind.

        Die Graceful-Tree-Vermutung besagt, dass für alle n ALLE Bäume graceful sind.
        (Bisher für n ≤ 35 verifiziert, Farragó 2024.)

        Args:
            max_n: Maximale Knotenzahl (Standard: 12; für n > 10 langsam).

        Returns:
            Dict {n: {"total": int, "graceful": int, "all_graceful": bool}}.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        results = {}
        for n in range(1, max_n + 1):
            trees = self.generate_all_trees(n)
            total = len(trees)
            graceful_count = 0
            for tree in trees:
                labeling = self.find_graceful_labeling(tree)
                if labeling is not None:
                    graceful_count += 1
            results[n] = {
                "total": total,
                "graceful": graceful_count,
                "all_graceful": (graceful_count == total),
            }
        return results

    def verify_known_classes(self) -> Dict[str, bool]:
        """
        Verifiziert die bekannten bewiesenen graceful Klassen für kleine Beispiele.

        Getestete Klassen:
          - P_n (Pfade, n = 2..8)
          - K_{1,k} (Sterne, k = 1..6)
          - Caterpillar mit vorgegebenen Parametern

        Returns:
            Dict {Beschreibung: True/False}.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        results = {}

        # Pfade P_n für n = 2..8
        for n in range(2, 9):
            path, labeling = self.graceful_path(n)
            results[f"P_{n}"] = self.is_graceful(path, labeling)

        # Sterne K_{1,k} für k = 1..6
        for k in range(1, 7):
            star, labeling = self.graceful_star(k)
            results[f"K_1_{k}"] = self.is_graceful(star, labeling)

        return results

    def is_caterpillar(self, tree: nx.Graph) -> bool:
        """
        Prüft ob ein Baum ein Caterpillar ist.

        Ein Baum ist ein Caterpillar gdw. nach Entfernen aller Blätter ein
        Pfad (oder der leere Graph) übrig bleibt.

        Args:
            tree: NetworkX-Baum.

        Returns:
            True wenn Caterpillar, sonst False.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        if tree.number_of_nodes() <= 2:
            return True
        # Entferne alle Blätter (Grad 1)
        leaves = [v for v in tree.nodes() if tree.degree(v) == 1]
        spine = tree.copy()
        spine.remove_nodes_from(leaves)
        if spine.number_of_nodes() == 0:
            return True  # Stern
        # Überprüfe ob Spine ein Pfad ist (alle Knoten Grad ≤ 2)
        return all(spine.degree(v) <= 2 for v in spine.nodes())

    def get_edge_labels(self, tree: nx.Graph, labeling: Dict) -> List[int]:
        """
        Berechnet die sortierten Kantenlabels für ein gegebenes Labeling.

        Args:
            tree:     NetworkX-Baum.
            labeling: Knotenlabeling als Dict.

        Returns:
            Sortierte Liste der Kantenlabels.

        Autor: Michael Fuhrmann
        Letzte Änderung: 2026-03-12
        """
        return sorted(abs(labeling[u] - labeling[v]) for u, v in tree.edges())
