# graph_theory.py – Graphentheorie und Kombinatorik

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Modul:** `src/graph_theory.py`

---

## Übersicht

Das Modul `graph_theory.py` implementiert Datenstrukturen und Algorithmen der Graphentheorie sowie kombinatorische Formeln:

- **Graph-Klasse** mit Adjazenzliste (gerichtet/ungerichtet, gewichtet)
- **Traversierungsalgorithmen:** BFS, DFS
- **Kürzeste Wege:** Dijkstra, Bellman-Ford, Floyd-Warshall
- **Minimale Spannbäume:** Kruskal (Union-Find), Prim (MinHeap)
- **Graphenanalyse:** Färbung, Euler-/Hamiltonpfade, Cliquen
- **Graphkonstruktoren:** $K_n$, $C_n$, $P_n$, $K_{m,n}$, Petersen, Gitter
- **Kombinatorik:** Binomialkoeffizient, Catalan, Bell, Stirling, Derangements

---

## Mathematischer Hintergrund

### Grundbegriffe

Ein **Graph** $G = (V, E)$ besteht aus:
- Knotenmenge $V$ (vertices)
- Kantenmenge $E \subseteq V \times V$ (edges)

**Wichtige Varianten:**

| Typ | Eigenschaft |
|-----|-------------|
| Ungerichtet | $(u,v) \in E \iff (v,u) \in E$ |
| Gerichtet (Digraph) | $(u,v) \in E$ impliziert nicht $(v,u) \in E$ |
| Gewichtet | Jede Kante $(u,v)$ hat Gewicht $w(u,v) \in \mathbb{R}$ |
| DAG | Directed Acyclic Graph – gerichtet und zyklenfrei |

**Grundgrößen:**

$$|E(K_n)| = \binom{n}{2} = \frac{n(n-1)}{2} \qquad \text{(vollständiger Graph)}$$

$$\text{Dichte: } d(G) = \frac{|E|}{\binom{|V|}{2}} \in [0, 1]$$

---

## Klasse `Graph`

### Interne Darstellung

```
adj: { vertex: set(neighbors) }
weights: { (u,v): float }
```

Die Adjazenzliste bietet:
- $O(1)$ Knotengrad-Abfrage
- $O(\deg(v))$ Nachbarschaftsabfrage
- $O(|V| + |E|)$ Speicher

### API-Übersicht

| Methode | Komplexität | Beschreibung |
|---------|------------|-------------|
| `add_vertex(v)` | $O(1)$ | Knoten hinzufügen |
| `add_edge(u,v,w)` | $O(1)$ | Kante hinzufügen |
| `remove_vertex(v)` | $O(\deg(v))$ | Knoten entfernen |
| `remove_edge(u,v)` | $O(1)$ | Kante entfernen |
| `degree(v)` | $O(1)$ | Knotengrad |
| `is_connected()` | $O(V+E)$ | Zusammenhang (BFS) |
| `connected_components()` | $O(V+E)$ | Alle Komponenten |
| `has_cycle()` | $O(V+E)$ | Zykluserkennung (DFS) |
| `is_tree()` | $O(V+E)$ | Baum-Test |
| `is_bipartite()` | $O(V+E)$ | 2-Färbbarkeit (BFS) |
| `adjacency_matrix()` | $O(V^2)$ | Als 2D-Liste |
| `complement()` | $O(V^2)$ | Komplementgraph |
| `subgraph(verts)` | $O(V+E)$ | Induzierter Teilgraph |

---

## Traversierungsalgorithmen

### Breitensuche (BFS)

**Idee:** Erkunde alle Knoten schichtweise nach Entfernung vom Start.

**Algorithmus:**
1. Start: Distanz $d[s] = 0$, alle anderen $\infty$
2. Warteschlange (FIFO): $Q = [s]$
3. Für jeden Knoten $u$ aus $Q$: relaxiere alle Nachbarn

```
d[v] := d[u] + 1 falls d[v] = ∞
```

**Zeitkomplexität:** $O(V + E)$

**Anwendung:** Kürzeste Wege in ungewichteten Graphen, Level-Order-Traversal.

### Tiefensuche (DFS)

**Idee:** Erkunde so tief wie möglich, bevor zurückgegangen wird.

**Zeitkomplexität:** $O(V + E)$

**Anwendung:** Zykluserkennung, topologische Sortierung, Zusammenhangskomponenten, Brücken/Artikulationspunkte.

---

## Kürzeste-Wege-Algorithmen

### Dijkstra

**Voraussetzung:** Nicht-negative Kantengewichte.

**Idee (Greedy):** Immer den Knoten mit kleinstem bekannten Abstand wählen und relaxieren.

$$d[v] = \min(d[v], \, d[u] + w(u,v))$$

**Zeitkomplexität:** $O((V+E) \log V)$ mit Binary MinHeap.

**Korrektheit:** Greedy-Choice-Property: Einmal finalisierter Abstand ist optimal.

### Bellman-Ford

**Voraussetzung:** Beliebige Kantengewichte (erkennt negative Zyklen).

**Idee:** $|V|-1$ Relaxierungsrunden über alle Kanten.

**Zeitkomplexität:** $O(V \cdot E)$

**Negativer Zyklus:** Falls nach $|V|-1$ Runden noch Relaxierungen möglich → negativer Zyklus vorhanden.

### Floyd-Warshall

**Voraussetzung:** Beliebige Gewichte, kein negativer Zyklus.

**Idee (DP):** $d[i][j] = \min(d[i][j], \, d[i][k] + d[k][j])$ für alle Zwischenknoten $k$.

**Zeitkomplexität:** $O(V^3)$

**Anwendung:** Alle kürzesten Wege in einem Schritt.

---

## Minimale Spannbäume

Ein **Minimaler Spannbaum** (MST) ist ein Teilgraph, der:
- alle Knoten enthält,
- zusammenhängend ist,
- azyklisch ist (Baum),
- das minimale Gesamtgewicht hat.

**MST-Eigenschaften:**
- $|V|-1$ Kanten
- Eindeutig wenn alle Gewichte verschieden

### Kruskal (Union-Find)

```
Sortiere Kanten aufsteigend nach Gewicht
Für jede Kante (u,v,w):
    Falls find(u) ≠ find(v): Kante hinzufügen, union(u,v)
```

**Zeitkomplexität:** $O(E \log E)$
**Korrektheit:** Schnitteigenschaft (Cut Property).

### Prim (MinHeap)

```
Starte von Knoten s, in_MST = {s}
Heap: alle Kanten von s
Solange Heap nicht leer:
    (w, u, v) = pop_min(Heap)
    Falls v ∉ in_MST: Kante hinzufügen, neue Kanten von v pushen
```

**Zeitkomplexität:** $O((V+E) \log V)$

---

## Topologische Sortierung (Kahn)

Für **DAGs**: Ordnet Knoten so, dass alle Kanten von links nach rechts zeigen.

**Algorithmus:**
1. Berechne Eingangsgrade $\deg^-(v)$
2. Initialisiere Queue mit allen $v$ mit $\deg^-(v) = 0$
3. Füge $v$ zur Reihenfolge, verringere $\deg^-$ der Nachbarn
4. Wenn $|V|$ Knoten verarbeitet: DAG; sonst: Zyklus vorhanden

**Zeitkomplexität:** $O(V + E)$

---

## Graphenanalyse

### Graphenfärbung

Die **chromatische Zahl** $\chi(G)$ ist die minimale Anzahl Farben für eine gültige Färbung (keine zwei benachbarten Knoten gleiche Farbe).

**Greedy-Obergrenze:** $\chi(G) \leq \Delta(G) + 1$ (Brook's Theorem: $\chi(G) \leq \Delta(G)$ für die meisten Graphen).

| Graph | $\chi(G)$ |
|-------|----------|
| Bipartit | 2 |
| Zyklus $C_n$ (gerade) | 2 |
| Zyklus $C_n$ (ungerade) | 3 |
| $K_n$ | $n$ |
| Planar | $\leq 4$ (Vierfarbensatz) |

### Euler-Wege und -Kreise

**Satz von Euler (1736):**
- **Euler-Kreis** (geschlossen) $\iff$ zusammenhängend + **alle** Knotengrade gerade
- **Euler-Pfad** (offen) $\iff$ zusammenhängend + genau **2** Knoten mit ungeradem Grad

### Hamilton-Problem

**Hamiltonpfad:** Besucht jeden Knoten **genau einmal**.
**NP-vollständig** – kein effizienter Algorithmus bekannt.

Implementiert via **Backtracking** (exponentielle Laufzeit, nur für kleine Graphen).

### Cliquen und Unabhängige Mengen

$$\omega(G) = \text{Kliquenzahl (größte vollständige Teilmenge)}$$
$$\alpha(G) = \text{Unabhängigkeitszahl} = \omega(\bar{G})$$

**Relation:** $\alpha(G) \cdot \omega(G) \geq n$ (Ramsey-Theorie)

---

## Bekannte Graphen

| Funktion | Beschreibung | Eigenschaften |
|----------|-------------|---------------|
| `complete_graph(n)` | $K_n$ | $\binom{n}{2}$ Kanten, $\chi = n$ |
| `cycle_graph(n)` | $C_n$ | $n$ Kanten, 2-regulär |
| `path_graph(n)` | $P_n$ | $n-1$ Kanten, Baum |
| `bipartite_complete_graph(m,n)` | $K_{m,n}$ | $mn$ Kanten, $\chi = 2$ |
| `petersen_graph()` | Petersen | 10 Knoten, 15 Kanten, 3-regulär, nicht hamiltonsch |
| `grid_graph(m,n)` | $m \times n$ Gitter | $2mn - m - n$ Kanten, bipartit |

### Petersen-Graph

Der Petersen-Graph ist ein wichtiges Gegenbeispiel:
- 3-regulär, nicht planar
- Hat **keinen** Hamiltonkreis (aber Hamiltonpfad)
- Chromatische Zahl: 3
- Enthält $K_{3,3}$ als Minoren (Kuratowski: nicht planar)

---

## Kombinatorik

### Binomialkoeffizient

$$\binom{n}{k} = \frac{n!}{k!(n-k)!}$$

**Pascalsches Dreieck:** $\binom{n}{k} = \binom{n-1}{k-1} + \binom{n-1}{k}$

Implementiert iterativ ohne Fakultät (verhindert Überlauf).

### Stirling-Zahlen 2. Art

$S(n,k)$: Anzahl der Partitionen einer $n$-Menge in $k$ nichtleere Teilmengen.

**Rekursion:** $S(n,k) = k \cdot S(n-1,k) + S(n-1,k-1)$

| $S(n,k)$ | $k=1$ | $k=2$ | $k=3$ | $k=4$ |
|----------|-------|-------|-------|-------|
| $n=1$ | 1 | | | |
| $n=2$ | 1 | 1 | | |
| $n=3$ | 1 | 3 | 1 | |
| $n=4$ | 1 | 7 | 6 | 1 |

### Bell-Zahlen

$B_n = \sum_{k=0}^{n} S(n,k)$ = Anzahl aller Partitionen einer $n$-Menge.

$B_0=1, B_1=1, B_2=2, B_3=5, B_4=15, B_5=52$

### Catalan-Zahlen

$$C_n = \frac{1}{n+1}\binom{2n}{n}$$

$C_0=1, C_1=1, C_2=2, C_3=5, C_4=14, C_5=42$

**Zählt u.a.:**
- Korrekt geklammerte Ausdrücke mit $n$ Klammerpaaren
- Binäre Bäume mit $n+1$ Blättern
- Monotone Gitterpfade ohne Diagonalüberschreitung

### Derangements (Subfakultät)

$!n$: Permutationen von $n$ Elementen ohne Fixpunkt.

**Rekursion:** $!n = (n-1)(!(n-1) + !(n-2))$

**Asymptotisch:** $!n \approx n!/e$

$!0=1, !1=0, !2=1, !3=2, !4=9, !5=44$

### Ganzzahlpartitionen

$p(n)$: Anzahl der Darstellungen von $n$ als Summe positiver ganzer Zahlen.

$p(1)=1, p(2)=2, p(3)=3, p(4)=5, p(5)=7, p(10)=42$

---

## Verwendungsbeispiele

```python
from src.graph_theory import (
    Graph, complete_graph, dijkstra, kruskal_mst,
    is_eulerian, catalan_number, bell_numbers
)

# Graph erstellen und analysieren
g = Graph(directed=False)
g.add_edge(0, 1, weight=3.0)
g.add_edge(1, 2, weight=1.0)
g.add_edge(0, 2, weight=5.0)

# Kürzeste Wege (Dijkstra)
result = dijkstra(g, 0)
print(result['distances'])  # {0: 0, 1: 3.0, 2: 4.0}
print(result['path'][2])    # [0, 1, 2]

# Minimaler Spannbaum
mst = kruskal_mst(g)
print(f"MST Kanten: {mst.edges()}")

# Euler-Eigenschaften
g_cycle = complete_graph(5)  # K_5
euler = is_eulerian(g_cycle)
print(f"Euler-Kreis: {euler['euler_circuit']}")   # True (alle Grad 4, gerade)

# Kombinatorik
print(catalan_number(5))    # 42
print(bell_numbers(5))      # [1, 1, 2, 5, 15, 52]
```

---

## Komplexitätsübersicht

| Algorithmus | Zeit | Raum |
|------------|------|------|
| BFS | $O(V+E)$ | $O(V)$ |
| DFS | $O(V+E)$ | $O(V)$ |
| Dijkstra | $O((V+E)\log V)$ | $O(V)$ |
| Bellman-Ford | $O(V \cdot E)$ | $O(V)$ |
| Floyd-Warshall | $O(V^3)$ | $O(V^2)$ |
| Kruskal-MST | $O(E \log E)$ | $O(V)$ |
| Prim-MST | $O((V+E)\log V)$ | $O(V)$ |
| Topol. Sortierung | $O(V+E)$ | $O(V)$ |
| Greedy-Färbung | $O(V+E)$ | $O(V)$ |
| Clique (Bron-Kerbosch) | $O(3^{n/3})$ | $O(V)$ |
| Hamiltonpfad | $O(n!)$ | $O(V)$ |

---

## Tests

Die Tests liegen in `/tests/test_graph_theory.py`:

- `TestGraphBasicOperations`: add/remove vertex/edge, degree, Adjazenzmatrix
- `TestGraphProperties`: Zusammenhang, Zyklen, Bipartitheit, Komplement
- `TestBFS`: Abstände, Besuchsreihenfolge, Vorgänger
- `TestDFS`: Entdeckungs-/Abschlusszeiten
- `TestDijkstra`: Gewichtete Wege, Pfadrekonstruktion, unerreichbare Knoten
- `TestBellmanFord`: Negative Gewichte, Zykluserkennung
- `TestFloydWarshall`: Alle-Paare-Distanzen
- `TestSpanningTrees`: MST-Korrektheit (Kruskal + Prim)
- `TestTopologicalSort`: DAG-Sortierung, Zyklus-Erkennung
- `TestGraphAnalysis`: Färbung, Euler, Hamilton, Cliquen
- `TestGraphConstructors`: $K_n$, $C_n$, $P_n$, Petersen, Gitter
- `TestCombinatorics`: Binomial, Stirling, Bell, Catalan, Derangements, Partitionen
- `TestEdgeCases`: Leerer Graph, Einzelknoten, Schleifen, große Graphen

---

## Referenzen

- Diestel, R. (2017). *Graph Theory* (5th ed.). Springer.
- Cormen, T. H. et al. (2022). *Introduction to Algorithms* (4th ed.). MIT Press.
- West, D. B. (2001). *Introduction to Graph Theory* (2nd ed.). Prentice Hall.
- Graham, R. L., Knuth, D. E., Patashnik, O. (1994). *Concrete Mathematics*. Addison-Wesley.
