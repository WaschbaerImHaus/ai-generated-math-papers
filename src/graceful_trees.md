# graceful_trees.py – Dokumentation

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `src/graceful_trees.py`

---

## Überblick

Das Modul untersucht die **Graceful-Tree-Vermutung** (auch *Ringel-Kotzig-Vermutung*), eine der bekanntesten offenen Vermutungen in der kombinatorischen Graphentheorie.

---

## Mathematischer Hintergrund

### Definition: Graceful Labeling

Sei $T$ ein Baum mit $n$ Knoten und $n-1$ Kanten. Ein **graceful Labeling** ist eine injektive Abbildung

$$f: V(T) \to \{0, 1, \ldots, n-1\}$$

sodass die induzierten **Kantenlabels**

$$|f(u) - f(v)| \quad \text{für alle Kanten } \{u,v\} \in E(T)$$

paarweise verschieden sind und genau die Menge $\{1, 2, \ldots, n-1\}$ ergeben.

### Die Graceful-Tree-Vermutung

> **Vermutung (Rosa 1967 / Ringel-Kotzig):**
> *Jeder Baum besitzt ein graceful Labeling.*

**Status**: Offen. Verifiziert für alle Bäume mit bis zu 35 Knoten (Farragó 2024).

### Bekannte bewiesene Klassen

| Klasse | Beweis | Konstruktion |
|--------|--------|--------------|
| Pfade $P_n$ | Trivial | Alternierende Zuweisung |
| Sterne $K_{1,k}$ | Trivial | Zentrum = 0, Blätter = 1..k |
| Caterpillar-Graphen | Gallian (1994) | Induktiv über den Spine |
| Lobster-Graphen | Teils bekannt | Komplex |
| Helmets, Fans, ... | Verschiedene | Explizite Formeln |

---

## Klasse: `GracefulTreeConjecture`

### Konstruktor

```python
GracefulTreeConjecture(max_backtrack_n: int = 10)
```

- `max_backtrack_n`: Maximale Knotenanzahl für Backtracking-Suche.

---

### Methoden

#### `is_graceful(tree, labeling) → bool`

Prüft ob ein gegebenes Labeling graceful ist.

**Algorithmus:**
1. Bijektivitätscheck: Labels $\subseteq \{0,\ldots,n-1\}$, keine Duplikate.
2. Kantenlabels berechnen: $\{|f(u)-f(v)| : \{u,v\} \in E\}$.
3. Gleichheit mit $\{1,\ldots,n-1\}$ prüfen.

**Komplexität:** $O(|E|)$ nach dem Bijektivitätscheck.

---

#### `find_graceful_labeling(tree) → Optional[Dict]`

Backtracking-Suche nach einem graceful Labeling.

**Heuristik:** Knoten werden nach absteigendem Grad sortiert (Knoten hohen Grades schränken früher ein → besseres Pruning).

**Pruning-Strategie:**
- Nach jeder Zuweisung werden Kantenlabels zu bereits beschrifteten Nachbarn geprüft.
- Konflikte (doppelte Kantenlabels) werden sofort erkannt → frühes Backtrack.

**Komplexität:** Worst case $O(n!)$, praktisch durch Pruning deutlich besser.

---

#### `generate_all_trees(n) → List[nx.Graph]`

Erzeugt alle nicht-isomorphen Bäume mit $n$ Knoten über den NetworkX-Generator `nonisomorphic_trees(n)`.

Die Anzahl nicht-isomorpher Bäume folgt der **OEIS-Sequenz A000055**:

| $n$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|-----|---|---|---|---|---|---|---|---|---|----|
| Anzahl | 1 | 1 | 1 | 2 | 3 | 6 | 11 | 23 | 47 | 106 |

---

#### `graceful_path(n) → Tuple[Graph, Dict]`

Explizites graceful Labeling für Pfad $P_n$:

$$f(v_i) = \begin{cases} i/2 & i \text{ gerade} \\ n-1-\lfloor i/2 \rfloor & i \text{ ungerade} \end{cases}$$

Beispiel für $P_5$: Labels $(0, 4, 1, 3, 2)$, Kanten: $4, 3, 2, 1$.

---

#### `graceful_star(k) → Tuple[Graph, Dict]`

Explizites graceful Labeling für $K_{1,k}$:

$$f(\text{Zentrum}) = 0, \quad f(\text{Blatt}_i) = i \quad (i=1,\ldots,k)$$

Kantenlabels: $\{1, 2, \ldots, k\}$ ✓

---

#### `graceful_wheel(n) → Tuple[Graph, Optional[Dict]]`

Sucht graceful Labeling für Wheel-Graph $W_n$ via Backtracking. Für $n \leq 10$ praktikabel.

---

#### `is_caterpillar(tree) → bool`

Prüft ob ein Baum ein **Caterpillar** ist:

> Ein Baum ist ein Caterpillar genau dann, wenn nach Entfernung aller Blätter (Knoten mit Grad 1) ein Pfad (oder der leere Graph) übrig bleibt.

---

#### `count_graceful_trees(max_n) → Dict`

Zählt für $n = 1, \ldots, \text{max\_n}$ die nicht-isomorphen Bäume und wie viele davon graceful sind.

**Erwartetes Ergebnis (Vermutung):** Für alle $n$ sind alle Bäume graceful.

---

#### `verify_known_classes() → Dict[str, bool]`

Verifiziert bekannte Klassen (Pfade $P_2$–$P_8$, Sterne $K_{1,1}$–$K_{1,6}$).

---

## Beispiel-Verwendung

```python
from graceful_trees import GracefulTreeConjecture
import networkx as nx

gtc = GracefulTreeConjecture()

# P_4 ist graceful:
path, labeling = gtc.graceful_path(4)
print(gtc.is_graceful(path, labeling))  # True

# Backtracking für beliebigen Baum:
tree = nx.path_graph(7)
labeling = gtc.find_graceful_labeling(tree)
print(labeling)  # {0: 0, 1: 6, 2: 1, 3: 5, 4: 2, 5: 4, 6: 3}

# Statistik:
stats = gtc.count_graceful_trees(max_n=8)
for n, info in stats.items():
    print(f"n={n}: {info['graceful']}/{info['total']} graceful")
```

---

## Literatur

- Rosa, A. (1967): On certain valuations of the vertices of a graph. *Theory of Graphs (Internat. Symposium Rome)*, 349–355.
- Gallian, J.A. (2024): A Dynamic Survey of Graph Labeling. *Electronic Journal of Combinatorics DS6*.
- Ringel, G. (1974): *Problem 25 in Theory of Graphs and Its Applications*.
- Kotzig, A. und Rosa, A. (1970): Magic valuations of finite graphs. *Canad. Math. Bull.*
