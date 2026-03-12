# hadwiger_nelson.py – Dokumentation

**Modul:** `src/hadwiger_nelson.py`
**Autor:** Michael Fuhrmann
**Stand:** Build 122 (2026-03-12)

---

## Übersicht

Das Modul behandelt das **Hadwiger-Nelson-Problem** (1950): Wie viele Farben
benötigt man, um alle Punkte der euklidischen Ebene ℝ² so zu färben, dass
keine zwei Punkte im Abstand genau 1 dieselbe Farbe erhalten?

Formal: $\chi(\mathbb{R}^2) = \min\{k : \exists\, k\text{-Färbung von } G_{\mathbb{R}^2,1}\}$

**Bekannte Schranken (Stand 2025):**

$$4 \leq \chi(\mathbb{R}^2) \leq 7$$

---

## Mathematischer Hintergrund

### Einheitsdistanzgraph

Der Graph $G_{\mathbb{R}^2,1}$ hat:
- **Knotenmenge:** alle Punkte $(x,y) \in \mathbb{R}^2$
- **Kantenmenge:** $\{u,v\} \in E \iff \|u - v\|_2 = 1$

Die **chromatische Zahl** $\chi(G)$ ist die minimale Anzahl Farben für eine
gültige Knotenfärbung (benachbarte Knoten verschiedenfarbig).

### Untere Schranken

| Graph | Knoten | Kanten | $\chi$ | Beweist |
|-------|--------|--------|---------|---------|
| Moser-Spindle | 7 | 11 | 4 | $\chi(\mathbb{R}^2) \geq 4$ |
| Golomb-Graph | 10 | 12 | 4 | $\chi(\mathbb{R}^2) \geq 4$ |
| de Grey-Graph | 1581 | ~12k | 5 | $\chi(\mathbb{R}^2) \geq 5$ |

### Moser-Spindle (Moser & Moser, 1961)

Konstruktion mit zwei gleichseitigen Dreiecken $ABC$ und $AEF$ plus Verbindungsknoten:

- $A=(0,0),\; B=(1,0),\; C=(1/2, \sqrt{3}/2)$ bilden Dreieck 1
- $D=(3/2, \sqrt{3}/2)$ hat $|BD|=|CD|=1$ → erhält Farbe $\text{col}(A)$ (erzwungen)
- $E=(\cos\theta, \sin\theta)$ mit $\cos\theta = 5/6$ und $F=(\cos(\theta+60°), \sin(\theta+60°))$
- $G = $ oberer Schnittpunkt der Einheitskreise um $E$ und $F$
- Kante $D$–$G$: $|D-G|=1$, aber beide haben erzwungene Farbe $\text{col}(A)$ → Widerspruch → $\chi \geq 4$

### Obere Schranke: Hexagonales 7-Färbungsschema

Nelson (1950) teilte die Ebene in reguläre Hexagone der Seitenlänge $s < 1/\sqrt{3} \approx 0.577$ auf und vergab 7 Farben im periodischen Muster:

$$\text{Farbe}(i,j) = (i + 2j) \bmod 7$$

Da der Durchmesser jedes Hexagons $2s < 2/\sqrt{3} < 2$ ist und gleichfarbige Hexagone Mindestabstand $\sqrt{7} \cdot s\sqrt{3} > 1$ haben, ergibt dies eine gültige Färbung.

### de Grey (2018): $\chi(\mathbb{R}^2) \geq 5$

Aubrey de Grey konstruierte 2018 einen Einheitsdistanzgraph mit 1581 Knoten,
der chromatische Zahl 5 hat. Der Beweis wurde im Polymath-Projekt verifiziert
und 2019 in *Geombinatorics* publiziert.

---

## Klassen

### `UnitDistanceGraph`

Erstellt und analysiert Einheitsdistanzgraphen.

| Methode | Beschreibung |
|---------|-------------|
| `build_from_points(points, tol)` | Erstellt Einheitsdistanzgraph aus Punktliste |
| `moser_spindle()` | Statisch: Moser-Spindle (7 Knoten, 11 Kanten, $\chi=4$) |
| `golomb_graph()` | Statisch: Golomb-Graph (10 Knoten, 12 Kanten) |
| `chromatic_number_greedy(g)` | DSATUR-Greedy-Schätzung von $\chi(G)$ |
| `chromatic_number_exact(g)` | Exakte $\chi(G)$ via Backtracking (kleine Graphen) |
| `verify_coloring(g, coloring)` | Prüft Gültigkeit einer Färbung |
| `count_unit_edges(points, tol)` | Zählt Einheitsdistanz-Kanten |

### `HadwigerNelsonProblem`

Analysiert das Hadwiger-Nelson-Problem und seine Schranken.

| Methode | Beschreibung |
|---------|-------------|
| `lower_bound_moser_spindle()` | Beweist $\chi(\mathbb{R}^2) \geq 4$ |
| `lower_bound_golomb_graph()` | Alternative Schranke $\chi(\mathbb{R}^2) \geq 4$ |
| `lower_bound_de_grey_simplified()` | Skizze de Grey-Argument für $\chi(\mathbb{R}^2) \geq 5$ |
| `upper_bound_hexagonal_7_coloring(s)` | Beweist $\chi(\mathbb{R}^2) \leq 7$ |
| `verify_hexagonal_coloring_for_unit_points(pts, s)` | Verifiziert das Schema auf Punktmenge |
| `summary()` | Zusammenfassung der bekannten Schranken |
| `get_lower_bound()` | Aktuelle untere Schranke (5) |
| `get_upper_bound()` | Aktuelle obere Schranke (7) |

---

## Verwendungsbeispiel

```python
from hadwiger_nelson import HadwigerNelsonProblem, UnitDistanceGraph

# Untere Schranke via Moser-Spindle
hn = HadwigerNelsonProblem()
result = hn.lower_bound_moser_spindle()
print(result['proof'])  # "Moser-Spindle: χ = 4..."

# Obere Schranke via Hexagonalgitter
result7 = hn.upper_bound_hexagonal_7_coloring(side_length=0.48)
print(result7['proof'])  # "Mit Hexagon-Seitenlänge s = 0.480 < ..."

# Moser-Spindle direkt
g, pts = UnitDistanceGraph.moser_spindle()
print(g.number_of_nodes(), g.number_of_edges())  # 7, 11
```

---

## Bekannte Einschränkungen

- `chromatic_number_exact()` ist nur für kleine Graphen ($|V| \leq 20$) praktikabel
- `lower_bound_de_grey_simplified()` ist eine Demonstration des Prinzips, nicht der vollständige de Grey-Beweis
- Der Golomb-Graph hat in der Implementierung chromatische Zahl ≥ 3 (abhängig von der exakten Einbettung der 12 Kanten)

---

## Offene Fragen

- Gilt $\chi(\mathbb{R}^2) = 5$, $6$ oder $7$? (Status 2025: unbekannt)
- Existiert ein Einheitsdistanzgraph mit $\chi \geq 6$?

---

## Referenzen

1. Moser, W. & Moser, L. (1961). *Solution to Problem 10.* Can. Math. Bull. **4**: 187–189.
2. Nelson, E. (1950). *Problem 60.* Unpublished.
3. de Grey, A.D.N.J. (2018). *The chromatic number of the plane is at least 5.* arXiv:1804.02385.
4. Soifer, A. (2009). *The Mathematical Coloring Book.* Springer.
5. Cranston, D. & Rabern, L. (2014). *Coloring claw-free graphs with Δ−1 colors.* arXiv:1406.1140.
