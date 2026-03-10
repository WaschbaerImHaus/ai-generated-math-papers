# topology.py – Topologie und Geometrie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Modul:** `src/topology.py`

---

## Übersicht

Das Modul `topology.py` implementiert grundlegende Konzepte der Topologie und Geometrie:

- **Metrische Räume** mit axiomatischer Überprüfung
- **Standardmetriken** (euklidisch, Manhattan, Chebyshev, Lp)
- **Topologische Eigenschaften** (Zusammenhang, Hausdorff-Abstand, Kompaktheit)
- **Parametrische Kurven** (Bogenlänge, Krümmung, Umlaufzahl)
- **Simpliziale Homologie** (Euler-Charakteristik, Betti-Zahlen)
- **Fraktale Dimensionen** (Box-Counting, Cantor-Menge, Sierpinski-Dreieck)

---

## Mathematischer Hintergrund

### Metrische Räume

Ein metrischer Raum $(X, d)$ ist eine Menge $X$ zusammen mit einer Abstandsfunktion $d: X \times X \to \mathbb{R}_{\geq 0}$, die folgende **vier Axiome** erfüllt:

| Axiom | Name | Formel |
|-------|------|--------|
| M1 | Nichtnegativität | $d(x,y) \geq 0$ |
| M2 | Definitheit | $d(x,y) = 0 \iff x = y$ |
| M3 | Symmetrie | $d(x,y) = d(y,x)$ |
| M4 | Dreiecksungleichung | $d(x,z) \leq d(x,y) + d(y,z)$ |

### Standardmetriken im $\mathbb{R}^n$

Die wichtigsten Metriken sind Spezialfälle der $L^p$-Norm:

$$d_p(x, y) = \left( \sum_{i=1}^{n} |x_i - y_i|^p \right)^{1/p}$$

| Metrik | $p$ | Formel | Beschreibung |
|--------|-----|--------|-------------|
| Manhattan | 1 | $\sum |x_i - y_i|$ | Taxicab-Abstand |
| Euklidisch | 2 | $\sqrt{\sum (x_i-y_i)^2}$ | Luftlinie |
| Chebyshev | $\infty$ | $\max |x_i - y_i|$ | Schachbrett |
| Diskret | – | $0$ oder $1$ | Identitätsprüfung |

**Norm-Äquivalenz in $\mathbb{R}^n$:**
$$d_\infty(x,y) \leq d_2(x,y) \leq d_1(x,y) \leq \sqrt{n} \cdot d_2(x,y)$$

---

## Klassen und Funktionen

### `MetricSpace(distance_func)`

Kapselt einen metrischen Raum mit einer Metrikfunktion.

```python
space = MetricSpace(euclidean_metric)
result = space.verify_metric_axioms([[0,0], [1,0], [0,1]])
# → {'non_negativity': True, 'identity': True, 'symmetry': True, 'triangle_inequality': True}
```

**Methoden:**

| Methode | Beschreibung |
|---------|-------------|
| `verify_metric_axioms(points)` | Prüft alle 4 Axiome numerisch |
| `open_ball(center, radius, points)` | $B(x,r) = \{y : d(x,y) < r\}$ |
| `is_open_set(subset, all_points)` | Prüft Offenheit in diskreter Menge |
| `is_cauchy_sequence(sequence)` | Prüft Cauchy-Bedingung |

### Topologische Funktionen

| Funktion | Beschreibung |
|----------|-------------|
| `is_connected(points, metric, radius)` | $\varepsilon$-Zusammenhang via Union-Find |
| `hausdorff_distance(A, B, metric)` | $H(A,B) = \max(\vec{d}(A,B), \vec{d}(B,A))$ |
| `compute_diameter(points, metric)` | $\text{diam}(A) = \sup_{x,y \in A} d(x,y)$ |
| `is_compact_discrete(points, metric, eps)` | Epsilon-Netz für endliche Mengen |

---

## Parametrische Kurven

### Theorie

Eine parametrische Kurve ist eine stetige Abbildung:
$$\gamma: [a,b] \to \mathbb{R}^n, \quad t \mapsto (\gamma_1(t), \ldots, \gamma_n(t))$$

#### Bogenlänge

$$L(\gamma) = \int_a^b |\gamma'(t)| \, dt = \int_a^b \sqrt{\sum_{i=1}^n \dot{\gamma}_i(t)^2} \, dt$$

Numerisch via **Trapezregel** mit $n$ Teilintervallen.

#### Krümmung

**2D:**
$$\kappa(t) = \frac{|x'y'' - y'x''|}{(x'^2 + y'^2)^{3/2}}$$

**3D:**
$$\kappa(t) = \frac{|\gamma' \times \gamma''|}{|\gamma'|^3}$$

#### Umlaufzahl

Für eine geschlossene Kurve $\gamma$ und Punkt $p \notin \gamma$:
$$n(\gamma, p) = \frac{1}{2\pi} \oint_\gamma d\theta = \frac{1}{2\pi i} \oint_\gamma \frac{dz}{z - p}$$

### Vorgefertigte Kurven

| Funktion | Parametrisierung | Beschreibung |
|----------|-----------------|-------------|
| `circle_curve(r, center)` | $(r\cos t, r\sin t)$ | Kreis |
| `lissajous_curve(a, b, δ)` | $(\sin(at+\delta), \sin(bt))$ | Lissajous-Figur |
| `helix_curve(r, pitch)` | $(r\cos t, r\sin t, \text{pitch}\cdot t/2\pi)$ | Schraubenlinie |

---

## Simpliziale Homologie

### Euler-Charakteristik

Die **Euler-Charakteristik** ist eine fundamentale topologische Invariante:

$$\chi = V - E + F$$

Bekannte Werte:

| Fläche | $\chi$ | Geschlecht $g$ |
|--------|--------|---------------|
| Sphäre $S^2$ | 2 | 0 |
| Torus $T^2$ | 0 | 1 |
| Doppeltorus | -2 | 2 |
| Kleinsche Flasche | 0 | – |
| Projektive Ebene $\mathbb{R}P^2$ | 1 | – |

**Euler-Polyeder-Formel** (Euler 1758): Für jeden konvexen Polyeder gilt $\chi = 2$.

### Geschlecht

$$\text{Orientierbar:} \quad \chi = 2 - 2g \quad \Rightarrow \quad g = \frac{2-\chi}{2}$$
$$\text{Nicht-orientierbar:} \quad \chi = 2 - g \quad \Rightarrow \quad g = 2 - \chi$$

### Betti-Zahlen für Graphen

$$\beta_0 = \text{Anzahl Zusammenhangskomponenten}$$
$$\beta_1 = |E| - |V| + \beta_0 \quad \text{(unabhängige Zyklen)}$$

---

## Fraktale Dimensionen

### Box-Counting-Dimension (Minkowski-Bouligand)

$$d_{\text{box}} = -\lim_{\varepsilon \to 0} \frac{\log N(\varepsilon)}{\log \varepsilon}$$

$N(\varepsilon)$: Anzahl der Boxen der Größe $\varepsilon$, die mindestens einen Punkt enthalten.

Numerisch via **lineare Regression** von $\log N$ gegen $\log(1/\varepsilon)$.

### Selbstähnlichkeitsdimension

Für selbstähnliche Fraktale mit $N$ Kopien der Skalierung $r$:

$$d = \frac{\log N}{\log(1/r)}$$

| Fraktal | $N$ | $r$ | Dimension |
|---------|-----|-----|-----------|
| Cantor-Menge | 2 | 1/3 | $\log 2 / \log 3 \approx 0.6309$ |
| Sierpinski-Dreieck | 3 | 1/2 | $\log 3 / \log 2 \approx 1.585$ |
| Koch-Kurve | 4 | 1/3 | $\log 4 / \log 3 \approx 1.262$ |

---

## Verwendungsbeispiele

```python
from src.topology import (
    euclidean_metric, MetricSpace, circle_curve,
    euler_characteristic_polygon, hausdorff_dimension_cantor
)

# Metrischer Raum
space = MetricSpace(euclidean_metric)
axioms = space.verify_metric_axioms([[0,0], [1,0], [0,1], [1,1]])
print(axioms)  # Alle True

# Kreiskurve
kreis = circle_curve(radius=2.0)
print(f"Umfang: {kreis.arc_length():.4f}")      # ≈ 12.5664 = 4π
print(f"Geschlossen: {kreis.is_closed()}")        # True
print(f"Krümmung: {kreis.curvature(0.0):.4f}")   # ≈ 0.5 = 1/r
print(f"Umlaufzahl: {kreis.winding_number([0,0])}")  # 1

# Würfel
chi = euler_characteristic_polygon(8, 12, 6)
print(f"Euler-Charakteristik des Würfels: {chi}")  # 2

# Cantor-Dimension
d = hausdorff_dimension_cantor()
print(f"Hausdorff-Dimension: {d:.4f}")  # 0.6309
```

---

## Tests

Die Tests liegen in `/tests/test_topology.py` und decken folgende Bereiche ab:

- Alle 5 Metriktypen mit Randwerten
- MetricSpace: Axiomprüfung, offene Kugeln, Cauchy-Folgen
- Topologische Eigenschaften: Zusammenhang, Hausdorff-Abstand, Durchmesser
- ParametricCurve: Auswertung, Bogenlänge, Krümmung, Abgeschlossenheit, Umlaufzahl
- Lissajous, Helix
- Simpliziale Homologie: Euler-Formel, Geschlecht, Betti-Zahlen
- Fraktale Dimensionen: Cantor, Sierpinski, Box-Counting

---

## Referenzen

- Munkres, J. R. (2000). *Topology* (2nd ed.). Prentice Hall.
- Falconer, K. (2003). *Fractal Geometry: Mathematical Foundations and Applications*. Wiley.
- Hatcher, A. (2002). *Algebraic Topology*. Cambridge University Press.
- do Carmo, M. P. (1976). *Differential Geometry of Curves and Surfaces*. Prentice Hall.
