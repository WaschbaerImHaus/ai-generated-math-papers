# Modul: Algebraische Geometrie (`algebraic_geometry.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Dieses Modul implementiert grundlegende Konzepte der algebraischen Geometrie: die Untersuchung
von Lösungsmengen polynomialer Gleichungssysteme (Varietäten) mit Methoden der kommutativen
Algebra. Es deckt affine und projektive Varietäten, den Hilbert'schen Nullstellensatz, elliptische
Kurven, Schnitttheorie und Morphismen algebraischer Varietäten ab.

## Mathematischer Hintergrund

### Affine Varietät

Eine affine algebraische Varietät $V(I)$ zum Ideal $I \subseteq k[x_1, \ldots, x_n]$ ist die Menge
aller gemeinsamen Nullstellen der Erzeuger:

$$V(I) = \{ \mathbf{a} \in k^n \mid f(\mathbf{a}) = 0 \text{ für alle } f \in I \}$$

### Hilbert-Basissatz

Jedes Ideal $I \subseteq k[x_1, \ldots, x_n]$ (mit $k$ einem Körper) ist endlich erzeugt:

$$I = \langle f_1, f_2, \ldots, f_m \rangle$$

Gröbner-Basen liefern ein algorithmisches Werkzeug zur Ideal-Arithmetik (Divisionsalgorithmus).

### Hilbert'scher Nullstellensatz

**Schwacher Nullstellensatz:** Für $f_1, \ldots, f_m \in \mathbb{C}[x_1, \ldots, x_n]$ gilt:

$$V(f_1, \ldots, f_m) = \emptyset \iff 1 \in \langle f_1, \ldots, f_m \rangle$$

**Starker Nullstellensatz:** Für ein Ideal $I$ und $f$ mit $f|_{V(I)} = 0$ gilt:

$$f^r \in I \quad \text{für ein } r \geq 1$$

### Elliptische Kurven

Eine elliptische Kurve in Weierstraß-Normalform ist:

$$E: y^2 = x^3 + ax + b, \quad \Delta = -16(4a^3 + 27b^2) \neq 0$$

Die $j$-Invariante klassifiziert elliptische Kurven bis auf Isomorphie:

$$j(E) = -1728 \cdot \frac{(4a)^3}{\Delta / (-16)} = 1728 \cdot \frac{4 \cdot (4a)^3}{4a^3 + 27b^2} \cdot \frac{1}{-16}$$

### Bézout'scher Schnittmultiplizitätssatz

Zwei projektive Kurven vom Grad $d_1$ bzw. $d_2$ ohne gemeinsame Komponenten schneiden sich
(über $\overline{k}$, gezählt mit Multiplizitäten) in genau

$$d_1 \cdot d_2 \text{ Punkten}$$

### Krull-Dimension

Die Krull-Dimension von $V(I)$ entspricht der Länge der längsten Primidealkette in $k[x_1,\ldots,x_n]/I$. Algorithmisch wird sie über die Gröbner-Basis und Hilbert-Polynom berechnet.

### Glattheitskriterium (Jacobi)

Ein Punkt $\mathbf{p} \in V(I)$ ist **singulär**, wenn alle partiellen Ableitungen verschwinden:

$$\operatorname{rank} J(\mathbf{p}) < \operatorname{codim} V(I), \quad J_{ij} = \frac{\partial f_i}{\partial x_j}\bigg|_{\mathbf{p}}$$

---

## Hilfsfunktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `_to_sympy_poly` | `(expr, variables) -> Poly` | Konvertiert Ausdruck in SymPy-Polynom |
| `_jacobian_matrix` | `(polys, variables) -> list` | Berechnet Jacobi-Matrix symbolisch |
| `_eval_at_point` | `(expr, variables, point) -> float` | Wertet symbolischen Ausdruck in einem Punkt aus |

---

## Klassen und Methoden

### `AffineVariety`

Repräsentiert eine affine algebraische Varietät $V(I) \subseteq k^n$, definiert durch ein Polynomideal.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(generators: list, variables: list)` | Erzeugt Varietät aus Polynomerzeuger-Liste und Variablenliste |
| `_get_groebner` | `() -> GroebnerBasis` | Berechnet intern Gröbner-Basis (gecacht) |
| `dimension` | `() -> int` | Berechnet Krull-Dimension über Hilbert-Polynom |
| `is_irreducible` | `() -> bool` | Prüft Irreduzibilität via Radikalideal-Vergleich |
| `singular_points` | `() -> list` | Findet singuläre Punkte mittels Jacobi-Kriterium |
| `intersect` | `(other: AffineVariety) -> AffineVariety` | Schnitt zweier Varietäten (Ideal-Vereinigung) |

### `ProjectiveVariety`

Projektive algebraische Varietät $V_+(I) \subseteq \mathbb{P}^n$, definiert durch homogene Polynome.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(homogeneous_poly: list, n_vars: int)` | Erzeugt projektive Varietät |
| `_poly_degree` | `(f) -> int` | Berechnet Grad eines homogenen Polynoms |
| `degree` | `() -> int` | Gibt den Grad der Varietät zurück |
| `is_projective_smooth` | `() -> bool` | Prüft Glattheit via projektives Jacobi-Kriterium |

### `HilbertBasisTheorem`

Demonstration und algorithmische Umsetzung des Hilbert-Basissatzes via Gröbner-Basen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `groebner_basis` | `(generators: list, variables: list) -> list` | Berechnet Gröbner-Basis (lex-Ordnung) |
| `ideal_membership` | `(poly, generators: list, variables: list) -> bool` | Prüft Ideal-Mitgliedschaft via Normalform |

### `NullstellensatzDemo`

Demonstriert und prüft den Hilbert'schen Nullstellensatz.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `weak_nullstellensatz` | `(generators: list, variables: list) -> dict` | Schwacher Nullstellensatz: prüft ob $V = \emptyset$ und ob $1 \in I$ |
| `strong_nullstellensatz` | `(f, generators: list, variables: list) -> dict` | Starker Nullstellensatz: prüft Radikalmitgliedschaft |

### `EllipticCurveGeometry`

Elliptische Kurve $y^2 = x^3 + ax + b$ mit algebraisch-geometrischen Methoden.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(a, b)` | Erzeugt elliptische Kurve mit Koeffizienten $a, b$ |
| `discriminant` | `() -> Expr` | Berechnet Diskriminante $\Delta = -16(4a^3 + 27b^2)$ |
| `j_invariant` | `() -> Expr` | Berechnet $j$-Invariante der Kurve |
| `is_smooth` | `() -> bool` | Prüft Glattheit: $\Delta \neq 0$ |
| `rational_points_mod_p` | `(p: int) -> list` | Berechnet alle rationalen Punkte $\pmod{p}$ |

### `IntersectionTheory`

Schnitttheorie nach Bézout und lokale Schnittmultiplizitäten.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `bezout_theorem` | `(d1: int, d2: int) -> int` | Bézout-Schnittzahl: $d_1 \cdot d_2$ |
| `intersection_multiplicity` | `(f, g, point: dict, variables: list) -> int` | Lokale Schnittmultiplizität in einem Punkt |

### `MorphismOfVarieties`

Morphismus $\varphi: V \to W$ zwischen algebraischen Varietäten.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(source: AffineVariety, target: AffineVariety, coord_map: list)` | Definiert Morphismus durch Koordinatenabbildung |
| `is_dominant` | `() -> bool` | Prüft Dominanz (dichtes Bild) |
| `pullback` | `(poly) -> object` | Berechnet Pullback $\varphi^*(f) = f \circ \varphi$ |

---

## Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `compute_groebner_basis` | `(gens: list, vars: list) -> list` | Berechnet Gröbner-Basis eines Ideals |
| `ideal_radical` | `(gens: list, vars: list) -> list` | Berechnet das Radikalideal $\sqrt{I}$ |
| `zariski_closure_demo` | `(points: list, vars: list) -> dict` | Demonstriert Zariski-Abschluss einer Punktmenge |
| `hilbert_polynomial` | `(gens: list, vars: list) -> Expr` | Berechnet Hilbert-Polynom des Quotientenrings |

---

## Beispiele

```python
from sympy import symbols
from algebraic_geometry import AffineVariety, EllipticCurveGeometry, IntersectionTheory

x, y = symbols('x y')

# Affine Varietät: Kreis x² + y² - 1 = 0
V = AffineVariety([x**2 + y**2 - 1], [x, y])
print("Dimension:", V.dimension())          # 1 (Kurve)
print("Irreduzibel:", V.is_irreducible())   # True

# Elliptische Kurve y² = x³ - x
E = EllipticCurveGeometry(-1, 0)
print("Diskriminante:", E.discriminant())   # -16·(4·(-1)³ + 27·0²) = 64
print("Glatt:", E.is_smooth())              # True
print("Punkte mod 5:", E.rational_points_mod_p(5))

# Bézout-Schnittzahl
it = IntersectionTheory()
print("Schnittpunkte Grad 3 × Grad 4:", it.bezout_theorem(3, 4))  # 12

# Hilbert-Basissatz
from algebraic_geometry import HilbertBasisTheorem
gb = HilbertBasisTheorem.groebner_basis([x**2 - y, x**3 - x], [x, y])
print("Gröbner-Basis:", gb)
```

---

## Technische Hinweise

- **Gröbner-Basen** werden intern mit SymPy (`sympy.polys.groebnertools`) berechnet.
- Die **Krull-Dimension** wird über das Hilbert-Polynom des assoziierten graduierten Rings approximiert.
- **Singuläre Punkte** werden numerisch mit dem Jacobi-Kriterium gefunden.
- Das **Jacobi-Kriterium** für projektive Glattheit prüft alle affinen Karten.
