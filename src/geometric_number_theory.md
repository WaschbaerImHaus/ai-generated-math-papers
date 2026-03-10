# Modul: Geometrische Zahlentheorie (`geometric_number_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul implementiert die Kernkonzepte der Geometrischen Zahlentheorie (Geometry of
Numbers): Gitter im $\mathbb{R}^n$, LLL-Gitterbasisreduktion (Lenstra-Lenstra-Lovász),
Minkowski-Sätze über konvexe Körper, binäre quadratische Formen mit Gauß'scher Reduktion,
Gaußkomposition von quadratischen Formen sowie Grundlagen der Kugelpackungstheorie.

## Mathematischer Hintergrund

### Gitter im $\mathbb{R}^n$

Ein **Gitter** $\Lambda \subset \mathbb{R}^n$ ist die Menge aller ganzzahligen Linearkombinationen
linear unabhängiger Basisvektoren $b_1, \ldots, b_m$:

$$\Lambda = \left\{\sum_{i=1}^m a_i b_i \;\middle|\; a_i \in \mathbb{Z}\right\}$$

Die **Gitterdeterminante** (Fundamentalbereichsvolumen) für $m = n$:

$$\det(\Lambda) = |\det(B)|, \quad B = (b_1 \mid \ldots \mid b_n)$$

Für $m < n$ (nicht-quadratisch): $\det(\Lambda) = \sqrt{\det(B^\top B)}$ (Gram-Determinante).

Das **duale Gitter** $\Lambda^* = \{x \in \mathbb{R}^n \mid \langle x, \lambda \rangle \in \mathbb{Z}\ \forall\, \lambda \in \Lambda\}$
besitzt die Basis $(B^{-1})^\top$ und erfüllt $\det(\Lambda^*) = 1 / \det(\Lambda)$.

### LLL-Algorithmus (Lenstra-Lenstra-Lovász, 1982)

Eine Basis $b_1, \ldots, b_n$ ist **LLL-reduziert** (mit Parameter $\delta$), wenn:

**Größenkondition:**
$$|\mu_{ij}| \leq \tfrac{1}{2} \quad \forall\, i > j, \qquad \mu_{ij} = \frac{\langle b_i, b_j^*\rangle}{\langle b_j^*, b_j^*\rangle}$$

**Lovász-Bedingung:**
$$\delta \|b_{k-1}^*\|^2 \leq \|b_k^* + \mu_{k,k-1} b_{k-1}^*\|^2$$

Die LLL-reduzierte Basis liefert eine Näherung des kürzesten Vektors:

$$\|b_1\| \leq 2^{(n-1)/4} \cdot \lambda_1(\Lambda)$$

Anwendungen: Kryptanalyse von Gitterproblemen, NTRU, Faktorisierung von Polynomen.

### Minkowski-Sätze

**Minkowskis Fundamentalsatz (1889):** Sei $K \subset \mathbb{R}^n$ konvex, symmetrisch
($K = -K$), dann:

$$\text{Vol}(K) > 2^n \cdot \det(\Lambda) \implies K \text{ enthält einen Gitterpunkt } \lambda \neq 0$$

**Minkowski-Schranke** für algebraische Zahlkörper $K$ vom Grad $n$ mit $r_2$ komplexen Paaren:

$$M_K = \frac{n!}{n^n} \cdot \left(\frac{4}{\pi}\right)^{r_2} \cdot \sqrt{|\Delta_K|}$$

Jede Idealklasse enthält einen ganzen Ideal-Repräsentanten mit Norm $\leq M_K$.

### Binäre quadratische Formen

Eine **binäre quadratische Form** ist $f(x, y) = ax^2 + bxy + cy^2$ mit $a, b, c \in \mathbb{Z}$.

**Diskriminante:** $D = b^2 - 4ac$

- $D < 0$, $a > 0$: positiv definit (alle Werte $> 0$ außer $(0,0)$)
- $D > 0$: indefinit

**Gauss-Reduktion** (Normalform für positiv definite Formen):

$$|b| \leq a \leq c, \quad b \geq 0 \text{ falls } |b| = a \text{ oder } a = c$$

Jede Äquivalenzklasse unter $GL_2(\mathbb{Z})$ enthält genau eine reduzierte Form.

**Klassenzahl** $h(D)$: Anzahl der Äquivalenzklassen reduzierter Formen zur Diskriminante $D$.

### Kugelpackungen

**Packungsdichte** $\Delta$ und **Kissing Number** $\tau(n)$ (Kusszahl):

| Dimension $n$ | Optimales Gitter | $\Delta$ | $\tau(n)$ |
|:---:|:---:|:---:|:---:|
| 2 | $A_2$ (Dreiecksgitter) | $\pi/(2\sqrt{3})$ | 6 |
| 3 | $D_3 = $ fcc | $\pi/(3\sqrt{2})$ | 12 |
| 8 | $E_8$ | $\pi^4/384$ | 240 |
| 24 | Leech $\Lambda_{24}$ | $\pi^{12}/479001600$ | 196560 |

---

## Klassen und Methoden

### `Lattice`

**Beschreibung:** Repräsentiert ein Gitter $\Lambda \subset \mathbb{R}^n$ durch eine
Basismatrix. Berechnet Determinante, Gram-Matrix, duales Gitter und approximiert
kürzeste Vektoren via LLL.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `__init__()` | `(basis_matrix: np.ndarray)` | Initialisiert mit Basismatrix (Zeilen = Basisvektoren) |
| `determinant()` | `() -> float` | $\det(\Lambda)$: direkt oder via Gram-Determinante |
| `fundamental_domain_volume()` | `() -> float` | $|\det(B)|$ — Volumen des Fundamentalparallelotops |
| `contains_point()` | `(point: np.ndarray, tol: float) -> bool` | Prüft ob $p = B \cdot a$ mit $a \in \mathbb{Z}^n$ |
| `shortest_vector_approx()` | `() -> np.ndarray` | LLL-Näherung des kürzesten Gittervektors |
| `successive_minima_approx()` | `(k: int) -> float` | Näherung für $\lambda_k(\Lambda)$ (k-tes Sukzessivminimum) |
| `gram_matrix()` | `() -> np.ndarray` | $G = B \cdot B^\top$ mit $G_{ij} = \langle b_i, b_j \rangle$ |
| `dual_lattice_basis()` | `() -> np.ndarray` | Basis des dualen Gitters $(B^{-1})^\top$ |

---

### `LLLReduction`

**Beschreibung:** Implementiert den LLL-Gitterbasisreduktionsalgorithmus mit
einstellbarem Lovász-Parameter $\delta$. Enthält Gram-Schmidt-Orthogonalisierung
und Reduktionsprüfung.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `__init__()` | `(delta: float = 0.75)` | Lovász-Parameter $\delta \in (1/4, 1)$ |
| `gram_schmidt_orthogonalization()` | `(basis: np.ndarray) -> tuple` | Gibt $(B^*, \mu)$: orthogonale Basis und $\mu_{ij}$-Koeffizienten |
| `reduce()` | `(basis_matrix: np.ndarray) -> np.ndarray` | Führt LLL-Reduktion durch (Größenreduktion + Lovász) |
| `is_lll_reduced()` | `(basis_matrix: np.ndarray) -> bool` | Prüft Größenkondition und Lovász-Bedingung |

---

### `QuadraticForm`

**Beschreibung:** Repräsentiert eine binäre quadratische Form $f(x,y) = ax^2 + bxy + cy^2$.
Berechnet Diskriminante, prüft Definitheit, reduziert auf Gauss-Normalform, findet
äquivalente Formen via $GL_2(\mathbb{Z})$-Transformationen.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `__init__()` | `(a: int, b: int, c: int)` | Initialisiert $f = ax^2 + bxy + cy^2$ |
| `discriminant()` | `() -> int` | $D = b^2 - 4ac$ |
| `evaluate()` | `(x: int, y: int) -> int` | Berechnet $f(x, y)$ |
| `is_positive_definite()` | `() -> bool` | $a > 0$ und $D < 0$ |
| `is_reduced()` | `() -> bool` | Gauss-Normalform: $|b| \leq a \leq c$ |
| `reduce()` | `() -> QuadraticForm` | Reduktion auf Gauss'sche Normalform |
| `equivalent_forms()` | `(max_coeff: int = 10) -> list` | $GL_2(\mathbb{Z})$-äquivalente Formen |
| `represents()` | `(n: int, max_search: int = 100) -> bool` | Prüft ob $f(x,y) = n$ lösbar |

---

### `GaussComposition`

**Beschreibung:** Implementiert die Gauß'sche Komposition von Äquivalenzklassen
binärer quadratischer Formen (Klassengruppe). Berechnet Klassenzahl und
Klassengruppen-Struktur für negative Diskriminanten.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `class_number()` | `(D: int) -> int` | $h(D)$: Anzahl reduzierter Formen zur Diskriminante $D < 0$ |
| `reduced_forms()` | `(D: int) -> list[QuadraticForm]` | Alle reduzierten Formen zur Diskriminante $D$ |
| `compose()` | `(f1: QuadraticForm, f2: QuadraticForm) -> QuadraticForm` | Gauß-Komposition: $[f_1] \cdot [f_2]$ in der Klassengruppe |
| `principal_form()` | `(D: int) -> QuadraticForm` | Hauptform $f_0$ (neutrales Element der Klassengruppe) |
| `order_in_class_group()` | `(f: QuadraticForm) -> int` | Ordnung von $[f]$ in der Klassengruppe |

---

### Standalone-Funktionen (Minkowski)

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `minkowski_convex_body_theorem()` | `(volume, det_lattice, dim) -> bool` | Prüft $\text{Vol}(K) > 2^n \det(\Lambda)$ |
| `minkowski_linear_forms()` | `(lambda_coeff, bounds) -> bool` | Überprüft ob $\prod c_i \geq 1$ (Simultane Approximation) |
| `minkowski_bound()` | `(discriminant, degree, num_complex_pairs) -> float` | $M_K = \frac{n!}{n^n} \left(\frac{4}{\pi}\right)^{r_2} \sqrt{|\Delta_K|}$ |

### Standalone-Funktionen (Kugelpackung)

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `sphere_packing_density()` | `(n: int) -> dict` | Bekannte optimale Packungsdichten für Dimension $n$ |
| `kissing_number()` | `(n: int) -> Optional[int]` | Bekannte Kusszahlen $\tau(n)$ |
| `leech_lattice_properties()` | `() -> dict` | Eigenschaften des Leech-Gitters $\Lambda_{24}$ |

---

## Beispiele

```python
import numpy as np
from geometric_number_theory import Lattice, LLLReduction, QuadraticForm
from geometric_number_theory import minkowski_convex_body_theorem, minkowski_bound

# Gitter: Z² mit gedrehter Basis
basis = np.array([[2.0, 1.0], [0.0, 3.0]])
lat = Lattice(basis)
print(lat.determinant())         # 6.0
print(lat.gram_matrix())         # [[5, 3], [3, 9]]

# Kürzester Vektor via LLL
sv = lat.shortest_vector_approx()
print(sv)   # Näherung des kürzesten Gittervektors

# LLL-Reduktion direkt
lll = LLLReduction(delta=0.75)
reduced = lll.reduce(np.array([[1.0, 1.0, 1.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0]]))

# Quadratische Formen
f = QuadraticForm(2, -1, 3)   # 2x² - xy + 3y²
print(f.discriminant())         # (-1)² - 4·2·3 = -23
print(f.is_positive_definite()) # True
f_red = f.reduce()              # Gauss-Normalform

# Klassenzahl h(-23)
from geometric_number_theory import GaussComposition
gc = GaussComposition()
print(gc.class_number(-23))     # 3

# Minkowski-Satz
print(minkowski_convex_body_theorem(volume=10.0, det_lattice=1.0, dim=2))  # True (10 > 4)

# Minkowski-Schranke für Q(√-5), Grad 2, r2=1, Δ=-20
print(minkowski_bound(-20, degree=2, num_complex_pairs=1))  # ≈ 2.86
```

## Verweise

- Verwandte Module: `linear_algebra.py`, `algebra.py`, `algorithmic_number_theory.py`
- Literatur:
  - Minkowski: *Geometrie der Zahlen* (1910)
  - Lenstra, Lenstra, Lovász: *Factoring Polynomials with Rational Coefficients* (1982)
  - Cassels: *An Introduction to the Geometry of Numbers* (Springer, 1997)
  - Conway & Sloane: *Sphere Packings, Lattices and Groups* (Springer, 1999)
  - Cox: *Primes of the Form x² + ny²* (Wiley, 1989)
