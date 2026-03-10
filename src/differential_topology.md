# Modul: Differentialtopologie (`differential_topology.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Dieses Modul implementiert grundlegende Konzepte der Differentialtopologie: glatte Mannigfaltigkeiten,
Differentialformen mit äußerer Ableitung, de Rham-Kohomologie, Morse-Theorie, Transversalitätsbedingungen,
Vektorbündel, charakteristische Klassen und die Topologie der Sphären $S^n$.

## Mathematischer Hintergrund

### Glatte Mannigfaltigkeiten

Eine glatte $n$-dimensionale **Mannigfaltigkeit** $M$ ist ein topologischer Raum, der lokal
homöomorph zu $\mathbb{R}^n$ ist, mit glatten Kartenwechseln. Tangential- und Kotangentialraum
haben jeweils Dimension $n$:

$$T_p M \cong \mathbb{R}^n, \quad T^*_p M \cong (\mathbb{R}^n)^*$$

### Differentialformen und äußere Ableitung

Eine **$k$-Differentialform** $\omega$ auf $U \subseteq \mathbb{R}^n$ ist ein glatter Schnitt im
$k$-ten äußeren Potenzprodukt $\bigwedge^k T^* U$. Das **äußere Produkt** (Wedge-Produkt):

$$(\alpha \wedge \beta)(v_1, \ldots, v_{p+q}) = \frac{1}{p!\, q!} \sum_{\sigma \in S_{p+q}} \operatorname{sgn}(\sigma)\, \alpha(v_{\sigma(1)}, \ldots) \beta(\ldots, v_{\sigma(p+q)})$$

Die **äußere Ableitung** $d: \Omega^k \to \Omega^{k+1}$ erfüllt $d^2 = 0$:

$$d\omega = \sum_{I} \sum_{j} \frac{\partial f_I}{\partial x_j} dx_j \wedge dx_{i_1} \wedge \cdots \wedge dx_{i_k}$$

### de Rham-Kohomologie

Die **de Rham-Kohomologiegruppen** messen Hindernisse für exakte Formen:

$$H^k_{\mathrm{dR}}(M) = \frac{\ker(d: \Omega^k \to \Omega^{k+1})}{\operatorname{im}(d: \Omega^{k-1} \to \Omega^k)}$$

Eine Form heißt **geschlossen**, falls $d\omega = 0$; **exakt**, falls $\omega = d\eta$.
Nach dem de Rham-Satz gilt $H^k_{\mathrm{dR}}(M) \cong H^k(M; \mathbb{R})$.

**Betti-Zahlen** für die $n$-Sphäre $S^n$ ($n \geq 1$):

$$b_k(S^n) = \begin{cases} 1 & k = 0 \text{ oder } k = n \\ 0 & \text{sonst} \end{cases}$$

**Euler-Charakteristik:**

$$\chi(M) = \sum_{k=0}^{n} (-1)^k b_k$$

### Morse-Theorie

Eine glatte Funktion $f: M \to \mathbb{R}$ ist eine **Morse-Funktion**, wenn alle kritischen Punkte
nicht-degeneriert sind ($\det H_f(p) \neq 0$). Der **Morse-Index** eines kritischen Punkts
ist die Anzahl negativer Eigenwerte der Hesse-Matrix:

$$H_f(p) = \left(\frac{\partial^2 f}{\partial x_i \partial x_j}\right)_{ij}$$

**Morse-Ungleichungen:** Mit $c_k$ = Anzahl der kritischen Punkte vom Index $k$:

$$c_k \geq b_k, \quad \sum_{k} (-1)^k c_k = \chi(M)$$

**Henkelzerlegung:** Jede Mannigfaltigkeit lässt sich aus Henkeln (Handlebodies) der Dimensionen
$0, 1, \ldots, n$ aufbauen; ein kritischer Punkt vom Index $\lambda$ entspricht einem $\lambda$-Henkel.

### Transversalität und Sard'sches Theorem

Zwei Untermannigfaltigkeiten $A, B \subseteq M$ sind **transversal** ($A \pitchfork B$), wenn:

$$T_p A + T_p B = T_p M \quad \text{für alle } p \in A \cap B$$

**Sard'sches Theorem:** Die Menge der kritischen Werte einer glatten Abbildung hat
Lebesgue-Maß Null; reguläre Werte liegen dicht.

**Schittdimension:** Wenn $A \pitchfork B$, dann $\dim(A \cap B) = \dim A + \dim B - \dim M$.

### Charakteristische Klassen

- **Chern-Klassen** $c_k(E) \in H^{2k}(M; \mathbb{Z})$ für komplexe Vektorbündel $E$
- **Pontryagin-Klassen** $p_k(E) \in H^{4k}(M; \mathbb{Z})$ für reelle Bündel
- **Euler-Klasse** $e(E) \in H^n(M; \mathbb{Z})$ für orientierte reelle Bündel vom Rang $n$
- **Stiefel-Whitney-Klassen** $w_k(E) \in H^k(M; \mathbb{Z}/2)$

### Satz von Poincaré-Hopf

Für ein glattes Vektorfeld $X$ auf einer kompakten Mannigfaltigkeit $M$:

$$\sum_p \operatorname{ind}(X, p) = \chi(M)$$

**Hairy-Ball-Theorem:** Jedes tangentiale Vektorfeld auf $S^{2k}$ hat mindestens eine Nullstelle.

### Whitney-Einbettungssatz

Jede glatte $n$-Mannigfaltigkeit lässt sich glatt in $\mathbb{R}^{2n}$ einbetten:

$$\dim_{\text{Einbettung}} = 2n$$

---

## Klassen und Methoden

### `SmoothManifold`

Repräsentiert eine abstrakte $n$-dimensionale glatte Mannigfaltigkeit.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(dimension: int, name: str = 'M')` | Erzeugt $n$-dimensionale Mannigfaltigkeit |
| `tangent_space_dim` | `() -> int` | Dimension des Tangentialraums: $n$ |
| `cotangent_space_dim` | `() -> int` | Dimension des Kotangentialraums: $n$ |
| `__repr__` | `() -> str` | Textuelle Darstellung |

### `DifferentialForm`

$k$-Differentialform auf einer offenen Menge des $\mathbb{R}^n$ mit SymPy-Ausdrücken.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(degree: int, coefficients: list, variables: list)` | Erzeugt $k$-Form mit Koeffizienten und Variablen |
| `wedge` | `(other: DifferentialForm) -> DifferentialForm` | Äußeres Produkt $\omega \wedge \eta$ |
| `exterior_derivative` | `() -> DifferentialForm` | Äußere Ableitung $d\omega$ |
| `is_closed` | `() -> bool` | Prüft ob $d\omega = 0$ (geschlossen) |
| `is_exact` | `() -> bool` | Prüft ob $\omega = d\eta$ für eine Form $\eta$ (exakt) |
| `__repr__` | `() -> str` | Textuelle Darstellung |

### `DeRhamCohomology`

Berechnet de Rham-Kohomologie und Betti-Zahlen für klassische Räume.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `closed_forms` | `(degree: int, forms: List[DifferentialForm]) -> list` | Filtert geschlossene $k$-Formen |
| `exact_forms` | `(degree: int, forms: List[DifferentialForm]) -> list` | Filtert exakte $k$-Formen |
| `betti_numbers_sphere` | `(n: int) -> List[int]` | Betti-Zahlen von $S^n$ |
| `betti_numbers_torus` | `(n: int) -> List[int]` | Betti-Zahlen des $n$-Torus $T^n$ |
| `euler_characteristic` | `(betti_numbers: List[int]) -> int` | Euler-Charakteristik aus Betti-Zahlen |

### `MorseTheory`

Morse-Theorie: kritische Punkte, Morse-Index und Henkelzerlegung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `morse_function_critical_points` | `(f_expr, variables: List) -> list` | Kritische Punkte ($\nabla f = 0$) via SymPy |
| `morse_index` | `(f_expr, variables: List, point: dict) -> int` | Morse-Index (neg. Eigenwerte der Hesse-Matrix) |
| `is_morse_function` | `(f_expr, variables: List) -> bool` | Prüft ob alle kritischen Punkte nicht-degeneriert sind |
| `morse_inequalities_check` | `(critical_points: list, betti_numbers: list) -> dict` | Überprüft Morse-Ungleichungen |
| `handle_decomposition_demo` | `(f_expr, variables: List) -> Dict` | Henkelzerlegung aus kritischen Punkten |

### `TransversalityTheory`

Transversalitätsbedingungen und Sard'sches Theorem.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `are_transverse` | `(dim_A, dim_B, dim_M, ...) -> bool` | Prüft Transversalitätsbedingung |
| `transverse_intersection_dimension` | `(dim_A, dim_B, dim_M) -> int` | Dimension $\dim(A \cap B) = \dim A + \dim B - \dim M$ |
| `sard_theorem_demo` | `(...) -> str` | Demonstriert Sard'sches Theorem |

### `VectorBundle`

Vektorbündel $E \to B$ über einer Basis $B$ mit Faser $\mathbb{R}^r$ (oder $\mathbb{C}^r$).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(base_dim, fiber_dim, total_space_name)` | Erzeugt Vektorbündel |
| `total_space_dim` | `() -> int` | $\dim E = \dim B + \text{Rang}$ |
| `rank` | `() -> int` | Rang des Bündels |
| `is_trivial` | `() -> bool` | Prüft Trivialität |
| `__repr__` | `() -> str` | Textuelle Darstellung |

### `CharacteristicClasses`

Charakteristische Klassen: Chern, Pontryagin, Euler, Stiefel-Whitney.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `chern_class_demo` | `(n: int) -> Dict` | Chern-Klassen $c_k$ für ein Vektorbündel über $\mathbb{CP}^n$ |
| `pontryagin_class_demo` | `(n: int) -> Dict` | Pontryagin-Klassen $p_k$ für reelle Bündel |
| `euler_class_demo` | `(dimension: int) -> Dict` | Euler-Klasse $e(E)$ |
| `stiefel_whitney_demo` | `(n: int) -> Dict` | Stiefel-Whitney-Klassen $w_k$ |

### `SphereTopology`

Topologie der $n$-Sphäre $S^n$: Homotopie-, Homologiegruppen, Orientierbarkeit.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(n: int)` | Erzeugt $S^n$ |
| `homotopy_groups_low` | `(n: Optional[int] = None) -> Dict` | Homotopiegruppen $\pi_k(S^n)$ für $k \leq n+3$ |
| `is_orientable` | `() -> bool` | Orientierbarkeit (immer `True` für $S^n$) |
| `homology_groups` | `(n: Optional[int] = None) -> Dict` | Homologiegruppen $H_k(S^n; \mathbb{Z})$ |
| `__repr__` | `() -> str` | Textuelle Darstellung |

---

## Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `degree_of_map` | `(f_values, g_values) -> int` | Abbildungsgrad einer glatten Abbildung $S^n \to S^n$ |
| `hairy_ball_theorem_demo` | `() -> Dict` | Demonstriert den Hairy-Ball-Satz für gerade $n$ |
| `poincare_hopf_theorem_demo` | `(euler_char, vector_field) -> Dict` | Poincaré-Hopf: $\sum \operatorname{ind} = \chi$ |
| `whitney_embedding_dimension` | `(manifold_dim: int) -> int` | Whitney-Einbettungsdimension $= 2n$ |
| `compute_jacobian` | `(f_exprs, variables) -> Matrix` | Jacobi-Matrix via SymPy |
| `implicit_function_theorem_check` | `(F, vars, point) -> dict` | Prüft Voraussetzungen des impliziten Funktionensatzes |

---

## Beispiele

```python
from sympy import symbols, sin, cos
from differential_topology import (
    SmoothManifold, DifferentialForm, DeRhamCohomology,
    MorseTheory, SphereTopology, whitney_embedding_dimension
)

# Glatte Mannigfaltigkeit
M = SmoothManifold(3, 'M')
print(M)  # Smooth 3-dimensional manifold M

# Differentialform: ω = x dx + y dy (1-Form in R²)
x, y = symbols('x y')
omega = DifferentialForm(1, [x, y], [x, y])
print("dω:", omega.exterior_derivative())     # Prüft auf 2-Form
print("Geschlossen:", omega.is_closed())      # Prüft ob dω=0

# de Rham-Kohomologie von S³
dr = DeRhamCohomology()
b = dr.betti_numbers_sphere(3)
print("Betti-Zahlen S³:", b)        # [1, 0, 0, 1]
print("χ(S³):", dr.euler_characteristic(b))  # 0

# Morse-Theorie: Höhenfunktion auf S²
x, y, z = symbols('x y z')
f = z  # Morse-Funktion (Nord- und Südpol sind kritische Punkte)
mt = MorseTheory()
pts = mt.morse_function_critical_points(
    x**2 + y**2 - 1,  # Randbedingung
    [x, y]
)

# Sphärentopologie
S4 = SphereTopology(4)
print("π_k(S⁴):", S4.homotopy_groups_low())
print("H_k(S⁴):", S4.homology_groups())
print("Orientierbar:", S4.is_orientable())  # True

# Whitney-Einbettung
print("Einbettungsdimension für n=5:", whitney_embedding_dimension(5))  # 10
```

---

## Technische Hinweise

- **Äußere Ableitung** wird symbolisch mit SymPy berechnet (partielle Ableitungen und Permutationsvorzeichen).
- **Kritische Punkte** werden über `sympy.solve` des Gradienten berechnet; bei nicht-elementaren Systemen kann dies langsam sein.
- **Morse-Index** wird numerisch über Eigenwerte der numerisch ausgewerteten Hesse-Matrix bestimmt.
- **Homotopiegruppen** $\pi_k(S^n)$ sind für $k \leq n+3$ tabellenbasiert; darüber hinaus sind sie häufig unbekannt oder sehr schwer berechenbar.
- **Charakteristische Klassen** werden algebraisch/demonstrativ zurückgegeben (keine vollständige Kohomologieberechnung).
