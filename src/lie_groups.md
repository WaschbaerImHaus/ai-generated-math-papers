# Modul: Lie-Gruppen und Lie-Algebren (`lie_groups.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Dieses Modul implementiert die Theorie der Lie-Gruppen und Lie-Algebren: Matrixgruppen wie
$SO(n)$, $SU(n)$ und $GL(n)$, die zugehörigen Lie-Algebren, die Exponentialabbildung,
die Baker-Campbell-Hausdorff-Formel sowie Darstellungstheorie und Cartan-Klassifikation
einfacher Lie-Algebren — ohne externe Lie-Bibliotheken.

## Mathematischer Hintergrund

### Lie-Gruppe

Eine **Lie-Gruppe** $G$ ist gleichzeitig eine differenzierbare Mannigfaltigkeit und eine Gruppe,
wobei Multiplikation und Inversion glatte Abbildungen sind. Wichtige Matrixgruppen:

$$GL(n, \mathbb{R}) = \{ A \in M_n(\mathbb{R}) \mid \det A \neq 0 \}$$

$$SO(n) = \{ A \in GL(n, \mathbb{R}) \mid A^\top A = I,\, \det A = 1 \}$$

$$SU(n) = \{ A \in GL(n, \mathbb{C}) \mid A^\dagger A = I,\, \det A = 1 \}$$

Dimensionen: $\dim SO(n) = \tfrac{n(n-1)}{2}$, $\dim SU(n) = n^2 - 1$.

### Lie-Algebra

Die **Lie-Algebra** $\mathfrak{g} = T_e G$ ist der Tangentialraum am Einselement $e$.
Sie wird durch die **Lie-Klammer** (Kommutator) zu einer Algebra:

$$[X, Y] = XY - YX$$

Die Lie-Algebren der klassischen Gruppen:
$$\mathfrak{so}(n) = \{ X \in M_n(\mathbb{R}) \mid X + X^\top = 0 \}, \quad \mathfrak{su}(n) = \{ X \in M_n(\mathbb{C}) \mid X + X^\dagger = 0,\, \operatorname{tr} X = 0 \}$$

### Exponentialabbildung

Die **Exponentialabbildung** verbindet Lie-Algebra und Lie-Gruppe:

$$\exp: \mathfrak{g} \to G, \quad X \mapsto e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}$$

Für Matrixgruppen: $e^X \in G$ wenn $X \in \mathfrak{g}$.

### Baker-Campbell-Hausdorff-Formel

Das Produkt zweier Gruppenelemente in Lie-Algebra-Koordinaten:

$$\log(e^X e^Y) = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}\bigl([X,[X,Y]] + [Y,[Y,X]]\bigr) + \cdots$$

### Killing-Form und Halbeinfachheit

Die **Killing-Form** auf $\mathfrak{g}$ ist die symmetrische Bilinearform:

$$\kappa(X, Y) = \operatorname{tr}(\operatorname{ad}_X \circ \operatorname{ad}_Y), \quad (\operatorname{ad}_X)(Z) = [X, Z]$$

Nach dem Satz von Cartan: $\mathfrak{g}$ ist **halbeinfach** $\iff$ $\kappa$ ist nicht-degeneriert.

### Cartan-Klassifikation

Einfache komplexe Lie-Algebren werden durch Dynkin-Diagramme klassifiziert:

| Typ | Algebra | Gruppe | Dimension |
|-----|---------|--------|-----------|
| $A_n$ | $\mathfrak{sl}(n+1)$ | $SU(n+1)$ | $n(n+2)$ |
| $B_n$ | $\mathfrak{so}(2n+1)$ | $SO(2n+1)$ | $n(2n+1)$ |
| $C_n$ | $\mathfrak{sp}(2n)$ | $Sp(2n)$ | $n(2n+1)$ |
| $D_n$ | $\mathfrak{so}(2n)$ | $SO(2n)$ | $n(2n-1)$ |

Ausnahmen: $G_2, F_4, E_6, E_7, E_8$.

### Casimir-Operator

Für eine Darstellung $\rho$ von $\mathfrak{g}$ mit Basis $\{T_a\}$:

$$C_2 = \sum_{a,b} g^{ab} \rho(T_a) \rho(T_b), \quad g_{ab} = \kappa(T_a, T_b)$$

Der Casimir-Operator vertauscht mit allen Darstellungsmatrizen (Schur-Lemma: $C_2 = \lambda \cdot I$).

---

## Klassen und Methoden

### `LieGroup` (abstrakte Basisklasse)

Abstrakte Basisklasse für alle Lie-Gruppen mit Pflichtmethoden.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `dimension` | `() -> int` | Dimension der Lie-Gruppe als Mannigfaltigkeit |
| `is_compact` | `() -> bool` | Gibt an, ob die Gruppe kompakt ist |
| `is_connected` | `() -> bool` | Gibt an, ob die Gruppe zusammenhängend ist |

### `MatrixLieGroup`

Allgemeine Matrixgruppe (Untergruppe von $GL(n)$).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(n: int, field: str = 'real')` | Erzeugt $n \times n$ Matrixgruppe über `'real'` oder `'complex'` |
| `dimension` | `() -> int` | Dimension der Gruppe |
| `is_compact` | `() -> bool` | Kompaktheit |
| `is_connected` | `() -> bool` | Zusammenhang |
| `identity` | `() -> np.ndarray` | Einheitsmatrix $I_n$ |
| `compose` | `(A, B) -> np.ndarray` | Gruppenoperation: Matrixprodukt $AB$ |
| `inverse` | `(A) -> np.ndarray` | Inverses $A^{-1}$ |
| `is_element` | `(A, tol=1e-10) -> bool` | Prüft Gruppenzugehörigkeit |

### `SO`

Spezielle orthogonale Gruppe $SO(n)$: Drehungen im $n$-dimensionalen reellen Raum.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(n: int)` | Erzeugt $SO(n)$ |
| `dimension` | `() -> int` | $\dim SO(n) = n(n-1)/2$ |
| `is_compact` | `() -> bool` | Immer `True` |
| `is_connected` | `() -> bool` | `True` für $n \geq 2$ |
| `is_element` | `(A, tol=1e-10) -> bool` | Prüft $A^\top A = I$ und $\det A = 1$ |
| `rotation_2d` | `(theta: float) -> np.ndarray` | 2D-Rotationsmatrix $R(\theta)$ |
| `rotation_3d_x` | `(theta: float) -> np.ndarray` | 3D-Drehung um $x$-Achse |
| `rotation_3d_y` | `(theta: float) -> np.ndarray` | 3D-Drehung um $y$-Achse |
| `rotation_3d_z` | `(theta: float) -> np.ndarray` | 3D-Drehung um $z$-Achse |

### `SU`

Spezielle unitäre Gruppe $SU(n)$: unitäre Transformationen mit Determinante 1.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(n: int)` | Erzeugt $SU(n)$ |
| `dimension` | `() -> int` | $\dim SU(n) = n^2 - 1$ |
| `is_compact` | `() -> bool` | Immer `True` |
| `is_connected` | `() -> bool` | Immer `True` |
| `is_element` | `(A, tol=1e-10) -> bool` | Prüft $A^\dagger A = I$ und $\det A = 1$ |
| `su2_from_angles` | `(theta, phi, psi) -> np.ndarray` | $SU(2)$-Element aus Euler-Winkeln |

### `GL`

Allgemeine lineare Gruppe $GL(n, \mathbb{R/C})$.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(n: int, field: str = 'real')` | Erzeugt $GL(n)$ |
| `is_element` | `(A, tol=1e-10) -> bool` | Prüft $\det A \neq 0$ |

### `LieAlgebra`

Lie-Algebra $\mathfrak{g}$, gegeben durch eine Basis von Matrizen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(basis_matrices: list)` | Erzeugt Lie-Algebra mit gegebener Matrixbasis |
| `bracket` | `(X, Y) -> np.ndarray` | Lie-Klammer $[X, Y] = XY - YX$ |
| `adjoint_representation` | `(X) -> np.ndarray` | Adjungierte Darstellung $\operatorname{ad}_X: Y \mapsto [X,Y]$ |
| `killing_form` | `(X, Y) -> complex` | Killing-Form $\kappa(X,Y) = \operatorname{tr}(\operatorname{ad}_X \operatorname{ad}_Y)$ |
| `is_semisimple` | `(tol=1e-10) -> bool` | Halbeinfachheit: Killing-Form nicht-degeneriert |

### `ExponentialMap`

Exponentialabbildung zwischen Lie-Algebra und Lie-Gruppe.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `exp_map` | `(X: np.ndarray) -> np.ndarray` | Matrixexponential $e^X$ via scipy |
| `log_map` | `(A: np.ndarray) -> np.ndarray` | Matrixlogarithmus $\log A$ |
| `baker_campbell_hausdorff` | `(X, Y, order: int = 4) -> np.ndarray` | BCH-Formel bis zur gewünschten Ordnung |

---

## Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `one_parameter_subgroup` | `(X: np.ndarray, t: float) -> np.ndarray` | Einparametrige Untergruppe $\exp(tX)$ |
| `classify_simple_lie_algebra` | `(n: int, type_char: str) -> dict` | Cartan-Klassifikation (Typ $A_n, B_n, C_n, D_n$) |
| `exceptional_lie_algebras_info` | `() -> dict` | Informationen zu $G_2, F_4, E_6, E_7, E_8$ |
| `cartan_matrix` | `(type_char: str, n: int) -> np.ndarray` | Cartan-Matrix für klassische Typen |
| `dynkin_diagram_info` | `(type_char: str, n: int) -> dict` | Dynkin-Diagramm Beschreibung |
| `fundamental_representation_su2` | `(j: float) -> dict` | Spin-$j$ Darstellung von $SU(2)$ |
| `casimir_element` | `(basis_matrices: list) -> np.ndarray` | Quadratischer Casimir-Operator |

---

## Beispiele

```python
import numpy as np
from lie_groups import SO, SU, LieAlgebra, ExponentialMap, classify_simple_lie_algebra

# SO(3): Drehgruppe im 3D-Raum
so3 = SO(3)
print("dim SO(3):", so3.dimension())          # 3
print("kompakt:", so3.is_compact())            # True

R = so3.rotation_3d_z(np.pi / 4)             # 45°-Drehung um z-Achse
print("Element von SO(3):", so3.is_element(R))# True

# SU(2): Quantenmechanische Drehgruppe
su2 = SU(2)
print("dim SU(2):", su2.dimension())          # 3
U = su2.su2_from_angles(0.5, 1.0, 0.3)
print("Element von SU(2):", su2.is_element(U))# True

# Lie-Algebra so(3): schiefsymmetrische Matrizen
L1 = np.array([[0,-1,0],[1,0,0],[0,0,0]], dtype=float)
L2 = np.array([[0,0,-1],[0,0,0],[1,0,0]], dtype=float)
L3 = np.array([[0,0,0],[0,0,-1],[0,1,0]], dtype=float)
algebra = LieAlgebra([L1, L2, L3])
print("[L1, L2]:", algebra.bracket(L1, L2))    # = L3
print("Halbeinfach:", algebra.is_semisimple())  # True

# Exponentialabbildung
em = ExponentialMap()
X = np.pi/2 * L3
R_exp = em.exp_map(X)                          # = Rotation um z um 90°
print("exp(π/2 · L3):", R_exp)

# Cartan-Klassifikation
info = classify_simple_lie_algebra(3, 'A')
print("A_3 = sl(4):", info)
```

---

## Technische Hinweise

- **Matrixexponential** wird via `scipy.linalg.expm` berechnet (Padé-Approximation).
- **Matrixlogarithmus** via `scipy.linalg.logm` (nur für invertierbare Matrizen nahe $I$ stabil).
- Die **BCH-Formel** wird als endliche Summe bis zur angegebenen Ordnung berechnet.
- Die **Killing-Form** wird über die adjungierte Darstellung $\operatorname{ad}_X$ berechnet.
- **Cartan-Matrizen** werden für die Typen $A_n, B_n, C_n, D_n$ explizit konstruiert.
