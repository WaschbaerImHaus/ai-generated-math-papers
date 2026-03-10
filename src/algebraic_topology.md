# Modul: Algebraische Topologie (`algebraic_topology.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Dieses Modul implementiert grundlegende und fortgeschrittene Konzepte der algebraischen Topologie:
simpliziale und singuläre Homologie, Kohomologie, Homotopiegruppen, Faserbündel,
Spektralsequenzen, K-Theorie und CW-Komplexe. Es dient sowohl als Lernwerkzeug als auch
zur algorithmischen Berechnung topologischer Invarianten.

## Mathematischer Hintergrund

### Simplizialkomplex und Homologie

Ein Simplizialkomplex $K$ besteht aus Simplizes (Punkte, Kanten, Dreiecke, ...).
Der **Randoperator** $\partial_k$ bildet $k$-Simplizes auf $(k-1)$-Simplizes ab:

$$\partial_k [v_0, \ldots, v_k] = \sum_{i=0}^{k} (-1)^i [v_0, \ldots, \hat{v}_i, \ldots, v_k]$$

Die $k$-te **Homologiegruppe** misst $k$-dimensionale Löcher:

$$H_k(K) = \frac{\ker(\partial_k)}{\operatorname{im}(\partial_{k+1})}$$

### Betti-Zahlen und Euler-Charakteristik

Die **Betti-Zahlen** $\beta_k = \operatorname{rank}(H_k)$ sind topologische Invarianten.
Die **Euler-Charakteristik** verbindet sie kombinatorisch:

$$\chi(X) = \sum_{k=0}^{\dim X} (-1)^k \cdot \beta_k = \sum_{k=0}^{\dim X} (-1)^k \cdot C_k$$

wobei $C_k$ die Anzahl der $k$-Zellen ist.

### Smith-Normalform

Zur Berechnung der Homologiegruppen wird die ganzzahlige Matrix $\partial_k$ in Smith-Normalform
transformiert: $D = U A V$ mit Diagonaleinträgen $d_1 \mid d_2 \mid \cdots \mid d_r$:

$$H_k \cong \mathbb{Z}^{\beta_k} \oplus \bigoplus_i \mathbb{Z}/d_i\mathbb{Z}$$

### Kohomologie und Cup-Produkt

Die Kohomologiegruppen $H^k(X; R)$ sind dual zur Homologie.
Das **Cup-Produkt** $\smile: H^p \otimes H^q \to H^{p+q}$ macht $H^*(X; R)$ zu einem Ring.
**Poincaré-Dualität** für geschlossene orientierbare $n$-Mannigfaltigkeiten:

$$H^k(M; \mathbb{Z}) \cong H_{n-k}(M; \mathbb{Z})$$

### Homotopiegruppen

Die $n$-te Homotopiegruppe $\pi_n(X, x_0)$ besteht aus Homotopieklassen von Abbildungen $S^n \to X$:

$$\pi_1(S^1) \cong \mathbb{Z}, \quad \pi_n(S^n) \cong \mathbb{Z} \text{ für } n \geq 1, \quad \pi_3(S^2) \cong \mathbb{Z}$$

### Hopf-Faserung

Die Hopf-Faserung ist ein Faserbündel $S^1 \to S^3 \to S^2$ mit der langen exakten Sequenz:

$$\cdots \to \pi_n(S^1) \to \pi_n(S^3) \to \pi_n(S^2) \to \pi_{n-1}(S^1) \to \cdots$$

### Bott-Periodizität (K-Theorie)

Die komplexe K-Theorie hat Periode 2 (Bott-Periodizität):

$$\tilde{K}(S^n) \cong \begin{cases} \mathbb{Z} & n \text{ gerade} \\ 0 & n \text{ ungerade} \end{cases}$$

---

## Hilfsfunktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `compute_smith_normal_form` | `(A: np.ndarray) -> (D, U, V)` | Smith-Normalform einer ganzzahligen Matrix |
| `homology_from_boundary_matrices` | `(boundary_ops: dict) -> dict` | Berechnet Homologiegruppen aus Randmatrizen |
| `classify_surface` | `(genus: int, orientable: bool) -> dict` | Klassifiziert Flächen nach Geschlecht und Orientierbarkeit |
| `van_kampen_free_product` | `(g1: str, g2: str) -> str` | Freies Produkt zweier Gruppen (Seifert-van-Kampen) |
| `lyndon_hochschild_serre_demo` | `() -> str` | Demonstriert Lyndon-Hochschild-Serre-Spektralsequenz |
| `classifying_space_demo` | `(group: str) -> str` | Beschreibt klassifizierenden Raum $BG$ einer Gruppe |

---

## Klassen und Methoden

### `SimplicialComplex`

Simplizialkomplex aus einer Liste von Simplizes mit automatischer Hinzunahme aller Seiten.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(simplices: List[tuple])` | Erzeugt Komplex; fügt alle Teilsimplizes (Seiten) hinzu |
| `_add_with_faces` | `(simplex: tuple)` | Fügt Simplex mit allen Seiten rekursiv ein |
| `faces` | `() -> List[tuple]` | Gibt alle Simplizes aller Dimensionen zurück |
| `_simplices_of_dim` | `(k: int) -> List[tuple]` | Gibt alle $k$-dimensionalen Simplizes zurück |
| `euler_characteristic` | `() -> int` | Berechnet $\chi = \sum (-1)^k C_k$ |
| `boundary_matrix` | `(p: int) -> np.ndarray` | Randmatrix $\partial_p$ als ganzzahlige Matrix |
| `betti_numbers` | `() -> Dict[int, int]` | Berechnet Betti-Zahlen $\beta_k$ via Smith-Normalform |

### `SimplicialHomology`

Berechnet simpliziale Homologiegruppen eines `SimplicialComplex`.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(complex: SimplicialComplex)` | Bindet Simplizialkomplex ein |
| `chain_groups` | `() -> Dict[int, int]` | Dimension der Kettengruppen $C_k$ |
| `boundary_operators` | `() -> Dict[int, np.ndarray]` | Alle Randmatrizen $\partial_k$ |
| `homology_groups` | `() -> Dict[int, Dict]` | Homologiegruppen mit Torsion und freiem Rang |
| `betti_numbers` | `() -> Dict[int, int]` | Betti-Zahlen $\beta_k$ |

### `SingularHomology`

Singuläre Homologie für klassische topologische Räume (analytisch/tabellenbasiert).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(space_name: str = "")` | Erzeugt singuläre Homologie-Berechner |
| `homology_of_sphere` | `(n: int) -> Dict[int, str]` | Homologie der $n$-Sphäre $S^n$ |
| `homology_of_torus` | `() -> Dict[int, str]` | Homologie des Torus $T^2$ |
| `homology_of_rp` | `(n: int) -> Dict[int, str]` | Homologie des projektiven Raums $\mathbb{RP}^n$ |
| `homology_of_cp` | `(n: int) -> Dict[int, str]` | Homologie des komplexen projektiven Raums $\mathbb{CP}^n$ |
| `mayer_vietoris_demo` | `() -> str` | Demonstriert die Mayer-Vietoris-Sequenz |

### `CohomologyRing`

Kohomologiering mit Cup-Produkt und Poincaré-Dualität.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(space_name: str = "")` | Erzeugt Kohomologie-Berechner |
| `cohomology_of_sphere` | `(n: int) -> Dict[int, str]` | Kohomologie der $n$-Sphäre |
| `cohomology_of_torus` | `() -> Dict[int, str]` | Kohomologie des Torus |
| `cup_product_demo` | `(space: str) -> str` | Demonstriert Cup-Produkt-Struktur |
| `poincare_duality_demo` | `(space: str, dim: int) -> str` | Demonstriert Poincaré-Dualität |

### `HomotopyGroups`

Homotopiegruppen, Seifert-van-Kampen-Theorem und Überlagerungsräume.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `()` | Erzeugt Homotopie-Berechner |
| `fundamental_group` | `(space: str) -> str` | Fundamentalgruppe $\pi_1$ bekannter Räume |
| `higher_homotopy_sphere` | `(n: int, k: int) -> str` | Homotopiegruppe $\pi_k(S^n)$ |
| `seifert_van_kampen_demo` | `() -> str` | Demonstriert Seifert-van-Kampen-Theorem |
| `covering_space_demo` | `() -> str` | Überlagerungsraum-Theorie |
| `long_exact_sequence_demo` | `() -> str` | Lange exakte Sequenz einer Faserung |

### `FiberBundle`

Faserbündel $F \to E \to B$ mit Hopf-Faserung und charakteristischen Klassen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(base: str, fiber: str, total: str)` | Definiert Faserbündel durch Basis, Faser, Totalraum |
| `hopf_fibration_demo` | `() -> Dict` | Hopf-Faserung $S^1 \hookrightarrow S^3 \twoheadrightarrow S^2$ |
| `vector_bundle_demo` | `() -> str` | Vektorbündel-Beispiele (Möbius-Band, Tangentialbündel) |
| `characteristic_classes_demo` | `() -> str` | Charakteristische Klassen des Bündels |
| `euler_class_demo` | `() -> str` | Euler-Klasse des Bündels |

### `SpectralSequence`

Spektralsequenzen (Serre, Leray-Hirsch).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(name: str)` | Erzeugt Spektralsequenz mit Name |
| `serre_spectral_sequence_demo` | `() -> str` | Serre-Spektralsequenz für Faserungen |
| `leray_hirsch_theorem_demo` | `() -> str` | Leray-Hirsch-Theorem |
| `e2_page_description` | `() -> str` | Beschreibung der $E_2$-Seite |

### `KTheory`

Komplexe K-Theorie mit Bott-Periodizität und Atiyah-Singer-Indexsatz.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `()` | Erzeugt K-Theorie-Berechner |
| `k_group_of_sphere` | `(n: int) -> Dict` | $\tilde{K}(S^n)$ via Bott-Periodizität |
| `bott_periodicity_demo` | `() -> str` | Demonstriert Bott-Periodizität |
| `chern_character_demo` | `() -> str` | Chern-Charakter $\mathrm{ch}: K(X) \to H^*(X;\mathbb{Q})$ |
| `atiyah_singer_index_theorem_demo` | `() -> str` | Atiyah-Singer-Indexsatz Demonstration |

### `CWComplex`

CW-Komplex mit zellulärer Homologie.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(cells: Dict[int, int])` | Erzeugt CW-Komplex aus Zellanzahlen pro Dimension |
| `euler_characteristic` | `() -> int` | Euler-Charakteristik $\chi = \sum (-1)^k C_k$ |
| `attaching_map_example` | `() -> str` | Beispiel für Anheftungsabbildungen |
| `cellular_homology_demo` | `() -> Dict[int, int]` | Zelluläre Homologiegruppen (vereinfacht) |

---

## Beispiele

```python
from algebraic_topology import SimplicialComplex, SimplicialHomology, SingularHomology, KTheory

# Torus als Simplizialkomplex (minimal: 7 Ecken, 21 Kanten, 14 Dreiecke)
# Vereinfachtes Beispiel: Dreieck (topologisch = D²)
K = SimplicialComplex([(0,1,2)])
print("Euler-Charakteristik:", K.euler_characteristic())  # 1
print("Betti-Zahlen:", K.betti_numbers())  # {0: 1, 1: 0, 2: 0}

# Homologie der Sphären
sh = SingularHomology()
print("H_k(S^2):", sh.homology_of_sphere(2))
# {0: 'Z', 1: '0', 2: 'Z'}

print("H_k(RP^3):", sh.homology_of_rp(3))
# {0: 'Z', 1: 'Z/2Z', 2: '0', 3: 'Z'}

# K-Theorie
kt = KTheory()
print("K(S^4):", kt.k_group_of_sphere(4))   # {'K0_tilde': 'Z', 'K1_tilde': '0'}
print("K(S^3):", kt.k_group_of_sphere(3))   # {'K0_tilde': '0', 'K1_tilde': 'Z'}
```

---

## Technische Hinweise

- Die **Smith-Normalform** wird durch ganzzahlige Zeilenoperationen berechnet (ohne externe Bibliotheken).
- **Singuläre Homologie** bekannter Räume wird über analytisch bekannte Ergebnisse zurückgegeben.
- **Homotopiegruppen** $\pi_k(S^n)$ werden tabellenbasiert für $k \leq n+3$ zurückgegeben.
- Die K-Theorie nutzt **Bott-Periodizität** mit Periode 2 für $\tilde{K}(S^n)$.
