# Modul: Kombinatorik (`combinatorics.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `combinatorics.py` implementiert Kombinatorik und Diskrete Mathematik in sieben thematischen Klassen. Es umfasst Permutationsgruppen, abzählende Kombinatorik (Stirling-, Euler-, Catalan-Zahlen), erzeugende Funktionen, Ramsey-Theorie, Graphen-Kombinatorik (chromatisches Polynom, Kirchhoff-Baum-Satz) und Kombinatorik auf Wörtern (Lyndon-Wörter, de-Bruijn-Folgen, Bearbeitungsabstand).

---

## Mathematischer Hintergrund

### Stirling-Zahlen

**Stirling-Zahlen 1. Art** $s(n,k)$ (vorzeichenbehaftet) zählen Permutationen von $\{1,\ldots,n\}$ mit genau $k$ Zyklen:
$$s(n,k) = s(n-1,k-1) - (n-1)\cdot s(n-1,k), \quad s(0,0) = 1$$

**Stirling-Zahlen 2. Art** $S(n,k)$ zählen Partitionen einer $n$-Menge in $k$ nicht-leere Blöcke:
$$S(n,k) = k \cdot S(n-1,k) + S(n-1,k-1), \quad S(0,0) = 1$$

### Euler-Zahlen

$A(n,k)$ = Anzahl der Permutationen von $\{1,\ldots,n\}$ mit genau $k$ Aufstiegen (Positionen $i$ mit $\pi(i) < \pi(i+1)$):
$$A(n,k) = (k+1) \cdot A(n-1,k) + (n-k) \cdot A(n-1,k-1)$$

### Partitionen und spezielle Zahlen

**Partitionszahl** $p(n)$: Anzahl der Möglichkeiten, $n$ als ungeordnete Summe positiver ganzer Zahlen zu schreiben. Euler'sche erzeugende Funktion:
$$\sum_{n=0}^\infty p(n) x^n = \prod_{k=1}^\infty \frac{1}{1-x^k}$$

**Motzkin-Zahl** $M_n$: Anzahl der Wege auf dem Gitter $\mathbb{Z}^2$ von $(0,0)$ nach $(n,0)$ mit Schritten $(1,1)$, $(1,0)$, $(1,-1)$, die nie unter die $x$-Achse fallen.

**Narayana-Zahl:** $N(n,k) = \frac{1}{n}\binom{n}{k}\binom{n}{k-1}$, Verfeinerung der Catalan-Zahlen.

**Wahlproblem (Ballot Problem):**
$$P(\text{Kandidat A liegt stets vorn}) = \frac{p - q}{p + q}, \quad p > q$$

### Ramsey-Theorie

**Ramsey-Zahl** $R(s,t)$: Kleinste natürliche Zahl $n$, sodass jede 2-Färbung der Kanten von $K_n$ einen roten $K_s$ oder blauen $K_t$ enthält.

$$R(s,t) \leq R(s-1,t) + R(s,t-1)$$

Bekannte Werte: $R(3,3)=6$, $R(3,4)=9$, $R(3,5)=14$, $R(4,4)=18$.

### Kirchhoff-Matrix-Satz (Graphen-Kombinatorik)

**Anzahl der Spannbäume** eines Graphen $G$ mit Laplace-Matrix $L$:
$$\tau(G) = \det(\tilde{L})$$

wobei $\tilde{L}$ die um eine Zeile/Spalte reduzierte Laplace-Matrix ist.

### Chromatisches Polynom

Für einen Graphen $G$ mit $n$ Knoten und einer Kante $\{u,v\}$:
$$P(G, k) = P(G - e, k) - P(G / e, k)$$

(Lösch-Kontraktion; $P(K_n, k) = k(k-1)\cdots(k-n+1)$)

### Erzeugende Funktionen

**Gewöhnliche erzeugende Funktion (OGF):**
$$F(x) = \sum_{n=0}^\infty a_n x^n$$

**Exponential erzeugende Funktion (EGF):**
$$F(x) = \sum_{n=0}^\infty a_n \frac{x^n}{n!}$$

Multiplikation zweier OGF entspricht der Faltung der Folgen: $(f * g)_n = \sum_{k=0}^n f_k g_{n-k}$.

---

## Klassen und Methoden

### `PermutationGroup`

Symmetrische Gruppe $S_n$ – Theorie der Permutationen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `all_permutations` | `(n: int) → list[list[int]]` | Alle $n!$ Permutationen von $\{0,\ldots,n-1\}$ |
| `permutation_sign` | `(perm: list[int]) → int` | Vorzeichen: $(-1)^{\text{Transpositionsanzahl}}$ |
| `cycle_decomposition` | `(perm: list[int]) → list[list[int]]` | Zerlegung in disjunkte Zyklen |
| `order_of_permutation` | `(perm: list[int]) → int` | Ordnung = kgV der Zykellängen |
| `compose_permutations` | `(p1, p2: list[int]) → list[int]` | Komposition $p_1 \circ p_2$ |
| `inverse_permutation` | `(p: list[int]) → list[int]` | Inverse Permutation $p^{-1}$ |
| `is_even` | `(perm: list[int]) → bool` | Prüft ob Permutation gerade ist |
| `fixed_points` | `(perm: list[int]) → list[int]` | Fixpunkte $\{i : p(i) = i\}$ |

### `EnumerativeCombinatorics`

Abzählende Kombinatorik – spezielle Zahlenfolgen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `multinomial` | `(n: int, ks: list[int]) → int` | Multinomialkoeffizient $\binom{n}{k_1,\ldots,k_r}$ |
| `stirling_first` | `(n, k: int) → int` | Stirling-Zahlen 1. Art $s(n,k)$ |
| `stirling_second` | `(n, k: int) → int` | Stirling-Zahlen 2. Art $S(n,k)$ |
| `euler_number_a` | `(n, k: int) → int` | Euler-Zahlen $A(n,k)$ (Aufstiege) |
| `partition_number` | `(n: int) → int` | Partitionszahl $p(n)$ |
| `motzkin_number` | `(n: int) → int` | Motzkin-Zahl $M_n$ |
| `narayana_number` | `(n, k: int) → int` | Narayana-Zahl $N(n,k)$ |
| `ballot_problem` | `(p, q: int) → float` | Wahrscheinlichkeit beim Wahlproblem |
| `inclusion_exclusion` | `(sets: list[set]) → int` | Inklusion-Exklusion $|\bigcup A_i|$ |

### `GeneratingFunctions`

Erzeugende Funktionen und Potenzreihenoperationen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `ordinary_gf` | `(sequence, n_terms) → list[float]` | OGF: erste $n$ Koeffizienten |
| `exponential_gf` | `(sequence, n_terms) → list[float]` | EGF: erste $n$ Koeffizienten |
| `fibonacci_gf` | `() → str` | OGF der Fibonacci-Folge als String |
| `catalan_gf` | `() → str` | OGF der Catalan-Zahlen als String |
| `power_series_multiply` | `(f, g, n_terms) → list[float]` | Faltung zweier Potenzreihen |
| `gf_from_recurrence` | `(a0, a1, recurrence_fn, n_terms) → list[float]` | GF aus Rekurrenzrelation |

### `RamseyTheory`

Ramsey-Theorie und Extremalkombinatorik.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `ramsey_number_R` | `(s, t: int) → int \| None` | Obere Schranke für $R(s,t)$ |
| `ramsey_graph_coloring` | `(n, k: int) → dict` | Zufällige $k$-Färbung von $K_n$ |
| `van_der_waerden_numbers` | `() → dict` | Bekannte van-der-Waerden-Zahlen $W(k,r)$ |
| `hales_jewett_theorem_demo` | `(n, k: int) → dict` | Hales-Jewett-Demo für $[n]^k$ |
| `pigeonhole_principle` | `(n_pigeons, n_holes: int) → dict` | Schubfachprinzip-Illustration |

### `GraphCombinatorics`

Kombinatorische Graphentheorie.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `chromatic_polynomial` | `(adj_matrix) → list[int]` | Chromatisches Polynom $P(G, k)$ |
| `tutte_polynomial_demo` | `(adj_matrix) → dict` | Tutte-Polynom (Demo) |
| `matching_number` | `(adj_matrix) → int` | Maximale Matchinggröße |
| `perfect_matching_count` | `(adj_matrix) → int` | Anzahl perfekter Matchings |
| `spanning_tree_count` | `(adj_matrix) → int` | Kirchhoff-Matrix-Satz: $\tau(G)$ |
| `eulerian_circuit_exists` | `(adj_matrix) → bool` | Euler-Kriterium (alle Grade gerade) |
| `hamiltonian_path_check` | `(adj_matrix) → bool` | Hamiltonpfad per Backtracking |

### `WordCombinatorics`

Kombinatorik auf Wörtern und Zeichenketten.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `lyndon_word_check` | `(word: str) → bool` | Prüft ob Wort ein Lyndon-Wort ist |
| `necklaces` | `(n, k: int) → int` | Anzahl $k$-färbiger Halsketten der Länge $n$ |
| `primitive_word` | `(word: str) → bool` | Prüft ob Wort primitiv (nicht Potenz) |
| `longest_common_subsequence` | `(s1, s2: str) → int` | LCS-Länge per Dynamic Programming |
| `edit_distance` | `(s1, s2: str) → int` | Levenshtein-Editierdistanz |
| `de_bruijn_sequence` | `(n, k: int) → str` | de-Bruijn-Folge: alle $k^n$ Teilwörter einmal |

### Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `mat_mul` | `(A, B) → list[list[int]]` | Matrixmultiplikation für ganzzahlige Matrizen |
| `mat_pow` | `(M, exp) → list[list[int]]` | schnelle Matrixpotenz per Quadrierung |
| `is_happy` | `(num: int) → bool` | Glückliche Zahl (iterative Quersummenquadrate) |

---

## Wichtige kombinatorische Identitäten

**Bell-Zahlen** (Gesamtanzahl Partitionen einer $n$-Menge):
$$B_n = \sum_{k=0}^n S(n,k)$$

**Catalan-Zahlen:**
$$C_n = \frac{1}{n+1}\binom{2n}{n}, \quad C_0=1, C_1=1, C_2=2, C_3=5, \ldots$$

**Multinomialkoeffizient:**
$$\binom{n}{k_1, k_2, \ldots, k_r} = \frac{n!}{k_1!\, k_2!\, \cdots\, k_r!}, \quad \sum k_i = n$$

**Inklusion-Exklusion:**
$$\left|\bigcup_{i=1}^n A_i\right| = \sum|A_i| - \sum|A_i \cap A_j| + \ldots + (-1)^{n+1}|A_1 \cap \ldots \cap A_n|$$

---

## Beispiele

```python
from combinatorics import PermutationGroup, EnumerativeCombinatorics, GraphCombinatorics, WordCombinatorics

# -- Permutationen --
pg = PermutationGroup(4)
perms = PermutationGroup.all_permutations(3)   # 6 Permutationen
sign = PermutationGroup.permutation_sign([1, 0, 2, 3])  # -1 (ungerade)
cycles = PermutationGroup.cycle_decomposition([1, 2, 0, 3])  # [[0,1,2], [3]]
order = PermutationGroup.order_of_permutation([1, 2, 0, 3])  # kgV(3,1) = 3

# -- Stirling-Zahlen --
ec = EnumerativeCombinatorics()
print(ec.stirling_second(4, 2))  # S(4,2) = 7
print(ec.partition_number(10))   # p(10) = 42

# -- Graphen-Kombinatorik: Spannbäume --
gc = GraphCombinatorics()
K4 = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]
print(gc.spanning_tree_count(K4))   # 16 Spannbäume von K₄

# -- de-Bruijn-Folge --
wc = WordCombinatorics()
db = wc.de_bruijn_sequence(2, 3)    # Länge 2³ = 8, enthält alle 2-Bit Teilwörter
print(db)

# -- Editierdistanz --
print(wc.edit_distance("kitten", "sitting"))  # 3
```

---

## Tests

**Testdatei:** `tests/test_combinatorics.py`
**Abdeckung:** Alle Klassen, Rekurrenzrelationen, Graphen-Spezialfälle, Ramsey-Bounds, Wörter-Edge-Cases

---

## Implementierungshinweise

- **Stirling-Zahlen:** Werden per `@lru_cache` rekursiv berechnet – für große $n$ ist der Cache entscheidend.
- **Chromatisches Polynom:** Verwendet Lösch-Kontraktion rekursiv; Komplexität $O(2^m)$ für $m$ Kanten.
- **Kirchhoff-Satz:** Die Determinante der reduzierten Laplace-Matrix wird per Gauß-Elimination berechnet.
- **de-Bruijn-Folge:** Martin-Algorithmus mit rekursiver Konstruktion in $O(k^n)$ Zeit.
- **LCS und Editierdistanz:** Standard-DP mit $O(mn)$ Zeit- und Platzkomplexität.
