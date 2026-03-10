# multilinear_algebra.py — Multilineare Algebra

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Datei:** `src/multilinear_algebra.py`

---

## Übersicht

Dieses Modul implementiert die fundamentalen Strukturen der **multilinearen Algebra**:

1. **Tensorprodukt** — Rang-1-Tensoren, Kronecker-Produkt, allgemeine Tensorklasse
2. **Äußere Algebra** — Wedge-Produkt, äußere Potenzen, Hodge-Operator
3. **Gram-Algebra** — Gram-Matrix, Gram-Determinante (Volumen)
4. **Symmetrische Algebra** — Symmetrische Potenzen, Polarisierung
5. **Multilineare Formen** — Symmetrie-/Alternierend-Tests, Signatur (Sylvester)

---

## 1. Tensorprodukt

### 1.1 Tensorprodukt zweier Vektoren

Das **äußere Produkt** (Tensorprodukt) zweier Vektoren $u \in \mathbb{R}^n$, $v \in \mathbb{R}^m$ ist:

$$
(u \otimes v)_{ij} = u_i \cdot v_j, \quad i = 1,\ldots,n,\; j = 1,\ldots,m
$$

Das Ergebnis ist eine $n \times m$-Matrix vom Rang 1.

```python
result = tensor_product_vectors([1, 2], [3, 4])
# [[3, 4], [6, 8]]
```

### 1.2 Kronecker-Produkt (Tensorprodukt von Matrizen)

Das **Kronecker-Produkt** $A \otimes B$ für $A \in \mathbb{R}^{m \times n}$, $B \in \mathbb{R}^{p \times q}$:

$$
A \otimes B = \begin{pmatrix}
a_{11} B & a_{12} B & \cdots \\
a_{21} B & a_{22} B & \cdots \\
\vdots & & \ddots
\end{pmatrix} \in \mathbb{R}^{mp \times nq}
$$

**Eigenschaften:**
- $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$
- $(A \otimes B)^T = A^T \otimes B^T$
- $\det(A \otimes B) = \det(A)^q \cdot \det(B)^p$ für quadratische Matrizen

### 1.3 Klasse `Tensor`

Allgemeiner $(p,q)$-Tensor mit $p$ kovarianten und $q$ kontravarianten Indizes:

$$
T^{i_1 \ldots i_q}_{j_1 \ldots j_p}
$$

Intern als `numpy`-Array der Form `[n₁,...,nₚ,m₁,...,mq]` gespeichert.

**Operationen:**

| Methode | Formel | Beschreibung |
|---------|--------|--------------|
| `contract(i, j)` | $T^{\ldots a \ldots}_{\ldots a \ldots}$ | Kontraktion (Spurbildung) |
| `symmetrize()` | $T_{(i_1\ldots i_n)} = \frac{1}{n!} \sum_\sigma T_{\sigma(i_1)\ldots\sigma(i_n)}$ | Symmetrisierung |
| `antisymmetrize()` | $T_{[i_1\ldots i_n]} = \frac{1}{n!} \sum_\sigma \text{sgn}(\sigma) T_{\sigma(i_1)\ldots}$ | Antisymmetrisierung |
| `__add__` | $S + T$ | Addition gleicher Typen |
| `__mul__` | $\lambda \cdot T$ | Skalarmultiplikation |

### 1.4 Tensorkontraktion

Für einen Tensor $T$ mit Rang $\geq 2$, Kontraktion über Achsen $i$ und $j$:

$$
(\text{tr}_{ij} T)_{\ldots} = \sum_a T_{\ldots a \ldots a \ldots}
$$

Implementiert via `numpy.trace(T, axis1=i, axis2=j)`.

---

## 2. Äußere (Graßmann-)Algebra

### 2.1 Wedge-Produkt zweier Vektoren

Das **äußere Produkt** $u \wedge v$ für $u, v \in \mathbb{R}^n$ ist ein Bivector mit Koeffizienten:

$$
(u \wedge v)_{ij} = u_i v_j - u_j v_i, \quad i < j
$$

Der Ergebnisvektor hat Länge $\binom{n}{2}$ und liegt im Raum $\Lambda^2(\mathbb{R}^n)$.

**Eigenschaften:**
- **Antisymmetrie:** $u \wedge v = -(v \wedge u)$
- **Nilpotenz:** $u \wedge u = 0$
- **Linearität:** $(au + bw) \wedge v = a(u \wedge v) + b(w \wedge v)$

```python
wedge_product([1,0,0], [0,1,0])
# [1.0, 0.0, 0.0]  →  e₁∧e₂ Koeffizient ist 1
```

### 2.2 Äußere Potenz $k$-ter Grad

Die **$k$-te äußere Potenz** $u_1 \wedge \cdots \wedge u_k$ berechnet sich durch Unterdeterminanten:

$$
(u_1 \wedge \cdots \wedge u_k)_{i_1\ldots i_k} = \det \begin{pmatrix} u_1[i_1] & \cdots & u_k[i_1] \\ \vdots & & \vdots \\ u_1[i_k] & \cdots & u_k[i_k] \end{pmatrix}
$$

Spezialfall $k = n$: Liefert die **Determinante** als einzigen Koeffizienten.

### 2.3 Klasse `ExteriorAlgebra`

Die **äußere Algebra** $\Lambda(V) = \bigoplus_{k=0}^n \Lambda^k(V)$ für $V = \mathbb{R}^n$:

$$
\dim(\Lambda^k(V)) = \binom{n}{k}, \quad \dim(\Lambda(V)) = 2^n
$$

k-Formen werden als Dictionaries `{(i₁,...,iₖ): Koeffizient}` mit $i_1 < \cdots < i_k$ dargestellt.

**Wedge-Produkt zweier Formen:**

$$
(\alpha \wedge \beta)(v_1,\ldots,v_{p+q}) = \frac{1}{p!\,q!} \sum_{\sigma \in S_{p+q}} \text{sgn}(\sigma)\, \alpha(v_{\sigma(1)},\ldots) \beta(v_{\sigma(p+1)},\ldots)
$$

**Hodge-Dual:**

Der Hodge-Sternoperator $\star: \Lambda^k(V) \to \Lambda^{n-k}(V)$ wirkt auf Basisformen:

$$
\star(e_{i_1} \wedge \cdots \wedge e_{i_k}) = \text{sgn}(\sigma)\, \sqrt{|\det g|}\, e_{j_1} \wedge \cdots \wedge e_{j_{n-k}}
$$

wobei $\{i_1,\ldots,i_k, j_1,\ldots,j_{n-k}\} = \{0,\ldots,n-1\}$ und $\sigma$ die Permutation zur Sortierung ist.

---

## 3. Gram-Matrix und Gram-Determinante

### 3.1 Gram-Matrix

Für Vektoren $v_1, \ldots, v_k \in \mathbb{R}^n$ ist die **Gram-Matrix**:

$$
G_{ij} = \langle v_i, v_j \rangle = \sum_l (v_i)_l (v_j)_l
$$

$G$ ist positiv semidefinit und positiv definit genau dann, wenn die Vektoren **linear unabhängig** sind.

### 3.2 Gram-Determinante

$$
\det(G) = |v_1 \wedge \cdots \wedge v_k|^2 \geq 0
$$

Das ist das **Quadrat des Volumens** des von den Vektoren aufgespannten Parallelotops.

Spezialfall $k = n$: $\det(G) = \det(V)^2$ für die Matrix $V$ mit Vektoren als Zeilen.

---

## 4. Symmetrische Algebra

### 4.1 Dimensionsformel

Der Raum $S^k(V)$ der symmetrischen $k$-Tensoren hat Dimension:

$$
\dim(S^k(V)) = \binom{n+k-1}{k} = \frac{(n+k-1)!}{k!\,(n-1)!}
$$

(Anzahl der Monomiale vom Grad $k$ in $n$ Variablen)

| $n$ \ $k$ | 0 | 1 | 2 | 3 | 4 |
|-----------|---|---|---|---|---|
| 2 | 1 | 2 | 3 | 4 | 5 |
| 3 | 1 | 3 | 6 | 10 | 15 |
| 4 | 1 | 4 | 10 | 20 | 35 |

### 4.2 Polarisierung

Die **Polarisierung** eines homogenen Polynoms $p(x)$ vom Grad $k$ liefert eine symmetrische $k$-lineare Form $P(v_1,\ldots,v_k)$ mit $P(v,\ldots,v) = p(v)$:

$$
P(v_1,\ldots,v_k) = \frac{1}{k!} \sum_{S \subseteq \{1,\ldots,k\}} (-1)^{k-|S|} p\!\left(\sum_{i \in S} v_i\right)
$$

Implementiert als Skalierung der Multiindex-Koeffizienten.

---

## 5. Multilineare Formen

### 5.1 Alternierend / Antisymmetrisch

Ein Tensor $T$ ist **alternierend**, wenn:

$$
T[\ldots, i, \ldots, j, \ldots] = -T[\ldots, j, \ldots, i, \ldots] \quad \forall\, i \neq j
$$

### 5.2 Symmetrisch

Ein Tensor $T$ ist **symmetrisch**, wenn:

$$
T[\ldots, i, \ldots, j, \ldots] = T[\ldots, j, \ldots, i, \ldots] \quad \forall\, i, j
$$

### 5.3 Sylvester'scher Trägheitssatz (Signatur)

Für eine symmetrische Bilinearform $B: V \times V \to \mathbb{R}$ mit darstellender Matrix $M$ ist die **Signatur** $(p, q)$ eine Invariante:

$$
p = \#\{\lambda_i > 0\}, \quad q = \#\{\lambda_i < 0\}
$$

wobei $\lambda_i$ die Eigenwerte von $M$ sind. Beispiele:

| Form | Signatur | Bedeutung |
|------|----------|-----------|
| Euklidisches Skalarprodukt | $(n, 0)$ | Positiv definit |
| Minkowski-Metrik | $(1, 3)$ | Pseudoeuklidisch |
| Lorentz-Metrik | $(3, 1)$ | Alternat. Konvention |
| Indifferente Form | $(p, q)$ mit $p,q>0$ | Indefinit |

---

## Verwendete Bibliotheken

| Bibliothek | Verwendung |
|------------|-----------|
| `numpy` | Tensoroperationen, Eigenwerte, Matrixalgebra |
| `itertools` | Permutationen, Kombinationen |
| `math` | Faktorielle, Binomialkoeffizienten |

---

## Abhängigkeiten

- Keine internen Projektabhängigkeiten
- Rein standalone nutzbar

---

## Tests

Alle Tests in `tests/test_multilinear_algebra.py`. Abdeckung: 35+ Tests.
