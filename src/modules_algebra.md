# modules_algebra.py — Moduln über Ringen

## Überblick

Dieses Modul implementiert die Theorie der **R-Moduln über Ringen**, die
Verallgemeinerung von Vektorräumen auf Ringe statt Körper.

Der wichtigste Anwendungsfall sind **endlich erzeugte abelsche Gruppen**,
die als ℤ-Moduln der Form $M = \mathbb{Z}^n / \operatorname{Im}(A)$
für eine ganzzahlige Matrix $A$ dargestellt werden.

---

## Mathematischer Hintergrund

### R-Modul

Ein **R-Modul** $M$ über einem Ring $R$ ist eine abelsche Gruppe $(M, +)$
zusammen mit einer Skalarmultiplikation:
$$R \times M \to M, \quad (r, m) \mapsto r \cdot m$$
die folgende Axiome erfüllt:
- $r(m_1 + m_2) = rm_1 + rm_2$ (Distributivgesetz)
- $(r_1 + r_2)m = r_1 m + r_2 m$
- $(r_1 r_2)m = r_1(r_2 m)$
- $1_R \cdot m = m$

Über $R = \mathbb{Z}$ sind Moduln dasselbe wie abelsche Gruppen.

---

## Funktionen

### `smith_normal_form(matrix)`

Berechnet die **Smith-Normalform** einer ganzzahligen Matrix:
$$A = U \cdot D \cdot V, \quad U, V \in \mathrm{GL}_n(\mathbb{Z}), \quad D = \mathrm{diag}(d_1, \ldots, d_r, 0, \ldots, 0)$$
mit der **Teilbarkeitskette** $d_1 \mid d_2 \mid \cdots \mid d_r$.

**Klassifikationssatz:**
$$\mathbb{Z}^n / \operatorname{Im}(A) \cong \mathbb{Z}_{d_1} \times \cdots \times \mathbb{Z}_{d_r} \times \mathbb{Z}^{n-r}$$

**Algorithmus:** Gauß-Elimination mit elementaren ganzzahligen Zeilen- und
Spaltenoperationen. In jedem Schritt wird das betragskleinste Nicht-Null-Element
als Pivot gewählt.

**Rückgabe:**
```python
{
    'invariant_factors': [2, 4],   # d_1 | d_2 | ...
    'rank': 2,                     # Rang der Matrix
    'U': [...],                    # linke Transformationsmatrix
    'D': [...],                    # Diagonalmatrix
    'V': [...],                    # rechte Transformationsmatrix
}
```

**Beispiel:**
```python
smith_normal_form([[2, 4], [6, 8]])
# → invariant_factors = [2, 4]
```

---

### `module_from_matrix(A)`

Konstruiert den $\mathbb{Z}$-Modul $M = \mathbb{Z}^n / \operatorname{Im}(A)$:

```python
module_from_matrix([[6]])
# → torsion_part = 'ℤ_6', free_part = '0'
```

---

### `structure_theorem_abelian_groups(n)`

**Klassifikationssatz endlich erzeugter abelscher Gruppen:**

Jede endliche abelsche Gruppe der Ordnung $n$ ist isomorph zu einem
Produkt zyklischer Gruppen:

$$G \cong \mathbb{Z}_{d_1} \times \cdots \times \mathbb{Z}_{d_k}$$
mit $d_1 \mid d_2 \mid \cdots \mid d_k$ und $\prod d_i = n$.

Die Anzahl nicht-isomorpher abelscher Gruppen der Ordnung $n = p_1^{a_1} \cdots p_r^{a_r}$
ist $p(a_1) \cdot p(a_2) \cdots p(a_r)$, wobei $p(k)$ die Partitionszahl ist.

**Beispiele:**
- $n = 4$: $\mathbb{Z}_4$, $\mathbb{Z}_2 \times \mathbb{Z}_2$
- $n = 12$: $\mathbb{Z}_{12}$, $\mathbb{Z}_2 \times \mathbb{Z}_6$
- $n = 8$: $\mathbb{Z}_8$, $\mathbb{Z}_4 \times \mathbb{Z}_2$, $\mathbb{Z}_2^3$

---

### `tensor_product_modules(M_gens, N_gens, R_mod)`

**Tensorprodukt** $M \otimes_\mathbb{Z} N$ für zyklische $\mathbb{Z}$-Moduln:

$$\mathbb{Z}_m \otimes_\mathbb{Z} \mathbb{Z}_n \cong \mathbb{Z}_{\gcd(m,n)}$$

**Beispiel:**
```python
tensor_product_modules([4], [6], 0)
# → order = 2  (ℤ_4 ⊗ ℤ_6 ≅ ℤ_2)
```

---

### `hom_module(M_order, N_order)`

**Hom-Modul** $\operatorname{Hom}_\mathbb{Z}(\mathbb{Z}_m, \mathbb{Z}_n)$:

$$|\operatorname{Hom}(\mathbb{Z}_m, \mathbb{Z}_n)| = \gcd(m, n)$$

Ein Homomorphismus $f: \mathbb{Z}_m \to \mathbb{Z}_n$ ist durch $f(1)$ bestimmt,
und es muss $m \cdot f(1) \equiv 0 \pmod{n}$ gelten.

**Beispiel:**
```python
hom_module(4, 6)
# → order = 2  (|Hom(ℤ_4, ℤ_6)| = gcd(4,6) = 2)
```

---

### `exact_sequence_check(maps, modules)`

Prüft ob eine kurze exakte Sequenz $0 \to A \xrightarrow{f} B \xrightarrow{g} C \to 0$
**exakt** ist, d.h. $\operatorname{Im}(f) = \ker(g)$ an jedem Punkt.

**Standard-Beispiel:**
$$0 \to \mathbb{Z} \xrightarrow{\cdot n} \mathbb{Z} \to \mathbb{Z}_n \to 0 \quad \text{(exakt)}$$

---

### `free_resolution(module_matrix)`

Berechnet eine **freie Auflösung** des Moduls $M = \mathbb{Z}^n / \operatorname{Im}(A)$:

$$0 \to \mathbb{Z}^m \xrightarrow{A} \mathbb{Z}^n \xrightarrow{\varepsilon} M \to 0$$

Die Länge der Auflösung ist 1 für endlich erzeugte $\mathbb{Z}$-Moduln
(da $\mathbb{Z}$ globale Dimension 1 hat).

---

## Klasse `Module`

```python
m = Module([[2, 0], [0, 3]])
m.rank()                # Freier Rang
m.torsion_submodule()   # Torsionsanteil
m.free_part()           # Freier Anteil
m.is_free()             # Keine Torsion?
m.is_finitely_generated()  # Immer True
```

---

## Tests

Die Testsuite (`tests/test_modules_algebra.py`) umfasst **42 Tests**:

| Testklasse | Tests | Abdeckung |
|-----------|-------|-----------|
| `TestSmithNormalForm` | 10 | SNF-Algorithmus |
| `TestModuleFromMatrix` | 3 | Modul-Konstruktion |
| `TestStructureTheoremAbelianGroups` | 9 | Klassifikation |
| `TestTensorProduct` | 4 | Tensorprodukt |
| `TestHomModule` | 5 | Hom-Modul |
| `TestExactSequence` | 3 | Exakte Sequenzen |
| `TestFreeResolution` | 2 | Freie Auflösung |
| `TestModuleClass` | 6 | Module-Klasse |

---

## Autor

**Kurt Ingwer** — @lastModified 2026-03-10
