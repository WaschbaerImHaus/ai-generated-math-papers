# group_theory.py — Gruppentheorie-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Build:** 11

---

## Überblick

Das Modul implementiert die klassische endliche Gruppentheorie vollständig in Python.
Es stellt Klassen für alle wichtigen Gruppenstrukturen sowie freie Funktionen für
zentrale gruppentheoretische Sätze bereit.

---

## Mathematische Grundlagen

### Definition einer Gruppe

Eine **Gruppe** $(G, \cdot)$ ist eine Menge $G$ mit einer binären Operation $\cdot$, die:

1. **Abgeschlossenheit**: $\forall a, b \in G: a \cdot b \in G$
2. **Assoziativität**: $\forall a, b, c \in G: (a \cdot b) \cdot c = a \cdot (b \cdot c)$
3. **Neutrales Element**: $\exists e \in G: e \cdot a = a \cdot e = a$
4. **Inverses Element**: $\forall a \in G: \exists a^{-1}: a \cdot a^{-1} = a^{-1} \cdot a = e$

### Satz von Lagrange

Für jede Untergruppe $H \leq G$ gilt:
$$|G| = |H| \cdot [G:H]$$
wobei $[G:H]$ der **Index** von $H$ in $G$ ist.

### Sylow-Sätze

Sei $|G| = p^k \cdot m$ mit $\gcd(p, m) = 1$. Dann:
- **Existenz**: $G$ besitzt eine Sylow-$p$-Untergruppe der Ordnung $p^k$
- **Konjugiertheit**: Alle Sylow-$p$-Untergruppen sind zueinander konjugiert
- **Anzahl**: $n_p \equiv 1 \pmod{p}$ und $n_p \mid m$

### Burnside-Lemma

Für eine Gruppenoperation $G \times X \to X$:
$$|X/G| = \frac{1}{|G|} \sum_{g \in G} |\operatorname{Fix}(g)|$$
wobei $\operatorname{Fix}(g) = \{x \in X \mid g \cdot x = x\}$.

### Hauptsatz über endliche abelsche Gruppen

Jede endliche abelsche Gruppe $G$ der Ordnung $n$ ist isomorph zu:
$$G \cong \mathbb{Z}_{p_1^{a_1}} \times \mathbb{Z}_{p_2^{a_2}} \times \cdots \times \mathbb{Z}_{p_k^{a_k}}$$
wobei $n = p_1^{a_1} \cdot p_2^{a_2} \cdots p_k^{a_k}$.

---

## Klassen

### `Group`

Allgemeine endliche Gruppe, dargestellt durch Elementmenge und binäre Operation.

```python
G = Group(elements, operation, identity=None, name="G")
```

**Methoden:**

| Methode | Rückgabe | Beschreibung |
|---------|----------|--------------|
| `order()` | `int` | $\|G\|$ |
| `is_abelian()` | `bool` | Prüft Kommutativität |
| `is_cyclic()` | `bool` | $\exists g: \langle g \rangle = G$ |
| `element_order(a)` | `int` | $\min\{n \geq 1: a^n = e\}$ |
| `inverse(a)` | `Any` | $a^{-1}$ |
| `power(a, n)` | `Any` | $a^n$ (auch negativ) |
| `subgroups()` | `list[Subgroup]` | Alle Untergruppen $H \leq G$ |
| `normal_subgroups()` | `list[Subgroup]` | Alle Normalteiler $N \trianglelefteq G$ |
| `center()` | `Subgroup` | $Z(G) = \{z \mid \forall g: zg=gz\}$ |
| `commutator_subgroup()` | `Subgroup` | $[G,G] = \langle [a,b] \rangle$ |
| `cosets(H, left)` | `list[list]` | Links-/Rechtsnebenklassen |
| `cayley_table()` | `np.ndarray` | Multiplikationstabelle |

### `Subgroup(Group)`

Untergruppe $H \leq G$ mit Referenz auf Obergruppe.

```python
H = Subgroup(elements, operation, parent=G, identity=e)
```

Zusätzliche Methode: `index()` → $[G:H] = |G|/|H|$

### `GroupHomomorphism`

Strukturerhaltende Abbildung $\varphi: G \to H$ mit $\varphi(ab) = \varphi(a)\varphi(b)$.

```python
phi = GroupHomomorphism(source, target, mapping, name="φ")
```

**Methoden:**

| Methode | Rückgabe | Beschreibung |
|---------|----------|--------------|
| `kernel()` | `Subgroup` | $\ker(\varphi) = \{g \mid \varphi(g) = e_H\}$ |
| `image()` | `Subgroup` | $\operatorname{im}(\varphi) = \{\varphi(g)\}$ |
| `is_injective()` | `bool` | $\ker(\varphi) = \{e\}$ |
| `is_surjective()` | `bool` | $\operatorname{im}(\varphi) = H$ |
| `is_isomorphism()` | `bool` | bijektiver Homomorphismus |

### `QuotientGroup(Group)`

Faktorgruppe $G/N$ für Normalteiler $N \trianglelefteq G$.

```python
Q = QuotientGroup(G, N)
```

Elemente sind Frozensets (Nebenklassen $gN$).

### `PermutationGroup(Group)`

Untergruppe von $S_n$. Elemente sind Tupel $\sigma = (\sigma(0), \ldots, \sigma(n-1))$.

```python
PG = PermutationGroup(elements, n, name="PG")
```

**Methoden:**

| Methode | Rückgabe | Beschreibung |
|---------|----------|--------------|
| `cycle_decomposition(σ)` | `list[list[int]]` | Disjunkte Zykel |
| `cycle_type(σ)` | `tuple[int,...]` | Sortierte Zyklenlängen |
| `sign(σ)` | `int` | $\pm 1$ (Vorzeichen) |
| `is_even_permutation(σ)` | `bool` | $\text{sgn}(\sigma) = +1$ |
| `orbit(x)` | `set[int]` | $\text{Orb}_G(x)$ |
| `stabilizer(x)` | `PermutationGroup` | $\text{Stab}_G(x)$ |

---

## Freie Funktionen

### `cyclic_group(n)` → `Group`
$\mathbb{Z}_n = (\{0,\ldots,n-1\}, + \bmod n)$

### `dihedral_group(n)` → `Group`
$D_n$: Symmetriegruppe des regulären $n$-Ecks, Ordnung $2n$.
Elemente: $(k, f)$ mit $k \in \{0,\ldots,n-1\}$, $f \in \{0,1\}$
($f=0$: Drehung $r^k$, $f=1$: Spiegelung $sr^k$)

### `quaternion_group()` → `Group`
$Q_8 = \{\pm 1, \pm i, \pm j, \pm k\}$
Relationen: $i^2 = j^2 = k^2 = ijk = -1$

### `symmetric_group(n)` → `PermutationGroup`
$S_n$ aller Permutationen, Ordnung $n!$

### `alternating_group(n)` → `PermutationGroup`
$A_n$ gerader Permutationen, Ordnung $n!/2$

### `lagrange_theorem_check(G, H)` → `dict`
Überprüft $|G| = |H| \cdot [G:H]$ und berechnet alle Nebenklassen.

### `sylow_theorems(G, p)` → `dict`
Berechnet alle Sylow-$p$-Untergruppen und prüft die drei Sylow-Bedingungen.

### `classify_abelian_group(n)` → `dict`
Hauptsatz: Alle Isomorphieklassen abelscher Gruppen der Ordnung $n$.

### `is_simple_group(G)` → `bool`
Prüft ob $G$ einfach ist. Für kleine Gruppen direkt, für große via Konjugationsklassen.

### `group_action(G, X, action)` → `dict`
Orbits, Stabilisatoren, Burnside-Lemma.

### `direct_product(G, H)` → `Group`
$G \times H$ mit komponentenweiser Operation, Ordnung $|G| \cdot |H|$.

### `semidirect_product(N, H, phi)` → `Group`
$N \rtimes_\varphi H$ mit Operation $(n_1, h_1)(n_2, h_2) = (n_1 \cdot \varphi(h_1)(n_2),\; h_1 h_2)$.

---

## Beispiele

```python
from group_theory import cyclic_group, dihedral_group, symmetric_group

# Zyklische Gruppe ℤ_6
Z6 = cyclic_group(6)
print(Z6.is_abelian())   # True
print(Z6.is_cyclic())    # True
print(Z6.order())        # 6

# Diedergruppe D_4 (Quadrat-Symmetrien)
D4 = dihedral_group(4)
print(D4.order())        # 8
print(D4.is_abelian())   # False

# Symmetrische Gruppe S_3
S3 = symmetric_group(3)
print(S3.order())        # 6
subs = S3.subgroups()
print(len(subs))         # 6 (Untergruppen von S_3)

# Lagrange-Satz
from group_theory import lagrange_theorem_check
H = next(h for h in subs if h.order() == 3)
result = lagrange_theorem_check(S3, H)
print(result['index'])   # 2

# Klassifikation
from group_theory import classify_abelian_group
r = classify_abelian_group(8)
print(r['description'])  # ['Z_8', 'Z_4 × Z_2', 'Z_2 × Z_2 × Z_2']
```

---

## Bekannte Gruppen und ihre Eigenschaften

| Gruppe | Ordnung | Abelsch | Einfach | Zyklisch |
|--------|---------|---------|---------|----------|
| $\mathbb{Z}_n$ (prim) | $n$ | ✓ | ✓ | ✓ |
| $S_n$, $n \geq 3$ | $n!$ | ✗ | ✗ | ✗ |
| $A_n$, $n \geq 5$ | $n!/2$ | ✗ | ✓ | ✗ |
| $D_n$, $n \geq 3$ | $2n$ | ✗ | ✗ | ✗ |
| $Q_8$ | $8$ | ✗ | ✗ | ✗ |
| $A_3 \cong \mathbb{Z}_3$ | $3$ | ✓ | ✓ | ✓ |
