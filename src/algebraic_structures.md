# algebraic_structures.py — Algebraische Strukturhierarchie

## Übersicht

Dieses Modul implementiert die vollständige Hierarchie algebraischer Strukturen
von der allgemeinsten (Magma) bis hin zur speziellsten Struktur (Gruppe), sowie
freie algebraische Strukturen und konkrete Beispiele.

## Algebraische Strukturhierarchie

```
Magma
  └── Halbgruppe (+ Assoziativität)
        └── Monoid (+ Neutrales Element)
              └── Gruppe (+ Inverse)
                    └── Abelsche Gruppe (+ Kommutativität)
```

## Klassen

### `Magma`

Die allgemeinste algebraische Struktur: eine Menge $M$ mit einer inneren binären Operation $*$.

**Axiom:** $\forall a, b \in M: a * b \in M$ (Abgeschlossenheit)

**Methoden:**
- `operation_table()` — Cayley-Tabelle als 2D-Liste
- `is_closed()` — Prüft Abgeschlossenheit
- `is_commutative()` — Prüft $a * b = b * a$

---

### `Semigroup(Magma)`

Zusätzliches Axiom: **Assoziativität**

$$\forall a, b, c \in S: (a * b) * c = a * (b * c)$$

**Methoden:**
- `is_associative()` — Prüft Assoziativgesetz (O(n³))
- `is_semigroup()` — Vollständige Halbgruppen-Prüfung

**Beispiele:** $(\mathbb{N}, +)$, $(\mathbb{N}, \cdot)$, $(\Sigma^+, \text{Konkat.})$

---

### `Monoid(Semigroup)`

Zusätzliches Axiom: **Neutrales Element**

$$\exists e \in M: \forall a \in M: e * a = a * e = a$$

**Methoden:**
- `identity_element()` — Sucht und gibt $e$ zurück
- `has_identity()` — Prüft Existenz von $e$
- `is_monoid()` — Vollständige Monoid-Prüfung
- `powers(element, n)` — Berechnet $[e, a, a^2, \ldots, a^n]$

**Beispiele:** $(\mathbb{N}_0, +, 0)$, $(\mathbb{N}, \cdot, 1)$, $(\Sigma^*, \text{Konkat.}, \varepsilon)$

---

### `GroupFromMagma(Monoid)`

Zusätzliches Axiom: **Inverse**

$$\forall a \in G \; \exists a^{-1} \in G: a * a^{-1} = a^{-1} * a = e$$

**Methoden:**
- `inverse(element)` — Findet $a^{-1}$
- `is_group()` — Vollständige Gruppen-Prüfung

---

## Freie Strukturen

### `free_monoid(alphabet)`

Der **freie Monoid** $(\Sigma^*, \cdot, \varepsilon)$ über einem Alphabet $\Sigma$:
- Elemente: alle endlichen Wörter (Tupel) über $\Sigma$
- Operation: Konkatenation
- Neutrales Element: leeres Wort $\varepsilon = ()$

**Universelle Eigenschaft:** Jede Funktion $f: \Sigma \to M$ in einen Monoid $M$
faktorisiert eindeutig über $\Sigma^*$.

### `free_semigroup(generators, max_length)`

Die **freie Halbgruppe** $(\Sigma^+, \cdot)$ — alle nicht-leeren Wörter bis `max_length`.

---

## Beispiel-Monoide

### `string_monoid(alphabet, max_length)`

Monoid der Strings über endlichem Alphabet mit abgeschnittener Konkatenation.

```python
m = string_monoid("ab", max_length=2)
# Elemente: '', 'a', 'b', 'aa', 'ab', 'ba', 'bb'
m.identity_element()  # → ''
```

### `matrix_monoid(n, field='Z2')`

Monoid der $n \times n$-Matrizen über $\mathbb{Z}/2\mathbb{Z}$ unter Multiplikation.

- $2^{n^2}$ Elemente (für $n=2$: 16 Matrizen)
- Neutrales Element: Einheitsmatrix $I_n$

### `transformation_monoid(n)`

Transformationsmonoid $T_n$: alle Funktionen $\{1,\ldots,n\} \to \{1,\ldots,n\}$.

- $n^n$ Elemente
- Operation: Komposition $f \circ g$
- Neutrales Element: Identitätsabbildung

---

## Hierarchie-Übersicht

```python
h = algebraic_structure_hierarchy()
# Gibt Dict zurück mit allen Strukturen, Axiomen, Beispielen
```

Hierarchie:
$$\text{Magma} \supset \text{Halbgruppe} \supset \text{Monoid} \supset \text{Gruppe} \supset \text{Abel. Gruppe}$$
$$\supset \text{Ring} \supset \text{Integritätsbereich} \supset \text{Körper}$$

---

## Verwendungsbeispiele

```python
from algebraic_structures import Magma, Semigroup, Monoid, GroupFromMagma

# ℤ/3ℤ als Gruppe
elements = [0, 1, 2]
op = lambda a, b: (a + b) % 3
g = GroupFromMagma(elements, op)
print(g.is_group())          # True
print(g.inverse(1))          # 2 (da 1+2=0)
print(g.powers(1, 5))        # [0, 1, 2, 0, 1, 2]

# Cayley-Tabelle ausgeben
table = g.operation_table()
```

## Mathematische Notizen

- Die Klassen sind durch Vererbung verbunden: jede Unterklasse **prüft** die Axiome der Oberklasse
- `is_associative()` hat Laufzeit $O(n^3)$ — bei großen Strukturen langsam
- Das neutrale Element ist eindeutig in Monoiden (Beweis: $e = e * e' = e'$)
- Inverse sind in Gruppen eindeutig (Beweis: $b = b*e = b*(a*c) = (b*a)*c = e*c = c$)
