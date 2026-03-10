# category_theory.py – Dokumentation

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Sprache:** Python 3.13

---

## Überblick

Das Modul `category_theory.py` implementiert die grundlegenden Konzepte der
**Kategorientheorie** – eines fundamentalen Zweigs der modernen Mathematik,
der universelle Strukturen und deren Beziehungen untersucht.

Kategorientheorie ist das "Gerüst", auf dem moderne Algebra, Topologie und
Informatik (Typsysteme, funktionale Programmierung) aufgebaut sind.

---

## Mathematische Grundlagen

### Definition einer Kategorie

Eine **Kategorie** $\mathcal{C}$ besteht aus:

1. **Objekten** $\text{Ob}(\mathcal{C})$ – beliebige abstrakte Entitäten
2. **Morphismen** $\text{Mor}_\mathcal{C}(A, B)$ für je zwei Objekte $A, B$
3. **Komposition**: Für $f: A \to B$ und $g: B \to C$ existiert $g \circ f: A \to C$
4. **Identitäten**: Für jedes $A$ ein $\text{id}_A: A \to A$

**Axiome:**

$$h \circ (g \circ f) = (h \circ g) \circ f \quad \text{(Assoziativität)}$$

$$\text{id}_B \circ f = f = f \circ \text{id}_A \quad \text{(Einheitsgesetze)}$$

### Beispielkategorien

| Kategorie | Objekte | Morphismen |
|-----------|---------|------------|
| **Set** | Mengen | Funktionen |
| **Grp** | Gruppen | Gruppenhomomorphismen |
| **Top** | Topologische Räume | Stetige Abbildungen |
| **Vect_k** | Vektorräume über $k$ | Lineare Abbildungen |

---

## Klassen

### `Object`

Repräsentiert ein Objekt $A \in \text{Ob}(\mathcal{C})$.

```python
obj = Object("A", data=some_set)
```

Objekte werden durch ihren Namen identifiziert (`__eq__` und `__hash__` basieren auf `name`).

---

### `Morphism`

Repräsentiert einen Morphismus $f: A \to B$.

```python
f = Morphism(source=A, target=B, name="f", func=lambda x: x+1)
```

**Komposition:**

$$g \circ f: A \to C \quad \text{für } f: A \to B,\ g: B \to C$$

```python
gof = g.compose(f)   # Gibt neuen Morphismus zurück
```

---

### `Category`

Verwaltet Objekte und Morphismen einer Kategorie.

```python
cat = Category("C")
cat.add_object(A)      # Erstellt automatisch id_A
cat.add_morphism(f)
gof = cat.compose(f, g)   # g ∘ f
id_A = cat.identity(A)
```

---

### `Functor`

Ein **kovarianter Funktor** $F: \mathcal{C} \to \mathcal{D}$ besteht aus:

$$F: \text{Ob}(\mathcal{C}) \to \text{Ob}(\mathcal{D})$$

$$F: \text{Mor}_\mathcal{C}(A,B) \to \text{Mor}_\mathcal{D}(F(A), F(B))$$

**Eigenschaften:**

$$F(\text{id}_A) = \text{id}_{F(A)} \quad \text{(Identitätserhaltung)}$$

$$F(g \circ f) = F(g) \circ F(f) \quad \text{(Kompositionserhaltung)}$$

```python
F = Functor("F", source_cat, target_cat, object_map, morphism_map)
F.preserves_identity()     # → bool
F.preserves_composition()  # → bool
```

---

### `NaturalTransformation`

Eine **natürliche Transformation** $\alpha: F \Rightarrow G$ zwischen Funktoren
$F, G: \mathcal{C} \to \mathcal{D}$ ist eine Familie von Morphismen:

$$\alpha_A: F(A) \to G(A) \quad \text{für jedes } A \in \text{Ob}(\mathcal{C})$$

**Natürlichkeitsbedingung** (die Naturquadrate kommutieren):

$$\alpha_B \circ F(f) = G(f) \circ \alpha_A \quad \text{für alle } f: A \to B$$

Diagrammatisch:

$$\begin{array}{ccc}
F(A) & \xrightarrow{F(f)} & F(B) \\
\downarrow_{\alpha_A} & & \downarrow_{\alpha_B} \\
G(A) & \xrightarrow{G(f)} & G(B)
\end{array}$$

```python
alpha = NaturalTransformation("α", F, G, components={A: alpha_A, B: alpha_B})
alpha.naturality_check()  # → bool
```

---

## Freie Funktionen

### `set_category(sets)`

Erstellt die **Set-Kategorie** für gegebene Python-Mengen:

$$\text{Set}: \quad \text{Ob} = \text{gegebene Mengen}, \quad \text{Mor}(A, B) = B^A$$

### `opposite_category(cat)`

Erstellt die **duale Kategorie** $\mathcal{C}^{op}$:

$$\text{Mor}_{\mathcal{C}^{op}}(A, B) = \text{Mor}_\mathcal{C}(B, A)$$

### `product_category(cat1, cat2)`

**Produktkategorie** $\mathcal{C}_1 \times \mathcal{C}_2$:

- Objekte: Paare $(A, B)$
- Morphismen: Paare $(f, g): (A, B) \to (A', B')$

### `is_isomorphism(f, cat)`

$f: A \to B$ ist ein **Isomorphismus**, falls $\exists g: B \to A$ mit:

$$g \circ f = \text{id}_A \quad \text{und} \quad f \circ g = \text{id}_B$$

### `is_epimorphism(f, cat)`

$f: A \to B$ ist ein **Epimorphismus** (rechts-kürzbar):

$$g \circ f = h \circ f \implies g = h$$

In **Set**: surjektive Funktionen.

### `is_monomorphism(f, cat)`

$f: A \to B$ ist ein **Monomorphismus** (links-kürzbar):

$$f \circ g = f \circ h \implies g = h$$

In **Set**: injektive Funktionen.

### `hom_set(A, B, cat)`

Gibt das **Hom-Set** zurück:

$$\text{Hom}_\mathcal{C}(A, B) = \{f \in \text{Mor}(\mathcal{C}) \mid f: A \to B\}$$

### `yoneda_functor(cat, obj)`

Das **Yoneda-Lemma** ist eines der fundamentalsten Ergebnisse der Kategorientheorie:

$$\text{Nat}(h_A, F) \cong F(A)$$

für jeden Funktor $F: \mathcal{C} \to \text{Set}$, wobei:

$$h_A = \text{Hom}_\mathcal{C}(A, -): \mathcal{C} \to \text{Set}$$

Die **Yoneda-Einbettung** $A \mapsto h_A$ ist vollständig und treu:

$$\mathcal{C} \hookrightarrow [\mathcal{C}^{op}, \text{Set}]$$

---

## Verwendungsbeispiel

```python
from category_theory import Object, Morphism, Category, Functor

# Einfache Kategorie erstellen
cat = Category("MeineKategorie")
A = Object("A")
B = Object("B")
C = Object("C")

cat.add_object(A)
cat.add_object(B)
cat.add_object(C)

# Morphismen
f = Morphism(A, B, "f", func=lambda x: x + 1)
g = Morphism(B, C, "g", func=lambda x: x * 2)
cat.add_morphism(f)
cat.add_morphism(g)

# Komposition g ∘ f: A → C
gof = cat.compose(f, g)
print(gof.func(3))  # (3+1)*2 = 8

# Identitätsmorphismus
id_A = cat.identity(A)
print(id_A.name)  # "id_A"
```

---

## Literatur

- Saunders Mac Lane: *Categories for the Working Mathematician* (1998)
- Awodey, Steve: *Category Theory* (2010)
- Riehl, Emily: *Category Theory in Context* (2016)
- nLab: https://ncatlab.org/nlab/show/category+theory
