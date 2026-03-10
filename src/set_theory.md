# set_theory.py – Mathematische Mengenlehre

**Modul:** `src/set_theory.py`
**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Tests:** `tests/test_set_theory.py` (138 Tests, alle grün)

---

## Überblick

Das Modul implementiert die vollständige mathematische Mengenlehre objektorientiert
und domänengetrieben. Es deckt alle fundamentalen Konzepte der axiomatischen
Mengenlehre ab – von einfachen Mengenoperationen bis hin zu transfiniten
Kardinal- und Ordinalzahlen sowie dem ZFC-Axiomensystem.

---

## 1. MathSet – Mathematische Menge

Eine Menge ist eine Zusammenfassung von Objekten (Elementen) zu einem Ganzen.
Intern als `frozenset` gespeichert, damit `MathSet`-Instanzen selbst als
Elemente anderer Mengen (z. B. Potenzmenge) verwendet werden können.

### Mengenoperationen

| Operation | Notation | Methode | Operator |
|-----------|----------|---------|----------|
| Vereinigung | $A \cup B$ | `union(B)` | `A \| B` |
| Durchschnitt | $A \cap B$ | `intersection(B)` | `A & B` |
| Differenz | $A \setminus B$ | `difference(B)` | `A - B` |
| Symmetr. Differenz | $A \triangle B$ | `symmetric_difference(B)` | `A ^ B` |
| Komplement | $A^c = U \setminus A$ | `complement(U)` | – |
| Potenzmenge | $\mathcal{P}(A)$ | `power_set()` | – |
| Kartesisches Produkt | $A \times B$ | `cartesian_product(B)` | – |

### Mengeneigenschaften

| Eigenschaft | Notation | Methode | Operator |
|-------------|----------|---------|----------|
| Teilmenge | $A \subseteq B$ | `is_subset(B)` | `A <= B` |
| Echte Teilmenge | $A \subsetneq B$ | `is_proper_subset(B)` | `A < B` |
| Disjunktheit | $A \cap B = \emptyset$ | `is_disjoint(B)` | – |
| Mächtigkeit | $|A|$ | `cardinality()` | `len(A)` |

### Gesetze der Mengenalgebra

**Kommutativgesetze:**
$$A \cup B = B \cup A \qquad A \cap B = B \cap A$$

**Assoziativgesetze:**
$$(A \cup B) \cup C = A \cup (B \cup C) \qquad (A \cap B) \cap C = A \cap (B \cap C)$$

**Distributivgesetze:**
$$A \cup (B \cap C) = (A \cup B) \cap (A \cup C)$$
$$A \cap (B \cup C) = (A \cap B) \cup (A \cap C)$$

**De Morganschen Gesetze:**
$$(A \cup B)^c = A^c \cap B^c \qquad (A \cap B)^c = A^c \cup B^c$$

### Potenzmenge

Für eine $n$-elementige Menge $A$ gilt:
$$|\mathcal{P}(A)| = 2^n$$

Die Potenzmenge enthält immer $\emptyset$ und $A$ selbst.

---

## 2. Relation – Binäre Relation

Eine binäre Relation $R \subseteq A \times B$ ist eine Menge geordneter Paare.

### Relationseigenschaften (für $R \subseteq A \times A$)

| Eigenschaft | Definition | Methode |
|-------------|-----------|---------|
| Reflexiv | $\forall a \in A: (a,a) \in R$ | `is_reflexive()` |
| Symmetrisch | $(a,b) \in R \Rightarrow (b,a) \in R$ | `is_symmetric()` |
| Antisymmetrisch | $(a,b),(b,a) \in R \Rightarrow a=b$ | `is_antisymmetric()` |
| Transitiv | $(a,b),(b,c) \in R \Rightarrow (a,c) \in R$ | `is_transitive()` |
| Äquivalenzrelation | reflexiv + symmetrisch + transitiv | `is_equivalence()` |
| Partielle Ordnung | reflexiv + antisymmetrisch + transitiv | `is_partial_order()` |
| Totalordnung | partielle Ordnung + $\forall a,b: (a,b) \in R$ oder $(b,a) \in R$ | `is_total_order()` |

### Äquivalenzklassen

Für eine Äquivalenzrelation $R$ auf $A$ ist die **Äquivalenzklasse** von $a$:
$$[a]_R = \{b \in A : (a,b) \in R\}$$

Die **Quotientenmenge** ist:
$$A/R = \{[a]_R : a \in A\}$$

Die Äquivalenzklassen bilden eine **Partition** von $A$.

### Abschlussoperationen

| Operation | Definition | Methode |
|-----------|-----------|---------|
| Transitiver Abschluss | Kleinste transitive Oberrelation | `transitive_closure()` |
| Reflexiver Abschluss | $R \cup \{(a,a) : a \in A\}$ | `reflexive_closure()` |
| Symmetrischer Abschluss | $R \cup R^{-1}$ | `symmetric_closure()` |
| Inverse Relation | $R^{-1} = \{(b,a) : (a,b) \in R\}$ | `inverse()` |
| Komposition | $(a,c) \in S \circ R \Leftrightarrow \exists b: (a,b) \in R \wedge (b,c) \in S$ | `composition(S)` |

Der transitive Abschluss wird via **Floyd-Warshall-Algorithmus** berechnet.

---

## 3. MathFunction – Mathematische Funktion

Eine Funktion $f: A \to B$ ist eine linkstotale, rechtseindeutige Relation:
- **Linkstotal:** $\forall a \in A \; \exists b \in B: f(a) = b$
- **Rechtseindeutig:** $f(a) = b \wedge f(a) = b' \Rightarrow b = b'$

### Funktionseigenschaften

| Eigenschaft | Definition | Methode |
|-------------|-----------|---------|
| Injektiv | $f(a) = f(b) \Rightarrow a = b$ | `is_injective()` |
| Surjektiv | $\forall b \in B \; \exists a \in A: f(a) = b$ | `is_surjective()` |
| Bijektiv | injektiv $\wedge$ surjektiv | `is_bijective()` |
| Identität | $f(a) = a \; \forall a \in A$ | `is_identity()` |

### Bild und Urbild

$$f(S) = \{f(a) : a \in S\} \quad \text{(Bild von } S \subseteq A\text{)}$$
$$f^{-1}(T) = \{a \in A : f(a) \in T\} \quad \text{(Urbild von } T \subseteq B\text{)}$$

### Komposition

$$(g \circ f)(a) = g(f(a)) \quad \text{für } f: A \to B, \; g: B \to C$$

---

## 4. Kardinalitätstheorie

### Gleichmächtigkeit

Zwei Mengen $A$ und $B$ sind **gleichmächtig** ($|A| = |B|$), wenn eine
Bijektion $h: A \to B$ existiert.

### Cantor-Schröder-Bernstein-Satz

**Satz:** Wenn $f: A \to B$ injektiv und $g: B \to A$ injektiv, dann existiert
eine Bijektion $h: A \to B$.

**Konstruktion:**

Definiere die Schröder-Bernstein-Kette:
$$A_0 = A \setminus g(B), \quad A_{n+1} = g(f(A_n))$$
$$C = \bigcup_{n=0}^\infty A_n$$

Dann:
$$h(a) = \begin{cases} f(a) & \text{falls } a \in C \\ g^{-1}(a) & \text{falls } a \notin C \end{cases}$$

### Abzählbarkeit

| Menge | Mächtig­keit | Abzählbar? | Beweis-Idee |
|-------|------------|-----------|------------|
| $\mathbb{N}$ | $\aleph_0$ | ja | Identität |
| $\mathbb{Z}$ | $\aleph_0$ | ja | Alternierend: $0,1,-1,2,-2,\ldots$ |
| $\mathbb{Q}$ | $\aleph_0$ | ja | Cantors Diagonalzählung |
| $\mathbb{R}$ | $2^{\aleph_0}$ | **nein** | Cantors Diagonalargument |
| $\mathcal{P}(\mathbb{N})$ | $2^{\aleph_0}$ | **nein** | Cantors Satz: $|A| < |\mathcal{P}(A)|$ |

### Cantors Diagonalargument

**Theorem:** $|\mathbb{N}| < |\mathcal{P}(\mathbb{N})|$

**Beweis:** Angenommen, es gibt eine Aufzählung $S_0, S_1, S_2, \ldots$ aller
Teilmengen von $\mathbb{N}$. Definiere:
$$D = \{n \in \mathbb{N} : n \notin S_n\}$$

Dann gilt für alle $n \in \mathbb{N}$:
$$n \in D \Leftrightarrow n \notin S_n \Rightarrow D \neq S_n$$

Widerspruch: $D$ ist eine Teilmenge von $\mathbb{N}$, taucht aber in keiner
Aufzählung auf. $\square$

### Beth-Zahlen und Aleph-Zahlen

$$\beth_0 = \aleph_0 = |\mathbb{N}|$$
$$\beth_{n+1} = 2^{\beth_n} = |\mathcal{P}^{(n+1)}(\mathbb{N})|$$

$$\beth_0 = \aleph_0, \quad \beth_1 = 2^{\aleph_0} = |\mathbb{R}|, \quad \beth_2 = 2^{2^{\aleph_0}} = |\mathcal{P}(\mathbb{R})|$$

**Kontinuumshypothese (CH):** $2^{\aleph_0} = \aleph_1$

CH ist **unabhängig von ZFC** (Gödel 1938, Cohen 1963).

---

## 5. Ordinalzahlen

### Definition

Ordinalzahlen beschreiben **Wohlordnungstypen**. Jede Wohlordnung ist isomorph
zu genau einer Ordinalzahl.

| Ordinalzahl | Bedeutung |
|------------|-----------|
| $0, 1, 2, \ldots$ | Endliche Ordinalzahlen |
| $\omega$ | Kleinste unendliche Ordinalzahl (Typ von $\mathbb{N}$) |
| $\omega + 1$ | Nachfolger von $\omega$ |
| $\omega \cdot 2 = \omega + \omega$ | Zwei Kopien von $\omega$ |
| $\omega^2$ | $\omega$ Kopien von $\omega$ |

### Nicht-Kommutativität der Ordinalzahlarithmetik

$$\omega + 1 \neq 1 + \omega$$

- $\omega + 1 = \{0, 1, 2, \ldots, \omega\}$ – hat ein größtes Element
- $1 + \omega = \omega$ – das $1$ links wird "geschluckt"

$$\omega \cdot 2 \neq 2 \cdot \omega$$

- $\omega \cdot 2 = \omega + \omega$ – zwei Kopien
- $2 \cdot \omega = \omega$ – endlich viele Kopien bleiben $\omega$

### Limeszahlen und Nachfolger-Ordinalzahlen

- **Nachfolger-Ordinalzahl:** $\alpha + 1$ für eine Ordinalzahl $\alpha$
- **Limeszahl:** Weder $0$ noch Nachfolger (z. B. $\omega, \omega \cdot 2, \omega^2$)

### Wohlordnungssatz

Jede Totalordnung auf einer endlichen Menge ist eine Wohlordnung
(jede nicht-leere Teilmenge hat ein kleinstes Element).

---

## 6. ZFC-Axiome

Das **Zermelo-Fraenkel-System mit Auswahlaxiom** (ZFC) ist das
Standard-Axiomensystem der modernen Mathematik.

| # | Name | Formale Aussage |
|---|------|----------------|
| 1 | Extensionalität | $\forall A \forall B [\forall x(x \in A \leftrightarrow x \in B) \to A = B]$ |
| 2 | Leere Menge | $\exists A \forall x \neg(x \in A)$ |
| 3 | Paarmenge | $\forall a \forall b \exists P \forall x(x \in P \leftrightarrow x=a \vee x=b)$ |
| 4 | Vereinigung | $\forall \mathcal{F} \exists U \forall x[x \in U \leftrightarrow \exists A(A \in \mathcal{F} \wedge x \in A)]$ |
| 5 | Potenzmenge | $\forall A \exists P \forall B[B \in P \leftrightarrow B \subseteq A]$ |
| 6 | Unendlichkeit | $\exists I[\emptyset \in I \wedge \forall x \in I(x \cup \{x\} \in I)]$ |
| 7 | Aussonderung | $\forall A \forall p \exists S \forall x[x \in S \leftrightarrow x \in A \wedge \varphi(x,p)]$ |
| 8 | Ersetzung | $\forall A[\forall x \in A \exists! y\, \varphi(x,y) \to \exists B \forall x \in A \exists y \in B\, \varphi(x,y)]$ |
| 9 | Fundierung | $\forall A[A \neq \emptyset \to \exists x \in A(x \cap A = \emptyset)]$ |
| 10 | Auswahlaxiom | $\forall \mathcal{F}[\emptyset \notin \mathcal{F} \to \exists f: \mathcal{F} \to \bigcup \mathcal{F} \; \forall A \in \mathcal{F}(f(A) \in A)]$ |

### Äquivalente zum Auswahlaxiom

Das Auswahlaxiom (AC) ist in ZF äquivalent zu:
- **Lemma von Zorn:** Jede Ketten-vollständige geordnete Menge hat ein maximales Element
- **Wohlordnungssatz (Zermelo):** Jede Menge kann wohlgeordnet werden
- **Tychonoff-Satz:** Produkt kompakter Räume ist kompakt
- **Hamel-Basis:** Jeder Vektorraum hat eine Basis

---

## 7. SetFamily – Mengenfamilien

Eine Mengenfamilie $\{A_i : i \in I\}$ ist eine indizierte Sammlung von Mengen.

### σ-Algebra

Eine Familie $\Sigma$ über dem Universum $U$ ist eine **σ-Algebra**, wenn:
1. $U \in \Sigma$
2. $A \in \Sigma \Rightarrow A^c \in \Sigma$ (Komplement-Abschluss)
3. $A_1, A_2, \ldots \in \Sigma \Rightarrow \bigcup_n A_n \in \Sigma$ (Abzählbare-Vereinigungs-Abschluss)

Kleinstes Beispiel: **Triviale σ-Algebra** $\{\emptyset, U\}$
Größtes Beispiel: **Diskrete σ-Algebra** $\mathcal{P}(U)$

### Partition

Eine Familie $\{A_i\}$ ist eine **Partition** von $U$, wenn:
- $A_i \neq \emptyset$ für alle $i$
- $\bigcup_i A_i = U$ (Überdeckung)
- $A_i \cap A_j = \emptyset$ für $i \neq j$ (paarweise Disjunktheit)

---

## 8. Filter

Ein **Filter** $\mathcal{F}$ auf $U$ ist eine nicht-leere Familie von
Teilmengen von $U$ mit:
1. $\emptyset \notin \mathcal{F}$ (Nicht-Trivialität)
2. $A \in \mathcal{F}, A \subseteq B \subseteq U \Rightarrow B \in \mathcal{F}$ (Aufwärtsabschluss)
3. $A, B \in \mathcal{F} \Rightarrow A \cap B \in \mathcal{F}$ (Schnitt-Abschluss)

Ein **Ultrafilter** ist ein maximaler Filter:
$$\forall A \subseteq U: A \in \mathcal{F} \text{ oder } A^c \in \mathcal{F}$$

### Hauptfilter

Der von $B \subseteq U$ erzeugte **Hauptfilter**:
$$\mathcal{F}_B = \{A \subseteq U : B \subseteq A\}$$

---

## Verwendungsbeispiele

```python
from set_theory import MathSet, Relation, MathFunction, zfc_axioms, Ordinal

# Mengen erstellen und verknüpfen
A = MathSet([1, 2, 3])
B = MathSet([2, 3, 4])
print(A | B)   # {1, 2, 3, 4}
print(A & B)   # {2, 3}
print(A ^ B)   # {1, 4}

# Potenzmenge
print(len(A.power_set()))  # 8 = 2^3

# Äquivalenzrelation
R = Relation({(1,1),(2,2),(3,3),(1,3),(3,1)}, A, A)
print(R.is_equivalence())     # True
print(R.equivalence_classes()) # [{1,3}, {2}]

# Bijektive Funktion
C = MathSet([4, 5, 6])
f = MathFunction({1:4, 2:5, 3:6}, A, C)
print(f.is_bijective())       # True
print(f.inverse_function()(4))  # 1

# Ordinalzahlen (nicht-kommutativ!)
omega = Ordinal('omega')
print(omega + Ordinal(1))     # ω+1
print(Ordinal(1) + omega)     # ω   (≠ ω+1!)

# ZFC-Axiome
axioms = zfc_axioms()
print(axioms['choice']['name'])  # Auswahlaxiom
```

---

## Komplexität

| Operation | Zeitkomplexität |
|-----------|----------------|
| Vereinigung/Durchschnitt | $O(\min(|A|, |B|))$ |
| Potenzmenge | $O(2^n)$ |
| Kartesisches Produkt | $O(|A| \cdot |B|)$ |
| Transitiver Abschluss (Floyd-Warshall) | $O(n^3)$ |
| Äquivalenzklassen | $O(|R|^2)$ |
| σ-Algebra-Abschluss | $O(2^{|U|})$ |
| Filter-Prüfung | $O(2^{|U|} \cdot |F|^2)$ |
