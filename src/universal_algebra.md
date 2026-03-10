# universal_algebra.py — Universelle Algebra

## Übersicht

Die **universelle Algebra** (auch allgemeine Algebra) untersucht algebraische
Strukturen auf einer abstrakten Ebene — unabhängig von spezifischen Realisierungen.
Statt Gruppen, Ringe oder Verbände einzeln zu betrachten, entwickelt sie eine
einheitliche Theorie für alle solchen Strukturen.

---

## Grundbegriffe

### Signatur (Typ)

Eine **Signatur** $\Sigma = (F, \text{ar})$ legt fest, welche Operationssymbole
mit welchen Aritäten vorhanden sind:

| Symbol | Arität | Bedeutung |
|--------|--------|-----------|
| $\cdot$ | 2 | Binäre Operation |
| $e$ | 0 | Konstante (ausgezeichnetes Element) |
| $^{-1}$ | 1 | Unäre Operation |

**Gruppen-Signatur:** $\{{\cdot}/2,\; e/0,\; {}^{-1}/1\}$

```python
sig = Signature({'mul': 2, 'e': 0, 'inv': 1})
sig.has_constants()  # True
sig.has_unary()      # True
```

---

### Algebra

Eine **$\Sigma$-Algebra** $(A, (f_i^A)_{i \in I})$ ist eine Menge $A$ mit
Operationsinterpretationen $f_i^A: A^{\text{ar}(i)} \to A$.

```python
# ℤ/3ℤ als Algebra
sig = Signature({'add': 2, 'zero': 0})
alg = Algebra([0, 1, 2], sig, {
    'add': lambda a, b: (a + b) % 3,
    'zero': lambda: 0,
})
alg.apply('add', 1, 2)  # → 0
```

**Unteralgebra:** $B \subseteq A$ heißt Unteralgebra, wenn $B$ unter allen
Operationen abgeschlossen ist: $\forall f \in F, b_1,\ldots,b_n \in B: f(b_1,\ldots,b_n) \in B$

**Erzeugte Unteralgebra:** $[G]$ = kleinste Unteralgebra die $G$ enthält
(Fixpunktiteration).

---

## Varietäten

Eine **Varietät** (equationale Klasse) ist eine Klasse von Algebren, die durch
Gleichungsaxiome definiert ist.

### Klassische Varietäten

| Varietät | Signatur | Axiome |
|----------|----------|--------|
| Gruppen | $\{\cdot, e, {}^{-1}\}$ | Assoziativität, Einheit, Inverse |
| Abelsche Gruppen | wie Gruppen | + Kommutativität |
| Ringe | $\{+, \cdot, 0, 1, -\}$ | Gruppe für +, Monoid für ·, Distributivität |
| Verbände | $\{\wedge, \vee\}$ | Kommutativität, Assoziativität, Absorption |

---

## Birkhoff-Charakterisierungssatz

> **Satz (Birkhoff, 1935):** Eine Klasse $K$ von $\Sigma$-Algebren ist genau dann
> eine Varietät (durch Gleichungen axiomatisierbar), wenn $K = HSP(K)$.

Die drei Operatoren:

$$H(K) = \{\varphi(A) \mid A \in K,\; \varphi \text{ surjektiver Homomorphismus}\}$$

$$S(K) = \{B \mid B \leq A \text{ für ein } A \in K\}$$

$$P(K) = \left\{\prod_{i \in I} A_i \mid \text{alle } A_i \in K\right\}$$

**Bedeutung:**
- Varietäten sind unter $H$, $S$, $P$ abgeschlossen
- Umgekehrt: HSP-abgeschlossene Klassen sind immer Varietäten
- Verknüpft strukturelle Eigenschaften mit equationaler Axiomatisierung

```python
result = birkhoff_theorem_demo()
# {'theorem_name': 'Birkhoff-Charakterisierungssatz (1935)',
#  'operators': {'H': '...', 'S': '...', 'P': '...'},
#  'examples': [...],  # Gruppen, Ringe, Verbände
#  ...}
```

---

## Freie Algebren

### Term-Algebra

Die **Term-Algebra** $T_\Sigma(X)$ über Generatoren $X$ ist die freie Algebra in
der Varietät aller $\Sigma$-Algebren. Ihre Elemente sind alle sinnvollen Terme:

- **Tiefe 0:** Generatoren $x \in X$ und Konstanten
- **Tiefe 1:** $f(t_1, \ldots, t_n)$ für $f \in F$ und Terme $t_i$ der Tiefe 0
- **Tiefe $d$:** rekursiv aufgebaut

**Universelle Eigenschaft:** Für jede $\Sigma$-Algebra $A$ und Funktion $f: X \to A$
gibt es genau einen Homomorphismus $\varphi: T_\Sigma(X) \to A$ mit $\varphi|_X = f$.

```python
sig = Signature({'mul': 2, 'e': 0})
ta = term_algebra(sig, ['x', 'y'])
# depth_0_terms: ['x', 'y', 'e']
# depth_1_terms: ['mul(x,x)', 'mul(x,y)', ..., 'mul(y,e)', ...]
```

---

## Kongruenzen und Quotientenalgebren

Eine **Kongruenz** $\theta$ auf $A$ ist eine Äquivalenzrelation, die mit allen
Operationen verträglich ist:

$$\forall f \in F: (a_1, b_1), \ldots, (a_n, b_n) \in \theta \Rightarrow (f(a_1,\ldots,a_n), f(b_1,\ldots,b_n)) \in \theta$$

Die **Quotientenalgebra** $A/\theta$ hat Äquivalenzklassen als Elemente:

$$f([a_1], \ldots, [a_n]) = [f(a_1, \ldots, a_n)]$$

(wohldefiniert durch Kongruenzeigenschaft)

---

## Subdirekte Irreduzibilität

Eine Algebra $A$ heißt **subdirekt irreduzibel**, wenn der Durchschnitt aller
nicht-trivialen Kongruenzen (≠ Gleichheitsrelation $\Delta$) selbst nicht-trivial ist.

> **Birkhoffs Darstellungssatz:** Jede Algebra ist subdirekt isomorph zu einem
> Produkt subdirekt irreduzibler Algebren.

Diese subdirekt irreduziblen Algebren sind die "Bausteine" jeder Varietät.

---

## Verwendungsbeispiele

```python
from universal_algebra import (
    Signature, Algebra, Variety,
    groups_variety, birkhoff_theorem_demo,
    congruence_relation, quotient_algebra
)

# Gruppen-Varietät
v = groups_variety()
print(v.axioms)  # Alle Gruppenaxiome

# Birkhoff-Satz Demonstration
demo = birkhoff_theorem_demo()
print(demo['statement'])

# Kongruenz prüfen
sig = Signature({'add': 2, 'zero': 0})
alg = Algebra([0,1,2], sig, {'add': lambda a,b: (a+b)%3, 'zero': lambda: 0})
is_cong = congruence_relation(alg, [(0, 1)])  # True/False
```
