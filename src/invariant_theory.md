# invariant_theory.py — Invariantentheorie

## Übersicht

Dieses Modul implementiert die klassische **Invariantentheorie**: das Studium von
Polynomen, die unter Gruppenwirkungen unverändert bleiben.

## Grundbegriffe

Sei $G$ eine endliche Gruppe, die auf dem Polynomring $k[x_1, \ldots, x_n]$ wirkt.
Ein Polynom $f$ heißt **$G$-invariant**, wenn für alle $g \in G$ gilt:

$$g \cdot f = f$$

Der **Invariantenring** $k[x_1, \ldots, x_n]^G$ ist die Menge aller $G$-invarianten Polynome.

---

## Reynolds-Operator

Der **Reynolds-Operator** projiziert ein beliebiges Polynom auf den invarianten Unterraum:

$$R(f) = \frac{1}{|G|} \sum_{g \in G} g \cdot f$$

**Eigenschaften:**
- $R \circ R = R$ (Projektor)
- $\text{Bild}(R) = k[x_1,\ldots,x_n]^G$
- Wohldefiniert wenn $\text{char}(k) \nmid |G|$

```python
from invariant_theory import reynolds_operator
import sympy as sp

x1, x2 = sp.symbols('x1 x2')
# S₂ wirkt durch Variablenvertauschung
s2_elements = [(0,1), (1,0)]
def action(poly, perm):
    return poly.subs([(x1,x2),(x2,x1)], simultaneous=True) if perm==(1,0) else poly

result = reynolds_operator(x1, s2_elements, action)
# → (x1 + x2) / 2
```

---

## Elementarsymmetrische Polynome

Die **elementarsymmetrischen Polynome** $e_1, \ldots, e_n$ in $x_1, \ldots, x_n$:

$$e_1 = x_1 + x_2 + \cdots + x_n$$
$$e_2 = \sum_{i < j} x_i x_j$$
$$e_k = \sum_{i_1 < i_2 < \cdots < i_k} x_{i_1} x_{i_2} \cdots x_{i_k}$$
$$e_n = x_1 x_2 \cdots x_n$$

Sie erzeugen den gesamten Invariantenring unter $S_n$.

---

## Newton-Identitäten

Die **Newton-Identitäten** verbinden die Potenzsummen $p_k = \sum x_i^k$ mit den
elementarsymmetrischen Polynomen $e_k$:

$$p_k - e_1 p_{k-1} + e_2 p_{k-2} - \cdots + (-1)^{k-1} e_{k-1} p_1 + (-1)^k k e_k = 0$$

für $k \leq n$. Sie erlauben die Umrechnung zwischen beiden Basen des Invariantenrings.

---

## Molien-Reihe

Die **Molien-Reihe** (Poincaré-Reihe) eines Invariantenrings:

$$M(t) = \frac{1}{|G|} \sum_{g \in G} \frac{1}{\det(I - t \cdot g)}$$

Der Koeffizient von $t^k$ gibt die **Dimension** des Raums der $G$-invarianten Polynome
vom Grad $k$ an.

**Beispiel für $S_2$:**
$$M(t) = \frac{1}{(1-t)(1-t^2)}$$
Koeffizienten: $1, 1, 2, 2, 3, 3, \ldots$ (Anzahl der Basis-Invarianten pro Grad)

---

## Hilbert-Basissatz

> **Satz (Hilbert, 1890):** Für jede endliche Gruppe $G$ ist der Invariantenring
> $k[x_1, \ldots, x_n]^G$ endlich erzeugt als $k$-Algebra.

**Beispiel: $S_2$ auf $\mathbb{R}[x, y]$**

Der Invariantenring $\mathbb{R}[x,y]^{S_2}$ wird erzeugt von:
- $e_1 = x + y$ (Grad 1)
- $e_2 = xy$ (Grad 2)

Jedes $S_2$-invariante Polynom lässt sich als Polynom in $e_1, e_2$ schreiben:

$$x^2 + y^2 = e_1^2 - 2e_2$$
$$x^3 + y^3 = e_1^3 - 3e_1 e_2$$
$$(x-y)^2 = e_1^2 - 4e_2$$

---

## Diskriminante als Invariante

Die **Diskriminante** ist ein $S_n$-invariantes Polynom:

$$\Delta(f) = \prod_{i < j} (x_i - x_j)^2$$

wobei $x_1, \ldots, x_n$ die Wurzeln von $f$ sind.

- Quadratisch: $\Delta = b^2 - 4ac$
- Kubisch: $\Delta = 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2$
- $\Delta > 0$: alle Wurzeln reell und verschieden
- $\Delta = 0$: mindestens eine Doppelwurzel
- $\Delta < 0$: komplexe Wurzeln vorhanden (für reelles $f$)

---

## Hauptsatz der Theorie symmetrischer Polynome

> **Hauptsatz:** Jedes symmetrische Polynom in $x_1, \ldots, x_n$ ist eindeutig
> als Polynom in den elementarsymmetrischen Polynomen $e_1, \ldots, e_n$ darstellbar.

```python
result = fundamental_theorem_symmetric_poly('x1**2 + x2**2', 2)
# → {'is_symmetric': True, 'representation_in_e': 'e1**2 - 2*e2', 'verified': True}
```

---

## Wichtige Hinweis: Simultane Substitution

Bei SymPy-Variablenvertauschungen (z.B. $x_1 \leftrightarrow x_2$) muss
`simultaneous=True` verwendet werden, da `.subs()` standardmäßig sequenziell
arbeitet und Ketteneffekte entstehen:

```python
# Falsch: x1 → x2 → x1 (Kette)
poly.subs([(x1, x2), (x2, x1)])  # Gibt 2*x1 statt x1+x2!

# Korrekt: simultane Substitution
poly.subs([(x1, x2), (x2, x1)], simultaneous=True)  # Korrekt
```
