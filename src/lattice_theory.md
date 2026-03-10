# lattice_theory.py — Verbandstheorie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Modul:** `src/lattice_theory.py`

---

## Überblick

Das Modul implementiert die grundlegenden Strukturen der **Verbandstheorie** (englisch: *Lattice Theory*):
Halbordnungen, Verbände und Boolesche Algebren sowie wichtige Sätze wie den
Birkhoff-Darstellungssatz und die Dedekind-Zahlen.

---

## Mathematische Grundlagen

### Halbordnung (Partial Order)

Eine **Halbordnung** $(M, \leq)$ ist eine Menge $M$ mit einer binären Relation $\leq$, die folgende Axiome erfüllt:

$$\text{Reflexivität:} \quad a \leq a$$
$$\text{Antisymmetrie:} \quad a \leq b \wedge b \leq a \Rightarrow a = b$$
$$\text{Transitivität:} \quad a \leq b \wedge b \leq c \Rightarrow a \leq c$$

**Hasse-Diagramm:** Die reduzierte Darstellung einer Halbordnung. Kante $a \to b$ nur wenn $a < b$ und kein $c$ mit $a < c < b$ existiert.

---

### Verband (Lattice)

Ein **Verband** $(L, \leq)$ ist eine Halbordnung, in der jedes Paar $(a, b)$ besitzt:

- **Infimum (Meet):** $a \wedge b = \max\{c \mid c \leq a \wedge c \leq b\}$
- **Supremum (Join):** $a \vee b = \min\{c \mid a \leq c \wedge b \leq c\}$

**Verbandsaxiome (äquivalente Formulierung):**

$$a \wedge a = a, \quad a \vee a = a \quad \text{(Idempotenz)}$$
$$a \wedge b = b \wedge a, \quad a \vee b = b \vee a \quad \text{(Kommutativität)}$$
$$a \wedge (b \wedge c) = (a \wedge b) \wedge c \quad \text{(Assoziativität)}$$
$$a \wedge (a \vee b) = a, \quad a \vee (a \wedge b) = a \quad \text{(Absorption)}$$

---

### Distributiver Verband

Ein Verband heißt **distributiv**, wenn gilt:

$$a \wedge (b \vee c) = (a \wedge b) \vee (a \wedge c)$$

(Äquivalent: $a \vee (b \wedge c) = (a \vee b) \wedge (a \vee c)$)

**Klassifikation:** $L$ ist distributiv $\Leftrightarrow$ $L$ enthält keinen Teilverband isomorph zu $M_3$ (Diamant) oder $N_5$ (Pentagon).

---

### Modularer Verband

Ein Verband heißt **modular**, wenn gilt:

$$a \leq c \Rightarrow a \vee (b \wedge c) = (a \vee b) \wedge c$$

Jeder distributive Verband ist modular, aber der Diamant $M_3$ ist modular ohne distributiv zu sein.

---

### Komplementierter Verband

Ein beschränkter Verband mit $\bot, \top$ heißt **komplementiert**, wenn:

$$\forall a \in L \; \exists \bar{a}: \quad a \wedge \bar{a} = \bot \wedge a \vee \bar{a} = \top$$

---

### Boolesche Algebra

Eine **Boolesche Algebra** ist ein distributiver, komplementierter, beschränkter Verband:

$$\text{Idempotenz}: \quad a \wedge a = a, \quad a \vee a = a$$
$$\text{Komplementarität}: \quad a \wedge \bar{a} = \bot, \quad a \vee \bar{a} = \top$$
$$\text{De Morgan}: \quad \overline{a \wedge b} = \bar{a} \vee \bar{b}, \quad \overline{a \vee b} = \bar{a} \wedge \bar{b}$$

**Standardbeispiele:**
- $\mathcal{P}(S)$ mit $\subseteq$, $\cap$, $\cup$, $\emptyset$, $S$
- $\{0,1\}^n$ mit bitweisem AND/OR
- Aussagenlogik

---

### Atome und Koatome

In einer Booleschen Algebra mit $\bot$:
- **Atom:** $a \neq \bot$ mit $\bot < a$ und kein $b$ mit $\bot < b < a$
- **Koatom:** $a \neq \top$ mit $a < \top$ und kein $b$ mit $a < b < \top$

---

## Birkhoff-Darstellungssatz

**Satz (Birkhoff, 1937):** Jeder endliche distributive Verband $L$ ist isomorph zur Menge der **Ordnungsideale** von $J_{\text{irr}}(L)$:

$$L \cong \mathcal{J}(J_{\text{irr}}(L))$$

wobei $J_{\text{irr}}(L) = \{j \in L \mid j \neq \bot,\; j = a \vee b \Rightarrow j = a \text{ oder } j = b\}$ die **join-irreduziblen** Elemente sind.

Ein **Ordnungsideal** $I \subseteq P$ ist nach unten abgeschlossen: $b \leq a \in I \Rightarrow b \in I$.

---

## Dedekind-Zahlen

Die **Dedekind-Zahl** $D(n)$ zählt die monotonen Booleschen Funktionen auf $n$ Variablen (äquivalent: Antiketten in $\mathcal{P}(\{1,\ldots,n\})$):

| $n$ | $D(n)$ |
|-----|--------|
| 0 | 2 |
| 1 | 3 |
| 2 | 6 |
| 3 | 20 |
| 4 | 168 |
| 5 | 7 581 |
| 6 | 7 828 354 |
| 7 | 2 414 682 040 998 |
| 8 | 56 130 437 228 687 557 907 788 |

$D(8)$ wurde 2023 berechnet (Dedekind-Problem gelöst nach über 40 Jahren).

---

## Jordan-Dedekind-Kettenbedingung

Ein Verband erfüllt die **Jordan-Dedekind-Kettenbedingung (JD)**, wenn alle maximalen Ketten zwischen je zwei Elementen $a \leq b$ gleich lang sind. Dies ist äquivalent zur Existenz einer **Rangfunktion**:

$$r: L \to \mathbb{N}_0 \quad \text{mit} \quad a \lessdot b \Rightarrow r(b) = r(a) + 1$$

(wobei $a \lessdot b$ bedeutet: $b$ überdeckt $a$ direkt).

---

## Wichtige Beispiele

### Teilerverband von $n$

Elemente: $\{d \mid d \mid n\}$, Ordnung: $d_1 \leq d_2 \Leftrightarrow d_1 \mid d_2$

$$d_1 \wedge d_2 = \gcd(d_1, d_2), \quad d_1 \vee d_2 = \mathrm{lcm}(d_1, d_2)$$

### Potenzmengenverband $\mathcal{P}(n)$

$$A \wedge B = A \cap B, \quad A \vee B = A \cup B, \quad \bot = \emptyset, \quad \top = \{0,\ldots,n-1\}$$

### Partitionsverband $\Pi_n$

Elemente: alle Partitionen von $\{1,\ldots,n\}$, geordnet nach **Verfeinerung**.
$|\Pi_n| = B(n)$ (Bell-Zahl). Beispiel: $|\Pi_3| = 5$, $|\Pi_4| = 15$.

---

## Stone-Dualität (endlicher Fall)

**Satz (Stone):** Jede endliche Boolesche Algebra $B$ ist isomorph zur Potenzmenge ihrer Atome:

$$B \cong \mathcal{P}(\mathrm{At}(B))$$

Insbesondere hat jede endliche Boolesche Algebra $2^k$ Elemente (für $k = |\mathrm{At}(B)|$).

---

## API-Referenz

### Klasse `PartialOrder`

| Methode | Beschreibung |
|---------|-------------|
| `is_partial_order()` | Prüft alle drei Axiome |
| `leq(a, b)` | $a \leq b$? |
| `hasse_diagram()` | Reduzierte Überdeckungsrelation |
| `minimal_elements()` | Elemente ohne kleineres |
| `maximal_elements()` | Elemente ohne größeres |
| `chain_length()` | Länge der längsten Kette |
| `is_chain()` | Total geordnet? |

### Klasse `Lattice(PartialOrder)`

| Methode | Beschreibung |
|---------|-------------|
| `meet(a, b)` | $a \wedge b$ |
| `join(a, b)` | $a \vee b$ |
| `top()` | $\top$ oder `None` |
| `bottom()` | $\bot$ oder `None` |
| `is_distributive()` | $a \wedge (b \vee c) = (a \wedge b) \vee (a \wedge c)$ |
| `is_modular()` | Modularitätsgesetz |
| `is_complemented()` | Jedes $a$ hat $\bar{a}$ |
| `complement(a)` | Findet $\bar{a}$ |

### Klasse `BooleanAlgebra(Lattice)`

| Methode | Beschreibung |
|---------|-------------|
| `is_boolean_algebra()` | Alle Axiome |
| `atoms()` | Elemente direkt über $\bot$ |
| `coatoms()` | Elemente direkt unter $\top$ |
| `stone_representation()` | $B \cong \mathcal{P}(\mathrm{At}(B))$ |

### Freie Funktionen

| Funktion | Beschreibung |
|----------|-------------|
| `divisibility_lattice(n)` | Teilerverband von $n$ |
| `power_set_lattice(n)` | $\mathcal{P}(\{0,\ldots,n-1\})$ |
| `partition_lattice(n)` | Partitionsverband $\Pi_n$ |
| `birkhoff_representation_theorem(lat)` | Birkhoff-Satz |
| `dedekind_numbers(n)` | $D(n)$ für $n \leq 6$ |
| `jordan_dedekind_chain_condition(lat)` | JD-Bedingung |
