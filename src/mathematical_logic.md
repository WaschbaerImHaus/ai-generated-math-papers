# mathematical_logic.py – Dokumentation

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Datei:** `src/mathematical_logic.py`

---

## Übersicht

Das Modul `mathematical_logic.py` implementiert Kernbereiche der mathematischen
und formalen Logik: Aussagenlogik, Resolution, SAT-Solving, Prädikatenlogik
erster Stufe, Beweiskalküle, Modallogik und Gödels Unvollständigkeitssätze.

---

## 1. Aussagenlogik (Propositional Logic)

### 1.1 Proposition

Atomare logische Aussage mit Namen und optionalem Wahrheitswert.

```
p = Proposition("p", True)
p.evaluate({"p": True})  # → True
```

### 1.2 LogicFormula

Zusammengesetzte Formel mit folgenden Operatoren:

| Operator | Symbol | Bedeutung |
|----------|--------|-----------|
| `AND`    | ∧      | Konjunktion |
| `OR`     | ∨      | Disjunktion |
| `NOT`    | ¬      | Negation |
| `IMPLIES`| →      | Implikation |
| `IFF`    | ↔      | Äquivalenz |
| `XOR`    | ⊕      | Exklusives Oder |
| `NAND`   | ↑      | Nicht-Und |
| `NOR`    | ↓      | Nicht-Oder |
| `ATOM`   | –      | Atomare Variable |

**Wahrheitstabellen der Konnektive:**

| p | q | p∧q | p∨q | p→q | p↔q | p⊕q |
|---|---|-----|-----|-----|-----|-----|
| T | T |  T  |  T  |  T  |  T  |  F  |
| T | F |  F  |  T  |  F  |  F  |  T  |
| F | T |  F  |  T  |  T  |  F  |  T  |
| F | F |  F  |  F  |  T  |  T  |  F  |

**Wichtige Äquivalenzen:**

$$p \rightarrow q \equiv \neg p \vee q$$

$$p \leftrightarrow q \equiv (p \rightarrow q) \wedge (q \rightarrow p)$$

$$p \oplus q \equiv (p \vee q) \wedge \neg(p \wedge q)$$

### 1.3 De Morganschen Gesetze

$$\neg(A \wedge B) \equiv \neg A \vee \neg B$$

$$\neg(A \vee B) \equiv \neg A \wedge \neg B$$

### 1.4 Normalformen

**Konjunktive Normalform (CNF):**

Eine Formel ist in CNF, wenn sie eine Konjunktion von Klauseln ist,
wobei jede Klausel eine Disjunktion von Literalen ist:

$$(l_{11} \vee \cdots \vee l_{1k_1}) \wedge (l_{21} \vee \cdots \vee l_{2k_2}) \wedge \cdots$$

Konversion: Implikationen eliminieren → NNF → AND über OR verteilen.

**Disjunktive Normalform (DNF):**

Eine Formel ist in DNF, wenn sie eine Disjunktion von Konjunktionen ist:

$$(l_{11} \wedge \cdots \wedge l_{1k_1}) \vee (l_{21} \wedge \cdots \wedge l_{2k_2}) \vee \cdots$$

### 1.5 Semantische Konzepte

**Tautologie:** Eine Formel ist wahr unter allen Belegungen.
$$\models \varphi \iff \forall \beta: \beta(\varphi) = \top$$

Beispiel: $p \vee \neg p$ (Satz vom ausgeschlossenen Dritten).

**Kontradiktion:** Eine Formel ist falsch unter allen Belegungen.
$$\models \neg \varphi \iff \forall \beta: \beta(\varphi) = \bot$$

Beispiel: $p \wedge \neg p$.

**Erfüllbarkeit:** Es existiert mindestens eine wahre Belegung.
$$\exists \beta: \beta(\varphi) = \top$$

**Logische Konsequenz:**
$$\Gamma \models \varphi \iff \text{jedes Modell von } \Gamma \text{ ist auch Modell von } \varphi$$

**Semantische Äquivalenz:**
$$\varphi \equiv \psi \iff \models (\varphi \leftrightarrow \psi)$$

---

## 2. Resolution und SAT-Solving

### 2.1 Resolutionsregel

Die Resolutionsregel für Klauseln:

$$\frac{A \vee l \quad B \vee \neg l}{A \vee B}$$

Literale werden als Strings dargestellt: `"p"` für $p$, `"-p"` für $\neg p$.

```python
c1 = frozenset({"p", "q"})
c2 = frozenset({"-p", "r"})
result = resolution_step(c1, c2)  # frozenset({"q", "r"})
```

### 2.2 Resolutionswiderlegung

Eine Klauselmenge ist unerfüllbar, wenn die leere Klausel $\square$ ableitbar ist:

$$\Gamma \vdash \square \iff \Gamma \text{ unerfüllbar}$$

Der Algorithmus sättigt die Klauselmenge durch alle möglichen Resolutionsschritte.

### 2.3 DPLL-Algorithmus

Davis-Putnam-Logemann-Loveland (1962):

1. **Unit-Propagation:** Erzwinge Belegung für Einheitsklauseln $\{l\}$
2. **Reine Literal-Elimination:** Literale die nur positiv/negativ vorkommen
3. **Branching:** Wähle Variable, probiere True und False (Backtracking)

$$\text{DPLL}(S) = \begin{cases} \top & \text{wenn } S = \emptyset \\ \bot & \text{wenn } \square \in S \\ \text{DPLL}(S[x:=T]) \vee \text{DPLL}(S[x:=F]) & \text{sonst} \end{cases}$$

### 2.4 Tseitin-Transformation

Wandelt beliebige Formeln in equisatisfiable CNF um (polynomial in der Formelgröße).

Für jede Teilformel $\varphi_i$ wird eine Hilfsvariable $t_i$ eingeführt:

- $t \leftrightarrow (A \wedge B)$: $\{(\neg t \vee A), (\neg t \vee B), (t \vee \neg A \vee \neg B)\}$
- $t \leftrightarrow (A \vee B)$: $\{(\neg t \vee A \vee B), (t \vee \neg A), (t \vee \neg B)\}$
- $t \leftrightarrow \neg A$: $\{(t \vee A), (\neg t \vee \neg A)\}$

---

## 3. Prädikatenlogik (First-Order Logic)

### 3.1 Terme

Terme sind induktiv definiert:
- Variablen: $x, y, z$
- Konstanten: $a, b, c$
- Funktionsanwendung: $f(t_1, \ldots, t_n)$

```python
x = Term.variable("x")
a = Term.constant("a")
fx = Term.function("f", x)   # f(x)
```

### 3.2 Prädikate und FOL-Formeln

$$P(t_1, \ldots, t_n), \quad \forall x.\varphi, \quad \exists x.\varphi$$

**Freie und gebundene Variablen:**

In $\forall x.P(x, y)$ ist $x$ gebunden und $y$ frei.

Eine geschlossene Formel (Satz) hat keine freien Variablen.

**Substitution:** $\varphi[x/t]$ ersetzt freie Vorkommen von $x$ durch $t$.

### 3.3 Pränex-Normalform (PNF)

Alle Quantoren stehen vorne:

$$Q_1 x_1 Q_2 x_2 \cdots Q_n x_n . \psi$$

Regeln zum Vorziehen von Quantoren:
- $\neg \forall x.\varphi \equiv \exists x.\neg \varphi$
- $\neg \exists x.\varphi \equiv \forall x.\neg \varphi$
- $(\forall x.\varphi) \wedge \psi \equiv \forall x.(\varphi \wedge \psi)$ (wenn $x \notin FV(\psi)$)

### 3.4 Skolemisierung

Eliminiert Existenzquantoren durch Skolem-Funktionen:

$$\forall x_1 \cdots \forall x_n . \exists y . \varphi \rightsquigarrow \forall x_1 \cdots \forall x_n . \varphi[y / f(x_1, \ldots, x_n)]$$

Ohne vorangehende Allquantoren: Skolem-Konstante $c$ statt Funktion.

### 3.5 Herbrand-Universum

Das Herbrand-Universum $H_\Sigma$ einer Signatur enthält alle Grundterme:

$$H_0 = \{c \mid c \text{ Konstante}\}$$
$$H_{i+1} = H_i \cup \{f(t_1, \ldots, t_n) \mid f \text{ n-stellig}, t_i \in H_i\}$$
$$H_\Sigma = \bigcup_{i \geq 0} H_i$$

---

## 4. Beweiskalküle

### 4.1 Natürliches Schließen (Gentzen, 1935)

Inferenzregeln (Auswahl):

**Einführungsregeln:**

$$\frac{A \quad B}{A \wedge B} \wedge I \qquad \frac{A}{A \vee B} \vee I_1 \qquad \frac{[A] \vdots B}{A \rightarrow B} \rightarrow I$$

**Eliminierungsregeln:**

$$\frac{A \wedge B}{A} \wedge E_1 \qquad \frac{A \rightarrow B \quad A}{B} \rightarrow E \qquad \frac{A \quad \neg A}{\bot} \neg E$$

**Modus Ponens** (→E):
$$\frac{A \rightarrow B \quad A}{B}$$

**Modus Tollens:**
$$\frac{A \rightarrow B \quad \neg B}{\neg A}$$

### 4.2 Hilbert-Kalkül

Axiomenschemata (für beliebige Formeln $A, B, C$):

- **K-Axiom:** $A \rightarrow (B \rightarrow A)$
- **S-Axiom:** $(A \rightarrow (B \rightarrow C)) \rightarrow ((A \rightarrow B) \rightarrow (A \rightarrow C))$
- **Doppelnegation:** $\neg\neg A \rightarrow A$

Einzige Schlussregel: Modus Ponens.

---

## 5. Modallogik

### 5.1 Kripke-Semantik

Ein **Kripke-Modell** $M = (W, R, V)$ besteht aus:
- $W$: Menge möglicher Welten
- $R \subseteq W \times W$: Erreichbarkeitsrelation
- $V$: Valuation $V(p) \subseteq W$

**Semantik der Modaloperatoren:**

$$M, w \models \Box A \iff \forall v: (wRv \Rightarrow M, v \models A)$$

$$M, w \models \Diamond A \iff \exists v: (wRv \wedge M, v \models A)$$

Beziehung: $\Diamond A \equiv \neg \Box \neg A$

### 5.2 Modalsysteme und Rahmen-Eigenschaften

| System | Zusatzaxiom | Eigenschaft von R |
|--------|-------------|-------------------|
| K      | –           | beliebig |
| D      | $\Box A \rightarrow \Diamond A$ | seriell |
| T      | $\Box A \rightarrow A$ | reflexiv |
| B      | $A \rightarrow \Box\Diamond A$ | reflexiv + symmetrisch |
| S4     | $\Box A \rightarrow \Box\Box A$ | reflexiv + transitiv |
| S5     | $\Diamond A \rightarrow \Box\Diamond A$ | reflexiv + transitiv + euklidisch |

**Korrespondenztheorem:** Axiom T korrespondiert mit Reflexivität,
Axiom 4 mit Transitivität, Axiom 5 mit Euklidizität.

---

## 6. Gödels Unvollständigkeitssätze

### 6.1 Gödel-Nummerierung

Jedem Symbol wird eine Zahl zugeordnet. Eine Formel $\langle s_1, s_2, \ldots, s_n \rangle$ erhält die Gödel-Zahl:

$$\lceil \varphi \rceil = 2^{c(s_1)} \cdot 3^{c(s_2)} \cdot 5^{c(s_3)} \cdots p_n^{c(s_n)}$$

wobei $p_i$ die $i$-te Primzahl und $c(s_i)$ der Code des Symbols ist.

### 6.2 Erster Unvollständigkeitssatz (Gödel, 1931)

**Satz:** Sei $F$ ein konsistentes, rekursiv aufzählbares formales System,
das die Peano-Arithmetik ausdrücken kann. Dann gibt es eine Aussage $G$, die:
- wahr ist (in der Standardinterpretation)
- in $F$ nicht beweisbar ist
- in $F$ nicht widerlegbar ist

**Konstruktion (Diagonallemma):**

Es gibt eine Aussage $G$ mit $F \vdash G \leftrightarrow \neg \text{Bew}_F(\lceil G \rceil)$.

$G$ sagt: „Ich bin in $F$ nicht beweisbar."

**Beweis durch Widerspruch:**
- Ist $G$ beweisbar in $F$: Dann gilt $\text{Bew}_F(\lceil G \rceil)$, also $\neg G$ – Widerspruch.
- Ist $\neg G$ beweisbar: Dann gilt $\text{Bew}_F(\lceil G \rceil)$, also $G$ – Widerspruch.

### 6.3 Zweiter Unvollständigkeitssatz

**Satz:** Kein konsistentes formales System $F$ kann seine eigene Konsistenz $\text{Con}(F)$ beweisen.

$$F \not\vdash \text{Con}(F) \quad (\text{sofern } F \text{ konsistent})$$

---

## API-Referenz (Kurzübersicht)

```python
# Aussagenlogik
Proposition(name, value=None)
LogicFormula(op, *args)
LogicFormula.atom(name)
formula.evaluate(assignment)
formula.to_cnf()
formula.to_dnf()
formula.simplify()
truth_table(formula, variables)
is_tautology(formula, variables)
is_satisfiable(formula, variables)
is_contradiction(formula, variables)
logical_consequence(premises, conclusion, variables)
logical_equivalence(f1, f2, variables)

# Resolution / SAT
resolution_step(clause1, clause2)
resolution_refutation(clauses)
dpll(clauses, assignment=None)
tseitin_transform(formula)

# Prädikatenlogik
Term.variable(name), Term.constant(name), Term.function(name, *args)
Predicate(name, arity)
FOLFormula.forall(var, formula), FOLFormula.exists(var, formula)
prenex_normal_form(formula)
skolemize(formula)
herbrand_universe(constants, function_symbols, depth=2)

# Beweiskalküle
NaturalDeductionProof()
HilbertSystem()
modus_ponens(major, minor)
modus_tollens(major, minor)

# Modallogik
KripkeFrame(worlds, accessibility)
ModalFormula.box(formula), ModalFormula.diamond(formula)
modal_system(system)

# Gödel
godel_numbering(formula_str)
godel_decode(number)
incompleteness_demonstration()
```
