# model_theory.py — Modelltheorie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Sprache:** Python 3.13

---

## Überblick

Die **Modelltheorie** ist ein Teilgebiet der mathematischen Logik, das die Beziehung
zwischen formalen Sprachen/Theorien und ihren Interpretationen (Modellen) untersucht.
Das Modul `model_theory.py` implementiert die zentralen Konzepte objektorientiert
und dokumentiert die mathematischen Hintergründe ausführlich.

---

## 1. Signatur und Struktur

### Signatur

Eine **Signatur** (auch: Sprache, Vokabular) $L = (\mathcal{F}, \mathcal{R}, \mathcal{C})$ legt
die nicht-logischen Symbole fest:

- **Funktionssymbole** $\mathcal{F}$: $f^{(n)}$ mit Arität $n$ (z.B. $+ : A^2 \to A$)
- **Relationssymbole** $\mathcal{R}$: $R^{(n)} \subseteq A^n$ (z.B. $< \subseteq A^2$)
- **Konstantensymbole** $\mathcal{C}$: nullstellige Funktionen $c \in A$

### L-Struktur

Eine **L-Struktur** $\mathcal{M} = (A, F^\mathcal{M}, R^\mathcal{M}, C^\mathcal{M})$ besteht aus:

- **Universum** (Grundmenge) $A \neq \emptyset$
- **Funktionsinterpretation**: $f^\mathcal{M} : A^n \to A$ für jedes $f \in \mathcal{F}$
- **Relationsinterpretation**: $R^\mathcal{M} \subseteq A^n$ für jedes $R \in \mathcal{R}$
- **Konstanteninterpretation**: $c^\mathcal{M} \in A$ für jedes $c \in \mathcal{C}$

---

## 2. Tarski-Wahrheitsbedingungen

Die **rekursive Wahrheitsdefinition** (Tarski 1933) legt fest, wann eine Formel $\varphi$
in einer Struktur $\mathcal{M}$ unter einer Variablenbelegung $s : \text{Var} \to A$ gilt:

$$\mathcal{M} \models R(t_1,\ldots,t_n)[s] \iff (t_1^{\mathcal{M}}[s],\ldots,t_n^{\mathcal{M}}[s]) \in R^\mathcal{M}$$

$$\mathcal{M} \models \lnot\varphi[s] \iff \mathcal{M} \not\models \varphi[s]$$

$$\mathcal{M} \models (\varphi \land \psi)[s] \iff \mathcal{M} \models \varphi[s] \text{ und } \mathcal{M} \models \psi[s]$$

$$\mathcal{M} \models (\varphi \lor \psi)[s] \iff \mathcal{M} \models \varphi[s] \text{ oder } \mathcal{M} \models \psi[s]$$

$$\mathcal{M} \models \exists x\,\varphi[s] \iff \exists a \in A: \mathcal{M} \models \varphi[s(x \mapsto a)]$$

$$\mathcal{M} \models \forall x\,\varphi[s] \iff \forall a \in A: \mathcal{M} \models \varphi[s(x \mapsto a)]$$

### Formelnotation im Modul

Da Python keine Unicode-Operatoren direkt parst, verwendet das Modul eine
Präfix-Notation:

| Mathematisch | Modul-Notation |
|:------------|:--------------|
| $\lnot \varphi$ | `NOT(phi)` |
| $\varphi \land \psi$ | `AND(phi,psi)` |
| $\varphi \lor \psi$ | `OR(phi,psi)` |
| $\varphi \to \psi$ | `IMPLIES(phi,psi)` |
| $\varphi \leftrightarrow \psi$ | `IFF(phi,psi)` |
| $\exists x\,\varphi$ | `EXISTS(x,phi)` |
| $\forall x\,\varphi$ | `FORALL(x,phi)` |
| $R(t_1,\ldots,t_n)$ | `R(t1,...,tn)` |
| $t_1 = t_2$ | `=(t1,t2)` oder `t1=t2` |

---

## 3. Elementare Einbettung

Eine Abbildung $f: \mathcal{M} \to \mathcal{N}$ heißt **elementare Einbettung** gdw.
für alle $L$-Formeln $\varphi(x_1,\ldots,x_n)$ und alle $a_1,\ldots,a_n \in A$:

$$\mathcal{M} \models \varphi[a_1,\ldots,a_n] \iff \mathcal{N} \models \varphi[f(a_1),\ldots,f(a_n)]$$

Insbesondere gilt:
- Jede elementare Einbettung ist ein **injektiver Homomorphismus**
- Elementare Einbettung $\Rightarrow$ Einbettung $\Rightarrow$ Homomorphismus (nicht umgekehrt)

**Elementare Äquivalenz** $\mathcal{M} \equiv \mathcal{N}$: Beide Strukturen erfüllen
dieselben $L$-Sätze (Formeln ohne freie Variablen).

---

## 4. Kompaktheitssatz

$$\boxed{T \text{ hat ein Modell} \iff \text{jede endliche Teiltheorie } T_0 \subseteq T \text{ hat ein Modell}}$$

**Historisch:** Gödel (1930, für abzählbare Sprachen), Maltsev (1936, allgemein).

**Anwendung — Nicht-Standard-Arithmetik:**

Sei $T = \text{Th}(\mathbb{N})$ und $c$ eine neue Konstante. Definiere:

$$T' = T \cup \{c > 0,\; c > 1,\; c > 2,\;\ldots\}$$

Jede endliche Teilmenge von $T'$ hat ein Modell (nimm $\mathbb{N}$ und setze $c$ groß genug).
Nach dem Kompaktheitssatz hat $T'$ ein Modell $\mathcal{M}^*$, in dem $c^{\mathcal{M}^*}$
größer als jede Standardzahl ist — ein **nicht-standard Element**.

---

## 5. Löwenheim-Skolem-Sätze

### Abwärts-Löwenheim-Skolem (Löwenheim 1915, Skolem 1920)

$$\text{Hat } T \text{ ein unendliches Modell, so hat } T \text{ ein abzählbares Modell.}$$

Allgemein: Für jede unendliche Kardinalität $\kappa \leq |T| + \aleph_0$
hat $T$ ein Modell der Kardinalität $\kappa$.

**Beispiel:** Die Theorie DLO (Dense Linear Order without Endpoints) hat
$\mathbb{R}$ (Kardinalität $2^{\aleph_0}$) als Modell. Nach dem Abwärts-L-S
hat DLO auch $(\mathbb{Q}, <)$ als abzählbares Modell.

### Aufwärts-Löwenheim-Skolem (Tarski-Vaught 1936)

$$\text{Hat } T \text{ ein unendliches Modell der Kardinalität } \kappa, \text{ so hat } T \text{ Modelle jeder Kardinalität } \lambda \geq \kappa.$$

**Beweis-Idee:** Füge $\lambda$ neue Konstanten $c_\alpha$ hinzu mit Axiomen
$c_\alpha \neq c_\beta$ ($\alpha \neq \beta$). Kompaktheit liefert das Modell.

### Skolem-Paradoxon (1922)

ZF-Mengenlehre beweist die Existenz überabzählbarer Mengen.
Nach Abwärts-L-S hat ZF aber ein **abzählbares Modell** $\mathcal{M}$.

**Auflösung:** "Überabzählbar" ist relativ zum Modell:

$$X \text{ ist überabzählbar in } \mathcal{M} \iff \nexists f \in \mathcal{M}: f: \omega^\mathcal{M} \twoheadrightarrow X$$

Die Bijektion existiert *außerhalb* von $\mathcal{M}$, aber nicht als *Element* von $\mathcal{M}$.

---

## 6. Typentheorie

Ein **$n$-Typ** über $T$ ist eine Menge $p(x_1,\ldots,x_n)$ von $L$-Formeln mit
$T \cup p$ konsistent.

- **Vollständiger Typ**: $p$ ist maximal konsistent (für jede Formel $\varphi$:
  $\varphi \in p$ oder $\lnot\varphi \in p$)
- **Realisierung**: $\bar{a} \in A^n$ realisiert $p$ gdw. $\mathcal{M} \models \varphi[\bar{a}]$ für alle $\varphi \in p$
- **Auslassung**: $p$ wird ausgelassen wenn kein Tupel $p$ realisiert

### Typraum $S_n(T)$

Der Raum aller vollständigen $n$-Typen über $T$ trägt die **Stone-Topologie**:

$$[{\varphi}] = \{p \in S_n(T) \mid \varphi \in p\}$$

Eigenschaften von $S_n(T)$:
- **Kompakt** (Kompaktheitssatz)
- **Hausdorff** (totalgetrennt)
- **Perfekt** für vollständige $T$ ohne isolierte Typen

### Auslassungstypen-Satz (Henkin 1954)

Sei $T$ vollständig abzählbar und $p(x)$ ein **nicht-isolierter** Typ
(kein $\theta$ mit $T \models \theta \to \varphi$ für alle $\varphi \in p$ und $T \models \exists x\,\theta$).
Dann existiert ein abzählbares Modell $\mathcal{M} \models T$, das $p$ **auslässt**.

---

## 7. κ-Kategorik und Vollständigkeit

### κ-Kategorik

Eine Theorie $T$ heißt **$\kappa$-kategorisch** gdw. alle Modelle von $T$ der
Kardinalität $\kappa$ isomorph sind (genau ein Modell bis auf Isomorphie).

| Theorie | $\aleph_0$-kategorisch | überabzählbar-kategorisch |
|:--------|:----------------------:|:-------------------------:|
| DLO     | ✓ (Cantor 1895)        | ✗                         |
| ACF$_p$ | ✗                      | ✓ (für alle $\lambda > \aleph_0$) |
| PA      | ✗                      | ✗                         |

### Vaught-Kriterium (1954)

$$T \text{ ohne endliche Modelle, } \kappa\text{-kategorisch für ein unendliches } \kappa \implies T \text{ vollständig.}$$

**Beweis:** Sind $\mathcal{M}_1, \mathcal{M}_2 \models T$ unendlich, so existieren
elementare Erweiterungen $\mathcal{N}_i$ der Kardinalität $\kappa$. Dann
$\mathcal{N}_1 \cong \mathcal{N}_2$ (Kategorik) $\Rightarrow \text{Th}(\mathcal{M}_1) = \text{Th}(\mathcal{M}_2)$.

### Morleys Kategorien-Satz (1965)

$$\boxed{T \text{ abzählbar, } \kappa\text{-kategorisch für ein überabzählbares } \kappa \implies T \; \lambda\text{-kategorisch für alle überabzählbaren } \lambda.}$$

Dies ist eines der tiefsten Ergebnisse der Modelltheorie und motivierte
Shelahs **Stabilitätstheorie** (Classification Theory, 1978).

**Stabilitätshierarchie (Shelah):**
- **$\omega$-stabil**: $|S_n(A)| \leq \aleph_0$ für alle abzählbaren $A$
- **Superstabil**: stabil in allen $\lambda \geq 2^{\aleph_0}$
- **Stabil**: stabil in bestimmten Kardinalitäten
- **Instabil**: nicht stabil (z.B. DLO, PA)

*Alle überabzählbar-kategorialen Theorien sind $\omega$-stabil.*

---

## 8. Quantorenelimination

Eine Theorie $T$ hat **Quantorenelimination (QE)** gdw. für jede $L$-Formel $\varphi$
eine quantorenfreie Formel $\psi$ existiert mit:

$$T \models \varphi \leftrightarrow \psi$$

### Wichtige QE-Theorien

**DLO** (Langford 1927, Back-and-Forth):
$$T \models \exists y\,(x_1 < y \land y < x_2) \;\leftrightarrow\; x_1 < x_2$$

**ACF** (Tarski 1948, Resultanten):
$$\text{ACF} \models \exists y\, p(y) = 0 \;\leftrightarrow\; \Delta(p) \neq 0 \text{ (für quadr. } p\text{)}$$

wobei $\Delta$ die Diskriminante bezeichnet.

**RCF** (Tarski 1951):
$$\text{RCF} \models \exists y\,(y^2 = x) \;\leftrightarrow\; x \geq 0$$

Tarskis QE für RCF beweist die **Entscheidbarkeit der reellen Arithmetik**
(im Gegensatz zur unentscheidbaren Arithmetik der natürlichen Zahlen, Gödel).

**Presburger-Arithmetik** $\text{Th}(\mathbb{Z}, +, <)$ (1929):
$$\text{Pres} \models \exists y\,(x = 2 \cdot y) \;\leftrightarrow\; 2 \mid x$$

---

## 9. API-Übersicht

### Klassen

| Klasse | Beschreibung |
|:-------|:------------|
| `Signature` | Formale Signatur $L = (\mathcal{F}, \mathcal{R}, \mathcal{C})$ |
| `Structure` | L-Struktur mit Tarski-Wahrheitsbedingungen |
| `ElementaryEmbedding` | Elementare Einbettung $f: \mathcal{M} \to \mathcal{N}$ |
| `Type` | p-Typ mit Vollständigkeit, Realisierung, Auslassung |

### Funktionen

| Funktion | Beschreibung |
|:---------|:------------|
| `elementary_equivalence(M, N, sentences)` | $\mathcal{M} \equiv \mathcal{N}$ |
| `is_definable(structure, subset, formula)` | Definierbarkeit $X = \varphi^\mathcal{M}$ |
| `compactness_theorem_demo()` | Kompaktheitssatz-Demo (Nicht-Standard-Arith.) |
| `nonstandard_arithmetic_demo()` | Nicht-Standard-Modell-Demo |
| `lowenheim_skolem_downward(axioms, size)` | Abwärts-L-S (DLO/ℚ-Demo) |
| `lowenheim_skolem_upward(axioms, size)` | Aufwärts-L-S |
| `skolem_paradox_explanation()` | Skolem-Paradoxon mit Auflösung |
| `type_space(theory, variables)` | Typraum $S_n(T)$ |
| `omitting_types_theorem_demo()` | Auslassungstypen-Satz |
| `theory_is_categorical(name, card)` | $\kappa$-Kategorik |
| `vaught_test(theory_name)` | Vaught-Vollständigkeitskriterium |
| `morley_theorem_demo()` | Morleys Kategorien-Satz |
| `quantifier_elimination_demo(theory)` | QE für DLO, ACF, RCF, Presburger |

---

## 10. Beispiele

### Struktur erstellen und Formel auswerten

```python
from src.model_theory import Signature, Structure

sig = Signature(functions={}, relations={'LT': 2}, constants=[])
universe = [0, 1, 2, 3, 4]
lt_rel = {(a, b) for a in universe for b in universe if a < b}

M = Structure(
    universe=universe,
    signature=sig,
    functions={},
    relations={'LT': lt_rel},
    constants={}
)

# Tarski: ∀x ¬(x < x) — Irreflexivität
print(M.satisfies('FORALL(x,NOT(LT(x,x)))'))  # True

# ∃y (y < 3) — es gibt etwas kleiner als 3
print(M.satisfies('EXISTS(y,LT(y,3))'))  # True
```

### Kompaktheitssatz-Demo

```python
from src.model_theory import compactness_theorem_demo

result = compactness_theorem_demo()
print(result['theorem'])      # Kompaktheitssatz
print(result['conclusion'])   # Nicht-Standard-Modell garantiert
```

### DLO — ℵ₀-Kategorik prüfen

```python
from src.model_theory import theory_is_categorical

info = theory_is_categorical('DLO', 'aleph_0')
print(info['is_categorical'])  # True
print(info['unique_model'])    # (ℚ, <)
```

### Quantorenelimination in RCF

```python
from src.model_theory import quantifier_elimination_demo

qe = quantifier_elimination_demo('RCF')
for ex in qe['examples']:
    print(f"Mit Quantor:  {ex['with_quantifier']}")
    print(f"Quantorenfrei: {ex['quantifier_free']}")
```

---

## 11. Weiterführende Literatur

- **Hodges, W.** (1993). *Model Theory*. Cambridge University Press.
- **Chang, C. C. & Keisler, H. J.** (1990). *Model Theory* (3rd ed.). North-Holland.
- **Marker, D.** (2002). *Model Theory: An Introduction*. Springer.
- **Shelah, S.** (1978). *Classification Theory*. North-Holland.
- **Tarski, A.** (1951). *A Decision Method for Elementary Algebra and Geometry*.
