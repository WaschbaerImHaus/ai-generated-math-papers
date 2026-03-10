# proof_theory_formal.py — Formale Beweistheorie

**Autor:** Kurt Ingwer
**Version:** 1.0
**Letzte Änderung:** 2026-03-10
**Abhängigkeiten:** Python 3.13, keine externen Bibliotheken

---

## Überblick

Das Modul implementiert Kernkonzepte der formalen Beweistheorie, die auf Gentzen (1934/1936) zurückgehen. Es deckt vier große Bereiche ab:

1. **Sequenzenkalkül LK** — formale Beweissysteme mit expliziten Regelanwendungen
2. **Natürliches Schließen (ND)** — beweisnaher intuitiver Stil
3. **Beweiskomplexität** — Vergleich der Stärke verschiedener Kalküle
4. **Ordinale Analyse & Reverse Mathematics** — quantitative Messung von Systemstärken

---

## Teil 1: Sequenzenkalkül LK (Gentzen 1934)

### Grundbegriffe

Ein **Sequent** hat die Form:

$$\Gamma \vdash \Delta$$

wobei $\Gamma$ (Antecedent) und $\Delta$ (Succedent) endliche Listen von Formeln sind. Die Bedeutung:

> Wenn alle Formeln in $\Gamma$ wahr sind, dann ist mindestens eine Formel in $\Delta$ wahr.

Formal ist $\Gamma \vdash \Delta$ äquivalent zu:

$$\neg(\Gamma_1 \wedge \cdots \wedge \Gamma_n) \vee \Delta_1 \vee \cdots \vee \Delta_m$$

### Klasse `Sequent`

```python
s = Sequent(["A", "B"], ["C"])
print(s)          # "A, B ⊢ C"
s.is_axiom()      # True wenn A∩D ≠ ∅ (Identitäts-Axiom)
```

Das **Identitäts-Axiom** (Init-Regel):

$$\overline{\Gamma, A \vdash A, \Delta} \quad (\text{Axiom})$$

### Klasse `ProofTree`

Beweisbäume sind rekursive Strukturen. Jeder Knoten enthält:
- `conclusion` — das bewiesene Sequent
- `rule` — die angewendete Regel
- `premises` — Liste der Teilbäume (Prämissen)

Ein Beweis ist **gültig** (`is_valid()`), wenn alle Blätter Axiome sind:

```python
blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
blatt.is_valid()  # True
```

### LK-Regeln (`LKRules`)

#### Strukturregeln

| Regel | Notation | Beschreibung |
|-------|----------|--------------|
| Weakening Left | $\frac{\Gamma \vdash \Delta}{A, \Gamma \vdash \Delta}$ WL | Formel im Antecedent hinzufügen |
| Weakening Right | $\frac{\Gamma \vdash \Delta}{\Gamma \vdash \Delta, A}$ WR | Formel im Succedent hinzufügen |
| Contraction Left | $\frac{A, A, \Gamma \vdash \Delta}{A, \Gamma \vdash \Delta}$ CL | Dopplung im Antecedent entfernen |
| Contraction Right | $\frac{\Gamma \vdash \Delta, A, A}{\Gamma \vdash \Delta, A}$ CR | Dopplung im Succedent entfernen |

#### Logische Regeln

**Konjunktion:**

$$\frac{\Gamma \vdash \Delta, A \quad \Gamma \vdash \Delta, B}{\Gamma \vdash \Delta, A \wedge B} \wedge\text{R} \qquad \frac{A, B, \Gamma \vdash \Delta}{A \wedge B, \Gamma \vdash \Delta} \wedge\text{L}$$

**Implikation:**

$$\frac{A, \Gamma \vdash \Delta, B}{\Gamma \vdash \Delta, A \to B} {\to}\text{R} \qquad \frac{\Gamma \vdash \Delta, A \quad B, \Sigma \vdash \Pi}{A \to B, \Gamma, \Sigma \vdash \Delta, \Pi} {\to}\text{L}$$

**Negation:**

$$\frac{\Gamma \vdash \Delta, A}{\neg A, \Gamma \vdash \Delta} \neg\text{L} \qquad \frac{A, \Gamma \vdash \Delta}{\Gamma \vdash \Delta, \neg A} \neg\text{R}$$

#### Schnitt-Regel (Cut)

$$\frac{\Gamma \vdash \Delta, A \quad A, \Sigma \vdash \Pi}{\Gamma, \Sigma \vdash \Delta, \Pi} \text{Cut}$$

Die Cut-Formel $A$ tritt in der Konklusion nicht mehr auf. Der **Hauptsatz** (Gentzen 1934) besagt, dass diese Regel eliminierbar ist.

---

## Teil 2: Der Hauptsatz (Schnittelimination)

### Satz (Gentzen 1934)

> Jeder LK-Beweis kann in einen **schnittfreien** Beweis umgewandelt werden.

```python
demo = cut_elimination_demo()
# demo["hauptsatz"]            — Beschreibung
# demo["beweis_mit_cut"]       — Beweis mit Cut
# demo["schnittfreier_beweis"] — Schnittfreier Beweis
# demo["folgerungen"]          — Konsequenzen des Satzes
```

### Folgerungen

**Teilformel-Eigenschaft (Subformula Property):** In schnittfreien Beweisen kommen in jedem Sequent nur **Teilformeln der Schlusssequenz** vor. Dies macht Beweise analysierbar und erklärt, warum Schnittelimination wichtig ist.

**Konsistenz:** Es gibt keinen schnittfreien Beweis von $\vdash$ (leeres Succedent), also gilt $\text{LK} \not\vdash \bot$.

**Komplexität:** Schnittelimination kann zu nicht-elementarer Beweisverlängerung führen (Statman 1979): Die Größe kann $(n+1)$-fach exponentiell wachsen.

```python
subformula_property(proof_tree)  # True wenn Eigenschaft gilt
```

---

## Teil 3: Natürliches Schließen (ND)

### Grundidee

Im Natürlichen Schließen (Gentzen 1935) beginnt man mit **Annahmen** und leitet schrittweise die Konklusion ab. Annahmen können "entladen" werden.

### Klasse `NDProof`

```python
beweis = NDProof()
id_a   = beweis.assume("A")          # Annahme A (ID: 0)
id_b   = beweis.assume("B")          # Annahme B (ID: 1)
id_und = beweis.apply_and_intro(id_a, id_b)  # A∧B (ID: 2)
id_a2  = beweis.apply_and_elim_left(id_und)  # A   (ID: 3)
```

### Einführungsregeln

| Regel | Mathematisch | Methode |
|-------|-------------|---------|
| $\wedge$I | $\frac{A \quad B}{A \wedge B}$ | `apply_and_intro(id_a, id_b)` |
| $\vee$I$_L$ | $\frac{A}{A \vee B}$ | `apply_or_intro_left(id_a, "B")` |
| $\to$I | $\frac{[A] \quad B}{A \to B}$ (A entladen) | `apply_implies_intro(id_a, id_b)` |
| $\neg$I | $\frac{[A] \quad \bot}{\neg A}$ (A entladen) | `apply_neg_intro(id_a, id_bot)` |

### Eliminierungsregeln

| Regel | Mathematisch | Methode |
|-------|-------------|---------|
| $\wedge$E$_L$ | $\frac{A \wedge B}{A}$ | `apply_and_elim_left(id_und)` |
| $\wedge$E$_R$ | $\frac{A \wedge B}{B}$ | `apply_and_elim_right(id_und)` |
| $\to$E (MP) | $\frac{A \to B \quad A}{B}$ | `apply_implies_elim(id_impl, id_a)` |
| $\neg$E | $\frac{A \quad \neg A}{\bot}$ | `apply_neg_elim(id_a, id_neg)` |
| EFQ | $\frac{\bot}{A}$ | `apply_ex_falso(id_bot, "A")` |

### Automatischer Beweis

```python
beweis = prove_simple("A→A")        # Identität
beweis = prove_simple("(A∧B)→A")    # Projektion
beweis = prove_simple("A→(A∨B)")    # Disjunktions-Einführung
```

---

## Teil 4: Hilbert-Kalkül H

### Axiomenschemata

Das Hilbert-System verwendet drei Schemata und Modus Ponens:

$$K: \quad A \to (B \to A)$$
$$S: \quad (A \to (B \to C)) \to ((A \to B) \to (A \to C))$$
$$DN: \quad \neg\neg A \to A$$

Mit K und S allein erhält man das intuitionistische System. Mit DN zusätzlich: klassische Logik.

### Beispiel: Beweis von $A \to A$

```
1. (A → ((A→A) → A)) → ((A → (A→A)) → (A → A))   [Schema S]
2. A → ((A→A) → A)                                  [Schema K]
3. (A → (A→A)) → (A → A)                            [MP aus 1,2]
4. A → (A → A)                                       [Schema K]
5. A → A                                             [MP aus 3,4]
```

### Deduktionstheorem

$$\Gamma, A \vdash B \iff \Gamma \vdash A \to B$$

```python
ergebnis = deduction_theorem(["A", "B"], "C")
# ergebnis["neue_konklusion"] == "B → C"
```

---

## Teil 5: Beweiskomplexität

### Hierarchie der Beweissysteme

Beweissysteme nach aufsteigender Stärke (Simulierbarkeit):

$$\text{Resolution} \leq \text{Cutting Planes} \leq \text{Frege} \leq \text{Extended Frege}$$

### Cook-Reckhow-Theorem (1979)

$$P = NP \iff \text{Es gibt kein superpolynomiales Beweissystem}$$

Mit anderen Worten: Wenn P = NP, dann gibt es ein polynomiales System für alle Tautologien.

### PHP — Taubenschlag-Prinzip

Das Taubenschlag-Prinzip $\text{PHP}_n$: $n+1$ Tauben in $n$ Löcher → keine injektive Abbildung.

**Kodierung:**
- Variable $p_{i,j}$: "Taube $i$ ist in Loch $j$"
- Klauseln Typ 1: $p_{i,1} \vee \cdots \vee p_{i,n}$ (jede Taube muss irgendwo)
- Klauseln Typ 2: $\neg p_{i_1,j} \vee \neg p_{i_2,j}$ (kein Loch doppelt)

**Satz (Ben-Sasson & Wigderson 2001):** Jeder Resolution-Beweis von $\text{PHP}_n$ hat Größe $2^{\Omega(n)}$.

**Wichtiger Kontrast:** Frege-Beweise für PHP haben Größe $O(n^3)$ — exponentieller Unterschied!

```python
ergebnis = pigeonhole_principle_hardness()
# ergebnis["haerteste_klasse"]["gap"] — Erklärung des Unterschieds
```

---

## Teil 6: Ordinale Analyse

### Beweistheoretische Ordinalzahlen

Die **beweistheoretische Ordinalzahl** $|T|$ misst die Stärke eines formalen Systems:

| System | Ordinalzahl | Bedeutung |
|--------|------------|-----------|
| PRA | $\omega^\omega$ | Primitive rekursive Arithmetik |
| PA | $\varepsilon_0$ | Peano-Arithmetik |
| ACA₀ | $\varepsilon_0$ | Arithmetische Komprehension |
| ATR₀ | $\Gamma_0$ | Arithmetische transfinite Rekursion |
| $\Pi_1^1$-CA₀ | $\Psi(\Omega_\omega)$ | Π₁¹-Komprehension |

### ε₀ (Epsilon-Null)

$\varepsilon_0$ ist die kleinste Ordinalzahl $\alpha$ mit $\omega^\alpha = \alpha$:

$$\varepsilon_0 = \sup\{\omega,\ \omega^\omega,\ \omega^{\omega^\omega},\ \ldots\}$$

**Cantor-Normalform:** Jede Ordinalzahl $< \varepsilon_0$ hat eine eindeutige Darstellung:

$$\alpha = \omega^{\beta_1} \cdot c_1 + \omega^{\beta_2} \cdot c_2 + \cdots + \omega^{\beta_n} \cdot c_n, \quad \beta_1 > \cdots > \beta_n,\ c_i \in \mathbb{N}^+$$

```python
demo = epsilon_zero_demo()
# demo["cantor_normalform"]["beispiele"] — Beispiele für kleine Zahlen
# demo["gentzen"]["satz"]               — "Con(PA) ⟺ WO(ε₀)"
```

### Gentzens Konsistenzbeweis (1936)

**Satz:** $\text{Con}(\text{PA})$ ist beweisbar in $\text{PRA} + \text{TI}(\varepsilon_0)$.

**Hauptidee:**
1. Weise jedem PA-Beweis eine Ordinalzahl $o(\pi) < \varepsilon_0$ zu
2. Beweis-Transformationen reduzieren $o(\pi)$ streng
3. Da $\varepsilon_0$ wohlgeordnet ist, terminiert die Reduktion
4. Terminierung impliziert: kein Beweis von $\bot$

```python
skizze = gentzen_consistency_proof_sketch()
# skizze["beweisidee"] — Schrittweise Erläuterung
```

---

## Teil 7: Reverse Mathematics

### Die Big Five (Simpson 1999)

Reverse Mathematics untersucht, welche Axiome zur Beweisbarkeit nötig sind. Die fünf Hauptsysteme:

| System | Ordinalzahl | Charakteristische Sätze |
|--------|------------|------------------------|
| RCA₀ | $\omega^\omega$ | Rekursive Funktionen, Grundarithmetik |
| WKL₀ | $\omega^\omega$ | Heine-Borel, König-Lemma |
| ACA₀ | $\varepsilon_0$ | Bolzano-Weierstraß, Ramsey |
| ATR₀ | $\Gamma_0$ | Determinierheit, Wohlordnungsvergleich |
| $\Pi_1^1$-CA₀ | $\Psi(\Omega_\omega)$ | Π₁¹-Mengen, hyperarithmetische Analysis |

```
RCA₀ ⊂ WKL₀ ⊂ ACA₀ ⊂ ATR₀ ⊂ Π₁¹-CA₀
```

### Beispielanalysen

```python
# Bolzano-Weierstraß ist äquivalent zu ACA₀ (über RCA₀)
reverse_math_example("bolzano_weierstrass")

# Heine-Borel ist äquivalent zu WKL₀
reverse_math_example("heine_borel")

# Ramsey RT²₂ ist äquivalent zu ACA₀
reverse_math_example("ramsey")

# König-Lemma ist äquivalent zu WKL₀
reverse_math_example("konig")
```

---

## Verwendung

```python
from src.proof_theory_formal import (
    Sequent, ProofTree, LKRules,
    NDProof, prove_simple,
    HilbertCalculus, deduction_theorem,
    cut_elimination_demo, subformula_property,
    proof_complexity_comparison,
    resolution_proof_size, pigeonhole_principle_hardness,
    proof_theoretic_ordinal, epsilon_zero_demo,
    gentzen_consistency_proof_sketch,
    big_five_systems, reverse_math_example,
)

# Sequent erstellen und prüfen
s = Sequent(["A", "B"], ["A"])
print(s.is_axiom())   # True

# Beweisbaum aufbauen
blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
mit_extra = LKRules.weakening_left(blatt, "B")
print(mit_extra.is_valid())  # True

# Natürliches Schließen
beweis = NDProof()
a = beweis.assume("A")
b = beweis.assume("B")
und = beweis.apply_and_intro(a, b)
print(beweis.get_formula(und))   # (A∧B)

# Ordinalzahlen
info = proof_theoretic_ordinal("PA")
print(info["ordinal"])   # ε₀
```

---

## Literatur

- Gentzen, G. (1934): Untersuchungen über das logische Schließen I, II. *Math. Z.*
- Gentzen, G. (1936): Die Widerspruchsfreiheit der reinen Zahlentheorie.
- Cook, S. & Reckhow, R. (1979): The relative efficiency of propositional proof systems.
- Ben-Sasson, E. & Wigderson, A. (2001): Short proofs are narrow — resolution made simple.
- Simpson, S. (1999): Subsystems of Second Order Arithmetic. Springer.
- Buss, S. (1998): Handbook of Proof Theory. Elsevier.
