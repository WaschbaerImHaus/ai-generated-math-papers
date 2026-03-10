# representation_theory.py — Darstellungstheorie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Modul:** `src/representation_theory.py`

---

## Überblick

Das Modul implementiert die Grundlagen der **Darstellungstheorie endlicher Gruppen**.
Eine Darstellung verbindet abstrakte Gruppenelemente mit konkreten invertierbaren Matrizen
und ermöglicht so algebraische Probleme durch lineare Algebra zu lösen.

---

## Mathematische Grundlagen

### Darstellung einer Gruppe

Eine **Darstellung** einer Gruppe $G$ über $\mathbb{C}$ ist ein Gruppenhomomorphismus

$$\rho: G \to \mathrm{GL}(n, \mathbb{C})$$

d.h. eine Abbildung, die jedem Gruppenelement $g \in G$ eine invertierbare $n \times n$-Matrix $\rho(g)$ zuordnet, sodass:

$$\rho(g \cdot h) = \rho(g) \cdot \rho(h) \quad \forall g, h \in G$$

Die Zahl $n$ heißt **Dimension** (oder Grad) der Darstellung.

### Charakter

Der **Charakter** $\chi_\rho: G \to \mathbb{C}$ einer Darstellung $\rho$ ist die Spurabbildung:

$$\chi_\rho(g) = \mathrm{Tr}(\rho(g))$$

Wichtige Eigenschaften:
- $\chi_\rho$ ist eine **Klassenfunktion**: $\chi_\rho(g) = \chi_\rho(hgh^{-1})$
- Äquivalente Darstellungen haben **gleiche Charaktere**
- $\chi_\rho(e) = \dim(\rho)$

### Irreduzibilität und Schur-Orthogonalität

Eine Darstellung $\rho$ heißt **irreduzibel**, wenn kein nichttrivialer invarianter Teilraum existiert.

Das **Schur-Orthogonalitätskriterium** besagt:

$$\langle \chi_i, \chi_j \rangle = \frac{1}{|G|} \sum_{g \in G} \chi_i(g) \overline{\chi_j(g)} = \delta_{ij}$$

Daraus folgt: $\rho$ ist irreduzibel genau dann wenn $\langle \chi_\rho, \chi_\rho \rangle = 1$.

---

## Schur-Lemma

**Satz (Schur):** Sei $T: V \to W$ ein Intertwiner (Morphismus), d.h. $T \circ \rho(g) = \sigma(g) \circ T$ für alle $g \in G$. Dann gilt:

1. Falls $\rho$ und $\sigma$ **inäquivalent** sind: $T = 0$
2. Falls $\rho = \sigma$ **irreduzibel** über $\mathbb{C}$: $T = \lambda I$ für ein $\lambda \in \mathbb{C}$

**Konsequenz:** Irreduzible Darstellungen sind "unzerlegbar" – kein Intertwiner kann sie vermischen.

---

## Maschke-Satz

**Satz (Maschke):** Sei $G$ eine endliche Gruppe, $k$ ein Körper mit $\mathrm{char}(k) \nmid |G|$ (insbesondere $\mathrm{char}(k) = 0$). Dann ist jede Darstellung von $G$ über $k$ **vollständig reduzibel**:

$$\rho \cong m_1 \rho_1 \oplus m_2 \rho_2 \oplus \cdots \oplus m_r \rho_r$$

wobei $\rho_1, \ldots, \rho_r$ irreduzibel sind und $m_i \geq 0$ die Multiplizitäten.

**Fails:** Für $\mathrm{char}(k) \mid |G|$ (modulare Darstellungstheorie) gilt der Satz nicht.

---

## Zerlegung in Irreduzible

Die **Multiplizität** $m_i$ der $i$-ten irreduziblen Darstellung in $\rho$ berechnet sich durch:

$$m_i = \langle \chi_\rho, \chi_i \rangle = \frac{1}{|G|} \sum_{g \in G} \chi_\rho(g) \overline{\chi_i(g)}$$

---

## Burnside-Lemma (Cauchy-Frobenius)

Wirkt eine Gruppe $G$ auf einer Menge $X$, dann ist die **Anzahl der Orbits**:

$$|X/G| = \frac{1}{|G|} \sum_{g \in G} |\mathrm{Fix}(g)|$$

wobei $\mathrm{Fix}(g) = \{x \in X \mid g \cdot x = x\}$.

---

## Charaktertafel

Die **Charaktertafel** einer Gruppe $G$ ist die Matrix $(\chi_i(C_j))$ mit:
- Zeilen: irreduzible Darstellungen $\rho_1, \ldots, \rho_r$
- Spalten: Konjugationsklassen $C_1, \ldots, C_r$

**Wichtige Relation:** $\sum_i \dim(\rho_i)^2 = |G|$

---

## Beispiel: $\mathbb{Z}/2\mathbb{Z}$

| $g$ | $\rho_1(g)$ | $\rho_2(g)$ |
|-----|------------|------------|
| $0$ | $[1]$ | $[1]$ |
| $1$ | $[1]$ | $[-1]$ |

Charaktertafel:

| | $\{0\}$ | $\{1\}$ |
|--|--------|--------|
| $\chi_1$ | $1$ | $1$ |
| $\chi_2$ | $1$ | $-1$ |

---

## Beispiel: $S_3$

$S_3$ hat $|S_3| = 6$ und 3 Konjugationsklassen:
- $C_1 = \{e\}$ (1 Element)
- $C_2 = \{(12), (13), (23)\}$ (3 Transpositionen)
- $C_3 = \{(123), (132)\}$ (2 Dreizyklen)

**Charaktertafel:**

| Darst. | $\dim$ | $C_1$ | $C_2$ | $C_3$ |
|--------|--------|-------|-------|-------|
| trivial | $1$ | $1$ | $1$ | $1$ |
| Signum | $1$ | $1$ | $-1$ | $1$ |
| Standard | $2$ | $2$ | $0$ | $-1$ |

Probe: $1^2 + 1^2 + 2^2 = 6 = |S_3|$ ✓

---

## API-Referenz

### Klasse `Representation`

| Methode | Beschreibung |
|---------|-------------|
| `dimension()` | Gibt $n$ zurück |
| `is_homomorphism()` | Prüft $\rho(gh) = \rho(g)\rho(h)$ |
| `character()` | Berechnet $\chi(g) = \mathrm{Tr}(\rho(g))$ |
| `is_irreducible()` | Prüft $\langle\chi,\chi\rangle = 1$ |
| `direct_sum(other)` | $\rho \oplus \sigma$ (block-diagonal) |
| `tensor_product(other)` | $\rho \otimes \sigma$ (Kronecker) |
| `is_equivalent(other)` | Charaktervergleich |

### Freie Funktionen

| Funktion | Beschreibung |
|----------|-------------|
| `trivial_representation(elems)` | $\rho(g) = [1]$ |
| `regular_representation(elems, table)` | Links-Multiplikation |
| `character_table(elems, classes, reps)` | Vollständige Charaktertafel |
| `schur_lemma(rho, sigma, T)` | Intertwiner-Analyse |
| `maschke_theorem_check(order, char)` | Maschke-Anwendbarkeit |
| `decompose_representation(rep, irreps)` | Multiplizitäten $m_i$ |
| `burnside_lemma(elems, action)` | Anzahl der Orbits |
| `z2_representations()` | Irreduzible von $\mathbb{Z}/2\mathbb{Z}$ |
| `s3_representations()` | Irreduzible von $S_3$ |
