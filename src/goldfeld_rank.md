# goldfeld_rank.py — Goldfeld-Rangvermutung und Selmer-Gruppen

## Übersicht

Dieses Modul implementiert die Analyse der Goldfeld-Rangvermutung für elliptische Kurven:

- **EllipticCurveRank**: Weierstraß-Kurven y² = x³ + Ax + B, Diskriminante, Frobenius-Spur
- **GoldfeldRankConjecture**: Durchschnittlicher Rang = 1/2, Bhargava-Shankar-Schranken
- **SelmerGroupAnalysis**: 2-Selmer-Gruppe, Torsionspunkte, Nagell-Lutz

---

## Klasse `EllipticCurveRank`

### Weierstraß-Form

Eine elliptische Kurve über ℚ in kurzer Weierstraß-Form:

$$E: y^2 = x^3 + Ax + B$$

**Diskriminante**:
$$\Delta = -16(4A^3 + 27B^2) \neq 0$$

Die Bedingung Δ ≠ 0 sichert, dass die Kurve nicht-singulär (glatt) ist.

### Frobenius-Spur

Für eine Primzahl p ohne schlechte Reduktion:

$$a_p = p + 1 - \#E(\mathbb{F}_p)$$

**THEOREM** (Hasse): $|a_p| \leq 2\sqrt{p}$

### L-Funktion

$$L(E, s) = \prod_{p \nmid \Delta} (1 - a_p p^{-s} + p^{1-2s})^{-1} \cdot \prod_{p \mid \Delta} (\cdots)$$

---

## Klasse `GoldfeldRankConjecture`

### Goldfeld-Vermutung

**CONJECTURE** (Goldfeld 1979):

$$\lim_{N \to \infty} \frac{1}{N^2} \sum_{\substack{|a|, |b| \leq N \\ \Delta \neq 0}} \operatorname{rg}(E_{a,b}(\mathbb{Q})) = \frac{1}{2}$$

> Der durchschnittliche Rang einer "zufälligen" elliptischen Kurve ist 1/2.

### Bhargava-Shankar-Schranken

**THEOREM** (Bhargava-Shankar 2013):

$$\text{avg. } \dim_{\mathbb{F}_2} S^{(2)}(E/\mathbb{Q}) = 3 \implies \text{avg. Rang} \leq \frac{3}{2}$$

**THEOREM** (Bhargava-Shankar 2015):

$$\text{avg. Rang} \leq \frac{7}{6} \approx 1{,}167$$

Dies ist die beste bekannte obere Schranke für den Durchschnittsrang.

### Heegner-Punkte und Klassenzahl h(-163) = 1

**THEOREM** (Baker-Heegner-Stark, 1966–1967):

Es gibt genau 9 imaginär-quadratische Körper ℚ(√-d) der Klassenzahl 1:
$$d \in \{1, 2, 3, 7, 11, 19, 43, 67, 163\}$$

**h(-163) = 1**: Der Körper ℚ(√-163) hat eindeutige Primfaktorzerlegung.

**THEOREM** (Goldfeld 1983): Für Heegner-Punkte gilt:
$$h(-D) \geq \frac{1}{7000} \log D$$

**THEOREM** (Gross-Zagier 1986):
Falls ord_{s=1} L(E, s) = 1, dann ist der Heegner-Punkt ein rationaler Punkt von unendlicher Ordnung, d.h. rang(E(ℚ)) = 1.

**THEOREM** (Kolyvagin 1988):
Falls ord_{s=1} L(E, s) ≤ 1, dann ist Ш(E/ℚ) endlich.

---

## Klasse `SelmerGroupAnalysis`

### 2-Selmer-Gruppe

Die 2-Selmer-Gruppe passt in die exakte Sequenz:

$$0 \to E(\mathbb{Q})/2E(\mathbb{Q}) \to S^{(2)}(E/\mathbb{Q}) \to \text{Ш}(E/\mathbb{Q})[2] \to 0$$

**Schlussfolgerung**: rang(E(ℚ)) ≤ dim S^{(2)} - dim Ш[2]

### Nagell-Lutz-Theorem

**THEOREM** (Nagell-Lutz):
Jeder ganzzahlige Torsionspunkt (x, y) mit y ≠ 0 erfüllt y² | Δ.

### Birch & Swinnerton-Dyer (BSD)

**CONJECTURE** (Birch & Swinnerton-Dyer, 1965):

$$\operatorname{rg}(E(\mathbb{Q})) = \operatorname{ord}_{s=1} L(E, s)$$

---

## Verwendungsbeispiele

```python
from goldfeld_rank import EllipticCurveRank, GoldfeldRankConjecture, SelmerGroupAnalysis

# Einzelne Kurve
curve = EllipticCurveRank(-1, 0)  # y² = x³ - x
print(curve.is_on_curve(0, 0))   # True
print(curve.ap(7))               # Frobenius-Spur bei p=7

# Goldfeld-Familie
goldfeld = GoldfeldRankConjecture(a_range=10, b_range=10)
goldfeld.build_curve_family()
print(goldfeld.average_rank_estimate())  # heuristischer Durchschnittsrang
print(goldfeld.heegner_point_example())  # h(-163)=1

# Selmer-Analyse
selmer = SelmerGroupAnalysis(curve)
print(selmer.torsion_subgroup_points())  # Torsionspunkte
print(selmer.selmer_rank_upper_bound())  # Selmer-Rang-Schranke
```

---

## Literatur

- Goldfeld (1979): Vermutung über Durchschnittsrang = 1/2
- Bhargava, Shankar (2013): avg. 2-Selmer ≤ 3/2
- Bhargava, Shankar (2015): avg. Rang ≤ 7/6
- Gross, Zagier (1986): Heegner-Punkt-Formel
- Kolyvagin (1988): Endlichkeit von Ш

**Autor:** Michael Fuhrmann
**Stand:** 2026-03-12
