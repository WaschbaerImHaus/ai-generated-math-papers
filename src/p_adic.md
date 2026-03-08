# p_adic.py – Dokumentation

**Datei:** `src/p_adic.py`
**Autor:** Kurt Ingwer
**Erstellt:** 2026-03-08
**Python:** 3.13+

---

## Übersicht

Das Modul `p_adic.py` implementiert die Theorie der **p-adischen Zahlen** – ein alternatives Zahlensystem, das auf einem völlig anderen Abstandsbegriff als die reellen Zahlen basiert.

Kernidee: Statt "kleine Unterschiede" durch kleine Absolutwerte zu messen, misst die p-adische Metrik, wie stark eine Zahl durch eine Primzahl p teilbar ist. Je stärker teilbar, desto "näher" bei Null.

---

## Mathematischer Hintergrund

### p-adische Bewertung

Für eine Primzahl p und eine ganze Zahl n ≠ 0:
```
v_p(n) = max{k ∈ ℕ : p^k | n}
```

Beispiele:
- v_2(12) = 2  (12 = 2² · 3)
- v_3(12) = 1  (12 = 2² · 3¹)
- v_5(12) = 0  (5 ∤ 12)
- v_p(0) = +∞  (0 ist durch jede Potenz teilbar)

**Bewertungsaxiome:**
- v_p(mn) = v_p(m) + v_p(n)
- v_p(m+n) ≥ min(v_p(m), v_p(n))  ← "ultrametrische" Eigenschaft

### p-adischer Betrag

```
|n|_p = p^{-v_p(n)}   (n ≠ 0)
|0|_p = 0
```

Der p-adische Betrag erfüllt die **ultrametrische Dreiecksungleichung**:
```
|m+n|_p ≤ max(|m|_p, |n|_p)
```
Diese ist stärker als die übliche Dreiecksungleichung!

### Satz von Ostrowski (1916)

Alle nicht-trivialen absoluten Werte auf ℚ sind äquivalent zu genau einem der:
- Dem gewöhnlichen Betrag |·|_∞
- Einem p-adischen Betrag |·|_p für eine Primzahl p

Die **Produktformel** verbindet alle Absolutwerte:
```
|n|_∞ · Π_{p prim} |n|_p = 1   für alle n ∈ ℤ, n ≠ 0
```

### p-adische Zahlen Q_p

Die p-adischen Zahlen Q_p sind die Vervollständigung von ℚ bezüglich |·|_p (analog wie ℝ die Vervollständigung bezüglich |·|_∞ ist).

Jede p-adische Zahl hat eine eindeutige Darstellung:
```
x = Σ_{k=v}^∞ a_k · p^k,   a_k ∈ {0, ..., p-1}
```

Besonderheit: Negative ganze Zahlen haben unendlich viele Einsen in der p-adischen Entwicklung! Z.B. in Q_2:
```
-1 = ...1111111_2  (alle Binärstellen sind 1)
```

---

## Funktionen und Klassen

### `p_adic_valuation(n, p) -> int | float`

Berechnet v_p(n). Gibt `float('inf')` für n=0 zurück.

```python
p_adic_valuation(12, 2)   # → 2
p_adic_valuation(12, 3)   # → 1
p_adic_valuation(0, 7)    # → inf
```

### `p_adic_norm(n, p) -> float`

Berechnet |n|_p = p^{-v_p(n)}.

```python
p_adic_norm(12, 2)   # → 0.25  (= 2^{-2})
p_adic_norm(12, 3)   # → 0.333 (= 3^{-1})
p_adic_norm(12, 5)   # → 1.0   (5 ∤ 12)
```

### `p_adic_distance(x, y, p) -> float`

d_p(x, y) = |x - y|_p. Ultrametrische Distanz.

```python
p_adic_distance(0, 8, 2)   # → 0.125 (= 2^{-3}, 8 ist nah bei 0!)
p_adic_distance(0, 3, 2)   # → 1.0   (3 ist weit von 0 in 2-adischer Metrik)
```

---

### Klasse `PAdicNumber`

Repräsentiert eine p-adische Zahl durch ihre Zifferentwicklung:
```
x = a_0 · p^0 + a_1 · p^1 + a_2 · p^2 + ...
```

#### `from_integer(n, p, precision=20)` (Klassenmethode)

Konvertiert eine ganze Zahl in p-adische Darstellung.

```python
x = PAdicNumber.from_integer(13, p=2, precision=8)
# 13 = 1·2^0 + 0·2^1 + 1·2^2 + 1·2^3 = 1101_2
# digits = [1, 0, 1, 1, 0, 0, 0, 0]
```

#### `from_fraction(num, den, p, precision=20)` (Klassenmethode)

Konvertiert einen Bruch num/den in p-adische Darstellung (via modularem Inversem).

#### `to_integer_approx(terms=10) -> int`

Konvertiert zurück zu einer ganzen Zahl modulo p^terms.

#### `norm() -> float`

Gibt die p-adische Norm zurück: p^{-v} wobei v die Position der ersten Nicht-Null-Ziffer ist.

#### Operatoren `+` und `*`

Addition und Multiplikation mit korrekter Übertragsbehandlung.

---

### `hensel_lift(f_coeffs, p, initial_root, n_lifts=5)`

Implementiert **Hensels Lemma** – das p-adische Newton-Verfahren.

**Grundaussage:** Sei f ein Polynom über ℤ. Wenn:
- f(a) ≡ 0 (mod p)  (Startwurzel)
- f'(a) ≢ 0 (mod p)  (Nicht-singulär)

Dann gibt es ein eindeutiges b ≡ a (mod p) mit f(b) ≡ 0 (mod p²).

**Liftformel:**
```
a_{k+1} = a_k - f(a_k) · [f'(a_k)]^{-1}  (mod p^{k+1})
```

**Beispiel:** x² ≡ -1 (mod 5^n)
```python
root = hensel_lift([1, 0, 1], p=5, initial_root=2, n_lifts=3)
# Startwurzel: 2 (da 2² = 4 ≡ -1 mod 5)
# Lift zu mod 625: root² ≡ -1 (mod 625)
```

---

### `p_adic_exp(x_digits, p, precision=10)`

p-adische Exponentialfunktion:
```
exp_p(x) = Σ_{n=0}^∞ x^n / n!
```

Konvergiert für |x|_p < p^{-1/(p-1)}, d.h.:
- Für p ≥ 3: v_p(x) ≥ 1
- Für p = 2: v_2(x) ≥ 2

### `p_adic_log(x_digits, p, precision=10)`

p-adischer Logarithmus:
```
log_p(1+y) = Σ_{n=1}^∞ (-1)^{n+1} y^n / n
```

Konvergiert für |x - 1|_p < 1 (d.h. x ≡ 1 mod p).

---

### `ostrowski_theorem_demo(n) -> dict`

Demonstriert den **Satz von Ostrowski** durch Berechnung der Produktformel:

```
|n|_∞ · Π_p |n|_p = 1
```

Rückgabe enthält:
- `abs_infinity`: gewöhnlicher Betrag |n|
- `p_adic_norms`: Dict {p: |n|_p} für kleine Primzahlen
- `prime_factorization`: Primfaktorzerlegung
- `product_formula_exact`: exaktes Produkt (sollte = 1 sein)
- `product_equals_one`: True falls Produktformel erfüllt

---

## Intuition und Anwendungen

### Warum p-adische Zahlen?

Das **Hasse-Minkowski-Prinzip** sagt (grob): Eine quadratische Form hat eine rationale Lösung genau dann, wenn sie reelle und p-adische Lösungen für alle p hat. D.h. lokale (p-adische) Aussagen implizieren globale Aussagen.

### Anwendungen in der modernen Mathematik

1. **Fermatsche Kurven:** p-adische Analysis ermöglicht Beweise über diophantische Gleichungen.
2. **p-adische L-Funktionen:** Iwasawa-Theorie, Kummer-Kongruenzen für Bernoulli-Zahlen.
3. **Kryptographie:** Gitterbasierte Kryptographie, p-adische Interpolation.
4. **Physik:** p-adische Strings, Quantengravitation.

### Seltsame Eigenschaft: "Jeder Punkt eines Balls ist Mittelpunkt"

In der ultrametrischen Topologie gilt: Für jeden Punkt y in einer offenen Kugel B(x, r) ist auch B(y, r) = B(x, r). Kugeln können nicht überlappen – sie sind entweder disjunkt oder eine enthält die andere.

---

## Numerische Hinweise

- **Endliche Approximation:** Alle Berechnungen sind modulo p^precision – echte p-adische Zahlen sind unendlich.
- **Negative Zahlen:** Im p-Komplement: -n = p^precision - n (mod p^precision).
- **Brüche:** Nur möglich, wenn der Nenner zu p koprim ist (sonst liegt ein p-adischer Pol vor).
- **Hensel-Lift:** Für singuläre Wurzeln (f'(a) ≡ 0 mod p) funktioniert das Lifting nicht.

---

## Literaturempfehlungen

- Neukirch, J.: *Algebraische Zahlentheorie*
- Gouvêa, F. Q.: *p-adic Numbers: An Introduction* (sehr zugänglich)
- Koblitz, N.: *p-adic Numbers, p-adic Analysis, and Zeta-Functions*
- Schikhof, W.: *Ultrametric Calculus*
