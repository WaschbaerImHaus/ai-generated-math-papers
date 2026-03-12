# agoh_giuga.py — Agoh-Vermutung, Giuga-Vermutung und Bernoulli-Zahlen

## Übersicht

Dieses Modul implementiert zwei eng verwandte Primzahl-Charakterisierungen:

- **BernoulliNumbers**: Exakte Bernoulli-Zahlen B_0..B_N via sympy
- **AgohConjecture**: n prim ↔ n·B_{n-1} ≡ −1 (mod n)
- **GiugaConjecture**: n prim ↔ Giuga-Bedingung für alle Primteiler
- **AgohGiuga**: Kombination beider Vermutungen + Äquivalenzanalyse

---

## Bernoulli-Zahlen

### Definition

Die Bernoulli-Zahlen sind definiert durch die erzeugende Funktion:

$$\frac{t}{e^t - 1} = \sum_{n=0}^{\infty} B_n \frac{t^n}{n!}$$

### Bekannte Werte (sympy-Konvention: B₁ = +1/2)

| k | B_k |
|---|-----|
| 0 | 1 |
| 1 | +1/2 |
| 2 | 1/6 |
| 3 | 0 |
| 4 | −1/30 |
| 5 | 0 |
| 6 | 1/42 |
| 8 | −1/30 |
| 10 | 5/66 |
| 12 | −691/2730 |
| 20 | −174611/330 |

**Hinweis**: Es gibt zwei Konventionen:
- Erste Konvention: B₁ = −1/2 (historisch, Agoh-Formulierung)
- Zweite Konvention: B₁ = +1/2 (sympy, modernere Quellen)

### Eigenschaften

1. **Verschwinden**: $B_{2k+1} = 0$ für alle $k \geq 1$
2. **Vorzeichen**: $(-1)^{k+1} B_{2k} > 0$ für $k \geq 1$

### Von-Staudt-Clausen-Theorem

**THEOREM** (von Staudt, Clausen 1840):

$$B_{2k} + \sum_{\substack{p \text{ prim} \\ (p-1) \mid 2k}} \frac{1}{p} \in \mathbb{Z}$$

Beispiel: $B_2 + \frac{1}{2} + \frac{1}{3} = \frac{1}{6} + \frac{1}{2} + \frac{1}{3} = 1 \in \mathbb{Z}$

---

## Agoh-Vermutung

### Formulierung

**CONJECTURE** (Agoh 1995):

$$n \text{ ist prim} \iff n \cdot B_{n-1} \equiv -1 \pmod{n}$$

### Interpretation

Sei $B_{n-1} = p/q$ (in niedrigsten Termen). Die Kongruenz bedeutet:

$$n \cdot \frac{p}{q} + 1 \equiv 0 \pmod{n} \iff n \mid (np + q)$$

Für $\gcd(q, n) = 1$ vereinfacht sich dies zu: der Zähler von $n \cdot B_{n-1} + 1$ ist durch $n$ teilbar.

### Numerische Verifikation

Die Vermutung ist verifiziert für alle $n < 10^{13800}$.

**Verbindung zur Wilson-Vermutung**: Beide charakterisieren Primzahlen über modular-algebraische Bedingungen.

---

## Giuga-Vermutung

### Formulierung

**CONJECTURE** (Giuga 1950):

$n$ ist prim genau dann, wenn für jeden Primteiler $p$ von $n$ gilt:

$$p \mid \left(\frac{n}{p} - 1\right) \quad \text{und} \quad (p-1) \mid \left(\frac{n}{p} - 1\right)$$

### Äquivalente Summenformulierung

**CONJECTURE** (Borwein et al. 1996):

$$n \text{ ist prim} \iff \sum_{k=1}^{n-1} k^{n-1} \equiv -1 \pmod{n}$$

### Giuga-Zahlen

Eine **Giuga-Zahl** wäre ein zusammengesetztes $n$, das die Giuga-Bedingung trotzdem erfüllt.

**CONJECTURE**: Es existieren keine Giuga-Zahlen (verifiziert bis $n = 10^{13800}$).

---

## Äquivalenz Agoh ↔ Giuga

**THEOREM** (Borwein, Borwein, Borwein, Girgensohn 1996):

$$\text{Agoh-Vermutung} \iff \text{Giuga-Vermutung}$$

Beide Vermutungen sind äquivalent, aber beide **offen** (unbewiesen).

---

## Vergleich mit anderen Primzahl-Charakterisierungen

| Satz | Formulierung | Status |
|------|--------------|--------|
| Wilson | $n$ prim $\iff$ $(n-1)! \equiv -1 \pmod{n}$ | **THEOREM** |
| Fermat (Klein) | $n$ prim $\Rightarrow$ $a^{n-1} \equiv 1 \pmod{n}$ für $\gcd(a,n)=1$ | **THEOREM** |
| Giuga | $n$ prim $\iff$ Giuga-Bedingung | **CONJECTURE** |
| Agoh | $n$ prim $\iff$ $n \cdot B_{n-1} \equiv -1 \pmod{n}$ | **CONJECTURE** |

---

## Verwendungsbeispiele

```python
from agoh_giuga import BernoulliNumbers, AgohConjecture, GiugaConjecture, AgohGiuga

# Bernoulli-Zahlen
bern = BernoulliNumbers(max_n=20)
print(bern.B(4))   # Fraction(-1, 30)
print(bern.B(10))  # Fraction(5, 66)

# Agoh-Bedingung
agoh = AgohConjecture()
print(agoh.check_agoh(7))   # True (7 ist prim)
print(agoh.check_agoh(9))   # False (9 = 3² zusammengesetzt)

# Giuga-Bedingung
giuga = GiugaConjecture()
print(giuga.is_giuga_prime(13))  # True
result = giuga.check_range(100)  # keine Giuga-Zahlen bis 100

# Kombinierte Analyse
ag = AgohGiuga()
result = ag.equivalence_check(30)
print(result["equivalence_holds_in_range"])  # True
```

---

## Hinweise zur Bernoulli-Konvention

sympy nutzt die **zweite Konvention** B₁ = +1/2. In der historischen Literatur
(insbesondere Agoh 1995) wird oft B₁ = −1/2 verwendet.

Die Agoh-Vermutung ist so formulierbar, dass sie von der Konvention unabhängig ist —
die entscheidende Eigenschaft ist die Primzahl-Charakterisierung über Kongruenzen.

---

## Literatur

- Bernoulli (1713): Ars Conjectandi — Einführung der Bernoulli-Zahlen
- Giuga (1950): Über eine Vermutung von Giuga
- Agoh (1995): On Giuga's Conjecture
- Borwein et al. (1996): Giuga's Conjecture on Primality
- von Staudt (1840), Clausen (1840): Staudt-Clausen-Theorem

**Autor:** Michael Fuhrmann
**Stand:** 2026-03-12
