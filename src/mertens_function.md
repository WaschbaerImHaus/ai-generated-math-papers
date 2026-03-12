# mertens_function.py — Mertens-Funktion und Liouville-Funktion

## Übersicht

Dieses Modul implementiert zwei fundamentale additive Funktionen der analytischen Zahlentheorie:

- **Möbius-Funktion** μ(n) via linearem Sieb
- **Mertens-Funktion** M(x) = Σ_{n≤x} μ(n)
- **Liouville-Funktion** λ(n) = (−1)^{Ω(n)}
- **Summatorische Liouville-Funktion** L(x) = Σ_{n≤x} λ(n)

---

## Klasse `MertensFunction`

### Mathematischer Hintergrund

Die **Möbius-Funktion** μ: ℕ → {−1, 0, 1} ist definiert als:

$$\mu(n) = \begin{cases} 1 & n = 1 \\ (-1)^k & n = p_1 p_2 \cdots p_k \text{ (quadratfrei, } k \text{ Primfaktoren)} \\ 0 & p^2 \mid n \text{ für ein Prim } p \end{cases}$$

Die **Mertens-Funktion** akkumuliert μ:

$$M(x) = \sum_{n=1}^{\lfloor x \rfloor} \mu(n)$$

### Bekannte Werte

| x | M(x) |
|---|------|
| 1 | 1 |
| 2 | 0 |
| 3 | -1 |
| 10 | -1 |
| 100 | 1 |
| 1000 | 2 |

### Theoretische Verbindungen

#### Verbindung zur Riemann-Hypothese

**THEOREM** (Unter RH): M(x) = O(x^{1/2 + ε}) für alle ε > 0.

Dies ist sogar äquivalent zur RH:

$$\text{RH} \iff M(x) = O(x^{1/2+\varepsilon}) \; \forall \varepsilon > 0$$

#### Explizite Formel

Formal (unter RH):

$$M(x) \sim \sum_{\substack{\rho \\ \zeta(\rho)=0}} \frac{x^\rho}{\rho \cdot \zeta'(\rho)}$$

Das Wachstumsverhalten ist durch die Realteile der nicht-trivialen Nullstellen bestimmt.

#### Mertens-Vermutung (WIDERLEGT)

**CONJECTURE** (Mertens 1897): $|M(x)| < \sqrt{x}$ für alle $x \geq 1$.

**THEOREM** (Odlyzko/te Riele 1985): Die Mertens-Vermutung ist **falsch**.

> Bis heute ist kein explizites Gegenbeispiel bekannt. Heuristisch liegt das kleinste Gegenbeispiel bei $x \approx \exp(1{,}59 \cdot 10^{40})$.

#### Pintz-Resultat

**THEOREM** (Pintz 1987):

$$\limsup_{x \to \infty} \frac{M(x)}{\sqrt{x}} > 1{,}06$$

### Implementierung

Der **lineare Sieb** berechnet μ(n) in O(N log log N) Zeit:

```python
mf = MertensFunction(10000)
mf.M(1000)   # = 2
mf.mobius(6) # = 1  (6 = 2·3, zwei Primfaktoren)
```

---

## Klasse `LiouvilleFunction`

### Mathematischer Hintergrund

**Ω(n)** = totale Anzahl Primfaktoren mit Vielfachheit:
- Ω(1) = 0, Ω(p) = 1, Ω(p²) = 2, Ω(12) = 3 (12 = 2²·3)

**Liouville-Funktion**:
$$\lambda(n) = (-1)^{\Omega(n)}$$

**Summatorische Funktion**:
$$L(x) = \sum_{n=1}^{\lfloor x \rfloor} \lambda(n)$$

### Eigenschaften

- λ ist **vollständig multiplikativ**: λ(mn) = λ(m)·λ(n) für alle m, n
- μ ist nur multiplikativ: μ(mn) = μ(m)μ(n) für gcd(m,n) = 1

**Faltungsbeziehung**:
$$\sum_{d^2 \mid n} \lambda(n/d^2) = \mu(n)$$

### Pólya-Vermutung (WIDERLEGT)

**CONJECTURE** (Pólya 1919): L(x) ≤ 0 für alle x ≥ 2.

**THEOREM** (Haselgrove 1960): Die Pólya-Vermutung ist **falsch**.

Kleinstes bekanntes Gegenbeispiel: x ≈ **906 316 571**.

---

## Verwendungsbeispiele

```python
# Mertens-Funktion
mf = MertensFunction(10000)
print(mf.M(10))     # -1
print(mf.M(100))    # 1
print(mf.M(1000))   # 2
print(mf.mobius(30)) # -1 (30 = 2·3·5)

# Ratio-Analyse
_, max_r = mf.max_ratio()
print(f"Max |M(x)/√x| = {max_r:.4f}")  # sollte < 1 bis 10000

# Liouville-Funktion
lf = LiouvilleFunction(1000)
print(lf.liouville(4))  # 1  (4=2², Ω=2)
print(lf.liouville(6))  # 1  (6=2·3, Ω=2)
print(lf.L(100))        # ...
```

---

## Literatur

- Mertens (1897): Vermutung über M(x)
- Odlyzko, te Riele (1985): Widerlegung der Mertens-Vermutung
- Pintz (1987): Ω-Resultat für M(x)/√x > 1.06
- Haselgrove (1960): Widerlegung der Pólya-Vermutung

**Autor:** Michael Fuhrmann
**Stand:** 2026-03-12
