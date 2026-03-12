# sun_tzu_squares.py — Sun-Tzu-Quadrate und Chinesischer Restsatz

## Übersicht

Dieses Modul implementiert:

- **ChineseRemainderTheorem**: Klassischer CRT nach Sun-Tzu (~3. Jh. n.Chr.)
- **SunTzuSquares**: Aufeinanderfolgende Quadrate mit vorgegebenen Resten
- **QuadraticResidueAnalysis**: Legendre-Symbol, Euler-Kriterium, quadratisches Reziprozitätsgesetz

---

## Chinesischer Restsatz (CRT)

### Formulierung

**THEOREM** (Sun-Tzu, ~3. Jh. n.Chr.):

Seien $m_1, \ldots, m_n$ paarweise teilerfremde Moduli. Dann hat das System

$$x \equiv r_1 \pmod{m_1}, \quad x \equiv r_2 \pmod{m_2}, \quad \ldots, \quad x \equiv r_n \pmod{m_n}$$

genau eine Lösung modulo $M = m_1 \cdots m_n$.

### Konstruktiver Beweis

$$M = \prod_{i=1}^n m_i, \quad M_i = \frac{M}{m_i}, \quad y_i = M_i^{-1} \pmod{m_i}$$

$$x = \left(\sum_{i=1}^n r_i \cdot M_i \cdot y_i\right) \bmod M$$

### Beispiel (Sun-Tzu)

$$x \equiv 2 \pmod{3}, \quad x \equiv 3 \pmod{5}, \quad x \equiv 2 \pmod{7} \implies x \equiv 23 \pmod{105}$$

---

## Klasse `SunTzuSquares`

### Problemstellung

Gesucht: Ganzzahl $a \geq 0$ mit

$$a^2 \equiv r_0 \pmod{m_0}, \quad (a+1)^2 \equiv r_1 \pmod{m_1}, \quad \ldots, \quad (a+n-1)^2 \equiv r_{n-1} \pmod{m_{n-1}}$$

### Lösungsmethode (CRT-Reduktion)

1. Berechne alle Quadratwurzeln: $(a+i)^2 \equiv r_i \pmod{m_i}$
   → $a+i \equiv \pm\sqrt{r_i} \pmod{m_i}$ (via Tonelli-Shanks)
   → $a \equiv (\pm\sqrt{r_i} - i) \pmod{m_i}$

2. CRT löst das simultane System für $a$.

### Sun-Tzu-Vermutung #285

**CONJECTURE #285**:

Für paarweise teilerfremde Moduli $m_1, \ldots, m_n$ und vorgegebene quadratische Reste $r_i$ (mit $r_i = \text{QR mod } m_i$) existieren stets $n$ aufeinanderfolgende Quadrate $a^2, (a+1)^2, \ldots, (a+n-1)^2$ mit $(a+i)^2 \equiv r_i \pmod{m_i}$.

---

## Tonelli-Shanks-Algorithmus

Berechnet $\sqrt{a} \pmod{p}$ für ungerade Primzahl $p$:

1. **Euler-Kriterium**: $a^{(p-1)/2} \equiv 1 \pmod{p}$ ↔ $a$ ist QR
2. **Spezialfall** $p \equiv 3 \pmod{4}$: $\sqrt{a} = a^{(p+1)/4} \bmod p$
3. **Allgemein**: Schreibe $p-1 = Q \cdot 2^S$, iteriere bis $t \equiv 1$

---

## Quadratische Reste und Legendre-Symbol

### Definition

$$\left(\frac{a}{p}\right) = \begin{cases} 0 & p \mid a \\ 1 & a \text{ ist QR mod } p \\ -1 & a \text{ ist kein QR mod } p \end{cases}$$

### Euler-Kriterium

**THEOREM** (Euler):

$$\left(\frac{a}{p}\right) \equiv a^{(p-1)/2} \pmod{p}$$

### Quadratisches Reziprozitätsgesetz

**THEOREM** (Gauss, 1796):

$$\left(\frac{p}{q}\right) \cdot \left(\frac{q}{p}\right) = (-1)^{\frac{p-1}{2} \cdot \frac{q-1}{2}}$$

für verschiedene ungerade Primzahlen $p \neq q$.

**Folgerung**:
- Falls $p \equiv 1 \pmod{4}$ oder $q \equiv 1 \pmod{4}$: $\left(\frac{p}{q}\right) = \left(\frac{q}{p}\right)$
- Falls $p \equiv q \equiv 3 \pmod{4}$: $\left(\frac{p}{q}\right) = -\left(\frac{q}{p}\right)$

---

## Verwendungsbeispiele

```python
from sun_tzu_squares import SunTzuSquares, ChineseRemainderTheorem, QuadraticResidueAnalysis

# CRT
x, M = ChineseRemainderTheorem.solve([2, 3, 2], [3, 5, 7])
# x = 23, M = 105

# Aufeinanderfolgende Quadrate
st = SunTzuSquares([3, 5, 7], [1, 1, 1])
a = st.find_a()  # kleinstes a mit a²≡1(3), (a+1)²≡1(5), (a+2)²≡1(7)

# Quadratwurzeln
roots = st.sqrt_mod_prime(4, 5)  # = [2, 3]

# Legendre-Symbol
sym = st.legendre_symbol(2, 7)  # = 1 (2 ist QR mod 7)

# Reziprozitätsgesetz
result = QuadraticResidueAnalysis.quadratic_reciprocity(3, 7)
print(result["reciprocity_holds"])  # True
```

---

## Literatur

- Sun-Tzu: Suànjīng (~3. Jh. n.Chr.) — ursprüngliches CRT
- Tonelli (1891), Shanks (1972): Quadratwurzeln mod p
- Gauss (1796): Quadratisches Reziprozitätsgesetz
- Ireland, Rosen: "A Classical Introduction to Modern Number Theory"

**Autor:** Michael Fuhrmann
**Stand:** 2026-03-12
