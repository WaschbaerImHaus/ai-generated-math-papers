# arbitrary_precision.py – Dokumentation

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Sprache:** Python 3.13 / mpmath 1.3.0

---

## Überblick

Das Modul `arbitrary_precision.py` stellt Funktionen für **beliebig genaue
numerische Berechnungen** bereit. Als Backend wird `mpmath` (Multiple Precision
Mathematics) verwendet, das Berechnungen mit einer wählbaren Anzahl von
Dezimalstellen ermöglicht – weit jenseits der float64-Genauigkeit (≈15 Stellen).

**Typische Anwendungen:**
- Verifikation der Riemann-Hypothese auf hohe Stellenzahl
- Berechnung mathematischer Konstanten (π, e, γ) auf Tausende von Stellen
- Exakte Bernoulli-Zahlen für zahlentheoretische Formeln
- Kettenbruchentwicklungen für Irrationalitätsnachweise

---

## Präzision

```python
import mpmath
mpmath.mp.dps = 50    # 50 Dezimalstellen
```

Die globale Präzision `mp.dps` (decimal places) steuert die Genauigkeit aller
mpmath-Operationen. Mehr Stellen bedeuten langsamere, aber genauere Berechnungen.

| `mp.dps` | Anwendung |
|----------|-----------|
| 15 | float64-äquivalent |
| 50 | Riemann-Nullstellen (Standardwert) |
| 100 | Hohe numerische Präzision |
| 1000+ | π/e-Berechnungen zu Demonstrationszwecken |

---

## Funktionen

### `set_precision(decimal_digits)`

Setzt die globale mpmath-Präzision.

```python
set_precision(100)
```

---

### `riemann_zeta_highprec(s, digits=50)`

Berechnet die **Riemann-Zeta-Funktion** $\zeta(s)$:

$$\zeta(s) = \sum_{n=1}^{\infty} \frac{1}{n^s} \quad \text{für } \text{Re}(s) > 1$$

Durch analytische Fortsetzung für alle $s \in \mathbb{C} \setminus \{1\}$ definiert.

**Spezielle Werte:**

$$\zeta(2) = \frac{\pi^2}{6} \approx 1.6449\ldots$$

$$\zeta(4) = \frac{\pi^4}{90} \approx 1.0823\ldots$$

$$\zeta(-1) = -\frac{1}{12} \quad \text{(analytische Fortsetzung)}$$

$$\zeta(0) = -\frac{1}{2}$$

**Funktionalgleichung:**

$$\zeta(s) = 2^s \pi^{s-1} \sin\!\left(\frac{\pi s}{2}\right) \Gamma(1-s)\, \zeta(1-s)$$

```python
z = riemann_zeta_highprec(2, digits=50)
# → (1.6449340668482264364724...+0j)
```

---

### `riemann_zero_verify(n, digits=100)`

Verifiziert die $n$-te nicht-triviale **Riemann-Nullstelle** $\rho_n$.

Die **Riemann-Hypothese** (eines der Millennium-Probleme) besagt:

$$\text{Re}(\rho_n) = \frac{1}{2} \quad \text{für alle nicht-trivialen Nullstellen } \rho_n$$

Bekannte Nullstellen:

| $n$ | $\rho_n = \frac{1}{2} + it_n$ |
|-----|-------------------------------|
| 1 | $\frac{1}{2} + 14.134725\ldots i$ |
| 2 | $\frac{1}{2} + 21.022040\ldots i$ |
| 3 | $\frac{1}{2} + 25.010858\ldots i$ |
| 4 | $\frac{1}{2} + 30.424876\ldots i$ |

```python
result = riemann_zero_verify(1, digits=100)
# {'n': 1, 'real_part': 0.5, 'imag_part': '14.134...', 'on_critical_line': True}
```

---

### `pi_highprec(digits=1000)`

Berechnet $\pi$ auf die gewünschte Stellenzahl.

mpmath nutzt den **Chudnovsky-Algorithmus** (ca. 14 Dezimalstellen pro Term):

$$\frac{1}{\pi} = \frac{12}{640320^{3/2}} \sum_{k=0}^{\infty} \frac{(-1)^k (6k)! (13591409 + 545140134k)}{(3k)! (k!)^3 640320^{3k}}$$

```python
pi_str = pi_highprec(digits=100)
# "3.14159265358979323846264338327950288419716939937510..."
```

---

### `e_highprec(digits=1000)`

Berechnet die **Euler-Zahl** $e$ auf die gewünschte Stellenzahl.

$$e = \sum_{n=0}^{\infty} \frac{1}{n!} = 1 + 1 + \frac{1}{2} + \frac{1}{6} + \frac{1}{24} + \cdots$$

```python
e_str = e_highprec(digits=100)
# "2.71828182845904523536028747135266249775724709369995..."
```

---

### `gamma_highprec(s, digits=50)`

Berechnet die **Gamma-Funktion** $\Gamma(s)$:

$$\Gamma(s) = \int_0^{\infty} t^{s-1} e^{-t}\, dt \quad \text{für } \text{Re}(s) > 0$$

**Rekursionsformel:**

$$\Gamma(s+1) = s \cdot \Gamma(s)$$

**Spezielle Werte:**

$$\Gamma(n) = (n-1)! \quad \text{für } n \in \mathbb{N}$$

$$\Gamma\!\left(\tfrac{1}{2}\right) = \sqrt{\pi}$$

$$\Gamma'(1) = -\gamma \quad \text{(Euler-Mascheroni-Konstante)}$$

```python
gamma_val = gamma_highprec(5, digits=50)  # = 24.0
gamma_half = gamma_highprec(0.5, digits=50)  # = sqrt(π)
```

---

### `bernoulli_highprec(n, digits=50)`

Berechnet die **Bernoulli-Zahl** $B_n$:

Erzeugende Funktion:

$$\frac{t}{e^t - 1} = \sum_{n=0}^{\infty} B_n \frac{t^n}{n!}$$

**Erste Werte:**

| $n$ | $B_n$ |
|-----|-------|
| 0 | $1$ |
| 1 | $-\frac{1}{2}$ |
| 2 | $\frac{1}{6}$ |
| 3 | $0$ |
| 4 | $-\frac{1}{30}$ |
| 6 | $\frac{1}{42}$ |

Für ungerade $n > 1$ gilt: $B_n = 0$.

**Zusammenhang mit Zeta-Funktion:**

$$\zeta(-n) = -\frac{B_{n+1}}{n+1} \quad \text{für } n \geq 0$$

```python
b2 = float(bernoulli_highprec(2, digits=30))  # = 1/6
b4 = float(bernoulli_highprec(4, digits=30))  # = -1/30
```

---

### `continued_fraction_expansion(x, n_terms=20)`

Berechnet die **Kettenbruchentwicklung** einer reellen Zahl:

$$x = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cfrac{1}{a_3 + \cdots}}}$$

Kurzschreibweise: $x = [a_0;\, a_1,\, a_2,\, a_3,\, \ldots]$

**Algorithmus:**
1. $a_i = \lfloor x_i \rfloor$
2. $x_{i+1} = \frac{1}{x_i - a_i}$

**Bekannte Kettenbrüche:**

| Zahl | Kettenbruch |
|------|-------------|
| $\pi$ | $[3; 7, 15, 1, 292, 1, 1, 1, 2, \ldots]$ |
| $e$ | $[2; 1, 2, 1, 1, 4, 1, 1, 6, \ldots]$ |
| $\sqrt{2}$ | $[1; 2, 2, 2, 2, \ldots]$ (periodisch) |
| $\varphi$ | $[1; 1, 1, 1, 1, \ldots]$ (goldener Schnitt) |

```python
cf_pi = continued_fraction_expansion(math.pi, n_terms=5)
# [3, 7, 15, 1, 292]
```

---

### Weitere Hilfsfunktionen

#### `euler_mascheroni_highprec(digits=50)`

Die **Euler-Mascheroni-Konstante**:

$$\gamma = \lim_{n \to \infty} \left(\sum_{k=1}^n \frac{1}{k} - \ln n\right) \approx 0.5772156649\ldots$$

#### `log_highprec(x, base=None, digits=50)`

Logarithmus mit hoher Präzision:

$$\ln(x) = \int_1^x \frac{dt}{t} \quad \text{oder} \quad \log_b(x) = \frac{\ln x}{\ln b}$$

#### `sqrt_highprec(x, digits=50)`

Quadratwurzel mit hoher Präzision via Newton-Iteration.

---

## Verwendungsbeispiel

```python
from arbitrary_precision import (
    riemann_zeta_highprec, riemann_zero_verify, pi_highprec,
    bernoulli_highprec, continued_fraction_expansion
)
import math

# π auf 50 Stellen
pi_50 = pi_highprec(50)
print(pi_50)

# Erste Riemann-Nullstelle verifizieren
zero1 = riemann_zero_verify(1, digits=50)
print(f"ρ₁ ≈ {zero1['real_part']} + {zero1['imag_part'][:10]}i")
print(f"Auf kritischer Linie: {zero1['on_critical_line']}")

# ζ(2) = π²/6
z2 = riemann_zeta_highprec(2, digits=50)
print(f"ζ(2) = {z2.real:.20f}")
print(f"π²/6 = {math.pi**2/6:.20f}")

# Kettenbruch von π
cf = continued_fraction_expansion(math.pi, n_terms=8)
print(f"π = [{cf[0]}; {cf[1:]}")
```

---

## Hinweise zur Präzision

- Mehr Dezimalstellen → deutlich langsamere Berechnung
- `riemann_zero_verify` mit `digits=100` dauert einige Sekunden
- Für `pi_highprec(digits=10000)` sind wenige Sekunden Rechenzeit normal
- Die Genauigkeit von `riemann_zero_verify` ist durch mpmath's interne
  `zetazero`-Funktion begrenzt (sehr hochwertige Implementierung)

---

## Literatur

- Borwein, Bailey, Girgensohn: *Experimentation in Mathematics* (2004)
- Titchmarsh: *The Theory of the Riemann Zeta-Function* (1986)
- mpmath-Dokumentation: https://mpmath.org/doc/current/
- Khinchin: *Continued Fractions* (1997)
