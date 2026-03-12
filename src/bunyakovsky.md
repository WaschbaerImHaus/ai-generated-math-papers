# bunyakovsky.py — Bunyakovsky-Vermutung

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Datum**: 2026-03-12

---

## Überblick

Die **Bunyakovsky-Vermutung** (Viktor Bunyakovsky, 1857) ist eine fundamentale
offene Frage über primzahlerzeugende Polynome:

> Sei $f(x) \in \mathbb{Z}[x]$ irreduzibel mit positivem Leitkoeffizienten.
> Falls $\gcd(\{f(n) : n \in \mathbb{N}\}) = 1$ (kein fester Primteiler),
> dann nimmt $f(n)$ **unendlich viele Primwerte** an.

---

## Bedingungen

### Bedingung 1: Irreduzibilität
$f(x)$ ist nicht als Produkt zweier nicht-konstanter ganzzahliger Polynome schreibbar.

### Bedingung 2: Kein fester Primteiler
$$\gcd(\{f(n) : n = 1, 2, 3, \ldots\}) = 1$$

**Gegenbeispiel**: $f(x) = x^2 - x = x(x-1)$ hat festen Primteiler 2,
da $f(n) = n(n-1)$ für alle $n$ gerade ist. → Vermutung **nicht anwendbar**.

### Bedingung 3: Positiver Leitkoeffizient
Damit $f(n) \to +\infty$ für $n \to \infty$ und Primwerte möglich sind.

---

## Bewiesene Spezialfälle

### Grad 1: Dirichlet-Theorem (1837, BEWIESEN)

Für $f(x) = ax + b$ mit $\gcd(a, b) = 1$:
> Es gibt unendlich viele Primzahlen der Form $an + b$ (Primzahlen in arithmetischen Progressionen).

**Beweis**: L-Funktionen $L(s, \chi)$ und Charakter-Orthogonalität.

---

## Offene Fälle (alle Conjectures)

| Polynom | Erste Primwerte | Status |
|---------|----------------|--------|
| $f(x) = x^2 + 1$ | $2, 5, 17, 37, 41, 53, 61, 73, 89, 97, \ldots$ | **Conjecture** |
| $f(x) = x^2 + x + 1$ | $3, 7, 13, 31, 43, 73, 157, \ldots$ | **Conjecture** |
| $f(x) = x^4 + 1$ | $2, 17, 82 = 2 \cdot 41, \ldots$ | **Conjecture** |
| $f(x) = x^2 + x + 41$ (Euler) | $41, 43, 47, 53, 61, \ldots$ | **Conjecture** |

---

## Hardy-Littlewood-Dichteformel

Für irreduzibles $f$ ohne feste Primteiler wird die Dichte heuristisch bestimmt:

$$\pi_f(x) \sim C_f \cdot \frac{x}{\ln x}$$

mit der **Hardy-Littlewood-Konstante**:

$$C_f = \prod_{p \text{ prim}} \frac{p - \rho_f(p)}{p - 1}$$

wobei $\rho_f(p) = \#\{a \in \{0, \ldots, p-1\} : f(a) \equiv 0 \pmod{p}\}$.

**Achtung**: Diese Formel ist für $\deg(f) \geq 2$ eine **Vermutung** (Conjecture).

---

## Bateman-Horn-Vermutung (1962, Conjecture)

Verallgemeinerung auf Systeme von Polynomen $f_1, \ldots, f_k$:

$$\#\{n \leq x : f_1(n), \ldots, f_k(n) \text{ alle prim}\} \sim C \cdot \frac{x}{(\ln x)^k}$$

mit einem Euler-Produkt-Korrekturfaktor $C$ über alle Primzahlen.

---

## Klassen

### `BunyakovskyConjecture`

| Methode | Beschreibung |
|---------|-------------|
| `fixed_prime_divisor(limit)` | Findet ersten festen Primteiler (falls vorhanden) |
| `gcd_of_values(limit)` | Berechnet $\gcd(\{f(n)\})$ |
| `find_prime_values(n_limit)` | Alle $n$ mit $f(n)$ prim |
| `prime_density(x)` | $\pi_f(x) / x$ |
| `hardy_littlewood_constant(limit)` | Numerische C_f-Berechnung |
| `bunyakovsky_criterion_check()` | Prüft alle Bedingungen |
| `compare_with_heuristic(x)` | Tatsächlich vs. H-L-Vorhersage |

### Hilfsfunktionen

| Funktion | Erstellt |
|---------|---------|
| `make_x_squared_plus_1()` | $f(x) = x^2 + 1$ |
| `make_x_squared_plus_x_plus_1()` | $f(x) = x^2 + x + 1$ |
| `make_linear(a, b)` | $f(x) = ax + b$ (Dirichlet, BEWIESEN) |

---

## Beispiel: f(x) = x² + 1

```python
f = make_x_squared_plus_1()
result = f.bunyakovsky_criterion_check()
# gcd = 1, keine festen Primteiler, Conjecture anwendbar
prime_vals = f.find_prime_values(1000)
# ~183 Primwerte von 1000 Kandidaten (Dichte ~18.3%)
```

Theoretische Dichte nach H-L: $C_f / \ln(x) \approx 1.33 / \ln(1000) \approx 19.4\%$ — sehr gut!
