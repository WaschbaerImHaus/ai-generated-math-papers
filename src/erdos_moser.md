# erdos_moser.py – Erdős-Moser-Vermutung

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Erstellt**: 2026-03-12
**Status**: CONJECTURE – offen

---

## Übersicht

Das Modul `erdos_moser.py` untersucht die **Erdős-Moser-Vermutung** aus der Zahlentheorie:

> **CONJECTURE (Erdős-Moser, 1953)**
> Die Gleichung $1^k + 2^k + \cdots + (m-1)^k = m^k$ hat für $k \geq 2$ keine natürlichzahlige Lösung $m \geq 2$.

Die einzige bekannte Lösung ist $(k=1, m=3)$: $1 + 2 = 3$.

---

## Mathematischer Hintergrund

### Die Gleichung

$$S_k(m-1) := \sum_{i=1}^{m-1} i^k = m^k$$

Für $k=1$ ist die Lösung $m=3$ analytisch beweisbar:

$$S_1(m-1) = \frac{m(m-1)}{2} = m \iff m-1 = 2 \iff m = 3$$

### Faulhabersche Formel und Bernoulli-Zahlen

Die Potenzsumme $S_k(N) = 1^k + \cdots + N^k$ hat eine exakte Formel:

$$S_k(N) = \frac{1}{k+1} \sum_{j=0}^{k} \binom{k+1}{j} B_j \, N^{k+1-j}$$

wobei $B_j$ die **Bernoulli-Zahlen** sind:

| $n$ | $B_n$ |
|-----|-------|
| 0 | $1$ |
| 1 | $+\tfrac{1}{2}$ (sympy-Konvention) |
| 2 | $\tfrac{1}{6}$ |
| 3 | $0$ |
| 4 | $-\tfrac{1}{30}$ |

**Wichtig**: Für ungerades $n \geq 3$ gilt $B_n = 0$.

### Moser's untere Schranke

Leo Moser bewies 1953:

> Falls $(k, m)$ mit $k \geq 2$ eine Lösung ist, dann gilt $m > 10^{10^6}$.

Gallot, Moree, Zudilin verbesserten dies 2011:

> $m > 2{,}7139 \times 10^{1{,}484 \times 10^9}$, und zusätzlich:
> - $k$ ist gerade
> - $2m-1 \mid k$
> - Alle Primteiler von $m$ sind Wieferich-Primzahlen bzgl. $m$

### Modulare Ausschlüsse

Für festes $k$ und einen Modulus $q$ kann man prüfen:

$$S_k(m-1) \equiv m^k \pmod{q}$$

Gilt die Kongruenz nicht, ist $m$ als Lösung ausgeschlossen. Je mehr Moduli testen, desto mehr Kandidaten können eliminiert werden.

---

## Klassenstruktur: `ErdosMoser`

### Attribute

| Attribut | Typ | Beschreibung |
|----------|-----|--------------|
| `MOSER_LOWER_BOUND_EXPONENT` | `int` | $10^6$ – Exponent in Moser's Schranke $m > 10^{10^6}$ |
| `_bernoulli_cache` | `dict` | Cache für sympy-Bernoulli-Zahlen |

### Methoden

#### Numerische Verifikation

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `is_solution` | `(k, m) → bool` | Prüft ob $S_k(m-1) = m^k$ gilt |
| `find_solutions` | `(k_max, m_max) → List[(k,m)]` | Exhaustive Suche nach Lösungen |
| `verify_k1_solution` | `() → dict` | Analytischer Beweis für $k=1$ |

#### Modulare Ausschlüsse

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `modular_residues` | `(k, q) → List[int]` | Zulässige Reste mod $q$ |
| `is_excluded_by_modulus` | `(k, m, q) → bool` | Ob $m$ durch $q$ ausgeschlossen |
| `find_excluding_moduli` | `(k, m, q_max) → List[int]` | Alle ausschließenden Moduli |
| `modular_exclusion_analysis` | `(k, m_range, q) → dict` | Effektivitäts-Analyse |

#### Bernoulli-Verbindung

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `bernoulli_number` | `(n) → Rational` | $B_n$ exakt via sympy |
| `faulhaber_formula` | `(N, k) → Rational` | $S_k(N)$ via Faulhaber |
| `faulhaber_equals_power` | `(k, m) → bool` | Symbolische Gleichheitsprüfung |
| `bernoulli_connection_summary` | `(k_max) → List[dict]` | Koeffizienten-Übersicht |

#### Moser's Schranke

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `mosers_lower_bound` | `() → dict` | Dokumentation der Schranken |
| `is_trivial_solution` | `(k, m) → bool` | Ob $(k,m) = (1,3)$ |
| `conjecture_status` | `() → dict` | Aktueller Beweisstand |

#### Summen-Analyse

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `power_sum_values` | `(k, m_max) → List[dict]` | Vergleich $S_k$ vs $m^k$ |
| `growth_comparison` | `(k, m_values) → dict` | Asymptotisches Wachstum |

---

## Verwendungsbeispiele

```python
from erdos_moser import ErdosMoser

em = ErdosMoser()

# Einzige bekannte Lösung
print(em.is_solution(1, 3))  # True: 1+2 = 3

# Kein Gegenbeispiel für k=2..5, m=2..500
solutions = em.find_solutions(k_max=5, m_max=500)
print(solutions)  # [(1, 3)]

# Analytischer Beweis für k=1
proof = em.verify_k1_solution()
print(proof['proof_steps'])

# Bernoulli-Zahlen
print(em.bernoulli_number(4))  # -1/30

# Faulhaber: S_3(5) = 1+8+27+64+125 = 225
print(em.faulhaber_formula(5, 3))  # 225

# Moser's Schranke
print(em.mosers_lower_bound())
```

---

## Asymptotisches Verhalten

Für große $m$ gilt:

$$S_k(m-1) \approx \frac{m^{k+1}}{k+1}$$

Da $\frac{m^{k+1}}{k+1} \gg m^k$ für großes $m$ und $k \geq 1$, müsste eine Lösung in einem sehr engen "Überkreuzungsbereich" liegen. Das macht Lösungen für große $m$ a priori unwahrscheinlich, ist aber kein Beweis.

---

## Wichtige Abhängigkeiten

| Bibliothek | Verwendung |
|-----------|------------|
| `sympy` | `bernoulli`, `binomial`, `Rational` – exakte symbolische Berechnungen |
| `numpy` | Implizit über Python-Summen |

---

## Beweisstand

| Aspekt | Status |
|--------|--------|
| $k=1$, $m=3$ ist Lösung | **BEWIESEN** (analytisch) |
| Keine andere Lösung für $k=2, m \leq 10^{10^6}$ | **BEWIESEN** (Moser 1953) |
| Keine Lösung für $k \geq 2$ allgemein | **OFFEN** – Conjecture |

---

## Literatur

- Leo Moser: *On the Diophantine equation $1^k + 2^k + \cdots + (m-1)^k = m^k$*, Scripta Math. 19, 1953
- Yves Gallot, Pieter Moree, Wadim Zudilin: *The Erdős-Moser equation $1^k + 2^k + \cdots + (m-1)^k = m^k$ revisited using continued fractions*, Math. Comp. 80, 2011
