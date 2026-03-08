# analytic_number_theory.py – Analytische Zahlentheorie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/analytic_number_theory.py`

---

## Überblick

Dieses Modul implementiert die Kernwerkzeuge der **analytischen Zahlentheorie** – dem Gebiet, das Analysis und Zahlentheorie verbindet, um Aussagen über Primzahlverteilungen zu beweisen.

**Mathematisches Fernziel:** Goldbach-Vermutung: Jede gerade Zahl > 2 = p + q

### Inhaltsübersicht

| Bereich | Inhalt |
|---------|--------|
| **Primzahlzählung** | π(x), Li(x), PNT-Vergleich |
| **Tschebyschow-Funktionen** | θ(x), ψ(x), Λ(n) |
| **Dirichlet-Charaktere** | Hauptcharakter χ₀, L-Funktionen |
| **Siebmethoden** | Semiprimen, Chen-Zerlegung, k-fast-Primen |
| **Riemanns Explizite Formel** | ψ(x) aus Nullstellen von ζ |
| **Primzahl-Lücken** | Cramér-Vermutung |

---

## Primzahlzählfunktion und Logarithmisches Integral

### `prime_counting_function(x: float) → int`

Berechnet **π(x)** = Anzahl der Primzahlen ≤ x exakt via Sieb des Eratosthenes.

**Verwendung:** `bytearray` statt `list[bool]` für Speichereffizienz (8× weniger Speicher).

**Komplexität:** Zeit O(x log log x), Speicher O(x).

**Der Primzahlsatz** (Hadamard & de la Vallée Poussin, 1896):
```
π(x) ~ x / ln(x)    (asymptotisch, relativer Fehler → 0)
```

**Bessere Approximation** (Gauß/Riemann):
```
π(x) ≈ Li(x) = ∫₂ˣ dt/ln(t)
```

---

### `logarithmic_integral(x: float, terms=100) → float`

Berechnet **Li(x)** = ∫₂ˣ dt/ln(t).

Li(x) ist die **beste einfache Approximation** für π(x).
Bekanntes Ergebnis: Li(x) > π(x) für alle bekannten x (bis ~10^316 gilt dies immer).

**Berechnung:**
```
li(x) = γ + ln|ln(x)| + Σ_{n=1}^∞ (ln x)^n / (n · n!)
Li(x) = li(x) - li(2)
```

γ ≈ 0.5772156649... ist die **Euler-Mascheroni-Konstante**.

### `prime_number_theorem_comparison(x_values) → list[dict]`

Vergleicht π(x), x/ln(x) und Li(x) für gegebene x-Werte und zeigt die relativen Fehler.

**Typische Genauigkeiten:**
| x | Fehler x/ln(x) | Fehler Li(x) |
|---|----------------|-------------|
| 10^3 | ~13% | ~0.4% |
| 10^6 | ~6% | ~0.01% |
| 10^9 | ~4% | ~0.0001% |

---

## Tschebyschow-Funktionen

Diese Funktionen sind für Beweise des Primzahlsatzes und der Riemann-Hypothese fundamentaler als π(x).

### `von_mangoldt_function(n: int) → float`

**Von-Mangoldt-Funktion Λ(n):**
```
Λ(n) = ln(p)   falls n = pᵏ für eine Primzahl p und k ≥ 1
Λ(n) = 0       sonst
```

**Bedeutung:** Λ kodiert die Primzahlverteilung mit dem "natürlichen" Gewicht ln(p). Sie ist das zentrale Gewicht in Riemanns Expliziter Formel.

**Beispiele:**
- Λ(2) = ln(2), Λ(4) = ln(2), Λ(8) = ln(2) (Potenzen von 2)
- Λ(3) = ln(3), Λ(9) = ln(3)
- Λ(6) = 0 (6 = 2·3, zwei verschiedene Primteiler)

### `chebyshev_theta(x: float) → float`

**Erste Tschebyschow-Funktion:**
```
θ(x) = Σ_{p ≤ x} ln(p)   (Summe über Primzahlen)
```

**Äquivalenz zum Primzahlsatz:**
```
θ(x) ~ x   genau dann, wenn π(x) ~ x/ln(x)
```

**Riemann-Hypothese impliziert:**
```
|θ(x) - x| < c·√x·(ln x)²
```

θ(x) macht **Sprünge der Höhe ln(p)** an jeder Primzahl p.

### `chebyshev_psi(x: float) → float`

**Zweite Tschebyschow-Funktion:**
```
ψ(x) = Σ_{n=1}^{x} Λ(n) = Σ_{pᵏ ≤ x} ln(p)
```

**Riemanns Explizite Formel** verbindet ψ(x) direkt mit den Nullstellen ρ von ζ:
```
ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - ½·ln(1-x^{-2})
```

Dieser Zusammenhang ist das Herzstück der analytischen Zahlentheorie.

### `chebyshev_psi_vs_x(x) → dict`

Vergleicht ψ(x) mit x. Unter der RH gilt: `|ψ(x) - x| = O(√x · ln²x)`.

---

## Dirichlet-Charaktere und L-Funktionen

### `is_coprime(a, b) → bool`

Prüft ob `gcd(a, b) = 1` (teilerfremd).

### `principal_dirichlet_character(n, modulus) → int`

**Hauptcharakter χ₀ modulo q:**
```
χ₀(n) = 1   falls gcd(n, q) = 1
χ₀(n) = 0   falls gcd(n, q) > 1
```

**Dirichlet-Charaktere** sind periodische multiplikative Funktionen χ: ℤ → ℂ:
- χ(mn) = χ(m)·χ(n) (multiplikativ)
- χ(n+q) = χ(n) (periodisch)
- χ(n) = 0 wenn gcd(n,q) > 1

### `dirichlet_l_function_partial(s, character_values, modulus, terms=5000) → complex`

**Dirichlet-L-Funktion** via Partialsumme für Re(s) > 1:
```
L(s, χ) = Σ_{n=1}^∞ χ(n)/n^s
```

**Bedeutung des Dirichlet-Theorems:**
- L(1, χ) ≠ 0 für χ ≠ χ₀
- → Unendlich viele Primzahlen in jeder Restklasse gcd(r,q) = 1

### `dirichlet_prime_counting_in_residue_class(x, residue, modulus) → int`

Zählt Primzahlen `p ≤ x` mit `p ≡ r (mod q)`.

**Asymptotik (Gleichverteilung):**
```
π(x; q, r) ~ π(x) / φ(q)
```

Primzahlen verteilen sich gleichmäßig auf alle teilerfremden Restklassen.

---

## Siebmethoden für Goldbach-verwandte Probleme

### `is_semiprime(n) → bool`

Prüft ob n ein **Semiprim** (Produkt genau zweier Primzahlen, mit Vielfachheit) ist.

**Semiprimen:** 4=2·2, 6=2·3, 9=3·3, 10=2·5, ...

**Relevanz für Chen-Theorem:** n = p + m, wobei m prim oder semiprim.

### `chen_decomposition(n) → Optional[tuple]`

Sucht eine **Chen-Zerlegung** n = p + m.

**Chen-Theorem (Chen Jingrun, 1966):**
Jede genügend große gerade Zahl n ist darstellbar als:
```
n = p + m,   wobei p prim und m prim oder semiprim
```

Dies ist das **stärkste bewiesene Ergebnis** in Richtung der Goldbach-Vermutung!

| Typ | Beschreibung |
|-----|-------------|
| "prime" | m ist prim (Goldbach-Zerlegung, "1+1") |
| "semiprime" | m ist semiprim (Chen-Zerlegung, "1+2") |

**Goldbach wäre "1+1" – Chen beweist "1+2".**

### `count_omega(n) → int`

Berechnet `ω(n)` = Anzahl der **verschiedenen** Primteiler von n.

```
ω(12) = ω(2²·3) = 2    (Primteiler: 2, 3)
ω(30) = ω(2·3·5) = 3
```

### `almost_prime(n, k) → bool`

Prüft ob n ein **k-fast-Prim** ist: `Ω(n) = k` (Gesamtanzahl Primfaktoren mit Vielfachheit).

| k | Bezeichnung |
|---|-------------|
| 1 | Primzahl |
| 2 | Semiprim |
| 3 | Triprime |

**Chen:** Jede gerade Zahl = Prim + 1-fast-Prim ODER Prim + 2-fast-Prim.

---

## Riemanns Explizite Formel

### `explicit_formula_psi(x, zeros_t, terms=10) → dict`

Berechnet die **Explizite Formel** für ψ(x):

**Riemanns Explizite Formel (1859):**
```
ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - ½·ln(1-x^{-2})
```

Die Summe läuft über alle nicht-trivialen Nullstellen ρ = 1/2 + it.

**Paarweise Vereinfachung** (für ρ = 1/2 + it und ρ̄ = 1/2 - it):
```
x^ρ/ρ + x^{ρ̄}/ρ̄ = 2·Re(x^ρ/ρ)
```

**Deutung:** Die Nullstellen der Zeta-Funktion erzeugen **Schwingungen** um den Hauptterm x. Mehr Nullstellen → bessere Approximation.

**Rückgabe:**
```python
{
    "x": x,
    "psi_exact": 5.99...,
    "psi_explicit_formula": 6.01...,
    "error": 0.02,
    "zeros_used": 10,
    "note": "Mehr Nullstellen → bessere Approximation"
}
```

---

## Primzahl-Lücken

### `prime_gaps(limit) → list[dict]`

Berechnet alle Primzahl-Lücken (Abstände aufeinanderfolgender Primzahlen) bis `limit`.

**Cramér-Vermutung:** Für große p gilt: `p_{n+1} - p_n = O((ln p_n)²)`

**Bekannte große Lücken:**
| Primzahl p | Nächste Prim | Lücke |
|-----------|-------------|-------|
| 23 | 29 | 6 |
| 89 | 97 | 8 |
| 113 | 127 | 14 |

**Rückgabe pro Lücke:**
```python
{
    "prime": 23,
    "next_prime": 29,
    "gap": 6,
    "cramer_bound": (ln 23)² ≈ 10.02   # Cramér-Schranke
}
```

### `maximal_prime_gap(limit) → dict`

Findet die größte Primzahl-Lücke bis zur Schranke.

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | log, sqrt, gcd |
| `cmath` | Komplexe Exponentialfunktion |
| `numpy` | Linspace (könnte entfernt werden) |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from analytic_number_theory import (
    prime_counting_function, logarithmic_integral,
    chebyshev_psi, von_mangoldt_function,
    chen_decomposition, explicit_formula_psi, prime_gaps
)

# Primzahlsatz-Vergleich
pi_1000 = prime_counting_function(1000)
li_1000 = logarithmic_integral(1000)
print(f"π(1000) = {pi_1000}")   # 168
print(f"Li(1000) = {li_1000:.2f}")  # ≈ 177.6 (etwas zu groß)

# Tschebyschow-Funktionen
psi_100 = chebyshev_psi(100)
print(f"ψ(100) = {psi_100:.4f}")    # ≈ 94.0
print(f"ψ(100)/100 = {psi_100/100:.4f}")  # ≈ 0.94 (nahe 1, PNT)

# Von-Mangoldt-Funktion
import math
print(von_mangoldt_function(8))   # ln(2) ≈ 0.693 (8 = 2³)
print(von_mangoldt_function(6))   # 0.0   (6 = 2·3, nicht Primzahlpotenz)

# Chen-Zerlegung
result = chen_decomposition(20)
print(result)  # (3, 17, 'prime') oder (7, 13, 'prime')

# Explizite Formel (bekannte Nullstellen)
known_zeros = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351]
formula = explicit_formula_psi(x=100, zeros_t=known_zeros, terms=5)
print(f"ψ(100) exakt: {formula['psi_exact']:.2f}")
print(f"ψ(100) Formel: {formula['psi_explicit_formula']:.2f}")

# Primzahl-Lücken
gaps = prime_gaps(100)
max_gap = max(gaps, key=lambda g: g["gap"])
print(f"Größte Lücke: {max_gap['gap']} zwischen {max_gap['prime']} und {max_gap['next_prime']}")
# → Lücke 8 zwischen 89 und 97
```

---

## Mathematische Verbindungen

### Primzahlsatz ↔ Zeta-Funktion
Der Primzahlsatz `π(x) ~ x/ln(x)` ist äquivalent dazu, dass ζ(s) auf der Geraden Re(s) = 1 keine Nullstellen hat – das war das zentrale Ergebnis von Hadamard und de la Vallée Poussin (1896).

### Riemann-Hypothese → stärkerer Primzahlsatz
Falls RH gilt:
```
|π(x) - Li(x)| = O(√x · ln x)
```
Das wäre die bestmögliche Schranke.

### Chen-Theorem → Weg zu Goldbach
Chen (1966) bewies sein Theorem mithilfe des **Selberg-Siebs**. Die Goldbach-Vermutung wäre bewiesen, wenn man das "1+2" zu "1+1" verbessern könnte. Das ist das Ziel aktiver mathematischer Forschung.
