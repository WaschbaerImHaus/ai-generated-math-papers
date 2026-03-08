# proof_theory.py – Beweistheorie-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/proof_theory.py`

---

## Überblick

Dieses Modul implementiert algorithmische Werkzeuge zur Untersuchung offener mathematischer Vermutungen. Es dient als Grundlage für das Fernziel, bedeutende mathematische Vermutungen zu beweisen oder zu widerlegen.

### Inhaltsübersicht

| Bereich | Inhalt |
|---------|--------|
| **Collatz-Vermutung** | Folgenberechnung, Stoppzeit, Verifikation |
| **Goldbach-Vermutung** | Zerlegungen, Verifikation |
| **Zwillingsprimzahlen** | Suche, Zählung, Konstante |
| **Riemann-Hypothese** | Partielle Zeta, Euler-Maclaurin, Nullstellensuche |
| **Formale Beweistechniken** | ProofByInduction-Klasse |
| **Primzahlalgorithmen** | Sieb, Miller-Rabin, Legendre, Jacobi, CRT |

---

## Collatz-Vermutung (3n+1-Problem)

**Vermutung (Lothar Collatz, 1937):** Für jede positive ganze Zahl n erreicht die folgende Folge nach endlich vielen Schritten die Zahl 1.

**Regel:**
```
f(n) = n/2      falls n gerade
f(n) = 3n + 1   falls n ungerade
```

**Status:** Offen. Verifiziert bis ca. 2^68 ≈ 2.95 × 10^20.

### `collatz_sequence(n) → list[int]`
Berechnet die vollständige Folge von n bis 1 (inklusiv).

### `collatz_stopping_time(n) → int`
Anzahl der Schritte bis 1 erreicht wird. Zahlen mit großen Stoppzeiten sind oft Kandidaten für Musteranalysen.

### `collatz_max_value(n) → int`
Gibt den Maximalwert zurück, den die Folge bei Startzahl n erreicht. Interessant: Für n=27 erreicht die Folge 9232 bevor sie fällt.

### `collatz_verify_range(limit) → dict`
Verifiziert die Vermutung empirisch für alle Zahlen 1 bis limit.

**Rückgabe:**
```python
{
    "verified": True,           # Kein Gegenbeispiel gefunden
    "counterexample": None,     # Erste Gegenbeispiel oder None
    "max_stopping_time": (n, steps),
    "max_value_seen": (n, maxval),
    "range_checked": limit
}
```

---

## Goldbach-Vermutung

**Vermutung (Christian Goldbach, 1742):** Jede gerade Zahl > 2 lässt sich als Summe zweier Primzahlen schreiben.

**Status:** Offen. Verifiziert bis 4 × 10^18.

**Beispiele:** 4 = 2+2, 6 = 3+3, 8 = 3+5, 10 = 3+7 = 5+5, ...

### `is_prime_fast(n) → bool`
Schneller deterministischer Primzahltest mit 6k±1-Optimierung.
Alle Primzahlen > 3 haben die Form 6k+1 oder 6k-1.

### `goldbach_decomposition(n) → Optional[tuple[int, int]]`
Findet **eine** Goldbach-Zerlegung von n (die mit dem kleinsten Primfaktor).

### `goldbach_all_decompositions(n) → list[tuple[int, int]]`
Gibt **alle** Goldbach-Zerlegungen `(p, q)` mit `p ≤ q` und `p + q = n` zurück.

### `goldbach_verify_range(limit) → dict`
Verifiziert die Goldbach-Vermutung für alle geraden Zahlen bis limit.

---

## Zwillingsprimzahl-Vermutung

**Vermutung:** Es gibt unendlich viele Zwillingsprimzahlpaare `(p, p+2)`.

**Bekannte Paare:** (3,5), (5,7), (11,13), (17,19), (29,31), ...

**Status:** Offen. Zhang (2013): Abstände < 246 kommen unendlich oft vor.

### `find_twin_primes(limit) → list[tuple[int, int]]`
Findet alle Zwillingsprimzahlpaare bis zur Schranke.

### `twin_prime_count(limit) → int`
Zählt Zwillingsprimzahlpaare bis zur Schranke.

**Hardy-Littlewood-Schätzung:**
```
π₂(x) ~ 2·C₂ · x / (ln x)²
```
mit der Zwillingsprimzahlkonstante C₂ ≈ 0.6601618...

### `twin_prime_constant() → float`
Berechnet die **Zwillingsprimzahlkonstante**:
```
C₂ = Π_{p prim, p>2} p(p-2)/(p-1)²  ≈ 0.6601618...
```
(Produkt über die ersten 1000 ungeraden Primzahlen)

---

## Riemann-Zeta-Funktion

### `riemann_zeta_partial(s, terms=10000) → complex`
Berechnet ζ(s) via **Partialsumme** `Σ 1/nˢ`.
**Warnung:** Konvergiert nur für `Re(s) > 1`.

### `riemann_zeta_euler_maclaurin(s, terms=100) → complex`
**Euler-Maclaurin-Beschleunigung** – konvergiert für `Re(s) > 0` (außer s=1):

```
ζ(s) ≈ Σ_{n=1}^{N} n^{-s} + N^{1-s}/(s-1) + N^{-s}/2 + Bernoulli-Korrekturen
```

Deutlich schneller als die pure Partialsumme.

### `check_riemann_hypothesis_numerically(imaginary_range, resolution=100) → dict`
Numerische Überprüfung der **Riemann-Hypothese** im angegebenen Bereich.

**Methode:** Berechnet |ζ(1/2 + it)| auf einem Gitter und sucht lokale Minima nahe Null.

**Bekannte Nullstellen** (Imaginärteile):
```
t₁ ≈ 14.1347, t₂ ≈ 21.0220, t₃ ≈ 25.0109, t₄ ≈ 30.4249, t₅ ≈ 32.9351, ...
```

---

## Klasse: `ProofByInduction`

Hilfsklasse für Beweise durch **vollständige Induktion**.

**Schema:**
```
1. Basisfall: P(n₀) beweisen
2. Induktionsschritt: P(n) → P(n+1) beweisen
3. Schluss: P(n) gilt für alle n ≥ n₀
```

### `verify_base_case(predicate, base) → bool`
Überprüft den Basisfall `P(base)`.

### `empirical_verify(predicate, start, end) → dict`
Empirische Überprüfung für `start ≤ n ≤ end`.

**⚠️ Wichtig:** Dies ist **kein formaler Beweis**! Empirische Verifikation kann nur Gegenbeispiele finden.

---

## Sieb des Eratosthenes

### `sieve_of_eratosthenes(limit) → list[int]`

Findet alle Primzahlen bis `limit`.

**Algorithmus (ca. 240 v. Chr.):**
```
1. Erstelle Bool-Array 2 bis limit, alle True
2. Für jede Primzahl p: markiere p², p²+p, ... als False
3. Übrig gebliebene True-Einträge sind Primzahlen
```

**Komplexität:** Zeit O(n log log n), Speicher O(n).

---

## Miller-Rabin-Primzahltest

### `miller_rabin_primality_test(n, rounds=20) → bool`

**Probabilistischer Primzahltest** – sehr effizient für große Zahlen.

**Grundlage:** Fermatscher Kleiner Satz: Wenn p prim, dann `a^{p-1} ≡ 1 (mod p)`.

**Deterministische Basen:**
```
{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
```
Mit diesen Basen: **deterministisch für alle n < 3.317 × 10^24**.

**Fehlerwahrscheinlichkeit** bei probabilistischer Nutzung: < 4^{-20} ≈ 10^{-12}.

**Algorithmus:**
1. Schreibe `n-1 = 2^r · d` (d ungerade)
2. Teste für jeden Zeugen `a`: Berechne `a^d mod n`
3. Quadriere r-1 mal und prüfe auf starken Primzeugen

---

## Zahlentheoretische Hilfsfunktionen

### `legendre_symbol(a, p) → int`

**Legendre-Symbol (a/p)** für ungerade Primzahl p:

```
(a/p) = 0   falls p | a
(a/p) = 1   falls a quadratischer Rest mod p   (∃x: x² ≡ a mod p)
(a/p) = -1  falls a quadratischer Nichtrest mod p
```

**Berechnung via Euler-Kriterium:**
```
(a/p) ≡ a^{(p-1)/2} (mod p)
```

### `jacobi_symbol(a, n) → int`

**Jacobi-Symbol (a/n)** – Verallgemeinerung des Legendre-Symbols auf zusammengesetzte n.

**Wichtig:** `(a/n) = 1` bedeutet **nicht**, dass a quadratischer Rest mod n ist!

Berechnung via **quadratischem Reziprozitätsgesetz**.

### `chinese_remainder_theorem(remainders, moduli) → int`

**Chinesischer Restsatz:** Löst simultane Kongruenzen:
```
x ≡ r₁ (mod m₁)
x ≡ r₂ (mod m₂)
...
x ≡ rₖ (mod mₖ)
```

**Voraussetzung:** Die Moduln mᵢ müssen paarweise koprim sein.

**Konstruktion (Bézout):**
```
M = m₁ · m₂ · ... · mₖ
Mᵢ = M / mᵢ
yᵢ = Mᵢ⁻¹ mod mᵢ
x = Σ rᵢ · Mᵢ · yᵢ (mod M)
```

**Anwendungen:** Kryptographie, Parallele Berechnungen, Beweis-Techniken.

---

## `conjecture_status_report() → dict`

Gibt eine Übersicht über alle bekannten offenen Vermutungen zurück.

**Millennium-Probleme (je 1 Million USD Preisgeld):**
- Riemann-Hypothese (1859) – offen
- P vs NP (1971) – offen
- Hodge-Vermutung (1950) – offen
- Yang-Mills (2000) – offen
- Navier-Stokes (2000) – offen
- BSD-Vermutung (1965) – offen
- Poincaré-Vermutung (1904) – **bewiesen** (Perelman 2003)

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `itertools` | Kombinatorische Hilfsfunktionen |
| `sympy` | Symbolische Mathematik |
| `numpy` | Numerische Arrays für Nullstellensuche |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from proof_theory import (
    collatz_sequence, goldbach_decomposition, find_twin_primes,
    miller_rabin_primality_test, chinese_remainder_theorem, sieve_of_eratosthenes
)

# Collatz-Folge für 27 (berühmt für ihre Länge)
seq = collatz_sequence(27)
print(f"Länge: {len(seq)}, Maximum: {max(seq)}")  # Länge: 112, Max: 9232

# Goldbach-Zerlegung
print(goldbach_decomposition(100))   # (3, 97) oder andere Zerlegung

# Zwillingsprimzahlen
twins = find_twin_primes(100)
print(twins)  # [(3,5),(5,7),(11,13),(17,19),(29,31),(41,43),(59,61),(71,73)]

# Miller-Rabin (deterministisch für n < 3.3×10^24)
print(miller_rabin_primality_test(999983))   # True (prim)
print(miller_rabin_primality_test(999984))   # False (= 2^5 × 3 × ...)

# Sieb
primes = sieve_of_eratosthenes(50)
print(primes)  # [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Chinesischer Restsatz: x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
x = chinese_remainder_theorem([2, 3, 2], [3, 5, 7])
print(x)   # 23 (= 23 mod 105)
```
