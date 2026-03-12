# legendre_brocard.py — Dokumentation

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Version**: 1.0

---

## Übersicht

Das Modul `legendre_brocard.py` implementiert zwei klassische, bis heute ungeklärte Vermutungen der analytischen Zahlentheorie:

1. **Legendre-Vermutung** (1798)
2. **Brocard-Vermutung** (1904, Primzahlquadrate)

> **Hinweis**: Die hier implementierte Brocard-Vermutung ist die Konjektur über Primzahlen zwischen aufeinanderfolgenden Primzahlquadraten — **nicht** die Brocard-Ramanujan-Gleichung $n! + 1 = m^2$ (die ist in `brocard_extension.py`).

---

## Mathematischer Hintergrund

### 1. Legendre-Vermutung

**Conjecture (Legendre, 1798)**:
Für jede positive ganze Zahl $n \geq 1$ existiert mindestens eine Primzahl $p$ mit:

$$n^2 < p < (n+1)^2$$

Die Intervalllänge beträgt $(n+1)^2 - n^2 = 2n+1$ und wächst linear. Nach dem Primzahlsatz liegen im Intervall asymptotisch:

$$\pi((n+1)^2) - \pi(n^2) \approx \frac{2n+1}{\ln(n^2)} \approx \frac{n}{\ln n} \to \infty$$

Dies macht die Vermutung für große $n$ sehr plausibel, ist aber kein Beweis.

**Verbindung zu Bertrand's Postulat** (Theorem, Tschebyschew 1852):
Für alle $n \geq 1$ existiert eine Primzahl $p$ mit $n < p \leq 2n$.
Die Legendre-Vermutung ist **stärker** als Bertrand's Postulat, impliziert es aber nicht direkt (andere Intervallstruktur).

**Beweisstatus**: Offen. Numerisch verifiziert bis $n = 10^{10}$ (Oliveira e Silva, 2014).

**Verbindung zur Riemann-Hypothese**:
Unter RH gilt $g_n = O(\sqrt{p_n} \ln p_n)$, was Legendre für hinreichend große $n$ fast zeigt.

---

### 2. Brocard-Vermutung (1904)

**Conjecture (Henri Brocard, 1904)**:
Für alle $n \geq 2$ liegen zwischen $p_n^2$ und $p_{n+1}^2$ mindestens 4 Primzahlen:

$$\pi(p_{n+1}^2) - \pi(p_n^2) \geq 4 \quad \forall n \geq 2$$

**Sonderfall $n=1$**: Zwischen $2^2=4$ und $3^2=9$ liegen nur die Primzahlen 5 und 7, also 2. Die Vermutung gilt daher erst ab $n=2$.

**Erstes Beispiel ($n=2$)**: Zwischen $3^2=9$ und $5^2=25$ liegen: 11, 13, 17, 19, 23 → 5 Primzahlen ✓

**Beweisstatus**: Offen. Verifiziert für alle Primzahlen $p_n < 10^9$.

---

## Klassen

### `LegendreConjecture`

Verifikation und statistische Analyse der Legendre-Vermutung.

#### Wichtige Methoden

| Methode | Beschreibung |
|---------|-------------|
| `interval_primes(n)` | Gibt `LegendreInterval` mit allen Primzahlen in $(n^2, (n+1)^2)$ zurück |
| `verify_range(n_start, n_end)` | Verifiziert Vermutung im Bereich, gibt Statistik zurück |
| `bertrand_connection(n)` | Zeigt Verbindung zu Bertrand's Postulat für gegebenes $n$ |
| `max_gap_between_squares(n_start, n_end)` | Maximale Primzahllücke innerhalb der Quadrat-Intervalle |
| `gap_histogram_data(n_start, n_end)` | Histogramm der Primzahlanzahl pro Intervall |

#### Rückgabe von `verify_range()`

```python
{
    "verified": True,                # Kein Gegenbeispiel gefunden
    "counterexamples": [],           # Liste von n ohne Primzahl
    "min_count": 2,                  # Minimale Anzahl (meist bei kleinen n)
    "max_count": 15,                 # Maximale Anzahl
    "min_n": 2,                      # n wo Minimum liegt
    "max_n": 999,                    # n wo Maximum liegt
    "total_verified": 1000,          # Anzahl geprüfter Intervalle
    "elapsed_seconds": 0.42,
}
```

#### Datenklasse `LegendreInterval`

```python
@dataclass
class LegendreInterval:
    n: int          # Basiswert
    n_sq: int       # n²
    np1_sq: int     # (n+1)²
    primes: List[int]  # Primzahlen im Intervall (exklusiv)
    count: int      # Anzahl
```

---

### `BrocardConjecture7`

Verifikation der Brocard-Vermutung (1904) über Primzahlen zwischen Primzahlquadraten.

**Klassenattribut**: `MINIMUM_REQUIRED = 4`

#### Wichtige Methoden

| Methode | Beschreibung |
|---------|-------------|
| `interval_primes(n)` | Primzahlen in $(p_n^2, p_{n+1}^2)$ |
| `verify_range(n_start, n_end)` | Verifikation mit Statistik |
| `minimum_tracking(n_start, n_end)` | Verfolgt laufendes Minimum |
| `statistics(n_start, n_end)` | Vollständige Verteilungsstatistik |

#### Datenklasse `BrocardInterval`

```python
@dataclass
class BrocardInterval:
    n: int          # Index der Primzahl pₙ
    p_n: int        # pₙ
    p_n1: int       # pₙ₊₁
    primes: List[int]  # Primzahlen in (pₙ², pₙ₊₁²)
    count: int
```

---

## Bekannte Werte

### Legendre-Intervalle (kleine n)

| n | Intervall | Primzahlen | Anzahl |
|---|-----------|-----------|--------|
| 1 | (1, 4) | 2, 3 | 2 |
| 2 | (4, 9) | 5, 7 | 2 |
| 3 | (9, 16) | 11, 13 | 2 |
| 4 | (16, 25) | 17, 19, 23 | 3 |
| 5 | (25, 36) | 29, 31 | 2 |
| 10 | (100, 121) | 101, 103, 107, 109, 113 | 5 |

### Brocard-Intervalle (n ≥ 2)

| n | $p_n$ | $p_{n+1}$ | Intervall | Primzahlen | Anzahl |
|---|-------|----------|----------|-----------|--------|
| 1 | 2 | 3 | (4, 9) | 5, 7 | 2 (*) |
| 2 | 3 | 5 | (9, 25) | 11,13,17,19,23 | 5 |
| 3 | 5 | 7 | (25, 49) | 29,31,37,41,43,47 | 6 |
| 4 | 7 | 11 | (49, 121) | viele | ≥4 |

(*) n=1 ist der einzige Fall mit < 4; die Vermutung gilt erst ab n=2.

---

## Verwendungsbeispiele

```python
from legendre_brocard import LegendreConjecture, BrocardConjecture7

# Legendre-Vermutung
lc = LegendreConjecture()

# Primzahlen zwischen 100 und 121
iv = lc.interval_primes(10)
print(iv.primes)  # [101, 103, 107, 109, 113]

# Verifikation für n=1..1000
result = lc.verify_range(1, 1000)
print(result["verified"])  # True

# Brocard-Vermutung
bc = BrocardConjecture7()

# Primzahlen zwischen 9 (=3²) und 25 (=5²)
iv = bc.interval_primes(2)
print(iv.primes)  # [11, 13, 17, 19, 23]
print(iv.count)   # 5 ≥ 4 ✓

# Verifikation
result = bc.verify_range(2, 500)
print(result["verified"])  # True
print(result["min_count"]) # ≥ 4
```

---

## Mathematische Invarianten (als Tests verifiziert)

1. Alle Primzahlen in `interval_primes(n)` sind tatsächlich prim (sympy.isprime)
2. Alle Primzahlen liegen im offenen Intervall (exklusiv Grenzen)
3. Intervalllänge von `LegendreInterval`: `np1_sq - n_sq == 2*n + 1`
4. Laufendes Minimum in `minimum_tracking()` ist monoton fallend
5. `BrocardConjecture7.MINIMUM_REQUIRED == 4`

---

## Abhängigkeiten

- `sympy`: `isprime`, `primerange`, `prime`, `nextprime`, `primepi`
- `numpy`: Statistische Berechnungen
- Python Standard: `math`, `time`, `collections`, `dataclasses`

---

## Bezug zu anderen Modulen

| Modul | Bezug |
|-------|-------|
| `andrica_cramer.py` | Andrica-Vermutung ist äquivalent zu Legendre für aufeinanderfolgende Primzahlen |
| `proof_theory.py` | Bertrand's Postulat (Theorem) |
| `analytic_number_theory.py` | Primzahlsatz, Primzahldichte |
| `brocard_extension.py` | Brocard-Ramanujan-Gleichung (ANDERE Vermutung!) |
