# andrica_cramer.py — Dokumentation

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Version**: 1.0

---

## Übersicht

Das Modul `andrica_cramer.py` implementiert zwei fundamentale, bis heute ungeklärte Vermutungen über Primzahllücken:

1. **Andrica-Vermutung** (1985, Dorin Andrica)
2. **Cramér-Vermutung** (1936, Harald Cramér) inkl. **Granville-Korrektur** (1995)

Beide Vermutungen machen Aussagen über die Größe von Primzahllücken $g_n = p_{n+1} - p_n$ und stehen in enger mathematischer Beziehung zueinander sowie zur Riemann-Hypothese und der Legendre-Vermutung.

---

## Mathematischer Hintergrund

### 1. Andrica-Vermutung

**Conjecture (Andrica, 1985)**:
Für alle $n \geq 1$ gilt:

$$A_n = \sqrt{p_{n+1}} - \sqrt{p_n} < 1$$

**Äquivalente Formulierungen**:
1. $\sqrt{p_{n+1}} < \sqrt{p_n} + 1$
2. $p_{n+1} < (\sqrt{p_n} + 1)^2 = p_n + 2\sqrt{p_n} + 1$
3. $g_n = p_{n+1} - p_n < 2\sqrt{p_n} + 1$

Formulierung 3 zeigt: Andrica behauptet, dass die Primzahllücken sublinear in $\sqrt{p_n}$ wachsen.

**Bekanntes Maximum**: Das Maximum liegt bei $n=4$ (d.h. $p_4 = 7$, $p_5 = 11$):

$$A_4 = \sqrt{11} - \sqrt{7} \approx 0.6708$$

Für alle bekannten $n$ bleibt $A_n$ deutlich unter 1. Numerisch verifiziert bis $4 \times 10^{18}$.

**Verbindung zur Legendre-Vermutung**:
Die Andrica-Vermutung für AUFEINANDERFOLGENDE Primzahlen ist äquivalent zur Legendre-Vermutung. Wenn $A_n < 1$, dann liegt $p_{n+1}$ im Intervall $(p_n, (\sqrt{p_n}+1)^2)$.

---

### 2. Cramér-Vermutung

**Conjecture (Harald Cramér, 1936)**:

$$\limsup_{n \to \infty} \frac{g_n}{(\ln p_n)^2} = 1$$

Das bedeutet: Die Primzahllücken wachsen höchstens wie $(\ln p_n)^2$, und es gibt unendlich viele $n$ wo die Lücke nahe an $(\ln p_n)^2$ ist.

**Cramér's probabilistisches Modell** (Heuristik):
Cramér modellierte Primzahlen als unabhängige Zufallsereignisse mit Wahrscheinlichkeit $1/\ln n$. Daraus folgt die Voraussage, dass $g_n / (\ln p_n)^2$ asymptotisch exponentialverteilt ist mit Parameter 1.

---

### 3. Granville-Korrektur (1995)

**Conjecture (Andrew Granville, 1995)**:
Basierend auf einer verfeinerten Heuristik (Hardy-Littlewood-Primzahltuples) argumentierte Granville, dass der tatsächliche Limes Superior **größer als 1** ist:

$$\limsup_{n \to \infty} \frac{g_n}{(\ln p_n)^2} = 2e^{-\gamma} \approx 1.1229$$

Dabei ist $\gamma \approx 0.5772156649...$ die **Euler-Mascheroni-Konstante**.

**Die Granville-Konstante**:

$$C_G = 2e^{-\gamma} = 2 \cdot e^{-0.5772...} \approx 1.12292...$$

Diese Korrektur berücksichtigt, dass gerade Zahlen nie prim sind (außer 2), was die Wahrscheinlichkeit großer Lücken leicht erhöht.

---

## Hierarchie der Stärke

$$\text{Cramér} \Rightarrow \text{Andrica} \Rightarrow \text{Legendre} \Rightarrow \text{Bertrand}$$

- **Cramér** impliziert **Andrica**: Für große $p$: $(\ln p)^2 \ll 2\sqrt{p}+1$
- **Andrica** ist äquivalent zu **Legendre** (für aufeinanderfolgende Primzahlen)
- **Legendre** impliziert **Bertrand** (Primzahl in $(n, 2n]$)

---

## Klassen

### `AndricaConjecture`

Verifikation und Analyse der Andrica-Vermutung.

**Klassenattribute**:
- `KNOWN_MAX_VALUE`: $\sqrt{11} - \sqrt{7} \approx 0.6708$
- `KNOWN_MAX_N`: 4

#### Wichtige Methoden

| Methode | Beschreibung |
|---------|-------------|
| `compute(n)` | Berechnet $A_n = \sqrt{p_{n+1}} - \sqrt{p_n}$ für Index $n$ |
| `verify_range(n_start, n_end)` | Verifiziert $A_n < 1$ für alle $n$ im Bereich |
| `max_value_analysis(n_end)` | Verfolgt laufendes Maximum und neue Rekorde |
| `gap_statistics(n_start, n_end)` | Lücken vs. Andrica-Schranke |

#### Rückgabe von `compute(n)`

```python
AndricaValue(
    n=4,
    p_n=7,
    p_n1=11,
    gap=4,
    value=0.6708203932499369,  # sqrt(11) - sqrt(7)
    holds=True                  # < 1
)
```

---

### `CramerConjecture`

Statistische Analyse der Cramér-Vermutung über normierte Primzahllücken.

**Klassenattribute**:
- `GRANVILLE_CONSTANT`: $2e^{-\gamma} \approx 1.1229$
- `EULER_MASCHERONI`: $\gamma \approx 0.5772$

#### Wichtige Methoden

| Methode | Beschreibung |
|---------|-------------|
| `compute(n)` | Berechnet $C_n = g_n / (\ln p_n)^2$ und Granville-Verhältnis |
| `prime_gap_statistics(n_start, n_end)` | Umfassende Primzahllückenstatistik |
| `granville_correction(n_start, n_end)` | Analyse der Granville-Korrektur |
| `maxima_analysis(n_end)` | Record-Maxima der Cramér-Werte |
| `distribution_analysis(n_start, n_end)` | Verteilungsanalyse mit Quantilen |

#### Rückgabe von `prime_gap_statistics()`

```python
{
    "max_cramer_value": 1.0506,    # Maximum von g_n/(log pₙ)²
    "max_cramer_n": 30,            # Index wo Maximum liegt
    "max_cramer_primes": (113, 127),
    "mean_cramer": 0.317,          # Durchschnitt
    "std_cramer": 0.248,
    "granville_constant": 1.12292,
    "granville_exceedances": [30, ...],  # n wo C_n > Granville
    "max_granville_ratio": 0.935,
    "gap_distribution": {2: 35, 4: 40, ...},
    "total_analyzed": 1000,
    "elapsed_seconds": 0.18,
}
```

#### Datenklasse `CramerValue`

```python
@dataclass
class CramerValue:
    n: int
    p_n: int
    p_n1: int
    gap: int
    log_sq: float       # (log pₙ)²
    cramer_val: float   # g_n / (log pₙ)²
    granville_ratio: float  # cramer_val / (2e^{-γ})
```

---

## Bekannte Werte

### Andrica-Maxima (erste Rekorde)

| n | $p_n$ | $p_{n+1}$ | $g_n$ | $A_n$ |
|---|-------|----------|-------|-------|
| 1 | 2 | 3 | 1 | 0.4142 |
| 2 | 3 | 5 | 2 | 0.5174 |
| 3 | 5 | 7 | 2 | 0.4365 |
| 4 | 7 | 11 | 4 | **0.6708** ← globales Maximum |
| 9 | 23 | 29 | 6 | 0.5781 |

### Cramér-Maxima (erste Rekorde bis n=100)

| n | $p_n$ | $g_n$ | $C_n = g_n/(\ln p_n)^2$ |
|---|-------|-------|----------------------|
| 1 | 2 | 1 | 2.082 |
| 4 | 7 | 4 | 2.738 |
| 9 | 23 | 6 | 1.745 |
| 30 | 113 | 14 | ≈1.05 |

---

## Verwendungsbeispiele

```python
from andrica_cramer import AndricaConjecture, CramerConjecture, GRANVILLE_CONSTANT

# Andrica-Vermutung
ac = AndricaConjecture()

# Bekanntes Maximum prüfen
av = ac.compute(4)
print(f"A_4 = {av.value:.6f}")  # 0.670820
print(f"Gilt: {av.holds}")      # True

# Verifikation für n=1..1000
result = ac.verify_range(1, 1000)
print(result["verified"])      # True
print(result["max_value"])     # ≈ 0.6708

# Laufendes Maximum
analysis = ac.max_value_analysis(200)
for event in analysis["new_max_events"][:5]:
    print(f"n={event['n']}: A_n={event['value']:.4f} (p={event['p_n']})")

# Cramér-Vermutung
cc = CramerConjecture()

print(f"Granville-Konstante: {GRANVILLE_CONSTANT:.6f}")  # 1.122922

# Statistik für n=1..10000
stats = cc.prime_gap_statistics(1, 10000)
print(f"Max Cramér-Wert: {stats['max_cramer_value']:.4f}")
print(f"Granville-Überschreitungen: {len(stats['granville_exceedances'])}")

# Verteilung
dist = cc.distribution_analysis(1, 1000)
print(f"Median: {dist['quantiles']['q50']:.4f}")
print(f"99%-Quantil: {dist['quantiles']['q99']:.4f}")
```

---

## Mathematische Invarianten (als Tests verifiziert)

1. `AndricaConjecture.KNOWN_MAX_VALUE == sqrt(11) - sqrt(7) ≈ 0.6708`
2. `GRANVILLE_CONSTANT == 2 * exp(-euler_mascheroni) ≈ 1.1229`
3. Alle `AndricaValue.holds` sind `True` für n=1..1000
4. Maximaler Andrica-Wert liegt nahe an 0.6708 (bekanntes globales Maximum)
5. Alle Primzahllücken liegen unter Andrica-Schranke: $g_n < 2\sqrt{p_n}+1$
6. Record-Maxima der Cramér-Werte sind streng monoton steigend
7. Für $p \geq 127$: $(\ln p)^2 < 2\sqrt{p}+1$ (Cramér ist stärker als Andrica)

---

## Abhängigkeiten

- `sympy`: `prime`, `nextprime`, `primerange`
- `numpy`: Statistische Berechnungen, Quantile, Histogramme
- Python Standard: `math`, `time`, `dataclasses`, `collections`

---

## Beweisstatus (Stand 2026)

| Vermutung | Status | Verifiziert bis |
|-----------|--------|----------------|
| Andrica | Offen | $4 \times 10^{18}$ |
| Cramér | Offen | Heuristisch |
| Granville-Korrektur | Offen | Heuristisch |
| Baker-Harman-Pintz 2001 | **Bewiesen** | $g_n \ll p_n^{0.525}$ |
| Unter RH | **Bedingt bewiesen** | $g_n \ll \sqrt{p_n} \ln p_n$ |

---

## Bezug zu anderen Modulen

| Modul | Bezug |
|-------|-------|
| `legendre_brocard.py` | Legendre-Vermutung: äquivalent zu Andrica für Primzahlpaare |
| `analytic_number_theory.py` | Primzahlsatz, $\pi(x)$, Zeta-Funktion |
| `complex_analysis.py` | Riemann-$\zeta$, deren Nullstellen bedingen Lückenschranken |
| `millennium_problems.py` | Riemann-Hypothese: impliziert $g_n \ll \sqrt{p} \ln p$ |
