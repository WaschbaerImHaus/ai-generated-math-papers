# twin_prime_analysis.py — Zwillingsprimzahl-Analyse

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Datum**: 2026-03-12

---

## Überblick

Dieses Modul analysiert **Zwillingsprimzahlpaare** $(p, p+2)$ umfassend:
Zählung, Dichte, Vergleich mit theoretischen Vorhersagen und moderne Durchbrüche.

---

## Mathematischer Hintergrund

### Zwillingsprimzahl-Vermutung (Conjecture — OFFEN)

> Es gibt unendlich viele Paare $(p, p+2)$, wobei beide Primzahlen sind.

Bekannte Paare: $(3,5), (5,7), (11,13), (17,19), (29,31), (41,43), \ldots$

**Status**: Trotz massiver numerischer Evidenz und modernen Durchbrüchen — **unbewiesen**.

### Hardy-Littlewood Vermutung B (Conjecture — OFFEN)

$$\pi_2(x) \sim 2 C_2 \cdot \frac{x}{(\ln x)^2}$$

mit der **Zwillingsprimzahl-Konstante**:
$$C_2 = \prod_{p \geq 3, \, p \text{ prim}} \frac{p(p-2)}{(p-1)^2} \approx 0.6601618158\ldots$$

### Brun-Konstante $B_2$ (BEWIESEN, Brun 1919)

$$B_2 = \sum_{\substack{(p, p+2) \text{ Zwillingsprim}}} \left(\frac{1}{p} + \frac{1}{p+2}\right) \approx 1.9021605831\ldots$$

Die Konvergenz dieser Reihe ist bewiesen — **unabhängig davon**, ob es unendlich viele Zwillingsprimpaare gibt.

---

## Moderne Durchbrüche

### Bounded Gaps zwischen Primzahlen

| Ergebnis | Autor | Jahr | Schranke |
|---------|-------|------|---------|
| $\liminf (p_{n+1} - p_n) < 70\,000\,000$ | Yitang Zhang | 2013 | **BEWIESEN** |
| $\liminf (p_{n+1} - p_n) < 600$ | James Maynard | 2013 | **BEWIESEN** |
| $\liminf (p_{n+1} - p_n) \leq 246$ | Polymath8b | 2014 | **BEWIESEN** |
| $\liminf (p_{n+1} - p_n) = 2$ | Vermutung | offen | **Conjecture** |

### GPY-Siebtechnik (Goldston-Pintz-Yıldırım, 2005)

Zeigt: Für jedes $\varepsilon > 0$ gibt es $\liminf (p_{n+1} - p_n) / \ln p_n < \varepsilon$.
Das bedeutet: Primzahlen häufen sich "beliebig eng" relativ zu ihrem Logarithmus.

**Schranke von 2 (Zwillingsprimzahlen) nicht erreichbar** wegen des **Selberg Parity Problems**.

### Parity Problem (Selberg)

Siebmethoden können nicht direkt zwischen:
- $n = p$ (prim) und $n = p_1 \cdot p_2$ (Semiprim)

unterscheiden, wenn beide die gleiche Siebstruktur haben. Dies "blockiert" den
direkten Beweis der Zwillingsprimzahl-Vermutung via Siebmethoden.

---

## Klassen

### `TwinPrimeAnalysis`

| Methode | Beschreibung |
|---------|-------------|
| `get_twin_pairs()` | Alle $(p, p+2)$ bis `limit` |
| `count_twin_primes(x)` | $\pi_2(x)$ = Anzahl der Paare bis $x$ |
| `hardy_littlewood_prediction(x)` | $2 C_2 x / (\ln x)^2$ |
| `compute_C2_numerically(limit)` | Eulerprodukt $\prod_{p \geq 3} p(p-2)/(p-1)^2$ |
| `prime_gaps_of_size(gap)` | Alle $p$ mit $p_{\text{next}} - p = \text{gap}$ |
| `gap_statistics(max_gap)` | Häufigkeit jeder Lückengröße |
| `cramer_model_prediction(x, gap)` | Cramér-Modell-Vorhersage |
| `zhang_maynard_summary()` | Historische Ergebnisse (Zhang, Maynard, Polymath) |

### `PrimeSieveGaps`

| Methode | Beschreibung |
|---------|-------------|
| `get_all_gaps()` | Alle $(p_n, p_{n+1} - p_n)$ |
| `gaps_by_size(size)` | Alle $p$ mit Lücke = size |
| `cramers_model_comparison()` | Tatsächlich vs. Cramér-Vorhersage |
| `erdos_rankin_lower_bound(x)` | Untere Schranke für max. Lücke |
| `maximal_gap_record()` | Größte Primzahllücke bis `limit` |

---

## Bekannte Werte

| $x$ | $\pi_2(x)$ | H-L Vorhersage |
|-----|-----------|----------------|
| $10^3$ | 35 | ~27 |
| $10^4$ | 205 | ~156 |
| $10^5$ | 1224 | ~827 |
| $10^6$ | 8169 | ~5765 |

---

## Cramér-Modell vs. Realität

Das **Cramér-Modell** (1936, heuristisch) behandelt Primzahlen als "unabhängige
Zufallsvariablen" mit Wahrscheinlichkeit $1/\ln n$.

**Vorhersage für mittlere Lücke** um $x$: $\approx \ln x$

**Vorhersage für maximale Lücke** bis $x$: $\approx (\ln x)^2$ (Cramér-Vermutung, **Conjecture**)

Tatsächlich ist die maximale Lücke bis $10^6$ gleich **148** (nach $p = 492113$),
während $(\ln 10^6)^2 \approx 191$. Das Cramér-Modell überschätzt leicht.
