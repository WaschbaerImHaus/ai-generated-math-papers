# normality_digits.py — Normalität von Zahlen: Ziffernanalyse und Borel-Theorie

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `src/normality_digits.py`

---

## Übersicht

Dieses Modul untersucht, ob mathematische Konstanten wie $\pi$, $e$ und $\sqrt{2}$
die Eigenschaft der **Normalität** besitzen — eine fundamentale zahlentheoretische
Eigenschaft, deren Nachweis für die meisten interessanten Konstanten noch aussteht.

---

## Was ist eine normale Zahl?

### Definition

Eine reelle Zahl $x$ heißt **normal zur Basis $b$**, wenn jede endliche Ziffernsequenz
der Länge $k$ in der $b$-adischen Entwicklung von $x$ mit der erwarteten Häufigkeit
$b^{-k}$ auftritt.

Eine Zahl heißt **absolut normal**, wenn sie normal zu **jeder** Basis $b \geq 2$ ist.

### Wichtige Sätze und Conjectures

| Aussage | Status |
|---|---|
| Fast alle reellen Zahlen sind normal (Borel 1909) | **BEWIESEN** |
| $C_{10} = 0.12345678910...$ ist normal zur Basis 10 | **BEWIESEN** (Champernowne 1933) |
| Copeland-Erdős-Konstante (Primzahlen) ist normal | **BEWIESEN** (1946) |
| $\pi$ ist normal | **CONJECTURE** (starke empirische Evidenz) |
| $e$ ist normal | **CONJECTURE** |
| $\sqrt{2}$ ist normal | **CONJECTURE** |
| $\ln(2)$ ist normal | **CONJECTURE** |

---

## Klasse `NormalNumberAnalysis`

Statistische Analyse von Ziffernfolgen auf Konsistenz mit Normalität.

### Unterstützte Konstanten

| Name | Konstante | Normalitätsstatus |
|---|---|---|
| `'pi'` | $\pi = 3.14159\ldots$ | CONJECTURE |
| `'e'` | $e = 2.71828\ldots$ | CONJECTURE |
| `'sqrt2'` | $\sqrt{2} = 1.41421\ldots$ | CONJECTURE |
| `'champernowne'` | $C_{10} = 0.12345\ldots$ | BEWIESEN (Basis 10) |

### Chi²-Goodness-of-Fit-Test

Testet die **Nullhypothese** $H_0$: Die Ziffern 0–9 sind gleichverteilt.

**Teststatistik**:
$$\chi^2 = \sum_{d=0}^{9} \frac{(O_d - E)^2}{E}, \quad E = \frac{n}{10}$$

Unter $H_0$ gilt asymptotisch: $\chi^2 \sim \chi^2(9)$ (9 Freiheitsgrade).

**Interpretation**:
- $p > 0.05$: Keine Evidenz gegen Normalität (H₀ nicht verworfen)
- $p < 0.05$: Schwache Evidenz gegen Gleichverteilung
- $p < 0.01$: Starke Evidenz gegen Gleichverteilung

### Runs-Test (Wald-Wolfowitz)

Testet die **Unabhängigkeit** aufeinanderfolgender Ziffern.

Ziffern werden binärisiert: $b_i = 1$ falls $d_i > 4$, sonst $0$.
Ein "Run" ist eine maximale Folge gleicher Symbole.

**Teststatistik**:
$$Z = \frac{R - \mu_R}{\sigma_R}, \quad \mu_R = \frac{2n_0 n_1}{n} + 1, \quad \sigma_R^2 = \frac{2n_0 n_1 (2n_0 n_1 - n)}{n^2(n-1)}$$

### Longest-Run-Test

Berechnet die **längste ununterbrochene Folge** gleicher Ziffern.

Theoretischer Erwartungswert in einer zufälligen Folge von $n$ Ziffern:
$$E[\text{longest run}] \approx \log_{10}(n)$$

### Methoden

| Methode | Beschreibung |
|---|---|
| `get_digits_of_pi(n)` | Erste $n$ Dezimalstellen von $\pi$ (mpmath) |
| `get_digits_of_e(n)` | Erste $n$ Dezimalstellen von $e$ |
| `get_digits_of_sqrt2(n)` | Erste $n$ Dezimalstellen von $\sqrt{2}$ |
| `get_digits_of_champernowne(n)` | Erste $n$ Stellen der Champernowne-Konstante |
| `digit_frequency(digits)` | Absolute und relative Häufigkeit 0–9 |
| `chi_squared_test(digits)` | $\chi^2$-Test auf Gleichverteilung |
| `runs_test(digits)` | Wald-Wolfowitz-Runs-Test |
| `longest_run_test(digits)` | Maximale Run-Längen-Analyse |
| `analyze_constant(name, n)` | Vollständige Analyse einer Konstante |

---

## Klasse `NormalityBounds`

Theoretische Schranken und Ergebnisse zur Borel-Normalitätstheorie.

### Satz von Borel (1909) — BEWIESEN

Die Menge der reellen Zahlen in $[0,1]$, die **nicht** absolut normal sind,
hat das **Lebesgue-Maß Null**.

**Beweisskizze**: Via Borel-Cantelli-Lemma und Gesetz der großen Zahlen.

**Wichtige Einschränkung**: Der Satz ist nicht-konstruktiv — er zeigt die Existenz
normaler Zahlen (fast alle!), sagt aber **nichts** über spezifische Zahlen wie $\pi$ aus.

### Diskrepanz

Die **Diskrepanz** misst die maximale Abweichung der empirischen Ziffernverteilung
von der theoretischen Gleichverteilung:

$$D_n = \max_{d \in \{0,\ldots,9\}} \left| \frac{\#\{i \leq n : a_i = d\}}{n} - \frac{1}{10} \right|$$

Für eine normale Zahl gilt: $D_n \to 0$ für $n \to \infty$.

**Erwartete Schranke** (Conjecture für $\pi$ etc.):
$$D_n = O\!\left(\sqrt{\frac{\log n}{n}}\right)$$

### Copeland-Erdős-Konstante

$$0.\underbrace{2}_{p_1}\underbrace{3}_{p_2}\underbrace{5}_{p_3}\underbrace{7}_{p_4}\underbrace{11}_{p_5}\underbrace{13}_{p_6}\ldots$$

**BEWIESEN** normal zur Basis 10 (Copeland & Erdős 1946).

Beweis nutzt den Primzahlsatz: Primzahlen sind hinreichend dicht.

### Methoden

| Methode | Beschreibung |
|---|---|
| `borel_measure_theorem()` | Beschreibung des Borel-Satzes |
| `compute_discrepancy(digits)` | Empirische Diskrepanz $D_n$ |
| `normality_bound_conjecture(n)` | Theoretische Schranken |
| `copeland_erdos_constant_digits(n)` | Ziffern der Copeland-Erdős-Konstante |

---

## Beispielverwendung

```python
from normality_digits import NormalNumberAnalysis, NormalityBounds

# Analyse von π
analysis = NormalNumberAnalysis(precision=1000)
result = analysis.analyze_constant('pi', n_digits=1000)
print(result['chi_squared_test']['p_value'])   # > 0.05 erwartet
print(result['chi_squared_test']['interpretation'])

# Champernowne (bewiesen normal)
result_champ = analysis.analyze_constant('champernowne', n_digits=500)
print(result_champ['normalcy_status'])  # "BEWIESEN normal zur Basis 10"

# Diskrepanz
bounds = NormalityBounds()
digits = analysis.get_digits_of_pi(500)
d = bounds.compute_discrepancy(digits)
print(f"Diskrepanz von π (500 Stellen): {d:.4f}")

# Borel-Theorie
borel = bounds.borel_measure_theorem()
print(borel['statement'])
```

---

## Warum ist Normalität schwer zu beweisen?

Die Vermutung, dass $\pi$ normal ist, besteht seit über 100 Jahren.
Das fundamentale Problem: Normalität ist eine **globale** Eigenschaft
(alle Stellen, asymptotisch), aber wir haben nur Zugang zu endlich vielen Stellen.

**Bekannte Hindernisse**:
1. Keine effektive Methode, die Verteilung in beliebig langen Blöcken zu kontrollieren
2. Analytische Methoden (Fourier, L-Funktionen) liefern nur Informationen über
   Summen, nicht über Gleichverteilung
3. Algebraische Eigenschaften von $\pi$, $e$ etc. (Transzendenz) helfen nicht direkt

**Aktueller Stand**: Alle numerischen Evidenzen bis Billiarden von Stellen
zeigen keine Abweichung von Normalität — aber kein Beweis existiert.

---

## Offene Probleme (Stand 2024)

1. Ist $\pi$ normal (zur Basis 10, zur Basis 2, absolut)? **(CONJECTURE)**
2. Ist $e$ normal? **(CONJECTURE)**
3. Ist $\sqrt{2}$ normal? **(CONJECTURE)**
4. Gibt es algebraisch irrationale Zahlen, die nicht normal sind? **(UNBEKANNT)**
5. Gilt $D_n(\pi) = O(\sqrt{\log n / n})$? **(CONJECTURE)**

---

## Quellen

- Borel, É. (1909). *Les probabilités dénombrables et leurs applications arithmétiques*
- Champernowne, D.G. (1933). *The construction of decimals normal in the scale of ten*
- Copeland, A.H. & Erdős, P. (1946). *Note on normal numbers*
- Bailey, D.H. & Crandall, R.E. (2002). *On the random character of fundamental constant expansions*
- Wald, A. & Wolfowitz, J. (1940). *On a test whether two samples are from the same population*
