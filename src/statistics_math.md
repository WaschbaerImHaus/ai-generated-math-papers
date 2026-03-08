# statistics_math.py – Statistik-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/statistics_math.py`

---

## Überblick

Dieses Modul implementiert deskriptive Statistik, Wahrscheinlichkeitsverteilungen, Hypothesentests, das Bayes-Theorem und Monte-Carlo-Simulationen. Es verwendet keine externen Statistikbibliotheken – alle Algorithmen sind in reinem Python implementiert.

### Inhaltsübersicht

| Kategorie | Funktionen |
|-----------|-----------|
| Deskriptive Statistik | `mean`, `median`, `mode`, `variance`, `std_dev`, `quartiles`, `iqr`, `skewness`, `kurtosis` |
| Normalverteilung | `normal_pdf`, `normal_cdf`, `normal_ppf` |
| Binomialverteilung | `binomial_pmf`, `binomial_cdf` |
| Poisson-Verteilung | `poisson_pmf`, `poisson_cdf` |
| Hypothesentests | `t_test_one_sample`, `t_test_two_sample`, `chi_square_test` |
| Bayes | `bayes_theorem` |
| Monte-Carlo | `monte_carlo_pi` |

---

## Hilfsfunktionen (intern)

### `_erf(x)` – Gauss'sche Fehlerfunktion
```
erf(x) = (2/√π) · ∫₀ˣ exp(-t²) dt
```
Approximation nach Abramowitz & Stegun 7.1.26. Fehler < 1.5×10⁻⁷.

### `_probit(p)` – Quantilfunktion der Standardnormalverteilung
Gibt z-Wert zurück mit `CDF(z) = p`. Rational-Approximation nach Peter Acklam (Fehler < 1.15×10⁻⁹).

### `_gamma_inc_reg(a, x)` – Regularisierte unvollständige Gammafunktion P(a,x)
- Für `x < a+1`: Reihenentwicklung
- Für `x ≥ a+1`: Kettenbruch (Lentz-Algorithmus)

### `_t_cdf(t, df)` – CDF der t-Verteilung
Verwendet regularisierte Betafunktion `I_x(df/2, 0.5)`.

### `_chi2_cdf(x, df)` – CDF der Chi-Quadrat-Verteilung
`P(X ≤ x) = P(df/2, x/2)` über regularisierte Gammafunktion.

---

## Deskriptive Statistik

### `mean(data) → float`
**Arithmetischer Mittelwert:**
```
μ = (1/n) · Σ xᵢ
```

### `median(data) → float`
**Median (Zentralwert):** Der mittlere Wert der sortierten Daten.
- Ungerade n: mittlerer Wert
- Gerade n: Durchschnitt der beiden mittleren Werte

**Eigenschaft:** Robuster gegen Ausreißer als der Mittelwert.

### `mode(data) → List[float]`
**Modalwert(e):** Der/die am häufigsten vorkommende(n) Wert(e).
Bei mehreren gleichhäufigen Werten werden alle zurückgegeben (multimodal).

### `variance(data, ddof=1) → float`
**Varianz:**
```
s² = (1/(n-ddof)) · Σ(xᵢ - μ)²
```

| ddof | Bedeutung |
|------|-----------|
| 0 | Populationsvarianz |
| 1 | Stichprobenvarianz (Bessel-Korrektur) |

Die Bessel-Korrektur (ddof=1) liefert einen **unverzerrten Schätzer** für die Populationsvarianz.

### `std_dev(data, ddof=1) → float`
**Standardabweichung:** `s = √variance(data, ddof)`

### `quartiles(data) → (Q1, Q2, Q3)`
**Quartile:** Teilen die sortierte Datenmenge in vier gleich große Teile.
- Q1 = Median der unteren Hälfte (25. Perzentil)
- Q2 = Median (50. Perzentil)
- Q3 = Median der oberen Hälfte (75. Perzentil)

### `iqr(data) → float`
**Interquartilsabstand (IQR):** `Q3 - Q1`
Robustes Streuungsmaß, resistent gegen Ausreißer.

### `skewness(data) → float`
**Schiefe (3. standardisiertes Moment):**
```
g₁ = (1/n · Σ(xᵢ-μ)³) / σ³
```

| Wert | Bedeutung |
|------|-----------|
| g₁ > 0 | Rechtsschief (langer rechter Schwanz) |
| g₁ = 0 | Symmetrisch |
| g₁ < 0 | Linksschief (langer linker Schwanz) |

### `kurtosis(data) → float`
**Exzess-Kurtosis (4. Moment - 3):**
```
g₂ = (1/n · Σ(xᵢ-μ)⁴) / σ⁴ - 3
```

| Wert | Typ | Bedeutung |
|------|-----|-----------|
| g₂ = 0 | Mesokurtisch | Normalverteilung |
| g₂ > 0 | Leptokurtisch | Spitzere Kurve, schwere Ränder |
| g₂ < 0 | Platykurtisch | Flachere Kurve |

---

## Normalverteilung N(μ, σ²)

### `normal_pdf(x, mu=0, sigma=1) → float`
**Wahrscheinlichkeitsdichtefunktion:**
```
f(x) = (1/(σ·√(2π))) · exp(-0.5·((x-μ)/σ)²)
```

Die "Glockenkurve" – die wichtigste Verteilung der Statistik (Zentraler Grenzwertsatz).

### `normal_cdf(x, mu=0, sigma=1) → float`
**Kumulative Verteilungsfunktion:** `P(X ≤ x)`
```
F(x) = 0.5 · (1 + erf((x-μ) / (σ·√2)))
```

### `normal_ppf(p, mu=0, sigma=1) → float`
**Quantilfunktion (Inverse CDF):** Gibt x zurück mit `P(X ≤ x) = p`
```
x = μ + σ·√2 · erfinv(2p-1)
```

---

## Binomialverteilung B(n, p)

Modelliert die Anzahl der Erfolge bei n unabhängigen Versuchen mit Erfolgswahrscheinlichkeit p.

### `binomial_pmf(k, n, p) → float`
**Wahrscheinlichkeitsmasse:**
```
P(X = k) = C(n,k) · pᵏ · (1-p)^{n-k}
```
Verwendet Log-Darstellung für numerische Stabilität bei großen n.

### `binomial_cdf(k, n, p) → float`
**Kumulative Wahrscheinlichkeit:** `P(X ≤ k) = Σᵢ₌₀ᵏ P(X=i)`

---

## Poisson-Verteilung Poi(λ)

Modelliert die Anzahl seltener Ereignisse in einem festen Intervall (Anrufe/Stunde, Fehler/Seite, ...).

### `poisson_pmf(k, lam) → float`
**Wahrscheinlichkeitsmasse:**
```
P(X = k) = (λᵏ · e⁻λ) / k!
```
Verwendet Log-Darstellung für numerische Stabilität.

### `poisson_cdf(k, lam) → float`
**Kumulative Wahrscheinlichkeit:** `P(X ≤ k)`

---

## Hypothesentests

### `t_test_one_sample(data, mu) → (t_stat, p_value)`

**Ein-Stichproben-t-Test:**
- H₀: Populationsmittelwert = μ
- H₁: Populationsmittelwert ≠ μ (zweiseitig)

**Teststatistik:**
```
t = (x̄ - μ) / (s / √n)
```
Freiheitsgrade: df = n-1. Zweiseitiger p-Wert: `2·P(T ≥ |t|)`

### `t_test_two_sample(data1, data2, equal_var=True) → (t_stat, p_value)`

**Zwei-Stichproben-t-Test:**
- H₀: μ₁ = μ₂

| `equal_var` | Test |
|-------------|------|
| True | Student-t-Test (gepoolte Varianz) |
| False | Welch-t-Test (Welch-Satterthwaite-Freiheitsgrade) |

**Gepoolte Varianz (Student):**
```
sp² = ((n₁-1)·s₁² + (n₂-1)·s₂²) / (n₁+n₂-2)
se = √(sp² · (1/n₁ + 1/n₂))
```

### `chi_square_test(observed, expected) → (chi2_stat, p_value)`

**Chi-Quadrat-Anpassungstest:** Vergleicht beobachtete vs. erwartete Häufigkeiten.
- H₀: Daten folgen der erwarteten Verteilung

**Teststatistik:**
```
χ² = Σ (Oᵢ - Eᵢ)² / Eᵢ
```
Freiheitsgrade: df = Anzahl Kategorien - 1.

---

## `bayes_theorem(prior, likelihood, evidence) → float`

**Bayes-Theorem:**
```
P(A|B) = P(B|A) · P(A) / P(B)
```

| Parameter | Symbol | Bedeutung |
|-----------|--------|-----------|
| `prior` | P(A) | A-priori-Wahrscheinlichkeit der Hypothese |
| `likelihood` | P(B\|A) | Likelihood der Beobachtung |
| `evidence` | P(B) | Gesamtwahrscheinlichkeit der Beobachtung |

---

## `monte_carlo_pi(n_samples, seed=None) → float`

**Monte-Carlo-Schätzung von π:**

Methode: Zufällige Punkte im Einheitsquadrat [0,1]² werden generiert.
Der Anteil innerhalb des Viertelkreises (x²+y²≤1) approximiert π/4.

```
π ≈ 4 · (Punkte im Kreis) / (Gesamtpunkte)
```

**Konvergenzrate:** O(1/√n) durch Gesetz der großen Zahlen.

---

## Verwendungsbeispiele

```python
from statistics_math import mean, std_dev, normal_pdf, t_test_one_sample, bayes_theorem

data = [2.1, 2.5, 2.3, 2.8, 2.4, 2.6, 2.2, 2.7]

print(mean(data))          # 2.45
print(std_dev(data))       # 0.232 (Stichprobe)

# Normalverteilung
print(normal_pdf(0))       # 0.3989 (Standardnormalverteilung im Maximum)
print(normal_cdf(1.96))    # ≈ 0.975 (97.5. Perzentil)

# Hypothesentest: Ist der Mittelwert = 2.5?
t, p = t_test_one_sample(data, mu=2.5)
print(f"t={t:.3f}, p={p:.3f}")   # Entscheidung ob H₀ abzulehnen

# Bayes-Theorem: P(Krank|Test positiv)
# Prior P(Krank) = 0.01, Sensitivität P(Pos|Krank) = 0.95, P(Pos) = 0.02
posterior = bayes_theorem(prior=0.01, likelihood=0.95, evidence=0.02)
print(f"P(Krank|positiv) = {posterior:.3f}")  # 0.475

# Monte-Carlo Pi
pi_approx = monte_carlo_pi(1_000_000, seed=42)
print(f"π ≈ {pi_approx:.5f}")
```
