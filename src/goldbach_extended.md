# goldbach_extended.py — Erweiterte Goldbach-Analysen

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Datum**: 2026-03-12

---

## Überblick

Dieses Modul implementiert erweiterte Goldbach-Analysen:

| Aussage | Status |
|---------|--------|
| Binäre Goldbach-Vermutung: jede gerade $n \geq 4 = p + q$ | **Conjecture — OFFEN** |
| Ternäre Goldbach-Vermutung: jede ungerade $n \geq 7 = p + q + r$ | **BEWIESEN** (Helfgott 2013) |
| Chen's Theorem: jede große gerade $n = p + P_2$ | **BEWIESEN** (Chen 1973) |
| Vinogradov-Satz: jede hinreichend große ungerade $n = p + q + r$ | **BEWIESEN** (Vinogradov 1937) |

---

## Binäre Goldbach-Vermutung

### Statement (Conjecture, OFFEN seit 1742)

> Jede gerade Zahl $n \geq 4$ ist die Summe zweier Primzahlen.

**Numerische Verifikation**: Bis $4 \times 10^{18}$ (Oliveira e Silva et al., 2012).

### Goldbach-Komet

$G(n) = \#\{(p, q) : p + q = n, \, p \leq q, \, p, q \text{ prim}\}$ für gerade $n$ zeigt
eine "kometenartige" Struktur, die bei $n$ mit vielen Primfaktoren höher liegt.

### Hardy-Littlewood-Heuristik (Conjecture A)

$$G(n) \sim 2\Pi_2 \cdot \prod_{\substack{p > 2 \\ p \mid n}} \frac{p-1}{p-2} \cdot \frac{n}{(\ln n)^2}$$

mit $\Pi_2 = \prod_{p > 2} \frac{p(p-2)}{(p-1)^2} \approx 0.6601618\ldots$ (Zwillingsprimkonstante).

---

## Ternäre Goldbach-Vermutung (BEWIESEN)

### Theorem (Helfgott 2013)

> Jede ungerade Zahl $n \geq 7$ ist die Summe dreier Primzahlen.

**Beweismethode**:
1. **Kreismethode** (Hardy-Littlewood): Darstellungsfunktion als Fourier-Integral
2. **Numerische Verifikation** für $n < 8.875 \times 10^{30}$ (Computer-gestützt)
3. **Analytischer Beweis** für $n > 8.875 \times 10^{30}$ via Exponentailsummen

**Wichtige Konsequenz**: Für $n = 5$: $5 = 2 + 3$ (nur 2 Summanden nötig — aber auch $5 = 2 + 2 + ?$ scheitert da $5-4=1$ nicht prim). Für $n \geq 7$: immer 3 Summanden möglich.

**Gegenbeispiel prüfen**: $9 \neq 2+3+4$, da $4$ nicht prim. Aber $9 = 2+2+5 \checkmark$ oder $9 = 3+3+3 \checkmark$.

---

## Chen's Theorem (BEWIESEN)

### Theorem (Jingrun Chen, 1973)

> Jede hinreichend große gerade Zahl $n$ ist darstellbar als $n = p + m$,
> wobei $p$ prim und $m$ entweder prim oder Produkt zweier Primzahlen ist.

$m$ ist ein sogenanntes **$P_2$-Semiprim** (höchstens 2 Primfaktoren).

**Bedeutung**: Dies ist das stärkste bisher bewiesene Ergebnis in Richtung der binären Goldbach-Vermutung.

---

## Vinogradov-Methode

### Theorem (Vinogradov, 1937 — BEWIESEN)

> Jede hinreichend große ungerade Zahl $N$ ist Summe dreier Primzahlen.

**Klassische Schranke**: "hinreichend groß" bedeutet $N > e^{e^{11.503}} \approx 10^{43000}$.
Helfgott (2013) reduzierte dies auf alle $N \geq 7$.

### Kreismethode

Die **Darstellungsfunktion**:
$$r_3(N) = \sum_{\substack{p_1 + p_2 + p_3 = N \\ p_i \text{ prim}}} \ln p_1 \cdot \ln p_2 \cdot \ln p_3$$

wird via Fourier-Inversion dargestellt:
$$r_3(N) = \int_0^1 S(\alpha)^3 \cdot e^{-2\pi i N \alpha} \, d\alpha$$

mit der **Exponentialsumme**:
$$S(\alpha) = \sum_{p \leq N} \Lambda(p) \cdot e^{2\pi i \alpha p}$$

### Major/Minor-Bögen

- **Major-Bögen** $\mathfrak{M}$: $\alpha$ nahe rationalen $a/q$ mit $q \leq Q = N^{1/3}$
  → $S(\alpha)$ ist groß, liefert **Hauptterm**
- **Minor-Bögen** $\mathfrak{m}$: Restliches $[0,1] \setminus \mathfrak{M}$
  → $S(\alpha)$ ist klein via Vinogradov-Schranken

**Hauptterm**:
$$\int_{\mathfrak{M}} S(\alpha)^3 e^{-2\pi i N\alpha} d\alpha \sim \mathfrak{S}(N) \cdot \frac{N^2}{2(\ln N)^3}$$

mit der **singulären Reihe** $\mathfrak{S}(N) \geq c > 0$ für ungerades $N$.

---

## Schnirelmann-Konstante

### Theorem (Schnirelmann, 1930 — BEWIESEN)

> Jede positive ganze Zahl $\geq 2$ ist die Summe von höchstens $s$ Primzahlen,
> wobei $s < \infty$ eine absolute Konstante ist.

Nach Helfgott (2013): $s \leq 4$ (4 Primzahlen reichen für alle $n \geq 2$).

---

## Klassen

### `GoldbachExtended`

| Methode | Beschreibung |
|---------|-------------|
| `goldbach_representations(n)` | Alle $(p,q)$ mit $p+q=n$, $p \leq q$ prim |
| `goldbach_count(n)` | $G(n)$ |
| `goldbach_comet(n_max)` | $(n, G(n))$ für alle geraden $n \leq n_{\max}$ |
| `hardy_littlewood_goldbach(n)` | H-L-Heuristik für $G(n)$ |
| `check_ternary_goldbach(n)` | Findet $(p,q,r)$ mit $p+q+r=n$ |
| `verify_ternary_goldbach_range(n_max)` | Verifikation für alle ungeraden $n \leq n_{\max}$ |
| `chen_theorem_check(n)` | Chen-Darstellung $n = p + P_2$ |
| `schnirelmann_density(n_limit)` | Minimale Primanzahl für kleine Zahlen |
| `goldbach_conjecture_status()` | Status-Übersicht aller Goldbachs |

### `VinogradovMethod`

| Methode | Beschreibung |
|---------|-------------|
| `exponential_sum(alpha, N)` | $S(\alpha) = \sum_{p \leq N} \Lambda(p) e^{2\pi i\alpha p}$ |
| `estimate_r3(N)` | Numerische Schätzung von $r_3(N)$ |
| `major_arc_definition(N)` | Definition der Major-Bögen |
| `singular_series(N)` | Singuläre Reihe $\mathfrak{S}(N)$ |
| `minor_arc_bound_explanation()` | Erklärung der Minor-Bogen-Abschätzung |
