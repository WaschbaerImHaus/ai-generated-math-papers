# zeta_odd_values.py – Dokumentation

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12

---

## Überblick

Dieses Modul analysiert die **Riemann-Zeta-Funktion an ungeraden ganzzahligen Stellen**
$\zeta(2k+1)$. Während die geraden Werte $\zeta(2k)$ durch Bernoulli-Zahlen vollständig
bekannt sind, ist die algebraische Natur von $\zeta(3), \zeta(5), \zeta(7), ...$ weitgehend
ein offenes Rätsel der Mathematik.

---

## Bekannte Resultate (Stand 2026)

| Wert | Numerisch | Irrationalitätsstatus |
|------|-----------|----------------------|
| $\zeta(3)$ | $1.2020569031...$ | **IRRATIONAL** (Apéry 1979) |
| $\zeta(5)$ | $1.0369277551...$ | **CONJECTURE** (unbewiesen) |
| $\zeta(7)$ | $1.0083492773...$ | **CONJECTURE** (unbewiesen) |
| $\zeta(9)$ | $1.0020083928...$ | **CONJECTURE** (unbewiesen) |
| $\zeta(11)$ | $1.0004941886...$ | **CONJECTURE** (unbewiesen) |

---

## Mathematische Hintergründe

### Gerade vs. ungerade Zeta-Werte

**Gerade Werte** (vollständig bekannt):
$$\zeta(2k) = \frac{(-1)^{k+1} (2\pi)^{2k} B_{2k}}{2 \cdot (2k)!}, \quad k \geq 1$$

**Ungerade Werte**: Keine vergleichbare geschlossene Formel bekannt.

### Apéry-Beweis für $\zeta(3)$ (1979)

Roger Apéry konstruierte eine Folge rationaler Zahlen $a_n/b_n$ mit:

$$\frac{a_n}{b_n} \to \zeta(3), \qquad \left|\zeta(3) - \frac{a_n}{b_n}\right| < \frac{1}{b_n^{1+\delta}}$$

für ein $\delta > 0$. Da $b_n$ ganzzahlig ist, folgt aus dem Irrationalitätskriterium,
dass $\zeta(3)$ irrational sein muss.

**Apéry-Rekursion** für den Nenner $b_n$:

$$b_n = \sum_{k=0}^n \binom{n}{k}^2 \binom{n+k}{k}^2$$

Konvergenzrate: $\left|\zeta(3) - a_n/b_n\right| \sim C \cdot (\sqrt{2}-1)^{4n}$

### Kettenbruch von $\zeta(3)$

$$\zeta(3) = [1; 4, 1, 18, 1, 1, 1, 4, 1, ...]$$

Der erste Koeffizient ist $a_0 = 1$, der zweite $a_1 = 4$.
Kein periodisches Muster bekannt.

### Hypergeometrische Darstellung von $\zeta(3)$ (Euler)

$$\zeta(3) = \frac{5}{2} \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{n^3 \binom{2n}{n}}$$

Konvergenzfaktor $\approx 1/4^n$ pro Term – deutlich schneller als die Dirichlet-Reihe.

---

## Bewiesene Theoreme

### Ball-Rivoal-Theorem (2001)

> **THEOREM** (Ball, Rivoal 2001 – vollständig bewiesen):
>
> Die Dimension des $\mathbb{Q}$-Vektorraums, der von
> $\{1, \zeta(3), \zeta(5), \ldots, \zeta(2n+1)\}$ aufgespannt wird, wächst mindestens wie
>
> $$\dim_{\mathbb{Q}} \text{Span}\{1, \zeta(3), ..., \zeta(2n+1)\} \geq \frac{(1+o(1)) \ln n}{1 + \ln 2}$$
>
> Insbesondere: **Unendlich viele $\zeta(2k+1)$ sind irrational.**

Die Irrationalität von $\zeta(5)$ selbst ist damit **nicht** bewiesen – nur, dass
unendlich viele der Werte irrational sind.

### Zudilin-Theorem (2004)

> **THEOREM** (Zudilin 2004 – vollständig bewiesen):
>
> Mindestens einer der vier Werte $\zeta(5), \zeta(7), \zeta(9), \zeta(11)$ ist irrational.

Zudilin identifiziert **nicht**, welcher der vier Werte irrational ist.

---

## Klasse: ZetaOddValues

### Methoden

| Methode | Beschreibung | Rückgabe |
|---------|-------------|---------|
| `zeta3()` | $\zeta(3)$ hochpräzise (50 Stellen) | `mpmath.mpf` |
| `zeta5()` | $\zeta(5)$ hochpräzise | `mpmath.mpf` |
| `zeta7()` | $\zeta(7)$ hochpräzise | `mpmath.mpf` |
| `compute_odd_zeta(k)` | $\zeta(2k+1)$ für $k \geq 1$ | `mpmath.mpf` |
| `apery_rational_sequence(n)` | Apéry-Approximationsfolge $(b_n, a_n/6)$ | `(int, Fraction)` |
| `apery_sequence_convergence(max_n)` | Konvergenzanalyse | `List[Tuple]` |
| `apery_continued_fraction(n)` | Kettenbruch-Koeffizienten von $\zeta(3)$ | `List[int]` |
| `evaluate_continued_fraction(coeffs)` | Wertet Kettenbruch aus | `float` |
| `zeta3_hypergeometric(n)` | Euler-hypergeometrische Formel | `float` |
| `zeta5_hypergeometric_approx(n)` | Dirichlet-Summe für $\zeta(5)$ | `float` |
| `zeta3_knopp_acceleration(n)` | Knopp-beschleunigte Berechnung | `float` |
| `ball_rivoal_theorem()` | Theorem-Details und Numerik | `dict` |
| `zudilin_theorem()` | Theorem-Details und Werte | `dict` |
| `verify_zeta3_digits(n)` | Verifikation von $\zeta(3)$ auf $n$ Stellen | `bool` |
| `apery_irrationality_measure()` | Irrationalitätsmaß $\mu(\zeta(3))$ | `dict` |
| `zeta_odd_table()` | Tabelle aller $\zeta(2k+1)$ bis $k=7$ | `List[Tuple]` |

---

## Irrationalitätsmaß von $\zeta(3)$

Das Irrationalitätsmaß $\mu(\alpha)$ gibt an, wie gut $\alpha$ durch Brüche approximierbar ist.
Für die Gleichung $|\alpha - p/q| < q^{-\mu}$ gilt:

| Konstante | Oberschranke $\mu$ | Quelle |
|-----------|-------------------|--------|
| $\zeta(3)$ | $\leq 5.513890$ | Rhin-Viola (2001) |
| $\pi$ | $\leq 7.103205$ | Salikhov (2008) |
| $e$ | $= 2$ | Hermite (1873) |
| Beliebige irrationale Zahl | $\geq 2$ | Dirichlet |

---

## Offene Fragen (Stand 2026)

- **Irrationalität von $\zeta(5), \zeta(7), ...$**: Einzeln unbewiesen
- **Transzendenz von $\zeta(3)$**: Unbewiesen (stärker als Irrationalität)
- **Algebraische Unabhängigkeit**: Sind $\zeta(3), \zeta(5), \zeta(7)$ algebraisch unabhängig?
- **Explizite Formel** für $\zeta(2k+1)$: Nicht bekannt

---

## Abhängigkeiten

- `mpmath`: Hochpräzise Berechnung von $\zeta(s)$
- `numpy`: Array-Operationen
- `fractions.Fraction`: Exakte rationale Arithmetik für Apéry-Folge
