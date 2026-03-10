# Modul: Transzendenztheorie (`transcendence_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul implementiert Algorithmen und Demonstrationen zur Transzendenztheorie: Klassifikation
reeller Zahlen als rational, algebraisch oder transzendent, Liouville-Zahlen und
Irrationalitätsmaße, Demonstrationen des Lindemann-Weierstrass-Satzes, des Gelfond-Schneider-Satzes
und von Bakers Theorem sowie die algebraische Zahl als vollständige Python-Klasse.
Optional verwendet es SymPy für exakte symbolische Berechnungen.

## Mathematischer Hintergrund

### Klassifikation von Zahlen

$$\mathbb{Z} \subset \mathbb{Q} \subset \overline{\mathbb{Q}} \subset \mathbb{R} \subset \mathbb{C}$$

- **Rational** $x \in \mathbb{Q}$: $x = p/q$, $p, q \in \mathbb{Z}$, $q \neq 0$
- **Algebraisch** $x \in \overline{\mathbb{Q}}$: Nullstelle eines Polynoms $f \in \mathbb{Q}[t]$
- **Transzendent**: weder rational noch algebraisch

Der **algebraische Grad** von $x$ ist der Grad des Minimalpolynoms (irreduzibel über $\mathbb{Q}$).

### Liouville-Zahlen und Irrationalitätsmaß

Das **Irrationalitätsmaß** (Liouville-Exponent) von $x$:

$$\mu(x) = \sup\left\{\mu : \left|x - \frac{p}{q}\right| < \frac{1}{q^\mu} \text{ unendlich viele Lösungen } \frac{p}{q} \in \mathbb{Q}\right\}$$

Klassische Schranken:
- $\mu(x) = 1$ für alle rationalen $x$
- $\mu(x) \geq 2$ für alle irrationalen $x$ (**Dirichlet**)
- $\mu(x) = 2$ für quadratische Irrationale (Hurwitz, scharf)
- $\mu(x) = \infty$ für **Liouville-Zahlen** (transzendent!)

**Liouville-Zahl:**

$$L = \sum_{k=1}^{\infty} 10^{-k!} = 0.1100010000000000000000010\ldots$$

**Hurwitz-Theorem:** $|x - p/q| < 1/(\sqrt{5}\, q^2)$ hat unendlich viele Lösungen;
$\sqrt{5}$ ist scharf (der goldene Schnitt).

**Bekannte Irrationalitätsmaße:**
| Zahl | $\mu$ | Quelle |
|------|-------|--------|
| $\pi$ | $\geq 7.103$ | Salikhov (2008) |
| $e$ | $= 2$ (exakt) | Klassisch |
| $\ln 2$ | $\geq 3.57$ | Rukhadze (1987) |
| $\sqrt{2}$ | $= 2$ (exakt) | Lagrange (periodischer Kettenbruch) |
| $\zeta(3)$ | $\geq 5.51$ | Rhin-Viola (2001) |

### Lindemann-Weierstrass-Satz (1882/1885)

> Sind $\alpha_1, \ldots, \alpha_n \in \overline{\mathbb{Q}}$ paarweise verschieden, so sind
> $e^{\alpha_1}, \ldots, e^{\alpha_n}$ linear unabhängig über $\overline{\mathbb{Q}}$.

Wichtige Folgerungen:
- $e$ ist transzendent (Hermite, 1873): $\alpha = 1 \in \overline{\mathbb{Q}} \Rightarrow e^1 = e \notin \overline{\mathbb{Q}}$
- $\pi$ ist transzendent (Lindemann, 1882): $e^{i\pi} = -1 \in \overline{\mathbb{Q}}$, $i\pi \neq 0 \Rightarrow i\pi \notin \overline{\mathbb{Q}} \Rightarrow \pi \notin \overline{\mathbb{Q}}$
- $\sin(1), \cos(1), \ln(2)$ sind transzendent

### Gelfond-Schneider-Satz (1934)

> Sind $a, b \in \overline{\mathbb{Q}}$ mit $a \notin \{0, 1\}$ und $b \notin \mathbb{Q}$,
> dann ist $a^b$ transzendent.

Löst Hilberts 7. Problem. Konsequenzen:
- $2^{\sqrt{2}}$ (Gelfond-Schneider-Konstante) ist transzendent
- $e^\pi$ (Gelfond-Konstante) ist transzendent: $e^\pi = (-1)^{-i}$
- $i^i = e^{-\pi/2}$ ist transzendent

### Baker-Theorem (1966)

> Logarithmische Linearkombinationen: Sind $\alpha_1, \ldots, \alpha_n$ algebraisch $\neq 0,1$
> und $\beta_0, \beta_1, \ldots, \beta_n$ algebraisch, nicht alle null, dann:
>
> $$\beta_0 + \beta_1 \log\alpha_1 + \cdots + \beta_n \log\alpha_n \neq 0$$
> (und wenn nicht null, dann transzendent)

Alan Baker erhielt dafür 1970 die Fields-Medaille.
Anwendungen: effektive Lösung diophantischer Gleichungen, elliptische Logarithmen.

### Kreisteilungspolynome

$$\Phi_n(x) = \prod_{\substack{1 \leq k \leq n \\ \gcd(k,n)=1}} (x - e^{2\pi i k/n}), \qquad x^n - 1 = \prod_{d \mid n} \Phi_d(x)$$

---

## Klassen und Methoden

### `AlgebraicNumber`

**Beschreibung:** Stellt eine algebraische Zahl $\alpha$ durch ihr Minimalpolynom dar.
Berechnet Konjugierte, Norm und prüft ob $\alpha$ eine Einheitswurzel ist.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `__init__()` | `(poly_coeffs: list)` | Koeffizienten $[a_0, a_1, \ldots, a_n]$ aufsteigend; z.B. $[-2, 0, 1]$ für $\sqrt{2}$ |
| `degree()` | `() -> int` | Algebraischer Grad (= Grad des Minimalpolynoms) |
| `conjugates()` | `() -> list` | Alle Konjugierten (Nullstellen des Minimalpolynoms, numerisch) |
| `norm()` | `() -> complex` | $N(\alpha) = \prod_\sigma \sigma(\alpha)$ (Produkt aller Konjugierten) |
| `is_unit_root()` | `(tol: float = 1e-8) -> bool` | Prüft ob $|\alpha| = 1$ und $\alpha^n = 1$ für ein $n$ |

---

### Standalone-Funktionen: Klassifikation

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `is_algebraic_integer_check()` | `(x: float, max_degree: int = 6) -> bool` | Sucht normiertes ganzzahliges Polynom mit $x$ als Nullstelle |
| `minimal_polynomial()` | `(x: float, max_degree: int = 6, max_coeff: int = 10) -> Optional[list]` | Findet das Minimalpolynom von $x$ (kleinster Grad) |
| `algebraic_degree()` | `(x: float, max_degree: int = 8) -> int` | Algebraischer Grad (1 = rational, $>1$ = irrational algebraisch, $d+1$ = wahrsch. transzendent) |
| `classify_number()` | `(x: float) -> str` | Gibt `'rational'`, `'algebraic'` oder `'probably_transcendental'` zurück |

---

### Standalone-Funktionen: Liouville-Zahlen

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `liouville_number()` | `(n_terms: int = 10) -> float` | Berechnet $L = \sum_{k=1}^{n_\text{terms}} 10^{-k!}$ |
| `is_liouville_like()` | `(x: float, precision: int = 50) -> bool` | Heuristik: Irrationalitätsexponent $\geq$ `precision` |
| `liouville_approximation_exponent()` | `(x: float, max_q: int = 1000) -> float` | Schätzt $\mu(x) = \sup\{-\log|x-p/q|/\log q\}$ numerisch |

---

### Standalone-Funktionen: Klassische Sätze (Demonstrationen)

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `lindemann_weierstrass_demo()` | `() -> dict` | Beispiele zu L-W: $e$, $\pi$, $e^{\sqrt{2}}$, Folgerungen |
| `e_is_transcendental_demo()` | `() -> dict` | Numerische Verifikation + Beweisskizze (Hermite 1873) |
| `pi_is_transcendental_demo()` | `() -> dict` | Euler'sche Identität → $\pi$ transzendent (Lindemann 1882) |
| `gelfond_schneider_demo()` | `() -> dict` | $2^{\sqrt{2}}$, $e^\pi$, $i^i$ als Gelfond-Schneider-Beispiele |
| `baker_theorem_demo()` | `() -> dict` | Baker-Theorem: log. Linearkombinationen, Fields-Medaille 1970 |

---

### Standalone-Funktionen: Irrationalitätsmaße

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `irrationality_measure()` | `(x: float, max_q: int = 10000) -> float` | Numerische Schätzung von $\mu(x)$ |
| `known_irrationality_measures()` | `() -> dict` | Bekannte exakte/untere Schranken für $\mu(\pi)$, $\mu(e)$, usw. |
| `diophantine_approximation_quality()` | `(x: float, max_q: int = 1000) -> dict` | Hurwitz-Analyse: Brüche mit $|x-p/q| < 1/(\sqrt{5}\,q^2)$ |
| `continued_fraction_irrationality()` | `(x: float, n_terms: int = 20) -> dict` | Kettenbruch $[a_0; a_1, \ldots]$ + Irrationalitätsmaß-Schätzung via $\mu \approx 2 + \log(a_{k+1})/\log(q_k)$ |

---

### Standalone-Funktionen: Algebraische Strukturen

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `cyclotomic_polynomial()` | `(n: int) -> list` | $\Phi_n(x)$ via SymPy oder Fallback für Primzahlen |
| `number_field_degree()` | `(poly_coeffs: list) -> int` | Grad des Zahlkörpers $[\mathbb{Q}(\alpha):\mathbb{Q}]$ |

---

## Beispiele

```python
from transcendence_theory import (
    classify_number, minimal_polynomial, algebraic_degree,
    liouville_number, irrationality_measure,
    lindemann_weierstrass_demo, gelfond_schneider_demo,
    AlgebraicNumber, cyclotomic_polynomial
)
import math

# Klassifikation
print(classify_number(0.5))           # 'rational'
print(classify_number(math.sqrt(2)))  # 'algebraic'
print(classify_number(math.pi))       # 'probably_transcendental'
print(classify_number(math.e))        # 'probably_transcendental'

# Minimalpolynom von √2
poly = minimal_polynomial(math.sqrt(2))  # [-2, 0, 1] (x² - 2)
print(algebraic_degree(math.sqrt(2)))    # 2

# Algebraischer Grad von ³√2
print(algebraic_degree(2 ** (1/3)))      # 3

# Liouville-Zahl
L = liouville_number(10)  # ≈ 0.1100010000000000000000010...
print(irrationality_measure(L, max_q=500))   # sehr groß (Liouville-artig)
print(irrationality_measure(math.pi))         # ≈ 7-8 (Salikhov: ≥ 7.103)

# Gelfond-Schneider
gs = gelfond_schneider_demo()
print(gs['beispiele']['2^sqrt(2)']['wert'])  # ≈ 2.6651...

# Algebraische Zahl: √2 (Minpoly x²-2)
alpha = AlgebraicNumber([-2, 0, 1])
print(alpha.degree())      # 2
print(alpha.conjugates())  # [√2, -√2]
print(alpha.norm())        # -2.0 (= Produkt der Nullstellen)

# Kreisteilungspolynom Φ₅(x) = x⁴ + x³ + x² + x + 1
print(cyclotomic_polynomial(5))  # [1, 1, 1, 1, 1]

# Kettenbruch von π
from transcendence_theory import continued_fraction_irrationality
result = continued_fraction_irrationality(math.pi, n_terms=15)
print(result['kettenbruch'])  # [3, 7, 15, 1, 292, ...]
```

## Verweise

- Verwandte Module: `complex_analysis.py`, `analytic_number_theory.py`, `proof_theory.py`, `algebra.py`
- Literatur:
  - Baker: *Transcendental Number Theory* (Cambridge, 1975)
  - Nesterenko & Philippon (Hrsg.): *Introduction to Algebraic Independence Theory* (Springer, 2001)
  - Waldschmidt: *Diophantine Approximation on Linear Algebraic Groups* (Springer, 2000)
  - Hermite: *Sur la fonction exponentielle* (1873)
  - Lindemann: *Über die Zahl π* (1882)
  - Gelfond: *On Hilbert's seventh problem* (1934)
