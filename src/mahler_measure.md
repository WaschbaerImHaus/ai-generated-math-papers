# mahler_measure.py – Mahler-Maß von Polynomen

**Autor:** Michael Fuhrmann
**Erstellt:** 2026-03-12
**Letzte Änderung:** 2026-03-12
**Build:** 122

---

## Überblick

Das Mahler-Maß ist ein fundamentales Maß für die arithmetische Komplexität von
Polynomen. Es verbindet Polynomtheorie, Transzendenztheorie und die Theorie
algebraischer Zahlen.

---

## Mathematischer Hintergrund

### Definition

Für ein Polynom $P(x) = a_n x^n + \cdots + a_0 \in \mathbb{C}[x]$ mit Wurzeln
$\alpha_1, \ldots, \alpha_n$ (gezählt mit Vielfachheit) ist das **Mahler-Maß**:

$$M(P) = \exp\left(\int_0^1 \log|P(e^{2\pi i t})|\, dt\right)$$

**Produktformel (Jensen-Formel):**
$$M(P) = |a_n| \cdot \prod_{|\alpha_i| > 1} |\alpha_i|$$

Das heißt: Multipliziere den Absolutwert des Leitkoeffizienten mit allen
Wurzeln außerhalb des Einheitskreises.

**Logarithmisches Mahler-Maß:**
$$m(P) = \log M(P) = \int_0^1 \log|P(e^{2\pi i t})|\, dt$$

### Wichtige Eigenschaften

| Eigenschaft | Formel |
|------------|--------|
| Multiplikativität | $M(PQ) = M(P) \cdot M(Q)$ |
| Skalierung | $M(cP) = |c| \cdot M(P)$ |
| Für Kreisteilungspolynome | $M(\Phi_n) = 1$ |
| Für $P(x) = x - \alpha$ | $M = \max(1, |\alpha|)$ |

---

## Lehrers Problem

Das **Lehmer-Dezim-Polynom** (1933):
$$L(x) = x^{10} + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1$$

hat Mahler-Maß $M(L) \approx 1.17628081825991...$

**Lehrers Problem (CONJECTURE, 1933, offen):**
> Gibt es ein ganzzahliges Polynom $P \in \mathbb{Z}[x]$, das kein
> Kreisteilungspolynom ist, mit $1 < M(P) < M(L)$?

Status: **OFFEN** – trotz intensiver rechnergestützter Suche bis Grad ~30.

---

## Kronecker-Theorem

$$M(P) = 1 \iff \text{alle Wurzeln von } P \text{ liegen auf dem Einheitskreis oder sind } 0$$

Äquivalent: $P$ ist ein Produkt aus Kreisteilungspolynomen $\Phi_n$ und $x^k$.

---

## Smyths Schranke

**Smyth (1971):** Für ganzzahlige, **nicht-reziproke** Polynome gilt:
$$M(P) \geq 1.3247179572... \quad \text{(kleinste Pisot-Zahl, Wurzel von } x^3 - x - 1\text{)}$$

Ein Polynom $P$ heißt **reziprok**, wenn $P(x) = \pm x^n P(1/x)$, d.h.
die Koeffizienten sind palindromisch oder antipalindromisch.

Das Lehmer-Polynom ist reziprok – deshalb gilt Smyths Schranke nicht für es.

---

## Salem-Zahlen

Eine **Salem-Zahl** ist eine algebraische Ganzzahl $\alpha > 1$, sodass:
- alle konjugierten Zahlen außer $\alpha$ und $1/\alpha$ auf dem Einheitskreis liegen

Salem-Zahlen sind Kandidaten für Lehrers Problem. Das Mahler-Maß des
Minimalpolynoms einer Salem-Zahl $\alpha$ ist gerade $\alpha$.

---

## Schur-Siegel-Smyth-Spurpfad

**Satz (Schur 1918, Siegel 1945):** Für eine total positive algebraische
Ganzzahl $\alpha$ vom Grad $d$ gilt:
$$\frac{\text{Spur}(\alpha)}{d} \geq 1.7719...$$

Diese Schranke ist eng mit Lehrers Problem verknüpft.

---

## Klassen

### `MahlerMeasure`

Hauptklasse für die Berechnung des Mahler-Maßes.

| Methode | Beschreibung |
|---------|-------------|
| `compute_product_formula()` | $M(P) = |a_n| \cdot \prod_{|\alpha_i|>1} |\alpha_i|$ |
| `compute_jensen_integral(n)` | Numerische Integration auf Einheitskreis |
| `logarithmic_mahler_measure()` | $m(P) = \log M(P)$ |
| `is_kronecker(tol)` | Prüft ob $M(P) = 1$ |
| `is_reciprocal()` | Prüft ob $P$ palindromisch ist |
| `smyth_bound_applies()` | Prüft Anwendbarkeit von Smyths Schranke |
| `lehmer_polynomial()` | Klassenmethode: Lehmer-Dezim-Polynom |
| `smyth_polynomial()` | Klassenmethode: $x^3 - x - 1$ |
| `cyclotomic(n)` | Klassenmethode: $\Phi_n(x)$ |

### `SchurSiegelSmythTrace`

Analyse von Spuren total positiver algebraischer Zahlen.

| Methode | Beschreibung |
|---------|-------------|
| `trace_ratio(coeffs)` | Berechnet Spur$(α)$/Grad |

---

## Verwendungsbeispiele

```python
from mahler_measure import MahlerMeasure

# Lehmer-Polynom
mm = MahlerMeasure.lehmer_polynomial()
print(mm.compute_product_formula())  # ≈ 1.17628

# Kreisteilungspolynom Φ₆(x) = x²-x+1
phi6 = MahlerMeasure.cyclotomic(6)
print(phi6.compute_product_formula())  # = 1.0
print(phi6.is_kronecker())             # True

# Smyth-Polynom x³-x-1
smyth = MahlerMeasure.smyth_polynomial()
print(smyth.compute_product_formula())  # ≈ 1.3247
print(smyth.smyth_bound_applies())      # True (nicht reziprok)

# Eigenes Polynom: 2x² - 5x + 3 = 2(x-1)(x-3/2)
mm_custom = MahlerMeasure([3, -5, 2])
print(mm_custom.compute_product_formula())  # = 2 * (3/2) = 3.0
```

---

## Offene Probleme

1. **Lehrers Problem** (CONJECTURE, offen seit 1933): Kleinste Mahler-Maß-Schranke
2. **Endlichkeit der Salem-Zahlen ≤ c**: Für jede Schranke $c > 1$ nur endlich viele?
3. **Mahler-Maß und L-Funktionen**: Zusammenhang $m(P)$ mit speziellen Werten von $L$-Funktionen (Boyd-Deninger-Konjecturen)

---

## Literatur

- Lehmer, D.H. (1933). *Factorization of certain cyclotomic functions.* Ann. Math. 34.
- Smyth, C.J. (1971). *On the product of the conjugates outside the unit circle.* Bull. LMS.
- Boyd, D.W. (1981). *Speculations concerning the range of Mahler's measure.*
- Everest, G., Ward, T. (1999). *Heights of Polynomials and Entropy in Algebraic Dynamics.*
