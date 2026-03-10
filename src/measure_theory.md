# Maßtheorie – Modul-Dokumentation

**Datei:** `src/measure_theory.py`
**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10

---

## Überblick

Das Modul implementiert die fundamentalen Konzepte der **Maßtheorie**, die als
mathematisches Fundament der modernen Analysis und Wahrscheinlichkeitstheorie dient.
Die Maßtheorie verallgemeinert die Begriffe von „Länge", „Fläche" und „Volumen" auf
beliebige abstrakte Mengen.

---

## Klasse `SigmaAlgebra`

### Mathematische Definition

Eine **σ-Algebra** (sigma-Algebra) $\mathcal{F}$ über einer Grundmenge $\Omega$ ist eine
Menge von Teilmengen, die folgende drei Axiome erfüllt:

$$
\begin{aligned}
&\text{(1)} \quad \Omega \in \mathcal{F} \\
&\text{(2)} \quad A \in \mathcal{F} \Rightarrow A^c = \Omega \setminus A \in \mathcal{F} \\
&\text{(3)} \quad A_1, A_2, \ldots \in \mathcal{F} \Rightarrow \bigcup_{i=1}^{\infty} A_i \in \mathcal{F}
\end{aligned}
$$

### Wichtige Spezialfälle

| σ-Algebra | Beschreibung |
|-----------|-------------|
| $\{\emptyset, \Omega\}$ | Triviale (gröbste) σ-Algebra |
| $\mathcal{P}(\Omega)$ | Potenzmenge (feinste σ-Algebra) |
| $\mathcal{B}(\mathbb{R})$ | Borel-σ-Algebra auf $\mathbb{R}$ (erzeugt von offenen Mengen) |

### Methoden

- `is_sigma_algebra()` — Prüft alle drei Axiome
- `is_closed_under_complement()` — Axiom (2)
- `is_closed_under_countable_union()` — Axiom (3) (endliche Näherung)
- `generated_sigma_algebra(generators)` — $\sigma(\mathcal{G})$ = kleinste σ-Algebra mit $\mathcal{G} \subseteq \mathcal{F}$
- `borel_sets_finite()` — Borel = Potenzmenge auf endlichen Mengen

---

## Klasse `Measure`

### Mathematische Definition

Ein **Maß** $\mu: \mathcal{F} \to [0, \infty]$ erfüllt:

$$
\begin{aligned}
&\text{(1)} \quad \mu(\emptyset) = 0 \\
&\text{(2) σ-Additivität:} \quad A_i \cap A_j = \emptyset \Rightarrow \mu\!\left(\bigsqcup_{i=1}^\infty A_i\right) = \sum_{i=1}^\infty \mu(A_i)
\end{aligned}
$$

### Spezielle Maße

| Maß | Definition |
|-----|-----------|
| Wahrscheinlichkeitsmaß | $\mu(\Omega) = 1$ |
| σ-endliches Maß | $\Omega = \bigcup_n A_n$ mit $\mu(A_n) < \infty$ |
| Zählmaß | $\mu(A) = |A|$ |

### Äußeres Maß

Das **äußere Maß** ist das Infimum aller Überdeckungen:

$$\mu^*(A) = \inf\!\left\{\sum_i \mu(A_i) \,\Big|\, A \subseteq \bigcup_i A_i,\; A_i \in \mathcal{F}\right\}$$

---

## Klasse `LebesgueMeasure`

Das **Lebesgue-Maß** $\lambda$ auf $\mathbb{R}$ ist das eindeutige translationsinvariante Maß mit:

$$\lambda([a, b]) = b - a$$

### Cantor-Menge

Die **Cantor-Menge** entsteht durch iteratives Entfernen des mittleren Drittels:

$$\lambda(\text{Cantor}) = \lim_{n \to \infty} \left(\frac{2}{3}\right)^n = 0$$

Obwohl die Cantor-Menge überabzählbar viele Elemente enthält, hat sie Maß **null**.

### Fette Cantor-Menge (Smith-Volterra-Cantor)

Im $n$-ten Schritt wird ein Intervall der Länge $\varepsilon / 4^n$ aus jedem verbleibenden Intervall entfernt:

$$\lambda(\text{Fette Cantor-Menge}) = 1 - \varepsilon$$

---

## Funktion `lebesgue_integral`

Das **Lebesgue-Integral** $\int_a^b f \, d\lambda$ wird numerisch durch Approximation auf feinem Gitter berechnet.

**Schlüsselprinzip:** Beim Lebesgue-Integral wird der **Wertebereich** partitioniert (nicht der Definitionsbereich wie beim Riemann-Integral):

$$\int f \, d\mu = \sum_{i} y_i \cdot \mu\!\left(\{x : f(x) \approx y_i\}\right)$$

Für stetige Funktionen stimmt dies mit dem Riemann-Integral überein.

---

## Funktion `riemann_vs_lebesgue`

Vergleich der Integrationsbegriffe:

| Funktion | Riemann-int.? | Lebesgue-int.? | Integralwert |
|----------|:-------------:|:--------------:|:------------:|
| Dirichlet $\mathbf{1}_\mathbb{Q}$ | Nein | Ja | 0 |
| $\sin(1/x)$ | Ja | Ja | ≈ -0.504 |
| $1/\sqrt{x}$ | Uneigentlich | Ja | 2 |

**Dirichlet-Funktion:**
$$f(x) = \begin{cases} 1 & x \in \mathbb{Q} \\ 0 & x \notin \mathbb{Q} \end{cases}$$
$\Rightarrow$ Riemann-Ober- und -Untersumme weichen auseinander (1 ≠ 0), aber $\lambda(\mathbb{Q}) = 0$, also $\int f \, d\lambda = 0$.

---

## Klasse `MeasurableFunction`

Eine Funktion $f: (\Omega, \mathcal{F}) \to (\mathbb{R}, \mathcal{B})$ ist **messbar**, wenn:

$$f^{-1}(B) \in \mathcal{F} \quad \text{für alle Borel-Mengen } B \subseteq \mathbb{R}$$

**Integral einer Treppenfunktion:**
$$\int f \, d\mu = \sum_{i} a_i \cdot \mu(A_i) \quad \text{für } f = \sum_i a_i \cdot \mathbf{1}_{A_i}$$

**$L^p$-Norm:**
$$\|f\|_p = \left(\int |f|^p \, d\mu\right)^{1/p}$$

---

## Konvergenzsätze

### Satz von der monotonen Konvergenz (Beppo Levi)

$$f_n \nearrow f \text{ messbar}, f_n \geq 0 \quad \Rightarrow \quad \lim_{n \to \infty} \int f_n \, d\mu = \int f \, d\mu$$

**Beispiel:** $f_n(x) = (1 - 1/n) \cdot x^2 \nearrow x^2$, $\int_0^1 f_n \, d\lambda = (1-1/n)/3 \to 1/3$.

### Satz von der dominierten Konvergenz (Lebesgue)

$$|f_n| \leq g \text{ (integrierbar)}, f_n \to f \text{ (p.ü.)} \quad \Rightarrow \quad \lim_{n \to \infty} \int f_n \, d\mu = \int f \, d\mu$$

**Beispiel:** $f_n(x) = \sin(nx)/n \to 0$, dominiert durch $g(x) = 1$.

### Satz von Fubini

$$\int_{\Omega_1 \times \Omega_2} f \, d(\mu_1 \otimes \mu_2) = \int_{\Omega_1} \left(\int_{\Omega_2} f(x,y) \, d\mu_2(y)\right) d\mu_1(x)$$

**Verifikation:** $\iint_{[0,1]^2} (x^2 + y^2) \, d(x, y) = 2/3$.

---

## Klasse `LpSpace`

**$L^p$-Räume** sind vollständige normierte Räume (Banach-Räume) für $p \geq 1$:

$$L^p(\mu) = \left\{f : \|f\|_p = \left(\int |f|^p \, d\mu\right)^{1/p} < \infty\right\}$$

Für $p = 2$ ist $L^2$ ein **Hilbert-Raum** mit Skalarprodukt:

$$\langle f, g \rangle = \int f \cdot g \, d\mu$$

### Hölder-Ungleichung

Für $\frac{1}{p} + \frac{1}{q} = 1$ gilt:

$$\int |f \cdot g| \, d\mu \leq \|f\|_p \cdot \|g\|_q$$

### Minkowski-Ungleichung (Dreiecksungleichung)

$$\|f + g\|_p \leq \|f\|_p + \|g\|_p$$

---

## Funktion `radon_nikodym_theorem_demo`

**Satz von Radon-Nikodym:**
Wenn $\nu \ll \mu$ (ν ist absolut stetig bezüglich μ), d.h.:

$$\mu(A) = 0 \Rightarrow \nu(A) = 0$$

dann existiert eine messbare Funktion $f \geq 0$ (Radon-Nikodym-Ableitung) mit:

$$\nu(A) = \int_A f \, d\mu \quad \text{für alle } A \in \mathcal{F}$$

Schreibweise: $f = \dfrac{d\nu}{d\mu}$ (verallgemeinerte „Ableitung" zweier Maße).

---

## Funktion `cantor_function_demo`

Die **Cantor-Funktion** (Teufelsleiter) $c: [0,1] \to [0,1]$ ist:

- **Stetig** und monoton wachsend
- **Fast überall** Ableitung 0: $c'(x) = 0$ für $\lambda$-fast alle $x$
- **Nicht absolut stetig**: kein $f \in L^1$ mit $c(x) = c(0) + \int_0^x f(t) \, dt$
- **Singulär**: Das Lebesgue-Stieltjes-Maß $\mu_c$ ist singulär, $\mu_c \perp \lambda$

**Konstruktion via ternärer Entwicklung:**
- Schreibe $x$ in Basis 3
- Bis zur ersten „1": ersetze „2" durch „1"
- An der ersten „1": restliche Stellen → 0
- Interpretiere als Binärzahl

$$c(1/3) = c(2/3) = 1/2, \quad c(0) = 0, \quad c(1) = 1$$
