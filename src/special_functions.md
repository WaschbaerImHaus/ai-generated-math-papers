# Spezielle Funktionen – Modul-Dokumentation

**Datei:** `src/special_functions.py`
**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10

---

## Überblick

Das Modul implementiert die wichtigsten **speziellen Funktionen** der mathematischen Physik
und Analysis. Diese Funktionen treten als Lösungen klassischer Differentialgleichungen auf
und bilden vollständige Orthogonalsysteme.

---

## Klasse `BesselFunctions`

### Bessel-Differentialgleichung

Die Bessel-Funktionen sind Lösungen der **Bessel-DGL** der Ordnung $\nu$:

$$x^2 y'' + x y' + (x^2 - \nu^2) y = 0$$

### Funktionstypen

| Funktion | Symbol | Eigenschaft |
|----------|--------|-------------|
| Bessel 1. Art | $J_\nu(x)$ | Endlich bei $x = 0$ für $\nu \geq 0$ |
| Bessel 2. Art (Neumann) | $Y_\nu(x)$ | Singulär bei $x = 0$ |
| Modifiziert 1. Art | $I_\nu(x)$ | Exponentiell wachsend |
| Modifiziert 2. Art | $K_\nu(x)$ | Exponentiell abklingend |

### Rekurrenzrelation

$$J_{\nu-1}(x) + J_{\nu+1}(x) = \frac{2\nu}{x} J_\nu(x)$$

### Integraldarstellung (ganzzahlige Ordnung)

$$J_n(x) = \frac{1}{\pi} \int_0^\pi \cos(n\tau - x \sin\tau) \, d\tau$$

### Asymptotik für große $x$

$$J_\nu(x) \approx \sqrt{\frac{2}{\pi x}} \cos\!\left(x - \frac{\nu\pi}{2} - \frac{\pi}{4}\right), \quad x \to \infty$$

---

## Klasse `LegendrePolynomials`

### Legendre-Differentialgleichung

$$\left(1 - x^2\right) y'' - 2x y' + n(n+1) y = 0, \quad x \in [-1, 1]$$

### Rodrigues-Formel

$$P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} \left(x^2 - 1\right)^n$$

### Niedrige Grade

$$P_0(x) = 1, \quad P_1(x) = x, \quad P_2(x) = \frac{3x^2 - 1}{2}, \quad P_3(x) = \frac{5x^3 - 3x}{2}$$

### Orthogonalität

$$\int_{-1}^{1} P_n(x) \, P_m(x) \, dx = \frac{2}{2n+1} \delta_{nm}$$

### Erzeugende Funktion

$$\frac{1}{\sqrt{1 - 2xt + t^2}} = \sum_{n=0}^{\infty} P_n(x) \, t^n, \quad |t| < 1$$

### Zugeordnete Legendre-Polynome

$$P_n^m(x) = (-1)^m \left(1 - x^2\right)^{m/2} \frac{d^m}{dx^m} P_n(x)$$

### Kugelflächenfunktionen

$$Y_n^m(\theta, \phi) \propto P_n^{|m|}(\cos\theta) \cdot \begin{cases} \cos(m\phi) & m \geq 0 \\ \sin(|m|\phi) & m < 0 \end{cases}$$

---

## Klasse `HypergeometricFunctions`

### Definition

Die **Gaußsche hypergeometrische Funktion** ist definiert durch die Potenzreihe:

$$\,_2F_1(a, b; c; z) = \sum_{n=0}^{\infty} \frac{(a)_n (b)_n}{(c)_n \, n!} z^n, \quad |z| < 1$$

wobei $(a)_n = a(a+1) \cdots (a+n-1)$ das **Pochhammer-Symbol** ist.

### Gaußsche hypergeometrische DGL

$$z(1-z) w'' + [c - (a+b+1)z] w' - ab \, w = 0$$

### Wichtige Sonderfälle

| Funktion | Hypergeometrische Darstellung |
|----------|------------------------------|
| $(1-z)^{-a}$ | $\,_1F_0(a; -; z)$ |
| $\ln(1+z)/z$ | $\,_2F_1(1, 1; 2; -z)$ |
| $P_n(x)$ | $\,_2F_1(-n, n+1; 1; \frac{1-x}{2})$ |
| $K(k)$ | $\frac{\pi}{2} \,_2F_1(\frac{1}{2}, \frac{1}{2}; 1; k^2)$ |
| $e^z$ | $\,_1F_1(a; a; z)$ |

---

## Klasse `AiryFunctions`

### Airy-Differentialgleichung

$$y'' = x \, y$$

Die beiden linear unabhängigen Lösungen sind $\text{Ai}(x)$ und $\text{Bi}(x)$.

### Eigenschaften

| Funktion | $x \to +\infty$ | $x \to -\infty$ |
|----------|:---------------:|:---------------:|
| $\text{Ai}(x)$ | $\to 0$ (abklingend) | Oszillierend |
| $\text{Bi}(x)$ | $\to \infty$ (wachsend) | Oszillierend |

### Asymptotik für großes positives $x$

$$\text{Ai}(x) \approx \frac{e^{-\zeta}}{2\sqrt{\pi} \, x^{1/4}}, \quad \zeta = \frac{2}{3} x^{3/2}$$

### Nullstellen von $\text{Ai}(x)$

Alle Nullstellen liegen auf der **negativen reellen Achse**. Asymptotisch:

$$a_n \approx -\left[\frac{3\pi(4n-1)}{8}\right]^{2/3}$$

### Physikalische Bedeutung

Airy-Funktionen treten in der **WKB-Näherung** der Quantenmechanik auf (lineare Potentiale)
und in der geometrischen Optik (Kaustiken, Beugung).

---

## Klasse `OrthogonalPolynomials`

Alle klassischen orthogonalen Polynome erfüllen eine **Drei-Term-Rekurrenz**:

$$p_{n+1}(x) = (\alpha_n x + \beta_n) p_n(x) - \gamma_n p_{n-1}(x)$$

### Übersichtstabelle

| Polynom | Symbol | Gewicht $w(x)$ | Intervall |
|---------|--------|:--------------:|:---------:|
| Hermite | $H_n(x)$ | $e^{-x^2}$ | $(-\infty, \infty)$ |
| Laguerre | $L_n(x)$ | $e^{-x}$ | $[0, \infty)$ |
| Chebyshev I | $T_n(x)$ | $1/\sqrt{1-x^2}$ | $[-1, 1]$ |
| Chebyshev II | $U_n(x)$ | $\sqrt{1-x^2}$ | $[-1, 1]$ |
| Legendre | $P_n(x)$ | $1$ | $[-1, 1]$ |
| Gegenbauer | $C_n^\alpha(x)$ | $(1-x^2)^{\alpha-1/2}$ | $[-1, 1]$ |

### Hermite-Polynome

$$H_0 = 1, \quad H_1 = 2x, \quad H_2 = 4x^2 - 2, \quad H_{n+1} = 2x H_n - 2n H_{n-1}$$

$$\int_{-\infty}^{\infty} H_n(x) H_m(x) e^{-x^2} dx = 2^n n! \sqrt{\pi} \, \delta_{nm}$$

### Chebyshev-Polynome

$$T_n(\cos\theta) = \cos(n\theta), \quad U_n(\cos\theta) = \frac{\sin((n+1)\theta)}{\sin\theta}$$

---

## Funktion `gamma_function_properties`

Die **Gamma-Funktion** $\Gamma: \mathbb{C} \setminus \{0, -1, -2, \ldots\} \to \mathbb{C}$:

$$\Gamma(z) = \int_0^\infty t^{z-1} e^{-t} \, dt, \quad \operatorname{Re}(z) > 0$$

**Funktionalgleichung:** $\Gamma(z+1) = z \cdot \Gamma(z)$

**Reflexionsformel (Euler):**
$$\Gamma(z) \cdot \Gamma(1-z) = \frac{\pi}{\sin(\pi z)}$$

**Duplikationsformel (Legendre):**
$$\Gamma(z) \cdot \Gamma\!\left(z + \tfrac{1}{2}\right) = \frac{\sqrt{\pi}}{2^{2z-1}} \Gamma(2z)$$

**Besondere Werte:** $\Gamma(1) = 1$, $\Gamma(1/2) = \sqrt{\pi}$, $\Gamma(n+1) = n!$

---

## Funktion `beta_function`

$$B(a, b) = \frac{\Gamma(a)\,\Gamma(b)}{\Gamma(a+b)} = \int_0^1 t^{a-1}(1-t)^{b-1} \, dt$$

**Symmetrie:** $B(a, b) = B(b, a)$

**Besondere Werte:** $B(1/2, 1/2) = \pi$, $B(1, 1) = 1$

---

## Funktion `digamma_function`

$$\psi(x) = \frac{d}{dx} \ln\Gamma(x) = \frac{\Gamma'(x)}{\Gamma(x)}$$

**Rekurrenz:** $\psi(x+1) = \psi(x) + \dfrac{1}{x}$

**Besondere Werte:** $\psi(1) = -\gamma$ (Euler-Mascheroni), $\psi(n+1) = H_n - \gamma$

---

## Funktion `riemann_zeta_special_values`

Die Riemann-Zeta-Funktion an geraden positiven Stellen:

$$\zeta(2) = \frac{\pi^2}{6}, \quad \zeta(4) = \frac{\pi^4}{90}, \quad \zeta(6) = \frac{\pi^6}{945}$$

**Bernoulli-Zusammenhang:**
$$\zeta(2n) = \frac{(-1)^{n+1} (2\pi)^{2n} B_{2n}}{2 \cdot (2n)!}$$

**Negative Argumente:** $\zeta(-1) = -\frac{1}{12}$ (Ramanujan), $\zeta(0) = -\frac{1}{2}$

---

## Funktion `elliptic_integrals`

**Vollständige elliptische Integrale:**

$$K(k) = \int_0^{\pi/2} \frac{d\theta}{\sqrt{1 - k^2 \sin^2\theta}}, \quad E(k) = \int_0^{\pi/2} \sqrt{1 - k^2 \sin^2\theta} \, d\theta$$

**Legendre-Relation:** $K(k) E(k') + K(k') E(k) - K(k) K(k') = \dfrac{\pi}{2}$

mit komplementärem Modulus $k' = \sqrt{1-k^2}$.

**Anwendung (Pendel):** Periode eines Pendels der Länge $L$ mit Amplitude $\theta_0$:

$$T = 4\sqrt{\frac{L}{g}} \, K\!\left(\sin\frac{\theta_0}{2}\right)$$

---

## Funktion `error_function_properties`

$$\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \, dt$$

**Zusammenhang mit Normalverteilung $\mathcal{N}(0,1)$:**

$$P(|X| \leq x) = \operatorname{erf}\!\left(\frac{x}{\sqrt{2}}\right)$$

| Sigma | Wahrscheinlichkeit |
|:-----:|:-----------------:|
| $1\sigma$ | $\approx 68.27\%$ |
| $2\sigma$ | $\approx 95.45\%$ |
| $3\sigma$ | $\approx 99.73\%$ |

**Asymptotik:** $\operatorname{erfc}(x) \approx \dfrac{e^{-x^2}}{x\sqrt{\pi}}$ für $x \to \infty$
