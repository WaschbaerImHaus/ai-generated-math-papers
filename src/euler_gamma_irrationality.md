# euler_gamma_irrationality.py – Dokumentation

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12

---

## Überblick

Dieses Modul untersucht die **Euler-Mascheroni-Konstante** $\gamma \approx 0.5772156649...$
und die damit verbundenen **Stieltjes-Konstanten** $\gamma_n$. Die Irrationalität von $\gamma$
ist eine der ältesten offenen Fragen der Zahlentheorie und steht bis heute (2026) unbewiesen.

---

## Mathematische Grundlagen

### Euler-Mascheroni-Konstante

**Definition (Grenzwert der harmonischen Reihe):**

$$\gamma = \lim_{n \to \infty} \left( \sum_{k=1}^n \frac{1}{k} - \ln n \right) \approx 0.57721566490153286...$$

**Verbindung zur Zeta-Funktion:**

$$\gamma = \lim_{s \to 1} \left( \zeta(s) - \frac{1}{s-1} \right)$$

**Verbindung zur Gamma-Funktion (Weierstraß-Produkt):**

$$\Gamma(z) = \frac{e^{-\gamma z}}{z} \prod_{n=1}^{\infty} \left(1 + \frac{z}{n}\right)^{-1} e^{z/n}$$

insbesondere: $\Gamma'(1) = \psi(1) \cdot \Gamma(1) = -\gamma$

### Reihenentwicklungen

**Vacca-Reihe (1910):**

$$\gamma = \sum_{k=1}^{\infty} (-1)^k \frac{\lfloor \log_2 k \rfloor}{k}$$

**Apéry-artige Darstellung:**

$$\gamma = 1 - \sum_{n=2}^{\infty} \frac{\zeta(n) - 1}{n}$$

Diese Reihen **beweisen nicht** die Irrationalität von $\gamma$, verdeutlichen aber seine
Verbindungen zu harmonischen Reihen.

### Irrationalitäts-Schranke (Papanikolaou 1997)

> **CONJECTURE** (kein Beweis der Irrationalität):
> Falls $\gamma = p/q$ mit $\gcd(p,q)=1$, dann gilt $q > 10^{242080}$.

Dies schließt nur sehr kleine Nenner aus. Die Irrationalität von $\gamma$ ist **offen** (Stand 2026).

---

## Klasse: EulerMascheroniConstant

### Initialisierung

```python
calc = EulerMascheroniConstant(precision=50)
```

`precision`: Anzahl der Dezimalstellen (1–100, Standard: 50).
Intern wird mpmath mit `precision + 10` Stellen konfiguriert.

### Methoden

| Methode | Beschreibung | Rückgabe |
|---------|-------------|---------|
| `compute()` | Hochpräziser Wert von $\gamma$ (Brent-McMillan) | `mpmath.mpf` |
| `compute_via_harmonic_series(n)` | Näherung via $H_n - \ln n$ | `float` |
| `compute_via_zeta_limit(ε)` | Näherung via $\zeta(1+\varepsilon) - 1/\varepsilon$ | `float` |
| `alternating_series_representation(n)` | Vacca-Reihe | `float` |
| `apery_like_representation(n)` | Formel via $\zeta(n)-1$ | `float` |
| `rational_approximation_bound()` | Papanikolaou-Schranke | `(str, int)` |
| `connection_to_gamma_function()` | Verbindung zu $\Gamma(z)$ | `dict` |
| `compute_to_n_digits(n)` | $\gamma$ als String mit $n$ Stellen | `str` |
| `verify_known_digits()` | Verifikation gegen 50 bekannte Stellen | `bool` |

### Numerische Konvergenzraten

| Methode | Fehler nach $n$ Termen |
|---------|----------------------|
| Harmonische Reihe | $O(1/n)$ – sehr langsam |
| Vacca-Reihe | $O(\log n / n)$ – langsam |
| $\zeta(n)$-Darstellung | $O(1/n \cdot 2^{-n})$ – schnell |
| mpmath (Brent-McMillan) | $O(\exp(-n))$ – optimal |

---

## Klasse: StieltjesConstants

### Definition

Die Stieltjes-Konstanten $\gamma_n$ erscheinen in der **Laurent-Entwicklung** von $\zeta(s)$ um den Pol $s=1$:

$$\zeta(s) = \frac{1}{s-1} + \sum_{n=0}^{\infty} \frac{(-1)^n \gamma_n}{n!} (s-1)^n$$

wobei $\gamma_0 = \gamma$ (Euler-Mascheroni-Konstante).

**Explizite Definition:**

$$\gamma_n = \lim_{N \to \infty} \left( \sum_{k=1}^N \frac{(\ln k)^n}{k} - \frac{(\ln N)^{n+1}}{n+1} \right)$$

### Numerische Werte

| $n$ | $\gamma_n$ (Näherung) |
|-----|----------------------|
| 0 | $+0.5772156649...$ (= $\gamma$) |
| 1 | $-0.0728158454...$ |
| 2 | $-0.0096903631...$ |
| 3 | $+0.0020538344...$ |
| 4 | $+0.0023228943...$ |

Das Vorzeichenmuster für $n \geq 2$ ist **unregelmäßig** – kein allgemeines Muster bekannt.

### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `compute(n)` | Berechnet $\gamma_n$ |
| `compute_all(max_n)` | Liste $[\gamma_0, ..., \gamma_{max_n}]$ |
| `laurent_expansion_zeta(s, terms)` | $\zeta(s)$ via Laurent-Entwicklung |
| `verify_gamma0_equals_euler()` | Prüft $\gamma_0 = \gamma$ |
| `alternating_sign_pattern(max_n)` | Vorzeichenliste |
| `size_growth(max_n)` | Absolutwerte $|\gamma_n|$ |

---

## Bekannte offene Fragen

- **Irrationalität von $\gamma$**: Offen seit Euler (1735), bis heute unbewiesen
- **Transzendenz von $\gamma$**: Noch offener als Irrationalität
- **Vorzeichenmuster der $\gamma_n$**: Kein allgemeines Muster bekannt
- **Wachstum $|\gamma_n|$**: Asymptotik bekannt, aber feinere Struktur offen

---

## Abhängigkeiten

- `mpmath`: Hochpräzise Arithmetik (Brent-McMillan-Algorithmus, `mpmath.euler`, `mpmath.stieltjes`)
- `sympy`: Symbolische Unterstützung
- `numpy`: Numerische Arrays
- `math`: Standardfunktionen
