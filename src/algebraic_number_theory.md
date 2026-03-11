# algebraic_number_theory.py — Dokumentation

**Autor:** Michael Fuhrmann
**Letzte Änderung:** 2026-03-11
**Build:** 84

---

## Überblick

Das Modul `algebraic_number_theory.py` implementiert die grundlegenden Konzepte der
**algebraischen Zahlentheorie** — dem Studium von Zahlkörpern und ihren ganzen Zahlen.

Es umfasst sieben Hauptbereiche:

1. [Zahlkörper-Grundlagen](#1-zahlkörper-grundlagen)
2. [Ideal-Arithmetik (Dedekind)](#2-ideal-arithmetik-dedekind)
3. [Klassengruppe](#3-klassengruppe)
4. [Einheitengruppe (Dirichlet)](#4-einheitengruppe-dirichlet)
5. [Quadratische Körper](#5-quadratische-körper)
6. [Zyklotomische Körper](#6-zyklotomische-körper)
7. [Lokale Theorie](#7-lokale-theorie)

---

## 1. Zahlkörper-Grundlagen

### Mathematischer Hintergrund

Ein **Zahlkörper** $K$ ist eine endliche Körpererweiterung von $\mathbb{Q}$:

$$[K : \mathbb{Q}] = n < \infty$$

Jeder Zahlkörper entsteht durch Adjunktion einer algebraischen Zahl $\alpha$:

$$K = \mathbb{Q}(\alpha) \cong \mathbb{Q}[x]/(f(x))$$

wobei $f \in \mathbb{Q}[x]$ das **Minimalpolynom** von $\alpha$ ist (irreduzibel, normiert).

### Klasse `NumberField`

```python
nf = NumberField([-2, 0, 1])   # x² - 2 = 0, erzeugt Q(√2)
```

**Koeffizientenkonvention:** Aufsteigend — `[a₀, a₁, ..., aₙ]` entspricht
$f(x) = a_0 + a_1 x + \ldots + a_n x^n$

#### `degree()` — Grad $[K:\mathbb{Q}]$

$$n = \deg(f)$$

#### `signature()` — Signatur $(r_1, r_2)$

Die Signatur beschreibt die Einbettungen von $K$ in $\mathbb{C}$:

- $r_1$ = Anzahl reeller Einbettungen $K \hookrightarrow \mathbb{R}$
- $r_2$ = Anzahl konjugierter Paare $K \hookrightarrow \mathbb{C}$ (nicht reell)
- **Beziehung:** $r_1 + 2r_2 = n$

| Körpertyp | $r_1$ | $r_2$ | Beispiel |
|-----------|--------|--------|---------|
| Reell-quadratisch ($d > 0$) | 2 | 0 | $\mathbb{Q}(\sqrt{2})$ |
| Imaginär-quadratisch ($d < 0$) | 0 | 1 | $\mathbb{Q}(\sqrt{-1})$ |
| Totalreell (Grad $n$) | $n$ | 0 | $\mathbb{Q}(\sqrt[n]{2})$ für bestimmte $n$ |

#### `ring_of_integers()` — Ganzheitsring $\mathcal{O}_K$

Der **Ring der ganzen Zahlen** $\mathcal{O}_K$ besteht aus allen $\alpha \in K$,
die ein normiertes Polynom mit ganzzahligen Koeffizienten erfüllen.

Für quadratische Körper $K = \mathbb{Q}(\sqrt{d})$ ($d$ quadratfrei):

$$\mathcal{O}_K = \begin{cases}
\mathbb{Z}\left[\frac{1+\sqrt{d}}{2}\right] & \text{falls } d \equiv 1 \pmod{4} \\
\mathbb{Z}[\sqrt{d}] & \text{falls } d \equiv 2,3 \pmod{4}
\end{cases}$$

#### `discriminant()` — Diskriminante $\Delta(K/\mathbb{Q})$

$$\Delta(K/\mathbb{Q}) = (-1)^{n(n-1)/2} \cdot N_{K/\mathbb{Q}}(f'(\alpha))$$

Für quadratische Körper:

$$\Delta = \begin{cases}
d & \text{falls } d \equiv 1 \pmod{4} \\
4d & \text{sonst}
\end{cases}$$

**Bedeutung:** $p \mid \Delta \iff p$ verzweigt in $K$.

---

## 2. Ideal-Arithmetik (Dedekind)

### Dedekind-Ringe

$\mathcal{O}_K$ ist ein **Dedekind-Ring**: noetherisch, ganzabgeschlossen, und
jedes Primideal $\neq 0$ ist maximal. Haupteigenschaft:

> **Eindeutige Primidealzerlegung:** Jedes Ideal $I \neq 0$ in $\mathcal{O}_K$ hat eine
> eindeutige Darstellung $I = \mathfrak{p}_1^{e_1} \cdots \mathfrak{p}_r^{e_r}$.

### Klasse `DedekindIdeal`

```python
ideal = DedekindIdeal([2, 3], field_degree=2)   # erzeugt von 2 und 3
```

#### `norm()` — Norm $N(I)$

$$N(I) = |\mathcal{O}_K/I|$$

Für das Hauptideal $(a)$ mit $a \in \mathbb{Z}$: $N((a)) = a^{[K:\mathbb{Q}]}$

#### `is_prime_ideal()` — Primalitätstest

Ein Ideal $\mathfrak{p}$ ist prim genau dann, wenn $\mathcal{O}_K/\mathfrak{p}$
ein Integritätsring ist.

### Funktion `ideal_factorization(p, poly_coeffs)`

Zerlegung einer rationalen Primzahl $p$ in $\mathcal{O}_K$:

$$p\mathcal{O}_K = \mathfrak{p}_1^{e_1} \cdot \mathfrak{p}_2^{e_2} \cdots \mathfrak{p}_r^{e_r}$$

Für quadratische Körper $\mathbb{Q}(\sqrt{d})$ und $p$ ungerade:

| Legendre-Symbol $(d/p)$ | Zerlegungstyp | Beschreibung |
|-------------------------|---------------|-------------|
| $-1$ | **Inert** | $p$ bleibt prim, $N(\mathfrak{p}) = p^2$ |
| $0$ | **Verzweigt** | $(p) = \mathfrak{p}^2$, $p \mid \Delta$ |
| $1$ | **Zerfällt** | $(p) = \mathfrak{p}\overline{\mathfrak{p}}$, $N(\mathfrak{p}) = p$ |

### Fundamentalsatz: `efg_equation_check(p, poly_coeffs)`

$$\sum_{i=1}^{r} e_i \cdot f_i = [K:\mathbb{Q}] = n$$

wobei:
- $e_i$ = Verzweigungsindex von $\mathfrak{p}_i$ über $p$
- $f_i$ = Trägheitsgrad: $[\mathcal{O}_K/\mathfrak{p}_i : \mathbb{F}_p]$
- $r$ = Anzahl der Primideale über $p$

---

## 3. Klassengruppe

### Mathematischer Hintergrund

Die **Klassengruppe** $\text{Cl}(K)$ misst, wie weit $\mathcal{O}_K$ von einem
Hauptidealbereich entfernt ist:

$$\text{Cl}(K) = \{\text{gebrochene Ideale}\}/\{\text{Hauptideale}\}$$

Die **Klassenzahl** $h(K) = |\text{Cl}(K)|$ ist endlich.

### `class_number(d)` — Klassenzahl von $\mathbb{Q}(\sqrt{d})$

Bekannte Werte (imaginär-quadratisch):

| $d$ | $h$ | Bemerkung |
|-----|-----|-----------|
| $-1$ | 1 | $\mathbb{Z}[i]$ = Gaußsche ganze Zahlen (PID) |
| $-3$ | 1 | $\mathbb{Z}[\omega]$ = Eisensteinsche Zahlen (PID) |
| $-5$ | 2 | erste nicht-triviale Klasse |
| $-23$ | 3 | |
| $-163$ | 1 | größter negativer Stark-Heegner-Körper |

**Satz (Stark-Heegner):** Die einzigen $d < 0$ mit $h(\mathbb{Q}(\sqrt{d})) = 1$ sind:
$$d \in \{-1, -2, -3, -7, -11, -19, -43, -67, -163\}$$

### `minkowski_bound(d)` — Minkowski-Schranke

$$M_K = \left(\frac{n!}{n^n}\right) \left(\frac{4}{\pi}\right)^{r_2} \sqrt{|\Delta|}$$

Für quadratische Körper:

$$M_K = \begin{cases}
\frac{2}{\pi}\sqrt{|\Delta|} & \text{imaginär-quadratisch} \\
\frac{1}{2}\sqrt{|\Delta|} & \text{reell-quadratisch}
\end{cases}$$

**Bedeutung:** Jede Idealklasse enthält ein Ideal mit $N(I) \leq M_K$.
Man muss nur Primideale mit Norm $\leq M_K$ untersuchen.

### `is_pid(d)` — Hauptidealbereichtest

$$\mathcal{O}_K \text{ ist PID} \iff h(K) = 1$$

---

## 4. Einheitengruppe (Dirichlet)

### Dirichletscher Einheitensatz

$$\mathcal{O}_K^\times \cong \mu(K) \times \mathbb{Z}^r$$

wobei:
- $\mu(K)$ = endliche Gruppe der Einheitswurzeln in $K$
- $r = r_1 + r_2 - 1$ = **Einheitenrang**

### `unit_rank(poly_coeffs)` — Einheitenrang

| Körpertyp | $r_1$ | $r_2$ | $r = r_1 + r_2 - 1$ |
|-----------|--------|--------|----------------------|
| $\mathbb{Q}(\sqrt{d})$, $d > 0$ | 2 | 0 | 1 |
| $\mathbb{Q}(\sqrt{d})$, $d < 0$ | 0 | 1 | 0 |
| $\mathbb{Q}(\zeta_p)$, $p$ prim | 0 | $(p-1)/2$ | $(p-3)/2$ |

### `fundamental_units_quadratic(d)` — Fundamentaleinheit

Für reell-quadratische Körper ($r = 1$) gibt es eine **Fundamentaleinheit** $\varepsilon$:

$$\mathcal{O}_{\mathbb{Q}(\sqrt{d})}^\times = \{\pm \varepsilon^n : n \in \mathbb{Z}\}$$

$\varepsilon$ löst die **Pell-Gleichung:**

$$x^2 - d y^2 = \pm 1$$

Bekannte Fundamentaleinheiten:

| $d$ | $\varepsilon$ | $N(\varepsilon)$ |
|-----|--------------|-----------------|
| 2 | $1 + \sqrt{2}$ | $-1$ |
| 3 | $2 + \sqrt{3}$ | $1$ |
| 5 | $(1+\sqrt{5})/2$ | $-1$ |
| 7 | $8 + 3\sqrt{7}$ | $1$ |

### `regulator_estimate(d)` — Regulator

$$R_K = \log|\varepsilon| \quad \text{(für reell-quadratische Körper)}$$

Der Regulator erscheint in der Dirichletschen Klassenzahlformel:

$$\lim_{s \to 1} (s-1) \zeta_K(s) = \frac{2^{r_1}(2\pi)^{r_2} h_K R_K}{w_K \sqrt{|\Delta|}}$$

---

## 5. Quadratische Körper

### Klasse `QuadraticField`

```python
K = QuadraticField(-5)   # Q(√-5)
K.discriminant()         # -20
K.signature()            # (0, 1)
K.class_number()         # 2
K.ring_of_integers_basis()  # ["1", "√-5"]
```

### `legendre_symbol_fast(a, p)` — Legendre-Symbol $(a/p)$

**Euler-Kriterium:**
$$\left(\frac{a}{p}\right) \equiv a^{(p-1)/2} \pmod{p}$$

$$\left(\frac{a}{p}\right) = \begin{cases}
0 & \text{falls } p \mid a \\
1 & \text{falls } a \text{ QR mod } p \\
-1 & \text{falls } a \text{ NR mod } p
\end{cases}$$

### `kronecker_symbol(a, n)` — Kronecker-Symbol

Verallgemeinerung des Jacobi-Symbols auf beliebige $n$:

$$\left(\frac{a}{n}\right) = \left(\frac{a}{-1}\right) \cdot \left(\frac{a}{2}\right)^{v_2(n)} \cdot \prod_{p \mid n, p \text{ odd}} \left(\frac{a}{p}\right)^{v_p(n)}$$

Für die Primzahl 2 gelten Spezialregeln:

$$\left(\frac{a}{2}\right) = \begin{cases}
0 & \text{falls } a \equiv 0 \pmod{2} \\
1 & \text{falls } a \equiv \pm 1 \pmod{8} \\
-1 & \text{falls } a \equiv \pm 3 \pmod{8}
\end{cases}$$

### `quadratic_reciprocity_law(p, q)` — Quadratisches Reziprozitätsgesetz

**Gaußsches QRG (1796):** Für verschiedene ungerade Primzahlen $p, q$:

$$\left(\frac{p}{q}\right)\left(\frac{q}{p}\right) = (-1)^{\frac{p-1}{2} \cdot \frac{q-1}{2}}$$

**Ergänzungssätze:**

$$\left(\frac{-1}{p}\right) = (-1)^{(p-1)/2} = \begin{cases}
1 & p \equiv 1 \pmod{4} \\
-1 & p \equiv 3 \pmod{4}
\end{cases}$$

$$\left(\frac{2}{p}\right) = (-1)^{(p^2-1)/8} = \begin{cases}
1 & p \equiv \pm 1 \pmod{8} \\
-1 & p \equiv \pm 3 \pmod{8}
\end{cases}$$

---

## 6. Zyklotomische Körper

### Klasse `CyclotomicField`

```python
cf = CyclotomicField(5)   # Q(ζ₅)
cf.degree()               # φ(5) = 4
cf.signature()            # (0, 2)
```

### Mathematischer Hintergrund

Der $n$-te zyklotomische Körper $\mathbb{Q}(\zeta_n)$ mit $\zeta_n = e^{2\pi i/n}$:

$$[\mathbb{Q}(\zeta_n) : \mathbb{Q}] = \varphi(n)$$

**Galoisgruppe:**

$$\text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \cong (\mathbb{Z}/n\mathbb{Z})^\times$$

durch $\sigma_a(\zeta_n) = \zeta_n^a$ für $\gcd(a,n) = 1$.

**Ganzheitsring:**

$$\mathcal{O}_{\mathbb{Q}(\zeta_n)} = \mathbb{Z}[\zeta_n]$$

$\mathbb{Z}$-Basis: $\{1, \zeta_n, \zeta_n^2, \ldots, \zeta_n^{\varphi(n)-1}\}$

### `splitting_of_prime(p)` — Zerlegung von $p$ in $\mathbb{Z}[\zeta_n]$

- **$p \mid n$:** Vollständig verzweigt: $(p) = \mathfrak{p}^{\varphi(n)}$
- **$p \nmid n$:** Sei $f = \text{ord}_n(p)$ (Ordnung von $p$ mod $n$).
  Dann: $(p) = \mathfrak{p}_1 \cdots \mathfrak{p}_r$ mit $r = \varphi(n)/f$, $e = 1$

Fundamentalsatz: $e \cdot f \cdot r = \varphi(n)$

### `kummer_lifting(p, n)` — Kummer-Theorie

Kummers Strategie für Fermats letzten Satz in $\mathbb{Z}[\zeta_p]$:

$$x^p + y^p = z^p \implies (x+y)(x+\zeta_p y)(x+\zeta_p^2 y)\cdots(x+\zeta_p^{p-1}y) = z^p$$

**Reguläre Primzahlen** ($p \nmid h(\mathbb{Q}(\zeta_p))$): $p < 37$ außer der leeren Menge.
**Irreguläre Primzahlen:** $37, 59, 67, 101, \ldots$

---

## 7. Lokale Theorie

### `p_adic_completion(poly_coeffs, p)` — $p$-adische Vervollständigung

**Produktformel:**
$$K \otimes_\mathbb{Q} \mathbb{Q}_p \cong \prod_{i=1}^r K_{\mathfrak{p}_i}$$

Jeder lokale Faktor $K_{\mathfrak{p}}$ ist ein lokaler Körper mit:
- Trägheitsgrad $f_i$ über $\mathbb{Q}_p$
- Verzweigungsindex $e_i$
- Residuenkörper $\mathbb{F}_{p^{f_i}}$

### `hensel_lifting_general(poly_coeffs, a0, p, steps)` — Hensel-Lifting

**Hensels Lemma:** Sei $f \in \mathbb{Z}_p[x]$, $a_0 \in \mathbb{Z}$ mit:
- $f(a_0) \equiv 0 \pmod{p}$
- $f'(a_0) \not\equiv 0 \pmod{p}$

Dann gibt es eindeutig $\hat{a} \in \mathbb{Z}_p$ mit $f(\hat{a}) = 0$ und $\hat{a} \equiv a_0 \pmod{p}$.

**Newton-Iteration** (quadratische Konvergenz):

$$a_{n+1} = a_n - \frac{f(a_n)}{f'(a_n)} \pmod{p^{2^n}}$$

### `local_norm_symbol(a, b, p)` — Hilbert-Symbol $(a,b)_p$

Das **Hilbert-Symbol** $(a,b)_p \in \{+1,-1\}$ ist definiert:

$$(a,b)_p = 1 \iff ax^2 + by^2 = z^2 \text{ hat eine nichttriviale Lösung in } \mathbb{Q}_p$$

**Eigenschaften:**
- Symmetrisch: $(a,b)_p = (b,a)_p$
- Bimultiplikativ: $(a,bc)_p = (a,b)_p \cdot (a,c)_p$
- **Produktformel:** $\prod_v (a,b)_v = 1$ (über alle Stellen $v$)

Für ungerades $p$ und $a, b \in \mathbb{Z}_p^\times$:

$$\left(\frac{a}{p}\right)^{v_p(b)} \cdot \left(\frac{b}{p}\right)^{v_p(a)} \cdot (-1)^{v_p(a) \cdot v_p(b) \cdot \frac{p-1}{2}}$$

---

## 8. Klassenzahlformel

### `dirichlet_class_number_formula(d)` — Dirichletsche Klassenzahlformel

**Für imaginär-quadratische Körper** ($d < 0$):

$$h(K) = \frac{w_K \sqrt{|\Delta|}}{2\pi} \cdot L(1, \chi_\Delta)$$

wobei:
- $w_K$ = Anzahl der Einheitswurzeln
- $L(1,\chi_\Delta) = \sum_{n=1}^\infty \chi_\Delta(n)/n$ mit $\chi_\Delta = (\Delta/\cdot)$ (Kronecker-Symbol)

**Für reell-quadratische Körper** ($d > 0$):

$$h(K) \cdot R_K = \frac{\sqrt{\Delta}}{2} \cdot L(1, \chi_\Delta)$$

Diese Formel verbindet die analytische ($L$-Funktion) mit der algebraischen
(Klassenzahl, Regulator) Seite der Zahlentheorie — ein zentrales Resultat der Arithmetik!

---

## Abhängigkeiten

| Bibliothek | Verwendung |
|------------|------------|
| `math` | Standardfunktionen (sqrt, log, pi) |
| `cmath` | Komplexe Nullstellen |
| `fractions` | Rationale Arithmetik (optional) |
| `functools` | `@lru_cache` für Memoization |
| `numpy` (optional) | Numerische Nullstellenberechnung via `np.roots` |
| `sympy` (optional) | Symbolische Berechnungen |

---

## Testabdeckung

Datei: `tests/test_algebraic_number_theory.py`

| Testklasse | Anzahl Tests | Abdeckung |
|------------|-------------|-----------|
| `TestHilfsfunktionen` | 15 | Quadratfrei, Euler-φ, Legendre, Pell |
| `TestNumberField` | 8 | Grad, Signatur, OI-Basis, Diskriminante |
| `TestIdealArithmetik` | 11 | Norm, Zerlegung, e/f/g-Satz |
| `TestKlassengruppe` | 13 | h(d), Minkowski, PID, Gruppenstruktur |
| `TestEinheitengruppe` | 9 | Rang, Fundamentaleinheit, Regulator |
| `TestQuadraticField` | 14 | Klasse, Diskriminante, Einheitswurzeln |
| `TestSymbole` | 12 | Legendre, Kronecker, Jacobi, Reziprozität |
| `TestCyclotomicField` | 12 | Grad, Signatur, Zerlegung, Kummer |
| `TestLokaleTheorie` | 9 | p-adisch, Hensel, Hilbert-Symbol |
| `TestKlassenzahlformel` | 7 | Formelwerte, L-Werte, Regulator |
| **Gesamt** | **110+** | **alle grün** |
