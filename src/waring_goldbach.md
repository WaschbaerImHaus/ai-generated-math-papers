# waring_goldbach.py – Waring-Goldbach-Problem

**Autor:** Michael Fuhrmann
**Erstellt:** 2026-03-12
**Letzte Änderung:** 2026-03-12
**Build:** 122

---

## Überblick

Das **Waring-Goldbach-Problem** (benannt nach Edward Waring und Christian
Goldbach) fragt: Wie viele Primzahlpotenzen $k$-ter Ordnung genügen, um jede
hinreichend große natürliche Zahl darzustellen?

---

## Mathematischer Hintergrund

### Warings Problem (klassisch, ohne Primzahlbedingung)

Waring (1770) behauptete:
- Jede positive ganze Zahl = Summe von ≤ 4 Quadraten (Lagrange 1770, bewiesen)
- Jede positive ganze Zahl = Summe von ≤ 9 Kubiken
- Jede positive ganze Zahl = Summe von ≤ 19 Viertenpotenzen

**Notation:**
- $g(k)$ = minimale Anzahl $k$-ter Potenzen für **alle** positiven ganzen Zahlen
- $G(k)$ = minimale Anzahl $k$-ter Potenzen für alle **hinreichend großen** Zahlen

### Waring-Goldbach (Primzahlpotenzen)

**Einschränkung:** Die Summanden müssen Potenzen von **Primzahlen** sein.

**Notation:**
- $G_{\text{prim}}(k)$ = minimale Anzahl $k$-ter Primzahlpotenzen für hinr. große $n$

---

## Bekannte Resultate

### k = 1: Goldbach und Vinogradov

**Goldbach-Vermutung (1742, CONJECTURE, offen):**
$$\text{Jede gerade Zahl} \geq 4 = p_1 + p_2 \quad (p_1, p_2 \text{ prim})$$

**Vinogradov (1937), vervollständigt Helfgott (2013):**
$$\text{Jede ungerade Zahl} \geq 7 = p_1 + p_2 + p_3 \quad (p_i \text{ prim})$$
(Numerisch verifiziert bis $4 \cdot 10^{18}$ für die Goldbach-Vermutung)

### k = 2: Quadrate (Hua 1938)

$$\text{Jede hinreichend große } n = p_1^2 + p_2^2 + p_3^2 + p_4^2 + p_5^2$$

**Hua's Theorem:** $G(2) \leq 5$

**Lokale Bedingungen:** Nicht jede Zahl ist als Summe von Primzahlquadraten
darstellbar. Die Darstellbarkeit hängt von $n \pmod{8}$ und $n \pmod{3}$ ab.

**Wichtiger Hinweis:** $n = 10$ ist **nicht** als Summe von Primzahlquadraten
darstellbar (die einzigen Primzahlquadrate $\leq 10$ sind $4 = 2^2$ und $9 = 3^2$,
aber $4+4=8, 4+9=13 \neq 10$).

### k = 3: Kuben (Hua 1938, Vaughan-Wooley)

$$G(3) \leq 9 \quad \text{(Hua)}, \quad G(3) \leq 7 \quad \text{(verbessert)}$$

### Allgemeine Schranken

| $k$ | $G(k)$ Hua (1938) | Verbessert |
|-----|------------------|------------|
| 1   | 3 (Vinogradov)   | 3          |
| 2   | 5                | 5          |
| 3   | 9                | 7          |
| 4   | 15               | 13         |
| 5   | 21               | 15         |
| 6   | 33               | 19         |

**Vinogradov-Wooley-Formel (Efficient Congruencing, Wooley 2019):**
$$G(k) \leq k(3\log k + 4.2 \log\log k + C)$$

---

## Methodik: Kreismethode (Hardy-Littlewood-Vinogradov)

Die **Kreismethode** (Hardy-Ramanujan 1918, verfeinert von Vinogradov) ist das
Hauptwerkzeug für Waring-Goldbach-Probleme.

### Grundidee

Für $r_s(n)$ = Anzahl der Darstellungen $n = p_1^k + \cdots + p_s^k$:

$$r_s(n) = \int_0^1 f(\alpha)^s e^{-2\pi i \alpha n}\, d\alpha$$

wobei $f(\alpha) = \sum_{p \leq N} e^{2\pi i \alpha p^k}$ (Exponentialsumme über Primzahlen).

### Singuläre Reihe

Die **Asymptotik** für $r_s(n)$ bei $s \geq G(k)$:
$$r_s(n) \sim \mathfrak{S}(n) \cdot \frac{\Gamma(1 + 1/k)^s}{\Gamma(s/k)} \cdot n^{s/k - 1} / (\log n)^s$$

wobei $\mathfrak{S}(n)$ die **singuläre Reihe** ist:
$$\mathfrak{S}(n) = \prod_p \left(1 + \sum_{m=1}^\infty \frac{\chi_k(p^m, n)}{p^{ms}}\right) > 0$$

---

## Klassen

### `WaringGoldbachK2`

Spezialisiert auf $k=2$ (Primzahlquadrate).

| Methode | Beschreibung |
|---------|-------------|
| `local_conditions(n)` | Lokale Bedingungen mod 8 und mod 24 |
| `represent_as_prime_squares(n, max_terms)` | Findet Darstellung $n = \sum p_i^2$ |
| `verify_hua_theorem(limit)` | Numerische Verifikation bis `limit` |
| `count_representations(n, num_terms)` | Zählt Darstellungen |
| `minimal_representation(n)` | Minimale Termsanzahl |
| `singular_series_heuristic(n)` | Heuristische Schätzung von $\mathfrak{S}(n)$ |

### `WaringGoldbachGeneral`

Allgemeines Waring-Goldbach für beliebige $k \geq 1$.

| Methode | Beschreibung |
|---------|-------------|
| `g_bound_hua()` | Hua (1938) Schranke für $G(k)$ |
| `g_bound_vinogradov_wooley()` | Vinogradov-Wooley Schranke |
| `represent_as_k_prime_powers(n, max_terms)` | Darstellung als $k$-te Primzahlpotenzen |
| `vinogradov_three_primes_check(n)` | Vinogradov/Helfgott für $k=1$ |
| `goldbach_conjecture_check(n)` | Goldbach-Vermutung für gerade $n$ |
| `waring_goldbach_statistics(limit)` | Statistiken über Darstellungen |
| `hua_theorem_statement()` | Hua-Theorem als String |
| `explicit_vinogradov_wooley_bound()` | Formel als String |

---

## Verwendungsbeispiele

```python
from waring_goldbach import WaringGoldbachK2, WaringGoldbachGeneral

# k=2: Primzahlquadrate
wg2 = WaringGoldbachK2(max_prime=200)

# 13 = 2² + 3² = 4 + 9
rep = wg2.represent_as_prime_squares(13, max_terms=2)
print(rep)   # [2, 3]
print([p**2 for p in rep])  # [4, 9]

# Wie viele Terme braucht 38?
k, primes = wg2.minimal_representation(38)
print(f"38 = Summe von {k} Primzahlquadraten: {primes}")

# k=3: Kubiken von Primzahlen
wg3 = WaringGoldbachGeneral(k=3, max_prime=100)
print(wg3.g_bound_hua())         # 9
print(wg3.hua_theorem_statement())

# Goldbach-Vermutung prüfen
wg1 = WaringGoldbachGeneral(k=1, max_prime=500)
result = wg1.goldbach_conjecture_check(100)
print(result['representation'])  # [3, 97] oder ähnlich

# Vinogradov für ungerade Zahlen
result = wg1.vinogradov_three_primes_check(15)
print(result['representation'])  # [3, 5, 7]
```

---

## Lokale Obstruktionen

### Warum ist k=2 subtiler?

Nicht alle natürlichen Zahlen sind als Summen von Primzahlquadraten darstellbar:

- $n = 10$: Nicht darstellbar (Primzahlquadrate $\leq 10$: nur $4, 9$; aber $4+4=8, 4+9=13$)
- $n = 1, 2, 3, 5, 6, 7, 10, 11, 14, 15, ...$: Nicht als Summe von Primzahlquadraten darstellbar

Darstellbare Zahlen: $4, 8, 9, 12, 13, 16, 17, 18, 21, 22, 25, 26, 27, 29, ...$

Die Menge der **nicht** darstellbaren Zahlen ist jedoch endlich für jeden festen $s \geq G(k)$.

---

## Offene Probleme

1. **Goldbach-Vermutung** (CONJECTURE, 1742, offen):
   Jede gerade Zahl $\geq 4 = p_1 + p_2$

2. **Exakter Wert von G(3)**:
   Hua zeigt $G(3) \leq 9$, aber der genaue Wert ist unbekannt.

3. **Linnik's Constant** für Waring-Goldbach:
   Ab welchem $n$ gilt die asymptotische Formel für $r_s(n)$?

---

## Literatur

- Hua, L.K. (1938). *On Waring's problem.* Quarterly J. Math. 9.
- Vinogradov, I.M. (1937). *Some theorems concerning the primes.* Mat. Sbornik.
- Wooley, T.D. (2019). *Nested efficient congruencing and relatives of Vinogradov's mean value theorem.* Proc. LMS.
- Vaughan, R.C. (1997). *The Hardy-Littlewood method.* Cambridge University Press.
- Helfgott, H.A. (2013). *Major arcs for Goldbach's theorem.* arXiv:1305.2897.
