# pillai_chowla.py – Dokumentation

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12

---

## Überblick

Dieses Modul untersucht die **Pillai-Chowla-Vermutung** über ggT-Eigenschaften von
Fakultäten sowie zugehörige klassische zahlentheoretische Sätze: Wilson-Theorem,
Mertens' 3. Theorem und das Sylvester-Produkt.

---

## Pillai-Chowla-Vermutung

### Schwache Form (CONJECTURE, unbewiesen)

$$\gcd(n!+1,\; (n+1)!+1) = 1 \quad \text{für alle } n \geq 1$$

Numerisch verifiziert für alle $n \leq 200$ (diese Implementierung).

### Stärkere Form (FALSCH für allgemeine $n \neq m$)

Die naiv vermutete stärkere Form $\gcd(n!+1, m!+1) = 1$ für alle $n \neq m$
gilt **nicht**! Bekannte Gegenbeispiele:

| $n$ | $m$ | $\gcd(n!+1,\; m!+1)$ |
|-----|-----|----------------------|
| 3 | 6 | 7 |
| 5 | 10 | 11 |
| 7 | 9 | 71 |

**Erklärung:** Sei $p$ prim mit $p \mid n!+1$. Nach Wilson-Lemma gilt $p > n$.
Falls zusätzlich $p \leq m$, dann gilt $p \mid m!$ (Primzahl erscheint im Faktorprodukt),
also $p \mid m! + 1$. Daher: $p = \gcd(n!+1, m!+1) > 1$.

---

## Wilson-Theorem und Verwandtes

### Wilson-Theorem (BEWIESENES THEOREM)

> **THEOREM** (Wilson 1770, Beweis Lagrange 1771):
>
> Eine natürliche Zahl $p > 1$ ist **genau dann prim**, wenn gilt:
> $$(p-1)! \equiv -1 \pmod{p}$$

**Verbindung zu Pillai-Chowla:**
Falls $p \mid n!+1$ mit $p$ prim und $n < p$, dann $(n! \equiv -1 \pmod{p})$.
Wilson-Theorem besagt dies nur für $n = p-1$. Für $n < p-1$ sind solche Primzahlen
"zufällig" und bilden eine seltene Klasse.

### Wilson-Primzahlen

> **DEFINITION:** $p$ heißt **Wilson-Primzahl**, wenn $p^2 \mid (p-1)!+1$.

Bekannte Wilson-Primzahlen (bis $2 \times 10^{13}$): **nur 5, 13, 563**.

**CONJECTURE:** Es gibt unendlich viele Wilson-Primzahlen (unbewiesen).

### Wilson-Lemma (BEWIESENES LEMMA)

> **LEMMA:** Sei $n \geq 1$. Jeder Primfaktor $p$ von $n!+1$ erfüllt $p > n$.
>
> *Beweis:* Sei $p \leq n$ prim. Dann $p \mid n!$, also $p \mid n!+1 - n! = 1$. Widerspruch.

---

## Sylvester-Produkt und verwandte Folgen

### Primorial-Definition

$$p_n\# = 2 \cdot 3 \cdot 5 \cdots p_n \quad \text{(Produkt der ersten } n \text{ Primzahlen)}$$

Verbindung zu Euklid-Beweis: $p_n\# + 1$ hat einen Primfaktor $> p_n$.

### Echte Sylvester-Folge

Rekursion: $a_1 = 2$, $a_{n+1} = a_n^2 - a_n + 1$

$$2, 3, 7, 43, 1807, 3263443, 10650056950807, ...$$

Eigenschaft: $\sum_{k=1}^n 1/a_k < 1$ (ägyptische Brüche, optimale Darstellung).

---

## Mertens' 3. Theorem

### Aussage (BEWIESENES THEOREM)

> **THEOREM** (Mertens 1874):
>
> $$\prod_{p \leq n} \frac{p}{p-1} \sim e^{\gamma} \ln n \quad \text{für } n \to \infty$$
>
> wobei $\gamma$ die Euler-Mascheroni-Konstante ist.

Äquivalente Form:
$$\prod_{p \leq n} \left(1 - \frac{1}{p}\right)^{-1} \sim e^{\gamma} \ln n$$

**Verbindung zu $\gamma$:** Das Mertens-Produkt ermöglicht, $\gamma$ numerisch aus der
Primzahlverteilung zu schätzen:
$$\gamma \approx \ln\left(\prod_{p \leq n} \frac{p}{p-1}\right) - \ln(\ln n)$$

Die Konvergenz ist langsam (logarithmisch), aber numerisch illustrativ.

---

## Klasse: PillaiChowla

### Initialisierung

```python
pc = PillaiChowla(max_n=200)
```

### Methoden

| Methode | Beschreibung | Rückgabe |
|---------|-------------|---------|
| `factorial_plus_one(n)` | $n!+1$ (gecacht) | `int` |
| `gcd_consecutive(n)` | $\gcd(n!+1, (n+1)!+1)$ | `int` |
| `verify_consecutive(max_n)` | Pillai-Chowla schwach bis `max_n` | `dict` |
| `gcd_arbitrary(n, m)` | $\gcd(n!+1, m!+1)$ | `int` |
| `verify_general(max_val)` | Alle Paare $n < m \leq \text{max\_val}$ | `dict` |
| `wilson_theorem_verify(up_to)` | Wilson für alle $p \leq \text{up\_to}$ | `dict` |
| `wilson_primes(up_to)` | Findet Wilson-Primzahlen $\leq \text{up\_to}$ | `List[int]` |
| `sylvester_product(n)` | $1 + \prod_{k=1}^n p_k$ | `int` |
| `sylvester_sequence(n)` | Echte Sylvester-Folge | `List[int]` |
| `primorial_plus_one(n)` | $p_n\# + 1$ | `int` |
| `mertens_third_theorem_numerical(n_vals)` | Numerische Mertens-Verifikation | `dict` |
| `mertens_connection_to_gamma(max_n)` | $\gamma$-Schätzung aus Mertens | `dict` |
| `prime_factors_of_factorial_plus_one(n)` | Primfaktoren von $n!+1$ | `List[int]` |
| `wilson_lemma_verify(max_n)` | Verifiziert Wilson-Lemma | `dict` |
| `summary_open_problems()` | Offene Probleme | `List[dict]` |

---

## Offene Probleme (Stand 2026)

| Problem | Status |
|---------|--------|
| Pillai-Chowla (schwach): $\gcd(n!+1,(n+1)!+1)=1$ | **CONJECTURE** – unbewiesen |
| Unendlich viele Wilson-Primzahlen | **CONJECTURE** – unbewiesen |
| Unendlich viele primoriale Primzahlen ($p_n\# + 1$) | **CONJECTURE** – unbewiesen |
| $\gcd(n!+1, m!+1) = 1$ für alle $n \neq m$ | **FALSCH** (Gegenbeispiele bekannt) |

---

## Abhängigkeiten

- `math`: `gcd`, `factorial`, `log`
- `sympy`: `isprime`, `primerange`, `primorial`
- `mpmath`: Euler-Konstante für Mertens-Vergleich
- `numpy`: Array-Operationen
- `fractions`: Exakte Arithmetik
