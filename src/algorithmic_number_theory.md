# Modul: Algorithmische Zahlentheorie (`algorithmic_number_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul implementiert klassische und moderne Algorithmen der algorithmischen Zahlentheorie:
Primzahlsiebe (Eratosthenes, Atkin), probabilistische und deterministische Primzahltests,
Faktorisierungsverfahren (Probedivision, Pollard-$\rho$, Fermat), diskrete Logarithmen
(BSGS, Pohlig-Hellman, Indexkalkül) sowie Kettenbrüche und diophantische Approximation.
Wichtige Basisfunktionen (`miller_rabin_primality_test`, `jacobi_symbol`) werden aus
`proof_theory.py` importiert, um Duplikate zu vermeiden.

## Mathematischer Hintergrund

### Primzahlsiebe

**Sieb des Eratosthenes** – Zeitkomplexität $O(n \log\log n)$:

$$\text{Markiere alle Vielfachen } p^2, p^2+p, p^2+2p, \ldots \text{ für jedes Prim } p \leq \sqrt{n}$$

**Sieb von Atkin** – Zeitkomplexität $O(n / \log\log n)$:
Verwendet drei quadratische Formen:

$$4x^2 + y^2 \equiv 1 \pmod{4}, \quad 3x^2 + y^2 \equiv 7 \pmod{12}, \quad 3x^2 - y^2 \equiv 11 \pmod{12}$$

### Primzahltests

**Fermat-Test:** Kleiner Satz von Fermat — für Primzahl $p$ gilt $a^{p-1} \equiv 1 \pmod{p}$.
Carmichael-Zahlen (z.B. 561) täuschen diesen Test für *alle* Basen.

**Solovay-Strassen:** Euler-Kriterium via Jacobi-Symbol:

$$a^{(n-1)/2} \equiv \left(\frac{a}{n}\right) \pmod{n}$$

Fehlerwahrscheinlichkeit $\leq (1/2)^k$ nach $k$ Runden.

**Lucas-Test:** $n$ ist prim, wenn eine primitive Wurzel $a$ modulo $n$ existiert:

$$a^{n-1} \equiv 1 \pmod{n} \quad \text{und} \quad a^{(n-1)/q} \not\equiv 1 \pmod{n} \text{ für alle Primteiler } q \mid n-1$$

**Korselt-Kriterium** (Carmichael-Zahlen): $n$ zusammengesetzt, quadratfrei,
$(p-1) \mid (n-1)$ für jeden Primteiler $p$ von $n$.

### Faktorisierung

**Pollard-$\rho$** (Floyd-Zykluserkennung): Erwartet $O(n^{1/4})$ Schritte:

$$f(x) = x^2 + c \pmod{n}, \quad d = \gcd(|x_\text{hase} - x_\text{schildkröte}|,\, n)$$

**Pollard's $p-1$:** Effektiv wenn $p-1$ $B$-glatt ist:

$$a \leftarrow 2^{\prod_{p \leq B} p^{\lfloor \log_B p \rfloor}} \pmod{n}, \quad d = \gcd(a-1, n)$$

**Fermat-Faktorisierung:** Sucht $n = a^2 - b^2 = (a+b)(a-b)$, effizient wenn $p \approx q$.

### Diskreter Logarithmus

Gesucht: $x$ mit $g^x \equiv h \pmod{p}$.

**Baby-Step-Giant-Step** – $O(\sqrt{p})$ Zeit und Speicher:

$$x = im - j, \quad m = \lceil\sqrt{p-1}\rceil$$
Baby Steps: $h \cdot g^j$; Giant Steps: $g^{im}$; Kollisionssuche.

**Pohlig-Hellman:** Nutzt CRT wenn $p-1 = q_1^{e_1} \cdots q_r^{e_r}$ (glatt):

$$x \pmod{q_i^{e_i}} \text{ via BSGS in Untergruppe der Ordnung } q_i^{e_i}$$

**Indexkalkül** – subexponentiell $L[1/2]$: Faktorbasis, Beziehungsphase, lineare Algebra.

### Kettenbrüche

$$x = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cdots}} = [a_0; a_1, a_2, \ldots]$$

Konvergente: $p_n = a_n p_{n-1} + p_{n-2}$, $q_n = a_n q_{n-1} + q_{n-2}$.
Anwendung: Faktorisierung via Kettenbruch-Methode (CFRAC).

---

## Klassen und Methoden

### `PrimeSieve`

**Beschreibung:** Enthält verschiedene Siebverfahren zur Primzahlgenerierung sowie
Analysefunktionen für Primzahllücken, Zwillingsprimzahlen und die Primzahlzählfunktion.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `sieve_of_eratosthenes()` | `(n: int) -> list[int]` | Alle Primzahlen $\leq n$ via Eratosthenes (delegiert an `proof_theory`) |
| `sieve_of_atkin()` | `(n: int) -> list[int]` | Moderneres Atkin-Sieb, effizienter für große $n$ |
| `prime_gaps()` | `(n: int) -> list[tuple]` | Liste $(p_1, p_2, p_2 - p_1)$ aufeinanderfolgender Primzahlpaare |
| `twin_primes()` | `(n: int) -> list[tuple[int,int]]` | Alle Zwillingsprimzahlpaare $(p, p+2)$ mit $p \leq n$ |
| `prime_counting_exact()` | `(n: int) -> int` | $\pi(n)$ — exakte Anzahl der Primzahlen $\leq n$ |

---

### `PrimalityTests`

**Beschreibung:** Sammlung probabilistischer und deterministischer Primzahltests.
Miller-Rabin aus `proof_theory` wird als deterministischer Fallback verwendet.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `aks_test_demo()` | `(n: int) -> dict` | AKS-Demo: perfekte Potenz-Test + deterministischer Miller-Rabin |
| `solovay_strassen()` | `(n: int, k: int = 10) -> bool` | Euler-Kriterium via Jacobi-Symbol, Fehler $\leq (1/2)^k$ |
| `fermat_test()` | `(n: int, a: int) -> bool` | Fermat-Test zur Basis $a$ (Carmichael-Zahlen täuschen!) |
| `lucas_primality_test()` | `(n: int) -> dict` | Sucht primitive Wurzel als Lucas-Zeugen |
| `is_carmichael_number()` | `(n: int) -> bool` | Korselt-Kriterium: quadratfrei und $(p-1) \mid (n-1)$ |

---

### `IntegerFactorization`

**Beschreibung:** Kombiniert mehrere Faktorisierungsverfahren. `factor_completely()` wählt
automatisch die geeignete Methode (Probedivision → Pollard-$\rho$ → Fermat).

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `trial_division()` | `(n: int) -> dict` | Probedivision bis $\sqrt{n}$, $O(\sqrt{n})$ |
| `pollard_rho()` | `(n: int, seed: int = 2) -> Optional[int]` | Floyd-Zykluserkennung, erwartet $O(n^{1/4})$ |
| `pollard_p_minus_1()` | `(n: int, B: int = 100) -> Optional[int]` | Effektiv für $B$-glatte $p-1$ |
| `fermat_factoring()` | `(n: int) -> Optional[tuple]` | $n = a^2 - b^2$, effizient bei $p \approx q$ |
| `factor_completely()` | `(n: int) -> dict` | Vollständige Primfaktorisierung (kombinierte Strategie) |
| `smooth_number_check()` | `(n: int, B: int) -> bool` | Prüft ob $n$ $B$-glatt ist (alle Primteiler $\leq B$) |

---

### `DiscreteLogarithm`

**Beschreibung:** Verschiedene Algorithmen zur Lösung von $g^x \equiv h \pmod{p}$.
Von Brute Force über Baby-Step-Giant-Step bis zum Pohlig-Hellman-Algorithmus.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `baby_step_giant_step()` | `(g: int, h: int, p: int) -> Optional[int]` | BSGS: $O(\sqrt{p})$ Zeit und Speicher |
| `pohlig_hellman()` | `(g: int, h: int, p: int) -> Optional[int]` | CRT-basiert, effizient wenn $p-1$ glatt |
| `index_calculus_demo()` | `(g: int, h: int, p: int, B: int = 20) -> dict` | Indexkalkül-Demo: Faktorbasis + Beziehungsphase |
| `discrete_log_brute_force()` | `(g: int, h: int, p: int) -> Optional[int]` | $O(p)$, nur für $p < 10^6$ |

---

### `ContinuedFractions`

**Beschreibung:** Kettenbruchentwicklung, Konvergenten-Berechnung und Anwendungen
in der diophantischen Approximation und Faktorisierung.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `to_continued_fraction()` | `(x: float, n_terms: int) -> list[int]` | Entwicklung $[a_0; a_1, a_2, \ldots]$ |
| `convergents()` | `(cf: list[int]) -> list[tuple]` | Konvergente $(p_k, q_k)$ via Rekurrenz |
| `from_continued_fraction()` | `(cf: list[int]) -> float` | Rekonstruktion der Zahl aus Kettenbruch |
| `best_rational_approximation()` | `(x: float, max_denom: int) -> tuple` | Bester Näherungsbruch $p/q$ mit $q \leq$ `max_denom` |
| `cfrac_factoring_demo()` | `(n: int) -> Optional[int]` | Kettenbruch-basierte Faktorisierung (CFRAC-Demo) |

---

### `DiophantineApproximation`

**Beschreibung:** Dirichlet-Approximation, Hurwitz-Satz und simultane diophantische
Näherungen für irrationale Zahlen.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `dirichlet_approximation()` | `(x: float, Q: int) -> tuple` | Bester $p/q$ mit $q \leq Q$, $|x - p/q| < 1/qQ$ |
| `hurwitz_theorem_check()` | `(x: float, max_q: int) -> list` | Brüche mit $|x - p/q| < 1/(\sqrt{5}\, q^2)$ (Hurwitz) |
| `simultaneous_approximation()` | `(alphas: list, Q: int) -> tuple` | Simultane Approximation mehrerer irrationaler Zahlen |
| `three_distance_theorem()` | `(alpha: float, n: int) -> list` | Drei-Abstands-Satz: Abstände $\{k\alpha\}$ auf dem Kreis |

---

## Beispiele

```python
from algorithmic_number_theory import (
    PrimeSieve, PrimalityTests, IntegerFactorization, DiscreteLogarithm
)

# Atkin-Sieb
ps = PrimeSieve()
primes = ps.sieve_of_atkin(100)
# [2, 3, 5, 7, 11, 13, ..., 97]

twins = ps.twin_primes(50)
# [(3, 5), (5, 7), (11, 13), (17, 19), (29, 31), (41, 43)]

# Solovay-Strassen
pt = PrimalityTests()
print(pt.solovay_strassen(997, k=15))   # True (prim)
print(pt.is_carmichael_number(561))     # True (561 = 3×11×17)

# Vollständige Faktorisierung
fact = IntegerFactorization()
result = fact.factor_completely(360)
# {'factors': {2: 3, 3: 2, 5: 1}, 'factorization_str': '2^3 × 3^2 × 5'}

# Diskreter Logarithmus
dl = DiscreteLogarithm()
x = dl.baby_step_giant_step(g=2, h=22, p=29)
# x=19: 2^19 = 524288 ≡ 22 (mod 29) ✓

# Pohlig-Hellman
x2 = dl.pohlig_hellman(g=3, h=13, p=17)
```

## Verweise

- Verwandte Module: `proof_theory.py`, `algebra.py`, `p_adic.py`
- Literatur:
  - Crandall & Pomerance: *Prime Numbers: A Computational Perspective* (Springer, 2005)
  - Cohen: *A Course in Computational Algebraic Number Theory* (Springer, 1993)
  - Lenstra: *Factoring Integers with Elliptic Curves* (1987)
  - Agrawal, Kayal, Saxena: *PRIMES is in P* (2004)
