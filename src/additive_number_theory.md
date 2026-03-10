# Modul: Additive Zahlentheorie (`additive_number_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul implementiert die zentralen Gebiete der additiven Zahlentheorie: die Theorie der
ganzzahligen Partitionen, das Waring-Problem, die Goldbach-Vermutung, arithmetische
Progressionen in den Primzahlen sowie Summenmengen-Theorie. Alle Algorithmen sind didaktisch
kommentiert und dienen als Lernmaterial für additive Kombinatorik und analytische Zahlentheorie.

## Mathematischer Hintergrund

### Partitionen

Die Partitionsfunktion $p(n)$ zählt die Anzahl der Möglichkeiten, $n$ als geordnete Summe
positiver ganzer Zahlen zu schreiben (Reihenfolge irrelevant):

$$p(4) = 5: \quad 4 = 3+1 = 2+2 = 2+1+1 = 1+1+1+1$$

**Euler'sche erzeugende Funktion:**

$$\sum_{n=0}^{\infty} p(n)\, x^n = \prod_{k=1}^{\infty} \frac{1}{1-x^k}$$

**Euler-Pentagonalzahltheorem (Rekurrenz):**

$$p(n) = \sum_{k \neq 0} (-1)^{k+1}\, p\!\left(n - \frac{k(3k-1)}{2}\right)$$

mit den verallgemeinerten Pentagonalzahlen $\omega(k) = \frac{k(3k-1)}{2}$.

**Ramanujan-Kongruenzen (1919):**

$$p(5k+4) \equiv 0 \pmod{5}, \quad p(7k+5) \equiv 0 \pmod{7}, \quad p(11k+6) \equiv 0 \pmod{11}$$

### Hardy-Ramanujan-Asymptotik

$$p(n) \sim \frac{1}{4n\sqrt{3}} \exp\!\left(\pi\sqrt{\frac{2n}{3}}\right)$$

**Rademacher-Formel (exakt konvergent, 1937):**

$$p(n) = \frac{1}{\pi\sqrt{2}} \sum_{k=1}^{\infty} A_k(n)\,\sqrt{k}\,
  \frac{d}{dn}\!\left[\frac{\sinh\!\left(\frac{\pi}{k}\sqrt{\frac{2}{3}\!\left(n-\frac{1}{24}\right)}\right)}
  {\sqrt{n-\frac{1}{24}}}\right]$$

### Waring-Problem

Für jeden Exponenten $k \geq 1$ existiert $g(k)$: die minimale Anzahl $k$-ter Potenzen,
die jede natürliche Zahl darstellen. Bekannte exakte Werte:

| $k$ | $g(k)$ | Satz |
|-----|--------|------|
| 2 | 4 | Lagrange (1770) |
| 3 | 9 | Wieferich & Kempner (1909) |
| 4 | 19 | Balasubramanian et al. (1986) |

Asymptotisches Analogon: $G(k) = $ minimale Anzahl für alle *hinreichend großen* $n$.
$G(4) = 15$ (Davenport, 1939).

**Legendre-Dreisatz:** $n = a^2 + b^2 + c^2$ hat eine Lösung genau dann wenn
$n \neq 4^a(8b+7)$.

### Goldbach-Vermutung

$$\forall\, n > 2,\; n \text{ gerade}: \quad n = p + q \quad (p, q \text{ prim})$$

Bis $n \leq 4 \cdot 10^{18}$ rechnerisch verifiziert, analytisch unbewiesen.

**Schwache Goldbach-Vermutung (Helfgott, 2013):**

$$\forall\, n > 5,\; n \text{ ungerade}: \quad n = p_1 + p_2 + p_3$$

**Chen-Theorem (1966):** Jede hinreichend große gerade Zahl $n = p + m$, wobei $m$
höchstens zwei Primfaktoren besitzt.

### Green-Tao-Satz

$$\text{Die Primzahlen enthalten arithmetische Progressionen beliebiger Länge.}$$

---

## Klassen und Methoden

### `PartitionFunction`

**Beschreibung:** Berechnet die Partitionsfunktion $p(n)$ und erzeugt Partitionen in
verschiedenen Varianten. Implementiert dynamische Programmierung, den Münz-Wechsel-Ansatz
und das Euler-Pentagonalzahltheorem.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `partition_count()` | `(n: int) -> int` | $p(n)$ via DP (Münz-Wechsel) |
| `generate_partitions()` | `(n: int) -> list[tuple]` | Alle Partitionen (absteigend sortiert) |
| `partition_with_parts()` | `(n: int, k: int) -> list[tuple]` | Partitionen in genau $k$ Teile |
| `partition_into_distinct_parts()` | `(n: int) -> list[tuple]` | Partitionen in paarweise verschiedene Teile |
| `euler_pentagonal_theorem()` | `(n: int) -> int` | $p(n)$ via Pentagonalzahl-Rekurrenz |
| `ramanujan_partition_congruences()` | `(p_vals: list[int]) -> dict` | Verifiziert alle drei Ramanujan-Kongruenzen |

---

### `HardyRamanujanCircleMethod`

**Beschreibung:** Implementiert asymptotische Formeln und die Kreismethode für $p(n)$.
Der Integrationskreis wird in Farey-Bögen um rationale Punkte $h/k$ aufgeteilt.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `asymptotic_partition()` | `(n: int) -> float` | Hardy-Ramanujan-Asymptotik $\frac{1}{4n\sqrt{3}} e^{\pi\sqrt{2n/3}}$ |
| `rademacher_formula()` | `(n: int, k_max: int = 10) -> float` | Rademacher-Reihe mit Dedekind-Summen und Kloosterman-Termen |
| `farey_arc_contributions()` | `(n: int, k_max: int = 5) -> list[dict]` | Amplituden und Phasen der Farey-Bogen-Beiträge |

---

### `WaringProblem`

**Beschreibung:** Berechnet $g(k)$- und $G(k)$-Werte und stellt Zahlen als Summen
von $k$-ten Potenzen dar. Enthält Lagrange-Vier-Quadrate-Darstellung und Kuben-Greedy.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `g_k()` | `(k: int) -> Optional[int]` | Gibt $g(k)$ zurück (bekannt für $k \leq 10$) |
| `G_k()` | `(k: int) -> Optional[int]` | Gibt $G(k)$ zurück (asymptotisch optimal) |
| `represent_as_squares()` | `(n: int, k: int = 4) -> Optional[tuple]` | $n = a_1^2 + \ldots + a_k^2$ via Backtracking |
| `represent_as_cubes()` | `(n: int) -> Optional[tuple]` | $n$ als Summe von Kuben (Greedy + Backtracking) |
| `four_squares_theorem()` | `(n: int) -> tuple[int,int,int,int]` | Lagrange-Darstellung $n = a^2+b^2+c^2+d^2$ |
| `three_squares_theorem_check()` | `(n: int) -> bool` | Legendre-Kriterium: $n \neq 4^a(8b+7)$ |

---

### `GoldbachAnalysis`

**Beschreibung:** Analysiert die Goldbach-Vermutung: findet alle Zerlegungen gerader
Zahlen in zwei Primzahlen, erzeugt Komet-Daten, prüft die ternäre Vermutung (Helfgott)
und das Chen-Theorem.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `goldbach_decompositions()` | `(n: int) -> list[tuple[int,int]]` | Alle Paare $(p, q)$ mit $p+q=n$, beide prim |
| `goldbach_count()` | `(n: int) -> int` | Anzahl $r(n)$ der Goldbach-Zerlegungen |
| `goldbach_comet_data()` | `(n_max: int) -> list[tuple[int,int]]` | Goldbach-Komet: $(n, r(n))$ für gerade $n \leq n_\text{max}$ |
| `ternary_goldbach_check()` | `(n: int) -> Optional[tuple]` | $n = p_1 + p_2 + p_3$ (schwache Vermutung, Helfgott 2013) |
| `chen_theorem_check()` | `(n: int) -> Optional[tuple]` | $n = p + m$ mit $\Omega(m) \leq 2$ (Chen 1966) |

---

### `ArithmeticProgressions`

**Beschreibung:** Untersucht arithmetische Progressionen in Primzahlen und ganzen
Zahlen. Implementiert Dirichlet-Progressionen, Van-der-Waerden-Zahlen und
Green-Tao-Suche.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `dirichlet_primes_in_ap()` | `(a: int, d: int, limit: int) -> list[int]` | Primzahlen der Form $a + kd$ (Dirichlet-Satz) |
| `prime_ap_of_length()` | `(k: int, limit: int) -> Optional[tuple]` | Primzahl-AP der Länge $k$ (Green-Tao-Suche) |
| `van_der_waerden_number()` | `(k: int, r: int) -> Optional[int]` | Bekannte Van-der-Waerden-Zahlen $W(k;r)$ |
| `check_van_der_waerden()` | `(coloring: list, k: int) -> bool` | Prüft ob eine $r$-Färbung eine einfarbige $k$-AP enthält |
| `consecutive_prime_aps()` | `(length: int, limit: int) -> list` | Konsekutive Primzahl-APs bis `limit` |

---

### `SumSetTheory`

**Beschreibung:** Implementiert Summenmengen-Operationen und additive Kombinatorik:
Cauchy-Davenport-Schranken, additive Energie, Freiman-Ruzsa-Koeffizient und
Sidon-Mengen.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `sumset()` | `(A: set, B: set) -> set` | $A + B = \{a+b : a \in A, b \in B\}$ |
| `difference_set()` | `(A: set, B: set) -> set` | $A - B = \{a-b : a \in A, b \in B\}$ |
| `iterated_sumset()` | `(A: set, k: int) -> set` | $kA = A + A + \ldots + A$ ($k$-fach) |
| `cauchy_davenport_bound()` | `(A: set, B: set, p: int) -> int` | $\min(p, |A|+|B|-1)$: untere Schranke für $|A+B|$ in $\mathbb{Z}/p\mathbb{Z}$ |
| `additive_energy()` | `(A: set, B: set) -> int` | $E(A,B) = |\{(a_1,a_2,b_1,b_2): a_1+b_1=a_2+b_2\}|$ |
| `freiman_ruzsa_coefficient()` | `(A: set) -> float` | Ruzsa-Verdopplungskoeffizient $|A+A| / |A|$ |
| `is_sidon_set()` | `(A: set) -> bool` | Prüft: alle Summen $a+b$ paarweise verschieden (Sidon/B₂-Menge) |
| `find_sidon_set()` | `(n: int, size: int) -> set` | Sucht Sidon-Menge der Größe `size` in $\{1,\ldots,n\}$ |

---

### Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `schur_number()` | `(r: int) -> Optional[int]` | Bekannte Schur-Zahlen $S(r)$ |
| `schur_coloring_check()` | `(coloring: list) -> bool` | Prüft ob Färbung summenfreie Klassen hat |
| `happy_ending_problem()` | `(points: list) -> Optional[tuple]` | Konvexe 4-Punkte-Teilmenge (Erdős-Szekeres) |
| `bertrand_postulate_verify()` | `(n: int) -> tuple[int, ...]` | Primzahlen in $(n, 2n)$ (Bertrand'sches Postulat) |

---

## Beispiele

```python
from additive_number_theory import PartitionFunction, GoldbachAnalysis, WaringProblem

# Partitionsfunktion
pf = PartitionFunction()
print(pf.partition_count(10))       # 42
print(pf.generate_partitions(4))    # [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]

# Ramanujan-Kongruenzen verifizieren
results = pf.ramanujan_partition_congruences([0, 1, 2, 3])
# p(4) = 5 ≡ 0 (mod 5), p(5+4)=p(9)=30 ≡ 0 (mod 5), ...

# Goldbach-Zerlegungen
ga = GoldbachAnalysis()
decomps = ga.goldbach_decompositions(28)
# [(5, 23), (11, 17)]

# Vier-Quadrate-Darstellung
wp = WaringProblem()
print(wp.four_squares_theorem(23))  # (4, 2, 2, 1) → 16+4+2+1 = 23 ✗
# → (3, 3, 2, 1) → 9+9+4+1 = 23 ✓

# Ternäre Goldbach-Vermutung
print(ga.ternary_goldbach_check(77))  # z.B. (3, 7, 67)
```

## Verweise

- Verwandte Module: `analytic_number_theory.py`, `proof_theory.py`, `algebra.py`
- Literatur:
  - Hardy & Ramanujan: *Asymptotic Formulae in Combinatory Analysis* (1918)
  - Rademacher: *On the Partition Function* (1937)
  - Helfgott: *Major arcs for Goldbach's theorem* (2013)
  - Green & Tao: *The primes contain arbitrarily long arithmetic progressions* (2004)
  - Nathanson: *Additive Number Theory* (Springer, 1996)
