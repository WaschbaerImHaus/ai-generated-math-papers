# erdos_selfridge.py – Erdős-Selfridge-Vermutung

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Erstellt**: 2026-03-12
**Status**: CONJECTURE – offen

---

## Übersicht

Das Modul `erdos_selfridge.py` untersucht die **Erdős-Selfridge-Vermutung** (1975):

> **CONJECTURE (Erdős-Selfridge, 1975)**
> Der Binomialkoeffizient $\binom{n}{k}$ ist niemals eine Primzahlpotenz $p^a$ (mit $p$ prim, $a \geq 1$) für $n \geq k + 2$ und $k \geq 2$.

---

## Mathematischer Hintergrund

### Binomialkoeffizient

$$\binom{n}{k} = \frac{n!}{k! \cdot (n-k)!} = \frac{n(n-1)\cdots(n-k+1)}{k!}$$

Die Vermutung sagt: Für $n \geq k+2$ und $k \geq 2$ ist dieser Ausdruck nie eine reine Primzahlpotenz.

### Triviale Fälle (außerhalb der Vermutung)

| Fall | Wert | Anmerkung |
|------|------|-----------|
| $k = 1$ | $\binom{n}{1} = n$ | Primzahlpotenz wenn $n = p^a$ – trivial, außerhalb der Vermutung |
| $k = n-1$ | $\binom{n}{n-1} = n$ | Symmetrisch zu $k=1$ |
| $n = k+1$ | $\binom{k+1}{k} = k+1$ | Zu klein: $n < k+2$ |

### Sylvesters Theorem (1892) – BEWIESEN

> Falls $n \geq 2k$, dann hat $\binom{n}{k}$ einen Primteiler $p > k$.

**Folgerung für die Vermutung**: Wenn $\binom{n}{k} = p^a$, dann wäre $p$ der einzige Primteiler. Mit Sylvester muss $p > k$ gelten. Aber das Produkt $n(n-1)\cdots(n-k+1)$ enthält viele Faktoren mit verschiedenen Primteilern, was $p^a$ für $a \geq 1$ typischerweise verhindert.

### Spezialfall $k = 2$

$$\binom{n}{2} = \frac{n(n-1)}{2} = p^a \iff n(n-1) = 2p^a$$

Da $\gcd(n, n-1) = 1$, müssten $n$ und $n-1$ die Form $\{2, p^a\}$ oder $\{1, 2p^a\}$ haben. Eine Lösung für $n \geq 4$ wurde nie gefunden.

### Granvilles Analyse

Andrew Granville untersuchte Spezialfälle:
- **$k$ prim**: Sylvester + Bertrand-Postulat schränkt stark ein
- **$n = p^a + k$**: p-adische Bewertungsargumente
- **asymptotisch**: Unter GRH für fast alle $(n, k)$ bewiesen

---

## Klassenstruktur: `ErdosSelfridge`

### Attribute

| Attribut | Typ | Beschreibung |
|----------|-----|--------------|
| `_prime_power_cache` | `dict` | Cache für berechnete Primzahlpotenz-Tests |

### Methoden

#### Primzahlpotenz-Prüfung

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `is_binomial_prime_power` | `(n, k) → Optional[(p,a)]` | Gibt $(p, a)$ wenn $\binom{n}{k} = p^a$, sonst `None` |
| `check_conjecture` | `(n, k) → dict` | Vollständige Prüfung mit Gegenbeispiel-Test |

#### Numerische Verifikation

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `verify_up_to` | `(n_max) → dict` | Exhaustive Prüfung bis $n = n_\max$ |
| `find_prime_power_binomials` | `(n_max, k_max) → List[dict]` | Alle $\binom{n}{k}$ die Primzahlpotenzen sind |

#### Sylvesters Theorem

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `sylvester_theorem_check` | `(n, k) → dict` | Verifiziert Primteiler $> k$ |
| `verify_sylvester` | `(n_max) → dict` | Massenverifikation des Theorems |

#### Spezielle Fälle

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `analyze_k1` | `(n_max) → List[dict]` | $\binom{n}{1} = n$ Primzahlpotenz-Fälle |
| `analyze_k2` | `(n_max) → List[dict]` | $\binom{n}{2}$ Analyse |
| `granville_special_cases` | `() → List[dict]` | Dokumentierte Spezialfälle |

#### Statistik

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `prime_factorization_statistics` | `(n_max) → dict` | Anzahl Primteiler von $\binom{n}{k}$ |
| `conjecture_status` | `() → dict` | Aktueller Beweisstand |

---

## Verwendungsbeispiele

```python
from erdos_selfridge import ErdosSelfridge, _is_prime_power, _binomial_coefficient

es = ErdosSelfridge()

# Prüfe C(9,3) = 84 – keine Primzahlpotenz
print(_binomial_coefficient(9, 3))     # 84
print(_is_prime_power(84))             # None (84 = 2²·3·7)
print(es.is_binomial_prime_power(9, 3))  # None ✓

# Verifikation bis n = 200
result = es.verify_up_to(200)
print(result['conjecture_holds'])  # True
print(result['counterexamples_found'])  # 0

# Sylvester Theorem: C(10,3) = 120 hat Primteiler > 3
sylv = es.sylvester_theorem_check(10, 3)
print(sylv['prime_factors_greater_k'])  # [5]

# Statistik: Wie viele C(n,k) sind Primzahlpotenzen?
stats = es.prime_factorization_statistics(n_max=50)
print(stats['one_prime_factor']['fraction'])  # sehr klein

# k=1 Analyse (trivial, außerhalb der Vermutung)
k1_cases = es.analyze_k1(n_max=20)
# Enthält n=4 (4=2²), n=8 (8=2³), etc.
```

---

## Wichtige Beispiele

### $\binom{9}{3} = 84$ – kein Gegenbeispiel

$$\binom{9}{3} = \frac{9 \cdot 8 \cdot 7}{6} = 84 = 2^2 \cdot 3 \cdot 7$$

$84$ hat **drei verschiedene Primteiler** → keine Primzahlpotenz. Die Vermutung hält.

### $\binom{4}{2} = 6$ – kein Gegenbeispiel

$$\binom{4}{2} = 6 = 2 \cdot 3 \quad \text{(keine Primzahlpotenz)}$$

### Sylvester für $\binom{6}{2} = 15$

$$\binom{6}{2} = 15 = 3 \cdot 5$$

$n = 6 \geq 2k = 4$, Primteiler $5 > k = 2$ ✓

---

## Beweisstand

| Aspekt | Status |
|--------|--------|
| Kein Gegenbeispiel für $n \leq 200$ | **NUMERISCH BESTÄTIGT** |
| Sylvester's Theorem ($n \geq 2k$) | **BEWIESEN** (1892) |
| Spezialfälle (Granville) | **TEILWEISE BEWIESEN** |
| Allgemeiner Fall $n \geq k + 2, k \geq 2$ | **OFFEN** – Conjecture |

---

## Wichtige Abhängigkeiten

| Bibliothek | Verwendung |
|-----------|------------|
| `sympy` | `isprime`, `factorint`, `perfect_power` für exakte Faktorisierung |
| `math` | `math.comb` für schnelle Binomialkoeffizienten |

---

## Literatur

- Paul Erdős, John Selfridge: *The product of consecutive integers is never a power*, Illinois J. Math. **19**, 1975
- James Joseph Sylvester: *On arithmetical series*, Messenger of Math. **21**, 1892
- Andrew Granville, Olivier Ramaré: *Explicit bounds on exponential sums and the scarcity of squarefree binomial coefficients*, Mathematika **43**, 1996
