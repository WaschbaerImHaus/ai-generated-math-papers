# artin_primitive_roots.py – Artin-Vermutung über primitive Wurzeln

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Erstellt**: 2026-03-12
**Status**: CONJECTURE – unter GRH bewiesen (Hooley 1967)

---

## Übersicht

Das Modul `artin_primitive_roots.py` untersucht die **Artin-Vermutung über primitive Wurzeln**:

> **CONJECTURE (Emil Artin, 1927)**
> Sei $a \in \mathbb{Z}$ mit $a \neq 0, \pm1$ und $a$ kein vollständiges Quadrat.
> Dann gibt es unendlich viele Primzahlen $p$, für die $a$ eine **primitive Wurzel modulo $p$** ist.
> Die Dichte dieser Primzahlen beträgt die **Artin-Konstante** $C_A \approx 0{,}3739558...$

---

## Mathematischer Hintergrund

### Primitive Wurzel

$a$ ist eine **primitive Wurzel modulo $p$** (Primzahl), wenn:

$$\text{ord}_p(a) = p - 1$$

d.h. $a$ erzeugt die gesamte multiplikative Gruppe $(\mathbb{Z}/p\mathbb{Z})^*$.

**Äquivalente Bedingung**: $a^{(p-1)/q} \not\equiv 1 \pmod{p}$ für jeden Primteiler $q$ von $p-1$.

### Artin-Konstante

$$C_A = \prod_{p \text{ prim}} \left(1 - \frac{1}{p(p-1)}\right) \approx 0{,}3739558136192...$$

Dieses unendliche Produkt konvergiert und gibt die erwartete Dichte der Primzahlen an, für die $a$ primitiv ist.

**Erste Faktoren**:

| $p$ | $1 - \frac{1}{p(p-1)}$ |
|-----|------------------------|
| 2 | $1 - \frac{1}{2} = 0{,}5$ |
| 3 | $1 - \frac{1}{6} \approx 0{,}833$ |
| 5 | $1 - \frac{1}{20} = 0{,}95$ |
| 7 | $1 - \frac{1}{42} \approx 0{,}976$ |

### Hooley's Resultat (1967)

Christopher Hooley bewies 1967 **unter der verallgemeinerten Riemann-Vermutung (GRH)**:

$$\#\{p \leq x : a \text{ primitiv mod } p\} \sim C_A \cdot \frac{x}{\ln x} \quad (x \to \infty)$$

Das ist Artin's Vermutung in quantitativer Form – aber nur bedingt (GRH nicht bewiesen!).

### Heath-Brown (1986) – unbedingtes Resultat

D. R. Heath-Brown bewies **ohne GRH**:

> Mindestens eine der drei Zahlen $\{2, 3, 5\}$ ist primitive Wurzel für unendlich viele Primzahlen $p$.

Dies ist das stärkste bekannte **unbedingte** Resultat zur Artin-Vermutung.

---

## Klassenstruktur: `ArtinPrimitiveRoots`

### Attribute

| Attribut | Typ | Beschreibung |
|----------|-----|--------------|
| `ARTIN_CONSTANT_APPROX` | `float` | $C_A \approx 0{,}3739558136...$ |
| `_prime_cache` | `dict` | Cache für Primzahlsiebe |

### Methoden

#### Primitive Wurzel Prüfen

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `is_primitive_root` | `(a, p) → bool` | Prüft $\text{ord}_p(a) = p-1$ |
| `primitive_root_check_detailed` | `(a, p) → dict` | Detailanalyse mit Faktorisierung |

#### Empirische Dichte

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `empirical_density` | `(a, limit) → dict` | Dichte mit Vergleich zu $C_A$ |
| `count_primitive_root_primes` | `(a, limit) → int` | Anzahl der Primzahlen mit $a$ primitiv |
| `first_primitive_root_primes` | `(a, count) → List[int]` | Erste $n$ Primzahlen |

#### Artin-Konstante

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `artin_constant` | `(num_primes) → float` | Approximation von $C_A$ |
| `artin_constant_convergence` | `(steps) → List[dict]` | Konvergenzanalyse |

#### Hooley & Heath-Brown

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `hooley_density` | `(a) → dict` | GRH-basierte Dichte für Basis $a$ |
| `heath_brown_result` | `() → dict` | Heath-Brown 1986 Theorem |

#### Basis-Analyse

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `analyze_base` | `(a, limit) → dict` | Vollanalyse einer Basis |
| `compare_bases` | `(bases, limit) → List[dict]` | Vergleich mehrerer Basen |
| `conjecture_status` | `() → dict` | Aktueller Beweisstand |

---

## Verwendungsbeispiele

```python
from artin_primitive_roots import ArtinPrimitiveRoots

artin = ArtinPrimitiveRoots()

# Primitive Wurzel prüfen
print(artin.is_primitive_root(2, 5))   # True  (ord_5(2) = 4 = φ(5))
print(artin.is_primitive_root(2, 7))   # False (ord_7(2) = 3 ≠ 6)
print(artin.is_primitive_root(3, 7))   # True  (ord_7(3) = 6 = φ(7))

# Artin-Konstante berechnen
c = artin.artin_constant(num_primes=1000)
print(f"C_A ≈ {c:.6f}")  # 0.373956...

# Empirische Dichte für Basis 2
density = artin.empirical_density(2, limit=10000)
print(f"Dichte: {density['empirical_density']:.4f}")  # ≈ 0.3739

# Erste 10 Primzahlen mit 2 als primitiver Wurzel
primes = artin.first_primitive_root_primes(2, count=10)
print(primes)  # [3, 5, 11, 13, 19, 29, 37, 53, 59, 61]

# Heath-Brown Theorem
hb = artin.heath_brown_result()
print(hb['unconditional'])  # True
```

---

## Bekannte Ergebnisse für Basis 2

Primzahlen $p$ für die 2 eine primitive Wurzel ist:
$3, 5, 11, 13, 19, 29, 37, 53, 59, 61, 67, 83, ...$

(Folge A001122 in der OEIS)

---

## Spezifische Basen $a = 2, 3, 5, 7$

| Basis $a$ | Artin gilt? | Anmerkung |
|-----------|------------|-----------|
| $a = 2$ | Ja | Kein Quadrat, $\neq \pm 1$ |
| $a = 3$ | Ja | Kein Quadrat, $\neq \pm 1$ |
| $a = 4 = 2^2$ | Nein | Vollständiges Quadrat! |
| $a = 5$ | Ja | Kein Quadrat, $\neq \pm 1$ |
| $a = 7$ | Ja | Kein Quadrat, $\neq \pm 1$ |
| $a = -1$ | Nein | Sonderfall $a = -1$ |

---

## Beweisstand

| Aspekt | Status |
|--------|--------|
| Unendlich viele primit. Wurzeln allgemein | **OFFEN** – Conjecture |
| Dichte = $C_A$ unter GRH | **BEWIESEN** (Hooley 1967) |
| Mindestens eine aus $\{2,3,5\}$ unbedingt | **BEWIESEN** (Heath-Brown 1986) |

---

## Wichtige Abhängigkeiten

| Bibliothek | Verwendung |
|-----------|------------|
| `sympy` | `isprime`, `factorint`, `n_order`, `totient` |
| `numpy` | Sieb des Eratosthenes (Boolean-Array) |

---

## Literatur

- Emil Artin: Brieflich an Helmut Hasse, 1927 (ursprüngliche Vermutung)
- Christopher Hooley: *On Artin's conjecture*, J. reine angew. Math. **225**, 1967
- D. R. Heath-Brown: *Artin's conjecture for primitive roots*, Q. J. Math. Oxford **37**, 1986
