# van_der_waerden.py — Dokumentation

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `van_der_waerden.py`

---

## Überblick

Das Modul implementiert die Berechnung und Verifikation von **Van-der-Waerden-Zahlen** $W(k;r)$ sowie das zugrunde liegende **AP-Färbungsproblem** als SAT-Instanz.

**Van der Waerdens Theorem (1927)**: Für alle $k, r \geq 1$ existiert eine endliche Zahl $W(k;r)$, sodass jede $r$-Färbung von $\{1, \ldots, W(k;r)\}$ eine monochromatische arithmetische Progression (AP) der Länge $k$ enthält.

---

## Klassen

### `ArithmeticProgressionColoring`

SAT-Kodierung des AP-Färbungsproblems. Kodiert die Frage:

> Gibt es eine $r$-Färbung von $\{1,\ldots,N\}$, die **keine** monochromatische $k$-term AP enthält?

#### SAT-Variablen

$$x_{i,c} = \text{True} \iff \text{Zahl } i \text{ erhält Farbe } c$$

Interne Variable-ID: $\text{var}(i,c) = (i-1) \cdot r + c$ (1-basiert, pysat-Konvention).

#### CNF-Klauseln

| Typ | Bedeutung | Klausel |
|-----|-----------|---------|
| 1 | Mindestens eine Farbe | $\bigvee_c x_{i,c}$ für alle $i$ |
| 2 | Höchstens eine Farbe | $\neg x_{i,c_1} \vee \neg x_{i,c_2}$ für $c_1 \neq c_2$ |
| 3 | Keine mono AP | $\neg x_{a,c} \vee \neg x_{a+d,c} \vee \ldots$ für alle APs |

#### Methoden

| Methode | Rückgabe | Beschreibung |
|---------|---------|--------------|
| `_var(i, c)` | `int` | SAT-Variablen-ID |
| `_compute_all_aps()` | `List[List[int]]` | Alle $k$-term APs in $\{1,\ldots,N\}$ |
| `build_cnf()` | `CNF` | CNF-Formel aufbauen |
| `solve()` | `Optional[List[int]]` | SAT-Lösung oder `None` |
| `extract_coloring_dict()` | `Optional[Dict[int,int]]` | Färbung als Dict |

---

### `VanDerWaerdenNumbers`

Hauptklasse für Berechnung und Verifikation von $W(k;r)$.

#### Bekannte Werte (Konvention: `(k, r)`)

| $(k, r)$ | $W(k;r)$ | Bedeutung |
|----------|---------|-----------|
| $(3, 2)$ | $9$ | Kleinste N: jede 2-Färbung hat mono 3-AP |
| $(4, 2)$ | $35$ | Kleinste N: jede 2-Färbung hat mono 4-AP |
| $(5, 2)$ | $178$ | Kleinste N: jede 2-Färbung hat mono 5-AP |
| $(3, 3)$ | $27$ | Kleinste N: jede 3-Färbung hat mono 3-AP |
| $(3, 4)$ | $293$ | Kleinste N: jede 4-Färbung hat mono 3-AP |

#### Berlekamp-Untergrenze

$$W(p+1; 2) > p \cdot 2^p \quad \text{für jede Primzahl } p$$

Berlekamp (1968) konstruierte explizit AP-freie 2-Färbungen von $\{1,\ldots,p \cdot 2^p\}$ mittels des Galois-Körpers $\text{GF}(2^p)$.

Beispiele:

| $p$ | Untergrenze $p \cdot 2^p$ | $W(p+1;2)$ |
|----|--------------------------|------------|
| $2$ | $8$ | $W(3;2) = 9 > 8$ ✓ |
| $3$ | $24$ | $W(4;2) = 35 > 24$ ✓ |
| $5$ | $160$ | $W(6;2) > 160$ |

#### Gowers-Obergrenze

$$W(k;2) \leq 2^{2^{k+9}}$$

Diese Schranke (Gowers 2001) ist astronomisch groß und dient nur zur theoretischen Vollständigkeit.

#### Berechnungsalgorithmus (`compute`)

Schrittweise SAT-Suche:

```
für N = k, k+1, k+2, ...:
    falls SAT(N, k, r) = UNSAT:
        W(k;r) = N  ← gefunden
    falls SAT(N, k, r) = SAT:
        weiter suchen
```

#### Methoden

| Methode | Rückgabe | Beschreibung |
|---------|---------|--------------|
| `lookup(k, r)` | `Optional[int]` | Tabellen-Lookup |
| `compute(k, r, max_N)` | `Optional[int]` | SAT-Berechnung |
| `verify_lower_bound(W, k, r)` | `bool` | Untergrenze verifizieren |
| `verify_upper_bound(W, k, r)` | `bool` | Obergrenze verifizieren |
| `find_ap_free_coloring(N, k, r)` | `Optional[List]` | SAT-Zeugen suchen |
| `is_ap_free_coloring(coloring, k)` | `bool` | Färbung prüfen |
| `is_arithmetic_progression(lst)` | `bool` | AP-Test |
| `detect_ap(sequence, k)` | `List[List[int]]` | APs in Folge finden |
| `berlekamp_lower_bound(p)` | `int` | Berlekamp-Schranke |
| `gowers_upper_bound(k)` | `int` | Gowers-Schranke |
| `count_ap_free_colorings(N, k, r)` | `int` | Brute-Force-Zählung |

---

## Abhängigkeiten

- `pysat` (SAT-Solver): `from pysat.solvers import Solver`, `from pysat.formula import CNF`
- `numpy` (numerische Hilfsfunktionen)
- `sympy` (nicht direkt, aber im Projekt verfügbar)

---

## Mathematische Hintergründe

### Van der Waerdens Theorem

**Satz** (van der Waerden 1927): Für alle $k, r \geq 1$ existiert $W(k;r) < \infty$.

**Beweisidee**: Doppelte Induktion über $r$ und $k$. Für feste $r$: Zeige $W(k;r) \leq f(k,r)$ für eine primitive-rekursive Funktion $f$.

### Konventionshinweis

In der Literatur gibt es zwei Notationen:

- **Diese Implementierung**: $W(k;r)$ mit $k =$ AP-Länge, $r =$ Farbenanzahl
- **Alternative**: $W(r;k)$ (umgekehrte Reihenfolge) — z.B. in manchen Tabellen

Die `KNOWN_VALUES`-Tabelle verwendet den Schlüssel `(k, r)`.

---

## Verwendungsbeispiele

```python
from van_der_waerden import VanDerWaerdenNumbers, ArithmeticProgressionColoring

vdw = VanDerWaerdenNumbers()

# Bekannten Wert nachschlagen
print(vdw.lookup(3, 2))  # 9

# W(3;2) via SAT berechnen
print(vdw.compute(3, 2, max_N=15))  # 9

# AP-freie Färbung von {1,...,8} für k=3, r=2 finden
coloring = vdw.find_ap_free_coloring(8, k=3, r=2)
print(coloring)  # z.B. [2, 2, 1, 1, 2, 2, 1, 1]

# Prüfen ob eine Folge eine AP ist
print(vdw.is_arithmetic_progression([1, 3, 5]))  # True
print(vdw.is_arithmetic_progression([1, 2, 4]))  # False

# Berlekamp-Untergrenze für Primzahl p=3
print(vdw.berlekamp_lower_bound(3))  # 24
```
