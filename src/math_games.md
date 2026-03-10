# Modul: Mathematische Spiele (`math_games.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul implementiert kombinatorische Spieltheorie und zahlentheoretische Spiele:
Nim und Sprague-Grundy-Theorie, zahlentheoretische Spiele (Euklid, Teiler, Primzahl,
Fibonacci-Nim, Chomp), magische und lateinische Quadrate sowie perfekte, befreundete
und figurate Zahlen. Der Fokus liegt auf der mathematischen Theorie hinter jedem Spiel.

## Mathematischer Hintergrund

### Nim und Sprague-Grundy-Theorie

**Nim:** $n$ Haufen mit $a_1, \ldots, a_n$ Steinen. Spieler nehmen abwechselnd aus
einem Haufen. Wer den letzten Stein nimmt, gewinnt (Normalspiel).

**Sprague-Grundy-Theorem:**

$$G(\text{Stellung}) = \bigoplus_{i} G(H_i) \quad \text{(Nim-Summe = bitweises XOR)}$$

$G = 0$: Zweitspieler gewinnt; $G \neq 0$: Erststieler gewinnt.

**Grundy-Wert** (Mex = Minimum Excludant):

$$G(n) = \text{mex}\{G(n - m) : m \in \text{Züge}\}, \qquad \text{mex}(S) = \min(\mathbb{N}_0 \setminus S)$$

**Spielsumme:**

$$G(G_1 + G_2) = G(G_1) \oplus G(G_2)$$

### Zahlentheoretische Spiele

**Euklid-Spiel:** Stapel $(a, b)$ mit $a \leq b$. Zug: subtrahiere positives Vielfaches
von $a$ von $b$. Spieler 1 gewinnt genau dann wenn $b \geq 2a$ (Paritätskontrolle).

**Divisor-Spiel:** Starte mit $n$, subtrahiere echten Teiler $d < n$. Gerade $n$:
Spieler 1 gewinnt (nimmt 1, den einzigen ungeraden echten Teiler).

**Fibonacci-Nim (Zeckendorf-Theorem):** Erststieler gewinnt genau dann wenn $n$
keine Fibonacci-Zahl ist. Gewinnstrategie: nimm den kleinsten Zeckendorf-Term.

**Chomp (Strategie-Stehlen):** Für $m, n \geq 2$ hat Spieler 1 immer eine Gewinnstrategie
(existentiell, nicht konstruktiv für allgemeine $m, n$):
> Annahme Spieler 2 hat Strategie $\sigma$. Spieler 1 wählt obere rechte Ecke →
> die resultierende Stellung ist für $\sigma$ mindestens so gut → Widerspruch.

### Magische Quadrate

**Magische Summe** eines $n \times n$-Quadrats mit Einträgen $1, \ldots, n^2$:

$$S(n) = \frac{n(n^2 + 1)}{2}$$

Lo-Shu (ältestes magisches Quadrat, $\approx 2800$ v. Chr.):

$$\begin{pmatrix} 2 & 7 & 6 \\ 9 & 5 & 1 \\ 4 & 3 & 8 \end{pmatrix}, \quad S(3) = 15$$

### Figurate Zahlen

| Name | Formel | Folge |
|------|--------|-------|
| Dreieckszahl | $T_n = n(n+1)/2$ | 1, 3, 6, 10, 15, … |
| Quadratzahl | $S_n = n^2$ | 1, 4, 9, 16, 25, … |
| Pentagonalzahl | $P_n = n(3n-1)/2$ | 1, 5, 12, 22, 35, … |
| Hexagonalzahl | $H_n = n(2n-1)$ | 1, 6, 15, 28, 45, … |

Allgemeine $s$-eckige Zahl:

$$P(s, n) = \frac{n\left[(s-2)n - (s-4)\right]}{2}$$

**Dreieckszahl-Kriterium:** $n = T_k \iff 8n + 1$ ist ein perfektes Quadrat.

### Perfekte und befreundete Zahlen

**Perfekte Zahl:** $\sigma(n) = 2n$ (Summe aller Teiler). Bekannte Werte: 6, 28, 496, 8128, …

$$\text{perfekt: } \sigma(n) = 2n, \quad \text{abundant: } \sigma(n) > 2n, \quad \text{defizient: } \sigma(n) < 2n$$

**Befreundete Zahlen:** $\sigma(a) - a = b$ und $\sigma(b) - b = a$.
Bekannte Paare: $(220, 284)$, $(1184, 1210)$, $(2620, 2924)$.

**Kanonenkugelproblem:** $\sum_{i=1}^n i^2 = k^2$ hat nur Lösungen $n=1$ und $n=24$
($S_{24} = 4900 = 70^2$).

---

## Klassen und Methoden

### `NimGame`

**Beschreibung:** Implementiert Nim und grundlegende Sprague-Grundy-Theorie. Berechnet
den Nim-Wert (XOR aller Haufengrößen), den optimalen Gewinnzug und allgemeine
Grundy-Werte für beliebige Zugmengen.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `nim_value()` | `(piles: list) -> int` | $G = a_1 \oplus a_2 \oplus \ldots \oplus a_n$ |
| `nim_winning_move()` | `(piles: list) -> Optional[tuple[int,int]]` | Gewinnzug: `(Haufen-Index, neue Größe)` oder `None` |
| `grundy_value()` | `(n: int, moves: set) -> int` | $G(n) = \text{mex}\{G(n-m) : m \in \text{Züge}\}$ (Bottom-Up) |
| `mex()` | `(values: set) -> int` | Minimum Excludant: kleinste nicht-negative Zahl $\notin$ `values` |
| `nim_sum()` | `(values: list) -> int` | Bitweises XOR aller Werte |

---

### `SpragueGrundyTheory`

**Beschreibung:** Verallgemeinerte Sprague-Grundy-Theorie für beliebige kombinatorische
Spiele. Berechnet Grundy-Werte via Memoization, Kayles-Werte und Subtraktionsspiele.

| Methode | Signatur | Beschreibung |
|---------|----------|-------------|
| `analyze_game()` | `(position: int, move_fn: Callable, max_depth: int = 20) -> int` | Grundy-Wert via rekursiver mex-Berechnung mit Cache |
| `game_sum_grundy()` | `(g1: int, g2: int) -> int` | $G(G_1 + G_2) = g_1 \oplus g_2$ |
| `kayles_grundy()` | `(n: int) -> list` | Kayles-Grundy-Werte $G(0), \ldots, G(n)$ (Züge: 1 oder 2 Kegel) |
| `subtraction_game_grundy()` | `(n: int, subtraction_set: set) -> list` | Grundy-Werte für Subtraktionsspiel |
| `take_away_game()` | `(n: int, moves: list) -> list` | Wrapper für `subtraction_game_grundy()` |

---

### Standalone-Funktionen: Zahlentheoretische Spiele

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `euclid_game()` | `(a: int, b: int) -> dict` | Euklid-Spiel: Gewinner + Spielverlauf-Simulation |
| `divisor_game()` | `(n: int) -> dict` | Teiler-Spiel: Grundy-Wert + Gewinner |
| `prime_game()` | `(n: int) -> dict` | Primzahlspiel: erlaubte Züge = Primzahlen ≤ $n$ |
| `fibonacci_nim()` | `(n: int) -> dict` | Fibonacci-Nim: Zeckendorf-Zerlegung + Strategie |
| `chomp_game_demo()` | `(m: int, n: int) -> dict` | Chomp: Strategie-Stehlen-Argument + bekannte Fälle |

---

### Standalone-Funktionen: Magische Quadrate

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `magic_square_check()` | `(matrix: list) -> bool` | Prüft Zeilen-, Spalten- und Diagonalensummen |
| `magic_square_3x3()` | `(start: int = 1) -> list` | Lo-Shu-Quadrat mit Startwert `start` |
| `magic_square_sum()` | `(n: int) -> int` | Magische Summe $S(n) = n(n^2+1)/2$ |
| `latin_square_check()` | `(matrix: list) -> bool` | Prüft: jede Zeile/Spalte enthält alle $n$ Werte genau einmal |

---

### Standalone-Funktionen: Perfekte, befreundete Zahlen

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `perfect_number_check()` | `(n: int) -> bool` | $\sigma(n) = 2n$ |
| `abundant_deficient_classify()` | `(n: int) -> str` | `'perfekt'`, `'abundant'` oder `'defizient'` |
| `friendly_numbers()` | `(n: int, max_search: int = 1000) -> list` | Sucht befreundete Zahlen: $\sigma(a)-a = b$, $\sigma(b)-b = a$ |

---

### Standalone-Funktionen: Figurate Zahlen

| Funktion | Signatur | Beschreibung |
|----------|----------|-------------|
| `triangular_number()` | `(n: int) -> int` | $T_n = n(n+1)/2$ |
| `square_number()` | `(n: int) -> int` | $S_n = n^2$ |
| `pentagonal_number()` | `(n: int) -> int` | $P_n = n(3n-1)/2$ |
| `hexagonal_number()` | `(n: int) -> int` | $H_n = n(2n-1)$ |
| `is_triangular()` | `(n: int) -> bool` | $n = T_k \iff 8n+1$ perfektes Quadrat |
| `polygonal_number()` | `(s: int, n: int) -> int` | $P(s,n) = n[(s-2)n-(s-4)]/2$ ($s$-eckige Zahl) |
| `cannonball_problem()` | `(n: int) -> dict` | $\sum_{i=1}^k i^2 = m^2$? Bekannte Lösungen $k=1$ und $k=24$ |

---

## Beispiele

```python
from math_games import (
    NimGame, SpragueGrundyTheory,
    euclid_game, divisor_game, fibonacci_nim, chomp_game_demo,
    magic_square_3x3, magic_square_check, magic_square_sum,
    perfect_number_check, friendly_numbers,
    triangular_number, polygonal_number, cannonball_problem
)

# Nim
piles = [3, 5, 7]
print(NimGame.nim_value(piles))           # 3 XOR 5 XOR 7 = 1 ≠ 0 → Spieler 1 gewinnt
print(NimGame.nim_winning_move(piles))    # z.B. (2, 6): Haufen 2 von 7 auf 6 setzen

# Grundy-Werte für Subtraktionsspiel {1, 2, 3}
g = NimGame.grundy_value(10, moves={1, 2, 3})
print(g)   # periodisch: 0,1,2,3,0,1,2,3,...

# Sprague-Grundy für allgemeines Spiel
sg = SpragueGrundyTheory()
move_fn = lambda n: [n - 1, n - 2] if n >= 2 else ([n - 1] if n == 1 else [])
print(sg.analyze_game(8, move_fn))  # Grundy-Wert für "1 oder 2 nehmen" bei n=8

# Euklid-Spiel
result = euclid_game(14, 21)
print(result['gewinner'])   # Spieler 1 (21 ≥ 2·14 → Nein; aber Analyse gibt Antwort)

# Fibonacci-Nim
print(fibonacci_nim(8))   # n=8 ist Fibonacci → Spieler 2 gewinnt
print(fibonacci_nim(10))  # 10=8+2 (Zeckendorf) → Spieler 1 gewinnt

# Chomp
print(chomp_game_demo(4, 5)['gewinner'])  # Spieler 1 (Strategie-Stehlen)

# Magisches Quadrat
sq = magic_square_3x3()
print(sq)                         # [[2,7,6],[9,5,1],[4,3,8]]
print(magic_square_check(sq))     # True
print(magic_square_sum(3))        # 15

# Perfekte Zahlen
print(perfect_number_check(28))   # True (1+2+4+7+14=28)
print(perfect_number_check(12))   # False (abundant)

# Befreundete Zahlen
print(friendly_numbers(220))     # [284]

# Figurate Zahlen
print(triangular_number(10))     # 55
print(polygonal_number(6, 5))    # 45 (5. Hexagonalzahl)
print(cannonball_problem(30)['loesungen'])
# [{'n':1,'summe':1,'quadratwurzel':1}, {'n':24,'summe':4900,'quadratwurzel':70}]
```

## Verweise

- Verwandte Module: `algebra.py`, `proof_theory.py`, `analytic_number_theory.py`
- Literatur:
  - Berlekamp, Conway, Guy: *Winning Ways for your Mathematical Plays* (A K Peters, 2001)
  - Ferguson: *Game Theory* (Academic Press, 1995)
  - Hardy & Wright: *An Introduction to the Theory of Numbers* (Oxford, 2008)
  - Conway & Guy: *The Book of Numbers* (Copernicus, 1996)
  - Gardner: *Mathematical Games* (Scientific American, div. Jahrgänge)
