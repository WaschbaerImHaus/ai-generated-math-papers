"""
math_games.py – Mathematische Spiele

Dieses Modul implementiert kombinatorische Spieltheorie und zahlentheoretische Spiele:
- Nim und Sprague-Grundy-Theorie
- Zahlentheoretische Spiele (Euklid, Teiler, Primzahl, Fibonacci-Nim, Chomp)
- Magische Quadrate und lateinische Quadrate
- Perfekte, befreundete und figurate Zahlen

@author    Michael Fuhrmann
@version   1.0.0
@timestamp 2026-03-10
"""

import math
import functools
from typing import Optional, Callable, Set, List, Tuple


# ---------------------------------------------------------------------------
# 1. NimGame – Nim und Sprague-Grundy
# ---------------------------------------------------------------------------

class NimGame:
    """
    Implementierung des Nim-Spiels und grundlegender Sprague-Grundy-Theorie.

    Nim: Zwei Spieler wechseln sich ab, aus einem oder mehreren Haufen Steine
    zu nehmen. Wer den letzten Stein nimmt, gewinnt (Normalspiel-Konvention).

    @author    Michael Fuhrmann
    @version   1.0.0
    @timestamp 2026-03-10
    """

    @staticmethod
    def nim_value(piles: list) -> int:
        """
        Berechnet den Nim-Grundwert (Nim-Wert) einer Stellung:
        G = a₁ XOR a₂ XOR … XOR aₙ

        G = 0 → Zweitspieler gewinnt (Verlierer-Position für Erststieler).
        G ≠ 0 → Erststieler gewinnt.

        @param piles Liste der Haufengrößen (nicht-negativ)
        @return      Nim-Wert (XOR aller Haufengrößen)

        @timestamp 2026-03-10
        """
        result = 0
        for pile in piles:
            result ^= pile  # bitweises XOR (Nim-Summe)
        return result

    @staticmethod
    def nim_winning_move(piles: list) -> Optional[Tuple[int, int]]:
        """
        Berechnet den Gewinnzug in einer Nim-Stellung.

        Methode: Finde einen Haufen i und eine Zielgröße t < piles[i],
        sodass nach dem Zug die Nim-Summe = 0 wird.

        @param piles Liste der Haufengrößen
        @return      (Haufen-Index, neue Haufengröße) oder None (wenn Verlierposition)

        @timestamp 2026-03-10
        """
        nim_sum = NimGame.nim_value(piles)
        if nim_sum == 0:
            return None  # Verlierposition – kein Gewinnzug vorhanden

        for i, pile in enumerate(piles):
            # Zielgröße: XOR des Haufens mit der Nim-Summe
            target = pile ^ nim_sum
            if target < pile:
                # Reduziere Haufen i von pile auf target (target < pile → gültiger Zug)
                return (i, target)
        return None  # Sollte nie eintreten wenn nim_sum ≠ 0

    @staticmethod
    def grundy_value(n: int, moves: set) -> int:
        """
        Berechnet den Grundy-Wert G(n) für ein Spiel mit Startzustand n
        und erlaubten Zugmengen moves.

        G(n) = mex{G(n - m) : m ∈ moves, m ≤ n}

        @param n     Aktueller Spielzustand (nicht-negativ)
        @param moves Menge erlaubter Züge (positive ganze Zahlen)
        @return      Grundy-Wert G(n)

        @timestamp 2026-03-10
        """
        # Tabellarische Berechnung (Bottom-Up) bis n
        g = [0] * (n + 1)
        for i in range(1, n + 1):
            reachable = set()
            for m in moves:
                if m <= i:
                    reachable.add(g[i - m])
            g[i] = NimGame.mex(reachable)
        return g[n]

    @staticmethod
    def mex(values: set) -> int:
        """
        Minimum Excludant (mex): Kleinste nicht-negative ganze Zahl,
        die nicht in `values` enthalten ist.

        mex({0,1,3}) = 2

        @param values Menge nicht-negativer ganzer Zahlen
        @return       Minimum Excludant

        @timestamp 2026-03-10
        """
        m = 0
        while m in values:
            m += 1
        return m

    @staticmethod
    def nim_sum(values: list) -> int:
        """
        Berechnet die Nim-Summe (bitweises XOR) einer Liste von Werten.

        @param values Liste von nicht-negativen ganzen Zahlen
        @return       Nim-Summe (XOR aller Werte)

        @timestamp 2026-03-10
        """
        result = 0
        for v in values:
            result ^= v
        return result


# ---------------------------------------------------------------------------
# 2. SpragueGrundyTheory – Sprague-Grundy-Theorie
# ---------------------------------------------------------------------------

class SpragueGrundyTheory:
    """
    Implementierung der Sprague-Grundy-Theorie für kombinatorische Spiele.

    Jedes endliche, deterministische Zwei-Spieler-Spiel ohne Zufallskomponente
    (normal convention) ist äquivalent zu einem Nim-Haufen der Größe G(Stellung).

    @author    Michael Fuhrmann
    @version   1.0.0
    @timestamp 2026-03-10
    """

    @staticmethod
    def analyze_game(position: int, move_fn: Callable, max_depth: int = 20) -> int:
        """
        Berechnet den Grundy-Wert einer Spielstellung rekursiv.

        move_fn(pos) → List[int]: Gibt alle erreichbaren Nachfolgestellungen zurück.

        @param position  Aktuelle Spielstellung (nicht-negative ganze Zahl)
        @param move_fn   Funktion position → [nachfolger1, nachfolger2, ...]
        @param max_depth Maximale Rekursionstiefe
        @return          Grundy-Wert

        @timestamp 2026-03-10
        """
        # Memoization-Cache für Effizienz
        cache = {}

        def _grundy(pos: int, depth: int) -> int:
            if pos in cache:
                return cache[pos]
            if depth == 0:
                return 0
            successors = move_fn(pos)
            if not successors:
                cache[pos] = 0  # Terminale Stellung (kein Zug möglich) → Grundy = 0
                return 0
            reachable = {_grundy(s, depth - 1) for s in successors}
            g = NimGame.mex(reachable)
            cache[pos] = g
            return g

        return _grundy(position, max_depth)

    @staticmethod
    def game_sum_grundy(g1: int, g2: int) -> int:
        """
        Berechnet den Grundy-Wert der Summe zweier Spiele:
        G(G₁ + G₂) = G(G₁) XOR G(G₂)

        @param g1 Grundy-Wert des ersten Spiels
        @param g2 Grundy-Wert des zweiten Spiels
        @return   Grundy-Wert der Spielsumme

        @timestamp 2026-03-10
        """
        return g1 ^ g2  # Nim-Summe der Grundy-Werte

    @staticmethod
    def kayles_grundy(n: int) -> list:
        """
        Berechnet die Grundy-Werte für Kayles bis Länge n.

        Kayles: Zwei Spieler wechseln sich ab, einen oder zwei benachbarte
        Kegel umzuwerfen (aus einer Reihe von n Kegeln).
        Erlaubte Züge: m ∈ {1, 2} (1 oder 2 Kegel entfernen).

        @param n Maximale Reihenlänge
        @return  Liste der Grundy-Werte G(0), G(1), ..., G(n)

        @timestamp 2026-03-10
        """
        g = [0] * (n + 1)
        for i in range(1, n + 1):
            reachable = set()
            if i >= 1:
                reachable.add(g[i - 1])  # 1 Kegel entfernen
            if i >= 2:
                reachable.add(g[i - 2])  # 2 Kegel entfernen
            g[i] = NimGame.mex(reachable)
        return g

    @staticmethod
    def subtraction_game_grundy(n: int, subtraction_set: set) -> list:
        """
        Berechnet Grundy-Werte für ein Subtraktionsspiel bis n.

        Subtraktionsspiel: Erlaubt ist, genau s Steine zu entfernen, wobei
        s ∈ subtraction_set. Wer den letzten Stein nimmt, gewinnt.

        @param n               Maximale Steinanzahl
        @param subtraction_set Erlaubte Subtraktionsmengen (positive ganze Zahlen)
        @return                Liste der Grundy-Werte G(0), ..., G(n)

        @timestamp 2026-03-10
        """
        g = [0] * (n + 1)
        for i in range(1, n + 1):
            reachable = set()
            for s in subtraction_set:
                if s <= i:
                    reachable.add(g[i - s])
            g[i] = NimGame.mex(reachable)
        return g

    @staticmethod
    def take_away_game(n: int, moves: list) -> list:
        """
        Grundy-Werte für ein Take-Away-Spiel: Aus einem Haufen von n Steinen
        darf man genau m ∈ moves Steine nehmen.

        @param n     Maximale Steinanzahl
        @param moves Liste erlaubter Zuggrößen
        @return      Grundy-Werte G(0), ..., G(n)

        @timestamp 2026-03-10
        """
        return SpragueGrundyTheory.subtraction_game_grundy(n, set(moves))


# ---------------------------------------------------------------------------
# 3. Zahlentheoretische Spiele
# ---------------------------------------------------------------------------

def euclid_game(a: int, b: int) -> dict:
    """
    Euklid-Spiel: Zwei Spieler haben zwei Stapel (a, b) mit a ≤ b.
    Zug: Subtrahiere ein positives Vielfaches des kleineren vom größeren.
    Wer den ersten Spieler in eine (0,0)-Stellung zwingen kann, gewinnt.

    Analyse: Wenn b ≥ 2a, kann Spieler 1 die Parität der Züge steuern → gewinnt.

    @param a Erste Zahl (kleiner oder gleich b)
    @param b Zweite Zahl
    @return  Wörterbuch mit Gewinner und Analyse

    @timestamp 2026-03-10
    """
    if a > b:
        a, b = b, a  # Sicherstellen a ≤ b

    moves = []
    pos = (a, b)
    first_player_wins = _euclid_game_wins(a, b)

    # Simuliere optimalen Spielverlauf
    current_a, current_b = a, b
    while current_a > 0 and current_b > 0:
        if current_a > current_b:
            current_a, current_b = current_b, current_a
        # Optimaler Zug: wähle k so, dass Nachfolger Verlierstellung ist
        q = current_b // current_a
        # Wenn b >= 2a, kann man entweder k=q oder k=q-1 wählen
        for k in [q, q - 1]:
            if k >= 1:
                new_b = current_b - k * current_a
                if new_b >= 0:
                    moves.append((current_a, current_b, k, (min(current_a, new_b), max(current_a, new_b))))
                    current_b = new_b
                    if current_a > current_b:
                        current_a, current_b = current_b, current_a
                    break
        else:
            break

    return {
        'a': a, 'b': b,
        'gewinner': 'Spieler 1' if first_player_wins else 'Spieler 2',
        'erklaerung': (
            f"b={b} >= 2a={2*a}: Spieler 1 gewinnt durch Paritätskontrolle."
            if b >= 2 * a else
            f"b={b} < 2a={2*a}: Eindeutiger Zug, rekursiv analysiert."
        ),
        'züge': moves[:10]
    }


def _euclid_game_wins(a: int, b: int) -> bool:
    """
    Iterative Hilfsfunktion: Bestimmt, ob der aktuelle Spieler im Euklid-Spiel gewinnt.

    Algorithmus (iterativ, kein Rekursionsproblem):
    - Falls b >= 2a: aktueller Spieler gewinnt (Paritätskontrolle möglich)
    - Falls b < 2a: einzig möglicher Zug ist b → b - a (Spieler wechselt)
      Wir verfolgen, ob die Anzahl der erzwungenen Züge gerade/ungerade ist.

    @param a Kleinere Zahl (a <= b)
    @param b Größere Zahl
    @return  True wenn der aktuelle Spieler gewinnt

    @timestamp 2026-03-10
    """
    # Sicherstellen: a <= b
    if a > b:
        a, b = b, a

    # Laufvariable: True = aktueller Spieler gewinnt
    current_wins = True  # Wird bei jedem erzwungenen Schritt invertiert

    while a > 0:
        if b >= 2 * a:
            # Aktueller Spieler kann Parität steuern → gewinnt wenn current_wins True
            return current_wins
        # Einzig möglicher Zug: b - a (Spieler wechselt → invertiere Perspektive)
        b, a = a, b - a
        if a > b:
            a, b = b, a
        current_wins = not current_wins  # Spieler wechselt

    # a == 0: Wer hier ist, hat verloren (kein Zug möglich)
    return not current_wins


def divisor_game(n: int) -> dict:
    """
    Teiler-Spiel: Zwei Spieler wechseln sich ab, einen echten Teiler d < n
    des aktuellen Werts zu subtrahieren. Wer auf 0 reduziert, gewinnt.

    Analyse: Wenn n gerade, gewinnt Spieler 1; wenn n ungerade, gewinnt Spieler 2.

    @param n Startzahl (positive ganze Zahl)
    @return  Wörterbuch mit Gewinner und Analyse

    @timestamp 2026-03-10
    """
    # Berechne Grundy-Werte für alle k von 0 bis n
    grundy = [0] * (n + 1)
    for k in range(2, n + 1):
        reachable = set()
        # Alle echten Teiler d < k
        for d in range(1, k):
            if k % d == 0:
                reachable.add(grundy[k - d])
        grundy[k] = NimGame.mex(reachable)

    first_wins = grundy[n] > 0
    return {
        'n': n,
        'grundy_wert': grundy[n],
        'gewinner': 'Spieler 1' if first_wins else 'Spieler 2',
        'erklaerung': (
            f"n={n} gerade → Spieler 1 nimmt 1 (einziger ungerader echter Teiler von n)."
            if n % 2 == 0 and n > 1 else
            f"n={n} ungerade → alle echten Teiler sind ungerade → gerade Situation."
        )
    }


def prime_game(n: int) -> dict:
    """
    Primzahlspiel: Zwei Spieler wechseln sich ab, eine Primzahl p ≤ n
    vom aktuellen Wert zu subtrahieren. Wer auf 0 reduziert, gewinnt.

    Analyse über Grundy-Werte.

    @param n Startzahl (positive ganze Zahl)
    @return  Wörterbuch mit Grundy-Wert und Gewinner

    @timestamp 2026-03-10
    """
    def _is_prime(k: int) -> bool:
        if k < 2:
            return False
        if k == 2:
            return True
        if k % 2 == 0:
            return False
        return all(k % i != 0 for i in range(3, int(math.sqrt(k)) + 1, 2))

    # Alle Primzahlen bis n als erlaubte Züge
    primes = {p for p in range(2, n + 1) if _is_prime(p)}

    # Grundy-Werte berechnen
    grundy = [0] * (n + 1)
    for k in range(1, n + 1):
        reachable = set()
        for p in primes:
            if p <= k:
                reachable.add(grundy[k - p])
        grundy[k] = NimGame.mex(reachable)

    return {
        'n': n,
        'grundy_wert': grundy[n],
        'gewinner': 'Spieler 1' if grundy[n] > 0 else 'Spieler 2',
        'erlaubte_zuege': sorted(primes),
        'grundy_tabelle': grundy
    }


def fibonacci_nim(n: int) -> dict:
    """
    Fibonacci-Nim: Im ersten Zug darf man 1 bis n-1 Steine nehmen.
    Danach darf man maximal doppelt so viele wie der Vorgänger nehmen.
    Wer den letzten Stein nimmt, gewinnt.

    Satz (Zeckendorf): Erststieler gewinnt genau dann, wenn n keine Fibonacci-Zahl ist.

    @param n Startanzahl der Steine
    @return  Wörterbuch mit Analyse

    @timestamp 2026-03-10
    """
    # Fibonacci-Zahlen bis n berechnen
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])
    fib_set = set(fibs)

    is_fibonacci = n in fib_set

    # Zeckendorf-Zerlegung: Schreibe n als Summe nicht-aufeinanderfolgender Fibonacci-Zahlen
    def zeckendorf(k: int) -> list:
        """Zeckendorf-Zerlegung von k."""
        fib_list = sorted([f for f in fib_set if f <= k], reverse=True)
        result = []
        remaining = k
        for f in fib_list:
            if f <= remaining:
                result.append(f)
                remaining -= f
        return result

    return {
        'n': n,
        'ist_fibonacci': is_fibonacci,
        'gewinner': 'Spieler 2' if is_fibonacci else 'Spieler 1',
        'erklaerung': (
            f"n={n} ist Fibonacci-Zahl → Spieler 2 gewinnt."
            if is_fibonacci else
            f"n={n} ist keine Fibonacci-Zahl → Spieler 1 gewinnt."
        ),
        'zeckendorf_zerlegung': zeckendorf(n),
        'strategie': (
            "Spieler 1 nimmt den kleinsten Zeckendorf-Term; dann ist die verbleibende "
            "Zahl eine Fibonacci-Zahl (Verlierposition für Spieler 2)."
            if not is_fibonacci else
            "Jeder Zug von Spieler 1 hinterlässt eine Nicht-Fibonacci-Zahl "
            "(Gewinnerposition für Spieler 2)."
        )
    }


def chomp_game_demo(m: int, n: int) -> dict:
    """
    Chomp (Schokoladenspiel) Demonstration:
    Eine m×n-Schokoladentafel: Spieler nehmen abwechselnd ein Quadrat (i,j)
    und alles rechts davon und darüber. Das vergiftete Feld (0,0) unten links
    darf nicht genommen werden. Wer (0,0) nehmen muss, verliert.

    Existenzsatz: Spieler 1 hat für m,n ≥ 2 immer eine Gewinnstrategie
    (Strategie-Stehlen-Argument).

    @param m Zeilen der Schokoladentafel
    @param n Spalten der Schokoladentafel
    @return  Wörterbuch mit Erklärung

    @timestamp 2026-03-10
    """
    return {
        'm': m, 'n': n,
        'gewinner': (
            'Spieler 2' if (m == 1 and n == 1)
            else 'Spieler 1'
        ),
        'erklaerung': (
            "1×1 Tafel: Spieler 1 muss (0,0) nehmen → verliert."
            if m == 1 and n == 1 else
            "Strategie-Stehlen: Da (m,n) ≠ (1,1), hat Spieler 1 eine Gewinnstrategie "
            "(existentiell, nicht konstruktiv für allgemeine m,n)."
        ),
        'strategie_stehlen': (
            "Annahme: Spieler 2 hat Gewinnstrategie σ. Spieler 1 wählt (m-1,n-1) "
            "(obere rechte Ecke). Der verbleibende Zustand ist für σ mindestens so gut "
            "→ Spieler 1 kann σ imitieren → Widerspruch."
        ),
        'explizit_bekannt': {
            '2x2': 'Spieler 1 nimmt (1,1) → 2×1-Tafel → Spieler 1 gewinnt',
            '3x3': 'Erststieler gewinnt, Gewinnzug bekannt',
            'nx1': 'Spieler 1 nimmt alle außer (0,0): Spieler 2 muss (0,0) nehmen',
        }
    }


# ---------------------------------------------------------------------------
# 4. Magische Quadrate und Zahlenrätsel
# ---------------------------------------------------------------------------

def magic_square_check(matrix: list) -> bool:
    """
    Prüft, ob eine quadratische Matrix ein magisches Quadrat ist.

    Ein n×n-magisches Quadrat enthält die Zahlen 1 bis n² genau einmal,
    und alle Zeilen-, Spalten- und Diagonalensummen sind gleich.

    @param matrix 2D-Liste (n×n) mit ganzen Zahlen
    @return       True wenn magisches Quadrat

    @timestamp 2026-03-10
    """
    n = len(matrix)
    if n == 0:
        return False
    # Quadratisch prüfen
    if any(len(row) != n for row in matrix):
        return False

    # Magische Summe berechnen
    magic_sum = magic_square_sum(n)

    # Zeilensummen prüfen
    for row in matrix:
        if sum(row) != magic_sum:
            return False

    # Spaltensummen prüfen
    for j in range(n):
        if sum(matrix[i][j] for i in range(n)) != magic_sum:
            return False

    # Hauptdiagonale
    if sum(matrix[i][i] for i in range(n)) != magic_sum:
        return False

    # Nebendiagonale
    if sum(matrix[i][n - 1 - i] for i in range(n)) != magic_sum:
        return False

    return True


def magic_square_3x3(start: int = 1) -> list:
    """
    Erzeugt das 3×3-magische Quadrat (Lo Shu, ~2800 v. Chr.) mit Einträgen
    start bis start+8.

    Standard-Anordnung (start=1):
    ```
    2  7  6
    9  5  1
    4  3  8
    ```
    Magische Summe = 15 (für start=1).

    @param start Startwert (Standard 1, dann Einträge 1-9)
    @return      3×3-Matrix als 2D-Liste

    @timestamp 2026-03-10
    """
    # Lo-Shu-Quadrat (normiert auf 1..9, dann verschoben)
    base = [
        [2, 7, 6],
        [9, 5, 1],
        [4, 3, 8]
    ]
    offset = start - 1
    return [[base[i][j] + offset for j in range(3)] for i in range(3)]


def magic_square_sum(n: int) -> int:
    """
    Berechnet die magische Summe eines n×n-Quadrats mit Einträgen 1 bis n².

    S(n) = n(n² + 1) / 2

    @param n Ordnung des magischen Quadrats
    @return  Magische Summe

    @timestamp 2026-03-10
    """
    return n * (n * n + 1) // 2


def latin_square_check(matrix: list) -> bool:
    """
    Prüft, ob eine Matrix ein lateinisches Quadrat ist.

    Ein n×n-lateinisches Quadrat enthält in jeder Zeile und Spalte
    genau die Elemente {0,1,…,n-1} (oder ein beliebiges n-elementiges Alphabet).

    @param matrix 2D-Liste (n×n)
    @return       True wenn lateinisches Quadrat

    @timestamp 2026-03-10
    """
    n = len(matrix)
    if n == 0:
        return False
    if any(len(row) != n for row in matrix):
        return False

    # Alle verwendeten Werte sammeln
    all_values = set()
    for row in matrix:
        all_values.update(row)

    if len(all_values) != n:
        return False

    # Jede Zeile enthält alle Werte genau einmal
    for row in matrix:
        if set(row) != all_values:
            return False

    # Jede Spalte enthält alle Werte genau einmal
    for j in range(n):
        col = {matrix[i][j] for i in range(n)}
        if col != all_values:
            return False

    return True


def perfect_number_check(n: int) -> bool:
    """
    Prüft, ob n eine perfekte Zahl ist: σ(n) = 2n, d.h. die Summe aller
    positiven Teiler gleich 2n.

    Bekannte perfekte Zahlen: 6, 28, 496, 8128, 33550336, ...

    @param n Positive ganze Zahl
    @return  True wenn n perfekt

    @timestamp 2026-03-10
    """
    if n <= 1:
        return False
    divisor_sum = sum(d for d in range(1, n + 1) if n % d == 0)
    return divisor_sum == 2 * n


def abundant_deficient_classify(n: int) -> str:
    """
    Klassifiziert n als:
    - 'perfekt'    : σ(n) = 2n  (z.B. 6, 28, 496)
    - 'abundant'   : σ(n) > 2n  (z.B. 12, 18, 20)
    - 'defizient'  : σ(n) < 2n  (z.B. 1, 2, 3, 4, 5)

    @param n Positive ganze Zahl
    @return  Klassifikationsstring

    @timestamp 2026-03-10
    """
    if n <= 0:
        raise ValueError("n muss positiv sein.")
    divisor_sum = sum(d for d in range(1, n + 1) if n % d == 0)
    if divisor_sum == 2 * n:
        return 'perfekt'
    elif divisor_sum > 2 * n:
        return 'abundant'
    else:
        return 'defizient'


def friendly_numbers(n: int, max_search: int = 1000) -> list:
    """
    Sucht alle befreundeten Zahlen bis max_search, die mit n befreundet sind.

    Zahlen a und b heißen befreundet (amicable), wenn σ(a) - a = b und σ(b) - b = a,
    d.h. die Summe der echten Teiler von a gleich b und umgekehrt.

    Bekannte Paare: (220, 284), (1184, 1210), (2620, 2924), ...

    @param n          Ausgangszahl
    @param max_search Maximale Suche bis zu dieser Grenze
    @return           Liste aller mit n befreundeten Zahlen

    @timestamp 2026-03-10
    """
    def proper_divisor_sum(k: int) -> int:
        """Summe der echten Teiler von k (ohne k selbst)."""
        return sum(d for d in range(1, k) if k % d == 0)

    friends = []
    s_n = proper_divisor_sum(n)

    # Direkte Prüfung: s(n) = m und s(m) = n
    if 1 < s_n <= max_search and s_n != n:
        if proper_divisor_sum(s_n) == n:
            friends.append(s_n)

    # Erweiterte Suche bis max_search
    for m in range(2, max_search + 1):
        if m == n:
            continue
        s_m = proper_divisor_sum(m)
        if s_m == n and proper_divisor_sum(n) == m:
            if m not in friends:
                friends.append(m)

    return sorted(friends)


# ---------------------------------------------------------------------------
# 5. Figurate Zahlen
# ---------------------------------------------------------------------------

def triangular_number(n: int) -> int:
    """
    Berechnet die n-te Dreieckszahl T_n = n(n+1)/2.

    Darstellung: Punkte, die ein gleichseitiges Dreieck bilden.
    T_1=1, T_2=3, T_3=6, T_4=10, ...

    @param n Positive ganze Zahl
    @return  n-te Dreieckszahl

    @timestamp 2026-03-10
    """
    if n < 0:
        raise ValueError("n muss nicht-negativ sein.")
    return n * (n + 1) // 2


def square_number(n: int) -> int:
    """
    Berechnet die n-te Quadratzahl S_n = n².

    @param n Nicht-negative ganze Zahl
    @return  n-te Quadratzahl

    @timestamp 2026-03-10
    """
    if n < 0:
        raise ValueError("n muss nicht-negativ sein.")
    return n * n


def pentagonal_number(n: int) -> int:
    """
    Berechnet die n-te Pentagonalzahl P_n = n(3n-1)/2.

    Darstellung: Punkte in Fünfeckform.
    P_1=1, P_2=5, P_3=12, P_4=22, ...

    @param n Positive ganze Zahl
    @return  n-te Pentagonalzahl

    @timestamp 2026-03-10
    """
    if n < 0:
        raise ValueError("n muss nicht-negativ sein.")
    return n * (3 * n - 1) // 2


def hexagonal_number(n: int) -> int:
    """
    Berechnet die n-te Hexagonalzahl H_n = n(2n-1).

    Darstellung: Punkte in Sechseckform.
    H_1=1, H_2=6, H_3=15, H_4=28, ...

    @param n Positive ganze Zahl
    @return  n-te Hexagonalzahl

    @timestamp 2026-03-10
    """
    if n < 0:
        raise ValueError("n muss nicht-negativ sein.")
    return n * (2 * n - 1)


def is_triangular(n: int) -> bool:
    """
    Prüft, ob n eine Dreieckszahl ist.

    n = T_k ⟺ 8n + 1 ist perfektes Quadrat.

    @param n Nicht-negative ganze Zahl
    @return  True wenn n eine Dreieckszahl ist

    @timestamp 2026-03-10
    """
    if n < 0:
        return False
    # 8n + 1 muss perfektes Quadrat sein
    discriminant = 8 * n + 1
    sqrt_d = int(math.isqrt(discriminant))
    return sqrt_d * sqrt_d == discriminant


def polygonal_number(s: int, n: int) -> int:
    """
    Berechnet die n-te s-eckige Zahl (polygonale Zahl).

    Formel: P(s, n) = n · ((s-2)·n - (s-4)) / 2

    Spezialfälle:
    - s=3: Dreieckszahlen
    - s=4: Quadratzahlen
    - s=5: Pentagonalzahlen
    - s=6: Hexagonalzahlen

    @param s Anzahl der Ecken (s ≥ 3)
    @param n Positiver Index
    @return  n-te s-eckige Zahl

    @timestamp 2026-03-10
    """
    if s < 3:
        raise ValueError("s muss mindestens 3 sein (Dreieck).")
    if n < 0:
        raise ValueError("n muss nicht-negativ sein.")
    return n * ((s - 2) * n - (s - 4)) // 2


def cannonball_problem(n: int) -> dict:
    """
    Kanonenkugelproblem: Für welche n ist die Summe 1² + 2² + … + n²
    eine perfekte Quadratzahl?

    Summe der Quadrate: S_n = n(n+1)(2n+1)/6
    Bekannte Lösung: n=24 → S_24 = 4900 = 70² (einzige Lösung > 1).

    @param n Maximale Anzahl von Lagen
    @return  Wörterbuch mit gefundenen Lösungen bis n

    @timestamp 2026-03-10
    """
    solutions = []
    for k in range(1, n + 1):
        # Summe der Quadrate 1² + 2² + ... + k²
        s = k * (k + 1) * (2 * k + 1) // 6
        sqrt_s = int(math.isqrt(s))
        if sqrt_s * sqrt_s == s:
            solutions.append({
                'n': k,
                'summe': s,
                'quadratwurzel': sqrt_s,
                'gleichung': f"1² + 2² + ... + {k}² = {s} = {sqrt_s}²"
            })

    return {
        'gesuchte_n': n,
        'loesungen': solutions,
        'mathematischer_satz': (
            "Die einzigen n mit Σᵢ₌₁ⁿ i² = k² sind n=1 (S=1=1²) und n=24 (S=4900=70²)."
        )
    }
