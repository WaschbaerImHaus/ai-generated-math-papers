"""
@file recursion_theory.py
@brief Rekursionstheorie / Berechenbarkeitstheorie: Turing-Maschinen, μ-rekursive Funktionen,
       Halteproblem, Reduktionen, arithmetische Hierarchie, Komplexitätstheorie und
       Kolmogorov-Komplexität.

@description
    Dieses Modul implementiert die wichtigsten Konzepte der Rekursionstheorie:

    - Deterministische Turing-Maschinen (DTM) mit Simulator
    - Beispiel-TMs: Palindrome, Binärinkrement, Unäre Addition
    - μ-rekursive Funktionen: Nullfunktion, Nachfolger, Projektion,
      primitive Rekursion, μ-Operator
    - Ackermann-Funktion (total rekursiv, nicht primitiv rekursiv)
    - Halteproblem-Unentscheidbarkeit (Cantor-Diagonalargument)
    - Satz von Rice
    - Arithmetische Hierarchie (Σₙ/Πₙ/Δₙ)
    - Post'sches Korrespondenzproblem
    - Komplexitätsklassen P, NP, PSPACE
    - Kolmogorov-Komplexität (obere Schranke via zlib)

    Alle Beweise sind als strukturierte Dicts zurückgegeben, damit sie
    in anderen Modulen weiterverarbeitet oder angezeigt werden können.

    Formeln im KaTeX-Format sind in der Dokumentationsdatei recursion_theory.md
    zu finden.

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

from __future__ import annotations

import zlib
import sys
from enum import Enum
from typing import Callable, Optional, Any


# ==============================================================================
# Turing-Maschinen
# ==============================================================================

class TapeSymbol(Enum):
    """
    Spezialsymbole für das TM-Band.
    BLANK steht für eine leere Bandzelle (entspricht dem Symbol '_').
    """
    BLANK = '_'


class Direction(Enum):
    """
    Bewegungsrichtungen des Lese-/Schreibkopfes einer Turing-Maschine.
    LEFT: einen Schritt nach links, RIGHT: einen Schritt nach rechts,
    STAY: Kopf bleibt an aktueller Position (Erweiterung).
    """
    LEFT  = 'L'
    RIGHT = 'R'
    STAY  = 'S'


class TuringMachine:
    """
    Deterministische Turing-Maschine M = (Q, Σ, Γ, δ, q₀, q_acc, q_rej).

    @description
        - Q          : endliche Zustandsmenge
        - Σ          : Eingabealphabet (Teilmenge von Γ, ohne BLANK)
        - Γ          : Bandalphabet (enthält BLANK = '_')
        - δ          : (state, symbol) → (new_state, write_symbol, direction)
        - q₀         : Startzustand
        - q_acc      : Akzeptierzustand
        - q_rej      : Verwerfzustand

    Das Band wird intern als Dictionary {position: symbol} repräsentiert,
    sodass es in beide Richtungen unbegrenzt ist.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        states: set[str],
        input_alphabet: set[str],
        tape_alphabet: set[str],
        transitions: dict[tuple, tuple],
        initial_state: str,
        accept_state: str,
        reject_state: str
    ) -> None:
        """
        Initialisiert eine deterministische Turing-Maschine.

        @param states          Menge aller Zustände, z.B. {'q0', 'q1', 'accept', 'reject'}
        @param input_alphabet  Eingabealphabet, z.B. {'a', 'b'}
        @param tape_alphabet   Bandalphabet, z.B. {'a', 'b', '_'}
        @param transitions     Übergangsfunktion als Dict:
                               {(zustand, symbol): (neuer_zustand, schreibsymbol, richtung)}
                               Richtung ist ein Direction-Enum-Wert.
        @param initial_state   Startzustand
        @param accept_state    Akzeptierzustand (TM hält und akzeptiert)
        @param reject_state    Verwerfzustand (TM hält und verwirft)
        """
        self.states = states
        self.input_alphabet = input_alphabet
        self.tape_alphabet = tape_alphabet
        self.transitions = transitions
        self.initial_state = initial_state
        self.accept_state = accept_state
        self.reject_state = reject_state

    def _init_tape(self, input_str: str) -> dict[int, str]:
        """
        Initialisiert das Band mit der Eingabe.
        Position 0 ist das erste Zeichen der Eingabe.

        @param input_str  Die Eingabezeichenkette
        @return           Dictionary {position: symbol}
        """
        tape: dict[int, str] = {}
        for i, ch in enumerate(input_str):
            tape[i] = ch
        return tape

    def _read(self, tape: dict[int, str], head: int) -> str:
        """
        Liest das Symbol an Position head.
        Gibt BLANK zurück, falls die Zelle nicht beschrieben wurde.

        @param tape  Das Band als Dict
        @param head  Aktuelle Kopfposition
        @return      Gelesenes Symbol
        """
        return tape.get(head, TapeSymbol.BLANK.value)

    def _tape_to_str(self, tape: dict[int, str], head: int) -> str:
        """
        Wandelt das Band in einen lesbaren String um.
        Zeigt den relevanten Bereich um den Kopf herum.

        @param tape  Das Band als Dict
        @param head  Aktuelle Kopfposition
        @return      String-Darstellung des Bandes
        """
        if not tape:
            return TapeSymbol.BLANK.value
        lo = min(min(tape.keys()), head) - 1
        hi = max(max(tape.keys()), head) + 1
        return ''.join(tape.get(i, TapeSymbol.BLANK.value) for i in range(lo, hi + 1))

    def run(self, input_str: str, max_steps: int = 10000) -> dict:
        """
        Führt die Turing-Maschine auf der Eingabe aus.

        @description
            Die TM hält, wenn sie den Akzeptier- oder Verwerfzustand erreicht
            oder das Schritt-Limit überschritten wird (Timeout).

        @param input_str  Eingabe-String
        @param max_steps  Maximale Anzahl von Berechnungsschritten
        @return           Dict mit Schlüsseln:
                          - 'accepted': True/False
                          - 'steps':    Anzahl ausgeführter Schritte
                          - 'tape':     Bandzustand als String
                          - 'halted':   True wenn TM normal gehalten hat
        """
        tape  = self._init_tape(input_str)
        head  = 0
        state = self.initial_state
        steps = 0

        while steps < max_steps:
            # Haltebedingung prüfen
            if state == self.accept_state:
                return {
                    'accepted': True,
                    'steps': steps,
                    'tape': self._tape_to_str(tape, head),
                    'halted': True
                }
            if state == self.reject_state:
                return {
                    'accepted': False,
                    'steps': steps,
                    'tape': self._tape_to_str(tape, head),
                    'halted': True
                }

            symbol = self._read(tape, head)
            key = (state, symbol)

            # Kein Übergang definiert → implizit verwerfen
            if key not in self.transitions:
                return {
                    'accepted': False,
                    'steps': steps,
                    'tape': self._tape_to_str(tape, head),
                    'halted': True
                }

            new_state, write_sym, direction = self.transitions[key]
            tape[head] = write_sym

            # Kopf bewegen
            if direction == Direction.LEFT:
                head -= 1
            elif direction == Direction.RIGHT:
                head += 1
            # Direction.STAY: Kopf bleibt

            state  = new_state
            steps += 1

        # Timeout erreicht
        return {
            'accepted': False,
            'steps': steps,
            'tape': self._tape_to_str(tape, head),
            'halted': False
        }

    def accepts(self, input_str: str, max_steps: int = 10000) -> Optional[bool]:
        """
        Prüft ob die TM die Eingabe akzeptiert.

        @param input_str  Eingabe-String
        @param max_steps  Maximale Schritte
        @return           True  → akzeptiert
                          False → verworfen
                          None  → Timeout (TM hat nicht gehalten)
        """
        result = self.run(input_str, max_steps)
        if not result['halted']:
            return None
        return result['accepted']

    def configuration_sequence(self, input_str: str, max_steps: int = 20) -> list[str]:
        """
        Gibt die Sequenz von Konfigurationen der TM zurück.

        @description
            Eine Konfiguration ist ein Tripel (Bandzustand, Kopfposition, Zustand).
            Darstellung: 'BANDSYMBOL...⟨ZUSTAND⟩BANDSYMBOL...'

        @param input_str  Eingabe-String
        @param max_steps  Maximale Anzahl Schritte
        @return           Liste von Konfigurationsbeschreibungen als Strings
        """
        tape   = self._init_tape(input_str)
        head   = 0
        state  = self.initial_state
        configs: list[str] = []

        def make_config() -> str:
            """Erzeugt einen lesbaren Konfigurations-String."""
            if not tape:
                symbols = {0: TapeSymbol.BLANK.value}
            else:
                lo = min(min(tape.keys()), head)
                hi = max(max(tape.keys()), head)
                symbols = {i: tape.get(i, TapeSymbol.BLANK.value) for i in range(lo - 1, hi + 2)}
            parts = []
            for pos in sorted(symbols.keys()):
                if pos == head:
                    parts.append(f'[{state}]{symbols[pos]}')
                else:
                    parts.append(symbols[pos])
            return ''.join(parts)

        configs.append(make_config())

        for _ in range(max_steps):
            if state in (self.accept_state, self.reject_state):
                break

            symbol = self._read(tape, head)
            key = (state, symbol)
            if key not in self.transitions:
                break

            new_state, write_sym, direction = self.transitions[key]
            tape[head] = write_sym

            if direction == Direction.LEFT:
                head -= 1
            elif direction == Direction.RIGHT:
                head += 1

            state = new_state
            configs.append(make_config())

        return configs


# ==============================================================================
# Beispiel-Turing-Maschinen
# ==============================================================================

def tm_recognizes_palindromes() -> TuringMachine:
    """
    Erzeugt eine TM, die Palindrome über dem Alphabet {a, b} erkennt.

    @description
        Algorithmus: Merke erstes Zeichen, lösche es, gehe ans Ende,
        vergleiche mit letztem Zeichen, lösche es, gehe zurück.
        Wiederhole bis Band leer.

        Zustände:
        - q0 : Lese und merke erstes Zeichen
        - q1a : Gehe nach rechts, suche das Ende (letztes Zeichen ist 'a')
        - q1b : Gehe nach rechts, suche das Ende (letztes Zeichen ist 'b')
        - q2a : Vergleiche letztes Zeichen mit gemerkt 'a'
        - q2b : Vergleiche letztes Zeichen mit gemerkt 'b'
        - q3  : Gehe zurück zur Mitte
        - accept, reject

    @return TuringMachine-Objekt
    """
    BLANK = TapeSymbol.BLANK.value
    L = Direction.LEFT
    R = Direction.RIGHT

    states = {'q0', 'q1a', 'q1b', 'q2a', 'q2b', 'q3', 'accept', 'reject'}
    input_alpha = {'a', 'b'}
    tape_alpha  = {'a', 'b', 'X', BLANK}

    trans: dict[tuple, tuple] = {
        # q0: Lese erstes Zeichen
        ('q0', 'a') : ('q1a', 'X', R),   # Merke 'a', schreibe X (markiert gelesen)
        ('q0', 'b') : ('q1b', 'X', R),   # Merke 'b', schreibe X
        ('q0', 'X') : ('q0',  'X', R),   # Überspringe bereits markierte Felder
        ('q0', BLANK): ('accept', BLANK, R),  # Leeres Band → Palindrom

        # q1a: Bewege Kopf nach rechts bis zum Ende (merke 'a')
        ('q1a', 'a') : ('q1a', 'a', R),
        ('q1a', 'b') : ('q1a', 'b', R),
        ('q1a', 'X') : ('q1a', 'X', R),
        ('q1a', BLANK): ('q2a', BLANK, L),  # Ende erreicht, gehe einen Schritt zurück

        # q1b: Bewege Kopf nach rechts bis zum Ende (merke 'b')
        ('q1b', 'a') : ('q1b', 'a', R),
        ('q1b', 'b') : ('q1b', 'b', R),
        ('q1b', 'X') : ('q1b', 'X', R),
        ('q1b', BLANK): ('q2b', BLANK, L),

        # q2a: Vergleiche letztes Zeichen mit 'a'
        ('q2a', 'a') : ('q3', 'X', L),    # Stimmt überein, markiere und gehe zurück
        ('q2a', 'b') : ('reject', 'b', L),
        ('q2a', 'X') : ('accept', 'X', L),  # Nur X übrig → Palindrom (ungerade Länge)
        ('q2a', BLANK): ('accept', BLANK, L),

        # q2b: Vergleiche letztes Zeichen mit 'b'
        ('q2b', 'b') : ('q3', 'X', L),
        ('q2b', 'a') : ('reject', 'a', L),
        ('q2b', 'X') : ('accept', 'X', L),
        ('q2b', BLANK): ('accept', BLANK, L),

        # q3: Gehe nach links zurück zum nächsten unmarkierten Zeichen
        ('q3', 'a') : ('q3', 'a', L),
        ('q3', 'b') : ('q3', 'b', L),
        ('q3', 'X') : ('q0', 'X', R),     # Linkestes X gefunden, gehe eine Position vor
        ('q3', BLANK): ('q0', BLANK, R),
    }

    return TuringMachine(
        states=states,
        input_alphabet=input_alpha,
        tape_alphabet=tape_alpha,
        transitions=trans,
        initial_state='q0',
        accept_state='accept',
        reject_state='reject'
    )


def tm_binary_increment() -> TuringMachine:
    """
    Erzeugt eine TM, die eine Binärzahl um 1 erhöht.

    @description
        Algorithmus: Gehe ans rechte Ende, addiere 1 von rechts nach links.
        Bei Übertrag setze Bit auf 0, sonst setze Bit auf 1 und halt.
        Beispiel: 011 → 100

        Zustände:
        - q0  : Gehe ans Ende der Binärzahl
        - q1  : Addiere 1 (mit Übertrag) von rechts nach links
        - accept, reject

    @return TuringMachine-Objekt
    """
    BLANK = TapeSymbol.BLANK.value
    L = Direction.LEFT
    R = Direction.RIGHT

    states = {'q0', 'q1', 'accept', 'reject'}
    input_alpha = {'0', '1'}
    tape_alpha  = {'0', '1', BLANK}

    trans: dict[tuple, tuple] = {
        # q0: Gehe nach rechts bis zum Ende
        ('q0', '0')  : ('q0', '0', R),
        ('q0', '1')  : ('q0', '1', R),
        ('q0', BLANK): ('q1', BLANK, L),   # Ende erreicht, beginne Inkrement

        # q1: Addiere 1 mit Übertrag
        ('q1', '1')  : ('q1', '0', L),    # 1+1=10: schreibe 0, Übertrag nach links
        ('q1', '0')  : ('accept', '1', R), # 0+1=1: schreibe 1, fertig
        ('q1', BLANK): ('accept', '1', R), # Übertrag über MSB: neues Bit schreiben
    }

    return TuringMachine(
        states=states,
        input_alphabet=input_alpha,
        tape_alphabet=tape_alpha,
        transitions=trans,
        initial_state='q0',
        accept_state='accept',
        reject_state='reject'
    )


def tm_unary_addition() -> TuringMachine:
    """
    Erzeugt eine TM, die unäre Addition berechnet: 1^n 0 1^m → 1^(n+m).

    @description
        Eingabeformat: n Einsen, dann eine 0, dann m Einsen.
        Ausgabe: (n+m) Einsen ohne Trenner.

        Algorithmus:
        1. Finde die 0 (Trenner)
        2. Ersetze 0 durch 1
        3. Lösche letzte 1

        Zustände:
        - q0  : Gehe nach rechts bis zur 0
        - q1  : Ersetze 0 durch 1
        - q2  : Gehe nach rechts bis zum Ende
        - q3  : Lösche die letzte 1
        - accept, reject

    @return TuringMachine-Objekt
    """
    BLANK = TapeSymbol.BLANK.value
    L = Direction.LEFT
    R = Direction.RIGHT

    states = {'q0', 'q1', 'q2', 'q3', 'accept', 'reject'}
    input_alpha = {'0', '1'}
    tape_alpha  = {'0', '1', BLANK}

    trans: dict[tuple, tuple] = {
        # q0: Suche den Trenner 0
        ('q0', '1')  : ('q0', '1', R),
        ('q0', '0')  : ('q1', '1', R),    # Trenner gefunden, ersetze durch 1

        # q1: Gehe ans Ende des zweiten Blocks
        ('q1', '1')  : ('q1', '1', R),
        ('q1', BLANK): ('q2', BLANK, L),  # Ende erreicht

        # q2: Lösche die letzte 1 (der Trenner wurde bereits als 1 gezählt)
        ('q2', '1')  : ('accept', BLANK, L),
    }

    return TuringMachine(
        states=states,
        input_alphabet=input_alpha,
        tape_alphabet=tape_alpha,
        transitions=trans,
        initial_state='q0',
        accept_state='accept',
        reject_state='reject'
    )


# ==============================================================================
# μ-rekursive Funktionen
# ==============================================================================

def zero_function(n: int) -> int:
    """
    Nullfunktion Z: ℕ → ℕ, Z(n) = 0.

    @description
        Eine der drei Grundfunktionen der μ-rekursiven Funktionen.
        Gibt für jede natürliche Zahl 0 zurück.

    @param n  Eingabe (beliebige natürliche Zahl)
    @return   Immer 0
    """
    return 0


def successor_function(n: int) -> int:
    """
    Nachfolgerfunktion S: ℕ → ℕ, S(n) = n + 1.

    @description
        Eine der drei Grundfunktionen der μ-rekursiven Funktionen.
        Entspricht dem Nachfolger in der Peano-Arithmetik.

    @param n  Eingabe
    @return   n + 1
    """
    return n + 1


def projection(args: tuple, i: int) -> int:
    """
    Projektionsfunktion Uᵢⁿ: ℕⁿ → ℕ, Uᵢⁿ(x₁,...,xₙ) = xᵢ.

    @description
        Eine der drei Grundfunktionen der μ-rekursiven Funktionen.
        Gibt das i-te Argument zurück (0-indiziert).

    @param args  Tupel der Eingaben (x₁, ..., xₙ)
    @param i     Index des gewünschten Arguments (0-basiert)
    @return      args[i]
    """
    return args[i]


def primitive_recursion(base: Callable, step: Callable) -> Callable:
    """
    Primitive Rekursion: Schema zum Aufbau neuer primitiv-rekursiver Funktionen.

    @description
        Gegeben:
          g(x̄)        → Basisfall (n-stellig)
          h(x̄, y, z)  → Schrittfunktion (n+2-stellig), z = bisheriges Ergebnis

        Definiert:
          f(x̄, 0)   = g(x̄)
          f(x̄, y+1) = h(x̄, y, f(x̄, y))

        Alle primitiv-rekursiven Funktionen entstehen durch endliche Anwendung
        von Komposition und primitivem Rekursionsschema aus {Z, S, Uᵢⁿ}.

    @param base  Funktion g, aufgerufen als g(*x_bar)
    @param step  Funktion h, aufgerufen als h(*x_bar, y, prev)
    @return      Neue Funktion f(*x_bar, y) → int
    """
    def f(*args: int) -> int:
        """
        Ausführung des primitiven Rekursionsschemas.
        args[-1] ist das Rekursionsargument y, args[:-1] sind die freien Variablen x̄.
        """
        *x_bar, y = args
        # Basisfall: f(x̄, 0) = g(x̄)
        result = base(*x_bar)
        # Schritt: f(x̄, k+1) = h(x̄, k, f(x̄, k))
        for k in range(y):
            result = step(*x_bar, k, result)
        return result

    return f


def mu_operator(predicate: Callable) -> Callable:
    """
    μ-Operator (Minimierungsoperator): μy[P(x̄, y) = 0].

    @description
        Gegeben: P(x̄, y) → {0, 1}  (oder int, wobei 0 = "wahr")

        μy[P(x̄, y)] = kleinste y ∈ ℕ mit P(x̄, y) = 0,
                       falls ein solches y existiert.

        Falls kein y existiert, divergiert die Funktion (PartialRecursion).
        Der μ-Operator macht primitive Rekursion zur totalen/partiellen
        μ-Rekursion (= Turing-Berechenbarkeit).

    @param predicate  Funktion P, aufgerufen als P(*x_bar, y) → int
    @return           Partielle Funktion mu_f(*x_bar, max_y=10000) → int
    """
    def mu_f(*args: int, max_y: int = 10_000) -> int:
        """
        Sucht das kleinste y mit predicate(*args, y) == 0.

        @param args   Freie Variablen x̄
        @param max_y  Suchgrenze (verhindert Endlosschleife)
        @return       Kleinstes y mit P(x̄, y) = 0
        @raises ValueError  Falls kein y innerhalb von max_y gefunden
        """
        for y in range(max_y):
            if predicate(*args, y) == 0:
                return y
        raise ValueError(
            f"μ-Operator: Kein y in [0, {max_y}) gefunden mit predicate(*{args}, y) = 0"
        )

    return mu_f


def ackermann_function(m: int, n: int) -> int:
    """
    Ackermann-Funktion A(m, n): total rekursiv, aber nicht primitiv rekursiv.

    @description
        Definition:
          A(0, n)   = n + 1
          A(m+1, 0) = A(m, 1)
          A(m+1, n+1) = A(m, A(m+1, n))

        Die Ackermann-Funktion wächst schneller als jede primitiv-rekursive Funktion.
        Für große m, n ist die Berechnung praktisch unmöglich:
          A(4, 2) = 2^(2^(2^(2^2))) - 3  (eine Zahl mit ~19.728 Stellen)

        Wichtige Werte:
          A(0, n) = n+1
          A(1, n) = n+2
          A(2, n) = 2n+3
          A(3, n) = 2^(n+3) - 3
          A(4, n) = Turm aus (n+3) Zweien - 3

        Implementiert mit iterativem Stack-Aufruf (verhindert Python RecursionError).

    @param m  Erste Eingabe (≥ 0)
    @param n  Zweite Eingabe (≥ 0)
    @return   A(m, n)
    @raises RecursionError  Für sehr große m (≥ 5) praktisch nicht berechenbar
    """
    # Iterative Implementierung mit explizitem Stack
    stack = [m]
    n_curr = n

    while stack:
        m_top = stack[-1]

        if m_top == 0:
            # A(0, n) = n + 1
            stack.pop()
            n_curr += 1
        elif n_curr == 0:
            # A(m+1, 0) = A(m, 1)
            stack[-1] -= 1
            n_curr = 1
        else:
            # A(m+1, n+1) = A(m, A(m+1, n))
            # → Berechne erst A(m+1, n), dann A(m, Ergebnis)
            stack[-1] -= 1       # Für äußeren Aufruf A(m, ...)
            stack.append(m_top)  # Für inneren Aufruf A(m+1, n)
            n_curr -= 1

    return n_curr


def ackermann_growth_demo() -> dict:
    """
    Demonstriert das explosive Wachstum der Ackermann-Funktion.

    @description
        Zeigt ausgewählte Werte und erklärt warum die Funktion
        nicht primitiv rekursiv ist (sie wächst schneller als
        jede primitiv-rekursive Funktion).

    @return  Dict mit berechneten Werten und Erklärungen
    """
    return {
        'A(0,0)': ackermann_function(0, 0),   # = 1
        'A(1,1)': ackermann_function(1, 1),   # = 3
        'A(2,2)': ackermann_function(2, 2),   # = 7
        'A(3,3)': ackermann_function(3, 3),   # = 61
        'A(3,4)': ackermann_function(3, 4),   # = 125
        'A(4,1)': 65533,                       # Zu groß für direkte Berechnung
        'A(4,2)_digits': 19728,                # Anzahl Dezimalstellen von A(4,2)
        'pattern_m0': 'A(0, n) = n + 1  (lineare +1)',
        'pattern_m1': 'A(1, n) = n + 2  (lineare +2)',
        'pattern_m2': 'A(2, n) = 2n + 3 (linear)',
        'pattern_m3': 'A(3, n) = 2^(n+3) - 3 (exponentiell)',
        'pattern_m4': 'A(4, n) = Turm aus (n+3) Zweien - 3 (tetrational)',
        'not_primitive_recursive': (
            'Für jede primitiv-rekursive Funktion f existiert ein m₀ mit '
            'A(m₀, n) > f(n) für alle hinreichend großen n. '
            'Daher kann A nicht primitiv-rekursiv sein.'
        )
    }


# ==============================================================================
# Halteproblem und Entscheidbarkeit
# ==============================================================================

def halting_problem_undecidability_proof() -> dict:
    """
    Cantor-Diagonalisierungs-Beweis der Unentscheidbarkeit des Halteproblems.

    @description
        HALT = {⟨M, w⟩ : M hält auf Eingabe w}

        Beweis durch Widerspruch (Diagonalargument):
        1. Annahme: HALT ist entscheidbar, also existiert TM H mit:
              H(⟨M, w⟩) = accept, falls M auf w hält
              H(⟨M, w⟩) = reject, sonst
        2. Konstruiere D(⟨M⟩):
              Wenn H(⟨M, M⟩) = accept (M hält auf sich selbst):
                  D schleift (hält NICHT)
              Wenn H(⟨M, M⟩) = reject (M hält NICHT auf sich selbst):
                  D hält (akzeptiert)
        3. Was passiert bei D(⟨D⟩)?
              Falls D auf ⟨D⟩ hält → H(⟨D,D⟩) = accept → D schleift  → Widerspruch
              Falls D auf ⟨D⟩ nicht hält → H(⟨D,D⟩) = reject → D hält → Widerspruch
        4. Beide Fälle ergeben Widerspruch → HALT ist unentscheidbar. □

    @return  Dict mit Beweis-Metadaten
    """
    return {
        'theorem': 'Das Halteproblem HALT ist unentscheidbar.',
        'proof_method': 'Cantor-Diagonalargument (Widerspruchsbeweis)',
        'problem': 'HALT = {⟨M, w⟩ : Turing-Maschine M hält auf Eingabe w}',
        'assumption': 'Annahme: HALT ist entscheidbar durch TM H.',
        'diagonal_machine': (
            'Definiere D(⟨M⟩): '
            'Falls H(⟨M, M⟩) = accept → schleife; '
            'Falls H(⟨M, M⟩) = reject → halte.'
        ),
        'contradiction': (
            'D(⟨D⟩): '
            'Falls D hält → H akzeptiert → D schleift (Widerspruch). '
            'Falls D schleift → H verwirft → D hält (Widerspruch).'
        ),
        'conclusion': 'HALT ist unentscheidbar. Kein Algorithmus kann allgemein bestimmen, ob eine TM hält.',
        'corollary': (
            'Folgerung: Es gibt unendlich viele unentscheidbare Probleme '
            '(via Many-One-Reduktionen von HALT).'
        ),
        'historical_note': (
            'Alan Turing bewies dies 1936 in "On Computable Numbers, '
            'with an Application to the Entscheidungsproblem".'
        )
    }


def rice_theorem_demo() -> dict:
    """
    Satz von Rice: Jede nicht-triviale semantische Eigenschaft von TMs ist unentscheidbar.

    @description
        Formale Aussage:
        Sei C eine Klasse von partiell berechenbaren Funktionen.
        Falls {} ⊊ C ⊊ {alle partiell berechenbaren Funktionen},
        dann ist die Sprache L_C = {⟨M⟩ : φ_M ∈ C} unentscheidbar.

        Nicht-trivial bedeutet:
        - Mindestens eine TM berechnet eine Funktion in C
        - Mindestens eine TM berechnet eine Funktion außerhalb von C

        Beispiele nicht-trivialer (also unentscheidbarer) Eigenschaften:
        - "M hält auf der leeren Eingabe"
        - "M berechnet eine totale Funktion"
        - "M akzeptiert mindestens eine Eingabe"
        - "M und M' berechnen die gleiche Funktion"

        Beispiele trivialer (also möglicherweise entscheidbarer) Eigenschaften:
        - "M hat genau 5 Zustände" (syntaktisch, nicht semantisch)

    @return  Dict mit Theorem-Information und Beispielen
    """
    return {
        'theorem': 'Satz von Rice (1951)',
        'formal_statement': (
            'Sei C eine nicht-triviale Klasse partiell berechenbarer Funktionen. '
            'Dann ist L_C = {⟨M⟩ : φ_M ∈ C} unentscheidbar.'
        ),
        'proof_idea': (
            'Reduktion von ATM (Akzeptanzproblem) auf L_C: '
            'Konstruiere aus ⟨M, w⟩ eine TM M\' mit: '
            'M\'(x) = simuliere M auf w; falls akzeptiert, simuliere f(x). '
            'Dann gilt φ_{M\'} ∈ C ⟺ M akzeptiert w.'
        ),
        'non_trivial_condition': (
            'C ist nicht-trivial: ∅ ⊊ C ⊊ Menge aller partiell berechenbaren Funktionen.'
        ),
        'undecidable_examples': [
            'Hält M auf leerer Eingabe?',
            'Ist die von M berechnete Funktion total?',
            'Akzeptiert M mindestens eine Eingabe?',
            'Berechnen M₁ und M₂ dieselbe Funktion?',
            'Ist L(M) regulär / kontextfrei / rekursiv aufzählbar?'
        ],
        'decidable_examples': [
            'Hat M genau k Zustände? (syntaktische, keine semantische Eigenschaft)',
        ],
        'consequence': (
            'Es ist prinzipiell unmöglich, einen allgemeinen Programmanalysator '
            'zu schreiben, der beliebige nicht-triviale Eigenschaften '
            'von Programmen korrekt bestimmt.'
        )
    }


def decidable_problems() -> list[dict]:
    """
    Liste klassischer entscheidbarer Probleme mit Begründung.

    @description
        Ein Problem ist entscheidbar, falls eine TM existiert,
        die auf jeder Eingabe hält und korrekt akzeptiert oder verwirft.

    @return  Liste von Dicts mit Schlüsseln: 'name', 'description', 'reason', 'complexity'
    """
    return [
        {
            'name': 'Wortproblem für reguläre Sprachen',
            'description': 'Gegeben DFA A und Wort w: Gilt w ∈ L(A)?',
            'reason': (
                'Simuliere den DFA auf w: linearer Durchlauf, '
                'deterministisch, hält immer.'
            ),
            'complexity': 'O(|w|)'
        },
        {
            'name': 'Äquivalenz regulärer Ausdrücke',
            'description': 'Gegeben RegEx r₁, r₂: Gilt L(r₁) = L(r₂)?',
            'reason': (
                'Beide in DFAs übersetzen, Produkt-Automat für symmetrische '
                'Differenz konstruieren, Erreichbarkeit leerer Sprache prüfen.'
            ),
            'complexity': 'PSPACE-vollständig'
        },
        {
            'name': 'Halteproblem für beschränkte TMs',
            'description': 'Hält TM M auf w innerhalb von t Schritten?',
            'reason': (
                'Simuliere M für genau t Schritte. Falls nicht gehalten, '
                'verwerfe. Simulation terminiert immer.'
            ),
            'complexity': 'O(t · |M|)'
        },
        {
            'name': 'Presburger-Arithmetik',
            'description': (
                'Entscheidbarkeit von Aussagen über ℕ mit Addition '
                '(ohne Multiplikation): ∃∀-Formeln über (ℕ, +, 0, 1).'
            ),
            'reason': (
                'Quantoren-Elimination (Cooper-Algorithmus): '
                'Jede Presburger-Formel ist äquivalent zu einer quantorenfreien.'
            ),
            'complexity': '2^(2^n) (mehrfach exponentiell)'
        },
        {
            'name': 'Wortproblem für kontextfreie Grammatiken',
            'description': 'Gegeben CFG G und Wort w: Gilt w ∈ L(G)?',
            'reason': (
                'CYK-Algorithmus (Cocke-Younger-Kasami) in O(|w|³ · |G|).'
            ),
            'complexity': 'O(n³)'
        },
        {
            'name': 'Leerheitsproblem für kontextfreie Grammatiken',
            'description': 'Gegeben CFG G: Gilt L(G) = ∅?',
            'reason': (
                'Prüfe ob Startsymbol aus einem Terminalterminal ableitbar. '
                'Markierungsalgorithmus auf dem Grammatikgraphen.'
            ),
            'complexity': 'Polynomiell'
        }
    ]


def undecidable_problems() -> list[dict]:
    """
    Liste klassischer unentscheidbarer Probleme mit Begründungsangaben.

    @return  Liste von Dicts mit Schlüsseln: 'name', 'description', 'reduction_from', 'year'
    """
    return [
        {
            'name': 'Halteproblem (HALT / ATM)',
            'description': 'Hält TM M auf Eingabe w?',
            'reduction_from': 'Direkter Beweis (Diagonalargument)',
            'year': 1936,
            'notes': 'Σ₁-vollständig (rekursiv aufzählbar, aber nicht entscheidbar)'
        },
        {
            'name': 'Äquivalenz kontextfreier Grammatiken',
            'description': 'Gegeben G₁, G₂: Gilt L(G₁) = L(G₂)?',
            'reduction_from': 'Post-Korrespondenzproblem',
            'year': 1961,
            'notes': 'Π₂-vollständig'
        },
        {
            'name': 'Post-Korrespondenzproblem (PCP)',
            'description': (
                'Gegeben Dominos (tᵢ, bᵢ): Existiert Sequenz i₁,...,iₙ '
                'mit t_{i₁}...t_{iₙ} = b_{i₁}...b_{iₙ}?'
            ),
            'reduction_from': 'Halteproblem (modifizierte PCP-Variante)',
            'year': 1946,
            'notes': 'Wichtiges Zwischenproblem für weitere Unentscheidbarkeitsbeweise'
        },
        {
            'name': "Hilberts 10. Problem",
            'description': (
                'Gegeben diophantische Gleichung P(x₁,...,xₙ) = 0: '
                'Hat sie eine ganzzahlige Lösung?'
            ),
            'reduction_from': 'Halteproblem (MRDP-Theorem)',
            'year': 1970,
            'notes': (
                'Matiyasevich-Robinson-Davis-Putnam: Jede r.e. Menge ist '
                'diophantisch darstellbar.'
            )
        },
        {
            'name': 'Collatz-Problem (offen)',
            'description': (
                'Gilt für jede natürliche Zahl n > 0: '
                'Die Collatz-Folge (3n+1 oder n/2) erreicht schließlich 1?'
            ),
            'reduction_from': 'Offen - möglicherweise unentscheidbar',
            'year': 1937,
            'notes': 'Bisher kein Beweis bekannt; gilt für alle n < 2^68 (numerisch)'
        },
        {
            'name': 'Totale TM (TOTAL)',
            'description': 'Gegeben TM M: Hält M auf allen Eingaben?',
            'reduction_from': 'Halteproblem (Π₂-Vollständigkeit)',
            'year': None,
            'notes': 'Π₂-vollständig (nicht in Σ₁ ∪ Π₁)'
        }
    ]


# ==============================================================================
# Reduktionen
# ==============================================================================

def many_one_reduction(A: str, B: str) -> dict:
    """
    Beschreibt eine Many-One-Reduktion A ≤_m B mit konkreten Beispielen.

    @description
        A ≤_m B bedeutet: Es gibt eine totale berechenbare Funktion f: Σ* → Σ* mit
            x ∈ A  ⟺  f(x) ∈ B
        für alle x ∈ Σ*.

        Eigenschaften:
        - Falls B entscheidbar → A entscheidbar
        - Falls A unentscheidbar → B unentscheidbar
        - Many-One-Reduktion ist transitiv und reflexiv
        - Schwächer als Turing-Reduktion (kein Oracle, keine Mehrfachnutzung)

    @param A  Name des Ausgangs-Problems
    @param B  Name des Ziel-Problems
    @return   Dict mit Erklärung der Reduktion
    """
    # Klassische Reduktionen als Nachschlagewerk
    known_reductions = {
        ('HALT', 'ATM'): {
            'description': 'HALT ≤_m ATM: Aus ⟨M, w⟩ konstruiere M\' die w ignoriert und M simuliert. M\' akzeptiert ⟺ M hält auf w.',
            'function': 'f(⟨M, w⟩) = ⟨M\', w\'⟩ wo M\' = "simuliere M auf w, akzeptiere falls halt"'
        },
        ('ATM', 'HALT'): {
            'description': 'ATM ≤_m HALT: Aus ⟨M, w⟩ konstruiere M\' die bei Ablehnung schleift.',
            'function': 'f(⟨M, w⟩) = ⟨M\', w⟩ wo M\' = "simuliere M; bei Akzeptanz halte, bei Ablehnung schleife"'
        },
        ('PCP', 'CFG-EQ'): {
            'description': 'PCP ≤_m Äquivalenzproblem für CFGs: Codiere PCP-Instanz als zwei CFGs.',
            'function': 'Komplexe Konstruktion über PCP-Grammatiken'
        }
    }

    key = (A, B)
    if key in known_reductions:
        base = known_reductions[key]
        return {
            'reduction': f'{A} ≤_m {B}',
            'type': 'Many-One-Reduktion',
            'description': base['description'],
            'reduction_function': base['function'],
            'implication': f'Wenn {B} entscheidbar wäre, wäre auch {A} entscheidbar.'
        }

    return {
        'reduction': f'{A} ≤_m {B}',
        'type': 'Many-One-Reduktion',
        'description': (
            f'Eine berechenbare Funktion f existiert mit: '
            f'x ∈ {A} ⟺ f(x) ∈ {B}.'
        ),
        'implication': (
            f'Falls {B} entscheidbar → {A} entscheidbar. '
            f'Falls {A} unentscheidbar → {B} unentscheidbar.'
        ),
        'note': 'Keine spezifische Reduktion für dieses Paar implementiert.'
    }


def turing_reduction(A: str, B: str) -> dict:
    """
    Beschreibt eine Turing-Reduktion A ≤_T B.

    @description
        A ≤_T B bedeutet: Es gibt eine Oracle-TM M^B die A entscheidet,
        wobei M^B das Oracle für B beliebig oft und mit beliebigen Fragen nutzen darf.

        Stärker als Many-One-Reduktion:
        - Oracle darf mehrfach aufgerufen werden
        - Ergebnis darf negiert werden (A ≤_T B ≠> co-A ≤_m B allgemein)
        - co-HALT ≤_T HALT (da: Oracle für HALT, Antwort negieren)
        - Aber: co-HALT ≤_m HALT? Nein (HALT ∈ Σ₁, co-HALT ∈ Π₁, Π₁ ≠ Σ₁)

        Turing-Grade (Degrees of Unsolvability):
        - 0  = Grad der entscheidbaren Probleme
        - 0' = Grad des Halteproblems (HALT)
        - 0'' = Grad von HALT^HALT (Halteproblem relativiert zu 0')

    @param A  Name des Ausgangs-Problems
    @param B  Name des Ziel-Problems
    @return   Dict mit Erklärung
    """
    return {
        'reduction': f'{A} ≤_T {B}',
        'type': 'Turing-Reduktion (Oracle-Reduktion)',
        'description': (
            f'Es gibt eine Oracle-TM M^{{{B}}} die {A} entscheidet, '
            f'indem sie das Oracle für {B} als Blackbox nutzt.'
        ),
        'difference_to_many_one': (
            'Anders als bei ≤_m darf das Oracle '
            'mehrfach aufgerufen und die Antwort kann beliebig verarbeitet werden.'
        ),
        'example': (
            'co-HALT ≤_T HALT: Oracle-TM fragt HALT und negiert die Antwort. '
            'Aber co-HALT ≤_m HALT gilt NICHT.'
        ),
        'implication': (
            f'Falls {B} entscheidbar → {A} entscheidbar. '
            f'Die Umkehrung gilt nicht allgemein bei ≤_T.'
        )
    }


def reduction_degrees() -> dict:
    """
    Beschreibt die Struktur der Turing-Grade (Degrees of Unsolvability).

    @description
        Die Turing-Grade bilden eine partielle Ordnung (≤_T).
        Wichtige Grade:
        - deg(0): Alle entscheidbaren Probleme (Rekursiven)
        - deg(0'): Das Halteproblem HALT (Σ₁-vollständig)
        - deg(0''): HALT' = Halteproblem eines HALT-Oracle
        - deg(0⁽ⁿ⁾): n-faches Sprungoperator-Bild von 0

        Post-Problem (offen): Gibt es einen Grad d mit 0 < d < 0'?
        (Beantwortet: Ja, durch Friedberg-Muchnik 1956/57 - "rekursiv aufzählbare Grade")

    @return  Dict mit Graden und Beispielen
    """
    return {
        'degree_0': {
            'name': 'Grad 0 (Rekursive Sprachen)',
            'description': 'Entscheidbare Probleme',
            'examples': [
                'Wortproblem für reguläre Sprachen',
                'Wortproblem für CFGs',
                'Presburger-Arithmetik'
            ]
        },
        'degree_0_prime': {
            'name': "Grad 0' (Halteproblem)",
            'description': 'Σ₁-vollständige Probleme',
            'examples': [
                'HALT = {⟨M, w⟩ : M hält auf w}',
                'ATM = {⟨M, w⟩ : M akzeptiert w}',
                'Leerheitsproblem für TMs'
            ]
        },
        'degree_0_double_prime': {
            'name': "Grad 0'' (Halteproblem relativiert zu HALT)",
            'description': 'Σ₂-vollständige Probleme',
            'examples': [
                'TOTAL = {⟨M⟩ : M hält auf allen Eingaben}',
                'INF = {⟨M⟩ : |L(M)| = ∞}'
            ]
        },
        'post_problem': {
            'description': (
                'Post 1944: Gibt es Grade zwischen 0 und 0\'? '
                'Antwort: Ja (Friedberg-Muchnik 1956/57). '
                'Solche Grade heißen "rekursiv aufzählbar, aber nicht entscheidbar, '
                'und nicht Turing-äquivalent zu HALT".'
            )
        }
    }


# ==============================================================================
# Arithmetische Hierarchie
# ==============================================================================

def arithmetical_hierarchy() -> dict:
    """
    Die Arithmetische Hierarchie: Klassifikation von Teilmengen von ℕ.

    @description
        Die arithmetische Hierarchie klassifiziert Mengen anhand der
        logischen Komplexität ihrer Definitionen in Peano-Arithmetik.

        Quantoren-Alternationen:
        Σ₀ = Π₀ = Δ₁ : Entscheidbare Mengen (rekursiv, keine unbeschränkten Quantoren)
        Σ₁ : ∃-Formeln (rekursiv aufzählbar, r.e.)
        Π₁ : ∀-Formeln (co-r.e., Komplement von Σ₁-Mengen)
        Δ₂ : Σ₂ ∩ Π₂
        Σ₂ : ∃∀-Formeln
        Π₂ : ∀∃-Formeln
        ...
        Σₙ₊₁ : ∃y : φ(x, y) mit φ ∈ Πₙ
        Πₙ₊₁ : ∀y : φ(x, y) mit φ ∈ Σₙ

        Alle arithmetisch definierbaren Mengen: ⋃ Σₙ = arithmetische Hierarchie

    @return  Dict mit Level-Beschreibungen und Beispielen
    """
    return {
        'Σ₀': {
            'also_called': 'Δ₁ = Π₀',
            'definition': 'Entscheidbare (rekursive) Mengen. Beschränkte Quantoren (∃x < n, ∀x < n).',
            'examples': [
                'Menge der Primzahlen',
                'Menge der geraden Zahlen',
                'Wortproblem für reguläre Sprachen'
            ],
            'completeness': None
        },
        'Σ₁': {
            'definition': '∃y : R(x, y) mit R entscheidbar. Rekursiv aufzählbare (r.e.) Mengen.',
            'examples': [
                'HALT = {⟨M, w⟩ : M hält auf w}  (Halteproblem)',
                'ATM = {⟨M, w⟩ : M akzeptiert w}',
                'Menge der Gödelzahlen beweisbarer Sätze (PA)'
            ],
            'completeness': 'ATM ist Σ₁-vollständig.'
        },
        'Π₁': {
            'definition': '∀y : R(x, y) mit R entscheidbar. Komplement von Σ₁-Mengen.',
            'examples': [
                'co-HALT = {⟨M, w⟩ : M schleift auf w}',
                'TOTAL = {⟨M⟩ : M hält auf allen Eingaben}',
                'INF = {⟨M⟩ : L(M) ist unendlich}'
            ],
            'completeness': 'co-HALT ist Π₁-vollständig.'
        },
        'Σ₂': {
            'definition': '∃y₁ ∀y₂ : R(x, y₁, y₂) mit R entscheidbar.',
            'examples': [
                'FIN = {⟨M⟩ : |L(M)| < ∞}  (Endlichkeitsproblem)',
                'COF = {⟨M⟩ : L(M) ist cofinit}',
                'Menge der erfüllbaren Π₁-Formeln'
            ],
            'completeness': 'FIN ist Σ₂-vollständig.'
        },
        'Π₂': {
            'definition': '∀y₁ ∃y₂ : R(x, y₁, y₂) mit R entscheidbar.',
            'examples': [
                'TOTAL = {⟨M⟩ : M hält auf allen Eingaben}  (auch Π₂)',
                'EQ_TM = {⟨M₁, M₂⟩ : L(M₁) = L(M₂)}'
            ],
            'completeness': 'TOTAL ist Π₂-vollständig.'
        },
        'Σ₃': {
            'definition': '∃∀∃-Quantorenstruktur.',
            'examples': [
                'REC = {⟨M⟩ : L(M) ist entscheidbar}  (Σ₃-vollständig)'
            ],
            'completeness': 'REC ist Σ₃-vollständig.'
        },
        'structure': {
            'inclusions': 'Σₙ ∪ Πₙ ⊆ Δₙ₊₁ ⊆ Σₙ₊₁ ∩ Πₙ₊₁',
            'proper': 'Alle Inklusionen sind echt (echte Hierarchie)',
            'union': 'Arithmetische Hierarchie = ⋃ₙ Σₙ (abzählbar additive Mengen)'
        }
    }


def sigma1_complete_problem() -> dict:
    """
    ATM (Akzeptanzproblem) ist Σ₁-vollständig.

    @description
        ATM = {⟨M, w⟩ : M akzeptiert w}

        Beweis in zwei Teilen:
        1. ATM ∈ Σ₁: Nichtdeterministische Simulation zeigt Akzeptanz.
           ⟨M, w⟩ ∈ ATM ⟺ ∃t: M akzeptiert w in t Schritten.
           "∃t" ist ein existentieller Quantor über die Berechnungshistorie.

        2. ATM ist Σ₁-hart: Für jedes L ∈ Σ₁ gilt L ≤_m ATM.
           Sei L ∈ Σ₁, also L = {x : M_L hält und akzeptiert x}.
           Many-One-Reduktion: f(x) = ⟨M_L, x⟩.
           Dann: x ∈ L ⟺ M_L akzeptiert x ⟺ ⟨M_L, x⟩ ∈ ATM.

    @return  Dict mit Beweis-Teilen
    """
    return {
        'problem': 'ATM = {⟨M, w⟩ : Turing-Maschine M akzeptiert Eingabe w}',
        'claim': 'ATM ist Σ₁-vollständig.',
        'part1_membership': {
            'claim': 'ATM ∈ Σ₁',
            'proof': (
                '⟨M, w⟩ ∈ ATM ⟺ ∃t ∈ ℕ: "M akzeptiert w in ≤ t Schritten". '
                'Der Teilausdruck ist entscheidbar (endliche Simulation), '
                'also ist ATM durch eine ∃-Formel definiert → ATM ∈ Σ₁.'
            )
        },
        'part2_hardness': {
            'claim': 'ATM ist Σ₁-hart (jedes r.e. Problem ist ≤_m ATM)',
            'proof': (
                'Sei L ∈ Σ₁, d.h. L = L(M_L) für eine TM M_L. '
                'Definiere f(x) = ⟨M_L, x⟩. '
                'Dann: x ∈ L ⟺ M_L akzeptiert x ⟺ f(x) ∈ ATM. '
                'f ist offensichtlich berechenbar → L ≤_m ATM.'
            )
        },
        'conclusion': (
            'ATM ist Σ₁-vollständig: '
            'Es ist das "schwerste" r.e. Problem bezüglich ≤_m.'
        )
    }


# ==============================================================================
# Post'sches Korrespondenzproblem
# ==============================================================================

def post_correspondence_problem(dominoes: list[tuple]) -> Optional[bool]:
    """
    Post'sches Korrespondenzproblem (PCP): Entscheidet für kleine Instanzen.

    @description
        Gegeben: Liste von Dominos [(t₁, b₁), (t₂, b₂), ..., (tₙ, bₙ)].
        Frage: Gibt es eine nicht-leere Sequenz i₁, i₂, ..., iₖ (mit Wiederholung)
               so dass: t_{i₁} t_{i₂} ... t_{iₖ} = b_{i₁} b_{i₂} ... b_{iₖ}?

        Algorithmus: Backtracking-Suche (BFS/DFS) mit Zustandsspeicher.
        Zustand = (oberer String, unterer String).

        Für kleine Instanzen entscheidbar, allgemein unentscheidbar.
        Gibt None zurück, falls Suchtiefe überschritten wird.

    @param dominoes  Liste von Tupeln (top, bottom) mit Strings
    @return          True falls Lösung gefunden, False falls keine Lösung,
                     None falls Timeout (Instanz zu groß)
    """
    if not dominoes:
        return False

    # BFS: Zustand = (aktueller_top_string, aktueller_bottom_string, sequenz_laenge)
    from collections import deque

    max_sequence_length = 8   # Maximale Sequenzlänge für Backtracking
    max_string_length   = 50  # Maximale Stringlänge (verhindert Explosion)

    # Initialer Zustand: leere Strings
    # Wir nutzen Präfixvergleich: top = aktueller oberer String, bottom = aktueller unterer
    queue: deque = deque()
    visited: set = set()

    # Starte mit jedem Domino als erstem Element
    for i, (top, bot) in enumerate(dominoes):
        state = (top, bot, 1)
        queue.append(state)

    while queue:
        top, bot, depth = queue.popleft()

        # Beide gleich lang und identisch → Lösung gefunden
        if top == bot and top != '':
            return True

        # Kürzeren String verlängern (gemeinsamer Präfix)
        # Normalisierung: Subtrahiere gemeinsamen Präfix
        common = min(len(top), len(bot))
        # Prüfe ob einer ein Präfix des anderen ist
        if len(top) <= len(bot):
            if not bot.startswith(top):
                continue  # Kein gemeinsamer Präfix möglich
            remainder = ('', bot[len(top):])
        else:
            if not top.startswith(bot):
                continue
            remainder = (top[len(bot):], '')

        # Tiefenbeschränkung
        if depth >= max_sequence_length:
            # Nicht definitiv keine Lösung, aber zu tief
            # Wir geben None zurück wenn wir am Ende des Timeouts sind
            continue

        # Zu langer String
        if len(top) > max_string_length or len(bot) > max_string_length:
            continue

        # Zustandsnormalisierung gegen Zyklen
        state_key = (remainder[0], remainder[1])
        if state_key in visited:
            continue
        visited.add(state_key)

        # Alle Dominos anfügen
        for top_add, bot_add in dominoes:
            new_top = remainder[0] + top_add
            new_bot = remainder[1] + bot_add
            queue.append((new_top, new_bot, depth + 1))

    # Keine Lösung innerhalb der Grenzen gefunden
    # Bei sehr kleinen Instanzen ist das tatsächlich False
    if len(dominoes) <= 3 and max_sequence_length >= 8:
        return False
    return None


# ==============================================================================
# Komplexitätstheorie (Grundlagen)
# ==============================================================================

def time_complexity_classes() -> dict:
    """
    Übersicht der wichtigsten Zeitkomplexitätsklassen.

    @description
        Die Komplexitätstheorie klassifiziert Probleme nach dem Ressourcenverbrauch
        (Zeit, Platz) zur Lösung auf einer TM.

        Wichtige Klassen:
        - P      : Polynomial time (deterministische TM)
        - NP     : Nondeterministic Polynomial time
        - co-NP  : Komplement von NP
        - PSPACE : Polynomial space
        - EXPTIME: Exponentialzeit
        - EXPSPACE: Exponentialplatz

        Bekannte Inklusionen: P ⊆ NP ∩ co-NP ⊆ PSPACE ⊆ EXPTIME

        Offen: P = NP? (Millennium-Problem, Clay Institute)

    @return  Dict mit Klassenbeschreibungen und Beispielen
    """
    return {
        'P': {
            'definition': 'Von deterministischer TM in polynomieller Zeit entscheidbar.',
            'formal': 'P = ⋃_{k≥1} DTIME(nᵏ)',
            'examples': [
                'Sortieren (O(n log n))',
                'Kürzeste Wege (Dijkstra: O(n² log n))',
                'Lineare Programmierung (Ellipsoidmethode)',
                'Primalitätstest (AKS: O(n^6))'
            ]
        },
        'NP': {
            'definition': (
                'Von nichtdeterministischer TM in polynomieller Zeit entscheidbar. '
                'Äquivalent: Lösungen sind in polynomieller Zeit verifizierbar.'
            ),
            'formal': 'NP = ⋃_{k≥1} NTIME(nᵏ)',
            'examples': [
                'SAT (Erfüllbarkeit Boolescher Formeln)',
                'Traveling Salesman (Entscheidungsversion)',
                'Clique, Independent Set, Vertex Cover',
                'Rucksackproblem (Knapsack)'
            ]
        },
        'co_NP': {
            'definition': (
                'Komplement von NP. '
                'L ∈ co-NP ⟺ Lᶜ ∈ NP.'
            ),
            'examples': [
                'TAUTOLOGIE (ist Formel Tautologie?)',
                'UNSAT (ist Formel unerfüllbar?)',
                'Kein Hamilton-Kreis vorhanden'
            ]
        },
        'PSPACE': {
            'definition': 'Polynomial space (deterministische TM).',
            'formal': 'PSPACE = ⋃_{k≥1} DSPACE(nᵏ)',
            'examples': [
                'QBF (Quantified Boolean Formula)',
                'Spiellösbarkeit (z.B. Geography-Spiel)',
                'Äquivalenz regulärer Ausdrücke'
            ],
            'note': 'NP ⊆ PSPACE (Savitch: PSPACE = NPSPACE)'
        },
        'EXPTIME': {
            'definition': 'Exponentielle Zeit.',
            'formal': 'EXPTIME = ⋃_{k≥1} DTIME(2^{nᵖ})',
            'examples': [
                'Schachspiel (exakte Lösung)',
                'ALT (Alternating polynomial time) = PSPACE',
                'Presburger-Arithmetik (Untere Schranke)'
            ]
        },
        'inclusions': {
            'chain': 'P ⊆ NP ∩ co-NP ⊆ NP ⊆ PSPACE ⊆ EXPTIME',
            'open_problems': [
                'P =? NP  (Millennium-Problem)',
                'NP =? co-NP',
                'P =? PSPACE',
                'NP ⊊ EXPTIME (bewiesen: echte Inklusion)'
            ]
        }
    }


def np_complete_problems() -> list[dict]:
    """
    Klassische NP-vollständige Probleme.

    @description
        Ein Problem L ist NP-vollständig, falls:
        1. L ∈ NP
        2. Für alle L' ∈ NP gilt L' ≤_p L (Polynomial-Zeit-Reduktion)

        Falls ein NP-vollständiges Problem in P läge, wäre P = NP.

    @return  Liste von Dicts mit Schlüsseln: 'name', 'description', 'reduction_from'
    """
    return [
        {
            'name': 'SAT (Boolean Satisfiability)',
            'description': (
                'Gegeben Boolesche Formel φ: Gibt es eine Belegung '
                'der Variablen, die φ wahr macht?'
            ),
            'reduction_from': 'NP-Definition (Cook-Levin-Satz)',
            'year': 1971,
            'note': 'Erstes bewiesenes NP-vollständiges Problem'
        },
        {
            'name': '3-SAT',
            'description': (
                'SAT eingeschränkt auf Formeln in 3-KNF '
                '(Konjunktive Normalform mit je 3 Literalen pro Klausel).'
            ),
            'reduction_from': 'SAT',
            'year': 1972
        },
        {
            'name': 'CLIQUE',
            'description': (
                'Gegeben Graph G und k: Enthält G eine Clique (vollständigen Teilgraphen) '
                'der Größe k?'
            ),
            'reduction_from': '3-SAT',
            'year': 1972
        },
        {
            'name': 'VERTEX COVER',
            'description': (
                'Gegeben Graph G und k: Gibt es eine Knotenüberdeckung '
                '(Teilmenge der Knoten, die alle Kanten abdeckt) der Größe ≤ k?'
            ),
            'reduction_from': 'CLIQUE (Komplement)',
            'year': 1972
        },
        {
            'name': 'HAMILTONIAN CIRCUIT',
            'description': (
                'Gegeben Graph G: Enthält G einen Hamiltonkreis '
                '(Kreis der jeden Knoten genau einmal besucht)?'
            ),
            'reduction_from': '3-SAT',
            'year': 1972
        },
        {
            'name': 'KNAPSACK (0/1 Rucksack)',
            'description': (
                'Gegeben n Gegenstände mit Gewichten wᵢ und Werten vᵢ, '
                'sowie Kapazität W: Gibt es eine Auswahl mit Gesamtgewicht ≤ W '
                'und Gesamtwert ≥ V?'
            ),
            'reduction_from': 'SUBSET-SUM',
            'year': 1972
        },
        {
            'name': 'GRAPH 3-COLORING',
            'description': (
                'Gegeben Graph G: Kann G mit 3 Farben gefärbt werden, '
                'so dass keine zwei benachbarten Knoten die gleiche Farbe haben?'
            ),
            'reduction_from': '3-SAT',
            'year': 1972
        }
    ]


def cook_levin_theorem_sketch() -> dict:
    """
    Skizze des Beweises des Cook-Levin-Satzes: SAT ist NP-vollständig.

    @description
        Beweis-Idee:
        1. SAT ∈ NP: Eine nichtdeterministische TM rät eine Belegung und verifiziert
           in polynomieller Zeit.
        2. SAT ist NP-hart: Für jede Sprache L ∈ NP konstruiere polynomielle Reduktion L ≤_p SAT.

        Kernidee der Reduktion (NP-TM → SAT):
        - Für eine NP-TM M und Eingabe w der Länge n erstelle eine Tabelle (Tableau):
          - Zeile t (0 ≤ t ≤ p(n)): Konfiguration der TM zum Zeitpunkt t
          - Spalte j (0 ≤ j ≤ p(n)): Symbol an Bandposition j
        - Füge Boolesche Variablen C[t,j,a] ein: "Zum Zeitpunkt t enthält Zelle j Symbol a"
        - Kodiere TM-Übergänge als SAT-Klauseln
        - Gesamtgröße der SAT-Formel: O(p(n)²) (polynomiell)

    @return  Dict mit Beweis-Struktur
    """
    return {
        'theorem': 'SAT ist NP-vollständig. (Cook 1971, Levin 1973)',
        'part1': {
            'claim': 'SAT ∈ NP',
            'proof': (
                'Nichtdeterministische TM rät eine Variablenbelegung und '
                'wertet die Formel in polynomieller Zeit aus.'
            )
        },
        'part2': {
            'claim': 'SAT ist NP-hart (∀L ∈ NP: L ≤_p SAT)',
            'tableau_construction': (
                'Gegeben NP-TM M und Eingabe w (|w| = n):\n'
                '1. Tableau T[0..p(n)][0..p(n)]: TM-Konfigurationen\n'
                '2. Variablen C[t, j, a]: "Zur Zeit t ist an Position j Symbol a"\n'
                '3. Klauseln für:\n'
                '   - Eindeutigkeit: Genau ein Symbol pro Zelle und Zeit\n'
                '   - Startbedingung: Konfiguration entspricht Eingabe\n'
                '   - Übergangsrelation: Korrekte TM-Übergänge\n'
                '   - Akzeptanzbedingung: Akzeptierzustand wird erreicht'
            ),
            'formula_size': 'O(p(n)² · |Γ| · |Q|) = polynomiell in n',
            'conclusion': (
                'x ∈ L ⟺ TM M akzeptiert x ⟺ SAT-Formel φ_x erfüllbar.'
            )
        }
    }


def polynomial_hierarchy() -> dict:
    """
    Die Polynomielle Hierarchie PH.

    @description
        Verallgemeinert NP und co-NP durch abwechselnde Quantoren:

        Σᵖ₀ = Πᵖ₀ = P
        Σᵖ₁ = NP  (∃-Orakel)
        Πᵖ₁ = co-NP  (∀-Orakel)
        Σᵖ₂ = NP^NP  (NP mit NP-Orakel)
        Πᵖ₂ = co-NP^NP
        ...
        Σᵖₙ₊₁ = NP^{Σᵖₙ}

        PH = ⋃ₙ Σᵖₙ

        Bekannte Resultate:
        - Falls P = NP → PH = P (Kollaps auf erste Ebene)
        - Falls Σᵖₙ = Πᵖₙ → PH kollapiert auf Ebene n
        - PH ⊆ PSPACE (PH kollabiert falls PH = PSPACE)

    @return  Dict mit Hierachie-Levels und Beispielen
    """
    return {
        'Σᵖ₀ = P': {
            'description': 'Polynomielle Zeit (deterministisch)',
            'examples': ['Sortieren', 'Kürzeste Wege', 'Primalitätstest']
        },
        'Σᵖ₁ = NP': {
            'description': '∃y: P(x, y) mit |y| = poly(|x|)',
            'examples': ['SAT', 'CLIQUE', 'Hamilton-Kreis']
        },
        'Πᵖ₁ = co-NP': {
            'description': '∀y: P(x, y)',
            'examples': ['TAUTOLOGIE', 'UNSAT', 'Kein Hamilton-Kreis']
        },
        'Σᵖ₂ = NP^NP': {
            'description': '∃y₁ ∀y₂: P(x, y₁, y₂)',
            'examples': [
                'MINSAT: Gibt es eine minimale erfüllende Belegung?',
                'Π₁-Theorie von (ℕ, +, ×) (Fragmente)'
            ]
        },
        'Πᵖ₂ = co-NP^NP': {
            'description': '∀y₁ ∃y₂: P(x, y₁, y₂)',
            'examples': ['MINSAT-Komplement']
        },
        'PH': {
            'definition': 'PH = ⋃ₙ Σᵖₙ (Vereinigung aller Ebenen)',
            'collapse_theorems': [
                'P = NP → PH = P',
                'Σᵖₙ = Πᵖₙ → PH kollapiert auf Ebene n',
                'PH = PSPACE → PH hat endliche Ebene (vermutlich nicht)'
            ],
            'open_problem': 'Ist PH = PSPACE? Ist PH = NP?'
        }
    }


# ==============================================================================
# Kolmogorov-Komplexität
# ==============================================================================

def kolmogorov_complexity_approx(string: str) -> dict:
    """
    Obere Schranke für die Kolmogorov-Komplexität K(s) via zlib-Kompression.

    @description
        Die Kolmogorov-Komplexität K(s) ist definiert als:
          K(s) = min{|p| : U(p) = s}
        wobei U eine universelle TM ist und p das kürzeste Programm.

        K(s) ist nicht berechenbar (da Halteproblem nicht entscheidbar).
        Obere Schranke: Jedes Kompressionsprogramm liefert eine obere Schranke.
        Hier: zlib (DEFLATE-Algorithmus) als Kompressor.

        Eigenschaften:
        - K(s) ≤ |s| + c  (triviale obere Schranke: "drucke s")
        - K(s) ≥ K(s) (nicht berechenbar, nur approximierbar)
        - Einfache Strings wie "aaaa..." haben niedrige K-Komplexität
        - Zufällige Strings haben K(s) ≈ |s| (inkompressibel)

    @param string  Eingabezeichenkette
    @return        Dict mit Komplexitätsabschätzung und Metadaten
    """
    raw_bytes      = string.encode('utf-8')
    original_len   = len(raw_bytes)
    compressed     = zlib.compress(raw_bytes, level=9)
    compressed_len = len(compressed)

    # Verhältnis: < 1 bedeutet gut komprimierbar (strukturiert)
    ratio = compressed_len / original_len if original_len > 0 else 1.0

    # Komplexitätsbewertung
    if ratio < 0.3:
        complexity_class = 'sehr niedrig (hochgradig strukturiert)'
    elif ratio < 0.6:
        complexity_class = 'niedrig (strukturiert)'
    elif ratio < 0.9:
        complexity_class = 'mittel'
    else:
        complexity_class = 'hoch (nahezu zufällig / inkompressibel)'

    return {
        'string_length': original_len,
        'compressed_length': compressed_len,
        'compression_ratio': round(ratio, 4),
        'upper_bound_bits': compressed_len * 8,
        'complexity_class': complexity_class,
        'theoretical_note': (
            'K(s) ist nicht berechenbar. '
            'zlib-Länge ist nur eine obere Schranke für K(s). '
            'Die wahre Kolmogorov-Komplexität kann kleiner sein.'
        ),
        'incompressible_definition': (
            f'String gilt als inkompressibel, falls K(s) ≥ |s| - c für kleine Konstante c. '
            f'Aktuell: compressed/original = {ratio:.3f}'
        )
    }


def incompressible_string_demo(n: int = 20) -> dict:
    """
    Demonstriert die Existenz inkompressibler Strings der Länge n.

    @description
        Zählargument: Es gibt 2ⁿ Binärstrings der Länge n, aber nur
        2⁰ + 2¹ + ... + 2^(n-1) = 2ⁿ - 1 < 2ⁿ kürzere Programme.
        Also hat mindestens ein String der Länge n keine kürzere Beschreibung.
        Tatsächlich gilt: Fast alle langen Strings sind inkompressibel.

        Die Menge der Strings mit K(s) < |s| - c hat höchstens
        Anteil 2^(-c) aller Strings der Länge n.

    @param n  Länge der betrachteten Strings
    @return   Dict mit Zählargument und Beispielen
    """
    import random

    # Strukturierter String (niedrige Komplexität)
    structured = 'a' * n
    # Zufälliger String (hohe Komplexität)
    random.seed(42)
    random_str = ''.join(random.choice('ab') for _ in range(n))

    structured_info = kolmogorov_complexity_approx(structured)
    random_info     = kolmogorov_complexity_approx(random_str)

    # Berechne: Wie viele Strings der Länge n sind kompressibel?
    total_strings_of_length_n  = 2 ** n
    shorter_programs_count     = 2 ** n - 1  # Alle Programme der Länge < n

    return {
        'n': n,
        'total_binary_strings': total_strings_of_length_n,
        'shorter_programs_exist': shorter_programs_count,
        'incompressible_exist': (
            f'Mindestens {total_strings_of_length_n - shorter_programs_count} '
            f'String(s) der Länge {n} sind inkompressibel.'
        ),
        'counting_argument': (
            f'Es gibt 2^{n} = {total_strings_of_length_n} Strings der Länge {n}, '
            f'aber nur {shorter_programs_count} < {total_strings_of_length_n} '
            f'Programme kürzerer Länge → Schubfachprinzip.'
        ),
        'density': (
            f'Anteil der Strings mit K(s) < |s| - c: ≤ 2^(-c). '
            f'Bei c=1: ≤ 50% sind um 1 Bit kompressibel.'
        ),
        'structured_example': {
            'string': structured[:20] + ('...' if n > 20 else ''),
            'upper_bound': structured_info['compressed_length'],
            'ratio': structured_info['compression_ratio'],
            'note': 'Viele Wiederholungen → sehr kompressibel → niedrige K-Komplexität'
        },
        'random_example': {
            'string': random_str[:20] + ('...' if n > 20 else ''),
            'upper_bound': random_info['compressed_length'],
            'ratio': random_info['compression_ratio'],
            'note': 'Zufälliger String → kaum kompressibel → hohe K-Komplexität'
        }
    }
