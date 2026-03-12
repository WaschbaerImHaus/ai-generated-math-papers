"""
@file parallel_compute.py
@brief Parallele symbolische und numerische Berechnungen via ProcessPoolExecutor.
@description
    Dieses Modul ermöglicht die parallelisierte Ausführung von mathematischen
    Berechnungen, die sonst sequenziell ausgeführt würden.

    Implementiert:
    - parallel_symbolic: Allgemeine Parallelisierung beliebiger Funktionen
    - parallel_integrate: Parallele symbolische Integration mehrerer Ausdrücke
    - parallel_limit: Parallele Grenzwertberechnung mehrerer Ausdrücke
    - parallel_derivative: Parallele Ableitungsberechnung mehrerer Ausdrücke

    Warum Parallelisierung?
        Symbolische Berechnungen (SymPy) sind oft CPU-gebunden und unabhängig
        voneinander. Mit ProcessPoolExecutor (Multi-Prozess, nicht Threads)
        können wir echte Parallelität nutzen – Python's GIL blockiert hier nicht,
        da jeder Prozess seinen eigenen Interpreter hat.

    Wichtige Einschränkungen:
        - Funktionen müssen pickle-bar sein (keine Lambda-Ausdrücke als func)
        - SymPy-Ausdrücke müssen als Strings oder Symbole übergeben werden
        - Overhead durch Prozess-Start: lohnt sich erst ab ~0.1s pro Aufgabe

@author Michael Fuhrmann
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Optional
import re as _re_pc
import sympy as sp
from sympy.parsing.sympy_parser import (
    parse_expr as _parse_expr_pc,
    standard_transformations as _std_tr_pc,
    implicit_multiplication_application as _impl_mult_pc,
)

# Erlaubte Transformationen für sicheres Parsing
_SAFE_TR_PC = _std_tr_pc + (_impl_mult_pc,)

# Whitelist-Muster: Nur mathematische Ausdrücke erlaubt
_SAFE_CHARS_PC = _re_pc.compile(r'^[a-zA-Z0-9\s\+\-\*\/\^\(\)\.\,\_\=\<\>]+$')


def _safe_sympify(s: str) -> sp.Basic:
    """
    @brief Parst einen Ausdruck sicher mit parse_expr statt sympify.
    @description
        Validiert Länge und Zeichenmenge, bevor der Ausdruck geparst wird.
        Verhindert Code-Injection über sympify().
    @param s: Mathematischer Ausdruck als String
    @return SymPy-Ausdruck
    @raises ValueError bei ungültigem Ausdruck
    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """
    if len(s) > 1000 or not _SAFE_CHARS_PC.match(s):
        raise ValueError(f"Ungültiger Ausdruck (unerlaubte Zeichen oder zu lang): '{s}'")
    return _parse_expr_pc(s, transformations=_SAFE_TR_PC)


# ===========================================================================
# HILFSFUNKTIONEN (müssen auf Modulebene sein, damit sie pickle-bar sind)
# ===========================================================================

def _integrate_single(args: tuple) -> tuple:
    """
    Berechnet ein einzelnes symbolisches Integral (Hilfsfunktion für Parallelisierung).

    Muss auf Modulebene definiert sein, damit ProcessPoolExecutor sie in
    Kindprozesse serialisieren (pickle) kann.

    @param args: Tupel (expr_str, x_str) mit SymPy-Ausdruck und Variablenname als Strings
    @return: Tupel (expr_str, integral_str) – Ergebnis als String
    @lastModified 2026-03-10
    """
    expr_str, x_str = args
    # Ausdruck und Variable in Kindprozess neu aufbauen (pickle-Sicherheit)
    x = sp.Symbol(x_str)
    expr = _safe_sympify(expr_str)
    result = sp.integrate(expr, x)
    return (expr_str, str(result))


def _limit_single(args: tuple) -> tuple:
    """
    Berechnet einen einzelnen Grenzwert (Hilfsfunktion für Parallelisierung).

    @param args: Tupel (expr_str, x_str, point_str, direction)
    @return: Tupel (expr_str, limit_str) – Ergebnis als String
    @lastModified 2026-03-10
    """
    expr_str, x_str, point_str, direction = args
    x = sp.Symbol(x_str)
    expr = _safe_sympify(expr_str)
    # Punkt parsen (kann 'oo', '-oo' oder eine Zahl sein)
    point = _safe_sympify(point_str)
    result = sp.limit(expr, x, point, direction)
    return (expr_str, str(result))


def _derivative_single(args: tuple) -> tuple:
    """
    Berechnet eine einzelne symbolische Ableitung (Hilfsfunktion für Parallelisierung).

    @param args: Tupel (expr_str, x_str, order)
    @return: Tupel (expr_str, derivative_str) – Ergebnis als String
    @lastModified 2026-03-10
    """
    expr_str, x_str, order = args
    x = sp.Symbol(x_str)
    expr = _safe_sympify(expr_str)
    result = sp.diff(expr, x, order)
    return (expr_str, str(result))


# ===========================================================================
# ÖFFENTLICHE API
# ===========================================================================

def parallel_symbolic(
    func: Callable,
    args_list: list[tuple],
    max_workers: int = 4
) -> dict[tuple, Any]:
    """
    Führt eine beliebige Funktion parallelisiert für eine Liste von Argumenten aus.

    Die Funktion muss pickle-bar sein (also auf Modulebene definiert, kein Lambda).
    Die Ergebnisse werden als Dictionary zurückgegeben, wobei die Schlüssel die
    ursprünglichen Argument-Tupel sind.

    Parallelisierungsstrategie:
        - ProcessPoolExecutor statt ThreadPoolExecutor (kein GIL-Problem)
        - as_completed() verarbeitet Ergebnisse sofort, sobald sie fertig sind
        - Bei Fehler in einem Future: Ausnahme wird in results[args] gespeichert

    Beispiel:
        >>> import sympy as sp
        >>> def quad(args): x_val, = args; return x_val ** 2
        >>> parallel_symbolic(quad, [(2,), (3,), (4,)], max_workers=2)
        {(2,): 4, (3,): 9, (4,): 16}

    @param func: Callable das genau ein Argument (Tupel) annimmt
    @param args_list: Liste von Argument-Tupeln (eines pro Aufruf)
    @param max_workers: Maximale Anzahl paralleler Prozesse (Standard: 4)
    @return: Dictionary {args_tuple: Ergebnis} für alle Eingaben
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    results: dict[tuple, Any] = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Alle Aufgaben einreichen, Future → args Zuordnung merken
        futures = {executor.submit(func, args): args for args in args_list}

        # Ergebnisse einsammeln, sobald sie fertig sind
        for future in as_completed(futures):
            args = futures[future]
            try:
                # Strukturiertes Ergebnis: {"success": True, "value": ..., "error": None}
                results[args] = {"success": True, "value": future.result(), "error": None}
            except Exception as exc:
                # Strukturierter Fehler: {"success": False, "value": None, "error": str}
                # Ermöglicht Callers, Fehler eindeutig von Ergebnissen zu unterscheiden.
                results[args] = {"success": False, "value": None, "error": str(exc)}

    return results


def parallel_integrate(
    exprs: list[str],
    x: str = 'x',
    max_workers: int = 4
) -> dict[str, str]:
    """
    Integriert mehrere symbolische Ausdrücke parallel.

    Jeder Ausdruck wird in einem separaten Prozess symbolisch integriert.
    Die Ausdrücke werden als Strings übergeben und zurückgegeben (pickle-sicher).

    Mathematisch:
        Berechnet ∫ expr dx für jeden Ausdruck in der Liste.
        Ergebnis ist das unbestimmte Integral (ohne Konstante C).

    Beispiel:
        >>> parallel_integrate(['x**2', 'sin(x)', 'exp(x)'], x='x')
        {'x**2': 'x**3/3', 'sin(x)': '-cos(x)', 'exp(x)': 'exp(x)'}

    @param exprs: Liste von SymPy-Ausdrücken als Strings
    @param x: Name der Integrationsvariable (Standard: 'x')
    @param max_workers: Maximale Anzahl paralleler Prozesse
    @return: Dictionary {ausdruck_str: integral_str}
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Argument-Tupel für _integrate_single vorbereiten
    args_list = [(expr_str, x) for expr_str in exprs]
    raw_results = parallel_symbolic(_integrate_single, args_list, max_workers)

    # Ergebnisse in benutzerfreundliches Format umwandeln
    result: dict[str, str] = {}
    for args, entry in raw_results.items():
        expr_str = args[0]
        # Strukturiertes Ergebnis auswerten: {"success": bool, "value": ..., "error": ...}
        if not entry["success"]:
            result[expr_str] = f"Fehler: {entry['error']}"
        else:
            # value ist Tupel (expr_str, integral_str)
            result[expr_str] = entry["value"][1]

    return result


def parallel_limit(
    exprs: list[str],
    x: str = 'x',
    point: str = '0',
    direction: str = '+',
    max_workers: int = 4
) -> dict[str, str]:
    """
    Berechnet Grenzwerte mehrerer symbolischer Ausdrücke parallel.

    Mathematisch:
        Berechnet lim_{x→point} expr für jeden Ausdruck in der Liste.

    Beispiel:
        >>> parallel_limit(['sin(x)/x', '(1-cos(x))/x**2'], x='x', point='0')
        {'sin(x)/x': '1', '(1-cos(x))/x**2': '1/2'}

    @param exprs: Liste von SymPy-Ausdrücken als Strings
    @param x: Name der Grenzwertvariable
    @param point: Grenzwertpunkt als String (z.B. '0', 'oo', 'pi')
    @param direction: Richtung '+' (von rechts), '-' (von links), '+-' (beidseitig)
    @param max_workers: Maximale Anzahl paralleler Prozesse
    @return: Dictionary {ausdruck_str: grenzwert_str}
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Argument-Tupel für _limit_single vorbereiten
    args_list = [(expr_str, x, point, direction) for expr_str in exprs]
    raw_results = parallel_symbolic(_limit_single, args_list, max_workers)

    result: dict[str, str] = {}
    for args, entry in raw_results.items():
        expr_str = args[0]
        # Strukturiertes Ergebnis auswerten: {"success": bool, "value": ..., "error": ...}
        if not entry["success"]:
            result[expr_str] = f"Fehler: {entry['error']}"
        else:
            result[expr_str] = entry["value"][1]

    return result


def parallel_derivative(
    exprs: list[str],
    x: str = 'x',
    order: int = 1,
    max_workers: int = 4
) -> dict[str, str]:
    """
    Berechnet Ableitungen mehrerer symbolischer Ausdrücke parallel.

    Mathematisch:
        Berechnet d^order/dx^order expr für jeden Ausdruck in der Liste.

    Beispiel:
        >>> parallel_derivative(['x**3', 'sin(x)', 'exp(2*x)'], order=2)
        {'x**3': '6*x', 'sin(x)': '-sin(x)', 'exp(2*x)': '4*exp(2*x)'}

    @param exprs: Liste von SymPy-Ausdrücken als Strings
    @param x: Name der Ableitungsvariable
    @param order: Ableitungsordnung (1 = erste Ableitung, 2 = zweite, etc.)
    @param max_workers: Maximale Anzahl paralleler Prozesse
    @return: Dictionary {ausdruck_str: ableitung_str}
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Argument-Tupel für _derivative_single vorbereiten
    args_list = [(expr_str, x, order) for expr_str in exprs]
    raw_results = parallel_symbolic(_derivative_single, args_list, max_workers)

    result: dict[str, str] = {}
    for args, entry in raw_results.items():
        expr_str = args[0]
        # Strukturiertes Ergebnis auswerten: {"success": bool, "value": ..., "error": ...}
        if not entry["success"]:
            result[expr_str] = f"Fehler: {entry['error']}"
        else:
            result[expr_str] = entry["value"][1]

    return result
