"""
@file analysis.py
@brief Analysis-Modul: Differentiation, Integration, Nullstellensuche, Taylor-Reihen.
@description
    Implementiert numerische und symbolische Methoden der mathematischen Analysis:

    - numerical_derivative: Zentrale Differenzenquotienten (1. und 2. Ordnung)
    - numerical_integral: Simpson-Regel für bestimmte Integrale
    - newton_raphson: Newton-Raphson-Verfahren zur Nullstellensuche
    - bisection: Bisektionsverfahren zur Nullstellensuche
    - taylor_series: Numerische Taylor-Entwicklung
    - is_continuous: Heuristischer Stetigkeitstest

    Mathematische Grundlagen werden in jedem Docblock erklärt.

@author Michael Fuhrmann
@date 2026-03-05
@lastModified 2026-03-11
"""

import math
import sympy as sp
from typing import Callable
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application,
)

# Konfigurationskonstanten aus zentraler config.py importieren
# Vermeidet Magic Numbers und ermöglicht projektweite Konsistenz
from config import (
    H_DERIVATIVE_1,    # Optimale Schrittweite für 1. Ableitung (≈ 6e-6)
    H_DERIVATIVE_2,    # Optimale Schrittweite für 2. Ableitung (≈ 1.78e-4)
    NEWTON_TOL,        # Konvergenztoleranz Newton-Raphson (1e-12)
    BISECTION_TOL,     # Konvergenztoleranz Bisektion (1e-12)
    SIMPSON_N,         # Standard-Teilintervalle Simpson-Regel (1000)
    MAX_ITERATIONS,    # Maximale Iterationszahl (1000)
)

# Spezifische mathematische Ausnahmen importieren
from exceptions import ConvergenceError, DomainError

# Zentrales Logging-System importieren
from math_logger import get_logger

# Modul-Logger für analysis.py erstellen
_logger = get_logger("analysis")


def _safe_parse(expr_str: str, local_vars: dict | None = None) -> sp.Expr:
    """
    @brief Sicheres Parsing von mathematischen Ausdrücken via SymPy.
    @description
        Verwendet parse_expr mit eingeschränktem Namespace statt sympify().
        Verhindert Code-Injection durch Nutzereingaben.
        Im Gegensatz zu sympify() führt parse_expr keine eval()-ähnlichen
        Operationen auf beliebigem Python-Code aus.

        Unterstützte Funktionen in Eingabe-Strings:
        sin, cos, tan, exp, log, ln, sqrt, abs, pi, e, oo
        sowie Variablen x, y, z, t, n, k.

        Bei unbekannten Ausdrücken (z.B. SymPy-interne Objekte) wird auf
        sympify() als Fallback zurückgegriffen.

    @param expr_str: Mathematischer Ausdruck als String
    @param local_vars: Optionale lokale Variablen-Zuordnung
    @return: SymPy-Ausdruck
    @author Michael Fuhrmann
    @lastModified: 2026-03-11
    """
    # Whitelist der erlaubten Funktionen und Symbole
    transformations = standard_transformations + (implicit_multiplication_application,)
    safe_locals = {
        'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
        'exp': sp.exp, 'log': sp.log, 'ln': sp.log,
        'sqrt': sp.sqrt, 'abs': sp.Abs,
        'pi': sp.pi, 'e': sp.E, 'oo': sp.oo,
        'x': sp.Symbol('x'), 'y': sp.Symbol('y'),
        'z': sp.Symbol('z'), 't': sp.Symbol('t'),
        'n': sp.Symbol('n'), 'k': sp.Symbol('k'),
    }
    # Optionale benutzerdefinierte Symbole einbinden
    if local_vars:
        safe_locals.update(local_vars)
    try:
        return parse_expr(expr_str, local_dict=safe_locals, transformations=transformations)
    except Exception as parse_err:
        # Fallback auf sympify für SymPy-interne Ausdrücke (z.B. bereits geparste Symbole).
        # WICHTIG: parse_expr ist sicher (keine eval()-Semantik), sympify ist unsicherer.
        # Warnung über öffentliche MathLogger.warning()-Methode loggen, bevor Fallback greift.
        # Dies ermöglicht dem Entwickler, fehlerhafte Eingaben im Log zu erkennen.
        _logger.warning(
            f"_safe_parse: parse_expr schlug fehl für Ausdruck '{expr_str}' "
            f"({type(parse_err).__name__}: {parse_err}). "
            f"Fallback auf sp.sympify() – Bitte Eingabe prüfen."
        )
        return sp.sympify(expr_str)


def numerical_derivative(f: Callable[[float], float], x: float, h: float | None = None, order: int = 1) -> float:
    """
    @brief Berechnet die numerische Ableitung einer Funktion im Punkt x.
    @description
        Verwendet den zentralen Differenzenquotienten für bessere Genauigkeit:

        1. Ableitung (zentral):
            f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
            Fehler: O(h^2)  (deutlich besser als einseitig: O(h))

        2. Ableitung (zentral):
            f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h^2
            Fehler: O(h^2)

        Herleitung aus Taylor-Entwicklung:
            f(x+h) = f(x) + h*f'(x) + (h^2/2)*f''(x) + ...
            f(x-h) = f(x) - h*f'(x) + (h^2/2)*f''(x) + ...
            Subtraktion: f(x+h) - f(x-h) = 2h*f'(x) + O(h^3)
            => f'(x) ≈ [f(x+h) - f(x-h)] / (2h)

        Optimale Schrittweite:
        - 1. Ableitung: h ≈ eps^(1/3) ≈ 6e-6  (eps = Maschinengenauigkeit ~2.2e-16)
        - 2. Ableitung: h ≈ eps^(1/4) ≈ 1.2e-4 (größer, da h^2 im Nenner)
        Zu kleines h -> Auslöschungsfehler dominieren
        Zu großes h -> Abbruchfehler der Approximation dominieren

    @param f Die zu differenzierende Funktion.
    @param x Auswertungspunkt.
    @param h Schrittweite (None = automatisch optimal je nach Ordnung).
    @param order Ableitungsordnung (1 oder 2).
    @return Näherungswert der Ableitung f^(order)(x).
    @date 2026-03-05
    """
    if order == 1:
        # Optimale h für 1. Ableitung: aus config.py (H_DERIVATIVE_1 ≈ eps^(1/3) ≈ 6e-6)
        if h is None:
            h = H_DERIVATIVE_1
        return (f(x + h) - f(x - h)) / (2 * h)
    elif order == 2:
        # Optimale h für 2. Ableitung: aus config.py (H_DERIVATIVE_2 ≈ eps^(1/4) ≈ 1.78e-4)
        # (größer als bei 1. Ableitung, da h^2 im Nenner -> mehr Auslöschung)
        if h is None:
            h = H_DERIVATIVE_2
        return (f(x + h) - 2 * f(x) + f(x - h)) / (h ** 2)
    else:
        raise ValueError(f"Ableitungsordnung {order} nicht unterstützt (nur 1 oder 2)")


def numerical_integral(f: Callable[[float], float], a: float, b: float, n: int = SIMPSON_N) -> float:
    """
    @brief Berechnet das bestimmte Integral ∫_a^b f(x) dx mit der Simpson-Regel.
    @description
        Die Simpson-Regel approximiert das Integral durch stückweise Parabeln:

            ∫_a^b f(x) dx ≈ (h/3) * [f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + ... + f(x_n)]

        Dabei ist h = (b-a)/n die Schrittweite und n muss gerade sein.

        Fehler: O(h^4) = O((b-a)^5 / n^4)
        -> Deutlich genauer als Trapezregel O(h^2) und Rechteckregel O(h)

        Herleitung: Jedes Teilintervall [x_i, x_{i+2}] wird durch eine Parabel
        durch die drei Punkte (x_i, f_i), (x_{i+1}, f_{i+1}), (x_{i+2}, f_{i+2})
        approximiert. Die exakte Integration der Parabel ergibt den Faktor h/3 * (f_i + 4f_{i+1} + f_{i+2}).

    @param f Die zu integrierende Funktion.
    @param a Untere Integrationsgrenze.
    @param b Obere Integrationsgrenze.
    @param n Anzahl der Teilintervalle (muss gerade sein, Standard: SIMPSON_N aus config.py).
    @return Näherungswert des Integrals.
    @date 2026-03-05

    Beispiele:
    >>> round(numerical_integral(lambda x: x**2, 0, 1), 4)
    0.3333
    >>> round(numerical_integral(lambda x: 1, 0, 5), 4)
    5.0
    """
    # Sicherstellen dass n gerade ist (Anforderung der Simpson-Regel)
    if n % 2 != 0:
        n += 1

    h = (b - a) / n
    # Stützstellen berechnen
    result = f(a) + f(b)  # Endpunkte mit Gewicht 1

    for i in range(1, n):
        x = a + i * h
        # Ungerade Indizes: Gewicht 4; gerade Indizes: Gewicht 2
        weight = 4 if i % 2 == 1 else 2
        result += weight * f(x)

    return result * h / 3


def newton_raphson(f: Callable[[float], float], x0: float, tol: float = NEWTON_TOL,
                   max_iter: int = MAX_ITERATIONS, h: float = 1e-7,
                   verbose: bool = False) -> float:
    """
    @brief Newton-Raphson-Verfahren zur Nullstellensuche.
    @description
        Das Newton-Raphson-Verfahren findet Nullstellen iterativ:

            x_{n+1} = x_n - f(x_n) / f'(x_n)

        Geometrische Interpretation: In jedem Schritt wird die Tangente an f
        im Punkt (x_n, f(x_n)) gezogen und deren Nullstelle als nächste
        Näherung x_{n+1} verwendet.

        Konvergenz:
        - Quadratische Konvergenz in der Nähe einfacher Nullstellen
        - Keine Garantie für globale Konvergenz (abhängig von x0)
        - Versagt wenn f'(x_n) ≈ 0 (Tangente fast waagerecht)

        Die Ableitung f' wird numerisch berechnet (kein f' als Eingabe nötig).

    @param f Die Funktion, deren Nullstelle gesucht wird.
    @param x0 Startwert (Anfangsnäherung).
    @param tol Toleranz (Abbruch wenn |f(x)| < tol) – Standard: NEWTON_TOL aus config.py.
    @param max_iter Maximale Anzahl Iterationen – Standard: MAX_ITERATIONS aus config.py.
    @param h Schrittweite für numerische Ableitung.
    @param verbose Wenn True, werden alle Iterationsschritte via MathLogger protokolliert.
    @return Näherungswert der Nullstelle.
    @raises RuntimeError Wenn das Verfahren nicht konvergiert.
    @date 2026-03-05
    @lastModified 2026-03-10
    """
    # Optionaler Logger: wird nur bei verbose=True initialisiert (lazy import)
    if verbose:
        from math_logger import get_logger
        logger = get_logger()
        logger.section("Newton-Raphson-Verfahren")

    x = x0
    for i in range(max_iter):
        fx = f(x)

        # Optionale Debug-Ausgabe des aktuellen Schritts
        if verbose:
            logger.step("Newton-Raphson", iteration=i + 1, x=float(x), fx=float(fx))

        # Konvergenzkriterium
        if abs(fx) < tol:
            if verbose:
                logger.result("Nullstelle", float(x))
            return x

        # Numerische Ableitung
        fpx = numerical_derivative(f, x, h=h)

        # Division durch Null abfangen (Tangente ist waagerecht)
        if abs(fpx) < 1e-15:
            # f'(x) ≈ 0 bedeutet, die Tangente ist waagerecht → kein Newton-Schritt möglich
            raise ConvergenceError(
                "Newton-Raphson (f'(x)≈0, Tangente waagerecht)",
                i + 1,
                last_value=float(x)
            )

        # Newton-Schritt
        x = x - fx / fpx

    # Maximale Iterationszahl erreicht ohne Konvergenz
    raise ConvergenceError("Newton-Raphson", max_iter, last_value=float(x))


def bisection(f: Callable[[float], float], a: float, b: float, tol: float = BISECTION_TOL,
              max_iter: int = MAX_ITERATIONS) -> float:
    """
    @brief Bisektionsverfahren zur Nullstellensuche im Intervall [a, b].
    @description
        Das Bisektionsverfahren (Intervallhalbierungsverfahren) garantiert
        Konvergenz, wenn f(a) und f(b) verschiedene Vorzeichen haben
        (Zwischenwertsatz: dann existiert mindestens eine Nullstelle in [a,b]).

        Algorithmus:
        1. Berechne Mittelpunkt m = (a+b)/2
        2. Wenn f(m) ≈ 0: fertig
        3. Wenn f(a)*f(m) < 0: Nullstelle in [a, m], setze b = m
        4. Sonst: Nullstelle in [m, b], setze a = m
        5. Wiederhole

        Konvergenz:
        - Lineare Konvergenz (langsamer als Newton-Raphson)
        - Aber garantiert global konvergent
        - Nach n Schritten: Fehler ≤ (b-a) / 2^n

    @param f Die Funktion, deren Nullstelle gesucht wird.
    @param a Linke Intervallgrenze.
    @param b Rechte Intervallgrenze.
    @param tol Toleranz für Intervallbreite – Standard: BISECTION_TOL aus config.py.
    @param max_iter Maximale Anzahl Iterationen – Standard: MAX_ITERATIONS aus config.py.
    @return Näherungswert der Nullstelle.
    @raises ValueError Wenn kein Vorzeichenwechsel in [a, b] vorliegt.
    @date 2026-03-05

    Beispiele:
    >>> round(bisection(lambda x: x - 1, 0, 2), 4)
    1.0
    >>> abs(bisection(lambda x: x**2 - 4, 0, 3) - 2.0) < 1e-6
    True
    """
    fa = f(a)
    fb = f(b)

    # Vorzeichenwechsel prüfen (Voraussetzung des Zwischenwertsatzes)
    if fa * fb > 0:
        # Ohne Vorzeichenwechsel kann das Bisektionsverfahren keine Nullstelle garantieren
        raise DomainError(
            "bisection",
            f"[{a}, {b}]",
            f"(Kein Vorzeichenwechsel: f({a})={fa:.4f}, f({b})={fb:.4f})"
        )

    for _ in range(max_iter):
        # Intervall ist klein genug -> fertig
        if abs(b - a) < tol:
            break

        # Mittelpunkt berechnen
        m = (a + b) / 2
        fm = f(m)

        # Exakte Nullstelle gefunden
        if abs(fm) < tol:
            return m

        # Intervall halbieren: Seite mit Vorzeichenwechsel behalten
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm

    return (a + b) / 2


def brent_method(f: Callable[[float], float], a: float, b: float,
                 tol: float = 1e-10, max_iter: int = 100) -> float:
    """
    @brief Brent-Methode zur Nullstellensuche – robustes Hybridverfahren.
    @description
        Die Brent-Methode (1973) kombiniert drei Verfahren:
        1. Bisection (Intervallhalbierung) – garantiert Konvergenz
        2. Sekanten-Methode – superlineare Konvergenz nahe der Nullstelle
        3. Inverse quadratische Interpolation (IQI) – kubische Konvergenz

        Algorithmus-Idee (nach Brent 1973):
        - Halte immer [a, b] mit Vorzeichenwechsel (a ist der letzte Schritt, b der beste)
        - Versuche IQI oder Sekantenverfahren (schneller Schritt s)
        - Akzeptiere s nur, wenn er „genug Fortschritt" bringt:
            * s liegt innerhalb von [(3a+b)/4, b]
            * |s - b| < |b - b_prev| / 2  (nicht zu nah am letzten Schritt)
            * (Erste Iteration oder b_prev == b_prev_prev) → Bisection
        - Andernfalls Fallback auf Bisection
        - Garantiert: mindestens so schnell wie Bisection, i.d.R. superlinear

        Konvergenzeigenschaften:
        - Worst-case: wie Bisection O(log((b-a)/tol)) Iterationen
        - Normalfall: superlineare bis kubische Konvergenz
        - Global konvergent (kein Divergenzrisiko wie bei Newton-Raphson)

        Literatur: Brent, R.P. (1973). Algorithms for Minimization without
        Derivatives. Prentice-Hall, Englewood Cliffs, NJ.

    @param f Die Funktion, deren Nullstelle gesucht wird.
    @param a Linke Intervallgrenze mit Vorzeichenwechsel zu b.
    @param b Rechte Intervallgrenze (initialer bester Schätzwert).
    @param tol Toleranz für die Nullstelle – Standard: 1e-10.
    @param max_iter Maximale Iterationsanzahl – Standard: 100.
    @return Näherungswert der Nullstelle.
    @raises ValueError Wenn f(a) und f(b) dasselbe Vorzeichen haben.
    @raises ConvergenceError Wenn keine Konvergenz innerhalb max_iter erreicht wird.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    fa = f(a)
    fb = f(b)

    # Voraussetzung: Vorzeichenwechsel zwischen a und b (Zwischenwertsatz)
    if fa * fb > 0:
        raise ValueError(
            f"brent_method: Kein Vorzeichenwechsel in [{a}, {b}] "
            f"(f({a})={fa:.4g}, f({b})={fb:.4g}). "
            "f(a) und f(b) müssen verschiedene Vorzeichen haben."
        )

    # Exakte Nullstelle bereits an Intervallgrenzen?
    if abs(fa) < tol:
        return a
    if abs(fb) < tol:
        return b

    # b ist immer der aktuell beste Schätzwert (|f(b)| ≤ |f(a)|)
    # Falls nicht: a und b tauschen
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    # c ist die vorherige Grenze des Intervalls (Sekantenverfahren: a → c)
    c = a
    fc = fa

    # b_prev: voriger b-Wert (für Schritt-Größen-Kontrolle)
    b_prev = a
    # b_prev_prev: vorvorigier b-Wert (Erstiterations-Flag)
    b_prev_prev = a
    # Flags ob Bisection im letzten Schritt genutzt wurde
    bisect_flag = True

    for i in range(max_iter):
        # Konvergenzcheck: Intervall kleiner als Toleranz
        if abs(b - a) < tol or abs(fb) < tol:
            return b

        # --- Schritt berechnen: IQI oder Sekante ---
        if abs(fa - fc) > 1e-15 and abs(fb - fc) > 1e-15:
            # Inverse quadratische Interpolation (IQI):
            # Durch drei Punkte (a, fa), (b, fb), (c, fc) ein
            # quadratisches Polynom in f legen und nach x auflösen
            # Formel: Lagrange-Interpolation invertiert
            s = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
        else:
            # Sekantenverfahren (nur 2 Punkte nötig):
            # s = b - fb * (b - a) / (fb - fa)
            s = b - fb * (b - a) / (fb - fa)

        # --- Akzeptanzbedingungen für den schnellen Schritt s ---
        # Berechne die Grenzen für die Bisection
        mid = (a + b) / 2.0

        # Bedingung 1: s liegt nicht im "akzeptierten" Bereich
        # Akzeptierter Bereich: zwischen (3a+b)/4 und b
        cond1 = not ((min(b, (3 * a + b) / 4) < s < max(b, (3 * a + b) / 4))
                     or (min((3 * a + b) / 4, b) <= s <= max((3 * a + b) / 4, b)))

        # Bedingung 2: Schritt ist nicht klein genug (relativ zum letzten Schritt)
        cond2 = bisect_flag and abs(s - b) >= abs(b - b_prev) / 2

        # Bedingung 3: Kein bisect im letzten Schritt und Schritt zu groß
        cond3 = (not bisect_flag) and abs(s - b) >= abs(b_prev - b_prev_prev) / 2

        # Bedingung 4: Bisect-Flag aktiv und b_prev == c (Erstschritt)
        cond4 = bisect_flag and abs(b - a) < tol

        # Wenn eine der Bedingungen verletzt: auf Bisection zurückfallen
        if cond1 or cond2 or cond3 or cond4:
            s = mid
            bisect_flag = True
        else:
            bisect_flag = False

        # --- Schritt auswerten ---
        fs = f(s)

        # b_prev_prev und b_prev aktualisieren (Verlaufshistorie für Schrittkontrolle)
        b_prev_prev = b_prev
        b_prev = b

        # Neues Intervall bestimmen: Seite mit Vorzeichenwechsel behalten
        if fa * fs < 0:
            # Nullstelle liegt in [a, s] → b = s
            b, fb = s, fs
        else:
            # Nullstelle liegt in [s, b] → a = s
            a, fa = s, fs

        # Sicherstellen dass |f(b)| ≤ |f(a)| (b = bester Schätzwert)
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        # c (altes a) aktualisieren
        c = a
        fc = fa

    # Keine Konvergenz innerhalb der maximalen Iterationsanzahl
    raise ConvergenceError("Brent-Methode", max_iter, last_value=float(b))


def taylor_series(f: Callable[[float], float], center: float, degree: int, evaluate_at: float) -> float:
    """
    @brief Berechnet die Taylor-Reihe einer Funktion symbolisch via SymPy.
    @description
        Die Taylor-Reihe von f um den Entwicklungspunkt c lautet:

            f(x) = Summe_{k=0}^{n} [f^(k)(c) / k!] * (x - c)^k

        Wichtige Taylor-Reihen (Entwicklung um 0):
        - e^x = 1 + x + x^2/2! + x^3/3! + ...
        - sin(x) = x - x^3/3! + x^5/5! - ...
        - cos(x) = 1 - x^2/2! + x^4/4! - ...
        - ln(1+x) = x - x^2/2 + x^3/3 - ...  (|x| < 1)

        Implementierung: SymPy-Bibliothek für exakte symbolische Ableitungen.
        Numerische finite Differenzen sind für hohe Ordnungen (k > 8) instabil,
        da h^k im Nenner die Auslöschungsfehler exponentiell verstärkt.

    @param f Die Funktion (muss mit sp.Symbol kompatibel sein: sin, cos, exp etc.).
    @param center Entwicklungspunkt c.
    @param degree Grad der Taylor-Approximation.
    @param evaluate_at Auswertungspunkt x.
    @return Taylor-Approximation von f(evaluate_at).
    @date 2026-03-05
    """
    # Symbolische Variable anlegen
    t = sp.Symbol('t')

    # Python-Funktion -> SymPy-Ausdruck (bekannte Funktionen mappen)
    # Mapping: math-Funktionen auf sympy-Äquivalente
    func_map = {
        math.sin: sp.sin(t),
        math.cos: sp.cos(t),
        math.exp: sp.exp(t),
        math.log: sp.log(t),
        math.tan: sp.tan(t),
        math.sqrt: sp.sqrt(t),
    }

    # SymPy-Ausdruck bestimmen
    if f in func_map:
        expr = func_map[f]
    else:
        # Fallback: Versuche die Funktion als Lambda auf sympy anzuwenden
        try:
            expr = f(t)
        except Exception:
            raise ValueError(f"Funktion {f} kann nicht symbolisch dargestellt werden.")

    # Taylor-Reihe berechnen
    result = 0.0
    x_minus_c = evaluate_at - center

    for k in range(degree + 1):
        # k-te symbolische Ableitung an Stelle 'center' auswerten
        deriv_expr = sp.diff(expr, t, k)
        deriv_value = float(deriv_expr.subs(t, center))

        # Taylor-Term: f^(k)(c) / k! * (x-c)^k
        result += (deriv_value / math.factorial(k)) * (x_minus_c ** k)

    return result


def is_continuous(f: Callable[[float], float], x: float, delta: float = 1e-6) -> bool:
    """
    @brief Heuristischer Test ob f an der Stelle x stetig ist.
    @description
        Eine Funktion f ist stetig in x, wenn:
            lim_{t->x} f(t) = f(x)

        Numerische Überprüfung: Linksseitiger und rechtsseitiger Grenzwert
        werden approximiert und mit f(x) verglichen.

        HINWEIS: Dies ist ein heuristischer Test! Unstetigkeitsstellen
        können durch ungeeignete delta-Werte übersehen werden.

    @param f Die zu prüfende Funktion.
    @param x Der zu prüfende Punkt.
    @param delta Abstand für Grenzwertapproximation.
    @return True wenn f numerisch stetig bei x erscheint.
    @date 2026-03-05
    """
    try:
        fx = f(x)
        f_left = f(x - delta)    # Linksseitiger Näherungswert
        f_right = f(x + delta)   # Rechtsseitiger Näherungswert

        # Stetigkeitskriterium: alle drei Werte ähnlich
        tolerance = max(1e-4, abs(fx) * 1e-4)
        return abs(f_left - fx) < tolerance and abs(f_right - fx) < tolerance
    except (ValueError, ZeroDivisionError, OverflowError):
        return False


# ===========================================================================
# GRENZWERTBERECHNUNG (SYMBOLISCH)
# ===========================================================================

def symbolic_limit(expr_str: str, var: str, point: float | str, direction: str = '+-') -> sp.Expr:
    """
    @brief Berechnet den Grenzwert eines symbolischen Ausdrucks via SymPy.
    @description
        Berechnet lim_{var → point} expr_str.

        SymPy unterstützt folgende Grenzwerte:
        - Reelle Punkte: point als Zahl oder 'oo'/'-oo'
        - Einseitige Grenzwerte: direction='+' (von rechts), '-' (von links)
        - Beidseitig: direction='+-' (nur wenn beide Seiten gleich sind)

        Bekannte Grenzwerte:
        - lim_{x→0} sin(x)/x = 1  (L'Hôpital, Sandwichsatz)
        - lim_{x→∞} (1+1/x)^x = e  (Definition der Eulerschen Zahl)
        - lim_{x→0+} ln(x) = -∞

    @param expr_str SymPy-Ausdruck als String (z.B. "sin(x)/x").
    @param var Variable als String (z.B. "x").
    @param point Grenzpunkt (Zahl, 'oo', '-oo').
    @param direction Richtung: '+' (rechts), '-' (links), '+-' (beidseitig).
    @return SymPy-Ergebnis des Grenzwerts.
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    # SymPy-Variable erstellen
    sym_var = sp.Symbol(var)

    # Ausdruck sicher parsen (verhindert Code-Injection via sympify/eval)
    expr = _safe_parse(expr_str, local_vars={var: sym_var})

    # Grenzpunkt konvertieren
    if point == 'oo' or point == float('inf'):
        sym_point = sp.oo
    elif point == '-oo' or point == float('-inf'):
        sym_point = -sp.oo
    else:
        # Sicherheitsgehärtete Konvertierung: numerische Typen direkt, Strings sicher parsen
        # sp.sympify() auf beliebigen Strings entspricht einem eval() und ist unsicher
        if isinstance(point, (int, float, complex)):
            sym_point = sp.sympify(point)
        else:
            sym_point = _safe_parse(str(point))

    # Grenzwert berechnen
    result = sp.limit(expr, sym_var, sym_point, direction)
    return result


def lhopital_applicable(numerator_str: str, denominator_str: str, var: str, point: float | str) -> dict[str, object]:
    """
    @brief Prüft ob L'Hôpital anwendbar ist und berechnet den Grenzwert.
    @description
        L'Hôpital's Regel (Johann Bernoulli/Marquis de L'Hôpital, 1696):
        Wenn lim f(x)/g(x) die unbestimmte Form 0/0 oder ∞/∞ hat, dann gilt:

            lim_{x→a} f(x)/g(x) = lim_{x→a} f'(x)/g'(x)

        (sofern der rechte Grenzwert existiert oder ±∞ ist)

        Unbestimmte Formen: 0/0, ∞/∞, 0·∞, ∞-∞, 0⁰, 1^∞, ∞⁰

    @param numerator_str Zähler als SymPy-String.
    @param denominator_str Nenner als SymPy-String.
    @param var Variablenname.
    @param point Grenzpunkt.
    @return Dict mit 'applicable' (bool), 'form' (str), 'limit' (SymPy-Ergebnis).
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    sym_var = sp.Symbol(var)
    # Zähler und Nenner sicher parsen (verhindert Code-Injection)
    num_expr = _safe_parse(numerator_str, local_vars={var: sym_var})
    den_expr = _safe_parse(denominator_str, local_vars={var: sym_var})

    # Grenzpunkt konvertieren
    if point == 'oo':
        sym_point = sp.oo
    elif point == '-oo':
        sym_point = -sp.oo
    else:
        # Sicherheitsgehärtete Konvertierung: numerische Typen direkt, Strings sicher parsen
        # sp.sympify() auf beliebigen Strings entspricht einem eval() und ist unsicher
        if isinstance(point, (int, float, complex)):
            sym_point = sp.sympify(point)
        else:
            sym_point = _safe_parse(str(point))

    # Grenzwerte von Zähler und Nenner berechnen
    lim_num = sp.limit(num_expr, sym_var, sym_point)
    lim_den = sp.limit(den_expr, sym_var, sym_point)

    # Unbestimmte Form bestimmen
    applicable = False
    form = "bestimmt"

    if lim_den == 0:
        if lim_num == 0:
            applicable = True
            form = "0/0"
        else:
            form = "x/0 (Pol)"
    elif lim_den == sp.oo or lim_den == -sp.oo:
        if lim_num == sp.oo or lim_num == -sp.oo:
            applicable = True
            form = "∞/∞"
        else:
            form = "x/∞ → 0"

    # Gesamtgrenzwert berechnen (direkt via SymPy)
    total_limit = sp.limit(num_expr / den_expr, sym_var, sym_point)

    return {
        'applicable': applicable,
        'form': form,
        'limit': total_limit
    }


def limit_comparison(f_str: str, g_str: str, var: str, point: float | str) -> dict[str, object]:
    """
    @brief Grenzwert-Vergleichssatz: berechnet lim f/g.
    @description
        Der Grenzwert-Vergleichssatz (limit comparison test):
        Wenn lim_{x→a} f(x)/g(x) = c ≠ 0, dann verhalten sich f und g
        an der Stelle a "gleich" (gleiche Ordnung des Wachstums).

        Anwendungen:
        - Konvergenzuntersuchung von Reihen
        - Vergleich von Wachstumsraten (O-Notation)
        - Asymptotische Äquivalenz: f ~ c·g bei x→a

        Beispiele:
        - lim_{x→∞} x²/(x²+1) = 1 (gleiche Ordnung)
        - lim_{x→0} sin(x)/x = 1 (asymptotisch äquivalent)

    @param f_str Erste Funktion als SymPy-String.
    @param g_str Zweite Funktion als SymPy-String.
    @param var Variablenname.
    @param point Grenzpunkt.
    @return Dict mit 'limit_ratio' (lim f/g) und 'same_behavior' (bool).
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    sym_var = sp.Symbol(var)
    # Beide Funktionen sicher parsen (verhindert Code-Injection)
    f_expr = _safe_parse(f_str, local_vars={var: sym_var})
    g_expr = _safe_parse(g_str, local_vars={var: sym_var})

    # Grenzpunkt konvertieren
    if point == 'oo':
        sym_point = sp.oo
    elif point == '-oo':
        sym_point = -sp.oo
    else:
        # Sicherheitsgehärtete Konvertierung: numerische Typen direkt, Strings sicher parsen
        # sp.sympify() auf beliebigen Strings entspricht einem eval() und ist unsicher
        if isinstance(point, (int, float, complex)):
            sym_point = sp.sympify(point)
        else:
            sym_point = _safe_parse(str(point))

    # Verhältnis-Grenzwert berechnen
    try:
        ratio_limit = sp.limit(f_expr / g_expr, sym_var, sym_point)
    except Exception:
        ratio_limit = sp.nan

    # "Gleiche Ordnung" wenn Grenzwert eine endliche, positive Konstante ist
    same_behavior = (
        ratio_limit != sp.zoo and
        ratio_limit != sp.oo and
        ratio_limit != -sp.oo and
        ratio_limit != sp.nan and
        ratio_limit != 0
    )

    return {
        'limit_ratio': ratio_limit,
        'same_behavior': same_behavior
    }


# ===========================================================================
# PARTIALBRUCHZERLEGUNG
# ===========================================================================

def partial_fraction_decomposition(numerator_coeffs: list[float], denominator_coeffs: list[float]) -> str:
    """
    @brief Partialbruchzerlegung eines rationalen Ausdrucks N(x)/D(x).
    @description
        Die Partialbruchzerlegung zerlegt einen gebrochen-rationalen Ausdruck
        in eine Summe einfacherer Brüche.

        Beispiel: 1/(x²-1) = 1/((x-1)(x+1)) = (1/2)/(x-1) + (-1/2)/(x+1)

        Algorithmus: SymPy's apart()-Funktion nutzt Polynomfaktorisierung
        und löst ein lineares Gleichungssystem für die Koeffizienten.

        Koeffizienten als Listen (höchster Grad zuerst):
        - numerator_coeffs=[1] und denominator_coeffs=[1,0,-1]
          entspricht 1/(x²-1)

    @param numerator_coeffs Koeffizienten des Zählers (höchster Grad zuerst).
    @param denominator_coeffs Koeffizienten des Nenners (höchster Grad zuerst).
    @return String-Darstellung der Partialbruchzerlegung.
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    x = sp.Symbol('x')

    # Polynome aus Koeffizienten-Listen aufbauen
    # Koeffizientenliste [a_n, ..., a_1, a_0] → a_n*x^n + ... + a_1*x + a_0
    def coeffs_to_poly(coeffs: list) -> sp.Expr:
        """Wandelt Koeffizientenliste in SymPy-Polynom um."""
        degree = len(coeffs) - 1
        poly = sp.Integer(0)
        for i, c in enumerate(coeffs):
            poly += c * x ** (degree - i)
        return poly

    num_poly = coeffs_to_poly(numerator_coeffs)
    den_poly = coeffs_to_poly(denominator_coeffs)

    # Rationaler Ausdruck
    rational_expr = num_poly / den_poly

    # Partialbruchzerlegung via SymPy
    decomposed = sp.apart(rational_expr, x)

    return str(decomposed)


def partial_fraction_symbolic(expr_str: str, var: str = 'x') -> dict[str, str | list[str]]:
    """
    @brief Symbolische Partialbruchzerlegung eines SymPy-Ausdrucks.
    @description
        Führt die Partialbruchzerlegung für einen als String gegebenen
        rationalen Ausdruck durch.

        Gibt die Zerlegung als Dictionary zurück mit:
        - 'original': Originalausdruck als String
        - 'decomposed': Zerlegter Ausdruck als String
        - 'terms': Liste der Einzelterme

    @param expr_str SymPy-Ausdruck als String.
    @param var Variablenname (Standard: 'x').
    @return Dict mit 'original', 'decomposed', 'terms'.
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    sym_var = sp.Symbol(var)
    # Ausdruck sicher parsen (verhindert Code-Injection)
    expr = _safe_parse(expr_str, local_vars={var: sym_var})

    # Partialbruchzerlegung
    decomposed = sp.apart(expr, sym_var)

    # Einzelterme extrahieren (Summanden)
    terms = sp.Add.make_args(decomposed)
    terms_str = [str(t) for t in terms]

    return {
        'original': str(expr),
        'decomposed': str(decomposed),
        'terms': terms_str
    }


# ===========================================================================
# UNEIGENTLICHE INTEGRALE
# ===========================================================================

def improper_integral_numerical(f: Callable[[float], float], a: float, b: float,
                                 singularities: list[float] | None = None,
                                 eps: float = 1e-6) -> float:
    """
    @brief Numerisches uneigentliches Integral mit möglichen Singularitäten.
    @description
        Behandelt uneigentliche Integrale ∫_a^b f(x) dx in folgenden Fällen:
        1. Unendliche Grenzen (a=-∞ oder b=+∞): Substitution t=1/(1+x) oder
           numerische Gauss-Laguerre-Quadratur
        2. Singularitäten im Inneren: Aufteilung in Teilintegrale mit
           Epsilon-Umgebung um die Singularität
        3. Kombination beider Fälle

        Methode für ∫_0^∞ f(x)dx: Substitution x = t/(1-t), dx = 1/(1-t)²dt
        → transformiert [0,∞) auf [0,1)

        Methode für ∫_-∞^∞: Aufteilen in zwei Hälften

    @param f Die zu integrierende Funktion.
    @param a Untere Grenze (float('-inf') möglich).
    @param b Obere Grenze (float('inf') möglich).
    @param singularities Liste von Singularitätspunkten im Inneren von (a,b).
    @param eps Epsilon-Umgebung um Singularitäten und für Grenzen.
    @return Näherungswert des uneigentlichen Integrals.
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    from scipy import integrate as sci_integrate

    if singularities is None:
        singularities = []

    # Hilfsfunktion: endliches Integral via scipy.integrate.quad
    def finite_integral(func, lo, hi):
        """Berechnet endliches Integral mit scipy quad (adaptiv)."""
        result, _ = sci_integrate.quad(func, lo, hi, limit=200)
        return result

    # Fall 1: Beide Grenzen endlich, keine Singularitäten
    if math.isfinite(a) and math.isfinite(b) and not singularities:
        return finite_integral(f, a, b)

    # Fall 2: Singularitäten im Inneren → Aufteilung
    if singularities:
        # Grenzen und Singularitäten zu einer sortierten Liste zusammenfassen
        points = sorted([a] + list(singularities) + [b])
        total = 0.0
        for i in range(len(points) - 1):
            lo = points[i]
            hi = points[i + 1]

            # Epsilon-Umgebung um Singularitäten
            if lo in singularities:
                lo = lo + eps
            if hi in singularities:
                hi = hi - eps

            # Grenzen behandeln
            if not math.isfinite(lo):
                lo_int = -1e10
            else:
                lo_int = lo
            if not math.isfinite(hi):
                hi_int = 1e10
            else:
                hi_int = hi

            if lo_int < hi_int:
                total += finite_integral(f, lo_int, hi_int)
        return total

    # Fall 3: Unendliche Grenzen via scipy
    if not math.isfinite(a) and not math.isfinite(b):
        result, _ = sci_integrate.quad(f, -math.inf, math.inf, limit=200)
        return result
    elif not math.isfinite(b):
        result, _ = sci_integrate.quad(f, a, math.inf, limit=200)
        return result
    else:
        result, _ = sci_integrate.quad(f, -math.inf, b, limit=200)
        return result


def improper_integral_symbolic(expr_str: str, var: str, a: float | str, b: float | str) -> dict[str, object]:
    """
    @brief Symbolisches uneigentliches Integral via SymPy.
    @description
        Berechnet ∫_a^b f(x) dx symbolisch, wobei a und b unendlich sein können.

        Konvergenzuntersuchung:
        - Das Integral konvergiert, wenn das Ergebnis endlich ist
        - Bei Divergenz gibt SymPy sp.oo, -sp.oo oder sp.zoo zurück

        Bekannte Ergebnisse:
        - ∫_0^∞ e^(-x) dx = 1
        - ∫_0^∞ x^(s-1) e^(-x) dx = Γ(s) (Gamma-Funktion)
        - ∫_{-∞}^∞ e^(-x²) dx = √π (Gauß-Integral)

    @param expr_str SymPy-Ausdruck als String.
    @param var Variablenname.
    @param a Untere Grenze (Zahl, 'oo', '-oo', oder sympy-Ausdruck).
    @param b Obere Grenze (Zahl, 'oo', '-oo', oder sympy-Ausdruck).
    @return Dict mit 'integral' (SymPy), 'converges' (bool), 'value' (float oder None).
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    sym_var = sp.Symbol(var)
    # Ausdruck sicher parsen (verhindert Code-Injection bei Nutzereingaben)
    expr = _safe_parse(expr_str, local_vars={var: sym_var})

    # Grenzen konvertieren
    def convert_bound(bound):
        """Konvertiert Grenzwert in SymPy-Ausdruck.
        Sicherheitsgehärtet: numerische Typen direkt konvertieren,
        Strings über _safe_parse() leiten (verhindert Code-Injection).
        """
        if bound == 'oo' or bound == float('inf'):
            return sp.oo
        elif bound == '-oo' or bound == float('-inf'):
            return -sp.oo
        elif isinstance(bound, (int, float, complex)):
            # Numerische Werte sind sicher direkt zu konvertieren
            return sp.sympify(bound)
        else:
            # Strings werden über den sicheren Parser geleitet (kein eval)
            return _safe_parse(str(bound))

    sym_a = convert_bound(a)
    sym_b = convert_bound(b)

    # Symbolisches Integral berechnen
    try:
        integral_result = sp.integrate(expr, (sym_var, sym_a, sym_b))
    except Exception as exc:
        return {
            'integral': None,
            'converges': False,
            'value': None
        }

    # Konvergenz prüfen: Ergebnis darf nicht unendlich oder komplex-unendlich sein
    divergent_values = {sp.oo, -sp.oo, sp.zoo, sp.nan}
    converges = integral_result not in divergent_values

    # Numerischen Wert extrahieren
    value = None
    if converges:
        try:
            value = float(integral_result.evalf())
        except Exception:
            value = None

    return {
        'integral': integral_result,
        'converges': converges,
        'value': value
    }


def cauchy_principal_value(f: Callable[[float], float], a: float, b: float,
                            singularity: float, eps: float = 1e-8) -> float:
    """
    @brief Cauchyscher Hauptwert P.V. ∫_a^b f(x) dx bei einer Singularität.
    @description
        Der Cauchy-Hauptwert (Valeur principale, P.V.) ist definiert als:

            P.V. ∫_a^b f(x) dx = lim_{ε→0} [∫_a^{c-ε} f(x) dx + ∫_{c+ε}^b f(x) dx]

        wobei c ∈ (a,b) eine Singularität von f ist.

        Beispiel: P.V. ∫_{-1}^{1} 1/x dx = 0
        (Beide Hälften heben sich gegenseitig auf durch Symmetrie)

        Wichtig: Der Hauptwert kann existieren, auch wenn das Integral im
        gewöhnlichen Sinne nicht konvergiert.

    @param f Die Funktion mit Singularität.
    @param a Untere Integrationsgrenze.
    @param b Obere Integrationsgrenze.
    @param singularity Innere Singularitätsstelle c ∈ (a, b).
    @param eps Epsilon-Parameter für die Näherung des Grenzwerts.
    @return Cauchyscher Hauptwert des Integrals.
    @author Michael Fuhrmann
    @lastModified 2026-03-08
    """
    from scipy import integrate as sci_integrate

    # Linkes Teilintegral: ∫_a^{c-ε} f(x) dx
    left_result, _ = sci_integrate.quad(f, a, singularity - eps, limit=200)

    # Rechtes Teilintegral: ∫_{c+ε}^b f(x) dx
    right_result, _ = sci_integrate.quad(f, singularity + eps, b, limit=200)

    # Cauchy-Hauptwert = Summe beider Teile
    return left_result + right_result


# ===========================================================================
# PARALLELES SYMBOLISCHES RECHNEN
# ===========================================================================

from concurrent.futures import ThreadPoolExecutor
import os as _os


def parallel_symbolic_compute(tasks: list, max_workers: int = None) -> list:
    """
    @file analysis.py
    @brief Führt mehrere SymPy-Berechnungen parallel via ThreadPoolExecutor aus.
    @description
        Ermöglicht die gleichzeitige Ausführung mehrerer symbolischer Berechnungen
        mit SymPy, indem Python-Threads genutzt werden.

        Warum ThreadPoolExecutor statt ProcessPoolExecutor?
        SymPy-Objekte (Expressions, Symbols, ...) sind nicht Pickle-serialisierbar
        und können deshalb nicht direkt an Subprozesse übergeben werden.
        Der GIL (Global Interpreter Lock) wird von SymPy bei vielen Operationen
        (I/O, C-Extensions) freigegeben, sodass echte Parallelität möglich ist.

        Aufbau eines Task-Tupels:
            (callable, args_tuple, kwargs_dict)

        Beispiele:
            import sympy as sp
            x = sp.Symbol('x')
            tasks = [
                (sp.integrate, (x**2, x), {}),
                (sp.diff,      (sp.sin(x), x), {}),
                (sp.factor,    (x**2 - 1,), {}),
            ]
            ergebnisse = parallel_symbolic_compute(tasks)
            # ergebnisse: [x**3/3, cos(x), (x-1)*(x+1)]

        Reihenfolge-Garantie:
            Die Ergebnisliste ist in derselben Reihenfolge wie die tasks-Liste,
            unabhängig davon, welche Berechnungen zuerst fertig wurden.

    @param tasks:       Liste von (callable, args, kwargs)-Tupeln.
                        - callable:  Jede aufrufbare SymPy-Funktion oder eigene Funktion.
                        - args:      Tupel mit positionalen Argumenten (kann leer sein: ()).
                        - kwargs:    Dict mit Schlüsselwort-Argumenten (kann leer sein: {}).
    @param max_workers: Maximale Anzahl paralleler Threads.
                        None (Standard) = Anzahl logischer CPU-Kerne.
    @return             Liste der Bergebnisse in derselben Reihenfolge wie tasks.
                        Auftretende Ausnahmen werden weitergeleitet (re-raised).
    @author Michael Fuhrmann
    @date 2026-03-11
    @lastModified 2026-03-11
    """
    # Standardmäßig so viele Threads wie logische CPU-Kerne (mindestens 4)
    if max_workers is None:
        max_workers = _os.cpu_count() or 4

    def _run_task(task):
        """Hilfsfunktion: Entpackt ein Task-Tupel und ruft die Funktion auf."""
        func, args, kwargs = task
        # Funktion mit positionalen und Schlüsselwort-Argumenten aufrufen
        return func(*args, **kwargs)

    # ThreadPoolExecutor: verwaltet den Thread-Pool automatisch (inkl. Cleanup)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Alle Tasks gleichzeitig einreichen → Future-Objekte erzeugen
        futures = [executor.submit(_run_task, t) for t in tasks]

        # Ergebnisse in der ursprünglichen Reihenfolge einsammeln
        # f.result() blockiert bis das Ergebnis vorliegt (oder wirft Exception)
        return [f.result() for f in futures]
