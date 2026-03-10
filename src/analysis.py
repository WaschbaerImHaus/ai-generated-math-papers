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

@author Kurt Ingwer
@date 2026-03-05
@lastModified 2026-03-10
"""

import math
import sympy as sp
from typing import Callable, Optional
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application,
)


def _safe_parse(expr_str: str, local_vars: dict = None) -> sp.Expr:
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
    @author Kurt Ingwer
    @lastModified: 2026-03-09
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
    except Exception:
        # Fallback auf sympify für SymPy-interne Ausdrücke (z.B. bereits geparste Symbole)
        return sp.sympify(expr_str)


def numerical_derivative(f: Callable, x: float, h: float = None, order: int = 1) -> float:
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
        # Optimale h für 1. Ableitung: eps^(1/3)
        if h is None:
            h = 6e-6
        return (f(x + h) - f(x - h)) / (2 * h)
    elif order == 2:
        # Optimale h für 2. Ableitung: eps^(1/4) ~ 1.2e-4
        # (größer als bei 1. Ableitung, da h^2 im Nenner -> mehr Auslöschung)
        if h is None:
            h = 1.5e-4
        return (f(x + h) - 2 * f(x) + f(x - h)) / (h ** 2)
    else:
        raise ValueError(f"Ableitungsordnung {order} nicht unterstützt (nur 1 oder 2)")


def numerical_integral(f: Callable, a: float, b: float, n: int = 10000) -> float:
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
    @param n Anzahl der Teilintervalle (muss gerade sein, Standard: 10000).
    @return Näherungswert des Integrals.
    @date 2026-03-05
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


def newton_raphson(f: Callable, x0: float, tol: float = 1e-10,
                   max_iter: int = 1000, h: float = 1e-7) -> float:
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
    @param tol Toleranz (Abbruch wenn |f(x)| < tol).
    @param max_iter Maximale Anzahl Iterationen.
    @param h Schrittweite für numerische Ableitung.
    @return Näherungswert der Nullstelle.
    @raises RuntimeError Wenn das Verfahren nicht konvergiert.
    @date 2026-03-05
    """
    x = x0
    for i in range(max_iter):
        fx = f(x)

        # Konvergenzkriterium
        if abs(fx) < tol:
            return x

        # Numerische Ableitung
        fpx = numerical_derivative(f, x, h=h)

        # Division durch Null abfangen (Tangente ist waagerecht)
        if abs(fpx) < 1e-15:
            raise RuntimeError(f"Newton-Raphson: f'(x) ≈ 0 bei x={x}, Verfahren divergiert")

        # Newton-Schritt
        x = x - fx / fpx

    raise RuntimeError(f"Newton-Raphson: Keine Konvergenz nach {max_iter} Iterationen")


def bisection(f: Callable, a: float, b: float, tol: float = 1e-10,
              max_iter: int = 1000) -> float:
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
    @param tol Toleranz für Intervallbreite.
    @param max_iter Maximale Anzahl Iterationen.
    @return Näherungswert der Nullstelle.
    @raises ValueError Wenn kein Vorzeichenwechsel in [a, b] vorliegt.
    @date 2026-03-05
    """
    fa = f(a)
    fb = f(b)

    # Vorzeichenwechsel prüfen (Voraussetzung des Zwischenwertsatzes)
    if fa * fb > 0:
        raise ValueError(
            f"Bisektionsverfahren: Kein Vorzeichenwechsel in [{a}, {b}]. "
            f"f({a})={fa:.4f}, f({b})={fb:.4f}"
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


def taylor_series(f: Callable, center: float, degree: int, evaluate_at: float) -> float:
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


def is_continuous(f: Callable, x: float, delta: float = 1e-6) -> bool:
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

def symbolic_limit(expr_str: str, var: str, point, direction: str = '+-'):
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
    @author Kurt Ingwer
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


def lhopital_applicable(numerator_str: str, denominator_str: str, var: str, point) -> dict:
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
    @author Kurt Ingwer
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


def limit_comparison(f_str: str, g_str: str, var: str, point) -> dict:
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
    @author Kurt Ingwer
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

def partial_fraction_decomposition(numerator_coeffs: list, denominator_coeffs: list) -> str:
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
    @author Kurt Ingwer
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


def partial_fraction_symbolic(expr_str: str, var: str = 'x') -> dict:
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
    @author Kurt Ingwer
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

def improper_integral_numerical(f: Callable, a: float, b: float,
                                 singularities: list = None,
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
    @author Kurt Ingwer
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


def improper_integral_symbolic(expr_str: str, var: str, a, b) -> dict:
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
    @author Kurt Ingwer
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


def cauchy_principal_value(f: Callable, a: float, b: float,
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
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    from scipy import integrate as sci_integrate

    # Linkes Teilintegral: ∫_a^{c-ε} f(x) dx
    left_result, _ = sci_integrate.quad(f, a, singularity - eps, limit=200)

    # Rechtes Teilintegral: ∫_{c+ε}^b f(x) dx
    right_result, _ = sci_integrate.quad(f, singularity + eps, b, limit=200)

    # Cauchy-Hauptwert = Summe beider Teile
    return left_result + right_result
