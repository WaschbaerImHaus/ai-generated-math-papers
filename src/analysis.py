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

@author Reisen macht Spass... mit Pia und Dirk e.Kfm.
@date 2026-03-05
"""

import math
import sympy as sp
from typing import Callable, Optional


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
