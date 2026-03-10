"""
@file config.py
@brief Globale Konfigurationskonstanten für das specialist-maths Projekt.
@description
    Diese Datei zentralisiert alle numerischen Parameter, Toleranzwerte
    und Algorithmus-Konfigurationen, die in verschiedenen Modulen
    verwendet werden. Durch die Zentralisierung können Werte an einer
    einzigen Stelle angepasst werden, was Konsistenz sicherstellt
    und doppelte Definitionen vermeidet.

    Enthält:
    - Maschinengenauigkeit und optimale Schrittweiten für Differentiation
    - Konvergenztoleranzen für iterative Verfahren
    - Algorithmus-Parameter (Iterationsgrenzen, Terme-Anzahl)
    - Lanczos-Koeffizienten für die Gamma-Funktion
    - Mathematische Konstanten (Goldener Schnitt)
    - Metadaten (Autor, Version)

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import math
import sympy as sp
from typing import Union

# ===========================================================================
# MASCHINENGENAUIGKEIT UND SCHRITTWEITEN
# ===========================================================================

# Maschinengenauigkeit für IEEE 754 double precision
# eps = 2^-52 ≈ 2.22e-16; wir verwenden 1e-15 als praktischen Wert
EPSILON: float = 1e-15

# Optimale Schrittweite für den zentralen Differenzenquotienten der 1. Ableitung
# Herleitung: Minimiere Rundungs- + Abbruchfehler → h ≈ eps^(1/3) ≈ 6.06e-6
# Zu klein: Auslöschung; zu groß: Taylorabbruchfehler dominiert
H_DERIVATIVE_1: float = EPSILON ** (1.0 / 3.0)  # ≈ 1e-5

# Optimale Schrittweite für die 2. Ableitung (größer wegen h² im Nenner)
# Herleitung: h ≈ eps^(1/4) ≈ 1.78e-4
# Der h²-Nenner verstärkt Auslöschungsfehler stärker als beim 1. Differenzenquotienten
H_DERIVATIVE_2: float = EPSILON ** (1.0 / 4.0)  # ≈ 1.78e-4

# ===========================================================================
# ITERATIONSPARAMETER
# ===========================================================================

# Standard-Maximale Anzahl an Iterationen für iterative Verfahren
# Gilt für Newton-Raphson, Bisektions-, Optimierungs- und Konvergenzalgorithmen
MAX_ITERATIONS: int = 1000

# Konvergenztoleranz für das Newton-Raphson-Verfahren
# Abbruch wenn |f(x_n)| < NEWTON_TOL oder |x_{n+1} - x_n| < NEWTON_TOL
# Strenger als Bisektion, da Newton quadratisch konvergiert
NEWTON_TOL: float = 1e-12

# Konvergenztoleranz für das Bisektionsverfahren
# Abbruch wenn Intervallbreite |b - a| < BISECTION_TOL
# Garantierte lineare Konvergenz: nach n Schritten Fehler ≤ (b-a)/2^n
BISECTION_TOL: float = 1e-12

# ===========================================================================
# INTEGRATIONSPARAMETER
# ===========================================================================

# Standard-Anzahl der Teilintervalle für die Simpson-Regel
# Muss gerade sein. Fehler: O((b-a)^5 / n^4) → bei n=1000 sehr klein
SIMPSON_N: int = 1000

# ===========================================================================
# PRIMZAHL-ALGORITHMEN
# ===========================================================================

# Deterministische Miller-Rabin-Zeugen für alle Zahlen < 3.317.044.064.679.887.385.961.981
# Diese Zeugenmengen garantieren korrekte Ergebnisse ohne Probabilistik.
# Quelle: Bach & Sorenson (1993), erweitert durch Pomerance et al.
# Für n < 3.215.031.751 reichen: [2, 3, 5, 7]
# Für n < 3.317.044.064.679.887.385.961.981: vollständige Liste unten
PRIME_MILLER_RABIN_WITNESSES: list[int] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

# ===========================================================================
# TAYLOR-REIHEN
# ===========================================================================

# Standard-Anzahl der Terme für Taylor-Reihenentwicklungen
# Zu viele Terme: symbolisch langsam; zu wenige: schlechte Approximation
TAYLOR_DEFAULT_TERMS: int = 10

# ===========================================================================
# RIEMANN-ZETA-FUNKTION / EULER-KNOPP-BESCHLEUNIGUNG
# ===========================================================================

# Anzahl der Terme für die Euler-Knopp-Beschleunigung der alternierenden Eta-Reihe
# η(s) = Σ (-1)^{k+1} / k^s, Euler-Knopp transformiert auf schnell konvergente Folge.
# 60 Terme liefern Maschinengenauigkeit (~16 Stellen) für Re(s) > 0
ZETA_EULER_KNOPP_TERMS: int = 60

# Anzahl der Terme für die Euler-Maclaurin-Summation der Zeta-Funktion
# Wird für Re(s) > 1 verwendet, wo die Reihe direkt konvergiert.
# 200 Terme liefern ausreichend genaue Werte für Standardanwendungen
ZETA_EULER_MACLAURIN_TERMS: int = 200

# ===========================================================================
# LANCZOS-APPROXIMATION DER GAMMA-FUNKTION
# ===========================================================================

# Lanczos-Parameter g = 7 (Verschiebung der Pol-Annäherung)
# Koeffizienten p_k für die Lanczos-Approximation:
#   Γ(z+1) ≈ √(2π) · [(z+g+0.5)/(e)]^(z+0.5) · A_g(z)
#   wobei A_g(z) = c_0 + Σ_{k=1}^{g} c_k / (z+k)
# Koeffizienten aus: Lanczos (1964), "A precision approximation of the gamma function"
LANCZOS_COEFFICIENTS: dict[str, int | list[float]] = {
    'g': 7,
    'coefficients': [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ]
}

# ===========================================================================
# MPMATH-KONFIGURATION (Hochpräzisions-Arithmetik)
# ===========================================================================

# Standard-Anzahl der Dezimalstellen für mpmath-Berechnungen
# 50 Dezimalstellen entsprechen ~166 Bit Genauigkeit
# Für Riemann-Hypothese-Berechnungen empfehlenswert: 100+
MPMATH_DEFAULT_DPS: int = 50

# ===========================================================================
# MATHEMATISCHE KONSTANTEN
# ===========================================================================

# Goldener Schnitt: φ = (1 + √5) / 2 ≈ 1.6180339887...
# Tritt in Fibonacci-Zahlen, Spiralen und Optimierungsalgorithmen auf
# Golden Section Search: Intervall wird im Verhältnis φ geteilt
GOLDEN_RATIO: float = (1.0 + math.sqrt(5.0)) / 2.0  # ≈ 1.618033988749895

# Reziproker Goldener Schnitt: 1/φ = φ - 1 = 2/(1+√5) ≈ 0.6180339887...
# Wird in Golden-Section-Search als Schrittfaktor verwendet:
# x1 = b - INV_GOLDEN_RATIO * (b-a); x2 = a + INV_GOLDEN_RATIO * (b-a)
INV_GOLDEN_RATIO: float = 2.0 / (1.0 + math.sqrt(5.0))  # ≈ 0.6180339887498949

# ===========================================================================
# PROJEKT-METADATEN
# ===========================================================================

# Autor aller Quelldateien im Projekt
AUTHOR: str = "Kurt Ingwer"

# Aktuelle Versionsnummer des Projekts (entspricht Build-Nummer)
VERSION: str = "11.0"


# ===========================================================================
# HILFSFUNKTIONEN – SYMBOLISCHE VEREINFACHUNG
# ===========================================================================

def smart_simplify(expr) -> "sp.Expr":
    """
    @brief Vereinfacht einen SymPy-Ausdruck intelligent mit mehreren Strategien.
    @description
        Wendet nacheinander sp.simplify() und sp.nsimplify() an, um den Ausdruck
        in seine einfachste Form zu bringen:

        1. sp.simplify(expr): Allgemeine algebraische und trigonometrische Vereinfachung.
           Beispiel: sin²(x) + cos²(x) → 1
        2. sp.nsimplify(expr, rational=True): Erkennt einfache rationale Zahlen und
           bekannte Konstanten (π, √2, e) aus Gleitkomma-Näherungen.
           Beispiel: 3.14159... → π

        Bei Fehler (z.B. nicht-symbolischer Ausdruck) wird der Ausdruck unverändert zurückgegeben.

        Anwendung: Alle Funktionen, die SymPy-Objekte zurückgeben, können das Ergebnis
        durch smart_simplify() lesbarer und kompakter machen.

    @param expr: SymPy-Ausdruck oder konvertierbarer Wert.
    @return: Vereinfachter SymPy-Ausdruck (oder unveränderter Eingabe bei Fehler).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Beispiele:
        >>> import sympy as sp
        >>> x = sp.Symbol('x')
        >>> smart_simplify(sp.sin(x)**2 + sp.cos(x)**2)
        1
        >>> smart_simplify(sp.Rational(6, 4))
        3/2
        >>> smart_simplify(3.141592653589793)
        pi
    """
    try:
        # Schritt 1: Allgemeine Vereinfachung (Algebra, Trigonometrie, Logarithmen)
        simplified = sp.simplify(expr)

        # Schritt 2: Rationale Vereinfachung (erkennt π, √n, e aus Floats)
        # rational=True: versucht Brüche statt Dezimalzahlen darzustellen
        simplified = sp.nsimplify(simplified, rational=True)

        return simplified
    except Exception:
        # Bei Fehler: unveränderter Ausdruck zurückgeben (z.B. für Nicht-SymPy-Typen)
        return expr


# ===========================================================================
# HEURISTIK-SCHWELLWERTE FÜR AUTO_COMPUTE
# ===========================================================================

# Maximaler Polynomgrad für symbolische Auswertung.
# Oberhalb dieses Werts wird auf numerische Berechnung umgeschaltet.
AUTO_COMPUTE_MAX_DEGREE: int = 10

# Maximale Anzahl symbolischer Operationen (SymPy count_ops) für symbolische Auswertung.
# Komplexere Ausdrücke (viele Operationen) werden numerisch ausgewertet.
AUTO_COMPUTE_MAX_OPS: int = 50


def auto_compute(expr, x_val: float) -> Union[float, complex]:
    """
    Wählt automatisch zwischen symbolischer (SymPy) und numerischer (numpy) Auswertung.

    Heuristik:
        - Wenige Variablen UND niedriger Grad (≤ AUTO_COMPUTE_MAX_DEGREE)
          UND wenige Operationen (≤ AUTO_COMPUTE_MAX_OPS) → SymPy symbolisch
        - Sonst → numpy numerisch (schneller für hochdimensionale/komplexe Ausdrücke)

    Warum diese Unterscheidung?
        - SymPy ist exakt, aber langsam für große Ausdrücke (symbolische Vereinfachung)
        - numpy.float64 ist ~1000× schneller, aber akkumuliert Rundungsfehler
        - Für Polynome kleinen Grades lohnt sich der Exaktheits-Vorteil von SymPy

    Beispiel:
        >>> import sympy as sp
        >>> x = sp.Symbol('x')
        >>> auto_compute(x**2 + 3*x + 1, 2.0)   # symbolisch (Grad 2 ≤ 10)
        (11+0j)
        >>> auto_compute(x**15, 2.0)              # numerisch (Grad 15 > 10)
        32768.0

    @param expr: SymPy-Ausdruck oder Callable (f(x_val))
    @param x_val: Auswertungspunkt (reelle Zahl)
    @return: Numerischer Wert als float oder complex
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    import numpy as np

    # Prüfe ob expr ein SymPy-Ausdruck ist
    if not isinstance(expr, sp.Expr):
        # Callable: direkt numerisch auswerten
        return float(expr(x_val))

    # Variablenmenge bestimmen
    free_vars = expr.free_symbols

    if len(free_vars) == 0:
        # Konstanter Ausdruck: symbolisch auswerten, kein x_val benötigt
        return complex(expr.evalf())

    if len(free_vars) > 1:
        # Mehrere Variablen: nur erste Variable mit x_val belegen, numerisch auswerten
        x_sym = list(free_vars)[0]
        return float(expr.subs(x_sym, x_val).evalf())

    # Genau eine Variable: Grad und Operationszahl prüfen
    x_sym = list(free_vars)[0]

    # Grad des Ausdrucks bestimmen (nur sinnvoll für Polynome)
    try:
        deg = sp.degree(expr, x_sym)
        if deg == sp.oo:
            # Nicht-Polynome (sin, exp, …) als "hohen Grad" behandeln
            deg = AUTO_COMPUTE_MAX_DEGREE + 1
    except Exception:
        deg = AUTO_COMPUTE_MAX_DEGREE + 1  # Unbekannter Grad → numerisch

    # Operationsanzahl messen (Komplexitätsmaß für SymPy-Ausdruck)
    ops = sp.count_ops(expr)

    # Entscheidung: symbolisch oder numerisch?
    if deg <= AUTO_COMPUTE_MAX_DEGREE and ops <= AUTO_COMPUTE_MAX_OPS:
        # Symbolisch: exakte Auswertung via SymPy-Substitution
        result = expr.subs(x_sym, x_val).evalf()
        return complex(result)
    else:
        # Numerisch: schnelle Auswertung via numpy lambdify
        # lambdify erzeugt eine reine Python/numpy-Funktion aus dem SymPy-Ausdruck
        f_num = sp.lambdify(x_sym, expr, modules=['numpy'])
        return float(f_num(x_val))
