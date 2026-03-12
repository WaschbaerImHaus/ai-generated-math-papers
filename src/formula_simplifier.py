"""
@file formula_simplifier.py
@brief Automatische Formel-Vereinfachung mit SymPy.
@description
    Dieses Modul stellt Werkzeuge zur automatischen Vereinfachung
    mathematischer Ausdrücke bereit. Es nutzt verschiedene SymPy-
    Vereinfachungsstrategien und wählt die kompakteste Darstellung.

    Strategien:
    - simplify():   Allgemeine algebraische Vereinfachung
    - nsimplify():  Annäherung an rationale/irrationale Zahlen
    - radsimp():    Vereinfachung durch Rationalisierung des Nenners
    - trigsimp():   Trigonometrische Vereinfachungen
    - expand():     Ausmultiplizieren von Ausdrücken

    Die Methode simplify_auto() wählt automatisch die kürzeste Darstellung.

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import re as _re_fs
import sympy as sp
from typing import Any
from sympy.parsing.sympy_parser import (
    parse_expr as _parse_expr_fs,
    standard_transformations as _std_tr_fs,
    implicit_multiplication_application as _impl_mult_fs,
)

# Erlaubte Transformationen für sicheres Parsing
_SAFE_TR_FS = _std_tr_fs + (_impl_mult_fs,)

# Whitelist-Muster: Nur mathematische Ausdrücke (Buchstaben, Ziffern, Operatoren)
_SAFE_CHARS_FS = _re_fs.compile(r'^[a-zA-Z0-9\s\+\-\*\/\^\(\)\.\,\_\=\<\>]+$')


def _safe_to_sympy(expr: Any) -> sp.Basic:
    """
    @brief Konvertiert Eingabe sicher in einen SymPy-Ausdruck.
    @description
        Wenn expr bereits ein SymPy-Ausdruck ist, wird er unverändert zurückgegeben.
        Strings werden mit parse_expr geparst (sicherer als sympify).
        Andere Typen (int, float) werden direkt konvertiert.
    @param expr: SymPy-Ausdruck, String oder numerischer Wert
    @return SymPy-Ausdruck
    @raises ValueError bei ungültigem String-Ausdruck
    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """
    if isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, str):
        if len(expr) > 1000 or not _SAFE_CHARS_FS.match(expr):
            raise ValueError(f"Ungültiger Ausdruck (unerlaubte Zeichen oder zu lang): '{expr}'")
        return _parse_expr_fs(expr, transformations=_SAFE_TR_FS)
    # int, float und ähnliche Typen direkt via sp.sympify konvertieren (kein String → sicher)
    return sp.sympify(expr)


class FormulaSimplifier:
    """
    @brief Klasse zur automatischen Vereinfachung mathematischer Ausdrücke.
    @description
        Bietet mehrere Vereinfachungsstrategien und wählt automatisch
        die am besten geeignete aus. Unterstützt auch LaTeX-Export
        und Batch-Verarbeitung von Ergebnis-Dictionaries.

    @example
        fs = FormulaSimplifier()
        x = sp.Symbol('x')
        result = fs.simplify_auto(sp.sin(x)**2 + sp.cos(x)**2)
        # Ergibt: 1

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def simplify(self, expr: Any) -> sp.Expr:
        """
        @brief Vereinfacht einen SymPy-Ausdruck mit mehreren Strategien.
        @description
            Wendet der Reihe nach simplify(), nsimplify() und radsimp() an.
            Gibt die erste wesentliche Vereinfachung zurück.

            Mathematische Vereinfachungsregeln:
            - sin²(x) + cos²(x) = 1 (Pythagoreischer Satz)
            - e^(ln x) = x
            - a/b + c/b = (a+c)/b (Bruchaddition)

        @param expr: SymPy-Ausdruck oder konvertierbarer Wert.
        @return: Vereinfachter SymPy-Ausdruck.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Eingabe sicher zu SymPy-Ausdruck konvertieren (parse_expr statt sympify)
        expr = _safe_to_sympy(expr)

        # Schrittweise Vereinfachungen anwenden
        simplified = sp.simplify(expr)

        # Versuche nsimplify für numerische Annäherung (z.B. 3.14159... → pi)
        try:
            ns = sp.nsimplify(simplified, rational=False, tolerance=1e-8)
            if len(str(ns)) < len(str(simplified)):
                simplified = ns
        except Exception:
            pass  # nsimplify kann bei manchen Ausdrücken scheitern

        # Versuche radsimp (Rationalisierung des Nenners)
        try:
            rs = sp.radsimp(simplified)
            if len(str(rs)) < len(str(simplified)):
                simplified = rs
        except Exception:
            pass

        return simplified

    def simplify_auto(self, expr: Any) -> sp.Expr:
        """
        @brief Wählt automatisch die kürzeste Vereinfachungsform.
        @description
            Probiert verschiedene Vereinfachungsstrategien aus und gibt
            die Form mit der kürzesten String-Darstellung zurück.

            Verglichene Strategien:
            1. Original-Ausdruck
            2. sp.simplify() – allgemeine Vereinfachung
            3. sp.nsimplify() – numerische Annäherung
            4. sp.trigsimp() – trigonometrische Identitäten
            5. sp.powsimp() – Potenz-Vereinfachungen
            6. sp.cancel() – Bruchkürzung
            7. sp.factor() – Faktorisierung
            8. sp.expand() – Ausmultiplizieren

        @param expr: SymPy-Ausdruck oder konvertierbarer Wert.
        @return: Kürzeste Darstellung des Ausdrucks.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Eingabe sicher zu SymPy konvertieren (parse_expr statt sympify)
        expr = _safe_to_sympy(expr)

        # Kandidaten-Liste aller Vereinfachungsformen
        candidates = [expr]

        # Alle Strategien versuchen und Kandidaten sammeln
        strategies = [
            sp.simplify,
            lambda e: sp.nsimplify(e, rational=False, tolerance=1e-8),
            sp.trigsimp,
            sp.powsimp,
            sp.cancel,
            sp.factor,
            sp.expand,
            sp.radsimp,
        ]

        for strategy in strategies:
            try:
                result = strategy(expr)
                if result is not None:
                    candidates.append(result)
            except Exception:
                pass  # Strategie kann bei bestimmten Ausdrücken versagen

        # Kürzeste String-Darstellung auswählen
        # Bei Gleichheit wird die einfachere Form bevorzugt (erste in Liste)
        best = min(candidates, key=lambda e: len(str(e)))
        return best

    def to_latex(self, expr: Any) -> str:
        """
        @brief Konvertiert einen SymPy-Ausdruck in einen LaTeX-String.
        @description
            Nutzt sp.latex() für die Konvertierung. Der erzeugte LaTeX-Code
            kann direkt in KaTeX oder MathJax verwendet werden.

            Beispiele:
            - x**2 → x^{2}
            - sqrt(x) → \\sqrt{x}
            - 1/3 → \\frac{1}{3}
            - Integral → \\int

        @param expr: SymPy-Ausdruck oder konvertierbarer Wert.
        @return: LaTeX-String (ohne Umgebungs-Delimitoren).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Eingabe sicher zu SymPy konvertieren (parse_expr statt sympify)
        expr = _safe_to_sympy(expr)

        # LaTeX erzeugen mit SymPy
        return sp.latex(expr)

    def simplify_all(self, results_dict: dict) -> dict:
        """
        @brief Vereinfacht alle Werte in einem Ergebnis-Dictionary.
        @description
            Iteriert über alle Schlüssel-Wert-Paare und versucht,
            SymPy-kompatible Werte zu vereinfachen. Nicht-SymPy-Werte
            (z.B. int, float, str) werden unverändert übernommen.

            Typischer Anwendungsfall: Ergebnis-Dictionaries aus
            Gleichungslösern, Integratoren oder Matrix-Berechnungen.

        @param results_dict: Dictionary mit mathematischen Ergebnissen.
        @return: Neues Dictionary mit vereinfachten Werten.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        simplified_dict = {}

        for key, value in results_dict.items():
            try:
                # Nur SymPy-kompatible Werte vereinfachen
                if isinstance(value, (sp.Basic, int, float)):
                    simplified_dict[key] = self.simplify_auto(value)
                elif isinstance(value, (list, tuple)):
                    # Listen/Tupel element-weise vereinfachen
                    simplified_list = []
                    for item in value:
                        try:
                            simplified_list.append(self.simplify_auto(item))
                        except Exception:
                            simplified_list.append(item)
                    simplified_dict[key] = type(value)(simplified_list)
                else:
                    # Sonstige Typen (str, bool, etc.) unverändert übernehmen
                    simplified_dict[key] = value
            except Exception:
                # Fallback: Original-Wert bei Fehler
                simplified_dict[key] = value

        return simplified_dict
