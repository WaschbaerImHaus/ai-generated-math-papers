"""
@file computation_strategy.py
@brief Automatische Auswahl zwischen symbolischer und numerischer Berechnung.
@description
    Dieses Modul entscheidet automatisch, ob ein mathematisches Problem
    symbolisch (exakt, mit SymPy) oder numerisch (schnell, mit NumPy/SciPy)
    gelöst werden soll.

    Entscheidungsregeln:
    - Polynome mit Grad ≤ 10 → symbolisch
    - Matrizen > 50×50 → numerisch
    - Integration mit exaktem Ergebnis gewünscht → symbolisch
    - Große Datenmengen → numerisch
    - Standardfall → symbolisch für kleine, numerisch für große Probleme

    Die Klasse AdaptiveSolver wählt automatisch die Strategie und
    löst das Problem entsprechend.

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import sympy as sp
import numpy as np
from typing import Any, Optional


# Schwellenwerte für die Strategiewahl
# Polynome bis Grad 10 können noch effizient symbolisch gelöst werden
_MAX_POLYNOMIAL_DEGREE_SYMBOLIC = 10

# Matrizen über 50×50 sind numerisch deutlich effizienter
_MAX_MATRIX_SIZE_SYMBOLIC = 50

# Maximale Anzahl von Integrationspunkten für symbolische Berechnung
_MAX_INTEGRATION_POINTS_SYMBOLIC = 100


def choose_strategy(
    problem_type: str,
    size: Any,
    precision_needed: bool = False
) -> str:
    """
    @brief Wählt automatisch die optimale Berechnungsstrategie.
    @description
        Entscheidet basierend auf Problemtyp und Größe, ob symbolisch
        oder numerisch gerechnet werden soll.

        Entscheidungsbaum:
        1. Wenn precision_needed=True → symbolisch (exakte Ergebnisse gewünscht)
        2. Polynome: Grad ≤ 10 → symbolisch, sonst → numerisch
        3. Matrizen: ≤ 50×50 → symbolisch, sonst → numerisch
        4. Integration: Mit precision_needed → symbolisch, sonst → numerisch
        5. Standardfall: numerisch für Geschwindigkeit

    @param problem_type: Art des Problems: 'polynomial', 'matrix', 'integration',
                         'derivative', 'ode', 'general'
    @param size: Größe des Problems (Grad für Polynome, n für n×n-Matrizen, usw.)
    @param precision_needed: Wenn True, wird symbolische Berechnung bevorzugt.
    @return: 'symbolic' oder 'numeric'
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Exaktheit explizit gewünscht → immer symbolisch
    if precision_needed:
        return 'symbolic'

    problem_type = problem_type.lower()

    # Polynome: kleiner Grad → symbolisch (SymPy kann exakte Nullstellen berechnen)
    if problem_type == 'polynomial':
        degree = int(size) if size is not None else 0
        if degree <= _MAX_POLYNOMIAL_DEGREE_SYMBOLIC:
            return 'symbolic'
        else:
            return 'numeric'

    # Matrizen: kleine Matrizen → symbolisch (exakte Determinante/Eigenwerte)
    elif problem_type == 'matrix':
        matrix_size = int(size) if size is not None else 0
        if matrix_size <= _MAX_MATRIX_SIZE_SYMBOLIC:
            return 'symbolic'
        else:
            return 'numeric'

    # Integration: ohne Präzisionsanforderung → numerisch (schneller)
    elif problem_type == 'integration':
        return 'numeric'

    # Ableitung: üblicherweise numerisch ausreichend
    elif problem_type == 'derivative':
        return 'numeric'

    # ODE: numerisch (Runge-Kutta etc. sind standard)
    elif problem_type == 'ode':
        return 'numeric'

    # Standardfall: numerisch für Geschwindigkeit
    else:
        return 'numeric'


class AdaptiveSolver:
    """
    @brief Adaptiver Löser der automatisch zwischen Strategien wählt.
    @description
        Wählt je nach Ausdruck und Problemgröße zwischen symbolischer
        und numerischer Berechnung. Nutzt choose_strategy() für die
        automatische Entscheidung.

    @example
        solver = AdaptiveSolver()
        x = sp.Symbol('x')
        result = solver.solve(x**2 - 4, method='auto')
        # Ergibt: [-2, 2] (symbolisch, da Grad 2 ≤ 10)

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self):
        """
        @brief Initialisiert den AdaptiveSolver.
        @description
            Setzt Standardkonfiguration für Schwellenwerte und
            bevorzugte Strategien.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Standardkonfiguration
        self._prefer_symbolic = False   # Im Zweifel numerisch für Geschwindigkeit
        self._last_strategy = None      # Merkt die zuletzt verwendete Strategie

    def solve(self, expr: Any, method: str = 'auto', variable: Optional[sp.Symbol] = None) -> Any:
        """
        @brief Löst einen mathematischen Ausdruck mit automatischer Strategiewahl.
        @description
            Bei method='auto': Strategie wird automatisch gewählt.
            Bei method='symbolic': Erzwingt symbolische Berechnung (SymPy).
            Bei method='numeric': Erzwingt numerische Berechnung (NumPy).

            Für Polynome wird nach Nullstellen gesucht.
            Für andere Ausdrücke wird sp.simplify() oder np.float64() genutzt.

        @param expr: SymPy-Ausdruck oder Python-Funktion.
        @param method: 'auto', 'symbolic' oder 'numeric'
        @param variable: SymPy-Symbol für die Auflösung (optional, Standard: x)
        @return: Ergebnis der Berechnung.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Strategie bestimmen
        if method == 'auto':
            # Automatische Strategiewahl
            strategy = self._determine_strategy(expr)
        elif method in ('symbolic', 'numeric'):
            strategy = method
        else:
            raise ValueError(f"Unbekannte Methode '{method}'. "
                             f"Verwende 'auto', 'symbolic' oder 'numeric'.")

        # Strategie merken für Abfrage
        self._last_strategy = strategy

        # Eingabe zu SymPy konvertieren falls nötig
        if not isinstance(expr, sp.Basic):
            try:
                expr = sp.sympify(expr)
            except Exception:
                # Fallback: direkte Rückgabe wenn nicht konvertierbar
                return expr

        # Variable bestimmen (Standard: x oder erste freie Variable)
        if variable is None:
            free_vars = list(expr.free_symbols)
            if free_vars:
                variable = sorted(free_vars, key=str)[0]
            else:
                variable = sp.Symbol('x')

        # Je nach Strategie lösen
        if strategy == 'symbolic':
            return self._solve_symbolic(expr, variable)
        else:
            return self._solve_numeric(expr, variable)

    def _determine_strategy(self, expr: Any) -> str:
        """
        @brief Bestimmt automatisch die optimale Strategie für einen Ausdruck.
        @description
            Analysiert den Ausdruck und wählt basierend auf Typ und Größe.

        @param expr: Zu analysierender Ausdruck.
        @return: 'symbolic' oder 'numeric'
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # SymPy-Ausdruck analysieren
        if isinstance(expr, sp.Basic):
            # Prüfe ob es ein Polynom ist
            free_vars = list(expr.free_symbols)
            if free_vars:
                var = sorted(free_vars, key=str)[0]
                try:
                    poly = sp.Poly(expr, var)
                    degree = poly.degree()
                    return choose_strategy('polynomial', degree)
                except (sp.PolynomialError, Exception):
                    pass

        return 'symbolic'  # Standardfall für SymPy-Ausdrücke

    def _solve_symbolic(self, expr: sp.Expr, variable: sp.Symbol) -> Any:
        """
        @brief Löst symbolisch mit SymPy.
        @param expr: SymPy-Ausdruck.
        @param variable: Auflösungsvariable.
        @return: Lösungsliste oder vereinfachter Ausdruck.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        try:
            # Versuche Gleichung zu lösen (Ausdruck = 0)
            solutions = sp.solve(expr, variable)
            if solutions is not None:
                return solutions
        except Exception:
            pass

        # Fallback: Vereinfachen
        return sp.simplify(expr)

    def _solve_numeric(self, expr: sp.Expr, variable: sp.Symbol) -> Any:
        """
        @brief Löst numerisch mit NumPy.
        @description
            Konvertiert den Ausdruck zu einer Python-Funktion und
            nutzt NumPy für numerische Berechnungen.

        @param expr: SymPy-Ausdruck.
        @param variable: Auflösungsvariable.
        @return: Numerische Näherungslösung.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        try:
            # Symbolisch auswerten für numerische Annäherung
            solutions = sp.nsolve(expr, variable, 0)
            return [float(solutions)]
        except Exception:
            pass

        # Fallback: symbolisch
        return sp.solve(expr, variable)

    @property
    def last_strategy(self) -> Optional[str]:
        """
        @brief Gibt die zuletzt verwendete Strategie zurück.
        @return: 'symbolic', 'numeric' oder None wenn noch nicht aufgerufen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return self._last_strategy
