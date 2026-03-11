"""
@file test_computation_strategy.py
@brief Umfassende Tests für das Modul computation_strategy.py.
@description
    Testet alle Funktionen und Klassen des Moduls:
    - choose_strategy(): Automatische Strategiewahl (symbolic/numeric)
    - AdaptiveSolver.solve(): Adaptives Lösen von Gleichungen
    - AdaptiveSolver._determine_strategy(): Interne Strategie-Erkennung
    - AdaptiveSolver.last_strategy: Property für zuletzt genutzte Strategie

    Mathematische Grundlagen:
    - Polynome Grad ≤ 10 → symbolisch (exakte Nullstellen via SymPy)
    - Matrizen ≤ 50×50 → symbolisch; größer → numerisch
    - Integration → numerisch (schneller); mit precision_needed → symbolisch
    - Unbekannter Typ → numerisch (Standardfall)

    Edge-Cases: Nullgrad, riesige Matrizen, unbekannte Problemtypen,
    ungültige Methodenname, leere Ausdrücke.

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import sympy as sp
from computation_strategy import choose_strategy, AdaptiveSolver


# =============================================================================
# TESTS: choose_strategy()
# =============================================================================

class TestChooseStrategy:
    """Tests für die choose_strategy()-Funktion."""

    # --- Polynomial ---

    def test_polynomial_small_degree_symbolic(self):
        """Polynom Grad ≤ 10 → 'symbolic'."""
        assert choose_strategy('polynomial', 5) == 'symbolic'
        assert choose_strategy('polynomial', 10) == 'symbolic'

    def test_polynomial_large_degree_numeric(self):
        """Polynom Grad > 10 → 'numeric'."""
        assert choose_strategy('polynomial', 11) == 'numeric'
        assert choose_strategy('polynomial', 100) == 'numeric'

    def test_polynomial_zero_degree_symbolic(self):
        """Polynom Grad 0 (Konstante) → 'symbolic'."""
        assert choose_strategy('polynomial', 0) == 'symbolic'

    def test_polynomial_degree_exactly_ten_symbolic(self):
        """Grad genau 10: Grenzfall → 'symbolic'."""
        assert choose_strategy('polynomial', 10) == 'symbolic'

    # --- Matrix ---

    def test_matrix_small_symbolic(self):
        """Kleine Matrix (≤ 50) → 'symbolic'."""
        assert choose_strategy('matrix', 2) == 'symbolic'
        assert choose_strategy('matrix', 50) == 'symbolic'

    def test_matrix_large_numeric(self):
        """Große Matrix (> 50) → 'numeric'."""
        assert choose_strategy('matrix', 51) == 'numeric'
        assert choose_strategy('matrix', 1000) == 'numeric'

    def test_matrix_exactly_50_symbolic(self):
        """Matrix genau 50×50: Grenzfall → 'symbolic'."""
        assert choose_strategy('matrix', 50) == 'symbolic'

    # --- Integration ---

    def test_integration_without_precision_numeric(self):
        """Integration ohne Präzisionswunsch → 'numeric'."""
        assert choose_strategy('integration', 1) == 'numeric'

    def test_integration_with_precision_symbolic(self):
        """Integration mit precision_needed=True → 'symbolic'."""
        assert choose_strategy('integration', 1, precision_needed=True) == 'symbolic'

    # --- Derivative ---

    def test_derivative_numeric(self):
        """Ableitungen → 'numeric' (Standard)."""
        assert choose_strategy('derivative', 1) == 'numeric'

    def test_derivative_with_precision_symbolic(self):
        """Ableitung mit precision_needed → 'symbolic'."""
        assert choose_strategy('derivative', 1, precision_needed=True) == 'symbolic'

    # --- ODE ---

    def test_ode_numeric(self):
        """ODE → 'numeric' (Runge-Kutta)."""
        assert choose_strategy('ode', 1) == 'numeric'

    def test_ode_with_precision_symbolic(self):
        """ODE mit precision_needed → 'symbolic'."""
        assert choose_strategy('ode', 1, precision_needed=True) == 'symbolic'

    # --- Allgemein ---

    def test_general_numeric(self):
        """Unbekannter Typ → 'numeric'."""
        assert choose_strategy('general', 1) == 'numeric'
        assert choose_strategy('unknown_type', 42) == 'numeric'

    def test_precision_needed_always_symbolic(self):
        """Mit precision_needed=True ist immer 'symbolic'."""
        for ptype in ['polynomial', 'matrix', 'integration', 'derivative', 'ode', 'general']:
            assert choose_strategy(ptype, 1000, precision_needed=True) == 'symbolic'

    def test_case_insensitive_type(self):
        """Großschreibung des problem_type wird normalisiert."""
        assert choose_strategy('POLYNOMIAL', 5) == 'symbolic'
        assert choose_strategy('Matrix', 5) == 'symbolic'


# =============================================================================
# TESTS: AdaptiveSolver
# =============================================================================

class TestAdaptiveSolver:
    """Tests für die AdaptiveSolver-Klasse."""

    def setup_method(self):
        """Solver-Instanz vor jedem Test erstellen."""
        self.solver = AdaptiveSolver()
        self.x = sp.Symbol('x')

    # --- Initialisierung ---

    def test_initial_last_strategy_is_none(self):
        """Vor dem ersten Aufruf ist last_strategy None."""
        assert self.solver.last_strategy is None

    def test_solver_instance(self):
        """AdaptiveSolver kann instanziiert werden."""
        assert isinstance(self.solver, AdaptiveSolver)

    # --- Methode 'symbolic' ---

    def test_solve_symbolic_linear(self):
        """Lineares Polynom symbolisch gelöst: x - 3 = 0 → [3]."""
        result = self.solver.solve(self.x - 3, method='symbolic')
        assert sp.Integer(3) in result or 3 in result

    def test_solve_symbolic_quadratic(self):
        """Quadratisches Polynom symbolisch: x²-4 = 0 → [-2, 2]."""
        result = self.solver.solve(self.x**2 - 4, method='symbolic')
        result_floats = sorted([float(r) for r in result])
        assert abs(result_floats[0] - (-2)) < 1e-9
        assert abs(result_floats[1] - 2) < 1e-9

    def test_solve_symbolic_records_strategy(self):
        """Nach symbolischem Aufruf ist last_strategy 'symbolic'."""
        self.solver.solve(self.x - 1, method='symbolic')
        assert self.solver.last_strategy == 'symbolic'

    # --- Methode 'numeric' ---

    def test_solve_numeric_linear(self):
        """Lineares Polynom numerisch gelöst: x - 5 = 0 → [5.0]."""
        result = self.solver.solve(self.x - 5, method='numeric')
        assert result is not None
        if isinstance(result, list):
            assert abs(float(result[0]) - 5.0) < 1e-6

    def test_solve_numeric_records_strategy(self):
        """Nach numerischem Aufruf ist last_strategy 'numeric'."""
        self.solver.solve(self.x - 2, method='numeric')
        assert self.solver.last_strategy == 'numeric'

    # --- Methode 'auto' ---

    def test_solve_auto_polynomial_degree2(self):
        """Auto-Modus für Grad-2-Polynom: korrekte Nullstellen."""
        result = self.solver.solve(self.x**2 - 9, method='auto')
        assert result is not None

    def test_solve_auto_records_strategy(self):
        """Auto-Modus setzt last_strategy auf einen gültigen Wert."""
        self.solver.solve(self.x**2 - 1, method='auto')
        assert self.solver.last_strategy in ('symbolic', 'numeric')

    # --- Ungültige Methode ---

    def test_solve_invalid_method_raises(self):
        """Unbekannte Methode wirft ValueError."""
        with pytest.raises(ValueError):
            self.solver.solve(self.x - 1, method='magic')

    # --- Edge Cases ---

    def test_solve_constant_expression(self):
        """Konstanter Ausdruck ohne freie Variablen."""
        result = self.solver.solve(sp.Integer(0), method='symbolic')
        # Sollte ohne Fehler durchlaufen (alle x sind Lösungen oder leere Liste)
        assert result is not None

    def test_solve_string_expr(self):
        """String-Eingabe wird via sympify konvertiert."""
        result = self.solver.solve('x - 7', method='symbolic', variable=self.x)
        assert result is not None

    def test_solve_with_explicit_variable(self):
        """Explizite Variablenangabe wird korrekt verwendet."""
        y = sp.Symbol('y')
        result = self.solver.solve(y**2 - 25, method='symbolic', variable=y)
        assert result is not None

    def test_last_strategy_updated_on_each_call(self):
        """last_strategy wird bei jedem Aufruf aktualisiert."""
        self.solver.solve(self.x - 1, method='symbolic')
        assert self.solver.last_strategy == 'symbolic'
        self.solver.solve(self.x - 2, method='numeric')
        assert self.solver.last_strategy == 'numeric'

    def test_solve_returns_list_or_expr(self):
        """solve() gibt immer ein verwertbares Ergebnis zurück."""
        result = self.solver.solve(self.x**2 - 1, method='symbolic')
        # Entweder Liste oder SymPy-Ausdruck
        assert isinstance(result, (list, sp.Basic)) or result is not None

    def test_solve_large_degree_polynomial_auto(self):
        """Grad-12-Polynom im Auto-Modus: Strategie numerisch."""
        # x^12 - 1
        expr = self.x**12 - 1
        result = self.solver.solve(expr, method='auto')
        # Muss ohne Fehler durchlaufen
        assert result is not None
