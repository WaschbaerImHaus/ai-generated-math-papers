"""
@file test_formula_simplifier.py
@brief Tests für das FormulaSimplifier-Modul.
@description
    Testet alle Methoden der FormulaSimplifier-Klasse:
    - simplify(): Allgemeine Vereinfachung
    - simplify_auto(): Automatische Auswahl der kürzesten Form
    - to_latex(): LaTeX-Konvertierung
    - simplify_all(): Batch-Vereinfachung

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import sympy as sp
import pytest
from formula_simplifier import FormulaSimplifier


class TestFormulaSimplifier:
    """Tests für die FormulaSimplifier-Klasse."""

    def setup_method(self):
        """Erstellt eine FormulaSimplifier-Instanz vor jedem Test."""
        self.fs = FormulaSimplifier()
        # Symbolische Variablen
        self.x = sp.Symbol('x')
        self.y = sp.Symbol('y')
        self.n = sp.Symbol('n', positive=True)

    def test_simplify_trig_identity(self):
        """Test: sin²(x) + cos²(x) wird zu 1 vereinfacht."""
        expr = sp.sin(self.x)**2 + sp.cos(self.x)**2
        result = self.fs.simplify(expr)
        assert result == sp.Integer(1), f"Erwartet 1, bekam {result}"

    def test_simplify_fraction(self):
        """Test: (x²-1)/(x-1) wird zu (x+1) vereinfacht."""
        expr = (self.x**2 - 1) / (self.x - 1)
        result = self.fs.simplify(expr)
        # Beide Formen sind korrekt: x+1 oder Äquivalentes
        assert sp.simplify(result - (self.x + 1)) == 0

    def test_simplify_exp_log(self):
        """Test: exp(log(x)) wird zu x vereinfacht (für positive x)."""
        x_pos = sp.Symbol('x', positive=True)
        expr = sp.exp(sp.log(x_pos))
        result = self.fs.simplify(expr)
        assert result == x_pos

    def test_simplify_integer_input(self):
        """Test: Ganzzahl-Eingabe wird korrekt verarbeitet."""
        result = self.fs.simplify(42)
        assert result == sp.Integer(42)

    def test_simplify_auto_chooses_shorter(self):
        """Test: simplify_auto() gibt kürzere Darstellung zurück."""
        expr = sp.sin(self.x)**2 + sp.cos(self.x)**2
        result = self.fs.simplify_auto(expr)
        # Das Ergebnis soll kürzer als der Originalausdruck sein
        assert len(str(result)) <= len(str(expr))

    def test_simplify_auto_constant(self):
        """Test: simplify_auto() vereinfacht Konstante korrekt."""
        result = self.fs.simplify_auto(sp.pi * 0)
        assert result == sp.Integer(0)

    def test_simplify_auto_power(self):
        """Test: simplify_auto() mit Potenzausdrücken."""
        expr = self.x**2 * self.x**3
        result = self.fs.simplify_auto(expr)
        # x^5 oder äquivalent
        assert sp.simplify(result - self.x**5) == 0

    def test_to_latex_simple(self):
        """Test: Einfacher Ausdruck → LaTeX."""
        expr = self.x**2 + 1
        latex = self.fs.to_latex(expr)
        # Muss x^2 enthalten
        assert 'x' in latex
        assert '2' in latex

    def test_to_latex_fraction(self):
        """Test: Bruch wird als \\frac{...}{...} dargestellt."""
        expr = sp.Rational(1, 3)
        latex = self.fs.to_latex(expr)
        assert 'frac' in latex

    def test_to_latex_sqrt(self):
        """Test: Wurzel wird als \\sqrt{...} dargestellt."""
        expr = sp.sqrt(self.x)
        latex = self.fs.to_latex(expr)
        assert 'sqrt' in latex

    def test_to_latex_string_input(self):
        """Test: String-Eingabe wird zu LaTeX konvertiert."""
        latex = self.fs.to_latex('x**2 + 1')
        assert 'x' in latex

    def test_simplify_all_dict(self):
        """Test: simplify_all() vereinfacht alle Werte im Dictionary."""
        expr1 = sp.sin(self.x)**2 + sp.cos(self.x)**2
        expr2 = self.x**2 * self.x
        data = {'a': expr1, 'b': expr2, 'c': 'text'}

        result = self.fs.simplify_all(data)

        # Trigonometrischer Ausdruck sollte vereinfacht werden
        assert sp.simplify(result['a'] - 1) == 0
        # Text-Wert unverändert
        assert result['c'] == 'text'

    def test_simplify_all_with_numbers(self):
        """Test: simplify_all() verarbeitet auch numerische Werte."""
        data = {'pi_approx': 3.14159, 'integer': 42}
        result = self.fs.simplify_all(data)
        # Zahlen sollten als SymPy-Ausdrücke vorliegen oder unverändert
        assert result['integer'] is not None

    def test_simplify_all_with_list(self):
        """Test: simplify_all() vereinfacht Listen element-weise."""
        data = {
            'roots': [sp.sin(self.x)**2 + sp.cos(self.x)**2, self.x**2 * self.x**0]
        }
        result = self.fs.simplify_all(data)
        # Erste Element sollte zu 1 vereinfacht werden
        assert sp.simplify(result['roots'][0] - 1) == 0

    def test_simplify_auto_rational(self):
        """Test: simplify_auto() erkennt rationale Approximationen."""
        # sqrt(4) = 2
        expr = sp.sqrt(4)
        result = self.fs.simplify_auto(expr)
        assert result == sp.Integer(2)

    def test_to_latex_integral(self):
        """Test: Integral wird korrekt in LaTeX konvertiert."""
        expr = sp.Integral(self.x**2, self.x)
        latex = self.fs.to_latex(expr)
        # Integral-Symbol muss vorhanden sein
        assert 'int' in latex.lower() or '\\int' in latex

    def test_simplify_nested_expression(self):
        """Test: Verschachtelter Ausdruck wird korrekt vereinfacht."""
        # (a+b)² = a²+2ab+b²
        a, b = sp.symbols('a b')
        expr = (a + b)**2 - a**2 - 2*a*b - b**2
        result = self.fs.simplify(expr)
        assert result == sp.Integer(0)
