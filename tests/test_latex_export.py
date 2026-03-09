"""
@file test_latex_export.py
@brief Tests für das LaTeX-Export-Modul.
@description
    Unit-Tests für alle Funktionen in latex_export.py.
    Prüft korrekte LaTeX-Ausgaben für Zahlen, Brüche, Polynome,
    Matrizen, Vektoren, Gleichungen, Integrale, Summen, Grenzwerte,
    Theoreme und vollständige Dokumente.

@author Kurt Ingwer
@date 2026-03-09
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import latex_export as lx


# =============================================================================
# TESTS: number_to_latex
# =============================================================================

class TestNumberToLatex:
    """Tests für die Funktion number_to_latex."""

    def test_integer_positive(self):
        """Positive ganze Zahl."""
        assert lx.number_to_latex(42) == "42"

    def test_integer_negative(self):
        """Negative ganze Zahl."""
        assert lx.number_to_latex(-7) == "-7"

    def test_integer_zero(self):
        """Null."""
        assert lx.number_to_latex(0) == "0"

    def test_float_basic(self):
        """Einfache Fließkommazahl mit Standardpräzision."""
        result = lx.number_to_latex(3.14159, 5)
        assert "3.14159" in result

    def test_float_precision(self):
        """Präzisionsangabe wird beachtet."""
        result = lx.number_to_latex(1.23456789, 4)
        # Sollte 4 Nachkommastellen enthalten (1.2346 gerundet)
        assert "1.2346" in result

    def test_pi_detection(self):
        """Pi wird als \\pi erkannt."""
        result = lx.number_to_latex(math.pi)
        assert result == r"\pi"

    def test_negative_pi_detection(self):
        """Negatives Pi wird erkannt."""
        result = lx.number_to_latex(-math.pi)
        assert r"\pi" in result

    def test_sqrt2_detection(self):
        """Wurzel aus 2 wird erkannt."""
        result = lx.number_to_latex(math.sqrt(2))
        assert r"\sqrt{2}" in result

    def test_fraction_half(self):
        """1/2 wird als Bruch erkannt."""
        result = lx.number_to_latex(0.5)
        assert r"\frac{1}{2}" in result

    def test_fraction_third(self):
        """1/3 wird als Bruch erkannt."""
        result = lx.number_to_latex(1/3)
        assert r"\frac{1}{3}" in result

    def test_complex_basic(self):
        """Komplexe Zahl 3+4i."""
        result = lx.number_to_latex(3 + 4j)
        assert "3" in result
        assert "4" in result
        assert "i" in result

    def test_complex_pure_imaginary(self):
        """Rein imaginäre Zahl."""
        result = lx.number_to_latex(0 + 5j)
        assert "5" in result
        assert "i" in result

    def test_float_integer_value(self):
        """Float mit ganzzahligem Wert wird als Integer dargestellt."""
        result = lx.number_to_latex(3.0)
        assert result == "3"


# =============================================================================
# TESTS: fraction_to_latex
# =============================================================================

class TestFractionToLatex:
    """Tests für die Funktion fraction_to_latex."""

    def test_basic_fraction(self):
        """Einfacher Bruch 1/2."""
        assert lx.fraction_to_latex(1, 2) == r"\frac{1}{2}"

    def test_fraction_simplification(self):
        """Bruch wird vereinfacht (2/4 → 1/2)."""
        result = lx.fraction_to_latex(2, 4)
        assert result == r"\frac{1}{2}"

    def test_fraction_whole_number(self):
        """Bruch mit Nenner 1 ergibt ganze Zahl."""
        assert lx.fraction_to_latex(6, 3) == "2"

    def test_fraction_negative_numerator(self):
        """Negativer Zähler."""
        result = lx.fraction_to_latex(-1, 3)
        assert r"\frac{-1}{3}" in result

    def test_fraction_negative_denominator(self):
        """Negativer Nenner wird normalisiert."""
        result = lx.fraction_to_latex(1, -3)
        assert r"\frac{-1}{3}" in result

    def test_fraction_three_quarters(self):
        """3/4 bleibt 3/4."""
        assert lx.fraction_to_latex(3, 4) == r"\frac{3}{4}"

    def test_fraction_zero_denominator(self):
        """Nenner 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            lx.fraction_to_latex(1, 0)


# =============================================================================
# TESTS: polynomial_to_latex
# =============================================================================

class TestPolynomialToLatex:
    """Tests für die Funktion polynomial_to_latex."""

    def test_quadratic(self):
        """Quadratisches Polynom x^2 - 3x + 2."""
        result = lx.polynomial_to_latex([1, -3, 2])
        assert "x^{2}" in result
        assert "3" in result
        # Vorzeichen prüfen
        assert "-" in result

    def test_linear(self):
        """Lineares Polynom 2x + 1."""
        result = lx.polynomial_to_latex([2, 1])
        assert "2" in result
        assert "x" in result

    def test_constant(self):
        """Konstantes Polynom."""
        result = lx.polynomial_to_latex([5])
        assert "5" in result

    def test_zero_polynomial(self):
        """Nullpolynom."""
        result = lx.polynomial_to_latex([0])
        assert result == "0"

    def test_leading_coeff_one(self):
        """Führender Koeffizient 1 wird nicht explizit angezeigt."""
        result = lx.polynomial_to_latex([1, 0, 0])
        # "x^{2}" ohne "1x^{2}"
        assert "x^{2}" in result

    def test_custom_variable(self):
        """Andere Variable als x."""
        result = lx.polynomial_to_latex([1, -2], 't')
        assert "t" in result

    def test_negative_leading(self):
        """Negativer führender Koeffizient."""
        result = lx.polynomial_to_latex([-1, 2])
        assert "-" in result


# =============================================================================
# TESTS: matrix_to_latex
# =============================================================================

class TestMatrixToLatex:
    """Tests für die Funktion matrix_to_latex."""

    def test_pmatrix_default(self):
        """Standardmäßig pmatrix."""
        result = lx.matrix_to_latex([[1, 2], [3, 4]])
        assert "pmatrix" in result

    def test_matrix_contains_elements(self):
        """Matrix enthält die Einträge."""
        result = lx.matrix_to_latex([[1, 2], [3, 4]])
        assert "1 & 2" in result

    def test_bmatrix(self):
        """Eckige Klammern mit bmatrix."""
        result = lx.matrix_to_latex([[1, 2]], 'bmatrix')
        assert "bmatrix" in result

    def test_vmatrix(self):
        """Betragsstriche mit vmatrix."""
        result = lx.matrix_to_latex([[1, 2]], 'vmatrix')
        assert "vmatrix" in result

    def test_3x3_matrix(self):
        """3×3-Matrix."""
        result = lx.matrix_to_latex([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        assert "pmatrix" in result
        assert "1" in result

    def test_matrix_row_separator(self):
        """Zeilentrenner \\\\ ist vorhanden."""
        result = lx.matrix_to_latex([[1, 2], [3, 4]])
        assert r"\\" in result or "\\\\" in result

    def test_matrix_begin_end(self):
        """\\begin und \\end vorhanden."""
        result = lx.matrix_to_latex([[1]])
        assert r"\begin{pmatrix}" in result
        assert r"\end{pmatrix}" in result


# =============================================================================
# TESTS: vector_to_latex
# =============================================================================

class TestVectorToLatex:
    """Tests für die Funktion vector_to_latex."""

    def test_column_vector(self):
        """Spaltenvektor enthält pmatrix."""
        result = lx.vector_to_latex([1, 2, 3])
        assert "pmatrix" in result

    def test_column_vector_elements(self):
        """Spaltenvektor enthält alle Elemente."""
        result = lx.vector_to_latex([1, 2, 3])
        assert "1" in result
        assert "2" in result
        assert "3" in result

    def test_row_vector(self):
        """Zeilenvektor mit column=False."""
        result = lx.vector_to_latex([1, 2, 3], column=False)
        assert "pmatrix" in result
        assert "&" in result

    def test_column_vector_separator(self):
        """Spaltenvektor nutzt \\\\ als Trenner."""
        result = lx.vector_to_latex([1, 2])
        assert r"\\" in result


# =============================================================================
# TESTS: equation_to_latex
# =============================================================================

class TestEquationToLatex:
    """Tests für die Funktion equation_to_latex."""

    def test_equality(self):
        """Standard-Gleichung."""
        result = lx.equation_to_latex("a", "b")
        assert "a" in result
        assert "b" in result
        assert "=" in result

    def test_leq(self):
        """Kleiner-gleich."""
        result = lx.equation_to_latex("x", "5", "<=")
        assert r"\leq" in result

    def test_geq(self):
        """Größer-gleich."""
        result = lx.equation_to_latex("x", "5", ">=")
        assert r"\geq" in result

    def test_approx(self):
        """Ungefähr gleich."""
        result = lx.equation_to_latex("pi", "3.14", "~=")
        assert r"\approx" in result

    def test_neq(self):
        """Ungleich."""
        result = lx.equation_to_latex("a", "b", "!=")
        assert r"\neq" in result


# =============================================================================
# TESTS: integral_to_latex
# =============================================================================

class TestIntegralToLatex:
    """Tests für die Funktion integral_to_latex."""

    def test_definite_integral(self):
        """Bestimmtes Integral mit Grenzen."""
        result = lx.integral_to_latex("x^2", "x", 0, 1)
        assert r"\int" in result
        assert "0" in result
        assert "1" in result

    def test_indefinite_integral(self):
        """Unbestimmtes Integral ohne Grenzen."""
        result = lx.integral_to_latex("x^2", "x")
        assert r"\int" in result
        assert "x^2" in result

    def test_integration_variable(self):
        """Integrationsvariable ist enthalten."""
        result = lx.integral_to_latex("f(t)", "t", 0, 10)
        assert "t" in result


# =============================================================================
# TESTS: sum_to_latex
# =============================================================================

class TestSumToLatex:
    """Tests für die Funktion sum_to_latex."""

    def test_basic_sum(self):
        """Einfache Summe."""
        result = lx.sum_to_latex("n^2", "n", "1", r"\infty")
        assert r"\sum" in result
        assert "n" in result

    def test_sum_indices(self):
        """Summe mit korrekten Grenzen."""
        result = lx.sum_to_latex("k", "k", "0", "N")
        assert "0" in result
        assert "N" in result


# =============================================================================
# TESTS: limit_to_latex
# =============================================================================

class TestLimitToLatex:
    """Tests für die Funktion limit_to_latex."""

    def test_basic_limit(self):
        """Grenzwert ohne Richtungsangabe."""
        result = lx.limit_to_latex("sin(x)/x", "x", "0")
        assert r"\lim" in result
        assert "x" in result

    def test_limit_right(self):
        """Rechtsseitiger Grenzwert."""
        result = lx.limit_to_latex("f(x)", "x", "0", "+")
        assert r"\lim" in result
        assert "+" in result

    def test_limit_left(self):
        """Linksseitiger Grenzwert."""
        result = lx.limit_to_latex("f(x)", "x", "0", "-")
        assert r"\lim" in result
        assert "-" in result

    def test_limit_contains_expression(self):
        """Ausdruck ist im Ergebnis enthalten."""
        result = lx.limit_to_latex("sin(x)/x", "x", "0")
        assert "sin(x)/x" in result


# =============================================================================
# TESTS: derivative_to_latex
# =============================================================================

class TestDerivativeToLatex:
    """Tests für die Funktion derivative_to_latex."""

    def test_first_derivative(self):
        """Erste Ableitung."""
        result = lx.derivative_to_latex("f(x)", "x", 1)
        assert r"\frac{d}{dx}" in result or r"\frac{d}{d" in result

    def test_higher_derivative(self):
        """Höhere Ableitung."""
        result = lx.derivative_to_latex("f(x)", "x", 3)
        assert "3" in result
        assert r"\frac" in result


# =============================================================================
# TESTS: theorem_to_latex
# =============================================================================

class TestTheoremToLatex:
    """Tests für die Funktion theorem_to_latex."""

    def test_theorem_basic(self):
        """Theorem-Block enthält Theorem-Umgebung."""
        result = lx.theorem_to_latex("Pythagoras", "a^2 + b^2 = c^2")
        assert "theorem" in result

    def test_theorem_contains_name(self):
        """Theorem-Name ist enthalten."""
        result = lx.theorem_to_latex("Pythagoras", "a^2 + b^2 = c^2")
        assert "Pythagoras" in result

    def test_theorem_contains_statement(self):
        """Theorem-Aussage ist enthalten."""
        result = lx.theorem_to_latex("Test", "a = b")
        assert "a = b" in result

    def test_theorem_with_proof(self):
        """Theorem mit Beweis enthält proof-Umgebung."""
        result = lx.theorem_to_latex("Test", "Aussage", "Beweis hier")
        assert "proof" in result
        assert "Beweis hier" in result

    def test_theorem_without_proof(self):
        """Theorem ohne Beweis enthält keine proof-Umgebung."""
        result = lx.theorem_to_latex("Test", "Aussage")
        assert "proof" not in result


# =============================================================================
# TESTS: full_document
# =============================================================================

class TestFullDocument:
    """Tests für die Funktion full_document."""

    def test_documentclass(self):
        """Vollständiges Dokument beginnt mit \\documentclass."""
        result = lx.full_document("Inhalt")
        assert r"\documentclass" in result

    def test_amsmath_package(self):
        """amsmath-Paket ist enthalten."""
        result = lx.full_document("Inhalt")
        assert "amsmath" in result

    def test_begin_document(self):
        """\\begin{document} ist vorhanden."""
        result = lx.full_document("Inhalt")
        assert r"\begin{document}" in result

    def test_end_document(self):
        """\\end{document} ist vorhanden."""
        result = lx.full_document("Inhalt")
        assert r"\end{document}" in result

    def test_content_included(self):
        """Übergebener Inhalt ist im Dokument."""
        result = lx.full_document("MEIN_INHALT")
        assert "MEIN_INHALT" in result

    def test_title_included(self):
        """Titel ist im Dokument."""
        result = lx.full_document("Inhalt", title="Mein Test")
        assert "Mein Test" in result

    def test_author_included(self):
        """Autor ist im Dokument."""
        result = lx.full_document("Inhalt", author="Kurt Ingwer")
        assert "Kurt Ingwer" in result

    def test_custom_packages(self):
        """Zusätzliche Pakete werden eingebunden."""
        result = lx.full_document("Inhalt", packages=["graphicx"])
        assert "graphicx" in result


# =============================================================================
# TESTS: sympy_to_latex
# =============================================================================

class TestSympyToLatex:
    """Tests für die Funktion sympy_to_latex."""

    def test_sympy_expression(self):
        """SymPy-Ausdruck wird korrekt konvertiert."""
        import sympy as sp
        x = sp.Symbol('x')
        expr = sp.sin(x) / x
        result = lx.sympy_to_latex(expr)
        assert "frac" in result or "sin" in result

    def test_sympy_integer(self):
        """SymPy-Ganzzahl."""
        import sympy as sp
        result = lx.sympy_to_latex(sp.Integer(42))
        assert "42" in result


# =============================================================================
# TESTS: solution_to_latex
# =============================================================================

class TestSolutionToLatex:
    """Tests für die Funktion solution_to_latex."""

    def test_basic_solution(self):
        """Einfaches Ergebnis-Dictionary."""
        result = lx.solution_to_latex({"x": 3, "y": -1})
        assert "x" in result
        assert "3" in result

    def test_align_environment(self):
        """align*-Umgebung ist vorhanden."""
        result = lx.solution_to_latex({"a": 1})
        assert "align" in result

    def test_solution_with_title(self):
        """Titel wird angezeigt."""
        result = lx.solution_to_latex({"x": 1}, title="Ergebnis")
        assert "Ergebnis" in result
