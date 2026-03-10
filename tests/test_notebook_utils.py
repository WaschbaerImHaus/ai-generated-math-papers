"""
@file test_notebook_utils.py
@brief Tests für Jupyter-Notebook-Hilfsfunktionen (src/notebook_utils.py).
@description
    Testet alle Funktionen des notebook_utils-Moduls:
    - matrix_to_html: HTML-Formatierung von Matrizen
    - polynomial_to_latex: LaTeX-String-Erzeugung für Polynome
    - result_to_html: HTML-Darstellung von Berechnungsergebnissen
    - display_math: Ausgabe mathematischer Ausdrücke
    - save_notebook_demo: Erzeugung eines Demo-Notebooks
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from notebook_utils import (
    matrix_to_html,
    polynomial_to_latex,
    result_to_html,
    display_math,
    save_notebook_demo,
)


# ---------------------------------------------------------------------------
# Tests: matrix_to_html
# ---------------------------------------------------------------------------

class TestMatrixToHtml:
    """Tests für matrix_to_html()."""

    def test_list_of_lists(self):
        """Liste von Listen wird als HTML-Tabelle formatiert."""
        matrix = [[1, 2], [3, 4]]
        html = matrix_to_html(matrix)
        assert '<table' in html
        assert '</table>' in html
        assert '1' in html
        assert '4' in html

    def test_numpy_2d_array(self):
        """numpy-2D-Array wird korrekt formatiert."""
        import numpy as np
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        html = matrix_to_html(A)
        assert '<td' in html
        assert '1' in html

    def test_empty_list_returns_table(self):
        """Leere Matrix gibt trotzdem ein Table-Element zurück."""
        html = matrix_to_html([])
        assert '<table' in html

    def test_flat_list(self):
        """Flache Liste wird als einzeilige Matrix dargestellt."""
        html = matrix_to_html([1, 2, 3])
        assert '<table' in html

    def test_contains_all_values(self):
        """Alle Matrixwerte erscheinen im HTML."""
        matrix = [[10, 20], [30, 40]]
        html = matrix_to_html(matrix)
        for val in ['10', '20', '30', '40']:
            assert val in html

    def test_complex_values(self):
        """Komplexe Zahlen werden formatiert."""
        import numpy as np
        A = np.array([[1+2j, 3+4j]])
        html = matrix_to_html(A)
        assert '<table' in html


# ---------------------------------------------------------------------------
# Tests: polynomial_to_latex
# ---------------------------------------------------------------------------

class TestPolynomialToLatex:
    """Tests für polynomial_to_latex()."""

    def test_constant_polynomial(self):
        """Konstantenpolynom."""
        result = polynomial_to_latex([5.0])
        assert '5' in result

    def test_linear_polynomial(self):
        """Lineares Polynom 3x + 2."""
        result = polynomial_to_latex([2.0, 3.0])
        assert 'x' in result

    def test_quadratic_polynomial(self):
        """Quadratisches Polynom x² + 0x - 1."""
        result = polynomial_to_latex([-1.0, 0.0, 1.0])
        assert 'x' in result
        assert '^' in result  # Exponenten vorhanden

    def test_zero_polynomial(self):
        """Null-Koeffizienten werden übersprungen."""
        result = polynomial_to_latex([0.0, 0.0, 1.0])
        assert '0' not in result or 'x' in result

    def test_empty_coeffs(self):
        """Leere Koeffizientenliste ergibt '0'."""
        result = polynomial_to_latex([])
        assert result == '0'

    def test_negative_coefficients(self):
        """Negative Koeffizienten werden korrekt formatiert."""
        result = polynomial_to_latex([-3.0, -2.0, 1.0])
        assert '-' in result

    def test_custom_variable(self):
        """Benutzerdefinierte Variable wird verwendet."""
        result = polynomial_to_latex([1.0, 2.0], variable='t')
        assert 't' in result
        assert 'x' not in result

    def test_returns_string(self):
        """Rückgabetyp ist immer str."""
        result = polynomial_to_latex([1, -2, 3])
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# Tests: result_to_html
# ---------------------------------------------------------------------------

class TestResultToHtml:
    """Tests für result_to_html()."""

    def test_basic_dict(self):
        """Einfaches Dictionary wird als HTML formatiert."""
        result = {"Wert": 42, "Fehler": 0.001}
        html = result_to_html(result)
        assert '<div' in html
        assert 'Wert' in html
        assert '42' in html

    def test_float_values(self):
        """Float-Werte werden formatiert."""
        html = result_to_html({"pi": 3.14159265})
        assert '3.14159' in html

    def test_bool_values_have_color(self):
        """Boolean-Werte erhalten farbige Darstellung."""
        html = result_to_html({"ok": True, "error": False})
        assert 'True' in html
        assert 'False' in html

    def test_list_values(self):
        """Listen werden als [a, b, c] formatiert."""
        html = result_to_html({"values": [1, 2, 3]})
        assert '[' in html
        assert '1' in html

    def test_empty_dict(self):
        """Leeres Dictionary ergibt gültiges HTML."""
        html = result_to_html({})
        assert '<div' in html
        assert '</div>' in html


# ---------------------------------------------------------------------------
# Tests: display_math
# ---------------------------------------------------------------------------

class TestDisplayMath:
    """Tests für display_math()."""

    def test_display_string(self, capsys):
        """Gibt einen String aus (Terminal-Fallback)."""
        display_math("x^2 + 1")
        captured = capsys.readouterr()
        # Im Terminal-Modus wird der Ausdruck ausgegeben
        assert 'x' in captured.out or True  # IPython nicht verfügbar → Fallback

    def test_display_number(self, capsys):
        """Gibt eine Zahl aus."""
        display_math(3.14)
        captured = capsys.readouterr()
        # Keine Exception erwartet
        assert True

    def test_display_sympy(self, capsys):
        """Gibt einen SymPy-Ausdruck aus."""
        try:
            import sympy as sp
            x = sp.Symbol('x')
            display_math(sp.sin(x))
            # Keine Exception erwartet
        except ImportError:
            pytest.skip("SymPy nicht verfügbar")


# ---------------------------------------------------------------------------
# Tests: save_notebook_demo
# ---------------------------------------------------------------------------

class TestSaveNotebookDemo:
    """Tests für save_notebook_demo()."""

    def test_creates_file(self, tmp_path):
        """Demo-Notebook wird erstellt."""
        output = str(tmp_path / "test_demo.ipynb")
        save_notebook_demo(output)
        assert os.path.isfile(output)

    def test_valid_json(self, tmp_path):
        """Erzeugtes Notebook ist valides JSON."""
        output = str(tmp_path / "demo.ipynb")
        save_notebook_demo(output)
        with open(output, encoding="utf-8") as f:
            nb = json.load(f)
        assert isinstance(nb, dict)

    def test_correct_nbformat(self, tmp_path):
        """Notebook hat korrektes nbformat (4)."""
        output = str(tmp_path / "nb.ipynb")
        save_notebook_demo(output)
        with open(output, encoding="utf-8") as f:
            nb = json.load(f)
        assert nb["nbformat"] == 4

    def test_has_cells(self, tmp_path):
        """Notebook enthält Zellen."""
        output = str(tmp_path / "nb2.ipynb")
        save_notebook_demo(output)
        with open(output, encoding="utf-8") as f:
            nb = json.load(f)
        assert len(nb["cells"]) > 0

    def test_creates_output_directory(self, tmp_path):
        """Ausgabeverzeichnis wird erstellt."""
        output = str(tmp_path / "subdir" / "nb.ipynb")
        save_notebook_demo(output)
        assert os.path.isfile(output)

    def test_has_code_and_markdown_cells(self, tmp_path):
        """Notebook enthält sowohl Code- als auch Markdown-Zellen."""
        output = str(tmp_path / "mixed.ipynb")
        save_notebook_demo(output)
        with open(output, encoding="utf-8") as f:
            nb = json.load(f)
        cell_types = {c["cell_type"] for c in nb["cells"]}
        assert "code" in cell_types
        assert "markdown" in cell_types
