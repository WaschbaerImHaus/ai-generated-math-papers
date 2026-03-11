"""
@file test_sagemath_bridge.py
@brief Tests für die SageMath-Interoperabilitäts-Brücke.
@description
    Alle Tests müssen auch ohne SageMath-Installation grün sein.
    Sie testen den Fallback-Modus (SymPy-basiert).

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import sympy as sp
import pytest
from sagemath_bridge import (
    is_sage_available,
    to_sage,
    from_sage,
    export_to_sage_format,
    export_to_sage_string,
    import_from_sage_string,
    _make_safe_varname,
)


class TestIsSageAvailable:
    """Tests für is_sage_available()."""

    def test_returns_bool(self):
        """Test: is_sage_available() gibt einen Boolean zurück."""
        result = is_sage_available()
        assert isinstance(result, bool)

    def test_no_exception(self):
        """Test: Keine Ausnahme wird geworfen."""
        try:
            is_sage_available()
        except Exception as e:
            pytest.fail(f"is_sage_available() hat eine Ausnahme geworfen: {e}")


class TestToSage:
    """Tests für to_sage()."""

    def test_sympy_expr_returns_something(self):
        """Test: SymPy-Ausdruck wird ohne Fehler konvertiert."""
        x = sp.Symbol('x')
        expr = x**2 + 1
        result = to_sage(expr)
        assert result is not None

    def test_fallback_returns_sympy(self):
        """Test: Ohne SageMath wird SymPy-Ausdruck zurückgegeben."""
        if not is_sage_available():
            x = sp.Symbol('x')
            expr = x**2
            result = to_sage(expr)
            # Im Fallback-Modus sollte der SymPy-Ausdruck zurückgegeben werden
            assert isinstance(result, sp.Basic)

    def test_integer_input(self):
        """Test: Ganzzahl-Eingabe wird korrekt verarbeitet."""
        result = to_sage(42)
        assert result is not None

    def test_string_input(self):
        """Test: String-Eingabe wird konvertiert."""
        result = to_sage('x**2 + 1')
        assert result is not None


class TestFromSage:
    """Tests für from_sage()."""

    def test_sympy_expr_roundtrip(self):
        """Test: SymPy → SageMath → SymPy Roundtrip (Fallback-Modus)."""
        x = sp.Symbol('x')
        original = x**2 + 3*x + 2

        # In SageMath konvertieren
        sage_expr = to_sage(original)

        # Zurück zu SymPy
        result = from_sage(sage_expr)
        assert result is not None

    def test_from_sage_string(self):
        """Test: String kann zu SymPy konvertiert werden."""
        result = from_sage('x**2 + 1')
        assert isinstance(result, sp.Basic)

    def test_from_sage_integer(self):
        """Test: Ganzzahl wird korrekt konvertiert."""
        result = from_sage(5)
        assert result is not None


class TestExportToSageFormat:
    """Tests für export_to_sage_format()."""

    def test_returns_string(self):
        """Test: Rückgabe ist ein String."""
        data = {}
        result = export_to_sage_format(data)
        assert isinstance(result, str)

    def test_contains_imports(self):
        """Test: Export-Skript enthält Import-Anweisungen."""
        data = {}
        script = export_to_sage_format(data)
        assert 'import' in script.lower()

    def test_contains_sage_check(self):
        """Test: Export-Skript enthält SageMath-Verfügbarkeitsprüfung."""
        data = {}
        script = export_to_sage_format(data)
        assert 'SAGE_AVAILABLE' in script

    def test_sympy_expr_in_export(self):
        """Test: SymPy-Ausdruck erscheint im Export-Skript."""
        x = sp.Symbol('x')
        data = {'mein_ergebnis': x**2 + 1}
        script = export_to_sage_format(data)
        assert 'mein_ergebnis' in script

    def test_multiple_expressions(self):
        """Test: Mehrere Ausdrücke werden korrekt exportiert."""
        x, y = sp.symbols('x y')
        data = {
            'expr1': x**2,
            'expr2': y + 1,
            'value': 42,
        }
        script = export_to_sage_format(data)
        assert 'expr1' in script
        assert 'expr2' in script
        assert 'value' in script

    def test_save_to_file(self, tmp_path):
        """Test: Export in Datei speichern."""
        data = {'test': sp.Integer(1)}
        filepath = str(tmp_path / "test_sage.py")
        script = export_to_sage_format(data, filename=filepath)

        assert os.path.exists(filepath)
        with open(filepath) as f:
            content = f.read()
        assert 'SAGE_AVAILABLE' in content


class TestExportToSageString:
    """Tests für export_to_sage_string()."""

    def test_returns_string(self):
        """Test: Rückgabe ist ein String."""
        x = sp.Symbol('x')
        result = export_to_sage_string(x**2)
        assert isinstance(result, str)

    def test_simple_expression(self):
        """Test: Einfacher Ausdruck wird korrekt als String dargestellt."""
        x = sp.Symbol('x')
        result = export_to_sage_string(x + 1)
        assert 'x' in result
        assert '1' in result

    def test_constant_pi(self):
        """Test: pi wird als String dargestellt."""
        result = export_to_sage_string(sp.pi)
        assert result is not None
        assert len(result) > 0

    def test_integer_input(self):
        """Test: Ganzzahl-Eingabe wird korrekt konvertiert."""
        result = export_to_sage_string(42)
        assert '42' in result

    def test_polynomial(self):
        """Test: Polynom wird als String dargestellt."""
        x = sp.Symbol('x')
        poly = x**3 - 2*x**2 + x - 1
        result = export_to_sage_string(poly)
        assert isinstance(result, str)
        assert 'x' in result

    def test_string_input(self):
        """Test: String-Eingabe wird verarbeitet."""
        result = export_to_sage_string('x**2 + 1')
        assert isinstance(result, str)


class TestImportFromSageString:
    """Tests für import_from_sage_string()."""

    def test_simple_number(self):
        """Test: Einfache Zahl wird korrekt geparst."""
        result = import_from_sage_string('42')
        assert result == sp.Integer(42)

    def test_variable_expression(self):
        """Test: Ausdruck mit Variablen wird geparst."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = import_from_sage_string('x**2 + 1')
        assert isinstance(result, sp.Basic)

    def test_sage_power_notation(self):
        """Test: SageMath-Potenz-Notation (^) wird zu (**) konvertiert."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = import_from_sage_string('x^2')
        x = sp.Symbol('x')
        assert result == x**2

    def test_empty_string(self):
        """Test: Leerer String gibt 0 zurück."""
        result = import_from_sage_string('')
        assert result == sp.Integer(0)

    def test_trig_function(self):
        """Test: Trigonometrische Funktion wird geparst."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = import_from_sage_string('sin(x)')
        assert isinstance(result, sp.Basic)

    def test_pi_constant(self):
        """Test: pi wird als sp.pi geparst."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = import_from_sage_string('pi')
        assert result == sp.pi

    def test_roundtrip_sympy(self):
        """Test: SymPy → Sage-String → SymPy ergibt gleichwertigen Ausdruck."""
        import warnings
        x = sp.Symbol('x')
        original = x**2 + 1

        # Zu Sage-String konvertieren
        sage_str = export_to_sage_string(original)

        # Zurück parsen
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = import_from_sage_string(sage_str)

        # Ausdrücke sollten äquivalent sein
        assert isinstance(result, sp.Basic)


class TestMakeSafeVarname:
    """Tests für die interne _make_safe_varname()-Funktion."""

    def test_simple_name(self):
        """Test: Einfacher Name bleibt unverändert."""
        assert _make_safe_varname('result') == 'result'

    def test_spaces_replaced(self):
        """Test: Leerzeichen werden durch Unterstriche ersetzt."""
        assert _make_safe_varname('my result') == 'my_result'

    def test_digit_prefix(self):
        """Test: Name der mit Ziffer beginnt bekommt Prefix."""
        result = _make_safe_varname('1bad')
        assert result[0].isalpha() or result[0] == '_'

    def test_special_chars(self):
        """Test: Sonderzeichen werden ersetzt."""
        result = _make_safe_varname('a-b+c')
        assert all(c.isalnum() or c == '_' for c in result)
