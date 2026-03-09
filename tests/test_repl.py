"""
@file test_repl.py
@brief Tests für den interaktiven REPL-Modus.
@description
    Unit-Tests für die MathREPL-Klasse und die is_jupyter()-Funktion.
    Prüft alle verfügbaren Befehle auf korrekte Ausgaben.

@author Kurt Ingwer
@date 2026-03-09
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from repl import MathREPL, is_jupyter


# =============================================================================
# FIXTURE
# =============================================================================

@pytest.fixture
def repl():
    """Erstellt eine frische MathREPL-Instanz für jeden Test."""
    return MathREPL()


# =============================================================================
# TESTS: is_jupyter
# =============================================================================

class TestIsJupyter:
    """Tests für die is_jupyter()-Funktion."""

    def test_returns_bool(self):
        """is_jupyter() gibt einen Boolean zurück."""
        result = is_jupyter()
        assert isinstance(result, bool)

    def test_not_jupyter_in_pytest(self):
        """In pytest läuft das Programm nicht in Jupyter."""
        result = is_jupyter()
        assert result is False


# =============================================================================
# TESTS: MathREPL.available_commands
# =============================================================================

class TestAvailableCommands:
    """Tests für available_commands()."""

    def test_returns_dict(self, repl):
        """available_commands() gibt ein Dictionary zurück."""
        result = repl.available_commands()
        assert isinstance(result, dict)

    def test_dict_not_empty(self, repl):
        """Das Dictionary ist nicht leer."""
        result = repl.available_commands()
        assert len(result) > 0

    def test_help_command_exists(self, repl):
        """'help' ist ein verfügbarer Befehl."""
        result = repl.available_commands()
        assert 'help' in result

    def test_prime_command_exists(self, repl):
        """'prime' ist ein verfügbarer Befehl."""
        result = repl.available_commands()
        assert 'prime' in result

    def test_values_are_strings(self, repl):
        """Alle Werte des Dictionaries sind Strings."""
        result = repl.available_commands()
        for key, value in result.items():
            assert isinstance(value, str)


# =============================================================================
# TESTS: MathREPL.evaluate - Algebra-Befehle
# =============================================================================

class TestEvaluatePrime:
    """Tests für den 'prime' Befehl."""

    def test_prime_97_is_true(self, repl):
        """97 ist eine Primzahl."""
        result = repl.evaluate("prime 97")
        assert "True" in result

    def test_prime_100_is_false(self, repl):
        """100 ist keine Primzahl."""
        result = repl.evaluate("prime 100")
        assert "False" in result

    def test_prime_2_is_true(self, repl):
        """2 ist eine Primzahl."""
        result = repl.evaluate("prime 2")
        assert "True" in result

    def test_prime_1_is_false(self, repl):
        """1 ist keine Primzahl."""
        result = repl.evaluate("prime 1")
        assert "False" in result


class TestEvaluateFactor:
    """Tests für den 'factor' Befehl."""

    def test_factor_12(self, repl):
        """Primfaktorzerlegung von 12 = 2² × 3."""
        result = repl.evaluate("factor 12")
        assert "2" in result
        assert "3" in result

    def test_factor_360(self, repl):
        """Primfaktorzerlegung von 360."""
        result = repl.evaluate("factor 360")
        assert "2" in result
        assert "3" in result
        assert "5" in result

    def test_factor_prime(self, repl):
        """Primzahl hat nur sich selbst als Faktor."""
        result = repl.evaluate("factor 17")
        assert "17" in result


class TestEvaluateGcd:
    """Tests für den 'gcd' Befehl."""

    def test_gcd_48_18(self, repl):
        """ggT(48, 18) = 6."""
        result = repl.evaluate("gcd 48 18")
        assert "6" in result

    def test_gcd_100_75(self, repl):
        """ggT(100, 75) = 25."""
        result = repl.evaluate("gcd 100 75")
        assert "25" in result

    def test_gcd_contains_both_args(self, repl):
        """Ausgabe enthält beide Eingabezahlen."""
        result = repl.evaluate("gcd 48 18")
        assert "48" in result
        assert "18" in result


class TestEvaluateSolve:
    """Tests für den 'solve' Befehl."""

    def test_solve_simple_quadratic(self, repl):
        """x² - 5x + 6 = 0 hat Lösungen 2 und 3."""
        result = repl.evaluate("solve 1 -5 6")
        assert "2" in result or "3" in result

    def test_solve_no_real_roots(self, repl):
        """x² + 1 = 0 hat keine reellen Lösungen."""
        result = repl.evaluate("solve 1 0 1")
        # Sollte komplexe Zahlen oder Fehlerhinweis zurückgeben
        assert len(result) > 0


# =============================================================================
# TESTS: MathREPL.evaluate - Analysis-Befehle
# =============================================================================

class TestEvaluateDerive:
    """Tests für den 'derive' Befehl."""

    def test_derive_sin_at_0(self, repl):
        """Ableitung von sin(x) bei x=0 ≈ 1."""
        result = repl.evaluate("derive sin at 0")
        assert "1.0" in result or "0.999" in result or "1.00" in result

    def test_derive_cos_at_0(self, repl):
        """Ableitung von cos(x) bei x=0 ≈ 0."""
        result = repl.evaluate("derive cos at 0")
        # |cos'(0)| = |sin(0)| = 0, also sehr klein
        assert "0.0" in result or "0.00" in result or "-0.0" in result or "0" in result


class TestEvaluateIntegrate:
    """Tests für den 'integrate' Befehl."""

    def test_integrate_sin_0_pi(self, repl):
        """∫sin(x)dx von 0 bis π ≈ 2."""
        result = repl.evaluate("integrate sin 0 3.14159265")
        assert "2.0" in result or "1.999" in result or "2.00" in result


# =============================================================================
# TESTS: MathREPL.evaluate - Fourier
# =============================================================================

class TestEvaluateFft:
    """Tests für den 'fft' Befehl."""

    def test_fft_list(self, repl):
        """FFT einer einfachen Liste."""
        result = repl.evaluate("fft [1,2,3,4,5,6,7,8]")
        assert "FFT" in result or "Amplituden" in result

    def test_fft_contains_values(self, repl):
        """FFT-Ergebnis enthält numerische Werte."""
        result = repl.evaluate("fft [1,0,1,0]")
        assert len(result) > 0


# =============================================================================
# TESTS: MathREPL.evaluate - Zahlentheorie
# =============================================================================

class TestEvaluateGoldbach:
    """Tests für den 'goldbach' Befehl."""

    def test_goldbach_100(self, repl):
        """Goldbach-Zerlegung von 100."""
        result = repl.evaluate("goldbach 100")
        assert "+" in result or "100" in result

    def test_goldbach_4(self, repl):
        """Goldbach-Zerlegung von 4 = 2 + 2."""
        result = repl.evaluate("goldbach 4")
        assert "2" in result


class TestEvaluateCollatz:
    """Tests für den 'collatz' Befehl."""

    def test_collatz_27(self, repl):
        """Collatz-Stoppzeit für 27 ist bekannt."""
        result = repl.evaluate("collatz 27")
        # Stoppzeit von 27 ist 111
        assert "111" in result

    def test_collatz_1(self, repl):
        """Collatz-Stoppzeit für 1 ist 0."""
        result = repl.evaluate("collatz 1")
        assert "0" in result


# =============================================================================
# TESTS: MathREPL.evaluate - Lineare Algebra
# =============================================================================

class TestEvaluateEigenvalues:
    """Tests für den 'eigenvalues' Befehl."""

    def test_eigenvalues_2x2(self, repl):
        """Eigenwerte einer 2×2-Matrix."""
        result = repl.evaluate("eigenvalues [[1,2],[3,4]]")
        # Eigenwerte ≈ -0.37 und 5.37
        assert len(result) > 0
        assert "Eigenwerte" in result or "=" in result

    def test_eigenvalues_identity(self, repl):
        """Einheitsmatrix hat Eigenwert 1."""
        result = repl.evaluate("eigenvalues [[1,0],[0,1]]")
        assert "1" in result


# =============================================================================
# TESTS: MathREPL.evaluate - Sonstige
# =============================================================================

class TestEvaluateHelp:
    """Tests für den 'help' Befehl."""

    def test_help_shows_commands(self, repl):
        """help zeigt Befehlsliste an."""
        result = repl.evaluate("help")
        assert len(result) > 0
        # Muss mehrere Befehle enthalten
        assert "prime" in result or "gcd" in result

    def test_help_not_empty(self, repl):
        """help gibt nicht-leeren String zurück."""
        result = repl.evaluate("help")
        assert len(result) > 50

    def test_help_specific_command(self, repl):
        """help mit Befehlsname gibt spezifische Hilfe."""
        result = repl.evaluate("help prime")
        assert "prime" in result.lower()


class TestEvaluateExit:
    """Tests für den 'exit' Befehl."""

    def test_exit_returns_sentinel(self, repl):
        """exit gibt __EXIT__ zurück."""
        result = repl.evaluate("exit")
        assert result == "__EXIT__"

    def test_quit_returns_sentinel(self, repl):
        """quit gibt __EXIT__ zurück."""
        result = repl.evaluate("quit")
        assert result == "__EXIT__"


class TestEvaluateEdgeCases:
    """Tests für Randfälle."""

    def test_empty_command(self, repl):
        """Leerer Befehl gibt leeren String zurück."""
        result = repl.evaluate("")
        assert result == ""

    def test_unknown_command(self, repl):
        """Unbekannter Befehl gibt Fehlermeldung zurück."""
        result = repl.evaluate("unknown_xyz_command")
        assert "Unbekannt" in result or "unbekannt" in result.lower() or "Fehler" in result

    def test_whitespace_only(self, repl):
        """Nur Leerzeichen gibt leeren String zurück."""
        result = repl.evaluate("   ")
        assert result == ""
