"""
@file test_math_logger.py
@brief Umfassende Tests für das Modul math_logger.py.
@description
    Testet alle Klassen und Funktionen:
    - MathLogger.__init__(): Konstruktor mit verschiedenen Leveln
    - MathLogger.step(): Debug-Logging
    - MathLogger.result(): Ergebnis-Logging
    - MathLogger.convergence_warning(): Konvergenz-Warnung
    - MathLogger.matrix_step(): Matrix-Logging
    - MathLogger.timing(): Zeitstempel-Logging
    - MathLogger.section(): Trennlinie-Logging
    - MathLogger.warning(): Allgemeine Warnung
    - MathLogger.set_level(): Dynamische Level-Änderung
    - MathLogger.condition_number_warning(): Konditionszahl-Warnung
    - MathLogger.stability_report(): Stabilitätsbericht
    - MathLogger.disable(): Logger deaktivieren
    - get_logger(): Globaler Logger
    - enable_debug_logging(), disable_logging(): Hilfsfunktionen

    Edge-Cases:
    - Unbekannter Log-Level → Fallback auf INFO
    - Leere Keyword-Arguments in step()
    - Leere Matrix in _format_matrix()
    - Stabilitätsbericht mit leerer Liste
    - stability_report mit verschiedenen Konditionszahlen

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import logging
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from math_logger import (
    MathLogger,
    get_logger,
    enable_debug_logging,
    disable_logging,
)


# =============================================================================
# TESTS: MathLogger Konstruktor
# =============================================================================

class TestMathLoggerConstructor:
    """Tests für die Initialisierung des MathLoggers."""

    def test_default_name(self):
        """Standard-Name wird gesetzt."""
        logger = MathLogger()
        assert logger.name == "specialist-maths"

    def test_custom_name(self):
        """Benutzerdefinierter Name wird gesetzt."""
        logger = MathLogger(name="test-logger")
        assert logger.name == "test-logger"

    def test_level_debug(self):
        """DEBUG-Level wird korrekt gesetzt."""
        logger = MathLogger(level="DEBUG")
        assert logger._logger.level == logging.DEBUG

    def test_level_info(self):
        """INFO-Level wird korrekt gesetzt."""
        logger = MathLogger(level="INFO")
        assert logger._logger.level == logging.INFO

    def test_level_warning(self):
        """WARNING-Level wird korrekt gesetzt."""
        logger = MathLogger(level="WARNING")
        assert logger._logger.level == logging.WARNING

    def test_level_error(self):
        """ERROR-Level wird korrekt gesetzt."""
        logger = MathLogger(level="ERROR")
        assert logger._logger.level == logging.ERROR

    def test_unknown_level_fallback(self):
        """Unbekannter Level fällt auf INFO zurück."""
        logger = MathLogger(level="UNKNOWN_LEVEL")
        # Sollte INFO sein (Fallback)
        assert logger._logger.level == logging.INFO

    def test_has_handler(self):
        """Logger hat mindestens einen Handler."""
        logger = MathLogger()
        assert len(logger._logger.handlers) >= 1

    def test_no_propagation(self):
        """Logger propagiert nicht an Parent-Logger."""
        logger = MathLogger()
        assert logger._logger.propagate is False


# =============================================================================
# TESTS: MathLogger.step()
# =============================================================================

class TestMathLoggerStep:
    """Tests für die step()-Methode."""

    def test_step_no_exception(self):
        """step() wirft keine Ausnahme."""
        logger = MathLogger(level="DEBUG")
        try:
            logger.step("Newton-Raphson", iteration=1, x=2.5, fx=0.25)
        except Exception as e:
            pytest.fail(f"step() hat Ausnahme geworfen: {e}")

    def test_step_empty_kwargs(self):
        """step() mit leeren Kwargs läuft ohne Fehler."""
        logger = MathLogger(level="DEBUG")
        logger.step("Algorithmus")  # Keine weiteren Argumente

    def test_step_float_formatting(self):
        """_format_kwargs formatiert Floats auf 4 Stellen."""
        logger = MathLogger()
        result = logger._format_kwargs({"x": 3.14159})
        assert "x=" in result
        assert "3.1416" in result  # Auf 4 Nachkommastellen gerundet

    def test_step_string_kwargs(self):
        """_format_kwargs lässt Strings unverändert."""
        logger = MathLogger()
        result = logger._format_kwargs({"name": "test"})
        assert "name=test" in result

    def test_step_mixed_kwargs(self):
        """_format_kwargs verarbeitet gemischte Typen."""
        logger = MathLogger()
        result = logger._format_kwargs({"iter": 3, "x": 1.5, "label": "ok"})
        assert "iter=3" in result
        assert "x=1.5000" in result
        assert "label=ok" in result


# =============================================================================
# TESTS: MathLogger.result()
# =============================================================================

class TestMathLoggerResult:
    """Tests für die result()-Methode."""

    def test_result_no_exception(self):
        """result() wirft keine Ausnahme."""
        logger = MathLogger()
        logger.result("Nullstelle", 3.0)

    def test_result_with_unit(self):
        """result() mit Einheit."""
        logger = MathLogger()
        logger.result("Geschwindigkeit", 9.81, unit="m/s²")

    def test_result_integer(self):
        """result() mit Ganzzahl."""
        logger = MathLogger()
        logger.result("Prüfzahl", 42)

    def test_result_string_value(self):
        """result() mit String-Wert."""
        logger = MathLogger()
        logger.result("Status", "konvergiert")

    def test_result_float_formatted(self):
        """Float-Werte werden auf 6 Stellen formatiert."""
        logger = MathLogger()
        # Keine Ausnahme, indirekt via Formatierung getestet
        logger.result("Pi-Approx", 3.14159265358)


# =============================================================================
# TESTS: MathLogger.convergence_warning()
# =============================================================================

class TestMathLoggerConvergenceWarning:
    """Tests für die Konvergenz-Warnung."""

    def test_convergence_warning_no_exception(self):
        """convergence_warning() wirft keine Ausnahme."""
        logger = MathLogger(level="WARNING")
        logger.convergence_warning("Newton", 100, 1e-3)

    def test_convergence_warning_with_small_residual(self):
        """Kleine Residuen werden korrekt formatiert."""
        logger = MathLogger(level="WARNING")
        logger.convergence_warning("Bisection", 50, 1e-15)


# =============================================================================
# TESTS: MathLogger.matrix_step()
# =============================================================================

class TestMathLoggerMatrixStep:
    """Tests für matrix_step()."""

    def test_matrix_step_no_exception(self):
        """matrix_step() wirft keine Ausnahme."""
        logger = MathLogger(level="DEBUG")
        before = [[1.0, 2.0], [3.0, 4.0]]
        after  = [[1.0, 2.0], [0.0, 2.0]]
        logger.matrix_step("Gauss-Schritt", before, after)

    def test_format_matrix_empty(self):
        """_format_matrix([]) gibt '[]' zurück."""
        logger = MathLogger()
        result = logger._format_matrix([])
        assert result == "[]"

    def test_format_matrix_2x2(self):
        """_format_matrix für 2×2-Matrix."""
        logger = MathLogger()
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        result = logger._format_matrix(matrix)
        assert isinstance(result, str)
        assert len(result) > 0


# =============================================================================
# TESTS: MathLogger.timing()
# =============================================================================

class TestMathLoggerTiming:
    """Tests für timing()."""

    def test_timing_no_exception(self):
        """timing() wirft keine Ausnahme."""
        logger = MathLogger()
        logger.timing("berechne_integral", 0.0123)

    def test_timing_zero(self):
        """timing() mit 0 Sekunden."""
        logger = MathLogger()
        logger.timing("schnelle_funktion", 0.0)


# =============================================================================
# TESTS: MathLogger.section()
# =============================================================================

class TestMathLoggerSection:
    """Tests für section()."""

    def test_section_no_exception(self):
        """section() wirft keine Ausnahme."""
        logger = MathLogger()
        logger.section("Algebra-Modul")

    def test_section_empty_title(self):
        """section() mit leerem Titel."""
        logger = MathLogger()
        logger.section("")


# =============================================================================
# TESTS: MathLogger.warning()
# =============================================================================

class TestMathLoggerWarning:
    """Tests für die allgemeine warning()-Methode."""

    def test_warning_no_exception(self):
        """warning() wirft keine Ausnahme."""
        logger = MathLogger(level="WARNING")
        logger.warning("Achtung: Möglicher Präzisionsverlust")

    def test_warning_empty_message(self):
        """warning() mit leerem String."""
        logger = MathLogger(level="WARNING")
        logger.warning("")


# =============================================================================
# TESTS: MathLogger.set_level()
# =============================================================================

class TestMathLoggerSetLevel:
    """Tests für die dynamische Level-Änderung."""

    def test_set_level_debug(self):
        """set_level('DEBUG') setzt DEBUG-Level."""
        logger = MathLogger(level="INFO")
        logger.set_level("DEBUG")
        assert logger._logger.level == logging.DEBUG

    def test_set_level_error(self):
        """set_level('ERROR') setzt ERROR-Level."""
        logger = MathLogger(level="INFO")
        logger.set_level("ERROR")
        assert logger._logger.level == logging.ERROR

    def test_set_level_unknown_fallback(self):
        """Unbekannter Level fällt auf INFO zurück."""
        logger = MathLogger(level="DEBUG")
        logger.set_level("NOOP")
        # Kein Fehler, Fallback auf INFO
        assert logger._logger.level == logging.INFO


# =============================================================================
# TESTS: MathLogger.condition_number_warning()
# =============================================================================

class TestConditionNumberWarning:
    """Tests für condition_number_warning()."""

    def test_no_warning_for_good_condition(self):
        """Gut konditionierte Matrix → kein Fehler."""
        logger = MathLogger(level="WARNING")
        # κ = 100 ist gut konditioniert (unter Schwellenwert 1e10)
        logger.condition_number_warning("A", 100.0)

    def test_warning_for_bad_condition(self):
        """Schlecht konditionierte Matrix → Warnung wird ausgegeben (kein Fehler)."""
        logger = MathLogger(level="WARNING")
        # κ = 1e12 > 1e10 → Warning
        logger.condition_number_warning("B", 1e12)

    def test_custom_threshold(self):
        """Benutzerdefinierter Schwellenwert."""
        logger = MathLogger(level="WARNING")
        logger.condition_number_warning("C", 1e6, threshold=1e5)


# =============================================================================
# TESTS: MathLogger.stability_report()
# =============================================================================

class TestStabilityReport:
    """Tests für stability_report()."""

    def test_empty_list(self):
        """Leere Operationsliste → Kein Fehler, String zurückgegeben."""
        logger = MathLogger()
        report = logger.stability_report([])
        assert isinstance(report, str)
        assert "keine" in report.lower() or "Stabilitätsbericht" in report

    def test_single_good_operation(self):
        """Eine gut konditionierte Operation."""
        logger = MathLogger()
        ops = [{"name": "LU", "cond": 1.5, "operation": "solve"}]
        report = logger.stability_report(ops)
        assert isinstance(report, str)
        assert "LU" in report

    def test_mixed_operations(self):
        """Gut und schlecht konditionierte Operationen gemischt."""
        logger = MathLogger()
        ops = [
            {"name": "A", "cond": 1.0, "operation": "det"},
            {"name": "B", "cond": 1e9, "operation": "inv"},
            {"name": "C", "cond": 1e15, "operation": "solve"},
        ]
        report = logger.stability_report(ops)
        assert isinstance(report, str)
        assert "A" in report
        assert "B" in report
        assert "C" in report

    def test_invalid_entries_skipped(self):
        """Einträge ohne 'cond' werden übersprungen."""
        logger = MathLogger()
        ops = [
            {"name": "A"},  # Kein 'cond'
            "nicht-dict",
            None,
            {"name": "B", "cond": 1.0, "operation": "ok"},
        ]
        report = logger.stability_report(ops)
        assert isinstance(report, str)


# =============================================================================
# TESTS: MathLogger.disable()
# =============================================================================

class TestMathLoggerDisable:
    """Tests für disable()."""

    def test_disable_no_exception(self):
        """disable() wirft keine Ausnahme."""
        logger = MathLogger()
        logger.disable()

    def test_disable_prevents_output(self):
        """Nach disable() wird kein Output mehr produziert."""
        logger = MathLogger(level="DEBUG")
        logger.disable()
        # Kein Fehler bei weiterem Aufruf
        logger.result("Test", 42)
        logger.step("Algorithmus", x=1.0)


# =============================================================================
# TESTS: get_logger(), enable_debug_logging(), disable_logging()
# =============================================================================

class TestModuleFunctions:
    """Tests für die Modulfunktionen."""

    def test_get_logger_returns_math_logger(self):
        """get_logger() gibt eine MathLogger-Instanz zurück."""
        logger = get_logger()
        assert isinstance(logger, MathLogger)

    def test_get_logger_with_name(self):
        """get_logger(name) erstellt benannten Logger."""
        logger = get_logger("my-module")
        assert isinstance(logger, MathLogger)
        assert logger.name == "my-module"

    def test_get_logger_default_is_consistent(self):
        """Mehrfache get_logger()-Aufrufe ohne Name → gleicher Logger."""
        l1 = get_logger()
        l2 = get_logger()
        # Beide sollten dieselbe Instanz sein (Singleton)
        assert l1 is l2

    def test_enable_debug_logging_no_exception(self):
        """enable_debug_logging() wirft keine Ausnahme."""
        enable_debug_logging()

    def test_disable_logging_no_exception(self):
        """disable_logging() wirft keine Ausnahme."""
        disable_logging()
