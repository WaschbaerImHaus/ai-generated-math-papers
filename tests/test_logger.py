"""
@file test_logger.py
@brief Tests für das Logging-Modul (math_logger.py).
@description
    Überprüft alle Funktionen und Klassen des MathLogger-Moduls auf
    korrektes Verhalten, fehlende Ausnahmen und Level-Steuerung.

    Getestete Funktionalität:
    - MathLogger.__init__: Initialisierung mit verschiedenen Parametern
    - MathLogger.step: Berechnungsschritt protokollieren
    - MathLogger.result: Endergebnis protokollieren
    - MathLogger.convergence_warning: Konvergenz-Warnung
    - MathLogger.matrix_step: Matrix-Umformungsschritt
    - MathLogger.timing: Ausführungszeit protokollieren
    - MathLogger.section: Abschnittstrennlinie
    - MathLogger.set_level: Log-Level zur Laufzeit ändern
    - MathLogger.disable: Logger deaktivieren
    - get_logger: Globalen Logger abrufen
    - enable_debug_logging: Debug-Level aktivieren
    - disable_logging: Logging deaktivieren

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import logging
import pytest
import time
import io

# src-Verzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from math_logger import MathLogger, get_logger, enable_debug_logging, disable_logging


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def capture_logger_output(logger: MathLogger) -> io.StringIO:
    """
    @brief Fängt die Konsolenausgabe eines Loggers in einem StringIO-Puffer auf.
    @description
        Tauscht vorhandene StreamHandler gegen einen neuen aus, der in StringIO schreibt.
        Ermöglicht Tests, die den tatsächlichen Log-Output prüfen.

    @param logger: MathLogger-Instanz
    @return StringIO-Objekt mit aufgezeichneter Ausgabe
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    buffer = io.StringIO()
    # Alle Handler durch StringIO-Handler ersetzen
    logger._logger.handlers.clear()
    handler = logging.StreamHandler(buffer)
    handler.setFormatter(logging.Formatter("%(levelname)s:%(message)s"))
    logger._logger.addHandler(handler)
    return buffer


# ===========================================================================
# TESTS: MathLogger.__init__
# ===========================================================================

class TestMathLoggerInit:
    """Tests für MathLogger-Initialisierung."""

    def test_default_init_succeeds(self):
        """Standard-Initialisierung ohne Argumente wirft keine Exception."""
        logger = MathLogger()
        assert logger is not None

    def test_custom_name(self):
        """Initialisierung mit benutzerdefiniertem Namen."""
        logger = MathLogger(name="test-logger")
        assert logger.name == "test-logger"

    def test_debug_level(self):
        """Initialisierung mit DEBUG-Level."""
        logger = MathLogger(level="DEBUG")
        assert logger._logger.level == logging.DEBUG

    def test_info_level(self):
        """Initialisierung mit INFO-Level (Standard)."""
        logger = MathLogger(level="INFO")
        assert logger._logger.level == logging.INFO

    def test_warning_level(self):
        """Initialisierung mit WARNING-Level."""
        logger = MathLogger(level="WARNING")
        assert logger._logger.level == logging.WARNING

    def test_error_level(self):
        """Initialisierung mit ERROR-Level."""
        logger = MathLogger(level="ERROR")
        assert logger._logger.level == logging.ERROR

    def test_unknown_level_defaults_to_info(self):
        """Unbekannter Level → Standard INFO."""
        logger = MathLogger(level="UNKNOWN")
        assert logger._logger.level == logging.INFO

    def test_log_to_file_creates_file(self, tmp_path):
        """log_to_file=True erstellt eine Log-Datei im angegebenen Verzeichnis."""
        logger = MathLogger(name="file-test", log_to_file=True, log_dir=str(tmp_path))
        # Mindestens ein Handler muss ein FileHandler sein
        file_handlers = [h for h in logger._logger.handlers if isinstance(h, logging.FileHandler)]
        assert len(file_handlers) == 1

    def test_multiple_loggers_independent(self):
        """Zwei MathLogger-Instanzen sind voneinander unabhängig."""
        logger1 = MathLogger(name="l1", level="DEBUG")
        logger2 = MathLogger(name="l2", level="ERROR")
        assert logger1._logger.level != logger2._logger.level


# ===========================================================================
# TESTS: MathLogger.step
# ===========================================================================

class TestMathLoggerStep:
    """Tests für MathLogger.step()."""

    def test_step_does_not_raise(self):
        """step() wirft keine Exception."""
        logger = MathLogger(level="DEBUG")
        logger.step("Newton-Raphson", iteration=1, x=2.5, fx=0.25)  # Kein Fehler erwartet

    def test_step_with_empty_kwargs(self):
        """step() ohne kwargs funktioniert."""
        logger = MathLogger(level="DEBUG")
        logger.step("Test-Algorithmus")  # Kein Fehler erwartet

    def test_step_with_float_values(self):
        """step() mit float-Werten formatiert korrekt."""
        logger = MathLogger(level="DEBUG")
        buffer = capture_logger_output(logger)
        logger.step("Test", x=1.23456789, fx=0.00001234)
        output = buffer.getvalue()
        assert "Test" in output
        # Float sollte auf 4 Dezimalstellen formatiert sein
        assert "1.2346" in output

    def test_step_with_integer_values(self):
        """step() mit ganzzahligen Werten."""
        logger = MathLogger(level="DEBUG")
        buffer = capture_logger_output(logger)
        logger.step("Test", iteration=5)
        output = buffer.getvalue()
        assert "iteration=5" in output

    def test_step_appears_at_debug_level(self):
        """step() erscheint nur wenn Level DEBUG ist."""
        # Logger mit DEBUG
        logger_debug = MathLogger(level="DEBUG")
        buf_debug = capture_logger_output(logger_debug)
        logger_debug.step("TestAlgo", iteration=1, x=1.0)

        # Logger mit INFO (step soll nicht erscheinen)
        logger_info = MathLogger(level="INFO")
        buf_info = capture_logger_output(logger_info)
        logger_info.step("TestAlgo", iteration=1, x=1.0)

        assert len(buf_debug.getvalue()) > 0
        assert len(buf_info.getvalue()) == 0

    def test_step_with_string_value(self):
        """step() mit String-Wert."""
        logger = MathLogger(level="DEBUG")
        buffer = capture_logger_output(logger)
        logger.step("Test", phase="Vorwärts")
        output = buffer.getvalue()
        assert "Vorwärts" in output


# ===========================================================================
# TESTS: MathLogger.result
# ===========================================================================

class TestMathLoggerResult:
    """Tests für MathLogger.result()."""

    def test_result_does_not_raise(self):
        """result() wirft keine Exception."""
        logger = MathLogger()
        logger.result("Nullstelle", 3.0)  # Kein Fehler erwartet

    def test_result_with_unit(self):
        """result() mit Einheit funktioniert."""
        logger = MathLogger()
        logger.result("Geschwindigkeit", 9.81, "m/s²")  # Kein Fehler erwartet

    def test_result_with_integer(self):
        """result() mit ganzzahligem Wert."""
        logger = MathLogger()
        buffer = capture_logger_output(logger)
        logger.result("Primzahl", 17)
        output = buffer.getvalue()
        assert "17" in output

    def test_result_with_float(self):
        """result() formatiert float auf 6 Dezimalstellen."""
        logger = MathLogger()
        buffer = capture_logger_output(logger)
        logger.result("Wert", 1.41421356)
        output = buffer.getvalue()
        assert "1.414214" in output

    def test_result_appears_at_info_level(self):
        """result() erscheint bei INFO-Level."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.result("Test", 42.0)
        assert len(buffer.getvalue()) > 0

    def test_result_contains_description(self):
        """result() enthält die Beschreibung im Output."""
        logger = MathLogger()
        buffer = capture_logger_output(logger)
        logger.result("Meine-Nullstelle", 3.14)
        assert "Meine-Nullstelle" in buffer.getvalue()


# ===========================================================================
# TESTS: MathLogger.convergence_warning
# ===========================================================================

class TestMathLoggerConvergenceWarning:
    """Tests für MathLogger.convergence_warning()."""

    def test_convergence_warning_does_not_raise(self):
        """convergence_warning() wirft keine Exception."""
        logger = MathLogger()
        logger.convergence_warning("Newton-Raphson", 50, 0.001)

    def test_convergence_warning_appears_at_warning_level(self):
        """Warnung erscheint bei WARNING-Level."""
        logger = MathLogger(level="WARNING")
        buffer = capture_logger_output(logger)
        logger.convergence_warning("TestAlgo", 10, 0.05)
        output = buffer.getvalue()
        assert len(output) > 0
        assert "TestAlgo" in output

    def test_convergence_warning_contains_iteration(self):
        """Warnung enthält Iterationsnummer."""
        logger = MathLogger(level="WARNING")
        buffer = capture_logger_output(logger)
        logger.convergence_warning("Algo", 42, 0.001)
        assert "42" in buffer.getvalue()

    def test_convergence_warning_contains_residual(self):
        """Warnung enthält das Residuum."""
        logger = MathLogger(level="WARNING")
        buffer = capture_logger_output(logger)
        logger.convergence_warning("Algo", 5, 1.2345e-3)
        output = buffer.getvalue()
        # Residuum als wissenschaftliche Notation
        assert "1.2345e-03" in output or "1.23" in output


# ===========================================================================
# TESTS: MathLogger.matrix_step
# ===========================================================================

class TestMathLoggerMatrixStep:
    """Tests für MathLogger.matrix_step()."""

    def test_matrix_step_does_not_raise(self):
        """matrix_step() wirft keine Exception."""
        logger = MathLogger(level="DEBUG")
        before = [[2.0, 1.0], [1.0, 3.0]]
        after = [[2.0, 1.0], [0.0, 2.5]]
        logger.matrix_step("Zeile 2 -= 0.5 × Zeile 1", before, after)

    def test_matrix_step_with_description(self):
        """matrix_step() mit Beschreibung."""
        logger = MathLogger(level="DEBUG")
        logger.matrix_step(
            "Test-Op",
            [[1, 0], [0, 1]],
            [[2, 0], [0, 2]],
            description="Skalar-Multiplikation"
        )

    def test_matrix_step_appears_at_debug_level(self):
        """matrix_step() erscheint bei DEBUG-Level."""
        logger = MathLogger(level="DEBUG")
        buffer = capture_logger_output(logger)
        before = [[1.0, 0.0], [0.0, 1.0]]
        after = [[2.0, 0.0], [0.0, 2.0]]
        logger.matrix_step("Test", before, after)
        assert len(buffer.getvalue()) > 0

    def test_matrix_step_not_at_info_level(self):
        """matrix_step() erscheint NICHT bei INFO-Level (ist DEBUG)."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.matrix_step("Test", [[1.0]], [[2.0]])
        assert len(buffer.getvalue()) == 0


# ===========================================================================
# TESTS: MathLogger.timing
# ===========================================================================

class TestMathLoggerTiming:
    """Tests für MathLogger.timing()."""

    def test_timing_does_not_raise(self):
        """timing() wirft keine Exception."""
        logger = MathLogger()
        logger.timing("meine_funktion", 0.0123)

    def test_timing_appears_at_info_level(self):
        """timing() erscheint bei INFO-Level."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.timing("test_func", 0.5678)
        output = buffer.getvalue()
        assert len(output) > 0
        assert "test_func" in output

    def test_timing_contains_elapsed_time(self):
        """timing() enthält die Zeitangabe."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.timing("func", 1.2345)
        assert "1.2345" in buffer.getvalue()


# ===========================================================================
# TESTS: MathLogger.section
# ===========================================================================

class TestMathLoggerSection:
    """Tests für MathLogger.section()."""

    def test_section_does_not_raise(self):
        """section() wirft keine Exception."""
        logger = MathLogger()
        logger.section("Mein Abschnitt")

    def test_section_appears_at_info_level(self):
        """section() erscheint bei INFO-Level."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.section("Test Abschnitt")
        output = buffer.getvalue()
        assert len(output) > 0

    def test_section_contains_title(self):
        """section() enthält den Titel (großgeschrieben)."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.section("Newton-Raphson")
        output = buffer.getvalue().upper()
        assert "NEWTON-RAPHSON" in output

    def test_section_contains_separator(self):
        """section() enthält Trennzeichen."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)
        logger.section("Test")
        output = buffer.getvalue()
        assert "=" in output


# ===========================================================================
# TESTS: MathLogger.set_level / disable
# ===========================================================================

class TestMathLoggerLevelControl:
    """Tests für Level-Steuerung."""

    def test_set_level_to_debug(self):
        """set_level('DEBUG') wechselt zu DEBUG."""
        logger = MathLogger(level="INFO")
        logger.set_level("DEBUG")
        assert logger._logger.level == logging.DEBUG

    def test_set_level_to_warning(self):
        """set_level('WARNING') wechselt zu WARNING."""
        logger = MathLogger(level="INFO")
        logger.set_level("WARNING")
        assert logger._logger.level == logging.WARNING

    def test_set_level_enables_debug_output(self):
        """Nach set_level('DEBUG') erscheinen Debug-Meldungen."""
        logger = MathLogger(level="INFO")
        buffer = capture_logger_output(logger)

        # Vor dem Wechsel: step() nicht sichtbar
        logger.step("Test1", x=1.0)
        assert "Test1" not in buffer.getvalue()

        # Nach dem Wechsel: step() sichtbar
        logger.set_level("DEBUG")
        logger.step("Test2", x=2.0)
        assert "Test2" in buffer.getvalue()

    def test_disable_suppresses_all_output(self):
        """disable() unterdrückt alle Ausgaben."""
        logger = MathLogger(level="DEBUG")
        buffer = capture_logger_output(logger)

        logger.disable()
        logger.step("Geheim", x=999.0)
        logger.result("Geheim", 999.0)
        logger.section("Geheimabschnitt")

        assert len(buffer.getvalue()) == 0


# ===========================================================================
# TESTS: Globale Hilfsfunktionen
# ===========================================================================

class TestGlobalFunctions:
    """Tests für get_logger(), enable_debug_logging(), disable_logging()."""

    def test_get_logger_without_name_returns_default(self):
        """get_logger() ohne Argument gibt globalen Logger zurück."""
        logger = get_logger()
        assert isinstance(logger, MathLogger)

    def test_get_logger_with_name_returns_new_logger(self):
        """get_logger(name) gibt neuen Logger mit diesem Namen zurück."""
        logger = get_logger("mein-logger")
        assert isinstance(logger, MathLogger)
        assert logger.name == "mein-logger"

    def test_get_logger_different_instances(self):
        """Verschiedene Namen → verschiedene Instanzen."""
        l1 = get_logger("a")
        l2 = get_logger("b")
        assert l1 is not l2

    def test_enable_debug_logging_does_not_raise(self):
        """enable_debug_logging() wirft keine Exception."""
        enable_debug_logging()  # Kein Fehler erwartet

    def test_disable_logging_does_not_raise(self):
        """disable_logging() wirft keine Exception."""
        disable_logging()  # Kein Fehler erwartet

    def test_enable_debug_sets_debug_level(self):
        """enable_debug_logging() setzt globalen Logger auf DEBUG."""
        enable_debug_logging()
        from math_logger import _default_logger
        assert _default_logger._logger.level == logging.DEBUG

    def test_disable_suppresses_global_output(self):
        """disable_logging() unterdrückt Ausgaben des globalen Loggers."""
        # Globalen Logger holen und Ausgabe zurücksetzen
        from math_logger import _default_logger
        buffer = capture_logger_output(_default_logger)
        disable_logging()
        _default_logger.result("Soll nicht erscheinen", 0.0)
        assert len(buffer.getvalue()) == 0


# ===========================================================================
# EDGE-CASE-TESTS
# ===========================================================================

class TestEdgeCases:
    """Edge-Cases und Sonderfälle."""

    def test_step_with_very_large_float(self):
        """step() mit sehr großem float-Wert."""
        logger = MathLogger(level="DEBUG")
        logger.step("Test", x=1e300)  # Kein Fehler erwartet

    def test_step_with_negative_float(self):
        """step() mit negativem float-Wert."""
        logger = MathLogger(level="DEBUG")
        logger.step("Test", x=-3.14159)  # Kein Fehler erwartet

    def test_result_with_none_value(self):
        """result() mit None-Wert."""
        logger = MathLogger()
        logger.result("Test", None)  # Kein Fehler erwartet

    def test_result_with_string_value(self):
        """result() mit String-Wert."""
        logger = MathLogger()
        logger.result("Status", "Konvergiert")  # Kein Fehler erwartet

    def test_result_with_list_value(self):
        """result() mit Listen-Wert."""
        logger = MathLogger()
        logger.result("Eigenwerte", [1.0, 2.0, 3.0])  # Kein Fehler erwartet

    def test_matrix_step_with_empty_matrix(self):
        """matrix_step() mit leerer Matrix."""
        logger = MathLogger(level="DEBUG")
        logger.matrix_step("Test", [], [])  # Kein Fehler erwartet

    def test_section_with_empty_title(self):
        """section() mit leerem Titel."""
        logger = MathLogger()
        logger.section("")  # Kein Fehler erwartet

    def test_multiple_log_calls_sequence(self):
        """Mehrere verschiedene Log-Methoden nacheinander."""
        logger = MathLogger(level="DEBUG")
        # Alle Methoden nacheinander – keine Exception darf auftreten
        logger.section("Start")
        logger.step("Algo", iteration=1, x=1.0, fx=0.5)
        logger.convergence_warning("Algo", 10, 0.01)
        logger.result("Ergebnis", 1.4142)
        logger.timing("algo_func", 0.001)
        logger.matrix_step("Op", [[1.0]], [[2.0]])

    def test_log_to_file_with_nonexistent_dir(self, tmp_path):
        """log_to_file mit noch nicht existierendem Unterverzeichnis."""
        new_dir = str(tmp_path / "deep" / "nested" / "dir")
        # Sollte das Verzeichnis erstellen und nicht abstürzen
        logger = MathLogger(name="test", log_to_file=True, log_dir=new_dir)
        assert logger is not None

    def test_timing_with_zero_elapsed(self):
        """timing() mit 0.0 Sekunden."""
        logger = MathLogger()
        logger.timing("instant_func", 0.0)  # Kein Fehler erwartet

    def test_convergence_warning_with_zero_residual(self):
        """convergence_warning() mit Residuum = 0."""
        logger = MathLogger()
        logger.convergence_warning("Test", 1, 0.0)  # Kein Fehler erwartet
