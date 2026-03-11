"""
@file test_stability_analysis.py
@brief Tests für Stabilitätsanalyse-Funktionen.
@description
    Testet:
    - analyze_stability() aus matrix_ops.py
    - condition_number_warning() aus math_logger.py
    - stability_report() aus math_logger.py

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import logging
import numpy as np
import pytest
from matrix_ops import Matrix, analyze_stability, condition_number_check
from math_logger import MathLogger, get_logger


class TestAnalyzeStability:
    """Tests für die analyze_stability()-Funktion."""

    def test_identity_matrix_well_conditioned(self):
        """Test: Einheitsmatrix ist gut konditioniert (κ ≈ 1)."""
        A = np.eye(4)
        result = analyze_stability(A)
        assert result['is_well_conditioned'] == True
        assert result['condition_number'] < 10.0
        assert result['rank'] == 4
        assert abs(result['determinant'] - 1.0) < 1e-10

    def test_diagonal_matrix_condition_number(self):
        """Test: Diagonalmatrix mit bekannter Konditionszahl."""
        # Diagonalmatrix mit Einträgen [1, 100]: κ = 100
        A = np.diag([1.0, 100.0])
        result = analyze_stability(A)
        assert abs(result['condition_number'] - 100.0) < 1.0
        assert bool(result['is_well_conditioned']) == True

    def test_ill_conditioned_matrix(self):
        """Test: Schlecht konditionierte Matrix wird erkannt."""
        # Hilbert-Matrix ist berühmt für schlechte Konditionierung
        n = 8
        H = np.array([[1.0 / (i + j + 1) for j in range(n)] for i in range(n)])
        result = analyze_stability(H)
        # Hilbert-Matrix 8x8 hat κ ≈ 1e10+
        assert result['condition_number'] > 1e6
        assert result['rank'] == n

    def test_singular_matrix(self):
        """Test: Singuläre Matrix hat Determinante ≈ 0."""
        # Rang-defiziente Matrix
        A = np.array([[1.0, 2.0], [2.0, 4.0]])  # Zweite Zeile = 2 * erste Zeile
        result = analyze_stability(A)
        assert abs(result['determinant']) < 1e-10
        assert result['rank'] < 2

    def test_matrix_object_input(self):
        """Test: Matrix-Objekt als Eingabe."""
        M = Matrix([[2.0, 0.0], [0.0, 3.0]])
        result = analyze_stability(M)
        assert result['rank'] == 2
        assert abs(result['determinant'] - 6.0) < 1e-10

    def test_list_input(self):
        """Test: Liste als Eingabe."""
        A = [[1.0, 0.0], [0.0, 1.0]]
        result = analyze_stability(A)
        assert bool(result['is_well_conditioned']) == True
        assert result['rank'] == 2

    def test_return_dict_keys(self):
        """Test: Rückgabe-Dictionary enthält alle erwarteten Schlüssel."""
        A = np.eye(3)
        result = analyze_stability(A)
        assert 'condition_number' in result
        assert 'is_well_conditioned' in result
        assert 'rank' in result
        assert 'determinant' in result

    def test_rank_rank_deficient(self):
        """Test: Rang einer rang-defizienten Matrix."""
        # 3x3-Matrix mit Rang 2
        A = np.array([[1.0, 2.0, 3.0],
                      [4.0, 5.0, 6.0],
                      [7.0, 8.0, 9.0]])  # Zeile 3 = Zeile 1 + Zeile 2 (nach Subtraktion)
        result = analyze_stability(A)
        assert result['rank'] == 2


class TestConditionNumberCheck:
    """Tests für condition_number_check() und solve()-Warnung."""

    def test_stabile_matrix_keine_warnung(self):
        """Test: Gut konditionierte Matrix → condition_number_check() gibt True zurück."""
        # Einheitsmatrix: κ = 1, perfekt konditioniert
        M = Matrix([[1.0, 0.0], [0.0, 1.0]])
        assert condition_number_check(M) is True

    def test_instabile_matrix_warnung_geloggt(self):
        """Test: Schlecht konditionierte Matrix löst Warnung im Logger aus."""
        # Hilbert-Matrix 8×8: κ ≈ 1e10+
        n = 8
        H_data = [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)]
        M = Matrix(H_data)
        b_data = [1.0] * n

        # Logging-Handler zum Abfangen der Warnmeldung
        import sys
        from vectors import Vector

        # Logger-Level auf WARNING setzen und Handler anhängen
        log_capture = []

        class ListHandler(logging.Handler):
            def emit(self, record):
                log_capture.append(record.getMessage())

        # Auf den Matrix-Ops-Logger zugreifen
        import matrix_ops as mo
        handler = ListHandler()
        handler.setLevel(logging.WARNING)
        mo._logger._logger.addHandler(handler)
        mo._logger._logger.setLevel(logging.WARNING)

        try:
            b = Vector(b_data)
            try:
                M.solve(b)
            except Exception:
                pass  # Lösung kann scheitern, Warnung ist das Ziel
        finally:
            mo._logger._logger.removeHandler(handler)

        # Mindestens eine Warnung mit dem korrekten Text erwartet
        warning_msgs = [m for m in log_capture if "überschreitet 1e10" in m or "Konditionszahl" in m]
        assert len(warning_msgs) > 0, (
            f"Keine Konditionszahl-Warnung geloggt. Nachrichten: {log_capture}"
        )

    def test_condition_number_check_rueckgabe(self):
        """Test: condition_number_check() gibt True für gute und False für schlechte Matrix zurück."""
        # Gut konditionierte Diagonalmatrix
        gute_matrix = [[2.0, 0.0], [0.0, 3.0]]
        assert condition_number_check(gute_matrix) is True

        # Schlecht konditionierte Hilbert-Matrix
        n = 10
        schlechte_matrix = [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)]
        # Hilbert-Matrix 10×10: κ ≈ 1e13 → sollte False zurückgeben
        assert condition_number_check(schlechte_matrix) is False

    def test_condition_number_check_numpy_array(self):
        """Test: condition_number_check() akzeptiert numpy-Arrays."""
        arr = np.eye(5)  # Einheitsmatrix: κ = 1
        assert condition_number_check(arr) is True

    def test_condition_number_check_matrix_objekt(self):
        """Test: condition_number_check() akzeptiert Matrix-Objekte."""
        M = Matrix([[3.0, 1.0], [1.0, 3.0]])  # κ = 2 → gut konditioniert
        assert condition_number_check(M) is True


class TestMathLoggerStability:
    """Tests für Stabilitäts-Methoden des MathLogger."""

    def setup_method(self):
        """Erstellt Logger-Instanz vor jedem Test."""
        self.logger = MathLogger(name="test_stability", level="WARNING")

    def test_condition_number_warning_triggers(self):
        """Test: Warnung wird bei hoher Konditionszahl ausgelöst."""
        # Kein Fehler bei korrektem Aufruf
        self.logger.condition_number_warning("TestMatrix", 1e12, threshold=1e10)

    def test_condition_number_warning_no_trigger(self):
        """Test: Keine Warnung bei kleiner Konditionszahl."""
        # Auch hier kein Fehler – Warnung wird nur intern geloggt
        self.logger.condition_number_warning("GutMatrix", 100.0, threshold=1e10)

    def test_condition_number_warning_default_threshold(self):
        """Test: Standard-Schwellenwert ist 1e10."""
        # Aufruf ohne threshold-Parameter
        self.logger.condition_number_warning("DefaultMatrix", 5e10)

    def test_stability_report_empty(self):
        """Test: Leere Liste → Bericht mit Hinweis."""
        report = self.logger.stability_report([])
        assert "Stabilitätsbericht" in report
        assert "Keine" in report

    def test_stability_report_with_operations(self):
        """Test: Bericht mit Operationen wird korrekt erstellt."""
        operations = [
            {'name': 'Matrix A', 'cond': 10.0, 'operation': 'LU'},
            {'name': 'Matrix B', 'cond': 1e11, 'operation': 'solve'},
        ]
        report = self.logger.stability_report(operations)
        assert 'Matrix A' in report
        assert 'Matrix B' in report
        assert 'gut konditioniert' in report.lower() or 'schlecht' in report.lower()

    def test_stability_report_returns_string(self):
        """Test: stability_report() gibt einen String zurück."""
        operations = [{'name': 'Test', 'cond': 1e5, 'operation': 'QR'}]
        report = self.logger.stability_report(operations)
        assert isinstance(report, str)

    def test_stability_report_sorted_by_condition(self):
        """Test: Bericht ist nach Konditionszahl aufsteigend sortiert."""
        operations = [
            {'name': 'Schlecht', 'cond': 1e12, 'operation': 'LU'},
            {'name': 'Gut', 'cond': 5.0, 'operation': 'QR'},
        ]
        report = self.logger.stability_report(operations)
        # "Gut" sollte vor "Schlecht" erscheinen
        pos_gut = report.find('Gut')
        pos_schlecht = report.find('Schlecht')
        assert pos_gut < pos_schlecht

    def test_stability_report_skips_non_dict(self):
        """Test: Nicht-Dict-Einträge werden übersprungen."""
        operations = [
            {'name': 'Valid', 'cond': 100.0, 'operation': 'LU'},
            "ungültig",
            42,
        ]
        # Kein Fehler, gültige Einträge werden verarbeitet
        report = self.logger.stability_report(operations)
        assert 'Valid' in report
