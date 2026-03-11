"""
@file test_benchmark_regression.py
@brief Benchmark-Regressions-Tests zur Laufzeitüberwachung.
@description
    Misst die Ausführungszeit kritischer Berechnungen und vergleicht
    sie mit einem gespeicherten Baseline-Wert. Bei Verschlechterung
    um mehr als 20% wird eine Warnung ausgegeben (kein Test-Fehler).

    Getestete Operationen:
    - prime_factorization(10**12): Primfaktorzerlegung großer Zahlen
    - FFT von 1024 Punkten: Fourier-Transformation
    - 100×100 Matrix-LU: LU-Zerlegung großer Matrizen

    Baseline wird in debugging/benchmark_baseline.json gespeichert.

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import time
import json
import warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pytest

# Pfad zur Baseline-Datei
BASELINE_FILE = os.path.join(
    os.path.dirname(__file__), '..', 'debugging', 'benchmark_baseline.json'
)

# Maximale Toleranz: 20% Verschlechterung → Warnung
REGRESSION_THRESHOLD = 1.20

# Harte Obergrenzen (Test schlägt fehl wenn überschritten)
HARD_LIMITS = {
    'prime_factorization_1e12': 1.0,    # max 1 Sekunde
    'fft_1024': 0.1,                     # max 0.1 Sekunde
    'matrix_lu_100x100': 0.5,            # max 0.5 Sekunde
}


def load_baseline() -> dict:
    """
    @brief Lädt die gespeicherte Baseline aus JSON-Datei.
    @return: Dictionary mit Baseline-Werten oder leeres Dict wenn nicht vorhanden.
    """
    if os.path.exists(BASELINE_FILE):
        try:
            with open(BASELINE_FILE, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return {}
    return {}


def save_baseline(new_values: dict) -> None:
    """
    @brief Speichert neue Baseline-Werte in JSON-Datei.
    @description
        Aktualisiert vorhandene Werte und schreibt die Datei.
    @param new_values: Neue Messwerte als Dictionary.
    """
    baseline = load_baseline()
    baseline.update(new_values)
    try:
        with open(BASELINE_FILE, 'w') as f:
            json.dump(baseline, f, indent=2)
    except IOError:
        pass  # Schreibfehler ignorieren in Tests


def check_regression(name: str, measured: float, baseline: dict) -> None:
    """
    @brief Prüft ob eine Messung stark von der Baseline abweicht.
    @description
        Warnt (via warnings.warn) wenn die gemessene Zeit um mehr als
        REGRESSION_THRESHOLD über der Baseline liegt.

    @param name: Name der Messung.
    @param measured: Gemessene Zeit in Sekunden.
    @param baseline: Dictionary mit Baseline-Werten.
    """
    if name in baseline:
        baseline_val = baseline[name]
        if measured > baseline_val * REGRESSION_THRESHOLD:
            pct = (measured / baseline_val - 1) * 100
            warnings.warn(
                f"Leistungsregression in '{name}': "
                f"{measured:.4f}s vs Baseline {baseline_val:.4f}s "
                f"(+{pct:.1f}%)",
                UserWarning,
                stacklevel=2
            )


class TestBenchmarkPrimeFactorization:
    """Benchmark-Test für Primfaktorzerlegung."""

    def test_prime_factorization_1e12_speed(self):
        """Test: prime_factorization(10**12) ist schneller als 1 Sekunde."""
        from algebra_numbertheory import prime_factorization

        baseline = load_baseline()

        # Messung
        start = time.perf_counter()
        result = prime_factorization(10**12)
        elapsed = time.perf_counter() - start

        # Ergebnis prüfen
        assert isinstance(result, dict), "Ergebnis muss ein Dictionary sein"
        assert len(result) > 0, "Muss mindestens einen Primfaktor haben"

        # 10^12 = 2^12 * 5^12
        assert 2 in result
        assert 5 in result

        # Hartes Limit prüfen
        assert elapsed < HARD_LIMITS['prime_factorization_1e12'], (
            f"prime_factorization(10**12) dauerte {elapsed:.3f}s "
            f"(Limit: {HARD_LIMITS['prime_factorization_1e12']}s)"
        )

        # Regression prüfen (nur Warnung)
        check_regression('prime_factorization_1e12', elapsed, baseline)

        # Baseline aktualisieren
        save_baseline({'prime_factorization_1e12': elapsed})


class TestBenchmarkFFT:
    """Benchmark-Test für FFT."""

    def test_fft_1024_speed(self):
        """Test: FFT von 1024 Punkten ist schneller als 0.1 Sekunden."""
        baseline = load_baseline()

        # Test-Signal: Sinussumme
        n = 1024
        t = np.linspace(0, 1, n)
        signal = np.sin(2 * np.pi * 10 * t) + 0.5 * np.sin(2 * np.pi * 50 * t)

        # Messung
        start = time.perf_counter()
        result = np.fft.fft(signal)
        elapsed = time.perf_counter() - start

        # Ergebnis prüfen
        assert len(result) == n, "FFT muss gleich viele Punkte wie Eingabe haben"
        assert result.dtype == complex, "FFT-Ergebnis muss komplex sein"

        # Hartes Limit
        assert elapsed < HARD_LIMITS['fft_1024'], (
            f"FFT({n}) dauerte {elapsed:.4f}s "
            f"(Limit: {HARD_LIMITS['fft_1024']}s)"
        )

        # Regression prüfen
        check_regression('fft_1024', elapsed, baseline)
        save_baseline({'fft_1024': elapsed})

    def test_fft_result_correctness(self):
        """Test: FFT-Ergebnis ist mathematisch korrekt."""
        # Einfacher Test: FFT eines reinen Kosinus
        n = 64
        freq = 4  # 4 Hz
        t = np.linspace(0, 1, n, endpoint=False)
        signal = np.cos(2 * np.pi * freq * t)

        fft_result = np.fft.fft(signal)
        magnitudes = np.abs(fft_result)

        # Peak sollte bei Frequenz 4 (Index 4) liegen
        peak_index = np.argmax(magnitudes[:n // 2])
        assert peak_index == freq, f"Peak bei {peak_index}, erwartet {freq}"


class TestBenchmarkMatrixLU:
    """Benchmark-Test für LU-Zerlegung."""

    def test_matrix_lu_100x100_speed(self):
        """Test: LU-Zerlegung einer 100×100-Matrix unter 0.5 Sekunden."""
        baseline = load_baseline()

        # Zufällige 100×100 Matrix (gut konditioniert)
        np.random.seed(42)  # Reproduzierbar
        A = np.random.randn(100, 100)
        A = A + 100 * np.eye(100)  # Diagonal-Dominanz für Stabilität

        # Messung
        start = time.perf_counter()
        P, L, U = _lu_decomposition_numpy(A)
        elapsed = time.perf_counter() - start

        # Ergebnis prüfen: A = P@L@U (scipy.linalg.lu-Konvention)
        reconstruction = P @ L @ U
        error = np.max(np.abs(reconstruction - A))
        assert error < 1e-8, f"LU-Zerlegung ungenau: max Fehler = {error:.2e}"

        # Hartes Limit
        assert elapsed < HARD_LIMITS['matrix_lu_100x100'], (
            f"LU(100x100) dauerte {elapsed:.3f}s "
            f"(Limit: {HARD_LIMITS['matrix_lu_100x100']}s)"
        )

        # Regression prüfen
        check_regression('matrix_lu_100x100', elapsed, baseline)
        save_baseline({'matrix_lu_100x100': elapsed})

    def test_matrix_lu_correctness(self):
        """Test: LU-Zerlegung ist mathematisch korrekt für kleine Matrix."""
        A = np.array([[2.0, 1.0, -1.0],
                      [-3.0, -1.0, 2.0],
                      [-2.0, 1.0, 2.0]])

        P, L, U = _lu_decomposition_numpy(A)

        # A = P @ L @ U (scipy-Konvention)
        reconstruction = P @ L @ U
        error = np.max(np.abs(reconstruction - A))
        assert error < 1e-10


def _lu_decomposition_numpy(A: np.ndarray) -> tuple:
    """
    @brief LU-Zerlegung via scipy.linalg.lu (Hilfsfunktion für Tests).
    @description
        scipy.linalg.lu gibt (P, L, U) mit A = P·L·U zurück,
        d.h. P ist so dass A = P·L·U (nicht P·A = L·U!).
        P ist hier eine Permutationsmatrix, nicht ihr Inverses.

    @param A: Eingangsmatrix.
    @return: Tuple (P, L, U) mit A = P·L·U.
    """
    from scipy.linalg import lu
    return lu(A)
