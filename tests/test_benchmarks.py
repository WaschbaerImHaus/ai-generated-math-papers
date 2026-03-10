"""
@file test_benchmarks.py
@brief Benchmark-Regressions-Tests für kritische Funktionen.
@description
    Misst die Laufzeit von 8 kritischen Funktionen und vergleicht sie mit
    einer gespeicherten Baseline. Warnt (schlägt aber nicht fehl) wenn die
    aktuelle Laufzeit > 150% der Baseline ist.

    Getestete Funktionen:
        1. prime_counting_function(1000)
        2. fourier_coefficients_delta(30)
        3. riemann_zeta(2.0) aus complex_analysis
        4. von_mangoldt_function in Schleife (n=1..500)
        5. FFT via numpy (n=4096 Punkte)
        6. auto_compute (symbolisch, Grad 5)
        7. auto_compute (numerisch, Grad 15)
        8. parallel_integrate (3 Ausdrücke)

    Baseline-Datei: debugging/benchmark_baseline.json

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import json
import time
import warnings
import pytest
import numpy as np

# Pfad zum src/-Verzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Pfad zur Baseline-Datei
BASELINE_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'debugging', 'benchmark_baseline.json'
)

# Warnung bei > 150% der Baseline-Zeit
REGRESSION_THRESHOLD = 1.5


def load_baseline() -> dict:
    """
    Lädt die Baseline-Zeiten aus der JSON-Datei.

    @return: Dictionary {funktionsname: baseline_sekunden} oder leer, falls keine Baseline
    @lastModified 2026-03-10
    """
    if not os.path.exists(BASELINE_PATH):
        return {}
    with open(BASELINE_PATH, 'r', encoding='utf-8') as f:
        return json.load(f)


def save_baseline(data: dict) -> None:
    """
    Speichert Benchmark-Zeiten als neue Baseline-Datei.

    @param data: Dictionary {funktionsname: sekunden}
    @lastModified 2026-03-10
    """
    os.makedirs(os.path.dirname(BASELINE_PATH), exist_ok=True)
    with open(BASELINE_PATH, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)


def check_regression(name: str, elapsed: float, baseline: dict) -> None:
    """
    Vergleicht aktuelle Laufzeit mit der Baseline und gibt ggf. Warnung aus.

    Warnt (via warnings.warn), schlägt aber nicht fehl. So blockieren
    Performance-Regressionen nicht den Build, machen aber auf Probleme aufmerksam.

    @param name: Name der Funktion
    @param elapsed: Aktuelle Laufzeit in Sekunden
    @param baseline: Dictionary mit Baseline-Werten
    @lastModified 2026-03-10
    """
    if name in baseline and baseline[name] > 0:
        ratio = elapsed / baseline[name]
        if ratio > REGRESSION_THRESHOLD:
            warnings.warn(
                f"PERFORMANCE-REGRESSION: {name} ist {ratio:.1f}× langsamer als Baseline "
                f"({elapsed:.4f}s vs. {baseline[name]:.4f}s Baseline). "
                f"Schwelle: {REGRESSION_THRESHOLD:.0%}.",
                UserWarning,
                stacklevel=2
            )


# ===========================================================================
# BENCHMARK-TESTS
# ===========================================================================

class TestBenchmarks:
    """Benchmark-Regressions-Tests für kritische Funktionen."""

    # Klassenvariable: Baseline einmal laden, nicht bei jedem Test
    _baseline: dict = {}
    _current_times: dict = {}

    @classmethod
    def setup_class(cls):
        """Baseline laden bevor Tests laufen."""
        cls._baseline = load_baseline()
        cls._current_times = {}

    @classmethod
    def teardown_class(cls):
        """Neue Baseline speichern (falls noch keine vorhanden oder deutlich schneller)."""
        if cls._current_times:
            if not cls._baseline:
                # Keine Baseline vorhanden: aktuelle Zeiten als Baseline speichern
                save_baseline(cls._current_times)
            else:
                # Baseline aktualisieren: Minimum aus alter und neuer Zeit
                updated = {}
                for name, t in cls._current_times.items():
                    if name in cls._baseline:
                        # Nur aktualisieren wenn deutlich schneller (10% Puffer)
                        updated[name] = min(cls._baseline[name], t)
                    else:
                        updated[name] = t
                save_baseline(updated)

    def _time_it(self, name: str, func, *args, **kwargs):
        """
        Misst die Laufzeit einer Funktion und prüft gegen Baseline.

        @param name: Bezeichner für die Funktion
        @param func: Aufrufbare Funktion
        @param args/kwargs: Argumente für func
        @return: Rückgabewert von func
        @lastModified 2026-03-10
        """
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start

        # Zeiten speichern
        TestBenchmarks._current_times[name] = elapsed

        # Gegen Baseline prüfen (Warnung, kein Fehler)
        check_regression(name, elapsed, self._baseline)

        return result

    # ------------------------------------------------------------------
    # BENCHMARK 1: Primzahlzählfunktion
    # ------------------------------------------------------------------

    def test_benchmark_prime_counting(self):
        """Benchmark: prime_counting_function(1000) – Sieb des Eratosthenes."""
        from analytic_number_theory import prime_counting_function

        result = self._time_it('prime_counting_function_1000', prime_counting_function, 1000)
        # Korrektheitsprüfung: π(1000) = 168
        assert result == 168, f"prime_counting_function(1000) = {result}, erwartet 168"

    # ------------------------------------------------------------------
    # BENCHMARK 2: Fourier-Koeffizienten der Delta-Funktion (FFT-basiert)
    # ------------------------------------------------------------------

    def test_benchmark_fourier_coefficients_delta(self):
        """Benchmark: fourier_coefficients_delta(30) – FFT-basierte Berechnung."""
        from modular_forms import fourier_coefficients_delta

        result = self._time_it(
            'fourier_coefficients_delta_30', fourier_coefficients_delta, 30
        )
        # τ(1) = 1, τ(2) = -24 (bekannte Werte)
        assert result[0] == 1,  f"τ(1) = {result[0]}, erwartet 1"
        assert result[1] == -24, f"τ(2) = {result[1]}, erwartet -24"

    # ------------------------------------------------------------------
    # BENCHMARK 3: Riemann-Zeta-Funktion
    # ------------------------------------------------------------------

    def test_benchmark_riemann_zeta(self):
        """Benchmark: riemann_zeta(2.0) – Euler-Knopp-Beschleunigung."""
        from complex_analysis import riemann_zeta

        result = self._time_it('riemann_zeta_2', riemann_zeta, 2.0)
        # ζ(2) = π²/6 ≈ 1.6449340668...
        expected = np.pi**2 / 6
        # Toleranz 1e-3: Die Implementierung nutzt Euler-Knopp mit begrenzten Termen
        assert abs(result - expected) < 1e-3, f"ζ(2) = {result}, erwartet ≈ {expected}"

    # ------------------------------------------------------------------
    # BENCHMARK 4: Von-Mangoldt-Funktion in Schleife
    # ------------------------------------------------------------------

    def test_benchmark_von_mangoldt_loop(self):
        """Benchmark: von_mangoldt_function für n = 1..500."""
        from analytic_number_theory import von_mangoldt_function

        def run_loop():
            return [von_mangoldt_function(n) for n in range(1, 501)]

        result = self._time_it('von_mangoldt_loop_500', run_loop)
        # Λ(2) = ln(2), Λ(4) = ln(2) (da 4 = 2²), Λ(6) = 0 (zwei Primteiler)
        assert abs(result[1] - np.log(2)) < 1e-10   # Λ(2) = ln(2)
        assert abs(result[5] - 0.0) < 1e-10          # Λ(6) = 0

    # ------------------------------------------------------------------
    # BENCHMARK 5: FFT mit numpy (4096 Punkte)
    # ------------------------------------------------------------------

    def test_benchmark_fft_numpy(self):
        """Benchmark: numpy FFT mit 4096 Punkten."""
        signal = np.random.rand(4096)

        def run_fft():
            spectrum = np.fft.fft(signal)
            recovered = np.fft.ifft(spectrum)
            return recovered

        result = self._time_it('numpy_fft_4096', run_fft)
        # Rundreise FFT → IFFT muss Signal rekonstruieren
        assert np.allclose(result.real, signal, atol=1e-10), "FFT-Rundreise fehlgeschlagen"

    # ------------------------------------------------------------------
    # BENCHMARK 6: auto_compute symbolisch (Grad 5)
    # ------------------------------------------------------------------

    def test_benchmark_auto_compute_symbolic(self):
        """Benchmark: auto_compute mit symbolischem Pfad (Grad 5 ≤ 10)."""
        import sympy as sp
        from config import auto_compute

        x = sp.Symbol('x')
        expr = x**5 + 3*x**3 - 2*x + 1  # Grad 5, wenige Ops → symbolisch

        result = self._time_it('auto_compute_symbolic_deg5', auto_compute, expr, 2.0)
        # x=2: 32 + 24 - 4 + 1 = 53
        assert abs(result.real - 53.0) < 1e-8, f"auto_compute(x^5+..., 2) = {result}, erwartet 53"

    # ------------------------------------------------------------------
    # BENCHMARK 7: auto_compute numerisch (Grad 15)
    # ------------------------------------------------------------------

    def test_benchmark_auto_compute_numeric(self):
        """Benchmark: auto_compute mit numerischem Pfad (Grad 15 > 10)."""
        import sympy as sp
        from config import auto_compute

        x = sp.Symbol('x')
        expr = x**15  # Grad 15 > AUTO_COMPUTE_MAX_DEGREE → numerisch

        result = self._time_it('auto_compute_numeric_deg15', auto_compute, expr, 2.0)
        # 2^15 = 32768
        assert abs(result - 32768.0) < 1.0, f"auto_compute(x^15, 2) = {result}, erwartet 32768"

    # ------------------------------------------------------------------
    # BENCHMARK 8: parallel_integrate (3 Ausdrücke)
    # ------------------------------------------------------------------

    def test_benchmark_parallel_integrate(self):
        """Benchmark: parallel_integrate mit 3 Ausdrücken (1 Worker, stabil)."""
        from parallel_compute import parallel_integrate

        exprs = ['x**2', 'sin(x)', 'exp(x)']
        results = self._time_it(
            'parallel_integrate_3', parallel_integrate, exprs, 'x', 1
        )
        # Alle Ausdrücke müssen integriert worden sein
        assert len(results) == 3
        assert 'x**2' in results
        assert 'sin(x)' in results
        assert 'exp(x)' in results
