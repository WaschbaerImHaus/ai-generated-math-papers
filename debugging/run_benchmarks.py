"""
@file run_benchmarks.py
@brief Eigenständiges Skript zum Ausführen und Speichern von Benchmarks.
@description
    Misst die Laufzeit kritischer Funktionen und speichert / vergleicht
    Ergebnisse mit der Baseline in debugging/benchmark_baseline.json.

    Kann direkt ausgeführt werden:
        python3 debugging/run_benchmarks.py

    Optionen:
        --reset    Baseline zurücksetzen (aktuelle Zeiten als neue Baseline)
        --verbose  Detaillierte Ausgabe

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import json
import time
import argparse
import numpy as np

# Pfad zum src/-Verzeichnis hinzufügen
SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src')
sys.path.insert(0, SRC_DIR)

# Pfad zur Baseline-Datei
BASELINE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'benchmark_baseline.json')

# Regression-Schwelle: 150%
REGRESSION_THRESHOLD = 1.5


def load_baseline() -> dict:
    """
    Lädt die gespeicherte Baseline aus der JSON-Datei.

    @return: Dictionary {name: sekunden} oder leeres Dict falls nicht vorhanden
    @lastModified 2026-03-10
    """
    if os.path.exists(BASELINE_PATH):
        with open(BASELINE_PATH, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}


def save_baseline(data: dict) -> None:
    """
    Speichert Benchmark-Ergebnisse als Baseline.

    @param data: Dictionary {name: sekunden}
    @lastModified 2026-03-10
    """
    with open(BASELINE_PATH, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)
    print(f"Baseline gespeichert: {BASELINE_PATH}")


def run_benchmark(name: str, func, *args, **kwargs) -> tuple[float, object]:
    """
    Führt eine Funktion aus und misst die Laufzeit.

    @param name: Bezeichner der Funktion (für Anzeige)
    @param func: Auszuführende Funktion
    @return: Tupel (laufzeit_sekunden, ergebnis)
    @lastModified 2026-03-10
    """
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return elapsed, result


def format_time(seconds: float) -> str:
    """
    Formatiert eine Zeit in lesbares Format.

    @param seconds: Zeit in Sekunden
    @return: Formatierter String (z.B. "1.23 ms" oder "456 μs")
    @lastModified 2026-03-10
    """
    if seconds >= 1.0:
        return f"{seconds:.3f} s"
    elif seconds >= 0.001:
        return f"{seconds * 1000:.2f} ms"
    else:
        return f"{seconds * 1_000_000:.1f} μs"


def main():
    """
    Hauptfunktion: Führt alle Benchmarks aus und gibt Ergebnisse aus.

    @lastModified 2026-03-10
    """
    # Argumente parsen
    parser = argparse.ArgumentParser(description='Benchmark-Runner für specialist-maths')
    parser.add_argument('--reset', action='store_true',
                        help='Baseline zurücksetzen (neue Baseline aus aktuellen Zeiten)')
    parser.add_argument('--verbose', action='store_true',
                        help='Detaillierte Ausgabe')
    args = parser.parse_args()

    print("=" * 70)
    print("  specialist-maths Benchmark-Suite")
    print("=" * 70)

    # Baseline laden (falls nicht --reset)
    baseline = {} if args.reset else load_baseline()
    if baseline:
        print(f"Baseline geladen: {BASELINE_PATH}")
    else:
        print("Keine Baseline vorhanden – neue Baseline wird erstellt.")
    print()

    # Ergebnisse sammeln
    current_times: dict[str, float] = {}
    regressions: list[str] = []

    # ------------------------------------------------------------------
    # BENCHMARK 1: Primzahlzählfunktion
    # ------------------------------------------------------------------
    print("1. prime_counting_function(1000) ...")
    from analytic_number_theory import prime_counting_function
    elapsed, result = run_benchmark('prime_counting_function_1000', prime_counting_function, 1000)
    current_times['prime_counting_function_1000'] = elapsed
    status = "OK" if result == 168 else f"FEHLER (erwartet 168, erhalten {result})"
    print(f"   Zeit: {format_time(elapsed)}, π(1000) = {result} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 2: Fourier-Koeffizienten der Delta-Funktion
    # ------------------------------------------------------------------
    print("2. fourier_coefficients_delta(30) ...")
    from modular_forms import fourier_coefficients_delta
    elapsed, result = run_benchmark(
        'fourier_coefficients_delta_30', fourier_coefficients_delta, 30
    )
    current_times['fourier_coefficients_delta_30'] = elapsed
    ok = (result[0] == 1 and result[1] == -24)
    status = "OK" if ok else f"FEHLER (τ(1)={result[0]}, τ(2)={result[1]})"
    print(f"   Zeit: {format_time(elapsed)}, τ(1)={result[0]}, τ(2)={result[1]} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 3: Riemann-Zeta-Funktion
    # ------------------------------------------------------------------
    print("3. riemann_zeta(2.0) ...")
    from complex_analysis import riemann_zeta
    elapsed, result = run_benchmark('riemann_zeta_2', riemann_zeta, 2.0)
    current_times['riemann_zeta_2'] = elapsed
    expected_zeta2 = np.pi**2 / 6
    err = abs(result - expected_zeta2)
    status = "OK" if err < 1e-6 else f"FEHLER (Fehler={err:.2e})"
    print(f"   Zeit: {format_time(elapsed)}, ζ(2) = {result:.10f} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 4: Von-Mangoldt-Funktion in Schleife
    # ------------------------------------------------------------------
    print("4. von_mangoldt_function(n) für n=1..500 ...")
    from analytic_number_theory import von_mangoldt_function
    def mangoldt_loop():
        return [von_mangoldt_function(n) for n in range(1, 501)]
    elapsed, result = run_benchmark('von_mangoldt_loop_500', mangoldt_loop)
    current_times['von_mangoldt_loop_500'] = elapsed
    ok = (abs(result[1] - np.log(2)) < 1e-10 and abs(result[5]) < 1e-10)
    status = "OK" if ok else "FEHLER"
    print(f"   Zeit: {format_time(elapsed)}, Λ(2)={result[1]:.6f}, Λ(6)={result[5]:.6f} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 5: numpy FFT (4096 Punkte)
    # ------------------------------------------------------------------
    print("5. numpy FFT (4096 Punkte) ...")
    np.random.seed(42)
    signal = np.random.rand(4096)
    def fft_roundtrip():
        return np.fft.ifft(np.fft.fft(signal)).real
    elapsed, result = run_benchmark('numpy_fft_4096', fft_roundtrip)
    current_times['numpy_fft_4096'] = elapsed
    ok = np.allclose(result, signal, atol=1e-10)
    status = "OK" if ok else "FEHLER (Rundreise nicht exakt)"
    print(f"   Zeit: {format_time(elapsed)}, FFT-Rundreise [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 6: auto_compute symbolisch
    # ------------------------------------------------------------------
    print("6. auto_compute(x^5 + 3x^3 - 2x + 1, x=2.0) ...")
    import sympy as sp
    from config import auto_compute
    x = sp.Symbol('x')
    expr_sym = x**5 + 3*x**3 - 2*x + 1
    elapsed, result = run_benchmark('auto_compute_symbolic_deg5', auto_compute, expr_sym, 2.0)
    current_times['auto_compute_symbolic_deg5'] = elapsed
    ok = (abs(result.real - 53.0) < 1e-8)
    status = "OK" if ok else f"FEHLER ({result} ≠ 53)"
    print(f"   Zeit: {format_time(elapsed)}, Ergebnis = {result.real:.1f} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 7: auto_compute numerisch
    # ------------------------------------------------------------------
    print("7. auto_compute(x^15, x=2.0) ...")
    expr_num = x**15
    elapsed, result = run_benchmark('auto_compute_numeric_deg15', auto_compute, expr_num, 2.0)
    current_times['auto_compute_numeric_deg15'] = elapsed
    ok = (abs(result - 32768.0) < 1.0)
    status = "OK" if ok else f"FEHLER ({result} ≠ 32768)"
    print(f"   Zeit: {format_time(elapsed)}, Ergebnis = {result:.0f} [{status}]")

    # ------------------------------------------------------------------
    # BENCHMARK 8: parallel_integrate (sequenziell, 1 Worker)
    # ------------------------------------------------------------------
    print("8. parallel_integrate(['x**2', 'sin(x)', 'exp(x)'], max_workers=1) ...")
    from parallel_compute import parallel_integrate
    elapsed, result = run_benchmark(
        'parallel_integrate_3', parallel_integrate, ['x**2', 'sin(x)', 'exp(x)'], 'x', 1
    )
    current_times['parallel_integrate_3'] = elapsed
    ok = (len(result) == 3)
    status = "OK" if ok else "FEHLER"
    print(f"   Zeit: {format_time(elapsed)}, {len(result)} Integrale [{status}]")

    # ------------------------------------------------------------------
    # REGRESSIONEN PRÜFEN
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("  Regressionsanalyse")
    print("=" * 70)

    has_regression = False
    for name, elapsed in current_times.items():
        if name in baseline and baseline[name] > 0:
            ratio = elapsed / baseline[name]
            baseline_str = format_time(baseline[name])
            current_str = format_time(elapsed)
            if ratio > REGRESSION_THRESHOLD:
                print(f"  [REGRESSION] {name}: {current_str} vs. {baseline_str} ({ratio:.1f}×)")
                regressions.append(name)
                has_regression = True
            elif ratio < 0.9:
                print(f"  [SCHNELLER]  {name}: {current_str} vs. {baseline_str} ({ratio:.2f}×)")
            else:
                if args.verbose:
                    print(f"  [OK]         {name}: {current_str} vs. {baseline_str} ({ratio:.2f}×)")
        else:
            if args.verbose:
                print(f"  [NEU]        {name}: {format_time(elapsed)} (kein Baseline-Vergleich)")

    if not has_regression:
        print("  Keine Regressionen gefunden.")

    # ------------------------------------------------------------------
    # BASELINE AKTUALISIEREN
    # ------------------------------------------------------------------
    print()
    if args.reset or not baseline:
        save_baseline(current_times)
    else:
        # Nur aktualisieren wenn neue Zeiten vorhanden
        updated = dict(baseline)
        for name, t in current_times.items():
            if name not in updated:
                updated[name] = t
            else:
                # Schnellere Zeit als neues Minimum übernehmen
                updated[name] = min(updated[name], t)
        save_baseline(updated)

    print("Benchmarks abgeschlossen.")
    return 0 if not has_regression else 1


if __name__ == '__main__':
    sys.exit(main())
