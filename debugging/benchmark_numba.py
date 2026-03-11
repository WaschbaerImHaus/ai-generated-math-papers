"""
benchmark_numba.py — Geschwindigkeitsvergleich Python vs. Numba-JIT.

Misst den Speedup-Faktor der Numba-beschleunigten Implementierungen
gegenüber den reinen Python-Versionen für:
  1. Sieb des Eratosthenes (n = 10^6)
  2. Euler-Knopp-Beschleunigung für η(s) (100 Aufrufe)

@author: Michael Fuhrmann
@date: 2026-03-11
@lastModified: 2026-03-11
"""

import sys
import os
import time
import math

# Projekt-Root zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.numba_jit import (
    sieve_primes_fast,
    sieve_numpy,
    eta_euler_fast,
    eta_euler_accelerated_jit,
    NUMBA_AVAILABLE,
    warmup,
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen: Python-Referenzimplementierungen
# ---------------------------------------------------------------------------

def python_sieve(limit: int) -> list:
    """
    Reine Python-Implementierung des Siebs für den Vergleich.

    @param limit: Obere Grenze (inklusiv)
    @return: Liste aller Primzahlen ≤ limit
    """
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    p = 2
    while p * p <= limit:
        if is_prime[p]:
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False
        p += 1
    return [n for n in range(2, limit + 1) if is_prime[n]]


def python_eta(s: complex, n: int = 60) -> complex:
    """
    Reine Python-Implementierung der Euler-Knopp-Eta für den Vergleich.

    @param s: Komplexe Zahl mit Re(s) > 0
    @param n: Anzahl Terme
    @return: η(s)
    """
    c = [1.0 / (complex(j + 1) ** s) for j in range(n + 1)]
    delta = list(c)
    result = complex(0.0)
    sign = 1
    power_half = 0.5
    for k in range(n):
        result += sign * power_half * delta[0]
        sign *= -1
        power_half *= 0.5
        for j in range(len(delta) - 1):
            delta[j] = delta[j + 1] - delta[j]
        delta.pop()
    return result


# ---------------------------------------------------------------------------
# Benchmark: Sieb
# ---------------------------------------------------------------------------

def benchmark_sieve(n: int = 1_000_000, repetitions: int = 3):
    """
    Vergleicht Python-Sieb vs. Numba-Sieb für n = 10^6.

    @param n: Siebgrenze
    @param repetitions: Anzahl Wiederholungen für stabilere Messung
    """
    print(f"\n{'='*60}")
    print(f"BENCHMARK: Sieb des Eratosthenes (n = {n:,})")
    print(f"{'='*60}")

    # --- Python-Sieb messen ---
    py_times = []
    for _ in range(repetitions):
        t0 = time.perf_counter()
        py_result = python_sieve(n)
        py_times.append(time.perf_counter() - t0)
    py_time = min(py_times)  # Bestes Ergebnis (keine GC-Störung)
    print(f"Python-Sieb:    {py_time*1000:.2f} ms  ({len(py_result)} Primzahlen)")

    # --- Numba-Sieb messen (ohne Warmup-Zeit) ---
    if NUMBA_AVAILABLE:
        numba_times = []
        for _ in range(repetitions):
            t0 = time.perf_counter()
            nb_result = sieve_primes_fast(n)
            numba_times.append(time.perf_counter() - t0)
        nb_time = min(numba_times)
        speedup = py_time / nb_time if nb_time > 0 else float("inf")
        print(f"Numba-Sieb:     {nb_time*1000:.2f} ms  ({len(nb_result)} Primzahlen)")
        print(f"Speedup:        {speedup:.1f}×")
        # Korrektheitsprüfung: letzte Primzahl muss übereinstimmen
        assert py_result[-1] == nb_result[-1], "FEHLER: Ergebnisse stimmen nicht überein!"
        print("Korrektheit:    OK (letzte Primzahl übereinstimmend)")
    else:
        print("Numba nicht verfügbar — nur Python-Variante gemessen.")


# ---------------------------------------------------------------------------
# Benchmark: Euler-Knopp η(s)
# ---------------------------------------------------------------------------

def benchmark_eta(repetitions: int = 100):
    """
    Vergleicht Python-eta vs. Numba-eta für 100 Aufrufe mit verschiedenen s.

    @param repetitions: Anzahl der Aufruf-Wiederholungen
    """
    print(f"\n{'='*60}")
    print(f"BENCHMARK: Euler-Knopp η(s) ({repetitions} Aufrufe)")
    print(f"{'='*60}")

    # Testeingaben: kritische Gerade und reelle Werte
    test_values = [
        0.5 + 14.135j,
        0.5 + 21.022j,
        2.0 + 0j,
        0.5 + 0.5j,
    ]

    # --- Python-eta messen ---
    t0 = time.perf_counter()
    for _ in range(repetitions):
        for s in test_values:
            python_eta(s, n=60)
    py_time = time.perf_counter() - t0
    total_calls = repetitions * len(test_values)
    print(f"Python-eta:     {py_time*1000:.2f} ms  ({total_calls} Aufrufe)")
    print(f"                {py_time/total_calls*1000:.3f} ms/Aufruf")

    # --- Numba-eta messen ---
    if NUMBA_AVAILABLE:
        t0 = time.perf_counter()
        for _ in range(repetitions):
            for s in test_values:
                eta_euler_fast(s, n=60)
        nb_time = time.perf_counter() - t0
        speedup = py_time / nb_time if nb_time > 0 else float("inf")
        print(f"Numba-eta:      {nb_time*1000:.2f} ms  ({total_calls} Aufrufe)")
        print(f"                {nb_time/total_calls*1000:.3f} ms/Aufruf")
        print(f"Speedup:        {speedup:.1f}×")

        # Korrektheitsprüfung: η(2) = π²/12
        expected = math.pi**2 / 12
        nb_val = eta_euler_fast(2+0j)
        py_val = python_eta(2+0j)
        diff = abs(nb_val.real - py_val.real)
        print(f"Korrektheit:    |Numba - Python| = {diff:.2e}  (OK wenn < 1e-12)")
    else:
        print("Numba nicht verfügbar — nur Python-Variante gemessen.")


# ---------------------------------------------------------------------------
# Warmup-Zeit messen
# ---------------------------------------------------------------------------

def measure_warmup():
    """
    Misst die JIT-Kompilierungszeit beim ersten Aufruf.
    """
    if not NUMBA_AVAILABLE:
        print("\nNumba nicht verfügbar — kein Warmup nötig.")
        return

    print(f"\n{'='*60}")
    print("WARMUP: JIT-Kompilierungszeit")
    print(f"{'='*60}")
    print("(Numba kompiliert beim ersten Aufruf — danach sofort)")

    t0 = time.perf_counter()
    warmup()
    warmup_time = time.perf_counter() - t0
    print(f"Warmup-Zeit:    {warmup_time*1000:.0f} ms")
    print("Nach Warmup: JIT-Funktionen sofort verfügbar.")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("specialist-maths — Numba-JIT Benchmark")
    print(f"Numba verfügbar: {NUMBA_AVAILABLE}")

    # Zuerst warmup damit JIT-Kompilierung nicht in Benchmarks fällt
    measure_warmup()

    # Benchmark 1: Sieb
    benchmark_sieve(n=1_000_000, repetitions=5)

    # Benchmark 2: Eta-Funktion
    benchmark_eta(repetitions=100)

    print(f"\n{'='*60}")
    print("Benchmark abgeschlossen.")
