"""
@file profile_performance.py
@brief Performance-Profiling für alle Hauptmodule.
@description
    Misst Ausführungszeiten der wichtigsten Funktionen und identifiziert Bottlenecks.
    Vergleicht Ausführungszeit mit und ohne lru_cache (Cache-Warmup vs. Cold-Start).

    Getestete Funktionen:
        - is_prime()            aus algebra_numbertheory.py (Cache: 10000 Einträge)
        - prime_factorization() aus algebra_numbertheory.py (Cache: 1000 Einträge)
        - euler_phi()           aus algebra_numbertheory.py (Cache: 1000 Einträge)
        - bernoulli_number()    aus p_adic.py               (Cache: 200 Einträge)
        - p_adic_valuation()    aus p_adic.py               (Cache: 5000 Einträge)
        - collatz_sequence()    aus proof_theory.py         (Cache: 500 Einträge)
        - is_prime_fast()       aus proof_theory.py         (Cache: 10000 Einträge)

    Ausgabe: Tabelle mit Funktion, Aufrufanzahl, Gesamtzeit, Zeit/Aufruf,
             sowie cProfile-Zusammenfassung der teuersten Funktionen.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import cProfile
import pstats
import io
import time
import sys
import os
import functools

# Projektpfad zum sys.path hinzufügen (Skript liegt in /debugging/, Module in /src/)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Module importieren
from algebra_numbertheory import is_prime, prime_factorization, euler_phi
from p_adic import p_adic_valuation, bernoulli_number
from proof_theory import collatz_sequence, is_prime_fast


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def benchmark(func, args_list: list, label: str) -> dict:
    """
    @brief Führt eine Funktion für eine Liste von Argumenten aus und misst die Zeit.
    @description
        Führt func(*args) für jedes args-Tupel in args_list aus.
        Misst:
            - Gesamtausführungszeit (Cold-Start, ohne Cache-Vorwärmung)
            - Gesamtausführungszeit (Warm-Start, Cache bereits befüllt)
            - Zeit pro Aufruf (Cold / Warm)

        Der erste Durchlauf (Cold) befüllt den lru_cache.
        Der zweite Durchlauf (Warm) liest ausschließlich aus dem Cache.

    @param func: Zu testende Funktion (mit @lru_cache dekoriert)
    @param args_list: Liste von Argument-Tupeln, z.B. [(17,), (100,), ...]
    @param label: Anzeigename für die Tabelle
    @return: Dictionary mit Timing-Ergebnissen
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n_calls = len(args_list)

    # --- Cold-Start: Cache ist leer (oder wurde geleert) ---
    # Cache leeren falls vorhanden (nur bei lru_cache-dekorierten Funktionen)
    if hasattr(func, 'cache_clear'):
        func.cache_clear()

    start_cold = time.perf_counter()
    for args in args_list:
        # Einheitliche Aufrufstruktur: Tupel-Entpackung
        if isinstance(args, tuple):
            func(*args)
        else:
            func(args)
    end_cold = time.perf_counter()
    time_cold = end_cold - start_cold

    # --- Warm-Start: Cache ist bereits befüllt ---
    start_warm = time.perf_counter()
    for args in args_list:
        if isinstance(args, tuple):
            func(*args)
        else:
            func(args)
    end_warm = time.perf_counter()
    time_warm = end_warm - start_warm

    # Cache-Statistik abrufen (falls lru_cache)
    cache_info = None
    if hasattr(func, 'cache_info'):
        cache_info = func.cache_info()

    return {
        'label': label,
        'n_calls': n_calls,
        'time_cold_s': time_cold,
        'time_warm_s': time_warm,
        'time_per_call_cold_us': (time_cold / n_calls) * 1e6,
        'time_per_call_warm_us': (time_warm / n_calls) * 1e6,
        'speedup': time_cold / time_warm if time_warm > 0 else float('inf'),
        'cache_info': cache_info,
    }


def print_table(results: list):
    """
    @brief Gibt die Benchmark-Ergebnisse als formatierte Tabelle aus.
    @param results: Liste von Dictionaries aus benchmark()
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Tabellenbreite berechnen
    col_widths = {
        'label':    28,
        'n_calls':   8,
        'cold_ms':  10,
        'warm_ms':  10,
        'cold_us':  12,
        'warm_us':  12,
        'speedup':  10,
    }

    # Trennlinie
    total_width = sum(col_widths.values()) + len(col_widths) * 3
    line = '-' * total_width

    # Kopfzeile
    print(line)
    print(
        f"{'Funktion':<{col_widths['label']}} | "
        f"{'Aufrufe':>{col_widths['n_calls']}} | "
        f"{'Cold [ms]':>{col_widths['cold_ms']}} | "
        f"{'Warm [ms]':>{col_widths['warm_ms']}} | "
        f"{'Cold/Aufruf[µs]':>{col_widths['cold_us']}} | "
        f"{'Warm/Aufruf[µs]':>{col_widths['warm_us']}} | "
        f"{'Speedup':>{col_widths['speedup']}}"
    )
    print(line)

    # Datenzeilen
    for r in results:
        print(
            f"{r['label']:<{col_widths['label']}} | "
            f"{r['n_calls']:>{col_widths['n_calls']}} | "
            f"{r['time_cold_s']*1000:>{col_widths['cold_ms']}.3f} | "
            f"{r['time_warm_s']*1000:>{col_widths['warm_ms']}.3f} | "
            f"{r['time_per_call_cold_us']:>{col_widths['cold_us']}.3f} | "
            f"{r['time_per_call_warm_us']:>{col_widths['warm_us']}.3f} | "
            f"{r['speedup']:>{col_widths['speedup']}.1f}x"
        )

    print(line)
    print()

    # Cache-Statistiken ausgeben
    print("Cache-Statistiken (hits / misses / maxsize / currsize):")
    for r in results:
        ci = r['cache_info']
        if ci is not None:
            print(
                f"  {r['label']:<28}: "
                f"hits={ci.hits:>8}, misses={ci.misses:>6}, "
                f"maxsize={ci.maxsize}, currsize={ci.currsize}"
            )
        else:
            print(f"  {r['label']:<28}: kein lru_cache")
    print()


# ===========================================================================
# BENCHMARK-KONFIGURATION
# ===========================================================================

def generate_test_cases():
    """
    @brief Erstellt Testfälle für alle Benchmark-Funktionen.
    @description
        Generiert Listen von Argument-Tupeln für jeden Benchmark.
        Enthält absichtlich viele Duplikate, um Cache-Effizienz zu messen.
    @return: Dictionary mit Testfälle-Listen je Funktion
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Testfälle für is_prime: Mischung aus kleinen, mittleren, großen Zahlen + Duplikate
    primes_test = (
        [(n,) for n in range(2, 500)]           # 498 Zahlen (einmalig)
        + [(n,) for n in range(2, 200)]          # 198 Duplikate (Cache-Test)
        + [(n,) for n in [97, 101, 103, 107, 109, 113] * 20]  # häufige Primzahlen
    )

    # Testfälle für prime_factorization: Zahlen mit bekannter Zerlegung
    factorization_test = (
        [(n,) for n in range(2, 200)]
        + [(n,) for n in range(2, 100)]          # Duplikate
        + [(360,), (720,), (1024,), (999,), (360,), (720,)] * 10
    )

    # Testfälle für euler_phi: φ(n) für n = 1..300
    euler_phi_test = (
        [(n,) for n in range(1, 301)]
        + [(n,) for n in range(1, 101)]          # Duplikate
    )

    # Testfälle für bernoulli_number: Gerade Indizes (ungerade B_n = 0 für n > 1)
    bernoulli_test = (
        [(n,) for n in range(0, 40, 2)]         # B_0 bis B_38
        + [(n,) for n in range(0, 20, 2)] * 5  # Häufige Duplikate
    )

    # Testfälle für p_adic_valuation: (n, p)-Paare mit p = 2, 3, 5, 7
    p_adic_val_test = (
        [(n, p) for n in range(1, 200) for p in [2, 3, 5]]    # 597 unique Paare
        + [(n, 2) for n in range(1, 100)]                        # Duplikate für p=2
        + [(12, 2), (12, 3), (12, 5)] * 50                      # häufige Wiederholung
    )

    # Testfälle für collatz_sequence: Zahlen 1..80 (kleine, da Folgen lang sein können)
    collatz_test = (
        [(n,) for n in range(1, 80)]
        + [(n,) for n in range(1, 40)]           # Duplikate
        + [(27,), (6,), (12,)] * 20              # Häufige Werte
    )

    # Testfälle für is_prime_fast: wie is_prime
    prime_fast_test = (
        [(n,) for n in range(2, 500)]
        + [(n,) for n in range(2, 200)]
        + [(97,), (101,), (103,)] * 30
    )

    return {
        'is_prime':            primes_test,
        'prime_factorization': factorization_test,
        'euler_phi':           euler_phi_test,
        'bernoulli_number':    bernoulli_test,
        'p_adic_valuation':    p_adic_val_test,
        'collatz_sequence':    collatz_test,
        'is_prime_fast':       prime_fast_test,
    }


# ===========================================================================
# CPROFILER-ANALYSE
# ===========================================================================

def run_cprofiler(func, args_list: list, label: str, top_n: int = 10):
    """
    @brief Führt cProfile für eine Funktion aus und gibt die teuersten Aufrufe aus.
    @param func: Zu profilierende Funktion
    @param args_list: Argument-Tupel-Liste
    @param label: Anzeigename
    @param top_n: Anzahl der anzuzeigenden teuersten Funktionen
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Cache leeren für saubere Profiling-Ergebnisse (Cold-Start)
    if hasattr(func, 'cache_clear'):
        func.cache_clear()

    # cProfile starten
    profiler = cProfile.Profile()
    profiler.enable()

    # Funktion ausführen
    for args in args_list:
        if isinstance(args, tuple):
            func(*args)
        else:
            func(args)

    profiler.disable()

    # Ergebnisse als String aufbereiten
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream)
    stats.sort_stats('cumulative')
    stats.print_stats(top_n)

    print(f"\n--- cProfile: {label} (Top {top_n} nach kumulierter Zeit) ---")
    print(stream.getvalue())


# ===========================================================================
# HAUPTPROGRAMM
# ===========================================================================

def main():
    """
    @brief Startet das komplette Performance-Profiling.
    @description
        1. Generiert Testfälle für alle Benchmark-Funktionen
        2. Führt Cold-Start und Warm-Start Benchmarks durch
        3. Gibt formatierte Ergebnistabelle aus
        4. Führt cProfile für die langsamsten Funktionen durch
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    print("=" * 80)
    print("  Performance-Profiling: specialist-maths")
    print("  Vergleich: Cold-Start (Cache leer) vs. Warm-Start (Cache befüllt)")
    print("=" * 80)
    print()

    # Testfälle generieren
    print("Generiere Testfälle...")
    test_cases = generate_test_cases()

    # Benchmarks definieren: (Funktion, Testfälle-Schlüssel, Anzeigename)
    benchmarks = [
        (is_prime,            'is_prime',            'is_prime (algebra_nt)'),
        (prime_factorization, 'prime_factorization', 'prime_factorization'),
        (euler_phi,           'euler_phi',            'euler_phi'),
        (bernoulli_number,    'bernoulli_number',     'bernoulli_number (p_adic)'),
        (p_adic_valuation,    'p_adic_valuation',     'p_adic_valuation'),
        (collatz_sequence,    'collatz_sequence',     'collatz_sequence'),
        (is_prime_fast,       'is_prime_fast',        'is_prime_fast (proof)'),
    ]

    # Alle Benchmarks ausführen
    print("Starte Benchmarks...\n")
    results = []
    for func, key, label in benchmarks:
        args = test_cases[key]
        try:
            result = benchmark(func, args, label)
            results.append(result)
            print(f"  ✓ {label} ({result['n_calls']} Aufrufe)")
        except Exception as exc:
            print(f"  ✗ {label}: FEHLER – {exc}")

    # Ergebnistabelle ausgeben
    print()
    print("ERGEBNISTABELLE:")
    print_table(results)

    # cProfile für die drei teuersten Funktionen (nach Cold-Start-Zeit)
    sorted_results = sorted(results, key=lambda r: r['time_cold_s'], reverse=True)
    top3 = sorted_results[:3]

    print("=" * 80)
    print("  cProfile-Analyse der drei teuersten Funktionen (Cold-Start)")
    print("=" * 80)

    # Mapping: Anzeigename → Funktion
    func_map = {label: func for func, _, label in benchmarks}

    for r in top3:
        label = r['label']
        key_map = {label_: key for _, key, label_ in benchmarks}
        key = key_map.get(label)
        if key is None:
            continue
        func = func_map.get(label)
        if func is None:
            continue
        args = test_cases[key]
        # Nur die ersten 200 Aufrufe für cProfile (Übersichtlichkeit)
        run_cprofiler(func, args[:200], label, top_n=8)

    print("=" * 80)
    print("  Profiling abgeschlossen.")
    print("=" * 80)


if __name__ == '__main__':
    main()
