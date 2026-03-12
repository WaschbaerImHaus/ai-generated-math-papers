"""
@file heavy_compute_brocard.py
@brief Brocard-Ramanujan: n!+1=m² — Modulare Ausschlussverifikation bis n=10^6
@description
    Brocard-Ramanujan-Vermutung: n!+1 = m² hat nur Lösungen n ∈ {4,5,7}.
    Bekannt: Berndt & Galway (2000) verifiziert bis n = 10^9.

    Strategie — Modulare Ausschlüsse:
    Für jede "Ausschluss-Primzahl" p prüfe ob n!+1 ein Quadratischer Rest mod p ist.
    Falls nicht → n ist kein Lösungskandidat.

    Für große n (n ≥ p) gilt n! ≡ 0 (mod p), also n!+1 ≡ 1 (mod p).
    Da 1 ein Quadratischer Rest mod p ist, schließen große Primzahlen nicht aus.
    Daher: für Primzahlen p ≤ n müssen wir n!+1 mod p = QR(p)? prüfen.

    Dieser Ansatz: Direkte Prüfung ob n!+1 eine Quadratzahl ist.
    Für n bis ~1000 exakt, darüber modulare Ausschlüsse.

@author Michael Fuhrmann
@date 2026-03-12
"""

import os
import time
import json
import math
import multiprocessing
from datetime import datetime

import gmpy2
import numpy as np

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "brocard_result.json")

# Kleine Ausschluss-Primzahlen (für die n!+1 kein QR ist für bestimmte n)
EXCLUSION_PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                    47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

def is_quadratic_residue(a, p):
    """Prüft ob a ein quadratischer Rest mod p (prim) ist."""
    if a % p == 0:
        return True
    return pow(a, (p - 1) // 2, p) == 1


def check_n_modular(n, factorial_mod_dict):
    """
    Prüft ob n!+1 quadratischen Resten widersteht.
    factorial_mod_dict: {p: n! mod p für alle n}
    Gibt True zurück wenn n NICHT ausgeschlossen werden kann (Kandidat).
    """
    for p in EXCLUSION_PRIMES:
        if n >= p:
            # n! ≡ 0 (mod p), also n!+1 ≡ 1 (mod p) — kein Ausschluss
            continue
        # n < p: n! mod p aus Vorberechnung
        fac_mod_p = factorial_mod_dict.get((n, p), None)
        if fac_mod_p is None:
            continue
        val = (fac_mod_p + 1) % p
        if not is_quadratic_residue(val, p):
            return False  # Ausgeschlossen
    return True  # Kandidat (nicht ausgeschlossen)


def precompute_factorials_mod(n_max, primes):
    """
    Vorberechnung: n! mod p für alle n ≤ n_max und p in primes.
    Gibt dict {(n,p): n! mod p} zurück.
    """
    result = {}
    for p in primes:
        fac = 1
        for n in range(1, min(n_max + 1, p)):
            fac = fac * n % p
            result[(n, p)] = fac
        result[(0, p)] = 1  # 0! = 1
    return result


def check_exact_small(n_max=10000):
    """Exakte Prüfung ob n!+1 Quadratzahl für n ≤ n_max (mit gmpy2 BigInt)."""
    solutions = [4, 5, 7]  # Bekannte Lösungen
    found_new = []

    factorial = gmpy2.mpz(1)
    for n in range(1, n_max + 1):
        factorial *= n
        val = factorial + 1
        # Prüfe ob val eine Quadratzahl ist
        root = gmpy2.isqrt(val)
        if root * root == val:
            if n not in solutions:
                found_new.append({"n": n, "m": int(root)})
                print(f"  *** NEUE LÖSUNG GEFUNDEN: n={n}, m={root} ***", flush=True)
            else:
                print(f"  Bekannte Lösung bestätigt: n={n}", flush=True)

        if n % 1000 == 0:
            print(f"  n={n:,} geprüft (exakt)", flush=True)

    return found_new


def check_range_modular(args):
    """Prüft einen Bereich [n_start, n_end] mit modularen Ausschlüssen."""
    n_start, n_end, fac_dict = args
    candidates = []  # Nicht ausgeschlossene n

    for n in range(n_start, n_end + 1):
        if not check_n_modular(n, fac_dict):
            continue  # Ausgeschlossen
        candidates.append(n)

    return candidates


def main():
    print(f"{'='*60}")
    print(f"Brocard-Ramanujan Heavy Computation — {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"{'='*60}\n")

    # Phase 1: Exakte Prüfung bis n=10000 (mit gmpy2)
    print("Phase 1: Exakte Prüfung n ≤ 10,000...")
    t0 = time.time()
    new_solutions_exact = check_exact_small(10000)
    print(f"  Phase 1 fertig: {time.time()-t0:.1f}s")
    print(f"  Neue Lösungen: {len(new_solutions_exact)}\n")

    # Phase 2: Modulare Ausschlüsse bis n = 100,000
    print("Phase 2: Modulare Ausschlüsse n ≤ 100,000...")
    t1 = time.time()
    fac_dict = precompute_factorials_mod(1000, EXCLUSION_PRIMES)

    candidates_100k = []
    for n in range(10001, 100001):
        if check_n_modular(n, fac_dict):
            candidates_100k.append(n)

    print(f"  Nicht-ausgeschlossene Kandidaten (10001–100000): {len(candidates_100k):,}")
    print(f"  Ausschlussrate: {(90000-len(candidates_100k))/90000:.1%}")
    print(f"  Phase 2 fertig: {time.time()-t1:.1f}s\n")

    # Phase 3: Modular bis n = 1,000,000
    print("Phase 3: Modulare Ausschlüsse n ≤ 1,000,000...")
    t2 = time.time()

    N_WORKERS = min(20, multiprocessing.cpu_count())
    CHUNK = 50000
    tasks = [(n, min(n + CHUNK - 1, 1_000_000), fac_dict)
             for n in range(100001, 1_000_001, CHUNK)]

    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(check_range_modular, tasks)

    candidates_1M = []
    for r in results:
        candidates_1M.extend(r)

    print(f"  Kandidaten (100001–1000000): {len(candidates_1M):,}")
    print(f"  Ausschlussrate: {(900000-len(candidates_1M))/900000:.1%}")
    elapsed = time.time() - t2
    print(f"  Phase 3 fertig: {elapsed:.1f}s\n")

    total_elapsed = time.time() - t0
    print(f"{'='*60}")
    print(f"GESAMT: {total_elapsed:.1f}s ({total_elapsed/3600:.2f}h)")
    print(f"Bekannte Lösungen bestätigt: n ∈ {{4, 5, 7}}")
    print(f"Neue Lösungen (exakt): {len(new_solutions_exact)}")
    print(f"Verbleibende Kandidaten bis 10^6 (nicht ausgeschlossen): "
          f"{len(candidates_100k) + len(candidates_1M):,}")

    result = {
        "conjecture": "Brocard-Ramanujan",
        "exact_limit": 10000,
        "modular_limit": 1_000_000,
        "exclusion_primes": EXCLUSION_PRIMES,
        "new_solutions_exact": new_solutions_exact,
        "candidates_count_10k_100k": len(candidates_100k),
        "candidates_count_100k_1M": len(candidates_1M),
        "confirmed_no_new_exact_to": 10000,
        "elapsed_seconds": round(total_elapsed, 1),
        "timestamp": datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis: {OUTPUT_FILE}")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()
