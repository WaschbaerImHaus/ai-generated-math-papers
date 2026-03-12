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
    t_start = time.time()
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

        if n % 5000 == 0:
            print(f"  n={n:,}/{n_max:,} geprüft | {time.time()-t_start:.0f}s", flush=True)

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

    # Phase 1: Exakte Prüfung bis n=100_000 (mit gmpy2 BigInt)
    # Berndt & Galway (2000) verifizierten bis n=10^9 — wir prüfen exakt bis 100K
    # Laufzeit-Schätzung: ~30-60 Min (O(n²·log n) wegen wachsender BigInts)
    print("Phase 1: Exakte Prüfung n ≤ 100,000 (gmpy2 BigInt)...")
    t0 = time.time()
    new_solutions_exact = check_exact_small(100_000)
    print(f"  Phase 1 fertig: {time.time()-t0:.1f}s")
    print(f"  Neue Lösungen: {len(new_solutions_exact)}\n")

    # Hinweis: Modulare Ausschlüsse (Phase 2/3) sind für n > max(EXCLUSION_PRIMES)=97
    # nutzlos, da n! ≡ 0 (mod p) für alle p ≤ n → n!+1 ≡ 1 (mod p) ist immer QR.
    # Daher: nur exakte Prüfung ist informativ.

    total_elapsed = time.time() - t0
    print(f"{'='*60}")
    print(f"GESAMT: {total_elapsed:.1f}s ({total_elapsed/3600:.2f}h)")
    print(f"Bekannte Lösungen bestätigt: n ∈ {{4, 5, 7}}")
    print(f"Neue Lösungen (exakt): {len(new_solutions_exact)}")
    if not new_solutions_exact:
        print(f"BESTÄTIGT: n!+1 ≠ m² für alle n ∈ [8, 100,000]")

    result = {
        "conjecture": "Brocard-Ramanujan",
        "exact_limit": 100_000,
        "method": "gmpy2 BigInt exakt (isqrt-Vergleich)",
        "known_solutions": [4, 5, 7],
        "new_solutions_exact": new_solutions_exact,
        "confirmed_no_new_exact_to": 100_000,
        "note": "Modulare Ausschlüsse nutzlos für n>97 (n!≡0 mod p → n!+1≡1 immer QR)",
        "elapsed_seconds": round(total_elapsed, 1),
        "timestamp": datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis: {OUTPUT_FILE}")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()
