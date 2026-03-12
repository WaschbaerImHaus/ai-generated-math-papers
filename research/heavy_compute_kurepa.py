"""
@file heavy_compute_kurepa.py
@brief Kurepa-Vermutung: !p ≢ 0 (mod p) — Verifikation bis p = 10_000_000
@description
    Kurepa-Vermutung (1950): Für jede Primzahl p gilt !p ≢ 0 (mod p),
    wobei !p = Σ_{k=0}^{p-1} k! die Linksfakultät ist.

    Strategie:
    - Wilson: (p-1)! ≡ -1 (mod p) → !p ≡ !(p-1) - 1 (mod p)
    - Numba JIT für innere Schleife: ~100x Speedup
    - Multiprocessing: 20 Kerne parallel

    Bekannt: Kein Gegenbeispiel bis p < 10^6 (Koukou & Prodinger u.a.)
    Ziel: Verifikation bis p = 10_000_000

@author Michael Fuhrmann
@date 2026-03-12
"""

import os
import sys
import time
import json
import multiprocessing
from datetime import datetime

import numba
import numpy as np
from sympy import primerange

# ─── Numba-beschleunigte Kernfunktion ────────────────────────────────────────

@numba.njit(cache=True)
def kurepa_check_range(primes_arr):
    """
    Prüft Kurepa-Vermutung für alle Primzahlen im Array.
    Gibt Array der Primzahlen zurück, für die !p ≡ 0 (mod p).
    Verwendet Wilson's Theorem: (p-1)! ≡ -1 (mod p)
    """
    violations = []
    for i in range(len(primes_arr)):
        p = primes_arr[i]
        # Berechne !p = Σ_{k=0}^{p-1} k! mod p
        # Inkrementell: fak = k! mod p, s = laufende Summe
        fak = 1  # 0! = 1
        s = 1    # Summe beginnt bei k=0: 0!=1
        for k in range(1, p):
            fak = fak * k % p
            s = (s + fak) % p
        if s == 0:
            violations.append(p)
    return violations


def check_chunk(args):
    """Verarbeitet einen Chunk von Primzahlen (für multiprocessing)."""
    chunk_primes, chunk_id = args
    arr = np.array(chunk_primes, dtype=np.int64)
    start = time.time()
    result = kurepa_check_range(arr)
    elapsed = time.time() - start
    print(f"  Chunk {chunk_id:3d}: {len(chunk_primes):8,} Primzahlen "
          f"({chunk_primes[0]:,}–{chunk_primes[-1]:,}) — {elapsed:.1f}s — "
          f"Verletzungen: {len(result)}", flush=True)
    return list(result)


def main():
    LIMIT = 10_000_000
    CHUNK_SIZE = 5_000  # Primzahlen pro Chunk
    N_WORKERS = min(20, multiprocessing.cpu_count())
    OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "kurepa_result.json")

    print(f"{'='*60}")
    print(f"Kurepa Heavy Computation — Start: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Grenze: p ≤ {LIMIT:,} | Kerne: {N_WORKERS} | Chunk: {CHUNK_SIZE:,}")
    print(f"{'='*60}")

    # Alle Primzahlen bis LIMIT generieren
    print("Generiere Primzahlen...", flush=True)
    t0 = time.time()
    primes = list(primerange(3, LIMIT + 1))  # p=2: !2=1≠0, trivial
    print(f"  {len(primes):,} Primzahlen bis {LIMIT:,} ({time.time()-t0:.1f}s)")

    # Warm-up: Numba JIT kompilieren
    print("Numba JIT warm-up...", flush=True)
    _ = kurepa_check_range(np.array([5, 7, 11], dtype=np.int64))
    print("  JIT bereit.")

    # Chunks aufteilen
    chunks = []
    for i in range(0, len(primes), CHUNK_SIZE):
        chunk = primes[i:i + CHUNK_SIZE]
        chunks.append((chunk, i // CHUNK_SIZE))

    print(f"\nVerarbeite {len(chunks):,} Chunks mit {N_WORKERS} Kernen...\n")
    t1 = time.time()

    all_violations = []
    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(check_chunk, chunks)
    for r in results:
        all_violations.extend(r)

    elapsed = time.time() - t1
    print(f"\n{'='*60}")
    print(f"Fertig: {elapsed:.0f}s ({elapsed/3600:.2f}h)")
    print(f"Geprüfte Primzahlen: {len(primes):,}")
    print(f"Verletzungen (Gegenbeispiele): {len(all_violations)}")
    if all_violations:
        print(f"ACHTUNG — Gegenbeispiele gefunden: {all_violations[:10]}")
    else:
        print("BESTÄTIGT: Kein Gegenbeispiel bis p =", LIMIT)

    # Ergebnis speichern
    result = {
        "conjecture": "Kurepa",
        "limit": LIMIT,
        "primes_checked": len(primes),
        "violations": all_violations,
        "confirmed": len(all_violations) == 0,
        "elapsed_seconds": round(elapsed, 1),
        "timestamp": datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis gespeichert: {OUTPUT_FILE}")
    return result


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()
