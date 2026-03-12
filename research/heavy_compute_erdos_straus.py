"""
@file heavy_compute_erdos_straus.py
@brief Erdős-Straus: 4/n = 1/x + 1/y + 1/z — Verifikation bis n = 10^9
@description
    Erdős-Straus-Vermutung (1948): Für alle n ≥ 2 existieren positive
    ganze Zahlen x,y,z mit 4/n = 1/x + 1/y + 1/z.

    Strategie (für Primzahlen n=p):
    - n ≡ 1 (mod 4): 4/p = 1/((p+3)/4) + ... (parametrische Formeln)
    - n ≡ 3 (mod 4): 4/p = 1/((p+1)/4) + 1/((p+1)/4·p) [BEWIESEN in erdos_straus_ext.py]
    - n ≡ 1 (mod 3): 4/p = 1/((p+2)/3) + ...
    - Fallback: Brute-Force Suche
    Für Nicht-Primzahlen: meist trivial (nutze Teiler)

    Ziel: Alle n bis 10^9 verifizieren (nur Primzahlen sind kritisch)
    Multiprocessing mit 20 Kernen.

@author Michael Fuhrmann
@date 2026-03-12
"""

import os
import time
import json
import multiprocessing
from datetime import datetime

import numba
import numpy as np
from sympy import primerange

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "erdos_straus_result.json")

# ─── Numba-beschleunigte Prüfung ─────────────────────────────────────────────

@numba.njit(cache=True)
def find_straus_decomp(p):
    """
    Findet x,y,z > 0 mit 4/p = 1/x + 1/y + 1/z für Primzahl p.
    Gibt 1 zurück wenn gefunden, 0 wenn nicht (Gegenbeispiel!).
    """
    # Parametrische Formeln (schnell)

    # p ≡ 3 mod 4: 4/p = 1/((p+1)/4) + 1/((p+1)/4·p)  (2-Term → 3-Term mit z=0 trivial)
    # Das ist schon ein Beweis — überspringen wir es hier für Geschwindigkeit

    # p ≡ 1 mod 4: 4/p = 1/((p+3)/4) + ...
    if p % 4 == 3:
        # Bewiesen: 4/p = 1/k + 1/(kp) mit k = (p+1)/4
        # Das ist immer gültig, also trivial ✓
        return 1

    # p ≡ 0 mod 3 (also p=3):
    if p == 3:
        return 1  # 4/3 = 1/1 + 1/3 + 1/... nein, 4/3 > 1... 4/3=1+1/3=1/1+1/3, brauchen 3 Terme: 1/1+1/4+1/12

    # p ≡ 1 mod 3: 4/p = 1/((p+2)/3) + 1/((p+2)p/3) wenn (p+2)%3==0
    if (p + 2) % 3 == 0:
        k = (p + 2) // 3
        # 4/p = 1/k + 3/(kp) = 1/k + 1/(kp/3) wenn 3|k... prüfe direkt
        # 1/k + 1/((kp+1)/... hmm, lass uns rechnerisch prüfen
        # 4/p = 1/k + r wobei r = 4/p - 1/k = (4k-p)/(pk)
        # Für k=(p+2)/3: 4k-p = 4(p+2)/3 - p = (4p+8-3p)/3 = (p+8)/3
        # r = (p+8)/(3pk)
        pass

    # p ≡ 1 mod 4: 4/p = 1/((p+3)/4) + 1/p·(p+3)/4/...
    # Systematische Formel: für p ≡ 1 mod 4, sei k=(p+3)/4: 4/p=1/k + (4k-p)/(pk)
    # = 1/k + 1/(pk/(4k-p)) wenn (4k-p) | pk
    if (p + 3) % 4 == 0:
        k = (p + 3) // 4
        # 4k - p = p+3 - p = 3... Nein: k = (p+3)/4, 4k = p+3, 4k-p = 3
        # r = 3/(pk) → 4/p = 1/k + 3/(pk) = 1/k + 1/(pk/3) wenn 3|pk
        # pk = p(p+3)/4. Wenn 3|p(p+3): 3|p oder 3|(p+3)
        # Wenn 3|p: p=3 (Primzahl) → separat. Wenn 3|(p+3): p≡0(mod3) → p=3
        # Also für p≡1(mod4) und p≢3(mod3): 4/p = 1/((p+3)/4) + ???
        # Brauchen weiteren Term. Prüfe Formel: 4/p = 1/p + 3/p = 1/p + 1/(p/3+..)
        pass

    # Brute-Force für restliche Fälle (sollte selten sein)
    # 4/p = 1/x + 1/y + 1/z, x ≤ y ≤ z, x ≥ ceil(3/4·p) = ceil(p·3/4)
    # x ≥ ceil(p/4)
    x_min = (p + 3) // 4  # ceil(p/4)
    x_max = p              # 1/x ≥ 1/p, also x ≤ p (da 4/p≥1/p trivial, aber 1/x≤4/p→x≥p/4)

    for x in range(x_min, x_max + 1):
        # Rest: 4/p - 1/x = (4x - p)/(px)
        num = 4 * x - p
        den = p * x
        if num <= 0:
            continue
        # Jetzt 1/y + 1/z = num/den, y ≤ z, y ≥ ceil(2·den/num/2) = ceil(2·den/num... hmm
        # y ≥ ceil(den/num) (da 1/y ≤ num/den → y ≥ den/num)
        # Genauer: 1/y ≥ num/(2·den) → y ≤ 2·den/num
        y_min = (den + num - 1) // num   # ceil(den/num)
        y_max = 2 * den // num

        for y in range(y_min, y_max + 1):
            # z: 1/z = num/den - 1/y = (num·y - den)/(den·y)
            rnum = num * y - den
            rden = den * y
            if rnum <= 0:
                continue
            if rden % rnum == 0:
                return 1  # Gefunden!

    return 0  # Kein Gegenbeispiel für diese p (sollte nie passieren)


@numba.njit(cache=True)
def check_prime_chunk(primes_arr):
    """Prüft Erdős-Straus für einen Array von Primzahlen."""
    violations = []
    for i in range(len(primes_arr)):
        p = primes_arr[i]
        if find_straus_decomp(p) == 0:
            violations.append(p)
    return violations


def process_chunk(args):
    chunk, chunk_id = args
    arr = np.array(chunk, dtype=np.int64)
    t = time.time()
    result = check_prime_chunk(arr)
    violations = list(result)
    elapsed = time.time() - t
    if chunk_id % 100 == 0:
        print(f"  Chunk {chunk_id:5d}: {len(chunk):8,} Primzahlen "
              f"({chunk[0]:,}–{chunk[-1]:,}) {elapsed:.1f}s", flush=True)
    return violations


def main():
    LIMIT = 10**9
    CHUNK_SIZE = 50_000
    N_WORKERS = min(20, multiprocessing.cpu_count())

    print(f"{'='*60}")
    print(f"Erdős-Straus Heavy Computation — {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Grenze: n ≤ {LIMIT:,} | Kerne: {N_WORKERS}")
    print(f"{'='*60}\n")

    # Nur Primzahlen sind kritisch (andere sind trivial über Teiler)
    print("Generiere Primzahlen bis 10^9...", flush=True)
    t0 = time.time()
    # Sieve of Eratosthenes für Primzahlen bis 10^9 ist zu viel RAM
    # Wir verarbeiten in Segmenten
    SEGMENT = 50_000_000  # 50M Primzahlen auf einmal

    all_violations = []
    chunk_id = 0
    total_primes = 0
    total_elapsed = 0

    # Warm-up JIT
    _ = check_prime_chunk(np.array([5, 7, 11, 13], dtype=np.int64))
    print("Numba JIT ready.\n")

    for seg_start in range(2, LIMIT + 1, SEGMENT):
        seg_end = min(seg_start + SEGMENT - 1, LIMIT)
        print(f"Segment {seg_start:,}–{seg_end:,}...", flush=True)
        t_seg = time.time()

        primes_seg = list(primerange(seg_start, seg_end + 1))
        total_primes += len(primes_seg)

        # Chunks
        chunks = [(primes_seg[i:i+CHUNK_SIZE], chunk_id + i//CHUNK_SIZE)
                  for i in range(0, len(primes_seg), CHUNK_SIZE)]
        chunk_id += len(chunks)

        with multiprocessing.Pool(N_WORKERS) as pool:
            results = pool.map(process_chunk, chunks)

        for r in results:
            all_violations.extend(r)

        seg_time = time.time() - t_seg
        total_elapsed += seg_time
        print(f"  → {len(primes_seg):,} Primzahlen, {seg_time:.0f}s, "
              f"Gesamtzeit: {total_elapsed:.0f}s ({total_elapsed/3600:.2f}h)", flush=True)

        # Schätze verbleibende Zeit
        progress = seg_end / LIMIT
        if progress > 0:
            eta = total_elapsed / progress * (1 - progress)
            print(f"  → Fortschritt: {progress:.1%}, ETA: {eta:.0f}s ({eta/3600:.2f}h)\n",
                  flush=True)

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Fertig: {elapsed:.0f}s ({elapsed/3600:.2f}h)")
    print(f"Geprüfte Primzahlen bis {LIMIT:,}: {total_primes:,}")
    print(f"Gegenbeispiele: {len(all_violations)}")
    if not all_violations:
        print(f"BESTÄTIGT: Erdős-Straus gilt für alle Primzahlen ≤ {LIMIT:,}")

    result = {
        "conjecture": "Erdos-Straus",
        "limit": LIMIT,
        "primes_checked": total_primes,
        "violations": all_violations[:100],
        "confirmed": len(all_violations) == 0,
        "elapsed_seconds": round(elapsed, 1),
        "timestamp": datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis gespeichert: {OUTPUT_FILE}")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()
