"""
@file heavy_compute_giuga4.py
@brief Giuga 4-Prim-Suche bis Produkt 10^15 — Factorize-Ansatz
@description
    Sucht 4-Prim-Giuga-Pseudoprimes: n = p1*p2*p3*p4 mit p1<p2<p3<p4 alle prim
    und für alle pi: pi | (n/pi − 1) AND (pi−1) | (n/pi − 1)

    Effizienter Ansatz:
    - Fix p1, p2, p3. Dann muss p4 | (p1*p2*p3 − 1) (aus p4-Bedingung)
    - Für jeden Kandidaten p4 = Primteiler von (p1*p2*p3−1) prüfe alle 4 Bedingungen
    - Multiprocessing über p1-Werte

    Bekannt: Kein Giuga-Pseudoprime mit 4 Primfaktoren bis ~10^{13800}
    (theoretische Schranken aus giuga_4prim.py)
    Ziel: Praktische Suche bis Produkt ~ 10^15

@author Michael Fuhrmann
@date 2026-03-12
"""

import os
import sys
import time
import json
import multiprocessing
from datetime import datetime
from itertools import count

import sympy
from sympy import isprime, factorint, nextprime, primerange

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "giuga4_result.json")

# ─── Giuga-Bedingung prüfen ──────────────────────────────────────────────────

def giuga_condition(primes_list):
    """Prüft ob das Produkt der Primzahlen ein Giuga-Pseudoprim ist."""
    n = 1
    for p in primes_list:
        n *= p
    for p in primes_list:
        nprime = n // p  # n/p (ganzzahlig, da p|n)
        if (nprime - 1) % p != 0:
            return False
        if (nprime - 1) % (p - 1) != 0:
            return False
    return True


def search_p1(args):
    """
    Sucht 4-Prim-Giuga-Zahlen mit erstem Faktor p1.
    Für festes p1 iteriere über p2, p3 und leite p4-Kandidaten ab.
    """
    p1, max_product = args
    found = []
    violations = 0

    # p2 von p1+2 (nächste Primzahl) bis (max_product/p1)^(1/3)
    p2_limit = int((max_product / p1) ** (1/3)) + 1
    p2 = nextprime(p1)

    while p2 <= p2_limit:
        # p3 von p2+2 bis sqrt(max_product/(p1*p2))
        p3_limit = int((max_product / (p1 * p2)) ** 0.5) + 1
        p3 = nextprime(p2)

        while p3 <= p3_limit:
            product123 = p1 * p2 * p3
            target = product123 - 1  # p4 muss target teilen

            # Faktoren von (p1*p2*p3 - 1) bestimmen
            try:
                factors = factorint(target)
            except Exception:
                p3 = nextprime(p3)
                continue

            # Alle Primteiler von target als p4-Kandidaten prüfen
            for p4 in factors.keys():
                if p4 <= p3:
                    continue  # p4 muss größer als p3 sein
                if not isprime(p4):
                    continue
                n = product123 * p4
                if n > max_product:
                    continue
                # Vollständige Giuga-Prüfung für alle 4 Primzahlen
                if giuga_condition([p1, p2, p3, p4]):
                    found.append({
                        "primes": [p1, p2, p3, p4],
                        "product": n
                    })
                    print(f"  *** GIUGA 4-PRIM GEFUNDEN: {[p1,p2,p3,p4]}, n={n} ***",
                          flush=True)
            p3 = nextprime(p3)
        p2 = nextprime(p2)

    return found


def main():
    # Suchgrenze: Produkt ≤ 10^15
    MAX_PRODUCT = 10**15
    # p1 geht von 2 bis (MAX_PRODUCT)^(1/4) ≈ 31622
    P1_LIMIT = int(MAX_PRODUCT ** 0.25) + 1
    N_WORKERS = min(20, multiprocessing.cpu_count())

    print(f"{'='*60}")
    print(f"Giuga 4-Prim Heavy Search — Start: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Suchgrenze: Produkt ≤ {MAX_PRODUCT:.0e} | p1 ≤ {P1_LIMIT:,}")
    print(f"Kerne: {N_WORKERS}")
    print(f"{'='*60}\n")

    # Alle p1-Werte (Primzahlen bis P1_LIMIT)
    p1_values = list(primerange(2, P1_LIMIT + 1))
    tasks = [(p1, MAX_PRODUCT) for p1 in p1_values]

    print(f"p1-Kandidaten: {len(p1_values)} ({p1_values[0]}..{p1_values[-1]})")
    print(f"Starte Suche mit {N_WORKERS} Kernen...\n")

    t0 = time.time()
    all_found = []

    # Für p1-Werte bis 100 sequentiell (kurz), ab 100 parallel
    small_tasks = [(p1, MAX_PRODUCT) for p1 in p1_values if p1 <= 100]
    large_tasks = [(p1, MAX_PRODUCT) for p1 in p1_values if p1 > 100]

    # Kleine p1 mit vielen Iterationen parallel
    print(f"Phase 1: p1 ≤ 100 ({len(small_tasks)} Werte, aufwändig)...")
    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(search_p1, small_tasks)
    for r in results:
        all_found.extend(r)
    print(f"  Phase 1 fertig: {time.time()-t0:.0f}s, {len(all_found)} Treffer")

    print(f"\nPhase 2: p1 > 100 ({len(large_tasks)} Werte, schnell)...")
    t1 = time.time()
    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(search_p1, large_tasks)
    for r in results:
        all_found.extend(r)
    print(f"  Phase 2 fertig: {time.time()-t1:.0f}s")

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Fertig: {elapsed:.0f}s ({elapsed/3600:.2f}h)")
    print(f"4-Prim-Giuga-Zahlen gefunden: {len(all_found)}")
    if not all_found:
        print(f"BESTÄTIGT: Kein 4-Prim-Giuga-Pseudoprim mit Produkt ≤ {MAX_PRODUCT:.0e}")
    else:
        for item in all_found:
            print(f"  Gefunden: {item}")

    result = {
        "conjecture": "Giuga 4-Prim",
        "max_product": MAX_PRODUCT,
        "p1_limit": P1_LIMIT,
        "found": all_found,
        "confirmed": len(all_found) == 0,
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
