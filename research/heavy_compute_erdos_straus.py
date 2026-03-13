"""
@file heavy_compute_erdos_straus.py
@brief Erdős-Straus: 4/p = 1/x + 1/y + 1/z — Verifikation bis n = 10^9
@description
    Erdős-Straus-Vermutung (1948): Für alle n ≥ 2 existieren positive
    ganze Zahlen x,y,z mit 4/n = 1/x + 1/y + 1/z.

    Algorithmus (O(√p) Worst-Case pro Primzahl):

    Schlüsselbeobachtung: 4/p = 1/x + n/D mit x = ceil(p/4), n = 4x-p, D = px.
    Für 1/y + 1/z = n/D gilt: (ny - D)(nz - D) = D².
    Sei d ein Teiler von D² (= p²x²), dann:
        y = (D + d) / n   [ganzzahlig gdw. n | (D + d)]
        z = (D + D²/d) / n = (D/d · (D + d)) / n
    Es reicht, Teiler d von D² mit n|(D+d) zu finden.

    Zerlegung d = p^a · e mit a ∈ {0,1,2} und e | x²:

    WICHTIG (Bug-Erkenntnis Build 179):
    Nicht nur Teiler von D=px, sondern von D²=p²x² müssen geprüft werden!
    Gegenbeispiel der alten Version: p=2521, d=848=2⁴×53.
        x=636=2²×3×53 → 848|x²=2⁴×9×53², aber 848∤x (nur 2²|x, nicht 2⁴).
    Fix: vollständige Enumeration aller Teiler von x² via Faktorisierung.

    Beweis der Korrektheit für z:
        d|D² → D²/d ganzzahlig. n|(D+d) → D≡−d (mod n) → D²≡d² (mod n).
        Also D+D²/d ≡ −d+d = 0 (mod n) → z = (D+D²/d)/n ganzzahlig. ✓

    Laufzeitanalyse:
        - p ≡ 3 mod 4: O(1) via bewiesene Formel
        - p ≡ 1 mod 4: O(√p) via Teilersuche
          Faktorisierung von x: O(√x)=O(√p); Teileranzahl d(x²): typisch <500

    Verifiziert: 0 Gegenbeispiele bis 1.000.000 (Python) und bis 10^9 (Ziel).

@author Michael Fuhrmann
@date 2026-03-13
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

@numba.njit(cache=False)
def find_straus_decomp(p):
    """
    Findet x,y,z > 0 mit 4/p = 1/x + 1/y + 1/z für Primzahl p.
    Gibt 1 zurück wenn gefunden, 0 wenn nicht (Gegenbeispiel!).

    Algorithmus: Suche x ab ceil(p/4). Für jedes x:
      - Berechne n=4x-p, D=px.
      - Faktorisiere x via Trial Division.
      - Erzeuge ALLE Teiler von x² (nicht nur x!).
      - Für a=0: prüfe alle e|x²  (d = e)
      - Für a=1: prüfe alle t|x   (d = p·t, kein Überlauf da t≤x)
      - Für a=2: prüfe e=1         (d = p², kleiner Spezialfall)

    @param p  Primzahl ≥ 2
    @return   1 wenn Lösung existiert, 0 sonst
    @lastModified 2026-03-13
    """
    # ── Spezialfälle ────────────────────────────────────────────────────
    if p == 2:
        return 1  # 4/2 = 1/1 + 1/2 + 1/2 ✓
    if p == 3:
        return 1  # 4/3 = 1/1 + 1/4 + 1/12 ✓

    # ── p ≡ 3 mod 4: bewiesene O(1)-Formel ─────────────────────────────
    # Sei k=(p+1)//4: 4/p = 1/k + 1/(2kp) + 1/(2kp) ✓
    if p % 4 == 3:
        return 1

    # ── p ≡ 1 mod 4: Teilersuche auf D² = p²x² ─────────────────────────
    x_min = (p + 3) // 4   # = ceil(p/4), ganzzahlig da p≡1 mod 4

    for dx in range(300):   # In Praxis reichen fast immer 1-5 Iterationen
        x = x_min + dx
        n = 4 * x - p       # Rest-Zähler: n ≥ 3, wächst um 4 pro Schritt
        D = p * x           # Rest-Nenner: D=px (int64-sicher für p≤10^9)

        # ── Faktorisiere x via Trial Division ───────────────────────────
        # Speichere Primfaktoren und Exponenten in fixen Arrays (numba-kompatibel)
        pf = np.zeros(20, dtype=np.int64)   # Primfaktoren von x
        pe = np.zeros(20, dtype=np.int64)   # Exponenten in x
        nf = np.int64(0)                    # Anzahl verschiedener Primfaktoren

        temp = x
        f = np.int64(2)
        while f * f <= temp:
            if temp % f == 0:
                e = np.int64(0)
                while temp % f == 0:
                    e += 1
                    temp //= f
                pf[nf] = f
                pe[nf] = e
                nf += 1
            f += 1
        if temp > 1:
            pf[nf] = temp
            pe[nf] = np.int64(1)
            nf += 1

        # ── Erzeuge alle Teiler von x² ──────────────────────────────────
        # x² = ∏ pf[i]^(2·pe[i]). Teiler: alle ∏ pf[i]^k_i mit 0≤k_i≤2·pe[i].
        # Algorithmus: starte mit [1], multipliziere sukzessive mit pf[i]^1..^(2e).
        # Maximale Teileranzahl von x² für x≤2.5×10^8: typisch <1000, max ~20000.
        MAX_DIVS = np.int64(3000)
        divs = np.zeros(3000, dtype=np.int64)
        divs[0] = np.int64(1)
        nd = np.int64(1)

        for i in range(nf):
            old_nd = nd       # Anzahl Teiler vor diesem Primfaktor
            pk = pf[i]        # Startet bei pf[i]^1
            for k in range(1, 2 * pe[i] + 1):
                for j in range(old_nd):
                    if nd < MAX_DIVS:
                        divs[nd] = divs[j] * pk
                        nd += 1
                pk *= pf[i]   # pf[i]^k → pf[i]^(k+1)

        # ── a=0: prüfe alle Teiler e von x² als d = e ───────────────────
        # d|x² ⊆ D²=p²x². Kein Überlaufrisiko: e≤x²≤(2.5×10^8)²=6×10^16<int64max.
        for j in range(nd):
            if (D + divs[j]) % n == 0:
                return 1

        # ── a=1: prüfe alle Teiler t von x als d = p·t ──────────────────
        # d=p·t|D²=p²x² (da t|x). Kein Überlauf: p·t ≤ p·x = D ≤ int64max/2.
        t = np.int64(1)
        while t * t <= x:
            if x % t == 0:
                t2 = x // t
                if (D + p * t) % n == 0:
                    return 1
                if (D + p * t2) % n == 0:
                    return 1
            t += 1

        # ── a=2, e=1: d = p² ────────────────────────────────────────────
        # d=p²|D²=p²x². p²+D ≤ 10^18 + 2.5×10^17 < int64max für p≤10^9.
        if (D + p * p) % n == 0:
            return 1

    # Sollte nie erreicht werden für Primzahlen ≤ 10^9
    return 0


@numba.njit(cache=False)
def check_prime_chunk(primes_arr):
    """
    Prüft Erdős-Straus für einen Array von Primzahlen.

    @param primes_arr  Array von int64-Primzahlen
    @return            Liste von Gegenbeispielen (leer wenn Vermutung gilt)
    @lastModified 2026-03-13
    """
    violations = []
    for i in range(len(primes_arr)):
        p = primes_arr[i]
        if find_straus_decomp(p) == 0:
            violations.append(p)
    return violations


def process_chunk(args):
    """
    Verarbeitet einen Chunk von Primzahlen mit Zeitnahme.

    @param args   Tupel (chunk: list[int], chunk_id: int)
    @return       Liste von Gegenbeispielen im Chunk
    @lastModified 2026-03-13
    """
    chunk, chunk_id = args
    arr = np.array(chunk, dtype=np.int64)
    t = time.time()
    result = check_prime_chunk(arr)
    violations = list(result)
    elapsed = time.time() - t

    # Fortschrittsausgabe alle 100 Chunks
    if chunk_id % 100 == 0:
        rate = len(chunk) / elapsed if elapsed > 0 else 0
        print(f"  Chunk {chunk_id:5d}: {len(chunk):8,} Primzahlen "
              f"({chunk[0]:,}–{chunk[-1]:,}) {elapsed:.2f}s "
              f"({rate/1e6:.1f}M prim/s)", flush=True)
    return violations


def main():
    """
    Hauptfunktion: Erdős-Straus für alle Primzahlen bis 10^9 verifizieren.

    @lastModified 2026-03-13
    """
    LIMIT       = 10**9
    CHUNK_SIZE  = 100_000      # Primzahlen pro Worker-Chunk
    N_WORKERS   = min(20, multiprocessing.cpu_count())
    SEGMENT     = 50_000_000   # Zahlenbereich pro Segment (Primzahl-Sieb)

    print(f"{'='*60}")
    print(f"Erdős-Straus Heavy Computation — {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Grenze: n ≤ {LIMIT:,} | Kerne: {N_WORKERS}")
    print(f"Algorithmus: O(√p) Teilersuche auf D²=p²x² (Build 179)")
    print(f"{'='*60}\n")

    # ── JIT-Warm-up ─────────────────────────────────────────────────────
    print("Numba JIT kompilieren...", flush=True)
    t_jit = time.time()
    # Warm-up inkl. bekannter Problemfälle
    warmup = np.array([2, 3, 5, 7, 11, 13, 17, 97, 193, 241,
                       2521, 1000003, 999999937], dtype=np.int64)
    result_warmup = check_prime_chunk(warmup)
    if list(result_warmup):
        print(f"FEHLER: Gegenbeispiele im Warm-up: {list(result_warmup)}")
        return
    print(f"Numba JIT ready. ({time.time()-t_jit:.1f}s)\n")

    all_violations = []
    chunk_id       = 0
    total_primes   = 0
    total_elapsed  = 0.0
    t0             = time.time()

    for seg_start in range(2, LIMIT + 1, SEGMENT):
        seg_end = min(seg_start + SEGMENT - 1, LIMIT)
        print(f"Segment {seg_start:,}–{seg_end:,}...", flush=True)
        t_seg = time.time()

        # Primzahlen im Segment sieben
        primes_seg  = list(primerange(seg_start, seg_end + 1))
        total_primes += len(primes_seg)

        # In Chunks aufteilen
        chunks = [(primes_seg[i:i+CHUNK_SIZE], chunk_id + i//CHUNK_SIZE)
                  for i in range(0, len(primes_seg), CHUNK_SIZE)]
        chunk_id += len(chunks)

        # Parallel auf Worker-Pool verteilen
        with multiprocessing.Pool(N_WORKERS) as pool:
            results = pool.map(process_chunk, chunks)

        for r in results:
            all_violations.extend(r)

        seg_time       = time.time() - t_seg
        total_elapsed += seg_time
        prim_rate      = len(primes_seg) / seg_time if seg_time > 0 else 0

        print(f"  → {len(primes_seg):,} Primzahlen, {seg_time:.1f}s "
              f"({prim_rate/1e6:.1f}M prim/s), "
              f"Gesamt: {total_elapsed:.1f}s ({total_elapsed/3600:.2f}h)", flush=True)

        # ETA-Schätzung
        progress = seg_end / LIMIT
        if 0 < progress < 1.0:
            eta = total_elapsed / progress * (1.0 - progress)
            print(f"  → Fortschritt: {progress:.1%}, "
                  f"ETA: {eta:.0f}s ({eta/3600:.2f}h)\n", flush=True)

    # ── Ergebnis ────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Fertig: {elapsed:.0f}s ({elapsed/3600:.2f}h)")
    print(f"Geprüfte Primzahlen bis {LIMIT:,}: {total_primes:,}")
    print(f"Gegenbeispiele: {len(all_violations)}")

    if not all_violations:
        print(f"✓ BESTÄTIGT: Erdős-Straus gilt für alle Primzahlen ≤ {LIMIT:,}")
    else:
        print(f"✗ GEGENBEISPIELE GEFUNDEN: {all_violations[:10]}")

    result = {
        "conjecture":      "Erdos-Straus",
        "limit":           LIMIT,
        "primes_checked":  total_primes,
        "violations":      all_violations[:100],
        "confirmed":       len(all_violations) == 0,
        "elapsed_seconds": round(elapsed, 1),
        "algorithm":       "O(sqrt(p)) divisor search on D^2=p^2*x^2 (Build 179)",
        "timestamp":       datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis gespeichert: {OUTPUT_FILE}")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()
