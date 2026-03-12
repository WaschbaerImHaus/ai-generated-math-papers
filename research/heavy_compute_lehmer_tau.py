"""
@file heavy_compute_lehmer_tau.py
@brief Lehmer τ-Vermutung: τ(n) ≠ 0 bis n = 50_000 — Python BigInt korrekt
@description
    Ramanujan τ-Funktion via Δ(q) = q·∏_{n≥1}(1−q^n)^24 = Σ τ(n)q^n.
    Vermutung (Lehmer 1947): τ(n) ≠ 0 für alle n ≥ 1.

    KORREKTE Implementierung mit Python-BigInts (beliebige Genauigkeit).
    τ-Werte wachsen als ~n^{11/2} (Deligne-Schranke: |τ(p)| ≤ 2p^{11/2}).
    Für n=1000: |τ(n)| ≈ 10^{33} → int64 NICHT ausreichend, Python int nötig.

    Strategie:
    1. Exakt für n ≤ 5000 mit Python-Listen (BigInt, Konvolutionsmethode)
    2. Modular-Check für n ≤ 50000: τ(n) mod p₁, ..., τ(n) mod pₖ berechnen
       Falls alle ≡ 0 → mögliche Nullstelle (dann exakt prüfen)

    Alternativ: sympy.ntheory.ram_sum oder direkte Formel.

@author Michael Fuhrmann
@date 2026-03-12
"""

import os
import time
import json
from datetime import datetime

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "lehmer_tau_result.json")

# Bekannte τ-Werte zur Verifikation (via Multiplikativität verifiziert)
# τ(100) = τ(4)*τ(25) = (-1472)*(-25499225) = 37534859200
# τ(1000) = τ(8)*τ(125) = 84480*(-359001100500) = -30328412970240000
KNOWN_TAU = {
    1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830,
    6: -6048, 7: -16744, 8: 84480, 9: -113643,
    10: -115920, 11: 534612, 12: -370944,
    100: 37534859200, 1000: -30328412970240000
}


def compute_tau_bigint(n_max):
    """
    Berechnet τ(1)...τ(n_max) mit Python-BigInts.
    Methode: Koeffizientenentwicklung von Δ(q) = q·∏(1-q^k)^24

    Repräsentiere Polynom als Python-Liste von BigInts.
    """
    print(f"  Berechne Δ(q)-Koeffizienten bis n={n_max:,} (Python BigInt)...",
          flush=True)
    t = time.time()

    size = n_max + 1
    # coeffs[i] = Koeffizient bei q^i in ∏_{k=1}^{n_max}(1-q^k)^24
    # Starte mit 1 (leeres Produkt)
    coeffs = [0] * size
    coeffs[0] = 1

    # Multipliziere inkrementell mit (1-q^k)^24 für k=1..n_max
    for k in range(1, n_max + 1):
        if k % 500 == 0:
            print(f"    k={k:,}/{n_max:,} ({time.time()-t:.0f}s)", flush=True)

        # Binomialkoeffizienten C(24,j) mit Vorzeichen
        # (1-x)^24 = Σ_{j=0}^{24} C(24,j)(-1)^j x^j
        binom24 = [
            1, -24, 276, -2024, 10626, -42504, 134596, -346104,
            735471, -1307504, 1961256, -2496144, 2704156, -2496144,
            1961256, -1307504, 735471, -346104, 134596, -42504,
            10626, -2024, 276, -24, 1
        ]

        # Neue Koeffizienten: new[i] = Σ_{j≥0} binom24[j] * coeffs[i - j*k]
        # Berechne in-place von hinten nach vorne
        for i in range(size - 1, -1, -1):
            s = coeffs[i]  # j=0 Term
            for j in range(1, 25):
                idx = i - j * k
                if idx < 0:
                    break
                s += binom24[j] * coeffs[idx]
            coeffs[i] = s

    # Δ(q) = q · (Ergebnis) → τ(n) = coeffs[n-1] für n ≥ 1
    tau = [coeffs[n - 1] for n in range(1, n_max + 1)]
    print(f"  Fertig: {time.time()-t:.1f}s")
    return tau


def main():
    N_MAX = 50_000  # ~85 min (O(N²) skaliert: 5000→50000 = 100× länger)

    print(f"{'='*60}")
    print(f"Lehmer τ Heavy Computation — {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Ziel: τ(n) ≠ 0 für alle n ≤ {N_MAX:,} (exakt, BigInt)")
    print(f"{'='*60}\n")

    t0 = time.time()
    tau = compute_tau_bigint(N_MAX)
    elapsed_compute = time.time() - t0

    # Bekannte Werte verifizieren
    print("\nVerifiziere bekannte Werte:")
    all_correct = True
    for n, expected in sorted(KNOWN_TAU.items()):
        if n <= N_MAX:
            got = tau[n - 1]
            ok = got == expected
            print(f"  τ({n:4d}) = {got} {'✓' if ok else f'✗ erwartet {expected}'}")
            if not ok:
                all_correct = False

    # Nullstellen suchen
    zeros = [n for n in range(1, N_MAX + 1) if tau[n - 1] == 0]

    # Deligne-Schranken prüfen: |τ(p)| ≤ 2·p^{11/2} für Primzahlen p
    from sympy import isprime
    deligne_violations = []
    for n in range(2, min(N_MAX + 1, 1000)):
        if isprime(n):
            bound = 2 * n**5.5
            if abs(tau[n - 1]) > bound:
                deligne_violations.append(n)

    print(f"\n{'='*60}")
    print(f"Ergebnis:")
    print(f"  Berechnungszeit: {elapsed_compute:.0f}s")
    print(f"  Geprüfte n: 1..{N_MAX:,}")
    print(f"  Bekannte Werte korrekt: {all_correct}")
    print(f"  Nullstellen: {zeros}")
    print(f"  Deligne-Verletzungen: {deligne_violations}")
    if not zeros:
        print(f"BESTÄTIGT: τ(n) ≠ 0 für alle n ≤ {N_MAX:,}")

    result = {
        "conjecture": "Lehmer-tau",
        "n_max": N_MAX,
        "method": "BigInt polynomial expansion",
        "zeros": zeros,
        "deligne_violations": deligne_violations,
        "known_values_correct": all_correct,
        "tau_sample": {str(n): tau[n-1] for n in [1,2,3,4,5,6,7,8,9,10,12,100,1000]
                       if n <= N_MAX},
        "elapsed_seconds": round(elapsed_compute, 1),
        "timestamp": datetime.now().isoformat(),
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Ergebnis: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
