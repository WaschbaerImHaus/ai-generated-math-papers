"""
numba_jit.py — Numba-beschleunigte Kernfunktionen für specialist-maths.

Enthält JIT-kompilierte Versionen der numerisch intensivsten Funktionen.
Falls Numba nicht verfügbar ist, werden transparente Python-Fallbacks genutzt.

Numba kompiliert beim ERSTEN Aufruf (~1-3 Sek), dann ist die Funktion schnell.

@author: Michael Fuhrmann
@date: 2026-03-11
@lastModified: 2026-03-11
"""
import numpy as np

# Numba-Import mit Fallback-Dekorator
try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    # Fallback: Dekoratoren, die nichts tun — API bleibt identisch
    def njit(*args, **kwargs):
        """Dummy-Dekorator wenn Numba nicht installiert ist."""
        def decorator(func):
            return func
        return decorator if args and callable(args[0]) else decorator
    prange = range
    NUMBA_AVAILABLE = False


@njit(cache=True)
def sieve_numpy(limit: int) -> np.ndarray:
    """
    Sieb des Eratosthenes via Numba-JIT — NumPy boolean Array.

    Zeitkomplexität: O(n log log n), aber 10–50× schneller als Python-Schleife
    durch Maschinencode-Kompilierung und Cache-freundlichen Speicherzugriff.
    Gibt ein boolean Array zurück: result[i] = True bedeutet i ist prim.

    @param limit: Obere Grenze (inklusiv)
    @return: numpy boolean Array der Länge limit+1
    @author: Michael Fuhrmann
    @date: 2026-03-11
    @lastModified: 2026-03-11
    """
    # Boolean-Array initialisieren: alle als prim annehmen
    is_prime = np.ones(limit + 1, dtype=np.bool_)
    is_prime[0] = False
    if limit >= 1:
        is_prime[1] = False

    # Sieben: alle Vielfachen als nicht-prim markieren
    p = 2
    while p * p <= limit:
        if is_prime[p]:
            # Alle Vielfachen ab p² als nicht-prim markieren
            multiple = p * p
            while multiple <= limit:
                is_prime[multiple] = False
                multiple += p
        p += 1

    return is_prime


def sieve_primes_fast(limit: int) -> list:
    """
    Gibt Liste aller Primzahlen bis limit zurück (Numba-beschleunigt).

    Wrapper um sieve_numpy(), der das boolean Array in eine int-Liste umwandelt.
    Bei limit < 2 wird sofort eine leere Liste zurückgegeben.

    @param limit: Obere Grenze (inklusiv)
    @return: Liste aller Primzahlen ≤ limit
    @author: Michael Fuhrmann
    @date: 2026-03-11
    @lastModified: 2026-03-11
    """
    if limit < 2:
        return []
    # Numba-JIT liefert boolean Array — np.where gibt Indizes der True-Einträge
    arr = sieve_numpy(limit)
    return list(np.where(arr)[0])


@njit(cache=True)
def eta_euler_accelerated_jit(s_real: float, s_imag: float, n: int) -> tuple:
    """
    Euler-Knopp-Beschleunigung für η(s) = Σ (-1)^{k-1}/k^s.

    Berechnet den Wert der Dirichlet-Eta-Funktion über die E₁-Transformation:
        S = Σ_{k=0}^{n-1} (-1)^k · (1/2^{k+1}) · Δ^k c_0
    mit c_j = 1/(j+1)^s (unsigned Terme, Vorzeichen extern).

    Gibt (real, imag) des Ergebnisses zurück.
    Numba unterstützt kein Python-complex, daher Real- und Imaginärteil getrennt.

    @param s_real: Realteil von s
    @param s_imag: Imaginärteil von s
    @param n: Anzahl Terme (60 = Maschinengenauigkeit, Fehler ~2^{-60})
    @return: (real, imag) Tuple von η(s)
    @author: Michael Fuhrmann
    @date: 2026-03-11
    @lastModified: 2026-03-11
    """
    # Unsigned Terme c[j] = 1/(j+1)^s als getrennter Real-/Imaginärteil
    # (j+1)^s = exp(s * ln(j+1)) → 1/(j+1)^s = exp(-s * ln(j+1))
    delta_r = np.zeros(n + 1, dtype=np.float64)
    delta_i = np.zeros(n + 1, dtype=np.float64)

    for j in range(n + 1):
        ln_j1 = np.log(j + 1.0)
        # Real-/Imaginärteil von exp(-s * ln(j+1)):
        #   Betrag    = exp(-s_real * ln)
        #   Phase     = -s_imag * ln
        delta_r[j] = np.exp(-s_real * ln_j1) * np.cos(-s_imag * ln_j1)
        delta_i[j] = np.exp(-s_real * ln_j1) * np.sin(-s_imag * ln_j1)

    # Euler-Knopp-Summation: S = Σ (-1)^k · 2^{-(k+1)} · Δ^k c_0
    result_r = 0.0
    result_i = 0.0
    sign = 1.0        # (-1)^k, startet mit +1
    power_half = 0.5  # 1/2^{k+1}, startet mit 0.5

    for k in range(n):
        # Aktuellen Δ^k c_0 akkumulieren
        result_r += sign * power_half * delta_r[0]
        result_i += sign * power_half * delta_i[0]
        sign *= -1.0
        power_half *= 0.5

        # Vorwärtsdifferenz: Δ^{k+1} c_j = Δ^k c_{j+1} - Δ^k c_j
        for j in range(n - k):
            delta_r[j] = delta_r[j + 1] - delta_r[j]
            delta_i[j] = delta_i[j + 1] - delta_i[j]

    return result_r, result_i


def eta_euler_fast(s: complex, n: int = 60) -> complex:
    """
    Numba-beschleunigte Euler-Knopp-Beschleunigung für η(s).

    Wrapper um eta_euler_accelerated_jit() für Python-complex-Kompatibilität.
    10–30× schneller als die reine Python-Implementierung in complex_analysis.py.

    @param s: Komplexe Zahl mit Re(s) > 0
    @param n: Anzahl Terme (Standard: 60 → Maschinengenauigkeit)
    @return: η(s) als complex
    @author: Michael Fuhrmann
    @date: 2026-03-11
    @lastModified: 2026-03-11
    """
    r, i = eta_euler_accelerated_jit(s.real, s.imag, n)
    return complex(r, i)


def warmup():
    """
    Initialisiert den Numba-JIT-Compiler durch einen ersten Dummy-Aufruf.

    Nach diesem Aufruf sind alle @njit-Funktionen vollständig kompiliert
    und werden ohne Verzögerung aufgerufen. Sollte beim Programmstart
    einmalig aufgerufen werden, wenn Performance kritisch ist.

    @author: Michael Fuhrmann
    @date: 2026-03-11
    @lastModified: 2026-03-11
    """
    if NUMBA_AVAILABLE:
        # Kleine Eingaben für schnellen Warmup (~1–3 Sek. beim ersten Start)
        sieve_numpy(100)
        eta_euler_accelerated_jit(2.0, 0.0, 10)
