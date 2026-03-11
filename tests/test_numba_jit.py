"""
Tests für das numba_jit-Modul.

Prüft Korrektheit der JIT-kompilierten und Python-Fallback-Implementierungen
von Sieb des Eratosthenes und Euler-Knopp-Beschleunigung für die Eta-Funktion.

@author: Michael Fuhrmann
@date: 2026-03-11
@lastModified: 2026-03-11
"""

import math
import pytest
import numpy as np


# Modul-Import — funktioniert mit und ohne Numba
from src.numba_jit import (
    sieve_primes_fast,
    eta_euler_fast,
    NUMBA_AVAILABLE,
    sieve_numpy,
    eta_euler_accelerated_jit,
)


# ===========================================================================
# Tests: sieve_primes_fast()
# ===========================================================================

class TestSievePrimesFast:
    """Tests für die Numba-beschleunigte Sieb-Implementierung."""

    def test_sieve_correctness(self):
        """
        Vergleicht sieve_primes_fast(100) mit der bekannten Liste aller
        Primzahlen bis 100 (25 Stück).
        """
        # Bekannte Primzahlen bis 100
        expected = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
            31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
            73, 79, 83, 89, 97
        ]
        result = sieve_primes_fast(100)
        assert result == expected, (
            f"Primzahlen bis 100 stimmen nicht: {result}"
        )

    def test_sieve_count_to_100(self):
        """Es gibt genau 25 Primzahlen ≤ 100."""
        result = sieve_primes_fast(100)
        assert len(result) == 25

    def test_sieve_small_limits(self):
        """Randfall: Grenzen 0, 1 und 2."""
        assert sieve_primes_fast(0) == []
        assert sieve_primes_fast(1) == []
        assert sieve_primes_fast(2) == [2]

    def test_sieve_limit_is_prime(self):
        """Wenn limit selbst eine Primzahl ist, muss sie enthalten sein."""
        result = sieve_primes_fast(97)
        assert 97 in result
        result2 = sieve_primes_fast(7)
        assert result2 == [2, 3, 5, 7]

    def test_sieve_limit_is_composite(self):
        """Wenn limit eine zusammengesetzte Zahl ist, darf sie nicht enthalten sein."""
        result = sieve_primes_fast(100)
        assert 100 not in result
        assert 99 not in result
        assert 98 not in result

    def test_sieve_large(self):
        """
        sieve_primes_fast(10000): letzte Primzahl muss 9973 sein.
        π(10000) = 1229 Primzahlen.
        """
        result = sieve_primes_fast(10000)
        assert result[-1] == 9973, (
            f"Letzte Primzahl bis 10000 soll 9973 sein, erhalten: {result[-1]}"
        )
        assert len(result) == 1229, (
            f"π(10000) = 1229, erhalten: {len(result)}"
        )

    def test_sieve_first_elements(self):
        """Die ersten 10 Primzahlen müssen korrekt sein."""
        result = sieve_primes_fast(30)
        assert result[:10] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    def test_sieve_numpy_returns_bool_array(self):
        """sieve_numpy() gibt ein numpy boolean Array zurück."""
        arr = sieve_numpy(20)
        assert arr.dtype == np.bool_
        assert arr.shape == (21,)
        # Index 2, 3, 5, 7, 11, 13, 17, 19 müssen True sein
        primes_20 = [2, 3, 5, 7, 11, 13, 17, 19]
        for p in primes_20:
            assert arr[p] == True, f"arr[{p}] soll True sein"
        # Index 0, 1, 4, 6, 8, 9 ... müssen False sein
        composites = [0, 1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20]
        for c in composites:
            assert arr[c] == False, f"arr[{c}] soll False sein"

    def test_sieve_type_is_list(self):
        """sieve_primes_fast() gibt eine Python-Liste zurück (nicht numpy-Array)."""
        result = sieve_primes_fast(50)
        assert isinstance(result, list)


# ===========================================================================
# Tests: eta_euler_fast()
# ===========================================================================

class TestEtaEulerFast:
    """Tests für die Numba-beschleunigte Euler-Knopp-Eta-Berechnung."""

    def test_eta_s2(self):
        """
        η(2) = π²/12 ≈ 0.8224670334241132.

        Aus der Identität: ζ(2) = π²/6, η(2) = (1 - 2^{1-2})·ζ(2) = π²/12.
        """
        expected = math.pi ** 2 / 12
        result = eta_euler_fast(2 + 0j)
        assert abs(result.real - expected) < 1e-10, (
            f"η(2) soll ≈ {expected:.10f} sein, erhalten: {result.real:.10f}"
        )
        assert abs(result.imag) < 1e-10, (
            f"η(2) soll reell sein, Imaginärteil: {result.imag}"
        )

    def test_eta_s1_ln2(self):
        """
        η(1) = ln(2) ≈ 0.6931471805599453 (alternierende harmonische Reihe).
        """
        expected = math.log(2)
        result = eta_euler_fast(1 + 0j)
        assert abs(result.real - expected) < 1e-10, (
            f"η(1) soll ln(2) ≈ {expected:.10f} sein, erhalten: {result.real:.10f}"
        )

    def test_eta_s3(self):
        """
        η(3) = 3/4 · ζ(3) ≈ 3/4 · 1.202056903... ≈ 0.901542677...
        """
        # Apéry-Konstante ζ(3) ≈ 1.2020569031595942
        zeta3 = 1.2020569031595942
        expected = 0.75 * zeta3  # η(3) = (1-2^{1-3}) · ζ(3) = 3/4 · ζ(3)
        result = eta_euler_fast(3 + 0j)
        assert abs(result.real - expected) < 1e-9, (
            f"η(3) soll ≈ {expected:.9f} sein, erhalten: {result.real:.9f}"
        )

    def test_eta_critical_line(self):
        """
        η(0.5 + 14.135j) — nahe der ersten nichttrivialen Nullstelle von ζ.

        Die erste Nullstelle von ζ liegt bei s ≈ 0.5 + 14.1347i.
        η(s) = (1 - 2^{1-s}) · ζ(s), also ist η(s) dort ≈ (1-2^{0.5-14.135i}) · 0.
        Das Ergebnis soll endlich und nicht NaN/Inf sein.
        """
        s = 0.5 + 14.135j
        result = eta_euler_fast(s)
        # Ergebnis muss endliche komplexe Zahl sein
        assert math.isfinite(result.real), f"Realteil ist nicht endlich: {result.real}"
        assert math.isfinite(result.imag), f"Imaginärteil ist nicht endlich: {result.imag}"

    def test_eta_symmetry_real_axis(self):
        """η(s) für reelle s > 0 gibt reelles Ergebnis (Imaginärteil ≈ 0)."""
        for s_val in [0.5, 1.5, 2.0, 4.0]:
            result = eta_euler_fast(complex(s_val))
            assert abs(result.imag) < 1e-12, (
                f"η({s_val}) soll reell sein, Imaginärteil: {result.imag:.2e}"
            )

    def test_eta_euler_jit_tuple_output(self):
        """eta_euler_accelerated_jit() gibt ein Tupel (float, float) zurück."""
        r, i = eta_euler_accelerated_jit(2.0, 0.0, 60)
        expected = math.pi ** 2 / 12
        assert abs(r - expected) < 1e-10
        assert abs(i) < 1e-10

    def test_eta_term_count_robustness(self):
        """η(2) soll auch mit unterschiedlichen n-Werten korrekt sein."""
        expected = math.pi ** 2 / 12
        for n_terms in [30, 40, 60]:
            result = eta_euler_fast(2 + 0j, n=n_terms)
            assert abs(result.real - expected) < 1e-8, (
                f"η(2) mit n={n_terms} zu ungenau: {result.real}"
            )


# ===========================================================================
# Tests: Modul-Metadaten
# ===========================================================================

class TestModuleMetadata:
    """Tests für Modulstatus und Konfiguration."""

    def test_numba_available(self):
        """
        NUMBA_AVAILABLE soll True sein, da Numba 0.64.0 installiert ist.
        Falls Numba nicht installiert wäre, wäre der Wert False — kein Fehler.
        """
        # Numba ist laut Aufgabenstellung installiert (0.64.0)
        assert NUMBA_AVAILABLE is True, (
            "Numba 0.64.0 soll installiert sein — NUMBA_AVAILABLE ist False. "
            "Bitte 'pip install numba' ausführen."
        )

    def test_numba_available_is_bool(self):
        """NUMBA_AVAILABLE muss ein bool sein."""
        assert isinstance(NUMBA_AVAILABLE, bool)

    def test_sieve_fast_callable(self):
        """sieve_primes_fast muss aufrufbar sein."""
        assert callable(sieve_primes_fast)

    def test_eta_fast_callable(self):
        """eta_euler_fast muss aufrufbar sein."""
        assert callable(eta_euler_fast)
