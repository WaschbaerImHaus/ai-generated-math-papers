"""
@file test_algebra_numbertheory.py
@brief Tests für algebra_numbertheory.py – is_prime, prime_factorization, euler_phi, RSA.
@description
    Überprüft alle zahlentheoretischen Funktionen mit Fokus auf:
    - Korrektheit der Ergebnisse
    - Caching-Verhalten (lru_cache)
    - euler_phi nutzt prime_factorization intern
    - Randfälle und Sonderwerte

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from algebra_numbertheory import is_prime, prime_factorization, euler_phi


class TestIsPrime:
    """Tests für is_prime()."""

    def test_prime_2(self):
        """Test: 2 ist eine Primzahl."""
        assert is_prime(2) is True

    def test_prime_3(self):
        """Test: 3 ist eine Primzahl."""
        assert is_prime(3) is True

    def test_prime_97(self):
        """Test: 97 ist eine Primzahl."""
        assert is_prime(97) is True

    def test_composite_4(self):
        """Test: 4 ist keine Primzahl."""
        assert is_prime(4) is False

    def test_composite_100(self):
        """Test: 100 ist keine Primzahl."""
        assert is_prime(100) is False

    def test_edge_1(self):
        """Test: 1 ist keine Primzahl."""
        assert is_prime(1) is False

    def test_edge_0(self):
        """Test: 0 ist keine Primzahl."""
        assert is_prime(0) is False

    def test_negative(self):
        """Test: Negative Zahl ist keine Primzahl."""
        assert is_prime(-7) is False

    def test_large_prime(self):
        """Test: Große Primzahl 9999991 wird korrekt erkannt."""
        assert is_prime(9999991) is True

    def test_cached_result(self):
        """Test: Wiederholter Aufruf nutzt Cache (schneller)."""
        is_prime.cache_clear()
        # Erste Berechnung
        start1 = time.perf_counter()
        is_prime(9999991)
        t1 = time.perf_counter() - start1

        # Zweite Berechnung (gecacht)
        start2 = time.perf_counter()
        is_prime(9999991)
        t2 = time.perf_counter() - start2

        # Cache-Treffer sollte deutlich schneller sein
        assert t2 < t1 or t2 < 1e-5  # Gecachter Aufruf < 10 Mikrosekunden


class TestPrimeFactorization:
    """Tests für prime_factorization()."""

    def test_factorize_12(self):
        """Test: 12 = 2² × 3."""
        assert prime_factorization(12) == {2: 2, 3: 1}

    def test_factorize_360(self):
        """Test: 360 = 2³ × 3² × 5."""
        assert prime_factorization(360) == {2: 3, 3: 2, 5: 1}

    def test_prime_7(self):
        """Test: Primzahl 7 hat Faktorisierung {7: 1}."""
        assert prime_factorization(7) == {7: 1}

    def test_power_of_2(self):
        """Test: 2^10 = 1024."""
        assert prime_factorization(1024) == {2: 10}

    def test_large_number(self):
        """Test: 10^12 = 2^12 × 5^12."""
        result = prime_factorization(10**12)
        assert 2 in result and 5 in result
        assert result[2] == 12
        assert result[5] == 12

    def test_lru_cache_present(self):
        """Test: prime_factorization hat @lru_cache-Dekorator."""
        assert hasattr(prime_factorization, 'cache_info')

    def test_cached_result(self):
        """Test: Wiederholter Aufruf nutzt lru_cache."""
        prime_factorization.cache_clear()

        prime_factorization(100003)  # Erzwinge ersten Aufruf (Berechnung)
        info_after_first = prime_factorization.cache_info()
        assert info_after_first.misses >= 1

        prime_factorization(100003)  # Zweiter Aufruf → Cache-Treffer
        info_after_second = prime_factorization.cache_info()
        assert info_after_second.hits >= 1

    def test_reconstruction(self):
        """Test: Aus den Faktoren kann n rekonstruiert werden."""
        for n in [12, 360, 1000, 7919]:
            factors = prime_factorization(n)
            reconstructed = 1
            for p, e in factors.items():
                reconstructed *= p ** e
            assert reconstructed == n, f"Rekonstruktion von {n} fehlgeschlagen"


class TestEulerPhi:
    """Tests für euler_phi() – Prüfung der Implementierung und des Cachings."""

    def test_phi_1(self):
        """Test: phi(1) = 1."""
        assert euler_phi(1) == 1

    def test_phi_prime_7(self):
        """Test: phi(7) = 6 (7 ist Primzahl, phi(p) = p-1)."""
        assert euler_phi(7) == 6

    def test_phi_12(self):
        """Test: phi(12) = 4."""
        assert euler_phi(12) == 4

    def test_phi_prime_power(self):
        """Test: phi(8) = 4 (phi(2^3) = 2^2 = 4)."""
        assert euler_phi(8) == 4

    def test_phi_100(self):
        """Test: phi(100) = 40."""
        assert euler_phi(100) == 40

    def test_phi_prime_large(self):
        """Test: phi(p) = p-1 für große Primzahl p=97."""
        assert euler_phi(97) == 96

    def test_phi_uses_prime_factorization(self):
        """Test: euler_phi() ruft prime_factorization() intern auf (gecacht)."""
        # Nach euler_phi()-Aufruf muss prime_factorization-Cache Treffer zeigen
        prime_factorization.cache_clear()
        euler_phi.cache_clear()

        # Für n=12: prime_factorization(12) wird intern aufgerufen
        euler_phi(12)

        # prime_factorization sollte gecacht worden sein (mindestens 1 Miss)
        info = prime_factorization.cache_info()
        assert info.misses >= 1, "prime_factorization wurde nicht intern aufgerufen"

    def test_phi_formula(self):
        """Test: phi(n) = n * prod(1 - 1/p) für n mit bekannten Primfaktoren."""
        # phi(30) = 30 * (1-1/2) * (1-1/3) * (1-1/5) = 30 * 1/2 * 2/3 * 4/5 = 8
        assert euler_phi(30) == 8

    def test_phi_multiplicative(self):
        """Test: phi(m*n) = phi(m)*phi(n) wenn ggT(m,n) = 1."""
        # phi(35) = phi(5) * phi(7) = 4 * 6 = 24
        assert euler_phi(35) == euler_phi(5) * euler_phi(7)

    def test_phi_lru_cache(self):
        """Test: euler_phi hat @lru_cache-Dekorator."""
        assert hasattr(euler_phi, 'cache_info')

    def test_phi_shared_cache_with_factorization(self):
        """Test: euler_phi und prime_factorization teilen den Cache effizient."""
        prime_factorization.cache_clear()
        euler_phi.cache_clear()

        # Zwei verschiedene Aufrufe mit demselben n
        result1 = euler_phi(360)
        result2 = euler_phi(360)

        # Zweiter Aufruf ist gecacht
        assert result1 == result2
        info = euler_phi.cache_info()
        assert info.hits >= 1
