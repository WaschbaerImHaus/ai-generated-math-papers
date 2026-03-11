"""
@file test_math_helpers.py
@brief Tests für das zentrale Hilfsfunktionen-Modul math_helpers.py.
@description
    Umfassende pytest-Tests für alle Funktionen in math_helpers.py.
    Deckt Standardfälle, Edge-Cases und mathematische Korrektheit ab.

    Testbereiche:
    - is_prime(): Primzahltest
    - prime_factorization(): Primfaktorzerlegung
    - primes_up_to(): Sieb des Eratosthenes
    - primes_generator(): Primzahlgenerator
    - first_n_primes(): Erste n Primzahlen
    - gcd(): Größter gemeinsamer Teiler
    - lcm(): Kleinstes gemeinsames Vielfaches
    - extended_gcd(): Erweiterter euklidischer Algorithmus
    - mod_inverse(): Modulares Inverses
    - euler_phi(): Euler-Phi-Funktion
    - legendre_symbol(): Legendre-Symbol
    - jacobi_symbol(): Jacobi-Symbol
    - chinese_remainder_theorem(): Chinesischer Restsatz

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from math_helpers import (
    is_prime, prime_factorization, primes_up_to, primes_generator,
    first_n_primes, gcd, lcm, extended_gcd, mod_inverse,
    euler_phi, legendre_symbol, jacobi_symbol, chinese_remainder_theorem
)


# =============================================================================
# TESTS: is_prime
# =============================================================================

class TestIsPrime:
    """Tests für die Primzahltest-Funktion."""

    def test_kleine_primzahlen(self):
        """Bekannte kleine Primzahlen müssen korrekt erkannt werden."""
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
            assert is_prime(p), f"{p} sollte eine Primzahl sein"

    def test_keine_primzahlen(self):
        """Zusammengesetzte Zahlen dürfen nicht als prim erkannt werden."""
        for n in [0, 1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 25, 100]:
            assert not is_prime(n), f"{n} sollte keine Primzahl sein"

    def test_negativ(self):
        """Negative Zahlen sind keine Primzahlen."""
        for n in [-1, -2, -7, -100]:
            assert not is_prime(n)

    def test_grosse_primzahl(self):
        """Große Primzahl korrekt erkennen."""
        assert is_prime(104729)   # bekannte Primzahl
        assert not is_prime(104730)

    def test_mersenne_primzahl(self):
        """M_7 = 127 ist eine Mersenne-Primzahl."""
        assert is_prime(127)

    def test_twin_primes(self):
        """Zwillingsprimzahlen-Paare."""
        twins = [(11, 13), (17, 19), (29, 31), (41, 43), (59, 61)]
        for p, q in twins:
            assert is_prime(p) and is_prime(q)


# =============================================================================
# TESTS: prime_factorization
# =============================================================================

class TestPrimeFactorization:
    """Tests für die Primfaktorzerlegung."""

    def test_einfache_primzahl(self):
        """Primzahl hat nur sich selbst als Faktor."""
        assert prime_factorization(7) == {7: 1}
        assert prime_factorization(13) == {13: 1}

    def test_primzahlpotenz(self):
        """Primzahlpotenzen korrekt zerlegen."""
        assert prime_factorization(8) == {2: 3}   # 2^3
        assert prime_factorization(27) == {3: 3}  # 3^3
        assert prime_factorization(32) == {2: 5}  # 2^5

    def test_zusammengesetzt(self):
        """Zusammengesetzte Zahlen korrekt zerlegen."""
        assert prime_factorization(12) == {2: 2, 3: 1}   # 4·3
        assert prime_factorization(360) == {2: 3, 3: 2, 5: 1}  # 8·9·5

    def test_produkt_verifizieren(self):
        """Rekonstruiertes Produkt muss n ergeben."""
        for n in [100, 360, 1000, 9999, 10000]:
            factors = prime_factorization(n)
            product = 1
            for p, e in factors.items():
                product *= p ** e
            assert product == n, f"Zerlegung von {n} inkonsistent"

    def test_kleine_zahlen(self):
        """Sonderfälle kleine Zahlen."""
        assert prime_factorization(2) == {2: 1}
        assert prime_factorization(4) == {2: 2}


# =============================================================================
# TESTS: primes_up_to
# =============================================================================

class TestPrimesUpTo:
    """Tests für das Sieb des Eratosthenes."""

    def test_bis_30(self):
        """Primzahlen bis 30 korrekt bestimmen."""
        expected = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        assert primes_up_to(30) == expected

    def test_bis_2(self):
        """Nur 2 bis 2."""
        assert primes_up_to(2) == [2]

    def test_bis_1(self):
        """Leere Liste bis 1."""
        assert primes_up_to(1) == []

    def test_anzahl_bis_100(self):
        """π(100) = 25 Primzahlen bis 100."""
        assert len(primes_up_to(100)) == 25

    def test_anzahl_bis_1000(self):
        """π(1000) = 168 Primzahlen bis 1000."""
        assert len(primes_up_to(1000)) == 168

    def test_alle_sind_prim(self):
        """Alle zurückgegebenen Zahlen müssen prim sein."""
        for p in primes_up_to(200):
            assert is_prime(p), f"{p} in Sieb ist keine Primzahl"


# =============================================================================
# TESTS: first_n_primes
# =============================================================================

class TestFirstNPrimes:
    """Tests für erste-n-Primzahlen-Funktion."""

    def test_erste_10(self):
        """Die ersten 10 Primzahlen."""
        expected = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        assert first_n_primes(10) == expected

    def test_erste_1(self):
        """Erste Primzahl ist 2."""
        assert first_n_primes(1) == [2]


# =============================================================================
# TESTS: gcd
# =============================================================================

class TestGcd:
    """Tests für den größten gemeinsamen Teiler."""

    def test_bekannte_werte(self):
        """ggT bekannter Zahlenpaare."""
        assert gcd(12, 8) == 4
        assert gcd(17, 5) == 1   # teilerfremd
        assert gcd(100, 75) == 25

    def test_symmetrie(self):
        """ggT ist symmetrisch: ggT(a,b) = ggT(b,a)."""
        assert gcd(48, 18) == gcd(18, 48)

    def test_mit_null(self):
        """ggT(a, 0) = a."""
        assert gcd(5, 0) == 5
        assert gcd(0, 7) == 7
        assert gcd(0, 0) == 0

    def test_negative_zahlen(self):
        """ggT ist immer nichtnegativ."""
        assert gcd(-12, 8) == 4
        assert gcd(-15, -10) == 5

    def test_identitaet(self):
        """ggT(a, a) = a."""
        for a in [1, 5, 12, 100]:
            assert gcd(a, a) == a


# =============================================================================
# TESTS: lcm
# =============================================================================

class TestLcm:
    """Tests für das kleinste gemeinsame Vielfache."""

    def test_bekannte_werte(self):
        """kgV bekannter Paare."""
        assert lcm(4, 6) == 12
        assert lcm(3, 7) == 21
        assert lcm(12, 8) == 24

    def test_teilerfremd(self):
        """kgV(a, b) = a·b wenn ggT(a,b)=1."""
        assert lcm(5, 7) == 35
        assert lcm(11, 13) == 143

    def test_mit_null(self):
        """kgV(0, b) = 0."""
        assert lcm(0, 5) == 0
        assert lcm(3, 0) == 0


# =============================================================================
# TESTS: extended_gcd
# =============================================================================

class TestExtendedGcd:
    """Tests für den erweiterten euklidischen Algorithmus."""

    def test_bezout_identitaet(self):
        """a·x + b·y = ggT(a, b)."""
        for a, b in [(12, 8), (17, 5), (100, 75), (48, 18)]:
            g, x, y = extended_gcd(a, b)
            assert a * x + b * y == g, f"Bézout-Identität für ({a},{b}) verletzt"
            assert g == gcd(a, b)

    def test_teilerfremd(self):
        """Für teilerfremde Zahlen ist ggT=1."""
        g, x, y = extended_gcd(17, 5)
        assert g == 1
        assert 17 * x + 5 * y == 1


# =============================================================================
# TESTS: mod_inverse
# =============================================================================

class TestModInverse:
    """Tests für das modulare Inverse."""

    def test_einfache_inverse(self):
        """Bekannte modulare Inverse."""
        # 3^{-1} mod 7 = 5, weil 3·5 = 15 ≡ 1 (mod 7)
        assert mod_inverse(3, 7) == 5
        # 2^{-1} mod 5 = 3, weil 2·3 = 6 ≡ 1 (mod 5)
        assert mod_inverse(2, 5) == 3

    def test_korrektheit(self):
        """a · a^{-1} ≡ 1 (mod m)."""
        for a, m in [(3, 7), (5, 11), (7, 13), (17, 100)]:
            if gcd(a, m) == 1:
                inv = mod_inverse(a, m)
                assert (a * inv) % m == 1

    def test_kein_inverses(self):
        """Fehler wenn ggT(a, m) ≠ 1."""
        with pytest.raises(ValueError):
            mod_inverse(6, 4)  # ggT(6, 4) = 2 ≠ 1


# =============================================================================
# TESTS: euler_phi
# =============================================================================

class TestEulerPhi:
    """Tests für Eulers Phi-Funktion."""

    def test_bekannte_werte(self):
        """Bekannte φ-Werte."""
        assert euler_phi(1) == 1
        assert euler_phi(2) == 1
        assert euler_phi(6) == 2
        assert euler_phi(12) == 4

    def test_primzahl(self):
        """φ(p) = p-1 für Primzahlen."""
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
            assert euler_phi(p) == p - 1, f"φ({p}) sollte {p-1} sein"

    def test_primzahlpotenz(self):
        """φ(p^k) = p^{k-1}·(p-1)."""
        assert euler_phi(4) == 2    # φ(2^2) = 2
        assert euler_phi(8) == 4    # φ(2^3) = 4
        assert euler_phi(9) == 6    # φ(3^2) = 6
        assert euler_phi(27) == 18  # φ(3^3) = 18

    def test_multiplikativ(self):
        """φ(mn) = φ(m)·φ(n) wenn ggT(m,n)=1."""
        # ggT(5, 7) = 1 → φ(35) = φ(5)·φ(7) = 4·6 = 24
        assert euler_phi(35) == euler_phi(5) * euler_phi(7)
        # ggT(4, 9) = 1 → φ(36) = φ(4)·φ(9) = 2·6 = 12
        assert euler_phi(36) == euler_phi(4) * euler_phi(9)

    def test_fehler_bei_null(self):
        """φ(0) soll ValueError werfen."""
        with pytest.raises(ValueError):
            euler_phi(0)


# =============================================================================
# TESTS: legendre_symbol
# =============================================================================

class TestLegendreSymbol:
    """Tests für das Legendre-Symbol."""

    def test_quadratische_reste(self):
        """Quadratische Reste mod 7: 1,2,4 (weil 1²=1, 2²=4, 3²=2 mod 7)."""
        # 1 ≡ 1² mod 7, also (1/7) = 1
        assert legendre_symbol(1, 7) == 1
        # 2 ≡ 3² mod 7, also (2/7) = 1
        assert legendre_symbol(2, 7) == 1
        # 4 ≡ 2² mod 7, also (4/7) = 1
        assert legendre_symbol(4, 7) == 1

    def test_nichtrest(self):
        """Quadratische Nichtreste mod 7: 3, 5, 6."""
        assert legendre_symbol(3, 7) == -1
        assert legendre_symbol(5, 7) == -1
        assert legendre_symbol(6, 7) == -1

    def test_teiler(self):
        """(p/p) = 0."""
        assert legendre_symbol(7, 7) == 0
        assert legendre_symbol(0, 5) == 0

    def test_euler_kriterium(self):
        """(a/p) ≡ a^{(p-1)/2} (mod p)."""
        p = 11
        for a in range(1, p):
            expected = pow(a, (p - 1) // 2, p)
            expected = 1 if expected == 1 else (-1 if expected == p - 1 else 0)
            assert legendre_symbol(a, p) == expected


# =============================================================================
# TESTS: jacobi_symbol
# =============================================================================

class TestJacobiSymbol:
    """Tests für das Jacobi-Symbol."""

    def test_stimmt_mit_legendre_fuer_primzahl_ueberein(self):
        """Für Primzahl n stimmt Jacobi mit Legendre überein."""
        for p in [3, 5, 7, 11, 13]:
            for a in range(0, p):
                assert jacobi_symbol(a, p) == legendre_symbol(a, p), \
                    f"Jacobi ≠ Legendre für a={a}, p={p}"

    def test_eins(self):
        """(1/n) = 1 für alle ungeraden n."""
        for n in [1, 3, 5, 7, 9, 15, 21]:
            assert jacobi_symbol(1, n) == 1

    def test_null(self):
        """(0/n) = 0 für n > 1."""
        for n in [3, 5, 7, 9]:
            assert jacobi_symbol(0, n) == 0

    def test_multiplikativ(self):
        """(a·b/n) = (a/n)·(b/n)."""
        n = 15
        for a in range(1, 8):
            for b in range(1, 8):
                assert jacobi_symbol(a * b, n) == jacobi_symbol(a, n) * jacobi_symbol(b, n), \
                    f"Multiplikativität verletzt für a={a}, b={b}, n={n}"


# =============================================================================
# TESTS: chinese_remainder_theorem
# =============================================================================

class TestChineseRemainderTheorem:
    """Tests für den Chinesischen Restsatz."""

    def test_einfaches_system(self):
        """x ≡ 2 (mod 3), x ≡ 3 (mod 5) → x = 8."""
        x, M = chinese_remainder_theorem([2, 3], [3, 5])
        assert x == 8
        assert M == 15
        # Nachprüfen
        assert x % 3 == 2
        assert x % 5 == 3

    def test_drei_gleichungen(self):
        """x ≡ 1 (mod 2), x ≡ 2 (mod 3), x ≡ 3 (mod 5) → x = 23."""
        x, M = chinese_remainder_theorem([1, 2, 3], [2, 3, 5])
        assert M == 30
        assert x % 2 == 1
        assert x % 3 == 2
        assert x % 5 == 3

    def test_loesung_konsistent(self):
        """Lösung muss alle Kongruenzen erfüllen."""
        remainders = [3, 5, 7]
        moduli = [11, 13, 17]
        x, M = chinese_remainder_theorem(remainders, moduli)
        for r, m in zip(remainders, moduli):
            assert x % m == r, f"CRT-Lösung erfüllt {x} ≡ {r} (mod {m}) nicht"

    def test_gesamtmodulus(self):
        """Gesamtmodulus ist das Produkt der Moduln."""
        _, M = chinese_remainder_theorem([1, 2], [3, 5])
        assert M == 15
