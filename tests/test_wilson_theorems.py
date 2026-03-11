"""
Testsuite für Wilson-Satz-Erweiterungen (Papers 9–12)
======================================================
Autor: Michael Fuhrmann
Erstellt: 2026-03-11
Letzte Änderung: 2026-03-11

Testet:
- Klassischer Wilson-Satz (Primzahlen bis 50)
- Wilson für Primzahlpotenzen (p=3,5,7 mit k=1,2,3)
- Wilson-Quotient und Wilson-Primzahlen (5, 13, 563)
- Allgemeiner Wilson für Gruppen Z/nZ (n=2,...,20)
- Harmonische Summe modulo p (Wolstenholme)
- Quadratreste: -1 als Quadrat modulo p
- Gausssche Verallgemeinerung
"""

import math
import pytest
from typing import List


# ---------------------------------------------------------------------------
# Hilfsfunktionen (ohne externe Bibliotheken, nur Python-Standardbibliothek)
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """
    Primzahltest via Probedivision.

    :param n: zu testende Zahl
    :returns: True wenn n eine Primzahl ist
    :rtype: bool
    :timestamp: 2026-03-11
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


def mod_inverse(a: int, m: int) -> int:
    """
    Modulares Inverses via erweitertem euklidischen Algorithmus.

    :param a: Zahl, deren Inverses gesucht wird
    :param m: Modulus
    :returns: a^{-1} mod m
    :rtype: int
    :timestamp: 2026-03-11
    """
    g, x, _ = extended_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"Kein Inverses: gcd({a}, {m}) = {g} != 1")
    return x % m


def extended_gcd(a: int, b: int):
    """
    Erweiterter euklidischer Algorithmus.

    :param a: erste Zahl
    :param b: zweite Zahl
    :returns: (gcd, x, y) mit gcd = a*x + b*y
    :rtype: tuple
    :timestamp: 2026-03-11
    """
    if a == 0:
        return b, 0, 1
    g, x, y = extended_gcd(b % a, a)
    return g, y - (b // a) * x, x


def factorial_mod(n: int, m: int) -> int:
    """
    Berechnet n! mod m.

    :param n: Zahl, deren Fakultät berechnet wird
    :param m: Modulus
    :returns: n! mod m
    :rtype: int
    :timestamp: 2026-03-11
    """
    result = 1
    for i in range(2, n + 1):
        result = (result * i) % m
    return result


def product_of_units(n: int) -> int:
    """
    Berechnet das Produkt aller Einheiten modulo n.
    Einheiten sind a mit 1 <= a <= n und gcd(a,n) = 1.

    :param n: Modulus
    :returns: Produkt aller Einheiten mod n
    :rtype: int
    :timestamp: 2026-03-11
    """
    result = 1
    for a in range(1, n):
        if math.gcd(a, n) == 1:
            result = (result * a) % n
    return result


def wilson_quotient(p: int) -> int:
    """
    Berechnet den Wilson-Quotienten W_p = ((p-1)! + 1) / p.
    Nur für Primzahlen p definiert.

    :param p: Primzahl
    :returns: W_p = ((p-1)! + 1) // p
    :rtype: int
    :timestamp: 2026-03-11
    """
    fact = math.factorial(p - 1)
    return (fact + 1) // p


def is_wilson_prime(p: int) -> bool:
    """
    Prüft, ob p eine Wilson-Primzahl ist (p^2 | (p-1)! + 1).

    :param p: Primzahl
    :returns: True wenn p eine Wilson-Primzahl ist
    :rtype: bool
    :timestamp: 2026-03-11
    """
    if not is_prime(p):
        return False
    fact = math.factorial(p - 1)
    return (fact + 1) % (p * p) == 0


def harmonic_sum_mod_p(p: int) -> int:
    """
    Berechnet die harmonische Summe 1 + 1/2 + ... + 1/(p-1) modulo p.
    Jeder Bruch 1/j ist das multiplikative Inverse von j modulo p.

    :param p: Primzahl >= 3
    :returns: harmonische Summe mod p
    :rtype: int
    :timestamp: 2026-03-11
    """
    total = 0
    for j in range(1, p):
        total = (total + mod_inverse(j, p)) % p
    return total


def is_quadratic_residue_minus1(p: int) -> bool:
    """
    Prüft, ob -1 ein quadratischer Rest modulo p ist.
    Nach Euler: -1 ist QR mod p genau dann wenn p ≡ 1 (mod 4).

    :param p: ungerade Primzahl
    :returns: True wenn -1 ein QR mod p ist
    :rtype: bool
    :timestamp: 2026-03-11
    """
    # Direkte Prüfung via Euler-Kriterium: (-1)^((p-1)/2) mod p
    exp = (p - 1) // 2
    return pow(-1, exp) == 1  # also: exp gerade, also p ≡ 1 mod 4


def gauss_wilson_expected(n: int) -> int:
    """
    Gibt den erwarteten Wert des Produkts der Einheiten modulo n zurück
    gemäß der Gaußschen Verallgemeinerung:
    -1 (d.h. n-1) wenn n = 1,2,4,p^k oder 2p^k; sonst 1.

    :param n: Modulus >= 1
    :returns: erwarteter Wert (0 oder n-1)
    :rtype: int
    :timestamp: 2026-03-11
    """
    if n == 1:
        return 0  # -1 mod 1 = 0
    if n == 2:
        return 1  # -1 mod 2 = 1
    if n == 4:
        return 3  # -1 mod 4 = 3
    # Prüfe n = p^k (ungerade Primzahl p)
    for p in range(3, n + 1, 2):
        if not is_prime(p):
            continue
        pk = p
        while pk <= n:
            if pk == n:
                return n - 1
            pk *= p
    # Prüfe n = 2*p^k
    for p in range(3, n + 1, 2):
        if not is_prime(p):
            continue
        pk = p
        while 2 * pk <= n:
            if 2 * pk == n:
                return n - 1
            pk *= p
    return 1


# ---------------------------------------------------------------------------
# Testklassen
# ---------------------------------------------------------------------------

class TestClassicalWilson:
    """
    Tests für den klassischen Satz von Wilson: (p-1)! ≡ -1 (mod p)
    genau dann wenn p prim ist.
    """

    def test_wilson_prime_2(self):
        """Wilson-Satz für die kleinste Primzahl p=2."""
        assert factorial_mod(1, 2) % 2 == 1  # 1! = 1 ≡ -1 ≡ 1 (mod 2)

    def test_wilson_prime_3(self):
        """Wilson-Satz für p=3: 2! = 2 ≡ -1 (mod 3)."""
        assert factorial_mod(2, 3) == 2  # 2! = 2 ≡ -1 mod 3

    def test_wilson_prime_5(self):
        """Wilson-Satz für p=5: 4! = 24 ≡ -1 (mod 5)."""
        assert factorial_mod(4, 5) == 4  # 24 mod 5 = 4 = -1 mod 5

    def test_wilson_prime_7(self):
        """Wilson-Satz für p=7: 6! = 720 ≡ -1 (mod 7)."""
        assert factorial_mod(6, 7) == 6

    def test_wilson_prime_11(self):
        """Wilson-Satz für p=11."""
        assert factorial_mod(10, 11) == 10

    def test_wilson_prime_13(self):
        """Wilson-Satz für p=13."""
        assert factorial_mod(12, 13) == 12

    def test_wilson_prime_17(self):
        """Wilson-Satz für p=17."""
        assert factorial_mod(16, 17) == 16

    def test_wilson_prime_19(self):
        """Wilson-Satz für p=19."""
        assert factorial_mod(18, 19) == 18

    def test_wilson_prime_23(self):
        """Wilson-Satz für p=23."""
        assert factorial_mod(22, 23) == 22

    def test_wilson_prime_29(self):
        """Wilson-Satz für p=29."""
        assert factorial_mod(28, 29) == 28

    def test_wilson_prime_31(self):
        """Wilson-Satz für p=31."""
        assert factorial_mod(30, 31) == 30

    def test_wilson_prime_37(self):
        """Wilson-Satz für p=37."""
        assert factorial_mod(36, 37) == 36

    def test_wilson_prime_41(self):
        """Wilson-Satz für p=41."""
        assert factorial_mod(40, 41) == 40

    def test_wilson_prime_43(self):
        """Wilson-Satz für p=43."""
        assert factorial_mod(42, 43) == 42

    def test_wilson_prime_47(self):
        """Wilson-Satz für p=47."""
        assert factorial_mod(46, 47) == 46

    def test_wilson_all_primes_up_to_50(self):
        """Alle Primzahlen bis 50 erfüllen den Wilson-Satz."""
        primes_up_to_50 = [p for p in range(2, 51) if is_prime(p)]
        for p in primes_up_to_50:
            assert factorial_mod(p - 1, p) == p - 1, \
                f"Wilson-Satz verletzt für p={p}: {factorial_mod(p-1,p)} != {p-1}"

    def test_wilson_composite_4(self):
        """Zusammengesetzte Zahl 4 erfüllt den Wilson-Satz nicht."""
        # (4-1)! = 6 ≡ 2 (mod 4) ≠ -1 ≡ 3 (mod 4)
        assert factorial_mod(3, 4) == 2

    def test_wilson_composite_6(self):
        """Zusammengesetzte Zahl 6 erfüllt den Wilson-Satz nicht."""
        assert factorial_mod(5, 6) == 0  # 5! = 120 ≡ 0 (mod 6)

    def test_wilson_composite_9(self):
        """Zusammengesetzte Zahl 9 erfüllt den Wilson-Satz nicht."""
        assert factorial_mod(8, 9) == 0  # 8! = 40320 ≡ 0 (mod 9)

    def test_wilson_composite_15(self):
        """Zusammengesetzte Zahl 15 erfüllt den Wilson-Satz nicht."""
        assert factorial_mod(14, 15) == 0


class TestWilsonPrimePowers:
    """
    Tests für Wilson-Satz bei Primzahlpotenzen.
    Für ungerades p und k>=1: Produkt der Einheiten mod p^k ≡ -1.
    """

    def test_wilson_p3_k1(self):
        """Wilson für p=3, k=1: Produkt der Einheiten mod 3."""
        # Einheiten mod 3: {1, 2}, Produkt = 2 ≡ -1 (mod 3)
        assert product_of_units(3) == 2

    def test_wilson_p3_k2(self):
        """Wilson für p=3, k=2: Produkt der Einheiten mod 9."""
        # Einheiten mod 9: {1,2,4,5,7,8}, Produkt ≡ -1 ≡ 8 (mod 9)
        assert product_of_units(9) == 8

    def test_wilson_p3_k3(self):
        """Wilson für p=3, k=3: Produkt der Einheiten mod 27."""
        # Produkt aller Einheiten mod 27 ≡ -1 ≡ 26 (mod 27)
        assert product_of_units(27) == 26

    def test_wilson_p5_k1(self):
        """Wilson für p=5, k=1: Produkt der Einheiten mod 5."""
        assert product_of_units(5) == 4  # ≡ -1 (mod 5)

    def test_wilson_p5_k2(self):
        """Wilson für p=5, k=2: Produkt der Einheiten mod 25."""
        assert product_of_units(25) == 24  # ≡ -1 (mod 25)

    def test_wilson_p5_k3(self):
        """Wilson für p=5, k=3: Produkt der Einheiten mod 125."""
        assert product_of_units(125) == 124  # ≡ -1 (mod 125)

    def test_wilson_p7_k1(self):
        """Wilson für p=7, k=1: Produkt der Einheiten mod 7."""
        assert product_of_units(7) == 6  # ≡ -1 (mod 7)

    def test_wilson_p7_k2(self):
        """Wilson für p=7, k=2: Produkt der Einheiten mod 49."""
        assert product_of_units(49) == 48  # ≡ -1 (mod 49)

    def test_wilson_p2_k1(self):
        """Wilson für p=2, k=1: Produkt mod 2 ist 1 (nicht -1 in üblichem Sinn)."""
        assert product_of_units(2) == 1  # Einzige Einheit: 1

    def test_wilson_p2_k2(self):
        """Wilson für p=2, k=2: Produkt der Einheiten mod 4 ist 3 ≡ -1 (mod 4)."""
        assert product_of_units(4) == 3  # {1,3}, Produkt = 3 ≡ -1 (mod 4)

    def test_wilson_p2_k3(self):
        """Wilson für p=2, k=3: Produkt der Einheiten mod 8 ist 1 (nicht -1!)."""
        # mod 8: Einheiten = {1,3,5,7}, Produkt = 105 ≡ 1 (mod 8)
        assert product_of_units(8) == 1

    def test_wilson_p2_k4(self):
        """Wilson für p=2, k=4: Produkt der Einheiten mod 16 ist 1."""
        assert product_of_units(16) == 1


class TestWilsonQuotient:
    """
    Tests für den Wilson-Quotienten W_p = ((p-1)! + 1) / p
    und Wilson-Primzahlen.
    """

    def test_wilson_quotient_p2(self):
        """Wilson-Quotient für p=2: W_2 = 1."""
        assert wilson_quotient(2) == 1

    def test_wilson_quotient_p3(self):
        """Wilson-Quotient für p=3: W_3 = 1."""
        assert wilson_quotient(3) == 1

    def test_wilson_quotient_p5(self):
        """Wilson-Quotient für p=5: W_5 = 5 (Wilson-Primzahl!)."""
        assert wilson_quotient(5) == 5

    def test_wilson_quotient_p7(self):
        """Wilson-Quotient für p=7: W_7 = 103."""
        assert wilson_quotient(7) == 103

    def test_wilson_quotient_p11(self):
        """Wilson-Quotient für p=11: W_11 = 329891."""
        assert wilson_quotient(11) == 329891

    def test_wilson_prime_5(self):
        """p=5 ist eine Wilson-Primzahl: 5^2 | 4! + 1 = 25."""
        assert is_wilson_prime(5) is True
        assert (math.factorial(4) + 1) % 25 == 0

    def test_wilson_prime_13(self):
        """p=13 ist eine Wilson-Primzahl: 13^2 | 12! + 1."""
        assert is_wilson_prime(13) is True
        assert (math.factorial(12) + 1) % (13 * 13) == 0

    def test_not_wilson_prime_3(self):
        """p=3 ist KEINE Wilson-Primzahl."""
        assert is_wilson_prime(3) is False

    def test_not_wilson_prime_7(self):
        """p=7 ist KEINE Wilson-Primzahl: W_7 = 103, 103 mod 7 = 5 ≠ 0."""
        assert is_wilson_prime(7) is False
        assert wilson_quotient(7) % 7 == 5

    def test_not_wilson_prime_11(self):
        """p=11 ist KEINE Wilson-Primzahl."""
        assert is_wilson_prime(11) is False

    def test_not_wilson_prime_17(self):
        """p=17 ist KEINE Wilson-Primzahl."""
        assert is_wilson_prime(17) is False

    def test_wilson_quotient_divisibility_by_p5(self):
        """W_5 ist durch 5 teilbar (Wilson-Primzahl-Kriterium)."""
        assert wilson_quotient(5) % 5 == 0

    def test_wilson_quotient_divisibility_by_p13(self):
        """W_13 ist durch 13 teilbar (Wilson-Primzahl-Kriterium)."""
        assert wilson_quotient(13) % 13 == 0


class TestGeneralWilsonAbelian:
    """
    Tests für den allgemeinen Wilson-Satz für Gruppen (Z/nZ)*.
    Produkt der Einheiten = -1 mod n genau für n = 1,2,4,p^k,2p^k.
    """

    def test_wilson_n2(self):
        """Produkt der Einheiten mod 2 ist 1 ≡ -1 (mod 2)."""
        assert product_of_units(2) == 1  # 1 ≡ -1 mod 2

    def test_wilson_n3(self):
        """Produkt der Einheiten mod 3 ist 2 ≡ -1 (mod 3)."""
        assert product_of_units(3) == 2

    def test_wilson_n4(self):
        """Produkt der Einheiten mod 4 ist 3 ≡ -1 (mod 4)."""
        assert product_of_units(4) == 3

    def test_wilson_n5(self):
        """Produkt der Einheiten mod 5 ist 4 ≡ -1 (mod 5)."""
        assert product_of_units(5) == 4

    def test_wilson_n6(self):
        """Produkt der Einheiten mod 6 = 2*3 ist 5 ≡ -1 (mod 6)."""
        # 6 = 2*3, also n = 2*p^1 mit p=3. Produkt = -1.
        assert product_of_units(6) == 5

    def test_wilson_n7(self):
        """Produkt der Einheiten mod 7 ist 6 ≡ -1 (mod 7)."""
        assert product_of_units(7) == 6

    def test_wilson_n8(self):
        """Produkt der Einheiten mod 8 ist 1 (nicht -1), da 8 = 2^3."""
        # 8 = 2^3: Produkt ist 1, nicht -1
        assert product_of_units(8) == 1

    def test_wilson_n9(self):
        """Produkt der Einheiten mod 9 = 3^2 ist 8 ≡ -1 (mod 9)."""
        assert product_of_units(9) == 8

    def test_wilson_n10(self):
        """Produkt der Einheiten mod 10 = 2*5 ist 9 ≡ -1 (mod 10)."""
        # 10 = 2*5^1, also n = 2*p^k. Produkt = -1.
        assert product_of_units(10) == 9

    def test_wilson_n12(self):
        """Produkt der Einheiten mod 12 ist 1 (nicht -1), da 12 = 4*3."""
        # 12 = 2^2 * 3: nicht von der Form 4, p^k oder 2p^k
        assert product_of_units(12) == 1

    def test_wilson_n15(self):
        """Produkt der Einheiten mod 15 ist 1 (nicht -1), da 15 = 3*5."""
        # 15 = 3*5: zwei verschiedene ungerade Primfaktoren → Produkt = 1
        assert product_of_units(15) == 1

    def test_wilson_n16(self):
        """Produkt der Einheiten mod 16 = 2^4 ist 1 (nicht -1)."""
        assert product_of_units(16) == 1

    def test_wilson_n18(self):
        """Produkt der Einheiten mod 18 = 2*3^2 ist 17 ≡ -1 (mod 18)."""
        # 18 = 2*9 = 2*3^2: n = 2*p^k. Produkt = -1.
        assert product_of_units(18) == 17

    def test_gauss_generalisation_all_n_2_to_20(self):
        """Gaußsche Verallgemeinerung für alle n von 2 bis 20."""
        for n in range(2, 21):
            actual = product_of_units(n)
            expected = gauss_wilson_expected(n)
            assert actual == expected, \
                f"Gaußsche Verallgemeinerung verletzt für n={n}: " \
                f"Produkt={actual}, erwartet={expected}"


class TestHarmonicSum:
    """
    Tests für das harmonische Lemma von Wolstenholme:
    1 + 1/2 + ... + 1/(p-1) ≡ 0 (mod p) für alle Primzahlen p >= 3.
    """

    def test_harmonic_sum_p3(self):
        """Harmonische Summe mod 3: 1 + 1/2 = 1 + 2 = 3 ≡ 0 (mod 3)."""
        assert harmonic_sum_mod_p(3) == 0

    def test_harmonic_sum_p5(self):
        """Harmonische Summe mod 5: 1 + 3 + 2 + 4 = 10 ≡ 0 (mod 5)."""
        assert harmonic_sum_mod_p(5) == 0

    def test_harmonic_sum_p7(self):
        """Harmonische Summe mod 7."""
        assert harmonic_sum_mod_p(7) == 0

    def test_harmonic_sum_p11(self):
        """Harmonische Summe mod 11."""
        assert harmonic_sum_mod_p(11) == 0

    def test_harmonic_sum_p13(self):
        """Harmonische Summe mod 13."""
        assert harmonic_sum_mod_p(13) == 0

    def test_harmonic_sum_all_primes_up_to_50(self):
        """Wolstenholmes Lemma für alle Primzahlen von 3 bis 50."""
        for p in range(3, 51):
            if is_prime(p):
                assert harmonic_sum_mod_p(p) == 0, \
                    f"Harmonische Summe nicht null für p={p}"


class TestQuadraticResidue:
    """
    Tests für die Charakterisierung von -1 als quadratischem Rest modulo p.
    -1 ist QR mod p genau dann wenn p ≡ 1 (mod 4).
    """

    def test_minus1_qr_p5(self):
        """5 ≡ 1 (mod 4): -1 ist QR mod 5."""
        # 2^2 = 4 ≡ -1 (mod 5)
        assert 5 % 4 == 1
        assert pow(2, 2, 5) == 4  # = -1 mod 5
        assert is_quadratic_residue_minus1(5) is True

    def test_minus1_qr_p13(self):
        """13 ≡ 1 (mod 4): -1 ist QR mod 13."""
        assert 13 % 4 == 1
        assert is_quadratic_residue_minus1(13) is True

    def test_minus1_nonqr_p3(self):
        """3 ≡ 3 (mod 4): -1 ist KEIN QR mod 3."""
        assert 3 % 4 == 3
        assert is_quadratic_residue_minus1(3) is False

    def test_minus1_nonqr_p7(self):
        """7 ≡ 3 (mod 4): -1 ist KEIN QR mod 7."""
        assert 7 % 4 == 3
        assert is_quadratic_residue_minus1(7) is False

    def test_minus1_nonqr_p11(self):
        """11 ≡ 3 (mod 4): -1 ist KEIN QR mod 11."""
        assert is_quadratic_residue_minus1(11) is False

    def test_minus1_qr_p17(self):
        """17 ≡ 1 (mod 4): -1 ist QR mod 17."""
        assert is_quadratic_residue_minus1(17) is True
        # Überprüfe: gibt es x mit x^2 ≡ -1 (mod 17)?
        found = any(pow(x, 2, 17) == 16 for x in range(1, 17))
        assert found

    def test_minus1_qr_all_primes_up_to_50(self):
        """Wilson-Kriterium für -1 als QR stimmt mit p ≡ 1 (mod 4) überein."""
        for p in range(3, 51):
            if is_prime(p):
                via_wilson = is_quadratic_residue_minus1(p)
                via_congruence = (p % 4 == 1)
                assert via_wilson == via_congruence, \
                    f"Diskrepanz für p={p}: Wilson={via_wilson}, Kongruenz={via_congruence}"

    def test_wilson_double_count_formula(self):
        """
        Überprüft die Doppelzählungsformel aus Paper 12:
        (p-1)! ≡ M^2 * (-1)^((p-1)/2) (mod p), M = ((p-1)/2)!
        """
        for p in [5, 7, 11, 13, 17, 19, 23]:
            m = (p - 1) // 2
            M = math.factorial(m)
            M2 = (M * M) % p
            sign = (-1) ** m % p
            lhs = math.factorial(p - 1) % p  # = p-1 = -1 mod p
            rhs = (M2 * sign) % p
            assert lhs == rhs, \
                f"Doppelzählungsformel verletzt für p={p}: lhs={lhs}, rhs={rhs}"


class TestFermatViaWilson:
    """
    Tests für den Beweis des kleinen Satzes von Fermat via Wilson.
    a^(p-1) ≡ 1 (mod p) für gcd(a,p) = 1.
    """

    def test_fermat_p5_a2(self):
        """Fermatscher Satz: 2^4 ≡ 1 (mod 5)."""
        assert pow(2, 4, 5) == 1

    def test_fermat_p7_a3(self):
        """Fermatscher Satz: 3^6 ≡ 1 (mod 7)."""
        assert pow(3, 6, 7) == 1

    def test_fermat_p11_a5(self):
        """Fermatscher Satz: 5^10 ≡ 1 (mod 11)."""
        assert pow(5, 10, 11) == 1

    def test_permutation_argument(self):
        """
        Überprüft das Permutationsargument aus dem Fermat-Beweis:
        Die Vielfachen {a, 2a, ..., (p-1)a} sind eine Permutation von {1,...,p-1} mod p.
        """
        p = 7
        a = 3
        multiples = sorted((j * a) % p for j in range(1, p))
        assert multiples == list(range(1, p))

    def test_product_cancellation(self):
        """
        Überprüft den Kürzungsschritt im Fermat-Beweis via Wilson.
        a^(p-1) * (p-1)! ≡ (p-1)! (mod p) impliziert a^(p-1) ≡ 1.
        """
        for p in [5, 7, 11, 13]:
            for a in range(1, p):
                fact = math.factorial(p - 1) % p  # = p-1 = -1
                lhs = (pow(a, p - 1, p) * fact) % p
                rhs = fact
                assert lhs == rhs, \
                    f"Kürzungsschritt verletzt für p={p}, a={a}"
