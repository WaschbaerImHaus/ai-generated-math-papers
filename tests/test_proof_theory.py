"""
Tests für das Beweistheorie-Modul (proof_theory.py).

Testet alle implementierten Funktionen zu offenen Vermutungen,
Primzahltheorie und formalen Beweishilfsmitteln.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import pytest
import sys
import os

# Projektpfad für Import hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from proof_theory import (
    # Collatz
    collatz_sequence,
    collatz_stopping_time,
    collatz_max_value,
    collatz_verify_range,
    # Goldbach
    is_prime_fast,
    goldbach_decomposition,
    goldbach_all_decompositions,
    goldbach_verify_range,
    # Zwillingsprimzahlen
    find_twin_primes,
    twin_prime_count,
    # Riemann
    riemann_zeta_partial,
    riemann_zeta_euler_maclaurin,
    check_riemann_hypothesis_numerically,
    # Formale Beweise
    ProofByInduction,
    sieve_of_eratosthenes,
    miller_rabin_primality_test,
    # Zahlentheorie
    legendre_symbol,
    jacobi_symbol,
    chinese_remainder_theorem,
    conjecture_status_report
)


# ===========================================================================
# COLLATZ-VERMUTUNG TESTS
# ===========================================================================

class TestCollatz:
    """Tests zur Collatz-Vermutung und verwandten Funktionen."""

    def test_collatz_sequence_known_values(self):
        """Bekannte Collatz-Folgen überprüfen."""
        # n=1: sofort bei 1
        assert collatz_sequence(1) == [1]

        # n=2: 2 → 1
        assert collatz_sequence(2) == [2, 1]

        # n=3: 3 → 10 → 5 → 16 → 8 → 4 → 2 → 1
        expected_3 = [3, 10, 5, 16, 8, 4, 2, 1]
        assert collatz_sequence(3) == expected_3

        # n=6: 6 → 3 → 10 → 5 → 16 → 8 → 4 → 2 → 1
        assert collatz_sequence(6) == [6, 3, 10, 5, 16, 8, 4, 2, 1]

    def test_collatz_always_ends_at_one(self):
        """Collatz-Folge endet für kleine Zahlen immer bei 1."""
        for n in range(1, 100):
            seq = collatz_sequence(n)
            assert seq[-1] == 1, f"Folge für n={n} endet nicht bei 1"

    def test_collatz_sequence_invalid_input(self):
        """Ungültige Eingaben werfen ValueError."""
        with pytest.raises(ValueError):
            collatz_sequence(0)
        with pytest.raises(ValueError):
            collatz_sequence(-5)

    def test_collatz_stopping_time(self):
        """Stoppzeiten für bekannte Zahlen prüfen."""
        assert collatz_stopping_time(1) == 0   # schon bei 1
        assert collatz_stopping_time(2) == 1   # 2 → 1
        assert collatz_stopping_time(3) == 7   # 3 → ... → 1 (7 Schritte)

    def test_collatz_stopping_time_27(self):
        """n=27 ist berühmt für seine lange Stoppzeit (111 Schritte)."""
        assert collatz_stopping_time(27) == 111

    def test_collatz_max_value(self):
        """Maximalwerte in Collatz-Folgen."""
        assert collatz_max_value(1) == 1
        assert collatz_max_value(3) == 16       # 3 → 10 → 5 → 16 → ...
        assert collatz_max_value(27) == 9232    # berühmter Maximalwert

    def test_collatz_verify_range(self):
        """Verifikation für kleinen Bereich."""
        result = collatz_verify_range(100)
        assert result["verified"] is True
        assert result["counterexample"] is None
        assert result["range_checked"] == 100

    def test_collatz_verify_range_finds_long_sequences(self):
        """Verifikation findet Zahlen mit langen Folgen."""
        result = collatz_verify_range(30)
        # n=27 hat Stoppzeit 111 – muss als maximale Stoppzeit erkannt werden
        # (innerhalb von 1..30 ist 27 die langsamste)
        assert result["max_stopping_time"][0] == 27


# ===========================================================================
# GOLDBACH-VERMUTUNG TESTS
# ===========================================================================

class TestGoldbach:
    """Tests zur Goldbach-Vermutung."""

    def test_is_prime_fast_small_primes(self):
        """Kleine Primzahlen korrekt erkannt."""
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
        for p in primes:
            assert is_prime_fast(p), f"{p} sollte prim sein"

    def test_is_prime_fast_composites(self):
        """Zusammengesetzte Zahlen korrekt abgelehnt."""
        composites = [0, 1, 4, 6, 8, 9, 10, 12, 15, 25, 49, 100]
        for c in composites:
            assert not is_prime_fast(c), f"{c} sollte nicht prim sein"

    def test_goldbach_decomposition_basic(self):
        """Goldbach-Zerlegungen für bekannte Fälle."""
        # 4 = 2 + 2
        p, q = goldbach_decomposition(4)
        assert p + q == 4 and is_prime_fast(p) and is_prime_fast(q)

        # 6 = 3 + 3
        p, q = goldbach_decomposition(6)
        assert p + q == 6

        # 100 hat mindestens eine Zerlegung
        result = goldbach_decomposition(100)
        assert result is not None
        p, q = result
        assert p + q == 100 and is_prime_fast(p) and is_prime_fast(q)

    def test_goldbach_decomposition_invalid(self):
        """Ungültige Eingaben werfen ValueError."""
        with pytest.raises(ValueError):
            goldbach_decomposition(2)    # nicht > 2
        with pytest.raises(ValueError):
            goldbach_decomposition(7)    # ungerade
        with pytest.raises(ValueError):
            goldbach_decomposition(1)

    def test_goldbach_all_decompositions(self):
        """Alle Zerlegungen einer Zahl."""
        # 4 hat genau eine: (2,2)
        decomps = goldbach_all_decompositions(4)
        assert len(decomps) == 1
        assert (2, 2) in decomps

        # 10 hat zwei: (3,7) und (5,5)
        decomps = goldbach_all_decompositions(10)
        assert len(decomps) == 2
        assert (3, 7) in decomps
        assert (5, 5) in decomps

    def test_goldbach_verify_range(self):
        """Goldbach für kleinen Bereich verifizieren."""
        result = goldbach_verify_range(100)
        assert result["verified"] is True
        assert result["counterexample"] is None


# ===========================================================================
# ZWILLINGSPRIMZAHLEN TESTS
# ===========================================================================

class TestTwinPrimes:
    """Tests zu Zwillingsprimzahlpaaren."""

    def test_find_twin_primes_first_pairs(self):
        """Erste bekannte Zwillingsprimzahlpaare."""
        twins = find_twin_primes(50)
        expected = [(3, 5), (5, 7), (11, 13), (17, 19), (29, 31), (41, 43)]
        for pair in expected:
            assert pair in twins, f"Zwillingspaar {pair} fehlt"

    def test_twin_prime_count_increases(self):
        """Anzahl der Paare wächst mit der Schranke."""
        count_100 = twin_prime_count(100)
        count_1000 = twin_prime_count(1000)
        assert count_1000 > count_100

    def test_twin_primes_are_prime(self):
        """Alle gefundenen Paare sind tatsächlich Primzahlen."""
        twins = find_twin_primes(200)
        for p, q in twins:
            assert is_prime_fast(p), f"{p} ist nicht prim"
            assert is_prime_fast(q), f"{q} ist nicht prim"
            assert q == p + 2, f"Differenz ist nicht 2: ({p}, {q})"


# ===========================================================================
# RIEMANN-ZETA TESTS
# ===========================================================================

class TestRiemannZeta:
    """Tests zur Riemann-Zeta-Funktion."""

    def test_zeta_s2_pi_squared_over_6(self):
        """ζ(2) = π²/6 ≈ 1.6449... (Basler Problem, Euler 1734)."""
        import math
        zeta_2 = riemann_zeta_euler_maclaurin(complex(2, 0), terms=500)
        expected = math.pi ** 2 / 6
        assert abs(zeta_2.real - expected) < 0.01, (
            f"ζ(2) = {zeta_2.real}, erwartet {expected}"
        )

    def test_zeta_s4(self):
        """ζ(4) = π⁴/90 ≈ 1.0823..."""
        import math
        zeta_4 = riemann_zeta_euler_maclaurin(complex(4, 0), terms=200)
        expected = math.pi ** 4 / 90
        assert abs(zeta_4.real - expected) < 0.01

    def test_zeta_partial_sum_converges(self):
        """Partialsumme konvergiert für Re(s) > 1."""
        z = riemann_zeta_partial(complex(3, 0), terms=5000)
        import math
        expected = math.pi ** 6 / 945  # ζ(6) = π⁶/945
        # Weniger streng, da ζ(3) ≈ 1.202 (Apérys Konstante)
        assert 1.0 < z.real < 1.5

    def test_zeta_partial_sum_raises_for_small_s(self):
        """Partialsumme wirft Fehler für Re(s) ≤ 1."""
        with pytest.raises(ValueError):
            riemann_zeta_partial(complex(0.5, 14.0))  # auf kritischer Geraden

    def test_riemann_hypothesis_check_first_zeros(self):
        """Erste bekannte Nullstellen auf kritischer Geraden bestätigen."""
        result = check_riemann_hypothesis_numerically(
            imaginary_range=(10.0, 35.0),
            resolution=500
        )
        # Bekannte Nullstellen t₁ ≈ 14.13, t₂ ≈ 21.02, t₃ ≈ 25.01
        known = result["known_zeros_in_range"]
        assert len(known) >= 3, "Mindestens 3 bekannte Nullstellen sollten im Bereich liegen"


# ===========================================================================
# SIEB DES ERATOSTHENES UND MILLER-RABIN TESTS
# ===========================================================================

class TestPrimality:
    """Tests zu Primzahltests und dem Sieb."""

    def test_sieve_of_eratosthenes_first_primes(self):
        """Sieb gibt korrekte Primzahlen zurück."""
        primes = sieve_of_eratosthenes(30)
        assert primes == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    def test_sieve_empty_for_small_limit(self):
        """Sieb gibt leere Liste für limit < 2 zurück."""
        assert sieve_of_eratosthenes(0) == []
        assert sieve_of_eratosthenes(1) == []

    def test_sieve_count_matches_known(self):
        """Anzahl der Primzahlen bis 100 ist 25 (bekannter Wert)."""
        primes = sieve_of_eratosthenes(100)
        assert len(primes) == 25

    def test_miller_rabin_primes(self):
        """Miller-Rabin erkennt Primzahlen korrekt."""
        known_primes = [2, 3, 5, 7, 11, 13, 17, 19, 97, 1009, 104729]
        for p in known_primes:
            assert miller_rabin_primality_test(p), f"{p} sollte prim sein"

    def test_miller_rabin_composites(self):
        """Miller-Rabin erkennt zusammengesetzte Zahlen."""
        composites = [4, 6, 8, 9, 15, 100, 561, 1105, 8911]
        for c in composites:
            assert not miller_rabin_primality_test(c), f"{c} sollte nicht prim sein"

    def test_miller_rabin_large_prime(self):
        """Miller-Rabin für große Primzahl (Mersenne-Primzahl: 2^31 - 1 = 2147483647)."""
        mersenne_31 = 2 ** 31 - 1
        assert miller_rabin_primality_test(mersenne_31)

    def test_sieve_and_miller_rabin_agree(self):
        """Sieb und Miller-Rabin liefern identische Ergebnisse bis 1000."""
        sieve_primes = set(sieve_of_eratosthenes(1000))
        miller_primes = set(n for n in range(2, 1001) if miller_rabin_primality_test(n))
        assert sieve_primes == miller_primes


# ===========================================================================
# ZAHLENTHEORETISCHE TESTS
# ===========================================================================

class TestNumberTheory:
    """Tests für zahlentheoretische Hilfsfunktionen."""

    def test_legendre_symbol_quadratic_residues(self):
        """Legendre-Symbol für quadratische Reste."""
        # Quadratische Reste mod 7: 1²=1, 2²=4, 3²=2 → {1, 2, 4}
        assert legendre_symbol(1, 7) == 1
        assert legendre_symbol(2, 7) == 1
        assert legendre_symbol(4, 7) == 1
        # Nichtrest: 3, 5, 6
        assert legendre_symbol(3, 7) == -1
        assert legendre_symbol(5, 7) == -1
        # Null: 7 teilt 7
        assert legendre_symbol(7, 7) == 0

    def test_legendre_symbol_invalid(self):
        """Legendre-Symbol wirft Fehler für nicht-prim oder gerade p."""
        with pytest.raises(ValueError):
            legendre_symbol(3, 4)    # 4 ist nicht prim
        with pytest.raises(ValueError):
            legendre_symbol(1, 2)    # 2 ist gerade

    def test_jacobi_symbol_basic(self):
        """Jacobi-Symbol für einfache Fälle."""
        # (a/p) = Jacobi = Legendre wenn p prim
        assert jacobi_symbol(2, 7) == 1    # wie Legendre
        assert jacobi_symbol(3, 7) == -1
        # (a/n) = 0 wenn gcd(a,n) > 1
        assert jacobi_symbol(3, 9) == 0

    def test_chinese_remainder_theorem_basic(self):
        """CRT löst einfache simultane Kongruenzen."""
        # x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
        # Lösung: x = 23
        result = chinese_remainder_theorem([2, 3, 2], [3, 5, 7])
        assert result == 23
        assert result % 3 == 2
        assert result % 5 == 3
        assert result % 7 == 2

    def test_chinese_remainder_theorem_single(self):
        """CRT mit einer Kongruenz."""
        result = chinese_remainder_theorem([5], [11])
        assert result == 5

    def test_chinese_remainder_theorem_invalid(self):
        """CRT wirft Fehler bei nicht-koprimen Moduln."""
        with pytest.raises(ValueError):
            chinese_remainder_theorem([1, 2], [4, 6])   # gcd(4,6) = 2 ≠ 1

    def test_chinese_remainder_theorem_mismatched_lengths(self):
        """CRT wirft Fehler bei falscher Länge."""
        with pytest.raises(ValueError):
            chinese_remainder_theorem([1, 2, 3], [3, 5])


# ===========================================================================
# BEWEISSTRUKTUR TESTS
# ===========================================================================

class TestProofByInduction:
    """Tests für die Induktions-Beweishilfe."""

    def test_verify_base_case_gauss(self):
        """Basisfall für Gaußsche Summenformel: Σk = n(n+1)/2."""
        gauss = lambda n: sum(range(1, n + 1)) == n * (n + 1) // 2
        assert ProofByInduction.verify_base_case(gauss, base=1) is True

    def test_empirical_verify_gauss(self):
        """Empirische Verifikation der Gaußschen Summe bis 1000."""
        gauss = lambda n: sum(range(1, n + 1)) == n * (n + 1) // 2
        result = ProofByInduction.empirical_verify(gauss, start=1, end=1000)
        assert result["verified"] is True
        assert result["counterexample"] is None

    def test_empirical_verify_finds_counterexample(self):
        """Empirische Verifikation findet Gegenbeispiel bei falscher Aussage."""
        # Falsche Aussage: n² > 2n für alle n ≥ 1 (falsch für n=1,2)
        false_predicate = lambda n: n ** 2 > 2 * n
        result = ProofByInduction.empirical_verify(false_predicate, start=1, end=10)
        assert result["verified"] is False
        assert result["counterexample"] is not None
        # Gegenbeispiele: n=1 (1 > 2 falsch), n=2 (4 > 4 falsch)
        assert result["counterexample"] <= 2


# ===========================================================================
# STATUS-BERICHT TESTS
# ===========================================================================

class TestConjectureStatus:
    """Tests für den Vermutungsstatus-Bericht."""

    def test_status_report_contains_all_millennium_problems(self):
        """Statusbericht enthält alle 7 Millennium-Probleme."""
        report = conjecture_status_report()
        millennium = report["millennium_problems"]
        assert len(millennium) == 7

    def test_poincare_marked_as_solved(self):
        """Poincaré-Vermutung als gelöst markiert."""
        report = conjecture_status_report()
        poincare = report["millennium_problems"]["Poincaré-Vermutung"]
        assert "bewiesen" in poincare["status"]

    def test_status_report_has_other_conjectures(self):
        """Statusbericht enthält weitere bekannte Vermutungen."""
        report = conjecture_status_report()
        other = report["other_conjectures"]
        assert "Goldbach-Vermutung" in other
        assert "Collatz-Vermutung" in other
        assert "Zwillingsprimzahl" in other
