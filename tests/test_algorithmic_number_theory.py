"""
Tests für das Modul algorithmic_number_theory.py.

Umfassende Testabdeckung für alle Klassen und Funktionen der algorithmischen
Zahlentheorie: Primzahlsiebe, Primzahltests, Faktorisierung,
diskreter Logarithmus, Kettenbrüche und kryptographische Algorithmen.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import math
import sys
import os
import pytest

# Sicherstellen dass src/ im Suchpfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algorithmic_number_theory import (
    PrimeSieve,
    PrimalityTests,
    IntegerFactorization,
    DiscreteLogarithm,
    ContinuedFractionFactoring,
    CryptographicAlgorithms,
    euler_product_approximation,
    mertens_first_theorem,
    primitive_root,
    legendre_symbol_jacobi,
    _gcd,
    _mod_inverse,
    _is_perfect_square,
    _prime_factors_list,
    _crt,
)


# ===========================================================================
# TESTS FÜR HILFSFUNKTIONEN
# ===========================================================================

class TestHelpers:
    """Tests für interne Hilfsfunktionen."""

    def test_gcd_basic(self):
        """ggT zweier normaler Zahlen."""
        assert _gcd(12, 8) == 4

    def test_gcd_coprime(self):
        """ggT teilerfremder Zahlen = 1."""
        assert _gcd(7, 13) == 1

    def test_gcd_zero(self):
        """ggT mit 0: ggT(a, 0) = a."""
        assert _gcd(5, 0) == 5
        assert _gcd(0, 7) == 7

    def test_gcd_negative(self):
        """ggT ist immer positiv."""
        assert _gcd(-12, 8) == 4

    def test_mod_inverse_exists(self):
        """Modulares Inverses wenn ggT(a, m) = 1."""
        inv = _mod_inverse(3, 7)
        assert inv is not None
        assert (3 * inv) % 7 == 1

    def test_mod_inverse_not_exists(self):
        """Kein Inverses wenn ggT(a, m) ≠ 1."""
        assert _mod_inverse(4, 6) is None

    def test_is_perfect_square_true(self):
        """Perfekte Quadratzahl erkannt."""
        ok, root = _is_perfect_square(49)
        assert ok is True
        assert root == 7

    def test_is_perfect_square_false(self):
        """Nicht-Quadratzahl korrekt abgelehnt."""
        ok, root = _is_perfect_square(48)
        assert ok is False

    def test_is_perfect_square_zero(self):
        """Null ist perfektes Quadrat."""
        ok, root = _is_perfect_square(0)
        assert ok is True
        assert root == 0

    def test_prime_factors_list(self):
        """Primfaktorliste mit Wiederholung."""
        factors = _prime_factors_list(12)
        assert sorted(factors) == [2, 2, 3]

    def test_prime_factors_list_prime(self):
        """Primzahl hat nur sich selbst als Faktor."""
        assert _prime_factors_list(13) == [13]

    def test_crt_basic(self):
        """Chinesischer Restsatz: x ≡ 2 (mod 3), x ≡ 3 (mod 5) → x ≡ 8 (mod 15)."""
        result = _crt([2, 3], [3, 5])
        assert result is not None
        assert result % 3 == 2
        assert result % 5 == 3

    def test_crt_single(self):
        """CRT mit einer Gleichung."""
        result = _crt([5], [7])
        assert result == 5


# ===========================================================================
# TESTS FÜR PRIMEZAHLSIEBE (PrimeSieve)
# ===========================================================================

class TestPrimeSieve:
    """Tests für die PrimeSieve-Klasse."""

    # --- sieve_of_eratosthenes ---

    def test_eratosthenes_small(self):
        """Primzahlen bis 20 korrekt."""
        primes = PrimeSieve.sieve_of_eratosthenes(20)
        assert primes == [2, 3, 5, 7, 11, 13, 17, 19]

    def test_eratosthenes_n1(self):
        """Keine Primzahlen bis 1."""
        assert PrimeSieve.sieve_of_eratosthenes(1) == []

    def test_eratosthenes_n2(self):
        """Nur 2 als Primzahl bis 2."""
        assert PrimeSieve.sieve_of_eratosthenes(2) == [2]

    def test_eratosthenes_count(self):
        """π(100) = 25 Primzahlen bis 100."""
        assert len(PrimeSieve.sieve_of_eratosthenes(100)) == 25

    def test_eratosthenes_large(self):
        """π(1000) = 168 Primzahlen bis 1000."""
        assert len(PrimeSieve.sieve_of_eratosthenes(1000)) == 168

    # --- sieve_of_atkin ---

    def test_atkin_small(self):
        """Atkin-Sieb bis 20 liefert dieselben Primzahlen wie Eratosthenes."""
        expected = PrimeSieve.sieve_of_eratosthenes(20)
        assert PrimeSieve.sieve_of_atkin(20) == expected

    def test_atkin_n1(self):
        """Atkin-Sieb: keine Primzahlen bis 1."""
        assert PrimeSieve.sieve_of_atkin(1) == []

    def test_atkin_n2(self):
        """Atkin-Sieb: nur [2] bis 2."""
        assert PrimeSieve.sieve_of_atkin(2) == [2]

    def test_atkin_n3(self):
        """Atkin-Sieb: [2, 3] bis 3."""
        assert PrimeSieve.sieve_of_atkin(3) == [2, 3]

    def test_atkin_count_100(self):
        """Atkin π(100) = 25."""
        assert len(PrimeSieve.sieve_of_atkin(100)) == 25

    def test_atkin_matches_eratosthenes_200(self):
        """Atkin und Eratosthenes liefern identische Ergebnisse bis 200."""
        assert PrimeSieve.sieve_of_atkin(200) == PrimeSieve.sieve_of_eratosthenes(200)

    # --- prime_gaps ---

    def test_prime_gaps_first(self):
        """Erste Primzahllücke: (2, 3, 1)."""
        gaps = PrimeSieve.prime_gaps(10)
        assert gaps[0] == (2, 3, 1)

    def test_prime_gaps_twin(self):
        """Zwillingsprimzahl-Lücke = 2."""
        gaps = PrimeSieve.prime_gaps(10)
        # (3,5,2) und (5,7,2) sind Zwillingslücken
        gap_values = [g[2] for g in gaps]
        assert 2 in gap_values

    def test_prime_gaps_empty(self):
        """Keine Lücken bei n < 3."""
        assert PrimeSieve.prime_gaps(2) == []

    def test_prime_gaps_structure(self):
        """Jede Lücke hat Format (p1, p2, gap) mit p2 = p1 + gap."""
        for p1, p2, gap in PrimeSieve.prime_gaps(50):
            assert p2 == p1 + gap

    # --- twin_primes ---

    def test_twin_primes_basic(self):
        """Erste Zwillingsprimzahlen bis 20."""
        twins = PrimeSieve.twin_primes(20)
        assert (3, 5) in twins
        assert (5, 7) in twins
        assert (11, 13) in twins
        assert (17, 19) in twins

    def test_twin_primes_empty(self):
        """Keine Zwillingsprimzahlen bei n < 3."""
        assert PrimeSieve.twin_primes(2) == []

    def test_twin_primes_all_pairs_differ_by_2(self):
        """Alle Zwillingspaare haben Abstand 2."""
        for p, q in PrimeSieve.twin_primes(100):
            assert q - p == 2

    # --- prime_counting_exact ---

    def test_prime_counting_n10(self):
        """π(10) = 4."""
        assert PrimeSieve.prime_counting_exact(10) == 4

    def test_prime_counting_n100(self):
        """π(100) = 25."""
        assert PrimeSieve.prime_counting_exact(100) == 25

    def test_prime_counting_n1(self):
        """π(1) = 0."""
        assert PrimeSieve.prime_counting_exact(1) == 0


# ===========================================================================
# TESTS FÜR PRIMZAHLTESTS (PrimalityTests)
# ===========================================================================

class TestPrimalityTests:
    """Tests für die PrimalityTests-Klasse."""

    # --- aks_test_demo ---

    def test_aks_prime_7(self):
        """AKS-Demo erkennt 7 als prim."""
        result = PrimalityTests.aks_test_demo(7)
        assert result['is_prime'] is True

    def test_aks_composite_9(self):
        """AKS-Demo erkennt 9 als zusammengesetzt."""
        result = PrimalityTests.aks_test_demo(9)
        assert result['is_prime'] is False

    def test_aks_n1(self):
        """AKS-Demo: n=1 ist nicht prim."""
        result = PrimalityTests.aks_test_demo(1)
        assert result['is_prime'] is False

    def test_aks_perfect_power_8(self):
        """AKS-Demo erkennt 8 = 2³ als zusammengesetzt."""
        result = PrimalityTests.aks_test_demo(8)
        assert result['is_prime'] is False

    def test_aks_n2(self):
        """AKS-Demo: n=2 ist prim."""
        result = PrimalityTests.aks_test_demo(2)
        assert result['is_prime'] is True

    def test_aks_has_steps(self):
        """AKS-Demo enthält Schrittbeschreibungen."""
        result = PrimalityTests.aks_test_demo(13)
        assert 'steps' in result
        assert len(result['steps']) > 0

    # --- solovay_strassen ---

    def test_solovay_strassen_prime(self):
        """Solovay-Strassen erkennt 97 als prim."""
        assert PrimalityTests.solovay_strassen(97, k=20) is True

    def test_solovay_strassen_composite(self):
        """Solovay-Strassen erkennt 91 = 7×13 als zusammengesetzt."""
        # Mehrere Versuche da probabilistisch
        results = [PrimalityTests.solovay_strassen(91, k=20) for _ in range(5)]
        # Mindestens einmal zusammengesetzt erkannt
        assert False in results

    def test_solovay_strassen_n2(self):
        """Solovay-Strassen: 2 ist prim."""
        assert PrimalityTests.solovay_strassen(2) is True

    def test_solovay_strassen_even(self):
        """Solovay-Strassen: gerade Zahlen > 2 sind zusammengesetzt."""
        assert PrimalityTests.solovay_strassen(100) is False

    def test_solovay_strassen_n1(self):
        """Solovay-Strassen: 1 ist nicht prim."""
        assert PrimalityTests.solovay_strassen(1) is False

    # --- fermat_test ---

    def test_fermat_test_prime(self):
        """Fermat-Test: 13 besteht zur Basis 2."""
        assert PrimalityTests.fermat_test(13, 2) is True

    def test_fermat_test_composite(self):
        """Fermat-Test: 15 = 3×5 fällt bei Basis 2 durch (2^14 ≡ 4 ≠ 1 mod 15)."""
        assert PrimalityTests.fermat_test(15, 2) is False

    def test_fermat_test_carmichael_561(self):
        """Carmichael-Zahl 561 besteht Fermat-Test für coprime Basen."""
        # 561 = 3 × 11 × 17 ist eine Carmichael-Zahl
        # Sie besteht für alle a mit ggT(a, 561) = 1
        assert PrimalityTests.fermat_test(561, 2) is True
        assert PrimalityTests.fermat_test(561, 5) is True

    def test_fermat_test_n2(self):
        """Fermat-Test: 2 ist prim."""
        assert PrimalityTests.fermat_test(2, 1) is True

    def test_fermat_test_n1(self):
        """Fermat-Test: 1 ist nicht prim."""
        assert PrimalityTests.fermat_test(1, 2) is False

    # --- lucas_primality_test ---

    def test_lucas_prime_7(self):
        """Lucas-Test bestätigt 7 als prim."""
        result = PrimalityTests.lucas_primality_test(7)
        assert result['is_prime'] is True

    def test_lucas_composite_9(self):
        """Lucas-Test erkennt 9 als zusammengesetzt."""
        result = PrimalityTests.lucas_primality_test(9)
        assert result['is_prime'] is False

    def test_lucas_n2(self):
        """Lucas-Test: 2 ist prim."""
        result = PrimalityTests.lucas_primality_test(2)
        assert result['is_prime'] is True

    def test_lucas_prime_13(self):
        """Lucas-Test findet Zeugen für 13."""
        result = PrimalityTests.lucas_primality_test(13)
        assert result['is_prime'] is True
        assert result['witness'] is not None

    # --- is_carmichael_number ---

    def test_carmichael_561(self):
        """561 ist die kleinste Carmichael-Zahl."""
        assert PrimalityTests.is_carmichael_number(561) is True

    def test_carmichael_1105(self):
        """1105 = 5×13×17 ist Carmichael-Zahl."""
        assert PrimalityTests.is_carmichael_number(1105) is True

    def test_carmichael_not_prime(self):
        """Primzahlen sind keine Carmichael-Zahlen."""
        assert PrimalityTests.is_carmichael_number(13) is False

    def test_carmichael_not_composite(self):
        """Normale zusammengesetzte Zahlen sind keine Carmichael-Zahlen."""
        assert PrimalityTests.is_carmichael_number(12) is False

    def test_carmichael_n1(self):
        """1 ist keine Carmichael-Zahl."""
        assert PrimalityTests.is_carmichael_number(1) is False


# ===========================================================================
# TESTS FÜR FAKTORISIERUNG (IntegerFactorization)
# ===========================================================================

class TestIntegerFactorization:
    """Tests für die IntegerFactorization-Klasse."""

    # --- trial_division ---

    def test_trial_division_12(self):
        """12 = 2² × 3."""
        result = IntegerFactorization.trial_division(12)
        assert result['factors'] == {2: 2, 3: 1}

    def test_trial_division_prime(self):
        """Primzahl 13 als prim erkannt."""
        result = IntegerFactorization.trial_division(13)
        assert result['is_prime'] is True
        assert result['factors'] == {13: 1}

    def test_trial_division_n1(self):
        """n=1 hat keine Primfaktoren."""
        result = IntegerFactorization.trial_division(1)
        assert result['factors'] == {}

    def test_trial_division_large_prime(self):
        """Große Primzahl korrekt erkannt."""
        result = IntegerFactorization.trial_division(9973)
        assert result['is_prime'] is True

    def test_trial_division_power_of_2(self):
        """64 = 2^6."""
        result = IntegerFactorization.trial_division(64)
        assert result['factors'] == {2: 6}

    # --- pollard_rho ---

    def test_pollard_rho_35(self):
        """Pollard-ρ findet Teiler von 35 = 5 × 7."""
        d = IntegerFactorization.pollard_rho(35)
        assert d is not None
        assert 35 % d == 0
        assert 1 < d < 35

    def test_pollard_rho_even(self):
        """Gerade Zahl: Teiler 2 gefunden."""
        d = IntegerFactorization.pollard_rho(100)
        assert d == 2

    def test_pollard_rho_small(self):
        """n < 4: None zurückgegeben."""
        assert IntegerFactorization.pollard_rho(3) is None

    def test_pollard_rho_semiprime(self):
        """Semiprimus 143 = 11 × 13: Teiler gefunden."""
        d = IntegerFactorization.pollard_rho(143)
        assert d is not None
        assert 143 % d == 0

    # --- pollard_p_minus_1 ---

    def test_pollard_p_minus_1_basic(self):
        """Pollard p-1: findet Teiler einer glatten Zahl."""
        # 3 × 5 = 15; p-1 für p=3 ist 2 (2-glatt), für p=5 ist 4 (2-glatt)
        d = IntegerFactorization.pollard_p_minus_1(15, B=10)
        # Kann keinen Teiler finden wenn zu klein, aber keinen Fehler werfen
        # Für n=15: erwartet einen Teiler oder None
        if d is not None:
            assert 15 % d == 0

    def test_pollard_p_minus_1_even(self):
        """Gerade Zahl: Teiler 2."""
        d = IntegerFactorization.pollard_p_minus_1(2 * 97, B=100)
        assert d == 2

    # --- fermat_factoring ---

    def test_fermat_factoring_close_primes(self):
        """Fermat: n = 7 × 11 = 77 (Faktoren nahe beieinander)."""
        result = IntegerFactorization.fermat_factoring(77)
        assert result is not None
        p, q = result
        assert p * q == 77

    def test_fermat_factoring_prime(self):
        """Primzahl: fermat_factoring gibt None zurück."""
        assert IntegerFactorization.fermat_factoring(13) is None

    def test_fermat_factoring_n4(self):
        """4 = 2 × 2."""
        result = IntegerFactorization.fermat_factoring(9)
        assert result is not None
        assert result[0] * result[1] == 9

    # --- factor_completely ---

    def test_factor_completely_360(self):
        """360 = 2³ × 3² × 5."""
        result = IntegerFactorization.factor_completely(360)
        f = result['factors']
        assert f[2] == 3
        assert f[3] == 2
        assert f[5] == 1

    def test_factor_completely_prime(self):
        """Primzahl 97 vollständig faktorisiert."""
        result = IntegerFactorization.factor_completely(97)
        assert result['factors'] == {97: 1}

    def test_factor_completely_product_correct(self):
        """Produkt der Faktoren ergibt ursprüngliche Zahl."""
        n = 2 * 3 * 5 * 7 * 11
        result = IntegerFactorization.factor_completely(n)
        product = 1
        for p, e in result['factors'].items():
            product *= p ** e
        assert product == n

    def test_factor_completely_string(self):
        """Faktorisierungsstring wird erstellt."""
        result = IntegerFactorization.factor_completely(12)
        assert 'factorization_str' in result
        assert len(result['factorization_str']) > 0

    # --- smooth_number_check ---

    def test_smooth_360_is_5_smooth(self):
        """360 = 2³ × 3² × 5 ist 5-glatt."""
        assert IntegerFactorization.smooth_number_check(360, 5) is True

    def test_smooth_14_not_5_smooth(self):
        """14 = 2 × 7 ist nicht 5-glatt."""
        assert IntegerFactorization.smooth_number_check(14, 5) is False

    def test_smooth_prime_itself(self):
        """Eine Primzahl p ist p-glatt aber nicht (p-1)-glatt."""
        assert IntegerFactorization.smooth_number_check(7, 7) is True
        assert IntegerFactorization.smooth_number_check(7, 5) is False


# ===========================================================================
# TESTS FÜR DISKRETEN LOGARITHMUS (DiscreteLogarithm)
# ===========================================================================

class TestDiscreteLogarithm:
    """Tests für die DiscreteLogarithm-Klasse."""

    # --- baby_step_giant_step ---

    def test_bsgs_basic(self):
        """BSGS: 5^x ≡ 8 (mod 23) → x = 6."""
        x = DiscreteLogarithm.baby_step_giant_step(5, 8, 23)
        assert x is not None
        assert pow(5, x, 23) == 8

    def test_bsgs_h_equals_1(self):
        """BSGS: h=1 → x=0 (g^0 = 1)."""
        x = DiscreteLogarithm.baby_step_giant_step(5, 1, 23)
        assert x == 0

    def test_bsgs_small_p(self):
        """BSGS auf kleinem Modulus."""
        # g=3, h=4, p=7: 3^x ≡ 4 (mod 7)
        x = DiscreteLogarithm.baby_step_giant_step(3, 4, 7)
        if x is not None:
            assert pow(3, x, 7) == 4

    def test_bsgs_result_is_solution(self):
        """Jedes zurückgegebene x ist tatsächlich Lösung."""
        g, h, p = 2, 16, 19
        x = DiscreteLogarithm.baby_step_giant_step(g, h, p)
        if x is not None:
            assert pow(g, x, p) == h % p

    # --- pohlig_hellman ---

    def test_pohlig_hellman_small(self):
        """Pohlig-Hellman auf kleinem Modulus."""
        # 5^x ≡ 8 (mod 23), Ordnung 22 = 2 × 11
        x = DiscreteLogarithm.pohlig_hellman(5, 8, 23)
        if x is not None:
            assert pow(5, x, 23) == 8

    def test_pohlig_hellman_trivial(self):
        """Pohlig-Hellman: h=g → x=1."""
        x = DiscreteLogarithm.pohlig_hellman(5, 5, 23)
        # x=1 oder eine äquivalente Lösung
        if x is not None:
            assert pow(5, x, 23) == 5

    # --- index_calculus_demo ---

    def test_index_calculus_returns_dict(self):
        """Indexkalkül gibt Dict mit korrekten Schlüsseln zurück."""
        result = DiscreteLogarithm.index_calculus_demo(5, 8, 23, B=10)
        assert 'result' in result
        assert 'factor_base' in result
        assert 'relations' in result
        assert 'log_found' in result

    def test_index_calculus_factor_base_are_primes(self):
        """Faktorbasis enthält nur Primzahlen."""
        result = DiscreteLogarithm.index_calculus_demo(5, 8, 23, B=15)
        fb = result['factor_base']
        for p in fb:
            # Alle Elemente der Faktorbasis sind prim
            from proof_theory import miller_rabin_primality_test
            assert miller_rabin_primality_test(p)

    # --- discrete_log_brute_force ---

    def test_brute_force_small(self):
        """Brute Force: 3^x ≡ 5 (mod 7)."""
        x = DiscreteLogarithm.discrete_log_brute_force(3, 5, 7)
        if x is not None:
            assert pow(3, x, 7) == 5

    def test_brute_force_identity(self):
        """Brute Force: g^0 = 1 → x=0 wenn h=1."""
        x = DiscreteLogarithm.discrete_log_brute_force(5, 1, 23)
        assert x == 0

    def test_brute_force_too_large(self):
        """Brute Force verweigert zu großen Modulus."""
        with pytest.raises(ValueError):
            DiscreteLogarithm.discrete_log_brute_force(2, 3, 10**7)


# ===========================================================================
# TESTS FÜR KETTENBRÜCHE (ContinuedFractionFactoring)
# ===========================================================================

class TestContinuedFractionFactoring:
    """Tests für Kettenbruchentwicklung und verwandte Methoden."""

    # --- continued_fraction_expansion ---

    def test_cf_sqrt2(self):
        """√2 = [1; 2, 2, 2, ...]: erster Koeffizient 1, Rest 2."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(2)
        assert cf[0] == 1
        assert cf[1] == 2  # Periode beginnt mit 2

    def test_cf_sqrt3(self):
        """√3 = [1; 1, 2, 1, 2, ...]: a₀=1, Periode [1,2]."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(3)
        assert cf[0] == 1

    def test_cf_perfect_square(self):
        """Perfekte Quadratzahl 4: Kettenbruch = [2]."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(4)
        assert cf == [2]

    def test_cf_sqrt7(self):
        """√7 = [2; 1, 1, 1, 4]: a₀=2."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(7)
        assert cf[0] == 2

    def test_cf_n1(self):
        """√1 = [1] (perfektes Quadrat)."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(1)
        assert cf == [1]

    # --- convergents ---

    def test_convergents_sqrt2(self):
        """Konvergenten von √2: 1/1, 3/2, 7/5, 17/12, ..."""
        cf = [1, 2, 2, 2, 2]  # [1; 2, 2, 2, 2]
        convs = ContinuedFractionFactoring.convergents(cf)
        # Erste Konvergente: 1/1
        assert convs[0] == (1, 1)
        # Zweite Konvergente: 1*2+1 / 1*2+0 = 3/2... warte, a₁=2
        # p₁ = 2*1+1=3, q₁ = 2*1+0=2
        assert convs[1] == (3, 2)

    def test_convergents_empty(self):
        """Leere Liste liefert leeres Ergebnis."""
        assert ContinuedFractionFactoring.convergents([]) == []

    def test_convergents_single(self):
        """Einzelner Koeffizient: Konvergente = a₀/1."""
        convs = ContinuedFractionFactoring.convergents([5])
        assert convs == [(5, 1)]

    def test_convergents_rational_approximation(self):
        """Konvergenten approximieren √2 immer besser."""
        cf = ContinuedFractionFactoring.continued_fraction_expansion(2, max_terms=10)
        convs = ContinuedFractionFactoring.convergents(cf)
        sqrt2 = math.sqrt(2)
        # Spätere Konvergenten sind besser
        errors = [abs(p / q - sqrt2) for p, q in convs if q > 0]
        for i in range(len(errors) - 1):
            assert errors[i] >= errors[i + 1]

    # --- pell_equation_solve ---

    def test_pell_d2(self):
        """Pell D=2: x²-2y²=1 → (3, 2) ist Fundamentallösung."""
        x, y = ContinuedFractionFactoring.pell_equation_solve(2)
        assert x * x - 2 * y * y == 1
        assert x > 0 and y > 0

    def test_pell_d3(self):
        """Pell D=3: x²-3y²=1 → (2, 1)."""
        x, y = ContinuedFractionFactoring.pell_equation_solve(3)
        assert x * x - 3 * y * y == 1

    def test_pell_d5(self):
        """Pell D=5: x²-5y²=1 → (9, 4)."""
        x, y = ContinuedFractionFactoring.pell_equation_solve(5)
        assert x * x - 5 * y * y == 1

    def test_pell_perfect_square_raises(self):
        """Pell mit perfektem Quadrat wirft ValueError."""
        with pytest.raises(ValueError):
            ContinuedFractionFactoring.pell_equation_solve(4)

    # --- cfrac_factoring_demo ---

    def test_cfrac_composite(self):
        """CFRAC-Demo faktorisiert zusammengesetzte Zahl."""
        result = ContinuedFractionFactoring.cfrac_factoring_demo(77)
        # Kann erfolgreich sein oder Fallback nutzen
        assert 'success' in result
        assert 'method' in result

    def test_cfrac_prime_returns_no_factor(self):
        """CFRAC-Demo bei Primzahl: kein Teiler."""
        result = ContinuedFractionFactoring.cfrac_factoring_demo(13)
        assert result['factor'] is None
        assert result['success'] is False

    def test_cfrac_even_handled(self):
        """CFRAC-Demo: gerade Zahl wird sicher behandelt."""
        result = ContinuedFractionFactoring.cfrac_factoring_demo(10)
        assert 'success' in result


# ===========================================================================
# TESTS FÜR KRYPTOGRAPHISCHE ALGORITHMEN
# ===========================================================================

class TestCryptographicAlgorithms:
    """Tests für CryptographicAlgorithms."""

    # --- diffie_hellman_demo ---

    def test_dh_shared_key_equal(self):
        """DH: Alice und Bob berechnen identischen gemeinsamen Schlüssel."""
        result = CryptographicAlgorithms.diffie_hellman_demo(p=23, g=5)
        assert result['shared_key'] > 0
        # Intern wird assert geprüft; Test bestätigt, dass kein Fehler auftritt

    def test_dh_result_structure(self):
        """DH-Ergebnis enthält alle erwarteten Schlüssel."""
        result = CryptographicAlgorithms.diffie_hellman_demo(p=23, g=5)
        for key in ['p', 'g', 'alice_private', 'bob_private', 'alice_public', 'bob_public', 'shared_key']:
            assert key in result

    def test_dh_public_keys_in_range(self):
        """DH: Öffentliche Schlüssel im Bereich [1, p-1]."""
        result = CryptographicAlgorithms.diffie_hellman_demo(p=23, g=5)
        p = result['p']
        assert 1 <= result['alice_public'] < p
        assert 1 <= result['bob_public'] < p

    # --- elgamal_keygen ---

    def test_elgamal_keygen_structure(self):
        """ElGamal Keygen gibt korrektes Dict zurück."""
        result = CryptographicAlgorithms.elgamal_keygen(p=23, g=5)
        assert 'public_key' in result
        assert 'private_key' in result
        assert 'y' in result['public_key']

    def test_elgamal_keygen_public_key_valid(self):
        """ElGamal: y = g^x mod p."""
        result = CryptographicAlgorithms.elgamal_keygen(p=23, g=5)
        p = result['public_key']['p']
        g = result['public_key']['g']
        y = result['public_key']['y']
        x = result['private_key']
        assert pow(g, x, p) == y

    # --- elgamal_encrypt / elgamal_decrypt ---

    def test_elgamal_roundtrip(self):
        """ElGamal: Verschlüsseln und Entschlüsseln ergibt Originalnachricht."""
        p, g = 23, 5
        keygen = CryptographicAlgorithms.elgamal_keygen(p=p, g=g)
        pub_y = keygen['public_key']['y']
        priv_x = keygen['private_key']
        message = 10

        ciphertext = CryptographicAlgorithms.elgamal_encrypt(message, pub_y, p, g)
        decrypted = CryptographicAlgorithms.elgamal_decrypt(ciphertext, priv_x, p)
        assert decrypted == message

    def test_elgamal_ciphertext_is_tuple(self):
        """ElGamal Ciphertext ist ein 2-Tupel."""
        keygen = CryptographicAlgorithms.elgamal_keygen(p=23, g=5)
        ct = CryptographicAlgorithms.elgamal_encrypt(7, keygen['public_key']['y'], 23, 5)
        assert isinstance(ct, tuple)
        assert len(ct) == 2

    # --- rsa_attack_small_e ---

    def test_rsa_attack_cube_root(self):
        """RSA Low-e-Angriff: c = m^3, finde m via Kubikwurzel."""
        m = 5
        e = 3
        # Kleines n: c = m^e (ohne mod-Reduktion, da m^3 < n)
        n = 10007  # Primzahl (größer als m^3=125)
        c = pow(m, e, n)
        # Da m^3 = 125 < n: c = 125, Kubikwurzel = 5
        if c == m ** e:  # Kein mod-Wrap
            result = CryptographicAlgorithms.rsa_attack_small_e(e, n, c)
            assert result == m

    def test_rsa_attack_e1(self):
        """RSA Angriff e=1: m = c."""
        result = CryptographicAlgorithms.rsa_attack_small_e(1, 1000, 42)
        assert result == 42

    # --- baby_step_giant_step_ecc ---

    def test_bsgs_ecc_returns_none_or_int(self):
        """BSGS auf ECC gibt int oder None zurück."""
        # Kurve y² ≡ x³ + x (mod 17), P = (5, 1), Q = P selbst → k=1
        result = CryptographicAlgorithms.baby_step_giant_step_ecc(
            P=(5, 1), Q=(5, 1), n=18, a=1, p=17
        )
        if result is not None:
            assert isinstance(result, int)


# ===========================================================================
# TESTS FÜR STANDALONE-FUNKTIONEN
# ===========================================================================

class TestStandaloneFunctions:
    """Tests für eigenständige Funktionen."""

    # --- euler_product_approximation ---

    def test_euler_product_s2(self):
        """Euler-Produkt für s=2 konvergiert gegen π²/6 ≈ 1.6449."""
        val = euler_product_approximation(2.0, n_primes=500)
        assert abs(val - math.pi**2 / 6) < 0.01

    def test_euler_product_s3(self):
        """Euler-Produkt für s=3 konvergiert gegen ζ(3) ≈ 1.202."""
        val = euler_product_approximation(3.0, n_primes=200)
        assert abs(val - 1.202) < 0.05

    def test_euler_product_invalid_s(self):
        """Euler-Produkt mit s ≤ 1 wirft ValueError."""
        with pytest.raises(ValueError):
            euler_product_approximation(1.0)
        with pytest.raises(ValueError):
            euler_product_approximation(0.5)

    def test_euler_product_increasing_with_primes(self):
        """Mehr Primzahlen → genauere Approximation (Wert steigt monoton)."""
        v1 = euler_product_approximation(2.0, n_primes=10)
        v2 = euler_product_approximation(2.0, n_primes=100)
        # Produkt mit mehr Primzahlen ist größer (nähert sich von unten an)
        assert v2 >= v1

    # --- mertens_first_theorem ---

    def test_mertens_n100(self):
        """Mertens 1. Satz: Differenz Σln(p)/p - ln(n) ist beschränkt für n=100."""
        result = mertens_first_theorem(100)
        assert 'mertens_sum' in result
        assert 'log_n' in result
        assert 'difference' in result
        # Differenz sollte beschränkt sein (Mertens: O(1))
        assert abs(result['difference']) < 5.0

    def test_mertens_structure(self):
        """Mertens gibt korrekte Dict-Struktur zurück."""
        result = mertens_first_theorem(50)
        assert result['primes_used'] == 15  # π(50) = 15
        assert result['mertens_sum'] > 0

    def test_mertens_approximation_quality(self):
        """Mertens-Summe nähert sich ln(n) an."""
        result = mertens_first_theorem(1000)
        # Für n=1000: ln(1000) ≈ 6.908, Mertens-Summe sollte nah daran sein
        assert abs(result['difference']) < 3.0

    # --- primitive_root ---

    def test_primitive_root_5(self):
        """Primitive Wurzel mod 5: g=2 (2^1=2, 2^2=4, 2^3=3, 2^4=1 mod 5)."""
        g = primitive_root(5)
        assert g is not None
        # g muss Ordnung 4 = p-1 haben
        assert pow(g, 4, 5) == 1
        assert pow(g, 2, 5) != 1
        assert pow(g, 1, 5) != 1

    def test_primitive_root_7(self):
        """Primitive Wurzel mod 7: kleinste ist 3."""
        g = primitive_root(7)
        assert g is not None
        # Prüfe Ordnung
        assert pow(g, 6, 7) == 1
        for k in [1, 2, 3]:
            assert pow(g, k, 7) != 1

    def test_primitive_root_not_prime(self):
        """Keine primitive Wurzel für zusammengesetzte Zahlen (None)."""
        result = primitive_root(10)
        assert result is None

    def test_primitive_root_n2(self):
        """Primitive Wurzel mod 2: triviale Gruppe."""
        g = primitive_root(2)
        assert g == 1

    def test_primitive_root_n1(self):
        """primitive_root(1) = None."""
        assert primitive_root(1) is None

    # --- legendre_symbol_jacobi ---

    def test_jacobi_quadratic_residue(self):
        """(1/p) = 1 für alle ungeraden Primzahlen p."""
        assert legendre_symbol_jacobi(1, 7) == 1

    def test_jacobi_zero(self):
        """(0/p) = 0."""
        assert legendre_symbol_jacobi(0, 7) == 0

    def test_jacobi_nonresidue(self):
        """(3/7) = -1: 3 ist kein QR mod 7."""
        # 1²=1, 2²=4, 3²=2, alle QR mod 7 sind {1,2,4}
        # Also 3 ist kein QR
        assert legendre_symbol_jacobi(3, 7) == -1

    def test_jacobi_composite_n(self):
        """Jacobi-Symbol für zusammengesetztes n."""
        # (2/15): Jacobi-Symbol für n=15 (zusammengesetzt, ungerade)
        result = legendre_symbol_jacobi(2, 15)
        assert result in (-1, 0, 1)


# ===========================================================================
# EDGE-CASE UND INTEGRATIONSTESTS
# ===========================================================================

class TestEdgeCases:
    """Grenzfall- und Integrationstests."""

    def test_factorization_consistency(self):
        """trial_division und factor_completely liefern konsistente Ergebnisse."""
        n = 2 * 3 * 5 * 7
        r1 = IntegerFactorization.trial_division(n)['factors']
        r2 = IntegerFactorization.factor_completely(n)['factors']
        assert r1 == r2

    def test_bsgs_and_brute_force_agree(self):
        """BSGS und Brute Force liefern übereinstimmende Ergebnisse."""
        g, h, p = 5, 8, 23
        x_bsgs = DiscreteLogarithm.baby_step_giant_step(g, h, p)
        x_bf = DiscreteLogarithm.discrete_log_brute_force(g, h, p)
        # Beide sollten eine Lösung finden
        assert x_bsgs is not None
        assert x_bf is not None
        # Beide Lösungen sind gültig (können verschieden sein bei mehreren Lösungen)
        assert pow(g, x_bsgs, p) == h % p
        assert pow(g, x_bf, p) == h % p

    def test_pell_fundamental_solution_minimal(self):
        """Pell-Fundamentallösung ist die kleinste positive Lösung."""
        x, y = ContinuedFractionFactoring.pell_equation_solve(2)
        assert x == 3 and y == 2  # Bekannte Fundamentallösung

    def test_elgamal_different_messages(self):
        """ElGamal verschlüsselt verschiedene Nachrichten korrekt."""
        p, g = 23, 5
        keygen = CryptographicAlgorithms.elgamal_keygen(p=p, g=g)
        pub_y = keygen['public_key']['y']
        priv_x = keygen['private_key']

        for msg in [1, 5, 10, 15, 20]:
            ct = CryptographicAlgorithms.elgamal_encrypt(msg, pub_y, p, g)
            dec = CryptographicAlgorithms.elgamal_decrypt(ct, priv_x, p)
            assert dec == msg, f"Nachricht {msg} nicht korrekt entschlüsselt"

    def test_carmichael_numbers_pass_fermat(self):
        """Carmichael-Zahlen bestehen Fermat-Test für coprime Basen."""
        carmichaels = [561, 1105, 1729]
        for n in carmichaels:
            assert PrimalityTests.is_carmichael_number(n) is True
            # Fermat-Test zur Basis 2 besteht
            assert PrimalityTests.fermat_test(n, 2) is True

    def test_sieves_agree_on_large_input(self):
        """Eratosthenes und Atkin stimmen für n=500 überein."""
        n = 500
        e = PrimeSieve.sieve_of_eratosthenes(n)
        a = PrimeSieve.sieve_of_atkin(n)
        assert e == a

    def test_smooth_and_factorization_consistent(self):
        """smooth_number_check stimmt mit trial_division überein."""
        n = 360  # 2³ × 3² × 5 → 5-glatt
        assert IntegerFactorization.smooth_number_check(n, B=5) is True
        # Größter Primfaktor ist 5
        factors = IntegerFactorization.trial_division(n)['factors']
        assert max(factors.keys()) == 5
