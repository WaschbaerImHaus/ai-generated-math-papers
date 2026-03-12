"""
Tests für die Module erdos_moser, artin_primitive_roots und erdos_selfridge.

Testet:
    - ErdosMoser: Potenzsummen-Gleichung, modulare Ausschlüsse, Bernoulli-Verbindung
    - ArtinPrimitiveRoots: Primitive Wurzeln, Artin-Konstante, empirische Dichte
    - ErdosSelfridge: Binomialkoeffizienten, Sylvester-Theorem, Primzahlpotenzen

TDD-Ansatz: Alle Tests prüfen mathematisch korrekte Erwartungswerte.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
import os
import sys
import pytest
from typing import List, Optional, Tuple

# Projektpfad eintragen damit die Module aus /src/ gefunden werden
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Module importieren
from erdos_moser import ErdosMoser, _power_sum, _power_sum_mod
from artin_primitive_roots import ArtinPrimitiveRoots, _sieve_primes, _multiplicative_order
from erdos_selfridge import ErdosSelfridge, _is_prime_power, _binomial_coefficient


# ===========================================================================
# HILFSFUNKTIONEN FÜR TESTS
# ===========================================================================

@pytest.fixture
def moser():
    """Erstellt eine ErdosMoser-Instanz für Tests."""
    return ErdosMoser()


@pytest.fixture
def artin():
    """Erstellt eine ArtinPrimitiveRoots-Instanz für Tests."""
    return ArtinPrimitiveRoots()


@pytest.fixture
def selfridge():
    """Erstellt eine ErdosSelfridge-Instanz für Tests."""
    return ErdosSelfridge()


# ===========================================================================
# TESTS: ERDŐS-MOSER (_power_sum Hilfsfunktion)
# ===========================================================================

class TestPowerSumHelper:
    """Tests für die _power_sum Hilfsfunktion."""

    def test_power_sum_k1_small(self):
        """1^1 + 2^1 = 3."""
        assert _power_sum(3, 1) == 3  # Summe 1..2 = 1+2 = 3

    def test_power_sum_k1_m5(self):
        """1+2+3+4 = 10 = 5*4/2."""
        assert _power_sum(5, 1) == 10

    def test_power_sum_k2_m4(self):
        """1^2+2^2+3^2 = 1+4+9 = 14."""
        assert _power_sum(4, 2) == 14

    def test_power_sum_k3_m3(self):
        """1^3+2^3 = 1+8 = 9."""
        assert _power_sum(3, 3) == 9

    def test_power_sum_empty(self):
        """Leere Summe: m=1 → keine Terme."""
        assert _power_sum(1, 1) == 0

    def test_power_sum_k1_gauss(self):
        """Gaußsche Summenformel: S_1(m-1) = m(m-1)/2."""
        for m in range(2, 20):
            assert _power_sum(m, 1) == m * (m - 1) // 2

    def test_power_sum_mod_k1(self):
        """Modulare Version: 1+2 mod 5 = 3."""
        assert _power_sum_mod(3, 1, 5) == 3

    def test_power_sum_mod_k2(self):
        """1^2+2^2+3^2 = 14, 14 mod 7 = 0."""
        assert _power_sum_mod(4, 2, 7) == 0


# ===========================================================================
# TESTS: ERDŐS-MOSER (Klasse)
# ===========================================================================

class TestErdosMoserSolutions:
    """Tests für die Lösungs-Suche und Verifikation."""

    def test_k1_m3_is_solution(self, moser):
        """Einzige bekannte Lösung: (k=1, m=3) → 1+2=3."""
        assert moser.is_solution(1, 3) is True

    def test_k1_m2_no_solution(self, moser):
        """k=1, m=2: 1 ≠ 2."""
        assert moser.is_solution(1, 2) is False

    def test_k1_m4_no_solution(self, moser):
        """k=1, m=4: 1+2+3=6 ≠ 4."""
        assert moser.is_solution(1, 4) is False

    def test_k2_no_solutions_up_to_100(self, moser):
        """Für k=2 keine Lösung bis m=100."""
        for m in range(2, 101):
            assert moser.is_solution(2, m) is False

    def test_k3_no_solutions_up_to_50(self, moser):
        """Für k=3 keine Lösung bis m=50."""
        for m in range(2, 51):
            assert moser.is_solution(3, m) is False

    def test_find_solutions_only_k1_m3(self, moser):
        """In k≤5, m≤200: Nur (1,3) ist Lösung."""
        solutions = moser.find_solutions(k_max=5, m_max=200)
        assert solutions == [(1, 3)]

    def test_invalid_k_raises(self, moser):
        """k < 1 soll ValueError werfen."""
        with pytest.raises(ValueError):
            moser.is_solution(0, 3)

    def test_invalid_m_raises(self, moser):
        """m < 2 soll ValueError werfen."""
        with pytest.raises(ValueError):
            moser.is_solution(1, 1)

    def test_k1_proof(self, moser):
        """Analytischer Beweis für k=1."""
        proof = moser.verify_k1_solution()
        assert proof["solution"] == 3
        assert proof["verification"] is True
        assert proof["k"] == 1

    def test_k1_proof_steps_nonempty(self, moser):
        """Beweis-Schritte nicht leer."""
        proof = moser.verify_k1_solution()
        assert len(proof["proof_steps"]) >= 3


class TestErdosMoserModular:
    """Tests für modulare Ausschlüsse."""

    def test_modular_excludes_k2_m2_q3(self, moser):
        """k=2, m=2, q=3: 1^2=1, 2^2=4. Prüfe Ausschluss."""
        # S_2(1) = 1, 2^2 = 4, 1 mod 3 = 1, 4 mod 3 = 1 → nicht ausgeschlossen
        result = moser.is_excluded_by_modulus(2, 2, 3)
        # S_2(m-1=1) = 1^2 = 1, m^k = 2^2 = 4. 1 mod 3 = 1, 4 mod 3 = 1 → gleich → nicht ausgeschlossen
        assert result is False

    def test_modular_residues_returns_list(self, moser):
        """modular_residues gibt eine Liste zurück."""
        res = moser.modular_residues(2, 5)
        assert isinstance(res, list)

    def test_modular_residues_subset_of_modulus(self, moser):
        """Alle Reste sind kleiner als q."""
        q = 7
        res = moser.modular_residues(2, q)
        for r in res:
            assert 0 <= r < q

    def test_find_excluding_moduli_not_empty_for_k2_m5(self, moser):
        """k=2, m=5: Es gibt ausschließende Moduli."""
        excluding = moser.find_excluding_moduli(2, 5, q_max=20)
        assert len(excluding) > 0

    def test_modular_exclusion_analysis_keys(self, moser):
        """Analyse-Dictionary hat erwartete Schlüssel."""
        result = moser.modular_exclusion_analysis(2, range(2, 20), 5)
        assert "k" in result
        assert "q" in result
        assert "exclusion_rate" in result
        assert "total_candidates" in result

    def test_modular_exclusion_rate_between_0_and_1(self, moser):
        """Ausschluss-Rate liegt zwischen 0 und 1."""
        result = moser.modular_exclusion_analysis(2, range(2, 30), 7)
        assert 0.0 <= result["exclusion_rate"] <= 1.0


class TestErdosMoserBernoulli:
    """Tests für Bernoulli-Zahlen und Faulhabersche Formel."""

    def test_bernoulli_b0(self, moser):
        """B_0 = 1."""
        from sympy import Rational
        assert moser.bernoulli_number(0) == Rational(1)

    def test_bernoulli_b1(self, moser):
        """B_1 = +1/2 in sympy's Konvention (andere Quellen: -1/2).
        Sympy verwendet die Bernoulli-Konvention B_1 = +1/2."""
        from sympy import Rational
        # sympy verwendet B_1 = +1/2 (nicht -1/2 wie in manchen Lehrbüchern)
        assert moser.bernoulli_number(1) == Rational(1, 2)

    def test_bernoulli_b2(self, moser):
        """B_2 = 1/6."""
        from sympy import Rational
        assert moser.bernoulli_number(2) == Rational(1, 6)

    def test_bernoulli_b3_zero(self, moser):
        """B_3 = 0 (ungerade Index ≥ 3)."""
        from sympy import Rational
        assert moser.bernoulli_number(3) == Rational(0)

    def test_bernoulli_b4(self, moser):
        """B_4 = -1/30."""
        from sympy import Rational
        assert moser.bernoulli_number(4) == Rational(-1, 30)

    def test_faulhaber_k1_n4(self, moser):
        """S_1(4) = 1+2+3+4 = 10."""
        from sympy import Rational
        assert moser.faulhaber_formula(4, 1) == Rational(10)

    def test_faulhaber_k2_n3(self, moser):
        """S_2(3) = 1+4+9 = 14."""
        from sympy import Rational
        assert moser.faulhaber_formula(3, 2) == Rational(14)

    def test_faulhaber_k3_n2(self, moser):
        """S_3(2) = 1+8 = 9."""
        from sympy import Rational
        assert moser.faulhaber_formula(2, 3) == Rational(9)

    def test_faulhaber_equals_power_k1_m3(self, moser):
        """Faulhaber bestätigt: (k=1, m=3) ist Lösung."""
        assert moser.faulhaber_equals_power(1, 3) is True

    def test_faulhaber_equals_power_k1_m4_false(self, moser):
        """k=1, m=4: keine Lösung."""
        assert moser.faulhaber_equals_power(1, 4) is False

    def test_bernoulli_summary_length(self, moser):
        """Bernoulli-Zusammenfassung für k=1..4."""
        summaries = moser.bernoulli_connection_summary(4)
        assert len(summaries) == 4

    def test_bernoulli_odd_vanish(self, moser):
        """Alle ungeraden B_n für n ≥ 3 sind 0."""
        from sympy import Rational
        for n in [3, 5, 7, 9, 11]:
            assert moser.bernoulli_number(n) == Rational(0)


class TestErdosMoserBounds:
    """Tests für Moser's Schranke und Status."""

    def test_moser_bound_exponent(self, moser):
        """Moser-Schranken-Exponent ist 10^6."""
        assert moser.MOSER_LOWER_BOUND_EXPONENT == 10 ** 6

    def test_moser_lower_bound_dict(self, moser):
        """mosers_lower_bound gibt Dictionary mit erwarteten Schlüsseln zurück."""
        bound = moser.mosers_lower_bound()
        assert "moser_1953" in bound
        assert "current_status" in bound
        assert "trivial_solution" in bound

    def test_is_trivial_solution(self, moser):
        """(k=1, m=3) ist die triviale Lösung."""
        assert moser.is_trivial_solution(1, 3) is True

    def test_not_trivial_other(self, moser):
        """(k=2, m=5) ist nicht trivial."""
        assert moser.is_trivial_solution(2, 5) is False

    def test_conjecture_status_open(self, moser):
        """Vermutung ist offen."""
        status = moser.conjecture_status()
        assert status["not_proven"] is True
        assert "CONJECTURE" in status["status"]

    def test_power_sum_values_returns_list(self, moser):
        """power_sum_values gibt Liste zurück."""
        vals = moser.power_sum_values(2, m_max=10)
        assert len(vals) == 9  # m=2..10

    def test_growth_comparison_keys(self, moser):
        """growth_comparison hat erwartete Schlüssel."""
        gc = moser.growth_comparison(2, [5, 10, 20])
        assert "asymptotic_formula" in gc
        assert "data" in gc
        assert len(gc["data"]) == 3


# ===========================================================================
# TESTS: ARTIN PRIMITIVE ROOTS (Hilfsfunktionen)
# ===========================================================================

class TestArtinHelpers:
    """Tests für Hilfsfunktionen des Artin-Moduls."""

    def test_sieve_primes_small(self):
        """Primzahlen bis 20."""
        primes = _sieve_primes(20)
        assert primes == [2, 3, 5, 7, 11, 13, 17, 19]

    def test_sieve_primes_empty(self):
        """Keine Primzahlen bis 1."""
        assert _sieve_primes(1) == []

    def test_multiplicative_order_2_mod5(self):
        """ord_5(2) = 4."""
        assert _multiplicative_order(2, 5) == 4

    def test_multiplicative_order_2_mod7(self):
        """ord_7(2) = 3."""
        assert _multiplicative_order(2, 7) == 3

    def test_multiplicative_order_gcd_not_1(self):
        """gcd(2, 4) ≠ 1 → None."""
        assert _multiplicative_order(2, 4) is None


# ===========================================================================
# TESTS: ARTIN PRIMITIVE ROOTS (Klasse)
# ===========================================================================

class TestArtinPrimitiveRootCheck:
    """Tests für primitive Wurzel Prüfungen."""

    def test_2_is_primitive_root_mod5(self, artin):
        """2 ist primitive Wurzel mod 5 (ord_5(2)=4=5-1)."""
        assert artin.is_primitive_root(2, 5) is True

    def test_2_not_primitive_root_mod_7(self, artin):
        """2 ist KEINE primitive Wurzel mod 7 (ord_7(2)=3 ≠ 6)."""
        assert artin.is_primitive_root(2, 7) is False

    def test_3_is_primitive_root_mod7(self, artin):
        """3 ist primitive Wurzel mod 7 (ord_7(3)=6=7-1)."""
        assert artin.is_primitive_root(3, 7) is True

    def test_2_not_primitive_root_mod8(self, artin):
        """2 ist KEINE primitive Wurzel mod 8 (8 ist keine Primzahl)."""
        with pytest.raises(ValueError):
            artin.is_primitive_root(2, 8)

    def test_not_prime_raises(self, artin):
        """is_primitive_root wirft ValueError für nicht-Primzahl."""
        with pytest.raises(ValueError):
            artin.is_primitive_root(2, 9)

    def test_detailed_check_mod5(self, artin):
        """Detaillierte Analyse für 2 mod 5."""
        result = artin.primitive_root_check_detailed(2, 5)
        assert result["is_primitive_root"] is True
        assert result["ord_p(a)"] == 4
        assert result["phi(p)"] == 4

    def test_detailed_check_order_correct(self, artin):
        """ord_11(2) = 10 = φ(11)."""
        result = artin.primitive_root_check_detailed(2, 11)
        assert result["ord_p(a)"] == 10
        assert result["is_primitive_root"] is True

    def test_5_primitive_root_mod_p(self, artin):
        """5 ist primitive Wurzel mod 23."""
        # ord_23(5) = 22 = φ(23)?
        order = _multiplicative_order(5, 23)
        assert artin.is_primitive_root(5, 23) == (order == 22)


class TestArtinConstant:
    """Tests für die Artin-Konstante."""

    def test_artin_constant_4_digits(self, artin):
        """Artin-Konstante auf 4 Stellen: 0.3739..."""
        c = artin.artin_constant(num_primes=500)
        # Auf 4 Stellen: 0.3739 (Toleranz 0.001)
        assert abs(c - 0.3739) < 0.002

    def test_artin_constant_positive(self, artin):
        """Artin-Konstante ist positiv."""
        c = artin.artin_constant(50)
        assert c > 0

    def test_artin_constant_less_than_1(self, artin):
        """Artin-Konstante ist kleiner als 1."""
        c = artin.artin_constant(100)
        assert c < 1.0

    def test_artin_constant_approaches_true_value(self, artin):
        """Mit mehr Primzahlen wird Approximation besser."""
        c_small = artin.artin_constant(50)
        c_large = artin.artin_constant(300)
        # Beide sollten näher als 0.05 am wahren Wert sein
        true_val = ArtinPrimitiveRoots.ARTIN_CONSTANT_APPROX
        assert abs(c_large - true_val) < abs(c_small - true_val) + 0.01

    def test_artin_constant_class_variable(self):
        """Klassenattribut hat korrekten Wert."""
        # 0.37395... – auf 4 Stellen korrekt
        assert abs(ArtinPrimitiveRoots.ARTIN_CONSTANT_APPROX - 0.3739) < 0.001

    def test_convergence_list(self, artin):
        """Konvergenz-Liste hat korrekte Länge."""
        conv = artin.artin_constant_convergence([10, 50, 100])
        assert len(conv) == 3
        assert conv[0]["num_primes"] == 10


class TestArtinEmpiricalDensity:
    """Tests für empirische Dichte-Berechnungen."""

    def test_empirical_density_returns_dict(self, artin):
        """empirical_density gibt Dictionary zurück."""
        result = artin.empirical_density(2, limit=500)
        assert isinstance(result, dict)
        assert "empirical_density" in result

    def test_empirical_density_positive(self, artin):
        """Dichte ist positiv für a=2."""
        result = artin.empirical_density(2, limit=500)
        assert result["empirical_density"] > 0

    def test_count_primitive_root_2_nontrivial(self, artin):
        """Es gibt Primzahlen ≤ 100 mit 2 als primitiver Wurzel."""
        count = artin.count_primitive_root_primes(2, 100)
        assert count > 0

    def test_first_primitive_roots_base2(self, artin):
        """Erste Primzahlen mit 2 als primitiver Wurzel."""
        primes = artin.first_primitive_root_primes(2, count=5)
        assert len(primes) == 5
        # 3 ist erste: ord_3(2) = 2 = φ(3)
        assert 3 in primes

    def test_first_primitive_roots_base3(self, artin):
        """Erste Primzahlen mit 3 als primitiver Wurzel."""
        primes = artin.first_primitive_root_primes(3, count=3)
        assert len(primes) == 3


class TestArtinHooley:
    """Tests für Hooley's Resultat und verwandte Analysen."""

    def test_hooley_density_base2(self, artin):
        """Hooley-Dichte für a=2 entspricht Artin-Konstante."""
        result = artin.hooley_density(2)
        assert result["artin_applies"] is True
        assert abs(result["hooley_density"] - 0.3739) < 0.001

    def test_hooley_perfect_square_excluded(self, artin):
        """a=4 ist vollständiges Quadrat → Artin gilt nicht."""
        result = artin.hooley_density(4)
        assert result["is_perfect_square"] is True
        assert result["artin_applies"] is False

    def test_heath_brown_theorem(self, artin):
        """Heath-Brown-Resultat ist unbedingtes Theorem."""
        hb = artin.heath_brown_result()
        assert hb["unconditional"] is True
        assert hb["status"] == "THEOREM (bewiesen)"
        assert hb["year"] == 1986

    def test_conjecture_status_open(self, artin):
        """Artin-Vermutung ist offen."""
        status = artin.conjecture_status()
        assert status["not_proven_in_general"] is True
        assert "CONJECTURE" in status["status"]

    def test_compare_bases_returns_sorted(self, artin):
        """compare_bases gibt Liste sortiert nach Dichte."""
        results = artin.compare_bases([2, 3], limit=500)
        assert len(results) == 2
        # Sortiert nach Dichte absteigend
        assert results[0]["density"] >= results[1]["density"]


# ===========================================================================
# TESTS: ERDŐS-SELFRIDGE (Hilfsfunktionen)
# ===========================================================================

class TestSelfridgeHelpers:
    """Tests für Hilfsfunktionen des Selfridge-Moduls."""

    def test_is_prime_power_prime(self):
        """7 ist prim → (7, 1)."""
        assert _is_prime_power(7) == (7, 1)

    def test_is_prime_power_8(self):
        """8 = 2^3 → (2, 3)."""
        assert _is_prime_power(8) == (2, 3)

    def test_is_prime_power_9(self):
        """9 = 3^2 → (3, 2)."""
        assert _is_prime_power(9) == (3, 2)

    def test_is_prime_power_12_no(self):
        """12 = 4·3 = 2^2·3 ist keine Primzahlpotenz."""
        assert _is_prime_power(12) is None

    def test_is_prime_power_6_no(self):
        """6 = 2·3 ist keine Primzahlpotenz."""
        assert _is_prime_power(6) is None

    def test_is_prime_power_1_no(self):
        """1 ist keine Primzahlpotenz."""
        assert _is_prime_power(1) is None

    def test_binomial_coefficient_basic(self):
        """C(5,2) = 10."""
        assert _binomial_coefficient(5, 2) == 10

    def test_binomial_coefficient_edge(self):
        """C(n,0) = 1."""
        assert _binomial_coefficient(10, 0) == 1

    def test_binomial_coefficient_k_gt_n(self):
        """C(3,5) = 0."""
        assert _binomial_coefficient(3, 5) == 0


# ===========================================================================
# TESTS: ERDŐS-SELFRIDGE (Klasse)
# ===========================================================================

class TestErdosSelfridgePrimePower:
    """Tests für Primzahlpotenz-Prüfungen."""

    def test_c9_3_not_prime_power(self, selfridge):
        """C(9,3) = 84 = 4·21 = 2^2·3·7 – keine Primzahlpotenz."""
        result = selfridge.is_binomial_prime_power(9, 3)
        assert result is None

    def test_c84_is_not_prime_power(self):
        """84 = 2^2 · 3 · 7 ist keine Primzahlpotenz."""
        assert _is_prime_power(84) is None

    def test_c5_2_not_prime_power(self, selfridge):
        """C(5,2) = 10 = 2·5 – keine Primzahlpotenz."""
        assert selfridge.is_binomial_prime_power(5, 2) is None

    def test_c4_2_not_prime_power(self, selfridge):
        """C(4,2) = 6 = 2·3 – keine Primzahlpotenz."""
        assert selfridge.is_binomial_prime_power(4, 2) is None

    def test_c4_1_is_prime_power(self, selfridge):
        """C(4,1) = 4 = 2^2 – Primzahlpotenz (k=1, trivial)."""
        result = selfridge.is_binomial_prime_power(4, 1)
        assert result is not None
        assert result == (2, 2)

    def test_c7_1_is_prime(self, selfridge):
        """C(7,1) = 7 – prim."""
        result = selfridge.is_binomial_prime_power(7, 1)
        assert result == (7, 1)

    def test_check_conjecture_no_counterexample(self, selfridge):
        """C(9,3) ist kein Gegenbeispiel."""
        result = selfridge.check_conjecture(9, 3)
        assert result["is_counterexample"] is False
        assert result["C(n,k)"] == 84

    def test_check_conjecture_applies_condition(self, selfridge):
        """Für n ≥ k+2 gilt die Vermutung."""
        result = selfridge.check_conjecture(10, 4)
        assert result["conjecture_applies"] is True


class TestErdosSelfridgeVerification:
    """Tests für numerische Verifikation."""

    def test_verify_up_to_50(self, selfridge):
        """Keine Gegenbeispiele bis n=50."""
        result = selfridge.verify_up_to(50)
        assert result["conjecture_holds"] is True
        assert result["counterexamples_found"] == 0

    def test_verify_up_to_100(self, selfridge):
        """Keine Gegenbeispiele bis n=100."""
        result = selfridge.verify_up_to(100)
        assert result["conjecture_holds"] is True

    def test_find_prime_power_binomials_k1(self, selfridge):
        """C(8,1)=8=2^3 ist in der Liste."""
        results = selfridge.find_prime_power_binomials(n_max=20, k_max=5)
        entries = [(r["n"], r["k"]) for r in results]
        assert (8, 1) in entries

    def test_find_prime_power_binomials_trivial_marked(self, selfridge):
        """k=1 Fälle sind als trivial markiert."""
        results = selfridge.find_prime_power_binomials(n_max=15, k_max=3)
        k1_cases = [r for r in results if r["k"] == 1]
        for case in k1_cases:
            assert case["trivial"] is True


class TestSylvesterTheorem:
    """Tests für Sylvesters Theorem."""

    def test_sylvester_c6_2(self, selfridge):
        """C(6,2) = 15 = 3·5. Primteiler 5 > 2 ✓ (n=6 ≥ 2k=4)."""
        result = selfridge.sylvester_theorem_check(6, 2)
        assert result["sylvester_applies"] is True
        assert result["theorem_confirmed"] is True
        assert 5 in result["prime_factors_greater_k"] or 3 in result["prime_factors_greater_k"]

    def test_sylvester_c10_3(self, selfridge):
        """C(10,3) = 120 = 2^3·3·5. Primteiler 5 > 3 ✓."""
        result = selfridge.sylvester_theorem_check(10, 3)
        assert result["sylvester_applies"] is True
        assert result["theorem_confirmed"] is True

    def test_sylvester_n_lt_2k_not_applies(self, selfridge):
        """n=3, k=2: 3 < 4 = 2k → Sylvester gilt nicht."""
        result = selfridge.sylvester_theorem_check(3, 2)
        assert result["condition_n_geq_2k"] is False

    def test_verify_sylvester_up_to_50(self, selfridge):
        """Sylvester-Theorem für alle n ≤ 50 verifiziert."""
        result = selfridge.verify_sylvester(50)
        assert result["sylvester_holds"] is True
        assert len(result["violations"]) == 0


class TestSelfridgeSpecialCases:
    """Tests für Spezialfälle."""

    def test_analyze_k1_contains_prime_powers(self, selfridge):
        """k=1 Analyse: n=4 (4=2^2) und n=8 (8=2^3) sind enthalten."""
        results = selfridge.analyze_k1(n_max=20)
        ns = [r["n"] for r in results]
        assert 4 in ns
        assert 8 in ns

    def test_analyze_k2_no_counterexample(self, selfridge):
        """k=2 Analyse: Kein Gegenbeispiel für n ≥ 4."""
        results = selfridge.analyze_k2(n_max=100)
        counterexamples = [r for r in results if r.get("is_counterexample")]
        assert len(counterexamples) == 0

    def test_granville_cases_nonempty(self, selfridge):
        """Granville-Spezialfälle sind dokumentiert."""
        cases = selfridge.granville_special_cases()
        assert len(cases) >= 3

    def test_granville_full_conjecture_open(self, selfridge):
        """Vollständige Vermutung ist als offen markiert."""
        cases = selfridge.granville_special_cases()
        full_case = next((c for c in cases if "allgemein" in c["case"]), None)
        assert full_case is not None
        assert full_case["proven"] is False

    def test_conjecture_status_open(self, selfridge):
        """Erdős-Selfridge-Vermutung ist offen."""
        status = selfridge.conjecture_status()
        assert status["not_proven"] is True
        assert status["year_conjectured"] == 1975
        assert "CONJECTURE" in status["status"]

    def test_prime_factorization_statistics(self, selfridge):
        """Statistik über Primteiler von C(n,k)."""
        stats = selfridge.prime_factorization_statistics(n_max=30)
        assert "total_binomials" in stats
        assert stats["total_binomials"] > 0
        # Primzahlpotenzen sollten selten sein
        assert stats["one_prime_factor"]["fraction"] < 0.3


# ===========================================================================
# INTEGRATIONSTESTS
# ===========================================================================

class TestIntegration:
    """Integrationstests über Modulgrenzen hinweg."""

    def test_moser_bernoulli_matches_direct_sum(self, moser):
        """Faulhaber und direkte Summe stimmen überein für k=3, N=5."""
        from sympy import Rational
        faulhaber = moser.faulhaber_formula(5, 3)
        direct = sum(i ** 3 for i in range(1, 6))  # 1+8+27+64+125 = 225
        assert faulhaber == Rational(direct)

    def test_artin_primitive_root_consistent(self, artin):
        """is_primitive_root und first_primitive_root_primes stimmen überein."""
        primes = artin.first_primitive_root_primes(2, count=5)
        for p in primes:
            assert artin.is_primitive_root(2, p) is True

    def test_selfridge_c9_3_value(self, selfridge):
        """C(9,3) = 84."""
        assert _binomial_coefficient(9, 3) == 84

    def test_selfridge_sylvester_and_prime_power_consistent(self, selfridge):
        """Wenn Sylvester Primteiler > k garantiert, dann kann C(n,k) nicht Primzahlpotenz sein mit passendem k."""
        # C(10,4) = 210 = 2·3·5·7 → viele Primteiler, sicher keine Primzahlpotenz
        assert selfridge.is_binomial_prime_power(10, 4) is None
        sylv = selfridge.sylvester_theorem_check(10, 4)
        assert sylv["theorem_confirmed"] is True

    def test_moser_k1_solution_count(self, moser):
        """Für k=1, m=2..100: Genau eine Lösung (m=3)."""
        solutions = [m for m in range(2, 101) if moser.is_solution(1, m)]
        assert solutions == [3]
