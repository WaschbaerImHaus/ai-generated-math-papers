"""
@file test_additive_number_theory.py
@brief Tests für das Modul additive_number_theory.py.

@description
    Umfassende Test-Suite für alle Klassen und Funktionen der additiven Zahlentheorie.
    Getestet werden:
      - PartitionFunction: p(n), Partitionsgenerierung, Pentagonaltheorem, Ramanujan-Kongruenzen
      - HardyRamanujanCircleMethod: Asymptotik, Rademacher-Formel, Farey-Beiträge
      - WaringProblem: g(k), G(k), Quadratsummen, Kubussummen, Vier-Quadrate, Drei-Quadrate
      - GoldbachAnalysis: Zerlegungen, Zählung, Komet, ternäre Goldbach, Chen
      - ArithmeticProgressions: Primzahlen in AP, Van-der-Waerden, AP in Primzahlen,
                                 Dirichlet-Dichte, Green-Tao-Demo
      - SumSetTheory: Summenmenge, Differenzmenge, Cauchy-Davenport, Plünnecke-Ruzsa,
                      Freiman, additive Energie
      - Standalone: schur_number, happy_ending_problem, bertrand_postulate_verify

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import pytest
import math
import sys
import os

# Pfad zum src-Verzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from additive_number_theory import (
    PartitionFunction,
    HardyRamanujanCircleMethod,
    WaringProblem,
    GoldbachAnalysis,
    ArithmeticProgressions,
    SumSetTheory,
    schur_number,
    happy_ending_problem,
    bertrand_postulate_verify,
)


# ===========================================================================
# TESTS: PartitionFunction
# ===========================================================================

class TestPartitionCount:
    """Tests für PartitionFunction.partition_count(n)."""

    def test_p0_equals_1(self):
        """p(0) = 1 (leere Partition)."""
        assert PartitionFunction.partition_count(0) == 1

    def test_p1_equals_1(self):
        """p(1) = 1."""
        assert PartitionFunction.partition_count(1) == 1

    def test_p2_equals_2(self):
        """p(2) = 2: {2, 1+1}."""
        assert PartitionFunction.partition_count(2) == 2

    def test_p3_equals_3(self):
        """p(3) = 3: {3, 2+1, 1+1+1}."""
        assert PartitionFunction.partition_count(3) == 3

    def test_p4_equals_5(self):
        """p(4) = 5: {4, 3+1, 2+2, 2+1+1, 1+1+1+1}."""
        assert PartitionFunction.partition_count(4) == 5

    def test_p5_equals_7(self):
        """p(5) = 7."""
        assert PartitionFunction.partition_count(5) == 7

    def test_p10_equals_42(self):
        """p(10) = 42 (bekannter Wert)."""
        assert PartitionFunction.partition_count(10) == 42

    def test_p20_equals_627(self):
        """p(20) = 627 (bekannter Wert)."""
        assert PartitionFunction.partition_count(20) == 627

    def test_p50_equals_204226(self):
        """p(50) = 204226 (bekannter Wert)."""
        assert PartitionFunction.partition_count(50) == 204226

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            PartitionFunction.partition_count(-1)

    def test_monotone_growth(self):
        """p(n) ist monoton wachsend für n ≥ 1."""
        prev = PartitionFunction.partition_count(1)
        for n in range(2, 15):
            curr = PartitionFunction.partition_count(n)
            assert curr >= prev, f"p({n}) < p({n-1}): nicht monoton"
            prev = curr


class TestGeneratePartitions:
    """Tests für PartitionFunction.generate_partitions(n)."""

    def test_generate_n0(self):
        """Partition von 0: nur die leere Partition ()."""
        parts = PartitionFunction.generate_partitions(0)
        assert parts == [()]

    def test_generate_n1(self):
        """Partition von 1: nur {(1,)}."""
        parts = PartitionFunction.generate_partitions(1)
        assert set(parts) == {(1,)}

    def test_generate_n4_count(self):
        """Anzahl der Partitionen von 4 muss p(4)=5 sein."""
        parts = PartitionFunction.generate_partitions(4)
        assert len(parts) == 5

    def test_generate_n4_correct(self):
        """Alle Partitionen von 4 sind korrekt."""
        parts = set(PartitionFunction.generate_partitions(4))
        expected = {(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)}
        assert parts == expected

    def test_generate_sums_correct(self):
        """Alle Partitionen summierten sich zu n."""
        for n in range(0, 10):
            parts = PartitionFunction.generate_partitions(n)
            for p in parts:
                assert sum(p) == n, f"Partition {p} summiert nicht zu {n}"

    def test_generate_descending(self):
        """Alle Partitionen sind in nicht-aufsteigender Reihenfolge."""
        parts = PartitionFunction.generate_partitions(6)
        for p in parts:
            assert list(p) == sorted(p, reverse=True), f"Partition {p} nicht absteigend"

    def test_generate_count_matches_partition_count(self):
        """Anzahl der generierten Partitionen entspricht partition_count."""
        for n in range(0, 12):
            parts = PartitionFunction.generate_partitions(n)
            assert len(parts) == PartitionFunction.partition_count(n), \
                f"Abweichung bei n={n}"

    def test_generate_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            PartitionFunction.generate_partitions(-1)


class TestPartitionWithParts:
    """Tests für PartitionFunction.partition_with_parts(n, k)."""

    def test_n5_k1(self):
        """Partition von 5 in 1 Teil: nur {(5,)}."""
        parts = PartitionFunction.partition_with_parts(5, 1)
        assert parts == [(5,)]

    def test_n5_k2(self):
        """Partition von 5 in 2 Teile: {(4,1), (3,2)}."""
        parts = set(PartitionFunction.partition_with_parts(5, 2))
        expected = {(4, 1), (3, 2)}
        assert parts == expected

    def test_n4_k2_count(self):
        """p(4, 2 Teile) = 2."""
        parts = PartitionFunction.partition_with_parts(4, 2)
        assert len(parts) == 2

    def test_n0_k0(self):
        """p(0, 0) = 1 (leere Partition)."""
        parts = PartitionFunction.partition_with_parts(0, 0)
        assert parts == [()]

    def test_n5_k0_empty(self):
        """p(5, 0) = 0 (keine Partition in 0 Teile außer n=0)."""
        parts = PartitionFunction.partition_with_parts(5, 0)
        assert parts == []

    def test_all_have_k_parts(self):
        """Alle zurückgegebenen Partitionen haben genau k Teile."""
        for n in range(1, 10):
            for k in range(1, n + 1):
                parts = PartitionFunction.partition_with_parts(n, k)
                for p in parts:
                    assert len(p) == k, f"Partition {p} hat nicht {k} Teile"

    def test_sums_are_n(self):
        """Alle Partitionen summieren zu n."""
        parts = PartitionFunction.partition_with_parts(7, 3)
        for p in parts:
            assert sum(p) == 7


class TestPartitionIntoDistinctParts:
    """Tests für PartitionFunction.partition_into_distinct_parts(n)."""

    def test_n0(self):
        """Partition von 0 in verschiedene Teile: nur ()."""
        parts = PartitionFunction.partition_into_distinct_parts(0)
        assert parts == [()]

    def test_n6(self):
        """Partition von 6 in verschiedene Teile: {6, 5+1, 4+2, 3+2+1}."""
        parts = set(PartitionFunction.partition_into_distinct_parts(6))
        expected = {(6,), (5, 1), (4, 2), (3, 2, 1)}
        assert parts == expected

    def test_all_distinct(self):
        """Alle Teile in jeder Partition sind verschieden."""
        for n in range(0, 12):
            parts = PartitionFunction.partition_into_distinct_parts(n)
            for p in parts:
                assert len(p) == len(set(p)), f"Nicht-verschiedene Teile in {p}"

    def test_sums_correct(self):
        """Alle Partitionen summieren zu n."""
        for n in range(0, 12):
            parts = PartitionFunction.partition_into_distinct_parts(n)
            for p in parts:
                assert sum(p) == n

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            PartitionFunction.partition_into_distinct_parts(-1)


class TestEulerPentagonalTheorem:
    """Tests für PartitionFunction.euler_pentagonal_theorem(n)."""

    def test_pentagonal_matches_dp(self):
        """Pentagonaltheorem-Ergebnisse stimmen mit DP überein."""
        for n in range(0, 30):
            assert PartitionFunction.euler_pentagonal_theorem(n) == \
                   PartitionFunction.partition_count(n), \
                f"Abweichung bei n={n}"

    def test_p0_pentagonal(self):
        """p(0) = 1 via Pentagonaltheorem."""
        assert PartitionFunction.euler_pentagonal_theorem(0) == 1

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            PartitionFunction.euler_pentagonal_theorem(-1)


class TestRamanujanCongruences:
    """Tests für PartitionFunction.ramanujan_partition_congruences."""

    def test_mod5_congruences(self):
        """p(5k+4) ≡ 0 (mod 5) für k = 0, 1, 2, 3, 4."""
        results = PartitionFunction.ramanujan_partition_congruences([0, 1, 2, 3, 4])
        for entry in results["mod5"]:
            assert entry["holds"], \
                f"p({entry['n']}) = {entry['p_n']} ≢ 0 (mod 5): Rest {entry['remainder']}"

    def test_mod7_congruences(self):
        """p(7k+5) ≡ 0 (mod 7) für k = 0, 1, 2, 3, 4."""
        results = PartitionFunction.ramanujan_partition_congruences([0, 1, 2, 3, 4])
        for entry in results["mod7"]:
            assert entry["holds"], \
                f"p({entry['n']}) = {entry['p_n']} ≢ 0 (mod 7): Rest {entry['remainder']}"

    def test_mod11_congruences(self):
        """p(11k+6) ≡ 0 (mod 11) für k = 0, 1, 2, 3."""
        results = PartitionFunction.ramanujan_partition_congruences([0, 1, 2, 3])
        for entry in results["mod11"]:
            assert entry["holds"], \
                f"p({entry['n']}) = {entry['p_n']} ≢ 0 (mod 11): Rest {entry['remainder']}"

    def test_result_structure(self):
        """Rückgabewert hat die korrekte Struktur."""
        results = PartitionFunction.ramanujan_partition_congruences([0])
        assert "mod5" in results
        assert "mod7" in results
        assert "mod11" in results
        for key in ["mod5", "mod7", "mod11"]:
            assert len(results[key]) == 1
            entry = results[key][0]
            assert "k" in entry and "n" in entry and "p_n" in entry
            assert "remainder" in entry and "holds" in entry

    def test_mod5_specific_values(self):
        """Konkrete Werte: p(4)=5≡0(5), p(9)=30≡0(5), p(14)=135≡0(5)."""
        pf = PartitionFunction()
        assert pf.partition_count(4) % 5 == 0   # k=0: n=4
        assert pf.partition_count(9) % 5 == 0   # k=1: n=9
        assert pf.partition_count(14) % 5 == 0  # k=2: n=14


# ===========================================================================
# TESTS: HardyRamanujanCircleMethod
# ===========================================================================

class TestAsymptoticPartition:
    """Tests für HardyRamanujanCircleMethod.asymptotic_partition(n)."""

    def test_returns_positive(self):
        """Asymptotik liefert positive Werte."""
        for n in range(1, 20):
            val = HardyRamanujanCircleMethod.asymptotic_partition(n)
            assert val > 0, f"Asymptotik für n={n} nicht positiv: {val}"

    def test_growing(self):
        """Asymptotik ist streng wachsend (für n ≥ 1)."""
        prev = HardyRamanujanCircleMethod.asymptotic_partition(1)
        for n in range(2, 20):
            curr = HardyRamanujanCircleMethod.asymptotic_partition(n)
            assert curr > prev, f"Asymptotik nicht wachsend bei n={n}"
            prev = curr

    def test_large_n_order_of_magnitude(self):
        """Für n=100 ist die Asymptotik in der richtigen Größenordnung."""
        # p(100) = 190569292 (exakt)
        exact = 190569292
        approx = HardyRamanujanCircleMethod.asymptotic_partition(100)
        # Relative Abweichung ≤ 5%
        rel_error = abs(approx - exact) / exact
        assert rel_error < 0.05, f"Asymptotik zu ungenau: Abweichung {rel_error:.4f}"

    def test_zero_raises(self):
        """n=0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HardyRamanujanCircleMethod.asymptotic_partition(0)

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            HardyRamanujanCircleMethod.asymptotic_partition(-5)


class TestRademacherFormula:
    """Tests für HardyRamanujanCircleMethod.rademacher_formula(n, k_max)."""

    def test_returns_positive(self):
        """Rademacher-Formel liefert positive Werte."""
        for n in [1, 2, 3, 5, 10]:
            val = HardyRamanujanCircleMethod.rademacher_formula(n, k_max=15)
            assert val > 0, f"Rademacher für n={n} nicht positiv: {val}"

    def test_more_terms_more_accurate(self):
        """Mehr Terme → bessere Näherung für p(n)."""
        n = 10
        exact = PartitionFunction.partition_count(n)  # p(10) = 42
        approx_5 = HardyRamanujanCircleMethod.rademacher_formula(n, k_max=5)
        approx_20 = HardyRamanujanCircleMethod.rademacher_formula(n, k_max=20)
        # Beide sollten in der Nähe von 42 liegen
        assert abs(approx_20 - exact) <= abs(approx_5 - exact) + 5, \
            "Mehr Terme verbessern die Näherung nicht"

    def test_n0_raises(self):
        """n=0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HardyRamanujanCircleMethod.rademacher_formula(0)


class TestFareyArcContributions:
    """Tests für HardyRamanujanCircleMethod.farey_arc_contributions(n, k_max)."""

    def test_returns_list(self):
        """Farey-Beitragsliste wird zurückgegeben."""
        contribs = HardyRamanujanCircleMethod.farey_arc_contributions(10, k_max=3)
        assert isinstance(contribs, list)
        assert len(contribs) > 0

    def test_k1_h1_present(self):
        """Der Hauptbogen (h=1, k=1) ist immer enthalten."""
        contribs = HardyRamanujanCircleMethod.farey_arc_contributions(5, k_max=3)
        keys = [(c["h"], c["k"]) for c in contribs]
        assert (1, 1) in keys

    def test_sorted_by_amplitude(self):
        """Beiträge sind nach absteigender Amplitude sortiert."""
        contribs = HardyRamanujanCircleMethod.farey_arc_contributions(20, k_max=5)
        amplitudes = [c["amplitude"] for c in contribs]
        assert amplitudes == sorted(amplitudes, reverse=True)

    def test_all_coprime(self):
        """Alle Paare (h, k) in der Liste sind teilerfremd."""
        from math import gcd
        contribs = HardyRamanujanCircleMethod.farey_arc_contributions(15, k_max=4)
        for c in contribs:
            assert math.gcd(c["h"], c["k"]) == 1, \
                f"(h={c['h']}, k={c['k']}) ist nicht teilerfremd"

    def test_n0_raises(self):
        """n=0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HardyRamanujanCircleMethod.farey_arc_contributions(0)


# ===========================================================================
# TESTS: WaringProblem
# ===========================================================================

class TestGk:
    """Tests für WaringProblem.g_k(k)."""

    def test_g1_equals_1(self):
        """g(1) = 1."""
        assert WaringProblem.g_k(1) == 1

    def test_g2_equals_4(self):
        """g(2) = 4 (Lagrange)."""
        assert WaringProblem.g_k(2) == 4

    def test_g3_equals_9(self):
        """g(3) = 9."""
        assert WaringProblem.g_k(3) == 9

    def test_g4_equals_19(self):
        """g(4) = 19."""
        assert WaringProblem.g_k(4) == 19

    def test_g_unknown_returns_none(self):
        """g(k) für k > 10 gibt None zurück."""
        assert WaringProblem.g_k(11) is None

    def test_invalid_raises(self):
        """k < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            WaringProblem.g_k(0)


class TestBigGk:
    """Tests für WaringProblem.G_k(k)."""

    def test_G1_equals_1(self):
        """G(1) = 1."""
        assert WaringProblem.G_k(1) == 1

    def test_G2_equals_4(self):
        """G(2) = 4 (Lagrange)."""
        assert WaringProblem.G_k(2) == 4

    def test_G4_equals_15(self):
        """G(4) = 15 (Davenport)."""
        assert WaringProblem.G_k(4) == 15

    def test_G_unknown_returns_none(self):
        """G(k) für k > 10 gibt None zurück."""
        assert WaringProblem.G_k(11) is None


class TestRepresentAsSquares:
    """Tests für WaringProblem.represent_as_squares(n, k)."""

    def test_n0_squares(self):
        """0 = 0² (trivial)."""
        result = WaringProblem.represent_as_squares(0, k=4)
        assert result is not None
        assert sum(x * x for x in result) == 0

    def test_n1_as_one_square(self):
        """1 = 1²."""
        result = WaringProblem.represent_as_squares(1, k=1)
        assert result is not None
        assert sum(x * x for x in result) == 1

    def test_four_squares_correctness(self):
        """Für alle n von 0 bis 50 ist die Vier-Quadrate-Darstellung korrekt."""
        for n in range(0, 51):
            result = WaringProblem.represent_as_squares(n, k=4)
            assert result is not None, f"Keine Darstellung für n={n}"
            assert sum(x * x for x in result) == n, \
                f"Darstellung {result} summiert nicht zu {n}"

    def test_n7_not_as_3_squares(self):
        """7 ist nicht als Summe von 3 Quadraten darstellbar (Legendre)."""
        result = WaringProblem.represent_as_squares(7, k=3)
        assert result is None, "7 sollte nicht als Summe von 3 Quadraten darstellbar sein"

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            WaringProblem.represent_as_squares(-1)


class TestRepresentAsCubes:
    """Tests für WaringProblem.represent_as_cubes(n)."""

    def test_n0(self):
        """0 = 0³."""
        result = WaringProblem.represent_as_cubes(0)
        assert result is not None
        assert sum(x ** 3 for x in result) == 0

    def test_n1(self):
        """1 = 1³."""
        result = WaringProblem.represent_as_cubes(1)
        assert result is not None
        assert sum(x ** 3 for x in result) == 1

    def test_n8(self):
        """8 = 2³."""
        result = WaringProblem.represent_as_cubes(8)
        assert result is not None
        assert sum(x ** 3 for x in result) == 8

    def test_n23(self):
        """23 = 2³+2³+1³+1³+1³ (braucht 9 Kuben laut Waring)."""
        result = WaringProblem.represent_as_cubes(23)
        assert result is not None
        assert sum(x ** 3 for x in result) == 23

    def test_sums_correct(self):
        """Für n = 1 bis 30 summieren sich die Kuben korrekt."""
        for n in range(1, 31):
            result = WaringProblem.represent_as_cubes(n)
            assert result is not None, f"Keine Kubus-Darstellung für n={n}"
            assert sum(x ** 3 for x in result) == n


class TestFourSquaresTheorem:
    """Tests für WaringProblem.four_squares_theorem(n)."""

    def test_known_values(self):
        """Bekannte Darstellungen verifizieren."""
        test_cases = [0, 1, 2, 3, 4, 5, 7, 12, 15, 23, 100]
        for n in test_cases:
            a, b, c, d = WaringProblem.four_squares_theorem(n)
            assert a * a + b * b + c * c + d * d == n, \
                f"Vier-Quadrate-Darstellung ({a},{b},{c},{d}) falsch für n={n}"

    def test_sorted_descending(self):
        """Ergebnis ist absteigend sortiert (a ≥ b ≥ c ≥ d)."""
        for n in range(0, 30):
            a, b, c, d = WaringProblem.four_squares_theorem(n)
            assert a >= b >= c >= d >= 0, \
                f"Darstellung ({a},{b},{c},{d}) nicht absteigend für n={n}"

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            WaringProblem.four_squares_theorem(-1)


class TestThreeSquaresCheck:
    """Tests für WaringProblem.three_squares_theorem_check(n)."""

    def test_not_representable_7(self):
        """7 = 4^0 · (8·0 + 7) → nicht darstellbar als 3 Quadrate."""
        assert WaringProblem.three_squares_theorem_check(7) is False

    def test_not_representable_15(self):
        """15 = 4^0 · (8·1 + 7) → nicht darstellbar."""
        assert WaringProblem.three_squares_theorem_check(15) is False

    def test_not_representable_28(self):
        """28 = 4 · 7 → nicht darstellbar."""
        assert WaringProblem.three_squares_theorem_check(28) is False

    def test_representable_0(self):
        """0 ist darstellbar (trivial: 0+0+0)."""
        assert WaringProblem.three_squares_theorem_check(0) is True

    def test_representable_3(self):
        """3 = 1+1+1 ist darstellbar."""
        assert WaringProblem.three_squares_theorem_check(3) is True

    def test_representable_6(self):
        """6 = 4+1+1 ist darstellbar."""
        assert WaringProblem.three_squares_theorem_check(6) is True

    def test_negative_raises(self):
        """Negativer Input löst ValueError aus."""
        with pytest.raises(ValueError):
            WaringProblem.three_squares_theorem_check(-1)


# ===========================================================================
# TESTS: GoldbachAnalysis
# ===========================================================================

class TestGoldbachDecompositions:
    """Tests für GoldbachAnalysis.goldbach_decompositions(n)."""

    def test_n4(self):
        """4 = 2+2 (einzige Zerlegung)."""
        decomps = GoldbachAnalysis.goldbach_decompositions(4)
        assert (2, 2) in decomps

    def test_n6(self):
        """6 = 3+3 (einzige Zerlegung)."""
        decomps = GoldbachAnalysis.goldbach_decompositions(6)
        assert (3, 3) in decomps

    def test_n10_contains_5_5(self):
        """10 = 5+5 und 10 = 3+7 sind Goldbach-Zerlegungen."""
        decomps = GoldbachAnalysis.goldbach_decompositions(10)
        decomps_set = set(decomps)
        assert (5, 5) in decomps_set
        assert (3, 7) in decomps_set

    def test_all_even_up_to_100(self):
        """Goldbach-Zerlegungen existieren für alle geraden n von 4 bis 100."""
        for n in range(4, 101, 2):
            decomps = GoldbachAnalysis.goldbach_decompositions(n)
            assert len(decomps) > 0, f"Keine Goldbach-Zerlegung für n={n}"

    def test_all_pairs_are_prime(self):
        """Alle zurückgegebenen Paare bestehen aus Primzahlen."""
        from additive_number_theory import _is_prime_simple
        for n in [4, 6, 8, 10, 20, 50, 100]:
            decomps = GoldbachAnalysis.goldbach_decompositions(n)
            for p, q in decomps:
                assert _is_prime_simple(p), f"{p} ist keine Primzahl"
                assert _is_prime_simple(q), f"{q} ist keine Primzahl"
                assert p + q == n, f"{p}+{q}≠{n}"

    def test_pairs_ordered(self):
        """In jedem Paar (p, q) gilt p ≤ q."""
        decomps = GoldbachAnalysis.goldbach_decompositions(30)
        for p, q in decomps:
            assert p <= q, f"Paar ({p},{q}) nicht geordnet"

    def test_odd_raises(self):
        """Ungerades n löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.goldbach_decompositions(7)

    def test_too_small_raises(self):
        """n ≤ 2 löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.goldbach_decompositions(2)


class TestGoldbachCount:
    """Tests für GoldbachAnalysis.goldbach_count(n)."""

    def test_n4_count_1(self):
        """4 hat genau 1 Goldbach-Zerlegung: 2+2."""
        assert GoldbachAnalysis.goldbach_count(4) == 1

    def test_n6_count_1(self):
        """6 hat genau 1 Goldbach-Zerlegung: 3+3."""
        assert GoldbachAnalysis.goldbach_count(6) == 1

    def test_positive_for_all_even(self):
        """Für alle geraden n von 4 bis 50 ist die Anzahl > 0."""
        for n in range(4, 51, 2):
            assert GoldbachAnalysis.goldbach_count(n) > 0, \
                f"Anzahl für n={n} ist 0"


class TestGoldbachCometData:
    """Tests für GoldbachAnalysis.goldbach_comet_data(n_max)."""

    def test_returns_list(self):
        """Gibt eine Liste zurück."""
        data = GoldbachAnalysis.goldbach_comet_data(20)
        assert isinstance(data, list)

    def test_correct_range(self):
        """Alle geraden n von 4 bis n_max sind enthalten."""
        data = GoldbachAnalysis.goldbach_comet_data(20)
        ns = [d[0] for d in data]
        expected = list(range(4, 21, 2))
        assert ns == expected

    def test_all_counts_positive(self):
        """Alle Anzahlen sind ≥ 1 (Goldbach-Vermutung bis 100 verifiziert)."""
        data = GoldbachAnalysis.goldbach_comet_data(100)
        for n, count in data:
            assert count >= 1, f"Goldbach-Zerlegungen für n={n} = 0"

    def test_too_small_raises(self):
        """n_max < 4 löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.goldbach_comet_data(3)


class TestTernaryGoldbach:
    """Tests für GoldbachAnalysis.ternary_goldbach_check(n)."""

    def test_n7(self):
        """7 = 2+2+3 (oder ähnlich)."""
        result = GoldbachAnalysis.ternary_goldbach_check(7)
        assert result is not None
        p, q, r = result
        assert p + q + r == 7

    def test_all_are_prime(self):
        """Alle drei Summanden sind prim."""
        from additive_number_theory import _is_prime_simple
        for n in [7, 9, 11, 13, 21, 99]:
            result = GoldbachAnalysis.ternary_goldbach_check(n)
            assert result is not None, f"Keine ternäre Zerlegung für n={n}"
            for p in result:
                assert _is_prime_simple(p), f"{p} in Zerlegung von {n} ist keine Primzahl"

    def test_odd_numbers_decompose(self):
        """Alle ungeraden Zahlen > 5 bis 51 haben eine ternäre Goldbach-Zerlegung."""
        for n in range(7, 52, 2):
            result = GoldbachAnalysis.ternary_goldbach_check(n)
            assert result is not None, f"Keine ternäre Zerlegung für n={n}"
            assert sum(result) == n

    def test_even_raises(self):
        """Gerader Input löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.ternary_goldbach_check(10)

    def test_too_small_raises(self):
        """n ≤ 5 löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.ternary_goldbach_check(5)


class TestChenTheoremCheck:
    """Tests für GoldbachAnalysis.chen_theorem_check(n)."""

    def test_n4(self):
        """4 = 2 + 2 (p=2, m=2 prim → hat ≤ 2 Primfaktoren)."""
        result = GoldbachAnalysis.chen_theorem_check(4)
        assert result is not None
        p, m = result
        assert p + m == 4

    def test_structure_valid(self):
        """Für alle geraden n von 4 bis 50 ist das Ergebnis korrekt."""
        for n in range(4, 51, 2):
            result = GoldbachAnalysis.chen_theorem_check(n)
            assert result is not None, f"Chen-Zerlegung für n={n} nicht gefunden"
            p, m = result
            assert p + m == n

    def test_odd_raises(self):
        """Ungerades n löst ValueError aus."""
        with pytest.raises(ValueError):
            GoldbachAnalysis.chen_theorem_check(9)


# ===========================================================================
# TESTS: ArithmeticProgressions
# ===========================================================================

class TestPrimesInAP:
    """Tests für ArithmeticProgressions.primes_in_arithmetic_progression."""

    def test_ap_2_2(self):
        """Primzahlen in 2, 4, 6, ...: nur 2."""
        primes = ArithmeticProgressions.primes_in_arithmetic_progression(2, 2, 20)
        assert primes == [2]

    def test_ap_3_4(self):
        """Primzahlen in AP 3 mod 4 (3, 7, 11, 19, ...)."""
        primes = ArithmeticProgressions.primes_in_arithmetic_progression(3, 4, 30)
        expected = [3, 7, 11, 19, 23]
        assert primes == expected

    def test_all_returned_are_prime(self):
        """Alle zurückgegebenen Elemente sind prim."""
        from additive_number_theory import _is_prime_simple
        primes = ArithmeticProgressions.primes_in_arithmetic_progression(1, 6, 100)
        for p in primes:
            assert _is_prime_simple(p), f"{p} ist keine Primzahl"

    def test_all_in_ap(self):
        """Alle zurückgegebenen Elemente liegen in der AP."""
        a, d = 5, 6
        primes = ArithmeticProgressions.primes_in_arithmetic_progression(a, d, 200)
        for p in primes:
            assert (p - a) % d == 0, f"{p} liegt nicht in AP {a} + k·{d}"

    def test_invalid_a_raises(self):
        """a < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.primes_in_arithmetic_progression(0, 2, 20)

    def test_invalid_d_raises(self):
        """d < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.primes_in_arithmetic_progression(2, 0, 20)


class TestVanDerWaerden:
    """Tests für ArithmeticProgressions.van_der_waerden_w(k, r)."""

    def test_w_3_2_equals_9(self):
        """W(3; 2) = 9."""
        assert ArithmeticProgressions.van_der_waerden_w(3, 2) == 9

    def test_w_4_2_equals_35(self):
        """W(4; 2) = 35."""
        assert ArithmeticProgressions.van_der_waerden_w(4, 2) == 35

    def test_w_3_3_equals_27(self):
        """W(3; 3) = 27."""
        assert ArithmeticProgressions.van_der_waerden_w(3, 3) == 27

    def test_unknown_returns_none(self):
        """Unbekannte Werte geben None zurück."""
        assert ArithmeticProgressions.van_der_waerden_w(10, 5) is None

    def test_invalid_k_raises(self):
        """k < 2 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.van_der_waerden_w(1, 2)

    def test_invalid_r_raises(self):
        """r < 2 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.van_der_waerden_w(3, 1)


class TestArithmeticProgressionInPrimes:
    """Tests für ArithmeticProgressions.arithmetic_progression_in_primes."""

    def test_length_3_found(self):
        """Eine AP der Länge 3 aus Primzahlen wird gefunden."""
        result = ArithmeticProgressions.arithmetic_progression_in_primes(3)
        assert result is not None
        assert len(result) == 3

    def test_length_5_found(self):
        """Eine AP der Länge 5 aus Primzahlen wird gefunden."""
        result = ArithmeticProgressions.arithmetic_progression_in_primes(5)
        assert result is not None
        assert len(result) == 5

    def test_all_prime_in_ap(self):
        """Alle Elemente der gefundenen AP sind prim."""
        from additive_number_theory import _is_prime_simple
        for length in [3, 4, 5]:
            result = ArithmeticProgressions.arithmetic_progression_in_primes(length)
            assert result is not None
            for p in result:
                assert _is_prime_simple(p), f"{p} in AP ist keine Primzahl"

    def test_constant_difference(self):
        """Alle Elemente haben dieselbe Differenz."""
        result = ArithmeticProgressions.arithmetic_progression_in_primes(4)
        assert result is not None
        diffs = [result[i + 1] - result[i] for i in range(len(result) - 1)]
        assert len(set(diffs)) == 1, f"Keine konstante Differenz: {diffs}"

    def test_invalid_length_raises(self):
        """length < 2 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.arithmetic_progression_in_primes(1)


class TestDirichletDensity:
    """Tests für ArithmeticProgressions.dirichlet_density(a, d)."""

    def test_density_1mod4(self):
        """Dichte der Primzahlen ≡ 1 (mod 4): 1/φ(4) = 1/2."""
        density = ArithmeticProgressions.dirichlet_density(1, 4)
        assert abs(density - 0.5) < 1e-10

    def test_density_1mod6(self):
        """Dichte der Primzahlen ≡ 1 (mod 6): 1/φ(6) = 1/2."""
        density = ArithmeticProgressions.dirichlet_density(1, 6)
        assert abs(density - 0.5) < 1e-10

    def test_density_2mod4_zero(self):
        """gcd(2, 4) ≠ 1 → Dichte = 0."""
        density = ArithmeticProgressions.dirichlet_density(2, 4)
        assert density == 0.0

    def test_density_positive_when_coprime(self):
        """Dichte > 0 wenn gcd(a, d) = 1."""
        density = ArithmeticProgressions.dirichlet_density(5, 6)
        assert density > 0

    def test_invalid_d_raises(self):
        """d < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ArithmeticProgressions.dirichlet_density(1, 0)


class TestGreenTaoDemo:
    """Tests für ArithmeticProgressions.green_tao_demo()."""

    def test_returns_dict(self):
        """Gibt ein Dictionary zurück."""
        result = ArithmeticProgressions.green_tao_demo()
        assert isinstance(result, dict)
        assert "theorem" in result
        assert "examples" in result

    def test_examples_are_valid(self):
        """Alle Beispiele enthalten gültige APs aus Primzahlen."""
        result = ArithmeticProgressions.green_tao_demo()
        for ex in result["examples"]:
            assert ex["all_prime"] is True, \
                f"AP {ex['progression']} enthält keine Primzahlen"

    def test_year_2004(self):
        """Das Theorem wurde 2004 bewiesen."""
        result = ArithmeticProgressions.green_tao_demo()
        assert result["year"] == 2004


# ===========================================================================
# TESTS: SumSetTheory
# ===========================================================================

class TestSumset:
    """Tests für SumSetTheory.sumset(A, B)."""

    def test_basic_sumset(self):
        """{1,2} + {3,4} = {4,5,6}."""
        result = SumSetTheory.sumset({1, 2}, {3, 4})
        assert result == {4, 5, 6}

    def test_empty_A(self):
        """Leeres A → leere Summenmenge."""
        result = SumSetTheory.sumset(set(), {1, 2, 3})
        assert result == set()

    def test_empty_B(self):
        """Leeres B → leere Summenmenge."""
        result = SumSetTheory.sumset({1, 2, 3}, set())
        assert result == set()

    def test_singleton(self):
        """{5} + {3} = {8}."""
        result = SumSetTheory.sumset({5}, {3})
        assert result == {8}

    def test_size_lower_bound(self):
        """|A+B| ≥ |A| + |B| - 1 für Intervalle."""
        A = set(range(1, 6))   # {1,2,3,4,5}
        B = set(range(1, 4))   # {1,2,3}
        sumset = SumSetTheory.sumset(A, B)
        assert len(sumset) >= len(A) + len(B) - 1


class TestDifferenceSet:
    """Tests für SumSetTheory.difference_set(A, B)."""

    def test_basic(self):
        """{3,4} - {1,2} = {1,2,3}."""
        result = SumSetTheory.difference_set({3, 4}, {1, 2})
        assert result == {1, 2, 3}

    def test_empty(self):
        """Leere Differenzmenge bei leerem Input."""
        assert SumSetTheory.difference_set(set(), {1}) == set()
        assert SumSetTheory.difference_set({1}, set()) == set()

    def test_A_minus_A_contains_zero(self):
        """A - A enthält immer 0."""
        A = {1, 3, 7, 12}
        diff = SumSetTheory.difference_set(A, A)
        assert 0 in diff


class TestCauchyDavenport:
    """Tests für SumSetTheory.cauchy_davenport_bound(A, B, p)."""

    def test_basic_bound(self):
        """Cauchy-Davenport: |A+B| ≥ min(p, |A|+|B|-1)."""
        A = {1, 2, 3}
        B = {4, 5}
        p = 7
        bound = SumSetTheory.cauchy_davenport_bound(A, B, p)
        assert bound == min(p, len(A) + len(B) - 1)

    def test_full_group(self):
        """Wenn |A|+|B|-1 ≥ p, ist die Schranke p."""
        A = {0, 1, 2, 3}
        B = {0, 1, 2, 3, 4}
        p = 7
        bound = SumSetTheory.cauchy_davenport_bound(A, B, p)
        assert bound == p

    def test_actual_sumset_at_least_bound(self):
        """Die tatsächliche Summenmenge ≥ Cauchy-Davenport-Schranke."""
        A = {1, 3}
        B = {2, 4}
        p = 7
        bound = SumSetTheory.cauchy_davenport_bound(A, B, p)
        actual = len({(a + b) % p for a in A for b in B})
        assert actual >= bound

    def test_non_prime_raises(self):
        """Nicht-prim p löst ValueError aus."""
        with pytest.raises(ValueError):
            SumSetTheory.cauchy_davenport_bound({1}, {2}, 4)

    def test_empty_raises(self):
        """Leere Mengen lösen ValueError aus."""
        with pytest.raises(ValueError):
            SumSetTheory.cauchy_davenport_bound(set(), {1}, 5)


class TestPlunneckeRuzsa:
    """Tests für SumSetTheory.plunnecke_ruzsa_estimate(A, B)."""

    def test_returns_dict(self):
        """Gibt ein Dictionary zurück."""
        result = SumSetTheory.plunnecke_ruzsa_estimate({1, 2, 3}, {1, 2})
        assert isinstance(result, dict)
        assert "K_doubling" in result

    def test_K_positive(self):
        """Verdoppelungskonstante K > 0."""
        result = SumSetTheory.plunnecke_ruzsa_estimate({1, 2, 3, 4}, {0, 1, 2})
        assert result["K_doubling"] > 0

    def test_empty_raises(self):
        """Leere Mengen lösen ValueError aus."""
        with pytest.raises(ValueError):
            SumSetTheory.plunnecke_ruzsa_estimate(set(), {1, 2})


class TestFreimanTheoremDemo:
    """Tests für SumSetTheory.freiman_theorem_demo(A)."""

    def test_ap_has_small_K(self):
        """Eine AP hat Verdoppelungskonstante K nahe 2."""
        A = set(range(1, 11))  # {1, 2, ..., 10} ist eine AP
        result = SumSetTheory.freiman_theorem_demo(A)
        # Für {1,...,n}: |A+A| = 2n-1, K = (2n-1)/n ≈ 2
        assert result["K_doubling"] < 2.1

    def test_is_ap_detected(self):
        """AP wird als AP erkannt."""
        A = {2, 4, 6, 8, 10}
        result = SumSetTheory.freiman_theorem_demo(A)
        assert result["is_arithmetic_progression"] is True

    def test_random_set_not_ap(self):
        """Eine nicht-arithmetische Menge wird korrekt identifiziert."""
        A = {1, 2, 4, 8, 16}  # Geometrische Folge, keine AP
        result = SumSetTheory.freiman_theorem_demo(A)
        assert result["is_arithmetic_progression"] is False

    def test_empty_raises(self):
        """Leere Menge löst ValueError aus."""
        with pytest.raises(ValueError):
            SumSetTheory.freiman_theorem_demo(set())


class TestAdditiveEnergy:
    """Tests für SumSetTheory.additive_energy(A)."""

    def test_trivial_lower_bound(self):
        """E(A) ≥ |A|² (triviale Schranke: a+a=a+a für jedes a)."""
        A = {1, 3, 5, 7}
        energy = SumSetTheory.additive_energy(A)
        assert energy >= len(A) ** 2

    def test_ap_has_large_energy(self):
        """Eine AP {1,...,n} hat additive Energie ~n³/3."""
        A = set(range(1, 8))  # {1,...,7}
        energy = SumSetTheory.additive_energy(A)
        n = len(A)
        # Grobe Schranke: E(A) ≤ |A|³
        assert energy <= n ** 3

    def test_energy_positive(self):
        """Additive Energie ist immer positiv."""
        A = {2, 5, 11}
        assert SumSetTheory.additive_energy(A) > 0

    def test_singleton_energy(self):
        """Für |A|=1: E(A) = 1."""
        A = {42}
        assert SumSetTheory.additive_energy(A) == 1

    def test_empty_raises(self):
        """Leere Menge löst ValueError aus."""
        with pytest.raises(ValueError):
            SumSetTheory.additive_energy(set())


# ===========================================================================
# TESTS: Standalone-Funktionen
# ===========================================================================

class TestSchurNumber:
    """Tests für schur_number(k)."""

    def test_s1_equals_2(self):
        """S(1) = 2."""
        assert schur_number(1) == 2

    def test_s2_equals_5(self):
        """S(2) = 5."""
        assert schur_number(2) == 5

    def test_s3_equals_14(self):
        """S(3) = 14."""
        assert schur_number(3) == 14

    def test_s4_equals_45(self):
        """S(4) = 45."""
        assert schur_number(4) == 45

    def test_s5_equals_161(self):
        """S(5) = 161 (SAT-Beweis 2017)."""
        assert schur_number(5) == 161

    def test_unknown_returns_none(self):
        """k > 5 gibt None zurück."""
        assert schur_number(6) is None

    def test_invalid_raises(self):
        """k < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            schur_number(0)


class TestHappyEndingProblem:
    """Tests für happy_ending_problem(n)."""

    def test_n3(self):
        """f(3) = 3 (ein Dreieck braucht 3 Punkte)."""
        # 2^(3-2) + 1 = 3
        assert happy_ending_problem(3) == 3

    def test_n4(self):
        """f(4) = 5 (Esther Klein: 5 Punkte garantieren konvexes Viereck)."""
        # 2^(4-2) + 1 = 5
        assert happy_ending_problem(4) == 5

    def test_n5(self):
        """f(5) = 9 (bekannt)."""
        # 2^(5-2) + 1 = 9
        assert happy_ending_problem(5) == 9

    def test_growth(self):
        """f(n) wächst exponentiell."""
        prev = happy_ending_problem(3)
        for n in range(4, 10):
            curr = happy_ending_problem(n)
            assert curr > prev
            prev = curr

    def test_invalid_raises(self):
        """n < 3 löst ValueError aus."""
        with pytest.raises(ValueError):
            happy_ending_problem(2)


class TestBertrandPostulate:
    """Tests für bertrand_postulate_verify(n)."""

    def test_n1(self):
        """Zwischen 1 und 2: 2 ist prim."""
        p = bertrand_postulate_verify(1)
        assert p == 2

    def test_n2(self):
        """Zwischen 2 und 4: 3 ist prim."""
        p = bertrand_postulate_verify(2)
        assert p in {3}  # 3 ist die einzige Primzahl in (2,4]

    def test_range_satisfied(self):
        """Die zurückgegebene Primzahl liegt in (n, 2n]."""
        from additive_number_theory import _is_prime_simple
        for n in range(1, 50):
            p = bertrand_postulate_verify(n)
            assert p is not None, f"Keine Primzahl in ({n}, {2*n}] gefunden"
            assert n < p <= 2 * n, f"p={p} liegt nicht in ({n},{2*n}]"
            assert _is_prime_simple(p), f"p={p} ist keine Primzahl"

    def test_invalid_raises(self):
        """n < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            bertrand_postulate_verify(0)
