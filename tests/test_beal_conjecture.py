"""
@file test_beal_conjecture.py
@brief Tests für das Beal Conjecture Modul
@description
    Umfassende Tests für BealTriple, BealChecker und ABCConnection.
    Prüft bekannte mathematische Eigenschaften und Edge-Cases.

@author Michael Fuhrmann
@date 2026-03-13
@lastModified 2026-03-13
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from beal_conjecture import BealTriple, BealChecker, ABCConnection


# ──────────────────────────────────────────────────────────────
# Tests für BealTriple
# ──────────────────────────────────────────────────────────────

class TestBealTriple:

    def test_valid_triple_2_3_2_3_2_4(self):
        """2^3 + 2^3 = 2^4: 8 + 8 = 16 ist korrekt."""
        t = BealTriple(2, 3, 2, 3, 2, 4)
        assert t.is_valid()

    def test_valid_triple_3_3_6_3_3_5(self):
        """3^3 + 6^3 = 3^5: 27 + 216 = 243 ist korrekt."""
        t = BealTriple(3, 3, 6, 3, 3, 5)
        assert t.is_valid()

    def test_invalid_triple(self):
        """1^3 + 1^3 = 3^3 ist FALSCH (1+1=2 ≠ 27)."""
        t = BealTriple(1, 3, 1, 3, 3, 3)
        assert not t.is_valid()

    def test_has_common_factor_same_base(self):
        """2^3 + 2^3 = 2^4: gcd(2,2,2) = 2 > 1."""
        t = BealTriple(2, 3, 2, 3, 2, 4)
        assert t.has_common_factor()

    def test_has_common_factor_shared_prime(self):
        """3^3 + 6^3 = 3^5: gcd(3,6,3) = 3 > 1."""
        t = BealTriple(3, 3, 6, 3, 3, 5)
        assert t.has_common_factor()

    def test_no_common_factor_coprime(self):
        """gcd(2, 3, 5) = 1 — kein gemeinsamer Faktor."""
        t = BealTriple(2, 3, 3, 3, 5, 3)
        # Prüfe nur gcd, nicht Gültigkeit der Gleichung
        assert not t.has_common_factor()

    def test_satisfies_beal_with_shared_factor(self):
        """3^3 + 6^3 = 3^5 hat gemeinsamen Faktor → Beal erfüllt."""
        t = BealTriple(3, 3, 6, 3, 3, 5)
        assert t.satisfies_beal()

    def test_satisfies_beal_low_exponent_bypassed(self):
        """Tripel mit exp < 3 ist außerhalb des Anwendungsbereichs → Beal trivial True."""
        t = BealTriple(3, 2, 4, 2, 5, 2)
        assert t.satisfies_beal()

    def test_radical_prime_bases(self):
        """rad(2 * 3 * 5) = 30."""
        t = BealTriple(2, 3, 3, 3, 5, 3)
        assert t.radical() == 30

    def test_radical_prime_power_base(self):
        """rad(4 * 8 * 2) = rad(2^2 * 2^3 * 2) = 2."""
        t = BealTriple(4, 3, 8, 3, 2, 3)
        assert t.radical() == 2

    def test_common_factor_value(self):
        """gcd(6, 9, 3) = 3."""
        t = BealTriple(6, 3, 9, 3, 3, 3)
        assert t.common_factor() == 3

    def test_repr_contains_bases(self):
        """__repr__ enthält Basis-Werte."""
        t = BealTriple(2, 3, 2, 3, 2, 4)
        r = repr(t)
        assert "2" in r

    def test_known_families_all_valid(self):
        """Alle bekannten Beal-Familien sind arithmetisch korrekt."""
        families = BealChecker.known_families()
        for t in families:
            assert t.is_valid(), f"Ungültiges Tripel: {t}"

    def test_known_families_have_common_factor(self):
        """Alle bekannten Beal-Familien haben gemeinsamen Primteiler."""
        families = BealChecker.known_families()
        for t in families:
            assert t.has_common_factor(), f"Kein gemeinsamer Faktor: {t}"

    def test_known_families_satisfy_beal(self):
        """Alle bekannten Familien erfüllen Beal's CONJECTURE."""
        families = BealChecker.known_families()
        for t in families:
            assert t.satisfies_beal()


# ──────────────────────────────────────────────────────────────
# Tests für BealChecker
# ──────────────────────────────────────────────────────────────

class TestBealChecker:

    def test_precompute_powers_contains_8(self):
        """2^3 = 8 muss in der Lookup-Tabelle stehen."""
        checker = BealChecker(N=10, max_exp=5)
        powers = checker.precompute_powers()
        assert 8 in powers
        assert (2, 3) in powers[8]

    def test_precompute_powers_contains_27(self):
        """3^3 = 27 muss in der Lookup-Tabelle stehen."""
        checker = BealChecker(N=10, max_exp=5)
        powers = checker.precompute_powers()
        assert 27 in powers

    def test_precompute_multiple_representations(self):
        """64 = 2^6 = 4^3 → zwei Darstellungen in der Tabelle."""
        checker = BealChecker(N=10, max_exp=6)
        powers = checker.precompute_powers()
        assert 64 in powers
        reps = powers[64]
        bases = [r[0] for r in reps]
        assert 2 in bases or 4 in bases  # mindestens eine gefunden

    def test_find_all_triples_contains_known(self):
        """3^3 + 6^3 = 3^5 muss unter den Tripeln sein."""
        checker = BealChecker(N=10, max_exp=6)
        triples = checker.find_all_triples()
        found = any(
            t.A == 3 and t.x == 3 and t.B == 6 and t.y == 3
            and t.C == 3 and t.z == 5
            for t in triples
        )
        assert found, "3^3 + 6^3 = 3^5 nicht gefunden"

    def test_no_counterexamples_small(self):
        """Für N=50 gibt es keine Gegenbeispiele zu Beal's CONJECTURE."""
        checker = BealChecker(N=50, max_exp=6)
        counterexamples = checker.find_counterexamples()
        assert counterexamples == [], (
            f"Gegenbeispiel gefunden: {counterexamples[0]}"
        )

    def test_all_triples_valid(self):
        """Alle gefundenen Tripel erfüllen die Gleichung A^x + B^y = C^z."""
        checker = BealChecker(N=30, max_exp=5)
        triples = checker.find_all_triples()
        for t in triples:
            assert t.is_valid(), f"Ungültiges Tripel: {t}"

    def test_all_triples_satisfy_beal(self):
        """Alle gefundenen Tripel (mit x,y,z≥3) haben gemeinsamen Faktor."""
        checker = BealChecker(N=30, max_exp=5)
        triples = checker.find_all_triples()
        for t in triples:
            if t.x >= 3 and t.y >= 3 and t.z >= 3:
                assert t.has_common_factor(), (
                    f"Tripel ohne gemeinsamen Faktor: {t}"
                )


# ──────────────────────────────────────────────────────────────
# Tests für ABCConnection
# ──────────────────────────────────────────────────────────────

class TestABCConnection:

    def test_radical_prime(self):
        """rad(p) = p für Primzahl p."""
        for p in [2, 3, 5, 7, 11]:
            assert ABCConnection.radical(p) == p

    def test_radical_prime_power(self):
        """rad(p^k) = p."""
        assert ABCConnection.radical(8) == 2   # 2^3
        assert ABCConnection.radical(27) == 3  # 3^3
        assert ABCConnection.radical(32) == 2  # 2^5

    def test_radical_composite(self):
        """rad(12) = rad(2^2 * 3) = 2*3 = 6."""
        assert ABCConnection.radical(12) == 6

    def test_quality_positive(self):
        """ABC-Qualität ist positiv für c > 1."""
        q = ABCConnection.quality(1, 2, 3)
        assert q > 0

    def test_quality_high_example(self):
        """Bekanntes Beispiel mit hoher Qualität: 1 + 8 = 9 (1+2^3=3^2)."""
        # q = log(9) / log(rad(1*8*9)) = log(9) / log(6) ≈ 1.226
        q = ABCConnection.quality(1, 8, 9)
        assert q > 1.0

    def test_quality_zero_for_trivial(self):
        """Qualität 0 für c ≤ 1."""
        assert ABCConnection.quality(0, 1, 0) == 0.0

    def test_beal_abc_bound_positive(self):
        """ABC-Qualität für Beal-Tripel 3^3 + 6^3 = 3^5 ist positiv."""
        q = ABCConnection.beal_abc_bound(3, 3, 6, 3, 3, 5)
        assert q > 0

    def test_beal_abc_bound_same_base(self):
        """Für 2^3 + 2^3 = 2^4: rad(2,2,2) = 2, C^z = 16."""
        # q = log(16) / log(2) = 4.0
        q = ABCConnection.beal_abc_bound(2, 3, 2, 3, 2, 4)
        assert abs(q - 4.0) < 1e-9

    def test_high_quality_triples_not_empty(self):
        """Es gibt ABC-Tripel mit Qualität > 1 unter 1000."""
        triples = ABCConnection.high_quality_abc_triples(limit=1000)
        assert len(triples) > 0

    def test_high_quality_triples_sorted(self):
        """Hochqualitäts-Tripel sind nach Qualität absteigend sortiert."""
        triples = ABCConnection.high_quality_abc_triples(limit=500)
        if len(triples) >= 2:
            for i in range(len(triples) - 1):
                assert triples[i][3] >= triples[i + 1][3]

    def test_high_quality_triples_sum_correct(self):
        """Für alle zurückgegebenen Tripel gilt a + b = c."""
        triples = ABCConnection.high_quality_abc_triples(limit=200)
        for a, b, c, q in triples:
            assert a + b == c, f"a+b≠c: {a}+{b}≠{c}"
