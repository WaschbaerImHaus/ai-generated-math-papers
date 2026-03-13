"""
@file test_frankl_union_closed.py
@brief Tests für das Frankl Union-Closed Conjecture Modul
@description
    Umfassende Tests für UnionClosedFamily, FranklChecker und UnionClosedGenerator.
    Alle Testfälle prüfen bekannte mathematische Eigenschaften.

@author Michael Fuhrmann
@date 2026-03-13
@lastModified 2026-03-13
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from frankl_union_closed import (
    UnionClosedFamily,
    FranklChecker,
    UnionClosedGenerator,
)


# ──────────────────────────────────────────────────────────────
# Hilfsfunktionen
# ──────────────────────────────────────────────────────────────

def make_family(list_of_sets):
    """Erstellt UnionClosedFamily aus einer Liste von gewöhnlichen Mengen."""
    return UnionClosedFamily([frozenset(s) for s in list_of_sets])


# ──────────────────────────────────────────────────────────────
# Tests für UnionClosedFamily.__init__ und Grundeigenschaften
# ──────────────────────────────────────────────────────────────

class TestUnionClosedFamilyInit:

    def test_empty_family(self):
        """Leere Familie hat 0 Elemente."""
        fam = UnionClosedFamily([])
        assert len(fam) == 0

    def test_single_set(self):
        """Familie mit einer Menge hat Länge 1."""
        fam = make_family([{1, 2}])
        assert len(fam) == 1

    def test_deduplication(self):
        """Duplikate werden entfernt."""
        fam = make_family([{1}, {1}, {2}])
        assert len(fam) == 2

    def test_universe_computed(self):
        """Grundmenge wird korrekt berechnet."""
        fam = make_family([{1, 2}, {2, 3}])
        assert fam.universe == frozenset({1, 2, 3})

    def test_from_bitmasks_single(self):
        """Einzelne Bitmaske 0b011 = {0,1} über n=3."""
        fam = UnionClosedFamily.from_bitmasks(3, [0b011])
        assert frozenset({0, 1}) in fam.sets

    def test_from_bitmasks_empty_set(self):
        """Bitmaske 0 entspricht der leeren Menge."""
        fam = UnionClosedFamily.from_bitmasks(3, [0])
        assert frozenset() in fam.sets

    def test_from_bitmasks_multiple(self):
        """Mehrere Bitmasken erzeugen korrekte Mengen."""
        fam = UnionClosedFamily.from_bitmasks(3, [0b001, 0b010, 0b011])
        assert frozenset({0}) in fam.sets
        assert frozenset({1}) in fam.sets
        assert frozenset({0, 1}) in fam.sets


# ──────────────────────────────────────────────────────────────
# Tests für is_union_closed
# ──────────────────────────────────────────────────────────────

class TestIsUnionClosed:

    def test_power_set_is_uc(self):
        """Potenzmenge ist vereinigungsabgeschlossen."""
        fam = make_family([set(), {1}, {2}, {1, 2}])
        assert fam.is_union_closed()

    def test_chain_is_uc(self):
        """Kette {{1},{1,2},{1,2,3}} ist uc-abgeschlossen."""
        fam = make_family([{1}, {1, 2}, {1, 2, 3}])
        assert fam.is_union_closed()

    def test_not_uc(self):
        """Familie {{1},{2}} ist NICHT uc-abgeschlossen ({1}∪{2}={1,2} fehlt)."""
        fam = make_family([{1}, {2}])
        assert not fam.is_union_closed()

    def test_single_set_is_uc(self):
        """Eine Einzelmenge ist immer uc-abgeschlossen."""
        fam = make_family([{1, 2, 3}])
        assert fam.is_union_closed()

    def test_empty_family_is_uc(self):
        """Leere Familie ist trivialerweise uc-abgeschlossen."""
        fam = UnionClosedFamily([])
        assert fam.is_union_closed()


# ──────────────────────────────────────────────────────────────
# Tests für frequency_vector und max_frequency_ratio
# ──────────────────────────────────────────────────────────────

class TestFrequencyVector:

    def test_frequency_all_sets_contain_element(self):
        """Wenn alle Mengen Element 1 enthalten, ist freq[0] = |F|."""
        fam = make_family([{1}, {1, 2}, {1, 3}])
        # Universe = {1,2,3}, sorted = [1,2,3]
        freq = fam.frequency_vector()
        assert freq[0] == 3  # Element 1 in allen 3 Mengen

    def test_frequency_length_equals_universe(self):
        """Länge des Häufigkeitsvektors = |Grundmenge|."""
        fam = make_family([{1, 2}, {2, 3}, {1, 3}])
        freq = fam.frequency_vector()
        assert len(freq) == 3

    def test_frequency_empty_family(self):
        """Leere Familie hat leeren Häufigkeitsvektor."""
        fam = UnionClosedFamily([])
        assert fam.frequency_vector() == []

    def test_max_frequency_ratio_power_set(self):
        """In der Potenzmenge von {0,1} kommt jedes Element in 2 von 4 Mengen vor."""
        fam = UnionClosedGenerator.power_set(2)
        ratio = fam.max_frequency_ratio()
        # 2^1 = 2 von 2^2 = 4 Mengen → Ratio = 0.5
        assert abs(ratio - 0.5) < 1e-10

    def test_max_frequency_ratio_chain(self):
        """Kette {{0},{0,1},{0,1,2}}: Element 0 in allen 3 Mengen → ratio = 1.0."""
        fam = UnionClosedGenerator.chain_family(3)
        assert abs(fam.max_frequency_ratio() - 1.0) < 1e-10


# ──────────────────────────────────────────────────────────────
# Tests für satisfies_frankl und frankl_witness
# ──────────────────────────────────────────────────────────────

class TestFranklProperty:

    def test_power_set_satisfies_frankl(self):
        """Potenzmenge erfüllt Frankl's CONJECTURE."""
        for n in range(1, 5):
            fam = UnionClosedGenerator.power_set(n)
            assert fam.satisfies_frankl(), f"Potenzmenge n={n} verletzt Frankl"

    def test_chain_satisfies_frankl(self):
        """Ketten erfüllen Frankl's CONJECTURE (Element 0 in allen Mengen)."""
        for n in range(1, 6):
            fam = UnionClosedGenerator.chain_family(n)
            assert fam.satisfies_frankl()

    def test_interval_family_satisfies_frankl(self):
        """Intervall-Familie erfüllt Frankl's CONJECTURE."""
        for n in range(1, 5):
            fam = UnionClosedGenerator.interval_family(n)
            assert fam.satisfies_frankl()

    def test_frankl_witness_found(self):
        """Für Kette der Länge 3 ist Zeuge vorhanden."""
        fam = UnionClosedGenerator.chain_family(3)
        witness = fam.frankl_witness()
        assert witness is not None
        # Zeuge muss in >= |F|/2 Mengen vorkommen
        count = sum(1 for s in fam.sets if witness in s)
        assert count >= len(fam) / 2

    def test_only_empty_set_is_trivial(self):
        """Familie {{}} (nur leere Menge) gilt als trivial erfüllt."""
        fam = make_family([set()])
        assert fam.satisfies_frankl()

    def test_frankl_witness_none_for_empty_universe(self):
        """Bei leerer Grundmenge gibt es keinen Zeugen."""
        fam = make_family([set()])
        witness = fam.frankl_witness()
        assert witness is None


# ──────────────────────────────────────────────────────────────
# Tests für FranklChecker
# ──────────────────────────────────────────────────────────────

class TestFranklChecker:

    def test_verify_n1_no_violations(self):
        """Für n=1 gibt es keine Verletzungen von Frankl's CONJECTURE."""
        total, violations, all_ok = FranklChecker.verify_up_to_n(1)
        assert all_ok
        assert violations == 0

    def test_verify_n2_no_violations(self):
        """Für n=2 gibt es keine Verletzungen von Frankl's CONJECTURE."""
        total, violations, all_ok = FranklChecker.verify_up_to_n(2)
        assert all_ok
        assert violations == 0

    def test_verify_n3_no_violations(self):
        """Für n=3 gibt es keine Verletzungen von Frankl's CONJECTURE."""
        total, violations, all_ok = FranklChecker.verify_up_to_n(3)
        assert all_ok
        assert violations == 0

    def test_verify_returns_positive_total(self):
        """Die Gesamtanzahl geprüfter Familien ist positiv."""
        total, _, _ = FranklChecker.verify_up_to_n(2)
        assert total > 0

    def test_gilmer_bound_positive(self):
        """Gilmer-Schranke ist positiv für |F| > 0."""
        assert FranklChecker.gilmer_bound(100) == pytest.approx(1.0)

    def test_chase_lovasz_bound(self):
        """Chase-Lovász-Schranke ist korrekt."""
        assert FranklChecker.chase_lovasz_bound(100) == pytest.approx(38.0)

    def test_gilmer_weaker_than_frankl(self):
        """Gilmer-Schranke (0.01) liegt unter Frankl-Ziel (0.5)."""
        for size in [10, 100, 1000]:
            assert FranklChecker.gilmer_bound(size) < 0.5 * size

    def test_chase_lovasz_weaker_than_frankl(self):
        """Chase-Lovász-Schranke (0.38) liegt unter Frankl-Ziel (0.5)."""
        for size in [10, 100, 1000]:
            assert FranklChecker.chase_lovasz_bound(size) < 0.5 * size


# ──────────────────────────────────────────────────────────────
# Tests für UnionClosedGenerator
# ──────────────────────────────────────────────────────────────

class TestUnionClosedGenerator:

    def test_power_set_size(self):
        """Potenzmenge von {0,...,n-1} hat 2^n Mengen."""
        for n in range(1, 5):
            fam = UnionClosedGenerator.power_set(n)
            assert len(fam) == 2 ** n

    def test_power_set_is_uc(self):
        """Generierte Potenzmenge ist uc-abgeschlossen."""
        for n in range(1, 5):
            fam = UnionClosedGenerator.power_set(n)
            assert fam.is_union_closed()

    def test_interval_family_contains_singletons(self):
        """Intervall-Familie über n=4 enthält alle Singleton-Intervalle."""
        fam = UnionClosedGenerator.interval_family(4)
        for i in range(4):
            assert frozenset({i}) in fam.sets

    def test_interval_family_uc_closure_is_uc(self):
        """
        Die Intervall-Familie selbst ist nur dann uc-abgeschlossen, wenn alle
        Vereinigungen wieder Intervalle sind. Für n=1,2 gilt das; für n≥3 nicht
        unbedingt (z.B. {0}∪{2}={0,2} ist kein Intervall).
        Wir prüfen stattdessen, dass der uc-Abschluss korrekt funktioniert.
        """
        for n in range(1, 4):
            fam = UnionClosedGenerator.interval_family(n)
            # uc-Abschluss der Intervall-Familie ist uc-abgeschlossen
            closed = UnionClosedGenerator.union_closure(fam.sets)
            assert closed.is_union_closed()

    def test_union_closure_of_single(self):
        """Vereinigungsabschluss einer Einzelmenge ist die Menge selbst."""
        fam = UnionClosedGenerator.union_closure([frozenset({1, 2})])
        assert len(fam) == 1
        assert frozenset({1, 2}) in fam.sets

    def test_union_closure_two_disjoint(self):
        """Abschluss von {{1},{2}} enthält auch {1,2}."""
        fam = UnionClosedGenerator.union_closure([frozenset({1}), frozenset({2})])
        assert frozenset({1, 2}) in fam.sets
        assert fam.is_union_closed()

    def test_random_uc_is_union_closed(self):
        """Zufällig erzeugte Familie ist uc-abgeschlossen."""
        for seed in range(5):
            fam = UnionClosedGenerator.random_union_closed(4, seed=seed)
            assert fam.is_union_closed()

    def test_chain_family_size(self):
        """Kette der Länge n hat n Elemente."""
        for n in range(1, 6):
            fam = UnionClosedGenerator.chain_family(n)
            assert len(fam) == n

    def test_chain_family_is_uc(self):
        """Kette ist uc-abgeschlossen."""
        for n in range(1, 6):
            fam = UnionClosedGenerator.chain_family(n)
            assert fam.is_union_closed()
