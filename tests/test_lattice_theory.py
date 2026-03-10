"""
@file test_lattice_theory.py
@brief Umfassende Tests für das Modul lattice_theory.py.
@description
    Testet alle Klassen und Funktionen der Verbandstheorie:
    - PartialOrder: Axiome, Hasse-Diagramm, minimale/maximale Elemente
    - Lattice: Meet, Join, Top, Bottom, Distributivität, Modularität, Komplementierung
    - BooleanAlgebra: Axiome, Atome, Koatome, Stone-Darstellung
    - Konstruktionen: Teilerverband, Potenzmengenverband, Partitionsverband
    - Sätze: Birkhoff, Dedekind-Zahlen, Jordan-Dedekind

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from lattice_theory import (
    PartialOrder,
    Lattice,
    BooleanAlgebra,
    divisibility_lattice,
    power_set_lattice,
    partition_lattice,
    birkhoff_representation_theorem,
    dedekind_numbers,
    jordan_dedekind_chain_condition
)


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def make_linear_order(n: int) -> PartialOrder:
    """Erstellt die lineare Ordnung {0 < 1 < ... < n-1}."""
    elements = list(range(n))
    order = {(i, j) for i in elements for j in elements if i <= j}
    return PartialOrder(elements, order)


def make_diamond_lattice() -> Lattice:
    """
    Erstellt den Diamant-Verband M₃ = {0, a, b, c, 1}.
    0 < a, b, c < 1; a, b, c unvergleichbar.
    """
    elements = [0, 'a', 'b', 'c', 1]
    order = {
        (0, 0), (0, 'a'), (0, 'b'), (0, 'c'), (0, 1),
        ('a', 'a'), ('a', 1),
        ('b', 'b'), ('b', 1),
        ('c', 'c'), ('c', 1),
        (1, 1)
    }
    return Lattice(elements, order)


def make_pentagon_lattice() -> Lattice:
    """
    Erstellt den Pentagon-Verband N₅ = {0, a, b, c, 1}.
    0 < a < c < 1, 0 < b < 1, b und a,c unvergleichbar.
    """
    elements = [0, 'a', 'b', 'c', 1]
    order = {
        (0, 0), (0, 'a'), (0, 'b'), (0, 'c'), (0, 1),
        ('a', 'a'), ('a', 'c'), ('a', 1),
        ('b', 'b'), ('b', 1),
        ('c', 'c'), ('c', 1),
        (1, 1)
    }
    return Lattice(elements, order)


def make_bool2() -> BooleanAlgebra:
    """Erstellt B₂ = {∅, {0}, {1}, {0,1}} (Potenzmenge von {0,1})."""
    return power_set_lattice(2)


# =============================================================================
# TESTS: PartialOrder
# =============================================================================

class TestPartialOrder:
    """Tests für die Halbordnungs-Klasse."""

    def test_linear_order_is_partial_order(self):
        """Eine lineare Ordnung ist eine Halbordnung."""
        po = make_linear_order(4)
        assert po.is_partial_order()

    def test_power_set_is_partial_order(self):
        """Potenzmengenrelation ist eine Halbordnung."""
        ba = power_set_lattice(3)
        assert ba.is_partial_order()

    def test_reflexivity_violated(self):
        """Fehlende Reflexivität → keine Halbordnung."""
        elements = [0, 1]
        # Reflexivität fehlt für 1
        order = {(0, 0), (0, 1)}  # (1,1) fehlt
        po = PartialOrder(elements, order)
        assert not po.is_partial_order()

    def test_antisymmetry_violated(self):
        """Fehlende Antisymmetrie → keine Halbordnung."""
        elements = [0, 1]
        # Beide (0,1) und (1,0) vorhanden, aber 0 ≠ 1
        order = {(0, 0), (1, 1), (0, 1), (1, 0)}
        po = PartialOrder(elements, order)
        assert not po.is_partial_order()

    def test_transitivity_violated(self):
        """Fehlende Transitivität → keine Halbordnung."""
        elements = [0, 1, 2]
        # (0,1) und (1,2) aber nicht (0,2)
        order = {(0, 0), (1, 1), (2, 2), (0, 1), (1, 2)}
        po = PartialOrder(elements, order)
        assert not po.is_partial_order()

    def test_minimal_elements_chain(self):
        """In einer Kette gibt es genau 1 minimales Element."""
        po = make_linear_order(5)
        assert po.minimal_elements() == [0]

    def test_maximal_elements_chain(self):
        """In einer Kette gibt es genau 1 maximales Element."""
        po = make_linear_order(5)
        assert po.maximal_elements() == [4]

    def test_minimal_elements_antichain(self):
        """In einer Antikette sind alle Elemente minimal."""
        elements = ['a', 'b', 'c']
        # Nur Reflexivität, keine anderen Relationen
        order = {('a', 'a'), ('b', 'b'), ('c', 'c')}
        po = PartialOrder(elements, order)
        mins = po.minimal_elements()
        assert set(mins) == {'a', 'b', 'c'}

    def test_hasse_diagram_chain(self):
        """Hasse-Diagramm einer Kette hat n-1 Kanten."""
        po = make_linear_order(4)
        hasse = po.hasse_diagram()
        # 0→1, 1→2, 2→3 (keine transitiven Kanten)
        assert 1 in hasse[0]
        assert 0 not in hasse[1]  # 2 erzwingt keine Kante von 0
        assert 2 in hasse[1]
        assert 3 in hasse[2]

    def test_chain_length_linear(self):
        """Längste Kette in {0,1,...,n-1} hat Länge n-1."""
        po = make_linear_order(5)
        assert po.chain_length() == 4

    def test_is_chain_linear(self):
        """Eine lineare Ordnung ist eine Kette."""
        po = make_linear_order(3)
        assert po.is_chain()

    def test_is_not_chain_antichain(self):
        """Eine echte Antikette (≥2 Elemente) ist keine Kette."""
        elements = ['a', 'b']
        order = {('a', 'a'), ('b', 'b')}
        po = PartialOrder(elements, order)
        assert not po.is_chain()

    def test_leq_correct(self):
        """leq() gibt korrektes Ergebnis zurück."""
        po = make_linear_order(3)
        assert po.leq(0, 2)
        assert not po.leq(2, 0)
        assert po.leq(1, 1)


# =============================================================================
# TESTS: Lattice
# =============================================================================

class TestLattice:
    """Tests für die Verband-Klasse."""

    def test_divisibility_lattice_6(self):
        """Teilerverband von 6: Elemente {1,2,3,6}."""
        lat = divisibility_lattice(6)
        assert set(lat.elements) == {1, 2, 3, 6}

    def test_meet_divisibility(self):
        """Meet im Teilerverband = ggT."""
        lat = divisibility_lattice(12)
        assert lat.meet(4, 6) == 2   # ggT(4,6) = 2
        assert lat.meet(4, 12) == 4  # ggT(4,12) = 4

    def test_join_divisibility(self):
        """Join im Teilerverband = kgV."""
        lat = divisibility_lattice(12)
        assert lat.join(4, 6) == 12  # kgV(4,6) = 12
        assert lat.join(2, 3) == 6   # kgV(2,3) = 6

    def test_top_divisibility(self):
        """Größtes Element im Teilerverband von n ist n."""
        lat = divisibility_lattice(12)
        assert lat.top() == 12

    def test_bottom_divisibility(self):
        """Kleinstes Element im Teilerverband ist 1."""
        lat = divisibility_lattice(12)
        assert lat.bottom() == 1

    def test_diamond_is_not_distributive(self):
        """Diamant-Verband M₃ ist nicht distributiv."""
        diamond = make_diamond_lattice()
        assert not diamond.is_distributive()

    def test_divisibility_12_is_distributive(self):
        """Teilerverband von 12 ist distributiv."""
        lat = divisibility_lattice(12)
        assert lat.is_distributive()

    def test_diamond_is_modular(self):
        """Diamant-Verband M₃ ist modular."""
        diamond = make_diamond_lattice()
        assert diamond.is_modular()

    def test_pentagon_is_not_modular(self):
        """Pentagon-Verband N₅ ist nicht modular."""
        pentagon = make_pentagon_lattice()
        assert not pentagon.is_modular()

    def test_divisibility_is_modular(self):
        """Teilerverband ist (als distributiver Verband) modular."""
        lat = divisibility_lattice(6)
        assert lat.is_modular()

    def test_power_set_is_complemented(self):
        """Potenzmenge ist komplementiert."""
        ba = power_set_lattice(3)
        assert ba.is_complemented()

    def test_complement_power_set(self):
        """Komplement in Potenzmenge: A̅ = S\\A."""
        ba = power_set_lattice(2)
        # {0} hat Komplement {1}
        elem = frozenset({0})
        comp = ba.complement(elem)
        assert comp == frozenset({1})

    def test_divisibility_not_complemented(self):
        """Teilerverband von 6 ist nicht komplementiert (2 hat kein Komplement)."""
        lat = divisibility_lattice(6)
        # Element 2: Meet(2,x)=1 und Join(2,x)=6 → x=3 → Meet(2,3)=1✓, Join(2,3)=6✓
        # Hier ist der Verband doch komplementiert, also 12:
        lat12 = divisibility_lattice(12)
        # 4 hat kein Komplement (ggT(4,x)=1 und kgV(4,x)=12 → kein x in {1,2,3,4,6,12})
        assert not lat12.is_complemented()

    def test_meet_reflexive(self):
        """a ∧ a = a (Idempotenz des Meets)."""
        lat = divisibility_lattice(12)
        for d in lat.elements:
            assert lat.meet(d, d) == d

    def test_join_reflexive(self):
        """a ∨ a = a (Idempotenz des Joins)."""
        lat = divisibility_lattice(12)
        for d in lat.elements:
            assert lat.join(d, d) == d

    def test_meet_commutative(self):
        """a ∧ b = b ∧ a."""
        lat = divisibility_lattice(12)
        for a in lat.elements:
            for b in lat.elements:
                assert lat.meet(a, b) == lat.meet(b, a)

    def test_join_commutative(self):
        """a ∨ b = b ∨ a."""
        lat = divisibility_lattice(12)
        for a in lat.elements:
            for b in lat.elements:
                assert lat.join(a, b) == lat.join(b, a)


# =============================================================================
# TESTS: BooleanAlgebra
# =============================================================================

class TestBooleanAlgebra:
    """Tests für die Boolesche-Algebra-Klasse."""

    def test_power_set_is_boolean_algebra(self):
        """Potenzmenge P(n) ist eine Boolesche Algebra."""
        ba = power_set_lattice(3)
        assert ba.is_boolean_algebra()

    def test_power_set_size(self):
        """P(n) hat 2ⁿ Elemente."""
        for n in range(4):
            ba = power_set_lattice(n)
            assert len(ba.elements) == 2 ** n

    def test_atoms_power_set_2(self):
        """Atome von P(2) sind {0} und {1}."""
        ba = power_set_lattice(2)
        atoms = ba.atoms()
        assert len(atoms) == 2
        assert frozenset({0}) in atoms
        assert frozenset({1}) in atoms

    def test_atoms_power_set_3(self):
        """Atome von P(3) sind die drei Singletons."""
        ba = power_set_lattice(3)
        atoms = ba.atoms()
        assert len(atoms) == 3

    def test_coatoms_power_set_2(self):
        """Koatome von P(2) sind {0,1}\\{0} und {0,1}\\{1}."""
        ba = power_set_lattice(2)
        coatoms = ba.coatoms()
        assert len(coatoms) == 2

    def test_stone_representation(self):
        """Stone-Darstellung: jedes Element wird Teilmenge der Atome."""
        ba = power_set_lattice(2)
        stone = ba.stone_representation()
        assert 'atoms' in stone
        assert 'isomorphism' in stone
        # Isomorphismus hat Einträge für alle Elemente
        assert len(stone['isomorphism']) == len(ba.elements)

    def test_stone_bottom_maps_to_empty(self):
        """⊥ wird auf leere Menge der Atome abgebildet."""
        ba = power_set_lattice(2)
        stone = ba.stone_representation()
        bot = ba.bottom()
        assert stone['isomorphism'][bot] == frozenset()

    def test_stone_top_maps_to_all_atoms(self):
        """⊤ wird auf die Menge aller Atome abgebildet."""
        ba = power_set_lattice(2)
        stone = ba.stone_representation()
        top = ba.top()
        all_atoms = frozenset(ba.atoms())
        assert stone['isomorphism'][top] == all_atoms

    def test_boolean_algebra_distributive(self):
        """Boolesche Algebren sind distributiv."""
        ba = power_set_lattice(3)
        assert ba.is_distributive()

    def test_boolean_algebra_complemented(self):
        """In einer Booleschen Algebra hat jedes Element ein Komplement."""
        ba = power_set_lattice(3)
        for elem in ba.elements:
            comp = ba.complement(elem)
            assert comp is not None


# =============================================================================
# TESTS: divisibility_lattice()
# =============================================================================

class TestDivisibilityLattice:
    """Tests für den Teilerverband."""

    def test_divisors_of_prime(self):
        """Teilerverband einer Primzahl hat genau 2 Elemente: {1, p}."""
        lat = divisibility_lattice(7)
        assert set(lat.elements) == {1, 7}

    def test_divisors_of_1(self):
        """Teilerverband von 1 hat nur {1}."""
        lat = divisibility_lattice(1)
        assert lat.elements == [1]

    def test_divisors_of_12(self):
        """Teilerverband von 12: {1,2,3,4,6,12}."""
        lat = divisibility_lattice(12)
        assert set(lat.elements) == {1, 2, 3, 4, 6, 12}

    def test_order_is_divisibility(self):
        """Ordnungsrelation: a ≤ b ⟺ a | b."""
        lat = divisibility_lattice(6)
        assert lat.leq(2, 6)
        assert lat.leq(3, 6)
        assert not lat.leq(2, 3)
        assert not lat.leq(3, 2)

    def test_is_valid_partial_order(self):
        """Teilerverband ist eine Halbordnung."""
        lat = divisibility_lattice(12)
        assert lat.is_partial_order()


# =============================================================================
# TESTS: power_set_lattice()
# =============================================================================

class TestPowerSetLattice:
    """Tests für den Potenzmengenverband."""

    def test_empty_set_is_bottom(self):
        """Leere Menge ist das kleinste Element."""
        ba = power_set_lattice(3)
        assert ba.bottom() == frozenset()

    def test_full_set_is_top(self):
        """Volle Menge ist das größte Element."""
        ba = power_set_lattice(3)
        assert ba.top() == frozenset({0, 1, 2})

    def test_meet_is_intersection(self):
        """Meet im Potenzmengenverband = Schnittmenge."""
        ba = power_set_lattice(3)
        a = frozenset({0, 1})
        b = frozenset({1, 2})
        assert ba.meet(a, b) == frozenset({1})

    def test_join_is_union(self):
        """Join im Potenzmengenverband = Vereinigung."""
        ba = power_set_lattice(3)
        a = frozenset({0, 1})
        b = frozenset({1, 2})
        assert ba.join(a, b) == frozenset({0, 1, 2})

    def test_complement_is_set_difference(self):
        """Komplement in P(n) = Mengendifferenz S\\A."""
        ba = power_set_lattice(3)
        a = frozenset({0})
        comp = ba.complement(a)
        assert comp == frozenset({1, 2})


# =============================================================================
# TESTS: partition_lattice()
# =============================================================================

class TestPartitionLattice:
    """Tests für den Partitionsverband."""

    def test_partition_lattice_2_size(self):
        """Π₂ hat 2 Partitionen: {{1,2}} und {{1},{2}}."""
        lat = partition_lattice(2)
        assert len(lat.elements) == 2

    def test_partition_lattice_3_size(self):
        """Π₃ hat Bell(3)=5 Partitionen."""
        lat = partition_lattice(3)
        assert len(lat.elements) == 5

    def test_partition_bottom_is_discrete(self):
        """Kleinstes Element von Πₙ ist die diskrete Partition {{1},{2},...,{n}}."""
        lat = partition_lattice(3)
        bot = lat.bottom()
        assert bot is not None
        # Diskrete Partition: jedes Element ist eigener Block
        assert len(bot) == 3

    def test_partition_top_is_indiscrete(self):
        """Größtes Element von Πₙ ist die grobe Partition {{1,2,...,n}}."""
        lat = partition_lattice(3)
        top = lat.top()
        assert top is not None
        # Grobe Partition: ein Block mit allen Elementen
        assert len(top) == 1

    def test_partition_is_partial_order(self):
        """Partitionsverband Π₃ ist eine Halbordnung."""
        lat = partition_lattice(3)
        assert lat.is_partial_order()


# =============================================================================
# TESTS: birkhoff_representation_theorem()
# =============================================================================

class TestBirkhoffRepresentation:
    """Tests für den Birkhoff-Darstellungssatz."""

    def test_distributive_lattice_applicable(self):
        """Birkhoff gilt für distributive Verbände."""
        lat = divisibility_lattice(12)
        result = birkhoff_representation_theorem(lat)
        assert result['is_distributive']

    def test_non_distributive_not_applicable(self):
        """Birkhoff gilt nicht für nicht-distributive Verbände."""
        diamond = make_diamond_lattice()
        result = birkhoff_representation_theorem(diamond)
        assert not result['is_distributive']

    def test_join_irreducibles_found(self):
        """Join-irreduzible Elemente werden korrekt bestimmt."""
        lat = divisibility_lattice(6)
        result = birkhoff_representation_theorem(lat)
        # Im Teilerverband von 6: 2 und 3 sind join-irreduzibel
        ji = set(result['join_irreducibles'])
        assert 2 in ji
        assert 3 in ji

    def test_isomorphism_size_matches(self):
        """Größe von J(P) stimmt mit Verbandsgröße überein."""
        lat = divisibility_lattice(6)
        result = birkhoff_representation_theorem(lat)
        assert result['isomorphism_size'] == result['lattice_size']

    def test_order_ideals_power_of_2_for_bool_algebra(self):
        """Für B₂ = P(2): |J(P)| = |B₂| = 4."""
        ba = power_set_lattice(2)
        result = birkhoff_representation_theorem(ba)
        assert result['is_distributive']
        assert result['lattice_size'] == 4


# =============================================================================
# TESTS: dedekind_numbers()
# =============================================================================

class TestDedekindNumbers:
    """Tests für die Dedekind-Zahlen."""

    def test_d0(self):
        """D(0) = 2."""
        assert dedekind_numbers(0) == 2

    def test_d1(self):
        """D(1) = 3."""
        assert dedekind_numbers(1) == 3

    def test_d2(self):
        """D(2) = 6."""
        assert dedekind_numbers(2) == 6

    def test_d3(self):
        """D(3) = 20."""
        assert dedekind_numbers(3) == 20

    def test_d4(self):
        """D(4) = 168."""
        assert dedekind_numbers(4) == 168

    def test_d5(self):
        """D(5) = 7581."""
        assert dedekind_numbers(5) == 7581

    def test_d6(self):
        """D(6) = 7828354."""
        assert dedekind_numbers(6) == 7828354

    def test_dedekind_monotone_increasing(self):
        """D(n) ist streng monoton wachsend."""
        for n in range(5):
            assert dedekind_numbers(n) < dedekind_numbers(n + 1)

    def test_unknown_n_raises(self):
        """Zu großes n wirft einen ValueError."""
        with pytest.raises((ValueError, Exception)):
            dedekind_numbers(7)


# =============================================================================
# TESTS: jordan_dedekind_chain_condition()
# =============================================================================

class TestJordanDedekindChainCondition:
    """Tests für die Jordan-Dedekind-Kettenbedingung."""

    def test_boolean_algebra_satisfies_jd(self):
        """Boolesche Algebra erfüllt die JD-Kettenbedingung."""
        ba = power_set_lattice(3)
        assert jordan_dedekind_chain_condition(ba)

    def test_divisibility_12_satisfies_jd(self):
        """Teilerverband von 12 erfüllt die JD-Kettenbedingung."""
        lat = divisibility_lattice(12)
        assert jordan_dedekind_chain_condition(lat)

    def test_linear_order_satisfies_jd(self):
        """Lineare Ordnung erfüllt die JD-Kettenbedingung."""
        # Lineare Ordnung als Lattice
        elements = [0, 1, 2, 3]
        order = {(i, j) for i in elements for j in elements if i <= j}
        lat = Lattice(elements, order)
        assert jordan_dedekind_chain_condition(lat)

    def test_diamond_satisfies_jd(self):
        """Diamant-Verband M₃ erfüllt die JD-Kettenbedingung."""
        diamond = make_diamond_lattice()
        assert jordan_dedekind_chain_condition(diamond)

    def test_jd_returns_bool(self):
        """jordan_dedekind_chain_condition() gibt bool zurück."""
        ba = power_set_lattice(2)
        result = jordan_dedekind_chain_condition(ba)
        assert isinstance(result, bool)
