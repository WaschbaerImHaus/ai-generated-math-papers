"""
@file test_algebraic_structures.py
@brief Tests für das algebraic_structures-Modul.
       Prüft Magma, Halbgruppe, Monoid, Gruppe und freie Strukturen.
@description
    Testet die Klassen und Funktionen aus algebraic_structures.py:
    - Magma: Cayley-Tabelle, Abgeschlossenheit, Kommutativität
    - Semigroup: Assoziativität, Halbgruppen-Prüfung
    - Monoid: Neutrales Element, Potenzen
    - GroupFromMagma: Inverse, Gruppen-Prüfung
    - Freie Strukturen: free_monoid, free_semigroup
    - Hierarchie: algebraic_structure_hierarchy
    - Beispiele: string_monoid, matrix_monoid, transformation_monoid

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
import sys
import os

# Projektpfad einfügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algebraic_structures import (
    Magma, Semigroup, Monoid, GroupFromMagma,
    free_monoid, free_semigroup, algebraic_structure_hierarchy,
    string_monoid, matrix_monoid, transformation_monoid
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen für Test-Algebren
# ---------------------------------------------------------------------------

def make_z3_addition():
    """Erstellt ℤ/3ℤ unter Addition (Gruppe der Ordnung 3)."""
    elements = [0, 1, 2]
    operation = lambda a, b: (a + b) % 3
    return elements, operation


def make_z4_addition():
    """Erstellt ℤ/4ℤ unter Addition (Gruppe der Ordnung 4)."""
    elements = [0, 1, 2, 3]
    operation = lambda a, b: (a + b) % 4
    return elements, operation


def make_nonassoc_magma():
    """Erstellt ein nicht-assoziatives Magma mit 2 Elementen."""
    elements = [0, 1]
    # a*b = 1 immer (nicht assoziativ: (0*0)*0=1*0=1, aber 0*(0*0)=0*1=1 hier doch assoziativ!)
    # Besseres Beispiel: a*b = a für a != b, a*a = 1-a
    table = [[1, 0], [0, 1]]  # a*b = 1-a wenn a=b, sonst a
    # Tatsächlich: a*b definiert durch table
    operation = lambda a, b: table[a][b]
    return elements, operation


# ---------------------------------------------------------------------------
# Tests für Magma
# ---------------------------------------------------------------------------

class TestMagma:
    """Tests für die Magma-Klasse."""

    def test_magma_creation(self):
        """Test: Magma kann erstellt werden."""
        elements, op = make_z3_addition()
        m = Magma(elements, op)
        assert len(m.elements) == 3

    def test_operation_table_z3(self):
        """Test: Cayley-Tabelle von ℤ/3ℤ hat korrekte Einträge."""
        elements, op = make_z3_addition()
        m = Magma(elements, op)
        table = m.operation_table()
        assert len(table) == 3
        assert table[0] == [0, 1, 2]  # 0+0=0, 0+1=1, 0+2=2
        assert table[1] == [1, 2, 0]  # 1+0=1, 1+1=2, 1+2=0
        assert table[2] == [2, 0, 1]  # 2+0=2, 2+1=0, 2+2=1

    def test_is_closed_z3(self):
        """Test: ℤ/3ℤ unter Addition ist abgeschlossen."""
        elements, op = make_z3_addition()
        m = Magma(elements, op)
        assert m.is_closed() is True

    def test_is_closed_not_closed(self):
        """Test: Nicht-abgeschlossene Menge wird erkannt."""
        # Operation gibt Element außerhalb der Menge zurück
        elements = [1, 2]
        op = lambda a, b: a + b  # 1+2=3, nicht in {1,2}
        m = Magma(elements, op)
        assert m.is_closed() is False

    def test_is_commutative_z3(self):
        """Test: ℤ/3ℤ unter Addition ist kommutativ."""
        elements, op = make_z3_addition()
        m = Magma(elements, op)
        assert m.is_commutative() is True

    def test_is_commutative_noncommutative(self):
        """Test: Nicht-kommutative Operation wird erkannt."""
        elements = [0, 1]
        op = lambda a, b: 1 - b  # nicht kommutativ: 0*1=0, aber 1*0=1
        m = Magma(elements, op)
        assert m.is_commutative() is False

    def test_op_direct_call(self):
        """Test: Direkte Operation über op()-Methode."""
        elements, op = make_z3_addition()
        m = Magma(elements, op)
        assert m.op(1, 2) == 0
        assert m.op(2, 2) == 1


# ---------------------------------------------------------------------------
# Tests für Semigroup
# ---------------------------------------------------------------------------

class TestSemigroup:
    """Tests für die Semigroup-Klasse."""

    def test_semigroup_is_associative_z3(self):
        """Test: ℤ/3ℤ ist assoziativ."""
        elements, op = make_z3_addition()
        sg = Semigroup(elements, op)
        assert sg.is_associative() is True

    def test_semigroup_is_semigroup_z3(self):
        """Test: ℤ/3ℤ ist eine Halbgruppe."""
        elements, op = make_z3_addition()
        sg = Semigroup(elements, op)
        assert sg.is_semigroup() is True

    def test_semigroup_not_associative(self):
        """Test: Nicht-assoziative Operation wird erkannt."""
        # Rock-Paper-Scissors: nicht assoziativ
        elements = ['rock', 'paper', 'scissors']
        wins = {
            ('rock', 'scissors'): 'rock', ('scissors', 'rock'): 'rock',
            ('paper', 'rock'): 'paper', ('rock', 'paper'): 'paper',
            ('scissors', 'paper'): 'scissors', ('paper', 'scissors'): 'scissors',
            ('rock', 'rock'): 'rock', ('paper', 'paper'): 'paper',
            ('scissors', 'scissors'): 'scissors',
        }
        op = lambda a, b: wins[(a, b)]
        sg = Semigroup(elements, op)
        # Rock-Paper-Scissors ist nicht assoziativ
        assert sg.is_associative() is False

    def test_semigroup_inherits_magma(self):
        """Test: Semigroup erbt Magma-Methoden."""
        elements, op = make_z3_addition()
        sg = Semigroup(elements, op)
        # Sollte Magma-Methoden haben
        assert hasattr(sg, 'operation_table')
        assert hasattr(sg, 'is_closed')
        assert hasattr(sg, 'is_commutative')


# ---------------------------------------------------------------------------
# Tests für Monoid
# ---------------------------------------------------------------------------

class TestMonoid:
    """Tests für die Monoid-Klasse."""

    def test_monoid_identity_z3(self):
        """Test: Neutrales Element von ℤ/3ℤ ist 0."""
        elements, op = make_z3_addition()
        m = Monoid(elements, op)
        assert m.identity_element() == 0

    def test_monoid_has_identity_z3(self):
        """Test: ℤ/3ℤ hat ein neutrales Element."""
        elements, op = make_z3_addition()
        m = Monoid(elements, op)
        assert m.has_identity() is True

    def test_monoid_no_identity(self):
        """Test: Struktur ohne neutrales Element wird erkannt."""
        # Positiv ganze Zahlen {1,2,3} unter Multiplikation haben kein Neutrales in {1,2,3}
        # außer wenn wir 1 einschließen
        elements = [2, 3, 4]
        op = lambda a, b: a * b  # Ergebnis oft außerhalb — kein Neutrales in {2,3,4}
        m = Monoid(elements, op)
        assert m.has_identity() is False

    def test_monoid_is_monoid_z3(self):
        """Test: ℤ/3ℤ ist ein Monoid."""
        elements, op = make_z3_addition()
        m = Monoid(elements, op)
        assert m.is_monoid() is True

    def test_monoid_powers_z3(self):
        """Test: Potenzen in ℤ/3ℤ (wiederholte Addition)."""
        elements, op = make_z3_addition()
        m = Monoid(elements, op)
        # Potenzen von 1: 0, 1, 2, 0, 1, ...
        powers = m.powers(1, 4)
        assert powers[0] == 0   # 1⁰ = e = 0
        assert powers[1] == 1   # 1¹ = 1
        assert powers[2] == 2   # 1² = 2
        assert powers[3] == 0   # 1³ = 0 (mod 3)
        assert powers[4] == 1   # 1⁴ = 1 (mod 3)

    def test_monoid_powers_length(self):
        """Test: powers() gibt Liste mit n+1 Elementen zurück."""
        elements, op = make_z3_addition()
        m = Monoid(elements, op)
        powers = m.powers(2, 5)
        assert len(powers) == 6  # e, a, a², ..., a⁵

    def test_monoid_powers_no_identity_raises(self):
        """Test: powers() ohne Neutrales Element wirft ValueError."""
        elements = [2, 4]
        op = lambda a, b: a * b % 6
        m = Monoid(elements, op)
        if not m.has_identity():
            with pytest.raises(ValueError):
                m.powers(2, 3)


# ---------------------------------------------------------------------------
# Tests für GroupFromMagma
# ---------------------------------------------------------------------------

class TestGroupFromMagma:
    """Tests für die GroupFromMagma-Klasse."""

    def test_group_inverse_z3(self):
        """Test: Inverse in ℤ/3ℤ."""
        elements, op = make_z3_addition()
        g = GroupFromMagma(elements, op)
        assert g.inverse(0) == 0  # 0 + 0 = 0
        assert g.inverse(1) == 2  # 1 + 2 = 0
        assert g.inverse(2) == 1  # 2 + 1 = 0

    def test_group_is_group_z3(self):
        """Test: ℤ/3ℤ ist eine Gruppe."""
        elements, op = make_z3_addition()
        g = GroupFromMagma(elements, op)
        assert g.is_group() is True

    def test_group_is_group_z4(self):
        """Test: ℤ/4ℤ ist eine Gruppe."""
        elements, op = make_z4_addition()
        g = GroupFromMagma(elements, op)
        assert g.is_group() is True
        assert g.inverse(3) == 1  # 3 + 1 = 4 ≡ 0 mod 4

    def test_group_not_group_without_inverses(self):
        """Test: Monoid ohne Inverse ist keine Gruppe."""
        # ℕ unter Addition: 1 hat kein Inverses in ℕ
        elements = [0, 1, 2, 3]
        # Kein Inverses: z.B. 3 hat kein Inverses in {0,1,2,3} unter +
        # Simuliere: 3+? = 0 gibt 4-3=1... aber 0+3=3 nicht 0
        # Also normaler ℤ/4ℤ hat Inverse → nehme andere Struktur
        elements = [1, 2, 3]
        op = lambda a, b: a * b % 4  # 2*2=0 nicht in Menge
        g = GroupFromMagma(elements, op)
        # Nicht abgeschlossen → kein Monoid → keine Gruppe
        assert g.is_group() is False


# ---------------------------------------------------------------------------
# Tests für freie Strukturen
# ---------------------------------------------------------------------------

class TestFreeStructures:
    """Tests für freie algebraische Strukturen."""

    def test_free_monoid_alphabet(self):
        """Test: Freier Monoid enthält Alphabet."""
        fm = free_monoid(['a', 'b'])
        assert 'a' in fm['alphabet']
        assert 'b' in fm['alphabet']

    def test_free_monoid_identity(self):
        """Test: Freier Monoid hat leeres Wort als neutrales Element."""
        fm = free_monoid(['a', 'b'])
        assert fm['identity'] == ()

    def test_free_monoid_operation(self):
        """Test: Freier Monoid-Operation ist Konkatenation."""
        fm = free_monoid(['a', 'b'])
        op = fm['operation']
        assert op(('a',), ('b',)) == ('a', 'b')
        assert op((), ('a',)) == ('a',)
        assert op(('a',), ()) == ('a',)

    def test_free_monoid_example_words(self):
        """Test: Beispielwörter enthalten Wörter bis Länge 2."""
        fm = free_monoid(['a', 'b'])
        words = fm['example_words']
        assert () in words          # Leerstring
        assert ('a',) in words      # Einzelbuchstaben
        assert ('b',) in words
        assert ('a', 'a') in words  # Wörter der Länge 2
        assert ('a', 'b') in words

    def test_free_semigroup_no_empty_word(self):
        """Test: Freie Halbgruppe enthält kein leeres Wort."""
        fs = free_semigroup(['a', 'b'])
        assert () not in fs['words']

    def test_free_semigroup_generators(self):
        """Test: Freie Halbgruppe enthält Generatoren."""
        fs = free_semigroup(['x', 'y'])
        assert ('x',) in fs['words']
        assert ('y',) in fs['words']

    def test_free_semigroup_max_length(self):
        """Test: Freie Halbgruppe respektiert max_length."""
        fs = free_semigroup(['a'], max_length=2)
        # Nur Wörter der Länge 1 und 2
        for word in fs['words']:
            assert len(word) <= 2

    def test_free_semigroup_count(self):
        """Test: Anzahl der Wörter ist korrekt."""
        # |Σ|=1, max_length=3: 1 + 1 + 1 = 3 Wörter
        fs = free_semigroup(['a'], max_length=3)
        assert fs['count'] == 3


# ---------------------------------------------------------------------------
# Tests für Hierarchie
# ---------------------------------------------------------------------------

class TestAlgebraicHierarchy:
    """Tests für die algebraische Strukturhierarchie."""

    def test_hierarchy_has_all_structures(self):
        """Test: Hierarchie enthält alle algebraischen Strukturen."""
        h = algebraic_structure_hierarchy()
        expected = ['Magma', 'Halbgruppe', 'Monoid', 'Gruppe', 'Abelsche Gruppe',
                    'Ring', 'Integritätsbereich', 'Körper']
        for name in expected:
            assert name in h, f"'{name}' fehlt in der Hierarchie"

    def test_hierarchy_magma_axioms(self):
        """Test: Magma hat Abgeschlossenheit als Axiom."""
        h = algebraic_structure_hierarchy()
        assert any('Abgeschlossenheit' in ax for ax in h['Magma']['axiome'])

    def test_hierarchy_koerper_no_extension(self):
        """Test: Körper hat keine weiteren Erweiterungen."""
        h = algebraic_structure_hierarchy()
        assert h['Körper']['erweitert_zu'] == []

    def test_hierarchy_gruppe_has_inverse_axiom(self):
        """Test: Gruppe hat Inversen-Axiom."""
        h = algebraic_structure_hierarchy()
        axioms = h['Gruppe']['axiome']
        assert any('Inverse' in ax or 'Invers' in ax for ax in axioms)


# ---------------------------------------------------------------------------
# Tests für Beispiel-Monoide
# ---------------------------------------------------------------------------

class TestExampleMonoids:
    """Tests für konkrete Monoid-Beispiele."""

    def test_string_monoid_creation(self):
        """Test: String-Monoid kann erstellt werden."""
        m = string_monoid("ab", max_length=1)
        # Elemente: '', 'a', 'b' (3 Elemente)
        assert len(m.elements) == 3
        assert '' in m.elements

    def test_string_monoid_identity(self):
        """Test: Leerstring ist neutrales Element des String-Monoids."""
        m = string_monoid("ab", max_length=2)
        assert m.identity_element() == ''

    def test_string_monoid_is_monoid(self):
        """Test: String-Monoid ist ein Monoid."""
        m = string_monoid("a", max_length=2)
        assert m.is_monoid() is True

    def test_matrix_monoid_2x2_z2(self):
        """Test: 2×2-Matrizen über ℤ/2ℤ bilden einen Monoid."""
        m = matrix_monoid(2, 'Z2')
        # 2^(2²) = 16 Matrizen
        assert len(m.elements) == 16

    def test_matrix_monoid_identity_is_unit_matrix(self):
        """Test: Einheitsmatrix ist neutrales Element des Matrizenmonoids."""
        m = matrix_monoid(2, 'Z2')
        e = m.identity_element()
        # Einheitsmatrix: ((1,0),(0,1))
        assert e == ((1, 0), (0, 1))

    def test_matrix_monoid_is_monoid(self):
        """Test: Matrizenmonoid ist ein Monoid."""
        m = matrix_monoid(2, 'Z2')
        assert m.is_monoid() is True

    def test_transformation_monoid_t2(self):
        """Test: T₂ hat 4 Elemente (2² = 4 Funktionen)."""
        m = transformation_monoid(2)
        assert len(m.elements) == 4

    def test_transformation_monoid_t3(self):
        """Test: T₃ hat 27 Elemente (3³ = 27 Funktionen)."""
        m = transformation_monoid(3)
        assert len(m.elements) == 27

    def test_transformation_monoid_identity(self):
        """Test: Identitätsabbildung ist neutrales Element von T_n."""
        m = transformation_monoid(3)
        e = m.identity_element()
        # Identität: f(i) = i für alle i ∈ {1,2,3}
        assert e == (1, 2, 3)

    def test_transformation_monoid_is_monoid(self):
        """Test: Transformationsmonoid ist ein Monoid."""
        m = transformation_monoid(2)
        assert m.is_monoid() is True

    def test_transformation_monoid_n5_raises(self):
        """Test: T_n für n>4 wirft ValueError."""
        with pytest.raises(ValueError):
            transformation_monoid(5)

    def test_matrix_monoid_unsupported_field_raises(self):
        """Test: Nicht-unterstützter Körper wirft NotImplementedError."""
        with pytest.raises(NotImplementedError):
            matrix_monoid(2, 'R')
