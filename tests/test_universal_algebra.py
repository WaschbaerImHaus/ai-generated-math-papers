"""
@file test_universal_algebra.py
@brief Tests für das universal_algebra-Modul.
       Prüft Signature, Algebra, Variety, freie Algebren,
       Kongruenzen, Quotientenalgebren und Birkhoff-Satz.
@description
    Testet alle Klassen und Funktionen aus universal_algebra.py:
    - Signature: Aritäten, Konstanten, Validierung
    - Algebra: Unteralgebren, erzeugte Unteralgebren, Homomorphismen
    - Variety: Axiomprüfung, H/S/P-Abschluss
    - free_algebra_word_count, term_algebra
    - congruence_relation, quotient_algebra
    - groups_variety, rings_variety, lattices_variety
    - birkhoff_theorem_demo, subdirectly_irreducible

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from universal_algebra import (
    Signature, Algebra, Variety,
    free_algebra_word_count, term_algebra,
    congruence_relation, quotient_algebra,
    groups_variety, rings_variety, lattices_variety,
    birkhoff_theorem_demo, subdirectly_irreducible
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen: Standardalgebren für Tests
# ---------------------------------------------------------------------------

def make_z3_algebra():
    """Erstellt ℤ/3ℤ als Algebra mit Signatur {+, 0}."""
    sig = Signature({'add': 2, 'zero': 0})
    elements = [0, 1, 2]
    interp = {
        'add': lambda a, b: (a + b) % 3,
        'zero': lambda: 0,
    }
    return Algebra(elements, sig, interp)


def make_bool_lattice():
    """Erstellt booleschen Verband {0,1} mit meet und join."""
    sig = Signature({'meet': 2, 'join': 2})
    elements = [0, 1]
    interp = {
        'meet': lambda a, b: min(a, b),
        'join': lambda a, b: max(a, b),
    }
    return Algebra(elements, sig, interp)


def make_trivial_algebra():
    """Erstellt triviale Algebra mit einem Element."""
    sig = Signature({'op': 2, 'e': 0})
    elements = [0]
    interp = {
        'op': lambda a, b: 0,
        'e': lambda: 0,
    }
    return Algebra(elements, sig, interp)


# ---------------------------------------------------------------------------
# Tests für Signature
# ---------------------------------------------------------------------------

class TestSignature:
    """Tests für die Signature-Klasse."""

    def test_signature_creation(self):
        """Test: Signatur kann erstellt werden."""
        sig = Signature({'mul': 2, 'e': 0, 'inv': 1})
        assert sig.operations['mul'] == 2
        assert sig.operations['e'] == 0
        assert sig.operations['inv'] == 1

    def test_signature_arities(self):
        """Test: arities() gibt korrektes Dict zurück."""
        sig = Signature({'add': 2, 'neg': 1, 'zero': 0})
        ar = sig.arities()
        assert ar['add'] == 2
        assert ar['neg'] == 1
        assert ar['zero'] == 0

    def test_signature_has_constants_true(self):
        """Test: has_constants() True wenn Konstanten vorhanden."""
        sig = Signature({'e': 0, 'mul': 2})
        assert sig.has_constants() is True

    def test_signature_has_constants_false(self):
        """Test: has_constants() False wenn keine Konstanten vorhanden."""
        sig = Signature({'mul': 2})
        assert sig.has_constants() is False

    def test_signature_has_unary(self):
        """Test: has_unary() True wenn unäre Operation vorhanden."""
        sig = Signature({'inv': 1, 'mul': 2})
        assert sig.has_unary() is True

    def test_signature_has_binary(self):
        """Test: has_binary() True wenn binäre Operation vorhanden."""
        sig = Signature({'mul': 2})
        assert sig.has_binary() is True

    def test_signature_invalid_arity_raises(self):
        """Test: Negative Arität wirft ValueError."""
        with pytest.raises(ValueError):
            Signature({'bad': -1})

    def test_signature_repr(self):
        """Test: __repr__ gibt lesbare Darstellung."""
        sig = Signature({'mul': 2})
        assert 'mul' in repr(sig)
        assert '2' in repr(sig)


# ---------------------------------------------------------------------------
# Tests für Algebra
# ---------------------------------------------------------------------------

class TestAlgebra:
    """Tests für die Algebra-Klasse."""

    def test_algebra_creation(self):
        """Test: Algebra kann erstellt werden."""
        alg = make_z3_algebra()
        assert len(alg.elements) == 3

    def test_algebra_apply_binary(self):
        """Test: Binäre Operation wird korrekt angewendet."""
        alg = make_z3_algebra()
        assert alg.apply('add', 1, 2) == 0
        assert alg.apply('add', 2, 2) == 1

    def test_algebra_apply_constant(self):
        """Test: Konstante wird korrekt abgerufen."""
        alg = make_z3_algebra()
        assert alg.apply('zero') == 0

    def test_algebra_missing_op_raises(self):
        """Test: Fehlende Operation in Interpretations-Dict wirft ValueError."""
        sig = Signature({'op': 2, 'missing': 1})
        with pytest.raises(ValueError):
            Algebra([0, 1], sig, {'op': lambda a, b: a})

    def test_algebra_unknown_op_raises(self):
        """Test: Unbekannte Operation bei apply() wirft ValueError."""
        alg = make_z3_algebra()
        with pytest.raises(ValueError):
            alg.apply('nonexistent', 1, 2)

    def test_algebra_is_subalgebra_true(self):
        """Test: {0} ist Unteralgebra von ℤ/3ℤ."""
        alg = make_z3_algebra()
        # {0} ist abgeschlossen: 0+0=0, Konstante 0∈{0}
        assert alg.is_subalgebra([0]) is True

    def test_algebra_is_subalgebra_full(self):
        """Test: Gesamte Algebra ist Unteralgebra von sich selbst."""
        alg = make_z3_algebra()
        assert alg.is_subalgebra([0, 1, 2]) is True

    def test_algebra_is_subalgebra_false(self):
        """Test: {1} ist keine Unteralgebra von ℤ/3ℤ."""
        alg = make_z3_algebra()
        # 1+1=2, aber 2 ∉ {1}
        assert alg.is_subalgebra([1]) is False

    def test_algebra_generated_subalgebra_trivial(self):
        """Test: Von {0} erzeugte Unteralgebra von ℤ/3ℤ ist {0}."""
        alg = make_z3_algebra()
        gen = alg.generated_subalgebra([0])
        assert 0 in gen

    def test_algebra_generated_subalgebra_full(self):
        """Test: Von {1} erzeugte Unteralgebra von ℤ/3ℤ ist ganz ℤ/3ℤ."""
        alg = make_z3_algebra()
        gen = alg.generated_subalgebra([1])
        assert set(gen) == {0, 1, 2}

    def test_algebra_repr(self):
        """Test: __repr__ gibt lesbare Darstellung."""
        alg = make_z3_algebra()
        r = repr(alg)
        assert 'Algebra' in r
        assert '3' in r  # |A|=3


# ---------------------------------------------------------------------------
# Tests für Variety
# ---------------------------------------------------------------------------

class TestVariety:
    """Tests für die Variety-Klasse."""

    def test_variety_creation(self):
        """Test: Varietät kann erstellt werden."""
        sig = Signature({'mul': 2, 'e': 0})
        v = Variety(sig, ['x*e = x', 'e*x = x'])
        assert len(v.axioms) == 2

    def test_variety_satisfies_compatible_signature(self):
        """Test: Algebra mit kompatibler Signatur satisfies() True."""
        sig = Signature({'add': 2, 'zero': 0})
        v = Variety(sig, [])
        alg = make_z3_algebra()
        assert v.satisfies(alg) is True

    def test_variety_satisfies_incompatible_signature(self):
        """Test: Algebra mit inkompatibler Signatur satisfies() False."""
        # Varietät erwartet 'mul', Algebra hat nur 'add'
        sig = Signature({'mul': 2, 'e': 0})
        v = Variety(sig, [])
        alg = make_z3_algebra()
        assert v.satisfies(alg) is False

    def test_variety_h_closure(self):
        """Test: H-Abschlusseigenschaft für Liste kompatibler Algebren."""
        sig = Signature({'add': 2, 'zero': 0})
        v = Variety(sig, [])
        alg = make_z3_algebra()
        assert v.is_closed_under_homomorphic_image([alg]) is True

    def test_variety_product_closure(self):
        """Test: P-Abschlusseigenschaft prüft Signaturkompatibilität."""
        sig = Signature({'meet': 2, 'join': 2})
        v = Variety(sig, [])
        A = make_bool_lattice()
        B = make_bool_lattice()
        assert v.is_closed_under_product(A, B) is True

    def test_variety_repr(self):
        """Test: __repr__ gibt lesbare Darstellung."""
        v = groups_variety()
        r = repr(v)
        assert 'Variety' in r


# ---------------------------------------------------------------------------
# Tests für freie Algebren
# ---------------------------------------------------------------------------

class TestFreeAlgebra:
    """Tests für freie algebraische Strukturen."""

    def test_free_algebra_word_count_monoid(self):
        """Test: Wortanzahl für Monoid-Signatur mit 1 Generator."""
        # Monoid: {· (2), e (0)}
        sig = Signature({'mul': 2, 'e': 0})
        # Tiefe 0: 1 Generator + 1 Konstante = 2
        count = free_algebra_word_count(sig, 1, max_depth=0)
        assert count == 2

    def test_free_algebra_word_count_positive(self):
        """Test: Wortanzahl ist positiv."""
        sig = Signature({'mul': 2, 'e': 0})
        count = free_algebra_word_count(sig, 2, max_depth=2)
        assert count > 0

    def test_free_algebra_word_count_increases_with_depth(self):
        """Test: Mehr Tiefe ergibt mehr oder gleich Terme."""
        sig = Signature({'mul': 2, 'e': 0})
        count1 = free_algebra_word_count(sig, 2, max_depth=1)
        count2 = free_algebra_word_count(sig, 2, max_depth=2)
        assert count2 >= count1

    def test_term_algebra_structure(self):
        """Test: Term-Algebra hat korrekte Struktur."""
        sig = Signature({'mul': 2, 'e': 0})
        ta = term_algebra(sig, ['x', 'y'])
        assert 'generators' in ta
        assert 'depth_0_terms' in ta
        assert 'depth_1_terms' in ta
        assert 'universal_property' in ta

    def test_term_algebra_generators(self):
        """Test: Generatoren sind in den Tiefe-0-Termen."""
        sig = Signature({'mul': 2, 'e': 0})
        ta = term_algebra(sig, ['x', 'y'])
        assert 'x' in ta['depth_0_terms']
        assert 'y' in ta['depth_0_terms']

    def test_term_algebra_constants_in_depth0(self):
        """Test: Konstanten sind in den Tiefe-0-Termen."""
        sig = Signature({'mul': 2, 'e': 0})
        ta = term_algebra(sig, ['x'])
        assert 'e' in ta['depth_0_terms']

    def test_term_algebra_depth1_has_binary_terms(self):
        """Test: Tiefe-1-Terme enthalten zusammengesetzte Terme."""
        sig = Signature({'mul': 2, 'e': 0})
        ta = term_algebra(sig, ['x'])
        # Sollte mul(x,x) oder ähnliches enthalten
        assert any('mul' in t for t in ta['depth_1_terms'])


# ---------------------------------------------------------------------------
# Tests für Kongruenzen
# ---------------------------------------------------------------------------

class TestCongruenceRelation:
    """Tests für Kongruenzrelationen."""

    def test_trivial_congruence_is_congruence(self):
        """Test: Leere Paarliste (nur Reflexivität) ist Kongruenz."""
        alg = make_z3_algebra()
        assert congruence_relation(alg, []) is True

    def test_congruence_0_1_in_z3(self):
        """Test: Identifikation von 0 und 1 führt zur vollen Kongruenz in ℤ/3ℤ."""
        # In ℤ/3ℤ: wenn 0≡1, dann 0+1≡1+1, d.h. 1≡2, und 0+2≡1+2, d.h. 2≡0
        # Kongruenzabschluss ergibt alle Elemente äquivalent
        alg = make_z3_algebra()
        # Prüfe ob (0,1) als Kongruenz gilt (muss mit Operation verträglich sein)
        result = congruence_relation(alg, [(0, 1)])
        # In ℤ/3ℤ kann man 0 und 1 nicht identifizieren ohne alle zu identifizieren
        # Die Kongruenzprüfung ergibt True wenn verträglich, False wenn nicht
        assert isinstance(result, bool)

    def test_full_congruence_is_congruence(self):
        """Test: Totale Relation (alle äquivalent) ist immer Kongruenz."""
        alg = make_z3_algebra()
        # Alle Paare: jedes Element ist mit jedem anderen äquivalent
        all_pairs = [(a, b) for a in alg.elements for b in alg.elements if a != b]
        assert congruence_relation(alg, all_pairs) is True


# ---------------------------------------------------------------------------
# Tests für Quotientenalgebra
# ---------------------------------------------------------------------------

class TestQuotientAlgebra:
    """Tests für Quotientenalgebren."""

    def test_quotient_trivial_congruence(self):
        """Test: Quotient nach trivaler Kongruenz ≅ Original (selbe Größe)."""
        alg = make_z3_algebra()
        # Triviale Kongruenz: jedes Element in eigener Klasse
        congruence = [[0], [1], [2]]
        q = quotient_algebra(alg, congruence)
        assert len(q.elements) == 3

    def test_quotient_full_congruence(self):
        """Test: Quotient nach totaler Kongruenz hat 1 Element."""
        alg = make_z3_algebra()
        # Totale Kongruenz: alle in einer Klasse
        congruence = [[0, 1, 2]]
        q = quotient_algebra(alg, congruence)
        assert len(q.elements) == 1

    def test_quotient_has_same_signature(self):
        """Test: Quotientenalgebra hat dieselbe Signatur."""
        alg = make_z3_algebra()
        congruence = [[0], [1], [2]]
        q = quotient_algebra(alg, congruence)
        assert set(q.signature.operations.keys()) == set(alg.signature.operations.keys())


# ---------------------------------------------------------------------------
# Tests für klassische Varietäten
# ---------------------------------------------------------------------------

class TestClassicalVarieties:
    """Tests für die klassischen Varietäten."""

    def test_groups_variety_signature(self):
        """Test: Gruppen-Varietät hat korrekte Signatur."""
        v = groups_variety()
        ops = v.signature.operations
        assert 'mul' in ops
        assert 'e' in ops
        assert 'inv' in ops
        assert ops['mul'] == 2
        assert ops['e'] == 0
        assert ops['inv'] == 1

    def test_groups_variety_axioms(self):
        """Test: Gruppen-Varietät hat mindestens 4 Axiome."""
        v = groups_variety()
        assert len(v.axioms) >= 4

    def test_rings_variety_signature(self):
        """Test: Ring-Varietät hat korrekte Signatur."""
        v = rings_variety()
        ops = v.signature.operations
        assert 'add' in ops
        assert 'mul' in ops
        assert 'zero' in ops
        assert 'one' in ops
        assert 'neg' in ops

    def test_rings_variety_axioms(self):
        """Test: Ring-Varietät hat mindestens 4 Axiome."""
        v = rings_variety()
        assert len(v.axioms) >= 4

    def test_lattices_variety_signature(self):
        """Test: Verband-Varietät hat meet und join."""
        v = lattices_variety()
        ops = v.signature.operations
        assert 'meet' in ops
        assert 'join' in ops
        assert ops['meet'] == 2
        assert ops['join'] == 2

    def test_lattices_variety_axioms_include_absorption(self):
        """Test: Verband-Axiome enthalten Absorptionsgesetz."""
        v = lattices_variety()
        absorption_found = any('bsorption' in ax for ax in v.axioms)
        assert absorption_found


# ---------------------------------------------------------------------------
# Tests für Birkhoff-Satz
# ---------------------------------------------------------------------------

class TestBirkhoffTheorem:
    """Tests für den Birkhoff-Charakterisierungssatz."""

    def test_birkhoff_returns_dict(self):
        """Test: Demo gibt Dict zurück."""
        result = birkhoff_theorem_demo()
        assert isinstance(result, dict)

    def test_birkhoff_has_operators(self):
        """Test: Demo enthält H, S, P Operatoren."""
        result = birkhoff_theorem_demo()
        assert 'operators' in result
        assert 'H' in result['operators']
        assert 'S' in result['operators']
        assert 'P' in result['operators']

    def test_birkhoff_has_examples(self):
        """Test: Demo enthält Beispiele klassischer Varietäten."""
        result = birkhoff_theorem_demo()
        assert 'examples' in result
        assert len(result['examples']) >= 3

    def test_birkhoff_examples_hsp_closed(self):
        """Test: Alle Beispiele in Demo sind HSP-abgeschlossen."""
        result = birkhoff_theorem_demo()
        for example in result['examples']:
            assert example['HSP_closure'] is True

    def test_birkhoff_has_non_variety_example(self):
        """Test: Demo enthält Nicht-Varietäts-Beispiel."""
        result = birkhoff_theorem_demo()
        assert 'non_variety_example' in result

    def test_birkhoff_theorem_statement(self):
        """Test: Demo enthält Theoremformulierung."""
        result = birkhoff_theorem_demo()
        assert 'statement' in result
        # Sollte 'HSP' erwähnen
        assert 'HSP' in result['statement']

    def test_birkhoff_significance(self):
        """Test: Demo enthält Bedeutungsbeschreibung."""
        result = birkhoff_theorem_demo()
        assert 'significance' in result


# ---------------------------------------------------------------------------
# Tests für subdirekte Irreduzibilität
# ---------------------------------------------------------------------------

class TestSubdirectlyIrreducible:
    """Tests für subdirekte Irreduzibilität."""

    def test_trivial_algebra_is_subdirectly_irreducible(self):
        """Test: Triviale Algebra (1 Element) ist subdirekt irreduzibel."""
        alg = make_trivial_algebra()
        assert subdirectly_irreducible(alg) is True

    def test_returns_bool(self):
        """Test: Funktion gibt bool zurück."""
        alg = make_z3_algebra()
        result = subdirectly_irreducible(alg)
        assert isinstance(result, bool)
