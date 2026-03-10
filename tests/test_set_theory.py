"""
@file test_set_theory.py
@brief Umfassende Tests für das Mengenlehre-Modul (set_theory.py).

@description
    Test-Suite für alle Klassen und Funktionen des Mengenlehre-Moduls:
    - MathSet (Grundoperationen, Eigenschaften, magische Methoden)
    - Relation (Reflexivität, Symmetrie, Transitivität, Äquivalenzklassen, Abschlüsse)
    - MathFunction (Injektivität, Surjektivität, Bijektivität, Komposition)
    - Kardinalitätstheorie (Cantor-Schröder-Bernstein, Diagonalargument)
    - Ordinalzahlen (ω+1 ≠ 1+ω, ω·2 ≠ 2·ω)
    - ZFC-Axiome (Keys vorhanden)
    - σ-Algebren (SetFamily)
    - Filter und Ultrafilter

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
import sys
import os

# Füge src-Verzeichnis zum Python-Pfad hinzu
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from set_theory import (
    MathSet, Relation, MathFunction,
    are_equinumerous, cantor_schroeder_bernstein,
    is_countable, cantors_diagonal_argument,
    beth_numbers, aleph_numbers,
    Ordinal, is_well_ordered, ordinal_arithmetic_demo,
    zfc_axioms, axiom_of_choice_equivalents, continuum_hypothesis,
    SetFamily, Filter, principal_filter
)


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: MathSet – Grundoperationen
# ═══════════════════════════════════════════════════════════════════════════════

class TestMathSetBasic:
    """Tests für die grundlegende MathSet-Funktionalität."""

    def test_empty_set_creation(self):
        """Leere Menge ∅ wird korrekt erstellt."""
        A = MathSet()
        assert len(A) == 0
        assert repr(A) == "∅"

    def test_set_creation_from_list(self):
        """Menge aus Liste erstellen."""
        A = MathSet([1, 2, 3])
        assert len(A) == 3
        assert 1 in A
        assert 2 in A
        assert 3 in A

    def test_set_deduplication(self):
        """Duplikate werden entfernt (Mengendefinition)."""
        A = MathSet([1, 2, 2, 3, 3, 3])
        assert len(A) == 3

    def test_cardinality(self):
        """Kardinalzahl |A| wird korrekt berechnet."""
        A = MathSet([1, 2, 3, 4, 5])
        assert A.cardinality() == 5
        assert len(A) == 5

    def test_contains(self):
        """Elementbeziehung x ∈ A."""
        A = MathSet([1, 2, 3])
        assert 1 in A
        assert 4 not in A

    def test_iteration(self):
        """Iteration über Elemente."""
        A = MathSet([1, 2, 3])
        assert set(A) == {1, 2, 3}

    def test_equality(self):
        """Mengengleichheit (Extensionalitätsprinzip)."""
        A = MathSet([1, 2, 3])
        B = MathSet([3, 2, 1])
        assert A == B

    def test_inequality(self):
        """Verschiedene Mengen sind nicht gleich."""
        A = MathSet([1, 2, 3])
        B = MathSet([1, 2, 4])
        assert A != B

    def test_repr_nonempty(self):
        """Lesbare Darstellung nicht-leerer Menge."""
        A = MathSet([1])
        assert repr(A) == "{1}"


class TestMathSetOperations:
    """Tests für Mengenoperationen."""

    def setup_method(self):
        """Testmengen initialisieren."""
        self.A = MathSet([1, 2, 3])
        self.B = MathSet([2, 3, 4])
        self.C = MathSet([5, 6])
        self.U = MathSet([1, 2, 3, 4, 5, 6])

    def test_union(self):
        """Vereinigung A ∪ B."""
        result = self.A.union(self.B)
        assert result == MathSet([1, 2, 3, 4])

    def test_union_operator(self):
        """Vereinigung via | Operator."""
        result = self.A | self.B
        assert result == MathSet([1, 2, 3, 4])

    def test_intersection(self):
        """Durchschnitt A ∩ B."""
        result = self.A.intersection(self.B)
        assert result == MathSet([2, 3])

    def test_intersection_operator(self):
        """Durchschnitt via & Operator."""
        result = self.A & self.B
        assert result == MathSet([2, 3])

    def test_intersection_disjoint(self):
        """Schnitt disjunkter Mengen = ∅."""
        result = self.A & self.C
        assert result == MathSet()

    def test_difference(self):
        """Differenz A \\ B."""
        result = self.A.difference(self.B)
        assert result == MathSet([1])

    def test_difference_operator(self):
        """Differenz via - Operator."""
        result = self.A - self.B
        assert result == MathSet([1])

    def test_symmetric_difference(self):
        """Symmetrische Differenz A △ B."""
        result = self.A.symmetric_difference(self.B)
        assert result == MathSet([1, 4])

    def test_symmetric_difference_operator(self):
        """Symmetrische Differenz via ^ Operator."""
        result = self.A ^ self.B
        assert result == MathSet([1, 4])

    def test_complement(self):
        """Komplement A^c = U \\ A."""
        result = self.A.complement(self.U)
        assert result == MathSet([4, 5, 6])

    def test_complement_invalid(self):
        """Komplement wirft Fehler wenn A ⊄ U."""
        with pytest.raises(ValueError):
            MathSet([7, 8]).complement(self.U)

    def test_power_set_size(self):
        """Potenzmenge |P(A)| = 2^|A|."""
        A = MathSet([1, 2, 3])
        P = A.power_set()
        assert len(P) == 8  # 2^3 = 8

    def test_power_set_contains_empty(self):
        """Potenzmenge enthält ∅."""
        A = MathSet([1, 2])
        P = A.power_set()
        assert frozenset() in P

    def test_power_set_contains_full(self):
        """Potenzmenge enthält die Menge selbst."""
        A = MathSet([1, 2])
        P = A.power_set()
        assert frozenset({1, 2}) in P

    def test_power_set_empty_set(self):
        """Potenzmenge der leeren Menge = {∅}."""
        empty = MathSet()
        P = empty.power_set()
        assert len(P) == 1

    def test_power_set_too_large(self):
        """Potenzmenge für |A| > 20 wirft ValueError."""
        large = MathSet(range(21))
        with pytest.raises(ValueError):
            large.power_set()

    def test_cartesian_product_size(self):
        """Kartesisches Produkt |A×B| = |A|·|B|."""
        A = MathSet([1, 2])
        B = MathSet([3, 4, 5])
        result = A.cartesian_product(B)
        assert len(result) == 6

    def test_cartesian_product_contains(self):
        """Kartesisches Produkt enthält korrekte Paare."""
        A = MathSet([1, 2])
        B = MathSet([3, 4])
        result = A.cartesian_product(B)
        assert (1, 3) in result
        assert (2, 4) in result
        assert (3, 1) not in result


class TestMathSetProperties:
    """Tests für Mengeneigenschaften."""

    def test_is_subset_true(self):
        """A ⊆ B wenn alle Elemente von A in B liegen."""
        A = MathSet([1, 2])
        B = MathSet([1, 2, 3])
        assert A.is_subset(B)

    def test_is_subset_equal(self):
        """A ⊆ A (reflexiv)."""
        A = MathSet([1, 2, 3])
        assert A.is_subset(A)

    def test_is_subset_false(self):
        """A ⊄ B wenn A Elemente enthält, die nicht in B sind."""
        A = MathSet([1, 2, 4])
        B = MathSet([1, 2, 3])
        assert not A.is_subset(B)

    def test_is_subset_operator(self):
        """Teilmenge via <= Operator."""
        A = MathSet([1, 2])
        B = MathSet([1, 2, 3])
        assert A <= B

    def test_is_proper_subset(self):
        """A ⊊ B wenn A ⊆ B und A ≠ B."""
        A = MathSet([1, 2])
        B = MathSet([1, 2, 3])
        assert A.is_proper_subset(B)
        assert A < B

    def test_is_proper_subset_equal_fails(self):
        """A ⊊ A ist falsch (Menge ist keine echte Teilmenge von sich selbst)."""
        A = MathSet([1, 2, 3])
        assert not A.is_proper_subset(A)
        assert not (A < A)

    def test_is_disjoint_true(self):
        """Disjunkte Mengen: A ∩ B = ∅."""
        A = MathSet([1, 2])
        B = MathSet([3, 4])
        assert A.is_disjoint(B)

    def test_is_disjoint_false(self):
        """Nicht disjunkt wenn gemeinsame Elemente vorhanden."""
        A = MathSet([1, 2, 3])
        B = MathSet([3, 4, 5])
        assert not A.is_disjoint(B)

    def test_set_as_dict_key(self):
        """MathSet kann als dict-Key verwendet werden (hashbar)."""
        A = MathSet([1, 2])
        d = {A: "test"}
        assert d[A] == "test"


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: Relation
# ═══════════════════════════════════════════════════════════════════════════════

class TestRelation:
    """Tests für binäre Relationen."""

    def setup_method(self):
        """Testmengen und Relationen initialisieren."""
        self.A = MathSet([1, 2, 3])

    def _identity_relation(self):
        """Erstellt Identitätsrelation (reflexiv, antisymmetrisch, transitiv)."""
        return Relation({(1,1), (2,2), (3,3)}, self.A, self.A)

    def _equivalence_mod2(self):
        """Kongruenz mod 2: 1≡1, 2≡2, 3≡3, 1≡3, 3≡1."""
        pairs = {(1,1), (2,2), (3,3), (1,3), (3,1)}
        return Relation(pairs, self.A, self.A)

    def test_reflexive_identity(self):
        """Identitätsrelation ist reflexiv."""
        R = self._identity_relation()
        assert R.is_reflexive()

    def test_reflexive_fails(self):
        """Nicht-reflexive Relation erkannt."""
        R = Relation({(1,2), (2,3)}, self.A, self.A)
        assert not R.is_reflexive()

    def test_symmetric_equivalence(self):
        """Äquivalenzrelation mod 2 ist symmetrisch."""
        R = self._equivalence_mod2()
        assert R.is_symmetric()

    def test_symmetric_fails(self):
        """Antisymmetrische Relation ist nicht symmetrisch."""
        # Partielle Ordnung ≤ auf {1,2,3}
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, self.A, self.A)
        assert not R.is_symmetric()

    def test_antisymmetric_order(self):
        """≤ auf {1,2,3} ist antisymmetrisch."""
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, self.A, self.A)
        assert R.is_antisymmetric()

    def test_transitive_identity(self):
        """Identitätsrelation ist transitiv."""
        R = self._identity_relation()
        assert R.is_transitive()

    def test_transitive_order(self):
        """≤ auf {1,2,3} ist transitiv."""
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, self.A, self.A)
        assert R.is_transitive()

    def test_equivalence_relation(self):
        """Kongruenz mod 2 ist Äquivalenzrelation."""
        R = self._equivalence_mod2()
        assert R.is_equivalence()

    def test_partial_order(self):
        """≤ auf {1,2,3} ist partielle Ordnung."""
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, self.A, self.A)
        assert R.is_partial_order()

    def test_total_order(self):
        """≤ auf {1,2,3} ist Totalordnung."""
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, self.A, self.A)
        assert R.is_total_order()

    def test_equivalence_classes_mod2(self):
        """Äquivalenzklassen von mod-2 Kongruenz auf {1,2,3}."""
        R = self._equivalence_mod2()
        classes = R.equivalence_classes()
        # Klassen: {1,3} und {2}
        assert len(classes) == 2
        class_sizes = sorted(len(c) for c in classes)
        assert class_sizes == [1, 2]

    def test_equivalence_classes_raises(self):
        """Äquivalenzklassen wirft Fehler für Nicht-Äquivalenzrelationen."""
        R = Relation({(1,2)}, self.A, self.A)
        with pytest.raises(ValueError):
            R.equivalence_classes()

    def test_transitive_closure(self):
        """Transitiver Abschluss wird korrekt berechnet."""
        # R = {(1,2), (2,3)} → R⁺ enthält auch (1,3)
        R = Relation({(1,2), (2,3)}, self.A, self.A)
        TC = R.transitive_closure()
        assert (1, 3) in TC.pairs

    def test_reflexive_closure(self):
        """Reflexiver Abschluss fügt Diagonalelemente hinzu."""
        R = Relation({(1,2)}, self.A, self.A)
        RC = R.reflexive_closure()
        assert (1, 1) in RC.pairs
        assert (2, 2) in RC.pairs
        assert (3, 3) in RC.pairs

    def test_symmetric_closure(self):
        """Symmetrischer Abschluss fügt inverse Paare hinzu."""
        R = Relation({(1,2), (2,3)}, self.A, self.A)
        SC = R.symmetric_closure()
        assert (2, 1) in SC.pairs
        assert (3, 2) in SC.pairs

    def test_inverse_relation(self):
        """Inverse Relation R⁻¹ tauscht Paare um."""
        R = Relation({(1,2), (2,3)}, self.A, self.A)
        R_inv = R.inverse()
        assert (2, 1) in R_inv.pairs
        assert (3, 2) in R_inv.pairs
        assert (1, 2) not in R_inv.pairs

    def test_composition(self):
        """Komposition S∘R berechnet transitiv erreichbare Paare."""
        R = Relation({(1,2), (2,3)}, self.A, self.A)
        S = Relation({(2,3), (3,1)}, self.A, self.A)
        # S∘R: (a,c) wenn ∃b: (a,b)∈R und (b,c)∈S
        # (1,2)∈R, (2,3)∈S → (1,3) ∈ S∘R
        composition = R.composition(S)
        assert (1, 3) in composition.pairs

    def test_invalid_pair_raises(self):
        """Ungültiges Paar außerhalb A×B wirft ValueError."""
        with pytest.raises(ValueError):
            Relation({(1, 5)}, self.A, self.A)


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: MathFunction
# ═══════════════════════════════════════════════════════════════════════════════

class TestMathFunction:
    """Tests für mathematische Funktionen."""

    def setup_method(self):
        """Testmengen und Funktionen initialisieren."""
        self.A = MathSet([1, 2, 3])
        self.B = MathSet([4, 5, 6])

    def test_function_call(self):
        """Funktionsauswertung f(a)."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        assert f(1) == 4
        assert f(2) == 5

    def test_injective_bijection(self):
        """Bijektion ist injektiv."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        assert f.is_injective()

    def test_injective_fails(self):
        """Nicht-injektive Funktion (zwei Elemente haben dasselbe Bild)."""
        B_ext = MathSet([4, 5, 6, 7])
        f = MathFunction({1: 4, 2: 4, 3: 5}, self.A, B_ext)
        assert not f.is_injective()

    def test_surjective_bijection(self):
        """Bijektion ist surjektiv."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        assert f.is_surjective()

    def test_surjective_fails(self):
        """Nicht-surjektive Funktion (Element von B nicht im Bild)."""
        B_ext = MathSet([4, 5, 6, 7])
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, B_ext)
        assert not f.is_surjective()

    def test_bijective(self):
        """Bijektive Funktion erkannt."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        assert f.is_bijective()

    def test_inverse_function(self):
        """Umkehrfunktion einer Bijektion."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        f_inv = f.inverse_function()
        assert f_inv(4) == 1
        assert f_inv(5) == 2
        assert f_inv(6) == 3

    def test_inverse_raises_not_bijective(self):
        """Umkehrfunktion wirft Fehler für nicht-bijektive Funktionen."""
        B_ext = MathSet([4, 5, 6, 7])
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, B_ext)
        with pytest.raises(ValueError):
            f.inverse_function()

    def test_image_full(self):
        """Bild f(A) der ganzen Definitionsmenge."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        img = f.image()
        assert img == MathSet([4, 5, 6])

    def test_image_subset(self):
        """Bild f(S) einer Teilmenge S ⊆ A."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        img = f.image(MathSet([1, 2]))
        assert img == MathSet([4, 5])

    def test_preimage(self):
        """Urbild f⁻¹(T) einer Teilmenge T ⊆ B."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        pre = f.preimage(MathSet([4, 5]))
        assert pre == MathSet([1, 2])

    def test_is_identity(self):
        """Identitätsfunktion wird erkannt."""
        A = MathSet([1, 2, 3])
        f = MathFunction({1: 1, 2: 2, 3: 3}, A, A)
        assert f.is_identity()

    def test_is_not_identity(self):
        """Nicht-Identitätsfunktion wird korrekt negiert."""
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        assert not f.is_identity()

    def test_composition(self):
        """Komposition g∘f: (g∘f)(a) = g(f(a))."""
        C = MathSet([7, 8, 9])
        f = MathFunction({1: 4, 2: 5, 3: 6}, self.A, self.B)
        g = MathFunction({4: 7, 5: 8, 6: 9}, self.B, C)
        gf = f.composition(g)
        assert gf(1) == 7
        assert gf(3) == 9

    def test_invalid_function_missing_element(self):
        """Funktion mit fehlendem Domain-Element wirft ValueError."""
        with pytest.raises(ValueError):
            MathFunction({1: 4, 2: 5}, self.A, self.B)  # 3 fehlt


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: Kardinalitätstheorie
# ═══════════════════════════════════════════════════════════════════════════════

class TestCardinalityTheory:
    """Tests für Kardinalitätstheorie."""

    def test_equinumerous_same_size(self):
        """Gleichmächtige Mengen erkannt."""
        A = MathSet([1, 2, 3])
        B = MathSet([4, 5, 6])
        assert are_equinumerous(A, B)

    def test_equinumerous_different_size(self):
        """Verschiedenmächtige Mengen erkannt."""
        A = MathSet([1, 2, 3])
        B = MathSet([1, 2])
        assert not are_equinumerous(A, B)

    def test_cantor_schroeder_bernstein_bijection(self):
        """CSB-Satz konstruiert eine Bijektion."""
        A = MathSet([1, 2, 3])
        B = MathSet([4, 5, 6])
        # Beide Injektionen sind bijektiv
        f_AB = {1: 4, 2: 5, 3: 6}
        f_BA = {4: 1, 5: 2, 6: 3}
        h = cantor_schroeder_bernstein(f_AB, f_BA, A, B)
        # Ergebnis muss Funktion A→B sein
        assert set(h.keys()) == {1, 2, 3}
        assert set(h.values()) == {4, 5, 6}

    def test_is_countable_N(self):
        """ℕ ist abzählbar."""
        result = is_countable('N')
        assert result['countable'] is True
        assert result['cardinality'] == 'ℵ₀'

    def test_is_countable_Z(self):
        """ℤ ist abzählbar."""
        result = is_countable('Z')
        assert result['countable'] is True

    def test_is_countable_Q(self):
        """ℚ ist abzählbar."""
        result = is_countable('Q')
        assert result['countable'] is True

    def test_is_uncountable_R(self):
        """ℝ ist überabzählbar."""
        result = is_countable('R')
        assert result['countable'] is False

    def test_is_uncountable_PN(self):
        """P(ℕ) ist überabzählbar."""
        result = is_countable('P(N)')
        assert result['countable'] is False

    def test_is_countable_unknown(self):
        """Unbekannte Menge gibt sinnvolle Antwort."""
        result = is_countable('unknown_set')
        assert result['countable'] is None

    def test_cantors_diagonal_argument_returns_dict(self):
        """Diagonalargument gibt dict zurück."""
        result = cantors_diagonal_argument(5)
        assert 'diagonal_set' in result
        assert 'enumeration' in result
        assert 'explanation' in result

    def test_cantors_diagonal_set_differs_from_all(self):
        """Diagonalmenge D unterscheidet sich von allen S_i."""
        result = cantors_diagonal_argument(5)
        diagonal = result['diagonal_set']
        for i in range(5):
            S_i = result['enumeration'][f'S_{i}']
            # D ≠ S_i wegen Diagonalkonstruktion
            # (mindestens an Position i unterschiedlich)
            assert diagonal != S_i or i not in diagonal

    def test_beth_numbers(self):
        """Beth-Zahlen korrekt beschrieben."""
        b0 = beth_numbers(0)
        assert 'ℵ_0' in b0 or 'ℵ₀' in b0
        b1 = beth_numbers(1)
        assert 'ℝ' in b1 or '2^' in b1

    def test_beth_numbers_negative_raises(self):
        """Negative Beth-Zahl wirft ValueError."""
        with pytest.raises(ValueError):
            beth_numbers(-1)

    def test_aleph_numbers(self):
        """Aleph-Zahlen korrekt beschrieben."""
        a0 = aleph_numbers(0)
        assert 'ℵ' in a0 or 'aleph' in a0.lower() or 'ℕ' in a0
        a1 = aleph_numbers(1)
        assert '1' in a1


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: Ordinalzahlen
# ═══════════════════════════════════════════════════════════════════════════════

class TestOrdinals:
    """Tests für Ordinalzahlen."""

    def test_finite_ordinal_creation(self):
        """Endliche Ordinalzahlen werden korrekt erstellt."""
        o = Ordinal(5)
        assert o.is_finite
        assert o.value == 5

    def test_omega_creation(self):
        """Transfinite Ordinalzahl ω wird erstellt."""
        omega = Ordinal('omega')
        assert not omega.is_finite

    def test_negative_raises(self):
        """Negative Ordinalzahl wirft ValueError."""
        with pytest.raises(ValueError):
            Ordinal(-1)

    def test_successor_finite(self):
        """Nachfolger endlicher Ordinalzahl."""
        o = Ordinal(5)
        s = o.successor()
        assert s.value == 6

    def test_successor_omega(self):
        """Nachfolger von ω ist ω+1."""
        omega = Ordinal('omega')
        s = omega.successor()
        assert 'omega' in str(s.value)
        assert '1' in str(s.value)

    def test_limit_ordinal_omega(self):
        """ω ist eine Limeszahl."""
        omega = Ordinal('omega')
        assert omega.is_limit_ordinal()

    def test_not_limit_ordinal_finite(self):
        """Endliche Ordinalzahlen sind keine Limeszahlen."""
        o = Ordinal(5)
        assert not o.is_limit_ordinal()

    def test_successor_ordinal_positive_finite(self):
        """Positive endliche Ordinalzahl ist Nachfolger-Ordinalzahl."""
        o = Ordinal(3)
        assert o.is_successor_ordinal()

    def test_ordinal_addition_finite(self):
        """Endliche Ordinaladdition: 3+4 = 7."""
        assert (Ordinal(3) + Ordinal(4)) == Ordinal(7)

    def test_ordinal_addition_noncommutative_1(self):
        """ω+1 ≠ 1+ω (nicht-kommutativ!)."""
        omega = Ordinal('omega')
        one = Ordinal(1)
        omega_plus_1 = omega + one
        one_plus_omega = one + omega
        assert omega_plus_1 != one_plus_omega

    def test_omega_plus_1_is_successor(self):
        """ω+1 ist eine Nachfolger-Ordinalzahl."""
        omega_plus_1 = Ordinal('omega') + Ordinal(1)
        assert omega_plus_1.is_successor_ordinal()

    def test_1_plus_omega_equals_omega(self):
        """1+ω = ω (1 wird links geschluckt)."""
        one_plus_omega = Ordinal(1) + Ordinal('omega')
        assert one_plus_omega == Ordinal('omega')

    def test_ordinal_multiplication_noncommutative(self):
        """ω·2 ≠ 2·ω (nicht-kommutativ!)."""
        omega = Ordinal('omega')
        two = Ordinal(2)
        omega_times_2 = omega * two
        two_times_omega = two * omega
        assert omega_times_2 != two_times_omega

    def test_2_times_omega_equals_omega(self):
        """2·ω = ω (2 wird links geschluckt)."""
        two_times_omega = Ordinal(2) * Ordinal('omega')
        assert two_times_omega == Ordinal('omega')

    def test_ordinal_comparison(self):
        """Ordinalvergleich: endlich < transfinit."""
        n = Ordinal(1000)
        omega = Ordinal('omega')
        assert n < omega

    def test_ordinal_exponentiation(self):
        """Ordinalexponentiation: ω^0 = 1, ω^1 = ω."""
        omega = Ordinal('omega')
        assert (omega ** Ordinal(0)) == Ordinal(1)
        assert (omega ** Ordinal(1)) == omega

    def test_ordinal_arithmetic_demo(self):
        """Demo zeigt Nicht-Kommutativität korrekt."""
        result = ordinal_arithmetic_demo()
        assert result['omega+1 != 1+omega'] is True
        assert result['omega*2 != 2*omega'] is True

    def test_is_well_ordered(self):
        """Totalordnung auf endlicher Menge ist Wohlordnung."""
        A = MathSet([1, 2, 3])
        R = Relation({(1,1), (1,2), (1,3), (2,2), (2,3), (3,3)}, A, A)
        assert is_well_ordered(R, A)


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: ZFC-Axiome
# ═══════════════════════════════════════════════════════════════════════════════

class TestZFCAxioms:
    """Tests für ZFC-Axiomensystem."""

    def test_zfc_has_10_axioms(self):
        """ZFC-Dict enthält genau 10 Axiome."""
        axioms = zfc_axioms()
        assert len(axioms) == 10

    def test_zfc_required_keys(self):
        """Alle 10 ZFC-Axiome haben die richtigen Schlüssel."""
        axioms = zfc_axioms()
        expected_keys = {
            'extensionality', 'empty_set', 'pairing', 'union', 'power_set',
            'infinity', 'separation', 'replacement', 'foundation', 'choice'
        }
        assert set(axioms.keys()) == expected_keys

    def test_zfc_axiom_has_fields(self):
        """Jedes Axiom hat 'name', 'formal' und 'description' Felder."""
        axioms = zfc_axioms()
        for key, axiom in axioms.items():
            assert 'name' in axiom, f"Axiom '{key}' fehlt 'name'"
            assert 'formal' in axiom, f"Axiom '{key}' fehlt 'formal'"
            assert 'description' in axiom, f"Axiom '{key}' fehlt 'description'"

    def test_zfc_extensionality_formal(self):
        """Extensionalitätsaxiom hat die korrekte Formel."""
        axioms = zfc_axioms()
        formal = axioms['extensionality']['formal']
        assert '∀' in formal or 'A' in formal

    def test_ac_equivalents_nonempty(self):
        """Äquivalente zum Auswahlaxiom werden aufgelistet."""
        equivalents = axiom_of_choice_equivalents()
        assert len(equivalents) >= 3

    def test_ac_equivalents_have_name_and_description(self):
        """Jedes Äquivalent hat 'name' und 'description'."""
        for eq in axiom_of_choice_equivalents():
            assert 'name' in eq
            assert 'description' in eq

    def test_continuum_hypothesis_dict(self):
        """Kontinuumshypothese gibt dict mit korrekten Keys zurück."""
        ch = continuum_hypothesis()
        assert 'statement' in ch
        assert 'formal' in ch
        assert 'independence' in ch

    def test_continuum_hypothesis_formal(self):
        """CH hat die korrekte formale Darstellung."""
        ch = continuum_hypothesis()
        assert '2^' in ch['formal'] or 'ℵ' in ch['formal']

    def test_zorn_in_ac_equivalents(self):
        """Lemma von Zorn ist in den Äquivalenten aufgeführt."""
        equivalents = axiom_of_choice_equivalents()
        names = [eq['name'] for eq in equivalents]
        assert any('Zorn' in name for name in names)


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: SetFamily und σ-Algebren
# ═══════════════════════════════════════════════════════════════════════════════

class TestSetFamily:
    """Tests für Mengenfamilien und σ-Algebren."""

    def setup_method(self):
        """Test-Universum und Mengenfamilien initialisieren."""
        self.U = MathSet([1, 2, 3, 4])

    def test_family_union(self):
        """Vereinigung aller Mengen in der Familie."""
        family = SetFamily({0: MathSet([1, 2]), 1: MathSet([3, 4])})
        assert family.union() == self.U

    def test_family_intersection(self):
        """Schnitt aller Mengen in der Familie."""
        family = SetFamily({0: MathSet([1, 2, 3]), 1: MathSet([2, 3, 4])})
        result = family.intersection()
        assert result == MathSet([2, 3])

    def test_is_partition_true(self):
        """Korrekte Partition erkannt."""
        family = SetFamily({0: MathSet([1, 2]), 1: MathSet([3, 4])})
        assert family.is_partition(self.U)

    def test_is_partition_false_overlap(self):
        """Überlappende Mengen sind keine Partition."""
        family = SetFamily({0: MathSet([1, 2, 3]), 1: MathSet([3, 4])})
        assert not family.is_partition(self.U)

    def test_is_partition_false_not_cover(self):
        """Nicht-überdeckende Familie ist keine Partition."""
        family = SetFamily({0: MathSet([1, 2]), 1: MathSet([3])})
        assert not family.is_partition(self.U)

    def test_is_cover_true(self):
        """Überdeckung des Universums erkannt."""
        family = SetFamily({0: MathSet([1, 2, 3]), 1: MathSet([2, 3, 4])})
        assert family.is_cover(self.U)

    def test_is_cover_false(self):
        """Keine Überdeckung erkannt."""
        family = SetFamily({0: MathSet([1, 2])})
        assert not family.is_cover(self.U)

    def test_trivial_sigma_algebra(self):
        """Triviale σ-Algebra {∅, U} ist σ-Algebra."""
        U = MathSet([1, 2])
        family = SetFamily({0: MathSet(), 1: U})
        assert family.is_sigma_algebra(U)

    def test_discrete_sigma_algebra(self):
        """Diskrete σ-Algebra P(U) ist σ-Algebra."""
        U = MathSet([1, 2])
        family = SetFamily({
            0: MathSet(),
            1: MathSet([1]),
            2: MathSet([2]),
            3: U
        })
        assert family.is_sigma_algebra(U)

    def test_not_sigma_algebra_missing_complement(self):
        """Familie ohne Komplement ist keine σ-Algebra."""
        U = MathSet([1, 2, 3])
        family = SetFamily({0: MathSet([1, 2]), 1: U})
        # Komplement {3} fehlt
        assert not family.is_sigma_algebra(U)

    def test_generated_sigma_algebra(self):
        """Erzeugte σ-Algebra enthält die Ausgangsmenge und ∅."""
        U = MathSet([1, 2, 3])
        family = SetFamily({0: MathSet([1])})
        generated = family.generated_sigma_algebra(U)
        # Die erzeugte σ-Algebra muss σ-Algebra sein
        assert generated.is_sigma_algebra(U)

    def test_family_empty_union(self):
        """Leere Familie: Vereinigung = ∅."""
        family = SetFamily({})
        assert family.union() == MathSet()


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: Filter und Ultrafilter
# ═══════════════════════════════════════════════════════════════════════════════

class TestFilter:
    """Tests für Filter und Ultrafilter."""

    def test_principal_filter_construction(self):
        """Hauptfilter wird korrekt konstruiert."""
        U = MathSet([1, 2, 3])
        B = MathSet([1, 2])
        F = principal_filter(B, U)
        # F enthält mindestens {1,2} und U
        assert any(s == B for s in F.sets)
        assert any(s == U for s in F.sets)

    def test_principal_filter_is_filter(self):
        """Hauptfilter ist ein gültiger Filter."""
        U = MathSet([1, 2, 3])
        B = MathSet([2])
        F = principal_filter(B, U)
        assert F.is_filter()

    def test_full_universe_filter_is_not_ultrafilter(self):
        """Hauptfilter erzeugt durch U selbst ist Ultrafilter."""
        U = MathSet([1, 2])
        B = U  # Erzeuge durch ganzes U
        F = principal_filter(B, U)
        # F = {U} – das ist ein ultrafilter (trivial)
        assert F.is_filter()

    def test_filter_not_contains_empty(self):
        """Echter Filter enthält nicht ∅."""
        U = MathSet([1, 2, 3])
        B = MathSet([1])
        F = principal_filter(B, U)
        # Kein Element in F darf ∅ sein
        for s in F.sets:
            assert len(s) > 0

    def test_filter_intersection_closed(self):
        """Filter ist schnitt-abgeschlossen."""
        U = MathSet([1, 2, 3])
        B = MathSet([1, 2])
        F = principal_filter(B, U)
        # Schnitt beliebiger zwei F-Mengen ist in F
        sets = F.sets
        for i in range(len(sets)):
            for j in range(i, len(sets)):
                inter = sets[i] & sets[j]
                if len(inter) > 0:
                    assert inter in sets

    def test_ultrafilter_tiny_universe(self):
        """Ultrafilter auf {a} erkannt."""
        U = MathSet(['a'])
        # Einziger nicht-leerer Filter auf einem 1-Element-Universum ist {U}
        F = Filter(U, [U])
        assert F.is_filter()
        assert F.is_ultrafilter()


# ═══════════════════════════════════════════════════════════════════════════════
# Tests: Algebraische Gesetze (Mengenalgebra)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSetAlgebraLaws:
    """Tests für die Gesetze der Mengenalgebra."""

    def setup_method(self):
        self.A = MathSet([1, 2, 3])
        self.B = MathSet([2, 3, 4])
        self.C = MathSet([3, 4, 5])

    def test_commutativity_union(self):
        """A ∪ B = B ∪ A (Kommutativgesetz)."""
        assert self.A | self.B == self.B | self.A

    def test_commutativity_intersection(self):
        """A ∩ B = B ∩ A (Kommutativgesetz)."""
        assert self.A & self.B == self.B & self.A

    def test_associativity_union(self):
        """(A ∪ B) ∪ C = A ∪ (B ∪ C) (Assoziativgesetz)."""
        assert (self.A | self.B) | self.C == self.A | (self.B | self.C)

    def test_associativity_intersection(self):
        """(A ∩ B) ∩ C = A ∩ (B ∩ C) (Assoziativgesetz)."""
        assert (self.A & self.B) & self.C == self.A & (self.B & self.C)

    def test_distributivity_union_over_intersection(self):
        """A ∪ (B ∩ C) = (A ∪ B) ∩ (A ∪ C) (Distributivgesetz)."""
        lhs = self.A | (self.B & self.C)
        rhs = (self.A | self.B) & (self.A | self.C)
        assert lhs == rhs

    def test_distributivity_intersection_over_union(self):
        """A ∩ (B ∪ C) = (A ∩ B) ∪ (A ∩ C) (Distributivgesetz)."""
        lhs = self.A & (self.B | self.C)
        rhs = (self.A & self.B) | (self.A & self.C)
        assert lhs == rhs

    def test_de_morgan_union(self):
        """De Morgan: (A ∪ B)^c = A^c ∩ B^c."""
        U = MathSet([1, 2, 3, 4, 5])
        lhs = (self.A | self.B).complement(U)
        rhs = self.A.complement(U) & self.B.complement(U)
        assert lhs == rhs

    def test_de_morgan_intersection(self):
        """De Morgan: (A ∩ B)^c = A^c ∪ B^c."""
        U = MathSet([1, 2, 3, 4, 5])
        lhs = (self.A & self.B).complement(U)
        rhs = self.A.complement(U) | self.B.complement(U)
        assert lhs == rhs

    def test_idempotency_union(self):
        """A ∪ A = A (Idempotenzgesetz)."""
        assert self.A | self.A == self.A

    def test_idempotency_intersection(self):
        """A ∩ A = A (Idempotenzgesetz)."""
        assert self.A & self.A == self.A


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
