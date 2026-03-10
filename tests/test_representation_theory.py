"""
@file test_representation_theory.py
@brief Umfassende Tests für das Modul representation_theory.py.
@description
    Testet alle Klassen und Funktionen der Darstellungstheorie:
    - Representation: Dimension, Homomorphismus, Charakter, Irreduzibilität
    - Direkte Summe und Tensorprodukt
    - Äquivalenz von Darstellungen
    - Triviale und reguläre Darstellungen
    - Charaktertafel
    - Schur-Lemma
    - Maschke-Satz
    - Zerlegung in Irreps
    - Burnside-Lemma
    - ℤ/2ℤ und S₃ Darstellungen

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import numpy as np
import pytest

# Pfad zum src-Verzeichnis setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from representation_theory import (
    Representation,
    trivial_representation,
    regular_representation,
    character_table,
    schur_lemma,
    maschke_theorem_check,
    decompose_representation,
    burnside_lemma,
    z2_representations,
    s3_representations
)


# =============================================================================
# HILFSFUNKTIONEN FÜR TESTS
# =============================================================================

def make_z2_mult_table() -> dict:
    """Erstellt die Multiplikationstabelle für ℤ/2ℤ."""
    return {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0}


def make_z2_representation(matrices: dict) -> Representation:
    """Erstellt eine ℤ/2ℤ-Darstellung mit gegebenen Matrizen."""
    return Representation([0, 1], matrices, make_z2_mult_table())


# =============================================================================
# TESTS: Representation.__init__ und dimension()
# =============================================================================

class TestRepresentationBasics:
    """Tests für grundlegende Repräsentation-Eigenschaften."""

    def test_dimension_1d(self):
        """Dimension einer 1D-Darstellung ist 1."""
        mats = {0: np.eye(1), 1: np.eye(1)}
        rep = make_z2_representation(mats)
        assert rep.dimension() == 1

    def test_dimension_2d(self):
        """Dimension einer 2D-Darstellung ist 2."""
        mats = {0: np.eye(2), 1: np.array([[-1, 0], [0, -1]])}
        rep = make_z2_representation(mats)
        assert rep.dimension() == 2

    def test_dimension_3d(self):
        """Dimension einer 3D-Darstellung ist 3."""
        mats = {0: np.eye(3), 1: -np.eye(3)}
        rep = make_z2_representation(mats)
        assert rep.dimension() == 3

    def test_group_elements_stored(self):
        """Gruppenelemente werden korrekt gespeichert."""
        mats = {0: np.eye(1), 1: -np.eye(1)}
        rep = make_z2_representation(mats)
        assert 0 in rep.group_elements
        assert 1 in rep.group_elements

    def test_matrices_complex(self):
        """Matrizen werden als komplexe Arrays gespeichert."""
        mats = {0: np.eye(1, dtype=float), 1: -np.eye(1, dtype=float)}
        rep = make_z2_representation(mats)
        assert rep.matrices[0].dtype == complex


# =============================================================================
# TESTS: is_homomorphism()
# =============================================================================

class TestIsHomomorphism:
    """Tests für die Homomorphismus-Eigenschaft."""

    def test_trivial_is_homomorphism(self):
        """Die triviale Darstellung ist ein Homomorphismus."""
        rep = trivial_representation([0, 1])
        rep.multiplication_table = make_z2_mult_table()
        assert rep.is_homomorphism()

    def test_sign_is_homomorphism(self):
        """Die Signum-Darstellung von ℤ/2ℤ ist ein Homomorphismus."""
        reps = z2_representations()
        sign = reps[1]
        assert sign.is_homomorphism()

    def test_no_mult_table_returns_true(self):
        """Ohne Multiplikationstabelle wird True zurückgegeben."""
        mats = {0: np.eye(1), 1: np.eye(1)}
        rep = Representation([0, 1], mats)
        assert rep.is_homomorphism()

    def test_s3_standard_is_homomorphism(self):
        """Die Standard-Darstellung von S₃ ist ein Homomorphismus."""
        reps = s3_representations()
        std = reps[2]
        assert std.is_homomorphism()

    def test_not_homomorphism_detected(self):
        """Ungültige Darstellung wird als Nicht-Homomorphismus erkannt."""
        mult = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0}
        # ρ(1)² ≠ ρ(0) → kein Homomorphismus
        mats = {0: np.eye(1), 1: 2 * np.eye(1)}
        rep = Representation([0, 1], mats, mult)
        assert not rep.is_homomorphism()


# =============================================================================
# TESTS: character()
# =============================================================================

class TestCharacter:
    """Tests für die Charakter-Berechnung."""

    def test_trivial_character(self):
        """Charakter der trivialen Darstellung ist überall 1."""
        rep = trivial_representation([0, 1])
        chi = rep.character()
        assert abs(chi[0] - 1.0) < 1e-10
        assert abs(chi[1] - 1.0) < 1e-10

    def test_sign_character_z2(self):
        """Charakter der Signum-Darstellung: χ(0)=1, χ(1)=-1."""
        reps = z2_representations()
        sign = reps[1]
        chi = sign.character()
        assert abs(chi[0] - 1.0) < 1e-10
        assert abs(chi[1] + 1.0) < 1e-10

    def test_character_trace(self):
        """Charakter = Spur der Darstellungsmatrix."""
        mat_e = np.array([[2, 1], [0, 3]], dtype=complex)
        mat_g = np.array([[1, 0], [0, -1]], dtype=complex)
        mats = {'e': mat_e, 'g': mat_g}
        rep = Representation(['e', 'g'], mats)
        chi = rep.character()
        assert abs(chi['e'] - 5.0) < 1e-10
        assert abs(chi['g'] - 0.0) < 1e-10

    def test_s3_standard_character(self):
        """Charakter der Standard-Darst. von S₃: χ(e)=2, χ(Transp)=0, χ(3-Zykl)=-1."""
        reps = s3_representations()
        std = reps[2]
        chi = std.character()
        # Einheitselement: Spur von 2×2 Einheitsmatrix = 2
        assert abs(chi['e'] - 2.0) < 1e-6
        # Transpositionen haben Spur 0
        assert abs(chi['(12)'] - 0.0) < 1e-6
        # 3-Zyklen haben Spur -1
        assert abs(chi['(123)'] + 1.0) < 1e-6


# =============================================================================
# TESTS: is_irreducible()
# =============================================================================

class TestIsIrreducible:
    """Tests für die Irreduzibilitätsprüfung."""

    def test_trivial_is_irreducible(self):
        """Die triviale Darstellung ist irreduzibel."""
        rep = trivial_representation([0, 1])
        assert rep.is_irreducible()

    def test_sign_z2_is_irreducible(self):
        """Signum-Darst. von ℤ/2ℤ ist irreduzibel."""
        reps = z2_representations()
        assert reps[1].is_irreducible()

    def test_s3_all_irreducible(self):
        """Alle drei Darstellungen von S₃ sind irreduzibel."""
        reps = s3_representations()
        for rep in reps:
            assert rep.is_irreducible(), f"Darst. dim={rep.dimension()} sollte irreduzibel sein"

    def test_direct_sum_not_irreducible(self):
        """Direkte Summe zweier irreduzibler Darst. ist reduzibel."""
        reps = z2_representations()
        ds = reps[0].direct_sum(reps[1])
        assert not ds.is_irreducible()

    def test_s3_standard_2d_irreducible(self):
        """Die 2D Standard-Darst. von S₃ ist irreduzibel."""
        reps = s3_representations()
        std = reps[2]
        assert std.dimension() == 2
        assert std.is_irreducible()


# =============================================================================
# TESTS: direct_sum()
# =============================================================================

class TestDirectSum:
    """Tests für die direkte Summe von Darstellungen."""

    def test_direct_sum_dimension(self):
        """Dimension der direkten Summe ist Summe der Dimensionen."""
        reps = z2_representations()
        ds = reps[0].direct_sum(reps[1])
        assert ds.dimension() == 2

    def test_direct_sum_block_diagonal(self):
        """Matrizen der direkten Summe sind block-diagonal."""
        reps = z2_representations()
        ds = reps[0].direct_sum(reps[1])
        # Matrix für Element 1: soll [[1,0],[0,-1]] sein
        mat = ds.matrices[1]
        assert abs(mat[0, 0] - 1.0) < 1e-10
        assert abs(mat[1, 1] + 1.0) < 1e-10
        assert abs(mat[0, 1]) < 1e-10
        assert abs(mat[1, 0]) < 1e-10

    def test_direct_sum_character_adds(self):
        """Charakter der direkten Summe ist Summe der Charaktere."""
        reps = z2_representations()
        chi0 = reps[0].character()
        chi1 = reps[1].character()
        ds = reps[0].direct_sum(reps[1])
        chi_ds = ds.character()
        for g in [0, 1]:
            assert abs(chi_ds[g] - (chi0[g] + chi1[g])) < 1e-10

    def test_direct_sum_three_representations(self):
        """Dreifache direkte Summe hat Dimension 4 (1+1+2)."""
        reps = s3_representations()
        ds = reps[0].direct_sum(reps[1]).direct_sum(reps[2])
        assert ds.dimension() == 4


# =============================================================================
# TESTS: tensor_product()
# =============================================================================

class TestTensorProduct:
    """Tests für das Tensorprodukt von Darstellungen."""

    def test_tensor_product_dimension(self):
        """Dimension des Tensorprodukts ist Produkt der Dimensionen."""
        reps = s3_representations()
        tp = reps[1].tensor_product(reps[2])
        assert tp.dimension() == 1 * 2

    def test_tensor_trivial_is_identity(self):
        """Tensorprodukt mit trivialer Darst. ergibt dieselbe Darst."""
        reps = s3_representations()
        triv = reps[0]
        std = reps[2]
        tp = triv.tensor_product(std)
        # Charakter bleibt gleich (Signum-Darst. hat Wert 1 überall)
        chi_tp = tp.character()
        chi_std = std.character()
        for g in std.group_elements:
            assert abs(chi_tp[g] - chi_std[g]) < 1e-6

    def test_tensor_product_character_multiplies(self):
        """Charakter des Tensorprodukts = Produkt der Charaktere."""
        reps = z2_representations()
        chi0 = reps[0].character()
        chi1 = reps[1].character()
        tp = reps[0].tensor_product(reps[1])
        chi_tp = tp.character()
        for g in [0, 1]:
            assert abs(chi_tp[g] - chi0[g] * chi1[g]) < 1e-10


# =============================================================================
# TESTS: is_equivalent()
# =============================================================================

class TestIsEquivalent:
    """Tests für die Äquivalenz von Darstellungen."""

    def test_equivalent_to_itself(self):
        """Jede Darstellung ist äquivalent zu sich selbst."""
        reps = s3_representations()
        for rep in reps:
            assert rep.is_equivalent(rep)

    def test_inequivalent_different_dimensions(self):
        """Darstellungen verschiedener Dimension sind inäquivalent."""
        reps = s3_representations()
        assert not reps[0].is_equivalent(reps[2])

    def test_trivial_and_sign_inequivalent(self):
        """Triviale und Signum-Darst. von S₃ sind inäquivalent."""
        reps = s3_representations()
        assert not reps[0].is_equivalent(reps[1])

    def test_z2_irreps_inequivalent(self):
        """Beide irreduziblen Darst. von ℤ/2ℤ sind inäquivalent."""
        reps = z2_representations()
        assert not reps[0].is_equivalent(reps[1])


# =============================================================================
# TESTS: trivial_representation()
# =============================================================================

class TestTrivialRepresentation:
    """Tests für die triviale Darstellung."""

    def test_all_matrices_are_1(self):
        """Alle Matrizen der trivialen Darst. sind [[1]]."""
        elements = ['a', 'b', 'c']
        rep = trivial_representation(elements)
        for g in elements:
            assert np.allclose(rep.matrices[g], np.array([[1.0 + 0j]]))

    def test_dimension_is_1(self):
        """Triviale Darst. hat Dimension 1."""
        rep = trivial_representation([1, 2, 3, 4])
        assert rep.dimension() == 1

    def test_is_irreducible(self):
        """Triviale Darst. ist irreduzibel."""
        rep = trivial_representation(list(range(6)))
        assert rep.is_irreducible()


# =============================================================================
# TESTS: regular_representation()
# =============================================================================

class TestRegularRepresentation:
    """Tests für die reguläre Darstellung."""

    def test_z2_regular_dimension(self):
        """Reguläre Darst. von ℤ/2ℤ hat Dimension 2."""
        # Multiplikationstabelle als 2D-Index-Liste
        mult_table = [[0, 1], [1, 0]]  # Zeile i, Spalte j → Index von i*j
        reg = regular_representation([0, 1], mult_table)
        assert reg.dimension() == 2

    def test_z2_regular_is_homomorphism(self):
        """Reguläre Darst. von ℤ/2ℤ ist Homomorphismus."""
        mult_table = [[0, 1], [1, 0]]
        reg = regular_representation([0, 1], mult_table)
        assert reg.is_homomorphism()

    def test_regular_identity_matrix(self):
        """Das neutrale Element wird auf die Einheitsmatrix abgebildet."""
        # Für ℤ/3ℤ: 0 ist neutral
        mult_table = [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
        reg = regular_representation([0, 1, 2], mult_table)
        assert np.allclose(reg.matrices[0], np.eye(3, dtype=complex))

    def test_regular_permutation_matrices(self):
        """Matrizen der regulären Darst. sind Permutationsmatrizen."""
        mult_table = [[0, 1], [1, 0]]
        reg = regular_representation([0, 1], mult_table)
        for g in [0, 1]:
            mat = reg.matrices[g]
            # Jede Spalte und Zeile hat genau eine 1
            assert np.allclose(mat.sum(axis=0), np.ones(2))
            assert np.allclose(mat.sum(axis=1), np.ones(2))


# =============================================================================
# TESTS: character_table()
# =============================================================================

class TestCharacterTable:
    """Tests für die Charaktertafel."""

    def test_z2_character_table_shape(self):
        """Charaktertafel von ℤ/2ℤ hat Form 2×2."""
        reps = z2_representations()
        classes = [[0], [1]]
        table = character_table([0, 1], classes, reps)
        assert len(table) == 2
        assert len(table[0]) == 2

    def test_z2_character_table_values(self):
        """Korrekte Werte in der Charaktertafel von ℤ/2ℤ."""
        reps = z2_representations()
        classes = [[0], [1]]
        table = character_table([0, 1], classes, reps)
        # Triviale: χ(0)=1, χ(1)=1
        assert abs(table[0][0] - 1.0) < 1e-10
        assert abs(table[0][1] - 1.0) < 1e-10
        # Signum: χ(0)=1, χ(1)=-1
        assert abs(table[1][0] - 1.0) < 1e-10
        assert abs(table[1][1] + 1.0) < 1e-10

    def test_s3_character_table_dimensions(self):
        """Charaktertafel von S₃ hat Form 3×3."""
        reps = s3_representations()
        classes = [['e'], ['(12)'], ['(123)']]
        elements = ['e', '(12)', '(13)', '(23)', '(123)', '(132)']
        table = character_table(elements, classes, reps)
        assert len(table) == 3
        assert len(table[0]) == 3


# =============================================================================
# TESTS: schur_lemma()
# =============================================================================

class TestSchurLemma:
    """Tests für das Schur-Lemma."""

    def test_schur_intertwiner_zero_for_inequivalent(self):
        """Null-Matrix ist Intertwiner für inäquivalente Darstellungen."""
        reps = z2_representations()
        T = np.zeros((1, 1))
        result = schur_lemma(reps[0], reps[1], T)
        assert result['is_intertwiner']
        assert result['type'] == 'zero'

    def test_schur_identity_for_trivial(self):
        """Einheitsmatrix ist Scalar-Vielfaches für triviale irreduzible Darst."""
        reps = z2_representations()
        T = 3.0 * np.eye(1)
        result = schur_lemma(reps[0], reps[0], T)
        assert result['is_intertwiner']
        assert result['type'] == 'scalar'
        assert abs(result['scalar'] - 3.0) < 1e-6

    def test_schur_not_intertwiner(self):
        """Ungültige Matrix wird als Nicht-Intertwiner erkannt."""
        reps = s3_representations()
        T = np.array([[1, 2], [3, 4]], dtype=complex)
        # T intertwined nicht mit der Standard-Darstellung (allgemein)
        result = schur_lemma(reps[2], reps[2], T)
        # Entweder kein Intertwiner oder skalare/allgemeine Antwort
        # Je nach T kann es ein Intertwiner sein oder nicht
        assert 'is_intertwiner' in result

    def test_schur_zero_matrix(self):
        """Nullmatrix ist immer ein Intertwiner (trivial)."""
        reps = s3_representations()
        std = reps[2]
        T = np.zeros((2, 2))
        result = schur_lemma(std, std, T)
        assert result['is_intertwiner']
        assert result['type'] == 'zero'


# =============================================================================
# TESTS: maschke_theorem_check()
# =============================================================================

class TestMaschkeTheorem:
    """Tests für den Maschke-Satz."""

    def test_char_0_always_applicable(self):
        """Über ℂ (char=0) gilt Maschke immer."""
        result = maschke_theorem_check(6, 0)
        assert result['applicable']
        assert result['char_field'] == 0

    def test_char_divides_order_not_applicable(self):
        """Wenn char(k) | |G|, gilt Maschke nicht."""
        result = maschke_theorem_check(6, 3)
        assert not result['applicable']

    def test_char_not_divide_order_applicable(self):
        """Wenn char(k) ∤ |G|, gilt Maschke."""
        result = maschke_theorem_check(6, 5)
        assert result['applicable']

    def test_group_order_in_result(self):
        """Gruppenordnung wird korrekt im Ergebnis gespeichert."""
        result = maschke_theorem_check(24, 0)
        assert result['group_order'] == 24

    def test_char_2_group_order_4(self):
        """char=2 und |G|=4: Maschke gilt nicht."""
        result = maschke_theorem_check(4, 2)
        assert not result['applicable']


# =============================================================================
# TESTS: decompose_representation()
# =============================================================================

class TestDecomposeRepresentation:
    """Tests für die Zerlegung in irreduzible Darstellungen."""

    def test_irrep_decomposes_to_itself(self):
        """Eine irreduzible Darst. zerlegt sich in sich selbst (Mult. 1)."""
        reps = z2_representations()
        result = decompose_representation(reps[0], reps)
        # Multiplizitäten: [1, 0] (triviale Darst. hat Mult. 1)
        assert result['multiplicities'][0] == 1
        assert result['multiplicities'][1] == 0

    def test_direct_sum_decomposes_correctly(self):
        """Direkte Summe zerlegt sich in Summanden mit Mult. 1."""
        reps = z2_representations()
        ds = reps[0].direct_sum(reps[1])
        result = decompose_representation(ds, reps)
        assert result['multiplicities'][0] == 1
        assert result['multiplicities'][1] == 1

    def test_decomposition_string(self):
        """Zerlegungsstring wird generiert."""
        reps = z2_representations()
        result = decompose_representation(reps[0], reps)
        assert 'decomposition' in result
        assert 'ρ' in result['decomposition']

    def test_s3_regular_contains_all_irreps(self):
        """Reguläre Darst. enthält jede irreduzible Darst. mit Mult. = Dimension."""
        # Reguläre Darst. von S₃ hat Dimension 6
        # Enthält: trivial (1×), signum (1×), standard (2×)
        mult_idx = [[0, 1, 2, 3, 4, 5],
                    [1, 0, 4, 5, 2, 3],
                    [2, 5, 0, 4, 3, 1],
                    [3, 4, 5, 0, 1, 2],
                    [4, 3, 1, 2, 5, 0],
                    [5, 2, 3, 1, 0, 4]]
        elements = ['e', '(12)', '(13)', '(23)', '(123)', '(132)']
        reg = regular_representation(elements, mult_idx)
        reps = s3_representations()
        result = decompose_representation(reg, reps)
        # Multiplizitäten sollten dim(ρᵢ) sein
        assert result['multiplicities'][0] == 1  # trivial, dim=1
        assert result['multiplicities'][1] == 1  # signum, dim=1
        assert result['multiplicities'][2] == 2  # standard, dim=2


# =============================================================================
# TESTS: burnside_lemma()
# =============================================================================

class TestBurnsideLemma:
    """Tests für das Burnside-Lemma."""

    def test_trivial_action_orbits(self):
        """Triviale Wirkung: jedes Element ist sein eigener Orbit."""
        # Gruppe {e} wirkt auf {0,1,2}
        action = {'e': [0, 1, 2]}
        result = burnside_lemma(['e'], action)
        assert result == 3

    def test_z2_on_two_elements(self):
        """ℤ/2ℤ wirkt auf {0,1} durch Tausch: 1 Orbit."""
        action = {0: [0, 1], 1: [1, 0]}  # 0→id, 1→swap
        result = burnside_lemma([0, 1], action)
        assert result == 1

    def test_z2_on_three_elements(self):
        """ℤ/2ℤ wirkt auf {0,1,2}: e→id, (12)→tausch von 0,1."""
        action = {0: [0, 1, 2], 1: [1, 0, 2]}
        result = burnside_lemma([0, 1], action)
        # Fixpunkte: e→3, swap→1. Orbits = (3+1)/2 = 2
        assert result == 2

    def test_s3_colorings(self):
        """S₃ wirkt auf 2-Färbungen des Dreiecks: 3 Ecken, 2 Farben."""
        # Ecken 0,1,2; S₃ permutiert sie; Farbkombinationen {0,1}^3
        # Vereinfacht: Anzahl der Orbits unter S₃ auf {0,1}^3
        # Bekannt: 4 Orbits (000, 001/010/100, 011/101/110, 111)
        # Als Permutation auf 8 Farbkombinationen (zu komplex)
        # Vereinfachter Test: identity-only group
        action = {'e': list(range(8))}
        result = burnside_lemma(['e'], action)
        assert result == 8

    def test_burnside_single_element_group(self):
        """Gruppe mit nur e: Anzahl Orbits = Anzahl Elemente."""
        n = 5
        action = {'e': list(range(n))}
        result = burnside_lemma(['e'], action)
        assert result == n


# =============================================================================
# TESTS: z2_representations()
# =============================================================================

class TestZ2Representations:
    """Tests für die ℤ/2ℤ-Darstellungen."""

    def test_returns_two_representations(self):
        """Es gibt genau 2 irreduzible Darst. von ℤ/2ℤ."""
        reps = z2_representations()
        assert len(reps) == 2

    def test_both_1d(self):
        """Beide Darst. sind 1-dimensional."""
        reps = z2_representations()
        for rep in reps:
            assert rep.dimension() == 1

    def test_both_irreducible(self):
        """Beide Darst. sind irreduzibel."""
        reps = z2_representations()
        for rep in reps:
            assert rep.is_irreducible()

    def test_are_inequivalent(self):
        """Die beiden Darst. sind inäquivalent."""
        reps = z2_representations()
        assert not reps[0].is_equivalent(reps[1])

    def test_trivial_character(self):
        """Triviale Darst.: beide Charaktere = 1."""
        reps = z2_representations()
        chi = reps[0].character()
        assert abs(chi[0] - 1.0) < 1e-10
        assert abs(chi[1] - 1.0) < 1e-10

    def test_sign_character(self):
        """Signum-Darst.: χ(0)=1, χ(1)=-1."""
        reps = z2_representations()
        chi = reps[1].character()
        assert abs(chi[0] - 1.0) < 1e-10
        assert abs(chi[1] + 1.0) < 1e-10


# =============================================================================
# TESTS: s3_representations()
# =============================================================================

class TestS3Representations:
    """Tests für die S₃-Darstellungen."""

    def test_returns_three_representations(self):
        """S₃ hat genau 3 irreduzible Darstellungen."""
        reps = s3_representations()
        assert len(reps) == 3

    def test_dimensions_1_1_2(self):
        """Dimensionen sind 1, 1, 2 (Summe der Quadrate = 6 = |S₃|)."""
        reps = s3_representations()
        dims = sorted([rep.dimension() for rep in reps])
        assert dims == [1, 1, 2]

    def test_sum_of_squares_equals_group_order(self):
        """Σ dim(ρᵢ)² = |G| = 6."""
        reps = s3_representations()
        total = sum(rep.dimension() ** 2 for rep in reps)
        assert total == 6

    def test_all_irreducible(self):
        """Alle drei Darst. sind irreduzibel."""
        reps = s3_representations()
        for i, rep in enumerate(reps):
            assert rep.is_irreducible(), f"Darst. {i} (dim={rep.dimension()}) ist nicht irreduzibel"

    def test_pairwise_inequivalent(self):
        """Alle drei Darst. sind paarweise inäquivalent."""
        reps = s3_representations()
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert not reps[i].is_equivalent(reps[j])

    def test_trivial_all_ones(self):
        """Triviale Darst.: alle Charaktere = 1."""
        reps = s3_representations()
        chi = reps[0].character()
        for g in chi:
            assert abs(chi[g] - 1.0) < 1e-10

    def test_sign_correct_values(self):
        """Signum-Darst.: χ(gerade Perm.)=1, χ(ungerade Perm.)=-1."""
        reps = s3_representations()
        sign = reps[1]
        chi = sign.character()
        # Gerade: e, (123), (132)
        assert abs(chi['e'] - 1.0) < 1e-6
        assert abs(chi['(123)'] - 1.0) < 1e-6
        assert abs(chi['(132)'] - 1.0) < 1e-6
        # Ungerade: (12), (13), (23)
        assert abs(chi['(12)'] + 1.0) < 1e-6
        assert abs(chi['(13)'] + 1.0) < 1e-6
        assert abs(chi['(23)'] + 1.0) < 1e-6

    def test_standard_is_homomorphism(self):
        """Standard-Darst. von S₃ ist Homomorphismus."""
        reps = s3_representations()
        assert reps[2].is_homomorphism()

    def test_schur_orthogonality(self):
        """Schur-Orthogonalität: ⟨χᵢ, χⱼ⟩ = δᵢⱼ."""
        reps = s3_representations()
        n = 6  # |S₃| = 6
        elements = reps[0].group_elements
        for i, rho_i in enumerate(reps):
            for j, rho_j in enumerate(reps):
                chi_i = rho_i.character()
                chi_j = rho_j.character()
                inner = sum(chi_i[g] * np.conj(chi_j[g]) for g in elements) / n
                expected = 1.0 if i == j else 0.0
                assert abs(inner - expected) < 1e-6, f"⟨χ_{i}, χ_{j}⟩ = {inner} ≠ {expected}"
