"""
@file test_group_theory.py
@brief Umfassende Tests für das Gruppentheorie-Modul (group_theory.py).
@description
    Testfälle für alle Klassen und Funktionen der Gruppentheorie:
    - Group: Grundoperationen, Untergruppen, Normalteiler, Zentrum
    - Subgroup: Lagrange-Index, Untergruppen-Test
    - GroupHomomorphism: Kern, Bild, Injektivität, Isomorphismus
    - QuotientGroup: Faktorgruppe, Ordnung, Nebenklassen
    - PermutationGroup: Zyklen, Vorzeichen, Orbits, Stabilisatoren
    - Freie Funktionen: Lagrange, Sylow, Burnside, Klassifikation

    Mindestens 35 Tests abgedeckt.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import pytest

# Pfad zum Quellverzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from group_theory import (
    Group, Subgroup, GroupHomomorphism, QuotientGroup, PermutationGroup,
    lagrange_theorem_check, sylow_theorems, classify_abelian_group,
    is_simple_group, group_action, direct_product, semidirect_product,
    symmetric_group, alternating_group, cyclic_group, dihedral_group,
    quaternion_group
)


# =============================================================================
# HILFS-FIXTURES
# =============================================================================

@pytest.fixture
def Z6():
    """Zyklische Gruppe ℤ_6."""
    return cyclic_group(6)


@pytest.fixture
def Z3():
    """Zyklische Gruppe ℤ_3."""
    return cyclic_group(3)


@pytest.fixture
def Z2():
    """Zyklische Gruppe ℤ_2."""
    return cyclic_group(2)


@pytest.fixture
def S3():
    """Symmetrische Gruppe S_3."""
    return symmetric_group(3)


@pytest.fixture
def A3():
    """Alternierende Gruppe A_3."""
    return alternating_group(3)


@pytest.fixture
def A5():
    """Alternierende Gruppe A_5 (einfach, Ordnung 60)."""
    return alternating_group(5)


@pytest.fixture
def D4():
    """Diedergruppe D_4 (Ordnung 8)."""
    return dihedral_group(4)


@pytest.fixture
def Q8():
    """Quaternionengruppe Q_8."""
    return quaternion_group()


# =============================================================================
# TEST 1-5: Zyklische Gruppe ℤ_6
# =============================================================================

class TestCyclicGroup:
    """Tests für zyklische Gruppen."""

    def test_z6_order(self, Z6):
        """Test 1: ℤ_6 hat Ordnung 6."""
        assert Z6.order() == 6

    def test_z6_is_abelian(self, Z6):
        """Test 2: ℤ_6 ist abelsch."""
        assert Z6.is_abelian() is True

    def test_z6_is_cyclic(self, Z6):
        """Test 3: ℤ_6 ist zyklisch (wird von 1 erzeugt)."""
        assert Z6.is_cyclic() is True

    def test_z6_identity(self, Z6):
        """Test 4: Neutrales Element von ℤ_6 ist 0."""
        assert Z6.identity == 0

    def test_z6_element_orders(self, Z6):
        """Test 5: Elementordnungen in ℤ_6 sind 1,6,3,2,3,6."""
        expected_orders = {0: 1, 1: 6, 2: 3, 3: 2, 4: 3, 5: 6}
        for elem, expected in expected_orders.items():
            assert Z6.element_order(elem) == expected, \
                f"ord({elem}) sollte {expected} sein"


# =============================================================================
# TEST 6-10: Symmetrische Gruppe S_3
# =============================================================================

class TestSymmetricGroup:
    """Tests für symmetrische Gruppen."""

    def test_s3_order(self, S3):
        """Test 6: S_3 hat Ordnung 6 (= 3!)."""
        assert S3.order() == 6

    def test_s3_not_abelian(self, S3):
        """Test 7: S_3 ist nicht-abelsch."""
        assert S3.is_abelian() is False

    def test_s3_not_cyclic(self, S3):
        """Test 8: S_3 ist nicht zyklisch."""
        assert S3.is_cyclic() is False

    def test_s3_identity(self, S3):
        """Test 9: Neutrales Element von S_3 ist (0,1,2)."""
        assert S3.identity == (0, 1, 2)

    def test_s4_order(self):
        """Test 10: S_4 hat Ordnung 24 (= 4!)."""
        S4 = symmetric_group(4)
        assert S4.order() == 24


# =============================================================================
# TEST 11-13: Alternierende Gruppe A_n
# =============================================================================

class TestAlternatingGroup:
    """Tests für alternierende Gruppen."""

    def test_a3_order(self, A3):
        """Test 11: A_3 hat Ordnung 3 (= 3!/2)."""
        assert A3.order() == 3

    def test_a3_is_abelian(self, A3):
        """Test 12: A_3 ≅ ℤ_3 ist abelsch."""
        assert A3.is_abelian() is True

    def test_a5_order(self, A5):
        """Test 13: A_5 hat Ordnung 60 (= 5!/2)."""
        assert A5.order() == 60


# =============================================================================
# TEST 14-16: Diedergruppe D_n
# =============================================================================

class TestDihedralGroup:
    """Tests für Diedergruppen."""

    def test_d4_order(self, D4):
        """Test 14: D_4 hat Ordnung 8 (= 2·4)."""
        assert D4.order() == 8

    def test_d3_order(self):
        """Test 15: D_3 hat Ordnung 6 (= 2·3)."""
        D3 = dihedral_group(3)
        assert D3.order() == 6

    def test_d4_not_abelian(self, D4):
        """Test 16: D_4 ist nicht-abelsch."""
        assert D4.is_abelian() is False


# =============================================================================
# TEST 17-19: Quaternionengruppe Q_8
# =============================================================================

class TestQuaternionGroup:
    """Tests für die Quaternionengruppe Q_8."""

    def test_q8_order(self, Q8):
        """Test 17: Q_8 hat Ordnung 8."""
        assert Q8.order() == 8

    def test_q8_not_abelian(self, Q8):
        """Test 18: Q_8 ist nicht-abelsch (ij = k ≠ -k = ji)."""
        assert Q8.is_abelian() is False

    def test_q8_elements(self, Q8):
        """Test 19: Q_8 enthält genau die 8 Elemente {±1, ±i, ±j, ±k}."""
        expected = {'1', '-1', 'i', '-i', 'j', '-j', 'k', '-k'}
        assert set(Q8.elements) == expected


# =============================================================================
# TEST 20-22: Lagrange-Satz
# =============================================================================

class TestLagrangeTheorem:
    """Tests für den Satz von Lagrange."""

    def test_lagrange_s3_index2(self, S3):
        """Test 20: In S_3 hat eine Untergruppe der Ordnung 2 Index 3."""
        subs = S3.subgroups()
        # Untergruppe der Ordnung 2 finden
        H2 = next(h for h in subs if h.order() == 2)
        result = lagrange_theorem_check(S3, H2)
        assert result['lagrange_holds'] is True
        assert result['index'] == 3
        assert result['order_G'] == 6
        assert result['order_H'] == 2

    def test_lagrange_s3_index3(self, S3):
        """Test 21: In S_3 hat eine Untergruppe der Ordnung 3 Index 2."""
        subs = S3.subgroups()
        H3 = next(h for h in subs if h.order() == 3)
        result = lagrange_theorem_check(S3, H3)
        assert result['lagrange_holds'] is True
        assert result['index'] == 2

    def test_lagrange_num_cosets(self, S3):
        """Test 22: Anzahl der Nebenklassen = Index."""
        subs = S3.subgroups()
        H2 = next(h for h in subs if h.order() == 2)
        result = lagrange_theorem_check(S3, H2)
        # 3 Nebenklassen für Index 3
        assert result['num_cosets'] == 3


# =============================================================================
# TEST 23-25: Sylow-Sätze
# =============================================================================

class TestSylowTheorems:
    """Tests für die Sylow-Sätze."""

    def test_sylow_s3_p2(self, S3):
        """Test 23: In S_3 gibt es n_2 = 3 Sylow-2-Untergruppen."""
        result = sylow_theorems(S3, 2)
        assert result['n_p'] == 3

    def test_sylow_s3_p3(self, S3):
        """Test 24: In S_3 gibt es n_3 = 1 Sylow-3-Untergruppe."""
        result = sylow_theorems(S3, 3)
        assert result['n_p'] == 1

    def test_sylow_congruence_check(self, S3):
        """Test 25: Sylow-Bedingung n_p ≡ 1 (mod p) ist erfüllt."""
        result_p2 = sylow_theorems(S3, 2)
        result_p3 = sylow_theorems(S3, 3)
        # n_2 = 3 ≡ 1 (mod 2), n_3 = 1 ≡ 1 (mod 3)
        assert result_p2['checks']['congruence_mod_p'] is True
        assert result_p3['checks']['congruence_mod_p'] is True


# =============================================================================
# TEST 26-28: Klassifikation abelscher Gruppen
# =============================================================================

class TestClassifyAbelianGroup:
    """Tests für den Hauptsatz über endliche abelsche Gruppen."""

    def test_classify_n4(self):
        """Test 26: n=4 ergibt {ℤ_4, ℤ_2×ℤ_2} — 2 Isomorphieklassen."""
        result = classify_abelian_group(4)
        assert result['count'] == 2
        # Klassen überprüfen: ℤ_4 und ℤ_2×ℤ_2
        descs = result['description']
        assert 'Z_4' in descs
        assert 'Z_2 × Z_2' in descs

    def test_classify_n6(self):
        """Test 27: n=6 ergibt 1 Isomorphieklasse (6=2·3, teilerfremd → ℤ_2×ℤ_3≅ℤ_6)."""
        result = classify_abelian_group(6)
        # ℤ_6 ≅ ℤ_2 × ℤ_3 (chinesischer Restsatz): genau 1 Klasse
        assert result['count'] == 1

    def test_classify_n8(self):
        """Test 28: n=8 ergibt 3 Isomorphieklassen: ℤ_8, ℤ_4×ℤ_2, ℤ_2×ℤ_2×ℤ_2."""
        result = classify_abelian_group(8)
        assert result['count'] == 3


# =============================================================================
# TEST 29-30: Einfache Gruppen
# =============================================================================

class TestSimpleGroup:
    """Tests für einfache Gruppen."""

    def test_a5_is_simple(self, A5):
        """Test 29: A_5 ist eine einfache Gruppe (Ordnung 60)."""
        assert is_simple_group(A5) is True

    def test_z6_not_simple(self, Z6):
        """Test 30: ℤ_6 ist nicht einfach (hat echte Normalteiler)."""
        assert is_simple_group(Z6) is False


# =============================================================================
# TEST 31-33: Gruppenoperation und Burnside
# =============================================================================

class TestGroupAction:
    """Tests für Gruppenoperationen und Burnside-Lemma."""

    def test_s3_transitive_on_3_points(self, S3):
        """Test 31: S_3 operiert transitiv auf {0, 1, 2}."""
        X = [0, 1, 2]
        action = lambda g, x: g[x]  # Permutation auf Punkt anwenden
        result = group_action(S3, X, action)
        assert result['is_transitive'] is True

    def test_s3_one_orbit(self, S3):
        """Test 32: S_3 hat genau 1 Orbit auf {0, 1, 2} (transitiv)."""
        X = [0, 1, 2]
        action = lambda g, x: g[x]
        result = group_action(S3, X, action)
        assert result['num_orbits'] == 1

    def test_burnside_trivial_action(self, Z2):
        """Test 33: Triviale Operation: jedes g fixiert alle x → |X/G| = |X|."""
        X = [0, 1, 2]
        # Triviale Operation: jedes Gruppenelement fixiert alle Punkte
        action = lambda g, x: x
        result = group_action(Z2, X, action)
        # Burnside: (1/2) · (3 + 3) = 3 Orbits
        assert result['num_orbits_burnside'] == 3


# =============================================================================
# TEST 34-35: GroupHomomorphism (Vorzeichen-Homomorphismus)
# =============================================================================

class TestGroupHomomorphism:
    """Tests für Gruppenhomomorphismen."""

    def test_sign_homomorphism_s3_to_z2(self, S3):
        """Test 34: Vorzeichen-Homomorphismus sgn: S_3 → ℤ_2 ist wohldefiniert."""
        Z2 = cyclic_group(2)

        # Vorzeichen-Abbildung: sgn(σ) = 0 wenn gerade, 1 wenn ungerade
        def sign_map(sigma: tuple) -> int:
            # Vorzeichen berechnen via Transpositions-Zählung
            n = len(sigma)
            visited = [False] * n
            num_transpositions = 0
            for start in range(n):
                if not visited[start]:
                    cycle_len = 0
                    current = start
                    while not visited[current]:
                        visited[current] = True
                        current = sigma[current]
                        cycle_len += 1
                    num_transpositions += cycle_len - 1
            return num_transpositions % 2  # 0 = gerade, 1 = ungerade

        phi = GroupHomomorphism(S3, Z2, sign_map, name="sgn")
        # Kern sind die geraden Permutationen (A_3)
        kern = phi.kernel()
        assert kern.order() == 3  # |A_3| = 3

    def test_sign_homomorphism_surjective(self, S3):
        """Test 35: Vorzeichen-Homomorphismus sgn: S_3 → ℤ_2 ist surjektiv."""
        Z2 = cyclic_group(2)

        def sign_map(sigma: tuple) -> int:
            n = len(sigma)
            visited = [False] * n
            num_transpositions = 0
            for start in range(n):
                if not visited[start]:
                    cycle_len = 0
                    current = start
                    while not visited[current]:
                        visited[current] = True
                        current = sigma[current]
                        cycle_len += 1
                    num_transpositions += cycle_len - 1
            return num_transpositions % 2

        phi = GroupHomomorphism(S3, Z2, sign_map, name="sgn")
        assert phi.is_surjective() is True


# =============================================================================
# ZUSÄTZLICHE TESTS (36+): Edge Cases und weitere Klassen
# =============================================================================

class TestAdditional:
    """Zusätzliche Tests für Edge-Cases und weitere Funktionalität."""

    def test_direct_product_order(self, Z2, Z3):
        """Test 36: |ℤ_2 × ℤ_3| = 6."""
        Z2xZ3 = direct_product(Z2, Z3)
        assert Z2xZ3.order() == 6

    def test_direct_product_abelian(self, Z2, Z3):
        """Test 37: ℤ_2 × ℤ_3 ist abelsch."""
        Z2xZ3 = direct_product(Z2, Z3)
        assert Z2xZ3.is_abelian() is True

    def test_permutation_cycle_type(self):
        """Test 38: Zyklentyp einer Transposition in S_4."""
        S4 = symmetric_group(4)
        # Transposition (0 1) → σ = (1, 0, 2, 3)
        transposition = (1, 0, 2, 3)
        ct = S4.cycle_type(transposition)
        # Zyklentyp: (2, 1, 1) = ein 2-Zykel, zwei 1-Zykel
        assert 2 in ct

    def test_permutation_sign_identity(self):
        """Test 39: Vorzeichen der Identitätspermutation ist +1."""
        S3 = symmetric_group(3)
        identity = (0, 1, 2)
        assert S3.sign(identity) == 1

    def test_permutation_sign_transposition(self):
        """Test 40: Vorzeichen einer Transposition ist -1."""
        S3 = symmetric_group(3)
        # Transposition (0 1): σ(0)=1, σ(1)=0, σ(2)=2
        transposition = (1, 0, 2)
        assert S3.sign(transposition) == -1

    def test_orbit_s3(self):
        """Test 41: Orbit von 0 unter S_3 ist {0, 1, 2}."""
        S3 = symmetric_group(3)
        orbit = S3.orbit(0)
        assert orbit == {0, 1, 2}

    def test_stabilizer_s3(self):
        """Test 42: Stab_S3(2) fixiert Position 2 (isomorph zu S_2)."""
        S3 = symmetric_group(3)
        stab = S3.stabilizer(2)
        # Permutationen die 2 festhalten: (0,1,2) und (1,0,2)
        assert stab.order() == 2

    def test_cayley_table_z3(self, Z3):
        """Test 43: Cayley-Tabelle von ℤ_3 hat die richtige Gestalt."""
        import numpy as np
        table = Z3.cayley_table()
        assert table.shape == (3, 3)
        # Diagonale: 0+0=0, 1+1=2, 2+2=1 → Index 0, 2, 1
        assert table[0, 0] == 0  # 0+0 = 0 → Index 0
        assert table[1, 1] == 2  # 1+1 = 2 → Index 2

    def test_center_z6(self, Z6):
        """Test 44: Zentrum von ℤ_6 ist ℤ_6 selbst (abelsche Gruppe)."""
        center = Z6.center()
        assert center.order() == 6

    def test_center_s3(self, S3):
        """Test 45: Zentrum von S_3 ist trivial {e} (S_3 ist nicht-abelsch)."""
        center = S3.center()
        assert center.order() == 1

    def test_commutator_subgroup_abelian(self, Z6):
        """Test 46: Kommutatorgruppe von ℤ_6 ist trivial {0}."""
        comm = Z6.commutator_subgroup()
        assert comm.order() == 1

    def test_commutator_subgroup_s3(self, S3):
        """Test 47: Kommutatorgruppe von S_3 ist A_3 (Ordnung 3)."""
        comm = S3.commutator_subgroup()
        assert comm.order() == 3

    def test_quotient_group_order(self, S3):
        """Test 48: |S_3/A_3| = 2 (da |S_3|=6, |A_3|=3)."""
        subs = S3.subgroups()
        # A_3 als Normalteiler der Ordnung 3
        A3 = next(h for h in subs if h.order() == 3)
        Q = QuotientGroup(S3, A3)
        assert Q.order() == 2

    def test_subgroup_index(self, S3):
        """Test 49: Index von A_3 in S_3 ist 2."""
        subs = S3.subgroups()
        A3_sub = next(h for h in subs if h.order() == 3)
        assert A3_sub.index() == 2

    def test_z5_is_simple(self):
        """Test 50: ℤ_5 ist einfach (p=5 ist prim)."""
        Z5 = cyclic_group(5)
        assert is_simple_group(Z5) is True

    def test_q8_center(self, Q8):
        """Test 51: Zentrum von Q_8 ist {1, -1} (Ordnung 2)."""
        center = Q8.center()
        assert center.order() == 2

    def test_q8_not_cyclic(self, Q8):
        """Test 52: Q_8 ist nicht zyklisch."""
        assert Q8.is_cyclic() is False

    def test_element_power(self, Z6):
        """Test 53: 2^3 = 6 ≡ 0 (mod 6) in ℤ_6."""
        # 2+2+2 = 6 ≡ 0 (mod 6)
        result = Z6.power(2, 3)
        assert result == 0

    def test_inverse_in_z6(self, Z6):
        """Test 54: Inverses von 2 in ℤ_6 ist 4 (da 2+4=6≡0)."""
        assert Z6.inverse(2) == 4

    def test_a5_simple_via_conjugacy_classes(self, A5):
        """Test 55: A_5 hat 5 Konjugationsklassen (Ordnungen 1,15,20,12,12)."""
        from group_theory import _conjugacy_classes
        classes = _conjugacy_classes(A5)
        # A_5 hat genau 5 Konjugationsklassen
        assert len(classes) == 5
        class_sizes = sorted([len(c) for c in classes])
        # Konjugationsklassen-Ordnungen: 1, 15, 20, 12, 12
        assert class_sizes == sorted([1, 15, 20, 12, 12])
