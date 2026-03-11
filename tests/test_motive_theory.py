r"""
@file test_motive_theory.py
@brief Umfassende Tests für das Modul motive_theory.py.
@description
    Testet alle Klassen und Funktionen der Motiventheorie:

    - Motive: Konstruktor, Betti-Zahlen, Euler-Charakteristik,
              Hodge-Zahlen, Tate-Twist, duales Motiv
    - ChowGroup: Konstruktor, Zykel hinzufügen, Schnittprodukt,
                 Picard-Rang
    - MotivicCohomology: Dimensionsberechnung, Beilinson-Regulator,
                         Milnor K-Theorie (endliche Körper, Zahlkörper)
    - StandardConjectures: Hard-Lefschetz, Künneth-Formel, Hodge-Vermutung
    - RealizationFunctor: Betti-, de-Rham-, l-adische Realisierung,
                          Vergleichsisomorphismen
    - motivic_l_function: Euler-Produkt für verschiedene Varietäten,
                          Konvergenzbereich

    Edge-Cases: negative Dimensionen, ungültige Typen, Nullpunkte,
    Grenzwerte, komplementäre Kodimensionen.

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import pytest

# Projektpfad einbinden
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from motive_theory import (
    Motive,
    ChowGroup,
    MotivicCohomology,
    StandardConjectures,
    RealizationFunctor,
    motivic_l_function,
    _is_prime,
    _euler_phi,
)
from exceptions import InvalidInputError


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

class TestIsPrime:
    """Tests für die interne _is_prime()-Funktion."""

    def test_small_primes(self):
        """Kleine Primzahlen werden korrekt erkannt."""
        assert _is_prime(2) is True
        assert _is_prime(3) is True
        assert _is_prime(5) is True
        assert _is_prime(7) is True
        assert _is_prime(11) is True

    def test_composites(self):
        """Zusammengesetzte Zahlen werden korrekt abgelehnt."""
        assert _is_prime(1) is False
        assert _is_prime(4) is False
        assert _is_prime(9) is False
        assert _is_prime(25) is False

    def test_zero_one_negative(self):
        """0, 1 und negative Zahlen sind keine Primzahlen."""
        assert _is_prime(0) is False
        assert _is_prime(1) is False
        assert _is_prime(-5) is False

    def test_larger_prime(self):
        """Größere Primzahl (97) korrekt erkannt."""
        assert _is_prime(97) is True


class TestEulerPhi:
    """Tests für die interne _euler_phi()-Funktion."""

    def test_phi_1(self):
        """φ(1) = 1."""
        assert _euler_phi(1) == 1

    def test_phi_prime(self):
        """φ(p) = p - 1 für Primzahl p."""
        assert _euler_phi(5) == 4
        assert _euler_phi(7) == 6

    def test_phi_prime_power(self):
        """φ(p²) = p(p-1)."""
        assert _euler_phi(4) == 2   # φ(2²) = 2
        assert _euler_phi(9) == 6   # φ(3²) = 6

    def test_phi_product(self):
        """φ(6) = 2 (φ multiplikativ: φ(2)·φ(3) = 1·2 = 2)."""
        assert _euler_phi(6) == 2


# =============================================================================
# MOTIVE-KLASSE
# =============================================================================

class TestMotive:
    """Tests für die Motive-Klasse."""

    # --- Konstruktor ---

    def test_basic_construction(self):
        """Einfache Erstellung eines Motivs ohne Fehler."""
        m = Motive('elliptic_curve', 1)
        assert m.variety_name == 'elliptic_curve'
        assert m.dimension == 1
        assert m.tate_twist == 0

    def test_dimension_zero_allowed(self):
        """Dimension 0 (Punkt) ist erlaubt."""
        m = Motive('point', 0)
        assert m.dimension == 0

    def test_negative_dimension_raises(self):
        """Negative Dimension wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            Motive('test', -1)

    def test_invalid_cohomology_type_raises(self):
        """Unbekannter Kohomologietyp wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            Motive('elliptic_curve', 1, cohomology_type='unknown')

    def test_valid_cohomology_types(self):
        """Alle gültigen Kohomologietypen funktionieren."""
        for ctype in ['betti', 'de_rham', 'l_adic', 'motivic', 'all']:
            m = Motive('p1', 1, cohomology_type=ctype)
            assert m.cohomology_type == ctype

    def test_tate_twist_stored(self):
        """Tate-Twist wird korrekt gespeichert."""
        m = Motive('elliptic_curve', 1, tate_twist=3)
        assert m.tate_twist == 3

    # --- Betti-Zahlen ---

    def test_betti_point(self):
        """Betti-Zahlen für einen Punkt: [1]."""
        m = Motive('point', 0)
        assert m.betti_numbers == [1]

    def test_betti_p1(self):
        """Betti-Zahlen für P¹: [1, 0, 1]."""
        m = Motive('p1', 1)
        assert m.betti_numbers == [1, 0, 1]

    def test_betti_p2(self):
        """Betti-Zahlen für P²: [1, 0, 1, 0, 1]."""
        m = Motive('p2', 2)
        assert m.betti_numbers == [1, 0, 1, 0, 1]

    def test_betti_elliptic_curve(self):
        """Betti-Zahlen für elliptische Kurve: [1, 2, 1]."""
        m = Motive('elliptic_curve', 1)
        assert m.betti_numbers == [1, 2, 1]

    def test_betti_k3(self):
        """Betti-Zahlen für K3-Fläche: [1, 0, 22, 0, 1]."""
        m = Motive('k3', 2)
        assert m.betti_numbers == [1, 0, 22, 0, 1]

    # --- Euler-Charakteristik ---

    def test_euler_point(self):
        """χ(Punkt) = 1."""
        m = Motive('point', 0)
        assert m.euler_characteristic() == 1

    def test_euler_p1(self):
        """χ(P¹) = 2."""
        m = Motive('p1', 1)
        assert m.euler_characteristic() == 2

    def test_euler_p2(self):
        """χ(P²) = 3."""
        m = Motive('p2', 2)
        assert m.euler_characteristic() == 3

    def test_euler_elliptic_curve(self):
        """χ(Elliptische Kurve) = 0."""
        m = Motive('elliptic_curve', 1)
        assert m.euler_characteristic() == 0

    def test_euler_k3(self):
        """χ(K3) = 24."""
        m = Motive('k3', 2)
        assert m.euler_characteristic() == 24

    # --- Hodge-Zahlen ---

    def test_hodge_elliptic_curve(self):
        """Hodge-Zahlen der elliptischen Kurve: h^{1,0}=h^{0,1}=1."""
        m = Motive('elliptic_curve', 1)
        hodge = m.hodge_numbers()
        assert hodge.get((1, 0)) == 1
        assert hodge.get((0, 1)) == 1

    def test_hodge_k3(self):
        """Hodge-Zahlen der K3-Fläche: h^{2,0}=h^{0,2}=1, h^{1,1}=20."""
        m = Motive('k3_surface', 2)
        hodge = m.hodge_numbers()
        assert hodge.get((2, 0)) == 1
        assert hodge.get((0, 2)) == 1
        assert hodge.get((1, 1)) == 20

    def test_hodge_p2_diagonal(self):
        """Hodge-Zahlen für P²: Nur diagonale Einträge h^{k,k}=1."""
        m = Motive('p2', 2)
        hodge = m.hodge_numbers()
        assert hodge.get((0, 0)) == 1
        assert hodge.get((1, 1)) == 1
        assert hodge.get((2, 2)) == 1

    def test_hodge_returns_dict(self):
        """hodge_numbers() gibt immer ein Dictionary zurück."""
        m = Motive('elliptic_curve', 1)
        hodge = m.hodge_numbers()
        assert isinstance(hodge, dict)

    # --- Tate-Twist ---

    def test_tate_twist_apply(self):
        """Tate-Twist M(n) erhöht den Twist um n."""
        m = Motive('elliptic_curve', 1, tate_twist=0)
        twisted = m.tate_twist_apply(2)
        assert twisted.tate_twist == 2
        assert twisted.dimension == m.dimension

    def test_tate_twist_negative(self):
        """Negativer Tate-Twist M(-1) ist erlaubt."""
        m = Motive('p1', 1, tate_twist=1)
        twisted = m.tate_twist_apply(-1)
        assert twisted.tate_twist == 0

    def test_tate_twist_returns_new_motive(self):
        """tate_twist_apply() gibt ein neues Motive-Objekt zurück."""
        m = Motive('elliptic_curve', 1)
        twisted = m.tate_twist_apply(1)
        assert twisted is not m

    # --- Duales Motiv ---

    def test_dual_returns_motive(self):
        """dual() gibt ein Motive-Objekt zurück."""
        m = Motive('elliptic_curve', 1)
        d = m.dual()
        assert isinstance(d, Motive)

    def test_dual_dimension_preserved(self):
        """Duales Motiv hat dieselbe Dimension."""
        m = Motive('p2', 2)
        assert m.dual().dimension == 2

    # --- repr ---

    def test_repr_contains_variety_name(self):
        """String-Repräsentation enthält den Varietätsnamen."""
        m = Motive('elliptic_curve', 1)
        r = repr(m)
        assert 'elliptic_curve' in r


# =============================================================================
# CHOW-GRUPPE
# =============================================================================

class TestChowGroup:
    """Tests für die ChowGroup-Klasse."""

    def setup_method(self):
        """Basisvarietät für jeden Test."""
        self.ec = Motive('elliptic_curve', 1)
        self.p2 = Motive('p2', 2)

    def test_construction_valid(self):
        """Gültige Chow-Gruppe kann erstellt werden."""
        cg = ChowGroup(self.ec, codimension=0)
        assert cg.codimension == 0

    def test_codimension_out_of_range_raises(self):
        """Kodimension außerhalb [0, d] wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            ChowGroup(self.ec, codimension=2)  # dim=1, codim=2 ungültig

    def test_codimension_negative_raises(self):
        """Negative Kodimension wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            ChowGroup(self.ec, codimension=-1)

    def test_add_cycle(self):
        """Zykel kann zur Gruppe hinzugefügt werden."""
        cg = ChowGroup(self.p2, codimension=1)
        cg.add_cycle(1, 'hyperplane')
        assert len(cg._generators) == 1

    def test_add_multiple_cycles(self):
        """Mehrere Zykel können hinzugefügt werden."""
        cg = ChowGroup(self.p2, codimension=1)
        cg.add_cycle(1, 'H')
        cg.add_cycle(2, '2H')
        assert len(cg._generators) == 2

    def test_intersection_pairing_complementary(self):
        """Schnittprodukt für komplementäre Kodimensionen (j+k=d)."""
        cg1 = ChowGroup(self.p2, codimension=1)
        cg2 = ChowGroup(self.p2, codimension=1)
        cg1.add_cycle(3, 'Z1')
        cg2.add_cycle(2, 'Z2')
        result = cg1.intersection_pairing(cg2)
        # codim 1 + codim 1 = dim 2 → Ergebnis ist Integer
        assert result == 6

    def test_intersection_pairing_non_complementary_none(self):
        """Schnittprodukt für nicht-komplementäre Kodimensionen → None."""
        cg1 = ChowGroup(self.p2, codimension=0)
        cg2 = ChowGroup(self.p2, codimension=1)
        result = cg1.intersection_pairing(cg2)
        assert result is None

    def test_picard_rank_elliptic_curve(self):
        """Picard-Rang für elliptische Kurve (CH^1): 1."""
        cg = ChowGroup(self.ec, codimension=1)
        assert cg.picard_group_rank() == 1

    def test_picard_rank_p_n(self):
        """Picard-Rang für P^n (CH^1): 1."""
        cg = ChowGroup(self.p2, codimension=1)
        assert cg.picard_group_rank() == 1

    def test_picard_rank_wrong_codim_none(self):
        """picard_group_rank() gibt None zurück, wenn codim ≠ 1."""
        cg = ChowGroup(self.p2, codimension=0)
        assert cg.picard_group_rank() is None

    def test_repr_contains_variety_name(self):
        """String-Repräsentation enthält Varietätsnamen."""
        cg = ChowGroup(self.ec, codimension=1)
        assert 'elliptic_curve' in repr(cg)


# =============================================================================
# MOTIVISCHE KOHOMOLOGIE
# =============================================================================

class TestMotivicCohomology:
    """Tests für die MotivicCohomology-Klasse."""

    def setup_method(self):
        """Motivische Kohomologie für häufig genutzte Varietäten."""
        self.mc_ec = MotivicCohomology(Motive('elliptic_curve', 1))
        self.mc_p2 = MotivicCohomology(Motive('p2', 2))

    # --- dimension() ---

    def test_dimension_chow_diagonal(self):
        """H^{2k,k}: entspricht CH^k(X) ⊗ Q → nutzt b_{2k}."""
        mc = MotivicCohomology(Motive('p2', 2))
        # H^{0,0}: b_0 = 1
        assert mc.dimension(0, 0) == 1
        # H^{4,2}: b_4 = 1 für P²
        assert mc.dimension(4, 2) == 1

    def test_dimension_off_diagonal_zero(self):
        """H^{p,q} für p ≠ 2q und p ≠ 2q-1 → 0."""
        mc = MotivicCohomology(Motive('p2', 2))
        assert mc.dimension(1, 0) == 0
        assert mc.dimension(3, 0) == 0

    def test_dimension_betti_related(self):
        """H^{2q-1,q}: nutzt Betti-Zahl b_{2q-1}."""
        mc = MotivicCohomology(Motive('elliptic_curve', 1))
        # H^{1,1}: b_1 = 2 für elliptische Kurve
        assert mc.dimension(1, 1) == 2

    # --- beilinson_regulator() ---

    def test_beilinson_regulator_elliptic_rank0(self):
        """Beilinson-Regulator (1,1) für elliptische Kurve: 1.0."""
        mc = MotivicCohomology(Motive('elliptic_curve', 1))
        r = mc.beilinson_regulator(1, 1)
        assert abs(r - 1.0) < 1e-12

    def test_beilinson_regulator_zero_dim(self):
        """Regulator für (p,q) mit dim=0 → 0.0."""
        mc = MotivicCohomology(Motive('p2', 2))
        # H^{1,0}: dim=0 → Regulator=0
        r = mc.beilinson_regulator(1, 0)
        assert abs(r) < 1e-12

    def test_beilinson_regulator_positive(self):
        """Beilinson-Regulator ist ≥ 0."""
        mc = MotivicCohomology(Motive('p2', 2))
        r = mc.beilinson_regulator(4, 2)
        assert r >= 0.0

    # --- milnor_k_theory() ---

    def test_milnor_k0_finite_field(self):
        """K^M_0(F_p) = Z (Beschreibung)."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(0, field_char=5)
        assert 'Z' in result['description']
        assert result['n'] == 0

    def test_milnor_k1_finite_field(self):
        """K^M_1(F_p) ≅ Z/(p-1)Z (Bass-Tate)."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(1, field_char=5)
        assert result['order'] == 4  # p-1 = 4

    def test_milnor_kn_finite_field_zero(self):
        """K^M_n(F_p) = 0 für n ≥ 2."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(3, field_char=7)
        assert result['order'] == 1
        assert result['rank'] == 0

    def test_milnor_k2_rationals(self):
        """K^M_2(Q) = Z/2Z (Moore, 1968)."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(2, field_char=0)
        assert result['order'] == 2

    def test_milnor_kn_rationals_large_n_zero(self):
        """K^M_n(Q) = 0 für n ≥ 3."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(5, field_char=0)
        assert result['order'] == 1

    def test_milnor_returns_dict(self):
        """milnor_k_theory() gibt immer ein Dictionary zurück."""
        mc = MotivicCohomology(Motive('point', 0))
        result = mc.milnor_k_theory(1)
        assert isinstance(result, dict)
        assert 'n' in result
        assert 'description' in result


# =============================================================================
# STANDARD-VERMUTUNGEN
# =============================================================================

class TestStandardConjectures:
    """Tests für StandardConjectures."""

    def setup_method(self):
        """Standardvermutungen für typische Varietäten."""
        self.sc_ec = StandardConjectures(Motive('elliptic_curve', 1))
        self.sc_k3 = StandardConjectures(Motive('k3', 2))
        self.sc_p2 = StandardConjectures(Motive('p2', 2))

    # --- Hard-Lefschetz ---

    def test_hard_lefschetz_ec(self):
        """Hard-Lefschetz für elliptische Kurve: Poincaré-Dualität erfüllt."""
        result = self.sc_ec.check_hard_lefschetz()
        assert result['poincare_duality'] is True
        assert result['hard_lefschetz'] == 'verified'

    def test_hard_lefschetz_k3(self):
        """Hard-Lefschetz für K3-Fläche: b_0=b_4=1 und b_1=b_3=0."""
        result = self.sc_k3.check_hard_lefschetz()
        assert result['poincare_duality'] is True

    def test_hard_lefschetz_returns_betti(self):
        """check_hard_lefschetz() enthält Betti-Zahlen."""
        result = self.sc_p2.check_hard_lefschetz()
        assert 'betti_numbers' in result
        assert isinstance(result['betti_numbers'], list)

    # --- Künneth-Formel ---

    def test_kunneth_euler_product(self):
        """Künneth: χ(X×Y) = χ(X)·χ(Y)."""
        m_p1 = Motive('p1', 1)
        m_p2 = Motive('p2', 2)
        sc = StandardConjectures(m_p1)
        result = sc.check_kunneth(m_p2)
        # χ(P¹)=2, χ(P²)=3 → χ(P¹×P²)=6
        assert result['euler_product_check'] is True
        assert result['euler_product'] == 6

    def test_kunneth_product_dimension(self):
        """Künneth: Dimensionen werden addiert."""
        m1 = Motive('p1', 1)
        m2 = Motive('p2', 2)
        sc = StandardConjectures(m1)
        result = sc.check_kunneth(m2)
        assert result['product_dimension'] == 3

    def test_kunneth_betti_list(self):
        """Künneth: Betti-Zahlen des Produkts sind eine Liste."""
        m1 = Motive('p1', 1)
        m2 = Motive('p1', 1)
        sc = StandardConjectures(m1)
        result = sc.check_kunneth(m2)
        assert isinstance(result['kunneth_betti'], list)

    def test_kunneth_conjecture_c_open(self):
        """Künneth-Vermutung C ist im Allgemeinen offen."""
        result = self.sc_p2.check_kunneth(Motive('p1', 1))
        assert result['conjecture_C_status'] == 'open_in_general'

    # --- Hodge-Vermutung ---

    def test_hodge_conjecture_millennium(self):
        """Hodge-Vermutung als Millennium-Problem markiert."""
        result = self.sc_p2.hodge_conjecture_evidence()
        assert 'OPEN' in result['millennium_status']

    def test_hodge_conjecture_elliptic_curve_proven(self):
        """Lefschetz (1,1)-Theorem für elliptische Kurve bewiesen."""
        result = self.sc_ec.hodge_conjecture_evidence()
        assert any('Lefschetz' in case for case in result['proven_cases'])

    def test_hodge_conjecture_returns_dict(self):
        """hodge_conjecture_evidence() gibt Dictionary zurück."""
        result = self.sc_p2.hodge_conjecture_evidence()
        assert isinstance(result, dict)
        assert 'pp_hodge_classes' in result


# =============================================================================
# REALISIERUNGSFUNKTOREN
# =============================================================================

class TestRealizationFunctor:
    """Tests für RealizationFunctor."""

    def setup_method(self):
        """Motiv für Tests."""
        self.m_ec = Motive('elliptic_curve', 1)
        self.m_p2 = Motive('p2', 2)

    # --- Konstruktor ---

    def test_valid_type_betti(self):
        """Realisierung 'betti' kann erstellt werden."""
        rf = RealizationFunctor(self.m_ec, 'betti')
        assert rf.realization_type == 'betti'

    def test_valid_type_de_rham(self):
        """Realisierung 'de_rham' kann erstellt werden."""
        rf = RealizationFunctor(self.m_ec, 'de_rham')
        assert rf.realization_type == 'de_rham'

    def test_valid_type_l_adic(self):
        """Realisierung 'l_adic' kann erstellt werden."""
        rf = RealizationFunctor(self.m_ec, 'l_adic')
        assert rf.realization_type == 'l_adic'

    def test_invalid_type_raises(self):
        """Ungültiger Realisierungstyp wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            RealizationFunctor(self.m_ec, 'motivic')

    # --- compute() ---

    def test_betti_compute_keys(self):
        """Betti-Realisierung enthält erwartete Schlüssel."""
        rf = RealizationFunctor(self.m_ec, 'betti')
        result = rf.compute()
        assert result['type'] == 'betti'
        assert 'vector_space_dimension' in result
        assert 'betti_numbers' in result
        assert 'hodge_numbers' in result

    def test_de_rham_compute_keys(self):
        """de-Rham-Realisierung enthält erwartete Schlüssel."""
        rf = RealizationFunctor(self.m_ec, 'de_rham')
        result = rf.compute()
        assert result['type'] == 'de_rham'
        assert 'hodge_filtration' in result

    def test_l_adic_compute_keys(self):
        """l-adische Realisierung enthält erwartete Schlüssel."""
        rf = RealizationFunctor(self.m_ec, 'l_adic')
        result = rf.compute()
        assert result['type'] == 'l_adic'
        assert 'galois_representation' in result

    def test_betti_dimension_ec(self):
        """Betti-Realisierung der elliptischen Kurve: dim = 1+2+1 = 4."""
        rf = RealizationFunctor(self.m_ec, 'betti')
        result = rf.compute()
        assert result['vector_space_dimension'] == 4

    def test_l_adic_galois_rep_ec(self):
        """Galois-Darstellung der elliptischen Kurve ist 2-dimensional."""
        rf = RealizationFunctor(self.m_ec, 'l_adic')
        result = rf.compute()
        assert '2-dimensional' in result['galois_representation']

    # --- Vergleichsisomorphismus ---

    def test_comparison_betti_de_rham(self):
        """Hodge-de-Rham-Vergleich ist kanonisch."""
        rf_b = RealizationFunctor(self.m_ec, 'betti')
        rf_dr = RealizationFunctor(self.m_ec, 'de_rham')
        comp = rf_b.comparison_isomorphism(rf_dr)
        assert comp['is_canonical'] is True
        assert 'Hodge' in comp['isomorphism']

    def test_comparison_betti_l_adic(self):
        """Étale-Betti-Vergleich ist kanonisch."""
        rf_b = RealizationFunctor(self.m_ec, 'betti')
        rf_l = RealizationFunctor(self.m_ec, 'l_adic')
        comp = rf_b.comparison_isomorphism(rf_l)
        assert comp['is_canonical'] is True

    def test_comparison_de_rham_l_adic(self):
        """p-adischer Hodge-Vergleich ist nicht kanonisch."""
        rf_dr = RealizationFunctor(self.m_ec, 'de_rham')
        rf_l = RealizationFunctor(self.m_ec, 'l_adic')
        comp = rf_dr.comparison_isomorphism(rf_l)
        assert comp['is_canonical'] is False

    def test_comparison_dimension_match(self):
        """Dimensionen stimmen bei gleichem Motiv immer überein."""
        rf_b = RealizationFunctor(self.m_ec, 'betti')
        rf_dr = RealizationFunctor(self.m_ec, 'de_rham')
        comp = rf_b.comparison_isomorphism(rf_dr)
        assert comp['dimension_match'] is True


# =============================================================================
# MOTIVISCHE L-FUNKTION
# =============================================================================

class TestMotivicLFunction:
    """Tests für motivic_l_function()."""

    def test_trivial_motive_s2_positive(self):
        """L(Q(0), 2) = ζ(2) ≈ π²/6 ≈ 1.6449 (positiver Wert)."""
        m = Motive('point', 0)
        result = motivic_l_function(m, s=2.0, prime_bound=100)
        assert result > 1.0

    def test_trivial_motive_large_s_converges(self):
        """L(Q(0), s) konvergiert für großes s gegen 1."""
        m = Motive('point', 0)
        result = motivic_l_function(m, s=10.0, prime_bound=50)
        # ζ(10) ≈ 1.000995
        assert abs(result - 1.0) < 0.01

    def test_elliptic_curve_l_function_positive(self):
        """Hasse-Weil-L-Funktion der elliptischen Kurve: L(E, s) > 0 für s >> 1."""
        m = Motive('elliptic_curve', 1)
        result = motivic_l_function(m, s=2.0, prime_bound=30)
        assert result > 0.0

    def test_p1_l_function_positive(self):
        """L(h(P¹), s) = ζ(s)ζ(s-1): positiv für s > 1."""
        m = Motive('p1', 1)
        result = motivic_l_function(m, s=3.0, prime_bound=30)
        assert result > 0.0

    def test_s_too_small_raises(self):
        """s ≤ 0.5 wirft InvalidInputError (außerhalb Konvergenzbereich)."""
        m = Motive('elliptic_curve', 1)
        with pytest.raises(InvalidInputError):
            motivic_l_function(m, s=0.3)

    def test_s_exactly_at_boundary_raises(self):
        """s = 0.5 (Grenzwert) wirft InvalidInputError."""
        m = Motive('point', 0)
        with pytest.raises(InvalidInputError):
            motivic_l_function(m, s=0.5)

    def test_prime_bound_affects_result(self):
        """Größerer prime_bound gibt genaueres Ergebnis (ζ(2) → π²/6)."""
        m = Motive('point', 0)
        r_small = motivic_l_function(m, s=2.0, prime_bound=20)
        r_large = motivic_l_function(m, s=2.0, prime_bound=200)
        # Mit mehr Primzahlen nähern wir uns π²/6 ≈ 1.6449
        assert abs(r_large - math.pi**2 / 6) < abs(r_small - math.pi**2 / 6)

    def test_general_motive_l_function(self):
        """Allgemeines Motiv mit Tate-Twist: L(M(1), s) korrekt berechnet."""
        m = Motive('abelian_1', 1, tate_twist=1)
        result = motivic_l_function(m, s=3.0, prime_bound=20)
        assert result > 0.0

    def test_l_function_returns_positive_float(self):
        """motivic_l_function() gibt immer eine positive Zahl zurück."""
        m = Motive('elliptic_curve', 1)
        for s in [1.5, 2.0, 3.0, 5.0]:
            result = motivic_l_function(m, s=s, prime_bound=20)
            assert isinstance(result, float)
            assert result > 0.0
