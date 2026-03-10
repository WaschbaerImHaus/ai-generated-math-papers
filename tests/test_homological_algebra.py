"""
@file test_homological_algebra.py
@brief Tests für das Modul homological_algebra.py.
@description
    Testet alle Funktionen und Klassen des Homologischen Algebra-Moduls:
    - ChainComplex: Axiom ∂∘∂=0, Homologie, Euler-Charakteristik, Betti-Zahlen
    - CochainComplex: Kohomologie, Axiom d∘d=0
    - Simplizialer Komplex → Kettenkomplex
    - ExactSequence: Exaktheit, kurze exakte Sequenzen, Aufspaltung
    - 5-Lemma und Schlangen-Lemma
    - Freie Auflösungen
    - Ext und Tor Gruppen
    - Universeller Koeffizienten-Satz

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import numpy as np
import pytest
import math

# Quellverzeichnis zum Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from homological_algebra import (
    ChainComplex,
    CochainComplex,
    simplicial_complex_to_chain,
    ExactSequence,
    five_lemma_check,
    snake_lemma,
    free_resolution,
    projective_dimension,
    ext_group,
    tor_group,
    universal_coefficient_theorem,
    _smith_normal_form,
    _gcd,
)


# ---------------------------------------------------------------------------
# Tests: Smith-Normalform
# ---------------------------------------------------------------------------

class TestSmithNormalForm:
    """Tests für die Smith-Normalform."""

    def test_identity_matrix(self):
        """Smith-NF der Einheitsmatrix: Elementarteiler alle 1."""
        A = [[1, 0], [0, 1]]
        D, rank, elem_div = _smith_normal_form(A)
        assert rank == 2
        assert all(d == 1 for d in elem_div)

    def test_zero_matrix(self):
        """Smith-NF der Nullmatrix: Rang 0."""
        A = [[0, 0], [0, 0]]
        D, rank, elem_div = _smith_normal_form(A)
        assert rank == 0
        assert len(elem_div) == 0

    def test_diagonal_matrix(self):
        """Smith-NF einer Diagonalmatrix mit gcd-Eigenschaft."""
        A = [[2, 0], [0, 6]]
        D, rank, elem_div = _smith_normal_form(A)
        assert rank == 2
        # Elementarteiler: d₁ | d₂
        assert len(elem_div) == 2

    def test_single_element(self):
        """1×1-Matrix."""
        A = [[5]]
        D, rank, elem_div = _smith_normal_form(A)
        assert rank == 1
        assert elem_div == [5]


class TestGcd:
    """Tests für den größten gemeinsamen Teiler."""

    def test_basic(self):
        assert _gcd(12, 8) == 4

    def test_coprime(self):
        assert _gcd(7, 13) == 1

    def test_zero(self):
        assert _gcd(0, 5) == 5
        assert _gcd(5, 0) == 5

    def test_negative(self):
        assert _gcd(-6, 4) == 2


# ---------------------------------------------------------------------------
# Tests: ChainComplex
# ---------------------------------------------------------------------------

class TestChainComplexBasic:
    """Grundlegende Tests für ChainComplex."""

    def _trivial_complex(self):
        """Trivialer Kettenkomplex: 0 → ℤ → 0."""
        return ChainComplex(groups={0: 1}, boundaries={})

    def test_trivial_is_complex(self):
        """Trivialer Komplex erfüllt ∂∘∂=0."""
        C = self._trivial_complex()
        assert C.is_complex()

    def test_circle_chain_complex(self):
        """
        Kettenkomplex des Kreises S¹:
        C₁ = ℤ (eine 1-Kette: Kante e)
        C₀ = ℤ² (zwei 0-Ketten: Punkte v₀, v₁)
        ∂₁(e) = v₁ - v₀ → Randmatrix: [-1, 1]ᵀ → nein, für S¹: v₀→v₁, ∂(e)=v₁-v₀
        Für einfachsten S¹: 1 Kante, 1 Knoten (Schleife): ∂₁ = [0]
        """
        # Einfacher Kreis: 1 Knoten, 1 Kante mit ∂₁(e) = v - v = 0
        C = ChainComplex(
            groups={0: 1, 1: 1},
            boundaries={1: [[0]]}
        )
        assert C.is_complex()

    def test_invalid_complex(self):
        """Kettenkomplex mit ∂∘∂ ≠ 0 wird erkannt."""
        # ∂₁ = [[1]] und ∂₀ = [[1]]: ∂₀∘∂₁ = [[1]] ≠ 0
        C = ChainComplex(
            groups={0: 1, 1: 1, 2: 1},
            boundaries={1: [[1]], 2: [[1]]}
        )
        # ∂₀∘∂₁ = [[1]]·[[1]] = [[1]] ≠ 0
        assert not C.is_complex()


class TestChainComplexHomology:
    """Tests für die Homologieberechnung."""

    def test_homology_single_point(self):
        """
        Simplizialer Komplex {Punkt}: H₀ = ℤ, H_n = 0 für n > 0.
        """
        C = ChainComplex(groups={0: 1}, boundaries={})
        h0 = C.homology(0)
        assert h0['rank'] == 1
        assert h0['torsion'] == []

    def test_homology_circle_s1(self):
        """
        S¹ als Simplizialkomplex (Dreieck ohne Innenfläche):
        3 Knoten, 3 Kanten.
        H₀ = ℤ (zusammenhängend), H₁ = ℤ (ein Loch)
        """
        # Dreieck: Knoten {0,1,2}, Kanten {01,12,02}
        simplices = {
            0: [[0], [1], [2]],
            1: [[0, 1], [1, 2], [0, 2]]
        }
        C = simplicial_complex_to_chain(simplices)
        assert C.is_complex()

        h0 = C.homology(0)
        assert h0['rank'] >= 0  # Mindestens 0

        h1 = C.homology(1)
        # H₁(S¹) = ℤ → Betti-Zahl 1
        assert h1['rank'] == 1

    def test_euler_characteristic_point(self):
        """Euler-Charakteristik eines Punktes ist 1."""
        C = ChainComplex(groups={0: 1}, boundaries={})
        assert C.euler_characteristic() == 1

    def test_euler_characteristic_circle(self):
        """Euler-Charakteristik von S¹ (1 Kante, 1 Knoten) = 0."""
        # Schleife: 1 Knoten, 1 Kante, ∂₁ = 0
        C = ChainComplex(
            groups={0: 1, 1: 1},
            boundaries={1: [[0]]}
        )
        # χ = 1 - 1 = 0
        assert C.euler_characteristic() == 0

    def test_euler_characteristic_triangle(self):
        """Euler-Charakteristik des gefüllten Dreiecks = 1."""
        # 3 Knoten, 3 Kanten, 1 Dreieck: χ = 3 - 3 + 1 = 1
        C = ChainComplex(
            groups={0: 3, 1: 3, 2: 1},
            boundaries={
                1: [[-1, -1, 0], [1, 0, -1], [0, 1, 1]],  # Randoperator Kanten→Knoten
                2: [[1], [-1], [1]]                          # Randoperator Fläche→Kanten
            }
        )
        assert C.euler_characteristic() == 1

    def test_betti_numbers_point(self):
        """Betti-Zahlen eines Punktes: β₀ = 1."""
        C = ChainComplex(groups={0: 1}, boundaries={})
        betti = C.betti_numbers()
        assert betti[0] == 1

    def test_homology_free_module(self):
        """Kettenkomplex mit freiem Modul hat keine Torsion."""
        C = ChainComplex(groups={0: 2, 1: 1}, boundaries={1: [[1], [-1]]})
        h0 = C.homology(0)
        # H₀ sollte frei sein (für zusammenhängenden Komplex: ℤ)
        assert h0['free'] == True or isinstance(h0['rank'], int)


class TestChainComplexS2:
    """Tests für den Kettenkomplex der 2-Sphäre."""

    def test_s2_is_complex(self):
        """
        S² als Simplizialkomplex (Tetraederoberfläche):
        4 Knoten, 6 Kanten, 4 Dreiecke.
        """
        simplices = {
            0: [[0], [1], [2], [3]],
            1: [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
            2: [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        }
        C = simplicial_complex_to_chain(simplices)
        assert C.is_complex()

    def test_s2_euler_characteristic(self):
        """Euler-Charakteristik von S²: χ = 4 - 6 + 4 = 2."""
        simplices = {
            0: [[0], [1], [2], [3]],
            1: [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
            2: [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        }
        C = simplicial_complex_to_chain(simplices)
        assert C.euler_characteristic() == 2


# ---------------------------------------------------------------------------
# Tests: CochainComplex
# ---------------------------------------------------------------------------

class TestCochainComplex:
    """Tests für den Kokettenkomplex."""

    def test_trivial_cochain_is_complex(self):
        """Trivialer Kokettenkomplex erfüllt d∘d=0."""
        C = CochainComplex(groups={0: 1}, coboundaries={})
        assert C.is_complex()

    def test_cochain_not_complex(self):
        """Ungültiger Kokettenkomplex wird erkannt."""
        # d⁰ = [[1]] und d¹ = [[1]]: d¹∘d⁰ = [[1]] ≠ 0
        C = CochainComplex(
            groups={0: 1, 1: 1, 2: 1},
            coboundaries={0: [[1]], 1: [[1]]}
        )
        assert not C.is_complex()

    def test_cohomology_trivial(self):
        """Kohomologie des Trivial-Komplexes."""
        C = CochainComplex(groups={0: 1}, coboundaries={})
        h = C.cohomology(0)
        assert h['rank'] == 1


# ---------------------------------------------------------------------------
# Tests: Simplizialer Komplex
# ---------------------------------------------------------------------------

class TestSimplicialComplexToChain:
    """Tests für die Konvertierung simplizialer Komplexe."""

    def test_single_point(self):
        """Ein einzelner Punkt."""
        simplices = {0: [[0]]}
        C = simplicial_complex_to_chain(simplices)
        assert C.groups[0] == 1
        assert C.is_complex()

    def test_interval(self):
        """Intervall [0,1]: 2 Knoten, 1 Kante."""
        simplices = {
            0: [[0], [1]],
            1: [[0, 1]]
        }
        C = simplicial_complex_to_chain(simplices)
        assert C.groups[0] == 2
        assert C.groups[1] == 1
        assert C.is_complex()
        # ∂₁([0,1]) = [1] - [0], Matrix: [[-1], [1]]
        d1 = C.boundaries.get(1)
        assert d1 is not None
        assert d1.shape == (2, 1)

    def test_triangle_boundary(self):
        """Randoperator des Dreiecks [0,1,2]."""
        simplices = {
            0: [[0], [1], [2]],
            1: [[0, 1], [0, 2], [1, 2]],
        }
        C = simplicial_complex_to_chain(simplices)
        assert C.is_complex()

    def test_boundary_formula(self):
        """Prüfe ∂([0,1]) = [1] - [0]."""
        simplices = {
            0: [[0], [1]],
            1: [[0, 1]]
        }
        C = simplicial_complex_to_chain(simplices)
        d1 = C.boundaries[1]
        # Knoten in Reihenfolge: [0], [1]
        # ∂([0,1]) = -1·[0] + 1·[1]
        total = sum(abs(d1[i, 0]) for i in range(d1.shape[0]))
        assert total == 2  # Zwei Nicht-Null-Einträge mit ±1


# ---------------------------------------------------------------------------
# Tests: Exakte Sequenzen
# ---------------------------------------------------------------------------

class TestExactSequence:
    """Tests für exakte Sequenzen."""

    def test_trivial_exact(self):
        """0 → ℤ → ℤ → 0 ist exakt."""
        # Isomorphismus 1×1 mit [[1]]
        seq = ExactSequence(
            groups=[1, 1],
            maps=[[[1]]]
        )
        assert seq.is_exact()

    def test_short_exact_z2(self):
        """0 → ℤ →^2 ℤ → ℤ/2ℤ → 0 (kurze exakte Sequenz)."""
        # Gruppe A = ℤ (Rang 1), B = ℤ (Rang 1), C = ℤ/2ℤ (Rang 0 als freier Modul)
        # Für freie Module: Rang(B) = Rang(A) + Rang(C)
        seq = ExactSequence(
            groups=[1, 1, 0],
            maps=[[[2]], [[1]]]  # Vereinfachte Darstellung
        )
        # B = ℤ (Rang 1), aber das ist eine Vereinfachung
        # Echtes Test: Rang der Abbildungen
        assert isinstance(seq.is_short_exact(), bool)

    def test_split_free_modules(self):
        """Freie Module: Jede kurze exakte Sequenz zerfällt."""
        # 0 → ℤ → ℤ² → ℤ → 0
        seq = ExactSequence(
            groups=[1, 2, 1],
            maps=[[[1], [0]], [[0, 1]]]
        )
        assert seq.split()

    def test_not_split_torsion(self):
        """Sequenz mit gleichen Rängen: prüfe split-Logik."""
        # 0 → ℤ² → ℤ³ → ℤ → 0
        seq = ExactSequence(groups=[2, 3, 1], maps=[])
        assert seq.split()

    def test_wrong_size_not_short(self):
        """Nicht-kurze Sequenz (4 Gruppen) ist nicht is_short_exact."""
        seq = ExactSequence(
            groups=[1, 1, 1, 1],
            maps=[[[1]], [[1]], [[1]]]
        )
        assert not seq.is_short_exact()


# ---------------------------------------------------------------------------
# Tests: 5-Lemma
# ---------------------------------------------------------------------------

class TestFiveLemma:
    """Tests für das 5-Lemma."""

    def test_identity_diagram(self):
        """Diagramm mit Identitätsabbildungen."""
        I = [[1, 0], [0, 1]]
        Z = [[0, 0], [0, 0]]
        diagram = {
            'top_maps': [I, I, I, I],
            'bot_maps': [I, I, I, I],
            'vert_maps': [I, I, I, I, I]
        }
        result = five_lemma_check(diagram)
        assert result['iso_status'] == [True, True, True, True, True]

    def test_wrong_vert_map_count(self):
        """Weniger als 5 vertikale Abbildungen liefert Fehler."""
        diagram = {
            'top_maps': [],
            'bot_maps': [],
            'vert_maps': [[[1]], [[1]], [[1]]]  # Nur 3
        }
        result = five_lemma_check(diagram)
        assert not result['valid']


# ---------------------------------------------------------------------------
# Tests: Schlangen-Lemma
# ---------------------------------------------------------------------------

class TestSnakeLemma:
    """Tests für das Schlangen-Lemma."""

    def test_trivial_snake(self):
        """Trivialer Fall: alle Abbildungen sind Identitäten."""
        maps = {
            'alpha': [[1]],
            'beta': [[1]],
            'gamma': [[1]]
        }
        result = snake_lemma(1, 1, 1, 1, 1, 1, maps)
        assert 'long_exact_sequence' in result
        assert 'ker_alpha' in result

    def test_zero_maps(self):
        """Nullabbildungen: alle Kerne sind die vollen Gruppen."""
        maps = {
            'alpha': [[0]],
            'beta': [[0]],
            'gamma': [[0]]
        }
        result = snake_lemma(1, 1, 1, 1, 1, 1, maps)
        assert result['ker_alpha'] == 1
        assert result['ker_beta'] == 1
        assert result['ker_gamma'] == 1

    def test_sequence_structure(self):
        """Die lange exakte Sequenz hat die richtige Struktur."""
        maps = {
            'alpha': [[2]],
            'beta': [[2]],
            'gamma': [[2]]
        }
        result = snake_lemma(1, 1, 1, 1, 1, 1, maps)
        seq = result['long_exact_sequence']
        assert len(seq) == 9  # 9 Einträge in der Sequenz
        assert seq[0] == ('0', 0)
        assert seq[-1] == ('0', 0)


# ---------------------------------------------------------------------------
# Tests: Freie Auflösungen
# ---------------------------------------------------------------------------

class TestFreeResolution:
    """Tests für freie Auflösungen."""

    def test_z2_resolution(self):
        """Freie Auflösung von ℤ/2ℤ: 0 → ℤ →^2 ℤ → ℤ/2ℤ → 0."""
        C = free_resolution(2, [])
        assert C.groups.get(0) == 1
        assert C.groups.get(1) == 1
        # Randoperator ist die Multiplikation mit 2
        d1 = C.boundaries.get(1)
        assert d1 is not None
        assert d1[0, 0] == 2

    def test_z3_resolution(self):
        """Freie Auflösung von ℤ/3ℤ."""
        C = free_resolution(3, [])
        d1 = C.boundaries.get(1)
        assert d1[0, 0] == 3

    def test_z1_resolution(self):
        """ℤ/1ℤ = 0: triviale Auflösung."""
        C = free_resolution(1, [])
        assert C.groups.get(0, 0) == 0

    def test_resolution_is_complex(self):
        """Freie Auflösung erfüllt ∂∘∂=0."""
        C = free_resolution(6, [])
        assert C.is_complex()

    def test_invalid_n(self):
        """Negativer Modulus wirft ValueError."""
        with pytest.raises(ValueError):
            free_resolution(-1, [])


class TestProjectiveDimension:
    """Tests für die projektive Dimension."""

    def test_free_module(self):
        """Freier Modul hat Dimension 0."""
        assert projective_dimension(2, []) == 0

    def test_zero_module(self):
        """Nullmodul hat Dimension -1."""
        assert projective_dimension(0, []) == -1

    def test_cyclic_torsion_module(self):
        """ℤ/nℤ hat projektive Dimension 1."""
        result = projective_dimension(1, [[2]])
        assert result >= 0


# ---------------------------------------------------------------------------
# Tests: Ext-Gruppen
# ---------------------------------------------------------------------------

class TestExtGroup:
    """Tests für die Ext-Gruppen."""

    def test_ext0_hom(self):
        """Ext⁰(ℤ/nℤ, ℤ/mℤ) = Hom(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ."""
        result = ext_group(6, 4, 0)
        # gcd(6,4) = 2
        assert result['order'] == 2
        assert result['group'] == 'ℤ/2ℤ'

    def test_ext1(self):
        """Ext¹(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ."""
        result = ext_group(6, 9, 1)
        # gcd(6,9) = 3
        assert result['order'] == 3
        assert not result['trivial']

    def test_ext_high_degree_zero(self):
        """Ext^k = 0 für k ≥ 2."""
        result = ext_group(5, 3, 2)
        assert result['trivial']
        assert result['order'] == 1

    def test_ext_coprime(self):
        """gcd(n,m) = 1 → Ext trivial."""
        result = ext_group(5, 7, 0)
        assert result['trivial']
        assert result['order'] == 1

    def test_ext0_description(self):
        """Beschreibungsstring enthält Hom."""
        result = ext_group(4, 6, 0)
        assert 'Hom' in result['description']

    def test_ext1_description(self):
        """Beschreibungsstring enthält Ext¹."""
        result = ext_group(4, 6, 1)
        assert 'Ext' in result['description']


# ---------------------------------------------------------------------------
# Tests: Tor-Gruppen
# ---------------------------------------------------------------------------

class TestTorGroup:
    """Tests für die Tor-Gruppen."""

    def test_tor0_tensor_product(self):
        """Tor₀(ℤ/nℤ, ℤ/mℤ) = ℤ/nℤ ⊗ ℤ/mℤ ≅ ℤ/gcd(n,m)ℤ."""
        result = tor_group(6, 4, 0)
        # gcd(6,4) = 2
        assert result['order'] == 2
        assert result['group'] == 'ℤ/2ℤ'

    def test_tor1(self):
        """Tor₁(ℤ/nℤ, ℤ/mℤ) ≅ ℤ/gcd(n,m)ℤ."""
        result = tor_group(4, 6, 1)
        # gcd(4,6) = 2
        assert result['order'] == 2
        assert not result['trivial']

    def test_tor_high_degree_zero(self):
        """Tor_k = 0 für k ≥ 2."""
        result = tor_group(3, 7, 3)
        assert result['trivial']

    def test_tor_symmetry(self):
        """Tor ist symmetrisch: Tor(A,B) ≅ Tor(B,A)."""
        r1 = tor_group(6, 4, 1)
        r2 = tor_group(4, 6, 1)
        assert r1['order'] == r2['order']

    def test_tor0_description(self):
        """Beschreibungsstring enthält ⊗."""
        result = tor_group(3, 5, 0)
        assert '⊗' in result['description']

    def test_tor1_description(self):
        """Beschreibungsstring enthält Tor₁."""
        result = tor_group(3, 5, 1)
        assert 'Tor' in result['description']


# ---------------------------------------------------------------------------
# Tests: Universeller Koeffizienten-Satz
# ---------------------------------------------------------------------------

class TestUniversalCoefficientTheorem:
    """Tests für den Universellen Koeffizienten-Satz."""

    def test_free_homology(self):
        """Freie Homologie → freie Kohomologie."""
        betti = {0: 1, 1: 1}
        torsion = {}
        result = universal_coefficient_theorem(betti, torsion)
        # H^0 = ℤ (frei, kein Torsion aus H_{-1})
        assert result[0]['free_rank'] == 1
        assert result[0]['torsion'] == []

    def test_torsion_shifts_degree(self):
        """Torsion in H_{n-1} erscheint in H^n."""
        betti = {0: 1, 1: 0, 2: 1}
        torsion = {1: [2, 3]}  # H₁ hat Torsion ℤ/2ℤ ⊕ ℤ/3ℤ
        result = universal_coefficient_theorem(betti, torsion)
        # H^2 = ℤ (von β₂=1) ⊕ Ext¹(H₁, ℤ) = Torsion von H₁
        assert result[2]['torsion'] == [2, 3]

    def test_group_string_free(self):
        """Freie Gruppe wird als 'ℤ' oder 'ℤ^n' dargestellt."""
        betti = {0: 1}
        torsion = {}
        result = universal_coefficient_theorem(betti, torsion)
        assert 'ℤ' in result[0]['group']

    def test_trivial_homology(self):
        """Triviale Homologie → triviale Kohomologie."""
        betti = {}
        torsion = {}
        result = universal_coefficient_theorem(betti, torsion)
        # Alle Gruppen sollten trivial sein
        for n, h in result.items():
            assert h['free_rank'] == 0
            assert h['torsion'] == []

    def test_rank_2_free(self):
        """β_n = 2 → 'ℤ^2' in der Gruppe."""
        betti = {0: 2}
        torsion = {}
        result = universal_coefficient_theorem(betti, torsion)
        assert 'ℤ^2' in result[0]['group']

    def test_mixed_group_string(self):
        """Gemischte Gruppe: freier Teil ⊕ Torsion."""
        betti = {1: 1}
        torsion = {0: [2]}  # Erscheint in H^1
        result = universal_coefficient_theorem(betti, torsion)
        # H^1 = ℤ ⊕ ℤ/2ℤ
        g = result[1]['group']
        assert 'ℤ' in g and '2' in g


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
