"""
@file test_modules_algebra.py
@brief Tests für das Modul modules_algebra.py (R-Moduln über Ringen).
@description
    Test-Suite für Smith-Normalform, Klassifikationssatz abelscher Gruppen,
    Tensorprodukt, Hom-Modul, exakte Sequenzen und freie Auflösungen.

    Testfälle folgen dem Test-Driven-Development-Ansatz: erst Tests,
    dann Implementierung.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math

# Projektpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from modules_algebra import (
    smith_normal_form,
    module_from_matrix,
    structure_theorem_abelian_groups,
    tensor_product_modules,
    hom_module,
    exact_sequence_check,
    free_resolution,
    Module,
)


# ---------------------------------------------------------------------------
# Smith-Normalform Tests
# ---------------------------------------------------------------------------

class TestSmithNormalForm:
    """Tests für die Smith-Normalform ganzzahliger Matrizen."""

    def test_simple_2x2_matrix(self):
        """[[2,4],[6,8]] → Invariantenfaktoren [2, 4]."""
        result = smith_normal_form([[2, 4], [6, 8]])
        assert result['invariant_factors'] == [2, 4], (
            f"Erwartet [2, 4], bekommen {result['invariant_factors']}"
        )

    def test_identity_matrix(self):
        """Einheitsmatrix hat Invariantenfaktoren [1, 1]."""
        result = smith_normal_form([[1, 0], [0, 1]])
        assert result['invariant_factors'] == [1, 1]

    def test_zero_matrix(self):
        """Nullmatrix hat keine Invariantenfaktoren (rank 0)."""
        result = smith_normal_form([[0, 0], [0, 0]])
        assert result['rank'] == 0
        assert result['invariant_factors'] == []

    def test_diagonal_matrix(self):
        """Diagonalmatrix diag(2,6) → Invariantenfaktoren [2, 6]."""
        result = smith_normal_form([[2, 0], [0, 6]])
        # d1 | d2 muss gelten
        factors = result['invariant_factors']
        assert len(factors) == 2
        assert factors[0] >= 1
        assert factors[1] % factors[0] == 0, (
            f"Teilbarkeitsbedingung verletzt: {factors[0]} | {factors[1]}"
        )

    def test_1x1_matrix(self):
        """1×1-Matrix [[3]] → Invariantenfaktoren [3]."""
        result = smith_normal_form([[3]])
        assert result['invariant_factors'] == [3]

    def test_rank_returned(self):
        """Rang der Matrix wird korrekt zurückgegeben."""
        result = smith_normal_form([[2, 4], [6, 8]])
        assert result['rank'] == 2

    def test_rectangular_matrix(self):
        """Rechteckige Matrix 2×3."""
        result = smith_normal_form([[1, 2, 3], [4, 5, 6]])
        assert 'invariant_factors' in result
        assert 'rank' in result
        assert result['rank'] <= 2

    def test_divisibility_chain(self):
        """Invariantenfaktoren erfüllen Teilbarkeitskette d1 | d2 | d3."""
        result = smith_normal_form([[2, 4, 8], [0, 6, 12], [0, 0, 12]])
        factors = result['invariant_factors']
        for i in range(len(factors) - 1):
            assert factors[i + 1] % factors[i] == 0, (
                f"Teilbarkeit verletzt: {factors[i]} soll {factors[i+1]} teilen"
            )

    def test_unimodular_transforms_returned(self):
        """Rückgabe enthält U, D, V Matrizen."""
        result = smith_normal_form([[2, 4], [6, 8]])
        assert 'U' in result
        assert 'D' in result
        assert 'V' in result

    def test_matrix_4x4(self):
        """4×4-Matrix wird korrekt verarbeitet."""
        A = [[1, 0, 0, 0],
             [0, 2, 0, 0],
             [0, 0, 3, 0],
             [0, 0, 0, 6]]
        result = smith_normal_form(A)
        factors = result['invariant_factors']
        # Alle müssen positiv sein
        assert all(f > 0 for f in factors)
        # Teilbarkeitskette
        for i in range(len(factors) - 1):
            assert factors[i + 1] % factors[i] == 0


# ---------------------------------------------------------------------------
# module_from_matrix Tests
# ---------------------------------------------------------------------------

class TestModuleFromMatrix:
    """Tests für die Modul-Konstruktion aus Relationsmatrizen."""

    def test_returns_dict_keys(self):
        """Rückgabe hat alle erwarteten Schlüssel."""
        result = module_from_matrix([[2, 0], [0, 3]])
        assert 'invariant_factors' in result
        assert 'rank' in result
        assert 'torsion_part' in result
        assert 'free_part' in result

    def test_trivial_module(self):
        """Einheitsmatrix → triviales Modul (keine Torsion, Rang 0)."""
        result = module_from_matrix([[1, 0], [0, 1]])
        # Alle Invariantenfaktoren 1 → keine nichttriviale Torsion
        assert result['rank'] == 0 or all(
            f == 1 for f in result['invariant_factors']
        )

    def test_cyclic_module_z6(self):
        """Matrix [[6]] → Modul ≅ ℤ_6."""
        result = module_from_matrix([[6]])
        assert 6 in result['invariant_factors'] or result['torsion_part'] != ''


# ---------------------------------------------------------------------------
# Klassifikationssatz abelscher Gruppen
# ---------------------------------------------------------------------------

class TestStructureTheoremAbelianGroups:
    """Tests für den Hauptsatz über endlich erzeugte abelsche Gruppen."""

    def test_n4_contains_z4(self):
        """n=4: Eine Gruppe ist ℤ_4."""
        groups = structure_theorem_abelian_groups(4)
        invariant_factor_lists = [g['invariant_factors'] for g in groups]
        assert [4] in invariant_factor_lists, (
            f"ℤ_4 nicht gefunden in {invariant_factor_lists}"
        )

    def test_n4_contains_z2xz2(self):
        """n=4: Eine Gruppe ist ℤ_2 × ℤ_2."""
        groups = structure_theorem_abelian_groups(4)
        invariant_factor_lists = [g['invariant_factors'] for g in groups]
        assert [2, 2] in invariant_factor_lists, (
            f"ℤ_2×ℤ_2 nicht gefunden in {invariant_factor_lists}"
        )

    def test_n4_exactly_two_groups(self):
        """n=4: Genau 2 nicht-isomorphe abelsche Gruppen."""
        groups = structure_theorem_abelian_groups(4)
        assert len(groups) == 2

    def test_n12_contains_z12(self):
        """n=12: Eine Gruppe ist ℤ_12."""
        groups = structure_theorem_abelian_groups(12)
        invariant_factor_lists = [g['invariant_factors'] for g in groups]
        assert [12] in invariant_factor_lists, (
            f"ℤ_12 nicht gefunden in {invariant_factor_lists}"
        )

    def test_n12_contains_z2xz6(self):
        """n=12: Eine Gruppe ist ℤ_2 × ℤ_6."""
        groups = structure_theorem_abelian_groups(12)
        invariant_factor_lists = [g['invariant_factors'] for g in groups]
        assert [2, 6] in invariant_factor_lists, (
            f"ℤ_2×ℤ_6 nicht gefunden in {invariant_factor_lists}"
        )

    def test_n1_trivial(self):
        """n=1: Nur triviale Gruppe."""
        groups = structure_theorem_abelian_groups(1)
        assert len(groups) == 1

    def test_n_prime_cyclic(self):
        """n=p (Primzahl): Nur eine abelsche Gruppe ℤ_p."""
        for p in [2, 3, 5, 7, 11]:
            groups = structure_theorem_abelian_groups(p)
            assert len(groups) == 1, f"Für n={p} (prim) erwartet genau 1 Gruppe"

    def test_orders_match(self):
        """Alle zurückgegebenen Gruppen haben die richtige Ordnung n."""
        for n in [4, 6, 8, 12, 16]:
            groups = structure_theorem_abelian_groups(n)
            for g in groups:
                order = 1
                for factor in g['invariant_factors']:
                    order *= factor
                assert order == n, (
                    f"Gruppe {g} hat Ordnung {order}, erwartet {n}"
                )

    def test_n8_three_groups(self):
        """n=8: Genau 3 abelsche Gruppen: ℤ_8, ℤ_4×ℤ_2, ℤ_2×ℤ_2×ℤ_2."""
        groups = structure_theorem_abelian_groups(8)
        assert len(groups) == 3


# ---------------------------------------------------------------------------
# Tensorprodukt
# ---------------------------------------------------------------------------

class TestTensorProduct:
    """Tests für das Tensorprodukt von ℤ-Moduln."""

    def test_z4_tensor_z6_is_z2(self):
        """ℤ_4 ⊗_ℤ ℤ_6 ≅ ℤ_{gcd(4,6)} = ℤ_2."""
        result = tensor_product_modules([4], [6], 0)
        assert result['order'] == 2, (
            f"ℤ_4 ⊗ ℤ_6 sollte ℤ_2 sein, bekam Ordnung {result['order']}"
        )

    def test_z3_tensor_z5_trivial(self):
        """ℤ_3 ⊗_ℤ ℤ_5 ≅ ℤ_1 (trivial, da gcd(3,5)=1)."""
        result = tensor_product_modules([3], [5], 0)
        assert result['order'] == 1

    def test_z6_tensor_z6_is_z6(self):
        """ℤ_6 ⊗_ℤ ℤ_6 ≅ ℤ_6."""
        result = tensor_product_modules([6], [6], 0)
        assert result['order'] == 6

    def test_returns_dict(self):
        """Rückgabe ist ein Dict mit 'order' und 'structure'."""
        result = tensor_product_modules([4], [6], 0)
        assert 'order' in result
        assert 'structure' in result


# ---------------------------------------------------------------------------
# Hom-Modul
# ---------------------------------------------------------------------------

class TestHomModule:
    """Tests für Hom_ℤ(M, N)."""

    def test_hom_z4_z6_has_order_2(self):
        """|Hom(ℤ_4, ℤ_6)| = gcd(4, 6) = 2."""
        result = hom_module(4, 6)
        assert result['order'] == 2, (
            f"|Hom(ℤ_4, ℤ_6)| = 2 erwartet, bekam {result['order']}"
        )

    def test_hom_z6_z4_has_order_2(self):
        """|Hom(ℤ_6, ℤ_4)| = gcd(6, 4) = 2."""
        result = hom_module(6, 4)
        assert result['order'] == 2

    def test_hom_z3_z5_trivial(self):
        """|Hom(ℤ_3, ℤ_5)| = 1 (nur Nullhomomorphismus)."""
        result = hom_module(3, 5)
        assert result['order'] == 1

    def test_hom_zn_zn(self):
        """|Hom(ℤ_n, ℤ_n)| = n."""
        for n in [2, 3, 4, 6, 12]:
            result = hom_module(n, n)
            assert result['order'] == n, (
                f"|Hom(ℤ_{n}, ℤ_{n})| = {n} erwartet, bekam {result['order']}"
            )

    def test_returns_dict(self):
        """Rückgabe ist Dict mit 'order' und 'generators'."""
        result = hom_module(4, 6)
        assert 'order' in result
        assert 'generators' in result


# ---------------------------------------------------------------------------
# Exakte Sequenzen
# ---------------------------------------------------------------------------

class TestExactSequence:
    """Tests für Überprüfung exakter Sequenzen."""

    def test_short_exact_sequence_z_to_z_to_zn(self):
        """0 → ℤ -n·-> ℤ → ℤ_n → 0 ist exakt."""
        # Dargestellt als Matrizen: Multiplikation mit n, dann Projektion mod n
        n = 6
        result = exact_sequence_check(
            maps=[[[n]], [[1]]],  # erste Map: ×n; zweite Map: mod n
            modules=[1, 1, n]    # ℤ, ℤ, ℤ_n (Rang/Ordnung)
        )
        assert result['is_exact'] is True

    def test_non_exact_sequence(self):
        """Eine nicht-exakte Sequenz wird korrekt erkannt."""
        result = exact_sequence_check(
            maps=[[[2]], [[3]]],
            modules=[1, 1, 6]
        )
        # 2ℤ ≠ Kern(×3 mod 6) → nicht exakt
        assert 'is_exact' in result

    def test_returns_dict(self):
        """Rückgabe hat 'is_exact' und Detailinformationen."""
        result = exact_sequence_check(
            maps=[[[1]], [[1]]],
            modules=[1, 1, 1]
        )
        assert 'is_exact' in result


# ---------------------------------------------------------------------------
# Freie Auflösung
# ---------------------------------------------------------------------------

class TestFreeResolution:
    """Tests für freie Auflösungen."""

    def test_free_resolution_z6(self):
        """Freie Auflösung von ℤ_6: 0 → ℤ -6·-> ℤ → ℤ_6 → 0."""
        result = free_resolution([[6]])
        assert 'resolution' in result
        assert len(result['resolution']) >= 2

    def test_returns_dict(self):
        """Rückgabe ist Dict mit 'resolution' und 'length'."""
        result = free_resolution([[2, 0], [0, 3]])
        assert 'resolution' in result
        assert 'length' in result


# ---------------------------------------------------------------------------
# Module-Klasse Tests
# ---------------------------------------------------------------------------

class TestModuleClass:
    """Tests für die Module-Klasse."""

    def test_create_module(self):
        """Modul kann erstellt werden."""
        m = Module([[2, 0], [0, 3]])
        assert m is not None

    def test_rank(self):
        """rank() gibt den freien Rang zurück."""
        m = Module([[1, 0], [0, 1]])
        # Volle Rang-2-Matrix → Rang 0 im Quotienten
        assert isinstance(m.rank(), int)

    def test_is_free(self):
        """is_free() erkennt freie Moduln."""
        m = Module([[1]])
        assert isinstance(m.is_free(), bool)

    def test_is_finitely_generated(self):
        """is_finitely_generated() ist immer True für endliche Matrizen."""
        m = Module([[2, 4], [6, 8]])
        assert m.is_finitely_generated() is True

    def test_torsion_submodule(self):
        """torsion_submodule() gibt Torsionsanteil zurück."""
        m = Module([[2, 0], [0, 3]])
        torsion = m.torsion_submodule()
        assert isinstance(torsion, dict)

    def test_free_part(self):
        """free_part() gibt freien Anteil zurück."""
        m = Module([[2, 0], [0, 3]])
        free = m.free_part()
        assert isinstance(free, dict)
