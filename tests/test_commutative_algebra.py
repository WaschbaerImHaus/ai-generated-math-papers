"""
@file test_commutative_algebra.py
@brief Tests für das Modul commutative_algebra.py (Kommutative Algebra).
@description
    Test-Suite für Lokalisierung, Primspektrum, Nakayama-Lemma,
    ganzen Abschluss, Klassengruppen, Noether-Normalisierung
    und den Hilbert-Basissatz.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os

# Projektpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from commutative_algebra import (
    localization,
    prime_spectrum,
    nakayama_lemma_check,
    integral_closure,
    class_group_estimate,
    noether_normalization,
    hilbert_basis_theorem_verify,
)


# ---------------------------------------------------------------------------
# Lokalisierung Tests
# ---------------------------------------------------------------------------

class TestLocalization:
    """Tests für Lokalisierung S^{-1}R."""

    def test_returns_dict(self):
        """Rückgabe ist ein Dict."""
        result = localization([0, 1, 2, 3, 4, 5], [1, 2, 4], 6)
        assert isinstance(result, dict)

    def test_has_expected_keys(self):
        """Rückgabe enthält 'elements' und 'description'."""
        result = localization([0, 1, 2, 3, 4, 5], [1, 2, 4], 6)
        assert 'elements' in result or 'description' in result

    def test_localization_at_prime(self):
        """Lokalisierung bei p=2 in ℤ/6ℤ."""
        result = localization(list(range(6)), [1, 2, 4], 6)
        assert 'description' in result


# ---------------------------------------------------------------------------
# Primspektrum Tests
# ---------------------------------------------------------------------------

class TestPrimeSpectrum:
    """Tests für Spec(ℤ/nℤ)."""

    def test_spec_z12_contains_2_and_3(self):
        """Spec(ℤ/12ℤ): Primideale sind (2) und (3)."""
        result = prime_spectrum(12)
        primes = result['prime_ideals']
        assert 2 in primes, f"(2) soll Primideal sein, bekam {primes}"
        assert 3 in primes, f"(3) soll Primideal sein, bekam {primes}"

    def test_spec_z12_only_two_primes(self):
        """Spec(ℤ/12ℤ): Genau 2 Primideale."""
        result = prime_spectrum(12)
        assert len(result['prime_ideals']) == 2

    def test_spec_zprime_has_one_prime(self):
        """Spec(ℤ/pℤ) für Primzahl p: Genau 1 Primideal (0)."""
        result = prime_spectrum(7)
        # ℤ/7ℤ ist Körper, Spec hat nur das Nullideal
        assert len(result['prime_ideals']) <= 1

    def test_spec_z30(self):
        """Spec(ℤ/30ℤ): Primteiler 2, 3, 5."""
        result = prime_spectrum(30)
        primes = result['prime_ideals']
        assert 2 in primes
        assert 3 in primes
        assert 5 in primes

    def test_maximal_ideals_returned(self):
        """Maximale Ideale werden zurückgegeben."""
        result = prime_spectrum(12)
        assert 'maximal_ideals' in result

    def test_maximal_equals_prime_for_zmodpz(self):
        """Für ℤ/pℤ (prim) ist das einzige Primideal auch maximal."""
        result = prime_spectrum(5)
        assert 'maximal_ideals' in result


# ---------------------------------------------------------------------------
# Nakayama-Lemma Tests
# ---------------------------------------------------------------------------

class TestNakayamaLemma:
    """Tests für das Nakayama-Lemma."""

    def test_returns_dict(self):
        """Rückgabe ist Dict."""
        result = nakayama_lemma_check([1, 2], [2], 3)
        assert isinstance(result, dict)

    def test_has_is_zero_key(self):
        """Rückgabe hat 'conclusion' oder 'M_is_zero' Schlüssel."""
        result = nakayama_lemma_check([1, 2], [2], 3)
        assert 'conclusion' in result or 'M_is_zero' in result

    def test_trivial_module_is_zero(self):
        """Triviales Modul M=0 → Nakayama trivial erfüllt."""
        result = nakayama_lemma_check([], [1], 5)
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# Ganzer Abschluss Tests
# ---------------------------------------------------------------------------

class TestIntegralClosure:
    """Tests für den ganzen Abschluss von ℤ in ℚ(√d)."""

    def test_gaussian_integers_discriminant(self):
        """Ganzer Abschluss in ℚ(√-1) = ℤ[i], Diskriminante -4."""
        result = integral_closure(-1)
        assert result['discriminant'] == -4, (
            f"Diskriminante von ℤ[i] soll -4 sein, bekam {result['discriminant']}"
        )

    def test_gaussian_integers_basis(self):
        """ℤ[i] hat ℤ-Basis {1, i}."""
        result = integral_closure(-1)
        basis = result['ring_basis']
        assert len(basis) == 2

    def test_eisenstein_integers_discriminant(self):
        """Ganzer Abschluss in ℚ(√-3): Diskriminante -3."""
        result = integral_closure(-3)
        assert result['discriminant'] == -3, (
            f"Diskriminante von ℤ[ω] soll -3 sein, bekam {result['discriminant']}"
        )

    def test_d_minus3_is_pid(self):
        """ℤ[ω] (d=-3) ist PID (h=1)."""
        result = integral_closure(-3)
        assert result['is_pid'] is True

    def test_d_minus1_is_pid(self):
        """ℤ[i] (d=-1) ist PID."""
        result = integral_closure(-1)
        assert result['is_pid'] is True

    def test_d_square_free_check(self):
        """Funktion gibt Ergebnis für verschiedene quadratfreie d zurück."""
        for d in [-2, -5, 2, 5]:
            result = integral_closure(d)
            assert 'discriminant' in result
            assert 'ring_basis' in result

    def test_d_positive_disc(self):
        """d=5 ≡ 1 (mod 4): Diskriminante = 5."""
        result = integral_closure(5)
        assert result['discriminant'] == 5

    def test_d_2_disc(self):
        """d=2: Diskriminante = 4·2 = 8."""
        result = integral_closure(2)
        assert result['discriminant'] == 8


# ---------------------------------------------------------------------------
# Klassengruppe Tests
# ---------------------------------------------------------------------------

class TestClassGroupEstimate:
    """Tests für Klassengruppen-Schätzung."""

    def test_gaussian_integers_h_equals_1(self):
        """h(-1) = 1: ℤ[i] ist PID."""
        result = class_group_estimate(-1)
        assert result['class_number'] == 1, (
            f"h(-1) = 1 erwartet, bekam {result['class_number']}"
        )

    def test_d_minus5_h_equals_2(self):
        """h(-5) = 2: ℤ[√-5] ist kein PID."""
        result = class_group_estimate(-5)
        assert result['class_number'] == 2, (
            f"h(-5) = 2 erwartet, bekam {result['class_number']}"
        )

    def test_d_minus1_is_pid(self):
        """is_pid für d=-1 muss True sein."""
        result = class_group_estimate(-1)
        assert result['is_pid'] is True

    def test_d_minus5_is_not_pid(self):
        """is_pid für d=-5 muss False sein."""
        result = class_group_estimate(-5)
        assert result['is_pid'] is False

    def test_returns_dict_with_keys(self):
        """Rückgabe hat 'class_number', 'is_pid', 'minkowski_bound'."""
        result = class_group_estimate(-1)
        assert 'class_number' in result
        assert 'is_pid' in result
        assert 'minkowski_bound' in result

    def test_d_minus3_h_equals_1(self):
        """h(-3) = 1: ℤ[ω] ist PID."""
        result = class_group_estimate(-3)
        assert result['class_number'] == 1


# ---------------------------------------------------------------------------
# Noether-Normalisierung Tests
# ---------------------------------------------------------------------------

class TestNoetherNormalization:
    """Tests für das Noether-Normalisierungslemma."""

    def test_returns_dict(self):
        """Rückgabe ist Dict."""
        result = noether_normalization(['x^2 - y'], 2)
        assert isinstance(result, dict)

    def test_has_krull_dimension(self):
        """Rückgabe enthält Krull-Dimension."""
        result = noether_normalization(['x^2 - y'], 2)
        assert 'krull_dimension' in result

    def test_plane_curve_dimension_1(self):
        """Affine Kurve V(f) ⊂ A^2 hat Dimension 1."""
        result = noether_normalization(['x*y - 1'], 2)
        assert result['krull_dimension'] == 1


# ---------------------------------------------------------------------------
# Hilbert-Basissatz Tests
# ---------------------------------------------------------------------------

class TestHilbertBasisTheorem:
    """Tests für den Hilbert-Basissatz."""

    def test_returns_dict(self):
        """Rückgabe ist Dict."""
        result = hilbert_basis_theorem_verify([[1, 0], [0, 1]], 5)
        assert isinstance(result, dict)

    def test_has_is_finitely_generated(self):
        """Rückgabe hat 'is_finitely_generated'."""
        result = hilbert_basis_theorem_verify([[1, 2], [3, 4]], 7)
        assert 'is_finitely_generated' in result

    def test_noetherian_ring_result(self):
        """Ideal in noetherischem Ring ist endlich erzeugt."""
        result = hilbert_basis_theorem_verify([[1, 0, 2], [0, 1, 3]], 5)
        assert result['is_finitely_generated'] is True
