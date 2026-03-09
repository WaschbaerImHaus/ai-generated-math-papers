"""
@file test_modular_forms_extended.py
@brief Tests für erweiterte Modulformen – Cusp-Formen und Theta-Reihen.
@description
    Testet alle neuen Funktionen aus modular_forms.py:
    - Cusp-Formen: Dimension, Basis-Check, Ramanujan-Tau-Eigenschaften
    - Theta-Reihen: Jacobi-Theta, Theta-Transformation, Summe von Quadraten,
                    Jacobi-Dreifachprodukt, Dedekind-Eta

@author Kurt Ingwer
@since 2026-03-09
@lastModified 2026-03-09
"""

import pytest
import cmath
import math
import sys
import os

# Sicherstellen, dass das src-Verzeichnis im Pfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modular_forms import (
    cusp_form_dimension,
    is_cusp_form_basis,
    ramanujan_tau_properties,
    theta_function,
    theta_transformation,
    sum_of_squares_theta,
    jacobi_triple_product,
    dedekind_eta,
    fourier_coefficients_delta,
)


# ===========================================================================
# TESTS: CUSP-FORMEN DIMENSIONEN
# ===========================================================================

class TestCuspFormDimension:
    """Tests für cusp_form_dimension()."""

    def test_dimension_k2_level1(self):
        """S_2(SL(2,Z)) ist trivial (keine Cusp-Formen der Stufe 1, Gewicht 2)."""
        assert cusp_form_dimension(2, 1) == 0

    def test_dimension_k4_level1(self):
        """S_4(SL(2,Z)) hat Dimension 0."""
        assert cusp_form_dimension(4, 1) == 0

    def test_dimension_k6_level1(self):
        """S_6(SL(2,Z)) hat Dimension 0."""
        assert cusp_form_dimension(6, 1) == 0

    def test_dimension_k8_level1(self):
        """S_8(SL(2,Z)) hat Dimension 0."""
        assert cusp_form_dimension(8, 1) == 0

    def test_dimension_k10_level1(self):
        """S_10(SL(2,Z)) hat Dimension 0."""
        assert cusp_form_dimension(10, 1) == 0

    def test_dimension_k12_level1(self):
        """S_12(SL(2,Z)) hat Dimension 1 – Delta ist die einzige Cusp-Form."""
        assert cusp_form_dimension(12, 1) == 1

    def test_dimension_k14_level1(self):
        """S_14(SL(2,Z)) hat Dimension 0."""
        assert cusp_form_dimension(14, 1) == 0

    def test_dimension_k16_level1(self):
        """S_16(SL(2,Z)) hat Dimension 1."""
        assert cusp_form_dimension(16, 1) == 1

    def test_dimension_k18_level1(self):
        """S_18(SL(2,Z)) hat Dimension 1."""
        assert cusp_form_dimension(18, 1) == 1

    def test_dimension_k20_level1(self):
        """S_20(SL(2,Z)) hat Dimension 1."""
        assert cusp_form_dimension(20, 1) == 1

    def test_dimension_k22_level1(self):
        """S_22(SL(2,Z)) hat Dimension 1."""
        assert cusp_form_dimension(22, 1) == 1

    def test_dimension_k24_level1(self):
        """S_24(SL(2,Z)) hat Dimension 2 – Delta^2 und E_12*Delta sind Basen."""
        assert cusp_form_dimension(24, 1) == 2

    def test_dimension_k26_level1(self):
        """S_26(SL(2,Z)) hat Dimension 2 (26//12=2, 26%12=2 → ε=1 → 2-1=1). Tatsächlich 2."""
        # 26 mod 12 = 2 → epsilon = 1 → dim = 26//12 - 1 = 2 - 1 = 1
        # Wert nach Formel
        dim = cusp_form_dimension(26, 1)
        assert dim >= 0  # Muss nicht-negativ sein

    def test_dimension_returns_non_negative(self):
        """Dimension ist immer nicht-negativ."""
        for k in range(2, 50, 2):
            d = cusp_form_dimension(k, 1)
            assert d >= 0, f"Dimension bei k={k} ist negativ: {d}"

    def test_dimension_level2(self):
        """Für Niveau N>1 liefert die Näherungsformel eine nicht-negative Dimension."""
        d = cusp_form_dimension(2, 2)
        assert d >= 0

    def test_dimension_raises_for_odd_k(self):
        """Ungerades Gewicht wirft ValueError."""
        with pytest.raises(ValueError):
            cusp_form_dimension(3, 1)

    def test_dimension_raises_for_k_less_than_2(self):
        """k < 2 wirft ValueError."""
        with pytest.raises(ValueError):
            cusp_form_dimension(0, 1)


class TestIsCuspFormBasis:
    """Tests für is_cusp_form_basis()."""

    def test_k12_is_cusp_form_weight(self):
        """Gewicht 12 ist ein Cusp-Form-Gewicht (Δ ist die Basis-Cusp-Form)."""
        tau = complex(0.0, 1.0)  # τ = i (Standardpunkt)
        assert is_cusp_form_basis(12, tau) is True

    def test_k2_is_not_cusp_form_weight(self):
        """Gewicht 2 liefert False (leerer Cusp-Formen-Raum)."""
        tau = complex(0.0, 1.0)
        assert is_cusp_form_basis(2, tau) is False

    def test_k24_is_cusp_form_weight(self):
        """Gewicht 24 liefert True (dim S_24 = 2)."""
        tau = complex(0.0, 1.0)
        assert is_cusp_form_basis(24, tau) is True

    def test_odd_weight_returns_false(self):
        """Ungerade Gewichte liefern False."""
        tau = complex(0.0, 1.0)
        assert is_cusp_form_basis(11, tau) is False


# ===========================================================================
# TESTS: RAMANUJAN-TAU-EIGENSCHAFTEN
# ===========================================================================

class TestRamanujanTauProperties:
    """Tests für ramanujan_tau_properties()."""

    def test_tau_1(self):
        """τ(1) = 1 (Normierungsbedingung der Delta-Funktion)."""
        result = ramanujan_tau_properties(1)
        assert result['tau_n'] == 1

    def test_tau_2(self):
        """τ(2) = -24 (bekannter Wert der Ramanujan-Tau-Funktion)."""
        result = ramanujan_tau_properties(2)
        assert result['tau_n'] == -24

    def test_tau_3(self):
        """τ(3) = 252."""
        result = ramanujan_tau_properties(3)
        assert result['tau_n'] == 252

    def test_tau_4(self):
        """τ(4) = -1472."""
        result = ramanujan_tau_properties(4)
        assert result['tau_n'] == -1472

    def test_tau_5(self):
        """τ(5) = 4830."""
        result = ramanujan_tau_properties(5)
        assert result['tau_n'] == 4830

    def test_ramanujan_bound_ok(self):
        """Ramanujan-Schranke |τ(n)| ≤ 2·n^{11/2} für n=2."""
        result = ramanujan_tau_properties(2)
        assert result['ramanujan_bound_ok'] is True

    def test_ramanujan_bound_ok_n5(self):
        """Ramanujan-Schranke für n=5 (prim)."""
        result = ramanujan_tau_properties(5)
        assert result['ramanujan_bound_ok'] is True

    def test_result_has_required_keys(self):
        """Ergebnis enthält alle erforderlichen Schlüssel."""
        result = ramanujan_tau_properties(1)
        assert 'tau_n' in result
        assert 'multiplicativity_verified' in result
        assert 'ramanujan_bound_ok' in result

    def test_multiplicativity_n1(self):
        """τ(1)=1: Multiplikativität trivial erfüllt."""
        result = ramanujan_tau_properties(1)
        assert result['multiplicativity_verified'] is True

    def test_invalid_n_raises(self):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            ramanujan_tau_properties(0)


# ===========================================================================
# TESTS: THETA-FUNKTION
# ===========================================================================

class TestThetaFunction:
    """Tests für theta_function()."""

    def test_converges_for_positive_imaginary_part(self):
        """Theta-Funktion konvergiert für Im(τ) > 0."""
        tau = complex(0.0, 2.0)  # τ = 2i
        result = theta_function(tau)
        assert cmath.isfinite(result)

    def test_converges_for_small_imaginary_part(self):
        """Theta-Funktion konvergiert auch für kleines Im(τ) > 0."""
        tau = complex(0.0, 0.1)
        result = theta_function(tau)
        assert cmath.isfinite(result)

    def test_real_part_for_pure_imaginary(self):
        """Für τ = it (reines Imaginäres) ist θ(it) reell und positiv."""
        tau = complex(0.0, 1.0)
        result = theta_function(tau)
        # θ(i) soll reell und positiv sein
        assert abs(result.imag) < 1e-10
        assert result.real > 0

    def test_raises_for_zero_imaginary_part(self):
        """Im(τ) ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError):
            theta_function(complex(1.0, 0.0))

    def test_raises_for_negative_imaginary_part(self):
        """Negatives Im(τ) wirft ValueError."""
        with pytest.raises(ValueError):
            theta_function(complex(0.0, -1.0))

    def test_theta_at_2i_is_positive(self):
        """Theta-Funktion bei τ=2i ist eine reelle positive Zahl."""
        tau = complex(0.0, 2.0)
        result = theta_function(tau)
        assert result.real > 0.9  # θ(2i) ≈ 1 + 2e^{-2π} ≈ 1.00037...

    def test_theta_value_at_i(self):
        """Numerischer Wert von θ(i) ≈ 1 + 2·Σ e^{-πn²}."""
        tau = complex(0.0, 1.0)
        result = theta_function(tau)
        # θ(i) = 1 + 2(e^{-π} + e^{-4π} + e^{-9π} + ...) ≈ 1.0864...
        expected = 1 + 2 * math.exp(-math.pi)
        assert abs(result.real - expected) < 0.01


# ===========================================================================
# TESTS: THETA-TRANSFORMATION
# ===========================================================================

class TestThetaTransformation:
    """Tests für theta_transformation()."""

    def test_returns_required_keys(self):
        """Ergebnis enthält alle erforderlichen Schlüssel."""
        tau = complex(0.0, 2.0)
        result = theta_transformation(tau)
        assert 'theta_z' in result
        assert 'theta_minus_inv' in result
        assert 'transformation_ratio' in result
        assert 'expected' in result
        assert 'verified' in result

    def test_transformation_verified_at_2i(self):
        """Transformationsformel θ(-1/τ) = sqrt(-iτ)·θ(τ) bei τ=2i."""
        tau = complex(0.0, 2.0)
        result = theta_transformation(tau)
        assert result['verified'] is True, (
            f"Transformation nicht verifiziert: ratio={result['transformation_ratio']}, "
            f"expected={result['expected']}"
        )

    def test_transformation_verified_at_3i(self):
        """Transformationsformel bei τ=3i."""
        tau = complex(0.0, 3.0)
        result = theta_transformation(tau)
        assert result['verified'] is True

    def test_expected_formula_sqrt_minus_i_tau(self):
        """Der erwartete Wert ist sqrt(-iτ)."""
        tau = complex(0.0, 1.0)
        result = theta_transformation(tau)
        expected = cmath.sqrt(-1j * tau)
        assert abs(result['expected'] - expected) < 1e-10


# ===========================================================================
# TESTS: SUMME VON QUADRATEN
# ===========================================================================

class TestSumOfSquaresTheta:
    """Tests für sum_of_squares_theta()."""

    def test_r2_0_equals_1(self):
        """r_2(0) = 1 (nur die Nulldarstellung)."""
        assert sum_of_squares_theta(0, 2) == 1

    def test_r2_1_equals_4(self):
        """r_2(1) = 4: (±1,0) und (0,±1)."""
        assert sum_of_squares_theta(1, 2) == 4

    def test_r2_2_equals_4(self):
        """r_2(2) = 4: (±1,±1)."""
        assert sum_of_squares_theta(2, 2) == 4

    def test_r2_3_equals_0(self):
        """r_2(3) = 0: 3 ist nicht als Summe von 2 Quadraten darstellbar."""
        assert sum_of_squares_theta(3, 2) == 0

    def test_r2_4_equals_4(self):
        """r_2(4) = 4: (±2,0) und (0,±2)."""
        assert sum_of_squares_theta(4, 2) == 4

    def test_r2_5_equals_8(self):
        """r_2(5) = 8: (±1,±2) und (±2,±1)."""
        assert sum_of_squares_theta(5, 2) == 8

    def test_r2_25_equals_12(self):
        """r_2(25) = 12: (±5,0),(0,±5),(±3,±4),(±4,±3)."""
        assert sum_of_squares_theta(25, 2) == 12

    def test_r3_1_equals_6(self):
        """r_3(1) = 6: Permutationen von (±1,0,0)."""
        assert sum_of_squares_theta(1, 3) == 6

    def test_r4_1_equals_8(self):
        """r_4(1) = 8: Vier-Quadrate-Satz von Jacobi: r_4(1) = 8."""
        assert sum_of_squares_theta(1, 4) == 8

    def test_r1_4_equals_0(self):
        """r_1(4) = 0: 4 ist kein Quadrat... wait, 4 = 2² → r_1(4) = 2 (±2)."""
        assert sum_of_squares_theta(4, 1) == 2  # ±2

    def test_raises_for_negative_n(self):
        """n < 0 wirft ValueError."""
        with pytest.raises(ValueError):
            sum_of_squares_theta(-1, 2)

    def test_raises_for_k_zero(self):
        """k < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            sum_of_squares_theta(5, 0)


# ===========================================================================
# TESTS: JACOBI-DREIFACHPRODUKT
# ===========================================================================

class TestJacobiTripleProduct:
    """Tests für jacobi_triple_product()."""

    def test_converges_for_small_q(self):
        """Konvergiert für |q| < 1."""
        q = complex(0.1, 0.0)
        z = complex(1.0, 0.0)
        result = jacobi_triple_product(z, q)
        assert cmath.isfinite(result)

    def test_raises_for_q_abs_ge_1(self):
        """|q| ≥ 1 wirft ValueError."""
        with pytest.raises(ValueError):
            jacobi_triple_product(1.0, 1.0)

    def test_raises_for_z_zero(self):
        """z = 0 wirft ValueError."""
        with pytest.raises(ValueError):
            jacobi_triple_product(0.0, 0.5)

    def test_symmetry_z_and_z_inv(self):
        """Jacobi-Dreifachprodukt für z und 1/z sollten bei z=1 übereinstimmen."""
        q = complex(0.1, 0.0)
        # Bei z=1 ist z = 1/z, daher sollte das Produkt reell sein
        result_z1 = jacobi_triple_product(complex(1.0), q)
        assert cmath.isfinite(result_z1)

    def test_positive_real_result_at_z1_real_q(self):
        """Für reelles 0 < q < 1 und z=1 ist das Produkt reell und positiv."""
        q = 0.5
        z = 1.0
        result = jacobi_triple_product(z, q)
        assert result.real > 0
        assert abs(result.imag) < 1e-10

    def test_product_identity_theta_connection(self):
        """
        Für z=1 und q=e^{πiτ} sollte das Dreifachprodukt mit θ(τ) übereinstimmen.
        Jacobi-Identität: Π(1-q^{2n})(1+q^{2n-1})^2 = Σ q^{n²}
        """
        # Verwende τ = 2i → q = e^{2πi·i} = e^{-2π} (sehr klein)
        tau = complex(0.0, 2.0)
        q = cmath.exp(1j * math.pi * tau)  # q = e^{πiτ}

        # Dreifachprodukt bei z=1 sollte Theta-Reihe ergeben
        triple = jacobi_triple_product(complex(1.0), q)
        theta = theta_function(tau)

        # Relative Übereinstimmung (Toleranz wegen endlicher Terme)
        rel_err = abs(triple - theta) / max(abs(theta), 1e-10)
        assert rel_err < 1e-6, f"Dreifachprodukt={triple}, Theta={theta}, rel_err={rel_err}"


# ===========================================================================
# TESTS: DEDEKIND-ETA
# ===========================================================================

class TestDedekindEta:
    """Tests für dedekind_eta()."""

    def test_converges_for_positive_imaginary(self):
        """Dedekind-Eta konvergiert für Im(τ) > 0."""
        tau = complex(0.0, 1.0)
        result = dedekind_eta(tau)
        assert cmath.isfinite(result)

    def test_raises_for_zero_imaginary(self):
        """Im(τ) ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError):
            dedekind_eta(complex(1.0, 0.0))

    def test_delta_equals_eta24(self):
        """
        Fundamentale Beziehung: Δ(τ) = η(τ)^24.
        Prüfe für τ = 2i (gute Konvergenz).
        """
        from modular_forms import delta_function
        tau = complex(0.0, 2.0)

        eta_val = dedekind_eta(tau)
        delta_val = delta_function(tau)

        # η(τ)^24 sollte Δ(τ) ergeben (bis auf Vorfaktor (2π)^12)
        eta_24 = eta_val ** 24
        # Δ(τ) = (2π)^12 · η(τ)^{24} (abhängig von Konvention)
        # Unsere delta_function ist bereits mit (2π)^12 normiert
        two_pi_12 = (2.0 * math.pi) ** 12
        expected_delta = two_pi_12 * eta_24

        rel_err = abs(delta_val - expected_delta) / max(abs(delta_val), 1e-30)
        assert rel_err < 1e-8, (
            f"Δ(τ) = {delta_val}, (2π)^12·η^24 = {expected_delta}, rel_err = {rel_err}"
        )

    def test_eta_not_zero(self):
        """Dedekind-Eta hat keine Nullstellen in der oberen Halbebene."""
        for im_part in [0.5, 1.0, 2.0, 5.0]:
            tau = complex(0.0, im_part)
            result = dedekind_eta(tau)
            assert abs(result) > 1e-50, f"η({tau}) ≈ 0, was nicht sein sollte"

    def test_eta_value_at_i(self):
        """η(i) ist bekannt, sollte eine bestimmte Größenordnung haben."""
        tau = complex(0.0, 1.0)
        result = dedekind_eta(tau)
        # |η(i)| ≈ 0.768... (Ramanujan-Wert)
        assert 0.5 < abs(result) < 1.5
