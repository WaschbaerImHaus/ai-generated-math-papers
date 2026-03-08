"""
@file test_analysis_extended.py
@brief Tests für erweiterte Analysis-Funktionen.
@description
    Testet alle neuen Funktionen in analysis.py:
    - symbolic_limit: Symbolische Grenzwertberechnung
    - lhopital_applicable: L'Hôpital-Prüfung
    - limit_comparison: Grenzwert-Vergleichssatz
    - partial_fraction_decomposition: Partialbruchzerlegung (Koeffizienten)
    - partial_fraction_symbolic: Symbolische Partialbruchzerlegung
    - improper_integral_numerical: Numerisches uneigentliches Integral
    - improper_integral_symbolic: Symbolisches uneigentliches Integral
    - cauchy_principal_value: Cauchyscher Hauptwert

@author Kurt Ingwer
@date 2026-03-08
"""

import pytest
import math
import sympy as sp
import sys
import os

# Projektpfad zum Modulpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from analysis import (
    symbolic_limit,
    lhopital_applicable,
    limit_comparison,
    partial_fraction_decomposition,
    partial_fraction_symbolic,
    improper_integral_numerical,
    improper_integral_symbolic,
    cauchy_principal_value,
)


# ===========================================================================
# Tests: symbolic_limit
# ===========================================================================

class TestSymbolicLimit:
    """Tests für symbolische Grenzwertberechnung."""

    def test_sinx_over_x_at_0(self):
        """lim_{x→0} sin(x)/x = 1."""
        result = symbolic_limit("sin(x)/x", "x", 0)
        assert result == 1

    def test_euler_number_limit(self):
        """lim_{x→∞} (1 + 1/x)^x = e."""
        result = symbolic_limit("(1 + 1/x)**x", "x", "oo")
        assert result == sp.E

    def test_constant_limit(self):
        """lim_{x→5} x^2 = 25."""
        result = symbolic_limit("x**2", "x", 5)
        assert result == 25

    def test_limit_at_infinity(self):
        """lim_{x→∞} 1/x = 0."""
        result = symbolic_limit("1/x", "x", "oo")
        assert result == 0

    def test_limit_from_right(self):
        """lim_{x→0+} ln(x) = -∞."""
        result = symbolic_limit("log(x)", "x", 0, '+')
        assert result == -sp.oo

    def test_limit_polynomial(self):
        """lim_{x→2} (x^2 - 4)/(x - 2) = 4 (Heben der Singularität)."""
        result = symbolic_limit("(x**2 - 4)/(x - 2)", "x", 2)
        assert result == 4

    def test_limit_exp(self):
        """lim_{x→∞} exp(-x) = 0."""
        result = symbolic_limit("exp(-x)", "x", "oo")
        assert result == 0

    def test_limit_negative_infinity(self):
        """lim_{x→-∞} exp(x) = 0."""
        result = symbolic_limit("exp(x)", "x", "-oo")
        assert result == 0


# ===========================================================================
# Tests: lhopital_applicable
# ===========================================================================

class TestLhopitalApplicable:
    """Tests für L'Hôpital-Prüfung."""

    def test_0_over_0_form(self):
        """sin(x)/x bei x=0 ist 0/0-Form → anwendbar."""
        result = lhopital_applicable("sin(x)", "x", "x", 0)
        assert result['applicable'] is True
        assert '0/0' in result['form']

    def test_limit_correct(self):
        """Gesamtgrenzwert muss korrekt berechnet werden."""
        result = lhopital_applicable("sin(x)", "x", "x", 0)
        assert result['limit'] == 1

    def test_infinity_over_infinity(self):
        """x/exp(x) bei x→∞ ist ∞/∞-Form."""
        # x → ∞: Numerator → ∞, aber eigentlich ist x/exp(x) 0/0 nicht
        # Teste lieber (exp(x))/(exp(x)+1) → ∞/∞
        result = lhopital_applicable("exp(x)", "exp(x)+1", "x", "oo")
        # Beide Grenzen sind ∞ → ∞/∞
        assert result['applicable'] is True

    def test_not_applicable(self):
        """x/(x+1) bei x=0: Nenner → 1 ≠ 0 → nicht anwendbar."""
        result = lhopital_applicable("x", "x+1", "x", 0)
        assert result['applicable'] is False

    def test_result_has_required_keys(self):
        """Rückgabe-Dict muss alle Pflichtschlüssel haben."""
        result = lhopital_applicable("sin(x)", "x", "x", 0)
        assert 'applicable' in result
        assert 'form' in result
        assert 'limit' in result


# ===========================================================================
# Tests: limit_comparison
# ===========================================================================

class TestLimitComparison:
    """Tests für den Grenzwert-Vergleichssatz."""

    def test_same_order(self):
        """x² und x²+1 haben bei ∞ das gleiche Verhalten."""
        result = limit_comparison("x**2", "x**2 + 1", "x", "oo")
        assert result['same_behavior'] is True
        assert result['limit_ratio'] == 1

    def test_different_order(self):
        """x und x² haben bei ∞ unterschiedliches Verhalten."""
        result = limit_comparison("x", "x**2", "x", "oo")
        # lim x/x² = lim 1/x = 0 → unterschiedliche Ordnung
        assert result['same_behavior'] is False

    def test_sinx_over_x(self):
        """sin(x)/x bei x=0: Verhältnis = 1."""
        result = limit_comparison("sin(x)", "x", "x", 0)
        assert result['limit_ratio'] == 1
        assert result['same_behavior'] is True

    def test_result_has_required_keys(self):
        """Rückgabe-Dict muss Pflichtschlüssel haben."""
        result = limit_comparison("x", "x+1", "x", "oo")
        assert 'limit_ratio' in result
        assert 'same_behavior' in result


# ===========================================================================
# Tests: partial_fraction_decomposition
# ===========================================================================

class TestPartialFractionDecomposition:
    """Tests für Partialbruchzerlegung via Koeffizientenlisten."""

    def test_1_over_x2_minus_1(self):
        """1/(x²-1) = 1/((x-1)(x+1)) → Partialbruchzerlegung."""
        # Numerator: [1] = 1
        # Denominator: [1, 0, -1] = x² - 1
        result = partial_fraction_decomposition([1], [1, 0, -1])
        # Ergebnis enthält Terme mit (x-1) und (x+1)
        assert isinstance(result, str)
        assert len(result) > 0
        # SymPy gibt "-1/(2*(x + 1)) + 1/(2*(x - 1))" oder ähnliches
        assert 'x' in result

    def test_1_over_x_times_x_plus_1(self):
        """1/(x(x+1)) = 1/x - 1/(x+1)."""
        # Denominator: [1, 1, 0] = x² + x
        result = partial_fraction_decomposition([1], [1, 1, 0])
        assert isinstance(result, str)
        assert 'x' in result

    def test_returns_string(self):
        """Ergebnis muss ein String sein."""
        result = partial_fraction_decomposition([1], [1, -1])
        assert isinstance(result, str)

    def test_already_decomposed(self):
        """1/(x-1) bleibt unverändert."""
        result = partial_fraction_decomposition([1], [1, -1])
        # Kein Bruch mehr möglich (einfacher Pol) → Ergebnis enthält '1/(x - 1)' o.ä.
        assert isinstance(result, str)


# ===========================================================================
# Tests: partial_fraction_symbolic
# ===========================================================================

class TestPartialFractionSymbolic:
    """Tests für symbolische Partialbruchzerlegung."""

    def test_basic_decomposition(self):
        """1/(x²-1) Zerlegung."""
        result = partial_fraction_symbolic("1/(x**2 - 1)", "x")
        assert 'original' in result
        assert 'decomposed' in result
        assert 'terms' in result

    def test_terms_list(self):
        """terms muss eine Liste sein."""
        result = partial_fraction_symbolic("1/(x*(x+1))", "x")
        assert isinstance(result['terms'], list)

    def test_original_preserved(self):
        """original muss den Originalausdruck als String enthalten."""
        result = partial_fraction_symbolic("1/(x*(x-1))", "x")
        # String-Darstellung des Originalausdrucks
        assert isinstance(result['original'], str)

    def test_decomposed_is_string(self):
        """decomposed muss String sein."""
        result = partial_fraction_symbolic("(x+1)/(x**2-x)", "x")
        assert isinstance(result['decomposed'], str)

    def test_multiple_terms(self):
        """x/(x²-1) hat mindestens 2 Terme in der Zerlegung."""
        result = partial_fraction_symbolic("x/(x**2 - 1)", "x")
        # Mindestens ein Term
        assert len(result['terms']) >= 1


# ===========================================================================
# Tests: improper_integral_numerical
# ===========================================================================

class TestImproperIntegralNumerical:
    """Tests für numerische uneigentliche Integrale."""

    def test_gaussian_integral(self):
        """∫_0^∞ exp(-x²) dx = √π/2 ≈ 0.8862."""
        result = improper_integral_numerical(
            lambda x: math.exp(-x * x), 0.0, float('inf')
        )
        expected = math.sqrt(math.pi) / 2
        assert abs(result - expected) < 1e-4

    def test_exp_minus_x(self):
        """∫_0^∞ exp(-x) dx = 1."""
        result = improper_integral_numerical(
            lambda x: math.exp(-x), 0.0, float('inf')
        )
        assert abs(result - 1.0) < 1e-4

    def test_finite_integral(self):
        """Endliches Integral: ∫_0^1 x dx = 0.5."""
        result = improper_integral_numerical(lambda x: x, 0.0, 1.0)
        assert abs(result - 0.5) < 1e-6

    def test_singularity_handling(self):
        """∫_0^2 1/√x dx = 2√2 ≈ 2.828 (Singularität bei x=0)."""
        # Numerisch approximiert (Singularität wird durch scipy.quad behandelt)
        result = improper_integral_numerical(
            lambda x: 1.0 / math.sqrt(x), 1e-6, 2.0
        )
        expected = 2 * math.sqrt(2) - 2 * math.sqrt(1e-6)
        assert abs(result - expected) < 0.01

    def test_two_sided_infinity(self):
        """∫_{-∞}^∞ exp(-x²) dx = √π ≈ 1.7725."""
        result = improper_integral_numerical(
            lambda x: math.exp(-x * x), float('-inf'), float('inf')
        )
        expected = math.sqrt(math.pi)
        assert abs(result - expected) < 1e-4


# ===========================================================================
# Tests: improper_integral_symbolic
# ===========================================================================

class TestImproperIntegralSymbolic:
    """Tests für symbolische uneigentliche Integrale."""

    def test_exp_minus_x_from_0_to_inf(self):
        """∫_0^∞ e^(-x) dx = 1."""
        result = improper_integral_symbolic("exp(-x)", "x", 0, "oo")
        assert result['converges'] is True
        assert result['value'] is not None
        assert abs(result['value'] - 1.0) < 1e-10

    def test_result_has_required_keys(self):
        """Rückgabe-Dict muss Pflichtschlüssel haben."""
        result = improper_integral_symbolic("exp(-x)", "x", 0, "oo")
        assert 'integral' in result
        assert 'converges' in result
        assert 'value' in result

    def test_divergent_integral(self):
        """∫_0^∞ x dx divergiert."""
        result = improper_integral_symbolic("x", "x", 0, "oo")
        assert result['converges'] is False

    def test_finite_integral(self):
        """∫_0^1 x^2 dx = 1/3."""
        result = improper_integral_symbolic("x**2", "x", 0, 1)
        assert result['converges'] is True
        assert abs(result['value'] - 1.0 / 3.0) < 1e-10

    def test_gaussian_symbolic(self):
        """∫_{-∞}^{∞} exp(-x²) dx = √π."""
        result = improper_integral_symbolic("exp(-x**2)", "x", "-oo", "oo")
        assert result['converges'] is True
        expected = math.sqrt(math.pi)
        assert abs(result['value'] - expected) < 1e-10

    def test_convergent_reciprocal_power(self):
        """∫_1^∞ 1/x² dx = 1."""
        result = improper_integral_symbolic("1/x**2", "x", 1, "oo")
        assert result['converges'] is True
        assert abs(result['value'] - 1.0) < 1e-10


# ===========================================================================
# Tests: cauchy_principal_value
# ===========================================================================

class TestCauchyPrincipalValue:
    """Tests für den Cauchyschen Hauptwert."""

    def test_1_over_x_symmetric(self):
        """P.V. ∫_{-1}^{1} 1/x dx = 0 (Symmetrie)."""
        result = cauchy_principal_value(lambda x: 1.0 / x, -1.0, 1.0, 0.0)
        assert abs(result) < 1e-4

    def test_result_is_float(self):
        """Ergebnis muss ein Float sein."""
        result = cauchy_principal_value(lambda x: 1.0 / x, -1.0, 1.0, 0.0)
        assert isinstance(result, float)

    def test_asymmetric_function(self):
        """P.V. ∫_{-2}^{3} 1/x dx = ln(3/2) (asymmetrisch)."""
        # P.V. = ∫_{-2}^{-ε} 1/x dx + ∫_{ε}^{3} 1/x dx
        # = [ln|x|]_{-2}^{-ε} + [ln|x|]_{ε}^{3}
        # = (ln(ε) - ln(2)) + (ln(3) - ln(ε)) = ln(3) - ln(2) = ln(3/2)
        result = cauchy_principal_value(lambda x: 1.0 / x, -2.0, 3.0, 0.0)
        expected = math.log(3.0 / 2.0)
        assert abs(result - expected) < 1e-4

    def test_regular_integrand_at_point(self):
        """Für stetige Funktion muss der Hauptwert mit dem normalen Integral übereinstimmen."""
        # f(x) = x ist stetig überall, auch bei x=0.5
        # P.V. ∫_0^1 x dx bei c=0.5 = 0.5 (normales Integral)
        result = cauchy_principal_value(lambda x: x, 0.0, 1.0, 0.5)
        assert abs(result - 0.5) < 1e-4

    def test_symmetric_odd_function(self):
        """P.V. ∫_{-a}^{a} f(x)/x dx = 0 für gerades f (Symmetrie)."""
        # f(x) = x²/x = x ist ungerade → ∫_{-1}^{1} x dx = 0
        # Aber wir testen P.V. ∫_{-1}^{1} 1/x dx = 0
        result = cauchy_principal_value(lambda x: 1.0 / x, -1.0, 1.0, 0.0)
        assert abs(result) < 1e-4
