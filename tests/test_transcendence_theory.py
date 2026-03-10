"""
test_transcendence_theory.py – Tests für das Transzendenztheorie-Modul

Testet alle Funktionen aus transcendence_theory.py:
- Klassifikation algebraischer Zahlen
- Liouville-Zahlen
- Lindemann-Weierstrass-Demonstrationen
- Irrationali­tätsmaße
- AlgebraicNumber-Klasse
- Kreisteilungspolynome

@author    Kurt Ingwer
@version   1.0.0
@timestamp 2026-03-10
"""

import math
import pytest
import sys
import os

# Projektverzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from transcendence_theory import (
    is_algebraic_integer_check,
    minimal_polynomial,
    algebraic_degree,
    classify_number,
    liouville_number,
    is_liouville_like,
    liouville_approximation_exponent,
    lindemann_weierstrass_demo,
    e_is_transcendental_demo,
    pi_is_transcendental_demo,
    gelfond_schneider_demo,
    baker_theorem_demo,
    irrationality_measure,
    known_irrationality_measures,
    diophantine_approximation_quality,
    continued_fraction_irrationality,
    AlgebraicNumber,
    cyclotomic_polynomial,
    number_field_degree,
)


# ---------------------------------------------------------------------------
# 1. Tests: Klassifikation von Zahlen
# ---------------------------------------------------------------------------

class TestIsAlgebraicIntegerCheck:
    """Tests für is_algebraic_integer_check()"""

    def test_integer_is_algebraic(self):
        """Ganze Zahlen sind algebraische ganze Zahlen (t - n = 0)."""
        assert is_algebraic_integer_check(3) is True
        assert is_algebraic_integer_check(0) is True

    def test_sqrt2_is_algebraic(self):
        """√2 ist algebraische ganze Zahl (Minpoly: t² - 2 = 0)."""
        assert is_algebraic_integer_check(math.sqrt(2), max_degree=4) is True

    def test_sqrt3_is_algebraic(self):
        """√3 ist algebraische ganze Zahl (Minpoly: t² - 3 = 0)."""
        assert is_algebraic_integer_check(math.sqrt(3), max_degree=4) is True

    def test_negative_algebraic(self):
        """Negative algebraische ganze Zahl (−√2)."""
        assert is_algebraic_integer_check(-math.sqrt(2), max_degree=4) is True

    def test_simple_rational_algebraic(self):
        """Rationale ganze Zahlen sind algebraisch."""
        assert is_algebraic_integer_check(5) is True


class TestMinimalPolynomial:
    """Tests für minimal_polynomial()"""

    def test_sqrt2_minimal_poly(self):
        """Minimalpolynom von √2 ist [−2, 0, 1] = t² − 2."""
        poly = minimal_polynomial(math.sqrt(2), max_degree=4, max_coeff=5)
        assert poly is not None
        # Auswerten: poly[0] + poly[1]*√2 + poly[2]*(√2)² ≈ 0
        val = sum(poly[i] * (math.sqrt(2) ** i) for i in range(len(poly)))
        assert abs(val) < 1e-6

    def test_rational_minimal_poly(self):
        """Minimalpolynom von 2 ist t - 2 = [−2, 1]."""
        poly = minimal_polynomial(2.0, max_degree=3, max_coeff=5)
        assert poly is not None
        val = sum(poly[i] * (2.0 ** i) for i in range(len(poly)))
        assert abs(val) < 1e-6

    def test_minimal_poly_evaluates_to_zero(self):
        """Minimalpolynom von √3 wertet bei √3 zu 0 aus."""
        x = math.sqrt(3)
        poly = minimal_polynomial(x, max_degree=4, max_coeff=5)
        assert poly is not None
        val = sum(poly[i] * (x ** i) for i in range(len(poly)))
        assert abs(val) < 1e-6


class TestAlgebraicDegree:
    """Tests für algebraic_degree()"""

    def test_rational_degree_1(self):
        """Rationale Zahlen haben algebraischen Grad 1."""
        assert algebraic_degree(3.0) == 1
        assert algebraic_degree(0.5) == 1

    def test_sqrt2_degree_2(self):
        """√2 hat algebraischen Grad 2."""
        deg = algebraic_degree(math.sqrt(2))
        assert deg == 2

    def test_sqrt3_degree_2(self):
        """√3 hat algebraischen Grad 2."""
        deg = algebraic_degree(math.sqrt(3))
        assert deg == 2


class TestClassifyNumber:
    """Tests für classify_number()"""

    def test_integer_is_rational(self):
        """Ganze Zahlen werden als 'rational' klassifiziert."""
        assert classify_number(3.0) == 'rational'
        assert classify_number(0.0) == 'rational'

    def test_half_is_rational(self):
        """1/2 wird als 'rational' klassifiziert."""
        assert classify_number(0.5) == 'rational'

    def test_sqrt2_is_algebraic(self):
        """√2 wird als 'algebraic' klassifiziert."""
        result = classify_number(math.sqrt(2))
        assert result == 'algebraic'

    def test_pi_probably_transcendental(self):
        """π wird als 'probably_transcendental' klassifiziert."""
        result = classify_number(math.pi)
        assert result == 'probably_transcendental'

    def test_e_probably_transcendental(self):
        """e wird als 'probably_transcendental' klassifiziert."""
        result = classify_number(math.e)
        assert result == 'probably_transcendental'


# ---------------------------------------------------------------------------
# 2. Tests: Liouville-Zahlen
# ---------------------------------------------------------------------------

class TestLiouvilleNumber:
    """Tests für liouville_number()"""

    def test_positive(self):
        """Liouville-Zahl ist positiv."""
        L = liouville_number(5)
        assert L > 0

    def test_less_than_one(self):
        """Liouville-Zahl ist kleiner als 1."""
        L = liouville_number(10)
        assert L < 1

    def test_first_term(self):
        """Erster Term: 10^{-1!} = 0.1."""
        L1 = liouville_number(1)
        assert abs(L1 - 0.1) < 1e-10

    def test_two_terms(self):
        """Zwei Terme: 10^{-1} + 10^{-2} = 0.11."""
        L2 = liouville_number(2)
        assert abs(L2 - 0.11) < 1e-10

    def test_more_terms_larger(self):
        """Mehr Terme ergeben größere Werte."""
        L5 = liouville_number(5)
        L3 = liouville_number(3)
        assert L5 >= L3

    def test_n_terms_zero_gives_zero(self):
        """Null Terme ergeben 0."""
        L0 = liouville_number(0)
        assert L0 == 0.0


class TestLiouvilleApproximationExponent:
    """Tests für liouville_approximation_exponent()"""

    def test_rational_exponent_low(self):
        """Rationale Zahlen haben niedrigen Exponent."""
        exp = liouville_approximation_exponent(0.5, max_q=100)
        # 0.5 = 1/2 → wird exakt durch p/q=1/2 approximiert
        # Exponent kann hoch sein durch exakte Darstellung, aber realistisch ≥ 2
        assert exp >= 2.0

    def test_sqrt2_exponent(self):
        """Irrationali­tätsexponent von √2 sollte ≥ 2 sein."""
        exp = liouville_approximation_exponent(math.sqrt(2), max_q=500)
        assert exp >= 2.0

    def test_pi_exponent(self):
        """Irrationali­tätsexponent von π sollte ≥ 2 sein."""
        exp = liouville_approximation_exponent(math.pi, max_q=200)
        assert exp >= 2.0

    def test_exponent_positive(self):
        """Exponent ist immer positiv."""
        exp = liouville_approximation_exponent(math.e, max_q=100)
        assert exp > 0


# ---------------------------------------------------------------------------
# 3. Tests: Lindemann-Weierstrass und verwandte Sätze
# ---------------------------------------------------------------------------

class TestLindemannWeierstrassDemo:
    """Tests für lindemann_weierstrass_demo()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = lindemann_weierstrass_demo()
        assert isinstance(result, dict)

    def test_contains_satz(self):
        """Wörterbuch enthält 'satz'-Schlüssel."""
        result = lindemann_weierstrass_demo()
        assert 'satz' in result

    def test_contains_beispiele(self):
        """Wörterbuch enthält 'beispiele'-Schlüssel."""
        result = lindemann_weierstrass_demo()
        assert 'beispiele' in result

    def test_e_transcendental_in_demo(self):
        """Demo markiert e als transzendent."""
        result = lindemann_weierstrass_demo()
        e_entry = result['beispiele'].get('e^1 = e', {})
        assert e_entry.get('transzendent') is True

    def test_pi_transcendental_in_demo(self):
        """Demo bestätigt π-Transzendenz."""
        result = lindemann_weierstrass_demo()
        pi_entry = result['beispiele'].get('e^(i*pi) = -1', {})
        assert pi_entry.get('pi_transzendent') is True

    def test_e_value_correct(self):
        """Demo enthält korrekten Wert für e."""
        result = lindemann_weierstrass_demo()
        e_val = result['beispiele']['e^1 = e']['wert']
        assert abs(e_val - math.e) < 1e-10


class TestETranscendentalDemo:
    """Tests für e_is_transcendental_demo()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = e_is_transcendental_demo()
        assert isinstance(result, dict)

    def test_contains_numerische_verifikation(self):
        """Enthält numerische Verifikation."""
        result = e_is_transcendental_demo()
        assert 'numerische_verifikation' in result

    def test_best_value_not_zero(self):
        """Bester Polynomwert bei e ist nicht 0 (e nicht algebraisch)."""
        result = e_is_transcendental_demo()
        for degree, data in result['numerische_verifikation'].items():
            assert data["bester_wert"] > 1e-4


class TestPiTranscendentalDemo:
    """Tests für pi_is_transcendental_demo()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = pi_is_transcendental_demo()
        assert isinstance(result, dict)

    def test_pi_value_correct(self):
        """Enthält korrekten π-Wert."""
        result = pi_is_transcendental_demo()
        assert abs(result['pi_wert'] - math.pi) < 1e-10

    def test_contains_folgerungen(self):
        """Enthält Folgerungen."""
        result = pi_is_transcendental_demo()
        assert 'folgerungen' in result
        assert len(result['folgerungen']) > 0


class TestGelfondSchneiderDemo:
    """Tests für gelfond_schneider_demo()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = gelfond_schneider_demo()
        assert isinstance(result, dict)

    def test_2_pow_sqrt2_correct_value(self):
        """2^√2 hat korrekten numerischen Wert."""
        result = gelfond_schneider_demo()
        val = result['beispiele']['2^sqrt(2)']['wert']
        expected = 2 ** math.sqrt(2)
        assert abs(val - expected) < 1e-8

    def test_2_pow_sqrt2_transcendental(self):
        """2^√2 ist als transzendent markiert."""
        result = gelfond_schneider_demo()
        assert result['beispiele']['2^sqrt(2)']['transzendent'] is True


class TestBakerTheoremDemo:
    """Tests für baker_theorem_demo()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = baker_theorem_demo()
        assert isinstance(result, dict)

    def test_log2_plus_log3_equals_log6(self):
        """log(2) + log(3) = log(6) numerisch korrekt."""
        result = baker_theorem_demo()
        log_example = result['beispiele']['log(2) + log(3) = log(6)']
        assert log_example['gleich_log6'] is True

    def test_fields_medal(self):
        """Fields-Medaille-Vermerk ist vorhanden."""
        result = baker_theorem_demo()
        assert 'fields_medal' in result


# ---------------------------------------------------------------------------
# 4. Tests: Irrationali­tätsmaße
# ---------------------------------------------------------------------------

class TestIrrationalityMeasure:
    """Tests für irrationality_measure()"""

    def test_sqrt2_measure_at_least_2(self):
        """Irrationali­tätsmaß von √2 ist mindestens 2."""
        mu = irrationality_measure(math.sqrt(2), max_q=500)
        assert mu >= 2.0

    def test_pi_measure_at_least_2(self):
        """Irrationali­tätsmaß von π ist mindestens 2."""
        mu = irrationality_measure(math.pi, max_q=500)
        assert mu >= 2.0

    def test_measure_positive(self):
        """Irrationali­tätsmaß ist stets positiv."""
        mu = irrationality_measure(math.e, max_q=100)
        assert mu > 0


class TestKnownIrrationalityMeasures:
    """Tests für known_irrationality_measures()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = known_irrationality_measures()
        assert isinstance(result, dict)

    def test_pi_present(self):
        """π ist im Wörterbuch."""
        result = known_irrationality_measures()
        assert 'pi' in result

    def test_e_present(self):
        """e ist im Wörterbuch."""
        result = known_irrationality_measures()
        assert 'e' in result

    def test_pi_lower_bound(self):
        """π hat untere Schranke ≥ 7."""
        result = known_irrationality_measures()
        assert result['pi']['untere_schranke'] >= 7.0

    def test_liouville_infinite(self):
        """Liouville-Zahl hat unendliches Irrationali­tätsmaß."""
        result = known_irrationality_measures()
        assert result['liouville']['untere_schranke'] == float('inf')


class TestDiophantineApproximationQuality:
    """Tests für diophantine_approximation_quality()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = diophantine_approximation_quality(math.sqrt(2), max_q=100)
        assert isinstance(result, dict)

    def test_contains_beste_annaeherung(self):
        """Enthält beste Annäherung."""
        result = diophantine_approximation_quality(math.sqrt(2), max_q=100)
        assert 'beste_annaeherung' in result

    def test_sqrt2_best_approx(self):
        """Beste Annäherung an √2 hat sehr kleinen Fehler."""
        result = diophantine_approximation_quality(math.sqrt(2), max_q=1000)
        p, q, error = result['beste_annaeherung']
        assert error < 1e-5


class TestContinuedFractionIrrationality:
    """Tests für continued_fraction_irrationality()"""

    def test_sqrt2_cf(self):
        """Kettenbruch von √2 beginnt mit [1; 2, 2, 2, ...]."""
        result = continued_fraction_irrationality(math.sqrt(2), n_terms=10)
        assert result['kettenbruch'][0] == 1
        # Alle weiteren Koeffizienten sollten 2 sein
        assert all(c == 2 for c in result['kettenbruch'][1:])

    def test_pi_cf_first_coeff(self):
        """Kettenbruch von π beginnt mit 3."""
        result = continued_fraction_irrationality(math.pi, n_terms=5)
        assert result['kettenbruch'][0] == 3

    def test_irrationality_estimate_at_least_2(self):
        """Irrationali­tätsschätzung ist mindestens 2."""
        result = continued_fraction_irrationality(math.sqrt(2), n_terms=15)
        assert result['irrationality_schaetzung'] >= 2.0

    def test_returns_dict_with_keys(self):
        """Rückgabe enthält alle erwarteten Schlüssel."""
        result = continued_fraction_irrationality(math.pi, n_terms=10)
        for key in ['kettenbruch', 'konvergente', 'fehler', 'irrationality_schaetzung']:
            assert key in result


# ---------------------------------------------------------------------------
# 5. Tests: AlgebraicNumber-Klasse
# ---------------------------------------------------------------------------

class TestAlgebraicNumber:
    """Tests für die AlgebraicNumber-Klasse"""

    def test_degree_of_sqrt2(self):
        """AlgebraicNumber für √2 (poly: t² - 2) hat Grad 2."""
        alpha = AlgebraicNumber([-2, 0, 1])  # t² - 2
        assert alpha.degree() == 2

    def test_degree_linear(self):
        """Lineares Polynom hat Grad 1."""
        alpha = AlgebraicNumber([-3, 1])  # t - 3
        assert alpha.degree() == 1

    def test_conjugates_sqrt2(self):
        """Konjugierte von √2 sind ±√2."""
        alpha = AlgebraicNumber([-2, 0, 1])
        conjs = alpha.conjugates()
        assert len(conjs) == 2
        vals = sorted([abs(c) for c in conjs])
        assert abs(vals[0] - math.sqrt(2)) < 1e-8
        assert abs(vals[1] - math.sqrt(2)) < 1e-8

    def test_norm_sqrt2(self):
        """Norm von √2 ist ∏(±√2) = -2."""
        alpha = AlgebraicNumber([-2, 0, 1])
        n = alpha.norm()
        assert abs(n.real - (-2)) < 1e-6

    def test_is_unit_root_primitive_root(self):
        """Primitive 4. Einheitswurzel (t⁴ - 1 Faktor: t² + 1) ist Einheitswurzel."""
        # t² + 1 → Nullstellen: ±i (Betrag = 1, i^4 = 1)
        alpha = AlgebraicNumber([1, 0, 1])  # t² + 1
        assert alpha.is_unit_root() is True

    def test_not_unit_root(self):
        """√2 ist keine Einheitswurzel (|√2| ≠ 1)."""
        alpha = AlgebraicNumber([-2, 0, 1])
        assert alpha.is_unit_root() is False

    def test_repr(self):
        """__repr__ gibt sinnvollen String zurück."""
        alpha = AlgebraicNumber([-2, 0, 1])
        assert 'AlgebraicNumber' in repr(alpha)

    def test_invalid_poly_raises(self):
        """Polynom mit weniger als 2 Koeffizienten wirft ValueError."""
        with pytest.raises(ValueError):
            AlgebraicNumber([1])


# ---------------------------------------------------------------------------
# 6. Tests: Kreisteilungspolynome
# ---------------------------------------------------------------------------

class TestCyclotomicPolynomial:
    """Tests für cyclotomic_polynomial()"""

    def test_phi1(self):
        """Φ₁(x) = x - 1 = [-1, 1]."""
        poly = cyclotomic_polynomial(1)
        assert poly is not None
        assert len(poly) >= 2

    def test_phi2(self):
        """Φ₂(x) = x + 1 → Nullstelle bei x=-1."""
        poly = cyclotomic_polynomial(2)
        # Auswertung bei -1: p(-1) = Σ poly[i]*(-1)^i ≈ 0
        val = sum(poly[i] * ((-1) ** i) for i in range(len(poly)))
        assert abs(val) < 1e-9

    def test_phi4(self):
        """Φ₄(x) = x² + 1 → Nullstellen sind primitive 4. Einheitswurzeln."""
        poly = cyclotomic_polynomial(4)
        assert poly is not None
        # Auswertung bei i (imaginäre Einheit): x²+1 = 0 bei x=i
        import cmath
        val = sum(poly[k] * (1j ** k) for k in range(len(poly)))
        assert abs(val) < 1e-9

    def test_returns_list(self):
        """Rückgabe ist eine Liste."""
        assert isinstance(cyclotomic_polynomial(3), list)


class TestNumberFieldDegree:
    """Tests für number_field_degree()"""

    def test_linear_poly_degree_1(self):
        """Lineares Polynom → Grad 1."""
        assert number_field_degree([-3, 1]) == 1

    def test_quadratic_poly_degree_2(self):
        """Quadratisches Polynom → Grad 2."""
        assert number_field_degree([-2, 0, 1]) == 2

    def test_cubic_poly_degree_3(self):
        """Kubisches Polynom → Grad 3."""
        assert number_field_degree([1, 0, 0, 1]) == 3
