"""
Tests für Mahler-Maß, Vollkommene Zahlen und Waring-Goldbach.

Testet alle drei Module:
    - mahler_measure.py: MahlerMeasure, SchurSiegelSmythTrace
    - perfect_numbers.py: SigmaFunction, EvenPerfectNumbers, OddPerfectNumberBounds
    - waring_goldbach.py: WaringGoldbachK2, WaringGoldbachGeneral

**Testabdeckung:**
    - Mahler-Maß (Lehmer-Polynom, Smyth-Polynom, Kreisteilungspolynome)
    - Gerade vollkommene Zahlen (erste 5: 6, 28, 496, 8128, 33550336)
    - Ungerade vollkommene Zahlen (Ausschluss, Touchard, Nielsen)
    - Sigma-Funktion (Multiplikativität, Abundanz)
    - Waring-Goldbach k=2 (10 = 1²+1²+2²+2², 30 als Quadrate)
    - Waring-Goldbach allgemein (Vinogradov, Goldbach-Vermutung)

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import pytest
import math
import sys
import os

# Sicherstelle, dass das src-Verzeichnis im Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from mahler_measure import MahlerMeasure, SchurSiegelSmythTrace
from perfect_numbers import SigmaFunction, EvenPerfectNumbers, OddPerfectNumberBounds
from waring_goldbach import WaringGoldbachK2, WaringGoldbachGeneral, _sieve_primes


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def approx_equal(a: float, b: float, tol: float = 1e-4) -> bool:
    """Prüft näherungsweise Gleichheit zweier Floats."""
    return abs(a - b) < tol


# ===========================================================================
# TESTS: MAHLER-MAß
# ===========================================================================

class TestMahlerMeasureBasic:
    """Grundlegende Tests für die MahlerMeasure-Klasse."""

    def test_lehmer_polynomial_mahler_measure_4_decimal(self):
        """Lehmer-Polynom M ≈ 1.17628 auf 4 Dezimalstellen."""
        mm = MahlerMeasure.lehmer_polynomial()
        m = mm.compute_product_formula()
        # Auf 4 Dezimalstellen: 1.1763 ± 0.0001
        assert abs(m - 1.17628) < 1e-4, f"Erwartet 1.17628, erhalten {m}"

    def test_lehmer_polynomial_degree(self):
        """Lehmer-Polynom hat Grad 10."""
        mm = MahlerMeasure.lehmer_polynomial()
        assert mm.degree == 10

    def test_lehmer_polynomial_coefficients(self):
        """Lehmer-Polynom hat korrekte Koeffizienten."""
        mm = MahlerMeasure.lehmer_polynomial()
        expected = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
        assert mm.coeffs == expected

    def test_lehmer_is_reciprocal(self):
        """Lehmer-Polynom ist reziprok (palindromisch)."""
        mm = MahlerMeasure.lehmer_polynomial()
        assert mm.is_reciprocal() is True

    def test_lehmer_smyth_bound_does_not_apply(self):
        """Smyth-Schranke gilt nicht für reziprokes Lehmer-Polynom."""
        mm = MahlerMeasure.lehmer_polynomial()
        assert mm.smyth_bound_applies() is False

    def test_smyth_polynomial_mahler_measure(self):
        """Smyth-Polynom x³-x-1 hat M ≈ 1.3247."""
        mm = MahlerMeasure.smyth_polynomial()
        m = mm.compute_product_formula()
        assert abs(m - MahlerMeasure.SMYTH_CONSTANT) < 1e-3, \
            f"Erwartet {MahlerMeasure.SMYTH_CONSTANT:.4f}, erhalten {m:.4f}"

    def test_smyth_polynomial_is_not_reciprocal(self):
        """Smyth-Polynom x³-x-1 ist nicht reziprok."""
        mm = MahlerMeasure.smyth_polynomial()
        assert mm.is_reciprocal() is False

    def test_smyth_bound_applies_to_smyth_polynomial(self):
        """Smyth-Schranke gilt für nicht-reziprokes Smyth-Polynom."""
        mm = MahlerMeasure.smyth_polynomial()
        assert mm.smyth_bound_applies() is True

    def test_linear_polynomial_mahler(self):
        """M(x - 2) = 2 (einzige Wurzel ist 2, außerhalb des Einheitskreises)."""
        mm = MahlerMeasure([-2, 1])  # x - 2
        m = mm.compute_product_formula()
        assert abs(m - 2.0) < 1e-6, f"Erwartet 2.0, erhalten {m}"

    def test_monic_polynomial_inside_unit_circle(self):
        """M(x - 1/2) = 1 (Wurzel bei 1/2, innerhalb des Einheitskreises)."""
        # Für ganzzahlige Polynome: M(2x - 1) = 2 (Leitkoeffizient 2)
        # Wurzel bei 1/2 (innen), Leitkoeffizient 2
        mm = MahlerMeasure([-1, 2])  # 2x - 1, Wurzel bei 1/2
        m = mm.compute_product_formula()
        assert abs(m - 2.0) < 1e-6

    def test_constant_polynomial(self):
        """M(c) = |c| für Konstanten (Grad-0-Polynom)."""
        mm = MahlerMeasure([3])
        m = mm.compute_product_formula()
        assert abs(m - 3.0) < 1e-6

    def test_product_formula_vs_jensen(self):
        """Produktformel und Jensen-Formel stimmen überein."""
        mm = MahlerMeasure.lehmer_polynomial()
        m_prod = mm.compute_product_formula()
        m_jensen = mm.compute_jensen_integral(num_points=50000)
        assert abs(m_prod - m_jensen) < 1e-3, \
            f"Produktformel {m_prod:.6f} ≠ Jensen {m_jensen:.6f}"

    def test_logarithmic_mahler_measure(self):
        """m(P) = log M(P)."""
        mm = MahlerMeasure.lehmer_polynomial()
        m = mm.compute_product_formula()
        log_m = mm.logarithmic_mahler_measure()
        assert abs(log_m - math.log(m)) < 1e-6


class TestMahlerMeasureCyclotomic:
    """Tests für Kreisteilungspolynome (M = 1 nach Kronecker)."""

    def test_cyclotomic_1_mahler_is_1(self):
        """Φ₁(x) = x-1 hat M = 1."""
        mm = MahlerMeasure.cyclotomic(1)
        m = mm.compute_product_formula()
        assert abs(m - 1.0) < 1e-6, f"Φ₁: M={m}"

    def test_cyclotomic_2_mahler_is_1(self):
        """Φ₂(x) = x+1 hat M = 1."""
        mm = MahlerMeasure.cyclotomic(2)
        m = mm.compute_product_formula()
        assert abs(m - 1.0) < 1e-6

    def test_cyclotomic_3_mahler_is_1(self):
        """Φ₃(x) = x²+x+1 hat M = 1."""
        mm = MahlerMeasure.cyclotomic(3)
        m = mm.compute_product_formula()
        assert abs(m - 1.0) < 1e-6

    def test_cyclotomic_6_mahler_is_1(self):
        """Φ₆(x) = x²-x+1 hat M = 1."""
        mm = MahlerMeasure.cyclotomic(6)
        m = mm.compute_product_formula()
        assert abs(m - 1.0) < 1e-6

    def test_cyclotomic_is_kronecker(self):
        """Kreisteilungspolynome erfüllen Kronecker-Theorem (M=1)."""
        for n in [1, 2, 3, 4, 5, 6, 7, 8, 12]:
            mm = MahlerMeasure.cyclotomic(n)
            assert mm.is_kronecker(tol=1e-5), f"Φ_{n} ist nicht Kronecker"

    def test_lehmer_not_kronecker(self):
        """Lehmer-Polynom ist kein Kronecker-Polynom (M > 1)."""
        mm = MahlerMeasure.lehmer_polynomial()
        assert mm.is_kronecker() is False

    def test_kronecker_theorem_x_polynomial(self):
        """M(x) = 1 (Nullpolynom-Ausnahme: Monome sind Kronecker)."""
        mm = MahlerMeasure([0, 1])  # x
        m = mm.compute_product_formula()
        assert abs(m - 1.0) < 1e-6


class TestMahlerMeasureMultiplication:
    """Tests für die Multiplikativität des Mahler-Maßes."""

    def test_mahler_multiplicativity(self):
        """M(P·Q) ≈ M(P)·M(Q) (Multiplikativitätseigenschaft)."""
        # P = x-2, Q = x-3: P·Q = x²-5x+6
        p = MahlerMeasure([-2, 1])   # M = 2
        q = MahlerMeasure([-3, 1])   # M = 3
        # Produkt: x²-5x+6 = [6, -5, 1]
        pq = MahlerMeasure([6, -5, 1])
        assert abs(pq.compute_product_formula() - 6.0) < 1e-5

    def test_mahler_scaling(self):
        """M(c·P) = |c|·M(P) für Skalar c."""
        mm = MahlerMeasure([-2, 1])  # x-2, M=2
        mm_scaled = MahlerMeasure([-6, 3])  # 3(x-2), M=6
        m1 = mm.compute_product_formula()
        m2 = mm_scaled.compute_product_formula()
        assert abs(m2 - 3 * m1) < 1e-5

    def test_invalid_zero_polynomial(self):
        """Nullpolynom löst ValueError aus."""
        with pytest.raises(ValueError):
            MahlerMeasure([0, 0, 0])

    def test_empty_polynomial(self):
        """Leeres Polynom löst ValueError aus."""
        with pytest.raises(ValueError):
            MahlerMeasure([])

    def test_leading_zero_removal(self):
        """Führende Nullen werden korrekt entfernt."""
        mm = MahlerMeasure([1, 2, 0])  # 0·x² + 2x + 1 = 2x + 1
        assert mm.degree == 1
        assert mm.coeffs == [1, 2]


class TestSchurSiegelSmyth:
    """Tests für den Schur-Siegel-Smyth-Spurpfad."""

    def test_trace_ratio_positive(self):
        """Spurverhältnis sollte positiv für total positive Polynome sein."""
        sss = SchurSiegelSmythTrace()
        # x²-3x+1: Wurzeln (3±√5)/2, beide positiv
        ratio = sss.trace_ratio([1, -3, 1])
        assert ratio is not None
        assert ratio > 0

    def test_trace_ratio_lehmer(self):
        """Spurverhältnis für Lehmer-Polynom."""
        sss = SchurSiegelSmythTrace()
        ratio = sss.trace_ratio(MahlerMeasure.LEHMER_COEFFS)
        assert ratio is not None

    def test_sss_constant(self):
        """SSS-Konstante liegt im richtigen Bereich."""
        assert 1.7 < SchurSiegelSmythTrace.SSS_CONSTANT < 1.9


# ===========================================================================
# TESTS: SIGMA-FUNKTION
# ===========================================================================

class TestSigmaFunction:
    """Tests für die σ-Funktion."""

    def test_sigma_prime(self):
        """σ(p) = p+1 für Primzahlen."""
        assert SigmaFunction.sigma(2) == 3
        assert SigmaFunction.sigma(3) == 4
        assert SigmaFunction.sigma(7) == 8
        assert SigmaFunction.sigma(13) == 14

    def test_sigma_prime_power(self):
        """σ(pᵏ) = (p^{k+1}-1)/(p-1)."""
        # σ(8) = σ(2³) = (2⁴-1)/(2-1) = 15
        assert SigmaFunction.sigma(8) == 15
        # σ(9) = σ(3²) = (3³-1)/(3-1) = 13
        assert SigmaFunction.sigma(9) == 13
        # σ(25) = σ(5²) = (5³-1)/(5-1) = 31
        assert SigmaFunction.sigma(25) == 31

    def test_sigma_1_is_1(self):
        """σ(1) = 1."""
        assert SigmaFunction.sigma(1) == 1

    def test_sigma_multiplicativity_check(self):
        """σ(mn) = σ(m)·σ(n) für teilerfremde m, n."""
        # 3 und 4 sind teilerfremd
        assert SigmaFunction.multiplicativity_check(3, 4) is True
        # 5 und 7 sind teilerfremd
        assert SigmaFunction.multiplicativity_check(5, 7) is True

    def test_sigma_multiplicativity_fails_for_non_coprime(self):
        """ValueError wenn m und n nicht teilerfremd."""
        with pytest.raises(ValueError):
            SigmaFunction.multiplicativity_check(4, 6)  # ggt(4,6) = 2

    def test_sigma_from_factorization(self):
        """σ aus Primfaktorzerlegung stimmt mit direktem σ überein."""
        # n = 12 = 2²·3: σ(12) = σ(4)·σ(3) = 7·4 = 28
        assert SigmaFunction.sigma_from_factorization({2: 2, 3: 1}) == 28
        assert SigmaFunction.sigma(12) == 28

    def test_sigma_invalid_input(self):
        """ValueError für n ≤ 0."""
        with pytest.raises(ValueError):
            SigmaFunction.sigma(0)
        with pytest.raises(ValueError):
            SigmaFunction.sigma(-5)

    def test_is_perfect_known_cases(self):
        """Bekannte vollkommene Zahlen werden erkannt."""
        for n in [6, 28, 496, 8128]:
            assert SigmaFunction.is_perfect(n), f"{n} sollte vollkommen sein"

    def test_is_not_perfect(self):
        """Nicht-vollkommene Zahlen werden korrekt abgelehnt."""
        for n in [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 100]:
            assert not SigmaFunction.is_perfect(n), f"{n} sollte nicht vollkommen sein"

    def test_abundant_numbers(self):
        """12, 18, 20 sind abundant."""
        for n in [12, 18, 20, 24]:
            assert SigmaFunction.is_abundant(n), f"{n} sollte abundant sein"

    def test_deficient_numbers(self):
        """Primzahlen und Primzahlpotenzen sind defizient."""
        for n in [2, 3, 4, 5, 7, 8, 9, 11, 13]:
            assert SigmaFunction.is_deficient(n), f"{n} sollte defizient sein"

    def test_abundancy_index_perfect(self):
        """Abundanzindex vollkommener Zahlen ist 2."""
        for n in [6, 28, 496]:
            idx = SigmaFunction.abundancy_index(n)
            assert abs(idx - 2.0) < 1e-10, f"σ({n})/{n} = {idx} ≠ 2"

    def test_perfect_numbers_up_to(self):
        """Findet alle vollkommenen Zahlen bis 10000."""
        perfects = SigmaFunction.perfect_numbers_up_to(10000)
        assert perfects == [6, 28, 496, 8128]

    def test_sigma_sequence_length(self):
        """sigma_sequence gibt korrekte Anzahl von Paaren zurück."""
        seq = SigmaFunction.sigma_sequence(10)
        assert len(seq) == 10
        assert seq[0] == (1, 1)
        assert seq[1] == (2, 3)


# ===========================================================================
# TESTS: GERADE VOLLKOMMENE ZAHLEN
# ===========================================================================

class TestEvenPerfectNumbers:
    """Tests für gerade vollkommene Zahlen (Euklid-Euler)."""

    def test_first_five_even_perfect_numbers(self):
        """Erste 5 gerade vollkommene Zahlen: 6, 28, 496, 8128, 33550336."""
        epn = EvenPerfectNumbers()
        first_five = epn.first_n_even_perfect(5)
        expected = [6, 28, 496, 8128, 33550336]
        assert first_five == expected, f"Erwartet {expected}, erhalten {first_five}"

    def test_first_even_perfect_is_6(self):
        """Erste gerade vollkommene Zahl ist 6 = 2¹(2²-1) = 2·3."""
        epn = EvenPerfectNumbers()
        first = epn.first_n_even_perfect(1)
        assert first == [6]

    def test_mersenne_number(self):
        """2^p - 1 Mersenne-Zahlen korrekt berechnet."""
        assert EvenPerfectNumbers.mersenne_number(2) == 3
        assert EvenPerfectNumbers.mersenne_number(3) == 7
        assert EvenPerfectNumbers.mersenne_number(5) == 31
        assert EvenPerfectNumbers.mersenne_number(7) == 127

    def test_is_mersenne_prime_known(self):
        """Bekannte Mersenne-Primzahlen werden erkannt."""
        for p in [2, 3, 5, 7, 13, 17, 19, 31]:
            assert EvenPerfectNumbers.is_mersenne_prime(p), \
                f"2^{p}-1 sollte Mersenne-Primzahl sein"

    def test_not_mersenne_prime(self):
        """Composite Mersenne-Zahlen werden abgelehnt."""
        # 2^11 - 1 = 2047 = 23 × 89 (nicht prim)
        assert not EvenPerfectNumbers.is_mersenne_prime(11)
        # 2^4 - 1 = 15 (nicht prim, aber 4 ist auch nicht prim)
        assert not EvenPerfectNumbers.is_mersenne_prime(4)

    def test_even_perfect_from_exponent_p2(self):
        """2^{2-1}(2^2-1) = 2·3 = 6."""
        n = EvenPerfectNumbers.even_perfect_from_exponent(2)
        assert n == 6

    def test_even_perfect_from_exponent_p3(self):
        """2^{3-1}(2^3-1) = 4·7 = 28."""
        n = EvenPerfectNumbers.even_perfect_from_exponent(3)
        assert n == 28

    def test_even_perfect_from_non_mersenne_prime(self):
        """even_perfect_from_exponent gibt None zurück wenn 2^p-1 nicht prim."""
        n = EvenPerfectNumbers.even_perfect_from_exponent(11)
        assert n is None

    def test_verify_euclid_euler_p2(self):
        """Verifikation des Euklid-Euler-Theorems für p=2."""
        result = EvenPerfectNumbers.verify_euclid_euler(2)
        assert result['p_is_prime'] is True
        assert result['mersenne_is_prime'] is True
        assert result['perfect_number'] == 6
        assert result['sigma_equals_2n'] is True

    def test_verify_euclid_euler_p3(self):
        """Verifikation des Euklid-Euler-Theorems für p=3."""
        result = EvenPerfectNumbers.verify_euclid_euler(3)
        assert result['perfect_number'] == 28
        assert result['sigma_equals_2n'] is True

    def test_euclid_euler_formula_str(self):
        """Formel-String enthält korrekte Werte."""
        s = EvenPerfectNumbers.euclid_euler_formula_str(2)
        assert '6' in s

    def test_digit_count_6(self):
        """Dezimalstellen-Zählung für n=6 korrekt."""
        count = EvenPerfectNumbers.even_perfect_digit_count(2)
        assert count == 1  # 6 hat 1 Stelle

    def test_digit_count_28(self):
        """Dezimalstellen-Zählung für n=28 korrekt."""
        count = EvenPerfectNumbers.even_perfect_digit_count(3)
        assert count == 2  # 28 hat 2 Stellen

    def test_first_8_even_perfect_numbers(self):
        """Erste 8 gerade vollkommene Zahlen entsprechen bekannten Werten."""
        epn = EvenPerfectNumbers()
        first_eight = epn.first_n_even_perfect(8)
        assert len(first_eight) == 8
        # Ersten 5 sind bekannt
        assert first_eight[:5] == [6, 28, 496, 8128, 33550336]

    def test_all_even_perfect_are_sigma_2n(self):
        """Alle ersten 4 geraden vollkommenen Zahlen erfüllen σ(n) = 2n."""
        for n in [6, 28, 496, 8128]:
            assert SigmaFunction.is_perfect(n), f"{n} sollte vollkommen sein"


# ===========================================================================
# TESTS: UNGERADE VOLLKOMMENE ZAHLEN
# ===========================================================================

class TestOddPerfectNumberBounds:
    """Tests für Schranken und Struktur ungerader vollkommener Zahlen."""

    def test_touchard_6_violates_congruence(self):
        """6 ≡ 6 (mod 12) ≠ 1 und ≡ 6 (mod 36) ≠ 9: Touchard verletzt."""
        result = OddPerfectNumberBounds.touchard_congruence(6)
        assert result['n_mod_12'] == 6
        assert not result['satisfies_n_equiv_1_mod_12']
        assert not result['satisfies_n_equiv_9_mod_36']
        assert not result['satisfies_touchard']

    def test_touchard_1_mod_12(self):
        """n=1 erfüllt Touchard: 1 ≡ 1 (mod 12)."""
        result = OddPerfectNumberBounds.touchard_congruence(1)
        assert result['satisfies_touchard'] is True

    def test_touchard_9_mod_36(self):
        """n=9 erfüllt Touchard: 9 ≡ 9 (mod 36)."""
        result = OddPerfectNumberBounds.touchard_congruence(9)
        assert result['satisfies_touchard'] is True

    def test_touchard_excludes_even_numbers_report(self):
        """Touchard-Test für gerade Zahlen (geraden vollkommenen)."""
        result = OddPerfectNumberBounds.touchard_congruence(28)
        # 28 % 12 = 4, 28 % 36 = 28 → Touchard verletzt für 28
        assert result['n_mod_12'] == 4
        assert not result['satisfies_touchard']

    def test_euler_form_check_odd(self):
        """Euler-Form-Check für kleine ungerade Zahlen."""
        # n=5: 5 = 5^1, Exponent 1 = 4·0+1, 5 ≡ 1 (mod 4)
        result = OddPerfectNumberBounds.euler_form_check(5)
        assert result['odd'] is True

    def test_euler_form_check_even_rejected(self):
        """Euler-Form-Check lehnt gerade Zahlen ab."""
        result = OddPerfectNumberBounds.euler_form_check(6)
        assert result['odd'] is False

    def test_nielsen_bound_small_number(self):
        """Kleine Zahlen erfüllen Nielsen-Schranke nicht (< 9 Primfaktoren)."""
        result = OddPerfectNumberBounds.nielsen_bound_check(15)  # 15 = 3·5
        assert result['num_distinct'] == 2
        assert not result['satisfies_nielsen']

    def test_nielsen_bound_constant(self):
        """Nielsen-Schranke ist 9."""
        assert OddPerfectNumberBounds.MIN_DISTINCT_PRIME_FACTORS == 9

    def test_chein_structure_prime(self):
        """Chein-Struktur für Primzahl (1 Primfaktor, Exponent 1)."""
        result = OddPerfectNumberBounds.chein_structure(7)
        assert result['exponents'] == [1]
        assert result['odd_exponent_count'] == 1
        assert result['has_exactly_one_odd_exponent'] is True

    def test_chein_multiplicity_minimum(self):
        """Mindest-Primfaktoranzahl mit Vielfachheit ist 101."""
        assert OddPerfectNumberBounds.MIN_PRIME_FACTORS_WITH_MULTIPLICITY == 101

    def test_is_excluded_simple_number(self):
        """Einfache Zahlen werden korrekt als nicht-vollkommen ausgeschlossen."""
        excluded, reasons = OddPerfectNumberBounds.is_excluded_as_odd_perfect(15)
        assert excluded is True  # σ(15) = 24 ≠ 30

    def test_summary_contains_conjecture_word(self):
        """Zusammenfassung enthält Conjecture-Hinweis."""
        summary = OddPerfectNumberBounds.summary_necessary_conditions()
        assert 'CONJECTURE' in summary or 'Conjecture' in summary.lower() or 'conjecture' in summary.lower()

    def test_lower_bound_brent_cohen(self):
        """Brent-Cohen-te Riele Schranke ist korrekt gesetzt."""
        assert OddPerfectNumberBounds.LOWER_BOUND_BRENT_COHEN == 10**300

    def test_touchard_n_equiv_1_mod_12_examples(self):
        """Verschiedene n ≡ 1 (mod 12) erfüllen Touchard."""
        for n in [1, 13, 25, 37, 49, 61]:
            result = OddPerfectNumberBounds.touchard_congruence(n)
            assert result['satisfies_touchard'], f"{n} ≡ 1 (mod 12) sollte Touchard erfüllen"

    def test_touchard_n_equiv_9_mod_36_examples(self):
        """Verschiedene n ≡ 9 (mod 36) erfüllen Touchard."""
        for n in [9, 45, 81, 117]:
            result = OddPerfectNumberBounds.touchard_congruence(n)
            assert result['satisfies_touchard'], f"{n} ≡ 9 (mod 36) sollte Touchard erfüllen"


# ===========================================================================
# TESTS: WARING-GOLDBACH k=2
# ===========================================================================

class TestWaringGoldbachK2:
    """Tests für Waring-Goldbach mit Primzahlquadraten (k=2)."""

    @pytest.fixture
    def wg2(self):
        """Waring-Goldbach k=2 Instanz."""
        return WaringGoldbachK2(max_prime=400)

    def test_hua_g2_constant(self):
        """G(2) ≤ 5 nach Hua (1938)."""
        assert WaringGoldbachK2.HUA_G2 == 5

    def test_13_is_sum_of_2_prime_squares(self):
        """13 = 2² + 3² = 4 + 9 (2 Primzahlquadrate).

        Hinweis: 10 ist NICHT als Summe von Primzahlquadraten darstellbar,
        da nur 4 und 9 ≤ 10 Primzahlquadrate sind (2²=4, 3²=9) und
        4+4=8, 4+9=13 ≠ 10. 1 ist keine Primzahl.
        """
        wg2 = WaringGoldbachK2(max_prime=20)
        rep = wg2.represent_as_prime_squares(13, max_terms=2)
        assert rep is not None, "13 sollte als Summe von 2 Primzahlquadraten darstellbar sein"
        squares = [p**2 for p in rep]
        assert sum(squares) == 13, f"Summe {sum(squares)} ≠ 13"
        assert len(rep) == 2

    def test_30_is_sum_of_prime_squares(self):
        """30 ist als Summe von Primzahlquadraten darstellbar."""
        wg2 = WaringGoldbachK2(max_prime=50)
        rep = wg2.represent_as_prime_squares(30, max_terms=5)
        assert rep is not None, "30 sollte als Summe von ≤5 Primzahlquadraten darstellbar sein"
        squares = [p**2 for p in rep]
        assert sum(squares) == 30

    def test_4_equals_2squared(self):
        """4 = 2² (1 Term)."""
        wg2 = WaringGoldbachK2(max_prime=20)
        rep = wg2.represent_as_prime_squares(4, max_terms=1)
        assert rep == [2]

    def test_9_equals_3squared(self):
        """9 = 3² (1 Term)."""
        wg2 = WaringGoldbachK2(max_prime=20)
        rep = wg2.represent_as_prime_squares(9, max_terms=1)
        assert rep == [3]

    def test_representation_sum_correct(self, wg2):
        """Darstellungen sind korrekt (Summe stimmt)."""
        for n in [4, 8, 9, 12, 13, 17, 20, 25, 30, 49, 50]:
            rep = wg2.represent_as_prime_squares(n, max_terms=5)
            if rep is not None:
                squares = [p**2 for p in rep]
                assert sum(squares) == n, f"Für n={n}: {rep} → Summe {sum(squares)} ≠ {n}"

    def test_local_conditions(self, wg2):
        """Lokale Bedingungen werden korrekt berechnet."""
        lc = wg2.local_conditions(29)
        assert lc['n'] == 29
        assert 'n_mod_8' in lc

    def test_count_representations_positive(self, wg2):
        """Anzahl der Darstellungen ist nicht-negativ."""
        count = wg2.count_representations(25, num_terms=2)
        assert count >= 0

    def test_25_is_5squared(self, wg2):
        """25 = 5² (1 Term)."""
        rep = wg2.represent_as_prime_squares(25, max_terms=1)
        assert rep == [5]

    def test_minimal_representation(self, wg2):
        """Minimale Darstellung hat korrekte Termsanzahl."""
        k, rep = wg2.minimal_representation(4)
        assert k == 1
        assert rep == [2]

    def test_5_prime_squares_sufficient(self, wg2):
        """Für viele Zahlen genügen ≤5 Primzahlquadrate."""
        successes = 0
        for n in range(4, 200):
            rep = wg2.represent_as_prime_squares(n, max_terms=5)
            if rep is not None:
                successes += 1
        # Mindestens 80% sollten darstellbar sein
        assert successes > 150, f"Nur {successes}/196 darstellbar"

    def test_invalid_n_zero(self, wg2):
        """n=0 ergibt None (kein gültiges Ergebnis)."""
        rep = wg2.represent_as_prime_squares(0, max_terms=5)
        assert rep is None


# ===========================================================================
# TESTS: WARING-GOLDBACH ALLGEMEIN
# ===========================================================================

class TestWaringGoldbachGeneral:
    """Tests für allgemeines Waring-Goldbach-Problem."""

    @pytest.fixture
    def wg3(self):
        """Waring-Goldbach k=3 Instanz."""
        return WaringGoldbachGeneral(k=3, max_prime=100)

    @pytest.fixture
    def wg1(self):
        """Waring-Goldbach k=1 Instanz (Goldbach/Vinogradov)."""
        return WaringGoldbachGeneral(k=1, max_prime=500)

    def test_invalid_k(self):
        """k < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            WaringGoldbachGeneral(k=0)

    def test_g_bounds_hua_k2(self):
        """G(2) = 5 nach Hua."""
        wg2 = WaringGoldbachGeneral(k=2, max_prime=100)
        assert wg2.g_bound_hua() == 5

    def test_g_bounds_hua_k3(self, wg3):
        """G(3) ≤ 9 nach Hua."""
        bound = wg3.g_bound_hua()
        assert bound == 9

    def test_g_bounds_hua_k4(self):
        """G(4) ≤ 15 nach Hua."""
        wg4 = WaringGoldbachGeneral(k=4, max_prime=50)
        assert wg4.g_bound_hua() == 15

    def test_vinogradov_wooley_bound_k2(self):
        """Vinogradov-Wooley Schranke für k=2 ist 5."""
        wg2 = WaringGoldbachGeneral(k=2, max_prime=100)
        assert wg2.g_bound_vinogradov_wooley() == 5

    def test_vinogradov_wooley_bound_k1(self, wg1):
        """Vinogradov-Wooley Schranke für k=1 ist 3."""
        assert wg1.g_bound_vinogradov_wooley() == 3

    def test_represent_cubes_of_primes(self, wg3):
        """Darstellung als Summe von Primzahlkuben."""
        # 8 = 2³ (1 Term)
        rep = wg3.represent_as_k_prime_powers(8, max_terms=1)
        assert rep == [2], f"8 = 2³, erhalten: {rep}"

    def test_represent_27_as_prime_cube(self, wg3):
        """27 = 3³ (1 Term)."""
        rep = wg3.represent_as_k_prime_powers(27, max_terms=1)
        assert rep == [3]

    def test_represent_sum_correct(self, wg3):
        """Darstellungen als k-te Potenzen sind korrekt."""
        for n in [8, 16, 27, 35, 100]:
            rep = wg3.represent_as_k_prime_powers(n, max_terms=9)
            if rep is not None:
                powers_sum = sum(p**3 for p in rep)
                assert powers_sum == n, f"Für n={n}: Summe={powers_sum} ≠ {n}"

    def test_vinogradov_three_primes_odd(self, wg1):
        """Vinogradov: Ungerade Zahl = p₁+p₂+p₃."""
        result = wg1.vinogradov_three_primes_check(15)
        # 15 = 3+5+7 oder andere Kombinationen
        assert result['is_odd'] is True
        assert result['found'] is True

    def test_vinogradov_even_rejected(self, wg1):
        """Vinogradov-Prüfung für gerade Zahlen gibt Hinweis."""
        result = wg1.vinogradov_three_primes_check(10)
        assert result['is_odd'] is False

    def test_goldbach_conjecture_check_basic(self, wg1):
        """Goldbach-Vermutung für kleine gerade Zahlen."""
        for n in [4, 6, 8, 10, 12, 20, 30, 100]:
            result = wg1.goldbach_conjecture_check(n)
            assert result['found'], f"Goldbach für n={n} nicht erfüllt"
            if result['representation']:
                p1, p2 = result['representation']
                assert p1 + p2 == n

    def test_goldbach_status_is_conjecture(self, wg1):
        """Goldbach-Status enthält 'CONJECTURE'."""
        result = wg1.goldbach_conjecture_check(10)
        assert 'CONJECTURE' in result['goldbach_status']

    def test_goldbach_invalid_input(self, wg1):
        """Goldbach für ungerade Zahlen gibt Fehler zurück."""
        result = wg1.goldbach_conjecture_check(7)
        assert 'error' in result

    def test_hua_theorem_statement(self, wg3):
        """Hua-Theorem-String für k=3 enthält G(3) ≤ 9."""
        stmt = wg3.hua_theorem_statement()
        assert '9' in stmt

    def test_explicit_vinogradov_wooley_string(self, wg3):
        """Vinogradov-Wooley String enthält Formel."""
        s = wg3.explicit_vinogradov_wooley_bound()
        assert 'G(3)' in s or 'G(' in s


# ===========================================================================
# TESTS: SIEB-HILFSFUNKTION
# ===========================================================================

class TestSievePrimes:
    """Tests für den Sieb des Eratosthenes."""

    def test_primes_up_to_10(self):
        """Primzahlen bis 10: [2, 3, 5, 7]."""
        assert _sieve_primes(10) == [2, 3, 5, 7]

    def test_primes_up_to_2(self):
        """Primzahlen bis 2: [2]."""
        assert _sieve_primes(2) == [2]

    def test_primes_up_to_1(self):
        """Keine Primzahlen bis 1."""
        assert _sieve_primes(1) == []

    def test_primes_up_to_0(self):
        """Keine Primzahlen bis 0."""
        assert _sieve_primes(0) == []

    def test_primes_count_up_to_100(self):
        """Bis 100 gibt es 25 Primzahlen."""
        assert len(_sieve_primes(100)) == 25

    def test_primes_all_prime(self):
        """Alle gefundenen Zahlen sind prim."""
        import sympy
        primes = _sieve_primes(50)
        for p in primes:
            assert sympy.isprime(p), f"{p} ist keine Primzahl"


# ===========================================================================
# INTEGRATIONSTESTS
# ===========================================================================

class TestIntegration:
    """Integrationstests über mehrere Module."""

    def test_perfect_numbers_sigma_consistency(self):
        """σ(n) = 2n ↔ is_perfect(n): Konsistenz."""
        from sympy import divisor_sigma
        for n in range(1, 10000):
            expected = int(divisor_sigma(n)) == 2 * n
            actual = SigmaFunction.is_perfect(n)
            if expected != actual:
                pytest.fail(f"Inkonsistenz bei n={n}: σ-Test={expected}, is_perfect={actual}")

    def test_even_perfect_euclid_euler_consistency(self):
        """Euklid-Euler: Alle ersten 4 geraden vollkommenen Zahlen haben Mersenne-Form."""
        epn = EvenPerfectNumbers()
        for p in [2, 3, 5, 7]:
            n = epn.even_perfect_from_exponent(p)
            assert n is not None
            assert SigmaFunction.is_perfect(n), f"n={n} (p={p}) ist nicht vollkommen"

    def test_mahler_larger_than_1_implies_not_kronecker(self):
        """M(P) > 1 impliziert P ist kein Kronecker-Polynom."""
        # Lehmer: M > 1
        mm = MahlerMeasure.lehmer_polynomial()
        m = mm.compute_product_formula()
        assert m > 1.0
        assert not mm.is_kronecker()

    def test_smyth_constant_greater_than_lehmer(self):
        """Smyth-Konstante > Lehmer-Maß (wie erwartet für nicht-reziproke Polynome)."""
        assert MahlerMeasure.SMYTH_CONSTANT > MahlerMeasure.LEHMER_MAHLER

    def test_waring_goldbach_primes_are_prime(self):
        """Alle zurückgegebenen 'Primzahlen' in Darstellungen sind tatsächlich prim."""
        import sympy
        wg2 = WaringGoldbachK2(max_prime=100)
        for n in [4, 9, 25, 49, 50]:
            rep = wg2.represent_as_prime_squares(n, max_terms=5)
            if rep is not None:
                for p in rep:
                    assert sympy.isprime(p), f"{p} in Darstellung von {n} ist keine Primzahl"
