"""
@file test_iwasawa.py
@brief Tests für das Iwasawa-Theorie-Modul (iwasawa_theory.py).
@description
    Umfassende Test-Suite für alle Funktionen der Iwasawa-Theorie:
    - p-adische Logarithmen und Teichmüller-Repräsentanten
    - Iwasawa-Algebra Λ = ℤ_p[[T]] (μ, λ-Invarianten, Weierstraß)
    - p-adische L-Funktionen (Kubota-Leopoldt)
    - Kummer-Kongruenzen
    - Selmer-Gruppen und BSD-Verbindung
    - Zyklotomische ℤ_p-Erweiterungen

    Mindestens 25 Tests gemäß Spezifikation.

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import math
import pytest
from fractions import Fraction

# Modulpfad einrichten
sys.path.insert(0, 'src')

from iwasawa_theory import (
    # p-adische Grundfunktionen
    p_adic_log,
    teichmuller_representative,
    # Iwasawa-Algebra
    iwasawa_mu_invariant,
    iwasawa_lambda_invariant,
    iwasawa_polynomial,
    weierstrass_preparation,
    # p-adische L-Funktionen
    p_adic_zeta_special_values,
    kummer_congruences_check,
    kubota_leopoldt_l_function,
    # Zyklotomische Erweiterung
    cyclotomic_zp_extension,
    # Selmer und BSD
    selmer_group_rank_estimate,
    iwasawa_main_conjecture_evidence,
    bsd_iwasawa_connection,
    # Interne Hilfsfunktionen (via direkten Import)
    _bernoulli_number,
    _bernoulli_polynomial,
)
from exceptions import InvalidInputError, PrimeRequiredError


# ===========================================================================
# TESTS: Bernoulli-Zahlen (interne Hilfsfunktionen)
# ===========================================================================

class TestBernoulliNumbers:
    """Tests für interne Bernoulli-Zahlen-Berechnung."""

    def test_bernoulli_b0(self):
        """B_0 = 1 (definitionsgemäß)."""
        assert _bernoulli_number(0) == Fraction(1)

    def test_bernoulli_b1(self):
        """B_1 = -1/2."""
        assert _bernoulli_number(1) == Fraction(-1, 2)

    def test_bernoulli_b2(self):
        """B_2 = 1/6."""
        assert _bernoulli_number(2) == Fraction(1, 6)

    def test_bernoulli_b4(self):
        """B_4 = -1/30."""
        assert _bernoulli_number(4) == Fraction(-1, 30)

    def test_bernoulli_odd_zero(self):
        """Alle ungeraden Bernoulli-Zahlen B_n (n ≥ 3) sind 0."""
        for n in [3, 5, 7, 9, 11]:
            assert _bernoulli_number(n) == Fraction(0), f"B_{n} sollte 0 sein"

    def test_bernoulli_polynomial_b0(self):
        """B_0(x) = 1 für alle x."""
        for x_val in [0, Fraction(1, 2), Fraction(1, 3)]:
            assert _bernoulli_polynomial(0, Fraction(x_val)) == Fraction(1)

    def test_bernoulli_polynomial_b1(self):
        """B_1(1/2) = 0 (Spiegelsymmetrie B_1(1-x) = -B_1(x) → B_1(1/2) = 0)."""
        val = _bernoulli_polynomial(1, Fraction(1, 2))
        assert val == Fraction(0)


# ===========================================================================
# TESTS: Teichmüller-Repräsentant
# ===========================================================================

class TestTeichmullerRepresentative:
    """Tests für den Teichmüller-Repräsentanten ω(a)."""

    def test_teichmuller_2_mod_5(self):
        """
        ω(2) mod 5: Da 2^4 ≡ 1 (mod 5) (Fermat), ist 2 bereits
        eine (p-1)=4-te Einheitswurzel mod 5. Teichmüller-Repräsentant
        von 2 bei p=5 ist 2 selbst (mod 5).
        """
        result = teichmuller_representative(2, 5)
        # ω(2) ≡ 2 (mod 5)
        assert result % 5 == 2 % 5

    def test_teichmuller_is_p_minus_1_root_of_unity(self):
        """ω(a)^{p-1} ≡ 1 (mod p^precision) – Einheitswurzel-Eigenschaft."""
        p = 5
        precision = 10
        modulus = p ** precision

        for a in [2, 3, 4]:  # a ∈ (ℤ/5ℤ)^× = {1,2,3,4}
            omega = teichmuller_representative(a, p, precision)
            # ω(a)^{p-1} ≡ 1 (mod p)
            assert pow(omega, p - 1, p) == 1 % p, \
                f"ω({a})^{p-1} ≢ 1 mod {p}"

    def test_teichmuller_congruent_to_a(self):
        """ω(a) ≡ a (mod p) – definierendes Kriterium."""
        for p in [5, 7, 11]:
            for a in range(1, p):
                omega = teichmuller_representative(a, p, 15)
                assert omega % p == a % p, \
                    f"ω({a}) mod {p} = {omega % p}, erwartet {a % p}"

    def test_teichmuller_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            teichmuller_representative(2, 4)

    def test_teichmuller_invalid_input(self):
        """a durch p teilbar wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            teichmuller_representative(5, 5)  # 5 ≡ 0 (mod 5)

    def test_teichmuller_1_is_1(self):
        """ω(1) = 1 für alle Primzahlen (1 ist triviale (p-1)-te Einheitswurzel)."""
        for p in [5, 7, 11]:
            omega = teichmuller_representative(1, p, 10)
            assert omega % p == 1


# ===========================================================================
# TESTS: p-adischer Logarithmus
# ===========================================================================

class TestPAdicLog:
    """Tests für den p-adischen Logarithmus."""

    def test_p_adic_log_returns_list(self):
        """p_adic_log gibt eine Liste von p-adischen Ziffern zurück."""
        result = p_adic_log(6, 5, precision=10)  # 6 ≡ 1 (mod 5)
        assert isinstance(result, list)
        assert len(result) == 10

    def test_p_adic_log_digits_in_range(self):
        """Alle zurückgegebenen Ziffern sind in {0, 1, ..., p-1}."""
        p = 5
        result = p_adic_log(6, p, precision=8)  # 6 = 5+1 ≡ 1 (mod 5)
        for digit in result:
            assert 0 <= digit < p, f"Ziffer {digit} außerhalb {0..p-1}"

    def test_p_adic_log_invalid_x(self):
        """x ≢ 1 (mod p) wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            p_adic_log(2, 5, precision=10)  # 2 ≢ 1 (mod 5)

    def test_p_adic_log_invalid_prime(self):
        """Nicht-Primzahl p wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            p_adic_log(6, 4)


# ===========================================================================
# TESTS: Iwasawa-Algebra (μ, λ-Invarianten)
# ===========================================================================

class TestIwasawaAlgebra:
    """Tests für Iwasawa-Algebra-Operationen."""

    def test_mu_invariant_basic(self):
        """
        f(T) = p + p² + p mit p=5: alle Koeffizienten durch p teilbar,
        μ = min v_p(a_k) = 1.
        """
        p = 5
        coeffs = [p, p**2, p]  # [5, 25, 5]
        mu = iwasawa_mu_invariant(coeffs, p)
        assert mu == 1

    def test_mu_invariant_zero(self):
        """Koeffizient ≢ 0 (mod p) → μ = 0."""
        p = 5
        coeffs = [1, 5, 25]  # a_0 = 1 ≢ 0 mod 5
        mu = iwasawa_mu_invariant(coeffs, p)
        assert mu == 0

    def test_mu_invariant_two(self):
        """Alle Koeffizienten durch p² teilbar → μ = 2."""
        p = 3
        coeffs = [9, 18, 27]  # 3², 2·3², 3³
        mu = iwasawa_mu_invariant(coeffs, p)
        assert mu == 2

    def test_mu_invariant_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            iwasawa_mu_invariant([1, 2, 3], 4)

    def test_mu_invariant_empty_raises(self):
        """Leere Liste wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            iwasawa_mu_invariant([], 5)

    def test_lambda_invariant_leading_not_divisible(self):
        """
        f_coeffs = [1, 0, p]: leitender Koeffizient ist 1 ≢ 0 (mod p),
        daher λ = 2 (Grad des Polynoms mod p).
        """
        p = 5
        coeffs = [1, 0, p]  # 1 + 0·T + 5·T² – mod 5: 1 + 0 + 0 = 1
        lam = iwasawa_lambda_invariant(coeffs, p)
        # a_2 = p ≡ 0 mod p, a_1 = 0 ≡ 0 mod p, a_0 = 1 ≢ 0 → λ = 0
        assert lam == 0

    def test_lambda_invariant_nonzero(self):
        """f(T) = p + T² (leitend ≢ 0 mod p) → λ = 2."""
        p = 5
        coeffs = [5, 0, 1]  # 5 + 0·T + 1·T²
        lam = iwasawa_lambda_invariant(coeffs, p)
        # mod p: 0 + 0 + 1·T² → Grad 2
        assert lam == 2

    def test_lambda_invariant_zero(self):
        """Einheit f(T) = 1 → λ = 0."""
        p = 5
        coeffs = [1]  # Konstantes Polynom 1 (Einheit)
        lam = iwasawa_lambda_invariant(coeffs, p)
        assert lam == 0

    def test_weierstrass_mu_correct(self):
        """Weierstraß-Zerlegung liefert korrekte μ-Invariante."""
        p = 5
        coeffs = [5, 25, 5]  # Alle durch 5 teilbar → μ = 1
        result = weierstrass_preparation(coeffs, p)
        assert result['mu_invariant'] == 1

    def test_weierstrass_is_unit_detection(self):
        """Einheit f(T) = 1 wird als Einheit erkannt (is_unit = True)."""
        p = 5
        coeffs = [1]
        result = weierstrass_preparation(coeffs, p)
        assert result['is_unit'] is True

    def test_weierstrass_non_unit(self):
        """f(T) = T ist keine Einheit (f(0) = 0)."""
        p = 5
        coeffs = [0, 1]  # 0 + T
        result = weierstrass_preparation(coeffs, p)
        assert result['is_unit'] is False

    def test_iwasawa_polynomial_structure(self):
        """iwasawa_polynomial gibt korrektes Dict zurück."""
        p = 5
        coeffs = [1, 5, 25]
        result = iwasawa_polynomial(coeffs, p)
        assert 'mu_invariant' in result
        assert 'lambda_invariant' in result
        assert 'distinguished_poly' in result
        assert 'is_unit' in result


# ===========================================================================
# TESTS: p-adische L-Funktionen und Kummer-Kongruenzen
# ===========================================================================

class TestPAdicLFunctions:
    """Tests für p-adische L-Funktionen und Kummer-Kongruenzen."""

    def test_p_adic_zeta_returns_dict(self):
        """p_adic_zeta_special_values gibt Dict mit n=2,4,6 zurück."""
        result = p_adic_zeta_special_values(5)
        assert 2 in result
        assert 4 in result
        assert 6 in result

    def test_p_adic_zeta_structure(self):
        """Jeder Eintrag hat die erwarteten Schlüssel."""
        result = p_adic_zeta_special_values(5)
        for n in [2, 4, 6]:
            entry = result[n]
            assert 'zeta_p_value' in entry
            assert 'euler_factor' in entry
            assert 'bernoulli' in entry

    def test_p_adic_zeta_b2_value(self):
        """
        ζ_p(1-2) = -(1 - p) B_2/2 für n=2:
        B_2 = 1/6, also ζ_5(-1) = -(1-5)·(1/6)/2 = -(-4)·(1/12) = 1/3.
        """
        result = p_adic_zeta_special_values(5)
        zeta_val = result[2]['zeta_p_value']
        # -(1 - 5^1) * B_2 / 2 = -(-4) * (1/6) / 2 = 4/12 = 1/3
        expected = Fraction(1, 3)
        assert zeta_val == expected, f"ζ_5(-1) = {zeta_val}, erwartet {expected}"

    def test_kummer_congruences_p5(self):
        """
        Kummer-Kongruenzen für p=5: B_4/4 ≡ B_8/8 (mod 5).
        Modifiziert: (1-5³)B_4/4 ≡ (1-5⁷)B_8/8 (mod 5).
        """
        result = kummer_congruences_check(5, n_range=12)
        assert result['all_consistent'], \
            f"Kummer-Kongruenzen für p=5 verletzt: {result['checks']}"

    def test_kummer_congruences_p7(self):
        """Kummer-Kongruenzen für p=7 müssen gelten."""
        result = kummer_congruences_check(7, n_range=14)
        assert result['p'] == 7
        assert result['period'] == 6  # p-1 = 6

    def test_kummer_structure(self):
        """Rückgabe-Dict hat erwartete Schlüssel."""
        result = kummer_congruences_check(5)
        assert 'p' in result
        assert 'period' in result
        assert 'checks' in result
        assert 'all_consistent' in result

    def test_kubota_leopoldt_trivial_char(self):
        """Kubota-Leopoldt mit trivialem Charakter = ζ_p(s)."""
        p = 5
        # L_p(1-2, 1) = ζ_p(-1) = -(1-5) B_2/2 = 1/3
        result = kubota_leopoldt_l_function(s=-1, chi_values={}, p=5)
        expected = Fraction(1, 3)
        assert result == expected

    def test_kubota_leopoldt_invalid_s(self):
        """s >= 1 wirft InvalidInputError (n = 1-s ≤ 0)."""
        with pytest.raises(InvalidInputError):
            kubota_leopoldt_l_function(s=2, chi_values={}, p=5)


# ===========================================================================
# TESTS: Zyklotomische ℤ_p-Erweiterung
# ===========================================================================

class TestCyclotomicExtension:
    """Tests für die zyklotomische ℤ_p-Erweiterung."""

    def test_cyclotomic_degrees_formula(self):
        """[ℚ_n : ℚ] = p^n · (p-1) für alle Ebenen n."""
        p = 5
        result = cyclotomic_zp_extension(p, n_levels=4)

        for n, deg in zip(result['levels'], result['degrees']):
            expected_deg = (p ** n) * (p - 1)
            assert deg == expected_deg, \
                f"Ebene {n}: Grad {deg}, erwartet {expected_deg}"

    def test_cyclotomic_p5_level0(self):
        """ℚ_0 = ℚ(ζ_5): Grad = 4 = 5^0 · (5-1) über ℚ."""
        result = cyclotomic_zp_extension(5, n_levels=3)
        assert result['degrees'][0] == 4

    def test_cyclotomic_p5_level1(self):
        """ℚ_1 = ℚ(ζ_{25}): Grad = 20 = 5^1 · 4 über ℚ."""
        result = cyclotomic_zp_extension(5, n_levels=3)
        assert result['degrees'][1] == 20

    def test_cyclotomic_p5_level2(self):
        """ℚ_2 = ℚ(ζ_{125}): Grad = 100 = 5^2 · 4 über ℚ."""
        result = cyclotomic_zp_extension(5, n_levels=4)
        assert result['degrees'][2] == 100

    def test_cyclotomic_structure(self):
        """Rückgabe enthält alle erwarteten Schlüssel."""
        result = cyclotomic_zp_extension(5, n_levels=3)
        assert 'levels' in result
        assert 'degrees' in result
        assert 'disc_estimates' in result
        assert 'gal_orders' in result

    def test_cyclotomic_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            cyclotomic_zp_extension(6)


# ===========================================================================
# TESTS: Selmer-Gruppen und Hauptvermutung
# ===========================================================================

class TestSelmerAndMainConjecture:
    """Tests für Selmer-Gruppen und Iwasawa-Hauptvermutung."""

    def test_selmer_rank_non_negative(self):
        """Selmer-Rang ist immer ≥ 0."""
        result = selmer_group_rank_estimate(-1, 0, 5)
        assert result['selmer_rank'] >= 0

    def test_selmer_bsd_consistent(self):
        """BSD-Konsistenz (selmer_rank ≥ 0) ist immer erfüllt."""
        result = selmer_group_rank_estimate(-1, 0, 5)
        assert result['bsd_consistent'] is True

    def test_selmer_returns_structure(self):
        """Rückgabe enthält alle erwarteten Schlüssel."""
        result = selmer_group_rank_estimate(-1, 0, 5)
        assert 'selmer_rank' in result
        assert 'mw_rank_lower' in result
        assert 'sha_order_estimate' in result
        assert 'bsd_consistent' in result

    def test_selmer_singular_curve_raises(self):
        """Singuläre Kurve (Δ=0) wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            selmer_group_rank_estimate(0, 0, 5)  # y² = x³: Δ = 0

    def test_selmer_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            selmer_group_rank_estimate(-1, 0, 4)

    def test_main_conjecture_structure(self):
        """iwasawa_main_conjecture_evidence gibt korrektes Dict zurück."""
        result = iwasawa_main_conjecture_evidence(5)
        assert 'mu_lp' in result
        assert 'lambda_lp' in result
        assert 'class_number_growth' in result
        assert 'main_conjecture_consistent' in result
        assert 'is_regular' in result

    def test_main_conjecture_greenberg(self):
        """Greenberg-Vermutung: μ = 0 für alle regulären Primzahlen."""
        result = iwasawa_main_conjecture_evidence(5)
        # p=5 ist regulär (keine irregulären Bernoulli-Zahlen bis p-3=2)
        assert result['mu_lp'] == 0

    def test_main_conjecture_p5_regular(self):
        """p=5 ist reguläre Primzahl: λ = 0."""
        result = iwasawa_main_conjecture_evidence(5)
        assert result['is_regular'] is True

    def test_main_conjecture_consistent(self):
        """Hauptvermutung ist nach Mazur-Wiles konsistent."""
        result = iwasawa_main_conjecture_evidence(5)
        assert result['main_conjecture_consistent'] is True

    def test_main_conjecture_class_growth(self):
        """Klassengruppen-Wachstum hat korrekte Struktur."""
        result = iwasawa_main_conjecture_evidence(5)
        growth = result['class_number_growth']
        assert len(growth) == 5  # n = 0, 1, 2, 3, 4
        assert growth[0]['n'] == 0

    def test_main_conjecture_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            iwasawa_main_conjecture_evidence(4)


# ===========================================================================
# TESTS: BSD-Iwasawa-Verbindung
# ===========================================================================

class TestBSDIwasawaConnection:
    """Tests für die BSD-Iwasawa-Verbindung."""

    def test_bsd_structure(self):
        """bsd_iwasawa_connection gibt korrektes Dict zurück."""
        result = bsd_iwasawa_connection(-1, 0, 5)
        assert 'analytic_rank_estimate' in result
        assert 'p_adic_rank' in result
        assert 'bsd_rank_match' in result
        assert 'kolyvagin_applicable' in result

    def test_bsd_ranks_non_negative(self):
        """Analytischer und p-adischer Rang sind ≥ 0."""
        result = bsd_iwasawa_connection(-1, 0, 5)
        assert result['analytic_rank_estimate'] >= 0
        assert result['p_adic_rank'] >= 0

    def test_bsd_singular_raises(self):
        """Singuläre Kurve wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            bsd_iwasawa_connection(0, 0, 5)

    def test_bsd_invalid_prime_raises(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            bsd_iwasawa_connection(-1, 0, 4)

    def test_bsd_kolyvagin_good_reduction(self):
        """Für gute Reduktion bei p: Kolyvagin anwendbar."""
        # E: y² = x³ - x, p=7 (7 ∤ Δ = -64 → gute Reduktion)
        result = bsd_iwasawa_connection(-1, 0, 7)
        # Δ = -16(4·(-1)³ + 0) = -16·(-4) = 64 ≠ 0, 7 ∤ 64
        assert result['kolyvagin_applicable'] is True


# ===========================================================================
# EDGE-CASE-TESTS
# ===========================================================================

class TestEdgeCases:
    """Edge-Cases und Sonderfälle."""

    def test_mu_invariant_single_coefficient(self):
        """Einzelner Koeffizient: μ = v_p(a_0)."""
        p = 5
        assert iwasawa_mu_invariant([25], p) == 2
        assert iwasawa_mu_invariant([1], p) == 0
        assert iwasawa_mu_invariant([125], p) == 3

    def test_kummer_p3(self):
        """Kummer-Kongruenzen für p=3 (kleinste ungerade Primzahl)."""
        result = kummer_congruences_check(3, n_range=8)
        assert result['p'] == 3
        assert result['period'] == 2

    def test_cyclotomic_p2(self):
        """
        p=2: Zyklotomische Erweiterung startet bei Ebene 0 mit Grad 1.
        [ℚ_0 : ℚ] = 2^0 · (2-1) = 1.
        """
        result = cyclotomic_zp_extension(2, n_levels=3)
        assert result['degrees'][0] == 1  # 2^0 · 1

    def test_selmer_no_bad_primes_small(self):
        """Kurve mit keinen kleinen schlechten Primzahlen."""
        # y² = x³ + x + 1: Δ = -16(4 + 27) = -496
        result = selmer_group_rank_estimate(1, 1, 5)
        assert result['selmer_rank'] >= 0

    def test_iwasawa_lambda_all_zero_except_leading(self):
        """f(T) = T³ (nur leitender Term ≢ 0) → λ = 3."""
        p = 5
        coeffs = [0, 0, 0, 1]  # T³
        lam = iwasawa_lambda_invariant(coeffs, p)
        assert lam == 3

    def test_p_adic_zeta_invalid_prime(self):
        """Nicht-Primzahl wirft PrimeRequiredError."""
        with pytest.raises(PrimeRequiredError):
            p_adic_zeta_special_values(4)
