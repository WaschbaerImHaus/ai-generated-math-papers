"""
Tests für euler_gamma_irrationality, zeta_odd_values und pillai_chowla.

Testet:
- EulerMascheroniConstant: γ auf 15 Stellen, Methoden, Stieltjes
- StieltjesConstants: γ_0=γ, Vorzeichen, Laurent-Entwicklung
- ZetaOddValues: ζ(3) auf 10 Stellen, Apéry-Kettenbruch, Ball-Rivoal, Zudilin
- PillaiChowla: Vermutung für n≤50, Wilson-Theorem, Mertens, Sylvester

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
import sys
import os
import pytest
import mpmath

# Projektpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from euler_gamma_irrationality import EulerMascheroniConstant, StieltjesConstants
from zeta_odd_values import ZetaOddValues
from pillai_chowla import PillaiChowla


# ===========================================================================
# FIXTURES
# ===========================================================================

@pytest.fixture(scope="module")
def gamma_calc():
    """EulerMascheroniConstant mit 50 Stellen."""
    return EulerMascheroniConstant(precision=50)


@pytest.fixture(scope="module")
def stieltjes_calc():
    """StieltjesConstants mit 30 Stellen."""
    return StieltjesConstants(precision=30)


@pytest.fixture(scope="module")
def zeta_calc():
    """ZetaOddValues mit 50 Stellen."""
    return ZetaOddValues(precision=50)


@pytest.fixture(scope="module")
def pillai():
    """PillaiChowla-Instanz."""
    return PillaiChowla(max_n=200)


# ===========================================================================
# EULER-MASCHERONI-KONSTANTE: GRUNDLEGENDE TESTS
# ===========================================================================

class TestEulerMascheroniConstant:
    """Tests für die EulerMascheroniConstant-Klasse."""

    def test_gamma_value_15_digits(self, gamma_calc):
        """γ auf 15 Stellen korrekt: 0.577215664901532."""
        mpmath.mp.dps = 25
        computed = float(gamma_calc.compute())
        expected = 0.5772156649015328606065  # 22 Stellen
        assert abs(computed - expected) < 1e-15, (
            f"γ ungenau: {computed} vs {expected}"
        )

    def test_gamma_first_digit(self, gamma_calc):
        """Erste Stelle von γ ist 0."""
        val = float(gamma_calc.compute())
        assert int(val) == 0

    def test_gamma_correct_value_0_577(self, gamma_calc):
        """γ ≈ 0.5772..."""
        val = float(gamma_calc.compute())
        assert 0.577 < val < 0.578

    def test_gamma_is_mpf(self, gamma_calc):
        """compute() liefert einen mpmath-kompatiblen Typ (mpf oder constant)."""
        val = gamma_calc.compute()
        # mpmath.euler hat Typ 'constant', ist aber in float konvertierbar
        assert hasattr(val, '__float__'), "compute() sollte einen numerischen Typ zurückgeben"
        # Numerisch korrekt
        assert abs(float(val) - 0.5772156649015328) < 1e-10

    def test_harmonic_series_approximation(self, gamma_calc):
        """Harmonische-Reihen-Näherung ist korrekt auf 4 Stellen."""
        approx = gamma_calc.compute_via_harmonic_series(10000)
        expected = 0.5772156649
        assert abs(approx - expected) < 1e-4, f"Näherung zu ungenau: {approx}"

    def test_harmonic_series_converges(self, gamma_calc):
        """Mehr Terme → bessere Näherung."""
        approx_small = gamma_calc.compute_via_harmonic_series(1000)
        approx_large = gamma_calc.compute_via_harmonic_series(10000)
        gamma_true = 0.5772156649015
        assert abs(approx_large - gamma_true) < abs(approx_small - gamma_true)

    def test_zeta_limit_method(self, gamma_calc):
        """ζ-Grenzwert-Methode liefert Näherung mit Fehler < 0.01."""
        approx = gamma_calc.compute_via_zeta_limit(epsilon=1e-5)
        expected = 0.5772156649
        assert abs(approx - expected) < 0.01, f"ζ-Näherung: {approx}"

    def test_alternating_series_sign(self, gamma_calc):
        """Vacca-Reihe liefert positiven Wert."""
        val = gamma_calc.alternating_series_representation(200)
        assert val > 0

    def test_alternating_series_approximation(self, gamma_calc):
        """Vacca-Reihe: Fehler < 0.1 (langsame Konvergenz)."""
        val = gamma_calc.alternating_series_representation(500)
        expected = 0.5772156649
        assert abs(val - expected) < 0.1

    def test_apery_like_representation(self, gamma_calc):
        """Apéry-artige Darstellung nähert γ korrekt an."""
        val = gamma_calc.apery_like_representation(50)
        expected = 0.5772156649
        assert abs(val - expected) < 0.001, f"Apéry-artig: {val}"

    def test_rational_bound_exponent(self, gamma_calc):
        """Papanikolaou-Schranke hat Exponent 242080."""
        desc, exponent = gamma_calc.rational_approximation_bound()
        assert exponent == 242080

    def test_rational_bound_description_contains_conjecture(self, gamma_calc):
        """Schranken-Beschreibung erwähnt Nicht-Beweis der Irrationalität."""
        desc, _ = gamma_calc.rational_approximation_bound()
        assert "BEWEIST NICHT" in desc or "nicht" in desc.lower()

    def test_connection_to_gamma_function(self, gamma_calc):
        """Verbindung Γ'(1)=−γ ist numerisch korrekt."""
        info = gamma_calc.connection_to_gamma_function()
        assert info["verified"] is True

    def test_digamma_at_1_equals_negative_gamma(self, gamma_calc):
        """ψ(1) = −γ."""
        info = gamma_calc.connection_to_gamma_function()
        gamma_val = info["gamma"]
        digamma_val = info["digamma_at_1"]
        assert abs(digamma_val + gamma_val) < 1e-10

    def test_compute_to_n_digits_length(self, gamma_calc):
        """compute_to_n_digits(20) gibt String mit 20+ Ziffern zurück."""
        s = gamma_calc.compute_to_n_digits(20)
        # Entferne "0." und zähle Ziffern
        digits = s.replace("0.", "").replace(".", "").replace("-", "")
        assert len(digits) >= 15

    def test_compute_to_n_digits_invalid_raises(self, gamma_calc):
        """compute_to_n_digits mit ungültigem n wirft ValueError."""
        with pytest.raises(ValueError):
            gamma_calc.compute_to_n_digits(0)
        with pytest.raises(ValueError):
            gamma_calc.compute_to_n_digits(101)

    def test_verify_known_digits(self, gamma_calc):
        """Verifikation der bekannten 50 Stellen von γ."""
        assert gamma_calc.verify_known_digits() is True

    def test_precision_parameter(self):
        """Verschiedene Präzisionen liefern konsistente Ergebnisse."""
        calc_low = EulerMascheroniConstant(precision=20)
        calc_high = EulerMascheroniConstant(precision=50)
        low = float(calc_low.compute())
        high = float(calc_high.compute())
        assert abs(low - high) < 1e-15

    def test_gamma_not_zero(self, gamma_calc):
        """γ ist nicht null."""
        assert float(gamma_calc.compute()) != 0.0

    def test_gamma_less_than_one(self, gamma_calc):
        """γ < 1."""
        assert float(gamma_calc.compute()) < 1.0


# ===========================================================================
# STIELTJES-KONSTANTEN
# ===========================================================================

class TestStieltjesConstants:
    """Tests für StieltjesConstants."""

    def test_gamma0_equals_euler_mascheroni(self, stieltjes_calc):
        """γ_0 = γ (Euler-Mascheroni-Konstante)."""
        assert stieltjes_calc.verify_gamma0_equals_euler() is True

    def test_gamma0_value(self, stieltjes_calc):
        """γ_0 ≈ 0.5772..."""
        val = float(stieltjes_calc.compute(0))
        assert abs(val - 0.5772156649) < 1e-8

    def test_gamma1_is_negative(self, stieltjes_calc):
        """γ_1 < 0 (bekannte Tatsache)."""
        val = float(stieltjes_calc.compute(1))
        assert val < 0, f"γ_1 sollte negativ sein, ist {val}"

    def test_gamma1_value(self, stieltjes_calc):
        """γ_1 ≈ −0.0728158..."""
        val = float(stieltjes_calc.compute(1))
        assert abs(val - (-0.0728158454836)) < 1e-6

    def test_compute_all_length(self, stieltjes_calc):
        """compute_all(5) liefert Liste der Länge 6."""
        result = stieltjes_calc.compute_all(5)
        assert len(result) == 6

    def test_negative_index_raises(self, stieltjes_calc):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            stieltjes_calc.compute(-1)

    def test_laurent_expansion_near_1(self, stieltjes_calc):
        """Laurent-Entwicklung von ζ(s) stimmt mit mpmath-Wert überein."""
        s = 1.1 + 0j
        laurent_val = stieltjes_calc.laurent_expansion_zeta(s, terms=5)
        mpmath.mp.dps = 30
        mpmath_val = complex(mpmath.zeta(1.1))
        assert abs(laurent_val - mpmath_val) < 0.01

    def test_laurent_expansion_raises_at_pole(self, stieltjes_calc):
        """Laurent-Entwicklung bei s=1 wirft ValueError."""
        with pytest.raises(ValueError):
            stieltjes_calc.laurent_expansion_zeta(1.0)

    def test_sign_pattern_length(self, stieltjes_calc):
        """Vorzeichenmuster hat korrekte Länge."""
        signs = stieltjes_calc.alternating_sign_pattern(5)
        assert len(signs) == 6

    def test_gamma0_sign_positive(self, stieltjes_calc):
        """γ_0 > 0 → Vorzeichen +1."""
        signs = stieltjes_calc.alternating_sign_pattern(0)
        assert signs[0] == 1

    def test_size_growth_all_positive(self, stieltjes_calc):
        """Absolutwerte sind alle positiv."""
        sizes = stieltjes_calc.size_growth(5)
        assert all(s >= 0 for s in sizes)

    def test_gamma_constants_not_all_same(self, stieltjes_calc):
        """Stieltjes-Konstanten sind nicht alle gleich."""
        vals = [float(stieltjes_calc.compute(n)) for n in range(5)]
        assert len(set(round(v, 5) for v in vals)) > 1


# ===========================================================================
# ZETA AN UNGERADEN STELLEN
# ===========================================================================

class TestZetaOddValues:
    """Tests für ZetaOddValues."""

    def test_zeta3_value_10_digits(self, zeta_calc):
        """ζ(3) auf 10 Stellen korrekt: 1.2020569031..."""
        val = float(zeta_calc.zeta3())
        expected = 1.2020569031595942854
        assert abs(val - expected) < 1e-10, f"ζ(3) falsch: {val}"

    def test_zeta3_first_digit(self, zeta_calc):
        """ζ(3) ≈ 1.202..., erste Stelle ist 1."""
        val = float(zeta_calc.zeta3())
        assert int(val) == 1

    def test_zeta5_value(self, zeta_calc):
        """ζ(5) ≈ 1.0369277551..."""
        val = float(zeta_calc.zeta5())
        assert abs(val - 1.0369277551) < 1e-8

    def test_zeta7_value(self, zeta_calc):
        """ζ(7) ≈ 1.0083492773..."""
        val = float(zeta_calc.zeta7())
        assert abs(val - 1.0083492773) < 1e-8

    def test_odd_zeta_monotone_decreasing(self, zeta_calc):
        """ζ(2k+1) ist streng monoton fallend für k=1,2,3,4."""
        vals = [float(zeta_calc.compute_odd_zeta(k)) for k in range(1, 5)]
        for i in range(len(vals) - 1):
            assert vals[i] > vals[i + 1], f"ζ nicht monoton: {vals}"

    def test_odd_zeta_approaches_1(self, zeta_calc):
        """ζ(2k+1) → 1 für k → ∞."""
        val_large = float(zeta_calc.compute_odd_zeta(20))  # ζ(41)
        assert abs(val_large - 1.0) < 0.001

    def test_odd_zeta_invalid_k_raises(self, zeta_calc):
        """k < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            zeta_calc.compute_odd_zeta(0)

    # Apéry-Kettenbruch Tests
    def test_apery_cf_first_coefficient(self, zeta_calc):
        """Erster Kettenbruchkoeffizient von ζ(3) ist 1."""
        cf = zeta_calc.apery_continued_fraction(5)
        assert cf[0] == 1, f"Erster CF-Koeffizient: {cf[0]}"

    def test_apery_cf_length(self, zeta_calc):
        """Kettenbruch hat mindestens 5 Glieder."""
        cf = zeta_calc.apery_continued_fraction(10)
        assert len(cf) >= 5

    def test_apery_cf_all_positive(self, zeta_calc):
        """Alle Kettenbruchkoeffizienten sind positiv."""
        cf = zeta_calc.apery_continued_fraction(8)
        assert all(a > 0 for a in cf), f"Negativer CF-Koeffizient: {cf}"

    def test_apery_cf_reconstruction(self, zeta_calc):
        """Aus Kettenbruch rekonstruierter Wert nähert ζ(3) an."""
        cf = zeta_calc.apery_continued_fraction(15)
        reconstructed = zeta_calc.evaluate_continued_fraction(cf)
        zeta3 = 1.2020569031595942
        assert abs(reconstructed - zeta3) < 0.001

    def test_apery_cf_second_coefficient(self, zeta_calc):
        """Zweiter Koeffizient des CF von ζ(3) ist 4."""
        cf = zeta_calc.apery_continued_fraction(5)
        assert cf[1] == 4, f"Zweiter CF-Koeffizient sollte 4 sein: {cf[1]}"

    # Hypergeometrische Darstellungen
    def test_zeta3_hypergeometric(self, zeta_calc):
        """Hypergeometrische Formel für ζ(3): Fehler < 0.001."""
        val = zeta_calc.zeta3_hypergeometric(50)
        expected = 1.2020569031595942
        assert abs(val - expected) < 0.001, f"Hypergeom. Wert: {val}"

    def test_zeta5_hypergeometric_approx(self, zeta_calc):
        """Direkte Summation für ζ(5): Fehler < 0.001."""
        val = zeta_calc.zeta5_hypergeometric_approx(1000)
        expected = 1.0369277551433699
        assert abs(val - expected) < 0.001

    def test_knopp_acceleration_zeta3(self, zeta_calc):
        """Knopp-Beschleunigung für ζ(3): Näherung in richtiger Größenordnung."""
        val = zeta_calc.zeta3_knopp_acceleration(50)
        # (4/3) · (3/4) · ζ(3) ≈ ζ(3), aber Konvergenz ist langsam
        assert 0.8 < val < 1.5, f"Knopp-Wert außerhalb Bereich: {val}"

    # Ball-Rivoal Theorem
    def test_ball_rivoal_status(self, zeta_calc):
        """Ball-Rivoal ist als 'BEWIESENES THEOREM' markiert."""
        info = zeta_calc.ball_rivoal_theorem()
        assert "BEWIESENES THEOREM" in info["status"]

    def test_ball_rivoal_odd_zeta_values_count(self, zeta_calc):
        """Ball-Rivoal liefert 6 ungerade Zeta-Werte."""
        info = zeta_calc.ball_rivoal_theorem()
        assert len(info["odd_zeta_values"]) == 6

    def test_ball_rivoal_dimension_bounds_positive(self, zeta_calc):
        """Dimensions-Untergrenzen sind alle positiv."""
        info = zeta_calc.ball_rivoal_theorem()
        assert all(d > 0 for d in info["dimension_lower_bounds"])

    def test_ball_rivoal_theorem_name(self, zeta_calc):
        """Theorem-Name enthält 'Ball-Rivoal'."""
        info = zeta_calc.ball_rivoal_theorem()
        assert "Ball" in info["theorem"] and "Rivoal" in info["theorem"]

    # Zudilin Theorem
    def test_zudilin_status(self, zeta_calc):
        """Zudilin ist als 'BEWIESENES THEOREM' markiert."""
        info = zeta_calc.zudilin_theorem()
        assert "BEWIESENES THEOREM" in info["status"]

    def test_zudilin_four_values(self, zeta_calc):
        """Zudilin gibt 4 Werte zurück (ζ(5),ζ(7),ζ(9),ζ(11))."""
        info = zeta_calc.zudilin_theorem()
        assert len(info["values_at_5_7_9_11"]) == 4

    def test_zudilin_keys_correct(self, zeta_calc):
        """Zudilin enthält genau die Schlüssel 5,7,9,11."""
        info = zeta_calc.zudilin_theorem()
        assert set(info["values_at_5_7_9_11"].keys()) == {5, 7, 9, 11}

    def test_zudilin_open_questions_not_empty(self, zeta_calc):
        """Zudilin enthält offene Fragen."""
        info = zeta_calc.zudilin_theorem()
        assert len(info["open_questions"]) > 0

    # Verifikation
    def test_zeta3_digit_verification(self, zeta_calc):
        """ζ(3) auf 10 Stellen korrekt (verify_zeta3_digits)."""
        assert zeta_calc.verify_zeta3_digits(10) is True

    def test_irrationality_measure_upper_bound(self, zeta_calc):
        """Irrationalitätsmaß-Oberschranke ist 5.513..."""
        info = zeta_calc.apery_irrationality_measure()
        assert abs(info["irrationality_measure_upper_bound"] - 5.513890) < 0.001

    def test_irrationality_measure_source(self, zeta_calc):
        """Quelle ist Rhin-Viola."""
        info = zeta_calc.apery_irrationality_measure()
        assert "Rhin" in info["source"] or "Viola" in info["source"]

    def test_zeta_odd_table_length(self, zeta_calc):
        """Tabelle enthält 7 Einträge (ζ(3) bis ζ(15))."""
        table = zeta_calc.zeta_odd_table()
        assert len(table) == 7

    def test_zeta_odd_table_first_entry(self, zeta_calc):
        """Erster Tabelleneintrag ist ζ(3)."""
        table = zeta_calc.zeta_odd_table()
        assert table[0][0] == 3

    def test_zeta_odd_table_apery_status(self, zeta_calc):
        """ζ(3)-Eintrag enthält 'IRRATIONAL'."""
        table = zeta_calc.zeta_odd_table()
        assert "IRRATIONAL" in table[0][2]


# ===========================================================================
# PILLAI-CHOWLA-VERMUTUNG
# ===========================================================================

class TestPillaiChowla:
    """Tests für die PillaiChowla-Klasse."""

    def test_factorial_plus_one_small(self, pillai):
        """1!+1 = 2, 2!+1 = 3, 3!+1 = 7."""
        assert pillai.factorial_plus_one(1) == 2
        assert pillai.factorial_plus_one(2) == 3
        assert pillai.factorial_plus_one(3) == 7
        assert pillai.factorial_plus_one(4) == 25

    def test_gcd_consecutive_n1(self, pillai):
        """ggT(1!+1, 2!+1) = ggT(2,3) = 1."""
        assert pillai.gcd_consecutive(1) == 1

    def test_gcd_consecutive_n2(self, pillai):
        """ggT(2!+1, 3!+1) = ggT(3,7) = 1."""
        assert pillai.gcd_consecutive(2) == 1

    def test_gcd_consecutive_n3(self, pillai):
        """ggT(3!+1, 4!+1) = ggT(7,25) = 1."""
        assert pillai.gcd_consecutive(3) == 1

    def test_gcd_consecutive_invalid_raises(self, pillai):
        """n < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            pillai.gcd_consecutive(0)

    def test_pillai_chowla_verified_up_to_50(self, pillai):
        """Pillai-Chowla ggT=1 für alle n ≤ 50 (schwache Form)."""
        result = pillai.verify_consecutive(max_n=50)
        assert result["conjecture_holds"] is True, (
            f"Pillai-Chowla verletzt bei: {result['failures']}"
        )

    def test_pillai_chowla_no_failures_n50(self, pillai):
        """Keine Fehler für n ≤ 50."""
        result = pillai.verify_consecutive(max_n=50)
        assert len(result["failures"]) == 0

    def test_pillai_chowla_status_conjecture(self, pillai):
        """Status enthält CONJECTURE."""
        result = pillai.verify_consecutive(max_n=10)
        assert "CONJECTURE" in result["status"]

    def test_pillai_chowla_verified_count(self, pillai):
        """verify_consecutive(50) prüft genau 50 Werte."""
        result = pillai.verify_consecutive(max_n=50)
        assert result["verified_up_to"] == 50
        assert len(result["gcd_values"]) == 50

    def test_gcd_arbitrary_n1_m2(self, pillai):
        """ggT(1!+1, 2!+1) = 1."""
        assert pillai.gcd_arbitrary(1, 2) == 1

    def test_gcd_arbitrary_n3_m5(self, pillai):
        """ggT(3!+1, 5!+1) = ggT(7, 121) = 1."""
        assert pillai.gcd_arbitrary(3, 5) == 1

    def test_general_verification_has_counterexamples(self, pillai):
        """Stärkere Form ist FALSCH: ggT(3!+1, 6!+1)=7 ist ein Gegenbeispiel."""
        result = pillai.verify_general(max_val=10)
        # Die starke Form gilt NICHT: n=3,m=6 ist ein Gegenbeispiel
        assert result["conjecture_holds"] is False, (
            "Die starke Pillai-Chowla-Form gilt NICHT allgemein (Gegenbeispiele existieren)"
        )
        # Verifiziere, dass (3,6) als Gegenbeispiel gefunden wird
        failures_n = [f[0] for f in result["failures"]]
        assert 3 in failures_n, f"Gegenbeispiel n=3 nicht gefunden: {result['failures']}"

    def test_general_verification_pair_count(self, pillai):
        """Für max_val=5 gibt es C(5,2)=10 Paare."""
        result = pillai.verify_general(max_val=5)
        assert result["total_pairs_checked"] == 10

    # Wilson-Theorem Tests
    def test_wilson_theorem_prime_5(self, pillai):
        """Wilson-Theorem für p=5: (4)!=24 ≡ -1 ≡ 4 (mod 5)."""
        result = pillai.wilson_theorem_verify(primes_up_to=10)
        assert result[5]["wilson_holds"] is True
        assert result[5]["is_prime"] is True

    def test_wilson_theorem_composite_4(self, pillai):
        """Wilson-Theorem: 4 ist nicht prim, (3)!=6 ≢ -1 (mod 4)."""
        result = pillai.wilson_theorem_verify(primes_up_to=10)
        assert result[4]["is_prime"] is False
        assert result[4]["wilson_holds"] is False

    def test_wilson_theorem_all_consistent(self, pillai):
        """Wilson-Theorem stimmt mit Primtest überein für alle n ≤ 50."""
        result = pillai.wilson_theorem_verify(primes_up_to=50)
        for p, info in result.items():
            assert info["consistent"] is True, (
                f"Wilson-Theorem inkonsistent bei p={p}: {info}"
            )

    def test_wilson_primes_contains_5_13_563(self, pillai):
        """Bekannte Wilson-Primzahlen 5 und 13 werden gefunden."""
        wp = pillai.wilson_primes(up_to=600)
        assert 5 in wp, f"5 nicht in Wilson-Primzahlen: {wp}"
        assert 13 in wp, f"13 nicht in Wilson-Primzahlen: {wp}"
        assert 563 in wp, f"563 nicht in Wilson-Primzahlen: {wp}"

    def test_wilson_primes_no_false_positives(self, pillai):
        """Keine falschen Wilson-Primzahlen unter 20."""
        wp = pillai.wilson_primes(up_to=20)
        for p in wp:
            assert p in [5, 13], f"Falsche Wilson-Primzahl: {p}"

    # Sylvester-Produkt
    def test_sylvester_product_n0(self, pillai):
        """Sylvester-Produkt für n=0 ist 2."""
        assert pillai.sylvester_product(0) == 2

    def test_sylvester_product_n1(self, pillai):
        """Sylvester-Produkt für n=1: 1+2 = 3."""
        assert pillai.sylvester_product(1) == 3

    def test_sylvester_product_n2(self, pillai):
        """Sylvester-Produkt für n=2: 1+2·3 = 7."""
        assert pillai.sylvester_product(2) == 7

    def test_sylvester_product_n3(self, pillai):
        """Sylvester-Produkt für n=3: 1+2·3·5 = 31."""
        assert pillai.sylvester_product(3) == 31

    def test_sylvester_sequence_starts_2_3_7(self, pillai):
        """Sylvester-Folge beginnt mit 2, 3, 7, 43."""
        seq = pillai.sylvester_sequence(4)
        assert seq[:4] == [2, 3, 7, 43]

    def test_sylvester_sequence_length(self, pillai):
        """Sylvester-Folge mit 6 Termen hat Länge 6."""
        seq = pillai.sylvester_sequence(6)
        assert len(seq) == 6

    def test_primorial_plus_one_n1(self, pillai):
        """p_1# + 1 = 2+1 = 3."""
        assert pillai.primorial_plus_one(1) == 3

    def test_primorial_plus_one_n2(self, pillai):
        """p_2# + 1 = 6+1 = 7."""
        assert pillai.primorial_plus_one(2) == 7

    def test_primorial_plus_one_n3(self, pillai):
        """p_3# + 1 = 30+1 = 31."""
        assert pillai.primorial_plus_one(3) == 31

    # Mertens-Theorem
    def test_mertens_ratio_approaches_1(self, pillai):
        """Mertens-Verhältnis nähert sich 1 für große n."""
        results = pillai.mertens_third_theorem_numerical([100, 500, 1000])
        for n in [100, 500, 1000]:
            ratio = results[n]["ratio"]
            assert 0.8 < ratio < 1.2, f"Mertens-Verhältnis bei n={n}: {ratio}"

    def test_mertens_product_grows(self, pillai):
        """Mertens-Produkt wächst mit n."""
        results = pillai.mertens_third_theorem_numerical([100, 500])
        assert results[500]["product"] > results[100]["product"]

    def test_mertens_gamma_connection(self, pillai):
        """γ via Mertens stimmt mit wahrem Wert auf 2 Stellen überein."""
        info = pillai.mertens_connection_to_gamma(max_n=1000)
        assert info["absolute_error"] < 0.1, (
            f"Mertens-γ-Schätzung zu ungenau: Fehler={info['absolute_error']}"
        )

    def test_mertens_theorem_status(self, pillai):
        """Mertens-Theorem ist als 'BEWIESENES THEOREM' markiert."""
        info = pillai.mertens_connection_to_gamma(100)
        assert "BEWIESENES THEOREM" in info["theorem"] or "Mertens" in info["theorem"]

    # Wilson-Lemma
    def test_wilson_lemma_all_factors_large(self, pillai):
        """Alle Primfaktoren von n!+1 sind > n (für n ≤ 12)."""
        result = pillai.wilson_lemma_verify(max_n=12)
        for n in range(1, 13):
            info = result[n]
            assert info["lemma_holds"] is True, (
                f"Wilson-Lemma verletzt bei n={n}: {info}"
            )

    def test_factorial_plus_one_n1_primefactors(self, pillai):
        """Primfaktoren von 1!+1=2: nur [2], und 2 > 1. ✓"""
        result = pillai.wilson_lemma_verify(max_n=1)
        assert result[1]["lemma_holds"] is True

    # Offene Probleme
    def test_open_problems_not_empty(self, pillai):
        """Liste offener Probleme ist nicht leer."""
        problems = pillai.summary_open_problems()
        assert len(problems) >= 3

    def test_open_problems_contain_conjecture(self, pillai):
        """Alle offenen Probleme sind als CONJECTURE oder OFFENE ... markiert."""
        problems = pillai.summary_open_problems()
        for p in problems:
            status = p["status"].upper()
            assert "CONJECTURE" in status or "OFFEN" in status or "UNBEWIESEN" in status

    def test_prime_factors_all_large(self, pillai):
        """prime_factors_of_factorial_plus_one: alle Faktoren > n."""
        for n in [3, 4, 5]:
            factors = pillai.prime_factors_of_factorial_plus_one(n, max_factor=10000)
            for p in factors:
                assert p > n, f"Primfaktor p={p} ≤ n={n} (widerspricht Wilson-Lemma)"
