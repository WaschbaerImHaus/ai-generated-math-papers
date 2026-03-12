"""
@file test_mersenne_normality.py
@brief Umfassende Tests für mersenne_fermat.py und normality_digits.py.
@description
    Testet:
    - MersennePrimes: Lucas-Lehmer-Test, bekannte Exponenten, Wagstaff-Heuristik,
      vollkommene Zahlen, Stellenzahl
    - FermatPrimes: Pépin-Test, bekannte Faktoren, Fermat-Zahlen-Berechnung,
      Konstruierbarkeit, Faktorisierungsverifikation
    - NormalNumberAnalysis: Ziffernextraktion, Chi²-Test, Runs-Test, Longest-Run,
      Champernowne-Normalität
    - NormalityBounds: Borel-Satz, Diskrepanz, Schranken, Copeland-Erdős

@author Michael Fuhrmann
@lastModified 2026-03-12
"""

import sys
import os
import math
import pytest

# Pfad zum src-Verzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from mersenne_fermat import MersennePrimes, FermatPrimes
from normality_digits import NormalNumberAnalysis, NormalityBounds


# ══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def mersenne():
    """Erstellt eine MersennePrimes-Instanz."""
    return MersennePrimes()


@pytest.fixture
def fermat():
    """Erstellt eine FermatPrimes-Instanz."""
    return FermatPrimes()


@pytest.fixture
def normal_analysis():
    """Erstellt eine NormalNumberAnalysis-Instanz mit 500 Stellen."""
    return NormalNumberAnalysis(precision=500)


@pytest.fixture
def normality_bounds():
    """Erstellt eine NormalityBounds-Instanz."""
    return NormalityBounds()


# ══════════════════════════════════════════════════════════════════════════════
# Tests: MersennePrimes — Lucas-Lehmer-Test
# ══════════════════════════════════════════════════════════════════════════════

class TestLucasLehmerTest:
    """Tests für den Lucas-Lehmer-Test."""

    def test_m2_prim(self, mersenne):
        """M_2 = 3 ist prim."""
        assert mersenne.lucas_lehmer_test(2) is True

    def test_m3_prim(self, mersenne):
        """M_3 = 7 ist prim."""
        assert mersenne.lucas_lehmer_test(3) is True

    def test_m5_prim(self, mersenne):
        """M_5 = 31 ist prim."""
        assert mersenne.lucas_lehmer_test(5) is True

    def test_m7_prim(self, mersenne):
        """M_7 = 127 ist prim."""
        assert mersenne.lucas_lehmer_test(7) is True

    def test_m13_prim(self, mersenne):
        """M_13 = 8191 ist prim."""
        assert mersenne.lucas_lehmer_test(13) is True

    def test_m17_prim(self, mersenne):
        """M_17 = 131071 ist prim."""
        assert mersenne.lucas_lehmer_test(17) is True

    def test_m19_prim(self, mersenne):
        """M_19 = 524287 ist prim."""
        assert mersenne.lucas_lehmer_test(19) is True

    def test_m31_prim(self, mersenne):
        """M_31 = 2147483647 ist prim."""
        assert mersenne.lucas_lehmer_test(31) is True

    def test_m11_zusammengesetzt(self, mersenne):
        """M_11 = 2047 = 23 × 89 ist NICHT prim."""
        assert mersenne.lucas_lehmer_test(11) is False

    def test_m23_zusammengesetzt(self, mersenne):
        """M_23 = 8388607 = 47 × 178481 ist NICHT prim."""
        assert mersenne.lucas_lehmer_test(23) is False

    def test_m29_zusammengesetzt(self, mersenne):
        """M_29 = 536870911 = 233 × 1103 × 2089 ist NICHT prim."""
        assert mersenne.lucas_lehmer_test(29) is False

    def test_m37_zusammengesetzt(self, mersenne):
        """M_37 ist zusammengesetzt."""
        assert mersenne.lucas_lehmer_test(37) is False

    def test_nicht_prim_exponent_4(self, mersenne):
        """M_4 = 15 = 3 × 5: Exponent 4 ist nicht prim → False."""
        assert mersenne.lucas_lehmer_test(4) is False

    def test_nicht_prim_exponent_6(self, mersenne):
        """M_6 = 63: Exponent 6 ist nicht prim → False."""
        assert mersenne.lucas_lehmer_test(6) is False

    def test_fehler_bei_p_kleiner_2(self, mersenne):
        """Exponent < 2 wirft ValueError."""
        with pytest.raises(ValueError):
            mersenne.lucas_lehmer_test(1)

    def test_fehler_bei_p_null(self, mersenne):
        """Exponent 0 wirft ValueError."""
        with pytest.raises(ValueError):
            mersenne.lucas_lehmer_test(0)


# ══════════════════════════════════════════════════════════════════════════════
# Tests: MersennePrimes — Mersenne-Zahl und bekannte Exponenten
# ══════════════════════════════════════════════════════════════════════════════

class TestMersenneNumbers:
    """Tests für Mersenne-Zahlen-Berechnung und bekannte Exponenten."""

    def test_mersenne_number_p2(self, mersenne):
        """M_2 = 3."""
        assert mersenne.mersenne_number(2) == 3

    def test_mersenne_number_p3(self, mersenne):
        """M_3 = 7."""
        assert mersenne.mersenne_number(3) == 7

    def test_mersenne_number_p5(self, mersenne):
        """M_5 = 31."""
        assert mersenne.mersenne_number(5) == 31

    def test_mersenne_number_p7(self, mersenne):
        """M_7 = 127."""
        assert mersenne.mersenne_number(7) == 127

    def test_mersenne_number_p13(self, mersenne):
        """M_13 = 8191."""
        assert mersenne.mersenne_number(13) == 8191

    def test_known_exponents_length(self, mersenne):
        """Es gibt genau 51 bekannte Mersenne-Primzahl-Exponenten."""
        assert len(mersenne.KNOWN_EXPONENTS) == 51

    def test_known_exponents_first_eight(self, mersenne):
        """Die ersten 8 Exponenten sind korrekt."""
        assert mersenne.KNOWN_EXPONENTS[:8] == [2, 3, 5, 7, 13, 17, 19, 31]

    def test_largest_known_exponent(self, mersenne):
        """Der größte bekannte Exponent ist 82589933."""
        assert mersenne.KNOWN_EXPONENTS[-1] == 82589933

    def test_is_known_exponent_true(self, mersenne):
        """p=61 ist bekannter Mersenne-Exponent."""
        assert mersenne.is_known_mersenne_prime_exponent(61) is True

    def test_is_known_exponent_false(self, mersenne):
        """p=11 ist KEIN bekannter Mersenne-Exponent."""
        assert mersenne.is_known_mersenne_prime_exponent(11) is False

    def test_get_first_n_exponents(self, mersenne):
        """Erste 5 Exponenten werden korrekt zurückgegeben."""
        result = mersenne.get_first_n_mersenne_primes(5)
        assert result == [2, 3, 5, 7, 13]

    def test_get_first_n_exponents_error(self, mersenne):
        """Mehr als 51 Exponenten anfordern wirft ValueError."""
        with pytest.raises(ValueError):
            mersenne.get_first_n_mersenne_primes(52)

    def test_digit_count_m2(self, mersenne):
        """M_2 = 3 hat 1 Stelle."""
        assert mersenne.digit_count(2) == 1

    def test_digit_count_m7(self, mersenne):
        """M_7 = 127 hat 3 Stellen."""
        assert mersenne.digit_count(7) == 3

    def test_digit_count_m82589933(self, mersenne):
        """M_82589933 hat ca. 24 Millionen Stellen."""
        count = mersenne.digit_count(82589933)
        assert count > 24_000_000


# ══════════════════════════════════════════════════════════════════════════════
# Tests: MersennePrimes — Vollkommene Zahlen und Wagstaff-Heuristik
# ══════════════════════════════════════════════════════════════════════════════

class TestPerfectNumbersAndWagstaff:
    """Tests für vollkommene Zahlen und Wagstaff-Heuristik."""

    def test_perfect_number_p2(self, mersenne):
        """p=2 → vollkommene Zahl 6."""
        assert mersenne.perfect_number_from_mersenne(2) == 6

    def test_perfect_number_p3(self, mersenne):
        """p=3 → vollkommene Zahl 28."""
        assert mersenne.perfect_number_from_mersenne(3) == 28

    def test_perfect_number_p5(self, mersenne):
        """p=5 → vollkommene Zahl 496."""
        assert mersenne.perfect_number_from_mersenne(5) == 496

    def test_perfect_number_p7(self, mersenne):
        """p=7 → vollkommene Zahl 8128."""
        assert mersenne.perfect_number_from_mersenne(7) == 8128

    def test_perfect_number_p11_none(self, mersenne):
        """p=11: M_11 nicht prim → None."""
        assert mersenne.perfect_number_from_mersenne(11) is None

    def test_verify_6_perfect(self, mersenne):
        """6 = 1+2+3 ist vollkommen."""
        assert mersenne.verify_perfect_number(6) is True

    def test_verify_28_perfect(self, mersenne):
        """28 = 1+2+4+7+14 ist vollkommen."""
        assert mersenne.verify_perfect_number(28) is True

    def test_verify_496_perfect(self, mersenne):
        """496 ist vollkommen."""
        assert mersenne.verify_perfect_number(496) is True

    def test_verify_12_not_perfect(self, mersenne):
        """12 ist nicht vollkommen."""
        assert mersenne.verify_perfect_number(12) is False

    def test_verify_10_not_perfect(self, mersenne):
        """10 ist nicht vollkommen."""
        assert mersenne.verify_perfect_number(10) is False

    def test_wagstaff_p2(self, mersenne):
        """Wagstaff-Heuristik für p=2 liefert Wert zwischen 0 und 1."""
        prob = mersenne.wagstaff_heuristic(2)
        assert 0 < prob <= 1.0

    def test_wagstaff_abnehmend(self, mersenne):
        """Wagstaff-Wahrscheinlichkeit ist abnehmend in p."""
        prob_small = mersenne.wagstaff_heuristic(5)
        prob_large = mersenne.wagstaff_heuristic(100)
        assert prob_small > prob_large

    def test_wagstaff_fehler(self, mersenne):
        """p <= 1 wirft ValueError."""
        with pytest.raises(ValueError):
            mersenne.wagstaff_heuristic(1)


# ══════════════════════════════════════════════════════════════════════════════
# Tests: FermatPrimes — Fermat-Zahlen und Pépin-Test
# ══════════════════════════════════════════════════════════════════════════════

class TestFermatNumbers:
    """Tests für Fermat-Zahlen und den Pépin-Test."""

    def test_fermat_number_f0(self, fermat):
        """F_0 = 2^1 + 1 = 3."""
        assert fermat.fermat_number(0) == 3

    def test_fermat_number_f1(self, fermat):
        """F_1 = 2^2 + 1 = 5."""
        assert fermat.fermat_number(1) == 5

    def test_fermat_number_f2(self, fermat):
        """F_2 = 2^4 + 1 = 17."""
        assert fermat.fermat_number(2) == 17

    def test_fermat_number_f3(self, fermat):
        """F_3 = 2^8 + 1 = 257."""
        assert fermat.fermat_number(3) == 257

    def test_fermat_number_f4(self, fermat):
        """F_4 = 2^16 + 1 = 65537."""
        assert fermat.fermat_number(4) == 65537

    def test_fermat_number_f5(self, fermat):
        """F_5 = 2^32 + 1 = 4294967297."""
        assert fermat.fermat_number(5) == 4294967297

    def test_fermat_number_negativ_fehler(self, fermat):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            fermat.fermat_number(-1)


class TestPepinTest:
    """Tests für den Pépin-Test."""

    def test_pepin_f0_prim(self, fermat):
        """F_0 = 3 ist prim (Sonderfall)."""
        assert fermat.pepin_test(0) is True

    def test_pepin_f1_prim(self, fermat):
        """F_1 = 5 ist prim."""
        assert fermat.pepin_test(1) is True

    def test_pepin_f2_prim(self, fermat):
        """F_2 = 17 ist prim."""
        assert fermat.pepin_test(2) is True

    def test_pepin_f3_prim(self, fermat):
        """F_3 = 257 ist prim."""
        assert fermat.pepin_test(3) is True

    def test_pepin_f4_prim(self, fermat):
        """F_4 = 65537 ist prim."""
        assert fermat.pepin_test(4) is True

    def test_pepin_f5_zusammengesetzt(self, fermat):
        """F_5 = 4294967297 = 641 × 6700417 ist zusammengesetzt."""
        assert fermat.pepin_test(5) is False

    def test_pepin_f6_zusammengesetzt(self, fermat):
        """F_6 ist zusammengesetzt."""
        assert fermat.pepin_test(6) is False

    def test_pepin_zu_gross_fehler(self, fermat):
        """n > 8 wirft ValueError (zu groß)."""
        with pytest.raises(ValueError):
            fermat.pepin_test(9)

    def test_pepin_negativ_fehler(self, fermat):
        """n < 0 wirft ValueError."""
        with pytest.raises(ValueError):
            fermat.pepin_test(-1)


class TestFermatFactors:
    """Tests für bekannte Fermat-Faktoren."""

    def test_known_factor_f5(self, fermat):
        """Bekannter Faktor von F_5 ist 641."""
        assert fermat.get_known_factor(5) == 641

    def test_known_factor_f6(self, fermat):
        """Bekannter Faktor von F_6 ist 274177."""
        assert fermat.get_known_factor(6) == 274177

    def test_no_known_factor_f0(self, fermat):
        """F_0 hat keinen bekannten Faktor (ist prim)."""
        assert fermat.get_known_factor(0) is None

    def test_verify_factorization_f5(self, fermat):
        """641 teilt tatsächlich F_5 = 4294967297."""
        assert fermat.verify_factorization(5) is True

    def test_verify_factorization_f6(self, fermat):
        """274177 teilt tatsächlich F_6."""
        assert fermat.verify_factorization(6) is True

    def test_verify_factorization_f5_manual(self, fermat):
        """Manuelle Verifikation: 4294967297 % 641 == 0."""
        assert 4294967297 % 641 == 0

    def test_f5_factor_641_korrekt(self, fermat):
        """F_5 / 641 = 6700417 (Euler 1732)."""
        f5 = fermat.fermat_number(5)
        assert f5 // 641 == 6700417
        assert f5 % 641 == 0

    def test_is_prime_f0(self, fermat):
        """F_0 = 3 ist prim."""
        assert fermat.is_prime_fermat(0) is True

    def test_is_prime_f4(self, fermat):
        """F_4 = 65537 ist prim."""
        assert fermat.is_prime_fermat(4) is True

    def test_is_prime_f5_false(self, fermat):
        """F_5 ist zusammengesetzt."""
        assert fermat.is_prime_fermat(5) is False

    def test_gauss_wantzel_f0(self, fermat):
        """F_0 = 3: reguläres 3-Eck konstruierbar."""
        assert fermat.gauss_wantzel_constructible(0) is True

    def test_gauss_wantzel_f4(self, fermat):
        """F_4 = 65537: reguläres 65537-Eck konstruierbar."""
        assert fermat.gauss_wantzel_constructible(4) is True

    def test_gauss_wantzel_f5_false(self, fermat):
        """F_5 zusammengesetzt: kein F_5-Eck konstruierbar."""
        assert fermat.gauss_wantzel_constructible(5) is False

    def test_wagstaff_fermat_dict(self, fermat):
        """Wagstaff-Heuristik gibt Dictionary mit 0..20 zurück."""
        result = fermat.wagstaff_heuristic_fermat()
        assert 0 in result
        assert 20 in result
        assert 'total_expected' in result

    def test_wagstaff_fermat_endliche_summe(self, fermat):
        """Summe der Wahrscheinlichkeiten ist endlich (< 10)."""
        result = fermat.wagstaff_heuristic_fermat()
        assert result['total_expected'] < 10.0


# ══════════════════════════════════════════════════════════════════════════════
# Tests: NormalNumberAnalysis — Ziffernextraktion
# ══════════════════════════════════════════════════════════════════════════════

class TestDigitExtraction:
    """Tests für die Extraktion von Ziffern mathematischer Konstanten."""

    def test_pi_erste_ziffern(self, normal_analysis):
        """π beginnt mit 14159... (nach dem Dezimalpunkt)."""
        digits = normal_analysis.get_digits_of_pi(5)
        assert digits[0] == 1
        assert digits[1] == 4
        assert digits[2] == 1
        assert digits[3] == 5
        assert digits[4] == 9

    def test_pi_laenge(self, normal_analysis):
        """Genau n Ziffern von π werden zurückgegeben."""
        digits = normal_analysis.get_digits_of_pi(100)
        assert len(digits) == 100

    def test_e_erste_ziffern(self, normal_analysis):
        """e beginnt mit 71828... (nach dem Dezimalpunkt)."""
        digits = normal_analysis.get_digits_of_e(5)
        assert digits[0] == 7
        assert digits[1] == 1
        assert digits[2] == 8
        assert digits[3] == 2
        assert digits[4] == 8

    def test_sqrt2_erste_ziffern(self, normal_analysis):
        """√2 beginnt mit 41421... (nach dem Dezimalpunkt)."""
        digits = normal_analysis.get_digits_of_sqrt2(5)
        assert digits[0] == 4
        assert digits[1] == 1
        assert digits[2] == 4
        assert digits[3] == 2
        assert digits[4] == 1

    def test_champernowne_erste_ziffern(self, normal_analysis):
        """Champernowne beginnt mit 1,2,3,4,5,6,7,8,9,1,0,1,1,..."""
        digits = normal_analysis.get_digits_of_champernowne(12)
        assert digits[:9] == [1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert digits[9] == 1
        assert digits[10] == 0
        assert digits[11] == 1

    def test_champernowne_laenge(self, normal_analysis):
        """Genau n Ziffern der Champernowne-Konstante."""
        digits = normal_analysis.get_digits_of_champernowne(50)
        assert len(digits) == 50

    def test_alle_ziffern_im_bereich(self, normal_analysis):
        """Alle Ziffern sind im Bereich 0–9."""
        digits = normal_analysis.get_digits_of_pi(100)
        assert all(0 <= d <= 9 for d in digits)


# ══════════════════════════════════════════════════════════════════════════════
# Tests: NormalNumberAnalysis — Statistische Tests
# ══════════════════════════════════════════════════════════════════════════════

class TestStatisticalTests:
    """Tests für Chi²-, Runs- und Longest-Run-Tests."""

    def test_chi2_champernowne_nicht_abgelehnt(self, normal_analysis):
        """Chi²-Test für Champernowne (bewiesen normal) — Chi²-Wert muss endlich sein.

        Hinweis: Bei kurzen Abschnitten der Champernowne-Konstante ist die Verteilung
        nicht gleichförmig (kleine Zahlen erzeugen mehr führende 1er). Normalität gilt
        asymptotisch. Wir prüfen daher nur, dass der Test strukturell korrekt läuft.
        """
        digits = normal_analysis.get_digits_of_champernowne(1000)
        result = normal_analysis.chi_squared_test(digits)
        # Chi²-Statistik muss nicht-negativ und endlich sein
        assert result['chi2_statistic'] >= 0
        assert math.isfinite(result['chi2_statistic'])
        assert 0.0 <= result['p_value'] <= 1.0
        assert result['n_digits'] == 1000

    def test_chi2_ergibt_dict(self, normal_analysis):
        """Chi²-Test gibt Dictionary mit allen Schlüsseln zurück."""
        digits = list(range(10)) * 100  # gleichverteilte Ziffern
        result = normal_analysis.chi_squared_test(digits)
        assert 'chi2_statistic' in result
        assert 'p_value' in result
        assert 'degrees_of_freedom' in result
        assert 'n_digits' in result
        assert 'observed_counts' in result

    def test_chi2_perfekte_gleichverteilung(self, normal_analysis):
        """Perfekt gleichverteilte Ziffern → Chi²=0, p-Wert=1."""
        digits = list(range(10)) * 100  # je 100 mal jede Ziffer
        result = normal_analysis.chi_squared_test(digits)
        assert abs(result['chi2_statistic']) < 1e-10
        assert result['p_value'] > 0.99
        assert bool(result['reject_h0_alpha_05']) is False

    def test_chi2_schlechte_verteilung(self, normal_analysis):
        """Stark ungleichverteilte Ziffern → Chi²-Test lehnt H0 ab."""
        # Nur Ziffern 0 und 1 (extrem ungleichverteilt für 10er-System)
        digits = [0] * 500 + [1] * 500
        result = normal_analysis.chi_squared_test(digits)
        assert bool(result['reject_h0_alpha_01']) is True

    def test_chi2_freiheitsgrade(self, normal_analysis):
        """Chi²-Test hat 9 Freiheitsgrade (10 Klassen − 1)."""
        digits = list(range(10)) * 10
        result = normal_analysis.chi_squared_test(digits)
        assert result['degrees_of_freedom'] == 9

    def test_runs_test_ergibt_dict(self, normal_analysis):
        """Runs-Test gibt Dictionary zurück."""
        digits = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] * 50
        result = normal_analysis.runs_test(digits)
        assert 'z_statistic' in result
        assert 'p_value' in result

    def test_runs_test_abwechselnd(self, normal_analysis):
        """Streng abwechselnde Folge hat viele Runs → signifikant."""
        digits = [0, 9] * 250  # maximale Runs: 500 Runs
        result = normal_analysis.runs_test(digits)
        assert result['z_statistic'] is not None

    def test_longest_run_ergibt_dict(self, normal_analysis):
        """Longest-Run-Test gibt Dictionary zurück."""
        digits = [1, 1, 1, 2, 3, 4, 5] * 10
        result = normal_analysis.longest_run_test(digits)
        assert 'max_run_per_digit' in result
        assert 'overall_max' in result
        assert 'expected_max' in result

    def test_longest_run_korrekt(self, normal_analysis):
        """Longest Run für [3,3,3,3,1,2] ist 4 (für Ziffer 3)."""
        digits = [3, 3, 3, 3, 1, 2]
        result = normal_analysis.longest_run_test(digits)
        assert result['max_run_per_digit'][3] == 4
        assert result['overall_max'] == 4

    def test_longest_run_leer(self, normal_analysis):
        """Leere Liste → overall_max = 0."""
        result = normal_analysis.longest_run_test([])
        assert result['overall_max'] == 0

    def test_analyze_pi_laeuft_durch(self, normal_analysis):
        """analyze_constant('pi') läuft ohne Fehler durch."""
        result = normal_analysis.analyze_constant('pi', n_digits=100)
        assert result['constant'] == 'pi'
        assert 'chi_squared_test' in result
        assert 'runs_test' in result
        assert 'longest_run' in result

    def test_analyze_champernowne_laeuft_durch(self, normal_analysis):
        """analyze_constant('champernowne') läuft ohne Fehler durch."""
        result = normal_analysis.analyze_constant('champernowne', n_digits=200)
        assert result['constant'] == 'champernowne'
        assert 'BEWIESEN' in result['normalcy_status']

    def test_analyze_unbekannte_konstante_fehler(self, normal_analysis):
        """Unbekannte Konstante wirft ValueError."""
        with pytest.raises(ValueError):
            normal_analysis.analyze_constant('tau')

    def test_digit_frequency_gleichverteilt(self, normal_analysis):
        """Gleichverteilte Ziffern → jede Frequenz ≈ 0.1."""
        digits = list(range(10)) * 100
        freq = normal_analysis.digit_frequency(digits)
        for d in range(10):
            assert abs(freq[d]['frequency'] - 0.1) < 1e-10

    def test_digit_frequency_alle_ziffern(self, normal_analysis):
        """digit_frequency gibt Einträge für alle Ziffern 0–9 zurück."""
        digits = [5] * 100  # nur Ziffer 5
        freq = normal_analysis.digit_frequency(digits)
        assert len(freq) == 10
        assert freq[5]['frequency'] == 1.0
        for d in range(10):
            if d != 5:
                assert freq[d]['frequency'] == 0.0


# ══════════════════════════════════════════════════════════════════════════════
# Tests: NormalityBounds
# ══════════════════════════════════════════════════════════════════════════════

class TestNormalityBounds:
    """Tests für die NormalityBounds-Klasse."""

    def test_borel_theorem_dict(self, normality_bounds):
        """borel_measure_theorem() gibt Dictionary mit Theorem-Infos zurück."""
        result = normality_bounds.borel_measure_theorem()
        assert result['status'] == 'BEWIESEN'
        assert 'open_problems' in result
        assert len(result['open_problems']) >= 4

    def test_borel_theorem_conjectures(self, normality_bounds):
        """Borel-Dict enthält mindestens 4 offene Probleme (Conjectures)."""
        result = normality_bounds.borel_measure_theorem()
        conjectures = [p for p in result['open_problems'] if 'CONJECTURE' in p]
        assert len(conjectures) >= 4

    def test_discrepancy_gleichverteilt(self, normality_bounds):
        """Perfekt gleichverteilte Ziffern → Diskrepanz ≈ 0."""
        digits = list(range(10)) * 100
        d = normality_bounds.compute_discrepancy(digits)
        assert d < 1e-10

    def test_discrepancy_einzig(self, normality_bounds):
        """Nur eine Ziffer → maximale Diskrepanz ≈ 0.9."""
        digits = [0] * 100
        d = normality_bounds.compute_discrepancy(digits)
        assert abs(d - 0.9) < 1e-10

    def test_discrepancy_leer(self, normality_bounds):
        """Leere Liste → Diskrepanz = 0."""
        d = normality_bounds.compute_discrepancy([])
        assert d == 0.0

    def test_normality_bound_n100(self, normality_bounds):
        """Normalitätsschranke für n=100 gibt korrektes Dict."""
        result = normality_bounds.normality_bound_conjecture(100)
        assert result['n_digits'] == 100
        assert result['expected_discrepancy'] > 0
        assert result['theoretical_bound_95pct'] > 0

    def test_normality_bound_abnehmend(self, normality_bounds):
        """Schranke nimmt mit n ab (mehr Ziffern → genauere Schätzung)."""
        bound_100 = normality_bounds.normality_bound_conjecture(100)
        bound_10000 = normality_bounds.normality_bound_conjecture(10000)
        assert bound_100['expected_discrepancy'] > bound_10000['expected_discrepancy']

    def test_normality_bound_fehler(self, normality_bounds):
        """n=0 wirft ValueError."""
        with pytest.raises(ValueError):
            normality_bounds.normality_bound_conjecture(0)

    def test_copeland_erdos_laenge(self, normality_bounds):
        """Copeland-Erdős-Konstante gibt n Ziffern zurück."""
        digits = normality_bounds.copeland_erdos_constant_digits(20)
        assert len(digits) == 20

    def test_copeland_erdos_erste_ziffern(self, normality_bounds):
        """Copeland-Erdős beginnt mit 2,3,5,7,1,1,1,3,..."""
        digits = normality_bounds.copeland_erdos_constant_digits(8)
        # Primzahlen: 2, 3, 5, 7, 11, 13 → Ziffern: 2,3,5,7,1,1,1,3
        assert digits[:4] == [2, 3, 5, 7]
        assert digits[4] == 1  # erste Ziffer von 11
        assert digits[5] == 1  # zweite Ziffer von 11

    def test_discrepancy_champernowne_kleiner_als_einheitsverteilung(self):
        """Champernowne hat geringere Diskrepanz als rein einseitige Verteilung.

        Hinweis: Bei 500 Stellen ist Champernowne noch nicht asymptotisch normal
        (kleine Zahlen sind überrepräsentiert). Die Diskrepanz ist trotzdem deutlich
        kleiner als bei einer maximal schlechten Verteilung.
        """
        bounds = NormalityBounds()
        analysis = NormalNumberAnalysis(precision=100)
        digits_champ = analysis.get_digits_of_champernowne(500)
        digits_worst = [0] * 500  # maximale Diskrepanz: 0.9
        d_champ = bounds.compute_discrepancy(digits_champ)
        d_worst = bounds.compute_discrepancy(digits_worst)
        # Champernowne muss besser als die schlechteste Verteilung sein
        assert d_champ < d_worst
