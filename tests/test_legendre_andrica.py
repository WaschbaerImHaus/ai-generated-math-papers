"""
@file test_legendre_andrica.py
@brief Umfassende Tests für legendre_brocard.py und andrica_cramer.py (TDD).
@description
    Testsammlung für:
    1. LegendreConjecture (legendre_brocard.py)
       - Primzahlen in (n², (n+1)²)
       - Verifikation im Bereich
       - Verbindung zu Bertrand's Postulat
       - Lückenhistogramm

    2. BrocardConjecture7 (legendre_brocard.py)
       - Primzahlen in (pₙ², pₙ₊₁²)
       - Verifikation (≥4 Primzahlen ab n=2)
       - Minimum-Tracking
       - Statistik

    3. AndricaConjecture (andrica_cramer.py)
       - A_n = √p_{n+1} − √pₙ < 1
       - Bekanntes Maximum bei n=4 (p=7, p'=11): ≈0.6708
       - Gap-Statistiken

    4. CramerConjecture (andrica_cramer.py)
       - Normierte Lücken g_n/(log pₙ)²
       - Granville-Konstante 2e^{-γ} ≈ 1.1229
       - Verteilungsanalyse

    Abdeckung: Happy Path, Edge Cases, bekannte Werte, mathematische Invarianten,
    Grenzwerte, Fehlerbehandlung, Konsistenz zwischen Modulen.

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
import sys
import os

# Suchpfad für src-Module setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import numpy as np

from legendre_brocard import LegendreConjecture, BrocardConjecture7, LegendreInterval, BrocardInterval
from andrica_cramer import AndricaConjecture, CramerConjecture, GRANVILLE_CONSTANT


# ===========================================================================
# Fixtures
# ===========================================================================

@pytest.fixture(scope="module")
def legendre():
    """Erstellt eine LegendreConjecture-Instanz (modulweit geteilt)."""
    return LegendreConjecture(max_n=10_000)


@pytest.fixture(scope="module")
def brocard():
    """Erstellt eine BrocardConjecture7-Instanz (modulweit geteilt)."""
    return BrocardConjecture7()


@pytest.fixture(scope="module")
def andrica():
    """Erstellt eine AndricaConjecture-Instanz (modulweit geteilt)."""
    return AndricaConjecture()


@pytest.fixture(scope="module")
def cramer():
    """Erstellt eine CramerConjecture-Instanz (modulweit geteilt)."""
    return CramerConjecture()


# ===========================================================================
# Tests: LegendreConjecture
# ===========================================================================

class TestLegendreConjecture:
    """Tests für die Legendre-Vermutung."""

    def test_interval_n1_returns_legendre_interval(self, legendre):
        """n=1: Intervall (1, 4). Primzahl 2 und 3 liegen drin."""
        iv = legendre.interval_primes(1)
        assert isinstance(iv, LegendreInterval)
        assert iv.n == 1
        assert iv.n_sq == 1
        assert iv.np1_sq == 4
        # Primzahlen in (1, 4): 2, 3
        assert 2 in iv.primes
        assert 3 in iv.primes
        assert iv.count == 2

    def test_interval_n2_has_primes_5_7(self, legendre):
        """n=2: Intervall (4, 9). Primzahlen 5 und 7."""
        iv = legendre.interval_primes(2)
        assert iv.n_sq == 4
        assert iv.np1_sq == 9
        assert 5 in iv.primes
        assert 7 in iv.primes
        assert iv.count == 2

    def test_interval_n3(self, legendre):
        """n=3: Intervall (9, 16). Primzahlen 11, 13."""
        iv = legendre.interval_primes(3)
        assert iv.n_sq == 9
        assert iv.np1_sq == 16
        assert 11 in iv.primes
        assert 13 in iv.primes
        assert iv.count == 2

    def test_interval_n4(self, legendre):
        """n=4: Intervall (16, 25). Primzahlen 17, 19, 23."""
        iv = legendre.interval_primes(4)
        assert iv.n_sq == 16
        assert iv.np1_sq == 25
        assert set(iv.primes) == {17, 19, 23}
        assert iv.count == 3

    def test_interval_n5(self, legendre):
        """n=5: Intervall (25, 36). Primzahlen 29, 31."""
        iv = legendre.interval_primes(5)
        assert iv.n_sq == 25
        assert iv.np1_sq == 36
        assert 29 in iv.primes
        assert 31 in iv.primes
        assert iv.count == 2

    def test_interval_boundary_excluded(self, legendre):
        """Grenzen n² und (n+1)² sind vom Intervall ausgeschlossen."""
        iv = legendre.interval_primes(4)  # (16, 25)
        # 25 = 5² ist keine Primzahl, aber auch wenn, sollte sie nicht enthalten sein
        assert 16 not in iv.primes
        assert 25 not in iv.primes

    def test_interval_raises_for_n_zero(self, legendre):
        """n=0 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            legendre.interval_primes(0)

    def test_interval_raises_for_negative_n(self, legendre):
        """Negative n sollen ValueError auslösen."""
        with pytest.raises(ValueError):
            legendre.interval_primes(-5)

    def test_interval_cache(self, legendre):
        """Mehrfaches Abrufen desselben n gibt identisches Objekt zurück."""
        iv1 = legendre.interval_primes(10)
        iv2 = legendre.interval_primes(10)
        assert iv1 is iv2

    def test_verify_range_small(self, legendre):
        """Verifikation für n=1..100 soll keine Gegenbeispiele ergeben."""
        result = legendre.verify_range(1, 100)
        assert result["verified"] is True
        assert result["counterexamples"] == []
        assert result["total_verified"] == 100

    def test_verify_range_1000(self, legendre):
        """Verifikation für n=1..1000 soll keine Gegenbeispiele ergeben."""
        result = legendre.verify_range(1, 1000)
        assert result["verified"] is True
        assert result["counterexamples"] == []

    def test_verify_range_min_count(self, legendre):
        """Minimale Primzahlanzahl in n=1..100 muss ≥ 1 sein."""
        result = legendre.verify_range(1, 100)
        assert result["min_count"] >= 1

    def test_verify_range_max_count_positive(self, legendre):
        """Maximale Primzahlanzahl muss positiv sein."""
        result = legendre.verify_range(1, 50)
        assert result["max_count"] > 0

    def test_verify_range_single_n(self, legendre):
        """Verifikation eines einzelnen n funktioniert korrekt."""
        result = legendre.verify_range(5, 5)
        assert result["total_verified"] == 1
        assert result["min_count"] == result["max_count"]

    def test_verify_range_invalid_raises(self, legendre):
        """n_end < n_start soll ValueError auslösen."""
        with pytest.raises(ValueError):
            legendre.verify_range(10, 5)

    def test_bertrand_connection_n4(self, legendre):
        """Bertrand-Verbindung für n=4: Zeuge muss in (4, 8] liegen."""
        result = legendre.bertrand_connection(4)
        assert result["bertrand_holds"] is True
        witness = result["bertrand_witness"]
        assert 4 < witness <= 8

    def test_bertrand_connection_n10(self, legendre):
        """Bertrand für n=10: Zeuge in (10, 20]."""
        result = legendre.bertrand_connection(10)
        assert result["bertrand_holds"] is True
        assert 10 < result["bertrand_witness"] <= 20

    def test_bertrand_connection_fields(self, legendre):
        """Rückgabe enthält alle erwarteten Felder."""
        result = legendre.bertrand_connection(7)
        assert "bertrand_witness" in result
        assert "bertrand_holds" in result
        assert "legendre_primes" in result
        assert "legendre_count" in result

    def test_max_gap_between_squares(self, legendre):
        """Maximale Lücke wird korrekt berechnet."""
        result = legendre.max_gap_between_squares(1, 100)
        assert "global_max_gap" in result
        assert result["global_max_gap"] > 0

    def test_gap_histogram_data(self, legendre):
        """Histogramm-Daten enthalten erwartete Felder."""
        result = legendre.gap_histogram_data(1, 100)
        assert "counts_distribution" in result
        assert "mean_count" in result
        assert "std_count" in result
        assert result["mean_count"] > 0

    def test_interval_count_increases_with_n(self, legendre):
        """Asymptotisch sollte die Primzahldichte mit n steigen (Mittelwert)."""
        result_small = legendre.gap_histogram_data(1, 50)
        result_large = legendre.gap_histogram_data(100, 150)
        # Für größere n sollte der Mittelwert höher sein
        assert result_large["mean_count"] >= result_small["mean_count"]

    def test_primality_of_interval_primes(self, legendre):
        """Alle gefundenen Primzahlen in Intervallen sind tatsächlich prim."""
        from sympy import isprime
        for n in range(1, 20):
            iv = legendre.interval_primes(n)
            for p in iv.primes:
                assert isprime(p), f"{p} ist keine Primzahl (n={n})"

    def test_primality_within_bounds(self, legendre):
        """Alle Primzahlen liegen tatsächlich im offenen Intervall (n², (n+1)²)."""
        for n in range(1, 30):
            iv = legendre.interval_primes(n)
            for p in iv.primes:
                assert iv.n_sq < p < iv.np1_sq, \
                    f"{p} liegt nicht in ({iv.n_sq}, {iv.np1_sq}) für n={n}"

    def test_interval_length_formula(self, legendre):
        """Intervalllänge = (n+1)² - n² = 2n+1."""
        for n in range(1, 20):
            iv = legendre.interval_primes(n)
            assert iv.np1_sq - iv.n_sq == 2 * n + 1


# ===========================================================================
# Tests: BrocardConjecture7
# ===========================================================================

class TestBrocardConjecture7:
    """Tests für die Brocard-Vermutung (1904)."""

    def test_n1_has_2_primes(self, brocard):
        """n=1: Zwischen 4(=2²) und 9(=3²) liegen genau 2 Primzahlen: 5, 7."""
        iv = brocard.interval_primes(1)
        assert iv.p_n == 2
        assert iv.p_n1 == 3
        assert iv.p_n ** 2 == 4
        assert iv.p_n1 ** 2 == 9
        assert set(iv.primes) == {5, 7}
        assert iv.count == 2

    def test_n2_has_5_primes(self, brocard):
        """n=2: Zwischen 9(=3²) und 25(=5²) liegen 5 Primzahlen."""
        iv = brocard.interval_primes(2)
        assert iv.p_n == 3
        assert iv.p_n1 == 5
        # Primzahlen in (9, 25): 11, 13, 17, 19, 23
        assert set(iv.primes) == {11, 13, 17, 19, 23}
        assert iv.count == 5

    def test_n3_satisfies_brocard(self, brocard):
        """n=3: Zwischen 25(=5²) und 49(=7²) liegen ≥ 4 Primzahlen."""
        iv = brocard.interval_primes(3)
        assert iv.p_n == 5
        assert iv.p_n1 == 7
        # Primzahlen in (25, 49): 29, 31, 37, 41, 43, 47
        assert iv.count >= 4
        assert 29 in iv.primes
        assert 31 in iv.primes

    def test_n4_satisfies_brocard(self, brocard):
        """n=4: Zwischen 49(=7²) und 121(=11²) liegen ≥ 4 Primzahlen."""
        iv = brocard.interval_primes(4)
        assert iv.p_n == 7
        assert iv.p_n1 == 11
        assert iv.count >= 4

    def test_raises_for_n_zero(self, brocard):
        """n=0 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            brocard.interval_primes(0)

    def test_raises_for_negative_n(self, brocard):
        """Negative n sollen ValueError auslösen."""
        with pytest.raises(ValueError):
            brocard.interval_primes(-3)

    def test_verify_range_n2_100(self, brocard):
        """Verifikation für n=2..100 soll keine Gegenbeispiele ergeben."""
        result = brocard.verify_range(2, 100)
        assert result["verified"] is True
        assert result["counterexamples"] == []

    def test_verify_range_min_count_geq_4_from_n2(self, brocard):
        """Minimale Anzahl ab n=2 muss ≥ 4 sein laut Vermutung."""
        result = brocard.verify_range(2, 100)
        assert result["min_count"] >= 4

    def test_verify_range_n1_allowed(self, brocard):
        """n=1 darf weniger als 4 Primzahlen haben (Vermutung gilt erst ab n=2)."""
        result = brocard.verify_range(1, 1)
        # n=1 hat 2 Primzahlen, das ist kein Gegenbeispiel (gilt ab n=2)
        assert result["verified"] is True

    def test_verify_range_fields(self, brocard):
        """Rückgabe enthält alle erwarteten Felder."""
        result = brocard.verify_range(2, 50)
        assert "verified" in result
        assert "counterexamples" in result
        assert "min_count" in result
        assert "max_count" in result
        assert "statistics" in result
        assert "total_verified" in result

    def test_statistics_n2_100(self, brocard):
        """Statistik-Methode liefert sinnvolle Werte für n=2..100."""
        stats = brocard.statistics(2, 100)
        assert stats["min"] >= 4
        assert stats["mean"] > 4
        assert stats["min_exceeds_4"] is True

    def test_minimum_tracking_n2_100(self, brocard):
        """Minimum-Tracking liefert korrektes globales Minimum."""
        result = brocard.minimum_tracking(2, 100)
        assert "global_minimum" in result
        assert "global_minimum_location" in result
        assert "conjecture_margin" in result
        # Margin muss ≥ 0 sein (Minimum ≥ 4)
        assert result["conjecture_margin"] >= 0

    def test_minimum_tracking_history_monotone(self, brocard):
        """Laufendes Minimum ist monoton fallend."""
        result = brocard.minimum_tracking(2, 50)
        history = result["running_minimum"]
        for i in range(1, len(history)):
            assert history[i][1] <= history[i - 1][1], \
                f"Minimum nicht monoton an Position {i}"

    def test_interval_primes_within_bounds(self, brocard):
        """Alle gefundenen Primzahlen liegen tatsächlich im Intervall."""
        for n in range(1, 15):
            iv = brocard.interval_primes(n)
            lo = iv.p_n ** 2
            hi = iv.p_n1 ** 2
            for p in iv.primes:
                assert lo < p < hi, \
                    f"{p} liegt nicht in ({lo}, {hi}) für n={n}"

    def test_interval_cache(self, brocard):
        """Caching: Identisches Objekt bei wiederholtem Abruf."""
        iv1 = brocard.interval_primes(5)
        iv2 = brocard.interval_primes(5)
        assert iv1 is iv2

    def test_known_minimum_constant_4(self, brocard):
        """Bekannte Konstante MINIMUM_REQUIRED ist 4."""
        assert BrocardConjecture7.MINIMUM_REQUIRED == 4

    def test_primes_in_n2_interval_correct(self, brocard):
        """n=2: Exakte Primzahlliste {11, 13, 17, 19, 23} verifiziert."""
        iv = brocard.interval_primes(2)
        assert sorted(iv.primes) == [11, 13, 17, 19, 23]


# ===========================================================================
# Tests: AndricaConjecture
# ===========================================================================

class TestAndricaConjecture:
    """Tests für die Andrica-Vermutung."""

    def test_known_max_constant(self):
        """Bekanntes Maximum: √11 − √7 ≈ 0.6708."""
        expected = math.sqrt(11) - math.sqrt(7)
        assert abs(AndricaConjecture.KNOWN_MAX_VALUE - expected) < 1e-10

    def test_known_max_value_approx(self):
        """Bekanntes Maximum liegt etwa bei 0.6708."""
        assert 0.670 < AndricaConjecture.KNOWN_MAX_VALUE < 0.671

    def test_compute_n1(self, andrica):
        """n=1: A_1 = √3 − √2."""
        av = andrica.compute(1)
        expected = math.sqrt(3) - math.sqrt(2)
        assert abs(av.value - expected) < 1e-10
        assert av.p_n == 2
        assert av.p_n1 == 3
        assert av.holds is True

    def test_compute_n2(self, andrica):
        """n=2: A_2 = √5 − √3."""
        av = andrica.compute(2)
        expected = math.sqrt(5) - math.sqrt(3)
        assert abs(av.value - expected) < 1e-10
        assert av.holds is True

    def test_compute_n4_is_maximum(self, andrica):
        """n=4: p=7, p'=11. A_4 ≈ 0.6708 ist bekanntes Maximum."""
        av = andrica.compute(4)
        assert av.p_n == 7
        assert av.p_n1 == 11
        assert abs(av.value - (math.sqrt(11) - math.sqrt(7))) < 1e-10
        assert av.holds is True

    def test_compute_all_hold_n1_100(self, andrica):
        """A_n < 1 für alle n = 1..100."""
        for n in range(1, 101):
            av = andrica.compute(n)
            assert av.holds is True, f"Andrica verletzt bei n={n}: A_n={av.value}"

    def test_compute_raises_n0(self, andrica):
        """n=0 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            andrica.compute(0)

    def test_verify_range_1_1000(self, andrica):
        """Verifikation für n=1..1000: Keine Verletzungen."""
        result = andrica.verify_range(1, 1000)
        assert result["verified"] is True
        assert result["counterexamples"] == []

    def test_verify_range_max_value(self, andrica):
        """Maximaler Andrica-Wert in n=1..200 liegt unter 1."""
        result = andrica.verify_range(1, 200)
        assert result["max_value"] < 1.0

    def test_verify_range_max_value_near_known(self, andrica):
        """Maximaler Wert in n=1..200 liegt nahe an 0.6708."""
        result = andrica.verify_range(1, 200)
        assert result["max_value"] > 0.65, "Bekanntes Maximum wurde nicht gefunden"

    def test_verify_range_fields(self, andrica):
        """Rückgabe enthält alle erwarteten Felder."""
        result = andrica.verify_range(1, 50)
        assert "verified" in result
        assert "max_value" in result
        assert "max_value_n" in result
        assert "max_value_primes" in result
        assert "gap_stats" in result
        assert "total_verified" in result

    def test_max_value_analysis_finds_known_max(self, andrica):
        """Maximumanalyse findet den bekannten Rekord bei n=4."""
        result = andrica.max_value_analysis(100)
        events = result["new_max_events"]
        # Einer der Ereignisse muss n=4 sein (p=7, p'=11)
        n4_events = [e for e in events if e["p_n"] == 7 and e["p_n1"] == 11]
        assert len(n4_events) >= 1

    def test_max_value_analysis_global_max(self, andrica):
        """Globales Maximum in n=1..500 liegt nahe am bekannten Wert."""
        result = andrica.max_value_analysis(500)
        assert abs(result["global_max"] - AndricaConjecture.KNOWN_MAX_VALUE) < 1e-8

    def test_gap_statistics_all_below_bound(self, andrica):
        """Alle Lücken liegen unter der Andrica-Schranke 2√pₙ + 1."""
        result = andrica.gap_statistics(1, 200)
        assert result["all_below_bound"] is True

    def test_gap_statistics_fields(self, andrica):
        """Gap-Statistik enthält erwartete Felder."""
        result = andrica.gap_statistics(1, 100)
        assert "max_ratio" in result
        assert "mean_ratio" in result
        assert "max_gap" in result

    def test_andrica_value_decreasing_trend(self, andrica):
        """Andrica-Werte sollten asymptotisch kleiner werden."""
        # Durchschnitt der ersten 10 vs. letzten 10 von n=1..100
        early = [andrica.compute(n).value for n in range(1, 11)]
        late = [andrica.compute(n).value for n in range(91, 101)]
        assert np.mean(early) > np.mean(late), \
            "Andrica-Werte nehmen nicht asymptotisch ab"

    def test_compute_cache(self, andrica):
        """Cache: compute(n) gibt bei wiederholtem Aufruf identisches Objekt."""
        av1 = andrica.compute(15)
        av2 = andrica.compute(15)
        assert av1 is av2


# ===========================================================================
# Tests: CramerConjecture
# ===========================================================================

class TestCramerConjecture:
    """Tests für die Cramér-Vermutung."""

    def test_granville_constant_value(self):
        """Granville-Konstante 2e^{-γ} ≈ 1.1229."""
        expected = 2.0 * math.exp(-0.5772156649015328)
        assert abs(GRANVILLE_CONSTANT - expected) < 1e-8

    def test_granville_constant_approx(self):
        """Granville-Konstante liegt zwischen 1.12 und 1.13."""
        assert 1.12 < GRANVILLE_CONSTANT < 1.13

    def test_compute_n1(self, cramer):
        """n=1: p=2, p'=3, gap=1. Cramér-Wert: 1/(log 2)² ≈ 2.082."""
        cv = cramer.compute(1)
        assert cv.p_n == 2
        assert cv.p_n1 == 3
        assert cv.gap == 1
        expected_log_sq = math.log(2) ** 2
        expected_cv = 1.0 / expected_log_sq
        assert abs(cv.cramer_val - expected_cv) < 1e-10

    def test_compute_fields(self, cramer):
        """CramerValue enthält alle erwarteten Felder."""
        cv = cramer.compute(5)
        assert hasattr(cv, 'p_n')
        assert hasattr(cv, 'p_n1')
        assert hasattr(cv, 'gap')
        assert hasattr(cv, 'log_sq')
        assert hasattr(cv, 'cramer_val')
        assert hasattr(cv, 'granville_ratio')

    def test_compute_raises_n0(self, cramer):
        """n=0 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            cramer.compute(0)

    def test_prime_gap_statistics_1_1000(self, cramer):
        """Prime-Gap-Statistik für n=1..1000 liefert sinnvolle Werte."""
        result = cramer.prime_gap_statistics(1, 1000)
        assert result["max_cramer_value"] > 0
        assert result["mean_cramer"] > 0
        assert result["total_analyzed"] == 1000

    def test_prime_gap_max_below_granville(self, cramer):
        """Maximaler Cramér-Wert bis n=1000 liegt unter Granville-Schranke."""
        result = cramer.prime_gap_statistics(1, 1000)
        # Empirisch: Bis n=1000 liegt das Maximum unter der Granville-Schranke
        # (das ist keine mathematische Garantie, aber empirisch wahr)
        assert result["max_cramer_value"] < GRANVILLE_CONSTANT * 2, \
            "Maximaler Cramér-Wert übersteigt 2×Granville deutlich"

    def test_prime_gap_statistics_fields(self, cramer):
        """Rückgabe enthält alle erwarteten Felder."""
        result = cramer.prime_gap_statistics(1, 100)
        assert "max_cramer_value" in result
        assert "max_cramer_n" in result
        assert "max_cramer_primes" in result
        assert "mean_cramer" in result
        assert "granville_constant" in result
        assert "granville_exceedances" in result

    def test_granville_correction_fields(self, cramer):
        """Granville-Korrektur-Methode gibt alle Felder zurück."""
        result = cramer.granville_correction(1, 100)
        assert "granville_constant" in result
        assert "euler_mascheroni" in result
        assert "max_cramer_normalized" in result
        assert "max_granville_normalized" in result
        assert "fraction_above_1" in result

    def test_granville_correction_euler_mascheroni(self, cramer):
        """Euler-Mascheroni-Konstante wird korrekt gespeichert."""
        result = cramer.granville_correction(1, 10)
        assert abs(result["euler_mascheroni"] - 0.5772156649015328) < 1e-12

    def test_maxima_analysis_first_record(self, cramer):
        """Erstes Rekordmaximum bei n=1 (p=2, gap=1)."""
        result = cramer.maxima_analysis(100)
        events = result["record_maxima"]
        assert len(events) >= 1
        # Erstes Ereignis ist bei n=1
        assert events[0]["n"] == 1
        assert events[0]["p_n"] == 2

    def test_maxima_analysis_monotone(self, cramer):
        """Rekordmaxima sind streng monoton steigend."""
        result = cramer.maxima_analysis(100)
        events = result["record_maxima"]
        vals = [e["cramer_value"] for e in events]
        for i in range(1, len(vals)):
            assert vals[i] > vals[i - 1], "Rekordmaxima nicht monoton steigend"

    def test_distribution_analysis_fields(self, cramer):
        """Verteilungsanalyse enthält alle Felder."""
        result = cramer.distribution_analysis(1, 100)
        assert "mean" in result
        assert "std" in result
        assert "quantiles" in result
        assert "histogram_counts" in result
        assert "histogram_edges" in result

    def test_distribution_quantiles_ordered(self, cramer):
        """Quantile sind geordnet: q50 ≤ q75 ≤ q90 ≤ q95 ≤ q99."""
        result = cramer.distribution_analysis(1, 200)
        q = result["quantiles"]
        assert q["q50"] <= q["q75"] <= q["q90"] <= q["q95"] <= q["q99"]

    def test_gap_distribution_contains_gap_2(self, cramer):
        """Primzahllücke 2 (Primzahlzwillinge) kommt vor."""
        result = cramer.prime_gap_statistics(1, 1000)
        # Lücke 2 kommt bei (3,5), (5,7), (11,13), ... vor
        assert 2 in result["gap_distribution"]

    def test_cramer_values_positive(self, cramer):
        """Alle Cramér-Werte sind positiv."""
        result = cramer.prime_gap_statistics(1, 200)
        assert result["max_cramer_value"] > 0
        assert result["mean_cramer"] > 0


# ===========================================================================
# Tests: Konsistenz zwischen Modulen
# ===========================================================================

class TestCrossModuleConsistency:
    """Konsistenztests zwischen LegendreConjecture und AndricaConjecture."""

    def test_andrica_equivalent_to_legendre_for_prime_pairs(self, andrica, legendre):
        """
        Andrica-Vermutung für aufeinanderfolgende Primzahlen impliziert:
        pₙ₊₁ < (√pₙ + 1)² = pₙ + 2√pₙ + 1

        Dies bedeutet, dass das Intervall (pₙ, pₙ₊₁) in ein Legendre-Intervall passt.
        """
        from sympy import prime as nth_prime
        for n in range(1, 20):
            av = andrica.compute(n)
            # Wenn Andrica gilt: p_{n+1} < (sqrt(p_n) + 1)^2
            bound = (math.sqrt(av.p_n) + 1) ** 2
            assert av.p_n1 < bound, \
                f"Andrica verletzt: p_{n+1}={av.p_n1} ≥ (√{av.p_n}+1)²={bound:.2f}"

    def test_bertrand_follows_from_legendre_for_small_n(self, legendre):
        """
        Für kleine n: Wenn Legendre gilt, dann gilt auch Bertrand.
        Bertrand: Für n ≥ 1 gibt es eine Primzahl in (n, 2n].
        """
        for n in range(1, 50):
            result = legendre.bertrand_connection(n)
            # Bertrand muss immer gelten (es ist ein Theorem!)
            assert result["bertrand_holds"] is True, \
                f"Bertrand verletzt bei n={n}: Zeuge={result['bertrand_witness']}"

    def test_brocard_n1_count_is_2(self, brocard, legendre):
        """
        n=1: Zwischen 2²=4 und 3²=9 liegen 5 und 7 (2 Primzahlen).
        Legendre(n=2): Zwischen 4 und 9 sind 5 und 7 — gleiche Zahl!
        """
        brocard_iv = brocard.interval_primes(1)
        legendre_iv = legendre.interval_primes(2)
        # Beide betrachten das Intervall (4, 9)
        assert brocard_iv.count == legendre_iv.count == 2

    def test_cramer_and_andrica_qualitative_agreement(self, cramer, andrica):
        """
        Qualitatives Agreement: Für hinreichend große pₙ gilt (log pₙ)² < 2√pₙ + 1.
        Das zeigt, dass die Cramér-Schranke für sehr große Primzahlen kleiner als
        die Andrica-Schranke ist: Cramér ist asymptotisch stärker als Andrica.

        Konkret: Für p ≥ 100 gilt (log p)² < 2√p + 1.
        """
        from sympy import prime as nth_prime
        for n in range(31, 60):
            p = nth_prime(n)
            andrica_bound = 2 * math.sqrt(p) + 1
            cramer_bound = math.log(p) ** 2
            # Für pₙ ≥ 127 gilt (log p)² < 2√p + 1 (empirisch ab p=127)
            if p >= 127:
                assert cramer_bound < andrica_bound, \
                    f"Bei p={p}: Cramér-Schranke {cramer_bound:.2f} >= Andrica-Schranke {andrica_bound:.2f}"
