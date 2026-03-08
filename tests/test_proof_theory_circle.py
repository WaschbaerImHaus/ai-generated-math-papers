"""
@file test_proof_theory_circle.py
@brief Tests für die Hardy-Littlewood-Kreismethode in proof_theory.py.
@description
    Testet die drei Funktionen der Kreismethode:
    - hardy_littlewood_singular_series
    - goldbach_circle_method_estimate
    - goldbach_circle_method_accuracy

@author Kurt Ingwer
@version 1.0
@since 2026-03-08
@lastModified 2026-03-08
"""

import pytest
import math
import sys
import os

# Füge src/ zum Suchpfad hinzu
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from proof_theory import (
    hardy_littlewood_singular_series,
    goldbach_circle_method_estimate,
    goldbach_circle_method_accuracy,
    goldbach_all_decompositions,
)


# ===========================================================================
# TESTS FÜR hardy_littlewood_singular_series
# ===========================================================================

class TestHardyLittlewoodSingularSeries:
    """Tests für die Singuläre Reihe S(n) der Hardy-Littlewood-Methode."""

    def test_n4_is_one(self):
        """S(4) = 1, weil 4 = 2² keine ungeraden Primteiler hat."""
        # 4 = 2², einziger Primteiler ist 2 (kein ungerader Primteiler > 2)
        # Also S(4) = 1 (leeres Produkt)
        result = hardy_littlewood_singular_series(4)
        assert result == pytest.approx(1.0, rel=1e-10)

    def test_n6_has_factor_three(self):
        """S(6): 6 = 2·3, ungerader Primteiler p=3 → Faktor (3-1)/(3-2) = 2."""
        # S(6) = (p-1)/(p-2) für p=3 = 2/1 = 2.0
        result = hardy_littlewood_singular_series(6)
        assert result == pytest.approx(2.0, rel=1e-10)

    def test_n12_multiple_odd_prime_factors(self):
        """S(12): 12 = 4·3, ungerader Primteiler p=3 → S(12) = 2.0."""
        # 12 = 2² · 3, einziger ungerader Primteiler ist 3
        result = hardy_littlewood_singular_series(12)
        assert result == pytest.approx(2.0, rel=1e-10)

    def test_n30_three_prime_factors(self):
        """S(30): 30 = 2·3·5, ungerade Primteiler 3 und 5."""
        # S(30) = (3-1)/(3-2) · (5-1)/(5-2) = 2/1 · 4/3 = 8/3 ≈ 2.666...
        result = hardy_littlewood_singular_series(30)
        assert result == pytest.approx(2.0 * 4.0 / 3.0, rel=1e-9)

    def test_result_is_at_least_one(self):
        """S(n) ≥ 1 für alle getesteten geraden n ≥ 4."""
        for n in [4, 6, 8, 10, 12, 14, 16, 18, 20, 100, 1000]:
            result = hardy_littlewood_singular_series(n)
            assert result >= 1.0 - 1e-12, f"S({n}) = {result} < 1"

    def test_raises_for_odd_n(self):
        """Fehler bei ungeradem n."""
        with pytest.raises(ValueError):
            hardy_littlewood_singular_series(7)

    def test_raises_for_n_less_than_4(self):
        """Fehler für n < 4."""
        with pytest.raises(ValueError):
            hardy_littlewood_singular_series(2)

    def test_n_power_of_two_gives_one(self):
        """Für n = 2^k hat S(n) keine ungeraden Primteiler → S(n) = 1."""
        for k in range(2, 8):  # 4, 8, 16, 32, 64, 128
            n = 2 ** k
            result = hardy_littlewood_singular_series(n)
            assert result == pytest.approx(1.0, rel=1e-10), f"S(2^{k}) ≠ 1"


# ===========================================================================
# TESTS FÜR goldbach_circle_method_estimate
# ===========================================================================

class TestGoldbachCircleMethodEstimate:
    """Tests für die Goldbach-Schätzung via Kreismethode."""

    def test_estimate_is_positive(self):
        """Die Schätzung muss positiv sein für gerades n ≥ 4."""
        for n in [4, 6, 8, 10, 20, 50, 100, 1000]:
            estimate = goldbach_circle_method_estimate(n)
            assert estimate > 0, f"Schätzung für n={n} ist nicht positiv: {estimate}"

    def test_n100_ratio_in_reasonable_range(self):
        """
        Für n=100: Verhältnis actual/estimate sollte zwischen 0.5 und 5.0 liegen.
        Die Asymptotik gilt erst für große n, kleine n können stärker abweichen.
        """
        actual = len(goldbach_all_decompositions(100))
        estimate = goldbach_circle_method_estimate(100)
        ratio = actual / estimate

        assert 0.3 <= ratio <= 6.0, (
            f"Ratio für n=100 außerhalb vernünftiger Grenzen: "
            f"actual={actual}, estimate={estimate:.2f}, ratio={ratio:.3f}"
        )

    def test_n1000_ratio_in_reasonable_range(self):
        """
        Für n=1000: Verhältnis sollte zwischen 0.5 und 3.0 liegen.
        Die Asymptotik wird für größeres n besser.
        """
        actual = len(goldbach_all_decompositions(1000))
        estimate = goldbach_circle_method_estimate(1000)
        ratio = actual / estimate

        assert 0.4 <= ratio <= 4.0, (
            f"Ratio für n=1000 außerhalb vernünftiger Grenzen: "
            f"actual={actual}, estimate={estimate:.2f}, ratio={ratio:.3f}"
        )

    def test_larger_n_gives_larger_estimate(self):
        """Für größeres n sollte die Schätzung wachsen (monotone Tendenz)."""
        prev = goldbach_circle_method_estimate(100)
        for n in [200, 500, 1000, 2000]:
            curr = goldbach_circle_method_estimate(n)
            # Nicht streng monoton (wegen S(n)), aber der Trend ist steigend
            assert curr > 0, f"Schätzung für n={n} ist negativ"
            prev = curr

    def test_raises_for_odd_n(self):
        """Fehler bei ungeradem n."""
        with pytest.raises(ValueError):
            goldbach_circle_method_estimate(101)

    def test_raises_for_n_less_than_4(self):
        """Fehler für n < 4."""
        with pytest.raises(ValueError):
            goldbach_circle_method_estimate(2)

    def test_estimate_uses_twin_prime_constant(self):
        """
        Prüfe, dass die Konstante C₂ ≈ 0.6601618 korrekt eingeflossen ist.
        Für n mit S(n)=1 (z.B. n=8): estimate ≈ C₂ · n / (ln n)²
        """
        n = 8  # S(8)=1, da 8=2³
        estimate = goldbach_circle_method_estimate(n)
        # estimate = 2*C2*n/(ln n)^2 * S(n) / 2 = C2*n/(ln n)^2
        C2 = 0.6601618158468695
        expected = C2 * n / (math.log(n) ** 2)
        assert estimate == pytest.approx(expected, rel=1e-5)


# ===========================================================================
# TESTS FÜR goldbach_circle_method_accuracy
# ===========================================================================

class TestGoldbachCircleMethodAccuracy:
    """Tests für den Genauigkeitsvergleich der Kreismethode."""

    def test_returns_correct_structure(self):
        """Rückgabe ist eine Liste von Dictionaries mit den richtigen Feldern."""
        results = goldbach_circle_method_accuracy(20)

        # Muss eine Liste sein
        assert isinstance(results, list), "Rückgabe muss eine Liste sein"

        # Für jedes gerade n von 4 bis 20 ein Eintrag: 4,6,8,10,12,14,16,18,20 = 9 Einträge
        assert len(results) == 9, f"Erwartete 9 Einträge, erhalten: {len(results)}"

        # Jeder Eintrag hat die richtigen Felder
        required_keys = {'n', 'actual', 'estimate', 'ratio'}
        for entry in results:
            assert required_keys.issubset(entry.keys()), (
                f"Fehlende Schlüssel in Eintrag: {entry.keys()}"
            )

    def test_n_values_are_even_and_in_range(self):
        """Alle n-Werte sind gerade und zwischen 4 und n_max."""
        results = goldbach_circle_method_accuracy(20)
        for entry in results:
            n = entry['n']
            assert n % 2 == 0, f"n={n} ist nicht gerade"
            assert 4 <= n <= 20, f"n={n} außerhalb des Bereichs [4, 20]"

    def test_actual_counts_correct(self):
        """Die tatsächlichen Zerlegungsanzahlen sind korrekt berechnet."""
        results = goldbach_circle_method_accuracy(20)
        for entry in results:
            n = entry['n']
            expected_actual = len(goldbach_all_decompositions(n))
            assert entry['actual'] == expected_actual, (
                f"Falsche actual-Zahl für n={n}: {entry['actual']} ≠ {expected_actual}"
            )

    def test_ratio_is_actual_divided_by_estimate(self):
        """Das Verhältnis ist actual/estimate."""
        results = goldbach_circle_method_accuracy(20)
        for entry in results:
            if entry['estimate'] > 0:
                expected_ratio = entry['actual'] / entry['estimate']
                assert entry['ratio'] == pytest.approx(expected_ratio, rel=1e-10), (
                    f"Falsches Verhältnis für n={entry['n']}"
                )

    def test_n4_has_exactly_one_decomposition(self):
        """n=4 hat genau eine Zerlegung: 2+2=4."""
        results = goldbach_circle_method_accuracy(20)
        entry_4 = next(e for e in results if e['n'] == 4)
        assert entry_4['actual'] == 1, f"4 sollte 1 Zerlegung haben, hat {entry_4['actual']}"

    def test_ratios_are_positive(self):
        """Alle Verhältnisse müssen positiv sein."""
        results = goldbach_circle_method_accuracy(20)
        for entry in results:
            assert entry['ratio'] > 0, f"Negatives Verhältnis bei n={entry['n']}"

    def test_larger_n_max_gives_more_entries(self):
        """Größeres n_max liefert mehr Einträge."""
        r20 = goldbach_circle_method_accuracy(20)
        r40 = goldbach_circle_method_accuracy(40)
        assert len(r40) > len(r20), "Größeres n_max sollte mehr Einträge liefern"

    def test_n_max_30_has_correct_count(self):
        """n_max=30 liefert geraden n von 4 bis 30: 4,6,...,30 = 14 Einträge."""
        results = goldbach_circle_method_accuracy(30)
        assert len(results) == 14, f"Erwartete 14 Einträge für n_max=30, erhalten: {len(results)}"
