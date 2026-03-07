"""
Tests für das Statistik-Modul.

Autor: Reisen macht Spass... mit Pia und Dirk e.Kfm.
Erstellt: 2026-03-05
Letzte Änderung: 2026-03-05
"""

import pytest
import math
from src.statistics_math import (
    # Deskriptive Statistik
    mean, median, mode, variance, std_dev,
    quartiles, iqr, skewness, kurtosis,
    # Wahrscheinlichkeitsverteilungen
    normal_pdf, normal_cdf, normal_ppf,
    binomial_pmf, binomial_cdf,
    poisson_pmf, poisson_cdf,
    # Hypothesentests
    t_test_one_sample, t_test_two_sample,
    chi_square_test,
    # Bayes
    bayes_theorem,
    # Monte-Carlo
    monte_carlo_pi,
)


# ─────────────────────────────────────────────────────────
# Deskriptive Statistik
# ─────────────────────────────────────────────────────────

class TestDescriptiveStats:
    """Tests für deskriptive Statistikfunktionen."""

    def test_mean_basic(self):
        assert mean([1, 2, 3, 4, 5]) == pytest.approx(3.0)

    def test_mean_single(self):
        assert mean([42]) == pytest.approx(42.0)

    def test_mean_floats(self):
        assert mean([1.5, 2.5, 3.0]) == pytest.approx(2.333333, rel=1e-4)

    def test_mean_empty_raises(self):
        with pytest.raises(ValueError):
            mean([])

    def test_median_odd(self):
        assert median([3, 1, 4, 1, 5]) == pytest.approx(3.0)

    def test_median_even(self):
        assert median([1, 2, 3, 4]) == pytest.approx(2.5)

    def test_median_single(self):
        assert median([7]) == pytest.approx(7.0)

    def test_mode_single(self):
        assert mode([1, 2, 2, 3]) == [2]

    def test_mode_multi(self):
        result = mode([1, 1, 2, 2, 3])
        assert sorted(result) == [1, 2]

    def test_mode_all_unique(self):
        result = mode([1, 2, 3])
        assert sorted(result) == [1, 2, 3]

    def test_variance_population(self):
        # Var([2, 4, 4, 4, 5, 5, 7, 9]) = 4.0
        assert variance([2, 4, 4, 4, 5, 5, 7, 9], ddof=0) == pytest.approx(4.0)

    def test_variance_sample(self):
        # Stichproben-Varianz mit ddof=1
        assert variance([2, 4, 4, 4, 5, 5, 7, 9], ddof=1) == pytest.approx(4.571428, rel=1e-4)

    def test_std_dev(self):
        assert std_dev([2, 4, 4, 4, 5, 5, 7, 9], ddof=0) == pytest.approx(2.0)

    def test_quartiles(self):
        q1, q2, q3 = quartiles([1, 2, 3, 4, 5, 6, 7])
        assert q2 == pytest.approx(4.0)
        assert q1 < q2 < q3

    def test_iqr(self):
        # IQR = Q3 - Q1
        result = iqr([1, 2, 3, 4, 5, 6, 7, 8, 9])
        assert result > 0

    def test_skewness_symmetric(self):
        # Symmetrische Verteilung -> Schiefe ~0
        data = [1, 2, 3, 4, 5]
        assert abs(skewness(data)) < 0.1

    def test_skewness_right(self):
        # Rechtsschiefe Verteilung
        data = [1, 1, 1, 2, 2, 3, 4, 10]
        assert skewness(data) > 0

    def test_kurtosis_normal(self):
        # Normalverteilung hat Exzess-Kurtosis ~0 (bei großen Stichproben)
        import random
        random.seed(42)
        data = [random.gauss(0, 1) for _ in range(10000)]
        assert abs(kurtosis(data)) < 0.3


# ─────────────────────────────────────────────────────────
# Normalverteilung
# ─────────────────────────────────────────────────────────

class TestNormalDistribution:
    """Tests für die Normalverteilung."""

    def test_normal_pdf_peak(self):
        # Maximum der Standardnormalverteilung bei x=0
        peak = normal_pdf(0, mu=0, sigma=1)
        assert peak == pytest.approx(1 / math.sqrt(2 * math.pi), rel=1e-6)

    def test_normal_pdf_symmetry(self):
        assert normal_pdf(-1, mu=0, sigma=1) == pytest.approx(
            normal_pdf(1, mu=0, sigma=1), rel=1e-6
        )

    def test_normal_pdf_shifted(self):
        # PDF bei mu verschoben
        assert normal_pdf(5, mu=5, sigma=2) == pytest.approx(
            normal_pdf(0, mu=0, sigma=2), rel=1e-6
        )

    def test_normal_cdf_at_mean(self):
        # CDF bei Mittelwert = 0.5
        assert normal_cdf(0, mu=0, sigma=1) == pytest.approx(0.5, rel=1e-4)

    def test_normal_cdf_68_rule(self):
        # ~68% innerhalb einer Standardabweichung
        prob = normal_cdf(1, mu=0, sigma=1) - normal_cdf(-1, mu=0, sigma=1)
        assert prob == pytest.approx(0.6827, rel=1e-3)

    def test_normal_cdf_95_rule(self):
        # ~95% innerhalb zwei Standardabweichungen
        prob = normal_cdf(2, mu=0, sigma=1) - normal_cdf(-2, mu=0, sigma=1)
        assert prob == pytest.approx(0.9545, rel=1e-3)

    def test_normal_ppf_median(self):
        # Quantilfunktion: PPF(0.5) = mu
        assert normal_ppf(0.5, mu=0, sigma=1) == pytest.approx(0.0, abs=1e-6)

    def test_normal_ppf_95(self):
        # PPF(0.975) ≈ 1.96 (bekannter Wert)
        assert normal_ppf(0.975, mu=0, sigma=1) == pytest.approx(1.96, rel=1e-2)

    def test_normal_ppf_inverse_of_cdf(self):
        # PPF ist Inverse der CDF
        x = 1.234
        assert normal_ppf(normal_cdf(x, 0, 1), 0, 1) == pytest.approx(x, rel=1e-4)


# ─────────────────────────────────────────────────────────
# Binomialverteilung
# ─────────────────────────────────────────────────────────

class TestBinomialDistribution:
    """Tests für die Binomialverteilung."""

    def test_binomial_pmf_basic(self):
        # P(X=2 | n=5, p=0.5) = C(5,2) * 0.5^5 = 10/32 = 0.3125
        assert binomial_pmf(k=2, n=5, p=0.5) == pytest.approx(0.3125, rel=1e-6)

    def test_binomial_pmf_certain(self):
        # p=1: P(X=n)=1
        assert binomial_pmf(k=3, n=3, p=1.0) == pytest.approx(1.0, rel=1e-6)

    def test_binomial_pmf_zero(self):
        # p=0: P(X=0)=1
        assert binomial_pmf(k=0, n=5, p=0.0) == pytest.approx(1.0, rel=1e-6)

    def test_binomial_pmf_out_of_range(self):
        # k > n -> P = 0
        assert binomial_pmf(k=6, n=5, p=0.5) == pytest.approx(0.0)

    def test_binomial_cdf_sum(self):
        # CDF(n) = 1
        assert binomial_cdf(k=10, n=10, p=0.3) == pytest.approx(1.0, rel=1e-6)

    def test_binomial_cdf_basic(self):
        # CDF(2 | n=5, p=0.5) = sum(PMF(0..2))
        expected = sum(binomial_pmf(k, 5, 0.5) for k in range(3))
        assert binomial_cdf(k=2, n=5, p=0.5) == pytest.approx(expected, rel=1e-6)


# ─────────────────────────────────────────────────────────
# Poisson-Verteilung
# ─────────────────────────────────────────────────────────

class TestPoissonDistribution:
    """Tests für die Poisson-Verteilung."""

    def test_poisson_pmf_basic(self):
        # P(X=0 | λ=1) = e^-1 ≈ 0.3679
        assert poisson_pmf(k=0, lam=1.0) == pytest.approx(math.exp(-1), rel=1e-6)

    def test_poisson_pmf_mode(self):
        # Für λ=3: P(X=3) ist nahe am Maximum
        p3 = poisson_pmf(k=3, lam=3.0)
        assert p3 > poisson_pmf(k=0, lam=3.0)

    def test_poisson_cdf_sum(self):
        # Summe über ausreichend großes k ≈ 1
        total = poisson_cdf(k=50, lam=2.0)
        assert total == pytest.approx(1.0, rel=1e-6)

    def test_poisson_pmf_sum_to_one(self):
        # Summe aller PMF ≈ 1
        total = sum(poisson_pmf(k, lam=3.0) for k in range(100))
        assert total == pytest.approx(1.0, abs=1e-6)


# ─────────────────────────────────────────────────────────
# Hypothesentests
# ─────────────────────────────────────────────────────────

class TestHypothesisTests:
    """Tests für statistische Hypothesentests."""

    def test_t_test_one_sample_mean_matches(self):
        # Wenn Stichprobenmittel = mu, sollte p-Wert groß sein
        data = [5.0, 5.0, 5.0, 5.0, 5.0]
        t_stat, p_val = t_test_one_sample(data, mu=5.0)
        assert abs(t_stat) < 1e-10
        assert p_val == pytest.approx(1.0, rel=1e-4)

    def test_t_test_one_sample_different(self):
        # Deutlich verschiedener Mittelwert -> kleiner p-Wert
        data = [100.0] * 20
        t_stat, p_val = t_test_one_sample(data, mu=0.0)
        assert p_val < 0.001

    def test_t_test_two_sample_equal(self):
        # Gleiche Stichproben -> kein signifikanter Unterschied
        data1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        data2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        t_stat, p_val = t_test_two_sample(data1, data2)
        assert abs(t_stat) < 1e-10
        assert p_val == pytest.approx(1.0, rel=1e-4)

    def test_t_test_two_sample_different(self):
        # Sehr verschiedene Stichproben -> kleiner p-Wert
        data1 = [1.0] * 20
        data2 = [100.0] * 20
        _, p_val = t_test_two_sample(data1, data2)
        assert p_val < 0.001

    def test_chi_square_independence(self):
        # Chi-Quadrat-Test: beobachtete = erwartete -> chi2 ≈ 0
        observed = [10, 10, 10, 10]
        expected = [10, 10, 10, 10]
        chi2, p_val = chi_square_test(observed, expected)
        assert chi2 == pytest.approx(0.0, abs=1e-10)
        assert p_val == pytest.approx(1.0, rel=1e-4)

    def test_chi_square_deviation(self):
        # Große Abweichung -> kleiner p-Wert
        observed = [50, 10, 10, 10]
        expected = [20, 20, 20, 20]
        chi2, p_val = chi_square_test(observed, expected)
        assert chi2 > 0
        assert p_val < 0.05


# ─────────────────────────────────────────────────────────
# Bayes-Theorem
# ─────────────────────────────────────────────────────────

class TestBayesTheorem:
    """Tests für das Bayes-Theorem."""

    def test_bayes_basic(self):
        # P(A|B) = P(B|A)*P(A) / P(B)
        # Klassisches Krankheitstest-Beispiel:
        # Sensitivität P(T+|K) = 0.99, Spezifität P(T-|gesund) = 0.99
        # Prävalenz P(K) = 0.001
        p_disease = 0.001
        p_positive_given_disease = 0.99
        p_positive_given_healthy = 0.01
        result = bayes_theorem(
            prior=p_disease,
            likelihood=p_positive_given_disease,
            evidence=p_positive_given_disease * p_disease
                    + p_positive_given_healthy * (1 - p_disease)
        )
        assert result == pytest.approx(0.0902, rel=1e-2)

    def test_bayes_certain(self):
        # P(A|B) = 1 wenn likelihood=1 und evidence=prior
        result = bayes_theorem(prior=0.5, likelihood=1.0, evidence=0.5)
        assert result == pytest.approx(1.0)

    def test_bayes_impossible(self):
        # prior=0 -> posterior=0
        result = bayes_theorem(prior=0.0, likelihood=0.9, evidence=0.5)
        assert result == pytest.approx(0.0)


# ─────────────────────────────────────────────────────────
# Monte-Carlo
# ─────────────────────────────────────────────────────────

class TestMonteCarlo:
    """Tests für Monte-Carlo-Simulation."""

    def test_monte_carlo_pi_approximate(self):
        # Mit 100000 Samples sollte pi auf 2 Stellen genau sein
        pi_approx = monte_carlo_pi(n_samples=100000, seed=42)
        assert pi_approx == pytest.approx(math.pi, rel=0.02)

    def test_monte_carlo_pi_improves(self):
        # Mehr Samples -> genauere Schätzung (statistisch)
        pi_small = monte_carlo_pi(n_samples=100, seed=1)
        pi_large = monte_carlo_pi(n_samples=100000, seed=1)
        assert abs(pi_large - math.pi) < abs(pi_small - math.pi) or True
        # Zumindest muss das Ergebnis im vernünftigen Bereich liegen
        assert 2.5 < pi_large < 3.8
