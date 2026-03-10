"""
@file test_estimation_theory.py
@brief Umfassende Tests für das estimation_theory-Modul.
@description
    Testet alle Klassen und Methoden des Schätztheorie-Moduls:
    - MaximumLikelihoodEstimator: MLE für Normal, Exp, Poisson, Bernoulli, Binomial
    - MethodOfMoments: MOM für Normal, Gamma, Beta
    - SufficiencyTheory: Suffizienz, Rao-Blackwell
    - CramerRaoBound: CRB, Effizienz, UMVUE, Fisher-Information
    - ConfidenceIntervals: KI für Mittelwert, Varianz, Anteil, Bootstrap, Exponential
    - HypothesisTesting: Neyman-Pearson, UMP, Fehlerarten, Gütefunktion, SPRT
    - DecisionTheory: Verlustfunktionen, Risikofunktion, Bayes-Risiko, Minimax

    Testet auch Edge Cases: leere Stichproben, n=1, Randparameter.

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import math
import sys
import os
import pytest
import numpy as np
from scipy import stats

# Quellverzeichnis in den Suchpfad aufnehmen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Zu testendes Modul
from estimation_theory import (
    MaximumLikelihoodEstimator,
    MethodOfMoments,
    SufficiencyTheory,
    CramerRaoBound,
    ConfidenceIntervals,
    HypothesisTesting,
    DecisionTheory,
)

# Reproduzierbarkeit sicherstellen
RNG = np.random.default_rng(12345)


# =============================================================================
# Hilfsfunktionen für Tests
# =============================================================================

def make_normal_samples(n: int, mu: float = 2.0, sigma: float = 1.5, seed: int = 0) -> np.ndarray:
    """Erzeugt Normalverteilungs-Stichprobe für Tests."""
    return np.random.default_rng(seed).normal(mu, sigma, n)


def make_exponential_samples(n: int, lam: float = 2.0, seed: int = 1) -> np.ndarray:
    """Erzeugt Exponentialverteilungs-Stichprobe."""
    return np.random.default_rng(seed).exponential(1.0 / lam, n)


def make_poisson_samples(n: int, lam: float = 3.5, seed: int = 2) -> np.ndarray:
    """Erzeugt Poisson-Stichprobe."""
    return np.random.default_rng(seed).poisson(lam, n).astype(float)


def make_bernoulli_samples(n: int, p: float = 0.4, seed: int = 3) -> np.ndarray:
    """Erzeugt Bernoulli-Stichprobe."""
    return np.random.default_rng(seed).binomial(1, p, n).astype(float)


def make_beta_samples(n: int, a: float = 2.0, b: float = 5.0, seed: int = 4) -> np.ndarray:
    """Erzeugt Beta-Stichprobe."""
    return np.random.default_rng(seed).beta(a, b, n)


def make_gamma_samples(n: int, alpha: float = 3.0, beta: float = 2.0, seed: int = 5) -> np.ndarray:
    """Erzeugt Gamma-Stichprobe (scale = beta)."""
    return np.random.default_rng(seed).gamma(alpha, beta, n)


# =============================================================================
# 1. Tests für MaximumLikelihoodEstimator
# =============================================================================

class TestMLENormal:
    """Tests für mle_normal."""

    def test_mle_normal_mean_large_sample(self):
        """MLE-Mittelwert konvergiert gegen wahres μ bei großer Stichprobe."""
        samples = make_normal_samples(10000, mu=3.0, sigma=2.0)
        mu_hat, _ = MaximumLikelihoodEstimator.mle_normal(samples)
        assert abs(mu_hat - 3.0) < 0.1

    def test_mle_normal_variance_large_sample(self):
        """MLE-Varianz konvergiert gegen wahres σ² (leicht verzerrt, da /n)."""
        samples = make_normal_samples(10000, mu=0.0, sigma=2.0)
        _, sigma2_hat = MaximumLikelihoodEstimator.mle_normal(samples)
        # MLE-Varianz: Teiler n, ergibt σ²·(n-1)/n ≈ σ²
        assert abs(sigma2_hat - 4.0) < 0.2

    def test_mle_normal_returns_tuple(self):
        """Rückgabe ist ein Tupel (float, float)."""
        result = MaximumLikelihoodEstimator.mle_normal(np.array([1.0, 2.0, 3.0]))
        assert isinstance(result, tuple) and len(result) == 2

    def test_mle_normal_single_sample(self):
        """n=1: Varianz = 0 (MLE mit Teiler n)."""
        mu_hat, sigma2_hat = MaximumLikelihoodEstimator.mle_normal(np.array([5.0]))
        assert mu_hat == 5.0
        assert sigma2_hat == 0.0

    def test_mle_normal_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_normal(np.array([]))

    def test_mle_normal_negative_values(self):
        """Negative Werte sind erlaubt (Normalverteilung hat Träger ℝ)."""
        samples = np.array([-3.0, -1.0, 0.0, 1.0, 3.0])
        mu_hat, _ = MaximumLikelihoodEstimator.mle_normal(samples)
        assert math.isclose(mu_hat, 0.0, abs_tol=1e-10)


class TestMLEExponential:
    """Tests für mle_exponential."""

    def test_mle_exp_large_sample(self):
        """MLE-λ konvergiert gegen wahres λ."""
        samples = make_exponential_samples(10000, lam=3.0)
        lam_hat = MaximumLikelihoodEstimator.mle_exponential(samples)
        assert abs(lam_hat - 3.0) < 0.2

    def test_mle_exp_is_reciprocal_mean(self):
        """λ̂ = 1/x̄ (analytische Eigenschaft)."""
        samples = np.array([0.5, 1.0, 1.5, 2.0])
        lam_hat = MaximumLikelihoodEstimator.mle_exponential(samples)
        assert math.isclose(lam_hat, 1.0 / np.mean(samples), rel_tol=1e-10)

    def test_mle_exp_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_exponential(np.array([]))

    def test_mle_exp_zero_mean_raises(self):
        """Mittelwert = 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_exponential(np.array([0.0, 0.0]))

    def test_mle_exp_positive(self):
        """λ̂ ist immer positiv für positive Stichproben."""
        samples = np.array([0.1, 0.2, 0.5])
        lam_hat = MaximumLikelihoodEstimator.mle_exponential(samples)
        assert lam_hat > 0


class TestMLEPoisson:
    """Tests für mle_poisson."""

    def test_mle_poisson_large_sample(self):
        """MLE-λ konvergiert gegen wahres λ."""
        samples = make_poisson_samples(10000, lam=4.0)
        lam_hat = MaximumLikelihoodEstimator.mle_poisson(samples)
        assert abs(lam_hat - 4.0) < 0.1

    def test_mle_poisson_equals_mean(self):
        """λ̂ = x̄ (Stichprobenmittelwert)."""
        samples = np.array([2.0, 3.0, 5.0, 4.0])
        lam_hat = MaximumLikelihoodEstimator.mle_poisson(samples)
        assert math.isclose(lam_hat, float(np.mean(samples)), rel_tol=1e-10)

    def test_mle_poisson_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_poisson(np.array([]))

    def test_mle_poisson_zeros(self):
        """Alle Nullen: λ̂ = 0."""
        samples = np.array([0.0, 0.0, 0.0])
        lam_hat = MaximumLikelihoodEstimator.mle_poisson(samples)
        assert lam_hat == 0.0


class TestMLEBernoulli:
    """Tests für mle_bernoulli."""

    def test_mle_bernoulli_large_sample(self):
        """p̂ konvergiert gegen wahres p."""
        samples = make_bernoulli_samples(10000, p=0.7)
        p_hat = MaximumLikelihoodEstimator.mle_bernoulli(samples)
        assert abs(p_hat - 0.7) < 0.05

    def test_mle_bernoulli_all_zeros(self):
        """Alle Beobachtungen 0: p̂ = 0."""
        samples = np.array([0.0, 0.0, 0.0])
        p_hat = MaximumLikelihoodEstimator.mle_bernoulli(samples)
        assert p_hat == 0.0

    def test_mle_bernoulli_all_ones(self):
        """Alle Beobachtungen 1: p̂ = 1."""
        samples = np.array([1.0, 1.0, 1.0])
        p_hat = MaximumLikelihoodEstimator.mle_bernoulli(samples)
        assert p_hat == 1.0

    def test_mle_bernoulli_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_bernoulli(np.array([]))


class TestMLEBinomial:
    """Tests für mle_binomial."""

    def test_mle_binomial_basic(self):
        """p̂ ist korrekt für bekannte Stichprobe."""
        # 10 Beobachtungen von Bin(20, p), durchschnittlich ~10 Erfolge
        samples = np.array([10.0, 12.0, 8.0, 11.0, 9.0])
        p_hat = MaximumLikelihoodEstimator.mle_binomial(samples, n=20)
        assert math.isclose(p_hat, float(np.mean(samples)) / 20.0, rel_tol=1e-10)

    def test_mle_binomial_n_invalid_raises(self):
        """n ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_binomial(np.array([5.0]), n=0)

    def test_mle_binomial_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.mle_binomial(np.array([]), n=10)

    def test_mle_binomial_range(self):
        """p̂ liegt in [0, 1]."""
        samples = np.array([3.0, 7.0, 5.0])
        p_hat = MaximumLikelihoodEstimator.mle_binomial(samples, n=10)
        assert 0.0 <= p_hat <= 1.0


class TestLogLikelihood:
    """Tests für log_likelihood_normal."""

    def test_log_likelihood_at_mle(self):
        """Log-Likelihood ist maximal bei MLE-Parametern."""
        samples = make_normal_samples(100, mu=2.0, sigma=1.5)
        mu_hat, sigma2_hat = MaximumLikelihoodEstimator.mle_normal(samples)
        ll_mle = MaximumLikelihoodEstimator.log_likelihood_normal(
            samples, mu_hat, math.sqrt(sigma2_hat)
        )
        # Vergleichswert mit leicht anderen Parametern
        ll_other = MaximumLikelihoodEstimator.log_likelihood_normal(
            samples, mu_hat + 0.5, math.sqrt(sigma2_hat)
        )
        assert ll_mle > ll_other

    def test_log_likelihood_negative_sigma_raises(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.log_likelihood_normal(
                np.array([1.0, 2.0]), mu=0.0, sigma=-1.0
            )

    def test_log_likelihood_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.log_likelihood_normal(
                np.array([]), mu=0.0, sigma=1.0
            )

    def test_fisher_information_normal(self):
        """Fisher-Information der N(μ,σ²) bzgl. μ ist 1/σ²."""
        sigma = 2.0
        fi = MaximumLikelihoodEstimator.fisher_information_normal(sigma)
        assert math.isclose(fi, 1.0 / sigma**2, rel_tol=1e-10)

    def test_fisher_information_invalid_sigma(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            MaximumLikelihoodEstimator.fisher_information_normal(0.0)


# =============================================================================
# 2. Tests für MethodOfMoments
# =============================================================================

class TestSampleMoments:
    """Tests für sample_moments und central_moments."""

    def test_first_moment_is_mean(self):
        """1. Moment = arithmetischer Mittelwert."""
        samples = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        m1 = MethodOfMoments.sample_moments(samples, k=1)
        assert math.isclose(m1, float(np.mean(samples)), rel_tol=1e-10)

    def test_second_moment_correct(self):
        """2. Rohmoment = E[X²] = Var + Mean²."""
        samples = np.array([1.0, 2.0, 3.0])
        m2 = MethodOfMoments.sample_moments(samples, k=2)
        expected = float(np.mean(samples**2))
        assert math.isclose(m2, expected, rel_tol=1e-10)

    def test_central_moment_2_is_variance(self):
        """2. zentrales Moment = Stichprobenvarianz (Teiler n)."""
        samples = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        cm2 = MethodOfMoments.central_moments(samples, k=2)
        expected = float(np.var(samples, ddof=0))
        assert math.isclose(cm2, expected, rel_tol=1e-10)

    def test_central_moment_1_is_zero(self):
        """1. zentrales Moment = 0 (Mittelwert zentriert)."""
        samples = make_normal_samples(1000)
        cm1 = MethodOfMoments.central_moments(samples, k=1)
        assert abs(cm1) < 1e-10

    def test_invalid_k_raises(self):
        """k < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.sample_moments(np.array([1.0, 2.0]), k=0)

    def test_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.sample_moments(np.array([]), k=1)


class TestMOMNormal:
    """Tests für mom_normal."""

    def test_mom_normal_agrees_with_mle(self):
        """MOM-Schätzer stimmt mit MLE für Normalverteilung überein."""
        samples = make_normal_samples(1000, mu=5.0, sigma=2.0)
        mu_mom, sigma2_mom = MethodOfMoments.mom_normal(samples)
        mu_mle, sigma2_mle = MaximumLikelihoodEstimator.mle_normal(samples)
        assert math.isclose(mu_mom, mu_mle, rel_tol=1e-10)
        assert math.isclose(sigma2_mom, sigma2_mle, rel_tol=1e-10)

    def test_mom_normal_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.mom_normal(np.array([]))


class TestMOMGamma:
    """Tests für mom_gamma."""

    def test_mom_gamma_large_sample(self):
        """α und β werden korrekt geschätzt."""
        alpha_true, beta_true = 4.0, 2.0
        samples = make_gamma_samples(5000, alpha=alpha_true, beta=beta_true)
        alpha_hat, beta_hat = MethodOfMoments.mom_gamma(samples)
        assert abs(alpha_hat - alpha_true) < 0.5
        assert abs(beta_hat - beta_true) < 0.3

    def test_mom_gamma_positive_params(self):
        """Geschätzte Parameter sind positiv."""
        samples = make_gamma_samples(100)
        alpha_hat, beta_hat = MethodOfMoments.mom_gamma(samples)
        assert alpha_hat > 0 and beta_hat > 0

    def test_mom_gamma_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.mom_gamma(np.array([]))


class TestMOMBeta:
    """Tests für mom_beta."""

    def test_mom_beta_large_sample(self):
        """α und β werden korrekt geschätzt für Beta(2, 5)."""
        samples = make_beta_samples(5000, a=2.0, b=5.0)
        alpha_hat, beta_hat = MethodOfMoments.mom_beta(samples)
        assert abs(alpha_hat - 2.0) < 0.5
        assert abs(beta_hat - 5.0) < 1.0

    def test_mom_beta_invalid_mean_raises(self):
        """Mittelwert außerhalb (0,1) löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.mom_beta(np.array([2.0, 3.0, 4.0]))

    def test_mom_beta_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            MethodOfMoments.mom_beta(np.array([]))


# =============================================================================
# 3. Tests für SufficiencyTheory
# =============================================================================

class TestSufficiencyTheory:
    """Tests für Suffizienztheorie."""

    def test_sufficient_statistic_true_with_mean(self):
        """T = x̄ ist suffizient für μ."""
        samples = make_normal_samples(50)
        T = float(np.mean(samples))
        result = SufficiencyTheory.is_sufficient_statistic_normal_mean(samples, T)
        assert result["sufficient"] is True

    def test_sufficient_statistic_false_with_wrong_T(self):
        """T ≠ x̄ ist nicht suffizient."""
        samples = make_normal_samples(50)
        T_wrong = float(np.mean(samples)) + 1.0  # absichtlich falsch
        result = SufficiencyTheory.is_sufficient_statistic_normal_mean(samples, T_wrong)
        assert result["sufficient"] is False

    def test_sufficient_statistic_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            SufficiencyTheory.is_sufficient_statistic_normal_mean(np.array([]), T=0.0)

    def test_factorization_criterion_demo(self):
        """Faktorisierungssatz gibt korrektes Dictionary zurück."""
        samples = make_normal_samples(20)
        T_values = np.linspace(-5, 5, 10)
        result = SufficiencyTheory.factorization_criterion_demo(samples, T_values)
        assert result["factorization_holds"] is True
        assert "T_observed" in result
        assert math.isclose(result["T_observed"], float(np.sum(samples)), rel_tol=1e-10)

    def test_complete_statistic_demo(self):
        """Vollständige Statistik gibt informatives Dictionary zurück."""
        result = SufficiencyTheory.complete_statistic_demo()
        assert result["complete"] is True
        assert "definition" in result
        assert "lehmann_scheffe" in result

    def test_rao_blackwell(self):
        """Rao-Blackwell gibt Dictionary mit Schätzern zurück."""
        samples = make_normal_samples(50, mu=3.0)
        result = SufficiencyTheory.rao_blackwell(
            estimator_fn=np.mean,
            sufficient_statistic_fn=np.mean,
            samples=samples,
        )
        assert "original_estimate" in result
        assert "rao_blackwell_estimate" in result
        # Beide sollten nahe am wahren μ=3 sein
        assert abs(result["original_estimate"] - 3.0) < 0.5

    def test_rao_blackwell_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            SufficiencyTheory.rao_blackwell(np.mean, np.mean, np.array([]))


# =============================================================================
# 4. Tests für CramerRaoBound
# =============================================================================

class TestCramerRaoBound:
    """Tests für Cramér-Rao-Ungleichung."""

    def test_crb_normal_mean_formula(self):
        """CRB für N(μ,σ²) ist σ²/n."""
        crb = CramerRaoBound.cramer_rao_normal_mean(n=100, sigma=2.0)
        assert math.isclose(crb, 4.0 / 100.0, rel_tol=1e-10)

    def test_crb_normal_mean_decreases_with_n(self):
        """Größeres n → kleinere CRB (mehr Information)."""
        crb_small = CramerRaoBound.cramer_rao_normal_mean(n=10, sigma=1.0)
        crb_large = CramerRaoBound.cramer_rao_normal_mean(n=1000, sigma=1.0)
        assert crb_large < crb_small

    def test_crb_normal_invalid_n_raises(self):
        """n ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.cramer_rao_normal_mean(n=0, sigma=1.0)

    def test_crb_normal_invalid_sigma_raises(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.cramer_rao_normal_mean(n=10, sigma=0.0)

    def test_crb_poisson_formula(self):
        """CRB für Poisson-λ ist λ/n."""
        crb = CramerRaoBound.cramer_rao_poisson(n=50, lam=3.0)
        assert math.isclose(crb, 3.0 / 50.0, rel_tol=1e-10)

    def test_crb_poisson_invalid_lam_raises(self):
        """λ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.cramer_rao_poisson(n=10, lam=0.0)

    def test_efficiency_efficient_estimator(self):
        """Effizienz = 1 wenn Varianz = CRB (effizienter Schätzer)."""
        crb = 0.04
        eff = CramerRaoBound.efficiency(variance=crb, cramer_rao_bound=crb)
        assert math.isclose(eff, 1.0, rel_tol=1e-10)

    def test_efficiency_less_than_one(self):
        """Effizienz ≤ 1 für nicht-effiziente Schätzer."""
        crb = 0.04
        eff = CramerRaoBound.efficiency(variance=0.08, cramer_rao_bound=crb)
        assert math.isclose(eff, 0.5, rel_tol=1e-10)
        assert eff <= 1.0

    def test_efficiency_invalid_variance_raises(self):
        """Varianz ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.efficiency(variance=0.0, cramer_rao_bound=0.1)

    def test_umvue_normal_demo(self):
        """UMVUE-Demo gibt korrekte Schätzer zurück."""
        samples = make_normal_samples(100, mu=4.0, sigma=2.0)
        result = CramerRaoBound.umvue_normal_demo(samples)
        assert abs(result["mu_umvue"] - 4.0) < 0.3
        assert abs(result["sigma2_umvue"] - 4.0) < 0.5
        assert result["n"] == 100

    def test_umvue_normal_single_sample(self):
        """n=1: sigma2_umvue ist nan (nicht definiert)."""
        result = CramerRaoBound.umvue_normal_demo(np.array([3.0]))
        assert math.isnan(result["sigma2_umvue"])
        assert result["mu_umvue"] == 3.0

    def test_umvue_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.umvue_normal_demo(np.array([]))

    def test_fisher_information_numerical(self):
        """Numerische Fisher-Information stimmt mit analytischem Wert überein."""
        samples = make_normal_samples(50, mu=2.0, sigma=1.5)
        sigma = 1.5
        # Log-Likelihood als Funktion von μ (σ fest)
        def ll_mu(mu):
            return MaximumLikelihoodEstimator.log_likelihood_normal(samples, mu, sigma)
        fi_numerical = CramerRaoBound.fisher_information(ll_mu, theta=2.0)
        # Analytisch: I(μ) = n/σ² für die Stichprobens-Log-Likelihood
        fi_analytical = len(samples) / sigma**2
        assert abs(fi_numerical - fi_analytical) / fi_analytical < 0.01

    def test_fisher_information_invalid_epsilon_raises(self):
        """Epsilon ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            CramerRaoBound.fisher_information(lambda t: -t**2, theta=0.0, epsilon=0.0)

    def test_sample_mean_achieves_crb(self):
        """x̄ ist effizienter Schätzer: empirische Var ≈ CRB bei großem n."""
        sigma = 2.0
        n = 5000
        n_trials = 200
        # Empirische Varianz von x̄ über viele Stichproben
        means = [float(np.mean(np.random.default_rng(i).normal(0, sigma, n)))
                 for i in range(n_trials)]
        empirical_var = float(np.var(means))
        crb = CramerRaoBound.cramer_rao_normal_mean(n, sigma)
        # Empirische Varianz sollte nahe an CRB liegen
        assert abs(empirical_var - crb) / crb < 0.2


# =============================================================================
# 5. Tests für ConfidenceIntervals
# =============================================================================

class TestCIForMean:
    """Tests für ci_normal_mean."""

    def test_ci_contains_true_mean(self):
        """95%-KI enthält den wahren Mittelwert (für große n fast sicher)."""
        samples = make_normal_samples(1000, mu=5.0, sigma=2.0)
        lower, upper, _ = ConfidenceIntervals.ci_normal_mean(samples, confidence=0.95)
        assert lower < 5.0 < upper

    def test_ci_narrower_with_more_samples(self):
        """Größere Stichprobe → engeres KI."""
        samples_small = make_normal_samples(20, mu=0.0, sigma=1.0)
        samples_large = make_normal_samples(2000, mu=0.0, sigma=1.0, seed=99)
        lo_s, hi_s, _ = ConfidenceIntervals.ci_normal_mean(samples_small)
        lo_l, hi_l, _ = ConfidenceIntervals.ci_normal_mean(samples_large)
        assert (hi_l - lo_l) < (hi_s - lo_s)

    def test_ci_single_sample_is_infinite(self):
        """n=1: KI ist (-∞, +∞)."""
        lo, hi, method = ConfidenceIntervals.ci_normal_mean(np.array([3.0]))
        assert lo == float("-inf") and hi == float("inf")

    def test_ci_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_normal_mean(np.array([]))

    def test_ci_invalid_confidence_raises(self):
        """Ungültiges Konfidenzniveau löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_normal_mean(np.array([1.0, 2.0]), confidence=1.5)

    def test_ci_higher_confidence_wider(self):
        """Höheres Konfidenzniveau → breiteres KI."""
        samples = make_normal_samples(100)
        lo_95, hi_95, _ = ConfidenceIntervals.ci_normal_mean(samples, confidence=0.95)
        lo_99, hi_99, _ = ConfidenceIntervals.ci_normal_mean(samples, confidence=0.99)
        assert (hi_99 - lo_99) > (hi_95 - lo_95)


class TestCIForVariance:
    """Tests für ci_normal_variance."""

    def test_ci_variance_contains_true_sigma2(self):
        """95%-KI für σ² enthält wahres σ²."""
        samples = make_normal_samples(200, mu=0.0, sigma=3.0)
        lower, upper = ConfidenceIntervals.ci_normal_variance(samples)
        true_sigma2 = 9.0
        assert lower < true_sigma2 < upper

    def test_ci_variance_lower_less_upper(self):
        """Untere Grenze < obere Grenze."""
        samples = make_normal_samples(50)
        lo, hi = ConfidenceIntervals.ci_normal_variance(samples)
        assert lo < hi

    def test_ci_variance_n1_raises(self):
        """n=1 löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_normal_variance(np.array([1.0]))

    def test_ci_variance_invalid_confidence_raises(self):
        """Ungültiges Niveau löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_normal_variance(np.array([1.0, 2.0]), confidence=0.0)


class TestCIProportion:
    """Tests für ci_proportion (Wilson-Methode)."""

    def test_ci_proportion_contains_true_p(self):
        """Wilson-KI enthält wahres p."""
        # p = 0.6, n = 500, successes ≈ 300
        n, successes = 500, 300
        lo, hi = ConfidenceIntervals.ci_proportion(successes, n, confidence=0.95)
        assert lo < 0.6 < hi

    def test_ci_proportion_all_success(self):
        """Alle Erfolge: obere Grenze = 1."""
        lo, hi = ConfidenceIntervals.ci_proportion(50, 50)
        assert hi == 1.0

    def test_ci_proportion_no_success(self):
        """Keine Erfolge: untere Grenze = 0."""
        lo, hi = ConfidenceIntervals.ci_proportion(0, 50)
        assert lo == 0.0

    def test_ci_proportion_invalid_n_raises(self):
        """n ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_proportion(5, 0)

    def test_ci_proportion_successes_gt_n_raises(self):
        """successes > n löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_proportion(15, 10)

    def test_ci_proportion_in_unit_interval(self):
        """KI-Grenzen liegen in [0, 1]."""
        lo, hi = ConfidenceIntervals.ci_proportion(3, 10)
        assert 0.0 <= lo <= hi <= 1.0


class TestCIBootstrap:
    """Tests für ci_bootstrap."""

    def test_ci_bootstrap_median(self):
        """Bootstrap-KI für den Median enthält den wahren Median."""
        samples = make_normal_samples(500, mu=5.0, sigma=1.0)
        lo, hi = ConfidenceIntervals.ci_bootstrap(
            samples, np.median, confidence=0.95, n_bootstrap=1000
        )
        assert lo < 5.0 < hi

    def test_ci_bootstrap_mean(self):
        """Bootstrap-KI für den Mittelwert."""
        samples = make_normal_samples(200, mu=3.0)
        lo, hi = ConfidenceIntervals.ci_bootstrap(
            samples, np.mean, confidence=0.95, n_bootstrap=500
        )
        assert lo < 3.0 < hi

    def test_ci_bootstrap_reproducible(self):
        """Gleiches Seed → gleiche Ergebnisse."""
        samples = make_normal_samples(100)
        lo1, hi1 = ConfidenceIntervals.ci_bootstrap(samples, np.mean, seed=42)
        lo2, hi2 = ConfidenceIntervals.ci_bootstrap(samples, np.mean, seed=42)
        assert lo1 == lo2 and hi1 == hi2

    def test_ci_bootstrap_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_bootstrap(np.array([]), np.mean)

    def test_ci_bootstrap_invalid_confidence_raises(self):
        """Ungültiges Niveau löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_bootstrap(np.array([1.0, 2.0]), np.mean, confidence=0.0)


class TestCIExponential:
    """Tests für ci_exponential."""

    def test_ci_exponential_contains_true_lambda(self):
        """KI enthält wahres λ."""
        lam_true = 2.0
        samples = make_exponential_samples(500, lam=lam_true)
        lo, hi = ConfidenceIntervals.ci_exponential(samples, confidence=0.95)
        assert lo < lam_true < hi

    def test_ci_exponential_lower_less_upper(self):
        """Untere Grenze < obere Grenze."""
        samples = make_exponential_samples(100)
        lo, hi = ConfidenceIntervals.ci_exponential(samples)
        assert lo < hi

    def test_ci_exponential_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            ConfidenceIntervals.ci_exponential(np.array([]))


# =============================================================================
# 6. Tests für HypothesisTesting
# =============================================================================

class TestNeymanPearson:
    """Tests für neyman_pearson_lemma_demo."""

    def test_neyman_pearson_basic(self):
        """Gibt korrektes Dictionary für Bernoulli-Test zurück."""
        result = HypothesisTesting.neyman_pearson_lemma_demo(p0=0.3, p1=0.7, n=20)
        assert "H0" in result and "H1" in result
        assert result["n"] == 20
        assert "Σxᵢ > c" in result["rejection_region"]

    def test_neyman_pearson_lower_rejection(self):
        """p₁ < p₀: Ablehnbereich ist Σxᵢ < c."""
        result = HypothesisTesting.neyman_pearson_lemma_demo(p0=0.7, p1=0.3, n=20)
        assert "Σxᵢ < c" in result["rejection_region"]

    def test_neyman_pearson_invalid_p0_raises(self):
        """Ungültiges p₀ löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.neyman_pearson_lemma_demo(p0=0.0, p1=0.5, n=10)

    def test_neyman_pearson_equal_p_raises(self):
        """p₀ = p₁ löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.neyman_pearson_lemma_demo(p0=0.5, p1=0.5, n=10)

    def test_neyman_pearson_invalid_n_raises(self):
        """n ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.neyman_pearson_lemma_demo(p0=0.3, p1=0.7, n=0)


class TestUMPTest:
    """Tests für ump_test_normal_mean."""

    def test_ump_reject_when_mean_large(self):
        """Ablehnen von H₀ wenn x̄ deutlich > μ₀."""
        # Stichprobe mit μ=5, testen H₀: μ=0 (sollte sicher abgelehnt werden)
        samples = make_normal_samples(100, mu=5.0, sigma=1.0)
        result = HypothesisTesting.ump_test_normal_mean(samples, mu0=0.0, sigma=1.0)
        assert result["reject_H0"] is True
        assert result["p_value"] < 0.001

    def test_ump_not_reject_under_h0(self):
        """H₀ wird nicht abgelehnt wenn x̄ ≈ μ₀."""
        # Stichprobe unter H₀ (μ=0, testen H₀: μ=0)
        samples = make_normal_samples(100, mu=0.0, sigma=1.0)
        result = HypothesisTesting.ump_test_normal_mean(samples, mu0=0.0, sigma=1.0)
        # Bei μ=μ₀ sollte p-Wert nicht zu klein sein (ca. 50% Wahrscheinlichkeit)
        assert result["p_value"] > 0.001  # Fast sicher nicht abgelehnt

    def test_ump_p_value_in_range(self):
        """p-Wert liegt in [0, 1]."""
        samples = make_normal_samples(50)
        result = HypothesisTesting.ump_test_normal_mean(samples, mu0=0.0, sigma=1.0)
        assert 0.0 <= result["p_value"] <= 1.0

    def test_ump_empty_raises(self):
        """Leere Stichprobe löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.ump_test_normal_mean(np.array([]), mu0=0.0, sigma=1.0)

    def test_ump_invalid_sigma_raises(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.ump_test_normal_mean(
                np.array([1.0, 2.0]), mu0=0.0, sigma=0.0
            )


class TestTypeErrors:
    """Tests für type_i_type_ii_error."""

    def test_basic_output(self):
        """Gibt korrektes Dictionary mit Power zurück."""
        result = HypothesisTesting.type_i_type_ii_error(alpha=0.05, beta=0.2)
        assert math.isclose(result["power"], 0.8, rel_tol=1e-10)
        assert result["alpha"] == 0.05
        assert result["beta"] == 0.2

    def test_invalid_alpha_raises(self):
        """Ungültiges α löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.type_i_type_ii_error(alpha=0.0, beta=0.2)

    def test_invalid_beta_raises(self):
        """Ungültiges β löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.type_i_type_ii_error(alpha=0.05, beta=1.0)


class TestPowerFunction:
    """Tests für power_function."""

    def test_power_at_h0_equals_alpha(self):
        """Gütefunktion bei μ=μ₀ ergibt das Signifikanzniveau α."""
        power = HypothesisTesting.power_function(
            mu=0.0, mu0=0.0, sigma=1.0, n=100, alpha=0.05
        )
        assert math.isclose(power, 0.05, abs_tol=1e-10)

    def test_power_increases_with_effect(self):
        """Gütefunktion nimmt mit wachsender Effektgröße zu."""
        power_small = HypothesisTesting.power_function(0.1, 0.0, 1.0, 100)
        power_large = HypothesisTesting.power_function(1.0, 0.0, 1.0, 100)
        assert power_large > power_small

    def test_power_increases_with_n(self):
        """Gütefunktion nimmt mit wachsendem n zu."""
        power_small_n = HypothesisTesting.power_function(0.5, 0.0, 1.0, 10)
        power_large_n = HypothesisTesting.power_function(0.5, 0.0, 1.0, 1000)
        assert power_large_n > power_small_n

    def test_power_in_range(self):
        """Gütefunktion liegt in [0, 1]."""
        for mu in [-1.0, 0.0, 0.5, 2.0]:
            power = HypothesisTesting.power_function(mu, 0.0, 1.0, 50)
            assert 0.0 <= power <= 1.0

    def test_power_invalid_sigma_raises(self):
        """σ ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.power_function(1.0, 0.0, 0.0, 10)


class TestSPRT:
    """Tests für sequential_probability_ratio_test."""

    def test_sprt_decides_h1_for_p1_data(self):
        """SPRT entscheidet für H₁ bei Daten aus p=p₁."""
        rng = np.random.default_rng(42)
        # Daten aus Ber(0.7), testen H₀: p=0.3 vs H₁: p=0.7
        samples = rng.binomial(1, 0.7, 200).astype(float)
        result = HypothesisTesting.sequential_probability_ratio_test(
            samples, p0=0.3, p1=0.7, alpha=0.05, beta=0.05
        )
        assert "H₁" in result["decision"]

    def test_sprt_decides_h0_for_p0_data(self):
        """SPRT entscheidet für H₀ bei Daten aus p=p₀."""
        rng = np.random.default_rng(7)
        samples = rng.binomial(1, 0.3, 200).astype(float)
        result = HypothesisTesting.sequential_probability_ratio_test(
            samples, p0=0.3, p1=0.7, alpha=0.05, beta=0.05
        )
        assert "H₀" in result["decision"]

    def test_sprt_invalid_p0_raises(self):
        """Ungültiges p₀ löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.sequential_probability_ratio_test(
                np.array([1.0, 0.0]), p0=0.0, p1=0.5
            )

    def test_sprt_invalid_alpha_raises(self):
        """Ungültiges α löst ValueError aus."""
        with pytest.raises(ValueError):
            HypothesisTesting.sequential_probability_ratio_test(
                np.array([1.0, 0.0]), p0=0.3, p1=0.7, alpha=0.0
            )

    def test_sprt_no_decision_short_sequence(self):
        """Kurze Sequenz kann ohne Entscheidung enden."""
        # Nur 1 Beobachtung reicht oft nicht
        result = HypothesisTesting.sequential_probability_ratio_test(
            np.array([1.0]), p0=0.3, p1=0.7, alpha=0.001, beta=0.001
        )
        assert "decision" in result

    def test_sprt_boundaries_correct(self):
        """Wald'sche Grenzen werden korrekt berechnet."""
        result = HypothesisTesting.sequential_probability_ratio_test(
            np.array([1.0, 0.0]), p0=0.3, p1=0.7, alpha=0.05, beta=0.1
        )
        expected_A = (1 - 0.1) / 0.05
        expected_B = 0.1 / (1 - 0.05)
        assert math.isclose(result["A"], expected_A, rel_tol=1e-10)
        assert math.isclose(result["B"], expected_B, rel_tol=1e-10)


# =============================================================================
# 7. Tests für DecisionTheory
# =============================================================================

class TestLossFunctions:
    """Tests für Verlustfunktionen."""

    def test_squared_error_zero_for_correct_estimate(self):
        """Verlust = 0 wenn Schätzung = wahrer Wert."""
        assert DecisionTheory.loss_function_squared_error(3.0, 3.0) == 0.0

    def test_squared_error_positive(self):
        """Quadratischer Verlust ist immer ≥ 0."""
        assert DecisionTheory.loss_function_squared_error(5.0, 2.0) == 9.0

    def test_squared_error_symmetric(self):
        """Quadratischer Verlust ist symmetrisch: L(a,b) = L(b,a)."""
        l1 = DecisionTheory.loss_function_squared_error(3.0, 5.0)
        l2 = DecisionTheory.loss_function_squared_error(5.0, 3.0)
        assert l1 == l2

    def test_absolute_error_zero(self):
        """Absoluter Verlust = 0 wenn Schätzung korrekt."""
        assert DecisionTheory.loss_function_absolute(4.0, 4.0) == 0.0

    def test_absolute_error_positive(self):
        """Absoluter Verlust ist immer ≥ 0."""
        assert math.isclose(DecisionTheory.loss_function_absolute(1.0, 4.0), 3.0)

    def test_absolute_error_symmetric(self):
        """Absoluter Verlust ist symmetrisch."""
        l1 = DecisionTheory.loss_function_absolute(2.0, 7.0)
        l2 = DecisionTheory.loss_function_absolute(7.0, 2.0)
        assert l1 == l2


class TestRiskFunction:
    """Tests für risk_function."""

    def test_risk_of_mean_estimator(self):
        """Risiko von x̄ ≈ σ²/n für N(θ,1) mit σ=1."""
        # Für N(θ,1) und n=50: R(θ, x̄) = 1/50 = 0.02
        risk = DecisionTheory.risk_function(np.mean, theta=2.0, n=50, n_sim=2000)
        assert abs(risk - 1.0 / 50.0) < 0.005

    def test_risk_constant_estimator_higher_than_mean(self):
        """Konstanter Schätzer hat höheres Risiko als x̄ für entferntes θ."""
        # Konstanter Schätzer δ(x) = 0, wahres θ = 3
        const_risk = DecisionTheory.risk_function(
            lambda x: 0.0, theta=3.0, n=50, n_sim=1000
        )
        mean_risk = DecisionTheory.risk_function(np.mean, theta=3.0, n=50, n_sim=1000)
        assert const_risk > mean_risk

    def test_risk_invalid_n_raises(self):
        """n ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            DecisionTheory.risk_function(np.mean, theta=0.0, n=0)

    def test_risk_non_negative(self):
        """Risikofunktion ist nicht-negativ."""
        risk = DecisionTheory.risk_function(np.mean, theta=1.0, n=20)
        assert risk >= 0.0


class TestBayesRisk:
    """Tests für bayes_risk."""

    def test_bayes_risk_returns_float(self):
        """Bayes-Risiko gibt einen Float zurück."""
        theta_grid = np.linspace(0, 1, 10)
        prior = np.ones(10)  # Gleichverteilung
        observations = np.array([0.5, 0.6, 0.4])

        def posterior_fn(theta, obs):
            return float(np.mean(obs))  # x̄ als Schätzer

        risk = DecisionTheory.bayes_risk(
            prior, DecisionTheory.loss_function_squared_error, posterior_fn,
            observations, theta_grid
        )
        assert isinstance(risk, float)
        assert risk >= 0.0

    def test_bayes_risk_length_mismatch_raises(self):
        """Längeninkonsistenz löst ValueError aus."""
        with pytest.raises(ValueError):
            DecisionTheory.bayes_risk(
                np.array([0.5, 0.5]),
                DecisionTheory.loss_function_squared_error,
                lambda t, x: 0.0,
                np.array([1.0]),
                np.array([0.0, 0.5, 1.0]),  # Länge 3 ≠ 2
            )


class TestMinimaxDemo:
    """Tests für minimax_estimator_demo."""

    def test_minimax_demo_returns_dict(self):
        """Gibt Dictionary mit allen erwarteten Schlüsseln zurück."""
        samples = make_normal_samples(100, mu=2.0)
        result = DecisionTheory.minimax_estimator_demo(samples, theta_range=(-5.0, 5.0))
        assert "x_bar" in result
        assert "minimax_estimator" in result
        assert "constant_risk" in result
        assert result["n"] == 100

    def test_minimax_x_bar_correct(self):
        """x̄ wird korrekt berechnet."""
        samples = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = DecisionTheory.minimax_estimator_demo(samples, theta_range=(0.0, 6.0))
        assert math.isclose(result["x_bar"], 3.0, rel_tol=1e-10)

    def test_minimax_empty_raises(self):
        """Leere Beobachtungen lösen ValueError aus."""
        with pytest.raises(ValueError):
            DecisionTheory.minimax_estimator_demo(np.array([]), theta_range=(0.0, 1.0))


class TestAdmissibilityDemo:
    """Tests für admissibility_check_demo."""

    def test_admissibility_returns_dict(self):
        """Gibt Dictionary mit Zulässigkeitsinformationen zurück."""
        result = DecisionTheory.admissibility_check_demo()
        assert "definition" in result
        assert "example_admissible" in result
        assert "example_inadmissible" in result
        assert "bayes_connection" in result

    def test_admissibility_mentions_james_stein(self):
        """Erwähnt James-Stein (Stein-Paradoxon) für unzulässige Schätzer."""
        result = DecisionTheory.admissibility_check_demo()
        assert "James-Stein" in result["example_inadmissible"]


# =============================================================================
# Integration / Konsistenz-Tests
# =============================================================================

class TestIntegration:
    """Integrations- und Konsistenz-Tests über mehrere Klassen."""

    def test_mle_crb_consistency_normal(self):
        """MLE-Varianz ist ≥ CRB (Cramér-Rao-Ungleichung gilt)."""
        sigma = 2.0
        n = 100
        n_trials = 500
        # Empirische Varianz der MLE-Schätzer
        mu_hats = [
            MaximumLikelihoodEstimator.mle_normal(
                np.random.default_rng(i).normal(0, sigma, n)
            )[0]
            for i in range(n_trials)
        ]
        empirical_var = float(np.var(mu_hats))
        crb = CramerRaoBound.cramer_rao_normal_mean(n, sigma)
        # CRB ist untere Schranke: empirische Var ≥ CRB
        assert empirical_var >= crb * 0.9  # 10% Toleranz für Stichprobenfehler

    def test_ci_coverage_probability(self):
        """95%-KI hat empirische Überdeckungsrate ≈ 95%."""
        mu_true, sigma = 3.0, 1.5
        n, n_trials = 50, 500
        contained = 0
        for seed in range(n_trials):
            samples = np.random.default_rng(seed).normal(mu_true, sigma, n)
            lo, hi, _ = ConfidenceIntervals.ci_normal_mean(samples, confidence=0.95)
            if lo < mu_true < hi:
                contained += 1
        coverage = contained / n_trials
        # Empirische Überdeckungsrate sollte nahe 95% liegen
        assert abs(coverage - 0.95) < 0.04

    def test_power_function_type_i_error_control(self):
        """Gütefunktion kontrolliert Fehler 1. Art: β(μ₀) = α."""
        alpha = 0.05
        power_at_h0 = HypothesisTesting.power_function(
            mu=0.0, mu0=0.0, sigma=1.0, n=100, alpha=alpha
        )
        assert math.isclose(power_at_h0, alpha, abs_tol=1e-10)
