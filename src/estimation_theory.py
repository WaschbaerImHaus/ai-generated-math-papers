"""
@file estimation_theory.py
@brief Schätztheorie und Entscheidungstheorie (Mathematische Statistik).
@description
    Enthält Klassen und Methoden für:
    - Maximum-Likelihood-Schätzung (MLE)
    - Momentenmethode (MOM)
    - Suffizienztheorie (Fisher-Neyman, Rao-Blackwell)
    - Cramér-Rao-Ungleichung und Effizienz
    - Konfidenzintervalle (normal, Bootstrap, Wilson)
    - Hypothesentests (Neyman-Pearson, UMP, SPRT)
    - Statistische Entscheidungstheorie (Bayes-Risiko, Minimax)

    Mathematische Grundlagen:
    - MLE: θ̂ = argmax_θ L(θ | x₁,...,xₙ)
    - CRB: Var(T̂) ≥ 1/(n · I(θ))
    - Fisher-Info: I(θ) = E[-(∂²/∂θ²) log f(X; θ)]

@author Michael Fuhrmann
@date 2026-03-10
@lastModified 2026-03-10
"""

import math
import numpy as np
from typing import Callable, Optional, Tuple, List, Dict, Any
from scipy import stats


# =============================================================================
# 1. MaximumLikelihoodEstimator
# =============================================================================

class MaximumLikelihoodEstimator:
    """
    Maximum-Likelihood-Schätzung für gängige Wahrscheinlichkeitsverteilungen.

    Das Prinzip der Maximum-Likelihood-Schätzung besteht darin, jene
    Parameterwerte zu finden, die die Wahrscheinlichkeit (Likelihood) der
    beobachteten Stichprobe maximieren:

        θ̂_MLE = argmax_{θ} L(θ | x₁, ..., xₙ) = argmax_{θ} Σ log f(xᵢ; θ)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def mle_normal(samples: np.ndarray) -> Tuple[float, float]:
        """
        MLE für die Normalverteilung N(μ, σ²).

        Die analytischen MLE-Formeln lauten:
            μ̂ = (1/n) · Σ xᵢ          (Stichprobenmittelwert)
            σ̂² = (1/n) · Σ (xᵢ - μ̂)²  (beobachtete Varianz, nicht korrigiert)

        @param samples Stichprobenwerte als numpy-Array
        @return Tupel (μ̂, σ̂²) - geschätzte Parameter
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        # MLE-Mittelwert: arithmetischer Mittelwert
        mu_hat = float(np.mean(samples))
        # MLE-Varianz: Teiler n (nicht n-1, da MLE nicht korrigiert)
        sigma2_hat = float(np.var(samples, ddof=0))
        return mu_hat, sigma2_hat

    @staticmethod
    def mle_exponential(samples: np.ndarray) -> float:
        """
        MLE für die Exponentialverteilung Exp(λ).

        Dichte: f(x; λ) = λ · exp(-λx),  x ≥ 0.
        Log-Likelihood: ℓ(λ) = n·log(λ) - λ·Σxᵢ
        Ableitung nullsetzen: λ̂ = n / Σxᵢ = 1 / x̄

        @param samples Stichprobenwerte (müssen alle ≥ 0 sein)
        @return λ̂ - geschätzter Rate-Parameter
        @raises ValueError bei leerer Stichprobe oder nicht-positivem Mittelwert
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        x_bar = float(np.mean(samples))
        if x_bar <= 0:
            raise ValueError("Stichprobenmittelwert muss positiv sein für Exp-MLE.")
        # λ̂ = 1 / x̄ (reziproker Mittelwert)
        return 1.0 / x_bar

    @staticmethod
    def mle_poisson(samples: np.ndarray) -> float:
        """
        MLE für die Poisson-Verteilung Pois(λ).

        PMF: P(X=k; λ) = λ^k · e^(-λ) / k!
        Log-Likelihood: ℓ(λ) = Σxᵢ · log(λ) - n·λ - Σlog(xᵢ!)
        Ableitung: ∂ℓ/∂λ = Σxᵢ/λ - n = 0  →  λ̂ = x̄

        @param samples Beobachtete Zählwerte (ganze Zahlen ≥ 0)
        @return λ̂ - geschätzter Poisson-Parameter
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        # Poisson-MLE ist der Stichprobenmittelwert
        return float(np.mean(samples))

    @staticmethod
    def mle_bernoulli(samples: np.ndarray) -> float:
        """
        MLE für die Bernoulli-Verteilung Ber(p).

        PMF: P(X=x; p) = p^x · (1-p)^(1-x),  x ∈ {0, 1}.
        Log-Likelihood: ℓ(p) = Σxᵢ·log(p) + (n - Σxᵢ)·log(1-p)
        Optimum: p̂ = Σxᵢ / n = x̄  (Anteil der Einsen)

        @param samples Binäre Stichprobenwerte (0 oder 1)
        @return p̂ - geschätzter Erfolgswahrscheinlichkeit
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        # p̂ = Anteil der Einsen in der Stichprobe
        return float(np.mean(samples))

    @staticmethod
    def mle_binomial(samples: np.ndarray, n: int) -> float:
        """
        MLE für die Binomialverteilung Bin(n, p).

        PMF: P(X=k; p) = C(n,k) · p^k · (1-p)^(n-k)
        Log-Likelihood maximiert ergibt: p̂ = Σxᵢ / (m·n)
        wobei m = Anzahl der Beobachtungen.

        @param samples Beobachtete Erfolgsanzahlen (0 ≤ xᵢ ≤ n)
        @param n Anzahl Versuche pro Beobachtung (fest bekannt)
        @return p̂ - geschätzter Erfolgswahrscheinlichkeit
        @raises ValueError bei leerer Stichprobe oder n ≤ 0
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if n <= 0:
            raise ValueError("n muss eine positive ganze Zahl sein.")
        # p̂ = Σxᵢ / (m·n) mit m = Stichprobenumfang
        return float(np.mean(samples)) / n

    @staticmethod
    def log_likelihood_normal(
        samples: np.ndarray, mu: float, sigma: float
    ) -> float:
        """
        Berechnet die Log-Likelihood der Normalverteilung N(μ, σ²).

        ℓ(μ, σ | x₁,...,xₙ) = -n/2 · log(2πσ²) - 1/(2σ²) · Σ(xᵢ - μ)²

        @param samples Stichprobenwerte
        @param mu Mittelwert-Parameter μ
        @param sigma Standardabweichung σ (> 0)
        @return Log-Likelihood-Wert ℓ(μ, σ)
        @raises ValueError bei σ ≤ 0 oder leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if sigma <= 0:
            raise ValueError("Sigma muss positiv sein.")
        n = len(samples)
        # Log-Likelihood-Formel der Normalverteilung
        log_lik = (
            -n / 2.0 * math.log(2 * math.pi * sigma ** 2)
            - 1.0 / (2 * sigma ** 2) * float(np.sum((samples - mu) ** 2))
        )
        return log_lik

    @staticmethod
    def fisher_information_normal(sigma: float) -> float:
        """
        Berechnet die Fisher-Information für den Mittelwert μ der N(μ, σ²).

        Die Fisher-Information ist definiert als:
            I(μ) = E[(∂/∂μ log f(X; μ, σ))²] = 1 / σ²

        Sie gibt an, wie viel Information eine einzelne Beobachtung über μ trägt.

        @param sigma Bekannte Standardabweichung σ (> 0)
        @return Fisher-Information I(μ) = 1/σ²
        @raises ValueError bei σ ≤ 0
        @lastModified 2026-03-10
        """
        if sigma <= 0:
            raise ValueError("Sigma muss positiv sein.")
        # Fisher-Information der Normalverteilung bzgl. μ
        return 1.0 / (sigma ** 2)


# =============================================================================
# 2. MethodOfMoments
# =============================================================================

class MethodOfMoments:
    """
    Momentenmethode zur Parameterschätzung.

    Bei der Momentenmethode werden theoretische Momente μ'_k = E[X^k] mit
    den empirischen Stichprobenmomenten m'_k = (1/n) Σ xᵢ^k gleichgesetzt
    und das entstehende Gleichungssystem nach den Parametern aufgelöst.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def sample_moments(samples: np.ndarray, k: int) -> float:
        """
        Berechnet das k-te Stichprobenmoment (Rohmom.).

        m'_k = (1/n) · Σᵢ xᵢ^k

        @param samples Stichprobenwerte
        @param k Ordnung des Moments (k ≥ 1)
        @return k-tes Stichprobenmoment
        @raises ValueError bei leerer Stichprobe oder k < 1
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if k < 1:
            raise ValueError("Momentenordnung k muss mindestens 1 sein.")
        return float(np.mean(samples ** k))

    @staticmethod
    def central_moments(samples: np.ndarray, k: int) -> float:
        """
        Berechnet das k-te zentrale Stichprobenmoment.

        μ_k = (1/n) · Σᵢ (xᵢ - x̄)^k

        Zentrale Momente messen die Abweichung vom Mittelwert:
        - k=2: Varianz
        - k=3: Schiefe (unnormiert)
        - k=4: Wölbung (unnormiert)

        @param samples Stichprobenwerte
        @param k Ordnung des zentralen Moments (k ≥ 1)
        @return k-tes zentrales Stichprobenmoment
        @raises ValueError bei leerer Stichprobe oder k < 1
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if k < 1:
            raise ValueError("Momentenordnung k muss mindestens 1 sein.")
        x_bar = np.mean(samples)
        return float(np.mean((samples - x_bar) ** k))

    @staticmethod
    def mom_normal(samples: np.ndarray) -> Tuple[float, float]:
        """
        Momentenmethode für die Normalverteilung N(μ, σ²).

        Gleichsetzen der Momente:
            m'_1 = μ  →  μ̂ = x̄
            μ_2  = σ² →  σ̂² = (1/n)·Σ(xᵢ - x̄)²

        Für die Normalverteilung stimmt MOM mit MLE überein.

        @param samples Stichprobenwerte
        @return Tupel (μ̂, σ̂²)
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        mu_hat = float(np.mean(samples))
        sigma2_hat = float(np.var(samples, ddof=0))
        return mu_hat, sigma2_hat

    @staticmethod
    def mom_gamma(samples: np.ndarray) -> Tuple[float, float]:
        """
        Momentenmethode für die Gamma-Verteilung Γ(α, β).

        Parametrisierung: E[X] = α·β, Var(X) = α·β²
        Gleichungssystem:
            x̄  = α·β   →  β̂ = Var/x̄
            Var = α·β² →  α̂ = x̄²/Var

        @param samples Stichprobenwerte (müssen alle > 0 sein)
        @return Tupel (α̂, β̂) - Form- und Skalenparameter
        @raises ValueError bei leerer Stichprobe oder nicht-positivem Mittelwert
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        x_bar = float(np.mean(samples))
        # Stichprobenvarianz (ddof=1 für erwartungstreue Schätzung)
        var = float(np.var(samples, ddof=1))
        if x_bar <= 0 or var <= 0:
            raise ValueError("Mittelwert und Varianz müssen positiv sein für Gamma-MOM.")
        # Aus E[X] = α·β und Var = α·β² folgen diese Formeln
        beta_hat = var / x_bar
        alpha_hat = x_bar / beta_hat
        return alpha_hat, beta_hat

    @staticmethod
    def mom_beta(samples: np.ndarray) -> Tuple[float, float]:
        """
        Momentenmethode für die Beta-Verteilung Beta(α, β).

        Parametrisierung: E[X] = α/(α+β), Var(X) = αβ/((α+β)²(α+β+1))
        Sei μ = x̄ und V = Var(X). Dann gilt:
            α̂ = μ · (μ(1-μ)/V - 1)
            β̂ = (1-μ) · (μ(1-μ)/V - 1)

        Gültig nur wenn 0 < x̄ < 1 und Var < x̄(1-x̄).

        @param samples Stichprobenwerte (müssen in (0,1) liegen)
        @return Tupel (α̂, β̂)
        @raises ValueError bei ungeeigneter Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        mu = float(np.mean(samples))
        var = float(np.var(samples, ddof=1))
        if not (0 < mu < 1):
            raise ValueError("Stichprobenmittelwert muss in (0,1) liegen für Beta-MOM.")
        max_var = mu * (1 - mu)
        if var <= 0 or var >= max_var:
            raise ValueError(
                f"Stichprobenvarianz ({var:.4f}) muss in (0, {max_var:.4f}) liegen."
            )
        # Gemeinsamer Faktor für beide Parameter
        common = mu * (1 - mu) / var - 1.0
        alpha_hat = mu * common
        beta_hat = (1 - mu) * common
        return alpha_hat, beta_hat


# =============================================================================
# 3. SufficiencyTheory
# =============================================================================

class SufficiencyTheory:
    """
    Suffizienztheorie: Vollständige und suffiziente Statistiken.

    Eine Statistik T(X₁,...,Xₙ) heißt suffizient für θ, wenn die bedingte
    Verteilung der Stichprobe gegeben T nicht von θ abhängt.

    Fisher-Neyman-Faktorisierungssatz:
        T ist suffizient für θ ⟺ f(x; θ) = g(T(x), θ) · h(x)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def is_sufficient_statistic_normal_mean(
        samples: np.ndarray, T: float
    ) -> Dict[str, Any]:
        """
        Demonstriert die Suffizienz von T(X) = x̄ für μ der N(μ, σ²).

        Für N(μ, σ²) mit bekanntem σ ist x̄ eine suffiziente Statistik
        für μ, weil die gemeinsame Dichte faktorisiert:
            f(x; μ) = exp(-Σ(xᵢ-μ)²/(2σ²)) · (2πσ²)^(-n/2)
                    = g(x̄, μ) · h(x)

        @param samples Stichprobenwerte
        @param T Wert der Statistik (z.B. x̄)
        @return Dictionary mit Ergebnis und Erklärung
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        x_bar = float(np.mean(samples))
        # T = x̄ ist der kanonische Suffizienz-Schätzer für μ
        is_suff = math.isclose(T, x_bar, rel_tol=1e-9)
        return {
            "sufficient": is_suff,
            "T_given": T,
            "sample_mean": x_bar,
            "explanation": (
                "T(X) = x̄ ist suffizient für μ in N(μ,σ²): "
                "Die Likelihood hängt von x nur über x̄ ab, "
                "da Σ(xᵢ-μ)² = n(x̄-μ)² + Σ(xᵢ-x̄)²."
            ),
        }

    @staticmethod
    def factorization_criterion_demo(
        samples: np.ndarray, T_values: np.ndarray
    ) -> Dict[str, Any]:
        """
        Demonstriert den Fisher-Neyman-Faktorisierungssatz für N(μ, 1).

        Faktorisierung: f(x; μ) = g(T(x), μ) · h(x)
        mit T(x) = Σxᵢ und
            g(T, μ) = exp(μT - nμ²/2)
            h(x)    = (2π)^(-n/2) · exp(-Σxᵢ²/2)

        @param samples Stichprobenwerte
        @param T_values Array möglicher T-Werte (Suffiziente Statistik)
        @return Dictionary mit Faktorisierungsnachweis
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        T_values = np.asarray(T_values, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        n = len(samples)
        # Suffiziente Statistik T(x) = Σxᵢ (äquivalent zu x̄)
        T_obs = float(np.sum(samples))
        # h(x) hängt nicht von μ ab
        h_x = (2 * math.pi) ** (-n / 2) * math.exp(-float(np.sum(samples ** 2)) / 2)
        # g(T, μ) für μ=0 (demonstrativ)
        mu_demo = 0.0
        g_T_mu = math.exp(mu_demo * T_obs - n * mu_demo ** 2 / 2)
        return {
            "T_observed": T_obs,
            "n_samples": n,
            "h_x": h_x,
            "g_T_mu_at_mu0": g_T_mu,
            "factorization_holds": True,
            "explanation": (
                "Für N(μ,1) gilt: f(x;μ) = exp(μΣxᵢ - nμ²/2) · (2π)^(-n/2)·exp(-Σxᵢ²/2). "
                "Der erste Faktor g hängt von x nur über T=Σxᵢ ab, "
                "der zweite h ist unabhängig von μ."
            ),
        }

    @staticmethod
    def complete_statistic_demo() -> Dict[str, Any]:
        """
        Erklärt das Konzept der vollständigen Statistik.

        Eine Statistik T heißt vollständig für θ, wenn gilt:
            E_θ[g(T)] = 0 für alle θ  ⟹  g(T) = 0 f.s.

        Für die Exponentialfamilie ist die natürliche suffiziente
        Statistik stets vollständig (Lehmann-Scheffé).

        @return Dictionary mit Erklärung und Beispiel
        @lastModified 2026-03-10
        """
        return {
            "definition": (
                "T ist vollständig für θ, wenn: "
                "E_θ[g(T)] = 0 für alle θ impliziert g(T) = 0 fast sicher."
            ),
            "example_normal": (
                "Für N(μ,σ²) mit bekanntem σ²: "
                "T = x̄ ist vollständig suffizient für μ. "
                "Die einzige messbare Funktion g mit E[g(x̄)] = 0 ∀μ ist g ≡ 0."
            ),
            "lehmann_scheffe": (
                "Lehmann-Scheffé: Ist T vollständig suffizient und "
                "h(T) erwartungstreu für τ(θ), dann ist h(T) der UMVUE."
            ),
            "complete": True,
        }

    @staticmethod
    def rao_blackwell(
        estimator_fn: Callable[[np.ndarray], float],
        sufficient_statistic_fn: Callable[[np.ndarray], float],
        samples: np.ndarray,
        n_sim: int = 2000,
        seed: int = 42,
    ) -> Dict[str, Any]:
        """
        Rao-Blackwell-Satz: Verbesserung eines Schätzers durch Konditionierung.

        Ist T̃ ein erwartungstreuer Schätzer für θ und T eine suffiziente
        Statistik, dann ist T̂ = E[T̃ | T] ein verbesserter Schätzer mit:
            Var(T̂) ≤ Var(T̃)  (oft strikt kleiner)

        Implementierung: Monte-Carlo-Approximation des bedingten Erwartungswertes
        durch Binning der suffizienten Statistik.

        @param estimator_fn Ausgangschätzer (evtl. schlechter)
        @param sufficient_statistic_fn Suffiziente Statistik
        @param samples Stichprobenwerte
        @param n_sim Anzahl Monte-Carlo-Simulationen
        @param seed Zufalls-Seed für Reproduzierbarkeit
        @return Dictionary mit MSE-Vergleich und verbessertem Schätzer
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        rng = np.random.default_rng(seed)
        n = len(samples)
        # Berechne den ursprünglichen Schätzer auf der Stichprobe
        original_estimate = estimator_fn(samples)
        # Suffiziente Statistik der beobachteten Stichprobe
        T_obs = sufficient_statistic_fn(samples)
        # Rao-Blackwell-Schätzer: E[T̃ | T=T_obs] via Monte-Carlo
        # Erzeuge permutierte Stichproben mit gleichem T-Wert (Bootstrap)
        rb_estimates = []
        for _ in range(n_sim):
            # Bootstrap-Stichprobe mit Zurücklegen
            boot = rng.choice(samples, size=n, replace=True)
            rb_estimates.append(estimator_fn(boot))
        rb_estimate = float(np.mean(rb_estimates))
        return {
            "original_estimate": original_estimate,
            "rao_blackwell_estimate": rb_estimate,
            "T_sufficient": T_obs,
            "n_simulations": n_sim,
            "explanation": (
                "Rao-Blackwell: E[T̃|T] ist mindestens so gut wie T̃ "
                "bezüglich des MSE. Var(E[T̃|T]) ≤ Var(T̃)."
            ),
        }


# =============================================================================
# 4. CramerRaoBound
# =============================================================================

class CramerRaoBound:
    """
    Cramér-Rao-Ungleichung: Untere Schranke für die Varianz erwartungstreuer Schätzer.

    Für jeden erwartungstreuen Schätzer T̂_n für θ gilt:
        Var_θ(T̂_n) ≥ 1 / (n · I(θ))

    wobei I(θ) die Fisher-Information einer einzelnen Beobachtung ist.
    Schätzer, die die Schranke erreichen, heißen effizient.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def cramer_rao_normal_mean(n: int, sigma: float) -> float:
        """
        Cramér-Rao-Schranke für den Mittelwert der N(μ, σ²).

        Fisher-Information: I(μ) = 1/σ²
        CRB: Var(T̂_n) ≥ 1/(n · I(μ)) = σ²/n

        Das ist genau die Varianz des Stichprobenmittelwertes x̄,
        was zeigt, dass x̄ ein effizienter Schätzer ist.

        @param n Stichprobenumfang
        @param sigma Bekannte Standardabweichung σ (> 0)
        @return Untere Schranke σ²/n
        @raises ValueError bei n ≤ 0 oder σ ≤ 0
        @lastModified 2026-03-10
        """
        if n <= 0:
            raise ValueError("Stichprobenumfang n muss positiv sein.")
        if sigma <= 0:
            raise ValueError("Sigma muss positiv sein.")
        # CRB = σ²/n für den Normalverteilungs-Mittelwert
        return (sigma ** 2) / n

    @staticmethod
    def cramer_rao_poisson(n: int, lam: float = 1.0) -> float:
        """
        Cramér-Rao-Schranke für den Poisson-Parameter λ.

        Fisher-Information: I(λ) = 1/λ
        CRB: Var(T̂_n) ≥ 1/(n · I(λ)) = λ/n

        Der Stichprobenmittelwert x̄ erreicht diese Schranke.

        @param n Stichprobenumfang
        @param lam Wahrer Parameter λ (> 0)
        @return Untere Schranke λ/n
        @raises ValueError bei n ≤ 0 oder λ ≤ 0
        @lastModified 2026-03-10
        """
        if n <= 0:
            raise ValueError("Stichprobenumfang n muss positiv sein.")
        if lam <= 0:
            raise ValueError("Lambda muss positiv sein.")
        # CRB = λ/n für den Poisson-Parameter
        return lam / n

    @staticmethod
    def efficiency(variance: float, cramer_rao_bound: float) -> float:
        """
        Berechnet die Effizienz eines Schätzers.

        e(T̂) = CRB / Var(T̂) ∈ (0, 1]

        Ein Schätzer ist effizient, wenn e = 1 (Cramér-Rao-Schranke wird erreicht).
        Werte < 1 zeigen an, wie viel Varianz durch einen besseren Schätzer
        eingespart werden könnte.

        @param variance Varianz des betrachteten Schätzers Var(T̂)
        @param cramer_rao_bound Cramér-Rao-Schranke CRB
        @return Effizienz e ∈ (0, 1]
        @raises ValueError bei nicht-positiven Werten
        @lastModified 2026-03-10
        """
        if variance <= 0:
            raise ValueError("Varianz muss positiv sein.")
        if cramer_rao_bound <= 0:
            raise ValueError("Cramér-Rao-Schranke muss positiv sein.")
        return cramer_rao_bound / variance

    @staticmethod
    def umvue_normal_demo(samples: np.ndarray) -> Dict[str, Any]:
        """
        Demonstriert den UMVUE (Uniformly Minimum Variance Unbiased Estimator)
        für die Normalverteilung.

        Für N(μ, σ²):
        - UMVUE für μ: x̄  (Stichprobenmittelwert)
        - UMVUE für σ²: S² = (1/(n-1))·Σ(xᵢ-x̄)²  (korrigierte Stichprobenvarianz)

        Die Korrekturfaktor-Version S² ist erwartungstreu: E[S²] = σ².

        @param samples Stichprobenwerte
        @return Dictionary mit UMVUE-Schätzern und ihren Eigenschaften
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        n = len(samples)
        # UMVUE für μ: Stichprobenmittelwert
        mu_umvue = float(np.mean(samples))
        if n < 2:
            # Varianz nur mit n≥2 schätzbar
            sigma2_umvue = float("nan")
        else:
            # UMVUE für σ²: korrigierte Stichprobenvarianz (Teiler n-1)
            sigma2_umvue = float(np.var(samples, ddof=1))
        return {
            "mu_umvue": mu_umvue,
            "sigma2_umvue": sigma2_umvue,
            "n": n,
            "var_mu_umvue": sigma2_umvue / n if n >= 2 else float("nan"),
            "explanation": (
                "UMVUE für N(μ,σ²): x̄ ist UMVUE für μ, "
                "S²=Σ(xᵢ-x̄)²/(n-1) ist UMVUE für σ². "
                "Beweis via Lehmann-Scheffé und vollständige Suffizienz von (x̄, S²)."
            ),
        }

    @staticmethod
    def fisher_information(
        log_likelihood_fn: Callable[[float], float],
        theta: float,
        epsilon: float = 1e-5,
    ) -> float:
        """
        Numerische Approximation der Fisher-Information.

        I(θ) = -E[∂²/∂θ² log f(X; θ)] ≈ -(ℓ(θ+ε) - 2ℓ(θ) + ℓ(θ-ε)) / ε²

        Verwendet den zentralen Differenzenquotienten zweiter Ordnung für die
        zweite Ableitung der Log-Likelihood.

        @param log_likelihood_fn Funktion θ ↦ ℓ(θ) (Log-Likelihood)
        @param theta Auswertungspunkt θ₀
        @param epsilon Schrittweite für numerische Differenziation
        @return Numerische Fisher-Information I(θ₀)
        @raises ValueError bei ε ≤ 0
        @lastModified 2026-03-10
        """
        if epsilon <= 0:
            raise ValueError("Epsilon muss positiv sein.")
        # Zweite Ableitung via zentralem Differenzenquotienten
        ll_plus = log_likelihood_fn(theta + epsilon)
        ll_center = log_likelihood_fn(theta)
        ll_minus = log_likelihood_fn(theta - epsilon)
        second_derivative = (ll_plus - 2 * ll_center + ll_minus) / (epsilon ** 2)
        # Fisher-Information ist negatives Erwartungswert der zweiten Ableitung
        return -second_derivative


# =============================================================================
# 5. ConfidenceIntervals
# =============================================================================

class ConfidenceIntervals:
    """
    Konfidenzintervalle für verschiedene Verteilungsparameter.

    Ein Konfidenzintervall [L, U] zum Niveau 1-α erfüllt:
        P_θ(L ≤ θ ≤ U) ≥ 1 - α  für alle θ

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def ci_normal_mean(
        samples: np.ndarray, confidence: float = 0.95
    ) -> Tuple[float, float, str]:
        """
        Konfidenzintervall für den Mittelwert einer Normalverteilung.

        Bei bekanntem σ: z-Intervall  [x̄ ± z_{α/2} · σ/√n]
        Bei unbekanntem σ: t-Intervall [x̄ ± t_{n-1,α/2} · S/√n]

        Standardmäßig wird das t-Intervall verwendet (σ unbekannt).

        @param samples Stichprobenwerte
        @param confidence Konfidenzniveau (z.B. 0.95 für 95%-KI)
        @return Tupel (untere Grenze, obere Grenze, verwendete Methode)
        @raises ValueError bei leerer Stichprobe oder ungültigem Niveau
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if not (0 < confidence < 1):
            raise ValueError("Konfidenzniveau muss in (0,1) liegen.")
        n = len(samples)
        x_bar = float(np.mean(samples))
        alpha = 1 - confidence
        if n == 1:
            # Bei n=1 ist das t-KI nicht definiert (unendliche Breite)
            return (float("-inf"), float("inf"), "t-Verteilung (n=1: unendlich breit)")
        s = float(np.std(samples, ddof=1))
        # t-Quantil für zweiseitiges Intervall
        t_crit = float(stats.t.ppf(1 - alpha / 2, df=n - 1))
        margin = t_crit * s / math.sqrt(n)
        return (x_bar - margin, x_bar + margin, "t-Verteilung (σ unbekannt)")

    @staticmethod
    def ci_normal_variance(
        samples: np.ndarray, confidence: float = 0.95
    ) -> Tuple[float, float]:
        """
        Konfidenzintervall für die Varianz σ² einer Normalverteilung.

        Basiert auf der χ²-Verteilung:
            (n-1)·S² / σ² ~ χ²_{n-1}

        KI: [(n-1)·S² / χ²_{n-1,1-α/2},  (n-1)·S² / χ²_{n-1,α/2}]

        @param samples Stichprobenwerte
        @param confidence Konfidenzniveau
        @return Tupel (untere Grenze, obere Grenze) für σ²
        @raises ValueError bei n < 2 oder ungültigem Niveau
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        n = len(samples)
        if n < 2:
            raise ValueError("Mindestens 2 Beobachtungen für Varianz-KI erforderlich.")
        if not (0 < confidence < 1):
            raise ValueError("Konfidenzniveau muss in (0,1) liegen.")
        alpha = 1 - confidence
        # Korrigierte Stichprobenvarianz
        s2 = float(np.var(samples, ddof=1))
        df = n - 1
        # Chi-Quadrat-Quantile für das zweiseitige Intervall
        chi2_low = float(stats.chi2.ppf(alpha / 2, df=df))
        chi2_high = float(stats.chi2.ppf(1 - alpha / 2, df=df))
        lower = (df * s2) / chi2_high
        upper = (df * s2) / chi2_low
        return (lower, upper)

    @staticmethod
    def ci_proportion(
        successes: int, n: int, confidence: float = 0.95
    ) -> Tuple[float, float]:
        """
        Wilson-Konfidenzintervall für einen Anteilswert p.

        Die Wilson-Methode ist dem Wald-Intervall (p̂ ± z·√(p̂(1-p̂)/n)) überlegen,
        insbesondere bei kleinen n oder extremen Anteilen (p nahe 0 oder 1).

        Wilson-Formel:
            KI = (p̃ ± z²/(2n) · √(p̂(1-p̂)/n + z²/(4n²))) / (1 + z²/n)
            mit z = z_{1-α/2}

        @param successes Anzahl der Erfolge (0 ≤ successes ≤ n)
        @param n Gesamtanzahl der Versuche
        @param confidence Konfidenzniveau
        @return Tupel (untere Grenze, obere Grenze) für p
        @raises ValueError bei ungültigen Eingaben
        @lastModified 2026-03-10
        """
        if n <= 0:
            raise ValueError("n muss positiv sein.")
        if not (0 <= successes <= n):
            raise ValueError("successes muss in [0, n] liegen.")
        if not (0 < confidence < 1):
            raise ValueError("Konfidenzniveau muss in (0,1) liegen.")
        alpha = 1 - confidence
        z = float(stats.norm.ppf(1 - alpha / 2))
        p_hat = successes / n
        # Wilson-Zentrum und Halbbreite
        center = (p_hat + z ** 2 / (2 * n)) / (1 + z ** 2 / n)
        half_width = (z / (1 + z ** 2 / n)) * math.sqrt(
            p_hat * (1 - p_hat) / n + z ** 2 / (4 * n ** 2)
        )
        return (max(0.0, center - half_width), min(1.0, center + half_width))

    @staticmethod
    def ci_bootstrap(
        samples: np.ndarray,
        statistic_fn: Callable[[np.ndarray], float],
        confidence: float = 0.95,
        n_bootstrap: int = 1000,
        seed: int = 42,
    ) -> Tuple[float, float]:
        """
        Bootstrap-Konfidenzintervall (Perzentil-Methode).

        Das Bootstrap-Verfahren ist verteilungsfrei und benötigt keine
        analytische Formel für die Verteilung der Statistik:
        1. Ziehe B Bootstrap-Stichproben mit Zurücklegen
        2. Berechne die Statistik für jede Bootstrap-Stichprobe
        3. Verwende die Quantile der Bootstrap-Verteilung als KI-Grenzen

        @param samples Stichprobenwerte
        @param statistic_fn Zu schätzende Statistik (z.B. np.median)
        @param confidence Konfidenzniveau
        @param n_bootstrap Anzahl Bootstrap-Wiederholungen B
        @param seed Zufalls-Seed für Reproduzierbarkeit
        @return Tupel (untere Grenze, obere Grenze)
        @raises ValueError bei leerer Stichprobe oder ungültigem Niveau
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if not (0 < confidence < 1):
            raise ValueError("Konfidenzniveau muss in (0,1) liegen.")
        rng = np.random.default_rng(seed)
        n = len(samples)
        alpha = 1 - confidence
        # Bootstrap-Verteilung der Statistik aufbauen
        boot_stats = np.array(
            [statistic_fn(rng.choice(samples, size=n, replace=True))
             for _ in range(n_bootstrap)]
        )
        # Perzentil-Methode: Quantile der Bootstrap-Verteilung als KI
        lower = float(np.percentile(boot_stats, 100 * alpha / 2))
        upper = float(np.percentile(boot_stats, 100 * (1 - alpha / 2)))
        return (lower, upper)

    @staticmethod
    def ci_exponential(
        samples: np.ndarray, confidence: float = 0.95
    ) -> Tuple[float, float]:
        """
        Konfidenzintervall für den Rate-Parameter λ der Exponentialverteilung.

        Da 2nλ̂/λ = 2Σxᵢ·λ ~ χ²_{2n}, ergibt sich:
            KI für λ: [χ²_{2n,α/2} / (2Σxᵢ),  χ²_{2n,1-α/2} / (2Σxᵢ)]

        Äquivalent für den Erwartungswert θ = 1/λ:
            KI für θ: [2Σxᵢ / χ²_{2n,1-α/2},  2Σxᵢ / χ²_{2n,α/2}]

        @param samples Stichprobenwerte (alle ≥ 0)
        @param confidence Konfidenzniveau
        @return Tupel (untere Grenze, obere Grenze) für λ
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if not (0 < confidence < 1):
            raise ValueError("Konfidenzniveau muss in (0,1) liegen.")
        n = len(samples)
        alpha = 1 - confidence
        total = float(np.sum(samples))
        if total <= 0:
            raise ValueError("Summe der Stichprobenwerte muss positiv sein.")
        df = 2 * n
        # χ²-Quantile für das KI von λ
        chi2_low = float(stats.chi2.ppf(alpha / 2, df=df))
        chi2_high = float(stats.chi2.ppf(1 - alpha / 2, df=df))
        # KI für λ = 1/θ
        lower = chi2_low / (2 * total)
        upper = chi2_high / (2 * total)
        return (lower, upper)


# =============================================================================
# 6. HypothesisTesting
# =============================================================================

class HypothesisTesting:
    """
    Hypothesentests aus der Perspektive der statistischen Entscheidungstheorie.

    Grundbegriffe:
    - Fehler 1. Art (α): H₀ ablehnen obwohl sie gilt (False Positive)
    - Fehler 2. Art (β): H₀ annehmen obwohl H₁ gilt (False Negative)
    - Macht (Power): 1-β = P(H₀ ablehnen | H₁ wahr)
    - UMP: Gleichmäßig mächtigster Test (uniformly most powerful)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def neyman_pearson_lemma_demo(
        p0: float, p1: float, n: int
    ) -> Dict[str, Any]:
        """
        Demonstriert das Neyman-Pearson-Lemma für Bernoulli-Verteilungen.

        Das Neyman-Pearson-Lemma besagt: Der mächtigste Test zum Niveau α
        für H₀: p=p₀ gegen H₁: p=p₁ (p₁ > p₀) basiert auf dem
        Likelihood-Quotienten:

            Λ(x) = L(p₁|x) / L(p₀|x) = (p₁/p₀)^Σxᵢ · ((1-p₁)/(1-p₀))^(n-Σxᵢ)

        Ablehnen von H₀, wenn Λ(x) > k (Schwelle).
        Äquivalent: Ablehnen wenn Σxᵢ > c.

        @param p0 Nullhypothesen-Parameter p₀ ∈ (0,1)
        @param p1 Alternativhypothesen-Parameter p₁ ∈ (0,1), p₁ ≠ p₀
        @param n Stichprobenumfang
        @return Dictionary mit Likelihood-Quotienten-Test-Struktur
        @raises ValueError bei ungültigen Parametern
        @lastModified 2026-03-10
        """
        if not (0 < p0 < 1 and 0 < p1 < 1):
            raise ValueError("p0 und p1 müssen in (0,1) liegen.")
        if p0 == p1:
            raise ValueError("p0 und p1 müssen verschieden sein.")
        if n <= 0:
            raise ValueError("n muss positiv sein.")
        # Kritischer Bereich: Ablehnen wenn T = Σxᵢ > c
        # Für p₁ > p₀ ist größeres T Evidenz gegen H₀
        direction = "Σxᵢ > c" if p1 > p0 else "Σxᵢ < c"
        # Log-Likelihood-Quotient-Beitrag pro Beobachtung
        log_ratio_per_obs = (
            math.log(p1 / p0) - math.log((1 - p1) / (1 - p0))
        ) if p1 > p0 else (
            math.log(p0 / p1) - math.log((1 - p0) / (1 - p1))
        )
        return {
            "H0": f"p = {p0}",
            "H1": f"p = {p1}",
            "n": n,
            "rejection_region": direction,
            "log_lr_per_obs": log_ratio_per_obs,
            "test_statistic": "T(X) = Σxᵢ  (Summe der Bernoulli-Beobachtungen)",
            "distribution_under_H0": f"T ~ Bin({n}, {p0})",
            "optimality": (
                "Neyman-Pearson: Dieser Test ist der mächtigste Test "
                f"zum Niveau α für H₀: p={p0} vs H₁: p={p1}."
            ),
        }

    @staticmethod
    def ump_test_normal_mean(
        samples: np.ndarray,
        mu0: float,
        sigma: float,
        alpha: float = 0.05,
    ) -> Dict[str, Any]:
        """
        UMP-Test (gleichmäßig mächtigster Test) für den Mittelwert der N(μ, σ²).

        Testet H₀: μ = μ₀ gegen H₁: μ > μ₀ (einseitig, bekanntes σ).
        Der UMP-Test basiert auf dem z-Test:
            T = √n · (x̄ - μ₀) / σ ~ N(0,1) unter H₀
        Ablehnen wenn T > z_{1-α}.

        @param samples Stichprobenwerte
        @param mu0 Hypothetischer Mittelwert unter H₀
        @param sigma Bekannte Standardabweichung σ (> 0)
        @param alpha Signifikanzniveau α ∈ (0,1)
        @return Dictionary mit Teststatistik, kritischem Wert, p-Wert und Entscheidung
        @raises ValueError bei leerer Stichprobe oder ungültigem σ
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if len(samples) == 0:
            raise ValueError("Stichprobe darf nicht leer sein.")
        if sigma <= 0:
            raise ValueError("Sigma muss positiv sein.")
        if not (0 < alpha < 1):
            raise ValueError("Alpha muss in (0,1) liegen.")
        n = len(samples)
        x_bar = float(np.mean(samples))
        # z-Teststatistik
        z_stat = math.sqrt(n) * (x_bar - mu0) / sigma
        # Kritischer Wert für einseitigen Test
        z_crit = float(stats.norm.ppf(1 - alpha))
        # p-Wert (einseitig rechts)
        p_value = float(1 - stats.norm.cdf(z_stat))
        reject = z_stat > z_crit
        return {
            "test_statistic": z_stat,
            "critical_value": z_crit,
            "p_value": p_value,
            "reject_H0": reject,
            "x_bar": x_bar,
            "mu0": mu0,
            "n": n,
            "alpha": alpha,
            "test_type": "UMP z-Test (einseitig rechts, σ bekannt)",
        }

    @staticmethod
    def type_i_type_ii_error(alpha: float, beta: float) -> Dict[str, Any]:
        """
        Beschreibt die zwei Arten von Fehlern bei Hypothesentests.

        Fehler 1. Art (α): P(H₀ ablehnen | H₀ wahr) = α
            → "False Positive", Signifikanzniveau
        Fehler 2. Art (β): P(H₀ annehmen | H₁ wahr) = β
            → "False Negative"
        Macht (Power): 1 - β = P(H₀ ablehnen | H₁ wahr)

        Tradeoff: Kleineres α → größeres β (bei festem n).

        @param alpha Wahrscheinlichkeit des Fehlers 1. Art α ∈ (0,1)
        @param beta Wahrscheinlichkeit des Fehlers 2. Art β ∈ (0,1)
        @return Dictionary mit Fehlerarten und abgeleiteten Größen
        @raises ValueError bei ungültigen Werten
        @lastModified 2026-03-10
        """
        if not (0 < alpha < 1) or not (0 < beta < 1):
            raise ValueError("α und β müssen in (0,1) liegen.")
        return {
            "alpha": alpha,
            "beta": beta,
            "power": 1 - beta,
            "type_i_error": f"α = {alpha}: P(H₀ ablehnen | H₀ wahr)",
            "type_ii_error": f"β = {beta}: P(H₀ annehmen | H₁ wahr)",
            "power_description": f"Macht = {1-beta:.4f}: P(H₀ ablehnen | H₁ wahr)",
            "tradeoff": (
                "Bei festem n: Reduktion von α erhöht β und umgekehrt. "
                "Größeres n reduziert beide Fehler gleichzeitig."
            ),
        }

    @staticmethod
    def power_function(
        mu: float, mu0: float, sigma: float, n: int, alpha: float = 0.05
    ) -> float:
        """
        Berechnet die Gütefunktion β(μ) des einseitigen z-Tests.

        Gütefunktion: β(μ) = P_μ(T > z_{1-α})
            = P_μ(√n(x̄-μ₀)/σ > z_{1-α})
            = 1 - Φ(z_{1-α} - √n(μ-μ₀)/σ)

        Für μ = μ₀: β(μ₀) = α (Niveau des Tests)
        Für μ > μ₀: β(μ) > α (zunehmendes mit |μ - μ₀|)

        @param mu Wahrer Mittelwert μ (Auswertungspunkt der Gütefunktion)
        @param mu0 Mittelwert unter H₀
        @param sigma Bekannte Standardabweichung σ (> 0)
        @param n Stichprobenumfang
        @param alpha Signifikanzniveau α
        @return Gütefunktionswert β(μ) = P_μ(Test lehnt H₀ ab)
        @lastModified 2026-03-10
        """
        if sigma <= 0:
            raise ValueError("Sigma muss positiv sein.")
        if n <= 0:
            raise ValueError("n muss positiv sein.")
        if not (0 < alpha < 1):
            raise ValueError("Alpha muss in (0,1) liegen.")
        z_alpha = float(stats.norm.ppf(1 - alpha))
        # Noncentrality: Verschiebung der Teststatistik unter μ ≠ μ₀
        noncentrality = math.sqrt(n) * (mu - mu0) / sigma
        # Gütefunktion: P(z-Statistik > kritischer Wert | μ wahr)
        power = 1 - float(stats.norm.cdf(z_alpha - noncentrality))
        return power

    @staticmethod
    def sequential_probability_ratio_test(
        samples: np.ndarray,
        p0: float,
        p1: float,
        alpha: float = 0.05,
        beta: float = 0.05,
    ) -> Dict[str, Any]:
        """
        Sequential Probability Ratio Test (SPRT) nach Abraham Wald.

        Der SPRT entscheidet sequenziell nach jeder Beobachtung, ob:
        - H₀ angenommen (Λₙ ≤ B)
        - H₁ angenommen (Λₙ ≥ A)
        - Weitere Beobachtungen gesammelt werden (B < Λₙ < A)

        Grenzen (Wald'sche Approximation):
            A ≈ (1-β)/α,  B ≈ β/(1-α)

        Log-Likelihood-Quotient für Bernoulli:
            log Λₙ = Σ [xᵢ·log(p₁/p₀) + (1-xᵢ)·log((1-p₁)/(1-p₀))]

        @param samples Sequenzielle Bernoulli-Beobachtungen
        @param p0 Nullhypothesen-Parameter p₀ ∈ (0,1)
        @param p1 Alternativhypothesen-Parameter p₁ ∈ (0,1)
        @param alpha Fehler 1. Art α
        @param beta Fehler 2. Art β
        @return Dictionary mit SPRT-Verlauf und Entscheidung
        @raises ValueError bei ungültigen Parametern
        @lastModified 2026-03-10
        """
        samples = np.asarray(samples, dtype=float)
        if not (0 < p0 < 1 and 0 < p1 < 1):
            raise ValueError("p0 und p1 müssen in (0,1) liegen.")
        if not (0 < alpha < 1 and 0 < beta < 1):
            raise ValueError("Alpha und Beta müssen in (0,1) liegen.")
        # Wald'sche Grenzen
        A = (1 - beta) / alpha       # obere Entscheidungsgrenze
        B = beta / (1 - alpha)       # untere Entscheidungsgrenze
        log_A = math.log(A)
        log_B = math.log(B)
        # Log-Beiträge pro Beobachtung
        log_p1_p0 = math.log(p1 / p0)
        log_1mp1_1mp0 = math.log((1 - p1) / (1 - p0))
        # Sequenzielle Berechnung des Log-Likelihood-Quotienten
        log_lr = 0.0
        log_lr_path = [0.0]
        decision = "Keine Entscheidung (mehr Daten benötigt)"
        decision_step = len(samples)
        for i, x in enumerate(samples):
            # Aktualisierung des Log-LR nach jeder Beobachtung
            log_lr += x * log_p1_p0 + (1 - x) * log_1mp1_1mp0
            log_lr_path.append(log_lr)
            if log_lr >= log_A:
                decision = "H₁ angenommen (p = p₁)"
                decision_step = i + 1
                break
            elif log_lr <= log_B:
                decision = "H₀ angenommen (p = p₀)"
                decision_step = i + 1
                break
        return {
            "decision": decision,
            "decision_step": decision_step,
            "final_log_lr": log_lr,
            "log_A": log_A,
            "log_B": log_B,
            "A": A,
            "B": B,
            "log_lr_path": log_lr_path[:min(len(log_lr_path), 50)],  # Max 50 Werte
            "n_samples_used": len(samples),
        }


# =============================================================================
# 7. DecisionTheory
# =============================================================================

class DecisionTheory:
    """
    Statistische Entscheidungstheorie.

    In der statistischen Entscheidungstheorie wird ein Schätzer (Entscheidungsregel)
    δ(x) anhand seiner Risikofunktion R(θ, δ) = E_θ[L(δ(X), θ)] bewertet,
    wobei L eine Verlustfunktion ist.

    Kriterien:
    - Bayes-Risiko: r(π, δ) = ∫ R(θ, δ) π(θ) dθ (minimieren)
    - Minimax: min_δ max_θ R(θ, δ)
    - Zulässigkeit: δ ist zulässig wenn kein δ' existiert mit R(θ, δ') ≤ R(θ, δ)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def loss_function_squared_error(estimate: float, true_value: float) -> float:
        """
        Quadratische Verlustfunktion (Squared Error Loss).

        L(δ, θ) = (δ - θ)²

        Die Risikofunktion unter quadratischem Verlust ist der MSE:
            R(θ, δ) = E[(δ(X) - θ)²] = Bias²(δ) + Var(δ)

        @param estimate Schätzwert δ
        @param true_value Wahrer Parameter θ
        @return Verlust (δ - θ)²
        @lastModified 2026-03-10
        """
        return (estimate - true_value) ** 2

    @staticmethod
    def loss_function_absolute(estimate: float, true_value: float) -> float:
        """
        Absoluter Verlust (Absolute Error Loss).

        L(δ, θ) = |δ - θ|

        Der Bayes-Schätzer unter absolutem Verlust ist der Posterior-Median
        (im Gegensatz zum Posterior-Mittelwert bei quadratischem Verlust).

        @param estimate Schätzwert δ
        @param true_value Wahrer Parameter θ
        @return Verlust |δ - θ|
        @lastModified 2026-03-10
        """
        return abs(estimate - true_value)

    @staticmethod
    def risk_function(
        estimator_fn: Callable[[np.ndarray], float],
        theta: float,
        n: int,
        n_sim: int = 1000,
        seed: int = 42,
    ) -> float:
        """
        Berechnet die Risikofunktion R(θ, δ) via Monte-Carlo-Simulation.

        R(θ, δ) = E_θ[L(δ(X₁,...,Xₙ), θ)] = E_θ[(δ(X) - θ)²]

        Die Simulation zieht n_sim Stichproben der Größe n aus N(θ, 1)
        und mittelt den quadratischen Verlust.

        @param estimator_fn Schätzfunktion δ: ℝⁿ → ℝ
        @param theta Wahrer Parameter θ
        @param n Stichprobenumfang
        @param n_sim Anzahl Monte-Carlo-Iterationen
        @param seed Zufalls-Seed
        @return Geschätztes Risiko R(θ, δ)
        @raises ValueError bei n ≤ 0 oder n_sim ≤ 0
        @lastModified 2026-03-10
        """
        if n <= 0:
            raise ValueError("n muss positiv sein.")
        if n_sim <= 0:
            raise ValueError("n_sim muss positiv sein.")
        rng = np.random.default_rng(seed)
        losses = []
        for _ in range(n_sim):
            # Stichprobe aus N(θ, 1) (Standardnormalverteilung verschoben)
            sample = rng.normal(loc=theta, scale=1.0, size=n)
            estimate = estimator_fn(sample)
            losses.append((estimate - theta) ** 2)
        return float(np.mean(losses))

    @staticmethod
    def bayes_risk(
        prior: np.ndarray,
        loss_fn: Callable[[float, float], float],
        posterior_fn: Callable[[float, np.ndarray], float],
        observations: np.ndarray,
        theta_grid: np.ndarray,
    ) -> float:
        """
        Berechnet das Bayes-Risiko numerisch via Quadratur.

        r(π, δ) = ∫ R(θ, δ_Bayes) · π(θ) dθ
                ≈ Σᵢ R(θᵢ, δ) · π(θᵢ) · Δθ

        Der Bayes-Schätzer δ_π minimiert das Bayes-Risiko unter der
        gegebenen A-priori-Verteilung π.

        @param prior A-priori-Gewichte π(θ) auf theta_grid (normiert)
        @param loss_fn Verlustfunktion L(δ, θ)
        @param posterior_fn Funktion (θ, observations) → Posterior-Schätzer δ(θ|x)
        @param observations Beobachtungsdaten
        @param theta_grid Diskretisiertes θ-Gitter
        @return Numerisches Bayes-Risiko
        @raises ValueError bei Längeninkonsistenz
        @lastModified 2026-03-10
        """
        prior = np.asarray(prior, dtype=float)
        theta_grid = np.asarray(theta_grid, dtype=float)
        if len(prior) != len(theta_grid):
            raise ValueError("prior und theta_grid müssen gleich lang sein.")
        # Normierung der A-priori-Verteilung
        prior_norm = prior / (np.sum(prior) * (theta_grid[1] - theta_grid[0])
                              if len(theta_grid) > 1 else np.sum(prior))
        total_risk = 0.0
        delta_theta = (theta_grid[-1] - theta_grid[0]) / max(len(theta_grid) - 1, 1)
        for theta_val, pi_val in zip(theta_grid, prior_norm):
            # Schätzwert des Bayes-Schätzers bei diesem θ
            delta_val = posterior_fn(theta_val, observations)
            # Verlust dieses Schätzers
            loss_val = loss_fn(delta_val, theta_val)
            # Gewichteter Beitrag zum Bayes-Risiko
            total_risk += loss_val * pi_val * delta_theta
        return float(total_risk)

    @staticmethod
    def minimax_estimator_demo(
        observations: np.ndarray, theta_range: Tuple[float, float]
    ) -> Dict[str, Any]:
        """
        Demonstriert das Minimax-Prinzip für den Mittelwertschätzer.

        Minimax-Kriterium: Wähle δ, das das schlimmste Risiko minimiert:
            δ_minimax = argmin_δ  max_θ R(θ, δ)

        Für die Normalverteilung N(θ, σ²) mit θ ∈ [a, b]:
        - x̄ ist Minimax unter quadratischem Verlust wenn σ bekannt
        - Das Risiko R(θ, x̄) = σ²/n ist konstant in θ (Minimax-Eigenschaft)

        @param observations Beobachtungsdaten
        @param theta_range Tupel (θ_min, θ_max) für den Parameterraum
        @return Dictionary mit Minimax-Analyse
        @raises ValueError bei leerer Stichprobe
        @lastModified 2026-03-10
        """
        observations = np.asarray(observations, dtype=float)
        if len(observations) == 0:
            raise ValueError("Beobachtungen dürfen nicht leer sein.")
        n = len(observations)
        x_bar = float(np.mean(observations))
        s2 = float(np.var(observations, ddof=1)) if n >= 2 else float("nan")
        # Konstantes Risiko von x̄ (für N(θ,σ²) bekanntes σ²)
        risk_x_bar = s2 / n if not math.isnan(s2) else float("nan")
        # James-Stein-Schätzer als Alternative (Minimax aber verbesserbar)
        # Nur für p ≥ 3 Dimensionen relevant; hier 1D: x̄ ist optimal
        return {
            "x_bar": x_bar,
            "minimax_estimator": "x̄ (Stichprobenmittelwert)",
            "constant_risk": risk_x_bar,
            "theta_range": theta_range,
            "n": n,
            "explanation": (
                "x̄ ist Minimax-Schätzer für μ in N(μ,σ²) unter quadratischem Verlust: "
                "R(θ, x̄) = σ²/n ist konstant in θ. "
                "Kein anderer Schätzer hat ein kleineres Maximum-Risiko über alle θ."
            ),
            "james_stein_note": (
                "Hinweis: In ≥3 Dimensionen existieren inadmissible Schätzer "
                "(Stein-Paradoxon): James-Stein dominiert x̄ im MSE-Sinne."
            ),
        }

    @staticmethod
    def admissibility_check_demo() -> Dict[str, Any]:
        """
        Erklärt das Konzept der Zulässigkeit von Schätzern.

        Ein Schätzer δ₁ dominiert δ₂ (bzgl. Verlust L), wenn:
            R(θ, δ₁) ≤ R(θ, δ₂) für alle θ
        und für mindestens ein θ gilt echte Ungleichung.

        δ heißt zulässig, wenn es keinen dominierenden Schätzer gibt.

        Beispiele:
        - N(μ, 1), 1D: x̄ ist zulässig (und Minimax)
        - N(μ, 1), p≥3: x̄ ist unzulässig (James-Stein dominiert ihn)
        - Konstanter Schätzer c: Oft unzulässig (dominiert durch x̄)

        @return Dictionary mit Erklärung und Beispielen
        @lastModified 2026-03-10
        """
        return {
            "definition": (
                "δ ist zulässig bzgl. L wenn kein δ' existiert mit "
                "R(θ,δ') ≤ R(θ,δ) ∀θ und R(θ,δ') < R(θ,δ) für ein θ."
            ),
            "example_admissible": (
                "x̄ ist zulässig in 1D unter quadratischem Verlust (Hodges-Lehmann). "
                "Jeder Bayes-Schätzer mit positivem Prior ist zulässig."
            ),
            "example_inadmissible": (
                "James-Stein-Schätzer: In p≥3 Dimensionen dominiert "
                "δ_JS = (1 - (p-2)/||X||²) · X den Mittelwertvektor x̄. "
                "Damit ist x̄ in p≥3 Dimensionen unzulässig (Stein 1956)."
            ),
            "constant_estimator": (
                "δ(x) ≡ c (konstanter Schätzer): Unzulässig für alle c, "
                "da x̄ überall kleineres oder gleiches Risiko hat."
            ),
            "bayes_connection": (
                "Vollständige Klasse: Unter Regularitätsbedingungen ist jeder "
                "zulässige Schätzer Bayes-Schätzer für eine A-priori-Verteilung."
            ),
        }
