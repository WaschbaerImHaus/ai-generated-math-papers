# Modul: Schätztheorie (`estimation_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `estimation_theory.py` implementiert die mathematische Statistik mit Schwerpunkt auf Schätz- und Entscheidungstheorie. Es umfasst Maximum-Likelihood-Schätzung (MLE), die Momentenmethode (MOM), Suffizienz und den Rao-Blackwell-Satz, die Cramér-Rao-Ungleichung, Konfidenzintervalle (normalverteilt, Bootstrap, Wilson), Hypothesentests (Neyman-Pearson-Lemma, UMP-Tests, SPRT) sowie statistische Entscheidungstheorie (Bayes-Risiko, Minimax).

---

## Mathematischer Hintergrund

### Maximum-Likelihood-Schätzung (MLE)

Der MLE-Schätzer maximiert die Likelihood-Funktion über den Parameterraum:

$$\hat{\theta}_{\text{MLE}} = \arg\max_\theta L(\theta \mid x_1, \ldots, x_n) = \arg\max_\theta \sum_{i=1}^n \log f(x_i; \theta)$$

**Analytische MLEs:**

| Verteilung | MLE |
|------------|-----|
| Normalverteilung $\mathcal{N}(\mu, \sigma^2)$ | $\hat{\mu} = \bar{x}$, $\hat{\sigma}^2 = \frac{1}{n}\sum(x_i - \bar{x})^2$ |
| Exponentialverteilung $\text{Exp}(\lambda)$ | $\hat{\lambda} = 1/\bar{x}$ |
| Poisson $\text{Poi}(\lambda)$ | $\hat{\lambda} = \bar{x}$ |
| Bernoulli $\text{Ber}(p)$ | $\hat{p} = \bar{x}$ |

### Cramér-Rao-Ungleichung

Die Varianz jedes **erwartungstreuen** Schätzers $\hat{\theta}$ für $\theta$ ist nach unten beschränkt durch die Cramér-Rao-Schranke (CRB):

$$\operatorname{Var}(\hat{\theta}) \geq \frac{1}{n \cdot I(\theta)}$$

**Fisher-Information:**
$$I(\theta) = \mathbb{E}\!\left[-\frac{\partial^2}{\partial \theta^2} \log f(X; \theta)\right] = \mathbb{E}\!\left[\left(\frac{\partial}{\partial \theta} \log f(X; \theta)\right)^2\right]$$

**Effizienz** eines Schätzers: $e(\hat{\theta}) = \frac{\text{CRB}}{\operatorname{Var}(\hat{\theta})} \in (0, 1]$. Ein Schätzer heißt **effizient**, wenn $e = 1$.

### Suffizienz und Rao-Blackwell

**Suffizienzkriterium (Fisher-Neyman):** $T(X)$ ist suffizienter Schätzer für $\theta$, wenn $f(x; \theta) = g(T(x), \theta) \cdot h(x)$ (Faktorisierung).

**Rao-Blackwell-Satz:** Sei $\hat{\theta}$ ein beliebiger Schätzer und $T$ suffizient. Dann ist $\tilde{\theta} = \mathbb{E}[\hat{\theta} \mid T]$ mindestens so gut:
$$\operatorname{Var}(\tilde{\theta}) \leq \operatorname{Var}(\hat{\theta})$$

**UMVUE** (Uniformly Minimum Variance Unbiased Estimator): Bester erwartungstreuer Schätzer – existiert genau dann, wenn ein vollständiger suffizienter Schätzer existiert (Lehmann-Scheffé-Satz).

### Konfidenzintervalle

**Normalverteiltes Mittel** ($\sigma^2$ bekannt): $\bar{x} \pm z_{\alpha/2} \cdot \sigma/\sqrt{n}$

**Normalverteiltes Mittel** ($\sigma^2$ unbekannt): $\bar{x} \pm t_{n-1,\alpha/2} \cdot s/\sqrt{n}$

**Varianz:** $\left[\frac{(n-1)s^2}{\chi^2_{n-1,\alpha/2}},\; \frac{(n-1)s^2}{\chi^2_{n-1,1-\alpha/2}}\right]$

**Wilson-Intervall** für Proportion $p$: löst das Überschreitungsproblem des Wald-Intervalls bei kleinen Stichproben.

**Bootstrap:** Nicht-parametrisches Konfidenzintervall durch Resampling (mit Zurücklegen), $B$ Replikationen.

### Neyman-Pearson-Lemma

Für den Test $H_0: \theta = \theta_0$ gegen $H_1: \theta = \theta_1$ ist der **Likelihood-Quotient-Test** (LRT) mit kritischer Region $\Lambda = L(\theta_1)/L(\theta_0) > c$ der **gleichmäßig stärkste Test** (UMP) zum Niveau $\alpha$:

$$\Lambda(\mathbf{x}) = \frac{L(\theta_1 \mid \mathbf{x})}{L(\theta_0 \mid \mathbf{x})} > c_\alpha$$

### Sequentieller Wahrscheinlichkeitsquotienten-Test (SPRT)

Wald'scher Test: Beobachte sequenziell bis $\Lambda_n$ eine Schranke $B = (1-\beta)/\alpha$ überschreitet ($H_1$ annehmen) oder unter $A = \beta/(1-\alpha)$ fällt ($H_0$ annehmen):

$$\Lambda_n = \prod_{i=1}^n \frac{f(x_i; \theta_1)}{f(x_i; \theta_0)}$$

### Statistische Entscheidungstheorie

**Verlustfunktion:** Bewertet den Schaden $L(\hat{\theta}, \theta)$ einer Entscheidung.

**Risikofunktion:** $R(\hat{\theta}, \theta) = \mathbb{E}_\theta[L(\hat{\theta}(X), \theta)]$

**Bayes-Risiko:** $r(\hat{\theta}) = \int R(\hat{\theta}, \theta) \, \pi(\theta) \, d\theta$ mit Priori-Verteilung $\pi$.

**Minimax-Schätzer:** Minimiert das schlechtmöglichste Risiko: $\min_{\hat{\theta}} \sup_\theta R(\hat{\theta}, \theta)$.

**Zulässigkeit:** $\hat{\theta}$ ist zulässig, wenn kein anderer Schätzer $\hat{\theta}'$ existiert mit $R(\hat{\theta}', \theta) \leq R(\hat{\theta}, \theta)$ für alle $\theta$ (echt kleiner für mindestens ein $\theta$).

---

## Klassen und Methoden

### `MaximumLikelihoodEstimator`

MLE für gängige Verteilungsfamilien.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `mle_normal` | `(samples) → (μ, σ²)` | MLE für Normalverteilung |
| `mle_exponential` | `(samples) → λ` | MLE für Exponentialverteilung |
| `mle_poisson` | `(samples) → λ` | MLE für Poisson-Verteilung |
| `mle_bernoulli` | `(samples) → p` | MLE für Bernoulli-Verteilung |
| `mle_binomial` | `(samples, n) → p` | MLE für Binomialverteilung bei bekanntem $n$ |
| `log_likelihood_normal` | `(samples, mu, sigma) → float` | Log-Likelihood-Wert $\ell(\mu,\sigma \mid x)$ |
| `fisher_information_normal` | `(sigma) → float` | Fisher-Information $I(\mu) = 1/\sigma^2$ |

### `MethodOfMoments`

Momentenmethode – Schätzung durch Gleichsetzen von Stichproben- und theoretischen Momenten.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `sample_moments` | `(samples, k) → float` | $k$-tes Stichprobenmoment $\hat{\mu}_k = \frac{1}{n}\sum x_i^k$ |
| `central_moments` | `(samples, k) → float` | $k$-tes zentrales Moment |
| `mom_normal` | `(samples) → (μ, σ²)` | MOM für Normalverteilung |
| `mom_gamma` | `(samples) → (α, β)` | MOM für Gamma-Verteilung |
| `mom_beta` | `(samples) → (α, β)` | MOM für Beta-Verteilung |

### `SufficiencyTheory`

Suffizienz, Vollständigkeit und Rao-Blackwell-Verbesserung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `is_sufficient_statistic_normal_mean` | `(samples, known_sigma) → dict` | Prüft Suffizienz von $\bar{X}$ für $\mu$ |
| `factorization_criterion_demo` | `(samples, ...) → dict` | Fisher-Neyman-Faktorisierung (Demo) |
| `complete_statistic_demo` | `() → dict` | Vollständigkeit und Lehmann-Scheffé-Satz |
| `rao_blackwell` | `(estimator_fn, sufficient_stat_fn, samples, ...) → dict` | Rao-Blackwell-Verbesserung |

### `CramerRaoBound`

Cramér-Rao-Ungleichung und Effizienz.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `cramer_rao_normal_mean` | `(n, sigma) → float` | CRB für $\hat{\mu}$ bei Normalverteilung: $\sigma^2/n$ |
| `cramer_rao_poisson` | `(n, lam=1.0) → float` | CRB für $\hat{\lambda}$ bei Poisson: $\lambda/n$ |
| `efficiency` | `(variance, cramer_rao_bound) → float` | Effizienz $e = \text{CRB}/\operatorname{Var}$ |
| `umvue_normal_demo` | `(samples) → dict` | UMVUE-Demo: $\bar{X}$ und $S^2$ |
| `fisher_information` | `(log_likelihood_fn, theta, ...) → float` | Numerische Fisher-Information |

### `ConfidenceIntervals`

Konfidenzintervalle für verschiedene Parameter und Methoden.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `ci_normal_mean` | `(samples, alpha, known_sigma) → (lo, hi)` | KI für $\mu$: $t$-Test oder $z$-Test |
| `ci_normal_variance` | `(samples, alpha) → (lo, hi)` | KI für $\sigma^2$: Chi-Quadrat-Verteilung |
| `ci_proportion` | `(n, k, alpha, method) → (lo, hi)` | KI für Anteil: Wald oder Wilson |
| `ci_bootstrap` | `(samples, estimator_fn, alpha, B) → (lo, hi)` | Bootstrap-Konfidenzintervall |
| `ci_exponential` | `(samples, alpha) → (lo, hi)` | KI für Exponentialrate $\lambda$ |

### `HypothesisTesting`

Hypothesentests auf Basis des Neyman-Pearson-Rahmens.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `neyman_pearson_lemma_demo` | `(samples, theta0, theta1, alpha, ...) → dict` | LRT-Demo nach Neyman-Pearson |
| `ump_test_normal_mean` | `(samples, mu0, sigma, alpha, alternative) → dict` | UMP-Test für $\mu$ (Normal, $\sigma$ bekannt) |
| `type_i_type_ii_error` | `(alpha, beta) → dict` | Fehler 1. und 2. Art, Trennschärfe |
| `power_function` | `(theta_vals, mu0, sigma, n, alpha, ...) → dict` | Trennschärfefunktion $\beta(\theta)$ |
| `sequential_probability_ratio_test` | `(samples, theta0, theta1, alpha, beta_err) → dict` | Wald-SPRT mit Stoppregeln |

### `DecisionTheory`

Statistische Entscheidungstheorie.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `loss_function_squared_error` | `(estimate, true_value) → float` | Quadratischer Verlust $L = (\hat{\theta}-\theta)^2$ |
| `loss_function_absolute` | `(estimate, true_value) → float` | Absoluter Verlust $L = |\hat{\theta}-\theta|$ |
| `risk_function` | `(estimator_fn, theta, sample_dist_fn, ...) → float` | Risikofunktion $R(\hat{\theta}, \theta)$ per Simulation |
| `bayes_risk` | `(estimator_fn, prior_samples, ...) → float` | Bayes-Risiko $r(\hat{\theta})$ |
| `minimax_estimator_demo` | `(...) → dict` | Minimax-Schätzer (Demo) |
| `admissibility_check_demo` | `() → dict` | Zulässigkeitsprüfung (Demo) |

---

## Beispiele

```python
import numpy as np
from estimation_theory import (MaximumLikelihoodEstimator, CramerRaoBound,
                                ConfidenceIntervals, HypothesisTesting)

# -- MLE für Normalverteilung --
rng = np.random.default_rng(42)
samples = rng.normal(loc=5.0, scale=2.0, size=100)

mle = MaximumLikelihoodEstimator()
mu_hat, var_hat = mle.mle_normal(samples)
print(f"μ̂ = {mu_hat:.3f}, σ̂² = {var_hat:.3f}")  # ≈ 5.0, 4.0

# -- Cramér-Rao-Schranke --
crb = CramerRaoBound()
bound = crb.cramer_rao_normal_mean(n=100, sigma=2.0)
print(f"CRB = {bound:.4f}")                        # 4/100 = 0.04
eff = crb.efficiency(variance=var_hat/100, cramer_rao_bound=bound)
print(f"Effizienz = {eff:.3f}")                    # ≈ 1.0 (MLE ist effizient)

# -- 95%-Konfidenzintervall --
ci = ConfidenceIntervals()
lo, hi = ci.ci_normal_mean(samples, alpha=0.05)
print(f"KI: [{lo:.3f}, {hi:.3f}]")

# -- Neyman-Pearson-Test --
ht = HypothesisTesting()
result = ht.ump_test_normal_mean(samples, mu0=4.5, sigma=2.0, alpha=0.05, alternative="greater")
print(result["decision"])  # "reject" oder "fail_to_reject"

# -- Bootstrap-KI --
lo_b, hi_b = ci.ci_bootstrap(samples, estimator_fn=np.mean, alpha=0.05, B=2000)
print(f"Bootstrap-KI: [{lo_b:.3f}, {hi_b:.3f}]")
```

---

## Tests

**Testdatei:** `tests/test_estimation_theory.py`
**Abdeckung:** Alle Klassen – MLE-Konsistenz, CRB-Verifikation, KI-Überdeckungswahrscheinlichkeit, Teststärke, Rao-Blackwell-Verbesserung, Entscheidungstheorie-Demos

---

## Implementierungshinweise

- **MLE-Konsistenz:** Für großes $n$ konvergieren alle MLEs gegen den wahren Parameter ($\sqrt{n}$-Konsistenz, asymptotische Normalität).
- **Fisher-Information numerisch:** Per zentralem Differenzenquotient der Log-Likelihood-Funktion.
- **Bootstrap:** Percentile-Methode; für bessere Genauigkeit ist BCa-Bootstrap (bias-corrected accelerated) möglich, aber nicht implementiert.
- **SPRT:** Wald'sche Näherungen $A \approx \beta/(1-\alpha)$, $B \approx (1-\beta)/\alpha$ für die Grenzen werden verwendet.
- **Minimax/Admissibility:** Als Demo implementiert; vollständige numerische Lösung erfordert Spieltheorie (vgl. `operations_research.py`).
