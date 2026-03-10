# Modul: Stochastik (`stochastic.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Dieses Modul implementiert grundlegende und fortgeschrittene stochastische Prozesse: diskrete
und zeitkontinuierliche Markov-Ketten, Brownsche Bewegung, Itô-Kalkül, stochastische
Differentialgleichungen, Ergodentheorie, Poisson-Prozesse und Gaußsche Prozesse — sowie
klassische Resultate wie den Zentralen Grenzwertsatz und die Gambler's-Ruin-Formel.

## Mathematischer Hintergrund

### Markov-Ketten

Eine **diskrete Markov-Kette** $(X_n)_{n \geq 0}$ erfüllt die Markov-Eigenschaft:

$$P(X_{n+1} = j \mid X_n = i, X_{n-1}, \ldots) = P_{ij}$$

Die **stationäre Verteilung** $\pi$ löst das lineare Gleichungssystem:

$$\pi P = \pi, \quad \sum_j \pi_j = 1$$

Für eine **ergodische** (irreduzible, aperiodische) Kette:

$$\lim_{n \to \infty} P^n_{ij} = \pi_j$$

**Mittlere erste Rückkehrzeit** zum Zustand $j$: $m_{jj} = 1/\pi_j$.

### Zeitkontinuierliche Markov-Ketten

Eine zeitkontinuierliche Kette wird durch die **Ratenmatrix** (Generator) $Q$ beschrieben:

$$Q_{ij} \geq 0 \text{ für } i \neq j, \quad Q_{ii} = -\sum_{j \neq i} Q_{ij}$$

Übergangswahrscheinlichkeiten via Matrixexponential:

$$P(t) = e^{Qt}$$

### Brownsche Bewegung

Ein **Wiener-Prozess** $W(t)$ erfüllt:
- $W(0) = 0$
- $W(t) - W(s) \sim \mathcal{N}(0, t-s)$ für $t > s$
- Unabhängige Zuwächse

Kovarianzstruktur und quadratische Variation:

$$\operatorname{Cov}(W(t), W(s)) = \min(t, s), \quad [W, W]_T = T \text{ (fast sicher)}$$

### Itô-Kalkül

Das **Itô-Integral** ist ein stochastisches Integral bezüglich $W(t)$:

$$\int_0^T f(t) \, dW(t) = \lim_{\|P\| \to 0} \sum_{k} f(t_k) \big[W(t_{k+1}) - W(t_k)\big]$$

**Itô-Isometrie:**

$$\mathbb{E}\left[\left(\int_0^T f \, dW\right)^2\right] = \int_0^T \mathbb{E}[f(t)^2] \, dt$$

**Itô-Formel** (stochastische Kettenregel) für $Y_t = g(t, X_t)$:

$$dY = \frac{\partial g}{\partial t} dt + \frac{\partial g}{\partial x} dX + \frac{1}{2} \frac{\partial^2 g}{\partial x^2} (dX)^2$$

wobei $(dW)^2 = dt$.

**Stratonovich-Integral** (mittelpunktsbasiert, Kettenregel wie gewöhnlich):

$$\int_0^T f \circ dW = \int_0^T f \, dW + \frac{1}{2}\int_0^T f'(t) \, dt$$

### Stochastische Differentialgleichungen (SDE)

Itô-SDE in allgemeiner Form:

$$dX_t = \mu(t, X_t) \, dt + \sigma(t, X_t) \, dW_t$$

**Euler-Maruyama-Verfahren** (Ordnung $1/2$):

$$X_{n+1} = X_n + \mu(t_n, X_n) \Delta t + \sigma(t_n, X_n) \Delta W_n$$

**Milstein-Verfahren** (Ordnung $1$, berücksichtigt $\sigma$-Ableitung):

$$X_{n+1} = X_n + \mu \Delta t + \sigma \Delta W_n + \frac{1}{2} \sigma \sigma' \left[(\Delta W_n)^2 - \Delta t\right]$$

**Geometrische Brownsche Bewegung** (Black-Scholes-Modell):

$$dS = \mu S \, dt + \sigma S \, dW \implies S(t) = S_0 \exp\!\left(\!\left(\mu - \frac{\sigma^2}{2}\right)t + \sigma W(t)\right)$$

### Ergodentheorie

**Birkhoff'sches Ergodentheorem:** Für eine maßerhaltende Abbildung $T$ und $f \in L^1$:

$$\frac{1}{N} \sum_{n=0}^{N-1} f(T^n x) \xrightarrow{N \to \infty} \int f \, d\mu \quad \text{(f.s.)}$$

**Lyapunov-Exponent** misst die exponentielle Divergenz von Trajektorien:

$$\lambda = \lim_{N \to \infty} \frac{1}{N} \sum_{n=0}^{N-1} \ln |f'(x_n)|$$

**Logistische Abbildung:** $x_{n+1} = r x_n (1 - x_n)$ — klassisches Beispiel für Chaos ($r > 3.57...$).

### Poisson-Prozess

Ein **Poisson-Prozess** $N(t)$ mit Rate $\lambda$ hat Zuwächse $\Delta N \sim \text{Poisson}(\lambda \Delta t)$:

$$P(N(t) = k) = \frac{(\lambda t)^k e^{-\lambda t}}{k!}, \quad \mathbb{E}[N(t)] = \lambda t, \quad \operatorname{Var}(N(t)) = \lambda t$$

### Gaußscher Prozess

Ein **Gaußscher Prozess** $f \sim \mathcal{GP}(m, k)$ ist vollständig durch Mittelwert $m(x)$
und Kovarianzfunktion (Kern) $k(x, x')$ bestimmt. RBF-Kern (Gauß-Kern):

$$k(x, x') = \sigma_f^2 \exp\!\left(-\frac{\|x - x'\|^2}{2\ell^2}\right)$$

### Gambler's Ruin

Ruinwahrscheinlichkeit für symmetrisches Random Walk (Gewinnwahrscheinlichkeit $p$, Startkapital $k$, Ziel $N$):

$$P(\text{Ruin} \mid X_0 = k) = \frac{(q/p)^k - (q/p)^N}{1 - (q/p)^N} \quad (p \neq 1/2)$$

---

## Klassen und Methoden

### `MarkovChain`

Diskrete homogene Markov-Kette mit endlichem Zustandsraum.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(transition_matrix: np.ndarray)` | Erzeugt Kette aus zeilenstochastischer Matrix $P$ |
| `stationary_distribution` | `() -> np.ndarray` | Berechnet $\pi$ mit $\pi P = \pi$ (linke Eigenvektoren) |
| `simulate` | `(n_steps, initial_state) -> np.ndarray` | Simuliert $n$ Schritte ab Startzustand |
| `n_step_distribution` | `(initial_dist, n) -> np.ndarray` | Verteilung nach $n$ Schritten: $\pi^{(0)} P^n$ |
| `is_ergodic` | `() -> bool` | Prüft Ergodizität (irreduzibel + aperiodisch) |
| `mean_first_passage_time` | `(i: int, j: int) -> float` | Mittlere Übergangszeit von $i$ nach $j$ |
| `absorption_probabilities` | `(...) -> np.ndarray` | Absorptionswahrscheinlichkeiten in absorbierende Zustände |

### `ContinuousTimeMarkovChain`

Zeitkontinuierliche Markov-Kette via Generator-Matrix und Matrixexponential.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(rate_matrix: np.ndarray)` | Erzeugt Kette aus Ratenmatrix $Q$ |
| `transition_matrix` | `(t: float) -> np.ndarray` | Übergangsmatrix $P(t) = e^{Qt}$ |
| `stationary_distribution` | `() -> np.ndarray` | Stationäre Verteilung (Kern von $Q^\top$) |
| `simulate` | `(t_max, initial_state) -> (times, states)` | Gillespie-Simulation |

### `BrownianMotion`

Simulation und Analyse der Brownschen Bewegung (Wiener-Prozess).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `simulate` | `(T, n_steps) -> (t, W)` | Simuliert Pfad via inkrementelle Gauß-Zuwächse |
| `covariance` | `(t: float, s: float) -> float` | $\operatorname{Cov}(W(t), W(s)) = \min(t,s)$ |
| `quadratic_variation` | `(T, n_steps) -> float` | Empirische quadratische Variation $\approx T$ |

### `ItoIntegral`

Itô- und Stratonovich-Integrale, Itô-Isometrie und Itô-Formel.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `ito_integral` | `(f: Callable, T, n_steps) -> float` | Numerisches Itô-Integral $\int_0^T f(t)\,dW$ |
| `stratonovich_integral` | `(f: Callable, T, n_steps) -> float` | Numerisches Stratonovich-Integral |
| `ito_formula_demo` | `(g: Callable, ...) -> dict` | Demonstriert Itô-Formel für $g(W_t)$ |
| `ito_isometry_check` | `(f: Callable, T, n_steps) -> dict` | Numerische Überprüfung der Itô-Isometrie |

### `StochasticDifferentialEquation`

Numerische Löser für SDEs: Euler-Maruyama, Milstein, geometrische Brownsche Bewegung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `euler_maruyama` | `(mu, sigma, X0, T, n_steps) -> (t, X)` | Euler-Maruyama-Verfahren (starke Ordnung 0.5) |
| `milstein` | `(mu, sigma, sigma_prime, X0, T, n_steps) -> (t, X)` | Milstein-Verfahren (starke Ordnung 1.0) |
| `geometric_brownian_motion` | `(mu, sigma, S0, T, n_steps) -> (t, S)` | Geometrische BM (Black-Scholes-Modell) |

### `ErgodicTheory`

Birkhoff'sches Ergodentheorem, Lyapunov-Exponent und logistische Abbildung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `birkhoff_ergodic_theorem_demo` | `(f, T_map, x0, N) -> dict` | Numerische Demonstration des Birkhoff-Theorems |
| `time_average` | `(f_values: np.ndarray) -> float` | Zeitdurchschnitt $\frac{1}{N}\sum f(x_n)$ |
| `is_mixing_demo` | `(sequence, lag=10) -> bool` | Prüft Mixing-Eigenschaft via Autokorrelation |
| `logistic_map` | `(r, x0, n_steps) -> np.ndarray` | Iteration der logistischen Abbildung |
| `lyapunov_exponent` | `(r, x0, n_steps) -> float` | Lyapunov-Exponent $\lambda$ der logistischen Abbildung |

### `PoissonProcess`

Poisson-Prozess mit Rate $\lambda$: Simulation, PMF und Momente.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(rate: float)` | Erzeugt Poisson-Prozess mit Rate $\lambda$ |
| `simulate` | `(T, n_paths=1) -> list` | Simuliert Ankunftszeitfolgen bis Zeit $T$ |
| `pmf` | `(k: int, t: float) -> float` | $P(N(t) = k) = \frac{(\lambda t)^k e^{-\lambda t}}{k!}$ |
| `mean` | `(t: float) -> float` | $\mathbb{E}[N(t)] = \lambda t$ |
| `variance` | `(t: float) -> float` | $\operatorname{Var}(N(t)) = \lambda t$ |

### `GaussianProcess`

Gaußscher Prozess mit RBF-Kern: Stichproben und Prädiktion.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(mean_fn: Callable, kernel_fn: Callable)` | Erzeugt GP mit Mittelwert- und Kernfunktion |
| `sample` | `(X: np.ndarray, n_samples=1) -> np.ndarray` | Zieht Pfade aus dem GP |
| `rbf_kernel` | `(x1, x2, length_scale, sigma_f) -> float` | RBF-Kern $k(x,x')$ |

---

## Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `random_walk_1d` | `(n_steps, p=0.5, x0=0) -> np.ndarray` | 1D Random Walk mit Schrittwahrscheinlichkeit $p$ |
| `gambler_ruin_probability` | `(p, k, N) -> float` | Ruinwahrscheinlichkeit bei Startkapital $k$ und Ziel $N$ |
| `central_limit_theorem_demo` | `(dist_fn, n_samples, n_trials) -> dict` | Numerische CLT-Demonstration |

---

## Beispiele

```python
import numpy as np
from stochastic import (
    MarkovChain, BrownianMotion, StochasticDifferentialEquation,
    PoissonProcess, ErgodicTheory, gambler_ruin_probability
)

# Markov-Kette: Wetter (Sonne / Regen)
P = np.array([[0.8, 0.2],
              [0.4, 0.6]])
mc = MarkovChain(P)
print("Stationäre Verteilung:", mc.stationary_distribution())  # [0.667, 0.333]
print("Ergodisch:", mc.is_ergodic())                           # True

# Brownsche Bewegung
bm = BrownianMotion()
t, W = bm.simulate(T=1.0, n_steps=1000)
print("Cov(W(0.3), W(0.7)):", bm.covariance(0.3, 0.7))  # 0.3

# Geometrische Brownsche Bewegung (Aktienkurs)
sde = StochasticDifferentialEquation()
mu, sigma, S0 = 0.05, 0.2, 100.0
t, S = sde.geometric_brownian_motion(mu, sigma, S0, T=1.0, n_steps=252)
print(f"Endpreis S(1): {S[-1]:.2f}")

# Poisson-Prozess
pp = PoissonProcess(rate=3.0)
print("P(N(2)=5):", pp.pmf(5, 2.0))   # Poisson(6)
print("E[N(2)]:", pp.mean(2.0))        # 6.0

# Logistische Abbildung und Lyapunov-Exponent
et = ErgodicTheory()
lam = et.lyapunov_exponent(r=3.9, x0=0.5, n_steps=10000)
print(f"Lyapunov-Exponent (r=3.9): {lam:.4f}")  # > 0 → Chaos

# Gambler's Ruin
p_ruin = gambler_ruin_probability(p=0.45, k=10, N=20)
print(f"Ruinwahrscheinlichkeit: {p_ruin:.4f}")
```

---

## Technische Hinweise

- **Stationäre Verteilung** der Markov-Kette wird über den linken Eigenvektor zum Eigenwert 1 berechnet.
- **Matrixexponential** für zeitkontinuierliche Ketten via `scipy.linalg.expm`.
- **Gillespie-Algorithmus** für zeitkontinuierliche Ketten: Wartezeiten sind exponentialverteilt.
- **Milstein-Verfahren** benötigt die Ableitung $\sigma'$ explizit als Parameter.
- **Gaußscher Prozess**: Stichproben via Cholesky-Zerlegung der Kovarianzmatrix.
- **Lyapunov-Exponent** wird als Ergodendurchschnitt von $\ln|f'(x_n)|$ berechnet.
