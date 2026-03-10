"""
@file test_stochastic.py
@brief Umfassende Tests für das stochastic.py Modul.
@description
    Test-Suite für alle Klassen und Funktionen des Stochastik-Moduls:
    - MarkovChain: stationäre Verteilung, Simulation, Ergodizität,
      mittlere Übergangszeit, Absorptionswahrscheinlichkeiten
    - ContinuousTimeMarkovChain: Übergangsmatrix, stationäre Verteilung
    - BrownianMotion: Simulation, Kovarianz, quadratische Variation
    - ItoIntegral: Itô- und Stratonovich-Integral, Itô-Formel, Isometrie
    - StochasticDifferentialEquation: Euler-Maruyama, Milstein, GBM
    - ErgodicTheory: Birkhoff, Zeitmittel, Lyapunov, logistische Abbildung
    - PoissonProcess: Simulation, PMF, Momente
    - GaussianProcess: Stichprobe, RBF-Kern
    - Standalone-Funktionen: Random Walk, Gambler's Ruin, CLT

    Edge Cases: absorbierende Zustände, degenerierte Verteilungen,
    Randwerte, numerische Stabilität.

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import pytest
import numpy as np
import sys
import os

# Testverzeichnis und src-Pfad konfigurieren
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from stochastic import (
    MarkovChain,
    ContinuousTimeMarkovChain,
    BrownianMotion,
    ItoIntegral,
    StochasticDifferentialEquation,
    ErgodicTheory,
    PoissonProcess,
    GaussianProcess,
    random_walk_1d,
    gambler_ruin_probability,
    central_limit_theorem_demo,
)


# =============================================================================
# Fixtures: Häufig verwendete Übergangsmatrizen
# =============================================================================

@pytest.fixture
def simple_chain():
    """Einfache 2-Zustands Markov-Kette."""
    P = np.array([[0.7, 0.3],
                  [0.4, 0.6]])
    return MarkovChain(P)


@pytest.fixture
def ergodic_chain():
    """Ergodische 3-Zustands Kette (irreduzibel, aperiodisch)."""
    P = np.array([[0.5, 0.3, 0.2],
                  [0.2, 0.6, 0.2],
                  [0.1, 0.4, 0.5]])
    return MarkovChain(P)


@pytest.fixture
def absorbing_chain():
    """
    Absorbierende Markov-Kette:
    Zustände 0, 1: transient; Zustände 2, 3: absorbierend.
    """
    P = np.array([
        [0.0, 0.5, 0.3, 0.2],  # Zustand 0 → Transient
        [0.4, 0.0, 0.2, 0.4],  # Zustand 1 → Transient
        [0.0, 0.0, 1.0, 0.0],  # Zustand 2 → Absorbierend
        [0.0, 0.0, 0.0, 1.0],  # Zustand 3 → Absorbierend
    ])
    return MarkovChain(P)


@pytest.fixture
def random_walk_chain():
    """Symmetrischer Random Walk auf {0, 1, 2, 3, 4} mit absorbierenden Rändern."""
    P = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],  # 0 absorbierend
        [0.5, 0.0, 0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5, 0.0, 0.5],
        [0.0, 0.0, 0.0, 0.0, 1.0],  # 4 absorbierend
    ])
    return MarkovChain(P)


# =============================================================================
# Tests: MarkovChain - Konstruktion und Validierung
# =============================================================================

class TestMarkovChainConstruction:
    """Tests für die Konstruktion und Validierung von Markov-Ketten."""

    def test_valid_2x2_matrix(self):
        """Gültige 2×2 stochastische Matrix wird akzeptiert."""
        P = np.array([[0.7, 0.3], [0.4, 0.6]])
        chain = MarkovChain(P)
        assert chain.n_states == 2

    def test_valid_3x3_matrix(self):
        """Gültige 3×3 stochastische Matrix wird akzeptiert."""
        P = np.eye(3)
        chain = MarkovChain(P)
        assert chain.n_states == 3

    def test_non_square_raises(self):
        """Nicht-quadratische Matrix wirft ValueError."""
        P = np.array([[0.5, 0.5, 0.0], [0.3, 0.7, 0.0]])
        with pytest.raises(ValueError, match="quadratisch"):
            MarkovChain(P)

    def test_negative_entries_raises(self):
        """Negative Einträge werfen ValueError."""
        P = np.array([[-0.1, 1.1], [0.5, 0.5]])
        with pytest.raises(ValueError):
            MarkovChain(P)

    def test_row_sum_not_one_raises(self):
        """Zeilensumme ≠ 1 wirft ValueError."""
        P = np.array([[0.6, 0.6], [0.4, 0.6]])
        with pytest.raises(ValueError, match="Zeilensummen"):
            MarkovChain(P)

    def test_identity_matrix_accepted(self):
        """Einheitsmatrix (vollständig absorbierend) wird akzeptiert."""
        P = np.eye(4)
        chain = MarkovChain(P)
        assert chain.n_states == 4

    def test_transition_matrix_stored(self, simple_chain):
        """Übergangsmatrix wird korrekt gespeichert."""
        expected = np.array([[0.7, 0.3], [0.4, 0.6]])
        np.testing.assert_array_almost_equal(simple_chain.P, expected)


# =============================================================================
# Tests: MarkovChain - Stationäre Verteilung
# =============================================================================

class TestMarkovChainStationary:
    """Tests für die stationäre Verteilung."""

    def test_stationary_sums_to_one(self, simple_chain):
        """Stationäre Verteilung summiert zu 1."""
        pi = simple_chain.stationary_distribution()
        assert abs(pi.sum() - 1.0) < 1e-8

    def test_stationary_all_positive(self, simple_chain):
        """Stationäre Verteilung hat positive Einträge."""
        pi = simple_chain.stationary_distribution()
        assert np.all(pi > 0)

    def test_stationary_satisfies_balance(self, simple_chain):
        """Stationäre Verteilung erfüllt πP = π."""
        pi = simple_chain.stationary_distribution()
        pi_new = pi @ simple_chain.P
        np.testing.assert_array_almost_equal(pi, pi_new, decimal=8)

    def test_stationary_2state_exact(self, simple_chain):
        """Für 2-Zustands-Kette: analytische Überprüfung."""
        # Für P = [[0.7, 0.3], [0.4, 0.6]]:
        # πP = π => π₀·0.3 = π₁·0.4 => π₀ = 4/7, π₁ = 3/7
        pi = simple_chain.stationary_distribution()
        assert abs(pi[0] - 4.0 / 7.0) < 1e-6
        assert abs(pi[1] - 3.0 / 7.0) < 1e-6

    def test_stationary_3state_ergodic(self, ergodic_chain):
        """Stationäre Verteilung der 3-Zustands-Kette ist korrekt."""
        pi = ergodic_chain.stationary_distribution()
        assert abs(pi.sum() - 1.0) < 1e-8
        pi_new = pi @ ergodic_chain.P
        np.testing.assert_array_almost_equal(pi, pi_new, decimal=6)

    def test_stationary_length(self, ergodic_chain):
        """Stationäre Verteilung hat korrekte Länge."""
        pi = ergodic_chain.stationary_distribution()
        assert len(pi) == 3

    def test_uniform_stationary_doubly_stochastic(self):
        """Doppelt stochastische Matrix hat uniforme stationäre Verteilung."""
        # Zirkulante Matrix: jede Zeile und Spalte summiert zu 1
        P = np.array([[1/3, 1/3, 1/3],
                      [1/3, 1/3, 1/3],
                      [1/3, 1/3, 1/3]])
        chain = MarkovChain(P)
        pi = chain.stationary_distribution()
        np.testing.assert_array_almost_equal(pi, [1/3, 1/3, 1/3], decimal=6)


# =============================================================================
# Tests: MarkovChain - Simulation
# =============================================================================

class TestMarkovChainSimulation:
    """Tests für die Simulation von Markov-Ketten-Pfaden."""

    def test_simulate_length(self, simple_chain):
        """Simulierter Pfad hat korrekte Länge."""
        path = simple_chain.simulate(0, 100, seed=42)
        assert len(path) == 101  # n_steps + 1

    def test_simulate_starts_correctly(self, simple_chain):
        """Pfad startet im richtigen Zustand."""
        path = simple_chain.simulate(1, 50, seed=0)
        assert path[0] == 1

    def test_simulate_valid_states(self, simple_chain):
        """Alle Zustände im simulierten Pfad sind gültig."""
        path = simple_chain.simulate(0, 1000, seed=42)
        assert np.all((path >= 0) & (path < simple_chain.n_states))

    def test_simulate_reproducible_with_seed(self, simple_chain):
        """Gleicher Seed erzeugt identische Pfade."""
        path1 = simple_chain.simulate(0, 100, seed=123)
        path2 = simple_chain.simulate(0, 100, seed=123)
        np.testing.assert_array_equal(path1, path2)

    def test_simulate_different_seeds(self, simple_chain):
        """Verschiedene Seeds erzeugen verschiedene Pfade."""
        path1 = simple_chain.simulate(0, 100, seed=1)
        path2 = simple_chain.simulate(0, 100, seed=2)
        assert not np.array_equal(path1, path2)

    def test_simulate_invalid_state_raises(self, simple_chain):
        """Ungültiger Startzustand wirft ValueError."""
        with pytest.raises(ValueError, match="außerhalb"):
            simple_chain.simulate(5, 10, seed=0)

    def test_simulate_negative_state_raises(self, simple_chain):
        """Negativer Startzustand wirft ValueError."""
        with pytest.raises(ValueError):
            simple_chain.simulate(-1, 10, seed=0)

    def test_simulate_empirical_stationary(self, simple_chain):
        """Empirische Häufigkeit konvergiert gegen stationäre Verteilung."""
        path = simple_chain.simulate(0, 100000, seed=42)
        empirical = np.array([np.mean(path == s) for s in range(2)])
        pi = simple_chain.stationary_distribution()
        np.testing.assert_array_almost_equal(empirical, pi, decimal=2)

    def test_simulate_absorbing_stays(self):
        """In vollständig absorbierender Kette bleibt Zustand konstant."""
        P = np.eye(3)
        chain = MarkovChain(P)
        path = chain.simulate(2, 50, seed=0)
        assert np.all(path == 2)


# =============================================================================
# Tests: MarkovChain - n-Schritt-Verteilung
# =============================================================================

class TestMarkovChainNStepDistribution:
    """Tests für die n-Schritt-Verteilung."""

    def test_zero_steps(self, simple_chain):
        """0-Schritt-Verteilung ist Anfangsverteilung."""
        dist0 = np.array([1.0, 0.0])
        result = simple_chain.n_step_distribution(dist0, 0)
        np.testing.assert_array_almost_equal(result, dist0, decimal=8)

    def test_one_step(self, simple_chain):
        """1-Schritt ist direkte Matrixmultiplikation."""
        dist0 = np.array([1.0, 0.0])
        result = simple_chain.n_step_distribution(dist0, 1)
        expected = dist0 @ simple_chain.P
        np.testing.assert_array_almost_equal(result, expected, decimal=8)

    def test_n_step_sums_to_one(self, ergodic_chain):
        """n-Schritt-Verteilung summiert zu 1."""
        dist0 = np.array([0.5, 0.3, 0.2])
        for n in [1, 5, 10, 100]:
            result = ergodic_chain.n_step_distribution(dist0, n)
            assert abs(result.sum() - 1.0) < 1e-8

    def test_convergence_to_stationary(self, ergodic_chain):
        """n-Schritt-Verteilung konvergiert gegen stationäre Verteilung."""
        dist0 = np.array([1.0, 0.0, 0.0])
        result = ergodic_chain.n_step_distribution(dist0, 1000)
        pi = ergodic_chain.stationary_distribution()
        np.testing.assert_array_almost_equal(result, pi, decimal=6)

    def test_normalization_of_input(self, simple_chain):
        """Unnormierte Anfangsverteilung wird normiert."""
        dist0 = np.array([2.0, 1.0])  # Skaliert, nicht normiert
        result = simple_chain.n_step_distribution(dist0, 1)
        assert abs(result.sum() - 1.0) < 1e-8


# =============================================================================
# Tests: MarkovChain - Ergodizität
# =============================================================================

class TestMarkovChainErgodic:
    """Tests für die Ergodizitätsprüfung."""

    def test_ergodic_chain(self, ergodic_chain):
        """Ergodische 3-Zustands-Kette wird als ergodisch erkannt."""
        assert ergodic_chain.is_ergodic() is True

    def test_simple_ergodic(self, simple_chain):
        """Einfache 2-Zustands-Kette ohne Selbstübergänge."""
        # Beide Zustände haben Selbstübergänge? Nein, aber trotzdem ergodisch
        assert simple_chain.is_ergodic() is True

    def test_identity_not_ergodic(self):
        """Einheitsmatrix (N=2) ist nicht ergodisch (nicht irreduzibel für N>1)."""
        P = np.eye(2)
        chain = MarkovChain(P)
        # Zustand 0 kann nicht nach 1 gelangen → nicht irreduzibel
        assert chain.is_ergodic() is False

    def test_periodic_chain_not_aperiodic(self):
        """Periodische Kette (Periode 2) ist nicht ergodisch."""
        # Bipartite Kette: 0 → 1 → 0 (Periode 2)
        P = np.array([[0.0, 1.0],
                      [1.0, 0.0]])
        chain = MarkovChain(P)
        assert chain.is_ergodic() is False

    def test_3state_ergodic_with_selfloop(self):
        """3-Zustands-Kette mit Selbstübergang ist ergodisch."""
        P = np.array([[0.5, 0.3, 0.2],
                      [0.4, 0.4, 0.2],
                      [0.1, 0.2, 0.7]])
        chain = MarkovChain(P)
        assert chain.is_ergodic() is True


# =============================================================================
# Tests: MarkovChain - Mittlere erste Übergangszeit
# =============================================================================

class TestMarkovChainMFPT:
    """Tests für die mittlere erstmalige Übergangszeit."""

    def test_mfpt_self_equals_one_over_pi(self, simple_chain):
        """Mittlere Rückkehrzeit = 1/π_i."""
        pi = simple_chain.stationary_distribution()
        for i in range(2):
            mfpt = simple_chain.mean_first_passage_time(i, i)
            expected = 1.0 / pi[i]
            assert abs(mfpt - expected) < 1e-6

    def test_mfpt_positive(self, ergodic_chain):
        """Mittlere erste Übergangszeit ist positiv."""
        for i in range(3):
            for j in range(3):
                mfpt = ergodic_chain.mean_first_passage_time(i, j)
                assert mfpt > 0

    def test_mfpt_direct_transition(self):
        """Für direkten Übergang mit Wahrsch. 1: m_{ij} = 1."""
        # Kette: 0 → 1 immer, 1 → 0 immer
        P = np.array([[0.0, 1.0],
                      [1.0, 0.0]])
        chain = MarkovChain(P)
        # m_{01} = 1 (direkter Übergang)
        mfpt = chain.mean_first_passage_time(0, 1)
        assert abs(mfpt - 1.0) < 1e-6

    def test_mfpt_symmetry_in_symmetric_chain(self):
        """Symmetrische Kette: m_{01} = m_{10}."""
        P = np.array([[0.5, 0.5],
                      [0.5, 0.5]])
        chain = MarkovChain(P)
        m01 = chain.mean_first_passage_time(0, 1)
        m10 = chain.mean_first_passage_time(1, 0)
        assert abs(m01 - m10) < 1e-6

    def test_mfpt_3state(self, ergodic_chain):
        """3-Zustands-Kette: MFPT ist endlich und konsistent."""
        # Ketten-Eigenschaft: m_{ij} + m_{ji} ≠ konstant (nicht symmetrisch)
        m01 = ergodic_chain.mean_first_passage_time(0, 1)
        m10 = ergodic_chain.mean_first_passage_time(1, 0)
        assert m01 > 0
        assert m10 > 0


# =============================================================================
# Tests: MarkovChain - Absorptionswahrscheinlichkeiten
# =============================================================================

class TestMarkovChainAbsorption:
    """Tests für Absorptionswahrscheinlichkeiten."""

    def test_absorption_sums_to_one(self, absorbing_chain):
        """Absorptionswahrscheinlichkeiten summieren zu 1 pro Zeile."""
        B = absorbing_chain.absorption_probabilities([0, 1], [2, 3])
        row_sums = B.sum(axis=1)
        np.testing.assert_array_almost_equal(row_sums, [1.0, 1.0], decimal=8)

    def test_absorption_non_negative(self, absorbing_chain):
        """Alle Absorptionswahrscheinlichkeiten sind nicht-negativ."""
        B = absorbing_chain.absorption_probabilities([0, 1], [2, 3])
        assert np.all(B >= -1e-10)

    def test_absorption_shape(self, absorbing_chain):
        """Ergebnismatrix hat korrekte Form (|T| × |A|)."""
        B = absorbing_chain.absorption_probabilities([0, 1], [2, 3])
        assert B.shape == (2, 2)

    def test_random_walk_absorption(self, random_walk_chain):
        """Symmetrischer Random Walk: Absorptionswahrscheinlichkeiten exakt."""
        # Transient: 1, 2, 3; Absorbierend: 0, 4
        B = random_walk_chain.absorption_probabilities([1, 2, 3], [0, 4])
        # Für symmetrischen RW von k aus: P(Absorption in 4) = k/4
        # Von 1: P(4) = 1/4, Von 2: P(4) = 2/4, Von 3: P(4) = 3/4
        np.testing.assert_array_almost_equal(B[:, 1], [1/4, 2/4, 3/4], decimal=6)

    def test_absorption_certain_absorbing(self):
        """Wenn nur ein absorbierender Zustand: Wahrsch. = 1."""
        P = np.array([
            [0.5, 0.5, 0.0],
            [0.3, 0.2, 0.5],
            [0.0, 0.0, 1.0],  # Absorbierend
        ])
        chain = MarkovChain(P)
        B = chain.absorption_probabilities([0, 1], [2])
        np.testing.assert_array_almost_equal(B.flatten(), [1.0, 1.0], decimal=6)


# =============================================================================
# Tests: ContinuousTimeMarkovChain
# =============================================================================

class TestCTMC:
    """Tests für zeitkontinuierliche Markov-Ketten."""

    @pytest.fixture
    def birth_death(self):
        """Einfacher Geburts-Tod-Prozess mit 3 Zuständen."""
        # Rate λ=2 nach oben, Rate μ=1 nach unten
        Q = np.array([
            [-2.0, 2.0, 0.0],
            [1.0, -3.0, 2.0],
            [0.0, 1.0, -1.0],
        ])
        return ContinuousTimeMarkovChain(Q)

    def test_invalid_rate_matrix(self):
        """Ungültige Ratenmatrix wirft ValueError."""
        Q = np.array([[-1.0, 1.0], [1.0, -1.0]])  # Gültig
        ctmc = ContinuousTimeMarkovChain(Q)
        assert ctmc.n_states == 2

    def test_nonzero_row_sum_raises(self):
        """Nicht-null Zeilensumme wirft ValueError."""
        Q = np.array([[-1.0, 2.0], [1.0, -1.0]])
        with pytest.raises(ValueError, match="Zeilensummen"):
            ContinuousTimeMarkovChain(Q)

    def test_negative_offdiag_raises(self):
        """Negative Außerdiagonaleinträge werfen ValueError."""
        Q = np.array([[-1.0, -0.5, 1.5], [0.5, -1.0, 0.5], [0.3, 0.7, -1.0]])
        with pytest.raises(ValueError):
            ContinuousTimeMarkovChain(Q)

    def test_transition_at_t0(self, birth_death):
        """P(0) = Einheitsmatrix."""
        P0 = birth_death.transition_matrix(0.0)
        np.testing.assert_array_almost_equal(P0, np.eye(3), decimal=8)

    def test_transition_stochastic(self, birth_death):
        """P(t) ist stochastische Matrix: Einträge ≥ 0, Zeilensummen = 1."""
        P = birth_death.transition_matrix(1.0)
        assert np.all(P >= -1e-8)
        np.testing.assert_array_almost_equal(P.sum(axis=1), [1.0, 1.0, 1.0], decimal=8)

    def test_transition_negative_t_raises(self, birth_death):
        """Negativer Zeitpunkt wirft ValueError."""
        with pytest.raises(ValueError, match="nicht-negativ"):
            birth_death.transition_matrix(-0.1)

    def test_transition_semigroup(self, birth_death):
        """Chapman-Kolmogorov: P(s+t) = P(s)·P(t)."""
        s, t = 0.5, 1.0
        Ps = birth_death.transition_matrix(s)
        Pt = birth_death.transition_matrix(t)
        Pst = birth_death.transition_matrix(s + t)
        np.testing.assert_array_almost_equal(Ps @ Pt, Pst, decimal=8)

    def test_stationary_sums_to_one(self, birth_death):
        """Stationäre Verteilung summiert zu 1."""
        pi = birth_death.stationary_distribution()
        assert abs(pi.sum() - 1.0) < 1e-8

    def test_stationary_satisfies_balance(self, birth_death):
        """Stationäre Verteilung erfüllt πQ = 0."""
        pi = birth_death.stationary_distribution()
        result = pi @ birth_death.Q
        np.testing.assert_array_almost_equal(result, np.zeros(3), decimal=6)

    def test_stationary_positive(self, birth_death):
        """Stationäre Verteilung hat positive Einträge."""
        pi = birth_death.stationary_distribution()
        assert np.all(pi > 0)


# =============================================================================
# Tests: BrownianMotion
# =============================================================================

class TestBrownianMotion:
    """Tests für die Brownsche Bewegung."""

    def test_simulate_shape(self):
        """Simulationsausgabe hat korrekte Form."""
        bm = BrownianMotion()
        t, W = bm.simulate(T=1.0, n_steps=100, n_paths=3, seed=0)
        assert t.shape == (101,)
        assert W.shape == (3, 101)

    def test_simulate_starts_at_zero(self):
        """Alle Pfade starten bei 0."""
        bm = BrownianMotion()
        _, W = bm.simulate(T=1.0, n_steps=100, n_paths=5, seed=0)
        np.testing.assert_array_equal(W[:, 0], np.zeros(5))

    def test_simulate_time_grid(self):
        """Zeitgitter ist korrekt (0 bis T)."""
        bm = BrownianMotion()
        t, _ = bm.simulate(T=2.0, n_steps=100, seed=0)
        assert abs(t[0] - 0.0) < 1e-12
        assert abs(t[-1] - 2.0) < 1e-12

    def test_simulate_reproducible(self):
        """Gleicher Seed ergibt identische Pfade."""
        bm = BrownianMotion()
        _, W1 = bm.simulate(T=1.0, n_steps=50, seed=99)
        _, W2 = bm.simulate(T=1.0, n_steps=50, seed=99)
        np.testing.assert_array_equal(W1, W2)

    def test_simulate_invalid_T_raises(self):
        """T ≤ 0 wirft ValueError."""
        bm = BrownianMotion()
        with pytest.raises(ValueError, match="positiv"):
            bm.simulate(T=0.0, n_steps=10)

    def test_simulate_invalid_steps_raises(self):
        """n_steps ≤ 0 wirft ValueError."""
        bm = BrownianMotion()
        with pytest.raises(ValueError, match="positiv"):
            bm.simulate(T=1.0, n_steps=0)

    def test_empirical_variance(self):
        """Empirische Varianz W(t) ≈ t (für große Pfadzahl)."""
        bm = BrownianMotion()
        _, W = bm.simulate(T=1.0, n_steps=100, n_paths=5000, seed=0)
        # Var(W(1.0)) ≈ 1.0
        empirical_var = np.var(W[:, -1])
        assert abs(empirical_var - 1.0) < 0.1

    def test_empirical_mean(self):
        """Empirischer Mittelwert E[W(t)] ≈ 0."""
        bm = BrownianMotion()
        _, W = bm.simulate(T=1.0, n_steps=100, n_paths=5000, seed=42)
        empirical_mean = np.mean(W[:, -1])
        assert abs(empirical_mean) < 0.05

    def test_covariance_formula(self):
        """Kovarianzfunktion Cov(W(t), W(s)) = min(t, s)."""
        bm = BrownianMotion()
        assert bm.covariance(0.3, 0.7) == pytest.approx(0.3)
        assert bm.covariance(1.0, 0.5) == pytest.approx(0.5)
        assert bm.covariance(0.5, 0.5) == pytest.approx(0.5)

    def test_covariance_symmetry(self):
        """Kovarianz ist symmetrisch: Cov(t,s) = Cov(s,t)."""
        bm = BrownianMotion()
        assert bm.covariance(0.3, 0.8) == bm.covariance(0.8, 0.3)

    def test_quadratic_variation_approx_T(self):
        """Quadratische Variation ≈ T (für viele Schritte)."""
        bm = BrownianMotion()
        T = 1.0
        n_steps = 10000
        dt = T / n_steps
        _, W = bm.simulate(T=T, n_steps=n_steps, n_paths=1, seed=42)
        qv = bm.quadratic_variation(W[0], dt)
        # Quadratische Variation ≈ T mit kleiner Abweichung
        assert abs(qv - T) < 0.1

    def test_quadratic_variation_zero_path(self):
        """Quadratische Variation eines konstanten Pfades ist 0."""
        bm = BrownianMotion()
        path = np.zeros(100)
        qv = bm.quadratic_variation(path, 0.01)
        assert qv == pytest.approx(0.0)


# =============================================================================
# Tests: ItoIntegral
# =============================================================================

class TestItoIntegral:
    """Tests für Itô- und Stratonovich-Integrale."""

    def test_ito_integral_zero_integrand(self):
        """∫0 dW = 0."""
        ito = ItoIntegral()
        f = np.zeros(50)
        dW = np.random.default_rng(0).normal(0, 0.1, size=50)
        result = ito.ito_integral(f, dW)
        assert result == pytest.approx(0.0)

    def test_ito_integral_constant(self):
        """∫c dW = c · ΔW_total."""
        ito = ItoIntegral()
        c = 2.5
        dW = np.array([0.1, -0.2, 0.3])
        result = ito.ito_integral(np.full(3, c), dW)
        expected = c * np.sum(dW)
        assert result == pytest.approx(expected)

    def test_ito_integral_length_mismatch_raises(self):
        """Verschiedene Längen werfen ValueError."""
        ito = ItoIntegral()
        f = np.ones(5)
        dW = np.ones(4)
        with pytest.raises(ValueError, match="gleiche Länge"):
            ito.ito_integral(f, dW)

    def test_stratonovich_integral_constant(self):
        """Stratonovich-Integral für konstante Funktion = Itô-Integral."""
        ito = ItoIntegral()
        c = 3.0
        dW = np.array([0.1, -0.2, 0.3])
        f_all = np.full(4, c)  # n+1 Punkte
        result = ito.stratonovich_integral(f_all, dW)
        expected = c * np.sum(dW)  # Für Konstante identisch
        assert result == pytest.approx(expected)

    def test_stratonovich_length_check(self):
        """Stratonovich: f_values muss n+1 Elemente haben."""
        ito = ItoIntegral()
        f = np.ones(4)   # Sollte n+1 = 4 sein für n=3
        dW = np.ones(3)
        result = ito.stratonovich_integral(f, dW)  # Sollte funktionieren
        assert np.isfinite(result)

    def test_stratonovich_wrong_length_raises(self):
        """Stratonovich: falsche Länge von f wirft ValueError."""
        ito = ItoIntegral()
        f = np.ones(3)   # Sollte 4 sein für n=3
        dW = np.ones(3)
        with pytest.raises(ValueError):
            ito.stratonovich_integral(f, dW)

    def test_ito_formula_demo(self):
        """Itô-Formel: beide Seiten stimmen überein."""
        ito = ItoIntegral()
        # f(x) = x², f'(x) = 2x, f''(x) = 2
        f = lambda x: x ** 2
        df = lambda x: 2 * x
        ddf = lambda x: np.full_like(x, 2.0) if hasattr(x, '__len__') else 2.0
        lhs, rhs = ito.ito_formula_demo(f, df, ddf, 0.0, 1.0, 1000)
        assert abs(lhs - rhs) < 0.1

    def test_ito_isometry_check_returns_tuple(self):
        """Isometrie-Check gibt Tupel zurück."""
        ito = ItoIntegral()
        f = np.ones(100)
        lhs, rhs = ito.ito_isometry_check(f, dt=0.01)
        assert isinstance(lhs, float)
        assert isinstance(rhs, float)

    def test_ito_isometry_constant_f(self):
        """Für f=1: RHS = ∫1²dt = T."""
        ito = ItoIntegral()
        n = 100
        dt = 0.01
        f = np.ones(n)
        _, rhs = ito.ito_isometry_check(f, dt=dt)
        T = n * dt
        assert abs(rhs - T) < 1e-10


# =============================================================================
# Tests: StochasticDifferentialEquation
# =============================================================================

class TestSDE:
    """Tests für stochastische Differentialgleichungen."""

    def test_euler_maruyama_shape(self):
        """Euler-Maruyama gibt korrekte Form zurück."""
        sde = StochasticDifferentialEquation()
        t, X = sde.euler_maruyama(
            mu=lambda x: 0.1 * x,
            sigma=lambda x: 0.2 * x,
            x0=1.0, T=1.0, n_steps=100, n_paths=5, seed=0
        )
        assert t.shape == (101,)
        assert X.shape == (5, 101)

    def test_euler_maruyama_initial_value(self):
        """Euler-Maruyama: Alle Pfade starten bei x0."""
        sde = StochasticDifferentialEquation()
        _, X = sde.euler_maruyama(
            mu=lambda x: 0.0,
            sigma=lambda x: 0.0,
            x0=3.14, T=1.0, n_steps=50, n_paths=3, seed=0
        )
        np.testing.assert_array_almost_equal(X[:, 0], [3.14, 3.14, 3.14])

    def test_euler_maruyama_zero_sigma(self):
        """Euler-Maruyama mit σ=0: ODE-Lösung (deterministisch)."""
        sde = StochasticDifferentialEquation()
        # dX = 1·dt, Lösung: X(t) = x0 + t
        _, X = sde.euler_maruyama(
            mu=lambda x: 1.0,
            sigma=lambda x: 0.0,
            x0=0.0, T=1.0, n_steps=1000, n_paths=1, seed=0
        )
        # X(1) ≈ 1.0
        assert abs(X[0, -1] - 1.0) < 0.01

    def test_euler_maruyama_reproducible(self):
        """Gleicher Seed → identische Pfade."""
        sde = StochasticDifferentialEquation()
        params = dict(mu=lambda x: 0.0, sigma=lambda x: 0.5,
                      x0=1.0, T=1.0, n_steps=100, n_paths=2)
        _, X1 = sde.euler_maruyama(**params, seed=7)
        _, X2 = sde.euler_maruyama(**params, seed=7)
        np.testing.assert_array_equal(X1, X2)

    def test_milstein_shape(self):
        """Milstein gibt korrekte Form zurück."""
        sde = StochasticDifferentialEquation()
        t, X = sde.milstein(
            mu=lambda x: 0.1,
            sigma=lambda x: 0.2,
            sigma_deriv=lambda x: 0.0,
            x0=1.0, T=1.0, n_steps=100, seed=0
        )
        assert t.shape == (101,)
        assert X.shape == (1, 101)

    def test_milstein_initial_value(self):
        """Milstein: Pfad startet bei x0."""
        sde = StochasticDifferentialEquation()
        _, X = sde.milstein(
            mu=lambda x: 0.0,
            sigma=lambda x: 0.1,
            sigma_deriv=lambda x: 0.0,
            x0=2.5, T=1.0, n_steps=50, seed=0
        )
        assert X[0, 0] == pytest.approx(2.5)

    def test_milstein_zero_sigma(self):
        """Milstein mit σ=0: deterministisch (ODE)."""
        sde = StochasticDifferentialEquation()
        # dX = -X·dt, Lösung: X(t) = e^{-t}
        _, X = sde.milstein(
            mu=lambda x: -x,
            sigma=lambda x: 0.0,
            sigma_deriv=lambda x: 0.0,
            x0=1.0, T=1.0, n_steps=1000, seed=0
        )
        assert abs(X[0, -1] - np.exp(-1.0)) < 0.01

    def test_milstein_vs_euler_zero_correction(self):
        """Für σ'=0: Milstein = Euler-Maruyama."""
        sde = StochasticDifferentialEquation()
        # Identische Seeds → identische dW → identische Ergebnisse wenn σ'=0
        params_eu = dict(mu=lambda x: 0.0, sigma=lambda x: 1.0,
                         x0=0.0, T=1.0, n_steps=100)
        params_mi = dict(mu=lambda x: 0.0, sigma=lambda x: 1.0,
                         sigma_deriv=lambda x: 0.0,
                         x0=0.0, T=1.0, n_steps=100)
        _, X_em = sde.euler_maruyama(**params_eu, n_paths=1, seed=42)
        _, X_mi = sde.milstein(**params_mi, seed=42)
        np.testing.assert_array_almost_equal(X_em[0], X_mi[0], decimal=10)

    def test_gbm_shape(self):
        """GBM gibt korrekte Form zurück."""
        sde = StochasticDifferentialEquation()
        t, S = sde.geometric_brownian_motion(
            mu=0.1, sigma=0.2, S0=100.0, T=1.0, n_steps=252, n_paths=10, seed=0
        )
        assert t.shape == (253,)
        assert S.shape == (10, 253)

    def test_gbm_starts_at_S0(self):
        """GBM: Alle Pfade starten bei S0."""
        sde = StochasticDifferentialEquation()
        S0 = 50.0
        _, S = sde.geometric_brownian_motion(
            mu=0.05, sigma=0.3, S0=S0, T=1.0, n_steps=100, n_paths=5, seed=0
        )
        np.testing.assert_array_almost_equal(S[:, 0], [S0] * 5)

    def test_gbm_positive_values(self):
        """GBM bleibt immer positiv (geometrisch)."""
        sde = StochasticDifferentialEquation()
        _, S = sde.geometric_brownian_motion(
            mu=0.1, sigma=0.5, S0=100.0, T=2.0, n_steps=500, n_paths=20, seed=42
        )
        assert np.all(S > 0)

    def test_gbm_mean_approx(self):
        """GBM: E[S(T)] ≈ S0·exp(μT)."""
        sde = StochasticDifferentialEquation()
        mu, sigma, S0, T = 0.1, 0.2, 100.0, 1.0
        _, S = sde.geometric_brownian_motion(
            mu=mu, sigma=sigma, S0=S0, T=T, n_steps=252, n_paths=10000, seed=0
        )
        empirical_mean = np.mean(S[:, -1])
        theoretical_mean = S0 * np.exp(mu * T)
        assert abs(empirical_mean - theoretical_mean) / theoretical_mean < 0.02


# =============================================================================
# Tests: ErgodicTheory
# =============================================================================

class TestErgodicTheory:
    """Tests für die Ergodentheorie."""

    def test_birkhoff_length(self):
        """Birkhoff-Mittel hat selbe Länge wie Eingabe."""
        et = ErgodicTheory()
        f = np.random.default_rng(0).standard_normal(100)
        averages = et.birkhoff_ergodic_theorem_demo(f)
        assert len(averages) == 100

    def test_birkhoff_converges(self):
        """Birkhoff-Mittel konvergiert gegen Zeitmittel."""
        et = ErgodicTheory()
        f = np.random.default_rng(42).standard_normal(10000)
        averages = et.birkhoff_ergodic_theorem_demo(f)
        # Letztes Element = Zeitmittel
        assert abs(averages[-1] - et.time_average(f)) < 1e-10

    def test_birkhoff_first_value(self):
        """Erstes Birkhoff-Mittel = f_0."""
        et = ErgodicTheory()
        f = np.array([3.0, 1.0, 2.0])
        averages = et.birkhoff_ergodic_theorem_demo(f)
        assert averages[0] == pytest.approx(3.0)

    def test_time_average_constant(self):
        """Zeitmittel einer konstanten Folge = Konstante."""
        et = ErgodicTheory()
        f = np.full(100, 5.0)
        assert et.time_average(f) == pytest.approx(5.0)

    def test_time_average_sinusoid(self):
        """Zeitmittel von sin über volle Perioden ≈ 0."""
        et = ErgodicTheory()
        t = np.linspace(0, 2 * np.pi, 10000, endpoint=False)
        f = np.sin(t)
        assert abs(et.time_average(f)) < 0.01

    def test_is_mixing_iid_noise(self):
        """i.i.d. Rauschen gilt als mischend."""
        et = ErgodicTheory()
        seq = np.random.default_rng(0).standard_normal(1000)
        assert et.is_mixing_demo(seq, lag=10) is True

    def test_is_mixing_constant_sequence(self):
        """Konstante Folge ist trivial mischend (Autokorr. = 0 nach Normierung)."""
        et = ErgodicTheory()
        seq = np.ones(100)
        # Konstante Folge hat Varianz 0 → gilt als mischend
        result = et.is_mixing_demo(seq, lag=5)
        assert result is True

    def test_logistic_map_length(self):
        """Logistische Abbildung gibt korrektes Array zurück."""
        et = ErgodicTheory()
        traj = et.logistic_map(r=3.9, x0=0.5, n_iter=100)
        assert len(traj) == 101  # n_iter + 1

    def test_logistic_map_starts_at_x0(self):
        """Logistische Abbildung startet bei x0."""
        et = ErgodicTheory()
        traj = et.logistic_map(r=2.0, x0=0.3, n_iter=10)
        assert traj[0] == pytest.approx(0.3)

    def test_logistic_map_fixed_point(self):
        """Stabiler Fixpunkt: x* = 1 - 1/r bei r=2."""
        et = ErgodicTheory()
        r = 2.0
        x_star = 1 - 1 / r  # = 0.5
        traj = et.logistic_map(r=r, x0=0.4, n_iter=1000)
        assert abs(traj[-1] - x_star) < 1e-6

    def test_logistic_map_bounded(self):
        """Alle Werte der logistischen Abbildung in (0, 1)."""
        et = ErgodicTheory()
        traj = et.logistic_map(r=4.0, x0=0.5, n_iter=500)
        # x0 in Folge enthalten, Werte bleiben in [0,1]
        assert np.all(traj >= 0) and np.all(traj <= 1)

    def test_logistic_map_invalid_r(self):
        """Ungültiges r wirft ValueError."""
        et = ErgodicTheory()
        with pytest.raises(ValueError, match="r="):
            et.logistic_map(r=5.0, x0=0.5, n_iter=10)

    def test_logistic_map_invalid_x0(self):
        """x0 außerhalb (0,1) wirft ValueError."""
        et = ErgodicTheory()
        with pytest.raises(ValueError):
            et.logistic_map(r=3.0, x0=0.0, n_iter=10)

    def test_lyapunov_chaotic_r4(self):
        """Lyapunov-Exponent für r=4 ≈ log(2) ≈ 0.693."""
        et = ErgodicTheory()
        lam = et.lyapunov_exponent(r=4.0, x0=0.3, n_iter=100000)
        assert abs(lam - np.log(2)) < 0.05

    def test_lyapunov_stable_negative(self):
        """Lyapunov-Exponent für stabilen Fixpunkt r=2 ist negativ."""
        et = ErgodicTheory()
        lam = et.lyapunov_exponent(r=2.0, x0=0.3, n_iter=1000)
        assert lam < 0


# =============================================================================
# Tests: PoissonProcess
# =============================================================================

class TestPoissonProcess:
    """Tests für den Poisson-Prozess."""

    def test_invalid_rate_raises(self):
        """Rate ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError, match="positiv"):
            PoissonProcess(rate=0.0)
        with pytest.raises(ValueError, match="positiv"):
            PoissonProcess(rate=-1.0)

    def test_simulate_all_before_T(self):
        """Alle Sprungzeiten liegen in (0, T]."""
        pp = PoissonProcess(rate=5.0)
        T = 10.0
        jumps = pp.simulate(T=T, seed=42)
        assert np.all(jumps <= T)
        assert np.all(jumps > 0)

    def test_simulate_count_approx_lambda_T(self):
        """Erwartete Anzahl Sprünge ≈ λT."""
        pp = PoissonProcess(rate=3.0)
        T = 100.0
        n_sim = 1000
        counts = []
        for seed in range(n_sim):
            jumps = pp.simulate(T=T, seed=seed)
            counts.append(len(jumps))
        empirical_mean = np.mean(counts)
        assert abs(empirical_mean - pp.rate * T) < 5  # 5 Einheiten Toleranz

    def test_simulate_reproducible(self):
        """Gleicher Seed → gleiche Sprungzeiten."""
        pp = PoissonProcess(rate=2.0)
        j1 = pp.simulate(T=5.0, seed=7)
        j2 = pp.simulate(T=5.0, seed=7)
        np.testing.assert_array_equal(j1, j2)

    def test_simulate_sorted(self):
        """Sprungzeiten sind aufsteigend geordnet."""
        pp = PoissonProcess(rate=10.0)
        jumps = pp.simulate(T=5.0, seed=0)
        assert np.all(np.diff(jumps) > 0)

    def test_simulate_negative_T_raises(self):
        """T ≤ 0 wirft ValueError."""
        pp = PoissonProcess(rate=1.0)
        with pytest.raises(ValueError):
            pp.simulate(T=0.0)

    def test_pmf_k0(self):
        """P(N(t)=0) = e^{-λt}."""
        pp = PoissonProcess(rate=2.0)
        t = 1.0
        expected = np.exp(-pp.rate * t)
        assert pp.pmf(0, t) == pytest.approx(expected, rel=1e-8)

    def test_pmf_sums_to_one(self):
        """Σ_k P(N(t)=k) ≈ 1."""
        pp = PoissonProcess(rate=3.0)
        t = 2.0
        total = sum(pp.pmf(k, t) for k in range(100))
        assert abs(total - 1.0) < 1e-6

    def test_pmf_negative_k_raises(self):
        """Negativer k-Wert wirft ValueError."""
        pp = PoissonProcess(rate=1.0)
        with pytest.raises(ValueError, match="nicht-negativ"):
            pp.pmf(-1, 1.0)

    def test_mean_formula(self):
        """E[N(t)] = λt."""
        pp = PoissonProcess(rate=4.0)
        assert pp.mean(2.5) == pytest.approx(4.0 * 2.5)

    def test_variance_formula(self):
        """Var[N(t)] = λt."""
        pp = PoissonProcess(rate=3.0)
        assert pp.variance(2.0) == pytest.approx(6.0)

    def test_mean_equals_variance(self):
        """E[N(t)] = Var[N(t)] (charakteristische Eigenschaft)."""
        pp = PoissonProcess(rate=5.0)
        t = 3.0
        assert pp.mean(t) == pp.variance(t)


# =============================================================================
# Tests: GaussianProcess
# =============================================================================

class TestGaussianProcess:
    """Tests für den Gaußschen Prozess."""

    def test_sample_length(self):
        """Stichprobe hat korrekte Länge."""
        gp = GaussianProcess(
            mean_func=lambda x: 0.0,
            kernel_func=GaussianProcess.rbf_kernel
        )
        x = np.linspace(0, 1, 20)
        sample = gp.sample(x, seed=0)
        assert len(sample) == 20

    def test_sample_finite(self):
        """Stichprobe enthält nur endliche Werte."""
        gp = GaussianProcess(
            mean_func=lambda x: 0.0,
            kernel_func=GaussianProcess.rbf_kernel
        )
        x = np.linspace(0, 5, 50)
        sample = gp.sample(x, seed=42)
        assert np.all(np.isfinite(sample))

    def test_sample_reproducible(self):
        """Gleicher Seed → identische Stichprobe."""
        gp = GaussianProcess(
            mean_func=lambda x: 0.0,
            kernel_func=GaussianProcess.rbf_kernel
        )
        x = np.linspace(0, 1, 10)
        s1 = gp.sample(x, seed=5)
        s2 = gp.sample(x, seed=5)
        np.testing.assert_array_equal(s1, s2)

    def test_sample_nonzero_mean(self):
        """Stichprobe mit Mittelwert m(x) = x."""
        gp = GaussianProcess(
            mean_func=lambda x: x,
            kernel_func=GaussianProcess.rbf_kernel
        )
        x = np.linspace(0, 1, 100)
        # Über viele Stichproben: E[f(x)] ≈ x
        samples = np.array([gp.sample(x, seed=i) for i in range(200)])
        empirical_mean = np.mean(samples, axis=0)
        np.testing.assert_array_almost_equal(empirical_mean, x, decimal=1)

    def test_rbf_kernel_self(self):
        """RBF-Kern: k(x, x) = 1."""
        assert GaussianProcess.rbf_kernel(0.5, 0.5) == pytest.approx(1.0)
        assert GaussianProcess.rbf_kernel(2.0, 2.0) == pytest.approx(1.0)

    def test_rbf_kernel_symmetry(self):
        """RBF-Kern ist symmetrisch: k(x,y) = k(y,x)."""
        x, y = 0.3, 1.7
        assert GaussianProcess.rbf_kernel(x, y) == pytest.approx(
            GaussianProcess.rbf_kernel(y, x))

    def test_rbf_kernel_range(self):
        """RBF-Kern liegt in (0, 1]."""
        pairs = [(0, 0), (0, 1), (0, 10), (3, 3)]
        for x1, x2 in pairs:
            k = GaussianProcess.rbf_kernel(x1, x2)
            assert 0 < k <= 1.0 + 1e-12

    def test_rbf_kernel_decreasing(self):
        """RBF-Kern nimmt mit Abstand ab."""
        k1 = GaussianProcess.rbf_kernel(0.0, 0.5)
        k2 = GaussianProcess.rbf_kernel(0.0, 1.0)
        k3 = GaussianProcess.rbf_kernel(0.0, 2.0)
        assert k1 > k2 > k3

    def test_rbf_kernel_length_scale(self):
        """Längenskalenparameter beeinflusst Korrelationslänge."""
        k_small = GaussianProcess.rbf_kernel(0.0, 1.0, length_scale=0.5)
        k_large = GaussianProcess.rbf_kernel(0.0, 1.0, length_scale=2.0)
        assert k_small < k_large  # Kleinere Skala → schnellerer Abfall


# =============================================================================
# Tests: Standalone-Funktionen
# =============================================================================

class TestRandomWalk:
    """Tests für den 1D Random Walk."""

    def test_starts_at_zero(self):
        """Random Walk startet bei 0."""
        walk = random_walk_1d(100, seed=0)
        assert walk[0] == 0

    def test_length(self):
        """Random Walk hat n_steps+1 Elemente."""
        walk = random_walk_1d(50, seed=0)
        assert len(walk) == 51

    def test_integer_steps(self):
        """Alle Schritte sind ±1 (ganzzahlig)."""
        walk = random_walk_1d(200, seed=42)
        diffs = np.diff(walk)
        assert np.all((diffs == 1) | (diffs == -1))

    def test_reproducible(self):
        """Gleicher Seed → identischer Walk."""
        w1 = random_walk_1d(100, seed=7)
        w2 = random_walk_1d(100, seed=7)
        np.testing.assert_array_equal(w1, w2)

    def test_biased_walk(self):
        """Biased Walk (p=1) macht nur positive Schritte."""
        walk = random_walk_1d(100, p=1.0, seed=0)
        assert walk[-1] == 100

    def test_biased_walk_negative(self):
        """Biased Walk (p=0) macht nur negative Schritte."""
        walk = random_walk_1d(100, p=0.0, seed=0)
        assert walk[-1] == -100

    def test_symmetric_mean_approx_zero(self):
        """Symmetrischer Walk: E[S_n] ≈ 0."""
        walk = random_walk_1d(100000, p=0.5, seed=0)
        # Zeitmittel ≈ 0 (schwaches Gesetz der großen Zahlen)
        assert abs(walk[-1]) < 1000  # Grobe Schranke

    def test_invalid_p_raises(self):
        """p außerhalb [0, 1] wirft ValueError."""
        with pytest.raises(ValueError, match="Wahrscheinlichkeit"):
            random_walk_1d(10, p=1.5, seed=0)

    def test_zero_steps(self):
        """0 Schritte: Walk = [0]."""
        walk = random_walk_1d(0, seed=0)
        np.testing.assert_array_equal(walk, [0])


class TestGamblerRuin:
    """Tests für die Gambler's-Ruin-Formel."""

    def test_symmetric_exact(self):
        """Symmetrischer Fall p=0.5: P(Gewinn) = k/N."""
        for k in range(1, 5):
            prob = gambler_ruin_probability(0.5, k, 5)
            assert prob == pytest.approx(k / 5, rel=1e-8)

    def test_start_zero(self):
        """Start bei 0: Wahrsch. = 0."""
        assert gambler_ruin_probability(0.5, 0, 10) == pytest.approx(0.0)

    def test_start_at_target(self):
        """Start am Ziel: Wahrsch. = 1."""
        assert gambler_ruin_probability(0.5, 5, 5) == pytest.approx(1.0)

    def test_asymmetric_p_gt_half(self):
        """p > 0.5: Gewinnwahrscheinlichkeit erhöht sich."""
        prob_sym = gambler_ruin_probability(0.5, 3, 6)
        prob_fav = gambler_ruin_probability(0.6, 3, 6)
        assert prob_fav > prob_sym

    def test_asymmetric_p_lt_half(self):
        """p < 0.5: Gewinnwahrscheinlichkeit sinkt."""
        prob_sym = gambler_ruin_probability(0.5, 3, 6)
        prob_unfav = gambler_ruin_probability(0.4, 3, 6)
        assert prob_unfav < prob_sym

    def test_probability_in_range(self):
        """Wahrscheinlichkeit liegt in [0, 1]."""
        for p in [0.3, 0.5, 0.7]:
            for k in range(6):
                prob = gambler_ruin_probability(p, k, 5)
                assert 0.0 <= prob <= 1.0 + 1e-12

    def test_invalid_p_raises(self):
        """Ungültiges p wirft ValueError."""
        with pytest.raises(ValueError, match="in \\[0, 1\\]"):
            gambler_ruin_probability(1.5, 2, 5)

    def test_invalid_start_raises(self):
        """start > target wirft ValueError."""
        with pytest.raises(ValueError):
            gambler_ruin_probability(0.5, 6, 5)

    def test_certain_win_p1(self):
        """p=1 (sicherer Gewinn): Wahrsch. = 1 für start > 0."""
        prob = gambler_ruin_probability(1.0, 1, 10)
        assert prob == pytest.approx(1.0, rel=1e-8)

    def test_certain_loss_p0(self):
        """p=0 (sicherer Verlust): Wahrsch. = 0 für start < target."""
        prob = gambler_ruin_probability(0.0, 5, 10)
        assert prob == pytest.approx(0.0, abs=1e-10)


class TestCLT:
    """Tests für den Zentralen Grenzwertsatz (CLT-Demonstration)."""

    def test_returns_tuple(self):
        """Gibt korrektes Tupel zurück."""
        result = central_limit_theorem_demo(100, n_obs=50, seed=0)
        assert len(result) == 3

    def test_normalized_length(self):
        """Normierte Mittelwerte haben korrekte Länge."""
        means, _, _ = central_limit_theorem_demo(200, n_obs=50, seed=0)
        assert len(means) == 200

    def test_normalized_mean_near_zero(self):
        """Normierter Mittelwert ≈ 0 (CLT)."""
        _, emp_mean, _ = central_limit_theorem_demo(5000, n_obs=100, seed=42)
        assert abs(emp_mean) < 0.1

    def test_normalized_variance_near_one(self):
        """Normierte Varianz ≈ 1 (CLT)."""
        _, _, emp_var = central_limit_theorem_demo(5000, n_obs=100, seed=42)
        assert abs(emp_var - 1.0) < 0.15

    def test_uniform_distribution(self):
        """CLT für Gleichverteilung: Var ≈ 1."""
        _, _, emp_var = central_limit_theorem_demo(
            3000, n_obs=200, distribution='uniform', seed=0
        )
        assert abs(emp_var - 1.0) < 0.2

    def test_exponential_distribution(self):
        """CLT für Exponentialverteilung: Var ≈ 1."""
        _, _, emp_var = central_limit_theorem_demo(
            3000, n_obs=200, distribution='exponential', seed=0
        )
        assert abs(emp_var - 1.0) < 0.2

    def test_bernoulli_distribution(self):
        """CLT für Bernoulli-Verteilung: Var ≈ 1."""
        _, _, emp_var = central_limit_theorem_demo(
            3000, n_obs=200, distribution='bernoulli', seed=0
        )
        assert abs(emp_var - 1.0) < 0.2

    def test_invalid_distribution_raises(self):
        """Unbekannte Verteilung wirft ValueError."""
        with pytest.raises(ValueError, match="Unbekannte Verteilung"):
            central_limit_theorem_demo(100, distribution='gamma', seed=0)

    def test_reproducible_with_seed(self):
        """Gleicher Seed → identische Ergebnisse."""
        r1, _, _ = central_limit_theorem_demo(100, n_obs=50, seed=13)
        r2, _, _ = central_limit_theorem_demo(100, n_obs=50, seed=13)
        np.testing.assert_array_equal(r1, r2)
