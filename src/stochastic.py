"""
@file stochastic.py
@brief Stochastik-Modul: Markov-Ketten, Brownsche Bewegung, Itô-Kalkül, SDEs,
       Ergodentheorie, Poisson-Prozess und Gaußsche Prozesse.
@description
    Dieses Modul implementiert grundlegende und fortgeschrittene stochastische
    Prozesse und Methoden:

    - MarkovChain: Diskrete Markov-Kette mit stationärer Verteilung, Simulation,
      mittlerer Übergangszeit und Absorptionswahrscheinlichkeiten
    - ContinuousTimeMarkovChain: Zeitkontinuierliche Markov-Kette via Matrixexponential
    - BrownianMotion: Simulation, Kovarianz, quadratische Variation
    - ItoIntegral: Itô- und Stratonovich-Integrale, Itô-Formel, Itô-Isometrie
    - StochasticDifferentialEquation: Euler-Maruyama, Milstein, geometrische BM
    - ErgodicTheory: Birkhoff-Ergodentheorem, Lyapunov-Exponent, logistische Abbildung
    - PoissonProcess: Simulation, PMF, Momente
    - GaussianProcess: Stichproben, RBF-Kern
    - Standalone-Funktionen: Random Walk, Gambler's Ruin, CLT-Demonstration

    Mathematische Grundlagen:
    - Übergangswahrscheinlichkeit: P(X_{n+1}=j | X_n=i) = P_{ij}
    - Stationäre Verteilung: πP = π, π·1 = 1
    - Brownsche Bewegung: W(t) ~ N(0, t), Cov(W(t),W(s)) = min(t,s)
    - Itô-Integral: ∫₀ᵀ f dW = lim Σ f(t_k)[W(t_{k+1}) - W(t_k)]
    - SDE: dX = μ(X)dt + σ(X)dW (Itô-Konvention)
    - Geometrische BM: S(t) = S₀ exp((μ - σ²/2)t + σW(t))

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
import scipy.linalg as la
from typing import Callable, Optional, Tuple, List


# =============================================================================
# Klasse: MarkovChain
# =============================================================================

class MarkovChain:
    """
    @class MarkovChain
    @brief Diskrete zeitliche Markov-Kette mit endlichem Zustandsraum.

    @description
        Modelliert einen stochastischen Prozess (X_n)_{n≥0} auf einem
        endlichen Zustandsraum S = {0, 1, ..., N-1}, wobei die Übergangs-
        wahrscheinlichkeiten durch die Matrix P gegeben sind:

          P_{ij} = P(X_{n+1} = j | X_n = i)

        Eigenschaften:
        - Jede Zeile von P summiert zu 1 (stochastische Matrix)
        - Stationäre Verteilung: π mit πP = π, Σπ_i = 1
        - Ergodisch ⟺ irreduzibel und aperiodisch

    @param transition_matrix Quadratische stochastische Übergangsmatrix (N×N)
    @lastModified 2026-03-10
    """

    def __init__(self, transition_matrix: np.ndarray) -> None:
        """
        @brief Initialisiert die Markov-Kette mit einer Übergangsmatrix.

        @param transition_matrix Stochastische N×N-Matrix; Zeilensummen ≈ 1
        @raises ValueError Wenn Matrix nicht quadratisch oder Zeilen nicht summieren
        @lastModified 2026-03-10
        """
        P = np.array(transition_matrix, dtype=float)

        # Überprüfung: quadratische Matrix
        if P.ndim != 2 or P.shape[0] != P.shape[1]:
            raise ValueError("Übergangsmatrix muss quadratisch sein.")

        # Überprüfung: alle Einträge nicht-negativ
        if np.any(P < -1e-10):
            raise ValueError("Übergangsmatrix darf keine negativen Einträge haben.")

        # Überprüfung: Zeilensummen ≈ 1 (stochastische Matrix)
        row_sums = P.sum(axis=1)
        if not np.allclose(row_sums, 1.0, atol=1e-8):
            raise ValueError(
                f"Zeilensummen müssen 1 sein, gefunden: {row_sums}"
            )

        self.P = P
        self.n_states = P.shape[0]

    def stationary_distribution(self) -> np.ndarray:
        """
        @brief Berechnet die stationäre Verteilung π via Eigenvektormethode.

        @description
            Sucht den linken Eigenvektor zum Eigenwert 1 der Übergangsmatrix.
            Die stationäre Verteilung π erfüllt:
              π = πP  ⟺  πᵀ = Pᵀ πᵀ (rechter Eigenvektor von Pᵀ)

            Vorgehen:
            1. Berechne Eigenwerte und -vektoren von Pᵀ
            2. Wähle Eigenvektor zum Eigenwert ≈ 1
            3. Normiere auf Summe 1

        @return Stationäre Verteilung als 1D-Array der Länge N
        @raises ValueError Wenn kein Eigenvektor zum Eigenwert 1 gefunden
        @lastModified 2026-03-10
        """
        # Eigenwerte und -vektoren der transponierten Matrix (linke EV)
        eigenvalues, eigenvectors = np.linalg.eig(self.P.T)

        # Suche Eigenwert nahe 1 (stochastische Matrizen haben immer λ=1)
        idx = np.argmin(np.abs(eigenvalues - 1.0))

        if np.abs(eigenvalues[idx] - 1.0) > 1e-6:
            raise ValueError("Kein Eigenwert 1 gefunden – keine stationäre Verteilung.")

        # Realteil nehmen (numerische Imaginarteile unterdrücken)
        pi = np.real(eigenvectors[:, idx])

        # Sicherstellen, dass alle Einträge positiv sind
        if np.all(pi <= 0):
            pi = -pi

        # Normieren: Summe = 1
        pi = pi / pi.sum()
        return pi

    def simulate(
        self,
        initial_state: int,
        n_steps: int,
        seed: Optional[int] = None
    ) -> np.ndarray:
        """
        @brief Simuliert einen Pfad der Markov-Kette.

        @description
            Erzeugt eine Folge (X_0, X_1, ..., X_{n_steps}) durch zufällige
            Übergänge gemäß der Übergangsmatrix P.

            Algorithmus: In jedem Schritt wird der nächste Zustand aus der
            Wahrscheinlichkeitsverteilung P[current, :] gezogen.

        @param initial_state Startzustand X_0 ∈ {0, ..., N-1}
        @param n_steps Anzahl Zeitschritte
        @param seed Zufallssaat für Reproduzierbarkeit
        @return Array der Zustände [X_0, X_1, ..., X_{n_steps}]
        @raises ValueError Wenn initial_state außerhalb des Zustandsraums
        @lastModified 2026-03-10
        """
        if not (0 <= initial_state < self.n_states):
            raise ValueError(
                f"initial_state={initial_state} außerhalb [0, {self.n_states-1}]."
            )

        rng = np.random.default_rng(seed)
        states = np.zeros(n_steps + 1, dtype=int)
        states[0] = initial_state

        for t in range(n_steps):
            current = states[t]
            # Ziehe nächsten Zustand gemäß Übergangswahrscheinlichkeiten
            states[t + 1] = rng.choice(self.n_states, p=self.P[current])

        return states

    def n_step_distribution(
        self,
        initial_dist: np.ndarray,
        n: int
    ) -> np.ndarray:
        """
        @brief Berechnet die n-Schritt-Verteilung P^n · π₀.

        @description
            Für eine Anfangsverteilung π₀ (Vektor) gilt nach n Schritten:
              π_n = π₀ · P^n

            Berechnung durch wiederholte Matrixmultiplikation oder
            schnelle Potenzierung.

        @param initial_dist Anfangsverteilung π₀ als 1D-Array (Summe = 1)
        @param n Anzahl der Schritte
        @return n-Schritt-Verteilung als 1D-Array
        @lastModified 2026-03-10
        """
        dist = np.array(initial_dist, dtype=float)

        # Normierung sicherstellen
        dist = dist / dist.sum()

        # P^n via numpy matrix power
        Pn = np.linalg.matrix_power(self.P, n)
        return dist @ Pn

    def is_ergodic(self) -> bool:
        """
        @brief Prüft, ob die Markov-Kette ergodisch ist.

        @description
            Eine Markov-Kette ist ergodisch, wenn sie:
            1. Irreduzibel ist: Jeder Zustand ist von jedem anderen erreichbar
               (zusammenhängender gerichteter Graph der Übergangsstruktur)
            2. Aperiodisch ist: Der ggT der Rückkehrzeiten jedes Zustands ist 1

            Überprüfung der Irreduzibilität via erreichbare Zustände (BFS/DFS):
            Berechne (I + P)^N > 0 (alle Einträge positiv).

            Überprüfung der Aperiodizität: Eine irreduzible Kette ist aperiodisch,
            wenn mindestens ein Zustand einen Selbstübergang P_{ii} > 0 hat,
            oder via Bestimmung der Periode durch Pfadlängenanalyse.

        @return True wenn die Kette ergodisch (irreduzibel + aperiodisch) ist
        @lastModified 2026-03-10
        """
        N = self.n_states

        # ── Irreduzibilität: (I+P)^N hat nur positive Einträge ──────────────
        # Wenn (I+P)^N_{ij} > 0 für alle i,j, dann ist j von i in ≤N Schritten
        # erreichbar. Dies prüft die starke Zusammenhangskomponente.
        A = np.eye(N) + self.P
        An = np.linalg.matrix_power(A, N)
        if not np.all(An > 1e-10):
            return False  # Nicht irreduzibel

        # ── Aperiodizität: Suche Zustand mit positivem Selbstübergang ───────
        # Einfachste Methode: Falls P_{ii} > 0 für irgendein i, ist die Kette
        # aperiodisch (Periode 1 für diesen Zustand, und bei Irreduzibilität
        # haben alle Zustände Periode 1).
        if np.any(np.diag(self.P) > 1e-10):
            return True

        # Aufwändigere Methode: Berechne Periode über Pfadlängen
        # Für jeden Zustand i: ggT aller n mit P^n_{ii} > 0
        # Approximation: Teste P^2, P^3, ..., P^{2N}
        period = N  # Obere Schranke
        Pk = self.P.copy()
        for k in range(2, 2 * N + 1):
            Pk = Pk @ self.P
            # Prüfe ob Diagonale von P^k positive Einträge hat
            for i in range(N):
                if Pk[i, i] > 1e-10:
                    import math
                    period = math.gcd(period, k)
        return period == 1

    def mean_first_passage_time(self, i: int, j: int) -> float:
        """
        @brief Berechnet die mittlere erstmalige Übergangszeit m_{ij}.

        @description
            Die mittlere erstmalige Übergangszeit von Zustand i nach j ist:
              m_{ij} = E[min{n ≥ 1 : X_n = j} | X_0 = i]

            Für i ≠ j berechnen wir das lineare Gleichungssystem:
              m_{ij} = 1 + Σ_{k≠j} P_{ik} · m_{kj}

            Für i = j (mittlere Rückkehrzeit):
              m_{ii} = 1 / π_i

        @param i Startzustand
        @param j Zielzustand
        @return Mittlere erstmalige Übergangszeit m_{ij}
        @lastModified 2026-03-10
        """
        N = self.n_states

        # Mittlere Rückkehrzeit: 1/π_j
        if i == j:
            pi = self.stationary_distribution()
            if pi[j] < 1e-14:
                return np.inf
            return 1.0 / pi[j]

        # Mittlere erstmalige Übergangszeit für i ≠ j:
        # Löse: m_k = 1 + Σ_{l≠j} P_{kl} · m_l  für k ≠ j, m_j = 0
        # Umformung: m - P_{-j,-j} · m = 1  (submatrix ohne Spalte/Zeile j)

        # Zustandsliste ohne Zielzustand j
        transient = [k for k in range(N) if k != j]
        n_t = len(transient)

        # Teilübergangsmatrix (ohne absorbierende Zeile/Spalte j)
        Q = self.P[np.ix_(transient, transient)]

        # Gleichungssystem: (I - Q) · m = 1
        A = np.eye(n_t) - Q
        b = np.ones(n_t)

        try:
            m_vec = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.inf

        # m_{ij} ist der Eintrag für Startzustand i
        idx = transient.index(i)
        return float(m_vec[idx])

    def absorption_probabilities(
        self,
        transient_states: List[int],
        absorbing_states: List[int]
    ) -> np.ndarray:
        """
        @brief Berechnet die Absorptionswahrscheinlichkeiten für absorbierende Ketten.

        @description
            In einer absorbierenden Markov-Kette mit transienten Zuständen T
            und absorbierenden Zuständen A gilt:
              B = (I - Q)^{-1} · R

            wobei:
            - Q = Teilmatrix der transient→transient Übergänge
            - R = Teilmatrix der transient→absorbing Übergänge
            - N = (I - Q)^{-1} = fundamentale Matrix
            - B = Absorptionswahrscheinlichkeiten (|T| × |A|-Matrix)

            B[i, j] = Wahrscheinlichkeit, von transientem Zustand i in
                      absorbierendem Zustand j absorbiert zu werden.

        @param transient_states Liste der transienten Zustände
        @param absorbing_states Liste der absorbierenden Zustände
        @return Matrix B der Absorptionswahrscheinlichkeiten (|T| × |A|)
        @lastModified 2026-03-10
        """
        t = transient_states
        a = absorbing_states

        # Teilübergangsmatrix transient → transient
        Q = self.P[np.ix_(t, t)]

        # Teilübergangsmatrix transient → absorbierend
        R = self.P[np.ix_(t, a)]

        # Fundamentale Matrix N = (I - Q)^{-1}
        I = np.eye(len(t))
        try:
            N_fund = np.linalg.inv(I - Q)
        except np.linalg.LinAlgError:
            raise ValueError("(I-Q) ist singulär – Prüfe ob Zustände wirklich transient sind.")

        # Absorptionswahrscheinlichkeiten B = N · R
        B = N_fund @ R
        return B


# =============================================================================
# Klasse: ContinuousTimeMarkovChain
# =============================================================================

class ContinuousTimeMarkovChain:
    """
    @class ContinuousTimeMarkovChain
    @brief Zeitkontinuierliche Markov-Kette (CTMC) mit Generatormatrix Q.

    @description
        Eine zeitkontinuierliche Markov-Kette wird durch ihre Generatormatrix
        (auch Ratenmatrix oder Q-Matrix) beschrieben:

          Q_{ij} ≥ 0  für i ≠ j  (Übergangsraten)
          Q_{ii} = -Σ_{j≠i} Q_{ij}  (Zeilensumme = 0)

        Die Übergangswahrscheinlichkeitsmatrix nach Zeit t ist:
          P(t) = exp(Qt)  (Matrixexponential)

        Stationäre Verteilung π: πQ = 0, π·1 = 1

    @param rate_matrix Generatormatrix Q (N×N, Zeilensummen = 0)
    @lastModified 2026-03-10
    """

    def __init__(self, rate_matrix: np.ndarray) -> None:
        """
        @brief Initialisiert die CTMC mit Generatormatrix Q.

        @param rate_matrix N×N-Matrix mit nicht-negativen Außerdiagonaleinträgen
                           und Zeilensummen = 0
        @raises ValueError Wenn Matrix ungültig
        @lastModified 2026-03-10
        """
        Q = np.array(rate_matrix, dtype=float)

        if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
            raise ValueError("Generatormatrix muss quadratisch sein.")

        # Überprüfung: Außerdiagonaleinträge ≥ 0
        N = Q.shape[0]
        for i in range(N):
            for j in range(N):
                if i != j and Q[i, j] < -1e-10:
                    raise ValueError(
                        f"Außerdiagonaleintrag Q[{i},{j}]={Q[i,j]} muss ≥ 0 sein."
                    )

        # Überprüfung: Zeilensummen ≈ 0
        row_sums = Q.sum(axis=1)
        if not np.allclose(row_sums, 0.0, atol=1e-8):
            raise ValueError(
                f"Zeilensummen der Generatormatrix müssen 0 sein, gefunden: {row_sums}"
            )

        self.Q = Q
        self.n_states = N

    def transition_matrix(self, t: float) -> np.ndarray:
        """
        @brief Berechnet die Übergangsmatrix P(t) = exp(Qt).

        @description
            Das Matrixexponential wird via scipy.linalg.expm berechnet,
            welches intern Padé-Approximation mit Skalierung verwendet.

            P(t) ist eine stochastische Matrix: P(t)_{ij} ≥ 0, Zeilensummen = 1.

        @param t Zeitpunkt t ≥ 0
        @return N×N Übergangsmatrix P(t)
        @raises ValueError Wenn t < 0
        @lastModified 2026-03-10
        """
        if t < 0:
            raise ValueError(f"Zeitpunkt t={t} muss nicht-negativ sein.")

        # Matrixexponential exp(Q·t)
        return la.expm(self.Q * t)

    def stationary_distribution(self) -> np.ndarray:
        """
        @brief Berechnet die stationäre Verteilung der CTMC.

        @description
            Die stationäre Verteilung π erfüllt: πQ = 0, Σπ_i = 1
            Dies ist äquivalent zu: Qᵀ π = 0

            Lösung: Linker Nullvektor von Q, normiert auf Summe 1.
            Berechnet via SVD des transponierten Systems.

        @return Stationäre Verteilung π als 1D-Array
        @lastModified 2026-03-10
        """
        # Nullvektor von Qᵀ via SVD (letzter rechter Singulärvektor)
        _, _, Vt = np.linalg.svd(self.Q.T)
        pi = np.abs(Vt[-1])  # Letzter Zeile von Vᵀ = Nullvektor

        # Normierung auf Summe 1
        pi = pi / pi.sum()
        return pi


# =============================================================================
# Klasse: BrownianMotion
# =============================================================================

class BrownianMotion:
    """
    @class BrownianMotion
    @brief Standard-Brownsche Bewegung (Wiener-Prozess) W(t).

    @description
        Die Brownsche Bewegung (oder Wiener-Prozess) ist ein stochastischer
        Prozess (W(t))_{t≥0} mit:
        - W(0) = 0 fast sicher
        - Unabhängige Inkremente: W(t) - W(s) ⊥ W(s) - W(r) für r ≤ s ≤ t
        - Normalverteilte Inkremente: W(t) - W(s) ~ N(0, t-s)
        - Stetige Pfade (fast sicher)

        Kovarianzfunktion: Cov(W(t), W(s)) = min(t, s)
        Quadratische Variation: [W]_t = t (fast sicher)

    @lastModified 2026-03-10
    """

    def simulate(
        self,
        T: float,
        n_steps: int,
        n_paths: int = 1,
        seed: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Simuliert Pfade der Brownschen Bewegung (Euler-Maruyama-Diskretisierung).

        @description
            Diskretisierung: W(t_{k+1}) = W(t_k) + √(Δt) · Z_k
            wobei Z_k ~ N(0,1) i.i.d.

            Dies entspricht dem Euler-Maruyama-Verfahren für dW = dW
            (triviale SDE ohne Drift und Diffusion = 1).

        @param T Endzeitpunkt (T > 0)
        @param n_steps Anzahl Zeitschritte (Diskretisierung)
        @param n_paths Anzahl unabhängiger Pfade (Standard: 1)
        @param seed Zufallssaat für Reproduzierbarkeit
        @return Tupel (t, W): Zeitgitter (n_steps+1,) und Pfade (n_paths, n_steps+1)
        @raises ValueError Wenn T ≤ 0 oder n_steps ≤ 0
        @lastModified 2026-03-10
        """
        if T <= 0:
            raise ValueError(f"Endzeitpunkt T={T} muss positiv sein.")
        if n_steps <= 0:
            raise ValueError(f"n_steps={n_steps} muss positiv sein.")

        dt = T / n_steps
        t = np.linspace(0, T, n_steps + 1)

        rng = np.random.default_rng(seed)

        # Inkremente: Z_k ~ N(0, √dt) für jeden Pfad
        # Form: (n_paths, n_steps)
        increments = rng.normal(0, np.sqrt(dt), size=(n_paths, n_steps))

        # Kumulierte Summe ergibt Brownsche Bewegung
        W = np.zeros((n_paths, n_steps + 1))
        W[:, 1:] = np.cumsum(increments, axis=1)

        return t, W

    def covariance(self, t: float, s: float) -> float:
        """
        @brief Berechnet die Kovarianzfunktion der Brownschen Bewegung.

        @description
            Für die Standard-Brownsche Bewegung gilt:
              Cov(W(t), W(s)) = E[W(t)·W(s)] = min(t, s)

            Dies folgt aus der Eigenschaft unabhängiger Inkremente:
            Für t ≥ s: E[W(t)·W(s)] = E[(W(t)-W(s))·W(s)] + E[W(s)²]
                                     = 0 + s = s = min(t,s)

        @param t Erster Zeitpunkt (t ≥ 0)
        @param s Zweiter Zeitpunkt (s ≥ 0)
        @return Kovarianz Cov(W(t), W(s)) = min(t, s)
        @lastModified 2026-03-10
        """
        return min(t, s)

    def quadratic_variation(
        self,
        path: np.ndarray,
        dt: float
    ) -> float:
        """
        @brief Berechnet die empirische quadratische Variation eines Pfades.

        @description
            Die quadratische Variation [W]_T eines Wiener-Prozesses ist:
              [W]_T = lim_{|Π|→0} Σ (W(t_{k+1}) - W(t_k))²  = T  (f.s.)

            Empirisch: Σ_{k=0}^{n-1} (W(t_{k+1}) - W(t_k))²

            Dies konvergiert gegen T (die Gesamtzeit), nicht gegen 0!
            (Im Gegensatz zu differenzierbaren Funktionen.)

        @param path Diskreter Pfad W(t_0), W(t_1), ..., W(t_n)
        @param dt Zeitschrittweite Δt (für Information, nicht für Berechnung)
        @return Empirische quadratische Variation Σ(ΔW)²
        @lastModified 2026-03-10
        """
        increments = np.diff(path)
        return float(np.sum(increments ** 2))


# =============================================================================
# Klasse: ItoIntegral
# =============================================================================

class ItoIntegral:
    """
    @class ItoIntegral
    @brief Implementierung von Itô- und Stratonovich-Integralen.

    @description
        Das Itô-Integral ist definiert als:
          ∫₀ᵀ f(t) dW(t) = lim Σ f(t_k) · [W(t_{k+1}) - W(t_k)]

        wobei f(t_k) am LINKEN Endpunkt des Intervalls ausgewertet wird.

        Das Stratonovich-Integral verwendet den MITTLEREN Punkt:
          ∫₀ᵀ f(t) ∘ dW(t) = lim Σ f((t_k+t_{k+1})/2) · [W(t_{k+1}) - W(t_k)]

        Itô-Isometrie: E[|∫₀ᵀ f dW|²] = E[∫₀ᵀ f² dt]
        Itô-Formel: df(W) = f'(W)dW + ½f''(W)dt

    @lastModified 2026-03-10
    """

    def ito_integral(
        self,
        f_values: np.ndarray,
        brownian_increments: np.ndarray
    ) -> float:
        """
        @brief Berechnet das Itô-Integral ∫f dW als Links-Riemann-Summe.

        @description
            Itô-Integral (Auswertung am linken Rand):
              I = Σ_{k=0}^{n-1} f(t_k) · ΔW_k

            Dabei gilt: len(f_values) = n (linke Endpunkte)
                        len(brownian_increments) = n (Inkremente ΔW_k)

        @param f_values Integrand f(t_k) an den linken Punkten t_0,...,t_{n-1}
        @param brownian_increments Inkremente ΔW_k = W(t_{k+1}) - W(t_k)
        @return Wert des Itô-Integrals
        @lastModified 2026-03-10
        """
        f = np.asarray(f_values, dtype=float)
        dW = np.asarray(brownian_increments, dtype=float)

        if len(f) != len(dW):
            raise ValueError(
                f"f_values (len={len(f)}) und brownian_increments (len={len(dW)}) "
                "müssen gleiche Länge haben."
            )

        # Links-Riemann-Summe: Itô-Konvention
        return float(np.dot(f, dW))

    def stratonovich_integral(
        self,
        f_values: np.ndarray,
        brownian_increments: np.ndarray
    ) -> float:
        """
        @brief Berechnet das Stratonovich-Integral ∫f ∘ dW (Mittelpunktsregel).

        @description
            Stratonovich-Integral (Auswertung am Mittelpunkt):
              S = Σ_{k=0}^{n-1} ½(f(t_k) + f(t_{k+1})) · ΔW_k

            Dabei gilt:
            - len(f_values) = n+1 (Werte an allen Gitterpunkten t_0,...,t_n)
            - len(brownian_increments) = n (Inkremente)

            Zusammenhang mit Itô-Integral:
              ∫f ∘ dW = ∫f dW + ½∫f' dt  (Itô-Stratonovich-Korrektur)

        @param f_values Integrand an ALLEN Gitterpunkten t_0,...,t_n (Länge n+1)
        @param brownian_increments Inkremente ΔW_k (Länge n)
        @return Wert des Stratonovich-Integrals
        @lastModified 2026-03-10
        """
        f = np.asarray(f_values, dtype=float)
        dW = np.asarray(brownian_increments, dtype=float)

        n = len(dW)
        if len(f) != n + 1:
            raise ValueError(
                f"f_values muss Länge {n+1} haben (n+1 Punkte für n Inkremente)."
            )

        # Trapezregel: Mittelpunkt-Approximation f(t_{k+1/2}) ≈ (f_k + f_{k+1})/2
        f_mid = 0.5 * (f[:-1] + f[1:])
        return float(np.dot(f_mid, dW))

    def ito_formula_demo(
        self,
        f: Callable[[float], float],
        df: Callable[[float], float],
        ddf: Callable[[float], float],
        x0: float,
        T: float,
        n_steps: int
    ) -> Tuple[float, float]:
        """
        @brief Demonstriert die Itô-Formel numerisch.

        @description
            Die Itô-Formel für f(W(t)) lautet:
              f(W(T)) - f(W(0)) = ∫₀ᵀ f'(W)dW + ½∫₀ᵀ f''(W)dt

            Diese Methode:
            1. Simuliert eine Brownsche Bewegung W(t)
            2. Berechnet die linke Seite: f(W(T)) - f(W(0))
            3. Berechnet die rechte Seite: Itô-Integral + ½·Riemann-Integral
            4. Gibt beide Seiten zurück (sollten übereinstimmen)

        @param f Funktion f(x)
        @param df Erste Ableitung f'(x)
        @param ddf Zweite Ableitung f''(x)
        @param x0 Startwert W(0) = x0 (üblicherweise 0)
        @param T Endzeitpunkt
        @param n_steps Anzahl Diskretisierungsschritte
        @return Tupel (linke_seite, rechte_seite) der Itô-Formel
        @lastModified 2026-03-10
        """
        dt = T / n_steps

        # Simuliere Brownsche Bewegung mit festem Seed für Demo
        bm = BrownianMotion()
        t_grid, W = bm.simulate(T, n_steps, n_paths=1, seed=42)
        W = W[0]  # Einzelner Pfad

        # Brownsche Inkremente
        dW = np.diff(W)

        # Linke Seite: f(W(T)) - f(x0)
        lhs = f(W[-1]) - f(x0)

        # Rechte Seite: Itô-Integral + ½·Riemann-Integral
        # Itô-Integral: ∫f'(W)dW (linke Auswertung)
        ito_part = self.ito_integral(df(W[:-1]), dW)

        # Riemann-Integral: ½∫f''(W)dt
        riemann_part = 0.5 * np.sum(ddf(W[:-1])) * dt

        rhs = ito_part + riemann_part
        return float(lhs), float(rhs)

    def ito_isometry_check(
        self,
        f_values: np.ndarray,
        dt: float
    ) -> Tuple[float, float]:
        """
        @brief Überprüft die Itô-Isometrie empirisch.

        @description
            Die Itô-Isometrie besagt:
              E[|∫₀ᵀ f dW|²] = E[∫₀ᵀ f² dt] = ∫₀ᵀ E[f²] dt

            Für einen einzelnen Pfad approximieren wir:
            - Linke Seite: (Itô-Integral)²
            - Rechte Seite: Σ f_k² · Δt

            Bei vielen unabhängigen Realisierungen konvergieren beide Seiten.

        @param f_values Werte des Integranden f(t_0),...,f(t_{n-1})
        @param dt Zeitschrittweite
        @return Tupel (quadrat_des_ito_integrals, riemann_summe_f_quadrat)
        @lastModified 2026-03-10
        """
        f = np.asarray(f_values, dtype=float)
        n = len(f)

        # Simuliere Brownsche Inkremente
        rng = np.random.default_rng(42)
        dW = rng.normal(0, np.sqrt(dt), size=n)

        # Linke Seite: (∫f dW)²
        ito_val = self.ito_integral(f, dW)
        lhs = ito_val ** 2

        # Rechte Seite: ∫f² dt ≈ Σ f_k² · Δt
        rhs = float(np.sum(f ** 2) * dt)

        return float(lhs), float(rhs)


# =============================================================================
# Klasse: StochasticDifferentialEquation
# =============================================================================

class StochasticDifferentialEquation:
    """
    @class StochasticDifferentialEquation
    @brief Numerische Lösung von SDEs der Form dX = μ(X)dt + σ(X)dW.

    @description
        Eine stochastische Differentialgleichung (SDE) in Itô-Form:
          dX(t) = μ(X(t), t) dt + σ(X(t), t) dW(t)

        Implementierte Verfahren:
        1. Euler-Maruyama: Ordnung O(√Δt) (stark), O(Δt) (schwach)
           X_{k+1} = X_k + μ(X_k)Δt + σ(X_k)ΔW_k
        2. Milstein: Ordnung O(Δt) (stark)
           X_{k+1} = X_k + μΔt + σΔW + ½σσ'((ΔW)² - Δt)
        3. Geometrische Brownsche Bewegung (exakte Lösung):
           S(t) = S₀ exp((μ - σ²/2)t + σW(t))

    @lastModified 2026-03-10
    """

    def euler_maruyama(
        self,
        mu: Callable[[float], float],
        sigma: Callable[[float], float],
        x0: float,
        T: float,
        n_steps: int,
        n_paths: int = 1,
        seed: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Euler-Maruyama-Verfahren für SDEs.

        @description
            Das Euler-Maruyama-Schema:
              X_{k+1} = X_k + μ(X_k)·Δt + σ(X_k)·ΔW_k

            mit ΔW_k ~ N(0, Δt) i.i.d.

            Starke Konvergenzordnung: O(√Δt)
            Schwache Konvergenzordnung: O(Δt)

        @param mu Driftfunktion μ(x)
        @param sigma Diffusionsfunktion σ(x)
        @param x0 Anfangswert X(0)
        @param T Endzeitpunkt
        @param n_steps Anzahl Zeitschritte
        @param n_paths Anzahl Simulationspfade
        @param seed Zufallssaat
        @return Tupel (t_grid, X): Zeitgitter und Pfade (n_paths × n_steps+1)
        @lastModified 2026-03-10
        """
        dt = T / n_steps
        t_grid = np.linspace(0, T, n_steps + 1)

        rng = np.random.default_rng(seed)

        # Initialzustand für alle Pfade
        X = np.zeros((n_paths, n_steps + 1))
        X[:, 0] = x0

        # Euler-Maruyama-Schritte
        for k in range(n_steps):
            # Brownsche Inkremente ΔW_k ~ N(0, √dt)
            dW = rng.normal(0, np.sqrt(dt), size=n_paths)
            x_curr = X[:, k]
            # Vektorisierte Berechnung (mu/sigma pfadweise)
            mu_vals = np.array([mu(xi) for xi in x_curr])
            sigma_vals = np.array([sigma(xi) for xi in x_curr])
            X[:, k + 1] = x_curr + mu_vals * dt + sigma_vals * dW

        return t_grid, X

    def milstein(
        self,
        mu: Callable[[float], float],
        sigma: Callable[[float], float],
        sigma_deriv: Callable[[float], float],
        x0: float,
        T: float,
        n_steps: int,
        seed: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Milstein-Verfahren für SDEs (höhere Konvergenzordnung).

        @description
            Das Milstein-Schema ergänzt Euler-Maruyama um einen Korrekturterm:
              X_{k+1} = X_k + μ(X_k)Δt + σ(X_k)ΔW_k
                       + ½σ(X_k)σ'(X_k)[(ΔW_k)² - Δt]

            Der Korrekturterm ½σσ'[(ΔW)²-Δt] berücksichtigt die quadratische
            Variation der Brownschen Bewegung (Itô-Korrektur).

            Starke Konvergenzordnung: O(Δt)  (besser als Euler-Maruyama O(√Δt))

        @param mu Driftfunktion μ(x)
        @param sigma Diffusionsfunktion σ(x)
        @param sigma_deriv Ableitung σ'(x) der Diffusionsfunktion
        @param x0 Anfangswert X(0)
        @param T Endzeitpunkt
        @param n_steps Anzahl Zeitschritte
        @param seed Zufallssaat
        @return Tupel (t_grid, X): Zeitgitter und Pfad (1 × n_steps+1)
        @lastModified 2026-03-10
        """
        dt = T / n_steps
        t_grid = np.linspace(0, T, n_steps + 1)

        rng = np.random.default_rng(seed)

        X = np.zeros(n_steps + 1)
        X[0] = x0

        for k in range(n_steps):
            dW = rng.normal(0, np.sqrt(dt))
            x_curr = X[k]
            mu_k = mu(x_curr)
            sigma_k = sigma(x_curr)
            sigma_d_k = sigma_deriv(x_curr)

            # Milstein-Korrekturterm: ½σσ'((ΔW)² - Δt)
            correction = 0.5 * sigma_k * sigma_d_k * (dW ** 2 - dt)
            X[k + 1] = x_curr + mu_k * dt + sigma_k * dW + correction

        return t_grid, X.reshape(1, -1)

    def geometric_brownian_motion(
        self,
        mu: float,
        sigma: float,
        S0: float,
        T: float,
        n_steps: int,
        n_paths: int = 1,
        seed: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Simuliert geometrische Brownsche Bewegung (Black-Scholes-Modell).

        @description
            Die geometrische Brownsche Bewegung (GBM) löst die SDE:
              dS = μ·S dt + σ·S dW

            Die exakte Lösung lautet:
              S(t) = S₀ · exp((μ - σ²/2)·t + σ·W(t))

            GBM wird in der Finanzmathematik (Black-Scholes) und in der
            Modellierung von Aktienpreisen verwendet.
            S(t) ist immer positiv (geometrisch, nicht arithmetisch).

        @param mu Driftrate (erwartet. logarithm. Rendite)
        @param sigma Volatilität (Standardabweichung der log-Rendite)
        @param S0 Anfangswert S(0) > 0
        @param T Endzeitpunkt (z.B. in Jahren)
        @param n_steps Anzahl Zeitschritte
        @param n_paths Anzahl Simulationspfade
        @param seed Zufallssaat
        @return Tupel (t_grid, S): Zeitgitter und Preispfade (n_paths × n_steps+1)
        @lastModified 2026-03-10
        """
        dt = T / n_steps
        t_grid = np.linspace(0, T, n_steps + 1)

        rng = np.random.default_rng(seed)

        # Simuliere Brownsche Bewegung
        Z = rng.normal(0, 1, size=(n_paths, n_steps))
        dW = np.sqrt(dt) * Z

        # Log-Renditen: Δlog(S) = (μ - σ²/2)Δt + σΔW
        log_returns = (mu - 0.5 * sigma ** 2) * dt + sigma * dW

        # Kumulierte Summe der Log-Renditen
        S = np.zeros((n_paths, n_steps + 1))
        S[:, 0] = S0
        S[:, 1:] = S0 * np.exp(np.cumsum(log_returns, axis=1))

        return t_grid, S


# =============================================================================
# Klasse: ErgodicTheory
# =============================================================================

class ErgodicTheory:
    """
    @class ErgodicTheory
    @brief Ergodentheorie: Zeitmittel, Birkhoff-Theorem, Lyapunov-Exponenten.

    @description
        Die Ergodentheorie untersucht das Langzeitverhalten dynamischer Systeme.

        Birkhoff'sches Ergodentheorem:
          Für einen ergodischen Prozess konvergiert das Zeitmittel gegen
          das Raumittel (Erwartungswert):
            (1/N) Σ_{k=0}^{N-1} f(X_k) → E[f(X)]  (f.s.)

        Lyapunov-Exponent λ misst die exponentielle Divergenz naher Trajektorien:
          λ = lim_{n→∞} (1/n) Σ log|f'(x_k)|

        Für λ > 0: Chaotisches Verhalten (exponentielle Sensitivität)

    @lastModified 2026-03-10
    """

    def birkhoff_ergodic_theorem_demo(
        self,
        f_values: np.ndarray
    ) -> np.ndarray:
        """
        @brief Demonstriert das Birkhoff'sche Ergodentheorem.

        @description
            Berechnet die laufenden Zeitmittel:
              A_n = (1/n) Σ_{k=0}^{n-1} f(X_k)

            Für ergodische Prozesse konvergiert A_n gegen E[f] (den Raummittelwert).
            Die Konvergenz ist für n→∞ garantiert (Birkhoff, 1931).

        @param f_values Beobachtungen f(X_0), f(X_1), ..., f(X_{N-1})
        @return Array der laufenden Zeitmittel A_1, A_2, ..., A_N
        @lastModified 2026-03-10
        """
        f = np.asarray(f_values, dtype=float)
        # Kumulierte Summen dividiert durch Index
        cumsum = np.cumsum(f)
        n = np.arange(1, len(f) + 1)
        return cumsum / n

    def time_average(self, f_values: np.ndarray) -> float:
        """
        @brief Berechnet das Zeitmittel einer Beobachtungsreihe.

        @description
            Zeitmittel: T̄ = (1/N) Σ_{k=0}^{N-1} f(X_k)

            Im ergodischen Fall: T̄ → E[f] für N → ∞.

        @param f_values Werte f(X_0), ..., f(X_{N-1})
        @return Zeitmittel (1/N)·Σf_k
        @lastModified 2026-03-10
        """
        return float(np.mean(f_values))

    def is_mixing_demo(self, sequence: np.ndarray, lag: int = 10) -> bool:
        """
        @brief Einfache Mischungsprüfung via Autokorrelation.

        @description
            Ein Prozess ist mischend, wenn zeitlich weit entfernte Ereignisse
            asymptotisch unabhängig werden:
              lim_{n→∞} P(A ∩ T^{-n}B) = P(A)·P(B)

            Approximation: Überprüfe, ob die Autokorrelation bei Lag=lag
            nahe 0 ist (schnelles Abklingen).

            Eine Sequenz gilt hier als "mischend", wenn:
            |Korr(X_k, X_{k+lag})| < 0.1 (empirisch)

        @param sequence Zeitreihe X_0, X_1, ..., X_{N-1}
        @param lag Verzögerung für die Autokorrelation (Standard: 10)
        @return True wenn Autokorrelation bei gegebenem Lag klein
        @lastModified 2026-03-10
        """
        seq = np.asarray(sequence, dtype=float)
        if len(seq) <= lag:
            return False

        # Autokorrelation bei gegebenem Lag
        x = seq - np.mean(seq)
        acf_0 = np.dot(x, x)
        if acf_0 < 1e-14:
            return True  # Konstante Sequenz

        acf_lag = np.dot(x[:-lag], x[lag:])
        autocorr = acf_lag / acf_0

        # Mischend wenn Autokorrelation < 0.1 in Betrag
        return bool(abs(autocorr) < 0.1)

    def logistic_map(
        self,
        r: float,
        x0: float,
        n_iter: int
    ) -> np.ndarray:
        """
        @brief Iteriert die logistische Abbildung f(x) = r·x·(1-x).

        @description
            Die logistische Abbildung:
              x_{n+1} = r · x_n · (1 - x_n),  x_0 ∈ (0, 1)

            Für verschiedene Parameter r zeigt sie unterschiedliches Verhalten:
            - 0 < r ≤ 1: Konvergenz zu 0
            - 1 < r ≤ 3: Stabile Fixpunkte
            - 3 < r ≤ 3.57: Periodenverdopplung (Bifurkation)
            - r ≈ 4: Chaotisches Verhalten (ergodisch, mischend)

            Der Fixpunkt x* = 1 - 1/r ist stabil für r < 3.

        @param r Wachstumsrate r ∈ (0, 4]
        @param x0 Anfangswert x_0 ∈ (0, 1)
        @param n_iter Anzahl Iterationen
        @return Array [x_0, x_1, ..., x_{n_iter}]
        @raises ValueError Wenn r außerhalb (0, 4] oder x0 außerhalb (0, 1)
        @lastModified 2026-03-10
        """
        if not (0 < r <= 4):
            raise ValueError(f"Parameter r={r} muss in (0, 4] liegen.")
        if not (0 < x0 < 1):
            raise ValueError(f"Startwert x0={x0} muss in (0, 1) liegen.")

        x = np.zeros(n_iter + 1)
        x[0] = x0

        for n in range(n_iter):
            x[n + 1] = r * x[n] * (1 - x[n])

        return x

    def lyapunov_exponent(
        self,
        r: float,
        x0: float,
        n_iter: int
    ) -> float:
        """
        @brief Berechnet den Lyapunov-Exponenten der logistischen Abbildung.

        @description
            Der Lyapunov-Exponent misst die Sensitivität gegenüber Anfangsbedingungen:
              λ = lim_{n→∞} (1/n) Σ_{k=0}^{n-1} log|f'(x_k)|

            Für f(x) = r·x·(1-x) gilt:
              f'(x) = r·(1 - 2x)

            Interpretation:
            - λ < 0: Stabiler Attraktor (Konvergenz)
            - λ = 0: Bifurkationspunkt
            - λ > 0: Chaotisches Verhalten (exponentielle Divergenz)

            Für r = 4 gilt: λ = log(2) ≈ 0.693

        @param r Wachstumsrate r ∈ (0, 4]
        @param x0 Anfangswert x_0 ∈ (0, 1)
        @param n_iter Anzahl Iterationen für Mittelwert
        @return Lyapunov-Exponent λ
        @lastModified 2026-03-10
        """
        trajectory = self.logistic_map(r, x0, n_iter)

        # Ableitungen: f'(x) = r(1 - 2x)
        derivs = r * (1 - 2 * trajectory[:-1])

        # Vermeide log(0) bei x = 0.5
        abs_derivs = np.abs(derivs)
        abs_derivs = np.where(abs_derivs < 1e-14, 1e-14, abs_derivs)

        # Lyapunov-Exponent als Zeitmittel der log-Ableitungen
        return float(np.mean(np.log(abs_derivs)))


# =============================================================================
# Klasse: PoissonProcess
# =============================================================================

class PoissonProcess:
    """
    @class PoissonProcess
    @brief Homogener Poisson-Prozess mit Intensität λ.

    @description
        Ein (homogener) Poisson-Prozess N(t) mit Rate λ > 0 hat:
        - N(0) = 0
        - Unabhängige Inkremente
        - N(t) - N(s) ~ Poisson(λ(t-s)) für 0 ≤ s < t
        - P(N(t) = k) = e^{-λt} (λt)^k / k!
        - E[N(t)] = Var[N(t)] = λt

        Simulationsalgorithmus: Zwischenankunftszeiten T_k ~ Exp(λ),
        Sprungzeiten S_k = T_1 + ... + T_k.

    @param rate Intensitätsparameter λ > 0
    @lastModified 2026-03-10
    """

    def __init__(self, rate: float) -> None:
        """
        @brief Initialisiert den Poisson-Prozess mit gegebener Rate λ.

        @param rate Intensität λ > 0
        @raises ValueError Wenn rate ≤ 0
        @lastModified 2026-03-10
        """
        if rate <= 0:
            raise ValueError(f"Rate λ={rate} muss positiv sein.")
        self.rate = rate

    def simulate(
        self,
        T: float,
        seed: Optional[int] = None
    ) -> np.ndarray:
        """
        @brief Simuliert die Sprungzeiten des Poisson-Prozesses auf [0, T].

        @description
            Algorithmus via Intertimes:
            1. Ziehe Zwischenankunftszeiten T_k ~ Exp(λ)
            2. Sprungzeiten: S_k = T_1 + ... + T_k
            3. Behalte alle S_k ≤ T

            Die Anzahl der Sprünge N(T) ~ Poisson(λT).

        @param T Endzeitpunkt
        @param seed Zufallssaat für Reproduzierbarkeit
        @return Array der Sprungzeiten [S_1, S_2, ..., S_{N(T)}]
        @lastModified 2026-03-10
        """
        if T <= 0:
            raise ValueError(f"T={T} muss positiv sein.")

        rng = np.random.default_rng(seed)

        # Erwartete Anzahl Sprünge + Sicherheitspuffer
        expected_jumps = int(self.rate * T)
        n_sample = max(10, int(2 * expected_jumps + 5 * np.sqrt(expected_jumps + 1)))

        # Ziehe Zwischenankunftszeiten ~ Exp(λ)
        inter_times = rng.exponential(1.0 / self.rate, size=n_sample)
        jump_times = np.cumsum(inter_times)

        # Behalte nur Sprünge vor T
        valid = jump_times[jump_times <= T]

        # Falls nicht genug Sprünge generiert wurden, iteriere
        while len(valid) == n_sample:
            more = rng.exponential(1.0 / self.rate, size=n_sample)
            more_jumps = np.cumsum(more) + jump_times[-1]
            new_valid = more_jumps[more_jumps <= T]
            valid = np.concatenate([valid, new_valid])
            if len(new_valid) == 0:
                break

        return valid

    def pmf(self, k: int, t: float) -> float:
        """
        @brief Berechnet P(N(t) = k) = e^{-λt}(λt)^k / k!

        @description
            Wahrscheinlichkeitsmassfunktion des Poisson-Prozesses:
              P(N(t) = k) = e^{-λt} · (λt)^k / k!

            Berechnung in Log-Skala für numerische Stabilität:
              log P = -λt + k·log(λt) - log(k!)

        @param k Anzahl der Ereignisse (k ≥ 0, ganzzahlig)
        @param t Zeitpunkt t > 0
        @return Wahrscheinlichkeit P(N(t) = k)
        @raises ValueError Wenn k < 0 oder t ≤ 0
        @lastModified 2026-03-10
        """
        if k < 0:
            raise ValueError(f"k={k} muss nicht-negativ sein.")
        if t <= 0:
            raise ValueError(f"t={t} muss positiv sein.")

        import math
        lam_t = self.rate * t

        # Log-Wahrscheinlichkeit für numerische Stabilität
        log_prob = -lam_t + k * np.log(lam_t) - math.lgamma(k + 1)
        return float(np.exp(log_prob))

    def mean(self, t: float) -> float:
        """
        @brief Erwartungswert E[N(t)] = λt.

        @param t Zeitpunkt t > 0
        @return E[N(t)] = λ·t
        @lastModified 2026-03-10
        """
        return self.rate * t

    def variance(self, t: float) -> float:
        """
        @brief Varianz Var[N(t)] = λt.

        @description
            Für den Poisson-Prozess gilt die Eigenschaft:
              E[N(t)] = Var[N(t)] = λt
            (Gleichheit von Erwartungswert und Varianz)

        @param t Zeitpunkt t > 0
        @return Var[N(t)] = λ·t
        @lastModified 2026-03-10
        """
        return self.rate * t


# =============================================================================
# Klasse: GaussianProcess
# =============================================================================

class GaussianProcess:
    """
    @class GaussianProcess
    @brief Gaußscher Prozess (GP) mit Mittelwertfunktion und Kernelfunktion.

    @description
        Ein Gaußscher Prozess ist eine Sammlung von Zufallsvariablen
        f(x_1), ..., f(x_n), deren gemeinsame Verteilung für jede endliche
        Auswahl von Eingaben multivariat normal ist:

          (f(x_1), ..., f(x_n)) ~ N(μ, K)

        wobei:
        - μ_i = m(x_i) (Mittelwertfunktion)
        - K_{ij} = k(x_i, x_j) (Kernfunktion/Kovarianzfunktion)

        GPs werden in der Bayes'schen nicht-parametrischen Statistik und
        im Maschinellen Lernen (Gaussian Process Regression) verwendet.

    @param mean_func Mittelwertfunktion m(x) → ℝ
    @param kernel_func Kernfunktion k(x, y) → ℝ (positiv semidefinit)
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        mean_func: Callable[[float], float],
        kernel_func: Callable[[float, float], float]
    ) -> None:
        """
        @brief Initialisiert den Gaußschen Prozess.

        @param mean_func Mittelwertfunktion m: ℝ → ℝ
        @param kernel_func Kernfunktion k: ℝ × ℝ → ℝ (positiv semidefinit)
        @lastModified 2026-03-10
        """
        self.mean_func = mean_func
        self.kernel_func = kernel_func

    def sample(
        self,
        x_points: np.ndarray,
        seed: Optional[int] = None
    ) -> np.ndarray:
        """
        @brief Zieht eine Stichprobe aus dem Gaußschen Prozess.

        @description
            Berechnet die Kovarianzmatrix K an den Punkten x_points:
              K_{ij} = k(x_i, x_j)

            Dann wird eine multivariate Normalverteilung gezogen:
              f ~ N(μ, K)

            Numerische Stabilität: Füge kleine Diagonalstörung hinzu
            (Jitter), um positive Definitheit sicherzustellen.

        @param x_points Eingabepunkte x_1, ..., x_n
        @param seed Zufallssaat
        @return Stichprobe f(x_1), ..., f(x_n) als 1D-Array
        @lastModified 2026-03-10
        """
        x = np.asarray(x_points, dtype=float)
        n = len(x)

        # Mittelwertvektor
        mu = np.array([self.mean_func(xi) for xi in x])

        # Kovarianzmatrix K_{ij} = k(x_i, x_j)
        K = np.array([[self.kernel_func(x[i], x[j]) for j in range(n)]
                      for i in range(n)])

        # Numerischer Jitter für positive Definitheit
        K += 1e-10 * np.eye(n)

        rng = np.random.default_rng(seed)

        # Cholesky-Zerlegung für stabiles Sampling
        try:
            L = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            # Fallback: Eigenwert-basiertes Sampling
            eigvals, eigvecs = np.linalg.eigh(K)
            eigvals = np.maximum(eigvals, 0)  # Numerisch negative → 0
            L = eigvecs @ np.diag(np.sqrt(eigvals))

        # Stichprobe: f = μ + L·z, z ~ N(0, I)
        z = rng.normal(0, 1, size=n)
        return mu + L @ z

    @staticmethod
    def rbf_kernel(
        x1: float,
        x2: float,
        length_scale: float = 1.0
    ) -> float:
        """
        @brief RBF/Gauß-Kern (Radial Basis Function Kernel).

        @description
            Der RBF-Kern (auch Gauß-Kern oder SE-Kern):
              k(x₁, x₂) = exp(-|x₁ - x₂|² / (2·ℓ²))

            Eigenschaften:
            - Stationär: k hängt nur von |x₁ - x₂| ab
            - Glatt: Unendlich oft differenzierbar
            - Universell: Kann beliebige stetige Funktionen approximieren
            - ℓ = length_scale steuert die Korrelationslänge

        @param x1 Erster Eingabepunkt
        @param x2 Zweiter Eingabepunkt
        @param length_scale Längenskalenparameter ℓ > 0 (Standard: 1.0)
        @return Kernwert k(x₁, x₂) ∈ (0, 1]
        @lastModified 2026-03-10
        """
        return float(np.exp(-0.5 * ((x1 - x2) / length_scale) ** 2))


# =============================================================================
# Standalone-Funktionen
# =============================================================================

def random_walk_1d(
    n_steps: int,
    p: float = 0.5,
    seed: Optional[int] = None
) -> np.ndarray:
    """
    @brief Simuliert einen 1D-Random-Walk S_n = Σ X_k.

    @description
        Der symmetrische 1D-Random-Walk:
          S_n = X_1 + X_2 + ... + X_n

        wobei X_k = +1 mit Wahrscheinlichkeit p
              X_k = -1 mit Wahrscheinlichkeit 1-p

        Für p = 0.5 (symmetrisch):
        - E[S_n] = 0
        - Var[S_n] = n
        - Skalierungslimit: S_{⌊nt⌋}/√n → W(t) (Brownsche Bewegung)

    @param n_steps Anzahl Schritte
    @param p Wahrscheinlichkeit für Schritt +1 (Standard: 0.5)
    @param seed Zufallssaat
    @return Array [S_0, S_1, ..., S_n] mit S_0 = 0
    @raises ValueError Wenn p ∉ [0, 1]
    @lastModified 2026-03-10
    """
    if not (0 <= p <= 1):
        raise ValueError(f"Wahrscheinlichkeit p={p} muss in [0, 1] liegen.")

    rng = np.random.default_rng(seed)

    # Schritte: +1 oder -1
    steps = rng.choice([1, -1], size=n_steps, p=[p, 1 - p])

    # Kumulierte Summe (mit S_0 = 0)
    walk = np.concatenate([[0], np.cumsum(steps)])
    return walk


def gambler_ruin_probability(
    p: float,
    start: int,
    target: int
) -> float:
    """
    @brief Berechnet die Ruinwahrscheinlichkeit beim Glücksspiel (Gambler's Ruin).

    @description
        Das Gambler's-Ruin-Problem:
        Ein Spieler beginnt mit k Einheiten, gewinnt 1 mit Wahrsch. p,
        verliert 1 mit Wahrsch. q = 1-p. Das Spiel endet bei 0 (Ruin)
        oder N (Gewinn).

        Ruinwahrscheinlichkeit (Verlust-Wahrscheinlichkeit) R(k):
        - Für p ≠ 0.5:
            R(k) = [1 - (q/p)^k] / [1 - (q/p)^N]
            (Wahrsch. N zu erreichen, nicht 0)
        - Für p = 0.5 (symmetrisch):
            R(k) = k/N  (Wahrsch. N zu erreichen)

        HINWEIS: Diese Funktion gibt die Wahrscheinlichkeit zurück,
        das Ziel N zu erreichen (= Gewinnwahrscheinlichkeit), NICHT zu ruinieren.

    @param p Gewinnwahrscheinlichkeit pro Runde ∈ [0, 1]
    @param start Startkapital k ∈ {0, 1, ..., target}
    @param target Zielkapital N > start
    @return Wahrscheinlichkeit, Ziel N zu erreichen (vor Ruin bei 0)
    @raises ValueError Wenn p ∉ [0, 1] oder start/target ungültig
    @lastModified 2026-03-10
    """
    if not (0 <= p <= 1):
        raise ValueError(f"p={p} muss in [0, 1] liegen.")
    if start < 0 or target <= 0 or start > target:
        raise ValueError(
            f"Ungültige Werte: start={start}, target={target}. "
            "Benötigt: 0 ≤ start ≤ target, target > 0."
        )

    k = start
    N = target
    q = 1 - p

    # Randfälle
    if k == 0:
        return 0.0  # Bereits ruiniert
    if k == N:
        return 1.0  # Bereits am Ziel

    # Randfall p = 0: sicherer Verlust → Gewinnwahrscheinlichkeit = 0
    if p < 1e-14:
        return 0.0

    # Randfall p = 1: sicherer Gewinn → Gewinnwahrscheinlichkeit = 1
    if p > 1 - 1e-14:
        return 1.0

    # Symmetrischer Fall p = q = 0.5
    if abs(p - 0.5) < 1e-12:
        return k / N

    # Allgemeiner Fall (p ≠ 0.5)
    r = q / p  # Verhältnis q/p
    return (1 - r ** k) / (1 - r ** N)


def central_limit_theorem_demo(
    n_samples: int,
    n_obs: int = 1000,
    distribution: str = 'uniform',
    seed: Optional[int] = None
) -> Tuple[np.ndarray, float, float]:
    """
    @brief Demonstriert den Zentralen Grenzwertsatz (CLT) numerisch.

    @description
        Der Zentrale Grenzwertsatz (CLT):
        Für i.i.d. Zufallsvariablen X_1, ..., X_n mit E[X_i] = μ, Var[X_i] = σ²:

          √n · (X̄_n - μ) / σ → N(0, 1)  (in Verteilung, n → ∞)

        Diese Funktion:
        1. Zieht n_samples Stichproben der Größe n_obs
        2. Berechnet die Stichprobenmittel X̄_1, ..., X̄_{n_samples}
        3. Normiert auf Z_i = √n·(X̄_i - μ)/σ
        4. Gibt die normierten Mittelwerte sowie empirischen Mittelwert/Varianz zurück

        Verfügbare Verteilungen: 'uniform' U(0,1), 'exponential' Exp(1), 'bernoulli' B(0.5)

    @param n_samples Anzahl der Stichprobenmittel (Realisierungen)
    @param n_obs Größe jeder Stichprobe n
    @param distribution Grundverteilung ('uniform', 'exponential', 'bernoulli')
    @param seed Zufallssaat
    @return Tupel (normalized_means, empirical_mean, empirical_variance)
    @raises ValueError Wenn distribution unbekannt
    @lastModified 2026-03-10
    """
    rng = np.random.default_rng(seed)

    # Parameter der Grundverteilungen (μ, σ)
    if distribution == 'uniform':
        # U(0, 1): μ = 0.5, σ = 1/√12
        samples = rng.uniform(0, 1, size=(n_samples, n_obs))
        mu_true = 0.5
        sigma_true = 1.0 / np.sqrt(12)
    elif distribution == 'exponential':
        # Exp(1): μ = 1, σ = 1
        samples = rng.exponential(1.0, size=(n_samples, n_obs))
        mu_true = 1.0
        sigma_true = 1.0
    elif distribution == 'bernoulli':
        # Bernoulli(0.5): μ = 0.5, σ = 0.5
        samples = rng.choice([0, 1], size=(n_samples, n_obs), p=[0.5, 0.5]).astype(float)
        mu_true = 0.5
        sigma_true = 0.5
    else:
        raise ValueError(
            f"Unbekannte Verteilung '{distribution}'. "
            "Wähle aus: 'uniform', 'exponential', 'bernoulli'."
        )

    # Stichprobenmittel
    sample_means = np.mean(samples, axis=1)

    # Normierung: Z = √n · (X̄ - μ) / σ
    normalized = np.sqrt(n_obs) * (sample_means - mu_true) / sigma_true

    empirical_mean = float(np.mean(normalized))
    empirical_var = float(np.var(normalized))

    return normalized, empirical_mean, empirical_var
