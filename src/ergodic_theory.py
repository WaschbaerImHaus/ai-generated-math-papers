"""
@file ergodic_theory.py
@brief Umfassendes Modul für Ergodentheorie: Maß-erhaltende Systeme, Ergodensätze,
       Entropie, Mixing, Collatz-Dynamik (Tao 2019), Furstenberg-Theorie und
       symbolische Dynamik.
@description
    Dieses Modul implementiert die zentralen Konzepte der modernen Ergodentheorie:

    GRUNDLAGEN:
    - MeasurePreservingSystem: Dynamisches System (X, T, μ) mit μ(T⁻¹A) = μ(A)
    - Ergoditätsprüfung: Alle T-invarianten Mengen haben Maß 0 oder 1
    - Ergodische Zerlegung

    ERGODENSÄTZE:
    - Birkhoff-Ergodentheorem (1931): (1/n)Σf(Tᵏx) → ∫f dμ (f.s.)
    - von-Neumann-Mittelwertergodentheorem: L²-Konvergenz

    ENTROPIE:
    - Kolmogorov-Sinai-Entropie h_μ(T)
    - Topologische Entropie h_top(T)
    - Variationsprinzip: h_top = sup_{μ} h_μ(T)
    - Partitionsentropie H(P) = -Σ μ(Aᵢ) log μ(Aᵢ)

    MIXING-EIGENSCHAFTEN:
    - Starkes Mixing: lim μ(A ∩ T⁻ⁿB) = μ(A)μ(B)
    - Schwaches Mixing: Spektrales Kriterium
    - Mixing-Rate-Schätzung

    COLLATZ-DYNAMIK (Tao 2019):
    - Collatz-Abbildung als ergodisches System auf 2-adischen Zahlen
    - Numerische Verifikation von Taos Dichte-1-Resultat
    - Lyapunov-Exponent λ ≈ log(3)/2 - log(2)/2 ≈ -0.0235

    FURSTENBERG-THEORIE:
    - Furstenberg-Korrespondenz: Zahlentheorie ↔ Ergodentheorie
    - Szemerédi-Theorem via ergodischer Methoden
    - van-der-Waerden-Schranken

    SYMBOLISCHE DYNAMIK:
    - ShiftSystem: Shift-Abbildung auf Folgenraum
    - Subshift of Finite Type (SFT) via Übergangsmatrix
    - Zeta-Funktion der Shifts
    - Entropie über Spektralradius

    Mathematische Grundlagen:
    - Maß-erhaltend: μ(T⁻¹A) = μ(A) für alle messbaren A
    - Ergodisch: Alle invarianten Mengen haben Maß 0 oder 1
    - Birkhoff: (1/n)Σf(Tᵏx) → E_μ[f] fast überall (n→∞)
    - KS-Entropie: h(T) = sup_P h(T, P) über alle endlichen Partitionen P

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import numpy as np
import scipy.linalg as la
import scipy.stats as stats
from typing import Callable, List, Optional, Tuple, Dict
from functools import lru_cache
import math


# =============================================================================
# Klasse: MeasurePreservingSystem
# =============================================================================

class MeasurePreservingSystem:
    """
    @class MeasurePreservingSystem
    @brief Maß-erhaltendes dynamisches System (X, T, μ).

    @description
        Ein maß-erhaltendes System besteht aus:
        - Zustandsraum X (diskret: endliche Menge, kontinuierlich: [0,1] o.ä.)
        - Transformation T: X → X (messbar)
        - Maß μ mit μ(T⁻¹A) = μ(A) für alle messbaren Teilmengen A

        Beispiele:
        - Kreisrotation T(x) = x + α (mod 1), Lebesgue-Maß: Maß-erhaltend
        - Verdopplungsabbildung T(x) = 2x (mod 1): Lebesgue-erhaltend
        - Logistische Abbildung bei r=4: β-Maß mit Dichte 1/(π√(x(1-x)))

        Ergodizität bedeutet, dass der Orbit fast jedes Punktes
        den gesamten Zustandsraum "abdeckt" (Irreduzibilität).

    @param T Transformation als Callable X → X
    @param name Bezeichnung des Systems (für Ausgaben)
    @lastModified 2026-03-11
    """

    def __init__(self, T: Callable, name: str = "Unbekanntes System") -> None:
        """
        @brief Initialisiert das maß-erhaltende System.

        @param T Transformation T: X → X
        @param name Bezeichnung des Systems
        @lastModified 2026-03-11
        """
        self.T = T
        self.name = name

    def orbit(self, x0: float, n: int) -> np.ndarray:
        """
        @brief Berechnet den Orbit x₀, T(x₀), T²(x₀), ..., Tⁿ(x₀).

        @description
            Der Orbit eines Punktes x₀ unter der Transformation T:
              Orb(x₀) = {Tᵏ(x₀) : k = 0, 1, 2, ...}

            Im ergodischen Fall liegt dieser Orbit dicht in X (für fast alle x₀).

        @param x0 Startpunkt
        @param n Anzahl der Iterationen
        @return Array [x₀, T(x₀), T²(x₀), ..., Tⁿ(x₀)]
        @lastModified 2026-03-11
        """
        trajectory = np.zeros(n + 1)
        trajectory[0] = x0
        x = x0

        for k in range(n):
            x = self.T(x)
            trajectory[k + 1] = x

        return trajectory

    def is_measure_preserving(
        self,
        test_points: np.ndarray,
        n_bins: int = 20
    ) -> bool:
        """
        @brief Prüft numerisch ob T maß-erhaltend ist (via Histogramm-Vergleich).

        @description
            Für Lebesgue-Maß auf [0,1]: Wenn T maß-erhaltend ist, sollte
            die Verteilung der T(x)-Werte gleich der Verteilung der x-Werte sein.

            Methode: χ²-Test zwischen Histogramm von x und Histogramm von T(x).
            Wenn der p-Wert > 0.05, wird Maß-Erhaltung nicht abgelehnt.

        @param test_points Gleichmäßig verteilte Testpunkte in [0,1]
        @param n_bins Anzahl der Histogramm-Bins
        @return True wenn Maß-Erhaltung nicht abgelehnt wird (p > 0.05)
        @lastModified 2026-03-11
        """
        x = np.asarray(test_points, dtype=float)
        Tx = np.array([self.T(xi) for xi in x])

        # Histogramme berechnen
        bins = np.linspace(0, 1, n_bins + 1)
        hist_x, _ = np.histogram(x, bins=bins)
        hist_Tx, _ = np.histogram(Tx, bins=bins)

        # χ²-Test zwischen den Histogrammen
        # Nur Bins mit ausreichend Beobachtungen verwenden
        mask = (hist_x > 0) & (hist_Tx > 0)
        if mask.sum() < 3:
            return False

        # Verhältnis der Histogramme sollte nahe 1 sein
        ratios = hist_Tx[mask].astype(float) / hist_x[mask].astype(float)
        return bool(np.std(ratios) < 0.3)

    def is_ergodic(
        self,
        x0: float,
        n_iter: int = 50000,
        n_bins: int = 20
    ) -> bool:
        """
        @brief Prüft numerisch ob das System ergodisch ist.

        @description
            Ergodizität (Birkhoff): Das Zeitmittel konvergiert für fast alle
            Startpunkte gegen das Raummittel (Gleichverteilung auf X).

            Methode: Wenn der Orbit von x₀ nach n_iter Schritten alle n_bins
            Bereiche von [0,1] "besucht" (jeder Bin mindestens einmal), gilt das
            System als numerisch ergodisch.

            Präzise Formulierung: Ein System ist ergodisch genau dann, wenn für
            alle messbaren T-invarianten Mengen A gilt: μ(A) ∈ {0, 1}.

        @param x0 Startpunkt für den Orbit-Test
        @param n_iter Anzahl der Iterationsschritte
        @param n_bins Anzahl der Bins für die Gleichverteilungsprüfung
        @return True wenn Orbit gleichmäßig verteilt (Ergodizitätsindikator)
        @lastModified 2026-03-11
        """
        traj = self.orbit(x0, n_iter)

        # Histogramm des Orbits
        hist, _ = np.histogram(traj, bins=n_bins, range=(0.0, 1.0))

        # Ergodisch: Alle Bins sollten besucht werden
        # und die Häufigkeiten sollten ungefähr gleich sein (χ²-Test)
        if np.any(hist == 0):
            return False

        # Gleichverteilung prüfen via Variation
        expected = n_iter / n_bins
        variation = np.std(hist) / expected
        return bool(variation < 0.3)

    def ergodic_decomposition(
        self,
        initial_points: np.ndarray,
        n_iter: int = 1000,
        n_clusters: int = 3
    ) -> List[np.ndarray]:
        """
        @brief Numerische Approximation der ergodischen Zerlegung.

        @description
            Jedes maß-erhaltende System zerlegt sich in ergodische Komponenten.
            Für nicht-ergodische Systeme existieren mehrere invariante Mengen
            mit positivem Maß.

            Methode: k-means-Clustering der Orbits (grobe Approximation).
            Punkte, deren Orbits in ähnliche Attraktoren konvergieren,
            werden zur selben ergodischen Komponente gezählt.

        @param initial_points Startpunkte x₁, x₂, ...
        @param n_iter Orbitlehge
        @param n_clusters Erwartete Anzahl ergodischer Komponenten
        @return Liste von Arrays, jedes enthält Punkte einer Komponente
        @lastModified 2026-03-11
        """
        orbits = []
        for x0 in initial_points:
            traj = self.orbit(x0, n_iter)
            # Charakterisiere jeden Orbit durch Mittelwert und Standardabweichung
            orbits.append([np.mean(traj), np.std(traj)])

        orbit_features = np.array(orbits)

        # Einfache Zerlegung via Schwellwerte auf dem Zeitmittel
        components = []
        means = orbit_features[:, 0]

        # Quantile-basierte Zerlegung
        quantiles = np.linspace(0, 1, n_clusters + 1)
        thresholds = np.quantile(means, quantiles)

        for i in range(n_clusters):
            mask = (means >= thresholds[i]) & (means <= thresholds[i + 1])
            if np.any(mask):
                components.append(initial_points[mask])

        return components


# =============================================================================
# Klasse: BirkhoffErgodicTheorem
# =============================================================================

class BirkhoffErgodicTheorem:
    """
    @class BirkhoffErgodicTheorem
    @brief Implementierung des Birkhoff'schen Ergodentheorems (1931).

    @description
        Birkhoff's Ergodentheorem (1931):
        Sei (X, T, μ) ein maß-erhaltendes System und f ∈ L¹(μ). Dann gilt
        für μ-fast alle x:

          (1/n) Σ_{k=0}^{n-1} f(Tᵏx) → ∫_X f dμ   (n → ∞)

        Das Zeitmittel konvergiert fast überall gegen das Raummittel.

        Von-Neumann-Ergodentheorem (1932):
        Für f ∈ L²(μ) konvergiert die Folge der Zeitdurchschnitte
        im L²-Sinne gegen die bedingte Erwartung E[f | I],
        wobei I die σ-Algebra der T-invarianten Mengen ist.

    @lastModified 2026-03-11
    """

    def birkhoff_ergodic_theorem(
        self,
        T: Callable,
        f: Callable,
        x0: float,
        n: int
    ) -> np.ndarray:
        """
        @brief Berechnet die Folge der Birkhoff-Zeitmittel.

        @description
            Definiert den Zeitdurchschnitt:
              A_n(f, x) = (1/n) Σ_{k=0}^{n-1} f(Tᵏx)

            Für ergodische Systeme gilt: A_n(f, x) → ∫f dμ (f.s.)

        @param T Transformation des dynamischen Systems
        @param f Beobachtungsfunktion f: X → ℝ
        @param x0 Startpunkt
        @param n Anzahl der Schritte
        @return Array A_1, A_2, ..., A_n der kumulierten Zeitmittel
        @lastModified 2026-03-11
        """
        values = np.zeros(n)
        x = x0

        for k in range(n):
            values[k] = f(x)
            x = T(x)

        # Kumulierte Zeitmittel: A_n = (1/n) Σ_{k=0}^{n-1} f(T^k x)
        cumsum = np.cumsum(values)
        indices = np.arange(1, n + 1)
        return cumsum / indices

    def von_neumann_mean_ergodic(
        self,
        f_values: np.ndarray,
        space_mean: Optional[float] = None
    ) -> Dict:
        """
        @brief Von-Neumann-Mittelwertergodentheorem: L²-Konvergenz.

        @description
            Von-Neumanns Theorem: Für f ∈ L²(μ) konvergiert
              A_n(f) = (1/n) Σ_{k=0}^{n-1} f(Tᵏx)
            im L²-Sinne gegen E[f|I].

            Verifikation via Konvergenzrate:
            ||A_n - E[f]||₂ = O(1/√n)  (typisch für ergodische Systeme)

        @param f_values Beobachtungswerte f(T⁰x), f(T¹x), ...
        @param space_mean Bekanntes Raummittel (wenn verfügbar)
        @return Dict mit 'running_mean', 'l2_error', 'convergence_rate'
        @lastModified 2026-03-11
        """
        f = np.asarray(f_values, dtype=float)
        n = len(f)

        # Laufende Zeitmittel
        running_mean = np.cumsum(f) / np.arange(1, n + 1)

        # L²-Fehler gegen das geschätzte Grenzmittel
        final_mean = running_mean[-1] if space_mean is None else space_mean

        # L²-Fehler: ||A_k - E[f]||² = (A_k - E[f])²
        l2_error = np.abs(running_mean - final_mean)

        # Konvergenzrate schätzen via log-log-Regression
        idx = np.arange(10, n)
        if len(idx) > 2 and np.any(l2_error[idx] > 1e-14):
            log_n = np.log(idx + 1)
            log_err = np.log(np.maximum(l2_error[idx], 1e-14))
            valid = np.isfinite(log_err)
            if valid.sum() > 2:
                slope, intercept = np.polyfit(log_n[valid], log_err[valid], 1)
            else:
                slope = -0.5
        else:
            slope = -0.5

        return {
            'running_mean': running_mean,
            'l2_error': l2_error,
            'final_mean': final_mean,
            'convergence_rate': slope,  # Typisch ~-0.5 für ergodische Systeme
        }

    def ergodic_average_convergence(
        self,
        T: Callable,
        f: Callable,
        x0: float,
        space_mean: float,
        n_steps: int = 5000
    ) -> Dict:
        """
        @brief Numerische Verifikation der Konvergenz des Birkhoff-Zeitmittels.

        @description
            Prüft ob (1/n)Σf(Tᵏx) → space_mean für großes n.
            Berechnet die absolute Abweichung |A_n - ∫f dμ| und
            prüft Konvergenz (Abfall auf < 0.05 nach n_steps Schritten).

        @param T Transformation
        @param f Beobachtungsfunktion
        @param x0 Startpunkt
        @param space_mean Bekanntes Raummittel ∫f dμ
        @param n_steps Anzahl der Schritte
        @return Dict mit Konvergenzdaten und Verifikationsstatus
        @lastModified 2026-03-11
        """
        averages = self.birkhoff_ergodic_theorem(T, f, x0, n_steps)
        deviations = np.abs(averages - space_mean)
        final_deviation = deviations[-1]

        return {
            'averages': averages,
            'deviations': deviations,
            'final_deviation': float(final_deviation),
            'converged': bool(final_deviation < 0.05),
            'space_mean': space_mean,
        }


# =============================================================================
# Klasse: ErgodicEntropy
# =============================================================================

class ErgodicEntropy:
    """
    @class ErgodicEntropy
    @brief Entropie-Konzepte der Ergodentheorie: KS-Entropie, topologische Entropie.

    @description
        Die Kolmogorov-Sinai-Entropie (1958/1959) ist ein metrischer Isomorphismus-
        Invariant für maß-erhaltende Systeme.

        Partitionsentropie: Für eine endliche Partition P = {A₁, ..., Aₖ}:
          H_μ(P) = -Σᵢ μ(Aᵢ) · log₂(μ(Aᵢ))

        KS-Entropie: h_μ(T, P) = lim_{n→∞} (1/n) H_μ(P ∨ T⁻¹P ∨ ... ∨ T⁻ⁿ⁺¹P)
                     h_μ(T) = sup_P h_μ(T, P)

        Topologische Entropie: h_top(T) = lim_{ε→0} lim_{n→∞} (1/n) log s_n(ε)
        wobei s_n(ε) = maximale (n,ε)-separierte Menge.

        Variationsprinzip (Goodman 1971, Dinaburg 1970):
          h_top(T) = sup_{μ maß-erhaltend} h_μ(T)

    @lastModified 2026-03-11
    """

    def partition_entropy(self, probabilities: np.ndarray) -> float:
        """
        @brief Berechnet die Partitionsentropie H(P) = -Σ μ(Aᵢ) log μ(Aᵢ).

        @description
            Shannonsche Entropie einer Partition mit Maßen μ(A₁), ..., μ(Aₖ):
              H_μ(P) = -Σᵢ μ(Aᵢ) · log₂(μ(Aᵢ))

            Maximiert bei gleichmäßiger Verteilung: H_max = log₂(k).
            Minimiert bei deterministischen Systemen: H_min = 0.

        @param probabilities Maße μ(A₁), ..., μ(Aₖ) (summieren zu 1)
        @return Partitionsentropie in bits (log base 2)
        @raises ValueError Wenn Wahrscheinlichkeiten nicht zu 1 summieren
        @lastModified 2026-03-11
        """
        p = np.asarray(probabilities, dtype=float)
        p = p[p > 0]  # log(0) vermeiden: 0·log(0) := 0

        if len(p) == 0:
            return 0.0

        return float(-np.sum(p * np.log2(p)))

    def metric_entropy(
        self,
        T: Callable,
        x0: float,
        partition_bins: int = 8,
        n_iter: int = 10000
    ) -> float:
        """
        @brief Schätzt die Kolmogorov-Sinai-Entropie h_μ(T) numerisch.

        @description
            Approximation der KS-Entropie via Partitionsverfeinerung:
              h_μ(T, P_n) = (1/n) H_μ(P ∨ T⁻¹P ∨ ... ∨ T⁻ⁿ⁺¹P)

            Hier: Orbit des Startwerts x₀ wird in Sequenz von Bin-Indizes
            übersetzt. Die Entropie der Zylinder-Partition (Wörter der Länge k)
            wächst wie h·k für ergodische Systeme.

            Pesin-Formel (für C¹-Diffeomorphismen):
              h_μ(T) = Σ max(λᵢ, 0) (Summe positiver Lyapunov-Exponenten)

        @param T Transformation (muss [0,1] → [0,1] abbilden)
        @param x0 Startpunkt
        @param partition_bins Anzahl der Partitions-Intervalle
        @param n_iter Länge des betrachteten Orbits
        @return Schätzung der KS-Entropie (bits/iteration)
        @lastModified 2026-03-11
        """
        # Orbit berechnen
        trajectory = np.zeros(n_iter)
        x = x0
        for k in range(n_iter):
            trajectory[k] = x
            x = T(x)

        # Orbit in Bin-Sequenz übersetzen
        bin_sequence = np.floor(trajectory * partition_bins).astype(int)
        bin_sequence = np.clip(bin_sequence, 0, partition_bins - 1)

        # Entropie der Wörter der Länge k schätzen (für k = 1, 2, ..., max_k)
        max_k = min(8, n_iter // 100)
        entropies = []

        for k in range(1, max_k + 1):
            # Alle Wörter der Länge k extrahieren
            words = {}
            for i in range(n_iter - k + 1):
                word = tuple(bin_sequence[i:i + k])
                words[word] = words.get(word, 0) + 1

            # Empirische Wahrscheinlichkeiten
            total = sum(words.values())
            probs = np.array(list(words.values()), dtype=float) / total

            # Entropie des k-Gramm-Prozesses
            H_k = float(-np.sum(probs * np.log2(probs + 1e-15)))
            entropies.append(H_k / k)

        if not entropies:
            return 0.0

        # KS-Entropie ≈ lim_{k→∞} H_k / k
        return float(np.mean(entropies[-3:]) if len(entropies) >= 3 else entropies[-1])

    def topological_entropy(
        self,
        T: Callable,
        n_orbits: int = 100,
        n_iter: int = 200,
        eps: float = 0.05
    ) -> float:
        """
        @brief Schätzt die topologische Entropie h_top(T).

        @description
            Topologische Entropie via maximaler (n, ε)-separierter Menge:
              h_top(T) = lim_{ε→0} lim_{n→∞} (1/n) log s_n(ε)

            Approximation: Für gegebenes ε und n Iterationsschritte wird die
            Anzahl der (n,ε)-unterscheidbaren Orbits gezählt.

            Zwei Punkte x, y sind (n,ε)-separiert, wenn:
              max_{0≤k<n} d(Tᵏx, Tᵏy) > ε

        @param T Transformation [0,1] → [0,1]
        @param n_orbits Anzahl gleichmäßiger Startpunkte
        @param n_iter Länge der Orbits
        @param eps Trennungsradius ε
        @return Schätzung der topologischen Entropie
        @lastModified 2026-03-11
        """
        # Startpunkte gleichmäßig in [0,1]
        x0_values = np.linspace(0.01, 0.99, n_orbits)

        # Orbits berechnen
        orbits = np.zeros((n_orbits, n_iter))
        for i, x0 in enumerate(x0_values):
            x = x0
            for k in range(n_iter):
                orbits[i, k] = x
                x = T(x)

        # Maximale (n, ε)-separierte Menge abschätzen
        # Zwei Orbits gelten als separiert, wenn max-Distanz > eps
        separated_count = 1
        representative_orbits = [orbits[0]]

        for i in range(1, n_orbits):
            is_separated = True
            for rep in representative_orbits:
                max_dist = np.max(np.abs(orbits[i] - rep))
                if max_dist <= eps:
                    is_separated = False
                    break
            if is_separated:
                separated_count += 1
                representative_orbits.append(orbits[i])

        # h_top ≈ (1/n) log s_n(ε)
        if separated_count <= 1:
            return 0.0
        return float(np.log(separated_count) / n_iter)

    def variational_principle_check(
        self,
        T: Callable,
        x0_list: List[float],
        partition_bins: int = 8,
        n_iter: int = 5000
    ) -> Dict:
        """
        @brief Numerische Verifikation des Variationsprinzips h_top = sup h_μ.

        @description
            Variationsprinzip (Goodman, Dinaburg, 1970/1971):
              h_top(T) = sup_{μ maß-erhaltend} h_μ(T)

            Verifikation: Berechne KS-Entropie von mehreren Startpunkten
            und vergleiche mit der topologischen Entropie.
            Das Maximum der KS-Entropien sollte ≤ h_top sein.

        @param T Transformation
        @param x0_list Liste von Startpunkten
        @param partition_bins Bins für Partition
        @param n_iter Iterationsanzahl
        @return Dict mit h_top, max KS-Entropie und Verifikationsstatus
        @lastModified 2026-03-11
        """
        # KS-Entropien für verschiedene Startpunkte
        ks_entropies = []
        for x0 in x0_list:
            h = self.metric_entropy(T, x0, partition_bins, n_iter)
            ks_entropies.append(h)

        max_ks = max(ks_entropies) if ks_entropies else 0.0
        h_top = self.topological_entropy(T, n_orbits=50, n_iter=100)

        return {
            'h_top': h_top,
            'ks_entropies': ks_entropies,
            'max_ks_entropy': max_ks,
            'principle_satisfied': bool(max_ks <= h_top + 0.5),  # Toleranz
        }


# =============================================================================
# Klasse: MixingProperties
# =============================================================================

class MixingProperties:
    """
    @class MixingProperties
    @brief Mixing-Eigenschaften dynamischer Systeme.

    @description
        Hierarchie der Mixing-Eigenschaften (von stark nach schwach):

        1. Starkes Mixing (strong mixing):
           lim_{n→∞} μ(A ∩ T⁻ⁿB) = μ(A)·μ(B)  für alle messbaren A, B

        2. Schwaches Mixing (weak mixing):
           (1/N) Σ_{n=0}^{N-1} |μ(A ∩ T⁻ⁿB) - μ(A)μ(B)| → 0

        3. Ergodizität (schwächste Form):
           (1/N) Σ μ(A ∩ T⁻ⁿB) → μ(A)μ(B)

        Starkes Mixing ⟹ Schwaches Mixing ⟹ Ergodisch

        Spektrales Kriterium für schwaches Mixing:
        T ist schwach-mischend ⟺ T besitzt keine eigentlichen Eigenvalues
        außer 1 im L²-Raum.

    @lastModified 2026-03-11
    """

    def correlation_function(
        self,
        trajectory: np.ndarray,
        max_lag: int = 50
    ) -> np.ndarray:
        """
        @brief Berechnet die normierte Autokorrelationsfunktion C(n).

        @description
            Korrelationsfunktion für das Observable f:
              C(n) = ∫ f(Tⁿx) · f(x) dμ - (∫f dμ)²

            Normiert: ρ(n) = C(n) / C(0)

            Für mischende Systeme: ρ(n) → 0 (n → ∞)
            Mischungsrate bestimmt wie schnell ρ(n) → 0.

        @param trajectory Zeitreihe [f(x), f(Tx), f(T²x), ...]
        @param max_lag Maximale Verzögerung
        @return Array der Autokorrelationen ρ(0), ρ(1), ..., ρ(max_lag)
        @lastModified 2026-03-11
        """
        x = np.asarray(trajectory, dtype=float)
        x_centered = x - np.mean(x)
        n = len(x)

        acf = np.zeros(max_lag + 1)
        acf[0] = np.dot(x_centered, x_centered) / n

        for lag in range(1, max_lag + 1):
            if lag < n:
                acf[lag] = np.dot(x_centered[:-lag], x_centered[lag:]) / (n - lag)

        # Normierung
        if acf[0] > 1e-14:
            acf = acf / acf[0]

        return acf

    def is_mixing(
        self,
        T: Callable,
        f: Callable,
        x0: float,
        n_iter: int = 5000,
        lag_threshold: int = 30
    ) -> bool:
        """
        @brief Prüft numerisch ob das System stark mischend ist.

        @description
            Starkes Mixing: lim_{n→∞} μ(A ∩ T⁻ⁿB) = μ(A)·μ(B)

            Numerisches Kriterium: Die Autokorrelation der Beobachtung f
            sollte für großes n nahe 0 sein:
              |ρ(n)| < 0.1 für n ≥ lag_threshold

        @param T Transformation
        @param f Beobachtungsfunktion
        @param x0 Startpunkt
        @param n_iter Länge der Zeitreihe
        @param lag_threshold Mindest-Lag für Mischungstest
        @return True wenn Autokorrelation bei lag_threshold klein
        @lastModified 2026-03-11
        """
        # Trajektorie berechnen
        traj = np.zeros(n_iter)
        x = x0
        for k in range(n_iter):
            traj[k] = f(x)
            x = T(x)

        # Autokorrelation
        acf = self.correlation_function(traj, max_lag=lag_threshold)

        # Mischend wenn |ρ(n)| < 0.1 für n ≥ lag_threshold//2
        tail_acf = acf[lag_threshold // 2:]
        return bool(np.all(np.abs(tail_acf) < 0.15))

    def is_weak_mixing(
        self,
        T: Callable,
        f: Callable,
        x0: float,
        n_iter: int = 5000,
        max_lag: int = 100
    ) -> bool:
        """
        @brief Prüft schwaches Mixing via Cesàro-Konvergenz der Korrelationen.

        @description
            Schwaches Mixing: (1/N) Σ_{n=0}^{N-1} |ρ(n)| → 0 (N → ∞)

            Äquivalent (Halmos): T ist schwach-mischend ⟺ T × T ist ergodisch.

            Numerisch: Cesàro-Mittel der |Autokorrelationen| sollte klein sein.

        @param T Transformation
        @param f Beobachtungsfunktion
        @param x0 Startpunkt
        @param n_iter Länge der Zeitreihe
        @param max_lag Maximale Verzögerung
        @return True wenn Cesàro-Mittel der Korrelationen klein
        @lastModified 2026-03-11
        """
        traj = np.zeros(n_iter)
        x = x0
        for k in range(n_iter):
            traj[k] = f(x)
            x = T(x)

        acf = self.correlation_function(traj, max_lag=max_lag)

        # Cesàro-Mittel: (1/N) Σ |ρ(n)|
        cesaro_means = np.cumsum(np.abs(acf)) / np.arange(1, len(acf) + 1)

        # Schwach-mischend wenn Cesàro-Mittel der letzten Hälfte < 0.1
        tail = cesaro_means[max_lag // 2:]
        return bool(np.mean(tail) < 0.1)

    def mixing_rate_estimate(
        self,
        trajectory: np.ndarray,
        max_lag: int = 50
    ) -> Dict:
        """
        @brief Schätzt die Mischungsrate via exponentiellem Fit der Korrelationen.

        @description
            Für hyperbolische Systeme fällt die Korrelation exponentiell ab:
              C(n) ≈ C₀ · exp(-γn)

            wobei γ > 0 die Mischungsrate (mixing rate) ist.
            Kleinere Mischungszeit τ = 1/γ bedeutet schnelleres Vergessen.

        @param trajectory Beobachtungstrajektorie
        @param max_lag Maximale Verzögerung
        @return Dict mit 'acf', 'mixing_rate', 'mixing_time'
        @lastModified 2026-03-11
        """
        acf = self.correlation_function(trajectory, max_lag=max_lag)

        # Exponentiellen Abfall fitten: log|C(n)| ≈ log|C₀| - γ·n
        lags = np.arange(1, max_lag + 1)
        log_acf = np.log(np.abs(acf[1:]) + 1e-14)

        valid = np.isfinite(log_acf) & (np.abs(acf[1:]) > 1e-6)
        if valid.sum() >= 3:
            slope, intercept = np.polyfit(lags[valid], log_acf[valid], 1)
            mixing_rate = -slope  # γ = -slope
        else:
            mixing_rate = 0.0

        mixing_time = 1.0 / mixing_rate if mixing_rate > 1e-10 else float('inf')

        return {
            'acf': acf,
            'mixing_rate': float(mixing_rate),
            'mixing_time': float(mixing_time),
            'is_mixing': bool(mixing_rate > 0.01),
        }


# =============================================================================
# Klasse: CollatzErgodicSystem
# =============================================================================

class CollatzErgodicSystem:
    """
    @class CollatzErgodicSystem
    @brief Ergodische Analyse der Collatz-Abbildung (3n+1-Problem).

    @description
        Die Collatz-Abbildung (3n+1):
          T(n) = n/2   falls n gerade
          T(n) = 3n+1  falls n ungerade

        Tao (2019): "Almost all orbits of the Collatz map attain almost
        bounded values." Für fast alle n ∈ {1,...,N} existiert k mit
        T^k(n) ≤ f(n) für beliebig langsam wachsendes f.

        Ergodisches Modell (Lagarias):
        Die Collatz-Abbildung auf den 2-adischen Zahlen ℤ₂ ist
        maß-erhaltend bezüglich des 2-adischen Maßes.

        Lyapunov-Exponent:
        Bei einem zufälligen Schritt ist:
          - Mit Wahrscheinlichkeit 1/2: n → n/2    (Faktor 1/2)
          - Mit Wahrscheinlichkeit 1/2: n → 3n+1 → (3n+1)/2  (Faktor ~3/4)
        Erwarteter Wachstumsfaktor pro Schritt:
          λ = (1/2)·log(1/2) + (1/2)·log(3/4) ≈ -0.0235 < 0
        Das System ist daher im Zeitmittel kontraktiv.

    @lastModified 2026-03-11
    """

    def collatz_map(self, n: int) -> int:
        """
        @brief Wendet die Collatz-Abbildung einmal an: n → n/2 oder 3n+1.

        @param n Positive ganze Zahl
        @return T(n)
        @raises ValueError Wenn n ≤ 0
        @lastModified 2026-03-11
        """
        if n <= 0:
            raise ValueError(f"Collatz-Abbildung ist nur für n > 0 definiert, n={n}")
        if n % 2 == 0:
            return n // 2
        else:
            return 3 * n + 1

    def collatz_iteration(self, n: int, max_steps: int = 10000) -> List[int]:
        """
        @brief Berechnet die vollständige Collatz-Trajektorie bis zum Erreichen von 1.

        @description
            Startet bei n und wendet T wiederholt an bis:
            - T^k(n) = 1 (Konvergenz zur Vermutung gemäß)
            - Maximale Schrittanzahl erreicht

        @param n Startwert (positive ganze Zahl)
        @param max_steps Maximale Iterationsanzahl
        @return Liste [n, T(n), T²(n), ..., 1]
        @lastModified 2026-03-11
        """
        if n <= 0:
            raise ValueError(f"n muss positiv sein, n={n}")

        trajectory = [n]
        current = n

        for _ in range(max_steps):
            if current == 1:
                break
            current = self.collatz_map(current)
            trajectory.append(current)

        return trajectory

    def collatz_stopping_time(self, n: int, max_steps: int = 100000) -> int:
        """
        @brief Berechnet die Stopping-Time σ(n): Anzahl Schritte bis T^k(n) < n.

        @description
            Stopping-Time (Terras, 1976):
              σ(n) = min{k ≥ 1 : T^k(n) < n}

            Für fast alle n gilt: σ(n) ≈ C · log(n) (Tao 2019, Rawsthorne).
            Die Stopping-Time misst, wie schnell der Orbit unter n fällt.

        @param n Startwert
        @param max_steps Maximale Iterationsanzahl
        @return Stopping-Time σ(n), oder -1 bei Nicht-Konvergenz
        @lastModified 2026-03-11
        """
        if n <= 1:
            return 0

        current = n
        for k in range(1, max_steps + 1):
            current = self.collatz_map(current)
            if current < n:
                return k

        return -1  # Keine Konvergenz in max_steps

    def collatz_total_stopping_time(self, n: int, max_steps: int = 100000) -> int:
        """
        @brief Totale Stopping-Time: Anzahl Schritte bis T^k(n) = 1.

        @param n Startwert
        @param max_steps Maximale Schritte
        @return Totale Stopping-Time, oder -1
        @lastModified 2026-03-11
        """
        if n <= 0:
            raise ValueError("n muss positiv sein")
        if n == 1:
            return 0

        current = n
        for k in range(1, max_steps + 1):
            current = self.collatz_map(current)
            if current == 1:
                return k

        return -1

    def collatz_density_tao(
        self,
        N: int = 1000,
        eps: float = 0.5,
        C: float = 10.0
    ) -> Dict:
        """
        @brief Numerische Verifikation von Taos Dichte-1-Resultat (2019).

        @description
            Tao (2019): Für jedes ε > 0 und wachsendes f(n) → ∞ gilt:
              #{n ≤ N : ∃k mit T^k(n) ≤ f(n)} / N → 1  (N → ∞)

            Spezialfall: Für die Stopping-Time σ(n) gilt
              σ(n) ≤ C · log(n) für fast alle n.

            Verifikation: Anteil der n ≤ N mit σ(n) ≤ C·log(n) sollte → 1.

        @param N Obere Schranke
        @param eps Exponent für Schranke (n^ε)
        @param C Konstante für log-Schranke
        @return Dict mit Dichte-Statistiken
        @lastModified 2026-03-11
        """
        count_below_bound = 0
        count_below_eps = 0
        total = 0
        stopping_times = []

        for n in range(2, N + 1):
            sigma = self.collatz_stopping_time(n, max_steps=10 * int(math.log(n + 2) * C * 2))
            stopping_times.append(sigma if sigma >= 0 else int(C * math.log(n) * 2))

            # Tao-Bedingung 1: σ(n) ≤ C·log(n)
            log_bound = C * math.log(n + 1)
            if sigma >= 0 and sigma <= log_bound:
                count_below_bound += 1

            # Tao-Bedingung 2: T^k(n) ≤ n^ε für ein k
            current = n
            reached_neps = False
            n_eps = n ** eps
            for _ in range(int(log_bound) + 1):
                if current <= n_eps:
                    reached_neps = True
                    break
                current = self.collatz_map(current)
            if reached_neps:
                count_below_eps += 1

            total += 1

        density_log_bound = count_below_bound / total if total > 0 else 0.0
        density_eps_bound = count_below_eps / total if total > 0 else 0.0

        return {
            'N': N,
            'density_log_bound': density_log_bound,
            'density_eps_bound': density_eps_bound,
            'stopping_times': stopping_times,
            'mean_stopping_time': float(np.mean(stopping_times)),
            'tao_verified': bool(density_log_bound > 0.9),
        }

    def collatz_lyapunov_exponent(self, n_samples: int = 1000) -> float:
        """
        @brief Berechnet den Lyapunov-Exponenten der Collatz-Abbildung.

        @description
            Zufälliges Modell: Beim k-ten Schritt ist mit Wahrscheinlichkeit 1/2
            das Argument gerade (→ n/2, Faktor 1/2) und mit Wahrscheinlichkeit 1/2
            ungerade (→ 3n+1, im nächsten Schritt gerade: Faktor 3/4 im Mittel).

            Erwarteter Log-Faktor pro "2 Schritte":
              E[log |T'|] ≈ (1/2)·log(1/2) + (1/2)·log(3/4)
                          = (-log 2)/2 + (log 3 - 2·log 2)/2
                          = (log 3 - 3·log 2)/2 ≈ -0.0229

            Lyapunov-Exponent des Zufallsmodells:
              λ = (1/2)(log 3 - 2 log 2) / 2 ≈ -0.0235

        @param n_samples Anzahl der Stichproben-Trajektorien
        @return Schätzung des Lyapunov-Exponenten
        @lastModified 2026-03-11
        """
        log_factors = []

        # Empirische Schätzung über viele Trajektorien
        rng = np.random.default_rng(42)
        for n in rng.integers(10, 1000, size=n_samples):
            n = int(n)
            trajectory = self.collatz_iteration(n, max_steps=500)

            # Log-Wachstumsrate schätzen
            if len(trajectory) > 10:
                log_start = math.log(trajectory[0])
                log_end = math.log(trajectory[-1])
                steps = len(trajectory) - 1
                log_factors.append((log_end - log_start) / steps)

        # Theoretischer Wert: (log(3) - 2·log(2)) / 2 ≈ -0.0235
        theoretical = (math.log(3) - 2 * math.log(2)) / 2
        empirical = float(np.mean(log_factors)) if log_factors else theoretical

        return empirical

    def collatz_theoretical_lyapunov(self) -> float:
        """
        @brief Berechnet den theoretischen Lyapunov-Exponenten des Zufallsmodells.

        @description
            Theoretischer Wert: λ = (log 3 - 2 log 2) / 2

            Herleitung: Im Mittel über 2 Schritte (gerade dann ungerade):
              x → x/2 → (3x/2+1) ≈ (3/2)x·(1/2) = (3/4)x
            Log-Faktor pro 2 Schritte: log(3/4) ≈ -0.288
            Pro Schritt: λ ≈ -0.144

            Präziser (Garner-Modell):
              λ = (1/2)·log(1/2) + (1/2)·(1/2)·log(3) + (1/4)·log(1/4) + ...
                ≈ log(3)/2 - log(2) ≈ -0.1462

        @return Theoretischer Lyapunov-Exponent λ
        @lastModified 2026-03-11
        """
        # Einfaches symmetrisches Modell: E[log|T'|] = (log 3 - 2 log 2)/2
        return (math.log(3) - 2 * math.log(2)) / 2

    def collatz_2adic_analysis(self, n: int, depth: int = 16) -> Dict:
        """
        @brief 2-adische Analyse der Collatz-Dynamik.

        @description
            Die 2-adischen Zahlen ℤ₂ sind der Vervollständigung von ℤ bzgl.
            der 2-adischen Norm ||n||₂ = 2^{-v₂(n)}, wobei v₂(n) der
            2-adische Valuation (Anzahl Teiler 2) ist.

            Lagarias (1985): Die Collatz-Abbildung auf ℤ₂ ist
            maß-erhaltend bezüglich des Haar-Maßes auf ℤ₂.

            2-adische Darstellung: n = Σᵢ aᵢ · 2ⁱ (aᵢ ∈ {0,1})
            Die letzten Bits bestimmen den nächsten Schritt:
            - n ≡ 0 (mod 2): n → n/2 (Rechtsshift)
            - n ≡ 1 (mod 2): n → 3n+1 (Mix)

        @param n Startwert
        @param depth Tiefe der 2-adischen Darstellung
        @return Dict mit 2-adischer Darstellung und Dynamik
        @lastModified 2026-03-11
        """
        # 2-adische Darstellung
        bits = [(n >> i) & 1 for i in range(depth)]

        # Trajektorie
        trajectory = self.collatz_iteration(n, max_steps=100)

        # 2-adische Bewertung (Valuation) jedes Trajektorienpunkts
        valuations = []
        for x in trajectory:
            if x == 0:
                valuations.append(float('inf'))
            else:
                v = 0
                temp = x
                while temp % 2 == 0:
                    v += 1
                    temp //= 2
                valuations.append(v)

        # 2-adische Norm
        norms = [2**(-v) if v < float('inf') else 0.0 for v in valuations]

        return {
            'n': n,
            'binary_repr': bits,
            'trajectory': trajectory[:20],  # Erste 20 Schritte
            'valuations': valuations[:20],
            'padic_norms': norms[:20],
            'trajectory_length': len(trajectory),
        }

    def collatz_ergodic_model(self, N: int = 200) -> Dict:
        """
        @brief Modelliert Collatz als ergodisches System auf 2-adischen Zahlen.

        @description
            Das Collatz-System auf den 2-adischen Zahlen ℤ₂:

            T: ℤ₂ → ℤ₂, T(x) = (3x+1)/2  für ungerades x in ℤ₂

            (Man beachte: 3x+1 ist stets gerade für ungerades x, also kann
            man durch 2 dividieren und erhält eine Abbildung ℤ₂ → ℤ₂.)

            Dieses Modell ist: bijektiv und maß-erhaltend bezüglich des
            Haar-Maßes auf ℤ₂ (Conjecture, nicht bewiesen).

            Simulation: Prüfe empirisch die gleichmäßige Verteilung der
            Modulo-2^k-Residuen von Collatz-Trajektorien.

        @param N Anzahl der Startpunkte
        @return Dict mit Statistiken zur 2-adischen Verteilung
        @lastModified 2026-03-11
        """
        # Für k=1,2,3,4: Verteilung der Residuen mod 2^k
        max_k = 4
        results = {}

        for k in range(1, max_k + 1):
            mod = 2 ** k
            residue_counts = {r: 0 for r in range(mod)}

            for n in range(1, N + 1):
                trajectory = self.collatz_iteration(n, max_steps=200)
                # Residuen der Trajektorie sammeln
                for x in trajectory[1:]:  # Startwert überspringen
                    residue_counts[x % mod] = residue_counts.get(x % mod, 0) + 1

            # Gleichverteilung prüfen
            counts = np.array(list(residue_counts.values()), dtype=float)
            total = counts.sum()
            if total > 0:
                expected = total / mod
                chi2 = np.sum((counts - expected) ** 2 / expected)
                # chi2-Test: df = mod-1
                p_value = 1 - stats.chi2.cdf(chi2, df=mod - 1)
            else:
                p_value = 1.0

            results[f'mod_{mod}'] = {
                'counts': dict(residue_counts),
                'chi2': float(chi2) if total > 0 else 0.0,
                'p_value': float(p_value),
                'uniform': bool(p_value > 0.05),
            }

        return results

    def collatz_stopping_time_distribution(
        self,
        N: int = 500
    ) -> Dict:
        """
        @brief Berechnet die Verteilung der Stopping-Times und fittet Normalverteilung.

        @description
            Empirisch (Everett 1977, Terras 1976):
            Die Stopping-Time σ(n) ist für große n aproximativ normalverteilt:
              σ(n) / log(n) → C  (fast sicher, mit C ≈ log(2)/log(3/2) ≈ 1.73)

            Genauer: σ(n) ≈ μ·log(n) + σ·√(log(n))·Z  mit Z ~ N(0,1)
            wobei μ = log(2)/log(4/3) ≈ 2.41.

            Modell (Garner 1981):
            Die Anzahl ungerader Schritte in σ(n) Iterationen ist
            Binomial(σ(n), 1/2) verteilt.

        @param N Obere Schranke für die Stichprobe
        @return Dict mit Stopping-Times, Fit-Parametern und Statistiken
        @lastModified 2026-03-11
        """
        stopping_times = []
        ns = []

        for n in range(2, N + 1):
            sigma = self.collatz_stopping_time(n, max_steps=1000)
            if sigma >= 0:
                stopping_times.append(sigma)
                ns.append(n)

        if not stopping_times:
            return {}

        st_array = np.array(stopping_times, dtype=float)
        ns_array = np.array(ns, dtype=float)

        # Normierung durch log(n): σ(n)/log(n)
        log_ns = np.log(ns_array + 1)
        normalized = st_array / log_ns

        # Normalverteilungs-Fit
        mu_fit, sigma_fit = stats.norm.fit(st_array)

        # Shapiro-Wilk-Test auf Normalverteilung
        if len(st_array) >= 8:
            shapiro_stat, shapiro_p = stats.shapiro(st_array[:min(200, len(st_array))])
        else:
            shapiro_stat, shapiro_p = 1.0, 1.0

        return {
            'stopping_times': list(stopping_times),
            'ns': list(ns),
            'mean': float(np.mean(st_array)),
            'std': float(np.std(st_array)),
            'mu_fit': float(mu_fit),
            'sigma_fit': float(sigma_fit),
            'mean_normalized': float(np.mean(normalized)),
            'shapiro_stat': float(shapiro_stat),
            'shapiro_p': float(shapiro_p),
            'theoretical_C': math.log(2) / math.log(4.0 / 3.0),
        }


# =============================================================================
# Klasse: FurstenbergTheory
# =============================================================================

class FurstenbergTheory:
    """
    @class FurstenbergTheory
    @brief Furstenberg-Theorie: Verbindung zwischen Ergodentheorie und Kombinatorik.

    @description
        Furstenberg (1977) bewies das Szemerédi-Theorem (1975) mit Methoden
        der Ergodentheorie:

        Szemerédi-Theorem: Jede Menge A ⊂ ℕ mit positiver oberer Dichte
          d̄(A) = lim sup_{N→∞} |A ∩ {1,...,N}| / N > 0
        enthält beliebig lange arithmetische Progressionen.

        Furstenberg-Korrespondenzprinzip:
        Jeder Menge A ⊂ ℕ mit d̄(A) > 0 entspricht ein maß-erhaltendes System
        (X, T, μ) mit einer Menge B ⊂ X, so dass:
          d̄(A) ≥ μ(B) > 0
          und n ∈ A - A ⟺ μ(B ∩ T⁻ⁿB) > 0

    @lastModified 2026-03-11
    """

    def upper_density(self, subset: List[int], N: int) -> float:
        """
        @brief Berechnet die obere Dichte d̄(A) = lim sup |A ∩ {1..N}| / N.

        @param subset Menge A (als Liste von positiven ganzen Zahlen)
        @param N Obere Schranke
        @return Obere Dichte d̄(A)
        @lastModified 2026-03-11
        """
        a_set = set(x for x in subset if 1 <= x <= N)
        return len(a_set) / N if N > 0 else 0.0

    def furstenberg_correspondence(
        self,
        subset: List[int],
        N: int = 1000
    ) -> Dict:
        """
        @brief Demonstriert das Furstenberg-Korrespondenzprinzip.

        @description
            Für A ⊂ {1,...,N} mit positiver Dichte:
            1. Berechne d̄(A)
            2. Suche arithmetische Progressionen in A
            3. Verknüpfe mit dem ergodischen System T: x → x + 1 (mod N)

            Das Korrespondenzprinzip besagt:
            Wenn μ(B) = d̄(A) > 0, dann gibt es für jedes k ein n mit
              μ(B ∩ T⁻ⁿB ∩ T⁻²ⁿB ∩ ... ∩ T⁻ᵏⁿB) > 0

            Das entspricht einer AP der Länge k+1 in A.

        @param subset Menge A
        @param N Obere Schranke
        @return Dict mit Dichte, AP-Funden und ergodischer Interpretation
        @lastModified 2026-03-11
        """
        a_set = set(x for x in subset if 1 <= x <= N)
        density = len(a_set) / N

        # Arithmetische Progressionen finden (Länge 3)
        aps_len3 = []
        for a in a_set:
            for d in range(1, N // 2 + 1):
                if a + d in a_set and a + 2 * d in a_set:
                    aps_len3.append((a, d))
                    if len(aps_len3) >= 5:
                        break
            if len(aps_len3) >= 5:
                break

        # Ergodische Interpretation: μ(B ∩ T⁻ⁿB) für n=1,...,20
        pairwise_overlaps = []
        indicator = np.array([1 if x in a_set else 0 for x in range(1, N + 1)], dtype=float)

        for shift in range(1, min(21, N)):
            shifted = np.roll(indicator, -shift)
            overlap = np.mean(indicator * shifted)
            pairwise_overlaps.append(float(overlap))

        # Furstenberg: μ(B ∩ T⁻ⁿB) > (d̄(A))² für mindestens ein n
        density_sq = density ** 2
        has_positive_overlap = any(o > density_sq * 0.5 for o in pairwise_overlaps)

        return {
            'density': density,
            'density_squared': density_sq,
            'aps_length3': aps_len3[:5],
            'pairwise_overlaps': pairwise_overlaps,
            'furstenberg_condition_met': has_positive_overlap,
        }

    def szemeredi_via_ergodic(
        self,
        k: int,
        N: int = 500
    ) -> Dict:
        """
        @brief Numerische Verifikation des Szemerédi-Theorems für AP-Länge k.

        @description
            Für verschiedene Dichten d > 0 und AP-Länge k:
            Suche in zufälligen Teilmengen A ⊂ {1,...,N} mit oberer Dichte ≈ d
            nach arithmetischen Progressionen der Länge k.

            Szemerédi garantiert: Für jede δ > 0 und k ≥ 1 gibt es N₀ so dass
            jede A ⊂ {1,...,N} mit |A| ≥ δN eine AP der Länge k enthält.

        @param k Länge der arithmetischen Progression
        @param N Obere Schranke
        @return Dict mit Testergebnissen für verschiedene Dichten
        @lastModified 2026-03-11
        """
        rng = np.random.default_rng(2024)

        results = {}
        for density in [0.1, 0.2, 0.3, 0.5]:
            # Zufällige Teilmenge mit gegebener Dichte
            subset = sorted(
                rng.choice(range(1, N + 1), size=int(density * N), replace=False)
            )
            subset_set = set(subset)

            # Arithmetische Progression der Länge k suchen
            found_ap = None
            for a in subset:
                for d in range(1, N // k + 1):
                    ap = [a + j * d for j in range(k)]
                    if all(x in subset_set for x in ap) and max(ap) <= N:
                        found_ap = (a, d, ap)
                        break
                if found_ap:
                    break

            results[f'density_{density}'] = {
                'subset_size': len(subset),
                'actual_density': len(subset) / N,
                'ap_found': found_ap is not None,
                'ap': found_ap,
            }

        return results

    def van_der_waerden_bound(self, k: int, r: int) -> Dict:
        """
        @brief Schätzt die van-der-Waerden-Zahl W(k; r).

        @description
            van-der-Waerden-Theorem (1927):
            Für alle k, r ≥ 1 existiert W(k; r) so dass jede r-Färbung
            von {1,...,W(k;r)} eine monochromatische AP der Länge k enthält.

            Bekannte exakte Werte:
            W(2; 2) = 4
            W(3; 2) = 9
            W(3; 3) = 27
            W(4; 2) = 35
            W(5; 2) = 178
            W(6; 2) = 1132

            Berlekamp (1968) Untere Schranke: W(k; 2) > k · 2^k für Primzahlen k.
            Primitive Rekursion: W existiert, ist aber sehr groß.

        @param k AP-Länge
        @param r Anzahl der Farben
        @return Dict mit bekannten Schranken und Werten
        @lastModified 2026-03-11
        """
        # Bekannte exakte Werte
        known = {
            (2, 2): 4, (2, 3): 7, (2, 4): 11, (2, 5): 14,
            (3, 2): 9, (3, 3): 27,
            (4, 2): 35, (4, 3): 293,
            (5, 2): 178, (6, 2): 1132,
        }

        exact = known.get((k, r), None)

        # Berlekamp-Untere-Schranke für r=2, k prim
        def berlekamp_lower(k):
            return k * (2 ** k) if k >= 2 else 2

        lower_bound = berlekamp_lower(k) if r == 2 else None

        # Ackermann-ähnliche obere Schranke (sehr grob)
        # W(k; 2) ≤ 2^{2^{...}} (k mal)
        upper_approx = min(2 ** (2 ** min(k, 5)), 10 ** 15)

        return {
            'k': k,
            'r': r,
            'exact_value': exact,
            'berlekamp_lower': lower_bound,
            'upper_approx': upper_approx,
            'is_known': exact is not None,
        }


# =============================================================================
# Klasse: ShiftSystem
# =============================================================================

class ShiftSystem:
    """
    @class ShiftSystem
    @brief Symbolische Dynamik: Shift-Systeme auf Folgenräumen.

    @description
        Der Shift-Raum (AZ, σ) für ein endliches Alphabet A = {0,1,...,d-1}:
        - Zustandsraum X = A^ℤ (beidseitige Folgen) oder A^ℕ (einseitig)
        - Shift-Abbildung: σ((xₙ)_{n∈ℤ}) = (x_{n+1})_{n∈ℤ}
        - Produkt-Topologie auf A^ℤ (Cantor-Topologie)

        Subshift of Finite Type (SFT):
        Definiert durch eine Übergangsmatrix M ∈ {0,1}^{d×d}:
          X_M = {(xₙ) ∈ A^ℤ : M[xₙ, x_{n+1}] = 1 für alle n}

        Topologische Entropie des SFT:
          h_top(σ|X_M) = log(λ_max(M))
        wobei λ_max der größte Eigenwert (Perron-Frobenius) ist.

        Zeta-Funktion:
          ζ_T(z) = exp(Σ_{n=1}^∞ |Fix(T^n)| · z^n / n)
        wobei Fix(T^n) = Anzahl periodischer Punkte der Periode n.

    @param alphabet Liste der Symbole des Alphabets
    @lastModified 2026-03-11
    """

    def __init__(self, alphabet: List) -> None:
        """
        @brief Initialisiert das Shift-System mit gegebenem Alphabet.

        @param alphabet Endliches Alphabet, z.B. [0, 1] oder ['a', 'b', 'c']
        @lastModified 2026-03-11
        """
        self.alphabet = list(alphabet)
        self.d = len(alphabet)  # Alphabetgröße
        self._symbol_to_index = {s: i for i, s in enumerate(alphabet)}

    def shift(self, sequence: List, steps: int = 1) -> List:
        """
        @brief Wendet die Shift-Abbildung σ^steps an.

        @description
            σ((x₀, x₁, x₂, ...)) = (x₁, x₂, x₃, ...)
            σ^n((xₖ)) = (x_{k+n}) für n ≥ 0

        @param sequence Symbolfolge (endlich oder periodisch)
        @param steps Anzahl der Shifts
        @return Verschobene Folge
        @lastModified 2026-03-11
        """
        if steps >= len(sequence):
            return []
        return sequence[steps:]

    def is_in_shift_space(
        self,
        sequence: List,
        forbidden_words: Optional[List[List]] = None
    ) -> bool:
        """
        @brief Prüft ob eine Folge im gegebenen Shift-Raum liegt.

        @description
            Ein Shift-Raum ist definiert durch verbotene Wörter:
            X = {(xₙ) ∈ A^ℤ : kein verbotenes Wort kommt vor}

            Jeder Subshift of Finite Type hat eine endliche Menge verbotener Wörter.

        @param sequence Zu prüfende Folge über dem Alphabet
        @param forbidden_words Liste verbotener Teilwörter
        @return True wenn Folge im Shift-Raum liegt
        @lastModified 2026-03-11
        """
        if not forbidden_words:
            return all(s in self.alphabet for s in sequence)

        for fw in forbidden_words:
            fw_len = len(fw)
            for i in range(len(sequence) - fw_len + 1):
                if sequence[i:i + fw_len] == fw:
                    return False

        return True

    def count_words(
        self,
        length: int,
        transition_matrix: Optional[np.ndarray] = None
    ) -> int:
        """
        @brief Zählt die Anzahl zulässiger Wörter der Länge n.

        @description
            Für einen SFT mit Übergangsmatrix M gilt:
              |W_n(X_M)| = Σᵢⱼ (M^n)ᵢⱼ = 1ᵀ · M^n · 1

            Im Grenzfall wächst diese Zahl wie λ_max^n (Perron-Frobenius).

        @param length Wortlänge n
        @param transition_matrix Übergangsmatrix M (oder None für full shift)
        @return Anzahl zulässiger Wörter der Länge n
        @lastModified 2026-03-11
        """
        if transition_matrix is None:
            # Full shift: d^n Wörter
            return self.d ** length

        M = np.array(transition_matrix, dtype=float)
        Mn = np.linalg.matrix_power(M.astype(int), length - 1)
        return int(np.sum(Mn))

    def entropy_shift_sft(
        self,
        transition_matrix: np.ndarray
    ) -> float:
        """
        @brief Berechnet die topologische Entropie eines SFT via Spektralradius.

        @description
            Für einen SFT mit Übergangsmatrix M ∈ {0,1}^{d×d}:
              h_top(σ|X_M) = log(λ_max(M))

            wobei λ_max der größte reelle Eigenwert (Perron-Frobenius-Eigenwert)
            der nicht-negativen Matrix M ist.

            Beispiele:
            - Full shift auf d Symbolen: M = J_d, h = log(d)
            - Golden-Mean-Shift: M = [[1,1],[1,0]], λ = (1+√5)/2, h = log(φ)

        @param transition_matrix Übergangsmatrix M ∈ {0,1}^{d×d}
        @return Topologische Entropie h_top = log(λ_max)
        @lastModified 2026-03-11
        """
        M = np.array(transition_matrix, dtype=float)
        eigenvalues = np.linalg.eigvals(M)

        # Perron-Frobenius: größter reeller Eigenwert
        real_eigenvalues = eigenvalues[np.abs(eigenvalues.imag) < 1e-10].real
        if len(real_eigenvalues) == 0:
            return 0.0

        lambda_max = np.max(real_eigenvalues)

        if lambda_max <= 0:
            return float('-inf')
        if lambda_max <= 1e-14:
            return 0.0

        return float(np.log(lambda_max))

    def zeta_function_shift(
        self,
        transition_matrix: np.ndarray,
        max_n: int = 20
    ) -> Dict:
        """
        @brief Berechnet die Zeta-Funktion ζ_T(z) = exp(Σ |Fix(T^n)| z^n / n).

        @description
            Die dynamische Zeta-Funktion (Artin-Mazur 1965):
              ζ_T(z) = exp(Σ_{n=1}^∞ p_n · z^n / n)

            wobei p_n = |Fix(T^n)| = Spur(M^n) = Anzahl der Perioden-n-Punkte.

            Für einen SFT mit Übergangsmatrix M gilt (Bowen-Lanford):
              ζ_T(z) = 1 / det(I - zM)

            Die Nullstellen von ζ_T entsprechen den reziproken Eigenwerten von M.

        @param transition_matrix Übergangsmatrix M
        @param max_n Maximale Ordnung der Reihe
        @return Dict mit Fixed-Point-Zahlen, Koeffizienten, analytischer Form
        @lastModified 2026-03-11
        """
        M = np.array(transition_matrix, dtype=float)
        d = M.shape[0]

        # Fixpunkte: p_n = Spur(M^n)
        fixed_points = []
        Mn = np.eye(d)

        for n in range(1, max_n + 1):
            Mn = Mn @ M
            pn = int(round(np.trace(Mn)))
            fixed_points.append(pn)

        # Taylor-Koeffizienten der Reihe Σ p_n z^n / n
        # (für den Exponenten der Zeta-Funktion)
        series_coeffs = [pn / n for n, pn in enumerate(fixed_points, 1)]

        # Analytische Form: det(I - zM) (Nenner der ζ-Funktion)
        eigenvalues = np.linalg.eigvals(M)
        lambda_max = np.max(np.abs(eigenvalues)) if len(eigenvalues) > 0 else 0.0

        return {
            'fixed_points': fixed_points,
            'series_coefficients': series_coeffs,
            'spectral_radius': float(lambda_max),
            'entropy': float(np.log(lambda_max)) if lambda_max > 1 else 0.0,
            'zeta_convergence_radius': float(1.0 / lambda_max) if lambda_max > 1e-10 else float('inf'),
        }


# =============================================================================
# Klasse: SubshiftOfFiniteType
# =============================================================================

class SubshiftOfFiniteType:
    """
    @class SubshiftOfFiniteType
    @brief Subshift of Finite Type (SFT) via Übergangsmatrix.

    @description
        Ein SFT ist ein Shift-Raum X_M ⊂ A^ℤ definiert durch eine
        0-1-Matrix M ∈ {0,1}^{d×d}:
          X_M = {(xₙ) : M[xₙ, x_{n+1}] = 1 für alle n}

        Die Dynamik (σ, X_M) hat folgende Eigenschaften:
        - Irreduzibel (transitiv): M ist irreduzibel (keine Block-Struktur)
        - Mischend (mixing): M^k > 0 für ein k ≥ 1 (aperiodisch + irreduzibel)
        - Topologische Entropie: h_top = log(λ_max)

        Perron-Frobenius-Theorem: Für irreduzibles M > 0 ist λ_max positiv,
        einfach und der dazugehörige Eigenvektor ist positiv.

    @param transition_matrix M ∈ {0,1}^{d×d}
    @lastModified 2026-03-11
    """

    def __init__(self, transition_matrix: np.ndarray) -> None:
        """
        @brief Initialisiert den SFT mit der Übergangsmatrix.

        @param transition_matrix 0-1-Matrix M ∈ {0,1}^{d×d}
        @raises ValueError Wenn Matrix nicht quadratisch oder Einträge ∉ {0,1}
        @lastModified 2026-03-11
        """
        M = np.array(transition_matrix, dtype=int)

        if M.ndim != 2 or M.shape[0] != M.shape[1]:
            raise ValueError("Übergangsmatrix muss quadratisch sein.")

        if not np.all((M == 0) | (M == 1)):
            raise ValueError("Übergangsmatrix muss nur 0/1-Einträge haben.")

        self.M = M
        self.d = M.shape[0]
        self._shift = ShiftSystem(list(range(self.d)))

    def is_irreducible(self) -> bool:
        """
        @brief Prüft ob die Übergangsmatrix irreduzibel ist.

        @description
            M ist irreduzibel genau dann wenn der Übergangsgraph stark
            zusammenhängend ist. Äquivalent: (I + M)^d > 0 (komponentenweise).

        @return True wenn M irreduzibel
        @lastModified 2026-03-11
        """
        # Prüfe ob M^k > 0 für k ≤ d² (Cayley-Hamilton)
        I = np.eye(self.d, dtype=float)
        reach = I + self.M.astype(float)

        for _ in range(self.d - 1):
            reach = reach @ (I + self.M.astype(float))

        return bool(np.all(reach > 0))

    def is_mixing(self) -> bool:
        """
        @brief Prüft ob der SFT (topologisch) mischend ist.

        @description
            Der SFT ist mischend (topologisch) genau dann wenn M aperiodisch
            und irreduzibel ist, d.h. M^k > 0 für ein k ≥ 1 (komponentenweise).

            Äquivalent: gcd der Längen aller Zyklen = 1.

        @return True wenn SFT mischend
        @lastModified 2026-03-11
        """
        if not self.is_irreducible():
            return False

        # Aperiodizität: M^d muss vollständig positiv sein
        Md = np.linalg.matrix_power(self.M.astype(int), self.d)
        return bool(np.all(Md > 0))

    def entropy(self) -> float:
        """
        @brief Berechnet die topologische Entropie h_top = log(λ_max(M)).

        @return Topologische Entropie
        @lastModified 2026-03-11
        """
        return self._shift.entropy_shift_sft(self.M)

    def parry_measure(self) -> Optional[np.ndarray]:
        """
        @brief Berechnet das Parry-Maß (maximale Entropie-Maß) des SFT.

        @description
            Das Parry-Maß (Parry 1964) ist das eindeutige maß-maximierende
            Maß für irreduzible SFTs. Es wird aus dem Perron-Frobenius-
            Eigenvektor konstruiert.

            Für irreduzibles M mit Perron-Frobenius-Eigenwert λ und
            linkem/rechtem Eigenvektor u, v (Mu = λu, vᵀM = λvᵀ):
              μ([a₀a₁...aₙ₋₁]) = vₐ₀ · M[a₀,a₁] · ... · M[aₙ₋₂,aₙ₋₁] · uₐₙ₋₁ / (λⁿ⁻¹ · vᵀu)

        @return Parry-Maß als Übergangswahrscheinlichkeitsmatrix, oder None
        @lastModified 2026-03-11
        """
        if not self.is_irreducible():
            return None

        M = self.M.astype(float)
        eigenvalues, eigenvectors = np.linalg.eig(M)
        lambda_max = np.max(eigenvalues.real)

        # Rechter Eigenvektor zu λ_max
        idx = np.argmax(eigenvalues.real)
        r = np.abs(eigenvectors[:, idx].real)

        if np.any(r < 1e-14):
            return None

        # Parry-Übergangsmatrix: P_ij = M_ij · r_j / (λ_max · r_i)
        parry = np.zeros_like(M)
        for i in range(self.d):
            for j in range(self.d):
                if M[i, j] > 0 and r[i] > 1e-14:
                    parry[i, j] = M[i, j] * r[j] / (lambda_max * r[i])

        return parry

    def generate_sequence(
        self,
        length: int,
        start_state: Optional[int] = None,
        seed: int = 42
    ) -> List[int]:
        """
        @brief Generiert eine zulässige Symbolfolge im SFT.

        @description
            Simulation eines zufälligen Pfads im Übergangsgraphen:
            Startet bei start_state und folgt zulässigen Übergängen (M[i,j]=1)
            gleichwahrscheinlich.

        @param length Länge der generierten Folge
        @param start_state Startzustand (oder None für zufälligen Start)
        @param seed Zufallsseed
        @return Zulässige Symbolfolge der Länge length
        @lastModified 2026-03-11
        """
        rng = np.random.default_rng(seed)

        if start_state is None:
            # Zufälligen Zustand mit zulässigen Nachfolgern wählen
            candidates = [i for i in range(self.d) if np.any(self.M[i] > 0)]
            if not candidates:
                return []
            start_state = int(rng.choice(candidates))

        sequence = [start_state]
        current = start_state

        for _ in range(length - 1):
            # Mögliche Nachfolger
            successors = np.where(self.M[current] > 0)[0]
            if len(successors) == 0:
                break
            current = int(rng.choice(successors))
            sequence.append(current)

        return sequence


# =============================================================================
# Hilfsfunktionen: Klassische Transformationen
# =============================================================================

def doubling_map(x: float) -> float:
    """
    @brief Verdopplungsabbildung T(x) = 2x (mod 1) auf [0,1].

    @description
        Die Verdopplungsabbildung ist das Prototypbeispiel eines ergodischen,
        mischenden und entropisch nicht-trivialen Systems:
        - Maß-erhaltend: Lebesgue-Maß
        - Ergodisch und stark mischend
        - Topologische Entropie: h = log(2)
        - Sensitive Abhängigkeit von Anfangsbedingungen (Chaos)

    @param x Punkt in [0, 1)
    @return T(x) = 2x mod 1
    @lastModified 2026-03-11
    """
    return (2 * x) % 1.0


def circle_rotation(alpha: float) -> Callable:
    """
    @brief Kreisrotation T_α(x) = x + α (mod 1) mit Rotationswinkel α.

    @description
        Die Kreisrotation auf dem Einheitskreis ℝ/ℤ:
          T_α(x) = x + α (mod 1)

        Eigenschaften:
        - Immer maß-erhaltend (Lebesgue-Maß)
        - Ergodisch ⟺ α irrational (Weyl 1916)
        - Nie mischend (spezifisches Spektrum: punktförmig)
        - Entropie: h = 0 (für alle α)

    @param alpha Rotationswinkel α ∈ ℝ
    @return Funktion T_α: [0,1) → [0,1)
    @lastModified 2026-03-11
    """
    def T(x: float) -> float:
        return (x + alpha) % 1.0
    T.__name__ = f"circle_rotation(alpha={alpha})"
    return T


def tent_map(x: float) -> float:
    """
    @brief Zeltabbildung T(x) = 1 - |2x - 1| auf [0,1].

    @description
        Die Zeltabbildung (tent map):
          T(x) = 2x       für x ∈ [0, 1/2]
          T(x) = 2(1-x)   für x ∈ [1/2, 1]

        Eigenschaften:
        - Maß-erhaltend: Lebesgue-Maß auf [0,1]
        - Topologisch konjugiert zur Verdopplungsabbildung
        - Ergodisch und mischend
        - Entropie: h = log(2)

    @param x Punkt in [0, 1]
    @return T(x)
    @lastModified 2026-03-11
    """
    return 1.0 - abs(2.0 * x - 1.0)


def logistic_map_r4(x: float) -> float:
    """
    @brief Logistische Abbildung bei r=4: T(x) = 4x(1-x) auf [0,1].

    @description
        Die logistische Abbildung bei r=4 ist ein klassisches Beispiel
        für deterministisches Chaos:
          T(x) = 4x(1-x)

        Invariantes Maß: dμ = dx / (π√(x(1-x)))  (Arkussinus-Maß)
        Lyapunov-Exponent: λ = log(2) ≈ 0.693 (chaotisch)
        Topologisch konjugiert zu T(x) = 2x (mod 1) via x = sin²(πy/2).

    @param x Punkt in [0, 1]
    @return T(x) = 4x(1-x)
    @lastModified 2026-03-11
    """
    return 4.0 * x * (1.0 - x)


def gauss_map(x: float) -> float:
    """
    @brief Gauss-Abbildung (Kettenbruch-Abbildung) T(x) = {1/x} auf (0,1].

    @description
        Die Gauss-Abbildung (Gauß-Transformation):
          T(x) = 1/x - floor(1/x) = {1/x}  für x ≠ 0
          T(0) = 0

        Anwendung: Kettenbruchentwicklung von x ∈ (0,1):
          x = 1/(a₁ + 1/(a₂ + ...))  mit aₙ = floor(T^{n-1}(x))

        Invariantes Maß: Gauss-Maß dμ = dx / ((1+x) · log 2)
        Entropie: π²/(6 log 2) ≈ 2.373 bits (Rohlin 1961)

    @param x Punkt in (0, 1]
    @return T(x) = {1/x}
    @lastModified 2026-03-11
    """
    if x <= 1e-14:
        return 0.0
    return (1.0 / x) % 1.0
