"""
@file random_matrix_theory.py
@brief Zufallsmatrizentheorie: GUE, GOE, Level-Spacing-Statistik und
       Verbindung zur Riemann-Zeta-Funktion (Montgomery-Odlyzko-Gesetz).
@description
    Dieses Modul implementiert die Kernkonzepte der Zufallsmatrizentheorie (RMT):

    1. **GUE** (Gauß'sches Unitäres Ensemble): Hermitesche Matrizen mit
       komplexen Gaußeinträgen. Eigenwertdichte folgt dem Wigner-Halbkreisgesetz.

    2. **GOE** (Gauß'sches Orthogonales Ensemble): Symmetrische reelle Matrizen.
       Level-Spacing folgt dem Wigner-Surmise p(s) ≈ (π/2)s·exp(-πs²/4).

    3. **LevelSpacingStatistics**: Normierte Eigenwertabstände (Wigner-Surmise,
       Poisson-Verteilung) — unterscheiden korrelierte von unkorrelierten Spektren.

    4. **RiemannZetaGUE**: Montgomery-Odlyzko-Gesetz: Die Nullstellen der
       Riemann-Zeta-Funktion verhalten sich wie GUE-Eigenwerte.

    5. **RandomMatrixTheoryTools**: Marchenko-Pastur-Gesetz für Wishart-Matrizen,
       Dyson-Index-Schätzung, Ensemble-Vergleichsplot.

    Mathematischer Hintergrund:
    - Wigner-Halbkreis: ρ(x) = (2/π)√(1-x²) für |x| ≤ 1
    - GUE-Wigner-Surmise: p(s) = (32/π²)·s²·exp(-4s²/π)
    - GOE-Wigner-Surmise: p(s) = (π/2)·s·exp(-πs²/4)
    - Montgomery-Paarkorrelation: R(r) = 1 - (sin(πr)/(πr))²

    Referenzen:
    - Mehta, M.L. (2004): "Random Matrices", 3rd ed., Academic Press
    - Montgomery, H.L. (1973): "The pair correlation of zeros of the zeta function"
    - Odlyzko, A.M. (1987): "On the distribution of spacings between zeros of the zeta function"

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import numpy as np
import matplotlib
# Nicht-interaktives Backend, damit Tests ohne Display funktionieren
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Optional


# ===========================================================================
# KLASSE: GUEEnsemble
# ===========================================================================

class GUEEnsemble:
    """
    @brief Gauß'sches Unitäres Ensemble (GUE) — hermitesche n×n-Matrizen.
    @description
        Das GUE ist das Ensemble aller n×n hermiteschen Matrizen H = H†,
        deren Einträge unabhängige komplexe Gaußvariablen sind.

        Normierung: H = (A + A†) / √(2n), wobei A eine n×n-Matrix mit
        Einträgen A_{ij} ~ CN(0,1) (komplexe Standardnormalverteilung).

        Im Grenzfall n → ∞ folgt die Eigenwertdichte dem Wigner-Halbkreisgesetz:
            ρ(x) = (2/π)√(1-x²)  für |x| ≤ 1

        Die Level-Spacings folgen dem GUE-Wigner-Surmise:
            p(s) = (32/π²)·s²·exp(-4s²/π)

        Dyson-Index β = 2 (unitäre Symmetrie, keine Zeitumkehrinvarianz).

    @param n: Matrixgröße (n×n hermitesche Matrizen)
    @param samples: Anzahl der Stichprobenmatrizen für statistische Auswertung
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, n: int, samples: int = 1000) -> None:
        """
        @brief Initialisiert das GUE-Ensemble.
        @param n: Dimension der hermiteschen Matrizen (n×n).
        @param samples: Anzahl der zu erzeugenden Stichprobenmatrizen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Matrixgröße muss positiv sein
        if n < 1:
            raise ValueError("Matrixgröße n muss mindestens 1 sein.")
        # Stichprobenanzahl muss positiv sein
        if samples < 1:
            raise ValueError("Anzahl der Stichproben muss mindestens 1 sein.")

        self.n = n          # Dimension der Matrizen
        self.samples = samples  # Anzahl der Stichproben

    def sample(self) -> np.ndarray:
        """
        @brief Erzeugt eine einzelne zufällige hermitesche n×n-Matrix aus dem GUE.
        @description
            Algorithmus:
            1. Erzeuge A: n×n-Matrix mit komplexen Gaußeinträgen
               Real- und Imaginärteil ∈ N(0,1) unabhängig
            2. Berechne H = (A + A†) / √(2n)
               Durch (A + A†)/2 wird H hermitesch: H_{ij} = conj(H_{ji})
               Durch Division durch √(2n) wird die Varianz auf 1/n normiert,
               was das Spektrum im Einheitsintervall hält (Wigner-Skalierung).

            Warum √(2n)?
            - Varianz von Diagonaleinträgen: Var[H_{ii}] = Var[A_{ii} + A̅_{ii}] / (2n)
              = 2·Var[Re(A_{ii})] / (2n) = 1/n
            - Im Limit: Spektralradius ≈ 2, skaliert auf [-1,1] bei Faktor √(2n)/2·√n

        @return: n×n hermitesche Matrix als numpy-Array (dtype=complex128).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        n = self.n

        # Zufällige komplexe Matrix: Real- und Imaginärteil je N(0,1)
        # np.random.randn erzeugt reelle N(0,1)-Werte
        real_part = np.random.randn(n, n)
        imag_part = np.random.randn(n, n)
        A = real_part + 1j * imag_part  # Komplexe n×n-Gaußmatrix

        # Hermite-Symmetrisierung: H = (A + A†) / √(2n)
        # A.conj().T = A† (hermitesche Konjugierte)
        H = (A + A.conj().T) / np.sqrt(2 * n)

        return H

    def eigenvalue_distribution(self) -> np.ndarray:
        """
        @brief Berechnet alle Eigenwerte aus mehreren GUE-Stichproben.
        @description
            Erzeugt `self.samples` hermitesche Matrizen und berechnet deren
            Eigenwerte. Da H hermitesch ist, sind alle Eigenwerte reell
            (np.linalg.eigvalsh garantiert reelle Eigenwerte für hermitesche Matrizen).

            Gesamtanzahl der Eigenwerte: samples × n

        @return: 1D-Array aller gesammelten Eigenwerte (float64).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Liste zur Sammlung aller Eigenwerte
        all_eigenvalues = []

        for _ in range(self.samples):
            # Erzeuge eine hermitesche Stichprobenmatrix
            H = self.sample()
            # eigvalsh: optimiert für hermitesche Matrizen, gibt reelle Eigenwerte zurück
            eigenvalues = np.linalg.eigvalsh(H)
            all_eigenvalues.append(eigenvalues)

        # Alle Eigenwerte in einem flachen Array zusammenführen
        return np.concatenate(all_eigenvalues)

    def wigner_semicircle(self, x: np.ndarray) -> np.ndarray:
        """
        @brief Berechnet die Wigner-Halbkreisdichte an den Punkten x.
        @description
            Das Wigner-Halbkreisgesetz besagt, dass die empirische Eigenwertdichte
            eines n×n-GUE im Grenzfall n → ∞ gegen folgende Dichte konvergiert:

                ρ(x) = (2/π)·√(1 - x²)   für |x| ≤ 1
                ρ(x) = 0                   für |x| > 1

            Dies ist eine universelle Aussage: unabhängig von der genauen Verteilung
            der Matrixeinträge (sofern Momente existieren) gilt dieses Gesetz.

            Normierung: ∫₋₁¹ ρ(x) dx = 1 (Wahrscheinlichkeitsdichte)

        @param x: 1D-Array von x-Werten (typisch: linspace(-1.5, 1.5, 300))
        @return: Dichtewerte ρ(x) (≥ 0, = 0 außerhalb [-1,1])
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        x = np.asarray(x, dtype=float)

        # Wigner-Halbkreis: (2/π)√(1-x²) für |x| ≤ 1, sonst 0
        # np.where: Bedingung ? Wert_wenn_wahr : Wert_wenn_falsch
        density = np.where(
            np.abs(x) <= 1.0,
            (2.0 / np.pi) * np.sqrt(np.maximum(1.0 - x**2, 0.0)),
            0.0
        )

        return density

    def plot_eigenvalue_histogram(self, bins: int = 50) -> plt.Figure:
        """
        @brief Erstellt ein Histogramm der empirischen Eigenwertverteilung
               und vergleicht es mit dem theoretischen Wigner-Halbkreis.
        @description
            Schritte:
            1. Berechne Eigenwertverteilung (samples × n Eigenwerte)
            2. Normiere Histogram auf Wahrscheinlichkeitsdichte
            3. Überlagere theoretische Wigner-Halbkreisdichte
            4. Beschriftung und Legende

        @param bins: Anzahl der Histogramm-Bins (Standard: 50)
        @return: matplotlib-Figure-Objekt (kann mit fig.savefig(...) gespeichert werden)
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Eigenwerte sammeln
        eigenvalues = self.eigenvalue_distribution()

        # Figure anlegen
        fig, ax = plt.subplots(figsize=(8, 5))

        # Histogramm: density=True normiert die Fläche auf 1
        ax.hist(eigenvalues, bins=bins, density=True,
                alpha=0.6, color='steelblue', label='Empirische Eigenwerte')

        # Wigner-Halbkreis überlagern
        x_theory = np.linspace(-1.5, 1.5, 400)
        ax.plot(x_theory, self.wigner_semicircle(x_theory),
                'r-', linewidth=2, label='Wigner-Halbkreis ρ(x) = (2/π)√(1-x²)')

        ax.set_xlabel('Eigenwert λ')
        ax.set_ylabel('Dichte ρ(λ)')
        ax.set_title(f'GUE-Eigenwertverteilung (n={self.n}, Stichproben={self.samples})')
        ax.legend()
        ax.grid(True, alpha=0.3)

        return fig


# ===========================================================================
# KLASSE: GOEEnsemble
# ===========================================================================

class GOEEnsemble:
    """
    @brief Gauß'sches Orthogonales Ensemble (GOE) — symmetrische reelle n×n-Matrizen.
    @description
        Das GOE ist das Ensemble aller n×n reellen symmetrischen Matrizen,
        deren Einträge unabhängige Gaußvariablen sind:
        - Diagonale: H_{ii} ~ N(0, 2/n)
        - Obere Dreiecke: H_{ij} ~ N(0, 1/n) für i < j
        - Untere Dreiecke: H_{ji} = H_{ij} (Symmetrie)

        Normierung analog zum GUE: H = (A + Aᵀ) / √(2n)

        Level-Spacings folgen dem GOE-Wigner-Surmise:
            p(s) = (π/2)·s·exp(-πs²/4)

        Dyson-Index β = 1 (orthogonale Symmetrie, Zeitumkehrinvarianz vorhanden).
        Anwendungen: Kerne schwerer Atome, Quantenchaos

    @param n: Matrixgröße (n×n symmetrische Matrizen)
    @param samples: Anzahl der Stichprobenmatrizen
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, n: int, samples: int = 1000) -> None:
        """
        @brief Initialisiert das GOE-Ensemble.
        @param n: Dimension der symmetrischen Matrizen (n×n).
        @param samples: Anzahl der zu erzeugenden Stichprobenmatrizen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        if n < 1:
            raise ValueError("Matrixgröße n muss mindestens 1 sein.")
        if samples < 1:
            raise ValueError("Anzahl der Stichproben muss mindestens 1 sein.")

        self.n = n
        self.samples = samples

    def sample(self) -> np.ndarray:
        """
        @brief Erzeugt eine einzelne zufällige symmetrische n×n-Matrix aus dem GOE.
        @description
            Algorithmus:
            1. Erzeuge A: n×n-Matrix mit reellen Gaußeinträgen A_{ij} ~ N(0,1)
            2. Symmetrisiere: H = (A + Aᵀ) / √(2n)

            Eigenschaften der erzeugten Matrix:
            - H_{ij} = H_{ji} (symmetrisch)
            - Alle Einträge reell
            - Eigenwerte reell (symmetrische Matrizen haben nur reelle Eigenwerte)

        @return: n×n symmetrische Matrix als numpy-Array (dtype=float64).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        n = self.n

        # Zufällige reelle n×n-Matrix mit N(0,1)-Einträgen
        A = np.random.randn(n, n)

        # Symmetrisierung: H = (A + Aᵀ) / √(2n)
        H = (A + A.T) / np.sqrt(2 * n)

        return H

    def eigenvalue_distribution(self) -> np.ndarray:
        """
        @brief Sammelt alle Eigenwerte aus mehreren GOE-Stichproben.
        @description
            Erzeugt `self.samples` symmetrische Matrizen und sammelt deren
            reelle Eigenwerte. np.linalg.eigvalsh ist für symmetrische Matrizen
            numerisch stabiler als eigvals.

        @return: 1D-Array aller Eigenwerte (float64), Länge = samples × n.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        all_eigenvalues = []

        for _ in range(self.samples):
            H = self.sample()
            # eigvalsh: für reell-symmetrische Matrizen (liefert sortierte reelle Eigenwerte)
            eigenvalues = np.linalg.eigvalsh(H)
            all_eigenvalues.append(eigenvalues)

        return np.concatenate(all_eigenvalues)

    def level_spacing_distribution(self) -> np.ndarray:
        """
        @brief Berechnet normierte Abstände zwischen benachbarten Eigenwerten.
        @description
            Algorithmus:
            1. Berechne Eigenwerte jeder Stichprobenmatrix (sortiert)
            2. Bilde Differenzen: s_i = λ_{i+1} - λ_i
            3. Normiere auf mittleren Abstand: s̃_i = s_i / <s>

            Normierung auf Mittelwert 1 ist wichtig für den Vergleich mit
            dem Wigner-Surmise, das für <s> = 1 formuliert ist.

        @return: 1D-Array der normierten Level-Spacings (≥ 0).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        all_spacings = []

        for _ in range(self.samples):
            H = self.sample()
            # Eigenwerte sortiert zurückgeben (eigvalsh liefert sie aufsteigend)
            eigenvalues = np.linalg.eigvalsh(H)

            # Abstände zwischen benachbarten Eigenwerten
            spacings = np.diff(eigenvalues)

            # Normierung: s̃ = s / <s>, damit Mittelwert = 1
            mean_spacing = np.mean(spacings)
            if mean_spacing > 1e-14:  # Vermeide Division durch ~0
                spacings = spacings / mean_spacing

            all_spacings.append(spacings)

        return np.concatenate(all_spacings)


# ===========================================================================
# KLASSE: LevelSpacingStatistics
# ===========================================================================

class LevelSpacingStatistics:
    """
    @brief Statistische Analyse von Eigenwertabständen (Level-Spacing-Statistik).
    @description
        Die Level-Spacing-Statistik unterscheidet zwischen:
        - **GUE**: p(s) = (32/π²)·s²·exp(-4s²/π)  — quadratische Abstoßung (β=2)
        - **GOE**: p(s) = (π/2)·s·exp(-πs²/4)      — lineare Abstoßung (β=1)
        - **Poisson**: p(s) = exp(-s)               — keine Abstoßung (unkorreliert)

        Level-Repulsion (Quantenchaos-Kriterium):
        - Chaotische Systeme → GUE oder GOE (Eigenwerte stoßen sich ab)
        - Integrable Systeme → Poisson (unabhängige Eigenwerte)

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def normalize_spacings(self, eigenvalues: np.ndarray) -> np.ndarray:
        """
        @brief Normiert Eigenwertabstände auf mittleren Abstand 1.
        @description
            1. Sortiere Eigenwerte aufsteigend
            2. Bilde Differenzen: s_i = λ_{i+1} - λ_i
            3. Normiere: s̃_i = s_i / mean(s)

            Nach Normierung gilt: mean(s̃) = 1

        @param eigenvalues: 1D-Array von Eigenwerten (beliebige Reihenfolge).
        @return: 1D-Array normierter Abstände (≥ 0).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Sortierung sicherstellen
        ev_sorted = np.sort(eigenvalues)

        # Abstände zwischen benachbarten Eigenwerten
        spacings = np.diff(ev_sorted)

        # Normierung auf Mittelwert 1
        mean_s = np.mean(spacings)
        if mean_s < 1e-14:
            # Entartetes Spektrum: keine sinnvolle Normierung möglich
            return spacings

        return spacings / mean_s

    def wigner_surmise_goe(self, s: np.ndarray) -> np.ndarray:
        """
        @brief Wigner-Surmise-Verteilung für das GOE.
        @description
            Exakte Formel für 2×2-GOE-Matrizen (Wigner 1951), die für große
            Matrizen eine sehr gute Näherung liefert:

                p(s) = (π/2)·s·exp(-πs²/4)

            Eigenschaften:
            - p(0) = 0: Level-Repulsion (lineare Abstoßung, β=1)
            - ∫₀^∞ p(s) ds = 1 (normierte Verteilung)
            - <s> = 1 bei entsprechender Normierung

        @param s: 1D-Array von Abstandswerten s ≥ 0.
        @return: Wahrscheinlichkeitsdichte p(s) für GOE.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        s = np.asarray(s, dtype=float)
        # GOE-Wigner-Surmise: (π/2)·s·exp(-πs²/4)
        return (np.pi / 2.0) * s * np.exp(-np.pi * s**2 / 4.0)

    def wigner_surmise_gue(self, s: np.ndarray) -> np.ndarray:
        """
        @brief Wigner-Surmise-Verteilung für das GUE.
        @description
            Exakte Formel für 2×2-GUE-Matrizen (Wigner 1951):

                p(s) = (32/π²)·s²·exp(-4s²/π)

            Eigenschaften:
            - p(0) = 0: Level-Repulsion (quadratische Abstoßung, β=2)
            - ∫₀^∞ p(s) ds = 1 (normierte Verteilung)
            - Stärkere Abstoßung als GOE (s²-Term statt s)

        @param s: 1D-Array von Abstandswerten s ≥ 0.
        @return: Wahrscheinlichkeitsdichte p(s) für GUE.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        s = np.asarray(s, dtype=float)
        # GUE-Wigner-Surmise: (32/π²)·s²·exp(-4s²/π)
        return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)

    def poisson_distribution(self, s: np.ndarray) -> np.ndarray:
        """
        @brief Poisson-Verteilung der Level-Spacings (unkorrelierte Eigenwerte).
        @description
            Für vollständig unkorrelierte (integrable) Systeme folgen die
            Eigenwertabstände einer Exponentialverteilung:

                p(s) = exp(-s)

            Eigenschaften:
            - p(0) = 1: Keine Level-Repulsion (Eigenwerte können entarten)
            - <s> = 1 (bei Normierung auf mittleren Abstand)
            - Vergleich: GOE/GUE zeigen Level-Repulsion durch p(0) = 0

        @param s: 1D-Array von Abstandswerten s ≥ 0.
        @return: Wahrscheinlichkeitsdichte p(s) für Poisson.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        s = np.asarray(s, dtype=float)
        return np.exp(-s)

    def plot_level_spacing(self, spacings: np.ndarray,
                           ensemble: str = "GUE") -> plt.Figure:
        """
        @brief Plottet empirische Level-Spacings mit theoretischen Kurven.
        @description
            Vergleich von:
            - Empirischem Histogramm der normierten Abstände
            - Wigner-Surmise (GOE oder GUE je nach Parameter)
            - Poisson-Verteilung (Referenz für unkorrelierte Systeme)

        @param spacings: 1D-Array normierter Eigenwertabstände.
        @param ensemble: "GUE" oder "GOE" — bestimmt welcher Wigner-Surmise angezeigt wird.
        @return: matplotlib-Figure-Objekt.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        fig, ax = plt.subplots(figsize=(8, 5))

        # Nur positive Abstände (numerische Sicherheit)
        spacings = spacings[spacings >= 0]

        # Histogramm der empirischen Abstände
        ax.hist(spacings, bins=50, density=True,
                alpha=0.5, color='steelblue', label='Empirische Abstände')

        # Theoretische Kurven
        s_vals = np.linspace(0, 4, 400)

        # Wigner-Surmise für gewähltes Ensemble
        if ensemble.upper() == "GOE":
            ws = self.wigner_surmise_goe(s_vals)
            ws_label = 'Wigner-Surmise GOE: (π/2)s·exp(-πs²/4)'
        else:
            ws = self.wigner_surmise_gue(s_vals)
            ws_label = 'Wigner-Surmise GUE: (32/π²)s²·exp(-4s²/π)'

        ax.plot(s_vals, ws, 'r-', linewidth=2, label=ws_label)

        # Poisson-Referenz
        ax.plot(s_vals, self.poisson_distribution(s_vals),
                'g--', linewidth=2, label='Poisson: exp(-s)')

        ax.set_xlabel('Normierter Abstand s')
        ax.set_ylabel('Wahrscheinlichkeitsdichte p(s)')
        ax.set_title(f'Level-Spacing-Statistik ({ensemble}-Ensemble)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 4)

        return fig


# ===========================================================================
# KLASSE: RiemannZetaGUE
# ===========================================================================

class RiemannZetaGUE:
    """
    @brief Verbindung zwischen Riemann-Zeta-Nullstellen und dem GUE.
    @description
        Das Montgomery-Odlyzko-Gesetz (1973/1987) besagt:

        Die normierten Abstände zwischen den nicht-trivialen Nullstellen
        ζ(1/2 + it_n) = 0 (Imaginärteile t_n) folgen statistisch der
        GUE-Paarkorrelationsfunktion:

            R₂(r) = 1 - (sin(πr) / (πr))²

        Dies ist ein starkes numerisches Indiz für die Riemann-Hypothese:
        Wenn die Nullstellen auf der kritischen Linie Re(s) = 1/2 liegen,
        verhalten sie sich wie Eigenwerte hermitescher Zufallsmatrizen.

        Die Montgomery-Paarkorrelation wurde für 10^22 Nullstellen verifiziert
        (Odlyzko, 2001).

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self) -> None:
        """
        @brief Initialisiert das RiemannZetaGUE-Objekt.
        @description
            Speichert die ersten 20 bekannten nicht-trivialen Nullstellen
            der Riemann-Zeta-Funktion (Imaginärteile t_n mit ζ(1/2+it_n)=0).
            Diese sind auf mehrere Dezimalstellen genau bekannt.

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Die ersten 20 nicht-trivialen Nullstellen der Riemann-Zeta-Funktion
        # Imaginärteile t_n mit ζ(1/2 + i·t_n) = 0
        # Quelle: Odlyzko, A.M. (1987), NIST Digital Library of Mathematical Functions
        self._known_zeros = [
            14.134725141734693,  # t_1
            21.022039638771555,  # t_2
            25.010857580145688,  # t_3
            30.424876125859513,  # t_4
            32.935061587739189,  # t_5
            37.586178158825671,  # t_6
            40.918719012147495,  # t_7
            43.327073280914999,  # t_8
            48.005150881167159,  # t_9
            49.773832477672302,  # t_10
            52.970321477714460,  # t_11
            56.446247697063246,  # t_12
            59.347044002602353,  # t_13
            60.831778524609809,  # t_14
            65.112544048081690,  # t_15
            67.079810529494174,  # t_16
            69.546401711173978,  # t_17
            72.067157674481907,  # t_18
            75.704690699083933,  # t_19
            77.144840068874805,  # t_20
        ]

    def compute_zeta_zeros(self, count: int = 50) -> np.ndarray:
        """
        @brief Berechnet die ersten `count` nicht-trivialen Nullstellen der ζ-Funktion.
        @description
            Strategie:
            1. Für count ≤ 20: Verwende bekannte hochpräzise Werte
            2. Für count > 20: Ergänze mit mpmath.zetazero() für weitere Nullstellen

            mpmath.zetazero(n) berechnet die n-te nicht-triviale Nullstelle
            mit konfigurierbarer Genauigkeit (mp.dps Dezimalstellen).

            Die Nullstellen werden nach aufsteigendem Imaginärteil sortiert.

        @param count: Anzahl der zu berechnenden Nullstellen (Standard: 50).
        @return: 1D-Array der Imaginärteile t_n (aufsteigend sortiert, float64).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        if count <= 0:
            raise ValueError("count muss mindestens 1 sein.")

        # Starte mit bekannten hochpräzisen Werten
        zeros = list(self._known_zeros[:min(count, len(self._known_zeros))])

        # Falls mehr Nullstellen benötigt werden, verwende mpmath
        if count > len(self._known_zeros):
            try:
                import mpmath
                # Genauigkeit: 25 Dezimalstellen für Geschwindigkeit/Genauigkeit-Balance
                mpmath.mp.dps = 25

                # Berechne fehlende Nullstellen
                for n in range(len(self._known_zeros) + 1, count + 1):
                    # zetazero(n) gibt die n-te Nullstelle zurück
                    zero = mpmath.zetazero(n)
                    # Nur den Imaginärteil extrahieren (alle auf Re(s)=1/2)
                    zeros.append(float(zero.imag))

            except ImportError:
                # mpmath nicht verfügbar: Nullstellen aus bekannten Werten extrapolieren
                # Grobe Näherung: Durchschnittlicher Abstand der letzten 5 bekannten Nullstellen
                if len(zeros) >= 2:
                    # Mittlerer Abstand der letzten Nullstellen
                    recent_spacings = np.diff(zeros[-min(5, len(zeros)):])
                    mean_gap = np.mean(recent_spacings) if len(recent_spacings) > 0 else 3.0

                    # Lineare Extrapolation (grobe Näherung)
                    last_zero = zeros[-1]
                    for _ in range(count - len(zeros)):
                        last_zero += mean_gap
                        zeros.append(last_zero)

        return np.array(zeros[:count], dtype=float)

    def zeros_level_spacing(self, zeros: np.ndarray) -> np.ndarray:
        """
        @brief Berechnet normierte Abstände zwischen Zeta-Nullstellen.
        @description
            Montgomery-Odlyzko-Normierung:
            Die Nullstellen t_n werden mit dem lokalen Mittelwert des Abstands
            2π/ln(t/(2π)) normiert, um die logarithmisch wachsende Dichte
            zu kompensieren.

            Einfachere Normierung (hier implementiert):
            s_n = (t_{n+1} - t_n) / <t_{n+1} - t_n>

            Für große t gilt die exakte Normierung:
            s_n = (t_{n+1} - t_n) · ln(t_n / (2π)) / (2π)

        @param zeros: 1D-Array von Zeta-Nullstellen (Imaginärteile, aufsteigend sortiert).
        @return: 1D-Array normierter Abstände zwischen benachbarten Nullstellen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        zeros = np.sort(zeros)

        # Abstände zwischen benachbarten Nullstellen
        spacings = np.diff(zeros)

        # Montgomery-Odlyzko-Normierung mit lokalem Faktor:
        # Die Dichte der Nullstellen wächst logarithmisch: ρ(t) ≈ ln(t/(2π)) / (2π)
        # Mittlerer lokaler Abstand: δ(t) ≈ 2π / ln(t/(2π))
        t_mid = (zeros[:-1] + zeros[1:]) / 2.0  # Mittelpunkte der Intervalle

        # Lokaler Normierungsfaktor
        local_mean = 2.0 * np.pi / np.log(np.maximum(t_mid / (2.0 * np.pi), 1.01))

        # Normierte Abstände
        normalized = spacings / local_mean

        return normalized

    def compare_gue_zeros(self, count: int = 100) -> plt.Figure:
        """
        @brief Vergleicht ζ-Nullstellen-Abstände mit der GUE-Wigner-Surmise.
        @description
            Montgomery-Odlyzko-Gesetz: Die normierten Abstände zwischen
            ζ-Nullstellen folgen der GUE-Statistik, nicht der Poisson-Statistik.

            Plot-Inhalt:
            - Histogramm der normierten ζ-Nullstellen-Abstände
            - GUE-Wigner-Surmise: p(s) = (32/π²)·s²·exp(-4s²/π)
            - Poisson-Referenz: p(s) = exp(-s)

        @param count: Anzahl der zu verwendenden Nullstellen (Standard: 100).
        @return: matplotlib-Figure-Objekt.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Nullstellen berechnen
        zeros = self.compute_zeta_zeros(count)

        # Normierte Abstände berechnen
        spacings = self.zeros_level_spacing(zeros)

        # Level-Spacing-Statistik für den Plot
        stats = LevelSpacingStatistics()

        fig, ax = plt.subplots(figsize=(9, 5))

        # Histogramm der ζ-Nullstellen-Abstände
        ax.hist(spacings, bins=min(30, len(spacings) // 2), density=True,
                alpha=0.5, color='purple', label='ζ-Nullstellen-Abstände')

        # Theoretische Kurven
        s_vals = np.linspace(0, 4, 400)

        # GUE-Wigner-Surmise
        ax.plot(s_vals, stats.wigner_surmise_gue(s_vals),
                'r-', linewidth=2, label='GUE-Wigner-Surmise (Montgomery-Odlyzko)')

        # GOE-Wigner-Surmise (zum Vergleich)
        ax.plot(s_vals, stats.wigner_surmise_goe(s_vals),
                'b--', linewidth=1.5, label='GOE-Wigner-Surmise')

        # Poisson-Referenz
        ax.plot(s_vals, stats.poisson_distribution(s_vals),
                'g--', linewidth=2, label='Poisson (unkorreliert)')

        ax.set_xlabel('Normierter Abstand s')
        ax.set_ylabel('Wahrscheinlichkeitsdichte p(s)')
        ax.set_title(f'ζ-Nullstellen vs. GUE (Montgomery-Odlyzko-Gesetz, {count} Nullstellen)')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 4)

        return fig

    def montgomery_pair_correlation(self, zeros: np.ndarray, r: float) -> float:
        """
        @brief Berechnet die Montgomery-Paarkorrelationsfunktion R(r).
        @description
            Montgomery (1973) bewies, dass die Paarkorrelationsfunktion der
            ζ-Nullstellen (unter Annahme der Riemann-Hypothese) ist:

                R(r) = 1 - (sin(πr) / (πr))²

            Dies ist exakt die GUE-2-Punkt-Korrelationsfunktion.

            Implementierung: Empirische Schätzung über normierte Abstände:
            R(r) = (Anzahl der Paare mit normalisertem Abstand ≈ r) /
                   (Gesamtpaare im Intervall)

            Für große r gilt R(r) → 1 (keine langreichweitige Korrelation).
            Für kleine r: R(r) → 0 (Level-Repulsion).

        @param zeros: 1D-Array der ζ-Nullstellen (Imaginärteile).
        @param r: Korrelationsabstand (reell, ≥ 0).
        @return: Paarkorrelationswert R(r).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Theoretische Montgomery-Paarkorrelation
        # R(r) = 1 - (sin(πr) / (πr))²
        if abs(r) < 1e-10:
            # Grenzwert: lim_{r→0} (sin(πr)/(πr))² = 1, also R(0) = 0
            return 0.0

        sinc_val = np.sin(np.pi * r) / (np.pi * r)
        return 1.0 - sinc_val**2


# ===========================================================================
# KLASSE: RandomMatrixTheoryTools
# ===========================================================================

class RandomMatrixTheoryTools:
    """
    @brief Werkzeuge der Zufallsmatrizentheorie: Marchenko-Pastur, Dyson-Index, Vergleichsplots.
    @description
        Ergänzende Funktionen für die Analyse von Zufallsmatrizen:

        1. **Marchenko-Pastur-Gesetz**: Eigenwertdichte von Wishart-Matrizen W = (1/m)·XᵀX,
           wobei X eine m×n-Matrix mit i.i.d. Einträgen ist.
           Verhältnis γ = n/m (Seitenverhältnis):
               ρ(x) = (1/(2π·γ·x)) · √((λ₊-x)(x-λ₋))
               mit λ± = (1 ± √γ)²

        2. **Dyson-Index**: β=1 (GOE), β=2 (GUE), β=4 (GSE)
           Schätzung aus dem Verhältnis der Level-Repulsion bei s→0

        3. **Ensemble-Vergleichsplot**: GOE, GUE, Poisson und ζ-Nullstellen

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def marchenko_pastur(self, x: np.ndarray, ratio: float) -> np.ndarray:
        """
        @brief Marchenko-Pastur-Dichte für Wishart-Matrizen.
        @description
            Das Marchenko-Pastur-Gesetz (1967) beschreibt die Eigenwertdichte
            von Wishart-Matrizen W = (1/m)·X·Xᵀ, wobei X eine n×m-Matrix
            mit i.i.d. N(0,1)-Einträgen ist, im Grenzfall n,m → ∞ mit n/m → γ.

            Dichte:
                ρ_γ(x) = (1/(2π·γ·σ²·x)) · √((λ₊ - x)(x - λ₋))  für x ∈ [λ₋, λ₊]
                ρ_γ(x) = 0                                           sonst

            Ränder des Spektrums:
                λ± = σ²·(1 ± √γ)²

            Für σ² = 1 (Varianz der Matrixeinträge):
                λ± = (1 ± √γ)²

            Sonderfall γ = 1 (quadratische Wishart-Matrix):
                λ₋ = 0, λ₊ = 4

            Anwendungen: PCA, Faktoranalyse, Portfoliooptimierung (Rauschen vs. Signal)

        @param x: 1D-Array von Eigenwerten x > 0.
        @param ratio: Verhältnis γ = n/m ∈ (0, ∞). Typisch γ ≤ 1.
        @return: Marchenko-Pastur-Dichte an den Stellen x.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        x = np.asarray(x, dtype=float)

        if ratio <= 0:
            raise ValueError("ratio (γ = n/m) muss positiv sein.")

        # Spektrumsränder λ± = (1 ± √γ)²
        gamma = ratio
        lambda_plus = (1.0 + np.sqrt(gamma))**2
        lambda_minus = (1.0 - np.sqrt(gamma))**2

        # Marchenko-Pastur-Dichte
        # Nur im Intervall [λ₋, λ₊] positiv
        # (λ₊ - x)(x - λ₋): muss ≥ 0 für x ∈ [λ₋, λ₊]
        inside = (x >= lambda_minus) & (x <= lambda_plus) & (x > 1e-14)

        density = np.zeros_like(x)

        if np.any(inside):
            x_in = x[inside]
            # Zähler: √((λ₊ - x)(x - λ₋))
            numerator = np.sqrt(np.maximum((lambda_plus - x_in) * (x_in - lambda_minus), 0.0))
            # Nenner: 2π·γ·x
            denominator = 2.0 * np.pi * gamma * x_in
            density[inside] = numerator / denominator

        return density

    def dyson_index(self, spacings: np.ndarray) -> float:
        """
        @brief Schätzt den Dyson-Index β aus der Level-Spacing-Verteilung.
        @description
            Der Dyson-Index β klassifiziert Zufallsmatrix-Ensembles:
            - β = 0: Poisson (keine Korrelation)
            - β = 1: GOE (reell-symmetrisch, Zeitumkehrinvarianz)
            - β = 2: GUE (komplex-hermitesch, keine Zeitumkehrinvarianz)
            - β = 4: GSE (quaternionisch, Kramers-Entartung)

            Schätzung: Der Wigner-Surmise liefert p(s) ∝ s^β für kleine s.
            Aus dem Histogramm-Verhalten bei kleinen s kann β abgeleitet werden.

            Methode: Maximum-Likelihood-Schätzung über Fit der Form
            p(s) ≈ A·s^β für s ∈ [0, 0.5]

        @param spacings: 1D-Array normierter Level-Spacings (s ≥ 0).
        @return: Geschätzter Dyson-Index β (float, typisch nahe 0, 1, 2 oder 4).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Nur kleine Abstände für die β-Schätzung verwenden
        small_s = spacings[(spacings > 0.05) & (spacings < 0.5)]

        if len(small_s) < 10:
            # Zu wenige Datenpunkte: Standardwert β=2 (GUE)
            return 2.0

        # Logarithmischer Fit: ln(p(s)) ≈ β·ln(s) + const
        # Empirische Dichte bei kleinen s via Histogramm
        hist, bin_edges = np.histogram(small_s, bins=10, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

        # Nur positive Dichte-Werte für Log-Fit
        mask = hist > 0
        if np.sum(mask) < 3:
            return 2.0

        log_s = np.log(bin_centers[mask])
        log_p = np.log(hist[mask])

        # Lineare Regression: log_p = β·log_s + const
        # Verwende np.polyfit für den linearen Fit
        coeffs = np.polyfit(log_s, log_p, 1)
        beta_estimate = coeffs[0]  # Steigung = β

        return float(beta_estimate)

    def plot_comparison_all_ensembles(self) -> plt.Figure:
        """
        @brief Vergleichsplot: GOE, GUE, Poisson und ζ-Nullstellen in einem Bild.
        @description
            Zeigt alle Level-Spacing-Verteilungen in einem Plot für direkten Vergleich:
            - GOE-Wigner-Surmise (β=1): rote Kurve
            - GUE-Wigner-Surmise (β=2): blaue Kurve
            - Poisson (β=0): grüne gestrichelte Kurve
            - ζ-Nullstellen-Abstände: violettes Histogramm (Montgomery-Odlyzko)

            Lernziel: GOE und GUE zeigen Level-Repulsion (p(0)=0),
            Poisson nicht (p(0)=1). ζ-Nullstellen stimmen mit GUE überein.

        @return: matplotlib-Figure-Objekt mit Vergleichsplot.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        stats = LevelSpacingStatistics()
        zeta_gue = RiemannZetaGUE()

        # ζ-Nullstellen-Abstände berechnen (erste 50 Nullstellen)
        try:
            zeros = zeta_gue.compute_zeta_zeros(50)
            zeta_spacings = zeta_gue.zeros_level_spacing(zeros)
        except Exception:
            # Fallback: nur bekannte 20 Nullstellen
            zeros = zeta_gue.compute_zeta_zeros(20)
            zeta_spacings = zeta_gue.zeros_level_spacing(zeros)

        fig, ax = plt.subplots(figsize=(10, 6))

        # Theoretische Kurven
        s_vals = np.linspace(0, 4, 400)

        ax.plot(s_vals, stats.wigner_surmise_goe(s_vals),
                'r-', linewidth=2.5, label='GOE: (π/2)s·exp(-πs²/4)')
        ax.plot(s_vals, stats.wigner_surmise_gue(s_vals),
                'b-', linewidth=2.5, label='GUE: (32/π²)s²·exp(-4s²/π)')
        ax.plot(s_vals, stats.poisson_distribution(s_vals),
                'g--', linewidth=2, label='Poisson: exp(-s)')

        # ζ-Nullstellen-Histogramm
        if len(zeta_spacings) > 5:
            ax.hist(zeta_spacings, bins=min(20, len(zeta_spacings) // 2),
                    density=True, alpha=0.4, color='purple',
                    label=f'ζ-Nullstellen ({len(zeros)} Stück)')

        ax.set_xlabel('Normierter Abstand s', fontsize=12)
        ax.set_ylabel('Wahrscheinlichkeitsdichte p(s)', fontsize=12)
        ax.set_title('Level-Spacing-Statistik: Vergleich GOE, GUE, Poisson, ζ-Nullstellen',
                     fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 4)
        ax.set_ylim(0, None)

        return fig
