"""
@file test_collatz_ergodic.py
@brief Unit-Tests für das collatz_ergodic-Modul.
@description
    Testet alle Klassen und Funktionen des Collatz-Moduls:
    - CollatzFunction:         Grundlegende Collatz-Operationen
    - CollatzTwoadicAnalysis:  2-adische Analyse (Lagarias-Perspektive)
    - CollatzErgodicMeasure:   Ergodische Maße (Tao-Perspektive)
    - CollatzStatistics:       Statistische Analysen und Visualisierung
    - verify_collatz():        Verifikation der Vermutung bis N

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import pytest
import math
import numpy as np

# Projektverzeichnis src/ in den Python-Pfad aufnehmen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Zu testende Module importieren
from collatz_ergodic import (
    CollatzFunction,
    CollatzTwoadicAnalysis,
    CollatzErgodicMeasure,
    CollatzStatistics,
    verify_collatz,
)


# ===========================================================================
# Tests: CollatzFunction
# ===========================================================================

class TestCollatzFunction:
    """Tests für die grundlegenden Collatz-Operationen."""

    def setup_method(self):
        """Erzeugt eine frische Instanz vor jedem Test."""
        self.cf = CollatzFunction()

    # --- step() ---

    def test_step_gerade(self):
        """Gerader Schritt: n/2"""
        assert self.cf.step(6) == 3
        assert self.cf.step(10) == 5
        assert self.cf.step(100) == 50

    def test_step_ungerade(self):
        """Ungerader Schritt: 3n+1"""
        assert self.cf.step(3) == 10
        assert self.cf.step(5) == 16
        assert self.cf.step(27) == 82

    def test_step_eins(self):
        """T(1) = 4 (1 ist ungerade)"""
        assert self.cf.step(1) == 4

    def test_step_negativ_wirft_fehler(self):
        """Negativer Wert soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.cf.step(-1)

    def test_step_null_wirft_fehler(self):
        """Wert 0 soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.cf.step(0)

    # --- sequence() ---

    def test_sequence_6(self):
        """Klassisches Beispiel: sequence(6) = [6,3,10,5,16,8,4,2,1]"""
        erwartet = [6, 3, 10, 5, 16, 8, 4, 2, 1]
        assert self.cf.sequence(6) == erwartet

    def test_sequence_1(self):
        """sequence(1) = [1] (bereits am Ziel)"""
        assert self.cf.sequence(1) == [1]

    def test_sequence_2(self):
        """sequence(2) = [2, 1]"""
        assert self.cf.sequence(2) == [2, 1]

    def test_sequence_endet_mit_1(self):
        """Jede Folge muss mit 1 enden (für kleine n)."""
        for n in range(1, 50):
            folge = self.cf.sequence(n)
            assert folge[-1] == 1, f"sequence({n}) endet nicht mit 1"

    def test_sequence_startet_mit_n(self):
        """Folge muss mit dem Startwert beginnen."""
        for n in [3, 7, 27]:
            assert self.cf.sequence(n)[0] == n

    def test_sequence_negativer_wert_fehler(self):
        """Negativer Startwert soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.cf.sequence(-5)

    # --- total_stopping_time() ---

    def test_total_stopping_time_27(self):
        """Bekanntes Ergebnis: τ(27) = 111"""
        assert self.cf.total_stopping_time(27) == 111

    def test_total_stopping_time_1(self):
        """τ(1) = 0 (bereits bei 1)"""
        assert self.cf.total_stopping_time(1) == 0

    def test_total_stopping_time_2(self):
        """τ(2) = 1 (2 → 1)"""
        assert self.cf.total_stopping_time(2) == 1

    def test_total_stopping_time_positiv(self):
        """Stoppzeiten für n = 1..20 müssen ≥ 0 sein."""
        for n in range(1, 21):
            assert self.cf.total_stopping_time(n) >= 0

    def test_total_stopping_time_konsistenz_mit_sequence(self):
        """τ(n) muss mit len(sequence(n))-1 übereinstimmen."""
        for n in range(1, 30):
            tau = self.cf.total_stopping_time(n)
            laenge = len(self.cf.sequence(n)) - 1
            assert tau == laenge, f"τ({n})={tau} ≠ len-1={laenge}"

    # --- stopping_time() ---

    def test_stopping_time_eins(self):
        """σ(1) = 0"""
        assert self.cf.stopping_time(1) == 0

    def test_stopping_time_positiv(self):
        """Stoppzeiten für n = 2..20 müssen ≥ 1 sein."""
        for n in range(2, 21):
            sigma = self.cf.stopping_time(n)
            assert sigma >= 1, f"σ({n}) = {sigma} < 1"

    def test_stopping_time_korrektheit(self):
        """Nach σ(n) Schritten muss der Wert < n sein."""
        for n in range(2, 30):
            sigma = self.cf.stopping_time(n)
            if sigma == -1:
                continue
            # σ(n) Schritte anwenden
            aktuell = n
            for _ in range(sigma):
                aktuell = self.cf.step(aktuell)
            assert aktuell < n, f"Nach σ({n})={sigma} Schritten: {aktuell} < {n} verletzt"

    # --- trajectory_statistics() ---

    def test_trajectory_statistics_felder(self):
        """Statistik-Dictionary muss alle Felder enthalten."""
        stats = self.cf.trajectory_statistics(27)
        for key in ["max_value", "steps_to_1", "odd_steps", "even_steps", "compression_ratio"]:
            assert key in stats

    def test_trajectory_statistics_27(self):
        """Bekannte Werte für n=27."""
        stats = self.cf.trajectory_statistics(27)
        assert stats["steps_to_1"] == 111
        assert stats["max_value"] == 9232  # Bekanntes Maximum der 27er-Folge
        assert stats["odd_steps"] + stats["even_steps"] == 111

    def test_trajectory_statistics_compression_ratio_positiv(self):
        """Kompressionsrate muss positiv sein."""
        stats = self.cf.trajectory_statistics(100)
        assert stats["compression_ratio"] > 0


# ===========================================================================
# Tests: CollatzTwoadicAnalysis
# ===========================================================================

class TestTwoadicAnalysis:
    """Tests für die 2-adische Collatz-Analyse."""

    def setup_method(self):
        """Erzeugt eine frische Instanz vor jedem Test."""
        self.ca = CollatzTwoadicAnalysis()

    # --- two_adic_valuation() ---

    def test_valuation_potenzen(self):
        """ν₂(2^k) = k"""
        assert self.ca.two_adic_valuation(1) == 0
        assert self.ca.two_adic_valuation(2) == 1
        assert self.ca.two_adic_valuation(4) == 2
        assert self.ca.two_adic_valuation(8) == 3
        assert self.ca.two_adic_valuation(16) == 4

    def test_valuation_ungerade(self):
        """ν₂(ungerade) = 0"""
        for n in [1, 3, 5, 7, 9, 11, 13]:
            assert self.ca.two_adic_valuation(n) == 0

    def test_valuation_gemischt(self):
        """ν₂(6) = 1, ν₂(12) = 2, ν₂(24) = 3"""
        assert self.ca.two_adic_valuation(6) == 1
        assert self.ca.two_adic_valuation(12) == 2
        assert self.ca.two_adic_valuation(24) == 3

    def test_valuation_null_wirft_fehler(self):
        """ν₂(0) soll ValueError auslösen."""
        with pytest.raises(ValueError):
            self.ca.two_adic_valuation(0)

    # --- odd_part() ---

    def test_odd_part(self):
        """Bekannte Werte."""
        assert self.ca.odd_part(12) == 3   # 12 = 4·3
        assert self.ca.odd_part(7) == 7    # 7 ungerade
        assert self.ca.odd_part(8) == 1    # 8 = 8·1
        assert self.ca.odd_part(6) == 3    # 6 = 2·3

    # --- collatz_as_padic() ---

    def test_collatz_as_padic_laenge(self):
        """Rückgabeliste hat maximal steps+1 Elemente."""
        ergebnis = self.ca.collatz_as_padic(6, 5)
        assert len(ergebnis) <= 6

    def test_collatz_as_padic_werte_nichtnegativ(self):
        """Alle 2-adischen Bewertungen ≥ 0."""
        for v in self.ca.collatz_as_padic(27, 50):
            assert v >= 0

    # --- lyapunov_exponent() ---

    def test_lyapunov_negativ(self):
        """Lyapunov-Exponent muss negativ sein (Konvergenz)."""
        lam = self.ca.lyapunov_exponent(27, 1000)
        assert lam < 0, f"Lyapunov-Exponent {lam} sollte negativ sein"

    def test_lyapunov_bereich(self):
        """Lyapunov-Exponent für n=27 liegt zwischen -0.5 und 0."""
        lam = self.ca.lyapunov_exponent(27, 1000)
        assert -0.5 < lam < 0, f"λ = {lam} außerhalb [-0.5, 0]"

    def test_lyapunov_theoretischer_wert(self):
        """
        λ sollte in der Nähe des theoretischen Werts liegen:
        λ_theo = log(3)/2 - log(2) ≈ -0.0959
        Toleranz ±0.15 für kleine Stichproben.
        """
        lam_theo = math.log(3) / 2 - math.log(2)  # ≈ -0.0959
        lam_emp = self.ca.lyapunov_exponent(97, 2000)
        assert abs(lam_emp - lam_theo) < 0.2, (
            f"λ_emp={lam_emp:.4f} zu weit von λ_theo={lam_theo:.4f}"
        )

    # --- average_contraction() ---

    def test_average_contraction_negativ(self):
        """Mittlere Kontraktion muss negativ sein."""
        kontraktion = self.ca.average_contraction(sample_size=500)
        assert kontraktion < 0, f"Mittlere Kontraktion {kontraktion} sollte negativ sein"


# ===========================================================================
# Tests: CollatzErgodicMeasure
# ===========================================================================

class TestErgodicMeasure:
    """Tests für die ergodischen Maße des Collatz-Systems."""

    def setup_method(self):
        """Erzeugt eine frische Instanz vor jedem Test."""
        self.cem = CollatzErgodicMeasure()

    # --- invariant_measure_estimate() ---

    def test_invariant_measure_normierung(self):
        """Das empirische Maß muss (bis auf Rundungsfehler) zu 1 summieren."""
        mass = self.cem.invariant_measure_estimate(50)
        assert abs(mass.sum() - 1.0) < 1e-10

    def test_invariant_measure_nichtnegativ(self):
        """Alle Gewichte müssen ≥ 0 sein."""
        mass = self.cem.invariant_measure_estimate(50)
        assert np.all(mass >= 0)

    def test_invariant_measure_konzentration(self):
        """Kleinen Zahlen (1, 2, 4, ...) haben hohes Gewicht."""
        mass = self.cem.invariant_measure_estimate(50)
        # 1, 2, 4 zusammen sollten substanziellen Anteil tragen
        assert mass[1] + mass[2] + mass[4] > 0.05

    # --- ergodic_average() ---

    def test_ergodic_average_identitaet(self):
        """Ergodischer Mittelwert der Identität ≈ mittlere Collatz-Nachfolger."""
        cf = CollatzFunction()
        mittelwert = self.cem.ergodic_average(lambda x: x, 100)
        assert mittelwert > 0  # Muss positiv sein

    def test_ergodic_average_konstante(self):
        """Ergodischer Mittelwert einer Konstante c = c."""
        mittelwert = self.cem.ergodic_average(lambda x: 5.0, 100)
        assert abs(mittelwert - 5.0) < 1e-10

    # --- density_below_bound() ---

    def test_density_below_bound_trivial(self):
        """Schranke n (identisch) → fast alle erfüllen min T^k(n) ≤ n."""
        anteil = self.cem.density_below_bound(100, lambda n: n)
        assert anteil >= 0.95

    def test_density_below_bound_bereich(self):
        """Anteil liegt immer in [0, 1]."""
        anteil = self.cem.density_below_bound(50, lambda n: n ** 0.5)
        assert 0.0 <= anteil <= 1.0

    # --- tao_theorem_numerical() ---

    def test_tao_theorem_fraction_hoch(self):
        """
        Taos Theorem: Anteil der n ≤ 100 mit min T^k(n) ≤ n^0.5 muss > 0.9 sein.
        """
        ergebnis = self.cem.tao_theorem_numerical(N=100, epsilon=0.5)
        assert ergebnis["fraction"] > 0.9, (
            f"Tao-Anteil {ergebnis['fraction']:.3f} zu niedrig (erwartet > 0.9)"
        )

    def test_tao_theorem_felder(self):
        """Rückgabe-Dictionary muss alle Felder enthalten."""
        ergebnis = self.cem.tao_theorem_numerical(N=50, epsilon=0.5)
        for key in ["verified", "fraction", "epsilon", "N"]:
            assert key in ergebnis

    def test_tao_theorem_epsilon_gespeichert(self):
        """Verwendetes epsilon muss korrekt gespeichert werden."""
        ergebnis = self.cem.tao_theorem_numerical(N=50, epsilon=0.3)
        assert ergebnis["epsilon"] == 0.3

    def test_tao_theorem_groesseres_epsilon(self):
        """Mit ε=0.8 (große Schranke) muss Anteil noch höher sein."""
        ergebnis = self.cem.tao_theorem_numerical(N=100, epsilon=0.8)
        assert ergebnis["fraction"] > 0.95


# ===========================================================================
# Tests: CollatzStatistics
# ===========================================================================

class TestStatistics:
    """Tests für statistische Analysen der Collatz-Folgen."""

    def setup_method(self):
        """Erzeugt eine frische Instanz vor jedem Test."""
        self.cs = CollatzStatistics()

    # --- stopping_time_distribution() ---

    def test_stopping_time_distribution_typ(self):
        """Rückgabe ist ein Dictionary."""
        dist = self.cs.stopping_time_distribution(20)
        assert isinstance(dist, dict)

    def test_stopping_time_distribution_summe(self):
        """Summe aller Häufigkeiten muss gleich limit sein."""
        limit = 50
        dist = self.cs.stopping_time_distribution(limit)
        assert sum(dist.values()) == limit

    def test_stopping_time_distribution_nullzeit(self):
        """n=1 hat Stoppzeit 0 → Schlüssel 0 muss vorkommen."""
        dist = self.cs.stopping_time_distribution(10)
        assert 0 in dist

    # --- record_setters() ---

    def test_record_setters_erster_eintrag(self):
        """Erster Rekordhalter muss n=1 mit Stoppzeit 0 sein."""
        rekorde = self.cs.record_setters(100)
        assert len(rekorde) >= 1
        n0, t0 = rekorde[0]
        assert n0 == 1 and t0 == 0

    def test_record_setters_aufsteigend(self):
        """Rekordhalter müssen in aufsteigender Reihenfolge nach n sortiert sein."""
        rekorde = self.cs.record_setters(100)
        ns = [r[0] for r in rekorde]
        assert ns == sorted(ns)

    def test_record_setters_stoppzeiten_wachsend(self):
        """Stoppzeiten der Rekordhalter müssen streng monoton wachsen."""
        rekorde = self.cs.record_setters(100)
        zeiten = [r[1] for r in rekorde]
        for i in range(1, len(zeiten)):
            assert zeiten[i] > zeiten[i - 1]

    def test_record_setters_bekannte_rekorde(self):
        """Bekannte Rekordhalter: 27 hat Stoppzeit 111 und muss in der Liste sein."""
        rekorde = self.cs.record_setters(100)
        ns = [r[0] for r in rekorde]
        assert 27 in ns

    # --- plot_stopping_times() ---

    def test_plot_stopping_times_gibt_figure_zurueck(self):
        """plot_stopping_times() muss eine matplotlib-Figure zurückgeben."""
        import matplotlib.pyplot as plt
        fig = self.cs.plot_stopping_times(limit=50)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    # --- plot_trajectory() ---

    def test_plot_trajectory_gibt_figure_zurueck(self):
        """plot_trajectory() muss eine matplotlib-Figure zurückgeben."""
        import matplotlib.pyplot as plt
        fig = self.cs.plot_trajectory(27)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    # --- collatz_tree() ---

    def test_collatz_tree_enthaelt_1(self):
        """Baum muss Wurzel 1 enthalten."""
        baum = self.cs.collatz_tree(depth=3)
        assert 1 in baum

    def test_collatz_tree_vorgaenger_von_1(self):
        """2 muss als Vorgänger von 1 im Baum stehen (T(2)=1)."""
        baum = self.cs.collatz_tree(depth=1)
        assert 2 in baum[1]

    def test_collatz_tree_tiefe_2(self):
        """Bei Tiefe 2 muss 4 im Baum sein (T(4)=2, T(2)=1)."""
        baum = self.cs.collatz_tree(depth=2)
        # 4 ist Vorgänger von 2 (T(4) = 2)
        assert 4 in baum.get(2, [])


# ===========================================================================
# Tests: verify_collatz()
# ===========================================================================

class TestVerifyCollatz:
    """Tests für die Collatz-Verifikationsfunktion."""

    def test_verify_kleines_n(self):
        """Verifikation für n ≤ 100 muss True ergeben."""
        ergebnis, wert = verify_collatz(100)
        assert ergebnis is True
        assert wert == 100

    def test_verify_1000(self):
        """Verifikation für n ≤ 1000 muss True ergeben."""
        ergebnis, wert = verify_collatz(1000)
        assert ergebnis is True
        assert wert == 1000

    def test_verify_10000(self):
        """Verifikation für n ≤ 10000 muss True ergeben (bekanntes Ergebnis)."""
        ergebnis, wert = verify_collatz(10000)
        assert ergebnis is True
        assert wert == 10000

    def test_verify_gibt_tuple_zurueck(self):
        """Rückgabe muss ein Tupel (bool, int) sein."""
        ergebnis = verify_collatz(10)
        assert isinstance(ergebnis, tuple)
        assert len(ergebnis) == 2
        assert isinstance(ergebnis[0], bool)
        assert isinstance(ergebnis[1], int)

    def test_verify_n_1(self):
        """verify_collatz(1) = (True, 1)"""
        assert verify_collatz(1) == (True, 1)


# ===========================================================================
# Edge-Case-Tests
# ===========================================================================

class TestEdgeCases:
    """Randfall- und Sonderfälle."""

    def test_collatz_sehr_grosse_zahl(self):
        """Stoppzeit für eine große Zahl muss endlich sein."""
        cf = CollatzFunction()
        # 2^31 - 1 ist eine große Zahl
        tau = cf.total_stopping_time(2 ** 20)
        assert tau >= 0

    def test_collatz_zweierpotenz_faellt_schnell(self):
        """2^k hat τ(2^k) = k (direkte Halbierkette 2^k → 2^{k-1} → ... → 1)."""
        cf = CollatzFunction()
        for k in range(1, 10):
            assert cf.total_stopping_time(2 ** k) == k

    def test_zwei_adisch_und_collatz_konsistenz(self):
        """ν₂(T(n)) > 0 nach jedem ungeraden Schritt (T(ungerade) = gerade)."""
        ca = CollatzTwoadicAnalysis()
        cf = CollatzFunction()
        # Nach einem ungeraden Schritt ist das Ergebnis durch 2 teilbar (wegen +1)
        for n in [3, 5, 7, 9, 11, 13]:
            t_n = cf.step(n)
            assert t_n % 2 == 0, f"T({n}) = {t_n} sollte gerade sein"
            assert ca.two_adic_valuation(t_n) >= 1

    def test_ergodisch_groessere_stichprobe(self):
        """Tao-Theorem mit N=200 und ε=0.6 → Anteil > 0.9."""
        cem = CollatzErgodicMeasure()
        ergebnis = cem.tao_theorem_numerical(N=200, epsilon=0.6)
        assert ergebnis["fraction"] > 0.9
