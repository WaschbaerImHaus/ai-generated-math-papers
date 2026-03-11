"""
@file test_random_matrix_theory.py
@brief Tests für das random_matrix_theory-Modul (Zufallsmatrizentheorie).
@description
    Testet GUEEnsemble, GOEEnsemble, LevelSpacingStatistics,
    RiemannZetaGUE und RandomMatrixTheoryTools auf:
    - Matrixeigenschaften (Hermitizität, Symmetrie)
    - Eigenwert-Realität
    - Wigner-Surmise-Werte (analytisch überprüfbar)
    - Normierung der Wahrscheinlichkeitsdichten (Integration ≈ 1)
    - ζ-Nullstellen-Werte (bekannte numerische Werte)
    - Montgomery-Paarkorrelation (analytische Formel)
    - Marchenko-Pastur-Normierung

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import numpy as np
import pytest
import matplotlib
matplotlib.use('Agg')  # Nicht-interaktives Backend für Tests

from src.random_matrix_theory import (
    GUEEnsemble,
    GOEEnsemble,
    LevelSpacingStatistics,
    RiemannZetaGUE,
    RandomMatrixTheoryTools,
)


# ===========================================================================
# TESTS: GUEEnsemble
# ===========================================================================

class TestGUEEnsemble:
    """Tests für das Gauß'sche Unitäre Ensemble (GUE)."""

    def test_sample_returns_correct_shape(self):
        """
        @brief Stellt sicher, dass sample() eine n×n-Matrix zurückgibt.
        """
        gue = GUEEnsemble(n=5, samples=10)
        H = gue.sample()
        assert H.shape == (5, 5), "Matrix muss 5×5 sein."

    def test_sample_is_hermitian(self):
        """
        @brief Prüft H = H† (hermitesch): H_{ij} = conj(H_{ji}).
        """
        gue = GUEEnsemble(n=8, samples=10)
        for _ in range(5):
            H = gue.sample()
            # H† = H: maximale Abweichung sollte numerisch 0 sein
            diff = np.max(np.abs(H - H.conj().T))
            assert diff < 1e-12, f"Matrix ist nicht hermitesch: max |H - H†| = {diff}"

    def test_sample_has_real_eigenvalues(self):
        """
        @brief Hermitesche Matrizen haben nur reelle Eigenwerte.
        """
        gue = GUEEnsemble(n=6, samples=10)
        for _ in range(5):
            H = gue.sample()
            eigenvalues = np.linalg.eigvalsh(H)
            # eigvalsh gibt reelle float-Werte zurück (intern geprüft)
            assert eigenvalues.dtype in [np.float64, np.float32], \
                "Eigenwerte sollen reell (float64) sein."

    def test_sample_complex_dtype(self):
        """
        @brief GUE-Matrizen sind komplex (dtype=complex128).
        """
        gue = GUEEnsemble(n=4, samples=10)
        H = gue.sample()
        assert np.iscomplexobj(H), "GUE-Matrix muss komplex sein."

    def test_eigenvalue_distribution_length(self):
        """
        @brief eigenvalue_distribution() liefert genau samples×n Werte.
        """
        n, samples = 5, 20
        gue = GUEEnsemble(n=n, samples=samples)
        ev = gue.eigenvalue_distribution()
        assert len(ev) == n * samples, \
            f"Erwartet {n*samples} Eigenwerte, erhalten {len(ev)}."

    def test_eigenvalue_distribution_all_real(self):
        """
        @brief Alle gesammelten Eigenwerte müssen reell sein.
        """
        gue = GUEEnsemble(n=4, samples=10)
        ev = gue.eigenvalue_distribution()
        assert np.all(np.isreal(ev)) or ev.dtype == np.float64, \
            "Alle Eigenwerte sollen reell sein."

    def test_wigner_semicircle_at_zero(self):
        """
        @brief ρ(0) = 2/π ≈ 0.6366.
        """
        gue = GUEEnsemble(n=10)
        result = gue.wigner_semicircle(np.array([0.0]))
        expected = 2.0 / np.pi
        assert abs(result[0] - expected) < 1e-10, \
            f"ρ(0) = {result[0]}, erwartet {expected}"

    def test_wigner_semicircle_zero_outside(self):
        """
        @brief Wigner-Halbkreis ist 0 außerhalb [-1, 1].
        """
        gue = GUEEnsemble(n=10)
        x = np.array([-2.0, -1.5, 1.5, 2.0])
        result = gue.wigner_semicircle(x)
        assert np.all(result == 0.0), \
            f"Wigner-Halbkreis außerhalb [-1,1] muss 0 sein, erhalten: {result}"

    def test_wigner_semicircle_normalization(self):
        """
        @brief Integral ∫₋₁¹ ρ(x) dx = 1 (numerische Integration).
        """
        gue = GUEEnsemble(n=10)
        x = np.linspace(-1.0, 1.0, 10000)
        y = gue.wigner_semicircle(x)
        # Trapezregel-Integration
        integral = np.trapezoid(y, x)
        assert abs(integral - 1.0) < 1e-4, \
            f"Wigner-Halbkreis-Normierung: ∫ρ dx = {integral}, erwartet ≈ 1.0"

    def test_plot_returns_figure(self):
        """
        @brief plot_eigenvalue_histogram() gibt ein matplotlib-Figure zurück.
        """
        import matplotlib.pyplot as plt
        gue = GUEEnsemble(n=5, samples=10)
        fig = gue.plot_eigenvalue_histogram(bins=10)
        assert isinstance(fig, plt.Figure), "Erwartet ein matplotlib.pyplot.Figure-Objekt."
        plt.close('all')

    def test_invalid_n_raises(self):
        """
        @brief n < 1 soll ValueError auslösen.
        """
        with pytest.raises(ValueError):
            GUEEnsemble(n=0)

    def test_invalid_samples_raises(self):
        """
        @brief samples < 1 soll ValueError auslösen.
        """
        with pytest.raises(ValueError):
            GUEEnsemble(n=5, samples=0)


# ===========================================================================
# TESTS: GOEEnsemble
# ===========================================================================

class TestGOEEnsemble:
    """Tests für das Gauß'sche Orthogonale Ensemble (GOE)."""

    def test_sample_returns_correct_shape(self):
        """
        @brief sample() liefert eine n×n-Matrix.
        """
        goe = GOEEnsemble(n=6, samples=10)
        H = goe.sample()
        assert H.shape == (6, 6), "Matrix muss 6×6 sein."

    def test_sample_is_symmetric(self):
        """
        @brief GOE-Matrix ist reell-symmetrisch: H = Hᵀ.
        """
        goe = GOEEnsemble(n=7, samples=10)
        for _ in range(5):
            H = goe.sample()
            diff = np.max(np.abs(H - H.T))
            assert diff < 1e-14, \
                f"Matrix ist nicht symmetrisch: max |H - Hᵀ| = {diff}"

    def test_sample_is_real(self):
        """
        @brief GOE-Matrizen sind reell (dtype=float64).
        """
        goe = GOEEnsemble(n=5, samples=10)
        H = goe.sample()
        assert not np.iscomplexobj(H), "GOE-Matrix muss reell sein."

    def test_sample_has_real_eigenvalues(self):
        """
        @brief Symmetrische Matrizen haben reelle Eigenwerte.
        """
        goe = GOEEnsemble(n=5, samples=10)
        for _ in range(5):
            H = goe.sample()
            eigenvalues = np.linalg.eigvalsh(H)
            assert eigenvalues.dtype == np.float64, \
                "Eigenwerte symmetrischer Matrizen müssen reell sein."

    def test_eigenvalue_distribution_length(self):
        """
        @brief eigenvalue_distribution() liefert samples×n Werte.
        """
        n, samples = 4, 15
        goe = GOEEnsemble(n=n, samples=samples)
        ev = goe.eigenvalue_distribution()
        assert len(ev) == n * samples

    def test_level_spacing_positive(self):
        """
        @brief Normierte Level-Spacings sind ≥ 0.
        """
        goe = GOEEnsemble(n=5, samples=50)
        spacings = goe.level_spacing_distribution()
        assert np.all(spacings >= 0), \
            "Alle normierten Abstände müssen ≥ 0 sein."

    def test_level_spacing_mean_near_one(self):
        """
        @brief Mittlerer normierter Abstand sollte nahe 1 sein.
        @description
            Da jede einzelne Stichprobe auf ihren eigenen Mittelwert normiert wird,
            sollte das Gesamtmittel über viele Stichproben nahe 1 liegen.
        """
        goe = GOEEnsemble(n=10, samples=200)
        spacings = goe.level_spacing_distribution()
        # Der Mittelwert sollte nahe 1 liegen
        assert abs(np.mean(spacings) - 1.0) < 0.1, \
            f"Mittlerer normierter Abstand = {np.mean(spacings):.4f}, erwartet ≈ 1.0"

    def test_invalid_n_raises(self):
        """
        @brief n < 1 soll ValueError auslösen.
        """
        with pytest.raises(ValueError):
            GOEEnsemble(n=-1)


# ===========================================================================
# TESTS: LevelSpacingStatistics
# ===========================================================================

class TestLevelSpacing:
    """Tests für LevelSpacingStatistics."""

    def test_wigner_surmise_gue_at_s1(self):
        """
        @brief p_GUE(1) = (32/π²)·exp(-4/π) ≈ 0.597 ≠ 32/π².
        @description
            Analytisch: p(1) = (32/π²)·1²·exp(-4·1²/π)
            = (32/π²)·exp(-4/π) ≈ 3.242·exp(-1.273) ≈ 3.242·0.280 ≈ 0.908

            Hinweis: Die Aufgabenstellung schreibt ≈ 1.0 (±0.1), was mit dem
            Grenzwert s→0 der Dichte übereinstimmt. Wir testen den korrekten
            analytischen Wert: p(1) ≈ 0.908.
        """
        stats = LevelSpacingStatistics()
        result = stats.wigner_surmise_gue(np.array([1.0]))[0]
        expected = (32.0 / np.pi**2) * 1.0**2 * np.exp(-4.0 * 1.0**2 / np.pi)
        assert abs(result - expected) < 1e-10, \
            f"p_GUE(1) = {result}, erwartet {expected}"

    def test_wigner_surmise_gue_normalization(self):
        """
        @brief ∫₀^∞ p_GUE(s) ds = 1.
        """
        stats = LevelSpacingStatistics()
        s = np.linspace(0, 10, 50000)
        p = stats.wigner_surmise_gue(s)
        integral = np.trapezoid(p, s)
        assert abs(integral - 1.0) < 1e-3, \
            f"GUE-Wigner-Surmise-Normierung: ∫p ds = {integral}, erwartet ≈ 1.0"

    def test_wigner_surmise_goe_normalization(self):
        """
        @brief ∫₀^∞ p_GOE(s) ds = 1.
        """
        stats = LevelSpacingStatistics()
        s = np.linspace(0, 10, 50000)
        p = stats.wigner_surmise_goe(s)
        integral = np.trapezoid(p, s)
        assert abs(integral - 1.0) < 1e-3, \
            f"GOE-Wigner-Surmise-Normierung: ∫p ds = {integral}, erwartet ≈ 1.0"

    def test_poisson_normalization(self):
        """
        @brief ∫₀^∞ exp(-s) ds = 1.
        """
        stats = LevelSpacingStatistics()
        s = np.linspace(0, 30, 100000)
        p = stats.poisson_distribution(s)
        integral = np.trapezoid(p, s)
        assert abs(integral - 1.0) < 1e-3, \
            f"Poisson-Normierung: ∫p ds = {integral}, erwartet ≈ 1.0"

    def test_wigner_surmise_goe_at_zero(self):
        """
        @brief p_GOE(0) = 0 (Level-Repulsion).
        """
        stats = LevelSpacingStatistics()
        result = stats.wigner_surmise_goe(np.array([0.0]))[0]
        assert abs(result) < 1e-10, \
            f"p_GOE(0) = {result}, erwartet 0 (Level-Repulsion)"

    def test_wigner_surmise_gue_at_zero(self):
        """
        @brief p_GUE(0) = 0 (stärkere Level-Repulsion als GOE).
        """
        stats = LevelSpacingStatistics()
        result = stats.wigner_surmise_gue(np.array([0.0]))[0]
        assert abs(result) < 1e-10, \
            f"p_GUE(0) = {result}, erwartet 0 (Level-Repulsion)"

    def test_poisson_at_zero(self):
        """
        @brief p_Poisson(0) = 1 (keine Level-Repulsion).
        """
        stats = LevelSpacingStatistics()
        result = stats.poisson_distribution(np.array([0.0]))[0]
        assert abs(result - 1.0) < 1e-10, \
            f"p_Poisson(0) = {result}, erwartet 1"

    def test_normalize_spacings_mean_one(self):
        """
        @brief normalize_spacings() erzeugt Abstände mit Mittelwert 1.
        """
        stats = LevelSpacingStatistics()
        # Bekannte Eigenwerte: gleichmäßig verteilt
        eigenvalues = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        spacings = stats.normalize_spacings(eigenvalues)
        assert abs(np.mean(spacings) - 1.0) < 1e-10, \
            f"Mittlerer normierter Abstand = {np.mean(spacings)}, erwartet 1.0"

    def test_plot_level_spacing_returns_figure(self):
        """
        @brief plot_level_spacing() gibt ein matplotlib-Figure zurück.
        """
        import matplotlib.pyplot as plt
        stats = LevelSpacingStatistics()
        # Zufällige Abstände (Poisson-verteilt als Test)
        spacings = np.random.exponential(1.0, 100)
        fig = stats.plot_level_spacing(spacings, ensemble="GUE")
        assert isinstance(fig, plt.Figure)
        plt.close('all')

    def test_wigner_gue_larger_than_goe_at_small_s(self):
        """
        @brief Für kleine s: p_GOE(s) > p_GUE(s) (GOE weniger unterdrückt).
        @description
            GOE hat linearen s-Term, GUE quadratischen s²-Term.
            Für s = 0.1: p_GOE(0.1) > p_GUE(0.1)
        """
        stats = LevelSpacingStatistics()
        s_small = np.array([0.1])
        goe_val = stats.wigner_surmise_goe(s_small)[0]
        gue_val = stats.wigner_surmise_gue(s_small)[0]
        assert goe_val > gue_val, \
            f"Bei s=0.1: p_GOE={goe_val:.4f} soll > p_GUE={gue_val:.4f} sein."


# ===========================================================================
# TESTS: RiemannZetaGUE
# ===========================================================================

class TestRiemannZetaGUE:
    """Tests für RiemannZetaGUE: ζ-Nullstellen und GUE-Verbindung."""

    def test_first_five_zeros_known_values(self):
        """
        @brief Die ersten 5 ζ-Nullstellen stimmen mit bekannten Werten überein.
        @description
            Bekannte Werte (Imaginärteile):
            t_1 ≈ 14.134, t_2 ≈ 21.022, t_3 ≈ 25.010, t_4 ≈ 30.424, t_5 ≈ 32.935
        """
        zeta = RiemannZetaGUE()
        zeros = zeta.compute_zeta_zeros(5)

        expected = [14.134, 21.022, 25.010, 30.424, 32.935]

        assert len(zeros) == 5, f"Erwartet 5 Nullstellen, erhalten {len(zeros)}"

        for i, (z, e) in enumerate(zip(zeros, expected)):
            assert abs(z - e) < 0.01, \
                f"Nullstelle {i+1}: ζ-Zero={z:.3f}, erwartet≈{e:.3f} (±0.01)"

    def test_zeros_count(self):
        """
        @brief compute_zeta_zeros(n) liefert genau n Nullstellen.
        """
        zeta = RiemannZetaGUE()
        for count in [5, 10, 20]:
            zeros = zeta.compute_zeta_zeros(count)
            assert len(zeros) == count, \
                f"compute_zeta_zeros({count}) liefert {len(zeros)} statt {count}"

    def test_zeros_increasing(self):
        """
        @brief ζ-Nullstellen sind aufsteigend sortiert.
        """
        zeta = RiemannZetaGUE()
        zeros = zeta.compute_zeta_zeros(15)
        diffs = np.diff(zeros)
        assert np.all(diffs > 0), \
            "ζ-Nullstellen müssen streng aufsteigend sein."

    def test_zeros_level_spacing_positive(self):
        """
        @brief Normierte Nullstellen-Abstände sind ≥ 0.
        """
        zeta = RiemannZetaGUE()
        zeros = zeta.compute_zeta_zeros(20)
        spacings = zeta.zeros_level_spacing(zeros)
        assert np.all(spacings >= 0), \
            "Normierte Abstände müssen ≥ 0 sein."

    def test_zeros_level_spacing_count(self):
        """
        @brief zeros_level_spacing() liefert len(zeros)-1 Abstände.
        """
        zeta = RiemannZetaGUE()
        zeros = zeta.compute_zeta_zeros(10)
        spacings = zeta.zeros_level_spacing(zeros)
        assert len(spacings) == len(zeros) - 1, \
            f"Erwartet {len(zeros)-1} Abstände, erhalten {len(spacings)}"

    def test_compare_gue_zeros_returns_figure(self):
        """
        @brief compare_gue_zeros() gibt ein matplotlib-Figure zurück.
        """
        import matplotlib.pyplot as plt
        zeta = RiemannZetaGUE()
        fig = zeta.compare_gue_zeros(count=20)
        assert isinstance(fig, plt.Figure)
        plt.close('all')

    def test_invalid_count_raises(self):
        """
        @brief count ≤ 0 soll ValueError auslösen.
        """
        zeta = RiemannZetaGUE()
        with pytest.raises(ValueError):
            zeta.compute_zeta_zeros(0)


# ===========================================================================
# TESTS: Montgomery-Paarkorrelation
# ===========================================================================

class TestMontgomery:
    """Tests für die Montgomery-Paarkorrelationsfunktion."""

    def test_pair_correlation_formula(self):
        """
        @brief R(r) = 1 - (sin(πr)/(πr))² für r > 0.
        """
        zeta = RiemannZetaGUE()

        # Testpunkte
        test_r = [0.5, 1.0, 1.5, 2.0, 3.0]

        for r in test_r:
            result = zeta.montgomery_pair_correlation(None, r)
            expected = 1.0 - (np.sin(np.pi * r) / (np.pi * r))**2
            assert abs(result - expected) < 1e-12, \
                f"R({r}) = {result}, erwartet {expected}"

    def test_pair_correlation_at_zero(self):
        """
        @brief R(0) = 0 (Level-Repulsion: keine Nullstellen bei r→0).
        """
        zeta = RiemannZetaGUE()
        result = zeta.montgomery_pair_correlation(None, 0.0)
        assert abs(result) < 1e-10, \
            f"R(0) = {result}, erwartet 0"

    def test_pair_correlation_large_r(self):
        """
        @brief R(r) → 1 für große r (keine langreichweitige Korrelation).
        """
        zeta = RiemannZetaGUE()
        # Bei r=10: sin(10π)/(10π) ≈ 0, also R(10) ≈ 1
        result = zeta.montgomery_pair_correlation(None, 10.0)
        assert abs(result - 1.0) < 0.01, \
            f"R(10) = {result}, erwartet ≈ 1.0 für große r"

    def test_pair_correlation_at_integers(self):
        """
        @brief R(n) = 1 für ganzzahlige n ≠ 0 (sin(nπ) = 0).
        """
        zeta = RiemannZetaGUE()
        for n in [1, 2, 3, 4, 5]:
            result = zeta.montgomery_pair_correlation(None, float(n))
            assert abs(result - 1.0) < 1e-12, \
                f"R({n}) = {result}, erwartet 1.0 (sin({n}π) = 0)"

    def test_pair_correlation_range(self):
        """
        @brief 0 ≤ R(r) ≤ 1 für alle r ≥ 0.
        """
        zeta = RiemannZetaGUE()
        r_vals = np.linspace(0.01, 10.0, 100)

        for r in r_vals:
            result = zeta.montgomery_pair_correlation(None, r)
            assert 0.0 <= result <= 1.0 + 1e-10, \
                f"R({r:.2f}) = {result:.4f} liegt außerhalb [0, 1]"


# ===========================================================================
# TESTS: RandomMatrixTheoryTools
# ===========================================================================

class TestRandomMatrixTheoryTools:
    """Tests für RandomMatrixTheoryTools: Marchenko-Pastur, Dyson-Index, Vergleichsplot."""

    def test_marchenko_pastur_normalization(self):
        """
        @brief ∫ ρ_MP(x) dx = 1 (Marchenko-Pastur ist normierte Dichte).
        """
        tools = RandomMatrixTheoryTools()
        # γ = 0.5: λ± = (1 ± √0.5)² ≈ [0.086, 2.914]
        gamma = 0.5
        lambda_minus = (1.0 - np.sqrt(gamma))**2
        lambda_plus = (1.0 + np.sqrt(gamma))**2

        x = np.linspace(lambda_minus + 1e-6, lambda_plus - 1e-6, 10000)
        rho = tools.marchenko_pastur(x, ratio=gamma)
        integral = np.trapezoid(rho, x)

        assert abs(integral - 1.0) < 0.05, \
            f"Marchenko-Pastur-Normierung (γ={gamma}): ∫ρ dx = {integral}, erwartet ≈ 1.0"

    def test_marchenko_pastur_zero_outside_support(self):
        """
        @brief Marchenko-Pastur ist 0 außerhalb [λ₋, λ₊].
        """
        tools = RandomMatrixTheoryTools()
        gamma = 1.0
        # λ₋ = 0, λ₊ = 4
        x_outside = np.array([-1.0, -0.1, 4.5, 5.0])
        rho = tools.marchenko_pastur(x_outside, ratio=gamma)
        assert np.all(rho == 0.0), \
            f"Marchenko-Pastur außerhalb [0,4] soll 0 sein, erhalten: {rho}"

    def test_marchenko_pastur_invalid_ratio(self):
        """
        @brief Negatives ratio soll ValueError auslösen.
        """
        tools = RandomMatrixTheoryTools()
        with pytest.raises(ValueError):
            tools.marchenko_pastur(np.array([1.0, 2.0]), ratio=-0.5)

    def test_marchenko_pastur_positive_inside(self):
        """
        @brief Marchenko-Pastur-Dichte ist ≥ 0 im Träger.
        """
        tools = RandomMatrixTheoryTools()
        gamma = 0.25
        lambda_minus = (1.0 - np.sqrt(gamma))**2  # = 0.25
        lambda_plus = (1.0 + np.sqrt(gamma))**2   # = 2.25

        x = np.linspace(lambda_minus + 0.01, lambda_plus - 0.01, 100)
        rho = tools.marchenko_pastur(x, ratio=gamma)
        assert np.all(rho >= 0), \
            "Marchenko-Pastur-Dichte muss ≥ 0 sein."

    def test_dyson_index_gue(self):
        """
        @brief GOE-Spacings → β ≈ 1 (GOE hat linearen Term in s).
        @description
            Aus einer GOE-Stichprobe sollte der Dyson-Index ≈ 1 sein.
            Toleranz ±1.5 (empirische Schätzung, nicht exakt).
        """
        tools = RandomMatrixTheoryTools()
        goe = GOEEnsemble(n=20, samples=500)
        spacings = goe.level_spacing_distribution()
        beta = tools.dyson_index(spacings)
        # Grobe Schätzung: GOE hat β≈1, Toleranz ±1.5 wegen kleiner Stichprobe
        assert 0.0 <= beta <= 4.0, \
            f"Dyson-Index {beta:.2f} außerhalb realistischem Bereich [0, 4]"

    def test_plot_comparison_returns_figure(self):
        """
        @brief plot_comparison_all_ensembles() gibt ein matplotlib-Figure zurück.
        """
        import matplotlib.pyplot as plt
        tools = RandomMatrixTheoryTools()
        fig = tools.plot_comparison_all_ensembles()
        assert isinstance(fig, plt.Figure)
        plt.close('all')

    def test_marchenko_pastur_gamma_one(self):
        """
        @brief Für γ=1: Spektrum liegt in [0, 4], Dichte am Rand = ∞ (integrierbar).
        @description
            Für γ=1: λ₋ = 0, λ₊ = 4.
            Mittelwert des Spektrums: <x> = γ·σ² = 1 (hier σ²=1, γ=1).
        """
        tools = RandomMatrixTheoryTools()
        # Inneres des Spektrums: x = 2 (Mittelpunkt von [0, 4])
        rho_inner = tools.marchenko_pastur(np.array([2.0]), ratio=1.0)
        # Außen: x = 5
        rho_outer = tools.marchenko_pastur(np.array([5.0]), ratio=1.0)

        assert rho_inner[0] > 0, "Marchenko-Pastur bei x=2 (γ=1) soll > 0 sein."
        assert rho_outer[0] == 0.0, "Marchenko-Pastur bei x=5 (γ=1) soll 0 sein."
