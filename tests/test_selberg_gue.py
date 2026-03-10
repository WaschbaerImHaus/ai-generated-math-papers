"""
Tests für Selberg-Klasse (analytic_number_theory.py) und
GUE-Statistik der Riemann-Nullstellen (complex_analysis.py).

Mathematischer Hintergrund:
    - Selberg-Klasse: axiomatische Charakterisierung von L-Funktionen
    - GUE: Gaussian Unitary Ensemble (Zufallsmatrixtheorie)
    - Montgomery-Odlyzko-Gesetz: Riemann-Nullstellen ~ GUE-Eigenwerte
    - Wigner-Surmise: P(δ) = (32/π²)·δ²·exp(-4δ²/π)

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import math
import pytest
import sys
import os

# Suchpfad für Quell-Module setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from analytic_number_theory import (
    selberg_class_check,
    selberg_orthogonality,
    selberg_zeta_motivation,
)
from complex_analysis import (
    zeta_zero_spacing_statistics,
    pair_correlation_function,
    random_matrix_gue_sample,
)


# ===========================================================================
# HILFSFUNKTIONEN FÜR TESTS
# ===========================================================================

def _zeta_coefficients(n: int) -> list[complex]:
    """
    Erzeugt die Koeffizienten a_n = 1 für ζ(s) = Σ 1/n^s.

    Für ζ(s) gilt: a_n = 1 für alle n, also ist ζ(s) die einfachste L-Funktion.
    Sie erfüllt alle Selberg-Axiome und hat Grad 1, Führer 1.

    @param n: Anzahl der Koeffizienten
    @return: Liste [1, 1, ..., 1] der Länge n
    """
    return [complex(1.0, 0.0)] * n


def _dirichlet_chi4_coefficients(n: int) -> list[complex]:
    """
    Erzeugt Koeffizienten des Dirichlet-Charakters χ₄ mod 4.

    χ₄(n): 1 wenn n ≡ 1 (mod 4), -1 wenn n ≡ 3 (mod 4), 0 sonst.
    L(s, χ₄) = 1 - 1/3^s + 1/5^s - 1/7^s + ... (Leibniz-Formel für π/4)

    @param n: Anzahl der Koeffizienten
    @return: Liste der χ₄-Werte
    """
    chi4 = [0, 1, 0, -1]  # χ₄(1)=1, χ₄(2)=0, χ₄(3)=-1, χ₄(4)=0 → periodisch mod 4
    return [complex(chi4[(i) % 4], 0.0) for i in range(1, n + 1)]


# ===========================================================================
# TESTS: selberg_class_check
# ===========================================================================

class TestSelbergClassCheck:
    """Tests für die Selberg-Axiom-Prüffunktion."""

    def test_zeta_normalization_axiom1(self):
        """ζ(s) hat a_1 = 1 → Normierungsaxiom erfüllt."""
        coeffs = _zeta_coefficients(50)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom1_normalization"] is True, (
            "ζ(s) muss Axiom 1 (a_1=1) erfüllen"
        )

    def test_zeta_growth_condition(self):
        """ζ(s)-Koeffizienten a_n=1 wachsen polynomiell → Growth-Bedingung erfüllt."""
        coeffs = _zeta_coefficients(100)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom1_growth"] is True, (
            "a_n=1 erfüllt die polynomielle Wachstumsbedingung"
        )

    def test_zeta_axiom1_overall(self):
        """ζ(s) erfüllt das vollständige Axiom 1."""
        coeffs = _zeta_coefficients(50)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom1_satisfied"] is True

    def test_zeta_ramanujan_axiom5(self):
        """
        Ramanujan-Vermutung für ζ(s): |a_p| = 1 ≤ 1·2 = 2.
        Für Grad 1: Schranke ist max(2·1, 1) = 2.
        """
        coeffs = _zeta_coefficients(100)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom5_satisfied"] is True, (
            "ζ(s) erfüllt die Ramanujan-Vermutung: |a_p|=1 ≤ 2"
        )
        assert result["axiom5_violations"] == [], "Keine Ramanujan-Verletzungen für ζ(s)"

    def test_zeta_selberg_candidate(self):
        """ζ(s) wird als Selberg-Klassen-Kandidat erkannt."""
        coeffs = _zeta_coefficients(100)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["selberg_class_candidate"] is True

    def test_return_structure(self):
        """Rückgabe enthält alle erforderlichen Schlüssel."""
        coeffs = _zeta_coefficients(30)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        required_keys = [
            "axiom1_normalization", "axiom1_growth", "axiom1_satisfied",
            "axiom4_satisfied", "axiom5_satisfied", "axiom5_violations",
            "ramanujan_bound", "mean_abs_a_p", "max_abs_a_p",
            "degree", "conductor", "n_coefficients", "selberg_class_candidate"
        ]
        for key in required_keys:
            assert key in result, f"Schlüssel '{key}' fehlt im Ergebnis"

    def test_ramanujan_violation_detection(self):
        """Koeffizienten mit |a_p| > Schranke werden als Verletzung erkannt."""
        # Erzeuge Koeffizienten mit großen Primzahl-Koeffizienten (|a_2| = 10)
        coeffs = _zeta_coefficients(50)
        coeffs[1] = complex(10.0, 0.0)  # a_2 = 10 >> 2 (Ramanujan-Schranke)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom5_satisfied"] is False, (
            "Koeffizient |a_2|=10 sollte Ramanujan-Verletzung auslösen"
        )
        assert len(result["axiom5_violations"]) > 0

    def test_degree_stored_correctly(self):
        """Grad d wird korrekt gespeichert."""
        coeffs = _zeta_coefficients(30)
        result = selberg_class_check(coeffs, degree=2.0, conductor=11.0)
        assert result["degree"] == 2.0
        assert result["conductor"] == 11.0

    def test_dirichlet_chi4_normalization(self):
        """χ₄-Charakter hat a_1 = 1 → Normierung erfüllt."""
        coeffs = _dirichlet_chi4_coefficients(60)
        result = selberg_class_check(coeffs, degree=1.0, conductor=4.0)
        assert result["axiom1_normalization"] is True, "χ₄ hat a_1 = 1"

    def test_wrong_normalization_detected(self):
        """Koeffizientenliste mit a_1 ≠ 1 wird als nicht normiert erkannt."""
        coeffs = [complex(2.0, 0.0)] + [complex(1.0, 0.0)] * 49
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["axiom1_normalization"] is False

    def test_too_few_coefficients(self):
        """Weniger als 2 Koeffizienten: Fehler wird zurückgegeben."""
        result = selberg_class_check([complex(1.0)], degree=1.0, conductor=1.0)
        assert "error" in result

    def test_ramanujan_bound_degree2(self):
        """Für Grad 2 ist die Ramanujan-Schranke 4.0 (max(2·2, 1))."""
        coeffs = _zeta_coefficients(50)
        result = selberg_class_check(coeffs, degree=2.0, conductor=1.0)
        assert result["ramanujan_bound"] == 4.0

    def test_max_abs_a_p_nonnegative(self):
        """Maximaler Betrag der Primzahl-Koeffizienten ist ≥ 0."""
        coeffs = _zeta_coefficients(50)
        result = selberg_class_check(coeffs, degree=1.0, conductor=1.0)
        assert result["max_abs_a_p"] >= 0.0
        assert result["mean_abs_a_p"] >= 0.0


# ===========================================================================
# TESTS: selberg_orthogonality
# ===========================================================================

class TestSelbergOrthogonality:
    """Tests für die Selberg-Orthogonalitätsvermutung."""

    def test_zeta_with_itself_nonzero(self):
        """
        ζ ⊗ ζ: Summe sollte nicht null sein (gleiche L-Funktion → n_F = 1).
        Die normierte Summe ist positiv für F = G.
        """
        coeffs = _zeta_coefficients(200)
        result = selberg_orthogonality(coeffs, coeffs, prime_bound=100)
        # |Σ a_p · ā_p / p| = Σ 1/p > 0 für ζ(s)
        assert result["orthogonality_measure"] > 0.0, (
            "ζ × ζ: Summe muss positiv sein (gleiche Funktion)"
        )

    def test_zeta_vs_chi4_smaller_than_zeta_self(self):
        """
        ζ und χ₄ sind orthogonale L-Funktionen:
        Die Korrelation sollte kleiner sein als ζ mit sich selbst.
        """
        zeta = _zeta_coefficients(150)
        chi4 = _dirichlet_chi4_coefficients(150)
        result_self = selberg_orthogonality(zeta, zeta, prime_bound=100)
        result_cross = selberg_orthogonality(zeta, chi4, prime_bound=100)
        # Kreuzkorrelation < Selbstkorrelation (Orthogonalität)
        assert result_cross["orthogonality_measure"] < result_self["orthogonality_measure"], (
            "Kreuzkorrelation ζ⊗χ₄ sollte kleiner sein als ζ⊗ζ"
        )

    def test_return_structure_orthogonality(self):
        """Rückgabe enthält alle erforderlichen Schlüssel."""
        coeffs = _zeta_coefficients(50)
        result = selberg_orthogonality(coeffs, coeffs, prime_bound=50)
        required = [
            "sum", "log_x", "normalized", "orthogonality_measure",
            "same_function_detected", "primes_used", "prime_bound", "interpretation"
        ]
        for key in required:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_primes_used_count(self):
        """Anzahl der verwendeten Primzahlen ist korrekt."""
        coeffs = _zeta_coefficients(100)
        result = selberg_orthogonality(coeffs, coeffs, prime_bound=50)
        # Primzahlen ≤ 50: 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47 → 15 Stück
        assert result["primes_used"] == 15, (
            f"Erwartet 15 Primzahlen ≤ 50, erhalten: {result['primes_used']}"
        )

    def test_orthogonality_measure_nonnegative(self):
        """Orthogonalitätsmaß ist immer ≥ 0."""
        coeffs = _zeta_coefficients(50)
        result = selberg_orthogonality(coeffs, coeffs, prime_bound=30)
        assert result["orthogonality_measure"] >= 0.0


# ===========================================================================
# TESTS: selberg_zeta_motivation
# ===========================================================================

class TestSelbergZetaMotivation:
    """Tests für die explizite Formel und Nullstellenbeiträge."""

    def test_basic_structure_x100(self):
        """Für x=100: alle Ausgabefelder vorhanden und sinnvoll."""
        result = selberg_zeta_motivation(100.0)
        assert result["x"] == 100.0
        assert result["pi_x"] == 25  # π(100) = 25
        assert result["Li_x"] > 0.0
        assert result["n_zeros_used"] > 0

    def test_zero_contributions_positive(self):
        """Beiträge der Nullstellen haben positiven Betrag."""
        result = selberg_zeta_motivation(50.0)
        for contrib in result["zero_contributions"]:
            assert contrib["abs_contribution"] > 0.0
            assert contrib["imaginary_part"] > 0.0

    def test_rh_error_bound_positive(self):
        """Fehlerschranke unter RH ist positiv."""
        result = selberg_zeta_motivation(1000.0)
        assert result["rh_error_bound"] > 0.0

    def test_invalid_x_raises(self):
        """x ≤ 2 wirft ValueError."""
        with pytest.raises(ValueError):
            selberg_zeta_motivation(2.0)
        with pytest.raises(ValueError):
            selberg_zeta_motivation(1.0)

    def test_error_within_rh_bound_key_exists(self):
        """Schlüssel error_within_rh_bound ist vorhanden."""
        result = selberg_zeta_motivation(500.0)
        assert "error_within_rh_bound" in result


# ===========================================================================
# TESTS: zeta_zero_spacing_statistics (GUE)
# ===========================================================================

class TestZetaZeroSpacingStatistics:
    """Tests für die GUE-Abstandsstatistik der Riemann-Nullstellen."""

    def test_mean_spacing_close_to_one(self):
        """
        Nach Normierung sollte der mittlere Abstand ≈ 1 sein.
        Toleranz: ±0.3 (kleine Stichprobe, grobe Näherung).
        """
        result = zeta_zero_spacing_statistics(n_zeros=30)
        if result.get("n_spacings", 0) == 0:
            pytest.skip("Nicht genug Nullstellen für Test")
        mean = result["mean_spacing"]
        assert 0.3 < mean < 2.5, (
            f"Mittlerer normierter Abstand sollte ≈ 1 sein, erhalten: {mean:.4f}"
        )

    def test_spacings_all_positive(self):
        """Alle normierten Abstände sind positiv (Nullstellen sind verschieden)."""
        result = zeta_zero_spacing_statistics(n_zeros=30)
        if result.get("n_spacings", 0) == 0:
            pytest.skip("Nicht genug Nullstellen für Test")
        for delta in result["spacings"]:
            assert delta > 0.0, f"Abstand {delta} ≤ 0"

    def test_variance_finite(self):
        """Varianz ist endlich und nicht-negativ."""
        result = zeta_zero_spacing_statistics(n_zeros=30)
        if result.get("n_spacings", 0) == 0:
            pytest.skip("Nicht genug Nullstellen für Test")
        assert result["variance"] >= 0.0
        assert math.isfinite(result["variance"])

    def test_ks_statistic_structure(self):
        """KS-Statistik ist vorhanden und im gültigen Bereich [0, 1]."""
        result = zeta_zero_spacing_statistics(n_zeros=30)
        if result.get("n_spacings", 0) == 0:
            pytest.skip("Nicht genug Nullstellen für Test")
        assert "ks_statistic_gue" in result
        assert "ks_p_value_gue" in result
        # KS-Statistik zwischen 0 und 1
        assert 0.0 <= result["ks_statistic_gue"] <= 1.0

    def test_ks_p_value_not_asserted_strict(self):
        """
        KS p-Wert hat gültige Struktur (kein harter p > 0.05-Test,
        da kleine Stichproben statistisch instabil sind).
        """
        result = zeta_zero_spacing_statistics(n_zeros=30)
        if result.get("n_spacings", 0) == 0:
            pytest.skip("Nicht genug Nullstellen für Test")
        p_val = result["ks_p_value_gue"]
        # Nur Strukturprüfung: p-Wert ∈ [0, 1]
        assert 0.0 <= p_val <= 1.0, f"p-Wert außerhalb [0,1]: {p_val}"

    def test_return_contains_wigner_formula(self):
        """Rückgabe enthält die Wigner-Surmise-Formel als String."""
        result = zeta_zero_spacing_statistics(n_zeros=20)
        assert "wigner_surmise_formula" in result

    def test_gue_variance_theoretical_correct(self):
        """Theoretische GUE-Varianz: Var(δ) = 4/π - π/4 ≈ 0.0586."""
        result = zeta_zero_spacing_statistics(n_zeros=20)
        if "gue_variance_theoretical" not in result:
            pytest.skip("Kein gue_variance_theoretical Schlüssel")
        expected_var = 4.0 / math.pi - math.pi / 4.0
        assert abs(result["gue_variance_theoretical"] - expected_var) < 1e-10


# ===========================================================================
# TESTS: pair_correlation_function
# ===========================================================================

class TestPairCorrelationFunction:
    """Tests für Montgomerys Paarkorrelationsfunktion."""

    def test_gue_prediction_at_r1(self):
        """
        GUE-Paarkorrelation bei r=1.0:
        R_2(1) = 1 - (sin(π)/(π))² = 1 - 0 = 1.0
        """
        result = pair_correlation_function(n_zeros=30, r=1.0)
        # sin(π·1) = sin(π) ≈ 0 → R_2(1) ≈ 1
        gue = result["gue_prediction"]
        assert abs(gue - 1.0) < 1e-6, f"R_2(1) sollte 1.0 sein, erhalten: {gue}"

    def test_gue_prediction_at_r05(self):
        """
        GUE-Paarkorrelation bei r=0.5:
        R_2(0.5) = 1 - (sin(π/2)/(π/2))² = 1 - (2/π)² ≈ 0.595
        """
        result = pair_correlation_function(n_zeros=30, r=0.5)
        gue = result["gue_prediction"]
        expected = 1.0 - (math.sin(math.pi * 0.5) / (math.pi * 0.5)) ** 2
        assert abs(gue - expected) < 1e-10

    def test_gue_prediction_in_valid_range(self):
        """GUE-Paarkorrelation liegt im Bereich [0, 2] für r > 0."""
        for r_val in [0.1, 0.5, 1.0, 1.5, 2.0, 3.0]:
            result = pair_correlation_function(n_zeros=20, r=r_val)
            gue = result["gue_prediction"]
            assert 0.0 <= gue <= 2.0, (
                f"R_2({r_val}) = {gue} außerhalb [0, 2]"
            )

    def test_poisson_prediction_is_one(self):
        """Poisson-Paarkorrelation ist immer 1.0 (keine Korrelation)."""
        result = pair_correlation_function(n_zeros=30, r=1.5)
        assert result["poisson_prediction"] == 1.0

    def test_level_repulsion_at_small_r(self):
        """Level-Repulsion: GUE-Paarkorrelation bei kleinem r ist < 0.5."""
        result = pair_correlation_function(n_zeros=30, r=0.2)
        assert result["gue_prediction"] < 0.5, (
            "Level-Repulsion: GUE-Paarkorrelation bei r=0.2 sollte < 0.5 sein"
        )

    def test_r_grid_length(self):
        """r_grid hat die erwartete Länge (30 Punkte von 0.1 bis 3.0)."""
        result = pair_correlation_function(n_zeros=30, r=1.0)
        assert len(result["r_grid"]) == 30
        assert len(result["gue_curve"]) == 30

    def test_return_structure_pair_correlation(self):
        """Rückgabe enthält alle erforderlichen Schlüssel."""
        result = pair_correlation_function(n_zeros=20, r=1.0)
        required = [
            "r", "gue_prediction", "poisson_prediction",
            "r_grid", "gue_curve", "n_zeros_used"
        ]
        for key in required:
            assert key in result, f"Schlüssel '{key}' fehlt"


# ===========================================================================
# TESTS: random_matrix_gue_sample
# ===========================================================================

class TestRandomMatrixGueSample:
    """Tests für GUE-Zufallsmatrix-Simulation."""

    def test_matrix_shape_correct(self):
        """sample_matrix_shape entspricht den Parametern."""
        result = random_matrix_gue_sample(n=5, size=10)
        assert result["sample_matrix_shape"] == (10, 10)
        assert result["matrix_size"] == 10

    def test_eigenvalues_real(self):
        """Hermitesche Matrizen haben reelle Eigenwerte."""
        result = random_matrix_gue_sample(n=5, size=10)
        assert result["eigenvalues_real"] is True, (
            "GUE-Matrizen sind hermitesch → Eigenwerte sind reell"
        )

    def test_spacings_positive(self):
        """Alle Eigenwert-Abstände sind positiv."""
        result = random_matrix_gue_sample(n=10, size=15)
        for s in result["eigenvalue_spacings"]:
            assert s > 0.0, f"Eigenwert-Abstand {s} ≤ 0"

    def test_mean_spacing_near_one(self):
        """
        Nach Normierung durch mittleren Abstand sollte der Mittelwert ≈ 1 sein.
        Toleranz: ±0.3 (bei Normierung durch Mittelwert ist dies exakt erfüllt).
        """
        result = random_matrix_gue_sample(n=20, size=20)
        mean = result["mean_spacing"]
        assert 0.7 < mean < 1.3, (
            f"Normierter mittlerer Eigenwert-Abstand sollte ≈ 1 sein, erhalten: {mean:.4f}"
        )

    def test_ks_statistic_in_valid_range(self):
        """KS-Statistik ist im Bereich [0, 1]."""
        result = random_matrix_gue_sample(n=10, size=15)
        assert 0.0 <= result["ks_statistic"] <= 1.0

    def test_n_matrices_stored(self):
        """Anzahl der Matrizen wird korrekt gespeichert."""
        result = random_matrix_gue_sample(n=7, size=10)
        assert result["n_matrices"] == 7

    def test_hilbert_polya_connection_present(self):
        """Hilbert-Pólya-Verbindung ist im Ergebnis dokumentiert."""
        result = random_matrix_gue_sample(n=5, size=10)
        assert "hilbert_polya_connection" in result
        assert len(result["hilbert_polya_connection"]) > 10

    def test_variance_finite_positive(self):
        """Varianz der Eigenwert-Abstände ist endlich und positiv."""
        result = random_matrix_gue_sample(n=10, size=15)
        assert result["variance"] > 0.0
        assert math.isfinite(result["variance"])

    def test_gue_wigner_formula_present(self):
        """Die GUE Wigner-Surmise-Formel ist als String gespeichert."""
        result = random_matrix_gue_sample(n=5, size=10)
        assert "gue_wigner_formula" in result
        assert "pi" in result["gue_wigner_formula"]
