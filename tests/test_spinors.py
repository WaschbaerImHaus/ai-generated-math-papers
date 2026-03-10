"""
@file test_spinors.py
@brief Umfassende Tests für src/spinors.py (Spinoren, Clifford-Algebren, Dirac-Gleichung).
@description
    Testet alle Funktionen des spinors-Moduls:

    - Gamma-Matrizen: Clifford-Relationen, Hermitizität, Spur
    - γ^5: Quadrat, Anti-Kommutator, Spur
    - Dirac-Spinoren: Dirac-Gleichung, Normierung
    - Spin-Gruppen-Elemente: Unitarität, Determinante, Spinor-Statistik
    - Weyl-Spinoren: Zerlegung, Projektor-Eigenschaft
    - Majorana-Bedingung
    - Gitter-Dirac-Operator: Hermitizität, Massenlücke

    ## Mathematische Grundlagen der Tests
    Die Clifford-Relation {γ^μ, γ^ν} = 2η^{μν}·I wird für alle (μ,ν)-Paare
    mit der Minkowski-Metrik η = diag(+1,-1,-1,-1) überprüft.

@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0.0
"""

import sys
import os
import pytest
import numpy as np

# Füge src-Verzeichnis zum Python-Pfad hinzu
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from spinors import (
    gamma_matrices,
    clifford_algebra_check,
    gamma5,
    dirac_spinor,
    dirac_equation_check,
    spin_group_element,
    weyl_spinors,
    majorana_condition_check,
    dirac_hamiltonian_1d,
    dirac_spectrum_1d,
)


# ===========================================================================
# Hilfsfunktionen für Tests
# ===========================================================================

def _minkowski_metric(n: int = 4) -> np.ndarray:
    """Gibt die Minkowski-Metrik diag(+1,-1,...,-1) zurück."""
    diag = [1.0] + [-1.0] * (n - 1)
    return np.diag(diag)


# ===========================================================================
# Tests: Gamma-Matrizen (dim=4, Dirac-Darstellung)
# ===========================================================================

class TestGammaMatrices4D:
    """Tests für 4×4-Gamma-Matrizen in der Dirac-Darstellung."""

    def test_gamma_matrices_returns_four_matrices(self):
        """
        @brief Gamma-Matrizen-Liste hat genau 4 Einträge.
        """
        gammas = gamma_matrices(4)
        assert len(gammas) == 4, "Es sollten 4 Gamma-Matrizen für dim=4 zurückgegeben werden"

    def test_gamma_matrices_shape_4x4(self):
        """
        @brief Alle Gamma-Matrizen sind 4×4-Matrizen.
        """
        gammas = gamma_matrices(4)
        for mu, g in enumerate(gammas):
            assert g.shape == (4, 4), f"γ^{mu} sollte eine 4×4-Matrix sein, got {g.shape}"

    def test_gamma0_gamma0_anticommutator(self):
        """
        @brief {γ^0, γ^0} = 2·I (positive Metrik-Komponente η^{00} = +1).

        Die Clifford-Relation für die zeitartige Komponente:
            {γ^0, γ^0} = 2·(γ^0)² = 2·η^{00}·I = +2·I
        """
        gammas = gamma_matrices(4)
        g0 = gammas[0]
        anticomm = g0 @ g0 + g0 @ g0   # {γ^0, γ^0} = 2·(γ^0)²
        expected = 2.0 * np.eye(4, dtype=complex)
        assert np.allclose(anticomm, expected, atol=1e-12), \
            f"{{γ^0, γ^0}} sollte 2·I sein, Fehler: {np.max(np.abs(anticomm - expected))}"

    def test_gamma1_gamma1_anticommutator(self):
        """
        @brief {γ^1, γ^1} = -2·I (negative Metrik-Komponente η^{11} = -1).

        Die Clifford-Relation für raumartige Komponenten:
            {γ^k, γ^k} = 2·η^{kk}·I = -2·I  für k=1,2,3
        """
        gammas = gamma_matrices(4)
        g1 = gammas[1]
        anticomm = g1 @ g1 + g1 @ g1   # {γ^1, γ^1} = 2·(γ^1)²
        expected = -2.0 * np.eye(4, dtype=complex)
        assert np.allclose(anticomm, expected, atol=1e-12), \
            f"{{γ^1, γ^1}} sollte -2·I sein, Fehler: {np.max(np.abs(anticomm - expected))}"

    def test_gamma2_gamma2_anticommutator(self):
        """
        @brief {γ^2, γ^2} = -2·I.
        """
        gammas = gamma_matrices(4)
        g2 = gammas[2]
        anticomm = 2.0 * (g2 @ g2)
        expected = -2.0 * np.eye(4, dtype=complex)
        assert np.allclose(anticomm, expected, atol=1e-12), \
            "{γ^2, γ^2} sollte -2·I sein"

    def test_gamma3_gamma3_anticommutator(self):
        """
        @brief {γ^3, γ^3} = -2·I.
        """
        gammas = gamma_matrices(4)
        g3 = gammas[3]
        anticomm = 2.0 * (g3 @ g3)
        expected = -2.0 * np.eye(4, dtype=complex)
        assert np.allclose(anticomm, expected, atol=1e-12), \
            "{γ^3, γ^3} sollte -2·I sein"

    def test_off_diagonal_anticommutators_zero(self):
        """
        @brief {γ^μ, γ^ν} = 0 für alle μ ≠ ν.

        Verschiedene Gamma-Matrizen anti-kommutieren:
            {γ^μ, γ^ν} = γ^μγ^ν + γ^νγ^μ = 0  für μ ≠ ν
        """
        gammas = gamma_matrices(4)
        for mu in range(4):
            for nu in range(4):
                if mu != nu:
                    anticomm = gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu]
                    assert np.allclose(anticomm, 0, atol=1e-12), \
                        f"{{γ^{mu}, γ^{nu}}} sollte 0 sein für μ≠ν"

    def test_gamma0_hermitian(self):
        """
        @brief γ^0 ist hermitesch: (γ^0)† = γ^0.
        """
        gammas = gamma_matrices(4)
        g0 = gammas[0]
        assert np.allclose(g0, g0.conj().T, atol=1e-12), \
            "γ^0 sollte hermitesch sein"

    def test_gamma_k_antihermitian(self):
        """
        @brief γ^k sind anti-hermitesch für k=1,2,3: (γ^k)† = -γ^k.
        """
        gammas = gamma_matrices(4)
        for k in range(1, 4):
            gk = gammas[k]
            assert np.allclose(gk, -gk.conj().T, atol=1e-12), \
                f"γ^{k} sollte anti-hermitesch sein: (γ^{k})† = -γ^{k}"


# ===========================================================================
# Tests: Clifford-Algebra-Überprüfung
# ===========================================================================

class TestCliffordAlgebraCheck:
    """Tests für die Clifford-Algebra-Überprüfungsfunktion."""

    def test_clifford_check_max_error_small(self):
        """
        @brief Maximaler Fehler bei korrekten Gamma-Matrizen < 1e-12.
        """
        gammas = gamma_matrices(4)
        result = clifford_algebra_check(gammas)
        assert result['max_error'] < 1e-12, \
            f"max_error sollte < 1e-12 sein, got {result['max_error']}"

    def test_clifford_check_is_clifford_algebra_true(self):
        """
        @brief is_clifford_algebra = True für korrekte Gamma-Matrizen.
        """
        gammas = gamma_matrices(4)
        result = clifford_algebra_check(gammas)
        assert result['is_clifford_algebra'], \
            "Korrekte Gamma-Matrizen sollten Clifford-Algebra-Bedingung erfüllen"

    def test_clifford_check_returns_anticommutators_dict(self):
        """
        @brief anticommutators enthält alle 16 Paare (μ,ν) mit μ,ν ∈ {0,1,2,3}.
        """
        gammas = gamma_matrices(4)
        result = clifford_algebra_check(gammas)
        assert len(result['anticommutators']) == 16, \
            "16 Paare (0..3)×(0..3) sollten in anticommutators sein"

    def test_clifford_check_diagonal_values(self):
        """
        @brief Diagonalpaare (μ,μ) haben korrekten Anti-Kommutator-Fehler ≈ 0.
        """
        gammas = gamma_matrices(4)
        result = clifford_algebra_check(gammas)
        for mu in range(4):
            assert result['anticommutators'][(mu, mu)] < 1e-12, \
                f"Fehler für ({mu},{mu}) sollte < 1e-12 sein"

    def test_clifford_check_with_wrong_matrices(self):
        """
        @brief Fehlerhafte Matrizen führen zu is_clifford_algebra = False.
        """
        # Zufällige Matrizen erfüllen die Clifford-Relation nicht
        rng = np.random.default_rng(42)
        wrong_gammas = [rng.random((4, 4)) + 1j * rng.random((4, 4)) for _ in range(4)]
        result = clifford_algebra_check(wrong_gammas)
        assert not result['is_clifford_algebra'], \
            "Zufällige Matrizen sollten Clifford-Relation NICHT erfüllen"


# ===========================================================================
# Tests: γ^5 (Chiralitäts-Matrix)
# ===========================================================================

class TestGamma5:
    """Tests für die Chiralitäts-Matrix γ^5."""

    def test_gamma5_square_is_identity(self):
        """
        @brief (γ^5)² = I.

        γ^5 ist ein Involutionsoperator: zweifache Anwendung ergibt Identität.
        """
        g5 = gamma5(4)
        g5_squared = g5 @ g5
        expected = np.eye(4, dtype=complex)
        assert np.allclose(g5_squared, expected, atol=1e-12), \
            f"(γ^5)² sollte I sein, Fehler: {np.max(np.abs(g5_squared - expected))}"

    def test_gamma5_anticommutes_with_all_gammas(self):
        """
        @brief {γ^5, γ^μ} = 0 für alle μ = 0,1,2,3.

        γ^5 anti-kommutiert mit allen 4 Gamma-Matrizen. Diese Eigenschaft
        ist fundamental für die Trennung chiraler Freiheitsgrade.
        """
        g5 = gamma5(4)
        gammas = gamma_matrices(4)
        for mu, gm in enumerate(gammas):
            anticomm = g5 @ gm + gm @ g5
            assert np.allclose(anticomm, 0, atol=1e-12), \
                f"{{γ^5, γ^{mu}}} sollte 0 sein, Fehler: {np.max(np.abs(anticomm))}"

    def test_gamma5_traceless(self):
        """
        @brief Tr(γ^5) = 0.
        """
        g5 = gamma5(4)
        trace = np.trace(g5)
        assert abs(trace) < 1e-12, \
            f"Tr(γ^5) sollte 0 sein, got {trace}"

    def test_gamma5_hermitian(self):
        """
        @brief (γ^5)† = γ^5 (hermitesch).
        """
        g5 = gamma5(4)
        assert np.allclose(g5, g5.conj().T, atol=1e-12), \
            "γ^5 sollte hermitesch sein"

    def test_gamma5_shape(self):
        """
        @brief γ^5 ist eine 4×4-Matrix.
        """
        g5 = gamma5(4)
        assert g5.shape == (4, 4), f"γ^5 sollte 4×4 sein, got {g5.shape}"


# ===========================================================================
# Tests: Dirac-Spinoren
# ===========================================================================

class TestDiracSpinor:
    """Tests für Dirac-Spinoren und die Dirac-Gleichung."""

    def test_dirac_spinor_satisfies_dirac_equation_rest(self):
        """
        @brief Ruhender Spinor (p=0) erfüllt Dirac-Gleichung mit Residuum < 1e-10.

        Für p=0 vereinfacht sich die Dirac-Gleichung zu:
            (m·γ^0 - m·I) u = m·(γ^0 - I)·u = 0
        """
        mass = 1.0
        momentum = np.array([0.0, 0.0, 0.0])
        for spin in [0, 1]:
            spinor = dirac_spinor(mass, momentum, spin)
            result = dirac_equation_check(spinor, mass, momentum)
            assert result['satisfies_dirac'], \
                f"Ruhender Spinor (spin={spin}) sollte Dirac-Gleichung erfüllen, " \
                f"Residuum: {result['residual']}"

    def test_dirac_spinor_satisfies_dirac_equation_moving(self):
        """
        @brief Bewegter Spinor erfüllt Dirac-Gleichung mit Residuum < 1e-10.
        """
        mass = 1.0
        momentum = np.array([1.0, 0.5, 0.3])
        for spin in [0, 1]:
            spinor = dirac_spinor(mass, momentum, spin)
            result = dirac_equation_check(spinor, mass, momentum)
            assert result['satisfies_dirac'], \
                f"Bewegter Spinor (spin={spin}) sollte Dirac-Gleichung erfüllen, " \
                f"Residuum: {result['residual']}"

    def test_dirac_spinor_is_4component(self):
        """
        @brief Dirac-Spinor hat genau 4 Komponenten.
        """
        spinor = dirac_spinor(1.0, np.array([0.0, 0.0, 0.0]), spin=0)
        assert spinor.shape == (4,), f"Spinor sollte 4-komponentig sein, got {spinor.shape}"

    def test_dirac_spinor_residual_below_threshold(self):
        """
        @brief Dirac-Gleichungs-Residuum ist < 1e-10 für verschiedene Massen.
        """
        for mass in [0.5, 1.0, 2.0, 5.0]:
            momentum = np.array([1.0, 0.0, 0.0])
            spinor = dirac_spinor(mass, momentum, spin=0)
            result = dirac_equation_check(spinor, mass, momentum)
            assert result['residual'] < 1e-10, \
                f"Residuum für m={mass} sollte < 1e-10 sein, got {result['residual']}"

    def test_dirac_equation_check_returns_energy(self):
        """
        @brief dirac_equation_check gibt korrekte Energie zurück.

        Für m=1, p=(0,0,0) gilt E = √(1² + 0) = 1.
        """
        mass = 1.0
        momentum = np.array([0.0, 0.0, 0.0])
        spinor = dirac_spinor(mass, momentum, spin=0)
        result = dirac_equation_check(spinor, mass, momentum)
        assert abs(result['energy'] - 1.0) < 1e-10, \
            f"Energie für m=1, p=0 sollte 1.0 sein, got {result['energy']}"

    def test_dirac_spinor_invalid_spin_raises(self):
        """
        @brief Ungültiger Spin-Wert wirft ValueError.
        """
        with pytest.raises(ValueError, match="spin"):
            dirac_spinor(1.0, np.array([0.0, 0.0, 0.0]), spin=2)


# ===========================================================================
# Tests: Spin-Gruppen-Elemente SU(2)
# ===========================================================================

class TestSpinGroupElement:
    """Tests für Spin(3) ≅ SU(2)-Gruppenelemente."""

    def test_spin_group_element_is_unitary(self):
        """
        @brief Spin-Gruppen-Elemente sind unitär: U†·U = I.
        """
        axis = np.array([0.0, 0.0, 1.0])
        U = spin_group_element(np.pi / 3, axis)
        product = U.conj().T @ U
        assert np.allclose(product, np.eye(2), atol=1e-12), \
            "Spin-Gruppen-Element sollte unitär sein (U†U = I)"

    def test_spin_group_element_determinant_one(self):
        """
        @brief |det(U)| = 1 (unitäre Matrix).

        SU(2)-Elemente haben det = +1, allgemein gilt |det| = 1.
        """
        for theta in [0.1, np.pi/4, np.pi/2, np.pi, 2*np.pi]:
            axis = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
            U = spin_group_element(theta, axis)
            det = np.linalg.det(U)
            assert abs(abs(det) - 1.0) < 1e-12, \
                f"|det(U)| sollte 1 sein für θ={theta:.2f}, got {abs(det)}"

    def test_spin_group_2pi_rotation_is_minus_identity(self):
        """
        @brief 2π-Rotation ergibt -I (Spinor-Statistik!).

        Das ist das fundamentale Merkmal von Spinoren:
        U(2π, n̂) = cos(π)·I + i·sin(π)·(σ·n̂) = -I + 0 = -I

        Gewöhnliche Vektoren kehren nach 2π in den Ausgangszustand zurück,
        Spinoren aber erst nach 4π!
        """
        axis = np.array([0.0, 0.0, 1.0])
        U = spin_group_element(2 * np.pi, axis)
        expected = -np.eye(2, dtype=complex)
        assert np.allclose(U, expected, atol=1e-12), \
            f"2π-Rotation sollte -I ergeben (Spinor-Statistik), got\n{U}"

    def test_spin_group_4pi_rotation_is_plus_identity(self):
        """
        @brief 4π-Rotation ergibt +I.

        Erst nach einer vollständigen 4π-Drehung kehrt der Spinor
        in seinen ursprünglichen Zustand zurück.
        """
        axis = np.array([0.0, 0.0, 1.0])
        U = spin_group_element(4 * np.pi, axis)
        expected = np.eye(2, dtype=complex)
        assert np.allclose(U, expected, atol=1e-12), \
            f"4π-Rotation sollte +I ergeben, got\n{U}"

    def test_spin_group_zero_rotation_is_identity(self):
        """
        @brief θ=0-Rotation ergibt I (Neutralelement der Gruppe).
        """
        axis = np.array([1.0, 0.0, 0.0])
        U = spin_group_element(0.0, axis)
        assert np.allclose(U, np.eye(2), atol=1e-12), \
            "θ=0 Rotation sollte Einheitsmatrix ergeben"


# ===========================================================================
# Tests: Weyl-Spinoren (chirale Zerlegung)
# ===========================================================================

class TestWeylSpinors:
    """Tests für die chirale Zerlegung von Dirac-Spinoren."""

    def test_weyl_spinors_sum_is_original(self):
        """
        @brief ψ_L + ψ_R = ψ (vollständige Zerlegung, P_L + P_R = I).
        """
        mass = 1.0
        psi = dirac_spinor(mass, np.array([0.5, 0.3, 0.1]), spin=0)
        result = weyl_spinors(psi)
        reconstructed = result['left'] + result['right']
        assert np.allclose(reconstructed, psi, atol=1e-12), \
            "ψ_L + ψ_R sollte ursprünglichen Spinor ψ ergeben"

    def test_weyl_spinors_left_projector_idempotent(self):
        """
        @brief P_L² = P_L (Projektor-Eigenschaft).

        Projektoren sind idempotent: zweifache Anwendung ändert nichts.
        """
        g5 = gamma5(4)
        I4 = np.eye(4, dtype=complex)
        P_L = (I4 - g5) / 2.0
        # P_L² = P_L
        P_L_sq = P_L @ P_L
        assert np.allclose(P_L_sq, P_L, atol=1e-12), \
            "P_L² sollte P_L sein (Projektor-Eigenschaft)"

    def test_weyl_spinors_right_projector_idempotent(self):
        """
        @brief P_R² = P_R (Projektor-Eigenschaft).
        """
        g5 = gamma5(4)
        I4 = np.eye(4, dtype=complex)
        P_R = (I4 + g5) / 2.0
        P_R_sq = P_R @ P_R
        assert np.allclose(P_R_sq, P_R, atol=1e-12), \
            "P_R² sollte P_R sein (Projektor-Eigenschaft)"

    def test_weyl_spinors_projectors_orthogonal(self):
        """
        @brief P_L · P_R = 0 (linke und rechte Projektion sind orthogonal).
        """
        g5 = gamma5(4)
        I4 = np.eye(4, dtype=complex)
        P_L = (I4 - g5) / 2.0
        P_R = (I4 + g5) / 2.0
        product = P_L @ P_R
        assert np.allclose(product, 0, atol=1e-12), \
            "P_L · P_R sollte 0 sein (Orthogonalität)"

    def test_weyl_spinors_returns_correct_keys(self):
        """
        @brief Rückgabe-Dictionary hat alle erforderlichen Schlüssel.
        """
        psi = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        result = weyl_spinors(psi)
        assert 'left' in result
        assert 'right' in result
        assert 'is_purely_left' in result
        assert 'is_purely_right' in result


# ===========================================================================
# Tests: Majorana-Bedingung
# ===========================================================================

class TestMajoranaCondition:
    """Tests für die Majorana-Fermion-Bedingung."""

    def test_majorana_zero_spinor_is_majorana(self):
        """
        @brief Null-Spinor ψ=0 erfüllt trivialerweise die Majorana-Bedingung.

        Für ψ=0: ψ^c = C·(γ^0·0*)^T = 0 = ψ. Residuum = 0.
        """
        psi_zero = np.zeros(4, dtype=complex)
        result = majorana_condition_check(psi_zero)
        assert result['is_majorana'], \
            "Null-Spinor sollte Majorana-Bedingung erfüllen"
        assert result['residual'] < 1e-12, \
            f"Residuum für Null-Spinor sollte 0 sein, got {result['residual']}"

    def test_majorana_check_returns_correct_keys(self):
        """
        @brief Rückgabe-Dictionary hat 'residual' und 'is_majorana'.
        """
        psi = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
        result = majorana_condition_check(psi)
        assert 'residual' in result
        assert 'is_majorana' in result

    def test_majorana_residual_is_float(self):
        """
        @brief residual ist ein reeller float-Wert.
        """
        psi = np.array([1.0, 2.0, 3.0, 4.0], dtype=complex)
        result = majorana_condition_check(psi)
        assert isinstance(result['residual'], float), \
            "residual sollte ein float sein"

    def test_majorana_residual_nonnegative(self):
        """
        @brief Residuum ist nicht-negativ (Norm ≥ 0).
        """
        rng = np.random.default_rng(123)
        # Erzeuge komplexen Spinor aus zwei reellen Normalverteilungen
        psi = rng.standard_normal(4) + 1j * rng.standard_normal(4)
        result = majorana_condition_check(psi)
        assert result['residual'] >= 0.0, \
            "Majorana-Residuum sollte nicht-negativ sein"


# ===========================================================================
# Tests: Gitter-Dirac-Operator
# ===========================================================================

class TestLatticeDirac:
    """Tests für den Wilson-Fermion-Hamiltonian auf dem 1D-Gitter."""

    def test_dirac_hamiltonian_shape(self):
        """
        @brief Hamilton-Matrix hat korrekte Form (2·n_sites × 2·n_sites).
        """
        n = 8
        H = dirac_hamiltonian_1d(n_sites=n, mass=1.0)
        assert H.shape == (2*n, 2*n), \
            f"Hamiltonian sollte (2·{n} × 2·{n}) = ({2*n}×{2*n}) sein, got {H.shape}"

    def test_dirac_hamiltonian_hermitian(self):
        """
        @brief Wilson-Fermion-Hamiltonian ist hermitesch: H† = H.

        Ein physikalischer Hamiltonian muss hermitesch sein, damit die
        Eigenwerte (Energien) reell sind.
        """
        H = dirac_hamiltonian_1d(n_sites=10, mass=1.0)
        assert np.allclose(H, H.conj().T, atol=1e-12), \
            "Hamiltonian sollte hermitesch sein (H† = H)"

    def test_dirac_spectrum_mass_gap_positive(self):
        """
        @brief Massenlücke > 0 für massive Fermionen.

        Der Wilson-Fermion-Hamiltonian ist gapped: min(|Eigenwerte|) > 0.
        """
        result = dirac_spectrum_1d(mass=1.0, n_sites=20)
        assert result['mass_gap'] > 0, \
            f"Massenlücke sollte > 0 sein, got {result['mass_gap']}"

    def test_dirac_spectrum_is_gapped_true(self):
        """
        @brief is_gapped = True für massive Wilson-Fermionen.
        """
        result = dirac_spectrum_1d(mass=1.0, n_sites=10)
        assert result['is_gapped'], \
            "Spektrum sollte gapped sein für m=1"

    def test_dirac_spectrum_returns_correct_keys(self):
        """
        @brief Rückgabe-Dictionary enthält alle erforderlichen Schlüssel.
        """
        result = dirac_spectrum_1d(mass=1.0, n_sites=10)
        for key in ('eigenvalues', 'mass_gap', 'is_gapped', 'n_zero_modes'):
            assert key in result, f"Schlüssel '{key}' fehlt im Rückgabe-Dictionary"

    def test_dirac_spectrum_eigenvalues_real(self):
        """
        @brief Eigenwerte eines hermiteschen Hamiltonians sind reell.

        numpy gibt für hermitesche Matrizen bei eigvalsh reelle Werte zurück.
        """
        result = dirac_spectrum_1d(mass=0.5, n_sites=12)
        # eigvalsh gibt reelle Werte zurück, Imaginärteil = 0
        assert result['eigenvalues'].dtype in (np.float64, np.float32), \
            f"Eigenwerte sollten reell sein (dtype), got {result['eigenvalues'].dtype}"

    def test_dirac_spectrum_count(self):
        """
        @brief Anzahl der Eigenwerte = 2·n_sites.
        """
        n = 15
        result = dirac_spectrum_1d(mass=1.0, n_sites=n)
        assert len(result['eigenvalues']) == 2 * n, \
            f"Anzahl Eigenwerte sollte 2·{n} = {2*n} sein"


# ===========================================================================
# Tests: 2D Gamma-Matrizen (Edge Case)
# ===========================================================================

class TestGammaMatrices2D:
    """Tests für 2×2-Gamma-Matrizen (2D-Darstellung)."""

    def test_gamma_matrices_2d_returns_two_matrices(self):
        """
        @brief Für dim=2 werden 2 Gamma-Matrizen zurückgegeben.
        """
        gammas = gamma_matrices(2)
        assert len(gammas) == 2, "Für dim=2 sollten 2 Gamma-Matrizen zurückgegeben werden"

    def test_gamma_matrices_2d_shape(self):
        """
        @brief 2D-Gamma-Matrizen sind 2×2.
        """
        gammas = gamma_matrices(2)
        for g in gammas:
            assert g.shape == (2, 2), f"2D-Gamma-Matrizen sollten 2×2 sein, got {g.shape}"

    def test_gamma_matrices_invalid_dim_raises(self):
        """
        @brief Ungültige Dimension wirft ValueError.
        """
        with pytest.raises(ValueError):
            gamma_matrices(3)


# ===========================================================================
# Tests: Gamma5 Fehlerbehandlung
# ===========================================================================

class TestGamma5EdgeCases:
    """Edge-Case-Tests für gamma5."""

    def test_gamma5_invalid_dim_raises(self):
        """
        @brief gamma5 mit dim≠4 wirft ValueError.
        """
        with pytest.raises(ValueError):
            gamma5(2)
