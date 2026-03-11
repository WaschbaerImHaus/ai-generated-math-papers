"""
@file test_symplectic_geometry.py
@brief Tests für das symplectic_geometry-Modul.
@description
    Testet alle Klassen und Funktionen aus symplectic_geometry.py:
    - SymplecticForm: evaluate, is_non_degenerate, is_closed, is_antisymmetric
    - HamiltonianSystem: vector_field, solve, poisson_bracket
    - SymplecticManifold: is_symplectic, darboux_theorem_local
    - harmonic_oscillator_hamiltonian
    - kepler_hamiltonian
    - double_pendulum_hamiltonian

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import os
import sys
import math

# src-Verzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pytest

from symplectic_geometry import (
    SymplecticForm,
    HamiltonianSystem,
    SymplecticManifold,
    harmonic_oscillator_hamiltonian,
    kepler_hamiltonian,
    double_pendulum_hamiltonian,
)


# ===========================================================================
# TESTS: SymplecticForm
# ===========================================================================

class TestSymplecticForm:
    """Tests für die SymplecticForm-Klasse."""

    def test_standard_form_dim1(self):
        """Standard-Symplektikform für n=1 ist [[0,1],[-1,0]]."""
        omega = SymplecticForm(1)
        expected = np.array([[0.0, 1.0], [-1.0, 0.0]])
        assert np.allclose(omega.matrix, expected), (
            "Standard-Symplektikmatrix für n=1 falsch"
        )

    def test_standard_form_dim2(self):
        """Standard-Symplektikform für n=2 hat korrekte Blockstruktur."""
        omega = SymplecticForm(2)
        J = omega.matrix
        assert J.shape == (4, 4), "Matrix muss 4×4 sein"
        # Überprüfe Blockstruktur: [[0,I],[-I,0]]
        assert np.allclose(J[:2, :2], np.zeros((2, 2))), "Oberer linker Block muss 0 sein"
        assert np.allclose(J[:2, 2:], np.eye(2)), "Oberer rechter Block muss I sein"
        assert np.allclose(J[2:, :2], -np.eye(2)), "Unterer linker Block muss -I sein"

    def test_evaluate_antisymmetry(self):
        """ω(u,v) = -ω(v,u) (Antisymmetrie)."""
        omega = SymplecticForm(2)
        u = np.array([1.0, 0.5, -1.0, 2.0])
        v = np.array([0.3, -1.0, 0.7, 1.2])
        assert np.isclose(omega.evaluate(u, v), -omega.evaluate(v, u)), (
            "Symplektische Form muss antisymmetrisch sein"
        )

    def test_evaluate_standard_case(self):
        """ω(e₁, e₂) = 0, ω(q₁, p₁) = 1 für Standardbasis."""
        omega = SymplecticForm(2)
        # Standardbasis-Vektoren im ℝ⁴: q₁=(1,0,0,0), p₁=(0,0,1,0)
        q1 = np.array([1.0, 0.0, 0.0, 0.0])
        p1 = np.array([0.0, 0.0, 1.0, 0.0])
        # ω(q₁, p₁) = 1
        assert np.isclose(omega.evaluate(q1, p1), 1.0), (
            "ω(q₁, p₁) muss 1 sein (kanonische Relation)"
        )

    def test_is_non_degenerate(self):
        """Standardsymplektikform ist nicht-entartet."""
        for n in [1, 2, 3]:
            omega = SymplecticForm(n)
            assert omega.is_non_degenerate(), f"n={n}: Standardform muss nicht-entartet sein"

    def test_is_closed(self):
        """Standardsymplektikform ist geschlossen."""
        omega = SymplecticForm(2)
        assert omega.is_closed(), "Standardsymplektikform muss geschlossen sein"

    def test_is_antisymmetric(self):
        """Symplektikmatrix J ist antisymmetrisch: Jᵀ = -J."""
        for n in [1, 2, 3]:
            omega = SymplecticForm(n)
            assert omega.is_antisymmetric(), (
                f"n={n}: Symplektikmatrix muss antisymmetrisch sein"
            )

    def test_wrong_vector_length_raises(self):
        """evaluate() mit falscher Vektorlänge wirft ValueError."""
        omega = SymplecticForm(2)
        u = np.array([1.0, 2.0])  # Zu kurz (braucht Länge 4)
        v = np.array([3.0, 4.0])
        with pytest.raises(ValueError):
            omega.evaluate(u, v)

    def test_custom_omega_matrix(self):
        """Benutzerdefinierte Symplektikmatrix wird akzeptiert."""
        J = np.array([[0.0, 2.0], [-2.0, 0.0]])
        omega = SymplecticForm(1, omega_matrix=J)
        u = np.array([1.0, 0.0])
        v = np.array([0.0, 1.0])
        assert np.isclose(omega.evaluate(u, v), 2.0), "ω(e₁, e₂) muss 2 sein"


# ===========================================================================
# TESTS: HamiltonianSystem
# ===========================================================================

class TestHamiltonianSystem:
    """Tests für das HamiltonianSystem."""

    def test_harmonic_oscillator_vector_field(self):
        """Harmonischer Oszillator: dq/dt = p, dp/dt = -ω²q."""
        omega = 1.0
        H = lambda q, p: harmonic_oscillator_hamiltonian(q, p, omega)
        sys = HamiltonianSystem(H, dim=1)

        q = np.array([1.0])
        p = np.array([0.0])
        dq, dp = sys.vector_field(q, p)

        # dq/dt = ∂H/∂p = p = 0
        assert np.isclose(dq[0], 0.0, atol=1e-4), f"dq/dt falsch: {dq[0]}"
        # dp/dt = -∂H/∂q = -ω²q = -1
        assert np.isclose(dp[0], -1.0, atol=1e-4), f"dp/dt falsch: {dp[0]}"

    def test_harmonic_oscillator_solve_returns_correct_shape(self):
        """solve() gibt korrekte Shapes für Trajektorie zurück."""
        H = lambda q, p: harmonic_oscillator_hamiltonian(q, p, 1.0)
        sys = HamiltonianSystem(H, dim=1)

        t, q_traj, p_traj = sys.solve(
            q0=[1.0], p0=[0.0],
            t_span=(0.0, 2 * math.pi),
            num_points=100
        )

        assert len(t) == 100, "t muss 100 Punkte haben"
        assert q_traj.shape == (1, 100), f"q_traj-Shape falsch: {q_traj.shape}"
        assert p_traj.shape == (1, 100), f"p_traj-Shape falsch: {p_traj.shape}"

    def test_harmonic_oscillator_energy_conserved(self):
        """Energie bleibt beim harmonischen Oszillator erhalten."""
        omega = 1.0
        H = lambda q, p: harmonic_oscillator_hamiltonian(q, p, omega)
        sys = HamiltonianSystem(H, dim=1)

        t, q_traj, p_traj = sys.solve(
            q0=[1.0], p0=[0.0],
            t_span=(0.0, 4.0),
            num_points=200
        )

        # Energie an jedem Punkt berechnen
        energies = [H(q_traj[:, i], p_traj[:, i]) for i in range(len(t))]
        energies = np.array(energies)
        # Energie muss näherungsweise konstant sein
        assert np.std(energies) < 1e-3, (
            f"Energie nicht erhalten (std = {np.std(energies):.6f})"
        )

    def test_poisson_bracket_canonical_relations(self):
        """Kanonische Poisson-Relationen: {qᵢ, pⱼ} = δᵢⱼ."""
        H = lambda q, p: harmonic_oscillator_hamiltonian(q, p, 1.0)
        sys = HamiltonianSystem(H, dim=1)

        # f(q,p) = q[0], g(q,p) = p[0]
        f_q = lambda q, p: float(q[0])
        g_p = lambda q, p: float(p[0])

        q = np.array([0.5])
        p = np.array([1.0])
        bracket = sys.poisson_bracket(f_q, g_p, q, p)

        assert np.isclose(bracket, 1.0, atol=1e-4), (
            f"{{q, p}} muss 1 sein, erhalten: {bracket}"
        )

    def test_poisson_bracket_antisymmetry(self):
        """Poisson-Klammer ist antisymmetrisch: {f,g} = -{g,f}."""
        H = lambda q, p: harmonic_oscillator_hamiltonian(q, p, 1.0)
        sys = HamiltonianSystem(H, dim=1)

        f = lambda q, p: float(q[0] ** 2)
        g = lambda q, p: float(p[0] ** 2)

        q = np.array([1.0])
        p = np.array([0.5])

        fg = sys.poisson_bracket(f, g, q, p)
        gf = sys.poisson_bracket(g, f, q, p)

        assert np.isclose(fg, -gf, atol=1e-4), (
            f"{{f,g}} = -{gf:.4f} aber {{g,f}} = {gf:.4f}: Antisymmetrie verletzt"
        )

    def test_poisson_bracket_energy_self(self):
        """{{H, H}} = 0 (Energie ist Erhaltungsgröße)."""
        omega = 1.5
        H_func = lambda q, p: harmonic_oscillator_hamiltonian(q, p, omega)
        sys = HamiltonianSystem(H_func, dim=1)

        q = np.array([1.0])
        p = np.array([0.5])
        bracket = sys.poisson_bracket(H_func, H_func, q, p)

        assert np.isclose(bracket, 0.0, atol=1e-4), (
            f"{{H, H}} muss 0 sein (Energieerhaltung), erhalten: {bracket}"
        )


# ===========================================================================
# TESTS: SymplecticManifold
# ===========================================================================

class TestSymplecticManifold:
    """Tests für die SymplecticManifold-Klasse."""

    def test_is_symplectic(self):
        """Standard-Mannigfaltigkeit ist symplektisch."""
        M = SymplecticManifold(2)
        assert M.is_symplectic(), "Standard-Mannigfaltigkeit muss symplektisch sein"

    def test_darboux_theorem_returns_dict(self):
        """darboux_theorem_local() gibt Dictionary zurück."""
        M = SymplecticManifold(2)
        result = M.darboux_theorem_local()
        assert isinstance(result, dict), "Rückgabe muss ein Dict sein"
        assert 'q_coords' in result
        assert 'p_coords' in result
        assert 'standard_form' in result

    def test_darboux_theorem_coordinate_count(self):
        """Darboux-Koordinaten haben die richtige Anzahl."""
        for n in [1, 2, 3]:
            M = SymplecticManifold(n)
            result = M.darboux_theorem_local()
            assert len(result['q_coords']) == n, f"n={n}: q-Koordinaten falsch"
            assert len(result['p_coords']) == n, f"n={n}: p-Koordinaten falsch"

    def test_volume_form_positive(self):
        """Symplektisches Volumen ist positiv."""
        M = SymplecticManifold(2)
        vol = M.volume_form()
        assert vol > 0, "Symplektisches Volumen muss positiv sein"

    def test_dimension(self):
        """Gesamtdimension ist 2n."""
        for n in [1, 2, 3, 4]:
            M = SymplecticManifold(n)
            assert M.total_dim == 2 * n, f"Dimension muss 2n={2*n} sein"


# ===========================================================================
# TESTS: Hamilton-Funktionen
# ===========================================================================

class TestHamiltonianFunctions:
    """Tests für die vordefinierten Hamilton-Funktionen."""

    def test_harmonic_oscillator_at_equilibrium(self):
        """H(0,0) = 0 für harmonischen Oszillator."""
        H = harmonic_oscillator_hamiltonian(np.array([0.0]), np.array([0.0]))
        assert np.isclose(H, 0.0), "H(0,0) muss 0 sein"

    def test_harmonic_oscillator_symmetry(self):
        """H(q,p) = H(-q,p) = H(q,-p) (Symmetrie)."""
        q = np.array([1.5])
        p = np.array([0.7])
        H_qp = harmonic_oscillator_hamiltonian(q, p)
        H_mq_p = harmonic_oscillator_hamiltonian(-q, p)
        H_q_mp = harmonic_oscillator_hamiltonian(q, -p)
        assert np.isclose(H_qp, H_mq_p), "H muss symmetrisch in q sein"
        assert np.isclose(H_qp, H_q_mp), "H muss symmetrisch in p sein"

    def test_harmonic_oscillator_positive(self):
        """H(q,p) ≥ 0 für harmonischen Oszillator."""
        for q_val in [-2.0, 0.0, 1.5]:
            for p_val in [-1.0, 0.0, 3.0]:
                H = harmonic_oscillator_hamiltonian(
                    np.array([q_val]), np.array([p_val])
                )
                assert H >= 0, f"H(q={q_val}, p={p_val}) = {H} muss ≥ 0 sein"

    def test_kepler_hamiltonian_far_field(self):
        """Kepler-Hamiltonian im Fernfeld: Potential ist klein."""
        q = np.array([100.0, 0.0, 0.0])
        p = np.array([0.0, 0.0, 0.0])
        H = kepler_hamiltonian(q, p, M=1.0, G=1.0)
        # Bei r=100: H ≈ -GM/r = -0.01
        assert abs(H - (-1.0 / 100.0)) < 1e-10, (
            f"Kepler-H im Fernfeld falsch: {H}"
        )

    def test_kepler_singularity(self):
        """kepler_hamiltonian wirft ValueError bei r=0."""
        with pytest.raises(ValueError):
            kepler_hamiltonian(np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))

    def test_double_pendulum_hamiltonian_runs(self):
        """double_pendulum_hamiltonian läuft ohne Fehler."""
        q = np.array([0.1, 0.2])
        p = np.array([0.0, 0.0])
        H = double_pendulum_hamiltonian(q, p)
        assert np.isfinite(H), "H muss endlich sein"
