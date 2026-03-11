"""
@file test_spinor_calculus.py
@brief Tests für das spinor_calculus-Modul.
@description
    Testet alle Klassen aus spinor_calculus.py:
    - PauliMatrices: sigma(), commutator(), anticommutator(), spin_rotation()
    - CliffordAlgebra: generators(), clifford_product(), grade_involution(),
                       reversion(), verify_clifford_relations()
    - Spinor: norm(), conjugate(), lorentz_transform()
    - DiracEquation: gamma_matrices(), free_particle_energy(), plane_wave_solution(),
                     verify_clifford_algebra()

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import pytest

from spinor_calculus import (
    PauliMatrices,
    CliffordAlgebra,
    Spinor,
    DiracEquation,
)


# ===========================================================================
# TESTS: PauliMatrices
# ===========================================================================

class TestPauliMatrices:
    """Tests für die Pauli-Matrizen."""

    def test_sigma_shapes(self):
        """Alle drei Pauli-Matrizen sind 2×2."""
        for i in [1, 2, 3]:
            sigma = PauliMatrices.sigma(i)
            assert sigma.shape == (2, 2), f"σ_{i} muss 2×2 sein"

    def test_sigma_hermitian(self):
        """Alle Pauli-Matrizen sind hermitesch: σᵢ† = σᵢ."""
        for i in [1, 2, 3]:
            sigma = PauliMatrices.sigma(i)
            assert np.allclose(sigma, sigma.conj().T, atol=1e-12), (
                f"σ_{i} muss hermitesch sein"
            )

    def test_sigma_traceless(self):
        """Alle Pauli-Matrizen sind spurlos: Tr(σᵢ) = 0."""
        for i in [1, 2, 3]:
            sigma = PauliMatrices.sigma(i)
            assert np.isclose(np.trace(sigma), 0.0, atol=1e-12), (
                f"Tr(σ_{i}) muss 0 sein"
            )

    def test_sigma_squared_is_identity(self):
        """σᵢ² = I für alle i."""
        I2 = np.eye(2, dtype=complex)
        for i in [1, 2, 3]:
            sigma = PauliMatrices.sigma(i)
            assert np.allclose(sigma @ sigma, I2, atol=1e-12), (
                f"σ_{i}² muss die Einheitsmatrix sein"
            )

    def test_sigma1_values(self):
        """σ₁ = [[0,1],[1,0]]."""
        sigma1 = PauliMatrices.sigma(1)
        expected = np.array([[0, 1], [1, 0]], dtype=complex)
        assert np.allclose(sigma1, expected), f"σ₁ falsch: {sigma1}"

    def test_sigma2_values(self):
        """σ₂ = [[0,-i],[i,0]]."""
        sigma2 = PauliMatrices.sigma(2)
        expected = np.array([[0, -1j], [1j, 0]], dtype=complex)
        assert np.allclose(sigma2, expected), f"σ₂ falsch: {sigma2}"

    def test_sigma3_values(self):
        """σ₃ = [[1,0],[0,-1]]."""
        sigma3 = PauliMatrices.sigma(3)
        expected = np.array([[1, 0], [0, -1]], dtype=complex)
        assert np.allclose(sigma3, expected), f"σ₃ falsch: {sigma3}"

    def test_commutator_12(self):
        """[σ₁, σ₂] = 2i·σ₃."""
        comm = PauliMatrices.commutator(1, 2)
        expected = 2j * PauliMatrices.sigma(3)
        assert np.allclose(comm, expected, atol=1e-12), (
            f"[σ₁, σ₂] muss 2i·σ₃ sein, erhalten: {comm}"
        )

    def test_commutator_23(self):
        """[σ₂, σ₃] = 2i·σ₁."""
        comm = PauliMatrices.commutator(2, 3)
        expected = 2j * PauliMatrices.sigma(1)
        assert np.allclose(comm, expected, atol=1e-12), (
            f"[σ₂, σ₃] muss 2i·σ₁ sein"
        )

    def test_commutator_31(self):
        """[σ₃, σ₁] = 2i·σ₂."""
        comm = PauliMatrices.commutator(3, 1)
        expected = 2j * PauliMatrices.sigma(2)
        assert np.allclose(comm, expected, atol=1e-12), (
            f"[σ₃, σ₁] muss 2i·σ₂ sein"
        )

    def test_anticommutator_same_index(self):
        """{{σᵢ, σᵢ}} = 2·I₂."""
        I2 = np.eye(2, dtype=complex)
        for i in [1, 2, 3]:
            anticomm = PauliMatrices.anticommutator(i, i)
            assert np.allclose(anticomm, 2 * I2, atol=1e-12), (
                f"{{σ_{i}, σ_{i}}} muss 2I sein"
            )

    def test_anticommutator_different_indices(self):
        """{{σᵢ, σⱼ}} = 0 für i ≠ j."""
        zero2 = np.zeros((2, 2), dtype=complex)
        for i, j in [(1, 2), (2, 3), (1, 3)]:
            anticomm = PauliMatrices.anticommutator(i, j)
            assert np.allclose(anticomm, zero2, atol=1e-12), (
                f"{{σ_{i}, σ_{j}}} muss 0 sein für i≠j"
            )

    def test_spin_rotation_identity(self):
        """Rotation um 0° gibt Einheitsmatrix."""
        U = PauliMatrices.spin_rotation(np.array([0.0, 0.0, 1.0]), 0.0)
        assert np.allclose(U, np.eye(2, dtype=complex), atol=1e-12), (
            "Rotation um 0° muss Einheitsmatrix geben"
        )

    def test_spin_rotation_360_degrees(self):
        """Rotation um 360° gibt -I (Spinor-Vorzeichen!)."""
        U = PauliMatrices.spin_rotation(np.array([0.0, 0.0, 1.0]), 2 * math.pi)
        assert np.allclose(U, -np.eye(2, dtype=complex), atol=1e-12), (
            "360°-Rotation muss -I geben (Spinor-Eigenschaft)"
        )

    def test_spin_rotation_unitary(self):
        """Spin-Rotationsmatrix ist unitär: U†U = I."""
        axis = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)
        angle = math.pi / 3
        U = PauliMatrices.spin_rotation(axis, angle)
        I2 = np.eye(2, dtype=complex)
        assert np.allclose(U.conj().T @ U, I2, atol=1e-12), (
            "Spin-Rotation muss unitär sein"
        )

    def test_invalid_sigma_index(self):
        """Ungültiger Index wirft ValueError."""
        with pytest.raises(ValueError):
            PauliMatrices.sigma(0)
        with pytest.raises(ValueError):
            PauliMatrices.sigma(4)


# ===========================================================================
# TESTS: CliffordAlgebra
# ===========================================================================

class TestCliffordAlgebra:
    """Tests für die Clifford-Algebra."""

    def test_generators_count(self):
        """Anzahl der Generatoren ist p + q."""
        for p, q in [(1, 0), (2, 0), (3, 0), (1, 3)]:
            cl = CliffordAlgebra(p, q)
            gammas = cl.generators()
            assert len(gammas) == p + q, (
                f"Cl({p},{q}): Anzahl Generatoren muss {p+q} sein, "
                f"erhalten: {len(gammas)}"
            )

    def test_generators_are_matrices(self):
        """Alle Generatoren sind quadratische Matrizen."""
        for p, q in [(1, 0), (2, 0), (3, 0)]:
            cl = CliffordAlgebra(p, q)
            gammas = cl.generators()
            for i, gamma in enumerate(gammas):
                assert gamma.ndim == 2, f"γ_{i} muss 2D sein"
                assert gamma.shape[0] == gamma.shape[1], f"γ_{i} muss quadratisch sein"

    def test_dirac_algebra_clifford_relations(self):
        """Cl(1,3): Clifford-Relationen {γ^μ, γ^ν} = 2g^{μν}·I erfüllt."""
        cl = CliffordAlgebra(1, 3)
        assert cl.verify_clifford_relations(), (
            "Dirac-Algebra Cl(1,3) muss Clifford-Relationen erfüllen"
        )

    def test_pauli_algebra_clifford_relations(self):
        """Cl(3,0): Clifford-Relationen erfüllt."""
        cl = CliffordAlgebra(3, 0)
        assert cl.verify_clifford_relations(), (
            "Pauli-Algebra Cl(3,0) muss Clifford-Relationen erfüllen"
        )

    def test_clifford_product_is_matrix_product(self):
        """Clifford-Produkt entspricht Matrixprodukt."""
        cl = CliffordAlgebra(3, 0)
        gammas = cl.generators()
        a = gammas[0]
        b = gammas[1]
        product = cl.clifford_product(a, b)
        assert np.allclose(product, a @ b, atol=1e-12), (
            "Clifford-Produkt muss Matrixprodukt sein"
        )

    def test_reversion_involution(self):
        """Reversion ist involutiv: ~~a = a."""
        cl = CliffordAlgebra(3, 0)
        gammas = cl.generators()
        a = gammas[0]
        rev_rev_a = cl.reversion(cl.reversion(a))
        assert np.allclose(rev_rev_a, a, atol=1e-10), (
            "Zweifache Reversion muss identisch sein"
        )


# ===========================================================================
# TESTS: Spinor
# ===========================================================================

class TestSpinor:
    """Tests für die Spinor-Klasse."""

    def test_norm_calculation(self):
        """Spinor-Norm wird korrekt berechnet."""
        psi = Spinor([1.0, 0.0, 0.0, 0.0])
        assert np.isclose(psi.norm(), 1.0), "Norm von (1,0,0,0) muss 1 sein"

    def test_norm_general(self):
        """Spinor-Norm für allgemeinen Spinor."""
        psi = Spinor([1.0, 1.0, 0.0, 0.0])
        assert np.isclose(psi.norm(), math.sqrt(2)), (
            "Norm von (1,1,0,0) muss √2 sein"
        )

    def test_conjugate_returns_array(self):
        """conjugate() gibt numpy-Array zurück."""
        psi = Spinor([1.0, 0.0, 0.0, 0.0])
        psi_bar = psi.conjugate()
        assert isinstance(psi_bar, np.ndarray), "Konjugierte muss numpy-Array sein"

    def test_lorentz_transform_identity(self):
        """Transformation mit Einheitsmatrix lässt Spinor unverändert."""
        components = np.array([1.0, 0.0, 0.5, -0.3], dtype=complex)
        psi = Spinor(components)
        I4 = np.eye(4, dtype=complex)
        psi_transformed = psi.lorentz_transform(I4)
        assert np.allclose(psi_transformed.components, components, atol=1e-12), (
            "Transformation mit I muss Spinor unverändert lassen"
        )


# ===========================================================================
# TESTS: DiracEquation
# ===========================================================================

class TestDiracEquation:
    """Tests für die Dirac-Gleichung."""

    def test_gamma_matrices_count(self):
        """gamma_matrices() gibt 4 Matrizen zurück."""
        dirac = DiracEquation(mass=1.0)
        gammas = dirac.gamma_matrices()
        assert len(gammas) == 4, "Dirac-Gleichung muss 4 Gamma-Matrizen haben"

    def test_gamma_matrices_shape(self):
        """Alle Gamma-Matrizen sind 4×4."""
        dirac = DiracEquation(mass=1.0)
        for i, gamma in enumerate(dirac.gamma_matrices()):
            assert gamma.shape == (4, 4), f"γ^{i} muss 4×4 sein"

    def test_verify_clifford_algebra(self):
        """Dirac-Gamma-Matrizen erfüllen Clifford-Algebra-Relationen."""
        dirac = DiracEquation(mass=1.0)
        assert dirac.verify_clifford_algebra(), (
            "Gamma-Matrizen müssen Clifford-Algebra-Relationen {γ^μ,γ^ν} = 2g^{μν} erfüllen"
        )

    def test_gamma0_hermitian(self):
        """γ⁰ ist hermitesch: (γ⁰)† = γ⁰."""
        dirac = DiracEquation()
        gamma0 = dirac.gamma_matrices()[0]
        assert np.allclose(gamma0, gamma0.conj().T, atol=1e-12), (
            "γ⁰ muss hermitesch sein"
        )

    def test_gamma_spatial_antihermitian(self):
        """γⁱ (i=1,2,3) sind antihermitesch: (γⁱ)† = -γⁱ."""
        dirac = DiracEquation()
        gammas = dirac.gamma_matrices()
        for i in [1, 2, 3]:
            gamma_i = gammas[i]
            assert np.allclose(gamma_i.conj().T, -gamma_i, atol=1e-12), (
                f"γ^{i} muss antihermitesch sein"
            )

    def test_free_particle_energy_massless(self):
        """Masseloses Teilchen: E = |p| (Lichtgeschwindigkeit)."""
        dirac = DiracEquation(mass=0.0)
        p = np.array([3.0, 4.0, 0.0])
        E = dirac.free_particle_energy(p)
        assert np.isclose(E, 5.0, atol=1e-12), (
            f"E = |p| = 5.0 erwartet, erhalten: {E}"
        )

    def test_free_particle_energy_rest(self):
        """Ruhendes Teilchen: E = m (Ruhemasse)."""
        mass = 2.5
        dirac = DiracEquation(mass=mass)
        E = dirac.free_particle_energy(np.array([0.0, 0.0, 0.0]))
        assert np.isclose(E, mass, atol=1e-12), (
            f"E = m = {mass} erwartet, erhalten: {E}"
        )

    def test_free_particle_energy_relativistic(self):
        """E² = p² + m² (Einstein-Dispersionsrelation)."""
        mass = 1.0
        dirac = DiracEquation(mass=mass)
        p = np.array([1.0, 2.0, 2.0])
        E = dirac.free_particle_energy(p)
        p_squared = np.sum(p ** 2)
        assert np.isclose(E ** 2, p_squared + mass ** 2, atol=1e-12), (
            "E² = p² + m² muss erfüllt sein"
        )

    def test_plane_wave_solution_shape(self):
        """Ebene-Welle-Lösung ist 4-komponentiger Spinor."""
        dirac = DiracEquation(mass=1.0)
        p = np.array([0.0, 0.0, 1.0])
        u = dirac.plane_wave_solution(p, spin='up')
        assert u.shape == (4,), f"Spinor muss 4 Komponenten haben, erhalten: {u.shape}"

    def test_plane_wave_both_spins(self):
        """Ebene-Welle-Lösung funktioniert für 'up' und 'down'."""
        dirac = DiracEquation(mass=1.0)
        p = np.array([1.0, 0.0, 0.0])
        u_up = dirac.plane_wave_solution(p, spin='up')
        u_down = dirac.plane_wave_solution(p, spin='down')
        assert u_up.shape == (4,), "spin-up Spinor muss 4 Komponenten haben"
        assert u_down.shape == (4,), "spin-down Spinor muss 4 Komponenten haben"
        # Die beiden Lösungen müssen verschieden sein
        assert not np.allclose(u_up, u_down), (
            "spin-up und spin-down Lösungen müssen verschieden sein"
        )
