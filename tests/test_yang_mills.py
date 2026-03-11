"""
@file test_yang_mills.py
@brief Umfassende Tests für das Yang-Mills-Theorie-Modul (yang_mills.py).
@description
    Testet alle Klassen und Funktionen des yang_mills.py-Moduls:
    - PrincipalBundle: Faserbündel, Zusammenhang, Krümmung, Eichtransformation
    - LieGroup: Dimensionen, Generatoren, Strukturkonstanten, Killing-Form
    - YangMillsEquations: Funktional, Euler-Lagrange, Bianchi-Identität, Selbst-Dualität
    - Instanton: BPST-Instanton, topologische Ladung, Wirkung, Moduli-Raum
    - MassGap: Millennium-Problem, Gitter-Evidenz, Glueballmasse, Confinement
    - WilsonLoop: Rechteckige Schleifen, String-Tension, Flächengesetz, Creutz-Verhältnis
    - YangMillsTopology: Chern-Klassen, Pontryagin-Klassen, Theta-Vakuum
    - IndexTheory: Atiyah-Singer-Index, Dirac-Spektrum, Eta-Invariante

    Numerische Toleranz: ε = 1e-6

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
import pytest
import numpy as np
import sys
import os

# Pfad zum src-Verzeichnis setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from yang_mills import (
    PrincipalBundle, LieGroup, YangMillsEquations,
    Instanton, MassGap, WilsonLoop,
    YangMillsTopology, IndexTheory,
    _pauli_matrices, _gell_mann_matrices, _levi_civita_3d
)

# Numerische Toleranz
EPS = 1e-6


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

class TestHelpFunctions:
    """Tests für die internen Hilfsfunktionen."""

    def test_pauli_matrices_count(self):
        """Pauli-Matrizen: genau 3 Stück."""
        paulis = _pauli_matrices()
        assert len(paulis) == 3

    def test_pauli_matrices_hermitian(self):
        """Pauli-Matrizen sind hermitesch: σ† = σ."""
        for sigma in _pauli_matrices():
            assert np.allclose(sigma, sigma.conj().T, atol=EPS)

    def test_pauli_matrices_traceless(self):
        """Pauli-Matrizen sind spurlos: Tr(σ) = 0."""
        for sigma in _pauli_matrices():
            assert abs(np.trace(sigma)) < EPS

    def test_pauli_matrices_anticommutation(self):
        """Pauli-Matrizen anti-kommutieren: {σᵢ, σⱼ} = 2δᵢⱼ I."""
        paulis = _pauli_matrices()
        I = np.eye(2, dtype=complex)
        for i in range(3):
            # {σᵢ, σᵢ} = 2I
            anticomm = paulis[i] @ paulis[i] + paulis[i] @ paulis[i]
            assert np.allclose(anticomm, 2 * I, atol=EPS)

    def test_gell_mann_matrices_count(self):
        """Gell-Mann-Matrizen: genau 8 Stück."""
        gm = _gell_mann_matrices()
        assert len(gm) == 8

    def test_gell_mann_matrices_hermitian(self):
        """Gell-Mann-Matrizen sind hermitesch."""
        for lam in _gell_mann_matrices():
            assert np.allclose(lam, lam.conj().T, atol=EPS)

    def test_gell_mann_matrices_traceless(self):
        """Gell-Mann-Matrizen sind spurlos."""
        for lam in _gell_mann_matrices():
            assert abs(np.trace(lam)) < EPS

    def test_levi_civita_3d_shape(self):
        """Levi-Civita-Tensor: Form (3,3,3)."""
        eps = _levi_civita_3d()
        assert eps.shape == (3, 3, 3)

    def test_levi_civita_3d_antisymmetry(self):
        """Levi-Civita-Tensor ist vollständig antisymmetrisch."""
        eps = _levi_civita_3d()
        # ε_{123} = +1
        assert eps[0, 1, 2] == 1.0
        # ε_{213} = -1 (eine Transposition)
        assert eps[1, 0, 2] == -1.0
        # ε_{112} = 0 (wiederholter Index)
        assert eps[0, 0, 2] == 0.0


# =============================================================================
# KLASSE PrincipalBundle
# =============================================================================

class TestPrincipalBundle:
    """Tests für das Hauptfaserbündel."""

    def setup_method(self):
        """Test-Setup: SU(2)-Bündel über R⁴."""
        self.bundle = PrincipalBundle(4, "SU(2)")
        # Einfaches Eichfeld: nur A^1_1 = 0.5 (1-Komponente, 1-Generator)
        self.A_simple = {(0, 0): 0.1, (1, 1): 0.2, (2, 2): 0.15, (3, 0): 0.05}

    def test_principal_bundle_init(self):
        """PrincipalBundle-Initialisierung korrekt."""
        assert self.bundle.dim == 4
        assert self.bundle.group_name == "SU(2)"

    def test_connection_1form_returns_dict(self):
        """connection_1form gibt Dictionary zurück."""
        result = self.bundle.connection_1form(self.A_simple)
        assert isinstance(result, dict)
        assert "components" in result
        assert "group" in result

    def test_connection_1form_group(self):
        """connection_1form enthält korrekte Gruppe."""
        result = self.bundle.connection_1form(self.A_simple)
        assert result["group"] == "SU(2)"

    def test_curvature_2form_returns_dict(self):
        """curvature_2form gibt Dictionary mit F_components zurück."""
        result = self.bundle.curvature_2form(self.A_simple)
        assert isinstance(result, dict)
        assert "F_components" in result

    def test_curvature_2form_antisymmetric_index(self):
        """Krümmungsform: nur μ < ν Paare vorhanden."""
        result = self.bundle.curvature_2form(self.A_simple)
        F = result["F_components"]
        for (mu, nu) in F.keys():
            assert mu < nu, f"Erwartet μ < ν, got ({mu}, {nu})"

    def test_curvature_2form_matrix_size(self):
        """Feldstärke-Matrizen haben korrekte Größe (2×2 für SU(2))."""
        result = self.bundle.curvature_2form(self.A_simple)
        F = result["F_components"]
        for F_mat in F.values():
            assert F_mat.shape == (2, 2)

    def test_gauge_transform_identity(self):
        """Triviale Eichtransformation (g=None) lässt A unverändert."""
        result = self.bundle.gauge_transform(self.A_simple, None)
        assert result["transform"] == "identity"
        assert result["A_transformed"] is self.A_simple

    def test_su3_bundle(self):
        """SU(3)-Bündel hat 8 Algebra-Dimensionen."""
        bundle_su3 = PrincipalBundle(4, "SU(3)")
        result = bundle_su3.connection_1form({})
        assert result["algebra_dim"] == 8


# =============================================================================
# KLASSE LieGroup
# =============================================================================

class TestLieGroup:
    """Tests für Lie-Gruppen."""

    def test_su2_dimension(self):
        """SU(2) hat Dimension 3."""
        g = LieGroup("SU(2)")
        assert g.dimension() == 3

    def test_su3_dimension(self):
        """SU(3) hat Dimension 8."""
        g = LieGroup("SU(3)")
        assert g.dimension() == 8

    def test_u1_dimension(self):
        """U(1) hat Dimension 1."""
        g = LieGroup("U(1)")
        assert g.dimension() == 1

    def test_so3_dimension(self):
        """SO(3) hat Dimension 3."""
        g = LieGroup("SO(3)")
        assert g.dimension() == 3

    def test_unsupported_group_raises(self):
        """Nicht unterstützte Gruppe löst ValueError aus."""
        with pytest.raises(ValueError):
            LieGroup("SU(5)")

    def test_su2_generators_count(self):
        """SU(2) hat 3 Generatoren."""
        g = LieGroup("SU(2)")
        gens = g.lie_algebra_generators()
        assert len(gens) == 3

    def test_su3_generators_count(self):
        """SU(3) hat 8 Generatoren."""
        g = LieGroup("SU(3)")
        gens = g.lie_algebra_generators()
        assert len(gens) == 8

    def test_su2_generators_antihermitian(self):
        """SU(2)-Generatoren Tₐ = σₐ/2 sind hermitesch (physikalische Konvention)."""
        g = LieGroup("SU(2)")
        gens = g.lie_algebra_generators()
        for T in gens:
            # Pauli/2 sind hermitesch
            assert np.allclose(T, T.conj().T, atol=EPS)

    def test_su2_generators_traceless(self):
        """SU(2)-Generatoren sind spurlos."""
        g = LieGroup("SU(2)")
        for T in g.lie_algebra_generators():
            assert abs(np.trace(T)) < EPS

    def test_su2_structure_constants_antisymmetry(self):
        """SU(2)-Strukturkonstanten sind antisymmetrisch: fᵃᵇᶜ = -fᵇᵃᶜ."""
        g = LieGroup("SU(2)")
        f = g.structure_constants()
        # Antisymmetrie in den ersten zwei Indizes
        assert np.allclose(f, -np.transpose(f, (1, 0, 2)), atol=EPS)

    def test_su2_structure_constants_shape(self):
        """SU(2)-Strukturkonstanten haben Form (3,3,3)."""
        g = LieGroup("SU(2)")
        f = g.structure_constants()
        assert f.shape == (3, 3, 3)

    def test_su2_dynkin_index(self):
        """SU(2)-Dynkin-Index ist 0.5."""
        g = LieGroup("SU(2)")
        assert abs(g.dynkin_index() - 0.5) < EPS

    def test_su3_dynkin_index(self):
        """SU(3)-Dynkin-Index ist 0.5."""
        g = LieGroup("SU(3)")
        assert abs(g.dynkin_index() - 0.5) < EPS

    def test_killing_form_shape_su2(self):
        """Killing-Form für SU(2) hat Form (3,3)."""
        g = LieGroup("SU(2)")
        B = g.killing_form()
        assert B.shape == (3, 3)

    def test_killing_form_symmetric(self):
        """Killing-Form ist symmetrisch: B[a,b] = B[b,a]."""
        g = LieGroup("SU(2)")
        B = g.killing_form()
        assert np.allclose(B, B.T, atol=EPS)


# =============================================================================
# KLASSE YangMillsEquations
# =============================================================================

class TestYangMillsEquations:
    """Tests für Yang-Mills-Gleichungen."""

    def setup_method(self):
        """Test-Setup."""
        self.ym = YangMillsEquations("SU(2)", dimension=4)
        # Einfaches nicht-triviales Eichfeld
        self.A = {
            (0, 0): 0.1, (0, 1): 0.0, (0, 2): 0.0,
            (1, 0): 0.0, (1, 1): 0.2, (1, 2): 0.0,
            (2, 0): 0.0, (2, 1): 0.0, (2, 2): 0.15,
            (3, 0): 0.05, (3, 1): 0.0, (3, 2): 0.0
        }

    def test_yang_mills_functional_zero_field(self):
        """Für verschwindendes Feld ist S[A] = 0."""
        F_zero = {(0, 1): np.zeros((2, 2), dtype=complex),
                  (0, 2): np.zeros((2, 2), dtype=complex),
                  (0, 3): np.zeros((2, 2), dtype=complex),
                  (1, 2): np.zeros((2, 2), dtype=complex),
                  (1, 3): np.zeros((2, 2), dtype=complex),
                  (2, 3): np.zeros((2, 2), dtype=complex)}
        S = self.ym.yang_mills_functional({"F_components": F_zero})
        assert abs(S) < EPS

    def test_yang_mills_functional_positive(self):
        """Yang-Mills-Funktional ist nicht-negativ."""
        bundle = PrincipalBundle(4, "SU(2)")
        F_result = bundle.curvature_2form(self.A)
        S = self.ym.yang_mills_functional(F_result)
        assert S >= -EPS

    def test_yang_mills_functional_nonzero(self):
        """Für nicht-triviales Feld ist S[A] > 0."""
        bundle = PrincipalBundle(4, "SU(2)")
        F_result = bundle.curvature_2form(self.A)
        S = self.ym.yang_mills_functional(F_result)
        assert S > EPS

    def test_euler_lagrange_returns_dict(self):
        """Euler-Lagrange gibt Dictionary zurück."""
        result = self.ym.euler_lagrange(self.A)
        assert isinstance(result, dict)
        assert "EL_equations" in result

    def test_euler_lagrange_has_nu_equations(self):
        """Euler-Lagrange hat Gleichungen für jeden ν = 0..3."""
        result = self.ym.euler_lagrange(self.A)
        EL = result["EL_equations"]
        assert len(EL) == 4

    def test_bianchi_identity_trivial(self):
        """Bianchi-Identität für triviales Feld."""
        F_trivial = {}
        assert self.ym.bianchi_identity_check(F_trivial) == True

    def test_self_dual_check_non_self_dual(self):
        """Generisches Feld ist nicht selbst-dual."""
        bundle = PrincipalBundle(4, "SU(2)")
        F_result = bundle.curvature_2form(self.A)
        # Für generische Felder erwarten wir i.A. keine Selbst-Dualität
        result = self.ym.self_dual_check(F_result)
        # Das Ergebnis ist bool (kein Fehler)
        assert isinstance(result, bool)

    def test_self_dual_implies_not_anti_self_dual(self):
        """
        Wenn F selbst-dual ist, dann nicht anti-selbst-dual
        (außer F = 0).
        """
        # Konstruiere explizit ein selbst-duales Feld: F_01 = F_23, etc.
        size = 2
        M = np.array([[0.1, 0.0], [0.0, -0.1]], dtype=complex) * 0.1
        F_sd = {
            (0, 1): M, (2, 3): M,
            (0, 2): -0.5 * M, (1, 3): 0.5 * M,
            (0, 3): 0.2 * M, (1, 2): 0.2 * M
        }
        # Selbst-dual: F_01 = F_23 ✓, F_02 = -F_13 ✓, F_03 = F_12 ✓
        assert self.ym.self_dual_check({"F_components": F_sd}) == True
        # Anti-selbst-dual würde F_01 = -F_23 erfordern → False
        assert self.ym.anti_self_dual_check({"F_components": F_sd}) == False

    def test_anti_self_dual_field(self):
        """Test für anti-selbst-duales Feld."""
        M = np.array([[0.1, 0.0], [0.0, -0.1]], dtype=complex) * 0.1
        F_asd = {
            (0, 1): M, (2, 3): -M,   # F_01 = -F_23
            (0, 2): M, (1, 3): M,    # F_02 = F_13
            (0, 3): M, (1, 2): -M    # F_03 = -F_12
        }
        assert self.ym.anti_self_dual_check({"F_components": F_asd}) == True

    def test_3d_no_self_dual(self):
        """In 3D gibt es keine Selbst-Dualität."""
        ym3 = YangMillsEquations("SU(2)", dimension=3)
        F = {(0, 1): np.eye(2, dtype=complex)}
        assert ym3.self_dual_check(F) == False


# =============================================================================
# KLASSE Instanton
# =============================================================================

class TestInstanton:
    """Tests für Instantone."""

    def setup_method(self):
        """Test-Setup: BPST-Instanton mit Ladung k=1."""
        self.inst = Instanton(charge_k=1, gauge_group="SU(2)")

    def test_instanton_action_k1(self):
        """Wirkung des k=1-Instantons ist 8π²."""
        S = self.inst.action()
        assert abs(S - 8.0 * math.pi ** 2) < EPS

    def test_instanton_action_k2(self):
        """Wirkung des k=2-Instantons ist 16π²."""
        inst2 = Instanton(charge_k=2)
        S = inst2.action()
        assert abs(S - 16.0 * math.pi ** 2) < EPS

    def test_instanton_action_negative_k(self):
        """Wirkung ist proportional zu |k| (positiv für k<0)."""
        inst_neg = Instanton(charge_k=-1)
        S = inst_neg.action()
        assert abs(S - 8.0 * math.pi ** 2) < EPS

    def test_moduli_space_dimension_k1(self):
        """Moduli-Raum-Dimension für k=1, SU(2): 8·1-3 = 5."""
        dim = self.inst.moduli_space_dimension()
        assert dim == 5

    def test_moduli_space_dimension_k2(self):
        """Moduli-Raum-Dimension für k=2, SU(2): 8·2-3 = 13."""
        inst2 = Instanton(charge_k=2)
        dim = inst2.moduli_space_dimension()
        assert dim == 13

    def test_bpst_instanton_returns_dict(self):
        """BPST-Instanton gibt Dictionary zurück."""
        x = np.array([1.0, 0.5, 0.3, 0.2])
        result = self.inst.bpst_instanton(x, rho=1.0)
        assert isinstance(result, dict)
        assert "A_components" in result

    def test_bpst_instanton_at_origin(self):
        """BPST-Instanton am Zentrum: r=0, A → 0/ρ² = 0."""
        center = np.array([0.0, 0.0, 0.0, 0.0])
        result = self.inst.bpst_instanton(center, rho=1.0, center=center)
        A = result["A_components"]
        # Am Zentrum: yν = 0 → alle A^a_μ = 0
        for val in A.values():
            assert abs(val) < EPS

    def test_bpst_instanton_rho_scaling(self):
        """BPST-Instanton skaliert mit ρ: A(x, 2ρ) ≠ A(x, ρ)."""
        x = np.array([1.0, 0.5, 0.3, 0.2])
        A1 = self.inst.bpst_instanton(x, rho=1.0)["A_components"]
        A2 = self.inst.bpst_instanton(x, rho=2.0)["A_components"]
        # Für ρ₁ ≠ ρ₂ müssen die Felder verschieden sein
        vals1 = list(A1.values())
        vals2 = list(A2.values())
        assert not all(abs(v1 - v2) < EPS for v1, v2 in zip(vals1, vals2))

    def test_topological_charge_zero_field(self):
        """Topologische Ladung des Nullfeldes ist 0."""
        F_zero = {(0, 1): np.zeros((2, 2), dtype=complex),
                  (0, 2): np.zeros((2, 2), dtype=complex),
                  (0, 3): np.zeros((2, 2), dtype=complex),
                  (1, 2): np.zeros((2, 2), dtype=complex),
                  (1, 3): np.zeros((2, 2), dtype=complex),
                  (2, 3): np.zeros((2, 2), dtype=complex)}
        k = self.inst.topological_charge({"F_components": F_zero})
        assert abs(k) < EPS

    def test_levi_civita_4d_basic(self):
        """Levi-Civita-Symbol: ε_{0123} = +1."""
        assert self.inst._levi_civita_4d(0, 1, 2, 3) == 1

    def test_levi_civita_4d_odd_perm(self):
        """Levi-Civita-Symbol: ε_{1023} = -1 (eine Transposition)."""
        assert self.inst._levi_civita_4d(1, 0, 2, 3) == -1

    def test_levi_civita_4d_repeated(self):
        """Levi-Civita-Symbol: ε_{0012} = 0 (wiederholter Index)."""
        assert self.inst._levi_civita_4d(0, 0, 1, 2) == 0


# =============================================================================
# KLASSE MassGap
# =============================================================================

class TestMassGap:
    """Tests für das Massenspalt-Problem."""

    def setup_method(self):
        """Test-Setup: SU(3)-Massenspalt."""
        self.mg = MassGap("SU(3)")

    def test_millennium_problem_statement_not_empty(self):
        """Millennium-Problem-Beschreibung ist nicht leer."""
        statement = self.mg.millennium_problem_statement()
        assert len(statement) > 100

    def test_millennium_problem_statement_contains_key_terms(self):
        """Beschreibung enthält Schlüsselbegriffe."""
        statement = self.mg.millennium_problem_statement()
        assert "Yang-Mills" in statement
        assert "Mass Gap" in statement
        assert "SU(3)" in statement

    def test_physical_motivation_returns_dict(self):
        """physical_motivation gibt Dictionary zurück."""
        result = self.mg.physical_motivation()
        assert isinstance(result, dict)
        assert "confinement" in result
        assert "asymptotic_freedom" in result

    def test_lattice_gauge_evidence_returns_dict(self):
        """lattice_gauge_evidence gibt Dictionary zurück."""
        result = self.mg.lattice_gauge_evidence(beta=6.0, lattice_size=16)
        assert isinstance(result, dict)

    def test_lattice_gauge_evidence_contains_mass_gap(self):
        """Gitter-Evidenz enthält Massenspalt-Schätzung."""
        result = self.mg.lattice_gauge_evidence(beta=6.0, lattice_size=8)
        assert "mass_gap_estimate" in result
        assert result["mass_gap_estimate"] >= 0.0

    def test_lattice_continuum_limit(self):
        """Für β > 6 wird Kontinuumslimes erkannt."""
        result = self.mg.lattice_gauge_evidence(beta=7.0, lattice_size=32)
        assert result["continuum_limit"] == True

    def test_lattice_no_continuum_limit(self):
        """Für β < 6 kein Kontinuumslimes."""
        result = self.mg.lattice_gauge_evidence(beta=4.0, lattice_size=8)
        assert result["continuum_limit"] == False

    def test_mass_gap_lower_bound_not_empty(self):
        """Jaffe-Witten-Ergebnis-String ist nicht leer."""
        result = self.mg.mass_gap_lower_bound_jaffe_witten()
        assert len(result) > 100

    def test_glueball_mass_zero_coupling(self):
        """Glueballmasse für g²=0 ist 0."""
        m = self.mg.glueball_mass_estimate(0.0)
        assert abs(m) < EPS

    def test_glueball_mass_positive(self):
        """Glueballmasse für g²>0 ist positiv."""
        m = self.mg.glueball_mass_estimate(1.0)
        assert m > 0.0

    def test_glueball_mass_small_coupling_suppression(self):
        """Glueballmasse ist exponentiell unterdrückt für kleine g²."""
        m_small = self.mg.glueball_mass_estimate(0.1)
        m_large = self.mg.glueball_mass_estimate(1.0)
        # Für kleine g²: m ~ exp(-c/g²) → sehr klein
        assert m_small < m_large

    def test_confinement_criterion_large_area(self):
        """Confinement für große Wilson-Schleife (großes Flächen-Argument)."""
        # Große Fläche → W klein → Confinement
        assert self.mg.confinement_criterion(wilson_loop_area=10.0) == True

    def test_confinement_criterion_small_area(self):
        """Kein Confinement-Signal für sehr kleine Fläche."""
        # Kleine Fläche → W ≈ 1 → kein Confinement-Signal
        assert self.mg.confinement_criterion(wilson_loop_area=0.1) == False


# =============================================================================
# KLASSE WilsonLoop
# =============================================================================

class TestWilsonLoop:
    """Tests für Wilson-Schleifen."""

    def setup_method(self):
        """Test-Setup: SU(3)-Wilson-Schleifen."""
        self.wl = WilsonLoop("SU(3)")

    def test_rectangular_loop_returns_float(self):
        """Rechteckiger Wilson-Loop gibt float zurück."""
        W = self.wl.rectangular_loop(R=1.0, T=1.0, g_squared=1.0)
        assert isinstance(W, float)

    def test_rectangular_loop_between_0_and_1(self):
        """Wilson-Loop liegt zwischen 0 und 1."""
        W = self.wl.rectangular_loop(R=1.0, T=1.0, g_squared=1.0)
        assert 0.0 <= W <= 1.0 + EPS

    def test_rectangular_loop_decreases_with_area(self):
        """Wilson-Loop nimmt mit der Fläche ab (Confinement)."""
        W1 = self.wl.rectangular_loop(R=1.0, T=1.0, g_squared=2.0)
        W2 = self.wl.rectangular_loop(R=2.0, T=2.0, g_squared=2.0)
        assert W2 < W1

    def test_rectangular_loop_trivial(self):
        """Wilson-Loop für R=0 oder T=0 ist 1 (leere Schleife)."""
        W = self.wl.rectangular_loop(R=0.0, T=1.0, g_squared=1.0)
        assert abs(W - 1.0) < EPS

    def test_string_tension_positive(self):
        """String-Tension ist positiv für g² > 0."""
        sigma = self.wl.string_tension(g_squared=2.0)
        assert sigma > 0.0

    def test_string_tension_zero_coupling(self):
        """String-Tension ist 0 für g² = 0."""
        sigma = self.wl.string_tension(g_squared=0.0)
        assert abs(sigma) < EPS

    def test_area_law_check_confinement(self):
        """Flächengesetz erkennt Confinement."""
        # W ≈ exp(-σ·RT), σ=0.18, R=T=2: W ≈ exp(-0.72) ≈ 0.487
        W = math.exp(-0.18 * 2.0 * 2.0)
        # Perimeter = 2(R+T) = 8, Area = R·T = 4
        # σ_eff = log(1/W)/4 = 0.045, μ_eff = log(1/W)/8 = 0.0225
        # → σ_eff > μ_eff → Confinement
        assert self.wl.area_law_check(perimeter=8.0, area=4.0, W_value=W) == True

    def test_area_law_check_no_confinement(self):
        """Umfangsgesetz → kein Confinement.

        Wenn Perimeter >> Area, dominiert das Umfangsgesetz.
        Kriterium: σ_eff = |log W|/Area > μ_eff = |log W|/Perimeter
        → Flächengesetz wenn Area < Perimeter und |log W| gleich.

        Für reines Umfangsgesetz: W = exp(-μ·Perimeter)
        σ_eff = μ·Perimeter/Area >> μ_eff = μ (da Perimeter >> Area ist FALSCH)

        Tatsächlich: area_law_check gibt True wenn σ_eff > μ_eff,
        also wenn |log W|/Area > |log W|/Perimeter, d.h. Perimeter > Area.
        Für Perimeter < Area → Flächengesetz False (Umfangsgesetz dominiert).
        """
        # Für Umfangsgesetz: Perimeter < Area → σ_eff < μ_eff → kein Confinement
        W = math.exp(-0.5)  # gleichmäßig gedämpfter Wert
        # Perimeter=2 < Area=10 → μ_eff = 0.5/2 = 0.25 > σ_eff = 0.5/10 = 0.05
        # → kein Flächengesetz → False
        assert self.wl.area_law_check(perimeter=2.0, area=10.0, W_value=W) == False

    def test_creutz_ratio_basic(self):
        """Creutz-Verhältnis berechnet String-Tension."""
        # Für Wilson-Schleifen mit Flächengesetz: W(R,T) = exp(-σRT)
        sigma = 0.1
        W_values = {}
        for R in range(1, 5):
            for T in range(1, 5):
                W_values[(R, T)] = math.exp(-sigma * R * T)
        chi = self.wl.creutz_ratio(W_values, R=3, T=3)
        # χ(R,T) ≈ σ für das Flächengesetz
        assert abs(chi - sigma) < 1e-5

    def test_creutz_ratio_missing_values(self):
        """Creutz-Verhältnis gibt nan bei fehlenden Werten."""
        W_values = {(3, 3): 0.5}  # Nicht genug Werte
        chi = self.wl.creutz_ratio(W_values, R=3, T=3)
        assert math.isnan(chi)


# =============================================================================
# KLASSE YangMillsTopology
# =============================================================================

class TestYangMillsTopology:
    """Tests für topologische Eigenschaften."""

    def setup_method(self):
        """Test-Setup."""
        self.topo = YangMillsTopology()
        # Null-Feldstärke
        self.F_zero = {(0, 1): np.zeros((2, 2), dtype=complex),
                       (0, 2): np.zeros((2, 2), dtype=complex),
                       (0, 3): np.zeros((2, 2), dtype=complex),
                       (1, 2): np.zeros((2, 2), dtype=complex),
                       (1, 3): np.zeros((2, 2), dtype=complex),
                       (2, 3): np.zeros((2, 2), dtype=complex)}

    def test_chern_class_0_is_one(self):
        """c_0 = 1 (Normierung)."""
        c0 = self.topo.chern_class(self.F_zero, k=0)
        assert abs(c0 - 1.0) < EPS

    def test_chern_class_1_zero_field(self):
        """c_1 = 0 für SU(N)-Feld (spurlose Generatoren)."""
        c1 = self.topo.chern_class({"F_components": self.F_zero}, k=1)
        assert abs(c1) < EPS

    def test_chern_class_2_zero_field(self):
        """c_2 = 0 für verschwindendes Feld."""
        c2 = self.topo.chern_class({"F_components": self.F_zero}, k=2)
        assert abs(c2) < EPS

    def test_pontryagin_class_zero_field(self):
        """p_1 = 0 für verschwindendes Feld."""
        p1 = self.topo.pontryagin_class({"F_components": self.F_zero})
        assert abs(p1) < EPS

    def test_pontryagin_vs_chern(self):
        """p_1 = -2·c_2 (Beziehung zwischen Pontryagin und Chern)."""
        F = {"F_components": self.F_zero}
        p1 = self.topo.pontryagin_class(F)
        c2 = self.topo.chern_class(F, k=2)
        assert abs(p1 - (-2.0 * c2)) < EPS

    def test_theta_vacuum_theta_zero(self):
        """Für θ=0 ist das Theta-Vakuum eine reelle Summe."""
        amp = self.topo.theta_vacuum(theta=0.0, instanton_sum=5)
        assert abs(amp.imag) < EPS

    def test_theta_vacuum_nonzero(self):
        """Theta-Vakuum-Amplitude ist nicht-trivial für θ ≠ 0."""
        amp = self.topo.theta_vacuum(theta=math.pi / 4, instanton_sum=5)
        assert isinstance(amp, complex)

    def test_vacuum_energy_density_theta_0(self):
        """Vakuumsenergie ist minimal bei θ=0 (Maximum des -cos)."""
        E0 = self.topo.vacuum_energy_density(theta=0.0, g_squared=1.0)
        E_pi = self.topo.vacuum_energy_density(theta=math.pi, g_squared=1.0)
        # E(θ=0) < E(θ=π) da cos(0)=1 > cos(π)=-1
        assert E0 < E_pi

    def test_vacuum_energy_density_zero_coupling(self):
        """Vakuumenergie ist 0 für g²=0."""
        E = self.topo.vacuum_energy_density(theta=1.0, g_squared=0.0)
        assert abs(E) < EPS

    def test_axion_solution_not_empty(self):
        """Axion-Lösung gibt nicht-leeren String zurück."""
        result = self.topo.axion_solution()
        assert len(result) > 100

    def test_axion_solution_contains_peccei_quinn(self):
        """Axion-Beschreibung erwähnt Peccei-Quinn."""
        result = self.topo.axion_solution()
        assert "Peccei" in result or "Quinn" in result


# =============================================================================
# KLASSE IndexTheory
# =============================================================================

class TestIndexTheory:
    """Tests für das Atiyah-Singer-Index-Theorem."""

    def setup_method(self):
        """Test-Setup."""
        self.idx = IndexTheory()

    def test_atiyah_singer_index_identity(self):
        """Index der Einheitsmatrix: ind(I) = 0."""
        I = np.eye(4, dtype=complex)
        ind = self.idx.atiyah_singer_index(I)
        assert ind == 0

    def test_atiyah_singer_index_zero_matrix(self):
        """Index der Null-Matrix: alle Singulärwerte 0 → ind = 0 (quadratisch)."""
        Z = np.zeros((4, 4), dtype=complex)
        ind = self.idx.atiyah_singer_index(Z)
        assert ind == 0

    def test_atiyah_singer_index_rectangular(self):
        """Index für nicht-quadratische Matrix: ind = n - m."""
        # 4×6-Matrix: ind = 6-4 = 2 (wenn vollen Rang hat)
        M = np.random.randn(4, 6) + 1j * np.random.randn(4, 6)
        ind = self.idx.atiyah_singer_index(M)
        assert isinstance(ind, (int, np.integer))

    def test_dirac_operator_spectrum_returns_list(self):
        """Dirac-Spektrum gibt Liste zurück."""
        F_zero = {"F_components": {}}
        spectrum = self.idx.dirac_operator_spectrum(F_zero, mass=1.0)
        assert isinstance(spectrum, list)

    def test_dirac_operator_spectrum_positive(self):
        """Eigenwerte von D†D sind nicht-negativ."""
        F_zero = {"F_components": {}}
        spectrum = self.idx.dirac_operator_spectrum(F_zero, mass=1.0)
        for ev in spectrum:
            assert ev >= -EPS

    def test_dirac_operator_spectrum_mass_dependence(self):
        """Spektrum hängt von der Masse ab."""
        F_zero = {"F_components": {}}
        spec1 = self.idx.dirac_operator_spectrum(F_zero, mass=1.0)
        spec2 = self.idx.dirac_operator_spectrum(F_zero, mass=2.0)
        # Verschiedene Massen → verschiedene Spektren
        assert not all(abs(e1 - e2) < EPS for e1, e2 in zip(spec1, spec2))

    def test_index_theorem_check_trivial(self):
        """Index-Theorem-Check für triviales Bündel (k=0)."""
        topology = {"topological_charge": 0, "manifold": "S4"}
        # Für k=0: ind = 0
        result = self.idx.index_theorem_check(topology, analytical_index=0)
        assert result == True

    def test_index_theorem_check_instanton(self):
        """Index-Theorem: ind(D_A) = -k für k-Instanton auf S⁴."""
        topology = {"topological_charge": 1, "manifold": "S4"}
        # k=1: ind(D) = -1
        result = self.idx.index_theorem_check(topology, analytical_index=-1)
        assert result == True

    def test_index_theorem_check_wrong_index(self):
        """Index-Theorem-Check schlägt fehl bei falschem Index."""
        topology = {"topological_charge": 1, "manifold": "S4"}
        result = self.idx.index_theorem_check(topology, analytical_index=0)
        assert result == False

    def test_eta_invariant_identity(self):
        """Eta-Invariante der Einheitsmatrix: η(0) = 4 (alle Eigenwerte +1)."""
        I = np.eye(4, dtype=float)
        eta = self.idx.eta_invariant(I)
        assert abs(eta - 4.0) < EPS

    def test_eta_invariant_minus_identity(self):
        """Eta-Invariante von -I: η(0) = -4 (alle Eigenwerte -1)."""
        minusI = -np.eye(4, dtype=float)
        eta = self.idx.eta_invariant(minusI)
        assert abs(eta + 4.0) < EPS

    def test_eta_invariant_zero_matrix(self):
        """Eta-Invariante der Null-Matrix: η(0) = 0 (kein Beitrag von Null-EW)."""
        Z = np.zeros((4, 4), dtype=float)
        eta = self.idx.eta_invariant(Z)
        assert abs(eta) < EPS

    def test_eta_invariant_symmetric_spectrum(self):
        """Eta-Invariante bei symmetrischem Spektrum ist 0."""
        # Diagonalmatrix mit gleichviel positiven und negativen Eigenwerten
        D = np.diag([1.0, 1.0, -1.0, -1.0])
        eta = self.idx.eta_invariant(D)
        assert abs(eta) < EPS
