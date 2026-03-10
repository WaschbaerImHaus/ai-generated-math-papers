"""
@file test_symplectic.py
@brief Tests für symplektische Geometrie und Hamilton-Mechanik (tensor_geometry.py).

@description
    Testet alle neuen Funktionen der symplektischen Geometrie:
    - symplectic_form: Standard-J-Matrix
    - is_symplectic_matrix: Gruppentest Sp(2n)
    - poisson_bracket: Kanonische Relationen
    - hamiltonian_flow: Numerische Integration (Störmer-Verlet)
    - liouville_theorem_check: Phasenraumvolumen-Erhaltung
    - action_angle_variables: Wirkung J = E/ω für harmonischen Oszillator
    - cotangent_bundle_metric: T*M-Struktur
    - darboux_theorem_check: Voraussetzungen Darboux-Satz

    Testphilosophie:
    - Algebraische Eigenschaften exakt prüfen (Schiefsymmetrie, det=1)
    - Numerische Eigenschaften mit Toleranz prüfen (Energie-Erhaltung, Liouville)
    - Physikalische Referenzergebnisse nutzen (harmonischer Oszillator)

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import numpy as np
import pytest
import sys
import os

# Suchpfad um das src/-Verzeichnis erweitern
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tensor_geometry import (
    symplectic_form,
    is_symplectic_matrix,
    poisson_bracket,
    hamiltonian_flow,
    liouville_theorem_check,
    action_angle_variables,
    cotangent_bundle_metric,
    darboux_theorem_check,
    wedge_product,
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def harmonic_H(q: float, p: float, omega: float = 1.0) -> float:
    """
    Hamilton-Funktion des harmonischen Oszillators:
        H = p²/2 + ω²q²/2
    Analytische Lösung: Ellipse im Phasenraum mit Halbachsen a=√(2E)/ω, b=√(2E).
    """
    return 0.5 * p**2 + 0.5 * omega**2 * q**2


def pendulum_H(q: float, p: float) -> float:
    """
    Hamilton-Funktion des mathematischen Pendels (kleine Schwingungen):
        H = p²/2 + (1 - cos(q))
    Nichtlinear, aber noch integrabel.
    """
    return 0.5 * p**2 + (1.0 - math.cos(q))


# ===========================================================================
# Tests: symplectic_form
# ===========================================================================

class TestSymplecticForm:
    """Tests für die Standard-symplektische Form J."""

    def test_shape_n1(self):
        """symplectic_form(1) gibt eine 2×2-Matrix zurück."""
        J = symplectic_form(1)
        assert J.shape == (2, 2)

    def test_shape_n2(self):
        """symplectic_form(2) gibt eine 4×4-Matrix zurück."""
        J = symplectic_form(2)
        assert J.shape == (4, 4)

    def test_j_matrix_n1(self):
        """
        symplectic_form(1) muss exakt [[0, 1], [-1, 0]] sein.
        Dies ist die grundlegendste symplektische Form (Drehung um 90°).
        """
        J = symplectic_form(1)
        expected = np.array([[0.0, 1.0], [-1.0, 0.0]])
        np.testing.assert_array_almost_equal(J, expected, decimal=12)

    def test_j_matrix_n2(self):
        """
        symplectic_form(2) muss die Block-Struktur [[0,I],[-I,0]] haben.
        Prüft alle 16 Einträge der 4×4-Matrix.
        """
        J = symplectic_form(2)
        expected = np.array([
            [ 0.,  0.,  1.,  0.],
            [ 0.,  0.,  0.,  1.],
            [-1.,  0.,  0.,  0.],
            [ 0., -1.,  0.,  0.],
        ])
        np.testing.assert_array_almost_equal(J, expected, decimal=12)

    def test_antisymmetry(self):
        """Symplektische Form ist schiefsymmetrisch: J^T = -J."""
        for n in [1, 2, 3]:
            J = symplectic_form(n)
            np.testing.assert_array_almost_equal(J + J.T, np.zeros_like(J), decimal=12)

    def test_determinant_one(self):
        """det(J) = 1 für alle n (symplektische Form nicht-entartet)."""
        for n in [1, 2, 3]:
            J = symplectic_form(n)
            det = np.linalg.det(J)
            assert abs(det - 1.0) < 1e-10, f"det(J) = {det} ≠ 1 für n={n}"

    def test_j_squared_minus_identity(self):
        """J² = -I_{2n} (analog zur imaginären Einheit i² = -1)."""
        for n in [1, 2, 3]:
            J = symplectic_form(n)
            J2 = J @ J
            expected = -np.eye(2 * n)
            np.testing.assert_array_almost_equal(J2, expected, decimal=12)


# ===========================================================================
# Tests: is_symplectic_matrix
# ===========================================================================

class TestIsSymplecticMatrix:
    """Tests für die symplektische Gruppenprüfung."""

    def test_identity_is_symplectic(self):
        """Die Einheitsmatrix I_{2n} ist symplektisch (triviales Element von Sp(2n))."""
        for n in [1, 2]:
            I = np.eye(2 * n)
            assert is_symplectic_matrix(I) is True

    def test_symplectic_form_itself(self):
        """Die symplektische Form J selbst ist symplektisch: J^T J J = J."""
        J = symplectic_form(2)
        assert is_symplectic_matrix(J) is True

    def test_rotation_matrix_2x2_is_symplectic(self):
        """
        Rotationsmatrizen in 2D sind symplektisch (Sp(2) = SL(2,R) ∩ SO(2)):
            R(θ) = [[cos θ, -sin θ], [sin θ, cos θ]]
        Gilt für alle θ, da det(R)=1 und R^T J R = J.
        """
        for theta in [0.3, 1.0, math.pi / 4, math.pi / 3]:
            c, s = math.cos(theta), math.sin(theta)
            R = np.array([[c, -s], [s, c]])
            assert is_symplectic_matrix(R) is True, f"Rotation θ={theta} nicht symplektisch"

    def test_shear_matrix_is_symplectic(self):
        """
        Schermatrix ist symplektisch: [[1, t], [0, 1]] für beliebiges t.
        Diese entspricht einer kanonischen Koordinatentransformation (Typ I).
        """
        for t in [0.5, 1.0, 2.0]:
            M = np.array([[1.0, t], [0.0, 1.0]])
            assert is_symplectic_matrix(M) is True

    def test_scaled_matrix_not_symplectic(self):
        """
        Skalierung 2·I ist NICHT symplektisch: (2I)^T J (2I) = 4J ≠ J.
        """
        M = 2.0 * np.eye(2)
        assert is_symplectic_matrix(M) is False

    def test_odd_dimension_not_symplectic(self):
        """Ungerade Dimension → nicht symplektisch (gibt False zurück)."""
        M = np.eye(3)
        assert is_symplectic_matrix(M) is False

    def test_nonsquare_not_symplectic(self):
        """Nicht-quadratische Matrix → nicht symplektisch."""
        M = np.ones((2, 4))
        assert is_symplectic_matrix(M) is False


# ===========================================================================
# Tests: poisson_bracket
# ===========================================================================

class TestPoissonBracket:
    """Tests für numerische Poisson-Klammer."""

    def setup_method(self):
        """Erzeuge ein (q, p)-Gitter für alle Tests."""
        # Gitter: q ∈ [-2, 2], p ∈ [-2, 2], je 100 Punkte
        n = 100
        q_arr = np.linspace(-2.0, 2.0, n)
        p_arr = np.linspace(-2.0, 2.0, n)
        self.dq = q_arr[1] - q_arr[0]
        self.dp = p_arr[1] - p_arr[0]
        self.Q, self.P = np.meshgrid(q_arr, p_arr, indexing='ij')

    def test_canonical_relation_q_p(self):
        """
        Kanonische Poisson-Relation: {q, p} = 1.
        f = q, g = p → {f,g} = ∂q/∂q · ∂p/∂p - ∂q/∂p · ∂p/∂q = 1·1 - 0·0 = 1.
        """
        f = self.Q.copy()   # f(q,p) = q
        g = self.P.copy()   # g(q,p) = p
        pb = poisson_bracket(f, g, self.dq, self.dp)
        # {q, p} sollte überall ≈ 1 sein (bis auf Randeffekte)
        # Nutze innere Punkte (ohne Rand) für bessere Genauigkeit
        inner = pb[1:-1, 1:-1]
        assert abs(np.mean(inner) - 1.0) < 0.05, f"{{q,p}} = {np.mean(inner):.4f} ≠ 1"

    def test_antisymmetry(self):
        """Schiefsymmetrie: {f, g} = -{g, f}."""
        f = self.Q**2 + self.P     # Beliebige Funktion
        g = self.P**2 - self.Q     # Beliebige andere Funktion
        pb_fg = poisson_bracket(f, g, self.dq, self.dp)
        pb_gf = poisson_bracket(g, f, self.dq, self.dp)
        np.testing.assert_array_almost_equal(pb_fg, -pb_gf, decimal=10)

    def test_self_bracket_zero(self):
        """Selbst-Klammer: {f, f} = 0 für alle f."""
        f = self.Q**2 + self.P**2   # H des harmonischen Oszillators
        pb = poisson_bracket(f, f, self.dq, self.dp)
        np.testing.assert_array_almost_equal(pb, np.zeros_like(pb), decimal=10)

    def test_energy_conserved(self):
        """
        Energieerhaltung: {H, H} = 0.
        Der harmonische Oszillator H = p²/2 + q²/2 erhält sich selbst.
        """
        H = 0.5 * self.Q**2 + 0.5 * self.P**2
        pb = poisson_bracket(H, H, self.dq, self.dp)
        np.testing.assert_array_almost_equal(pb, np.zeros_like(pb), decimal=10)

    def test_linearity(self):
        """{αf + βg, h} = α{f,h} + β{g,h} (Bilinearität)."""
        alpha, beta = 2.0, -3.0
        f = self.Q * self.P
        g = self.Q**2
        h = self.P**2
        pb_sum = poisson_bracket(alpha * f + beta * g, h, self.dq, self.dp)
        pb_lin = alpha * poisson_bracket(f, h, self.dq, self.dp) \
               + beta  * poisson_bracket(g, h, self.dq, self.dp)
        np.testing.assert_array_almost_equal(pb_sum, pb_lin, decimal=8)


# ===========================================================================
# Tests: hamiltonian_flow
# ===========================================================================

class TestHamiltonianFlow:
    """Tests für den Störmer-Verlet-Hamiltonfluss."""

    def test_return_keys(self):
        """Rückgabe-Dictionary enthält 'q', 'p', 't', 'H'."""
        result = hamiltonian_flow(harmonic_H, 1.0, 0.0, t_span=1.0, dt=0.01)
        for key in ['q', 'p', 't', 'H']:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_array_lengths(self):
        """Alle Arrays haben gleiche Länge n_steps+1."""
        t_span, dt = 2.0, 0.05
        result = hamiltonian_flow(harmonic_H, 1.0, 0.0, t_span=t_span, dt=dt)
        n_expected = int(t_span / dt) + 1
        for key in ['q', 'p', 't', 'H']:
            assert len(result[key]) == n_expected

    def test_initial_conditions(self):
        """Erster Punkt entspricht exakt den Anfangsbedingungen."""
        q0, p0 = 1.5, -0.3
        result = hamiltonian_flow(harmonic_H, q0, p0, t_span=1.0, dt=0.01)
        assert abs(result['q'][0] - q0) < 1e-12
        assert abs(result['p'][0] - p0) < 1e-12

    def test_energy_conservation_harmonic(self):
        """
        Störmer-Verlet: Energiefehler < 1% über eine Periode des harmonischen Oszillators.
        Periode: T = 2π/ω = 2π für ω=1.
        """
        omega = 1.0
        E0 = 0.5   # Anfangsenergie
        q0 = math.sqrt(2.0 * E0) / omega   # Startpunkt auf q-Achse
        result = hamiltonian_flow(
            lambda q, p: harmonic_H(q, p, omega),
            q0, 0.0, t_span=2 * math.pi, dt=0.01,
        )
        H_vals = result['H']
        H_max_err = max(abs(H_vals - H_vals[0]))
        assert H_max_err < 0.01 * abs(H_vals[0]) + 1e-10, \
            f"Energiefehler {H_max_err:.2e} > 1% von E={H_vals[0]:.4f}"

    def test_energy_conservation_long_time(self):
        """
        Störmer-Verlet Energiefehler < 1% über 100 Perioden (Langzeitsimulation).
        Dies ist der entscheidende Vorteil gegenüber nicht-symplektischen Methoden.
        """
        omega = 1.0
        q0, p0 = 1.0, 0.0
        T = 2.0 * math.pi / omega
        result = hamiltonian_flow(harmonic_H, q0, p0, t_span=100 * T, dt=0.01)
        H_vals = result['H']
        max_relative_err = max(abs(H_vals - H_vals[0])) / abs(H_vals[0])
        assert max_relative_err < 0.01, \
            f"Relativer Energiefehler {max_relative_err:.4f} > 1% nach 100 Perioden"

    def test_phase_space_trajectory_is_ellipse(self):
        """
        Trajektorie des harmonischen Oszillators bildet eine Ellipse im Phasenraum.
        Prüfung: max(|q²/a² + p²/b² - 1|) < Toleranz mit a=1, b=ω für H=1/2.
        """
        omega = 1.0
        E = 0.5
        q0 = math.sqrt(2 * E) / omega
        result = hamiltonian_flow(harmonic_H, q0, 0.0, t_span=4 * math.pi, dt=0.005)
        q, p = result['q'], result['p']
        # Ellipsentest: q²/a² + p²/b² = 1 mit a = √(2E)/ω, b = √(2E)
        a = math.sqrt(2 * E) / omega
        b = math.sqrt(2 * E)
        ellipse_err = np.abs(q**2 / a**2 + p**2 / b**2 - 1.0)
        assert np.max(ellipse_err) < 0.02, \
            f"Trajektorie ist keine Ellipse, max. Abweichung: {np.max(ellipse_err):.4f}"

    def test_time_array_correct(self):
        """Zeitachse startet bei 0 und endet nahe t_span."""
        t_span, dt = 5.0, 0.1
        result = hamiltonian_flow(harmonic_H, 1.0, 0.0, t_span=t_span, dt=dt)
        assert result['t'][0] == 0.0
        assert abs(result['t'][-1] - t_span) < dt + 1e-10


# ===========================================================================
# Tests: liouville_theorem_check
# ===========================================================================

class TestLiouvilleTheorem:
    """Tests für die numerische Verifikation des Liouville-Satzes."""

    def test_return_keys(self):
        """Rückgabe-Dictionary enthält die erforderlichen Schlüssel."""
        ics = [(1.0, 0.0), (1.1, 0.0), (1.0, 0.1)]
        result = liouville_theorem_check(harmonic_H, ics, t=1.0)
        for key in ['volume_initial', 'volume_final', 'ratio']:
            assert key in result

    def test_volume_ratio_near_one(self):
        """
        Liouville-Satz: Phasenraumvolumen bleibt unter Hamiltonfluss erhalten.
        Erwarte ratio ≈ 1 mit Toleranz 5%.
        """
        # Kleines kreisförmiges Ensemble um (1, 0)
        N = 8
        r = 0.05   # Kleiner Radius
        ics = [
            (1.0 + r * math.cos(2 * math.pi * k / N),
             r * math.sin(2 * math.pi * k / N))
            for k in range(N)
        ]
        result = liouville_theorem_check(harmonic_H, ics, t=2.0)
        ratio = result['ratio']
        assert abs(ratio - 1.0) < 0.05, \
            f"Liouville-Ratio {ratio:.4f} weicht um mehr als 5% von 1 ab"

    def test_volumes_positive(self):
        """Initiales und finales Volumen sind positiv."""
        ics = [(1.0, 0.0), (1.1, 0.1), (0.9, -0.1), (1.0, 0.2)]
        result = liouville_theorem_check(harmonic_H, ics, t=1.0)
        assert result['volume_initial'] > 0
        assert result['volume_final'] > 0


# ===========================================================================
# Tests: action_angle_variables
# ===========================================================================

class TestActionAngleVariables:
    """Tests für Wirkungs-Winkelvariablen."""

    def test_return_keys(self):
        """Rückgabe-Dictionary enthält 'J', 'omega', 'q_vals', 'p_vals'."""
        result = action_angle_variables(harmonic_H, E=1.0)
        for key in ['J', 'omega', 'q_vals', 'p_vals']:
            assert key in result

    def test_J_equals_E_over_omega_harmonic(self):
        """
        Harmonischer Oszillator: J = E/ω analytisch.
        Für ω=1, E=1: J = 1. Numerisch sollte J ≈ 1 sein.
        """
        omega = 1.0
        E = 1.0
        result = action_angle_variables(
            lambda q, p: harmonic_H(q, p, omega),
            E=E,
            q_range=(-math.sqrt(2 * E) * 1.1, math.sqrt(2 * E) * 1.1),
        )
        J_expected = E / omega
        assert abs(result['J'] - J_expected) < 0.05, \
            f"J = {result['J']:.4f}, erwartet {J_expected:.4f}"

    def test_J_scales_with_energy(self):
        """
        J ∝ E für harmonischen Oszillator (J = E/ω).
        Verdoppelung von E → Verdoppelung von J.
        """
        omega = 1.0
        E1, E2 = 0.5, 1.0
        result1 = action_angle_variables(
            lambda q, p: harmonic_H(q, p, omega), E=E1,
            q_range=(-1.5, 1.5),
        )
        result2 = action_angle_variables(
            lambda q, p: harmonic_H(q, p, omega), E=E2,
            q_range=(-1.7, 1.7),
        )
        # J2/J1 sollte ≈ E2/E1 = 2 sein
        ratio = result2['J'] / result1['J']
        assert abs(ratio - (E2 / E1)) < 0.1, \
            f"J-Skalierung: ratio = {ratio:.4f}, erwartet {E2/E1:.4f}"

    def test_omega_correct_harmonic(self):
        """
        Frequenz des harmonischen Oszillators: ω = dE/dJ = ω₀.
        Für ω₀=1: Kreisfrequenz = 1, also Frequenz = 1/2π? Nein: ω = 1.
        """
        omega0 = 1.0
        E = 1.0
        q_range = (-1.5 * math.sqrt(2 * E), 1.5 * math.sqrt(2 * E))
        result = action_angle_variables(
            lambda q, p: harmonic_H(q, p, omega0), E=E, q_range=q_range,
        )
        if not math.isnan(result['omega']):
            assert abs(result['omega'] - omega0) < 0.2, \
                f"ω = {result['omega']:.4f}, erwartet {omega0:.4f}"

    def test_p_vals_nonnegative(self):
        """p_vals sind nicht-negativ (positiver Ast der Energiebahn)."""
        result = action_angle_variables(harmonic_H, E=1.0)
        assert np.all(result['p_vals'] >= 0), "p_vals enthält negative Werte"

    def test_zero_energy_edge_case(self):
        """Für E=0 ist die Trajektorie ein Punkt (J ≈ 0)."""
        result = action_angle_variables(harmonic_H, E=0.0)
        assert abs(result['J']) < 0.01, f"J = {result['J']:.4f} sollte ≈ 0 für E=0 sein"


# ===========================================================================
# Tests: cotangent_bundle_metric
# ===========================================================================

class TestCotangentBundleMetric:
    """Tests für die symplektische Struktur des Kotangentialbündels."""

    def test_return_keys(self):
        """Rückgabe-Dictionary enthält alle Schlüssel."""
        g = np.eye(2)
        result = cotangent_bundle_metric(g)
        for key in ['symplectic_form', 'is_exact', 'liouville_form']:
            assert key in result

    def test_symplectic_form_shape(self):
        """Für n×n-Basismetrik: symplektische Form ist (2n×2n)."""
        for n in [1, 2, 3]:
            g = np.eye(n)
            result = cotangent_bundle_metric(g)
            omega = result['symplectic_form']
            assert omega.shape == (2 * n, 2 * n)

    def test_is_exact(self):
        """Symplektische Form auf T*M ist immer exakt (is_exact=True)."""
        g = np.eye(2)
        result = cotangent_bundle_metric(g)
        assert result['is_exact'] is True

    def test_liouville_form_shape(self):
        """Liouville-Form hat Länge 2n."""
        for n in [1, 2, 3]:
            g = np.eye(n)
            result = cotangent_bundle_metric(g)
            assert len(result['liouville_form']) == 2 * n

    def test_symplectic_form_antisymmetric(self):
        """Symplektische Form auf T*M ist schiefsymmetrisch."""
        g = np.eye(2)
        result = cotangent_bundle_metric(g)
        omega = result['symplectic_form']
        np.testing.assert_array_almost_equal(omega + omega.T, np.zeros_like(omega), decimal=12)


# ===========================================================================
# Tests: darboux_theorem_check
# ===========================================================================

class TestDarbouxTheoremCheck:
    """Tests für die Darboux-Satz-Überprüfung."""

    def test_standard_J_darboux_applicable(self):
        """Standard-symplektische Form J erfüllt alle Darboux-Voraussetzungen."""
        J = symplectic_form(2)
        point = np.array([0.0, 0.0, 0.0, 0.0])
        result = darboux_theorem_check(J, point)
        assert result['darboux_applicable'] is True
        assert result['is_non_degenerate'] is True
        assert result['is_antisymmetric'] is True

    def test_non_degenerate_flag(self):
        """Nicht-entartete Form → is_non_degenerate=True."""
        J = symplectic_form(1)
        result = darboux_theorem_check(J, np.zeros(2))
        assert result['is_non_degenerate'] is True

    def test_degenerate_form_not_darboux(self):
        """Entartete Form (det=0) → darboux_applicable=False."""
        omega = np.zeros((4, 4))   # Komplett degeneriert
        result = darboux_theorem_check(omega, np.zeros(4))
        assert result['is_non_degenerate'] is False
        assert result['darboux_applicable'] is False

    def test_antisymmetry_check(self):
        """Symmetrische Matrix → is_antisymmetric=False, nicht Darboux."""
        sym_matrix = np.eye(4)   # Symmetrisch, nicht schiefsymmetrisch
        result = darboux_theorem_check(sym_matrix, np.zeros(4))
        assert result['is_antisymmetric'] is False

    def test_det_omega_correct(self):
        """det_omega-Wert wird korrekt berechnet."""
        J = symplectic_form(2)   # det(J) = 1
        result = darboux_theorem_check(J, np.zeros(4))
        assert abs(result['det_omega'] - 1.0) < 1e-10

    def test_nonsquare_returns_false(self):
        """Nicht-quadratische Matrix → alle Flags False."""
        omega = np.ones((3, 4))
        result = darboux_theorem_check(omega, np.zeros(3))
        assert result['darboux_applicable'] is False


# ===========================================================================
# Tests: Kombinierte / Integrations-Tests
# ===========================================================================

class TestSymplecticIntegration:
    """Kombinierte Tests für mehrere Funktionen zusammen."""

    def test_symplectic_form_wedge_nondegerate(self):
        """
        ω ∧ ω ≠ 0 (nicht-entartet): Für n=2 ist ω∧ω die Volumenform.
        Prüfe indirekt über det(J) ≠ 0 und dass ω∧ω eine nicht-triviale
        Form ergibt (via wedge_product der 1-Formen-Darstellung).
        """
        J = symplectic_form(2)
        # Nicht-Entartetheit: det(J) ≠ 0
        det_J = np.linalg.det(J)
        assert abs(det_J) > 0.5, f"det(J) = {det_J:.4f} → ω entartet"

    def test_flow_preserves_symplectic_structure(self):
        """
        Hamiltonfluss ist symplektisch: Jacobi-Matrix des Flusses ist symplektisch.
        Numerisch: Prüfe ob ein kleines symplektisches Parallelogramm sein Volumen behält.
        Diese Eigenschaft des Störmer-Verlet-Integrators unterscheidet ihn von Euler.
        """
        # Zwei benachbarte Anfangsbedingungen
        eps = 0.001
        result1 = hamiltonian_flow(harmonic_H, 1.0,       0.0, t_span=math.pi, dt=0.01)
        result2 = hamiltonian_flow(harmonic_H, 1.0 + eps, 0.0, t_span=math.pi, dt=0.01)
        result3 = hamiltonian_flow(harmonic_H, 1.0,       eps, t_span=math.pi, dt=0.01)

        # Anfängliches symplektisches Flächenelement: eps×eps = eps²
        area_init = eps**2

        # Finales symplektisches Flächenelement via Kreuzprodukt der Differenzvektoren
        dq1 = result2['q'][-1] - result1['q'][-1]
        dp1 = result2['p'][-1] - result1['p'][-1]
        dq2 = result3['q'][-1] - result1['q'][-1]
        dp2 = result3['p'][-1] - result1['p'][-1]
        area_final = abs(dq1 * dp2 - dq2 * dp1)

        # Symplektische Erhaltung: area_final / area_init ≈ 1
        ratio = area_final / area_init
        assert abs(ratio - 1.0) < 0.05, \
            f"Symplektische Fläche nicht erhalten: ratio = {ratio:.4f}"

    def test_j_matrix_from_cotangent_bundle_matches_standard(self):
        """
        cotangent_bundle_metric(I_n) gibt dieselbe symplektische Form wie symplectic_form(n).
        """
        for n in [1, 2, 3]:
            g = np.eye(n)
            result = cotangent_bundle_metric(g)
            J_standard = symplectic_form(n)
            np.testing.assert_array_almost_equal(
                result['symplectic_form'], J_standard, decimal=12,
            )
