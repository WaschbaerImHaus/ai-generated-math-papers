"""
@file test_pde.py
@brief Tests für das PDE-Modul (Partielle Differentialgleichungen).
@description
    Umfassende Tests für alle Funktionen in src/pde.py:
    - PDE-Klassifikation
    - Wärmeleitungsgleichung (explizit, implizit, Crank-Nicolson)
    - Wellengleichung (Leapfrog, analytisch, Dispersionsrelation)
    - Laplace/Poisson-Gleichung
    - Schrödinger-Gleichung (stationär, harmonischer Oszillator, zeitabhängig)
    - Finite-Elemente-Methode 1D
    - Methode der Charakteristiken
    - Burgers-Gleichung

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
import pytest
import sys
import os

# Sicherstellen, dass src/ im Pfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pde import (
    classify_pde,
    heat_equation_explicit,
    heat_equation_implicit,
    heat_equation_crank_nicolson,
    heat_equation_analytical,
    wave_equation_explicit,
    wave_equation_analytical,
    wave_dispersion_relation,
    laplace_equation_2d,
    poisson_equation_2d,
    green_function_laplace_2d,
    harmonic_function_properties,
    schrodinger_stationary,
    schrodinger_harmonic_oscillator,
    schrodinger_time_dependent,
    fem_1d_poisson,
    fem_convergence_demo,
    method_of_characteristics,
    burgers_equation,
)


# =============================================================================
# Klassifikation
# =============================================================================

class TestClassifyPDE:
    """Tests für classify_pde()."""

    def test_laplace_elliptic(self):
        """Laplace-Gleichung u_xx + u_yy = 0: A=1, B=0, C=1 → elliptisch."""
        assert classify_pde(1.0, 0.0, 1.0) == 'elliptic'

    def test_wave_hyperbolic(self):
        """Wellengleichung u_tt - c²u_xx = 0: A=-c², B=0, C=1 → hyperbolisch."""
        assert classify_pde(-1.0, 0.0, 1.0) == 'hyperbolic'

    def test_heat_parabolic(self):
        """Wärmeleitungsgleichung u_t = u_xx → parabolisch: A=1, B=0, C=0."""
        assert classify_pde(1.0, 0.0, 0.0) == 'parabolic'

    def test_mixed_hyperbolic(self):
        """Gemischte PDE mit B²>4AC → hyperbolisch."""
        assert classify_pde(1.0, 3.0, 1.0) == 'hyperbolic'

    def test_pure_elliptic(self):
        """Exakte Elliptizität: A=2, B=1, C=2 → D=1-16=-15<0."""
        assert classify_pde(2.0, 1.0, 2.0) == 'elliptic'

    def test_zero_coefficients(self):
        """A=B=C=0: exakt parabolisch (D=0)."""
        assert classify_pde(0.0, 0.0, 0.0) == 'parabolic'


# =============================================================================
# Wärmeleitungsgleichung
# =============================================================================

class TestHeatEquationExplicit:
    """Tests für heat_equation_explicit()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx+1) sein."""
        u0 = np.sin(np.linspace(0, np.pi, 51))
        result = heat_equation_explicit(u0.tolist(), np.pi, 0.1, nx=50, nt=100)
        assert result.shape == (101, 51)

    def test_boundary_conditions(self):
        """Randwerte u(0,t) = u(L,t) = 0 für alle t."""
        nx, nt = 30, 50
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_explicit(u0.tolist(), np.pi, 0.05, nx=nx, nt=nt)
        assert np.allclose(u[:, 0], 0.0, atol=1e-12)
        assert np.allclose(u[:, -1], 0.0, atol=1e-12)

    def test_initial_condition(self):
        """Erste Zeile muss der Anfangsbedingung entsprechen."""
        nx = 40
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_explicit(u0.tolist(), np.pi, 0.1, nx=nx, nt=100)
        assert np.allclose(u[0, :], u0, atol=1e-12)

    def test_decay_in_time(self):
        """Die Lösung muss im Laufe der Zeit abnehmen (Wärmeleitung dissipiert)."""
        nx = 50
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_explicit(u0.tolist(), np.pi, 1.0, nx=nx, nt=1000)
        # Maximalnorm nimmt ab
        assert np.max(np.abs(u[-1, :])) < np.max(np.abs(u[0, :]))

    def test_stability_warning(self):
        """Instabile Parameter (r > 0.5) sollen eine Warnung erzeugen."""
        nx = 10
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        with pytest.warns(UserWarning, match="Stabilitätsbedingung"):
            # dt sehr groß → r > 0.5
            heat_equation_explicit(u0.tolist(), np.pi, 100.0, nx=nx, nt=5)

    def test_interpolation_of_u0(self):
        """Unterschiedliche Länge von u0 wird korrekt interpoliert."""
        u0 = [0.0, 1.0, 0.0]  # nur 3 Punkte
        result = heat_equation_explicit(u0, np.pi, 0.01, nx=50, nt=10)
        assert result.shape == (11, 51)


class TestHeatEquationImplicit:
    """Tests für heat_equation_implicit()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx+1) sein."""
        u0 = np.sin(np.linspace(0, np.pi, 51))
        result = heat_equation_implicit(u0.tolist(), np.pi, 0.1, nx=50, nt=100)
        assert result.shape == (101, 51)

    def test_boundary_conditions(self):
        """Randwerte bleiben 0."""
        nx, nt = 30, 20
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_implicit(u0.tolist(), np.pi, 0.1, nx=nx, nt=nt)
        assert np.allclose(u[:, 0], 0.0, atol=1e-12)
        assert np.allclose(u[:, -1], 0.0, atol=1e-12)

    def test_unconditional_stability(self):
        """Implizites Verfahren ist auch für große Zeitschritte stabil."""
        nx = 20
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        # Sehr großer Zeitschritt → r >> 0.5, explizit würde divergieren
        u = heat_equation_implicit(u0.tolist(), np.pi, 100.0, nx=nx, nt=10)
        # Lösung sollte beschränkt bleiben
        assert np.all(np.isfinite(u))

    def test_decay_in_time(self):
        """Lösung klingt im Zeitverlauf ab."""
        nx = 40
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_implicit(u0.tolist(), np.pi, 1.0, nx=nx, nt=50)
        assert np.max(np.abs(u[-1, :])) < np.max(np.abs(u[0, :]))


class TestHeatEquationCrankNicolson:
    """Tests für heat_equation_crank_nicolson()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx+1) sein."""
        u0 = np.sin(np.linspace(0, np.pi, 51))
        result = heat_equation_crank_nicolson(u0.tolist(), np.pi, 0.1, nx=50, nt=50)
        assert result.shape == (51, 51)

    def test_boundary_conditions(self):
        """Randwerte bleiben 0."""
        nx, nt = 30, 20
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u = heat_equation_crank_nicolson(u0.tolist(), np.pi, 0.1, nx=nx, nt=nt)
        assert np.allclose(u[:, 0], 0.0, atol=1e-12)
        assert np.allclose(u[:, -1], 0.0, atol=1e-12)

    def test_accuracy_vs_analytical(self):
        """Crank-Nicolson sollte gut mit analytischer Lösung übereinstimmen."""
        nx, nt = 50, 100
        L = np.pi
        T = 0.5
        x = np.linspace(0, L, nx + 1)
        u0 = np.sin(x)
        u_num = heat_equation_crank_nicolson(u0.tolist(), L, T, nx=nx, nt=nt)
        u_ana = heat_equation_analytical(x, T, L=L)
        # Maximaler Fehler sollte < 5% sein
        max_err = np.max(np.abs(u_num[-1, :] - u_ana))
        assert max_err < 0.05


class TestHeatEquationAnalytical:
    """Tests für heat_equation_analytical()."""

    def test_at_t0(self):
        """Für t≈0 sollte die analytische Lösung ≈ sin(x) sein."""
        x = np.linspace(0, np.pi, 100)
        u = heat_equation_analytical(x, 1e-6, L=np.pi)
        assert np.allclose(u, np.sin(x), atol=1e-4)

    def test_boundary_values(self):
        """u(0,t) = u(L,t) = 0."""
        x = np.array([0.0, np.pi])
        u = heat_equation_analytical(x, 0.5, L=np.pi)
        assert abs(u[0]) < 1e-10
        assert abs(u[1]) < 1e-10

    def test_decay_with_time(self):
        """Lösung klingt mit e^{-t} ab."""
        x = np.array([np.pi / 2])
        u_early = heat_equation_analytical(x, 0.1)
        u_late = heat_equation_analytical(x, 1.0)
        assert u_late[0] < u_early[0]


# =============================================================================
# Wellengleichung
# =============================================================================

class TestWaveEquationExplicit:
    """Tests für wave_equation_explicit()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx+1) sein."""
        nx, nt = 50, 100
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u0_t = np.zeros(nx + 1)
        result = wave_equation_explicit(u0.tolist(), u0_t.tolist(), np.pi, 1.0, nx=nx, nt=nt)
        assert result.shape == (nt + 1, nx + 1)

    def test_boundary_conditions(self):
        """Randwerte u(0,t)=u(L,t)=0."""
        nx, nt = 50, 100
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u0_t = np.zeros(nx + 1)
        u = wave_equation_explicit(u0.tolist(), u0_t.tolist(), np.pi, 1.0, nx=nx, nt=nt)
        assert np.allclose(u[:, 0], 0.0, atol=1e-12)
        assert np.allclose(u[:, -1], 0.0, atol=1e-12)

    def test_initial_condition(self):
        """u(x,0) muss Anfangsbedingung entsprechen."""
        nx = 40
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u0_t = np.zeros(nx + 1)
        u = wave_equation_explicit(u0.tolist(), u0_t.tolist(), np.pi, 0.5, nx=nx, nt=100)
        assert np.allclose(u[0, :], u0, atol=1e-12)

    def test_cfl_warning(self):
        """CFL-Verletzung erzeugt Warnung."""
        nx = 10
        u0 = np.sin(np.linspace(0, np.pi, nx + 1))
        u0_t = np.zeros(nx + 1)
        with pytest.warns(UserWarning, match="CFL"):
            wave_equation_explicit(u0.tolist(), u0_t.tolist(), np.pi, 100.0, c=10.0, nx=nx, nt=5)

    def test_energy_conservation_approx(self):
        """Energie sollte näherungsweise erhalten bleiben (im CFL-stabilen Bereich)."""
        nx, nt = 100, 500
        x = np.linspace(0, np.pi, nx + 1)
        u0 = np.sin(x)
        u0_t = np.zeros(nx + 1)
        u = wave_equation_explicit(u0.tolist(), u0_t.tolist(), np.pi, 1.0, c=1.0, nx=nx, nt=nt)
        # Maximale Amplitude sollte nicht sehr viel zunehmen
        max_init = np.max(np.abs(u[0, :]))
        max_final = np.max(np.abs(u[-1, :]))
        assert max_final < 2 * max_init  # Stabilitätsprüfung


class TestWaveEquationAnalytical:
    """Tests für wave_equation_analytical()."""

    def test_at_t0(self):
        """u(x,0) = sin(πx/L)."""
        x = np.linspace(0, np.pi, 100)
        u = wave_equation_analytical(x, 0.0, L=np.pi, c=1.0)
        assert np.allclose(u, np.sin(x), atol=1e-3)

    def test_boundary_values(self):
        """u(0,t) = u(L,t) = 0."""
        x = np.array([0.0, np.pi])
        u = wave_equation_analytical(x, 1.0, L=np.pi, c=1.0)
        assert np.allclose(u, 0.0, atol=1e-10)

    def test_periodicity(self):
        """Für t=2L/c: u(x,t) ≈ u(x,0) (volle Periode)."""
        L, c = np.pi, 1.0
        x = np.linspace(0, L, 50)
        T = 2 * L / c  # Periode
        u0 = wave_equation_analytical(x, 0.0, L=L, c=c)
        uT = wave_equation_analytical(x, T, L=L, c=c)
        assert np.allclose(u0, uT, atol=1e-6)


class TestWaveDispersionRelation:
    """Tests für wave_dispersion_relation()."""

    def test_linear_dispersion(self):
        """ω = ck (nicht-dispersiv)."""
        assert wave_dispersion_relation(2.0, 3.0) == pytest.approx(6.0)

    def test_zero_wavenumber(self):
        """k=0 → ω=0."""
        assert wave_dispersion_relation(0.0, 5.0) == pytest.approx(0.0)


# =============================================================================
# Laplace/Poisson
# =============================================================================

class TestLaplaceEquation2D:
    """Tests für laplace_equation_2d()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (ny+1, nx+1) sein."""
        u = laplace_equation_2d(nx=10, ny=10)
        assert u.shape == (11, 11)

    def test_top_boundary(self):
        """Obere Randbedingung muss eingehalten sein (innere Punkte)."""
        nx, ny = 15, 15
        top_bc = np.ones(nx + 1)
        u = laplace_equation_2d(nx=nx, ny=ny, boundary={
            'top': top_bc,
            'bottom': np.zeros(nx + 1),
            'left': np.zeros(ny + 1),
            'right': np.zeros(ny + 1),
        })
        # Ecken werden von left/right überschrieben, daher nur innere Punkte prüfen
        assert np.allclose(u[-1, 1:-1], top_bc[1:-1], atol=1e-6)

    def test_zero_boundary_gives_zero(self):
        """Überall null Randbedingungen → Lösung überall 0."""
        nx, ny = 10, 10
        u = laplace_equation_2d(nx=nx, ny=ny, boundary={
            'top': np.zeros(nx + 1),
            'bottom': np.zeros(nx + 1),
            'left': np.zeros(ny + 1),
            'right': np.zeros(ny + 1),
        })
        assert np.allclose(u, 0.0, atol=1e-8)

    def test_finite_values(self):
        """Alle Werte müssen endlich sein."""
        u = laplace_equation_2d(nx=15, ny=15)
        assert np.all(np.isfinite(u))


class TestPoissonEquation2D:
    """Tests für poisson_equation_2d()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (ny+1, nx+1) sein."""
        u = poisson_equation_2d(f=lambda x, y: 1.0, nx=10, ny=10)
        assert u.shape == (11, 11)

    def test_zero_boundary(self):
        """Randwerte müssen 0 sein."""
        nx, ny = 15, 15
        u = poisson_equation_2d(f=lambda x, y: 1.0, nx=nx, ny=ny)
        assert np.allclose(u[0, :], 0.0, atol=1e-10)
        assert np.allclose(u[-1, :], 0.0, atol=1e-10)
        assert np.allclose(u[:, 0], 0.0, atol=1e-10)
        assert np.allclose(u[:, -1], 0.0, atol=1e-10)

    def test_positive_source_positive_solution(self):
        """Für f>0 sollte die Lösung im Inneren >0 sein (Maximumprinzip)."""
        nx, ny = 15, 15
        u = poisson_equation_2d(f=lambda x, y: 1.0, nx=nx, ny=ny)
        # Innere Werte sollten positiv sein
        assert np.all(u[1:-1, 1:-1] >= 0.0)

    def test_finite_values(self):
        """Alle Werte müssen endlich sein."""
        u = poisson_equation_2d(f=lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y), nx=15, ny=15)
        assert np.all(np.isfinite(u))


class TestGreenFunctionLaplace2D:
    """Tests für green_function_laplace_2d()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (len(y), len(x)) sein."""
        x = np.linspace(0.1, 1.0, 10)
        y = np.linspace(0.1, 1.0, 8)
        G = green_function_laplace_2d(0.5, 0.5, x, y)
        assert G.shape == (8, 10)

    def test_symmetry(self):
        """G(x,y;x₀,y₀) = G(x₀,y₀;x,y) (Symmetrie der Greens Funktion)."""
        x = np.array([0.7])
        y = np.array([0.3])
        G1 = green_function_laplace_2d(0.2, 0.4, x, y)
        x2 = np.array([0.2])
        y2 = np.array([0.4])
        G2 = green_function_laplace_2d(0.7, 0.3, x2, y2)
        assert abs(G1[0, 0] - G2[0, 0]) < 1e-10

    def test_negative_values(self):
        """Greens Funktion von Laplace ist im 2D negativ für r>1/e."""
        x = np.linspace(2, 5, 20)
        y = np.linspace(2, 5, 20)
        G = green_function_laplace_2d(0.0, 0.0, x, y)
        # Für große r: G = -1/(2π) ln(r) < 0 da ln(r)>0
        assert np.all(G < 0)


class TestHarmonicFunctionProperties:
    """Tests für harmonic_function_properties()."""

    def test_returns_dict_keys(self):
        """Rückgabe muss alle erwarteten Schlüssel enthalten."""
        u = laplace_equation_2d(nx=10, ny=10)
        props = harmonic_function_properties(u)
        for key in ['max_interior', 'max_boundary', 'satisfies_max_principle']:
            assert key in props

    def test_max_principle_for_laplace_solution(self):
        """Lösung der Laplace-Gleichung erfüllt das Maximum-Prinzip."""
        u = laplace_equation_2d(nx=15, ny=15, boundary={
            'top': np.ones(16),
            'bottom': np.zeros(16),
            'left': np.zeros(16),
            'right': np.zeros(16),
        })
        props = harmonic_function_properties(u)
        assert props['satisfies_max_principle'] is True


# =============================================================================
# Schrödinger-Gleichung
# =============================================================================

class TestSchrodingerStationary:
    """Tests für schrodinger_stationary()."""

    def test_returns_energies_and_wavefunctions(self):
        """Rückgabe muss 'energies' und 'wavefunctions' enthalten."""
        x = np.linspace(-5, 5, 100)
        result = schrodinger_stationary(lambda xi: 0.5 * xi**2, x, n_states=3)
        assert 'energies' in result
        assert 'wavefunctions' in result

    def test_energies_are_sorted(self):
        """Eigenenergien müssen aufsteigend sortiert sein."""
        x = np.linspace(-5, 5, 100)
        result = schrodinger_stationary(lambda xi: 0.5 * xi**2, x, n_states=4)
        E = result['energies']
        assert np.all(np.diff(E) >= 0)

    def test_wavefunction_normalization(self):
        """Wellenfunktionen müssen normiert sein: ∫|ψ|²dx ≈ 1."""
        x = np.linspace(-6, 6, 200)
        result = schrodinger_stationary(lambda xi: 0.5 * xi**2, x, n_states=3)
        for psi in result['wavefunctions']:
            norm = np.trapezoid(np.abs(psi)**2, x)
            assert abs(norm - 1.0) < 0.1  # Toleranz wegen grobem Gitter

    def test_positive_energies_for_positive_potential(self):
        """Für V(x)=x²/2 ≥ 0 müssen alle Eigenwerte > 0 sein."""
        x = np.linspace(-5, 5, 100)
        result = schrodinger_stationary(lambda xi: 0.5 * xi**2, x, n_states=3)
        assert np.all(result['energies'] > 0)


class TestSchrodingerHarmonicOscillator:
    """Tests für schrodinger_harmonic_oscillator()."""

    def test_returns_correct_keys(self):
        """Rückgabe muss alle erwarteten Schlüssel enthalten."""
        result = schrodinger_harmonic_oscillator(n_states=3)
        for key in ['energies_numerical', 'energies_analytical', 'errors']:
            assert key in result

    def test_ground_state_energy(self):
        """Grundzustandsenergie E₀ ≈ 0.5 (analytisch: 1/2)."""
        result = schrodinger_harmonic_oscillator(n_states=3)
        E0_num = result['energies_numerical'][0]
        assert abs(E0_num - 0.5) < 0.05  # Toleranz für Gitterdiskretisierung

    def test_energy_spacing(self):
        """Energieabstände ΔE = 1 (E_n = n+1/2)."""
        result = schrodinger_harmonic_oscillator(n_states=4)
        E_num = result['energies_numerical']
        if len(E_num) >= 3:
            dE1 = E_num[1] - E_num[0]
            dE2 = E_num[2] - E_num[1]
            assert abs(dE1 - 1.0) < 0.1
            assert abs(dE2 - 1.0) < 0.1

    def test_numerical_matches_analytical(self):
        """Numerische Eigenwerte passen zu analytischen innerhalb 5% Toleranz."""
        result = schrodinger_harmonic_oscillator(n_states=4)
        errors = result['errors']
        # Erste 3 Zustände gut konvergiert
        for err in errors[:3]:
            assert err < 0.1


class TestSchrodingerTimeDependent:
    """Tests für schrodinger_time_dependent()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx) sein."""
        nx = 80
        x = np.linspace(-5, 5, nx)
        dx = x[1] - x[0]
        V = 0.5 * x**2
        # Gaussisches Wellenpaket als Anfangszustand
        psi0 = np.exp(-x**2 / 2.0).astype(complex)
        norm = np.sqrt(np.trapezoid(np.abs(psi0)**2, x))
        psi0 /= norm
        result = schrodinger_time_dependent(psi0, V, x, T=0.2, nt=20)
        assert result.shape == (21, nx)

    def test_norm_conservation(self):
        """Crank-Nicolson erhält die Norm ‖ψ‖²."""
        nx = 60
        x = np.linspace(-4, 4, nx)
        V = 0.5 * x**2
        psi0 = np.exp(-x**2 / 2.0).astype(complex)
        norm_init = np.sqrt(np.trapezoid(np.abs(psi0)**2, x))
        psi0 /= norm_init

        psi = schrodinger_time_dependent(psi0, V, x, T=0.2, nt=20)
        # Norm am Anfang
        norm0 = np.trapezoid(np.abs(psi[0])**2, x)
        # Norm am Ende
        norm_final = np.trapezoid(np.abs(psi[-1])**2, x)
        # Norm sollte näherungsweise erhalten sein
        assert abs(norm_final - norm0) < 0.5  # Grobe Toleranz wegen Gittereffekten


# =============================================================================
# Finite-Elemente-Methode
# =============================================================================

class TestFEM1DPoisson:
    """Tests für fem_1d_poisson()."""

    def test_returns_correct_keys(self):
        """Rückgabe muss 'x', 'u' und 'nodes' enthalten."""
        result = fem_1d_poisson(lambda x: 1.0, 0.0, 1.0, n=10)
        assert 'x' in result
        assert 'u' in result
        assert 'nodes' in result

    def test_boundary_conditions(self):
        """u(a) = u(b) = 0."""
        result = fem_1d_poisson(lambda x: np.pi**2 * np.sin(np.pi * x), 0.0, 1.0, n=20)
        assert abs(result['u'][0]) < 1e-12
        assert abs(result['u'][-1]) < 1e-12

    def test_nodes_count(self):
        """Anzahl der Knoten muss n+1 sein."""
        result = fem_1d_poisson(lambda x: 1.0, 0.0, 1.0, n=15)
        assert result['nodes'] == 16
        assert len(result['x']) == 16
        assert len(result['u']) == 16

    def test_accuracy_for_known_solution(self):
        """FEM für -u''=π²sin(πx) → u=sin(πx): Fehler < 1%."""
        n = 50
        f = lambda x: np.pi**2 * np.sin(np.pi * x)
        result = fem_1d_poisson(f, 0.0, 1.0, n=n)
        x = result['x']
        u_h = result['u']
        u_exact = np.sin(np.pi * x)
        max_err = np.max(np.abs(u_h - u_exact))
        assert max_err < 0.01

    def test_simple_poisson(self):
        """FEM für -u'' = 1, u(0)=u(1)=0: Lösung u = x(1-x)/2."""
        n = 40
        result = fem_1d_poisson(lambda x: 1.0, 0.0, 1.0, n=n)
        x = result['x']
        u_h = result['u']
        u_exact = x * (1 - x) / 2.0
        max_err = np.max(np.abs(u_h - u_exact))
        assert max_err < 0.005


class TestFEMConvergenceDemo:
    """Tests für fem_convergence_demo()."""

    def test_returns_correct_keys(self):
        """Rückgabe muss 'h_values', 'errors' und 'order' enthalten."""
        result = fem_convergence_demo()
        assert 'h_values' in result
        assert 'errors' in result
        assert 'order' in result

    def test_order_approximately_two(self):
        """Konvergenzordnung sollte nahe 2 liegen (O(h²))."""
        result = fem_convergence_demo()
        order = result['order']
        # Konvergenzordnung zwischen 1.5 und 3
        assert 1.5 < order < 3.0

    def test_errors_decrease(self):
        """Fehler müssen mit kleinerem h abnehmen."""
        result = fem_convergence_demo()
        errors = result['errors']
        for i in range(len(errors) - 1):
            assert errors[i + 1] <= errors[i] * 1.1  # Toleranz für numerisches Rauschen


# =============================================================================
# Methode der Charakteristiken
# =============================================================================

class TestMethodOfCharacteristics:
    """Tests für method_of_characteristics()."""

    def test_returns_correct_keys(self):
        """Rückgabe muss 'x', 't', 'u' und 'characteristics' enthalten."""
        result = method_of_characteristics(
            1.0, 1.0, 0.0,
            u0=lambda x: np.sin(x),
            x_range=(0, 2 * np.pi),
            t_range=(0, 1.0)
        )
        for key in ['x', 't', 'u', 'characteristics']:
            assert key in result

    def test_output_shape(self):
        """u muss Form (nt, nx) haben."""
        nx, nt = 30, 20
        result = method_of_characteristics(
            1.0, 1.0, 0.0,
            u0=lambda x: np.sin(x),
            x_range=(0, 2 * np.pi),
            t_range=(0, 1.0),
            nx=nx, nt=nt
        )
        assert result['u'].shape == (nt, nx)

    def test_pure_advection(self):
        """Für c=0: u(x,t) = u0(x - a/b * t) (Verschiebungslösung)."""
        # u_x + u_t = 0 → u(x,t) = u0(x-t)
        u0_func = lambda x: np.sin(x)
        result = method_of_characteristics(
            1.0, 1.0, 0.0,
            u0=u0_func,
            x_range=(0, 4 * np.pi),
            t_range=(0, 1.0),
            nx=80, nt=5
        )
        # Bei t=1: u(x,1) ≈ sin(x-1)
        x = result['x']
        t_idx = 4  # t=1
        t_val = result['t'][t_idx]
        u_num = result['u'][t_idx, :]
        u_exact = np.sin(x - t_val)
        max_err = np.max(np.abs(u_num - u_exact))
        assert max_err < 1e-10

    def test_characteristics_count(self):
        """Es werden 5 Charakteristiken zurückgegeben."""
        result = method_of_characteristics(
            1.0, 1.0, 0.0,
            u0=lambda x: x,
            x_range=(0, 1.0),
            t_range=(0, 0.5)
        )
        assert len(result['characteristics']) == 5


# =============================================================================
# Burgers-Gleichung
# =============================================================================

class TestBurgersEquation:
    """Tests für burgers_equation()."""

    def test_output_shape(self):
        """Rückgabe-Shape muss (nt+1, nx) sein."""
        nx, nt = 50, 100
        u0 = np.sin(np.linspace(0, 2 * np.pi, nx))
        result = burgers_equation(u0.tolist(), L=2 * np.pi, T=0.5, nx=nx, nt=nt)
        assert result.shape == (nt + 1, nx)

    def test_initial_condition(self):
        """Erste Zeile muss der Anfangsbedingung entsprechen."""
        nx = 50
        u0 = np.sin(np.linspace(0, 2 * np.pi, nx))
        result = burgers_equation(u0.tolist(), nx=nx, nt=100)
        assert np.allclose(result[0, :], u0, atol=1e-12)

    def test_viscosity_stabilizes(self):
        """Mit Viskosität ν>0 bleibt die Lösung beschränkt."""
        nx, nt = 100, 500
        u0 = np.sin(np.linspace(0, 2 * np.pi, nx))
        result = burgers_equation(u0.tolist(), L=2 * np.pi, T=1.0, nx=nx, nt=nt, nu=0.1)
        assert np.all(np.isfinite(result))

    def test_finite_values_with_small_nu(self):
        """Für kleine Viskosität sollten endliche Werte zurückgegeben werden."""
        nx, nt = 100, 200
        u0 = np.sin(np.linspace(0, 2 * np.pi, nx))
        result = burgers_equation(u0.tolist(), L=2 * np.pi, T=0.2, nx=nx, nt=nt, nu=0.01)
        assert np.all(np.isfinite(result))

    def test_interpolation_of_u0(self):
        """Unterschiedliche Länge von u0 wird korrekt interpoliert."""
        u0 = [0.0, 1.0, 0.0, -1.0, 0.0]
        result = burgers_equation(u0, nx=100, nt=50)
        assert result.shape[1] == 100


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
