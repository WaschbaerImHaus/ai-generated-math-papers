"""
@file test_differential_geometry.py
@brief Tests für das Differentialgeometrie-Modul (differential_geometry.py).
@description
    Umfassende Tests für:
    - ParametricCurve: Frenet-Serret-Rahmen, Krümmung, Torsion, Bogenlänge
    - ParametricSurface: Fundamentalformen, Gaußsche/Mittlere Krümmung
    - GeodesicComputation: Geodätische auf Sphäre, allgemeine Geodäten
    - RiemannianGeometry: Schnittkrümmung, Ricci-Tensor, Gauß-Bonnet
    - ClassicalSurfaces: Sphäre, Torus, Zylinder, Helikoid, Katenoid, Sattel
    - CurveTheory2D: Krümmung, Evolute, isoperimetrische Ungleichung

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import pytest
import numpy as np
import sympy as sp

from differential_geometry import (
    ParametricCurve,
    ParametricSurface,
    GeodesicComputation,
    RiemannianGeometry,
    ClassicalSurfaces,
    CurveTheory2D,
)


# ============================================================================
# Hilfsfunktionen für Tests
# ============================================================================

def assert_close(a: float, b: float, tol: float = 1e-4, msg: str = "") -> None:
    """Prüft |a - b| < tol."""
    assert abs(a - b) < tol, f"{msg}: |{a} - {b}| = {abs(a-b)} >= {tol}"


# ============================================================================
# Fixtures: Häufig verwendete Kurven
# ============================================================================

@pytest.fixture
def helix_curve():
    """Helix r(t) = (cos t, sin t, t). κ = 1/2, τ = 1/2."""
    t = sp.Symbol('t')
    return ParametricCurve([sp.cos(t), sp.sin(t), t], t, (0, 2 * math.pi))


@pytest.fixture
def circle_curve():
    """Einheitskreis r(t) = (cos t, sin t, 0). κ = 1, τ = 0."""
    t = sp.Symbol('t')
    return ParametricCurve([sp.cos(t), sp.sin(t), sp.Integer(0)], t, (0, 2 * math.pi))


@pytest.fixture
def line_curve():
    """Gerade r(t) = (t, t, t). κ = 0, τ = 0."""
    t = sp.Symbol('t')
    return ParametricCurve([t, t, t], t, (0, 1))


@pytest.fixture
def parabola_curve():
    """Parabel r(t) = (t, t², 0). 2D-Kurve."""
    t = sp.Symbol('t')
    return ParametricCurve([t, t**2, sp.Integer(0)], t, (-1, 1))


# ============================================================================
# Tests: ParametricCurve – Helix
# ============================================================================

class TestParametricCurveHelix:
    """Tests für die Helix r(t) = (cos t, sin t, t)."""

    def test_velocity_helix(self, helix_curve):
        """Geschwindigkeitsvektor der Helix ist (-sin t, cos t, 1)."""
        v = helix_curve.velocity(0.0)
        # Bei t=0: v = (0, 1, 1)
        assert_close(v[0], 0.0, tol=1e-8, msg="v_x bei t=0")
        assert_close(v[1], 1.0, tol=1e-8, msg="v_y bei t=0")
        assert_close(v[2], 1.0, tol=1e-8, msg="v_z bei t=0")

    def test_speed_helix(self, helix_curve):
        """Geschwindigkeit der Helix ist konstant |r'| = √2."""
        for t_val in [0.0, 0.5, 1.0, math.pi]:
            s = helix_curve.speed(t_val)
            assert_close(s, math.sqrt(2), tol=1e-8, msg=f"speed bei t={t_val}")

    def test_arc_length_helix(self, helix_curve):
        """Bogenlänge der Helix von 0 bis 2π ist 2π√2."""
        L = helix_curve.arc_length(0, 2 * math.pi)
        expected = 2 * math.pi * math.sqrt(2)
        assert_close(L, expected, tol=1e-5, msg="Bogenlänge Helix")

    def test_curvature_helix(self, helix_curve):
        """Krümmung der Helix ist konstant κ = 1/2."""
        for t_val in [0.0, 0.5, 1.0, math.pi]:
            kappa = helix_curve.curvature(t_val)
            assert_close(kappa, 0.5, tol=1e-4, msg=f"Krümmung Helix bei t={t_val}")

    def test_torsion_helix(self, helix_curve):
        """Torsion der Helix ist konstant τ = 1/2."""
        for t_val in [0.0, 0.5, 1.0]:
            tau = helix_curve.torsion(t_val)
            assert_close(tau, 0.5, tol=1e-3, msg=f"Torsion Helix bei t={t_val}")

    def test_unit_tangent_helix(self, helix_curve):
        """Tangenteneinheitsvektor der Helix hat Länge 1."""
        for t_val in [0.0, math.pi / 4, math.pi]:
            T = helix_curve.unit_tangent(t_val)
            assert_close(_norm(T), 1.0, tol=1e-8, msg=f"|T| bei t={t_val}")

    def test_principal_normal_helix(self, helix_curve):
        """Hauptnormalenvektor der Helix hat Länge 1."""
        N = helix_curve.principal_normal(1.0)
        assert_close(_norm(N), 1.0, tol=1e-5, msg="|N| Helix")

    def test_binormal_helix(self, helix_curve):
        """Binormalvektor der Helix hat Länge 1."""
        B = helix_curve.binormal(1.0)
        assert_close(_norm(B), 1.0, tol=1e-5, msg="|B| Helix")

    def test_frenet_frame_orthonormality(self, helix_curve):
        """Frenet-Serret-Rahmen {T, N, B} ist orthonormal."""
        frame = helix_curve.frenet_serret_frame(1.0)
        T, N, B = frame["T"], frame["N"], frame["B"]
        # Orthonormalität
        assert_close(np.dot(T, T), 1.0, tol=1e-5, msg="T·T=1")
        assert_close(np.dot(N, N), 1.0, tol=1e-5, msg="N·N=1")
        assert_close(np.dot(B, B), 1.0, tol=1e-5, msg="B·B=1")
        assert_close(np.dot(T, N), 0.0, tol=1e-5, msg="T·N=0")
        assert_close(np.dot(T, B), 0.0, tol=1e-5, msg="T·B=0")
        assert_close(np.dot(N, B), 0.0, tol=1e-5, msg="N·B=0")

    def test_frenet_frame_kappa_tau(self, helix_curve):
        """Frenet-Serret-Rahmen enthält korrekte κ und τ."""
        frame = helix_curve.frenet_serret_frame(1.0)
        assert_close(frame["kappa"], 0.5, tol=1e-4, msg="κ im Frame")
        assert_close(frame["tau"], 0.5, tol=1e-3, msg="τ im Frame")

    def test_frenet_serret_formulas_return(self, helix_curve):
        """frenet_serret_formulas() gibt korrektes Dictionary zurück."""
        formulas = helix_curve.frenet_serret_formulas()
        assert "T_prime" in formulas
        assert "N_prime" in formulas
        assert "B_prime" in formulas
        assert "κ" in formulas["T_prime"] or "kappa" in formulas["T_prime"].lower() or "κ" in formulas["T_prime"]


# ============================================================================
# Tests: ParametricCurve – Einheitskreis
# ============================================================================

class TestParametricCurveCircle:
    """Tests für den Einheitskreis r(t) = (cos t, sin t)."""

    def test_curvature_circle(self, circle_curve):
        """Krümmung des Einheitskreises ist κ = 1."""
        for t_val in [0.0, math.pi / 4, math.pi / 2, math.pi]:
            kappa = circle_curve.curvature(t_val)
            assert_close(kappa, 1.0, tol=1e-4, msg=f"κ Kreis bei t={t_val}")

    def test_torsion_circle(self, circle_curve):
        """Torsion des (ebenen) Kreises ist τ = 0."""
        for t_val in [0.0, 0.5, 1.0]:
            tau = circle_curve.torsion(t_val)
            assert_close(tau, 0.0, tol=1e-3, msg=f"τ Kreis bei t={t_val}")

    def test_arc_length_circle(self, circle_curve):
        """Bogenlänge des Einheitskreises von 0 bis 2π ist 2π."""
        L = circle_curve.arc_length(0, 2 * math.pi)
        assert_close(L, 2 * math.pi, tol=1e-5, msg="Bogenlänge Kreis")

    def test_speed_circle(self, circle_curve):
        """Geschwindigkeit des Einheitskreises ist |r'| = 1."""
        for t_val in [0.0, math.pi / 4, math.pi]:
            s = circle_curve.speed(t_val)
            assert_close(s, 1.0, tol=1e-8, msg=f"speed Kreis bei t={t_val}")

    def test_unit_tangent_circle(self, circle_curve):
        """T(0) = (0, 1, 0) für den Einheitskreis."""
        T = circle_curve.unit_tangent(0.0)
        assert_close(T[0], 0.0, tol=1e-6, msg="T_x")
        assert_close(T[1], 1.0, tol=1e-6, msg="T_y")


# ============================================================================
# Tests: ParametricCurve – Gerade (κ = 0)
# ============================================================================

class TestParametricCurveLine:
    """Edge Case: Gerade mit κ = 0."""

    def test_curvature_line_zero(self, line_curve):
        """Krümmung einer Geraden ist κ = 0."""
        for t_val in [0.0, 0.5, 1.0]:
            kappa = line_curve.curvature(t_val)
            assert_close(kappa, 0.0, tol=1e-8, msg=f"κ=0 Gerade bei t={t_val}")

    def test_torsion_line_zero(self, line_curve):
        """Torsion einer Geraden ist τ = 0."""
        tau = line_curve.torsion(0.5)
        assert_close(tau, 0.0, tol=1e-6, msg="τ=0 Gerade")

    def test_arc_length_line(self, line_curve):
        """Bogenlänge der Geraden (t,t,t) von 0 bis 1 ist √3."""
        L = line_curve.arc_length(0, 1)
        assert_close(L, math.sqrt(3), tol=1e-5, msg="Bogenlänge Gerade")


# ============================================================================
# Tests: ParametricCurve – 2D-Kurven
# ============================================================================

class TestParametricCurve2D:
    """Tests für 2D-Kurven (dim=2 oder 3 mit z=0)."""

    def test_curvature_parabola_at_vertex(self, parabola_curve):
        """Krümmung der Parabel y=t² am Scheitel t=0 ist κ = 2."""
        kappa = parabola_curve.curvature(0.0)
        assert_close(kappa, 2.0, tol=1e-3, msg="κ Parabel an t=0")

    def test_curvature_parabola_positive(self, parabola_curve):
        """Krümmung der Parabel ist positiv."""
        for t_val in [-0.5, 0.0, 0.5]:
            kappa = parabola_curve.curvature(t_val)
            assert kappa >= 0, f"κ < 0 bei t={t_val}"


# ============================================================================
# Tests: ParametricSurface – Sphäre
# ============================================================================

class TestParametricSurfaceSphere:
    """Tests für die Einheitssphäre."""

    @pytest.fixture
    def sphere(self):
        """Einheitssphäre als ParametricSurface."""
        def r(theta, phi):
            return np.array([
                math.sin(theta) * math.cos(phi),
                math.sin(theta) * math.sin(phi),
                math.cos(theta),
            ])
        return ParametricSurface(r, (0, math.pi), (0, 2 * math.pi))

    def test_unit_normal_on_sphere(self, sphere):
        """Einheitsnormale auf der Sphäre zeigt vom Ursprung weg."""
        n = sphere.unit_normal(math.pi / 2, 0.0)
        # Bei (θ=π/2, φ=0): Punkt ist (1,0,0), Normale ist auch (1,0,0)
        assert_close(_norm(n), 1.0, tol=1e-5, msg="|n|=1 Sphäre")
        assert_close(n[0], 1.0, tol=1e-4, msg="n_x Sphäre Äquator")

    def test_first_fundamental_form_sphere(self, sphere):
        """Erste Fundamentalform der Einheitssphäre: E=1, F=0, G=sin²θ."""
        theta = math.pi / 2
        phi = 0.0
        E, F, G = sphere.first_fundamental_form(theta, phi)
        assert_close(E, 1.0, tol=1e-4, msg="E Sphäre")
        assert_close(F, 0.0, tol=1e-4, msg="F Sphäre")
        assert_close(G, 1.0, tol=1e-4, msg="G=sin²(π/2)=1")

    def test_gaussian_curvature_sphere(self, sphere):
        """Gaußsche Krümmung der Einheitssphäre ist K = 1."""
        for theta in [math.pi / 4, math.pi / 2, 3 * math.pi / 4]:
            K = sphere.gaussian_curvature(theta, 0.0)
            assert_close(K, 1.0, tol=0.05, msg=f"K=1 Sphäre bei θ={theta:.2f}")

    def test_mean_curvature_sphere(self, sphere):
        """Betrag der mittleren Krümmung der Einheitssphäre ist |H| = 1.

        Das Vorzeichen hängt von der Orientierung der Normalen ab:
        Innenorientierung → H = -1, Außenorientierung → H = +1.
        """
        H = sphere.mean_curvature(math.pi / 2, 0.0)
        assert_close(abs(H), 1.0, tol=0.1, msg="|H|=1 Einheitssphäre")

    def test_principal_curvatures_sphere(self, sphere):
        """Hauptkrümmungen der Einheitssphäre: |κ₁| = |κ₂| = 1.

        Das Vorzeichen hängt von der Normalenorientierung ab.
        """
        k1, k2 = sphere.principal_curvatures(math.pi / 2, 0.0)
        assert_close(abs(k1), 1.0, tol=0.1, msg="|κ₁|=1 Sphäre")
        assert_close(abs(k2), 1.0, tol=0.1, msg="|κ₂|=1 Sphäre")

    def test_area_element_sphere_equator(self, sphere):
        """Flächenelement der Einheitssphäre bei θ=π/2 ist sin(π/2)=1."""
        dA = sphere.area_element(math.pi / 2, 0.0)
        assert_close(dA, 1.0, tol=1e-4, msg="dA Sphäre am Äquator")

    def test_total_area_sphere(self, sphere):
        """Gesamtfläche der Einheitssphäre ist 4π ≈ 12.566."""
        A = sphere.total_area(n_pts=25)
        assert_close(A, 4 * math.pi, tol=0.2, msg="Fläche Einheitssphäre ≈ 4π")


# ============================================================================
# Tests: ParametricSurface – Torus
# ============================================================================

class TestParametricSurfaceTorus:
    """Tests für den Torus (R=2, r=1)."""

    @pytest.fixture
    def torus(self):
        """Torus mit R=2, r=1."""
        R, r = 2.0, 1.0
        def r_torus(u, v):
            return np.array([
                (R + r * math.cos(v)) * math.cos(u),
                (R + r * math.cos(v)) * math.sin(u),
                r * math.sin(v),
            ])
        return ParametricSurface(r_torus, (0, 2 * math.pi), (0, 2 * math.pi))

    def test_gaussian_curvature_torus_outer(self, torus):
        """Gaußsche Krümmung des Torus an der Außenseite (v=0) ist K > 0."""
        K_outer = torus.gaussian_curvature(0.0, 0.0)
        assert K_outer > 0, f"K={K_outer} sollte >0 an der Außenseite sein"

    def test_gaussian_curvature_torus_inner(self, torus):
        """Gaußsche Krümmung des Torus an der Innenseite (v=π) ist K < 0."""
        K_inner = torus.gaussian_curvature(0.0, math.pi)
        assert K_inner < 0, f"K={K_inner} sollte <0 an der Innenseite sein"

    def test_gaussian_curvature_torus_side(self, torus):
        """Gaußsche Krümmung des Torus bei v=π/2 ist K = 0."""
        K_side = torus.gaussian_curvature(0.0, math.pi / 2)
        assert_close(K_side, 0.0, tol=0.1, msg="K≈0 an der Seite des Torus")

    def test_mean_curvature_torus_nonzero(self, torus):
        """Mittlere Krümmung des Torus ist im Allgemeinen ≠ 0."""
        H = torus.mean_curvature(0.0, 0.0)
        # H > 0 an der Außenseite
        assert H != 0.0, "H sollte nicht 0 sein"

    def test_first_fundamental_form_torus_F_zero(self, torus):
        """F=0 für orthogonale Koordinaten des Torus."""
        E, F, G = torus.first_fundamental_form(0.0, 0.0)
        assert_close(F, 0.0, tol=1e-4, msg="F=0 Torus")

    def test_total_area_torus(self, torus):
        """Gesamtfläche des Torus (R=2, r=1) ist 4π²Rr = 4π²·2·1 ≈ 78.957."""
        A = torus.total_area(n_pts=20)
        expected = 4 * math.pi ** 2 * 2.0 * 1.0
        assert_close(A, expected, tol=2.0, msg="Fläche Torus")


# ============================================================================
# Tests: ParametricSurface – Zylinder (K = 0)
# ============================================================================

class TestParametricSurfaceCylinder:
    """Tests für den Zylinder (K = 0, lokal flach)."""

    @pytest.fixture
    def cylinder(self):
        """Einheitszylinder."""
        def r_cyl(u, v):
            return np.array([math.cos(u), math.sin(u), v])
        return ParametricSurface(r_cyl, (0, 2 * math.pi), (0, 2.0))

    def test_gaussian_curvature_cylinder_zero(self, cylinder):
        """Gaußsche Krümmung des Zylinders ist K = 0."""
        for u in [0.0, math.pi / 4, math.pi / 2, math.pi]:
            for v in [0.5, 1.0]:
                K = cylinder.gaussian_curvature(u, v)
                assert_close(K, 0.0, tol=1e-3, msg=f"K=0 Zylinder bei u={u:.2f}")


# ============================================================================
# Tests: ClassicalSurfaces
# ============================================================================

class TestClassicalSurfaces:
    """Tests für ClassicalSurfaces."""

    @pytest.fixture
    def cs(self):
        return ClassicalSurfaces()

    def test_sphere_returns_surface(self, cs):
        """sphere() gibt ParametricSurface zurück."""
        S = cs.sphere(radius=1.0)
        assert isinstance(S, ParametricSurface)

    def test_sphere_gaussian_curvature(self, cs):
        """Sphäre mit R=2: K = 1/R² = 0.25."""
        S = cs.sphere(radius=2.0)
        K = S.gaussian_curvature(math.pi / 2, 0.0)
        assert_close(K, 0.25, tol=0.05, msg="K=1/4 Sphäre R=2")

    def test_torus_returns_surface(self, cs):
        """torus() gibt ParametricSurface zurück."""
        T = cs.torus(R=2.0, r=1.0)
        assert isinstance(T, ParametricSurface)

    def test_cylinder_gaussian_zero(self, cs):
        """Zylinder hat K = 0."""
        C = cs.cylinder(radius=1.0)
        K = C.gaussian_curvature(0.5, 0.5)
        assert_close(K, 0.0, tol=1e-3, msg="K=0 Zylinder")

    def test_helicoid_returns_surface(self, cs):
        """helicoid() gibt ParametricSurface zurück."""
        H = cs.helicoid()
        assert isinstance(H, ParametricSurface)

    def test_helicoid_mean_curvature_zero(self, cs):
        """Helikoid ist Minimalfläche: H ≈ 0."""
        H_surf = cs.helicoid()
        for u, v in [(0.0, 0.5), (1.0, 0.3), (-1.0, -0.3)]:
            H = H_surf.mean_curvature(u, v)
            assert_close(H, 0.0, tol=1e-3, msg=f"H=0 Helikoid bei u={u}, v={v}")

    def test_catenoid_returns_surface(self, cs):
        """catenoid() gibt ParametricSurface zurück."""
        C = cs.catenoid()
        assert isinstance(C, ParametricSurface)

    def test_catenoid_mean_curvature_zero(self, cs):
        """Katenoid ist Minimalfläche: H ≈ 0."""
        C = cs.catenoid()
        for u in [0.0, math.pi / 2, math.pi]:
            H = C.mean_curvature(u, 0.0)
            assert_close(H, 0.0, tol=1e-3, msg=f"H=0 Katenoid bei u={u:.2f}")

    def test_saddle_returns_surface(self, cs):
        """saddle_surface() gibt ParametricSurface zurück."""
        S = cs.saddle_surface()
        assert isinstance(S, ParametricSurface)

    def test_saddle_gaussian_curvature_negative(self, cs):
        """Sattelfläche hat K < 0."""
        S = cs.saddle_surface()
        K = S.gaussian_curvature(0.0, 0.0)
        assert K < 0, f"K={K} sollte negativ sein"

    def test_minimal_surface_check_helicoid(self, cs):
        """minimal_surface_check erkennt Helikoid als Minimalfläche."""
        H_surf = cs.helicoid()
        result = cs.minimal_surface_check(
            H_surf, (-2.0, 2.0), (-0.8, 0.8), n_samples=4, tol=1e-2
        )
        assert result["is_minimal"], f"Helikoid sollte Minimalfläche sein: {result}"

    def test_minimal_surface_check_catenoid(self, cs):
        """minimal_surface_check erkennt Katenoid als Minimalfläche."""
        C = cs.catenoid()
        result = cs.minimal_surface_check(
            C, (0.1, 2 * math.pi - 0.1), (-1.5, 1.5), n_samples=4, tol=1e-2
        )
        assert result["is_minimal"], f"Katenoid sollte Minimalfläche sein: {result}"

    def test_minimal_surface_check_sphere_not_minimal(self, cs):
        """minimal_surface_check erkennt Sphäre als KEINE Minimalfläche."""
        S = cs.sphere(radius=1.0)
        result = cs.minimal_surface_check(
            S, (0.2, math.pi - 0.2), (0.1, 2 * math.pi - 0.1),
            n_samples=4, tol=0.1
        )
        # Die Sphäre hat H ≈ 1 ≠ 0
        assert not result["is_minimal"], "Sphäre ist keine Minimalfläche"

    def test_cylinder_not_minimal(self, cs):
        """Zylinder ist keine Minimalfläche (H ≠ 0)."""
        C = cs.cylinder(radius=1.0)
        result = cs.minimal_surface_check(
            C, (0.1, 2 * math.pi - 0.1), (0.1, 1.9), n_samples=4, tol=0.1
        )
        assert not result["is_minimal"], "Zylinder ist keine Minimalfläche"


# ============================================================================
# Tests: GeodesicComputation
# ============================================================================

class TestGeodesicComputation:
    """Tests für geodätische Berechnungen."""

    @pytest.fixture
    def gc(self):
        return GeodesicComputation()

    def test_geodesic_sphere_runs(self, gc):
        """geodesic_sphere() terminiert erfolgreich."""
        result = gc.geodesic_sphere(math.pi / 2, 0.0, [1.0, 0.0], t_max=math.pi)
        assert result["success"], "Geodätenintegration sollte erfolgreich sein"

    def test_geodesic_sphere_returns_arrays(self, gc):
        """geodesic_sphere() gibt theta, phi, t Arrays zurück."""
        result = gc.geodesic_sphere(math.pi / 2, 0.0, [1.0, 0.0], t_max=1.0)
        assert "theta" in result
        assert "phi" in result
        assert "t" in result
        assert len(result["theta"]) > 0

    def test_geodesic_sphere_great_circle(self, gc):
        """Geodätische auf S² bleibt auf dem Äquator (θ'=0 für Äquatordirektion)."""
        # Startpunkt am Äquator, Bewegung entlang des Äquators (nur φ ändert sich)
        result = gc.geodesic_sphere(math.pi / 2, 0.0, [0.0, 1.0], t_max=1.0)
        theta_values = result["theta"]
        # θ sollte nahe π/2 bleiben
        for theta in theta_values:
            assert abs(theta - math.pi / 2) < 0.1, f"θ={theta} weicht zu weit von π/2 ab"

    def test_geodesic_general_flat_space(self, gc):
        """Geodätische im flachen Raum ist eine Gerade."""
        def flat_metric(x):
            return np.eye(2)

        result = gc.geodesic_equations_general(
            flat_metric,
            initial_pos=[0.0, 0.0],
            initial_vel=[1.0, 0.0],
            t_max=1.0,
        )
        assert result["success"]
        # Endposition sollte nahe (1, 0) sein
        pos = result["positions"]
        final_pos = pos[-1]
        assert_close(final_pos[0], 1.0, tol=1e-4, msg="x-Koordinate Geodätische flach")
        assert_close(final_pos[1], 0.0, tol=1e-4, msg="y-Koordinate Geodätische flach")

    def test_conjugate_points_sphere(self, gc):
        """conjugate_points_demo('sphere') gibt korrekte Struktur zurück."""
        result = gc.conjugate_points_demo("sphere")
        assert "description" in result
        assert "focal_point_distance" in result
        assert_close(result["focal_point_distance"], math.pi, tol=1e-10,
                     msg="Konjugierter Abstand auf S² ist π")

    def test_conjugate_points_plane(self, gc):
        """Ebene hat keine konjugierten Punkte."""
        result = gc.conjugate_points_demo("plane")
        assert "description" in result
        assert "keine konjugierten" in result["description"].lower() or \
               "keine" in result["description"].lower() or \
               "K=0" in result["description"]

    def test_exponential_map_surface(self, gc):
        """exponential_map_surface() gibt Punkt auf Einheitssphäre zurück."""
        p = np.array([1.0, 0.0, 0.0])
        v = np.array([0.0, 0.1, 0.0])  # kleiner Tangentialvektor
        q = gc.exponential_map_surface(p, v, t=0.5)
        # Auf Einheitssphäre: |q| = 1
        assert_close(_norm(q), 1.0, tol=0.1, msg="|exp_p(v)| ≈ 1")


# ============================================================================
# Tests: RiemannianGeometry
# ============================================================================

class TestRiemannianGeometry:
    """Tests für Riemannsche Geometrie."""

    @pytest.fixture
    def rg(self):
        return RiemannianGeometry()

    def test_ricci_flat_space(self, rg):
        """Ricci-Tensor im flachen Raum ist null."""
        R_flat = np.zeros((3, 3, 3, 3))
        metric = np.eye(3)
        Ric = rg.ricci_curvature(R_flat, metric)
        assert np.allclose(Ric, 0.0), "Ricci-Tensor flach = 0"

    def test_scalar_curvature_flat(self, rg):
        """Skalarkrümmung im flachen Raum ist 0."""
        R_flat = np.zeros((3, 3, 3, 3))
        metric = np.eye(3)
        Ric = rg.ricci_curvature(R_flat, metric)
        R_scalar = rg.scalar_curvature(Ric, metric)
        assert_close(R_scalar, 0.0, tol=1e-10, msg="Skalarkrümmung flach=0")

    def test_sectional_curvature_flat(self, rg):
        """Schnittkrümmung im flachen Raum ist 0."""
        R_flat = np.zeros((3, 3, 3, 3))
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        K = rg.sectional_curvature(R_flat, v1, v2)
        assert_close(K, 0.0, tol=1e-10, msg="K flach=0")

    def test_sectional_curvature_parallel_vectors(self, rg):
        """Schnittkrümmung für parallele Vektoren ist nicht definiert (0 zurückgeben)."""
        R_flat = np.zeros((3, 3, 3, 3))
        v = np.array([1.0, 0.0, 0.0])
        K = rg.sectional_curvature(R_flat, v, v)
        assert_close(K, 0.0, tol=1e-10, msg="K für parallele Vektoren = 0")

    def test_covariant_derivative_demo_structure(self, rg):
        """covariant_derivative_demo() gibt korrektes Dictionary zurück."""
        def flat(x): return np.eye(3)
        V = np.array([1.0, 0.0, 0.0])
        result = rg.covariant_derivative_demo(flat, V)
        assert "formula" in result
        assert "∇" in result["formula"]

    def test_parallel_transport_demo_structure(self, rg):
        """parallel_transport_demo() gibt korrektes Dictionary zurück."""
        curve = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
        vector = np.array([0.0, 1.0, 0.0])
        metric = np.eye(3)
        result = rg.parallel_transport_demo(curve, vector, metric)
        assert "initial_vector" in result
        assert "transported_vector" in result
        assert "holonomy_formula" in result

    def test_gauss_bonnet_sphere(self, rg):
        """Gauß-Bonnet für Sphäre: χ=2, ∫K dA = 4π."""
        result = rg.gauss_bonnet_theorem_demo("sphere")
        assert result["euler_characteristic"] == 2
        assert_close(result["integral_K"], 4 * math.pi, tol=1e-10,
                     msg="∫K dA = 4π für Sphäre")
        assert result["satisfied"]

    def test_gauss_bonnet_torus(self, rg):
        """Gauß-Bonnet für Torus: χ=0, ∫K dA = 0."""
        result = rg.gauss_bonnet_theorem_demo("torus")
        assert result["euler_characteristic"] == 0
        assert_close(result["integral_K"], 0.0, tol=1e-10,
                     msg="∫K dA = 0 für Torus")
        assert result["satisfied"]

    def test_gauss_bonnet_plane(self, rg):
        """Gauß-Bonnet für Ebene: K=0."""
        result = rg.gauss_bonnet_theorem_demo("plane")
        assert result["satisfied"]
        assert_close(result["K_constant"], 0.0, tol=1e-10, msg="K=0 Ebene")

    def test_gauss_bonnet_consistency(self, rg):
        """Gauß-Bonnet: ∫K dA = 2πχ für alle Flächen."""
        for surface in ["sphere", "torus", "plane"]:
            result = rg.gauss_bonnet_theorem_demo(surface)
            assert result["satisfied"], f"Gauß-Bonnet verletzt für {surface}"
            if result["euler_characteristic"] is not None:
                expected = 2 * math.pi * result["euler_characteristic"]
                assert_close(
                    result["integral_K"], expected, tol=1e-8,
                    msg=f"∫K dA = 2πχ für {surface}"
                )


# ============================================================================
# Tests: CurveTheory2D
# ============================================================================

class TestCurveTheory2D:
    """Tests für ebene Kurventheorie."""

    @pytest.fixture
    def ct(self):
        return CurveTheory2D()

    @pytest.fixture
    def x(self):
        return sp.Symbol('x')

    # --- curvature_from_cartesian ---

    def test_curvature_parabola(self, ct, x):
        """Krümmung y=x² an x=0 ist κ = 2."""
        kappa = ct.curvature_from_cartesian(x ** 2, x, 0.0)
        assert_close(kappa, 2.0, tol=1e-8, msg="κ=2 Parabel an x=0")

    def test_curvature_straight_line(self, ct, x):
        """Krümmung einer Geraden y=2x ist κ = 0."""
        kappa = ct.curvature_from_cartesian(2 * x, x, 0.0)
        assert_close(kappa, 0.0, tol=1e-10, msg="κ=0 Gerade")

    def test_curvature_circle_function(self, ct, x):
        """Krümmung des Einheitskreises y=√(1-x²) an x=0 ist κ = 1."""
        # y = sqrt(1 - x²), Krümmung am Hochpunkt (x=0) ist 1/R = 1
        y_expr = sp.sqrt(1 - x ** 2)
        kappa = ct.curvature_from_cartesian(y_expr, x, 0.0)
        assert_close(kappa, 1.0, tol=1e-6, msg="κ=1 Halbkreis an x=0")

    def test_curvature_positive(self, ct, x):
        """Krümmung ist immer ≥ 0."""
        for expr in [x ** 2, x ** 3, sp.sin(x)]:
            kappa = ct.curvature_from_cartesian(expr, x, 0.5)
            assert kappa >= 0, f"κ < 0 für {expr}"

    def test_curvature_sine(self, ct, x):
        """Krümmung von y=sin(x) an x=0 ist κ = 1 (da sin''(0)=-0, sin'(0)=1, Krümmung=0)."""
        # y=sin(x): f'=cos, f''=-sin; κ = |f''|/(1+f'²)^{3/2}
        # An x=0: f'=1, f''=0 → κ = 0
        kappa = ct.curvature_from_cartesian(sp.sin(x), x, 0.0)
        assert_close(kappa, 0.0, tol=1e-8, msg="κ=0 sin(x) an x=0")

    def test_curvature_sine_at_pi_half(self, ct, x):
        """Krümmung von y=sin(x) an x=π/2 ist κ = 1/(1+0)^{3/2} = 1."""
        kappa = ct.curvature_from_cartesian(sp.sin(x), x, math.pi / 2)
        assert_close(kappa, 1.0, tol=1e-6, msg="κ=1 sin(x) an x=π/2")

    # --- evolute ---

    def test_evolute_circle(self, ct):
        """Evolute des Einheitskreises r(t)=(cos t, sin t) ist der Ursprung (0,0)."""
        t = sp.Symbol('t')
        E = ct.evolute([sp.cos(t), sp.sin(t)], t)
        # Für den Einheitskreis: Krümmungsradius R=1, Evolute = Mittelpunkt = (0,0)
        # E_x = cos t - sin t · 1/κ · N_x etc.
        # Symbolisch prüfen bei t=0
        Ex_at_0 = float(E[0].subs(t, 0))
        Ey_at_0 = float(E[1].subs(t, 0))
        assert_close(Ex_at_0, 0.0, tol=1e-8, msg="Evolute Kreis E_x(0)=0")
        assert_close(Ey_at_0, 0.0, tol=1e-8, msg="Evolute Kreis E_y(0)=0")

    def test_evolute_returns_list(self, ct):
        """evolute() gibt eine Liste von 2 SymPy-Ausdrücken zurück."""
        t = sp.Symbol('t')
        E = ct.evolute([sp.cos(t), sp.sin(t)], t)
        assert len(E) == 2
        assert isinstance(E[0], sp.Expr)
        assert isinstance(E[1], sp.Expr)

    # --- involute ---

    def test_involute_returns_list(self, ct):
        """involute() gibt eine Liste von 2 Ausdrücken zurück."""
        t = sp.Symbol('t')
        I = ct.involute([sp.cos(t), sp.sin(t)], t, c=0.0)
        assert len(I) == 2

    def test_involute_line_from_circle(self, ct):
        """Evolvente des Kreises ist eine Spirale (nicht-triviale Kurve)."""
        t = sp.Symbol('t')
        I = ct.involute([sp.cos(t), sp.sin(t)], t, c=0.0)
        # Ausdruck sollte von t abhängen (nicht konstant)
        # Prüfen: Ausdruck ist kein konstanter Wert
        assert t in I[0].free_symbols or t in I[1].free_symbols

    # --- isoperimetric_inequality_check ---

    def test_isoperimetric_circle(self, ct):
        """Kreis erfüllt die isoperimetrische Gleichheit: Q = 1."""
        r = 1.0
        L = 2 * math.pi * r
        A = math.pi * r ** 2
        result = ct.isoperimetric_inequality_check(L, A)
        assert result["satisfied"]
        assert_close(result["ratio"], 1.0, tol=1e-6, msg="Kreis: Q=1")
        assert result["is_circle"]

    def test_isoperimetric_square(self, ct):
        """Quadrat mit Seite a=1: L=4, A=1. Q = 4π/16 ≈ 0.785 < 1."""
        a = 1.0
        L = 4 * a
        A = a ** 2
        result = ct.isoperimetric_inequality_check(L, A)
        assert result["satisfied"]
        assert result["ratio"] < 1.0, "Quadrat hat Q < 1"
        assert not result["is_circle"]

    def test_isoperimetric_satisfied_always(self, ct):
        """Isoperimetrische Ungleichung ist für beliebige L, A mit 4πA ≤ L² erfüllt."""
        # Verschiedene Formen
        for L, A in [(10, 5), (20, 20), (2 * math.pi, math.pi)]:
            result = ct.isoperimetric_inequality_check(L, A)
            if 4 * math.pi * A <= L ** 2:
                assert result["satisfied"], f"Soll erfüllt sein: L={L}, A={A}"

    def test_isoperimetric_inequality_values(self, ct):
        """Isoperimetrische Ungleichung: lhs und rhs korrekt berechnet."""
        L, A = 10.0, 5.0
        result = ct.isoperimetric_inequality_check(L, A)
        assert_close(result["lhs"], 4 * math.pi * A, tol=1e-10, msg="lhs=4πA")
        assert_close(result["rhs"], L ** 2, tol=1e-10, msg="rhs=L²")

    # --- four_vertex_theorem_demo ---

    def test_four_vertex_theorem(self, ct):
        """Vier-Scheitel-Satz: Ellipse hat ≥ 4 Scheitelpunkte."""
        result = ct.four_vertex_theorem_demo()
        assert result["satisfies_theorem"]
        assert result["num_vertices_found"] >= 4, \
            f"Nur {result['num_vertices_found']} Scheitel gefunden"

    def test_four_vertex_returns_dict(self, ct):
        """four_vertex_theorem_demo() gibt Dictionary mit korrekten Schlüsseln zurück."""
        result = ct.four_vertex_theorem_demo()
        assert "theorem" in result
        assert "kappa_max" in result
        assert "kappa_min" in result
        assert result["kappa_max"] > result["kappa_min"]

    # --- total_curvature_closed_curve ---

    def test_total_curvature_unit_circle(self, ct):
        """Gesamtkrümmung des Einheitskreises ist 2π."""
        def unit_circle(t):
            return [math.cos(t), math.sin(t)]
        total = ct.total_curvature_closed_curve(unit_circle, (0, 2 * math.pi))
        assert_close(abs(total), 2 * math.pi, tol=0.01, msg="Gesamtkrümmung Kreis=2π")

    def test_total_curvature_ellipse(self, ct):
        """Gesamtkrümmung der Ellipse ist 2π."""
        def ellipse(t):
            return [2.0 * math.cos(t), math.sin(t)]
        total = ct.total_curvature_closed_curve(ellipse, (0, 2 * math.pi))
        assert_close(abs(total), 2 * math.pi, tol=0.01, msg="Gesamtkrümmung Ellipse=2π")


# ============================================================================
# Tests: Edge Cases
# ============================================================================

class TestEdgeCases:
    """Edge Cases und Grenzfälle."""

    def test_curve_at_parameter_boundary(self):
        """Kurve an den Randpunkten des Parameterbereichs."""
        t = sp.Symbol('t')
        curve = ParametricCurve([t, t ** 2, sp.Integer(0)], t, (0, 1))
        # Kein Fehler an den Rändern
        v0 = curve.velocity(0.0)
        v1 = curve.velocity(1.0)
        assert v0 is not None
        assert v1 is not None

    def test_surface_degenerate_point(self):
        """Fläche an einem degenerierten Punkt (Pol der Sphäre)."""
        def sphere(theta, phi):
            return np.array([
                math.sin(theta) * math.cos(phi),
                math.sin(theta) * math.sin(phi),
                math.cos(theta),
            ])
        S = ParametricSurface(sphere, (0, math.pi), (0, 2 * math.pi))
        # Am Nordpol θ=0 (entarteter Punkt), kein Crash
        try:
            K = S.gaussian_curvature(0.01, 0.0)  # Nähe des Nordpols
            # K kann ungenau sein, aber kein Fehler
            assert isinstance(K, float)
        except Exception as e:
            pytest.fail(f"Fehler an degenerierten Punkt: {e}")

    def test_zero_perimeter_isoperimetric(self):
        """Isoperimetrische Ungleichung mit L=0."""
        ct = CurveTheory2D()
        result = ct.isoperimetric_inequality_check(0.0, 0.0)
        assert isinstance(result["ratio"], float)

    def test_negative_area_isoperimetric(self):
        """Isoperimetrische Ungleichung mit negativer Fläche."""
        ct = CurveTheory2D()
        result = ct.isoperimetric_inequality_check(10.0, -1.0)
        # Negativer Flächeninhalt → 4πA < 0 ≤ L² → erfüllt
        assert result["satisfied"]

    def test_principal_curvatures_relation(self):
        """K = κ₁ · κ₂ und H = (κ₁+κ₂)/2 für die Sphäre.

        Die Relationen gelten unabhängig vom Vorzeichen (Normalenorientierung).
        """
        def sphere(theta, phi):
            return np.array([
                math.sin(theta) * math.cos(phi),
                math.sin(theta) * math.sin(phi),
                math.cos(theta),
            ])
        S = ParametricSurface(sphere, (0, math.pi), (0, 2 * math.pi))
        u, v = math.pi / 2, 0.0
        k1, k2 = S.principal_curvatures(u, v)
        K = S.gaussian_curvature(u, v)
        H = S.mean_curvature(u, v)
        # K = κ₁ · κ₂ (Vorzeichen konsistent da beide aus gleicher Fundamentalform)
        assert_close(k1 * k2, K, tol=0.1, msg="K = κ₁·κ₂")
        # H = (κ₁ + κ₂) / 2 (Vorzeichen konsistent)
        assert_close((k1 + k2) / 2, H, tol=0.1, msg="H = (κ₁+κ₂)/2")

    def test_frenet_formula_keys(self):
        """frenet_serret_formulas() enthält alle erwarteten Schlüssel."""
        t = sp.Symbol('t')
        curve = ParametricCurve([sp.cos(t), sp.sin(t), t], t)
        formulas = curve.frenet_serret_formulas()
        for key in ["T_prime", "N_prime", "B_prime", "description"]:
            assert key in formulas, f"Schlüssel '{key}' fehlt in frenet_serret_formulas"

    def test_geodesic_sphere_angle_zero_direction(self):
        """Geodätische mit Nullrichtung auf S² – keine Exception."""
        gc = GeodesicComputation()
        try:
            result = gc.geodesic_sphere(math.pi / 2, 0.0, [0.0, 0.0], t_max=0.5)
            assert "theta" in result
        except Exception as e:
            pytest.fail(f"Exception bei Nullrichtung: {e}")

    def test_ricci_tensor_shape(self):
        """Ricci-Tensor hat korrekte Form (n×n)."""
        rg = RiemannianGeometry()
        R_flat = np.zeros((3, 3, 3, 3))
        metric = np.eye(3)
        Ric = rg.ricci_curvature(R_flat, metric)
        assert Ric.shape == (3, 3), f"Ricci-Form {Ric.shape} ≠ (3,3)"

    def test_arc_length_zero_interval(self):
        """Bogenlänge über Intervall der Länge 0 ist 0."""
        t = sp.Symbol('t')
        curve = ParametricCurve([sp.cos(t), sp.sin(t), t], t)
        L = curve.arc_length(1.0, 1.0)
        assert_close(L, 0.0, tol=1e-10, msg="Bogenlänge leeres Intervall = 0")

    def test_curvature_radius_one_circle(self):
        """Kreis mit Radius R=2: κ = 1/R = 0.5."""
        t = sp.Symbol('t')
        R = 2.0
        curve = ParametricCurve([R * sp.cos(t), R * sp.sin(t), sp.Integer(0)], t)
        kappa = curve.curvature(0.0)
        assert_close(kappa, 1.0 / R, tol=1e-4, msg="κ=1/R Kreis mit R=2")


# ============================================================================
# Tests: Numerische Konsistenz
# ============================================================================

class TestNumericalConsistency:
    """Tests zur Überprüfung numerischer Konsistenz."""

    def test_bogenlange_additiv(self):
        """Bogenlänge ist additiv: L(a,c) = L(a,b) + L(b,c)."""
        t = sp.Symbol('t')
        curve = ParametricCurve([sp.cos(t), sp.sin(t), t], t)
        a, b, c = 0.0, 1.0, 2.0
        L_ac = curve.arc_length(a, c)
        L_ab = curve.arc_length(a, b)
        L_bc = curve.arc_length(b, c)
        assert_close(L_ac, L_ab + L_bc, tol=1e-5, msg="Additivität der Bogenlänge")

    def test_gaussian_curvature_sphere_consistent(self):
        """K der Sphäre ist an mehreren Punkten konsistent K=1."""
        def sphere(theta, phi):
            return np.array([math.sin(theta)*math.cos(phi),
                             math.sin(theta)*math.sin(phi), math.cos(theta)])
        S = ParametricSurface(sphere, (0, math.pi), (0, 2 * math.pi))
        K_values = [
            S.gaussian_curvature(theta, phi)
            for theta in [math.pi / 4, math.pi / 2, 3 * math.pi / 4]
            for phi in [0.0, math.pi / 2, math.pi]
        ]
        for K in K_values:
            assert_close(K, 1.0, tol=0.05, msg=f"K={K} Einheitssphäre")

    def test_torsion_planar_curve_zero(self):
        """Ebene Kurve hat τ = 0."""
        t = sp.Symbol('t')
        # Beliebige ebene Kurve: z = 0
        curve = ParametricCurve([t ** 2, t ** 3, sp.Integer(0)], t)
        tau = curve.torsion(1.0)
        assert_close(tau, 0.0, tol=1e-4, msg="τ=0 für ebene Kurve")

    def test_speed_velocity_relation(self):
        """speed(t) = |velocity(t)|."""
        t = sp.Symbol('t')
        curve = ParametricCurve([sp.cos(t), sp.sin(t), t], t)
        for t_val in [0.0, 0.5, 1.0]:
            s = curve.speed(t_val)
            v = curve.velocity(t_val)
            assert_close(s, float(np.linalg.norm(v)), tol=1e-10,
                         msg=f"speed=|velocity| bei t={t_val}")


# ============================================================================
# Hilfsfunktion für Tests
# ============================================================================

def _norm(v) -> float:
    """Norm-Hilfsfunktion für Tests."""
    return float(np.linalg.norm(v))
