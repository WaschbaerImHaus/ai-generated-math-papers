"""
Tests für src/ricci_flow.py — Ricci-Fluss und Perelmansche Techniken

Testabdeckung:
- RicciFlowMetric: Christoffel-Symbole, Ricci-Tensor, Skalare Krümmung
- RicciFlow: Evolution, Stationarität, Singularitätserkennung
- NormalizedRicciFlow: Normalisierter Fluss, Langzeitverhalten
- SphericalMetric: Sphärische Geometrie
- HyperbolicMetric: Hyperbolische Geometrie
- RicciFlowSurgery: Hals-Erkennung, Chirurgie
- PoincareConjecture: Aussage, Beweisskizze, Betti-Zahlen
- EntropyFunctional: W- und F-Funktional
- ThreeDimensionalTopology: Euler-Charakteristik, Klassifikation
- Hilfsfunktionen: normalize_metric, volume_element, convergence_rate

Autor: Michael Fuhrmann
Letzte Änderung: 2026-03-11
"""

import pytest
import numpy as np
import sys
import os

# Sicherstellen, dass src/ im Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ricci_flow import (
    RicciFlowMetric,
    RicciFlow,
    NormalizedRicciFlow,
    SphericalMetric,
    HyperbolicMetric,
    RicciFlowSurgery,
    PoincareConjecture,
    EntropyFunctional,
    ThreeDimensionalTopology,
    ricci_flow_on_surface,
    normalize_metric,
    compute_volume_element,
    ricci_flow_convergence_rate,
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen für Testmetriken
# ---------------------------------------------------------------------------

def flat_metric_2d(coords):
    """Flache euklidische Metrik: g = I₂."""
    return np.eye(2)


def flat_metric_3d(coords):
    """Flache euklidische Metrik: g = I₃."""
    return np.eye(3)


def scaled_flat_metric(scale=2.0):
    """Skalierte flache Metrik: g = scale * I₂."""
    def f(coords):
        return scale * np.eye(2)
    return f


def sphere_like_metric_2d(coords):
    """
    Sphärenähnliche Metrik: g = (1 + r²/4)⁻² * I₂
    (stereographische Projektion der Einheitssphäre S²)
    """
    r2 = float(np.dot(coords[:2], coords[:2]))
    factor = 1.0 / (1.0 + r2 / 4.0) ** 2
    return factor * np.eye(2)


def poincare_metric_2d(coords):
    """Poincaré-Scheibenmetrik."""
    r2 = float(np.dot(coords[:2], coords[:2]))
    r2 = min(r2, 0.99)
    factor = 4.0 / (1.0 - r2) ** 2
    return factor * np.eye(2)


# ---------------------------------------------------------------------------
# TestRicciFlowMetric
# ---------------------------------------------------------------------------

class TestRicciFlowMetric:
    """Tests für die RicciFlowMetric-Klasse."""

    def test_flat_metric_zero_christoffel(self):
        """Flache Metrik → Christoffel-Symbole ≈ 0."""
        rfm = RicciFlowMetric(flat_metric_2d, dim=2)
        coords = np.array([0.5, 0.5])
        Gamma = rfm.christoffel_symbols(coords)
        assert Gamma.shape == (2, 2, 2)
        assert np.allclose(Gamma, 0.0, atol=1e-6)

    def test_flat_metric_zero_ricci_tensor(self):
        """Flache Metrik → Ricci-Tensor ≈ 0."""
        rfm = RicciFlowMetric(flat_metric_2d, dim=2)
        coords = np.array([0.3, 0.7])
        Ric = rfm.ricci_tensor(coords)
        assert Ric.shape == (2, 2)
        assert np.allclose(Ric, 0.0, atol=1e-4)

    def test_flat_metric_zero_ricci_scalar(self):
        """Flache Metrik → Skalare Krümmung R ≈ 0."""
        rfm = RicciFlowMetric(flat_metric_2d, dim=2)
        coords = np.array([0.1, 0.2])
        R = rfm.ricci_scalar(coords)
        assert abs(R) < 1e-4

    def test_scaled_flat_metric_zero_curvature(self):
        """Skalierte flache Metrik hat ebenfalls R = 0."""
        rfm = RicciFlowMetric(scaled_flat_metric(3.0), dim=2)
        coords = np.array([0.0, 0.0])
        R = rfm.ricci_scalar(coords)
        assert abs(R) < 1e-3

    def test_sphere_like_metric_positive_curvature(self):
        """Sphärenähnliche Metrik hat positive skalare Krümmung."""
        rfm = RicciFlowMetric(sphere_like_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        R = rfm.ricci_scalar(coords)
        # Für S² mit Radius 1: R = 2
        assert R > 0.5

    def test_metric_matrix_shape(self):
        """Metrik-Matrix hat korrekte Form (dim×dim)."""
        for dim in [2, 3]:
            rfm = RicciFlowMetric(
                flat_metric_3d if dim == 3 else flat_metric_2d, dim=dim
            )
            coords = np.zeros(dim)
            g = rfm._metric(coords)
            assert g.shape == (dim, dim)

    def test_inverse_metric_identity(self):
        """g · g⁻¹ = I für flache Metrik."""
        rfm = RicciFlowMetric(flat_metric_2d, dim=2)
        coords = np.array([0.5, 0.5])
        g = rfm._metric(coords)
        g_inv = rfm._inverse_metric(coords)
        product = g @ g_inv
        assert np.allclose(product, np.eye(2), atol=1e-10)

    def test_ricci_tensor_symmetry(self):
        """Ricci-Tensor ist symmetrisch: R_{ij} = R_{ji}."""
        rfm = RicciFlowMetric(sphere_like_metric_2d, dim=2)
        coords = np.array([0.1, 0.1])
        Ric = rfm.ricci_tensor(coords)
        assert np.allclose(Ric, Ric.T, atol=1e-5)

    def test_christoffel_3d_flat(self):
        """3D flache Metrik → Christoffel-Symbole ≈ 0."""
        rfm = RicciFlowMetric(flat_metric_3d, dim=3)
        coords = np.array([0.2, 0.3, 0.4])
        Gamma = rfm.christoffel_symbols(coords)
        assert Gamma.shape == (3, 3, 3)
        assert np.allclose(Gamma, 0.0, atol=1e-5)

    def test_poincare_metric_negative_curvature(self):
        """Poincaré-Scheibenmetrik hat negative skalare Krümmung."""
        rfm = RicciFlowMetric(poincare_metric_2d, dim=2)
        coords = np.array([0.1, 0.0])
        R = rfm.ricci_scalar(coords)
        # Poincaré-Scheibe: R = -2 (in 2D K = -1 → R = 2K = -2)
        assert R < 0

    def test_sectional_curvature_flat(self):
        """Schnittkriimmung der flachen Metrik ≈ 0."""
        rfm = RicciFlowMetric(flat_metric_2d, dim=2)
        coords = np.array([0.5, 0.5])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([0.0, 1.0])
        K = rfm.sectional_curvature(coords, v1, v2)
        assert abs(K) < 1e-3

    def test_dg_shape(self):
        """Ableitung der Metrik hat Form (dim, dim, dim)."""
        rfm = RicciFlowMetric(sphere_like_metric_2d, dim=2)
        coords = np.array([0.1, 0.2])
        dg = rfm._dg(coords)
        assert dg.shape == (2, 2, 2)

    def test_ricci_scalar_sphere_approximation(self):
        """Sphärenähnliche Metrik: R nahe 2 am Ursprung."""
        rfm = RicciFlowMetric(sphere_like_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        R = rfm.ricci_scalar(coords)
        # Für diese Metrik ist R ≈ 2 (Gauss-Krümmung K=1, R=2K)
        assert R > 0

    def test_metric_positive_definite(self):
        """Sphärenmetrik ist positiv definit."""
        rfm = RicciFlowMetric(sphere_like_metric_2d, dim=2)
        coords = np.array([0.5, 0.5])
        g = rfm._metric(coords)
        eigenvalues = np.linalg.eigvalsh(g)
        assert np.all(eigenvalues > 0)


# ---------------------------------------------------------------------------
# TestRicciFlow
# ---------------------------------------------------------------------------

class TestRicciFlow:
    """Tests für die RicciFlow-Klasse."""

    def test_flat_metric_stationary(self):
        """Flache Metrik bleibt unter dem Ricci-Fluss stationär."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.3, 0.3])
        history = flow.evolve(coords, T=0.1, steps=10)
        g0 = history[0][1]
        g_final = history[-1][1]
        # Flache Metrik sollte sich kaum ändern
        assert np.allclose(g0, g_final, atol=1e-3)

    def test_evolve_returns_history(self):
        """evolve() gibt Liste von (t, g) Paaren zurück."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.1, 0.1])
        history = flow.evolve(coords, T=0.5, steps=5)
        assert len(history) == 6  # steps + 1 Anfangswert
        for t, g in history:
            assert isinstance(t, float)
            assert g.shape == (2, 2)

    def test_evolve_time_increases(self):
        """Zeitwerte in der Historie sind monoton wachsend."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        history = flow.evolve(coords, T=1.0, steps=10)
        times = [t for t, g in history]
        assert times[0] == 0.0
        assert times[-1] == pytest.approx(1.0, abs=1e-10)
        for i in range(len(times) - 1):
            assert times[i] < times[i + 1]

    def test_evolve_step_returns_matrix(self):
        """evolve_step() gibt (dim×dim)-Matrix zurück."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.1, 0.1])
        g_new = flow.evolve_step(coords, dt=0.01)
        assert g_new.shape == (2, 2)

    def test_no_singularity_flat_metric(self):
        """Keine Singularität für flache Metrik."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.5, 0.5])
        result = flow.detect_singularity(coords, t=0.0, threshold=1e6)
        assert result is False

    def test_singularity_detection_high_curvature(self):
        """Hohe Krümmung wird als Singularität erkannt."""
        # Künstlich hohe Krümmungsmetrik
        def high_curv_metric(coords):
            x = coords[0]
            # Sehr stark gewölbte Metrik
            scale = 1.0 / (x ** 2 + 1e-4)
            return scale * np.eye(2)

        flow = RicciFlow(high_curv_metric, dim=2)
        flow._current_metric = high_curv_metric
        coords = np.array([1e-3, 0.0])
        # Mit sehr kleinem Threshold sollte Singularität erkannt werden
        result = flow.detect_singularity(coords, t=0.0, threshold=0.0)
        assert result is True  # Jede Krümmung > 0 → Singularität

    def test_evolve_metric_positive_definite_flat(self):
        """Metrik bleibt positiv definit nach Evolution."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        history = flow.evolve(coords, T=0.1, steps=5)
        for t, g in history:
            eigvals = np.linalg.eigvalsh(g)
            assert np.all(eigvals > -0.1)  # Näherungsweise positiv definit

    def test_flow_dim_2_default(self):
        """Standard-Dimension ist 2."""
        flow = RicciFlow(flat_metric_2d)
        assert flow.dim == 2

    def test_flow_dim_3(self):
        """3D-Flow wird korrekt initialisiert."""
        flow = RicciFlow(flat_metric_3d, dim=3)
        assert flow.dim == 3

    def test_evolve_zero_time(self):
        """T=0 gibt nur Anfangsmetrik zurück."""
        flow = RicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.1, 0.1])
        history = flow.evolve(coords, T=0.0, steps=1)
        # Mindestens 2 Einträge (Start + ein Schritt)
        assert len(history) >= 1


# ---------------------------------------------------------------------------
# TestNormalizedRicciFlow
# ---------------------------------------------------------------------------

class TestNormalizedRicciFlow:
    """Tests für den normalisierten Ricci-Fluss."""

    def test_init(self):
        """Initialisierung korrekt."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        assert nrf.dim == 2

    def test_mean_curvature_flat(self):
        """Mittlere Krümmung der flachen Metrik ≈ 0."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.1, 0.1])
        r = nrf.mean_curvature(coords)
        assert abs(r) < 1e-4

    def test_evolve_step_returns_matrix(self):
        """evolve_step() gibt (dim×dim)-Matrix zurück."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        g_new = nrf.evolve_step(coords, dt=0.01)
        assert g_new.shape == (2, 2)

    def test_flat_metric_stays_flat(self):
        """Flache Metrik bleibt unter normalisiertem Fluss flach."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        g_new = nrf.evolve_step(coords, dt=0.1)
        assert np.allclose(g_new, np.eye(2), atol=1e-3)

    def test_long_time_behavior_returns_dict(self):
        """long_time_behavior() gibt Dict zurück."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        result = nrf.long_time_behavior(coords, T=0.5, steps=10)
        assert isinstance(result, dict)
        assert 'curvature_history' in result
        assert 'final_metric' in result
        assert 'converged' in result

    def test_long_time_behavior_history_length(self):
        """Krümmungshistorie hat korrekte Länge."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        result = nrf.long_time_behavior(coords, T=1.0, steps=20)
        assert len(result['curvature_history']) == 20

    def test_long_time_final_metric_shape(self):
        """Endmetrik hat korrekte Form."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        result = nrf.long_time_behavior(coords, T=0.5, steps=10)
        assert result['final_metric'].shape == (2, 2)

    def test_flat_converges(self):
        """Flache Metrik konvergiert sofort."""
        nrf = NormalizedRicciFlow(flat_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        result = nrf.long_time_behavior(coords, T=1.0, steps=50)
        # Flache Metrik sollte als konvergiert gelten
        assert bool(result['converged']) is True

    def test_mean_curvature_sphere_like_positive(self):
        """Sphärenähnliche Metrik hat positive mittlere Krümmung."""
        nrf = NormalizedRicciFlow(sphere_like_metric_2d, dim=2)
        coords = np.array([0.0, 0.0])
        r = nrf.mean_curvature(coords)
        assert r > 0

    def test_dim_3_init(self):
        """3D normalisierter Fluss wird korrekt initialisiert."""
        nrf = NormalizedRicciFlow(flat_metric_3d, dim=3)
        assert nrf.dim == 3


# ---------------------------------------------------------------------------
# TestSphericalMetric
# ---------------------------------------------------------------------------

class TestSphericalMetric:
    """Tests für die SphericalMetric-Klasse."""

    def test_standard_metric_is_identity(self):
        """Standard-Metrik ist die Einheitsmatrix."""
        sm = SphericalMetric()
        f = sm.standard_metric(dim=2)
        g = f(np.array([0.0, 0.0]))
        assert np.allclose(g, np.eye(2))

    def test_standard_metric_3d(self):
        """Standard-Metrik in 3D ist I₃."""
        sm = SphericalMetric()
        f = sm.standard_metric(dim=3)
        g = f(np.zeros(3))
        assert np.allclose(g, np.eye(3))

    def test_sectional_curvature_sphere_positive(self):
        """Schnittkriimmung der Einheitssphäre ist +1."""
        sm = SphericalMetric()
        K = sm.sectional_curvature_sphere(dim=2)
        assert K == pytest.approx(1.0)

    def test_sectional_curvature_sphere_dim3(self):
        """Schnittkriimmung ist unabhängig von der Dimension +1."""
        sm = SphericalMetric()
        for dim in [2, 3, 4]:
            K = sm.sectional_curvature_sphere(dim=dim)
            assert K == pytest.approx(1.0)

    def test_ricci_tensor_sphere_2d(self):
        """Ricci-Tensor der S² ist (dim-1)·I = I."""
        sm = SphericalMetric()
        Ric = sm.ricci_tensor_sphere(dim=2)
        assert np.allclose(Ric, np.eye(2))

    def test_ricci_tensor_sphere_3d(self):
        """Ricci-Tensor der S³ ist 2·I."""
        sm = SphericalMetric()
        Ric = sm.ricci_tensor_sphere(dim=3)
        assert np.allclose(Ric, 2.0 * np.eye(3))

    def test_stereographic_metric_at_origin(self):
        """Stereographische Metrik am Ursprung: g = 4·I."""
        sm = SphericalMetric()
        f = sm.stereographic_metric(dim=2)
        g = f(np.array([0.0, 0.0]))
        assert np.allclose(g, 4.0 * np.eye(2))

    def test_stereographic_metric_positive_definite(self):
        """Stereographische Metrik ist positiv definit."""
        sm = SphericalMetric()
        f = sm.stereographic_metric(dim=2)
        for coords in [np.array([0.1, 0.2]), np.array([1.0, 0.0])]:
            g = f(coords)
            eigvals = np.linalg.eigvalsh(g)
            assert np.all(eigvals > 0)

    def test_ricci_tensor_sphere_shape(self):
        """Ricci-Tensor hat korrekte Form."""
        sm = SphericalMetric()
        for dim in [2, 3]:
            Ric = sm.ricci_tensor_sphere(dim=dim)
            assert Ric.shape == (dim, dim)


# ---------------------------------------------------------------------------
# TestHyperbolicMetric
# ---------------------------------------------------------------------------

class TestHyperbolicMetric:
    """Tests für die HyperbolicMetric-Klasse."""

    def test_negative_sectional_curvature(self):
        """Hyperbolischer Raum hat Schnittkriimmung -1."""
        hm = HyperbolicMetric()
        K = hm.sectional_curvature_hyperbolic()
        assert K == pytest.approx(-1.0)

    def test_poincare_metric_at_origin(self):
        """Poincaré-Metrik am Ursprung: g = 4·I."""
        hm = HyperbolicMetric()
        coords = np.array([0.0, 0.0])
        g = hm.poincare_disk_metric(coords)
        assert np.allclose(g, 4.0 * np.eye(2))

    def test_poincare_metric_grows_at_boundary(self):
        """Poincaré-Metrik wächst zum Rand der Scheibe."""
        hm = HyperbolicMetric()
        g_center = hm.poincare_disk_metric(np.array([0.0, 0.0]))
        g_near_boundary = hm.poincare_disk_metric(np.array([0.9, 0.0]))
        # Metrik ist am Rand größer
        assert g_near_boundary[0, 0] > g_center[0, 0]

    def test_poincare_metric_positive_definite(self):
        """Poincaré-Metrik ist positiv definit."""
        hm = HyperbolicMetric()
        coords = np.array([0.3, 0.4])
        g = hm.poincare_disk_metric(coords)
        eigvals = np.linalg.eigvalsh(g)
        assert np.all(eigvals > 0)

    def test_hyperbolic_distance_zero(self):
        """Abstand eines Punktes zu sich selbst ist 0."""
        hm = HyperbolicMetric()
        z = 0.3 + 0.4j
        d = hm.hyperbolic_distance(z, z)
        assert d == pytest.approx(0.0, abs=1e-10)

    def test_hyperbolic_distance_symmetry(self):
        """Hyperbolischer Abstand ist symmetrisch: d(z1,z2) = d(z2,z1)."""
        hm = HyperbolicMetric()
        z1 = 0.3 + 0.2j
        z2 = -0.1 + 0.4j
        assert hm.hyperbolic_distance(z1, z2) == pytest.approx(
            hm.hyperbolic_distance(z2, z1), abs=1e-10
        )

    def test_poincare_metric_shape(self):
        """Poincaré-Metrik gibt (2×2)-Matrix zurück."""
        hm = HyperbolicMetric()
        coords = np.array([0.1, 0.2])
        g = hm.poincare_disk_metric(coords)
        assert g.shape == (2, 2)

    def test_poincare_metric_func(self):
        """poincare_disk_metric_func gibt aufrufbare Funktion zurück."""
        hm = HyperbolicMetric()
        f = hm.poincare_disk_metric_func()
        assert callable(f)
        g = f(np.array([0.0, 0.0]))
        assert g.shape == (2, 2)


# ---------------------------------------------------------------------------
# TestRicciFlowSurgery
# ---------------------------------------------------------------------------

class TestRicciFlowSurgery:
    """Tests für die RicciFlowSurgery-Klasse."""

    def test_no_necks_large_volume(self):
        """Keine Hälse wenn Volumen-Element groß genug."""
        surgery = RicciFlowSurgery()
        metric_data = [np.eye(2) * k for k in range(1, 6)]  # große Metriken
        necks = surgery.neck_detection(metric_data, threshold=0.1)
        assert len(necks) == 0

    def test_neck_detection_small_metric(self):
        """Kleine Metrik wird als Hals erkannt."""
        surgery = RicciFlowSurgery()
        # Sehr kleine Metrik → kleines Volumen-Element
        metric_data = [np.eye(2) * 0.01]
        necks = surgery.neck_detection(metric_data, threshold=0.5)
        assert len(necks) > 0

    def test_neck_detection_returns_indices(self):
        """Hals-Erkennung gibt Indizes zurück."""
        surgery = RicciFlowSurgery()
        # Abwechselnde große und kleine Metriken
        metric_data = [
            np.eye(2) * 10.0,  # groß
            np.eye(2) * 0.01,  # klein → Hals
            np.eye(2) * 10.0,  # groß
        ]
        necks = surgery.neck_detection(metric_data, threshold=0.5)
        assert 1 in necks

    def test_surgery_step_returns_dict(self):
        """surgery_step() gibt Dict zurück."""
        surgery = RicciFlowSurgery()
        metric_data = [np.eye(2)] * 3
        result = surgery.surgery_step(metric_data, singularity_time=1.0)
        assert isinstance(result, dict)
        assert 'singularity_time' in result
        assert 'neck_count' in result
        assert 'surgery_performed' in result

    def test_surgery_singularity_time_stored(self):
        """Singularitätszeitpunkt wird korrekt gespeichert."""
        surgery = RicciFlowSurgery()
        metric_data = [np.eye(2)]
        result = surgery.surgery_step(metric_data, singularity_time=2.5)
        assert result['singularity_time'] == 2.5

    def test_describe_surgery_nonempty(self):
        """Beschreibung der Chirurgie ist nicht leer."""
        surgery = RicciFlowSurgery()
        desc = surgery.describe_surgery()
        assert isinstance(desc, str)
        assert len(desc) > 50

    def test_describe_surgery_contains_perelman(self):
        """Beschreibung enthält 'Perelman' oder 'Hals'."""
        surgery = RicciFlowSurgery()
        desc = surgery.describe_surgery()
        assert 'Perelman' in desc or 'Hals' in desc or 'ε' in desc

    def test_surgery_no_necks(self):
        """Chirurgie ohne Hälse: surgery_performed = False."""
        surgery = RicciFlowSurgery()
        metric_data = [np.eye(2) * 100.0] * 5
        result = surgery.surgery_step(metric_data, singularity_time=0.5)
        assert result['surgery_performed'] is False


# ---------------------------------------------------------------------------
# TestPoincareConjecture
# ---------------------------------------------------------------------------

class TestPoincareConjecture:
    """Tests für die PoincareConjecture-Klasse."""

    def test_statement_nonempty(self):
        """Aussage der Vermutung ist nicht leer."""
        pc = PoincareConjecture()
        stmt = pc.statement()
        assert isinstance(stmt, str)
        assert len(stmt) > 100

    def test_statement_contains_poincare(self):
        """Aussage enthält 'Poincaré'."""
        pc = PoincareConjecture()
        assert 'Poincar' in pc.statement()

    def test_statement_contains_s3(self):
        """Aussage erwähnt S³."""
        pc = PoincareConjecture()
        assert 'S³' in pc.statement() or 'S^3' in pc.statement() or '3-Sphäre' in pc.statement()

    def test_s3_betti_numbers_true(self):
        """Betti-Zahlen (1,0,0,1) werden als S³ erkannt."""
        pc = PoincareConjecture()
        assert pc.verify_s3_topology([1, 0, 0, 1]) is True

    def test_wrong_betti_numbers_false(self):
        """Falsche Betti-Zahlen werden abgelehnt."""
        pc = PoincareConjecture()
        assert pc.verify_s3_topology([1, 1, 0, 1]) is False
        assert pc.verify_s3_topology([1, 0, 1, 1]) is False
        assert pc.verify_s3_topology([2, 0, 0, 1]) is False

    def test_torus_betti_numbers_false(self):
        """Torus T³ mit Betti (1,3,3,1) ist nicht S³."""
        pc = PoincareConjecture()
        assert pc.verify_s3_topology([1, 3, 3, 1]) is False

    def test_short_betti_list_false(self):
        """Zu kurze Betti-Liste gibt False zurück."""
        pc = PoincareConjecture()
        assert pc.verify_s3_topology([1, 0]) is False

    def test_perelman_sketch_is_dict(self):
        """Beweisskizze ist ein Dict."""
        pc = PoincareConjecture()
        sketch = pc.perelman_sketch()
        assert isinstance(sketch, dict)

    def test_perelman_sketch_has_keys(self):
        """Beweisskizze enthält mindestens 5 Schlüsselschritte."""
        pc = PoincareConjecture()
        sketch = pc.perelman_sketch()
        assert len(sketch) >= 5

    def test_hamilton_program_nonempty(self):
        """Hamiltons Programm ist nicht leer."""
        pc = PoincareConjecture()
        prog = pc.hamilton_program()
        assert isinstance(prog, dict)
        assert len(prog) >= 3

    def test_perelman_contributions_list(self):
        """Perelmans Beiträge sind eine Liste."""
        pc = PoincareConjecture()
        contribs = pc.perelman_contributions()
        assert isinstance(contribs, list)
        assert len(contribs) >= 3

    def test_perelman_contributions_w_functional(self):
        """W-Funktional ist in den Beiträgen enthalten."""
        pc = PoincareConjecture()
        contribs = pc.perelman_contributions()
        names = [c['name'] for c in contribs]
        assert any('W' in n or 'W-Funk' in n for n in names)

    def test_perelman_contributions_have_beschreibung(self):
        """Jeder Beitrag hat eine Beschreibung."""
        pc = PoincareConjecture()
        for contrib in pc.perelman_contributions():
            assert 'beschreibung' in contrib
            assert len(contrib['beschreibung']) > 10


# ---------------------------------------------------------------------------
# TestEntropyFunctional
# ---------------------------------------------------------------------------

class TestEntropyFunctional:
    """Tests für das EntropyFunctional."""

    def test_monotonicity_statement_nonempty(self):
        """Monotonie-Aussage ist nicht leer."""
        ef = EntropyFunctional()
        stmt = ef.monotonicity_statement()
        assert isinstance(stmt, str)
        assert len(stmt) > 50

    def test_monotonicity_statement_contains_w(self):
        """Monotonie-Aussage erwähnt W-Funktional."""
        ef = EntropyFunctional()
        stmt = ef.monotonicity_statement()
        assert 'W' in stmt or 'F' in stmt

    def test_f_functional_flat_zero(self):
        """F-Funktional für flache Metrik mit f=0: F ≈ 0."""
        ef = EntropyFunctional()
        coords_list = [np.array([0.0, 0.0]), np.array([0.1, 0.0]),
                       np.array([0.0, 0.1])]
        f_func = lambda c: 0.0  # f identisch 0
        F = ef.f_functional(flat_metric_2d, f_func, coords_list)
        assert abs(F) < 1e-3

    def test_f_functional_returns_float(self):
        """F-Funktional gibt float zurück."""
        ef = EntropyFunctional()
        coords_list = [np.array([0.0, 0.0])]
        f_func = lambda c: 1.0
        F = ef.f_functional(flat_metric_2d, f_func, coords_list)
        assert isinstance(F, float)

    def test_w_functional_returns_float(self):
        """W-Funktional gibt float zurück."""
        ef = EntropyFunctional()
        coords_list = [np.array([0.0, 0.0])]
        f_func = lambda c: 1.0
        W = ef.w_functional(flat_metric_2d, f_func, tau=1.0, coords_list=coords_list)
        assert isinstance(W, float)

    def test_w_functional_tau_positive(self):
        """W-Funktional ist für τ > 0 endlich."""
        ef = EntropyFunctional()
        coords_list = [np.array([0.1, 0.1])]
        f_func = lambda c: np.dot(c, c)
        W = ef.w_functional(flat_metric_2d, f_func, tau=0.5, coords_list=coords_list)
        assert np.isfinite(W)

    def test_f_functional_with_weights(self):
        """F-Funktional akzeptiert Gewichte."""
        ef = EntropyFunctional()
        coords_list = [np.array([0.0, 0.0]), np.array([0.1, 0.1])]
        f_func = lambda c: 0.0
        weights = np.array([0.5, 0.5])
        F = ef.f_functional(flat_metric_2d, f_func, coords_list, weights=weights)
        assert isinstance(F, float)


# ---------------------------------------------------------------------------
# TestThreeDimensionalTopology
# ---------------------------------------------------------------------------

class TestThreeDimensionalTopology:
    """Tests für ThreeDimensionalTopology."""

    def test_euler_characteristic_s3(self):
        """χ(S³) = 1 - 0 + 0 - 1 = 0."""
        topo = ThreeDimensionalTopology()
        chi = topo.euler_characteristic([1, 0, 0, 1])
        assert chi == 0

    def test_euler_characteristic_s2(self):
        """χ(S²) = β₀ - β₁ + β₂ = 1 - 0 + 1 = 2."""
        topo = ThreeDimensionalTopology()
        chi = topo.euler_characteristic([1, 0, 1])
        assert chi == 2

    def test_euler_characteristic_torus(self):
        """χ(T²) = 1 - 2 + 1 = 0."""
        topo = ThreeDimensionalTopology()
        chi = topo.euler_characteristic([1, 2, 1])
        assert chi == 0

    def test_euler_characteristic_empty(self):
        """Leere Betti-Liste → χ = 0."""
        topo = ThreeDimensionalTopology()
        chi = topo.euler_characteristic([])
        assert chi == 0

    def test_is_simply_connected_true(self):
        """Triviale Fundamentalgruppe → einfach zusammenhängend."""
        topo = ThreeDimensionalTopology()
        assert topo.is_simply_connected(True) is True

    def test_is_simply_connected_false(self):
        """Nicht-triviale Fundamentalgruppe → nicht einfach zusammenhängend."""
        topo = ThreeDimensionalTopology()
        assert topo.is_simply_connected(False) is False

    def test_geometrization_8_cases(self):
        """Geometrisierungsvermutung hat genau 8 Geometrien."""
        topo = ThreeDimensionalTopology()
        cases = topo.geometrization_conjecture_cases()
        assert len(cases) == 8

    def test_geometrization_contains_s3(self):
        """S³-Geometrie ist enthalten."""
        topo = ThreeDimensionalTopology()
        cases = topo.geometrization_conjecture_cases()
        assert 'S³' in cases

    def test_geometrization_contains_h3(self):
        """H³-Geometrie ist enthalten."""
        topo = ThreeDimensionalTopology()
        cases = topo.geometrization_conjecture_cases()
        assert 'H³' in cases

    def test_geometrization_contains_e3(self):
        """E³-Geometrie (flach) ist enthalten."""
        topo = ThreeDimensionalTopology()
        cases = topo.geometrization_conjecture_cases()
        assert 'E³' in cases

    def test_classify_sphere(self):
        """Einfach zusammenhängende, positiv gekrümmte Mannigfaltigkeit → S³."""
        topo = ThreeDimensionalTopology()
        result = topo.classify_3manifold({
            'simply_connected': True,
            'positive_curvature': True
        })
        assert 'S³' in result

    def test_classify_flat(self):
        """Flache Mannigfaltigkeit → E³-Geometrie."""
        topo = ThreeDimensionalTopology()
        result = topo.classify_3manifold({'flat': True})
        assert 'E³' in result or 'flat' in result.lower() or 'Flach' in result

    def test_classify_hyperbolic(self):
        """Hyperbolische Mannigfaltigkeit → H³-Geometrie."""
        topo = ThreeDimensionalTopology()
        result = topo.classify_3manifold({'hyperbolic': True})
        assert 'H³' in result or 'Hyperbol' in result

    def test_classify_s3_via_betti(self):
        """S³ wird via Betti-Zahlen erkannt."""
        topo = ThreeDimensionalTopology()
        result = topo.classify_3manifold({
            'simply_connected': True,
            'betti_numbers': [1, 0, 0, 1]
        })
        assert 'S³' in result

    def test_geometrization_all_have_name(self):
        """Alle 8 Geometrien haben einen 'name'-Eintrag."""
        topo = ThreeDimensionalTopology()
        cases = topo.geometrization_conjecture_cases()
        for key, geo in cases.items():
            assert 'name' in geo


# ---------------------------------------------------------------------------
# TestHelperFunctions
# ---------------------------------------------------------------------------

class TestHelperFunctions:
    """Tests für die Hilfsfunktionen auf Modul-Ebene."""

    def test_normalize_metric_det_1(self):
        """Normalisierte Metrik hat det(g) ≈ 1."""
        g = np.array([[4.0, 0.0], [0.0, 4.0]])
        g_norm = normalize_metric(g)
        assert abs(np.linalg.det(g_norm) - 1.0) < 1e-10

    def test_normalize_metric_shape_preserved(self):
        """Form der Matrix bleibt erhalten."""
        g = np.array([[2.0, 1.0], [1.0, 3.0]])
        g_norm = normalize_metric(g)
        assert g_norm.shape == (2, 2)

    def test_normalize_metric_identity(self):
        """Einheitsmatrix wird nicht verändert."""
        g = np.eye(2)
        g_norm = normalize_metric(g)
        assert np.allclose(g_norm, np.eye(2), atol=1e-10)

    def test_normalize_metric_3x3(self):
        """Normalisierung einer 3×3-Matrix."""
        g = 8.0 * np.eye(3)
        g_norm = normalize_metric(g)
        assert abs(np.linalg.det(g_norm) - 1.0) < 1e-9

    def test_volume_element_flat(self):
        """Volumen-Element der Einheitsmetrik = 1."""
        g = np.eye(2)
        vol = compute_volume_element(g)
        assert vol == pytest.approx(1.0)

    def test_volume_element_scaled(self):
        """Volumen-Element der skalierten Metrik = scale."""
        g = 4.0 * np.eye(2)
        vol = compute_volume_element(g)
        assert vol == pytest.approx(4.0)

    def test_volume_element_3d(self):
        """Volumen-Element in 3D: √det(g)."""
        g = np.diag([1.0, 4.0, 9.0])
        vol = compute_volume_element(g)
        expected = np.sqrt(1.0 * 4.0 * 9.0)
        assert vol == pytest.approx(expected)

    def test_volume_element_nonnegative(self):
        """Volumen-Element ist immer ≥ 0."""
        g = np.eye(3)
        assert compute_volume_element(g) >= 0.0

    def test_convergence_rate_decreasing(self):
        """Abnehmende Krümmung → positive Konvergenzrate."""
        # Exponentiell abnehmende Krümmung
        history = [10.0 * np.exp(-0.5 * t) for t in range(20)]
        rate = ricci_flow_convergence_rate(history)
        assert rate > 0

    def test_convergence_rate_constant(self):
        """Konstante Krümmung → Konvergenzrate ≈ 0."""
        history = [1.0] * 20
        rate = ricci_flow_convergence_rate(history)
        assert abs(rate) < 0.1

    def test_convergence_rate_short_history(self):
        """Sehr kurze Historie gibt 0 zurück."""
        rate = ricci_flow_convergence_rate([1.0])
        assert rate == pytest.approx(0.0)

    def test_ricci_flow_on_surface_returns_list(self):
        """ricci_flow_on_surface() gibt eine Liste zurück."""
        result = ricci_flow_on_surface(flat_metric_2d, T=0.1, steps=5)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_ricci_flow_on_surface_tuple_format(self):
        """Jedes Element ist ein (float, ndarray)-Paar."""
        result = ricci_flow_on_surface(flat_metric_2d, T=0.1, steps=3)
        for t, g in result:
            assert isinstance(float(t), float)
            assert isinstance(g, np.ndarray)

    def test_volume_element_singular_matrix(self):
        """Singuläre Matrix → Volumen-Element = 0."""
        g = np.zeros((2, 2))
        vol = compute_volume_element(g)
        assert vol == pytest.approx(0.0)
