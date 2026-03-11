"""
@file test_tensor_geometry.py
@brief Tests für das Modul tensor_geometry.py (Tensor- und Differentialgeometrie).
@description
    Umfassende Test-Suite für alle Klassen und Funktionen des Tensor-Moduls.
    Testet:
    - Tensor-Klasse: Arithmetik, Kontraktion, Symmetrie, äußeres Produkt
    - MetricTensor-Klasse: Inneres Produkt, Norm, Determinante, Volumenelement
    - Christoffel-Symbole: Flacher Raum ≈ 0, Sphäre bekannte Werte
    - Krümmungstensoren: Riemann, Ricci-Tensor, Ricci-Skalar
    - Gaußsche Krümmung: Sphäre K=1/r², Ebene K=0
    - Geodäten: Rückgabe-Struktur, Energieerhaltung
    - Mannigfaltigkeiten: Sphäre, Torus, Hyperbolische Ebene, Sattelform
    - Differentialformen: Wedge-Produkt, Hodge-Stern
    - Hilfsfunktionen: Levi-Civita, Kronecker-Delta, Paralleltransport

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import numpy as np
import pytest
import sys
import os
import time

# Projektpfad einbinden
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from tensor_geometry import (
    Tensor,
    MetricTensor,
    christoffel_symbols,
    riemann_tensor,
    ricci_tensor,
    ricci_scalar,
    gaussian_curvature,
    geodesic_equation,
    geodesic_distance,
    sphere_metric,
    torus_metric,
    hyperbolic_plane_metric,
    saddle_metric,
    flat_metric,
    wedge_product,
    exterior_derivative,
    hodge_star,
    lie_derivative_vector,
    einstein_tensor,
    schwarzschild_metric,
    check_vacuum_solution,
    levi_civita_symbol,
    kronecker_delta,
    parallel_transport,
    clear_christoffel_cache,
    get_christoffel_cache_size,
    _christoffel_cache,
)


# ===========================================================================
# Tests: Tensor-Klasse
# ===========================================================================

class TestTensor:
    """Tests für die Tensor-Klasse."""

    def test_tensor_creation_scalar(self):
        """Skalar als (0,0)-Tensor erstellen."""
        t = Tensor(5.0, contravariant=0, covariant=0)
        assert float(t.data) == pytest.approx(5.0)
        assert t.rank == 0

    def test_tensor_creation_vector(self):
        """Vektor als (1,0)-Tensor."""
        v = Tensor([1.0, 2.0, 3.0], contravariant=1, covariant=0)
        assert v.shape == (3,)
        assert v.rank == 1
        assert v.contravariant == 1
        assert v.covariant == 0

    def test_tensor_creation_matrix(self):
        """Matrix als (1,1)-Tensor."""
        T = Tensor([[1, 2], [3, 4]], contravariant=1, covariant=1)
        assert T.shape == (2, 2)
        assert T.rank == 2

    def test_tensor_addition(self):
        """Tensoraddition zweier (0,2)-Tensoren."""
        A = Tensor([[1, 0], [0, 1]], covariant=2)
        B = Tensor([[2, 1], [1, 2]], covariant=2)
        C = A + B
        expected = np.array([[3, 1], [1, 3]])
        np.testing.assert_allclose(C.data, expected)
        assert C.covariant == 2

    def test_tensor_addition_wrong_shape(self):
        """Tensoraddition mit inkompatiblen Formen wirft Fehler."""
        A = Tensor([[1, 0], [0, 1]], covariant=2)
        B = Tensor([1, 2, 3], contravariant=1)
        with pytest.raises(ValueError):
            _ = A + B

    def test_tensor_subtraction(self):
        """Tensorsubtraktion."""
        A = Tensor([[3, 1], [1, 3]], covariant=2)
        B = Tensor([[1, 0], [0, 1]], covariant=2)
        C = A - B
        expected = np.array([[2, 1], [1, 2]])
        np.testing.assert_allclose(C.data, expected)

    def test_tensor_scalar_multiplication(self):
        """Skalarmultiplikation."""
        T = Tensor([1.0, 2.0, 3.0], contravariant=1)
        result = T * 3.0
        np.testing.assert_allclose(result.data, [3.0, 6.0, 9.0])

    def test_tensor_rmul(self):
        """Rechtsseitige Skalarmultiplikation."""
        T = Tensor([1.0, 2.0], contravariant=1)
        result = 2.0 * T
        np.testing.assert_allclose(result.data, [2.0, 4.0])

    def test_outer_product_rank(self):
        """Äußeres Produkt: Rang = Summe der Einzel-Ränge."""
        u = Tensor([1.0, 2.0], contravariant=1)       # Rang (1,0)
        v = Tensor([3.0, 4.0], contravariant=1)       # Rang (1,0)
        T = u.outer_product(v)
        # Ergebnis: Rang (2,0)
        assert T.contravariant == 2
        assert T.covariant == 0
        assert T.shape == (2, 2)

    def test_outer_product_values(self):
        """Äußeres Produkt: Werte korrekt."""
        u = Tensor([1.0, 2.0], contravariant=1)
        v = Tensor([3.0, 4.0], contravariant=1)
        T = u.outer_product(v)
        expected = np.array([[3.0, 4.0], [6.0, 8.0]])
        np.testing.assert_allclose(T.data, expected)

    def test_outer_product_mixed_rank(self):
        """Äußeres Produkt mit gemischten Rängen."""
        A = Tensor([[1, 0], [0, 1]], contravariant=1, covariant=1)  # (1,1)
        v = Tensor([1.0, 0.0], contravariant=1)                      # (1,0)
        T = A.outer_product(v)
        assert T.contravariant == 2
        assert T.covariant == 1

    def test_contract_trace(self):
        """Kontraktion eines (1,1)-Tensors ergibt die Spur."""
        T = Tensor([[1, 2], [3, 4]], contravariant=1, covariant=1)
        result = T.contract(0, 0)
        # Spur: 1 + 4 = 5
        assert float(result.data) == pytest.approx(5.0)

    def test_contract_rank_reduction(self):
        """Kontraktion reduziert den Rang um (1,1)."""
        T = Tensor([[1, 2], [3, 4]], contravariant=1, covariant=1)
        result = T.contract(0, 0)
        assert result.rank == 0

    def test_contract_no_index_error(self):
        """Kontraktion ohne passende Indizes wirft Fehler."""
        T = Tensor([1.0, 2.0], contravariant=1, covariant=0)
        with pytest.raises(ValueError):
            T.contract(0, 0)

    def test_trace_2x2(self):
        """Spur einer 2×2-Matrix."""
        T = Tensor([[1, 2], [3, 4]], contravariant=1, covariant=1)
        assert T.trace() == pytest.approx(5.0)

    def test_trace_identity(self):
        """Spur der 3×3-Einheitsmatrix = 3."""
        I = Tensor(np.eye(3), contravariant=1, covariant=1)
        assert I.trace() == pytest.approx(3.0)

    def test_trace_wrong_rank(self):
        """Spur für nicht-(1,1)-Tensor wirft Fehler."""
        T = Tensor([1, 2, 3], contravariant=1, covariant=0)
        with pytest.raises(ValueError):
            T.trace()

    def test_symmetrize(self):
        """Symmetrisierung: Ergebnis ist symmetrisch."""
        T = Tensor([[1, 2], [3, 4]], covariant=2)
        sym = T.symmetrize()
        assert sym.is_symmetric()
        # Wert: ½([[1,2],[3,4]] + [[1,3],[2,4]]) = [[1, 2.5], [2.5, 4]]
        np.testing.assert_allclose(sym.data, [[1.0, 2.5], [2.5, 4.0]])

    def test_antisymmetrize(self):
        """Antisymmetrisierung: Ergebnis ist antisymmetrisch."""
        T = Tensor([[0, 2], [-2, 0]], covariant=2)
        anti = T.antisymmetrize()
        assert anti.is_antisymmetric()

    def test_symmetrize_plus_antisymmetrize(self):
        """Zerlegung: sym + antisym = Original."""
        T = Tensor([[1, 3], [7, 2]], covariant=2)
        sym = T.symmetrize()
        anti = T.antisymmetrize()
        reconstructed = sym.data + anti.data
        np.testing.assert_allclose(reconstructed, T.data, atol=1e-12)

    def test_is_symmetric_identity(self):
        """Einheitsmatrix als (0,2)-Tensor ist symmetrisch."""
        I = Tensor(np.eye(3), covariant=2)
        assert I.is_symmetric()

    def test_is_symmetric_asymmetric(self):
        """Nicht-symmetrischer Tensor: is_symmetric = False."""
        T = Tensor([[1, 2], [3, 4]], covariant=2)
        assert not T.is_symmetric()

    def test_is_antisymmetric(self):
        """Antisymmetrischer Tensor korrekt erkannt."""
        A = Tensor([[0, 1, -2], [-1, 0, 3], [2, -3, 0]], covariant=2)
        assert A.is_antisymmetric()

    def test_is_antisymmetric_false(self):
        """Symmetrischer Tensor ist nicht antisymmetrisch."""
        I = Tensor(np.eye(2), covariant=2)
        assert not I.is_antisymmetric()

    def test_repr(self):
        """__repr__ gibt String zurück."""
        T = Tensor([1.0, 2.0], contravariant=1)
        r = repr(T)
        assert "Tensor(1,0)" in r


# ===========================================================================
# Tests: MetricTensor-Klasse
# ===========================================================================

class TestMetricTensor:
    """Tests für die MetricTensor-Klasse."""

    def test_creation_euclidean(self):
        """Euklidische Metrik erstellen."""
        g = MetricTensor([[1, 0], [0, 1]])
        assert g.n == 2
        np.testing.assert_allclose(g.g, np.eye(2))

    def test_creation_non_square_error(self):
        """Nicht-quadratische Matrix wirft Fehler."""
        with pytest.raises(ValueError):
            MetricTensor([[1, 0, 0], [0, 1, 0]])

    def test_determinant_identity(self):
        """Determinante der Einheitsmatrix = 1."""
        g = MetricTensor(np.eye(3))
        assert g.determinant() == pytest.approx(1.0)

    def test_determinant_scaled(self):
        """Determinante der 2r²-Diagonalmetrik."""
        g = MetricTensor([[4.0, 0], [0, 4.0]])
        assert g.determinant() == pytest.approx(16.0)

    def test_volume_element_identity(self):
        """Volumenelement der Einheitsmetrik = 1."""
        g = MetricTensor(np.eye(2))
        assert g.volume_element() == pytest.approx(1.0)

    def test_volume_element_scaled(self):
        """Volumenelement der skalierten Metrik."""
        g = MetricTensor([[9.0, 0], [0, 4.0]])
        # √(9·4) = √36 = 6
        assert g.volume_element() == pytest.approx(6.0)

    def test_inner_product_euclidean(self):
        """Euklidisches inneres Produkt = Standard-Skalarprodukt."""
        g = MetricTensor(np.eye(3))
        u = [1.0, 2.0, 3.0]
        v = [4.0, 5.0, 6.0]
        # Skalarprodukt: 1*4 + 2*5 + 3*6 = 32
        assert g.inner_product(u, v) == pytest.approx(32.0)

    def test_inner_product_orthogonal(self):
        """Orthogonale Vektoren haben inneres Produkt 0."""
        g = MetricTensor(np.eye(2))
        u = [1.0, 0.0]
        v = [0.0, 1.0]
        assert g.inner_product(u, v) == pytest.approx(0.0)

    def test_norm_euclidean(self):
        """Euklidische Norm eines Einheitsvektors = 1."""
        g = MetricTensor(np.eye(3))
        v = [1.0, 0.0, 0.0]
        assert g.norm(v) == pytest.approx(1.0)

    def test_norm_pythagorean(self):
        """Norm eines 3-4-5-Vektors = 5."""
        g = MetricTensor(np.eye(2))
        v = [3.0, 4.0]
        assert g.norm(v) == pytest.approx(5.0)

    def test_raise_index(self):
        """Indexhebung mit Einheitsmetrik ändert Tensor nicht."""
        g = MetricTensor(np.eye(3))
        T = Tensor([1.0, 2.0, 3.0], covariant=1)
        raised = g.raise_index(T, 0)
        np.testing.assert_allclose(raised.data, T.data, atol=1e-10)

    def test_lower_index(self):
        """Indexsenkung mit Einheitsmetrik ändert Tensor nicht."""
        g = MetricTensor(np.eye(3))
        T = Tensor([1.0, 2.0, 3.0], contravariant=1)
        lowered = g.lower_index(T, 0)
        np.testing.assert_allclose(lowered.data, T.data, atol=1e-10)

    def test_is_flat(self):
        """Konstante Metrik gilt immer als flach."""
        g = MetricTensor(np.eye(2))
        assert g.is_flat([0.0, 0.0])


# ===========================================================================
# Tests: Klassische Mannigfaltigkeiten
# ===========================================================================

class TestManifolds:
    """Tests für die vordefinierten Mannigfaltigkeiten."""

    def test_flat_metric_2d(self):
        """Flache 2D-Metrik gibt Einheitsmatrix zurück."""
        metric = flat_metric(2)
        g = metric([0.5, 0.3])
        np.testing.assert_allclose(g, np.eye(2))

    def test_flat_metric_3d(self):
        """Flache 3D-Metrik gibt 3×3-Einheitsmatrix zurück."""
        metric = flat_metric(3)
        g = metric([1.0, 2.0, 3.0])
        np.testing.assert_allclose(g, np.eye(3))

    def test_sphere_metric_shape(self):
        """Sphären-Metrik hat Form (2,2)."""
        metric = sphere_metric(r=1.0)
        g = metric([math.pi / 2, 0.0])
        assert g.shape == (2, 2)

    def test_sphere_metric_at_equator(self):
        """Sphären-Metrik am Äquator (θ=π/2): g = diag(r², r²)."""
        r = 2.0
        metric = sphere_metric(r=r)
        theta = math.pi / 2.0
        g = metric([theta, 0.0])
        expected = np.array([[r**2, 0], [0, r**2 * math.sin(theta)**2]])
        np.testing.assert_allclose(g, expected, atol=1e-12)

    def test_sphere_metric_at_pole(self):
        """Sphären-Metrik am Pol (θ=0): g_φφ = 0 (Koordinaten-Singularität)."""
        metric = sphere_metric(r=1.0)
        g = metric([0.0, 0.0])
        # g_φφ = sin²(0) = 0
        assert g[1, 1] == pytest.approx(0.0, abs=1e-12)

    def test_sphere_metric_symmetry(self):
        """Sphären-Metrik ist eine Diagonalmatrix."""
        metric = sphere_metric()
        g = metric([1.0, 0.5])
        assert g[0, 1] == pytest.approx(0.0)
        assert g[1, 0] == pytest.approx(0.0)

    def test_torus_metric_shape(self):
        """Torus-Metrik hat Form (2,2)."""
        metric = torus_metric(R=2.0, r=1.0)
        g = metric([0.0, 0.0])
        assert g.shape == (2, 2)

    def test_torus_metric_diagonal(self):
        """Torus-Metrik ist diagonal."""
        metric = torus_metric(R=2.0, r=1.0)
        g = metric([0.5, 1.0])
        assert g[0, 1] == pytest.approx(0.0)
        assert g[1, 0] == pytest.approx(0.0)

    def test_torus_metric_outer(self):
        """Torus-Metrik außen (θ=0): g_φφ = (R+r)²."""
        R, r = 3.0, 1.0
        metric = torus_metric(R=R, r=r)
        g = metric([0.0, 0.0])
        assert g[1, 1] == pytest.approx((R + r)**2)

    def test_torus_metric_inner(self):
        """Torus-Metrik innen (θ=π): g_φφ = (R-r)²."""
        R, r = 3.0, 1.0
        metric = torus_metric(R=R, r=r)
        g = metric([math.pi, 0.0])
        assert g[1, 1] == pytest.approx((R - r)**2, rel=1e-10)

    def test_hyperbolic_metric_shape(self):
        """Hyperbolische Metrik hat Form (2,2)."""
        metric = hyperbolic_plane_metric()
        g = metric([0.0, 1.0])
        assert g.shape == (2, 2)

    def test_hyperbolic_metric_diagonal(self):
        """Hyperbolische Metrik ist diagonal."""
        metric = hyperbolic_plane_metric()
        g = metric([0.0, 2.0])
        assert g[0, 1] == pytest.approx(0.0)

    def test_hyperbolic_metric_values(self):
        """Hyperbolische Metrik bei y=2: g = diag(1/4, 1/4)."""
        metric = hyperbolic_plane_metric()
        g = metric([0.0, 2.0])
        np.testing.assert_allclose(g, np.array([[0.25, 0], [0, 0.25]]))

    def test_saddle_metric_shape(self):
        """Sattelflächen-Metrik hat Form (2,2)."""
        metric = saddle_metric()
        g = metric([1.0, 1.0])
        assert g.shape == (2, 2)

    def test_saddle_metric_at_origin(self):
        """Sattelflächen-Metrik am Ursprung = Einheitsmatrix."""
        metric = saddle_metric()
        g = metric([0.0, 0.0])
        np.testing.assert_allclose(g, np.eye(2))

    def test_saddle_metric_symmetric(self):
        """Sattelflächen-Metrik ist symmetrisch."""
        metric = saddle_metric()
        g = metric([1.0, 2.0])
        np.testing.assert_allclose(g, g.T)


# ===========================================================================
# Tests: Christoffel-Symbole
# ===========================================================================

class TestChristoffelSymbols:
    """Tests für die Christoffel-Symbole."""

    def test_flat_metric_christoffel_near_zero(self):
        """Flache Metrik: alle Christoffel-Symbole ≈ 0."""
        metric = flat_metric(2)
        Gamma = christoffel_symbols(metric, [0.5, 0.5])
        np.testing.assert_allclose(Gamma, np.zeros((2, 2, 2)), atol=1e-8)

    def test_flat_metric_3d_christoffel_near_zero(self):
        """Flache 3D-Metrik: alle Christoffel-Symbole ≈ 0."""
        metric = flat_metric(3)
        Gamma = christoffel_symbols(metric, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(Gamma, np.zeros((3, 3, 3)), atol=1e-8)

    def test_christoffel_symmetry(self):
        """Christoffel-Symbole sind symmetrisch: Γ^k_{ij} = Γ^k_{ji}."""
        metric = sphere_metric(r=1.0)
        Gamma = christoffel_symbols(metric, [math.pi / 3, 0.5])
        n = Gamma.shape[0]
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    assert Gamma[k, i, j] == pytest.approx(Gamma[k, j, i], abs=1e-6)

    def test_christoffel_shape(self):
        """Christoffel-Symbole haben Form (n,n,n)."""
        metric = sphere_metric()
        Gamma = christoffel_symbols(metric, [1.0, 0.5])
        assert Gamma.shape == (2, 2, 2)

    def test_sphere_christoffel_theta_phi_phi(self):
        """
        Sphäre: Γ^θ_{φφ} = -sin(θ)cos(θ).
        Indizes (0=θ, 1=φ): Gamma[0,1,1] ≈ -sin(θ)cos(θ).
        """
        theta = math.pi / 3.0
        metric = sphere_metric(r=1.0)
        Gamma = christoffel_symbols(metric, [theta, 0.0])
        expected = -math.sin(theta) * math.cos(theta)
        assert Gamma[0, 1, 1] == pytest.approx(expected, abs=1e-5)

    def test_sphere_christoffel_phi_theta_phi(self):
        """
        Sphäre: Γ^φ_{θφ} = cos(θ)/sin(θ) = cot(θ).
        Indizes: Gamma[1,0,1] ≈ cot(θ).
        """
        theta = math.pi / 3.0
        metric = sphere_metric(r=1.0)
        Gamma = christoffel_symbols(metric, [theta, 0.0])
        expected = math.cos(theta) / math.sin(theta)
        assert Gamma[1, 0, 1] == pytest.approx(expected, abs=1e-5)


# ===========================================================================
# Tests: Krümmungstensoren
# ===========================================================================

class TestCurvature:
    """Tests für Riemann-Tensor, Ricci-Tensor und Ricci-Skalar."""

    def test_riemann_flat_near_zero(self):
        """Flacher Raum: Riemann-Tensor ≈ 0."""
        metric = flat_metric(2)
        R = riemann_tensor(metric, [0.5, 0.5])
        np.testing.assert_allclose(R, np.zeros((2, 2, 2, 2)), atol=1e-5)

    def test_riemann_tensor_shape(self):
        """Riemann-Tensor hat Form (n,n,n,n)."""
        metric = sphere_metric()
        R = riemann_tensor(metric, [math.pi / 2, 0.0])
        assert R.shape == (2, 2, 2, 2)

    def test_riemann_antisymmetry_ij(self):
        """Riemann-Tensor: R^l_{kij} = -R^l_{kji}."""
        metric = sphere_metric(r=1.0)
        R = riemann_tensor(metric, [math.pi / 2, 0.5])
        n = 2
        for l in range(n):
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        assert R[l, k, i, j] == pytest.approx(-R[l, k, j, i], abs=1e-4)

    def test_ricci_tensor_flat_near_zero(self):
        """Flacher Raum: Ricci-Tensor ≈ 0."""
        metric = flat_metric(2)
        Ric = ricci_tensor(metric, [0.5, 0.5])
        np.testing.assert_allclose(Ric, np.zeros((2, 2)), atol=1e-5)

    def test_ricci_tensor_shape(self):
        """Ricci-Tensor hat Form (n,n)."""
        metric = sphere_metric()
        Ric = ricci_tensor(metric, [math.pi / 2, 0.0])
        assert Ric.shape == (2, 2)

    def test_ricci_tensor_symmetry(self):
        """Ricci-Tensor ist symmetrisch."""
        metric = sphere_metric(r=1.0)
        Ric = ricci_tensor(metric, [math.pi / 2, 0.0])
        np.testing.assert_allclose(Ric, Ric.T, atol=1e-5)

    def test_ricci_scalar_flat(self):
        """Flacher Raum: Ricci-Skalar ≈ 0."""
        metric = flat_metric(2)
        R = ricci_scalar(metric, [0.5, 0.5])
        assert R == pytest.approx(0.0, abs=1e-5)

    def test_ricci_scalar_sphere_radius_1(self):
        """2-Sphäre mit r=1: Ricci-Skalar R ≈ 2.0 am Äquator."""
        metric = sphere_metric(r=1.0)
        # Äquator: θ = π/2 (Singularität vermeiden)
        R = ricci_scalar(metric, [math.pi / 2, 0.0])
        # Erwarteter Wert: R = 2/r² = 2
        assert R == pytest.approx(2.0, rel=0.1)

    def test_ricci_scalar_sphere_radius_2(self):
        """2-Sphäre mit r=2: Ricci-Skalar R ≈ 0.5."""
        metric = sphere_metric(r=2.0)
        R = ricci_scalar(metric, [math.pi / 2, 0.0])
        # Erwarteter Wert: R = 2/r² = 2/4 = 0.5
        assert R == pytest.approx(0.5, rel=0.15)


# ===========================================================================
# Tests: Gaußsche Krümmung
# ===========================================================================

class TestGaussianCurvature:
    """Tests für die Gaußsche Krümmung."""

    def test_flat_metric_gaussian_zero(self):
        """Flache Metrik: K ≈ 0."""
        metric = flat_metric(2)
        K = gaussian_curvature(metric, [0.5, 0.5])
        assert K == pytest.approx(0.0, abs=1e-5)

    def test_sphere_r1_gaussian_curvature(self):
        """2-Sphäre mit r=1: K ≈ 1 am Äquator."""
        metric = sphere_metric(r=1.0)
        K = gaussian_curvature(metric, [math.pi / 2, 0.0])
        # K = 1/r² = 1
        assert K == pytest.approx(1.0, rel=0.15)

    def test_sphere_r2_gaussian_curvature(self):
        """2-Sphäre mit r=2: K ≈ 0.25."""
        metric = sphere_metric(r=2.0)
        K = gaussian_curvature(metric, [math.pi / 2, 0.0])
        # K = 1/r² = 0.25
        assert K == pytest.approx(0.25, rel=0.15)

    def test_gaussian_wrong_dimension(self):
        """Gaußsche Krümmung für 3D-Metrik wirft Fehler."""
        metric = flat_metric(3)
        with pytest.raises(ValueError):
            gaussian_curvature(metric, [0.0, 0.0, 0.0])

    def test_saddle_gaussian_negative(self):
        """Sattelfläche: K < 0 am Ursprung (hyperbolisch)."""
        metric = saddle_metric()
        K = gaussian_curvature(metric, [0.01, 0.01])
        # K = -1/(1+x²+y²)² ≈ -1 nahe Ursprung
        assert K < 0


# ===========================================================================
# Tests: Geodäten
# ===========================================================================

class TestGeodesics:
    """Tests für Geodätenberechnung."""

    def test_geodesic_equation_returns_dict(self):
        """geodesic_equation gibt Dict mit 'trajectory', 'times', 'length' zurück."""
        metric = flat_metric(2)
        result = geodesic_equation(metric, [0.0, 0.0], [1.0, 0.0], t_max=0.1, n_steps=10)
        assert isinstance(result, dict)
        assert "trajectory" in result
        assert "times" in result
        assert "length" in result

    def test_geodesic_trajectory_shape(self):
        """Trajektorie hat Form (n_steps+1, dim)."""
        metric = flat_metric(2)
        n_steps = 20
        result = geodesic_equation(
            metric, [0.0, 0.0], [1.0, 0.0], t_max=0.5, n_steps=n_steps
        )
        assert result["trajectory"].shape == (n_steps + 1, 2)

    def test_geodesic_flat_straight_line(self):
        """Geodäte in flachem Raum = gerade Linie."""
        metric = flat_metric(2)
        v0 = [1.0, 0.0]
        t_max = 1.0
        result = geodesic_equation(metric, [0.0, 0.0], v0, t_max=t_max, n_steps=100)
        traj = result["trajectory"]
        # y-Koordinate sollte nahezu 0 bleiben
        np.testing.assert_allclose(traj[:, 1], 0.0, atol=1e-8)
        # x-Koordinate: linear von 0 bis ~t_max
        assert traj[-1, 0] == pytest.approx(t_max, rel=0.01)

    def test_geodesic_times_monotone(self):
        """Zeitpunkte sind monoton wachsend."""
        metric = flat_metric(2)
        result = geodesic_equation(metric, [0.0, 0.0], [1.0, 1.0], t_max=1.0, n_steps=50)
        times = result["times"]
        assert all(times[i] <= times[i + 1] for i in range(len(times) - 1))

    def test_geodesic_start_point(self):
        """Startpunkt der Trajektorie = x0."""
        metric = flat_metric(2)
        x0 = [3.0, -1.0]
        result = geodesic_equation(metric, x0, [1.0, 0.0], t_max=0.5, n_steps=10)
        np.testing.assert_allclose(result["trajectory"][0], x0)

    def test_geodesic_length_positive(self):
        """Geodätische Länge ist nicht-negativ."""
        metric = sphere_metric(r=1.0)
        result = geodesic_equation(
            metric, [math.pi / 2, 0.0], [0.0, 1.0], t_max=0.5, n_steps=50
        )
        assert result["length"] >= 0.0

    def test_geodesic_distance_same_point(self):
        """Geodätischer Abstand vom Punkt zu sich selbst ≈ 0."""
        metric = flat_metric(2)
        d = geodesic_distance(metric, [1.0, 1.0], [1.0, 1.0])
        assert d == pytest.approx(0.0, abs=1e-10)

    def test_geodesic_distance_flat(self):
        """Geodätischer Abstand in flachem Raum = euklidische Distanz."""
        metric = flat_metric(2)
        x0 = [0.0, 0.0]
        x1 = [3.0, 4.0]
        d = geodesic_distance(metric, x0, x1, n_steps=200)
        # Euklidische Distanz: √(9+16) = 5
        assert d == pytest.approx(5.0, rel=0.001)

    def test_geodesic_distance_symmetry(self):
        """Geodätischer Abstand ist symmetrisch: d(x0,x1) = d(x1,x0)."""
        metric = flat_metric(2)
        d1 = geodesic_distance(metric, [1.0, 2.0], [4.0, 6.0])
        d2 = geodesic_distance(metric, [4.0, 6.0], [1.0, 2.0])
        assert d1 == pytest.approx(d2, rel=1e-10)


# ===========================================================================
# Tests: Differentialformen
# ===========================================================================

class TestDifferentialForms:
    """Tests für Wedge-Produkt, Hodge-Stern und Lie-Ableitung."""

    def test_wedge_product_antisymmetric(self):
        """Wedge-Produkt ist antisymmetrisch: α∧β = -(β∧α)."""
        alpha = np.array([1.0, 2.0, 3.0])
        beta = np.array([4.0, 5.0, 6.0])
        W1 = wedge_product(alpha, beta)
        W2 = wedge_product(beta, alpha)
        np.testing.assert_allclose(W1, -W2)

    def test_wedge_product_diagonal_zero(self):
        """Wedge-Produkt: Diagonale ist immer 0 (α∧α = 0)."""
        alpha = np.array([1.0, 2.0, 3.0])
        W = wedge_product(alpha, alpha)
        np.testing.assert_allclose(np.diag(W), 0.0, atol=1e-12)

    def test_wedge_product_shape(self):
        """Wedge-Produkt zweier 1-Formen ergibt n×n-Matrix."""
        alpha = np.array([1.0, 2.0])
        beta = np.array([3.0, 4.0])
        W = wedge_product(alpha, beta)
        assert W.shape == (2, 2)

    def test_wedge_product_values(self):
        """Wedge-Produkt: konkrete Werte für 2D."""
        alpha = np.array([1.0, 0.0])
        beta = np.array([0.0, 1.0])
        W = wedge_product(alpha, beta)
        # (α∧β)_{01} = α_0 β_1 - α_1 β_0 = 1·1 - 0·0 = 1
        assert W[0, 1] == pytest.approx(1.0)
        assert W[1, 0] == pytest.approx(-1.0)

    def test_exterior_derivative_shape(self):
        """Äußere Ableitung einer 1-Form ergibt n×n-Matrix."""
        form = np.array([1.0, 2.0])
        result = exterior_derivative(form, [0.0, 0.0])
        assert result.shape == (2, 2)

    def test_hodge_star_2d_1form(self):
        """Hodge-Stern einer 1-Form in 2D dreht den Vektor."""
        g = MetricTensor(np.eye(2))
        alpha = np.array([1.0, 0.0])
        result = hodge_star(alpha, g)
        # *(1, 0) = (0, 1) · √|det(g)| = (0, 1)
        assert result.shape == (2,)
        np.testing.assert_allclose(result, [0.0, 1.0], atol=1e-12)

    def test_hodge_star_2d_1form_y_direction(self):
        """Hodge-Stern von (0,1) in 2D = (-1, 0)."""
        g = MetricTensor(np.eye(2))
        alpha = np.array([0.0, 1.0])
        result = hodge_star(alpha, g)
        np.testing.assert_allclose(result, [-1.0, 0.0], atol=1e-12)

    def test_lie_derivative_linear(self):
        """Lie-Ableitung von f(x,y) = x+y in Richtung (1,0) = 1."""
        def f(point):
            return point[0] + point[1]

        result = lie_derivative_vector(f, [1.0, 0.0], [0.5, 0.5])
        assert result == pytest.approx(1.0, abs=1e-8)

    def test_lie_derivative_constant(self):
        """Lie-Ableitung einer Konstanten = 0."""
        def f(point):
            return 5.0

        result = lie_derivative_vector(f, [1.0, 1.0], [0.0, 0.0])
        assert result == pytest.approx(0.0, abs=1e-8)

    def test_lie_derivative_quadratic(self):
        """Lie-Ableitung von f(x,y) = x² in Richtung (1,0) bei x=2: 2x = 4."""
        def f(point):
            return point[0] ** 2

        result = lie_derivative_vector(f, [1.0, 0.0], [2.0, 0.0])
        assert result == pytest.approx(4.0, rel=1e-5)


# ===========================================================================
# Tests: Einstein-Gleichungen
# ===========================================================================

class TestEinsteinEquations:
    """Tests für Einstein-Tensor und Schwarzschild-Metrik."""

    def test_einstein_tensor_flat_near_zero(self):
        """Flacher Raum: Einstein-Tensor ≈ 0."""
        metric = flat_metric(2)
        G = einstein_tensor(metric, [0.5, 0.5])
        np.testing.assert_allclose(G, np.zeros((2, 2)), atol=1e-5)

    def test_einstein_tensor_shape(self):
        """Einstein-Tensor hat Form (n,n)."""
        metric = sphere_metric(r=1.0)
        G = einstein_tensor(metric, [math.pi / 2, 0.0])
        assert G.shape == (2, 2)

    def test_schwarzschild_metric_shape(self):
        """Schwarzschild-Metrik hat Form (2,2)."""
        metric = schwarzschild_metric(M=1.0)
        g = metric([0.0, 10.0])  # weit entfernt
        assert g.shape == (2, 2)

    def test_schwarzschild_metric_far_field(self):
        """Schwarzschild-Metrik im Fernfeld (r → ∞): g ≈ Minkowski."""
        metric = schwarzschild_metric(M=1.0, c=1.0, G=1.0)
        # r_s = 2; bei r=10000 >> r_s: nahezu Minkowski
        g = metric([0.0, 10000.0])
        # g_tt ≈ -c² = -1, g_rr ≈ 1
        assert g[0, 0] == pytest.approx(-1.0, rel=1e-3)
        assert g[1, 1] == pytest.approx(1.0, rel=1e-3)

    def test_schwarzschild_metric_singularity_avoidance(self):
        """Schwarzschild-Metrik vermeidet Singularität bei r <= r_s."""
        metric = schwarzschild_metric(M=1.0)
        # r_s = 2; r = 1 < r_s: sollte keinen Fehler werfen
        g = metric([0.0, 1.0])
        assert g is not None
        assert not np.any(np.isnan(g))

    def test_check_vacuum_flat(self):
        """Flacher Raum erfüllt Vakuumgleichungen (G ≈ 0)."""
        metric = flat_metric(2)
        points = [[0.5, 0.5], [1.0, 1.0], [2.0, 2.0]]
        result = check_vacuum_solution(metric, points)
        assert result["is_vacuum"] is True
        assert result["points_checked"] == 3

    def test_check_vacuum_returns_dict(self):
        """check_vacuum_solution gibt Dict mit allen Schlüsseln zurück."""
        metric = flat_metric(2)
        result = check_vacuum_solution(metric, [[0.5, 0.5]])
        assert "is_vacuum" in result
        assert "max_deviation" in result
        assert "points_checked" in result


# ===========================================================================
# Tests: Hilfsfunktionen
# ===========================================================================

class TestHelpers:
    """Tests für Levi-Civita-Symbol, Kronecker-Delta und Paralleltransport."""

    def test_levi_civita_1d(self):
        """Levi-Civita-Symbol in 1D: ε_0 = 1."""
        eps = levi_civita_symbol(1)
        assert eps.shape == (1,)
        assert eps[0] == pytest.approx(1.0)

    def test_levi_civita_2d(self):
        """Levi-Civita-Symbol in 2D: ε_{01} = 1, ε_{10} = -1."""
        eps = levi_civita_symbol(2)
        assert eps.shape == (2, 2)
        assert eps[0, 1] == pytest.approx(1.0)
        assert eps[1, 0] == pytest.approx(-1.0)
        assert eps[0, 0] == pytest.approx(0.0)
        assert eps[1, 1] == pytest.approx(0.0)

    def test_levi_civita_3d_positive(self):
        """Levi-Civita in 3D: ε_{012} = 1 (gerade Permutation)."""
        eps = levi_civita_symbol(3)
        assert eps[0, 1, 2] == pytest.approx(1.0)

    def test_levi_civita_3d_negative(self):
        """Levi-Civita in 3D: ε_{021} = -1 (ungerade Permutation)."""
        eps = levi_civita_symbol(3)
        assert eps[0, 2, 1] == pytest.approx(-1.0)

    def test_levi_civita_3d_zero(self):
        """Levi-Civita in 3D: ε_{001} = 0 (wiederholter Index)."""
        eps = levi_civita_symbol(3)
        assert eps[0, 0, 1] == pytest.approx(0.0)

    def test_levi_civita_3d_all_permutations(self):
        """Levi-Civita in 3D: alle 6 Permutationen korrekt."""
        eps = levi_civita_symbol(3)
        assert eps[0, 1, 2] == pytest.approx(1.0)
        assert eps[1, 2, 0] == pytest.approx(1.0)
        assert eps[2, 0, 1] == pytest.approx(1.0)
        assert eps[0, 2, 1] == pytest.approx(-1.0)
        assert eps[2, 1, 0] == pytest.approx(-1.0)
        assert eps[1, 0, 2] == pytest.approx(-1.0)

    def test_kronecker_delta_2d(self):
        """Kronecker-Delta in 2D = Einheitsmatrix."""
        delta = kronecker_delta(2)
        np.testing.assert_allclose(delta, np.eye(2))

    def test_kronecker_delta_3d(self):
        """Kronecker-Delta in 3D = 3×3-Einheitsmatrix."""
        delta = kronecker_delta(3)
        np.testing.assert_allclose(delta, np.eye(3))

    def test_kronecker_delta_shape(self):
        """Kronecker-Delta hat Form (n,n)."""
        delta = kronecker_delta(5)
        assert delta.shape == (5, 5)

    def test_parallel_transport_flat_constant(self):
        """Paralleltransport in flachem Raum: Vektor bleibt konstant."""
        metric = flat_metric(2)
        v0 = [1.0, 0.0]
        # Kurve: gerade Linie
        curve = [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0]]
        transported = parallel_transport(metric, v0, curve)
        # Im flachen Raum ändert sich der Vektor nicht (Γ = 0)
        for v in transported:
            np.testing.assert_allclose(v, v0, atol=1e-8)

    def test_parallel_transport_output_length(self):
        """Paralleltransport gibt Vektor an jedem Kurvenpunkt zurück."""
        metric = flat_metric(2)
        curve = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]
        transported = parallel_transport(metric, [1.0, 0.0], curve)
        assert len(transported) == len(curve)

    def test_parallel_transport_start_preserved(self):
        """Paralleltransport: Anfangsvektor unverändert."""
        metric = flat_metric(2)
        v0 = [2.0, 3.0]
        curve = [[0.0, 0.0], [1.0, 1.0]]
        transported = parallel_transport(metric, v0, curve)
        np.testing.assert_allclose(transported[0], v0)


# ===========================================================================
# Edge-Cases und Sonderfälle
# ===========================================================================

class TestEdgeCases:
    """Tests für Randwerte und Sonderfälle."""

    def test_tensor_single_element(self):
        """Tensor mit einem Element."""
        T = Tensor([42.0], contravariant=1)
        assert T.shape == (1,)
        assert float(T.data[0]) == pytest.approx(42.0)

    def test_tensor_large_values(self):
        """Tensor mit sehr großen Werten."""
        T = Tensor([[1e15, 0], [0, 1e15]], covariant=2)
        assert T.is_symmetric()

    def test_tensor_negative_values(self):
        """Tensor mit negativen Werten."""
        T = Tensor([[-1.0, -2.0], [-3.0, -4.0]], contravariant=1, covariant=1)
        assert T.trace() == pytest.approx(-5.0)

    def test_metric_near_singular(self):
        """MetricTensor mit kleiner Determinante (fast singulär)."""
        g = MetricTensor([[1.0, 0.999], [0.999, 1.0]])
        assert abs(g.determinant()) > 1e-10  # Nicht exakt singulär

    def test_geodesic_zero_velocity(self):
        """Geodäte mit Null-Anfangsgeschwindigkeit bleibt am Startpunkt."""
        metric = flat_metric(2)
        x0 = [1.0, 2.0]
        result = geodesic_equation(metric, x0, [0.0, 0.0], t_max=1.0, n_steps=20)
        # Mit v=0 bleibt der Punkt (nahezu) stationär
        np.testing.assert_allclose(result["trajectory"][-1], x0, atol=1e-10)

    def test_sphere_metric_at_half_pi(self):
        """Sphäre bei θ=π/2: Metrik wohldefiniert."""
        metric = sphere_metric(r=1.0)
        g = metric([math.pi / 2, math.pi / 4])
        assert not np.any(np.isnan(g))
        assert g[0, 0] == pytest.approx(1.0)
        assert g[1, 1] == pytest.approx(1.0)

    def test_levi_civita_4d_shape(self):
        """Levi-Civita in 4D hat Form (4,4,4,4)."""
        eps = levi_civita_symbol(4)
        assert eps.shape == (4, 4, 4, 4)

    def test_levi_civita_4d_positive(self):
        """Levi-Civita in 4D: ε_{0123} = 1."""
        eps = levi_civita_symbol(4)
        assert eps[0, 1, 2, 3] == pytest.approx(1.0)

    def test_outer_product_rank_sum(self):
        """Äußeres Produkt: Rang ist Summe der Einzel-Ränge (allgemein)."""
        A = Tensor([[1, 0], [0, 1]], contravariant=1, covariant=1)  # (1,1)
        B = Tensor([[1, 0], [0, 1]], covariant=2)                   # (0,2)
        C = A.outer_product(B)
        assert C.contravariant == 1
        assert C.covariant == 3
        assert C.rank == 4

    def test_geodesic_distance_triangle_inequality(self):
        """Geodätischer Abstand erfüllt Dreiecksungleichung."""
        metric = flat_metric(2)
        a = [0.0, 0.0]
        b = [3.0, 0.0]
        c = [0.0, 4.0]
        d_ab = geodesic_distance(metric, a, b)
        d_bc = geodesic_distance(metric, b, c)
        d_ac = geodesic_distance(metric, a, c)
        # Dreiecksungleichung: d(a,c) ≤ d(a,b) + d(b,c)
        assert d_ac <= d_ab + d_bc + 1e-10

    def test_schwarzschild_radius_natural_units(self):
        """Schwarzschild-Radius in natürlichen Einheiten: r_s = 2M."""
        # In nat. Einheiten (G=1, c=1): r_s = 2GM/c² = 2M
        M = 3.0
        metric = schwarzschild_metric(M=M, c=1.0, G=1.0)
        # Weit außen (r >> r_s = 6): g_tt ≈ -1
        g_far = metric([0.0, 1000.0])
        assert g_far[0, 0] == pytest.approx(-1.0, rel=0.01)


# ===========================================================================
# Tests: Christoffel-Symbol-Caching (Build 51)
# ===========================================================================

class TestChristoffelCaching:
    """Tests für das Caching der Christoffel-Symbole."""

    def setup_method(self):
        """Cache vor jedem Test leeren für Isolation."""
        clear_christoffel_cache()

    def test_clear_cache_returns_count(self):
        """Test: clear_christoffel_cache() gibt die Anzahl gelöschter Einträge zurück."""
        # Einen Eintrag in den Cache schreiben
        metric = flat_metric(2)
        christoffel_symbols(metric, [0.0, 0.0])
        count = clear_christoffel_cache()
        assert count >= 1

    def test_get_cache_size_empty(self):
        """Test: get_christoffel_cache_size() gibt 0 nach dem Leeren zurück."""
        clear_christoffel_cache()
        assert get_christoffel_cache_size() == 0

    def test_cache_hit_on_repeated_call(self):
        """Test: Wiederholter Aufruf mit identischen Parametern trifft den Cache."""
        metric = flat_metric(2)
        point = [1.0, 2.0]

        clear_christoffel_cache()
        assert get_christoffel_cache_size() == 0

        # Erster Aufruf: Cache miss
        christoffel_symbols(metric, point)
        assert get_christoffel_cache_size() == 1

        # Zweiter Aufruf: Cache hit (Größe bleibt gleich)
        christoffel_symbols(metric, point)
        assert get_christoffel_cache_size() == 1

    def test_cache_grows_for_different_points(self):
        """Test: Cache wächst bei verschiedenen Punkten."""
        metric = flat_metric(2)
        clear_christoffel_cache()

        christoffel_symbols(metric, [0.0, 0.0])
        christoffel_symbols(metric, [1.0, 0.0])
        christoffel_symbols(metric, [0.0, 1.0])

        assert get_christoffel_cache_size() == 3

    def test_cache_result_is_same_object(self):
        """Test: Gecachtes Ergebnis ist dasselbe numpy-Array-Objekt."""
        metric = flat_metric(2)
        point = [0.5, 0.5]
        clear_christoffel_cache()

        result1 = christoffel_symbols(metric, point)
        result2 = christoffel_symbols(metric, point)

        # Exakt dasselbe Objekt (nicht nur gleiche Werte)
        assert result1 is result2

    def test_cache_speedup(self):
        """Test: Gecachter Aufruf ist schneller als erster Aufruf."""
        metric = sphere_metric()
        point = [1.0471, 0.5]  # ≈ π/3, 0.5
        clear_christoffel_cache()

        # Erster Aufruf (Berechnung)
        start1 = time.perf_counter()
        christoffel_symbols(metric, point)
        t1 = time.perf_counter() - start1

        # Zweiter Aufruf (Cache)
        start2 = time.perf_counter()
        christoffel_symbols(metric, point)
        t2 = time.perf_counter() - start2

        # Cache sollte deutlich schneller sein
        assert t2 < t1 or t2 < 1e-5  # Gecacht < 10 Mikrosekunden

    def test_cache_key_includes_step_size(self):
        """Test: Verschiedene Schrittweiten h erzeugen verschiedene Cache-Einträge."""
        metric = flat_metric(2)
        point = [1.0, 1.0]
        clear_christoffel_cache()

        christoffel_symbols(metric, point, h=1e-5)
        christoffel_symbols(metric, point, h=1e-6)

        # Zwei verschiedene h-Werte → zwei Cache-Einträge
        assert get_christoffel_cache_size() == 2
