r"""
@file test_geometric_flows.py
@brief Test-Suite für geometric_flows.py — Ricci-Fluss und Poincaré-Vermutung.
@description
    Umfassende Tests für alle Klassen in geometric_flows.py:
    - RiemannianMetric: Christoffel-Symbole, Riemann- und Ricci-Tensor, Ricci-Skalar
    - ThreeSphere: Metriktensor, topologische Invarianten
    - RicciFlow: Evolutionsgleichungen, Euler- und RK4-Schritte
    - RicciFlowSingularity: Aufblasraten, kanonische Nachbarschaften, κ-Nichtkollaps
    - PerelmanEntropy: W- und F-Funktionale, Monotoniecheck
    - RicciFlowWithSurgery: Singularitätserkennung, Chirurgie
    - PoincaréConjecture: Theorem-Aussagen, acht Geometrien, Beweis-Outline
    - GeometrizationConjecture: Klassifikation, JSJ-Zerlegung
    - UniformizationTheorem: Flächen-Klassifikation, Gauss-Bonnet, Ricci-Fluss 2D

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import numpy as np
import pytest

# Zu testende Module
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from geometric_flows import (
    RiemannianMetric,
    ThreeSphere,
    RicciFlow,
    RicciFlowSingularity,
    PerelmanEntropy,
    RicciFlowWithSurgery,
    PoincaréConjecture,
    GeometrizationConjecture,
    UniformizationTheorem,
    _finite_diff_gradient,
    _finite_diff_hessian,
)


# ===========================================================================
# HILFSFUNKTIONS-TESTS
# ===========================================================================

class TestFiniteDiff:
    r"""Tests für numerische Differentiationshelfer."""

    def test_gradient_linear(self):
        r"""Gradient einer linearen Funktion f(x) = 2x₀ + 3x₁."""
        f = lambda x: 2.0 * x[0] + 3.0 * x[1]
        x = np.array([1.0, 1.0])
        grad = _finite_diff_gradient(f, x)
        assert abs(grad[0] - 2.0) < 1e-7, "∂f/∂x₀ = 2"
        assert abs(grad[1] - 3.0) < 1e-7, "∂f/∂x₁ = 3"

    def test_gradient_quadratic(self):
        r"""Gradient von f(x) = x₀² + x₁² ist [2x₀, 2x₁]."""
        f = lambda x: x[0]**2 + x[1]**2
        x = np.array([3.0, 4.0])
        grad = _finite_diff_gradient(f, x)
        assert abs(grad[0] - 6.0) < 1e-5
        assert abs(grad[1] - 8.0) < 1e-5

    def test_hessian_quadratic(self):
        r"""Hesse-Matrix von f(x) = x₀² + x₁² ist 2·I."""
        f = lambda x: x[0]**2 + x[1]**2
        x = np.array([1.0, 1.0])
        H = _finite_diff_hessian(f, x)
        assert H.shape == (2, 2)
        assert abs(H[0, 0] - 2.0) < 1e-4
        assert abs(H[1, 1] - 2.0) < 1e-4
        assert abs(H[0, 1]) < 1e-4  # Keine Kreuzterme


# ===========================================================================
# RIEMANNISCHE METRIK
# ===========================================================================

class TestRiemannianMetric:
    r"""Tests für die RiemannianMetric-Klasse."""

    def test_euclidean_metric_identity(self):
        r"""Euklidische Metrik (Standard) ist die Einheitsmatrix."""
        m = RiemannianMetric(dim=3)
        g = m.metric_at(np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_almost_equal(g, np.eye(3))

    def test_custom_metric(self):
        r"""Benutzerdefinierter Metrik-Tensor wird korrekt zurückgegeben."""
        g_func = lambda coords: np.diag([2.0, 3.0, 4.0])
        m = RiemannianMetric(dim=3, metric_tensor_func=g_func)
        g = m.metric_at(np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_almost_equal(g, np.diag([2.0, 3.0, 4.0]))

    def test_christoffel_euclidean_zero(self):
        r"""Christoffel-Symbole der euklidischen Metrik sind 0."""
        m = RiemannianMetric(dim=2)
        G = m.christoffel_symbols(np.array([1.0, 1.0]))
        # Für flache Metrik: alle Γ = 0
        assert G.shape == (2, 2, 2)
        np.testing.assert_array_almost_equal(G, np.zeros((2, 2, 2)), decimal=6)

    def test_ricci_tensor_euclidean_zero(self):
        r"""Ricci-Tensor der euklidischen Metrik ist 0."""
        m = RiemannianMetric(dim=2)
        Ric = m.ricci_tensor(np.array([0.5, 0.5]))
        assert Ric.shape == (2, 2)
        np.testing.assert_array_almost_equal(Ric, np.zeros((2, 2)), decimal=4)

    def test_ricci_scalar_euclidean_zero(self):
        r"""Ricci-Skalar der euklidischen Metrik ist 0."""
        m = RiemannianMetric(dim=2)
        R = m.ricci_scalar(np.array([0.5, 0.5]))
        assert abs(R) < 1e-4, f"Erwartete R≈0, erhalten R={R}"

    def test_volume_form_euclidean(self):
        r"""Volumenelement der euklidischen Metrik ist 1."""
        m = RiemannianMetric(dim=3)
        vol = m.volume_form(np.array([0.0, 0.0, 0.0]))
        assert abs(vol - 1.0) < 1e-10

    def test_volume_form_scaled_metric(self):
        r"""Volumenelement für g = 4·I ist √(4³) = 8."""
        g_func = lambda coords: 4.0 * np.eye(3)
        m = RiemannianMetric(dim=3, metric_tensor_func=g_func)
        vol = m.volume_form(np.array([0.0, 0.0, 0.0]))
        assert abs(vol - 8.0) < 1e-8, f"Erwartete 8, erhalten {vol}"

    def test_sectional_curvature_flat(self):
        r"""Schnittkrümmung der euklidischen Metrik ist 0."""
        m = RiemannianMetric(dim=3)
        u = np.array([1.0, 0.0, 0.0])
        v = np.array([0.0, 1.0, 0.0])
        K = m.sectional_curvature((u, v), np.array([0.0, 0.0, 0.0]))
        assert abs(K) < 1e-4, f"Schnittkrümmung der flachen Metrik sollte 0 sein, erhalten {K}"

    def test_sectional_curvature_dependent_vectors_raises(self):
        r"""Abhängige Vektoren lösen ValueError aus."""
        m = RiemannianMetric(dim=2)
        u = np.array([1.0, 0.0])
        v = np.array([2.0, 0.0])  # Linear abhängig von u
        with pytest.raises(ValueError):
            m.sectional_curvature((u, v), np.array([0.0, 0.0]))

    def test_riemann_tensor_euclidean_zero(self):
        r"""Riemann-Tensor der euklidischen Metrik ist 0."""
        m = RiemannianMetric(dim=2)
        R = m.riemann_curvature_tensor(np.array([0.5, 0.5]))
        assert R.shape == (2, 2, 2, 2)
        np.testing.assert_array_almost_equal(R, np.zeros((2, 2, 2, 2)), decimal=3)


# ===========================================================================
# DREI-SPHÄRE S³
# ===========================================================================

class TestThreeSphere:
    r"""Tests für die ThreeSphere-Klasse."""

    def test_default_radius(self):
        r"""Standard-Radius ist 1."""
        s3 = ThreeSphere()
        assert s3.radius == 1.0

    def test_custom_radius(self):
        r"""Benutzerdefinierter Radius wird gesetzt."""
        s3 = ThreeSphere(radius=2.5)
        assert abs(s3.radius - 2.5) < 1e-10

    def test_negative_radius_raises(self):
        r"""Negativer Radius löst ValueError aus."""
        with pytest.raises(ValueError):
            ThreeSphere(radius=-1.0)

    def test_zero_radius_raises(self):
        r"""Null-Radius löst ValueError aus."""
        with pytest.raises(ValueError):
            ThreeSphere(radius=0.0)

    def test_metric_tensor_shape(self):
        r"""Metrischer Tensor hat Form (3,3)."""
        s3 = ThreeSphere()
        g = s3.metric_tensor(np.array([0.5, 1.0, 2.0]))
        assert g.shape == (3, 3)

    def test_metric_tensor_diagonal(self):
        r"""Metrischer Tensor in Hopf-Koordinaten ist diagonal."""
        s3 = ThreeSphere()
        g = s3.metric_tensor(np.array([math.pi/4, 0.0, 0.0]))
        # Off-Diagonalelemente sollten 0 sein
        assert abs(g[0, 1]) < 1e-12
        assert abs(g[0, 2]) < 1e-12
        assert abs(g[1, 2]) < 1e-12

    def test_metric_positive_definite(self):
        r"""Metrischer Tensor ist positiv-definit für η ≠ 0, π/2."""
        s3 = ThreeSphere()
        g = s3.metric_tensor(np.array([math.pi/4, 1.0, 2.0]))
        eigenvalues = np.linalg.eigvalsh(g)
        assert all(ev > 0 for ev in eigenvalues), f"Nicht positiv-definit: {eigenvalues}"

    def test_metric_radius_scaling(self):
        r"""Metrischer Tensor skaliert mit R²."""
        s3 = ThreeSphere(radius=2.0)
        coords = np.array([math.pi/4, 0.0, 0.0])
        g = s3.metric_tensor(coords)
        # g[0,0] = R² = 4
        assert abs(g[0, 0] - 4.0) < 1e-10

    def test_fundamental_group(self):
        r"""Fundamentalgruppe von S³ ist trivial."""
        s3 = ThreeSphere()
        fg = s3.fundamental_group()
        assert "{1}" in fg or "trivial" in fg.lower()

    def test_homology_groups(self):
        r"""Homologiegruppen von S³ sind korrekt."""
        s3 = ThreeSphere()
        H = s3.homology_groups()
        assert H["H_0"] == "Z"
        assert H["H_1"] == "0"
        assert H["H_2"] == "0"
        assert H["H_3"] == "Z"

    def test_homotopy_groups(self):
        r"""Homotopiegruppen von S³ sind korrekt."""
        s3 = ThreeSphere()
        pi = s3.homotopy_groups()
        assert pi["pi_1"] == "0"
        assert pi["pi_2"] == "0"
        assert pi["pi_3"] == "Z"

    def test_is_simply_connected(self):
        r"""S³ ist einfach-zusammenhängend."""
        s3 = ThreeSphere()
        assert s3.is_simply_connected() is True


# ===========================================================================
# RICCI-FLUSS
# ===========================================================================

class TestRicciFlow:
    r"""Tests für die RicciFlow-Klasse."""

    def setup_method(self):
        r"""Erstellt eine euklidische Metrik für Tests."""
        self.metric = RiemannianMetric(dim=2)
        self.flow = RicciFlow(self.metric, dim=2)

    def test_evolution_equation_string(self):
        r"""Evolutionsgleichung enthält erwartete Symbole."""
        eq = self.flow.evolution_equation()
        assert "R_{ij}" in eq or "Ric" in eq or "-2" in eq

    def test_normalized_evolution_equation_string(self):
        r"""Normalisierte Gleichung enthält 'normalisiert' oder 'avg'."""
        eq = self.flow.normalized_evolution_equation()
        assert "normalisiert" in eq.lower() or "avg" in eq.lower() or "R_avg" in eq

    def test_short_time_existence_mentions_deturck(self):
        r"""Kurzzeit-Existenzsatz erwähnt DeTurck."""
        desc = self.flow.short_time_existence()
        assert "DeTurck" in desc or "1983" in desc

    def test_euler_step_flat_metric_unchanged(self):
        r"""Euler-Schritt einer flachen Metrik ändert die Metrik nicht."""
        coords = np.array([0.5, 0.5])
        g_new = self.flow.euler_step(dt=0.01, coords=coords)
        # Für flache Metrik: Ric = 0 → g_new = g_old
        np.testing.assert_array_almost_equal(g_new, np.eye(2), decimal=4)

    def test_euler_step_returns_correct_shape(self):
        r"""Euler-Schritt gibt Matrix der richtigen Form zurück."""
        coords = np.array([0.5, 0.5])
        g_new = self.flow.euler_step(dt=0.01, coords=coords)
        assert g_new.shape == (2, 2)

    def test_rk4_step_flat_metric_unchanged(self):
        r"""RK4-Schritt einer flachen Metrik ändert die Metrik nicht."""
        coords = np.array([0.5, 0.5])
        g_new = self.flow.rk4_step(dt=0.01, coords=coords)
        np.testing.assert_array_almost_equal(g_new, np.eye(2), decimal=3)

    def test_rk4_step_returns_correct_shape(self):
        r"""RK4-Schritt gibt Matrix der richtigen Form zurück."""
        coords = np.array([0.5, 0.5])
        g_new = self.flow.rk4_step(dt=0.01, coords=coords)
        assert g_new.shape == (2, 2)


# ===========================================================================
# RICCI-FLUSS-SINGULARITÄTEN
# ===========================================================================

class TestRicciFlowSingularity:
    r"""Tests für die RicciFlowSingularity-Klasse."""

    def test_valid_singularity_types(self):
        r"""Alle drei Singularitätstypen können instanziiert werden."""
        for stype in ["neck-pinch", "cap", "cigar"]:
            s = RicciFlowSingularity(stype)
            assert s.singularity_type == stype

    def test_invalid_type_raises(self):
        r"""Ungültiger Typ löst ValueError aus."""
        with pytest.raises(ValueError):
            RicciFlowSingularity("unknown-type")

    def test_blow_up_rate_neck_pinch(self):
        r"""Aufblasrate für Nack-Pinch ist ~ 1/(T-t)."""
        s = RicciFlowSingularity("neck-pinch")
        rate = s.blow_up_rate(t=0.5, T=1.0)
        assert abs(rate - 2.0) < 1e-10  # 1/(1.0-0.5) = 2.0

    def test_blow_up_rate_increases_near_t(self):
        r"""Aufblasrate wächst wenn t → T."""
        s = RicciFlowSingularity("neck-pinch")
        rate1 = s.blow_up_rate(t=0.0, T=1.0)
        rate2 = s.blow_up_rate(t=0.9, T=1.0)
        assert rate2 > rate1, "Aufblasrate sollte nahe T größer sein"

    def test_blow_up_rate_past_t_raises(self):
        r"""t ≥ T löst ValueError aus."""
        s = RicciFlowSingularity("neck-pinch")
        with pytest.raises(ValueError):
            s.blow_up_rate(t=1.0, T=1.0)

    def test_canonical_neighborhood_neck_pinch(self):
        r"""Kanonische Nachbarschaft des Nack-Pinch erwähnt ε-Hals."""
        s = RicciFlowSingularity("neck-pinch")
        desc = s.canonical_neighborhood()
        assert "Hals" in desc or "neck" in desc.lower() or "ε-neck" in desc.lower()

    def test_canonical_neighborhood_cigar_excludes(self):
        r"""Zigarren-Singularität wird durch κ-Nichtkollaps ausgeschlossen."""
        s = RicciFlowSingularity("cigar")
        desc = s.canonical_neighborhood()
        assert "ausgeschlossen" in desc or "κ" in desc or "Nichtkollaps" in desc

    def test_kappa_noncollapsed_neck_pinch(self):
        r"""Nack-Pinch erfüllt κ-Nichtkollaps."""
        s = RicciFlowSingularity("neck-pinch")
        result = s.kappa_noncollapsed_check(r=1.0, kappa=0.1)
        assert result is True

    def test_kappa_noncollapsed_cigar_false(self):
        r"""Zigarren-Singularität verletzt κ-Nichtkollaps."""
        s = RicciFlowSingularity("cigar")
        result = s.kappa_noncollapsed_check(r=1.0, kappa=0.1)
        assert result is False

    def test_kappa_noncollapsed_invalid_r_raises(self):
        r"""r ≤ 0 löst ValueError aus."""
        s = RicciFlowSingularity("neck-pinch")
        with pytest.raises(ValueError):
            s.kappa_noncollapsed_check(r=0.0, kappa=0.1)


# ===========================================================================
# PERELMAN-ENTROPIE
# ===========================================================================

class TestPerelmanEntropy:
    r"""Tests für die PerelmanEntropy-Klasse."""

    def setup_method(self):
        r"""Initialisiert Entropie-Objekt mit euklidischer Metrik."""
        self.metric = RiemannianMetric(dim=1)
        self.f_func = lambda coords: float(np.sum(coords**2))  # f(x) = |x|²
        self.entropy = PerelmanEntropy(self.metric, self.f_func, tau=1.0)

    def test_invalid_tau_raises(self):
        r"""τ ≤ 0 löst ValueError aus."""
        metric = RiemannianMetric(dim=1)
        with pytest.raises(ValueError):
            PerelmanEntropy(metric, lambda x: 0.0, tau=0.0)

    def test_tau_stored_correctly(self):
        r"""τ wird korrekt gespeichert."""
        assert abs(self.entropy.tau - 1.0) < 1e-10

    def test_w_functional_returns_float(self):
        r"""W-Funktional gibt einen float zurück."""
        result = self.entropy.w_functional(self.f_func, num_points=10)
        assert isinstance(result, float)

    def test_f_functional_returns_float(self):
        r"""F-Funktional gibt einen float zurück."""
        result = self.entropy.f_functional(self.f_func, num_points=10)
        assert isinstance(result, float)

    def test_mu_functional_returns_float(self):
        r"""μ-Funktional gibt einen float zurück."""
        result = self.entropy.mu_functional(tau=1.0)
        assert isinstance(result, float)

    def test_nu_functional_returns_float(self):
        r"""ν-Funktional gibt einen float zurück."""
        result = self.entropy.nu_functional()
        assert isinstance(result, float)

    def test_nu_leq_mu(self):
        r"""ν(g) ≤ μ(g,τ) für τ = 1.0 (Infimum-Eigenschaft)."""
        nu = self.entropy.nu_functional()
        mu = self.entropy.mu_functional(tau=1.0)
        assert nu <= mu + 1e-8, f"ν={nu} sollte ≤ μ={mu} sein"

    def test_monotonicity_check_valid(self):
        r"""Monotonie-Check für t1 < t2 gibt bool zurück."""
        result = self.entropy.monotonicity_check(t1=0.0, t2=1.0)
        assert isinstance(result, bool)

    def test_monotonicity_check_invalid_raises(self):
        r"""t1 ≥ t2 löst ValueError aus."""
        with pytest.raises(ValueError):
            self.entropy.monotonicity_check(t1=1.0, t2=0.5)


# ===========================================================================
# RICCI-FLUSS MIT CHIRURGIE
# ===========================================================================

class TestRicciFlowWithSurgery:
    r"""Tests für die RicciFlowWithSurgery-Klasse."""

    def setup_method(self):
        r"""Initialisiert Ricci-Fluss-mit-Chirurgie-Objekt."""
        self.metric = RiemannianMetric(dim=3)
        self.flow = RicciFlowWithSurgery(self.metric)

    def test_initial_surgery_count_zero(self):
        r"""Anfangszähler für Chirurgien ist 0."""
        assert self.flow.surgery_count == 0

    def test_no_singularity_early_time(self):
        r"""Keine Singularität bei frühem Zeitpunkt t=0.1."""
        result = self.flow.detect_singularity(t=0.1)
        assert result is None

    def test_singularity_detected_near_blowup(self):
        r"""Singularität wird kurz vor Blow-up-Zeit erkannt."""
        result = self.flow.detect_singularity(t=0.99)
        assert result is not None
        assert "pinch" in result.lower() or "neck" in result.lower()

    def test_singularity_detected_at_blowup(self):
        r"""Singularität wird bei Blow-up-Zeit erkannt."""
        result = self.flow.detect_singularity(t=1.0)
        assert result is not None

    def test_perform_surgery_increments_counter(self):
        r"""Chirurgie erhöht den Zähler."""
        self.flow.perform_surgery("neck-pinch")
        assert self.flow.surgery_count == 1

    def test_perform_surgery_returns_dict(self):
        r"""Chirurgie gibt ein Dictionary zurück."""
        result = self.flow.perform_surgery("neck-pinch")
        assert isinstance(result, dict)
        assert "surgery_number" in result
        assert "description" in result

    def test_perform_multiple_surgeries(self):
        r"""Mehrere Chirurgien werden korrekt protokolliert."""
        self.flow.perform_surgery("neck-pinch")
        self.flow.perform_surgery("cap")
        assert self.flow.surgery_count == 2
        assert len(self.flow.surgery_log) == 2

    def test_canonical_neighborhood_theorem_string(self):
        r"""Kanonisches Nachbarschaftstheorem erwähnt Perelman."""
        desc = self.flow.canonical_neighborhood_theorem()
        assert "Perelman" in desc or "12.1" in desc

    def test_extinction_time_positive(self):
        r"""Auslöschungszeit ist positiv."""
        T = self.flow.extinction_time()
        assert T is not None
        assert T > 0


# ===========================================================================
# POINCARÉ-VERMUTUNG
# ===========================================================================

class TestPoincaréConjecture:
    r"""Tests für die PoincaréConjecture-Klasse."""

    def setup_method(self):
        r"""Initialisiert das Poincaré-Objekt."""
        self.pc = PoincaréConjecture()

    def test_theorem_statement_contains_poincare(self):
        r"""Theorem-Aussage enthält 'Poincaré'."""
        stmt = self.pc.theorem_statement()
        assert "Poincaré" in stmt or "Perelman" in stmt

    def test_theorem_statement_mentions_s3(self):
        r"""Theorem-Aussage erwähnt S³."""
        stmt = self.pc.theorem_statement()
        assert "S³" in stmt or "S^3" in stmt or "3-Sphäre" in stmt

    def test_geometrization_mentions_thurston(self):
        r"""Geometrisierungsvermutung erwähnt Thurston."""
        stmt = self.pc.geometrization_statement()
        assert "Thurston" in stmt

    def test_geometrization_mentions_8_geometries(self):
        r"""Geometrisierungsvermutung erwähnt 8 Geometrien."""
        stmt = self.pc.geometrization_statement()
        assert "8" in stmt or "acht" in stmt.lower()

    def test_eight_geometries_count(self):
        r"""Liste der 8 Thurston-Geometrien hat genau 8 Einträge."""
        geos = self.pc.eight_geometries()
        assert len(geos) == 8

    def test_eight_geometries_contain_s3(self):
        r"""S³ ist in den acht Geometrien enthalten."""
        geos = self.pc.eight_geometries()
        combined = " ".join(geos)
        assert "S³" in combined or "S^3" in combined

    def test_eight_geometries_contain_h3(self):
        r"""H³ ist in den acht Geometrien enthalten."""
        geos = self.pc.eight_geometries()
        combined = " ".join(geos)
        assert "H³" in combined or "H^3" in combined

    def test_proof_outline_has_steps(self):
        r"""Beweis-Outline enthält mindestens 4 Schritte."""
        steps = self.pc.proof_outline()
        assert len(steps) >= 4

    def test_proof_outline_mentions_entropy(self):
        r"""Beweis-Outline erwähnt Entropie."""
        steps = self.pc.proof_outline()
        combined = " ".join(steps)
        assert "Entropie" in combined or "entropy" in combined.lower() or "W-" in combined

    def test_perelman_papers_count(self):
        r"""Perelman hat genau 3 arXiv-Papers."""
        papers = self.pc.perelman_papers()
        assert len(papers) == 3

    def test_perelman_papers_arxiv_ids(self):
        r"""Alle Papers haben arXiv-IDs."""
        papers = self.pc.perelman_papers()
        arxiv_ids = [p["arxiv"] for p in papers]
        assert "math/0211159" in arxiv_ids
        assert "math/0303109" in arxiv_ids
        assert "math/0307245" in arxiv_ids

    def test_perelman_papers_dates(self):
        r"""Papers haben das richtige Jahr 2002 bzw. 2003."""
        papers = self.pc.perelman_papers()
        years = [p["date"][:4] for p in papers]
        assert "2002" in years
        assert "2003" in years

    def test_fields_medal_declined_mentions_2006(self):
        r"""Fields-Medaille-Notiz erwähnt 2006."""
        note = self.pc.fields_medal_declined()
        assert "2006" in note

    def test_fields_medal_declined_mentions_perelman(self):
        r"""Fields-Medaille-Notiz erwähnt Perelman."""
        note = self.pc.fields_medal_declined()
        assert "Perelman" in note


# ===========================================================================
# GEOMETRISIERUNGSVERMUTUNG
# ===========================================================================

class TestGeometrizationConjecture:
    r"""Tests für die GeometrizationConjecture-Klasse."""

    def setup_method(self):
        r"""Initialisiert das Geometrisierungs-Objekt."""
        self.gc = GeometrizationConjecture()

    def test_classify_s3(self):
        r"""S³ wird als sphärische Geometrie klassifiziert."""
        result = self.gc.thurston_classification("S3")
        assert "S³" in result or "sphärisch" in result.lower()

    def test_classify_hyperbolic(self):
        r"""Hyperbolische Mannigfaltigkeit wird korrekt klassifiziert."""
        result = self.gc.thurston_classification("hyperbolic")
        assert "H³" in result or "hyperbol" in result.lower()

    def test_classify_flat_torus(self):
        r"""Flacher Torus wird als E³-Geometrie klassifiziert."""
        result = self.gc.thurston_classification("flat-torus")
        assert "E³" in result or "flach" in result.lower() or "euklidisch" in result.lower()

    def test_classify_unknown_type(self):
        r"""Unbekannter Typ gibt Fallback-String zurück."""
        result = self.gc.thurston_classification("some-unknown-manifold")
        assert "JSJ" in result or "Unbekannt" in result or "unbekannt" in result.lower()

    def test_prime_decomposition_s3(self):
        r"""Primzerlegung von S³ ist [S³]."""
        parts = self.gc.prime_decomposition("S3")
        assert len(parts) == 1
        assert "S³" in parts[0]

    def test_jso_decomposition_string(self):
        r"""JSJ-Zerlegung gibt beschreibenden String zurück."""
        desc = self.gc.jaco_shalen_johannsen_decomposition()
        assert "JSJ" in desc or "Jaco" in desc

    def test_jso_decomposition_mentions_tori(self):
        r"""JSJ-Zerlegung erwähnt Tori als Zerlegungsflächen."""
        desc = self.gc.jaco_shalen_johannsen_decomposition()
        assert "Tori" in desc or "tori" in desc.lower() or "T²" in desc

    def test_hyperbolic_pieces_mentions_mostow(self):
        r"""Hyperbolische Stücke erwähnen Mostow-Rigidität."""
        desc = self.gc.hyperbolic_pieces()
        assert "Mostow" in desc


# ===========================================================================
# UNIFORMISIERUNGSSATZ
# ===========================================================================

class TestUniformizationTheorem:
    r"""Tests für die UniformizationTheorem-Klasse."""

    def test_genus_0_classify(self):
        r"""Genus 0 → S²."""
        ut = UniformizationTheorem(genus=0)
        result = ut.classify_surface()
        assert "S²" in result

    def test_genus_1_classify(self):
        r"""Genus 1 → Torus T²."""
        ut = UniformizationTheorem(genus=1)
        result = ut.classify_surface()
        assert "T²" in result or "Torus" in result.lower()

    def test_genus_2_classify(self):
        r"""Genus 2 → Σ₂ (hyperbolisch)."""
        ut = UniformizationTheorem(genus=2)
        result = ut.classify_surface()
        assert "Σ_2" in result or "hyperbol" in result.lower()

    def test_negative_genus_raises(self):
        r"""Negativer Genus löst ValueError aus."""
        with pytest.raises(ValueError):
            UniformizationTheorem(genus=-1)

    def test_euler_characteristic_genus_0(self):
        r"""χ(S²) = 2."""
        ut = UniformizationTheorem(genus=0)
        assert ut.euler_char == 2

    def test_euler_characteristic_genus_1(self):
        r"""χ(T²) = 0."""
        ut = UniformizationTheorem(genus=1)
        assert ut.euler_char == 0

    def test_euler_characteristic_genus_2(self):
        r"""χ(Σ₂) = -2."""
        ut = UniformizationTheorem(genus=2)
        assert ut.euler_char == -2

    def test_universal_cover_s2(self):
        r"""Universelle Überlagerung von S² ist S²."""
        ut = UniformizationTheorem(genus=0)
        uc = ut.universal_cover()
        assert "S²" in uc

    def test_universal_cover_torus(self):
        r"""Universelle Überlagerung von T² ist ℝ²."""
        ut = UniformizationTheorem(genus=1)
        uc = ut.universal_cover()
        assert "ℝ²" in uc or "R²" in uc or "Euklidisch" in uc or "euklidisch" in uc.lower()

    def test_universal_cover_genus_2(self):
        r"""Universelle Überlagerung von Σ₂ ist H²."""
        ut = UniformizationTheorem(genus=2)
        uc = ut.universal_cover()
        assert "H²" in uc or "hyperbol" in uc.lower()

    def test_gauss_bonnet_s2(self):
        r"""Gauss-Bonnet für S²: ∫K dA = 4π."""
        ut = UniformizationTheorem(genus=0)
        val = ut.gauss_bonnet(K_avg=1.0)
        assert abs(val - 4.0 * math.pi) < 1e-10

    def test_gauss_bonnet_torus(self):
        r"""Gauss-Bonnet für T²: ∫K dA = 0."""
        ut = UniformizationTheorem(genus=1)
        val = ut.gauss_bonnet(K_avg=0.0)
        assert abs(val) < 1e-10

    def test_gauss_bonnet_genus_2(self):
        r"""Gauss-Bonnet für Σ₂: ∫K dA = -4π."""
        ut = UniformizationTheorem(genus=2)
        val = ut.gauss_bonnet(K_avg=-1.0)
        assert abs(val - (-4.0 * math.pi)) < 1e-10

    def test_ricci_flow_2d_initial(self):
        r"""Ricci-Fluss 2D bei t=0 beschreibt Anfangsmetrik."""
        ut = UniformizationTheorem(genus=0)
        desc = ut.ricci_flow_2d(t=0.0)
        assert "t=0" in desc or "Anfang" in desc.lower()

    def test_ricci_flow_2d_convergence_s2(self):
        r"""Ricci-Fluss auf S² konvergiert zur runden Metrik."""
        ut = UniformizationTheorem(genus=0)
        desc = ut.ricci_flow_2d(t=10.0)
        assert "konvergiert" in desc.lower() or "runde" in desc.lower() or "Hamilton" in desc

    def test_ricci_flow_2d_convergence_hyperbolic(self):
        r"""Ricci-Fluss auf Σ₂ konvergiert zur hyperbolischen Metrik."""
        ut = UniformizationTheorem(genus=2)
        desc = ut.ricci_flow_2d(t=10.0)
        assert "hyperbol" in desc.lower() or "K=-1" in desc

    def test_ricci_flow_2d_negative_t_raises(self):
        r"""Negativer Zeitpunkt löst ValueError aus."""
        ut = UniformizationTheorem(genus=1)
        with pytest.raises(ValueError):
            ut.ricci_flow_2d(t=-1.0)
