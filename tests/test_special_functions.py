"""
Tests für das Modul spezieller Funktionen (src/special_functions.py).

Testet Bessel-, Legendre-, hypergeometrische, Airy-Funktionen,
orthogonale Polynome und klassische Sonderfälle.

Autor: Kurt Ingwer
Letzte Änderung: 2026-03-10
"""

import pytest
import numpy as np
from scipy import special as sc

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from special_functions import (
    BesselFunctions, LegendrePolynomials, HypergeometricFunctions,
    AiryFunctions, OrthogonalPolynomials, gamma_function_properties,
    beta_function, digamma_function, riemann_zeta_special_values,
    elliptic_integrals, error_function_properties,
)


# ---------------------------------------------------------------------------
# Bessel-Funktionen Tests
# ---------------------------------------------------------------------------

class TestBesselFunctions:
    """Tests für Bessel-Funktionen."""

    def setup_method(self):
        self.bf = BesselFunctions()

    def test_J0_at_zero(self):
        """J₀(0) = 1."""
        assert abs(self.bf.J(0, 0) - 1.0) < 1e-10

    def test_J1_at_zero(self):
        """J₁(0) = 0."""
        assert abs(self.bf.J(1, 0)) < 1e-10

    def test_J0_known_values(self):
        """J₀(2.4048) ≈ 0 (erste Nullstelle)."""
        # Erste Nullstelle von J₀ ≈ 2.4048
        first_zero = self.bf.zeros(0, 1)[0]
        assert abs(self.bf.J(0, first_zero)) < 1e-6

    def test_J_recurrence(self):
        """Rekurrenzrelation J_{ν-1} + J_{ν+1} = (2ν/x)·J_ν."""
        assert self.bf.recurrence(2.0, 3.0)
        assert self.bf.recurrence(1.0, 5.0)

    def test_Y_singular_at_zero(self):
        """Y₀(0) = -∞."""
        assert self.bf.Y(0, 0) == float('-inf')

    def test_Y_positive_for_large_x(self):
        """Y₀(x) existiert für x > 0."""
        val = self.bf.Y(0, 5.0)
        assert np.isfinite(val)

    def test_I_at_zero(self):
        """I₀(0) = 1."""
        assert abs(self.bf.I(0, 0) - 1.0) < 1e-10

    def test_K_decreasing(self):
        """K₀(x) ist abnehmend für x > 0."""
        k1 = self.bf.K(0, 1.0)
        k2 = self.bf.K(0, 2.0)
        assert k1 > k2 > 0

    def test_zeros_first_n(self):
        """Erste 3 Nullstellen von J₀ sind positiv und wachsend."""
        zeros = self.bf.zeros(0, 3)
        assert len(zeros) == 3
        assert all(z > 0 for z in zeros)
        assert zeros[0] < zeros[1] < zeros[2]

    def test_bessel_ode_string(self):
        """ODE-String enthält 'y'."""
        ode = self.bf.bessel_ode(1.5)
        assert 'y' in ode
        assert '1.5' in ode

    def test_J_negative_order(self):
        """J_{-n}(x) = (-1)^n J_n(x)."""
        x = 3.0
        # J_{-1}(x) = -J_1(x)
        assert abs(self.bf.J(-1, x) + self.bf.J(1, x)) < 1e-10


# ---------------------------------------------------------------------------
# Legendre-Polynome Tests
# ---------------------------------------------------------------------------

class TestLegendrePolynomials:
    """Tests für Legendre-Polynome."""

    def setup_method(self):
        self.lp = LegendrePolynomials()

    def test_P0(self):
        """P₀(x) = 1."""
        assert abs(self.lp.P(0, 0.5) - 1.0) < 1e-10
        assert abs(self.lp.P(0, -0.3) - 1.0) < 1e-10

    def test_P1(self):
        """P₁(x) = x."""
        x = 0.7
        assert abs(self.lp.P(1, x) - x) < 1e-10

    def test_P2(self):
        """P₂(x) = (3x²-1)/2."""
        x = 0.5
        expected = (3 * x ** 2 - 1) / 2
        assert abs(self.lp.P(2, x) - expected) < 1e-10

    def test_P3(self):
        """P₃(x) = (5x³-3x)/2."""
        x = 0.3
        expected = (5 * x ** 3 - 3 * x) / 2
        assert abs(self.lp.P(3, x) - expected) < 1e-8

    def test_orthogonality_different_n_m(self):
        """∫₋₁¹ P_n·P_m dx = 0 für n ≠ m."""
        result = self.lp.orthogonality_check(2, 3)
        assert abs(result) < 1e-8

    def test_orthogonality_same_n(self):
        """∫₋₁¹ P_n² dx = 2/(2n+1)."""
        n = 3
        result = self.lp.orthogonality_check(n, n)
        expected = 2.0 / (2 * n + 1)
        assert abs(result - expected) < 1e-8

    def test_rodrigues_matches_direct(self):
        """Rodrigues-Formel stimmt mit direkter Berechnung überein."""
        for n in range(5):
            for x in [-0.5, 0.0, 0.5]:
                assert abs(self.lp.rodrigues(n, x) - self.lp.P(n, x)) < 1e-10

    def test_generating_function(self):
        """Erzeugende Funktion stimmt mit Partialsumme überein."""
        x, t = 0.5, 0.3
        G = self.lp.generating_function(x, t)
        # Erste Terme der Reihe
        approx = sum(self.lp.P(n, x) * t ** n for n in range(20))
        assert abs(G - approx) < 1e-6

    def test_generating_function_invalid_t(self):
        """|t| ≥ 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            self.lp.generating_function(0.5, 1.5)

    def test_associated_legendre(self):
        """P_1^1(x) = -(1-x²)^{1/2}."""
        x = 0.5
        # P_1^1(x) = -√(1-x²)
        expected = -np.sqrt(1 - x ** 2)
        result = self.lp.P_assoc(1, 1, x)
        assert abs(result - expected) < 1e-8


# ---------------------------------------------------------------------------
# Hypergeometrische Funktionen Tests
# ---------------------------------------------------------------------------

class TestHypergeometricFunctions:
    """Tests für hypergeometrische Funktionen."""

    def setup_method(self):
        self.hf = HypergeometricFunctions()

    def test_F21_at_zero(self):
        """₂F₁(a,b;c;0) = 1."""
        assert abs(self.hf.F21(1, 2, 3, 0) - 1.0) < 1e-10

    def test_F21_geometric_series(self):
        """₂F₁(1,1;1;z) = 1/(1-z)."""
        z = 0.3
        result = self.hf.F21(1, 1, 1, z)
        expected = 1.0 / (1 - z)
        assert abs(result - expected) < 1e-8

    def test_F10_binomial(self):
        """₁F₀(a;;z) = (1-z)^{-a}."""
        a, z = 2.0, 0.4
        result = self.hf.F10(a, z)
        expected = (1 - z) ** (-a)
        assert abs(result - expected) < 1e-10

    def test_F10_invalid_z(self):
        """|z| ≥ 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            self.hf.F10(2.0, 1.5)

    def test_F11_at_zero(self):
        """₁F₁(a;c;0) = 1."""
        assert abs(self.hf.F11(1, 2, 0) - 1.0) < 1e-10

    def test_F11_exponential(self):
        """₁F₁(a;a;z) = e^z."""
        a, z = 2.0, 1.5
        result = self.hf.F11(a, a, z)
        expected = np.exp(z)
        assert abs(result - expected) < 1e-8

    def test_gauss_ode_string(self):
        """ODE-String der hypergeometrischen DGL."""
        ode = self.hf.Gauss_hypergeometric_ode()
        assert 'w' in ode
        assert "z(1-z)" in ode

    def test_special_cases_dict(self):
        """special_cases gibt Dict mit mindestens 5 Einträgen zurück."""
        cases = self.hf.special_cases()
        assert isinstance(cases, dict)
        assert len(cases) >= 5

    def test_F21_legendre_connection(self):
        """₂F₁(-n,n+1;1;(1-x)/2) = P_n(x)."""
        lp = LegendrePolynomials()
        n = 2
        x = 0.5
        z = (1 - x) / 2
        result = self.hf.F21(-n, n + 1, 1, z)
        expected = lp.P(n, x)
        assert abs(result - expected) < 1e-6


# ---------------------------------------------------------------------------
# Airy-Funktionen Tests
# ---------------------------------------------------------------------------

class TestAiryFunctions:
    """Tests für Airy-Funktionen."""

    def setup_method(self):
        self.af = AiryFunctions()

    def test_Ai_at_zero(self):
        """Ai(0) = Γ(2/3)/(3^{2/3}·Γ(1)) ≈ 0.3551."""
        result = self.af.Ai(0)
        expected = float(sc.airy(0)[0])
        assert abs(result - expected) < 1e-10

    def test_Bi_at_zero(self):
        """Bi(0) ist endlich."""
        result = self.af.Bi(0)
        assert np.isfinite(result)

    def test_Ai_decreasing_for_positive_x(self):
        """Ai(x) nimmt für positive x ab."""
        vals = [self.af.Ai(x) for x in [0, 1, 2, 3]]
        assert vals[0] > vals[1] > vals[2] > vals[3]

    def test_Bi_increasing_for_positive_x(self):
        """Bi(x) nimmt für positive x zu."""
        vals = [self.af.Bi(x) for x in [0, 1, 2, 3]]
        assert vals[0] < vals[1] < vals[2] < vals[3]

    def test_Ai_zeros(self):
        """Erste Nullstellen von Ai liegen auf negativer reeller Achse."""
        zeros = self.af.zeros_Ai(3)
        assert len(zeros) == 3
        assert all(z < 0 for z in zeros)

    def test_Ai_zeros_ordered(self):
        """Nullstellen sind nach Betrag wachsend geordnet (alle negativ)."""
        zeros = self.af.zeros_Ai(5)
        # sc.ai_zeros gibt negative Werte vom betragsmäßig kleinsten zum größten
        assert all(abs(zeros[i]) < abs(zeros[i + 1]) for i in range(len(zeros) - 1))

    def test_airy_ode(self):
        """Ai erfüllt die Airy-DGL: y'' - x·y = 0."""
        # Numerische Verifikation an x=1
        x = 1.0
        h = 1e-5
        Ai_minus = self.af.Ai(x - h)
        Ai_val = self.af.Ai(x)
        Ai_plus = self.af.Ai(x + h)
        # y'' ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
        second_deriv = (Ai_plus - 2 * Ai_val + Ai_minus) / h ** 2
        # Airy-DGL: y'' = x·y → residual ≈ 0
        residual = abs(second_deriv - x * Ai_val)
        assert residual < 1e-4

    def test_Ai_prime_at_zero(self):
        """Ai'(0) ≈ -0.2588."""
        result = self.af.Ai_prime(0)
        expected = float(sc.airy(0)[1])
        assert abs(result - expected) < 1e-10

    def test_asymptotic_Ai_accuracy(self):
        """Asymptotik stimmt für großes x näherungsweise."""
        x = 10.0
        exact = self.af.Ai(x)
        asymp = self.af.asymptotic_Ai(x)
        # Relativer Fehler sollte < 1%
        rel_error = abs(exact - asymp) / abs(exact)
        assert rel_error < 0.01

    def test_asymptotic_invalid_x(self):
        """Asymptotik für x ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            self.af.asymptotic_Ai(-1.0)


# ---------------------------------------------------------------------------
# Orthogonale Polynome Tests
# ---------------------------------------------------------------------------

class TestOrthogonalPolynomials:
    """Tests für orthogonale Polynome."""

    def setup_method(self):
        self.op = OrthogonalPolynomials()

    def test_hermite_H0(self):
        """H₀(x) = 1."""
        assert abs(self.op.hermite(0, 2.5) - 1.0) < 1e-10

    def test_hermite_H1(self):
        """H₁(x) = 2x."""
        x = 1.5
        assert abs(self.op.hermite(1, x) - 2 * x) < 1e-8

    def test_hermite_H2(self):
        """H₂(x) = 4x² - 2."""
        x = 1.0
        expected = 4 * x ** 2 - 2
        assert abs(self.op.hermite(2, x) - expected) < 1e-8

    def test_hermite_orthogonality(self):
        """H_n und H_m sind orthogonal für n ≠ m."""
        result = self.op.orthogonality_verify('hermite', 2, 3)
        assert abs(result) < 1e-5

    def test_laguerre_L0(self):
        """L₀(x) = 1."""
        assert abs(self.op.laguerre(0, 2.0) - 1.0) < 1e-10

    def test_laguerre_L1(self):
        """L₁(x) = 1 - x."""
        x = 2.5
        assert abs(self.op.laguerre(1, x) - (1 - x)) < 1e-8

    def test_laguerre_orthogonality(self):
        """L_n und L_m sind orthogonal für n ≠ m."""
        result = self.op.orthogonality_verify('laguerre', 1, 2)
        assert abs(result) < 1e-5

    def test_chebyshev_T_at_cos(self):
        """T_n(cos θ) = cos(nθ)."""
        theta = 0.7
        x = np.cos(theta)
        n = 4
        result = self.op.chebyshev_T(n, x)
        expected = np.cos(n * theta)
        assert abs(result - expected) < 1e-8

    def test_chebyshev_T0(self):
        """T₀(x) = 1."""
        assert abs(self.op.chebyshev_T(0, 0.5) - 1.0) < 1e-10

    def test_chebyshev_T1(self):
        """T₁(x) = x."""
        x = 0.7
        assert abs(self.op.chebyshev_T(1, x) - x) < 1e-10

    def test_chebyshev_U0(self):
        """U₀(x) = 1."""
        assert abs(self.op.chebyshev_U(0, 0.5) - 1.0) < 1e-10

    def test_chebyshev_U1(self):
        """U₁(x) = 2x."""
        x = 0.4
        assert abs(self.op.chebyshev_U(1, x) - 2 * x) < 1e-10

    def test_gegenbauer_reduces_to_legendre(self):
        """C_n^{1/2}(x) = P_n(x) (Gegenbauer → Legendre für α=1/2)."""
        lp = LegendrePolynomials()
        n = 3
        x = 0.5
        gegenbauer_val = self.op.gegenbauer(n, 0.5, x)
        legendre_val = lp.P(n, x)
        assert abs(gegenbauer_val - legendre_val) < 1e-8


# ---------------------------------------------------------------------------
# Gamma-Funktion Tests
# ---------------------------------------------------------------------------

class TestGammaFunction:
    """Tests für die Gamma-Funktion."""

    def test_gamma_1(self):
        """Γ(1) = 1."""
        result = gamma_function_properties()
        assert abs(result['Gamma_1'] - 1.0) < 1e-10

    def test_gamma_factorial(self):
        """Γ(n+1) = n!"""
        result = gamma_function_properties()
        assert abs(result['Gamma_3'] - 2.0) < 1e-10  # Γ(3) = 2!
        assert abs(result['Gamma_4'] - 6.0) < 1e-10  # Γ(4) = 3!

    def test_gamma_half(self):
        """Γ(1/2) = √π."""
        result = gamma_function_properties()
        assert abs(result['Gamma_half'] - np.sqrt(np.pi)) < 1e-10

    def test_functional_equation(self):
        """Γ(z+1) = z·Γ(z)."""
        result = gamma_function_properties()
        assert result['functional_equation_check']['matches']

    def test_reflection_formula(self):
        """Reflexionsformel: Γ(z)·Γ(1-z) = π/sin(πz)."""
        result = gamma_function_properties()
        assert result['reflection_formula']['matches']

    def test_duplication_formula(self):
        """Duplikationsformel für Gamma."""
        result = gamma_function_properties()
        assert result['duplication_formula']['matches']


# ---------------------------------------------------------------------------
# Beta-Funktion Tests
# ---------------------------------------------------------------------------

class TestBetaFunction:
    """Tests für die Beta-Funktion."""

    def test_beta_symmetry(self):
        """B(a,b) = B(b,a)."""
        assert abs(beta_function(2, 3) - beta_function(3, 2)) < 1e-10

    def test_beta_known_value(self):
        """B(1,1) = 1."""
        assert abs(beta_function(1, 1) - 1.0) < 1e-10

    def test_beta_half_half(self):
        """B(1/2, 1/2) = π."""
        result = beta_function(0.5, 0.5)
        assert abs(result - np.pi) < 1e-8

    def test_beta_gamma_relation(self):
        """B(a,b) = Γ(a)Γ(b)/Γ(a+b)."""
        a, b = 2.0, 3.0
        result = beta_function(a, b)
        expected = float(sc.gamma(a) * sc.gamma(b) / sc.gamma(a + b))
        assert abs(result - expected) < 1e-10

    def test_beta_invalid_args(self):
        """a ≤ 0 oder b ≤ 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            beta_function(-1, 2)
        with pytest.raises(ValueError):
            beta_function(1, 0)


# ---------------------------------------------------------------------------
# Digamma-Funktion Tests
# ---------------------------------------------------------------------------

class TestDigammaFunction:
    """Tests für die Digamma-Funktion."""

    def test_digamma_1(self):
        """ψ(1) = -γ (Euler-Mascheroni-Konstante)."""
        gamma_euler = 0.5772156649015328
        result = digamma_function(1.0)
        assert abs(result + gamma_euler) < 1e-8

    def test_digamma_recurrence(self):
        """ψ(x+1) = ψ(x) + 1/x."""
        x = 2.5
        psi_x = digamma_function(x)
        psi_x1 = digamma_function(x + 1)
        assert abs(psi_x1 - (psi_x + 1.0 / x)) < 1e-8

    def test_digamma_integer_2(self):
        """ψ(2) = 1 - γ."""
        gamma_euler = 0.5772156649015328
        result = digamma_function(2.0)
        assert abs(result - (1.0 - gamma_euler)) < 1e-8


# ---------------------------------------------------------------------------
# Riemann-Zeta Tests
# ---------------------------------------------------------------------------

class TestRiemannZeta:
    """Tests für die Riemann-Zeta-Funktion (spezielle Werte)."""

    def test_zeta_2(self):
        """ζ(2) = π²/6."""
        result = riemann_zeta_special_values()
        assert result['verifications']['zeta_2_vs_pi_squared_over_6']['match']

    def test_zeta_4(self):
        """ζ(4) = π⁴/90."""
        result = riemann_zeta_special_values()
        assert result['verifications']['zeta_4_vs_pi4_over_90']['match']

    def test_zeta_bernoulli(self):
        """Bernoulli-Zusammenhang im Dict vorhanden."""
        result = riemann_zeta_special_values()
        assert 'B_2' in result['bernoulli_connection']
        assert abs(result['bernoulli_connection']['B_2'] - 1.0 / 6.0) < 1e-10


# ---------------------------------------------------------------------------
# Elliptische Integrale Tests
# ---------------------------------------------------------------------------

class TestEllipticIntegrals:
    """Tests für elliptische Integrale."""

    def test_K_at_zero(self):
        """K(0) = π/2."""
        result = elliptic_integrals(0.0)
        assert abs(result['K_k'] - np.pi / 2) < 1e-8

    def test_E_at_zero(self):
        """E(0) = π/2."""
        result = elliptic_integrals(0.0)
        assert abs(result['E_k'] - np.pi / 2) < 1e-8

    def test_legendre_relation(self):
        """Legendre-Relation: K·E' + K'·E - K·K' = π/2."""
        result = elliptic_integrals(0.5)
        assert result['legendre_relation_holds']

    def test_numerical_vs_scipy(self):
        """Numerische Integration stimmt mit scipy überein."""
        result = elliptic_integrals(0.7)
        assert abs(result['K_k'] - result['K_numerical']) < 1e-5
        assert abs(result['E_k'] - result['E_numerical']) < 1e-5

    def test_K_increases_with_k(self):
        """K(k) ist monoton wachsend in k."""
        K1 = elliptic_integrals(0.3)['K_k']
        K2 = elliptic_integrals(0.7)['K_k']
        assert K1 < K2

    def test_invalid_k(self):
        """k ≥ 1 oder k < 0 löst ValueError aus."""
        with pytest.raises(ValueError):
            elliptic_integrals(1.0)
        with pytest.raises(ValueError):
            elliptic_integrals(-0.1)


# ---------------------------------------------------------------------------
# Fehlerfunktion Tests
# ---------------------------------------------------------------------------

class TestErrorFunction:
    """Tests für die Fehlerfunktion."""

    def test_erf_at_zero(self):
        """erf(0) = 0."""
        result = error_function_properties()
        assert abs(result['values']['erf_0']) < 1e-10

    def test_erf_at_infinity(self):
        """erf(∞) = 1."""
        result = error_function_properties()
        assert abs(result['values']['erf_inf'] - 1.0) < 1e-6

    def test_erfc_at_zero(self):
        """erfc(0) = 1."""
        result = error_function_properties()
        assert abs(result['values']['erfc_0'] - 1.0) < 1e-10

    def test_erf_erfc_complement(self):
        """erf(x) + erfc(x) = 1."""
        result = error_function_properties()
        erf1 = result['values']['erf_1']
        erfc1 = result['values']['erfc_1']
        assert abs(erf1 + erfc1 - 1.0) < 1e-10

    def test_symmetry(self):
        """erf ist antisymmetrisch: erf(-x) = -erf(x)."""
        result = error_function_properties()
        assert result['symmetry_holds']

    def test_sigma_probabilities(self):
        """1-Sigma: P(|X|≤1) ≈ 0.6827."""
        result = error_function_properties()
        prob_1sigma = result['sigma_probabilities']['1_sigma']
        assert abs(prob_1sigma - 0.6827) < 0.001

    def test_3_sigma_probability(self):
        """3-Sigma: P(|X|≤3) ≈ 0.9973."""
        result = error_function_properties()
        prob_3sigma = result['sigma_probabilities']['3_sigma']
        assert abs(prob_3sigma - 0.9973) < 0.001

    def test_erfinv_roundtrip(self):
        """erfinv(erf(x)) = x."""
        result = error_function_properties()
        assert abs(result['erfinv_values']['erfinv_erf_roundtrip'] - 1.23) < 1e-6

    def test_asymptotic_error(self):
        """Asymptotik für x=5 hat relativen Fehler < 10%."""
        result = error_function_properties()
        assert result['asymptotic']['relative_error'] < 0.1
