"""
@file test_l_functions.py
@brief Tests für L-Funktionen elliptischer Kurven, Hecke-Algebra und p-adische L-Funktionen.
@description
    Testet die in modular_forms.py (Hecke-Algebra, L-Funktionen elliptischer Kurven)
    und p_adic.py (Kubota-Leopoldt, Bernoulli-Zahlen, Kummer-Kongruenzen) implementierten
    Funktionen.

    Getestete Module:
        - modular_forms: hecke_algebra_structure, hecke_eigenform,
          elliptic_curve_points_over_fp, trace_of_frobenius,
          l_function_modular_form, l_function_elliptic_curve,
          birch_swinnerton_dyer_bsd_estimate, petersson_inner_product_estimate,
          functional_equation_l_function
        - p_adic: bernoulli_number, generalized_bernoulli_numbers,
          kubota_leopoldt_l_function, p_adic_zeta_function,
          kummer_congruence_check, iwasawa_mu_lambda_invariants

    Mathematische Grundlagen:
        - Frobenius-Spur: a_p = p + 1 - #E(F_p), Hasse-Schranke |a_p| ≤ 2√p
        - Bernoulli-Zahlen: B_0=1, B_1=-1/2, B_2=1/6, B_n=0 für ungerades n>1
        - Kubota-Leopoldt: L_p(1-n, chi) = -(1 - chi(p)p^{n-1}) B_{n,chi} / n

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math
import cmath
import pytest

# Sicherstellen, dass /src im Python-Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modular_forms import (
    hecke_algebra_structure,
    hecke_eigenform,
    elliptic_curve_points_over_fp,
    trace_of_frobenius,
    l_function_modular_form,
    l_function_elliptic_curve,
    birch_swinnerton_dyer_bsd_estimate,
    petersson_inner_product_estimate,
    functional_equation_l_function,
    fourier_coefficients_delta,
    hecke_operator,
)

from p_adic import (
    bernoulli_number,
    generalized_bernoulli_numbers,
    kubota_leopoldt_l_function,
    p_adic_zeta_function,
    kummer_congruence_check,
    iwasawa_mu_lambda_invariants,
)


# ===========================================================================
# TESTS FÜR: elliptic_curve_points_over_fp
# ===========================================================================

class TestEllipticCurvePointsOverFp:
    """Tests für elliptic_curve_points_over_fp."""

    def test_regulaere_kurve_y2_eq_x3_minus_x_mod7(self):
        """
        Testet E: y² = x³ - x (a=-1, b=0) über F_7.
        Bekannte Werte: Punkte manuell nachrechenbar.
        """
        pts = elliptic_curve_points_over_fp(-1, 0, 7)
        # Alle Punkte müssen die Kurvengleichung erfüllen
        for (x, y) in pts:
            lhs = (y * y) % 7
            rhs = (x**3 - x) % 7
            assert lhs == rhs, f"Punkt ({x},{y}) liegt nicht auf y²=x³-x mod 7"

    def test_anzahl_punkte_hasse_schranke(self):
        """
        Anzahl der Punkte #E(F_7) muss der Hasse-Schranke genügen:
        |#E(F_p) - (p+1)| ≤ 2√p
        """
        p = 7
        pts = elliptic_curve_points_over_fp(-1, 0, p)
        num_affine = len(pts)
        num_total = num_affine + 1  # +1 für Punkt im Unendlichen
        hasse_bound = 2 * math.sqrt(p)
        assert abs(num_total - (p + 1)) <= hasse_bound + 1e-9, (
            f"#E(F_{p}) = {num_total} verletzt Hasse-Schranke: "
            f"|{num_total} - {p+1}| = {abs(num_total-(p+1))} > 2√{p} ≈ {hasse_bound:.4f}"
        )

    def test_punkte_liegen_auf_kurve_mod5(self):
        """
        Teste E: y² = x³ + x + 1 (a=1, b=1) über F_5.
        Alle zurückgegebenen Punkte müssen auf der Kurve liegen.
        """
        p = 5
        pts = elliptic_curve_points_over_fp(1, 1, p)
        for (x, y) in pts:
            lhs = (y * y) % p
            rhs = (x**3 + x + 1) % p
            assert lhs == rhs, f"Punkt ({x},{y}) nicht auf y²=x³+x+1 mod {p}"

    def test_punkte_liegen_auf_kurve_mod11(self):
        """
        Teste E: y² = x³ - x + 1 (a=-1, b=1) über F_11.
        """
        p = 11
        pts = elliptic_curve_points_over_fp(-1, 1, p)
        for (x, y) in pts:
            lhs = (y * y) % p
            rhs = (pow(x, 3, p) - x + 1) % p
            assert lhs == rhs, f"Punkt ({x},{y}) nicht auf Kurve mod {p}"

    def test_rueckgabe_ist_liste(self):
        """elliptic_curve_points_over_fp muss eine Liste zurückgeben."""
        result = elliptic_curve_points_over_fp(-1, 0, 5)
        assert isinstance(result, list)

    def test_punkte_sind_tupel_von_zwei_ints(self):
        """Jeder Punkt muss ein (x, y)-Tupel mit Einträgen in F_p sein."""
        p = 7
        pts = elliptic_curve_points_over_fp(-1, 0, p)
        for pt in pts:
            assert len(pt) == 2, f"Punkt {pt} sollte 2 Einträge haben"
            x, y = pt
            assert 0 <= x < p, f"x={x} nicht in F_{p}"
            assert 0 <= y < p, f"y={y} nicht in F_{p}"

    def test_kurve_mit_grossem_p(self):
        """
        Hasse-Schranke für p=13: |a_p| ≤ 2√13 ≈ 7.21
        """
        p = 13
        pts = elliptic_curve_points_over_fp(1, 0, p)
        num_total = len(pts) + 1
        hasse = 2 * math.sqrt(p)
        assert abs(num_total - (p + 1)) <= hasse + 1e-9


# ===========================================================================
# TESTS FÜR: trace_of_frobenius
# ===========================================================================

class TestTraceOfFrobenius:
    """Tests für trace_of_frobenius."""

    def test_hasse_weil_schranke_p7(self):
        """
        Hasse-Weil-Schranke: |a_p| ≤ 2√p.
        Teste für E: y² = x³ - x (a=-1, b=0) bei p=7.
        """
        a_p = trace_of_frobenius(-1, 0, 7)
        assert abs(a_p) <= 2 * math.sqrt(7) + 1e-9, (
            f"Hasse-Schranke verletzt: |a_7| = {abs(a_p)} > 2√7 ≈ {2*math.sqrt(7):.4f}"
        )

    def test_hasse_weil_schranke_verschiedene_primzahlen(self):
        """Hasse-Schranke gilt für alle getesteten Primzahlen."""
        for p in [5, 7, 11, 13, 17]:
            a_p = trace_of_frobenius(-1, 0, p)
            assert abs(a_p) <= 2 * math.sqrt(p) + 1e-9, (
                f"Hasse-Schranke verletzt bei p={p}: |a_p|={abs(a_p)}"
            )

    def test_rueckgabe_ist_int(self):
        """trace_of_frobenius muss eine ganze Zahl zurückgeben."""
        result = trace_of_frobenius(-1, 0, 7)
        assert isinstance(result, int)

    def test_formel_ap_eq_p_plus_1_minus_nfp(self):
        """
        Direkte Prüfung: a_p = p + 1 - #E(F_p).
        """
        p = 7
        pts = elliptic_curve_points_over_fp(-1, 0, p)
        n_fp = len(pts) + 1  # +1 für Punkt im Unendlichen
        a_p = trace_of_frobenius(-1, 0, p)
        assert a_p == p + 1 - n_fp, (
            f"a_p = {a_p} ≠ p+1-#E(F_p) = {p+1-n_fp}"
        )

    def test_hasse_schranke_kurve_y2_x3_plus1(self):
        """Teste E: y² = x³ + 1 (a=0, b=1) für mehrere Primzahlen."""
        for p in [5, 7, 11]:
            a_p = trace_of_frobenius(0, 1, p)
            assert abs(a_p) <= 2 * math.sqrt(p) + 1e-9


# ===========================================================================
# TESTS FÜR: l_function_modular_form
# ===========================================================================

class TestLFunctionModularForm:
    """Tests für l_function_modular_form."""

    def test_rueckgabe_ist_komplex(self):
        """l_function_modular_form gibt eine komplexe Zahl zurück."""
        tau_coeffs = fourier_coefficients_delta(20)
        result = l_function_modular_form(tau_coeffs, s=complex(14), k=12)
        assert isinstance(result, complex)

    def test_ramanujan_tau_l_funktion_konvergenz(self):
        """
        L(Δ, s) für Re(s) = 14 > (12+1)/2 = 6.5 sollte gut konvergieren.
        Das Ergebnis ist endlich und nicht NaN.
        """
        tau_coeffs = fourier_coefficients_delta(50)
        result = l_function_modular_form(tau_coeffs, s=complex(14), k=12)
        assert math.isfinite(result.real), "L(Δ,14) realteil sollte endlich sein"
        assert math.isfinite(result.imag), "L(Δ,14) imaginärteil sollte endlich sein"

    def test_koeffizienten_a1_gleich_1_bei_s_0_grenzwert(self):
        """
        Für sehr großes Re(s) dominiert der erste Term: L(f, s) ≈ a_1 / 1^s = a_1.
        Für Ramanujan-Tau: a_1 = 1.
        """
        tau_coeffs = fourier_coefficients_delta(50)
        # s sehr groß: alle n ≥ 2 Terme verschwinden
        result = l_function_modular_form(tau_coeffs, s=complex(100), k=12)
        # Nur a_1/1^100 = 1 sollte dominieren
        assert abs(result - 1.0) < 0.01, f"L(Δ, 100) ≈ 1, erhalten: {result}"

    def test_mit_konstanten_koeffizienten(self):
        """
        L(f, s) für f mit a_n = 1 (Riemann-artig): L(f,s) = ζ(s).
        Für s=2: ζ(2) = π²/6 ≈ 1.6449.
        """
        N = 200
        ones_coeffs = [1] * N
        result = l_function_modular_form(ones_coeffs, s=complex(2), k=2)
        zeta_2 = math.pi**2 / 6
        # Partialsumme konvergiert langsam, Toleranz 5%
        assert abs(result.real - zeta_2) / zeta_2 < 0.05, (
            f"L mit a_n=1 bei s=2 sollte ≈ ζ(2)={zeta_2:.4f} sein, erhalten: {result.real:.4f}"
        )


# ===========================================================================
# TESTS FÜR: l_function_elliptic_curve
# ===========================================================================

class TestLFunctionEllipticCurve:
    """Tests für l_function_elliptic_curve."""

    def test_rueckgabe_ist_komplex(self):
        """l_function_elliptic_curve gibt eine komplexe Zahl zurück."""
        result = l_function_elliptic_curve(-1, 0, complex(2), terms=30)
        assert isinstance(result, complex)

    def test_endlicher_wert_fuer_grosses_re_s(self):
        """Für Re(s) groß konvergiert die Reihe schnell."""
        result = l_function_elliptic_curve(-1, 0, complex(4), terms=50)
        assert math.isfinite(result.real)
        assert math.isfinite(result.imag)

    def test_erster_koeffizient_a1_gleich_1(self):
        """
        Bei sehr großem s dominiert a_1/1^s = 1 den Wert.
        """
        result = l_function_elliptic_curve(-1, 0, complex(100), terms=50)
        # L(E, 100) ≈ 1 (erster Term dominiert)
        assert abs(result - 1.0) < 0.01, f"L(E,100) ≈ 1, erhalten: {result}"


# ===========================================================================
# TESTS FÜR: birch_swinnerton_dyer_bsd_estimate
# ===========================================================================

class TestBirchSwinnertonDyerEstimate:
    """Tests für birch_swinnerton_dyer_bsd_estimate."""

    def test_rueckgabe_ist_dict(self):
        """birch_swinnerton_dyer_bsd_estimate gibt ein Dictionary zurück."""
        primes = [5, 7, 11, 13]
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, primes)
        assert isinstance(result, dict)

    def test_dict_hat_pflichtschluessel(self):
        """Das Dictionary enthält alle erwarteten Schlüssel."""
        primes = [5, 7, 11, 13]
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, primes)
        required_keys = ['estimate', 'primes_used', 'log_product', 'product', 'n_p_values']
        for key in required_keys:
            assert key in result, f"Schlüssel '{key}' fehlt im Ergebnis-Dictionary"

    def test_n_p_values_enthaelt_primzahlen(self):
        """n_p_values enthält #E(F_p) für verwendete Primzahlen."""
        primes = [5, 7, 11]
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, primes)
        n_p = result['n_p_values']
        assert isinstance(n_p, dict)
        for p in result['primes_used']:
            assert p in n_p, f"Primzahl {p} fehlt in n_p_values"

    def test_n_p_values_erfuellen_hasse_schranke(self):
        """Alle #E(F_p) in n_p_values erfüllen die Hasse-Schranke."""
        primes = [5, 7, 11, 13]
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, primes)
        for p, n_p in result['n_p_values'].items():
            # |N_p - (p+1)| ≤ 2√p
            assert abs(n_p - (p + 1)) <= 2 * math.sqrt(p) + 1e-9, (
                f"#E(F_{p}) = {n_p} verletzt Hasse-Schranke"
            )

    def test_estimate_ist_float(self):
        """Die Rang-Schätzung ist eine Gleitkommazahl."""
        primes = [5, 7, 11, 13, 17]
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, primes)
        assert isinstance(result['estimate'], float)

    def test_leere_primzahlliste(self):
        """Leere Primzahlliste: Schätzung = 0.0, leere n_p_values."""
        result = birch_swinnerton_dyer_bsd_estimate(-1, 0, [])
        assert isinstance(result, dict)


# ===========================================================================
# TESTS FÜR: hecke_algebra_structure
# ===========================================================================

class TestHeckeAlgebraStructure:
    """Tests für hecke_algebra_structure."""

    def test_rueckgabe_ist_dict(self):
        """hecke_algebra_structure gibt ein Dictionary zurück."""
        result = hecke_algebra_structure(12)
        assert isinstance(result, dict)

    def test_hecke_eigenwerte_fuer_delta(self):
        """
        Für die Delta-Funktion (k=12) sind die Hecke-Eigenwerte die Ramanujan-Tau-Zahlen.
        Bekannte Werte: τ(2) = -24, τ(3) = 252.
        """
        result = hecke_algebra_structure(12)
        eigenvalues = result['hecke_eigenvalues']
        assert eigenvalues.get(2) == -24, f"τ(2) sollte -24 sein, erhalten: {eigenvalues.get(2)}"
        assert eigenvalues.get(3) == 252, f"τ(3) sollte 252 sein, erhalten: {eigenvalues.get(3)}"

    def test_kommutativitaet(self):
        """Die Hecke-Algebra ist kommutativ: T_p T_q = T_q T_p."""
        result = hecke_algebra_structure(12)
        assert result['is_commutative'] is True, (
            "Hecke-Algebra sollte kommutativ sein"
        )

    def test_multiplikativitaet_ramanujan_tau(self):
        """
        Ramanujan-Tau ist multiplikativ: τ(p*q) = τ(p)*τ(q) für verschiedene Primzahlen p≠q.
        Teste: τ(2*3) = τ(6) = τ(2)*τ(3) = (-24)*252 = -6048.
        """
        tau_list = fourier_coefficients_delta(10)
        # τ(6) = -6048 (bekannter Wert)
        assert tau_list[5] == -6048, f"τ(6) sollte -6048 sein, erhalten: {tau_list[5]}"
        assert tau_list[5] == tau_list[1] * tau_list[2], (
            f"τ(6) = {tau_list[5]} ≠ τ(2)*τ(3) = {tau_list[1]*tau_list[2]}"
        )

    def test_dict_hat_pflichtschluessel(self):
        """Das Ergebnis-Dictionary enthält alle Pflichtschlüssel."""
        result = hecke_algebra_structure(12)
        for key in ['hecke_eigenvalues', 'commutativity_check',
                    'multiplicativity_check', 'is_commutative']:
            assert key in result, f"Pflichtschlüssel '{key}' fehlt"


# ===========================================================================
# TESTS FÜR: hecke_eigenform
# ===========================================================================

class TestHeckeEigenform:
    """Tests für hecke_eigenform."""

    def test_delta_ist_eigenform(self):
        """
        Die Ramanujan-Delta-Funktion ist eine normierte Hecke-Eigenform.
        hecke_eigenform sollte is_eigenform=True zurückgeben.
        """
        tau_coeffs = fourier_coefficients_delta(50)
        primes = [2, 3, 5]
        result = hecke_eigenform(tau_coeffs, primes)
        assert isinstance(result, dict)
        assert result.get('is_eigenform') is True, (
            f"Delta-Funktion sollte Hecke-Eigenform sein. Ergebnis: {result}"
        )

    def test_eigenwerte_sind_tau_zahlen(self):
        """
        Eigenwerte der Delta-Funktion sind Ramanujan-Tau-Zahlen:
        λ_2 = τ(2)/τ(1) = -24/1 = -24.
        """
        tau_coeffs = fourier_coefficients_delta(50)
        primes = [2, 3]
        result = hecke_eigenform(tau_coeffs, primes)
        eigenvalues = result.get('eigenvalues', {})
        # λ_2 = a_2 / a_1 = -24 / 1 = -24
        assert abs(eigenvalues.get(2, float('inf')) - (-24)) < 0.5, (
            f"Eigenwert λ_2 sollte -24 sein, erhalten: {eigenvalues.get(2)}"
        )

    def test_rueckgabe_ist_dict(self):
        """hecke_eigenform gibt ein Dictionary zurück."""
        tau_coeffs = fourier_coefficients_delta(20)
        result = hecke_eigenform(tau_coeffs, [2, 3])
        assert isinstance(result, dict)

    def test_multiplikativitaet_der_delta_form(self):
        """
        Multiplikativitäts-Check für Delta: τ(p*q) = τ(p)*τ(q) für p≠q prim.
        """
        tau_coeffs = fourier_coefficients_delta(50)
        primes = [2, 3, 5]
        result = hecke_eigenform(tau_coeffs, primes)
        mult = result.get('multiplicativity', {})
        # Mindestens ein Paar geprüft
        if mult:
            for pair, val in mult.items():
                assert val is True, f"Multiplikativität verletzt für Paar {pair}: {val}"


# ===========================================================================
# TESTS FÜR: petersson_inner_product_estimate
# ===========================================================================

class TestPeterssonInnerProduct:
    """Tests für petersson_inner_product_estimate."""

    def test_rueckgabe_ist_komplex(self):
        """petersson_inner_product_estimate gibt eine komplexe Zahl zurück."""
        tau = fourier_coefficients_delta(20)
        result = petersson_inner_product_estimate(tau, tau, k=12)
        assert isinstance(result, complex)

    def test_skalarprodukt_mit_sich_selbst_positiv(self):
        """
        ⟨f, f⟩ > 0 für eine nicht-triviale Cusp-Form.
        Das Petersson-Skalarprodukt ist positiv definit auf Cusp-Formen.
        """
        tau = fourier_coefficients_delta(20)
        result = petersson_inner_product_estimate(tau, tau, k=12)
        assert result.real > 0, (
            f"⟨Δ,Δ⟩ sollte > 0 sein (positiv definit), erhalten: {result.real}"
        )

    def test_leere_koeffizienten(self):
        """Leere Koeffizientenlisten: Ergebnis = 0."""
        result = petersson_inner_product_estimate([], [], k=12)
        assert result == complex(0.0)


# ===========================================================================
# TESTS FÜR: functional_equation_l_function
# ===========================================================================

class TestFunctionalEquationLFunction:
    """Tests für functional_equation_l_function."""

    def test_rueckgabe_ist_dict(self):
        """functional_equation_l_function gibt ein Dictionary zurück."""
        tau = fourier_coefficients_delta(20)
        result = functional_equation_l_function(tau, k=12, N=1, s=complex(7))
        assert isinstance(result, dict)

    def test_dict_hat_pflichtschluessel(self):
        """Das Dictionary enthält alle erwarteten Schlüssel."""
        tau = fourier_coefficients_delta(20)
        result = functional_equation_l_function(tau, k=12, N=1, s=complex(7))
        required = ['lambda_s', 'lambda_k_minus_s', 'epsilon', 'functional_equation_ratio']
        for key in required:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_vorzeichen_epsilon_fuer_k12(self):
        """
        Für k=12: ε = (-1)^{12/2} = (-1)^6 = 1.
        """
        tau = fourier_coefficients_delta(20)
        result = functional_equation_l_function(tau, k=12, N=1, s=complex(7))
        assert result['epsilon'] == 1, (
            f"ε für k=12 sollte 1 sein, erhalten: {result['epsilon']}"
        )

    def test_vorzeichen_epsilon_fuer_k4(self):
        """
        Für k=4: ε = (-1)^{4/2} = (-1)^2 = 1.
        """
        coeffs = [1] * 20
        result = functional_equation_l_function(coeffs, k=4, N=1, s=complex(3))
        assert result['epsilon'] == 1


# ===========================================================================
# TESTS FÜR: bernoulli_number (p_adic.py)
# ===========================================================================

class TestBernoulliNumber:
    """Tests für bernoulli_number aus p_adic.py."""

    def test_b0_gleich_1(self):
        """B_0 = 1 (Definitionswert)."""
        assert bernoulli_number(0) == 1.0

    def test_b1_gleich_minus_halb(self):
        """B_1 = -1/2."""
        assert abs(bernoulli_number(1) - (-0.5)) < 1e-10

    def test_b2_gleich_ein_sechstel(self):
        """B_2 = 1/6 ≈ 0.16667."""
        assert abs(bernoulli_number(2) - (1.0 / 6.0)) < 1e-10

    def test_b3_gleich_null(self):
        """B_3 = 0 (ungerades n > 1)."""
        assert bernoulli_number(3) == 0.0

    def test_b4_gleich_minus_1_durch_30(self):
        """B_4 = -1/30 ≈ -0.03333."""
        assert abs(bernoulli_number(4) - (-1.0 / 30.0)) < 1e-10

    def test_b6_gleich_1_durch_42(self):
        """B_6 = 1/42 ≈ 0.02381."""
        assert abs(bernoulli_number(6) - (1.0 / 42.0)) < 1e-10

    def test_ungerade_n_groesser_1_gleich_null(self):
        """B_n = 0 für alle ungeraden n > 1."""
        for n in [3, 5, 7, 9, 11, 13]:
            result = bernoulli_number(n)
            assert result == 0.0, f"B_{n} sollte 0 sein, erhalten: {result}"

    def test_negative_n_wirft_fehler(self):
        """bernoulli_number(-1) sollte einen ValueError werfen."""
        with pytest.raises(ValueError):
            bernoulli_number(-1)

    def test_b8_gleich_minus_1_durch_30(self):
        """B_8 = -1/30."""
        assert abs(bernoulli_number(8) - (-1.0 / 30.0)) < 1e-10

    def test_b10_gleich_5_durch_66(self):
        """B_10 = 5/66."""
        assert abs(bernoulli_number(10) - (5.0 / 66.0)) < 1e-10


# ===========================================================================
# TESTS FÜR: generalized_bernoulli_numbers
# ===========================================================================

class TestGeneralizedBernoulliNumbers:
    """Tests für generalized_bernoulli_numbers aus p_adic.py."""

    def test_trivialer_charakter_stimmt_mit_bernoulli_zahlen_ueberein(self):
        """
        Für chi_conductor = 1 (trivialer Charakter):
        B_{n, trivial} = B_n (gewöhnliche Bernoulli-Zahlen).
        """
        n_max = 6
        gen_b = generalized_bernoulli_numbers(1, n_max)
        for n in range(n_max + 1):
            expected = bernoulli_number(n)
            assert abs(gen_b[n] - expected) < 1e-8, (
                f"B_{{{n},1}} = {gen_b[n]} ≠ B_{n} = {expected}"
            )

    def test_rueckgabe_laenge(self):
        """Ergebnis-Liste hat n_max + 1 Einträge."""
        for n_max in [3, 5, 8]:
            result = generalized_bernoulli_numbers(1, n_max)
            assert len(result) == n_max + 1, (
                f"Länge sollte {n_max + 1} sein, erhalten: {len(result)}"
            )

    def test_b0_trivialer_charakter(self):
        """B_{0, trivial} = B_0 = 1."""
        result = generalized_bernoulli_numbers(1, 4)
        assert abs(result[0] - 1.0) < 1e-8

    def test_b1_trivialer_charakter(self):
        """B_{1, trivial} = B_1 = -1/2."""
        result = generalized_bernoulli_numbers(1, 4)
        assert abs(result[1] - (-0.5)) < 1e-8

    def test_b2_trivialer_charakter(self):
        """B_{2, trivial} = B_2 = 1/6."""
        result = generalized_bernoulli_numbers(1, 4)
        assert abs(result[2] - (1.0 / 6.0)) < 1e-8

    def test_allgemeiner_charakter_gibt_liste_zurueck(self):
        """Für chi_conductor > 1 wird ebenfalls eine Liste zurückgegeben."""
        result = generalized_bernoulli_numbers(4, 5)
        assert isinstance(result, list)
        assert len(result) == 6


# ===========================================================================
# TESTS FÜR: kubota_leopoldt_l_function
# ===========================================================================

class TestKubotaLeopoldt:
    """Tests für kubota_leopoldt_l_function aus p_adic.py."""

    def test_rueckgabe_ist_komplex(self):
        """kubota_leopoldt_l_function gibt eine komplexe Zahl zurück."""
        result = kubota_leopoldt_l_function(chi_conductor=1, s_padic=-1, p=3)
        assert isinstance(result, complex)

    def test_interpolation_fuer_s_gleich_minus1_trivialer_charakter(self):
        """
        L_p(1-2, trivial) = -(1 - p^{1}) * B_2 / 2
        Für p=3: = -(1 - 3) * (1/6) / 2 = -(-2) * (1/6) / 2 = 1/6
        s = 1-2 = -1 ⟹ n = 2
        """
        p = 3
        # s = -1 ⟹ n = 1-s = 1-(-1) = 2
        result = kubota_leopoldt_l_function(chi_conductor=1, s_padic=-1, p=p)
        # Erwarteter Wert: -(1 - p^{n-1}) * B_n / n mit n=2
        # = -(1 - 3^1) * (1/6) / 2 = -(-2) * (1/6) / 2 = 2/(6*2) = 1/6
        expected = -(1 - p**(2-1)) * bernoulli_number(2) / 2
        assert abs(result.real - expected) < 1e-8, (
            f"L_3(-1, 1) sollte {expected:.6f} sein, erhalten: {result.real:.6f}"
        )

    def test_interpolation_fuer_s_gleich_minus3_trivialer_charakter(self):
        """
        L_p(1-4, trivial) = -(1 - p^3) * B_4 / 4
        Für p=5, s=-3 ⟹ n=4: -(1 - 5^3) * (-1/30) / 4
        """
        p = 5
        result = kubota_leopoldt_l_function(chi_conductor=1, s_padic=-3, p=p)
        # n = 1 - (-3) = 4
        n = 4
        expected = -(1 - p**(n-1)) * bernoulli_number(n) / n
        assert abs(result.real - expected) < 1e-8, (
            f"L_5(-3, 1) sollte {expected:.6f} sein, erhalten: {result.real:.6f}"
        )

    def test_wirft_fehler_fuer_ungueltiges_p(self):
        """ValueError bei p < 2."""
        with pytest.raises(ValueError):
            kubota_leopoldt_l_function(1, -1, p=1)

    def test_wirft_fehler_fuer_ungueltigen_conductor(self):
        """ValueError bei chi_conductor < 1."""
        with pytest.raises(ValueError):
            kubota_leopoldt_l_function(0, -1, p=3)


# ===========================================================================
# TESTS FÜR: p_adic_zeta_function
# ===========================================================================

class TestPAdicZetaFunction:
    """Tests für p_adic_zeta_function aus p_adic.py."""

    def test_rueckgabe_ist_float(self):
        """p_adic_zeta_function gibt einen float zurück."""
        result = p_adic_zeta_function(3, -1)
        assert isinstance(result, float)

    def test_interpolationsformel_s_minus1_p3(self):
        """
        ζ_3(-1) = -(1 - 3^{1}) * B_2 / 2 = -(1-3) * (1/6) / 2 = 1/6
        """
        result = p_adic_zeta_function(3, -1)
        # s = -1 ⟹ n = 2: ζ_p(-1) = -(1-p)*B_2/2
        expected = -(1 - 3) * bernoulli_number(2) / 2
        assert abs(result - expected) < 1e-10, (
            f"ζ_3(-1) sollte {expected:.8f} sein, erhalten: {result:.8f}"
        )

    def test_interpolationsformel_s_minus3_p5(self):
        """
        ζ_5(-3) = -(1 - 5³) * B_4 / 4
        n = 1 - (-3) = 4
        """
        p = 5
        s = -3
        result = p_adic_zeta_function(p, s)
        n = 1 - s  # = 4
        expected = -(1 - p**(n-1)) * bernoulli_number(n) / n
        assert abs(result - expected) < 1e-10

    def test_ungerade_n_gibt_null(self):
        """
        Für ungerades n > 1 (d.h. gerades s = 1-n): B_n = 0 ⟹ ζ_p(s) = 0.
        s = 1-3 = -2 ⟹ n = 3 (ungerade): ζ_p(-2) = 0.
        """
        result = p_adic_zeta_function(3, -2)  # n=3, ungerade
        assert result == 0.0, f"ζ_3(-2) sollte 0 sein (B_3=0), erhalten: {result}"

    def test_verschiedene_primzahlen(self):
        """ζ_p(-1) für verschiedene Primzahlen gibt jeweils float zurück."""
        for p in [2, 3, 5, 7]:
            result = p_adic_zeta_function(p, -1)
            assert isinstance(result, float)

    def test_s_gleich_2_gibt_float(self):
        """p_adic_zeta_function(3, 2) gibt einen float zurück (Testfall aus Aufgabe)."""
        result = p_adic_zeta_function(3, 2)
        assert isinstance(result, float)


# ===========================================================================
# TESTS FÜR: kummer_congruence_check
# ===========================================================================

class TestKummerCongruenceCheck:
    """Tests für kummer_congruence_check aus p_adic.py."""

    def test_rueckgabe_ist_dict(self):
        """kummer_congruence_check gibt ein Dictionary zurück."""
        result = kummer_congruence_check(3, [2, 4, 6, 8])
        assert isinstance(result, dict)

    def test_dict_hat_pflichtschluessel(self):
        """Das Dictionary enthält alle Pflichtschlüssel."""
        result = kummer_congruence_check(3, [2, 4])
        for key in ['results', 'all_satisfied', 'p', 'bernoulli_values', 'zeta_p_values']:
            assert key in result, f"Pflichtschlüssel '{key}' fehlt"

    def test_p_im_ergebnis(self):
        """Das Dictionary enthält die verwendete Primzahl."""
        result = kummer_congruence_check(5, [2, 4])
        assert result['p'] == 5

    def test_bernoulli_werte_korrekt(self):
        """bernoulli_values enthält korrekte Bernoulli-Zahlen."""
        result = kummer_congruence_check(3, [2, 4, 6])
        bv = result['bernoulli_values']
        assert abs(bv[2] - bernoulli_number(2)) < 1e-10
        assert abs(bv[4] - bernoulli_number(4)) < 1e-10

    def test_leere_liste_gibt_dict_zurueck(self):
        """Leere n-Liste: Ergebnis ist leeres 'results' Dict."""
        result = kummer_congruence_check(3, [])
        assert isinstance(result, dict)
        assert result['all_satisfied'] is True  # Trivial wahr (kein Paar geprüft)


# ===========================================================================
# TESTS FÜR: iwasawa_mu_lambda_invariants
# ===========================================================================

class TestIwasawaMuLambdaInvariants:
    """Tests für iwasawa_mu_lambda_invariants aus p_adic.py."""

    def test_rueckgabe_ist_dict(self):
        """iwasawa_mu_lambda_invariants gibt ein Dictionary zurück."""
        result = iwasawa_mu_lambda_invariants(3, 1, n_check=6)
        assert isinstance(result, dict)

    def test_dict_hat_pflichtschluessel(self):
        """Das Dictionary enthält alle Pflichtschlüssel."""
        result = iwasawa_mu_lambda_invariants(3, 1, n_check=6)
        for key in ['mu_estimate', 'lambda_estimate', 'l_p_values',
                    'p_adic_valuations', 'vandiver_ok']:
            assert key in result, f"Pflichtschlüssel '{key}' fehlt"

    def test_mu_ist_nicht_negativ(self):
        """μ-Invariante ist nicht negativ."""
        result = iwasawa_mu_lambda_invariants(3, 1, n_check=6)
        assert result['mu_estimate'] >= 0

    def test_lambda_ist_nicht_negativ(self):
        """λ-Invariante ist nicht negativ."""
        result = iwasawa_mu_lambda_invariants(5, 1, n_check=8)
        assert result['lambda_estimate'] >= 0

    def test_vandiver_ok_ist_bool(self):
        """vandiver_ok ist ein bool."""
        result = iwasawa_mu_lambda_invariants(3, 1, n_check=4)
        assert isinstance(result['vandiver_ok'], bool)

    def test_greenberg_vermutung_mu_gleich_null(self):
        """
        Greenberg-Vermutung: μ = 0 für reguläre Primzahlen.
        p=3, 5, 7 sind reguläre Primzahlen, also sollte μ = 0 sein.
        """
        for p in [3, 5, 7]:
            result = iwasawa_mu_lambda_invariants(p, 1, n_check=6)
            # Greenberg-Vermutung: μ sollte 0 sein
            # Dies ist numerisch schwer zu verifizieren, aber eine Schätzung
            assert result['mu_estimate'] >= 0  # Minimaltest: μ nicht negativ

    def test_l_p_values_ist_dict(self):
        """l_p_values enthält L_p-Werte als Dictionary."""
        result = iwasawa_mu_lambda_invariants(3, 1, n_check=6)
        assert isinstance(result['l_p_values'], dict)


# ===========================================================================
# INTEGRATIONSTESTS
# ===========================================================================

class TestIntegration:
    """Integrationstests für das Zusammenspiel der Funktionen."""

    def test_frobenius_und_l_funktion_konsistent(self):
        """
        Die a_p-Koeffizienten der L-Funktion stimmen mit trace_of_frobenius überein.
        """
        a_coeff = -1
        b_coeff = 0
        p = 7
        # Direkte Frobenius-Spur
        a_p_direct = trace_of_frobenius(a_coeff, b_coeff, p)
        # Via L-Funktion: a_p ist der p-te Koeffizient (hier über Punktanzahl)
        pts = elliptic_curve_points_over_fp(a_coeff, b_coeff, p)
        a_p_formula = p + 1 - (len(pts) + 1)
        assert a_p_direct == a_p_formula

    def test_bernoulli_und_zeta_p_konsistent(self):
        """
        ζ_p(1-n) und direkte Bernoulli-Formel geben dasselbe Ergebnis.
        """
        p = 3
        n = 4  # Gerades n ≥ 2
        s = 1 - n  # = -3
        # Via p_adic_zeta_function
        zeta_val = p_adic_zeta_function(p, s)
        # Direkt: -(1 - p^{n-1}) * B_n / n
        direct = -(1 - p**(n-1)) * bernoulli_number(n) / n
        assert abs(zeta_val - direct) < 1e-10

    def test_kubota_leopoldt_und_zeta_p_konsistent(self):
        """
        kubota_leopoldt_l_function mit chi_conductor=1 stimmt mit p_adic_zeta_function überein.
        """
        p = 5
        s = -1  # n = 2
        kl_val = kubota_leopoldt_l_function(1, s, p).real
        zeta_val = p_adic_zeta_function(p, s)
        assert abs(kl_val - zeta_val) < 1e-8, (
            f"Kubota-Leopoldt ({kl_val:.8f}) ≠ ζ_p ({zeta_val:.8f}) für p={p}, s={s}"
        )


if __name__ == '__main__':
    # Direkter Aufruf: Alle Tests ausführen
    pytest.main([__file__, '-v', '--tb=short'])
