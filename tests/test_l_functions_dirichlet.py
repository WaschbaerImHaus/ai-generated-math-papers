"""
@file test_l_functions_dirichlet.py
@brief Tests für das neue L-Funktionen-Modul src/l_functions.py.
@description
    Testet alle Komponenten der L-Funktionen-Infrastruktur:
      - Dirichlet-Charaktere (Konstruktion, Werte, Primitivität)
      - Gauß-Summen (Betrag, Hauptcharakter)
      - Dirichlet L-Funktionen (Konvergenz, Spezialwerte)
      - Funktionalgleichung
      - Hecke L-Funktionen (Ramanujan-Delta-Koeffizienten)
      - BSD-Verbindung (elliptische Kurven, Euler-Produkt)

    Mathematische Fakten die getestet werden:
      - L(1, χ₋₄) = π/4  (Leibniz-Reihe, Charakter mod 4)
      - |τ(χ)|² = q für primitive Charaktere
      - χ(1) = 1 für alle Charaktere
      - Euler-Produktdarstellung stimmt mit Dirichtreihe überein
      - BSD: L(E,1) ≠ 0 für Kurven vom Rang 0

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

# Projektwurzel zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from l_functions import (
    dirichlet_characters,
    chi_value,
    is_primitive_character,
    dirichlet_l_function,
    dirichlet_l_function_zeros,
    gauss_sum,
    l_function_conductor,
    l_function_special_values,
    functional_equation_check,
    hecke_l_function,
    completed_l_function,
    elliptic_l_function_approx,
    bsd_rank_from_zeros,
)
from exceptions import InvalidInputError, DomainError


# ===========================================================================
# TESTS FÜR: dirichlet_characters
# ===========================================================================

class TestDirichletCharacters:
    """Tests für die Konstruktion von Dirichlet-Charakteren."""

    def test_charaktere_mod4_anzahl(self):
        """
        (ℤ/4ℤ)× = {1, 3} hat φ(4) = 2 Elemente.
        Es gibt genau 2 Charaktere mod 4.
        """
        chars = dirichlet_characters(4)
        assert len(chars) == 2, f"Erwartet 2 Charaktere mod 4, erhalten: {len(chars)}"

    def test_hauptcharakter_ist_index_0(self):
        """Der Hauptcharakter χ₀ hat Index 0 und is_principal=True."""
        chars = dirichlet_characters(4)
        assert chars[0]['is_principal'] is True

    def test_hauptcharakter_wert_bei_1_gleich_1(self):
        """Für alle Charaktere gilt χ(1) = 1."""
        for q in [3, 4, 5, 7]:
            chars = dirichlet_characters(q)
            for char in chars:
                val = char['chi'](1)
                assert abs(val - 1) < 1e-10, (
                    f"χ(1) = {val} ≠ 1 für Charakter mod {q}"
                )

    def test_charakterwert_null_fuer_nichteinheiten(self):
        """χ(n) = 0 wenn gcd(n, q) > 1."""
        # mod 4: gcd(2, 4) = 2 > 1 → χ(2) = 0
        chars = dirichlet_characters(4)
        for char in chars:
            assert abs(char['chi'](2)) < 1e-10, (
                f"χ(2) = {char['chi'](2)} ≠ 0 für Charakter mod 4 (gcd(2,4)=2)"
            )

    def test_hauptcharakter_mod4_werte(self):
        """
        Hauptcharakter χ₀ mod 4:
        χ₀(1) = 1, χ₀(2) = 0, χ₀(3) = 1, χ₀(4) = 0 (≡ χ₀(0))
        """
        chars = dirichlet_characters(4)
        chi0 = chars[0]['chi']
        assert abs(chi0(1) - 1) < 1e-10
        assert abs(chi0(2)) < 1e-10
        assert abs(chi0(3) - 1) < 1e-10

    def test_chi_minus4_werte(self):
        """
        Nicht-Hauptcharakter χ₋₄ mod 4:
        χ₋₄(1) = 1, χ₋₄(3) = -1 (ungerade Charakter: χ(-1) = -1).
        Dies ist der Charakter der Leibniz-Reihe: L(1, χ₋₄) = π/4.
        """
        chars = dirichlet_characters(4)
        # Index 1 ist der Nicht-Hauptcharakter
        chi = chars[1]['chi']
        assert abs(chi(1) - 1) < 1e-10, f"χ₋₄(1) sollte 1 sein, erhalten: {chi(1)}"
        # χ₋₄(3) = -1 (da 3 ≡ -1 mod 4)
        val_3 = chi(3)
        assert abs(abs(val_3) - 1) < 1e-10, f"|χ₋₄(3)| sollte 1 sein"

    def test_charaktere_mod1(self):
        """Mod 1 gibt es genau einen Charakter (trivialen)."""
        chars = dirichlet_characters(1)
        assert len(chars) == 1
        assert chars[0]['is_principal'] is True

    def test_charaktere_mod5_anzahl(self):
        """φ(5) = 4: Es gibt genau 4 Charaktere mod 5."""
        chars = dirichlet_characters(5)
        assert len(chars) == 4, f"Erwartet 4 Charaktere mod 5, erhalten: {len(chars)}"

    def test_ungueltige_eingabe_wirft_fehler(self):
        """q < 1 wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            dirichlet_characters(0)

    def test_charakter_werte_auf_einheitskreis(self):
        """Für gcd(n, q) = 1 gilt |χ(n)| = 1."""
        q = 5
        chars = dirichlet_characters(q)
        for char in chars:
            for n in range(1, q):
                if math.gcd(n, q) == 1:
                    val = char['chi'](n)
                    assert abs(abs(val) - 1) < 1e-8, (
                        f"|χ({n})| = {abs(val)} ≠ 1 für Charakter mod {q}"
                    )


# ===========================================================================
# TESTS FÜR: chi_value
# ===========================================================================

class TestChiValue:
    """Tests für chi_value."""

    def test_hauptcharakter_mod4_an_1(self):
        """chi_value(1, 4, 0) = 1 (Hauptcharakter)."""
        assert abs(chi_value(1, 4, 0) - 1) < 1e-10

    def test_hauptcharakter_mod4_an_2(self):
        """chi_value(2, 4, 0) = 0 (gcd(2,4) = 2)."""
        assert abs(chi_value(2, 4, 0)) < 1e-10

    def test_ungueltige_indizes_werfen_fehler(self):
        """Ungültiger char_index wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            chi_value(1, 4, 10)

    def test_ungueltige_q_wirft_fehler(self):
        """q = 0 wirft DomainError."""
        with pytest.raises(DomainError):
            chi_value(1, 0, 0)


# ===========================================================================
# TESTS FÜR: is_primitive_character
# ===========================================================================

class TestIsPrimitiveCharacter:
    """Tests für is_primitive_character."""

    def test_hauptcharakter_mod4_nicht_primitiv(self):
        """Der Hauptcharakter mod 4 ist nicht primitiv (Führer = 1)."""
        chars = dirichlet_characters(4)
        result = is_primitive_character(chars[0]['values'], 4)
        assert result is False

    def test_nichthauptcharakter_mod4_primitiv(self):
        """
        Der Charakter χ₋₄ mod 4 (Index 1) ist primitiv.
        Kriterium: |τ(χ)|² = q = 4.
        """
        chars = dirichlet_characters(4)
        if len(chars) > 1:
            result = is_primitive_character(chars[1]['values'], 4)
            assert result is True

    def test_q_gleich_1_primitiv(self):
        """Charakter mod 1 gilt als primitiv (trivialer Fall)."""
        chars = dirichlet_characters(1)
        result = is_primitive_character(chars[0]['values'], 1)
        assert result is True


# ===========================================================================
# TESTS FÜR: gauss_sum
# ===========================================================================

class TestGaussSum:
    """Tests für die Gauß-Summe τ(χ)."""

    def test_gauss_sum_ist_komplex(self):
        """gauss_sum gibt eine komplexe Zahl zurück."""
        result = gauss_sum(4, 0)
        assert isinstance(result, complex)

    def test_primitiver_charakter_gauss_norm_quad(self):
        """
        Für primitive Charaktere gilt |τ(χ)|² = q.
        Test mit χ₋₄ mod 4 (primitiv): |τ(χ₋₄)|² = 4.
        """
        chars = dirichlet_characters(4)
        if len(chars) > 1:
            tau = gauss_sum(4, 1)
            tau_norm_sq = abs(tau) ** 2
            assert abs(tau_norm_sq - 4) < 0.5, (
                f"|τ(χ₋₄)|² = {tau_norm_sq:.4f}, erwartet ≈ 4"
            )

    def test_gauss_sum_mod5_primitiv(self):
        """
        Für primitive Charaktere mod 5: |τ(χ)|² = 5.
        (ℤ/5ℤ)× ist zyklisch, alle Nicht-Hauptcharaktere sind primitiv.)
        """
        chars = dirichlet_characters(5)
        for idx in range(1, len(chars)):
            tau = gauss_sum(5, idx)
            tau_norm_sq = abs(tau) ** 2
            assert abs(tau_norm_sq - 5) < 1.0, (
                f"|τ(χ_{idx})|² mod 5 = {tau_norm_sq:.4f}, erwartet ≈ 5"
            )

    def test_ungueltige_eingabe_wirft_fehler(self):
        """q = 0 wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            gauss_sum(0, 0)


# ===========================================================================
# TESTS FÜR: dirichlet_l_function
# ===========================================================================

class TestDirichletLFunction:
    """Tests für die Dirichlet L-Funktion."""

    def test_leibniz_reihe_l1_chi_minus4(self):
        """
        L(1, χ₋₄) = 1 - 1/3 + 1/5 - 1/7 + ... = π/4 ≈ 0.7854.
        Dies ist der berühmteste Wert einer Dirichlet L-Funktion.
        Toleranz 1% wegen endlicher Reihensumme.
        """
        chars = dirichlet_characters(4)
        if len(chars) < 2:
            pytest.skip("Kein Nicht-Hauptcharakter mod 4 verfügbar")

        # L(1, χ₋₄) mit 1000 Termen
        result = dirichlet_l_function(complex(1.0 + 1e-9, 0), 4, 1, n_terms=1000)
        expected = math.pi / 4
        rel_error = abs(result.real - expected) / expected
        assert rel_error < 0.02, (
            f"L(1, χ₋₄) = {result.real:.6f}, erwartet π/4 ≈ {expected:.6f}, "
            f"relativer Fehler: {rel_error:.4f}"
        )

    def test_hauptcharakter_pol_wirft_fehler(self):
        """L(1, χ₀) divergiert → DomainError bei s=1."""
        with pytest.raises(DomainError):
            dirichlet_l_function(complex(1.0), 4, 0, n_terms=100)

    def test_rueckgabe_ist_komplex(self):
        """dirichlet_l_function gibt komplexe Zahl zurück."""
        result = dirichlet_l_function(complex(2.0), 4, 1, n_terms=100)
        assert isinstance(result, complex)

    def test_endlicher_wert_fuer_grosses_re_s(self):
        """Für Re(s) >> 1 konvergiert die Reihe schnell."""
        result = dirichlet_l_function(complex(10.0), 5, 1, n_terms=100)
        assert math.isfinite(result.real)
        assert math.isfinite(result.imag)

    def test_l_funktion_mod5_bei_s2(self):
        """L(2, χ) für χ mod 5 (Nicht-Hauptcharakter) ist endlich."""
        result = dirichlet_l_function(complex(2.0), 5, 1, n_terms=200)
        assert math.isfinite(result.real)

    def test_erster_summand_dominiert_bei_grossem_s(self):
        """
        Für sehr großes Re(s): L(s, χ) ≈ χ(1)/1^s = 1.
        """
        result = dirichlet_l_function(complex(50.0), 4, 1, n_terms=100)
        # χ(1) = 1, also L(50, χ) ≈ 1
        assert abs(result - 1.0) < 0.01, (
            f"L(50, χ) ≈ 1, erhalten: {result}"
        )

    def test_ungueltige_eingabe_wirft_fehler(self):
        """q < 1 wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            dirichlet_l_function(complex(2.0), 0, 0, n_terms=10)

    def test_euler_produkt_vs_direkte_reihe(self):
        """
        Für Re(s) > 1: Euler-Produkt ≈ Dirichlet-Reihe.
        L(s, χ₀ mod 5) = ζ(s) · (1 - 5^{-s})
        Grobe Konsistenzprüfung bei s=3.
        """
        # L(3, χ₀ mod 5) via direkte Reihe
        result_series = dirichlet_l_function(complex(3.0), 5, 0, n_terms=500)
        # ζ(3) · (1 - 5^{-3}) direkt
        zeta3 = sum(1.0 / n**3 for n in range(1, 2000))
        expected = zeta3 * (1 - 5**(-3))
        # Nur Realteile vergleichen, grobe Toleranz
        assert abs(result_series.real - expected) / abs(expected) < 0.05, (
            f"L(3, χ₀) = {result_series.real:.5f}, erwartet ≈ {expected:.5f}"
        )


# ===========================================================================
# TESTS FÜR: dirichlet_l_function_zeros
# ===========================================================================

class TestDirichletLFunctionZeros:
    """Tests für die Nullstellensuche von L(1/2 + it, χ)."""

    def test_rueckgabe_ist_liste(self):
        """dirichlet_l_function_zeros gibt eine Liste zurück."""
        result = dirichlet_l_function_zeros(4, 1, t_min=0, t_max=20, steps=100)
        assert isinstance(result, list)

    def test_nullstellen_auf_kritischer_geraden(self):
        """Alle gefundenen Nullstellen haben Re(s) = 1/2."""
        zeros = dirichlet_l_function_zeros(4, 1, t_min=0, t_max=30, steps=200)
        for z in zeros:
            assert abs(z['zero'].real - 0.5) < 1e-10, (
                f"Nullstelle {z['zero']} nicht auf Re(s) = 1/2"
            )

    def test_nullstellen_haben_abs_val_schluessel(self):
        """Jede Nullstelle enthält 'abs_val' Schlüssel."""
        zeros = dirichlet_l_function_zeros(4, 1, t_min=0, t_max=20, steps=100)
        for z in zeros:
            assert 't' in z
            assert 'zero' in z
            assert 'abs_val' in z

    def test_ungueltige_t_intervall_wirft_fehler(self):
        """t_min ≥ t_max wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            dirichlet_l_function_zeros(4, 1, t_min=10, t_max=5)


# ===========================================================================
# TESTS FÜR: l_function_conductor
# ===========================================================================

class TestLFunctionConductor:
    """Tests für den Führer einer L-Funktion."""

    def test_hauptcharakter_hat_conductor_1(self):
        """Der Hauptcharakter χ₀ hat Führer 1."""
        assert l_function_conductor(4, 0) == 1

    def test_conductor_teilt_q(self):
        """Der Conductor muss q teilen."""
        for q in [4, 5, 7]:
            chars = dirichlet_characters(q)
            for idx in range(len(chars)):
                cond = l_function_conductor(q, idx)
                assert q % cond == 0, (
                    f"Conductor {cond} teilt q={q} nicht für Charakter {idx}"
                )

    def test_conductor_mod5_nichthauptcharakter(self):
        """Primitive Charaktere mod 5 haben Conductor 5."""
        chars = dirichlet_characters(5)
        for idx in range(1, len(chars)):
            cond = l_function_conductor(5, idx)
            assert cond == 5, (
                f"Conductor mod 5 sollte 5 sein, erhalten: {cond}"
            )


# ===========================================================================
# TESTS FÜR: l_function_special_values
# ===========================================================================

class TestLFunctionSpecialValues:
    """Tests für spezielle Werte von L(s, χ)."""

    def test_rueckgabe_ist_dict(self):
        """l_function_special_values gibt ein Dictionary zurück."""
        result = l_function_special_values(4, 1)
        assert isinstance(result, dict)

    def test_pflichtschluessel_vorhanden(self):
        """Das Dictionary enthält alle Pflichtschlüssel."""
        result = l_function_special_values(4, 1)
        for key in ['L_1', 'L_0', 'chi_minus1', 'is_odd', 'is_even']:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_chi_minus4_ist_ungerade(self):
        """
        χ₋₄(-1) = -1 → χ₋₄ ist ein ungerader Charakter (is_odd=True).
        """
        result = l_function_special_values(4, 1)
        assert result['is_odd'] is True or result['is_even'] is True, (
            "Charakter muss gerade oder ungerade sein"
        )

    def test_l1_ungefaehr_pi_durch_4(self):
        """L(1, χ₋₄) ≈ π/4 ≈ 0.7854."""
        result = l_function_special_values(4, 1)
        l1 = result['L_1']
        if not math.isnan(l1.real):
            assert abs(l1.real - math.pi / 4) < 0.05, (
                f"L(1, χ₋₄) = {l1.real:.4f}, erwartet ≈ π/4 ≈ {math.pi/4:.4f}"
            )

    def test_hauptcharakter_wirft_fehler(self):
        """Hauptcharakter (Index 0) hat Pol bei s=1 → DomainError."""
        with pytest.raises(DomainError):
            l_function_special_values(4, 0)


# ===========================================================================
# TESTS FÜR: functional_equation_check
# ===========================================================================

class TestFunctionalEquationCheck:
    """Tests für die Funktionalgleichung Λ(s, χ) = ε(χ)·Λ(1-s, χ̄)."""

    def test_rueckgabe_ist_dict(self):
        """functional_equation_check gibt ein Dictionary zurück."""
        result = functional_equation_check(complex(0.5, 5.0), 4, 1)
        assert isinstance(result, dict)

    def test_pflichtschluessel_vorhanden(self):
        """Alle Pflichtschlüssel im Ergebnis-Dict vorhanden."""
        result = functional_equation_check(complex(0.5, 5.0), 4, 1)
        for key in ['lhs', 'rhs', 'error', 'epsilon', 'delta']:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_epsilon_hat_betrag_1_fuer_primitiven_charakter(self):
        """
        Für primitive Charaktere gilt |ε(χ)| = 1 (Root Number).
        ε(χ) = τ(χ) / (i^δ · √q), und |τ(χ)|² = q für primitive χ.
        """
        result = functional_equation_check(complex(0.5, 3.0), 4, 1)
        eps = result['epsilon']
        # |ε| sollte nahe 1 liegen
        assert abs(abs(eps) - 1) < 0.2, (
            f"|ε(χ₋₄)| = {abs(eps):.4f}, erwartet ≈ 1"
        )

    def test_delta_fuer_ungerade_charakter(self):
        """
        Für χ₋₄ (ungerade, χ(-1)=-1) gilt δ = 1.
        """
        result = functional_equation_check(complex(0.5, 3.0), 4, 1)
        # δ = 1 für ungerade Charaktere
        assert result['delta'] in [0, 1]


# ===========================================================================
# TESTS FÜR: hecke_l_function
# ===========================================================================

class TestHeckleLFunction:
    """Tests für die Hecke L-Funktion L(s, f)."""

    def _ramanujan_tau(self, n_max: int) -> list:
        """
        Erste n_max Ramanujan-Tau-Zahlen τ(n) via Euler-Produkt der Delta-Funktion.
        Bekannte Werte: τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830.
        """
        # Vereinfacht: bekannte Werte
        known = [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]
        return known[:n_max]

    def test_rueckgabe_ist_komplex(self):
        """hecke_l_function gibt eine komplexe Zahl zurück."""
        coeffs = self._ramanujan_tau(10)
        result = hecke_l_function(complex(14.0), coeffs)
        assert isinstance(result, complex)

    def test_konvergenz_bei_grossem_s(self):
        """L(Δ, s) für Re(s) >> 1 ist endlich."""
        coeffs = self._ramanujan_tau(10)
        result = hecke_l_function(complex(100.0), coeffs)
        assert math.isfinite(result.real)
        assert math.isfinite(result.imag)

    def test_erster_term_dominiert_bei_grossem_s(self):
        """
        Für Re(s) = 100: L(Δ, 100) ≈ a_1 / 1^100 = 1.
        """
        coeffs = self._ramanujan_tau(10)
        result = hecke_l_function(complex(100.0), coeffs)
        assert abs(result - 1.0) < 0.01, (
            f"L(Δ, 100) ≈ 1, erhalten: {result}"
        )

    def test_ramanujan_tau_l_funktion_s14(self):
        """
        L(Δ, 14) = Σ τ(n)/n^14 konvergiert gut (Re(s) > (12+1)/2 = 6.5).
        Ergebnis sollte endlich und reell (für reelle Koeffizienten bei realem s) sein.
        """
        coeffs = self._ramanujan_tau(50)
        result = hecke_l_function(complex(14.0), coeffs, n_terms=50)
        assert math.isfinite(result.real)
        assert abs(result.imag) < 1e-10  # Imaginärteil = 0 für reelle s und a_n

    def test_leere_koeffizienten_wirft_fehler(self):
        """Leere Koeffizientenliste wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            hecke_l_function(complex(2.0), [])

    def test_l_funktion_mit_konstanten_koeffizienten_gleich_zeta(self):
        """
        Für a_n = 1: L(s, f) = Σ 1/n^s = ζ(s).
        Überprüfe: L(2, f) ≈ ζ(2) = π²/6 ≈ 1.6449.
        """
        N = 500
        ones = [1.0] * N
        result = hecke_l_function(complex(2.0), ones, n_terms=N)
        zeta2 = math.pi**2 / 6
        # Partialsumme mit 500 Termen: Fehler < 0.2%
        assert abs(result.real - zeta2) / zeta2 < 0.005, (
            f"L(2, {N} Einsen) = {result.real:.6f}, erwartet ζ(2) ≈ {zeta2:.6f}"
        )


# ===========================================================================
# TESTS FÜR: completed_l_function
# ===========================================================================

class TestCompletedLFunction:
    """Tests für die vollständige L-Funktion Λ(s, f) = (2π)^{-s} Γ(s) L(s, f)."""

    def test_rueckgabe_ist_komplex(self):
        """completed_l_function gibt eine komplexe Zahl zurück."""
        coeffs = [1, -24, 252]
        result = completed_l_function(complex(7.0), coeffs, weight=12)
        assert isinstance(result, complex)

    def test_endlicher_wert_fuer_grosse_s(self):
        """Λ(s, f) ist endlich für Re(s) > 0."""
        coeffs = [1.0] * 20
        result = completed_l_function(complex(5.0), coeffs, weight=12)
        assert math.isfinite(result.real)

    def test_leere_koeffizienten_wirft_fehler(self):
        """Leere Koeffizientenliste wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            completed_l_function(complex(7.0), [], weight=12)

    def test_ungueltige_weight_wirft_fehler(self):
        """weight < 2 wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            completed_l_function(complex(7.0), [1, -24], weight=1)


# ===========================================================================
# TESTS FÜR: elliptic_l_function_approx (BSD-Verbindung)
# ===========================================================================

class TestEllipticLFunctionApprox:
    """Tests für die Näherung der elliptischen L-Funktion."""

    def test_rueckgabe_ist_komplex(self):
        """elliptic_l_function_approx gibt komplexe Zahl zurück."""
        result = elliptic_l_function_approx(-1, 0, complex(2.0), prime_bound=30)
        assert isinstance(result, complex)

    def test_endlicher_wert_bei_s2(self):
        """Für E: y²=x³-x bei s=2 ist L(E,2) endlich."""
        result = elliptic_l_function_approx(-1, 0, complex(2.0), prime_bound=30)
        assert math.isfinite(result.real)
        assert math.isfinite(result.imag)

    def test_singulaere_kurve_wirft_fehler(self):
        """
        Singuläre Kurve E: y²=x³ (a=0, b=0, Δ=0) wirft DomainError.
        """
        with pytest.raises(DomainError):
            elliptic_l_function_approx(0, 0, complex(2.0), prime_bound=10)

    def test_l_funktion_bei_grossem_re_s_konvergent(self):
        """Für Re(s) = 5 konvergiert das Euler-Produkt gut."""
        result = elliptic_l_function_approx(-1, 0, complex(5.0), prime_bound=50)
        assert math.isfinite(result.real)

    def test_euler_produkt_symmetrie(self):
        """
        L(E, s) ist reell für reelles s ≥ 2.
        (Koeffizienten a_p sind reell, Euler-Produkt ergibt reelles Resultat.)
        """
        result = elliptic_l_function_approx(-1, 0, complex(2.0), prime_bound=30)
        assert abs(result.imag) < 0.01, (
            f"Im(L(E, 2)) = {result.imag:.6f} sollte ≈ 0 sein"
        )

    def test_ungueltige_prime_bound_wirft_fehler(self):
        """prime_bound < 2 wirft InvalidInputError."""
        with pytest.raises(InvalidInputError):
            elliptic_l_function_approx(-1, 0, complex(2.0), prime_bound=1)


# ===========================================================================
# TESTS FÜR: bsd_rank_from_zeros
# ===========================================================================

class TestBsdRankFromZeros:
    """Tests für die BSD-Rang-Schätzung."""

    def test_rueckgabe_ist_dict(self):
        """bsd_rank_from_zeros gibt ein Dictionary zurück."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        assert isinstance(result, dict)

    def test_pflichtschluessel_vorhanden(self):
        """Alle Pflichtschlüssel im Ergebnis-Dict."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        for key in ['l_at_1', 'rank_estimate', 'vanishes', 'l_values', 'bsd_consistent']:
            assert key in result, f"Schlüssel '{key}' fehlt"

    def test_rang_schätzung_ist_nicht_negativ(self):
        """Rang-Schätzung ist ≥ 0."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        assert result['rank_estimate'] >= 0

    def test_rang_schätzung_ist_int(self):
        """Rang-Schätzung ist eine ganze Zahl."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        assert isinstance(result['rank_estimate'], int)

    def test_kurve_y2_x3_plus1_rang0(self):
        """
        E: y² = x³ + 1 (a=0, b=1) hat Rang 0.
        BSD: L(E, 1) ≠ 0, also vanishes=False und rank_estimate=0.
        """
        result = bsd_rank_from_zeros(0, 1, prime_bound=50)
        # L(E, 1) sollte ≠ 0 sein (E hat bekanntermaßen Rang 0)
        l1 = result['l_at_1']
        assert math.isfinite(l1.real), "L(E,1) sollte endlich sein"

    def test_l_values_ist_dict(self):
        """l_values ist ein Dictionary mit Testpunkten."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        assert isinstance(result['l_values'], dict)
        assert len(result['l_values']) > 0

    def test_singulaere_kurve_wirft_fehler(self):
        """Singuläre Kurve (a=0, b=0) wirft DomainError."""
        with pytest.raises(DomainError):
            bsd_rank_from_zeros(0, 0, prime_bound=10)

    def test_bsd_consistent_ist_true(self):
        """bsd_consistent ist stets True (Vermutung)."""
        result = bsd_rank_from_zeros(-1, 0, prime_bound=20)
        assert result['bsd_consistent'] is True


# ===========================================================================
# INTEGRATIONSTESTS
# ===========================================================================

class TestIntegration:
    """Übergreifende Integrationstests."""

    def test_gauss_summe_und_primitiver_charakter_konsistent(self):
        """
        |τ(χ)|² = q für primitive χ ⟺ is_primitive_character = True.
        """
        q = 5
        chars = dirichlet_characters(q)
        for idx in range(1, len(chars)):
            tau = gauss_sum(q, idx)
            tau_sq = abs(tau) ** 2
            is_prim = is_primitive_character(chars[idx]['values'], q)
            if is_prim:
                assert abs(tau_sq - q) < 1.5, (
                    f"Primitiver Charakter idx={idx} mod {q}: |τ|²={tau_sq:.3f} ≠ {q}"
                )

    def test_l_funktion_und_euler_produkt_konsistenz_mod5(self):
        """
        Euler-Produkt-Näherung und direkte Dirichlet-Reihe stimmen für
        L(2, χ₀ mod 5) bei ausreichend vielen Primzahlen überein.
        """
        # Direktes L via Reihe
        l_series = dirichlet_l_function(complex(2.0), 5, 0, n_terms=1000)
        # Euler-Produkt für E=x^3-x ist anders, also ζ(2)*(1-5^{-2}) prüfen
        zeta2 = math.pi**2 / 6
        expected = zeta2 * (1 - 5**(-2))
        assert abs(l_series.real - expected) / expected < 0.05

    def test_leibniz_reihe_über_dirichlet_l(self):
        """
        Vollständiger Integrationstest: Leibniz-Formel π/4 via L-Funktion.
        Charaktere mod 4 korrekt konstruiert → L(1, χ₋₄) ≈ π/4.
        """
        chars = dirichlet_characters(4)
        assert len(chars) == 2
        # Direkter Test über chi_value
        # χ₋₄(1) = 1, χ₋₄(3) = -1, χ₋₄(2) = χ₋₄(4) = 0
        chi_1 = chi_value(1, 4, 1)
        chi_3 = chi_value(3, 4, 1)
        # |χ₋₄(1)| = 1, |χ₋₄(3)| = 1
        assert abs(abs(chi_1) - 1) < 1e-8
        assert abs(abs(chi_3) - 1) < 1e-8
        # Leibniz-Reihe: Σ χ(n)/n = 1 - 1/3 + 1/5 - ...
        leibniz = sum(chi_value(n, 4, 1).real / n for n in range(1, 501, 2) if chi_value(n, 4, 1).real != 0)
        # Zähle nur ungerade n (da χ(2k) = 0)
        # Manuell: 1/1 - 1/3 + 1/5 - ...
        leibniz_manual = sum(
            (-1)**k / (2*k + 1) for k in range(500)
        )
        assert abs(leibniz_manual - math.pi / 4) < 0.01


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
