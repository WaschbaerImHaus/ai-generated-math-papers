"""
@file test_langlands_program.py
@brief Umfassende pytest-Tests für das Modul langlands_program.py.
@description
    Testet alle 8 Klassen und 3 Hilfsfunktionen des Langlands-Programm-Moduls:
    - SatakeIsomorphism: Satake-Parameter, sphärische Funktion, Hecke-Algebra
    - LocalLanglandsGL1: Klassenkörpertheorie, Artin-Reziprozität
    - LocalLanglandsGL2: Hauptreihen, spezielle, supercuspidale Darstellungen
    - GlobalLanglands: Modularitätssatz, Lifting-Theorem
    - LanglandsFunctoriality: sym^k-Lifts, Basiswechsel, automorphe Induktion
    - ArtinReciprocity: Artin-L-Funktion, Artin-Führer
    - TraceFormula: geometrische/spektrale Seite der Spurformel
    - LDualGroup: Dualgruppen, L-Gruppen, Funktorialitäts-Checks
    - Hilfsfunktionen: local_langlands_parameter, epsilon_factor_normalization,
      langlands_correspondence_table

    Mathematische Grundlagen:
    - Satake-Isomorphismus: H(G//K) ≅ ℂ[X*(T̂)^W]
    - Lokale Langlands-Korrespondenz: irred. π ↔ WD-Darstellungen
    - Globale Korrespondenz (GL_2/Q): Modulform ↔ Galois-Darstellung
    - Artin-Reziprozität: Art_p: Q_p^× → Gal(Q_p^ab/Q_p)

    Alle Tests sind unabhängig voneinander und laufen parallelisierbar
    unter pytest-xdist (keine geteilten Zustände).

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
import sys
import os
import pytest

# Pfad zu den Quelldateien
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from langlands_program import (
    # Hilfsfunktionen
    local_langlands_parameter,
    epsilon_factor_normalization,
    langlands_correspondence_table,
    # Klassen
    SatakeIsomorphism,
    LocalLanglandsGL1,
    LocalLanglandsGL2,
    GlobalLanglands,
    LanglandsFunctoriality,
    ArtinReciprocity,
    TraceFormula,
    LDualGroup,
)


# ===========================================================================
# HILFSFUNKTION: local_langlands_parameter
# ===========================================================================

class TestLocalLanglandsParameter:
    """Testet die Hilfsfunktion local_langlands_parameter."""

    def test_rückgabe_dict_schlüssel(self):
        """Ergebnis-Dict enthält alle erwarteten Schlüssel."""
        result = local_langlands_parameter(5, 2.0)
        expected_keys = {"p", "a_p", "alpha", "beta", "|alpha|", "discriminant", "type",
                         "ramanujan_satisfied"}
        assert expected_keys == set(result.keys())

    def test_elliptischer_fall_typ(self):
        """Für a_p² < 4p ergibt sich der elliptische Typ."""
        # p=5, a_p=2: disc = 4 - 20 = -16 < 0 → elliptisch
        result = local_langlands_parameter(5, 2.0)
        assert result["type"] == "elliptic"

    def test_elliptischer_fall_alpha_beta_konjugiert(self):
        """Im elliptischen Fall sind α und β komplex-konjugiert."""
        result = local_langlands_parameter(5, 2.0)
        alpha = result["alpha"]
        beta = result["beta"]
        assert abs(alpha.real - beta.real) < 1e-12
        assert abs(alpha.imag + beta.imag) < 1e-12

    def test_elliptischer_fall_alpha_times_beta_gleich_p(self):
        """Im elliptischen Fall gilt α·β = p (unitäre Normierung Gewicht 2)."""
        p = 7
        a_p = 0.0
        result = local_langlands_parameter(p, a_p)
        alpha = result["alpha"]
        beta = result["beta"]
        product = alpha * beta
        assert abs(product.real - p) < 1e-8
        assert abs(product.imag) < 1e-8

    def test_elliptischer_fall_ramanujan_erfüllt(self):
        """Ramanujan-Vermutung: |α_p| = √p für elliptischen Fall."""
        p = 5
        a_p = 2.0
        result = local_langlands_parameter(p, a_p)
        assert abs(result["|alpha|"] - math.sqrt(p)) < 1e-8
        assert result["ramanujan_satisfied"] is True

    def test_parabolischer_fall_typ(self):
        """Für a_p² = 4p ergibt sich der parabolische Typ."""
        p = 4
        a_p = 4.0  # disc = 16 - 16 = 0
        result = local_langlands_parameter(p, a_p)
        assert result["type"] == "parabolic"

    def test_parabolischer_fall_alpha_gleich_beta(self):
        """Im parabolischen Fall gilt α = β."""
        p = 4
        a_p = 4.0
        result = local_langlands_parameter(p, a_p)
        assert abs(result["alpha"] - result["beta"]) < 1e-12

    def test_hyperbolischer_fall_typ(self):
        """Für a_p² > 4p ergibt sich der hyperbolische Typ."""
        p = 2
        a_p = 4.0  # disc = 16 - 8 = 8 > 0
        result = local_langlands_parameter(p, a_p)
        assert result["type"] == "hyperbolic"

    def test_alpha_plus_beta_gleich_a_p(self):
        """α + β = a_p gilt in allen Fällen (Vieta)."""
        p = 5
        a_p = 3.0
        result = local_langlands_parameter(p, a_p)
        summe = result["alpha"] + result["beta"]
        assert abs(summe.real - a_p) < 1e-8
        assert abs(summe.imag) < 1e-8

    def test_diskriminante_korrekt(self):
        """Diskriminante = a_p² - 4p wird korrekt berechnet."""
        p = 3
        a_p = 1.0
        result = local_langlands_parameter(p, a_p)
        expected_disc = 1.0 - 12.0
        assert abs(result["discriminant"] - expected_disc) < 1e-12


# ===========================================================================
# HILFSFUNKTION: epsilon_factor_normalization
# ===========================================================================

class TestEpsilonFactorNormalization:
    """Testet die Normierung des ε-Faktors."""

    def test_bei_s_halb_ergibt_eins(self):
        """Bei s=1/2 gilt N^{1/2-s} = N^0 = 1, also ε = root_number = 1."""
        result = epsilon_factor_normalization(27, 0.5)
        assert abs(result - 1.0) < 1e-10

    def test_bei_s_null_ergibt_sqrt_n(self):
        """Bei s=0 gilt ε = N^{1/2} = √N (root_number=1)."""
        N = 16
        result = epsilon_factor_normalization(N, 0.0)
        expected = math.sqrt(N)
        assert abs(result.real - expected) < 1e-8

    def test_bei_s_eins_ergibt_eins_durch_sqrt_n(self):
        """Bei s=1 gilt ε = N^{-1/2} = 1/√N."""
        N = 25
        result = epsilon_factor_normalization(N, 1.0)
        expected = 1.0 / math.sqrt(N)
        assert abs(result.real - expected) < 1e-8

    def test_leiter_eins_ergibt_immer_eins(self):
        """Bei Führer N=1 ist ε = 1 für alle s."""
        for s in [0.0, 0.5, 1.0, 2.0]:
            result = epsilon_factor_normalization(1, s)
            assert abs(result - 1.0) < 1e-10

    def test_rückgabe_ist_komplex(self):
        """Der ε-Faktor hat den Typ komplex."""
        result = epsilon_factor_normalization(11, 0.5 + 1j)
        assert isinstance(result, complex)


# ===========================================================================
# HILFSFUNKTION: langlands_correspondence_table
# ===========================================================================

class TestLanglandsCorrespondenceTable:
    """Testet die Tabellenfunktion für Langlands-Parameter."""

    def test_leere_liste(self):
        """Leere Eingabe ergibt leere Tabelle."""
        result = langlands_correspondence_table([])
        assert result == []

    def test_rückgabe_liste_von_dicts(self):
        """Ergebnis ist eine Liste von Dicts."""
        result = langlands_correspondence_table([(5, 2)])
        assert isinstance(result, list)
        assert isinstance(result[0], dict)

    def test_tabellen_schlüssel(self):
        """Jeder Tabelleneintrag enthält alle erforderlichen Schlüssel."""
        result = langlands_correspondence_table([(5, 2)])
        entry = result[0]
        for key in ("p", "a_p", "alpha_p", "beta_p", "|alpha_p|", "sqrt_p", "ramanujan", "type"):
            assert key in entry, f"Schlüssel '{key}' fehlt"

    def test_mehrere_primzahlen(self):
        """Tabelle kann mehrere Einträge enthalten."""
        # E: y² = x³ - x, Leiter N=32
        ap_data = [(5, 2), (13, 6), (17, 2), (29, -2)]
        result = langlands_correspondence_table(ap_data)
        assert len(result) == 4

    def test_sqrt_p_korrekt(self):
        """Das Feld sqrt_p entspricht √p."""
        result = langlands_correspondence_table([(5, 2)])
        assert abs(result[0]["sqrt_p"] - math.sqrt(5)) < 1e-12

    def test_a_p_3_mod_4_Null(self):
        """Für p≡3(4) und E: y²=x³-x gilt a_p=0 (Frobenius=komp. Konjugation)."""
        # p=3 ≡ 3 (mod 4), a_p=0
        result = langlands_correspondence_table([(3, 0)])
        entry = result[0]
        alpha = entry["alpha_p"]
        beta = entry["beta_p"]
        # α + β = 0
        assert abs(alpha + beta) < 1e-10

    def test_ramanujan_erfüllt_für_gute_reduktion(self):
        """Ramanujan ist True, wenn |α_p| = √p."""
        result = langlands_correspondence_table([(5, 2)])
        assert result[0]["ramanujan"] is True


# ===========================================================================
# KLASSE: SatakeIsomorphism
# ===========================================================================

class TestSatakeIsomorphism:
    """Testet den Satake-Isomorphismus für GL_n."""

    def test_initialisierung_n2(self):
        """Initialisierung mit n=2 setzt Attribut korrekt."""
        sat = SatakeIsomorphism(n=2)
        assert sat.n == 2

    def test_satake_parameter_aus_hecke_elliptisch(self):
        """Satake-Parameter (α,β) für a_p mit |a_p| < 2√p·p^{(k-1)/2}."""
        sat = SatakeIsomorphism(n=2)
        # Für Gewicht 2, p=5, a_p=2: norm = √5
        alpha, beta = sat.satake_parameters_from_hecke(a_p=2.0, p=5, weight=2)
        # α·β = p^{(k-1)} = p für k=2
        assert abs(alpha * beta - 5.0) < 1e-8

    def test_satake_parameter_summe_gleich_a_p(self):
        """α + β = a_p (normiert)."""
        sat = SatakeIsomorphism(n=2)
        a_p = 3.0
        alpha, beta = sat.satake_parameters_from_hecke(a_p=a_p, p=7, weight=2)
        # Normiert: a_p = 2·norm·cos(θ), aber α+β = a_p (direkte Prüfung)
        assert abs((alpha + beta).real - a_p) < 1e-8

    def test_sphärische_funktion_k0_ergibt_eins(self):
        """φ_{α,β}(p^0) = 1 (Normierung der sphärischen Funktion)."""
        sat = SatakeIsomorphism(n=2)
        # k=0: Summe Σ_{j=0}^{0} α^0·β^0 = 1
        alpha = complex(math.sqrt(5), 0)
        beta = complex(math.sqrt(5), 0)
        val = sat.spherical_function(alpha, beta, p=5, k=0)
        assert abs(val - 1.0) < 1e-10

    def test_sphärische_funktion_k1_gleich_a_p(self):
        """φ_{α,β}(p^1) = α + β = a_p (Hecke-Eigenwert)."""
        sat = SatakeIsomorphism(n=2)
        alpha = complex(2.0, 1.0)
        beta = complex(2.0, -1.0)
        val = sat.spherical_function(alpha, beta, p=5, k=1)
        expected = alpha + beta
        assert abs(val - expected) < 1e-10

    def test_sphärische_funktion_alpha_gleich_beta(self):
        """Entarteter Fall α=β: φ = (k+1)·α^k."""
        sat = SatakeIsomorphism(n=2)
        alpha = complex(2.0, 0.0)
        beta = complex(2.0, 0.0)
        k = 3
        val = sat.spherical_function(alpha, beta, p=5, k=k)
        expected = (k + 1) * (alpha ** k)
        assert abs(val - expected) < 1e-10

    def test_hecke_algebra_generators_struktur(self):
        """hecke_algebra_generators gibt Dict mit Generatoren zurück."""
        sat = SatakeIsomorphism(n=2)
        result = sat.hecke_algebra_generators(p=5)
        assert "group" in result
        assert "generators" in result
        assert result["group"] == "GL_2"
        assert len(result["generators"]) == 2  # T_{p^1}, T_{p^2}

    def test_hecke_algebra_weyl_gruppe(self):
        """Die Weyl-Gruppe von GL_n ist S_n."""
        sat = SatakeIsomorphism(n=3)
        result = sat.hecke_algebra_generators(p=7)
        assert result["weyl_group"] == "S_3"

    def test_satake_transform_linearität(self):
        """Satake-Transform: f → S(f) ist linear in den Werten."""
        sat = SatakeIsomorphism(n=2)
        alpha = complex(math.sqrt(5), 0)
        beta = complex(math.sqrt(5), 0)
        # f1 = {0: 1.0}, f2 = {0: 2.0}: S(f2) = 2·S(f1)
        f1 = {0: 1.0 + 0j}
        f2 = {0: 2.0 + 0j}
        t1 = sat.satake_transform(f1, alpha, beta)
        t2 = sat.satake_transform(f2, alpha, beta)
        assert abs(t2 - 2 * t1) < 1e-10


# ===========================================================================
# KLASSE: LocalLanglandsGL1
# ===========================================================================

class TestLocalLanglandsGL1:
    """Testet die lokale Langlands-Korrespondenz für GL_1."""

    def test_initialisierung_primzahl(self):
        """Initialisierung setzt prime korrekt."""
        llgl1 = LocalLanglandsGL1(p=5)
        assert llgl1.p == 5

    def test_class_field_theory_map_unramifiziert(self):
        """Für p ∤ mod ist die Abbildung unramifiziert."""
        llgl1 = LocalLanglandsGL1(p=5)
        # Trivial-Charakter: χ(n)=1 für alle n
        chi_values = {5: 1.0 + 0j, 7: 1.0 + 0j}
        result = llgl1.class_field_theory_map(chi_values, mod=7)
        assert result["is_ramified"] is False

    def test_class_field_theory_map_ramifiziert(self):
        """Für p | mod ist die Abbildung ramifiziert."""
        llgl1 = LocalLanglandsGL1(p=5)
        chi_values = {5: 1.0 + 0j}
        result = llgl1.class_field_theory_map(chi_values, mod=25)
        assert result["is_ramified"] is True

    def test_class_field_theory_map_schlüssel(self):
        """Rückgabe-Dict enthält alle erwarteten Schlüssel."""
        llgl1 = LocalLanglandsGL1(p=3)
        chi_values = {3: 1.0, 7: 1.0}
        result = llgl1.class_field_theory_map(chi_values, mod=7)
        for key in ("prime", "modulus", "is_ramified", "conductor_exponent",
                    "chi_p_frobenius", "rho_frobenius_eigenvalue", "artin_conductor"):
            assert key in result

    def test_frobenius_eigenwert_für_trivialcharakter(self):
        """Für den Trivial-Charakter ist der Frobenius-Eigenwert 1."""
        llgl1 = LocalLanglandsGL1(p=7)
        chi_values = {7: 1.0 + 0j}
        result = llgl1.class_field_theory_map(chi_values, mod=11)
        assert abs(result["rho_frobenius_eigenvalue"] - 1.0) < 1e-12

    def test_artin_reciprocity_local_a1(self):
        """Für a=1 ist die Valuation v_p(1)=0 und Artin-Symbol = 1."""
        llgl1 = LocalLanglandsGL1(p=3)
        result = llgl1.artin_reciprocity_local(a=1, n=5)
        assert result["v_p(a)"] == 0
        assert result["artin_symbol_mod_n"] == 1

    def test_artin_reciprocity_local_a_p(self):
        """Für a=p ist die Valuation v_p(p)=1."""
        p = 5
        llgl1 = LocalLanglandsGL1(p=p)
        result = llgl1.artin_reciprocity_local(a=p, n=7)
        assert result["v_p(a)"] == 1

    def test_artin_reciprocity_local_fehler_bei_null(self):
        """ValueError wenn a=0 übergeben wird."""
        llgl1 = LocalLanglandsGL1(p=3)
        with pytest.raises(ValueError):
            llgl1.artin_reciprocity_local(a=0, n=5)

    def test_local_norm_residue_symbol_p_ungerade_einheit(self):
        """Normrestsymbol für p ungerade und a = Einheit."""
        llgl1 = LocalLanglandsGL1(p=5)
        # a=4 = 2², Einheit mod 5, n=2 → Legendre(4, 5) = 1 (4≡4, 2²≡4 mod5 → QR)
        symbol = llgl1.local_norm_residue_symbol(a=4, p=5, n=2)
        assert symbol in {-1, 0, 1}

    def test_local_norm_residue_symbol_rückgabe_int(self):
        """Das Normrestsymbol gibt einen Integer zurück."""
        llgl1 = LocalLanglandsGL1(p=7)
        result = llgl1.local_norm_residue_symbol(a=9, p=7, n=2)
        assert isinstance(result, int)


# ===========================================================================
# KLASSE: LocalLanglandsGL2
# ===========================================================================

class TestLocalLanglandsGL2:
    """Testet die lokale Langlands-Korrespondenz für GL_2."""

    def test_initialisierung(self):
        """Initialisierung setzt prime korrekt."""
        llgl2 = LocalLanglandsGL2(p=7)
        assert llgl2.p == 7

    def test_principal_series_schlüssel(self):
        """principal_series gibt Dict mit allen Schlüsseln zurück."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.principal_series(chi1_p=1.0+0j, chi2_p=1.0+0j)
        for key in ("representation", "chi1_p", "chi2_p", "hecke_eigenvalue_a_p",
                    "central_character_omega_p", "is_irreducible", "satake_parameters"):
            assert key in result

    def test_principal_series_rep_typ(self):
        """Darstellungstyp ist 'principal_series'."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.principal_series(chi1_p=1.0, chi2_p=2.0)
        assert result["representation"] == "principal_series"

    def test_principal_series_hecke_eigenwert(self):
        """Hecke-Eigenwert a_p = χ_1(p) + χ_2(p)."""
        llgl2 = LocalLanglandsGL2(p=5)
        chi1 = 2.0 + 0j
        chi2 = 3.0 + 0j
        result = llgl2.principal_series(chi1_p=chi1, chi2_p=chi2)
        assert abs(result["hecke_eigenvalue_a_p"] - (chi1 + chi2)) < 1e-12

    def test_principal_series_zentralcharakter(self):
        """Zentralcharakter ω(p) = χ_1(p)·χ_2(p)."""
        llgl2 = LocalLanglandsGL2(p=5)
        chi1 = 2.0 + 0j
        chi2 = 3.0 + 0j
        result = llgl2.principal_series(chi1_p=chi1, chi2_p=chi2)
        assert abs(result["central_character_omega_p"] - chi1 * chi2) < 1e-12

    def test_special_representation_schlüssel(self):
        """special_representation gibt Dict mit allen Schlüsseln zurück."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.special_representation(chi_p=1.0+0j)
        for key in ("representation", "chi_p", "weil_deligne_eigenvalues",
                    "monodromy_N", "conductor_exponent"):
            assert key in result

    def test_special_representation_steinberg_typ(self):
        """Darstellungstyp ist 'special'."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.special_representation(chi_p=1.0)
        assert result["representation"] == "special"

    def test_special_representation_führer_exponent(self):
        """Steinberg-Darstellung hat Führer-Exponent 1."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.special_representation(chi_p=1.0)
        assert result["conductor_exponent"] == 1

    def test_special_representation_weil_deligne_eigenwerte(self):
        """WD-Eigenwerte sind χ(p)·√p und χ(p)/√p."""
        p = 5
        llgl2 = LocalLanglandsGL2(p=p)
        chi_p = 1.0
        result = llgl2.special_representation(chi_p=chi_p)
        rho1, rho2 = result["weil_deligne_eigenvalues"]
        assert abs(rho1 - chi_p * math.sqrt(p)) < 1e-10
        assert abs(rho2 - chi_p / math.sqrt(p)) < 1e-10

    def test_supercuspidal_marker_typ(self):
        """supercuspidal_marker gibt Typ 'supercuspidal' zurück."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.supercuspidal_marker()
        assert result["representation"] == "supercuspidal"

    def test_supercuspidal_l_faktor_eins(self):
        """L-Faktor für supercuspidale Darstellung ist '1'."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.supercuspidal_marker()
        assert result["l_factor"] == "1"

    def test_local_factor_l_principal_series_bei_großem_s(self):
        """L(s, π) für Hauptreihe mit α=β=1 bei großem Re(s) nahe 1."""
        llgl2 = LocalLanglandsGL2(p=5)
        # α=β=p^{-1/2}: unitäre Normierung (Ramanujan)
        alpha = complex(1.0 / math.sqrt(5), 0)
        beta = complex(1.0 / math.sqrt(5), 0)
        s = 2.0
        val = llgl2.local_factor_l("principal_series", s, alpha, beta)
        # Ergebnis ist endlich (kein Pol für Re(s) > 1)
        assert math.isfinite(abs(val))
        assert abs(val) > 0

    def test_local_factor_l_special_bei_alpha_null(self):
        """L(s, Sp) = 1 wenn α = 0 (Pol vermieden)."""
        llgl2 = LocalLanglandsGL2(p=5)
        val = llgl2.local_factor_l("special", s=2.0, alpha=0.0+0j)
        assert abs(val - 1.0) < 1e-12

    def test_local_factor_l_supercuspidal_immer_eins(self):
        """L(s, π_sc) = 1 für alle s."""
        llgl2 = LocalLanglandsGL2(p=5)
        for s in [0.5, 1.0, 2.0]:
            val = llgl2.local_factor_l("supercuspidal", s=s)
            assert abs(val - 1.0) < 1e-12

    def test_local_factor_l_unbekannter_typ_fehler(self):
        """ValueError bei unbekanntem Darstellungstyp."""
        llgl2 = LocalLanglandsGL2(p=5)
        with pytest.raises(ValueError):
            llgl2.local_factor_l("unbekannt", s=1.0)

    def test_local_epsilon_factor_principal_series_eins(self):
        """ε(1/2, π) = 1 für unramifizierte Hauptreihe."""
        llgl2 = LocalLanglandsGL2(p=5)
        eps = llgl2.local_epsilon_factor("principal_series")
        assert abs(eps - 1.0) < 1e-12

    def test_local_epsilon_factor_special_eins(self):
        """ε(1/2, Sp(χ_triv)) = 1 für triviale Steinberg."""
        llgl2 = LocalLanglandsGL2(p=5)
        eps = llgl2.local_epsilon_factor("special")
        assert abs(eps - 1.0) < 1e-12

    def test_local_epsilon_factor_supercuspidal_betrag_eins(self):
        """|ε(1/2, π_sc)| ≈ 1 (unitäre ε-Faktoren)."""
        llgl2 = LocalLanglandsGL2(p=5)
        eps = llgl2.local_epsilon_factor("supercuspidal", psi_level=1)
        # |ε| sollte nahe 1 sein (Gauß-Summe / p^{1/2} für Niveau 1)
        assert abs(abs(eps) - 1.0) < 1e-10

    def test_weil_deligne_to_admissible_principal_series(self):
        """WD → admissible für Hauptreihe gibt korrekten Typ zurück."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.weil_deligne_to_admissible(
            "principal_series", alpha=1.0+0j, beta=1.0+0j
        )
        assert result["type"] == "principal_series"

    def test_weil_deligne_to_admissible_special(self):
        """WD → admissible für Steinberg gibt Führer-Exponent 1."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.weil_deligne_to_admissible(
            "special", alpha=1.0+0j, beta=1.0+0j
        )
        assert result["conductor_exponent"] == 1

    def test_weil_deligne_to_admissible_supercuspidal(self):
        """WD → admissible für supercuspidal gibt Führer ≥ 2."""
        llgl2 = LocalLanglandsGL2(p=5)
        result = llgl2.weil_deligne_to_admissible(
            "supercuspidal", alpha=1.0+0j, beta=1.0+0j
        )
        assert result["conductor_exponent"] >= 2

    def test_weil_deligne_unbekannter_typ_fehler(self):
        """ValueError bei unbekanntem WD-Typ."""
        llgl2 = LocalLanglandsGL2(p=5)
        with pytest.raises(ValueError):
            llgl2.weil_deligne_to_admissible("falsch", alpha=0+0j, beta=0+0j)


# ===========================================================================
# KLASSE: GlobalLanglands
# ===========================================================================

class TestGlobalLanglands:
    """Testet die globale Langlands-Korrespondenz (Modularitätssatz)."""

    def test_modular_curve_wiles_check_rückgabe_schlüssel(self):
        """modular_curve_wiles_check gibt Dict mit Pflichtfeldern zurück."""
        gl = GlobalLanglands()
        a_coeff = {5: 2, 7: 0, 11: -2, 13: 6}
        result = gl.modular_curve_wiles_check(a_coeff, N=32)
        for key in ("conductor", "prime_checks", "all_ramanujan_satisfied",
                    "modularity_theorem", "conclusion"):
            assert key in result

    def test_modular_curve_wiles_check_ramanujan_schranke(self):
        """Ramanujan-Schranke |a_p| ≤ 2√p wird korrekt geprüft."""
        gl = GlobalLanglands()
        # a_5=2: |2| ≤ 2√5 ≈ 4.47 → erfüllt
        a_coeff = {5: 2}
        result = gl.modular_curve_wiles_check(a_coeff, N=32)
        assert result["prime_checks"][5]["satisfies_ramanujan"] is True

    def test_modular_curve_wiles_check_überschreitend_false(self):
        """Wenn |a_p| > 2√p, ist Ramanujan nicht erfüllt."""
        gl = GlobalLanglands()
        a_coeff = {5: 100}  # |100| >> 2√5
        result = gl.modular_curve_wiles_check(a_coeff, N=1)
        assert result["prime_checks"][5]["satisfies_ramanujan"] is False
        assert result["all_ramanujan_satisfied"] is False

    def test_modularity_lifting_theorem_demo_rückgabe(self):
        """modularity_lifting_theorem_demo gibt evidence dict zurück."""
        gl = GlobalLanglands()
        # Elliptische Kurve E: y² = x³ - x, Führer 32
        E_ap = {3: 0, 5: 2, 7: 0, 11: -2, 13: 6}
        result = gl.modularity_lifting_theorem_demo(E_conductor=32, E_ap=E_ap)
        assert isinstance(result, dict)
        assert "theorem" in result
        assert "conclusion" in result

    def test_modularity_lifting_theorem_demo_residual_mod_ell(self):
        """Residuelle Darstellungen mod 3, 5, 7 werden berechnet."""
        gl = GlobalLanglands()
        E_ap = {5: 2, 7: 0}
        result = gl.modularity_lifting_theorem_demo(E_conductor=32, E_ap=E_ap)
        assert "residual_mod_ell" in result
        residues = result["residual_mod_ell"]
        # Für ℓ=3,5,7 müssen Einträge existieren
        for ell in [3, 5, 7]:
            assert ell in residues

    def test_modularity_lifting_theorem_demo_schlussfolgerung_enthält_n(self):
        """Schlussfolgerung enthält den Führer N."""
        gl = GlobalLanglands()
        E_ap = {5: 2}
        result = gl.modularity_lifting_theorem_demo(E_conductor=32, E_ap=E_ap)
        assert "32" in result["conclusion"]

    def test_langlands_conjecture_evidence_gl1_bewiesen(self):
        """Für GL_1 ist die Vermutung bewiesen (Klassenkörpertheorie)."""
        gl = GlobalLanglands()
        result = gl.langlands_conjecture_evidence(rep_dim=1)
        assert result["evidence"]["status"] == "Bewiesen"

    def test_langlands_conjecture_evidence_gl2_bewiesen(self):
        """Für GL_2 ist die Vermutung bewiesen."""
        gl = GlobalLanglands()
        result = gl.langlands_conjecture_evidence(rep_dim=2)
        assert "Bewiesen" in result["evidence"]["status"]

    def test_langlands_conjecture_evidence_groß_n_offen(self):
        """Für n ≥ 5 ist die Vermutung noch offen (für Zahlkörper)."""
        gl = GlobalLanglands()
        result = gl.langlands_conjecture_evidence(rep_dim=10)
        assert "Offen" in result["evidence"]["status"] or "offen" in result["evidence"]["status"]

    def test_langlands_conjecture_evidence_funktionskörper_lafforgue(self):
        """Lafforgue-Referenz erscheint für alle Dimensionen."""
        gl = GlobalLanglands()
        result = gl.langlands_conjecture_evidence(rep_dim=5)
        assert "Lafforgue" in result["function_field"] or "Lafforgue" in result["function_field"]


# ===========================================================================
# KLASSE: LanglandsFunctoriality
# ===========================================================================

class TestLanglandsFunctoriality:
    """Testet die Langlands-Funktorialitäts-Lifts."""

    def test_symmetric_power_lift_k1_identität(self):
        """sym^1-Lift ist die Identität: a_p(Sym^1) = a_p."""
        lf = LanglandsFunctoriality()
        a_p = {5: 2, 7: -2}
        result = lf.symmetric_power_lift(a_p=a_p, p=5, k=1)
        # Für k=1: sym^1-Eigenwert = α+β = a_p
        for prime, ap in a_p.items():
            data = result["results"][prime]
            sym1 = data["sym_k_eigenvalue"]
            assert abs(sym1.real - ap) < 1e-6, (
                f"sym^1 sollte a_p={ap} ergeben, bekam {sym1}"
            )

    def test_symmetric_power_lift_k2_gl3(self):
        """sym²-Lift: Zielgruppe ist GL_3."""
        lf = LanglandsFunctoriality()
        result = lf.symmetric_power_lift(a_p={5: 2}, p=5, k=2)
        assert result["results"][5]["sym_k_group"] == "GL_3"

    def test_symmetric_power_lift_k2_status_bewiesen(self):
        """sym²-Lift (Gelbart-Jacquet) ist bewiesen."""
        lf = LanglandsFunctoriality()
        result = lf.symmetric_power_lift(a_p={5: 2}, p=5, k=2)
        assert "bewiesen" in result["status"]

    def test_symmetric_power_lift_k5_offen(self):
        """sym^5-Lift ist noch offen."""
        lf = LanglandsFunctoriality()
        result = lf.symmetric_power_lift(a_p={5: 2}, p=5, k=5)
        assert "offen" in result["status"]

    def test_symmetric_power_lift_anzahl_satake_parameter(self):
        """sym^k-Lift hat k+1 Satake-Parameter."""
        lf = LanglandsFunctoriality()
        k = 3
        result = lf.symmetric_power_lift(a_p={5: 2}, p=5, k=k)
        params = result["results"][5]["sym_k_satake_params"]
        assert len(params) == k + 1

    def test_base_change_n1_ergibt_a_p(self):
        """Basiswechsel für n=1 (triviale Erweiterung): a_q = a_p."""
        lf = LanglandsFunctoriality()
        a_p = {5: 2, 7: -2}
        result = lf.base_change_GL2(a_p=a_p, p=5, n=1)
        for prime, ap in a_p.items():
            data = result["results"][prime]
            # power_sums[1] = a_p
            assert abs(data["power_sums"][1].real - ap) < 1e-8

    def test_base_change_struktur(self):
        """base_change_GL2 gibt Dict mit theorem, results zurück."""
        lf = LanglandsFunctoriality()
        result = lf.base_change_GL2(a_p={5: 2}, p=5, n=2)
        assert "theorem" in result
        assert "results" in result

    def test_base_change_n2_potenzsummen_korrekt(self):
        """Für n=2: a_q(BC) = α² + β² = a_p² - 2p (Newton-Identität)."""
        lf = LanglandsFunctoriality()
        p = 5
        a_p = 2
        result = lf.base_change_GL2(a_p={p: a_p}, p=p, n=2)
        bc_aq = result["results"][p]["base_change_a_q"]
        # α²+β² = (α+β)² - 2αβ = a_p² - 2p
        expected = a_p ** 2 - 2 * p
        assert abs(bc_aq.real - expected) < 1e-6

    def test_automorphic_induction_n2_a_p_reell(self):
        """Automorphe Induktion für [K:Q]=2: a_p ist reell (2·Re(χ(p)))."""
        lf = LanglandsFunctoriality()
        chi_values = {5: complex(1.0, 0.5)}
        result = lf.automorphic_induction(chi_values=chi_values, p=5, K_degree=2)
        # a_p = 2·Re(χ(p))
        expected_ap = 2 * chi_values[5].real
        assert abs(result["a_p_induced"].real - expected_ap) < 1e-10

    def test_automorphic_induction_schlüssel(self):
        """automorphic_induction gibt Dict mit allen Pflichtfeldern zurück."""
        lf = LanglandsFunctoriality()
        result = lf.automorphic_induction(chi_values={5: 1.0+0j}, p=5, K_degree=2)
        for key in ("lift", "K_degree", "prime_p", "induced_satake_params",
                    "a_p_induced", "status"):
            assert key in result

    def test_langlands_dual_group_gl_n_selbstdual(self):
        """Dualgruppe von GL_n ist GL_n(ℂ) (selbstdual)."""
        lf = LanglandsFunctoriality()
        result = lf.langlands_dual_group("GL_n")
        assert result["is_self_dual"] is True
        assert "GL_n" in result["dual_group"]

    def test_langlands_dual_group_sp2n_ergibt_so(self):
        """Dualgruppe von Sp_{2n} ist SO_{2n+1}(ℂ)."""
        lf = LanglandsFunctoriality()
        result = lf.langlands_dual_group("Sp_2n")
        assert "SO" in result["dual_group"]


# ===========================================================================
# KLASSE: ArtinReciprocity
# ===========================================================================

class TestArtinReciprocity:
    """Testet die Artin-Reziprozität und Artin-L-Funktionen."""

    def test_artin_l_function_trivialcharakter_approximiert_zeta(self):
        """Für den Trivialcharakter approximiert L(s,1) die Riemann-Zeta-Funktion."""
        ar = ArtinReciprocity()
        primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
        chi_trivial = {p: 1.0 for p in primes}
        s = 2.0
        l_val = ar.artin_l_function(chi_trivial, s=s, primes=primes)
        # ζ(2) = π²/6 ≈ 1.6449...
        # Partielles Euler-Produkt konvergiert von unten
        assert l_val.real > 1.0
        # Grobe Schranke: partielles Produkt < ζ(2)+0.5
        assert l_val.real < 2.5

    def test_artin_l_function_gibt_eins_ohne_primes(self):
        """L-Funktion mit leerer Primzahl-Liste ergibt 1."""
        ar = ArtinReciprocity()
        result = ar.artin_l_function({}, s=2.0, primes=[])
        assert abs(result - 1.0) < 1e-12

    def test_artin_conductor_formula_keine_ramifizierten_primes(self):
        """Führer ist 1 wenn keine ramifizierten Primzahlen vorliegen."""
        ar = ArtinReciprocity()
        result = ar.artin_conductor_formula(ramified_primes={}, chi_order=1)
        assert result == 1

    def test_artin_conductor_formula_einfache_ramifizierung(self):
        """Einfache Ramifizierung bei p gibt Führer p^{tame_exp+swan_exp}."""
        ar = ArtinReciprocity()
        # Tame: chi_order=1, inertia_fixed_dim=0, wild=0 → exp=1
        result = ar.artin_conductor_formula(
            ramified_primes={5: {"inertia_fixed_dim": 0, "wild_part": 0}},
            chi_order=1
        )
        assert result == 5

    def test_artin_conductor_formula_mehrere_primes(self):
        """Mehrere ramifizierte Primzahlen: Führer = Produkt der p-Potenzen."""
        ar = ArtinReciprocity()
        result = ar.artin_conductor_formula(
            ramified_primes={
                2: {"inertia_fixed_dim": 0, "wild_part": 0},
                3: {"inertia_fixed_dim": 0, "wild_part": 0}
            },
            chi_order=1
        )
        # 2^1 · 3^1 = 6
        assert result == 6

    def test_global_artin_map_demo_schlüssel(self):
        """global_artin_map_demo gibt vollständiges Dict zurück."""
        ar = ArtinReciprocity()
        result = ar.global_artin_map_demo(n=5)
        for key in ("extension", "galois_group", "galois_elements",
                    "frobenius_at_primes", "kronecker_weber"):
            assert key in result

    def test_global_artin_map_demo_galois_gruppe_zyklotomisch(self):
        """Galois-Gruppe von Q(ζ_5)/Q hat φ(5)=4 Elemente."""
        ar = ArtinReciprocity()
        result = ar.global_artin_map_demo(n=5)
        galois_elements = result["galois_elements"]
        assert len(galois_elements) == 4  # φ(5) = 4
        # Elemente sind {1, 2, 3, 4} (relativ prim zu 5)
        assert set(galois_elements) == {1, 2, 3, 4}

    def test_global_artin_map_demo_frobenius_bei_2(self):
        """Frobenius_2 für Q(ζ_5)/Q: Frob_2 = 2 mod 5."""
        ar = ArtinReciprocity()
        result = ar.global_artin_map_demo(n=5)
        frob_at_2 = result["frobenius_at_primes"].get(2)
        assert frob_at_2 is not None
        assert frob_at_2["Frob_p"] == 2 % 5


# ===========================================================================
# KLASSE: TraceFormula
# ===========================================================================

class TestTraceFormula:
    """Testet die Arthur-Selberg-Spurformel (numerische Demo)."""

    def test_geometric_side_demo_leer(self):
        """Leere Konjugationsklassen-Liste ergibt total 0."""
        tf = TraceFormula()
        result = tf.geometric_side_demo([])
        assert abs(result["geometric_side_total"]) < 1e-12

    def test_geometric_side_demo_identität(self):
        """Identitäts-Beitrag: Vol(Γ\\H) / (4π)."""
        tf = TraceFormula()
        vol_quotient = math.pi / 3  # SL_2(Z)
        classes = [{"type": "identity", "vol_quotient": vol_quotient}]
        result = tf.geometric_side_demo(classes)
        expected = vol_quotient / (4 * math.pi)
        assert abs(result["geometric_side_total"].real - expected) < 1e-10

    def test_geometric_side_demo_hyperbolisch(self):
        """Hyperbolischer Beitrag: log(N_P) / (√N_P - 1/√N_P) > 0."""
        tf = TraceFormula()
        N_P = 4.0
        classes = [{"type": "hyperbolic", "norm": N_P}]
        result = tf.geometric_side_demo(classes)
        expected = math.log(N_P) / (math.sqrt(N_P) - 1 / math.sqrt(N_P))
        assert abs(result["geometric_side_total"].real - expected) < 1e-10

    def test_geometric_side_demo_num_classes(self):
        """num_classes gibt die Anzahl der übergebenen Klassen zurück."""
        tf = TraceFormula()
        classes = [
            {"type": "identity", "vol_quotient": 1.0},
            {"type": "hyperbolic", "norm": 4.0},
        ]
        result = tf.geometric_side_demo(classes)
        assert result["num_classes"] == 2

    def test_spectral_side_demo_leer(self):
        """Leere Eigenwert-Liste ergibt spektrale Seite 0."""
        tf = TraceFormula()
        result = tf.spectral_side_demo([])
        assert result["spectral_side"] == 0.0

    def test_spectral_side_demo_positive_eigenwerte(self):
        """Eigenwerte > 1/4 erzeugen positive spektrale Summe."""
        tf = TraceFormula()
        # Bekannte Maass-Form-Eigenwerte für SL_2(Z)
        eigenvalues = [91.14, 190.13, 281.66]
        result = tf.spectral_side_demo(eigenvalues)
        spectral_val = result["spectral_side_approx"]
        assert isinstance(spectral_val, float)
        assert spectral_val > 0

    def test_spectral_side_demo_weyl_schätzung(self):
        """Weyl-Gesetz: N(T) ~ Area(Γ\\H)·T/(4π) = T/12."""
        tf = TraceFormula()
        T = 100.0
        eigenvalues = [T]
        result = tf.spectral_side_demo(eigenvalues)
        weyl = result["weyl_estimate_N(T)"]
        # Für SL_2(Z): Area = π/3, Weyl ~ π/3 · T / (4π) = T/12
        assert abs(weyl - T / 12) < 1.0  # grobe Prüfung

    def test_trace_formula_identity_check_schlüssel(self):
        """trace_formula_identity_check gibt vollständiges Dict zurück."""
        tf = TraceFormula()
        result = tf.trace_formula_identity_check(N=1, weight=2)
        for key in ("N", "weight", "index_in_SL2Z", "area_fundamental_domain",
                    "dim_S_k_approx", "weyl_estimate"):
            assert key in result

    def test_trace_formula_identity_check_n1_index1(self):
        """Für N=1 ist der Index [SL_2(Z):Γ_0(1)] = 1."""
        tf = TraceFormula()
        result = tf.trace_formula_identity_check(N=1, weight=2)
        assert result["index_in_SL2Z"] == 1

    def test_trace_formula_identity_check_fläche_wächst_mit_n(self):
        """Fundamentalbereichs-Fläche wächst mit N."""
        tf = TraceFormula()
        r1 = tf.trace_formula_identity_check(N=1, weight=2)
        r2 = tf.trace_formula_identity_check(N=11, weight=2)
        assert r2["area_fundamental_domain"] > r1["area_fundamental_domain"]


# ===========================================================================
# KLASSE: LDualGroup
# ===========================================================================

class TestLDualGroup:
    """Testet die Langlands-Dualgruppen-Klasse."""

    def test_dual_group_gl_n(self):
        """dual_group('GL_n') gibt 'GL_n(C)' zurück."""
        ldg = LDualGroup()
        assert ldg.dual_group("GL_n") == "GL_n(C)"

    def test_dual_group_sl_n(self):
        """dual_group('SL_n') gibt 'PGL_n(C)' zurück."""
        ldg = LDualGroup()
        assert ldg.dual_group("SL_n") == "PGL_n(C)"

    def test_dual_group_sp_2n(self):
        """dual_group('Sp_2n') gibt 'SO_{2n+1}(C)' zurück (Typ C_n ↔ B_n)."""
        ldg = LDualGroup()
        result = ldg.dual_group("Sp_2n")
        assert result == "SO_{2n+1}(C)"

    def test_dual_group_so_2n_selbstdual(self):
        """dual_group('SO_2n') gibt 'SO_{2n}(C)' zurück (Typ D_n selbstdual)."""
        ldg = LDualGroup()
        result = ldg.dual_group("SO_2n")
        assert result == "SO_{2n}(C)"

    def test_dual_group_g2_selbstdual(self):
        """dual_group('G_2') ist G_2(C) (außergewöhnliche Gruppe, selbstdual)."""
        ldg = LDualGroup()
        assert ldg.dual_group("G_2") == "G_2(C)"

    def test_dual_group_e8_selbstdual(self):
        """dual_group('E_8') ist E_8(C) (selbstdual)."""
        ldg = LDualGroup()
        assert ldg.dual_group("E_8") == "E_8(C)"

    def test_dual_group_unbekannt_fallback(self):
        """Unbekannte Gruppe: Fallback gibt '^G(C)' zurück."""
        ldg = LDualGroup()
        result = ldg.dual_group("H_5")
        assert "H_5" in result  # Fallback enthält den Gruppen-Namen

    def test_l_group_gl_n_split(self):
        """L-Gruppe von GL_n: Galois-Wirkung ist trivial (split über Q)."""
        ldg = LDualGroup()
        result = ldg.l_group("GL_n")
        assert result["is_split_over_Q"] is True
        assert result["galois_action_on_dual"] == "trivial"

    def test_l_group_struktur(self):
        """l_group gibt Dict mit allen erwarteten Schlüsseln zurück."""
        ldg = LDualGroup()
        result = ldg.l_group("GL_2")
        for key in ("group", "connected_dual", "L_group", "galois_action_on_dual",
                    "is_split_over_Q", "langlands_parameters"):
            assert key in result

    def test_l_group_u_n_nicht_split(self):
        """U_n ist nicht über Q split (nicht-triviale Galois-Wirkung)."""
        ldg = LDualGroup()
        result = ldg.l_group("U_n")
        assert result["is_split_over_Q"] is False

    def test_check_functoriality_sym2_bewiesen(self):
        """sym²-Lift GL_2 → GL_3 ist bewiesen (Gelbart-Jacquet 1978)."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("GL_2", "GL_3", "sym2")
        assert result["status"] == "bewiesen"
        assert "Gelbart" in result["reference"]

    def test_check_functoriality_sym3_bewiesen(self):
        """sym³-Lift GL_2 → GL_4 ist bewiesen (Kim-Shahidi 2002)."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("GL_2", "GL_4", "sym3")
        assert result["status"] == "bewiesen"

    def test_check_functoriality_sym5_offen(self):
        """sym^5-Lift GL_2 → GL_6 ist noch offen."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("GL_2", "GL_6", "sym5")
        assert result["status"] == "offen"

    def test_check_functoriality_base_change_bewiesen(self):
        """Basiswechsel GL_2 → GL_2 ist bewiesen (Arthur-Clozel 1989)."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("GL_2", "GL_2", "base_change")
        assert "bewiesen" in result["status"]

    def test_check_functoriality_unbekannt_fallback(self):
        """Unbekannte Funktorialität gibt Status 'unbekannt / offen' zurück."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("G_2", "E_8", "mystery_lift")
        assert "unbekannt" in result["status"] or "offen" in result["status"]

    def test_check_functoriality_schlüssel(self):
        """check_functoriality gibt Dict mit source, target, lift_type zurück."""
        ldg = LDualGroup()
        result = ldg.check_functoriality("GL_2", "GL_3", "sym2")
        assert result["source"] == "GL_2"
        assert result["target"] == "GL_3"
        assert result["lift_type"] == "sym2"

    def test_pgl_n_dual_ist_sl_n(self):
        """dual_group('PGL_n') ist 'SL_n(C)' (dual zu SL_n → PGL_n)."""
        ldg = LDualGroup()
        assert ldg.dual_group("PGL_n") == "SL_n(C)"

    def test_sp4_dual_ist_so5(self):
        """dual_group('Sp_4') ist 'SO_5(C)' (konkrete sp4 Dual-Gruppe)."""
        ldg = LDualGroup()
        assert ldg.dual_group("Sp_4") == "SO_5(C)"


# ===========================================================================
# INTEGRATIONSTEST: Zusammenspiel der Module
# ===========================================================================

class TestIntegration:
    """Integrationstests, die mehrere Klassen zusammen testen."""

    def test_satake_parameter_konsistenz_mit_local_parameter(self):
        """SatakeIsomorphism und local_langlands_parameter liefern konsistente Ergebnisse."""
        sat = SatakeIsomorphism(n=2)
        p = 5
        a_p = 2.0
        # Via local_langlands_parameter
        params = local_langlands_parameter(p, a_p)
        alpha_lp = params["alpha"]
        # Via SatakeIsomorphism
        alpha_sat, _ = sat.satake_parameters_from_hecke(a_p=a_p, p=p, weight=2)
        # Beide |α_p| sollten √p sein
        assert abs(abs(alpha_lp) - math.sqrt(p)) < 1e-8
        assert abs(abs(alpha_sat) - math.sqrt(p)) < 1e-8

    def test_local_l_faktor_und_global_modularity(self):
        """LocalLanglandsGL2 und GlobalLanglands: Konsistenz der Ramanujan-Schranke."""
        p = 5
        a_p = 2.0
        llgl2 = LocalLanglandsGL2(p=p)
        gl = GlobalLanglands()
        # Satake-Parameter für die Hauptreihe
        alpha = complex(a_p / 2, math.sqrt(4 * p - a_p ** 2) / 2)
        beta = alpha.conjugate()
        # L-Faktor bei s=2
        l_val = llgl2.local_factor_l("principal_series", s=2.0, alpha=alpha, beta=beta)
        assert math.isfinite(abs(l_val))
        # Ramanujan global
        result = gl.modular_curve_wiles_check({p: int(a_p)}, N=32)
        assert result["prime_checks"][p]["satisfies_ramanujan"] is True

    def test_artin_l_function_trivial_char_und_epsilon_bei_halb(self):
        """Artin-L und ε-Faktor: Bei s=1/2 und N=1 ist ε=1."""
        ar = ArtinReciprocity()
        # Trivial-Charakter
        primes = [2, 3, 5, 7, 11, 13]
        chi_trivial = {p: 1.0 for p in primes}
        l_val = ar.artin_l_function(chi_trivial, s=2.0, primes=primes)
        assert l_val.real > 1.0
        # ε-Faktor bei s=1/2 und N=1 ist exakt 1
        eps = epsilon_factor_normalization(conductor=1, s=0.5)
        assert abs(eps - 1.0) < 1e-10

    def test_symmetrische_potenz_und_dualgruppe_konsistenz(self):
        """sym^k-Lift GL_2 → GL_{k+1} stimmt mit Zielgruppe überein."""
        lf = LanglandsFunctoriality()
        ldg = LDualGroup()
        for k in [2, 3, 4]:
            lift_result = lf.symmetric_power_lift(a_p={5: 2}, p=5, k=k)
            target_group = lift_result["results"][5]["sym_k_group"]
            assert target_group == f"GL_{k + 1}"
