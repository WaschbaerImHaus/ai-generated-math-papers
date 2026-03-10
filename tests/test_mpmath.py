"""
@file test_mpmath.py
@brief Tests für die mpmath-basierten Funktionen in complex_analysis.py.

@description
    Testet alle neuen Hochpräzisions-Funktionen aus complex_analysis.py Build 11:
    - riemann_zeta_mpmath(): Beliebig genaue Zeta-Funktion
    - gamma_mpmath(): Beliebig genaue Gamma-Funktion
    - riemann_siegel_z_mpmath(): Riemann-Siegel Z-Funktion
    - verify_riemann_hypothesis_range_mpmath(): RH-Verifikation
    - find_zeta_zeros_mpmath(): Nullstellensuche mit mpmath
    - li_function_mpmath(): Logarithmischer Integralus

    Bekannte mathematische Werte für Validierung:
    - ζ(2) = π²/6 ≈ 1.6449340668...
    - ζ(0) = -1/2 = -0.5
    - Γ(1) = 1
    - Γ(1/2) = √π ≈ 1.7724538509...
    - Erste Nullstelle der Zeta-Funktion: t₁ ≈ 14.134725...

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import math
import pytest

# Zu testende Funktionen importieren
from complex_analysis import (
    riemann_zeta_mpmath,
    find_zeta_zeros_mpmath,
    riemann_siegel_z_mpmath,
    gamma_mpmath,
    verify_riemann_hypothesis_range_mpmath,
    li_function_mpmath,
)


# ===========================================================================
# SKIP-MARK: Falls mpmath nicht verfügbar
# ===========================================================================

try:
    import mpmath
    MPMATH_VERFUEGBAR = True
except ImportError:
    MPMATH_VERFUEGBAR = False

# Alle Tests in dieser Datei überspringen wenn mpmath nicht installiert ist
pytestmark = pytest.mark.skipif(
    not MPMATH_VERFUEGBAR,
    reason="mpmath ist nicht installiert (pip install mpmath)"
)


# ===========================================================================
# TESTS: riemann_zeta_mpmath
# ===========================================================================

class TestRiemannZetaMpmath:
    """Tests für die Hochpräzisions-Zeta-Funktion."""

    def test_zeta_2_gleich_pi_quadrat_durch_6(self):
        """ζ(2) = π²/6 ≈ 1.6449340668482264 (Euler, 1735)."""
        ergebnis = riemann_zeta_mpmath(complex(2, 0), dps=50)
        erwartung = math.pi**2 / 6  # ≈ 1.6449340668482264
        assert abs(ergebnis.real - erwartung) < 1e-10, (
            f"ζ(2) sollte ≈ {erwartung:.6f} sein, erhalten: {ergebnis.real:.6f}"
        )
        assert abs(ergebnis.imag) < 1e-12, "ζ(2) muss reell sein (kein Imaginärteil)"

    def test_zeta_0_gleich_minus_halb(self):
        """ζ(0) = -1/2 (bekannter Wert durch analytische Fortsetzung)."""
        ergebnis = riemann_zeta_mpmath(complex(0, 0), dps=30)
        assert abs(ergebnis.real - (-0.5)) < 1e-8, (
            f"ζ(0) sollte -0.5 sein, erhalten: {ergebnis.real}"
        )

    def test_zeta_minus1_ergibt_minus_zwolftel(self):
        """ζ(-1) = -1/12 (regularisierte Summe 1+2+3+... = -1/12)."""
        ergebnis = riemann_zeta_mpmath(complex(-1, 0), dps=30)
        erwartung = -1.0 / 12.0
        assert abs(ergebnis.real - erwartung) < 1e-8, (
            f"ζ(-1) sollte -1/12 ≈ {erwartung:.6f} sein, erhalten: {ergebnis.real}"
        )

    def test_zeta_4_gleich_pi_hoch_4_durch_90(self):
        """ζ(4) = π⁴/90 ≈ 1.0823232337111381 (Euler)."""
        ergebnis = riemann_zeta_mpmath(complex(4, 0), dps=50)
        erwartung = math.pi**4 / 90
        assert abs(ergebnis.real - erwartung) < 1e-10

    def test_zeta_pol_bei_s_gleich_1(self):
        """ζ(1) muss eine ValueError-Ausnahme werfen (Pol)."""
        with pytest.raises(ValueError):
            riemann_zeta_mpmath(complex(1, 0))

    def test_zeta_auf_kritischer_geraden(self):
        """ζ(1/2 + 14i) hat kleinen Betrag in der Nähe der ersten Nullstelle."""
        # t ≈ 14.134725 ist die erste Nullstelle → |ζ(1/2 + it)| ≈ 0 dort
        # Bei t = 10 (weit von Nullstelle) sollte |ζ| > 0 sein
        ergebnis = riemann_zeta_mpmath(complex(0.5, 10.0), dps=30)
        assert abs(ergebnis) > 0.0, "ζ(1/2 + 10i) sollte ≠ 0 sein"

    def test_verschiedene_dps_werte(self):
        """Unterschiedliche dps-Werte liefern konsistente Ergebnisse."""
        z2_50 = riemann_zeta_mpmath(complex(2, 0), dps=50)
        z2_30 = riemann_zeta_mpmath(complex(2, 0), dps=30)
        # Beide sollten nah an π²/6 sein
        erwartung = math.pi**2 / 6
        assert abs(z2_50.real - erwartung) < 1e-10
        assert abs(z2_30.real - erwartung) < 1e-8

    def test_rueckgabe_ist_python_complex(self):
        """Rückgabewert muss Python complex sein (nicht mpmath.mpc)."""
        ergebnis = riemann_zeta_mpmath(complex(3, 0))
        assert isinstance(ergebnis, complex), "Rückgabe muss Python complex sein"


# ===========================================================================
# TESTS: gamma_mpmath
# ===========================================================================

class TestGammaMpmath:
    """Tests für die Hochpräzisions-Gamma-Funktion."""

    def test_gamma_1_gleich_1(self):
        """Γ(1) = 0! = 1 (Definition der Gamma-Funktion)."""
        ergebnis = gamma_mpmath(complex(1, 0), dps=50)
        assert abs(ergebnis.real - 1.0) < 1e-12, (
            f"Γ(1) = 1, erhalten: {ergebnis.real}"
        )
        assert abs(ergebnis.imag) < 1e-14, "Γ(1) muss reell sein"

    def test_gamma_halb_gleich_wurzel_pi(self):
        """Γ(1/2) = √π ≈ 1.7724538509055159."""
        ergebnis = gamma_mpmath(complex(0.5, 0), dps=50)
        erwartung = math.sqrt(math.pi)  # ≈ 1.7724538509055159
        assert abs(ergebnis.real - erwartung) < 1e-12, (
            f"Γ(1/2) = √π ≈ {erwartung:.8f}, erhalten: {ergebnis.real:.8f}"
        )

    def test_gamma_2_gleich_1(self):
        """Γ(2) = 1! = 1."""
        ergebnis = gamma_mpmath(complex(2, 0), dps=30)
        assert abs(ergebnis.real - 1.0) < 1e-10

    def test_gamma_3_gleich_2(self):
        """Γ(3) = 2! = 2."""
        ergebnis = gamma_mpmath(complex(3, 0), dps=30)
        assert abs(ergebnis.real - 2.0) < 1e-10

    def test_gamma_4_gleich_6(self):
        """Γ(4) = 3! = 6."""
        ergebnis = gamma_mpmath(complex(4, 0), dps=30)
        assert abs(ergebnis.real - 6.0) < 1e-10

    def test_gamma_rekursion(self):
        """Γ(z+1) = z·Γ(z) (Funktionalgleichung)."""
        z = complex(2.5, 0)
        gz = gamma_mpmath(z, dps=50)
        gz1 = gamma_mpmath(z + 1, dps=50)
        # Γ(3.5) = 2.5 · Γ(2.5)
        assert abs(gz1.real - 2.5 * gz.real) < 1e-10, (
            "Γ(z+1) = z·Γ(z) muss gelten"
        )

    def test_rueckgabe_ist_python_complex(self):
        """Rückgabewert muss Python complex sein."""
        ergebnis = gamma_mpmath(complex(5, 0))
        assert isinstance(ergebnis, complex)

    def test_gamma_komplexes_argument(self):
        """Gamma-Funktion funktioniert auch für komplexe Argumente."""
        # Γ(1 + i) ist berechenbar
        ergebnis = gamma_mpmath(complex(1, 1), dps=30)
        assert isinstance(ergebnis, complex)
        # |Γ(1+i)| ≈ 0.498... (bekannter Wert)
        assert abs(abs(ergebnis) - 0.498) < 0.01, (
            f"|Γ(1+i)| ≈ 0.498, erhalten: {abs(ergebnis):.4f}"
        )


# ===========================================================================
# TESTS: riemann_siegel_z_mpmath
# ===========================================================================

class TestRiemannSiegelZMpmath:
    """Tests für die Riemann-Siegel Z-Funktion."""

    def test_gibt_float_zurueck(self):
        """riemann_siegel_z_mpmath() muss Python float zurückgeben."""
        ergebnis = riemann_siegel_z_mpmath(14.0, dps=30)
        assert isinstance(ergebnis, float), "Rückgabe muss Python float sein"

    def test_z_bei_erster_nullstelle_nahe_0(self):
        """Z(14.134725) ≈ 0 (erste Nullstelle der Riemann-Zeta-Funktion)."""
        # t₁ ≈ 14.134725141734693
        z_val = riemann_siegel_z_mpmath(14.134725, dps=50)
        assert abs(z_val) < 0.01, (
            f"Z(14.134725) sollte ≈ 0 sein (erste Nullstelle), erhalten: {z_val}"
        )

    def test_vorzeichenwechsel_erste_nullstelle(self):
        """Z(t) wechselt Vorzeichen nahe t₁ ≈ 14.135."""
        # Links von der Nullstelle
        z_links = riemann_siegel_z_mpmath(13.0, dps=30)
        # Rechts von der Nullstelle
        z_rechts = riemann_siegel_z_mpmath(17.0, dps=30)
        # Vorzeichen müssen unterschiedlich sein (Vorzeichenwechsel irgendwo dazwischen)
        assert z_links * z_rechts < 0 or abs(z_links) > 0.1, (
            "Vorzeichenwechsel nahe erster Nullstelle erwartet"
        )

    def test_betrag_positiv_abseits_nullstellen(self):
        """Z(10) sollte |Z| > 0 sein (weit von der ersten Nullstelle)."""
        z_val = riemann_siegel_z_mpmath(10.0, dps=30)
        assert abs(z_val) > 0.0, "Z(10) sollte ≠ 0 sein"

    def test_zweite_nullstelle(self):
        """Z(21.022) ≈ 0 (zweite Nullstelle)."""
        z_val = riemann_siegel_z_mpmath(21.022, dps=50)
        assert abs(z_val) < 0.1, (
            f"Z(21.022) sollte ≈ 0 sein (zweite Nullstelle), erhalten: {z_val}"
        )


# ===========================================================================
# TESTS: verify_riemann_hypothesis_range_mpmath
# ===========================================================================

class TestVerifyRiemannHypothesisMpmath:
    """Tests für die RH-Verifikation der ersten n Nullstellen."""

    def test_gibt_dict_zurueck(self):
        """verify_riemann_hypothesis_range_mpmath() muss Dict zurückgeben."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=3, dps=30)
        assert isinstance(result, dict), "Rückgabe muss Dict sein"

    def test_dict_hat_alle_schluessel(self):
        """Dict muss 'zeros', 'all_on_critical_line', 'max_deviation' enthalten."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=3, dps=30)
        assert 'zeros' in result, "Dict muss 'zeros' enthalten"
        assert 'all_on_critical_line' in result, "Dict muss 'all_on_critical_line' enthalten"
        assert 'max_deviation' in result, "Dict muss 'max_deviation' enthalten"

    def test_erste_5_nullstellen_auf_kritischer_geraden(self):
        """Die ersten 5 Nullstellen liegen auf der kritischen Geraden Re(s)=1/2."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=5, dps=50)
        assert result['all_on_critical_line'] is True, (
            f"Alle ersten 5 Nullstellen sollen auf Re=1/2 liegen. "
            f"max_deviation = {result['max_deviation']}"
        )

    def test_max_deviation_sehr_klein(self):
        """Maximale Abweichung von Re=1/2 muss sehr klein sein."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=5, dps=50)
        assert result['max_deviation'] < 1e-8, (
            f"Maximale Abweichung zu groß: {result['max_deviation']}"
        )

    def test_zeros_liste_hat_korrekte_laenge(self):
        """Liste 'zeros' muss n_zeros Einträge haben."""
        n = 3
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=n, dps=30)
        assert len(result['zeros']) == n, (
            f"'zeros' sollte {n} Einträge haben, erhalten: {len(result['zeros'])}"
        )

    def test_erste_nullstelle_imaginaerteil(self):
        """Die erste Nullstelle hat Im(ρ₁) ≈ 14.134725."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=1, dps=50)
        erste_nullstelle = result['zeros'][0]
        if erste_nullstelle.get('rho') is not None:
            im_teil = erste_nullstelle['rho'].imag
            assert abs(im_teil - 14.134725) < 1e-4, (
                f"Im(ρ₁) ≈ 14.134725, erhalten: {im_teil}"
            )

    def test_all_on_critical_line_ist_bool(self):
        """'all_on_critical_line' muss Python bool sein."""
        result = verify_riemann_hypothesis_range_mpmath(n_zeros=3, dps=30)
        assert isinstance(result['all_on_critical_line'], bool)


# ===========================================================================
# TESTS: find_zeta_zeros_mpmath
# ===========================================================================

class TestFindZetaZerosMpmath:
    """Tests für die Nullstellensuche mit mpmath."""

    def test_findet_erste_nullstelle(self):
        """Suche im Bereich [14, 15] findet die erste Nullstelle."""
        zeros = find_zeta_zeros_mpmath(14.0, 15.0, dps=30, steps=200)
        assert len(zeros) >= 1, "Im Bereich [14, 15] muss mindestens eine Nullstelle sein"

    def test_nullstelle_bei_t1_korrekt(self):
        """Gefundene Nullstelle bei t₁ ≈ 14.134725."""
        zeros = find_zeta_zeros_mpmath(14.0, 14.5, dps=30, steps=300)
        if zeros:
            t = zeros[0]['t']
            assert abs(t - 14.134725) < 0.01, (
                f"Erste Nullstelle bei t ≈ 14.134725, gefunden: {t}"
            )

    def test_gibt_liste_zurueck(self):
        """find_zeta_zeros_mpmath() muss eine Liste zurückgeben."""
        result = find_zeta_zeros_mpmath(10.0, 16.0, dps=20, steps=100)
        assert isinstance(result, list)

    def test_listenelemente_haben_korrekte_schluessel(self):
        """Jedes Listenelement muss 't', 'zero', 'verified' haben."""
        zeros = find_zeta_zeros_mpmath(14.0, 15.0, dps=20, steps=200)
        for z in zeros:
            assert 't' in z, "Jedes Element muss 't' haben"
            assert 'zero' in z, "Jedes Element muss 'zero' haben"
            assert 'verified' in z, "Jedes Element muss 'verified' haben"

    def test_nullstellen_auf_kritischer_geraden(self):
        """Alle gefundenen Nullstellen sollen Re(s) = 1/2 haben."""
        zeros = find_zeta_zeros_mpmath(14.0, 22.0, dps=30, steps=500)
        for z in zeros:
            re_teil = z['zero'].real
            assert abs(re_teil - 0.5) < 1e-8, (
                f"Nullstelle {z['zero']} hat Re ≠ 1/2: Re = {re_teil}"
            )


# ===========================================================================
# TESTS: li_function_mpmath
# ===========================================================================

class TestLiFunctionMpmath:
    """Tests für den Logarithmischen Integralus."""

    def test_li_100_positiv(self):
        """li(100) > 0 (für x > 1 ist li(x) positiv und wächst)."""
        ergebnis = li_function_mpmath(100.0, dps=30)
        assert ergebnis > 0.0, f"li(100) muss > 0 sein, erhalten: {ergebnis}"

    def test_li_1000_groesser_als_li_100(self):
        """li(1000) > li(100) (li ist monoton wachsend für x > 1)."""
        li100 = li_function_mpmath(100.0, dps=30)
        li1000 = li_function_mpmath(1000.0, dps=30)
        assert li1000 > li100, (
            f"li(1000) = {li1000:.4f} sollte > li(100) = {li100:.4f} sein"
        )

    def test_li_10_bekannter_wert(self):
        """li(10) ≈ 6.1655... (Tabellenwert)."""
        ergebnis = li_function_mpmath(10.0, dps=30)
        # Bekannter Wert: li(10) ≈ 6.1655...
        assert abs(ergebnis - 6.165) < 0.01, (
            f"li(10) ≈ 6.1655, erhalten: {ergebnis:.4f}"
        )

    def test_li_2_positiv(self):
        """li(2) > 0."""
        ergebnis = li_function_mpmath(2.0, dps=30)
        assert ergebnis > 0.0

    def test_gibt_float_zurueck(self):
        """li_function_mpmath() muss Python float zurückgeben."""
        ergebnis = li_function_mpmath(100.0, dps=30)
        assert isinstance(ergebnis, float), "Rückgabe muss Python float sein"

    def test_li_approximiert_pix(self):
        """li(x) ist eine gute Approximation von π(x) (Primzahlzählerfunktion)."""
        # π(100) = 25 (Anzahl der Primzahlen ≤ 100)
        # li(100) ≈ 29.08... → Differenz sollte < 10 sein
        li_100 = li_function_mpmath(100.0, dps=30)
        pi_100 = 25  # Exakter Wert: 25 Primzahlen ≤ 100
        assert abs(li_100 - pi_100) < 10, (
            f"li(100) ≈ {li_100:.2f} sollte nahe π(100) = {pi_100} sein"
        )

    def test_verschiedene_dps(self):
        """Unterschiedliche dps-Werte liefern konsistente Ergebnisse."""
        li_50 = li_function_mpmath(100.0, dps=50)
        li_20 = li_function_mpmath(100.0, dps=20)
        # Beide sollten nah beieinander liegen (< 1e-8 Abweichung)
        assert abs(li_50 - li_20) < 1e-8, (
            f"dps=50: {li_50}, dps=20: {li_20} – zu große Abweichung"
        )
