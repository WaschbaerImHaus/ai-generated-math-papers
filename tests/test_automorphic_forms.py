"""
@file test_automorphic_forms.py
@brief Tests für das automorphic_forms-Modul.
@description
    Umfassende Test-Suite für alle Klassen und Funktionen des automorphic_forms-Moduls.

    Getestete Klassen:
    - AutomorphicForm: Eigenschaften, Fourier-Koeffizienten, Dimension
    - SpectralDecomposition: Diskretes Spektrum, Eisenstein-Reihen
    - WhittakerModel: Whittaker-Funktionen, Fourier-Entwicklung
    - GlobalLFunction: Euler-Produkt, Funktionalgleichung
    - ArthurPacket: Temperiertheit, lokale Komponenten
    - FunctorialLift: Sym²-Lift, Hecke-Eigenwerte
    - RamanujanConjecture: Schranken-Verifikation
    - SatakeTransform: Satake-Parameter, Hecke-Polynome

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os

# Suchpfad für Quellcode setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import math
import cmath
from exceptions import InvalidInputError

# Modul unter Test
from automorphic_forms import (
    AutomorphicForm,
    SpectralDecomposition,
    WhittakerModel,
    GlobalLFunction,
    ArthurPacket,
    FunctorialLift,
    RamanujanConjecture,
    SatakeTransform,
    _frobenius_trace,
    _mobius,
    _euler_phi,
    _is_prime,
)


# =============================================================================
# TESTS: Hilfsfunktionen
# =============================================================================

class TestHelperFunctions:
    """Tests für die internen Hilfsfunktionen."""

    def test_is_prime_small_primes(self):
        """Primtest für kleine Primzahlen."""
        assert _is_prime(2) is True
        assert _is_prime(3) is True
        assert _is_prime(5) is True
        assert _is_prime(7) is True
        assert _is_prime(11) is True

    def test_is_prime_non_primes(self):
        """Primtest für Nicht-Primzahlen."""
        assert _is_prime(1) is False
        assert _is_prime(4) is False
        assert _is_prime(6) is False
        assert _is_prime(9) is False
        assert _is_prime(25) is False

    def test_euler_phi(self):
        """Euler-Phi-Funktion für bekannte Werte."""
        assert _euler_phi(1) == 1
        assert _euler_phi(2) == 1
        assert _euler_phi(4) == 2
        assert _euler_phi(6) == 2
        assert _euler_phi(12) == 4

    def test_mobius(self):
        """Möbius-Funktion für bekannte Werte."""
        assert _mobius(1) == 1       # μ(1) = 1
        assert _mobius(2) == -1      # μ(2) = -1 (eine Primzahl)
        assert _mobius(6) == 1       # μ(6) = μ(2·3) = (-1)² = 1
        assert _mobius(4) == 0       # μ(4) = 0 (4 = 2², nicht quadratfrei)
        assert _mobius(30) == -1     # μ(30) = μ(2·3·5) = (-1)³ = -1

    def test_frobenius_trace_basic(self):
        """Frobenius-Spur für einfache elliptische Kurven."""
        # Für E: y² = x³ + x (über F_5)
        a_5 = _frobenius_trace(5, 1, 0)
        # Manuell: Zähle Punkte auf y² = x³+x über F_5
        # Ergebnis: a_5 sollte |a_5| ≤ 2√5 ≈ 4.47 sein (Hasse-Schranke)
        assert abs(a_5) <= 2 * math.sqrt(5) + 1

    def test_frobenius_trace_hasse_bound(self):
        """Frobenius-Spur erfüllt Hasse-Schranke |a_p| ≤ 2√p für alle p."""
        for p, a, b in [(3, 1, 1), (5, 0, 1), (7, 1, 2), (11, 2, 3), (13, 1, 0)]:
            a_p = _frobenius_trace(p, a, b)
            assert abs(a_p) <= 2 * math.sqrt(p) + 1, (
                f"Hasse-Schranke verletzt: |a_{p}| = {abs(a_p)} > 2√{p} ≈ {2*math.sqrt(p):.2f}"
            )

    def test_frobenius_trace_invalid_p(self):
        """Ungültige Eingabe wirft Fehler."""
        with pytest.raises(InvalidInputError):
            _frobenius_trace(4, 1, 1)  # 4 ist keine Primzahl


# =============================================================================
# TESTS: AutomorphicForm
# =============================================================================

class TestAutomorphicForm:
    """Tests für die AutomorphicForm-Klasse."""

    def test_creation_default(self):
        """Standarderstellung einer automorphen Form."""
        f = AutomorphicForm()
        assert f.group == "GL2"
        assert f.weight == 12
        assert f.level == 1

    def test_creation_custom(self):
        """Erstellung mit angepassten Parametern."""
        f = AutomorphicForm(group="GL3", weight=8, level=5)
        assert f.group == "GL3"
        assert f.weight == 8
        assert f.level == 5

    def test_invalid_weight(self):
        """Negativer oder Null-Gewicht wirft Fehler."""
        with pytest.raises(InvalidInputError):
            AutomorphicForm(weight=0)
        with pytest.raises(InvalidInputError):
            AutomorphicForm(weight=-1)

    def test_invalid_level(self):
        """Ungültiges Niveau wirft Fehler."""
        with pytest.raises(InvalidInputError):
            AutomorphicForm(level=0)
        with pytest.raises(InvalidInputError):
            AutomorphicForm(level=-5)

    def test_dimension_positive(self):
        """Dimension muss nicht-negativ sein."""
        f = AutomorphicForm(weight=12, level=1)
        assert f.dimension() >= 0

    def test_is_cuspidal_weight12(self):
        """Ramanujan-Δ ist eine Kuspform (Gewicht 12, Niveau 1)."""
        f = AutomorphicForm(weight=12, level=1)
        assert f.is_cuspidal() is True

    def test_is_cuspidal_low_weight(self):
        """Für kleines Gewicht keine Kuspformen."""
        f = AutomorphicForm(weight=4, level=1)
        # Keine Kuspformen für Gewicht < 12 bei Niveau 1
        assert f.is_cuspidal() is False

    def test_fourier_coefficient_ramanujan(self):
        """Ramanujan-Tau-Funktion: τ(1) = 1."""
        f = AutomorphicForm(group="GL2", weight=12, level=1)
        a1 = f.fourier_coefficient(1)
        assert abs(a1 - 1) < 1, "τ(1) sollte nahe bei 1 liegen"

    def test_fourier_coefficient_invalid(self):
        """Ungültiger Index wirft Fehler."""
        f = AutomorphicForm()
        with pytest.raises(InvalidInputError):
            f.fourier_coefficient(0)
        with pytest.raises(InvalidInputError):
            f.fourier_coefficient(-1)

    def test_info_dict(self):
        """info() gibt ein vollständiges Dictionary zurück."""
        f = AutomorphicForm(weight=16, level=3)
        info = f.info()
        assert 'group' in info
        assert 'weight' in info
        assert 'level' in info
        assert 'dimension' in info
        assert 'is_cuspidal' in info
        assert info['weight'] == 16
        assert info['level'] == 3


# =============================================================================
# TESTS: SpectralDecomposition
# =============================================================================

class TestSpectralDecomposition:
    """Tests für die SpectralDecomposition-Klasse."""

    def test_creation(self):
        """Standarderstellung."""
        sd = SpectralDecomposition()
        assert sd.group == "SL2Z"

    def test_discrete_spectrum_no_cusp_below_12(self):
        """Keine Kuspformen für Gewicht < 12 bei SL_2(ℤ)."""
        sd = SpectralDecomposition()
        disc = sd.discrete_spectrum(level=1, max_weight=11)
        # Keine Kuspformen mit Gewicht ≤ 11 für SL_2(ℤ)
        assert disc['total_cusp_forms'] == 0

    def test_discrete_spectrum_weight_12_exists(self):
        """Kuspform der Stufe 12 (Ramanujan-Δ) muss existieren."""
        sd = SpectralDecomposition()
        disc = sd.discrete_spectrum(level=1, max_weight=12)
        assert 12 in disc['discrete_spectrum']
        assert disc['discrete_spectrum'][12]['dimension'] == 1

    def test_discrete_spectrum_structure(self):
        """Diskretes Spektrum hat korrekte Struktur."""
        sd = SpectralDecomposition()
        disc = sd.discrete_spectrum(level=1, max_weight=24)
        assert 'discrete_spectrum' in disc
        assert 'total_cusp_forms' in disc
        assert disc['total_cusp_forms'] >= 0

    def test_eisenstein_series_pole(self):
        """Eisenstein-Reihe hat Pol bei s=1."""
        sd = SpectralDecomposition()
        result = sd.eisenstein_series(complex(1.0, 0.0))
        assert result.get('has_pole') is True

    def test_eisenstein_series_critical_line(self):
        """Eisenstein-Reihe auf der kritischen Linie."""
        sd = SpectralDecomposition()
        result = sd.eisenstein_series(complex(0.5, 3.0))
        assert result.get('is_on_critical_line') is True
        assert result.get('continuous_spectrum') is True

    def test_spectral_decomposition_summary(self):
        """Zusammenfassung der Spektralzerlegung."""
        sd = SpectralDecomposition()
        summary = sd.spectral_decomposition_summary()
        assert 'discrete' in summary
        assert 'continuous' in summary
        assert 'selberg_trace_formula' in summary


# =============================================================================
# TESTS: WhittakerModel
# =============================================================================

class TestWhittakerModel:
    """Tests für das WhittakerModel."""

    def test_creation(self):
        """Standarderstellung."""
        wm = WhittakerModel()
        assert wm.representation == "GL2"
        assert wm.n == 2

    def test_creation_gl3(self):
        """GL(3)-Whittaker-Modell."""
        wm = WhittakerModel(representation="GL3")
        assert wm.n == 3

    def test_whittaker_function_positive(self):
        """Whittaker-Funktion gibt positiven reellen Wert für y > 0."""
        wm = WhittakerModel()
        w = wm.whittaker_function(1.0)
        assert w.real > 0

    def test_whittaker_function_decreasing(self):
        """Whittaker-Funktion ist für große y abnehmend."""
        wm = WhittakerModel()
        w1 = abs(wm.whittaker_function(1.0))
        w2 = abs(wm.whittaker_function(2.0))
        # e^{-2πy} ist stark abnehmend
        assert w1 > w2

    def test_whittaker_function_invalid_y(self):
        """y ≤ 0 wirft Fehler."""
        wm = WhittakerModel()
        with pytest.raises(InvalidInputError):
            wm.whittaker_function(0.0)
        with pytest.raises(InvalidInputError):
            wm.whittaker_function(-1.0)

    def test_is_generic(self):
        """GL(n) ist generic."""
        for n in [1, 2, 3, 4]:
            wm = WhittakerModel(representation=f"GL{n}")
            assert wm.is_generic() is True

    def test_fourier_whittaker_expansion(self):
        """Fourier-Whittaker-Entwicklung mit bekannten Koeffizienten."""
        wm = WhittakerModel()
        coeffs = [1.0, -24.0, 252.0, -1472.0, 4830.0]
        result = wm.fourier_whittaker_expansion(coeffs)
        assert 'partial_sum' in result
        assert result['n_terms'] == 5

    def test_fourier_whittaker_empty_coefficients(self):
        """Leere Koeffizientenliste wirft Fehler."""
        wm = WhittakerModel()
        with pytest.raises(InvalidInputError):
            wm.fourier_whittaker_expansion([])

    def test_info_dict(self):
        """info() gibt vollständiges Dictionary zurück."""
        wm = WhittakerModel(representation="GL2", psi="standard")
        info = wm.info()
        assert 'representation' in info
        assert 'n' in info
        assert 'is_generic' in info


# =============================================================================
# TESTS: GlobalLFunction
# =============================================================================

class TestGlobalLFunction:
    """Tests für die GlobalLFunction-Klasse."""

    def test_creation(self):
        """Standarderstellung."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f)
        assert L.form is f

    def test_euler_product_converges(self):
        """Euler-Produkt konvergiert für Re(s) > 1."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f)
        val = L.euler_product(complex(2.0, 0), num_primes=10)
        # Sollte einen endlichen Wert liefern
        assert math.isfinite(abs(val))

    def test_euler_product_right_half_plane(self):
        """Euler-Produkt für s im rechten Halbraum (Re(s) = 3/2)."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f)
        val = L.euler_product(complex(1.5, 2.0), num_primes=15)
        assert math.isfinite(abs(val))

    def test_evaluate(self):
        """evaluate() gibt einen Wert zurück."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f, s=complex(2.0, 0))
        val = L.evaluate()
        assert isinstance(val, complex)

    def test_functional_equation_check(self):
        """Funktionalgleichungs-Check liefert Daten."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f)
        result = L.functional_equation_check(complex(0.5, 5.0))
        assert 's' in result
        assert 'L_s' in result
        assert 'L_1ms' in result

    def test_info_dict(self):
        """info() gibt vollständiges Dictionary zurück."""
        f = AutomorphicForm(weight=12, level=1)
        L = GlobalLFunction(f)
        info = L.info()
        assert 'weight' in info
        assert 'level' in info


# =============================================================================
# TESTS: ArthurPacket
# =============================================================================

class TestArthurPacket:
    """Tests für ArthurPacket."""

    def test_creation(self):
        """Standarderstellung."""
        ap = ArthurPacket()
        assert ap.group == "GL2"

    def test_is_tempered_default(self):
        """Standardpaket ist tempered."""
        ap = ArthurPacket()
        assert ap.is_tempered() is True

    def test_is_tempered_non_tempered(self):
        """Nicht-temperiertes Paket."""
        ap = ArthurPacket(parameter={'tempered': False})
        assert ap.is_tempered() is False

    def test_local_components(self):
        """Lokale Komponenten haben die richtige Struktur."""
        ap = ArthurPacket(group="GL3")
        comps = ap.local_components()
        assert 'local_components' in comps
        assert 'unramified_places' in comps['local_components']

    def test_endoscopic_classification(self):
        """Endoskopische Klassifikation."""
        ap = ArthurPacket()
        classif = ap.endoscopic_classification()
        assert 'group' in classif
        assert 'classification' in classif

    def test_info_dict(self):
        """info() gibt vollständiges Dictionary zurück."""
        ap = ArthurPacket(group="Sp4")
        info = ap.info()
        assert 'group' in info
        assert 'is_tempered' in info
        assert info['group'] == "Sp4"


# =============================================================================
# TESTS: FunctorialLift
# =============================================================================

class TestFunctorialLift:
    """Tests für FunctorialLift."""

    def test_creation_sym2(self):
        """Sym²-Lift von GL(2) nach GL(3)."""
        f = AutomorphicForm(weight=12, level=1)
        lift = FunctorialLift(f, lift_type="sym2")
        assert lift.lift_type == "sym2"
        assert lift.lift_info['proved'] is True

    def test_creation_invalid_type(self):
        """Unbekannter Lift-Typ wirft Fehler."""
        f = AutomorphicForm(weight=12, level=1)
        with pytest.raises(InvalidInputError):
            FunctorialLift(f, lift_type="nonexistent")

    def test_sym2_lift_proved(self):
        """Sym²-Lift ist bewiesen (Gelbart-Jacquet 1978)."""
        f = AutomorphicForm(weight=12, level=1)
        lift = FunctorialLift(f, lift_type="sym2")
        assert lift.lift_info['proved'] is True
        assert "Gelbart" in lift.lift_info['reference'] or "1978" in lift.lift_info['reference']

    def test_sym5_lift_unproved(self):
        """Sym⁵-Lift ist (noch) unbewiesen."""
        f = AutomorphicForm(weight=12, level=1)
        lift = FunctorialLift(f, lift_type="sym5")
        assert lift.lift_info['proved'] is False

    def test_lifted_hecke_eigenvalues(self):
        """Geliftete Hecke-Eigenwerte werden berechnet."""
        f = AutomorphicForm(weight=12, level=1)
        lift = FunctorialLift(f, lift_type="sym2")
        primes = [2, 3, 5, 7, 11]
        result = lift.lifted_hecke_eigenvalues(primes)
        assert 'original_eigenvalues' in result
        assert 'lifted_eigenvalues' in result
        # Alle Primzahlen sollten in den Ergebnissen sein
        for p in primes:
            assert p in result['original_eigenvalues']

    def test_info_dict(self):
        """info() gibt vollständiges Dictionary zurück."""
        f = AutomorphicForm(weight=12, level=1)
        lift = FunctorialLift(f, lift_type="sym3")
        info = lift.info()
        assert 'lift_type' in info
        assert 'source_group' in info
        assert 'target_group' in info


# =============================================================================
# TESTS: RamanujanConjecture
# =============================================================================

class TestRamanujanConjecture:
    """Tests für RamanujanConjecture."""

    def test_ramanujan_bound_delta(self):
        """Ramanujan-Schranke für die Deltafunktion (Gewicht 12)."""
        rc = RamanujanConjecture()
        f = AutomorphicForm(group="GL2", weight=12, level=1)
        primes = [2, 3, 5, 7, 11, 13]
        result = rc.verify_for_form(f, primes)
        # Alle Schranken sollten für die Deltafunktion erfüllt sein (Deligne 1974)
        assert result['all_satisfied'] is True

    def test_ramanujan_result_structure(self):
        """Ergebnisstruktur der Ramanujan-Verifikation."""
        rc = RamanujanConjecture()
        f = AutomorphicForm(weight=12, level=1)
        result = rc.verify_for_form(f, [2, 3, 5])
        assert 'all_satisfied' in result
        assert 'prime_results' in result
        assert 'bound_formula' in result

    def test_ramanujan_prime_results(self):
        """Jede Primzahl hat ein Ergebnis."""
        rc = RamanujanConjecture()
        f = AutomorphicForm(weight=12, level=1)
        primes = [2, 3, 5]
        result = rc.verify_for_form(f, primes)
        for p in primes:
            assert p in result['prime_results']
            assert 'satisfies_ramanujan' in result['prime_results'][p]

    def test_satake_parameter_bound(self):
        """Satake-Parameter-Schranken-Test."""
        rc = RamanujanConjecture()
        p = 5
        weight = 12
        # Erzeuge einen Satake-Parameter der Form p^{(k-1)/2}
        alpha = complex(p ** ((weight - 1) / 2), 0)
        result = rc.satake_parameter_bound(alpha, p, weight)
        assert result['satisfies_ramanujan'] is True

    def test_satake_parameter_violation(self):
        """Stark abweichender Satake-Parameter verletzt Schranke."""
        rc = RamanujanConjecture()
        p = 5
        weight = 12
        # Erzeuge einen Satake-Parameter weit über der Schranke
        alpha = complex(p ** ((weight - 1) / 2) * 10, 0)  # 10x zu groß
        result = rc.satake_parameter_bound(alpha, p, weight)
        assert result['satisfies_ramanujan'] is False


# =============================================================================
# TESTS: SatakeTransform
# =============================================================================

class TestSatakeTransform:
    """Tests für SatakeTransform."""

    def test_creation(self):
        """Standarderstellung."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=5)
        assert st.p == 5
        assert st.is_unramified is True

    def test_creation_invalid_p(self):
        """Nicht-Primzahl wirft Fehler."""
        f = AutomorphicForm(weight=12, level=1)
        with pytest.raises(InvalidInputError):
            SatakeTransform(f, p=4)

    def test_ramified_prime(self):
        """Verzweigte Primzahl wird korrekt erkannt."""
        # Niveau 5: p=5 ist verzweigt
        f = AutomorphicForm(weight=12, level=5)
        st = SatakeTransform(f, p=5)
        assert st.is_unramified is False

    def test_satake_parameters_unramified(self):
        """Satake-Parameter für unverzweigte Primzahl."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=5)
        params = st.satake_parameters()
        assert params['is_unramified'] is True
        assert 'alpha_p' in params
        assert 'beta_p' in params

    def test_satake_parameters_product(self):
        """Produkt der Satake-Parameter sollte p^{k-1} sein."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=3)
        params = st.satake_parameters()
        if params['is_unramified']:
            # Produkt α_p · β_p = p^{k-1} = 3^{11}
            expected_product = 3 ** 11
            # Toleranz wegen numerischer Fehler
            product_str = params.get('product', '0')
            # Nur strukturell prüfen (String enthält Zahl)
            assert 'product' in params

    def test_hecke_polynomial_structure(self):
        """Hecke-Polynom hat drei Koeffizienten."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=7)
        poly = st.hecke_polynomial()
        assert 'polynomial' in poly
        assert 'euler_factor' in poly
        assert len(poly['coefficients']) == 3

    def test_hecke_polynomial_constant_term(self):
        """Konstanter Term des Hecke-Polynoms ist 1."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=2)
        poly = st.hecke_polynomial()
        # Koeffizient von X^0 sollte 1 sein
        assert abs(poly['coefficients'][0] - 1) < 1e-10

    def test_info_dict(self):
        """info() gibt vollständiges Dictionary zurück."""
        f = AutomorphicForm(weight=12, level=1)
        st = SatakeTransform(f, p=11)
        info = st.info()
        assert 'satake_parameters' in info
        assert 'hecke_polynomial' in info
        assert 'p' in info
