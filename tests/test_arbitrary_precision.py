"""
@file test_arbitrary_precision.py
@brief Tests für das Arbitrary-Precision-Modul (arbitrary_precision.py).
@description
    Testfälle für alle Funktionen der beliebig genauen Arithmetik:
    - set_precision: Präzision setzen
    - riemann_zeta_highprec: ζ-Funktion mit bekannten Werten
    - riemann_zero_verify: Riemann-Nullstellen auf kritischer Linie
    - pi_highprec: π mit bekannten Dezimalstellen
    - e_highprec: e (Euler-Zahl) mit bekannten Dezimalstellen
    - gamma_highprec: Gamma-Funktion an bekannten Punkten
    - bernoulli_highprec: Bernoulli-Zahlen mit bekannten Werten
    - continued_fraction_expansion: Kettenbrüche für bekannte Zahlen
    - Hilfsfunktionen: euler_mascheroni, log, sqrt

    Mindestens 15 Tests abgedeckt.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest

# Pfad zum Quellverzeichnis hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import mpmath

from arbitrary_precision import (
    set_precision,
    riemann_zeta_highprec,
    riemann_zero_verify,
    pi_highprec,
    e_highprec,
    gamma_highprec,
    bernoulli_highprec,
    continued_fraction_expansion,
    euler_mascheroni_highprec,
    log_highprec,
    sqrt_highprec
)


# =============================================================================
# TESTS: Präzision
# =============================================================================

class TestPrecision:
    """Tests für die Präzisionseinstellung."""

    def test_set_precision(self):
        """Präzision wird korrekt gesetzt."""
        set_precision(50)
        assert mpmath.mp.dps == 50

    def test_set_precision_high(self):
        """Hohe Präzision wird unterstützt."""
        set_precision(500)
        assert mpmath.mp.dps == 500

    def test_set_precision_low(self):
        """Niedrige Präzision wird gesetzt."""
        set_precision(10)
        assert mpmath.mp.dps == 10


# =============================================================================
# TESTS: Riemann-Zeta-Funktion
# =============================================================================

class TestRiemannZeta:
    """Tests für die Riemann-Zeta-Funktion mit hoher Präzision."""

    def test_zeta_2(self):
        """ζ(2) = π²/6 ≈ 1.6449340668..."""
        result = riemann_zeta_highprec(2, digits=30)
        expected = math.pi ** 2 / 6
        assert abs(result.real - expected) < 1e-10

    def test_zeta_4(self):
        """ζ(4) = π⁴/90 ≈ 1.0823232337..."""
        result = riemann_zeta_highprec(4, digits=30)
        expected = math.pi ** 4 / 90
        assert abs(result.real - expected) < 1e-10

    def test_zeta_negative_one(self):
        """ζ(-1) = -1/12 (Ramanujan-Regularisierung)."""
        result = riemann_zeta_highprec(-1, digits=30)
        expected = -1 / 12
        assert abs(result.real - expected) < 1e-10

    def test_zeta_zero(self):
        """ζ(0) = -1/2."""
        result = riemann_zeta_highprec(0, digits=30)
        assert abs(result.real - (-0.5)) < 1e-10

    def test_zeta_complex(self):
        """ζ(s) für komplexes s ergibt komplexen Wert."""
        s = complex(2, 1)
        result = riemann_zeta_highprec(s, digits=20)
        assert isinstance(result, complex)

    def test_zeta_high_precision(self):
        """ζ(2) mit 50 Stellen ist nahe π²/6."""
        result = riemann_zeta_highprec(2, digits=50)
        expected = math.pi ** 2 / 6
        assert abs(result.real - expected) < 1e-12


# =============================================================================
# TESTS: Riemann-Nullstellen
# =============================================================================

class TestRiemannZeroVerify:
    """Tests für die Verifikation der Riemann-Nullstellen."""

    def test_first_zero_on_critical_line(self):
        """Erste Riemann-Nullstelle liegt auf der kritischen Linie."""
        result = riemann_zero_verify(1, digits=30)
        assert result['on_critical_line'] is True

    def test_first_zero_real_part(self):
        """Realteil der ersten Nullstelle ist 1/2."""
        result = riemann_zero_verify(1, digits=30)
        assert abs(result['real_part'] - 0.5) < 1e-10

    def test_first_zero_imaginary_part(self):
        """Imaginärteil der ersten Nullstelle ist ≈ 14.1347..."""
        result = riemann_zero_verify(1, digits=30)
        # Extrahiere Imaginärteil aus String
        im = float(result['imag_part'].replace('e', 'E').split('±')[0])
        assert abs(im - 14.1347) < 0.001

    def test_second_zero_on_critical_line(self):
        """Zweite Riemann-Nullstelle liegt auf der kritischen Linie."""
        result = riemann_zero_verify(2, digits=30)
        assert result['on_critical_line'] is True

    def test_zero_dict_structure(self):
        """Rückgabe-Dict hat die erwartete Struktur."""
        result = riemann_zero_verify(1, digits=20)
        assert 'n' in result
        assert 'real_part' in result
        assert 'imag_part' in result
        assert 'on_critical_line' in result
        assert result['n'] == 1

    def test_zero_invalid_index(self):
        """Ungültiger Index wirft ValueError."""
        with pytest.raises(ValueError):
            riemann_zero_verify(0, digits=20)


# =============================================================================
# TESTS: Konstanten (π, e)
# =============================================================================

class TestConstants:
    """Tests für hochpräzise Konstanten."""

    def test_pi_starts_with_3(self):
        """π beginnt mit 3."""
        result = pi_highprec(digits=50)
        assert result.startswith("3.")

    def test_pi_known_digits(self):
        """π stimmt mit bekannten Stellen überein."""
        result = pi_highprec(digits=20)
        pi_str = str(math.pi)
        # Erste 10 Stellen vergleichen
        assert result[:8] == pi_str[:8]

    def test_e_starts_with_2(self):
        """e beginnt mit 2."""
        result = e_highprec(digits=50)
        assert result.startswith("2.")

    def test_e_known_digits(self):
        """e stimmt mit bekannten Stellen überein."""
        result = e_highprec(digits=20)
        e_str = str(math.e)
        # Erste 8 Stellen vergleichen
        assert result[:8] == e_str[:8]

    def test_pi_length(self):
        """π-String hat mindestens die gewünschte Länge."""
        result = pi_highprec(digits=50)
        # Mindestens 50 Zeichen (inkl. "3.")
        assert len(result.replace('.', '')) >= 50

    def test_euler_mascheroni(self):
        """Euler-Mascheroni-Konstante ≈ 0.5772..."""
        result = euler_mascheroni_highprec(digits=20)
        gamma_val = float(result)
        assert abs(gamma_val - 0.5772156649) < 1e-8


# =============================================================================
# TESTS: Gamma-Funktion
# =============================================================================

class TestGammaFunction:
    """Tests für die Gamma-Funktion mit hoher Präzision."""

    def test_gamma_integer(self):
        """Γ(n) = (n-1)! für positive ganze Zahlen."""
        # Γ(5) = 4! = 24
        result = gamma_highprec(5, digits=30)
        assert abs(result.real - 24.0) < 1e-10

    def test_gamma_one_half(self):
        """Γ(1/2) = √π."""
        result = gamma_highprec(0.5, digits=30)
        expected = math.sqrt(math.pi)
        assert abs(result.real - expected) < 1e-10

    def test_gamma_one(self):
        """Γ(1) = 1."""
        result = gamma_highprec(1, digits=30)
        assert abs(result.real - 1.0) < 1e-10

    def test_gamma_complex(self):
        """Γ(s) für komplexes s."""
        result = gamma_highprec(complex(1, 1), digits=20)
        assert isinstance(result, complex)


# =============================================================================
# TESTS: Bernoulli-Zahlen
# =============================================================================

class TestBernoulliNumbers:
    """Tests für Bernoulli-Zahlen mit hoher Präzision."""

    def test_bernoulli_0(self):
        """B_0 = 1."""
        result = float(bernoulli_highprec(0, digits=20))
        assert abs(result - 1.0) < 1e-10

    def test_bernoulli_1(self):
        """B_1 = -1/2."""
        result = float(bernoulli_highprec(1, digits=20))
        assert abs(result - (-0.5)) < 1e-10

    def test_bernoulli_2(self):
        """B_2 = 1/6."""
        result = float(bernoulli_highprec(2, digits=20))
        assert abs(result - (1.0 / 6)) < 1e-10

    def test_bernoulli_odd_zero(self):
        """B_n = 0 für ungerade n > 1."""
        result = float(bernoulli_highprec(3, digits=20))
        assert abs(result) < 1e-10

        result5 = float(bernoulli_highprec(5, digits=20))
        assert abs(result5) < 1e-10

    def test_bernoulli_negative_index(self):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            bernoulli_highprec(-1, digits=20)

    def test_bernoulli_4(self):
        """B_4 = -1/30."""
        result = float(bernoulli_highprec(4, digits=20))
        assert abs(result - (-1.0 / 30)) < 1e-10


# =============================================================================
# TESTS: Kettenbruchentwicklung
# =============================================================================

class TestContinuedFraction:
    """Tests für die Kettenbruchentwicklung."""

    def test_cf_integer(self):
        """Ganze Zahl hat Kettenbruch [n]."""
        result = continued_fraction_expansion(5.0, n_terms=10)
        assert result[0] == 5

    def test_cf_pi_starts(self):
        """π = [3; 7, 15, 1, 292, ...]."""
        result = continued_fraction_expansion(math.pi, n_terms=5)
        assert result[0] == 3
        assert result[1] == 7

    def test_cf_e_starts(self):
        """e = [2; 1, 2, 1, 1, 4, ...]."""
        result = continued_fraction_expansion(math.e, n_terms=6)
        assert result[0] == 2
        assert result[1] == 1

    def test_cf_sqrt2_starts(self):
        """√2 = [1; 2, 2, 2, ...] (periodisch)."""
        result = continued_fraction_expansion(math.sqrt(2), n_terms=5)
        assert result[0] == 1
        # Die restlichen Koeffizienten sollten 2 sein
        for coef in result[1:4]:
            assert coef == 2

    def test_cf_half(self):
        """1/2 = [0; 2] als Kettenbruch."""
        result = continued_fraction_expansion(0.5, n_terms=5)
        assert result[0] == 0
        assert result[1] == 2


# =============================================================================
# TESTS: Hilfsfunktionen
# =============================================================================

class TestHelperFunctions:
    """Tests für weitere Hilfsfunktionen."""

    def test_log_natural(self):
        """Natürlicher Logarithmus von e ist 1."""
        result = float(log_highprec(math.e, digits=20))
        assert abs(result - 1.0) < 1e-10

    def test_log_base10(self):
        """log₁₀(100) = 2."""
        result = float(log_highprec(100, base=10, digits=20))
        assert abs(result - 2.0) < 1e-10

    def test_log_invalid(self):
        """Logarithmus für x ≤ 0 wirft ValueError."""
        with pytest.raises(ValueError):
            log_highprec(-1, digits=20)

    def test_sqrt_four(self):
        """√4 = 2."""
        result = float(sqrt_highprec(4, digits=20))
        assert abs(result - 2.0) < 1e-10

    def test_sqrt_two(self):
        """√2 ≈ 1.41421356..."""
        result = float(sqrt_highprec(2, digits=20))
        assert abs(result - math.sqrt(2)) < 1e-10

    def test_sqrt_invalid(self):
        """√x für x < 0 wirft ValueError."""
        with pytest.raises(ValueError):
            sqrt_highprec(-1, digits=20)


# ---------------------------------------------------------------------------
# NEUE FUNKTIONEN (Build 61): zeta_zeros_mpmath, verify_riemann_hypothesis_mpmath, pi_mpmath
# ---------------------------------------------------------------------------

class TestZetaZerosMpmath:
    """Tests für zeta_zeros_mpmath() (Build 61)."""

    def test_first_zero_imaginary_part(self):
        """Die erste Riemann-Nullstelle hat Im(s₁) ≈ 14.134..."""
        from arbitrary_precision import zeta_zeros_mpmath
        zeros = zeta_zeros_mpmath(1, prec=20)
        assert len(zeros) == 1
        import mpmath
        t1 = float(mpmath.im(zeros[0]))
        assert abs(t1 - 14.1347) < 0.001

    def test_first_zero_real_part(self):
        """Die erste Riemann-Nullstelle hat Re(s₁) = 0.5."""
        from arbitrary_precision import zeta_zeros_mpmath
        zeros = zeta_zeros_mpmath(1, prec=20)
        import mpmath
        sigma1 = float(mpmath.re(zeros[0]))
        assert abs(sigma1 - 0.5) < 1e-6

    def test_three_zeros_count(self):
        """zeta_zeros_mpmath(3) gibt 3 Nullstellen zurück."""
        from arbitrary_precision import zeta_zeros_mpmath
        zeros = zeta_zeros_mpmath(3, prec=15)
        assert len(zeros) == 3

    def test_invalid_n(self):
        """n < 1 wirft ValueError."""
        from arbitrary_precision import zeta_zeros_mpmath
        with pytest.raises(ValueError):
            zeta_zeros_mpmath(0, prec=10)

    def test_zeros_increasing_imaginary(self):
        """Imaginärteile der Nullstellen sind aufsteigend."""
        from arbitrary_precision import zeta_zeros_mpmath
        import mpmath
        zeros = zeta_zeros_mpmath(5, prec=15)
        imag_parts = [float(mpmath.im(z)) for z in zeros]
        assert imag_parts == sorted(imag_parts)


class TestVerifyRiemannHypothesisMpmath:
    """Tests für verify_riemann_hypothesis_mpmath() (Build 61)."""

    def test_verify_first_5_zeros(self):
        """Die ersten 5 Nullstellen liegen auf Re(s) = 1/2."""
        from arbitrary_precision import verify_riemann_hypothesis_mpmath
        result = verify_riemann_hypothesis_mpmath(5, prec=30)
        assert result['n_checked'] == 5
        assert result['all_on_critical_line'] is True
        assert result['max_deviation'] < 1e-10

    def test_result_keys(self):
        """Ergebnis enthält alle erwarteten Schlüssel."""
        from arbitrary_precision import verify_riemann_hypothesis_mpmath
        result = verify_riemann_hypothesis_mpmath(3, prec=20)
        assert 'n_checked' in result
        assert 'all_on_critical_line' in result
        assert 'max_deviation' in result
        assert 'zeros' in result
        assert 'real_parts' in result

    def test_real_parts_are_half(self):
        """Alle Realteile der ersten 3 Nullstellen sind ≈ 0.5."""
        from arbitrary_precision import verify_riemann_hypothesis_mpmath
        result = verify_riemann_hypothesis_mpmath(3, prec=30)
        for rp in result['real_parts']:
            assert abs(rp - 0.5) < 1e-8

    def test_invalid_n_zeros(self):
        """n_zeros < 1 wirft ValueError."""
        from arbitrary_precision import verify_riemann_hypothesis_mpmath
        with pytest.raises(ValueError):
            verify_riemann_hypothesis_mpmath(0, prec=20)


class TestPiMpmath:
    """Tests für pi_mpmath() (Build 61)."""

    def test_pi_starts_with_3(self):
        """π beginnt mit 3."""
        from arbitrary_precision import pi_mpmath
        pi_str = pi_mpmath(prec=10)
        assert pi_str.startswith("3.")

    def test_pi_known_digits(self):
        """π stimmt mit bekannten Stellen überein."""
        from arbitrary_precision import pi_mpmath
        pi_str = pi_mpmath(prec=15)
        # π = 3.14159265358979...
        pi_float = float(pi_str)
        assert abs(pi_float - math.pi) < 1e-12

    def test_pi_high_precision(self):
        """π mit 50 Stellen hat korrekte Länge."""
        from arbitrary_precision import pi_mpmath
        pi_str = pi_mpmath(prec=50)
        # String sollte mindestens 50 Zeichen haben
        assert len(pi_str) >= 50

    def test_pi_invalid_prec(self):
        """prec < 1 wirft ValueError."""
        from arbitrary_precision import pi_mpmath
        with pytest.raises(ValueError):
            pi_mpmath(prec=0)
