"""
@file test_config.py
@brief Unit-Tests für das Konfigurationsmodul config.py.
@description
    Überprüft alle globalen Konstanten in config.py auf:
    - Vorhandensein (AttributeError wenn fehlend)
    - Korrekten Datentyp
    - Sinnvolle Werte (Plausibilitätsprüfungen)
    - Mathematische Korrektheit (GOLDEN_RATIO, INV_GOLDEN_RATIO)

    Teststruktur (Test-Driven Development):
    - TestEpsilonAndStepSizes: Maschinengenauigkeit und Schrittweiten
    - TestIterationParameters: MAX_ITERATIONS, Toleranzen
    - TestAlgorithmParameters: SIMPSON_N, TAYLOR_DEFAULT_TERMS, Zeta-Parameter
    - TestMillerRabinWitnesses: Primzahl-Zeugen
    - TestLanczosCoefficients: Gamma-Funktion Koeffizienten
    - TestMathematicalConstants: Goldener Schnitt
    - TestMetadata: Autor, Version, MPMATH-Präzision

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math
import pytest

# Projektverzeichnis zum Suchpfad hinzufügen, damit src importierbar ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import config


# ===========================================================================
# TESTS: MASCHINENGENAUIGKEIT UND SCHRITTWEITEN
# ===========================================================================

class TestEpsilonAndStepSizes:
    """Testet EPSILON, H_DERIVATIVE_1, H_DERIVATIVE_2."""

    def test_epsilon_exists(self):
        """EPSILON muss als Attribut in config vorhanden sein."""
        assert hasattr(config, 'EPSILON'), "EPSILON fehlt in config.py"

    def test_epsilon_is_float(self):
        """EPSILON muss ein float sein."""
        assert isinstance(config.EPSILON, float), \
            f"EPSILON sollte float sein, ist aber {type(config.EPSILON)}"

    def test_epsilon_smaller_than_1e14(self):
        """EPSILON muss kleiner als 1e-14 sein (Maschinengenauigkeit-Anforderung)."""
        assert config.EPSILON < 1e-14, \
            f"EPSILON={config.EPSILON} ist nicht kleiner als 1e-14"

    def test_epsilon_positive(self):
        """EPSILON muss positiv sein."""
        assert config.EPSILON > 0, "EPSILON muss eine positive Zahl sein"

    def test_h_derivative_1_exists(self):
        """H_DERIVATIVE_1 muss vorhanden sein."""
        assert hasattr(config, 'H_DERIVATIVE_1'), "H_DERIVATIVE_1 fehlt in config.py"

    def test_h_derivative_1_is_float(self):
        """H_DERIVATIVE_1 muss ein float sein."""
        assert isinstance(config.H_DERIVATIVE_1, float)

    def test_h_derivative_1_reasonable_range(self):
        """H_DERIVATIVE_1 sollte im Bereich [1e-8, 1e-3] liegen (optimaler Bereich)."""
        assert 1e-8 <= config.H_DERIVATIVE_1 <= 1e-3, \
            f"H_DERIVATIVE_1={config.H_DERIVATIVE_1} liegt außerhalb des sinnvollen Bereichs"

    def test_h_derivative_1_approx_eps_third(self):
        """H_DERIVATIVE_1 ≈ EPSILON^(1/3) laut numerischer Theorie."""
        expected = config.EPSILON ** (1.0 / 3.0)
        # Toleranz 20%: verschiedene EPSILON-Werte möglich
        assert abs(config.H_DERIVATIVE_1 - expected) / expected < 0.2, \
            f"H_DERIVATIVE_1={config.H_DERIVATIVE_1} weicht stark von EPSILON^(1/3)={expected} ab"

    def test_h_derivative_2_exists(self):
        """H_DERIVATIVE_2 muss vorhanden sein."""
        assert hasattr(config, 'H_DERIVATIVE_2'), "H_DERIVATIVE_2 fehlt in config.py"

    def test_h_derivative_2_is_float(self):
        """H_DERIVATIVE_2 muss ein float sein."""
        assert isinstance(config.H_DERIVATIVE_2, float)

    def test_h_derivative_2_larger_than_h1(self):
        """H_DERIVATIVE_2 > H_DERIVATIVE_1: 2. Ableitung braucht größeres h."""
        assert config.H_DERIVATIVE_2 > config.H_DERIVATIVE_1, \
            "H_DERIVATIVE_2 sollte größer als H_DERIVATIVE_1 sein"

    def test_h_derivative_2_approx_eps_fourth(self):
        """H_DERIVATIVE_2 ≈ EPSILON^(1/4) laut numerischer Theorie."""
        expected = config.EPSILON ** (1.0 / 4.0)
        # Toleranz 20%
        assert abs(config.H_DERIVATIVE_2 - expected) / expected < 0.2, \
            f"H_DERIVATIVE_2={config.H_DERIVATIVE_2} weicht stark von EPSILON^(1/4)={expected} ab"


# ===========================================================================
# TESTS: ITERATIONSPARAMETER
# ===========================================================================

class TestIterationParameters:
    """Testet MAX_ITERATIONS, NEWTON_TOL, BISECTION_TOL."""

    def test_max_iterations_exists(self):
        """MAX_ITERATIONS muss vorhanden sein."""
        assert hasattr(config, 'MAX_ITERATIONS'), "MAX_ITERATIONS fehlt in config.py"

    def test_max_iterations_is_int(self):
        """MAX_ITERATIONS muss ein int sein."""
        assert isinstance(config.MAX_ITERATIONS, int), \
            f"MAX_ITERATIONS sollte int sein, ist aber {type(config.MAX_ITERATIONS)}"

    def test_max_iterations_positive(self):
        """MAX_ITERATIONS muss > 0 sein."""
        assert config.MAX_ITERATIONS > 0, \
            f"MAX_ITERATIONS={config.MAX_ITERATIONS} muss positiv sein"

    def test_max_iterations_at_least_100(self):
        """MAX_ITERATIONS sollte mindestens 100 betragen (genug für Konvergenz)."""
        assert config.MAX_ITERATIONS >= 100, \
            f"MAX_ITERATIONS={config.MAX_ITERATIONS} ist zu klein für gute Konvergenz"

    def test_newton_tol_exists(self):
        """NEWTON_TOL muss vorhanden sein."""
        assert hasattr(config, 'NEWTON_TOL'), "NEWTON_TOL fehlt in config.py"

    def test_newton_tol_is_float(self):
        """NEWTON_TOL muss ein float sein."""
        assert isinstance(config.NEWTON_TOL, float)

    def test_newton_tol_positive(self):
        """NEWTON_TOL muss positiv sein."""
        assert config.NEWTON_TOL > 0

    def test_newton_tol_smaller_than_1e8(self):
        """NEWTON_TOL muss kleiner als 1e-8 sein (hohe Genauigkeit)."""
        assert config.NEWTON_TOL < 1e-8, \
            f"NEWTON_TOL={config.NEWTON_TOL} ist zu groß für genaue Newton-Iteration"

    def test_bisection_tol_exists(self):
        """BISECTION_TOL muss vorhanden sein."""
        assert hasattr(config, 'BISECTION_TOL'), "BISECTION_TOL fehlt in config.py"

    def test_bisection_tol_is_float(self):
        """BISECTION_TOL muss ein float sein."""
        assert isinstance(config.BISECTION_TOL, float)

    def test_bisection_tol_positive(self):
        """BISECTION_TOL muss positiv sein."""
        assert config.BISECTION_TOL > 0

    def test_bisection_tol_smaller_than_1e8(self):
        """BISECTION_TOL muss kleiner als 1e-8 sein."""
        assert config.BISECTION_TOL < 1e-8


# ===========================================================================
# TESTS: INTEGRATIONSPARAMETER
# ===========================================================================

class TestIntegrationParameters:
    """Testet SIMPSON_N."""

    def test_simpson_n_exists(self):
        """SIMPSON_N muss vorhanden sein."""
        assert hasattr(config, 'SIMPSON_N'), "SIMPSON_N fehlt in config.py"

    def test_simpson_n_is_int(self):
        """SIMPSON_N muss ein int sein."""
        assert isinstance(config.SIMPSON_N, int)

    def test_simpson_n_positive(self):
        """SIMPSON_N muss positiv sein."""
        assert config.SIMPSON_N > 0

    def test_simpson_n_even(self):
        """SIMPSON_N muss gerade sein (Anforderung der Simpson-Regel)."""
        assert config.SIMPSON_N % 2 == 0, \
            f"SIMPSON_N={config.SIMPSON_N} muss gerade sein für die Simpson-Regel"

    def test_simpson_n_at_least_100(self):
        """SIMPSON_N sollte mindestens 100 betragen für gute Genauigkeit."""
        assert config.SIMPSON_N >= 100


# ===========================================================================
# TESTS: TAYLOR-REIHEN
# ===========================================================================

class TestTaylorParameters:
    """Testet TAYLOR_DEFAULT_TERMS."""

    def test_taylor_default_terms_exists(self):
        """TAYLOR_DEFAULT_TERMS muss vorhanden sein."""
        assert hasattr(config, 'TAYLOR_DEFAULT_TERMS'), \
            "TAYLOR_DEFAULT_TERMS fehlt in config.py"

    def test_taylor_default_terms_is_int(self):
        """TAYLOR_DEFAULT_TERMS muss ein int sein."""
        assert isinstance(config.TAYLOR_DEFAULT_TERMS, int)

    def test_taylor_default_terms_positive(self):
        """TAYLOR_DEFAULT_TERMS muss positiv sein."""
        assert config.TAYLOR_DEFAULT_TERMS > 0

    def test_taylor_default_terms_at_least_5(self):
        """TAYLOR_DEFAULT_TERMS sollte mindestens 5 betragen."""
        assert config.TAYLOR_DEFAULT_TERMS >= 5


# ===========================================================================
# TESTS: ZETA-FUNKTION PARAMETER
# ===========================================================================

class TestZetaParameters:
    """Testet ZETA_EULER_KNOPP_TERMS und ZETA_EULER_MACLAURIN_TERMS."""

    def test_euler_knopp_terms_exists(self):
        """ZETA_EULER_KNOPP_TERMS muss vorhanden sein."""
        assert hasattr(config, 'ZETA_EULER_KNOPP_TERMS')

    def test_euler_knopp_terms_is_int(self):
        """ZETA_EULER_KNOPP_TERMS muss ein int sein."""
        assert isinstance(config.ZETA_EULER_KNOPP_TERMS, int)

    def test_euler_knopp_terms_at_least_30(self):
        """ZETA_EULER_KNOPP_TERMS sollte mindestens 30 für Maschinengenauigkeit sein."""
        assert config.ZETA_EULER_KNOPP_TERMS >= 30

    def test_euler_maclaurin_terms_exists(self):
        """ZETA_EULER_MACLAURIN_TERMS muss vorhanden sein."""
        assert hasattr(config, 'ZETA_EULER_MACLAURIN_TERMS')

    def test_euler_maclaurin_terms_is_int(self):
        """ZETA_EULER_MACLAURIN_TERMS muss ein int sein."""
        assert isinstance(config.ZETA_EULER_MACLAURIN_TERMS, int)

    def test_euler_maclaurin_terms_positive(self):
        """ZETA_EULER_MACLAURIN_TERMS muss positiv sein."""
        assert config.ZETA_EULER_MACLAURIN_TERMS > 0

    def test_euler_maclaurin_more_than_knopp(self):
        """Euler-Maclaurin benötigt mehr Terme als Euler-Knopp für gleiche Genauigkeit."""
        assert config.ZETA_EULER_MACLAURIN_TERMS > config.ZETA_EULER_KNOPP_TERMS


# ===========================================================================
# TESTS: MILLER-RABIN ZEUGEN
# ===========================================================================

class TestMillerRabinWitnesses:
    """Testet PRIME_MILLER_RABIN_WITNESSES."""

    def test_witnesses_exists(self):
        """PRIME_MILLER_RABIN_WITNESSES muss vorhanden sein."""
        assert hasattr(config, 'PRIME_MILLER_RABIN_WITNESSES')

    def test_witnesses_is_list(self):
        """PRIME_MILLER_RABIN_WITNESSES muss eine Liste sein."""
        assert isinstance(config.PRIME_MILLER_RABIN_WITNESSES, list), \
            "PRIME_MILLER_RABIN_WITNESSES muss eine Liste sein"

    def test_witnesses_length_at_least_4(self):
        """Mindestens 4 Zeugen für deterministischen Test kleiner Zahlen."""
        assert len(config.PRIME_MILLER_RABIN_WITNESSES) >= 4, \
            f"Nur {len(config.PRIME_MILLER_RABIN_WITNESSES)} Zeugen, mindestens 4 erwartet"

    def test_witnesses_all_positive_ints(self):
        """Alle Zeugen müssen positive ganze Zahlen sein."""
        for w in config.PRIME_MILLER_RABIN_WITNESSES:
            assert isinstance(w, int) and w > 0, \
                f"Zeuge {w} ist kein positiver Integer"

    def test_witnesses_all_prime(self):
        """Alle Standard-Zeugen sollten Primzahlen sein (Best Practice)."""
        def is_small_prime(n):
            """Trivialer Primzahltest für kleine Zahlen."""
            if n < 2:
                return False
            if n == 2:
                return True
            if n % 2 == 0:
                return False
            for i in range(3, int(n**0.5) + 1, 2):
                if n % i == 0:
                    return False
            return True

        for w in config.PRIME_MILLER_RABIN_WITNESSES:
            assert is_small_prime(w), \
                f"Zeuge {w} ist keine Primzahl (nicht Standard-Zeuge)"

    def test_witnesses_contains_2(self):
        """2 muss in den Zeugen enthalten sein (wichtigster Zeuge)."""
        assert 2 in config.PRIME_MILLER_RABIN_WITNESSES

    def test_witnesses_contains_3(self):
        """3 muss in den Zeugen enthalten sein."""
        assert 3 in config.PRIME_MILLER_RABIN_WITNESSES

    def test_witnesses_sorted(self):
        """Zeugen sollten aufsteigend sortiert sein."""
        assert config.PRIME_MILLER_RABIN_WITNESSES == sorted(
            config.PRIME_MILLER_RABIN_WITNESSES
        ), "PRIME_MILLER_RABIN_WITNESSES sollte aufsteigend sortiert sein"


# ===========================================================================
# TESTS: LANCZOS-KOEFFIZIENTEN
# ===========================================================================

class TestLanczosCoefficients:
    """Testet LANCZOS_COEFFICIENTS."""

    def test_lanczos_exists(self):
        """LANCZOS_COEFFICIENTS muss vorhanden sein."""
        assert hasattr(config, 'LANCZOS_COEFFICIENTS')

    def test_lanczos_is_dict(self):
        """LANCZOS_COEFFICIENTS muss ein Dict sein."""
        assert isinstance(config.LANCZOS_COEFFICIENTS, dict)

    def test_lanczos_has_g_key(self):
        """LANCZOS_COEFFICIENTS muss den Schlüssel 'g' enthalten."""
        assert 'g' in config.LANCZOS_COEFFICIENTS, \
            "LANCZOS_COEFFICIENTS fehlt Schlüssel 'g'"

    def test_lanczos_has_coefficients_key(self):
        """LANCZOS_COEFFICIENTS muss den Schlüssel 'coefficients' enthalten."""
        assert 'coefficients' in config.LANCZOS_COEFFICIENTS, \
            "LANCZOS_COEFFICIENTS fehlt Schlüssel 'coefficients'"

    def test_lanczos_g_is_7(self):
        """Lanczos-Parameter g sollte 7 sein (Standardwert)."""
        assert config.LANCZOS_COEFFICIENTS['g'] == 7

    def test_lanczos_coefficients_is_list(self):
        """Die Koeffizientenliste muss eine Liste sein."""
        assert isinstance(config.LANCZOS_COEFFICIENTS['coefficients'], list)

    def test_lanczos_coefficients_length(self):
        """Für g=7 werden g+2=9 Koeffizienten benötigt (c_0 bis c_8)."""
        g = config.LANCZOS_COEFFICIENTS['g']
        coeffs = config.LANCZOS_COEFFICIENTS['coefficients']
        assert len(coeffs) == g + 2, \
            f"Für g={g} werden {g+2} Koeffizienten erwartet, {len(coeffs)} gefunden"

    def test_lanczos_first_coefficient_near_1(self):
        """Der erste Koeffizient c_0 ≈ 1.0 (Normierungsbedingung)."""
        c0 = config.LANCZOS_COEFFICIENTS['coefficients'][0]
        assert abs(c0 - 1.0) < 0.01, \
            f"Erster Lanczos-Koeffizient c_0={c0} sollte nahe 1.0 sein"


# ===========================================================================
# TESTS: MPMATH
# ===========================================================================

class TestMpmathConfig:
    """Testet MPMATH_DEFAULT_DPS."""

    def test_mpmath_dps_exists(self):
        """MPMATH_DEFAULT_DPS muss vorhanden sein."""
        assert hasattr(config, 'MPMATH_DEFAULT_DPS')

    def test_mpmath_dps_is_int(self):
        """MPMATH_DEFAULT_DPS muss ein int sein."""
        assert isinstance(config.MPMATH_DEFAULT_DPS, int)

    def test_mpmath_dps_at_least_15(self):
        """MPMATH_DEFAULT_DPS sollte mindestens 15 sein (doppelte Genauigkeit)."""
        assert config.MPMATH_DEFAULT_DPS >= 15

    def test_mpmath_dps_reasonable(self):
        """MPMATH_DEFAULT_DPS sollte ≤ 1000 sein (praktische Obergrenze)."""
        assert config.MPMATH_DEFAULT_DPS <= 1000


# ===========================================================================
# TESTS: MATHEMATISCHE KONSTANTEN
# ===========================================================================

class TestMathematicalConstants:
    """Testet GOLDEN_RATIO und INV_GOLDEN_RATIO."""

    def test_golden_ratio_exists(self):
        """GOLDEN_RATIO muss vorhanden sein."""
        assert hasattr(config, 'GOLDEN_RATIO')

    def test_golden_ratio_is_float(self):
        """GOLDEN_RATIO muss ein float sein."""
        assert isinstance(config.GOLDEN_RATIO, float)

    def test_golden_ratio_value(self):
        """GOLDEN_RATIO ≈ 1.6180339887... (Goldener Schnitt)."""
        expected = (1.0 + math.sqrt(5.0)) / 2.0
        assert abs(config.GOLDEN_RATIO - expected) < 1e-12, \
            f"GOLDEN_RATIO={config.GOLDEN_RATIO} stimmt nicht mit (1+√5)/2={expected} überein"

    def test_golden_ratio_property(self):
        """Goldener Schnitt Eigenschaft: φ² = φ + 1."""
        phi = config.GOLDEN_RATIO
        assert abs(phi**2 - (phi + 1)) < 1e-12, \
            f"φ²={phi**2} ≠ φ+1={phi+1}: Goldener-Schnitt-Eigenschaft verletzt"

    def test_inv_golden_ratio_exists(self):
        """INV_GOLDEN_RATIO muss vorhanden sein."""
        assert hasattr(config, 'INV_GOLDEN_RATIO')

    def test_inv_golden_ratio_is_float(self):
        """INV_GOLDEN_RATIO muss ein float sein."""
        assert isinstance(config.INV_GOLDEN_RATIO, float)

    def test_inv_golden_ratio_value(self):
        """INV_GOLDEN_RATIO = 2/(1+√5) ≈ 0.6180339887..."""
        expected = 2.0 / (1.0 + math.sqrt(5.0))
        assert abs(config.INV_GOLDEN_RATIO - expected) < 1e-12

    def test_golden_ratio_times_inv_is_1(self):
        """φ · (1/φ) = 1 (Reziproke Eigenschaft)."""
        product = config.GOLDEN_RATIO * config.INV_GOLDEN_RATIO
        assert abs(product - 1.0) < 1e-12, \
            f"φ · (1/φ) = {product} ≠ 1"

    def test_inv_golden_ratio_alternative(self):
        """1/φ = φ - 1 (weitere Eigenschaft des Goldenen Schnitts)."""
        assert abs(config.INV_GOLDEN_RATIO - (config.GOLDEN_RATIO - 1.0)) < 1e-12


# ===========================================================================
# TESTS: METADATEN
# ===========================================================================

class TestMetadata:
    """Testet AUTHOR und VERSION."""

    def test_author_exists(self):
        """AUTHOR muss vorhanden sein."""
        assert hasattr(config, 'AUTHOR')

    def test_author_is_string(self):
        """AUTHOR muss ein String sein."""
        assert isinstance(config.AUTHOR, str)

    def test_author_value(self):
        """AUTHOR muss 'Kurt Ingwer' sein."""
        assert config.AUTHOR == "Kurt Ingwer", \
            f"AUTHOR='{config.AUTHOR}' ist nicht 'Kurt Ingwer'"

    def test_version_exists(self):
        """VERSION muss vorhanden sein."""
        assert hasattr(config, 'VERSION')

    def test_version_is_string(self):
        """VERSION muss ein String sein."""
        assert isinstance(config.VERSION, str)

    def test_version_not_empty(self):
        """VERSION darf nicht leer sein."""
        assert len(config.VERSION) > 0

    def test_version_format(self):
        """VERSION sollte im Format 'X.Y' sein (z.B. '11.0')."""
        parts = config.VERSION.split('.')
        assert len(parts) >= 1, "VERSION muss mindestens eine Zahl enthalten"
        # Erster Teil muss eine Ganzzahl sein
        assert parts[0].isdigit(), f"Hauptversion '{parts[0]}' muss eine Zahl sein"


# ===========================================================================
# TESTS: VOLLSTÄNDIGKEIT (Smoke-Tests)
# ===========================================================================

class TestCompleteness:
    """Stellt sicher, dass alle erwarteten Konstanten vorhanden sind."""

    REQUIRED_CONSTANTS = [
        'EPSILON',
        'H_DERIVATIVE_1',
        'H_DERIVATIVE_2',
        'MAX_ITERATIONS',
        'NEWTON_TOL',
        'BISECTION_TOL',
        'SIMPSON_N',
        'PRIME_MILLER_RABIN_WITNESSES',
        'TAYLOR_DEFAULT_TERMS',
        'ZETA_EULER_KNOPP_TERMS',
        'ZETA_EULER_MACLAURIN_TERMS',
        'LANCZOS_COEFFICIENTS',
        'MPMATH_DEFAULT_DPS',
        'GOLDEN_RATIO',
        'INV_GOLDEN_RATIO',
        'AUTHOR',
        'VERSION',
    ]

    @pytest.mark.parametrize("const_name", REQUIRED_CONSTANTS)
    def test_constant_exists(self, const_name):
        """Jede erforderliche Konstante muss in config.py vorhanden sein."""
        assert hasattr(config, const_name), \
            f"Konstante '{const_name}' fehlt in config.py"
