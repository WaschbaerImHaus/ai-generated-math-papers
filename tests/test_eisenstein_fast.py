"""
@file test_eisenstein_fast.py
@brief Tests für die optimierte Eisenstein-Reihe eisenstein_series_fast().
@description
    Testet die Korrektheit und Konsistenz von eisenstein_series_fast()
    im Vergleich zur Referenzimplementierung eisenstein_series():

    - Übereinstimmung für gerades k (k=4, k=6, k=8, k=10)
    - Fallback für ungerades k auf die Standardimplementierung
    - Verschiedene Punkte τ in der oberen Halbebene
    - Edge-Cases: τ nahe dem Rand, große Imaginärteile

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import pytest
import cmath

# Sicherstellen dass src/ im Python-Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from modular_forms import eisenstein_series, eisenstein_series_fast


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def complex_close(a: complex, b: complex, atol: float = 1e-6) -> bool:
    """
    Prüft ob zwei komplexe Zahlen nahe beieinander liegen.
    Vergleicht Realteil und Imaginärteil separat mit gegebener Toleranz.
    """
    return (abs(a.real - b.real) < atol) and (abs(a.imag - b.imag) < atol)


# =============================================================================
# TESTS: Haupttest für eisenstein_series_fast
# =============================================================================

class TestEisensteinSeriesFast:
    """Testet eisenstein_series_fast() gegen die Referenzimplementierung."""

    # Standardpunkt in der oberen Halbebene für die meisten Tests
    z_standard = 0.5 + 1.5j

    def test_e4_matches_reference_standard_point(self):
        """G_4 an z=0.5+1.5j: fast-Version muss der Referenz entsprechen (Abweichung < 1e-6)."""
        z = self.z_standard
        ref = eisenstein_series(4, z)
        fast = eisenstein_series_fast(4, z)

        # Betrag der Abweichung prüfen
        assert abs(fast - ref) < 1e-6, (
            f"G_4({z}): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e6_matches_reference_standard_point(self):
        """G_6 an z=0.5+1.5j: fast-Version muss der Referenz entsprechen (Abweichung < 1e-6)."""
        z = self.z_standard
        ref = eisenstein_series(6, z)
        fast = eisenstein_series_fast(6, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_6({z}): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e4_matches_reference_imaginary_unit(self):
        """G_4 am Punkt τ=i: fast-Version muss der Referenz entsprechen."""
        z = 1j  # τ = i ist der Standardpunkt der oberen Halbebene
        ref = eisenstein_series(4, z)
        fast = eisenstein_series_fast(4, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_4(i): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e6_matches_reference_imaginary_unit(self):
        """G_6 am Punkt τ=i: fast-Version muss der Referenz entsprechen."""
        z = 1j
        ref = eisenstein_series(6, z)
        fast = eisenstein_series_fast(6, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_6(i): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e8_matches_reference(self):
        """G_8 an z=0.5+1.5j: fast-Version muss der Referenz entsprechen."""
        z = self.z_standard
        ref = eisenstein_series(8, z)
        fast = eisenstein_series_fast(8, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_8({z}): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e10_matches_reference(self):
        """G_10 an z=0.5+1.5j: fast-Version muss der Referenz entsprechen."""
        z = self.z_standard
        ref = eisenstein_series(10, z)
        fast = eisenstein_series_fast(10, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_10({z}): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e4_at_pure_imaginary(self):
        """G_4 an τ=2i (rein imaginär): Vergleich mit Referenz."""
        z = 2j
        ref = eisenstein_series(4, z)
        fast = eisenstein_series_fast(4, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_4(2i): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e6_at_pure_imaginary(self):
        """G_6 an τ=2i: Vergleich mit Referenz."""
        z = 2j
        ref = eisenstein_series(6, z)
        fast = eisenstein_series_fast(6, z)

        assert abs(fast - ref) < 1e-6, (
            f"G_6(2i): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_e4_at_rho(self):
        """G_4 am Fixpunkt ρ = e^(2πi/3) = -1/2 + i√3/2: Vergleich mit Referenz."""
        import math
        # ρ ist ein besonderer Punkt: E_4(ρ) = 0 für die normierte Eisenstein-Reihe
        z = complex(-0.5, math.sqrt(3) / 2)
        ref = eisenstein_series(4, z)
        fast = eisenstein_series_fast(4, z)

        assert abs(fast - ref) < 1e-5, (
            f"G_4(ρ): Referenz={ref}, Fast={fast}, Abweichung={abs(fast - ref)}"
        )

    def test_real_part_sign_consistency(self):
        """Realteil von eisenstein_series_fast muss dasselbe Vorzeichen haben wie Referenz."""
        z = 0.3 + 1.2j
        ref = eisenstein_series(4, z)
        fast = eisenstein_series_fast(4, z)

        # Vorzeichen des Realteils muss übereinstimmen (grobe Konsistenzprüfung)
        if abs(ref.real) > 1e-10:  # Nur wenn Realteil signifikant
            assert (ref.real > 0) == (fast.real > 0), (
                f"Vorzeichenunterschied im Realteil: ref={ref.real}, fast={fast.real}"
            )

    def test_imaginary_part_sign_consistency(self):
        """Imaginärteil von eisenstein_series_fast muss dasselbe Vorzeichen haben wie Referenz."""
        z = 0.7 + 0.9j
        ref = eisenstein_series(6, z)
        fast = eisenstein_series_fast(6, z)

        if abs(ref.imag) > 1e-10:
            assert (ref.imag > 0) == (fast.imag > 0), (
                f"Vorzeichenunterschied im Imaginärteil: ref={ref.imag}, fast={fast.imag}"
            )


# =============================================================================
# TESTS: Fallback-Verhalten
# =============================================================================

class TestEisensteinSeriesFastFallback:
    """Testet Fallback-Verhalten für ungerades k und k < 4."""

    z_test = 0.5 + 1.5j

    def test_k3_falls_back_to_standard(self):
        """Für k=3 (ungerade) muss eisenstein_series_fast auf Referenz zurückfallen."""
        # k=3 ist ungerade → Symmetrie gilt nicht → Fallback
        # ABER: eisenstein_series() wirft ValueError für k < 4
        with pytest.raises(ValueError):
            eisenstein_series_fast(3, self.z_test)

    def test_k2_falls_back_to_standard(self):
        """Für k=2 muss eisenstein_series_fast auf Referenz zurückfallen (k < 4)."""
        # k=2 ist zwar gerade, aber k < 4 → kein Fallback, Referenz wirft ValueError
        with pytest.raises(ValueError):
            eisenstein_series_fast(2, self.z_test)

    def test_k5_falls_back(self):
        """Für k=5 (ungerade) muss eisenstein_series_fast den Fallback nutzen."""
        # Ungerades k → Fallback → gleiche Exception wie eisenstein_series()
        with pytest.raises(ValueError):
            eisenstein_series_fast(5, self.z_test)

    def test_negative_imaginary_part_raises(self):
        """τ mit Im(τ) ≤ 0 muss eine ValueError auslösen."""
        z_invalid = 0.5 - 1.0j  # Im(τ) < 0: nicht in oberer Halbebene
        with pytest.raises(ValueError):
            eisenstein_series_fast(4, z_invalid)


# =============================================================================
# TESTS: Numerische Konsistenz
# =============================================================================

class TestEisensteinSeriesFastNumerical:
    """Testet numerische Eigenschaften von eisenstein_series_fast."""

    def test_returns_complex_type(self):
        """eisenstein_series_fast muss eine komplexe Zahl zurückgeben."""
        z = 0.5 + 1.5j
        result = eisenstein_series_fast(4, z)
        assert isinstance(result, complex), f"Erwartet complex, erhalten: {type(result)}"

    def test_result_is_finite(self):
        """eisenstein_series_fast darf keine NaN oder Inf-Werte zurückgeben."""
        z = 0.5 + 1.5j
        result = eisenstein_series_fast(4, z)
        # Realteil und Imaginärteil müssen endlich sein
        assert cmath.isfinite(result), f"Ergebnis ist nicht endlich: {result}"

    def test_multiple_points_finite(self):
        """eisenstein_series_fast muss für verschiedene τ-Werte endliche Ergebnisse liefern."""
        import math
        test_points = [
            0.5 + 1.0j,
            0.0 + 1.0j,    # rein imaginär
            0.0 + 2.0j,    # größerer Imaginärteil
            0.5 + 0.5j,    # kleinerer Imaginärteil (nahe Rand)
            -0.4 + 1.1j,   # negativer Realteil
        ]

        for z in test_points:
            result = eisenstein_series_fast(4, z)
            assert cmath.isfinite(result), (
                f"Nicht-endliches Ergebnis bei τ={z}: {result}"
            )

    def test_tolerance_k4_multiple_points(self):
        """G_4: Toleranztest für mehrere Punkte gleichzeitig."""
        test_points = [
            0.5 + 1.5j,
            0.1 + 1.0j,
            0.0 + 3.0j,
        ]
        for z in test_points:
            ref = eisenstein_series(4, z)
            fast = eisenstein_series_fast(4, z)
            assert abs(fast - ref) < 1e-6, (
                f"G_4({z}): zu große Abweichung {abs(fast - ref):.2e}"
            )

    def test_tolerance_k6_multiple_points(self):
        """G_6: Toleranztest für mehrere Punkte gleichzeitig."""
        test_points = [
            0.5 + 1.5j,
            0.3 + 2.0j,
            0.0 + 1.5j,
        ]
        for z in test_points:
            ref = eisenstein_series(6, z)
            fast = eisenstein_series_fast(6, z)
            assert abs(fast - ref) < 1e-6, (
                f"G_6({z}): zu große Abweichung {abs(fast - ref):.2e}"
            )

    def test_n_terms_effect(self):
        """Mehr Terme müssen zu besserer Übereinstimmung führen."""
        z = 0.5 + 1.5j
        # Mit weniger Termen → mehr Abweichung
        fast_10 = eisenstein_series_fast(4, z, n_terms=10)
        fast_50 = eisenstein_series_fast(4, z, n_terms=50)
        ref = eisenstein_series(4, z, n_terms=50)

        # n_terms=50 sollte näher an der Referenz (ebenfalls n_terms=50) sein
        error_10 = abs(fast_10 - ref)
        error_50 = abs(fast_50 - ref)
        # fast_50 sollte nicht schlechter als fast_10 sein
        # (bei identischer n_terms in Referenz: error_50 ≈ 0)
        assert error_50 <= error_10 + 1e-10, (
            f"Mehr Terme sollten weniger Fehler haben: error_10={error_10:.2e}, error_50={error_50:.2e}"
        )
