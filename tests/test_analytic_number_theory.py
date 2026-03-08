"""
Tests für das analytische Zahlentheorie-Modul (analytic_number_theory.py).

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from analytic_number_theory import (
    prime_counting_function,
    logarithmic_integral,
    prime_number_theorem_comparison,
    von_mangoldt_function,
    chebyshev_theta,
    chebyshev_psi,
    chebyshev_psi_vs_x,
    is_coprime,
    principal_dirichlet_character,
    dirichlet_l_function_partial,
    dirichlet_prime_counting_in_residue_class,
    is_semiprime,
    chen_decomposition,
    count_omega,
    almost_prime,
    prime_gaps,
    maximal_prime_gap
)


# ===========================================================================
# PRIMZAHLZÄHLFUNKTION TESTS
# ===========================================================================

class TestPrimeCountingFunction:
    """Tests für π(x)."""

    def test_small_known_values(self):
        """Bekannte Werte von π(x) für kleine x."""
        # π(2) = 1, π(3) = 2, π(5) = 3, π(10) = 4
        assert prime_counting_function(2) == 1
        assert prime_counting_function(3) == 2
        assert prime_counting_function(5) == 3
        assert prime_counting_function(10) == 4
        assert prime_counting_function(1) == 0

    def test_pi_100_equals_25(self):
        """π(100) = 25 (bekannter Wert)."""
        assert prime_counting_function(100) == 25

    def test_pi_1000_equals_168(self):
        """π(1000) = 168 (bekannter Wert)."""
        assert prime_counting_function(1000) == 168

    def test_pi_10000_equals_1229(self):
        """π(10000) = 1229 (bekannter Wert)."""
        assert prime_counting_function(10000) == 1229

    def test_pi_monotone(self):
        """π(x) ist monoton nicht-fallend."""
        prev = 0
        for x in range(1, 100):
            curr = prime_counting_function(x)
            assert curr >= prev, f"π nicht monoton bei x={x}"
            prev = curr


# ===========================================================================
# LOGARITHMISCHES INTEGRAL TESTS
# ===========================================================================

class TestLogarithmicIntegral:
    """Tests für Li(x)."""

    def test_li_better_than_simple_approximation(self):
        """Li(x) approximiert π(x) besser als x/ln(x) für x ≥ 1000.

        Hinweis: Für kleine x (z.B. x=100) kann x/ln(x) zufällig genauer sein.
        Die Überlegenheit von Li(x) ist asymptotisch und ab ~x=1000 stabil.
        """
        for x in [1000, 10000]:
            pi_x = prime_counting_function(x)
            li_x = logarithmic_integral(x)
            simple = x / math.log(x)
            error_li = abs(li_x - pi_x)
            error_simple = abs(simple - pi_x)
            assert error_li < error_simple, (
                f"Li({x}) sollte besser als x/ln(x) sein: "
                f"Fehler Li={error_li:.2f}, einfach={error_simple:.2f}"
            )

    def test_li_overcounts(self):
        """Li(x) > π(x) für alle bekannten x (bis ca. 10^316)."""
        for x in [10, 100, 1000]:
            assert logarithmic_integral(x) > prime_counting_function(x), (
                f"Li({x}) sollte > π({x}) sein"
            )

    def test_li_invalid_input(self):
        """Li(x) wirft Fehler für x ≤ 1."""
        with pytest.raises(ValueError):
            logarithmic_integral(0.5)
        with pytest.raises(ValueError):
            logarithmic_integral(1.0)


# ===========================================================================
# VON-MANGOLDT-FUNKTION TESTS
# ===========================================================================

class TestVonMangoldt:
    """Tests für Λ(n)."""

    def test_lambda_at_primes(self):
        """Λ(p) = ln(p) für Primzahlen p."""
        for p, expected in [(2, math.log(2)), (3, math.log(3)),
                            (5, math.log(5)), (7, math.log(7))]:
            result = von_mangoldt_function(p)
            assert abs(result - expected) < 1e-12, f"Λ({p}) = {result}, erwartet {expected}"

    def test_lambda_at_prime_powers(self):
        """Λ(p^k) = ln(p) für Primzahlpotenzen."""
        # Λ(4) = Λ(2²) = ln(2)
        assert abs(von_mangoldt_function(4) - math.log(2)) < 1e-12
        # Λ(8) = Λ(2³) = ln(2)
        assert abs(von_mangoldt_function(8) - math.log(2)) < 1e-12
        # Λ(9) = Λ(3²) = ln(3)
        assert abs(von_mangoldt_function(9) - math.log(3)) < 1e-12

    def test_lambda_at_composites(self):
        """Λ(n) = 0 für n mit verschiedenen Primteilern."""
        for n in [6, 10, 12, 15, 21, 30]:
            assert von_mangoldt_function(n) == 0.0, f"Λ({n}) sollte 0 sein"

    def test_lambda_at_1(self):
        """Λ(1) = 0."""
        assert von_mangoldt_function(1) == 0.0


# ===========================================================================
# TSCHEBYSCHOW-FUNKTIONEN TESTS
# ===========================================================================

class TestChebyshev:
    """Tests für θ(x) und ψ(x)."""

    def test_theta_small_values(self):
        """θ(x) für kleine bekannte x."""
        # θ(2) = ln(2)
        assert abs(chebyshev_theta(2) - math.log(2)) < 1e-12
        # θ(3) = ln(2) + ln(3) = ln(6)
        assert abs(chebyshev_theta(3) - math.log(6)) < 1e-12
        # θ(5) = ln(2) + ln(3) + ln(5) = ln(30)
        assert abs(chebyshev_theta(5) - math.log(30)) < 1e-12

    def test_psi_geq_theta(self):
        """ψ(x) ≥ θ(x) (wegen Primzahlpotenzen)."""
        for x in [10, 50, 100]:
            assert chebyshev_psi(x) >= chebyshev_theta(x), (
                f"ψ({x}) sollte ≥ θ({x}) sein"
            )

    def test_psi_approaches_x(self):
        """ψ(x)/x → 1 (Primzahlsatz)."""
        result = chebyshev_psi_vs_x(1000)
        assert 0.9 < result["ratio"] < 1.1, (
            f"ψ(1000)/1000 sollte nahe 1 sein, erhalten: {result['ratio']}"
        )

    def test_theta_lt_x(self):
        """θ(x) < x für alle x (bekanntes Ergebnis)."""
        for x in [10, 100, 1000]:
            assert chebyshev_theta(x) < x, f"θ({x}) sollte < {x} sein"


# ===========================================================================
# DIRICHLET-CHARAKTERE UND L-FUNKTIONEN TESTS
# ===========================================================================

class TestDirichlet:
    """Tests für Dirichlet-Charaktere und L-Funktionen."""

    def test_principal_character_coprime(self):
        """χ₀(n) = 1 wenn gcd(n, q) = 1."""
        # Modulo 5: χ₀(1)=χ₀(2)=χ₀(3)=χ₀(4)=1, χ₀(5)=0
        for n in [1, 2, 3, 4, 6, 7, 8, 9]:
            assert principal_dirichlet_character(n, 5) == 1

    def test_principal_character_not_coprime(self):
        """χ₀(n) = 0 wenn gcd(n, q) > 1."""
        for n in [5, 10, 15, 20]:
            assert principal_dirichlet_character(n, 5) == 0

    def test_l_function_reduces_to_zeta(self):
        """L(s, χ₀ mod 1) sollte ζ(s) sein."""
        # Für Modul 1 ist χ₀ der triviale Charakter, χ₀(n)=1 für alle n
        # → L(s, χ₀ mod 1) = ζ(s)
        from complex_analysis import _zeta_euler_maclaurin
        s = complex(3, 0)
        chi_values = [1]  # χ₀(1) = 1, periodisch mit Periode 1
        l_val = dirichlet_l_function_partial(s, chi_values, modulus=1, terms=2000)
        zeta_val = _zeta_euler_maclaurin(s, terms=2000)
        assert abs(l_val - zeta_val) < 0.01, (
            f"L(3, χ₀) sollte ζ(3) sein: {l_val} vs {zeta_val}"
        )

    def test_l_function_invalid_s(self):
        """L(s, χ) wirft Fehler für Re(s) ≤ 1."""
        with pytest.raises(ValueError):
            dirichlet_l_function_partial(complex(1, 0), [1, -1], modulus=2)

    def test_dirichlet_equidistribution(self):
        """Primzahlen sind gleichmäßig auf Restklassen verteilt."""
        # π(1000; 4, 1) ≈ π(1000; 4, 3) ≈ π(1000)/2 ≈ 84
        count_1_mod_4 = dirichlet_prime_counting_in_residue_class(1000, 1, 4)
        count_3_mod_4 = dirichlet_prime_counting_in_residue_class(1000, 3, 4)
        # Beide sollten nahe 84 liegen (π(1000) = 168, φ(4) = 2)
        assert 60 <= count_1_mod_4 <= 110, f"π(1000;4,1) = {count_1_mod_4}"
        assert 60 <= count_3_mod_4 <= 110, f"π(1000;4,3) = {count_3_mod_4}"
        # Chebyshev-Bias: Normalerweise count_3 > count_1
        # (Chebyshev-Phänomen, aber nicht immer)


# ===========================================================================
# SIEBMETHODEN UND CHEN-THEOREM TESTS
# ===========================================================================

class TestSieveMethods:
    """Tests für Siebmethoden und Chen-Zerlegung."""

    def test_is_semiprime_known_values(self):
        """Bekannte Semiprimen korrekt erkannt."""
        semiprimes = [4, 6, 9, 10, 14, 15, 21, 22, 25, 26]
        for n in semiprimes:
            assert is_semiprime(n), f"{n} sollte ein Semiprim sein"

    def test_is_semiprime_primes_not_semiprime(self):
        """Primzahlen sind keine Semiprimen."""
        for p in [2, 3, 5, 7, 11, 13]:
            assert not is_semiprime(p), f"{p} (prim) sollte kein Semiprim sein"

    def test_is_semiprime_composites_not_semiprime(self):
        """Zahlen mit ≥ 3 Primfaktoren sind keine Semiprimen."""
        assert not is_semiprime(8)    # 2³: Ω(8) = 3
        assert not is_semiprime(30)   # 2·3·5: ω(30) = 3
        assert not is_semiprime(12)   # 2²·3: Ω(12) = 3

    def test_chen_decomposition_goldbach_case(self):
        """Chen-Zerlegung findet Goldbach (1+1) wenn möglich."""
        # 12 = 5 + 7 (beide prim)
        result = chen_decomposition(12)
        assert result is not None
        p, m, typ = result
        assert p + m == 12

    def test_chen_decomposition_semiprime_case(self):
        """Chen-Zerlegung findet 1+2-Zerlegung."""
        # Für gerade Zahlen: mindestens eine Chen-Zerlegung existiert
        for n in [4, 6, 8, 10, 20, 50, 100]:
            result = chen_decomposition(n)
            assert result is not None, f"Chen-Zerlegung für {n} nicht gefunden"
            p, m, typ = result
            assert p + m == n
            assert typ in ('prime', 'semiprime')

    def test_chen_decomposition_invalid(self):
        """Chen-Zerlegung wirft Fehler für ungerade oder ≤ 2."""
        with pytest.raises(ValueError):
            chen_decomposition(7)
        with pytest.raises(ValueError):
            chen_decomposition(2)

    def test_count_omega_known(self):
        """ω(n) korrekt berechnet."""
        assert count_omega(1) == 0
        assert count_omega(2) == 1
        assert count_omega(4) == 1    # 2²: nur Primteiler 2
        assert count_omega(6) == 2    # 2·3
        assert count_omega(30) == 3   # 2·3·5

    def test_almost_prime_primes(self):
        """Primzahlen sind 1-fast-prim."""
        for p in [2, 3, 5, 7, 11, 13]:
            assert almost_prime(p, 1), f"{p} sollte 1-fast-prim sein"

    def test_almost_prime_semiprimes(self):
        """Semiprimen sind 2-fast-prim."""
        for n in [4, 6, 9, 10, 15]:
            assert almost_prime(n, 2), f"{n} sollte 2-fast-prim sein"


# ===========================================================================
# PRIMZAHL-LÜCKEN TESTS
# ===========================================================================

class TestPrimeGaps:
    """Tests für Primzahl-Lücken-Analyse."""

    def test_first_gap_is_1(self):
        """Erste Lücke (2→3) hat Größe 1."""
        gaps = prime_gaps(10)
        assert gaps[0]["gap"] == 1
        assert gaps[0]["prime"] == 2

    def test_all_gaps_positive(self):
        """Alle Lücken sind positiv."""
        gaps = prime_gaps(100)
        for g in gaps:
            assert g["gap"] > 0

    def test_maximal_gap_in_range(self):
        """Größte Lücke bis 100 ist 8 (zwischen 89 und 97)."""
        max_g = maximal_prime_gap(100)
        assert max_g["gap"] == 8, f"Größte Lücke bis 100 ist 8, erhalten: {max_g['gap']}"
        assert max_g["prime"] == 89

    def test_cramer_bound_plausible(self):
        """Cramér-Vermutung: Lücken wachsen höchstens wie (ln p)² für große p.

        Für kleine p (z.B. p=2, Lücke=1 aber (ln 2)²=0.48) gilt die
        asymptotische Schranke noch nicht. Der Test prüft daher nur für p ≥ 10.
        """
        gaps = prime_gaps(200)
        for g in gaps:
            if g["prime"] < 10:
                continue   # Cramér-Schranke ist asymptotisch, nicht für p<10
            # Großzügige Schranke: gap ≤ 4·(ln p)² (enthält Cramér-Konstante C)
            assert g["gap"] <= 4 * g["cramer_bound"], (
                f"Lücke {g['gap']} bei p={g['prime']} überschreitet 4·(ln p)²={4*g['cramer_bound']:.1f}"
            )
