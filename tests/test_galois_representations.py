"""
@file test_galois_representations.py
@brief Tests für das Galois-Darstellungen-Modul (galois_representations.py).
@description
    Umfassende pytest-Tests für alle Klassen und Funktionen des Moduls,
    inklusive mathematisch korrekter Testfälle aus der Theorie der
    Galois-Darstellungen und des Langlands-Programms.

    Testbereiche:
    - GaloisRepresentation (abstrakte Klasse via TateModule)
    - TateModuleRepresentation mit echten a_p-Werten für y²=x³-x
    - CyclotomicCharacter
    - DirichletCharacterRepresentation
    - SymmetricPowerRepresentation
    - LanglandsCorrespondence
    - Hilfsfunktionen: p_adic_matrix, galois_group_order, artin_conductor, weil_deligne

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
import math
import cmath

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import numpy as np

from galois_representations import (
    GaloisRepresentation,
    TateModuleRepresentation,
    CyclotomicCharacter,
    CyclotomicCharacterPower,
    DirichletCharacterRepresentation,
    SymmetricPowerRepresentation,
    LanglandsCorrespondence,
    p_adic_representation_matrix,
    galois_group_order,
    artin_conductor,
    weil_deligne_representation_demo,
    compute_ap_for_y2_x3_minus_x,
    _is_prime,
    _primes_up_to,
    _sym_power_trace,
)

# ===========================================================================
# FIXTURE: Elliptische Kurve y² = x³ - x (Conductor 32)
# ===========================================================================

# Bekannte a_p-Werte für y²=x³-x:
# p≡3 (mod 4) → a_p=0
# p≡1 (mod 4) → a_p aus Gaußscher Primzahl
AP_VALUES_E_MINUS_X = {
    5: -2,    # 5≡1 mod 4
    13: 6,    # 13≡1 mod 4
    17: 2,    # 17≡1 mod 4
    29: 6,    # 29≡1 mod 4
    37: 2,    # 37≡1 mod 4
    41: 10,   # 41≡1 mod 4
    53: -6,   # 53≡1 mod 4
    3: 0,     # 3≡3 mod 4
    7: 0,     # 7≡3 mod 4
    11: 0,    # 11≡3 mod 4
    19: 0,    # 19≡3 mod 4
    23: 0,    # 23≡3 mod 4
    31: 0,    # 31≡3 mod 4
    43: 0,    # 43≡3 mod 4
    47: 0,    # 47≡3 mod 4
}


@pytest.fixture
def tate_module_e():
    """Tate-Modul der elliptischen Kurve y²=x³-x mit ℓ=3."""
    return TateModuleRepresentation(AP_VALUES_E_MINUS_X, conductor=32, ell=3)


@pytest.fixture
def cyclotomic_5():
    """Zyklotomischer Charakter mod 5."""
    return CyclotomicCharacter(5)


@pytest.fixture
def dirichlet_chi4():
    """Dirichlet-Charakter mod 4: χ(1)=1, χ(3)=-1, χ(0)=χ(2)=0."""
    chi_values = {1: 1, 2: 0, 3: -1, 0: 0}
    return DirichletCharacterRepresentation(chi_values, conductor=4, ell=3)


@pytest.fixture
def langlands():
    """LanglandsCorrespondence-Instanz."""
    return LanglandsCorrespondence()


# ===========================================================================
# ABSCHNITT 1: HILFSFUNKTIONEN
# ===========================================================================

class TestHelpers:
    """Tests für interne Hilfsfunktionen."""

    def test_is_prime_small(self):
        """Kleine Primzahlen korrekt erkannt."""
        assert _is_prime(2)
        assert _is_prime(3)
        assert _is_prime(5)
        assert _is_prime(7)
        assert _is_prime(11)

    def test_is_prime_composite(self):
        """Zusammengesetzte Zahlen korrekt abgelehnt."""
        assert not _is_prime(1)
        assert not _is_prime(4)
        assert not _is_prime(9)
        assert not _is_prime(25)
        assert not _is_prime(100)

    def test_is_prime_edge_cases(self):
        """Edge-Cases: 0, 1, negative Zahlen."""
        assert not _is_prime(0)
        assert not _is_prime(-1)
        assert not _is_prime(-7)

    def test_primes_up_to_10(self):
        """Primzahlen bis 10."""
        assert _primes_up_to(10) == [2, 3, 5, 7]

    def test_primes_up_to_1(self):
        """Primzahlen bis 1: leer."""
        assert _primes_up_to(1) == []

    def test_primes_up_to_count(self):
        """Anzahl der Primzahlen bis 100 ist 25."""
        assert len(_primes_up_to(100)) == 25

    def test_sym_power_trace_k0(self):
        """Sym^0: Spur = 1 (triviale Darstellung)."""
        result = _sym_power_trace(2+0j, 3+0j, 0)
        assert abs(result - 1.0) < 1e-10

    def test_sym_power_trace_k1(self):
        """Sym^1: Spur = α + β (= ursprüngliche Darstellung)."""
        result = _sym_power_trace(2+0j, 3+0j, 1)
        assert abs(result - 5.0) < 1e-10

    def test_sym_power_trace_k2(self):
        """Sym^2: Spur = α² + αβ + β² = (α+β)² - αβ."""
        # α=2, β=3: α²+αβ+β²=4+6+9=19
        result = _sym_power_trace(2+0j, 3+0j, 2)
        assert abs(result - 19.0) < 1e-10

    def test_compute_ap_congruent_3_mod4(self):
        """a_p = 0 für p ≡ 3 (mod 4) bei y²=x³-x."""
        for p in [3, 7, 11, 19, 23]:
            a_p = compute_ap_for_y2_x3_minus_x(p)
            assert a_p == 0, f"Erwartet a_{p}=0, erhalten {a_p}"

    def test_compute_ap_p5(self):
        """a_5 = -2 für y²=x³-x."""
        assert compute_ap_for_y2_x3_minus_x(5) == -2

    def test_compute_ap_p13(self):
        """a_13 = 6 für y²=x³-x."""
        assert compute_ap_for_y2_x3_minus_x(13) == 6

    def test_compute_ap_hasse_bound(self):
        """Alle a_p erfüllen Hasse-Schranke |a_p| ≤ 2√p."""
        for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
            a_p = compute_ap_for_y2_x3_minus_x(p)
            assert abs(a_p) <= 2 * math.sqrt(p) + 1e-9, \
                f"Hasse-Schranke verletzt bei p={p}: |a_p|={abs(a_p)} > 2√{p}"


# ===========================================================================
# ABSCHNITT 2: GaloisRepresentation (Abstrakte Klasse)
# ===========================================================================

class TestGaloisRepresentationBase:
    """Tests für die abstrakte Basisklasse via TateModule."""

    def test_invalid_ell_not_prime(self):
        """ValueError bei nicht-primer ℓ."""
        with pytest.raises(ValueError, match="Primzahl"):
            TateModuleRepresentation({}, 1, ell=4)

    def test_invalid_dimension_zero(self):
        """ValueError bei Dimension 0 (direkt via Klasse nicht möglich, via Subklasse)."""
        # TateModule hat immer dim=2, teste ell=1 (nicht prim)
        with pytest.raises(ValueError):
            TateModuleRepresentation({}, 1, ell=1)

    def test_repr_format(self, tate_module_e):
        """__repr__ gibt sinnvollen String zurück."""
        r = repr(tate_module_e)
        assert "GaloisRepresentation" in r or "TateModule" in r

    def test_ell_attribute(self, tate_module_e):
        """ell-Attribut korrekt gesetzt."""
        assert tate_module_e.ell == 3

    def test_dim_attribute(self, tate_module_e):
        """dim-Attribut für TateModule ist 2."""
        assert tate_module_e.dim == 2


# ===========================================================================
# ABSCHNITT 3: TateModuleRepresentation
# ===========================================================================

class TestTateModuleRepresentation:
    """Tests für den Tate-Modul einer elliptischen Kurve."""

    def test_trace_frobenius_good_prime(self, tate_module_e):
        """a_p aus dict für gute Primzahl."""
        assert tate_module_e.trace_frobenius(5) == -2
        assert tate_module_e.trace_frobenius(13) == 6

    def test_trace_frobenius_p_congruent_3_mod4(self, tate_module_e):
        """a_p = 0 für p ≡ 3 mod 4 bei y²=x³-x."""
        for p in [3, 7, 11, 19, 23]:
            assert tate_module_e.trace_frobenius(p) == 0, \
                f"Erwartet 0 für p={p}, erhalten {tate_module_e.trace_frobenius(p)}"

    def test_trace_frobenius_bad_prime(self, tate_module_e):
        """a_p = 0 für schlechte Primzahl (p=2, Führer=32)."""
        # p=2 ist Verzweigungsprimzahl (Führer 32 = 2^5), nicht im dict
        assert tate_module_e.trace_frobenius(2) == 0

    def test_trace_frobenius_unknown_prime(self, tate_module_e):
        """Unbekannte Primzahl (nicht im dict): gibt 0 zurück."""
        assert tate_module_e.trace_frobenius(997) == 0

    def test_characteristic_polynomial_p5(self, tate_module_e):
        """char. Polynom bei p=5: x² + 2x + 5 (Koeff: 1, -a_p=-(-2)=2, p=5)."""
        coeffs = tate_module_e.characteristic_polynomial_frobenius(5)
        assert coeffs[0] == 1.0   # Leitkoeffizient
        assert coeffs[1] == 2.0   # -a_p = -(-2) = 2
        assert coeffs[2] == 5.0   # p = 5

    def test_characteristic_polynomial_p13(self, tate_module_e):
        """char. Polynom bei p=13: 1, -6, 13."""
        coeffs = tate_module_e.characteristic_polynomial_frobenius(13)
        assert coeffs[1] == -6.0

    def test_is_unramified_good_prime(self, tate_module_e):
        """Gute Primzahl p=5: unverzweigt."""
        assert tate_module_e.is_unramified(5)

    def test_is_unramified_bad_prime(self, tate_module_e):
        """Schlechte Primzahl p=2 (teilt 32): verzweigt."""
        assert not tate_module_e.is_unramified(2)

    def test_is_unramified_ell(self, tate_module_e):
        """p = ℓ = 3: verzweigt (ℓ-adische Darstellung ist bei ℓ ramifiziert)."""
        assert not tate_module_e.is_unramified(3)

    def test_frobenius_eigenvalues_hasse(self, tate_module_e):
        """Frobenius-Eigenwerte liegen auf Kreis mit Radius √p (Hasse)."""
        for p in [5, 13, 17, 29, 37]:
            if tate_module_e.is_unramified(p):
                alpha, beta = tate_module_e.frobenius_eigenvalues(p)
                assert abs(abs(alpha) - math.sqrt(p)) < 1e-8, \
                    f"|α| ≠ √{p} bei p={p}: |α|={abs(alpha)}"
                assert abs(abs(beta) - math.sqrt(p)) < 1e-8, \
                    f"|β| ≠ √{p} bei p={p}: |β|={abs(beta)}"

    def test_frobenius_eigenvalues_product(self, tate_module_e):
        """Produkt der Frobenius-Eigenwerte = p (Determinante = Normen-Charakter)."""
        for p in [5, 13, 17]:
            alpha, beta = tate_module_e.frobenius_eigenvalues(p)
            assert abs(alpha * beta - p) < 1e-8, \
                f"α·β ≠ {p} bei p={p}"

    def test_frobenius_eigenvalues_sum(self, tate_module_e):
        """Summe der Frobenius-Eigenwerte = a_p."""
        for p in [5, 13, 17]:
            a_p = AP_VALUES_E_MINUS_X[p]
            alpha, beta = tate_module_e.frobenius_eigenvalues(p)
            assert abs((alpha + beta) - a_p) < 1e-8, \
                f"α+β ≠ a_p={a_p} bei p={p}"

    def test_conductor(self, tate_module_e):
        """Führer N=32 korrekt gespeichert."""
        assert tate_module_e.conductor == 32

    def test_invalid_conductor(self):
        """ValueError bei Führer < 1."""
        with pytest.raises(ValueError, match="Führer"):
            TateModuleRepresentation({}, conductor=0, ell=3)

    def test_l_function_euler_factor_good_prime(self, tate_module_e):
        """Euler-Faktor bei guter Primzahl p=5, s=2: endlicher Wert."""
        factor = tate_module_e.l_function_euler_factor(5, 2.0)
        assert cmath.isfinite(factor)
        assert abs(factor) > 0

    def test_l_function_euler_factor_bad_prime(self, tate_module_e):
        """Euler-Faktor bei verzweigter Primzahl p=2: Rückgabe 1."""
        factor = tate_module_e.l_function_euler_factor(2, 2.0)
        assert abs(factor - 1.0) < 1e-10


# ===========================================================================
# ABSCHNITT 4: CyclotomicCharacter
# ===========================================================================

class TestCyclotomicCharacter:
    """Tests für den zyklotomischen Charakter."""

    def test_evaluate_basic(self, cyclotomic_5):
        """χ_5(Frob_p) = p mod 5."""
        assert cyclotomic_5.evaluate(7) == 7 % 5   # = 2
        assert cyclotomic_5.evaluate(11) == 11 % 5  # = 1
        assert cyclotomic_5.evaluate(13) == 13 % 5  # = 3

    def test_evaluate_p_mod_ell(self):
        """p mod ℓ für verschiedene Charaktere."""
        for ell in [2, 3, 5, 7, 11]:
            chi = CyclotomicCharacter(ell)
            for p in [3, 5, 7, 11, 13, 17, 19]:
                if p != ell:
                    assert chi.evaluate(p) == p % ell

    def test_is_odd(self, cyclotomic_5):
        """Zyklotomischer Charakter ist immer ungerade."""
        assert cyclotomic_5.is_odd() is True

    def test_is_odd_all_ell(self):
        """Ungerade für alle Primzahlen ℓ."""
        for ell in [2, 3, 5, 7, 11, 13]:
            chi = CyclotomicCharacter(ell)
            assert chi.is_odd()

    def test_invalid_ell(self):
        """ValueError bei nicht-primer ℓ."""
        with pytest.raises(ValueError, match="Primzahl"):
            CyclotomicCharacter(4)

    def test_evaluate_nonzero(self):
        """p mod ℓ ≠ 0 für p ≠ ℓ (p teilt nicht ℓ, da beide prim)."""
        chi = CyclotomicCharacter(7)
        for p in [2, 3, 5, 11, 13]:
            assert chi.evaluate(p) != 0

    def test_tensor_power_k1(self):
        """χ^1(Frob_p) = p mod ℓ."""
        chi = CyclotomicCharacter(7)
        power = chi.tensor_power(1)
        assert power.evaluate(11) == pow(11, 1, 7)

    def test_tensor_power_k2(self):
        """χ^2(Frob_p) = p² mod ℓ."""
        chi = CyclotomicCharacter(7)
        power = chi.tensor_power(2)
        assert power.evaluate(11) == pow(11, 2, 7)

    def test_repr(self, cyclotomic_5):
        """__repr__ enthält ell."""
        r = repr(cyclotomic_5)
        assert "5" in r


# ===========================================================================
# ABSCHNITT 5: DirichletCharacterRepresentation
# ===========================================================================

class TestDirichletCharacterRepresentation:
    """Tests für die Dirichlet-Charakter-Darstellung."""

    def test_trace_frobenius_chi4(self, dirichlet_chi4):
        """χ_4(3) = -1, χ_4(5)=χ_4(1)=1."""
        assert dirichlet_chi4.trace_frobenius(3) == -1   # 3 mod 4 = 3
        assert dirichlet_chi4.trace_frobenius(5) == 1    # 5 mod 4 = 1
        assert dirichlet_chi4.trace_frobenius(7) == -1   # 7 mod 4 = 3
        assert dirichlet_chi4.trace_frobenius(13) == 1   # 13 mod 4 = 1

    def test_trace_frobenius_ramified(self, dirichlet_chi4):
        """χ(p) = 0 wenn p | Führer (p=2 bei mod-4-Charakter)."""
        assert dirichlet_chi4.trace_frobenius(2) == 0

    def test_is_unramified_good_prime(self, dirichlet_chi4):
        """p=5 (∤ 4, ≠ 3): unverzweigt."""
        assert dirichlet_chi4.is_unramified(5)

    def test_is_unramified_ramified(self, dirichlet_chi4):
        """p=2 (| 4): verzweigt."""
        assert not dirichlet_chi4.is_unramified(2)

    def test_dimension_1(self, dirichlet_chi4):
        """Dirichlet-Charakter ist 1-dimensional."""
        assert dirichlet_chi4.dim == 1

    def test_conductor(self, dirichlet_chi4):
        """Führer 4 korrekt gespeichert."""
        assert dirichlet_chi4.conductor == 4

    def test_invalid_conductor(self):
        """ValueError bei Führer 0."""
        with pytest.raises(ValueError, match="Führer"):
            DirichletCharacterRepresentation({}, conductor=0, ell=3)

    def test_periodicity(self):
        """Charakter ist periodisch: χ(p) = χ(p mod q)."""
        chi_values = {1: 1, 2: 0, 3: -1, 0: 0}
        chi = DirichletCharacterRepresentation(chi_values, conductor=4, ell=5)
        # 17 ≡ 1 (mod 4) → χ(17) = χ(1) = 1
        assert chi.trace_frobenius(17) == 1
        # 19 ≡ 3 (mod 4) → χ(19) = χ(3) = -1
        assert chi.trace_frobenius(19) == -1

    def test_l_function_euler_factor_dim1(self, dirichlet_chi4):
        """Euler-Faktor für 1-dim. Darstellung: (1 - χ(p)·p^{-s})^{-1}."""
        p, s = 5, 2.0
        chi_p = dirichlet_chi4.trace_frobenius(p)  # = 1
        expected = 1.0 / (1.0 - chi_p * p**(-s))
        computed = dirichlet_chi4.l_function_euler_factor(p, s)
        assert abs(computed - expected) < 1e-10


# ===========================================================================
# ABSCHNITT 6: SymmetricPowerRepresentation
# ===========================================================================

class TestSymmetricPowerRepresentation:
    """Tests für die symmetrische Potenz einer 2-dim. Darstellung."""

    def test_dimension_k0(self, tate_module_e):
        """Sym^0 hat Dimension 1."""
        sym0 = SymmetricPowerRepresentation(tate_module_e, 0)
        assert sym0.dimension() == 1
        assert sym0.dim == 1

    def test_dimension_k1(self, tate_module_e):
        """Sym^1 hat Dimension 2 (= ursprüngliche Darstellung)."""
        sym1 = SymmetricPowerRepresentation(tate_module_e, 1)
        assert sym1.dimension() == 2

    def test_dimension_k2(self, tate_module_e):
        """Sym^2 hat Dimension 3."""
        sym2 = SymmetricPowerRepresentation(tate_module_e, 2)
        assert sym2.dimension() == 3

    def test_dimension_k4(self, tate_module_e):
        """Sym^4 hat Dimension 5."""
        sym4 = SymmetricPowerRepresentation(tate_module_e, 4)
        assert sym4.dimension() == 5

    def test_dimension_general(self, tate_module_e):
        """Sym^k hat Dimension k+1 für k=0,...,6."""
        for k in range(7):
            sym_k = SymmetricPowerRepresentation(tate_module_e, k)
            assert sym_k.dimension() == k + 1

    def test_trace_sym0(self, tate_module_e):
        """Sym^0: Frobenius-Spur = 1."""
        sym0 = SymmetricPowerRepresentation(tate_module_e, 0)
        trace = sym0.trace_frobenius(5)
        assert abs(trace - 1.0) < 1e-10

    def test_trace_sym1_equals_ap(self, tate_module_e):
        """Sym^1: Frobenius-Spur = a_p (ursprüngliche Darstellung)."""
        sym1 = SymmetricPowerRepresentation(tate_module_e, 1)
        for p in [5, 13, 17]:
            trace = sym1.trace_frobenius(p)
            a_p = AP_VALUES_E_MINUS_X[p]
            assert abs(trace - a_p) < 1e-8, \
                f"Sym^1 Spur ≠ a_p={a_p} bei p={p}"

    def test_trace_sym2_formula(self, tate_module_e):
        """Sym^2: Spur = α² + αβ + β²."""
        sym2 = SymmetricPowerRepresentation(tate_module_e, 2)
        p = 5
        alpha, beta = tate_module_e.frobenius_eigenvalues(p)
        expected = alpha**2 + alpha*beta + beta**2
        computed = sym2.trace_frobenius(p)
        assert abs(computed - expected) < 1e-8

    def test_invalid_base_dim(self):
        """ValueError wenn Basis-Darstellung nicht 2-dimensional."""
        chi_values = {1: 1, 3: -1}
        chi = DirichletCharacterRepresentation(chi_values, 4, 3)
        with pytest.raises(ValueError, match="2-dimensional"):
            SymmetricPowerRepresentation(chi, 2)

    def test_invalid_k_negative(self, tate_module_e):
        """ValueError für k < 0."""
        with pytest.raises(ValueError, match="≥ 0"):
            SymmetricPowerRepresentation(tate_module_e, -1)

    def test_ell_inherited(self, tate_module_e):
        """ℓ von der Basis-Darstellung geerbt."""
        sym2 = SymmetricPowerRepresentation(tate_module_e, 2)
        assert sym2.ell == tate_module_e.ell


# ===========================================================================
# ABSCHNITT 7: LanglandsCorrespondence
# ===========================================================================

class TestLanglandsCorrespondence:
    """Tests für die Langlands-Korrespondenz-Klasse."""

    def test_automorphic_l_function_converges(self, langlands, tate_module_e):
        """L(ρ, 2) konvergiert für Re(s)=2 > 1."""
        result = langlands.automorphic_l_function(tate_module_e, 2.0, num_primes=20)
        assert cmath.isfinite(result)
        assert abs(result) > 0

    def test_automorphic_l_function_positive(self, langlands, tate_module_e):
        """L(ρ, s) > 0 für reelles s >> 1 und reelle Koeffizienten."""
        result = langlands.automorphic_l_function(tate_module_e, 3.0, num_primes=15)
        # Für s=3 sollte das Euler-Produkt einen positiven reellen Wert annähern
        assert abs(result) > 0.1

    def test_check_ramanujan_conjecture_satisfied(self, langlands, tate_module_e):
        """Ramanujan-Schranke |a_p| ≤ 2√p erfüllt für alle guten p."""
        primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
        result = langlands.check_ramanujan_conjecture(tate_module_e, primes)
        assert result["all_satisfied"], \
            f"Ramanujan verletzt: {result['violations']}"
        assert result["num_checked"] > 0

    def test_ramanujan_max_ratio(self, langlands, tate_module_e):
        """Maximales |a_p|/(2√p) ≤ 1."""
        primes = [5, 13, 17, 29, 37, 41]
        result = langlands.check_ramanujan_conjecture(tate_module_e, primes)
        assert result["max_ratio"] <= 1.0 + 1e-9

    def test_ramanujan_skips_bad_primes(self, langlands, tate_module_e):
        """Verzweigte Primzahlen werden bei Ramanujan-Check übersprungen."""
        # p=2 ist verzweigt (Führer=32), wird übersprungen
        result = langlands.check_ramanujan_conjecture(tate_module_e, [2])
        assert result["num_checked"] == 0

    def test_check_functional_equation_returns_dict(self, langlands, tate_module_e):
        """Funktionalgleichungs-Check gibt korrektes Dict zurück."""
        result = langlands.check_functional_equation(tate_module_e, 1.5)
        assert "lambda_s" in result
        assert "lambda_1ms" in result
        assert "epsilon_estimate" in result

    def test_functional_equation_lambda_finite(self, langlands, tate_module_e):
        """Vollständige L-Funktion Λ ist endlich für s=1.5."""
        result = langlands.check_functional_equation(tate_module_e, 1.5)
        assert cmath.isfinite(result["lambda_s"])


# ===========================================================================
# ABSCHNITT 8: p-adische Darstellungsmatrix
# ===========================================================================

class TestPAdicRepresentationMatrix:
    """Tests für p_adic_representation_matrix."""

    def test_shape_2x2(self):
        """2×2-Matrix für n=2."""
        A = p_adic_representation_matrix(2, 3, 2)
        assert A.shape == (2, 2)

    def test_shape_3x3(self):
        """3×3-Matrix für n=3."""
        A = p_adic_representation_matrix(3, 5, 2)
        assert A.shape == (3, 3)

    def test_invertible_mod_p(self):
        """Matrix ist invertierbar mod p."""
        A = p_adic_representation_matrix(2, 3, 2)
        det_mod_p = int(round(np.linalg.det(A.astype(float)))) % 3
        assert det_mod_p != 0

    def test_entries_in_range(self):
        """Einträge liegen in {0, ..., p^k - 1}."""
        p, k = 3, 2
        A = p_adic_representation_matrix(2, p, k)
        modulus = p ** k
        assert np.all(A >= 0)
        assert np.all(A < modulus)

    def test_invalid_n(self):
        """ValueError für n < 1."""
        with pytest.raises(ValueError, match="n muss"):
            p_adic_representation_matrix(0, 3, 2)

    def test_invalid_p(self):
        """ValueError wenn p keine Primzahl."""
        with pytest.raises(ValueError, match="Primzahl"):
            p_adic_representation_matrix(2, 4, 2)

    def test_invalid_precision(self):
        """ValueError für precision < 1."""
        with pytest.raises(ValueError, match="precision"):
            p_adic_representation_matrix(2, 3, 0)

    def test_1x1_matrix(self):
        """1×1-Matrix: Skalar in ℤ/p^k."""
        A = p_adic_representation_matrix(1, 5, 3)
        assert A.shape == (1, 1)


# ===========================================================================
# ABSCHNITT 9: Galois-Gruppen-Ordnung
# ===========================================================================

class TestGaloisGroupOrder:
    """Tests für galois_group_order."""

    def test_x2_minus_2(self):
        """Gal(ℚ(√2)/ℚ) ≅ ℤ/2ℤ hat Ordnung 2."""
        order = galois_group_order("x**2 - 2")
        assert order == 2

    def test_x2_minus_3(self):
        """Gal(ℚ(√3)/ℚ) ≅ ℤ/2ℤ hat Ordnung 2."""
        order = galois_group_order("x**2 - 3")
        assert order == 2

    def test_x3_minus_2(self):
        """Gal(ℚ(∛2)/ℚ) ≅ S₃ hat Ordnung 6."""
        order = galois_group_order("x**3 - 2")
        assert order == 6

    def test_x4_minus_2(self):
        """Gal(ℚ(⁴√2)/ℚ) ≅ D₄ hat Ordnung 8."""
        order = galois_group_order("x**4 - 2")
        assert order == 8

    def test_linear_polynomial(self):
        """Lineares Polynom: Galois-Gruppe trivial, Ordnung 1."""
        order = galois_group_order("x - 3")
        assert order == 1

    def test_invalid_polynomial(self):
        """ValueError bei ungültigem Polynom."""
        with pytest.raises((ValueError, Exception)):
            galois_group_order("not_a_polynomial!!!")


# ===========================================================================
# ABSCHNITT 10: Artin-Leiter
# ===========================================================================

class TestArtinConductor:
    """Tests für artin_conductor."""

    def test_tate_module_returns_conductor(self, tate_module_e):
        """TateModule gibt exakten Führer N=32 zurück."""
        primes = [2, 3, 5, 7]
        result = artin_conductor(tate_module_e, primes)
        assert result == 32

    def test_dirichlet_returns_conductor(self, dirichlet_chi4):
        """DirichletCharacter gibt Führer q=4 zurück."""
        result = artin_conductor(dirichlet_chi4, [2, 3, 5])
        assert result == 4

    def test_general_rep_product_ramified(self, tate_module_e):
        """Artin-Leiter für TateModule mit anderem Führer."""
        rep = TateModuleRepresentation({5: -2}, conductor=15, ell=7)
        result = artin_conductor(rep, [3, 5, 7, 11])
        assert result == 15  # Führer direkt


# ===========================================================================
# ABSCHNITT 11: Weil-Deligne-Darstellung
# ===========================================================================

class TestWeilDeligneRepresentation:
    """Tests für weil_deligne_representation_demo."""

    def test_returns_dict(self):
        """Rückgabe ist ein Dict mit korrekten Schlüsseln."""
        result = weil_deligne_representation_demo(5)
        assert "prime" in result
        assert "frobenius_matrix" in result
        assert "monodromy_matrix" in result
        assert "type" in result
        assert "description" in result
        assert "nilpotent_check" in result

    def test_prime_field(self):
        """prime-Feld enthält die eingegebene Primzahl."""
        result = weil_deligne_representation_demo(7)
        assert result["prime"] == 7

    def test_frobenius_matrix_shape(self):
        """Frobenius-Matrix ist 2×2."""
        result = weil_deligne_representation_demo(5)
        assert result["frobenius_matrix"].shape == (2, 2)

    def test_monodromy_nilpotent(self):
        """N² = 0 (Monodromie ist nilpotent)."""
        result = weil_deligne_representation_demo(5)
        N = result["monodromy_matrix"]
        N_squared = N @ N
        assert np.allclose(N_squared, np.zeros((2, 2))), \
            f"N² ≠ 0: N²={N_squared}"
        assert result["nilpotent_check"] is True

    def test_good_reduction_p5(self):
        """p=5 hat gute Reduktion bei y²=x³-x."""
        result = weil_deligne_representation_demo(5)
        assert result["type"] == "good"

    def test_good_reduction_monodromy_zero(self):
        """Bei guter Reduktion ist N = 0."""
        result = weil_deligne_representation_demo(5)
        N = result["monodromy_matrix"]
        assert np.allclose(N, np.zeros((2, 2)))

    def test_additive_reduction_p2(self):
        """p=2 ist Verzweigungsprimzahl (additive Reduktion)."""
        result = weil_deligne_representation_demo(2)
        assert result["type"] == "additive"

    def test_invalid_prime(self):
        """ValueError für nicht-prime p."""
        with pytest.raises(ValueError, match="Primzahl"):
            weil_deligne_representation_demo(4)

    def test_description_contains_prime(self):
        """description enthält die Primzahl."""
        result = weil_deligne_representation_demo(13)
        assert "13" in result["description"]

    def test_frobenius_eigenvalues_good_p13(self):
        """Bei p=13 (gute Reduktion): Frobenius-Eigenwerte auf √13-Kreis."""
        result = weil_deligne_representation_demo(13)
        if result["type"] == "good":
            F = result["frobenius_matrix"]
            # Diagonal-Einträge sind die Eigenwerte α, β
            alpha = F[0, 0]
            beta = F[1, 1]
            sqrt_p = math.sqrt(13)
            assert abs(abs(alpha) - sqrt_p) < 0.5 or abs(abs(beta) - sqrt_p) < 0.5
