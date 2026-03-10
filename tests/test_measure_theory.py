"""
Tests für das Maßtheorie-Modul (src/measure_theory.py).

Testet σ-Algebren, Maße, Lebesgue-Integral, Konvergenzsätze und L^p-Räume.

Autor: Kurt Ingwer
Letzte Änderung: 2026-03-10
"""

import pytest
import numpy as np
from fractions import Fraction

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from measure_theory import (
    SigmaAlgebra, Measure, LebesgueMeasure, lebesgue_integral,
    riemann_vs_lebesgue, MeasurableFunction, monotone_convergence_theorem_demo,
    dominated_convergence_theorem_demo, fubini_theorem_demo, product_measure,
    LpSpace, radon_nikodym_theorem_demo, cantor_function_demo,
)


# ---------------------------------------------------------------------------
# SigmaAlgebra Tests
# ---------------------------------------------------------------------------

class TestSigmaAlgebra:
    """Tests für die SigmaAlgebra-Klasse."""

    def setup_method(self):
        """Erstellt Test-σ-Algebren."""
        self.omega = frozenset([1, 2, 3])
        # Minimale σ-Algebra: {∅, Ω}
        self.trivial = SigmaAlgebra(self.omega, [frozenset(), self.omega])
        # Potenzmenge (vollständige σ-Algebra)
        self.power_set = SigmaAlgebra(self.omega, [
            frozenset(), frozenset([1]), frozenset([2]), frozenset([3]),
            frozenset([1, 2]), frozenset([1, 3]), frozenset([2, 3]), self.omega
        ])

    def test_trivial_sigma_algebra(self):
        """Die triviale σ-Algebra {∅, Ω} ist gültig."""
        assert self.trivial.is_sigma_algebra()

    def test_power_set_is_sigma_algebra(self):
        """Die Potenzmenge ist eine σ-Algebra."""
        assert self.power_set.is_sigma_algebra()

    def test_missing_omega_is_not_sigma_algebra(self):
        """σ-Algebra ohne Ω ist ungültig."""
        bad = SigmaAlgebra(self.omega, [frozenset(), frozenset([1])])
        assert not bad.is_sigma_algebra()

    def test_closed_under_complement(self):
        """Triviale σ-Algebra ist unter Komplement abgeschlossen."""
        assert self.trivial.is_closed_under_complement()

    def test_not_closed_under_complement(self):
        """Menge ohne Komplement ist nicht unter Komplement abgeschlossen."""
        bad = SigmaAlgebra(self.omega, [frozenset(), self.omega, frozenset([1])])
        assert not bad.is_closed_under_complement()

    def test_closed_under_union(self):
        """Potenzmenge ist unter Vereinigung abgeschlossen."""
        assert self.power_set.is_closed_under_countable_union()

    def test_generated_sigma_algebra(self):
        """Erzeugte σ-Algebra enthält Generatoren und ist abgeschlossen."""
        gen = SigmaAlgebra(self.omega, [])
        generated = gen.generated_sigma_algebra([frozenset([1])])
        assert generated.is_sigma_algebra()
        assert frozenset([1]) in generated.sets
        assert frozenset([2, 3]) in generated.sets  # Komplement

    def test_generated_sigma_algebra_with_two_generators(self):
        """Erzeugte σ-Algebra mit zwei Generatoren."""
        gen = SigmaAlgebra(self.omega, [])
        generated = gen.generated_sigma_algebra([frozenset([1]), frozenset([2])])
        assert generated.is_sigma_algebra()

    def test_borel_sets_finite(self):
        """Borel-σ-Algebra auf endlicher Menge = Potenzmenge."""
        borel = self.trivial.borel_sets_finite()
        assert borel.is_sigma_algebra()
        # Potenzmenge von {1,2,3} hat 2³ = 8 Elemente
        assert len(borel.sets) == 8

    def test_sigma_algebra_repr(self):
        """Repr gibt sinnvolle Ausgabe."""
        r = repr(self.trivial)
        assert 'SigmaAlgebra' in r

    def test_empty_set_always_in_sigma_algebra(self):
        """Leere Menge ist immer in der σ-Algebra."""
        assert frozenset() in self.trivial.sets


# ---------------------------------------------------------------------------
# Measure Tests
# ---------------------------------------------------------------------------

class TestMeasure:
    """Tests für die Measure-Klasse."""

    def setup_method(self):
        """Erstellt Test-Maße."""
        omega = frozenset([1, 2, 3])
        # Potenzmenge
        from itertools import combinations
        all_subsets = [frozenset()] + [
            frozenset(c) for r in range(1, 4) for c in combinations([1, 2, 3], r)
        ]
        self.F = SigmaAlgebra(omega, all_subsets)

        # Wahrscheinlichkeitsmaß: P({1})=1/3, P({2})=1/3, P({3})=1/3
        mu_vals = {A: len(A) / 3.0 for A in all_subsets}
        self.prob_measure = Measure(self.F, mu_vals)

        # Nicht-normiertes Maß: μ({i}) = i
        mu_vals2 = {A: sum(A) for A in all_subsets}
        self.measure2 = Measure(self.F, mu_vals2)

    def test_is_measure(self):
        """Zählmaß ist ein gültiges Maß."""
        assert self.prob_measure.is_measure()

    def test_empty_set_measure_zero(self):
        """Maß der leeren Menge ist 0."""
        assert self.prob_measure.measure_of(frozenset()) == 0.0

    def test_probability_measure(self):
        """Wahrscheinlichkeitsmaß: P(Ω) = 1."""
        assert self.prob_measure.is_probability_measure()

    def test_non_probability_measure(self):
        """Nicht-normiertes Maß ist kein Wahrscheinlichkeitsmaß."""
        assert not self.measure2.is_probability_measure()

    def test_sigma_finite(self):
        """Endliches Maß ist σ-endlich."""
        assert self.prob_measure.is_sigma_finite()

    def test_measure_of_singleton(self):
        """Maß einer Einpunktmenge."""
        result = self.prob_measure.measure_of(frozenset([1]))
        assert abs(result - 1.0 / 3.0) < 1e-10

    def test_measure_of_omega(self):
        """Maß der Grundmenge."""
        result = self.prob_measure.measure_of(frozenset([1, 2, 3]))
        assert abs(result - 1.0) < 1e-10

    def test_measure_of_non_measurable_raises(self):
        """Nicht-messbare Menge löst ValueError aus."""
        # Erstelle kleine σ-Algebra ohne {1}
        omega = frozenset([1, 2])
        F_small = SigmaAlgebra(omega, [frozenset(), omega])
        mu_small = Measure(F_small, {omega: 1.0})
        with pytest.raises(ValueError):
            mu_small.measure_of(frozenset([1]))

    def test_outer_measure(self):
        """Äußeres Maß ist kleiner oder gleich echtem Maß."""
        outer = self.prob_measure.outer_measure(frozenset([1]))
        inner = self.prob_measure.measure_of(frozenset([1]))
        assert outer <= inner + 1e-10

    def test_outer_measure_monotone(self):
        """Äußeres Maß ist monoton: A ⊆ B → μ*(A) ≤ μ*(B)."""
        mu_1 = self.prob_measure.outer_measure(frozenset([1]))
        mu_12 = self.prob_measure.outer_measure(frozenset([1, 2]))
        assert mu_1 <= mu_12 + 1e-10


# ---------------------------------------------------------------------------
# LebesgueMeasure Tests
# ---------------------------------------------------------------------------

class TestLebesgueMeasure:
    """Tests für das Lebesgue-Maß."""

    def setup_method(self):
        self.lm = LebesgueMeasure()

    def test_interval_measure(self):
        """λ([0,1]) = 1."""
        assert abs(self.lm.measure_interval(0, 1) - 1.0) < 1e-10

    def test_interval_measure_general(self):
        """λ([a,b]) = b-a."""
        assert abs(self.lm.measure_interval(2.5, 5.7) - 3.2) < 1e-10

    def test_degenerate_interval(self):
        """λ([a,a]) = 0 (Punktmenge)."""
        assert self.lm.measure_interval(3, 3) == 0.0

    def test_invalid_interval_raises(self):
        """a > b löst ValueError aus."""
        with pytest.raises(ValueError):
            self.lm.measure_interval(3, 1)

    def test_union_disjoint_intervals(self):
        """Disjunkte Intervalle: λ([0,1] ∪ [2,3]) = 2."""
        result = self.lm.measure_union_intervals([(0, 1), (2, 3)])
        assert abs(result - 2.0) < 1e-10

    def test_union_overlapping_intervals(self):
        """Überlappende Intervalle korrekt berechnet."""
        result = self.lm.measure_union_intervals([(0, 2), (1, 3)])
        assert abs(result - 3.0) < 1e-10

    def test_union_empty_list(self):
        """Leere Liste hat Maß 0."""
        assert self.lm.measure_union_intervals([]) == 0.0

    def test_cantor_set_measure_zero(self):
        """Cantor-Menge hat Maß 0."""
        result = self.lm.cantor_set_measure()
        assert abs(result) < 1e-6

    def test_fat_cantor_set_measure(self):
        """Fette Cantor-Menge hat Maß 1 - ε."""
        epsilon = 0.3
        result = self.lm.fat_cantor_set_measure(epsilon)
        assert abs(result - (1.0 - epsilon)) < 1e-10

    def test_fat_cantor_invalid_epsilon(self):
        """epsilon außerhalb (0,1) löst ValueError aus."""
        with pytest.raises(ValueError):
            self.lm.fat_cantor_set_measure(1.5)


# ---------------------------------------------------------------------------
# Lebesgue-Integral Tests
# ---------------------------------------------------------------------------

class TestLebesgueIntegral:
    """Tests für das Lebesgue-Integral."""

    def test_constant_function(self):
        """∫₀¹ 2 dλ = 2."""
        result = lebesgue_integral(lambda x: 2.0, 0, 1)
        assert abs(result - 2.0) < 0.01

    def test_linear_function(self):
        """∫₀¹ x dλ = 1/2."""
        result = lebesgue_integral(lambda x: x, 0, 1)
        assert abs(result - 0.5) < 0.01

    def test_quadratic_function(self):
        """∫₀¹ x² dλ = 1/3."""
        result = lebesgue_integral(lambda x: x ** 2, 0, 1)
        assert abs(result - 1.0 / 3.0) < 0.01

    def test_sine_function(self):
        """∫₀^π sin(x) dλ = 2."""
        result = lebesgue_integral(np.sin, 0, np.pi)
        assert abs(result - 2.0) < 0.01

    def test_exponential_function(self):
        """∫₀¹ e^x dλ = e - 1."""
        result = lebesgue_integral(np.exp, 0, 1)
        assert abs(result - (np.e - 1)) < 0.01


# ---------------------------------------------------------------------------
# Riemann vs Lebesgue Tests
# ---------------------------------------------------------------------------

class TestRiemannVsLebesgue:
    """Tests für den Vergleich Riemann- vs. Lebesgue-Integral."""

    def test_dirichlet_function(self):
        """Dirichlet-Funktion: Lebesgue-int.=0, nicht Riemann-int."""
        result = riemann_vs_lebesgue('dirichlet')
        assert result['lebesgue_integrable'] is True
        assert result['riemann_integrable'] is False
        assert result['lebesgue_integral'] == 0.0

    def test_bounded_discontinuous(self):
        """Beschränkte Funktion mit Unstetigkeitsstelle."""
        result = riemann_vs_lebesgue('bounded_discontinuous')
        assert result['riemann_integrable'] is True
        assert result['lebesgue_integrable'] is True

    def test_improper_integral(self):
        """Uneigentliches Integral."""
        result = riemann_vs_lebesgue('improper')
        assert result['lebesgue_integrable'] is True
        assert abs(result['riemann_integral'] - 2.0) < 1e-10
        assert abs(result['lebesgue_integral'] - 2.0) < 0.1

    def test_invalid_type_raises(self):
        """Unbekannter Typ löst ValueError aus."""
        with pytest.raises(ValueError):
            riemann_vs_lebesgue('unknown_type')


# ---------------------------------------------------------------------------
# MeasurableFunction Tests
# ---------------------------------------------------------------------------

class TestMeasurableFunction:
    """Tests für messbare Funktionen."""

    def setup_method(self):
        omega = frozenset([1, 2, 3])
        from itertools import combinations
        all_subsets = [frozenset()] + [
            frozenset(c) for r in range(1, 4) for c in combinations([1, 2, 3], r)
        ]
        self.F = SigmaAlgebra(omega, all_subsets)
        mu_vals = {A: len(A) / 3.0 for A in all_subsets}
        self.mu = Measure(self.F, mu_vals)

    def test_constant_function_integral(self):
        """Integral einer konstanten Funktion f ≡ 2: ∫f dμ = 2·μ(Ω) = 2."""
        values = {1: 2.0, 2: 2.0, 3: 2.0}
        f = MeasurableFunction(self.F, values)
        result = f.integral(self.mu)
        assert abs(result - 2.0) < 1e-10

    def test_lp_norm_l1(self):
        """L¹-Norm für f ≡ 1 auf [0,1] = 1."""
        x_vals = list(np.linspace(0, 1, 100))
        f_vals = [1.0] * 100
        lp = LpSpace(1, self.mu)
        result = lp.norm(f_vals, x_vals)
        assert abs(result - 1.0) < 0.01

    def test_lp_norm_l2(self):
        """L²-Norm für f(x)=x auf [0,1] = 1/√3."""
        x_vals = list(np.linspace(0, 1, 10000))
        f_vals = list(np.array(x_vals))
        lp = LpSpace(2, self.mu)
        result = lp.norm(f_vals, x_vals)
        expected = 1.0 / np.sqrt(3)
        assert abs(result - expected) < 0.01

    def test_lp_norm_linf(self):
        """L∞-Norm = Maximum."""
        x_vals = list(np.linspace(0, 1, 100))
        f_vals = list(np.linspace(0, 3, 100))
        lp = LpSpace(float('inf'), self.mu)
        result = lp.norm(f_vals, x_vals)
        assert abs(result - 3.0) < 0.1


# ---------------------------------------------------------------------------
# Konvergenzsätze Tests
# ---------------------------------------------------------------------------

class TestConvergenceTheorems:
    """Tests für die großen Konvergenzsätze."""

    def test_monotone_convergence_theorem(self):
        """Satz von der monotonen Konvergenz: Integrale konvergieren."""
        result = monotone_convergence_theorem_demo()
        assert result['monotone_increasing'] is True
        assert result['converges_to_limit'] is True
        assert abs(result['exact_limit_integral'] - 1.0 / 3.0) < 1e-10

    def test_monotone_convergence_values(self):
        """Konvergenzwerte werden kleiner."""
        result = monotone_convergence_theorem_demo()
        convergence = result['convergence']
        # Letzte Werte kleiner als erste
        assert convergence[-1] < convergence[0]

    def test_dominated_convergence_theorem(self):
        """Satz von der dominierten Konvergenz: Konvergenz zum Nullintegral."""
        result = dominated_convergence_theorem_demo()
        assert result['converges'] is True
        assert result['limit_integral'] == 0.0

    def test_dominated_convergence_monotone_decrease(self):
        """Integrale nehmen ab."""
        result = dominated_convergence_theorem_demo()
        vals = result['convergence_to_zero']
        # Allgemein abnehmend (nicht streng monoton wegen Vorzeichen)
        assert vals[-1] < vals[0]

    def test_fubini_theorem(self):
        """Satz von Fubini: alle drei Methoden stimmen überein."""
        result = fubini_theorem_demo()
        assert result['fubini_holds'] is True
        assert abs(result['exact_value'] - 2.0 / 3.0) < 1e-10

    def test_fubini_accuracy(self):
        """Fubini-Näherung mit Fehler < 0.01."""
        result = fubini_theorem_demo()
        assert result['error_double'] < 0.01
        assert result['error_iterated_xy'] < 0.01
        assert result['error_iterated_yx'] < 0.01


# ---------------------------------------------------------------------------
# Produktmaß Tests
# ---------------------------------------------------------------------------

class TestProductMeasure:
    """Tests für das Produktmaß."""

    def setup_method(self):
        """Erstellt zwei einfache Maßräume."""
        omega1 = frozenset(['a', 'b'])
        F1 = SigmaAlgebra(omega1, [frozenset(), frozenset(['a']), frozenset(['b']), omega1])
        mu1_vals = {
            frozenset(): 0.0,
            frozenset(['a']): 0.4,
            frozenset(['b']): 0.6,
            omega1: 1.0,
        }
        self.mu1 = Measure(F1, mu1_vals)

        omega2 = frozenset([1, 2])
        F2 = SigmaAlgebra(omega2, [frozenset(), frozenset([1]), frozenset([2]), omega2])
        mu2_vals = {
            frozenset(): 0.0,
            frozenset([1]): 0.3,
            frozenset([2]): 0.7,
            omega2: 1.0,
        }
        self.mu2 = Measure(F2, mu2_vals)

    def test_product_measure_total(self):
        """Produktmaß: Gesamtmaß = μ₁(Ω₁) · μ₂(Ω₂) = 1."""
        result = product_measure(self.mu1, self.mu2)
        assert abs(result['total_measure'] - 1.0) < 1e-10

    def test_product_measure_rectangle(self):
        """Produktmaß eines Rechtecks: (μ₁⊗μ₂)(A×B) = μ₁(A)·μ₂(B)."""
        result = product_measure(self.mu1, self.mu2)
        A = frozenset(['a'])
        B = frozenset([1])
        val = result['product_measure_values'].get((A, B), None)
        assert val is not None
        assert abs(val - 0.4 * 0.3) < 1e-10

    def test_product_measure_empty_set(self):
        """Produktmaß der leeren Menge ist 0."""
        result = product_measure(self.mu1, self.mu2)
        val = result['product_measure_values'].get((frozenset(), frozenset()), 0.0)
        assert val == 0.0


# ---------------------------------------------------------------------------
# L^p-Raum Tests
# ---------------------------------------------------------------------------

class TestLpSpace:
    """Tests für L^p-Räume."""

    def setup_method(self):
        omega = frozenset([1])
        F = SigmaAlgebra(omega, [frozenset(), omega])
        mu = Measure(F, {omega: 1.0})
        self.mu = mu
        self.x = list(np.linspace(0, 1, 1000))

    def test_l2_inner_product(self):
        """L²-Skalarprodukt: ⟨f,g⟩ = ∫f·g dx."""
        lp = LpSpace(2, self.mu)
        f = list(np.ones(1000))
        g = list(np.ones(1000))
        result = lp.inner_product(f, g, self.x)
        assert abs(result - 1.0) < 0.01

    def test_l2_inner_product_orthogonal(self):
        """sin und cos sind orthogonal auf [0,π]."""
        x = list(np.linspace(0, np.pi, 10000))
        f = list(np.sin(np.array(x)))
        g = list(np.cos(np.array(x)))
        lp = LpSpace(2, self.mu)
        result = lp.inner_product(f, g, x)
        assert abs(result) < 0.01

    def test_inner_product_only_for_l2(self):
        """Skalarprodukt nur für p=2 definiert."""
        lp = LpSpace(3, self.mu)
        with pytest.raises(ValueError):
            lp.inner_product([1, 2], [1, 2], [0, 1])

    def test_holder_inequality(self):
        """Hölder-Ungleichung: ∫|fg| ≤ ‖f‖_p · ‖g‖_q."""
        lp = LpSpace(2, self.mu)
        x = list(np.linspace(0, 1, 1000))
        f = list(np.array(x) ** 2)
        g = list(np.sqrt(np.array(x)))
        assert lp.holder_inequality_check(f, g, x, 2.0, 2.0)

    def test_minkowski_inequality(self):
        """Minkowski-Ungleichung: ‖f+g‖_p ≤ ‖f‖_p + ‖g‖_p."""
        lp = LpSpace(2, self.mu)
        x = list(np.linspace(0, 1, 1000))
        f = list(np.sin(np.pi * np.array(x)))
        g = list(np.cos(np.pi * np.array(x)))
        assert lp.minkowski_inequality_check(f, g, x)

    def test_invalid_p_raises(self):
        """p < 1 löst ValueError aus."""
        with pytest.raises(ValueError):
            LpSpace(0.5, self.mu)

    def test_l1_norm(self):
        """L¹-Norm = ∫|f| dx."""
        lp = LpSpace(1, self.mu)
        x = list(np.linspace(0, 1, 10000))
        f = list(np.ones(10000))
        result = lp.norm(f, x)
        assert abs(result - 1.0) < 0.01


# ---------------------------------------------------------------------------
# Radon-Nikodym Tests
# ---------------------------------------------------------------------------

class TestRadonNikodym:
    """Tests für den Radon-Nikodym-Satz."""

    def test_radon_nikodym_all_verified(self):
        """Alle Verifikationen von ν(A) = ∫_A f dμ sind korrekt."""
        result = radon_nikodym_theorem_demo()
        assert result['all_verified'] is True

    def test_radon_nikodym_density_values(self):
        """Radon-Nikodym-Dichte hat korrekte Werte."""
        result = radon_nikodym_theorem_demo()
        density = result['density_values']
        assert abs(density[1] - 0.5) < 1e-10
        assert abs(density[2] - 1.0) < 1e-10
        assert abs(density[3] - 1.5) < 1e-10

    def test_radon_nikodym_empty_set(self):
        """ν(∅) = 0."""
        result = radon_nikodym_theorem_demo()
        empty_key = tuple(sorted(frozenset()))
        empty_data = result['verification'][empty_key]
        assert abs(empty_data['nu(A)']) < 1e-10
        assert abs(empty_data['integral_f_dmu']) < 1e-10


# ---------------------------------------------------------------------------
# Cantor-Funktion Tests
# ---------------------------------------------------------------------------

class TestCantorFunction:
    """Tests für die Cantor-Funktion (Teufelsleiter)."""

    def test_cantor_function_boundary_values(self):
        """c(0) = 0 und c(1) = 1."""
        result = cantor_function_demo()
        assert abs(result['sample_values']['0'] - 0.0) < 1e-10
        assert abs(result['sample_values']['1'] - 1.0) < 1e-10

    def test_cantor_function_half(self):
        """c(1/2) = 1/2 (Mittelpunkt)."""
        result = cantor_function_demo()
        assert abs(result['sample_values']['1/2'] - 0.5) < 1e-3

    def test_cantor_function_singular(self):
        """Cantor-Funktion ist singulär."""
        result = cantor_function_demo()
        assert result['is_singular'] is True

    def test_cantor_function_zero_derivative(self):
        """Ableitung ist fast überall 0 (≥ 80% der Punkte bei diskreter Approximation)."""
        result = cantor_function_demo()
        assert result['fraction_with_zero_derivative'] > 0.8

    def test_cantor_function_properties_list(self):
        """Eigenschaften-Liste ist nicht leer."""
        result = cantor_function_demo()
        assert len(result['properties']) >= 4

    def test_cantor_function_one_third(self):
        """c(1/3) = 1/2 (Ende des ersten entfernten Intervalls)."""
        result = cantor_function_demo()
        # c(1/3) sollte 1/2 sein (Cantor-Funktion ist auf [1/3,2/3] konstant = 1/2)
        assert abs(result['sample_values']['1/3'] - 0.5) < 1e-3

    def test_cantor_function_two_thirds(self):
        """c(2/3) = 1/2 (Anfang des ersten entfernten Intervalls)."""
        result = cantor_function_demo()
        assert abs(result['sample_values']['2/3'] - 0.5) < 1e-3
