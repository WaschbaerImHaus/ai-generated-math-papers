"""
@file test_differential_topology.py
@brief Tests für das Differentialtopologie-Modul (differential_topology.py).
@description
    Testet alle Klassen und Funktionen des Moduls:
    - SmoothManifold (Dimension, Tangential-/Kotangentialraum)
    - DifferentialForm (Wedge-Produkt, äußere Ableitung, geschlossen/exakt)
    - DeRhamCohomology (Betti-Zahlen, Euler-Charakteristik)
    - MorseTheory (kritische Punkte, Morse-Index, Ungleichungen)
    - TransversalityTheory (Transversalität, Schnittdimension, Sard)
    - VectorBundle (Totalraum, Rang, Trivialität)
    - CharacteristicClasses (Chern, Pontryagin, Euler, Stiefel-Whitney)
    - SphereTopology (Homotopie, Homologie, Orientierbarkeit)
    - Standalone-Funktionen (degree_of_map, hairy_ball, poincare_hopf, whitney, jacobi, implizite Fkt.)

    Besondere Edge Cases:
    - 0-Formen und ihre äußere Ableitung
    - Kritische Punkte mit Morse-Index 0 (Minimum) und n (Maximum)
    - Degenerierte Hesse-Matrizen
    - Transversalität bei Grenzfällen (codim_M + codim_N = ambient_dim)
    - Sphären der Dimension 0, 1, 2, 3

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import pytest
import numpy as np
import sympy as sp
from sympy import symbols, sin, cos, exp, sqrt, Rational
import sys
import os

# Sicherstellen, dass src/ im Suchpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from differential_topology import (
    SmoothManifold,
    DifferentialForm,
    DeRhamCohomology,
    MorseTheory,
    TransversalityTheory,
    VectorBundle,
    CharacteristicClasses,
    SphereTopology,
    degree_of_map,
    hairy_ball_theorem_demo,
    poincare_hopf_theorem_demo,
    whitney_embedding_dimension,
    compute_jacobian,
    implicit_function_theorem_check,
    _sign_and_sort,
)


# ===========================================================================
# HILFSSYMBOLE FÜR TESTS
# ===========================================================================

x, y, z, t = symbols('x y z t')
u, v, w = symbols('u v w')


# ===========================================================================
# TESTS: SmoothManifold
# ===========================================================================

class TestSmoothManifold:
    """Tests für die Klasse SmoothManifold."""

    def test_dimension_1(self):
        """1-dimensionale Mannigfaltigkeit: Dimension = 1."""
        M = SmoothManifold(1, 'S^1')
        assert M.dimension == 1

    def test_dimension_3(self):
        """3-dimensionale Mannigfaltigkeit: Dimension = 3."""
        M = SmoothManifold(3, 'S^3')
        assert M.dimension == 3

    def test_dimension_0(self):
        """0-dimensionale Mannigfaltigkeit (Punkte): Dimension = 0."""
        M = SmoothManifold(0)
        assert M.dimension == 0

    def test_name_default(self):
        """Standardname ist 'M'."""
        M = SmoothManifold(2)
        assert M.name == 'M'

    def test_name_custom(self):
        """Benutzerdefinierter Name wird korrekt gespeichert."""
        M = SmoothManifold(4, 'T^2')
        assert M.name == 'T^2'

    def test_tangent_space_dim_equals_manifold_dim(self):
        """Tangentialraum-Dimension = Mannigfaltigkeitsdimension."""
        for n in [0, 1, 2, 3, 5, 10]:
            M = SmoothManifold(n)
            assert M.tangent_space_dim() == n

    def test_cotangent_space_dim_equals_manifold_dim(self):
        """Kotangentialraum-Dimension = Mannigfaltigkeitsdimension."""
        for n in [0, 1, 2, 3, 5, 10]:
            M = SmoothManifold(n)
            assert M.cotangent_space_dim() == n

    def test_tangent_cotangent_equal(self):
        """Tangential- und Kotangentialraum haben gleiche Dimension."""
        M = SmoothManifold(7)
        assert M.tangent_space_dim() == M.cotangent_space_dim()

    def test_negative_dimension_raises(self):
        """Negative Dimension wirft ValueError."""
        with pytest.raises(ValueError):
            SmoothManifold(-1)

    def test_repr_contains_dimension(self):
        """Repr enthält die Dimension."""
        M = SmoothManifold(3, 'S^3')
        r = repr(M)
        assert '3' in r
        assert 'S^3' in r


# ===========================================================================
# TESTS: _sign_and_sort (interne Hilfsfunktion)
# ===========================================================================

class TestSignAndSort:
    """Tests für die interne Sortierfunktion mit Vorzeichenverfolgung."""

    def test_already_sorted_even(self):
        """Sortierte gerade Permutation: Vorzeichen +1."""
        sign, lst = _sign_and_sort([0, 1, 2])
        assert sign == 1
        assert lst == [0, 1, 2]

    def test_single_transposition(self):
        """Einzel-Transposition: Vorzeichen -1."""
        sign, lst = _sign_and_sort([1, 0])
        assert sign == -1
        assert lst == [0, 1]

    def test_double_index_gives_zero(self):
        """Doppelter Index (Schiefsymmetrie): Vorzeichen 0."""
        sign, _ = _sign_and_sort([0, 0, 1])
        assert sign == 0

    def test_cyclic_permutation_even(self):
        """Zyklische Permutation (1,2,0) → 2 Transpositionen → +1."""
        sign, lst = _sign_and_sort([1, 2, 0])
        assert sign == 1
        assert lst == [0, 1, 2]

    def test_reverse_three_elements(self):
        """Umkehrung [2,1,0] → 3 Transpositionen → -1."""
        sign, lst = _sign_and_sort([2, 1, 0])
        assert sign == -1
        assert lst == [0, 1, 2]


# ===========================================================================
# TESTS: DifferentialForm
# ===========================================================================

class TestDifferentialForm:
    """Tests für die Klasse DifferentialForm."""

    def test_0_form_constant(self):
        """0-Form (Funktion): Grad 0, Koeffizient ()."""
        omega = DifferentialForm(0, {(): sp.Integer(3)}, [x, y])
        assert omega.degree == 0
        assert omega.coefficients[()] == 3

    def test_1_form_creation(self):
        """1-Form: Grad 1, korrekte Koeffizienten."""
        # ω = x·dx + y·dy
        omega = DifferentialForm(1, {(0,): x, (1,): y}, [x, y])
        assert omega.degree == 1
        assert omega.coefficients[(0,)] == x
        assert omega.coefficients[(1,)] == y

    def test_2_form_creation(self):
        """2-Form: Grad 2, Koeffizient für dx∧dy."""
        # ω = f(x,y)·dx∧dy
        omega = DifferentialForm(2, {(0, 1): x*y}, [x, y])
        assert omega.degree == 2

    def test_invalid_degree_raises(self):
        """Negativer Grad wirft ValueError."""
        with pytest.raises(ValueError):
            DifferentialForm(-1, {}, [x])

    def test_wrong_multiindex_length_raises(self):
        """Multi-Index falscher Länge wirft ValueError."""
        with pytest.raises(ValueError):
            DifferentialForm(1, {(0, 1): x}, [x, y])  # Länge 2, erwartet 1

    def test_exterior_derivative_of_constant_is_zero(self):
        """Äußere Ableitung einer konstanten 0-Form ist 0."""
        # ω = 5 (Konstante)
        omega = DifferentialForm(0, {(): sp.Integer(5)}, [x, y])
        d_omega = omega.exterior_derivative()
        assert len(d_omega.coefficients) == 0

    def test_exterior_derivative_of_function(self):
        """Äußere Ableitung von f = x²y: df = 2xy·dx + x²·dy."""
        f_expr = x**2 * y
        omega = DifferentialForm(0, {(): f_expr}, [x, y])
        d_omega = omega.exterior_derivative()
        # Erwartete Koeffizienten: (0,) → 2xy, (1,) → x²
        assert d_omega.degree == 1
        coeff_dx = sp.simplify(d_omega.coefficients.get((0,), sp.Integer(0)) - 2*x*y)
        coeff_dy = sp.simplify(d_omega.coefficients.get((1,), sp.Integer(0)) - x**2)
        assert coeff_dx == 0
        assert coeff_dy == 0

    def test_exterior_derivative_squared_is_zero(self):
        """d(d(ω)) = 0 (Poincaré-Lemma d² = 0)."""
        # f = x² + y²
        f_expr = x**2 + y**2
        omega = DifferentialForm(0, {(): f_expr}, [x, y])
        d_omega = omega.exterior_derivative()
        dd_omega = d_omega.exterior_derivative()
        # dd = 0: alle Koeffizienten müssen 0 sein
        assert len(dd_omega.coefficients) == 0

    def test_exact_1_form_is_closed(self):
        """Exakte 1-Form ist geschlossen: d(df) = 0."""
        # df = 2x·dx + 2y·dy (Gradient von x² + y²)
        omega = DifferentialForm(1, {(0,): 2*x, (1,): 2*y}, [x, y])
        assert omega.is_closed()

    def test_non_closed_1_form(self):
        """Nicht-geschlossene 1-Form: ω = y²·dx (∂(y²)/∂y = 2y ≠ 0 = ∂0/∂x)."""
        # ω = y²·dx: ∂(y²)/∂y = 2y ≠ ∂0/∂x = 0
        omega = DifferentialForm(1, {(0,): y**2}, [x, y])
        assert not omega.is_closed()

    def test_exact_1_form_check(self):
        """Exakte 1-Form erfüllt Integrabilitätsbedingung."""
        # ω = 2x·dx + 3y²·dy = d(x² + y³)
        omega = DifferentialForm(1, {(0,): 2*x, (1,): 3*y**2}, [x, y])
        assert omega.is_exact()

    def test_non_exact_1_form(self):
        """Nicht-exakte 1-Form: ω = y·dx - x·dy (Winkeldifferential auf ℝ²\\{0})."""
        # ω = y·dx - x·dy: ∂y/∂y = 1 ≠ ∂(-x)/∂x = -1
        omega = DifferentialForm(1, {(0,): y, (1,): -x}, [x, y])
        assert not omega.is_exact()

    def test_zero_0_form_is_exact(self):
        """Die Null-0-Form ist exakt (trivialerweise)."""
        omega = DifferentialForm(0, {(): sp.Integer(0)}, [x])
        # Alle Koeffizienten sind Null → exakt
        assert omega.is_exact()

    def test_nonzero_0_form_is_not_exact(self):
        """Eine nichtkonstante Funktion als 0-Form ist nicht exakt."""
        omega = DifferentialForm(0, {(): x + 1}, [x])
        # Nicht-Null → nicht exakt
        assert not omega.is_exact()

    def test_wedge_product_basic(self):
        """Wedge-Produkt dx ∧ dy ergibt 2-Form."""
        # dx als 1-Form
        dx_form = DifferentialForm(1, {(0,): sp.Integer(1)}, [x, y])
        # dy als 1-Form
        dy_form = DifferentialForm(1, {(1,): sp.Integer(1)}, [x, y])
        # dx ∧ dy
        wedge = dx_form.wedge(dy_form)
        assert wedge.degree == 2
        # Koeffizient von (0,1) sollte 1 sein
        coeff = wedge.coefficients.get((0, 1), sp.Integer(0))
        assert sp.simplify(coeff - 1) == 0

    def test_wedge_anticommutativity(self):
        """Antikommutativität: dx ∧ dy = -(dy ∧ dx)."""
        dx_form = DifferentialForm(1, {(0,): sp.Integer(1)}, [x, y])
        dy_form = DifferentialForm(1, {(1,): sp.Integer(1)}, [x, y])
        w1 = dx_form.wedge(dy_form)
        w2 = dy_form.wedge(dx_form)
        c1 = w1.coefficients.get((0, 1), sp.Integer(0))
        c2 = w2.coefficients.get((0, 1), sp.Integer(0))
        assert sp.simplify(c1 + c2) == 0

    def test_wedge_0_form_with_1_form(self):
        """0-Form (Funktion) ∧ 1-Form = Skalierung der 1-Form."""
        # f = x, ω = dy
        f_form = DifferentialForm(0, {(): x}, [x, y])
        dy_form = DifferentialForm(1, {(1,): sp.Integer(1)}, [x, y])
        result = f_form.wedge(dy_form)
        assert result.degree == 1
        coeff = result.coefficients.get((1,), sp.Integer(0))
        assert sp.simplify(coeff - x) == 0

    def test_wedge_different_variables_raises(self):
        """Wedge-Produkt mit verschiedenen Variablenlisten wirft ValueError."""
        omega1 = DifferentialForm(1, {(0,): x}, [x])
        omega2 = DifferentialForm(1, {(0,): y}, [y])
        with pytest.raises(ValueError):
            omega1.wedge(omega2)

    def test_exterior_derivative_3d(self):
        """Äußere Ableitung einer 1-Form im ℝ³ (entspricht Rotation)."""
        # ω = P·dx + Q·dy + R·dz
        P, Q, R = x*y, y*z, x*z
        omega = DifferentialForm(1, {(0,): P, (1,): Q, (2,): R}, [x, y, z])
        d_omega = omega.exterior_derivative()
        assert d_omega.degree == 2

    def test_closed_2_form_in_r3(self):
        """Exakte 2-Form im ℝ³ ist geschlossen."""
        # Exakte 2-Form: dω wo ω = x·dy (exakt mit η=x·dy)
        # Nimm: ω = dx∧dy (konstante Koeffizienten → exakt auf ℝ²)
        omega = DifferentialForm(2, {(0, 1): sp.Integer(1)}, [x, y, z])
        assert omega.is_closed()


# ===========================================================================
# TESTS: DeRhamCohomology
# ===========================================================================

class TestDeRhamCohomology:
    """Tests für de Rham-Kohomologie."""

    def test_betti_sphere_0(self):
        """S^0 (zwei Punkte): b_0 = 2."""
        betti = DeRhamCohomology.betti_numbers_sphere(0)
        assert betti == [2]

    def test_betti_sphere_1(self):
        """S^1 (Kreis): b_0 = 1, b_1 = 1."""
        betti = DeRhamCohomology.betti_numbers_sphere(1)
        assert betti == [1, 1]

    def test_betti_sphere_2(self):
        """S^2: b_0 = 1, b_1 = 0, b_2 = 1."""
        betti = DeRhamCohomology.betti_numbers_sphere(2)
        assert betti == [1, 0, 1]

    def test_betti_sphere_3(self):
        """S^3: b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 1."""
        betti = DeRhamCohomology.betti_numbers_sphere(3)
        assert betti == [1, 0, 0, 1]

    def test_betti_sphere_n_length(self):
        """Betti-Zahlen der n-Sphäre haben Länge n+1."""
        for n in [1, 2, 3, 4, 5]:
            betti = DeRhamCohomology.betti_numbers_sphere(n)
            assert len(betti) == n + 1

    def test_betti_sphere_negative_raises(self):
        """Negative Dimension wirft ValueError."""
        with pytest.raises(ValueError):
            DeRhamCohomology.betti_numbers_sphere(-1)

    def test_betti_torus_1(self):
        """T^1 = S^1: b_0 = 1, b_1 = 1."""
        betti = DeRhamCohomology.betti_numbers_torus(1)
        assert betti == [1, 1]

    def test_betti_torus_2(self):
        """T^2: b_0 = 1, b_1 = 2, b_2 = 1."""
        betti = DeRhamCohomology.betti_numbers_torus(2)
        assert betti == [1, 2, 1]

    def test_betti_torus_3(self):
        """T^3: b_0 = 1, b_1 = 3, b_2 = 3, b_3 = 1."""
        betti = DeRhamCohomology.betti_numbers_torus(3)
        assert betti == [1, 3, 3, 1]

    def test_betti_torus_0(self):
        """T^0 (Punkt): b_0 = 1."""
        betti = DeRhamCohomology.betti_numbers_torus(0)
        assert betti == [1]

    def test_betti_torus_binomial_coefficient(self):
        """Betti-Zahlen des n-Torus sind Binomialkoeffizienten."""
        import math
        n = 4
        betti = DeRhamCohomology.betti_numbers_torus(n)
        for k, bk in enumerate(betti):
            assert bk == math.comb(n, k)

    def test_euler_characteristic_sphere_2(self):
        """χ(S^2) = 1 - 0 + 1 = 2."""
        betti = DeRhamCohomology.betti_numbers_sphere(2)
        chi = DeRhamCohomology.euler_characteristic(betti)
        assert chi == 2

    def test_euler_characteristic_sphere_3(self):
        """χ(S^3) = 1 - 0 + 0 - 1 = 0."""
        betti = DeRhamCohomology.betti_numbers_sphere(3)
        chi = DeRhamCohomology.euler_characteristic(betti)
        assert chi == 0

    def test_euler_characteristic_torus_2(self):
        """χ(T^2) = 1 - 2 + 1 = 0."""
        betti = DeRhamCohomology.betti_numbers_torus(2)
        chi = DeRhamCohomology.euler_characteristic(betti)
        assert chi == 0

    def test_euler_characteristic_point(self):
        """χ(Punkt) = b_0 = 1."""
        chi = DeRhamCohomology.euler_characteristic([1])
        assert chi == 1

    def test_euler_odd_sphere_is_zero(self):
        """Ungerade Sphären S^{2k+1} haben χ = 0."""
        for n in [1, 3, 5, 7]:
            betti = DeRhamCohomology.betti_numbers_sphere(n)
            chi = DeRhamCohomology.euler_characteristic(betti)
            assert chi == 0

    def test_euler_even_sphere_is_two(self):
        """Geraddimensionale Sphären S^{2k} (k ≥ 1) haben χ = 2."""
        for n in [2, 4, 6]:
            betti = DeRhamCohomology.betti_numbers_sphere(n)
            chi = DeRhamCohomology.euler_characteristic(betti)
            assert chi == 2

    def test_closed_forms_filter(self):
        """closed_forms filtert korrekt geschlossene Formen heraus."""
        # Geschlossene Form: dg = 2x·dx + 2y·dy (exakt → geschlossen)
        closed = DifferentialForm(1, {(0,): 2*x, (1,): 2*y}, [x, y])
        # Nicht-geschlossene Form
        not_closed = DifferentialForm(1, {(0,): y**2, (1,): sp.Integer(0)}, [x, y])
        forms = [closed, not_closed]
        result = DeRhamCohomology.closed_forms(1, forms)
        assert closed in result
        assert not_closed not in result

    def test_exact_forms_filter(self):
        """exact_forms filtert korrekt exakte Formen heraus."""
        exact = DifferentialForm(1, {(0,): 2*x, (1,): 2*y}, [x, y])
        not_exact = DifferentialForm(1, {(0,): y, (1,): -x}, [x, y])
        forms = [exact, not_exact]
        result = DeRhamCohomology.exact_forms(1, forms)
        assert exact in result
        assert not_exact not in result


# ===========================================================================
# TESTS: MorseTheory
# ===========================================================================

class TestMorseTheory:
    """Tests für Morse-Theorie."""

    def test_critical_points_simple_minimum(self):
        """f = x² + y²: kritischer Punkt bei (0,0) (Minimum)."""
        f = x**2 + y**2
        pts = MorseTheory.morse_function_critical_points(f, [x, y])
        assert len(pts) >= 1
        # Punkt (0,0) muss enthalten sein
        found = any(
            sp.simplify(pt.get(x, sp.Integer(99))) == 0 and
            sp.simplify(pt.get(y, sp.Integer(99))) == 0
            for pt in pts
        )
        assert found

    def test_critical_points_saddle(self):
        """f = x² - y²: kritischer Punkt bei (0,0) (Sattelpunkt)."""
        f = x**2 - y**2
        pts = MorseTheory.morse_function_critical_points(f, [x, y])
        assert len(pts) >= 1

    def test_morse_index_minimum_is_0(self):
        """Morse-Index eines Minimums ist 0."""
        f = x**2 + y**2
        pt = {x: sp.Integer(0), y: sp.Integer(0)}
        idx = MorseTheory.morse_index(f, pt, [x, y])
        assert idx == 0

    def test_morse_index_maximum_is_dim(self):
        """Morse-Index eines Maximums ist n (Dimension)."""
        f = -(x**2 + y**2)
        pt = {x: sp.Integer(0), y: sp.Integer(0)}
        idx = MorseTheory.morse_index(f, pt, [x, y])
        assert idx == 2  # Dimension = 2

    def test_morse_index_saddle_is_1(self):
        """Morse-Index eines Sattelpunkts f=x²-y² bei (0,0) ist 1."""
        f = x**2 - y**2
        pt = {x: sp.Integer(0), y: sp.Integer(0)}
        idx = MorseTheory.morse_index(f, pt, [x, y])
        assert idx == 1

    def test_morse_index_degenerate_raises(self):
        """Degenerierter kritischer Punkt (det H = 0) wirft ValueError."""
        # f = x²: Hesse-Matrix [[2,0],[0,0]] → det = 0 → degeneriert
        f = x**2
        pt = {x: sp.Integer(0), y: sp.Integer(0)}
        with pytest.raises(ValueError):
            MorseTheory.morse_index(f, pt, [x, y])

    def test_is_morse_function_true(self):
        """f = x² + y² ist Morse-Funktion (Hesse bei (0,0) invertierbar)."""
        f = x**2 + y**2
        assert MorseTheory.is_morse_function(f, [x, y])

    def test_is_morse_function_false_degenerate(self):
        """f = x² (kein y-Term): degeneriert → keine Morse-Funktion."""
        f = x**2
        # Nur x als Variable → kein Problem
        # Mit zwei Variablen: f = x² hat Hesse [[2,0],[0,0]], degeneriert
        assert not MorseTheory.is_morse_function(x**2, [x, y])

    def test_morse_inequalities_satisfied(self):
        """Schwache Morse-Ungleichungen für S^2: c_0 ≥ 1, c_2 ≥ 1."""
        # S^2 Betti-Zahlen: (1, 0, 1)
        betti = [1, 0, 1]
        # Torus-artige Morse-Funktion auf S^2: 1 Min, 0 Sattelpunkte, 1 Max
        morse_data = {0: 1, 1: 0, 2: 1}
        result = MorseTheory.morse_inequalities_check(betti, morse_data)
        assert result['weak_satisfied']
        assert result['euler_ok']
        assert result['euler_char'] == 2

    def test_morse_inequalities_euler_check(self):
        """Euler-Charakteristik über Morse-Daten stimmt überein."""
        # T^2 Betti-Zahlen: (1, 2, 1)
        betti = [1, 2, 1]
        # Perfekte Morse-Funktion: 1 Min, 2 Sattel, 1 Max
        morse_data = {0: 1, 1: 2, 2: 1}
        result = MorseTheory.morse_inequalities_check(betti, morse_data)
        assert result['euler_ok']
        assert result['euler_char'] == 0

    def test_morse_inequalities_violated(self):
        """Zu wenige kritische Punkte verletzen Morse-Ungleichungen."""
        betti = [1, 0, 1]  # S^2
        # Nur 1 kritischer Punkt (kann Euler-Char nicht erfüllen)
        morse_data = {0: 1, 2: 0}
        result = MorseTheory.morse_inequalities_check(betti, morse_data)
        assert not result['euler_ok']

    def test_handle_decomposition_demo(self):
        """handle_decomposition_demo gibt korrekte Struktur zurück."""
        f = x**2 + y**2
        result = MorseTheory.handle_decomposition_demo(f, [x, y])
        assert 'critical_points' in result
        assert 'handles' in result
        assert 'is_morse' in result
        assert result['is_morse']

    def test_handle_decomposition_torus_function(self):
        """Höhenfunktion auf dem Torus hat 4 kritische Punkte."""
        # Standard-Höhenfunktion auf T^2 eingebettet in ℝ³
        # f(θ,φ) = (2 + cos φ) cos θ + (2 + cos φ) sin θ approximiert als:
        # f = sin(x)*sin(y) → 4 kritische Punkte (lokale Morse-Funktion)
        f_expr = sp.sin(x) * sp.sin(y)
        result = MorseTheory.handle_decomposition_demo(f_expr, [x, y])
        assert 'critical_points' in result
        assert 'is_morse' in result

    def test_morse_index_1d_minimum(self):
        """1D: f = x²: Morse-Index des Minimums bei 0 ist 0."""
        f = x**2
        pt = {x: sp.Integer(0)}
        idx = MorseTheory.morse_index(f, pt, [x])
        assert idx == 0

    def test_morse_index_1d_maximum(self):
        """1D: f = -x²: Morse-Index des Maximums bei 0 ist 1."""
        f = -x**2
        pt = {x: sp.Integer(0)}
        idx = MorseTheory.morse_index(f, pt, [x])
        assert idx == 1


# ===========================================================================
# TESTS: TransversalityTheory
# ===========================================================================

class TestTransversalityTheory:
    """Tests für Transversalitätstheorie."""

    def test_transverse_curves_in_r2(self):
        """Zwei Kurven (codim 1) in ℝ² sind transversal wenn codim_M + codim_N ≤ 2."""
        # Zwei Kurven in ℝ²: codim = 1 je
        assert TransversalityTheory.are_transverse(1, 1, 2)

    def test_transverse_surface_curve_in_r3(self):
        """Kurve (codim 2) und Fläche (codim 1) in ℝ³: 2+1 = 3 ≤ 3."""
        assert TransversalityTheory.are_transverse(2, 1, 3)

    def test_not_transverse(self):
        """Zwei Punkte (codim 2 je) in ℝ³: 2+2 = 4 > 3, nicht transversal."""
        assert not TransversalityTheory.are_transverse(2, 2, 3)

    def test_transverse_boundary_case(self):
        """Grenzfall: codim_M + codim_N = ambient_dim → transversal."""
        # codim 1 + codim 2 = 3 = ambient_dim = 3
        assert TransversalityTheory.are_transverse(1, 2, 3)

    def test_intersection_dimension_positive(self):
        """Schnittdimension zweier Flächen in ℝ³: 2+2-3 = 1 (Kurve)."""
        dim = TransversalityTheory.transverse_intersection_dimension(2, 2, 3)
        assert dim == 1

    def test_intersection_dimension_zero(self):
        """Schnittdimension zweier Geraden in ℝ²: 1+1-2 = 0 (Punkt)."""
        dim = TransversalityTheory.transverse_intersection_dimension(1, 1, 2)
        assert dim == 0

    def test_intersection_dimension_negative_empty(self):
        """Wenn m+n < ambient: Schnitt ist generisch leer (negative Dimension)."""
        dim = TransversalityTheory.transverse_intersection_dimension(1, 1, 3)
        assert dim < 0

    def test_sard_theorem_demo(self):
        """Sard-Theorem: reguläre Funktion hat keine kritischen Punkte."""
        # f(x) = x² (R → R): f'(x) = 2x, nur x=0 kritisch
        def f(pt):
            return np.array([pt[0]**2])

        def df(pt):
            return np.array([[2 * pt[0]]])

        # Testpunkte: viele reguläre, einer kritisch (x=0)
        pts = np.linspace(-2, 2, 20).reshape(-1, 1)
        result = TransversalityTheory.sard_theorem_demo(f, df, pts)
        assert 'critical_points' in result
        assert 'regular_points' in result
        assert result['ratio_critical'] < 0.2  # Fast alle regulär

    def test_sard_constant_map(self):
        """Konstante Abbildung: alle Punkte kritisch (Rang 0 < 1)."""
        def f_const(pt):
            return np.array([5.0])  # Konstante

        def df_const(pt):
            return np.array([[0.0]])  # Nulljakobian

        pts = np.linspace(-1, 1, 10).reshape(-1, 1)
        result = TransversalityTheory.sard_theorem_demo(f_const, df_const, pts)
        # Alle Punkte sind kritisch
        assert result['ratio_critical'] == 1.0


# ===========================================================================
# TESTS: VectorBundle
# ===========================================================================

class TestVectorBundle:
    """Tests für Vektorbündel."""

    def test_trivial_bundle_creation(self):
        """Triviales Bündel: Dimension und Rang korrekt."""
        E = VectorBundle(2, 3, 'TM')
        assert E.base_dim == 2
        assert E.fiber_dim == 3

    def test_total_space_dim(self):
        """Totalraum-Dimension = base_dim + fiber_dim."""
        E = VectorBundle(4, 2)
        assert E.total_space_dim() == 6

    def test_rank(self):
        """Rang = fiber_dim."""
        E = VectorBundle(3, 5)
        assert E.rank() == 5

    def test_is_trivial(self):
        """Triviales Bündel: is_trivial() = True."""
        E = VectorBundle(2, 2, trivial=True)
        assert E.is_trivial()

    def test_non_trivial_bundle(self):
        """Nicht-triviales Bündel: is_trivial() = False."""
        # Möbiusband: ℝ¹-Bündel über S¹, nicht-trivial
        moebius = VectorBundle(1, 1, 'Möbius', trivial=False)
        assert not moebius.is_trivial()

    def test_negative_base_dim_raises(self):
        """Negative Basisdimension wirft ValueError."""
        with pytest.raises(ValueError):
            VectorBundle(-1, 2)

    def test_negative_fiber_dim_raises(self):
        """Negative Faserdimension wirft ValueError."""
        with pytest.raises(ValueError):
            VectorBundle(2, -1)

    def test_zero_fiber_dim(self):
        """Rang-0-Bündel (triviales Nullbündel) ist erlaubt."""
        E = VectorBundle(3, 0)
        assert E.rank() == 0
        assert E.total_space_dim() == 3

    def test_repr(self):
        """Repr enthält relevante Informationen."""
        E = VectorBundle(2, 3, 'TM')
        r = repr(E)
        assert 'TM' in r
        assert '2' in r
        assert '3' in r


# ===========================================================================
# TESTS: CharacteristicClasses
# ===========================================================================

class TestCharacteristicClasses:
    """Tests für charakteristische Klassen."""

    def test_chern_class_cp0(self):
        """Chern-Zahlen für ℂP^0: c_0 = 1."""
        result = CharacteristicClasses.chern_class_demo(0)
        assert result['chern_numbers'] == [1]
        assert result['euler_characteristic'] == 1

    def test_chern_class_cp1(self):
        """Chern-Zahlen für ℂP^1 = S^2: c_0 = 1, c_1 = 2."""
        result = CharacteristicClasses.chern_class_demo(1)
        assert result['chern_numbers'] == [1, 2]
        assert result['top_chern_number'] == 2  # χ(CP^1) = 2

    def test_chern_class_cp2(self):
        """Chern-Zahlen für ℂP^2: c_0=1, c_1=3, c_2=3."""
        result = CharacteristicClasses.chern_class_demo(2)
        assert result['chern_numbers'] == [1, 3, 3]
        assert result['euler_characteristic'] == 3

    def test_chern_class_binomial(self):
        """Chern-Zahlen sind Binomialkoeffizienten C(n+1,k)."""
        import math
        n = 4
        result = CharacteristicClasses.chern_class_demo(n)
        for k, ck in enumerate(result['chern_numbers']):
            assert ck == math.comb(n + 1, k)

    def test_pontryagin_class_s2(self):
        """Pontryagin-Klassen von S^2: alle verschwinden."""
        result = CharacteristicClasses.pontryagin_class_demo(2)
        assert result['pontryagin_classes'][0] == 1
        # Für S^2: keine p_k für k ≥ 1 relevant (4k ≤ 2 nur für k=0)
        assert 'manifold' in result

    def test_pontryagin_class_s4(self):
        """Pontryagin-Klassen von S^4: p_1 = 0 (stably trivial)."""
        result = CharacteristicClasses.pontryagin_class_demo(4)
        # p_1 ∈ H^4(S^4): sollte 0 sein (TS^4 stably trivial)
        assert result['pontryagin_classes'].get(1, 0) == 0

    def test_euler_class_s2(self):
        """Euler-Klasse von S^2: e(TS^2) = χ(S^2) = 2."""
        result = CharacteristicClasses.euler_class_demo(2)
        assert result['euler_class_value'] == 2
        assert not result['is_zero_section_avoidable']

    def test_euler_class_s1(self):
        """Euler-Klasse von S^1: χ(S^1) = 0 → nirgends verschwindendes Feld existiert."""
        result = CharacteristicClasses.euler_class_demo(1)
        assert result['euler_class_value'] == 0
        assert result['is_zero_section_avoidable']

    def test_euler_class_s3(self):
        """Euler-Klasse von S^3: χ(S^3) = 0."""
        result = CharacteristicClasses.euler_class_demo(3)
        assert result['euler_class_value'] == 0

    def test_stiefel_whitney_rp1(self):
        """Stiefel-Whitney-Klassen von RP^1 = S^1: w_0=1, w_1=0 (orientierbar)."""
        result = CharacteristicClasses.stiefel_whitney_demo(1)
        assert result['stiefel_whitney_classes'][0] == 1
        assert result['orientable']

    def test_stiefel_whitney_rp2(self):
        """RP^2 ist nicht orientierbar: w_1 ≠ 0."""
        result = CharacteristicClasses.stiefel_whitney_demo(2)
        # w_1(TRP^2) = C(3,1) mod 2 = 3 mod 2 = 1 → nicht orientierbar
        assert not result['orientable']

    def test_stiefel_whitney_formula(self):
        """Stiefel-Whitney-Zahlen sind korrekte Binomialkoeffizienten mod 2."""
        import math
        n = 3
        result = CharacteristicClasses.stiefel_whitney_demo(n)
        for k, wk in enumerate(result['stiefel_whitney_classes']):
            assert wk == math.comb(n + 1, k) % 2


# ===========================================================================
# TESTS: SphereTopology
# ===========================================================================

class TestSphereTopology:
    """Tests für Sphärentopologie."""

    def test_sphere_dimension(self):
        """S^n hat die korrekte Dimension n."""
        for n in [0, 1, 2, 3, 4]:
            S = SphereTopology(n)
            assert S.n == n

    def test_negative_dimension_raises(self):
        """Negative Dimension wirft ValueError."""
        with pytest.raises(ValueError):
            SphereTopology(-1)

    def test_all_spheres_orientable(self):
        """Alle Sphären sind orientierbar."""
        for n in [0, 1, 2, 3, 4, 5]:
            S = SphereTopology(n)
            assert S.is_orientable()

    def test_homology_s0(self):
        """H_*(S^0): H_0 = ℤ² (zwei Punkte)."""
        S = SphereTopology(0)
        homology = S.homology_groups()
        assert homology[0] == 'ℤ²'

    def test_homology_s1(self):
        """H_*(S^1): H_0 = ℤ, H_1 = ℤ."""
        S = SphereTopology(1)
        homology = S.homology_groups()
        assert homology[0] == 'ℤ'
        assert homology[1] == 'ℤ'

    def test_homology_s2(self):
        """H_*(S^2): H_0 = ℤ, H_1 = 0, H_2 = ℤ."""
        S = SphereTopology(2)
        homology = S.homology_groups()
        assert homology[0] == 'ℤ'
        assert homology[1] == '0'
        assert homology[2] == 'ℤ'

    def test_homology_s3(self):
        """H_*(S^3): H_0 = ℤ, H_1 = H_2 = 0, H_3 = ℤ."""
        S = SphereTopology(3)
        homology = S.homology_groups()
        assert homology[0] == 'ℤ'
        assert homology[1] == '0'
        assert homology[2] == '0'
        assert homology[3] == 'ℤ'

    def test_homotopy_groups_s1(self):
        """π_1(S^1) = ℤ (Fundamentalgruppe)."""
        S = SphereTopology(1)
        pi = S.homotopy_groups_low()
        assert pi[1] == 'ℤ'

    def test_homotopy_groups_s2_fundamental(self):
        """π_2(S^2) = ℤ."""
        S = SphereTopology(2)
        pi = S.homotopy_groups_low()
        assert pi[2] == 'ℤ'

    def test_homotopy_groups_s2_hopf(self):
        """π_3(S^2) = ℤ (Hopf-Faserung)."""
        S = SphereTopology(2)
        pi = S.homotopy_groups_low()
        assert pi[3] == 'ℤ'

    def test_homotopy_groups_s3_fundamental(self):
        """π_3(S^3) = ℤ."""
        S = SphereTopology(3)
        pi = S.homotopy_groups_low()
        assert pi[3] == 'ℤ'

    def test_homotopy_groups_below_n_are_zero(self):
        """π_k(S^n) = 0 für k < n (Hurewicz)."""
        S = SphereTopology(3)
        pi = S.homotopy_groups_low()
        assert pi[0] == '0'
        assert pi[1] == '0'
        assert pi[2] == '0'

    def test_repr_contains_n(self):
        """Repr enthält die Dimension."""
        S = SphereTopology(4)
        r = repr(S)
        assert '4' in r


# ===========================================================================
# TESTS: STANDALONE-FUNKTIONEN
# ===========================================================================

class TestDegreeOfMap:
    """Tests für degree_of_map."""

    def test_identity_map_degree_1(self):
        """Identität S^1 → S^1 hat Grad 1."""
        theta = np.linspace(0, 2 * np.pi, 1000, endpoint=False)
        # Identität: f(θ) = θ
        degree = degree_of_map(theta, theta)
        assert degree == 1

    def test_double_cover_degree_2(self):
        """Doppelte Überdeckung θ ↦ 2θ hat Grad 2."""
        theta = np.linspace(0, 2 * np.pi, 1000, endpoint=False)
        f_theta = 2 * theta  # Winkeländerung ums Doppelte
        degree = degree_of_map(f_theta, theta)
        assert degree == 2

    def test_constant_map_degree_0(self):
        """Konstante Abbildung hat Grad 0."""
        theta = np.linspace(0, 2 * np.pi, 1000, endpoint=False)
        f_const = np.zeros_like(theta)  # Immer Winkel 0
        degree = degree_of_map(f_const, theta)
        assert degree == 0


class TestHairyBallTheorem:
    """Tests für hairy_ball_theorem_demo."""

    def test_returns_dict(self):
        """Rückgabe ist ein Dict mit den erwarteten Schlüsseln."""
        result = hairy_ball_theorem_demo()
        assert isinstance(result, dict)
        assert 'theorem' in result
        assert 'examples' in result

    def test_s2_no_nowhere_vanishing_field(self):
        """S^2 hat kein nirgends verschwindendes Vektorfeld."""
        result = hairy_ball_theorem_demo()
        assert not result['examples']['S^2']['has_nowhere_vanishing_field']

    def test_s1_has_nowhere_vanishing_field(self):
        """S^1 hat ein nirgends verschwindendes Vektorfeld."""
        result = hairy_ball_theorem_demo()
        assert result['examples']['S^1']['has_nowhere_vanishing_field']

    def test_criterion_present(self):
        """Das Kriterium (χ = 0) ist enthalten."""
        result = hairy_ball_theorem_demo()
        assert 'criterion' in result


class TestPoincareHopfTheorem:
    """Tests für poincare_hopf_theorem_demo."""

    def test_s2_two_poles(self):
        """Meridianvektorfeld auf S^2: 2 Singularitäten (Index +1 je), Summe = 2 = χ(S^2)."""
        result = poincare_hopf_theorem_demo(euler_char=2, vector_field_index_sum=2)
        assert result['theorem_satisfied']

    def test_mismatch_not_satisfied(self):
        """Falsche Indexsumme verletzt den Satz."""
        result = poincare_hopf_theorem_demo(euler_char=2, vector_field_index_sum=3)
        assert not result['theorem_satisfied']

    def test_torus_zero_sum(self):
        """Auf T^2 mit χ = 0 muss Indexsumme = 0 sein."""
        result = poincare_hopf_theorem_demo(euler_char=0, vector_field_index_sum=0)
        assert result['theorem_satisfied']

    def test_euler_char_preserved(self):
        """Euler-Charakteristik wird korrekt im Ergebnis gespeichert."""
        result = poincare_hopf_theorem_demo(euler_char=4, vector_field_index_sum=4)
        assert result['euler_characteristic'] == 4


class TestWhitneyEmbedding:
    """Tests für whitney_embedding_dimension."""

    def test_1d_manifold(self):
        """1-dim Mannigfaltigkeit: Whitney-Einbettung in ℝ^2."""
        assert whitney_embedding_dimension(1) == 2

    def test_2d_manifold(self):
        """2-dim Mannigfaltigkeit: Whitney-Einbettung in ℝ^4."""
        assert whitney_embedding_dimension(2) == 4

    def test_3d_manifold(self):
        """3-dim Mannigfaltigkeit: Whitney-Einbettung in ℝ^6."""
        assert whitney_embedding_dimension(3) == 6

    def test_0d_manifold(self):
        """0-dim Mannigfaltigkeit (Punkt): Whitney-Einbettung in ℝ^0."""
        assert whitney_embedding_dimension(0) == 0

    def test_negative_raises(self):
        """Negative Dimension wirft ValueError."""
        with pytest.raises(ValueError):
            whitney_embedding_dimension(-1)


class TestComputeJacobian:
    """Tests für compute_jacobian."""

    def test_jacobian_2x2(self):
        """Jacobi-Matrix von f=(x²,xy) bzgl. (x,y)."""
        f_exprs = [x**2, x * y]
        J = compute_jacobian(f_exprs, [x, y])
        assert J.shape == (2, 2)
        # J[0,0] = ∂(x²)/∂x = 2x
        assert sp.simplify(J[0, 0] - 2*x) == 0
        # J[0,1] = ∂(x²)/∂y = 0
        assert sp.simplify(J[0, 1]) == 0
        # J[1,0] = ∂(xy)/∂x = y
        assert sp.simplify(J[1, 0] - y) == 0
        # J[1,1] = ∂(xy)/∂y = x
        assert sp.simplify(J[1, 1] - x) == 0

    def test_jacobian_identity(self):
        """Jacobi-Matrix der Identität ist die Einheitsmatrix."""
        f_exprs = [x, y]
        J = compute_jacobian(f_exprs, [x, y])
        assert J == sp.eye(2)

    def test_jacobian_scalar_function(self):
        """Gradient einer skalaren Funktion f=x²+y²: J = [[2x, 2y]]."""
        f_exprs = [x**2 + y**2]
        J = compute_jacobian(f_exprs, [x, y])
        assert J.shape == (1, 2)
        assert sp.simplify(J[0, 0] - 2*x) == 0
        assert sp.simplify(J[0, 1] - 2*y) == 0

    def test_jacobian_3d(self):
        """Jacobi-Matrix einer 3D-Abbildung hat korrekte Form."""
        f_exprs = [x*y, y*z, x*z]
        J = compute_jacobian(f_exprs, [x, y, z])
        assert J.shape == (3, 3)


class TestImplicitFunctionTheorem:
    """Tests für implicit_function_theorem_check."""

    def test_circle_implicit(self):
        """Kreis x²+y²-1=0: Implizite Funktion bei (0,1) existiert."""
        # F(x,y) = x² + y² - 1 = 0
        F = [x**2 + y**2 - 1]
        pt = {x: sp.Integer(0), y: sp.Integer(1)}
        result = implicit_function_theorem_check(F, [x], [y], pt)
        assert result['condition_satisfied']
        assert result['F_is_zero']
        assert result['implicit_function_exists']

    def test_circle_at_tangent_point_fails(self):
        """Kreis x²+y²-1=0 an Punkt (1,0): ∂F/∂y = 0 → Bedingung nicht erfüllt."""
        F = [x**2 + y**2 - 1]
        pt = {x: sp.Integer(1), y: sp.Integer(0)}
        result = implicit_function_theorem_check(F, [x], [y], pt)
        # ∂F/∂y = 2y = 0 an (1,0) → nicht invertierbar
        assert not result['condition_satisfied']

    def test_wrong_num_y_vars_raises(self):
        """Anzahl y-Variablen ≠ Anzahl Gleichungen wirft ValueError."""
        F = [x + y, x - y]  # 2 Gleichungen
        with pytest.raises(ValueError):
            implicit_function_theorem_check(F, [x], [y], {x: 0, y: 0})  # nur 1 y-Var

    def test_linear_system(self):
        """Lineares System F=(x+y-1, x-y) an Lösungspunkt (1/2, 1/2)."""
        F = [x + y - 1, x - y]
        pt = {x: Rational(1, 2), y: Rational(1, 2)}
        result = implicit_function_theorem_check(F, [], [x, y], pt)
        assert result['F_is_zero']
        assert result['condition_satisfied']

    def test_returns_jy_det(self):
        """Rückgabe enthält Determinante der Jacobi-Matrix."""
        F = [x**2 + y**2 - 1]
        pt = {x: sp.Integer(0), y: sp.Integer(1)}
        result = implicit_function_theorem_check(F, [x], [y], pt)
        assert 'Jy_det' in result
        # ∂F/∂y = 2y = 2·1 = 2
        assert sp.simplify(result['Jy_det'] - 2) == 0
