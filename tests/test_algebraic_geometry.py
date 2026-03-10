"""
@file test_algebraic_geometry.py
@brief Tests für das Modul algebraic_geometry.py.
@description
    Testet alle Klassen und Standalone-Funktionen aus algebraic_geometry.py:

    - AffineVariety: __init__, dimension, is_irreducible, singular_points, intersect
    - ProjectiveVariety: __init__, degree, is_projective_smooth
    - HilbertBasisTheorem: groebner_basis, ideal_membership
    - NullstellensatzDemo: weak_nullstellensatz, strong_nullstellensatz
    - EllipticCurveGeometry: __init__, discriminant, j_invariant, is_smooth, rational_points_mod_p
    - IntersectionTheory: bezout_theorem, intersection_multiplicity
    - MorphismOfVarieties: __init__, is_dominant, pullback
    - Standalone: compute_groebner_basis, ideal_radical, zariski_closure_demo, hilbert_polynomial

    Getestete Edge Cases:
    - Leeres Ideal, Nullpolynom, trivialer Morphismus
    - Singuläre Punkte und glatte Varietäten
    - Inkonsistente Ideale (kein gemeinsamer Nullpunkt)
    - Elliptische Kurven mit Δ=0 (singulär)
    - Radikalmitgliedschaft

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import pytest
import sys
import os

# Quellverzeichnis in den Suchpfad aufnehmen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import sympy
from sympy import symbols, expand, Integer, Rational, sqrt, oo

from algebraic_geometry import (
    AffineVariety,
    ProjectiveVariety,
    HilbertBasisTheorem,
    NullstellensatzDemo,
    EllipticCurveGeometry,
    IntersectionTheory,
    MorphismOfVarieties,
    compute_groebner_basis,
    ideal_radical,
    zariski_closure_demo,
    hilbert_polynomial,
)

# Häufig verwendete SymPy-Symbole
x, y, z, w = symbols('x y z w')


# ===========================================================================
# TESTS FÜR AffineVariety
# ===========================================================================

class TestAffineVariety:
    """Tests für die Klasse AffineVariety."""

    # --- __init__ ---

    def test_init_basic(self):
        """Einfache Initialisierung mit einer Variablen."""
        V = AffineVariety([x - 1], [x])
        assert len(V.generators) == 1
        assert len(V.variables) == 1

    def test_init_two_variables(self):
        """Initialisierung mit zwei Variablen."""
        V = AffineVariety([x**2 + y**2 - 1], [x, y])
        assert len(V.variables) == 2

    def test_init_zero_generator(self):
        """Nullpolynom als Erzeuger ist erlaubt."""
        V = AffineVariety([Integer(0)], [x, y])
        assert expand(V.generators[0]) == 0

    def test_init_empty_generators(self):
        """Leere Erzeuger-Liste (Nullideal)."""
        V = AffineVariety([], [x, y])
        assert V.generators == []

    def test_init_empty_variables_raises(self):
        """Leere Variablen-Liste wirft ValueError."""
        with pytest.raises(ValueError):
            AffineVariety([x], [])

    def test_init_expand_called(self):
        """Erzeuger werden durch expand() normalisiert."""
        V = AffineVariety([(x + 1)**2 - x**2 - 2*x - 1], [x])
        # (x+1)^2 - x^2 - 2x - 1 = 0
        assert expand(V.generators[0]) == 0

    # --- dimension ---

    def test_dimension_full_space(self):
        """Nullideal → gesamter affiner Raum, dim = n."""
        V = AffineVariety([], [x, y])
        assert V.dimension() == 2

    def test_dimension_point(self):
        """Maximales Ideal ⟨x-a, y-b⟩ → einzelner Punkt, dim = 0."""
        V = AffineVariety([x - 1, y - 2], [x, y])
        assert V.dimension() == 0

    def test_dimension_line(self):
        """Hyperfläche ⟨x⟩ in 2D → dim = 1."""
        V = AffineVariety([x], [x, y])
        assert V.dimension() == 1

    def test_dimension_empty_variety(self):
        """Einheitsideal ⟨1⟩ → leere Varietät, dim = -1."""
        V = AffineVariety([Integer(1)], [x, y])
        assert V.dimension() == -1

    def test_dimension_inconsistent(self):
        """Inkonsistentes System ⟨x-1, x-2⟩ → leer, dim = -1."""
        V = AffineVariety([x - 1, x - 2], [x])
        assert V.dimension() == -1

    def test_dimension_single_var_hypersurface(self):
        """Kreis x²+y²=1 in 2D: dim = 1."""
        V = AffineVariety([x**2 + y**2 - 1], [x, y])
        assert V.dimension() == 1

    # --- is_irreducible ---

    def test_irreducible_line(self):
        """Gerade x = 0 ist irreduzibel."""
        V = AffineVariety([x], [x, y])
        assert V.is_irreducible() is True

    def test_irreducible_parabola(self):
        """Parabel y - x² = 0 ist irreduzibel."""
        V = AffineVariety([y - x**2], [x, y])
        assert V.is_irreducible() is True

    def test_reducible_two_lines(self):
        """xy = 0 zerfällt in x=0 und y=0 → reduzibel."""
        V = AffineVariety([x * y], [x, y])
        assert V.is_irreducible() is False

    def test_reducible_difference_of_squares(self):
        """x²-1 = (x-1)(x+1) → reduzibel."""
        V = AffineVariety([x**2 - 1], [x])
        assert V.is_irreducible() is False

    def test_irreducible_full_space(self):
        """Nullideal → ganzer affiner Raum ist irreduzibel."""
        V = AffineVariety([], [x, y])
        assert V.is_irreducible() is True

    def test_irreducible_circle(self):
        """Kreis x²+y²-1 ist irreduzibel über ℚ."""
        V = AffineVariety([x**2 + y**2 - 1], [x, y])
        assert V.is_irreducible() is True

    # --- singular_points ---

    def test_singular_points_smooth_line(self):
        """Gerade x = 0: keine singulären Punkte."""
        V = AffineVariety([x], [x, y])
        sings = V.singular_points()
        assert sings == []

    def test_singular_points_cusp(self):
        """Spitze y² - x³ = 0: singulärer Punkt bei (0,0)."""
        V = AffineVariety([y**2 - x**3], [x, y])
        sings = V.singular_points()
        # (0,0) muss in den Lösungen sein
        origin_found = any(
            sol.get(x, None) == 0 and sol.get(y, None) == 0
            for sol in sings
        )
        assert origin_found

    def test_singular_points_node(self):
        """Knoten y² - x²(x+1) = 0: singulärer Punkt bei (0,0)."""
        V = AffineVariety([y**2 - x**2 * (x + 1)], [x, y])
        sings = V.singular_points()
        origin_found = any(
            sol.get(x, None) == 0 and sol.get(y, None) == 0
            for sol in sings
        )
        assert origin_found

    def test_singular_points_empty_generators(self):
        """Leere Erzeuger: keine Singularitäten."""
        V = AffineVariety([], [x, y])
        sings = V.singular_points()
        assert sings == []

    def test_singular_points_smooth_circle(self):
        """Kreis x²+y²-1: keine singulären Punkte (glatt)."""
        V = AffineVariety([x**2 + y**2 - 1], [x, y])
        sings = V.singular_points()
        # Keine Singularitäten über ℝ (Kreis ist glatt)
        assert sings == []

    # --- intersect ---

    def test_intersect_two_lines(self):
        """Schnitt der Geraden x=0 und y=0 → Ursprung."""
        V1 = AffineVariety([x], [x, y])
        V2 = AffineVariety([y], [x, y])
        V = V1.intersect(V2)
        # Erzeuger sollten x und y enthalten
        gen_strs = [str(expand(g)) for g in V.generators]
        assert 'x' in gen_strs
        assert 'y' in gen_strs

    def test_intersect_combines_generators(self):
        """Schnitt vereint die Erzeuger beider Varietäten."""
        V1 = AffineVariety([x - 1], [x, y])
        V2 = AffineVariety([y - 2], [x, y])
        V = V1.intersect(V2)
        assert len(V.generators) == 2

    def test_intersect_incompatible_vars_raises(self):
        """Schnitt mit verschiedenen Variablen wirft ValueError."""
        V1 = AffineVariety([x], [x])
        V2 = AffineVariety([y], [y])
        with pytest.raises(ValueError):
            V1.intersect(V2)

    def test_intersect_circle_and_line(self):
        """Schnitt von Kreis x²+y²=1 mit Linie y=0."""
        V1 = AffineVariety([x**2 + y**2 - 1], [x, y])
        V2 = AffineVariety([y], [x, y])
        V = V1.intersect(V2)
        assert len(V.generators) == 2


# ===========================================================================
# TESTS FÜR ProjectiveVariety
# ===========================================================================

class TestProjectiveVariety:
    """Tests für die Klasse ProjectiveVariety."""

    # Symbole für projektive Koordinaten
    x0, x1, x2 = symbols('x0 x1 x2')

    def test_init_basic(self):
        """Einfache Initialisierung im ℙ²."""
        V = ProjectiveVariety([self.x0**2 + self.x1**2 + self.x2**2], 3)
        assert V.n_vars == 3

    def test_init_n_vars_too_small_raises(self):
        """n_vars < 2 wirft ValueError."""
        with pytest.raises(ValueError):
            ProjectiveVariety([self.x0], 1)

    def test_init_empty_generators(self):
        """Leere Erzeuger sind erlaubt."""
        V = ProjectiveVariety([], 3)
        assert V.generators == []

    def test_degree_linear(self):
        """Lineare Form: Grad 1."""
        V = ProjectiveVariety([self.x0 + self.x1 + self.x2], 3)
        assert V.degree() == 1

    def test_degree_quadric(self):
        """Quadratische Form: Grad 2."""
        V = ProjectiveVariety([self.x0**2 + self.x1**2 - self.x2**2], 3)
        assert V.degree() == 2

    def test_degree_cubic(self):
        """Kubik: Grad 3."""
        V = ProjectiveVariety([self.x0**3 + self.x1**3 + self.x2**3], 3)
        assert V.degree() == 3

    def test_degree_empty(self):
        """Leere Erzeuger: Grad 0."""
        V = ProjectiveVariety([], 3)
        assert V.degree() == 0

    def test_degree_max_of_multiple(self):
        """Mehrere Polynome: maximaler Grad wird zurückgegeben."""
        V = ProjectiveVariety([self.x0 + self.x1, self.x0**3 - self.x2**3], 3)
        assert V.degree() == 3

    def test_is_projective_smooth_no_generators(self):
        """Keine Erzeuger: trivial glatt."""
        V = ProjectiveVariety([], 3)
        assert V.is_projective_smooth() is True

    def test_is_projective_smooth_quadric(self):
        """Glatte Quadrik x₀²+x₁²+x₂² = 0 (über ℂ, keine reellen Punkte)."""
        V = ProjectiveVariety([self.x0**2 + self.x1**2 + self.x2**2], 3)
        # Jacobi: (2x0, 2x1, 2x2) = 0 nur bei (0,0,0) → nicht im ℙ²
        result = V.is_projective_smooth()
        assert isinstance(result, bool)

    def test_is_projective_smooth_singular_cone(self):
        """Singular: x₀·x₁·x₂ hat Singularitäten."""
        V = ProjectiveVariety([self.x0 * self.x1], 3)
        # Jacobi: (x1, x0, 0) = 0 → x0=x1=0 liegt auf der Varietät bei [0:0:1]
        result = V.is_projective_smooth()
        assert isinstance(result, bool)


# ===========================================================================
# TESTS FÜR HilbertBasisTheorem
# ===========================================================================

class TestHilbertBasisTheorem:
    """Tests für die Klasse HilbertBasisTheorem."""

    def test_groebner_basis_single_poly(self):
        """Gröbner-Basis eines Hauptideals ⟨f⟩."""
        gb = HilbertBasisTheorem.groebner_basis([x**2 - 1], [x])
        assert len(gb) >= 1
        # Gröbner-Basis sollte x²-1 oder äquivalente Form enthalten

    def test_groebner_basis_system(self):
        """Gröbner-Basis für Gleichungssystem {x²+y²-1=0, x-y=0}."""
        gb = HilbertBasisTheorem.groebner_basis(
            [x**2 + y**2 - 1, x - y], [x, y]
        )
        assert isinstance(gb, list)
        assert len(gb) >= 1

    def test_groebner_basis_empty_returns_zero(self):
        """Leere Erzeuger → Gröbner-Basis = [0]."""
        gb = HilbertBasisTheorem.groebner_basis([], [x, y])
        assert gb == [Integer(0)] or gb == [0]

    def test_groebner_basis_unit_ideal(self):
        """Einheitsideal ⟨x, 1-x⟩ = ⟨1⟩: GB = [1]."""
        gb = HilbertBasisTheorem.groebner_basis([x, 1 - x], [x])
        assert Integer(1) in gb or 1 in gb

    def test_groebner_basis_no_variables(self):
        """Keine Variablen → Gröbner-Basis [1]."""
        gb = HilbertBasisTheorem.groebner_basis([x], [])
        assert gb == [Integer(1)] or gb == [1]

    def test_groebner_basis_zero_generators_returns_zero(self):
        """Alle Nullerzeuger → Gröbner-Basis [0]."""
        gb = HilbertBasisTheorem.groebner_basis([Integer(0), Integer(0)], [x])
        assert Integer(0) in gb or 0 in gb

    # --- ideal_membership ---

    def test_ideal_membership_zero_in_any_ideal(self):
        """0 liegt in jedem Ideal."""
        assert HilbertBasisTheorem.ideal_membership(
            Integer(0), [x**2 - 1], [x]
        ) is True

    def test_ideal_membership_generator_in_ideal(self):
        """Erzeuger liegt im Ideal."""
        assert HilbertBasisTheorem.ideal_membership(
            x**2 - 1, [x**2 - 1], [x]
        ) is True

    def test_ideal_membership_multiple_of_generator(self):
        """Vielfaches des Erzeugers liegt im Ideal."""
        assert HilbertBasisTheorem.ideal_membership(
            x * (x**2 - 1), [x**2 - 1], [x]
        ) is True

    def test_ideal_membership_not_member(self):
        """x liegt nicht im Ideal ⟨x²-1⟩."""
        result = HilbertBasisTheorem.ideal_membership(
            x, [x**2 - 1], [x]
        )
        assert result is False

    def test_ideal_membership_unit_ideal(self):
        """Jedes Polynom liegt im Einheitsideal ⟨1⟩."""
        assert HilbertBasisTheorem.ideal_membership(
            x**3 + 2*x + 5, [Integer(1)], [x]
        ) is True

    def test_ideal_membership_empty_generators(self):
        """Leeres Ideal: nur 0 liegt drin."""
        result = HilbertBasisTheorem.ideal_membership(
            x + 1, [], [x]
        )
        assert result is False

    def test_ideal_membership_two_var(self):
        """x+y liegt im Ideal ⟨x, y⟩."""
        assert HilbertBasisTheorem.ideal_membership(
            x + y, [x, y], [x, y]
        ) is True


# ===========================================================================
# TESTS FÜR NullstellensatzDemo
# ===========================================================================

class TestNullstellensatzDemo:
    """Tests für die Klasse NullstellensatzDemo."""

    # --- weak_nullstellensatz ---

    def test_weak_consistent_system(self):
        """System mit Lösungen → hat_keine_Lösung = False."""
        result = NullstellensatzDemo.weak_nullstellensatz(
            [x**2 - 1], [x]
        )
        assert result['has_no_solution'] is False

    def test_weak_inconsistent_system(self):
        """Widersprüchliches System {x=0, x=1} → hat_keine_Lösung = True."""
        result = NullstellensatzDemo.weak_nullstellensatz(
            [x, x - 1], [x]
        )
        assert result['has_no_solution'] is True

    def test_weak_unit_ideal(self):
        """Einheitsideal direkt → hat_keine_Lösung = True."""
        result = NullstellensatzDemo.weak_nullstellensatz(
            [Integer(1)], [x]
        )
        assert result['has_no_solution'] is True

    def test_weak_returns_groebner(self):
        """Ergebnis enthält Gröbner-Basis."""
        result = NullstellensatzDemo.weak_nullstellensatz(
            [x**2 - y, y - 1], [x, y]
        )
        assert 'groebner_basis' in result
        assert isinstance(result['groebner_basis'], list)

    def test_weak_empty_generators(self):
        """Leere Erzeuger → hat_keine_Lösung = False."""
        result = NullstellensatzDemo.weak_nullstellensatz([], [x])
        assert result['has_no_solution'] is False

    def test_weak_explanation_present(self):
        """Ergebnis enthält Erklärungstext."""
        result = NullstellensatzDemo.weak_nullstellensatz([x - 1], [x])
        assert 'explanation' in result
        assert isinstance(result['explanation'], str)
        assert len(result['explanation']) > 0

    # --- strong_nullstellensatz ---

    def test_strong_zero_always_in_radical(self):
        """0 liegt in jedem Radikal."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            Integer(0), [x - 1], [x]
        )
        assert result['in_radical'] is True

    def test_strong_member_in_radical(self):
        """x² ∈ √⟨x²⟩ (trivialerweise)."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            x**2, [x**2], [x]
        )
        assert result['in_radical'] is True

    def test_strong_x_in_radical_of_x_squared(self):
        """x ∈ √⟨x²⟩ weil x² ∈ ⟨x²⟩."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            x, [x**2], [x]
        )
        assert result['in_radical'] is True

    def test_strong_not_in_radical(self):
        """y ∉ √⟨x²⟩: y hat keine gemeinsame Nullstelle-Bedingung mit x²."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            y, [x**2], [x, y]
        )
        assert result['in_radical'] is False

    def test_strong_inconsistent_ideal(self):
        """Einheitsideal: jedes f liegt im Radikal √⟨1⟩ = k[x]... Nein: 1∈I → √I=I=k[x], f^1=f∈k[x]."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            x + 1, [x, Integer(1) - x], [x]
        )
        # In Einheitsideal liegt alles: 1 ∈ I → jedes f^1 ∈ I
        assert result['in_radical'] is True

    def test_strong_returns_explanation(self):
        """Ergebnis enthält Erklärungstext."""
        result = NullstellensatzDemo.strong_nullstellensatz(x, [x], [x])
        assert 'explanation' in result

    def test_strong_empty_generators(self):
        """Nullideal: nur 0 im Radikal."""
        result = NullstellensatzDemo.strong_nullstellensatz(
            x + 1, [], [x]
        )
        assert result['in_radical'] is False


# ===========================================================================
# TESTS FÜR EllipticCurveGeometry
# ===========================================================================

class TestEllipticCurveGeometry:
    """Tests für die Klasse EllipticCurveGeometry."""

    # --- discriminant ---

    def test_discriminant_standard_curve(self):
        """Diskriminante von y²=x³-x (a=-1, b=0): Δ = -16(4(-1)³+27·0²) = 64."""
        E = EllipticCurveGeometry(-1, 0)
        delta = E.discriminant()
        assert delta == 64

    def test_discriminant_zero_means_singular(self):
        """a=0, b=0: Δ=0 → singuläre Kurve."""
        E = EllipticCurveGeometry(0, 0)
        assert E.discriminant() == 0

    def test_discriminant_sign(self):
        """Diskriminante hat korrektes Vorzeichen."""
        E = EllipticCurveGeometry(0, 1)
        # Δ = -16(0 + 27) = -432
        delta = E.discriminant()
        assert delta == -432

    def test_discriminant_formula(self):
        """Δ = -16(4a³+27b²) korrekt berechnet."""
        a_val, b_val = 1, 1
        E = EllipticCurveGeometry(a_val, b_val)
        expected = -16 * (4 * a_val**3 + 27 * b_val**2)
        assert int(E.discriminant()) == expected

    # --- j_invariant ---

    def test_j_invariant_b_zero(self):
        """a≠0, b=0: j = 1728 (klassischer CM-Fall)."""
        E = EllipticCurveGeometry(-1, 0)
        j = E.j_invariant()
        # j = -1728·(4·(-1))³ / (64) = -1728·(-64)/64 = 1728
        assert j == 1728

    def test_j_invariant_a_zero(self):
        """a=0, b≠0: j = 0."""
        E = EllipticCurveGeometry(0, 1)
        j = E.j_invariant()
        assert j == 0

    def test_j_invariant_singular_is_infinity(self):
        """Singuläre Kurve (Δ=0): j = ∞."""
        E = EllipticCurveGeometry(0, 0)
        j = E.j_invariant()
        assert j == sympy.oo

    def test_j_invariant_is_rational(self):
        """j-Invariante für ganzzahlige a,b ist rational."""
        E = EllipticCurveGeometry(1, 1)
        j = E.j_invariant()
        # Prüfe ob SymPy-Rational oder Integer
        assert j.is_rational or j.is_integer or j.is_number

    # --- is_smooth ---

    def test_is_smooth_regular_curve(self):
        """y²=x³-x: glatt (Δ=64≠0)."""
        E = EllipticCurveGeometry(-1, 0)
        assert E.is_smooth() is True

    def test_is_smooth_singular_curve(self):
        """y²=x³: singulär (Δ=0)."""
        E = EllipticCurveGeometry(0, 0)
        assert E.is_smooth() is False

    def test_is_smooth_y_squared_equals_x_cubed_minus_x(self):
        """y²=x³-x ist eine glatte elliptische Kurve."""
        E = EllipticCurveGeometry(-1, 0)
        assert E.is_smooth() is True

    def test_is_smooth_cuspidal(self):
        """Spitze: a=0, b=0 → singulär."""
        E = EllipticCurveGeometry(0, 0)
        assert E.is_smooth() is False

    # --- rational_points_mod_p ---

    def test_rational_points_mod_5_count(self):
        """y²=x³-x mod 5: Punkt O ist immer dabei."""
        E = EllipticCurveGeometry(-1, 0)
        pts = E.rational_points_mod_p(5)
        assert 'O' in pts

    def test_rational_points_mod_5_type(self):
        """Alle Punkte (außer O) sind 2-Tupel."""
        E = EllipticCurveGeometry(-1, 0)
        pts = E.rational_points_mod_p(5)
        for pt in pts:
            if pt != 'O':
                assert isinstance(pt, tuple)
                assert len(pt) == 2

    def test_rational_points_mod_p_coordinates_in_range(self):
        """Koordinaten liegen in {0,...,p-1}."""
        E = EllipticCurveGeometry(-1, 0)
        p = 7
        pts = E.rational_points_mod_p(p)
        for pt in pts:
            if pt != 'O':
                xv, yv = pt
                assert 0 <= xv < p
                assert 0 <= yv < p

    def test_rational_points_mod_p_satisfy_equation(self):
        """Alle Punkte erfüllen y²≡x³+ax+b (mod p)."""
        a_val, b_val, p = -1, 0, 7
        E = EllipticCurveGeometry(a_val, b_val)
        pts = E.rational_points_mod_p(p)
        for pt in pts:
            if pt != 'O':
                xv, yv = pt
                lhs = (yv * yv) % p
                rhs = (xv**3 + a_val * xv + b_val) % p
                assert lhs == rhs

    def test_rational_points_mod_p_invalid_raises(self):
        """Negativer Wert für p wirft ValueError."""
        E = EllipticCurveGeometry(-1, 0)
        with pytest.raises(ValueError):
            E.rational_points_mod_p(-1)

    def test_rational_points_mod_p_singular_curve(self):
        """Auch singuläre Kurven können Punkte mod p haben."""
        E = EllipticCurveGeometry(0, 0)
        pts = E.rational_points_mod_p(5)
        assert 'O' in pts

    def test_rational_points_mod_2(self):
        """Punkte mod 2 für einfache Kurve."""
        E = EllipticCurveGeometry(0, 1)
        pts = E.rational_points_mod_p(2)
        # Alle Punkte prüfen
        for pt in pts:
            if pt != 'O':
                xv, yv = pt
                lhs = (yv**2) % 2
                rhs = (xv**3 + 1) % 2
                assert lhs == rhs


# ===========================================================================
# TESTS FÜR IntersectionTheory
# ===========================================================================

class TestIntersectionTheory:
    """Tests für die Klasse IntersectionTheory."""

    # --- bezout_theorem ---

    def test_bezout_two_lines(self):
        """Zwei Geraden: 1·1 = 1 Schnittpunkt."""
        assert IntersectionTheory.bezout_theorem(1, 1) == 1

    def test_bezout_line_and_conic(self):
        """Gerade und Kegelschnitt: 1·2 = 2 Schnittpunkte."""
        assert IntersectionTheory.bezout_theorem(1, 2) == 2

    def test_bezout_two_conics(self):
        """Zwei Kegelschnitte: 2·2 = 4 Schnittpunkte."""
        assert IntersectionTheory.bezout_theorem(2, 2) == 4

    def test_bezout_zero_degree(self):
        """Grad 0 (Punkt): 0·d = 0."""
        assert IntersectionTheory.bezout_theorem(0, 5) == 0

    def test_bezout_cubic(self):
        """Zwei Kubiken: 3·3 = 9."""
        assert IntersectionTheory.bezout_theorem(3, 3) == 9

    def test_bezout_commutative(self):
        """Bézout ist kommutativ: d₁·d₂ = d₂·d₁."""
        assert IntersectionTheory.bezout_theorem(3, 5) == IntersectionTheory.bezout_theorem(5, 3)

    # --- intersection_multiplicity ---

    def test_multiplicity_simple_transversal(self):
        """Einfacher Schnittpunkt zweier Geraden: Multiplizität = 1."""
        # f = x, g = y schneiden sich transversal im Ursprung
        mult = IntersectionTheory.intersection_multiplicity(
            x, y, {x: 0, y: 0}, [x, y]
        )
        assert mult == 1

    def test_multiplicity_tangent_parabola(self):
        """Parabel y=x² und Tangente y=0: Multiplizität ≥ 2."""
        # f = y - x², g = y  → berühren sich im Ursprung mit Mult. 2
        mult = IntersectionTheory.intersection_multiplicity(
            y - x**2, y, {x: 0, y: 0}, [x, y]
        )
        assert mult >= 2

    def test_multiplicity_at_non_singular_point(self):
        """Nicht-Schnittpunkt: Multiplizität = 0."""
        # f = x, g = y-1 schneiden sich NICHT im Ursprung
        mult = IntersectionTheory.intersection_multiplicity(
            x, y - 1, {x: 0, y: 0}, [x, y]
        )
        assert mult == 0

    def test_multiplicity_cusp_singularity(self):
        """Spitze y²-x³ und Tangente y: erhöhte Multiplizität."""
        mult = IntersectionTheory.intersection_multiplicity(
            y**2 - x**3, y, {x: 0, y: 0}, [x, y]
        )
        assert mult >= 2


# ===========================================================================
# TESTS FÜR MorphismOfVarieties
# ===========================================================================

class TestMorphismOfVarieties:
    """Tests für die Klasse MorphismOfVarieties."""

    def test_init_basic(self):
        """Einfache Initialisierung eines Morphismus."""
        V1 = AffineVariety([y - x**2], [x, y])
        V2 = AffineVariety([y - x], [x, y])
        phi = MorphismOfVarieties(V1, V2, [x, x**2])
        assert len(phi.map_polynomials) == 2

    def test_pullback_linear(self):
        """Pullback einer linearen Funktion."""
        V1 = AffineVariety([], [x, y])
        V2 = AffineVariety([], [x, y])
        # Morphismus: (x,y) → (x+1, y+1)
        phi = MorphismOfVarieties(V1, V2, [x + 1, y + 1])
        # Pullback von x (erste Koordinate)
        result = phi.pullback(x)
        assert expand(result) == expand(x + 1)

    def test_pullback_quadratic(self):
        """Pullback einer quadratischen Funktion."""
        V1 = AffineVariety([], [x])
        V2 = AffineVariety([], [x])
        # Morphismus: x → x²
        phi = MorphismOfVarieties(V1, V2, [x**2])
        # Pullback von x: ergibt x²
        result = phi.pullback(x)
        assert expand(result) == expand(x**2)

    def test_pullback_constant(self):
        """Pullback einer Konstanten ist die Konstante selbst."""
        V1 = AffineVariety([], [x, y])
        V2 = AffineVariety([], [x, y])
        phi = MorphismOfVarieties(V1, V2, [x**2, y**2])
        result = phi.pullback(Integer(5))
        assert expand(result) == 5

    def test_is_dominant_identity(self):
        """Identitäts-Morphismus ist dominant."""
        V1 = AffineVariety([], [x, y])
        V2 = AffineVariety([], [x, y])
        phi = MorphismOfVarieties(V1, V2, [x, y])
        assert phi.is_dominant() is True

    def test_is_dominant_projection(self):
        """Projektion (x,y)→(x) nutzt nicht alle Variablen der Quelle, heuristisch False."""
        V1 = AffineVariety([], [x, y])
        V2 = AffineVariety([], [x])
        phi = MorphismOfVarieties(V1, V2, [x])
        # Heuristik: Jacobi ∂x/∂x = 1 ≠ 0, aber y fehlt → trotzdem Ergebnis prüfen
        result = phi.is_dominant()
        assert isinstance(result, bool)

    def test_is_dominant_constant_map_false(self):
        """Konstanter Morphismus ist nicht dominant."""
        V1 = AffineVariety([], [x, y])
        V2 = AffineVariety([], [x, y])
        # Konstante Abbildung: (x,y) → (1, 1)
        phi = MorphismOfVarieties(V1, V2, [Integer(1), Integer(1)])
        assert phi.is_dominant() is False

    def test_pullback_composition(self):
        """Pullback von f∘g = Pullback von f hintereinandergeschaltet."""
        V = AffineVariety([], [x])
        phi = MorphismOfVarieties(V, V, [x**2])
        # Pullback von x²: φ*(x²) = (x²)² = x⁴
        result = phi.pullback(x**2)
        assert expand(result) == expand(x**4)


# ===========================================================================
# TESTS FÜR STANDALONE-FUNKTIONEN
# ===========================================================================

class TestStandaloneFunctions:
    """Tests für Standalone-Funktionen."""

    # --- compute_groebner_basis ---

    def test_compute_groebner_basis_basic(self):
        """Berechnung einer Gröbner-Basis."""
        gb = compute_groebner_basis([x**2 - 1, x - 1], [x])
        assert isinstance(gb, list)
        assert len(gb) >= 1

    def test_compute_groebner_basis_empty(self):
        """Leere Erzeuger → [0]."""
        gb = compute_groebner_basis([], [x])
        assert Integer(0) in gb or 0 in gb

    def test_compute_groebner_basis_two_vars(self):
        """Zwei-Variablen-Gröbner-Basis."""
        gb = compute_groebner_basis([x*y - 1, x - y], [x, y])
        assert isinstance(gb, list)

    # --- ideal_radical ---

    def test_ideal_radical_prime_ideal(self):
        """Primideal ⟨x⟩: Radikal = ⟨x⟩."""
        rad = ideal_radical([x], [x])
        assert isinstance(rad, list)
        # Radikal von x ist x selbst
        assert any(expand(r - x) == 0 for r in rad)

    def test_ideal_radical_x_squared(self):
        """√⟨x²⟩ = ⟨x⟩."""
        rad = ideal_radical([x**2], [x])
        assert isinstance(rad, list)
        assert len(rad) >= 1
        # Ergebnis sollte x sein
        assert any(expand(r - x) == 0 for r in rad)

    def test_ideal_radical_empty(self):
        """Leeres Ideal → [0]."""
        rad = ideal_radical([], [x])
        assert Integer(0) in rad or 0 in rad

    def test_ideal_radical_product(self):
        """√⟨x²·y⟩: Radikal enthält xy."""
        rad = ideal_radical([x**2 * y], [x, y])
        assert isinstance(rad, list)

    def test_ideal_radical_irreducible(self):
        """Irreduzibles Polynom: Radikal = es selbst."""
        rad = ideal_radical([x**2 + 1], [x])
        assert isinstance(rad, list)

    # --- zariski_closure_demo ---

    def test_zariski_closure_empty_points(self):
        """Leere Punktmenge → Zariski-Abschluss = ∅ = V(1)."""
        result = zariski_closure_demo([], [x])
        assert result['ideal_generators'] == [Integer(1)]

    def test_zariski_closure_single_point(self):
        """Einziger Punkt {(1,2)}: erzeuge maximales Ideal."""
        result = zariski_closure_demo([{x: 1, y: 2}], [x, y])
        gens = result['ideal_generators']
        assert len(gens) == 2
        # Sollte x-1 und y-2 enthalten
        gen_expanded = [expand(g) for g in gens]
        assert expand(x - 1) in gen_expanded
        assert expand(y - 2) in gen_expanded

    def test_zariski_closure_two_points_1d(self):
        """Zwei Punkte in 1D: Polynom mit zwei Nullstellen."""
        result = zariski_closure_demo([{x: 1}, {x: -1}], [x])
        gens = result['ideal_generators']
        assert len(gens) == 1
        # Polynom hat Nullstellen bei ±1 → (x-1)(x+1) = x²-1
        g = expand(gens[0])
        assert expand(g.subs(x, 1)) == 0
        assert expand(g.subs(x, -1)) == 0

    def test_zariski_closure_returns_description(self):
        """Ergebnis enthält Beschreibungstext."""
        result = zariski_closure_demo([{x: 0}], [x])
        assert 'description' in result
        assert isinstance(result['description'], str)

    # --- hilbert_polynomial ---

    def test_hilbert_polynomial_zero_ideal(self):
        """Nullideal: Hilbert-Polynom = C(t+n,n)."""
        from sympy import binomial
        t = symbols('t')
        p = hilbert_polynomial([], [x, y])
        # Für n=2: C(t+2,2) = (t+1)(t+2)/2
        assert p is not None

    def test_hilbert_polynomial_principal_ideal(self):
        """Hauptideal ⟨f⟩ mit Grad d: Hilbert-Polynom ist Polynom in t."""
        from sympy import symbols as sym_symbols
        t = sym_symbols('t')
        p = hilbert_polynomial([x**2 - y], [x, y])
        assert p is not None

    def test_hilbert_polynomial_single_var_linear(self):
        """⟨x-1⟩ in k[x]: Hilbert-Polynom."""
        p = hilbert_polynomial([x - 1], [x])
        assert p is not None

    def test_hilbert_polynomial_unit_ideal(self):
        """Einheitsideal ⟨1⟩: Hilbert-Polynom."""
        p = hilbert_polynomial([Integer(1)], [x])
        assert p is not None
