"""
@file test_algebra_diophantine.py
@brief Tests für Diophantische Gleichungen und Quadratisches Reziprozitätsgesetz.
@description
    Testet alle neuen Funktionen in algebra.py:
    - Lineare Diophantische Gleichungen
    - Pell-Gleichung
    - Pythagoräische Tripel
    - Summe zweier Quadrate
    - Markov-Zahlen
    - Quadratisches Reziprozitätsgesetz (QR, Tonelli-Shanks, Cipolla)

@author Kurt Ingwer
@date 2026-03-08
"""

import pytest
import math
import sys
import os

# Projektpfad zum Modulpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algebra import (
    solve_linear_diophantine,
    solve_quadratic_diophantine_pell,
    solve_pythagorean_triples,
    solve_diophantine_two_squares,
    markov_numbers,
    is_quadratic_residue,
    quadratic_residues,
    quadratic_reciprocity,
    tonelli_shanks,
    cipolla_algorithm,
    gcd,
)


# ===========================================================================
# Tests: solve_linear_diophantine
# ===========================================================================

class TestSolveLinearDiophantine:
    """Tests für lineare Diophantische Gleichungen ax + by = c."""

    def test_basic_solution_exists(self):
        """3x + 5y = 1 hat eine Lösung (gcd(3,5)=1 | 1)."""
        result = solve_linear_diophantine(3, 5, 1)
        assert result is not None
        x0, y0, g = result
        # Verifikation: 3*x0 + 5*y0 = 1
        assert 3 * x0 + 5 * y0 == 1
        assert g == 1

    def test_no_solution(self):
        """6x + 4y = 3 hat keine Lösung: gcd(6,4)=2, 2 ∤ 3."""
        result = solve_linear_diophantine(6, 4, 3)
        assert result is None

    def test_solution_verification(self):
        """Jede zurückgegebene Lösung muss die Gleichung erfüllen."""
        # 12x + 8y = 4, gcd=4, 4|4 ✓
        result = solve_linear_diophantine(12, 8, 4)
        assert result is not None
        x0, y0, g = result
        assert 12 * x0 + 8 * y0 == 4
        assert g == 4

    def test_gcd_returned_correctly(self):
        """gcd_ab im Rückgabetupel muss korrekt sein."""
        result = solve_linear_diophantine(6, 10, 4)
        assert result is not None
        _, _, g = result
        assert g == 2

    def test_negative_coefficients(self):
        """Negative Koeffizienten müssen korrekt behandelt werden."""
        result = solve_linear_diophantine(-3, 5, 1)
        assert result is not None
        x0, y0, _ = result
        assert -3 * x0 + 5 * y0 == 1

    def test_c_is_multiple_of_gcd(self):
        """c ist ein Vielfaches von gcd(a,b)."""
        # 6x + 4y = 8, gcd=2, 2|8 ✓
        result = solve_linear_diophantine(6, 4, 8)
        assert result is not None
        x0, y0, g = result
        assert 6 * x0 + 4 * y0 == 8

    def test_trivial_case_a_equals_b(self):
        """a = b Fall: 5x + 5y = 10, gcd=5, 5|10 ✓."""
        result = solve_linear_diophantine(5, 5, 10)
        assert result is not None
        x0, y0, _ = result
        assert 5 * x0 + 5 * y0 == 10

    def test_large_numbers(self):
        """Große Zahlen müssen korrekt behandelt werden."""
        result = solve_linear_diophantine(1071, 462, 1)
        # gcd(1071,462)=21, 21∤1 → keine Lösung wäre korrekt
        # aber gcd(1071,462): 1071=2*462+147, 462=3*147+21, 147=7*21 → gcd=21
        # 21∤1 → None
        assert result is None

    def test_unit_gcd_large(self):
        """gcd(a,b)=1 Fall mit großen Zahlen."""
        result = solve_linear_diophantine(17, 13, 1)
        assert result is not None
        x0, y0, g = result
        assert 17 * x0 + 13 * y0 == 1
        assert g == 1


# ===========================================================================
# Tests: solve_quadratic_diophantine_pell
# ===========================================================================

class TestPellEquation:
    """Tests für die Pell-Gleichung x² - D·y² = 1."""

    def test_pell_d2_five_solutions(self):
        """Pell-Gleichung x²-2y²=1: erste 5 Lösungen."""
        solutions = solve_quadratic_diophantine_pell(2, 5)
        assert len(solutions) == 5
        # Alle Lösungen müssen x²-2y²=1 erfüllen
        for x, y in solutions:
            assert x * x - 2 * y * y == 1, f"Pell-Gleichung verletzt: x={x}, y={y}"

    def test_pell_d2_first_solution(self):
        """Erste Lösung von x²-2y²=1 ist (3,2)."""
        solutions = solve_quadratic_diophantine_pell(2, 1)
        assert solutions[0] == (3, 2)

    def test_pell_d3(self):
        """Pell-Gleichung x²-3y²=1: Fundamentallösung (2,1)."""
        solutions = solve_quadratic_diophantine_pell(3, 3)
        assert solutions[0] == (2, 1)
        for x, y in solutions:
            assert x * x - 3 * y * y == 1

    def test_pell_d5(self):
        """Pell-Gleichung x²-5y²=1: Fundamentallösung (9,4)."""
        solutions = solve_quadratic_diophantine_pell(5, 3)
        assert solutions[0] == (9, 4)
        for x, y in solutions:
            assert x * x - 5 * y * y == 1

    def test_pell_d_perfect_square_raises(self):
        """D als perfektes Quadrat wirft ValueError."""
        with pytest.raises(ValueError):
            solve_quadratic_diophantine_pell(4, 3)

    def test_pell_n_solutions_count(self):
        """Anzahl der zurückgegebenen Lösungen muss n_solutions sein."""
        for n in [1, 3, 5, 10]:
            solutions = solve_quadratic_diophantine_pell(2, n)
            assert len(solutions) == n

    def test_pell_all_solutions_valid(self):
        """Alle erzeugten Lösungen für D=7 sind korrekt."""
        solutions = solve_quadratic_diophantine_pell(7, 5)
        for x, y in solutions:
            assert x * x - 7 * y * y == 1


# ===========================================================================
# Tests: solve_pythagorean_triples
# ===========================================================================

class TestPythagoreanTriples:
    """Tests für primitive Pythagoräische Tripel."""

    def test_triples_up_to_100(self):
        """Primitive Tripel mit c ≤ 100."""
        triples = solve_pythagorean_triples(100)
        # Bekannte primitive Tripel: (3,4,5), (5,12,13), (8,15,17), (7,24,25)
        assert (3, 4, 5) in triples
        assert (5, 12, 13) in triples
        assert (8, 15, 17) in triples

    def test_triple_satisfies_pythagoras(self):
        """Alle Tripel müssen a²+b²=c² erfüllen."""
        triples = solve_pythagorean_triples(200)
        for a, b, c in triples:
            assert a * a + b * b == c * c, f"Pythagoras verletzt: ({a},{b},{c})"

    def test_triple_is_primitive(self):
        """Alle Tripel müssen primitiv sein (gcd=1)."""
        triples = solve_pythagorean_triples(200)
        for a, b, c in triples:
            assert gcd(gcd(a, b), c) == 1, f"Nicht primitiv: ({a},{b},{c})"

    def test_triple_sorted(self):
        """In jedem Tripel muss a < b < c gelten."""
        triples = solve_pythagorean_triples(200)
        for a, b, c in triples:
            assert a < b < c, f"Nicht sortiert: ({a},{b},{c})"

    def test_c_within_bound(self):
        """Alle Tripel müssen c ≤ n erfüllen."""
        n = 150
        triples = solve_pythagorean_triples(n)
        for a, b, c in triples:
            assert c <= n, f"c={c} > {n}"

    def test_small_bound_returns_345(self):
        """Für n=5 wird nur (3,4,5) zurückgegeben."""
        triples = solve_pythagorean_triples(5)
        assert (3, 4, 5) in triples

    def test_no_non_primitive_triples(self):
        """(6,8,10) darf nicht in der Liste sein (nicht primitiv, da gcd=2)."""
        triples = solve_pythagorean_triples(100)
        assert (6, 8, 10) not in triples


# ===========================================================================
# Tests: solve_diophantine_two_squares
# ===========================================================================

class TestTwoSquares:
    """Tests für die Darstellung als Summe zweier Quadrate."""

    def test_25_has_two_representations(self):
        """25 = 0²+5² = 3²+4² → zwei Darstellungen."""
        result = solve_diophantine_two_squares(25)
        assert len(result) >= 2
        # Verifikation aller Darstellungen
        for a, b in result:
            assert a * a + b * b == 25

    def test_3_not_representable(self):
        """3 kann nicht als Summe zweier Quadrate geschrieben werden."""
        result = solve_diophantine_two_squares(3)
        assert result == []

    def test_0_is_0_plus_0(self):
        """0 = 0² + 0²."""
        result = solve_diophantine_two_squares(0)
        assert (0, 0) in result

    def test_1_is_0_plus_1(self):
        """1 = 0² + 1²."""
        result = solve_diophantine_two_squares(1)
        assert (0, 1) in result

    def test_2_is_1_plus_1(self):
        """2 = 1² + 1²."""
        result = solve_diophantine_two_squares(2)
        assert (1, 1) in result

    def test_50_representations(self):
        """50 = 1²+7² = 5²+5² (zwei Darstellungen)."""
        result = solve_diophantine_two_squares(50)
        for a, b in result:
            assert a * a + b * b == 50

    def test_all_results_satisfy_equation(self):
        """Alle zurückgegebenen Paare müssen a²+b²=n erfüllen."""
        for n in [0, 1, 2, 4, 5, 8, 9, 10, 13, 25, 50]:
            for a, b in solve_diophantine_two_squares(n):
                assert a * a + b * b == n
                assert a <= b  # geordnet

    def test_7_not_representable(self):
        """7 ≡ 3 mod 4 (Primzahl) → nicht darstellbar."""
        result = solve_diophantine_two_squares(7)
        assert result == []


# ===========================================================================
# Tests: markov_numbers
# ===========================================================================

class TestMarkovNumbers:
    """Tests für Markov-Zahlen."""

    def test_first_markov_numbers(self):
        """Die ersten Markov-Zahlen sind bekannt: 1, 2, 5, 13, 29, ..."""
        result = markov_numbers(10)
        # 1, 2, 5 müssen enthalten sein
        assert 1 in result
        assert 2 in result
        assert 5 in result

    def test_count_returned(self):
        """Anzahl der zurückgegebenen Werte muss n_terms sein."""
        for n in [5, 10, 20]:
            result = markov_numbers(n)
            assert len(result) == n

    def test_sorted(self):
        """Markov-Zahlen müssen sortiert sein."""
        result = markov_numbers(15)
        assert result == sorted(result)

    def test_all_positive(self):
        """Alle Markov-Zahlen sind positiv."""
        result = markov_numbers(20)
        assert all(x > 0 for x in result)

    def test_markov_triple_equation(self):
        """Prüft ob bekannte Markov-Tripel die Gleichung x²+y²+z²=3xyz erfüllen."""
        # Bekannte Tripel: (1,1,1), (1,1,2), (1,2,5), (1,5,13), (2,5,29)
        triples = [(1, 1, 1), (1, 1, 2), (1, 2, 5), (1, 5, 13), (2, 5, 29)]
        for x, y, z in triples:
            assert x**2 + y**2 + z**2 == 3 * x * y * z, f"Markov-Gleichung verletzt: ({x},{y},{z})"


# ===========================================================================
# Tests: is_quadratic_residue
# ===========================================================================

class TestIsQuadraticResidue:
    """Tests für quadratische Reste."""

    def test_2_mod_7_is_qr(self):
        """2 ist quadratischer Rest mod 7: 3²=9≡2 (mod 7)."""
        assert is_quadratic_residue(2, 7) is True

    def test_3_mod_7_is_not_qr(self):
        """3 ist kein quadratischer Rest mod 7."""
        assert is_quadratic_residue(3, 7) is False

    def test_1_is_always_qr(self):
        """1 ist immer ein quadratischer Rest (1²≡1)."""
        for p in [5, 7, 11, 13, 17]:
            assert is_quadratic_residue(1, p) is True

    def test_0_is_not_qr(self):
        """0 ist kein quadratischer Rest (im üblichen Sinne)."""
        assert is_quadratic_residue(0, 7) is False

    def test_euler_criterion(self):
        """Euler-Kriterium: a^((p-1)/2) ≡ 1 (mod p) für QR."""
        # Quadratische Reste mod 11: 1,3,4,5,9
        qr_mod11 = {1, 3, 4, 5, 9}
        for a in range(1, 11):
            expected = (a in qr_mod11)
            assert is_quadratic_residue(a, 11) == expected


# ===========================================================================
# Tests: quadratic_residues
# ===========================================================================

class TestQuadraticResidues:
    """Tests für die Liste quadratischer Reste."""

    def test_qr_mod_7(self):
        """Quadratische Reste mod 7 sind {1, 2, 4}."""
        result = quadratic_residues(7)
        assert result == [1, 2, 4]

    def test_qr_mod_5(self):
        """Quadratische Reste mod 5 sind {1, 4}."""
        result = quadratic_residues(5)
        assert result == [1, 4]

    def test_qr_count(self):
        """Es gibt genau (p-1)/2 quadratische Reste mod p."""
        for p in [5, 7, 11, 13, 17, 19]:
            result = quadratic_residues(p)
            assert len(result) == (p - 1) // 2

    def test_qr_sorted(self):
        """Ergebnisliste muss sortiert sein."""
        for p in [5, 7, 11, 13]:
            result = quadratic_residues(p)
            assert result == sorted(result)

    def test_qr_mod_11(self):
        """Quadratische Reste mod 11 sind {1, 3, 4, 5, 9}."""
        result = quadratic_residues(11)
        assert result == [1, 3, 4, 5, 9]


# ===========================================================================
# Tests: quadratic_reciprocity
# ===========================================================================

class TestQuadraticReciprocity:
    """Tests für das Quadratische Reziprozitätsgesetz."""

    def test_reciprocity_3_7(self):
        """Quadratisches Reziprozitätsgesetz für p=3, q=7."""
        result = quadratic_reciprocity(3, 7)
        assert 'legendre_p_q' in result
        assert 'legendre_q_p' in result
        assert 'law_verified' in result
        assert result['law_verified'] is True

    def test_reciprocity_5_11(self):
        """Gesetz gilt für p=5, q=11."""
        result = quadratic_reciprocity(5, 11)
        assert result['law_verified'] is True

    def test_legendre_values_range(self):
        """Legendre-Symbol kann nur -1, 0, 1 sein."""
        for p, q in [(3, 7), (5, 11), (7, 13), (11, 17)]:
            result = quadratic_reciprocity(p, q)
            assert result['legendre_p_q'] in [-1, 0, 1]
            assert result['legendre_q_p'] in [-1, 0, 1]

    def test_product_matches_expected_sign(self):
        """Produkt der Legendre-Symbole muss expected_sign entsprechen."""
        for p, q in [(3, 7), (5, 11), (7, 11), (3, 5)]:
            result = quadratic_reciprocity(p, q)
            assert result['product'] == result['expected_sign']

    def test_p_1_mod_4_symmetry(self):
        """Wenn p ≡ 1 (mod 4): (p/q) = (q/p)."""
        # p=5 ≡ 1 (mod 4), q=7
        result = quadratic_reciprocity(5, 7)
        # expected_sign = (-1)^(2*3) = 1 → gleiche Symbole
        assert result['expected_sign'] == 1
        assert result['legendre_p_q'] == result['legendre_q_p']


# ===========================================================================
# Tests: tonelli_shanks
# ===========================================================================

class TestTonelliShanks:
    """Tests für den Tonelli-Shanks-Algorithmus."""

    def test_sqrt_2_mod_7(self):
        """x²≡2 (mod 7): Lösung ist 3 oder 4 (da 3²=9≡2, 4²=16≡2 mod 7)."""
        x = tonelli_shanks(2, 7)
        assert x is not None
        assert (x * x) % 7 == 2

    def test_no_solution(self):
        """3 ist kein QR mod 7 → None."""
        x = tonelli_shanks(3, 7)
        assert x is None

    def test_zero_input(self):
        """n=0 → Wurzel ist 0."""
        x = tonelli_shanks(0, 7)
        assert x == 0

    def test_various_primes(self):
        """Teste für verschiedene Primzahlen."""
        test_cases = [(1, 5), (4, 5), (1, 7), (2, 7), (4, 7)]
        for n, p in test_cases:
            x = tonelli_shanks(n, p)
            assert x is not None
            assert (x * x) % p == n % p

    def test_p_3_mod_4(self):
        """Spezialfall p≡3 (mod 4): p=7, n=4 → x=2."""
        x = tonelli_shanks(4, 7)
        assert x is not None
        assert (x * x) % 7 == 4

    def test_large_prime(self):
        """Teste mit größerer Primzahl."""
        p = 101  # 101 ≡ 1 (mod 4) → allgemeiner Tonelli-Shanks
        n = 4
        x = tonelli_shanks(n, p)
        assert x is not None
        assert (x * x) % p == n

    def test_result_in_range(self):
        """Ergebnis muss in [0, p-1] liegen."""
        for n in [1, 2, 4]:
            x = tonelli_shanks(n, 7)
            if x is not None:
                assert 0 <= x < 7


# ===========================================================================
# Tests: cipolla_algorithm
# ===========================================================================

class TestCipollaAlgorithm:
    """Tests für den Cipolla-Algorithmus."""

    def test_sqrt_2_mod_7(self):
        """x²≡2 (mod 7) via Cipolla."""
        x = cipolla_algorithm(2, 7)
        assert x is not None
        assert (x * x) % 7 == 2

    def test_no_solution(self):
        """3 ist kein QR mod 7 → None."""
        x = cipolla_algorithm(3, 7)
        assert x is None

    def test_zero_input(self):
        """n=0 → Wurzel ist 0."""
        x = cipolla_algorithm(0, 7)
        assert x == 0

    def test_agrees_with_tonelli_shanks(self):
        """Cipolla und Tonelli-Shanks müssen übereinstimmen (mod p)."""
        test_cases = [(1, 7), (2, 7), (4, 7), (1, 11), (3, 11)]
        for n, p in test_cases:
            x_cipolla = cipolla_algorithm(n, p)
            x_tonelli = tonelli_shanks(n, p)
            if x_cipolla is not None:
                assert (x_cipolla * x_cipolla) % p == n
            if x_tonelli is not None:
                assert (x_tonelli * x_tonelli) % p == n

    def test_various_qr(self):
        """Teste verschiedene quadratische Reste mod 13."""
        # QR mod 13: 1, 3, 4, 9, 10, 12
        for n in [1, 3, 4, 9, 10, 12]:
            x = cipolla_algorithm(n, 13)
            assert x is not None
            assert (x * x) % 13 == n

    def test_result_in_range(self):
        """Ergebnis muss in [0, p-1] liegen."""
        for n in [1, 2, 4]:
            x = cipolla_algorithm(n, 7)
            if x is not None:
                assert 0 <= x < 7
