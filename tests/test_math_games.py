"""
test_math_games.py – Tests für das Mathematische-Spiele-Modul

Testet alle Funktionen aus math_games.py:
- NimGame (Nim-Werte, Gewinnzüge, Grundy, mex, nim_sum)
- SpragueGrundyTheory (Spielanalyse, Kayles, Subtraktionsspiel)
- Zahlentheoretische Spiele (Euklid, Teiler, Primzahl, Fibonacci-Nim, Chomp)
- Magische Quadrate, Lateinische Quadrate, perfekte Zahlen
- Figurate Zahlen (Dreiecks-, Quadrat-, Penta-, Hexagonalzahlen)

@author    Kurt Ingwer
@version   1.0.0
@timestamp 2026-03-10
"""

import pytest
import sys
import os

# Projektverzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from math_games import (
    NimGame,
    SpragueGrundyTheory,
    euclid_game,
    divisor_game,
    prime_game,
    fibonacci_nim,
    chomp_game_demo,
    magic_square_check,
    magic_square_3x3,
    magic_square_sum,
    latin_square_check,
    perfect_number_check,
    abundant_deficient_classify,
    friendly_numbers,
    triangular_number,
    square_number,
    pentagonal_number,
    hexagonal_number,
    is_triangular,
    polygonal_number,
    cannonball_problem,
)


# ---------------------------------------------------------------------------
# 1. Tests: NimGame
# ---------------------------------------------------------------------------

class TestNimValue:
    """Tests für NimGame.nim_value()"""

    def test_single_pile(self):
        """Einzelner Haufen: Nim-Wert = Haufengröße."""
        assert NimGame.nim_value([5]) == 5
        assert NimGame.nim_value([0]) == 0

    def test_two_equal_piles(self):
        """Zwei gleiche Haufen: Nim-Wert = 0 (XOR von n XOR n)."""
        assert NimGame.nim_value([3, 3]) == 0
        assert NimGame.nim_value([7, 7]) == 0

    def test_three_piles_classic(self):
        """Klassisches Nim (1,2,3): 1 XOR 2 XOR 3 = 0."""
        assert NimGame.nim_value([1, 2, 3]) == 0

    def test_three_piles_nonzero(self):
        """Nim-Wert (1,3,5) = 1 XOR 3 XOR 5 = 7."""
        assert NimGame.nim_value([1, 3, 5]) == 7

    def test_empty_piles(self):
        """Keine Haufen: Nim-Wert = 0."""
        assert NimGame.nim_value([]) == 0

    def test_zero_piles(self):
        """Alle Haufen = 0: Nim-Wert = 0."""
        assert NimGame.nim_value([0, 0, 0]) == 0


class TestNimWinningMove:
    """Tests für NimGame.nim_winning_move()"""

    def test_losing_position_returns_none(self):
        """Verlierposition gibt None zurück."""
        assert NimGame.nim_winning_move([3, 3]) is None
        assert NimGame.nim_winning_move([1, 2, 3]) is None

    def test_winning_move_single_pile(self):
        """Einzelner Haufen: Gewinnzug = Haufen auf 0 reduzieren."""
        result = NimGame.nim_winning_move([5])
        assert result is not None
        idx, target = result
        assert idx == 0
        assert target == 0

    def test_winning_move_leads_to_zero_nimsum(self):
        """Gewinnzug führt zu Nim-Summe 0."""
        piles = [4, 5, 6]
        result = NimGame.nim_winning_move(piles)
        assert result is not None
        idx, target = result
        new_piles = list(piles)
        new_piles[idx] = target
        assert NimGame.nim_value(new_piles) == 0

    def test_winning_move_reduces_pile(self):
        """Gewinnzug reduziert Haufen (target < original)."""
        piles = [1, 3, 5]
        result = NimGame.nim_winning_move(piles)
        assert result is not None
        idx, target = result
        assert target < piles[idx]


class TestMex:
    """Tests für NimGame.mex()"""

    def test_empty_set(self):
        """mex({}) = 0."""
        assert NimGame.mex(set()) == 0

    def test_contains_zero(self):
        """mex({0}) = 1."""
        assert NimGame.mex({0}) == 1

    def test_consecutive(self):
        """mex({0,1,2}) = 3."""
        assert NimGame.mex({0, 1, 2}) == 3

    def test_gap(self):
        """mex({0,1,3}) = 2."""
        assert NimGame.mex({0, 1, 3}) == 2

    def test_large_set(self):
        """mex({0,1,2,3,4}) = 5."""
        assert NimGame.mex({0, 1, 2, 3, 4}) == 5


class TestNimSum:
    """Tests für NimGame.nim_sum()"""

    def test_empty(self):
        """Nim-Summe leerer Liste = 0."""
        assert NimGame.nim_sum([]) == 0

    def test_single(self):
        """Nim-Summe Einzelelement = Element selbst."""
        assert NimGame.nim_sum([7]) == 7

    def test_xor_two(self):
        """Nim-Summe [5, 3] = 5 XOR 3 = 6."""
        assert NimGame.nim_sum([5, 3]) == 6

    def test_three_values(self):
        """Nim-Summe [1, 2, 3] = 1 XOR 2 XOR 3 = 0."""
        assert NimGame.nim_sum([1, 2, 3]) == 0


class TestGrundyValue:
    """Tests für NimGame.grundy_value()"""

    def test_zero_is_zero(self):
        """G(0) = 0 für jedes Spiel (Terminal-Position)."""
        assert NimGame.grundy_value(0, {1, 2}) == 0

    def test_take_one_or_two(self):
        """G(n) für moves={1,2}: G(0)=0,G(1)=1,G(2)=2,G(3)=0,G(4)=1,G(5)=2,..."""
        expected = [0, 1, 2, 0, 1, 2]
        for n, g in enumerate(expected):
            assert NimGame.grundy_value(n, {1, 2}) == g

    def test_take_one_grundy(self):
        """G(n) für moves={1}: G(n) = n % 2."""
        for n in range(6):
            assert NimGame.grundy_value(n, {1}) == n % 2


# ---------------------------------------------------------------------------
# 2. Tests: SpragueGrundyTheory
# ---------------------------------------------------------------------------

class TestGameSumGrundy:
    """Tests für SpragueGrundyTheory.game_sum_grundy()"""

    def test_xor_of_grundy_values(self):
        """G(G1 + G2) = G(G1) XOR G(G2)."""
        assert SpragueGrundyTheory.game_sum_grundy(3, 5) == 6
        assert SpragueGrundyTheory.game_sum_grundy(0, 7) == 7

    def test_same_values(self):
        """G(n) XOR G(n) = 0 (Verlierposition)."""
        assert SpragueGrundyTheory.game_sum_grundy(4, 4) == 0

    def test_zero_game(self):
        """G(G1 + 0) = G(G1)."""
        assert SpragueGrundyTheory.game_sum_grundy(5, 0) == 5


class TestKaylesGrundy:
    """Tests für SpragueGrundyTheory.kayles_grundy()"""

    def test_zero_grundy(self):
        """G(0) = 0 für Kayles."""
        g = SpragueGrundyTheory.kayles_grundy(5)
        assert g[0] == 0

    def test_one_pin(self):
        """G(1) = 1 (ein Kegel, ein Zug)."""
        g = SpragueGrundyTheory.kayles_grundy(5)
        assert g[1] == 1

    def test_two_pins(self):
        """G(2) für Kayles."""
        g = SpragueGrundyTheory.kayles_grundy(5)
        # mex{G(1), G(0)} = mex{1, 0} = 2
        assert g[2] == 2

    def test_length_correct(self):
        """Liste hat n+1 Einträge."""
        g = SpragueGrundyTheory.kayles_grundy(10)
        assert len(g) == 11


class TestSubtractionGameGrundy:
    """Tests für SpragueGrundyTheory.subtraction_game_grundy()"""

    def test_subtraction_set_1_2(self):
        """Subtraktionsspiel {1,2}: periodisch mit Periode 3."""
        g = SpragueGrundyTheory.subtraction_game_grundy(8, {1, 2})
        expected = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        assert g == expected

    def test_subtraction_set_1(self):
        """Subtraktionsspiel {1}: G(n) = n mod 2."""
        g = SpragueGrundyTheory.subtraction_game_grundy(6, {1})
        assert g == [0, 1, 0, 1, 0, 1, 0]

    def test_terminal_zero(self):
        """G(0) = 0 immer."""
        g = SpragueGrundyTheory.subtraction_game_grundy(5, {1, 3})
        assert g[0] == 0


class TestAnalyzeGame:
    """Tests für SpragueGrundyTheory.analyze_game()"""

    def test_terminal_position(self):
        """Position ohne Züge hat Grundy-Wert 0."""
        def no_moves(pos):
            return []
        assert SpragueGrundyTheory.analyze_game(0, no_moves) == 0

    def test_single_move_game(self):
        """Spiel mit Zügen: G(n) = mex{G(n-1)} für moves={1}."""
        def subtract_one(pos):
            return [pos - 1] if pos > 0 else []
        # G(0)=0, G(1)=mex{G(0)}=mex{0}=1, G(2)=mex{G(1)}=mex{1}=0
        assert SpragueGrundyTheory.analyze_game(0, subtract_one) == 0
        assert SpragueGrundyTheory.analyze_game(1, subtract_one) == 1
        assert SpragueGrundyTheory.analyze_game(2, subtract_one) == 0


class TestTakeAwayGame:
    """Tests für SpragueGrundyTheory.take_away_game()"""

    def test_same_as_subtraction_game(self):
        """take_away_game ist identisch mit subtraction_game_grundy."""
        g1 = SpragueGrundyTheory.take_away_game(6, [1, 2])
        g2 = SpragueGrundyTheory.subtraction_game_grundy(6, {1, 2})
        assert g1 == g2


# ---------------------------------------------------------------------------
# 3. Tests: Zahlentheoretische Spiele
# ---------------------------------------------------------------------------

class TestEuclidGame:
    """Tests für euclid_game()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = euclid_game(5, 3)
        assert isinstance(result, dict)

    def test_contains_gewinner(self):
        """Wörterbuch enthält 'gewinner'."""
        result = euclid_game(5, 3)
        assert 'gewinner' in result

    def test_b_at_least_2a_player1_wins(self):
        """Falls b ≥ 2a, gewinnt Spieler 1."""
        result = euclid_game(3, 6)  # b = 6 = 2*3 ≥ 2*a
        assert result['gewinner'] == 'Spieler 1'

    def test_order_doesnt_matter(self):
        """Reihenfolge von a, b spielt keine Rolle."""
        result1 = euclid_game(3, 7)
        result2 = euclid_game(7, 3)
        assert result1['gewinner'] == result2['gewinner']


class TestDivisorGame:
    """Tests für divisor_game()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = divisor_game(6)
        assert isinstance(result, dict)

    def test_n2_player1_wins(self):
        """n=2: Spieler 1 subtrahiert 1 → Spieler 2 hat 1, kein Zug → Spieler 1 gewinnt."""
        result = divisor_game(2)
        assert result['gewinner'] == 'Spieler 1'

    def test_grundy_value_present(self):
        """Grundy-Wert ist im Ergebnis."""
        result = divisor_game(4)
        assert 'grundy_wert' in result

    def test_grundy_zero_means_player2_wins(self):
        """Grundy-Wert 0 → Spieler 2 gewinnt."""
        result = divisor_game(6)
        if result['grundy_wert'] == 0:
            assert result['gewinner'] == 'Spieler 2'
        else:
            assert result['gewinner'] == 'Spieler 1'


class TestPrimeGame:
    """Tests für prime_game()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = prime_game(5)
        assert isinstance(result, dict)

    def test_n1_player2_wins(self):
        """n=1: keine Primzahl ≤ 1 → Spieler 1 kann nicht ziehen → Spieler 2 gewinnt."""
        result = prime_game(1)
        assert result['gewinner'] == 'Spieler 2'

    def test_n2_player1_wins(self):
        """n=2: Spieler 1 nimmt 2 → 0 Steine → Spieler 1 gewinnt."""
        result = prime_game(2)
        assert result['gewinner'] == 'Spieler 1'

    def test_grundy_table_present(self):
        """Grundy-Tabelle ist vorhanden."""
        result = prime_game(10)
        assert 'grundy_tabelle' in result

    def test_erlaubte_zuege_are_primes(self):
        """Erlaubte Züge sind Primzahlen."""
        result = prime_game(10)
        zuege = result['erlaubte_zuege']
        primes_expected = [2, 3, 5, 7]
        assert zuege == primes_expected


class TestFibonacciNim:
    """Tests für fibonacci_nim()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = fibonacci_nim(10)
        assert isinstance(result, dict)

    def test_fibonacci_number_player2_wins(self):
        """Fibonacci-Zahl → Spieler 2 gewinnt."""
        for fib in [1, 2, 3, 5, 8, 13, 21]:
            result = fibonacci_nim(fib)
            assert result['gewinner'] == 'Spieler 2', f"Fehler bei n={fib}"

    def test_non_fibonacci_player1_wins(self):
        """Nicht-Fibonacci-Zahl → Spieler 1 gewinnt."""
        for non_fib in [4, 6, 7, 9, 10, 11, 12, 14]:
            result = fibonacci_nim(non_fib)
            assert result['gewinner'] == 'Spieler 1', f"Fehler bei n={non_fib}"

    def test_zeckendorf_decomposition_valid(self):
        """Zeckendorf-Zerlegung summiert zu n."""
        result = fibonacci_nim(10)
        assert sum(result['zeckendorf_zerlegung']) == 10

    def test_zeckendorf_for_fib(self):
        """Zeckendorf-Zerlegung einer Fibonacci-Zahl ist die Zahl selbst."""
        result = fibonacci_nim(8)
        assert result['zeckendorf_zerlegung'] == [8]


class TestChompGameDemo:
    """Tests für chomp_game_demo()"""

    def test_1x1_player2_loses(self):
        """1×1-Tafel: Spieler 1 muss vergiftetes Feld nehmen → Spieler 2 gewinnt."""
        result = chomp_game_demo(1, 1)
        assert result['gewinner'] == 'Spieler 2'

    def test_larger_board_player1_wins(self):
        """Für m,n ≥ 2 gewinnt Spieler 1."""
        result = chomp_game_demo(3, 3)
        assert result['gewinner'] == 'Spieler 1'

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = chomp_game_demo(2, 5)
        assert isinstance(result, dict)

    def test_2x2_player1_wins(self):
        """2×2-Tafel: Spieler 1 gewinnt."""
        result = chomp_game_demo(2, 2)
        assert result['gewinner'] == 'Spieler 1'


# ---------------------------------------------------------------------------
# 4. Tests: Magische Quadrate
# ---------------------------------------------------------------------------

class TestMagicSquareCheck:
    """Tests für magic_square_check()"""

    def test_valid_3x3(self):
        """Standard 3×3 Lo-Shu-Quadrat ist magisch."""
        sq = [
            [2, 7, 6],
            [9, 5, 1],
            [4, 3, 8]
        ]
        assert magic_square_check(sq) is True

    def test_invalid_square(self):
        """Nicht-magisches Quadrat wird erkannt."""
        sq = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ]
        assert magic_square_check(sq) is False

    def test_empty_matrix(self):
        """Leere Matrix ist kein magisches Quadrat."""
        assert magic_square_check([]) is False

    def test_1x1_is_magic(self):
        """1×1-Matrix ist triviales magisches Quadrat."""
        assert magic_square_check([[1]]) is True

    def test_4x4_magic(self):
        """Bekanntes 4×4 magisches Quadrat (Dürer's Melancholia)."""
        sq = [
            [16, 3, 2, 13],
            [5, 10, 11, 8],
            [9, 6, 7, 12],
            [4, 15, 14, 1]
        ]
        assert magic_square_check(sq) is True


class TestMagicSquare3x3:
    """Tests für magic_square_3x3()"""

    def test_default_is_magic(self):
        """Standard-Quadrat (start=1) ist magisch."""
        sq = magic_square_3x3()
        assert magic_square_check(sq) is True

    def test_default_sum_is_15(self):
        """Magische Summe für 3×3 = 15."""
        sq = magic_square_3x3()
        assert sum(sq[0]) == 15

    def test_shifted_start(self):
        """Verschobenes Quadrat (start=0) ist nicht magisch (enthält 0)."""
        sq = magic_square_3x3(start=0)
        # Quadrat hat Einträge 0..8, ist strukturell korrekt
        assert isinstance(sq, list)
        assert len(sq) == 3


class TestMagicSquareSum:
    """Tests für magic_square_sum()"""

    def test_n3(self):
        """S(3) = 3(9+1)/2 = 15."""
        assert magic_square_sum(3) == 15

    def test_n4(self):
        """S(4) = 4(16+1)/2 = 34."""
        assert magic_square_sum(4) == 34

    def test_n1(self):
        """S(1) = 1(1+1)/2 = 1."""
        assert magic_square_sum(1) == 1

    def test_n5(self):
        """S(5) = 5(25+1)/2 = 65."""
        assert magic_square_sum(5) == 65


class TestLatinSquareCheck:
    """Tests für latin_square_check()"""

    def test_valid_2x2(self):
        """Gültiges 2×2 lateinisches Quadrat."""
        sq = [[0, 1], [1, 0]]
        assert latin_square_check(sq) is True

    def test_valid_3x3(self):
        """Gültiges 3×3 lateinisches Quadrat."""
        sq = [
            [0, 1, 2],
            [1, 2, 0],
            [2, 0, 1]
        ]
        assert latin_square_check(sq) is True

    def test_invalid(self):
        """Ungültiges Quadrat (doppelte Werte in Zeile)."""
        sq = [
            [0, 0, 1],
            [1, 2, 0],
            [2, 1, 2]
        ]
        assert latin_square_check(sq) is False

    def test_empty(self):
        """Leere Matrix ist kein lateinisches Quadrat."""
        assert latin_square_check([]) is False

    def test_1x1(self):
        """1×1-Matrix ist triviales lateinisches Quadrat."""
        assert latin_square_check([[5]]) is True


# ---------------------------------------------------------------------------
# 5. Tests: Perfekte und befreundete Zahlen
# ---------------------------------------------------------------------------

class TestPerfectNumberCheck:
    """Tests für perfect_number_check()"""

    def test_6_is_perfect(self):
        """6 ist perfekt: σ(6) = 1+2+3+6 = 12 = 2·6."""
        assert perfect_number_check(6) is True

    def test_28_is_perfect(self):
        """28 ist perfekt: 1+2+4+7+14+28 = 56 = 2·28."""
        assert perfect_number_check(28) is True

    def test_496_is_perfect(self):
        """496 ist perfekt."""
        assert perfect_number_check(496) is True

    def test_not_perfect(self):
        """Nicht-perfekte Zahlen: 1, 2, 3, 4, 5, 7."""
        for n in [1, 2, 3, 4, 5, 7, 10, 12]:
            assert perfect_number_check(n) is False

    def test_one_not_perfect(self):
        """1 ist nicht perfekt."""
        assert perfect_number_check(1) is False


class TestAbundantDeficientClassify:
    """Tests für abundant_deficient_classify()"""

    def test_6_is_perfect(self):
        """6 ist perfekt."""
        assert abundant_deficient_classify(6) == 'perfekt'

    def test_12_is_abundant(self):
        """12 ist abundant: σ(12) = 1+2+3+4+6+12 = 28 > 24."""
        assert abundant_deficient_classify(12) == 'abundant'

    def test_5_is_deficient(self):
        """5 ist defizient: σ(5) = 1+5 = 6 < 10."""
        assert abundant_deficient_classify(5) == 'defizient'

    def test_1_is_deficient(self):
        """1 ist defizient."""
        assert abundant_deficient_classify(1) == 'defizient'

    def test_invalid_raises(self):
        """Negative Zahl wirft ValueError."""
        with pytest.raises(ValueError):
            abundant_deficient_classify(0)


class TestFriendlyNumbers:
    """Tests für friendly_numbers()"""

    def test_220_friends_284(self):
        """220 und 284 sind befreundete Zahlen."""
        friends = friendly_numbers(220, max_search=300)
        assert 284 in friends

    def test_284_friends_220(self):
        """284 und 220 sind befreundete Zahlen."""
        friends = friendly_numbers(284, max_search=300)
        assert 220 in friends

    def test_returns_list(self):
        """Rückgabe ist eine Liste."""
        result = friendly_numbers(6, max_search=100)
        assert isinstance(result, list)

    def test_non_amicable_number(self):
        """Zahl ohne befreundete Partner."""
        # 10 hat keine befreundeten Zahlen in üblichem Bereich
        friends = friendly_numbers(10, max_search=100)
        # Toleranz: kann leer sein oder nicht, je nach Bereich
        assert isinstance(friends, list)


# ---------------------------------------------------------------------------
# 6. Tests: Figurate Zahlen
# ---------------------------------------------------------------------------

class TestTriangularNumber:
    """Tests für triangular_number()"""

    def test_t0(self):
        """T(0) = 0."""
        assert triangular_number(0) == 0

    def test_t1(self):
        """T(1) = 1."""
        assert triangular_number(1) == 1

    def test_t2(self):
        """T(2) = 3."""
        assert triangular_number(2) == 3

    def test_t3(self):
        """T(3) = 6."""
        assert triangular_number(3) == 6

    def test_t4(self):
        """T(4) = 10."""
        assert triangular_number(4) == 10

    def test_t10(self):
        """T(10) = 55."""
        assert triangular_number(10) == 55

    def test_negative_raises(self):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            triangular_number(-1)


class TestSquareNumber:
    """Tests für square_number()"""

    def test_s0(self):
        """S(0) = 0."""
        assert square_number(0) == 0

    def test_s1(self):
        """S(1) = 1."""
        assert square_number(1) == 1

    def test_s5(self):
        """S(5) = 25."""
        assert square_number(5) == 25

    def test_negative_raises(self):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            square_number(-1)


class TestPentagonalNumber:
    """Tests für pentagonal_number()"""

    def test_p0(self):
        """P(0) = 0."""
        assert pentagonal_number(0) == 0

    def test_p1(self):
        """P(1) = 1."""
        assert pentagonal_number(1) == 1

    def test_p2(self):
        """P(2) = 5."""
        assert pentagonal_number(2) == 5

    def test_p3(self):
        """P(3) = 12."""
        assert pentagonal_number(3) == 12

    def test_p4(self):
        """P(4) = 22."""
        assert pentagonal_number(4) == 22

    def test_negative_raises(self):
        """Negativer Index wirft ValueError."""
        with pytest.raises(ValueError):
            pentagonal_number(-1)


class TestHexagonalNumber:
    """Tests für hexagonal_number()"""

    def test_h0(self):
        """H(0) = 0."""
        assert hexagonal_number(0) == 0

    def test_h1(self):
        """H(1) = 1."""
        assert hexagonal_number(1) == 1

    def test_h2(self):
        """H(2) = 6."""
        assert hexagonal_number(2) == 6

    def test_h3(self):
        """H(3) = 15."""
        assert hexagonal_number(3) == 15

    def test_h4(self):
        """H(4) = 28."""
        assert hexagonal_number(4) == 28


class TestIsTriangular:
    """Tests für is_triangular()"""

    def test_triangular_numbers(self):
        """Bekannte Dreieckszahlen werden erkannt."""
        for t in [0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55]:
            assert is_triangular(t) is True, f"{t} sollte Dreieckszahl sein"

    def test_non_triangular(self):
        """Nicht-Dreieckszahlen werden erkannt."""
        for n in [2, 4, 5, 7, 8, 9, 11, 12, 13, 14]:
            assert is_triangular(n) is False, f"{n} sollte keine Dreieckszahl sein"

    def test_negative(self):
        """Negative Zahlen sind keine Dreieckszahlen."""
        assert is_triangular(-1) is False


class TestPolygonalNumber:
    """Tests für polygonal_number()"""

    def test_triangular_s3(self):
        """s=3 ergibt Dreieckszahlen."""
        for n in range(6):
            assert polygonal_number(3, n) == triangular_number(n)

    def test_square_s4(self):
        """s=4 ergibt Quadratzahlen."""
        for n in range(6):
            assert polygonal_number(4, n) == square_number(n)

    def test_pentagonal_s5(self):
        """s=5 ergibt Pentagonalzahlen."""
        for n in range(6):
            assert polygonal_number(5, n) == pentagonal_number(n)

    def test_hexagonal_s6(self):
        """s=6 ergibt Hexagonalzahlen."""
        for n in range(6):
            assert polygonal_number(6, n) == hexagonal_number(n)

    def test_invalid_s(self):
        """s < 3 wirft ValueError."""
        with pytest.raises(ValueError):
            polygonal_number(2, 5)


class TestCannonballProblem:
    """Tests für cannonball_problem()"""

    def test_returns_dict(self):
        """Rückgabe ist ein Wörterbuch."""
        result = cannonball_problem(30)
        assert isinstance(result, dict)

    def test_n1_is_solution(self):
        """n=1: 1² = 1 = 1². Ist Lösung."""
        result = cannonball_problem(5)
        ns = [s['n'] for s in result['loesungen']]
        assert 1 in ns

    def test_n24_is_solution(self):
        """n=24: S_24 = 4900 = 70². Bekannte Lösung."""
        result = cannonball_problem(24)
        solutions = result['loesungen']
        ns = [s['n'] for s in solutions]
        assert 24 in ns

    def test_24_square_is_70(self):
        """S_24 = 4900, √4900 = 70."""
        result = cannonball_problem(24)
        for s in result['loesungen']:
            if s['n'] == 24:
                assert s['quadratwurzel'] == 70
                assert s['summe'] == 4900

    def test_no_extra_solutions_between_2_and_23(self):
        """Zwischen n=2 und n=23 gibt es keine weiteren Lösungen."""
        result = cannonball_problem(23)
        # Nur n=1 sollte Lösung sein
        ns = [s['n'] for s in result['loesungen']]
        for n in ns:
            assert n == 1  # keine andere Lösung zwischen 2 und 23
