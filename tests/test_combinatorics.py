"""
@file test_combinatorics.py
@brief Umfassende Tests für combinatorics.py – Permutationsgruppen, abzählende Kombinatorik,
       erzeugende Funktionen, Ramsey-Theorie, Graphen-Kombinatorik, Kombinatorik auf Wörtern.
@author Kurt Ingwer
@lastModified 2026-03-10

Testabdeckung:
    - PermutationGroup: alle_permutationen, Vorzeichen, Zyklen, Ordnung, Komposition, Inverse
    - EnumerativeCombinatorics: Multinomial, Stirling 1./2. Art, Euler-Zahlen, Partition, Motzkin
    - GeneratingFunctions: OGF, EGF, Faltung, Rekurrenz
    - RamseyTheory: R(3,3)=6, Van-der-Waerden, Schubfachprinzip
    - GraphCombinatorics: chromatisches Polynom, Spannbäume, Matchings, Euler, Hamilton
    - CombinatoricsOnWords: Lyndon-Wörter, Halsketten, LCS, Edit-Distanz, De-Bruijn
    - Standalone: Fibonacci, Lucas, Tribonacci, Sylvester, Glückszahlen, Collatz-Baum
"""

import pytest
import math
from src.combinatorics import (
    PermutationGroup,
    EnumerativeCombinatorics,
    GeneratingFunctions,
    RamseyTheory,
    GraphCombinatorics,
    CombinatoricsOnWords,
    fibonacci,
    lucas_numbers,
    tribonacci,
    sylvester_sequence,
    happy_numbers,
    collatz_tree,
)


# ============================================================
# Tests: PermutationGroup
# ============================================================

class TestPermutationGroup:
    """Tests für die PermutationGroup-Klasse."""

    def test_init(self):
        """Symmetrische Gruppe S_3 initialisieren."""
        pg = PermutationGroup(3)
        assert pg.n == 3
        assert pg.identity == [0, 1, 2]

    def test_all_permutations_n0(self):
        """Leere Gruppe S_0 hat genau eine Permutation: []."""
        perms = PermutationGroup.all_permutations(0)
        assert perms == [[]]

    def test_all_permutations_n1(self):
        """S_1 hat genau 1 Permutation."""
        perms = PermutationGroup.all_permutations(1)
        assert perms == [[0]]

    def test_all_permutations_n3_count(self):
        """S_3 hat 3! = 6 Permutationen."""
        perms = PermutationGroup.all_permutations(3)
        assert len(perms) == 6

    def test_all_permutations_n4_count(self):
        """S_4 hat 4! = 24 Permutationen."""
        perms = PermutationGroup.all_permutations(4)
        assert len(perms) == 24

    def test_all_permutations_unique(self):
        """Alle Permutationen von S_3 sind verschieden."""
        perms = PermutationGroup.all_permutations(3)
        unique = set(tuple(p) for p in perms)
        assert len(unique) == 6

    def test_permutation_sign_identity(self):
        """Identitätspermutation hat Vorzeichen +1."""
        assert PermutationGroup.permutation_sign([0, 1, 2]) == 1

    def test_permutation_sign_single_transposition(self):
        """Einfache Transposition hat Vorzeichen -1."""
        # Tausche 0 und 1: [1, 0, 2] – eine Inversion
        assert PermutationGroup.permutation_sign([1, 0, 2]) == -1

    def test_permutation_sign_3_cycle(self):
        """3-Zyklus (0→1→2→0) ist gerade (Vorzeichen +1)."""
        # [1, 2, 0]: 0→1, 1→2, 2→0
        assert PermutationGroup.permutation_sign([1, 2, 0]) == 1

    def test_permutation_sign_empty(self):
        """Leere Permutation hat Vorzeichen +1."""
        assert PermutationGroup.permutation_sign([]) == 1

    def test_cycle_decomposition_identity(self):
        """Identität hat nur 1-elementige Zyklen (Fixpunkte)."""
        cycles = PermutationGroup.cycle_decomposition([0, 1, 2])
        # Alle Zyklen haben Länge 1
        assert all(len(c) == 1 for c in cycles)
        assert len(cycles) == 3

    def test_cycle_decomposition_transposition(self):
        """[1, 0, 2] hat einen 2-Zyklus und einen Fixpunkt."""
        cycles = PermutationGroup.cycle_decomposition([1, 0, 2])
        lengths = sorted(len(c) for c in cycles)
        assert lengths == [1, 2]

    def test_cycle_decomposition_3_cycle(self):
        """[1, 2, 0] ist ein 3-Zyklus."""
        cycles = PermutationGroup.cycle_decomposition([1, 2, 0])
        assert len(cycles) == 1
        assert len(cycles[0]) == 3

    def test_order_identity(self):
        """Ordnung der Identität ist 1."""
        assert PermutationGroup.order_of_permutation([0, 1, 2]) == 1

    def test_order_transposition(self):
        """Ordnung einer Transposition ist 2."""
        assert PermutationGroup.order_of_permutation([1, 0, 2]) == 2

    def test_order_3_cycle(self):
        """Ordnung eines 3-Zyklus ist 3."""
        assert PermutationGroup.order_of_permutation([1, 2, 0]) == 3

    def test_order_mixed(self):
        """Ordnung von [1, 0, 3, 2] (zwei Transpositionen) ist kgV(2,2)=2."""
        assert PermutationGroup.order_of_permutation([1, 0, 3, 2]) == 2

    def test_order_empty(self):
        """Leere Permutation hat Ordnung 1."""
        assert PermutationGroup.order_of_permutation([]) == 1

    def test_compose_permutations_identity(self):
        """Komposition mit Identität ergibt dieselbe Permutation."""
        p = [2, 0, 1]
        identity = [0, 1, 2]
        assert PermutationGroup.compose_permutations(p, identity) == p
        assert PermutationGroup.compose_permutations(identity, p) == p

    def test_compose_permutations_two_transpositions(self):
        """Komposition [1,0,2] ∘ [0,2,1] ergibt [1,2,0]."""
        p1 = [1, 0, 2]  # Tausch 0↔1
        p2 = [0, 2, 1]  # Tausch 1↔2
        result = PermutationGroup.compose_permutations(p1, p2)
        assert result == [1, 2, 0]

    def test_compose_permutations_error(self):
        """Komposition unterschiedlich langer Permutationen wirft ValueError."""
        with pytest.raises(ValueError):
            PermutationGroup.compose_permutations([0, 1], [0, 1, 2])

    def test_inverse_permutation_identity(self):
        """Inverse der Identität ist die Identität."""
        assert PermutationGroup.inverse_permutation([0, 1, 2]) == [0, 1, 2]

    def test_inverse_permutation_transposition(self):
        """Transposition ist selbstinvers."""
        p = [1, 0, 2]
        assert PermutationGroup.inverse_permutation(p) == p

    def test_inverse_permutation_3_cycle(self):
        """Inverse des 3-Zyklus [1,2,0] ist [2,0,1]."""
        assert PermutationGroup.inverse_permutation([1, 2, 0]) == [2, 0, 1]

    def test_inverse_inverse_is_original(self):
        """Doppelte Inverse ergibt die Originalpermutation."""
        p = [2, 0, 3, 1]
        inv = PermutationGroup.inverse_permutation(p)
        double_inv = PermutationGroup.inverse_permutation(inv)
        assert double_inv == p

    def test_is_even_identity(self):
        """Identität ist gerade."""
        assert PermutationGroup.is_even([0, 1, 2]) is True

    def test_is_even_transposition(self):
        """Einfache Transposition ist ungerade."""
        assert PermutationGroup.is_even([1, 0, 2]) is False

    def test_fixed_points_identity(self):
        """Identität hat alle Elemente als Fixpunkte."""
        assert PermutationGroup.fixed_points([0, 1, 2]) == [0, 1, 2]

    def test_fixed_points_transposition(self):
        """[1, 0, 2] hat Fixpunkt 2."""
        assert PermutationGroup.fixed_points([1, 0, 2]) == [2]

    def test_fixed_points_3_cycle(self):
        """3-Zyklus [1, 2, 0] hat keine Fixpunkte."""
        assert PermutationGroup.fixed_points([1, 2, 0]) == []

    def test_fixed_points_empty(self):
        """Leere Permutation hat keine Fixpunkte."""
        assert PermutationGroup.fixed_points([]) == []


# ============================================================
# Tests: EnumerativeCombinatorics
# ============================================================

class TestEnumerativeCombinatorics:
    """Tests für abzählende Kombinatorik."""

    def test_multinomial_basic(self):
        """Multinomialkoeffizient 4!/(2!·2!) = 6."""
        assert EnumerativeCombinatorics.multinomial(4, [2, 2]) == 6

    def test_multinomial_three_groups(self):
        """Multinomialkoeffizient 6!/(1!·2!·3!) = 60."""
        assert EnumerativeCombinatorics.multinomial(6, [1, 2, 3]) == 60

    def test_multinomial_single_group(self):
        """Multinomialkoeffizient n!/(n!) = 1."""
        assert EnumerativeCombinatorics.multinomial(5, [5]) == 1

    def test_multinomial_error(self):
        """Falsche Summe wirft ValueError."""
        with pytest.raises(ValueError):
            EnumerativeCombinatorics.multinomial(5, [2, 2])

    def test_stirling_first_base_cases(self):
        """Stirling 1. Art: s(0,0)=1, s(n,0)=0, s(0,k)=0."""
        assert EnumerativeCombinatorics.stirling_first(0, 0) == 1
        assert EnumerativeCombinatorics.stirling_first(3, 0) == 0
        assert EnumerativeCombinatorics.stirling_first(0, 2) == 0

    def test_stirling_first_values(self):
        """Bekannte Stirling-Zahlen 1. Art prüfen."""
        # s(3,1) = 2 (zwei 3-Zyklen: (123) und (132))
        assert EnumerativeCombinatorics.stirling_first(3, 1) == 2
        # s(3,2) = -3
        assert EnumerativeCombinatorics.stirling_first(3, 2) == -3
        # s(3,3) = 1
        assert EnumerativeCombinatorics.stirling_first(3, 3) == 1

    def test_stirling_first_row_sum(self):
        """Summe über k von |s(n,k)| = n! (Anzahl Permutationen in S_n)."""
        n = 4
        total = sum(abs(EnumerativeCombinatorics.stirling_first(n, k))
                    for k in range(n + 1))
        assert total == math.factorial(n)

    def test_stirling_second_base_cases(self):
        """Stirling 2. Art: S(0,0)=1, S(n,0)=0, S(0,k)=0."""
        assert EnumerativeCombinatorics.stirling_second(0, 0) == 1
        assert EnumerativeCombinatorics.stirling_second(3, 0) == 0
        assert EnumerativeCombinatorics.stirling_second(0, 2) == 0

    def test_stirling_second_values(self):
        """Bekannte Stirling-Zahlen 2. Art prüfen."""
        assert EnumerativeCombinatorics.stirling_second(3, 1) == 1
        assert EnumerativeCombinatorics.stirling_second(3, 2) == 3
        assert EnumerativeCombinatorics.stirling_second(3, 3) == 1
        assert EnumerativeCombinatorics.stirling_second(4, 2) == 7

    def test_stirling_second_bell_sum(self):
        """Summe S(n,k) über k = B_n (n-te Bell-Zahl)."""
        # B_4 = 15
        n = 4
        bell_4 = sum(EnumerativeCombinatorics.stirling_second(n, k) for k in range(n + 1))
        assert bell_4 == 15

    def test_euler_number_a_base(self):
        """Euler-Zahlen: A(0,0)=1, A(n,k)=0 außerhalb."""
        assert EnumerativeCombinatorics.euler_number_a(0, 0) == 1
        assert EnumerativeCombinatorics.euler_number_a(1, 0) == 1
        assert EnumerativeCombinatorics.euler_number_a(2, 0) == 1
        assert EnumerativeCombinatorics.euler_number_a(2, 1) == 1

    def test_euler_number_a_row_sum(self):
        """Summe A(n,k) über k = n! (alle Permutationen)."""
        n = 4
        total = sum(EnumerativeCombinatorics.euler_number_a(n, k) for k in range(n))
        assert total == math.factorial(n)

    def test_partition_number_small(self):
        """Partitionszahlen für kleine n."""
        expected = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        for i, e in enumerate(expected):
            assert EnumerativeCombinatorics.partition_number(i) == e

    def test_partition_number_zero(self):
        """p(0) = 1 (leere Partition)."""
        assert EnumerativeCombinatorics.partition_number(0) == 1

    def test_motzkin_number_small(self):
        """Motzkin-Zahlen: M_0=1, M_1=1, M_2=2, M_3=4, M_4=9."""
        assert EnumerativeCombinatorics.motzkin_number(0) == 1
        assert EnumerativeCombinatorics.motzkin_number(1) == 1
        assert EnumerativeCombinatorics.motzkin_number(2) == 2
        assert EnumerativeCombinatorics.motzkin_number(3) == 4
        assert EnumerativeCombinatorics.motzkin_number(4) == 9

    def test_narayana_number_catalan_sum(self):
        """Σ N(n,k) k=1..n = C_n (Catalan-Zahl)."""
        # C_5 = 42
        n = 5
        total = sum(EnumerativeCombinatorics.narayana_number(n, k) for k in range(1, n + 1))
        assert total == 42

    def test_narayana_number_basic(self):
        """N(3,1)=1, N(3,2)=3, N(3,3)=1."""
        assert EnumerativeCombinatorics.narayana_number(3, 1) == 1
        assert EnumerativeCombinatorics.narayana_number(3, 2) == 3
        assert EnumerativeCombinatorics.narayana_number(3, 3) == 1

    def test_ballot_problem_basic(self):
        """P = (3-2)/(3+2) = 0.2 für p=3, q=2."""
        result = EnumerativeCombinatorics.ballot_problem(3, 2)
        assert abs(result - 0.2) < 1e-12

    def test_ballot_problem_tie(self):
        """P = 0 wenn p ≤ q."""
        assert EnumerativeCombinatorics.ballot_problem(2, 3) == 0.0
        assert EnumerativeCombinatorics.ballot_problem(3, 3) == 0.0

    def test_ballot_problem_large_margin(self):
        """P = 1 wenn q = 0 (A erhält alle Stimmen)."""
        result = EnumerativeCombinatorics.ballot_problem(10, 0)
        assert abs(result - 1.0) < 1e-12

    def test_inclusion_exclusion_two_sets(self):
        """|A ∪ B| = |A| + |B| - |A ∩ B|."""
        A = {1, 2, 3, 4}
        B = {3, 4, 5, 6}
        result = EnumerativeCombinatorics.inclusion_exclusion([A, B])
        assert result == len(A | B)

    def test_inclusion_exclusion_three_sets(self):
        """Inklusion-Exklusion für drei Mengen."""
        A = {1, 2, 3}
        B = {2, 3, 4}
        C = {3, 4, 5}
        result = EnumerativeCombinatorics.inclusion_exclusion([A, B, C])
        assert result == len(A | B | C)

    def test_inclusion_exclusion_empty(self):
        """Leere Mengenliste: |∅| = 0."""
        assert EnumerativeCombinatorics.inclusion_exclusion([]) == 0

    def test_inclusion_exclusion_disjoint(self):
        """Disjunkte Mengen: |A ∪ B| = |A| + |B|."""
        A = {1, 2}
        B = {3, 4}
        assert EnumerativeCombinatorics.inclusion_exclusion([A, B]) == 4


# ============================================================
# Tests: GeneratingFunctions
# ============================================================

class TestGeneratingFunctions:
    """Tests für erzeugende Funktionen."""

    def test_ordinary_gf_basic(self):
        """OGF gibt die Koeffizienten unverändert zurück."""
        seq = [1, 2, 3, 4, 5]
        result = GeneratingFunctions.ordinary_gf(seq, 5)
        assert result == [1, 2, 3, 4, 5]

    def test_ordinary_gf_padding(self):
        """OGF füllt mit Nullen auf wenn Folge zu kurz."""
        result = GeneratingFunctions.ordinary_gf([1, 2], 5)
        assert result == [1, 2, 0, 0, 0]

    def test_ordinary_gf_truncate(self):
        """OGF kürzt wenn Folge zu lang."""
        result = GeneratingFunctions.ordinary_gf([1, 2, 3, 4, 5], 3)
        assert result == [1, 2, 3]

    def test_exponential_gf_basic(self):
        """EGF: a_n/n! für n=0,1,2."""
        # Seq = [1, 1, 2]: EGF = [1/0!, 1/1!, 2/2!] = [1, 1, 1]
        result = GeneratingFunctions.exponential_gf([1, 1, 2], 3)
        assert abs(result[0] - 1.0) < 1e-12
        assert abs(result[1] - 1.0) < 1e-12
        assert abs(result[2] - 1.0) < 1e-12

    def test_fibonacci_gf(self):
        """Fibonacci OGF-String korrekt."""
        gf = GeneratingFunctions.fibonacci_gf()
        assert "1 - x - x^2" in gf

    def test_catalan_gf(self):
        """Catalan OGF-String korrekt."""
        gf = GeneratingFunctions.catalan_gf()
        assert "sqrt" in gf or "√" in gf or "1 - 4x" in gf

    def test_power_series_multiply_identity(self):
        """Multiplikation mit [1] ergibt die ursprüngliche Reihe."""
        f = [1.0, 2.0, 3.0]
        g = [1.0]  # Konstante 1
        result = GeneratingFunctions.power_series_multiply(f, g, 3)
        assert result == pytest.approx([1.0, 2.0, 3.0])

    def test_power_series_multiply_convolution(self):
        """Faltung: (1+x)^2 = 1 + 2x + x^2."""
        f = [1.0, 1.0]  # 1 + x
        g = [1.0, 1.0]  # 1 + x
        result = GeneratingFunctions.power_series_multiply(f, g, 3)
        assert result == pytest.approx([1.0, 2.0, 1.0])

    def test_gf_from_recurrence_fibonacci(self):
        """Fibonacci-Rekurrenz: a_n = a_{n-1} + a_{n-2}."""
        fib_rec = lambda seq, n: seq[n - 1] + seq[n - 2]
        result = GeneratingFunctions.gf_from_recurrence(0, 1, fib_rec, 8)
        assert result == [0, 1, 1, 2, 3, 5, 8, 13]


# ============================================================
# Tests: RamseyTheory
# ============================================================

class TestRamseyTheory:
    """Tests für Ramsey-Theorie."""

    def test_ramsey_R33(self):
        """R(3,3) = 6 (bekannteste Ramsey-Zahl)."""
        assert RamseyTheory.ramsey_number_R(3, 3) == 6

    def test_ramsey_R34(self):
        """R(3,4) = 9."""
        assert RamseyTheory.ramsey_number_R(3, 4) == 9

    def test_ramsey_R44(self):
        """R(4,4) = 18."""
        assert RamseyTheory.ramsey_number_R(4, 4) == 18

    def test_ramsey_symmetry(self):
        """R(s,t) = R(t,s)."""
        assert RamseyTheory.ramsey_number_R(3, 4) == RamseyTheory.ramsey_number_R(4, 3)

    def test_ramsey_R11(self):
        """R(1,1) = 1."""
        assert RamseyTheory.ramsey_number_R(1, 1) == 1

    def test_ramsey_unknown(self):
        """Unbekannte Ramsey-Zahlen geben None zurück."""
        result = RamseyTheory.ramsey_number_R(6, 6)
        assert result is None

    def test_ramsey_graph_coloring_edges(self):
        """K_4 hat C(4,2) = 6 Kanten."""
        coloring = RamseyTheory.ramsey_graph_coloring(4, 2)
        assert len(coloring) == 6

    def test_ramsey_graph_coloring_valid_colors(self):
        """Alle Farben liegen in {0, ..., k-1}."""
        coloring = RamseyTheory.ramsey_graph_coloring(5, 3)
        assert all(0 <= c < 3 for c in coloring.values())

    def test_van_der_waerden_numbers(self):
        """W(3;2) = 9 ist ein bekannter Wert."""
        numbers = RamseyTheory.van_der_waerden_numbers()
        assert numbers[(3, 2)] == 9
        assert numbers[(4, 2)] == 35

    def test_hales_jewett_demo(self):
        """Hales-Jewett Demo gibt Dictionary zurück."""
        result = RamseyTheory.hales_jewett_theorem_demo(3, 2)
        assert result["dimension"] == 3
        assert result["alphabet_size"] == 2
        assert result["cube_size"] == 8

    def test_pigeonhole_principle_basic(self):
        """5 Tauben in 4 Löcher: mindestens 2 in einem Loch."""
        result = RamseyTheory.pigeonhole_principle(5, 4)
        assert result["at_least_two_in_one_hole"] is True
        assert result["guaranteed_min_in_one_hole"] == 2

    def test_pigeonhole_principle_equal(self):
        """Gleich viele Tauben wie Löcher: kein Loch muss mehr als 1 haben."""
        result = RamseyTheory.pigeonhole_principle(4, 4)
        assert result["at_least_two_in_one_hole"] is False

    def test_pigeonhole_principle_zero_holes_error(self):
        """0 Löcher wirft ValueError."""
        with pytest.raises(ValueError):
            RamseyTheory.pigeonhole_principle(5, 0)


# ============================================================
# Tests: GraphCombinatorics
# ============================================================

class TestGraphCombinatorics:
    """Tests für Graphen-Kombinatorik."""

    def _k3_adj(self):
        """Adjazenzmatrix von K_3 (vollständiger Graph, 3 Knoten)."""
        return [[0, 1, 1], [1, 0, 1], [1, 1, 0]]

    def _k4_adj(self):
        """Adjazenzmatrix von K_4 (vollständiger Graph, 4 Knoten)."""
        n = 4
        return [[0 if i == j else 1 for j in range(n)] for i in range(n)]

    def _path_adj(self, n):
        """Adjazenzmatrix eines Pfadgraphen P_n."""
        adj = [[0] * n for _ in range(n)]
        for i in range(n - 1):
            adj[i][i + 1] = 1
            adj[i + 1][i] = 1
        return adj

    def _cycle_adj(self, n):
        """Adjazenzmatrix eines Kreisgraphen C_n."""
        adj = self._path_adj(n)
        adj[0][n - 1] = 1
        adj[n - 1][0] = 1
        return adj

    def test_chromatic_polynomial_empty_graph(self):
        """Leerer Graph (kein Knoten): P = [1]."""
        result = GraphCombinatorics.chromatic_polynomial([])
        assert result == [1]

    def test_chromatic_polynomial_single_node(self):
        """Ein Knoten ohne Kanten: P(k) = k, also [0, 1]."""
        result = GraphCombinatorics.chromatic_polynomial([[0]])
        assert result == [0, 1]

    def test_chromatic_polynomial_k3(self):
        """K_3: P(k) = k(k-1)(k-2) = k³ - 3k² + 2k → [-0, 2, -3, 1]."""
        result = GraphCombinatorics.chromatic_polynomial(self._k3_adj())
        # Auswerten bei k=3: 3! = 6
        k = 3
        value = sum(result[i] * (k ** i) for i in range(len(result)))
        assert value == 6  # 3 Farben für K_3: 3·2·1 = 6 Färbungen

    def test_chromatic_polynomial_k3_k2(self):
        """K_3 kann mit 2 Farben nicht gefärbt werden: P(2) = 0."""
        result = GraphCombinatorics.chromatic_polynomial(self._k3_adj())
        k = 2
        value = sum(result[i] * (k ** i) for i in range(len(result)))
        assert value == 0

    def test_chromatic_polynomial_path_p3(self):
        """P_3 (Pfad 0-1-2): P(k) = k(k-1)^2."""
        adj = self._path_adj(3)
        result = GraphCombinatorics.chromatic_polynomial(adj)
        # P(3) = 3·4 = 12
        k = 3
        value = sum(result[i] * (k ** i) for i in range(len(result)))
        assert value == 12

    def test_spanning_tree_count_k3(self):
        """K_3 hat 3 Spannbäume (Cayley-Formel: n^{n-2} = 3^1 = 3)."""
        result = GraphCombinatorics.spanning_tree_count(self._k3_adj())
        assert result == 3

    def test_spanning_tree_count_k4(self):
        """K_4 hat 4^2 = 16 Spannbäume (Cayley-Formel)."""
        result = GraphCombinatorics.spanning_tree_count(self._k4_adj())
        assert result == 16

    def test_spanning_tree_count_path(self):
        """Pfadgraph hat genau 1 Spannbaum (er ist selbst ein Baum)."""
        result = GraphCombinatorics.spanning_tree_count(self._path_adj(4))
        assert result == 1

    def test_spanning_tree_count_single(self):
        """Einzelknoten hat 1 Spannbaum."""
        assert GraphCombinatorics.spanning_tree_count([[0]]) == 1

    def test_tutte_polynomial_demo(self):
        """Tutte-Polynom Demo gibt korrekte Spannbaumanzahl."""
        result = GraphCombinatorics.tutte_polynomial_demo(self._k3_adj())
        assert result["T_1_1"] == 3  # K_3 hat 3 Spannbäume

    def test_matching_number_k3(self):
        """K_3 hat Matchingzahl 1 (kein perfektes Matching bei ungerader Knotenanzahl)."""
        result = GraphCombinatorics.matching_number(self._k3_adj())
        assert result == 1

    def test_matching_number_path_p4(self):
        """P_4: maximales Matching hat Größe 2."""
        result = GraphCombinatorics.matching_number(self._path_adj(4))
        assert result == 2

    def test_matching_number_empty(self):
        """Leerer Graph: Matchingzahl 0."""
        assert GraphCombinatorics.matching_number([]) == 0

    def test_perfect_matching_count_k4(self):
        """K_4 hat 3 perfekte Matchings."""
        result = GraphCombinatorics.perfect_matching_count(self._k4_adj())
        assert result == 3

    def test_perfect_matching_count_odd(self):
        """Kein perfektes Matching bei ungerader Knotenanzahl."""
        result = GraphCombinatorics.perfect_matching_count(self._k3_adj())
        assert result == 0

    def test_eulerian_circuit_k3(self):
        """K_3 hat einen Euler-Kreis (alle Grad 2, gerade)."""
        assert GraphCombinatorics.eulerian_circuit_exists(self._k3_adj()) is True

    def test_eulerian_circuit_path(self):
        """Pfadgraph P_3 hat keinen Euler-Kreis (Endknoten haben Grad 1)."""
        assert GraphCombinatorics.eulerian_circuit_exists(self._path_adj(3)) is False

    def test_eulerian_circuit_empty(self):
        """Leerer Graph: Euler-Kreis trivialerweise vorhanden."""
        assert GraphCombinatorics.eulerian_circuit_exists([]) is True

    def test_hamiltonian_path_k3(self):
        """K_3 hat einen Hamiltonpfad."""
        assert GraphCombinatorics.hamiltonian_path_check(self._k3_adj()) is True

    def test_hamiltonian_path_single_node(self):
        """Einzelknoten hat einen Hamiltonpfad (er selbst)."""
        assert GraphCombinatorics.hamiltonian_path_check([[0]]) is True

    def test_hamiltonian_path_empty(self):
        """Leerer Graph: trivial True."""
        assert GraphCombinatorics.hamiltonian_path_check([]) is True

    def test_hamiltonian_path_disconnected(self):
        """Nicht zusammenhängender Graph hat keinen Hamiltonpfad."""
        # Graph: 4 Knoten, keine Kanten → kein Hamiltonpfad
        adj = [[0, 0, 0, 0]] * 4
        assert GraphCombinatorics.hamiltonian_path_check(adj) is False


# ============================================================
# Tests: CombinatoricsOnWords
# ============================================================

class TestCombinatoricsOnWords:
    """Tests für Kombinatorik auf Wörtern."""

    def test_lyndon_word_basic(self):
        """'ab' ist ein Lyndon-Wort (a < b, keine kleinere Rotation)."""
        assert CombinatoricsOnWords.lyndon_word_check("ab") is True

    def test_lyndon_word_abcd(self):
        """'abcd' ist ein Lyndon-Wort."""
        assert CombinatoricsOnWords.lyndon_word_check("abcd") is True

    def test_lyndon_word_not_aab(self):
        """'aab' ist kein Lyndon-Wort: Rotation 'aba' > 'aab', aber 'baa' > 'aab'...
        Nochmals: 'ab' + 'a' → Rotation 1: 'aba', Rotation 2: 'baa'. 'aab' < 'aba' ✓,
        'aab' < 'baa' ✓ → doch Lyndon-Wort? Nein: 'aab' ist Rot0, Rot1='aba', Rot2='baa'.
        'aab' < alle → tatsächlich Lyndon."""
        # Korrekte Prüfung: 'aab' → Rot1='aba' > 'aab' ✓, Rot2='baa' > 'aab' ✓ → Lyndon
        assert CombinatoricsOnWords.lyndon_word_check("aab") is True

    def test_lyndon_word_not_aba(self):
        """'aba' ist kein Lyndon-Wort: Rotation 'aab' < 'aba'."""
        assert CombinatoricsOnWords.lyndon_word_check("aba") is False

    def test_lyndon_word_empty(self):
        """Leeres Wort ist kein Lyndon-Wort."""
        assert CombinatoricsOnWords.lyndon_word_check("") is False

    def test_lyndon_word_single(self):
        """Einzelzeichen ist ein Lyndon-Wort (keine echten Rotationen)."""
        assert CombinatoricsOnWords.lyndon_word_check("a") is True

    def test_lyndon_word_repeated(self):
        """'aa' ist kein Lyndon-Wort (Rotation ergibt gleiches Wort)."""
        assert CombinatoricsOnWords.lyndon_word_check("aa") is False

    def test_necklaces_k2_n3(self):
        """3 Perlen, 2 Farben: 4 Halsketten (000,001,011,111)."""
        assert CombinatoricsOnWords.necklaces(3, 2) == 4

    def test_necklaces_k2_n4(self):
        """4 Perlen, 2 Farben: 6 Halsketten."""
        assert CombinatoricsOnWords.necklaces(4, 2) == 6

    def test_necklaces_k1_any(self):
        """1 Farbe: immer genau 1 Halskette."""
        assert CombinatoricsOnWords.necklaces(5, 1) == 1

    def test_necklaces_n0(self):
        """0 Perlen: 1 Halskette."""
        assert CombinatoricsOnWords.necklaces(0, 2) == 1

    def test_primitive_word_basic(self):
        """'ab' ist primitiv."""
        assert CombinatoricsOnWords.primitive_word("ab") is True

    def test_primitive_word_repeated(self):
        """'abab' ist nicht primitiv ('ab' ^ 2)."""
        assert CombinatoricsOnWords.primitive_word("abab") is False

    def test_primitive_word_single(self):
        """Einzelzeichen ist primitiv."""
        assert CombinatoricsOnWords.primitive_word("a") is True

    def test_lcs_basic(self):
        """LCS von 'ABCBDAB' und 'BDCAB' hat Länge 4."""
        result = CombinatoricsOnWords.longest_common_subsequence("ABCBDAB", "BDCAB")
        assert result == 4

    def test_lcs_identical(self):
        """LCS identischer Strings ist die Länge selbst."""
        s = "hello"
        assert CombinatoricsOnWords.longest_common_subsequence(s, s) == len(s)

    def test_lcs_empty(self):
        """LCS mit leerem String ist 0."""
        assert CombinatoricsOnWords.longest_common_subsequence("abc", "") == 0
        assert CombinatoricsOnWords.longest_common_subsequence("", "abc") == 0

    def test_lcs_no_common(self):
        """LCS ohne gemeinsame Zeichen ist 0."""
        assert CombinatoricsOnWords.longest_common_subsequence("abc", "xyz") == 0

    def test_edit_distance_basic(self):
        """Edit-Distanz von 'kitten' und 'sitting' ist 3."""
        result = CombinatoricsOnWords.edit_distance("kitten", "sitting")
        assert result == 3

    def test_edit_distance_identical(self):
        """Edit-Distanz identischer Strings ist 0."""
        assert CombinatoricsOnWords.edit_distance("hello", "hello") == 0

    def test_edit_distance_empty(self):
        """Edit-Distanz zu leerem String ist Länge des anderen."""
        assert CombinatoricsOnWords.edit_distance("abc", "") == 3
        assert CombinatoricsOnWords.edit_distance("", "abc") == 3

    def test_edit_distance_single_insert(self):
        """Edit-Distanz durch Einfügen: 'ab' → 'abc' = 1."""
        assert CombinatoricsOnWords.edit_distance("ab", "abc") == 1

    def test_edit_distance_single_delete(self):
        """Edit-Distanz durch Löschen: 'abc' → 'ab' = 1."""
        assert CombinatoricsOnWords.edit_distance("abc", "ab") == 1

    def test_edit_distance_symmetric(self):
        """Edit-Distanz ist symmetrisch."""
        s1, s2 = "algorithm", "altruistic"
        assert CombinatoricsOnWords.edit_distance(s1, s2) == CombinatoricsOnWords.edit_distance(s2, s1)

    def test_de_bruijn_n1_k2(self):
        """B(2,1): enthält alle 1-mere aus {0,1}."""
        seq = CombinatoricsOnWords.de_bruijn_sequence(1, 2)
        # Länge 2^1 = 2, enthält '0' und '1'
        assert len(seq) == 2
        assert '0' in seq and '1' in seq

    def test_de_bruijn_n2_k2(self):
        """B(2,2): enthält alle 2-mere aus {0,1}² zyklisch."""
        seq = CombinatoricsOnWords.de_bruijn_sequence(2, 2)
        # Länge 2^2 = 4
        assert len(seq) == 4
        # Alle 2-mere {00, 01, 10, 11} müssen zyklisch vorkommen
        cyclic = seq + seq
        kmers = {cyclic[i:i+2] for i in range(len(seq))}
        assert {"00", "01", "10", "11"} == kmers

    def test_de_bruijn_n3_k2_length(self):
        """B(2,3): Länge 2^3 = 8."""
        seq = CombinatoricsOnWords.de_bruijn_sequence(3, 2)
        assert len(seq) == 8


# ============================================================
# Tests: Standalone-Funktionen
# ============================================================

class TestStandaloneFunctions:
    """Tests für standalone kombinatorische Funktionen."""

    def test_fibonacci_base_cases(self):
        """F_0=0, F_1=1."""
        assert fibonacci(0) == 0
        assert fibonacci(1) == 1

    def test_fibonacci_small(self):
        """Fibonacci-Zahlen für kleine n."""
        expected = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        for i, e in enumerate(expected):
            assert fibonacci(i) == e

    def test_fibonacci_large(self):
        """Fibonacci(20) = 6765."""
        assert fibonacci(20) == 6765

    def test_fibonacci_negative_error(self):
        """Negative Eingabe wirft ValueError."""
        with pytest.raises(ValueError):
            fibonacci(-1)

    def test_lucas_numbers_basic(self):
        """Lucas-Zahlen: L_0=2, L_1=1, L_2=3, L_3=4, L_4=7."""
        result = lucas_numbers(5)
        assert result == [2, 1, 3, 4, 7]

    def test_lucas_numbers_empty(self):
        """lucas_numbers(0) = []."""
        assert lucas_numbers(0) == []

    def test_lucas_numbers_single(self):
        """lucas_numbers(1) = [2]."""
        assert lucas_numbers(1) == [2]

    def test_tribonacci_basic(self):
        """Tribonacci: T_0=0, T_1=0, T_2=1, T_3=1, T_4=2, T_5=4."""
        result = tribonacci(6)
        assert result == [0, 0, 1, 1, 2, 4]

    def test_tribonacci_empty(self):
        """tribonacci(0) = []."""
        assert tribonacci(0) == []

    def test_tribonacci_recurrence(self):
        """Prüfe Rekurrenz: T_n = T_{n-1} + T_{n-2} + T_{n-3}."""
        result = tribonacci(10)
        for i in range(3, 10):
            assert result[i] == result[i-1] + result[i-2] + result[i-3]

    def test_sylvester_sequence_basic(self):
        """Sylvester: a_0=2, a_1=3, a_2=7, a_3=43."""
        result = sylvester_sequence(4)
        assert result == [2, 3, 7, 43]

    def test_sylvester_sequence_empty(self):
        """sylvester_sequence(0) = []."""
        assert sylvester_sequence(0) == []

    def test_sylvester_sequence_property(self):
        """Ägyptische Brüche: 1 = Σ 1/(a_0 · a_1 · ... · a_{k-1}) · (a_{k}-1)/a_k."""
        # Einfachere Prüfung: Rekurrenz a_n = a_{n-1}*(a_{n-1}-1)+1
        result = sylvester_sequence(5)
        for i in range(1, 5):
            expected = result[i-1] * (result[i-1] - 1) + 1
            assert result[i] == expected

    def test_happy_numbers_basic(self):
        """Glückszahlen bis 20: 1, 7, 10, 13, 19."""
        result = happy_numbers(20)
        assert 1 in result
        assert 7 in result
        assert 10 in result
        assert 13 in result
        assert 19 in result
        # Unglückliche Zahlen nicht enthalten
        assert 2 not in result
        assert 3 not in result
        assert 4 not in result

    def test_happy_numbers_empty(self):
        """Keine Glückszahlen bis 0."""
        assert happy_numbers(0) == []

    def test_happy_numbers_singleton(self):
        """1 ist die erste Glückszahl."""
        assert happy_numbers(1) == [1]

    def test_collatz_tree_basic(self):
        """Collatz-Baum: 1 ist Vorgänger von 2 (2→1), und 2 → 1."""
        tree = collatz_tree(10)
        # 2 führt zu 1 (2/2=1), also 2 ist Vorgänger von 1
        assert 2 in tree[1]

    def test_collatz_tree_structure(self):
        """Collatz-Baum enthält alle Zahlen von 1 bis n_max als Schlüssel."""
        tree = collatz_tree(15)
        for i in range(1, 16):
            assert i in tree

    def test_collatz_tree_4_successor(self):
        """4 → 2 (4/2=2): 4 ist Vorgänger von 2."""
        tree = collatz_tree(10)
        assert 4 in tree[2]

    def test_collatz_tree_odd(self):
        """5 → 16 (3*5+1=16): 5 sollte Vorgänger von 16 sein, aber 16 > n_max=10."""
        tree = collatz_tree(10)
        # 16 ist nicht im Baum, also hat 5 keine nachverfolgten Vorgänger-Beziehungen
        # zu Zahlen > 10
        assert 16 not in tree
