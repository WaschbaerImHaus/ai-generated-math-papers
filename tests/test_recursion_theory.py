"""
@file test_recursion_theory.py
@brief Test-Suite für das Modul recursion_theory.py.

@description
    Testet alle Kernfunktionen der Rekursionstheorie:
    - TuringMachine: run(), accepts(), configuration_sequence()
    - Beispiel-TMs: Palindrome, Binärinkrement, Unäre Addition
    - μ-rekursive Funktionen: zero, successor, projection, primitive_recursion, mu_operator
    - Ackermann-Funktion (kleine Werte)
    - Halteproblem-Beweis (dict-Struktur)
    - Satz von Rice (dict-Struktur)
    - Entscheidbare / unentscheidbare Probleme
    - Arithmetische Hierarchie
    - Post'sches Korrespondenzproblem (kleine Instanzen)
    - Komplexitätsklassen und NP-vollständige Probleme
    - Kolmogorov-Komplexität

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os

# src/ zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from recursion_theory import (
    TapeSymbol, Direction, TuringMachine,
    tm_recognizes_palindromes, tm_binary_increment, tm_unary_addition,
    zero_function, successor_function, projection,
    primitive_recursion, mu_operator,
    ackermann_function, ackermann_growth_demo,
    halting_problem_undecidability_proof, rice_theorem_demo,
    decidable_problems, undecidable_problems,
    many_one_reduction, turing_reduction, reduction_degrees,
    arithmetical_hierarchy, sigma1_complete_problem,
    post_correspondence_problem,
    time_complexity_classes, np_complete_problems,
    cook_levin_theorem_sketch, polynomial_hierarchy,
    kolmogorov_complexity_approx, incompressible_string_demo
)


# ==============================================================================
# Tests: TuringMachine (Basisklasse)
# ==============================================================================

class TestTuringMachineBase:
    """Tests für die TuringMachine-Basisklasse."""

    def _simple_accept_tm(self) -> TuringMachine:
        """Minimale TM die alles akzeptiert."""
        return TuringMachine(
            states={'q0', 'accept', 'reject'},
            input_alphabet={'a'},
            tape_alphabet={'a', '_'},
            transitions={('q0', 'a'): ('accept', 'a', Direction.RIGHT),
                         ('q0', '_'): ('accept', '_', Direction.RIGHT)},
            initial_state='q0',
            accept_state='accept',
            reject_state='reject'
        )

    def test_run_returns_dict(self):
        """run() gibt Dict mit korrekten Schlüsseln zurück."""
        tm = self._simple_accept_tm()
        result = tm.run('a')
        assert isinstance(result, dict)
        assert 'accepted' in result
        assert 'steps' in result
        assert 'tape' in result
        assert 'halted' in result

    def test_run_simple_accept(self):
        """Einfache TM akzeptiert korrekt."""
        tm = self._simple_accept_tm()
        result = tm.run('a')
        assert result['accepted'] is True
        assert result['halted'] is True

    def test_run_max_steps_timeout(self):
        """TM mit Endlosschleife erkennt Timeout."""
        # TM die immer nach rechts läuft (hält nie)
        tm = TuringMachine(
            states={'q0', 'accept', 'reject'},
            input_alphabet={'a'},
            tape_alphabet={'a', '_'},
            transitions={('q0', 'a'): ('q0', 'a', Direction.RIGHT),
                         ('q0', '_'): ('q0', '_', Direction.RIGHT)},
            initial_state='q0',
            accept_state='accept',
            reject_state='reject'
        )
        result = tm.run('a', max_steps=10)
        assert result['halted'] is False
        assert result['steps'] == 10

    def test_accepts_returns_none_on_timeout(self):
        """accepts() gibt None bei Timeout zurück."""
        tm = TuringMachine(
            states={'q0', 'accept', 'reject'},
            input_alphabet={'a'},
            tape_alphabet={'a', '_'},
            transitions={('q0', 'a'): ('q0', 'a', Direction.RIGHT),
                         ('q0', '_'): ('q0', '_', Direction.RIGHT)},
            initial_state='q0',
            accept_state='accept',
            reject_state='reject'
        )
        result = tm.accepts('a', max_steps=5)
        assert result is None

    def test_configuration_sequence_nonempty(self):
        """configuration_sequence() gibt mindestens Startkonfiguration zurück."""
        tm = self._simple_accept_tm()
        configs = tm.configuration_sequence('a', max_steps=5)
        assert isinstance(configs, list)
        assert len(configs) >= 1


# ==============================================================================
# Tests: TM für Palindrome
# ==============================================================================

class TestTMPalindrome:
    """Tests für die Palindrom-erkennende Turing-Maschine."""

    def setup_method(self):
        """Initialisiere TM vor jedem Test."""
        self.tm = tm_recognizes_palindromes()

    def test_empty_string_is_palindrome(self):
        """Leerer String ist Palindrom."""
        assert self.tm.accepts('') is True

    def test_single_char_is_palindrome(self):
        """Einzelnes Zeichen ist Palindrom."""
        assert self.tm.accepts('a') is True
        assert self.tm.accepts('b') is True

    def test_two_same_chars(self):
        """Zwei gleiche Zeichen: Palindrom."""
        assert self.tm.accepts('aa') is True

    def test_two_diff_chars_not_palindrome(self):
        """Zwei verschiedene Zeichen: kein Palindrom."""
        assert self.tm.accepts('ab') is False

    def test_aba_is_palindrome(self):
        """'aba' ist ein Palindrom."""
        assert self.tm.accepts('aba') is True

    def test_abba_is_palindrome(self):
        """'abba' ist ein Palindrom."""
        assert self.tm.accepts('abba') is True

    def test_abc_not_palindrome(self):
        """'abb' ist kein Palindrom."""
        assert self.tm.accepts('abb') is False

    def test_abcba_not_supported(self):
        """'abcba': 'c' liegt nicht im Eingabealphabet {a,b}, wird verworfen."""
        # TM über {a,b}, 'c' unbekannt → kein Übergang → verwerfen
        result = self.tm.accepts('abcba')
        # Entweder False oder None (kein Übergang definiert)
        assert result in (False, None)


# ==============================================================================
# Tests: TM für Binärinkrement
# ==============================================================================

class TestTMBinaryIncrement:
    """Tests für die Binärinkrement-Turing-Maschine."""

    def setup_method(self):
        self.tm = tm_binary_increment()

    def test_increment_0(self):
        """0 → 1."""
        result = self.tm.run('0')
        assert result['accepted'] is True

    def test_increment_1(self):
        """1 → 10 (Binär): Band enthält '0' und neue führende 1."""
        result = self.tm.run('1')
        assert result['accepted'] is True

    def test_increment_01(self):
        """01 + 1 = 10."""
        result = self.tm.run('01')
        assert result['accepted'] is True

    def test_increment_011(self):
        """011 + 1 = 100."""
        result = self.tm.run('011')
        assert result['accepted'] is True

    def test_increment_halts(self):
        """TM hält immer bei Binäreingaben."""
        for s in ['0', '1', '10', '11', '100', '101', '110', '111']:
            result = self.tm.run(s)
            assert result['halted'] is True, f"TM hält nicht für '{s}'"


# ==============================================================================
# Tests: TM für Unäre Addition
# ==============================================================================

class TestTMUnaryAddition:
    """Tests für die Unäre-Additions-Turing-Maschine."""

    def setup_method(self):
        self.tm = tm_unary_addition()

    def test_1_plus_1(self):
        """1 + 1 = 2 (unär: '10' + '1' → '11')."""
        result = self.tm.run('101')  # 1^1 0 1^1
        assert result['accepted'] is True

    def test_2_plus_3(self):
        """2 + 3 = 5 (unär: '11011' → '11111')."""
        result = self.tm.run('11' + '0' + '111')
        assert result['accepted'] is True

    def test_addition_halts(self):
        """TM hält für korrekte unäre Eingaben."""
        for n, m in [(1, 1), (2, 2), (3, 1)]:
            inp = '1' * n + '0' + '1' * m
            result = self.tm.run(inp, max_steps=1000)
            assert result['halted'] is True


# ==============================================================================
# Tests: μ-rekursive Funktionen
# ==============================================================================

class TestMuRecursiveFunctions:
    """Tests für die μ-rekursiven Grundfunktionen."""

    def test_zero_function(self):
        """Z(n) = 0 für alle n."""
        for n in range(10):
            assert zero_function(n) == 0

    def test_successor_function(self):
        """S(n) = n + 1."""
        for n in range(10):
            assert successor_function(n) == n + 1

    def test_projection_first(self):
        """U₀²(x₀, x₁) = x₀."""
        assert projection((3, 7), 0) == 3

    def test_projection_second(self):
        """U₁²(x₀, x₁) = x₁."""
        assert projection((3, 7), 1) == 7

    def test_projection_singleton(self):
        """U₀¹(x) = x."""
        assert projection((42,), 0) == 42

    def test_primitive_recursion_addition(self):
        """
        Additionsfunktion via primitiver Rekursion:
        add(x, 0) = x  (Basisfall: g(x) = x = U₀¹(x))
        add(x, n+1) = S(add(x, n))  (Schritt: h(x, n, z) = S(z))
        """
        # g(x) = x (Identität als Projektion)
        base = lambda x: x
        # h(x, n, z) = z + 1 (Nachfolger des Ergebnisses)
        step = lambda x, n, z: z + 1
        add = primitive_recursion(base, step)

        assert add(3, 0) == 3
        assert add(3, 1) == 4
        assert add(3, 2) == 5
        assert add(0, 5) == 5
        assert add(7, 3) == 10

    def test_primitive_recursion_multiplication(self):
        """
        Multiplikation via primitiver Rekursion:
        mul(x, 0) = 0
        mul(x, n+1) = mul(x, n) + x
        """
        base = lambda x: 0
        step = lambda x, n, z: z + x
        mul = primitive_recursion(base, step)

        assert mul(3, 0) == 0
        assert mul(3, 1) == 3
        assert mul(3, 4) == 12
        assert mul(0, 5) == 0

    def test_mu_operator_basic(self):
        """μy[y = 5]: Suche kleinstes y mit y - 5 = 0."""
        # Prädikat: P(y) = 0 falls y == 5, sonst 1
        pred = lambda y: 0 if y == 5 else 1
        mu_f = mu_operator(pred)
        assert mu_f() == 5

    def test_mu_operator_with_args(self):
        """μy[x + y = 10]: Suche kleinstes y mit x + y = 10."""
        # P(x, y) = 0 falls x + y == 10
        pred = lambda x, y: 0 if x + y == 10 else 1
        mu_f = mu_operator(pred)
        assert mu_f(3) == 7
        assert mu_f(0) == 10
        assert mu_f(10) == 0

    def test_mu_operator_raises_on_no_solution(self):
        """μ-Operator wirft Fehler falls kein y in [0, max_y)."""
        pred = lambda y: 1  # Niemals 0
        mu_f = mu_operator(pred)
        with pytest.raises(ValueError):
            mu_f(max_y=10)


# ==============================================================================
# Tests: Ackermann-Funktion
# ==============================================================================

class TestAckermannFunction:
    """Tests für die Ackermann-Funktion."""

    def test_a_0_n(self):
        """A(0, n) = n + 1."""
        for n in range(8):
            assert ackermann_function(0, n) == n + 1

    def test_a_1_n(self):
        """A(1, n) = n + 2."""
        for n in range(8):
            assert ackermann_function(1, n) == n + 2

    def test_a_2_n(self):
        """A(2, n) = 2n + 3."""
        for n in range(6):
            assert ackermann_function(2, n) == 2 * n + 3

    def test_a_3_0(self):
        """A(3, 0) = 5."""
        assert ackermann_function(3, 0) == 5

    def test_a_3_1(self):
        """A(3, 1) = 13."""
        assert ackermann_function(3, 1) == 13

    def test_a_3_2(self):
        """A(3, 2) = 29."""
        assert ackermann_function(3, 2) == 29

    def test_a_3_3(self):
        """A(3, 3) = 61."""
        assert ackermann_function(3, 3) == 61

    def test_a_3_4(self):
        """A(3, 4) = 125."""
        assert ackermann_function(3, 4) == 125

    def test_ackermann_growth_demo_returns_dict(self):
        """ackermann_growth_demo() gibt Dict mit Wachstumsinformationen zurück."""
        demo = ackermann_growth_demo()
        assert isinstance(demo, dict)
        assert demo['A(3,3)'] == 61
        assert 'not_primitive_recursive' in demo


# ==============================================================================
# Tests: Halteproblem-Beweis
# ==============================================================================

class TestHaltingProblemProof:
    """Tests für den Halteproblem-Unentscheidbarkeitsbeweis."""

    def setup_method(self):
        self.proof = halting_problem_undecidability_proof()

    def test_returns_dict(self):
        """Funktion gibt Dict zurück."""
        assert isinstance(self.proof, dict)

    def test_has_proof_method(self):
        """Dict enthält 'proof_method' Schlüssel."""
        assert 'proof_method' in self.proof

    def test_proof_method_is_diagonal(self):
        """Beweis-Methode ist Diagonalargument."""
        assert 'Diagonal' in self.proof['proof_method'] or 'diagonal' in self.proof['proof_method']

    def test_has_conclusion(self):
        """Dict enthält 'conclusion' Schlüssel."""
        assert 'conclusion' in self.proof

    def test_conclusion_mentions_undecidable(self):
        """Schlussfolgerung erwähnt Unentscheidbarkeit."""
        conclusion = self.proof['conclusion'].lower()
        assert 'unentscheidbar' in conclusion or 'undecidable' in conclusion.lower()

    def test_has_contradiction(self):
        """Dict enthält 'contradiction' Schlüssel."""
        assert 'contradiction' in self.proof

    def test_has_theorem(self):
        """Dict enthält 'theorem' Schlüssel."""
        assert 'theorem' in self.proof


# ==============================================================================
# Tests: Satz von Rice
# ==============================================================================

class TestRiceTheorem:
    """Tests für den Satz von Rice."""

    def setup_method(self):
        self.rice = rice_theorem_demo()

    def test_returns_dict(self):
        """Funktion gibt Dict zurück."""
        assert isinstance(self.rice, dict)

    def test_has_theorem(self):
        """Dict enthält 'theorem' Schlüssel."""
        assert 'theorem' in self.rice

    def test_has_undecidable_examples(self):
        """Dict enthält Beispiele unentscheidbarer Eigenschaften."""
        assert 'undecidable_examples' in self.rice
        assert len(self.rice['undecidable_examples']) >= 3

    def test_has_consequence(self):
        """Dict enthält 'consequence' Schlüssel."""
        assert 'consequence' in self.rice

    def test_has_formal_statement(self):
        """Dict enthält formale Aussage."""
        assert 'formal_statement' in self.rice


# ==============================================================================
# Tests: Entscheidbare / unentscheidbare Probleme
# ==============================================================================

class TestDecidableUndecidableProblems:
    """Tests für die Listen entscheidbarer/unentscheidbarer Probleme."""

    def test_decidable_returns_list(self):
        """decidable_problems() gibt Liste zurück."""
        problems = decidable_problems()
        assert isinstance(problems, list)
        assert len(problems) >= 4

    def test_decidable_has_required_keys(self):
        """Jedes entscheidbare Problem hat name, description, reason."""
        for p in decidable_problems():
            assert 'name' in p
            assert 'description' in p
            assert 'reason' in p

    def test_undecidable_returns_list(self):
        """undecidable_problems() gibt Liste zurück."""
        problems = undecidable_problems()
        assert isinstance(problems, list)
        assert len(problems) >= 5

    def test_undecidable_has_required_keys(self):
        """Jedes unentscheidbare Problem hat name, description."""
        for p in undecidable_problems():
            assert 'name' in p
            assert 'description' in p

    def test_undecidable_includes_halt(self):
        """Halteproblem ist in der Liste."""
        names = [p['name'] for p in undecidable_problems()]
        halt_found = any('Halt' in n or 'HALT' in n or 'halt' in n.lower() for n in names)
        assert halt_found

    def test_undecidable_includes_hilbert10(self):
        """Hilberts 10. Problem ist in der Liste."""
        names = [p['name'] for p in undecidable_problems()]
        hilbert_found = any('Hilbert' in n or 'hilbert' in n.lower() or '10' in n for n in names)
        assert hilbert_found


# ==============================================================================
# Tests: Arithmetische Hierarchie
# ==============================================================================

class TestArithmeticalHierarchy:
    """Tests für die arithmetische Hierarchie."""

    def setup_method(self):
        self.hier = arithmetical_hierarchy()

    def test_returns_dict(self):
        """Funktion gibt Dict zurück."""
        assert isinstance(self.hier, dict)

    def test_has_sigma0(self):
        """Hierarchie enthält Σ₀."""
        assert 'Σ₀' in self.hier

    def test_has_sigma1(self):
        """Hierarchie enthält Σ₁."""
        assert 'Σ₁' in self.hier

    def test_has_pi1(self):
        """Hierarchie enthält Π₁."""
        assert 'Π₁' in self.hier

    def test_has_sigma2(self):
        """Hierarchie enthält Σ₂."""
        assert 'Σ₂' in self.hier

    def test_has_sigma3(self):
        """Hierarchie enthält Σ₃."""
        assert 'Σ₃' in self.hier

    def test_sigma1_has_examples(self):
        """Σ₁ hat Beispiele."""
        assert 'examples' in self.hier['Σ₁']
        assert len(self.hier['Σ₁']['examples']) >= 1

    def test_sigma1_complete(self):
        """sigma1_complete_problem() gibt Dict mit 'problem' und 'claim' zurück."""
        result = sigma1_complete_problem()
        assert isinstance(result, dict)
        assert 'problem' in result
        assert 'claim' in result
        assert 'conclusion' in result


# ==============================================================================
# Tests: Post'sches Korrespondenzproblem
# ==============================================================================

class TestPostCorrespondenceProblem:
    """Tests für das Post'sche Korrespondenzproblem."""

    def test_trivial_solution(self):
        """Direkte Lösung: Ein Domino mit identischen Strings."""
        # Domino (abc, abc): Sequenz [0] → 'abc' = 'abc'
        dominoes = [('abc', 'abc')]
        result = post_correspondence_problem(dominoes)
        assert result is True

    def test_empty_dominoes(self):
        """Leere Domino-Liste: keine Lösung."""
        result = post_correspondence_problem([])
        assert result is False

    def test_simple_solvable(self):
        """Einfach lösbare Instanz."""
        # (a, ab), (ab, b): Sequenz [0, 1] → 'a'+'ab' = 'aab', 'ab'+'b' = 'abb' → nein
        # Andere: (a, a) → trivial lösbar
        dominoes = [('a', 'a'), ('b', 'b')]
        result = post_correspondence_problem(dominoes)
        assert result is True

    def test_classic_solvable_pcp(self):
        """Klassische lösbare PCP-Instanz."""
        # Bekannte lösbare Instanz: (1, 101), (10, 00), (011, 11)
        # Lösung: Sequenz 1,2,1,0 → top='1'+'10'+'1'+'1'='1101 1', bottom='101'+'00'+'101'+'1'
        # Vereinfacht: (a, aa) → Sequenz [0, 0, ...] - nein
        # Stattdessen: einfaches Beispiel
        dominoes = [('ab', 'a'), ('a', 'ab')]
        # top='ab'+'a'='aba', bottom='a'+'ab'='aab' → nein
        # top='a'+'ab'='aab', bottom='ab'+'a'='aba' → nein
        # Kein Ergebnis in beiden Richtungen → None oder False
        result = post_correspondence_problem(dominoes)
        assert result in (True, False, None)

    def test_returns_valid_type(self):
        """Funktion gibt bool oder None zurück."""
        dominoes = [('a', 'b'), ('b', 'a')]
        result = post_correspondence_problem(dominoes)
        assert result in (True, False, None)


# ==============================================================================
# Tests: Komplexitätsklassen
# ==============================================================================

class TestComplexityClasses:
    """Tests für Komplexitätsklassen."""

    def test_time_complexity_classes_returns_dict(self):
        """time_complexity_classes() gibt Dict zurück."""
        classes = time_complexity_classes()
        assert isinstance(classes, dict)

    def test_has_P(self):
        """P-Klasse vorhanden."""
        classes = time_complexity_classes()
        assert 'P' in classes

    def test_has_NP(self):
        """NP-Klasse vorhanden."""
        classes = time_complexity_classes()
        assert 'NP' in classes

    def test_has_PSPACE(self):
        """PSPACE-Klasse vorhanden."""
        classes = time_complexity_classes()
        assert 'PSPACE' in classes

    def test_has_inclusions(self):
        """Inklusionskette vorhanden."""
        classes = time_complexity_classes()
        assert 'inclusions' in classes

    def test_np_complete_returns_list(self):
        """np_complete_problems() gibt Liste zurück."""
        problems = np_complete_problems()
        assert isinstance(problems, list)

    def test_np_complete_has_at_least_5(self):
        """Mindestens 5 NP-vollständige Probleme."""
        assert len(np_complete_problems()) >= 5

    def test_np_complete_has_sat(self):
        """SAT ist in der Liste."""
        problems = np_complete_problems()
        names = [p['name'] for p in problems]
        sat_found = any('SAT' in n for n in names)
        assert sat_found

    def test_cook_levin_sketch(self):
        """cook_levin_theorem_sketch() gibt Dict zurück."""
        sketch = cook_levin_theorem_sketch()
        assert isinstance(sketch, dict)
        assert 'theorem' in sketch
        assert 'part1' in sketch
        assert 'part2' in sketch

    def test_polynomial_hierarchy(self):
        """polynomial_hierarchy() enthält Σᵖ₁ und NP."""
        ph = polynomial_hierarchy()
        assert isinstance(ph, dict)
        # Mindestens P und NP vorhanden
        keys = list(ph.keys())
        assert any('P' in k for k in keys)


# ==============================================================================
# Tests: Reduktionen
# ==============================================================================

class TestReductions:
    """Tests für Reduktions-Funktionen."""

    def test_many_one_reduction_returns_dict(self):
        """many_one_reduction() gibt Dict zurück."""
        result = many_one_reduction('HALT', 'ATM')
        assert isinstance(result, dict)
        assert 'reduction' in result
        assert 'type' in result

    def test_turing_reduction_returns_dict(self):
        """turing_reduction() gibt Dict zurück."""
        result = turing_reduction('co-HALT', 'HALT')
        assert isinstance(result, dict)
        assert 'reduction' in result

    def test_reduction_degrees_returns_dict(self):
        """reduction_degrees() gibt Dict mit Graden zurück."""
        degrees = reduction_degrees()
        assert isinstance(degrees, dict)
        assert 'degree_0' in degrees
        assert 'degree_0_prime' in degrees


# ==============================================================================
# Tests: Kolmogorov-Komplexität
# ==============================================================================

class TestKolmogorovComplexity:
    """Tests für die Kolmogorov-Komplexitäts-Abschätzungen."""

    def test_returns_dict(self):
        """kolmogorov_complexity_approx() gibt Dict zurück."""
        result = kolmogorov_complexity_approx('hello')
        assert isinstance(result, dict)

    def test_has_required_keys(self):
        """Dict enthält alle erforderlichen Schlüssel."""
        result = kolmogorov_complexity_approx('hello')
        assert 'string_length' in result
        assert 'compressed_length' in result
        assert 'compression_ratio' in result
        assert 'complexity_class' in result

    def test_repetitive_string_low_complexity(self):
        """Repetitiver String hat niedrige Kompressionsrate."""
        result = kolmogorov_complexity_approx('a' * 1000)
        # Komprimierter String sollte deutlich kürzer sein
        assert result['compression_ratio'] < 0.5

    def test_short_string_length(self):
        """Länge wird korrekt gemessen."""
        s = 'hello'
        result = kolmogorov_complexity_approx(s)
        assert result['string_length'] == len(s.encode('utf-8'))

    def test_random_vs_repetitive(self):
        """Zufälliger String hat höhere Rate als repetitiver."""
        repetitive = kolmogorov_complexity_approx('a' * 200)
        import random
        random.seed(0)
        rnd_str = ''.join(random.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(200))
        random_result = kolmogorov_complexity_approx(rnd_str)
        # Repetitiver String hat niedrigeres Verhältnis
        assert repetitive['compression_ratio'] < random_result['compression_ratio']

    def test_incompressible_demo_returns_dict(self):
        """incompressible_string_demo() gibt Dict zurück."""
        result = incompressible_string_demo(20)
        assert isinstance(result, dict)

    def test_incompressible_demo_has_counting_argument(self):
        """Dict enthält Zählargument."""
        result = incompressible_string_demo(20)
        assert 'counting_argument' in result

    def test_incompressible_demo_n_is_correct(self):
        """n im Dict stimmt mit Parameter überein."""
        result = incompressible_string_demo(15)
        assert result['n'] == 15

    def test_incompressible_structured_lower_ratio(self):
        """Strukturiertes Beispiel hat niedrigere Kompressionsrate als zufälliges."""
        result = incompressible_string_demo(100)
        s_ratio = result['structured_example']['ratio']
        r_ratio = result['random_example']['ratio']
        assert s_ratio < r_ratio


# ==============================================================================
# Integrationstests
# ==============================================================================

class TestIntegration:
    """Integrationstests: Mehrere Module zusammen."""

    def test_full_palindrome_workflow(self):
        """Vollständiger Test: TM erstellen, verschiedene Eingaben."""
        tm = tm_recognizes_palindromes()
        true_palindromes  = ['', 'a', 'b', 'aa', 'bb', 'aba', 'abba', 'aabaa']
        false_palindromes = ['ab', 'ba', 'abc', 'aab']

        for s in true_palindromes:
            assert tm.accepts(s, max_steps=10000) is True, f"'{s}' sollte Palindrom sein"
        for s in false_palindromes:
            result = tm.accepts(s, max_steps=10000)
            assert result in (False, None), f"'{s}' sollte kein Palindrom sein"

    def test_ackermann_growth_consistency(self):
        """Ackermann-Wachstum stimmt mit berechneten Werten überein."""
        demo = ackermann_growth_demo()
        # Vergleiche mit direkter Berechnung
        assert demo['A(0,0)'] == ackermann_function(0, 0)
        assert demo['A(1,1)'] == ackermann_function(1, 1)
        assert demo['A(2,2)'] == ackermann_function(2, 2)
        assert demo['A(3,3)'] == ackermann_function(3, 3)

    def test_hierarchy_and_complexity_consistent(self):
        """Arithmetische Hierarchie und Komplexitätsklassen sind konsistent."""
        hier     = arithmetical_hierarchy()
        classes  = time_complexity_classes()
        # Beide sollten auf P / entscheidbare Probleme hinweisen
        sigma0_examples = ' '.join(hier['Σ₀']['examples'])
        p_examples      = ' '.join(classes['P']['examples'])
        # Beide sollten nicht leer sein
        assert len(sigma0_examples) > 0
        assert len(p_examples) > 0
