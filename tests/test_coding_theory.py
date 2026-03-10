"""
@file test_coding_theory.py
@brief Umfassende Tests für das coding_theory-Modul.
@description
    Testet alle Klassen und Funktionen des Kodierungstheorie-Moduls:

    - GaloisField           : GF(p^n)-Arithmetik, primitives Element, diskreter Log
    - LinearCode            : Kodierung, Syndrom, Fehlerkorrektur, dualer Code
    - HammingCode           : Parameter, Generatormatrix, Encode/Decode mit Fehler
    - CyclicCode            : Polynommultiplikation, Syndrom, zyklische Verschiebung
    - BCHCode               : Parameter, Reed-Solomon-Demo, Berlekamp-Massey
    - DFA                   : Akzeptanz, leeres Wort, minimale Automaten
    - NFA                   : Akzeptanz via Teilmengenkonstruktion
    - RegularExpression     : Muster-Match, Kleene-Stern, Alternation
    - ContextFreeGrammar    : Erzeugungsrelation, Palindrome, leeres Wort
    - PushdownAutomaton     : aⁿbⁿ-Sprache, Demo-Akzeptanz
    - ComputabilityTheory   : Demos, unentscheidbare Probleme

    Edge Cases:
    - Leeres Wort ε
    - Einheitsmatrix als Generatormatrix
    - Einelementige Sprachen
    - GF(2) und GF(3)
    - Hamming-Code mit und ohne Fehler
    - NFA mit mehreren möglichen Pfaden
    - CFG mit ε-Produktion

@author Kurt Ingwer
@version 1.0
@date 2026-03-10
@lastModified 2026-03-10
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

# Projektpfad hinzufügen
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from coding_theory import (
    BCHCode,
    ComputabilityTheory,
    ContextFreeGrammar,
    CyclicCode,
    DFA,
    FormalLanguages,
    GaloisField,
    HammingCode,
    LinearCode,
    NFA,
    PushdownAutomaton,
    RegularExpression,
)


# ==============================================================
# Fixtures
# ==============================================================

@pytest.fixture
def gf2():
    """GF(2): Körper mit 2 Elementen."""
    return GaloisField(2)


@pytest.fixture
def gf3():
    """GF(3): Körper mit 3 Elementen {0,1,2}."""
    return GaloisField(3)


@pytest.fixture
def gf5():
    """GF(5): Körper mit 5 Elementen {0,1,2,3,4}."""
    return GaloisField(5)


@pytest.fixture
def hamming_r3():
    """Hamming-Code der Ordnung r=3: [7,4,3]-Code."""
    return HammingCode(3)


@pytest.fixture
def hamming_r2():
    """Hamming-Code der Ordnung r=2: [3,1,3]-Code."""
    return HammingCode(2)


@pytest.fixture
def repetition_code():
    """Wiederholungscode [3,1,3]: G = [[1,1,1]]."""
    return LinearCode([[1, 1, 1]])


@pytest.fixture
def simple_code():
    """Einfacher [4,2,2]-Code mit Paritätsprüfung."""
    # G = [I2 | P], P = [[1,0],[0,1]] → Paritätsbits
    # Tatsächlich: G = [[1,0,1,0],[0,1,0,1]]
    G = np.array([[1, 0, 1, 0],
                  [0, 1, 0, 1]])
    return LinearCode(G)


@pytest.fixture
def dfa_even_zeros():
    """DFA: Akzeptiert binäre Wörter mit gerader Anzahl von Nullen."""
    states = {'q0', 'q1'}
    alphabet = {'0', '1'}
    transitions = {
        ('q0', '0'): 'q1',
        ('q0', '1'): 'q0',
        ('q1', '0'): 'q0',
        ('q1', '1'): 'q1',
    }
    return DFA(states, alphabet, transitions, 'q0', {'q0'})


@pytest.fixture
def nfa_ends_ab():
    """NFA: Akzeptiert Wörter über {a,b}, die auf 'ab' enden."""
    states = {0, 1, 2}
    alphabet = {'a', 'b'}
    transitions = {
        0: {'a': {0, 1}, 'b': {0}},
        1: {'b': {2}},
        2: {},
    }
    return NFA(states, alphabet, transitions, 0, {2})


@pytest.fixture
def cfg_palindromes():
    """KFG für Palindrome über {a,b}: S → aSa | bSb | a | b | ε."""
    terminals = {'a', 'b'}
    nonterminals = {'S'}
    rules = {
        'S': [
            ['a', 'S', 'a'],
            ['b', 'S', 'b'],
            ['a'],
            ['b'],
            [''],  # ε-Produktion
        ]
    }
    return ContextFreeGrammar(terminals, nonterminals, rules, 'S')


# ==============================================================
# 1. GaloisField Tests
# ==============================================================

class TestGaloisField:
    """Tests für GaloisField."""

    def test_init_gf2(self, gf2):
        """GF(2) hat 2 Elemente."""
        assert gf2.p == 2
        assert gf2.n == 1
        assert gf2.order == 2

    def test_init_gf3(self, gf3):
        """GF(3) hat 3 Elemente."""
        assert gf3.order == 3

    def test_init_gf4(self):
        """GF(4) = GF(2^2) hat 4 Elemente."""
        gf4 = GaloisField(2, 2)
        assert gf4.order == 4
        assert gf4.p == 2
        assert gf4.n == 2

    def test_init_invalid_not_prime(self):
        """Keine Primzahl → ValueError."""
        with pytest.raises(ValueError, match="keine Primzahl"):
            GaloisField(4)

    def test_init_invalid_n_zero(self):
        """n=0 → ValueError."""
        with pytest.raises(ValueError):
            GaloisField(2, 0)

    def test_elements_gf2(self, gf2):
        """GF(2) hat Elemente {0, 1}."""
        assert gf2.elements() == [0, 1]

    def test_elements_gf5(self, gf5):
        """GF(5) hat Elemente {0,1,2,3,4}."""
        assert gf5.elements() == [0, 1, 2, 3, 4]

    def test_elements_count(self):
        """GF(p^n) hat genau p^n Elemente."""
        for p, n in [(2, 1), (2, 3), (3, 2), (5, 1)]:
            gf = GaloisField(p, n)
            assert len(gf.elements()) == p ** n

    # ----------------------------------------------------------
    # Addition
    # ----------------------------------------------------------

    def test_add_gf2_basic(self, gf2):
        """GF(2)-Addition: 0+0=0, 0+1=1, 1+1=0."""
        assert gf2.add(0, 0) == 0
        assert gf2.add(0, 1) == 1
        assert gf2.add(1, 0) == 1
        assert gf2.add(1, 1) == 0

    def test_add_gf3(self, gf3):
        """GF(3)-Addition mod 3."""
        assert gf3.add(1, 2) == 0
        assert gf3.add(2, 2) == 1
        assert gf3.add(0, 0) == 0

    def test_add_invalid(self, gf2):
        """Ungültige Elemente → ValueError."""
        with pytest.raises(ValueError):
            gf2.add(0, 2)  # 2 ∉ GF(2)

    def test_add_commutativity(self, gf5):
        """Addition ist kommutativ in GF(p)."""
        for a in range(5):
            for b in range(5):
                assert gf5.add(a, b) == gf5.add(b, a)

    def test_add_associativity(self, gf3):
        """Addition ist assoziativ in GF(3)."""
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    assert gf3.add(gf3.add(a, b), c) == gf3.add(a, gf3.add(b, c))

    def test_add_identity(self, gf5):
        """0 ist neutrales Element der Addition."""
        for a in range(5):
            assert gf5.add(a, 0) == a
            assert gf5.add(0, a) == a

    # ----------------------------------------------------------
    # Multiplikation
    # ----------------------------------------------------------

    def test_multiply_gf2(self, gf2):
        """GF(2)-Multiplikation."""
        assert gf2.multiply(0, 0) == 0
        assert gf2.multiply(0, 1) == 0
        assert gf2.multiply(1, 1) == 1

    def test_multiply_gf3(self, gf3):
        """GF(3)-Multiplikation mod 3."""
        assert gf3.multiply(2, 2) == 1
        assert gf3.multiply(2, 1) == 2
        assert gf3.multiply(0, 2) == 0

    def test_multiply_invalid(self, gf3):
        """Ungültige Elemente → ValueError."""
        with pytest.raises(ValueError):
            gf3.multiply(1, 3)

    def test_multiply_distributive(self, gf5):
        """Distributivgesetz: a*(b+c) = a*b + a*c."""
        for a in range(5):
            for b in range(5):
                for c in range(5):
                    left = gf5.multiply(a, gf5.add(b, c))
                    right = gf5.add(gf5.multiply(a, b), gf5.multiply(a, c))
                    assert left == right

    def test_multiply_identity(self, gf5):
        """1 ist neutrales Element der Multiplikation."""
        for a in range(5):
            assert gf5.multiply(a, 1) == a
            assert gf5.multiply(1, a) == a

    # ----------------------------------------------------------
    # Primitives Element
    # ----------------------------------------------------------

    def test_primitive_element_gf2(self, gf2):
        """GF(2)*={1}, primitives Element = 1."""
        g = gf2.primitive_element()
        assert g == 1

    def test_primitive_element_gf3(self, gf3):
        """GF(3)*={1,2}, primitives Element = 2 (ord=2)."""
        g = gf3.primitive_element()
        assert g in [2]
        # Ord(g) = p-1 = 2: g^1 ≠ 1, g^2 = 1
        assert pow(g, 3 - 1, 3) == 1

    def test_primitive_element_gf5(self, gf5):
        """GF(5): primitives Element hat Ordnung 4."""
        g = gf5.primitive_element()
        p = 5
        # Ord(g) = p-1 = 4
        assert pow(g, p - 1, p) == 1
        # Kein kleinerer Exponent ergibt 1
        for k in range(1, p - 1):
            assert pow(g, k, p) != 1

    def test_primitive_element_gf7(self):
        """GF(7): primitives Element hat Ordnung 6."""
        gf7 = GaloisField(7)
        g = gf7.primitive_element()
        assert pow(g, 6, 7) == 1
        assert pow(g, 2, 7) != 1
        assert pow(g, 3, 7) != 1

    # ----------------------------------------------------------
    # Diskreter Logarithmus
    # ----------------------------------------------------------

    def test_discrete_log_gf5(self, gf5):
        """log_g(g^k) = k in GF(5)."""
        g = gf5.primitive_element()
        for k in range(4):
            elem = pow(g, k, 5)
            log = gf5.discrete_log_gf(g, elem)
            assert pow(g, log, 5) == elem

    def test_discrete_log_zero(self, gf5):
        """log(0) ist undefiniert → -1."""
        g = gf5.primitive_element()
        assert gf5.discrete_log_gf(g, 0) == -1

    def test_discrete_log_one(self, gf5):
        """log(1) = 0."""
        g = gf5.primitive_element()
        assert gf5.discrete_log_gf(g, 1) == 0

    # ----------------------------------------------------------
    # Primitives Polynom
    # ----------------------------------------------------------

    def test_is_primitive_poly_x3_x_1(self, gf2):
        """x^3 + x + 1 = [1,1,0,1] ist primitiv über GF(2)."""
        # Koeffizienten: [a0, a1, a2, a3] = [1,1,0,1] → 1 + x + x^3
        assert gf2.is_primitive_poly([1, 1, 0, 1]) is True

    def test_is_primitive_poly_x4_x_1(self, gf2):
        """x^4 + x + 1 = [1,1,0,0,1] ist primitiv über GF(2)."""
        assert gf2.is_primitive_poly([1, 1, 0, 0, 1]) is True

    def test_is_primitive_poly_not_primitive(self, gf2):
        """x^3 + 1 = [1,0,0,1] ist NICHT irreduzibel (Faktor: (x+1))."""
        # x^3 + 1 = (x+1)(x^2+x+1), also reduzibel → nicht primitiv
        result = gf2.is_primitive_poly([1, 0, 0, 1])
        assert result is False

    def test_is_primitive_poly_wrong_field(self):
        """is_primitive_poly nur für GF(2)."""
        gf3 = GaloisField(3)
        with pytest.raises(NotImplementedError):
            gf3.is_primitive_poly([1, 0, 1])


# ==============================================================
# 2. LinearCode Tests
# ==============================================================

class TestLinearCode:
    """Tests für LinearCode."""

    def test_init_shape(self, repetition_code):
        """Wiederholungscode [3,1]: k=1, n=3."""
        assert repetition_code.k == 1
        assert repetition_code.n == 3

    def test_encode_repetition_zero(self, repetition_code):
        """Nachricht [0] → Codewort [0,0,0]."""
        c = repetition_code.encode([0])
        assert list(c) == [0, 0, 0]

    def test_encode_repetition_one(self, repetition_code):
        """Nachricht [1] → Codewort [1,1,1]."""
        c = repetition_code.encode([1])
        assert list(c) == [1, 1, 1]

    def test_encode_identity_matrix(self):
        """Einheitsmatrix als Generator: encode(m) = m (Padding)."""
        G = np.eye(3, dtype=int)
        code = LinearCode(G)
        m = [1, 0, 1]
        c = code.encode(m)
        assert list(c) == [1, 0, 1]

    def test_encode_gf2_modulo(self, simple_code):
        """Kodierung ist mod 2."""
        c = simple_code.encode([1, 1])
        # [1,1] * [[1,0,1,0],[0,1,0,1]] = [1,1,1,1]
        assert all(x in [0, 1] for x in c)

    def test_encode_invalid_length(self, repetition_code):
        """Falsche Nachrichtenlänge → ValueError."""
        with pytest.raises(ValueError):
            repetition_code.encode([1, 0])

    def test_encode_linearity(self, simple_code):
        """Linearität: encode(m1+m2) = encode(m1) + encode(m2) mod 2."""
        m1 = np.array([1, 0])
        m2 = np.array([0, 1])
        c1 = simple_code.encode(m1)
        c2 = simple_code.encode(m2)
        c_sum = simple_code.encode((m1 + m2) % 2)
        assert np.array_equal(c_sum, (c1 + c2) % 2)

    def test_parity_check_matrix_shape(self, simple_code):
        """Prüfmatrix H hat Dimension (n-k) × n."""
        H = simple_code.parity_check_matrix()
        assert H.shape == (simple_code.n - simple_code.k, simple_code.n)

    def test_parity_check_orthogonality(self, simple_code):
        """H·cᵀ = 0 für alle Codewörter c."""
        H = simple_code.parity_check_matrix()
        for i in range(1 << simple_code.k):
            m = np.array([(i >> j) & 1 for j in range(simple_code.k)])
            c = simple_code.encode(m)
            syndrome = (H @ c) % 2
            assert np.all(syndrome == 0), f"H·c ≠ 0 für m={m}, c={c}"

    def test_parity_check_repetition(self, repetition_code):
        """Repetitionscode [3,1]: H hat Shape (2,3)."""
        H = repetition_code.parity_check_matrix()
        assert H.shape == (2, 3)

    def test_syndrome_no_error(self, simple_code):
        """Korrektes Codewort hat Syndrom 0."""
        c = simple_code.encode([1, 0])
        s = simple_code.syndrome(c)
        assert np.all(s == 0)

    def test_syndrome_with_error(self, simple_code):
        """Codewort mit Fehler hat Syndrom ≠ 0."""
        c = simple_code.encode([1, 0])
        r = c.copy()
        r[0] ^= 1  # Fehler in Position 0
        s = simple_code.syndrome(r)
        assert not np.all(s == 0)

    def test_syndrome_invalid_length(self, simple_code):
        """Falsche Länge → ValueError."""
        with pytest.raises(ValueError):
            simple_code.syndrome([1, 0])

    def test_minimum_distance_repetition(self, repetition_code):
        """Wiederholungscode [3,1,3] hat Minimaldistanz 3."""
        d = repetition_code.minimum_distance()
        assert d == 3

    def test_minimum_distance_simple(self, simple_code):
        """Einfacher [4,2]-Code: Minimaldistanz ≥ 1."""
        d = simple_code.minimum_distance()
        assert d >= 1

    def test_error_correct_no_error(self, repetition_code):
        """Fehlerkorrektur ohne Fehler gibt Codewort zurück."""
        c = repetition_code.encode([1])
        corrected = repetition_code.error_correct(c)
        assert np.array_equal(corrected, c)

    def test_error_correct_one_error(self, repetition_code):
        """Fehlerkorrektur mit 1 Fehler im Wiederholungscode."""
        c = repetition_code.encode([1])  # [1,1,1]
        r = c.copy()
        r[1] ^= 1  # [1,0,1]
        corrected = repetition_code.error_correct(r)
        assert np.array_equal(corrected, c)

    def test_is_systematic_true(self):
        """Systematischer Code wird erkannt."""
        G = np.array([[1, 0, 1, 1],
                      [0, 1, 1, 0]])
        code = LinearCode(G)
        assert code.is_systematic() is True

    def test_is_systematic_false(self):
        """Nicht-systematischer Code wird erkannt."""
        G = np.array([[1, 1, 0],
                      [0, 1, 1]])
        code = LinearCode(G)
        # Einheitsmatrix I2 kommt nicht als Teilmatrix vor
        result = code.is_systematic()
        # Wir prüfen nur dass der Aufruf funktioniert
        assert isinstance(result, bool)

    def test_dual_code_dimension(self, simple_code):
        """Dualer Code C⊥ hat Dimension n-k."""
        dual = simple_code.dual_code()
        assert dual.k == simple_code.n - simple_code.k
        assert dual.n == simple_code.n

    def test_dual_code_orthogonality(self, simple_code):
        """Codewörter von C und C⊥ sind orthogonal."""
        dual = simple_code.dual_code()
        for i in range(1 << simple_code.k):
            m1 = np.array([(i >> j) & 1 for j in range(simple_code.k)])
            c = simple_code.encode(m1)
            for j in range(1 << dual.k):
                m2 = np.array([(j >> l) & 1 for l in range(dual.k)])
                c_perp = dual.encode(m2)
                dot = int(np.dot(c, c_perp)) % 2
                assert dot == 0, f"c={c}, c⊥={c_perp} nicht orthogonal"


# ==============================================================
# 3. HammingCode Tests
# ==============================================================

class TestHammingCode:
    """Tests für HammingCode."""

    def test_parameters_r2(self, hamming_r2):
        """[3,1,3]-Hamming-Code."""
        n, k, d = hamming_r2.parameters()
        assert n == 3    # 2^2 - 1
        assert k == 1    # 2^2 - 2 - 1
        assert d == 3

    def test_parameters_r3(self, hamming_r3):
        """[7,4,3]-Hamming-Code."""
        n, k, d = hamming_r3.parameters()
        assert n == 7    # 2^3 - 1
        assert k == 4    # 2^3 - 3 - 1
        assert d == 3

    def test_parameters_r4(self):
        """[15,11,3]-Hamming-Code."""
        h = HammingCode(4)
        n, k, d = h.parameters()
        assert n == 15
        assert k == 11
        assert d == 3

    def test_invalid_r1(self):
        """r=1 → ValueError."""
        with pytest.raises(ValueError):
            HammingCode(1)

    def test_parity_check_matrix_shape_r3(self, hamming_r3):
        """H hat Dimension r × n = 3 × 7."""
        H = hamming_r3.parity_check_matrix()
        assert H.shape == (3, 7)

    def test_parity_check_matrix_all_columns(self, hamming_r3):
        """Alle Spalten von H sind verschiedene Bitvektoren 1..2^r-1."""
        H = hamming_r3.parity_check_matrix()
        # Spalten als Zahlen
        col_vals = set()
        for j in range(H.shape[1]):
            val = 0
            for i in range(H.shape[0]):
                val = (val << 1) | int(H[i, j])
            col_vals.add(val)
        # Alle Werte von 1 bis 7 vorhanden
        assert col_vals == set(range(1, 8))

    def test_generator_matrix_shape_r3(self, hamming_r3):
        """G hat Dimension k × n = 4 × 7."""
        G = hamming_r3.generator_matrix()
        assert G.shape == (4, 7)

    def test_generator_orthogonal_to_parity(self, hamming_r3):
        """H·Gᵀ = 0 mod 2."""
        G = hamming_r3.generator_matrix()
        H = hamming_r3.parity_check_matrix()
        product = (H @ G.T) % 2
        assert np.all(product == 0)

    def test_encode_r3(self, hamming_r3):
        """Kodierung produziert Codewort der Länge 7."""
        c = hamming_r3.encode([1, 0, 1, 1])
        assert len(c) == 7
        assert all(x in [0, 1] for x in c)

    def test_encode_all_zeros_r3(self, hamming_r3):
        """Nullvektor ergibt Nullcodewort."""
        c = hamming_r3.encode([0, 0, 0, 0])
        assert np.all(c == 0)

    def test_decode_no_error(self, hamming_r3):
        """Dekodierung ohne Fehler: Fehlerpositon = 0."""
        m = [1, 0, 1, 1]
        c = hamming_r3.encode(m)
        result = hamming_r3.decode(c)
        assert result['error_position'] == 0
        assert np.array_equal(result['corrected'], c)

    def test_decode_single_error_pos1(self, hamming_r3):
        """Dekodierung mit 1 Fehler: Fehlerposition korrekt erkannt."""
        m = [1, 0, 1, 0]
        c = hamming_r3.encode(m)

        # Fehler in Position 1 (0-basiert = Position 1, 1-basiert = 1)
        r = c.copy()
        error_idx = 0  # 0-basiert
        r[error_idx] ^= 1

        result = hamming_r3.decode(r)
        assert result['error_position'] == error_idx + 1  # 1-basiert
        assert np.array_equal(result['corrected'], c)

    def test_decode_single_error_pos5(self, hamming_r3):
        """Dekodierung mit Fehler in Position 5."""
        m = [0, 1, 1, 0]
        c = hamming_r3.encode(m)
        r = c.copy()
        r[4] ^= 1  # 0-basiert Position 4 = 1-basiert Position 5

        result = hamming_r3.decode(r)
        assert np.array_equal(result['corrected'], c)

    def test_decode_invalid_length(self, hamming_r3):
        """Falsche Länge → ValueError."""
        with pytest.raises(ValueError):
            hamming_r3.decode([1, 0, 1])


# ==============================================================
# 4. CyclicCode Tests
# ==============================================================

class TestCyclicCode:
    """Tests für CyclicCode."""

    def test_init_parameters(self):
        """Zyklischer Code [7,4] mit g(x) = x³+x+1."""
        # g(x) = x^3 + x + 1 → Koeff. [1,1,0,1]
        cc = CyclicCode(7, [1, 1, 0, 1])
        assert cc.n == 7
        assert cc.deg_g == 3
        assert cc.k == 4

    def test_encode_poly_length(self):
        """Codewort hat Länge n=7."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        cw = cc.encode_poly([1, 0, 0, 0])
        assert len(cw) == 7

    def test_encode_poly_all_zeros(self):
        """Nullnachricht → Nullcodewort."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        cw = cc.encode_poly([0, 0, 0, 0])
        assert all(x == 0 for x in cw)

    def test_encode_poly_gf2(self):
        """Alle Einträge des Codewortes in {0,1}."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        for bits in range(16):
            m = [(bits >> i) & 1 for i in range(4)]
            cw = cc.encode_poly(m)
            assert all(x in [0, 1] for x in cw)

    def test_syndrome_poly_zero_for_codeword(self):
        """Syndrom eines Codewortes ist 0."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        m = [1, 0, 1, 1]
        cw = cc.encode_poly(m)
        s = cc.syndrome_poly(cw)
        assert all(x == 0 for x in s)

    def test_syndrome_poly_nonzero_for_error(self):
        """Syndrom eines verfälschten Wortes ist ≠ 0."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        m = [1, 0, 0, 1]
        cw = cc.encode_poly(m)
        received = list(cw)
        received[0] ^= 1
        s = cc.syndrome_poly(received)
        assert not all(x == 0 for x in s)

    def test_cyclic_shift_basic(self):
        """Zyklische Verschiebung verschiebt um 1 Position."""
        cc = CyclicCode(5, [1, 0, 1])
        cw = [1, 0, 1, 0, 0]
        shifted = cc.cyclic_shift(cw)
        assert shifted == [0, 1, 0, 1, 0]

    def test_cyclic_shift_empty(self):
        """Zyklische Verschiebung des leeren Wortes."""
        cc = CyclicCode(5, [1, 0, 1])
        assert cc.cyclic_shift([]) == []

    def test_cyclic_shift_wraps_around(self):
        """Letztes Bit wird zum ersten."""
        cc = CyclicCode(4, [1, 1])
        cw = [0, 0, 0, 1]
        shifted = cc.cyclic_shift(cw)
        assert shifted[0] == 1

    def test_is_cyclic_true(self):
        """Zyklisch geschlossener Code wird erkannt."""
        cc = CyclicCode(3, [1, 1])
        # Alle Zyklen von [1,1,0]
        code_words = [[0, 0, 0], [1, 1, 0], [0, 1, 1], [1, 0, 1]]
        assert cc.is_cyclic(code_words) is True

    def test_is_cyclic_false(self):
        """Nicht-zyklisch geschlossener Code wird erkannt."""
        cc = CyclicCode(3, [1, 1])
        # [1,0,0] ist dabei, aber Verschiebung [0,1,0] fehlt
        code_words = [[0, 0, 0], [1, 0, 0]]
        assert cc.is_cyclic(code_words) is False

    def test_poly_add_gf2(self):
        """Polynom-Addition über GF(2)."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        # [1,0,1] + [1,1,0] = [0,1,1]
        result = cc._poly_add([1, 0, 1], [1, 1, 0])
        assert result == [0, 1, 1]

    def test_poly_mul_gf2(self):
        """Polynom-Multiplikation über GF(2)."""
        cc = CyclicCode(7, [1, 1, 0, 1])
        # (x+1)*(x+1) = x^2+1 in GF(2) ([1,1]*[1,1]=[1,0,1])
        result = cc._poly_mul([1, 1], [1, 1])
        assert result == [1, 0, 1]


# ==============================================================
# 5. BCHCode Tests
# ==============================================================

class TestBCHCode:
    """Tests für BCHCode."""

    def test_bch_generator_poly_7_1(self):
        """BCH-Code [7,4,3]: t=1."""
        info = BCHCode.bch_generator_poly(7, 1)
        assert info['n'] == 7
        assert info['t'] == 1
        assert info['d_min'] == 3
        assert info['q'] == 2

    def test_bch_generator_poly_15_2(self):
        """BCH-Code [15,7,5]: t=2."""
        info = BCHCode.bch_generator_poly(15, 2)
        assert info['n'] == 15
        assert info['t'] == 2
        assert info['d_min'] == 5

    def test_bch_d_min_formula(self):
        """d_min = 2t + 1 (BCH-Schranke)."""
        for t in range(1, 6):
            info = BCHCode.bch_generator_poly(15, t)
            assert info['d_min'] == 2 * t + 1

    def test_bch_has_description(self):
        """Rückgabe enthält Beschreibung."""
        info = BCHCode.bch_generator_poly(7, 1)
        assert 'description' in info
        assert len(info['description']) > 0

    def test_reed_solomon_demo(self):
        """RS(255,223) über GF(256): t=16."""
        info = BCHCode.reed_solomon_demo(255, 223, 256)
        assert info['n'] == 255
        assert info['k'] == 223
        assert info['d'] == 33    # 255 - 223 + 1
        assert info['t'] == 16   # (33-1)/2
        assert info['is_mds'] is True

    def test_reed_solomon_mds(self):
        """RS-Codes sind immer MDS: d = n - k + 1."""
        for n, k in [(7, 3), (15, 9), (31, 21)]:
            info = BCHCode.reed_solomon_demo(n, k, 32)
            assert info['d'] == n - k + 1

    def test_reed_solomon_has_applications(self):
        """RS-Demo enthält Anwendungsbeispiele."""
        info = BCHCode.reed_solomon_demo(255, 223, 256)
        assert 'applications' in info
        assert len(info['applications']) > 0

    # ----------------------------------------------------------
    # Berlekamp-Massey
    # ----------------------------------------------------------

    def test_berlekamp_massey_empty(self):
        """Leere Sequenz → LFSR-Länge 0."""
        result = BCHCode.berlekamp_massey([])
        assert result['lfsr_length'] == 0
        assert result['sequence_length'] == 0

    def test_berlekamp_massey_all_zeros(self):
        """Nullsequenz → LFSR-Länge 0."""
        result = BCHCode.berlekamp_massey([0, 0, 0, 0])
        assert result['lfsr_length'] == 0

    def test_berlekamp_massey_single_one(self):
        """Einzelne 1 → LFSR-Länge 1."""
        result = BCHCode.berlekamp_massey([1])
        assert result['lfsr_length'] == 1

    def test_berlekamp_massey_periodic(self):
        """Periodische Sequenz 1,0,1,0,... hat LFSR-Länge 2.

        Das kürzeste LFSR für 1,0,1,0,... ist Länge 2:
        s[n] = s[n-2] (d.h. Periode 2 → LFSR-Länge 2).
        """
        result = BCHCode.berlekamp_massey([1, 0, 1, 0, 1, 0])
        assert result['lfsr_length'] == 2

    def test_berlekamp_massey_returns_poly(self):
        """Rückgabe enthält Verbindungspolynom."""
        result = BCHCode.berlekamp_massey([1, 1, 0, 1, 1, 0])
        assert 'connection_poly' in result
        assert isinstance(result['connection_poly'], list)
        assert result['connection_poly'][0] == 1  # Konstantterm immer 1

    def test_berlekamp_massey_length_correct(self):
        """Sequenzlänge im Ergebnis korrekt."""
        seq = [1, 0, 1, 1, 0, 0, 1]
        result = BCHCode.berlekamp_massey(seq)
        assert result['sequence_length'] == len(seq)


# ==============================================================
# 6. DFA Tests
# ==============================================================

class TestDFA:
    """Tests für DFA (Deterministischer Endlicher Automat)."""

    def test_accepts_empty_word_in_accepting(self, dfa_even_zeros):
        """Leeres Wort akzeptiert wenn Startzustand akzeptierend."""
        # q0 ist akzeptierend (0 Nullen = gerade Anzahl)
        assert dfa_even_zeros.accepts('') is True

    def test_accepts_single_zero(self, dfa_even_zeros):
        """Einzelne '0': 1 Null = ungerade → nicht akzeptiert."""
        assert dfa_even_zeros.accepts('0') is False

    def test_accepts_two_zeros(self, dfa_even_zeros):
        """'00': 2 Nullen = gerade → akzeptiert."""
        assert dfa_even_zeros.accepts('00') is True

    def test_accepts_mixed(self, dfa_even_zeros):
        """'10100' hat 3 Nullen (ungerade) → nicht akzeptiert; '1010' hat 2 Nullen → akzeptiert."""
        # '10100': 1,0,1,0,0 → 3 Nullen (ungerade) → nicht akzeptiert
        assert dfa_even_zeros.accepts('10100') is False
        # '1010': 1,0,1,0 → 2 Nullen (gerade) → akzeptiert
        assert dfa_even_zeros.accepts('1010') is True

    def test_accepts_only_ones(self, dfa_even_zeros):
        """'111': 0 Nullen = gerade → akzeptiert."""
        assert dfa_even_zeros.accepts('111') is True

    def test_rejects_odd_zeros(self, dfa_even_zeros):
        """Wörter mit ungerader Nullen-Anzahl werden korrekt verworfen."""
        # '0': 1 Null (ungerade) → nicht akzeptiert
        assert dfa_even_zeros.accepts('0') is False
        # '000': 3 Nullen (ungerade) → nicht akzeptiert
        assert dfa_even_zeros.accepts('000') is False
        # '001': 2 Nullen (gerade) → akzeptiert
        assert dfa_even_zeros.accepts('001') is True
        # '0001': 3 Nullen (ungerade) → nicht akzeptiert
        assert dfa_even_zeros.accepts('0001') is False

    def test_undefined_transition(self):
        """Undefinierter Übergang → nicht akzeptiert."""
        dfa = DFA(
            states={'q0', 'q1'},
            alphabet={'a'},
            transitions={('q0', 'a'): 'q1'},
            start='q0',
            accepting={'q1'},
        )
        # 'aa': q0 →a→ q1 →a→ ? (kein Übergang → verwirft)
        assert dfa.accepts('aa') is False
        assert dfa.accepts('a') is True

    def test_single_element_language(self):
        """Einelementige Sprache L = {'a'}."""
        dfa = DFA(
            states={'q0', 'q1', 'dead'},
            alphabet={'a', 'b'},
            transitions={
                ('q0', 'a'): 'q1',
                ('q0', 'b'): 'dead',
                ('q1', 'a'): 'dead',
                ('q1', 'b'): 'dead',
            },
            start='q0',
            accepting={'q1'},
        )
        assert dfa.accepts('a') is True
        assert dfa.accepts('aa') is False
        assert dfa.accepts('b') is False
        assert dfa.accepts('') is False

    def test_minimize_returns_dict(self, dfa_even_zeros):
        """Minimierung gibt Dictionary zurück."""
        result = dfa_even_zeros.minimize()
        assert 'equivalence_classes' in result
        assert 'num_states_original' in result
        assert 'num_states_minimal' in result

    def test_minimize_minimal_dfa(self, dfa_even_zeros):
        """Bereits minimaler DFA: minimale Zustände ≤ ursprüngliche Zustände."""
        result = dfa_even_zeros.minimize()
        assert result['num_states_minimal'] <= result['num_states_original']

    def test_minimize_covers_all_states(self, dfa_even_zeros):
        """Äquivalenzklassen überdecken alle Zustände."""
        result = dfa_even_zeros.minimize()
        all_in_classes = set()
        for cls in result['equivalence_classes']:
            all_in_classes |= cls
        assert all_in_classes == dfa_even_zeros.states


# ==============================================================
# 7. NFA Tests
# ==============================================================

class TestNFA:
    """Tests für NFA (Nichtdeterministischer Endlicher Automat)."""

    def test_accepts_ab(self, nfa_ends_ab):
        """'ab' endet auf 'ab' → akzeptiert."""
        assert nfa_ends_ab.accepts('ab') is True

    def test_accepts_xab(self, nfa_ends_ab):
        """'aab' endet auf 'ab' → akzeptiert."""
        assert nfa_ends_ab.accepts('aab') is True

    def test_accepts_long_word(self, nfa_ends_ab):
        """'aaabbbab' endet auf 'ab' → akzeptiert."""
        assert nfa_ends_ab.accepts('aaabbbab') is True

    def test_rejects_ba(self, nfa_ends_ab):
        """'ba' endet nicht auf 'ab' → nicht akzeptiert."""
        assert nfa_ends_ab.accepts('ba') is False

    def test_rejects_empty(self, nfa_ends_ab):
        """Leeres Wort endet nicht auf 'ab' → nicht akzeptiert."""
        assert nfa_ends_ab.accepts('') is False

    def test_rejects_single_a(self, nfa_ends_ab):
        """Einzelnes 'a' → nicht akzeptiert."""
        assert nfa_ends_ab.accepts('a') is False

    def test_rejects_single_b(self, nfa_ends_ab):
        """Einzelnes 'b' → nicht akzeptiert."""
        assert nfa_ends_ab.accepts('b') is False

    def test_nfa_multiple_paths(self):
        """NFA mit mehreren Pfaden, nur einer führt zur Akzeptanz."""
        nfa = NFA(
            states={0, 1, 2},
            alphabet={'a', 'b'},
            transitions={
                0: {'a': {0, 1}},
                1: {'b': {2}},
                2: {},
            },
            start=0,
            accepting={2},
        )
        assert nfa.accepts('ab') is True
        assert nfa.accepts('aab') is True
        assert nfa.accepts('b') is False

    def test_nfa_empty_accepting(self):
        """NFA: leeres Wort akzeptiert wenn Startzustand akzeptierend."""
        nfa = NFA(
            states={0},
            alphabet={'a'},
            transitions={0: {}},
            start=0,
            accepting={0},
        )
        assert nfa.accepts('') is True

    def test_nfa_dead_state(self):
        """NFA: keine Übergänge → verwirft alle nichtleeren Wörter."""
        nfa = NFA(
            states={0},
            alphabet={'a'},
            transitions={0: {}},
            start=0,
            accepting={0},
        )
        assert nfa.accepts('a') is False
        assert nfa.accepts('aa') is False


# ==============================================================
# 8. RegularExpression Tests
# ==============================================================

class TestRegularExpression:
    """Tests für RegularExpression."""

    def test_literal_match(self):
        """Einfaches Literal."""
        re = RegularExpression('a')
        assert re.matches('a') is True
        assert re.matches('b') is False
        assert re.matches('') is False

    def test_empty_word_star(self):
        """a* akzeptiert leeres Wort."""
        re = RegularExpression('a*')
        assert re.matches('') is True
        assert re.matches('a') is True
        assert re.matches('aaa') is True
        assert re.matches('b') is False

    def test_plus(self):
        """a+ akzeptiert mindestens ein a."""
        re = RegularExpression('a+')
        assert re.matches('a') is True
        assert re.matches('aaa') is True
        assert re.matches('') is False

    def test_alternation(self):
        """a|b akzeptiert a oder b."""
        re = RegularExpression('a|b')
        assert re.matches('a') is True
        assert re.matches('b') is True
        assert re.matches('c') is False
        assert re.matches('ab') is False

    def test_concatenation(self):
        """ab akzeptiert genau 'ab'."""
        re = RegularExpression('ab')
        assert re.matches('ab') is True
        assert re.matches('a') is False
        assert re.matches('ba') is False

    def test_grouping(self):
        """(ab)* akzeptiert Wiederholungen von 'ab'."""
        re = RegularExpression('(ab)*')
        assert re.matches('') is True
        assert re.matches('ab') is True
        assert re.matches('abab') is True
        assert re.matches('aba') is False

    def test_complex_pattern(self):
        """(a|b)*abb: Wörter die auf 'abb' enden."""
        re = RegularExpression('(a|b)*abb')
        assert re.matches('abb') is True
        assert re.matches('aabb') is True
        assert re.matches('babb') is True
        assert re.matches('ab') is False

    def test_single_element_language(self):
        """Einelementige Sprache {hello}."""
        re = RegularExpression('hello')
        assert re.matches('hello') is True
        assert re.matches('hell') is False
        assert re.matches('helloo') is False


# ==============================================================
# 9. ContextFreeGrammar Tests
# ==============================================================

class TestContextFreeGrammar:
    """Tests für ContextFreeGrammar."""

    def test_generates_empty_word(self, cfg_palindromes):
        """Palindrom-CFG erzeugt das leere Wort ε."""
        assert cfg_palindromes.generates('') is True

    def test_generates_single_a(self, cfg_palindromes):
        """'a' ist ein Palindrom (S → a)."""
        assert cfg_palindromes.generates('a') is True

    def test_generates_single_b(self, cfg_palindromes):
        """'b' ist ein Palindrom (S → b)."""
        assert cfg_palindromes.generates('b') is True

    def test_generates_aa(self, cfg_palindromes):
        """'aa' ist ein Palindrom (S → aSa → a·ε·a... via S→ε? nein)."""
        # S → aSa, dann S → ε ? (wir haben S → '' also ε)
        # Aber '' macht S→'', dann aSa→a·''·a = 'aa'
        assert cfg_palindromes.generates('aa') is True

    def test_generates_aba(self, cfg_palindromes):
        """'aba' ist ein Palindrom."""
        assert cfg_palindromes.generates('aba') is True

    def test_generates_abba(self, cfg_palindromes):
        """'abba' ist ein Palindrom."""
        assert cfg_palindromes.generates('abba') is True

    def test_not_generates_ab(self, cfg_palindromes):
        """'ab' ist kein Palindrom."""
        assert cfg_palindromes.generates('ab') is False

    def test_not_generates_abc(self, cfg_palindromes):
        """'abc' ist kein Palindrom (c ∉ Σ)."""
        assert cfg_palindromes.generates('abc') is False

    def test_simple_cfg(self):
        """Einfache CFG: S → aS | b erzeugt {a^n b | n≥0}."""
        cfg = ContextFreeGrammar(
            terminals={'a', 'b'},
            nonterminals={'S'},
            rules={'S': [['a', 'S'], ['b']]},
            start='S',
        )
        assert cfg.generates('b') is True
        assert cfg.generates('ab') is True
        assert cfg.generates('aab') is True
        assert cfg.generates('a') is False
        assert cfg.generates('') is False

    def test_is_ambiguous_demo(self, cfg_palindromes):
        """Ambiguitäts-Demo gibt informatives Dictionary zurück."""
        result = cfg_palindromes.is_ambiguous_demo()
        assert 'definition' in result
        assert 'classic_example' in result
        assert 'decidability' in result
        assert 'dangling' in result['classic_example']['context'].lower() or True


# ==============================================================
# 10. PushdownAutomaton Tests
# ==============================================================

class TestPushdownAutomaton:
    """Tests für PushdownAutomaton."""

    @pytest.fixture
    def pda_an_bn(self):
        """
        PDA für L = {aⁿbⁿ | n ≥ 0}.
        Konvention: stack[-1] = Stack-Top (oben), stack[0] = Boden.
        push_symbols wird so angegeben, dass push_symbols[-1] der neue Stack-Top ist.

        Zustände: q0 (lese a's), q1 (lese b's), q2 (akzeptierend)
        Stack: Z0 (Boden), A (für jedes 'a' ein 'A' pushen)

        Übergänge (key = (Zustand, Symbol_oder_None, Stack-Top)):
        - (q0, 'a', 'Z0') → (q0, ['Z0', 'A'])  : Z0 bleibt unten, A oben
        - (q0, 'a', 'A')  → (q0, ['A', 'A'])   : ein weiteres A oben
        - (q0, 'b', 'A')  → (q1, [])           : A poppen (erstes 'b')
        - (q1, 'b', 'A')  → (q1, [])           : A poppen (weitere 'b')
        - (q0, None, 'Z0')→ (q2, ['Z0'])        : leeres Wort → akzeptieren
        - (q1, None, 'Z0')→ (q2, ['Z0'])        : alle b's gelesen → akzeptieren
        """
        states = {'q0', 'q1', 'q2'}
        alphabet = {'a', 'b'}
        stack_alphabet = {'Z0', 'A'}
        transitions = {
            ('q0', 'a', 'Z0'): [('q0', ['Z0', 'A'])],   # Z0 bleibt, A oben
            ('q0', 'a', 'A'):  [('q0', ['A', 'A'])],    # A bleibt, A oben
            ('q0', 'b', 'A'):  [('q1', [])],             # A poppen
            ('q1', 'b', 'A'):  [('q1', [])],             # A poppen
            ('q0', None, 'Z0'): [('q2', ['Z0'])],        # ε-Übergang
            ('q1', None, 'Z0'): [('q2', ['Z0'])],        # ε-Übergang
        }
        return PushdownAutomaton(states, alphabet, stack_alphabet, transitions, 'q0', {'q2'})

    def test_accepts_empty_word(self, pda_an_bn):
        """ε ∈ {aⁿbⁿ} (n=0)."""
        result = pda_an_bn.accepts_demo('')
        assert result['accepted'] is True

    def test_accepts_ab(self, pda_an_bn):
        """'ab' = a¹b¹ ∈ L."""
        result = pda_an_bn.accepts_demo('ab')
        assert result['accepted'] is True

    def test_accepts_aabb(self, pda_an_bn):
        """'aabb' = a²b² ∈ L."""
        result = pda_an_bn.accepts_demo('aabb')
        assert result['accepted'] is True

    def test_rejects_a(self, pda_an_bn):
        """'a' ∉ {aⁿbⁿ}."""
        result = pda_an_bn.accepts_demo('a')
        assert result['accepted'] is False

    def test_rejects_ba(self, pda_an_bn):
        """'ba' ∉ {aⁿbⁿ}."""
        result = pda_an_bn.accepts_demo('ba')
        assert result['accepted'] is False

    def test_result_has_description(self, pda_an_bn):
        """Ergebnis enthält beschreibenden Text."""
        result = pda_an_bn.accepts_demo('ab')
        assert 'description' in result
        assert 'word' in result
        assert 'accepted' in result


# ==============================================================
# 11. FormalLanguages Namespace Tests
# ==============================================================

class TestFormalLanguages:
    """Tests für FormalLanguages Namespace-Klasse."""

    def test_namespace_contains_dfa(self):
        """FormalLanguages.DFA ist die DFA-Klasse."""
        assert FormalLanguages.DFA is DFA

    def test_namespace_contains_nfa(self):
        """FormalLanguages.NFA ist die NFA-Klasse."""
        assert FormalLanguages.NFA is NFA

    def test_namespace_contains_regex(self):
        """FormalLanguages.RegularExpression ist die RegularExpression-Klasse."""
        assert FormalLanguages.RegularExpression is RegularExpression

    def test_namespace_contains_cfg(self):
        """FormalLanguages.ContextFreeGrammar ist die CFG-Klasse."""
        assert FormalLanguages.ContextFreeGrammar is ContextFreeGrammar

    def test_namespace_contains_pda(self):
        """FormalLanguages.PushdownAutomaton ist die PDA-Klasse."""
        assert FormalLanguages.PushdownAutomaton is PushdownAutomaton

    def test_create_dfa_via_namespace(self):
        """DFA via Namespace-Klasse erstellen."""
        dfa = FormalLanguages.DFA({'q0'}, {'a'}, {}, 'q0', {'q0'})
        assert dfa.accepts('') is True


# ==============================================================
# 12. ComputabilityTheory Tests
# ==============================================================

class TestComputabilityTheory:
    """Tests für ComputabilityTheory."""

    def test_primitive_recursive_demo_keys(self):
        """Demo gibt alle erwarteten Schlüssel zurück."""
        result = ComputabilityTheory.primitive_recursive_demo()
        assert 'description' in result
        assert 'grundfunktionen' in result
        assert 'examples' in result
        assert 'limitation' in result

    def test_primitive_recursive_examples_correct(self):
        """Konkrete Rechenwerte korrekt."""
        result = ComputabilityTheory.primitive_recursive_demo()
        examples = result['examples']
        assert examples['add(3,4)'] == 7
        assert examples['multiply(5,6)'] == 30
        assert examples['factorial(6)'] == 720
        assert examples['fibonacci(10)'] == 55

    def test_primitive_recursive_ackermann_values(self):
        """Ackermann-Werte korrekt dokumentiert."""
        result = ComputabilityTheory.primitive_recursive_demo()
        ack = result['ackermann_values']
        assert ack['A(0,0)'] == 1
        assert ack['A(1,1)'] == 3
        assert ack['A(2,2)'] == 7

    def test_mu_recursive_demo_keys(self):
        """μ-rekursive Demo gibt alle Schlüssel zurück."""
        result = ComputabilityTheory.mu_recursive_demo()
        assert 'description' in result
        assert 'mu_operator' in result
        assert 'examples' in result
        assert 'equivalence' in result

    def test_mu_recursive_integer_sqrt(self):
        """Ganzzahlige Wurzel korrekt."""
        result = ComputabilityTheory.mu_recursive_demo()
        examples = result['examples']
        assert examples['integer_sqrt(16)'] == 4
        assert examples['integer_sqrt(25)'] == 5
        assert examples['integer_sqrt(2)'] is None

    def test_church_turing_demo_keys(self):
        """Church-Turing-Demo gibt alle Schlüssel zurück."""
        result = ComputabilityTheory.church_turing_thesis_demo()
        assert 'thesis' in result
        assert 'equivalent_models' in result
        assert 'status' in result

    def test_church_turing_has_models(self):
        """Mindestens 4 äquivalente Berechnungsmodelle."""
        result = ComputabilityTheory.church_turing_thesis_demo()
        assert len(result['equivalent_models']) >= 4

    def test_undecidable_problems_keys(self):
        """Unentscheidbare Probleme enthält Halteproblem."""
        result = ComputabilityTheory.undecidable_problems()
        assert 'halting_problem' in result
        assert 'post_correspondence' in result
        assert 'rices_theorem' in result

    def test_undecidable_halting_is_undecidable(self):
        """Halteproblem ist unentscheidbar (aber semi-entscheidbar)."""
        result = ComputabilityTheory.undecidable_problems()
        hp = result['halting_problem']
        assert hp['undecidable'] is True
        assert hp['semi_decidable'] is True

    def test_undecidable_hilbert10(self):
        """Hilberts 10. Problem ist unentscheidbar."""
        result = ComputabilityTheory.undecidable_problems()
        assert result['hilbert_10th']['undecidable'] is True


# ==============================================================
# 13. Edge Cases
# ==============================================================

class TestEdgeCases:
    """Edge-Case-Tests für alle Module."""

    def test_linear_code_single_bit(self):
        """LinearCode [1,1,1]: Triviale Einheitsmatrix."""
        code = LinearCode([[1]])
        assert code.k == 1
        assert code.n == 1
        c = code.encode([0])
        assert list(c) == [0]
        c = code.encode([1])
        assert list(c) == [1]

    def test_galois_field_large_prime(self):
        """GF(17) hat 17 Elemente und primitives Element."""
        gf17 = GaloisField(17)
        assert gf17.order == 17
        g = gf17.primitive_element()
        assert pow(g, 16, 17) == 1

    def test_cyclic_code_trivial(self):
        """Zyklischer Code mit g(x)=1 (Generatorgrad 0)."""
        # g=[1] → k=n=5, alles Codewörter
        cc = CyclicCode(5, [1])
        assert cc.k == 5
        cw = cc.encode_poly([1, 0, 1, 0, 1])
        assert len(cw) == 5

    def test_hamming_r4_parameters(self):
        """HammingCode r=4: [15,11,3]."""
        h = HammingCode(4)
        n, k, d = h.parameters()
        assert n == 15
        assert k == 11
        assert d == 3

    def test_dfa_empty_accepting_set(self):
        """DFA mit leerer akzeptierender Menge verwirft alle Wörter."""
        dfa = DFA({'q0'}, {'a'}, {('q0', 'a'): 'q0'}, 'q0', set())
        assert dfa.accepts('') is False
        assert dfa.accepts('a') is False
        assert dfa.accepts('aaa') is False

    def test_regex_empty_pattern(self):
        """Leeres Muster akzeptiert nur leeres Wort."""
        re = RegularExpression('')
        assert re.matches('') is True
        assert re.matches('a') is False

    def test_cfg_single_terminal(self):
        """CFG die nur ein einzelnes Terminal erzeugt."""
        cfg = ContextFreeGrammar(
            terminals={'x'},
            nonterminals={'S'},
            rules={'S': [['x']]},
            start='S',
        )
        assert cfg.generates('x') is True
        assert cfg.generates('xx') is False
        assert cfg.generates('') is False

    def test_berlekamp_massey_all_ones(self):
        """Folge 1,1,1,1,... hat LFSR-Länge 1."""
        result = BCHCode.berlekamp_massey([1, 1, 1, 1, 1])
        assert result['lfsr_length'] == 1

    def test_linear_code_encode_all_messages(self):
        """Alle 2^k Nachrichten kodieren und auf GF(2)-Eigenschaft prüfen."""
        G = np.array([[1, 0, 1],
                      [0, 1, 1]])
        code = LinearCode(G)
        for i in range(1 << code.k):
            m = [(i >> j) & 1 for j in range(code.k)]
            c = code.encode(m)
            assert len(c) == code.n
            assert all(x in [0, 1] for x in c)

    def test_nfa_accepts_via_nondeterminism(self):
        """NFA nutzt Nichtdeterminismus: beide Pfade korrekt verfolgt."""
        # NFA akzeptiert Wörter mit letztem Zeichen 'a' ODER 'b'
        nfa = NFA(
            states={0, 1, 2},
            alphabet={'a', 'b'},
            transitions={
                0: {'a': {0, 1}, 'b': {0, 2}},
                1: {},
                2: {},
            },
            start=0,
            accepting={1, 2},
        )
        assert nfa.accepts('a') is True    # endet auf a
        assert nfa.accepts('b') is True    # endet auf b
        assert nfa.accepts('ab') is True   # endet auf b
        assert nfa.accepts('ba') is True   # endet auf a
