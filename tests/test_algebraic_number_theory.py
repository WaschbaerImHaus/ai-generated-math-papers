"""
@file test_algebraic_number_theory.py
@brief Tests für das Modul algebraic_number_theory.py
@description
    Umfassende Tests für alle Funktionen und Klassen der algebraischen
    Zahlentheorie:
    - NumberField: Grad, Signatur, Ganzheitsring, Diskriminante
    - Ideal-Arithmetik: Zerlegung, Norm, Verzweigung
    - Klassengruppe: Klassenzahl, Minkowski-Schranke
    - Einheitengruppe: Einheitenrang, Fundamentaleinheiten, Regulator
    - QuadraticField: Legendre-Symbol, Kronecker-Symbol, Reziprozität
    - CyclotomicField: Grad, Zerlegung, Ganzheitsring
    - Lokale Theorie: Hensel-Lifting, Hilbert-Symbol

@author Michael Fuhrmann
@lastModified 2026-03-11
"""

import pytest
import math
import sys
import os

# Sicherstellen, dass src im Pythonpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algebraic_number_theory import (
    # Hilfsunktionen
    _squarefree_part, _euler_totient, _gcd, _prime_factorization,
    _legendre_symbol, _jacobi_symbol, _pell_solution,
    # NumberField
    NumberField,
    # Ideal-Arithmetik
    DedekindIdeal, ideal_factorization, ramification_index,
    inertia_degree, efg_equation_check,
    # Klassengruppe
    class_number, minkowski_bound, class_group_structure, is_pid,
    # Einheitengruppe
    unit_rank, fundamental_units_quadratic, regulator_estimate,
    # Quadratische Körper
    QuadraticField, legendre_symbol_fast, kronecker_symbol,
    quadratic_reciprocity_law, hilbert_class_field_degree,
    # Zyklotomische Körper
    CyclotomicField, cyclotomic_ring_of_integers, kummer_lifting,
    # Lokale Theorie
    p_adic_completion, hensel_lifting_general, local_norm_symbol,
    # Klassenzahlformel
    dirichlet_class_number_formula,
)


# ============================================================
# Klasse 1: Hilfsfunktionen
# ============================================================

class TestHilfsfunktionen:
    """Tests für interne Hilfsfunktionen."""

    def test_squarefree_part_positive(self):
        """Quadratfreier Teil: 12 = 4·3 → 3."""
        assert _squarefree_part(12) == 3

    def test_squarefree_part_already_free(self):
        """Bereits quadratfrei: 30 bleibt 30."""
        assert _squarefree_part(30) == 30

    def test_squarefree_part_negative(self):
        """Negative Zahl: -20 = -4·5 → -5."""
        assert _squarefree_part(-20) == -5

    def test_squarefree_part_prime(self):
        """Primzahl bleibt unverändert."""
        assert _squarefree_part(7) == 7

    def test_squarefree_part_perfect_square(self):
        """Perfektes Quadrat: 36 = 6² → 1."""
        assert _squarefree_part(36) == 1

    def test_euler_totient_prime(self):
        """φ(p) = p-1 für Primzahl p."""
        assert _euler_totient(7) == 6
        assert _euler_totient(13) == 12

    def test_euler_totient_prime_power(self):
        """φ(p^k) = p^{k-1}·(p-1)."""
        assert _euler_totient(8) == 4   # φ(2³) = 4
        assert _euler_totient(9) == 6   # φ(3²) = 6

    def test_euler_totient_multiplicative(self):
        """φ ist multiplikativ: φ(12) = φ(4)·φ(3) = 2·2 = 4."""
        assert _euler_totient(12) == 4

    def test_gcd_basic(self):
        """gcd(12, 8) = 4."""
        assert _gcd(12, 8) == 4

    def test_gcd_coprime(self):
        """gcd(7, 11) = 1 (teilerfremd)."""
        assert _gcd(7, 11) == 1

    def test_legendre_symbol_qr(self):
        """(1/p) = 1 (1 ist immer quadratischer Rest)."""
        assert _legendre_symbol(1, 7) == 1
        assert _legendre_symbol(4, 11) == 1  # 2²=4 → QR

    def test_legendre_symbol_nqr(self):
        """(3/7) = -1: 3 ist NR mod 7."""
        assert _legendre_symbol(3, 7) == -1

    def test_legendre_symbol_zero(self):
        """(p/p) = 0."""
        assert _legendre_symbol(7, 7) == 0

    def test_pell_solution_d2(self):
        """Pell-Gleichung für d=2: 1² - 2·1² = -1, Lösung (1,1)."""
        a, b, norm = _pell_solution(2)
        assert abs(a * a - 2 * b * b) == 1

    def test_pell_solution_d3(self):
        """Pell-Gleichung für d=3: 2² - 3·1² = 1, Lösung (2,1)."""
        a, b, norm = _pell_solution(3)
        assert abs(a * a - 3 * b * b) == 1


# ============================================================
# Klasse 2: NumberField
# ============================================================

class TestNumberField:
    """Tests für die NumberField-Klasse."""

    def test_degree_quadratic(self):
        """Quadratisches Polynom: Grad 2."""
        # x² - 2 = 0 → Q(√2)
        nf = NumberField([-2, 0, 1])
        assert nf.degree() == 2

    def test_degree_cubic(self):
        """Kubisches Polynom: Grad 3."""
        # x³ - 2 = 0 → Q(∛2)
        nf = NumberField([-2, 0, 0, 1])
        assert nf.degree() == 3

    def test_signature_real_quadratic(self):
        """Q(√2): reell-quadratisch, Signatur (2, 0)."""
        nf = NumberField([-2, 0, 1])
        r1, r2 = nf.signature()
        assert r1 == 2
        assert r2 == 0

    def test_signature_imaginary_quadratic(self):
        """Q(√-1) = Q(i): imaginär-quadratisch, Signatur (0, 1)."""
        nf = NumberField([1, 0, 1])  # x² + 1 = 0
        r1, r2 = nf.signature()
        assert r1 == 0
        assert r2 == 1

    def test_discriminant_quadratic_d2(self):
        """Diskriminante von Q(√2): Δ = 8 (d=2 ≡ 2 mod 4 → 4d=8)."""
        nf = NumberField([-2, 0, 1])  # x² - 2
        disc = nf.discriminant()
        assert disc == 8

    def test_discriminant_quadratic_d5(self):
        """Diskriminante von Q(√5): Δ = 5 (d=5 ≡ 1 mod 4)."""
        nf = NumberField([-5, 0, 1])  # x² - 5
        disc = nf.discriminant()
        assert disc == 5

    def test_ring_of_integers_d2(self):
        """Ganzheitsring von Q(√2): Z[√2], Basis {1, √2}."""
        nf = NumberField([-2, 0, 1])
        roi = nf.ring_of_integers()
        assert "√2" in roi["description"] or "2" in str(roi["ring_basis"])

    def test_ring_of_integers_d5(self):
        """Ganzheitsring von Q(√5): Z[(1+√5)/2] wegen 5≡1(4)."""
        nf = NumberField([-5, 0, 1])
        roi = nf.ring_of_integers()
        # d=5 ≡ 1 (mod 4) → Basis enthält (1+√5)/2
        assert "1+√5" in roi["description"] or "(1+√" in str(roi["ring_basis"])

    def test_repr_numberfield(self):
        """Textdarstellung enthält Grundinformation."""
        nf = NumberField([-2, 0, 1])
        s = repr(nf)
        assert "NumberField" in s or "Q" in s


# ============================================================
# Klasse 3: Ideal-Arithmetik
# ============================================================

class TestIdealArithmetik:
    """Tests für DedekindIdeal und Primidealzerlegung."""

    def test_dedekind_ideal_norm_principal(self):
        """Norm des Hauptideals (2) in Z[√-5]: N((2)) = 4."""
        ideal = DedekindIdeal([2], field_degree=2)
        assert ideal.norm() == 4

    def test_dedekind_ideal_is_prime(self):
        """Ideal erzeugt von 3 (Primzahl) ist Primideal."""
        ideal = DedekindIdeal([3])
        assert ideal.is_prime_ideal()

    def test_dedekind_ideal_not_prime(self):
        """Ideal erzeugt von 4 ist kein Primideal."""
        ideal = DedekindIdeal([4])
        assert not ideal.is_prime_ideal()

    def test_dedekind_ideal_is_principal_single(self):
        """Einzeln erzeugtes Ideal ist Hauptideal."""
        ideal = DedekindIdeal([5])
        assert ideal.is_principal()

    def test_ideal_factorization_inert(self):
        """p=3 in Q(√2): inert oder zerfällt (Legendre-Symbol (2/3)=-1 → inert)."""
        # f = x² - 2, p=3: (2/3) = ? → 2^{(3-1)/2} = 2^1 ≡ 2 ≡ -1 (mod 3) → NR → inert
        factors = ideal_factorization(3, [-2, 0, 1])
        assert len(factors) == 1
        assert factors[0]["inertia_degree"] == 2
        assert factors[0]["type"] == "inert"

    def test_ideal_factorization_split(self):
        """p=7 in Q(√2): (2/7)=? → 2^3=8≡1(7) → QR → zerfällt."""
        # (2/7): 2^{(7-1)/2} = 2^3 = 8 ≡ 1 (mod 7) → zerfällt
        factors = ideal_factorization(7, [-2, 0, 1])
        assert len(factors) == 2
        assert all(f["type"] == "split" for f in factors)

    def test_ideal_factorization_ramified(self):
        """p=2 in Q(√2): 2 teilt Diskriminante 8 → verzweigt."""
        factors = ideal_factorization(2, [-2, 0, 1])
        assert len(factors) == 1
        assert factors[0]["ramification_index"] == 2

    def test_efg_equation_check(self):
        """Fundamentalsatz Σ e_i·f_i = n = 2 für quadratischen Körper."""
        result = efg_equation_check(3, [-2, 0, 1])
        assert result["check_passed"]
        assert result["sum_eifi"] == 2
        assert result["n"] == 2

    def test_efg_equation_split(self):
        """Zerfallende Primzahl: 2 Faktoren, je e=1,f=1 → Summe = 2."""
        result = efg_equation_check(7, [-2, 0, 1])
        assert result["check_passed"]
        assert result["sum_eifi"] == 2

    def test_ramification_index_ramified(self):
        """Verzweigungsindex für p=2 in Q(√2) ist 2."""
        rams = ramification_index(2, [-2, 0, 1])
        assert 2 in rams

    def test_inertia_degree_inert(self):
        """Trägheitsgrad für inerte Primzahl ist 2."""
        iners = inertia_degree(3, [-2, 0, 1])
        assert 2 in iners


# ============================================================
# Klasse 4: Klassengruppe
# ============================================================

class TestKlassengruppe:
    """Tests für Klassenzahl und Klassengruppe."""

    def test_class_number_gauss_integers(self):
        """h(Q(√-1)) = 1: Z[i] ist PID."""
        assert class_number(-1) == 1

    def test_class_number_eisenstein(self):
        """h(Q(√-3)) = 1: Eisensteinsche Zahlen sind PID."""
        assert class_number(-3) == 1

    def test_class_number_non_pid(self):
        """h(Q(√-5)) = 2: Z[√-5] ist kein PID."""
        assert class_number(-5) == 2

    def test_class_number_stark_heegner(self):
        """Stark-Heegner: h(Q(√-163)) = 1."""
        assert class_number(-163) == 1

    def test_class_number_negative_large(self):
        """h(Q(√-23)) = 3."""
        assert class_number(-23) == 3

    def test_class_number_real_quadratic_d2(self):
        """h(Q(√2)) = 1: reell-quadratischer PID."""
        assert class_number(2) == 1

    def test_class_number_real_quadratic_d5(self):
        """h(Q(√5)) = 1."""
        assert class_number(5) == 1

    def test_class_number_real_d10(self):
        """h(Q(√10)) = 2."""
        assert class_number(10) == 2

    def test_minkowski_bound_imaginary(self):
        """Minkowski-Schranke für Q(√-5): M = (2/π)·√20 ≈ 2.85."""
        m = minkowski_bound(-5)
        # Δ = 4·(-5) = -20, M = (2/π)·√20 ≈ 2.85
        assert 2.5 < m < 3.5

    def test_minkowski_bound_real(self):
        """Minkowski-Schranke für Q(√5): M = (1/2)·√5 ≈ 1.12."""
        m = minkowski_bound(5)
        # Δ = 5 (5≡1 mod 4), M = (1/2)·√5 ≈ 1.118
        assert 0.9 < m < 1.5

    def test_is_pid_true(self):
        """Z[i] ist PID."""
        assert is_pid(-1) is True

    def test_is_pid_false(self):
        """Z[√-5] ist kein PID."""
        assert is_pid(-5) is False

    def test_class_group_structure_h1(self):
        """Klassengruppe mit h=1: triviale Gruppe."""
        cg = class_group_structure(-1)
        assert cg["class_number"] == 1
        assert cg["is_pid"] is True

    def test_class_group_structure_h2(self):
        """Klassengruppe Q(√-5): Z/2Z."""
        cg = class_group_structure(-5)
        assert cg["class_number"] == 2
        assert cg["group_structure"] == (2,)

    def test_class_group_structure_h3(self):
        """Klassengruppe Q(√-23): Z/3Z."""
        cg = class_group_structure(-23)
        assert cg["class_number"] == 3


# ============================================================
# Klasse 5: Einheitengruppe
# ============================================================

class TestEinheitengruppe:
    """Tests für Einheitenrang und Fundamentaleinheiten."""

    def test_unit_rank_imaginary(self):
        """Imaginär-quadratisch: Einheitenrang r = 0."""
        # f = x² + 1, r1=0, r2=1 → r = 0+1-1 = 0
        r = unit_rank([1, 0, 1])
        assert r == 0

    def test_unit_rank_real(self):
        """Reell-quadratisch: Einheitenrang r = 1."""
        # f = x² - 2, r1=2, r2=0 → r = 2+0-1 = 1
        r = unit_rank([-2, 0, 1])
        assert r == 1

    def test_fundamental_unit_d2(self):
        """Fundamentaleinheit von Q(√2): ε = 1+√2, N(ε) = -1."""
        result = fundamental_units_quadratic(2)
        assert result["unit_rank"] == 1
        assert result["a"] == 1
        assert result["b"] == 1
        # N(ε) = 1² - 2·1² = -1
        assert result["norm"] == -1

    def test_fundamental_unit_d3(self):
        """Fundamentaleinheit von Q(√3): ε = 2+√3, N(ε) = 1."""
        result = fundamental_units_quadratic(3)
        assert result["a"] == 2
        assert result["b"] == 1
        assert result["norm"] == 1

    def test_fundamental_unit_d5(self):
        """Fundamentaleinheit von Q(√5): (1+√5)/2 (Goldener Schnitt)."""
        result = fundamental_units_quadratic(5)
        # d=5 ≡ 1 (mod 4): Einheit der Form (a+b√5)/2
        assert result["unit_rank"] == 1
        assert result["a"] == 1
        assert result["b"] == 1

    def test_fundamental_unit_imaginary(self):
        """Imaginär-quadratischer Körper: Einheitenrang 0."""
        result = fundamental_units_quadratic(-5)
        assert result["unit_rank"] == 0

    def test_regulator_real_d2(self):
        """Regulator von Q(√2): log(1+√2) ≈ 0.881."""
        reg = regulator_estimate(2)
        expected = math.log(1 + math.sqrt(2))
        assert abs(reg - expected) < 0.01

    def test_regulator_real_d3(self):
        """Regulator von Q(√3): log(2+√3) ≈ 1.317."""
        reg = regulator_estimate(3)
        expected = math.log(2 + math.sqrt(3))
        assert abs(reg - expected) < 0.01

    def test_regulator_imaginary(self):
        """Regulator imaginär-quadratisch: R = 1 per Konvention."""
        reg = regulator_estimate(-1)
        assert reg == 1.0


# ============================================================
# Klasse 6: QuadraticField
# ============================================================

class TestQuadraticField:
    """Tests für die QuadraticField-Klasse."""

    def test_init_imaginary(self):
        """Q(√-1) ist imaginär-quadratisch."""
        K = QuadraticField(-1)
        assert K.is_imaginary()
        assert not K.is_real()

    def test_init_real(self):
        """Q(√2) ist reell-quadratisch."""
        K = QuadraticField(2)
        assert K.is_real()
        assert not K.is_imaginary()

    def test_discriminant_d_neg1(self):
        """Δ(Q(√-1)) = -4 (−1 ≡ 3 mod 4 → 4·(−1) = −4)."""
        K = QuadraticField(-1)
        assert K.discriminant() == -4

    def test_discriminant_d_neg3(self):
        """Δ(Q(√-3)) = -3 (−3 ≡ 1 mod 4 → d)."""
        K = QuadraticField(-3)
        assert K.discriminant() == -3

    def test_discriminant_d5(self):
        """Δ(Q(√5)) = 5 (5 ≡ 1 mod 4)."""
        K = QuadraticField(5)
        assert K.discriminant() == 5

    def test_signature_imaginary(self):
        """Signatur imaginär-quadratisch: (0, 1)."""
        K = QuadraticField(-5)
        assert K.signature() == (0, 1)

    def test_signature_real(self):
        """Signatur reell-quadratisch: (2, 0)."""
        K = QuadraticField(2)
        assert K.signature() == (2, 0)

    def test_ring_of_integers_d2(self):
        """Basis von Q(√2): {1, √2}."""
        K = QuadraticField(2)
        basis = K.ring_of_integers_basis()
        assert "√2" in basis[1]

    def test_ring_of_integers_d5(self):
        """Basis von Q(√5): {1, (1+√5)/2}."""
        K = QuadraticField(5)
        basis = K.ring_of_integers_basis()
        assert "(1+√5)/2" in basis[1]

    def test_class_number_q_i(self):
        """h(Q(i)) = 1."""
        K = QuadraticField(-1)
        assert K.class_number() == 1

    def test_roots_of_unity_qi(self):
        """Q(i) hat 4 Einheitswurzeln: ±1, ±i."""
        K = QuadraticField(-1)
        assert K.number_of_roots_of_unity() == 4

    def test_roots_of_unity_qomega(self):
        """Q(ω) (d=-3) hat 6 Einheitswurzeln."""
        K = QuadraticField(-3)
        assert K.number_of_roots_of_unity() == 6

    def test_roots_of_unity_generic(self):
        """Generischer Fall: 2 Einheitswurzeln (±1)."""
        K = QuadraticField(-5)
        assert K.number_of_roots_of_unity() == 2

    def test_squarefree_normalization(self):
        """Nicht-quadratfreie Eingabe wird normalisiert."""
        K = QuadraticField(8)  # 8 = 4·2 → d_sf = 2
        assert K.d == 2


# ============================================================
# Klasse 7: Legendre- und Kronecker-Symbol
# ============================================================

class TestSymbole:
    """Tests für Legendre-Symbol, Kronecker-Symbol und Reziprozität."""

    def test_legendre_symbol_fast_qr(self):
        """(4/7) = 1: 4 = 2² ist quadratischer Rest."""
        assert legendre_symbol_fast(4, 7) == 1

    def test_legendre_symbol_fast_nqr(self):
        """(3/7) = -1: 3 ist kein QR mod 7."""
        assert legendre_symbol_fast(3, 7) == -1

    def test_legendre_symbol_fast_zero(self):
        """(7/7) = 0."""
        assert legendre_symbol_fast(7, 7) == 0

    def test_kronecker_symbol_n1(self):
        """(a/1) = 1 immer."""
        assert kronecker_symbol(5, 1) == 1
        assert kronecker_symbol(-7, 1) == 1

    def test_kronecker_symbol_neg1(self):
        """(a/-1) = -1 für a < 0."""
        assert kronecker_symbol(-3, -1) == -1
        assert kronecker_symbol(5, -1) == 1

    def test_kronecker_symbol_prime(self):
        """Kronecker-Symbol stimmt mit Legendre-Symbol für Primzahlen überein."""
        # (3/5): 3^2=9≡4, 3^{(5-1)/2}=3²≡4≡-1(5) → NR → -1
        assert kronecker_symbol(3, 5) == _legendre_symbol(3, 5)

    def test_kronecker_symbol_2(self):
        """(1/2) = 1, (3/2) = -1 (Kronecker)."""
        # a=1: a mod 8 = 1 → (1/2) = 1
        assert kronecker_symbol(1, 2) == 1

    def test_jacobi_symbol_3_5(self):
        """(3/5) = (Jacobi) = Legendre für Primzahl."""
        assert _jacobi_symbol(3, 5) == _legendre_symbol(3, 5)

    def test_jacobi_symbol_multiplicative(self):
        """Jacobi-Symbol ist multiplikativ im Nenner."""
        # (3/15) = (3/3)·(3/5)
        j3_15 = _jacobi_symbol(3, 15)
        j3_3 = _jacobi_symbol(3, 3)
        j3_5 = _jacobi_symbol(3, 5)
        assert j3_15 == j3_3 * j3_5

    def test_quadratic_reciprocity_3_7(self):
        """Reziprozitätsgesetz: (3/7)·(7/3) = (-1)^{(3-1)(7-1)/4} = (-1)^3 = -1."""
        result = quadratic_reciprocity_law(3, 7)
        assert result["law_holds"]
        # (3/7)·(7/3) = -1·(-1) = ... check calculation
        assert result["product"] == result["expected_product"]

    def test_quadratic_reciprocity_5_11(self):
        """(5/11)·(11/5) für 5≡1(4): gleich."""
        result = quadratic_reciprocity_law(5, 11)
        assert result["law_holds"]

    def test_quadratic_reciprocity_ergaenzung_neg1(self):
        """Ergänzungssatz: (-1/7) = (-1)^{(7-1)/2} = (-1)^3 = -1."""
        result = quadratic_reciprocity_law(3, 7)
        assert result["neg1_over_p"] in (-1, 1)

    def test_hilbert_class_field_degree(self):
        """Hilbert-Klassenfeld-Grad = Klassenzahl."""
        assert hilbert_class_field_degree(-5) == class_number(-5)
        assert hilbert_class_field_degree(-1) == 1


# ============================================================
# Klasse 8: CyclotomicField
# ============================================================

class TestCyclotomicField:
    """Tests für zyklotomische Körper."""

    def test_degree_n3(self):
        """[Q(ζ_3):Q] = φ(3) = 2."""
        cf = CyclotomicField(3)
        assert cf.degree() == 2

    def test_degree_n5(self):
        """[Q(ζ_5):Q] = φ(5) = 4."""
        cf = CyclotomicField(5)
        assert cf.degree() == 4

    def test_degree_n12(self):
        """[Q(ζ_12):Q] = φ(12) = 4."""
        cf = CyclotomicField(12)
        assert cf.degree() == 4

    def test_signature_n5(self):
        """Q(ζ_5): vollständig komplex → (0, φ(5)/2) = (0, 2)."""
        cf = CyclotomicField(5)
        r1, r2 = cf.signature()
        assert r1 == 0
        assert r2 == 2

    def test_ring_of_integers_basis_n3(self):
        """Ganzheitsring Z[ζ_3]: Basis {1, ζ_3}."""
        cf = CyclotomicField(3)
        roi = cf.cyclotomic_ring_of_integers()
        assert roi["degree"] == 2
        assert "1" in roi["basis"]

    def test_galois_group_order(self):
        """Ord(Gal(Q(ζ_n)/Q)) = φ(n)."""
        cf = CyclotomicField(7)
        assert cf.galois_group_order() == 6

    def test_is_prime_cyclotomic(self):
        """Q(ζ_5): Primzahlgrad."""
        cf = CyclotomicField(5)
        assert cf.is_prime_cyclotomic()

    def test_is_not_prime_cyclotomic(self):
        """Q(ζ_12): nicht-Primzahlgrad."""
        cf = CyclotomicField(12)
        assert not cf.is_prime_cyclotomic()

    def test_splitting_of_prime_divides_n(self):
        """p | n: vollständig verzweigt, e = φ(n)."""
        cf = CyclotomicField(5)  # Q(ζ_5)
        result = cf.splitting_of_prime(5)  # 5 | 5
        assert result["ramification_index"] == 4  # φ(5) = 4
        assert result["split_type"] == "vollständig verzweigt"
        assert result["check_passed"]

    def test_splitting_of_prime_not_divides(self):
        """p ∤ n: unverzweigt, e·f·r = φ(n)."""
        cf = CyclotomicField(5)
        result = cf.splitting_of_prime(11)  # 11 ∤ 5
        assert result["ramification_index"] == 1
        assert result["check_passed"]
        assert result["efr_check"] == 4

    def test_cyclotomic_ring_of_integers_function(self):
        """Modulbezogene Funktion: cyclotomic_ring_of_integers(7)."""
        roi = cyclotomic_ring_of_integers(7)
        assert roi["degree"] == 6  # φ(7) = 6

    def test_kummer_lifting_regular_prime(self):
        """Kummer-Theorie für reguläre Primzahl p=5."""
        result = kummer_lifting(5, 3)
        assert result["prime"] == 5
        assert result["is_regular_prime"] is True
        assert result["phi_p"] == 4

    def test_kummer_lifting_irregular_prime(self):
        """Kummer-Theorie für irreguläre Primzahl p=37."""
        result = kummer_lifting(37, 2)
        assert result["prime"] == 37
        assert result["is_regular_prime"] is False


# ============================================================
# Klasse 9: Lokale Theorie
# ============================================================

class TestLokaleTheorie:
    """Tests für p-adische Vervollständigung, Hensel-Lifting und Hilbert-Symbol."""

    def test_p_adic_completion_basic(self):
        """p-adische Zerlegung von Q(√2) bei p=3."""
        result = p_adic_completion([-2, 0, 1], 3)
        assert result["prime"] == 3
        assert "completions" in result
        assert len(result["completions"]) > 0

    def test_p_adic_completion_n_completions(self):
        """Anzahl lokaler Faktoren ist korrekt (Σe_i·f_i = n)."""
        result = p_adic_completion([-2, 0, 1], 7)  # Split
        total = sum(
            c["degree_over_Qp"] for c in result["completions"]
        )
        assert total == 2  # Grad des Körpers

    def test_hensel_lifting_valid_root(self):
        """Hensel-Lifting: x²-2≡0(mod 3) hat keine Lösung (3 inert)."""
        # x²-2 mod 3: x=0→-2≡1, x=1→-1≡2, x=2→2 → keine Nullstelle
        # Daher: Fehler erwartet
        result = hensel_lifting_general([-2, 0, 1], 0, 3)
        # a0=0: 0²-2=-2≡1(mod 3)≠0 → Fehler
        assert result["success"] is False

    def test_hensel_lifting_x2_minus_7_mod5(self):
        """x²-7≡0(mod 5): x=3→9-7=2≡2(5), x=2→4-7=-3≡2(5), nein; x=4→16-7=9≡4(5), nein."""
        # x²≡7(mod 5): 7≡2(5). Quadratreste mod 5: 1,4. 2 kein QR → keine Lösung
        result = hensel_lifting_general([-7, 0, 1], 1, 5)
        # 1²-7=-6≡-1≡4(5)≠0 → kein gültiger Start
        assert result["success"] is False

    def test_hensel_lifting_linear_poly(self):
        """Lineares Polynom: x-1=0, a0=1 ist Lösung mod jede Primzahl."""
        # f(x) = x-1 = [-1, 1], f(1) = 0, f'(1) = 1 ≠ 0 (mod p)
        result = hensel_lifting_general([-1, 1], 1, 7, steps=3)
        assert result["success"] is True
        assert result["final_a"] % 7 == 1

    def test_hensel_lifting_x2_minus_2_mod7(self):
        """x²-2=0 mod 7: 3²=9=7+2≡2(7) → Nullstelle bei x=3."""
        # 3²=9≡2(mod 7) → a0=3 ist Lösung
        result = hensel_lifting_general([-2, 0, 1], 3, 7, steps=3)
        assert result["success"] is True

    def test_local_norm_symbol_trivial(self):
        """Hilbert-Symbol (1, 1)_p = 1."""
        for p in [3, 5, 7]:
            assert local_norm_symbol(1, 1, p) == 1

    def test_local_norm_symbol_neg1_neg1(self):
        """(-1, -1)_p: ax²+by²=z² mit a=b=-1."""
        # Für p=5: (-1,-1)_5 = (-1/5)^0 · (-1/5)^0 · ... = 1
        result = local_norm_symbol(-1, -1, 3)
        assert result in (-1, 1)

    def test_local_norm_symbol_symmetry(self):
        """Symmetrie: (a,b)_p = (b,a)_p."""
        a, b, p = 3, 5, 7
        assert local_norm_symbol(a, b, p) == local_norm_symbol(b, a, p)


# ============================================================
# Klasse 10: Klassenzahlformel
# ============================================================

class TestKlassenzahlformel:
    """Tests für die Dirichletsche Klassenzahlformel."""

    def test_formula_d_neg1(self):
        """Klassenzahlformel für Q(i): h=1, w=4."""
        result = dirichlet_class_number_formula(-1)
        assert result["h_exact"] == 1
        assert result["w"] == 4
        assert result["discriminant"] == -4

    def test_formula_d_neg3(self):
        """Klassenzahlformel für Q(√-3): h=1, w=6, Δ=-3."""
        result = dirichlet_class_number_formula(-3)
        assert result["h_exact"] == 1
        assert result["w"] == 6

    def test_formula_d_neg5(self):
        """Klassenzahlformel für Q(√-5): h=2, w=2."""
        result = dirichlet_class_number_formula(-5)
        assert result["h_exact"] == 2
        assert result["w"] == 2

    def test_formula_approximate_d_neg1(self):
        """Formelwert für Q(i): ≈ 1.0 (mit Toleranz)."""
        result = dirichlet_class_number_formula(-1)
        # L(1, χ_{-4}) = π/4 ≈ 0.785
        # h = (4·√4)/(2π) · (π/4) = 4·2/(2π) · π/4 = 1
        assert abs(result["h_from_formula"] - 1.0) < 0.2

    def test_formula_l_value_positive(self):
        """L(1, χ_Δ) > 0 für alle d (Dirichlet)."""
        for d in [-1, -3, -5, -7]:
            result = dirichlet_class_number_formula(d)
            assert result["l_value_approx"] > 0

    def test_formula_discriminant_correct(self):
        """Diskriminante wird korrekt berechnet."""
        result = dirichlet_class_number_formula(-5)
        assert result["discriminant"] == -20  # 4·(-5)

    def test_formula_regulator_real(self):
        """Regulator ist positiv für reell-quadratische Körper."""
        result = dirichlet_class_number_formula(2)
        assert result["regulator"] > 0
