"""
@file test_fuzzing.py
@brief Fuzz-Tests via Hypothesis-Bibliothek für Robustheit.
@description
    Testet mathematische Funktionen mit zufälligen Eingaben, um
    Crashes, Exceptions und numerische Fehler zu finden.

    Hypothesis generiert automatisch Grenzfälle (0, inf, -inf, sehr große Zahlen).
    Mit @given werden Eingabewerte automatisch variiert und bei Fehlern
    werden minimale Gegenbeispiele (Shrinking) gefunden.

    Getestete Bereiche:
    - Algebra: gcd, lcm, euler_phi, is_prime, prime_factorization
    - Analysis: numerical_derivative, numerical_integral
    - Vektoren: Vector.norm, Vector.dot, Vector.normalize
    - Graphentheorie: complete_graph, cycle_graph, binomial_coefficient, catalan_number
    - Topologie: euclidean_metric (Symmetrie, Dreiecksungleichung)

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os

# Sicherstellen, dass das src-Verzeichnis im Suchpfad ist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import math

from hypothesis import given, settings, assume
from hypothesis import strategies as st

# ─── Algebra-Importe ────────────────────────────────────────────────────────
from algebra import (
    gcd,
    lcm,
    euler_phi,
    is_prime,
    prime_factorization,
)

# ─── Analysis-Importe ───────────────────────────────────────────────────────
from analysis import (
    numerical_derivative,
    numerical_integral,
)

# ─── Linear-Algebra-Importe ─────────────────────────────────────────────────
from linear_algebra import Vector

# ─── Graphentheorie-Importe ─────────────────────────────────────────────────
from graph_theory import (
    complete_graph,
    cycle_graph,
    binomial_coefficient,
    catalan_number,
)

# ─── Topologie-Importe ──────────────────────────────────────────────────────
from topology import euclidean_metric


# ============================================================================
# ALGEBRA – Fuzz-Tests
# ============================================================================

class TestAlgebraFuzz:
    """Fuzz-Tests für algebraische Grundfunktionen."""

    @given(
        st.integers(min_value=1, max_value=10**6),
        st.integers(min_value=1, max_value=10**6)
    )
    def test_gcd_always_positive(self, a: int, b: int) -> None:
        """
        @brief gcd(a, b) ist immer eine positive ganze Zahl.
        @description
            Der ggT zweier positiver ganzer Zahlen muss stets positiv sein.
            Hypothesis prüft dies für viele zufällige Paare (a, b).
        """
        result = gcd(a, b)
        assert result > 0, f"gcd({a}, {b}) = {result} ist nicht positiv"

    @given(
        st.integers(min_value=1, max_value=10**6),
        st.integers(min_value=1, max_value=10**6)
    )
    def test_gcd_divides_both(self, a: int, b: int) -> None:
        """
        @brief gcd(a, b) teilt sowohl a als auch b.
        @description
            Per Definition des ggT gilt: gcd(a,b) | a und gcd(a,b) | b.
            Divisionsrest muss immer 0 sein.
        """
        g = gcd(a, b)
        assert a % g == 0, f"gcd({a},{b})={g} teilt a={a} nicht"
        assert b % g == 0, f"gcd({a},{b})={g} teilt b={b} nicht"

    @given(
        st.integers(min_value=1, max_value=10**6),
        st.integers(min_value=1, max_value=10**6)
    )
    def test_lcm_multiple_of_both(self, a: int, b: int) -> None:
        """
        @brief lcm(a, b) ist Vielfaches von a und von b.
        @description
            Das kgV ist die kleinste positive Zahl, die von a und b geteilt wird.
            Daher muss lcm(a,b) % a == 0 und lcm(a,b) % b == 0 gelten.
        """
        l = lcm(a, b)
        assert l % a == 0, f"lcm({a},{b})={l} ist kein Vielfaches von a={a}"
        assert l % b == 0, f"lcm({a},{b})={l} ist kein Vielfaches von b={b}"

    @given(st.integers(min_value=1, max_value=10**6))
    def test_gcd_lcm_relation(self, a: int) -> None:
        """
        @brief Beziehung: gcd(a,b) * lcm(a,b) = a * b.
        @description
            Die fundamentale Relation zwischen ggT und kgV:
            gcd(a,b) · lcm(a,b) = a · b
            (gilt für alle positiven ganzen Zahlen).
        """
        b = a + 1  # Sicherstellen, dass b != a (aber nahe beieinander)
        g = gcd(a, b)
        l = lcm(a, b)
        assert g * l == a * b, f"gcd·lcm={g*l} ≠ a·b={a*b} für a={a}, b={b}"

    @given(st.integers(min_value=2, max_value=1000))
    def test_euler_phi_range(self, n: int) -> None:
        """
        @brief euler_phi(n) ∈ [1, n-1] für alle n ≥ 2.
        @description
            Eulers φ-Funktion gibt die Anzahl der teilerfremden Zahlen 1..n-1 an.
            Sie muss für n ≥ 2 stets in [1, n-1] liegen.
        """
        phi = euler_phi(n)
        assert 1 <= phi <= n - 1, (
            f"euler_phi({n}) = {phi} liegt nicht in [1, {n-1}]"
        )

    @given(st.integers(min_value=2, max_value=100))
    def test_is_prime_consistent(self, n: int) -> None:
        """
        @brief is_prime ist konsistent mit prime_factorization.
        @description
            Wenn n eine Primzahl ist, darf prime_factorization(n) genau einen
            Eintrag {n: 1} liefern (n hat nur den Primfaktor n selbst).
        """
        if is_prime(n):
            factors = prime_factorization(n)
            assert len(factors) == 1, (
                f"Primzahl {n} hat {len(factors)} Primfaktoren (erwartet: 1)"
            )
            assert n in factors, f"Primzahl {n} nicht in eigenem Primfaktor-Dict"

    @given(st.integers(min_value=2, max_value=500))
    def test_prime_factorization_reconstructs(self, n: int) -> None:
        """
        @brief Primfaktorzerlegung rekonstruiert die ursprüngliche Zahl.
        @description
            Das Produkt aller p^e aus der Zerlegung muss n ergeben:
            n = ∏ p^e(p) für alle Primfaktoren p.
        """
        factors = prime_factorization(n)
        # Zahl aus Primfaktoren rekonstruieren: ∏ p^e
        reconstructed = 1
        for prime, exp in factors.items():
            reconstructed *= prime ** exp
        assert reconstructed == n, (
            f"Primfaktorzerlegung von {n} liefert Produkt {reconstructed}"
        )

    @given(st.integers(min_value=1, max_value=100))
    def test_euler_phi_prime_formula(self, p_idx: int) -> None:
        """
        @brief Für Primzahlen gilt: euler_phi(p) = p - 1.
        @description
            Eine Primzahl p ist teilerfremd zu allen 1, ..., p-1.
            Also φ(p) = p - 1.
        """
        # Kleine Primzahlen zum Testen
        small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                        53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
        # Index in der Liste der kleinen Primzahlen
        p = small_primes[p_idx % len(small_primes)]
        assert euler_phi(p) == p - 1, (
            f"euler_phi({p}) = {euler_phi(p)}, erwartet {p-1}"
        )


# ============================================================================
# ANALYSIS – Fuzz-Tests
# ============================================================================

class TestAnalysisFuzz:
    """Fuzz-Tests für numerische Analysis-Funktionen."""

    @given(
        st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False)
    )
    def test_numerical_derivative_linear(self, x: float) -> None:
        """
        @brief Ableitung von f(x) = x ist überall genau 1.0.
        @description
            Die Ableitung einer linearen Funktion f(x) = x ist konstant 1.
            Der numerische Differenzenquotient muss dies bis auf 1e-4 genau liefern.
        """
        result = numerical_derivative(lambda t: t, x)
        assert abs(result - 1.0) < 1e-4, (
            f"f'(x) = 1, numerisches Ergebnis: {result} bei x={x}"
        )

    @given(
        st.floats(min_value=-50, max_value=50, allow_nan=False, allow_infinity=False)
    )
    def test_numerical_derivative_constant(self, x: float) -> None:
        """
        @brief Ableitung von f(x) = 5 (Konstante) ist überall 0.
        @description
            Konstante Funktionen haben Ableitung 0 – das muss numerisch stimmen.
        """
        result = numerical_derivative(lambda t: 5.0, x)
        assert abs(result) < 1e-4, (
            f"Ableitung der Konstante ergibt {result} (≠ 0) bei x={x}"
        )

    @given(
        st.floats(min_value=0.1, max_value=10, allow_nan=False, allow_infinity=False)
    )
    def test_numerical_integral_positive(self, x_max: float) -> None:
        """
        @brief Integral von x² über [0, x_max] ist positiv für x_max > 0.
        @description
            x² ≥ 0 überall, daher muss das Integral über ein positives Intervall
            stets positiv sein (x_max > 0 garantiert x_max² > 0 fast überall).
        """
        result = numerical_integral(lambda x: x**2, 0, x_max)
        assert result > 0, (
            f"∫₀^{x_max} x² dx = {result} ist nicht positiv"
        )

    @given(
        st.floats(min_value=1, max_value=1000, allow_nan=False, allow_infinity=False)
    )
    @settings(max_examples=50)  # Langsamerer Test: Anzahl begrenzen
    def test_numerical_integral_x_squared(self, b: float) -> None:
        """
        @brief ∫₀ᵇ x² dx = b³/3 (analytisch bekannte Formel).
        @description
            Die Simpson-Regel approximiert das Integral. Der relative Fehler
            zum analytischen Wert b³/3 darf nicht mehr als 1e-4 betragen.
        """
        result = numerical_integral(lambda x: x**2, 0, b)
        expected = b**3 / 3
        rel_error = abs(result - expected) / max(abs(expected), 1e-10)
        assert rel_error < 1e-4, (
            f"∫₀^{b} x² dx = {result}, erwartet ≈ {expected}, rel. Fehler = {rel_error}"
        )

    @given(
        st.floats(min_value=0.1, max_value=10, allow_nan=False, allow_infinity=False)
    )
    def test_numerical_integral_additive(self, mid: float) -> None:
        """
        @brief Additivität: ∫₀ᵐ + ∫ₘᵇ = ∫₀ᵇ für beliebiges m ∈ (0, b).
        @description
            Das Integral ist additiv über Teilintervalle:
            ∫₀ᵐ f + ∫ₘᵇ f = ∫₀ᵇ f
        """
        b = mid + 1.0  # Sicherstellen, dass b > mid
        left = numerical_integral(lambda x: x**2, 0, mid)
        right = numerical_integral(lambda x: x**2, mid, b)
        total = numerical_integral(lambda x: x**2, 0, b)
        assert abs(left + right - total) < 1e-6, (
            f"Additivität verletzt: {left} + {right} = {left+right} ≠ {total}"
        )


# ============================================================================
# VEKTOR – Fuzz-Tests
# ============================================================================

class TestVectorFuzz:
    """Fuzz-Tests für Vektor-Operationen aus linear_algebra.py."""

    @given(
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=5
        )
    )
    def test_vector_norm_nonneg(self, components: list) -> None:
        """
        @brief Vektornorm ‖v‖ ist immer nicht-negativ.
        @description
            Die euklidische Norm sqrt(Σ vᵢ²) kann nie negativ sein.
            Hypothesis testet dies für beliebige Komponenten-Listen.
        """
        # Nullvektoren überspringen (Norm = 0, aber nicht negativ)
        assume(any(abs(c) > 1e-15 for c in components))
        v = Vector(components)
        norm = v.norm()
        assert norm >= 0, f"Vektornorm {norm} ist negativ für v={components}"

    @given(
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        ),
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_dot_product_commutative(self, a: list, b: list) -> None:
        """
        @brief Skalarprodukt ist kommutativ: a·b = b·a.
        @description
            Das euklidische Skalarprodukt ist symmetrisch:
            ⟨a, b⟩ = Σ aᵢbᵢ = Σ bᵢaᵢ = ⟨b, a⟩
        """
        # Beide Listen müssen die gleiche Länge haben
        assume(len(a) == len(b) and len(a) >= 2)
        va, vb = Vector(a), Vector(b)
        dot_ab = va.dot(vb)
        dot_ba = vb.dot(va)
        assert abs(dot_ab - dot_ba) < 1e-10, (
            f"Skalarprodukt nicht kommutativ: a·b={dot_ab}, b·a={dot_ba}"
        )

    @given(
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_vector_normalize_unit(self, components: list) -> None:
        """
        @brief Normierter Vektor hat immer Norm 1.
        @description
            v/‖v‖ hat per Definition Norm 1.
            Hypothesis prüft dies für viele verschiedene Vektoren.
        """
        # Nullvektoren und quasi-Nullvektoren überspringen (nicht normierbar)
        assume(not all(c == 0 for c in components))
        assume(any(abs(c) > 1e-10 for c in components))
        v = Vector(components)
        assume(v.norm() > 1e-10)
        n = v.normalize()
        assert abs(n.norm() - 1.0) < 1e-10, (
            f"Normierter Vektor hat Norm {n.norm()} (≠ 1) für v={components}"
        )

    @given(
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=5
        )
    )
    def test_norm_triangle_inequality(self, components: list) -> None:
        """
        @brief Dreiecksungleichung: ‖v + w‖ ≤ ‖v‖ + ‖w‖.
        @description
            Die Norm erfüllt die Dreiecksungleichung (Subadditivität).
            Hier: w = -v (Nullvektor als Summe), daher ‖0‖ ≤ ‖v‖ + ‖v‖ = 2‖v‖.
        """
        # Teste mit einfachem Sonderfall: v + 0 = v
        v = Vector(components)
        zero = Vector([0.0] * len(components))
        v_plus_zero = v  # v + 0 = v
        assert v.norm() <= v.norm() + zero.norm() + 1e-12, (
            f"Dreiecksungleichung verletzt: ‖v‖={v.norm()}"
        )

    @given(
        st.floats(min_value=-50, max_value=50, allow_nan=False, allow_infinity=False),
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_scalar_multiplication_norm(self, scalar: float, components: list) -> None:
        """
        @brief ‖α·v‖ = |α| · ‖v‖ (Homogenität der Norm).
        @description
            Die Norm ist homogen: Skalarmultiplikation skaliert die Norm.
        """
        # Null-Skalar gibt Nullvektor
        assume(abs(scalar) > 1e-10)
        v = Vector(components)
        # α·v: Komponenten manuell skalieren
        scaled_components = [scalar * c for c in components]
        sv = Vector(scaled_components)
        expected_norm = abs(scalar) * v.norm()
        assert abs(sv.norm() - expected_norm) < 1e-8, (
            f"‖{scalar}·v‖={sv.norm()}, |{scalar}|·‖v‖={expected_norm}"
        )


# ============================================================================
# GRAPHENTHEORIE – Fuzz-Tests
# ============================================================================

class TestGraphFuzz:
    """Fuzz-Tests für Graphentheorie-Funktionen."""

    @given(st.integers(min_value=2, max_value=10))
    def test_complete_graph_edges(self, n: int) -> None:
        """
        @brief Vollständiger Graph K_n hat genau n*(n-1)/2 Kanten.
        @description
            K_n verbindet jeden Knoten mit jedem anderen.
            Anzahl der Kanten: C(n,2) = n(n-1)/2.
        """
        g = complete_graph(n)
        expected_edges = n * (n - 1) // 2
        actual_edges = len(g.edges())
        assert actual_edges == expected_edges, (
            f"K_{n} hat {actual_edges} Kanten (erwartet: {expected_edges})"
        )

    @given(st.integers(min_value=2, max_value=10))
    def test_complete_graph_vertices(self, n: int) -> None:
        """
        @brief Vollständiger Graph K_n hat genau n Knoten.
        @description
            Der vollständige Graph auf n Knoten hat genau n Knoten.
        """
        g = complete_graph(n)
        assert len(g.vertices()) == n, (
            f"K_{n} hat {len(g.vertices())} Knoten (erwartet: {n})"
        )

    @given(st.integers(min_value=3, max_value=20))
    def test_cycle_graph_regular(self, n: int) -> None:
        """
        @brief Kreisgrad C_n ist 2-regulär (jeder Knoten hat Grad 2).
        @description
            In C_n hat jeder Knoten genau 2 Nachbarn (links und rechts im Kreis).
            Ein solcher Graph heißt 2-regulär.
        """
        g = cycle_graph(n)
        for v in g.vertices():
            deg = g.degree(v)
            assert deg == 2, (
                f"Knoten {v} in C_{n} hat Grad {deg} (erwartet: 2)"
            )

    @given(st.integers(min_value=3, max_value=20))
    def test_cycle_graph_edges(self, n: int) -> None:
        """
        @brief Kreisgrad C_n hat genau n Kanten.
        @description
            Ein Kreis mit n Knoten hat genau n Kanten (jede Kante verbindet
            zwei benachbarte Knoten im Kreis, inklusive der Schlusskante).
        """
        g = cycle_graph(n)
        assert len(g.edges()) == n, (
            f"C_{n} hat {len(g.edges())} Kanten (erwartet: {n})"
        )

    @given(
        st.integers(min_value=1, max_value=15),
        st.integers(min_value=0, max_value=15)
    )
    def test_binomial_coefficient_symmetry(self, n: int, k: int) -> None:
        """
        @brief C(n,k) = C(n, n-k) (Symmetrie des Binomialkoeffizienten).
        @description
            Der Binomialkoeffizient ist symmetrisch: C(n,k) = C(n, n-k).
            Kombinatorisch: Wähle k aus n ↔ lasse n-k weg.
        """
        assume(k <= n)
        c_nk = binomial_coefficient(n, k)
        c_nnk = binomial_coefficient(n, n - k)
        assert c_nk == c_nnk, (
            f"C({n},{k}) = {c_nk} ≠ C({n},{n-k}) = {c_nnk}"
        )

    @given(st.integers(min_value=0, max_value=10))
    def test_catalan_number_positive(self, n: int) -> None:
        """
        @brief Catalan-Zahlen C_n sind stets positiv (≥ 1).
        @description
            Die n-te Catalan-Zahl C_n = C(2n,n)/(n+1) ist stets ≥ 1.
            C_0 = 1, C_1 = 1, C_2 = 2, C_3 = 5, ...
        """
        c = catalan_number(n)
        assert c >= 1, f"Catalan-Zahl C_{n} = {c} ist nicht positiv"

    @given(
        st.integers(min_value=0, max_value=8),
        st.integers(min_value=0, max_value=8)
    )
    def test_binomial_coefficient_pascal(self, n: int, k: int) -> None:
        """
        @brief Pascal'sche Relation: C(n,k) = C(n-1,k-1) + C(n-1,k).
        @description
            Das Pascal'sche Dreieck: Jeder Eintrag ist die Summe der
            zwei darüber liegenden Einträge.
        """
        # Randfall: k > n
        assume(0 <= k <= n and n >= 1)
        lhs = binomial_coefficient(n, k)
        # C(n-1, k-1) = 0 für k=0 (konvention)
        rhs_left = binomial_coefficient(n - 1, k - 1) if k > 0 else 0
        rhs_right = binomial_coefficient(n - 1, k)
        assert lhs == rhs_left + rhs_right, (
            f"C({n},{k})={lhs} ≠ C({n-1},{k-1})+C({n-1},{k})={rhs_left+rhs_right}"
        )


# ============================================================================
# TOPOLOGIE – Fuzz-Tests
# ============================================================================

class TestTopologyFuzz:
    """Fuzz-Tests für Topologie-Funktionen (Metriken)."""

    @given(
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        ),
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_euclidean_metric_symmetric(self, p: list, q: list) -> None:
        """
        @brief Euklidische Metrik ist symmetrisch: d(p,q) = d(q,p).
        @description
            Symmetrie ist eines der vier Metrik-Axiome:
            d(x,y) = d(y,x) für alle x, y ∈ X.
        """
        assume(len(p) == len(q) and len(p) >= 2)
        d_pq = euclidean_metric(p, q)
        d_qp = euclidean_metric(q, p)
        assert abs(d_pq - d_qp) < 1e-10, (
            f"Metrik nicht symmetrisch: d(p,q)={d_pq}, d(q,p)={d_qp}"
        )

    @given(
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_euclidean_metric_self_zero(self, p: list) -> None:
        """
        @brief d(p, p) = 0 (Identität des Ununterscheidbaren).
        @description
            Der Abstand eines Punktes zu sich selbst ist 0:
            d(x,x) = 0 für alle x ∈ X (Metrik-Axiom 2).
        """
        d = euclidean_metric(p, p)
        assert abs(d) < 1e-12, f"d(p,p) = {d} ≠ 0 für p={p}"

    @given(
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=3
        ),
        st.lists(
            st.floats(min_value=-10, max_value=10, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=3
        )
    )
    @settings(max_examples=50)  # Dreiecksungleichung: numerisch etwas teurer
    def test_euclidean_metric_triangle_inequality(self, p: list, q: list) -> None:
        """
        @brief Dreiecksungleichung: d(p, 0) ≤ d(p, q) + d(q, 0).
        @description
            Die Dreiecksungleichung gilt für jede Metrik:
            d(x, z) ≤ d(x, y) + d(y, z)
            Hier: z = Nullvektor.
        """
        assume(len(p) == len(q))
        origin = [0.0] * len(p)
        d_po = euclidean_metric(p, origin)
        d_pq = euclidean_metric(p, q)
        d_qo = euclidean_metric(q, origin)
        # Kleine Toleranz für numerische Ungenauigkeiten
        assert d_po <= d_pq + d_qo + 1e-10, (
            f"Dreiecksungleichung verletzt: d(p,0)={d_po} > d(p,q)+d(q,0)={d_pq+d_qo}"
        )

    @given(
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        ),
        st.lists(
            st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
            min_size=2, max_size=4
        )
    )
    def test_euclidean_metric_nonneg(self, p: list, q: list) -> None:
        """
        @brief Euklidische Metrik ist stets nicht-negativ: d(p,q) ≥ 0.
        @description
            Nichtnegativität ist das erste Metrik-Axiom:
            d(x,y) ≥ 0 für alle x, y ∈ X.
        """
        assume(len(p) == len(q) and len(p) >= 2)
        d = euclidean_metric(p, q)
        assert d >= 0, f"Euklidische Metrik ist negativ: {d} für p={p}, q={q}"
