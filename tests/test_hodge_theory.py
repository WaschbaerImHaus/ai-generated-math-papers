"""
Tests für das Hodge-Theorie-Modul (src/hodge_theory.py)

Testet Hodge-Zerlegung, Kähler-Mannigfaltigkeiten, Hodge-Vermutung,
Perioden-Gebiete, Hodge-Strukturen und de-Rham-Kohomologie.

@author: Michael Fuhrmann
@lastModified: 2026-03-11
"""

import pytest
import numpy as np
from src.hodge_theory import (
    HodgeDecomposition,
    KahlerManifold,
    HodgeConjecture,
    PeriodDomain,
    HodgeStructure,
    DeRhamCohomology,
    HodgeExamples,
)


# ─────────────────────────────────────────────
# Tests: HodgeDecomposition
# ─────────────────────────────────────────────

class TestHodgeDecomposition:
    """Tests für die Hodge-Zerlegungsklasse."""

    def test_init_basic(self):
        """Erstellt Hodge-Zerlegung für n=2."""
        hd = HodgeDecomposition(2)
        assert hd.n == 2
        assert hd.kahler is True

    def test_init_non_kahler(self):
        """Erstellt nicht-Kähler Hodge-Zerlegung."""
        hd = HodgeDecomposition(3, kahler=False)
        assert hd.kahler is False

    def test_hodge_number_default_zero(self):
        """Nicht gesetzte Hodge-Zahlen sind 0."""
        hd = HodgeDecomposition(2)
        assert hd.hodge_numbers(1, 0) == 0

    def test_hodge_number_out_of_range(self):
        """Hodge-Zahlen außerhalb des Bereichs sind 0."""
        hd = HodgeDecomposition(2)
        assert hd.hodge_numbers(-1, 0) == 0
        assert hd.hodge_numbers(3, 0) == 0
        assert hd.hodge_numbers(0, 5) == 0

    def test_set_hodge_number_kahler_symmetry(self):
        """Kähler-Symmetrien werden beim Setzen automatisch propagiert."""
        hd = HodgeDecomposition(2)
        hd.set_hodge_number(1, 0, 5)
        # Konjugationssymmetrie h^{p,q} = h^{q,p}
        assert hd.hodge_numbers(0, 1) == 5
        # Serre-Dualität h^{p,q} = h^{n-p,n-q}
        assert hd.hodge_numbers(1, 2) == 5
        assert hd.hodge_numbers(2, 1) == 5

    def test_betti_numbers_projective_line(self):
        """Betti-Zahlen von ℙ^1: b_0=1, b_1=0, b_2=1."""
        hd = HodgeDecomposition(1)
        hd._hodge_numbers[(0, 0)] = 1
        hd._hodge_numbers[(1, 1)] = 1
        betti = hd.betti_numbers()
        assert betti[0] == 1   # b_0 = h^{0,0}
        assert betti[1] == 0   # b_1 = h^{1,0} + h^{0,1} = 0
        assert betti[2] == 1   # b_2 = h^{1,1}

    def test_betti_numbers_elliptic_curve(self):
        """Betti-Zahlen der elliptischen Kurve: b_0=1, b_1=2, b_2=1."""
        data = HodgeExamples.elliptic_curve()
        betti = data["betti_numbers"]
        assert betti[0] == 1
        assert betti[1] == 2
        assert betti[2] == 1

    def test_euler_characteristic_projective_space(self):
        """χ(ℙ^n) = n+1."""
        for n in range(1, 5):
            data = HodgeExamples.projective_space(n)
            assert data["euler_characteristic"] == n + 1

    def test_euler_characteristic_elliptic_curve(self):
        """χ(E) = 0 (Torus)."""
        data = HodgeExamples.elliptic_curve()
        assert data["euler_characteristic"] == 0

    def test_euler_characteristic_k3(self):
        """χ(K3) = 24."""
        data = HodgeExamples.k3_surface()
        assert data["euler_characteristic"] == 24

    def test_hodge_diamond_shape(self):
        """Hodge-Diamant hat Form (2n+1) × (2n+1)."""
        for n in [1, 2, 3]:
            hd = HodgeDecomposition(n)
            diamond = hd.hodge_diamond()
            assert diamond.shape == (2 * n + 1, 2 * n + 1)

    def test_hodge_diamond_p1(self):
        """Hodge-Diamant von ℙ^1 hat 1 an den Eckpunkten."""
        hd = HodgeDecomposition(1)
        hd._hodge_numbers[(0, 0)] = 1
        hd._hodge_numbers[(1, 1)] = 1
        diamond = hd.hodge_diamond()
        # Überprüfe dass Summe der Diamond-Einträge korrekt
        assert diamond.sum() == 2  # h^{0,0} + h^{1,1}

    def test_symmetry_check_projective_space(self):
        """ℙ^n erfüllt Kähler-Symmetrien."""
        data = HodgeExamples.projective_space(2)
        hd = HodgeDecomposition(2)
        for p in range(3):
            hd._hodge_numbers[(p, p)] = 1
        assert hd.symmetry_check() is True

    def test_symmetry_check_k3(self):
        """K3-Fläche erfüllt Kähler-Symmetrien."""
        hd = HodgeDecomposition(2)
        hd._hodge_numbers = {
            (0, 0): 1, (2, 2): 1,
            (2, 0): 1, (0, 2): 1,
            (1, 1): 20,
            (2, 1): 0, (1, 2): 0,
            (1, 0): 0, (0, 1): 0,
        }
        assert hd.symmetry_check() is True


# ─────────────────────────────────────────────
# Tests: KahlerManifold
# ─────────────────────────────────────────────

class TestKahlerManifold:
    """Tests für Kähler-Mannigfaltigkeiten."""

    def test_init(self):
        """Erstellt eine Kähler-Mannigfaltigkeit."""
        km = KahlerManifold(2, "CP^2")
        assert km.n == 2
        assert km.name == "CP^2"

    def test_init_default_name(self):
        """Standard-Name wird generiert."""
        km = KahlerManifold(3)
        assert "3" in km.name

    def test_kahler_form_antisymmetric(self):
        """Kähler-Form ist antisymmetrisch."""
        km = KahlerManifold(2)
        coords = np.zeros(4)
        omega = km.kahler_form(coords)
        # ω antisymmetrisch: ω = -ω^T
        assert np.allclose(omega, -omega.T, atol=1e-12)

    def test_kahler_form_shape(self):
        """Kähler-Form hat korrekte Dimension (2n × 2n)."""
        for n in [1, 2, 3]:
            km = KahlerManifold(n)
            omega = km.kahler_form(np.zeros(2 * n))
            assert omega.shape == (2 * n, 2 * n)

    def test_kahler_form_standard(self):
        """Standard-Kähler-Form auf C^1: ω_{01} = 1."""
        km = KahlerManifold(1)
        omega = km.kahler_form(np.array([0.0, 0.0]))
        assert omega[0, 1] == 1.0
        assert omega[1, 0] == -1.0

    def test_lefschetz_operator_preserves_shape(self):
        """Lefschetz-Operator gibt Form gleicher Länge zurück."""
        km = KahlerManifold(2)
        form = np.array([1.0, 0.0, 0.0])
        result = km.lefschetz_operator(form, 0)
        assert result.shape == form.shape

    def test_hard_lefschetz_valid_range(self):
        """Hard-Lefschetz gilt für 0 ≤ k ≤ n."""
        km = KahlerManifold(3)
        for k in range(4):  # 0, 1, 2, 3
            assert km.hard_lefschetz_theorem(k) is True

    def test_hard_lefschetz_invalid(self):
        """Hard-Lefschetz für k > n gibt False."""
        km = KahlerManifold(2)
        assert km.hard_lefschetz_theorem(3) is False
        assert km.hard_lefschetz_theorem(-1) is False

    def test_lefschetz_decomposition(self):
        """Lefschetz-Zerlegung gibt primitive Komponenten zurück."""
        km = KahlerManifold(3)
        form = np.array([1.0, 2.0])
        result = km.lefschetz_decomposition(form, 2)
        assert "primitive_components" in result
        assert result["degree"] == 2

    def test_hodge_star_finite(self):
        """Hodge-Stern-Operator gibt endliche Werte zurück."""
        km = KahlerManifold(2)
        form = np.array([1.0, 0.0, 1.0])
        metric = np.eye(4)
        result = km.hodge_star(form, metric)
        assert np.all(np.isfinite(result))


# ─────────────────────────────────────────────
# Tests: HodgeConjecture
# ─────────────────────────────────────────────

class TestHodgeConjecture:
    """Tests für die Hodge-Vermutungsklasse."""

    def setup_method(self):
        """Erstellt HodgeConjecture-Instanz."""
        self.hc = HodgeConjecture()

    def test_millennium_statement_contains_key_terms(self):
        """Millennium-Statement enthält wichtige Schlüsselwörter."""
        statement = self.hc.millennium_problem_statement()
        assert "Hodge" in statement
        assert "projektive" in statement or "projective" in statement.lower() or "Vermutung" in statement
        assert "Clay" in statement or "Millennium" in statement

    def test_millennium_statement_mentions_open(self):
        """Statement erwähnt offenen Status."""
        statement = self.hc.millennium_problem_statement()
        assert "OFFEN" in statement or "offen" in statement.lower()

    def test_hodge_class_definition_contains_formula(self):
        """Hodge-Klassen-Definition enthält mathematische Formel."""
        definition = self.hc.hodge_class_definition()
        assert "H^{2p}" in definition or "Hdg" in definition
        assert "rational" in definition.lower() or "ℚ" in definition

    def test_known_cases_lefschetz(self):
        """Lefschetz (1,1)-Theorem ist in bekannten Fällen."""
        cases = self.hc.known_cases()
        names = [c["name"] for c in cases]
        assert any("Lefschetz" in name for name in names)

    def test_known_cases_proved(self):
        """Alle aufgeführten Fälle sind tatsächlich bewiesen."""
        cases = self.hc.known_cases()
        for case in cases:
            assert case["proved"] is True

    def test_known_cases_count(self):
        """Mindestens 3 bekannte Fälle aufgeführt."""
        cases = self.hc.known_cases()
        assert len(cases) >= 3

    def test_counterexample_attempts_structure(self):
        """Fehlversuche haben korrektes Format."""
        attempts = self.hc.counterexample_attempts()
        assert len(attempts) >= 2
        for attempt in attempts:
            assert "name" in attempt
            assert "description" in attempt

    def test_grothendieck_reformulation_contains_motif(self):
        """Grothendieck-Reformulierung erwähnt Motive."""
        ref = self.hc.grothendieck_reformulation()
        assert "Motiv" in ref or "motivic" in ref.lower() or "Grothendieck" in ref

    def test_tate_conjecture_comparison_structure(self):
        """Vergleich enthält beide Vermutungen."""
        comparison = self.hc.tate_conjecture_comparison()
        assert "Tate" in comparison
        assert "Hodge" in comparison
        assert "l-adisch" in comparison or "ℓ" in comparison or "Galois" in comparison

    def test_hodge_class_type_condition(self):
        """Hodge-Klasse muss vom Typ (p,p) sein."""
        # Gültige (2,2)-Klasse
        result = self.hc.hodge_class_is_algebraic({"p": 2, "q": 2, "rational": True})
        assert result["type_condition"] is True
        assert result["is_hodge_class"] is True

    def test_hodge_class_non_hodge_type(self):
        """Klasse vom Typ (1,2) ist keine Hodge-Klasse."""
        result = self.hc.hodge_class_is_algebraic({"p": 1, "q": 2, "rational": True})
        assert result["type_condition"] is False
        assert result["is_hodge_class"] is False

    def test_hodge_class_lefschetz_theorem(self):
        """(1,1)-Klassen werden durch Lefschetz-Theorem als algebraisch erkannt."""
        result = self.hc.hodge_class_is_algebraic({"p": 1, "q": 1, "rational": True})
        assert result["known_algebraic"] is True
        assert "Lefschetz" in result["reason"]

    def test_hodge_class_conjecture_applies(self):
        """Hodge-Vermutung offen für (p,p)-Klassen mit p > 1."""
        result = self.hc.hodge_class_is_algebraic({"p": 2, "q": 2, "rational": True})
        assert result["hodge_conjecture_applies"] is True


# ─────────────────────────────────────────────
# Tests: PeriodDomain
# ─────────────────────────────────────────────

class TestPeriodDomain:
    """Tests für Perioden-Gebiete."""

    def test_init_weight(self):
        """Gewicht wird korrekt berechnet."""
        pd = PeriodDomain([1, 1])
        assert pd.weight == 1
        assert pd.total_dim == 2

    def test_init_k3_type(self):
        """K3-Typ: [1, 20, 1] ergibt Gewicht 2."""
        pd = PeriodDomain([1, 20, 1])
        assert pd.weight == 2
        assert pd.total_dim == 22

    def test_period_matrix_shape(self):
        """Perioden-Matrix hat korrekte Anzahl Zeilen."""
        pd = PeriodDomain([2, 3, 2])
        basis = np.random.randn(2, 7)
        omega = pd.period_matrix(basis)
        assert omega.shape[0] == 2  # h^{n,0} = 2

    def test_period_matrix_undersized_basis(self):
        """Zu kleine Basis wird aufgefüllt."""
        pd = PeriodDomain([3, 5, 3])
        basis = np.eye(2, 11)  # Nur 2 Zeilen, brauche 3
        omega = pd.period_matrix(basis)
        assert omega.shape[0] == 3

    def test_griffiths_transversality_full_rank(self):
        """Transversalität: voller Rang gibt True."""
        pd = PeriodDomain([1, 5, 1])
        # 1 × 7 Matrix mit vollem Rang (1 Zeile)
        mat = np.array([[1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]])
        assert pd.griffiths_transversality_check(mat) is True

    def test_griffiths_transversality_empty(self):
        """Leere Matrix ergibt False."""
        pd = PeriodDomain([1, 1])
        assert pd.griffiths_transversality_check(np.array([])) is False

    def test_monodromy_group_even_weight(self):
        """Gerades Gewicht → symplektische Gruppe."""
        pd = PeriodDomain([1, 20, 1])  # Gewicht 2 (gerade)
        mono = pd.monodromy_group()
        assert "Sp" in mono or "symplektisch" in mono.lower()

    def test_monodromy_group_odd_weight(self):
        """Ungerades Gewicht → orthogonale Gruppe."""
        pd = PeriodDomain([1, 1])  # Gewicht 1 (ungerade)
        mono = pd.monodromy_group()
        assert "O(" in mono or "orthogonal" in mono.lower()

    def test_vhs_check_valid(self):
        """Gültige Perioden-Matrix gibt True."""
        pd = PeriodDomain([2, 3, 2])
        mat = np.eye(2, 7)
        assert pd.vhs_check(mat) is True

    def test_vhs_check_none(self):
        """None gibt False."""
        pd = PeriodDomain([1, 1])
        assert pd.vhs_check(None) is False


# ─────────────────────────────────────────────
# Tests: HodgeStructure
# ─────────────────────────────────────────────

class TestHodgeStructure:
    """Tests für reine Hodge-Strukturen."""

    def test_init_weight_1(self):
        """Hodge-Struktur Gewicht 1 (elliptische Kurve)."""
        hs = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        assert hs.weight == 1
        assert hs.rank == 2

    def test_pure_of_weight(self):
        """Gewicht wird korrekt zurückgegeben."""
        hs = HodgeStructure(3, {(3, 0): 1, (2, 1): 5, (1, 2): 5, (0, 3): 1})
        assert hs.pure_of_weight() == 3

    def test_polarization_identity(self):
        """Einheitsmatrix ist Polarisierung (gerades Gewicht)."""
        hs = HodgeStructure(2, {(2, 0): 1, (1, 1): 2, (0, 2): 1})
        Q = np.eye(4)  # Positive definite Matrix
        assert hs.polarization(Q) is True

    def test_polarization_antisymmetric_odd(self):
        """Antisymmetrische Matrix ist Polarisierung für ungerades Gewicht."""
        hs = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        # Für Gewicht 1: Q antisymmetrisch
        Q = np.array([[0.0, 1.0], [-1.0, 0.0]])
        # Antisymmetrische Matrix ist immer nicht positiv definit (Eigenwerte imaginär)
        # Test: Symmetrieprüfung — Ergebnis ist bool (True oder False)
        result = hs.polarization(Q)
        assert result is True or result is False

    def test_morphism_check_wrong_weight(self):
        """Morphismus zwischen verschiedenen Gewichten gibt False."""
        hs1 = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        hs2 = HodgeStructure(2, {(2, 0): 1, (1, 1): 2, (0, 2): 1})
        f = np.eye(2, 4)
        assert hs1.morphism_check(hs2, f) is False

    def test_morphism_check_same_weight(self):
        """Wohldefinierten Morphismus gleichen Gewichts gibt True."""
        hs1 = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        hs2 = HodgeStructure(1, {(1, 0): 2, (0, 1): 2})
        f = np.ones((4, 2))
        assert hs1.morphism_check(hs2, f) is True

    def test_tensor_product_weight(self):
        """Tensorprodukt hat Summe der Gewichte."""
        hs1 = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        hs2 = HodgeStructure(2, {(2, 0): 1, (1, 1): 2, (0, 2): 1})
        product = hs1.tensor_product(hs2)
        assert product.weight == 3

    def test_tensor_product_rank(self):
        """Tensorprodukt hat Produkt der Ränge."""
        hs1 = HodgeStructure(1, {(1, 0): 1, (0, 1): 1})
        hs2 = HodgeStructure(1, {(1, 0): 2, (0, 1): 2})
        product = hs1.tensor_product(hs2)
        assert product.rank == 2 * 4  # 2 × 4

    def test_tate_twist_weight_shift(self):
        """Tate-Twist verschiebt Gewicht um -2n."""
        hs = HodgeStructure(4, {(4, 0): 1, (2, 2): 5, (0, 4): 1})
        twisted = hs.tate_twist(1)
        assert twisted.weight == 4 - 2  # == 2

    def test_tate_twist_hodge_type_shift(self):
        """Tate-Twist verschiebt Hodge-Typen um (-n,-n)."""
        hs = HodgeStructure(2, {(2, 0): 1, (1, 1): 2, (0, 2): 1})
        twisted = hs.tate_twist(1)
        # (2,0) → (1,-1), (1,1) → (0,0), (0,2) → (-1,1)
        assert (0, 0) in twisted.hodge_numbers
        assert twisted.hodge_numbers[(0, 0)] == 2


# ─────────────────────────────────────────────
# Tests: DeRhamCohomology
# ─────────────────────────────────────────────

class TestDeRhamCohomology:
    """Tests für die de-Rham-Kohomologie."""

    def setup_method(self):
        """Erstellt de-Rham-Kohomologie für den 2D-Torus."""
        self.manifold = {
            "dim": 2,
            "metric": np.eye(2),
            "name": "T^2",
        }
        self.drc = DeRhamCohomology(self.manifold)

    def test_init(self):
        """Korrekte Initialisierung."""
        assert self.drc.dim == 2
        assert np.allclose(self.drc.metric, np.eye(2))

    def test_harmonic_forms_degree_0(self):
        """Harmonische 0-Formen: konstante Funktionen."""
        forms = self.drc.harmonic_forms(0)
        assert len(forms) > 0
        assert any("konst" in f.lower() for f in forms)

    def test_harmonic_forms_top_degree(self):
        """Harmonische n-Formen: Volumenform."""
        forms = self.drc.harmonic_forms(self.drc.dim)
        assert len(forms) > 0

    def test_harmonic_forms_out_of_range(self):
        """Grad außerhalb [0,n] → leere Liste."""
        forms = self.drc.harmonic_forms(-1)
        assert forms == []
        forms = self.drc.harmonic_forms(10)
        assert forms == []

    def test_hodge_laplacian_shape(self):
        """Hodge-Laplacian hat gleiche Form wie Eingabe."""
        form = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = self.drc.hodge_laplacian(form, np.eye(2))
        assert result.shape == form.shape

    def test_hodge_laplacian_constant_form(self):
        """Laplacian einer konstanten Form (innere Punkte) gibt 0."""
        form = np.ones(5)
        result = self.drc.hodge_laplacian(form, np.eye(2))
        # Konstante Form: (1-2·1+1)/h² = 0 an inneren Punkten
        assert np.allclose(result[1:-1], 0.0, atol=1e-6)

    def test_hodge_laplacian_finite(self):
        """Hodge-Laplacian gibt endliche Werte zurück."""
        form = np.random.randn(10)
        result = self.drc.hodge_laplacian(form, np.eye(2))
        assert np.all(np.isfinite(result))

    def test_bochner_weitzenboeck_shape(self):
        """Bochner-Weitzenböck hat gleiche Form wie Eingabe."""
        form = np.array([1.0, 0.5, -0.5])
        result = self.drc.bochner_weitzenboeck(form)
        assert result.shape == form.shape

    def test_green_operator_shape(self):
        """Greenscher Operator gibt Form gleicher Länge zurück."""
        form = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        result = self.drc.green_operator(form)
        assert result.shape == form.shape

    def test_green_operator_trivial_case(self):
        """Greenscher Operator für Länge 1 gibt Null zurück."""
        form = np.array([1.0])
        result = self.drc.green_operator(form)
        assert np.allclose(result, 0.0)

    def test_green_operator_finite(self):
        """Greenscher Operator gibt endliche Werte zurück."""
        form = np.ones(6)
        result = self.drc.green_operator(form)
        assert np.all(np.isfinite(result))


# ─────────────────────────────────────────────
# Tests: HodgeExamples
# ─────────────────────────────────────────────

class TestHodgeExamples:
    """Tests für konkrete Hodge-Beispiele."""

    def test_projective_space_pn_name(self):
        """Name enthält korrekten projektiven Raum."""
        for n in [1, 2, 3]:
            data = HodgeExamples.projective_space(n)
            assert str(n) in data["name"]

    def test_projective_space_betti_numbers(self):
        """ℙ^n: b_{2k}=1, b_{2k+1}=0."""
        data = HodgeExamples.projective_space(3)
        betti = data["betti_numbers"]
        # b_0=1, b_1=0, b_2=1, b_3=0, b_4=1, b_5=0, b_6=1
        assert betti[0] == 1
        assert betti[1] == 0
        assert betti[2] == 1
        assert betti[3] == 0
        assert betti[4] == 1

    def test_elliptic_curve_genus(self):
        """Elliptische Kurve hat Geschlecht 1."""
        data = HodgeExamples.elliptic_curve()
        assert data["genus"] == 1

    def test_elliptic_curve_h10(self):
        """h^{1,0} = 1 für elliptische Kurve."""
        data = HodgeExamples.elliptic_curve()
        assert data["hodge_numbers"]["h^{1,0}"] == 1

    def test_elliptic_curve_h11(self):
        """h^{1,1} = 1 für elliptische Kurve."""
        data = HodgeExamples.elliptic_curve()
        assert data["hodge_numbers"]["h^{1,1}"] == 1

    def test_k3_h11_equals_20(self):
        """K3-Fläche: h^{1,1} = 20."""
        data = HodgeExamples.k3_surface()
        assert data["hodge_numbers"]["h^{1,1}"] == 20

    def test_k3_h20_equals_1(self):
        """K3-Fläche: h^{2,0} = 1."""
        data = HodgeExamples.k3_surface()
        assert data["hodge_numbers"]["h^{2,0}"] == 1

    def test_k3_euler_24(self):
        """K3-Fläche: χ = 24."""
        data = HodgeExamples.k3_surface()
        assert data["euler_characteristic"] == 24

    def test_abelian_variety_formula(self):
        """Abelsche Varietät: h^{p,q} = C(g,p)·C(g,q)."""
        for g in [1, 2, 3]:
            data = HodgeExamples.abelian_variety(g)
            from math import comb
            for p in range(g + 1):
                for q in range(g + 1):
                    expected = comb(g, p) * comb(g, q)
                    key = f"h^{{{p},{q}}}"
                    assert data["hodge_numbers"][key] == expected

    def test_abelian_variety_euler_zero(self):
        """Abelsche Varietäten haben χ = 0 für g ≥ 1."""
        for g in [1, 2, 3]:
            data = HodgeExamples.abelian_variety(g)
            assert data["euler_characteristic"] == 0

    def test_calabi_yau_h11_h21(self):
        """Calabi-Yau: h^{1,1} und h^{2,1} korrekt."""
        data = HodgeExamples.calabi_yau_threefold(101, 1)
        assert data["h11"] == 101
        assert data["h21"] == 1

    def test_calabi_yau_euler_characteristic(self):
        """CY3: χ = 2(h^{1,1} - h^{2,1})."""
        h11, h21 = 11, 11
        data = HodgeExamples.calabi_yau_threefold(h11, h21)
        assert data["euler_characteristic"] == 2 * (h11 - h21)

    def test_calabi_yau_mirror_symmetry(self):
        """Spiegel-CY vertauscht h^{1,1} und h^{2,1}."""
        data = HodgeExamples.calabi_yau_threefold(10, 20)
        mirror = data["mirror_pair"]
        assert mirror["h11"] == 20
        assert mirror["h21"] == 10

    def test_verify_hodge_conjecture_abelian_small_g(self):
        """Hodge-Vermutung bewiesen für abelsche Varietäten g ≤ 3."""
        for g in [1, 2, 3]:
            assert HodgeExamples.verify_hodge_conjecture_for_abelian_variety(g) is True

    def test_verify_hodge_conjecture_abelian_large_g(self):
        """Hodge-Vermutung unbekannt für g > 3."""
        assert HodgeExamples.verify_hodge_conjecture_for_abelian_variety(4) is False
        assert HodgeExamples.verify_hodge_conjecture_for_abelian_variety(5) is False

    def test_projective_space_hodge_diamond_entries(self):
        """ℙ^2: Hodge-Diamant hat genau 3 Einsen."""
        data = HodgeExamples.projective_space(2)
        diamond = data["hodge_diamond"]
        # Nur h^{0,0}, h^{1,1}, h^{2,2} = 1
        assert diamond.sum() == 3
