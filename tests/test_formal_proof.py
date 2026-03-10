"""
@file test_formal_proof.py
@brief Tests für die formale Beweisinfrastruktur (formal_proof.py).
@description
    Umfassende Tests für alle Klassen und Funktionen der formalen
    Beweisinfrastruktur. Abgedeckte Bereiche:

    - ProofStep: symbolische und numerische Verifikation
    - ProofChain: Aufbau, Zusammenfassung, Textausgabe, HTML-Ausgabe
    - InductionProof: Basisfall, empirische Verifikation, symbolischer Beweis
    - ContradictionProof: Negation, Widerspruch, QED
    - TheoremRegistry: Registrierung, Abfrage, Klassiker laden
    - Vorgefertigte Beweise: Gauß, √2, Primzahlen, Fermat, Riemann, Goldbach

    Alle Tests sind unabhängig voneinander ausführbar.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import pytest
import sympy as sp

# Pfad zum src/-Verzeichnis hinzufügen, damit Importe funktionieren
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from formal_proof import (
    ProofStatus,
    ProofStep,
    ProofChain,
    InductionProof,
    ContradictionProof,
    TheoremRegistry,
    prove_gauss_sum,
    prove_sqrt2_irrational,
    prove_infinitely_many_primes,
    prove_fermat_little_theorem,
    riemann_hypothesis_evidence,
    goldbach_evidence,
)


# ===========================================================================
# TESTS: ProofStatus
# ===========================================================================

class TestProofStatus:
    """
    @brief Tests für den ProofStatus-Enum.
    @description Stellt sicher, dass alle Status-Werte korrekt definiert sind.
    @date 2026-03-10
    """

    def test_all_statuses_exist(self):
        """Alle sechs Status-Werte müssen vorhanden sein."""
        assert ProofStatus.UNVERIFIED.value == "unverified"
        assert ProofStatus.VERIFIED_SYMBOLIC.value == "symbolic"
        assert ProofStatus.VERIFIED_NUMERIC.value == "numeric"
        assert ProofStatus.ASSUMED.value == "assumed"
        assert ProofStatus.FAILED.value == "failed"
        assert ProofStatus.PENDING.value == "pending"

    def test_statuses_are_distinct(self):
        """Alle Status-Werte müssen sich voneinander unterscheiden."""
        values = [s.value for s in ProofStatus]
        # Keine Duplikate: Menge hat gleiche Länge wie Liste
        assert len(values) == len(set(values))


# ===========================================================================
# TESTS: ProofStep
# ===========================================================================

class TestProofStep:
    """
    @brief Tests für die ProofStep-Dataclass.
    @description Testet Erstellung, symbolische und numerische Verifikation.
    @date 2026-03-10
    """

    def test_default_status_is_unverified(self):
        """Ein neuer ProofStep hat standardmäßig UNVERIFIED-Status."""
        step = ProofStep("x > 0", "Annahme")
        assert step.status == ProofStatus.UNVERIFIED

    def test_verify_symbolic_x_squared_minus_x_squared(self):
        """
        Symbolische Verifikation: x² - x² = 0 muss True ergeben.
        Das ist der einfachste symbolische Beweis: beide Seiten sind gleich.
        """
        x = sp.Symbol('x')
        step = ProofStep("x² - x² = 0", "Selbstäquivalenz")
        # lhs = x², rhs = x² → Differenz = 0
        result = step.verify_symbolic(x**2, x**2)
        assert result is True
        assert step.status == ProofStatus.VERIFIED_SYMBOLIC

    def test_verify_symbolic_x_plus_1_not_equal_x(self):
        """
        Symbolische Verifikation: x+1 ≠ x muss False ergeben.
        Die Differenz (x+1) - x = 1 ≠ 0.
        """
        x = sp.Symbol('x')
        step = ProofStep("x+1 = x", "Falsche Aussage")
        result = step.verify_symbolic(x + 1, x)
        assert result is False
        assert step.status == ProofStatus.FAILED

    def test_verify_symbolic_algebraic_identity(self):
        """
        Symbolische Verifikation einer algebraischen Identität:
        (a+b)² = a² + 2ab + b²
        """
        a, b = sp.symbols('a b')
        step = ProofStep("(a+b)² = a² + 2ab + b²", "Binomische Formel")
        lhs = (a + b)**2
        rhs = a**2 + 2*a*b + b**2
        result = step.verify_symbolic(lhs, rhs)
        assert result is True
        assert step.status == ProofStatus.VERIFIED_SYMBOLIC

    def test_verify_symbolic_gauss_induction_step(self):
        """
        Symbolischer Induktionsschritt für Gaußsche Summe:
        n(n+1)/2 + (n+1) = (n+1)(n+2)/2
        """
        n = sp.Symbol('n', positive=True, integer=True)
        step = ProofStep("Induktionsschritt Gauß", "IV + nächster Term")
        lhs = sp.Rational(1, 2) * n * (n + 1) + (n + 1)
        rhs = sp.Rational(1, 2) * (n + 1) * (n + 2)
        result = step.verify_symbolic(lhs, rhs)
        assert result is True

    def test_verify_numeric_correct_values(self):
        """
        Numerische Verifikation mit korrekten Werten muss True ergeben.
        Prüft f(x) = x² für x = 0, 1, 2, 3.
        """
        step = ProofStep("f(x) = x²", "Quadratfunktion")
        result = step.verify_numeric(
            func=lambda x: x * x,
            test_values=[0, 1, 2, 3],
            expected=[0, 1, 4, 9]
        )
        assert result is True
        assert step.status == ProofStatus.VERIFIED_NUMERIC

    def test_verify_numeric_wrong_expected(self):
        """
        Numerische Verifikation mit falschen Erwartungswerten muss False ergeben.
        """
        step = ProofStep("f(x) = x + 1", "Inkrement")
        result = step.verify_numeric(
            func=lambda x: x + 1,
            test_values=[1, 2, 3],
            expected=[2, 3, 99]  # Letzter Wert ist falsch
        )
        assert result is False
        assert step.status == ProofStatus.FAILED

    def test_verify_numeric_mismatched_lengths(self):
        """
        Numerische Verifikation mit unterschiedlich langen Listen muss False ergeben.
        """
        step = ProofStep("Test", "Test")
        result = step.verify_numeric(
            func=lambda x: x,
            test_values=[1, 2, 3],
            expected=[1, 2]  # Zu kurz
        )
        assert result is False

    def test_numeric_checks_populated_after_verify(self):
        """
        Nach verify_numeric soll numeric_checks die geprüften Werte enthalten.
        """
        step = ProofStep("Identität", "f(x) = x")
        step.verify_numeric(lambda x: x, [1, 2, 3], [1, 2, 3])
        assert len(step.numeric_checks) == 3

    def test_str_representation(self):
        """
        __str__ muss Status, Aussage und Begründung enthalten.
        """
        step = ProofStep("Eine Aussage", "Eine Begründung", ProofStatus.ASSUMED)
        text = str(step)
        assert "assumed" in text
        assert "Eine Aussage" in text
        assert "Eine Begründung" in text


# ===========================================================================
# TESTS: ProofChain
# ===========================================================================

class TestProofChain:
    """
    @brief Tests für die ProofChain-Klasse.
    @description Testet Aufbau, Statusberechnung, Textausgabe und HTML.
    @date 2026-03-10
    """

    def test_empty_chain_has_no_steps(self):
        """Eine neue ProofChain hat keine Schritte und ist nicht vollständig."""
        proof = ProofChain("Test-Theorem")
        assert len(proof.steps) == 0
        assert proof.is_complete is False

    def test_assume_adds_step(self):
        """assume() soll einen Schritt mit ASSUMED-Status hinzufügen."""
        proof = ProofChain("Theorem")
        proof.assume("n ∈ ℕ")
        assert len(proof.steps) == 1
        assert proof.steps[0].status == ProofStatus.ASSUMED

    def test_step_adds_step(self):
        """step() ohne symbolische Ausdrücke soll einen UNVERIFIED-Schritt hinzufügen."""
        proof = ProofChain("Theorem")
        proof.step("Schritt 1", "Begründung 1")
        assert len(proof.steps) == 1
        assert proof.steps[0].status == ProofStatus.UNVERIFIED

    def test_step_with_symbolic_verification(self):
        """step() mit lhs und rhs soll symbolisch verifizieren."""
        x = sp.Symbol('x')
        proof = ProofChain("(x+1)² = x² + 2x + 1")
        proof.step(
            "(x+1)² = x² + 2x + 1",
            "Binomische Formel",
            (x + 1)**2,
            x**2 + 2*x + 1
        )
        assert proof.steps[0].status == ProofStatus.VERIFIED_SYMBOLIC

    def test_conclude_sets_is_complete(self):
        """conclude() soll is_complete auf True setzen."""
        proof = ProofChain("Test")
        proof.conclude("Beweis abgeschlossen")
        assert proof.is_complete is True

    def test_chain_has_steps_after_assume_step_conclude(self):
        """
        Nach assume(), step() und conclude() müssen 3 Schritte vorhanden sein.
        """
        proof = ProofChain("Multi-Schritt-Test")
        proof.assume("Voraussetzung A")
        proof.step("Schlussfolgerung B", "Aus A folgt B")
        proof.conclude("Ergebnis C", "q.e.d.")
        assert len(proof.steps) == 3

    def test_method_chaining_works(self):
        """Method-Chaining muss self zurückgeben."""
        proof = ProofChain("Chaining-Test")
        result = proof.assume("X").step("Y", "aus X").conclude("Z")
        # Das Ergebnis der letzten Methode ist die ProofChain
        assert result is proof

    def test_status_summary_counts_correctly(self):
        """status_summary() muss alle Zähler korrekt berechnen."""
        x = sp.Symbol('x')
        proof = ProofChain("Zähl-Test")
        proof.assume("A")                                         # ASSUMED
        proof.step("B", "weil", x, x)                            # VERIFIED_SYMBOLIC
        proof.step("C", "unbekannt")                              # UNVERIFIED
        proof.conclude("D")                                       # ASSUMED

        summary = proof.status_summary()
        assert summary['total_steps'] == 4
        assert summary['assumed'] == 2       # assume + conclude
        assert summary['symbolic'] == 1
        assert summary['unverified'] == 1
        assert summary['verified'] == 1      # nur symbolic zählt

    def test_status_summary_strength_symbolic(self):
        """Wenn symbolische Verifikation vorhanden: strength = 'symbolic'."""
        x = sp.Symbol('x')
        proof = ProofChain("Stärke-Test")
        proof.step("x=x", "Identität", x, x)
        summary = proof.status_summary()
        assert summary['strength'] == 'symbolic'

    def test_status_summary_strength_numeric(self):
        """Wenn nur numerische Verifikation: strength = 'numeric'."""
        proof = ProofChain("Numerisch-Test")
        step = ProofStep("f(1)=1", "Test", ProofStatus.VERIFIED_NUMERIC)
        proof.steps.append(step)
        summary = proof.status_summary()
        assert summary['strength'] == 'numeric'

    def test_status_summary_strength_informal(self):
        """Wenn nur ASSUMED/UNVERIFIED: strength = 'informal'."""
        proof = ProofChain("Informal")
        proof.assume("Annahme")
        summary = proof.status_summary()
        assert summary['strength'] == 'informal'

    def test_status_summary_strength_failed(self):
        """Wenn FAILED vorhanden: strength = 'failed'."""
        proof = ProofChain("Fail-Test")
        step = ProofStep("Falsch", "Fehler", ProofStatus.FAILED)
        proof.steps.append(step)
        summary = proof.status_summary()
        assert summary['strength'] == 'failed'

    def test_to_text_returns_nonempty_string(self):
        """to_text() muss einen nicht-leeren String zurückgeben."""
        proof = ProofChain("Test-Theorem", "Eine Beschreibung")
        proof.assume("Annahme A")
        proof.conclude("Schluss B")
        text = proof.to_text()
        assert isinstance(text, str)
        assert len(text) > 0

    def test_to_text_contains_theorem(self):
        """to_text() muss das Theorem enthalten."""
        proof = ProofChain("Mein besonderes Theorem")
        text = proof.to_text()
        assert "Mein besonderes Theorem" in text

    def test_to_html_returns_div(self):
        """to_html() muss einen String mit '<div' enthalten."""
        proof = ProofChain("HTML-Test")
        proof.assume("Voraussetzung")
        html = proof.to_html()
        assert "<div" in html

    def test_to_html_contains_proof_container(self):
        """to_html() muss 'proof-container' CSS-Klasse enthalten."""
        proof = ProofChain("Container-Test")
        html = proof.to_html()
        assert "proof-container" in html

    def test_to_html_color_coding(self):
        """to_html() muss verschiedene CSS-Klassen für verschiedene Status verwenden."""
        x = sp.Symbol('x')
        proof = ProofChain("Farbkodierung")
        proof.assume("Annahme")                       # proof-assumed
        proof.step("Symbolisch", "ok", x, x)          # proof-symbolic
        html = proof.to_html()
        assert "proof-assumed" in html
        assert "proof-symbolic" in html

    def test_repr_format(self):
        """__repr__ muss Theorem und Verifikationsstand enthalten."""
        proof = ProofChain("Repr-Theorem")
        repr_str = repr(proof)
        assert "Repr-Theorem" in repr_str
        assert "verifiziert" in repr_str

    def test_numeric_step_verifies_range(self):
        """numeric_step() soll einen Schritt mit numerischer Verifikation hinzufügen."""
        proof = ProofChain("Numerische Schritte")
        # P(n): n * 2 == n + n für alle n im Bereich
        proof.numeric_step(
            "2n = n + n für alle n ∈ [0, 10]",
            "Distributivgesetz",
            func=lambda n: n * 2 == n + n,
            test_range=range(0, 11)
        )
        assert len(proof.steps) == 1
        assert proof.steps[0].status == ProofStatus.VERIFIED_NUMERIC


# ===========================================================================
# TESTS: InductionProof
# ===========================================================================

class TestInductionProof:
    """
    @brief Tests für die InductionProof-Klasse.
    @description Testet Basisfall-Verifikation und Induktionsschritte.
    @date 2026-03-10
    """

    def test_base_case_gauss_sum_n1(self):
        """
        Basisfall der Gaußschen Summenformel für n=1:
        Σ(k=1..1) k = 1 = 1*(1+1)/2 = 1 ✓
        """
        proof = InductionProof(
            "Σk = n(n+1)/2",
            predicate=lambda n: sum(range(1, n + 1)) == n * (n + 1) // 2,
            n0=1
        )
        result = proof.verify_base_case()
        assert result is True
        assert proof.base_case_verified is True

    def test_base_case_adds_step(self):
        """verify_base_case() muss einen Schritt zur Kette hinzufügen."""
        proof = InductionProof(
            "Test-Induktion",
            predicate=lambda n: n >= 0,
            n0=0
        )
        initial_steps = len(proof.steps)
        proof.verify_base_case()
        assert len(proof.steps) == initial_steps + 1

    def test_base_case_fails_for_wrong_predicate(self):
        """Basisfall für ein falsches Prädikat muss False zurückgeben."""
        proof = InductionProof(
            "Immer falsch",
            predicate=lambda n: False,  # Immer falsch
            n0=0
        )
        result = proof.verify_base_case()
        assert result is False
        assert proof.base_case_verified is False

    def test_empirical_verification_gauss_sum_50(self):
        """
        Empirische Verifikation der Gaußschen Summe für n=1..50 muss True sein.
        """
        proof = InductionProof(
            "Gaußsche Summe",
            predicate=lambda n: sum(range(1, n + 1)) == n * (n + 1) // 2,
            n0=1
        )
        result = proof.verify_induction_empirically(50)
        assert result is True
        assert proof.induction_step_verified is True

    def test_empirical_verification_fails_for_wrong_predicate(self):
        """Empirische Verifikation für ein falsches Prädikat muss False sein."""
        # P(n): n² == n+1 ist fast immer falsch
        proof = InductionProof(
            "Falsche Aussage",
            predicate=lambda n: n * n == n + 1,
            n0=2
        )
        result = proof.verify_induction_empirically(10)
        assert result is False

    def test_symbolic_induction_gauss(self):
        """
        Symbolischer Beweis der Gaußschen Summe:
        Basisfall: 1 == 1*(1+1)/2
        Schritt: n(n+1)/2 + (n+1) == (n+1)(n+2)/2
        """
        n = sp.Symbol('n', positive=True, integer=True)
        proof = InductionProof(
            "Gauß symbolisch",
            predicate=lambda m: sum(range(1, m+1)) == m*(m+1)//2,
            n0=1
        )
        # Basisfall: 1 = 1·2/2
        base_lhs = sp.Integer(1)
        base_rhs = sp.Rational(1 * 2, 2)
        # Induktionsschritt
        step_lhs = sp.Rational(1, 2) * n * (n + 1) + (n + 1)
        step_rhs = sp.Rational(1, 2) * (n + 1) * (n + 2)

        result = proof.prove_by_symbolic_induction(
            base_lhs, base_rhs, step_lhs, step_rhs, n
        )
        assert result is True
        assert proof.base_case_verified is True
        assert proof.induction_step_verified is True

    def test_symbolic_induction_fails_for_wrong_step(self):
        """
        Symbolischer Beweis mit falscher rechter Seite muss False zurückgeben.
        """
        n = sp.Symbol('n')
        proof = InductionProof("Falsch", predicate=lambda m: True, n0=0)
        result = proof.prove_by_symbolic_induction(
            sp.Integer(0), sp.Integer(0),  # Basisfall stimmt
            n + 1, n + 2,                  # Schritt stimmt NICHT
            n
        )
        assert result is False


# ===========================================================================
# TESTS: ContradictionProof
# ===========================================================================

class TestContradictionProof:
    """
    @brief Tests für ContradictionProof (Widerspruchsbeweise).
    @description Testet Negation, Widerspruchsableitung und Abschluss.
    @date 2026-03-10
    """

    def test_negate_adds_step(self):
        """negate() muss einen ASSUMED-Schritt mit Negation hinzufügen."""
        proof = ContradictionProof("A ist wahr")
        proof.negate("Angenommen: A ist falsch")
        assert len(proof.steps) == 1
        assert proof.steps[0].status == ProofStatus.ASSUMED
        assert proof.negation == "Angenommen: A ist falsch"

    def test_derive_contradiction_sets_flag(self):
        """derive_contradiction() muss contradiction_found auf True setzen."""
        proof = ContradictionProof("Test-Theorem")
        proof.derive_contradiction("Widerspruch mit Axiom X", "Axiom X gilt immer")
        assert proof.contradiction_found is True

    def test_qed_completes_proof(self):
        """qed() muss is_complete auf True setzen."""
        proof = ContradictionProof("Test")
        proof.negate("Negation")
        proof.derive_contradiction("Widerspruch", "Begründung")
        proof.qed()
        assert proof.is_complete is True

    def test_qed_adds_theorem_to_conclusion(self):
        """qed() muss das Theorem in der Schlussfolgerung erwähnen."""
        proof = ContradictionProof("√2 ist irrational")
        proof.negate("√2 ist rational")
        proof.derive_contradiction("ggT > 1", "Beide Seiten gerade")
        proof.qed()
        # Die Schlussfolgerung muss das Theorem enthalten
        last_step = proof.steps[-1]
        assert "√2 ist irrational" in last_step.statement

    def test_full_contradiction_proof_structure(self):
        """
        Ein vollständiger Widerspruchsbeweis muss korrekte Struktur haben:
        negate → steps → derive_contradiction → qed
        """
        proof = ContradictionProof("∞ viele Primzahlen")
        proof.negate("Endlich viele Primzahlen")
        proof.step("Produkt + 1 konstruieren", "Euklid")
        proof.derive_contradiction("Neue Primzahl gefunden", "Kein Teiler in Liste")
        proof.qed()

        assert proof.is_complete is True
        assert proof.contradiction_found is True
        # Mindestens: negate, step, widerspruch, qed
        assert len(proof.steps) >= 4


# ===========================================================================
# TESTS: TheoremRegistry
# ===========================================================================

class TestTheoremRegistry:
    """
    @brief Tests für die TheoremRegistry-Klasse.
    @description Testet Registrierung, Abfrage und Laden klassischer Sätze.
    @date 2026-03-10
    """

    def setup_method(self):
        """Setzt die Registry vor jedem Test zurück."""
        # Registry leeren für isolierte Tests
        TheoremRegistry._registry.clear()

    def test_register_and_get(self):
        """Registrierter Beweis muss über get() abrufbar sein."""
        proof = ProofChain("Test-Satz")
        TheoremRegistry.register("test_theorem", proof)
        retrieved = TheoremRegistry.get("test_theorem")
        assert retrieved is proof

    def test_get_nonexistent_returns_none(self):
        """get() für unbekannten Namen muss None zurückgeben."""
        result = TheoremRegistry.get("nicht_vorhanden")
        assert result is None

    def test_list_theorems_after_register(self):
        """list_theorems() muss registrierte Namen enthalten."""
        TheoremRegistry.register("satz_a", ProofChain("A"))
        TheoremRegistry.register("satz_b", ProofChain("B"))
        theorems = TheoremRegistry.list_theorems()
        assert "satz_a" in theorems
        assert "satz_b" in theorems

    def test_list_theorems_is_sorted(self):
        """list_theorems() muss alphabetisch sortiert sein."""
        TheoremRegistry.register("zebra", ProofChain("Z"))
        TheoremRegistry.register("apfel", ProofChain("A"))
        theorems = TheoremRegistry.list_theorems()
        assert theorems == sorted(theorems)

    def test_load_classics_registers_at_least_four(self):
        """load_classics() muss mindestens 4 Sätze registrieren."""
        TheoremRegistry.load_classics()
        theorems = TheoremRegistry.list_theorems()
        assert len(theorems) >= 4

    def test_load_classics_contains_gauss_sum(self):
        """load_classics() muss 'gauss_sum' enthalten."""
        TheoremRegistry.load_classics()
        assert TheoremRegistry.get("gauss_sum") is not None

    def test_load_classics_contains_sqrt2(self):
        """load_classics() muss 'sqrt2_irrational' enthalten."""
        TheoremRegistry.load_classics()
        assert TheoremRegistry.get("sqrt2_irrational") is not None

    def test_load_classics_contains_primes(self):
        """load_classics() muss 'infinitely_many_primes' enthalten."""
        TheoremRegistry.load_classics()
        assert TheoremRegistry.get("infinitely_many_primes") is not None

    def test_load_classics_contains_fermat(self):
        """load_classics() muss 'fermat_little' enthalten."""
        TheoremRegistry.load_classics()
        assert TheoremRegistry.get("fermat_little") is not None


# ===========================================================================
# TESTS: prove_gauss_sum()
# ===========================================================================

class TestProveGaussSum:
    """
    @brief Tests für den vorgefertigten Gauß-Summen-Beweis.
    @description Verifiziert Typ, Vollständigkeit und Korrektheit.
    @date 2026-03-10
    """

    def test_returns_induction_proof(self):
        """prove_gauss_sum() muss eine InductionProof-Instanz zurückgeben."""
        proof = prove_gauss_sum()
        assert isinstance(proof, InductionProof)

    def test_is_complete(self):
        """Der Beweis muss vollständig sein (is_complete = True)."""
        proof = prove_gauss_sum()
        assert proof.is_complete is True

    def test_has_steps(self):
        """Der Beweis muss mindestens 4 Schritte enthalten."""
        proof = prove_gauss_sum()
        assert len(proof.steps) >= 4

    def test_theorem_in_title(self):
        """Das Theorem muss die Gaußsche Formel im Namen enthalten."""
        proof = prove_gauss_sum()
        assert "n(n+1)/2" in proof.theorem or "Σ" in proof.theorem

    def test_has_verified_steps(self):
        """Der Beweis muss mindestens einen verifizierten Schritt haben."""
        proof = prove_gauss_sum()
        summary = proof.status_summary()
        assert summary['verified'] > 0

    def test_base_case_verified(self):
        """Der Basisfall muss verifiziert sein."""
        proof = prove_gauss_sum()
        assert proof.base_case_verified is True

    def test_strength_is_symbolic_or_numeric(self):
        """Die Beweisstärke muss symbolic oder numeric sein."""
        proof = prove_gauss_sum()
        summary = proof.status_summary()
        assert summary['strength'] in ('symbolic', 'numeric')

    def test_to_text_nonempty(self):
        """to_text() muss einen nicht-leeren String zurückgeben."""
        proof = prove_gauss_sum()
        text = proof.to_text()
        assert len(text) > 0

    def test_to_html_contains_div(self):
        """to_html() muss '<div' enthalten."""
        proof = prove_gauss_sum()
        html = proof.to_html()
        assert "<div" in html


# ===========================================================================
# TESTS: prove_sqrt2_irrational()
# ===========================================================================

class TestProveSqrt2Irrational:
    """
    @brief Tests für den Irrationalitätsbeweis von √2.
    @date 2026-03-10
    """

    def test_returns_contradiction_proof(self):
        """prove_sqrt2_irrational() muss eine ContradictionProof-Instanz liefern."""
        proof = prove_sqrt2_irrational()
        assert isinstance(proof, ContradictionProof)

    def test_is_complete(self):
        """Der Beweis muss vollständig sein."""
        proof = prove_sqrt2_irrational()
        assert proof.is_complete is True

    def test_contradiction_found(self):
        """Ein Widerspruch muss explizit gefunden worden sein."""
        proof = prove_sqrt2_irrational()
        assert proof.contradiction_found is True

    def test_theorem_mentions_sqrt2(self):
        """Das Theorem muss √2 erwähnen."""
        proof = prove_sqrt2_irrational()
        assert "√2" in proof.theorem or "sqrt" in proof.theorem.lower()

    def test_has_multiple_steps(self):
        """Der Beweis muss mehrere Schritte enthalten."""
        proof = prove_sqrt2_irrational()
        assert len(proof.steps) >= 4

    def test_negation_is_set(self):
        """Die Negation muss gesetzt sein."""
        proof = prove_sqrt2_irrational()
        assert len(proof.negation) > 0

    def test_to_html_contains_div(self):
        """to_html() muss '<div' enthalten."""
        proof = prove_sqrt2_irrational()
        assert "<div" in proof.to_html()


# ===========================================================================
# TESTS: prove_infinitely_many_primes()
# ===========================================================================

class TestProveInfinitelyManyPrimes:
    """
    @brief Tests für Euklids Beweis unendlich vieler Primzahlen.
    @date 2026-03-10
    """

    def test_returns_contradiction_proof(self):
        """prove_infinitely_many_primes() muss ContradictionProof zurückgeben."""
        proof = prove_infinitely_many_primes()
        assert isinstance(proof, ContradictionProof)

    def test_is_complete(self):
        """Der Beweis muss vollständig sein."""
        proof = prove_infinitely_many_primes()
        assert proof.is_complete is True

    def test_contradiction_found(self):
        """Ein Widerspruch muss gefunden worden sein."""
        proof = prove_infinitely_many_primes()
        assert proof.contradiction_found is True

    def test_theorem_mentions_primes(self):
        """Das Theorem muss Primzahlen erwähnen."""
        proof = prove_infinitely_many_primes()
        assert "Prim" in proof.theorem or "prim" in proof.theorem.lower()

    def test_has_euclid_steps(self):
        """Der Beweis muss mindestens 4 Schritte haben (Euklids Schema)."""
        proof = prove_infinitely_many_primes()
        assert len(proof.steps) >= 4


# ===========================================================================
# TESTS: prove_fermat_little_theorem()
# ===========================================================================

class TestProveFermatLittleTheorem:
    """
    @brief Tests für Fermats kleinen Satz.
    @date 2026-03-10
    """

    def test_returns_proof_chain(self):
        """prove_fermat_little_theorem() muss eine ProofChain zurückgeben."""
        proof = prove_fermat_little_theorem(p=5)
        assert isinstance(proof, ProofChain)

    def test_is_complete(self):
        """Der Beweis muss vollständig sein."""
        proof = prove_fermat_little_theorem(p=5)
        assert proof.is_complete is True

    def test_theorem_contains_p(self):
        """Das Theorem muss die verwendete Primzahl p enthalten."""
        proof = prove_fermat_little_theorem(p=7)
        assert "7" in proof.theorem

    def test_has_numeric_verified_steps(self):
        """Der Beweis muss numerisch verifizierte Schritte haben."""
        proof = prove_fermat_little_theorem(p=5)
        summary = proof.status_summary()
        assert summary['numeric'] > 0

    def test_no_failed_steps_for_prime(self):
        """Für eine echte Primzahl dürfen keine FAILED-Schritte auftreten."""
        proof = prove_fermat_little_theorem(p=11)
        summary = proof.status_summary()
        assert summary['failed'] == 0

    def test_all_a_values_verified(self):
        """Für p=5 müssen alle 5 Werte a=0..4 verifiziert sein."""
        proof = prove_fermat_little_theorem(p=5)
        # Mindestens p numerisch verifizierte Schritte (plus Annahmen)
        numeric_steps = sum(
            1 for s in proof.steps
            if s.status == ProofStatus.VERIFIED_NUMERIC
        )
        assert numeric_steps >= 5


# ===========================================================================
# TESTS: riemann_hypothesis_evidence()
# ===========================================================================

class TestRiemannHypothesisEvidence:
    """
    @brief Tests für die Riemann-Hypothesen-Evidenz.
    @date 2026-03-10
    """

    def test_returns_proof_chain(self):
        """riemann_hypothesis_evidence() muss eine ProofChain zurückgeben."""
        proof = riemann_hypothesis_evidence()
        assert isinstance(proof, ProofChain)

    def test_is_not_complete(self):
        """Die Riemann-Hypothese ist nicht vollständig bewiesen."""
        proof = riemann_hypothesis_evidence()
        assert proof.is_complete is False

    def test_has_pending_step(self):
        """Muss mindestens einen PENDING-Schritt enthalten."""
        proof = riemann_hypothesis_evidence()
        pending_steps = [s for s in proof.steps if s.status == ProofStatus.PENDING]
        assert len(pending_steps) >= 1

    def test_has_verified_numeric_steps(self):
        """Muss numerisch verifizierte Schritte enthalten."""
        proof = riemann_hypothesis_evidence()
        numeric_steps = [s for s in proof.steps if s.status == ProofStatus.VERIFIED_NUMERIC]
        assert len(numeric_steps) >= 1

    def test_theorem_mentions_riemann_or_nullstellen(self):
        """Das Theorem muss Riemann oder Nullstellen erwähnen."""
        proof = riemann_hypothesis_evidence()
        assert any(
            word in proof.theorem
            for word in ["Riemann", "Nullstell", "ζ", "Re(s)"]
        )


# ===========================================================================
# TESTS: goldbach_evidence()
# ===========================================================================

class TestGoldbachEvidence:
    """
    @brief Tests für die Goldbach-Evidenz.
    @date 2026-03-10
    """

    def test_returns_proof_chain_100(self):
        """goldbach_evidence(100) muss eine ProofChain zurückgeben."""
        proof = goldbach_evidence(n_max=100)
        assert isinstance(proof, ProofChain)

    def test_is_complete_100(self):
        """goldbach_evidence(100) muss vollständig sein."""
        proof = goldbach_evidence(n_max=100)
        assert proof.is_complete is True

    def test_has_verified_numeric_status_100(self):
        """
        goldbach_evidence(100) muss mindestens einen VERIFIED_NUMERIC-Schritt haben.
        """
        proof = goldbach_evidence(n_max=100)
        numeric_steps = [s for s in proof.steps if s.status == ProofStatus.VERIFIED_NUMERIC]
        assert len(numeric_steps) >= 1

    def test_strength_is_numeric_100(self):
        """Die Beweisstärke muss 'numeric' sein (kein symbolischer Beweis)."""
        proof = goldbach_evidence(n_max=100)
        summary = proof.status_summary()
        assert summary['strength'] == 'numeric'

    def test_no_failed_steps_100(self):
        """Für n_max=100 dürfen keine FAILED-Schritte auftreten."""
        proof = goldbach_evidence(n_max=100)
        summary = proof.status_summary()
        assert summary['failed'] == 0

    def test_goldbach_evidence_10000(self):
        """
        goldbach_evidence(10000): Vollständige Verifikation bis 10000.
        Größerer Test — prüft Effizienz und Korrektheit.
        """
        proof = goldbach_evidence(n_max=10000)
        assert isinstance(proof, ProofChain)
        assert proof.is_complete is True
        summary = proof.status_summary()
        assert summary['strength'] == 'numeric'
        assert summary['failed'] == 0

    def test_theorem_mentions_goldbach(self):
        """Das Theorem muss Goldbach oder Primzahl erwähnen."""
        proof = goldbach_evidence(n_max=50)
        assert any(
            word in proof.theorem
            for word in ["Goldbach", "Prim", "gerade"]
        )


# ===========================================================================
# EDGE-CASE TESTS
# ===========================================================================

class TestEdgeCases:
    """
    @brief Tests für Randfälle und Sondersituationen.
    @description Prüft robustes Verhalten bei ungewöhnlichen Eingaben.
    @date 2026-03-10
    """

    def test_empty_proof_chain_summary(self):
        """Eine leere ProofChain hat sinnvolle Standardwerte."""
        proof = ProofChain("Leer")
        summary = proof.status_summary()
        assert summary['total_steps'] == 0
        assert summary['verified'] == 0
        assert summary['is_complete'] is False
        assert summary['strength'] == 'informal'

    def test_proof_step_with_exception_in_func(self):
        """verify_numeric muss mit fehlerhafter Funktion umgehen können."""
        step = ProofStep("Fehlertest", "Test")
        def bad_func(x):
            raise ValueError("Fehler!")
        result = step.verify_numeric(bad_func, [1], [1])
        assert result is False
        assert step.status == ProofStatus.FAILED

    def test_contradiction_proof_without_derive(self):
        """ContradictionProof ohne derive_contradiction hat contradiction_found=False."""
        proof = ContradictionProof("Theorem ohne Widerspruch")
        proof.negate("Negation")
        assert proof.contradiction_found is False

    def test_induction_proof_n0_zero(self):
        """InductionProof mit n0=0 muss Basisfall P(0) prüfen."""
        # Summe aller Zahlen 1..0 = 0 = 0*(0+1)/2 = 0
        proof = InductionProof(
            "Summe ist 0 für n=0",
            predicate=lambda n: sum(range(1, n + 1)) == n * (n + 1) // 2,
            n0=0
        )
        result = proof.verify_base_case()
        assert result is True

    def test_proof_chain_description_in_html(self):
        """to_html() muss die Beschreibung enthalten, wenn angegeben."""
        proof = ProofChain("Theorem", description="Meine spezielle Beschreibung")
        html = proof.to_html()
        assert "Meine spezielle Beschreibung" in html

    def test_proof_chain_description_in_text(self):
        """to_text() muss die Beschreibung enthalten."""
        proof = ProofChain("Theorem", description="Text-Beschreibung XYZ")
        text = proof.to_text()
        assert "Text-Beschreibung XYZ" in text

    def test_fermat_theorem_for_p2(self):
        """Fermats kleiner Satz für p=2 (kleinste Primzahl)."""
        proof = prove_fermat_little_theorem(p=2)
        assert proof.is_complete is True
        summary = proof.status_summary()
        assert summary['failed'] == 0

    def test_fermat_theorem_for_p13(self):
        """Fermats kleiner Satz für p=13."""
        proof = prove_fermat_little_theorem(p=13)
        assert proof.is_complete is True
        summary = proof.status_summary()
        assert summary['failed'] == 0

    def test_induction_proof_inherits_from_proof_chain(self):
        """InductionProof muss von ProofChain erben."""
        proof = InductionProof("Test", predicate=lambda n: True, n0=0)
        assert isinstance(proof, ProofChain)

    def test_contradiction_proof_inherits_from_proof_chain(self):
        """ContradictionProof muss von ProofChain erben."""
        proof = ContradictionProof("Test")
        assert isinstance(proof, ProofChain)

    def test_theorem_registry_overwrite(self):
        """Registrierung unter gleichem Namen überschreibt alten Eintrag."""
        TheoremRegistry._registry.clear()
        proof1 = ProofChain("Version 1")
        proof2 = ProofChain("Version 2")
        TheoremRegistry.register("test", proof1)
        TheoremRegistry.register("test", proof2)
        assert TheoremRegistry.get("test") is proof2

    def test_verify_symbolic_with_complex_expression(self):
        """
        Symbolische Verifikation einer komplexeren Identität:
        sin²(x) + cos²(x) = 1
        """
        x = sp.Symbol('x')
        step = ProofStep("Pythagoreische Identität", "Trigonometrie")
        lhs = sp.sin(x)**2 + sp.cos(x)**2
        rhs = sp.Integer(1)
        result = step.verify_symbolic(lhs, rhs)
        assert result is True
        assert step.status == ProofStatus.VERIFIED_SYMBOLIC

    def test_goldbach_evidence_small_n(self):
        """goldbach_evidence(10) muss für sehr kleine n funktionieren."""
        proof = goldbach_evidence(n_max=10)
        assert isinstance(proof, ProofChain)
        assert proof.is_complete is True
