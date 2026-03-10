"""
Tests für das Modul mathematical_logic.py
==========================================
Testet Aussagenlogik, Resolution, SAT, Prädikatenlogik,
Modallogik und Gödel-Nummerierung.

Autor: Kurt Ingwer
Letzte Änderung: 2026-03-10
"""

import pytest
from src.mathematical_logic import (
    # Aussagenlogik
    Proposition, LogicFormula,
    truth_table, is_tautology, is_satisfiable, is_contradiction,
    logical_consequence, logical_equivalence,
    # Resolution / SAT
    resolution_step, resolution_refutation, dpll, tseitin_transform,
    # Prädikatenlogik
    Term, Predicate, FOLFormula,
    prenex_normal_form, skolemize, herbrand_universe,
    # Beweiskalküle
    NaturalDeductionProof, HilbertSystem, modus_ponens, modus_tollens,
    # Modallogik
    KripkeFrame, ModalFormula, modal_system,
    # Gödel
    godel_numbering, godel_decode, incompleteness_demonstration,
)


# =============================================================================
# Hilfsfunktionen zum Erstellen von Formeln
# =============================================================================

def p():
    """Atomare Variable p."""
    return LogicFormula.atom("p")

def q():
    """Atomare Variable q."""
    return LogicFormula.atom("q")

def r():
    """Atomare Variable r."""
    return LogicFormula.atom("r")


# =============================================================================
# 1. PROPOSITION TESTS
# =============================================================================

class TestProposition:
    """Tests für die Proposition-Klasse."""

    def test_proposition_creation(self):
        """Atomare Aussage wird korrekt erstellt."""
        prop = Proposition("p", True)
        assert prop.name == "p"
        assert prop.value is True

    def test_proposition_evaluate_from_assignment(self):
        """Aussage wird korrekt aus Belegung ausgewertet."""
        prop = Proposition("p")
        assert prop.evaluate({"p": True}) is True
        assert prop.evaluate({"p": False}) is False

    def test_proposition_evaluate_from_value(self):
        """Aussage mit vorgegebenem Wert wird korrekt ausgewertet."""
        prop = Proposition("p", True)
        assert prop.evaluate({}) is True

    def test_proposition_str(self):
        """String-Darstellung ist der Name."""
        prop = Proposition("Regen")
        assert str(prop) == "Regen"

    def test_proposition_missing_variable_raises(self):
        """Fehlende Variable in Belegung wirft KeyError."""
        prop = Proposition("x")
        with pytest.raises(KeyError):
            prop.evaluate({"y": True})


# =============================================================================
# 2. LOGICFORMULA – GRUNDOPERATOREN
# =============================================================================

class TestLogicFormulaBasic:
    """Tests für grundlegende Operatoren der LogicFormula."""

    def test_atom_true(self):
        """Atom mit True-Belegung ist wahr."""
        f = p()
        assert f.evaluate({"p": True}) is True

    def test_atom_false(self):
        """Atom mit False-Belegung ist falsch."""
        f = p()
        assert f.evaluate({"p": False}) is False

    def test_not_true(self):
        """Negation von True ist False."""
        f = LogicFormula("NOT", p())
        assert f.evaluate({"p": True}) is False

    def test_not_false(self):
        """Negation von False ist True."""
        f = LogicFormula("NOT", p())
        assert f.evaluate({"p": False}) is True

    def test_and_tt(self):
        """True AND True ist True."""
        f = LogicFormula("AND", p(), q())
        assert f.evaluate({"p": True, "q": True}) is True

    def test_and_tf(self):
        """True AND False ist False."""
        f = LogicFormula("AND", p(), q())
        assert f.evaluate({"p": True, "q": False}) is False

    def test_or_ff(self):
        """False OR False ist False."""
        f = LogicFormula("OR", p(), q())
        assert f.evaluate({"p": False, "q": False}) is False

    def test_or_tf(self):
        """True OR False ist True."""
        f = LogicFormula("OR", p(), q())
        assert f.evaluate({"p": True, "q": False}) is True

    def test_implies_ff(self):
        """False → False ist True (Vakuumwahrheit)."""
        f = LogicFormula("IMPLIES", p(), q())
        assert f.evaluate({"p": False, "q": False}) is True

    def test_implies_tf(self):
        """True → False ist False."""
        f = LogicFormula("IMPLIES", p(), q())
        assert f.evaluate({"p": True, "q": False}) is False

    def test_implies_ft(self):
        """False → True ist True."""
        f = LogicFormula("IMPLIES", p(), q())
        assert f.evaluate({"p": False, "q": True}) is True

    def test_iff_tt(self):
        """True ↔ True ist True."""
        f = LogicFormula("IFF", p(), q())
        assert f.evaluate({"p": True, "q": True}) is True

    def test_iff_tf(self):
        """True ↔ False ist False."""
        f = LogicFormula("IFF", p(), q())
        assert f.evaluate({"p": True, "q": False}) is False

    def test_xor_tt(self):
        """True XOR True ist False."""
        f = LogicFormula("XOR", p(), q())
        assert f.evaluate({"p": True, "q": True}) is False

    def test_xor_tf(self):
        """True XOR False ist True."""
        f = LogicFormula("XOR", p(), q())
        assert f.evaluate({"p": True, "q": False}) is True

    def test_nand_tt(self):
        """True NAND True ist False."""
        f = LogicFormula("NAND", p(), q())
        assert f.evaluate({"p": True, "q": True}) is False

    def test_nor_ff(self):
        """False NOR False ist True."""
        f = LogicFormula("NOR", p(), q())
        assert f.evaluate({"p": False, "q": False}) is True

    def test_invalid_operator_raises(self):
        """Ungültiger Operator wirft ValueError."""
        with pytest.raises(ValueError):
            LogicFormula("INVALID", p())

    def test_str_atom(self):
        """Atom-String-Darstellung."""
        assert str(p()) == "p"

    def test_str_not(self):
        """NOT-String-Darstellung."""
        assert str(LogicFormula("NOT", p())) == "¬p"

    def test_str_and(self):
        """AND-String-Darstellung."""
        f = LogicFormula("AND", p(), q())
        assert "p" in str(f) and "q" in str(f) and "∧" in str(f)

    def test_str_or(self):
        """OR-String-Darstellung."""
        f = LogicFormula("OR", p(), q())
        assert "∨" in str(f)

    def test_str_implies(self):
        """IMPLIES-String-Darstellung."""
        f = LogicFormula("IMPLIES", p(), q())
        assert "→" in str(f)

    def test_get_variables(self):
        """Variablen werden korrekt extrahiert."""
        f = LogicFormula("AND", p(), LogicFormula("OR", q(), r()))
        variables = f.get_variables()
        assert variables == {"p", "q", "r"}


# =============================================================================
# 3. WAHRHEITSTABELLEN UND PRÜFFUNKTIONEN
# =============================================================================

class TestTruthTableAndChecks:
    """Tests für Wahrheitstabellen und logische Prüfungen."""

    def test_truth_table_size(self):
        """Wahrheitstabelle hat 2^n Zeilen."""
        f = LogicFormula("AND", p(), q())
        table = truth_table(f, ["p", "q"])
        assert len(table) == 4

    def test_truth_table_and(self):
        """AND-Wahrheitstabelle ist korrekt."""
        f = LogicFormula("AND", p(), q())
        table = truth_table(f, ["p", "q"])
        # Nur TT → True
        true_rows = [row for row in table if row["result"]]
        assert len(true_rows) == 1
        assert true_rows[0]["p"] is True and true_rows[0]["q"] is True

    def test_tautology_excluded_middle(self):
        """A ∨ ¬A ist Tautologie (ausgeschlossenes Drittes)."""
        f = LogicFormula("OR", p(), LogicFormula("NOT", p()))
        assert is_tautology(f, ["p"]) is True

    def test_not_tautology(self):
        """p ∧ q ist keine Tautologie."""
        f = LogicFormula("AND", p(), q())
        assert is_tautology(f, ["p", "q"]) is False

    def test_contradiction_basic(self):
        """A ∧ ¬A ist Kontradiktion."""
        f = LogicFormula("AND", p(), LogicFormula("NOT", p()))
        assert is_contradiction(f, ["p"]) is True

    def test_not_contradiction(self):
        """A ∨ B ist keine Kontradiktion."""
        f = LogicFormula("OR", p(), q())
        assert is_contradiction(f, ["p", "q"]) is False

    def test_satisfiable_true(self):
        """p ist erfüllbar."""
        f = p()
        assert is_satisfiable(f, ["p"]) is True

    def test_satisfiable_false(self):
        """A ∧ ¬A ist nicht erfüllbar."""
        f = LogicFormula("AND", p(), LogicFormula("NOT", p()))
        assert is_satisfiable(f, ["p"]) is False

    def test_logical_consequence_mp(self):
        """Modus Ponens als logische Konsequenz."""
        # p, p→q ⊨ q
        prem1 = p()
        prem2 = LogicFormula("IMPLIES", p(), q())
        conc = q()
        assert logical_consequence([prem1, prem2], conc, ["p", "q"]) is True

    def test_logical_consequence_fails(self):
        """Ungültige Schlussfolgerung wird erkannt."""
        # p ⊭ q
        assert logical_consequence([p()], q(), ["p", "q"]) is False

    def test_logical_equivalence_de_morgan_and(self):
        """De Morgan: ¬(p ∧ q) ≡ ¬p ∨ ¬q."""
        f1 = LogicFormula("NOT", LogicFormula("AND", p(), q()))
        f2 = LogicFormula("OR", LogicFormula("NOT", p()), LogicFormula("NOT", q()))
        assert logical_equivalence(f1, f2, ["p", "q"]) is True

    def test_logical_equivalence_de_morgan_or(self):
        """De Morgan: ¬(p ∨ q) ≡ ¬p ∧ ¬q."""
        f1 = LogicFormula("NOT", LogicFormula("OR", p(), q()))
        f2 = LogicFormula("AND", LogicFormula("NOT", p()), LogicFormula("NOT", q()))
        assert logical_equivalence(f1, f2, ["p", "q"]) is True

    def test_logical_equivalence_implication(self):
        """p → q ≡ ¬p ∨ q."""
        f1 = LogicFormula("IMPLIES", p(), q())
        f2 = LogicFormula("OR", LogicFormula("NOT", p()), q())
        assert logical_equivalence(f1, f2, ["p", "q"]) is True

    def test_double_negation_equivalence(self):
        """¬¬p ≡ p."""
        f1 = LogicFormula("NOT", LogicFormula("NOT", p()))
        f2 = p()
        assert logical_equivalence(f1, f2, ["p"]) is True


# =============================================================================
# 4. CNF / DNF KONVERSION
# =============================================================================

class TestNormalForms:
    """Tests für CNF- und DNF-Konversion."""

    def test_cnf_atom(self):
        """Atom in CNF ist Atom selbst."""
        cnf = p().to_cnf()
        assert cnf.evaluate({"p": True}) is True

    def test_cnf_preserves_semantics_simple(self):
        """CNF von p∧q hat gleiche Wahrheitswerte."""
        f = LogicFormula("AND", p(), q())
        cnf = f.to_cnf()
        for pv in [True, False]:
            for qv in [True, False]:
                assert cnf.evaluate({"p": pv, "q": qv}) == f.evaluate({"p": pv, "q": qv})

    def test_cnf_preserves_semantics_implies(self):
        """CNF von p→q hat gleiche Wahrheitswerte."""
        f = LogicFormula("IMPLIES", p(), q())
        cnf = f.to_cnf()
        for pv in [True, False]:
            for qv in [True, False]:
                assignment = {"p": pv, "q": qv}
                assert cnf.evaluate(assignment) == f.evaluate(assignment)

    def test_dnf_atom(self):
        """Atom in DNF ist Atom selbst."""
        dnf = p().to_dnf()
        assert dnf.evaluate({"p": False}) is False

    def test_dnf_preserves_semantics(self):
        """DNF von p∨(q∧r) hat gleiche Wahrheitswerte."""
        inner = LogicFormula("AND", q(), r())
        f = LogicFormula("OR", p(), inner)
        dnf = f.to_dnf()
        for pv in [True, False]:
            for qv in [True, False]:
                for rv in [True, False]:
                    assignment = {"p": pv, "q": qv, "r": rv}
                    assert dnf.evaluate(assignment) == f.evaluate(assignment)

    def test_simplify_double_negation(self):
        """¬¬p vereinfacht zu p."""
        f = LogicFormula("NOT", LogicFormula("NOT", p()))
        simplified = f.simplify()
        assert simplified.op == "ATOM"
        assert simplified.args[0] == "p"


# =============================================================================
# 5. RESOLUTION UND DPLL
# =============================================================================

class TestResolutionAndDPLL:
    """Tests für Resolution und DPLL SAT-Solver."""

    def test_resolution_step_basic(self):
        """Resolutionsschritt mit komplementären Literalen."""
        # {p, q} und {-p, r} → {q, r}
        c1 = frozenset({"p", "q"})
        c2 = frozenset({"-p", "r"})
        result = resolution_step(c1, c2)
        assert result == frozenset({"q", "r"})

    def test_resolution_step_empty(self):
        """Resolution ergibt leere Klausel bei komplementären Einheitsklauseln."""
        c1 = frozenset({"p"})
        c2 = frozenset({"-p"})
        result = resolution_step(c1, c2)
        assert result == frozenset()

    def test_resolution_step_none(self):
        """Kein Resolutionsschritt wenn keine komplementären Literale."""
        c1 = frozenset({"p", "q"})
        c2 = frozenset({"r", "s"})
        result = resolution_step(c1, c2)
        assert result is None

    def test_resolution_refutation_unsatisfiable(self):
        """Unerfüllbare Klauselmenge wird erkannt."""
        # {p} und {-p}
        clauses = [frozenset({"p"}), frozenset({"-p"})]
        assert resolution_refutation(clauses) is True

    def test_resolution_refutation_satisfiable(self):
        """Erfüllbare Klauselmenge wird korrekt erkannt."""
        # {p, q} und {-p, q} → q ist immer wahr → erfüllbar
        clauses = [frozenset({"p", "q"}), frozenset({"-p", "q"})]
        assert resolution_refutation(clauses) is False

    def test_dpll_satisfiable(self):
        """DPLL findet Belegung für erfüllbare Formel."""
        # {p} ∧ {q} → p=True, q=True
        clauses = [frozenset({"p"}), frozenset({"q"})]
        result = dpll(clauses)
        assert result is not None
        assert result.get("p") is True
        assert result.get("q") is True

    def test_dpll_unsatisfiable(self):
        """DPLL erkennt Unerfüllbarkeit."""
        # {p} ∧ {-p}
        clauses = [frozenset({"p"}), frozenset({"-p"})]
        assert dpll(clauses) is None

    def test_dpll_complex(self):
        """DPLL löst komplexere Formel (3-SAT Instanz)."""
        # (p ∨ q ∨ r) ∧ (-p ∨ q) ∧ (-q ∨ r) ∧ (-r)
        clauses = [
            frozenset({"p", "q", "r"}),
            frozenset({"-p", "q"}),
            frozenset({"-q", "r"}),
            frozenset({"-r"}),
        ]
        result = dpll(clauses)
        # Prüfe ob Belegung wirklich erfüllt
        if result is not None:
            def eval_clause(clause, assignment):
                for lit in clause:
                    var = lit.lstrip("-")
                    is_neg = lit.startswith("-")
                    val = assignment.get(var, False)
                    if (val and not is_neg) or (not val and is_neg):
                        return True
                return False
            for clause in clauses:
                assert eval_clause(clause, result)

    def test_tseitin_transform_result_is_list(self):
        """Tseitin-Transformation gibt Liste von Klauseln zurück."""
        f = LogicFormula("AND", p(), q())
        clauses = tseitin_transform(f)
        assert isinstance(clauses, list)
        assert len(clauses) > 0
        assert all(isinstance(c, frozenset) for c in clauses)

    def test_tseitin_satisfiable(self):
        """Tseitin-CNF ist erfüllbar wenn Original erfüllbar."""
        f = LogicFormula("OR", p(), q())
        clauses = tseitin_transform(f)
        result = dpll(clauses)
        assert result is not None


# =============================================================================
# 6. PRÄDIKATENLOGIK
# =============================================================================

class TestFirstOrderLogic:
    """Tests für Prädikatenlogik-Klassen."""

    def test_term_variable(self):
        """Variablen-Term wird korrekt erstellt."""
        x = Term.variable("x")
        assert x.is_variable is True
        assert x.name == "x"
        assert x.free_variables() == {"x"}

    def test_term_constant(self):
        """Konstanten-Term wird korrekt erstellt."""
        a = Term.constant("a")
        assert a.is_variable is False
        assert a.free_variables() == set()

    def test_term_function(self):
        """Funktions-Term wird korrekt erstellt."""
        x = Term.variable("x")
        fx = Term.function("f", x)
        assert fx.is_function is True
        assert fx.free_variables() == {"x"}

    def test_term_substitute(self):
        """Variablensubstitution in Termen."""
        x = Term.variable("x")
        a = Term.constant("a")
        result = x.substitute("x", a)
        assert result.name == "a"
        assert result.is_variable is False

    def test_term_str_function(self):
        """String-Darstellung eines Funktionsterms."""
        x = Term.variable("x")
        fx = Term.function("f", x)
        assert str(fx) == "f(x)"

    def test_predicate_creation(self):
        """Prädikat wird korrekt erstellt."""
        P = Predicate("P", 1)
        assert P.name == "P"
        assert P.arity == 1

    def test_predicate_apply(self):
        """Prädikat auf Term anwenden erstellt FOL-Formel."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = P.apply(x)
        assert formula.op == "PRED"

    def test_predicate_wrong_arity(self):
        """Falsche Arität wirft ValueError."""
        P = Predicate("P", 2)
        x = Term.variable("x")
        with pytest.raises(ValueError):
            P.apply(x)  # Nur 1 Argument, braucht 2

    def test_fol_formula_free_variables_forall(self):
        """∀x.P(x) hat keine freien Variablen."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = FOLFormula.forall("x", P.apply(x))
        assert formula.free_variables() == set()

    def test_fol_formula_free_variables_exists(self):
        """∃x.P(x,y) hat y als freie Variable."""
        P = Predicate("P", 2)
        x = Term.variable("x")
        y = Term.variable("y")
        inner = FOLFormula("PRED", P, [x, y])
        formula = FOLFormula.exists("x", inner)
        assert "y" in formula.free_variables()
        assert "x" not in formula.free_variables()

    def test_fol_formula_is_closed(self):
        """∀x.P(x) ist geschlossen."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = FOLFormula.forall("x", P.apply(x))
        assert formula.is_closed() is True

    def test_fol_formula_str_forall(self):
        """String-Darstellung von ∀-Formel."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = FOLFormula.forall("x", P.apply(x))
        s = str(formula)
        assert "∀" in s and "x" in s

    def test_herbrand_universe_no_functions(self):
        """Herbrand-Universum ohne Funktionssymbole sind nur Konstanten."""
        universe = herbrand_universe(["a", "b"], [], depth=2)
        assert "a" in universe
        assert "b" in universe
        assert len(universe) == 2

    def test_herbrand_universe_with_function(self):
        """Herbrand-Universum mit Funktion enthält zusammengesetzte Terme."""
        universe = herbrand_universe(["a"], [("f", 1)], depth=2)
        assert "a" in universe
        assert "f(a)" in universe

    def test_prenex_normal_form_basic(self):
        """PNF einer einfachen Formel wird berechnet."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = FOLFormula.forall("x", P.apply(x))
        pnf = prenex_normal_form(formula)
        assert pnf.op == "FORALL"

    def test_skolemize_exists(self):
        """Existenzquantor wird durch Skolem-Konstante ersetzt."""
        P = Predicate("P", 1)
        x = Term.variable("x")
        formula = FOLFormula.exists("x", P.apply(x))
        skolemized = skolemize(formula)
        # Nach Skolemisierung kein EXISTS mehr
        assert skolemized.op != "EXISTS"

    def test_skolemize_forall_exists(self):
        """∀x.∃y.P(x,y) wird skolemisiert zu ∀x.P(x, sk0(x))."""
        P = Predicate("P", 2)
        x = Term.variable("x")
        y = Term.variable("y")
        inner = FOLFormula("PRED", P, [x, y])
        formula = FOLFormula.forall("x", FOLFormula.exists("y", inner))
        skolemized = skolemize(formula)
        assert skolemized.op == "FORALL"


# =============================================================================
# 7. BEWEISKALKÜLE
# =============================================================================

class TestNaturalDeduction:
    """Tests für Natürliches Schließen."""

    def test_assume_returns_id(self):
        """Annahme gibt eine Zeilen-ID zurück."""
        proof = NaturalDeductionProof()
        line_id = proof.assume(p())
        assert isinstance(line_id, int)
        assert line_id >= 0

    def test_assume_adds_to_open(self):
        """Annahme wird zu offenen Annahmen hinzugefügt."""
        proof = NaturalDeductionProof()
        line_id = proof.assume(p())
        assert line_id in proof.open_assumptions

    def test_apply_and_intro(self):
        """∧I Regel erstellt AND-Formel."""
        proof = NaturalDeductionProof()
        id1 = proof.assume(p())
        id2 = proof.assume(q())
        id3 = proof.apply_rule("∧I", id1, id2)
        formula = proof.lines[id3]["formula"]
        assert isinstance(formula, LogicFormula)
        assert formula.op == "AND"

    def test_apply_mp(self):
        """→E (Modus Ponens) leitet Konsequens ab."""
        proof = NaturalDeductionProof()
        impl = LogicFormula("IMPLIES", p(), q())
        id1 = proof.assume(impl)
        id2 = proof.assume(p())
        id3 = proof.apply_rule("→E", id1, id2)
        formula = proof.lines[id3]["formula"]
        # Konsequens sollte q sein
        assert str(formula) == "q"

    def test_discharge_removes_assumption(self):
        """Entlasten einer Annahme schließt sie."""
        proof = NaturalDeductionProof()
        id1 = proof.assume(p())
        proof.discharge(id1)
        assert id1 not in proof.open_assumptions

    def test_is_valid_with_open_assumptions(self):
        """Beweis mit offenen Annahmen ist ungültig."""
        proof = NaturalDeductionProof()
        proof.assume(p())  # Nicht entlastet
        assert proof.is_valid() is False

    def test_is_valid_no_open_assumptions(self):
        """Beweis ohne offene Annahmen ist gültig."""
        proof = NaturalDeductionProof()
        assert proof.is_valid() is True


class TestModusPonensAndTollens:
    """Tests für MP und MT."""

    def test_modus_ponens_success(self):
        """MP: A→B, A ⊢ B."""
        major = LogicFormula("IMPLIES", p(), q())
        minor = p()
        result = modus_ponens(major, minor)
        assert result is not None
        assert str(result) == "q"

    def test_modus_ponens_wrong_major(self):
        """MP: Kein IMPLIES → None."""
        major = LogicFormula("AND", p(), q())
        result = modus_ponens(major, p())
        assert result is None

    def test_modus_ponens_mismatch(self):
        """MP: Antezedens stimmt nicht → None."""
        major = LogicFormula("IMPLIES", p(), q())
        result = modus_ponens(major, q())
        assert result is None

    def test_modus_tollens_success(self):
        """MT: A→B, ¬B ⊢ ¬A."""
        major = LogicFormula("IMPLIES", p(), q())
        minor = LogicFormula("NOT", q())
        result = modus_tollens(major, minor)
        assert result is not None
        assert result.op == "NOT"
        assert str(result.args[0]) == "p"

    def test_modus_tollens_wrong_minor(self):
        """MT: minor kein NOT → None."""
        major = LogicFormula("IMPLIES", p(), q())
        result = modus_tollens(major, q())
        assert result is None

    def test_hilbert_system_axioms(self):
        """Hilbert-System hat die richtigen Axiome."""
        hs = HilbertSystem()
        assert len(hs.AXIOMS) == 3
        assert "A→(B→A)" in hs.AXIOMS

    def test_hilbert_prove_returns_steps(self):
        """Hilbert-Beweis gibt Schritte zurück."""
        hs = HilbertSystem()
        steps = hs.prove("p→p")
        assert isinstance(steps, list)
        assert len(steps) > 0


# =============================================================================
# 8. MODALLOGIK
# =============================================================================

class TestKripkeFrame:
    """Tests für Kripke-Rahmen."""

    def _simple_frame(self):
        """Einfacher Kripke-Rahmen: w1 → w2, w2 → w1, w1 → w1, w2 → w2."""
        return KripkeFrame(
            ["w1", "w2"],
            {"w1": ["w1", "w2"], "w2": ["w1", "w2"]}
        )

    def test_worlds(self):
        """Welten werden korrekt gespeichert."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"]})
        assert "w1" in frame.worlds
        assert "w2" in frame.worlds

    def test_accessible_from(self):
        """Erreichbare Welten werden korrekt zurückgegeben."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"]})
        assert frame.accessible_from("w1") == ["w2"]

    def test_reflexive_true(self):
        """Reflexiver Rahmen wird erkannt."""
        frame = self._simple_frame()
        assert frame.is_reflexive() is True

    def test_reflexive_false(self):
        """Nicht-reflexiver Rahmen wird erkannt."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"], "w2": []})
        assert frame.is_reflexive() is False

    def test_symmetric_true(self):
        """Symmetrischer Rahmen wird erkannt."""
        frame = self._simple_frame()
        assert frame.is_symmetric() is True

    def test_symmetric_false(self):
        """Nicht-symmetrischer Rahmen wird erkannt."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"], "w2": ["w2"]})
        assert frame.is_symmetric() is False

    def test_transitive(self):
        """Transitiver Rahmen wird erkannt."""
        frame = self._simple_frame()
        assert frame.is_transitive() is True

    def test_serial_true(self):
        """Serieller Rahmen: jede Welt hat Nachfolger."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"], "w2": ["w1"]})
        assert frame.is_serial() is True

    def test_serial_false(self):
        """Nicht-serieller Rahmen: leere Erreichbarkeitsliste."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"], "w2": []})
        assert frame.is_serial() is False

    def test_modal_system_s5(self):
        """Vollständiger Rahmen (alle Welten erreichbar) ist S5."""
        frame = self._simple_frame()
        system = frame.modal_system()
        assert system == "S5"

    def test_modal_system_k(self):
        """Leerer Rahmen ist K."""
        frame = KripkeFrame(["w1"], {"w1": []})
        system = frame.modal_system()
        assert system == "K"


class TestModalFormula:
    """Tests für Modalformeln."""

    def _setup(self):
        """Erstellt einfaches Kripke-Modell."""
        frame = KripkeFrame(
            ["w1", "w2", "w3"],
            {"w1": ["w2", "w3"], "w2": ["w2"], "w3": ["w3"]}
        )
        # p gilt in w2 und w3, nicht in w1
        valuation = {"p": {"w2", "w3"}}
        return frame, valuation

    def test_atom_evaluate(self):
        """Atomare Formel wird korrekt ausgewertet."""
        frame, valuation = self._setup()
        f = ModalFormula.atom("p")
        assert f.evaluate(frame, valuation, "w2") is True
        assert f.evaluate(frame, valuation, "w1") is False

    def test_box_true(self):
        """□p gilt wenn p in allen erreichbaren Welten gilt."""
        frame, valuation = self._setup()
        f = ModalFormula.box(ModalFormula.atom("p"))
        # w1 kann w2 und w3 erreichen, p gilt in beiden → □p gilt in w1
        assert f.evaluate(frame, valuation, "w1") is True

    def test_box_false(self):
        """□p gilt nicht wenn p in einer erreichbaren Welt nicht gilt."""
        frame = KripkeFrame(
            ["w1", "w2"],
            {"w1": ["w1", "w2"], "w2": ["w2"]}
        )
        valuation = {"p": {"w2"}}  # p gilt nur in w2, nicht w1
        f = ModalFormula.box(ModalFormula.atom("p"))
        assert f.evaluate(frame, valuation, "w1") is False

    def test_diamond_true(self):
        """◇p gilt wenn p in mindestens einer erreichbaren Welt gilt."""
        frame, valuation = self._setup()
        f = ModalFormula.diamond(ModalFormula.atom("p"))
        assert f.evaluate(frame, valuation, "w1") is True

    def test_diamond_false(self):
        """◇p gilt nicht wenn p in keiner erreichbaren Welt gilt."""
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w2"], "w2": []})
        valuation = {"p": {"w1"}}  # p gilt nur in w1, w1 kann nur w2 erreichen
        f = ModalFormula.diamond(ModalFormula.atom("p"))
        assert f.evaluate(frame, valuation, "w1") is False

    def test_box_empty_accessible(self):
        """□p gilt trivialerweise wenn keine Welt erreichbar."""
        frame = KripkeFrame(["w1"], {"w1": []})
        valuation = {"p": set()}  # p gilt nirgends
        f = ModalFormula.box(ModalFormula.atom("p"))
        # Vakuumwahrheit: alle 0 erreichbaren Welten erfüllen p
        assert f.evaluate(frame, valuation, "w1") is True

    def test_not_formula(self):
        """NOT-Modalformel wird korrekt ausgewertet."""
        frame, valuation = self._setup()
        f = ModalFormula("NOT", ModalFormula.atom("p"))
        assert f.evaluate(frame, valuation, "w1") is True
        assert f.evaluate(frame, valuation, "w2") is False

    def test_modal_system_description_k(self):
        """modal_system gibt Beschreibung für K zurück."""
        desc = modal_system("K")
        assert "K" in desc
        assert "□" in desc

    def test_modal_system_description_s5(self):
        """modal_system gibt Beschreibung für S5 zurück."""
        desc = modal_system("S5")
        assert "S5" in desc

    def test_modal_system_all_known(self):
        """Alle bekannten Systeme haben eine Beschreibung."""
        for sys in ["K", "D", "T", "B", "S4", "S5"]:
            desc = modal_system(sys)
            assert isinstance(desc, str)
            assert len(desc) > 0

    def test_modal_system_unknown(self):
        """Unbekanntes System gibt Fehlermeldung zurück."""
        desc = modal_system("X99")
        assert "Unbekannt" in desc or "X99" in desc

    def test_str_box(self):
        """Box-Darstellung enthält □."""
        f = ModalFormula.box(ModalFormula.atom("p"))
        assert "□" in str(f)

    def test_str_diamond(self):
        """Diamond-Darstellung enthält ◇."""
        f = ModalFormula.diamond(ModalFormula.atom("p"))
        assert "◇" in str(f)


# =============================================================================
# 9. GÖDEL-NUMMERIERUNG
# =============================================================================

class TestGodelNumbering:
    """Tests für Gödel-Nummerierung."""

    def test_godel_numbering_returns_int(self):
        """Gödel-Nummerierung gibt eine ganze Zahl zurück."""
        n = godel_numbering("p")
        assert isinstance(n, int)
        assert n > 0

    def test_godel_numbering_deterministic(self):
        """Gleiche Formel → gleiche Gödel-Zahl."""
        n1 = godel_numbering("p→q")
        n2 = godel_numbering("p→q")
        assert n1 == n2

    def test_godel_numbering_different_formulas(self):
        """Verschiedene Formeln → verschiedene Gödel-Zahlen."""
        n1 = godel_numbering("p")
        n2 = godel_numbering("q")
        assert n1 != n2

    def test_godel_numbering_empty(self):
        """Leere Formel gibt 1 zurück (leeres Produkt)."""
        n = godel_numbering("")
        assert n == 1

    def test_godel_decode_returns_string(self):
        """Dekodierung gibt einen String zurück."""
        n = godel_numbering("p")
        decoded = godel_decode(n)
        assert isinstance(decoded, str)

    def test_incompleteness_demo_returns_dict(self):
        """Unvollständigkeits-Demo gibt Dict zurück."""
        result = incompleteness_demonstration()
        assert isinstance(result, dict)

    def test_incompleteness_demo_has_keys(self):
        """Demo enthält alle wichtigen Schlüssel."""
        result = incompleteness_demonstration()
        assert "theorem" in result
        assert "statement" in result
        assert "construction" in result
        assert "implications" in result

    def test_incompleteness_demo_construction(self):
        """Konstruktion enthält alle 4 Schritte."""
        result = incompleteness_demonstration()
        construction = result["construction"]
        assert "step1" in construction
        assert "step2" in construction
        assert "step3" in construction
        assert "step4" in construction

    def test_incompleteness_demo_implications(self):
        """Demo enthält mindestens eine Implikation."""
        result = incompleteness_demonstration()
        assert len(result["implications"]) > 0

    def test_incompleteness_demo_second_theorem(self):
        """Zweiter Unvollständigkeitssatz ist enthalten."""
        result = incompleteness_demonstration()
        assert "second_incompleteness" in result


# =============================================================================
# 10. INTEGRATIONSTESTS
# =============================================================================

class TestIntegration:
    """Integrationstests für Zusammenspiel der Komponenten."""

    def test_distributivity_law_cnf(self):
        """Distributivgesetz: p∧(q∨r) ≡ (p∧q)∨(p∧r) via CNF-Semantik."""
        f = LogicFormula("AND", p(), LogicFormula("OR", q(), r()))
        cnf = f.to_cnf()
        for pv in [True, False]:
            for qv in [True, False]:
                for rv in [True, False]:
                    assignment = {"p": pv, "q": qv, "r": rv}
                    assert cnf.evaluate(assignment) == f.evaluate(assignment)

    def test_contrapositive(self):
        """Kontraposition: p→q ≡ ¬q→¬p."""
        f1 = LogicFormula("IMPLIES", p(), q())
        f2 = LogicFormula("IMPLIES", LogicFormula("NOT", q()), LogicFormula("NOT", p()))
        assert logical_equivalence(f1, f2, ["p", "q"]) is True

    def test_hypothetical_syllogism(self):
        """Hypothetischer Syllogismus: p→q, q→r ⊢ p→r."""
        prem1 = LogicFormula("IMPLIES", p(), q())
        prem2 = LogicFormula("IMPLIES", q(), r())
        conc = LogicFormula("IMPLIES", p(), r())
        assert logical_consequence([prem1, prem2], conc, ["p", "q", "r"]) is True

    def test_resolution_and_dpll_agreement(self):
        """Resolution und DPLL stimmen bei unerfüllbarer Formel überein."""
        # p ∧ ¬p
        clauses = [frozenset({"p"}), frozenset({"-p"})]
        assert resolution_refutation(clauses) is True
        assert dpll(clauses) is None

    def test_modal_axiom_t(self):
        """T-Axiom: □p → p gilt in reflexivem Rahmen."""
        # Reflexiver Rahmen: jede Welt ist sich selbst erreichbar
        frame = KripkeFrame(["w1", "w2"], {"w1": ["w1", "w2"], "w2": ["w1", "w2"]})
        valuation = {"p": {"w1"}}  # p gilt nur in w1

        f_box_p = ModalFormula.box(ModalFormula.atom("p"))
        f_p = ModalFormula.atom("p")

        # T-Axiom: wenn □p gilt, dann gilt auch p
        for w in frame.worlds:
            box_val = f_box_p.evaluate(frame, valuation, w)
            p_val = f_p.evaluate(frame, valuation, w)
            # □p → p: wenn □p wahr, dann muss p wahr sein
            if box_val:
                assert p_val is True

    def test_fol_closed_formula(self):
        """∀x.∃y.R(x,y) ist geschlossene Formel."""
        R = Predicate("R", 2)
        x = Term.variable("x")
        y = Term.variable("y")
        inner = FOLFormula("PRED", R, [x, y])
        formula = FOLFormula.forall("x", FOLFormula.exists("y", inner))
        assert formula.is_closed() is True

    def test_cnf_used_in_dpll(self):
        """CNF-Konversion und DPLL arbeiten korrekt zusammen (via Tseitin)."""
        # Erstelle Formel, konvertiere via Tseitin, löse mit DPLL
        f = LogicFormula("AND", p(), q())
        clauses = tseitin_transform(f)
        result = dpll(clauses)
        # Formel ist erfüllbar
        assert result is not None
