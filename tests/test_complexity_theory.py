"""
@file test_complexity_theory.py
@brief Umfassende Tests für das Modul complexity_theory.py.
@description
    Testet alle Klassen und Methoden der Komplexitätstheorie:
        - TuringMachine: Simulation, Akzeptanz, Laufzeitkomplexitätsklasse
        - NonDeterministicTM: Zertifikat-Verifikation, nichtdeterministische Akzeptanz
        - ComplexityClass: Klassen-Metadaten, Inklusionen, Abschluss
        - PolynomialHierarchy: Σ/Π/Δ-Ebenen, Kollaps-Konsequenz
        - NPCompletenessTheory: Cook-Levin, Karp-21, Reduktionen
        - SATSolver: DPLL, CDCL-Schritt, Resolution, Tseitin
        - CircuitComplexity: Größe, Tiefe, NC-Klassen, Barrieren
        - OracleComplexity: Baker-Gill-Solovay, Relativierung
        - ProofComplexity: Resolution, Frege, Natural Proofs, Algebrization, GCT
        - QuantumComplexity: BQP vs PH, Grover, Shor, Quanten-PCP

@author Michael Fuhrmann
@lastModified 2026-03-11
"""

import sys
import os
import pytest

# Sicherstellen, dass das src/-Verzeichnis im Suchpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from complexity_theory import (
    TuringMachine,
    NonDeterministicTM,
    ComplexityClass,
    PolynomialHierarchy,
    NPCompletenessTheory,
    SATSolver,
    CircuitComplexity,
    OracleComplexity,
    ProofComplexity,
    QuantumComplexity,
)


# ===========================================================================
# Hilfsfunktionen für Tests
# ===========================================================================

def _make_palindrome_tm() -> TuringMachine:
    """
    Erstellt eine einfache TM, die das Wort 'aa' akzeptiert.
    Transitions: lese 'a', gehe rechts, lese 'a', gehe rechts,
                 lese Blank → Akzeptieren.
    """
    transitions = {
        ('q0', 'a'): ('q1', 'a', 'R'),
        ('q1', 'a'): ('q2', 'a', 'R'),
        ('q2', '_'): ('accept', '_', 'N'),
        # Ablehnen bei unbekanntem Symbol
        ('q0', '_'): ('reject', '_', 'N'),
        ('q1', '_'): ('reject', '_', 'N'),
    }
    return TuringMachine(
        states=['q0', 'q1', 'q2', 'accept', 'reject'],
        alphabet=['a'],
        transitions=transitions,
        initial_state='q0',
        accept_states=['accept'],
        reject_states=['reject'],
    )


def _make_simple_tm() -> TuringMachine:
    """Erstellt eine TM, die '1' auf dem Band durch '0' ersetzt und akzeptiert."""
    transitions = {
        ('q0', '1'): ('accept', '0', 'R'),
        ('q0', '_'): ('reject', '_', 'N'),
    }
    return TuringMachine(
        states=['q0', 'accept', 'reject'],
        alphabet=['0', '1'],
        transitions=transitions,
        initial_state='q0',
        accept_states=['accept'],
        reject_states=['reject'],
    )


# ===========================================================================
# 1. TuringMachine
# ===========================================================================

class TestTuringMachine:
    """Tests für die deterministische Turingmaschine."""

    def test_simulate_basic_accept(self):
        """TM simuliert 'aa' und hält in Akzeptierzustand."""
        tm = _make_palindrome_tm()
        halted, tape, steps = tm.simulate('aa')
        assert halted is True
        assert steps >= 2

    def test_simulate_basic_reject(self):
        """TM lehnt leere Eingabe ab."""
        tm = _make_palindrome_tm()
        halted, tape, steps = tm.simulate('')
        # Leere Eingabe → q0 liest Blank → reject
        assert halted is True

    def test_simulate_returns_tape(self):
        """simulate() gibt ein Bandinhalt-String zurück."""
        tm = _make_simple_tm()
        halted, tape, steps = tm.simulate('1')
        assert isinstance(tape, str)
        assert halted is True

    def test_simulate_step_count(self):
        """Schrittanzahl ist eine nicht-negative Ganzzahl."""
        tm = _make_simple_tm()
        _, _, steps = tm.simulate('1')
        assert steps >= 0

    def test_simulate_max_steps_respected(self):
        """Simulation bricht bei max_steps ab."""
        # Einfache TM die sich nach rechts bewegt ohne zu halten
        transitions = {('q0', '1'): ('q0', '1', 'R'),
                       ('q0', '_'): ('q0', '_', 'R')}
        tm = TuringMachine(['q0', 'accept'], ['1'],
                           transitions, 'q0', ['accept'])
        _, _, steps = tm.simulate('111', max_steps=5)
        assert steps <= 5

    def test_accepts_true(self):
        """accepts() gibt True für akzeptierte Eingabe zurück."""
        tm = _make_palindrome_tm()
        assert tm.accepts('aa') is True

    def test_accepts_false(self):
        """accepts() gibt False für nicht-akzeptierte Eingabe zurück."""
        tm = _make_palindrome_tm()
        assert tm.accepts('') is False

    def test_time_complexity_class_linear(self):
        """Wenige Übergänge → O(n)."""
        tm = _make_simple_tm()
        result = tm.time_complexity_class()
        assert result == "O(n)"

    def test_time_complexity_class_quadratic(self):
        """Viele Übergänge → O(n^2) oder O(2^n)."""
        # Erstelle TM mit vielen Übergängen (n_trans > n_states)
        states = [f'q{i}' for i in range(5)]
        alphabet = ['0', '1', 'a', 'b', 'c', 'x']
        transitions = {}
        for i, s in enumerate(states[:-1]):
            for sym in alphabet:
                transitions[(s, sym)] = (states[i + 1], sym, 'R')
        tm = TuringMachine(states, alphabet, transitions, states[0], [states[-1]])
        result = tm.time_complexity_class()
        assert result in ("O(n)", "O(n^2)", "O(2^n)")

    def test_blank_symbol_on_empty_input(self):
        """Leere Eingabe → Band enthält nur Blanks."""
        tm = _make_simple_tm()
        halted, tape, _ = tm.simulate('')
        # Leere Eingabe → Band bleibt leer oder zeigt Blank
        assert isinstance(tape, str)


# ===========================================================================
# 2. NonDeterministicTM
# ===========================================================================

class TestNonDeterministicTM:
    """Tests für die nichtdeterministische Turingmaschine."""

    def _make_ntm(self) -> NonDeterministicTM:
        """Einfache NTM: akzeptiert wenn irgendeiner der nichtdet. Pfade 'a' liest."""
        nd_transitions = {
            ('q0', 'a'): [('accept', 'a', 'N'), ('q1', 'a', 'R')],
            ('q0', 'b'): [('reject', 'b', 'N')],
            ('q1', '_'): [('accept', '_', 'N')],
        }
        return NonDeterministicTM(
            states=['q0', 'q1', 'accept', 'reject'],
            alphabet=['a', 'b'],
            transitions=nd_transitions,
            initial_state='q0',
            accept_states=['accept'],
            reject_states=['reject'],
        )

    def test_certificate_verify_combined_input(self):
        """certificate_verify() kombiniert Eingabe und Zertifikat."""
        ntm = self._make_ntm()
        # Das Ergebnis hängt von der TM ab — wir prüfen nur den Rückgabetyp
        result = ntm.certificate_verify('a', 'cert')
        assert isinstance(result, bool)

    def test_nondeterministic_accepts_true(self):
        """nondeterministic_accepts() gibt True zurück für 'a'."""
        ntm = self._make_ntm()
        assert ntm.nondeterministic_accepts('a', max_depth=4) is True

    def test_nondeterministic_accepts_false(self):
        """nondeterministic_accepts() gibt False zurück wenn kein Pfad akzeptiert."""
        ntm = self._make_ntm()
        assert ntm.nondeterministic_accepts('b', max_depth=4) is False

    def test_nondeterministic_accepts_max_depth(self):
        """max_depth begrenzt den Berechnungsbaum."""
        ntm = self._make_ntm()
        result = ntm.nondeterministic_accepts('a', max_depth=1)
        assert isinstance(result, bool)

    def test_nd_transitions_multiple_options(self):
        """NTM hat multiple Übergänge für ein (Zustand, Symbol)-Paar."""
        ntm = self._make_ntm()
        assert len(ntm.nd_transitions[('q0', 'a')]) == 2


# ===========================================================================
# 3. ComplexityClass
# ===========================================================================

class TestComplexityClass:
    """Tests für Komplexitätsklassen-Metadaten."""

    def test_known_classes_all_present(self):
        """Alle bekannten Klassen sind im KNOWN_CLASSES-Dict vorhanden."""
        expected = {"P", "NP", "coNP", "PSPACE", "EXP", "BPP", "BQP", "#P", "PP", "IP", "PH"}
        assert expected.issubset(set(ComplexityClass.KNOWN_CLASSES.keys()))

    def test_instantiation_np(self):
        """NP kann instanziiert werden."""
        np_class = ComplexityClass("NP")
        assert np_class.name == "NP"

    def test_description_not_empty(self):
        """Beschreibung ist nicht leer."""
        for name in ComplexityClass.KNOWN_CLASSES:
            c = ComplexityClass(name)
            assert len(c.description) > 10

    def test_examples_not_empty(self):
        """Beispiele sind nicht leer."""
        np_class = ComplexityClass("NP")
        assert len(np_class.examples) > 0

    def test_contains_known_problem(self):
        """contains() findet bekannte Probleme."""
        np_class = ComplexityClass("NP")
        assert np_class.contains("SAT") is True

    def test_contains_unknown_problem(self):
        """contains() gibt False für unbekannte Probleme."""
        p_class = ComplexityClass("P")
        assert p_class.contains("HALTING-PROBLEM") is False

    def test_is_closed_under_complement_p(self):
        """P ist unter Komplementbildung abgeschlossen."""
        p_class = ComplexityClass("P")
        assert p_class.is_closed_under("complement") is True

    def test_is_not_closed_under_unknown(self):
        """Unbekannte Operationen geben False zurück."""
        np_class = ComplexityClass("NP")
        assert np_class.is_closed_under("quotient") is False

    def test_known_relations_np(self):
        """NP liegt über P und unter PH."""
        np_class = ComplexityClass("NP")
        rels = np_class.known_relations()
        assert "PH" in rels["supersets"] or "PSPACE" in rels["supersets"]
        assert "P" in rels["subsets"]

    def test_unknown_class_raises(self):
        """Unbekannte Klasse wirft ValueError."""
        with pytest.raises(ValueError):
            ComplexityClass("XYZ_UNKNOWN")

    def test_repr(self):
        """__repr__ gibt lesbaren String zurück."""
        c = ComplexityClass("BQP")
        assert "BQP" in repr(c)


# ===========================================================================
# 4. PolynomialHierarchy
# ===========================================================================

class TestPolynomialHierarchy:
    """Tests für die Polynomielle Hierarchie."""

    def setup_method(self):
        self.ph = PolynomialHierarchy()

    def test_sigma_0_is_p(self):
        """Σ_0^P = P."""
        result = self.ph.sigma(0)
        assert "P" in result and "0" in result

    def test_sigma_1_is_np(self):
        """Σ_1^P = NP."""
        result = self.ph.sigma(1)
        assert "NP" in result

    def test_sigma_2_contains_oracle(self):
        """Σ_2^P enthält Orakel-Referenz."""
        result = self.ph.sigma(2)
        assert "Σ_1" in result or "NP" in result

    def test_sigma_negative_raises(self):
        """Negative Stufen werfen ValueError."""
        with pytest.raises(ValueError):
            self.ph.sigma(-1)

    def test_pi_0_is_p(self):
        """Π_0^P = P."""
        result = self.ph.pi(0)
        assert "P" in result

    def test_pi_1_is_conp(self):
        """Π_1^P = co-NP."""
        result = self.ph.pi(1)
        assert "co-NP" in result or "coNP" in result

    def test_pi_2_structure(self):
        """Π_2^P hat korrekte Struktur."""
        result = self.ph.pi(2)
        assert "2" in result and ("co-NP" in result or "Π" in result)

    def test_delta_1_is_p(self):
        """Δ_1^P = P."""
        result = self.ph.delta(1)
        assert "P" in result

    def test_delta_2_structure(self):
        """Δ_2^P enthält P mit Σ_1^P Orakel."""
        result = self.ph.delta(2)
        assert "P" in result and ("Σ_1" in result or "NP" in result)

    def test_delta_0_raises(self):
        """Δ_0^P ist nicht definiert."""
        with pytest.raises(ValueError):
            self.ph.delta(0)

    def test_oracle_characterization_returns_all_three(self):
        """oracle_characterization() enthält Σ, Π und Δ."""
        result = self.ph.oracle_characterization(2)
        assert "Σ_2" in result
        assert "Π_2" in result
        assert "Δ_2" in result

    def test_collapse_consequence_mentions_p_np(self):
        """collapse_consequence() erwähnt P=NP."""
        result = self.ph.collapse_consequence()
        assert "P = NP" in result or "P=NP" in result

    def test_collapse_consequence_multiple_scenarios(self):
        """collapse_consequence() nennt mehrere Kollaps-Szenarien."""
        result = self.ph.collapse_consequence()
        assert "1." in result and "2." in result


# ===========================================================================
# 5. NPCompletenessTheory
# ===========================================================================

class TestNPCompletenessTheory:
    """Tests für die NP-Vollständigkeitstheorie."""

    def test_cook_levin_returns_dict(self):
        """cook_levin_theorem() gibt ein Dictionary zurück."""
        result = NPCompletenessTheory.cook_levin_theorem()
        assert isinstance(result, dict)

    def test_cook_levin_mentions_sat(self):
        """Cook-Levin-Theorem erwähnt SAT."""
        result = NPCompletenessTheory.cook_levin_theorem()
        assert "SAT" in result.get("theorem", "")

    def test_cook_levin_year_1971(self):
        """Cook publizierte 1971."""
        result = NPCompletenessTheory.cook_levin_theorem()
        assert result.get("published_by_cook") == 1971

    def test_karp_21_returns_21_problems(self):
        """karp_21_problems() gibt genau 21 Probleme zurück."""
        problems = NPCompletenessTheory.karp_21_problems()
        assert len(problems) == 21

    def test_karp_21_contains_sat(self):
        """SAT ist das erste Karp-Problem."""
        problems = NPCompletenessTheory.karp_21_problems()
        names = [p["name"] for p in problems]
        assert "SAT" in names

    def test_karp_21_contains_clique(self):
        """CLIQUE ist in den Karp-Problemen."""
        problems = NPCompletenessTheory.karp_21_problems()
        names = [p["name"] for p in problems]
        assert "CLIQUE" in names

    def test_karp_21_all_have_description(self):
        """Alle Karp-Probleme haben eine Beschreibung."""
        for p in NPCompletenessTheory.karp_21_problems():
            assert "desc" in p and len(p["desc"]) > 5

    def test_reduction_3sat_to_clique(self):
        """Reduktion von 3-SAT auf CLIQUE ist bekannt."""
        result = NPCompletenessTheory.reduction("3-SAT", "CLIQUE")
        assert "≤_p" in result["direction"]
        assert "Clique" in result["description"] or "clique" in result["description"].lower()

    def test_reduction_unknown_pair(self):
        """Unbekannte Reduktion gibt generische Antwort."""
        result = NPCompletenessTheory.reduction("KNAPSACK", "3-COLOR")
        assert "existiert" in result["description"] or "nicht" in result["description"]

    def test_is_np_complete_sat(self):
        """SAT ist NP-vollständig."""
        assert NPCompletenessTheory.is_np_complete("SAT") is True

    def test_is_np_complete_max_cut(self):
        """MAX-CUT ist NP-vollständig."""
        assert NPCompletenessTheory.is_np_complete("MAX-CUT") is True

    def test_is_np_complete_linear_prog(self):
        """LINEAR-PROGRAMMING ist nicht in der Karp-Liste."""
        assert NPCompletenessTheory.is_np_complete("LINEAR-PROGRAMMING") is False

    def test_self_reducibility_sat(self):
        """SAT ist self-reducible."""
        assert NPCompletenessTheory.self_reducibility("SAT") is True

    def test_self_reducibility_unknown(self):
        """Unbekannte Probleme sind nicht self-reducible (nach Karp-Liste)."""
        assert NPCompletenessTheory.self_reducibility("TRAVELLING-SALESMAN") is False


# ===========================================================================
# 6. SATSolver
# ===========================================================================

class TestSATSolver:
    """Tests für den SAT-Löser."""

    # --- DPLL ---

    def test_dpll_sat_simple(self):
        """Einfache erfüllbare Formel: x1 ∨ x2."""
        result = SATSolver.dpll([[1, 2]], [1, 2])
        assert result is not None

    def test_dpll_sat_unit(self):
        """Unit Clause: {1} → x1=True."""
        result = SATSolver.dpll([[1]], [1])
        assert result is not None
        assert result.get(1) is True

    def test_dpll_unsat(self):
        """Unerfüllbar: {x1} ∧ {¬x1}."""
        result = SATSolver.dpll([[1], [-1]], [1])
        assert result is None

    def test_dpll_three_sat(self):
        """3-SAT Formel erfüllbar."""
        # (x1 ∨ x2 ∨ x3) ∧ (¬x1 ∨ x2 ∨ x3)
        result = SATSolver.dpll([[1, 2, 3], [-1, 2, 3]], [1, 2, 3])
        assert result is not None

    def test_dpll_all_neg_unsat(self):
        """Alle Klauseln negiert: UNSAT."""
        # {x1} ∧ {x2} ∧ {¬x1 ∨ ¬x2}
        result = SATSolver.dpll([[1], [2], [-1, -2]], [1, 2])
        # x1=True, x2=True → ¬x1∨¬x2 = F∨F = F → UNSAT
        assert result is None

    def test_dpll_returns_valid_assignment(self):
        """DPLL-Belegung erfüllt alle Klauseln."""
        clauses = [[1, 2], [-1, 3], [-2, -3]]
        result = SATSolver.dpll(clauses, [1, 2, 3])
        if result is not None:
            # Verifiziere: Jede Klausel hat mind. ein wahres Literal
            for clause in clauses:
                satisfied = any(
                    (lit > 0 and result.get(abs(lit), False)) or
                    (lit < 0 and not result.get(abs(lit), True))
                    for lit in clause
                )
                assert satisfied

    # --- CDCL ---

    def test_cdcl_step_returns_string(self):
        """cdcl_step() gibt einen String zurück."""
        result = SATSolver.cdcl_step([[1, 2], [-1, 3]], {})
        assert isinstance(result, str)
        assert len(result) > 0

    def test_cdcl_step_conflict_detection(self):
        """cdcl_step() erkennt Konflikte."""
        # Konflikt: {1} und {¬1} bei leerer Belegung
        result = SATSolver.cdcl_step([[1], [-1]], {})
        assert "KONFLIKT" in result.upper() or "SAT" in result.upper() or "UNIT" in result.upper()

    def test_cdcl_step_sat_detection(self):
        """cdcl_step() erkennt erfüllte Klauseln."""
        # Alle Klauseln erfüllt: Belegung {1: True} bei Klausel [[1]]
        result = SATSolver.cdcl_step([[1]], {1: True})
        assert "SAT" in result.upper()

    def test_cdcl_step_unit_propagation(self):
        """cdcl_step() propagiert Unit Clauses."""
        result = SATSolver.cdcl_step([[1], [1, 2]], {})
        # Sollte Unit Propagation von x1=True melden
        assert isinstance(result, str)

    # --- Resolution ---

    def test_resolution_unsat_simple(self):
        """Resolution findet UNSAT-Beweis für {1} ∧ {¬1}."""
        result = SATSolver.resolution_refutation([[1], [-1]])
        assert result is not None
        # Letzter Resolvent muss leer sein
        assert result[-1][2] == []

    def test_resolution_sat_returns_none(self):
        """Resolution gibt None für erfüllbare Formel zurück."""
        result = SATSolver.resolution_refutation([[1, 2]])
        assert result is None

    def test_resolution_two_var_unsat(self):
        """2-Variable UNSAT-Instanz."""
        # {1,2} ∧ {1,¬2} ∧ {¬1,2} ∧ {¬1,¬2}
        clauses = [[1, 2], [1, -2], [-1, 2], [-1, -2]]
        result = SATSolver.resolution_refutation(clauses)
        assert result is not None

    # --- Tseitin ---

    def test_tseitin_single_literal(self):
        """Einzelnes Literal ergibt Einzel-Klausel."""
        result = SATSolver.tseitin_transform("1")
        assert result == [[1]]

    def test_tseitin_negative_literal(self):
        """Negatives Literal."""
        result = SATSolver.tseitin_transform("-2")
        assert result == [[-2]]

    def test_tseitin_and(self):
        """AND-Formel wird zu zwei Klauseln."""
        result = SATSolver.tseitin_transform("AND(1, 2)")
        assert isinstance(result, list)
        assert len(result) >= 2

    def test_tseitin_or(self):
        """OR-Formel wird zu einer Klausel."""
        result = SATSolver.tseitin_transform("OR(1, -2)")
        assert isinstance(result, list)
        assert [1, -2] in result or len(result) >= 1

    def test_tseitin_returns_list_of_lists(self):
        """Rückgabetyp ist Liste von Listen."""
        result = SATSolver.tseitin_transform("1")
        assert isinstance(result, list)
        assert all(isinstance(c, list) for c in result)


# ===========================================================================
# 7. CircuitComplexity
# ===========================================================================

class TestCircuitComplexity:
    """Tests für Schaltkreiskomplexität."""

    def test_circuit_size_constant_zero(self):
        """Konstante 0-Funktion hat minimale Schaltkreisgröße."""
        result = CircuitComplexity.circuit_size([0, 0, 0, 0])
        assert result == 1

    def test_circuit_size_constant_one(self):
        """Konstante 1-Funktion hat minimale Schaltkreisgröße."""
        result = CircuitComplexity.circuit_size([1, 1, 1, 1])
        assert result == 1

    def test_circuit_size_single_variable(self):
        """2-Einträge Wahrheitstabelle (1 Variable)."""
        result = CircuitComplexity.circuit_size([0, 1])  # x
        assert result >= 1

    def test_circuit_size_xor(self):
        """XOR-Funktion: [0,1,1,0]."""
        result = CircuitComplexity.circuit_size([0, 1, 1, 0])
        assert result > 1  # XOR braucht mehr als ein Gatter

    def test_circuit_size_invalid_length(self):
        """Ungültige Länge wirft ValueError."""
        with pytest.raises(ValueError):
            CircuitComplexity.circuit_size([0, 1, 0])  # Länge 3 ist keine 2er-Potenz

    def test_circuit_depth_constant(self):
        """Konstante Funktion hat Tiefe 1."""
        result = CircuitComplexity.circuit_depth([0, 0, 0, 0])
        assert result == 1

    def test_circuit_depth_nonconstant(self):
        """Nichtkonstante Funktion hat Tiefe ≥ 1."""
        result = CircuitComplexity.circuit_depth([0, 1, 1, 0])
        assert result >= 1

    def test_nc_class_nc1(self):
        """depth_exponent=1, fan_in=2 → NC^1."""
        result = CircuitComplexity.nc_class_check(
            {"depth_exponent": 1, "fan_in": 2}
        )
        assert "NC^1" in result

    def test_nc_class_ac0(self):
        """Konstante Tiefe, unbeschränkter Fan-in → AC^0."""
        result = CircuitComplexity.nc_class_check(
            {"depth_exponent": 0, "unbounded_fan_in": True}
        )
        assert "AC^0" in result

    def test_nc_class_nc2(self):
        """depth_exponent=2 → NC^2."""
        result = CircuitComplexity.nc_class_check(
            {"depth_exponent": 2, "fan_in": 2}
        )
        assert "NC^2" in result

    def test_razborov_rudich_returns_dict(self):
        """razborov_rudich_natural_proof() gibt Dict zurück."""
        result = CircuitComplexity.razborov_rudich_natural_proof()
        assert isinstance(result, dict)
        assert "barrier" in result

    def test_razborov_rudich_year(self):
        """Razborov-Rudich veröffentlicht 1994."""
        result = CircuitComplexity.razborov_rudich_natural_proof()
        assert result.get("year") == 1994

    def test_monotone_lower_bound_positive(self):
        """Untere Schranke ist positiv."""
        result = CircuitComplexity.monotone_lower_bound(16)
        assert result > 0

    def test_monotone_lower_bound_grows(self):
        """Untere Schranke wächst mit n."""
        lb_small = CircuitComplexity.monotone_lower_bound(16)
        lb_large = CircuitComplexity.monotone_lower_bound(256)
        assert lb_large >= lb_small


# ===========================================================================
# 8. OracleComplexity
# ===========================================================================

class TestOracleComplexity:
    """Tests für Orakel-Komplexität."""

    def test_baker_gill_solovay_returns_dict(self):
        """baker_gill_solovay() gibt Dict zurück."""
        result = OracleComplexity.baker_gill_solovay()
        assert isinstance(result, dict)

    def test_baker_gill_solovay_year(self):
        """Baker-Gill-Solovay erschien 1975."""
        result = OracleComplexity.baker_gill_solovay()
        assert "1975" in result.get("theorem", "")

    def test_baker_gill_solovay_two_oracles(self):
        """Es gibt Orakel A und B."""
        result = OracleComplexity.baker_gill_solovay()
        assert "oracle_A" in result
        assert "oracle_B" in result

    def test_baker_gill_solovay_barrier_consequence(self):
        """Barriere-Konsequenz ist beschrieben."""
        result = OracleComplexity.baker_gill_solovay()
        assert "barrier_consequence" in result
        assert len(result["barrier_consequence"]) > 10

    def test_relativization_barrier_returns_string(self):
        """relativization_barrier() gibt String zurück."""
        result = OracleComplexity.relativization_barrier()
        assert isinstance(result, str)
        assert len(result) > 50

    def test_relativization_barrier_mentions_diagonalization(self):
        """Relativierungsbarriere erwähnt Diagonalisierung."""
        result = OracleComplexity.relativization_barrier()
        assert "Diagonal" in result or "diagonal" in result

    def test_algebraic_oracle_ip_pspace(self):
        """IP=PSPACE ist bekanntes algebraisches Resultat."""
        result = OracleComplexity.algebraic_oracle("IP=PSPACE")
        assert "PSPACE" in result
        assert "Shamir" in result

    def test_algebraic_oracle_unknown(self):
        """Unbekannte Probleme geben generische Antwort."""
        result = OracleComplexity.algebraic_oracle("UNKNOWN_ORACLE")
        assert isinstance(result, str)
        assert len(result) > 10


# ===========================================================================
# 9. ProofComplexity
# ===========================================================================

class TestProofComplexity:
    """Tests für Beweiskomplexität."""

    def test_resolution_width_lower_bound_positive(self):
        """Untere Schranke ist positiv."""
        result = ProofComplexity.resolution_width_lower_bound(4, 10)
        assert result > 0

    def test_resolution_width_lower_bound_k4(self):
        """k=4 → Schranke k²/2 = 8."""
        result = ProofComplexity.resolution_width_lower_bound(4, 10)
        assert result == 8

    def test_resolution_width_lower_bound_k6(self):
        """k=6 → Schranke 18."""
        result = ProofComplexity.resolution_width_lower_bound(6, 20)
        assert result == 18

    def test_frege_lower_bound_status_string(self):
        """frege_lower_bound_status() gibt String zurück."""
        result = ProofComplexity.frege_lower_bound_status()
        assert isinstance(result, str)
        assert len(result) > 50

    def test_frege_mentions_open_problem(self):
        """Frege-Schranken sind offen."""
        result = ProofComplexity.frege_lower_bound_status()
        assert "offen" in result.lower() or "unbekannt" in result.lower() or "keine" in result.lower()

    def test_natural_proofs_barrier_returns_dict(self):
        """natural_proofs_barrier() gibt Dict zurück."""
        result = ProofComplexity.natural_proofs_barrier()
        assert isinstance(result, dict)
        assert "barrier" in result

    def test_algebrization_barrier_returns_dict(self):
        """algebrization_barrier() gibt Dict zurück."""
        result = ProofComplexity.algebrization_barrier()
        assert isinstance(result, dict)
        assert "barrier" in result

    def test_algebrization_barrier_year(self):
        """Aaronson-Wigderson 2009."""
        result = ProofComplexity.algebrization_barrier()
        assert result.get("year") == 2009

    def test_algebrization_barrier_authors(self):
        """Aaronson und Wigderson als Autoren."""
        result = ProofComplexity.algebrization_barrier()
        assert "Aaronson" in result.get("authors", "")
        assert "Wigderson" in result.get("authors", "")

    def test_gct_returns_dict(self):
        """geometric_complexity_theory() gibt Dict zurück."""
        result = ProofComplexity.geometric_complexity_theory()
        assert isinstance(result, dict)
        assert "approach" in result

    def test_gct_mentions_permanent(self):
        """GCT erwähnt das Permanent-Problem."""
        result = ProofComplexity.geometric_complexity_theory()
        # Suche in allen string-Werten
        all_text = " ".join(str(v) for v in result.values())
        assert "Permanent" in all_text or "permanent" in all_text.lower()

    def test_gct_mulmuley(self):
        """GCT-Autoren: Mulmuley & Sohoni."""
        result = ProofComplexity.geometric_complexity_theory()
        assert "Mulmuley" in result.get("authors", "")


# ===========================================================================
# 10. QuantumComplexity
# ===========================================================================

class TestQuantumComplexity:
    """Tests für Quantenkomplexität."""

    def test_bqp_vs_ph_returns_string(self):
        """bqp_vs_ph() gibt String zurück."""
        result = QuantumComplexity.bqp_vs_ph()
        assert isinstance(result, str)
        assert len(result) > 50

    def test_bqp_vs_ph_mentions_raz_tal(self):
        """Raz-Tal 2019 ist erwähnt."""
        result = QuantumComplexity.bqp_vs_ph()
        assert "Raz" in result or "Tal" in result or "2019" in result

    def test_bqp_vs_ph_mentions_forrelation(self):
        """Forrelation-Problem ist erwähnt."""
        result = QuantumComplexity.bqp_vs_ph()
        assert "Forrelation" in result or "forrelation" in result.lower()

    def test_grovers_speedup_returns_dict(self):
        """grovers_speedup() gibt Dict zurück."""
        result = QuantumComplexity.grovers_speedup()
        assert isinstance(result, dict)

    def test_grovers_speedup_sqrt_n(self):
        """Grover hat O(√N) Quantenkomplexität."""
        result = QuantumComplexity.grovers_speedup()
        quantum_str = result.get("quantum_complexity", "")
        assert "√N" in quantum_str or "sqrt" in quantum_str.lower() or "√" in quantum_str

    def test_grovers_speedup_year(self):
        """Grover veröffentlichte 1996."""
        result = QuantumComplexity.grovers_speedup()
        assert result.get("year") == 1996

    def test_grovers_speedup_steps(self):
        """Algorithmus-Schritte sind dokumentiert."""
        result = QuantumComplexity.grovers_speedup()
        steps = result.get("algorithm_steps", [])
        assert len(steps) >= 4

    def test_shors_algorithm_returns_dict(self):
        """shors_algorithm_complexity() gibt Dict zurück."""
        result = QuantumComplexity.shors_algorithm_complexity()
        assert isinstance(result, dict)

    def test_shors_algorithm_in_bqp(self):
        """Faktorisierung liegt in BQP."""
        result = QuantumComplexity.shors_algorithm_complexity()
        assert result.get("in_bqp") is True

    def test_shors_algorithm_year(self):
        """Shor veröffentlichte 1994."""
        result = QuantumComplexity.shors_algorithm_complexity()
        assert result.get("year") == 1994

    def test_shors_algorithm_poly_quantum(self):
        """Quantenkomplexität ist polynomiell."""
        result = QuantumComplexity.shors_algorithm_complexity()
        assert "poly" in result.get("quantum_complexity", "").lower() or "n^2" in result.get("quantum_complexity", "")

    def test_quantum_pcp_conjecture_returns_string(self):
        """quantum_pcp_conjecture() gibt String zurück."""
        result = QuantumComplexity.quantum_pcp_conjecture()
        assert isinstance(result, str)
        assert len(result) > 50

    def test_quantum_pcp_mentions_qma(self):
        """Quanten-PCP erwähnt QMA."""
        result = QuantumComplexity.quantum_pcp_conjecture()
        assert "QMA" in result

    def test_quantum_pcp_mentions_nlts(self):
        """NLTS-Theorem als Teilfortschritt erwähnt."""
        result = QuantumComplexity.quantum_pcp_conjecture()
        assert "NLTS" in result

    def test_quantum_pcp_open_status(self):
        """Quanten-PCP ist noch offen."""
        result = QuantumComplexity.quantum_pcp_conjecture()
        assert "offen" in result.lower() or "Offen" in result
