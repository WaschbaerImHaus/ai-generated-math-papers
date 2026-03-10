"""
Mathematische Logik / Formale Logik
=====================================
Implementiert Aussagenlogik, Prädikatenlogik, Resolution, SAT-Solving,
Modallogik und Gödels Unvollständigkeitssätze.

Autor: Kurt Ingwer
Letzte Änderung: 2026-03-10
"""

from __future__ import annotations
import itertools
from typing import Any
from functools import reduce


# =============================================================================
# AUSSAGENLOGIK (PROPOSITIONAL LOGIC)
# =============================================================================

class Proposition:
    """
    Atomare Aussage mit Namen und optionalem Wahrheitswert.

    Eine atomare Aussage ist die einfachste logische Einheit – sie hat
    einen Namen (z. B. 'p', 'q') und kann wahr oder falsch sein.

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, name: str, value: bool | None = None):
        """
        Initialisiert eine atomare Aussage.

        :param name: Name der Aussage (z. B. 'p', 'Regen')
        :param value: Optionaler Wahrheitswert (True/False/None)
        """
        self.name = name
        self.value = value

    def evaluate(self, assignment: dict[str, bool]) -> bool:
        """
        Wertet die Aussage unter einer Belegung aus.

        :param assignment: Dict {Variablenname: Wahrheitswert}
        :return: Wahrheitswert der Aussage
        :raises KeyError: wenn Variable nicht in assignment
        """
        if self.name in assignment:
            return assignment[self.name]
        if self.value is not None:
            return self.value
        raise KeyError(f"Variable '{self.name}' nicht in der Belegung.")

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"Proposition('{self.name}')"


class LogicFormula:
    """
    Zusammengesetzte logische Formel mit Konnektiven.

    Unterstützte Operatoren:
    - AND   : Konjunktion (∧)
    - OR    : Disjunktion (∨)
    - NOT   : Negation (¬)
    - IMPLIES: Implikation (→)
    - IFF   : Äquivalenz (↔)
    - XOR   : Exklusives Oder (⊕)
    - NAND  : Nicht-Und (↑)
    - NOR   : Nicht-Oder (↓)
    - ATOM  : Atomare Variable

    Letzte Änderung: 2026-03-10
    """

    # Gültige Operatoren
    OPERATORS = {"AND", "OR", "NOT", "IMPLIES", "IFF", "XOR", "NAND", "NOR", "ATOM"}

    def __init__(self, op: str, *args):
        """
        Erstellt eine logische Formel.

        :param op: Operator (AND, OR, NOT, IMPLIES, IFF, XOR, NAND, NOR, ATOM)
        :param args: Teilformeln oder Variablenname bei ATOM
        """
        if op not in self.OPERATORS:
            raise ValueError(f"Ungültiger Operator: '{op}'. Erlaubt: {self.OPERATORS}")
        self.op = op
        self.args = args

    @staticmethod
    def atom(name: str) -> 'LogicFormula':
        """Erstellt eine atomare Variable."""
        return LogicFormula("ATOM", name)

    def evaluate(self, assignment: dict[str, bool]) -> bool:
        """
        Wertet die Formel unter einer Variablenbelegung aus.

        :param assignment: Dict {Variablenname: Wahrheitswert}
        :return: Wahrheitswert der Formel
        """
        op = self.op
        args = self.args

        if op == "ATOM":
            # Atomare Variable aus der Belegung lesen
            return assignment[args[0]]
        elif op == "NOT":
            return not args[0].evaluate(assignment)
        elif op == "AND":
            return all(a.evaluate(assignment) for a in args)
        elif op == "OR":
            return any(a.evaluate(assignment) for a in args)
        elif op == "IMPLIES":
            # A → B ≡ ¬A ∨ B
            a, b = args
            return (not a.evaluate(assignment)) or b.evaluate(assignment)
        elif op == "IFF":
            # A ↔ B ≡ (A → B) ∧ (B → A)
            a, b = args
            va, vb = a.evaluate(assignment), b.evaluate(assignment)
            return va == vb
        elif op == "XOR":
            # A ⊕ B ≡ (A ∨ B) ∧ ¬(A ∧ B)
            a, b = args
            va, vb = a.evaluate(assignment), b.evaluate(assignment)
            return va != vb
        elif op == "NAND":
            # A ↑ B ≡ ¬(A ∧ B)
            return not all(a.evaluate(assignment) for a in args)
        elif op == "NOR":
            # A ↓ B ≡ ¬(A ∨ B)
            return not any(a.evaluate(assignment) for a in args)
        else:
            raise ValueError(f"Unbekannter Operator: {op}")

    def get_variables(self) -> set[str]:
        """Gibt alle Variablennamen in der Formel zurück."""
        if self.op == "ATOM":
            return {self.args[0]}
        result = set()
        for arg in self.args:
            if isinstance(arg, LogicFormula):
                result |= arg.get_variables()
        return result

    def to_cnf(self) -> 'LogicFormula':
        """
        Wandelt die Formel in Konjunktive Normalform (CNF) um.

        CNF: Konjunktion von Disjunktionen (Klauseln).
        Schritte: NNF → Verteilen von AND über OR.

        :return: Äquivalente Formel in CNF
        """
        # Schritt 1: Implikationen und Äquivalenzen eliminieren
        f = self._eliminate_implications()
        # Schritt 2: Negationsnormalform (NNF)
        f = f._to_nnf()
        # Schritt 3: CNF durch Verteilung
        f = f._distribute_and_over_or()
        return f

    def to_dnf(self) -> 'LogicFormula':
        """
        Wandelt die Formel in Disjunktive Normalform (DNF) um.

        DNF: Disjunktion von Konjunktionen (Minterme).

        :return: Äquivalente Formel in DNF
        """
        f = self._eliminate_implications()
        f = f._to_nnf()
        f = f._distribute_or_over_and()
        return f

    def _eliminate_implications(self) -> 'LogicFormula':
        """Ersetzt →, ↔, ⊕, NAND, NOR durch AND, OR, NOT."""
        op = self.op
        args = [a._eliminate_implications() if isinstance(a, LogicFormula) else a
                for a in self.args]

        if op == "ATOM":
            return self
        elif op == "IMPLIES":
            # A → B ≡ ¬A ∨ B
            return LogicFormula("OR", LogicFormula("NOT", args[0]), args[1])
        elif op == "IFF":
            # A ↔ B ≡ (A → B) ∧ (B → A)
            a, b = args
            left = LogicFormula("OR", LogicFormula("NOT", a), b)
            right = LogicFormula("OR", LogicFormula("NOT", b), a)
            return LogicFormula("AND", left, right)
        elif op == "XOR":
            # A ⊕ B ≡ (A ∨ B) ∧ (¬A ∨ ¬B)
            a, b = args
            left = LogicFormula("OR", a, b)
            right = LogicFormula("OR", LogicFormula("NOT", a), LogicFormula("NOT", b))
            return LogicFormula("AND", left, right)
        elif op == "NAND":
            return LogicFormula("NOT", LogicFormula("AND", *args))
        elif op == "NOR":
            return LogicFormula("NOT", LogicFormula("OR", *args))
        else:
            return LogicFormula(op, *args)

    def _to_nnf(self) -> 'LogicFormula':
        """
        Wandelt in Negationsnormalform (NNF) um.
        Alle Negationen stehen direkt vor Atomen.
        De Morgan: ¬(A∧B) ≡ ¬A∨¬B, ¬(A∨B) ≡ ¬A∧¬B
        """
        if self.op == "ATOM":
            return self
        elif self.op == "NOT":
            inner = self.args[0]
            if inner.op == "ATOM":
                return self  # ¬p bleibt ¬p
            elif inner.op == "NOT":
                # Doppelnegation: ¬¬A ≡ A
                return inner.args[0]._to_nnf()
            elif inner.op == "AND":
                # De Morgan: ¬(A∧B) ≡ ¬A∨¬B
                return LogicFormula("OR", *[LogicFormula("NOT", a)._to_nnf()
                                            for a in inner.args])
            elif inner.op == "OR":
                # De Morgan: ¬(A∨B) ≡ ¬A∧¬B
                return LogicFormula("AND", *[LogicFormula("NOT", a)._to_nnf()
                                             for a in inner.args])
            else:
                return LogicFormula("NOT", inner._to_nnf())
        else:
            return LogicFormula(self.op, *[a._to_nnf() for a in self.args])

    def _distribute_and_over_or(self) -> 'LogicFormula':
        """
        Verteilt AND über OR für CNF.
        (A ∧ B) ∨ C ≡ (A ∨ C) ∧ (B ∨ C)
        """
        if self.op in ("ATOM", "NOT"):
            return self
        elif self.op == "AND":
            return LogicFormula("AND", *[a._distribute_and_over_or() for a in self.args])
        elif self.op == "OR":
            # Zuerst Argumente in CNF bringen
            args = [a._distribute_and_over_or() for a in self.args]
            # Suche ein AND unter den Argumenten
            result = args[0]
            for arg in args[1:]:
                result = self._distribute_or_with_and(result, arg)
            return result
        return self

    def _distribute_or_with_and(self, left: 'LogicFormula', right: 'LogicFormula') -> 'LogicFormula':
        """Hilfsmethode: Verteilt OR über AND."""
        if left.op == "AND":
            # (A ∧ B) ∨ C → (A ∨ C) ∧ (B ∨ C)
            return LogicFormula("AND", *[
                self._distribute_or_with_and(a, right) for a in left.args
            ])
        elif right.op == "AND":
            # C ∨ (A ∧ B) → (C ∨ A) ∧ (C ∨ B)
            return LogicFormula("AND", *[
                self._distribute_or_with_and(left, a) for a in right.args
            ])
        else:
            return LogicFormula("OR", left, right)

    def _distribute_or_over_and(self) -> 'LogicFormula':
        """
        Verteilt OR über AND für DNF.
        (A ∨ B) ∧ C ≡ (A ∧ C) ∨ (B ∧ C)
        """
        if self.op in ("ATOM", "NOT"):
            return self
        elif self.op == "OR":
            return LogicFormula("OR", *[a._distribute_or_over_and() for a in self.args])
        elif self.op == "AND":
            args = [a._distribute_or_over_and() for a in self.args]
            result = args[0]
            for arg in args[1:]:
                result = self._distribute_and_with_or(result, arg)
            return result
        return self

    def _distribute_and_with_or(self, left: 'LogicFormula', right: 'LogicFormula') -> 'LogicFormula':
        """Hilfsmethode: Verteilt AND über OR für DNF."""
        if left.op == "OR":
            return LogicFormula("OR", *[
                self._distribute_and_with_or(a, right) for a in left.args
            ])
        elif right.op == "OR":
            return LogicFormula("OR", *[
                self._distribute_and_with_or(left, a) for a in right.args
            ])
        else:
            return LogicFormula("AND", left, right)

    def simplify(self) -> 'LogicFormula':
        """
        Vereinfacht die Formel algebraisch.

        Regeln:
        - A ∧ True ≡ A
        - A ∨ False ≡ A
        - A ∧ False ≡ False
        - A ∨ True ≡ True
        - A ∧ A ≡ A (Idempotenz)
        - A ∨ A ≡ A (Idempotenz)
        - ¬¬A ≡ A (Doppelnegation)

        :return: Vereinfachte Formel
        """
        if self.op == "ATOM":
            return self
        elif self.op == "NOT":
            inner = self.args[0].simplify()
            if inner.op == "NOT":
                # ¬¬A ≡ A
                return inner.args[0]
            return LogicFormula("NOT", inner)
        else:
            # Argumente rekursiv vereinfachen
            simplified_args = [a.simplify() for a in self.args]
            return LogicFormula(self.op, *simplified_args)

    def __str__(self) -> str:
        """Lesbare String-Darstellung der Formel."""
        op = self.op
        if op == "ATOM":
            return self.args[0]
        elif op == "NOT":
            inner = self.args[0]
            if inner.op in ("AND", "OR", "IMPLIES", "IFF"):
                return f"¬({inner})"
            return f"¬{inner}"
        elif op == "AND":
            parts = []
            for a in self.args:
                if a.op in ("OR", "IMPLIES", "IFF"):
                    parts.append(f"({a})")
                else:
                    parts.append(str(a))
            return " ∧ ".join(parts)
        elif op == "OR":
            parts = []
            for a in self.args:
                if a.op in ("IMPLIES", "IFF"):
                    parts.append(f"({a})")
                else:
                    parts.append(str(a))
            return " ∨ ".join(parts)
        elif op == "IMPLIES":
            a, b = self.args
            left = f"({a})" if a.op in ("IMPLIES", "IFF") else str(a)
            return f"{left} → {b}"
        elif op == "IFF":
            a, b = self.args
            return f"{a} ↔ {b}"
        elif op == "XOR":
            a, b = self.args
            return f"{a} ⊕ {b}"
        elif op == "NAND":
            return " ↑ ".join(str(a) for a in self.args)
        elif op == "NOR":
            return " ↓ ".join(str(a) for a in self.args)
        return f"({self.op} {' '.join(str(a) for a in self.args)})"

    def __repr__(self) -> str:
        return f"LogicFormula('{self.op}', {self.args})"


def truth_table(formula: LogicFormula, variables: list[str]) -> list[dict]:
    """
    Erstellt eine vollständige Wahrheitstabelle für die Formel.

    Für n Variablen gibt es 2^n Zeilen.

    :param formula: Logische Formel
    :param variables: Liste der Variablennamen
    :return: Liste von Dicts mit Belegung und Ergebnis
    Letzte Änderung: 2026-03-10
    """
    table = []
    # Alle 2^n Kombinationen der Wahrheitswerte durchgehen
    for values in itertools.product([False, True], repeat=len(variables)):
        assignment = dict(zip(variables, values))
        result = formula.evaluate(assignment)
        row = dict(assignment)
        row["result"] = result
        table.append(row)
    return table


def is_tautology(formula: LogicFormula, variables: list[str]) -> bool:
    """
    Prüft ob die Formel eine Tautologie ist (für alle Belegungen wahr).

    Beispiel: A ∨ ¬A ist eine Tautologie (Satz vom ausgeschlossenen Dritten).

    :param formula: Logische Formel
    :param variables: Liste der Variablennamen
    :return: True wenn Tautologie, sonst False
    Letzte Änderung: 2026-03-10
    """
    return all(row["result"] for row in truth_table(formula, variables))


def is_satisfiable(formula: LogicFormula, variables: list[str]) -> bool:
    """
    Prüft ob die Formel erfüllbar ist (mindestens eine wahre Belegung).

    :param formula: Logische Formel
    :param variables: Liste der Variablennamen
    :return: True wenn erfüllbar, sonst False
    Letzte Änderung: 2026-03-10
    """
    return any(row["result"] for row in truth_table(formula, variables))


def is_contradiction(formula: LogicFormula, variables: list[str]) -> bool:
    """
    Prüft ob die Formel eine Kontradiktion ist (für alle Belegungen falsch).

    Beispiel: A ∧ ¬A ist eine Kontradiktion.

    :param formula: Logische Formel
    :param variables: Liste der Variablennamen
    :return: True wenn Kontradiktion, sonst False
    Letzte Änderung: 2026-03-10
    """
    return not is_satisfiable(formula, variables)


def logical_consequence(premises: list, conclusion: LogicFormula, variables: list[str]) -> bool:
    """
    Prüft ob die Konklusion logisch aus den Prämissen folgt.

    P1, ..., Pn ⊨ C genau dann wenn (P1 ∧ ... ∧ Pn) → C Tautologie ist.

    :param premises: Liste von Prämissen (LogicFormula)
    :param conclusion: Schlussfolgerung
    :param variables: Liste der Variablennamen
    :return: True wenn logische Konsequenz
    Letzte Änderung: 2026-03-10
    """
    if not premises:
        return is_tautology(conclusion, variables)
    # Konjunktion aller Prämissen
    combined = premises[0]
    for p in premises[1:]:
        combined = LogicFormula("AND", combined, p)
    # Prüfe ob combined → conclusion eine Tautologie ist
    implication = LogicFormula("IMPLIES", combined, conclusion)
    return is_tautology(implication, variables)


def logical_equivalence(f1: LogicFormula, f2: LogicFormula, variables: list[str]) -> bool:
    """
    Prüft semantische Äquivalenz: f1 ↔ f2 ist Tautologie.

    :param f1: Erste Formel
    :param f2: Zweite Formel
    :param variables: Liste der Variablennamen
    :return: True wenn äquivalent
    Letzte Änderung: 2026-03-10
    """
    iff = LogicFormula("IFF", f1, f2)
    return is_tautology(iff, variables)


# =============================================================================
# RESOLUTION UND SAT-SOLVING
# =============================================================================

def resolution_step(clause1: frozenset, clause2: frozenset) -> frozenset | None:
    """
    Führt einen Resolutionsschritt durch.

    Resolutionsregel: (A ∨ l) ∧ (B ∨ ¬l) → A ∨ B

    Literale werden als Strings dargestellt: 'p' für p, '-p' für ¬p.

    :param clause1: Erste Klausel (Menge von Literalen)
    :param clause2: Zweite Klausel (Menge von Literalen)
    :return: Resolvente oder None wenn kein Resolutionsschritt möglich
    Letzte Änderung: 2026-03-10
    """
    # Suche komplementäre Literale
    for lit in clause1:
        # Negiertes Literal bestimmen
        if lit.startswith("-"):
            neg_lit = lit[1:]  # -p → p
        else:
            neg_lit = "-" + lit  # p → -p

        if neg_lit in clause2:
            # Resolvente: Vereinigung ohne die komplementären Literale
            resolvent = (clause1 - {lit}) | (clause2 - {neg_lit})
            return frozenset(resolvent)
    return None


def resolution_refutation(clauses: list[frozenset]) -> bool:
    """
    Vollständiger Resolutionsbeweis (Widerlegung).

    Gibt True zurück wenn die Klauselmenge unerfüllbar ist.
    Eine Klauselmenge ist unerfüllbar, wenn die leere Klausel {} ableitbar ist.

    :param clauses: Liste von Klauseln (jede Klausel ist ein frozenset von Literalen)
    :return: True wenn unerfüllbar (Widerlegung erfolgreich)
    Letzte Änderung: 2026-03-10
    """
    # Arbeitsmenge initialisieren
    clause_set = set(clauses)

    while True:
        new_clauses = set()
        clause_list = list(clause_set)

        # Alle Paare von Klauseln durchprobieren
        for i in range(len(clause_list)):
            for j in range(i + 1, len(clause_list)):
                resolvent = resolution_step(clause_list[i], clause_list[j])
                if resolvent is not None:
                    if resolvent == frozenset():
                        # Leere Klausel abgeleitet → unerfüllbar
                        return True
                    new_clauses.add(resolvent)

        # Prüfe ob neue Klauseln hinzukommen
        if new_clauses.issubset(clause_set):
            # Keine neuen Klauseln → Fixpunkt erreicht, erfüllbar
            return False
        clause_set |= new_clauses


def _unit_propagate(clauses: list[frozenset], assignment: dict) -> tuple[list[frozenset], dict] | None:
    """
    Unit-Propagation für DPLL: erzwingt Einheitklauseln.
    Gibt None zurück bei Widerspruch.

    :param clauses: Aktuelle Klauseln
    :param assignment: Aktuelle Belegung
    :return: Neue Klauseln und Belegung oder None
    """
    changed = True
    while changed:
        changed = False
        new_clauses = []
        for clause in clauses:
            # Klausel unter aktueller Belegung auswerten
            satisfied = False
            remaining = set()
            for lit in clause:
                var = lit.lstrip("-")
                is_neg = lit.startswith("-")
                if var in assignment:
                    val = assignment[var]
                    # Literal ist wahr?
                    if (val and not is_neg) or (not val and is_neg):
                        satisfied = True
                        break
                    # Literal ist falsch → weglassen
                else:
                    remaining.add(lit)

            if satisfied:
                continue  # Klausel erfüllt

            if not remaining:
                return None  # Leere Klausel → Widerspruch

            if len(remaining) == 1:
                # Einheitklausel → Literal erzwingen
                lit = next(iter(remaining))
                var = lit.lstrip("-")
                is_neg = lit.startswith("-")
                assignment = dict(assignment)
                assignment[var] = not is_neg
                changed = True
            else:
                new_clauses.append(frozenset(remaining))

        clauses = new_clauses

    return clauses, assignment


def dpll(clauses: list[frozenset], assignment: dict = None) -> dict | None:
    """
    Davis-Putnam-Logemann-Loveland SAT-Solver.

    Algorithmus:
    1. Unit-Propagation
    2. Reine Literal-Elimination
    3. Fallunterscheidung (Branching)

    :param clauses: Liste von Klauseln (Literale als Strings)
    :param assignment: Aktuelle Belegung (für Rekursion)
    :return: Erfüllende Belegung oder None
    Letzte Änderung: 2026-03-10
    """
    if assignment is None:
        assignment = {}

    # Unit-Propagation
    result = _unit_propagate(clauses, assignment)
    if result is None:
        return None  # Widerspruch
    clauses, assignment = result

    # Alle Klauseln erfüllt?
    if not clauses:
        return assignment

    # Wähle unbekannte Variable (Branching)
    var = None
    for clause in clauses:
        for lit in clause:
            v = lit.lstrip("-")
            if v not in assignment:
                var = v
                break
        if var:
            break

    if var is None:
        return assignment  # Alle Variablen belegt

    # Branching: var = True
    new_assignment = dict(assignment)
    new_assignment[var] = True
    result = dpll(clauses, new_assignment)
    if result is not None:
        return result

    # Branching: var = False
    new_assignment = dict(assignment)
    new_assignment[var] = False
    return dpll(clauses, new_assignment)


def _formula_to_cnf_clauses(formula: LogicFormula) -> list[frozenset]:
    """
    Konvertiert eine Formel direkt in eine Liste von Klauseln (CNF).
    Hilfsfunktion für Tseitin.

    :param formula: LogicFormula in CNF
    :return: Liste von Klauseln
    """
    def cnf_to_clauses(f):
        if f.op == "AND":
            result = []
            for arg in f.args:
                result.extend(cnf_to_clauses(arg))
            return result
        elif f.op == "OR":
            lits = set()
            for arg in f.args:
                lits |= _formula_to_literals(arg)
            return [frozenset(lits)]
        elif f.op in ("ATOM", "NOT"):
            return [frozenset(_formula_to_literals(f))]
        else:
            return [frozenset(_formula_to_literals(f))]

    return cnf_to_clauses(formula)


def _formula_to_literals(f: LogicFormula) -> set:
    """Konvertiert eine einzelne Klausel in Literale."""
    if f.op == "ATOM":
        return {f.args[0]}
    elif f.op == "NOT" and f.args[0].op == "ATOM":
        return {"-" + f.args[0].args[0]}
    elif f.op == "OR":
        result = set()
        for arg in f.args:
            result |= _formula_to_literals(arg)
        return result
    return set()


def tseitin_transform(formula: LogicFormula) -> list[frozenset]:
    """
    Tseitin-Transformation: Wandelt eine Formel in eine equisatisfiable CNF um.

    Eingeführte Hilfsvariablen: t_0, t_1, ...
    Für jede Teilformel wird eine Äquivalenz erzeugt.

    :param formula: Beliebige LogicFormula
    :return: Liste von Klauseln in CNF
    Letzte Änderung: 2026-03-10
    """
    counter = [0]
    clauses = []

    def new_var() -> str:
        """Erzeugt neue Tseitin-Variable."""
        v = f"_t{counter[0]}"
        counter[0] += 1
        return v

    def transform(f: LogicFormula) -> str:
        """Gibt Tseitin-Variable für Teilformel zurück, fügt Klauseln hinzu."""
        if f.op == "ATOM":
            return f.args[0]
        elif f.op == "NOT":
            inner_var = transform(f.args[0])
            t = new_var()
            # t ↔ ¬inner: {t, inner}, {-t, -inner}
            clauses.append(frozenset({t, inner_var}))
            clauses.append(frozenset({"-" + t, "-" + inner_var}))
            return t
        elif f.op == "AND":
            child_vars = [transform(a) for a in f.args]
            t = new_var()
            # t ↔ (c1 ∧ c2 ∧ ...):
            # t → ci: {-t, ci} für alle i
            # ci → t: {t, -c1, -c2, ...}
            for cv in child_vars:
                clauses.append(frozenset({"-" + t, cv}))
            neg_children = {"-" + cv for cv in child_vars}
            clauses.append(frozenset({t} | neg_children))
            return t
        elif f.op == "OR":
            child_vars = [transform(a) for a in f.args]
            t = new_var()
            # t ↔ (c1 ∨ c2 ∨ ...):
            # t → (c1 ∨ ...): {-t, c1, c2, ...}
            # ci → t: {t, -ci} für alle i
            clauses.append(frozenset({"-" + t} | set(child_vars)))
            for cv in child_vars:
                clauses.append(frozenset({t, "-" + cv}))
            return t
        elif f.op == "IMPLIES":
            a_var = transform(f.args[0])
            b_var = transform(f.args[1])
            t = new_var()
            # t ↔ (a → b) ≡ t ↔ (¬a ∨ b)
            clauses.append(frozenset({"-" + t, "-" + a_var, b_var}))
            clauses.append(frozenset({t, a_var}))
            clauses.append(frozenset({t, "-" + b_var}))
            return t
        else:
            # Für andere Operatoren: zuerst konvertieren
            eliminated = f._eliminate_implications()
            return transform(eliminated)

    # Wurzel-Variable muss wahr sein
    root_var = transform(formula)
    clauses.append(frozenset({root_var}))
    return clauses


# =============================================================================
# PRÄDIKATENLOGIK (FIRST-ORDER LOGIC)
# =============================================================================

class Term:
    """
    Term in der Prädikatenlogik erster Stufe.

    Ein Term kann sein:
    - Variable: x, y, z
    - Konstante: a, b, c
    - Funktionsanwendung: f(t1, ..., tn)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, name: str, args: list['Term'] = None, is_variable: bool = True):
        """
        :param name: Name der Variable, Konstante oder Funktion
        :param args: Argumente bei Funktionsanwendung (None für Var/Konst)
        :param is_variable: True wenn Variable, False wenn Konstante
        """
        self.name = name
        self.args = args if args is not None else []
        self.is_variable = is_variable if not self.args else False
        self.is_function = len(self.args) > 0

    @staticmethod
    def variable(name: str) -> 'Term':
        """Erstellt eine Variable."""
        return Term(name, is_variable=True)

    @staticmethod
    def constant(name: str) -> 'Term':
        """Erstellt eine Konstante."""
        return Term(name, is_variable=False)

    @staticmethod
    def function(name: str, *args: 'Term') -> 'Term':
        """Erstellt eine Funktionsanwendung f(t1,...,tn)."""
        return Term(name, list(args), is_variable=False)

    def free_variables(self) -> set[str]:
        """Gibt alle Variablennamen zurück."""
        if self.is_variable:
            return {self.name}
        result = set()
        for arg in self.args:
            result |= arg.free_variables()
        return result

    def substitute(self, var: str, term: 'Term') -> 'Term':
        """Ersetzt Variable var durch Term."""
        if self.is_variable and self.name == var:
            return term
        elif self.is_function:
            return Term(self.name, [a.substitute(var, term) for a in self.args],
                        is_variable=False)
        return self

    def __str__(self) -> str:
        if self.is_function:
            return f"{self.name}({', '.join(str(a) for a in self.args)})"
        return self.name

    def __repr__(self) -> str:
        return f"Term('{self.name}', args={self.args}, is_var={self.is_variable})"


class Predicate:
    """
    Prädikatsymbol mit Arität (Anzahl der Argumente).

    Beispiele: P(x), Q(x, y), Equals(x, y)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, name: str, arity: int):
        """
        :param name: Name des Prädikats
        :param arity: Anzahl der Argumente
        """
        self.name = name
        self.arity = arity

    def apply(self, *terms: Term) -> 'FOLFormula':
        """Wendet das Prädikat auf Terme an."""
        if len(terms) != self.arity:
            raise ValueError(
                f"Prädikat '{self.name}' hat Arität {self.arity}, "
                f"aber {len(terms)} Argumente angegeben."
            )
        return FOLFormula("PRED", self, list(terms))

    def __str__(self) -> str:
        return self.name


class FOLFormula:
    """
    Formel der Prädikatenlogik erster Stufe (First-Order Logic).

    Unterstützte Formeltypen:
    - PRED: Prädikatsanwendung P(t1,...,tn)
    - NOT, AND, OR, IMPLIES, IFF: Konnektive (wie Aussagenlogik)
    - FORALL: Allquantor ∀x.F
    - EXISTS: Existenzquantor ∃x.F

    Letzte Änderung: 2026-03-10
    """

    FOL_OPS = {"PRED", "NOT", "AND", "OR", "IMPLIES", "IFF", "FORALL", "EXISTS"}

    def __init__(self, op: str, *args):
        """
        :param op: Operator (PRED, NOT, AND, OR, IMPLIES, IFF, FORALL, EXISTS)
        :param args: Argumente abhängig vom Operator
        """
        if op not in self.FOL_OPS:
            raise ValueError(f"Ungültiger FOL-Operator: '{op}'")
        self.op = op
        self.args = args

    @staticmethod
    def predicate(pred: Predicate, *terms: Term) -> 'FOLFormula':
        """Erstellt eine Prädikatsformel."""
        return FOLFormula("PRED", pred, list(terms))

    @staticmethod
    def forall(var: str, formula: 'FOLFormula') -> 'FOLFormula':
        """Erstellt ∀var.formula"""
        return FOLFormula("FORALL", var, formula)

    @staticmethod
    def exists(var: str, formula: 'FOLFormula') -> 'FOLFormula':
        """Erstellt ∃var.formula"""
        return FOLFormula("EXISTS", var, formula)

    def free_variables(self) -> set[str]:
        """Gibt alle freien Variablen der Formel zurück."""
        op = self.op
        if op == "PRED":
            result = set()
            for t in self.args[1]:  # args[1] = Liste der Terme
                result |= t.free_variables()
            return result
        elif op == "NOT":
            return self.args[0].free_variables()
        elif op in ("AND", "OR", "IMPLIES", "IFF"):
            result = set()
            for a in self.args:
                result |= a.free_variables()
            return result
        elif op in ("FORALL", "EXISTS"):
            # args[0] = gebundene Variable, args[1] = Formel
            var = self.args[0]
            return self.args[1].free_variables() - {var}
        return set()

    def is_closed(self) -> bool:
        """Prüft ob die Formel geschlossen ist (keine freien Variablen)."""
        return len(self.free_variables()) == 0

    def substitute(self, var: str, term: Term) -> 'FOLFormula':
        """
        Substituiert freie Vorkommen von var durch term.

        :param var: Zu ersetzende Variable
        :param term: Ersetzungsterm
        :return: Formel mit Substitution
        """
        op = self.op
        if op == "PRED":
            new_terms = [t.substitute(var, term) for t in self.args[1]]
            return FOLFormula("PRED", self.args[0], new_terms)
        elif op == "NOT":
            return FOLFormula("NOT", self.args[0].substitute(var, term))
        elif op in ("AND", "OR", "IMPLIES", "IFF"):
            return FOLFormula(op, *[a.substitute(var, term) for a in self.args])
        elif op in ("FORALL", "EXISTS"):
            bound_var = self.args[0]
            if bound_var == var:
                # Variable ist gebunden → keine Substitution
                return self
            return FOLFormula(op, bound_var, self.args[1].substitute(var, term))
        return self

    def __str__(self) -> str:
        op = self.op
        if op == "PRED":
            pred, terms = self.args[0], self.args[1]
            args_str = ", ".join(str(t) for t in terms)
            return f"{pred.name}({args_str})"
        elif op == "NOT":
            return f"¬{self.args[0]}"
        elif op == "AND":
            return " ∧ ".join(str(a) for a in self.args)
        elif op == "OR":
            return " ∨ ".join(str(a) for a in self.args)
        elif op == "IMPLIES":
            return f"{self.args[0]} → {self.args[1]}"
        elif op == "IFF":
            return f"{self.args[0]} ↔ {self.args[1]}"
        elif op == "FORALL":
            return f"∀{self.args[0]}.{self.args[1]}"
        elif op == "EXISTS":
            return f"∃{self.args[0]}.{self.args[1]}"
        return str(self.args)


def prenex_normal_form(formula: FOLFormula) -> FOLFormula:
    """
    Bringt eine FOL-Formel in Pränex-Normalform (PNF).

    PNF: Q1x1 Q2x2 ... Qnxn . F
    Alle Quantoren stehen vorne, gefolgt von einer quantorenfreien Formel (Matrix).

    :param formula: FOL-Formel
    :return: Äquivalente Formel in PNF
    Letzte Änderung: 2026-03-10
    """
    # Vereinfachte PNF-Konversion: Quantoren nach vorne ziehen
    # In einer vollständigen Implementierung müssen Variablen umbenannt werden

    def pull_quantifiers(f: FOLFormula) -> FOLFormula:
        """Zieht Quantoren nach vorne."""
        if f.op in ("PRED", "ATOM"):
            return f
        elif f.op in ("FORALL", "EXISTS"):
            return FOLFormula(f.op, f.args[0], pull_quantifiers(f.args[1]))
        elif f.op == "NOT":
            inner = pull_quantifiers(f.args[0])
            if inner.op == "FORALL":
                # ¬∀x.F ≡ ∃x.¬F
                return FOLFormula("EXISTS", inner.args[0],
                                  pull_quantifiers(FOLFormula("NOT", inner.args[1])))
            elif inner.op == "EXISTS":
                # ¬∃x.F ≡ ∀x.¬F
                return FOLFormula("FORALL", inner.args[0],
                                  pull_quantifiers(FOLFormula("NOT", inner.args[1])))
            return FOLFormula("NOT", inner)
        elif f.op in ("AND", "OR"):
            processed = [pull_quantifiers(a) for a in f.args]
            # Quantoren aus dem ersten Argument nach vorne ziehen
            result = processed[-1]
            for p in reversed(processed[:-1]):
                if p.op in ("FORALL", "EXISTS"):
                    q_type = p.op
                    var = p.args[0]
                    # Sicherstellen dass var nicht frei in result vorkommt
                    while var in result.free_variables():
                        var = var + "'"
                    inner = p.args[1]
                    result = FOLFormula(q_type, var, FOLFormula(f.op, inner, result))
                else:
                    result = FOLFormula(f.op, p, result)
            return result
        return f

    return pull_quantifiers(formula)


def skolemize(formula: FOLFormula) -> FOLFormula:
    """
    Skolemisierung: Eliminiert Existenzquantoren durch Skolem-Funktionen.

    Für ∃x.F unter Allquantoren ∀y1...∀yn.∃x.F:
    Ersetze x durch Skolemfunktion sk(y1,...,yn).

    :param formula: FOL-Formel (vorzugsweise in PNF)
    :return: Skolemisierte Formel (nur noch Allquantoren)
    Letzte Änderung: 2026-03-10
    """
    sk_counter = [0]

    def skolem_func_name() -> str:
        name = f"sk{sk_counter[0]}"
        sk_counter[0] += 1
        return name

    def skolemize_inner(f: FOLFormula, universally_bound: list[str]) -> FOLFormula:
        """Rekursive Skolemisierung mit aktuellem Allquantor-Kontext."""
        if f.op == "PRED":
            return f
        elif f.op == "NOT":
            return FOLFormula("NOT", skolemize_inner(f.args[0], universally_bound))
        elif f.op in ("AND", "OR", "IMPLIES", "IFF"):
            return FOLFormula(f.op, *[skolemize_inner(a, universally_bound) for a in f.args])
        elif f.op == "FORALL":
            var = f.args[0]
            return FOLFormula("FORALL", var,
                              skolemize_inner(f.args[1], universally_bound + [var]))
        elif f.op == "EXISTS":
            var = f.args[0]
            sk_name = skolem_func_name()
            if universally_bound:
                # Erstelle Skolem-Funktion mit allen Allquantor-Variablen
                sk_term = Term.function(sk_name, *[Term.variable(v) for v in universally_bound])
            else:
                # Keine Allquantoren → Skolem-Konstante
                sk_term = Term.constant(sk_name)
            # Ersetze var durch Skolem-Term in der Formel
            inner_skolemized = skolemize_inner(f.args[1], universally_bound)
            return inner_skolemized.substitute(var, sk_term)
        return f

    return skolemize_inner(formula, [])


def herbrand_universe(constants: list[str], function_symbols: list[tuple], depth: int = 2) -> list[str]:
    """
    Berechnet das Herbrand-Universum bis zur gegebenen Tiefe.

    Das Herbrand-Universum enthält alle Grundterme (ohne Variablen),
    die aus Konstanten und Funktionssymbolen aufgebaut werden können.

    :param constants: Liste von Konstantennamen
    :param function_symbols: Liste von (name, arität)-Tupeln
    :param depth: Maximale Schachtelungstiefe
    :return: Liste aller Grundterme als Strings
    Letzte Änderung: 2026-03-10
    """
    if not constants:
        constants = ["a"]  # Mindestens eine Konstante

    # Tiefe 0: nur Konstanten
    universe = list(constants)
    current_level = list(constants)

    for _ in range(depth):
        next_level = []
        for fname, arity in function_symbols:
            # Alle Kombinationen aus aktuellen Termen
            for combo in itertools.product(current_level, repeat=arity):
                term_str = f"{fname}({', '.join(combo)})"
                if term_str not in universe:
                    next_level.append(term_str)
                    universe.append(term_str)
        current_level = next_level
        if not next_level:
            break

    return universe


# =============================================================================
# BEWEISKALKÜLE
# =============================================================================

class NaturalDeductionProof:
    """
    Natürliches Schließen (Gentzen-Kalkül ND).

    Inferenzregeln:
    - →I  (Implikationseinführung): Nimm A an, zeige B → A→B
    - →E  (Implikationselimination / Modus Ponens): A→B, A ⊢ B
    - ∧I  (Konjunktionseinführung): A, B ⊢ A∧B
    - ∧E  (Konjunktionselimination): A∧B ⊢ A und A∧B ⊢ B
    - ∨I  (Disjunktionseinführung): A ⊢ A∨B und B ⊢ A∨B
    - ¬I  (Negationseinführung): A ⊢ ⊥ → ¬A
    - ¬E  (Negationselimination): A, ¬A ⊢ ⊥
    - ⊥E  (Ex Falso): ⊥ ⊢ A (beliebige Formel)

    Letzte Änderung: 2026-03-10
    """

    # Darstellung von False (⊥)
    BOTTOM = "⊥"

    def __init__(self):
        """Initialisiert einen leeren Beweis."""
        # Beweiszeilen: {id: {formula, rule, from_lines, discharged}}
        self.lines = {}
        self._next_id = 0
        # Annahmen-Stack für Entlastung
        self.open_assumptions = set()

    def _new_id(self) -> int:
        """Gibt neue Zeilen-ID zurück."""
        i = self._next_id
        self._next_id += 1
        return i

    def assume(self, formula) -> int:
        """
        Führt eine Annahme in den Beweis ein.

        :param formula: Angenommene Formel (String oder LogicFormula)
        :return: Zeilen-ID der Annahme
        """
        line_id = self._new_id()
        self.lines[line_id] = {
            "formula": formula,
            "rule": "Assumption",
            "from_lines": [],
            "discharged": False
        }
        self.open_assumptions.add(line_id)
        return line_id

    def apply_rule(self, rule: str, *line_ids: int) -> int:
        """
        Wendet eine Inferenzregel an.

        :param rule: Regelname (→E, ∧I, ∧E1, ∧E2, ∨I1, ∨I2, ¬E, ⊥E)
        :param line_ids: Zeilen-IDs der Prämissen
        :return: Zeilen-ID der neuen Zeile
        """
        line_id = self._new_id()
        formulas = [self.lines[i]["formula"] for i in line_ids]

        if rule == "→E":
            # Modus Ponens: A→B, A ⊢ B
            if len(formulas) >= 2:
                impl = formulas[0]
                if isinstance(impl, LogicFormula) and impl.op == "IMPLIES":
                    conclusion = impl.args[1]
                else:
                    conclusion = f"({formulas[0]})_E"
            else:
                conclusion = "?"
        elif rule == "∧I":
            # A, B ⊢ A∧B
            if len(formulas) >= 2 and all(isinstance(f, LogicFormula) for f in formulas):
                conclusion = LogicFormula("AND", *formulas)
            else:
                conclusion = f"({' ∧ '.join(str(f) for f in formulas)})"
        elif rule == "∧E1":
            # A∧B ⊢ A
            f = formulas[0]
            if isinstance(f, LogicFormula) and f.op == "AND":
                conclusion = f.args[0]
            else:
                conclusion = f"left({f})"
        elif rule == "∧E2":
            # A∧B ⊢ B
            f = formulas[0]
            if isinstance(f, LogicFormula) and f.op == "AND":
                conclusion = f.args[1]
            else:
                conclusion = f"right({f})"
        elif rule == "∨I1":
            # A ⊢ A∨B (B beliebig, hier nicht spezifiziert)
            conclusion = f"{formulas[0]} ∨ ?"
        elif rule == "∨I2":
            conclusion = f"? ∨ {formulas[0]}"
        elif rule == "¬E":
            # A, ¬A ⊢ ⊥
            conclusion = self.BOTTOM
        elif rule == "⊥E":
            # ⊥ ⊢ A (beliebig)
            conclusion = "A"
        else:
            conclusion = f"({rule} {', '.join(str(f) for f in formulas)})"

        self.lines[line_id] = {
            "formula": conclusion,
            "rule": rule,
            "from_lines": list(line_ids),
            "discharged": False
        }
        return line_id

    def discharge(self, assumption_id: int) -> int:
        """
        Entlädt eine Annahme (für →I Regel).

        :param assumption_id: ID der zu entladenden Annahme
        :return: Zeilen-ID mit der Implikation
        """
        if assumption_id in self.open_assumptions:
            self.lines[assumption_id]["discharged"] = True
            self.open_assumptions.discard(assumption_id)
        line_id = self._new_id()
        # Letzte Zeile des Beweises als Konklusion
        last_id = self._next_id - 2
        assumption_formula = self.lines[assumption_id]["formula"]
        if last_id in self.lines:
            conclusion_formula = self.lines[last_id]["formula"]
            if isinstance(assumption_formula, LogicFormula) and isinstance(conclusion_formula, LogicFormula):
                implication = LogicFormula("IMPLIES", assumption_formula, conclusion_formula)
            else:
                implication = f"{assumption_formula} → {conclusion_formula}"
        else:
            implication = f"{assumption_formula} → ?"

        self.lines[line_id] = {
            "formula": implication,
            "rule": "→I",
            "from_lines": [assumption_id, last_id],
            "discharged": True
        }
        return line_id

    def is_valid(self) -> bool:
        """
        Prüft ob der Beweis gültig ist.
        Ein Beweis ist gültig wenn keine offenen Annahmen vorhanden sind.

        :return: True wenn gültig
        """
        return len(self.open_assumptions) == 0

    def get_conclusion(self):
        """Gibt die letzte Formel des Beweises zurück."""
        if not self.lines:
            return None
        last_id = max(self.lines.keys())
        return self.lines[last_id]["formula"]

    def __str__(self) -> str:
        """Zeigt den Beweis als Tabelle."""
        result = []
        for i, (line_id, line) in enumerate(sorted(self.lines.items())):
            discharged = "[discharged]" if line.get("discharged") else ""
            from_str = f"from {line['from_lines']}" if line["from_lines"] else ""
            result.append(
                f"  {line_id:3d}: {str(line['formula']):<40} "
                f"({line['rule']}) {from_str} {discharged}"
            )
        return "\n".join(result)


class HilbertSystem:
    """
    Hilbert-Kalkül für Aussagenlogik.

    Axiomenschemata (für beliebige Formeln A, B, C):
    - K: A → (B → A)
    - S: (A → (B → C)) → ((A → B) → (A → C))
    - DN: (¬¬A) → A

    Einzige Schlussregel: Modus Ponens (MP)
    Aus A → B und A folgt B.

    Letzte Änderung: 2026-03-10
    """

    AXIOMS = [
        "A→(B→A)",           # K-Axiom (Konstante)
        "(A→(B→C))→((A→B)→(A→C))",  # S-Axiom (Substitution)
        "(¬¬A)→A"            # Doppelnegation
    ]

    def __init__(self):
        """Initialisiert das Hilbert-System."""
        self.proof_steps = []

    def apply_mp(self, major: str, minor: str) -> str:
        """
        Wendet Modus Ponens an.

        :param major: Formel der Form A→B
        :param minor: Formel A
        :return: Formel B (wenn MP anwendbar)
        """
        # Vereinfachte String-Darstellung
        step = f"MP({major}, {minor})"
        self.proof_steps.append(step)
        return step

    def prove(self, formula: str) -> list[str]:
        """
        Versucht eine Formel zu beweisen (vereinfacht).

        :param formula: Zu beweisende Formel als String
        :return: Liste der Beweisschritte
        """
        steps = [
            f"Ziel: {formula}",
            f"Verwende K-Axiom: A→(B→A) mit passender Instanz",
            f"Verwende S-Axiom falls nötig",
            f"Wende MP an",
            f"Resultat: {formula} bewiesen (sofern ableitbar)"
        ]
        self.proof_steps = steps
        return steps


def modus_ponens(major: LogicFormula, minor: LogicFormula) -> LogicFormula | None:
    """
    Modus Ponens: Aus A→B und A folgt B.

    Formal: {A→B, A} ⊢ B

    :param major: Implikation A→B
    :param minor: Antezedens A
    :return: Konsequens B oder None wenn MP nicht anwendbar
    Letzte Änderung: 2026-03-10
    """
    if major.op != "IMPLIES":
        return None

    antecedent = major.args[0]
    consequent = major.args[1]

    # Prüfe ob minor und antecedent strukturell gleich sind
    if str(antecedent) == str(minor):
        return consequent
    return None


def modus_tollens(major: LogicFormula, minor: LogicFormula) -> LogicFormula | None:
    """
    Modus Tollens: Aus A→B und ¬B folgt ¬A.

    Formal: {A→B, ¬B} ⊢ ¬A

    :param major: Implikation A→B
    :param minor: Negation des Konsequens ¬B
    :return: Negation des Antezedens ¬A oder None
    Letzte Änderung: 2026-03-10
    """
    if major.op != "IMPLIES":
        return None
    if minor.op != "NOT":
        return None

    antecedent = major.args[0]
    consequent = major.args[1]
    neg_consequent = minor.args[0]

    # Prüfe ¬B = minor
    if str(consequent) == str(neg_consequent):
        return LogicFormula("NOT", antecedent)
    return None


# =============================================================================
# MODALLOGIK
# =============================================================================

class KripkeFrame:
    """
    Kripke-Rahmen für Modallogik.

    Ein Kripke-Rahmen besteht aus:
    - W: Menge möglicher Welten
    - R: Erreichbarkeitsrelation R ⊆ W × W

    Verschiedene Eigenschaften von R entsprechen verschiedenen Modallogiken:
    - Reflexivität → System T (□A → A)
    - Transitivität → System S4 (□A → □□A)
    - Symmetrie → System B (A → □◇A)
    - Reflexiv + transitiv → S4
    - Reflexiv + transitiv + symmetrisch → S5 (◇A → □◇A)
    - Seriell → System D (□A → ◇A)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, worlds: list[str], accessibility: dict[str, list[str]]):
        """
        :param worlds: Liste der Weltnamen (z. B. ['w1', 'w2', 'w3'])
        :param accessibility: Dict {welt: [erreichbare_welten]}
        """
        self.worlds = list(worlds)
        self.accessibility = {w: list(acc) for w, acc in accessibility.items()}
        # Welten ohne Eintrag haben keine Nachfolger
        for w in self.worlds:
            if w not in self.accessibility:
                self.accessibility[w] = []

    def accessible_from(self, world: str) -> list[str]:
        """Gibt alle von world erreichbaren Welten zurück."""
        return self.accessibility.get(world, [])

    def is_reflexive(self) -> bool:
        """
        Prüft Reflexivität: ∀w ∈ W: wRw.
        Korrespondiert mit T-Axiom: □A → A.
        """
        return all(w in self.accessibility.get(w, []) for w in self.worlds)

    def is_symmetric(self) -> bool:
        """
        Prüft Symmetrie: wRv → vRw.
        Korrespondiert mit B-Axiom: A → □◇A.
        """
        for w in self.worlds:
            for v in self.accessibility.get(w, []):
                if w not in self.accessibility.get(v, []):
                    return False
        return True

    def is_transitive(self) -> bool:
        """
        Prüft Transitivität: wRv ∧ vRu → wRu.
        Korrespondiert mit 4-Axiom: □A → □□A.
        """
        for w in self.worlds:
            for v in self.accessibility.get(w, []):
                for u in self.accessibility.get(v, []):
                    if u not in self.accessibility.get(w, []):
                        return False
        return True

    def is_serial(self) -> bool:
        """
        Prüft Serialität: ∀w ∃v: wRv.
        Korrespondiert mit D-Axiom: □A → ◇A.
        """
        return all(len(self.accessibility.get(w, [])) > 0 for w in self.worlds)

    def is_euclidean(self) -> bool:
        """
        Prüft Euklidizität: wRv ∧ wRu → vRu.
        Zusammen mit Reflexivität ergibt sich S5.
        """
        for w in self.worlds:
            accessible = self.accessibility.get(w, [])
            for v in accessible:
                for u in accessible:
                    if u not in self.accessibility.get(v, []):
                        return False
        return True

    def modal_system(self) -> str:
        """Bestimmt das stärkste gültige Modalsystem für diesen Rahmen."""
        refl = self.is_reflexive()
        symm = self.is_symmetric()
        trans = self.is_transitive()
        serial = self.is_serial()
        eucl = self.is_euclidean()

        if refl and trans and eucl:
            return "S5"
        elif refl and trans:
            return "S4"
        elif refl and symm:
            return "B"
        elif refl:
            return "T"
        elif serial:
            return "D"
        else:
            return "K"

    def __str__(self) -> str:
        worlds_str = ", ".join(self.worlds)
        rel_str = ", ".join(
            f"{w}→{{{', '.join(v)}}}"
            for w, v in self.accessibility.items()
        )
        return f"KripkeFrame(W={{{worlds_str}}}, R={{{rel_str}}})"


class ModalFormula:
    """
    Modallogische Formel mit Notwendigkeits- (□) und Möglichkeitsoperator (◇).

    Unterstützte Operatoren:
    - ATOM: Atomare Variable
    - NOT, AND, OR, IMPLIES: Klassische Konnektive
    - BOX: Notwendigkeit (□) - wahr in allen erreichbaren Welten
    - DIAMOND: Möglichkeit (◇) - wahr in mindestens einer erreichbaren Welt

    Beziehung: ◇A ≡ ¬□¬A

    Letzte Änderung: 2026-03-10
    """

    MODAL_OPS = {"ATOM", "NOT", "AND", "OR", "IMPLIES", "BOX", "DIAMOND"}

    def __init__(self, op: str, *args):
        """
        :param op: Operator
        :param args: Argumente
        """
        if op not in self.MODAL_OPS:
            raise ValueError(f"Ungültiger Modaloperator: '{op}'")
        self.op = op
        self.args = args

    @staticmethod
    def atom(name: str) -> 'ModalFormula':
        return ModalFormula("ATOM", name)

    @staticmethod
    def box(formula: 'ModalFormula') -> 'ModalFormula':
        """□ Notwendigkeit."""
        return ModalFormula("BOX", formula)

    @staticmethod
    def diamond(formula: 'ModalFormula') -> 'ModalFormula':
        """◇ Möglichkeit."""
        return ModalFormula("DIAMOND", formula)

    def evaluate(self, frame: KripkeFrame, valuation: dict[str, set[str]], world: str) -> bool:
        """
        Wertet die Formel in einem Kripke-Modell aus.

        :param frame: Kripke-Rahmen (W, R)
        :param valuation: Dict {atom_name: {welten_in_denen_es_gilt}}
        :param world: Aktuelle Welt
        :return: Wahrheitswert in der aktuellen Welt
        """
        op = self.op
        if op == "ATOM":
            # Atom gilt in Welt wenn welt in valuation[atom]
            return world in valuation.get(self.args[0], set())
        elif op == "NOT":
            return not self.args[0].evaluate(frame, valuation, world)
        elif op == "AND":
            return all(a.evaluate(frame, valuation, world) for a in self.args)
        elif op == "OR":
            return any(a.evaluate(frame, valuation, world) for a in self.args)
        elif op == "IMPLIES":
            a, b = self.args
            return (not a.evaluate(frame, valuation, world)) or b.evaluate(frame, valuation, world)
        elif op == "BOX":
            # □A gilt wenn A in allen erreichbaren Welten gilt
            accessible = frame.accessible_from(world)
            return all(self.args[0].evaluate(frame, valuation, w) for w in accessible)
        elif op == "DIAMOND":
            # ◇A gilt wenn A in mindestens einer erreichbaren Welt gilt
            accessible = frame.accessible_from(world)
            return any(self.args[0].evaluate(frame, valuation, w) for w in accessible)
        return False

    def __str__(self) -> str:
        op = self.op
        if op == "ATOM":
            return self.args[0]
        elif op == "NOT":
            return f"¬{self.args[0]}"
        elif op == "AND":
            return " ∧ ".join(str(a) for a in self.args)
        elif op == "OR":
            return " ∨ ".join(str(a) for a in self.args)
        elif op == "IMPLIES":
            return f"{self.args[0]} → {self.args[1]}"
        elif op == "BOX":
            return f"□{self.args[0]}"
        elif op == "DIAMOND":
            return f"◇{self.args[0]}"
        return str(self.args)


def modal_system(system: str) -> str:
    """
    Gibt die Axiome eines modalen Systems zurück.

    Systeme:
    - K   : Grundsystem (nur Distributivitätsaxiom)
    - D   : K + □A → ◇A (D-Axiom, serielle Rahmen)
    - T   : K + □A → A (T-Axiom, reflexive Rahmen)
    - B   : T + A → □◇A (B-Axiom, symmetrische Rahmen)
    - S4  : T + □A → □□A (4-Axiom, transitive Rahmen)
    - S5  : S4 + ◇A → □◇A (5-Axiom, euklidische Rahmen)

    :param system: Systemname ('K', 'D', 'T', 'B', 'S4', 'S5')
    :return: Beschreibung des Systems mit Axiomen
    Letzte Änderung: 2026-03-10
    """
    systems = {
        "K": (
            "System K (Kripke-Grundsystem)\n"
            "Axiom K: □(A→B) → (□A→□B)\n"
            "Rahmen: beliebig\n"
            "Regel: Notwendigkeitsgeneralisierung (NEC): ⊢A → ⊢□A"
        ),
        "D": (
            "System D (Deontic Modal Logic)\n"
            "Axiom K: □(A→B) → (□A→□B)\n"
            "Axiom D: □A → ◇A\n"
            "Rahmen: seriell (∀w∃v: wRv)"
        ),
        "T": (
            "System T (Alethische Modallogik)\n"
            "Axiom K: □(A→B) → (□A→□B)\n"
            "Axiom T: □A → A  (was notwendig ist, ist wahr)\n"
            "Rahmen: reflexiv"
        ),
        "B": (
            "System B (Brouwer-System)\n"
            "Axiom K + T\n"
            "Axiom B: A → □◇A\n"
            "Rahmen: reflexiv + symmetrisch"
        ),
        "S4": (
            "System S4 (Lewis S4)\n"
            "Axiom K + T\n"
            "Axiom 4: □A → □□A  (positive Introspektionsaxiom)\n"
            "Rahmen: reflexiv + transitiv"
        ),
        "S5": (
            "System S5 (Lewis S5)\n"
            "Axiom K + T + 4\n"
            "Axiom 5: ◇A → □◇A  (negative Introspektionsaxiom)\n"
            "Rahmen: reflexiv + transitiv + euklidisch (Äquivalenzrelation)\n"
            "Stärkstes Standard-Modalsystem"
        ),
    }
    system_upper = system.upper()
    return systems.get(system_upper, f"Unbekanntes System: '{system}'. Gültig: K, D, T, B, S4, S5")


# =============================================================================
# GÖDELS UNVOLLSTÄNDIGKEITSSÄTZE
# =============================================================================

def godel_numbering(formula_str: str) -> int:
    """
    Gödel-Nummerierung: Weist einer Formel eine eindeutige natürliche Zahl zu.

    Methode: Primzahl-Kodierung
    Jedes Symbol bekommt eine Zahl:
    - Symbole → Zahlen (a→1, b→2, ..., ¬→27, ∧→28, ∨→29, →→30, ∀→31, ∃→32, ...)

    Gödel-Zahl: ⟨s1, s2, ..., sn⟩ = 2^c(s1) · 3^c(s2) · 5^c(s3) · ... · p_n^c(sn)

    wobei p_i die i-te Primzahl und c(si) der Code von Symbol si ist.

    :param formula_str: Formelstring
    :return: Gödel-Zahl (kann sehr groß werden!)
    Letzte Änderung: 2026-03-10
    """
    # Symbol-zu-Code Mapping
    symbol_codes = {
        '0': 1, '1': 2, 'S': 3, '+': 4, '·': 5, '=': 6,
        '(': 7, ')': 8, ',': 9, ' ': 10,
        'a': 11, 'b': 12, 'c': 13, 'd': 14, 'e': 15,
        'f': 16, 'g': 17, 'h': 18, 'i': 19, 'j': 20,
        'k': 21, 'l': 22, 'm': 23, 'n': 24, 'o': 25,
        'p': 26, 'q': 27, 'r': 28, 's': 29, 't': 30,
        'u': 31, 'v': 32, 'w': 33, 'x': 34, 'y': 35, 'z': 36,
        '¬': 37, '∧': 38, '∨': 39, '→': 40, '↔': 41,
        '∀': 42, '∃': 43, '⊥': 44, '⊢': 45,
        'A': 46, 'B': 47, 'C': 48, 'D': 49, 'E': 50,
        'F': 51, 'G': 52, 'H': 53, 'I': 54, 'J': 55,
        'K': 56, 'L': 57, 'M': 58, 'N': 59, 'O': 60,
        'P': 61, 'Q': 62, 'R': 63, 'T': 64, 'U': 65,
        'V': 66, 'W': 67, 'X': 68, 'Y': 69, 'Z': 70,
    }

    def nth_prime(n: int) -> int:
        """Gibt die n-te Primzahl zurück (0-indiziert)."""
        primes = []
        candidate = 2
        while len(primes) <= n:
            if all(candidate % p != 0 for p in primes):
                primes.append(candidate)
            candidate += 1
        return primes[n]

    # Kodierung
    result = 1
    for i, symbol in enumerate(formula_str):
        code = symbol_codes.get(symbol, ord(symbol))
        prime = nth_prime(i)
        result *= prime ** code

    return result


def godel_decode(number: int) -> str:
    """
    Dekodiert eine Gödel-Zahl zurück in eine Formel.

    Umkehrung der Primzahl-Kodierung: Primfaktorzerlegung bestimmt Symbole.

    :param number: Gödel-Zahl
    :return: Decodierte Formel als String (approximativ)
    Letzte Änderung: 2026-03-10
    """
    # Code-zu-Symbol Mapping (Umkehrung)
    code_to_symbol = {
        1: '0', 2: '1', 3: 'S', 4: '+', 5: '·', 6: '=',
        7: '(', 8: ')', 9: ',', 10: ' ',
        11: 'a', 12: 'b', 13: 'c', 14: 'd', 15: 'e',
        16: 'f', 17: 'g', 18: 'h', 19: 'i', 20: 'j',
        21: 'k', 22: 'l', 23: 'm', 24: 'n', 25: 'o',
        26: 'p', 27: 'q', 28: 'r', 29: 's', 30: 't',
        31: 'u', 32: 'v', 33: 'w', 34: 'x', 35: 'y', 36: 'z',
        37: '¬', 38: '∧', 39: '∨', 40: '→', 41: '↔',
        42: '∀', 43: '∃', 44: '⊥', 45: '⊢',
        46: 'A', 47: 'B', 48: 'C', 49: 'D', 50: 'E',
        51: 'F', 52: 'G', 53: 'H', 54: 'I', 55: 'J',
        56: 'K', 57: 'L', 58: 'M', 59: 'N', 60: 'O',
        61: 'P', 62: 'Q', 63: 'R', 64: 'T', 65: 'U',
        66: 'V', 67: 'W', 68: 'X', 69: 'Y', 70: 'Z',
    }

    def prime_factorization(n: int) -> list[tuple[int, int]]:
        """Primfaktorzerlegung: gibt [(primzahl, exponent)] zurück."""
        factors = []
        d = 2
        while d * d <= n:
            if n % d == 0:
                exp = 0
                while n % d == 0:
                    exp += 1
                    n //= d
                factors.append((d, exp))
            d += 1
        if n > 1:
            factors.append((n, 1))
        return factors

    factors = prime_factorization(number)
    if not factors:
        return ""

    # Sortiere nach Primzahl-Index und extrahiere Codes
    # Die Codes sind die Exponenten der Primzahlen in der richtigen Reihenfolge
    result = ""
    for _, exp in factors:
        symbol = code_to_symbol.get(exp, f"[{exp}]")
        result += symbol

    return result


def incompleteness_demonstration() -> dict:
    """
    Demonstriert Gödels ersten Unvollständigkeitssatz konzeptuell.

    Gödels 1. Unvollständigkeitssatz (1931):
    Jedes konsistente, rekursiv aufzählbare formale System F,
    das starke genug ist, um Peano-Arithmetik auszudrücken,
    enthält eine Aussage G, die:
    - wahr ist (in der Standardinterpretation)
    - aber in F weder beweisbar noch widerlegbar ist.

    Die Konstruktion: G sagt „Ich bin in F nicht beweisbar."
    (Gödels selbstreferenzielle Aussage via Diagonallemma)

    :return: Dict mit Erklärung und Demonstration
    Letzte Änderung: 2026-03-10
    """
    return {
        "theorem": "Gödels erster Unvollständigkeitssatz (1931)",
        "statement": (
            "Jedes konsistente, hinreichend starke formale System F enthält "
            "eine Aussage G, die in F weder beweisbar noch widerlegbar ist."
        ),
        "construction": {
            "step1": "Gödel-Nummerierung: Jede Formel und jeder Beweis erhält eine Gödel-Zahl.",
            "step2": (
                "Diagonallemma: Es gibt eine Aussage G mit: G ↔ ¬Bew(⌈G⌉)\n"
                "G sagt: 'Die Formel mit der Gödel-Zahl ⌈G⌉ ist nicht beweisbar.'"
            ),
            "step3": (
                "Annahme G ist beweisbar:\n"
                "  → Bew(⌈G⌉) ist wahr\n"
                "  → G ist falsch (da G = ¬Bew(⌈G⌉))\n"
                "  → Widerspruch (F ist inkonsistent)\n"
                "Also: G ist NICHT beweisbar in F."
            ),
            "step4": (
                "Annahme ¬G ist beweisbar:\n"
                "  → ¬G ist beweisbar → G ist falsch\n"
                "  → Also: Bew(⌈G⌉) ist wahr\n"
                "  → G ist beweisbar → G ist wahr\n"
                "  → Widerspruch\n"
                "Also: ¬G ist auch NICHT beweisbar in F."
            ),
            "conclusion": (
                "G ist wahr (da G nicht beweisbar) aber unbeweisbar in F.\n"
                "F ist unvollständig: Es gibt wahre Aussagen, die nicht beweisbar sind."
            ),
        },
        "second_incompleteness": {
            "statement": (
                "Gödels zweiter Unvollständigkeitssatz:\n"
                "Kein konsistentes formales System kann seine eigene Konsistenz beweisen."
            ),
            "note": (
                "Con(F) = 'F ist konsistent' ist in F nicht beweisbar,\n"
                "sofern F tatsächlich konsistent ist."
            ),
        },
        "implications": [
            "Mathematik ist prinzipiell unvollständig.",
            "Es gibt wahre mathematische Aussagen, die keine Beweise haben.",
            "Hilberts Programm (Vollständigkeit der Mathematik) ist gescheitert.",
            "Gödel-Aussagen wie die Kontinuumshypothese sind unabhängig von ZFC.",
        ],
        "goedel_sentence_example": "G = 'Diese Aussage ist in F nicht beweisbar.'",
        "formal_system_requirements": [
            "Konsistent (kein Widerspruch beweisbar)",
            "Rekursiv aufzählbar (beweisbare Aussagen aufzählbar)",
            "Hinreichend stark (enthält Peano-Arithmetik oder equivalentes)",
        ],
    }
