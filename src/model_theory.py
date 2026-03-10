"""
@file model_theory.py
@brief Modelltheorie: Strukturen, Erfüllbarkeit, Kompaktheitssatz, Löwenheim-Skolem,
       Typentheorie, Kategorien und Quantorenelimination.

@description
    Dieses Modul implementiert die zentralen Konzepte der mathematischen Modelltheorie:

    1. Signature        – Formale Signatur L = (F, R, C) mit Arität
    2. Structure        – L-Struktur M = (A, F^M, R^M, C^M) mit Tarski-Wahrheitsbedingungen
    3. ElementaryEmbedding – Elementare Einbettung f: M → N
    4. Kompaktheitssatz – Demonstration mit Nicht-Standard-Arithmetik
    5. Löwenheim-Skolem – Auf- und Abwärts-Sätze
    6. Type             – p-Typen, Vollständigkeit, Realisierung, Auslassung
    7. Kategorien       – κ-Kategorik, Vaught-Kriterium, Morley-Satz
    8. Quantorenelimination – DLO, Th(ℚ,<), ACF

    Formelsprache (vereinfachtes Parsen):
      Atomare Formeln:   "R(x,y)", "x=y", "=(x,y)"
      Negation:          "NOT(phi)"
      Konjunktion:       "AND(phi,psi)"
      Disjunktion:       "OR(phi,psi)"
      Implikation:       "IMPLIES(phi,psi)"
      Äquivalenz:        "IFF(phi,psi)"
      Existenzquantor:   "EXISTS(x,phi)"
      Allquantor:        "FORALL(x,phi)"

@author Kurt Ingwer
@lastModified 2026-03-10
"""

from __future__ import annotations
import itertools
import re
from typing import Any, Callable, Optional


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Signatur
# ═══════════════════════════════════════════════════════════════════════════════

class Signature:
    """
    Formale Signatur L = (F, R, C).

    Eine Signatur legt die nicht-logischen Symbole einer Sprache fest:
    - Funktionssymbole F mit Arität (f: A^n → A)
    - Relationssymbole  R mit Arität (R ⊆ A^n)
    - Konstantensymbole C (= nullstellige Funktionen)

    Beispiel Gruppentheorie: F={'*': 2, 'inv': 1}, R={}, C={'e': 0}
    Beispiel Ordnungstheorie: F={}, R={'<': 2, '≤': 2}, C={}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        functions: dict[str, int],   # Funktionssymbol → Arität
        relations: dict[str, int],   # Relationssymbol → Arität
        constants: list[str]         # Konstantensymbole (Arität 0)
    ):
        """
        Initialisiert die Signatur.

        @param functions  Dict {Symbolname: Arität} für Funktionen
        @param relations  Dict {Symbolname: Arität} für Relationen
        @param constants  Liste der Konstantensymbole
        """
        # Funktionssymbole mit ihren Aritäten
        self.functions: dict[str, int] = dict(functions)
        # Relationssymbole mit ihren Aritäten
        self.relations: dict[str, int] = dict(relations)
        # Konstantensymbole als Menge
        self.constants: list[str] = list(constants)

    def all_symbols(self) -> dict[str, str]:
        """
        Gibt alle Symbole der Signatur mit Typ zurück.

        @return Dict {Symbol: Typ} wobei Typ ∈ {'function','relation','constant'}
        """
        result: dict[str, str] = {}
        for sym in self.functions:
            result[sym] = 'function'
        for sym in self.relations:
            result[sym] = 'relation'
        for sym in self.constants:
            result[sym] = 'constant'
        return result

    def __repr__(self) -> str:
        return (f"Signature(F={self.functions}, R={self.relations}, "
                f"C={self.constants})")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. L-Struktur
# ═══════════════════════════════════════════════════════════════════════════════

class Structure:
    """
    L-Struktur M = (A, F^M, R^M, C^M).

    Eine L-Struktur interpretiert alle Symbole der Signatur L über
    einer nichtleeren Grundmenge (Universum) A:
    - f^M: A^n → A für jedes n-stellige Funktionssymbol f
    - R^M ⊆ A^n für jedes n-stellige Relationssymbol R
    - c^M ∈ A für jede Konstante c

    Tarski-Wahrheitsbedingungen (rekursiv):
      M ⊨ R(t₁,…,tₙ)[s]  ⟺  (t₁^M[s],…,tₙ^M[s]) ∈ R^M
      M ⊨ ¬φ[s]          ⟺  M ⊭ φ[s]
      M ⊨ φ∧ψ[s]         ⟺  M⊨φ[s]  und  M⊨ψ[s]
      M ⊨ φ∨ψ[s]         ⟺  M⊨φ[s]  oder M⊨ψ[s]
      M ⊨ ∃x φ[s]        ⟺  ∃a∈A: M⊨φ[s(x↦a)]
      M ⊨ ∀x φ[s]        ⟺  ∀a∈A: M⊨φ[s(x↦a)]

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        universe: list,                       # Grundmenge A ≠ ∅
        signature: Signature,                  # Formale Signatur
        functions: dict[str, Callable],        # f^M: aufrufbare Python-Funktionen
        relations: dict[str, set],             # R^M: Menge von Tupeln
        constants: dict[str, Any]              # c^M: Elemente aus A
    ):
        """
        Initialisiert die L-Struktur.

        @param universe   Nichtleere Liste der Grundmengenelemente
        @param signature  Signatur L
        @param functions  Dict {Symbolname: callable} für Funktionsinterpretation
        @param relations  Dict {Symbolname: set of tuples} für Relationsinterpretation
        @param constants  Dict {Symbolname: element} für Konstanteninterpretation
        @raises ValueError wenn das Universum leer ist
        """
        if not universe:
            raise ValueError("Das Universum einer L-Struktur muss nichtleer sein.")
        # Universum als Liste gespeichert (Reihenfolge für ∃/∀-Auswertung)
        self.universe: list = list(universe)
        # Signatur der Struktur
        self.signature: Signature = signature
        # Interpretation der Funktionssymbole
        self.functions: dict[str, Callable] = dict(functions)
        # Interpretation der Relationssymbole (Menge von Tupeln)
        self.relations: dict[str, set] = {k: set(v) for k, v in relations.items()}
        # Interpretation der Konstantensymbole
        self.constants: dict[str, Any] = dict(constants)

    # ──────────────────────────────────────────────────────────────────────────
    # Termauswertung
    # ──────────────────────────────────────────────────────────────────────────

    def interpret_term(self, term: str, assignment: dict) -> Any:
        """
        Wertet einen Term t unter Variablenbelegung s aus: t^M[s].

        Unterstützte Termformen:
        - Variablen:     "x", "y" (in assignment)
        - Konstanten:    "c", "e" (in self.constants)
        - Zahlliterale:  "0", "1", "42" (werden als int interpretiert)
        - Funktionsaufrufe: "f(t₁,…,tₙ)" (rekursiv)

        @param term        Term als String
        @param assignment  Variablenbelegung s: {Var → Element}
        @return            Auswertung t^M[s] ∈ A
        @raises KeyError   wenn Variable oder Konstante unbekannt
        """
        term = term.strip()

        # Zahlliteral (int oder float)
        try:
            return int(term)
        except ValueError:
            pass
        try:
            return float(term)
        except ValueError:
            pass

        # Variable in der aktuellen Belegung
        if term in assignment:
            return assignment[term]

        # Konstantensymbol der Struktur
        if term in self.constants:
            return self.constants[term]

        # Funktionsanwendung: "f(arg1, arg2, ...)"
        match = re.fullmatch(r'(\w+)\((.+)\)', term)
        if match:
            fname = match.group(1)
            args_str = match.group(2)
            # Argumente korrekt aufteilen (berücksichtigt verschachtelte Klammern)
            args = _split_args(args_str)
            evaluated_args = [self.interpret_term(a, assignment) for a in args]
            if fname in self.functions:
                return self.functions[fname](*evaluated_args)
            raise KeyError(f"Unbekanntes Funktionssymbol: '{fname}'")

        raise KeyError(f"Term '{term}' nicht auswertbar: weder Variable noch Konstante.")

    # ──────────────────────────────────────────────────────────────────────────
    # Formelauswertung
    # ──────────────────────────────────────────────────────────────────────────

    def satisfies_atomic(self, formula: str, assignment: dict) -> bool:
        """
        Prüft M ⊨ φ[s] für eine atomare Formel φ.

        Atomare Formeln:
        - Gleichheit:   "=(t1,t2)"  oder  "t1=t2"
        - Relation:     "R(t1,...,tn)"

        @param formula    Atomare Formel als String
        @param assignment Variablenbelegung
        @return           True gdw. M ⊨ φ[s]
        """
        formula = formula.strip()

        # Gleichheitsatom in Präfixnotation: "=(t1,t2)"
        if formula.startswith('=('):
            inner = formula[2:-1]
            parts = _split_args(inner)
            if len(parts) == 2:
                lhs = self.interpret_term(parts[0], assignment)
                rhs = self.interpret_term(parts[1], assignment)
                return lhs == rhs

        # Gleichheitsatom in Infixnotation: "t1=t2"
        # (nur wenn kein Klammerpaar außen vorhanden)
        if '=' in formula and not formula.startswith('('):
            # Prüfen ob es wirklich eine Gleichung und kein Funktionsaufruf ist
            eq_pos = formula.find('=')
            # Sicherstellen, dass das = nicht innerhalb von Klammern liegt
            depth = 0
            for i, ch in enumerate(formula):
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                elif ch == '=' and depth == 0:
                    lhs = formula[:i].strip()
                    rhs = formula[i+1:].strip()
                    lv = self.interpret_term(lhs, assignment)
                    rv = self.interpret_term(rhs, assignment)
                    return lv == rv

        # Relationales Atom: "R(t1,...,tn)"
        match = re.fullmatch(r'(\w+)\((.+)\)', formula)
        if match:
            rname = match.group(1)
            args_str = match.group(2)
            if rname in self.signature.relations:
                args = _split_args(args_str)
                evaluated = tuple(self.interpret_term(a, assignment) for a in args)
                return evaluated in self.relations.get(rname, set())

        # Nullstellige Relation (Satz): "P" ohne Argumente
        if formula in self.signature.relations:
            return () in self.relations.get(formula, set())

        # Fallback: unbekannte atomare Formel → False
        return False

    def satisfies(self, formula: str, assignment: dict = None) -> bool:
        """
        Prüft M ⊨ φ[s] nach den Tarski-Wahrheitsbedingungen (rekursiv).

        Unterstützte Formelkonstruktoren:
        - NOT(phi)          – Negation ¬φ
        - AND(phi,psi)      – Konjunktion φ∧ψ
        - OR(phi,psi)       – Disjunktion φ∨ψ
        - IMPLIES(phi,psi)  – Implikation φ→ψ
        - IFF(phi,psi)      – Bikonditional φ↔ψ
        - EXISTS(x,phi)     – Existenzquantor ∃xφ
        - FORALL(x,phi)     – Allquantor ∀xφ
        - Atomare Formeln   – R(t1,...), t1=t2, =(t1,t2)

        @param formula    Formel als String in der oben beschriebenen Syntax
        @param assignment Variablenbelegung s (None → leere Belegung)
        @return           True gdw. M ⊨ φ[s]
        """
        if assignment is None:
            assignment = {}
        formula = formula.strip()

        # ── Negation: NOT(phi) ────────────────────────────────────────────────
        if formula.startswith('NOT(') and formula.endswith(')'):
            inner = formula[4:-1]
            return not self.satisfies(inner, assignment)

        # ── Binäre Konnektive: OP(phi,psi) ───────────────────────────────────
        for prefix, op in [('AND(', 'AND'), ('OR(', 'OR'),
                           ('IMPLIES(', 'IMPLIES'), ('IFF(', 'IFF')]:
            if formula.startswith(prefix) and formula.endswith(')'):
                inner = formula[len(prefix):-1]
                parts = _split_args(inner)
                if len(parts) == 2:
                    left_val = self.satisfies(parts[0], assignment)
                    right_val = self.satisfies(parts[1], assignment)
                    if op == 'AND':
                        return left_val and right_val
                    elif op == 'OR':
                        return left_val or right_val
                    elif op == 'IMPLIES':
                        return (not left_val) or right_val
                    elif op == 'IFF':
                        return left_val == right_val

        # ── Existenzquantor: EXISTS(x,phi) ────────────────────────────────────
        if formula.startswith('EXISTS(') and formula.endswith(')'):
            inner = formula[7:-1]
            parts = _split_args(inner)
            if len(parts) >= 2:
                var = parts[0].strip()
                # Alles nach dem ersten Komma ist die Teilformel
                subformula = ','.join(parts[1:])
                # ∃xφ gilt gdw. ∃a∈A: M⊨φ[s(x↦a)]
                for element in self.universe:
                    new_assignment = dict(assignment)
                    new_assignment[var] = element
                    if self.satisfies(subformula, new_assignment):
                        return True
                return False

        # ── Allquantor: FORALL(x,phi) ─────────────────────────────────────────
        if formula.startswith('FORALL(') and formula.endswith(')'):
            inner = formula[7:-1]
            parts = _split_args(inner)
            if len(parts) >= 2:
                var = parts[0].strip()
                subformula = ','.join(parts[1:])
                # ∀xφ gilt gdw. ∀a∈A: M⊨φ[s(x↦a)]
                for element in self.universe:
                    new_assignment = dict(assignment)
                    new_assignment[var] = element
                    if not self.satisfies(subformula, new_assignment):
                        return False
                return True

        # ── Atomare Formel ────────────────────────────────────────────────────
        return self.satisfies_atomic(formula, assignment)

    def is_model_of(self, theory: list[str]) -> bool:
        """
        Prüft ob M ⊨ T: M erfüllt alle Sätze der Theorie T.

        @param theory Liste von Formeln (Sätze ohne freie Variablen)
        @return       True gdw. M ⊨ φ für alle φ ∈ T
        """
        return all(self.satisfies(axiom, {}) for axiom in theory)

    def elementary_diagram(self) -> list[str]:
        """
        Elementares Diagramm Diag^el(M):
        Alle atomaren Sätze und deren Negationen, die in M wahr bzw. falsch sind.

        Das elementare Diagramm enthält für alle Relationen R und alle Tupel
        (a₁,…,aₙ) aus dem Universum den Satz "R(a₁,…,aₙ)" (falls wahr)
        oder "NOT(R(a₁,…,aₙ))" (falls falsch).

        @return Liste von atomaren Sätzen (als Strings)
        """
        diagram: list[str] = []

        # Gleichheitsdiagramm: a=a für alle Elemente
        for a in self.universe:
            diagram.append(f"=({a},{a})")

        # Für jede Relation alle Tupel prüfen
        for rname, arity in self.signature.relations.items():
            for tup in itertools.product(self.universe, repeat=arity):
                atom = f"{rname}({','.join(str(x) for x in tup)})"
                if tup in self.relations.get(rname, set()):
                    diagram.append(atom)
                else:
                    diagram.append(f"NOT({atom})")

        return diagram

    def cardinality(self) -> int:
        """
        Gibt die Kardinalität |A| des Universums zurück.

        @return Anzahl der Elemente im Universum
        """
        return len(self.universe)

    def __repr__(self) -> str:
        return f"Structure(universe={self.universe}, sig={self.signature})"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Elementare Einbettung
# ═══════════════════════════════════════════════════════════════════════════════

class ElementaryEmbedding:
    """
    Elementare Einbettung f: M → N zwischen L-Strukturen.

    Eine Abbildung f: A → B ist eine elementare Einbettung gdw. für alle
    Formeln φ(x₁,…,xₙ) und alle a₁,…,aₙ ∈ A gilt:
        M ⊨ φ[a₁,…,aₙ]  ⟺  N ⊨ φ[f(a₁),…,f(aₙ)]

    Insbesondere ist jede elementare Einbettung ein injektiver Homomorphismus.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, source: Structure, target: Structure, mapping: dict):
        """
        Initialisiert die Einbettung.

        @param source  Quellstruktur M
        @param target  Zielstruktur N
        @param mapping Dict {a → f(a)} für alle a ∈ A
        """
        self.source = source
        self.target = target
        # Abbildung f: A → B
        self.mapping: dict = dict(mapping)

    def is_homomorphism(self) -> bool:
        """
        Prüft ob f ein Strukturhomomorphismus ist.

        Ein Homomorphismus f: M → N erhält:
        - Relationen: (a₁,…,aₙ) ∈ R^M ⇒ (f(a₁),…,f(aₙ)) ∈ R^N
        - Funktionen: f(f^M(a₁,…,aₙ)) = f^N(f(a₁),…,f(aₙ))
        - Konstanten: f(c^M) = c^N

        @return True gdw. f ein Homomorphismus ist
        """
        # Prüfe Relationstreue: R-Tupel werden auf R-Tupel abgebildet
        for rname, arity in self.source.signature.relations.items():
            src_rel = self.source.relations.get(rname, set())
            tgt_rel = self.target.relations.get(rname, set())
            for tup in src_rel:
                mapped = tuple(self.mapping.get(x, x) for x in tup)
                if mapped not in tgt_rel:
                    return False

        # Prüfe Konstantentreue: c^M wird auf c^N abgebildet
        for cname in self.source.signature.constants:
            src_val = self.source.constants.get(cname)
            tgt_val = self.target.constants.get(cname)
            if src_val is not None and tgt_val is not None:
                if self.mapping.get(src_val) != tgt_val:
                    return False

        return True

    def is_embedding(self) -> bool:
        """
        Prüft ob f eine injektive Einbettung ist (Homomorphismus + injektiv).

        Eine Einbettung (starker Homomorphismus) erhält zusätzlich:
        - Rückwärtsrichtung: (f(a₁),…,f(aₙ)) ∈ R^N ⇒ (a₁,…,aₙ) ∈ R^M

        @return True gdw. f eine Einbettung ist
        """
        # Muss zunächst Homomorphismus sein
        if not self.is_homomorphism():
            return False

        # Injektivität prüfen: f(a) ≠ f(b) für a ≠ b
        mapped_values = list(self.mapping.values())
        if len(mapped_values) != len(set(str(v) for v in mapped_values)):
            return False

        # Rückwärtsrichtung: N-Tupel im Bild müssen in R^M liegen
        img = set(self.mapping.values())
        for rname, arity in self.source.signature.relations.items():
            src_rel = self.source.relations.get(rname, set())
            tgt_rel = self.target.relations.get(rname, set())
            # Invertiere Abbildung
            inv_map = {v: k for k, v in self.mapping.items()}
            for tup in tgt_rel:
                # Prüfe nur Tupel, die vollständig im Bild von f liegen
                if all(x in img for x in tup):
                    preimage = tuple(inv_map[x] for x in tup)
                    if preimage not in src_rel:
                        return False

        return True

    def is_elementary(self, formulas: list[str]) -> bool:
        """
        Prüft ob f für die gegebene Formelliste eine elementare Einbettung ist.

        Für jede Formel φ(x₁,…,xₙ) und alle (a₁,…,aₙ) ∈ A^n:
        M ⊨ φ[a̅]  ⟺  N ⊨ φ[f(a̅)]

        @param formulas Liste von Formeln zum Testen
        @return         True gdw. f für alle gegebenen Formeln elementar ist
        """
        # Iteriere über alle möglichen Variablenbelegungen mit Elementen aus A
        for formula in formulas:
            # Bestimme freie Variablen (vereinfachte Heuristik: einbuchstabige klein)
            free_vars = list(set(re.findall(r'\b([a-z])\b', formula)))
            if not free_vars:
                # Satz (kein freie Var): direkt vergleichen
                src_val = self.source.satisfies(formula, {})
                tgt_val = self.target.satisfies(formula, {})
                if src_val != tgt_val:
                    return False
                continue

            # Teste für alle Belegungen der freien Variablen
            for combo in itertools.product(self.source.universe, repeat=len(free_vars)):
                src_assign = dict(zip(free_vars, combo))
                # Abgebildete Belegung
                tgt_assign = {v: self.mapping.get(a, a) for v, a in src_assign.items()}
                src_val = self.source.satisfies(formula, src_assign)
                tgt_val = self.target.satisfies(formula, tgt_assign)
                if src_val != tgt_val:
                    return False

        return True


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Elementare Äquivalenz
# ═══════════════════════════════════════════════════════════════════════════════

def elementary_equivalence(M: Structure, N: Structure, sentences: list[str]) -> bool:
    """
    Prüft M ≡ N: M und N sind elementar äquivalent.

    Zwei L-Strukturen M und N heißen elementar äquivalent (M ≡ N) gdw.
    für alle L-Sätze φ gilt: M ⊨ φ ⟺ N ⊨ φ.

    Da wir nur endliche Formellisten prüfen können, ist das Ergebnis
    relativ zur gegebenen Formelmenge zu verstehen.

    @param M         Erste L-Struktur
    @param N         Zweite L-Struktur
    @param sentences Liste von L-Sätzen (ohne freie Variablen)
    @return          True gdw. M und N alle gegebenen Sätze gleich beurteilen

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    for sentence in sentences:
        if M.satisfies(sentence, {}) != N.satisfies(sentence, {}):
            return False
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Definierbarkeit
# ═══════════════════════════════════════════════════════════════════════════════

def is_definable(structure: Structure, subset: list, formula_template: str) -> bool:
    """
    Prüft ob eine Teilmenge X ⊆ A durch eine Formel φ(x) definierbar ist.

    X ist durch φ(x) definierbar gdw. für alle a ∈ A gilt:
        a ∈ X  ⟺  M ⊨ φ(a)

    @param structure       L-Struktur M
    @param subset          Teilmenge X ⊆ A (als Liste)
    @param formula_template Formel φ mit genau einer freien Variable 'x'
    @return                True gdw. X = {a ∈ A | M ⊨ φ[x↦a]}

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    subset_set = set(subset)
    # Berechne die durch φ definierte Menge: {a ∈ A | M ⊨ φ[x↦a]}
    defined_set = set()
    for element in structure.universe:
        if structure.satisfies(formula_template, {'x': element}):
            defined_set.add(element)

    return subset_set == defined_set


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Kompaktheitssatz
# ═══════════════════════════════════════════════════════════════════════════════

def compactness_theorem_demo() -> dict:
    """
    Demonstration des Kompaktheitssatzes.

    Kompaktheitssatz (Gödel 1930, Maltsev 1936):
    Eine Theorie T hat ein Modell genau dann wenn jede endliche Teiltheorie
    T₀ ⊆ T ein Modell hat.

    Konstruktion eines Nicht-Standard-Modells der Arithmetik:
    - Beginne mit Th(ℕ) und füge neue Konstante c hinzu.
    - Füge Axiome hinzu: c > 0, c > 1, c > 2, ...
    - Jede endliche Teilmenge der erweiterten Theorie T' hat ein Modell
      (nimm ℕ und interpretiere c als genügend große natürliche Zahl).
    - Nach dem Kompaktheitssatz hat T' ein Modell – dies ist ein
      Nicht-Standard-Modell, in dem c kein Standard-Element ist.

    @return Dict mit Erläuterung und Beispielmodellen

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Axiome der Nicht-Standard-Arithmetik: c > n für n = 0,1,...,N-1
    N = 5
    nonstandard_axioms = [f"GT(c,{n})" for n in range(N)]

    # Signatur der Arithmetik mit Ordnung
    sig = Signature(
        functions={'succ': 1},
        relations={'GT': 2, 'EQ': 2},
        constants=['zero', 'c']
    )

    # Standard-Modell ℕ₅ = {0,1,2,...,4}: endliche Teiltheorien haben Modelle
    finite_models = []
    for k in range(1, N + 1):
        # Endliche Teiltheorie T_k = {c>0, c>1, ..., c>k-1}
        universe_k = list(range(k + 1))
        # c wird als k interpretiert (das größte Element)
        rel_GT_k = {(a, b) for a in universe_k for b in universe_k if a > b}
        rel_EQ_k = {(a, a) for a in universe_k}
        M_k = Structure(
            universe=universe_k,
            signature=sig,
            functions={'succ': lambda x, u=universe_k: min(x + 1, max(u))},
            relations={'GT': rel_GT_k, 'EQ': rel_EQ_k},
            constants={'zero': 0, 'c': k}
        )
        # Prüfe ob M_k die ersten k Axiome erfüllt
        satisfied = all(
            M_k.satisfies(f"GT(c,{n})", {}) for n in range(k)
        )
        finite_models.append({
            'finite_subtheory_size': k,
            'model_universe': universe_k,
            'c_interpretation': k,
            'satisfies_subtheory': satisfied
        })

    return {
        'theorem': 'Kompaktheitssatz',
        'statement': (
            'T hat ein Modell ⟺ jede endliche Teiltheorie T₀ ⊆ T hat ein Modell.'
        ),
        'demonstration': (
            'Nicht-Standard-Arithmetik: Durch Hinzufügen von c>0, c>1, c>2,... '
            'zu Th(ℕ) erhält man eine Theorie T\', deren jede endliche Teilmenge '
            'ein endliches Modell hat. Nach dem Kompaktheitssatz hat T\' ein Modell, '
            'in dem c kein Standard-Element ist (unendlich großes Element).'
        ),
        'nonstandard_axioms': nonstandard_axioms,
        'finite_subtheory_models': finite_models,
        'conclusion': (
            'Da alle endlichen Teiltheorien Modelle haben, garantiert der '
            'Kompaktheitssatz die Existenz eines Nicht-Standard-Modells.'
        )
    }


def nonstandard_arithmetic_demo() -> dict:
    """
    Nicht-Standard-Modell der Arithmetik.

    Durch den Kompaktheitssatz existiert ein Modell der Peano-Arithmetik (PA),
    das ein Element ω* enthält mit ω* > n für alle Standard-Naturzahlen n.

    Dies zeigt die Nicht-Axiomatisierbarkeit von ℕ durch endliche Axiomensätze
    und die Existenz von Modellen jeder unendlichen Kardinalität für Th(ℕ).

    Wir simulieren ein endliches "Nicht-Standard-Modell" mit einem ausgezeichneten
    maximalen Element, das die Rolle von ω* spielt.

    @return Dict mit Erläuterung und Struktur des Nicht-Standard-Modells

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Simuliertes Nicht-Standard-Modell: {0,1,2,3,ω*}
    # ω* > 0,1,2,3 und ω* ist "unendlich"
    standard_part = [0, 1, 2, 3]
    nonstandard_element = 'omega_star'
    universe = standard_part + [nonstandard_element]

    def is_standard(x):
        return isinstance(x, int)

    def gt_nonstandard(a, b):
        """Ordnung im Nicht-Standard-Modell: ω* > alle Standard-Elemente."""
        if a == 'omega_star' and b != 'omega_star':
            return True
        if a != 'omega_star' and b == 'omega_star':
            return False
        if is_standard(a) and is_standard(b):
            return a > b
        return False  # omega_star < omega_star: False

    gt_relation = {(a, b) for a in universe for b in universe if gt_nonstandard(a, b)}

    # Axiome die ω* erfüllt
    nonstandard_axioms = [f"GT(omega_star, {n})" for n in standard_part]

    sig = Signature(
        functions={},
        relations={'GT': 2},
        constants=['omega_star'] + [str(n) for n in standard_part]
    )
    consts = {str(n): n for n in standard_part}
    consts['omega_star'] = 'omega_star'

    M_nonstandard = Structure(
        universe=universe,
        signature=sig,
        functions={},
        relations={'GT': gt_relation},
        constants=consts
    )

    # Prüfe ob ω* größer als alle Standard-Elemente ist
    gt_checks = {
        f"omega_star > {n}": M_nonstandard.satisfies(f"GT(omega_star,{n})", {})
        for n in standard_part
    }

    return {
        'theorem': 'Nicht-Standard-Arithmetik (via Kompaktheitssatz)',
        'universe': universe,
        'nonstandard_element': nonstandard_element,
        'nonstandard_axioms': nonstandard_axioms,
        'nonstandard_element_checks': gt_checks,
        'all_axioms_satisfied': all(gt_checks.values()),
        'explanation': (
            'Das Element ω* ist größer als jede Standard-Zahl. '
            'Durch den Kompaktheitssatz existiert ein solches Element '
            'in einem Modell der vollständigen Arithmetik Th(ℕ). '
            'Solche Modelle haben jede unendliche Kardinalität.'
        )
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Löwenheim-Skolem-Sätze
# ═══════════════════════════════════════════════════════════════════════════════

def lowenheim_skolem_downward(theory_axioms: list[str], model_size: int) -> dict:
    """
    Abwärts-Löwenheim-Skolem-Satz (Löwenheim 1915, Skolem 1920).

    Satz: Besitzt eine (abzählbare) Theorie T ein unendliches Modell,
    so besitzt T auch ein abzählbares Modell.

    Allgemein: Hat T ein Modell der Kardinalität λ ≥ |T|+ℵ₀,
    so hat T ein Modell jeder Kardinalität κ mit ℵ₀ ≤ κ ≤ λ.

    Demonstration: Vollständig geordneter Körper.
    - ℝ (Kardinalität 2^ℵ₀) ist ein Modell der Theorie der dichten linearen Ordnung.
    - Nach dem Abwärts-L-S hat DLO auch ein abzählbares Modell → ℚ!
    - Das zeigt: ℚ ≡ ℝ bezüglich DLO-Axiomen (aber nicht isomorph!).

    @param theory_axioms Liste der Theorie-Axiome (Beschreibungen)
    @param model_size    Gewünschte Modellgröße (für endliche Demo)
    @return              Dict mit Erläuterung und Modell

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Demo: DLO (Dense Linear Order without Endpoints) über ℚ∩[0,1]
    # (endliche Approximation mit Brüchen)
    from fractions import Fraction

    # Erzeuge ein abzählbares dicht geordnetes Modell: Brüche k/n für kleine k,n
    dlo_universe = sorted(set(
        Fraction(k, d) for d in range(1, model_size + 1)
        for k in range(1, d)
    ))

    if not dlo_universe:
        dlo_universe = [Fraction(1, 2)]

    sig_dlo = Signature(
        functions={},
        relations={'LT': 2},
        constants=[]
    )
    lt_rel = {(a, b) for a in dlo_universe for b in dlo_universe if a < b}

    M_dlo = Structure(
        universe=dlo_universe,
        signature=sig_dlo,
        functions={},
        relations={'LT': lt_rel},
        constants={}
    )

    # DLO-Axiome (als natürlichsprachliche Beschreibungen):
    dlo_axioms_descriptions = [
        'Irreflexivität: ∀x ¬(x < x)',
        'Transitivität: ∀x∀y∀z (x<y ∧ y<z → x<z)',
        'Dichtheit: ∀x∀y (x<y → ∃z(x<z ∧ z<y))',
        'Kein kleinstes Element: ∀x ∃y (y<x)',
        'Kein größtes Element: ∀x ∃y (x<y)',
    ]

    # Prüfe Dichtheit exemplarisch: zwischen je zwei Elementen gibt es ein weiteres
    density_check = True
    for a, b in itertools.islice(itertools.combinations(dlo_universe, 2), 20):
        # Suche c mit a < c < b
        found_between = any(a < c < b for c in dlo_universe)
        if not found_between and b - a > Fraction(1, model_size):
            density_check = False
            break

    return {
        'theorem': 'Abwärts-Löwenheim-Skolem',
        'statement': (
            'Hat T ein unendliches Modell, so hat T ein abzählbares Modell. '
            '(Allg.: für jede κ mit ℵ₀ ≤ κ ≤ |M|)'
        ),
        'example_theory': 'DLO (Dense Linear Order without Endpoints)',
        'large_model': 'ℝ (Kardinalität 2^ℵ₀)',
        'small_model': 'ℚ (abzählbar, ≡ ℝ bezüglich DLO)',
        'finite_approximation_size': len(dlo_universe),
        'universe_sample': [str(x) for x in dlo_universe[:10]],
        'dlo_axioms': dlo_axioms_descriptions,
        'density_check_passed': density_check,
        'theory_axioms_given': theory_axioms,
        'corollary': (
            'ℚ ≡ ℝ (DLO): Beide Strukturen sind elementar äquivalent, '
            'aber nicht isomorph (ℝ vollständig, ℚ nicht).'
        )
    }


def lowenheim_skolem_upward(theory_axioms: list[str], base_model_size: int) -> dict:
    """
    Aufwärts-Löwenheim-Skolem-Satz (Tarski-Vaught 1936).

    Satz: Besitzt eine Theorie T ein unendliches Modell der Kardinalität κ,
    so besitzt T Modelle jeder Kardinalität λ ≥ κ.

    Beweis-Idee: Füge λ neue Konstantensymbole {c_α | α < λ} hinzu
    und neue Axiome c_α ≠ c_β (für α ≠ β). Jede endliche Teilmenge
    der erweiterten Theorie hat ein Modell (da das Ausgangsmodell unendlich ist).
    Nach dem Kompaktheitssatz hat die gesamte Theorie ein Modell der Größe ≥ λ.

    @param theory_axioms    Axiome der Ausgangstheorie (Beschreibungen)
    @param base_model_size  Größe des Ausgangsmodells (als Untergrenze)
    @return                 Dict mit Erläuterung

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Demonstration: ℕ als abzählbares Modell, Konstruktion größerer Modelle
    target_sizes = [base_model_size, base_model_size * 2, base_model_size * 10]

    models_info = []
    for size in target_sizes:
        models_info.append({
            'cardinality': size,
            'description': f'Modell der Kardinalität {size} (konstruiert durch Hinzufügen neuer Elemente)',
            'method': 'Ultrapotenz oder Kompaktheit + neue Konstanten'
        })

    return {
        'theorem': 'Aufwärts-Löwenheim-Skolem (Tarski-Vaught)',
        'statement': (
            'Hat T ein unendliches Modell der Kardinalität κ, '
            'so hat T Modelle jeder Kardinalität λ ≥ κ.'
        ),
        'base_model_size': base_model_size,
        'proof_sketch': (
            'Füge λ neue Konstantensymbole c_α hinzu mit Axiomen c_α ≠ c_β. '
            'Jede endliche Teilmenge hat ein Modell (Ausgangsmodell ist unendlich). '
            'Kompaktheitssatz: Die erweiterte Theorie hat ein Modell der Größe ≥ λ.'
        ),
        'target_models': models_info,
        'theory_axioms_given': theory_axioms,
        'consequence': (
            'Keine abzählbare Theorie mit einem unendlichen Modell ist kategorisch '
            'bezüglich aller Kardinalitäten → "Skolems Paradoxon" und L-S-Schranken.'
        )
    }


def skolem_paradox_explanation() -> dict:
    """
    Das Skolem-Paradoxon (1922).

    Paradoxon:
    - ZF-Mengenlehre postuliert die Existenz überabzählbarer Mengen (z.B. ℝ).
    - Der Abwärts-Löwenheim-Skolem-Satz besagt: Hat ZF ein Modell,
      so hat ZF auch ein abzählbares Modell M.
    - In M gibt es keine Bijektion zwischen ω^M und (ℝ)^M – obwohl M selbst abzählbar ist!

    Auflösung:
    - "Überabzählbar" ist relativ zum Modell.
    - Im abzählbaren Modell M fehlt die Bijektion als Element von M,
      obwohl sie außerhalb von M existiert (von außen: M ist abzählbar).
    - Der Begriff "überabzählbar" bedeutet: "keine Bijektion existiert innerhalb des Modells".
    - Absolute Überabzählbarkeit erfordert den "richtigen" Mengenbegriff von außen.

    @return Dict mit Erläuterung des Paradoxons und seiner Auflösung

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    return {
        'name': 'Skolem-Paradoxon',
        'year': 1922,
        'paradox': (
            'ZF-Mengenlehre hat (nach Abwärts-L-S) ein abzählbares Modell M, '
            'obwohl ZF die Existenz überabzählbarer Mengen beweist (Cantors Satz: |P(ω)| > |ω|).'
        ),
        'resolution': {
            'key_insight': 'Überabzählbarkeit ist relativ zum Modell',
            'explanation': (
                'Im abzählbaren Modell M ist (ℝ)^M eine Menge, aber keine Bijektion '
                'ω^M → (ℝ)^M liegt als Element in M. Von außen betrachtet ist M selbst '
                'abzählbar und (ℝ)^M ist abzählbar – aber diese Bijektion ist kein '
                'Element von M. Das Modell M "sieht" keine Bijektion und urteilt daher '
                'korrekt: "(ℝ)^M ist überabzählbar" (relativ zu M).'
            ),
            'technical': (
                'Überabzählbarkeit von X in M: In M existiert kein f: ω → X surjektiv. '
                'Die Bijektion mag extern existieren, aber nicht als Element von M. '
                'Das ist kein Widerspruch, sondern zeigt die Relativität mathematischer Begriffe.'
            )
        },
        'philosophical_implications': [
            'Mathematische Wahrheit ist relativ zu einem Modell/einer Interpretation.',
            'Kein abzählbares Axiomensystem kann ℕ oder ℝ bis auf Isomorphie charakterisieren.',
            'Der Mengenbegriff ist von innen und außen verschieden (innere vs. äußere Abzählbarkeit).',
            'Motiviert die Suche nach starken Axiomen (z.B. große Kardinalzahlaxiome).'
        ]
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Typentheorie
# ═══════════════════════════════════════════════════════════════════════════════

class Type:
    """
    p-Typ über einer Theorie T.

    Ein n-Typ über T ist eine konsistente Menge Γ(x₁,…,xₙ) von Formeln,
    sodass T ∪ Γ konsistent ist. Ein vollständiger n-Typ ist maximal konsistent.

    In einer Struktur M wird ein Typ p(x) durch ein Tupel a̅ ∈ A^n realisiert,
    wenn M ⊨ φ[a̅] für alle φ ∈ p.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, formulas: list[str], theory: list[str]):
        """
        Initialisiert den Typ.

        @param formulas Liste von Formeln φ(x₁,…,xₙ) (Typformeln)
        @param theory   Hintergrundtheorie T
        """
        # Typformeln: die Menge der Formeln im Typ
        self.formulas: list[str] = list(formulas)
        # Hintergrundtheorie
        self.theory: list[str] = list(theory)

    def is_complete(self, all_sentences: list[str]) -> bool:
        """
        Prüft ob der Typ vollständig (maximal) ist.

        Ein Typ p ist vollständig gdw. für jede Formel φ aus der Liste gilt:
        entweder φ ∈ p oder ¬φ ∈ p (bzw. NOT(φ) in unserer Notation).

        @param all_sentences Alle zu prüfenden Formeln
        @return              True gdw. p für alle φ entschieden ist
        """
        formula_set = set(self.formulas)
        for sentence in all_sentences:
            # Prüfe ob φ oder NOT(φ) im Typ enthalten ist
            has_positive = sentence in formula_set
            has_negative = f"NOT({sentence})" in formula_set
            if not (has_positive or has_negative):
                return False
        return True

    def is_realized(self, structure: Structure) -> bool:
        """
        Prüft ob der Typ in der Struktur M realisiert wird.

        p(x) wird in M realisiert gdw. ∃a ∈ A: M ⊨ φ[x↦a] für alle φ ∈ p.

        @param structure L-Struktur M
        @return          True gdw. ∃a ∈ A das alle Typformeln erfüllt
        """
        for element in structure.universe:
            # Prüfe ob dieses Element alle Typformeln erfüllt
            if all(structure.satisfies(formula, {'x': element})
                   for formula in self.formulas):
                return True
        return False

    def is_omitted(self, structure: Structure) -> bool:
        """
        Prüft ob der Typ in der Struktur M ausgelassen (omitted) wird.

        p(x) wird ausgelassen gdw. kein Element von M alle Formeln von p erfüllt.

        @param structure L-Struktur M
        @return          True gdw. p in M nicht realisiert ist
        """
        return not self.is_realized(structure)

    def __repr__(self) -> str:
        return f"Type(formulas={self.formulas})"


def type_space(theory: list[str], variables: list[str]) -> dict:
    """
    Typraum S_n(T): Beschreibung aller vollständigen n-Typen über T.

    Der Typraum S_n(T) trägt die Stone-Topologie:
    - Basismengen: [φ] = {p ∈ S_n(T) | φ ∈ p} für jede Formel φ
    - S_n(T) ist kompakt (Kompaktheitssatz)
    - S_n(T) ist Hausdorff (da totalgetrennt)
    - S_n(T) ist perfekt (kein isolierter Punkt bei vollständigen T)

    @param theory    Vollständige Theorie T
    @param variables Liste der freien Variablen x₁,…,xₙ
    @return          Dict mit Beschreibung des Typraums

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n = len(variables)
    return {
        'space': f'S_{n}(T)',
        'variables': variables,
        'theory': theory,
        'topology': 'Stone-Topologie',
        'basis': {
            'description': 'Basismengen [φ] = {p ∈ S_n(T) | φ ∈ p}',
            'example': '[x=y] enthält alle Typen, die x=y fordern'
        },
        'properties': {
            'compact': True,
            'hausdorff': True,
            'reason': 'Kompaktheitssatz ⇒ kompakt; totalgetrennt ⇒ Hausdorff'
        },
        'n_tuples': n,
        'description': (
            f'S_{n}(T) ist der Raum aller vollständigen {n}-Typen über T. '
            f'Jedes Tupel (a₁,…,a_{n}) aus einem Modell M von T realisiert '
            f'genau einen Typ tp^M(a₁,…,a_{n}) ∈ S_{n}(T).'
        )
    }


def omitting_types_theorem_demo() -> dict:
    """
    Auslassungstypen-Satz (Henkin 1954, Omitting Types Theorem).

    Satz: Sei T eine vollständige abzählbare Theorie und p(x) ein nicht-isolierter
    n-Typ über T (d.h. kein θ(x) mit T⊨θ→φ für alle φ∈p und T⊨∃xθ).
    Dann gibt es ein abzählbares Modell M von T, das p auslässt.

    Konsequenz:
    - Nicht alle Typen müssen in jedem Modell realisiert sein.
    - Atomare Modelle: realisieren nur isolierte Typen.
    - Saturierte Modelle: realisieren alle Typen.

    @return Dict mit Erläuterung

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Beispiel: Theorie der dichten linearen Ordnung (DLO)
    # Typ p(x) = {x > n | n ∈ ℕ} ist nicht isoliert (kein oberes Maximum)
    # Im Standardmodell ℚ wird dieser Typ nicht realisiert (ℚ ist Archimedisch)
    # aber in einem Nicht-Standard-Modell schon.

    example_type_formulas = [f"GT(x,{n})" for n in range(5)]

    # Konstruiere eine Struktur (endliche Approximation), die diesen Typ auslässt
    sig = Signature(functions={}, relations={'GT': 2}, constants=[])
    universe = list(range(10))
    gt_rel = {(a, b) for a in universe for b in universe if a > b}
    M_omitting = Structure(
        universe=universe,
        signature=sig,
        functions={},
        relations={'GT': gt_rel},
        constants={}
    )

    example_type = Type(formulas=example_type_formulas, theory=['DLO-Axiome'])

    # In diesem endlichen Modell ist der Typ p(x) = {x>0,x>1,...,x>4} realisiert durch 5,6,...
    is_realized_here = example_type.is_realized(M_omitting)

    # Modell das den Typ auslässt: Universum = {0,1,2,3,4}
    universe_small = list(range(5))
    gt_rel_small = {(a, b) for a in universe_small for b in universe_small if a > b}
    M_omits = Structure(
        universe=universe_small,
        signature=sig,
        functions={},
        relations={'GT': gt_rel_small},
        constants={}
    )
    is_omitted_here = example_type.is_omitted(M_omits)

    return {
        'theorem': 'Auslassungstypen-Satz (Omitting Types Theorem)',
        'statement': (
            'Sei T vollständig abzählbar und p(x) ein nicht-isolierter Typ. '
            'Dann existiert ein abzählbares Modell M von T, das p auslässt.'
        ),
        'example_type': example_type_formulas,
        'realizing_structure': {
            'universe_size': len(universe),
            'realizes_type': is_realized_here,
            'witnesses': [a for a in universe
                          if all(M_omitting.satisfies(f"GT(x,{n})", {'x': a})
                                 for n in range(5))]
        },
        'omitting_structure': {
            'universe_size': len(universe_small),
            'omits_type': is_omitted_here,
        },
        'model_types': {
            'atomic_model': 'Realisiert nur isolierte Typen (minimal)',
            'saturated_model': 'Realisiert alle Typen seiner Kardinalität (maximal)',
            'prime_model': 'Einbettbar in jedes Modell der gleichen Theorie'
        }
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Kategorien und Vollständigkeit
# ═══════════════════════════════════════════════════════════════════════════════

def theory_is_categorical(theory_name: str, cardinality: str) -> dict:
    """
    Untersucht die κ-Kategorik einer Theorie.

    Eine Theorie T heißt κ-kategorisch gdw. T (bis auf Isomorphie)
    genau ein Modell der Kardinalität κ besitzt.

    Bekannte Beispiele:
    - DLO ist ℵ₀-kategorisch (Cantor 1895):
      Alle abzählbaren dichten linearen Ordnungen ohne Endpunkte ≅ (ℚ,<).
    - Th(ℤ,+,0) ist nicht ℵ₀-kategorisch (verschiedene abzählbare Modelle).
    - ACF_p (algebraisch abgeschlossene Körper Charakteristik p) ist
      κ-kategorisch für alle überabzählbaren κ.
    - DLO ist nicht überabzählbar-kategorisch (verschiedene überabzählbare Dichten).

    @param theory_name  Name der Theorie (z.B. 'DLO', 'ACF0', 'PA')
    @param cardinality  Kardinalität als String (z.B. 'aleph_0', 'aleph_1', 'uncountable')
    @return             Dict mit Kategorik-Information

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Bekannte Kategorik-Tabelle
    known_categoricity: dict[str, dict[str, bool | str]] = {
        'DLO': {
            'aleph_0': True,
            'uncountable': False,
            'explanation': (
                'ℵ₀-kategorisch nach Cantor: Alle abzählbaren dichten linearen '
                'Ordnungen ohne Endpunkte sind isomorph zu (ℚ,<). '
                'Nicht überabzählbar-kategorisch: (ℝ,<) und (ℝ\\ℚ∪ℚ·√2,<) nicht isomorph.'
            ),
            'unique_countable_model': '(ℚ, <)'
        },
        'ACF0': {
            'aleph_0': False,
            'uncountable': True,
            'explanation': (
                'Nicht ℵ₀-kategorisch: Verschiedene abzählbare algebraisch abgeschlossene '
                'Körper der Char. 0 (verschiedene Transzendenzgrade über ℚ). '
                'Überabzählbar kategorisch: Kardinalität bestimmt Transzendenzgrad eindeutig.'
            ),
            'unique_countable_model': 'Nein (Transzendenzgrad variiert)'
        },
        'ACFp': {
            'aleph_0': False,
            'uncountable': True,
            'explanation': (
                'Wie ACF0: nicht ℵ₀-kategorisch, aber κ-kategorisch für alle überabzählbaren κ.'
            ),
            'unique_countable_model': 'Nein'
        },
        'PA': {
            'aleph_0': False,
            'uncountable': False,
            'explanation': (
                'Peano-Arithmetik ist weder ℵ₀-kategorisch noch überabzählbar-kategorisch. '
                'Nicht-Standard-Modelle existieren in allen unendlichen Kardinalitäten.'
            ),
            'unique_countable_model': 'Nein (Nicht-Standard-Modelle)'
        },
        'RCF': {
            'aleph_0': False,
            'uncountable': False,
            'explanation': (
                'Reell abgeschlossene Körper (RCF): nicht kategorisch in irgendeiner '
                'Kardinalität, aber vollständig (durch Quantorenelimination).'
            ),
            'unique_countable_model': 'Nein'
        }
    }

    if theory_name in known_categoricity:
        info = known_categoricity[theory_name]
        cat_in_card = info.get(cardinality, 'unbekannt')
        return {
            'theory': theory_name,
            'cardinality': cardinality,
            'is_categorical': cat_in_card,
            'explanation': info.get('explanation', ''),
            'unique_model': info.get('unique_countable_model', 'unbekannt')
        }
    else:
        return {
            'theory': theory_name,
            'cardinality': cardinality,
            'is_categorical': 'unbekannt',
            'explanation': f'Theorie {theory_name!r} nicht in der Datenbank.',
            'unique_model': 'unbekannt'
        }


def vaught_test(theory_name: str) -> dict:
    """
    Vaught-Kriterium für Vollständigkeit (Vaught 1954).

    Satz: Eine Theorie T ohne endliche Modelle, die κ-kategorisch für
    mindestens eine unendliche Kardinalität κ ist, ist vollständig.

    Beweis-Idee: Sind M₁, M₂ ⊨ T (unendlich) und T κ-kategorisch,
    so gibt es elementar überabzählbare N₁ ≅ N₂ (beide Kardinalität κ).
    Da elementare Erweiterungen dieselbe vollständige Theorie haben:
    Th(M₁) = Th(N₁) = Th(N₂) = Th(M₂) → T ist vollständig.

    @param theory_name Name der Theorie
    @return            Dict mit Vollständigkeitsinformation

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Bekannte Vollständigkeits-Information
    known_complete: dict[str, dict] = {
        'DLO': {
            'complete': True,
            'has_finite_models': False,
            'categorical_in': 'ℵ₀',
            'proof': 'ℵ₀-kategorisch + keine endlichen Modelle → vollständig (Vaught).',
            'decidable': True
        },
        'ACF0': {
            'complete': True,
            'has_finite_models': False,
            'categorical_in': 'alle überabzählbaren κ',
            'proof': 'Überabzählbar kategorisch + keine endlichen Modelle → vollständig.',
            'decidable': True
        },
        'PA': {
            'complete': False,
            'has_finite_models': False,
            'categorical_in': 'keine',
            'proof': 'Gödels Unvollständigkeitssatz: PA ist unvollständig.',
            'decidable': False
        },
        'RCF': {
            'complete': True,
            'has_finite_models': False,
            'categorical_in': 'keine (aber vollständig durch QE)',
            'proof': 'Tarski: RCF ist vollständig durch Quantorenelimination.',
            'decidable': True
        }
    }

    if theory_name in known_complete:
        info = known_complete[theory_name]
        vaught_applies = (
            not info['has_finite_models'] and
            info['categorical_in'] not in ['keine', '']
        )
        return {
            'theory': theory_name,
            'vaught_criterion_applies': vaught_applies,
            'complete': info['complete'],
            'has_finite_models': info['has_finite_models'],
            'categorical_in': info['categorical_in'],
            'proof': info['proof'],
            'decidable': info['decidable']
        }
    else:
        return {
            'theory': theory_name,
            'vaught_criterion_applies': 'unbekannt',
            'complete': 'unbekannt',
            'explanation': f'Theorie {theory_name!r} nicht in der Datenbank.'
        }


def morley_theorem_demo() -> dict:
    """
    Morleys Kategorien-Satz (1965).

    Satz (Morley 1965): Sei T eine abzählbare vollständige Theorie.
    Ist T κ-kategorisch für ein überabzählbares κ,
    so ist T λ-kategorisch für alle überabzählbaren λ.

    Dies ist eines der tiefsten Ergebnisse der Modelltheorie.
    Es führte zur Entwicklung der Stabilitätstheorie (Shelah).

    Stabilitätshierarchie (Shelah):
    - ω-stabil: Typraum über jedem abzählbaren A ist abzählbar
    - superstabil: stabil in allen λ ≥ 2^ℵ₀
    - stabil: stabil in gewissen Kardinalitäten
    - instabil: nicht stabil in irgendwelchen Kardinalitäten

    @return Dict mit Erläuterung des Morley-Satzes

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    return {
        'theorem': 'Morleys Kategorien-Satz',
        'year': 1965,
        'statement': (
            'Eine abzählbare vollständige Theorie T, die für ein überabzählbares κ '
            'κ-kategorisch ist, ist λ-kategorisch für alle überabzählbaren λ.'
        ),
        'examples': {
            'ACF0': {
                'categorical_in_uncountable': True,
                'all_uncountable_categorical': True,
                'reason': 'Überabzählbare ACF0-Modelle durch Transzendenzgrad bestimmt.'
            },
            'DLO': {
                'categorical_in_uncountable': False,
                'all_uncountable_categorical': False,
                'reason': 'DLO nur ℵ₀-kategorisch, nicht überabzählbar.'
            }
        },
        'shelah_stability': {
            'omega_stable': ['ACF0', 'ACFp'],
            'superstable': ['Th(ℤ, +)'],
            'stable_not_superstable': ['Th(ℤ, +, <)'],
            'unstable': ['DLO', 'Th(ℚ, +, <)', 'PA'],
            'explanation': (
                'ω-stabile Theorien sind genau die Morleys Kategorien-Satz anwendbaren Theorien '
                '(alle überabzählbar-kategorialen Theorien sind ω-stabil).'
            )
        },
        'significance': (
            'Der Morley-Satz motivierte Shelahs Stabilitätstheorie (Classification Theory), '
            'die die Modelltheorie revolutionierte und zeigt, wann Theorien "zahm" (stabil) '
            'oder "wild" (instabil) sind.'
        )
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 10. Quantorenelimination
# ═══════════════════════════════════════════════════════════════════════════════

def quantifier_elimination_demo(theory: str) -> dict:
    """
    Demonstration der Quantorenelimination.

    Eine Theorie T hat Quantorenelimination gdw. jede Formel φ äquivalent zu
    einer quantorenfreien Formel ψ ist (relativ zu T):
        T ⊨ φ ↔ ψ,    ψ quantorenfrei.

    Bekannte Theorien mit Quantorenelimination:
    1. DLO (dichte lineare Ordnung ohne Endpunkte, Langford 1927)
       - ∃y(x₁ < y ∧ y < x₂)  ≡  x₁ < x₂
       - Beweis durch Back-and-Forth-Methode
    2. Th(ℚ, <) = DLO (vollständig durch QE)
    3. ACF (algebraisch abgeschlossene Körper, Tarski 1948)
       - Eliminierung durch Resultanten
       - ∃y(p(y)=0) eliminierbar mit Grundlage Nullstellensatz
    4. RCF (reell abgeschlossene Körper, Tarski 1951)
       - ∃y(p(y)=0 ∧ y > 0) eliminierbar
       - Begründet Tarskis Entscheidbarkeit der reellen Arithmetik

    @param theory  Theoriename ('DLO', 'ACF', 'RCF', 'Presburger')
    @return        Dict mit QE-Information und Beispielen

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    qe_theories: dict[str, dict] = {
        'DLO': {
            'has_qe': True,
            'full_name': 'Dense Linear Order without Endpoints',
            'language': '(<)',
            'examples': [
                {
                    'with_quantifier': '∃y(x₁ < y ∧ y < x₂)',
                    'quantifier_free': 'x₁ < x₂',
                    'explanation': 'Ein Zwischenpunkt existiert gdw. x₁ < x₂ (Dichtheit).'
                },
                {
                    'with_quantifier': '∃y(y < x)',
                    'quantifier_free': '⊤ (wahr)',
                    'explanation': 'Kein kleinstes Element: links von jedem Punkt gibt es einen.'
                }
            ],
            'method': 'Back-and-Forth (Ehrenfeucht-Fraïssé-Spiele)',
            'consequence': 'DLO ist vollständig und entscheidbar.',
            'historical': 'Langford (1927), Cantor (1895 für ℵ₀-Kategorik)'
        },
        'ACF': {
            'has_qe': True,
            'full_name': 'Algebraically Closed Fields',
            'language': '(+, ·, 0, 1)',
            'examples': [
                {
                    'with_quantifier': '∃y(y² - x = 0)',
                    'quantifier_free': '⊤ (in ACF: jedes Element hat eine Quadratwurzel)',
                    'explanation': 'Algebraisch abgeschlossen: jedes Polynom hat eine Nullstelle.'
                },
                {
                    'with_quantifier': '∃y(p(y) = 0)',
                    'quantifier_free': 'Diskriminanten-Bedingung über Koeffizienten',
                    'explanation': 'Resultante eliminiert Existenzquantor über Polynomwurzel.'
                }
            ],
            'method': 'Resultanten, Nullstellensatz',
            'consequence': 'ACFp ist vollständig (für jede Charakteristik p) und entscheidbar.',
            'historical': 'Tarski (1948)'
        },
        'RCF': {
            'has_qe': True,
            'full_name': 'Real Closed Fields',
            'language': '(+, ·, 0, 1, <)',
            'examples': [
                {
                    'with_quantifier': '∃y(y² = x)',
                    'quantifier_free': 'x ≥ 0',
                    'explanation': 'Quadratwurzel existiert gdw. Argument nicht-negativ.'
                },
                {
                    'with_quantifier': '∃y(y³ + py + q = 0)',
                    'quantifier_free': 'Diskriminante 4p³ + 27q² ≤ 0',
                    'explanation': 'Kubische Gleichung hat reelle Lösung gdw. Diskriminante ≤ 0.'
                }
            ],
            'method': 'Cylindrical Algebraic Decomposition (Collins 1975)',
            'consequence': 'Reelle Arithmetik ist entscheidbar (Tarskis Theorem).',
            'historical': 'Tarski (1951)'
        },
        'Presburger': {
            'has_qe': True,
            'full_name': 'Presburger Arithmetic Th(ℤ, +, <)',
            'language': '(+, <, 0, 1)',
            'examples': [
                {
                    'with_quantifier': '∃y(x = 2·y)',
                    'quantifier_free': '2 | x  (x ist gerade)',
                    'explanation': 'Existenz von y mit 2y=x ⟺ 2 teilt x.'
                },
                {
                    'with_quantifier': '∃y(x = ny)',
                    'quantifier_free': 'n | x',
                    'explanation': 'Teilbarkeit ist quantorenfrei ausdrückbar.'
                }
            ],
            'method': 'Omega-Test (Pugh 1991), Presburger-Elimination',
            'consequence': 'Presburger-Arithmetik ist vollständig und entscheidbar (doppelt-exponentiell).',
            'historical': 'Presburger (1929)'
        }
    }

    if theory in qe_theories:
        info = qe_theories[theory]
        return {
            'theory': theory,
            'has_quantifier_elimination': info['has_qe'],
            'full_name': info['full_name'],
            'language': info['language'],
            'examples': info['examples'],
            'elimination_method': info['method'],
            'consequence': info['consequence'],
            'historical_reference': info['historical'],
            'definition': (
                'T hat Quantorenelimination gdw. für jede Formel φ eine '
                'quantorenfreie Formel ψ existiert mit T ⊨ φ ↔ ψ.'
            )
        }
    else:
        return {
            'theory': theory,
            'has_quantifier_elimination': 'unbekannt',
            'explanation': f'Theorie {theory!r} nicht in der Datenbank.',
            'known_qe_theories': list(qe_theories.keys())
        }


# ═══════════════════════════════════════════════════════════════════════════════
# Hilfsfunktionen
# ═══════════════════════════════════════════════════════════════════════════════

def _split_args(s: str) -> list[str]:
    """
    Teilt einen Argumentstring korrekt an Kommas auf,
    unter Berücksichtigung von verschachtelten Klammern.

    Beispiel: "f(x,y),z" → ["f(x,y)", "z"]

    @param s  Argumentstring
    @return   Liste der einzelnen Argumente (getrimmt)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    args: list[str] = []
    depth = 0
    current: list[str] = []

    for ch in s:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == ',' and depth == 0:
            # Trennstelle auf Tiefe 0: neues Argument beginnt
            args.append(''.join(current).strip())
            current = []
        else:
            current.append(ch)

    # Letztes Argument hinzufügen
    if current:
        args.append(''.join(current).strip())

    return args
