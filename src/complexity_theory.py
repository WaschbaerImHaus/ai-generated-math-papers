"""
@file complexity_theory.py
@brief Komplexitätstheorie — P vs NP, Turingmaschinen, NP-Vollständigkeit,
       Schaltkreiskomplexität, Orakel-Barrieren, Quantenkomplexität.
@description
    Dieses Modul implementiert die wichtigsten Konzepte der Berechnungskomplexitäts-
    theorie. Es dient als Lernwerkzeug und kann konkrete Berechnungen (z.B. DPLL-SAT,
    Turingmaschinen-Simulation) sowie theoretische Metadaten (Komplexitätsklassen,
    bekannte Reduktionen, Barrieren) liefern.

    Klassen:
        TuringMachine            — deterministische TM mit Band-Simulation
        NonDeterministicTM       — zertifikat-basierte NTM
        ComplexityClass          — Metadaten über Komplexitätsklassen (P, NP, …)
        PolynomialHierarchy      — Polynomielle Hierarchie Σ/Π/Δ
        NPCompletenessTheory     — Cook-Levin, Karp-21, Reduktionen
        SATSolver                — DPLL, CDCL-Schritt, Resolution, Tseitin
        CircuitComplexity        — Schaltkreisgröße/-tiefe, NC, Razborov-Rudich
        OracleComplexity         — Baker-Gill-Solovay, Relativierung
        ProofComplexity          — Resolution, Frege, Natural Proofs, Algebrization, GCT
        QuantumComplexity        — BQP, Grover, Shor, Quanten-PCP

@author Michael Fuhrmann
@lastModified 2026-03-11
"""

from __future__ import annotations

import itertools
import math
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _is_tautology(clause: List[int], assignment: Dict[int, bool]) -> bool:
    """
    Prüft, ob eine Klausel (Liste von Literalen) unter einer (Teil-)Belegung
    erfüllt ist.  Literal l > 0 entspricht Variable l, l < 0 entspricht ¬|l|.
    """
    for lit in clause:
        var = abs(lit)
        if var in assignment:
            val = assignment[var]
            if (lit > 0 and val) or (lit < 0 and not val):
                return True
    return False


def _unit_propagate(
    clauses: List[List[int]], assignment: Dict[int, bool]
) -> Tuple[List[List[int]], Dict[int, bool], bool]:
    """
    Führt Unit Propagation durch:  Wenn eine Klausel nur noch ein ungebundenes
    Literal hat, muss dieses True sein.  Gibt (neue_klauseln, neue_belegung,
    konflikt) zurück.
    """
    changed = True
    assignment = dict(assignment)  # Kopie
    while changed:
        changed = False
        new_clauses: List[List[int]] = []
        for clause in clauses:
            # Unerfüllte Literale in dieser Klausel ermitteln
            free = [l for l in clause if abs(l) not in assignment]
            if _is_tautology(clause, assignment):
                continue  # Klausel bereits erfüllt
            if not free:
                return clauses, assignment, True  # Konflikt: leere Klausel
            if len(free) == 1:
                # Unit: einziges freies Literal erzwingen
                lit = free[0]
                assignment[abs(lit)] = lit > 0
                changed = True
            else:
                new_clauses.append(clause)
        if changed:
            clauses = new_clauses
    return clauses, assignment, False


# ===========================================================================
# 1. Turingmaschinen
# ===========================================================================

class TuringMachine:
    """
    Deterministische Turingmaschine mit endlichem Zustandsübergangs-
    diagramm.  Das Band ist beidseitig unendlich (als dict implementiert).

    Übergangsformat:
        transitions[(zustand, symbol)] = (neuer_zustand, schreib_symbol, richtung)
        Richtung: 'L' (links), 'R' (rechts), 'N' (kein Schritt)
    """

    BLANK = '_'   # Blank-Symbol

    def __init__(
        self,
        states: List[str],
        alphabet: List[str],
        transitions: Dict[Tuple[str, str], Tuple[str, str, str]],
        initial_state: str,
        accept_states: List[str],
        reject_states: Optional[List[str]] = None,
    ):
        """
        @param states         Menge aller Zustände
        @param alphabet       Bandalphabet (ohne Blank)
        @param transitions    Übergangsfunktion als Dictionary
        @param initial_state  Startzustand
        @param accept_states  Akzeptierzustände
        @param reject_states  Ablehnungszustände (optional; Rest = unbekannt)
        """
        self.states = set(states)
        self.alphabet = set(alphabet) | {self.BLANK}
        self.transitions = transitions
        self.initial_state = initial_state
        self.accept_states = set(accept_states)
        self.reject_states = set(reject_states) if reject_states else set()

    def simulate(
        self, tape_input: str, max_steps: int = 10_000
    ) -> Tuple[bool, str, int]:
        """
        Simuliert die TM auf einem Eingabeband.

        @param tape_input  Eingabestring (wird auf Band geschrieben)
        @param max_steps   Maximale Schrittzahl (Schutz vor Endlosschleifen)
        @return (halted, output_tape, steps)
                halted      = True wenn Akzeptier-/Ablehnzustand erreicht
                output_tape = Bandinhalt ohne führende/abschließende Blanks
                steps       = Anzahl ausgeführter Schritte
        """
        # Band initialisieren: Position → Symbol
        band: Dict[int, str] = {i: c for i, c in enumerate(tape_input)}
        kopf = 0          # Kopfposition
        zustand = self.initial_state
        steps = 0

        while steps < max_steps:
            # Haltezustand erreicht?
            if zustand in self.accept_states or zustand in self.reject_states:
                break

            # Aktuelles Symbol lesen (Blank wenn nicht vorhanden)
            symbol = band.get(kopf, self.BLANK)

            # Übergang nachschlagen
            key = (zustand, symbol)
            if key not in self.transitions:
                # Keine Übergangsregel → implizites Ablehnen
                break

            neuer_zustand, schreib_symbol, richtung = self.transitions[key]
            band[kopf] = schreib_symbol
            zustand = neuer_zustand

            # Kopf bewegen
            if richtung == 'R':
                kopf += 1
            elif richtung == 'L':
                kopf -= 1
            # 'N' = kein Schritt

            steps += 1

        halted = (
            zustand in self.accept_states or zustand in self.reject_states
        )

        # Bandinhalt als String (ohne führende/abschließende Blanks)
        if band:
            min_pos = min(band.keys())
            max_pos = max(band.keys())
            tape_str = ''.join(
                band.get(i, self.BLANK) for i in range(min_pos, max_pos + 1)
            ).strip(self.BLANK) or self.BLANK
        else:
            tape_str = self.BLANK

        return halted, tape_str, steps

    def accepts(self, input_str: str, max_steps: int = 10_000) -> bool:
        """
        Gibt True zurück, wenn die TM den Eingabestring akzeptiert.

        @param input_str  Zu prüfender Eingabestring
        @param max_steps  Maximale Schrittzahl
        @return True wenn die TM in einem Akzeptierzustand hält
        """
        _, _, steps = self.simulate(input_str, max_steps)
        # Haltezustand nach der Simulation erneut prüfen
        band: Dict[int, str] = {i: c for i, c in enumerate(input_str)}
        kopf = 0
        zustand = self.initial_state
        s = 0

        while s < max_steps:
            if zustand in self.accept_states:
                return True
            if zustand in self.reject_states:
                return False
            symbol = band.get(kopf, self.BLANK)
            key = (zustand, symbol)
            if key not in self.transitions:
                return False
            neuer_zustand, schreib_symbol, richtung = self.transitions[key]
            band[kopf] = schreib_symbol
            zustand = neuer_zustand
            if richtung == 'R':
                kopf += 1
            elif richtung == 'L':
                kopf -= 1
            s += 1

        return False  # Zeitlimit überschritten

    def time_complexity_class(self) -> str:
        """
        Gibt eine heuristische Laufzeitkomplexitätsklasse zurück.
        Echte automatische Komplexitätsanalyse ist unentscheidbar (Rice's Theorem);
        diese Methode gibt eine konservative obere Schranke basierend auf der
        Anzahl der Zustände und der Übergangsregel-Struktur zurück.

        @return String wie "O(n)", "O(n^2)" usw.
        """
        n_states = len(self.states)
        n_trans = len(self.transitions)
        # Heuristik: Wenige Übergänge → lineare Zeit; viele → quadratisch
        if n_trans <= n_states:
            return "O(n)"
        elif n_trans <= n_states * n_states:
            return "O(n^2)"
        else:
            return "O(2^n)"


class NonDeterministicTM(TuringMachine):
    """
    Nichtdeterministische Turingmaschine (NTM).
    Erweiterung um zertifikat-basierte Verifikation (NP-Charakterisierung):
    Ein NTM akzeptiert, wenn *irgendein* nichtdeterministischer Pfad akzeptiert.
    Hier simulieren wir dies durch erschöpfende Suche über Zertifikate.
    """

    def __init__(
        self,
        states: List[str],
        alphabet: List[str],
        transitions: Dict[Tuple[str, str], List[Tuple[str, str, str]]],
        initial_state: str,
        accept_states: List[str],
        reject_states: Optional[List[str]] = None,
    ):
        """
        @param transitions  Übergangs*relation*: ein Schlüssel kann mehrere
                            mögliche Übergänge haben (Liste statt Einzelwert).
        """
        # Wir speichern die nichtdeterministische Relation separat
        self.nd_transitions: Dict[
            Tuple[str, str], List[Tuple[str, str, str]]
        ] = transitions
        # Für die Elternklasse: erster Übergang als deterministischer Fallback
        det_trans = {k: v[0] for k, v in transitions.items() if v}
        super().__init__(
            states, alphabet, det_trans, initial_state,
            accept_states, reject_states
        )

    def certificate_verify(self, input_str: str, certificate: str) -> bool:
        """
        Überprüft ein Zertifikat (Lösungszeuge) in Polynomialzeit.
        Das Zertifikat wird an die Eingabe angehängt (getrennt durch '#'),
        und die deterministische TM prüft, ob die Kombination akzeptiert wird.

        Diese Methode modelliert die NP-Definition:
            x ∈ L ⟺ ∃y : |y| ≤ p(|x|) ∧ V(x,y) = 1

        @param input_str   Eigentliche Eingabe x
        @param certificate Zeuge y
        @return True wenn das Zertifikat gültig ist
        """
        combined = input_str + '#' + certificate
        return self.accepts(combined)

    def nondeterministic_accepts(
        self, input_str: str, max_depth: int = 8
    ) -> bool:
        """
        Simuliert Nichtdeterminismus durch BFS über alle möglichen
        Berechnungspfade bis zur Tiefe max_depth.

        @param input_str  Eingabestring
        @param max_depth  Maximale Tiefe des Berechnungsbaums
        @return True wenn mindestens ein Pfad akzeptiert
        """
        # Zustand: (band, kopf_position, zustand)
        initial_band: Dict[int, str] = {i: c for i, c in enumerate(input_str)}
        frontier = [(initial_band, 0, self.initial_state)]

        for _ in range(max_depth):
            if not frontier:
                break
            next_frontier = []
            for band, kopf, zustand in frontier:
                if zustand in self.accept_states:
                    return True
                if zustand in self.reject_states:
                    continue
                symbol = band.get(kopf, self.BLANK)
                key = (zustand, symbol)
                moeglichkeiten = self.nd_transitions.get(key, [])
                for (nz, ns, richtung) in moeglichkeiten:
                    neues_band = dict(band)
                    neues_band[kopf] = ns
                    neuer_kopf = kopf + (1 if richtung == 'R' else -1 if richtung == 'L' else 0)
                    next_frontier.append((neues_band, neuer_kopf, nz))
            frontier = next_frontier

        # Prüfen ob irgendein Frontier-Zustand ein Akzeptierzustand ist
        return any(z in self.accept_states for _, _, z in frontier)


# ===========================================================================
# 2. Komplexitätsklassen
# ===========================================================================

class ComplexityClass:
    """
    Metadaten und Relationen zwischen bekannten Komplexitätsklassen.

    Die bekannten Inklusions-Ketten (alle bewiesen, sofern nicht anders vermerkt):
        L ⊆ NL = co-NL ⊆ P ⊆ NP ∩ co-NP ⊆ NP ⊆ PH ⊆ PSPACE = co-PSPACE ⊆ EXP
        P ⊆ BPP ⊆ BQP (vermutet: BQP ⊄ PH, aber nicht bewiesen)
        #P ist schwerer als NP (Toda 1991: PH ⊆ P^{#P})
    """

    # Beschreibungen der wichtigsten Komplexitätsklassen
    KNOWN_CLASSES: Dict[str, Dict[str, Any]] = {
        "P": {
            "description": (
                "Probleme, die von einer deterministischen TM in Polynomialzeit "
                "gelöst werden können. Formell: ∃k: L ∈ DTIME(n^k)."
            ),
            "examples": ["PRIMES (AKS 2002)", "LINEAR PROGRAMMING", "MATCHING"],
            "closed_under": ["complement", "union", "intersection", "concatenation"],
        },
        "NP": {
            "description": (
                "Probleme, die von einer nichtdeterministischen TM in Polynomialzeit "
                "gelöst werden. Äquivalent: Lösungen sind in Polynomialzeit verifizierbar."
            ),
            "examples": ["SAT", "CLIQUE", "HAMILTONIAN-CYCLE", "SUBSET-SUM", "3-COLORING"],
            "closed_under": ["union", "intersection"],
        },
        "coNP": {
            "description": (
                "Komplement von NP. Probleme, deren Verneinung in NP liegt. "
                "Beispiel: TAUTOLOGY (ist jede Belegung erfüllend?)."
            ),
            "examples": ["TAUTOLOGY", "GRAPH-NON-ISOMORPHISM (vermutet)", "UNSAT"],
            "closed_under": ["union", "intersection"],
        },
        "PSPACE": {
            "description": (
                "Probleme, die mit polynomiellem Speicher lösbar sind. "
                "PSPACE = co-PSPACE. Enthält PH und NP."
            ),
            "examples": ["TQBF (True Quantified Boolean Formula)", "GAME-OF-LIFE", "REGEX-EQUIVALENCE"],
            "closed_under": ["complement", "union", "intersection"],
        },
        "EXP": {
            "description": (
                "Probleme in DTIME(2^{p(n)}) für ein Polynom p. "
                "Enthält PSPACE. Enthält unentscheidbare Probleme NICHT."
            ),
            "examples": ["GENERALIZED-CHESS", "GENERALIZED-GO"],
            "closed_under": ["complement", "union", "intersection"],
        },
        "BPP": {
            "description": (
                "Bounded-error Probabilistic Polynomial time. "
                "Randomisierte Algorithmen mit Fehlerwahrscheinlichkeit ≤ 1/3. "
                "Vermutet: BPP = P (Derandomisierung)."
            ),
            "examples": ["PRIMALITY (vor AKS)", "POLYNOMIAL-IDENTITY-TESTING"],
            "closed_under": ["complement", "union", "intersection"],
        },
        "BQP": {
            "description": (
                "Bounded-error Quantum Polynomial time. "
                "Quantenalgorithmen mit Fehlerwahrscheinlichkeit ≤ 1/3. "
                "Enthält: Faktorisierung (Shor), Suche O(√N) (Grover)."
            ),
            "examples": ["INTEGER-FACTORING", "DISCRETE-LOG", "ELEMENT-DISTINCTNESS"],
            "closed_under": ["complement", "union", "intersection"],
        },
        "#P": {
            "description": (
                "Zählklasse: Anzahl der akzeptierenden NTM-Pfade. "
                "Toda 1991: PH ⊆ P^{#P}. #P-vollständig: #SAT, PERMANENT."
            ),
            "examples": ["#SAT", "PERMANENT", "#MATCHING"],
            "closed_under": ["union"],
        },
        "PP": {
            "description": (
                "Probabilistic Polynomial time (Mehrheitsentscheid). "
                "PP ⊇ NP ∪ co-NP. PP = PostBQP (Aaronson 2005)."
            ),
            "examples": ["MAJORITY-SAT"],
            "closed_under": ["complement", "union", "intersection"],
        },
        "IP": {
            "description": (
                "Interactive Proof System. Shamir 1992: IP = PSPACE. "
                "Verifier = polynomielle Zeit, unbeschränkte Runden."
            ),
            "examples": ["TQBF", "GRAPH-NON-ISOMORPHISM"],
            "closed_under": ["complement"],
        },
        "PH": {
            "description": (
                "Polynomielle Hierarchie: ⋃_k Σ_k^P. "
                "Vermutet: PH ist echt unendlich. Wenn P=NP dann PH = P."
            ),
            "examples": ["Σ_2^P: Π_2^P-Probleme", "Π_2-SAT"],
            "closed_under": ["union", "intersection"],
        },
    }

    # Bekannte Inklusionen (bewiesene Relationen)
    _RELATIONS: List[Tuple[str, str]] = [
        ("P", "NP"), ("P", "coNP"), ("P", "BPP"),
        ("NP", "PH"), ("coNP", "PH"),
        ("PH", "PSPACE"), ("PSPACE", "EXP"),
        ("NP", "PP"), ("coNP", "PP"),
        ("PP", "PSPACE"),
        ("BPP", "PSPACE"),
        ("BQP", "PSPACE"),
        ("NP", "#P"),   # #P ist mindestens so schwer wie NP
    ]

    def __init__(self, name: str):
        """
        @param name  Name der Komplexitätsklasse (z.B. "NP")
        """
        if name not in self.KNOWN_CLASSES:
            raise ValueError(
                f"Unbekannte Klasse '{name}'. Bekannte: {list(self.KNOWN_CLASSES)}"
            )
        self.name = name
        self._info = self.KNOWN_CLASSES[name]

    @property
    def description(self) -> str:
        """Gibt die Beschreibung der Klasse zurück."""
        return self._info["description"]

    @property
    def examples(self) -> List[str]:
        """Gibt Beispielprobleme zurück."""
        return self._info["examples"]

    def contains(self, problem: str) -> bool:
        """
        Prüft, ob ein (durch Namen bekanntes) Problem in dieser Klasse liegt.

        @param problem  Problemname (z.B. "SAT", "PRIMES")
        @return True wenn das Problem als Beispiel bekannt ist
        """
        return any(problem.upper() in ex.upper() for ex in self._info["examples"])

    def is_closed_under(self, operation: str) -> bool:
        """
        Gibt an, ob die Klasse unter einer Mengenoperation abgeschlossen ist.

        @param operation  "complement", "union" oder "intersection"
        @return True wenn abgeschlossen
        """
        return operation in self._info.get("closed_under", [])

    def known_relations(self) -> Dict[str, List[str]]:
        """
        Gibt bekannte Inklusionsbeziehungen zurück.

        @return Dict mit Schlüsseln "supersets" und "subsets"
        """
        supersets = [b for (a, b) in self._RELATIONS if a == self.name]
        subsets   = [a for (a, b) in self._RELATIONS if b == self.name]
        return {"supersets": supersets, "subsets": subsets}

    def __repr__(self) -> str:
        return f"ComplexityClass({self.name!r})"


class PolynomialHierarchy:
    """
    Polynomielle Hierarchie (PH).

    Definitionen:
        Σ_0^P = Π_0^P = Δ_0^P = P
        Σ_{k+1}^P = NP^{Σ_k^P}     (NP mit Σ_k^P-Orakel)
        Π_{k+1}^P = co-NP^{Σ_k^P}
        Δ_{k+1}^P = P^{Σ_k^P}

    Bekannte Fakten:
        - P = NP ⟹ PH kollabiert auf Σ_0^P = P
        - NP = co-NP ⟹ PH kollabiert auf Σ_1^P = NP
        - Für alle k: Σ_k^P ⊆ Σ_{k+1}^P (vermutet: echt)
    """

    def sigma(self, k: int) -> str:
        """
        Gibt die Beschreibung von Σ_k^P zurück.

        @param k  Stufe (k ≥ 0)
        @return String-Beschreibung
        """
        if k < 0:
            raise ValueError("k muss ≥ 0 sein.")
        if k == 0:
            return "Σ_0^P = P  (Polynomialzeit, keine Quantoren)"
        if k == 1:
            return "Σ_1^P = NP  (∃-Quantor, Zeuge polynomiell)"
        orakel = self.sigma(k - 1).split("=")[0].strip()
        return (
            f"Σ_{k}^P = NP^{{Σ_{k-1}^P}}  "
            f"(NP mit {orakel}-Orakel; entspricht ∃∀∃…-Quantorpräfix der Länge {k})"
        )

    def pi(self, k: int) -> str:
        """
        Gibt die Beschreibung von Π_k^P zurück.

        @param k  Stufe (k ≥ 0)
        @return String-Beschreibung
        """
        if k < 0:
            raise ValueError("k muss ≥ 0 sein.")
        if k == 0:
            return "Π_0^P = P  (Polynomialzeit, keine Quantoren)"
        if k == 1:
            return "Π_1^P = co-NP  (∀-Quantor)"
        return (
            f"Π_{k}^P = co-NP^{{Σ_{k-1}^P}}  "
            f"(Komplement von Σ_{k}^P; ∀∃∀…-Quantorpräfix der Länge {k})"
        )

    def delta(self, k: int) -> str:
        """
        Gibt die Beschreibung von Δ_k^P zurück.

        @param k  Stufe (k ≥ 1)
        @return String-Beschreibung
        """
        if k < 1:
            raise ValueError("k muss ≥ 1 sein.")
        if k == 1:
            return "Δ_1^P = P"
        return (
            f"Δ_{k}^P = P^{{Σ_{k-1}^P}}  "
            f"(Deterministisch polynomiell mit Σ_{k-1}^P-Orakel)"
        )

    def oracle_characterization(self, k: int) -> str:
        """
        Beschreibt die Orakel-Charakterisierung der k-ten Stufe.

        @param k  Stufe
        @return Beschreibung
        """
        return (
            f"Stufe {k} der Polynomiellen Hierarchie:\n"
            f"  {self.sigma(k)}\n"
            f"  {self.pi(k)}\n"
            f"  {self.delta(k) if k >= 1 else '(nicht definiert)'}\n"
            f"  Relation: Σ_{k}^P ∪ Π_{k}^P ⊆ Δ_{{k+1}}^P ⊆ Σ_{{k+1}}^P"
        )

    def collapse_consequence(self) -> str:
        """
        Gibt eine Beschreibung der Kollaps-Konsequenz zurück.

        @return String mit den wichtigsten Kollaps-Szenarien
        """
        return (
            "Kollaps der Polynomiellen Hierarchie:\n"
            "  1. P = NP  ⟹  PH = P = Σ_0^P  (vollständiger Kollaps)\n"
            "  2. NP = co-NP  ⟹  PH = NP = Σ_1^P  (Kollaps auf Stufe 1)\n"
            "  3. Σ_k^P = Π_k^P  ⟹  PH = Σ_k^P  (Kollaps auf Stufe k)\n"
            "  4. Σ_k^P = Σ_{k+1}^P  ⟹  PH = Σ_k^P  (ab Stufe k konstant)\n"
            "Folgerung: Wenn PH unendlich ist (vermutet), dann ist P ≠ NP."
        )


# ===========================================================================
# 3. NP-Vollständigkeit
# ===========================================================================

class NPCompletenessTheory:
    """
    Theorie der NP-Vollständigkeit:
        - Cook-Levin-Theorem (1971/1973): SAT ist NP-vollständig
        - Karp's 21 Probleme (1972): erste vollständige Liste
        - Polynomialzeit-Reduktionen
        - Selbstreduktion
    """

    # Karp's 21 NP-vollständige Probleme (historisch)
    _KARP_21: List[Dict[str, str]] = [
        {"name": "SAT",             "reduction_from": "Definition (Cook-Levin)", "year": 1971,
         "desc": "Gibt es eine Belegung, die eine boolesche Formel erfüllt?"},
        {"name": "3-SAT",           "reduction_from": "SAT",             "year": 1972,
         "desc": "SAT eingeschränkt auf Klauseln mit genau 3 Literalen."},
        {"name": "CLIQUE",          "reduction_from": "3-SAT",           "year": 1972,
         "desc": "Gibt es eine k-Clique im Graphen?"},
        {"name": "SET-PACKING",     "reduction_from": "CLIQUE",          "year": 1972,
         "desc": "k disjunkte Teilmengen aus einer Familie?"},
        {"name": "VERTEX-COVER",    "reduction_from": "CLIQUE",          "year": 1972,
         "desc": "k Knoten, die alle Kanten abdecken?"},
        {"name": "SET-COVER",       "reduction_from": "VERTEX-COVER",    "year": 1972,
         "desc": "k Mengen, die das Universum überdecken?"},
        {"name": "FEEDBACK-VERTEX-SET", "reduction_from": "VERTEX-COVER","year": 1972,
         "desc": "k Knoten entfernen, um azyklischen Graph zu erhalten?"},
        {"name": "FEEDBACK-ARC-SET","reduction_from": "FEEDBACK-VERTEX-SET","year": 1972,
         "desc": "k Kanten entfernen, um azyklischen gerichteten Graph zu erhalten?"},
        {"name": "DIRECTED-HAMILTONIAN-CYCLE","reduction_from": "3-SAT", "year": 1972,
         "desc": "Hamiltonkreis im gerichteten Graphen?"},
        {"name": "UNDIRECTED-HAMILTONIAN-CYCLE","reduction_from": "DIRECTED-HAMILTONIAN-CYCLE","year": 1972,
         "desc": "Hamiltonkreis im ungerichteten Graphen?"},
        {"name": "0-1-INTEGER-PROGRAMMING","reduction_from": "SAT",      "year": 1972,
         "desc": "0-1-ILP hat eine Lösung?"},
        {"name": "3-COLOR",         "reduction_from": "3-SAT",           "year": 1972,
         "desc": "Graph mit 3 Farben färbbar?"},
        {"name": "CLIQUE-COVER",    "reduction_from": "3-COLOR",         "year": 1972,
         "desc": "Graph ist Vereinigung von k Cliquen?"},
        {"name": "EXACT-COVER",     "reduction_from": "3-SAT",           "year": 1972,
         "desc": "Exakte Überdeckung des Universums durch k Mengen?"},
        {"name": "3D-MATCHING",     "reduction_from": "EXACT-COVER",     "year": 1972,
         "desc": "Perfektes 3-dimensionales Matching?"},
        {"name": "STEINER-TREE",    "reduction_from": "EXACT-COVER",     "year": 1972,
         "desc": "Steiner-Baum mit Gewicht ≤ k?"},
        {"name": "HITTING-SET",     "reduction_from": "SET-COVER",       "year": 1972,
         "desc": "Menge, die alle Teilmengen schneidet?"},
        {"name": "KNAPSACK",        "reduction_from": "EXACT-COVER",     "year": 1972,
         "desc": "Rucksackproblem: Wert ≥ k bei Gewicht ≤ W?"},
        {"name": "JOB-SEQUENCING", "reduction_from": "KNAPSACK",        "year": 1972,
         "desc": "Jobs einplanen mit Strafe ≤ k?"},
        {"name": "PARTITION",       "reduction_from": "KNAPSACK",        "year": 1972,
         "desc": "Menge in zwei Hälften gleicher Summe teilen?"},
        {"name": "MAX-CUT",         "reduction_from": "PARTITION",       "year": 1972,
         "desc": "Schnitt mit Gewicht ≥ k?"},
    ]

    _NP_COMPLETE_SET = {p["name"] for p in _KARP_21}

    @classmethod
    def cook_levin_theorem(cls) -> Dict[str, Any]:
        """
        Gibt historische und inhaltliche Informationen zum Cook-Levin-Theorem zurück.

        @return Dictionary mit Theorem-Details
        """
        return {
            "theorem": "SAT ist NP-vollständig.",
            "published_by_cook": 1971,
            "published_by_levin": 1973,
            "key_idea": (
                "Jede NTM-Berechnung in Polynomialzeit kann als boolesche Formel "
                "kodiert werden (Tableau-Konstruktion). Die Formel ist erfüllbar "
                "genau dann, wenn die NTM akzeptiert."
            ),
            "tableau_construction": (
                "Für Eingabe x, |x|=n, p(n) = Zeitschranke:\n"
                "  - Variablen: cell[i][j][s] = 'In Schritt i ist Position j im Zustand s'\n"
                "  - Klauseln kodieren: Initialisierung, Übergänge, Akzeptanz\n"
                "  - Größe: O(p(n)^2 · |Q| · |Γ|)"
            ),
            "consequence": (
                "Jedes NP-Problem kann in Polynomialzeit auf SAT reduziert werden. "
                "Daher ist SAT das 'schwerste' Problem in NP."
            ),
            "significance": (
                "Grundstein der Komplexitätstheorie. Zeigt: Wenn SAT ∈ P, dann P = NP."
            ),
        }

    @classmethod
    def karp_21_problems(cls) -> List[Dict[str, Any]]:
        """
        Gibt Karps 21 originale NP-vollständige Probleme zurück.

        @return Liste von Dictionaries mit Problemdetails
        """
        return list(cls._KARP_21)

    @classmethod
    def reduction(cls, problem_a: str, problem_b: str) -> Dict[str, str]:
        """
        Beschreibt die bekannte Polynomialzeitreduktion von problem_a nach problem_b.

        @param problem_a  Quellproblem (wird reduziert)
        @param problem_b  Zielproblem
        @return Dict mit Beschreibung der Reduktion
        """
        a = problem_a.upper()
        b = problem_b.upper()

        # Bekannte direkte Reduktionen
        known: Dict[Tuple[str, str], str] = {
            ("3-SAT", "CLIQUE"): (
                "Für jede Klausel (l1∨l2∨l3) füge 3 Knoten ein. "
                "Verbinde Knoten aus verschiedenen Klauseln, wenn sie konsistent sind (kein l und ¬l). "
                "Formel erfüllbar ⟺ k-Clique (k = Klauselanzahl)."
            ),
            ("CLIQUE", "VERTEX-COVER"): (
                "G hat k-Clique ⟺ Komplementgraph Ḡ hat (n-k)-Vertex-Cover. "
                "Konstruktion: Ḡ enthält genau die Kanten, die G nicht enthält."
            ),
            ("SAT", "3-SAT"): (
                "Klausel mit m > 3 Literalen wird durch Einführung neuer Hilfsvariablen "
                "in (m-2) Klauseln der Länge 3 transformiert (lineare Aufblähung)."
            ),
            ("3-SAT", "3-COLOR"): (
                "Gadget-Konstruktion: Wahrheitskonstanten T/F als Dreieck kodieren. "
                "Pro Variable xi: xi und ¬xi als Knoten mit Kante. "
                "Pro Klausel: 6-Knoten-Gadget, das nur mit gültigem Zeuge 3-färbbar ist."
            ),
            ("VERTEX-COVER", "SET-COVER"): (
                "Jede Kante e=(u,v) entspricht einem Element. "
                "Jeder Knoten v entspricht der Menge aller anliegenden Kanten. "
                "k-Vertex-Cover ⟺ k-Set-Cover."
            ),
        }

        result = known.get((a, b)) or known.get((b, a))
        if result:
            return {
                "from": problem_a,
                "to": problem_b,
                "direction": f"{problem_a} ≤_p {problem_b}",
                "description": result,
                "complexity": "Polynomialzeit",
            }
        else:
            return {
                "from": problem_a,
                "to": problem_b,
                "direction": f"{problem_a} ≤_p {problem_b}",
                "description": (
                    "Diese spezifische Reduktion ist nicht in der Datenbank. "
                    "Verwende Cook-Levin als generellen Weg: Beide Probleme sind NP-vollständig, "
                    "daher existiert eine Polynomialzeit-Reduktion (nicht explizit gespeichert)."
                ),
                "complexity": "Polynomialzeit (existiert, nicht explizit)",
            }

    @classmethod
    def is_np_complete(cls, problem_name: str) -> bool:
        """
        Gibt an, ob ein Problem als NP-vollständig bekannt ist (aus Karp's Liste).

        @param problem_name  Problemname
        @return True wenn NP-vollständig
        """
        return problem_name.upper() in cls._NP_COMPLETE_SET

    @classmethod
    def self_reducibility(cls, problem: str) -> bool:
        """
        Gibt an, ob ein Problem self-reducible ist.
        NP-vollständige Probleme sind typischerweise self-reducible:
        Die Optimierungsvariante reduziert sich auf die Entscheidungsvariante.

        @param problem  Problemname
        @return True für bekannte NP-vollständige Probleme
        """
        # Alle NP-vollständigen Probleme in der Karp-Liste sind self-reducible
        return cls.is_np_complete(problem)


# ===========================================================================
# 4. SAT-Löser
# ===========================================================================

class SATSolver:
    """
    Implementierungen verschiedener SAT-Algorithmen:
        - DPLL (Davis-Putnam-Logemann-Loveland, 1962)
        - CDCL-Schritt (Conflict-Driven Clause Learning)
        - Resolution-Widerlegung
        - Tseitin-Transformation (beliebige Formel → CNF)
    """

    @staticmethod
    def dpll(
        clauses: List[List[int]], variables: List[int]
    ) -> Optional[Dict[int, bool]]:
        """
        DPLL-Algorithmus für SAT.

        Eingabe: Klauselmenge in CNF (Konjunktive Normalform).
        Literale: positive Ganzzahl = Variable, negativ = negierte Variable.
        Beispiel: [[1, -2], [2, 3]] = (x1 ∨ ¬x2) ∧ (x2 ∨ x3)

        @param clauses    Liste von Klauseln (jede Klausel = Liste von Literalen)
        @param variables  Liste der Variablen-IDs
        @return Erfüllende Belegung als {var: bool} oder None wenn UNSAT
        """

        def solve(
            cls: List[List[int]], assignment: Dict[int, bool], remaining: List[int]
        ) -> Optional[Dict[int, bool]]:
            # Unit Propagation
            cls, assignment, conflict = _unit_propagate(cls, assignment)
            if conflict:
                return None  # Konflikt → Backtrack

            # Alle Klauseln erfüllt?
            active = [
                c for c in cls if not _is_tautology(c, assignment)
            ]
            if not active:
                return assignment  # SAT gefunden

            # Leere Klausel = Konflikt
            for c in active:
                free = [l for l in c if abs(l) not in assignment]
                if not free:
                    return None

            # Pure-Literal-Elimination
            all_lits: List[int] = [l for c in active for l in c]
            for lit in set(all_lits):
                if -lit not in all_lits and abs(lit) not in assignment:
                    assignment[abs(lit)] = lit > 0
                    all_lits = [l for c in active for l in c if abs(l) not in assignment]

            # Nächste Variable auswählen (First Unset)
            unset = [v for v in remaining if v not in assignment]
            if not unset:
                # Alles belegt — prüfen ob alle Klauseln erfüllt
                if all(_is_tautology(c, assignment) for c in cls):
                    return assignment
                return None

            var = unset[0]
            rest = unset[1:]

            # Branching: erst True, dann False
            for val in [True, False]:
                new_asgn = dict(assignment)
                new_asgn[var] = val
                result = solve(cls, new_asgn, rest)
                if result is not None:
                    return result

            return None  # UNSAT

        return solve(clauses, {}, list(variables))

    @staticmethod
    def cdcl_step(
        clauses: List[List[int]], assignment: Dict[int, bool]
    ) -> str:
        """
        Führt einen CDCL-Schritt (Conflict-Driven Clause Learning) aus.
        Gibt eine Beschreibung der durchgeführten Aktion zurück.

        CDCL ist der Kern moderner SAT-Löser (MiniSAT, zChaff, Glucose).

        @param clauses    Aktuelle Klauselmenge in CNF
        @param assignment Aktuelle Variablenbelegung
        @return Beschreibung des Schritts ("Unit Prop", "Konflikt", "Entscheidung", "SAT")
        """
        # Schritt 1: Unit Propagation
        new_cls, new_asgn, conflict = _unit_propagate(clauses, assignment)

        if conflict:
            # Schritt 2: Konfliktanalyse (vereinfacht)
            return (
                "KONFLIKT erkannt (leere Klausel). "
                "CDCL würde jetzt:\n"
                "  1. Konfliktklausel analysieren (Konfliktgraph)\n"
                "  2. Neue Lernklausel via 1-UIP (Unique Implication Point) schneiden\n"
                "  3. Nicht-chronologisches Backtracking zum Entscheidungslevel des zweiten Literals\n"
                "  4. Lernklausel zur Klauseldatenbank hinzufügen"
            )

        # Schritt 3: Alle Klauseln erfüllt?
        active = [c for c in new_cls if not _is_tautology(c, new_asgn)]
        if not active:
            return "SAT: Alle Klauseln erfüllt. Belegung ist eine Lösung."

        if new_asgn != assignment:
            propagated = set(new_asgn.keys()) - set(assignment.keys())
            return (
                f"UNIT PROPAGATION: Variable(n) {propagated} wurden erzwungen. "
                f"Neue Belegung hat {len(new_asgn)} Variablen gesetzt."
            )

        # Schritt 4: Entscheidung (VSIDS-Heuristik vereinfacht)
        all_vars = {abs(l) for c in active for l in c}
        unset_vars = all_vars - set(new_asgn.keys())
        if unset_vars:
            chosen = min(unset_vars)  # Einfachste Heuristik: kleinste Variable
            return (
                f"ENTSCHEIDUNG: Variable {chosen} = True gesetzt "
                f"(VSIDS-Heuristik würde häufigste Konflikt-Variable wählen). "
                f"Entscheidungslevel erhöht."
            )

        return "UNBEKANNTER ZUSTAND (kein Fortschritt möglich)."

    @staticmethod
    def resolution_refutation(
        clauses: List[List[int]]
    ) -> Optional[List[Tuple[List[int], List[int], List[int]]]]:
        """
        Versucht, durch Resolutionsschritte die leere Klausel abzuleiten
        (= UNSAT-Beweis).

        Resolution: Aus (A ∨ x) und (B ∨ ¬x) folgt (A ∨ B).

        @param clauses  Klauselmenge in CNF
        @return Liste von (Klausel1, Klausel2, Resolvent)-Tupeln oder None wenn SAT
        """
        clause_set = [frozenset(c) for c in clauses]
        derivations: List[Tuple[List[int], List[int], List[int]]] = []
        new: List[frozenset] = []

        for _ in range(1000):  # Iterationslimit
            pairs = list(itertools.combinations(clause_set, 2))
            added = False
            for c1, c2 in pairs:
                # Resolutionsvariable suchen
                for lit in c1:
                    if -lit in c2:
                        # Resolvent berechnen
                        resolvent = (c1 - {lit}) | (c2 - {-lit})
                        resolvent_f = frozenset(resolvent)

                        # Tautologie prüfen
                        is_taut = any(-l in resolvent_f for l in resolvent_f)
                        if is_taut:
                            continue

                        if resolvent_f not in clause_set:
                            derivations.append(
                                (sorted(c1), sorted(c2), sorted(resolvent_f))
                            )
                            if not resolvent_f:
                                return derivations  # Leere Klausel = UNSAT bewiesen
                            clause_set.append(resolvent_f)
                            new.append(resolvent_f)
                            added = True
                        break

            if not added:
                return None  # Keine neuen Resolventen → SAT (keine Widerlegung)

        return None  # Limit erreicht

    @staticmethod
    def tseitin_transform(formula_str: str) -> List[List[int]]:
        """
        Tseitin-Transformation: Wandelt eine aussagenlogische Formel in CNF um.
        Die Transformation ist linear in der Formellänge (im Gegensatz zur
        exponentiellen naiven CNF-Umformung).

        Unterstützte Formate:
            - Literale: "1", "-2", "3"
            - AND: "AND(a,b)"
            - OR: "OR(a,b)"
            - NOT: "NOT(a)"

        Für einfache Variablen-IDs (Integerliterale) wird die Formel direkt
        als CNF interpretiert.

        @param formula_str  Formelstring, z.B. "AND(OR(1,-2),OR(-1,3))"
        @return Liste von Klauseln (CNF)
        """
        # Vereinfachte Implementierung: Parsen als Liste von Literal-Listen
        # Echte Tseitin-Transformation für komplexe Formeln
        formula_str = formula_str.strip()

        # Fall 1: Direkte Klauselform "[[1,-2],[3,4]]"
        if formula_str.startswith("[[") or formula_str.startswith("[("):
            try:
                import ast
                result = ast.literal_eval(formula_str)
                return [list(c) for c in result]
            except Exception:
                pass

        # Fall 2: Einfaches Literal
        try:
            lit = int(formula_str)
            return [[lit]]  # Einzel-Literal als Einzel-Klausel
        except ValueError:
            pass

        # Fall 3: AND(sub1, sub2) — Konjunktion
        if formula_str.upper().startswith("AND("):
            inner = formula_str[4:-1]
            # Einfaches Splitting (kein Nesting unterstützt in dieser Kurzform)
            parts = inner.split(",", 1)
            if len(parts) == 2:
                cls1 = SATSolver.tseitin_transform(parts[0].strip())
                cls2 = SATSolver.tseitin_transform(parts[1].strip())
                return cls1 + cls2

        # Fall 4: OR(sub1, sub2) — Disjunktion
        if formula_str.upper().startswith("OR("):
            inner = formula_str[3:-1]
            parts = inner.split(",", 1)
            if len(parts) == 2:
                try:
                    l1 = int(parts[0].strip())
                    l2 = int(parts[1].strip())
                    return [[l1, l2]]
                except ValueError:
                    pass

        # Fall 5: NOT(lit)
        if formula_str.upper().startswith("NOT("):
            inner = formula_str[4:-1].strip()
            try:
                lit = int(inner)
                return [[-lit]]
            except ValueError:
                pass

        # Fallback: Gesamtformel als einzelne Klausel
        return [[1]]  # Triviale Fallback-Klausel


# ===========================================================================
# 5. Schaltkreiskomplexität
# ===========================================================================

class CircuitComplexity:
    """
    Schaltkreiskomplexität: Größe und Tiefe boolescher Schaltkreise.

    Bekannte Resultate:
        - AC^0 ⊊ TC^0 ⊆ NC^1 ⊆ L ⊆ NL ⊆ NC^2 ⊆ NC ⊆ P
        - PARITY ∉ AC^0 (Furst-Saxe-Sipser 1984, Håstad 1987)
        - CLIQUE benötigt exponentielle Monotone-Schaltkreisgröße (Razborov 1985)
        - Natural Proofs Barriere (Razborov-Rudich 1994)
    """

    @staticmethod
    def circuit_size(boolean_function_tt: List[int]) -> int:
        """
        Schätzt die Schaltkreisgröße für eine boolesche Funktion, gegeben durch
        ihre Wahrheitstabelle.

        Methode: Obere Schranke via DNF/CNF-Konstruktion.
        Exakte Minimierung ist NP-schwer (Schaltkreis-Minimierungsproblem).

        @param boolean_function_tt  Wahrheitstabelle als Liste von 0/1 (Länge muss 2^n sein)
        @return Obere Schranke auf die Schaltkreisgröße (Gatter-Anzahl)
        """
        n_rows = len(boolean_function_tt)
        if n_rows == 0:
            return 0

        # n = Anzahl der Eingabevariablen
        n = int(math.log2(n_rows))
        if 2 ** n != n_rows:
            raise ValueError(
                f"Länge der Wahrheitstabelle ({n_rows}) muss eine 2er-Potenz sein."
            )

        # Anzahl der 1-Einträge (Minterme für DNF)
        ones = sum(boolean_function_tt)
        zeros = n_rows - ones

        # Obere Schranke: Minimum von DNF-Größe und CNF-Größe
        # DNF: ones Terme à n Literale + (ones-1) AND-Gatter + OR-Gatter
        # CNF: zeros Terme à n Literale
        if ones == 0:
            return 1  # Konstante 0: ein Gatter
        if zeros == 0:
            return 1  # Konstante 1: ein Gatter

        dnf_gates = ones * n + (ones - 1) + 1  # Literale + ANDs + OR
        cnf_gates = zeros * n + (zeros - 1) + 1
        return min(dnf_gates, cnf_gates)

    @staticmethod
    def circuit_depth(boolean_function_tt: List[int]) -> int:
        """
        Schätzt die Schaltkreistiefe für eine boolesche Funktion.

        Obere Schranke: O(log n) für NC^1-Formeln, O(log^2 n) für NC^2.
        Triviale obere Schranke: Tiefe 2 für DNF/CNF (2-Level-Schaltkreis).

        @param boolean_function_tt  Wahrheitstabelle
        @return Tiefe des Schaltkreises (2-Level-DNF-Schranke)
        """
        n_rows = len(boolean_function_tt)
        if n_rows == 0:
            return 0
        n = int(math.log2(n_rows))
        if 2 ** n != n_rows:
            raise ValueError("Länge muss 2er-Potenz sein.")

        ones = sum(boolean_function_tt)
        if ones == 0 or ones == n_rows:
            return 1  # Konstante Funktion: Tiefe 1

        # Balanciertes Baumtiefe für NC^1-Simulation: O(log n)
        # 2-Level DNF/CNF: Tiefe 2 (praktische obere Schranke)
        return 2

    @staticmethod
    def nc_class_check(circuit_description: Dict[str, Any]) -> str:
        """
        Bestimmt die NC-Klasse eines Schaltkreises anhand seiner Parameter.

        NC^k: Polynomielle Größe, O(log^k n) Tiefe, Fan-in 2.

        @param circuit_description  Dict mit 'depth_exponent', 'size', 'n', 'fan_in'
        @return Klassenname "NC^1", "NC^2", "AC^0" etc.
        """
        depth_exp = circuit_description.get("depth_exponent", 2)
        fan_in = circuit_description.get("fan_in", 2)
        unbounded = circuit_description.get("unbounded_fan_in", False)

        if unbounded and depth_exp == 0:
            return "AC^0 (konstante Tiefe, unbeschränkter Fan-in)"
        if unbounded and depth_exp >= 1:
            return f"AC^{depth_exp} (Tiefe O(log^{depth_exp} n), unbeschränkter Fan-in)"
        if fan_in == 2:
            if depth_exp == 1:
                return "NC^1 (Tiefe O(log n), Fan-in 2)"
            elif depth_exp == 2:
                return "NC^2 (Tiefe O(log^2 n), Fan-in 2)"
            else:
                return f"NC^{depth_exp} (Tiefe O(log^{depth_exp} n), Fan-in 2)"
        return "Nicht klassifizierbar mit den gegebenen Parametern."

    @staticmethod
    def razborov_rudich_natural_proof() -> Dict[str, Any]:
        """
        Beschreibt die Natural-Proofs-Barriere (Razborov-Rudich 1994).

        Ein Beweis ist "natürlich" wenn:
        1. Konstruktiv: Gegeben Schaltkreis der Größe S(n), in Zeit poly(2^n) entscheidbar
        2. Nützlich: Trennt P/poly von NP
        3. Largeness: Die Eigenschaft gilt für einen positiven Anteil aller Funktionen

        Wenn sichere Pseudozufallsgeneratoren existieren (Standardannahme),
        können natürliche Beweise keine unteren Schranken beweisen.

        @return Dict mit Barriere-Details
        """
        return {
            "barrier": "Natural Proofs Barriere",
            "authors": "Alexander Razborov & Steven Rudich",
            "year": 1994,
            "formal_definition": {
                "constructive": (
                    "Eine kombinatorische Eigenschaft C ist konstruktiv, wenn "
                    "gegeben ein Schaltkreis C_n: {0,1}^n→{0,1} der Größe S(n), "
                    "in Zeit 2^{O(n)} entschieden werden kann ob C_n ∈ C."
                ),
                "useful": (
                    "C ist nützlich für eine untere Schranke gegen P/poly, wenn: "
                    "Alle Funktionen in C haben super-polynomielle Schaltkreiskomplexität."
                ),
                "large": (
                    "C ist groß, wenn |C_n| / 2^{2^n} ≥ 1/2^{n^{O(1)}} "
                    "(positiver Bruchteil aller n-Bit-Funktionen)."
                ),
            },
            "consequence": (
                "Wenn natürliche Beweise P≠NP beweisen könnten und sichere PRGs existieren, "
                "würde daraus folgen, dass diese PRGs von polynomiellen Schaltkreisen gebrochen "
                "werden können — Widerspruch zur Sicherheitsannahme."
            ),
            "examples_of_natural_proofs": [
                "Razborovs Monotone-Schaltkreis-Beweis (Approximationsmethode)",
                "Håstads Switching-Lemma (AC^0 untere Schranken)",
            ],
            "what_is_not_ruled_out": [
                "Algebraische Geometrie-Ansätze (Mulmuley-Sohoni GCT)",
                "Algebrization-resistente Beweise",
                "Beweise ohne Largeness-Eigenschaft",
            ],
        }

    @staticmethod
    def monotone_lower_bound(n: int) -> int:
        """
        Razborovs 1985 Ergebnis: Die monotone Schaltkreisgröße für das
        k-CLIQUE-Problem ist Ω(n^{k/4}) für k ≤ n^{1/4}.

        Dies zeigt: Monotone Schaltkreise sind exponentiell schwächer als
        allgemeine Schaltkreise (da CLIQUE in P liegt).

        @param n  Graphgröße (Anzahl der Knoten)
        @return Untere Schranke (als Ganzzahl approximiert)
        """
        # Für k = n^{1/4} ergibt sich Ω(2^{n^{1/4}})
        # Wir verwenden k = max(3, floor(n^(1/4)))
        k = max(3, int(n ** 0.25))
        # Razborovs untere Schranke: n^{k/4}
        lower = int(n ** (k / 4))
        return lower


# ===========================================================================
# 6. Orakel-Komplexität
# ===========================================================================

class OracleComplexity:
    """
    Orakel-Trennungen und die Relativierungsbarriere für P vs. NP.

    Baker-Gill-Solovay 1975:
        - ∃ Orakel A: P^A = NP^A  (Trennungs-Beweis nicht möglich mit Diagonalisierung allein)
        - ∃ Orakel B: P^B ≠ NP^B  (Nachweis unmöglich über Relativierung)

    Konsequenz: Kein Beweis für P≠NP, der sich relativiert, kann funktionieren.
    """

    @staticmethod
    def baker_gill_solovay() -> Dict[str, Any]:
        """
        Beschreibt das Baker-Gill-Solovay-Theorem (1975).

        @return Dict mit Theorem-Details und Beweisideen
        """
        return {
            "theorem": "Baker-Gill-Solovay 1975",
            "statement": (
                "Es existieren Orakel A und B, so dass:\n"
                "  P^A = NP^A  und  P^B ≠ NP^B"
            ),
            "oracle_A": {
                "description": "A = PSPACE-vollständiges Problem (z.B. TQBF)",
                "consequence": "P^A = NP^A = PSPACE^A = PSPACE",
                "intuition": (
                    "Mit einem PSPACE-Orakel können sowohl deterministische als auch "
                    "nichtdeterministische Maschinen PSPACE lösen. "
                    "Daher kollabiert die Hierarchie."
                ),
            },
            "oracle_B": {
                "description": "B = zufällig gewähltes Orakel (mit W-Maß 1)",
                "consequence": "P^B ≠ NP^B fast sicher",
                "intuition": (
                    "NP^B enthält das Problem 'Gibt es x ∈ B mit |x|=n?', "
                    "welches O(2^n) Abfragen erfordert, aber von keiner det. polynomiellen "
                    "Maschine mit polynomiell vielen B-Abfragen gelöst werden kann."
                ),
            },
            "barrier_consequence": (
                "Jede Beweistechnik für P≠NP, die sich relativiert "
                "(d.h. für beliebige Orakel funktioniert), scheitert notwendigerweise: "
                "Denn es gibt Orakel, für die P=NP gilt."
            ),
        }

    @staticmethod
    def relativization_barrier() -> str:
        """
        Erklärt die Relativierungsbarriere als Beweishindernis.

        @return Beschreibungstext
        """
        return (
            "Relativierungsbarriere (Baker-Gill-Solovay 1975):\n\n"
            "Definition: Ein Beweis 'relativiert', wenn er für beliebige Orakel O "
            "die analoge Aussage P^O = NP^O (bzw. ≠) beweist.\n\n"
            "Problem: Da sowohl Orakel A mit P^A = NP^A als auch B mit P^B ≠ NP^B existieren,\n"
            "kann kein relativierender Beweis die Frage P = NP? entscheiden.\n\n"
            "Beispiele für Beweistechniken, die relativieren (und damit scheitern):\n"
            "  - Diagonalisierungsargumente (Cantor, Turing)\n"
            "  - Zeitkonstruierbarkeitsargumente\n"
            "  - Kleene-Rekursions-basierte Argumente\n\n"
            "Nicht-relativierende Techniken (potenzielle Wege):\n"
            "  - Algebrization (schwache Form, reicht auch nicht)\n"
            "  - Geometric Complexity Theory (Mulmuley)\n"
            "  - Schaltkreis-untere-Schranken (bisher erfolgreich für eingeschränkte Modelle)"
        )

    @staticmethod
    def algebraic_oracle(problem: str) -> str:
        """
        Beschreibt algebraische Orakel-Resultate (LFKN, Shamir IP=PSPACE).

        @param problem  Problem oder Protokoll (z.B. "IP=PSPACE", "MIP*=RE")
        @return Beschreibung
        """
        known = {
            "IP=PSPACE": (
                "Shamir 1992: IP = PSPACE.\n"
                "Beweis über arithmetische Protokolle (LFKN 1990 für #SAT ∈ IP):\n"
                "  - Arithmetisierung: Boolesche Formeln → Polynome über endlichem Körper\n"
                "  - Sumcheck-Protokoll: Prüfer verifiziert Polynomauswertungen\n"
                "  - Beliebige PSPACE-Probleme haben interaktive Beweise\n"
                "Nicht-relativierend: ∃ Orakel O mit IP^O ≠ PSPACE^O."
            ),
            "MIP*=RE": (
                "Ji-Natarajan-Vidick-Wright-Yuen 2020: MIP* = RE.\n"
                "MIP* = Multi-Prover Interactive Proofs mit verschränkten Quantenproofs.\n"
                "RE = rekursiv aufzählbare Sprachen (enthält unentscheidbare Probleme!).\n"
                "Konsequenz: Connes' Einbettungsvermutung ist falsch."
            ),
            "PCP=NP": (
                "Arora-Lund-Motwani-Sudan-Szegedy 1992: PCP-Theorem.\n"
                "NP = PCP(O(log n), O(1)): Jeder NP-Beweis kann mit O(log n) Zufallsbits\n"
                "und O(1) Lesebits auf Korrektheit geprüft werden.\n"
                "Implikation: APX-Vollständigkeit von MAX-SAT (Inapproximierbarkeit)."
            ),
        }
        key = problem.upper().replace(" ", "").replace("_", "")
        for k, v in known.items():
            if k.replace("=", "").replace("*", "") in key.replace("=", "").replace("*", ""):
                return v
        return (
            f"Für '{problem}' ist kein algebraisches Orakel-Resultat gespeichert.\n"
            "Bekannte Resultate: IP=PSPACE, MIP*=RE, PCP=NP."
        )


# ===========================================================================
# 7. Beweiskomplexität
# ===========================================================================

class ProofComplexity:
    """
    Beweiskomplexität: Länge/Größe von Beweisen in formalen Systemen.

    Bekannte Resultate:
        - Resolution: Exponentiell für Pigeonhole-Prinzip (Haken 1985)
        - Frege-Systeme: Untere Schranken unbekannt
        - Natural Proofs Barriere: Razborov-Rudich 1994
        - Algebrization Barriere: Aaronson-Wigderson 2009
        - GCT: Mulmuley-Sohoni (ohne bekannte vollständige Barriere)
    """

    @staticmethod
    def resolution_width_lower_bound(k: int, n: int) -> int:
        """
        Ben-Sasson und Wigderson 1999:
        Die Resolution-Breite für die k-Clique-Formel über n Knoten ist Ω(k²).

        Resolution-Breite W(φ): Maximale Klauselgröße in einem minimalen Beweis.

        @param k  Clique-Größe
        @param n  Graphgröße (Knoten)
        @return Untere Schranke auf die Resolution-Breite
        """
        # Ben-Sasson-Wigderson: W(CLIQUE_n,k) ≥ k^2 / 2
        return (k * k) // 2

    @staticmethod
    def frege_lower_bound_status() -> str:
        """
        Gibt den aktuellen Stand zu Frege-System-Schranken zurück.

        @return Statusbeschreibung
        """
        return (
            "Frege-Systeme (auch: propositionale Beweissysteme):\n\n"
            "Status 2026: Keine superpolynomiellen unteren Schranken bekannt.\n\n"
            "Frege-System F:\n"
            "  - Arbeitet mit aussagenlogischen Formeln (beliebige Tiefe)\n"
            "  - Axiome: Tautologien des aussagenlogischen Kalküls\n"
            "  - Regeln: Modus Ponens, Substitution\n\n"
            "Extended Frege (EF = Frege + INTRO-Variablen):\n"
            "  - Polynomiell äquivalent zu allen 'natürlichen' Beweissystemen\n"
            "  - Optimal unter Annahme NP ⊄ coNP/poly\n\n"
            "Warum so schwierig:\n"
            "  1. Frege-Beweise können Tautologien kurz darstellen (z.B. PHP in Frege poly)\n"
            "  2. Bekannte Techniken (Haken-Stil für Resolution) greifen nicht\n"
            "  3. Fehlendes kombinatorisches Verständnis der Frege-Beweisstruktur\n\n"
            "Verbindung zu P vs NP:\n"
            "  Wenn NP ⊄ coNP/poly, dann hat kein Beweissystem polynomielle Beweislängen\n"
            "  für alle Tautologien. Extended Frege wäre dann nicht optimal."
        )

    @staticmethod
    def natural_proofs_barrier() -> Dict[str, Any]:
        """
        Beschreibt die Natural-Proofs-Barriere ausführlich.

        @return Dict mit Barriere-Details
        """
        return CircuitComplexity.razborov_rudich_natural_proof()

    @staticmethod
    def algebrization_barrier() -> Dict[str, Any]:
        """
        Beschreibt die Algebrization-Barriere (Aaronson-Wigderson 2009).

        Algebrization ist stärker als Relativierung:
        Man erlaubt dem Orakel, als Polynom über einem endlichen Körper
        ausgewertet zu werden.

        @return Dict mit Barriere-Details
        """
        return {
            "barrier": "Algebrization",
            "authors": "Scott Aaronson & Avi Wigderson",
            "year": 2009,
            "definition": (
                "Ein Beweis 'algebriert', wenn er für beliebige Orakel A und deren\n"
                "algebraische Erweiterung Ã (als Polynom über endlichem Körper) gilt."
            ),
            "strength": (
                "Stärker als Relativierung: Erfasst auch algebraische Techniken wie\n"
                "  - IP=PSPACE-Beweis (Shamir)\n"
                "  - Sumcheck-Protokoll\n"
                "  - LFKN-Protokoll für #SAT\n"
                "ABER: GCT (Mulmuley) algebriert möglicherweise nicht."
            ),
            "ruled_out": [
                "Alle diagonalisierenden Beweise",
                "Relativierende Beweise",
                "Algebraische Protokoll-basierte Beweise (IP=PSPACE-Stil)",
            ],
            "not_ruled_out": [
                "Geometric Complexity Theory (Mulmuley-Sohoni)",
                "Schaltkreis-Methoden ohne Polynominterpolation",
                "Kombinatorische Geometrie-Argumente",
            ],
            "separation": (
                "∃ Orakel A: NP^A ⊄ P/poly^{Ã}  und  ∃ Orakel B: P^B = NP^B\n"
                "Daher: Kein algebrierender Beweis kann P≠NP zeigen."
            ),
        }

    @staticmethod
    def geometric_complexity_theory() -> Dict[str, Any]:
        """
        Beschreibt den Geometric Complexity Theory (GCT)-Ansatz.

        GCT von Mulmuley & Sohoni (2001–) versucht, Schranken für algebraische
        Komplexitätsmaße über algebraische Geometrie und Darstellungstheorie zu beweisen.

        @return Dict mit GCT-Übersicht
        """
        return {
            "approach": "Geometric Complexity Theory (GCT)",
            "authors": "Ketan Mulmuley & Milind Sohoni",
            "started": 2001,
            "main_conjecture": (
                "Permanente eines n×n-Matrixpolynoms erfordert super-polynomiellen\n"
                "algebraischen Schaltkreis im Vergleich zur Determinante. "
                "(Valiant's VNP ≠ VP Hypothese, algebraisches Analogon zu P≠NP)"
            ),
            "key_idea": (
                "Komplexitätslücken zwischen Permanent und Determinante entsprechen\n"
                "Lücken zwischen Darstellungen von GL_n-Gruppenorbits in algebraischer Geometrie.\n"
                "Konkreter: Wenn perm_n nicht als det_m für m=poly(n) dargestellt werden kann,\n"
                "dann muss es 'fehlende' symmetrische Funktionen (Plethysmen) geben."
            ),
            "progress": [
                "2001: Grundrahmen und Verbindung zu Darstellungstheorie",
                "2008: GCT I+II — Formalisierung der Obstruktionen",
                "2012-2020: Partiell bewiesene Spezialfälle (beschränkte Tiefe)",
                "Aktuell: Berechnung algebraisch-geometrischer Invarianten (sehr schwer)",
            ],
            "barriers_overcome": (
                "GCT ist möglicherweise nicht-relativierend und nicht-algebrisierend,\n"
                "da es fundamentell andere mathematische Strukturen nutzt."
            ),
            "criticism": [
                "Sehr abstrakt und schwer zugänglich",
                "Keine konkreten Komplexitätsschranken bisher bewiesen",
                "Fehlende Brücke zwischen Darstellungstheorie und Komplexität",
            ],
        }


# ===========================================================================
# 8. Quantenkomplexität
# ===========================================================================

class QuantumComplexity:
    """
    Quantenkomplexitätstheorie: BQP, Quantenalgorithmen, Quantenbeweissysteme.

    Bekannte Fakten:
        - P ⊆ BQP ⊆ PSPACE
        - BQP ⊆ PP (Bennett-Bernstein-Brassard-Vazirani 1997)
        - INTEGER-FACTORING ∈ BQP (Shor 1994)
        - Grover 1996: Unstrukturierte Suche in O(√N) Quantenoperationen
        - Raz-Tal 2019: BQP ⊄ PH relativ zu einem Orakel (Orakel-Separation!)
    """

    @staticmethod
    def bqp_vs_ph() -> str:
        """
        Beschreibt den Raz-Tal 2019 Orakel-Separation BQP vs. PH.

        @return Beschreibungstext
        """
        return (
            "Raz-Tal 2019: Orakel-Separation BQP ⊄ PH\n\n"
            "Theorem (Raz-Tal): ∃ Orakel A mit L_A ∈ BQP^A \\ PH^A.\n\n"
            "Verwendetes Problem: 'Forrelation' (Aaronson 2010)\n"
            "  - Eingabe: Zwei Tabellen f,g: {0,1}^n → {±1}\n"
            "  - Aufgabe: Ist |Σ_{x,y} f(x)·(-1)^{x·y}·g(y)| / 2^n ≥ 1/100?\n"
            "  - BQP: 1 Quantenabfrage genügt (Hadamard-Transformation)\n"
            "  - PH: Benötigt Ω̃(n^{1/4}) Abfragen (Razborov-Rudich-Typ-Schranke)\n\n"
            "Bedeutung:\n"
            "  - Erste Orakel-Separation zwischen BQP und der gesamten PH\n"
            "  - Zeigt: Quantencomputer haben strukturelle Vorteile über die PH\n"
            "  - ABER: Relativierendes Ergebnis → sagt nichts über reale Separation aus\n\n"
            "Vorher bekannt:\n"
            "  - Bernstein-Vazirani 1993: BQP ⊄ BPP relativ zu Orakel\n"
            "  - Simon 1994: Subexponentieller Quantenvorteil"
        )

    @staticmethod
    def grovers_speedup() -> Dict[str, Any]:
        """
        Beschreibt Grovers Quantensuchalgorithmus (1996).

        @return Dict mit Algorithmus-Details
        """
        return {
            "algorithm": "Grover's Suchalgorithmus",
            "author": "Lov Grover",
            "year": 1996,
            "problem": (
                "Gegeben: Unstrukturierte Datenbank mit N Einträgen, "
                "genau 1 davon erfüllt Eigenschaft f(x)=1.\n"
                "Aufgabe: Finde diesen Eintrag."
            ),
            "classical_complexity": "Θ(N) Abfragen (deterministisch und randomisiert)",
            "quantum_complexity": "O(√N) Abfragen (optimal, BBBV 1994)",
            "speedup": "Quadratischer Speedup (nicht exponentiell)",
            "algorithm_steps": [
                "1. Initialisiere gleichverteilten Superpositionszustand |ψ⟩ = (1/√N)·Σ_x|x⟩",
                "2. Wende Orakeloperator U_f an: |x⟩ → (-1)^{f(x)}|x⟩ (invertiert Zielzustand)",
                "3. Wende Diffusionsoperator an: D = 2|ψ⟩⟨ψ| - I (Inversion um Mittelwert)",
                "4. Wiederhole Schritte 2-3 genau ⌊π/4·√N⌋ mal",
                "5. Messung ergibt Zielzustand mit Wahrscheinlichkeit ≥ 1-1/N",
            ],
            "geometric_interpretation": (
                "Rotation im 2D-Unterraum aufgespannt von |Ziel⟩ und |Nicht-Ziel⟩.\n"
                "Pro Iteration: Rotation um Winkel 2θ, wo sin(θ)=1/√N.\n"
                "Nach k Iterationen: Amplitude des Ziels = sin((2k+1)·θ).\n"
                "Optimal bei k ≈ π/(4θ) ≈ π√N/4."
            ),
            "applications": [
                "Beschleunigung von NP-Problemen: O(√{2^n}) statt O(2^n)",
                "Quantenverifizierung von NP-Zertifikaten",
                "Kollisionssuche in O(N^{1/3}) (Brassard et al.)",
                "Amplitude Amplification (Allgemeinform)",
            ],
            "limitations": (
                "Grover-Speedup bricht keine NP-Hürden: 3-SAT mit Grover: O(1.41^n) statt O(2^n).\n"
                "BBBV-Theorem: Jeder Quantenalgorithmus braucht Ω(√N) Orakelabfragen."
            ),
        }

    @staticmethod
    def shors_algorithm_complexity() -> Dict[str, Any]:
        """
        Beschreibt Shors Quantenfaktorisierungsalgorithmus (1994).

        @return Dict mit Algorithmus-Details
        """
        return {
            "algorithm": "Shor's Faktorisierungsalgorithmus",
            "author": "Peter Shor",
            "year": 1994,
            "problem": "INTEGER-FACTORING: Gegeben N = p·q (Semiprüm), finde p und q.",
            "classical_best": "GNFS: exp(O(n^{1/3}·(log n)^{2/3})) — subexponentiell",
            "quantum_complexity": "O(n^2 · log n · log log n) Quantenoperationen — polynomiell!",
            "in_bqp": True,
            "consequence_for_cryptography": (
                "RSA, Diffie-Hellman, ElGamal wären gebrochen durch skalierbaren Quantencomputer.\n"
                "Post-Quanten-Kryptographie (CRYSTALS-KYBER, CRYSTALS-DILITHIUM) als Ersatz."
            ),
            "algorithm_outline": [
                "1. Wähle zufällig a mit 1 < a < N, gcd(a,N)=1",
                "2. Finde Ordnung r von a mod N: a^r ≡ 1 (mod N)",
                "   → Quantenteil: Quantum Phase Estimation auf U|y⟩=|ay mod N⟩",
                "   → QFT (Quantum Fourier Transform) mit O(n^2) Gattern",
                "3. Wenn r gerade und a^{r/2} ≢ -1 (mod N):",
                "   → gcd(a^{r/2}±1, N) ergibt einen echten Teiler",
            ],
            "quantum_components": {
                "QFT": "Quantum Fourier Transform — O(n^2) Gatter, O(n) Tiefe",
                "QPE": "Quantum Phase Estimation — bestimmt Eigenphase auf n Bits genau",
                "Modular_exp": "Modulare Exponentiation — teuerster Teil: O(n^3) Elementaroperationen",
            },
            "implications_for_p_vs_np": (
                "Faktorisierung ∈ BQP zeigt: BQP enthält Probleme, die vermutlich nicht in P liegen.\n"
                "Faktorisierung liegt in NP ∩ co-NP (wahrscheinlich nicht NP-vollständig).\n"
                "Shor beweist NICHT P=NP, sondern zeigt Quantenvorteil für spez. Probleme."
            ),
        }

    @staticmethod
    def quantum_pcp_conjecture() -> str:
        """
        Beschreibt die Quantum PCP Conjecture.

        @return Beschreibungstext
        """
        return (
            "Quantum PCP Conjecture (QPCP-Vermutung):\n\n"
            "Klassisches PCP-Theorem: NP = PCP(O(log n), O(1)).\n"
            "Jeder NP-Beweis kann mit O(log n) Zufallsbits und O(1) Bits verifiziert werden.\n\n"
            "Quantum PCP Vermutung (Aharonov-Arad-Landau-Vazirani):\n"
            "  QMA = QPCP(O(log n), O(1))\n"
            "  Jeder QMA-Beweis (Quantenzustand) kann durch Messung von O(1) Qubits\n"
            "  mit polynomiell vielen klassischen Zufallsbits verifiziert werden.\n\n"
            "Status: Offen (einer der wichtigsten offenen Sätze der Quantenkomplexität)\n\n"
            "Was bekannt ist:\n"
            "  - LOCAL HAMILTONIAN ∈ QMA-vollständig (Kitaev 1999)\n"
            "  - 2-Local Hamiltonian ist QMA-vollständig (Kempe-Kitaev-Regev 2006)\n"
            "  - Ω(√n)-Schranke für 1D-Gitter: Nächster Schritt wäre Ω(n)\n\n"
            "Implikationen wenn QPCP wahr:\n"
            "  - Approximation von Grundzustandsenergien ist QMA-hart\n"
            "  - Starke Verbindung zur Festkörperphysik (Haldane-Gap, NLTS-Theorem)\n"
            "  - Grundlage für Quanten-Inapproximierbarkeit\n\n"
            "NLTS-Theorem (Anshu-Breuckmann-Nirkhe 2022):\n"
            "  Teilfortschritt: Es existieren Hamiltonians, deren Tieftemperaturzustände\n"
            "  nicht durch O(1)-tiefe Quantenschaltkreise approximiert werden können.\n"
            "  Dies ist eine notwendige (aber nicht hinreichende) Bedingung für QPCP."
        )
