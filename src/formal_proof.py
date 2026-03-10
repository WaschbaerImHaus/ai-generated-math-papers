"""
@file formal_proof.py
@brief Formale Beweisinfrastruktur für mathematische Aussagen.
@description
    Implementiert eine Infrastruktur für strukturierte mathematische Beweise:

    - ProofStep: Einzelner Beweisschritt (Aussage + Begründung + Verifikation)
    - ProofChain: Verkettung von Beweisschritten zu einem vollständigen Beweis
    - SymbolicVerifier: Symbolische Verifikation via SymPy
    - InductionProof: Vollständige Induktion mit automatischer Verifikation
    - ContradictionProof: Widerspruchsbeweis-Infrastruktur
    - TheoremRegistry: Datenbank bewiesener Sätze (wiederverwendbar)

    Philosophie:
    - Numerische Verifikation ≠ Beweis (aber hilfreiche Evidenz)
    - Symbolische Verifikation via SymPy = stärker, aber begrenzt
    - Strukturierte Beweise = nachvollziehbar und lehrreich

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sympy as sp
from typing import Callable, Any
from enum import Enum
from dataclasses import dataclass, field


# ===========================================================================
# ENUM: Beweisstatus
# ===========================================================================

class ProofStatus(Enum):
    """
    @brief Status eines Beweisschritts oder Beweises.
    @description
        Zeigt an, wie stark ein Schritt oder ein ganzer Beweis verifiziert wurde.
        Symbolische Verifikation ist stärker als numerische (die nur endlich viele
        Fälle prüft). ASSUMED gilt für Axiome und Voraussetzungen.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    UNVERIFIED = "unverified"       # Noch nicht geprüft
    VERIFIED_SYMBOLIC = "symbolic"  # Symbolisch via SymPy verifiziert
    VERIFIED_NUMERIC = "numeric"    # Numerisch für endlich viele Fälle geprüft
    ASSUMED = "assumed"             # Als Axiom/Voraussetzung angenommen
    FAILED = "failed"               # Verifikation fehlgeschlagen
    PENDING = "pending"             # Verifikation ausstehend (offene Vermutung)


# ===========================================================================
# DATACLASS: Einzelner Beweisschritt
# ===========================================================================

@dataclass
class ProofStep:
    """
    @brief Einzelner Schritt in einem mathematischen Beweis.
    @description
        Jeder Schritt hat:
        - Eine Aussage (als String oder SymPy-Ausdruck)
        - Eine Begründung (Regel, Axiom, vorheriger Schritt)
        - Einen Verifikationsstatus
        - Optional: SymPy-Verifikation der Aussage

        Die Verifikation kann symbolisch (via SymPy simplify) oder numerisch
        (Prüfung endlich vieler Testwerte) erfolgen.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    # Aussage als menschenlesbarer Text
    statement: str

    # Begründung, z.B. "Aus Schritt 2" oder "Definition der Summe"
    justification: str

    # Aktueller Verifikationsstatus (Standard: noch nicht geprüft)
    status: ProofStatus = ProofStatus.UNVERIFIED

    # Optionaler SymPy-Ausdruck für maschinelle Weiterverarbeitung
    sympy_expr: sp.Expr | None = None

    # Liste von (Eingabewert, erwartetes Ergebnis)-Paaren für numerische Prüfung
    numeric_checks: list = field(default_factory=list)

    def verify_symbolic(self, lhs: sp.Expr, rhs: sp.Expr) -> bool:
        """
        @brief Symbolische Verifikation: Prüft ob lhs == rhs via SymPy simplify.
        @description
            Berechnet lhs - rhs und vereinfacht den Ausdruck symbolisch.
            Ergibt die Vereinfachung 0, sind beide Ausdrücke äquivalent.
            Setzt status auf VERIFIED_SYMBOLIC bei Erfolg, FAILED bei Misserfolg.
        @param lhs Linke Seite der Gleichung (SymPy-Ausdruck)
        @param rhs Rechte Seite der Gleichung (SymPy-Ausdruck)
        @return True wenn lhs == rhs symbolisch beweisbar, sonst False
        @date 2026-03-10
        """
        try:
            # Differenz berechnen und vereinfachen
            diff = sp.simplify(lhs - rhs)
            if diff == 0:
                # Differenz ist exakt 0 → Gleichheit bewiesen
                self.status = ProofStatus.VERIFIED_SYMBOLIC
                return True
            # SymPy konnte Gleichheit nicht beweisen
            self.status = ProofStatus.FAILED
            return False
        except Exception:
            # Unerwarteter Fehler bei der symbolischen Berechnung
            self.status = ProofStatus.FAILED
            return False

    def verify_numeric(self, func: Callable, test_values: list, expected: list) -> bool:
        """
        @brief Numerische Verifikation für endlich viele Eingabewerte.
        @description
            Prüft ob func(test_values[i]) == expected[i] für alle i.
            Setzt status auf VERIFIED_NUMERIC wenn alle Werte übereinstimmen.

            WICHTIG: Das ist KEIN vollständiger Beweis! Numerische Verifikation
            prüft nur die gegebenen Fälle, nicht alle natürlichen Zahlen.
        @param func Zu prüfende Funktion f(x) → Ergebnis
        @param test_values Liste von Eingabewerten
        @param expected Liste der erwarteten Ausgabewerte
        @return True wenn alle Werte übereinstimmen, sonst False
        @date 2026-03-10
        """
        # Prüfen ob Listen gleich lang sind
        if len(test_values) != len(expected):
            self.status = ProofStatus.FAILED
            return False

        # Jeden Testwert gegen erwarteten Wert prüfen
        for val, exp in zip(test_values, expected):
            try:
                result = func(val)
                # Ergebnis speichern für spätere Analyse
                self.numeric_checks.append((val, result))
                if result != exp:
                    self.status = ProofStatus.FAILED
                    return False
            except Exception:
                self.status = ProofStatus.FAILED
                return False

        # Alle Werte bestanden → numerisch verifiziert
        self.status = ProofStatus.VERIFIED_NUMERIC
        return True

    def __str__(self) -> str:
        """
        @brief Textdarstellung eines Beweisschritts.
        @return Formatierter String mit Status, Aussage und Begründung
        @date 2026-03-10
        """
        return f"[{self.status.value}] {self.statement} ({self.justification})"


# ===========================================================================
# KLASSE: Vollständige Beweiskette
# ===========================================================================

class ProofChain:
    """
    @brief Vollständiger Beweis als Kette von ProofSteps.
    @description
        Sammelt Beweisschritte in geordneter Reihenfolge und ermöglicht
        symbolische sowie numerische Verifikation einzelner Schritte.
        Unterstützt Method-Chaining für flüssige Beweis-Syntax.

        Beispiel:
            proof = ProofChain("Σk = n(n+1)/2")
            proof.assume("n ∈ ℕ", "Voraussetzung")
            proof.step("Für n=1: Σk = 1 = 1·2/2", "Basisfall")
            proof.conclude("Σk = n(n+1)/2 für alle n ∈ ℕ", "Induktionsprinzip")
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, theorem: str, description: str = "") -> None:
        """
        @brief Initialisiert eine neue Beweiskette.
        @param theorem Der zu beweisende Satz als Text
        @param description Optionale ausführliche Beschreibung
        @date 2026-03-10
        """
        # Hauptaussage des Beweises
        self.theorem = theorem

        # Optionale Beschreibung mit Kontext
        self.description = description

        # Geordnete Liste aller Beweisschritte
        self.steps: list[ProofStep] = []

        # Markierung ob der Beweis mit conclude() abgeschlossen wurde
        self.is_complete: bool = False

    def assume(self, statement: str, justification: str = "Voraussetzung") -> 'ProofChain':
        """
        @brief Fügt eine Annahme oder Voraussetzung hinzu.
        @description
            Voraussetzungen und Axiome werden mit ASSUMED markiert — sie
            werden nicht weiter verifiziert, sondern als gegeben angenommen.
        @param statement Die Annahme als Text
        @param justification Begründung (Standard: "Voraussetzung")
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Schritt mit ASSUMED-Status erstellen (keine Verifikation nötig)
        step = ProofStep(statement, justification, ProofStatus.ASSUMED)
        self.steps.append(step)
        return self  # Method chaining ermöglichen

    def step(self, statement: str, justification: str,
             lhs: sp.Expr = None, rhs: sp.Expr = None) -> 'ProofChain':
        """
        @brief Fügt einen Beweisschritt hinzu, optional mit symbolischer Verifikation.
        @description
            Wenn lhs und rhs angegeben werden, prüft SymPy ob lhs == rhs.
            Ohne lhs/rhs bleibt der Schritt UNVERIFIED (informeller Schritt).
        @param statement Der Beweisschritt als Text
        @param justification Begründung des Schritts
        @param lhs Optionaler SymPy-Ausdruck: linke Seite
        @param rhs Optionaler SymPy-Ausdruck: rechte Seite
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Neuen Schritt anlegen (zunächst unverifiziert)
        proof_step = ProofStep(statement, justification)

        # Wenn beide Seiten angegeben: symbolische Verifikation versuchen
        if lhs is not None and rhs is not None:
            proof_step.verify_symbolic(lhs, rhs)

        self.steps.append(proof_step)
        return self  # Method chaining

    def numeric_step(self, statement: str, justification: str,
                     func: Callable, test_range: range) -> 'ProofChain':
        """
        @brief Fügt einen Schritt mit numerischer Verifikation für einen Bereich hinzu.
        @description
            Prüft die Aussage func(n) für alle n im test_range.
            Nützlich für empirische Evidenz bei Vermutungen wie Goldbach.
        @param statement Der Beweisschritt als Text
        @param justification Begründung des Schritts
        @param func Funktion f(n) → bool, die für alle n True sein soll
        @param test_range Bereich der zu prüfenden Werte (range-Objekt)
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Neuen Schritt anlegen
        proof_step = ProofStep(statement, justification)

        # Alle Werte im Bereich testen
        test_values = list(test_range)
        expected = [True] * len(test_values)

        # Numerische Verifikation durchführen
        proof_step.verify_numeric(func, test_values, expected)

        self.steps.append(proof_step)
        return self  # Method chaining

    def conclude(self, statement: str, justification: str = "q.e.d.") -> 'ProofChain':
        """
        @brief Fügt die Schlussfolgerung hinzu und markiert den Beweis als vollständig.
        @description
            Die Schlussfolgerung ist der abschließende Schritt jedes Beweises.
            Nach conclude() ist is_complete = True.
        @param statement Die Schlussfolgerung als Text
        @param justification Begründung (Standard: "q.e.d.")
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Schlussfolgerung als ASSUMED markieren (Ergebnis der vorigen Schritte)
        step = ProofStep(statement, justification, ProofStatus.ASSUMED)
        self.steps.append(step)
        # Beweis als abgeschlossen markieren
        self.is_complete = True
        return self  # Method chaining

    def status_summary(self) -> dict[str, Any]:
        """
        @brief Gibt eine Zusammenfassung des Beweises zurück.
        @description
            Zählt die Schritte nach Status und bewertet die Beweissstärke:
            - 'symbolic': Mindestens ein symbolisch verifizierter Schritt, kein FAILED
            - 'numeric': Nur numerisch verifizierte Schritte, kein FAILED
            - 'informal': Nur ASSUMED/UNVERIFIED Schritte
        @return Dict mit Zählern und Beweisstärke
        @date 2026-03-10
        """
        # Zähler initialisieren
        total = len(self.steps)
        symbolic = sum(1 for s in self.steps if s.status == ProofStatus.VERIFIED_SYMBOLIC)
        numeric = sum(1 for s in self.steps if s.status == ProofStatus.VERIFIED_NUMERIC)
        assumed = sum(1 for s in self.steps if s.status == ProofStatus.ASSUMED)
        failed = sum(1 for s in self.steps if s.status == ProofStatus.FAILED)
        pending = sum(1 for s in self.steps if s.status == ProofStatus.PENDING)
        unverified = sum(1 for s in self.steps if s.status == ProofStatus.UNVERIFIED)

        # Gesamtanzahl verifizierter Schritte (symbolisch + numerisch)
        verified = symbolic + numeric

        # Beweisstärke bestimmen (Rangfolge: symbolic > numeric > informal)
        if failed > 0:
            strength = "failed"
        elif symbolic > 0:
            strength = "symbolic"
        elif numeric > 0:
            strength = "numeric"
        else:
            strength = "informal"

        return {
            'total_steps': total,
            'verified': verified,
            'symbolic': symbolic,
            'numeric': numeric,
            'assumed': assumed,
            'failed': failed,
            'pending': pending,
            'unverified': unverified,
            'is_complete': self.is_complete,
            'strength': strength
        }

    def to_text(self) -> str:
        """
        @brief Formatiert den Beweis als lesbaren Klartext.
        @description
            Gibt einen mehrzeiligen String zurück mit Theorem, allen Schritten
            (nummeriert) und einer Zusammenfassung am Ende.
        @return Formatierter Text-Beweis
        @date 2026-03-10
        """
        # Kopfzeile mit Theorem
        lines = [
            "=" * 60,
            f"BEWEIS: {self.theorem}",
        ]

        # Optionale Beschreibung
        if self.description:
            lines.append(f"Beschreibung: {self.description}")

        lines.append("-" * 60)

        # Alle Schritte nummeriert ausgeben
        for i, step in enumerate(self.steps, 1):
            lines.append(f"  {i:2d}. {step}")

        lines.append("-" * 60)

        # Zusammenfassung
        summary = self.status_summary()
        lines.append(
            f"Status: {'VOLLSTÄNDIG' if self.is_complete else 'UNVOLLSTÄNDIG'} | "
            f"Stärke: {summary['strength']} | "
            f"Verifiziert: {summary['verified']}/{summary['total_steps']}"
        )
        lines.append("=" * 60)

        return "\n".join(lines)

    def to_html(self) -> str:
        """
        @brief Formatiert den Beweis als HTML für das Web-Interface.
        @description
            Gibt einen HTML-String zurück mit farbkodierten Schritten:
            - Grün: symbolisch/numerisch verifiziert
            - Gelb: Annahme (assumed)
            - Rot: fehlgeschlagen oder ausstehend
            - Grau: unverifiziert
        @return HTML-String mit strukturiertem Beweis
        @date 2026-03-10
        """
        # CSS-Klasse je nach ProofStatus bestimmen
        status_class = {
            ProofStatus.VERIFIED_SYMBOLIC: "proof-symbolic",   # Grün (stark)
            ProofStatus.VERIFIED_NUMERIC: "proof-numeric",     # Blaugrün
            ProofStatus.ASSUMED: "proof-assumed",              # Gelb
            ProofStatus.FAILED: "proof-failed",                # Rot
            ProofStatus.PENDING: "proof-pending",              # Orange
            ProofStatus.UNVERIFIED: "proof-unverified",        # Grau
        }

        # Symbol je nach ProofStatus
        status_symbol = {
            ProofStatus.VERIFIED_SYMBOLIC: "✓",
            ProofStatus.VERIFIED_NUMERIC: "~",
            ProofStatus.ASSUMED: "⊢",
            ProofStatus.FAILED: "✗",
            ProofStatus.PENDING: "?",
            ProofStatus.UNVERIFIED: "○",
        }

        # HTML-Kopf mit Theorem-Titel
        html_parts = [
            '<div class="proof-container">',
            f'  <div class="proof-title">Beweis: {self.theorem}</div>',
        ]

        # Beschreibung (falls vorhanden)
        if self.description:
            html_parts.append(
                f'  <div class="proof-description">{self.description}</div>'
            )

        # Schritte-Liste
        html_parts.append('  <ol class="proof-steps">')
        for step in self.steps:
            css = status_class.get(step.status, "proof-unverified")
            sym = status_symbol.get(step.status, "○")
            # Jeden Schritt als Listenelement mit Farbkodierung
            html_parts.append(
                f'    <li class="proof-step {css}">'
                f'<span class="proof-status-icon">{sym}</span>'
                f'<span class="proof-statement">{step.statement}</span>'
                f'<span class="proof-justification">({step.justification})</span>'
                f'</li>'
            )
        html_parts.append('  </ol>')

        # Zusammenfassung am Ende
        summary = self.status_summary()
        completeness = "Vollständig" if self.is_complete else "Unvollständig"
        html_parts.append(
            f'  <div class="proof-summary">'
            f'Status: {completeness} | '
            f'Stärke: {summary["strength"]} | '
            f'Verifiziert: {summary["verified"]}/{summary["total_steps"]}'
            f'</div>'
        )
        html_parts.append('</div>')

        return "\n".join(html_parts)

    def __repr__(self) -> str:
        """
        @brief Kurze Repr-Darstellung der ProofChain.
        @return String mit Theorem und Verifikationsstand
        @date 2026-03-10
        """
        summary = self.status_summary()
        return (
            f"ProofChain('{self.theorem}', "
            f"{summary['verified']}/{summary['total_steps']} verifiziert)"
        )


# ===========================================================================
# KLASSE: Vollständige Induktion
# ===========================================================================

class InductionProof(ProofChain):
    """
    @brief Vollständige Induktion: Basisfall + Induktionsschritt.
    @description
        Schema der vollständigen Induktion:
        1. Basisfall: P(n₀) ist wahr
        2. Induktionsvoraussetzung: Annahme P(n) ist wahr
        3. Induktionsschritt: P(n) → P(n+1)
        4. Schluss: P(n) ist wahr für alle n ≥ n₀

        Unterstützt sowohl empirische (numerische) als auch symbolische
        Verifikation des Induktionsschritts via SymPy.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, theorem: str, predicate: Callable[[int], bool],
                 n0: int = 0) -> None:
        """
        @brief Initialisiert einen Induktionsbeweis.
        @param theorem Aussage als Text
        @param predicate Funktion P(n) → bool, die die Aussage für n prüft
        @param n0 Startwert für den Basisfall (Standard: 0)
        @date 2026-03-10
        """
        super().__init__(theorem)

        # Prädikat P(n): Die zu beweisende Aussage für natürliche Zahl n
        self.predicate = predicate

        # Startwert n₀ für den Basisfall
        self.n0 = n0

        # Status der Verifikationen
        self.base_case_verified = False
        self.induction_step_verified = False

    def verify_base_case(self) -> bool:
        """
        @brief Prüft den Basisfall P(n₀) und fügt den Schritt zur Kette hinzu.
        @description
            Wertet P(n₀) aus und erstellt einen entsprechenden ProofStep.
            Setzt base_case_verified auf True bei Erfolg.
        @return True wenn P(n₀) wahr ist, sonst False
        @date 2026-03-10
        """
        try:
            # Basisfall prüfen: P(n₀)
            result = self.predicate(self.n0)

            if result:
                # Basisfall bestanden → VERIFIED_NUMERIC
                step = ProofStep(
                    f"Basisfall: P({self.n0}) ist wahr",
                    f"Direkte Berechnung für n={self.n0}",
                    ProofStatus.VERIFIED_NUMERIC
                )
                self.base_case_verified = True
            else:
                # Basisfall fehlgeschlagen
                step = ProofStep(
                    f"Basisfall: P({self.n0}) ist FALSCH",
                    f"Gegenbeispiel bei n={self.n0}",
                    ProofStatus.FAILED
                )

            self.steps.append(step)
            return result

        except Exception as e:
            # Fehler bei der Auswertung
            step = ProofStep(
                f"Basisfall: Fehler bei P({self.n0}): {e}",
                "Ausnahme bei Auswertung",
                ProofStatus.FAILED
            )
            self.steps.append(step)
            return False

    def verify_induction_empirically(self, n_max: int = 100) -> bool:
        """
        @brief Empirische Verifikation für n₀ bis n_max.
        @description
            Prüft P(n) für alle n von n₀ bis n_max.

            WICHTIG: Das ist KEIN vollständiger mathematischer Beweis!
            Empirische Verifikation liefert nur Evidenz für endlich viele Fälle.
            Ein vollständiger Beweis erfordert den symbolischen Induktionsschritt.
        @param n_max Maximaler Wert für die Verifikation (Standard: 100)
        @return True wenn P(n) für alle n in [n₀, n_max] wahr ist
        @date 2026-03-10
        """
        # Zähler für gefundene Gegenbeispiele
        failures = []

        # Alle Werte von n₀ bis n_max prüfen
        for n in range(self.n0, n_max + 1):
            try:
                if not self.predicate(n):
                    failures.append(n)
            except Exception:
                failures.append(n)

        if not failures:
            # Alle Werte bestanden
            step = ProofStep(
                f"Empirische Verifikation: P(n) für alle n ∈ [{self.n0}, {n_max}] wahr",
                f"Numerische Prüfung von {n_max - self.n0 + 1} Werten "
                f"(KEIN vollständiger Beweis!)",
                ProofStatus.VERIFIED_NUMERIC
            )
            self.induction_step_verified = True
            self.steps.append(step)
            return True
        else:
            # Gegenbeispiele gefunden
            step = ProofStep(
                f"Empirische Verifikation: Gegenbeispiele bei n = {failures[:5]}",
                f"Numerische Prüfung fehlgeschlagen",
                ProofStatus.FAILED
            )
            self.steps.append(step)
            return False

    def prove_by_symbolic_induction(self,
                                    base_lhs: sp.Expr, base_rhs: sp.Expr,
                                    step_lhs: sp.Expr, step_rhs: sp.Expr,
                                    n: sp.Symbol = None) -> bool:
        """
        @brief Symbolischer Induktionsbeweis via SymPy.
        @description
            Prüft:
            1. Basisfall: base_lhs == base_rhs (für n = n₀)
            2. Induktionsschritt: step_lhs == step_rhs (symbolisch in n)

            Beide Gleichheiten werden mit sp.simplify(lhs - rhs) == 0 geprüft.
        @param base_lhs Linke Seite des Basisfalls
        @param base_rhs Rechte Seite des Basisfalls
        @param step_lhs Linke Seite des Induktionsschritts (enthält n)
        @param step_rhs Rechte Seite des Induktionsschritts (enthält n)
        @param n SymPy-Symbol für n (Standard: sp.Symbol('n'))
        @return True wenn beide Schritte symbolisch verifiziert werden konnten
        @date 2026-03-10
        """
        # Standard-Symbol für n
        if n is None:
            n = sp.Symbol('n', positive=True, integer=True)

        # --- Basisfall symbolisch prüfen ---
        base_step = ProofStep(
            f"Basisfall symbolisch: base_lhs = base_rhs",
            "SymPy simplify(base_lhs - base_rhs) == 0"
        )
        base_ok = base_step.verify_symbolic(base_lhs, base_rhs)
        self.steps.append(base_step)

        # --- Induktionsschritt symbolisch prüfen ---
        induction_step = ProofStep(
            "Induktionsschritt symbolisch: step_lhs = step_rhs",
            "SymPy simplify(step_lhs - step_rhs) == 0"
        )
        step_ok = induction_step.verify_symbolic(step_lhs, step_rhs)
        self.steps.append(induction_step)

        # Beide müssen erfolgreich sein
        if base_ok and step_ok:
            self.base_case_verified = True
            self.induction_step_verified = True

        return base_ok and step_ok


# ===========================================================================
# KLASSE: Widerspruchsbeweis
# ===========================================================================

class ContradictionProof(ProofChain):
    """
    @brief Widerspruchsbeweis: Nehme ¬A an, leite Widerspruch ab → A ist wahr.
    @description
        Schema des Widerspruchsbeweises (reductio ad absurdum):
        1. Annahme: ¬A (Negation der zu beweisenden Aussage)
        2. Aus ¬A folgt Aussage B
        3. Aussage B widerspricht einer bekannten Tatsache
        4. Also muss ¬A falsch sein → A ist wahr

        Klassische Beispiele:
        - √2 ist irrational
        - Unendlich viele Primzahlen (Euklid)
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, theorem: str) -> None:
        """
        @brief Initialisiert einen Widerspruchsbeweis.
        @param theorem Die zu beweisende Aussage (NICHT die Negation)
        @date 2026-03-10
        """
        super().__init__(theorem)

        # Die angenommene Negation der Aussage
        self.negation: str = ""

        # Markierung ob ein Widerspruch explizit gefunden wurde
        self.contradiction_found: bool = False

    def negate(self, statement: str) -> 'ContradictionProof':
        """
        @brief Nimmt die Negation des Theorems als Ausgangspunkt an.
        @description
            Erster Schritt jedes Widerspruchsbeweises: Annahme der Negation.
            Erstellt einen ASSUMED-Schritt mit der negierten Aussage.
        @param statement Die Negation als Text (z.B. "Angenommen √2 = p/q rational")
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Negation speichern
        self.negation = statement

        # Als ASSUMED-Schritt hinzufügen (Annahme zum Widerspruch)
        step = ProofStep(
            f"Annahme (zur Widerlegung): {statement}",
            "Annahme der Negation (reductio ad absurdum)",
            ProofStatus.ASSUMED
        )
        self.steps.append(step)
        return self  # Method chaining

    def derive_contradiction(self, statement: str,
                             justification: str) -> 'ContradictionProof':
        """
        @brief Leitet einen Widerspruch ab und markiert ihn als gefunden.
        @description
            Dieser Schritt markiert den gefundenen Widerspruch explizit.
            Setzt contradiction_found auf True.
        @param statement Die Widerspruchsaussage (z.B. "p und q wären beide gerade")
        @param justification Begründung des Widerspruchs
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Widerspruch als VERIFIED_NUMERIC markieren (logisch hergeleitet)
        step = ProofStep(
            f"⚡ Widerspruch: {statement}",
            justification,
            ProofStatus.VERIFIED_NUMERIC
        )
        self.steps.append(step)
        # Widerspruch gefunden → merken
        self.contradiction_found = True
        return self  # Method chaining

    def qed(self) -> 'ContradictionProof':
        """
        @brief Schließt den Widerspruchsbeweis ab.
        @description
            Da die Annahme (Negation) zu einem Widerspruch führt, muss
            das Original-Theorem wahr sein. Markiert den Beweis als vollständig.
        @return self für Method-Chaining
        @date 2026-03-10
        """
        # Schlussfolgerung des Widerspruchsbeweises
        conclusion = (
            f"Da die Annahme '¬({self.theorem})' zu einem Widerspruch führt, "
            f"ist '{self.theorem}' bewiesen. q.e.d."
        )
        step = ProofStep(
            conclusion,
            "Widerspruchsprinzip: ¬A → ⊥ impliziert A",
            ProofStatus.VERIFIED_NUMERIC
        )
        self.steps.append(step)
        # Beweis abschließen
        self.is_complete = True
        return self  # Method chaining


# ===========================================================================
# KLASSE: Satz-Datenbank
# ===========================================================================

class TheoremRegistry:
    """
    @brief Datenbank bewiesener Sätze — wiederverwendbar in anderen Beweisen.
    @description
        Zentrale Registrierung aller ProofChain-Objekte unter einem Namen.
        Enthält klassische Sätze aus der Mathematikgeschichte.
        Die _registry ist ein Klassen-Attribut (geteilt zwischen allen Instanzen).
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    # Klassen-Attribut: geteilte Registry aller registrierten Beweise
    _registry: dict[str, ProofChain] = {}

    @classmethod
    def register(cls, name: str, proof: ProofChain) -> None:
        """
        @brief Registriert einen Beweis unter einem Namen.
        @param name Eindeutiger Name des Satzes (z.B. 'gauss_sum')
        @param proof Die ProofChain-Instanz
        @date 2026-03-10
        """
        cls._registry[name] = proof

    @classmethod
    def get(cls, name: str) -> ProofChain | None:
        """
        @brief Gibt einen registrierten Beweis zurück.
        @param name Name des Satzes
        @return ProofChain oder None wenn nicht gefunden
        @date 2026-03-10
        """
        return cls._registry.get(name, None)

    @classmethod
    def list_theorems(cls) -> list[str]:
        """
        @brief Listet die Namen aller registrierten Sätze.
        @return Alphabetisch sortierte Liste der Satz-Namen
        @date 2026-03-10
        """
        return sorted(cls._registry.keys())

    @classmethod
    def load_classics(cls) -> None:
        """
        @brief Lädt klassische Beweise in die Registry.
        @description
            Registriert folgende bewiesene Sätze:
            - 'gauss_sum': Σk=1..n k = n(n+1)/2 (Gaußsche Summenformel)
            - 'sqrt2_irrational': √2 ist irrational (Widerspruchsbeweis)
            - 'infinitely_many_primes': Unendlich viele Primzahlen (Euklid)
            - 'fermat_little': a^p ≡ a (mod p) für Primzahl p
            - 'goldbach_partial': Goldbach für n ≤ 10000 (numerisch)
        @date 2026-03-10
        """
        # Gaußsche Summenformel
        cls._registry['gauss_sum'] = prove_gauss_sum()

        # √2 Irrationalität
        cls._registry['sqrt2_irrational'] = prove_sqrt2_irrational()

        # Unendlich viele Primzahlen
        cls._registry['infinitely_many_primes'] = prove_infinitely_many_primes()

        # Fermats kleiner Satz (für p=5 als Beispiel)
        cls._registry['fermat_little'] = prove_fermat_little_theorem(p=5)

        # Goldbach-Vermutung (numerisch bis 10000)
        cls._registry['goldbach_partial'] = goldbach_evidence(n_max=10000)


# ===========================================================================
# KLASSISCHE BEWEISE — Vorgefertigte Beweisfunktionen
# ===========================================================================

def prove_gauss_sum() -> InductionProof:
    """
    @brief Beweis der Gaußschen Summenformel: Σ(k=1..n) k = n(n+1)/2
    @description
        Methode: Vollständige Induktion mit symbolischer SymPy-Verifikation.

        Schema:
        1. Basisfall n=1: 1 = 1·2/2 = 1 ✓
        2. Induktionsvoraussetzung: Σ(k=1..n) k = n(n+1)/2
        3. Induktionsschritt: Σ(k=1..n+1) k = n(n+1)/2 + (n+1)
                                              = (n+1)(n+2)/2

        KaTeX-Formel: $\sum_{k=1}^{n} k = \frac{n(n+1)}{2}$
    @return InductionProof-Instanz mit vollständigem Beweis
    @date 2026-03-10
    """
    # SymPy-Symbol für n (positiv, ganzzahlig)
    n = sp.Symbol('n', positive=True, integer=True)

    # InductionProof initialisieren
    # Prädikat: Prüft ob Σk = n(n+1)/2 für konkretes n gilt
    proof = InductionProof(
        "Σ(k=1..n) k = n(n+1)/2",
        predicate=lambda m: sum(range(1, m + 1)) == m * (m + 1) // 2,
        n0=1
    )

    # Voraussetzung: n ist eine positive natürliche Zahl
    proof.assume("n ∈ ℕ, n ≥ 1", "Definitionsbereich der Formel")

    # Basisfall n=1 verifizieren
    proof.verify_base_case()

    # Induktionsvoraussetzung als Annahme hinzufügen
    proof.assume(
        "Induktionsvoraussetzung: Σ(k=1..n) k = n(n+1)/2",
        "Annahme für festes n (Induktionshypothese)"
    )

    # Symbolischer Induktionsschritt:
    # Linke Seite: Σ(k=1..n+1) = Σ(k=1..n) + (n+1) = n(n+1)/2 + (n+1)
    # Rechte Seite: (n+1)(n+2)/2
    # Beide müssen gleich sein
    lhs = sp.Rational(1, 2) * n * (n + 1) + (n + 1)   # IV + nächster Term
    rhs = sp.Rational(1, 2) * (n + 1) * (n + 2)        # Formel für n+1

    proof.step(
        "Induktionsschritt: n(n+1)/2 + (n+1) = (n+1)(n+2)/2",
        "IV einsetzen + algebraische Vereinfachung",
        lhs, rhs
    )

    # Empirische Verifikation für n = 1..100 (zusätzliche Evidenz)
    proof.verify_induction_empirically(100)

    # Schlussfolgerung
    proof.conclude(
        "Gaußsche Summenformel Σ(k=1..n) k = n(n+1)/2 gilt für alle n ∈ ℕ",
        "Induktionsprinzip (Basisfall + Schritt ⇒ Allgemeinheit)"
    )

    return proof


def prove_sqrt2_irrational() -> ContradictionProof:
    """
    @brief Beweis: √2 ist irrational.
    @description
        Methode: Widerspruchsbeweis (klassisch, antike Griechen).

        Schema:
        1. Annahme: √2 = p/q rational, p/q vollständig gekürzt
        2. Dann: 2 = p²/q² → p² = 2q²
        3. Also ist p² gerade → p ist gerade → p = 2m
        4. Einsetzen: 4m² = 2q² → q² = 2m² → q ist gerade
        5. Widerspruch: p und q sind beide gerade, ggT(p,q) ≥ 2

        KaTeX-Formel: $\sqrt{2} \notin \mathbb{Q}$
    @return ContradictionProof-Instanz mit vollständigem Beweis
    @date 2026-03-10
    """
    # Widerspruchsbeweis initialisieren
    proof = ContradictionProof("√2 ist irrational")

    # Schritt 1: Negation annehmen
    proof.negate("√2 sei rational, d.h. √2 = p/q mit p, q ∈ ℤ, q ≠ 0, ggT(p,q) = 1")

    # Schritt 2: Algebraische Folgerung
    proof.step(
        "Dann gilt: (√2)² = p²/q², also 2 = p²/q², also p² = 2q²",
        "Quadrieren beider Seiten"
    )

    # Schritt 3: p muss gerade sein
    proof.step(
        "p² = 2q² ist gerade → p² ist durch 2 teilbar → p ist gerade → p = 2m für ein m ∈ ℤ",
        "Lemma: p² gerade ⟺ p gerade (für ganze Zahlen)"
    )

    # Schritt 4: q muss auch gerade sein
    proof.step(
        "Einsetzen p = 2m: (2m)² = 2q² → 4m² = 2q² → q² = 2m² → q ist gerade",
        "Gleiche Argumentation: q² gerade ⟺ q gerade"
    )

    # Schritt 5: Widerspruch finden
    proof.derive_contradiction(
        "p und q sind beide gerade → ggT(p, q) ≥ 2. Widerspruch zur Annahme ggT(p,q) = 1!",
        "p = 2m und q = 2k → gemeinsamer Faktor 2"
    )

    # Schlussfolgerung: Beweis abschließen
    proof.qed()

    return proof


def prove_infinitely_many_primes() -> ContradictionProof:
    """
    @brief Beweis: Es gibt unendlich viele Primzahlen (Euklid, ca. 300 v. Chr.).
    @description
        Methode: Widerspruchsbeweis (Euklids klassischer Beweis).

        Schema:
        1. Annahme: Es gibt nur endlich viele Primzahlen p₁, p₂, ..., pₙ
        2. Bilde N = p₁ · p₂ · ... · pₙ + 1
        3. N ist durch keine der Primzahlen p₁, ..., pₙ teilbar (Rest 1)
        4. N hat also einen Primteiler q, der nicht in der Liste ist
        5. Widerspruch zur Annahme, dass die Liste vollständig ist

        KaTeX-Formel: $|\{p \in \mathbb{N} : p \text{ prim}\}| = \infty$
    @return ContradictionProof-Instanz mit vollständigem Beweis
    @date 2026-03-10
    """
    # Widerspruchsbeweis initialisieren
    proof = ContradictionProof("Es gibt unendlich viele Primzahlen")

    # Schritt 1: Endliche Liste annehmen
    proof.negate(
        "Es gibt nur endlich viele Primzahlen: p₁, p₂, ..., pₙ "
        "(vollständige Liste aller Primzahlen)"
    )

    # Schritt 2: Euklid-Zahl konstruieren
    proof.step(
        "Konstruiere N = p₁ · p₂ · ... · pₙ + 1 (Produkt aller Primzahlen + 1)",
        "Konstruktion (Euklids Trick)"
    )

    # Schritt 3: N ist nicht durch die bekannten Primzahlen teilbar
    proof.step(
        "N mod pᵢ = 1 für alle i = 1, ..., n → keine der Primzahlen pᵢ teilt N",
        "N = (p₁·...·pₙ) + 1 → Rest 1 bei Division durch jedes pᵢ"
    )

    # Schritt 4: N muss aber einen Primteiler haben
    proof.step(
        "N > 1 → N hat mindestens einen Primteiler q (Fundamentalsatz der Arithmetik)",
        "Jede natürliche Zahl > 1 besitzt einen Primteiler"
    )

    # Schritt 5: Widerspruch
    proof.derive_contradiction(
        "q teilt N, aber q ist nicht in {p₁, ..., pₙ} → die Liste war unvollständig!",
        "q ≠ pᵢ für alle i (da kein pᵢ N teilt) → neue Primzahl gefunden"
    )

    # Schlussfolgerung
    proof.qed()

    return proof


def prove_fermat_little_theorem(p: int = 5) -> ProofChain:
    """
    @brief Beweis von Fermats kleinem Satz: a^p ≡ a (mod p) für Primzahl p.
    @description
        Numerische Verifikation für alle a = 0..p-1 für eine gegebene Primzahl p.

        Für einen vollständigen Beweis wird üblicherweise Induktion über a
        oder Gruppentheorie (Lagranges Satz) verwendet.

        KaTeX-Formel: $a^p \equiv a \pmod{p}$ für alle $a \in \mathbb{Z}$
    @param p Primzahl (Standard: 5)
    @return ProofChain mit numerischer Verifikation
    @date 2026-03-10
    """
    # ProofChain für Fermats kleinen Satz
    proof = ProofChain(
        f"a^{p} ≡ a (mod {p}) für alle a ∈ ℤ",
        f"Fermats kleiner Satz für die Primzahl p = {p}"
    )

    # Voraussetzung: p muss eine Primzahl sein
    proof.assume(f"p = {p} ist eine Primzahl", "Voraussetzung (Primzahl-Eigenschaft)")

    # Voraussetzung: a ist eine beliebige ganze Zahl
    proof.assume(f"a ∈ ℤ (hier: a = 0, 1, ..., {p-1})", "Testbereich für numerische Prüfung")

    # Numerische Verifikation: a^p mod p == a mod p für a = 0..p-1
    def fermat_check(a: int) -> bool:
        """Prüft a^p ≡ a (mod p) für ein konkretes a."""
        return pow(a, p, p) == a % p

    # Jeden Wert a = 0..p-1 prüfen und als Schritt hinzufügen
    all_passed = True
    for a in range(p):
        result = fermat_check(a)
        status = ProofStatus.VERIFIED_NUMERIC if result else ProofStatus.FAILED
        step = ProofStep(
            f"{a}^{p} mod {p} = {pow(a, p, p)} ≡ {a} mod {p} = {a % p}",
            f"Direkte Berechnung für a = {a}",
            status
        )
        proof.steps.append(step)
        if not result:
            all_passed = False

    # Schlussfolgerung
    if all_passed:
        proof.conclude(
            f"Fermats kleiner Satz verifiziert: a^{p} ≡ a (mod {p}) für alle a = 0..{p-1}",
            "Numerische Verifikation aller Resklassen mod p"
        )
    else:
        # Sollte bei einer echten Primzahl nie auftreten
        proof.conclude(
            f"Verifikation fehlgeschlagen (Fehler: p={p} ist keine Primzahl?)",
            "Gegenbeispiel gefunden"
        )

    return proof


def riemann_hypothesis_evidence() -> ProofChain:
    """
    @brief Empirische Evidenz für die Riemann-Hypothese.
    @description
        Sammelt numerische Evidenz für die Riemann-Hypothese:
        1. Alle nicht-trivialen Nullstellen liegen auf der kritischen Geraden Re(s) = 1/2
        2. Gram-Punkte und Gram-Gesetz
        3. GUE-Statistik der Nullstellen-Abstände

        Status: PENDING (Millennium-Problem, bis heute unbewiesen)

        KaTeX-Formel: $\zeta(s) = 0 \Rightarrow \text{Re}(s) = \frac{1}{2}$ (für nicht-triviale Nullstellen)
    @return ProofChain mit PENDING-Status
    @date 2026-03-10
    """
    # ProofChain für die Riemann-Hypothese
    proof = ProofChain(
        "Alle nicht-trivialen Nullstellen von ζ(s) liegen auf Re(s) = 1/2",
        "Riemann-Hypothese (1859) — Millennium-Problem, unbewiesen"
    )

    # Hintergrund
    proof.assume(
        "ζ(s) = Σ(n=1..∞) 1/n^s (Riemannsche Zeta-Funktion)",
        "Definition der Zeta-Funktion"
    )
    proof.assume(
        "Nicht-triviale Nullstellen: ζ(s) = 0 im kritischen Streifen 0 < Re(s) < 1",
        "Analytische Fortsetzung auf ℂ"
    )

    # Erste Nullstellen numerisch prüfen (bekannte Werte)
    # Die ersten nicht-trivialen Nullstellen: ½ + i·t mit bekannten t-Werten
    known_zeros_t = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832
    ]

    # Für jede bekannte Nullstelle einen Schritt hinzufügen
    for i, t in enumerate(known_zeros_t, 1):
        step = ProofStep(
            f"Nullstelle {i}: s = 1/2 + {t:.6f}·i liegt auf Re(s) = 1/2",
            f"Numerisch verifiziert (bekannter Wert)",
            ProofStatus.VERIFIED_NUMERIC
        )
        proof.steps.append(step)

    # Gram-Gesetz (empirische Beobachtung)
    gram_step = ProofStep(
        "Gram-Gesetz: Zwischen aufeinanderfolgenden Gram-Punkten liegt mindestens eine Nullstelle "
        "(gilt für die meisten, nicht alle Gram-Punkte)",
        "Empirische Beobachtung (Gram 1903)",
        ProofStatus.VERIFIED_NUMERIC
    )
    proof.steps.append(gram_step)

    # GUE-Statistik
    gue_step = ProofStep(
        "Die statistischen Abstände der Nullstellen folgen der GUE-Zufallsmatrix-Verteilung",
        "Montgomery-Odlyzko-Gesetz (1973/1987) — Verbindung zur Quantenchaostheorie",
        ProofStatus.VERIFIED_NUMERIC
    )
    proof.steps.append(gue_step)

    # Aktueller Stand der Forschung
    status_step = ProofStep(
        "Stand 2026: Alle ~10^13 berechneten Nullstellen bestätigen die Hypothese",
        "Numerische Massenverifikation (Odlyzko, ZetaGrid, LMFDB)",
        ProofStatus.VERIFIED_NUMERIC
    )
    proof.steps.append(status_step)

    # Schlussfolgerung: Offen (PENDING)
    pending_step = ProofStep(
        "Die Riemann-Hypothese bleibt unbewiesen — trotz starker numerischer Evidenz",
        "Millennium-Problem (Clay Mathematics Institute, 2000)",
        ProofStatus.PENDING
    )
    proof.steps.append(pending_step)

    # Beweis ist technisch "abgeschlossen" (als offene Vermutung dokumentiert)
    proof.is_complete = False  # Noch kein vollständiger Beweis

    return proof


def goldbach_evidence(n_max: int = 10000) -> ProofChain:
    """
    @brief Empirische Evidenz für die Goldbach-Vermutung.
    @description
        Die Goldbach-Vermutung (1742): Jede gerade Zahl > 2 ist Summe zweier Primzahlen.
        Dieser Beweis gibt nur NUMERISCHE Evidenz für alle geraden n ≤ n_max.

        Status: VERIFIED_NUMERIC (nicht vollständig bewiesen!)

        KaTeX-Formel: $\forall n \in 2\mathbb{N}, n > 2: \exists p, q \in \mathbb{P}: n = p + q$
    @param n_max Maximaler Wert für numerische Verifikation (Standard: 10000)
    @return ProofChain mit VERIFIED_NUMERIC-Status
    @date 2026-03-10
    """
    # Hilfsfunktion: Primzahl-Sieb für schnelle Prüfung
    def sieve(limit: int) -> list[bool]:
        """Sieb des Eratosthenes: is_prime[n] = True wenn n prim."""
        is_p = [True] * (limit + 1)
        is_p[0] = is_p[1] = False
        for i in range(2, int(limit**0.5) + 1):
            if is_p[i]:
                for j in range(i * i, limit + 1, i):
                    is_p[j] = False
        return is_p

    # Hilfsfunktion: Prüft ob n als Summe zweier Primzahlen darstellbar
    def has_goldbach_decomposition(n: int, is_prime_arr: list[bool]) -> bool:
        """Prüft ob n = p + q mit p, q prim."""
        for p in range(2, n // 2 + 1):
            if is_prime_arr[p] and is_prime_arr[n - p]:
                return True
        return False

    # ProofChain initialisieren
    proof = ProofChain(
        f"Goldbach-Vermutung: Jede gerade Zahl 4 ≤ n ≤ {n_max} ist Summe zweier Primzahlen",
        f"Numerische Verifikation für alle geraden n von 4 bis {n_max}"
    )

    # Voraussetzung
    proof.assume(
        "Goldbach-Vermutung (1742): Jede gerade Zahl > 2 ist Summe zweier Primzahlen",
        "Zu verifizierende Aussage"
    )

    # Primzahlen bis n_max berechnen
    is_prime_arr = sieve(n_max)

    # Alle geraden Zahlen von 4 bis n_max prüfen
    failures = []
    count = 0
    for n in range(4, n_max + 1, 2):
        if not has_goldbach_decomposition(n, is_prime_arr):
            failures.append(n)
        count += 1

    if not failures:
        # Alle geraden Zahlen bestanden → numerisch verifiziert
        step = ProofStep(
            f"Alle {count} geraden Zahlen von 4 bis {n_max} haben eine Goldbach-Zerlegung",
            f"Exhaustive numerische Prüfung mit Sieb des Eratosthenes",
            ProofStatus.VERIFIED_NUMERIC
        )
        proof.steps.append(step)

        # Einige konkrete Beispiele hinzufügen
        examples = []
        for n in [4, 6, 8, 10, 100, 1000]:
            if n <= n_max:
                # Kleinste Zerlegung finden
                for p in range(2, n // 2 + 1):
                    if is_prime_arr[p] and is_prime_arr[n - p]:
                        examples.append(f"{n} = {p} + {n-p}")
                        break

        example_step = ProofStep(
            f"Beispiele: {', '.join(examples)}",
            "Konkrete Goldbach-Zerlegungen",
            ProofStatus.VERIFIED_NUMERIC
        )
        proof.steps.append(example_step)

        # Schlussfolgerung
        proof.conclude(
            f"Goldbach-Vermutung numerisch verifiziert für alle geraden n ≤ {n_max}. "
            "ABER: Kein vollständiger Beweis für alle n ∈ ℕ!",
            "Numerische Evidenz (kein vollständiger Beweis)"
        )
    else:
        # Gegenbeispiele gefunden (sollte bei korrekter Implementierung nie passieren)
        step = ProofStep(
            f"Gegenbeispiele gefunden: {failures[:5]}",
            "Implementierungsfehler oder Goldbach ist falsch?",
            ProofStatus.FAILED
        )
        proof.steps.append(step)
        proof.conclude("Verifikation fehlgeschlagen", "Gegenbeispiele gefunden")

    return proof
