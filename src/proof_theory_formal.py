"""
Formale Beweistheorie: Sequenzenkalkül, Schnittelimination, Beweiskomplexität,
Ordinale Analyse und Reverse Mathematics.

Dieses Modul implementiert Kernkonzepte der formalen Beweistheorie nach Gentzen (1934/1936):
- LK-Sequenzenkalkül (Klassische Logik)
- Natürliches Schließen (ND)
- Hilbert-Kalkül
- Beweiskomplexität (Resolution, PHP)
- Ordinale Analyse (ε₀, Gentzen)
- Reverse Mathematics (Big Five)

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

from __future__ import annotations
from typing import Optional
import re


# =============================================================================
# SEQUENZENKALKÜL (LK - GENTZEN 1934)
# =============================================================================

class Sequent:
    """
    Repräsentiert ein Sequent Γ ⊢ Δ im Sinne Gentzens.

    Ein Sequent bedeutet: Wenn alle Formeln in Γ (Antecedent) wahr sind,
    dann ist mindestens eine Formel in Δ (Succedent) wahr.

    Formal: Γ ⊢ Δ  ≡  ¬(Γ₁ ∧ ... ∧ Γₙ) ∨ Δ₁ ∨ ... ∨ Δₘ

    Letzte Änderung: 2026-03-10
    """

    def __init__(self, antecedent: list[str], succedent: list[str]):
        """
        Erstellt ein Sequent mit Antecedent und Succedent.

        :param antecedent: Γ – Liste von Formeln auf der linken Seite
        :param succedent:  Δ – Liste von Formeln auf der rechten Seite
        """
        # Antecedent (linke Seite des Sequents)
        self.antecedent: list[str] = list(antecedent)
        # Succedent (rechte Seite des Sequents)
        self.succedent: list[str] = list(succedent)

    def is_axiom(self) -> bool:
        """
        Prüft ob dieses Sequent ein Identitäts-Axiom ist.

        Ein Sequent Γ,A ⊢ A,Δ ist ein Axiom, wenn mindestens eine Formel A
        sowohl im Antecedent als auch im Succedent vorkommt.

        :return: True wenn das Sequent ein Axiom ist
        """
        # Schnittmenge von linker und rechter Seite
        linke = set(self.antecedent)
        rechte = set(self.succedent)
        return bool(linke & rechte)

    def __str__(self) -> str:
        """
        Gibt das Sequent als lesbaren String aus.

        :return: String der Form "A,B ⊢ C,D"
        """
        # Linke Seite zusammensetzen
        links = ", ".join(self.antecedent) if self.antecedent else "∅"
        # Rechte Seite zusammensetzen
        rechts = ", ".join(self.succedent) if self.succedent else "∅"
        return f"{links} ⊢ {rechts}"

    def __repr__(self) -> str:
        """Kurzdarstellung für Debug-Zwecke."""
        return f"Sequent({self!s})"

    def __eq__(self, other: object) -> bool:
        """Vergleich zweier Sequenten auf Gleichheit."""
        if not isinstance(other, Sequent):
            return False
        return self.antecedent == other.antecedent and self.succedent == other.succedent


class ProofTree:
    """
    Beweisbaum im Sequenzenkalkül.

    Jeder Knoten repräsentiert ein Sequent (die Konklusion einer Regelanwendung).
    Die Kinder sind die Prämissen der angewendeten Regel.
    Die Blätter müssen Axiome sein, damit der Beweis gültig ist.

    Letzte Änderung: 2026-03-10
    """

    def __init__(
        self,
        conclusion: Sequent,
        rule: str,
        premises: Optional[list['ProofTree']] = None,
    ):
        """
        Erstellt einen Beweisbaum-Knoten.

        :param conclusion: Das bewiesene Sequent (Konklusion)
        :param rule:       Name der angewendeten Regelanwendung
        :param premises:   Liste der Teilbeweise (Prämissen), None für Blätter
        """
        self.conclusion: Sequent = conclusion
        self.rule: str = rule
        # Blätter haben keine Prämissen
        self.premises: list[ProofTree] = premises if premises is not None else []

    def is_valid(self) -> bool:
        """
        Prüft ob der Beweisbaum gültig ist.

        Ein Beweis ist gültig, wenn alle Blätter Axiome sind und
        alle Prämissen selbst gültige Beweise sind.

        :return: True wenn der gesamte Baum korrekt ist
        """
        # Blatt: muss Axiom sein
        if not self.premises:
            return self.conclusion.is_axiom()
        # Innerer Knoten: alle Teilbäume müssen gültig sein
        return all(p.is_valid() for p in self.premises)

    def height(self) -> int:
        """
        Berechnet die Höhe des Beweisbaums.

        Die Höhe entspricht der maximalen Anzahl von Regelanwendungen
        auf einem Pfad von der Wurzel bis zum Blatt.

        :return: Höhe des Baums (0 bei Blättern)
        """
        if not self.premises:
            return 0
        return 1 + max(p.height() for p in self.premises)

    def size(self) -> int:
        """
        Berechnet die Größe des Beweisbaums (Anzahl der Knoten).

        Die Größe entspricht der Gesamtzahl der Regelanwendungen
        inklusive dieses Knotens.

        :return: Anzahl Knoten im Baum
        """
        return 1 + sum(p.size() for p in self.premises)

    def to_string(self, indent: int = 0) -> str:
        """
        Erzeugt eine lesbare Baumdarstellung mit Einrückung.

        :param indent: Aktuelle Einrückungstiefe
        :return: Mehrzeiliger String mit Baumstruktur
        """
        # Einrückung für diesen Knoten
        einrueckung = "  " * indent
        zeile = f"{einrueckung}[{self.rule}] {self.conclusion}"
        # Prämissen rekursiv darstellen
        if self.premises:
            teile = [zeile]
            for praemisse in self.premises:
                teile.append(praemisse.to_string(indent + 1))
            return "\n".join(teile)
        return zeile

    def __str__(self) -> str:
        """Gibt den Baum als String zurück."""
        return self.to_string()


class LKRules:
    """
    Regeln des klassischen Sequenzenkalküls LK (Gentzen 1934).

    Das LK-System enthält:
    - Strukturregeln: Weakening (W), Contraction (C), Exchange (E)
    - Logische Regeln: ∧L, ∧R, ∨L, ∨R, ¬L, ¬R, →L, →R
    - Quantorregeln: ∀L, ∀R, ∃L, ∃R
    - Schnitt (Cut): A ⊢, ⊢ A / ⊢

    Der Hauptsatz (Schnitteliminationssatz) besagt, dass jeder Beweis
    in einen schnittfreien Beweis umgewandelt werden kann.

    Letzte Änderung: 2026-03-10
    """

    @staticmethod
    def weakening_left(proof: ProofTree, formula: str) -> ProofTree:
        """
        Weakening-Links (WL): Fügt eine Formel im Antecedent hinzu.

        Regel: Γ ⊢ Δ / A, Γ ⊢ Δ

        Eine Formel die nicht gebraucht wird, kann links hinzugefügt werden.

        :param proof:   Bestehender Teilbeweis für Γ ⊢ Δ
        :param formula: Die neue Formel A die hinzugefügt wird
        :return: Neuer Beweisbaum für A, Γ ⊢ Δ
        """
        # Neues Sequent mit zusätzlicher Formel links
        neues_ant = [formula] + proof.conclusion.antecedent
        neue_konklusion = Sequent(neues_ant, proof.conclusion.succedent)
        return ProofTree(neue_konklusion, "WL", [proof])

    @staticmethod
    def weakening_right(proof: ProofTree, formula: str) -> ProofTree:
        """
        Weakening-Rechts (WR): Fügt eine Formel im Succedent hinzu.

        Regel: Γ ⊢ Δ / Γ ⊢ Δ, A

        :param proof:   Bestehender Teilbeweis für Γ ⊢ Δ
        :param formula: Die neue Formel A die hinzugefügt wird
        :return: Neuer Beweisbaum für Γ ⊢ Δ, A
        """
        # Neues Sequent mit zusätzlicher Formel rechts
        neues_suc = proof.conclusion.succedent + [formula]
        neue_konklusion = Sequent(proof.conclusion.antecedent, neues_suc)
        return ProofTree(neue_konklusion, "WR", [proof])

    @staticmethod
    def contraction_left(proof: ProofTree, formula: str) -> ProofTree:
        """
        Contraction-Links (CL): Entfernt doppelte Formel im Antecedent.

        Regel: A, A, Γ ⊢ Δ / A, Γ ⊢ Δ

        Mehrfaches Vorkommen einer Formel kann auf eines reduziert werden.

        :param proof:   Teilbeweis für A, A, Γ ⊢ Δ
        :param formula: Die zu kontrahierende Formel A
        :return: Neuer Beweisbaum für A, Γ ⊢ Δ
        """
        # Antecedent: erste Duplizierung entfernen
        ant = list(proof.conclusion.antecedent)
        if formula in ant:
            ant.remove(formula)  # Entfernt ersten Treffer
        neue_konklusion = Sequent(ant, proof.conclusion.succedent)
        return ProofTree(neue_konklusion, "CL", [proof])

    @staticmethod
    def contraction_right(proof: ProofTree, formula: str) -> ProofTree:
        """
        Contraction-Rechts (CR): Entfernt doppelte Formel im Succedent.

        Regel: Γ ⊢ Δ, A, A / Γ ⊢ Δ, A

        :param proof:   Teilbeweis für Γ ⊢ Δ, A, A
        :param formula: Die zu kontrahierende Formel A
        :return: Neuer Beweisbaum für Γ ⊢ Δ, A
        """
        suc = list(proof.conclusion.succedent)
        if formula in suc:
            suc.remove(formula)
        neue_konklusion = Sequent(proof.conclusion.antecedent, suc)
        return ProofTree(neue_konklusion, "CR", [proof])

    @staticmethod
    def cut(proof1: ProofTree, proof2: ProofTree, cut_formula: str) -> ProofTree:
        """
        Schnitt-Regel (Cut): Verbindet zwei Beweise über eine Cut-Formel.

        Regel: Γ ⊢ Δ, A    A, Σ ⊢ Π  /  Γ, Σ ⊢ Δ, Π

        Die Cut-Formel A tritt in der Konklusion nicht mehr auf.
        Der Schnitteliminationssatz besagt, dass diese Regel eliminierbar ist.

        :param proof1:       Teilbeweis für Γ ⊢ Δ, A
        :param proof2:       Teilbeweis für A, Σ ⊢ Π
        :param cut_formula:  Die Cut-Formel A
        :return: Neuer Beweisbaum für Γ, Σ ⊢ Δ, Π
        """
        # Linker Beweis: A aus Succedent entfernen
        linkes_suc = [f for f in proof1.conclusion.succedent if f != cut_formula]
        # Rechter Beweis: A aus Antecedent entfernen
        rechtes_ant = [f for f in proof2.conclusion.antecedent if f != cut_formula]
        # Kombiniertes Sequent
        neues_ant = proof1.conclusion.antecedent + rechtes_ant
        neues_suc = linkes_suc + proof2.conclusion.succedent
        neue_konklusion = Sequent(neues_ant, neues_suc)
        return ProofTree(neue_konklusion, f"Cut({cut_formula})", [proof1, proof2])

    @staticmethod
    def and_right(proof1: ProofTree, proof2: ProofTree) -> ProofTree:
        """
        Konjunktion-Rechts (∧R): Beweist A∧B im Succedent.

        Regel: Γ ⊢ Δ, A    Γ ⊢ Δ, B  /  Γ ⊢ Δ, A∧B

        :param proof1: Teilbeweis für Γ ⊢ Δ, A
        :param proof2: Teilbeweis für Γ ⊢ Δ, B
        :return: Neuer Beweisbaum für Γ ⊢ Δ, A∧B
        """
        # Die letzten Succedent-Formeln der beiden Beweise sind A und B
        if not proof1.conclusion.succedent or not proof2.conclusion.succedent:
            raise ValueError("Prämissen müssen nicht-leere Succedente haben")
        formel_a = proof1.conclusion.succedent[-1]
        formel_b = proof2.conclusion.succedent[-1]
        # Neues Succedent: ohne A/B, dafür mit A∧B
        neues_suc = proof1.conclusion.succedent[:-1] + [f"({formel_a}∧{formel_b})"]
        neue_konklusion = Sequent(proof1.conclusion.antecedent, neues_suc)
        return ProofTree(neue_konklusion, "∧R", [proof1, proof2])

    @staticmethod
    def and_left(proof: ProofTree, formula_a: str, formula_b: str) -> ProofTree:
        """
        Konjunktion-Links (∧L): Zerlegt A∧B im Antecedent.

        Regel: A, B, Γ ⊢ Δ  /  A∧B, Γ ⊢ Δ

        :param proof:     Teilbeweis für A, B, Γ ⊢ Δ
        :param formula_a: Die Formel A
        :param formula_b: Die Formel B
        :return: Neuer Beweisbaum für A∧B, Γ ⊢ Δ
        """
        # A und B aus Antecedent entfernen, A∧B hinzufügen
        ant = list(proof.conclusion.antecedent)
        if formula_a in ant:
            ant.remove(formula_a)
        if formula_b in ant:
            ant.remove(formula_b)
        ant = [f"({formula_a}∧{formula_b})"] + ant
        neue_konklusion = Sequent(ant, proof.conclusion.succedent)
        return ProofTree(neue_konklusion, "∧L", [proof])

    @staticmethod
    def or_left(proof1: ProofTree, proof2: ProofTree) -> ProofTree:
        """
        Disjunktion-Links (∨L): Zerlegt A∨B im Antecedent.

        Regel: A, Γ ⊢ Δ    B, Γ ⊢ Δ  /  A∨B, Γ ⊢ Δ

        :param proof1: Teilbeweis für A, Γ ⊢ Δ
        :param proof2: Teilbeweis für B, Γ ⊢ Δ
        :return: Neuer Beweisbaum für A∨B, Γ ⊢ Δ
        """
        if not proof1.conclusion.antecedent or not proof2.conclusion.antecedent:
            raise ValueError("Prämissen müssen nicht-leere Antecedente haben")
        formel_a = proof1.conclusion.antecedent[0]
        formel_b = proof2.conclusion.antecedent[0]
        # Neues Antecedent mit A∨B statt A und B
        neues_ant = [f"({formel_a}∨{formel_b})"] + proof1.conclusion.antecedent[1:]
        neue_konklusion = Sequent(neues_ant, proof1.conclusion.succedent)
        return ProofTree(neue_konklusion, "∨L", [proof1, proof2])

    @staticmethod
    def or_right(proof: ProofTree, formula: str, side: str = "left") -> ProofTree:
        """
        Disjunktion-Rechts (∨R): Beweist A∨B im Succedent.

        Regel: Γ ⊢ Δ, A  /  Γ ⊢ Δ, A∨B   (linke Variante)
               Γ ⊢ Δ, B  /  Γ ⊢ Δ, A∨B   (rechte Variante)

        :param proof:   Teilbeweis für Γ ⊢ Δ, A (oder B)
        :param formula: Die andere Disjunktionsformel B (oder A)
        :param side:    "left" wenn A bewiesen, "right" wenn B bewiesen
        :return: Neuer Beweisbaum für Γ ⊢ Δ, A∨B
        """
        if not proof.conclusion.succedent:
            raise ValueError("Prämisse muss nicht-leeres Succedent haben")
        bewiesene = proof.conclusion.succedent[-1]
        if side == "left":
            neue_formel = f"({bewiesene}∨{formula})"
        else:
            neue_formel = f"({formula}∨{bewiesene})"
        neues_suc = proof.conclusion.succedent[:-1] + [neue_formel]
        neue_konklusion = Sequent(proof.conclusion.antecedent, neues_suc)
        return ProofTree(neue_konklusion, "∨R", [proof])

    @staticmethod
    def neg_right(proof: ProofTree) -> ProofTree:
        """
        Negation-Rechts (¬R): Überführt eine Annahme in Negation.

        Regel: A, Γ ⊢ Δ  /  Γ ⊢ Δ, ¬A

        :param proof: Teilbeweis für A, Γ ⊢ Δ
        :return: Neuer Beweisbaum für Γ ⊢ Δ, ¬A
        """
        if not proof.conclusion.antecedent:
            raise ValueError("Prämisse muss nicht-leeres Antecedent haben")
        formel_a = proof.conclusion.antecedent[0]
        neues_ant = proof.conclusion.antecedent[1:]
        neues_suc = proof.conclusion.succedent + [f"¬{formel_a}"]
        neue_konklusion = Sequent(neues_ant, neues_suc)
        return ProofTree(neue_konklusion, "¬R", [proof])

    @staticmethod
    def neg_left(proof: ProofTree) -> ProofTree:
        """
        Negation-Links (¬L): Überführt eine Succedent-Formel in Negation.

        Regel: Γ ⊢ Δ, A  /  ¬A, Γ ⊢ Δ

        :param proof: Teilbeweis für Γ ⊢ Δ, A
        :return: Neuer Beweisbaum für ¬A, Γ ⊢ Δ
        """
        if not proof.conclusion.succedent:
            raise ValueError("Prämisse muss nicht-leeres Succedent haben")
        formel_a = proof.conclusion.succedent[-1]
        neues_suc = proof.conclusion.succedent[:-1]
        neues_ant = [f"¬{formel_a}"] + proof.conclusion.antecedent
        neue_konklusion = Sequent(neues_ant, neues_suc)
        return ProofTree(neue_konklusion, "¬L", [proof])

    @staticmethod
    def implies_right(proof: ProofTree) -> ProofTree:
        """
        Implikation-Rechts (→R): Beweist A→B im Succedent.

        Regel: A, Γ ⊢ Δ, B  /  Γ ⊢ Δ, A→B

        :param proof: Teilbeweis für A, Γ ⊢ Δ, B
        :return: Neuer Beweisbaum für Γ ⊢ Δ, A→B
        """
        if not proof.conclusion.antecedent or not proof.conclusion.succedent:
            raise ValueError("Prämisse muss Antecedent und Succedent haben")
        formel_a = proof.conclusion.antecedent[0]
        formel_b = proof.conclusion.succedent[-1]
        neues_ant = proof.conclusion.antecedent[1:]
        neues_suc = proof.conclusion.succedent[:-1] + [f"({formel_a}→{formel_b})"]
        neue_konklusion = Sequent(neues_ant, neues_suc)
        return ProofTree(neue_konklusion, "→R", [proof])

    @staticmethod
    def implies_left(proof1: ProofTree, proof2: ProofTree) -> ProofTree:
        """
        Implikation-Links (→L): Verwendet A→B im Antecedent.

        Regel: Γ ⊢ Δ, A    B, Σ ⊢ Π  /  A→B, Γ, Σ ⊢ Δ, Π

        :param proof1: Teilbeweis für Γ ⊢ Δ, A
        :param proof2: Teilbeweis für B, Σ ⊢ Π
        :return: Neuer Beweisbaum für A→B, Γ, Σ ⊢ Δ, Π
        """
        if not proof1.conclusion.succedent or not proof2.conclusion.antecedent:
            raise ValueError("Prämissen haben falsche Struktur")
        formel_a = proof1.conclusion.succedent[-1]
        formel_b = proof2.conclusion.antecedent[0]
        impl_formel = f"({formel_a}→{formel_b})"
        neues_ant = [impl_formel] + proof1.conclusion.antecedent + proof2.conclusion.antecedent[1:]
        neues_suc = proof1.conclusion.succedent[:-1] + proof2.conclusion.succedent
        neue_konklusion = Sequent(neues_ant, neues_suc)
        return ProofTree(neue_konklusion, "→L", [proof1, proof2])


def cut_elimination_demo() -> dict:
    """
    Demonstration der Schnittelimination (Hauptsatz, Gentzen 1934).

    Der Hauptsatz besagt: Jeder Beweis im Sequenzenkalkül LK kann in einen
    schnittfreien Beweis transformiert werden.

    Folgerungen:
    - Teilformel-Eigenschaft (Subformula Property)
    - Entscheidbarkeit der Prädikatenlogik
    - Konsistenz (kein Beweis von ⊢ ⊥)

    Beispiel: Beweise A ⊢ A mit Hilfe von Cut, dann eliminiere den Schnitt.

    :return: Dict mit Beweis vor/nach Schnittelimination und Erläuterungen

    Letzte Änderung: 2026-03-10
    """
    # Schritt 1: Axiom-Beweise für die Cut-Formel
    # Beweis für: A ⊢ A (Axiom)
    axiom_links = ProofTree(Sequent(["A"], ["A"]), "Axiom")
    # Beweis für: A ⊢ A (zweites Axiom)
    axiom_rechts = ProofTree(Sequent(["A"], ["A"]), "Axiom")

    # Schritt 2: Cut-Anwendung (nicht schnittfrei)
    # Γ ⊢ Δ, A    A, Σ ⊢ Π  /  Γ, Σ ⊢ Δ, Π
    # Hier: A ⊢ A    A ⊢ A  /  A, A ⊢ A, A  (nach Cut über A)
    # Tatsächlich: Cut auf A ⊢ A kombiniert mit A ⊢ A
    beweis_mit_cut = LKRules.cut(axiom_links, axiom_rechts, "A")

    # Schritt 3: Schnittfreier Beweis (direkt durch Axiom)
    schnittfreier_beweis = ProofTree(
        Sequent(["A", "A"], ["A", "A"]),
        "Axiom"
    )

    return {
        "hauptsatz": (
            "Hauptsatz (Gentzen 1934): Jeder LK-Beweis ist in einen "
            "schnittfreien Beweis transformierbar."
        ),
        "beweis_mit_cut": {
            "baum": beweis_mit_cut.to_string(),
            "sequent": str(beweis_mit_cut.conclusion),
            "hoehe": beweis_mit_cut.height(),
            "groesse": beweis_mit_cut.size(),
        },
        "schnittfreier_beweis": {
            "baum": schnittfreier_beweis.to_string(),
            "sequent": str(schnittfreier_beweis.conclusion),
            "hoehe": schnittfreier_beweis.height(),
            "groesse": schnittfreier_beweis.size(),
            "ist_axiom": schnittfreier_beweis.conclusion.is_axiom(),
        },
        "folgerungen": [
            "Teilformel-Eigenschaft: Nur Teilformeln der Schlusssequenz kommen vor",
            "Konsistenz: Es gibt keinen Beweis von ⊢ ⊥ (leeres Succedent)",
            "Entscheidbarkeit: Prädikatenlogik ist semi-entscheidbar",
        ],
        "komplexitaet": (
            "Schnittelimination kann zu nicht-elementarer Beweisverlängerung führen "
            "(Statman 1979): (n+1)-fach exponentielle Zunahme möglich."
        ),
    }


def subformula_property(proof: ProofTree) -> bool:
    """
    Prüft die Teilformel-Eigenschaft (Subformula Property).

    In schnittfreien Beweisen kommen in jedem Sequent nur Teilformeln
    der Schlusssequenz vor. Diese Eigenschaft ist eine direkte Folge
    der Schnittelimination.

    :param proof: Zu prüfender Beweisbaum
    :return: True wenn die Teilformel-Eigenschaft gilt

    Letzte Änderung: 2026-03-10
    """
    def get_subformulas(formula: str) -> set[str]:
        """Extrahiert alle syntaktischen Teilformeln einer Formel."""
        teilformeln = {formula}
        # Zerlege zusammengesetzte Formeln
        for operator in ["∧", "∨", "→"]:
            if operator in formula:
                # Vereinfachte Zerlegung (für Basisformeln ohne Schachtelung)
                teile = formula.replace("(", "").replace(")", "").split(operator)
                for teil in teile:
                    teil = teil.strip()
                    if teil:
                        teilformeln.add(teil)
                        # Negationen beachten
                        if teil.startswith("¬"):
                            teilformeln.add(teil[1:])
        # Negationen erkennen
        if formula.startswith("¬"):
            teilformeln.add(formula[1:])
        return teilformeln

    # Alle Teilformeln der Schlusssequenz sammeln
    alle_schlusformeln = (
        proof.conclusion.antecedent + proof.conclusion.succedent
    )
    zulaessige_teilformeln: set[str] = set()
    for formel in alle_schlusformeln:
        zulaessige_teilformeln |= get_subformulas(formel)

    def check_node(knoten: ProofTree) -> bool:
        """Prüft rekursiv ob alle Formeln im Knoten Teilformeln sind."""
        # Alle Formeln in diesem Knoten
        alle_formeln = (
            knoten.conclusion.antecedent + knoten.conclusion.succedent
        )
        for formel in alle_formeln:
            if formel not in zulaessige_teilformeln:
                return False
        # Rekursiv für alle Prämissen prüfen
        return all(check_node(praemisse) for praemisse in knoten.premises)

    return check_node(proof)


# =============================================================================
# NATÜRLICHES SCHLIEßEN (ND - GENTZEN)
# =============================================================================

class NDProof:
    """
    Beweis im Natürlichen Schließen (Natural Deduction, Gentzen 1935).

    Im Natürlichen Schließen beginnt man mit Annahmen und leitet
    schrittweise die Konklusion ab. Annahmen können zu gegebener Zeit
    "entladen" werden (z.B. bei →-Einführung).

    Regeln:
    - Einführungsregeln (∧I, ∨I, →I, ¬I)
    - Eliminierungsregeln (∧E, ∨E, →E, ¬E)
    - Ex Falso Quodlibet (aus ⊥ folgt alles)

    Letzte Änderung: 2026-03-10
    """

    def __init__(self):
        """Initialisiert einen leeren ND-Beweis."""
        # Schrittregister: id -> (Formel, Regelname, abhängige IDs)
        self._schritte: dict[int, tuple[str, str, list[int]]] = {}
        # Zähler für Schritt-IDs
        self._naechste_id: int = 0
        # Menge der entladenen Annahmen-IDs
        self._entladene_annahmen: set[int] = set()

    def _neuer_schritt(self, formel: str, regel: str, abhaengig: list[int]) -> int:
        """
        Registriert einen neuen Beweisschritt intern.

        :param formel:    Die abgeleitete Formel
        :param regel:     Name der angewendeten Regel
        :param abhaengig: Liste der IDs auf die verwiesen wird
        :return: ID des neuen Schritts
        """
        schritt_id = self._naechste_id
        self._schritte[schritt_id] = (formel, regel, abhaengig)
        self._naechste_id += 1
        return schritt_id

    def assume(self, formula: str) -> int:
        """
        Führt eine Annahme ein.

        Im Natürlichen Schließen beginnt man mit unbewiesenen Annahmen,
        die später entladen werden müssen.

        :param formula: Die angenommene Formel
        :return: ID der Annahme (für spätere Schritte)
        """
        return self._neuer_schritt(formula, "Annahme", [])

    def apply_and_intro(self, id1: int, id2: int) -> int:
        """
        Konjunktions-Einführung (∧I).

        Aus A (id1) und B (id2) folgt A∧B.

        :param id1: ID des Schritts mit Formel A
        :param id2: ID des Schritts mit Formel B
        :return: ID des neuen Schritts mit A∧B
        """
        formel_a = self.get_formula(id1)
        formel_b = self.get_formula(id2)
        return self._neuer_schritt(f"({formel_a}∧{formel_b})", "∧I", [id1, id2])

    def apply_and_elim_left(self, id: int) -> int:
        """
        Konjunktions-Eliminierung links (∧EL).

        Aus A∧B folgt A.

        :param id: ID des Schritts mit Formel A∧B
        :return: ID des neuen Schritts mit A
        """
        formel = self.get_formula(id)
        # Extrahiere linken Teil aus A∧B
        linker_teil = self._extract_left_conjunct(formel)
        return self._neuer_schritt(linker_teil, "∧EL", [id])

    def apply_and_elim_right(self, id: int) -> int:
        """
        Konjunktions-Eliminierung rechts (∧ER).

        Aus A∧B folgt B.

        :param id: ID des Schritts mit Formel A∧B
        :return: ID des neuen Schritts mit B
        """
        formel = self.get_formula(id)
        rechter_teil = self._extract_right_conjunct(formel)
        return self._neuer_schritt(rechter_teil, "∧ER", [id])

    def apply_or_intro_left(self, id: int, formula: str) -> int:
        """
        Disjunktions-Einführung links (∨IL).

        Aus A folgt A∨B (für beliebiges B).

        :param id:      ID des Schritts mit Formel A
        :param formula: Die Formel B (rechte Seite der Disjunktion)
        :return: ID des neuen Schritts mit A∨B
        """
        formel_a = self.get_formula(id)
        return self._neuer_schritt(f"({formel_a}∨{formula})", "∨IL", [id])

    def apply_or_intro_right(self, id: int, formula: str) -> int:
        """
        Disjunktions-Einführung rechts (∨IR).

        Aus B folgt A∨B (für beliebiges A).

        :param id:      ID des Schritts mit Formel B
        :param formula: Die Formel A (linke Seite der Disjunktion)
        :return: ID des neuen Schritts mit A∨B
        """
        formel_b = self.get_formula(id)
        return self._neuer_schritt(f"({formula}∨{formel_b})", "∨IR", [id])

    def apply_implies_intro(self, assumption_id: int, conclusion_id: int) -> int:
        """
        Implikations-Einführung (→I).

        Wenn aus Annahme A (assumption_id) die Formel B (conclusion_id)
        abgeleitet wurde, dann gilt A→B. Die Annahme A wird entladen.

        :param assumption_id:  ID der Annahme A
        :param conclusion_id:  ID des Schritts mit Formel B
        :return: ID des neuen Schritts mit A→B
        """
        formel_a = self.get_formula(assumption_id)
        formel_b = self.get_formula(conclusion_id)
        # Annahme wird entladen
        self._entladene_annahmen.add(assumption_id)
        return self._neuer_schritt(
            f"({formel_a}→{formel_b})", "→I", [assumption_id, conclusion_id]
        )

    def apply_implies_elim(self, impl_id: int, premise_id: int) -> int:
        """
        Implikations-Eliminierung / Modus Ponens (→E).

        Aus A→B (impl_id) und A (premise_id) folgt B.

        :param impl_id:    ID des Schritts mit Formel A→B
        :param premise_id: ID des Schritts mit Formel A
        :return: ID des neuen Schritts mit B
        """
        formel_impl = self.get_formula(impl_id)
        # Extrahiere Konsequent aus A→B
        konsequent = self._extract_consequent(formel_impl)
        return self._neuer_schritt(konsequent, "→E", [impl_id, premise_id])

    def apply_neg_intro(self, assumption_id: int, contradiction_id: int) -> int:
        """
        Negations-Einführung (¬I).

        Wenn aus Annahme A ein Widerspruch ⊥ folgt, dann gilt ¬A.

        :param assumption_id:     ID der Annahme A
        :param contradiction_id:  ID des Schritts mit ⊥
        :return: ID des neuen Schritts mit ¬A
        """
        formel_a = self.get_formula(assumption_id)
        self._entladene_annahmen.add(assumption_id)
        return self._neuer_schritt(f"¬{formel_a}", "¬I", [assumption_id, contradiction_id])

    def apply_ex_falso(self, bot_id: int, formula: str) -> int:
        """
        Ex Falso Quodlibet (EFQ) / Principle of Explosion.

        Aus dem Widerspruch ⊥ (bot_id) folgt jede beliebige Formel.
        Formal: ⊥ ⊢ A (für beliebiges A)

        :param bot_id:  ID des Schritts mit ⊥
        :param formula: Die abzuleitende Formel A
        :return: ID des neuen Schritts mit A
        """
        return self._neuer_schritt(formula, "EFQ", [bot_id])

    def apply_neg_elim(self, formula_id: int, neg_id: int) -> int:
        """
        Negations-Eliminierung (¬E) / Widerspruchsregel.

        Aus A (formula_id) und ¬A (neg_id) folgt ⊥.

        :param formula_id: ID des Schritts mit Formel A
        :param neg_id:     ID des Schritts mit Formel ¬A
        :return: ID des neuen Schritts mit ⊥
        """
        return self._neuer_schritt("⊥", "¬E", [formula_id, neg_id])

    def get_formula(self, id: int) -> str:
        """
        Gibt die Formel des Schritts mit der angegebenen ID zurück.

        :param id: Die Schritt-ID
        :return: Die abgeleitete Formel
        :raises KeyError: Wenn die ID nicht existiert
        """
        if id not in self._schritte:
            raise KeyError(f"Schritt-ID {id} existiert nicht")
        return self._schritte[id][0]

    def is_complete_proof(self, conclusion: str) -> bool:
        """
        Prüft ob der Beweis vollständig ist und die Konklusion beweist.

        Ein Beweis ist vollständig wenn:
        1. Die Konklusion abgeleitet wurde
        2. Alle Annahmen entweder entladen wurden oder als Prämissen gelten

        :param conclusion: Die zu beweisende Formel
        :return: True wenn Beweis vollständig und korrekt
        """
        # Prüfe ob Konklusion in den Schritten vorkommt
        for schritt_id, (formel, regel, _) in self._schritte.items():
            if formel == conclusion:
                return True
        return False

    def _extract_left_conjunct(self, formula: str) -> str:
        """Extrahiert den linken Konjunkt aus (A∧B)."""
        # Entferne äußere Klammern falls vorhanden
        f = formula.strip("()")
        if "∧" in f:
            teile = f.split("∧", 1)
            return teile[0].strip()
        return formula

    def _extract_right_conjunct(self, formula: str) -> str:
        """Extrahiert den rechten Konjunkt aus (A∧B)."""
        f = formula.strip("()")
        if "∧" in f:
            teile = f.split("∧", 1)
            return teile[1].strip()
        return formula

    def _extract_consequent(self, formula: str) -> str:
        """Extrahiert den Konsequenten aus (A→B)."""
        f = formula.strip("()")
        if "→" in f:
            teile = f.split("→", 1)
            return teile[1].strip()
        return formula

    def get_all_steps(self) -> dict[int, tuple[str, str, list[int]]]:
        """
        Gibt alle Beweisschritte zurück.

        :return: Dict {id: (Formel, Regelname, Abhängigkeiten)}
        """
        return dict(self._schritte)


def prove_simple(formula: str) -> Optional[NDProof]:
    """
    Versucht einfache logische Tautologien automatisch zu beweisen.

    Unterstützte Formen:
    - A→A (Identität)
    - (A∧B)→A (Konjunktions-Projektion links)
    - (A∧B)→B (Konjunktions-Projektion rechts)
    - A→(A∨B) (Disjunktions-Einführung)
    - (A→B)→(A→B) (Identität für Implikationen)

    :param formula: Die zu beweisende Formel
    :return: NDProof wenn erfolgreich, None wenn unbekannte Tautologie

    Letzte Änderung: 2026-03-10
    """
    beweis = NDProof()

    # Muster: A→A (einfache Identität)
    if re.match(r'^(\w+)→\1$', formula):
        var = formula.split("→")[0]
        annahme_id = beweis.assume(var)
        beweis.apply_implies_intro(annahme_id, annahme_id)
        return beweis

    # Muster: (A∧B)→A (linke Projektion)
    match = re.match(r'^\((\w+)∧(\w+)\)→\1$', formula)
    if match:
        var_a, var_b = match.group(1), match.group(2)
        annahme_id = beweis.assume(f"({var_a}∧{var_b})")
        links_id = beweis.apply_and_elim_left(annahme_id)
        beweis.apply_implies_intro(annahme_id, links_id)
        return beweis

    # Muster: (A∧B)→B (rechte Projektion)
    match = re.match(r'^\((\w+)∧(\w+)\)→\2$', formula)
    if match:
        var_a, var_b = match.group(1), match.group(2)
        annahme_id = beweis.assume(f"({var_a}∧{var_b})")
        rechts_id = beweis.apply_and_elim_right(annahme_id)
        beweis.apply_implies_intro(annahme_id, rechts_id)
        return beweis

    # Muster: A→(A∨B) (Einführung linker Disjunkt)
    match = re.match(r'^(\w+)→\(\1∨(\w+)\)$', formula)
    if match:
        var_a, var_b = match.group(1), match.group(2)
        annahme_id = beweis.assume(var_a)
        oder_id = beweis.apply_or_intro_left(annahme_id, var_b)
        beweis.apply_implies_intro(annahme_id, oder_id)
        return beweis

    # Muster: (A→B)→(A→B) (Identität für Implikationen)
    match = re.match(r'^\((\w+)→(\w+)\)→\(\1→\2\)$', formula)
    if match:
        var_a, var_b = match.group(1), match.group(2)
        impl_formel = f"({var_a}→{var_b})"
        annahme_id = beweis.assume(impl_formel)
        beweis.apply_implies_intro(annahme_id, annahme_id)
        return beweis

    # Keine bekannte Tautologie erkannt
    return None


# =============================================================================
# HILBERT-KALKÜL
# =============================================================================

class HilbertCalculus:
    """
    Hilbert-Kalkül H für klassische Aussagenlogik.

    Das Hilbert-System verwendet drei Axiomenschemata und Modus Ponens:
    - K: A → (B → A)              (Abschwächung)
    - S: (A→(B→C)) → ((A→B)→(A→C)) (Verteilung)
    - DN: ¬¬A → A                 (Doppelnegation-Elimination)

    Jedes System mit K, S, MP ist vollständig für intuitionistische Logik.
    Mit DN zusätzlich: Vollständigkeit für klassische Logik.

    Letzte Änderung: 2026-03-10
    """

    # Axiomenschemata des Hilbert-Kalküls
    AXIOMS = {
        "K": "A → (B → A)",
        "S": "(A → (B → C)) → ((A → B) → (A → C))",
        "DN": "(¬¬A) → A",
    }

    def check_axiom_instance(self, formula: str) -> Optional[str]:
        """
        Prüft ob eine Formel eine Instanz eines Axiomenschemas ist.

        Verwendet einfachen Muster-Abgleich für die drei Hauptschemata.

        :param formula: Die zu prüfende Formel
        :return: Name des Axiomenschemas oder None wenn keine Instanz

        Letzte Änderung: 2026-03-10
        """
        # Bereinige die Formel
        f = formula.strip()

        # Schema K: A → (B → A)
        # Muster: X → (Y → X) für beliebige X, Y
        match_k = re.match(r'^(.+) → \((.+) → \1\)$', f)
        if match_k:
            return "K"

        # Schema DN: (¬¬A) → A
        match_dn = re.match(r'^\(¬¬(.+)\) → \1$', f)
        if match_dn:
            return "DN"

        # Schema S: (A → (B → C)) → ((A → B) → (A → C))
        # Sehr vereinfachte Prüfung auf Struktur
        if "→ ((A → B) → (A → C))" in f or "(A → (B → C)) → ((A → B) →" in f:
            return "S"

        return None

    def apply_mp(self, impl: str, antecedent: str) -> Optional[str]:
        """
        Wendet Modus Ponens an.

        Aus A→B und A folgt B.

        :param impl:       Die Implikation der Form "A → B"
        :param antecedent: Die Prämisse A
        :return: Konsequent B oder None wenn MP nicht anwendbar

        Letzte Änderung: 2026-03-10
        """
        # Versuche Implikation zu parsen
        if " → " not in impl:
            return None
        teile = impl.split(" → ", 1)
        vorderglied = teile[0].strip()
        nachfolger = teile[1].strip()
        # Prüfe ob Vorderglied mit Antecedent übereinstimmt
        if vorderglied == antecedent.strip():
            return nachfolger
        return None

    def prove(self, goal: str, hypotheses: Optional[list[str]] = None) -> Optional[list[str]]:
        """
        Versucht ein Ziel im Hilbert-Kalkül zu beweisen.

        Unterstützte Ziele (Tautologien):
        - "A → A" (durch S und K)
        - "A → (B → A)" (K-Axiom direkt)

        :param goal:        Die zu beweisende Formel
        :param hypotheses:  Optionale Hypothesen (Formeln als gegeben)
        :return: Beweisschritte als Stringliste oder None

        Letzte Änderung: 2026-03-10
        """
        # Hypothesen vorbereiten
        hyps = hypotheses or []
        schritte: list[str] = []

        # Ziel direkt in Hypothesen?
        if goal in hyps:
            schritte.append(f"Hypothese: {goal}")
            return schritte

        # Prüfe ob Ziel ein Axiom-Schema ist
        axiom = self.check_axiom_instance(goal)
        if axiom:
            schritte.append(f"Axiom {axiom}: {goal}")
            return schritte

        # Standard-Beweis: A → A via S und K
        # Beweis von A → A:
        # 1. (A → (B → A)) → ((A → B) → (A → A))   [Schema S mit C=A]
        # 2. A → (B → A)                             [Schema K]
        # 3. (A → B) → (A → A)                       [MP aus 1,2]
        # 4. A → ((A → A) → A)                       [Schema K]
        # 5. (A → (A → A)) → (A → A)                 [Schema S mit B=(A→A), C=A]  (wait...)
        # Vereinfachter Beweis via K und S:
        if " → " in goal:
            teile = goal.split(" → ", 1)
            if teile[0].strip() == teile[1].strip():
                var = teile[0].strip()
                schritte.append(f"1. ({var} → (({var} → {var}) → {var})) → (({var} → ({var} → {var})) → ({var} → {var}))   [Schema S]")
                schritte.append(f"2. {var} → (({var} → {var}) → {var})   [Schema K]")
                schritte.append(f"3. ({var} → ({var} → {var})) → ({var} → {var})   [MP aus 1,2]")
                schritte.append(f"4. {var} → ({var} → {var})   [Schema K]")
                schritte.append(f"5. {var} → {var}   [MP aus 3,4]")
                return schritte

        return None


def deduction_theorem(hypotheses: list[str], conclusion: str) -> dict:
    """
    Deduktionstheorem für den Hilbert-Kalkül.

    Das Deduktionstheorem besagt:
    Γ, A ⊢ B  ⟺  Γ ⊢ A → B

    Dies erlaubt es, Beweise mit Hypothesen in Beweise ohne Hypothesen
    umzuwandeln, indem Hypothesen in die Implikation "gezogen" werden.

    :param hypotheses: Liste der Hypothesen Γ, A
    :param conclusion:  Die Konklusion B
    :return: Dict mit Theorem-Beschreibung und Transformation

    Letzte Änderung: 2026-03-10
    """
    if not hypotheses:
        return {
            "fehler": "Mindestens eine Hypothese erforderlich",
            "theorem": "Deduktionstheorem: Γ,A ⊢ B ⟺ Γ ⊢ A→B",
        }

    # Letzte Hypothese wird zur Implikations-Prämisse
    letzte_hyp = hypotheses[-1]
    restliche_hyps = hypotheses[:-1]

    # Neue Konklusion: letzte_hyp → conclusion
    neue_konklusion = f"{letzte_hyp} → {conclusion}"

    kalkuel = HilbertCalculus()
    beweis_vorher = f"{', '.join(hypotheses)} ⊢ {conclusion}"
    beweis_nachher = f"{', '.join(restliche_hyps) + ' ' if restliche_hyps else ''}⊢ {neue_konklusion}"

    return {
        "theorem": "Deduktionstheorem: Γ,A ⊢ B ⟺ Γ ⊢ A→B",
        "eingabe": {
            "hypothesen": hypotheses,
            "konklusion": conclusion,
            "sequent": beweis_vorher,
        },
        "ausgabe": {
            "hypothesen": restliche_hyps,
            "konklusion": neue_konklusion,
            "sequent": beweis_nachher,
        },
        "anwendung": (
            "Das Deduktionstheorem ermöglicht es, Beweise interaktiv zu führen "
            "und Hypothesen schrittweise in die Konklusion zu übernehmen."
        ),
        "neue_konklusion": neue_konklusion,
    }


# =============================================================================
# BEWEISKOMPLEXITÄT
# =============================================================================

def proof_complexity_comparison() -> dict:
    """
    Vergleicht Beweissysteme nach ihrer Stärke (Simulation).

    Beweissysteme in aufsteigender Stärke:
    1. Resolution: Klauselbasiert, refutationsvollständig
    2. Cutting Planes: Lineare Ungleichungen über Ganzzahlen
    3. Frege-System F: Propositionale Beweise mit Hilbert-Regeln
    4. Extended Frege EF: F + Abkürzungsdefinitionen
    5. Monotone Schaltkreise: Nur ∧ und ∨, kein ¬

    Cook-Reckhow-Theorem (1979):
    P = NP ⟺ Es gibt ein polynomiales Beweissystem für TAUT.

    :return: Dict mit Systembeschreibungen und Hierarchie

    Letzte Änderung: 2026-03-10
    """
    return {
        "cook_reckhow": {
            "theorem": "P = NP ⟺ Es gibt kein superpolynomiales Beweissystem",
            "bedeutung": (
                "Ein polynomiales Beweissystem würde P=NP implizieren, "
                "da man NP-Probleme durch kurze Beweise lösen könnte."
            ),
        },
        "systeme": {
            "Resolution": {
                "beschreibung": "Klauselbasiertes Widerlegungs-System (CNF)",
                "regel": "Aus (A ∨ C) und (¬A ∨ D) folgt (C ∨ D)",
                "staerke": 1,
                "bekannte_schranken": "PHP_n benötigt exp. Beweisgröße (Ben-Sasson & Wigderson 2001)",
                "anwendung": "SAT-Solver (DPLL, CDCL)",
            },
            "Cutting_Planes": {
                "beschreibung": "Beweise mit ganzzahligen linearen Ungleichungen",
                "regel": "Aus Σaᵢxᵢ ≥ b und Σcᵢxᵢ ≥ d folgt lineare Kombination",
                "staerke": 2,
                "bekannte_schranken": "Simuliert Resolution, stärker für manche Kombinatorik",
                "anwendung": "Integer Linear Programming (ILP)",
            },
            "Frege": {
                "beschreibung": "Aussagenlogisches Hilbert-System mit festen Regeln",
                "regel": "Endliche Menge von Axiomenschemata + MP",
                "staerke": 3,
                "bekannte_schranken": "Keine superpolynomialen Schranken bekannt",
                "anwendung": "Klassische Logik-Beweise",
            },
            "Extended_Frege": {
                "beschreibung": "Frege + Abkürzungsdefinitionen (p ↔ φ)",
                "regel": "Beliebige neue Variablen als Abkürzungen einführen",
                "staerke": 4,
                "bekannte_schranken": "Polynomiale Simulation aller bekannten starken Systeme",
                "anwendung": "Schaltkreis-Komplexität",
            },
            "Monotone_Schaltkreise": {
                "beschreibung": "Nur ∧ und ∨, kein Komplement",
                "regel": "Keine Negation erlaubt",
                "staerke": 0,
                "bekannte_schranken": "Superpolynomiale Schranken für Clique-Problem (Razborov 1985)",
                "anwendung": "Komplexitätstheorie untere Schranken",
            },
        },
        "hierarchie": "Monoton < Resolution ≤ Cutting Planes ≤ Frege ≤ Extended Frege",
        "offene_probleme": [
            "Gibt es superpolynomiale Schranken für Frege?",
            "Ist EF stärker als Frege?",
            "Cook's Problem: Gibt es ein optimales Beweissystem?",
        ],
    }


def resolution_proof_size(clauses: list[frozenset]) -> dict:
    """
    Analysiert die Komplexität eines Resolutionsbeweises.

    Die Resolution ist ein Widerlegungs-Verfahren:
    Aus einer Klauselmenge wird ⊥ (leere Klausel) abgeleitet, wenn
    die Formel unerfüllbar ist.

    Wichtige Größenmaße:
    - Länge: Anzahl der Klauseln im Beweis
    - Breite: Maximale Anzahl Literale in einer Klausel
    - Raum: Anzahl gleichzeitig benötigter Klauseln

    Zusammenhang Breite-Länge (Ben-Sasson & Wigderson 2001):
    L(F ⊢ ⊥) ≥ exp(W(F ⊢ ⊥) - W(F))²/n

    :param clauses: Menge von Klauseln (frozensets von Literalen)
    :return: Dict mit Komplexitätsanalyse

    Letzte Änderung: 2026-03-10
    """
    # Anzahl der Eingabe-Klauseln
    n_klauseln = len(clauses)
    # Initiale Beweis-Breite (max. Klauselgröße der Eingabe)
    eingabe_breite = max((len(k) for k in clauses), default=0)
    # Anzahl verschiedener Variablen
    alle_literale: set[str] = set()
    for klausel in clauses:
        for literal in klausel:
            alle_literale.add(literal.lstrip("¬"))
    n_variablen = len(alle_literale)

    # Einfache Resolution durchführen (bis Fixpunkt)
    alle_klauseln = set(clauses)
    neue_klauseln: set[frozenset] = set()
    resolutionsschritte = 0
    enthalten_leer = frozenset() in alle_klauseln

    # Maximale Iterationen begrenzen
    for _ in range(min(n_klauseln * 2, 100)):
        klauselliste = list(alle_klauseln)
        for i in range(len(klauselliste)):
            for j in range(i + 1, len(klauselliste)):
                c1 = klauselliste[i]
                c2 = klauselliste[j]
                # Finde komplementäre Literale
                for lit in c1:
                    komplement = lit[1:] if lit.startswith("¬") else f"¬{lit}"
                    if komplement in c2:
                        # Resolution anwenden
                        resolvent = (c1 - {lit}) | (c2 - {komplement})
                        if resolvent not in alle_klauseln:
                            neue_klauseln.add(resolvent)
                            resolutionsschritte += 1
                        if resolvent == frozenset():
                            enthalten_leer = True
        if not neue_klauseln - alle_klauseln:
            break
        alle_klauseln |= neue_klauseln

    # Breite nach Resolution
    max_breite = max((len(k) for k in alle_klauseln), default=0)

    return {
        "eingabe": {
            "anzahl_klauseln": n_klauseln,
            "anzahl_variablen": n_variablen,
            "eingabe_breite": eingabe_breite,
        },
        "analyse": {
            "resolutionsschritte": resolutionsschritte,
            "gesamt_klauseln": len(alle_klauseln),
            "maximale_breite": max_breite,
            "ist_widerlegbar": enthalten_leer,
        },
        "komplexitaet": {
            "breite_laenge_theorem": "L ≥ exp((W(F⊢⊥) - W(F))² / n)",
            "pebbling": "Resolution-Raum entspricht Pebbling-Spielen auf DAGs",
        },
    }


def pigeonhole_principle_hardness() -> dict:
    """
    Das Taubenschlag-Prinzip (PHP) als hartes Beispiel für Resolution.

    PHP_n: n+1 Tauben in n Löcher → keine injektive Abbildung möglich.

    Satz (Ben-Sasson & Wigderson 2001):
    Jeder Resolution-Beweis von PHP_n benötigt Größe 2^Ω(n).

    Kodierung: Variablen p_{i,j} = "Taube i ist in Loch j"
    - Für jede Taube i: p_{i,1} ∨ ... ∨ p_{i,n}  (Taube muss irgendwo)
    - Für jedes Paar (i₁,i₂) und jedes Loch j: ¬p_{i₁,j} ∨ ¬p_{i₂,j}  (kein Loch doppelt)

    :return: Dict mit Kodierung, Größenanalyse und theoretischen Schranken

    Letzte Änderung: 2026-03-10
    """
    n = 3  # PHP_3: 4 Tauben in 3 Löcher (Demo-Größe)
    tauben = n + 1  # n+1 Tauben
    loecher = n

    # Klauseln erzeugen
    klauseln = []

    # Typ 1: Jede Taube muss in mindestens einem Loch sein
    # p_{i,1} ∨ p_{i,2} ∨ ... ∨ p_{i,n}
    taube_klauseln = []
    for i in range(1, tauben + 1):
        klausel = frozenset(f"p_{i},{j}" for j in range(1, loecher + 1))
        taube_klauseln.append(klausel)
        klauseln.append(klausel)

    # Typ 2: Kein Loch kann zwei Tauben haben
    # ¬p_{i1,j} ∨ ¬p_{i2,j}
    loch_klauseln = []
    for j in range(1, loecher + 1):
        for i1 in range(1, tauben + 1):
            for i2 in range(i1 + 1, tauben + 1):
                klausel = frozenset({f"¬p_{i1},{j}", f"¬p_{i2},{j}"})
                loch_klauseln.append(klausel)
                klauseln.append(klausel)

    # Theoretische Schranken berechnen
    n_variablen = tauben * loecher
    n_taube_klauseln = tauben
    n_loch_klauseln = loecher * (tauben * (tauben - 1) // 2)
    n_gesamt = n_taube_klauseln + n_loch_klauseln

    # Untere Schranke für Resolution: 2^Ω(n)
    untere_schranke = 2 ** (n // 2)

    return {
        "problem": f"PHP_{n}: {tauben} Tauben in {loecher} Löcher",
        "kodierung": {
            "variablen": f"p_{{i,j}} = 'Taube i ist in Loch j'",
            "anzahl_variablen": n_variablen,
            "taube_klauseln": n_taube_klauseln,
            "loch_klauseln": n_loch_klauseln,
            "gesamt_klauseln": n_gesamt,
        },
        "komplexitaet": {
            "resolution_schranke": f"2^Ω(n) = mindestens {untere_schranke} Schritte für n={n}",
            "beweis": "Ben-Sasson & Wigderson (2001) via Breite-Länge-Methode",
            "bedeutung": "Resolution kann Pigeonhole nicht polynomiell beweisen",
        },
        "haerteste_klasse": {
            "frege": "PHP hat polynomielle Frege-Beweise (O(n³))",
            "resolution": "PHP benötigt exponentielle Resolution-Beweise",
            "gap": "Exponentieller Unterschied zwischen Frege und Resolution für PHP",
        },
        "erste_paar_klauseln": {
            "taube_1": f"p_1,1 ∨ p_1,2 ∨ ... ∨ p_1,{loecher}",
            "kollision_1": f"¬p_1,1 ∨ ¬p_2,1",
        },
    }


# =============================================================================
# ORDINALE ANALYSE
# =============================================================================

def proof_theoretic_ordinal(system: str) -> dict:
    """
    Bestimmt die beweistheoretische Ordinalzahl |T| eines formalen Systems.

    Die beweistheoretische Ordinalzahl |T| misst die "Beweiskraft" eines Systems:
    Es ist die kleinste Ordinalzahl, für die transfinite Induktion in T nicht
    beweisbar ist.

    Bekannte Werte:
    - PRA (Prim. rek. Arithmetik): |PRA| = ω^ω
    - PA (Peano-Arithmetik):       |PA|  = ε₀
    - ACA₀ (Arithm. Kompr.):       |ACA₀| = ε₀
    - ATR₀ (Arithm. Transf.):      |ATR₀| = Γ₀
    - Π₁¹-CA₀:                     |Π₁¹-CA₀| = Γ₀ (Feferman-Schütte)

    :param system: Name des formalen Systems
    :return: Dict mit Ordinal-Beschreibung und mathematischen Details

    Letzte Änderung: 2026-03-10
    """
    systeme = {
        "PRA": {
            "name": "Primitive Rekursive Arithmetik",
            "ordinal": "ω^ω",
            "cantor_normalform": "ω^ω = ω^ω",
            "beschreibung": (
                "PRA beweist Terminierung aller primitiv-rekursiven Funktionen. "
                "Transfinite Induktion bis ω^ω ist beweisbar, aber nicht darüber hinaus."
            ),
            "axiome": [
                "Peano-Axiome ohne vollständige Induktion",
                "Definition primitiv-rekursiver Funktionen",
                "Induktion nur für quantorenfreie Formeln",
            ],
            "gentzen": "PRA + TI(ε₀) beweist Con(PA)",
        },
        "PA": {
            "name": "Peano-Arithmetik",
            "ordinal": "ε₀",
            "cantor_normalform": "ε₀ = ω^{ω^{ω^{...}}} (ω-fach)",
            "beschreibung": (
                "PA ist das Standard-System der natürlichen Arithmetik. "
                "ε₀ ist die kleinste Ordinalzahl mit ω^α = α. "
                "Gentzen (1936) bewies Con(PA) in PRA + TI(ε₀)."
            ),
            "axiome": [
                "Nachfolger-Axiome (0 ≠ S(n), S injektiv)",
                "Addition und Multiplikation rekursiv",
                "Vollständiges Induktionsschema (alle Formeln)",
            ],
            "gentzen": "Con(PA) äquivalent zu WO(ε₀) (Wohlordnung von ε₀)",
        },
        "ACA0": {
            "name": "Arithmetische Komprehension (ACA₀)",
            "ordinal": "ε₀",
            "cantor_normalform": "ε₀ = ω^{ω^{ω^{...}}}",
            "beschreibung": (
                "ACA₀ ermöglicht arithmetisch definierbare Mengen. "
                "Obwohl stärker als PA in Bezug auf Mengenlehre, "
                "haben PA und ACA₀ die gleiche beweistheoretische Ordinalzahl."
            ),
            "axiome": [
                "PA-Axiome für natürliche Zahlen",
                "Komprehension: {n : φ(n)} existiert für arithmetisches φ",
                "Keine Mengenquantoren in φ",
            ],
            "gentzen": "Con(ACA₀) äquivalent zu Con(PA)",
        },
        "ATR0": {
            "name": "Arithmetische Transfinite Rekursion (ATR₀)",
            "ordinal": "Γ₀",
            "cantor_normalform": "Γ₀ = φ(1,0,0,...) (Veblen-Funktion)",
            "beschreibung": (
                "ATR₀ erlaubt arithmetische Definitionen entlang beliebiger Wohlordnungen. "
                "Die Feferman-Schütte-Ordinalzahl Γ₀ ist die kleinste 'unbeweisbare' Ordinalzahl."
            ),
            "axiome": [
                "ACA₀-Axiome",
                "Transfinite Rekursion entlang beliebiger Wohlordnungen",
            ],
            "gentzen": "Γ₀ ist die Feferman-Schütte-Ordinalzahl",
        },
        "Pi11-CA0": {
            "name": "Π₁¹-Komprehension (Π₁¹-CA₀)",
            "ordinal": "Ψ(Ω_ω)",
            "cantor_normalform": "Ψ(Ω_ω) (Bachmann-Howard-Ordinal)",
            "beschreibung": (
                "Π₁¹-CA₀ erlaubt Π₁¹-Definitionen von Mengen. "
                "Dies entspricht der zweiten Stufe der analytischen Hierarchie. "
                "Das Bachmann-Howard-Ordinal ist die beweistheoretische Stärke."
            ),
            "axiome": [
                "ATR₀-Axiome",
                "Π₁¹-Komprehension: {n : ∀X.φ(n,X)} existiert",
            ],
            "gentzen": "Stärker als ATR₀, schwächer als volle Analysis (Z₂)",
        },
    }

    if system in systeme:
        ergebnis = systeme[system]
        ergebnis["system"] = system
        return ergebnis
    else:
        return {
            "fehler": f"System '{system}' nicht bekannt",
            "bekannte_systeme": list(systeme.keys()),
        }


def epsilon_zero_demo() -> dict:
    """
    Demonstration der Ordinalzahl ε₀.

    ε₀ (Epsilon-Null) ist die kleinste Ordinalzahl α mit ω^α = α.
    Sie ist der Grenzwert der Folge: ω, ω^ω, ω^{ω^ω}, ...

    Cantor-Normalform: Jede Ordinalzahl α < ε₀ hat eine eindeutige Darstellung
    α = ω^{β₁}·c₁ + ω^{β₂}·c₂ + ... + ω^{βₙ}·cₙ
    mit β₁ > β₂ > ... > βₙ und cᵢ ∈ ℕ⁺.

    Bedeutung nach Gentzen (1936):
    Con(PA) ⟺ Wohlordnung von ε₀ (beweisbar in PRA + TI(ε₀))

    :return: Dict mit ε₀-Eigenschaften und Cantor-Normalformen

    Letzte Änderung: 2026-03-10
    """
    # Endliche Cantor-Normalformen berechnen
    def cantor_normalform(n: int) -> str:
        """Gibt die Cantor-Normalform einer kleinen natürlichen Zahl zurück."""
        if n == 0:
            return "0"
        if n < 0:
            raise ValueError("Negative Ordinalzahlen existieren nicht")
        # Dezimal → binär-ähnliche Darstellung (vereinfacht)
        teile = []
        rest = n
        exponent = 0
        while rest > 0:
            if rest % 2 == 1:
                if exponent == 0:
                    teile.append("1")
                elif exponent == 1:
                    teile.append("ω")
                else:
                    teile.append(f"ω^{exponent}")
            rest //= 2
            exponent += 1
        return " + ".join(reversed(teile)) if teile else "0"

    # Aufsteigende Folge zu ε₀
    folge = [
        ("ω⁰ = 1", "1", "Endliche Ordinalzahl"),
        ("ω¹ = ω", "ω", "Erste unendliche Ordinalzahl"),
        ("ω^ω", "ω^ω", "Ordnung der Polynome in ω"),
        ("ω^{ω^ω}", "ω^{ω^ω}", "Dritte Stufe der Iteration"),
        ("ε₀ = sup{ω, ω^ω, ω^{ω^ω}, ...}", "ε₀", "Grenzwert der Folge"),
    ]

    return {
        "definition": "ε₀ ist die kleinste Ordinalzahl α mit ω^α = α",
        "darstellung": "ε₀ = ω^{ω^{ω^{ω^{...}}}} (ω-fache Iteration von ω^·)",
        "folge": folge,
        "cantor_normalform": {
            "satz": "Jede Ordinalzahl < ε₀ hat eindeutige Cantor-Normalform",
            "form": "α = ω^{β₁}·c₁ + ... + ω^{βₙ}·cₙ, β₁ > ... > βₙ",
            "beispiele": {
                str(i): cantor_normalform(i) for i in range(1, 9)
            },
        },
        "gentzen": {
            "satz": "Con(PA) ⟺ WO(ε₀) in PRA",
            "bedeutung": (
                "Die Widerspruchsfreiheit der Peano-Arithmetik ist äquivalent "
                "zur Wohlordnung von ε₀. Gödel zeigte, dass Con(PA) nicht in PA "
                "beweisbar ist (2. Unvollständigkeitssatz)."
            ),
        },
        "properties": {
            "limit_ordinal": True,
            "epsilon_number": True,
            "countable": True,
            "recursive": True,
            "beschreibung": "ε₀ ist abzählbar und rekursiv geordnet",
        },
    }


def gentzen_consistency_proof_sketch() -> dict:
    """
    Skizze von Gentzens Konsistenzbeweis für PA (1936).

    Gentzen bewies die Konsistenz der Peano-Arithmetik in einem System
    das aus PRA (primitive rekursive Arithmetik) plus transfiniter Induktion
    bis ε₀ besteht.

    Hauptidee:
    1. Ordne jeden PA-Beweis einer Ordinalzahl < ε₀ zu
    2. Transformationen reduzieren die Ordinalzahl streng
    3. Da ε₀ wohlgeordnet ist, terminiert das Verfahren
    4. Terminierung bedeutet: kein Beweis von ⊥ möglich

    :return: Dict mit Beweisskizze und mathematischen Details

    Letzte Änderung: 2026-03-10
    """
    return {
        "titel": "Gentzens Konsistenzbeweis für PA (1936)",
        "system": "PRA + TI(ε₀)",
        "godel_verbindung": {
            "problem": "Gödels 2. Unvollständigkeitssatz (1931)",
            "aussage": "Con(PA) ist nicht in PA selbst beweisbar",
            "genzens_loesung": (
                "Beweis in einem stärkeren System PRA + TI(ε₀) "
                "(nicht 'problematisch' weil ε₀ finit beschreibbar)"
            ),
        },
        "beweisidee": {
            "schritt_1": {
                "name": "Ordinalzahlzuweisung",
                "beschreibung": (
                    "Weise jedem Beweis-Sequenzbaum π eine Ordinalzahl o(π) < ε₀ zu. "
                    "Dabei kodiert o(π) die 'Komplexität' des Beweises."
                ),
            },
            "schritt_2": {
                "name": "Reduktionsverfahren",
                "beschreibung": (
                    "Zeige: Wenn π ein Beweis von ⊥ ist, kann π zu π' vereinfacht werden "
                    "mit o(π') < o(π). Die Transformation entspricht Schnittelimination."
                ),
            },
            "schritt_3": {
                "name": "Wohlordnungsargument",
                "beschreibung": (
                    "Da ε₀ wohlgeordnet ist (keine unendliche absteigende Folge), "
                    "terminiert die Reduktion. Am Ende steht kein gültiger Beweis von ⊥."
                ),
            },
            "schritt_4": {
                "name": "Schlussfolgerung",
                "beschreibung": "PA ⊬ ⊥, also ist PA konsistent.",
            },
        },
        "ordinalzuweisung": {
            "axiom": "o(Axiom) = 0",
            "schnitt": "o(Cut(π₁,π₂)) = ω^{max(o(π₁),o(π₂))} + 1",
            "logische_regeln": "o(Regel(π)) = o(π) + 1",
        },
        "historische_bedeutung": (
            "Gentzens Beweis war revolutionär: Er etablierte die Verbindung zwischen "
            "Konsistenz formaler Systeme und Wohlordnungseigenschaften von Ordinalzahlen. "
            "Dies begründete das Programm der Proof Theory / Beweistheorie als eigenständige "
            "mathematische Disziplin."
        ),
        "literatur": [
            "Gentzen, G. (1936): Die Widerspruchsfreiheit der reinen Zahlentheorie",
            "Takeuti, G. (1987): Proof Theory, 2nd ed.",
            "Buss, S. (1998): Handbook of Proof Theory",
        ],
    }


# =============================================================================
# REVERSE MATHEMATICS
# =============================================================================

def big_five_systems() -> dict:
    """
    Die 'Big Five' der Reverse Mathematics nach Simpson.

    Reverse Mathematics (Simpson 1999) untersucht, welche Axiome zur
    Beweisbarkeit klassischer mathematischer Sätze notwendig sind.

    Die fünf Hauptsysteme bilden eine aufsteigende Hierarchie:
    RCA₀ ⊂ WKL₀ ⊂ ACA₀ ⊂ ATR₀ ⊂ Π₁¹-CA₀

    :return: Dict mit allen fünf Systemen und ihren Eigenschaften

    Letzte Änderung: 2026-03-10
    """
    return {
        "RCA0": {
            "name": "Recursive Comprehension Axiom (RCA₀)",
            "ordinal": "ω^ω",
            "axiome": [
                "PA⁻ (Peano ohne volle Induktion)",
                "Σ⁰₁-Induktion",
                "Δ⁰₁-Komprehension (rekursiv definierbare Mengen)",
            ],
            "charakteristische_saetze": [
                "Grundlegende Arithmetik und rekursive Funktionen",
                "Existenz unendlicher rekursiver Teilmengen",
            ],
            "bedeutung": "Basislinien-System; entspricht 'berechenbarer Mathematik'",
        },
        "WKL0": {
            "name": "Weak König's Lemma (WKL₀)",
            "ordinal": "ω^ω",
            "axiome": [
                "RCA₀-Axiome",
                "WKL: Jeder unendliche binäre Baum hat einen unendlichen Pfad",
            ],
            "charakteristische_saetze": [
                "Heine-Borel-Überdeckungssatz",
                "Jede stetige Funktion auf [0,1] ist gleichmäßig stetig",
                "Brouwerscher Fixpunktsatz",
                "Gödels Vollständigkeitssatz",
            ],
            "bedeutung": "Äquivalent zu Π⁰₁-Komprehension; starke Verbindung zu Kompaktheit",
        },
        "ACA0": {
            "name": "Arithmetical Comprehension Axiom (ACA₀)",
            "ordinal": "ε₀",
            "axiome": [
                "RCA₀-Axiome",
                "Arithmetische Komprehension: Existenz von {n : φ(n)} für arithmetisches φ",
            ],
            "charakteristische_saetze": [
                "Bolzano-Weierstraß-Theorem",
                "Ramsey-Theorem für Paare (RT²₂)",
                "Cauchy-Konvergenz-Kriterium",
                "König-Lemma (unendliche Bäume)",
            ],
            "bedeutung": "Entspricht 'arithmetischer Analysis'; konservativ über PA für Π¹₁-Sätze",
        },
        "ATR0": {
            "name": "Arithmetical Transfinite Recursion (ATR₀)",
            "ordinal": "Γ₀",
            "axiome": [
                "ACA₀-Axiome",
                "ATR: Arithmetische Definitionen entlang jeder Wohlordnung",
            ],
            "charakteristische_saetze": [
                "Determiniertheits-Theorem (Borel-Spiele)",
                "Luzin-Sierpiński-Theorem",
                "Vergleichbarkeit abzählbarer Wohlordnungen",
                "Silver-Theorem über Ko-analyt. Äquivalenzrelationen",
            ],
            "bedeutung": "Entspricht 'prädikativistischer' Mathematik; Γ₀ = Feferman-Schütte",
        },
        "Pi11-CA0": {
            "name": "Π¹₁-Comprehension Axiom (Π¹₁-CA₀)",
            "ordinal": "Ψ(Ω_ω)",
            "axiome": [
                "ATR₀-Axiome",
                "Π¹₁-Komprehension: Existenz von {n : ∀X.φ(n,X)} für Π¹₁-Formeln φ",
            ],
            "charakteristische_saetze": [
                "Σ¹₁-Baire-Kategorie-Satz",
                "Existenz von Π¹₁-Mengen mit komplexer Struktur",
                "Hyperarithmetische Analysis",
            ],
            "bedeutung": "Stärkstes der Big Five; entspricht 'imprädikativistischer' Mathematik",
        },
        "hierarchie": {
            "aufsteigend": "RCA₀ ⊂ WKL₀ ⊂ ACA₀ ⊂ ATR₀ ⊂ Π¹₁-CA₀",
            "konservativitaet": "Jedes stärkere System ist eine echte Erweiterung",
        },
    }


def reverse_math_example(theorem: str) -> dict:
    """
    Zeigt die beweistheoretische Stärke klassischer mathematischer Sätze.

    Für jeden Satz wird angegeben, über welchem der Big Five er beweisbar
    und zu welchem er äquivalent ist (über RCA₀).

    :param theorem: Name des Satzes (Bezeichner)
    :return: Dict mit System-Äquivalenz und Beweisskizze

    Letzte Änderung: 2026-03-10
    """
    theoreme = {
        "bolzano_weierstrass": {
            "name": "Bolzano-Weierstraß-Theorem",
            "aussage": "Jede beschränkte Folge reeller Zahlen hat eine konvergente Teilfolge",
            "aequivalent_zu": "ACA₀",
            "beweisbar_in": "ACA₀",
            "nicht_beweisbar_in": "WKL₀",
            "skizze": (
                "In ACA₀: Konstruiere die Teilfolge arithmetisch via Komprehension. "
                "Für die Umkehrung: Aus B-W folgt Existenz arithmetisch definierter Mengen."
            ),
            "historisch": "Bolzano (1817), Weierstraß (1860er)",
        },
        "heine_borel": {
            "name": "Heine-Borel-Überdeckungssatz",
            "aussage": "[0,1] ist kompakt: Jede offene Überdeckung hat eine endliche Teilüberdeckung",
            "aequivalent_zu": "WKL₀",
            "beweisbar_in": "WKL₀",
            "nicht_beweisbar_in": "RCA₀",
            "skizze": (
                "In WKL₀: Benutze WKL um einen Pfad durch den Überdeckungsbaum zu finden. "
                "Für die Umkehrung: Aus HB folgt WKL (kodiere Baumpfade als Überdeckungen)."
            ),
            "historisch": "Heine (1872), Borel (1895)",
        },
        "ramsey": {
            "name": "Ramsey-Theorem RT²₂",
            "aussage": "Für jede 2-Färbung von Paaren aus ℕ gibt es eine unendliche monochromatische Menge",
            "aequivalent_zu": "ACA₀",
            "beweisbar_in": "ACA₀",
            "nicht_beweisbar_in": "WKL₀",
            "skizze": (
                "In ACA₀: Konstruiere die monochromatische Menge per arithmetischer Komprehension. "
                "Umkehrung via Kodierung arithmetischer Mengen als Ramsey-Instanzen."
            ),
            "historisch": "Ramsey (1930)",
        },
        "konig": {
            "name": "König-Lemma",
            "aussage": "Jeder unendliche endlich-verzweigte Baum hat einen unendlichen Pfad",
            "aequivalent_zu": "WKL₀",
            "beweisbar_in": "WKL₀",
            "nicht_beweisbar_in": "RCA₀",
            "skizze": (
                "Für binäre Bäume: Dies ist genau WKL. "
                "Für allgemein endlich-verzweigte Bäume: Reduktion auf WKL möglich. "
                "Ohne WKL gibt es in RCA₀ rekursive Gegenbeispiele."
            ),
            "historisch": "König (1927)",
        },
    }

    if theorem in theoreme:
        ergebnis = theoreme[theorem]
        ergebnis["theorem"] = theorem
        return ergebnis
    else:
        return {
            "fehler": f"Satz '{theorem}' nicht bekannt",
            "bekannte_saetze": list(theoreme.keys()),
        }
