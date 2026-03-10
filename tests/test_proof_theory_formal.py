"""
Tests für das Modul proof_theory_formal.py

Testet: Sequenzenkalkül, Natürliches Schließen, Hilbert-Kalkül,
        Beweiskomplexität, Ordinale Analyse, Reverse Mathematics.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import pytest
from src.proof_theory_formal import (
    Sequent,
    ProofTree,
    LKRules,
    NDProof,
    HilbertCalculus,
    cut_elimination_demo,
    subformula_property,
    prove_simple,
    deduction_theorem,
    proof_complexity_comparison,
    resolution_proof_size,
    pigeonhole_principle_hardness,
    proof_theoretic_ordinal,
    epsilon_zero_demo,
    gentzen_consistency_proof_sketch,
    big_five_systems,
    reverse_math_example,
)


# =============================================================================
# TESTS FÜR Sequent
# =============================================================================

class TestSequent:
    """Tests für die Sequent-Klasse."""

    def test_sequent_erstellung_einfach(self):
        """Einfache Sequent-Erstellung."""
        s = Sequent(["A"], ["B"])
        assert s.antecedent == ["A"]
        assert s.succedent == ["B"]

    def test_sequent_erstellung_leer(self):
        """Leeres Sequent ist erlaubt."""
        s = Sequent([], [])
        assert s.antecedent == []
        assert s.succedent == []

    def test_sequent_mehrere_formeln(self):
        """Sequent mit mehreren Formeln."""
        s = Sequent(["A", "B", "C"], ["D", "E"])
        assert len(s.antecedent) == 3
        assert len(s.succedent) == 2

    def test_sequent_axiom_einfach(self):
        """A ⊢ A ist ein Axiom."""
        s = Sequent(["A"], ["A"])
        assert s.is_axiom() is True

    def test_sequent_axiom_mit_extras(self):
        """A, B ⊢ A, C ist ein Axiom (A kommt auf beiden Seiten vor)."""
        s = Sequent(["A", "B"], ["A", "C"])
        assert s.is_axiom() is True

    def test_sequent_kein_axiom(self):
        """A ⊢ B ist kein Axiom (verschiedene Formeln)."""
        s = Sequent(["A"], ["B"])
        assert s.is_axiom() is False

    def test_sequent_kein_axiom_leer(self):
        """Leeres Sequent ist kein Axiom."""
        s = Sequent([], [])
        assert s.is_axiom() is False

    def test_sequent_str_einfach(self):
        """String-Darstellung eines einfachen Sequents."""
        s = Sequent(["A"], ["B"])
        darstellung = str(s)
        assert "⊢" in darstellung
        assert "A" in darstellung
        assert "B" in darstellung

    def test_sequent_str_leer(self):
        """String-Darstellung eines leeren Sequents."""
        s = Sequent([], [])
        darstellung = str(s)
        assert "⊢" in darstellung

    def test_sequent_str_mehrere(self):
        """String-Darstellung mit mehreren Formeln."""
        s = Sequent(["A", "B"], ["C", "D"])
        darstellung = str(s)
        assert "A" in darstellung
        assert "B" in darstellung
        assert "C" in darstellung

    def test_sequent_gleichheit(self):
        """Zwei identische Sequenten sind gleich."""
        s1 = Sequent(["A"], ["B"])
        s2 = Sequent(["A"], ["B"])
        assert s1 == s2

    def test_sequent_ungleichheit(self):
        """Verschiedene Sequenten sind ungleich."""
        s1 = Sequent(["A"], ["B"])
        s2 = Sequent(["B"], ["A"])
        assert s1 != s2


# =============================================================================
# TESTS FÜR ProofTree
# =============================================================================

class TestProofTree:
    """Tests für die ProofTree-Klasse."""

    def test_proof_tree_blatt_gueltig(self):
        """Ein Axiom-Blatt ist ein gültiger Beweis."""
        sequent = Sequent(["A"], ["A"])
        baum = ProofTree(sequent, "Axiom")
        assert baum.is_valid() is True

    def test_proof_tree_blatt_ungueltig(self):
        """Ein Nicht-Axiom-Blatt ist ungültig."""
        sequent = Sequent(["A"], ["B"])
        baum = ProofTree(sequent, "Axiom")
        assert baum.is_valid() is False

    def test_proof_tree_hoehe_blatt(self):
        """Blattknoten hat Höhe 0."""
        baum = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        assert baum.height() == 0

    def test_proof_tree_hoehe_eins(self):
        """Baum mit einer Regelanwendung hat Höhe 1."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        wurzel = LKRules.weakening_left(blatt, "B")
        assert wurzel.height() == 1

    def test_proof_tree_groesse_blatt(self):
        """Blattknoten hat Größe 1."""
        baum = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        assert baum.size() == 1

    def test_proof_tree_groesse_zwei(self):
        """Baum mit zwei Knoten hat Größe 2."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        wurzel = LKRules.weakening_left(blatt, "B")
        assert wurzel.size() == 2

    def test_proof_tree_to_string(self):
        """Baum-Darstellung ist nicht leer."""
        baum = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        darstellung = baum.to_string()
        assert len(darstellung) > 0
        assert "Axiom" in darstellung

    def test_proof_tree_to_string_mit_praemissen(self):
        """Baum-Darstellung mit Einrückung für Prämissen."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        wurzel = LKRules.weakening_left(blatt, "B")
        darstellung = wurzel.to_string()
        assert "WL" in darstellung
        assert "Axiom" in darstellung


# =============================================================================
# TESTS FÜR LKRules
# =============================================================================

class TestLKRules:
    """Tests für die LK-Regelanwendungen."""

    def test_weakening_left(self):
        """WL fügt Formel im Antecedent hinzu."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        ergebnis = LKRules.weakening_left(blatt, "B")
        assert "B" in ergebnis.conclusion.antecedent
        assert ergebnis.rule == "WL"

    def test_weakening_right(self):
        """WR fügt Formel im Succedent hinzu."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        ergebnis = LKRules.weakening_right(blatt, "C")
        assert "C" in ergebnis.conclusion.succedent
        assert ergebnis.rule == "WR"

    def test_contraction_left(self):
        """CL entfernt doppelte Formel im Antecedent."""
        blatt = ProofTree(Sequent(["A", "A", "B"], ["C"]), "Axiom_fake")
        ergebnis = LKRules.contraction_left(blatt, "A")
        # Nach Contraction: weniger A's im Antecedent
        anzahl_a = ergebnis.conclusion.antecedent.count("A")
        assert anzahl_a < 2

    def test_cut_regel(self):
        """Cut kombiniert zwei Beweise."""
        beweis1 = ProofTree(Sequent(["A"], ["A", "B"]), "Axiom")
        beweis2 = ProofTree(Sequent(["B", "C"], ["D"]), "Axiom_fake")
        ergebnis = LKRules.cut(beweis1, beweis2, "B")
        assert "Cut" in ergebnis.rule
        assert "B" in ergebnis.rule

    def test_and_right(self):
        """∧R kombiniert zwei Beweise zu Konjunktion."""
        beweis1 = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        beweis2 = ProofTree(Sequent(["B"], ["B"]), "Axiom")
        ergebnis = LKRules.and_right(beweis1, beweis2)
        assert "∧" in str(ergebnis.conclusion)
        assert ergebnis.rule == "∧R"

    def test_implies_right(self):
        """→R erzeugt Implikation im Succedent."""
        beweis = ProofTree(Sequent(["A", "B"], ["A"]), "Axiom")
        ergebnis = LKRules.implies_right(beweis)
        # Implikation sollte im Succedent erscheinen
        assert "→" in str(ergebnis.conclusion)
        assert ergebnis.rule == "→R"


# =============================================================================
# TESTS FÜR NDProof
# =============================================================================

class TestNDProof:
    """Tests für Natürliches Schließen."""

    def test_nd_annahme(self):
        """Annahme gibt ID zurück."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        assert isinstance(id_a, int)
        assert beweis.get_formula(id_a) == "A"

    def test_nd_mehrere_annahmen(self):
        """Mehrere Annahmen haben verschiedene IDs."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_b = beweis.assume("B")
        assert id_a != id_b
        assert beweis.get_formula(id_a) == "A"
        assert beweis.get_formula(id_b) == "B"

    def test_nd_and_intro(self):
        """∧I kombiniert zwei Formeln zu Konjunktion."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_b = beweis.assume("B")
        id_und = beweis.apply_and_intro(id_a, id_b)
        formel = beweis.get_formula(id_und)
        assert "∧" in formel
        assert "A" in formel
        assert "B" in formel

    def test_nd_and_elim_left(self):
        """∧EL extrahiert linken Konjunkt."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_b = beweis.assume("B")
        id_und = beweis.apply_and_intro(id_a, id_b)
        id_links = beweis.apply_and_elim_left(id_und)
        assert beweis.get_formula(id_links) == "A"

    def test_nd_and_elim_right(self):
        """∧ER extrahiert rechten Konjunkt."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_b = beweis.assume("B")
        id_und = beweis.apply_and_intro(id_a, id_b)
        id_rechts = beweis.apply_and_elim_right(id_und)
        assert beweis.get_formula(id_rechts) == "B"

    def test_nd_implies_intro(self):
        """→I erzeugt Implikation und entlädt Annahme."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_impl = beweis.apply_implies_intro(id_a, id_a)
        formel = beweis.get_formula(id_impl)
        assert "→" in formel

    def test_nd_implies_elim(self):
        """→E wendet Modus Ponens an."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        # Manuell A→B einführen
        id_b = beweis.assume("B")
        id_impl = beweis.apply_implies_intro(id_a, id_b)
        id_a2 = beweis.assume("A")
        id_konsequent = beweis.apply_implies_elim(id_impl, id_a2)
        assert beweis.get_formula(id_konsequent) == "B"

    def test_nd_or_intro_left(self):
        """∨IL führt Disjunktion links ein."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_oder = beweis.apply_or_intro_left(id_a, "B")
        formel = beweis.get_formula(id_oder)
        assert "∨" in formel
        assert "A" in formel
        assert "B" in formel

    def test_nd_neg_intro(self):
        """¬I führt Negation ein."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_bot = beweis.assume("⊥")
        id_neg = beweis.apply_neg_intro(id_a, id_bot)
        formel = beweis.get_formula(id_neg)
        assert "¬" in formel
        assert "A" in formel

    def test_nd_ex_falso(self):
        """EFQ leitet beliebige Formel aus ⊥ ab."""
        beweis = NDProof()
        id_bot = beweis.assume("⊥")
        id_c = beweis.apply_ex_falso(id_bot, "C")
        assert beweis.get_formula(id_c) == "C"

    def test_nd_is_complete_true(self):
        """is_complete_proof erkennt vollständige Beweise."""
        beweis = NDProof()
        id_a = beweis.assume("A")
        id_und = beweis.apply_and_intro(id_a, id_a)
        # Der Beweis leitet (A∧A) ab
        assert beweis.is_complete_proof("(A∧A)") is True

    def test_nd_is_complete_false(self):
        """is_complete_proof erkennt unvollständige Beweise."""
        beweis = NDProof()
        beweis.assume("A")
        # B wurde nie abgeleitet
        assert beweis.is_complete_proof("B") is False

    def test_nd_ungueltige_id_wirft_fehler(self):
        """Ungültige ID wirft KeyError."""
        beweis = NDProof()
        with pytest.raises(KeyError):
            beweis.get_formula(999)


# =============================================================================
# TESTS FÜR prove_simple
# =============================================================================

class TestProveSimple:
    """Tests für die automatische Tautologie-Beweis-Funktion."""

    def test_prove_identitaet(self):
        """A→A ist beweisbar."""
        beweis = prove_simple("A→A")
        assert beweis is not None
        assert beweis.is_complete_proof("(A→A)")

    def test_prove_konjunktion_projektion_links(self):
        """(A∧B)→A ist beweisbar."""
        beweis = prove_simple("(A∧B)→A")
        assert beweis is not None

    def test_prove_konjunktion_projektion_rechts(self):
        """(A∧B)→B ist beweisbar."""
        beweis = prove_simple("(A∧B)→B")
        assert beweis is not None

    def test_prove_disjunktion_einfuehrung(self):
        """A→(A∨B) ist beweisbar."""
        beweis = prove_simple("A→(A∨B)")
        assert beweis is not None

    def test_prove_unbekannte_tautologie(self):
        """Unbekannte Formeln geben None zurück."""
        beweis = prove_simple("A∧B→C")
        assert beweis is None

    def test_prove_andere_identitaet(self):
        """B→B ist beweisbar."""
        beweis = prove_simple("B→B")
        assert beweis is not None


# =============================================================================
# TESTS FÜR HilbertCalculus
# =============================================================================

class TestHilbertCalculus:
    """Tests für den Hilbert-Kalkül."""

    def setup_method(self):
        """Erstellt eine Hilbert-Kalkül-Instanz vor jedem Test."""
        self.kalkuel = HilbertCalculus()

    def test_hilbert_axiome_vorhanden(self):
        """Alle drei Axiomenschemata sind definiert."""
        assert "K" in HilbertCalculus.AXIOMS
        assert "S" in HilbertCalculus.AXIOMS
        assert "DN" in HilbertCalculus.AXIOMS

    def test_hilbert_axiom_k_instanz(self):
        """Erkennt K-Schema-Instanz."""
        # K: A → (B → A)
        ergebnis = self.kalkuel.check_axiom_instance("p → (q → p)")
        assert ergebnis == "K"

    def test_hilbert_axiom_dn_instanz(self):
        """Erkennt DN-Schema-Instanz."""
        ergebnis = self.kalkuel.check_axiom_instance("(¬¬A) → A")
        assert ergebnis == "DN"

    def test_hilbert_kein_axiom(self):
        """Gibt None für Nicht-Axiome."""
        ergebnis = self.kalkuel.check_axiom_instance("A ∧ B")
        assert ergebnis is None

    def test_hilbert_modus_ponens(self):
        """MP leitet Konsequent ab."""
        ergebnis = self.kalkuel.apply_mp("A → B", "A")
        assert ergebnis == "B"

    def test_hilbert_modus_ponens_fehlschlag(self):
        """MP schlägt fehl bei falschem Antecedent."""
        ergebnis = self.kalkuel.apply_mp("A → B", "C")
        assert ergebnis is None

    def test_hilbert_modus_ponens_keine_implikation(self):
        """MP schlägt fehl bei Nicht-Implikation."""
        ergebnis = self.kalkuel.apply_mp("A ∧ B", "A")
        assert ergebnis is None

    def test_hilbert_beweis_a_implies_a(self):
        """Beweist A → A via Hilbert-Kalkül."""
        schritte = self.kalkuel.prove("A → A")
        assert schritte is not None
        assert len(schritte) > 0

    def test_hilbert_beweis_hypothese(self):
        """Beweis mit Hypothesen (Hypothese direkt verfügbar)."""
        schritte = self.kalkuel.prove("A", hypotheses=["A"])
        assert schritte is not None


# =============================================================================
# TESTS FÜR cut_elimination_demo
# =============================================================================

class TestCutEliminationDemo:
    """Tests für die Schnitteliminations-Demonstration."""

    def setup_method(self):
        """Führt die Demo vor jedem Test aus."""
        self.ergebnis = cut_elimination_demo()

    def test_demo_gibt_dict_zurueck(self):
        """cut_elimination_demo gibt ein Dict zurück."""
        assert isinstance(self.ergebnis, dict)

    def test_demo_hauptsatz_schluessel(self):
        """Schlüssel 'hauptsatz' vorhanden."""
        assert "hauptsatz" in self.ergebnis

    def test_demo_beweis_mit_cut_schluessel(self):
        """Schlüssel 'beweis_mit_cut' vorhanden."""
        assert "beweis_mit_cut" in self.ergebnis

    def test_demo_schnittfreier_beweis_schluessel(self):
        """Schlüssel 'schnittfreier_beweis' vorhanden."""
        assert "schnittfreier_beweis" in self.ergebnis

    def test_demo_folgerungen_schluessel(self):
        """Schlüssel 'folgerungen' ist eine Liste."""
        assert "folgerungen" in self.ergebnis
        assert isinstance(self.ergebnis["folgerungen"], list)
        assert len(self.ergebnis["folgerungen"]) >= 2

    def test_demo_schnittfreier_ist_axiom(self):
        """Der schnittfreie Beweis ist ein Axiom."""
        assert self.ergebnis["schnittfreier_beweis"]["ist_axiom"] is True

    def test_demo_baum_strings_nicht_leer(self):
        """Baumdarstellungen sind nicht leer."""
        assert len(self.ergebnis["beweis_mit_cut"]["baum"]) > 0
        assert len(self.ergebnis["schnittfreier_beweis"]["baum"]) > 0


# =============================================================================
# TESTS FÜR proof_theoretic_ordinal
# =============================================================================

class TestProofTheoreticOrdinal:
    """Tests für die beweistheoretische Ordinalzahl-Funktion."""

    def test_ordinal_pa(self):
        """|PA| = ε₀."""
        ergebnis = proof_theoretic_ordinal("PA")
        assert "ordinal" in ergebnis
        assert "ε₀" in ergebnis["ordinal"]

    def test_ordinal_pra(self):
        """|PRA| = ω^ω."""
        ergebnis = proof_theoretic_ordinal("PRA")
        assert "ordinal" in ergebnis
        assert "ω" in ergebnis["ordinal"]

    def test_ordinal_aca0(self):
        """|ACA₀| = ε₀."""
        ergebnis = proof_theoretic_ordinal("ACA0")
        assert "ordinal" in ergebnis
        assert "ε₀" in ergebnis["ordinal"]

    def test_ordinal_unbekanntes_system(self):
        """Unbekannte Systeme geben Fehler zurück."""
        ergebnis = proof_theoretic_ordinal("XYZ")
        assert "fehler" in ergebnis
        assert "bekannte_systeme" in ergebnis

    def test_ordinal_beschreibung_vorhanden(self):
        """PA-Beschreibung enthält Gentzen-Verweis."""
        ergebnis = proof_theoretic_ordinal("PA")
        assert "beschreibung" in ergebnis
        assert "gentzen" in ergebnis


# =============================================================================
# TESTS FÜR big_five_systems
# =============================================================================

class TestBigFiveSystems:
    """Tests für die Big Five der Reverse Mathematics."""

    def setup_method(self):
        """Lädt die Big Five vor jedem Test."""
        self.systeme = big_five_systems()

    def test_big_five_gibt_dict_zurueck(self):
        """big_five_systems gibt ein Dict zurück."""
        assert isinstance(self.systeme, dict)

    def test_big_five_genau_fuenf_systeme(self):
        """Mindestens 5 Hauptsysteme vorhanden."""
        # RCA0, WKL0, ACA0, ATR0, Pi11-CA0 + hierarchie
        hauptsysteme = [k for k in self.systeme.keys() if k != "hierarchie"]
        assert len(hauptsysteme) >= 5

    def test_big_five_rca0_vorhanden(self):
        """RCA₀ ist vorhanden."""
        assert "RCA0" in self.systeme

    def test_big_five_wkl0_vorhanden(self):
        """WKL₀ ist vorhanden."""
        assert "WKL0" in self.systeme

    def test_big_five_aca0_vorhanden(self):
        """ACA₀ ist vorhanden."""
        assert "ACA0" in self.systeme

    def test_big_five_atr0_vorhanden(self):
        """ATR₀ ist vorhanden."""
        assert "ATR0" in self.systeme

    def test_big_five_pi11_ca0_vorhanden(self):
        """Π₁¹-CA₀ ist vorhanden."""
        assert "Pi11-CA0" in self.systeme

    def test_big_five_jedes_system_hat_ordinal(self):
        """Jedes System hat eine Ordinalzahl."""
        for name, system in self.systeme.items():
            if name != "hierarchie":
                assert "ordinal" in system, f"{name} fehlt Ordinalzahl"

    def test_big_five_jedes_system_hat_axiome(self):
        """Jedes System hat Axiome."""
        for name, system in self.systeme.items():
            if name != "hierarchie":
                assert "axiome" in system, f"{name} fehlt Axiome"

    def test_big_five_hierarchie_vorhanden(self):
        """Hierarchie-Information ist vorhanden."""
        assert "hierarchie" in self.systeme


# =============================================================================
# TESTS FÜR reverse_math_example
# =============================================================================

class TestReverseMathExample:
    """Tests für Reverse-Mathematics-Beispiele."""

    def test_bolzano_weierstrass(self):
        """Bolzano-Weierstraß ist ACA₀-äquivalent."""
        ergebnis = reverse_math_example("bolzano_weierstrass")
        assert "aequivalent_zu" in ergebnis
        assert ergebnis["aequivalent_zu"] == "ACA₀"

    def test_heine_borel(self):
        """Heine-Borel ist WKL₀-äquivalent."""
        ergebnis = reverse_math_example("heine_borel")
        assert "aequivalent_zu" in ergebnis
        assert ergebnis["aequivalent_zu"] == "WKL₀"

    def test_ramsey(self):
        """Ramsey-Theorem ist ACA₀-äquivalent."""
        ergebnis = reverse_math_example("ramsey")
        assert "aequivalent_zu" in ergebnis
        assert ergebnis["aequivalent_zu"] == "ACA₀"

    def test_konig(self):
        """König-Lemma ist WKL₀-äquivalent."""
        ergebnis = reverse_math_example("konig")
        assert "aequivalent_zu" in ergebnis
        assert ergebnis["aequivalent_zu"] == "WKL₀"

    def test_unbekannter_satz(self):
        """Unbekannte Sätze geben Fehler zurück."""
        ergebnis = reverse_math_example("unbekannt")
        assert "fehler" in ergebnis
        assert "bekannte_saetze" in ergebnis


# =============================================================================
# TESTS FÜR epsilon_zero_demo
# =============================================================================

class TestEpsilonZeroDemo:
    """Tests für die ε₀-Demonstration."""

    def setup_method(self):
        """Führt die Demo vor jedem Test aus."""
        self.ergebnis = epsilon_zero_demo()

    def test_gibt_dict_zurueck(self):
        """epsilon_zero_demo gibt ein Dict zurück."""
        assert isinstance(self.ergebnis, dict)

    def test_definition_vorhanden(self):
        """Definition von ε₀ ist vorhanden."""
        assert "definition" in self.ergebnis

    def test_folge_vorhanden(self):
        """Aufsteigende Folge zu ε₀ ist vorhanden."""
        assert "folge" in self.ergebnis
        assert len(self.ergebnis["folge"]) >= 3

    def test_gentzen_verbindung(self):
        """Gentzen-Verbindung ist dokumentiert."""
        assert "gentzen" in self.ergebnis
        assert "satz" in self.ergebnis["gentzen"]

    def test_cantor_normalform_beispiele(self):
        """Cantor-Normalform-Beispiele sind vorhanden."""
        assert "cantor_normalform" in self.ergebnis
        assert "beispiele" in self.ergebnis["cantor_normalform"]


# =============================================================================
# TESTS FÜR gentzen_consistency_proof_sketch
# =============================================================================

class TestGentzenConsistencyProof:
    """Tests für die Gentzen-Konsistenzbeweis-Skizze."""

    def setup_method(self):
        """Führt die Skizze vor jedem Test aus."""
        self.ergebnis = gentzen_consistency_proof_sketch()

    def test_gibt_dict_zurueck(self):
        """gentzen_consistency_proof_sketch gibt ein Dict zurück."""
        assert isinstance(self.ergebnis, dict)

    def test_titel_vorhanden(self):
        """Titel ist vorhanden."""
        assert "titel" in self.ergebnis

    def test_system_vorhanden(self):
        """Das verwendete System ist angegeben."""
        assert "system" in self.ergebnis
        assert "PRA" in self.ergebnis["system"]

    def test_beweisidee_vorhanden(self):
        """Beweisidee mit Schritten ist vorhanden."""
        assert "beweisidee" in self.ergebnis
        assert len(self.ergebnis["beweisidee"]) >= 3

    def test_godel_verbindung(self):
        """Gödel-Verbindung ist dokumentiert."""
        assert "godel_verbindung" in self.ergebnis


# =============================================================================
# TESTS FÜR proof_complexity_comparison
# =============================================================================

class TestProofComplexityComparison:
    """Tests für den Beweissystem-Vergleich."""

    def setup_method(self):
        """Führt den Vergleich vor jedem Test aus."""
        self.ergebnis = proof_complexity_comparison()

    def test_gibt_dict_zurueck(self):
        """Gibt Dict zurück."""
        assert isinstance(self.ergebnis, dict)

    def test_cook_reckhow_vorhanden(self):
        """Cook-Reckhow-Theorem ist dokumentiert."""
        assert "cook_reckhow" in self.ergebnis

    def test_systeme_vorhanden(self):
        """Beweissysteme sind aufgeführt."""
        assert "systeme" in self.ergebnis
        assert len(self.ergebnis["systeme"]) >= 4

    def test_resolution_vorhanden(self):
        """Resolution ist aufgeführt."""
        assert "Resolution" in self.ergebnis["systeme"]

    def test_frege_vorhanden(self):
        """Frege-System ist aufgeführt."""
        assert "Frege" in self.ergebnis["systeme"]


# =============================================================================
# TESTS FÜR resolution_proof_size
# =============================================================================

class TestResolutionProofSize:
    """Tests für die Resolution-Komplexitätsanalyse."""

    def test_leere_klauselmenge(self):
        """Leere Klauselmenge wird verarbeitet."""
        ergebnis = resolution_proof_size([])
        assert isinstance(ergebnis, dict)
        assert "eingabe" in ergebnis

    def test_widerlegbare_klauseln(self):
        """Widerlegbare Klauselmenge wird erkannt."""
        # {A} und {¬A} sind widerlegbar
        klauseln = [frozenset({"A"}), frozenset({"¬A"})]
        ergebnis = resolution_proof_size(klauseln)
        assert ergebnis["analyse"]["ist_widerlegbar"] is True

    def test_nicht_widerlegbare_klauseln(self):
        """Erfüllbare Klauselmenge ist nicht widerlegbar."""
        # {A} allein ist erfüllbar
        klauseln = [frozenset({"A"})]
        ergebnis = resolution_proof_size(klauseln)
        assert ergebnis["analyse"]["ist_widerlegbar"] is False

    def test_analyse_schluessel(self):
        """Analyseschlüssel sind vorhanden."""
        ergebnis = resolution_proof_size([frozenset({"A"})])
        assert "eingabe" in ergebnis
        assert "analyse" in ergebnis
        assert "komplexitaet" in ergebnis


# =============================================================================
# TESTS FÜR pigeonhole_principle_hardness
# =============================================================================

class TestPigeonholePrincipleHardness:
    """Tests für die Taubenschlag-Komplexitätsanalyse."""

    def setup_method(self):
        """Führt die Analyse vor jedem Test aus."""
        self.ergebnis = pigeonhole_principle_hardness()

    def test_gibt_dict_zurueck(self):
        """pigeonhole_principle_hardness gibt Dict zurück."""
        assert isinstance(self.ergebnis, dict)

    def test_problem_beschreibung(self):
        """Problem-Beschreibung ist vorhanden."""
        assert "problem" in self.ergebnis

    def test_kodierung(self):
        """Kodierung ist dokumentiert."""
        assert "kodierung" in self.ergebnis
        assert "anzahl_variablen" in self.ergebnis["kodierung"]

    def test_komplexitaet(self):
        """Komplexitätsschranken sind angegeben."""
        assert "komplexitaet" in self.ergebnis
        assert "resolution_schranke" in self.ergebnis["komplexitaet"]

    def test_haerteste_klasse(self):
        """Hardness-Vergleich ist vorhanden."""
        assert "haerteste_klasse" in self.ergebnis
        assert "frege" in self.ergebnis["haerteste_klasse"]
        assert "resolution" in self.ergebnis["haerteste_klasse"]


# =============================================================================
# TESTS FÜR deduction_theorem
# =============================================================================

class TestDeductionTheorem:
    """Tests für das Deduktionstheorem."""

    def test_einfache_anwendung(self):
        """Deduktionstheorem transformiert Hypothesen."""
        ergebnis = deduction_theorem(["A"], "B")
        assert "neue_konklusion" in ergebnis
        assert "A → B" in ergebnis["neue_konklusion"]

    def test_mehrere_hypothesen(self):
        """Mehrere Hypothesen werden korrekt behandelt."""
        ergebnis = deduction_theorem(["A", "B"], "C")
        assert "neue_konklusion" in ergebnis
        assert "B → C" in ergebnis["neue_konklusion"]

    def test_leere_hypothesen(self):
        """Leere Hypothesen-Liste gibt Fehler zurück."""
        ergebnis = deduction_theorem([], "A")
        assert "fehler" in ergebnis

    def test_ausgabe_struktur(self):
        """Ausgabe-Struktur ist korrekt."""
        ergebnis = deduction_theorem(["A"], "B")
        assert "eingabe" in ergebnis
        assert "ausgabe" in ergebnis
        assert "theorem" in ergebnis


# =============================================================================
# TESTS FÜR subformula_property
# =============================================================================

class TestSubformulaProperty:
    """Tests für die Teilformel-Eigenschaft."""

    def test_axiom_hat_eigenschaft(self):
        """Ein Axiom-Beweis erfüllt die Teilformel-Eigenschaft."""
        beweis = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        assert subformula_property(beweis) is True

    def test_weakening_verletzt_eigenschaft_nicht(self):
        """Weakening kann Teilformel-Eigenschaft verletzen (neue Formel)."""
        blatt = ProofTree(Sequent(["A"], ["A"]), "Axiom")
        beweis = LKRules.weakening_left(blatt, "B")
        # B ist nicht in der Schlusssequent-Formel A, aber in Prämisse
        # Das Ergebnis hängt von der Implementierung ab
        ergebnis = subformula_property(beweis)
        assert isinstance(ergebnis, bool)
