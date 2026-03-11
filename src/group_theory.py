"""
@file group_theory.py
@brief Gruppentheorie: Endliche Gruppen, Homomorphismen, Permutationsgruppen,
       Sylow-Sätze, Burnside-Lemma, Klassifikation abelscher Gruppen.
@description
    Dieses Modul implementiert die grundlegenden Strukturen und Algorithmen
    der Gruppentheorie:

    Klassen:
    - Group               – Allgemeine endliche Gruppe (Cayley-Tabelle)
    - Subgroup            – Untergruppe H ≤ G
    - GroupHomomorphism   – Homomorphismus φ: G → H
    - QuotientGroup       – Faktorgruppe G/N für N ◁ G
    - PermutationGroup    – Untergruppe von S_n

    Freie Funktionen:
    - lagrange_theorem_check()   – Lagrange-Satz und Nebenklassen
    - sylow_theorems()           – Sylow-Sätze für eine Primzahl p
    - classify_abelian_group()   – Hauptsatz über endliche abelsche Gruppen
    - is_simple_group()          – Einfachheit einer Gruppe
    - group_action()             – Gruppenoperation, Orbits, Burnside
    - direct_product()           – Direktes Produkt G × H
    - semidirect_product()       – Halbdirektes Produkt N ⋊_φ H
    - symmetric_group()          – S_n (alle Permutationen)
    - alternating_group()        – A_n (gerade Permutationen)
    - cyclic_group()             – ℤ_n (zyklische Gruppe)
    - dihedral_group()           – D_n (Diedergruppe, Ordnung 2n)
    - quaternion_group()         – Q_8 (Quaternionengruppe)

    Mathematische Grundlagen:
    - Cayley-Tabelle definiert die Gruppenoperation vollständig
    - Lagrange: |G| = |H| · [G:H] für H ≤ G
    - Sylow: Für p^k | |G| existieren Sylow-p-Untergruppen der Ordnung p^k
    - Burnside: |X/G| = (1/|G|) · Σ_{g∈G} |Fix(g)|

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import math
import itertools
from typing import Any, Callable, Optional
from collections import defaultdict

import numpy as np
import sys
import os

# Pfad für eigene Module setzen
sys.path.insert(0, os.path.dirname(__file__))
from algebra_core import gcd
from algebra_numbertheory import is_prime, prime_factorization


# =============================================================================
# HILFS-FUNKTIONEN
# =============================================================================

def _factorint(n: int) -> dict[int, int]:
    """
    @brief Primfaktorzerlegung von n als Dictionary {p: e}.
    @param n Positive ganze Zahl
    @return Dictionary mit Primfaktoren als Schlüssel und Exponenten als Werte
    @lastModified 2026-03-10
    """
    # prime_factorization() gibt bereits {p: e} zurück
    return prime_factorization(n)


def _divisors(n: int) -> list[int]:
    """
    @brief Alle Teiler von n (sortiert aufsteigend).
    @param n Positive ganze Zahl
    @return Sortierte Liste aller Teiler
    @lastModified 2026-03-10
    """
    divs = []
    for i in range(1, int(math.isqrt(n)) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)


# =============================================================================
# KLASSE: Group
# =============================================================================

class Group:
    """
    @brief Endliche Gruppe G = (elements, operation).
    @description
        Repräsentiert eine endliche Gruppe durch ihre Elementmenge und
        eine Binäroperation, die als Cayley-Tabelle oder Lambda gespeichert wird.

        Eine Gruppe (G, ·) erfüllt:
        1. Abgeschlossenheit: ∀ a,b ∈ G: a·b ∈ G
        2. Assoziativität:    ∀ a,b,c ∈ G: (a·b)·c = a·(b·c)
        3. Neutrales Element: ∃ e ∈ G: e·a = a·e = a
        4. Inverses Element:  ∀ a ∈ G: ∃ a⁻¹ ∈ G: a·a⁻¹ = a⁻¹·a = e

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        elements: list[Any],
        operation: Callable[[Any, Any], Any],
        identity: Optional[Any] = None,
        name: str = "G"
    ) -> None:
        """
        @brief Konstruktor der Gruppe.
        @param elements Liste aller Gruppenelemente
        @param operation Binäre Operation (a, b) → a·b
        @param identity Neutrales Element (wird automatisch ermittelt, wenn None)
        @param name Bezeichner der Gruppe (für Ausgaben)
        @lastModified 2026-03-10
        """
        self.elements = list(elements)
        self.operation = operation
        self.name = name

        # Neutrales Element bestimmen
        if identity is not None:
            self.identity = identity
        else:
            self.identity = self._find_identity()

        # Interne Cayley-Tabelle als Dictionary cachen
        self._table: Optional[dict] = None

    def _find_identity(self) -> Any:
        """
        @brief Sucht das neutrale Element e mit e·a = a·e = a für alle a.
        @return Das neutrale Element
        @raises ValueError Wenn kein neutrales Element gefunden wird
        @lastModified 2026-03-10
        """
        for e in self.elements:
            # Prüfen ob e für alle Elemente das neutrale Element ist
            if all(
                self.operation(e, a) == a and self.operation(a, e) == a
                for a in self.elements
            ):
                return e
        raise ValueError(f"Kein neutrales Element in {self.name} gefunden!")

    def _build_table(self) -> dict:
        """
        @brief Erstellt die vollständige Cayley-Tabelle als Dictionary.
        @return Dictionary (a, b) → a·b für alle a, b ∈ G
        @lastModified 2026-03-10
        """
        if self._table is None:
            self._table = {}
            for a in self.elements:
                for b in self.elements:
                    self._table[(a, b)] = self.operation(a, b)
        return self._table

    def op(self, a: Any, b: Any) -> Any:
        """
        @brief Wendet die Gruppenoperation an: a · b.
        @param a Erstes Element
        @param b Zweites Element
        @return Produkt a · b
        @lastModified 2026-03-10
        """
        return self.operation(a, b)

    def order(self) -> int:
        """
        @brief Gibt die Gruppenordnung |G| zurück.
        @return Anzahl der Elemente in G
        @lastModified 2026-03-10
        """
        return len(self.elements)

    def element_order(self, a: Any) -> int:
        """
        @brief Berechnet die Ordnung ord(a) = min{n ≥ 1 : a^n = e}.
        @param a Gruppenelement
        @return Ordnung des Elements
        @lastModified 2026-03-10
        """
        current = a
        for n in range(1, len(self.elements) + 1):
            if current == self.identity:
                return n
            current = self.operation(current, a)
        raise ValueError(f"Element {a} hat keine endliche Ordnung in {self.name}!")

    def inverse(self, a: Any) -> Any:
        """
        @brief Berechnet das inverse Element a⁻¹ mit a · a⁻¹ = e.
        @param a Gruppenelement
        @return Inverses Element a⁻¹
        @raises ValueError Wenn kein Inverses gefunden wird
        @lastModified 2026-03-10
        """
        e = self.identity
        for b in self.elements:
            if self.operation(a, b) == e and self.operation(b, a) == e:
                return b
        raise ValueError(f"Kein Inverses für {a} in {self.name} gefunden!")

    def power(self, a: Any, n: int) -> Any:
        """
        @brief Berechnet a^n (n-faches Produkt).
        @param a Gruppenelement
        @param n Exponent (kann negativ sein)
        @return a^n
        @lastModified 2026-03-10
        """
        if n == 0:
            return self.identity
        if n < 0:
            # Negatives Potenzieren: (a⁻¹)^|n|
            a = self.inverse(a)
            n = -n
        # Schnelle Exponentiation via wiederholtem Quadrieren
        result = self.identity
        base = a
        while n > 0:
            if n % 2 == 1:
                result = self.operation(result, base)
            base = self.operation(base, base)
            n //= 2
        return result

    def is_abelian(self) -> bool:
        """
        @brief Prüft ob G kommutativ (abelsch) ist: ∀ a,b: a·b = b·a.
        @return True wenn G abelsch ist
        @lastModified 2026-03-10
        """
        for a in self.elements:
            for b in self.elements:
                if self.operation(a, b) != self.operation(b, a):
                    return False
        return True

    def subgroups(self) -> list["Subgroup"]:
        """
        @brief Berechnet alle Untergruppen H ≤ G.
        @description
            Eine Teilmenge H ⊆ G ist Untergruppe wenn:
            1. e ∈ H (neutrales Element)
            2. a,b ∈ H → a·b ∈ H (Abgeschlossenheit)
            3. a ∈ H → a⁻¹ ∈ H (Inversenabgeschlossenheit)

            Aus dem Satz von Lagrange: |H| teilt |G|.
        @return Liste aller Untergruppen
        @lastModified 2026-03-10
        """
        n = len(self.elements)
        result = []
        e = self.identity

        # Alle Teilmengen mit neutralem Element und Lagrange-konformer Größe
        for size in _divisors(n):
            # Alle Teilmengen der Größe 'size' die e enthalten
            rest = [x for x in self.elements if x != e]
            for combo in itertools.combinations(rest, size - 1):
                candidate = [e] + list(combo)
                if self._is_subgroup(candidate):
                    # Als Subgroup-Objekt einfügen, falls nicht bereits vorhanden
                    cand_set = set(candidate)
                    already = any(
                        set(sg.elements) == cand_set for sg in result
                    )
                    if not already:
                        sub_op = lambda a, b, c=candidate: self.operation(a, b)
                        result.append(Subgroup(candidate, self.operation, self, e))
        return result

    def _is_subgroup(self, candidate: list[Any]) -> bool:
        """
        @brief Prüft ob eine Teilmenge eine gültige Untergruppe bildet.
        @param candidate Liste der Kandidaten-Elemente
        @return True wenn Untergruppe
        @lastModified 2026-03-10
        """
        cand_set = set(candidate)
        # Neutrales Element muss enthalten sein
        if self.identity not in cand_set:
            return False
        # Abgeschlossenheit unter Operation und Inversenbildung
        for a in candidate:
            for b in candidate:
                if self.operation(a, b) not in cand_set:
                    return False
        return True

    def normal_subgroups(self) -> list["Subgroup"]:
        """
        @brief Berechnet alle Normalteiler N ◁ G.
        @description
            N ist Normalteiler wenn: ∀ g ∈ G, n ∈ N: g·n·g⁻¹ ∈ N
            Äquivalent: gN = Ng für alle g ∈ G (links = rechts Nebenklassen)
        @return Liste aller Normalteiler
        @lastModified 2026-03-10
        """
        all_subs = self.subgroups()
        normal = []
        for H in all_subs:
            if self._is_normal(H.elements):
                normal.append(H)
        return normal

    def _is_normal(self, h_elements: list[Any]) -> bool:
        """
        @brief Prüft ob eine Untergruppe ein Normalteiler ist.
        @param h_elements Elemente der potenziellen Untergruppe
        @return True wenn Normalteiler
        @lastModified 2026-03-10
        """
        h_set = set(h_elements)
        for g in self.elements:
            g_inv = self.inverse(g)
            for h in h_elements:
                # Konjugiertes Element: g·h·g⁻¹ muss in H liegen
                conj = self.operation(self.operation(g, h), g_inv)
                if conj not in h_set:
                    return False
        return True

    def center(self) -> "Subgroup":
        """
        @brief Berechnet das Zentrum Z(G) = {z ∈ G | ∀g ∈ G: zg = gz}.
        @return Zentrum als Untergruppe
        @lastModified 2026-03-10
        """
        center_elements = []
        for z in self.elements:
            # z ist im Zentrum wenn es mit allen Elementen kommutiert
            if all(
                self.operation(z, g) == self.operation(g, z)
                for g in self.elements
            ):
                center_elements.append(z)
        return Subgroup(center_elements, self.operation, self, self.identity)

    def commutator_subgroup(self) -> "Subgroup":
        """
        @brief Berechnet die Kommutatorgruppe [G,G] = ⟨[a,b] | a,b ∈ G⟩.
        @description
            Kommutatoren: [a,b] = a⁻¹·b⁻¹·a·b
            Die Kommutatorgruppe ist der kleinste Normalteiler N mit G/N abelsch.
        @return Kommutatoruntergruppe [G,G]
        @lastModified 2026-03-10
        """
        # Alle Kommutatoren [a,b] = a⁻¹ b⁻¹ a b berechnen
        commutators = set()
        for a in self.elements:
            for b in self.elements:
                a_inv = self.inverse(a)
                b_inv = self.inverse(b)
                # [a,b] = a⁻¹ b⁻¹ a b
                comm = self.operation(
                    self.operation(a_inv, b_inv),
                    self.operation(a, b)
                )
                commutators.add(comm)

        # Abschluss unter der Gruppenoperation bilden (erzeugter Abschluss)
        generated = set(commutators)
        generated.add(self.identity)

        changed = True
        while changed:
            changed = False
            new_elements = set()
            for a in generated:
                for b in generated:
                    prod = self.operation(a, b)
                    if prod not in generated:
                        new_elements.add(prod)
                        changed = True
            generated.update(new_elements)

        gen_list = sorted(generated, key=lambda x: self.elements.index(x))
        return Subgroup(gen_list, self.operation, self, self.identity)

    def is_cyclic(self) -> bool:
        """
        @brief Prüft ob G zyklisch ist: G = ⟨g⟩ für ein g ∈ G.
        @description
            G ist zyklisch wenn ein Element der Ordnung |G| existiert.
        @return True wenn G zyklisch
        @lastModified 2026-03-10
        """
        n = len(self.elements)
        for g in self.elements:
            if self.element_order(g) == n:
                return True
        return False

    def generators(self) -> list[Any]:
        """
        @brief Findet minimale Erzeugendensysteme von G.
        @description
            Ein Element g ist Generator wenn ⟨g⟩ = G (für zyklische Gruppen).
            Für nicht-zyklische Gruppen: Menge von Elementen, die G erzeugen.
        @return Liste der Generatoren (Elemente der Ordnung |G| oder minimales Erzeugungssystem)
        @lastModified 2026-03-10
        """
        n = len(self.elements)
        # Zuerst Elemente der vollen Ordnung suchen (zyklisch erzeugend)
        full_gens = [g for g in self.elements if self.element_order(g) == n]
        if full_gens:
            return full_gens

        # Minimales Erzeugendensystem für nicht-zyklische Gruppen
        # Einfacher Ansatz: Alle Nicht-Identitäts-Elemente
        return [g for g in self.elements if g != self.identity]

    def cayley_table(self) -> np.ndarray:
        """
        @brief Erstellt die Cayley-Tabelle als numpy-Array.
        @description
            Zeile i, Spalte j: Index des Produkts elements[i] · elements[j]
            in der Elementliste.
        @return 2D numpy-Array mit Indizes der Produkte
        @lastModified 2026-03-10
        """
        n = len(self.elements)
        # Lookup-Dictionary: Element → Index
        index = {e: i for i, e in enumerate(self.elements)}
        table = np.zeros((n, n), dtype=int)
        for i, a in enumerate(self.elements):
            for j, b in enumerate(self.elements):
                prod = self.operation(a, b)
                table[i, j] = index[prod]
        return table

    def cosets(self, H: "Subgroup", left: bool = True) -> list[list[Any]]:
        """
        @brief Berechnet alle Links- oder Rechtsnebenklassen von H in G.
        @param H Untergruppe
        @param left True für Linksnebenklassen gH, False für Rechtsnebenklassen Hg
        @return Liste der Nebenklassen
        @lastModified 2026-03-10
        """
        h_set = set(H.elements)
        covered = set()
        coset_list = []

        for g in self.elements:
            if g not in covered:
                if left:
                    # Linksnebenklasse gH = {g·h | h ∈ H}
                    coset = [self.operation(g, h) for h in H.elements]
                else:
                    # Rechtsnebenklasse Hg = {h·g | h ∈ H}
                    coset = [self.operation(h, g) for h in H.elements]
                coset_list.append(coset)
                covered.update(coset)

        return coset_list

    def __repr__(self) -> str:
        return f"Group({self.name}, order={self.order()})"

    def __len__(self) -> int:
        return self.order()


# =============================================================================
# KLASSE: Subgroup
# =============================================================================

class Subgroup(Group):
    """
    @brief Untergruppe H ≤ G mit Einbettung in die Obergruppe.
    @description
        Eine Untergruppe erbt die Gruppenstruktur und enthält zusätzlich
        eine Referenz auf die Obergruppe G sowie die Einbettungsabbildung.

        Kriterium (Untergruppen-Test):
        H ≤ G ⟺ H ≠ ∅ und ∀ a,b ∈ H: a·b⁻¹ ∈ H

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        elements: list[Any],
        operation: Callable[[Any, Any], Any],
        parent: Group,
        identity: Optional[Any] = None,
        name: str = "H"
    ) -> None:
        """
        @brief Konstruktor der Untergruppe.
        @param elements Elemente der Untergruppe
        @param operation Gruppenoperation (von der Obergruppe geerbt)
        @param parent Obergruppe G
        @param identity Neutrales Element
        @param name Bezeichner
        @lastModified 2026-03-10
        """
        super().__init__(elements, operation, identity, name)
        self.parent = parent

    def index(self) -> int:
        """
        @brief Berechnet den Index [G:H] = |G|/|H|.
        @description
            Nach dem Satz von Lagrange gilt: |G| = [G:H] · |H|
        @return Index der Untergruppe in der Obergruppe
        @lastModified 2026-03-10
        """
        return self.parent.order() // self.order()

    def __repr__(self) -> str:
        return f"Subgroup({self.name}, order={self.order()}, parent={self.parent.name})"


# =============================================================================
# KLASSE: GroupHomomorphism
# =============================================================================

class GroupHomomorphism:
    """
    @brief Gruppenhomomorphismus φ: G → H mit φ(ab) = φ(a)φ(b).
    @description
        Ein Gruppenhomomorphismus ist eine strukturerhaltende Abbildung
        zwischen zwei Gruppen. Es gilt:
        - φ(e_G) = e_H  (Neutralelemente werden abgebildet)
        - φ(a⁻¹) = φ(a)⁻¹  (Inverse werden abgebildet)
        - kern(φ) = {g ∈ G | φ(g) = e_H} ist Normalteiler in G
        - bild(φ) = {φ(g) | g ∈ G} ist Untergruppe von H

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        source: Group,
        target: Group,
        mapping: Callable[[Any], Any],
        name: str = "φ"
    ) -> None:
        """
        @brief Konstruktor des Homomorphismus.
        @param source Quellgruppe G
        @param target Zielgruppe H
        @param mapping Abbildungsfunktion g → φ(g)
        @param name Bezeichner des Homomorphismus
        @raises ValueError Wenn die Homomorphismus-Eigenschaft verletzt ist
        @lastModified 2026-03-10
        """
        self.source = source
        self.target = target
        self.mapping = mapping
        self.name = name

        # Homomorphismus-Eigenschaft überprüfen
        if not self._verify():
            raise ValueError(
                f"Abbildung {name} ist kein Homomorphismus: "
                f"φ(ab) ≠ φ(a)φ(b) für manche a,b ∈ G"
            )

    def _verify(self) -> bool:
        """
        @brief Überprüft die Homomorphismus-Eigenschaft φ(ab) = φ(a)φ(b).
        @return True wenn Homomorphismus-Eigenschaft erfüllt
        @lastModified 2026-03-10
        """
        G = self.source
        H = self.target
        for a in G.elements:
            for b in G.elements:
                # φ(a·b) muss gleich φ(a)·φ(b) sein
                lhs = self.mapping(G.operation(a, b))
                rhs = H.operation(self.mapping(a), self.mapping(b))
                if lhs != rhs:
                    return False
        return True

    def __call__(self, g: Any) -> Any:
        """
        @brief Wertet φ(g) aus.
        @param g Element aus G
        @return Bild φ(g) in H
        @lastModified 2026-03-10
        """
        return self.mapping(g)

    def kernel(self) -> Subgroup:
        """
        @brief Berechnet den Kern ker(φ) = {g ∈ G | φ(g) = e_H}.
        @description
            Der Kern ist immer ein Normalteiler in G.
            G/ker(φ) ≅ bild(φ) (Homomorphiesatz)
        @return Kern als Untergruppe von G
        @lastModified 2026-03-10
        """
        e_H = self.target.identity
        kern_elements = [g for g in self.source.elements if self.mapping(g) == e_H]
        return Subgroup(kern_elements, self.source.operation, self.source,
                        self.source.identity, "ker(φ)")

    def image(self) -> Subgroup:
        """
        @brief Berechnet das Bild im(φ) = {φ(g) | g ∈ G}.
        @description
            Das Bild ist immer eine Untergruppe von H.
        @return Bild als Untergruppe von H
        @lastModified 2026-03-10
        """
        image_elements = list({self.mapping(g) for g in self.source.elements})
        return Subgroup(image_elements, self.target.operation, self.target,
                        self.target.identity, "im(φ)")

    def is_injective(self) -> bool:
        """
        @brief Prüft ob φ injektiv ist (ker(φ) = {e}).
        @description
            Ein Homomorphismus ist injektiv genau dann wenn ker(φ) = {e_G}.
        @return True wenn injektiv
        @lastModified 2026-03-10
        """
        return len(self.kernel().elements) == 1

    def is_surjective(self) -> bool:
        """
        @brief Prüft ob φ surjektiv ist (im(φ) = H).
        @return True wenn surjektiv
        @lastModified 2026-03-10
        """
        return len(self.image().elements) == len(self.target.elements)

    def is_isomorphism(self) -> bool:
        """
        @brief Prüft ob φ ein Isomorphismus ist (bijektiver Homomorphismus).
        @return True wenn Isomorphismus
        @lastModified 2026-03-10
        """
        return self.is_injective() and self.is_surjective()

    def __repr__(self) -> str:
        return (f"GroupHomomorphism({self.name}: {self.source.name} → "
                f"{self.target.name})")


# =============================================================================
# KLASSE: QuotientGroup
# =============================================================================

class QuotientGroup(Group):
    """
    @brief Faktorgruppe G/N für einen Normalteiler N ◁ G.
    @description
        Die Faktorgruppe besteht aus den Nebenklassen gN.
        Operation: (gN) · (hN) = (g·h)N

        Nach dem Homomorphiesatz gilt:
        G/ker(φ) ≅ im(φ) für jeden Homomorphismus φ: G → H

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, G: Group, N: Subgroup) -> None:
        """
        @brief Konstruktor der Faktorgruppe.
        @param G Übergruppe
        @param N Normalteiler N ◁ G
        @raises ValueError Wenn N kein Normalteiler ist
        @lastModified 2026-03-10
        """
        # Prüfen ob N wirklich Normalteiler ist
        if not G._is_normal(N.elements):
            raise ValueError(f"{N.name} ist kein Normalteiler von {G.name}!")

        self.G = G
        self.N = N

        # Nebenklassen als Frozensets (hashbar, für Operation als Dict-Key nutzbar)
        cosets = G.cosets(N, left=True)
        # Nebenklassen als Frozensets repräsentieren
        coset_frozensets = [frozenset(c) for c in cosets]

        # Nebenklassen-Operation: (gN)·(hN) = (gh)N
        def coset_op(c1: frozenset, c2: frozenset) -> frozenset:
            # Repräsentanten aus den Nebenklassen wählen
            g = next(iter(c1))
            h = next(iter(c2))
            # Produkt g·h
            gh = G.operation(g, h)
            # Nebenklasse von gh finden
            for cs in coset_frozensets:
                if gh in cs:
                    return cs
            raise RuntimeError(f"Nebenklasse für {gh} nicht gefunden!")

        # Neutrales Element: Nebenklasse von e (= N selbst)
        n_set = frozenset(N.elements)

        super().__init__(
            coset_frozensets,
            coset_op,
            identity=n_set,
            name=f"{G.name}/{N.name}"
        )

    def __repr__(self) -> str:
        return f"QuotientGroup({self.G.name}/{self.N.name}, order={self.order()})"


# =============================================================================
# KLASSE: PermutationGroup
# =============================================================================

class PermutationGroup(Group):
    """
    @brief Permutationsgruppe als Untergruppe von S_n.
    @description
        Elemente sind Tupel σ = (σ(0), σ(1), ..., σ(n-1)) (0-indiziert).
        Die Gruppenoperation ist die Komposition von Permutationen:
        (σ∘τ)(i) = σ(τ(i))

        Für eine Permutation σ ∈ S_n gilt:
        - Zyklentyp: Zerlegung in disjunkte Zykel
        - Vorzeichen sgn(σ) = (-1)^(Anzahl gerader Zykel)
        - Gerade Permutationen bilden A_n ≤ S_n

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        elements: list[tuple],
        n: int,
        name: str = "PermGroup"
    ) -> None:
        """
        @brief Konstruktor der Permutationsgruppe.
        @param elements Liste der Permutationen als Tupel (0-indiziert)
        @param n Grad der Permutationsgruppe
        @param name Bezeichner
        @lastModified 2026-03-10
        """
        self.n = n

        # Permutations-Komposition: (σ∘τ)(i) = σ(τ(i))
        def perm_op(sigma: tuple, tau: tuple) -> tuple:
            return tuple(sigma[tau[i]] for i in range(n))

        # Identitätspermutation: id(i) = i
        identity = tuple(range(n))

        super().__init__(elements, perm_op, identity, name)

    def cycle_decomposition(self, sigma: tuple) -> list[list[int]]:
        """
        @brief Zerlegt eine Permutation in disjunkte Zykel.
        @param sigma Permutation als Tupel
        @return Liste der Zykel (als Listen von Indizes)
        @lastModified 2026-03-10
        """
        visited = [False] * self.n
        cycles = []
        for start in range(self.n):
            if not visited[start]:
                cycle = []
                current = start
                while not visited[current]:
                    visited[current] = True
                    cycle.append(current)
                    current = sigma[current]
                if len(cycle) > 1 or len(cycles) == 0:
                    cycles.append(cycle)
                elif len(cycle) == 1:
                    # Fixpunkte (1-Zykel) separat behandeln
                    cycles.append(cycle)
        return cycles

    def cycle_type(self, sigma: tuple) -> tuple[int, ...]:
        """
        @brief Berechnet den Zyklentyp von σ (sortierte Zyklenlängen).
        @param sigma Permutation als Tupel
        @return Sortiertes Tupel der Zyklenlängen (absteigend)
        @lastModified 2026-03-10
        """
        cycles = self.cycle_decomposition(sigma)
        lengths = sorted([len(c) for c in cycles], reverse=True)
        return tuple(lengths)

    def sign(self, sigma: tuple) -> int:
        """
        @brief Berechnet das Vorzeichen (Signum) von σ.
        @description
            sgn(σ) = (-1)^(Anzahl der Transpositionen in der Zerlegung)
                   = (-1)^(n - Anzahl der Zykel)
            Gerade Permutation: sgn = +1, ungerade: sgn = -1
        @param sigma Permutation als Tupel
        @return +1 oder -1
        @lastModified 2026-03-10
        """
        cycles = self.cycle_decomposition(sigma)
        # Anzahl Transpositionen = Summe (Zyklenlänge - 1)
        num_transpositions = sum(len(c) - 1 for c in cycles)
        return (-1) ** num_transpositions

    def is_even_permutation(self, sigma: tuple) -> bool:
        """
        @brief Prüft ob σ eine gerade Permutation ist (sgn(σ) = +1).
        @param sigma Permutation als Tupel
        @return True wenn gerade Permutation
        @lastModified 2026-03-10
        """
        return self.sign(sigma) == 1

    def orbit(self, element: int) -> set[int]:
        """
        @brief Berechnet die Bahn (Orbit) eines Elements unter G.
        @description
            Orb_G(x) = {σ(x) | σ ∈ G}
        @param element Element aus {0, ..., n-1}
        @return Menge aller Bilder von element unter G
        @lastModified 2026-03-10
        """
        return {sigma[element] for sigma in self.elements}

    def stabilizer(self, element: int) -> "PermutationGroup":
        """
        @brief Berechnet den Stabilisator Stab_G(x) = {σ ∈ G | σ(x) = x}.
        @description
            Bahnformel: |Orb(x)| · |Stab(x)| = |G|
        @param element Zu stabilisierendes Element
        @return Stabilisator als Permutationsgruppe
        @lastModified 2026-03-10
        """
        stab_elements = [sigma for sigma in self.elements if sigma[element] == element]
        return PermutationGroup(stab_elements, self.n, f"Stab({element})")

    def __repr__(self) -> str:
        return f"PermutationGroup({self.name}, n={self.n}, order={self.order()})"


# =============================================================================
# FREIE FUNKTIONEN
# =============================================================================

def lagrange_theorem_check(G: Group, H: Subgroup) -> dict:
    """
    @brief Überprüft den Satz von Lagrange: |G| = |H| · [G:H].
    @description
        Satz von Lagrange: Für jede Untergruppe H ≤ G gilt:
        |G| = |H| · [G:H]
        wobei [G:H] = |G|/|H| der Index von H in G ist.

        Die Funktion berechnet alle Linksnebenklassen gH und verifiziert,
        dass sie eine Partition von G bilden.

    @param G Endliche Gruppe
    @param H Untergruppe von G
    @return Dictionary mit Ergebnissen:
            - 'order_G': |G|
            - 'order_H': |H|
            - 'index': [G:H]
            - 'left_cosets': Liste der Linksnebenklassen
            - 'right_cosets': Liste der Rechtsnebenklassen
            - 'lagrange_holds': True wenn |G| = |H| · [G:H]
            - 'is_normal': True wenn H Normalteiler ist
    @lastModified 2026-03-10
    """
    order_G = G.order()
    order_H = H.order()
    index = order_G // order_H

    # Links- und Rechtsnebenklassen berechnen
    left_cosets = G.cosets(H, left=True)
    right_cosets = G.cosets(H, left=False)

    # Lagrange-Bedingung prüfen
    lagrange_holds = (order_G == order_H * index) and (len(left_cosets) == index)

    # Normalteiler-Test: Links- = Rechtsnebenklassen?
    left_sets = [frozenset(c) for c in left_cosets]
    right_sets = [frozenset(c) for c in right_cosets]
    is_normal = set(left_sets) == set(right_sets)

    return {
        'order_G': order_G,
        'order_H': order_H,
        'index': index,
        'left_cosets': left_cosets,
        'right_cosets': right_cosets,
        'lagrange_holds': lagrange_holds,
        'is_normal': is_normal,
        'num_cosets': len(left_cosets)
    }


def sylow_theorems(G: Group, p: int) -> dict:
    """
    @brief Wendet die Sylow-Sätze für die Primzahl p auf G an.
    @description
        Sei |G| = p^k · m mit gcd(p, m) = 1. Dann gilt:
        1. Sylow-Satz: G besitzt eine Sylow-p-Untergruppe (Ordnung p^k)
        2. Sylow-Satz: Alle Sylow-p-Untergruppen sind konjugiert
        3. Sylow-Satz: n_p ≡ 1 (mod p) und n_p | m

    @param G Endliche Gruppe
    @param p Primzahl
    @return Dictionary mit:
            - 'p': Primzahl
            - 'p_power': k (maximale p-Potenz in |G|)
            - 'index_m': m = |G|/p^k
            - 'sylow_subgroups': Liste der Sylow-p-Untergruppen
            - 'n_p': Anzahl der Sylow-p-Untergruppen
            - 'checks': {'congruence': n_p ≡ 1 (mod p), 'divisibility': n_p | m}
    @lastModified 2026-03-10
    """
    order_G = G.order()

    # Maximale p-Potenz in |G| bestimmen
    k = 0
    temp = order_G
    while temp % p == 0:
        k += 1
        temp //= p
    m = order_G // (p ** k)  # m = |G| / p^k, gcd(p, m) = 1

    # Sylow-p-Untergruppen der Ordnung p^k finden
    sylow_order = p ** k
    all_subs = G.subgroups()
    sylow_subs = [H for H in all_subs if H.order() == sylow_order]

    n_p = len(sylow_subs)

    # Sylow-Bedingungen überprüfen
    checks = {
        'congruence_mod_p': (n_p % p == 1),  # n_p ≡ 1 (mod p)
        'divides_m': (m % n_p == 0) if n_p > 0 else True,  # n_p | m
        'exists': n_p > 0  # Existenz der Sylow-UG
    }

    # Konjugiertheit der Sylow-UG verifizieren (für n_p > 1)
    if n_p > 1:
        # Prüfen ob H1 und H2 konjugiert sind: ∃ g: gH1g⁻¹ = H2
        checks['conjugate'] = _are_conjugate(G, sylow_subs[0], sylow_subs[1])
    else:
        checks['conjugate'] = True  # Trivialerweise konjugiert

    return {
        'p': p,
        'p_power': k,
        'sylow_order': sylow_order,
        'index_m': m,
        'sylow_subgroups': sylow_subs,
        'n_p': n_p,
        'checks': checks
    }


def _are_conjugate(G: Group, H1: Subgroup, H2: Subgroup) -> bool:
    """
    @brief Prüft ob H1 und H2 konjugiert in G sind (∃ g ∈ G: gH1g⁻¹ = H2).
    @param G Übergruppe
    @param H1 Erste Untergruppe
    @param H2 Zweite Untergruppe
    @return True wenn konjugiert
    @lastModified 2026-03-10
    """
    h2_set = frozenset(H2.elements)
    for g in G.elements:
        g_inv = G.inverse(g)
        # Konjugierte Untergruppe: gH1g⁻¹
        conj = frozenset(
            G.operation(G.operation(g, h), g_inv)
            for h in H1.elements
        )
        if conj == h2_set:
            return True
    return False


def classify_abelian_group(n: int) -> dict:
    """
    @brief Klassifiziert alle abelschen Gruppen der Ordnung n (bis auf Isomorphie).
    @description
        Hauptsatz über endliche abelsche Gruppen:
        Jede endliche abelsche Gruppe G der Ordnung n ist isomorph zu einem
        direkten Produkt zyklischer Gruppen:
        G ≅ ℤ_{p1^a1} × ℤ_{p2^a2} × ... × ℤ_{pk^ak}
        wobei n = p1^a1 · p2^a2 · ... · pk^ak (Primfaktorzerlegung).

        Die Anzahl der Isomorphieklassen entspricht dem Produkt der
        Partitionen der Exponenten in der Primfaktorzerlegung.

    @param n Gruppenordnung
    @return Dictionary mit:
            - 'n': Ordnung
            - 'factorization': Primfaktorzerlegung von n
            - 'isomorphism_classes': Liste der Isomorphieklassen (als Tupel-Listen)
            - 'count': Anzahl der Isomorphieklassen
            - 'is_cyclic': True wenn ℤ_n die einzige Klasse ist
    @lastModified 2026-03-10
    """
    factors = _factorint(n)

    # Partitionen eines Exponenten e berechnen (p^e → Partitionen von e)
    def partitions(e: int) -> list[list[int]]:
        """Alle Partitionen von e als Liste."""
        if e == 0:
            return [[]]
        result = []
        def gen_parts(remaining: int, max_part: int, current: list[int]) -> None:
            if remaining == 0:
                result.append(current[:])
                return
            for part in range(min(remaining, max_part), 0, -1):
                current.append(part)
                gen_parts(remaining - part, part, current)
                current.pop()
        gen_parts(e, e, [])
        return result

    # Für jeden Primfaktor: alle Partitionen des Exponenten
    prime_partitions: list[tuple[int, list[list[int]]]] = []
    for p, e in sorted(factors.items()):
        parts = partitions(e)
        prime_partitions.append((p, parts))

    # Kreuzprodukt aller Partitionen
    all_classes = []
    for combo in itertools.product(*[parts for _, parts in prime_partitions]):
        # Kombination: für jeden Primfaktor eine Partition
        class_groups = []
        for i, (p, _) in enumerate(prime_partitions):
            for exp in combo[i]:
                class_groups.append(f"Z_{p**exp}")
        all_classes.append(class_groups)

    # Zyklisch wenn alle Partitionen [e] (maximale Potenz) sind
    is_cyclic_check = len(all_classes) == 1

    return {
        'n': n,
        'factorization': factors,
        'isomorphism_classes': all_classes,
        'count': len(all_classes),
        'is_cyclic': is_cyclic_check,
        'description': [' × '.join(cls) if cls else 'trivial' for cls in all_classes]
    }


def is_simple_group(G: Group) -> bool:
    """
    @brief Prüft ob G eine einfache Gruppe ist.
    @description
        G ist einfach ⟺ die einzigen Normalteiler sind {e} und G selbst.

        Bekannte einfache Gruppen:
        - ℤ_p für Primzahlen p (abelsche einfache Gruppen)
        - A_n für n ≥ 5 (alternierende Gruppen)
        - Die 26 sporadischen einfachen Gruppen (z.B. Monster-Gruppe)

        Optimierung für große Gruppen: Konjugationsklassen-basierter Test.
        Ein Normalteiler ist stets eine Vereinigung von Konjugationsklassen
        (inklusive {e} und G). Wenn keine nicht-triviale solche Vereinigung
        eine Untergruppe bildet, ist G einfach.

    @param G Endliche Gruppe
    @return True wenn G einfach ist
    @lastModified 2026-03-10
    """
    order_G = G.order()

    # Triviale Fälle
    if order_G == 1:
        return False  # Triviale Gruppe gilt per Konvention nicht als einfach

    # Für kleine Gruppen: direkter Normalteiler-Test
    if order_G <= 30:
        normal_subs = G.normal_subgroups()
        trivial_orders = {1, order_G}
        return all(N.order() in trivial_orders for N in normal_subs)

    # Für größere Gruppen: Konjugationsklassen-Methode
    # Normalteiler = Vereinigung von Konjugationsklassen (enthält immer e und G)
    # Schritt 1: Konjugationsklassen berechnen
    conjugacy_classes = _conjugacy_classes(G)

    # Schritt 2: Prüfen ob eine echte nicht-triviale Vereinigung von Klassen
    # eine Untergruppe (also Normalteiler) bildet
    # Die Klassen von e und G müssen immer enthalten sein
    e_class = next(cls for cls in conjugacy_classes if G.identity in cls)

    # Alle Teilmengen der Konjugationsklassen (ohne leere Menge und ganz G)
    classes_without_e = [cls for cls in conjugacy_classes if G.identity not in cls]

    for r in range(1, len(classes_without_e) + 1):
        for combo in itertools.combinations(classes_without_e, r):
            # Kandidat: {e} ∪ Vereinigung der gewählten Klassen
            candidate = list(e_class)
            for cls in combo:
                candidate.extend(cls)
            # Ordnung muss |G| teilen
            if order_G % len(candidate) != 0:
                continue
            # Prüfen ob Kandidat Untergruppe (= Normalteiler) bildet
            if G._is_subgroup(candidate) and 1 < len(candidate) < order_G:
                return False

    return True


def _conjugacy_classes(G: Group) -> list[list[Any]]:
    """
    @brief Berechnet alle Konjugationsklassen von G.
    @description
        Die Konjugationsklasse von g ist: Cl(g) = {x·g·x⁻¹ | x ∈ G}
        Konjugationsklassen bilden eine Partition von G.
    @param G Endliche Gruppe
    @return Liste der Konjugationsklassen
    @lastModified 2026-03-10
    """
    remaining = set(range(len(G.elements)))
    elem_to_idx = {e: i for i, e in enumerate(G.elements)}
    classes = []

    while remaining:
        # Repräsentant wählen
        idx = next(iter(remaining))
        g = G.elements[idx]

        # Konjugationsklasse von g berechnen: {x·g·x⁻¹ | x ∈ G}
        conj_class = set()
        for x in G.elements:
            x_inv = G.inverse(x)
            conj = G.operation(G.operation(x, g), x_inv)
            conj_class.add(elem_to_idx[conj])

        cls_elements = [G.elements[i] for i in conj_class]
        classes.append(cls_elements)
        remaining -= conj_class

    return classes


def group_action(G: Group, X: list[Any], action: Callable[[Any, Any], Any]) -> dict:
    """
    @brief Analysiert eine Gruppenoperation G × X → X.
    @description
        Eine Gruppenoperation G × X → X erfüllt:
        1. Identität: e · x = x für alle x ∈ X
        2. Verträglichkeit: (g·h) · x = g · (h · x) für alle g,h ∈ G, x ∈ X

        Wichtige Konzepte:
        - Orbit: Orb(x) = {g·x | g ∈ G}
        - Stabilisator: Stab(x) = {g ∈ G | g·x = x}
        - Bahnformel: |Orb(x)| · |Stab(x)| = |G|
        - Burnside-Lemma: |X/G| = (1/|G|) · Σ_{g∈G} |Fix(g)|
          wobei Fix(g) = {x ∈ X | g·x = x}

    @param G Gruppe
    @param X Menge auf der G operiert
    @param action Funktion (g, x) → g·x
    @return Dictionary mit:
            - 'orbits': Liste der Orbits
            - 'stabilizers': Dict x → Stab(x)
            - 'fixed_points': Dict g → Fix(g)
            - 'num_orbits_burnside': Anzahl Orbits via Burnside
            - 'is_transitive': True wenn |X/G| = 1
    @lastModified 2026-03-10
    """
    # Orbits berechnen (Union-Find-ähnlich)
    covered = set()
    orbits = []
    for x in X:
        x_key = x if not isinstance(x, list) else tuple(x)
        if x_key not in covered:
            orbit = []
            for g in G.elements:
                gx = action(g, x)
                gx_key = gx if not isinstance(gx, list) else tuple(gx)
                if gx_key not in {(o if not isinstance(o, list) else tuple(o)) for o in orbit}:
                    orbit.append(gx)
                covered.add(gx_key)
            orbits.append(orbit)

    # Stabilisatoren berechnen
    stabilizers = {}
    for x in X:
        stab = [g for g in G.elements if action(g, x) == x]
        stabilizers[str(x)] = stab

    # Fixpunkte für Burnside-Lemma: Fix(g) = {x ∈ X | g·x = x}
    fixed_points = {}
    burnside_sum = 0
    for g in G.elements:
        fix_g = [x for x in X if action(g, x) == x]
        fixed_points[str(g)] = fix_g
        burnside_sum += len(fix_g)

    # Burnside-Lemma: |X/G| = (1/|G|) · Σ |Fix(g)|
    num_orbits_burnside = burnside_sum // G.order()

    return {
        'orbits': orbits,
        'num_orbits': len(orbits),
        'stabilizers': stabilizers,
        'fixed_points': fixed_points,
        'burnside_sum': burnside_sum,
        'num_orbits_burnside': num_orbits_burnside,
        'is_transitive': len(orbits) == 1
    }


def direct_product(G: Group, H: Group) -> Group:
    """
    @brief Bildet das direkte Produkt G × H.
    @description
        G × H = {(g, h) | g ∈ G, h ∈ H}
        mit komponentenweiser Operation:
        (g1, h1) · (g2, h2) = (g1·g2, h1·h2)

        Ordnung: |G × H| = |G| · |H|
        G × H ist abelsch ⟺ G und H sind abelsch.

    @param G Erste Gruppe
    @param H Zweite Gruppe
    @return Direktes Produkt als Group-Objekt
    @lastModified 2026-03-10
    """
    # Elemente des Produkts: alle Paare (g, h)
    elements = [(g, h) for g in G.elements for h in H.elements]

    # Komponentenweise Operation
    def product_op(a: tuple, b: tuple) -> tuple:
        return (G.operation(a[0], b[0]), H.operation(a[1], b[1]))

    # Neutrales Element: (e_G, e_H)
    identity = (G.identity, H.identity)

    return Group(elements, product_op, identity, name=f"{G.name}×{H.name}")


def semidirect_product(N: Group, H: Group, phi: dict) -> Group:
    """
    @brief Bildet das halbdirekte Produkt N ⋊_φ H.
    @description
        Das halbdirekte Produkt N ⋊_φ H hat:
        - Elementen: {(n, h) | n ∈ N, h ∈ H}
        - Operation: (n1, h1) · (n2, h2) = (n1 · φ(h1)(n2), h1 · h2)

        Dabei ist φ: H → Aut(N) ein Homomorphismus.
        Für φ trivial: N ⋊ H = N × H (direktes Produkt).

    @param N Normalteiler-Untergruppe
    @param H Komplement-Untergruppe
    @param phi Dictionary: h → (Automorphismus von N als Dict n → φ(h)(n))
    @return Halbdirektes Produkt als Group-Objekt
    @lastModified 2026-03-10
    """
    # Elemente: alle Paare (n, h)
    elements = [(n, h) for n in N.elements for h in H.elements]

    # Halbdirektes Produkt Operation: (n1, h1)·(n2, h2) = (n1·φ(h1)(n2), h1·h2)
    def semi_op(a: tuple, b: tuple) -> tuple:
        n1, h1 = a
        n2, h2 = b
        # φ(h1) auf n2 anwenden
        phi_h1 = phi[h1]  # Automorphismus als Dictionary n → φ(h1)(n)
        phi_n2 = phi_h1[n2]
        # n1 · φ(h1)(n2)
        new_n = N.operation(n1, phi_n2)
        # h1 · h2
        new_h = H.operation(h1, h2)
        return (new_n, new_h)

    identity = (N.identity, H.identity)

    return Group(elements, semi_op, identity, name=f"{N.name}⋊{H.name}")


def symmetric_group(n: int) -> PermutationGroup:
    """
    @brief Erzeugt die symmetrische Gruppe S_n aller Permutationen von {0,...,n-1}.
    @description
        S_n = Gruppe aller Bijektionen {0,...,n-1} → {0,...,n-1}
        Ordnung: |S_n| = n!
        S_n ist nicht-abelsch für n ≥ 3.

    @param n Grad der Permutationsgruppe
    @return S_n als PermutationGroup
    @lastModified 2026-03-10
    """
    # Alle Permutationen von {0, ..., n-1} generieren
    all_perms = list(itertools.permutations(range(n)))
    return PermutationGroup(all_perms, n, name=f"S_{n}")


def alternating_group(n: int) -> PermutationGroup:
    """
    @brief Erzeugt die alternierende Gruppe A_n (gerade Permutationen).
    @description
        A_n = {σ ∈ S_n | sgn(σ) = +1}
        Ordnung: |A_n| = n!/2 für n ≥ 2
        A_n ist einfach für n ≥ 5 (klassisches Resultat der Gruppentheorie).
        A_3 ≅ ℤ_3 (zyklische Gruppe der Ordnung 3).

    @param n Grad
    @return A_n als PermutationGroup
    @lastModified 2026-03-10
    """
    # Alle geraden Permutationen aus S_n auswählen
    Sn = symmetric_group(n)
    even_perms = [sigma for sigma in Sn.elements if Sn.is_even_permutation(sigma)]
    return PermutationGroup(even_perms, n, name=f"A_{n}")


def cyclic_group(n: int) -> Group:
    """
    @brief Erzeugt die zyklische Gruppe ℤ_n = ({0,1,...,n-1}, + mod n).
    @description
        ℤ_n ist die einzige abelsche Gruppe der Ordnung n (bis auf Isomorphie)
        wenn n prim ist. Für zusammengesetzte n gibt es mehrere abelsche Gruppen.

        Eigenschaften:
        - Abelsch: ∀ a,b: a+b ≡ b+a (mod n)
        - Zyklisch: erzeugt von 1
        - Generator: alle k mit gcd(k, n) = 1 sind Generatoren
        - |generators| = φ(n) (Eulersche Phi-Funktion)

    @param n Gruppenordnung
    @return ℤ_n als Group-Objekt
    @lastModified 2026-03-10
    """
    elements = list(range(n))
    # Addition modulo n als Gruppenoperation
    operation = lambda a, b: (a + b) % n
    return Group(elements, operation, identity=0, name=f"Z_{n}")


def dihedral_group(n: int) -> Group:
    """
    @brief Erzeugt die Diedergruppe D_n (Symmetriegruppe des regulären n-Ecks).
    @description
        D_n hat Ordnung 2n und besteht aus:
        - n Drehungen: r^k für k = 0,...,n-1 (dargestellt als (k, 0))
        - n Spiegelungen: s·r^k für k = 0,...,n-1 (dargestellt als (k, 1))

        Gruppenstruktur:
        - r^n = e (Drehung um 360°)
        - s^2 = e (Spiegelung ist Involution)
        - s·r = r^{n-1}·s (Kommutationsrelation)

        Operation (k1, f1) · (k2, f2):
        - Wenn f1 = 0 (Drehung): (k1+k2 mod n, f2)
        - Wenn f1 = 1 (Spiegelung): (k1-k2 mod n, 1-f2)... korrektes Dihedral

    @param n Ordnung des Polygons (D_n hat Ordnung 2n)
    @return D_n als Group-Objekt
    @lastModified 2026-03-10
    """
    # Elemente als Paare (k, f): k ∈ {0,...,n-1}, f ∈ {0,1}
    # (k, 0) = r^k (Drehung), (k, 1) = s·r^k (Spiegelung)
    elements = [(k, f) for f in range(2) for k in range(n)]

    def dihedral_op(a: tuple, b: tuple) -> tuple:
        """
        Diedergruppen-Multiplikation.
        r = Drehung um 2π/n, s = Spiegelung
        Relation: s·r = r^{-1}·s
        """
        k1, f1 = a
        k2, f2 = b
        if f1 == 0:
            # Drehung · Drehung/Spiegelung: r^k1 · r^k2·s^f2 = r^{k1+k2}·s^f2
            return ((k1 + k2) % n, f2)
        else:
            # Spiegelung · Drehung/Spiegelung: s·r^k1 · r^k2·s^f2
            # = s · r^{k1+k2} · s^f2
            # Relation: s·r^m = r^{-m}·s → s·r^m·s = r^{-m}
            # (k1, 1)·(k2, 0) = (k1-k2 mod n, 1)
            # (k1, 1)·(k2, 1) = (k1-k2 mod n, 0)
            return ((k1 - k2) % n, 1 - f2)

    identity = (0, 0)  # e = r^0

    return Group(elements, dihedral_op, identity, name=f"D_{n}")


def quaternion_group() -> Group:
    """
    @brief Erzeugt die Quaternionengruppe Q_8 = {±1, ±i, ±j, ±k}.
    @description
        Q_8 ist eine nicht-abelsche Gruppe der Ordnung 8 mit den Relationen:
        i² = j² = k² = -1
        ij = k,  ji = -k
        jk = i,  kj = -i
        ki = j,  ik = -j
        (-1)² = 1

        Q_8 ist die einzige nicht-abelsche Gruppe der Ordnung 8,
        die keine Untergruppe isomorph zu ℤ_2 × ℤ_2 hat.

        Darstellung als String-Elemente: {'1', '-1', 'i', '-i', 'j', '-j', 'k', '-k'}

    @return Q_8 als Group-Objekt
    @lastModified 2026-03-10
    """
    # Elemente als Strings
    elements = ['1', '-1', 'i', '-i', 'j', '-j', 'k', '-k']

    # Multiplikationstabelle für Q_8
    q8_table = {
        ('1', '1'): '1',    ('1', '-1'): '-1',  ('1', 'i'): 'i',    ('1', '-i'): '-i',
        ('1', 'j'): 'j',    ('1', '-j'): '-j',  ('1', 'k'): 'k',    ('1', '-k'): '-k',
        ('-1', '1'): '-1',  ('-1', '-1'): '1',  ('-1', 'i'): '-i',  ('-1', '-i'): 'i',
        ('-1', 'j'): '-j',  ('-1', '-j'): 'j',  ('-1', 'k'): '-k',  ('-1', '-k'): 'k',
        ('i', '1'): 'i',    ('i', '-1'): '-i',  ('i', 'i'): '-1',   ('i', '-i'): '1',
        ('i', 'j'): 'k',    ('i', '-j'): '-k',  ('i', 'k'): '-j',   ('i', '-k'): 'j',
        ('-i', '1'): '-i',  ('-i', '-1'): 'i',  ('-i', 'i'): '1',   ('-i', '-i'): '-1',
        ('-i', 'j'): '-k',  ('-i', '-j'): 'k',  ('-i', 'k'): 'j',   ('-i', '-k'): '-j',
        ('j', '1'): 'j',    ('j', '-1'): '-j',  ('j', 'i'): '-k',   ('j', '-i'): 'k',
        ('j', 'j'): '-1',   ('j', '-j'): '1',   ('j', 'k'): 'i',    ('j', '-k'): '-i',
        ('-j', '1'): '-j',  ('-j', '-1'): 'j',  ('-j', 'i'): 'k',   ('-j', '-i'): '-k',
        ('-j', 'j'): '1',   ('-j', '-j'): '-1', ('-j', 'k'): '-i',  ('-j', '-k'): 'i',
        ('k', '1'): 'k',    ('k', '-1'): '-k',  ('k', 'i'): 'j',    ('k', '-i'): '-j',
        ('k', 'j'): '-i',   ('k', '-j'): 'i',   ('k', 'k'): '-1',   ('k', '-k'): '1',
        ('-k', '1'): '-k',  ('-k', '-1'): 'k',  ('-k', 'i'): '-j',  ('-k', '-i'): 'j',
        ('-k', 'j'): 'i',   ('-k', '-j'): '-i', ('-k', 'k'): '1',   ('-k', '-k'): '-1',
    }

    def q8_op(a: str, b: str) -> str:
        return q8_table[(a, b)]

    return Group(elements, q8_op, identity='1', name="Q_8")
