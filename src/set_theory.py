"""
@file set_theory.py
@brief Vollständiges Mengenlehre-Modul: Mengen, Relationen, Funktionen,
       Kardinalitätstheorie, Ordinalzahlen, ZFC-Axiome, Mengenfamilien, Filter.

@description
    Dieses Modul implementiert die grundlegenden Konzepte der mathematischen Mengenlehre
    objektorientiert und domänengetrieben. Es deckt folgende Bereiche ab:

    1. MathSet          – Mathematische Menge (union, intersection, power_set, …)
    2. Relation         – Binäre Relation R ⊆ A×B mit Abschlussoperationen
    3. MathFunction     – Funktion f: A→B als spezielle Relation
    4. Kardinalitätstheorie – Cantor-Schröder-Bernstein, Diagonalargument
    5. Ordinal          – Ordinalzahlen (endlich + transfinit: ω, ω+1, …)
    6. ZFC / AC         – Axiomensystem-Beschreibungen
    7. SetFamily        – Familie von Mengen, σ-Algebren
    8. Filter           – Mengentheoretischer Filter / Ultrafilter

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

from __future__ import annotations
import itertools
from typing import Any, Callable, Optional


# ═══════════════════════════════════════════════════════════════════════════════
# 1. MathSet – Mathematische Menge
# ═══════════════════════════════════════════════════════════════════════════════

class MathSet:
    """
    Mathematische Menge mit allen Standardoperationen der Mengenlehre.

    Intern als frozenset gespeichert, damit MathSet-Objekte selbst als
    Elemente anderer Mengen (z.B. Potenzmenge) verwendbar sind (hashbar).

    Beispiel:
        A = MathSet([1, 2, 3])
        B = MathSet([2, 3, 4])
        print(A | B)   # {1, 2, 3, 4}

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, elements=None):
        """
        Initialisiert die Menge.

        @param elements Iterierbares Objekt mit den Elementen oder None für ∅.
        @type  elements iterable | None
        """
        if elements is None:
            # Leere Menge ∅
            self._data: frozenset = frozenset()
        else:
            self._data = frozenset(elements)

    # ── Grundoperationen ────────────────────────────────────────────────────

    def union(self, other: 'MathSet') -> 'MathSet':
        """
        Vereinigung zweier Mengen: A ∪ B = {x : x∈A oder x∈B}.

        @param other Die zweite Menge B.
        @return      Neue MathSet-Instanz mit A ∪ B.
        """
        return MathSet(self._data | other._data)

    def intersection(self, other: 'MathSet') -> 'MathSet':
        """
        Durchschnitt zweier Mengen: A ∩ B = {x : x∈A und x∈B}.

        @param other Die zweite Menge B.
        @return      Neue MathSet-Instanz mit A ∩ B.
        """
        return MathSet(self._data & other._data)

    def difference(self, other: 'MathSet') -> 'MathSet':
        """
        Mengendifferenz: A \\ B = {x : x∈A und x∉B}.

        @param other Die zweite Menge B.
        @return      Neue MathSet-Instanz mit A \\ B.
        """
        return MathSet(self._data - other._data)

    def symmetric_difference(self, other: 'MathSet') -> 'MathSet':
        """
        Symmetrische Differenz: A △ B = (A \\ B) ∪ (B \\ A).

        @param other Die zweite Menge B.
        @return      Neue MathSet-Instanz mit A △ B.
        """
        return MathSet(self._data ^ other._data)

    def complement(self, universe: 'MathSet') -> 'MathSet':
        """
        Komplement bezüglich eines Universums U: A^c = U \\ A.

        @param universe Das Universum U, das A enthält.
        @return         Neue MathSet-Instanz mit U \\ A.
        @raises ValueError wenn A keine Teilmenge von U ist.
        """
        if not self._data <= universe._data:
            raise ValueError("A muss eine Teilmenge des Universums U sein.")
        return MathSet(universe._data - self._data)

    def power_set(self) -> 'MathSet':
        """
        Potenzmenge: P(A) = {S : S ⊆ A}.

        Für eine n-elementige Menge enthält P(A) genau 2^n Elemente.

        @return MathSet, dessen Elemente selbst frozenset-Instanzen sind.
        @raises ValueError wenn |A| > 20 (Schutz vor Speicherüberlauf).
        """
        if len(self._data) > 20:
            raise ValueError("Potenzmenge nur für |A| ≤ 20 berechenbar (2^20 = 1 048 576 Elemente).")
        # Erzeuge alle Teilmengen mittels itertools.chain.from_iterable
        elements = list(self._data)
        subsets = []
        for r in range(len(elements) + 1):
            for combo in itertools.combinations(elements, r):
                subsets.append(frozenset(combo))
        # Die Elemente der Potenzmenge sind selbst frozensets (hashbar)
        return MathSet(subsets)

    def cartesian_product(self, other: 'MathSet') -> 'MathSet':
        """
        Kartesisches Produkt: A × B = {(a,b) : a∈A, b∈B}.

        @param other Die zweite Menge B.
        @return      MathSet von Tupeln (a, b).
        """
        pairs = set(itertools.product(self._data, other._data))
        return MathSet(pairs)

    # ── Mengeneigenschaften ─────────────────────────────────────────────────

    def is_subset(self, other: 'MathSet') -> bool:
        """
        Prüft Teilmengenbeziehung: A ⊆ B ⟺ ∀x∈A: x∈B.

        @param other Menge B.
        @return      True wenn A ⊆ B.
        """
        return self._data <= other._data

    def is_proper_subset(self, other: 'MathSet') -> bool:
        """
        Prüft echte Teilmengenbeziehung: A ⊊ B ⟺ A ⊆ B und A ≠ B.

        @param other Menge B.
        @return      True wenn A ⊊ B.
        """
        return self._data < other._data

    def is_disjoint(self, other: 'MathSet') -> bool:
        """
        Prüft ob A und B disjunkt sind: A ∩ B = ∅.

        @param other Menge B.
        @return      True wenn A ∩ B = ∅.
        """
        return self._data.isdisjoint(other._data)

    def cardinality(self) -> int:
        """
        Mächtigkeit (Kardinalzahl) der Menge: |A|.

        @return Anzahl der Elemente in A.
        """
        return len(self._data)

    # ── Magische Methoden ───────────────────────────────────────────────────

    def __or__(self, other: 'MathSet') -> 'MathSet':
        """A | B  →  Vereinigung A ∪ B."""
        return self.union(other)

    def __and__(self, other: 'MathSet') -> 'MathSet':
        """A & B  →  Durchschnitt A ∩ B."""
        return self.intersection(other)

    def __sub__(self, other: 'MathSet') -> 'MathSet':
        """A - B  →  Differenz A \\ B."""
        return self.difference(other)

    def __xor__(self, other: 'MathSet') -> 'MathSet':
        """A ^ B  →  Symmetrische Differenz A △ B."""
        return self.symmetric_difference(other)

    def __le__(self, other: 'MathSet') -> bool:
        """A <= B  →  Teilmenge A ⊆ B."""
        return self.is_subset(other)

    def __lt__(self, other: 'MathSet') -> bool:
        """A < B  →  Echte Teilmenge A ⊊ B."""
        return self.is_proper_subset(other)

    def __eq__(self, other: object) -> bool:
        """A == B  →  Mengengleichheit (Extensionalitätsaxiom)."""
        if isinstance(other, MathSet):
            return self._data == other._data
        return NotImplemented

    def __contains__(self, x: Any) -> bool:
        """x in A  →  Elementbeziehung x ∈ A."""
        return x in self._data

    def __len__(self) -> int:
        """len(A)  →  Kardinalzahl |A|."""
        return len(self._data)

    def __iter__(self):
        """Iteration über die Elemente von A."""
        return iter(self._data)

    def __hash__(self) -> int:
        """Hashwert (benötigt für Verwendung als dict-Key oder in Mengen)."""
        return hash(self._data)

    def __repr__(self) -> str:
        """Lesbare Darstellung der Menge in Mengennotation."""
        if not self._data:
            return "∅"
        # Sortiere wenn möglich (für reproduzierbare Ausgabe)
        try:
            sorted_elems = sorted(self._data)
        except TypeError:
            sorted_elems = list(self._data)
        return "{" + ", ".join(repr(e) for e in sorted_elems) + "}"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Relation – Binäre Relation
# ═══════════════════════════════════════════════════════════════════════════════

class Relation:
    """
    Binäre Relation R ⊆ A × B.

    Eine binäre Relation ist eine Menge von geordneten Paaren (a, b),
    wobei a aus der Definitionsmenge (Domain) und b aus der Zielmenge
    (Codomain) stammt.

    Beispiel:
        A = MathSet([1, 2, 3])
        R = Relation({(1,1), (2,2), (3,3)}, A, A)
        print(R.is_reflexive())   # True

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, pairs: set, domain: MathSet, codomain: MathSet):
        """
        Initialisiert die Relation.

        @param pairs    Menge von Tupeln (a, b) mit a∈domain, b∈codomain.
        @param domain   Definitionsmenge A.
        @param codomain Zielmenge B.
        @raises ValueError wenn ein Paar außerhalb von A×B liegt.
        """
        self.domain = domain
        self.codomain = codomain
        # Validierung: alle Paare müssen in A×B liegen
        for a, b in pairs:
            if a not in domain:
                raise ValueError(f"Element {a!r} nicht in der Definitionsmenge.")
            if b not in codomain:
                raise ValueError(f"Element {b!r} nicht in der Zielmenge.")
        self.pairs: frozenset = frozenset(pairs)

    # ── Relationseigenschaften ──────────────────────────────────────────────

    def is_reflexive(self) -> bool:
        """
        Reflexivität: ∀a∈A: (a,a)∈R.

        Nur definiert wenn Domain = Codomain (homogene Relation).

        @return True wenn R reflexiv ist.
        """
        # Homogene Relation vorausgesetzt
        return all((a, a) in self.pairs for a in self.domain)

    def is_symmetric(self) -> bool:
        """
        Symmetrie: ∀(a,b)∈R: (b,a)∈R.

        @return True wenn R symmetrisch ist.
        """
        return all((b, a) in self.pairs for (a, b) in self.pairs)

    def is_antisymmetric(self) -> bool:
        """
        Antisymmetrie: (a,b)∈R ∧ (b,a)∈R ⟹ a=b.

        @return True wenn R antisymmetrisch ist.
        """
        for (a, b) in self.pairs:
            if a != b and (b, a) in self.pairs:
                return False
        return True

    def is_transitive(self) -> bool:
        """
        Transitivität: (a,b)∈R ∧ (b,c)∈R ⟹ (a,c)∈R.

        @return True wenn R transitiv ist.
        """
        for (a, b) in self.pairs:
            for (b2, c) in self.pairs:
                if b == b2 and (a, c) not in self.pairs:
                    return False
        return True

    def is_equivalence(self) -> bool:
        """
        Äquivalenzrelation: reflexiv ∧ symmetrisch ∧ transitiv.

        @return True wenn R eine Äquivalenzrelation ist.
        """
        return self.is_reflexive() and self.is_symmetric() and self.is_transitive()

    def is_partial_order(self) -> bool:
        """
        Partielle Ordnung: reflexiv ∧ antisymmetrisch ∧ transitiv.

        @return True wenn R eine partielle Ordnung ist.
        """
        return self.is_reflexive() and self.is_antisymmetric() and self.is_transitive()

    def is_total_order(self) -> bool:
        """
        Totalordnung: partielle Ordnung + alle Elemente sind vergleichbar.

        ∀a,b∈A: (a,b)∈R oder (b,a)∈R.

        @return True wenn R eine Totalordnung ist.
        """
        if not self.is_partial_order():
            return False
        # Totalitätsbedingung: alle Paare vergleichbar
        for a in self.domain:
            for b in self.domain:
                if (a, b) not in self.pairs and (b, a) not in self.pairs:
                    return False
        return True

    # ── Äquivalenzklassen ───────────────────────────────────────────────────

    def equivalence_classes(self) -> list[MathSet]:
        """
        Berechnet die Äquivalenzklassen [a]_R = {b∈A : (a,b)∈R}.

        Setzt voraus, dass R eine Äquivalenzrelation ist.

        @return Liste der Äquivalenzklassen (jede als MathSet).
        @raises ValueError wenn R keine Äquivalenzrelation ist.
        """
        if not self.is_equivalence():
            raise ValueError("R ist keine Äquivalenzrelation.")

        # Union-Find-ähnlicher Ansatz: gruppiere Elemente in Klassen
        classes: dict[Any, set] = {}
        for a in self.domain:
            # Suche ob a bereits in einer Klasse liegt
            representative = None
            for rep in classes:
                if (a, rep) in self.pairs:
                    representative = rep
                    break
            if representative is None:
                classes[a] = {a}
            else:
                classes[representative].add(a)

        return [MathSet(cls) for cls in classes.values()]

    def quotient_set(self) -> 'MathSet':
        """
        Quotientenmenge A/R = Menge aller Äquivalenzklassen.

        @return MathSet, dessen Elemente frozensets (Äquivalenzklassen) sind.
        """
        classes = self.equivalence_classes()
        # Jede Äquivalenzklasse wird als frozenset repräsentiert
        return MathSet(frozenset(cls._data) for cls in classes)

    # ── Abschlussoperationen ────────────────────────────────────────────────

    def transitive_closure(self) -> 'Relation':
        """
        Transitiver Abschluss R⁺ via Floyd-Warshall-Algorithmus.

        Floyd-Warshall: R_new = R ∪ {(a,c) : ∃b: (a,b)∈R ∧ (b,c)∈R}
        Iteriert bis Fixpunkt erreicht.

        @return Neue Relation mit transitivem Abschluss.
        """
        # Starte mit aktuellen Paaren
        closure = set(self.pairs)
        # Iteriere bis kein neues Paar mehr hinzukommt (Fixpunkt)
        changed = True
        while changed:
            changed = False
            new_pairs = set()
            for (a, b) in closure:
                for (b2, c) in closure:
                    if b == b2 and (a, c) not in closure:
                        new_pairs.add((a, c))
                        changed = True
            closure |= new_pairs
        return Relation(closure, self.domain, self.codomain)

    def reflexive_closure(self) -> 'Relation':
        """
        Reflexiver Abschluss: R ∪ {(a,a) : a∈A}.

        @return Neue Relation mit reflexivem Abschluss.
        """
        closure = set(self.pairs)
        for a in self.domain:
            closure.add((a, a))
        return Relation(closure, self.domain, self.codomain)

    def symmetric_closure(self) -> 'Relation':
        """
        Symmetrischer Abschluss: R ∪ R⁻¹ = R ∪ {(b,a) : (a,b)∈R}.

        @return Neue Relation mit symmetrischem Abschluss.
        """
        closure = set(self.pairs)
        for (a, b) in self.pairs:
            closure.add((b, a))
        return Relation(closure, self.domain, self.codomain)

    def composition(self, other: 'Relation') -> 'Relation':
        """
        Komposition S∘R: (a,c) ∈ S∘R ⟺ ∃b: (a,b)∈R ∧ (b,c)∈S.

        Notation: self = R (Definitionsbereich), other = S.

        @param other Relation S (domain muss zu codomain von R passen).
        @return      Komposition S∘R.
        """
        result = set()
        for (a, b) in self.pairs:
            for (b2, c) in other.pairs:
                if b == b2:
                    result.add((a, c))
        return Relation(result, self.domain, other.codomain)

    def inverse(self) -> 'Relation':
        """
        Inverse Relation R⁻¹ = {(b,a) : (a,b)∈R}.

        @return Neue Relation R⁻¹.
        """
        return Relation({(b, a) for (a, b) in self.pairs},
                        self.codomain, self.domain)

    def __repr__(self) -> str:
        """Lesbare Darstellung der Relation."""
        sorted_pairs = sorted(self.pairs) if self.pairs else []
        return f"Relation({set(sorted_pairs)})"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. MathFunction – Mathematische Funktion
# ═══════════════════════════════════════════════════════════════════════════════

class MathFunction:
    """
    Mathematische Funktion f: A → B als spezielle Relation.

    Eine Funktion ist eine linkstotale, rechtseindeutige Relation:
    - Linkstotal: ∀a∈A ∃b∈B: f(a) = b
    - Rechtseindeutig: f(a)=b ∧ f(a)=b' ⟹ b=b'

    Beispiel:
        A = MathSet([1, 2, 3])
        B = MathSet([2, 4, 6])
        f = MathFunction({1: 2, 2: 4, 3: 6}, A, B)
        print(f.is_bijective())  # True

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, mapping: dict, domain: MathSet, codomain: MathSet):
        """
        Initialisiert die Funktion.

        @param mapping  dict {a: f(a)} für alle a∈domain.
        @param domain   Definitionsmenge A.
        @param codomain Zielmenge B (nicht zwingend = Bild).
        @raises ValueError wenn mapping keine gültige Funktion f: A→B ist.
        """
        self.domain = domain
        self.codomain = codomain
        self.mapping: dict = dict(mapping)

        # Validierung: Linkstotalität (jedes a∈A hat ein Bild)
        for a in domain:
            if a not in self.mapping:
                raise ValueError(f"Element {a!r} aus A hat kein Bild (nicht linkstotal).")

        # Validierung: Bilder liegen in B
        for a, b in self.mapping.items():
            if a not in domain:
                raise ValueError(f"Schlüssel {a!r} nicht in der Definitionsmenge A.")
            if b not in codomain:
                raise ValueError(f"Bild {b!r} nicht in der Zielmenge B.")

    def __call__(self, a: Any) -> Any:
        """
        Funktionsauswertung f(a).

        @param a Element aus A.
        @return  f(a) ∈ B.
        """
        if a not in self.domain:
            raise ValueError(f"{a!r} nicht in der Definitionsmenge A.")
        return self.mapping[a]

    # ── Funktionseigenschaften ──────────────────────────────────────────────

    def is_injective(self) -> bool:
        """
        Injektivität (Eineindeutigkeit): f(a) = f(b) ⟹ a = b.

        @return True wenn f injektiv ist.
        """
        # Injektiv ⟺ alle Bilder verschieden
        images = list(self.mapping.values())
        return len(images) == len(set(images))

    def is_surjective(self) -> bool:
        """
        Surjektivität (Aufsurjektivität): ∀b∈B ∃a∈A: f(a) = b.

        @return True wenn f surjektiv ist.
        """
        return set(self.mapping.values()) == set(self.codomain)

    def is_bijective(self) -> bool:
        """
        Bijektivität: injektiv ∧ surjektiv.

        @return True wenn f bijektiv (umkehrbar eindeutig) ist.
        """
        return self.is_injective() and self.is_surjective()

    def inverse_function(self) -> 'MathFunction':
        """
        Umkehrfunktion f⁻¹: B → A (nur bei Bijektivität).

        @return MathFunction f⁻¹ mit f⁻¹(f(a)) = a.
        @raises ValueError wenn f nicht bijektiv ist.
        """
        if not self.is_bijective():
            raise ValueError("Umkehrfunktion existiert nur für bijektive Funktionen.")
        inverse_mapping = {b: a for a, b in self.mapping.items()}
        return MathFunction(inverse_mapping, self.codomain, self.domain)

    def composition(self, g: 'MathFunction') -> 'MathFunction':
        """
        Komposition g∘f: A → C, wobei f: A→B und g: B→C.

        (g∘f)(a) = g(f(a))

        @param g Funktion g: B→C (Codomain von f muss Domain von g sein).
        @return  Komposition g∘f.
        """
        composed = {a: g.mapping[self.mapping[a]] for a in self.domain}
        return MathFunction(composed, self.domain, g.codomain)

    def image(self, subset: Optional[MathSet] = None) -> MathSet:
        """
        Bild (Image) von f oder einer Teilmenge S ⊆ A.

        f(A) = {f(a) : a∈A} oder f(S) = {f(a) : a∈S}.

        @param subset Teilmenge S ⊆ A (optional, Standard: ganzes A).
        @return       Bild f(S) als MathSet.
        """
        if subset is None:
            subset = self.domain
        return MathSet(self.mapping[a] for a in subset if a in self.mapping)

    def preimage(self, subset: MathSet) -> MathSet:
        """
        Urbild (Preimage) einer Teilmenge T ⊆ B.

        f⁻¹(T) = {a∈A : f(a)∈T}.

        @param subset Teilmenge T ⊆ B.
        @return       Urbild f⁻¹(T) als MathSet.
        """
        return MathSet(a for a in self.domain if self.mapping[a] in subset)

    def is_identity(self) -> bool:
        """
        Prüft ob f die Identitätsfunktion ist: ∀a∈A: f(a) = a.

        @return True wenn f = id_A.
        """
        return all(self.mapping[a] == a for a in self.domain)

    def __repr__(self) -> str:
        """Lesbare Darstellung der Funktion."""
        return f"MathFunction({self.mapping})"


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Kardinalitätstheorie
# ═══════════════════════════════════════════════════════════════════════════════

def are_equinumerous(A: MathSet, B: MathSet) -> bool:
    """
    Prüft Gleichmächtigkeit |A| = |B| für endliche Mengen.

    Für endliche Mengen: Bijektion existiert ⟺ |A| = |B|.

    @param A Erste Menge.
    @param B Zweite Menge.
    @return  True wenn |A| = |B|.
    @lastModified 2026-03-10
    """
    return len(A) == len(B)


def cantor_schroeder_bernstein(
    f_AB: dict, f_BA: dict, A: MathSet, B: MathSet
) -> dict:
    """
    Konstruiert eine Bijektion aus zwei Injektionen (Cantor-Schröder-Bernstein-Satz).

    Satz (CSB): Wenn f: A→B injektiv und g: B→A injektiv, dann existiert
    eine Bijektion h: A→B.

    Konstruktion:
    - Berechne A₀ = A \\ g(B), dann A_{n+1} = g(f(A_n))
    - C = ⋃ A_n (Schroeder-Bernstein-Kette)
    - h(a) = f(a) falls a∈C, sonst h(a) = g⁻¹(a)

    @param f_AB dict {a: f(a)} für die Injektion f: A→B.
    @param f_BA dict {b: g(b)} für die Injektion g: B→A.
    @param A    Definitionsmenge.
    @param B    Zielmenge.
    @return     dict {a: h(a)} – die konstruierte Bijektion h: A→B.
    @lastModified 2026-03-10
    """
    # g⁻¹: Partielle Umkehrabbildung von g (nur auf g(B) definiert)
    g_inv = {v: k for k, v in f_BA.items()}

    # Berechne die "Schroeder-Bernstein-Kette" C = ⋃ A_n
    # A_0 = A \\ g(B) = Elemente von A, die kein Urbild unter g haben
    gB = set(f_BA.values())  # Bild von g in A
    A_set = set(A)

    # Starte mit A_0
    chain = A_set - gB
    C = set(chain)

    # Iteriere: A_{n+1} = g(f(A_n)) – wende f dann g an
    while True:
        # Wende f auf aktuelle Kette an, dann g
        new_elements = set()
        for a in chain:
            if a in f_AB:
                b = f_AB[a]
                if b in f_BA:
                    a_new = f_BA[b]
                    if a_new not in C:
                        new_elements.add(a_new)
        if not new_elements:
            break
        C |= new_elements
        chain = new_elements

    # Konstruiere Bijektion h
    h = {}
    for a in A:
        if a in C:
            # a liegt in der Kette → benutze f
            h[a] = f_AB[a]
        else:
            # a nicht in Kette → benutze g⁻¹
            if a in g_inv:
                h[a] = g_inv[a]
            else:
                # Fallback: f benutzen (sollte nicht vorkommen bei korrekten Injektionen)
                h[a] = f_AB[a]

    return h


def is_countable(description: str) -> dict:
    """
    Gibt eine Erklärung zurück, warum eine klassische Menge abzählbar oder
    überabzählbar ist.

    Unterstützte Beschreibungen: 'N', 'Z', 'Q', 'R', 'R^2', 'P(N)'.

    @param description Name der Menge ('N', 'Z', 'Q', 'R', 'R^2', 'P(N)').
    @return            dict mit 'set', 'countable', 'cardinality', 'proof_sketch'.
    @lastModified 2026-03-10
    """
    explanations = {
        'N': {
            'set': 'ℕ (natürliche Zahlen)',
            'countable': True,
            'cardinality': 'ℵ₀',
            'proof_sketch': (
                "ℕ ist per Definition abzählbar. Die Identität n↦n liefert "
                "die Bijektion ℕ→ℕ. |ℕ| = ℵ₀ (kleinstes unendliche Kardinalzahl)."
            )
        },
        'Z': {
            'set': 'ℤ (ganze Zahlen)',
            'countable': True,
            'cardinality': 'ℵ₀',
            'proof_sketch': (
                "Bijektion: 0↦0, 1↦1, -1↦2, 2↦3, -2↦4, … (alternierend).\n"
                "Formal: f(n) = 2n für n≥0, f(n) = -2n-1 für n<0."
            )
        },
        'Q': {
            'set': 'ℚ (rationale Zahlen)',
            'countable': True,
            'cardinality': 'ℵ₀',
            'proof_sketch': (
                "Cantors Diagonalzählung: Ordne p/q in einem 2D-Gitter an "
                "(Zeile=p, Spalte=q) und zähle diagonal durch. "
                "Alle Brüche werden erreicht → ℚ abzählbar."
            )
        },
        'R': {
            'set': 'ℝ (reelle Zahlen)',
            'countable': False,
            'cardinality': '2^ℵ₀ = ℶ₁ (überabzählbar)',
            'proof_sketch': (
                "Cantors Diagonalargument: Angenommen ℝ wäre abzählbar. "
                "Dann gäbe es eine Aufzählung r₁, r₂, r₃, … aller reellen Zahlen. "
                "Konstruiere x: x_n ≠ n-te Dezimalstelle von r_n. "
                "Dann x ≠ r_n für alle n – Widerspruch!"
            )
        },
        'R^2': {
            'set': 'ℝ² (Ebene)',
            'countable': False,
            'cardinality': '2^ℵ₀ = ℶ₁',
            'proof_sketch': (
                "Es gilt |ℝ²| = |ℝ|·|ℝ| = (2^ℵ₀)² = 2^ℵ₀. "
                "Bijektion ℝ²→ℝ via verschränkte Dezimalentwicklung: "
                "(0.a₁a₂…, 0.b₁b₂…) ↦ 0.a₁b₁a₂b₂…"
            )
        },
        'P(N)': {
            'set': 'P(ℕ) (Potenzmenge der natürlichen Zahlen)',
            'countable': False,
            'cardinality': '2^ℵ₀ = ℶ₁',
            'proof_sketch': (
                "Cantors Satz: |P(A)| > |A| für jede Menge A. "
                "Für A=ℕ: |P(ℕ)| = 2^ℵ₀ > ℵ₀. "
                "Bijektion P(ℕ)→{0,1}^ℕ→ℝ via charakteristische Funktionen."
            )
        }
    }

    if description not in explanations:
        return {
            'set': description,
            'countable': None,
            'cardinality': 'unbekannt',
            'proof_sketch': f"Keine Informationen für '{description}' vorhanden. "
                           "Unterstützt: N, Z, Q, R, R^2, P(N)."
        }

    return explanations[description]


def cantors_diagonal_argument(n: int = 5) -> dict:
    """
    Demonstriert Cantors Diagonalargument: |ℕ| < |P(ℕ)| = |ℝ|.

    Zeigt explizit ein binäres Wort, das nicht in einer hypothetischen
    Aufzählung aller Teilmengen von {0,…,n-1} enthalten sein kann.

    Vorgehen:
    1. Erzeuge n zufällige (deterministische) Teilmengen S₀,…,S_{n-1} von {0,…,n-1}
    2. Bilde Diagonale: D_i = 1 wenn i∉S_i, sonst D_i = 0
    3. D unterscheidet sich von jeder S_i an der i-ten Stelle

    @param n Anzahl der Beispielmengen (Standard: 5).
    @return  dict mit 'enumeration', 'diagonal_set', 'explanation'.
    @lastModified 2026-03-10
    """
    import random
    random.seed(42)  # Reproduzierbarkeit

    universe = list(range(n))

    # Erzeuge n deterministische Beispielmengen
    enumeration = []
    for i in range(n):
        # Jede Menge ist durch i*7 mod 2^n bestimmt (deterministisch)
        bits = (i * 7) % (2 ** n)
        subset = frozenset(j for j in universe if bits & (1 << j))
        enumeration.append(subset)

    # Diagonalargument: D = {i : i ∉ S_i}
    diagonal_set = frozenset(i for i in universe if i not in enumeration[i])

    # Erkläre warum D nicht in der Aufzählung ist
    differences = {}
    for i, S in enumerate(enumeration):
        if i in diagonal_set:
            differences[i] = f"D enthält {i}, aber S_{i} enthält {i} nicht → D ≠ S_{i}"
        else:
            differences[i] = f"D enthält {i} nicht, aber S_{i} enthält {i} → D ≠ S_{i}"

    return {
        'universe': set(universe),
        'enumeration': {f'S_{i}': set(S) for i, S in enumerate(enumeration)},
        'diagonal_set': set(diagonal_set),
        'differences': differences,
        'explanation': (
            f"Cantors Diagonalargument mit {n} Mengen:\n"
            f"D = {{i ∈ {{0,…,{n-1}}} : i ∉ S_i}} = {set(diagonal_set)}\n"
            f"D unterscheidet sich von jeder S_i an der Stelle i.\n"
            f"Daher: Keine Aufzählung kann alle Teilmengen von ℕ erfassen.\n"
            f"Folgerung: |P(ℕ)| > |ℕ|, also |ℝ| > ℵ₀."
        )
    }


def beth_numbers(n: int) -> str:
    """
    Beth-Zahlen ℶ_n (beth numbers).

    Definition:
    - ℶ_0 = ℵ_0 (Mächtigkeit von ℕ)
    - ℶ_{n+1} = 2^{ℶ_n} (Mächtigkeit der Potenzmenge)

    @param n Index der Beth-Zahl (n ≥ 0).
    @return  Beschreibung von ℶ_n als String.
    @raises  ValueError wenn n < 0.
    @lastModified 2026-03-10
    """
    if n < 0:
        raise ValueError("Beth-Zahl nur für n ≥ 0 definiert.")

    descriptions = {
        0: "ℶ_0 = ℵ_0 = |ℕ| (abzählbar unendlich)",
        1: "ℶ_1 = 2^ℵ_0 = |ℝ| = |P(ℕ)| (Mächtigkeit des Kontinuums)",
        2: "ℶ_2 = 2^(2^ℵ_0) = |P(ℝ)| = |Menge aller Funktionen ℝ→ℝ|",
        3: "ℶ_3 = 2^ℶ_2 = |P(P(ℝ))| (dritte Potenzmengenebene)",
    }

    if n in descriptions:
        return descriptions[n]
    return f"ℶ_{n} = 2^ℶ_{n-1} (Mächtigkeit der {n}-fachen Potenzmenge von ℕ)"


def aleph_numbers(n: int) -> str:
    """
    Aleph-Zahlen ℵ_n (aleph numbers).

    Beschreibt die n-te unendliche Kardinalzahl und ihre Beziehung zu ℶ_n.

    @param n Index der Aleph-Zahl (n ≥ 0).
    @return  Beschreibung von ℵ_n als String.
    @raises  ValueError wenn n < 0.
    @lastModified 2026-03-10
    """
    if n < 0:
        raise ValueError("Aleph-Zahl nur für n ≥ 0 definiert.")

    descriptions = {
        0: (
            "ℵ_0: Kleinste unendliche Kardinalzahl.\n"
            "= |ℕ| = |ℤ| = |ℚ| (alle abzählbaren Mengen).\n"
            "Kontinuumshypothese (CH): ℶ_1 = ℵ_1 (unbewiesen in ZFC)."
        ),
        1: (
            "ℵ_1: Nächste unendliche Kardinalzahl nach ℵ_0.\n"
            "= Mächtigkeit der Menge aller abzählbaren Ordinalzahlen.\n"
            "CH behauptet: ℵ_1 = 2^ℵ_0 = |ℝ| = ℶ_1."
        ),
        2: (
            "ℵ_2: Nächste unendliche Kardinalzahl nach ℵ_1.\n"
            "Generalisierte CH: ℵ_{n+1} = 2^ℵ_n für alle n."
        ),
    }

    if n in descriptions:
        return descriptions[n]
    return (
        f"ℵ_{n}: Die {n}-te unendliche Kardinalzahl in der Aleph-Hierarchie.\n"
        f"Beziehung zu Beth-Zahlen: ℶ_{n} ≥ ℵ_{n} (mit GCH: ℶ_{n} = ℵ_{n})."
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Ordinalzahlen
# ═══════════════════════════════════════════════════════════════════════════════

# Interne Repräsentation transfiniter Ordinalzahlen
# Wir kodieren endliche Ordinalzahlen als int und transfinite als Tupel:
# ω   →  ('omega',)
# ω+n →  ('omega', '+', n)
# ω·k →  ('omega', '*', k)
# ω²  →  ('omega', '^', 2)
# usw.

class Ordinal:
    """
    Ordinalzahl als Wohlordnungstyp.

    Unterstützt endliche Ordinalzahlen (als int) sowie transfinite Ordinalzahlen
    (ω, ω+1, ω+2, …, ω·2, ω², …) als symbolische Ausdrücke.

    Ordinaladdition ist NICHT kommutativ: 1+ω = ω ≠ ω+1.

    Darstellung:
        Ordinal(0), Ordinal(1), …, Ordinal('omega'), Ordinal('omega+1'), …

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    # Bekannte symbolische Konstanten
    OMEGA = 'omega'

    def __init__(self, value):
        """
        Initialisiert die Ordinalzahl.

        @param value int für endliche Ordinalzahlen; str für transfinite
                     ('omega', 'omega+1', 'omega*2', 'omega^2', usw.).
        """
        if isinstance(value, int):
            if value < 0:
                raise ValueError("Ordinalzahlen sind nicht-negativ.")
            self.value = value
            self.is_finite = True
        elif isinstance(value, str):
            self.value = value
            self.is_finite = False
        elif isinstance(value, Ordinal):
            self.value = value.value
            self.is_finite = value.is_finite
        else:
            raise TypeError(f"Ungültiger Ordinalzahl-Typ: {type(value)}")

    def successor(self) -> 'Ordinal':
        """
        Nachfolger-Ordinalzahl: α + 1.

        @return Nachfolger-Ordinalzahl.
        """
        if self.is_finite:
            return Ordinal(self.value + 1)
        # Transfinit: ω → ω+1, ω+n → ω+(n+1), etc.
        v = self.value
        if v == 'omega':
            return Ordinal('omega+1')
        elif v.startswith('omega+'):
            n = int(v[6:])
            return Ordinal(f'omega+{n + 1}')
        else:
            return Ordinal(v + '+1')

    def is_limit_ordinal(self) -> bool:
        """
        Prüft ob dies eine Limeszahl ist.

        Limeszahl: weder 0 noch Nachfolger-Ordinalzahl.
        Beispiele: ω, ω·2, ω², ε₀.

        @return True wenn Limeszahl.
        """
        if self.is_finite:
            return False  # Endliche Ordinalzahlen sind nie Limeszahlen
        v = self.value
        # ω, ω*k, ω^k sind Limeszahlen; ω+n für n>0 sind Nachfolger
        if v == 'omega':
            return True
        if v.startswith('omega+'):
            return False  # ω+n ist Nachfolger von ω+(n-1)
        if v.startswith('omega*') or v.startswith('omega^'):
            return True
        return True  # Im Zweifel transfinit = Limeszahl

    def is_successor_ordinal(self) -> bool:
        """
        Prüft ob dies eine Nachfolger-Ordinalzahl ist.

        @return True wenn Nachfolger-Ordinalzahl (also nicht 0 und nicht Limeszahl).
        """
        if self.is_finite:
            return self.value > 0
        return not self.is_limit_ordinal()

    def __lt__(self, other: 'Ordinal') -> bool:
        """
        Ordinalvergleich α < β.

        @param other Vergleichsordinalzahl β.
        @return      True wenn α < β.
        """
        if not isinstance(other, Ordinal):
            return NotImplemented

        # Endlich < Transfinit immer
        if self.is_finite and not other.is_finite:
            return True
        if not self.is_finite and other.is_finite:
            return False

        # Beide endlich: normaler Vergleich
        if self.is_finite and other.is_finite:
            return self.value < other.value

        # Beide transfinit: vereinfachte Ordnung
        order = ['omega', 'omega*2', 'omega*3', 'omega^2', 'omega^3']
        sv = self.value
        ov = other.value

        # ω+n < ω+m ⟺ n < m
        def parse_omega_n(s):
            """Extrahiert n aus 'omega+n', 0 für 'omega'."""
            if s == 'omega':
                return 0
            if s.startswith('omega+'):
                return int(s[6:])
            return None

        sn = parse_omega_n(sv)
        on = parse_omega_n(ov)

        if sn is not None and on is not None:
            return sn < on

        # Grobe Ordnung für andere transfinite Ausdrücke
        s_idx = order.index(sv) if sv in order else -1
        o_idx = order.index(ov) if ov in order else -1
        if s_idx >= 0 and o_idx >= 0:
            return s_idx < o_idx

        return str(sv) < str(ov)

    def __eq__(self, other: object) -> bool:
        """Gleichheit zweier Ordinalzahlen."""
        if isinstance(other, Ordinal):
            return self.value == other.value
        return NotImplemented

    def __le__(self, other: 'Ordinal') -> bool:
        """α ≤ β."""
        return self == other or self < other

    def __add__(self, other: 'Ordinal') -> 'Ordinal':
        """
        Ordinaladdition α + β (NICHT kommutativ!).

        Regeln:
        - α + 0 = α
        - n + ω = ω für endliches n (links von ω wird "geschluckt")
        - ω + n = ω+n
        - ω + ω = ω·2

        @param other Summand β.
        @return      Summe α + β als Ordinalzahl.
        """
        if not isinstance(other, Ordinal):
            return NotImplemented

        # α + 0 = α
        if other.is_finite and other.value == 0:
            return Ordinal(self.value)

        # 0 + β = β
        if self.is_finite and self.value == 0:
            return Ordinal(other.value)

        # Beide endlich
        if self.is_finite and other.is_finite:
            return Ordinal(self.value + other.value)

        # n + ω = ω (endliches n wird links geschluckt)
        if self.is_finite and other.value == 'omega':
            return Ordinal('omega')

        # ω + n = ω+n
        if self.value == 'omega' and other.is_finite:
            return Ordinal(f'omega+{other.value}')

        # ω+n + m = ω+(n+m)
        if isinstance(self.value, str) and self.value.startswith('omega+') and other.is_finite:
            n = int(self.value[6:])
            return Ordinal(f'omega+{n + other.value}')

        # ω + ω = ω·2
        if self.value == 'omega' and other.value == 'omega':
            return Ordinal('omega*2')

        # Allgemeiner Fall: symbolische Darstellung
        return Ordinal(f'{self.value}+{other.value}')

    def __mul__(self, other: 'Ordinal') -> 'Ordinal':
        """
        Ordinalmultiplikation α · β (NICHT kommutativ!).

        Regeln:
        - α · 0 = 0
        - n · ω = ω für endliches n > 0
        - ω · n = ω·n
        - ω · ω = ω²

        @param other Faktor β.
        @return      Produkt α · β.
        """
        if not isinstance(other, Ordinal):
            return NotImplemented

        # α · 0 = 0
        if other.is_finite and other.value == 0:
            return Ordinal(0)

        # α · 1 = α
        if other.is_finite and other.value == 1:
            return Ordinal(self.value)

        # 0 · β = 0
        if self.is_finite and self.value == 0:
            return Ordinal(0)

        # Beide endlich
        if self.is_finite and other.is_finite:
            return Ordinal(self.value * other.value)

        # n · ω = ω für n > 0 (endliches n links wird geschluckt)
        if self.is_finite and self.value > 0 and other.value == 'omega':
            return Ordinal('omega')

        # ω · n = ω*n
        if self.value == 'omega' and other.is_finite:
            return Ordinal(f'omega*{other.value}')

        # ω · ω = ω²
        if self.value == 'omega' and other.value == 'omega':
            return Ordinal('omega^2')

        return Ordinal(f'({self.value})*({other.value})')

    def __pow__(self, other: 'Ordinal') -> 'Ordinal':
        """
        Ordinalexponentiation α^β.

        Regeln:
        - α^0 = 1
        - α^1 = α
        - ω^n → 'omega^n'

        @param other Exponent β.
        @return      Potenz α^β.
        """
        if not isinstance(other, Ordinal):
            return NotImplemented

        if other.is_finite and other.value == 0:
            return Ordinal(1)
        if other.is_finite and other.value == 1:
            return Ordinal(self.value)
        if self.is_finite and other.is_finite:
            return Ordinal(self.value ** other.value)

        if self.value == 'omega' and other.is_finite:
            return Ordinal(f'omega^{other.value}')

        return Ordinal(f'({self.value})^({other.value})')

    def __repr__(self) -> str:
        """Lesbare Darstellung der Ordinalzahl."""
        if self.is_finite:
            return str(self.value)
        # Ersetze 'omega' durch 'ω' für schönere Ausgabe
        return self.value.replace('omega', 'ω')

    def __hash__(self) -> int:
        """Hashwert für Verwendung in Mengen."""
        return hash(self.value)


def is_well_ordered(relation: Relation, domain: MathSet) -> bool:
    """
    Prüft ob eine Relation eine Wohlordnung auf einer Menge ist.

    Wohlordnung: Totalordnung, bei der jede nicht-leere Teilmenge
    ein kleinstes Element hat.

    Für endliche Mengen: Totalordnung ⟹ Wohlordnung.

    @param relation Zu prüfende Relation.
    @param domain   Trägermenge.
    @return         True wenn Wohlordnung.
    @lastModified 2026-03-10
    """
    # Prüfe Totalordnung (notwendige Bedingung)
    if not relation.is_total_order():
        return False

    # Für endliche Mengen gilt: Totalordnung ⟹ Wohlordnung
    # (jede endliche nicht-leere Teilmenge hat ein kleinstes Element)
    # Für unendliche Mengen müsste dies tiefer geprüft werden
    return True


def ordinal_arithmetic_demo() -> dict:
    """
    Demonstriert die Nicht-Kommutativität der Ordinalzahlarithmetik.

    Zeigt: ω+1 ≠ 1+ω und ω·2 ≠ 2·ω.

    @return dict mit Beispielen und Erklärungen.
    @lastModified 2026-03-10
    """
    omega = Ordinal('omega')
    one = Ordinal(1)
    two = Ordinal(2)

    # Ordinaladdition
    omega_plus_1 = omega + one     # ω+1 (echter Nachfolger von ω)
    one_plus_omega = one + omega   # ω   (1 wird links geschluckt)

    # Ordinalmultiplikation
    omega_times_2 = omega * two    # ω·2 (zwei Kopien von ω)
    two_times_omega = two * omega  # ω   (endliches Vielfaches von ω)

    return {
        'omega+1': str(omega_plus_1),
        '1+omega': str(one_plus_omega),
        'omega+1 != 1+omega': omega_plus_1 != one_plus_omega,
        'omega*2': str(omega_times_2),
        '2*omega': str(two_times_omega),
        'omega*2 != 2*omega': omega_times_2 != two_times_omega,
        'explanation': (
            "Ordinaladdition ist NICHT kommutativ:\n"
            "• ω+1 = {0,1,2,…,ω} (Wohlordnungstyp mit größtem Element)\n"
            "• 1+ω = {0,1,2,…} = ω (das '1' links wird 'geschluckt')\n\n"
            "Ordinalmultiplikation ist NICHT kommutativ:\n"
            "• ω·2 = ω+ω = {0,1,…,ω,ω+1,…} (zwei Kopien von ω)\n"
            "• 2·ω = ω (endlich viele Kopien werden zu einer einzigen ω verschmolzen)"
        )
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Axiomensysteme
# ═══════════════════════════════════════════════════════════════════════════════

def zfc_axioms() -> dict:
    """
    Gibt die ZFC-Axiome (Zermelo-Fraenkel mit Auswahlaxiom) zurück.

    ZFC ist das Standard-Axiomensystem der modernen Mathematik.
    Es besteht aus 10 Axiomen (oder Axiomenschemata).

    @return dict {axiom_name: {'formal': ..., 'description': ...}}.
    @lastModified 2026-03-10
    """
    return {
        'extensionality': {
            'name': 'Extensionalitätsaxiom',
            'formal': '∀A ∀B [∀x(x∈A ↔ x∈B) → A=B]',
            'description': (
                "Zwei Mengen sind gleich, wenn sie dieselben Elemente enthalten. "
                "Die Identität einer Menge wird vollständig durch ihre Elemente bestimmt."
            )
        },
        'empty_set': {
            'name': 'Leermengenaxiom',
            'formal': '∃A ∀x ¬(x∈A)',
            'description': (
                "Es existiert eine Menge ohne Elemente: die leere Menge ∅. "
                "Aus dem Extensionalitätsaxiom folgt deren Eindeutigkeit."
            )
        },
        'pairing': {
            'name': 'Paarmengenaxiom',
            'formal': '∀a ∀b ∃P ∀x(x∈P ↔ x=a ∨ x=b)',
            'description': (
                "Für je zwei Mengen a, b existiert eine Menge {a,b}. "
                "Ermöglicht die Konstruktion geordneter Paare (a,b) = {{a},{a,b}}."
            )
        },
        'union': {
            'name': 'Vereinigungsaxiom',
            'formal': '∀F ∃U ∀x[x∈U ↔ ∃A(A∈F ∧ x∈A)]',
            'description': (
                "Für jede Mengenfamilie F existiert ihre Vereinigung ⋃F. "
                "Insbesondere: A∪B = ⋃{A,B}."
            )
        },
        'power_set': {
            'name': 'Potenzmengenaxiom',
            'formal': '∀A ∃P ∀B[B∈P ↔ B⊆A]',
            'description': (
                "Für jede Menge A existiert ihre Potenzmenge P(A) = {B : B⊆A}. "
                "Liefert die Grundlage für überabzählbare Mengen: |P(ℕ)| > |ℕ|."
            )
        },
        'infinity': {
            'name': 'Unendlichkeitsaxiom',
            'formal': '∃I[∅∈I ∧ ∀x∈I(x∪{x}∈I)]',
            'description': (
                "Es existiert eine induktive Menge I: ∅∈I und x∈I ⟹ x∪{x}∈I. "
                "Die minimale solche Menge entspricht den natürlichen Zahlen ℕ."
            )
        },
        'separation': {
            'name': 'Aussonderungsschema (Separation)',
            'formal': '∀A ∀p ∃S ∀x[x∈S ↔ x∈A ∧ φ(x,p)]',
            'description': (
                "Für jede Menge A und Eigenschaft φ existiert S = {x∈A : φ(x)}. "
                "Verhindert Russells Paradox (kein unbeschränktes Aussonderungsprinzip)."
            )
        },
        'replacement': {
            'name': 'Ersetzungsschema',
            'formal': '∀A[∀x∈A ∃!y φ(x,y) → ∃B ∀x∈A ∃y∈B φ(x,y)]',
            'description': (
                "Das Bild einer Menge unter einer (klassenmäßigen) Funktion ist eine Menge. "
                "Ermöglicht die Konstruktion von V_{ω+ω} und höheren Stufen der Von-Neumann-Hierarchie."
            )
        },
        'foundation': {
            'name': 'Fundierungsaxiom (Regularitätsaxiom)',
            'formal': '∀A[A≠∅ → ∃x∈A(x∩A=∅)]',
            'description': (
                "Jede nicht-leere Menge A enthält ein ∈-minimales Element. "
                "Verhindert zirkuläre Mengen (A∈A) und unendlich absteigende ∈-Ketten."
            )
        },
        'choice': {
            'name': 'Auswahlaxiom (Axiom of Choice, AC)',
            'formal': '∀F[∅∉F → ∃f: F→⋃F ∀A∈F(f(A)∈A)]',
            'description': (
                "Für jede Familie nicht-leerer Mengen existiert eine Auswahlfunktion. "
                "Äquivalent zu: Wohlordnungssatz, Lemma von Zorn, Tychonoff-Satz. "
                "Unabhängig von ZF (Gödel 1938: konsistent; Cohen 1963: unabhängig)."
            )
        }
    }


def axiom_of_choice_equivalents() -> list:
    """
    Listet wichtige Äquivalente zum Auswahlaxiom (AC) auf.

    Alle folgenden Aussagen sind in ZF äquivalent zu AC.

    @return Liste von dicts mit 'name' und 'description'.
    @lastModified 2026-03-10
    """
    return [
        {
            'name': "Lemma von Zorn",
            'description': (
                "Jede nicht-leere geordnete Menge, in der jede total geordnete "
                "Teilkette eine obere Schranke hat, enthält ein maximales Element."
            )
        },
        {
            'name': "Wohlordnungssatz (Zermelo)",
            'description': (
                "Jede Menge kann wohlgeordnet werden. "
                "Insbesondere: ℝ kann wohlgeordnet werden (keine explizite Konstruktion bekannt)."
            )
        },
        {
            'name': "Tychonoff-Satz",
            'description': (
                "Jedes Produkt kompakter topologischer Räume ist kompakt "
                "(in der Produkttopologie)."
            )
        },
        {
            'name': "Cantor-Bernstein-Schröder ist ohne AC beweisbar",
            'description': (
                "Der CSB-Satz benötigt kein AC. "
                "Aber: Trichotomie der Kardinalzahlen (|A|<|B| oder |A|=|B| oder |A|>|B|) "
                "ist äquivalent zu AC."
            )
        },
        {
            'name': "Basis-Satz (Hamel-Basis)",
            'description': (
                "Jeder Vektorraum über einem Körper hat eine Basis (Hamel-Basis). "
                "Äquivalent zu AC."
            )
        },
        {
            'name': "Satz vom maximalen Ideal",
            'description': (
                "Jeder echte Ideal in einem Ring mit 1 liegt in einem maximalen Ideal."
            )
        },
        {
            'name': "Ultrafilter-Lemma",
            'description': (
                "Jeder Filter auf einer Menge kann zu einem Ultrafilter erweitert werden. "
                "(Schwächer als AC, äquivalent zu Tychonoff für Hausdorff-Räume)."
            )
        }
    ]


def continuum_hypothesis() -> dict:
    """
    Erklärt die Kontinuumshypothese (KH / CH).

    Die KH besagt: 2^ℵ₀ = ℵ₁, d.h. es gibt keine Kardinalzahl κ mit
    ℵ₀ < κ < 2^ℵ₀ = |ℝ|.

    Gödel (1938): KH ist konsistent mit ZFC.
    Cohen (1963): ¬KH ist konsistent mit ZFC.
    Folge: KH ist unentscheidbar in ZFC (unabhängig).

    @return dict mit 'statement', 'formal', 'independence', 'generalized'.
    @lastModified 2026-03-10
    """
    return {
        'statement': (
            "Kontinuumshypothese (CH): Es gibt keine Kardinalzahl κ mit "
            "ℵ₀ < κ < 2^ℵ₀. Anders: |ℝ| = ℵ₁."
        ),
        'formal': '2^ℵ₀ = ℵ₁',
        'independence': (
            "Kurt Gödel (1938): CH ist konsistent mit ZFC (konstruierbare Universen L). "
            "Paul Cohen (1963): ¬CH ist konsistent mit ZFC (Forcing-Methode). "
            "Folge: CH ist in ZFC weder beweisbar noch widerlegbar – sie ist unabhängig."
        ),
        'generalized': (
            "Verallgemeinerte Kontinuumshypothese (GCH): 2^ℵ_α = ℵ_{α+1} für alle α. "
            "Äquivalent: ℶ_α = ℵ_α für alle Ordinalzahlen α. "
            "GCH impliziert AC und ist ebenfalls unabhängig von ZFC."
        ),
        'beth_aleph_relation': (
            "Unter GCH gilt ℶ_n = ℵ_n: "
            "ℶ_0 = ℵ_0, ℶ_1 = ℵ_1, ℶ_2 = ℵ_2, … "
            "Ohne GCH gilt nur: ℵ_n ≤ ℶ_n."
        )
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Mengenfamilien und σ-Algebren
# ═══════════════════════════════════════════════════════════════════════════════

class SetFamily:
    """
    Familie von Mengen {A_i : i∈I} mit Indexmenge I.

    Eine Mengenfamilie ist eine Abbildung I → P(U) von einer Indexmenge
    in die Potenzmenge des Universums.

    Unterstützt: Vereinigung, Schnitt, Partition-Prüfung, σ-Algebra-Prüfung.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, sets: dict):
        """
        Initialisiert die Mengenfamilie.

        @param sets dict {index: MathSet} – Indexierte Familie von Mengen.
        """
        self.sets = {k: (v if isinstance(v, MathSet) else MathSet(v))
                     for k, v in sets.items()}

    def union(self) -> MathSet:
        """
        Vereinigung aller Mengen: ⋃_{i∈I} A_i.

        @return Vereinigung aller indizierten Mengen.
        """
        if not self.sets:
            return MathSet()
        result = set()
        for A in self.sets.values():
            result |= set(A)
        return MathSet(result)

    def intersection(self) -> MathSet:
        """
        Schnitt aller Mengen: ⋂_{i∈I} A_i.

        @return Schnitt aller indizierten Mengen (∅ wenn Familie leer).
        """
        if not self.sets:
            return MathSet()
        it = iter(self.sets.values())
        result = set(next(it))
        for A in it:
            result &= set(A)
        return MathSet(result)

    def is_partition(self, universe: MathSet) -> bool:
        """
        Prüft ob die Familie eine Partition des Universums ist.

        Partition: Überdeckung + paarweise disjunkt + keine leeren Mengen.

        @param universe Das Universum U.
        @return         True wenn Partition von U.
        """
        # 1. Keine leeren Mengen
        if any(len(A) == 0 for A in self.sets.values()):
            return False

        # 2. Überdeckung: ⋃ A_i = U
        if self.union() != universe:
            return False

        # 3. Paarweise disjunkt
        sets_list = list(self.sets.values())
        for i in range(len(sets_list)):
            for j in range(i + 1, len(sets_list)):
                if not sets_list[i].is_disjoint(sets_list[j]):
                    return False

        return True

    def is_cover(self, universe: MathSet) -> bool:
        """
        Prüft ob die Familie eine Überdeckung des Universums ist.

        Überdeckung: ⋃_{i∈I} A_i ⊇ U (d.h. = U wenn A_i ⊆ U).

        @param universe Das Universum U.
        @return         True wenn Überdeckung.
        """
        return universe <= self.union()

    def is_sigma_algebra(self, universe: MathSet) -> bool:
        """
        Prüft ob die Familie eine σ-Algebra auf dem Universum ist.

        σ-Algebra Σ über U erfüllt:
        1. U ∈ Σ
        2. A ∈ Σ ⟹ A^c ∈ Σ (abgeschlossen unter Komplementbildung)
        3. A₁,A₂,… ∈ Σ ⟹ ⋃A_n ∈ Σ (abgeschlossen unter abzählbarer Vereinigung)

        Für endliche Familien: abzählbare Vereinigung = endliche Vereinigung.

        @param universe Das Universum U.
        @return         True wenn σ-Algebra.
        """
        sets_list = list(self.sets.values())

        # Eigenschaft 1: U ∈ Σ
        if universe not in sets_list:
            return False

        # Eigenschaft 2: Abgeschlossenheit unter Komplement
        for A in sets_list:
            A_complement = A.complement(universe)
            if A_complement not in sets_list:
                return False

        # Eigenschaft 3: Abgeschlossenheit unter endlicher Vereinigung
        # (für endliche Familien ausreichend als Approximation)
        for i, A in enumerate(sets_list):
            for j, B in enumerate(sets_list):
                A_union_B = A.union(B)
                if A_union_B not in sets_list:
                    return False

        return True

    def generated_sigma_algebra(self, universe: MathSet) -> 'SetFamily':
        """
        Erzeugt die von der Familie erzeugte σ-Algebra σ(F).

        σ(F) ist die kleinste σ-Algebra, die alle Mengen aus F enthält.

        Algorithmus: Schließe iterativ unter Komplement und endlicher Vereinigung ab.

        @param universe Das Universum U.
        @return         SetFamily, die σ(F) repräsentiert.
        """
        # Starte mit ∅, U und allen gegebenen Mengen
        generated = {MathSet(), universe}
        for A in self.sets.values():
            generated.add(A)

        # Iterativer Abschluss unter Komplement und Vereinigung
        changed = True
        while changed:
            changed = False
            current = list(generated)

            # Komplement-Abschluss
            for A in current:
                comp = A.complement(universe)
                if comp not in generated:
                    generated.add(comp)
                    changed = True

            # Vereinigungs-Abschluss
            current = list(generated)
            for A in current:
                for B in current:
                    union = A.union(B)
                    if union not in generated:
                        generated.add(union)
                        changed = True

        # Erstelle indizierte Familie
        return SetFamily({i: A for i, A in enumerate(generated)})

    def __repr__(self) -> str:
        """Lesbare Darstellung der Mengenfamilie."""
        return f"SetFamily({{{', '.join(f'{k}: {v}' for k, v in self.sets.items())}}})"


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Filter
# ═══════════════════════════════════════════════════════════════════════════════

class Filter:
    """
    Mengentheoretischer Filter auf einer Menge.

    Ein Filter F auf U ist eine nicht-leere Familie von Teilmengen von U mit:
    1. ∅ ∉ F  (Nicht-Trivialität)
    2. A∈F, A⊆B ⊆ U ⟹ B∈F  (Aufwärtsabgeschlossenheit)
    3. A,B∈F ⟹ A∩B∈F  (Schnitt-Abgeschlossenheit)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, universe: MathSet, sets: list[MathSet]):
        """
        Initialisiert den Filter.

        @param universe Das Universum U.
        @param sets     Liste von Teilmengen von U, die den Filter bilden sollen.
        """
        self.universe = universe
        self.sets = list(sets)

    def is_filter(self) -> bool:
        """
        Prüft ob die gegebenen Mengen einen Filter bilden.

        Bedingungen:
        1. ∅ ∉ F
        2. Aufwärtsabgeschlossenheit: A∈F, A⊆B ⊆ U ⟹ B∈F
        3. Schnitt-Abgeschlossenheit: A,B∈F ⟹ A∩B∈F

        @return True wenn Filtereigenschaften erfüllt.
        """
        # Bedingung 1: ∅ nicht in F
        empty = MathSet()
        if empty in self.sets:
            return False

        # Bedingung 3: Schnitt-Abgeschlossenheit (einfacher zu prüfen)
        for A in self.sets:
            for B in self.sets:
                intersection = A.intersection(B)
                if intersection not in self.sets:
                    return False

        # Bedingung 2: Aufwärtsabgeschlossenheit
        # Prüfe: A∈F und A⊆B⊆U ⟹ B∈F
        # Für endliche Universen: prüfe alle Teilmengen B mit A⊆B
        for A in self.sets:
            # Erzeuge alle Obermengen von A, die Teilmengen von U sind
            u_elements = list(self.universe)
            a_elements = set(A)
            extra = [x for x in u_elements if x not in a_elements]
            for r in range(len(extra) + 1):
                for combo in itertools.combinations(extra, r):
                    B = MathSet(a_elements | set(combo))
                    if B not in self.sets:
                        return False

        return True

    def is_ultrafilter(self) -> bool:
        """
        Prüft ob der Filter ein Ultrafilter ist.

        Ultrafilter: Maximaler Filter, d.h. für jede Teilmenge A⊆U gilt:
        entweder A∈F oder A^c∈F (aber nicht beide, da F kein ∅ enthält).

        @return True wenn Ultrafilter.
        """
        if not self.is_filter():
            return False

        # Ultrafilter-Bedingung: ∀A⊆U: A∈F oder (U\A)∈F
        u_elements = list(self.universe)
        n = len(u_elements)
        for r in range(n + 1):
            for combo in itertools.combinations(u_elements, r):
                A = MathSet(combo)
                A_comp = A.complement(self.universe)
                if A not in self.sets and A_comp not in self.sets:
                    return False

        return True

    def __repr__(self) -> str:
        """Lesbare Darstellung des Filters."""
        return f"Filter(universe={self.universe}, sets={self.sets})"


def principal_filter(base_set: MathSet, universe: MathSet) -> Filter:
    """
    Konstruiert den Hauptfilter (Principal Filter) erzeugt durch base_set B.

    Hauptfilter: F_B = {A ⊆ U : B ⊆ A}
    (Alle Obermengenvon B innerhalb von U)

    @param base_set Die erzeugende Menge B.
    @param universe Das Universum U.
    @return         Hauptfilter F_B als Filter-Objekt.
    @lastModified 2026-03-10
    """
    u_elements = list(universe)
    b_elements = set(base_set)
    extra = [x for x in u_elements if x not in b_elements]

    # Alle Obermengen von B in U
    filter_sets = []
    for r in range(len(extra) + 1):
        for combo in itertools.combinations(extra, r):
            A = MathSet(b_elements | set(combo))
            filter_sets.append(A)

    return Filter(universe, filter_sets)
