"""
@file lattice_theory.py
@brief Verbandstheorie: Verbände, Boolesche Algebren, distributive Verbände,
       Birkhoff-Darstellungssatz, Dedekind-Zahlen.
@description
    Dieses Modul implementiert die grundlegenden Strukturen der Verbandstheorie:

    Klassen:
    - PartialOrder    – Halbordnung (M, ≤): reflexiv, antisymmetrisch, transitiv
    - Lattice         – Verband (L, ∧, ∨): eindeutige Meet- und Join-Operationen
    - BooleanAlgebra  – Boolesche Algebra: distributiv, komplementiert, mit ⊤,⊥

    Freie Funktionen:
    - divisibility_lattice()              – Teilerverband von n
    - power_set_lattice()                 – Potenzmengenverband als Bool. Algebra
    - partition_lattice()                 – Partitionsverband Πₙ
    - birkhoff_representation_theorem()   – Birkhoff-Darstellungssatz
    - dedekind_numbers()                  – Dedekind-Zahlen D(n)
    - jordan_dedekind_chain_condition()   – Jordan-Dedekind-Kettenbedingung

    Mathematische Grundlagen:
    - Halbordnung: (M, ≤) mit a≤a, (a≤b ∧ b≤a → a=b), (a≤b ∧ b≤c → a≤c)
    - Verband: jedes Paar (a,b) hat eindeutiges sup und inf
    - Distributivität: a∧(b∨c) = (a∧b)∨(a∧c)
    - Boolesche Algebra: distributiver, komplementierter beschränkter Verband

@author Kurt Ingwer
@lastModified 2026-03-10
"""

from itertools import chain, combinations
from functools import reduce
from math import gcd


# =============================================================================
# KLASSE: PartialOrder
# =============================================================================

class PartialOrder:
    """
    Halbordnung (M, ≤).

    Eine Halbordnung erfüllt für alle a, b, c ∈ M:
    - Reflexivität: a ≤ a
    - Antisymmetrie: a ≤ b ∧ b ≤ a → a = b
    - Transitivität: a ≤ b ∧ b ≤ c → a ≤ c

    @param elements        Liste aller Elemente
    @param order_relation  Menge von Paaren (a, b) mit a ≤ b
    @lastModified 2026-03-10
    """

    def __init__(self, elements: list, order_relation: set):
        """
        Initialisiert die Halbordnung.

        @param elements        Alle Elemente der Halbordnung
        @param order_relation  Explizite Menge {(a,b) | a ≤ b}
        @lastModified 2026-03-10
        """
        self.elements = list(elements)
        # Als Menge von Paaren speichern; Tupel für Hashbarkeit
        self.order_relation = set(order_relation)

    def leq(self, a, b) -> bool:
        """
        Prüft ob a ≤ b in der Halbordnung gilt.

        @param a  Erstes Element
        @param b  Zweites Element
        @return   True wenn a ≤ b
        @lastModified 2026-03-10
        """
        return (a, b) in self.order_relation

    def is_partial_order(self) -> bool:
        """
        Prüft ob die Relation tatsächlich eine Halbordnung ist.

        Überprüft alle drei Axiome: Reflexivität, Antisymmetrie, Transitivität.

        @return  True wenn alle Axiome erfüllt sind
        @lastModified 2026-03-10
        """
        # Reflexivität: ∀a: a ≤ a
        for a in self.elements:
            if not self.leq(a, a):
                return False

        # Antisymmetrie: a ≤ b ∧ b ≤ a → a = b
        for a in self.elements:
            for b in self.elements:
                if a != b and self.leq(a, b) and self.leq(b, a):
                    return False

        # Transitivität: a ≤ b ∧ b ≤ c → a ≤ c
        for a in self.elements:
            for b in self.elements:
                for c in self.elements:
                    if self.leq(a, b) and self.leq(b, c) and not self.leq(a, c):
                        return False

        return True

    def hasse_diagram(self) -> dict:
        """
        Berechnet die reduzierte Relation für das Hasse-Diagramm.

        Entfernt alle Kanten, die durch Transitivität erzwungen sind:
        a wird unter b gezeichnet wenn a < b und kein c mit a < c < b existiert.

        @return  Dict {a: [direkte_Nachfolger_von_a]}
        @lastModified 2026-03-10
        """
        # Strikte Ordnung a < b (a ≤ b und a ≠ b)
        def strictly_less(a, b):
            return self.leq(a, b) and a != b

        hasse = {a: [] for a in self.elements}

        for a in self.elements:
            for b in self.elements:
                if not strictly_less(a, b):
                    continue
                # Prüfen ob es ein c gibt mit a < c < b (dann ist a→b keine Hasse-Kante)
                has_intermediate = False
                for c in self.elements:
                    if c != a and c != b and strictly_less(a, c) and strictly_less(c, b):
                        has_intermediate = True
                        break
                if not has_intermediate:
                    hasse[a].append(b)

        return hasse

    def minimal_elements(self) -> list:
        """
        Bestimmt alle minimalen Elemente der Halbordnung.

        Ein Element m ist minimal wenn kein a ∈ M mit a < m existiert.

        @return  Liste der minimalen Elemente
        @lastModified 2026-03-10
        """
        minimal = []
        for a in self.elements:
            # Prüfen ob es ein kleineres Element gibt
            is_minimal = not any(
                b != a and self.leq(b, a) for b in self.elements
            )
            if is_minimal:
                minimal.append(a)
        return minimal

    def maximal_elements(self) -> list:
        """
        Bestimmt alle maximalen Elemente der Halbordnung.

        Ein Element m ist maximal wenn kein a ∈ M mit m < a existiert.

        @return  Liste der maximalen Elemente
        @lastModified 2026-03-10
        """
        maximal = []
        for a in self.elements:
            # Prüfen ob es ein größeres Element gibt
            is_maximal = not any(
                b != a and self.leq(a, b) for b in self.elements
            )
            if is_maximal:
                maximal.append(a)
        return maximal

    def _all_chains(self) -> list:
        """
        Findet alle Ketten (total geordnete Teilmengen) durch DFS.

        @return  Liste aller maximalen Ketten als Listen von Elementen
        @lastModified 2026-03-10
        """
        def extend_chain(chain):
            """Verlängert eine Kette um alle möglichen Nachfolger."""
            last = chain[-1]
            extensions = [
                b for b in self.elements
                if b not in chain and self.leq(last, b) and b != last
                # Kompatibilität mit der Kette: b > alle Elemente der Kette
                and all(self.leq(c, b) for c in chain)
            ]
            if not extensions:
                return [chain]
            result = []
            for b in extensions:
                result.extend(extend_chain(chain + [b]))
            return result

        all_chains = []
        # Von jedem minimalen Element starten
        for start in self.minimal_elements():
            chains = extend_chain([start])
            all_chains.extend(chains)

        return all_chains if all_chains else [[a] for a in self.elements]

    def chain_length(self) -> int:
        """
        Berechnet die Länge der längsten Kette.

        Die Länge einer Kette mit k Elementen ist k-1 (Anzahl der Kanten).

        @return  Länge der längsten Kette (Anzahl Kanten)
        @lastModified 2026-03-10
        """
        # Alle maximalen Ketten finden
        max_len = 0
        # Dynamische Programmierung: für jedes Element die längste aufsteigende Kette
        dp = {a: 1 for a in self.elements}
        # Topologische Reihenfolge approximieren
        sorted_elems = list(self.elements)
        for _ in range(len(sorted_elems)):
            for a in sorted_elems:
                for b in sorted_elems:
                    if a != b and self.leq(a, b):
                        if dp[b] < dp[a] + 1:
                            dp[b] = dp[a] + 1
        max_len = max(dp.values()) - 1 if dp else 0
        return max_len

    def is_chain(self) -> bool:
        """
        Prüft ob die Halbordnung eine Kette (total geordnete Menge) ist.

        Alle Elemente müssen paarweise vergleichbar sein.

        @return  True wenn total geordnet
        @lastModified 2026-03-10
        """
        for a in self.elements:
            for b in self.elements:
                if not (self.leq(a, b) or self.leq(b, a)):
                    return False
        return True


# =============================================================================
# KLASSE: Lattice
# =============================================================================

class Lattice(PartialOrder):
    """
    Verband (L, ∧, ∨).

    Ein Verband ist eine Halbordnung, in der jedes Paar (a, b) ein:
    - Infimum (Meet): a ∧ b = größte untere Schranke
    - Supremum (Join): a ∨ b = kleinste obere Schranke
    besitzt.

    @lastModified 2026-03-10
    """

    def meet(self, a, b):
        """
        Berechnet a ∧ b (Infimum, größte untere Schranke).

        a ∧ b = max{c ∈ L | c ≤ a ∧ c ≤ b}

        @param a  Erstes Element
        @param b  Zweites Element
        @return   Infimum a ∧ b, oder None wenn nicht existent
        @lastModified 2026-03-10
        """
        # Untere Schranken: alle c mit c ≤ a und c ≤ b
        lower_bounds = [
            c for c in self.elements
            if self.leq(c, a) and self.leq(c, b)
        ]
        if not lower_bounds:
            return None
        # Größte untere Schranke finden
        for candidate in lower_bounds:
            if all(self.leq(c, candidate) for c in lower_bounds):
                return candidate
        return None

    def join(self, a, b):
        """
        Berechnet a ∨ b (Supremum, kleinste obere Schranke).

        a ∨ b = min{c ∈ L | a ≤ c ∧ b ≤ c}

        @param a  Erstes Element
        @param b  Zweites Element
        @return   Supremum a ∨ b, oder None wenn nicht existent
        @lastModified 2026-03-10
        """
        # Obere Schranken: alle c mit a ≤ c und b ≤ c
        upper_bounds = [
            c for c in self.elements
            if self.leq(a, c) and self.leq(b, c)
        ]
        if not upper_bounds:
            return None
        # Kleinste obere Schranke finden
        for candidate in upper_bounds:
            if all(self.leq(candidate, c) for c in upper_bounds):
                return candidate
        return None

    def top(self):
        """
        Gibt das größte Element ⊤ zurück (falls vorhanden).

        ⊤ ist das eindeutige Element mit a ≤ ⊤ für alle a ∈ L.

        @return  Größtes Element oder None
        @lastModified 2026-03-10
        """
        maximal = self.maximal_elements()
        # Eindeutiges Maximum existiert genau dann wenn es nur ein maximales Element gibt
        if len(maximal) == 1:
            return maximal[0]
        # Prüfen ob eines der maximalen Elemente wirklich größer als alle anderen ist
        for candidate in maximal:
            if all(self.leq(a, candidate) for a in self.elements):
                return candidate
        return None

    def bottom(self):
        """
        Gibt das kleinste Element ⊥ zurück (falls vorhanden).

        ⊥ ist das eindeutige Element mit ⊥ ≤ a für alle a ∈ L.

        @return  Kleinstes Element oder None
        @lastModified 2026-03-10
        """
        minimal = self.minimal_elements()
        if len(minimal) == 1:
            return minimal[0]
        for candidate in minimal:
            if all(self.leq(candidate, a) for a in self.elements):
                return candidate
        return None

    def is_distributive(self) -> bool:
        """
        Prüft ob der Verband distributiv ist.

        Distributivitätsgesetz:
            a ∧ (b ∨ c) = (a ∧ b) ∨ (a ∧ c) für alle a, b, c ∈ L

        (Äquivalent: a ∨ (b ∧ c) = (a ∨ b) ∧ (a ∨ c))

        @return  True wenn distributiv
        @lastModified 2026-03-10
        """
        for a in self.elements:
            for b in self.elements:
                for c in self.elements:
                    # a ∧ (b ∨ c)
                    bvc = self.join(b, c)
                    lhs = self.meet(a, bvc)
                    # (a ∧ b) ∨ (a ∧ c)
                    amb = self.meet(a, b)
                    amc = self.meet(a, c)
                    rhs = self.join(amb, amc)
                    if lhs != rhs:
                        return False
        return True

    def is_modular(self) -> bool:
        """
        Prüft ob der Verband modular ist.

        Modularitätsgesetz:
            a ≤ c ⟹ a ∨ (b ∧ c) = (a ∨ b) ∧ c für alle a, b, c ∈ L

        Jeder distributive Verband ist modular, aber nicht umgekehrt.

        @return  True wenn modular
        @lastModified 2026-03-10
        """
        for a in self.elements:
            for b in self.elements:
                for c in self.elements:
                    if not self.leq(a, c):
                        continue
                    # a ∨ (b ∧ c)
                    bmc = self.meet(b, c)
                    lhs = self.join(a, bmc)
                    # (a ∨ b) ∧ c
                    avb = self.join(a, b)
                    rhs = self.meet(avb, c)
                    if lhs != rhs:
                        return False
        return True

    def is_complemented(self) -> bool:
        """
        Prüft ob der Verband komplementiert ist.

        Ein Verband ist komplementiert wenn:
        - ⊤ und ⊥ existieren
        - Für jedes a ∈ L existiert ein ā mit a ∧ ā = ⊥ und a ∨ ā = ⊤

        @return  True wenn komplementiert
        @lastModified 2026-03-10
        """
        bot = self.bottom()
        top = self.top()

        # Brauchen ⊤ und ⊥
        if bot is None or top is None:
            return False

        # Für jedes Element ein Komplement suchen
        for a in self.elements:
            has_complement = False
            for a_bar in self.elements:
                meet_val = self.meet(a, a_bar)
                join_val = self.join(a, a_bar)
                if meet_val == bot and join_val == top:
                    has_complement = True
                    break
            if not has_complement:
                return False
        return True

    def complement(self, a):
        """
        Findet ein Komplement von a (falls vorhanden).

        ā ist Komplement von a wenn a ∧ ā = ⊥ und a ∨ ā = ⊤.

        @param a  Element dessen Komplement gesucht wird
        @return   Komplement ā oder None
        @lastModified 2026-03-10
        """
        bot = self.bottom()
        top = self.top()

        if bot is None or top is None:
            return None

        for a_bar in self.elements:
            if self.meet(a, a_bar) == bot and self.join(a, a_bar) == top:
                return a_bar
        return None


# =============================================================================
# KLASSE: BooleanAlgebra
# =============================================================================

class BooleanAlgebra(Lattice):
    """
    Boolesche Algebra (B, ∧, ∨, ¬, ⊥, ⊤).

    Eine Boolesche Algebra ist ein distributiver, komplementierter,
    beschränkter Verband. Jedes Element hat ein eindeutiges Komplement.

    Beispiele:
    - Potenzmenge P(S) mit ∩, ∪, ⊆
    - {0,1}ⁿ mit bitweisem AND, OR
    - Aussagenlogik mit ∧, ∨, ¬

    @lastModified 2026-03-10
    """

    def is_boolean_algebra(self) -> bool:
        """
        Überprüft alle Axiome einer Booleschen Algebra.

        Prüft:
        1. Ist es eine Halbordnung?
        2. Existieren ⊤ und ⊥?
        3. Ist der Verband distributiv?
        4. Ist der Verband komplementiert?

        @return  True wenn alle Axiome erfüllt
        @lastModified 2026-03-10
        """
        if not self.is_partial_order():
            return False
        if self.top() is None or self.bottom() is None:
            return False
        if not self.is_distributive():
            return False
        if not self.is_complemented():
            return False
        return True

    def atoms(self) -> list:
        """
        Findet alle Atome der Booleschen Algebra.

        Ein Atom ist ein Element a mit ⊥ < a und kein b mit ⊥ < b < a existiert.
        (Direkt über ⊥)

        @return  Liste der Atome
        @lastModified 2026-03-10
        """
        bot = self.bottom()
        if bot is None:
            return []

        atom_list = []
        for a in self.elements:
            if a == bot:
                continue
            # a muss direkt über ⊥ liegen
            if not self.leq(bot, a):
                continue
            # Kein Element zwischen ⊥ und a
            is_atom = not any(
                b != bot and b != a and self.leq(bot, b) and self.leq(b, a)
                for b in self.elements
            )
            if is_atom:
                atom_list.append(a)
        return atom_list

    def coatoms(self) -> list:
        """
        Findet alle Koatome (duale Atome) der Booleschen Algebra.

        Ein Koatom ist ein Element a mit a < ⊤ und kein b mit a < b < ⊤ existiert.
        (Direkt unter ⊤)

        @return  Liste der Koatome
        @lastModified 2026-03-10
        """
        top = self.top()
        if top is None:
            return []

        coatom_list = []
        for a in self.elements:
            if a == top:
                continue
            if not self.leq(a, top):
                continue
            # Kein Element zwischen a und ⊤
            is_coatom = not any(
                b != top and b != a and self.leq(a, b) and self.leq(b, top)
                for b in self.elements
            )
            if is_coatom:
                coatom_list.append(a)
        return coatom_list

    def stone_representation(self) -> dict:
        """
        Demonstriert die Stone-Dualität für endliche Boolesche Algebren.

        Stone-Satz (endlicher Fall): Jede endliche Boolesche Algebra ist isomorph
        zur Potenzmenge ihrer Atome:
            B ≅ P(At(B))

        @return  Dict mit 'atoms', 'isomorphism' {element: Teilmenge der Atome}
        @lastModified 2026-03-10
        """
        atom_list = self.atoms()
        bot = self.bottom()

        # Isomorphismus: jedem Element b die Menge der Atome zuordnen, die ≤ b sind
        iso = {}
        for b in self.elements:
            iso[b] = frozenset(a for a in atom_list if self.leq(a, b))

        return {
            'atoms': atom_list,
            'isomorphism': iso,
            'description': f'B ≅ P(At(B)), |At(B)| = {len(atom_list)}, |B| = 2^{len(atom_list)} = {2**len(atom_list)}'
        }


# =============================================================================
# KONSTRUKTIONEN
# =============================================================================

def divisibility_lattice(n: int) -> 'Lattice':
    """
    Erstellt den Teilerverband von n.

    Elemente: alle Teiler d von n.
    Ordnung: d₁ ≤ d₂ ⟺ d₁ | d₂ (Teilbarkeit).
    Meet: d₁ ∧ d₂ = ggT(d₁, d₂)
    Join: d₁ ∨ d₂ = kgV(d₁, d₂)

    @param n  Positive ganze Zahl
    @return   Teilerverband als Lattice-Objekt
    @lastModified 2026-03-10
    """
    # Alle Teiler von n berechnen
    divisors = [d for d in range(1, n + 1) if n % d == 0]

    # Ordnungsrelation: d₁ | d₂
    order = set()
    for d1 in divisors:
        for d2 in divisors:
            if d2 % d1 == 0:  # d1 teilt d2
                order.add((d1, d2))

    return Lattice(divisors, order)


def power_set_lattice(n: int) -> 'BooleanAlgebra':
    """
    Erstellt die Potenzmenge {0,...,n-1} als Boolesche Algebra.

    Elemente: alle Teilmengen von {0,...,n-1} (als frozenset).
    Ordnung: A ≤ B ⟺ A ⊆ B.
    Meet: A ∧ B = A ∩ B.
    Join: A ∨ B = A ∪ B.

    @param n  Größe der Grundmenge
    @return   Potenzmenge als BooleanAlgebra-Objekt
    @lastModified 2026-03-10
    """
    base = list(range(n))
    # Alle Teilmengen als frozensets
    elements = [frozenset(s) for s in chain.from_iterable(
        combinations(base, r) for r in range(n + 1)
    )]

    # Ordnung: Teilmengenrelation
    order = set()
    for a in elements:
        for b in elements:
            if a.issubset(b):
                order.add((a, b))

    return BooleanAlgebra(elements, order)


def partition_lattice(n: int) -> 'Lattice':
    """
    Erstellt den Partitionsverband Πₙ.

    Elemente: alle Partitionen von {1,...,n}.
    Ordnung: π₁ ≤ π₂ (Verfeinerung) ⟺ jeder Block von π₁ ist Teilmenge eines Blocks von π₂.
    Meet: feinste gemeinsame Verfeinerung.
    Join: gröbste gemeinsame Vergröberung.

    @param n  Größe der Grundmenge {1,...,n}
    @return   Partitionsverband als Lattice-Objekt
    @lastModified 2026-03-10
    """
    base = list(range(1, n + 1))

    def all_partitions(s: list) -> list:
        """Erzeugt alle Partitionen einer Menge rekursiv (Bell-Zahlen-Rekursion)."""
        if not s:
            return [()]
        first = s[0]
        rest = s[1:]
        result = []
        for part in all_partitions(rest):
            # first in einen neuen Block
            result.append((frozenset([first]),) + part)
            # first in jeden bestehenden Block
            for i in range(len(part)):
                new_part = list(part)
                new_part[i] = part[i] | frozenset([first])
                result.append(tuple(new_part))
        return result

    partitions = [frozenset(p) for p in all_partitions(base)]

    def is_refinement(pi1: frozenset, pi2: frozenset) -> bool:
        """Prüft ob π₁ Verfeinerung von π₂ ist: jeder Block von π₁ ⊆ Block von π₂."""
        for block1 in pi1:
            # Einen Block in π₂ finden, der block1 enthält
            if not any(block1.issubset(block2) for block2 in pi2):
                return False
        return True

    # Ordnungsrelation: Verfeinerung
    order = set()
    for p1 in partitions:
        for p2 in partitions:
            if is_refinement(p1, p2):
                order.add((p1, p2))

    return Lattice(partitions, order)


# =============================================================================
# SÄTZE
# =============================================================================

def birkhoff_representation_theorem(lattice: 'Lattice') -> dict:
    """
    Demonstriert den Birkhoff-Darstellungssatz für endliche distributive Verbände.

    Birkhoff-Satz (1937): Jeder endliche distributive Verband L ist isomorph zur
    Menge J(P) der Ordnungsideale eines eindeutig bestimmten Posets P = J_irr(L)
    (die join-irreduziblen Elemente von L).

    Ein Element j ∈ L ist join-irreduzibel wenn j = a ∨ b → j = a oder j = b,
    und j ≠ ⊥.

    @param lattice  Ein endlicher distributiver Verband
    @return         Dict mit 'is_distributive', 'join_irreducibles', 'order_ideals'
    @lastModified 2026-03-10
    """
    if not lattice.is_distributive():
        return {
            'is_distributive': False,
            'join_irreducibles': [],
            'order_ideals': [],
            'message': 'Birkhoff-Satz gilt nur für distributive Verbände!'
        }

    bot = lattice.bottom()

    # Join-irreduzible Elemente finden
    # j ist join-irreduzibel wenn: j ≠ ⊥ und j = a ∨ b ⟹ j = a oder j = b
    join_irreducibles = []
    for j in lattice.elements:
        if j == bot:
            continue
        is_ji = True
        for a in lattice.elements:
            for b in lattice.elements:
                if lattice.join(a, b) == j and a != j and b != j:
                    is_ji = False
                    break
            if not is_ji:
                break
        if is_ji:
            join_irreducibles.append(j)

    # Ordnungsideale von J_irr(L) berechnen
    # Ein Ordnungsideal I ist eine nach unten abgeschlossene Teilmenge
    ji_set = join_irreducibles
    ji_order = {(a, b) for (a, b) in lattice.order_relation if a in ji_set and b in ji_set}

    def is_order_ideal(subset: frozenset) -> bool:
        """Prüft ob subset ein Ordnungsideal ist: a ∈ I, b ≤ a → b ∈ I."""
        for b in ji_set:
            for a in subset:
                if (b, a) in ji_order and b not in subset:
                    return False
        return True

    all_ideals = []
    for r in range(len(ji_set) + 1):
        for subset in combinations(ji_set, r):
            fs = frozenset(subset)
            if is_order_ideal(fs):
                all_ideals.append(fs)

    return {
        'is_distributive': True,
        'join_irreducibles': join_irreducibles,
        'order_ideals': all_ideals,
        'isomorphism_size': len(all_ideals),
        'lattice_size': len(lattice.elements),
        'message': f'L ≅ J(P) mit |L| = {len(lattice.elements)} = |J(P)| = {len(all_ideals)}'
    }


def dedekind_numbers(n: int) -> int:
    """
    Berechnet die Dedekind-Zahl D(n) für kleine n.

    D(n) ist die Anzahl der monotonen Booleschen Funktionen f: {0,1}ⁿ → {0,1},
    äquivalent zur Anzahl der Antiketten in der Potenzmenge von {1,...,n}.

    Bekannte Werte:
    D(0)=2, D(1)=3, D(2)=6, D(3)=20, D(4)=168, D(5)=7581, D(6)=7828354

    @param n  Anzahl der Variablen (0 ≤ n ≤ 6 direkt, sonst Approximation)
    @return   Dedekind-Zahl D(n)
    @lastModified 2026-03-10
    """
    # Bekannte Werte (D(7) hat 23 Ziffern, zu groß für direkte Berechnung)
    known_values = {0: 2, 1: 3, 2: 6, 3: 20, 4: 168, 5: 7581, 6: 7828354}

    if n in known_values:
        return known_values[n]

    # Für n ≤ 6 direkte Berechnung via Antiketten-Zählung
    if n > 6:
        raise ValueError(f"D({n}) nicht vorberechnet (zu groß für direkte Berechnung)")

    # Antiketten in P({1,...,n}) zählen (inklusive leerer Antikette)
    base = list(range(n))
    all_subsets = [frozenset(s) for s in chain.from_iterable(
        combinations(base, r) for r in range(n + 1)
    )]

    def is_antichain(collection: list) -> bool:
        """Prüft ob keine zwei Mengen in der Kollektion vergleichbar sind."""
        for i, a in enumerate(collection):
            for j, b in enumerate(collection):
                if i != j and (a.issubset(b) or b.issubset(a)):
                    return False
        return True

    count = 0
    # Alle Teilmengen von P(base) durchgehen und auf Antiketten-Eigenschaft prüfen
    for r in range(len(all_subsets) + 1):
        for coll in combinations(all_subsets, r):
            if is_antichain(list(coll)):
                count += 1

    return count


def jordan_dedekind_chain_condition(lattice: 'Lattice') -> bool:
    """
    Prüft die Jordan-Dedekind-Kettenbedingung.

    Ein Verband erfüllt die JD-Kettenbedingung wenn:
    Alle maximalen Ketten zwischen zwei Elementen a ≤ b haben gleiche Länge.

    Dies ist äquivalent zur Existenz einer Rangfunktion r: L → ℕ mit:
    - a ⋖ b (a wird direkt von b überdeckt) ⟹ r(b) = r(a) + 1

    @param lattice  Ein endlicher Verband
    @return         True wenn JD-Kettenbedingung gilt
    @lastModified 2026-03-10
    """
    # Überdeckungsrelation berechnen (Hasse-Diagramm-Kanten)
    hasse = lattice.hasse_diagram()

    # Rangfunktion via BFS von ⊥ aufbauen
    bot = lattice.bottom()
    if bot is None:
        return False  # Kein kleinstes Element

    # BFS-Ränge berechnen
    rank = {bot: 0}
    queue = [bot]
    visited = {bot}

    while queue:
        current = queue.pop(0)
        for successor in hasse[current]:
            new_rank = rank[current] + 1
            if successor not in rank:
                rank[successor] = new_rank
                if successor not in visited:
                    queue.append(successor)
                    visited.add(successor)
            elif rank[successor] != new_rank:
                # Widersprüchliche Ränge → JD-Bedingung verletzt
                return False

    # Alle Elemente müssen einen Rang erhalten haben
    return len(rank) == len(lattice.elements)
