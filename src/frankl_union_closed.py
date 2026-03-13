"""
@file frankl_union_closed.py
@brief Frankl Union-Closed Conjecture — Computational Experiments
@description
    Frankl's Conjecture (1979): For any union-closed family F of finite sets
    (not consisting only of the empty set), there exists an element that
    belongs to at least |F|/2 of the sets.

    Status: OPEN (CONJECTURE). Best known bounds:
    - Gilmer 2022: existiert Element in ≥ 0.01 * |F| Mengen
    - Chase-Lovász 2023: existiert Element in ≥ 0.38 * |F| Mengen

    Computational approach:
    - Verify for all union-closed families on small ground sets (|U| ≤ 5)
    - Compute frequency vectors and find most frequent elements
    - Exhaustive search via bitset representation

@author Michael Fuhrmann
@date 2026-03-13
@lastModified 2026-03-13
"""

import itertools
import random
from typing import FrozenSet, List, Optional, Tuple
from collections import Counter


class UnionClosedFamily:
    """
    Repräsentiert eine vereinigungsabgeschlossene Familie von Mengen.

    Eine Familie F ist vereinigungsabgeschlossen (union-closed), wenn für alle
    A, B ∈ F gilt: A ∪ B ∈ F.

    Frankl's CONJECTURE (1979): Falls F ≠ {∅}, existiert ein Element das
    in mindestens |F|/2 Mengen vorkommt.
    """

    def __init__(self, sets: List[FrozenSet[int]]):
        """
        Initialisiert die Familie mit einer Liste von frozensets.

        @param sets          Liste von frozensets (Mengen der Familie)
        @lastModified 2026-03-13
        """
        # Duplikate entfernen, Reihenfolge beibehalten (als Liste)
        seen = set()
        unique = []
        for s in sets:
            fs = frozenset(s)
            if fs not in seen:
                seen.add(fs)
                unique.append(fs)
        self.sets: List[FrozenSet[int]] = unique
        # Grundmenge: alle vorkommenden Elemente
        self.universe: FrozenSet[int] = frozenset(
            e for s in self.sets for e in s
        )

    @classmethod
    def from_bitmasks(cls, n: int, bitmasks: List[int]) -> "UnionClosedFamily":
        """
        Erstellt eine Familie aus Bitmasken über der Grundmenge {0, ..., n-1}.

        Bit i der Maske entspricht dem Element i.

        @param n         Größe der Grundmenge
        @param bitmasks  Liste ganzer Zahlen als Bitmasken
        @return          UnionClosedFamily-Instanz
        @lastModified 2026-03-13
        """
        sets = []
        for mask in bitmasks:
            # Jedes gesetzte Bit i entspricht dem Element i in der Menge
            s = frozenset(i for i in range(n) if mask & (1 << i))
            sets.append(s)
        return cls(sets)

    def is_union_closed(self) -> bool:
        """
        Prüft, ob die Familie vereinigungsabgeschlossen ist.

        Für alle A, B ∈ F muss A ∪ B ∈ F gelten.

        @return True falls vereinigungsabgeschlossen, sonst False
        @lastModified 2026-03-13
        """
        set_collection = frozenset(self.sets)
        for a in self.sets:
            for b in self.sets:
                union = a | b
                if union not in set_collection:
                    return False
        return True

    def frequency_vector(self) -> List[int]:
        """
        Berechnet den Häufigkeitsvektor der Familie.

        freq[i] = Anzahl der Mengen in F, die das Element i enthalten.
        Die Reihenfolge entspricht der sortierten Grundmenge.

        @return Liste mit Häufigkeiten für jedes Element der Grundmenge
        @lastModified 2026-03-13
        """
        if not self.universe:
            return []
        sorted_elements = sorted(self.universe)
        freq = []
        for elem in sorted_elements:
            count = sum(1 for s in self.sets if elem in s)
            freq.append(count)
        return freq

    def max_frequency_ratio(self) -> float:
        """
        Berechnet das maximale Verhältnis max(freq[i]) / |F|.

        Frankl's CONJECTURE besagt: dieses Verhältnis ist ≥ 1/2
        für alle nicht-trivialen Familien.

        @return maximales Häufigkeitsverhältnis (0.0 falls Familie leer)
        @lastModified 2026-03-13
        """
        if not self.sets:
            return 0.0
        freq = self.frequency_vector()
        if not freq:
            # Nur die leere Menge ist enthalten
            return 0.0
        return max(freq) / len(self.sets)

    def satisfies_frankl(self) -> bool:
        """
        Prüft, ob Frankl's CONJECTURE für diese Familie erfüllt ist.

        Bedingung: max(freq[i]) ≥ |F| / 2, d.h. max_frequency_ratio ≥ 0.5.

        Hinweis: Familien, die nur die leere Menge enthalten, sind ausgenommen.

        @return True falls Frankl erfüllt
        @lastModified 2026-03-13
        """
        # Nur-leere-Menge-Familie ist ausgenommen (Bedingung F ≠ {∅})
        non_empty = [s for s in self.sets if len(s) > 0]
        if not non_empty:
            return True  # Trivialfall: keine nichtleere Menge
        return self.max_frequency_ratio() >= 0.5

    def frankl_witness(self) -> Optional[int]:
        """
        Gibt ein Element zurück, das in ≥ |F|/2 Mengen vorkommt (Zeuge).

        @return Zeugen-Element oder None falls kein Zeuge existiert
        @lastModified 2026-03-13
        """
        if not self.universe or not self.sets:
            return None
        sorted_elements = sorted(self.universe)
        threshold = len(self.sets) / 2.0
        for elem in sorted_elements:
            count = sum(1 for s in self.sets if elem in s)
            if count >= threshold:
                return elem
        return None

    def __len__(self) -> int:
        """Gibt die Anzahl der Mengen in der Familie zurück."""
        return len(self.sets)

    def __repr__(self) -> str:
        """String-Darstellung der Familie."""
        return f"UnionClosedFamily({self.sets})"


class FranklChecker:
    """
    Verifiziert Frankl's CONJECTURE exhaustiv für kleine Familien.

    Für Grundmengen der Größe n ≤ 5 können alle vereinigungsabgeschlossenen
    Familien enumeriert und geprüft werden.
    """

    @staticmethod
    def all_union_closed_families(n: int) -> List[UnionClosedFamily]:
        """
        Enumeriert alle vereinigungsabgeschlossenen Familien über {0,...,n-1}.

        Methode: Iteriere über alle 2^(2^n) Teilmengen der Potenzmenge und
        behalte nur die vereinigungsabgeschlossenen.
        Praktisch nur für n ≤ 4 durchführbar.

        @param n  Grundmengengröße (empfohlen: ≤ 4)
        @return   Liste aller uc-Familien über {0,...,n-1}
        @lastModified 2026-03-13
        """
        # Potenzmenge als Liste von frozensets
        ground = list(range(n))
        all_subsets = [
            frozenset(combo)
            for r in range(n + 1)
            for combo in itertools.combinations(ground, r)
        ]
        num_subsets = len(all_subsets)  # = 2^n

        result = []
        # Iteriere über alle 2^(2^n) Teilfamilien der Potenzmenge
        for mask in range(1, 1 << num_subsets):
            # Baue Familie aus den durch mask ausgewählten Teilmengen
            family_sets = [
                all_subsets[i]
                for i in range(num_subsets)
                if mask & (1 << i)
            ]
            fam = UnionClosedFamily(family_sets)
            if fam.is_union_closed():
                result.append(fam)
        return result

    @staticmethod
    def verify_up_to_n(n: int) -> Tuple[int, int, bool]:
        """
        Prüft Frankl's CONJECTURE für alle uc-Familien über Grundmengen bis Größe n.

        Familien die nur aus der leeren Menge bestehen sind ausgenommen.

        @param n  Maximale Grundmengengröße
        @return   Tupel (total_families, violations, all_ok)
        @lastModified 2026-03-13
        """
        total = 0
        violations = 0
        for size in range(1, n + 1):
            families = FranklChecker.all_union_closed_families(size)
            for fam in families:
                # Nur nicht-triviale Familien prüfen
                non_empty_sets = [s for s in fam.sets if len(s) > 0]
                if not non_empty_sets:
                    continue
                total += 1
                if not fam.satisfies_frankl():
                    violations += 1
        all_ok = violations == 0
        return (total, violations, all_ok)

    @staticmethod
    def gilmer_bound(family_size: int) -> float:
        """
        Gilmer 2022: Existiert Element in ≥ 0.01 * |F| Mengen.

        @param family_size  Größe der Familie |F|
        @return             untere Schranke nach Gilmer
        @lastModified 2026-03-13
        """
        return 0.01 * family_size

    @staticmethod
    def chase_lovasz_bound(family_size: int) -> float:
        """
        Chase-Lovász 2023: Existiert Element in ≥ 0.38 * |F| Mengen.

        @param family_size  Größe der Familie |F|
        @return             untere Schranke nach Chase-Lovász
        @lastModified 2026-03-13
        """
        return 0.38 * family_size


class UnionClosedGenerator:
    """
    Erzeugt spezielle Klassen von vereinigungsabgeschlossenen Familien.

    Diese dienen als Testfälle und Beispiele für Frankl's CONJECTURE.
    """

    @staticmethod
    def power_set(n: int) -> UnionClosedFamily:
        """
        Erzeugt die Potenzmenge 2^{0,...,n-1}.

        Die Potenzmenge ist trivialerweise vereinigungsabgeschlossen.
        Frankl ist erfüllt, da jedes Element in genau 2^(n-1) von 2^n Mengen
        vorkommt, also in der Hälfte.

        @param n  Grundmengengröße
        @return   Potenzmenge als UnionClosedFamily
        @lastModified 2026-03-13
        """
        ground = list(range(n))
        all_subsets = [
            frozenset(combo)
            for r in range(n + 1)
            for combo in itertools.combinations(ground, r)
        ]
        return UnionClosedFamily(all_subsets)

    @staticmethod
    def interval_family(n: int) -> UnionClosedFamily:
        """
        Erzeugt die Intervall-Familie: alle Mengen {i, i+1, ..., j} für 0 ≤ i ≤ j < n.

        Hinweis: Die Familie der ganzzahligen Intervalle ist NICHT zwingend
        vereinigungsabgeschlossen für n ≥ 3, da {0} ∪ {2} = {0,2} kein Intervall ist.
        Diese Methode gibt die Intervall-Mengen zurück (ohne uc-Abschluss).
        Für einen echten uc-Abschluss: UnionClosedGenerator.union_closure(fam.sets).

        @param n  Grundmengengröße
        @return   Familie der Intervall-Mengen als UnionClosedFamily
        @lastModified 2026-03-13
        """
        sets = []
        for i in range(n):
            for j in range(i, n):
                sets.append(frozenset(range(i, j + 1)))
        return UnionClosedFamily(sets)

    @staticmethod
    def random_union_closed(n: int, seed: int = 42) -> UnionClosedFamily:
        """
        Erzeugt eine zufällige vereinigungsabgeschlossene Familie.

        Methode: Wähle eine Zufallsfamilie und bilde deren Vereinigungsabschluss.

        @param n     Grundmengengröße
        @param seed  Zufallsseed für Reproduzierbarkeit
        @return      Zufällige uc-Familie
        @lastModified 2026-03-13
        """
        rng = random.Random(seed)
        ground = list(range(n))
        # Zufällige Auswahl von Teilmengen als Startfamilie
        all_subsets = [
            frozenset(combo)
            for r in range(1, n + 1)
            for combo in itertools.combinations(ground, r)
        ]
        k = max(1, rng.randint(1, max(1, len(all_subsets) // 2)))
        start_sets = rng.sample(all_subsets, min(k, len(all_subsets)))
        return UnionClosedGenerator.union_closure(start_sets)

    @staticmethod
    def union_closure(sets: List[FrozenSet[int]]) -> "UnionClosedFamily":
        """
        Berechnet den Vereinigungsabschluss einer gegebenen Familie.

        Iteriert, bis keine neuen Mengen mehr durch Vereinigung entstehen.

        @param sets  Ausgangsfamilie (Liste von frozensets)
        @return      Vereinigungsabgeschlossene Hülle
        @lastModified 2026-03-13
        """
        closed: set = set(frozenset(s) for s in sets)
        changed = True
        while changed:
            changed = False
            current = list(closed)
            for a in current:
                for b in current:
                    union = a | b
                    if union not in closed:
                        closed.add(union)
                        changed = True
        return UnionClosedFamily(list(closed))

    @staticmethod
    def chain_family(n: int) -> UnionClosedFamily:
        """
        Erzeugt eine Ketten-Familie: {{0}, {0,1}, {0,1,2}, ..., {0,...,n-1}}.

        Ketten sind vereinigungsabgeschlossen.

        @param n  Länge der Kette
        @return   Ketten-Familie als UnionClosedFamily
        @lastModified 2026-03-13
        """
        sets = [frozenset(range(i + 1)) for i in range(n)]
        return UnionClosedFamily(sets)
