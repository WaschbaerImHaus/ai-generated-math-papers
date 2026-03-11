"""
@file algebraic_structures.py
@brief Algebraische Strukturhierarchie: Magma → Halbgruppe → Monoid → Gruppe.
       Freie algebraische Strukturen, Cayley-Tabellen, Transformationsmonoide.
@description
    Dieses Modul implementiert die grundlegende Hierarchie algebraischer Strukturen:

    Klassen:
    - Magma              – Menge mit binärer Operation (keine Axiome)
    - Semigroup          – Magma mit Assoziativität
    - Monoid             – Halbgruppe mit neutralem Element
    - GroupFromMagma     – Monoid mit Inversen

    Freie Strukturen:
    - free_monoid()              – Freier Monoid über Alphabet (Wörter mit Konkatenation)
    - free_semigroup()           – Freie Halbgruppe (nicht-leere Wörter)
    - algebraic_structure_hierarchy() – Hierarchieübersicht als Dict

    Beispiele:
    - string_monoid()            – Monoid der Strings über endlichem Alphabet
    - matrix_monoid()            – n×n-Matrizen über ℤ/2ℤ unter Multiplikation
    - transformation_monoid()    – Alle Funktionen {1,...,n}→{1,...,n}

    Mathematische Grundlagen:
    - Magma: (M, ∗) mit a∗b ∈ M (Abgeschlossenheit)
    - Halbgruppe: Magma + (a∗b)∗c = a∗(b∗c)
    - Monoid: Halbgruppe + ∃e: e∗a = a∗e = a
    - Gruppe: Monoid + ∀a ∃a⁻¹: a∗a⁻¹ = e

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

from typing import Callable, Any, Optional
import itertools
import sys
import os

# Eigene Module aus demselben Verzeichnis importieren
sys.path.insert(0, os.path.dirname(__file__))


class Magma:
    """
    Magma (M, ∗): Menge mit einer inneren binären Operation.

    Das ist die allgemeinste algebraische Struktur — es werden keine Axiome
    (Assoziativität, Neutralelement, Inverse) gefordert. Lediglich die
    Abgeschlossenheit (∀a,b∈M: a∗b∈M) gilt per Definition.

    Beispiele: Nicht-assoziative Algebren, Quasigruppen.
    """

    def __init__(self, elements: list, operation: Callable[[Any, Any], Any]):
        """
        Initialisiert ein Magma mit einer Elementliste und binärer Operation.

        @param elements  Liste der Elemente der Menge M
        @param operation Callable: (a, b) → a∗b (binäre Operation)
        """
        # Elementliste als Liste speichern (Reihenfolge für Cayley-Tabelle)
        self.elements = list(elements)
        self.operation = operation
        # Index-Mapping für schnellen Lookup
        self._index = {e: i for i, e in enumerate(self.elements)}

    def op(self, a: Any, b: Any) -> Any:
        """
        Verknüpft zwei Elemente mit der binären Operation.

        @param a Erstes Element
        @param b Zweites Element
        @return  Ergebnis a∗b
        """
        return self.operation(a, b)

    def operation_table(self) -> list:
        """
        Erstellt die Cayley-Tabelle (Verknüpfungstafel) des Magmas.

        Die Cayley-Tabelle zeigt für alle Paare (a, b) das Ergebnis a∗b.
        Zeile i, Spalte j enthält elements[i] ∗ elements[j].

        @return 2D-Liste: table[i][j] = elements[i] ∗ elements[j]
        """
        n = len(self.elements)
        # Zweidimensionale Liste mit allen Verknüpfungen aufbauen
        table = []
        for a in self.elements:
            row = []
            for b in self.elements:
                row.append(self.op(a, b))
            table.append(row)
        return table

    def is_closed(self) -> bool:
        """
        Prüft ob die Operation abgeschlossen ist: ∀a,b∈M: a∗b∈M.

        Bei endlichen Mengen mit expliziter Elementliste ist dies genau dann
        erfüllt, wenn alle Ergebnisse der Cayley-Tabelle in der Elementmenge liegen.

        @return True wenn abgeschlossen, False sonst
        """
        elem_set = set(self.elements)
        for a in self.elements:
            for b in self.elements:
                result = self.op(a, b)
                if result not in elem_set:
                    return False
        return True

    def is_commutative(self) -> bool:
        """
        Prüft ob die Operation kommutativ ist: ∀a,b∈M: a∗b = b∗a.

        Ein kommutatives Magma wird auch 'kommutativer Quasigruppenoid' genannt.
        Für Gruppen ergibt dies eine abelsche Gruppe.

        @return True wenn kommutativ, False sonst
        """
        for a in self.elements:
            for b in self.elements:
                if self.op(a, b) != self.op(b, a):
                    return False
        return True


class Semigroup(Magma):
    """
    Halbgruppe (S, ∗): Magma mit Assoziativität.

    Zusätzliches Axiom gegenüber Magma:
    Assoziativität: ∀a,b,c∈S: (a∗b)∗c = a∗(b∗c)

    Beispiele: (ℕ, +), (ℕ, ·), (Σ⁺, Konkatenation), Matrizen unter Multiplikation.
    """

    def is_associative(self) -> bool:
        """
        Prüft das Assoziativgesetz: ∀a,b,c: (a∗b)∗c = a∗(b∗c).

        Laufzeit: O(n³) für n Elemente — bei großen Strukturen langsam.

        @return True wenn assoziativ, False sobald Gegenbeispiel gefunden
        """
        for a in self.elements:
            for b in self.elements:
                for c in self.elements:
                    # Linke Seite: (a∗b)∗c
                    left = self.op(self.op(a, b), c)
                    # Rechte Seite: a∗(b∗c)
                    right = self.op(a, self.op(b, c))
                    if left != right:
                        return False
        return True

    def is_semigroup(self) -> bool:
        """
        Prüft ob diese Struktur eine Halbgruppe ist.

        Eine Halbgruppe ist ein abgeschlossenes, assoziatives Magma.

        @return True wenn alle Halbgruppen-Axiome erfüllt sind
        """
        return self.is_closed() and self.is_associative()


class Monoid(Semigroup):
    """
    Monoid (M, ∗, e): Halbgruppe mit neutralem Element e.

    Zusätzliches Axiom gegenüber Halbgruppe:
    Neutrales Element: ∃e∈M: ∀a∈M: e∗a = a∗e = a

    Beispiele: (ℕ₀, +, 0), (ℕ, ·, 1), (Σ*, Konkatenation, ε), (Matrizen, ·, I).
    """

    def identity_element(self) -> Any:
        """
        Sucht und gibt das neutrale Element e zurück.

        Das neutrale Element e erfüllt: e∗a = a∗e = a für alle a∈M.
        Es ist eindeutig, falls es existiert.

        @return Neutrales Element oder None falls keines existiert
        """
        for e in self.elements:
            # Prüfe ob e neutral für alle anderen Elemente ist
            is_identity = True
            for a in self.elements:
                if self.op(e, a) != a or self.op(a, e) != a:
                    is_identity = False
                    break
            if is_identity:
                return e
        return None

    def has_identity(self) -> bool:
        """
        Prüft ob ein neutrales Element existiert.

        @return True wenn neutrales Element vorhanden, False sonst
        """
        return self.identity_element() is not None

    def is_monoid(self) -> bool:
        """
        Prüft ob diese Struktur ein Monoid ist.

        Ein Monoid ist eine Halbgruppe mit neutralem Element.

        @return True wenn alle Monoid-Axiome erfüllt sind
        """
        return self.is_semigroup() and self.has_identity()

    def powers(self, element: Any, n: int) -> list:
        """
        Berechnet die Potenzen [e, a, a², ..., aⁿ] eines Elements.

        Potenzen werden iterativ durch Linksmultiplikation berechnet:
        a⁰ = e, a^k = a^{k-1} ∗ a

        @param element Element a ∈ M
        @param n       Maximale Potenz
        @return        Liste [e, a, a², ..., aⁿ]
        @raises        ValueError wenn kein neutrales Element vorhanden
        """
        e = self.identity_element()
        if e is None:
            raise ValueError("Kein neutrales Element — powers() erfordert Monoid")
        # Potenzen iterativ aufbauen
        result = [e]
        current = e
        for _ in range(n):
            current = self.op(current, element)
            result.append(current)
        return result


class GroupFromMagma(Monoid):
    """
    Gruppe (G, ∗, e): Monoid, bei dem jedes Element ein Inverses besitzt.

    Zusätzliches Axiom gegenüber Monoid:
    Inverse: ∀a∈G ∃a⁻¹∈G: a∗a⁻¹ = a⁻¹∗a = e

    Diese Klasse erweitert Monoid und benutzt die Magma-Infrastruktur,
    um Gruppen aus expliziten Elementlisten und Operationen zu erzeugen.
    """

    def inverse(self, element: Any) -> Any:
        """
        Findet das inverse Element a⁻¹ für ein gegebenes a.

        Durchsucht alle Elemente nach b mit a∗b = b∗a = e.

        @param element Element a ∈ G
        @return        Inverses a⁻¹ oder None falls nicht vorhanden
        @raises        ValueError wenn kein neutrales Element existiert
        """
        e = self.identity_element()
        if e is None:
            raise ValueError("Kein neutrales Element — inverse() erfordert Monoid")
        # Suche b mit a∗b = e und b∗a = e
        for b in self.elements:
            if self.op(element, b) == e and self.op(b, element) == e:
                return b
        return None

    def is_group(self) -> bool:
        """
        Prüft ob diese Struktur eine Gruppe ist.

        Eine Gruppe ist ein Monoid, in dem jedes Element ein Inverses hat.

        @return True wenn alle Gruppen-Axiome erfüllt sind
        """
        if not self.is_monoid():
            return False
        # Prüfe ob jedes Element ein Inverses besitzt
        for a in self.elements:
            if self.inverse(a) is None:
                return False
        return True


# ---------------------------------------------------------------------------
# Freie algebraische Strukturen
# ---------------------------------------------------------------------------

def free_monoid(alphabet: list) -> dict:
    """
    Freier Monoid über einem Alphabet Σ: Menge aller endlichen Wörter mit Konkatenation.

    Der freie Monoid (Σ*, ·, ε) ist fundamental in der theoretischen Informatik
    und formalen Sprachtheorie. Er ist charakterisiert durch die universelle
    Eigenschaft: jede Funktion f: Σ → M in einen Monoid M faktorisiert eindeutig
    über den freien Monoid.

    @param alphabet  Liste von Symbolen (die Generatoren)
    @return          Dict mit 'alphabet', 'identity', 'operation', 'description'
    """
    return {
        'alphabet': list(alphabet),
        # Das neutrale Element ist das leere Wort ε
        'identity': (),
        # Die Operation ist die Konkatenation von Tupeln
        'operation': lambda a, b: a + b,
        'description': (
            f"Freier Monoid Σ* über Σ={alphabet}. "
            "Elemente: alle endlichen Wörter (Tupel) über Σ. "
            "Operation: Konkatenation. Neutrales Element: ε = ()."
        ),
        # Beispielwörter bis Länge 2
        'example_words': (
            [()]
            + [(x,) for x in alphabet]
            + [(x, y) for x in alphabet for y in alphabet]
        ),
    }


def free_semigroup(generators: list, max_length: int = 3) -> dict:
    """
    Freie Halbgruppe über Generatoren: alle nicht-leeren Wörter bis max_length.

    Im Gegensatz zum freien Monoid enthält die freie Halbgruppe kein leeres Wort.
    (Σ⁺, ·) ist die Restriktion von Σ* auf nicht-leere Wörter.

    @param generators  Liste von Generatoren/Symbolen
    @param max_length  Maximale Wortlänge für Aufzählung (default: 3)
    @return            Dict mit Generatoren, Wörtern und Beschreibung
    """
    # Alle nicht-leeren Wörter bis max_length aufzählen
    words = []
    for length in range(1, max_length + 1):
        for word in itertools.product(generators, repeat=length):
            words.append(word)

    return {
        'generators': list(generators),
        'max_length': max_length,
        'words': words,
        'count': len(words),
        'operation': lambda a, b: a + b,
        'description': (
            f"Freie Halbgruppe Σ⁺ über Generatoren {generators} "
            f"(Wörter bis Länge {max_length}). "
            "Operation: Konkatenation. Kein neutrales Element."
        ),
    }


def algebraic_structure_hierarchy() -> dict:
    """
    Gibt die Hierarchie algebraischer Strukturen als Dictionary zurück.

    Die Hierarchie folgt dem Schema:
    Magma ⊃ Halbgruppe ⊃ Monoid ⊃ Gruppe ⊃ Abelsche Gruppe
    und verzweigt sich bei Ring → Integritätsbereich → Körper.

    @return Dict mit Strukturname als Key, Axiomen und Beziehungen als Value
    """
    return {
        'Magma': {
            'axiome': ['Abgeschlossenheit: ∀a,b∈M: a∗b∈M'],
            'examples': ['Oktaven (nicht-assoziativ)', 'Quasigruppen'],
            'erweitert_zu': ['Halbgruppe'],
        },
        'Halbgruppe': {
            'axiome': [
                'Abgeschlossenheit',
                'Assoziativität: (a∗b)∗c = a∗(b∗c)',
            ],
            'examples': ['(ℕ,+)', '(ℕ,·)', 'Σ⁺ (nichtleere Wörter)'],
            'erweitert_zu': ['Monoid'],
        },
        'Monoid': {
            'axiome': [
                'Abgeschlossenheit',
                'Assoziativität',
                'Neutrales Element: ∃e: e∗a=a∗e=a',
            ],
            'examples': ['(ℕ₀,+,0)', '(ℕ,·,1)', 'Σ* (alle Wörter)', 'Matrizen (·,I)'],
            'erweitert_zu': ['Gruppe'],
        },
        'Gruppe': {
            'axiome': [
                'Abgeschlossenheit',
                'Assoziativität',
                'Neutrales Element',
                'Inverse: ∀a ∃a⁻¹: a∗a⁻¹=e',
            ],
            'examples': ['(ℤ,+,0)', 'GL(n,ℝ)', 'Sₙ', 'Diedergruppe Dₙ'],
            'erweitert_zu': ['Abelsche Gruppe', 'Ring (additive Gruppe)'],
        },
        'Abelsche Gruppe': {
            'axiome': [
                'Alle Gruppenaxiome',
                'Kommutativität: a∗b=b∗a',
            ],
            'examples': ['(ℤ,+)', '(ℝ,+)', '(ℤ/nℤ,+)'],
            'erweitert_zu': ['Ring'],
        },
        'Ring': {
            'axiome': [
                '(R,+) abelsche Gruppe',
                '(R,·) Monoid',
                'Distributivgesetze: a·(b+c)=a·b+a·c',
            ],
            'examples': ['ℤ', 'ℤ[x]', 'Matrizenringe', 'ℤ/nℤ'],
            'erweitert_zu': ['Integritätsbereich'],
        },
        'Integritätsbereich': {
            'axiome': [
                'Kommutativer Ring',
                'Nullteilerfreiheit: a·b=0 ⟹ a=0 oder b=0',
            ],
            'examples': ['ℤ', 'ℤ[x]', 'Gaussian integers ℤ[i]'],
            'erweitert_zu': ['Körper'],
        },
        'Körper': {
            'axiome': [
                'Kommutativer Ring',
                'Jedes Nichtnull-Element hat multiplikatives Inverses',
            ],
            'examples': ['ℚ', 'ℝ', 'ℂ', 'ℤ/pℤ (p prim)', 'GF(2^n)'],
            'erweitert_zu': [],
        },
    }


# ---------------------------------------------------------------------------
# Konkrete Monoid-Beispiele
# ---------------------------------------------------------------------------

def string_monoid(alphabet: str = "ab", max_length: int = 2) -> Monoid:
    """
    Monoid der Strings über einem endlichen Alphabet mit Konkatenation.

    Elemente: alle Strings über 'alphabet' bis Länge max_length (inkl. Leerstring).
    Operation: Stringkonkatenation.
    Neutrales Element: Leerstring ''.

    Da der echte freie Monoid unendlich ist, wird hier auf max_length beschränkt.
    Die Abgeschlossenheit gilt nur, wenn max_length ≥ 2·(Länge der längsten Wörter).

    @param alphabet   Zeichenkette der erlaubten Zeichen (default: "ab")
    @param max_length Maximale Stringlänge (default: 2)
    @return           Monoid-Instanz
    """
    # Alle Strings bis max_length erzeugen
    elements = ['']  # Leerstring = neutrales Element
    for length in range(1, max_length + 1):
        for chars in itertools.product(alphabet, repeat=length):
            elements.append(''.join(chars))

    # Abgeschnittene Konkatenation (für Abgeschlossenheit auf max_length)
    def concat(a: str, b: str) -> str:
        """Konkatenation, abgeschnitten auf max_length."""
        result = a + b
        # Ergebnis auf max_length kürzen damit Abgeschlossenheit gewährleistet
        return result[:max_length] if len(result) > max_length else result

    return Monoid(elements, concat)


def matrix_monoid(n: int, field: str = 'Z2') -> Monoid:
    """
    Monoid der n×n-Matrizen über ℤ/2ℤ (GF(2)) unter Multiplikation.

    Elemente: alle n×n-Matrizen mit Einträgen in {0,1}.
    Operation: Matrixmultiplikation modulo 2.
    Neutrales Element: Einheitsmatrix.

    Achtung: Anzahl der Elemente = 2^(n²) — wächst schnell!
    Für n≥3 wird die Elementmenge sehr groß.

    @param n     Matrizengröße n×n
    @param field Körper (nur 'Z2' unterstützt)
    @return      Monoid-Instanz
    """
    if field != 'Z2':
        raise NotImplementedError("Nur ℤ/2ℤ (field='Z2') implementiert")
    if n > 3:
        raise ValueError("n > 3 nicht sinnvoll (2^(n²) Matrizen werden zu viele)")

    # Alle n×n-Matrizen über {0,1} erzeugen
    # Jede Matrix als Tupel von Tupeln (hashbar)
    size = n * n
    elements = []
    for bits in itertools.product([0, 1], repeat=size):
        matrix = tuple(
            tuple(bits[i * n + j] for j in range(n))
            for i in range(n)
        )
        elements.append(matrix)

    def mat_mult(A: tuple, B: tuple) -> tuple:
        """Matrixmultiplikation modulo 2."""
        result = []
        for i in range(n):
            row = []
            for j in range(n):
                # Skalarprodukt der i-ten Zeile von A mit j-ter Spalte von B, mod 2
                val = sum(A[i][k] * B[k][j] for k in range(n)) % 2
                row.append(val)
            result.append(tuple(row))
        return tuple(result)

    return Monoid(elements, mat_mult)


def transformation_monoid(n: int) -> Monoid:
    """
    Transformationsmonoid T_n: alle Funktionen {1,...,n} → {1,...,n}.

    Elemente: alle nⁿ Abbildungen von {1,...,n} in sich selbst (als Tupel).
    Operation: Komposition f∘g (f nach g).
    Neutrales Element: Identitätsabbildung (1, 2, ..., n).

    T_n ist ein wichtiger Monoid mit n^n Elementen.
    Für n=3: 27 Elemente, für n=4: 256 Elemente.

    @param n Größe der Menge {1,...,n}
    @return  Monoid-Instanz
    """
    if n > 4:
        raise ValueError("n > 4 nicht sinnvoll (n^n Elemente werden zu viele)")

    base = list(range(1, n + 1))
    # Alle Funktionen {1,...,n}→{1,...,n} als Tupel (f(1), f(2), ..., f(n))
    elements = list(itertools.product(base, repeat=n))

    def compose(f: tuple, g: tuple) -> tuple:
        """
        Komposition f∘g: (f∘g)(i) = f(g(i)).
        Elemente sind 1-indiziert, Tupel 0-indiziert.
        """
        return tuple(f[g[i - 1] - 1] for i in range(1, n + 1))

    return Monoid(elements, compose)
