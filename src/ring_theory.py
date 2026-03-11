"""
@file ring_theory.py
@brief Ringtheorie – Ringe, Ideale, Quotientenringe und Polynomringe über endlichen Körpern.
@description
    Implementiert die grundlegenden algebraischen Strukturen der Ringtheorie:

    - Ring: Eine Menge mit Addition und Multiplikation (abstrakt über ℤ/nℤ realisiert)
    - Ideal: Teilmengen I ⊆ R mit r·a ∈ I für alle r ∈ R
    - QuotientRing: R/I als Ring der Nebenklassen
    - PolynomialRingModP: ℤ_p[x]/(f(x)) – endliche Körper der Ordnung p^n

    Zentrale Konzepte:
        - Primideale:   R/P ist Integritätsbereich  ⟺  P ist Primideal
        - Maximale Ideale: R/M ist Körper           ⟺  M ist maximales Ideal
        - Hauptideal:   Wird von einem einzigen Element erzeugt
        - Noethersche Ringe: Jede aufsteigende Idealkette stabilisiert
        - Krull-Dimension: Länge der längsten Primidealkette

    Polynomringe:
        - ℤ_p[x]/(f) ist genau dann ein Körper, wenn f irreduzibel über ℤ_p ist
        - Berlekamp-Algorithmus für Faktorisierung über endlichen Körpern

    Zahlenringe:
        - ℤ[√d] = {a + b√d : a,b ∈ ℤ} quadratischer Zahlkörper
        - Gauß'sche Zahlen ℤ[i]: PIR (Hauptidealring), euklidisch

    Anwendungen:
        - Kryptographie (Berlekamp, elliptische Kurven über GF(p^n))
        - Kodierungstheorie (zyklische Codes = Ideale in ℤ_p[x]/(x^n-1))
        - Algebraische Zahlentheorie (Klassenzahl, Einheitengruppe)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import math
from typing import Optional, Union
from math import gcd


# ===========================================================================
# HILFSFUNKTIONEN – POLYNOMOPERATIONEN ÜBER ℤ_p
# ===========================================================================

def _poly_mod(poly: list[int], p: int) -> list[int]:
    """
    Reduziert alle Koeffizienten eines Polynoms modulo p und entfernt führende Nullen.

    @param poly: Koeffizientenliste [a_n, a_{n-1}, ..., a_0] (höchster Grad zuerst)
    @param p: Primzahl (Charakteristik des Körpers)
    @return: Reduzierte Koeffizientenliste
    @lastModified: 2026-03-10
    """
    # Alle Koeffizienten modulo p reduzieren
    result = [c % p for c in poly]
    # Führende Nullen entfernen (Normalform)
    while len(result) > 1 and result[0] == 0:
        result.pop(0)
    return result


def _poly_degree(poly: list[int]) -> int:
    """
    Bestimmt den Grad eines Polynoms.

    @param poly: Koeffizientenliste [a_n, ..., a_0]
    @return: Grad des Polynoms (Länge - 1), -1 für Nullpolynom
    @lastModified: 2026-03-10
    """
    # Normalisierung: führende Nullen ignorieren
    p = [c for c in poly]
    while len(p) > 1 and p[0] == 0:
        p.pop(0)
    # Nullpolynom hat Grad -1 per Konvention
    if len(p) == 1 and p[0] == 0:
        return -1
    return len(p) - 1


def _poly_add_mod(f: list[int], g: list[int], p: int) -> list[int]:
    """
    Addiert zwei Polynome in ℤ_p[x].

    @param f: Erstes Polynom (höchster Grad zuerst)
    @param g: Zweites Polynom (höchster Grad zuerst)
    @param p: Primzahlmodul
    @return: Summe f + g in ℤ_p[x]
    @lastModified: 2026-03-10
    """
    # Gleiche Länge durch Auffüllen mit Nullen
    n = max(len(f), len(g))
    f_pad = [0] * (n - len(f)) + list(f)
    g_pad = [0] * (n - len(g)) + list(g)
    result = [(a + b) % p for a, b in zip(f_pad, g_pad)]
    return _poly_mod(result, p)


def _poly_mul_mod(f: list[int], g: list[int], p: int) -> list[int]:
    """
    Multipliziert zwei Polynome in ℤ_p[x].

    @param f: Erstes Polynom
    @param g: Zweites Polynom
    @param p: Primzahlmodul
    @return: Produkt f·g in ℤ_p[x]
    @lastModified: 2026-03-10
    """
    df = _poly_degree(f)
    dg = _poly_degree(g)
    # Nullpolynom * beliebig = Nullpolynom
    if df < 0 or dg < 0:
        return [0]
    # Koeffizienten des Produkts berechnen (Faltung)
    result = [0] * (df + dg + 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            result[i + j] = (result[i + j] + a * b) % p
    return _poly_mod(result, p)


def _poly_divmod_mod(f: list[int], g: list[int], p: int) -> tuple[list[int], list[int]]:
    """
    Polynomdivision mit Rest in ℤ_p[x]: f = q·g + r.

    @param f: Dividend
    @param g: Divisor (darf nicht Nullpolynom sein)
    @param p: Primzahlmodul
    @return: Tupel (Quotient q, Rest r)
    @raises ValueError: Wenn g das Nullpolynom ist
    @lastModified: 2026-03-10
    """
    if _poly_degree(g) < 0:
        raise ValueError("Division durch das Nullpolynom ist nicht erlaubt.")

    f = _poly_mod(f, p)
    g = _poly_mod(g, p)
    dg = _poly_degree(g)

    # Inverse des Leitkoeffizienten von g in ℤ_p
    lc_g_inv = pow(g[0], p - 2, p)  # Fermats kleiner Satz: a^{p-1} ≡ 1 (mod p)

    # Quotient-Array: Grad(f) - Grad(g) + 1 Koeffizienten (Platzhalter 0)
    df = _poly_degree(f)
    if df < dg:
        # f hat kleineren Grad als g → Quotient = 0, Rest = f
        return [0], f

    quotient_len = df - dg + 1
    quotient = [0] * quotient_len
    remainder = list(f)

    for step in range(quotient_len):
        # Aktuelle führende Stelle im Rest
        # Die führende Stelle ist remainder[step] nach step Eliminierungsschritten
        # Wir arbeiten mit dem aktuellen remainder direkt
        # Leitkoeffizient des aktuellen Rests normalisiert
        lc = remainder[0] if remainder else 0
        coeff = (lc * lc_g_inv) % p
        quotient[step] = coeff

        if coeff == 0:
            # Kein Beitrag: führende Null entfernen
            remainder.pop(0)
            continue

        # Subtrahiere coeff * g vom Rest, ausgerichtet am aktuellen Leitterm
        for i, c in enumerate(g):
            remainder[i] = (remainder[i] - coeff * c) % p

        # Führende Null entfernen (wird durch Subtraktion zu 0)
        remainder.pop(0)

    # Rest normalisieren
    if not remainder:
        remainder = [0]

    return _poly_mod(quotient, p), _poly_mod(remainder, p)


def _poly_eval_mod(poly: list[int], x: int, p: int) -> int:
    """
    Wertet ein Polynom in ℤ_p an der Stelle x aus (Horner-Schema).

    @param poly: Koeffizientenliste (höchster Grad zuerst)
    @param x: Auswertungspunkt
    @param p: Primzahlmodul
    @return: f(x) mod p
    @lastModified: 2026-03-10
    """
    result = 0
    for c in poly:
        result = (result * x + c) % p
    return result


# ===========================================================================
# KLASSE: Ring (ℤ/nℤ)
# ===========================================================================

class Ring:
    """
    Abstrakte Darstellung des Restklassenrings ℤ/nℤ.

    Ein Ring (R, +, ·) ist eine Menge mit zwei Verknüpfungen:
        - (R, +) ist eine abelsche Gruppe (kommutativ, mit 0)
        - (R, ·) ist ein Monoid (assoziativ, mit 1)
        - Distributivgesetze: a(b+c) = ab + ac und (a+b)c = ac + bc

    Beispiele:
        ℤ/2ℤ = {0, 1}           – Körper (2 Elemente)
        ℤ/6ℤ = {0,1,2,3,4,5}   – nicht Integritätsbereich (2·3 = 0)
        ℤ/7ℤ = {0,1,...,6}      – Körper (7 Primzahl)

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, n: int) -> None:
        """
        Erstellt den Ring ℤ/nℤ.

        @param n: Modulus (muss ≥ 2 sein)
        @raises ValueError: Wenn n < 2
        @lastModified: 2026-03-10
        """
        if n < 2:
            raise ValueError(f"Modulus muss ≥ 2 sein, erhalten: {n}")
        # Modulus speichern (Ordnung des Rings)
        self.n = n
        # Alle Elemente {0, 1, ..., n-1}
        self.elements: list[int] = list(range(n))

    def __repr__(self) -> str:
        return f"Ring(ℤ/{self.n}ℤ)"

    def add(self, a: int, b: int) -> int:
        """
        Addition in ℤ/nℤ: a + b mod n.

        @param a: Erstes Element
        @param b: Zweites Element
        @return: (a + b) mod n
        @lastModified: 2026-03-10
        """
        return (a + b) % self.n

    def mul(self, a: int, b: int) -> int:
        """
        Multiplikation in ℤ/nℤ: a · b mod n.

        @param a: Erstes Element
        @param b: Zweites Element
        @return: (a · b) mod n
        @lastModified: 2026-03-10
        """
        return (a * b) % self.n

    def is_commutative(self) -> bool:
        """
        Prüft ob der Ring kommutativ ist (a·b = b·a für alle a,b).

        ℤ/nℤ ist immer kommutativ, da die ganzzahlige Multiplikation kommutativ ist.

        @return: True (ℤ/nℤ ist stets kommutativ)
        @lastModified: 2026-03-10
        """
        # ℤ/nℤ ist als Restklassenring von ℤ immer kommutativ
        return True

    def characteristic(self) -> int:
        """
        Bestimmt die Charakteristik des Rings.

        Die Charakteristik char(R) ist die kleinste positive ganze Zahl n mit
        n·1 = 0 (n-fache Addition des Einselements ergibt 0).
        Gibt 0 zurück, wenn keine solche Zahl existiert.

        Für ℤ/nℤ gilt: char(ℤ/nℤ) = n.

        @return: Charakteristik des Rings
        @lastModified: 2026-03-10
        """
        # Für ℤ/nℤ: n·1 = n ≡ 0 (mod n), also char = n
        return self.n

    def units(self) -> list[int]:
        """
        Bestimmt alle Einheiten (invertierbare Elemente) des Rings.

        Ein Element a ∈ R ist eine Einheit, wenn es ein b ∈ R gibt mit a·b = 1.
        In ℤ/nℤ: a ist Einheit ⟺ gcd(a, n) = 1 (Satz von Euler).

        @return: Liste aller Einheiten in ℤ/nℤ
        @lastModified: 2026-03-10
        """
        # a ist Einheit genau dann wenn ggT(a, n) = 1
        return [a for a in self.elements if a > 0 and math.gcd(a, self.n) == 1]

    def zero_divisors(self) -> list[int]:
        """
        Bestimmt alle Nullteiler des Rings (ohne 0 selbst).

        Ein Element a ≠ 0 ist ein Nullteiler, wenn es ein b ≠ 0 gibt mit a·b = 0.
        In ℤ/nℤ: a ist Nullteiler ⟺ gcd(a, n) > 1 und a ≠ 0.

        @return: Liste aller Nullteiler (ohne 0)
        @lastModified: 2026-03-10
        """
        divisors = []
        for a in self.elements:
            if a == 0:
                continue
            # Suche b ≠ 0 mit a·b ≡ 0 (mod n)
            for b in self.elements:
                if b != 0 and (a * b) % self.n == 0:
                    divisors.append(a)
                    break
        return divisors

    def is_integral_domain(self) -> bool:
        """
        Prüft ob der Ring ein Integritätsbereich ist.

        Ein kommutativer Ring ohne Nullteiler ist ein Integritätsbereich.
        Für ℤ/nℤ gilt: Integritätsbereich ⟺ n ist eine Primzahl.

        @return: True wenn ℤ/nℤ ein Integritätsbereich ist
        @lastModified: 2026-03-10
        """
        # Integritätsbereich ⟺ keine Nullteiler ⟺ n prim
        return len(self.zero_divisors()) == 0

    def is_field(self) -> bool:
        """
        Prüft ob der Ring ein Körper ist.

        Ein Körper ist ein kommutativer Ring, in dem jedes Element ≠ 0 invertierbar ist.
        Für ℤ/nℤ gilt: Körper ⟺ n ist eine Primzahl.

        @return: True wenn ℤ/nℤ ein Körper ist
        @lastModified: 2026-03-10
        """
        # Körper ⟺ jedes Nicht-Null-Element ist Einheit ⟺ n prim
        units = self.units()
        non_zero = [a for a in self.elements if a != 0]
        return set(units) == set(non_zero)

    def is_euclidean_domain(self) -> bool:
        """
        Prüft ob der Ring ein euklidischer Bereich ist.

        Ein Integritätsbereich R heißt euklidisch, wenn es eine euklidische Norm
        φ: R\\{0} → ℕ gibt, sodass Division mit Rest möglich ist.
        Jeder Körper ist trivialerweise euklidisch.
        Für ℤ/nℤ: euklidisch ⟺ n prim (d.h. ℤ/nℤ ist Körper).

        @return: True wenn der Ring euklidisch ist
        @lastModified: 2026-03-10
        """
        # ℤ/pℤ (p prim) ist Körper, also euklidisch
        return self.is_field()

    def is_pid(self) -> bool:
        """
        Prüft ob der Ring ein Hauptidealring (PID) ist.

        Ein Hauptidealring ist ein Integritätsbereich, in dem jedes Ideal ein
        Hauptideal ist (von einem Element erzeugt wird).
        Jeder euklidische Bereich ist ein PID.
        Für ℤ/pℤ: PID ⟺ p prim.

        @return: True wenn der Ring ein Hauptidealring ist
        @lastModified: 2026-03-10
        """
        # Körper sind PIDs (jedes Ideal ist {0} oder R, beide Hauptideale)
        return self.is_field()

    def is_ufd(self) -> bool:
        """
        Prüft ob der Ring ein faktorieller Ring (UFD) ist.

        In einem UFD lässt sich jedes Element eindeutig (bis auf Reihenfolge und
        Einheiten) als Produkt irreduzibler Elemente schreiben.
        Jeder PID ist ein UFD.

        @return: True wenn der Ring ein UFD ist
        @lastModified: 2026-03-10
        """
        # Jeder PID ist UFD, Körper sind PIDs → Körper sind UFDs
        return self.is_pid()


# ===========================================================================
# KLASSE: Ideal
# ===========================================================================

class Ideal:
    """
    Ein Ideal I in einem Ring R = ℤ/nℤ.

    Ein Ideal I ⊆ R erfüllt:
        1. (I, +) ist eine Untergruppe von (R, +)
        2. Für alle r ∈ R und a ∈ I gilt: r·a ∈ I (Linksideal)
           und a·r ∈ I (Rechtsideal) – für kommutative Ringe identisch

    In ℤ/nℤ sind alle Ideale Hauptideale der Form (d) = {0, d, 2d, ...}
    wobei d ein Teiler von n ist.

    Wichtige Typen:
        - Primideal P:    R/P ist Integritätsbereich
        - Maximales Ideal M: R/M ist Körper
        - Hauptideal:     wird von einem Element erzeugt

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, ring: Ring, elements: list[int]) -> None:
        """
        Erstellt ein Ideal mit den angegebenen Elementen.

        @param ring: Der Grundring R = ℤ/nℤ
        @param elements: Liste der Elemente im Ideal
        @raises ValueError: Wenn die Menge kein gültiges Ideal ist
        @lastModified: 2026-03-10
        """
        self.ring = ring
        # Normalisierung: Elemente mod n, sortiert, ohne Duplikate
        self.elements: list[int] = sorted(set(x % ring.n for x in elements))

        # Validierung: Ideal muss die 0 enthalten
        if 0 not in self.elements:
            raise ValueError("Ein Ideal muss die Null (0) enthalten.")

        # Validierung: Abgeschlossenheit unter Addition
        for a in self.elements:
            for b in self.elements:
                s = (a + b) % ring.n
                if s not in self.elements:
                    raise ValueError(
                        f"Nicht unter Addition abgeschlossen: {a} + {b} = {s} ∉ Ideal"
                    )

        # Validierung: Absorptionseigenschaft r·a ∈ I für alle r ∈ R, a ∈ I
        for r in ring.elements:
            for a in self.elements:
                ra = (r * a) % ring.n
                if ra not in self.elements:
                    raise ValueError(
                        f"Absorptionseigenschaft verletzt: {r}·{a} = {ra} ∉ Ideal"
                    )

    def __repr__(self) -> str:
        return f"Ideal({self.elements} in {self.ring})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Ideal):
            return False
        return self.ring.n == other.ring.n and set(self.elements) == set(other.elements)

    def __le__(self, other: 'Ideal') -> bool:
        """Enthaltensein: I ⊆ J."""
        return all(x in other.elements for x in self.elements)

    def is_prime(self) -> bool:
        """
        Prüft ob das Ideal ein Primideal ist.

        P ist Primideal ⟺ R/P ist Integritätsbereich
                        ⟺ a·b ∈ P  →  a ∈ P oder b ∈ P

        In ℤ/nℤ: Das Ideal (d) ist prim ⟺ n/d eine Primzahl ist.

        @return: True wenn das Ideal ein Primideal ist
        @lastModified: 2026-03-10
        """
        n = self.ring.n
        # Primideal-Bedingung: ab ∈ I → a ∈ I oder b ∈ I
        ideal_set = set(self.elements)
        for a in self.ring.elements:
            for b in self.ring.elements:
                if (a * b) % n in ideal_set:
                    # a·b ∈ I: dann muss a ∈ I oder b ∈ I gelten
                    if a not in ideal_set and b not in ideal_set:
                        return False
        return True

    def is_maximal(self) -> bool:
        """
        Prüft ob das Ideal ein maximales Ideal ist.

        M ist maximal ⟺ R/M ist Körper
                      ⟺ kein Ideal J mit M ⊊ J ⊊ R existiert

        In ℤ/nℤ: Jedes maximale Ideal ist prim (da ℤ/nℤ noethersch ist).

        @return: True wenn das Ideal maximal ist
        @lastModified: 2026-03-10
        """
        n = self.ring.n
        # Maximales Ideal: kein echtes Zwischenideal M ⊊ J ⊊ R
        # In ℤ/nℤ: Ideale entsprechen Teilern von n, also (d) für d | n
        # (d) ist maximal ⟺ n/d ist prim
        ideal_set = set(self.elements)
        ring_set = set(self.ring.elements)

        # Suche echtes Zwischenideal
        for d in range(1, n):
            if n % d == 0:
                # Erzeuge Ideal (d) in ℤ/nℤ
                candidate = set(range(0, n, d))
                # Echtes Zwischenideal: ideal_set ⊊ candidate ⊊ ring_set
                if ideal_set < candidate < ring_set:
                    return False
        return True

    def is_principal(self) -> bool:
        """
        Prüft ob das Ideal ein Hauptideal ist.

        Ein Ideal I ist ein Hauptideal, wenn es ein a ∈ R gibt mit I = (a) = {r·a : r ∈ R}.
        In ℤ/nℤ sind alle Ideale Hauptideale (ℤ/nℤ ist noethersch und prim).

        @return: True wenn das Ideal ein Hauptideal ist
        @lastModified: 2026-03-10
        """
        n = self.ring.n
        ideal_set = set(self.elements)
        # Prüfe für jedes Element, ob es das Ideal erzeugt
        for a in self.elements:
            # Erzeuge (a) = {0·a, 1·a, ..., (n-1)·a} mod n
            generated = set((r * a) % n for r in range(n))
            if generated == ideal_set:
                return True
        return False


# ===========================================================================
# KLASSE: QuotientRing
# ===========================================================================

class QuotientRing(Ring):
    """
    Quotientenring R/I für einen Ring R = ℤ/nℤ und ein Ideal I.

    Die Elemente von R/I sind Nebenklassen a + I = {a + x : x ∈ I}.
    Die Arithmetik wird komponentenweise definiert:
        (a + I) + (b + I) = (a+b) + I
        (a + I) · (b + I) = (a·b) + I

    Eigenschaften:
        - R/P ist Integritätsbereich  ⟺  P ist Primideal
        - R/M ist Körper              ⟺  M ist maximales Ideal

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, ring: Ring, ideal: Ideal) -> None:
        """
        Erstellt den Quotientenring R/I.

        @param ring: Grundring R = ℤ/nℤ
        @param ideal: Ideal I ⊆ R
        @lastModified: 2026-03-10
        """
        self.base_ring = ring
        self.ideal = ideal

        # Nebenklassen berechnen: eine Repräsentanten-Menge
        # Wähle für jede Nebenklasse den kleinsten Vertreter
        ideal_set = set(ideal.elements)
        seen_classes: set[frozenset] = set()
        representatives = []

        for a in ring.elements:
            # Nebenklasse von a: {a + x : x ∈ I} mod n
            coset = frozenset((a + x) % ring.n for x in ideal_set)
            if coset not in seen_classes:
                seen_classes.add(coset)
                # Kleinsten Vertreter wählen
                representatives.append(min(coset))

        # Superklasse initialisieren mit Anzahl der Nebenklassen
        # (entspricht |R/I| = |R| / |I|)
        n_quotient = ring.n // len(ideal.elements)
        super().__init__(n_quotient)

        # Repräsentanten und Zuordnung speichern
        self._representatives = sorted(representatives)
        self._ideal_size = len(ideal.elements)

    def __repr__(self) -> str:
        return f"QuotientRing({self.base_ring}/{self.ideal})"

    def coset_rep(self, a: int) -> int:
        """
        Gibt den kanonischen Vertreter der Nebenklasse a + I zurück.

        @param a: Element aus dem Grundring
        @return: Kleinster nicht-negativer Vertreter der Nebenklasse
        @lastModified: 2026-03-10
        """
        ideal_set = set(self.ideal.elements)
        # Nebenklasse von a: alle (a + x) mod n für x ∈ I
        coset = {(a + x) % self.base_ring.n for x in ideal_set}
        return min(coset)

    def is_field(self) -> bool:
        """
        R/I ist Körper genau dann wenn I ein maximales Ideal ist.

        @return: True wenn R/I ein Körper ist
        @lastModified: 2026-03-10
        """
        return self.ideal.is_maximal()

    def is_integral_domain(self) -> bool:
        """
        R/I ist Integritätsbereich genau dann wenn I ein Primideal ist.

        @return: True wenn R/I ein Integritätsbereich ist
        @lastModified: 2026-03-10
        """
        return self.ideal.is_prime()


# ===========================================================================
# KLASSE: PolynomialRingModP
# ===========================================================================

class PolynomialRingModP:
    """
    Polynomring ℤ_p[x]/(f(x)) – Restklassenring modulo einem Polynom.

    Ist f(x) irreduzibel über ℤ_p, so ist ℤ_p[x]/(f(x)) ein endlicher Körper
    der Ordnung p^deg(f), auch bezeichnet als GF(p^n) oder F_{p^n}.

    Beispiel: ℤ_2[x]/(x²+x+1) ist der Körper GF(4) mit 4 Elementen:
        {0, 1, x, x+1}  (Darstellung als Polynome mod x²+x+1 über ℤ_2)

    Arithmetik:
        - Addition: komponentenweise mod p
        - Multiplikation: Polynomprodukt, dann mod f(x) und mod p

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, modulus_poly: list[int], p: int) -> None:
        """
        Erstellt ℤ_p[x]/(f(x)).

        @param modulus_poly: Koeffizienten von f(x) (höchster Grad zuerst), z.B. [1,1,1] für x²+x+1
        @param p: Primzahl (Charakteristik)
        @raises ValueError: Wenn p keine Primzahl ist
        @lastModified: 2026-03-10
        """
        # Primzahlprüfung
        if p < 2 or not all(p % i != 0 for i in range(2, int(p**0.5) + 1)):
            raise ValueError(f"{p} ist keine Primzahl.")

        self.p = p
        # Polynom mod p normalisieren
        self.modulus = _poly_mod(modulus_poly, p)
        self.degree = _poly_degree(self.modulus)

        if self.degree < 1:
            raise ValueError("Das Moduluspolynom muss Grad ≥ 1 haben.")

        # Ordnung des Körpers: p^deg(f)
        self.order = p ** self.degree

    def __repr__(self) -> str:
        return f"PolynomialRingModP(ℤ_{self.p}[x]/({self._poly_str(self.modulus)}))"

    def _poly_str(self, poly: list[int]) -> str:
        """Lesbare Darstellung eines Polynoms."""
        if not poly or (len(poly) == 1 and poly[0] == 0):
            return "0"
        terms = []
        deg = len(poly) - 1
        for i, c in enumerate(poly):
            if c == 0:
                continue
            power = deg - i
            if power == 0:
                terms.append(str(c))
            elif power == 1:
                terms.append(f"{c}x" if c != 1 else "x")
            else:
                terms.append(f"{c}x^{power}" if c != 1 else f"x^{power}")
        return " + ".join(terms) if terms else "0"

    def reduce(self, poly: list[int]) -> list[int]:
        """
        Reduziert ein Polynom modulo f(x) in ℤ_p[x].

        @param poly: Zu reduzierendes Polynom
        @return: Rest bei Division durch f(x) (Normalform)
        @lastModified: 2026-03-10
        """
        _, rem = _poly_divmod_mod(poly, self.modulus, self.p)
        return rem

    def add(self, f: list[int], g: list[int]) -> list[int]:
        """
        Addition in ℤ_p[x]/(f(x)).

        @param f: Erstes Polynom (bereits reduziert)
        @param g: Zweites Polynom (bereits reduziert)
        @return: f + g mod Moduluspolynom, mod p
        @lastModified: 2026-03-10
        """
        return self.reduce(_poly_add_mod(f, g, self.p))

    def mul(self, f: list[int], g: list[int]) -> list[int]:
        """
        Multiplikation in ℤ_p[x]/(f(x)).

        @param f: Erstes Polynom
        @param g: Zweites Polynom
        @return: f·g mod Moduluspolynom, mod p
        @lastModified: 2026-03-10
        """
        return self.reduce(_poly_mul_mod(f, g, self.p))

    def is_field(self) -> bool:
        """
        Prüft ob ℤ_p[x]/(f) ein Körper ist.

        ℤ_p[x]/(f) ist Körper ⟺ f irreduzibel über ℤ_p.

        @return: True wenn der Ring ein Körper ist
        @lastModified: 2026-03-10
        """
        return is_irreducible_mod_p(self.modulus, self.p)

    def elements(self) -> list[list[int]]:
        """
        Gibt alle Elemente von ℤ_p[x]/(f) zurück.

        Jedes Element ist ein Polynom vom Grad < deg(f) mit Koeffizienten in ℤ_p.
        Anzahl der Elemente: p^deg(f).

        @return: Liste aller Polynome (als Koeffizientenlisten)
        @lastModified: 2026-03-10
        """
        from itertools import product as iproduct
        result = []
        # Alle Polynome a_{n-1}x^{n-1} + ... + a_0 mit a_i ∈ {0,...,p-1}
        for coeffs in iproduct(range(self.p), repeat=self.degree):
            # Führende Nullen entfernen für kanonische Darstellung
            elem = _poly_mod(list(coeffs), self.p)
            result.append(elem)
        return result


# ===========================================================================
# FUNKTION: principal_ideal
# ===========================================================================

def principal_ideal(ring: Ring, generator: int) -> 'Ideal':
    """
    Erzeugt das Hauptideal (a) = {r·a : r ∈ R} in ℤ/nℤ.

    Das Hauptideal wird von einem einzigen Element a erzeugt.
    In ℤ/nℤ: (a) = {0, a, 2a, ..., (k-1)a} wobei k = n/gcd(a,n).

    @param ring: Grundring R = ℤ/nℤ
    @param generator: Erzeuger a des Ideals
    @return: Das Hauptideal (a)
    @lastModified: 2026-03-10
    """
    n = ring.n
    g = generator % n
    # (a) = {r·a mod n : r ∈ ℤ/nℤ}
    elements = sorted(set((r * g) % n for r in range(n)))
    return Ideal(ring, elements)


# ===========================================================================
# FUNKTION: ideal_sum
# ===========================================================================

def ideal_sum(I: Ideal, J: Ideal) -> Ideal:
    """
    Berechnet die Summe zweier Ideale: I + J = {a + b : a ∈ I, b ∈ J}.

    Die Summe I + J ist das kleinste Ideal, das beide Ideale enthält.
    In ℤ/nℤ: (d₁) + (d₂) = (gcd(d₁, d₂)).

    @param I: Erstes Ideal
    @param J: Zweites Ideal
    @return: Summe I + J als Ideal
    @raises ValueError: Wenn I und J in verschiedenen Ringen liegen
    @lastModified: 2026-03-10
    """
    if I.ring.n != J.ring.n:
        raise ValueError("Ideale müssen im selben Ring liegen.")
    n = I.ring.n
    # I + J = {a + b mod n : a ∈ I, b ∈ J}
    elements = sorted(set((a + b) % n for a in I.elements for b in J.elements))
    return Ideal(I.ring, elements)


# ===========================================================================
# FUNKTION: ideal_product
# ===========================================================================

def ideal_product(I: Ideal, J: Ideal) -> Ideal:
    """
    Berechnet das Produkt zweier Ideale: I·J = kleinster Ideal ⊇ {a·b : a∈I, b∈J}.

    Das Produkt I·J wird erzeugt von allen endlichen Summen aᵢ·bᵢ mit aᵢ ∈ I, bᵢ ∈ J.
    In der Praxis für endliche Ringe: schließe unter Addition ab.

    @param I: Erstes Ideal
    @param J: Zweites Ideal
    @return: Produkt I·J als Ideal
    @raises ValueError: Wenn I und J in verschiedenen Ringen liegen
    @lastModified: 2026-03-10
    """
    if I.ring.n != J.ring.n:
        raise ValueError("Ideale müssen im selben Ring liegen.")
    n = I.ring.n
    # Basis: alle paarweisen Produkte {a·b : a ∈ I, b ∈ J}
    products = set((a * b) % n for a in I.elements for b in J.elements)
    # Unter Addition abschließen (erzeugte Untergruppe)
    changed = True
    while changed:
        changed = False
        for a in list(products):
            for b in list(products):
                s = (a + b) % n
                if s not in products:
                    products.add(s)
                    changed = True
    return Ideal(I.ring, sorted(products))


# ===========================================================================
# FUNKTION: ideal_intersection
# ===========================================================================

def ideal_intersection(I: Ideal, J: Ideal) -> Ideal:
    """
    Berechnet den Schnitt zweier Ideale: I ∩ J.

    Der Schnitt zweier Ideale ist wieder ein Ideal.
    In ℤ/nℤ: (d₁) ∩ (d₂) = (kgV(d₁, d₂)).

    @param I: Erstes Ideal
    @param J: Zweites Ideal
    @return: Schnittideal I ∩ J
    @raises ValueError: Wenn I und J in verschiedenen Ringen liegen
    @lastModified: 2026-03-10
    """
    if I.ring.n != J.ring.n:
        raise ValueError("Ideale müssen im selben Ring liegen.")
    # Mengentheoretischer Schnitt
    elements = sorted(set(I.elements) & set(J.elements))
    return Ideal(I.ring, elements)


# ===========================================================================
# FUNKTION: chinese_remainder_ring
# ===========================================================================

def chinese_remainder_ring(n: int, moduli: list[int]) -> dict:
    """
    Chinesischer Restsatz für Ringe: ℤ/nℤ ≅ Π ℤ/pᵢ^aᵢℤ.

    Der CRT besagt: Wenn n = p₁^a₁ · p₂^a₂ · ... · pₖ^aₖ (Primfaktorzerlegung),
    dann gilt: ℤ/nℤ ≅ ℤ/p₁^a₁ℤ × ℤ/p₂^a₂ℤ × ... × ℤ/pₖ^aₖℤ

    als Ringisomorphismus (komponentenweise Addition und Multiplikation).

    Die Abbildung φ: ℤ/nℤ → Π ℤ/mᵢℤ ist gegeben durch
        φ(a) = (a mod m₁, a mod m₂, ..., a mod mₖ)

    @param n: Modulus des Ausgangsrings
    @param moduli: Liste der Teilmoduln mᵢ (müssen paarweise teilerfremd sein)
    @return: Dictionary mit CRT-Informationen und Isomorphismus-Tabelle
    @raises ValueError: Wenn die Moduln nicht paarweise teilerfremd sind
    @lastModified: 2026-03-10
    """
    # Prüfe paarweise Teilerfremdheit
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            if math.gcd(moduli[i], moduli[j]) != 1:
                raise ValueError(
                    f"Moduln müssen paarweise teilerfremd sein: "
                    f"gcd({moduli[i]}, {moduli[j]}) = {math.gcd(moduli[i], moduli[j])}"
                )

    # Prüfe ob das Produkt gleich n ist
    product = 1
    for m in moduli:
        product *= m
    if product != n:
        raise ValueError(
            f"Produkt der Moduln ({product}) muss gleich n ({n}) sein."
        )

    # CRT-Isomorphismus-Tabelle aufbauen
    # Für jedes a ∈ ℤ/nℤ: (a mod m₁, a mod m₂, ..., a mod mₖ)
    iso_map = {}
    for a in range(n):
        iso_map[a] = tuple(a % m for m in moduli)

    # Inverse: Rekonstruktion via CRT-Koeffizienten
    # Berechne Ni = n/mi und deren Inverse mod mi
    crt_coefficients = []
    for m in moduli:
        Ni = n // m
        # Ni^{-1} mod mi via erweitertem euklidischem Algorithmus
        Ni_inv = pow(Ni, -1, m)
        crt_coefficients.append((Ni, Ni_inv))

    return {
        'n': n,
        'moduli': moduli,
        'isomorphism': iso_map,
        'crt_coefficients': crt_coefficients,
        'description': f"ℤ/{n}ℤ ≅ " + " × ".join(f"ℤ/{m}ℤ" for m in moduli)
    }


# ===========================================================================
# FUNKTION: is_irreducible_mod_p
# ===========================================================================

def is_irreducible_mod_p(poly_coeffs: list[int], p: int) -> bool:
    """
    Prüft ob f(x) ∈ ℤ_p[x] irreduzibel ist.

    Ein Polynom f ist irreduzibel, wenn es keine nicht-triviale Faktorisierung
    f = g·h gibt mit 0 < deg(g), deg(h) < deg(f).

    Algorithmus für deg ≤ 10 (Brute-Force über alle Polynome bis Grad ⌊deg/2⌋):
        1. Prüfe auf Nullstellen in ℤ_p (lineare Faktoren)
        2. Prüfe auf quadratische Faktoren (falls deg ≥ 4)
        3. Allgemein: teste alle Polynome vom Grad 1 bis ⌊deg/2⌋

    @param poly_coeffs: Koeffizienten von f(x) (höchster Grad zuerst)
    @param p: Primzahl
    @return: True wenn f irreduzibel über ℤ_p ist
    @lastModified: 2026-03-10
    """
    from itertools import product as iproduct

    poly = _poly_mod(poly_coeffs, p)
    deg = _poly_degree(poly)

    # Konstante oder lineare Polynome sind (trivialerweise) irreduzibel
    if deg <= 1:
        return deg == 1

    # Schritt 1: Prüfe auf Nullstellen (lineare Faktoren)
    for x in range(p):
        if _poly_eval_mod(poly, x, p) == 0:
            return False  # (x - x₀) ist Faktor

    # Schritt 2: Prüfe auf Faktoren vom Grad 2 bis deg//2
    max_factor_deg = deg // 2
    for factor_deg in range(2, max_factor_deg + 1):
        # Iteriere über alle monischen Polynome vom Grad factor_deg über ℤ_p
        for lower_coeffs in iproduct(range(p), repeat=factor_deg):
            # Monisches Polynom: x^d + a_{d-1}x^{d-1} + ... + a_0
            candidate = [1] + list(lower_coeffs)
            candidate = _poly_mod(candidate, p)
            # Prüfe ob candidate f teilt (Rest = 0)
            _, rem = _poly_divmod_mod(poly, candidate, p)
            if _poly_degree(rem) < 0 or (len(rem) == 1 and rem[0] == 0):
                return False  # Echter Faktor gefunden

    return True


# ===========================================================================
# FUNKTION: polynomial_gcd_mod_p
# ===========================================================================

def polynomial_gcd_mod_p(f: list[int], g: list[int], p: int) -> list[int]:
    """
    Berechnet den größten gemeinsamen Teiler ggT(f, g) in ℤ_p[x].

    Verwendet den euklidischen Algorithmus für Polynome:
        ggT(f, g) = ggT(g, f mod g)
    bis der Rest 0 ist. Normalisiert auf monisches Polynom (Leitkoeffizient = 1).

    @param f: Erstes Polynom (höchster Grad zuerst)
    @param g: Zweites Polynom
    @param p: Primzahl
    @return: Monischer ggT in ℤ_p[x]
    @lastModified: 2026-03-10
    """
    f = _poly_mod(f, p)
    g = _poly_mod(g, p)

    # Euklidischer Algorithmus für Polynome
    while not (_poly_degree(g) < 0 or (len(g) == 1 and g[0] == 0)):
        _, rem = _poly_divmod_mod(f, g, p)
        f = g
        g = rem

    # Auf monisches Polynom normalisieren (Leitkoeffizient = 1)
    if _poly_degree(f) >= 0 and f[0] != 0:
        lc_inv = pow(f[0], p - 2, p)
        f = [(c * lc_inv) % p for c in f]

    return _poly_mod(f, p)


# ===========================================================================
# FUNKTION: polynomial_factorization_mod_p (Berlekamp-Algorithmus)
# ===========================================================================

def polynomial_factorization_mod_p(poly: list[int], p: int) -> dict:
    """
    Faktorisiert f(x) ∈ ℤ_p[x] mit dem Berlekamp-Algorithmus.

    Der Berlekamp-Algorithmus ist ein effizienter Algorithmus zur Faktorisierung
    quadratfreier Polynome über endlichen Körpern ℤ_p (p Primzahl):

    Schritte:
        1. Squarefree-Zerlegung: f = f₁ · f₂² · f₃³ · ... (Trennung vielfacher Faktoren)
        2. Berlekamp-Matrix Q aufstellen: Q[i][j] = Koeffizient von x^{i·p} mod f
        3. Kern(Q - I) berechnen → Berlekamp-Subalgebra
        4. Faktoren via GCD mit Translationspolynomen h - c (c ∈ ℤ_p) bestimmen

    @param poly: Koeffizienten von f(x) (höchster Grad zuerst)
    @param p: Primzahl
    @return: {'factors': list[list[int]], 'multiplicities': list[int]}
    @lastModified: 2026-03-10
    """
    poly = _poly_mod(poly, p)
    n = _poly_degree(poly)

    if n <= 0:
        return {'factors': [poly], 'multiplicities': [1]}

    # Schritt 1: Squarefree-Zerlegung
    # f'(x) = formale Ableitung
    def formal_derivative(f: list[int]) -> list[int]:
        """Formale Ableitung: d/dx(Σ aᵢxⁱ) = Σ i·aᵢ·x^{i-1}."""
        deg = len(f) - 1
        if deg == 0:
            return [0]
        result = []
        for i, c in enumerate(f[:-1]):
            power = deg - i
            result.append((power * c) % p)
        return _poly_mod(result, p)

    # Squarefree: f / gcd(f, f')
    df = formal_derivative(poly)
    d_is_zero = _poly_degree(df) < 0 or (len(df) == 1 and df[0] == 0)

    if d_is_zero:
        # Charakteristik p teilt alle Exponenten: f = g(x^p) → g(x) mit Berlekamp
        # Vereinfachte Behandlung: p-te Wurzel der Koeffizienten
        # Für kleine p: Brute-force über alle Faktoren
        return _factorize_brute_force(poly, p)

    g = polynomial_gcd_mod_p(poly, df, p)
    g_deg = _poly_degree(g)

    if g_deg == 0:
        # f ist bereits quadratfrei
        squarefree_part = poly
        multiplicity_factor = 1
    else:
        # Trenne quadratfreien Teil
        squarefree_part, _ = _poly_divmod_mod(poly, g, p)
        multiplicity_factor = 1

    # Schritt 2: Berlekamp-Algorithmus auf den quadratfreien Teil anwenden
    factors_of_sqf = _berlekamp(squarefree_part, p)

    # Multiplizitäten bestimmen
    all_factors = []
    all_mults = []
    remaining = poly

    for fac in factors_of_sqf:
        if _poly_degree(fac) < 1:
            continue
        mult = 0
        # Zähle wie oft fac in remaining vorkommt
        while True:
            q, r = _poly_divmod_mod(remaining, fac, p)
            r_deg = _poly_degree(r)
            if r_deg < 0 or (len(r) == 1 and r[0] == 0):
                remaining = q
                mult += 1
            else:
                break
        if mult > 0:
            all_factors.append(fac)
            all_mults.append(mult)

    # Falls nichts gefunden (f irreduzibel), f selbst zurückgeben
    if not all_factors:
        all_factors = [poly]
        all_mults = [1]

    return {'factors': all_factors, 'multiplicities': all_mults}


def _berlekamp(poly: list[int], p: int) -> list[list[int]]:
    """
    Berlekamp-Algorithmus: Faktorisierung quadratfreier Polynome über ℤ_p.

    Konstruiert die Berlekamp-Matrix Q, deren Kern die Faktorisierung liefert.

    @param poly: Quadratfreies monisches Polynom
    @param p: Primzahl
    @return: Liste der irreduziblen Faktoren
    @lastModified: 2026-03-10
    """
    n = _poly_degree(poly)
    if n <= 0:
        return [poly]
    if n == 1:
        return [poly]

    # Berlekamp-Matrix Q aufbauen: Q[i] = x^{i·p} mod f(x), als Zeilenvektor
    Q = []
    for i in range(n):
        # Berechne x^{i·p} mod f via wiederholtes Quadrieren
        base = [0] * (i + 1)  # x^i: [0,...,0,1,0,...,0]
        base[0] = 1
        base = base + [0] * (n - i - 1)  # Grad n-1
        # Eigentlich: Basis = x^i als Polynom vom Grad n-1 (Koeffizient bei x^i = 1)
        basis = [0] * n
        basis[n - 1 - i] = 1  # Koeffizient von x^i = 1
        # Berechne basis^p mod poly via schneller Exponentiation
        vec = _poly_pow_mod(basis, p, poly, p)
        # Als Vektor [a_0, a_1, ..., a_{n-1}] (Koeffizienten von x^0 bis x^{n-1})
        row = [0] * n
        for j, c in enumerate(vec):
            row[n - 1 - j] = c  # Umkehren: vec[0] = Leitkoeffizient
        Q.append(row)

    # Q - I berechnen (um Kern von Q - I zu finden)
    QmI = [[(Q[i][j] - (1 if i == j else 0)) % p for j in range(n)] for i in range(n)]

    # Kern(Q - I) via Gauß-Elimination über ℤ_p finden
    basis_vectors = _kernel_mod_p(QmI, n, p)

    # Anzahl der irreduziblen Faktoren = dim(Kern(Q-I))
    num_factors = len(basis_vectors)

    if num_factors == 1:
        # f ist irreduzibel
        return [poly]

    # Faktoren aus dem Kern extrahieren
    factors = [poly]
    for h_coeffs in basis_vectors[1:]:  # Ersten Vektor (konstant) überspringen
        # h(x) als Polynom
        h = [h_coeffs[n - 1 - i] for i in range(n - 1, -1, -1)]
        h = _poly_mod(h, p)
        # GCD von f mit (h - c) für c ∈ ℤ_p
        new_factors = []
        for factor in factors:
            for c in range(p):
                # h - c
                hc = list(h)
                if hc:
                    hc[-1] = (hc[-1] - c) % p
                else:
                    hc = [(-c) % p]
                hc = _poly_mod(hc, p)
                g = polynomial_gcd_mod_p(factor, hc, p)
                if _poly_degree(g) > 0 and _poly_degree(g) < _poly_degree(factor):
                    q, _ = _poly_divmod_mod(factor, g, p)
                    new_factors.append(g)
                    new_factors.append(_poly_mod(q, p))
                    break
            else:
                new_factors.append(factor)
        factors = [f for f in new_factors if _poly_degree(f) > 0]

    return factors if factors else [poly]


def _poly_pow_mod(base: list[int], exp: int, mod: list[int], p: int) -> list[int]:
    """
    Schnelle Polynom-Exponentiation: base^exp mod mod_poly in ℤ_p[x].

    @param base: Basispolynom
    @param exp: Exponent
    @param mod: Moduluspolynom
    @param p: Primzahl
    @return: base^exp mod mod_poly in ℤ_p[x]
    @lastModified: 2026-03-10
    """
    result = [1]  # Neutralelement für Multiplikation
    base = _poly_mod(base, p)
    # Reduziere base direkt
    _, base = _poly_divmod_mod(base, mod, p)

    while exp > 0:
        if exp % 2 == 1:
            result = _poly_mul_mod(result, base, p)
            _, result = _poly_divmod_mod(result, mod, p)
        base = _poly_mul_mod(base, base, p)
        _, base = _poly_divmod_mod(base, mod, p)
        exp //= 2
    return result


def _kernel_mod_p(matrix: list[list[int]], n: int, p: int) -> list[list[int]]:
    """
    Berechnet den Kern einer Matrix über ℤ_p via Gauß-Elimination.

    @param matrix: n×n Matrix über ℤ_p
    @param n: Matrixgröße
    @param p: Primzahl
    @return: Basis des Kerns (Liste der Basisvektoren)
    @lastModified: 2026-03-10
    """
    # Arbeite mit erweiterter Matrix [A | I] für Zeilenoperationen
    M = [list(row) for row in matrix]
    pivots = []  # Pivot-Spalten

    row = 0
    for col in range(n):
        # Suche Pivot in dieser Spalte
        pivot_row = None
        for r in range(row, n):
            if M[r][col] != 0:
                pivot_row = r
                break
        if pivot_row is None:
            continue

        # Zeilen tauschen
        M[row], M[pivot_row] = M[pivot_row], M[row]
        pivots.append(col)

        # Pivot normalisieren
        inv = pow(M[row][col], p - 2, p)
        M[row] = [(c * inv) % p for c in M[row]]

        # Eliminierung
        for r in range(n):
            if r != row and M[r][col] != 0:
                factor = M[r][col]
                M[r] = [(M[r][j] - factor * M[row][j]) % p for j in range(n)]

        row += 1

    # Freie Variablen (nicht-Pivot-Spalten)
    free_cols = [c for c in range(n) if c not in pivots]

    # Kernbasis: für jede freie Variable einen Basisvektor
    kernel = []

    # Konstanter Vektor (1-Vektor als trivialer Kern, wenn dim > 0)
    if not free_cols:
        return [[1] + [0] * (n - 1)]  # Nur triviale Lösung: nicht im Kern

    for fc in free_cols:
        vec = [0] * n
        vec[fc] = 1
        for i, pc in enumerate(pivots):
            # Lösung für Pivot-Variable aus freier Variable
            # Zeile i: x_{pc} = -M_rref[i][fc] * vec[fc]
            val = 0
            for j in free_cols:
                val = (val - M[i][j] * vec[j]) % p
            vec[pc] = val % p
        kernel.append(vec)

    # Füge den konstanten Vektor hinzu (immer im Kern, da Berlekamp 1 enthält)
    const_vec = [0] * n
    const_vec[0] = 1  # Konstante 1
    return [const_vec] + kernel


def _factorize_brute_force(poly: list[int], p: int) -> dict:
    """
    Brute-Force Faktorisierung für kleine Polynomgrade und Primzahlen.

    @param poly: Zu faktorisierendes Polynom
    @param p: Primzahl
    @return: {'factors': list, 'multiplicities': list}
    @lastModified: 2026-03-10
    """
    from itertools import product as iproduct
    poly = _poly_mod(poly, p)
    deg = _poly_degree(poly)
    factors = []
    mults = []
    remaining = poly

    for factor_deg in range(1, deg):
        for lower_coeffs in iproduct(range(p), repeat=factor_deg):
            candidate = _poly_mod([1] + list(lower_coeffs), p)
            mult = 0
            while True:
                q, r = _poly_divmod_mod(remaining, candidate, p)
                if _poly_degree(r) < 0 or (len(r) == 1 and r[0] == 0):
                    remaining = q
                    mult += 1
                else:
                    break
            if mult > 0:
                factors.append(candidate)
                mults.append(mult)

    # Verbleibender Teil (irreduzibel)
    if _poly_degree(remaining) > 0:
        factors.append(remaining)
        mults.append(1)

    if not factors:
        return {'factors': [poly], 'multiplicities': [1]}
    return {'factors': factors, 'multiplicities': mults}


# ===========================================================================
# FUNKTION: quotient_ring_field_check
# ===========================================================================

def quotient_ring_field_check(poly: list[int], p: int) -> dict:
    """
    Überprüft ob ℤ_p[x]/(f) ein Körper ist und bestimmt seine Eigenschaften.

    ℤ_p[x]/(f) ist Körper ⟺ f irreduzibel über ℤ_p.
    Die Ordnung des endlichen Körpers GF(p^n) ist p^n mit n = deg(f).

    @param poly: Koeffizienten von f(x) (höchster Grad zuerst)
    @param p: Primzahl
    @return: Dictionary mit Eigenschaften des Quotientenrings
    @lastModified: 2026-03-10
    """
    poly = _poly_mod(poly, p)
    deg = _poly_degree(poly)

    is_irred = is_irreducible_mod_p(poly, p)
    order = p ** deg if deg > 0 else p

    return {
        'polynomial': poly,
        'prime': p,
        'degree': deg,
        'is_irreducible': is_irred,
        'is_field': is_irred,
        'field_order': order if is_irred else None,
        'description': (
            f"GF({p}^{deg}) = GF({order})" if is_irred
            else f"ℤ_{p}[x]/({poly}) ist kein Körper (f nicht irreduzibel)"
        )
    }


# ===========================================================================
# FUNKTION: noetherian_check
# ===========================================================================

def noetherian_check(ideals: list) -> bool:
    """
    Prüft die noethersche Bedingung für eine Liste von Idealen.

    Ein Ring ist noethersch, wenn jede aufsteigende Kette von Idealen
    I₁ ⊆ I₂ ⊆ I₃ ⊆ ... nach endlich vielen Schritten stationär wird:
    es gibt ein N mit I_N = I_{N+1} = I_{N+2} = ...

    Testet empirisch für die gegebene Liste ob eine nicht-stationäre
    aufsteigende Kette vorliegt.

    @param ideals: Liste von Ideal-Objekten (aufsteigende Kette)
    @return: True wenn die Kette stationär wird (noethersche Bedingung erfüllt)
    @lastModified: 2026-03-10
    """
    if len(ideals) <= 1:
        return True

    # Prüfe ob es sich um eine aufsteigende Kette handelt
    # und ob sie stationär wird
    for i in range(len(ideals) - 1):
        I = ideals[i]
        J = ideals[i + 1]
        # I muss Teilmenge von J sein (aufsteigende Kette)
        if not I <= J:
            return True  # Keine aufsteigende Kette → noethersche Bedingung trivial

    # Prüfe ob die Kette stationär wird
    # Finde letztes echtes Inklusionsschritt
    stabilizes = False
    for i in range(len(ideals) - 1):
        I = ideals[i]
        J = ideals[i + 1]
        if set(I.elements) == set(J.elements):
            # Ab hier stationär
            stabilizes = True
            # Alle weiteren müssen gleich sein
            for k in range(i + 1, len(ideals)):
                if set(ideals[k].elements) != set(J.elements):
                    return False
            break

    # Eine endliche aufsteigende Kette wird immer stationär
    return True


# ===========================================================================
# FUNKTION: krull_dimension_estimate
# ===========================================================================

def krull_dimension_estimate(prime_ideals_chain: list) -> int:
    """
    Schätzt die Krull-Dimension anhand einer Primidealkette.

    Die Krull-Dimension eines Rings R ist die Länge der längsten
    strikten Kette von Primidealen:
        p₀ ⊊ p₁ ⊊ p₂ ⊊ ... ⊊ pₙ

    Beispiele:
        - Körper: dim = 0 (nur Primideal ist (0))
        - ℤ oder k[x]: dim = 1
        - k[x,y]: dim = 2

    @param prime_ideals_chain: Liste von Ideal-Objekten, die eine Primidealkette bilden
    @return: Länge der Kette (= Krull-Dimension wenn Kette maximal)
    @lastModified: 2026-03-10
    """
    if not prime_ideals_chain:
        return 0

    # Prüfe ob alle Ideale Primideale sind
    valid_chain_length = 0
    for i, ideal in enumerate(prime_ideals_chain):
        if not ideal.is_prime():
            break
        valid_chain_length += 1

    # Prüfe Striktheit der Inklusion
    strict_inclusions = 0
    for i in range(min(valid_chain_length - 1, len(prime_ideals_chain) - 1)):
        I = prime_ideals_chain[i]
        J = prime_ideals_chain[i + 1]
        if I <= J and set(I.elements) != set(J.elements):
            strict_inclusions += 1

    return strict_inclusions


# ===========================================================================
# FUNKTION: ring_of_integers
# ===========================================================================

def ring_of_integers(d: int) -> dict:
    """
    Bestimmt Eigenschaften des Rings der ganzen Zahlen ℤ[√d] (quadratischer Zahlkörper).

    ℤ[√d] = {a + b√d : a, b ∈ ℤ} mit der Norm N(a + b√d) = a² - d·b².

    Bekannte Fälle:
        d = -1: ℤ[i]  (Gauß'sche Zahlen, Klassenzahl 1, PID)
        d = -3: ℤ[ω]  (Eisenstein-Zahlen, Klassenzahl 1, PID)
        d = -2: ℤ[√-2] (Klassenzahl 1, PID)
        d =  2: ℤ[√2]  (PID, fundamentale Einheit 1 + √2)
        d =  3: ℤ[√3]  (PID, fundamentale Einheit 2 + √3)
        d =  5: ℤ[√5]  (nicht PID, Klassenzahl 2)

    Diskriminante:
        Δ = 4d     wenn d ≡ 2,3 (mod 4)
        Δ = d      wenn d ≡ 1   (mod 4)

    @param d: Quadratfreie ganze Zahl ≠ 0, 1 (Radikand)
    @return: Dictionary mit Ringeigenschaften
    @raises ValueError: Wenn d = 0 oder d ein vollständiges Quadrat ist
    @lastModified: 2026-03-10
    """
    if d == 0:
        raise ValueError("d = 0 ergibt keinen Zahlkörper.")
    if d == 1:
        raise ValueError("d = 1 ergibt ℤ (trivialer Fall).")

    # Prüfe ob d quadratfrei ist
    abs_d = abs(d)
    for prime in range(2, int(abs_d**0.5) + 1):
        if abs_d % (prime * prime) == 0:
            raise ValueError(f"d = {d} ist nicht quadratfrei (durch {prime}² teilbar).")

    # Diskriminante berechnen
    if d % 4 == 1:
        discriminant = d
        ring_of_integers_basis = f"{{1, (1+√{d})/2}}"
    else:
        discriminant = 4 * d
        ring_of_integers_basis = f"{{1, √{d}}}"

    # Bekannte Klassenzahlen und PID-Status (Tabelle für kleine |d|)
    # Negative d: imaginär-quadratische Zahlkörper
    # Positive d: reell-quadratische Zahlkörper
    pid_negative = {-1, -2, -3, -7, -11, -19, -43, -67, -163}
    # Reell-quadratische Zahlkörper ℤ[√d] mit Klassenzahl 1 (PIDs) für kleine d:
    # d=5 hat Klassenzahl 2 und ist KEIN PID. d=10 ebenfalls nicht.
    pid_positive_small = {2, 3, 6, 7, 11, 13, 14, 17, 19, 21}

    is_pid = (d in pid_negative) or (d > 0 and d in pid_positive_small)

    # Einheitengruppe
    if d < 0:
        # Imaginär-quadratisch: endliche Einheitengruppe
        if d == -1:
            units_desc = "{±1, ±i}"
            fundamental_unit = None
        elif d == -3:
            units_desc = "{±1, ±ω, ±ω²}  (6. Einheitswurzeln)"
            fundamental_unit = None
        else:
            units_desc = "{±1}"
            fundamental_unit = None
    else:
        # Reell-quadratisch: unendliche Einheitengruppe (Pell-Gleichung)
        units_desc = "{±(a + b√d)^n : n ∈ ℤ, a² - db² = 1}"
        # Fundamentale Einheit: kleinste Lösung der Pell-Gleichung x² - d·y² = 1
        fundamental_unit = _find_fundamental_unit(d)

    # Beschreibung des Rings
    if d == -1:
        ring_name = "ℤ[i] (Gauß'sche Zahlen)"
    elif d == -3:
        ring_name = "ℤ[ω] (Eisenstein-Zahlen, ω = e^{2πi/3})"
    else:
        ring_name = f"ℤ[√{d}]"

    return {
        'd': d,
        'ring_name': ring_name,
        'discriminant': discriminant,
        'ring_of_integers_basis': ring_of_integers_basis,
        'is_pid': is_pid,
        'units_description': units_desc,
        'fundamental_unit': fundamental_unit,
        'is_imaginary': d < 0,
        'norm_formula': f"N(a + b√{d}) = a² - {d}b²"
    }


def _find_fundamental_unit(d: int) -> Optional[tuple[int, int]]:
    """
    Findet die fundamentale Einheit in ℤ[√d] (kleinste Lösung der Pell-Gleichung).

    Löst x² - d·y² = ±1 durch Kettenbruchentwicklung von √d.

    @param d: Positives quadratfreies d
    @return: (a, b) mit a + b√d als fundamentale Einheit, oder None
    @lastModified: 2026-03-10
    """
    if d <= 0:
        return None

    # Kettenbruch-Algorithmus für √d
    sqrt_d = math.isqrt(d)
    if sqrt_d * sqrt_d == d:
        return None  # d ist vollständiges Quadrat

    # Suche Lösung von x² - d·y² = ±1 via Kettenbuch
    # Konvergenten des Kettenbruchs von √d
    m, d0, a0 = 0, 1, sqrt_d
    p_prev, p_curr = 1, a0
    q_prev, q_curr = 0, 1

    for _ in range(1000):
        m = d0 * a0 - m
        d0 = (d - m * m) // d0
        if d0 == 0:
            break
        a0 = (sqrt_d + m) // d0

        p_prev, p_curr = p_curr, a0 * p_curr + p_prev
        q_prev, q_curr = q_curr, a0 * q_curr + q_prev

        # Prüfe ob (p_curr, q_curr) Lösung von Pell ist
        if p_curr * p_curr - d * q_curr * q_curr in (1, -1):
            return (p_curr, q_curr)

    return None


# ===========================================================================
# FUNKTION: gaussian_integers_gcd
# ===========================================================================

def gaussian_integers_gcd(a: complex, b: complex) -> complex:
    """
    Berechnet den ggT zweier Gauß'scher Zahlen in ℤ[i] via euklidischem Algorithmus.

    ℤ[i] = {a + bi : a, b ∈ ℤ} ist ein euklidischer Bereich mit der Norm
    N(a + bi) = a² + b².

    Euklidische Division: Zu α, β ∈ ℤ[i], β ≠ 0, existieren γ, ρ ∈ ℤ[i] mit
        α = γ·β + ρ  und  N(ρ) < N(β)

    γ wird berechnet als gerundetes α/β ∈ ℚ[i].

    @param a: Erste Gauß'sche Zahl (als komplexe Zahl mit ganzzahligem Real- und Imaginärteil)
    @param b: Zweite Gauß'sche Zahl
    @return: Größter gemeinsamer Teiler (normiert: Real > 0 oder Real = 0, Imag > 0)
    @lastModified: 2026-03-10
    """
    def gauss_norm(z: complex) -> int:
        """Gauß'sche Norm N(a+bi) = a² + b²."""
        return round(z.real) ** 2 + round(z.imag) ** 2

    def gauss_round(z: complex) -> complex:
        """Rundet eine komplexe Zahl auf nächsten Gauß'schen Integer."""
        return complex(round(z.real), round(z.imag))

    # Euklidischer Algorithmus für Gauß'sche Zahlen
    a = complex(round(a.real), round(a.imag))
    b = complex(round(b.real), round(b.imag))

    while gauss_norm(b) > 0:
        # Berechne Quotient γ = round(a/b)
        gamma = gauss_round(a / b)
        # Rest ρ = a - γ·b
        a, b = b, a - gamma * b

    # Normalisierung: Einheiten in ℤ[i] sind {1, -1, i, -i}
    # Wähle Normalform mit positivem Realteil, oder wenn 0 dann positivem Imaginärteil
    result = complex(round(a.real), round(a.imag))

    # Normiere auf assoziiertes Element mit Re > 0 oder (Re = 0 und Im > 0)
    for unit in [1, -1, 1j, -1j]:
        candidate = result * unit
        cr = round(candidate.real)
        ci = round(candidate.imag)
        if cr > 0 or (cr == 0 and ci > 0):
            return complex(cr, ci)

    return result
