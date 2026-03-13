"""
@file beal_conjecture.py
@brief Beal Conjecture — Computational Verification and ABC-Connection
@description
    Beal's Conjecture (1993): If A^x + B^y = C^z where A, B, C, x, y, z are
    positive integers with x, y, z >= 3, then A, B, C have a common prime factor.

    Note: Fermat's Last Theorem is a special case (x=y=z=n ≥ 3, A,B,C paarweise
    teilerfremd würde einen Widerspruch zu FLT ergeben).

    Connection to ABC-conjecture: Wenn Beal für teilerfremde A,B,C falsch wäre,
    erhielte man ein ABC-Tripel mit hoher Qualität q = log(C^z)/log(rad(A^x*B^y*C^z)).
    Die ABC-Vermutung schränkt solche Qualitäten ein.

    Computational approach:
    - Precompute alle k^e für k ≤ N, e ∈ {3,...,max_exp} mit Primfaktorisierungen
    - Build lookup table: Wert → [(base, exp), ...]
    - Für jedes Paar (A^x, B^y): prüfe ob A^x + B^y = C^z für ein C,z
    - Verifiziere gcd(A,B,C) > 1 wann immer Gleichheit gilt

@author Michael Fuhrmann
@date 2026-03-13
@lastModified 2026-03-13
"""

import math
from typing import Dict, List, Optional, Tuple
from sympy import factorint, gcd, Integer


class BealTriple:
    """
    Repräsentiert ein Beal-Tripel der Form A^x + B^y = C^z.

    Frankl's NOTE: Dies ist ein CONJECTURE (Vermutung), kein bewiesener Satz.
    Alle Tripel mit x, y, z >= 3 sollten gcd(A, B, C) > 1 haben.
    """

    def __init__(self, A: int, x: int, B: int, y: int, C: int, z: int):
        """
        Initialisiert ein Beal-Tripel A^x + B^y = C^z.

        @param A  Basis der ersten Potenz
        @param x  Exponent der ersten Potenz (sollte >= 3 für Beal-Relevanz)
        @param B  Basis der zweiten Potenz
        @param y  Exponent der zweiten Potenz
        @param C  Basis der dritten Potenz
        @param z  Exponent der Summe
        @lastModified 2026-03-13
        """
        self.A = A
        self.x = x
        self.B = B
        self.y = y
        self.C = C
        self.z = z

    def is_valid(self) -> bool:
        """
        Prüft ob A^x + B^y = C^z arithmetisch korrekt ist.

        @return True falls die Gleichung stimmt
        @lastModified 2026-03-13
        """
        try:
            lhs = self.A ** self.x + self.B ** self.y
            rhs = self.C ** self.z
            return lhs == rhs
        except (OverflowError, ValueError):
            return False

    def has_common_factor(self) -> bool:
        """
        Prüft ob gcd(A, B, C) > 1, d.h. A, B, C einen gemeinsamen Primteiler haben.

        @return True falls gemeinsamer Faktor existiert
        @lastModified 2026-03-13
        """
        g = int(gcd(gcd(Integer(self.A), Integer(self.B)), Integer(self.C)))
        return g > 1

    def satisfies_beal(self) -> bool:
        """
        Prüft ob dieses Tripel Beal's CONJECTURE erfüllt.

        Beal's CONJECTURE (kein bewiesener Satz!): Falls x, y, z >= 3,
        dann muss gcd(A, B, C) > 1 gelten.

        @return True falls Beal erfüllt (oder nicht anwendbar wegen kleiner Exponenten)
        @lastModified 2026-03-13
        """
        # Beal gilt nur für x, y, z >= 3
        if self.x < 3 or self.y < 3 or self.z < 3:
            return True  # Nicht im Anwendungsbereich
        return self.has_common_factor()

    def radical(self) -> int:
        """
        Berechnet rad(A · B · C) = Produkt aller verschiedenen Primteiler von A·B·C.

        Das Radikal ist ein zentrales Konzept in der ABC-Vermutung.

        @return Radikal rad(A·B·C)
        @lastModified 2026-03-13
        """
        # Alle Primteiler von A, B, C sammeln (ohne Vielfachheit)
        primes: set = set()
        for val in [self.A, self.B, self.C]:
            if val > 1:
                primes.update(factorint(val).keys())
        result = 1
        for p in primes:
            result *= p
        return result

    def common_factor(self) -> int:
        """
        Gibt gcd(A, B, C) zurück.

        @return Größter gemeinsamer Teiler der drei Basen
        @lastModified 2026-03-13
        """
        return int(gcd(gcd(Integer(self.A), Integer(self.B)), Integer(self.C)))

    def __repr__(self) -> str:
        """String-Darstellung des Tripels."""
        return (
            f"BealTriple({self.A}^{self.x} + {self.B}^{self.y} = "
            f"{self.C}^{self.z}, "
            f"valid={self.is_valid()}, gcd={self.common_factor()})"
        )


def _radical_of_product(values: List[int]) -> int:
    """
    Hilfsfunktion: Berechnet das Radikal des Produkts einer Liste von Zahlen.

    rad(a1 * a2 * ... * ak) = Produkt aller verschiedenen Primteiler.

    @param values  Liste positiver ganzer Zahlen
    @return        Radikal des Produkts
    @lastModified 2026-03-13
    """
    primes: set = set()
    for v in values:
        if v > 1:
            primes.update(factorint(int(v)).keys())
    result = 1
    for p in primes:
        result *= p
    return result


class BealChecker:
    """
    Verifiziert Beal's CONJECTURE für Basen A, B, C bis zur Schranke N.

    Verwendet Vorberechnungs-Lookup-Tabellen für effiziente Suche.
    """

    def __init__(self, N: int = 100, max_exp: int = 10):
        """
        Initialisiert den Checker mit Suchschranken.

        @param N        Obere Schranke für Basen A, B, C ≤ N
        @param max_exp  Maximaler Exponent (mindestens 3)
        @lastModified 2026-03-13
        """
        self.N = N
        self.max_exp = max(3, max_exp)
        # Wird lazy berechnet
        self._powers: Optional[Dict[int, List[Tuple[int, int]]]] = None

    def precompute_powers(self) -> Dict[int, List[Tuple[int, int]]]:
        """
        Berechnet alle k^e für k ≤ N, e ∈ {3, ..., max_exp}.

        Gibt ein Dict zurück: Wert → [(base, exp), ...]
        Mehrere Darstellungen desselben Wertes (z.B. 64 = 2^6 = 4^3 = 8^2)
        werden alle gespeichert.

        @return Dict von Werten auf ihre (Basis, Exponent)-Darstellungen
        @lastModified 2026-03-13
        """
        powers: Dict[int, List[Tuple[int, int]]] = {}
        for base in range(1, self.N + 1):
            for exp in range(3, self.max_exp + 1):
                value = base ** exp
                if value not in powers:
                    powers[value] = []
                powers[value].append((base, exp))
        self._powers = powers
        return powers

    def _get_powers(self) -> Dict[int, List[Tuple[int, int]]]:
        """Gibt die Lookup-Tabelle zurück (berechnet sie falls nötig)."""
        if self._powers is None:
            self.precompute_powers()
        return self._powers  # type: ignore

    def find_all_triples(self) -> List[BealTriple]:
        """
        Sucht alle Beal-Tripel A^x + B^y = C^z mit A, B, C ≤ N, x, y, z ≥ 3.

        Methode:
        1. Für jedes Paar (A^x, B^y) berechne die Summe S = A^x + B^y
        2. Prüfe ob S in der Lookup-Tabelle als C^z vorhanden ist
        3. Füge gültige Tripel hinzu

        @return Liste aller gefundenen Beal-Tripel
        @lastModified 2026-03-13
        """
        powers = self._get_powers()
        triples = []
        power_values = list(powers.keys())

        for val_ax in power_values:
            for val_by in power_values:
                target = val_ax + val_by
                if target in powers:
                    # Tripel gefunden: alle Kombinationen von Darstellungen
                    for (A, x) in powers[val_ax]:
                        for (B, y) in powers[val_by]:
                            for (C, z) in powers[target]:
                                # Nur physikalisch gültige Tripel (alle Basen ≤ N)
                                if A <= self.N and B <= self.N and C <= self.N:
                                    triple = BealTriple(A, x, B, y, C, z)
                                    if triple.is_valid():
                                        triples.append(triple)
        return triples

    def find_counterexamples(self) -> List[BealTriple]:
        """
        Sucht nach Gegenbeispielen zu Beal's CONJECTURE.

        Ein Gegenbeispiel wäre: A^x + B^y = C^z mit x, y, z ≥ 3 und gcd(A,B,C) = 1.

        Erwartetes Ergebnis (gemäß CONJECTURE): leere Liste.

        @return Liste der Gegenbeispiele (erwartet: leer)
        @lastModified 2026-03-13
        """
        triples = self.find_all_triples()
        counterexamples = []
        for triple in triples:
            # x, y, z >= 3 ist schon durch precompute_powers garantiert
            if not triple.has_common_factor():
                counterexamples.append(triple)
        return counterexamples

    @staticmethod
    def known_families() -> List[BealTriple]:
        """
        Bekannte Familien gültiger Beal-Tripel (alle mit gemeinsamem Faktor).

        Beispiele:
        - 2^3 + 2^3 = 2^4:      8 + 8 = 16, gcd(2,2,2) = 2 ✓
        - 2^6 + 2^6 = 2^7:      64 + 64 = 128, gcd(2,2,2) = 2 ✓
        - 3^3 + 6^3 = 3^5:      27 + 216 = 243, gcd(3,6,3) = 3 ✓
        - 6^3 + 10^3 = 2^3*...  (komplexere Beispiele)
        - 14^3 + 17^3 = 31^... (nur wenn Gleichung stimmt)

        @return Liste bekannter gültiger Beal-Tripel
        @lastModified 2026-03-13
        """
        candidates = [
            BealTriple(2, 3, 2, 3, 2, 4),    # 8 + 8 = 16 = 2^4
            BealTriple(2, 6, 2, 6, 2, 7),    # 64 + 64 = 128 = 2^7
            BealTriple(3, 3, 6, 3, 3, 5),    # 27 + 216 = 243 = 3^5
            BealTriple(6, 3, 10, 3, 2, 3),   # würde geprüft (invalid falls falsch)
            BealTriple(2, 9, 8, 3, 4, 5),    # 512 + 512 = 1024 = 4^5? prüfen
            BealTriple(9, 3, 3, 6, 6, 4),    # 729 + 729 = 1458? prüfen
        ]
        # Nur arithmetisch gültige Tripel zurückgeben
        return [t for t in candidates if t.is_valid()]


class ABCConnection:
    """
    Verbindung zwischen Beal's CONJECTURE und der ABC-Vermutung.

    Die ABC-Vermutung (Oesterlé-Masser 1985) besagt:
    Für jedes ε > 0 existiert C_ε > 0, sodass für alle teilerfremden a + b = c
    gilt: c < C_ε * rad(a·b·c)^(1+ε).

    Equivalently: q(a,b,c) = log(c) / log(rad(abc)) < 1 + ε für fast alle Tripel.
    """

    @staticmethod
    def radical(n: int) -> int:
        """
        Berechnet das Radikal rad(n) = Produkt aller verschiedenen Primteiler von n.

        @param n  Positive ganze Zahl
        @return   rad(n)
        @lastModified 2026-03-13
        """
        if n <= 1:
            return n
        factors = factorint(int(n))
        result = 1
        for p in factors.keys():
            result *= p
        return result

    @staticmethod
    def quality(a: int, b: int, c: int) -> float:
        """
        Berechnet die ABC-Qualität q(a, b, c) = log(c) / log(rad(a·b·c)).

        Oesterlé-Masser-CONJECTURE: q < 1 + ε für fast alle teilerfremden Tripel
        mit a + b = c. Qualitäten > 1 sind sehr selten (höchste bekannte: ~1.6299).

        @param a  Erster Summand (a + b = c, a,b,c teilerfremd)
        @param b  Zweiter Summand
        @param c  Summe c = a + b
        @return   ABC-Qualität (0.0 falls c ≤ 1)
        @lastModified 2026-03-13
        """
        if c <= 1:
            return 0.0
        rad = _radical_of_product([a, b, c])
        if rad <= 1:
            return 0.0
        return math.log(c) / math.log(rad)

    @staticmethod
    def beal_abc_bound(A: int, x: int, B: int, y: int, C: int, z: int) -> float:
        """
        Berechnet die ABC-Qualität für ein Beal-Tripel A^x + B^y = C^z.

        Falls A, B, C teilerfremd wären (Gegenbeispiel zu Beal), würde gelten:
        q = z·log(C) / log(rad(A^x · B^y · C^z)) = z·log(C) / log(rad(A·B·C))

        ABC-CONJECTURE → q < 1 + ε → dies schränkt mögliche Exponenten ein.

        Da rad(A^x * B^y * C^z) = rad(A * B * C), vereinfacht sich die Formel.

        @param A  Basis der ersten Potenz
        @param x  Exponent der ersten Potenz
        @param B  Basis der zweiten Potenz
        @param y  Exponent der zweiten Potenz
        @param C  Basis der dritten Potenz (Summe)
        @param z  Exponent der Summe
        @return   ABC-Qualität des entsprechenden Tripels
        @lastModified 2026-03-13
        """
        if C <= 1 or z < 1:
            return 0.0
        # Radikal von A^x * B^y * C^z = rad(A) * rad(B) * rad(C) (gemeinsame Primes einmal)
        rad = _radical_of_product([A, B, C])
        if rad <= 1:
            return 0.0
        # ABC-Qualität: c = C^z, rad(abc) = rad(A·B·C)
        c_val = C ** z
        if c_val <= 1:
            return 0.0
        return math.log(c_val) / math.log(rad)

    @staticmethod
    def high_quality_abc_triples(limit: int = 1000) -> List[Tuple[int, int, int, float]]:
        """
        Sucht ABC-Tripel a + b = c mit hoher Qualität q > 1.

        Diese Tripel sind mathematisch bemerkenswert und selten.

        @param limit  Obere Schranke für c
        @return       Liste von (a, b, c, quality) mit q > 1, sortiert nach Qualität
        @lastModified 2026-03-13
        """
        results = []
        for c in range(3, limit + 1):
            for a in range(1, c):
                b = c - a
                if b <= 0 or a <= 0:
                    continue
                # Nur teilerfremde Tripel (gcd(a,b) = gcd(a,c) = 1)
                if int(gcd(Integer(a), Integer(b))) != 1:
                    continue
                q = ABCConnection.quality(a, b, c)
                if q > 1.0:
                    results.append((a, b, c, q))
        # Nach Qualität absteigend sortieren
        results.sort(key=lambda t: t[3], reverse=True)
        return results
