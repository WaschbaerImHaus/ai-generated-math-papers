"""
@file coding_theory.py
@brief Kodierungstheorie und Formale Sprachen – Galois-Felder, lineare Codes,
       zyklische Codes, BCH-Codes, Automaten, kontextfreie Grammatiken.
@description
    Dieses Modul implementiert die Kernkonzepte der algebraischen Kodierungstheorie
    sowie der Theorie formaler Sprachen und Automaten.

    Enthaltene Klassen:
    - GaloisField           : Galois-Körper GF(p^n), Arithmetik, primitives Element
    - LinearCode            : Lineare [n,k,d]-Codes über GF(2)
    - HammingCode           : Hamming-Codes [2^r-1, 2^r-r-1, 3]
    - CyclicCode            : Zyklische Codes mit Generatorpolynom
    - BCHCode               : BCH- und Reed-Solomon-Codes (Demo), Berlekamp-Massey
    - FormalLanguages       : DFA, NFA, reguläre Ausdrücke, CFG, Kellerautomat
    - ComputabilityTheory   : Berechenbarkeitstheorie (ergänzend zu recursion_theory.py)

    Mathematische Grundlagen:
        Linearer Code [n,k,d]:
            n = Codewortlänge, k = Informationsbits, d = Minimaldistanz
            Codewort c = m·G  (m: Nachricht, G: Generatormatrix k×n)
            Prüfmatrix H: H·cᵀ = 0 für alle Codewörter c
            Syndrom s = H·rᵀ (r: empfangenes Wort; s = 0 → fehlerfrei)

        Hamming-Code [2^r-1, 2^r-r-1, 3]:
            Kann 1 Fehler korrigieren (perfekter Code)
            Prüfmatrix H: Spalten = alle 2^r-1 Bitvektoren ≠ 0

        Zyklischer Code:
            Jede zyklische Verschiebung eines Codewortes ist wieder ein Codewort
            Beschrieben durch Generatorpolynom g(x) | (x^n - 1)

        Galois-Körper GF(p^n):
            Eindeutiger Körper mit p^n Elementen (p Primzahl)
            Multiplikative Gruppe ist zyklisch (primitives Element g: g^(p^n-1) = 1)

    HINWEIS: Hamming-Distanz, Singleton-Bound, (7,4)-Hamming-Code und Paritätscode
    sind bereits in information_theory.py (Klasse ErrorCorrection) implementiert.
    Dieses Modul ergänzt mit algebraischer Kodierungstheorie und formalen Sprachen.

@author Kurt Ingwer
@version 1.0
@date 2026-03-10
@lastModified 2026-03-10
"""

from __future__ import annotations

import re as _re
from itertools import product as _product
from typing import Any, Callable

import numpy as np


# ============================================================
# 1. GaloisField – Galois-Körper GF(p^n)
# ============================================================

class GaloisField:
    """
    @class GaloisField
    @brief Galois-Körper GF(p^n) mit p Primzahl und n ≥ 1.
    @description
        Für n=1: GF(p) = Z/pZ – Restklassenkörper mod p.
        Arithmetik erfolgt direkt mod p (Addition und Multiplikation).

        Für n>1 (GF(p^n)) werden die n ersten Elemente {0,1,...,p^n-1} als
        Repräsentanten der Polynomreste über GF(p)[x]/(f(x)) betrachtet.
        Die vollständige Arithmetik über Erweiterungskörper ist als Demo
        implementiert; komplette Erweiterungskörper-Multiplikation ist für
        n=1 exakt, für n>1 als konzeptioneller Hinweis gedacht.

        Formeln:
            |GF(p^n)| = p^n
            Multiplikative Gruppe: GF(p^n)* ist zyklisch der Ordnung p^n - 1
            Primitives Element g: ord(g) = p^n - 1

    @param p  Primzahl (Charakteristik des Körpers)
    @param n  Erweiterungsgrad (n=1 → Primkörper GF(p))

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, p: int, n: int = 1) -> None:
        """
        @brief Initialisiert den Galois-Körper GF(p^n).
        @param p  Primzahl
        @param n  Erweiterungsgrad (Standard: 1)
        @raises ValueError  Wenn p keine Primzahl oder n < 1
        @lastModified 2026-03-10
        """
        if not self._is_prime(p):
            raise ValueError(f"p={p} ist keine Primzahl.")
        if n < 1:
            raise ValueError(f"Erweiterungsgrad n muss ≥ 1 sein, ist aber {n}.")
        self.p = p          # Charakteristik
        self.n = n          # Erweiterungsgrad
        self.order = p ** n  # Körpergröße |GF(p^n)|

    # ----------------------------------------------------------
    # Hilfsmethoden
    # ----------------------------------------------------------

    @staticmethod
    def _is_prime(num: int) -> bool:
        """
        @brief Prüft, ob eine Zahl eine Primzahl ist.
        @param num  Zu prüfende Zahl
        @return True wenn Primzahl, sonst False
        @lastModified 2026-03-10
        """
        if num < 2:
            return False
        if num == 2:
            return True
        if num % 2 == 0:
            return False
        for i in range(3, int(num ** 0.5) + 1, 2):
            if num % i == 0:
                return False
        return True

    # ----------------------------------------------------------
    # Öffentliche Methoden
    # ----------------------------------------------------------

    def elements(self) -> list[int]:
        """
        @brief Gibt alle Elemente des Körpers GF(p^n) zurück.
        @description
            Für GF(p): {0, 1, ..., p-1}
            Für GF(p^n): {0, 1, ..., p^n-1} als Repräsentanten
        @return Liste aller Körperelemente
        @lastModified 2026-03-10
        """
        return list(range(self.order))

    def add(self, a: int, b: int) -> int:
        """
        @brief Addition in GF(p): (a + b) mod p.
        @description
            In GF(p^n) ist die Addition komponentenweise (Polynomaddition mod p).
            Für n=1 entspricht dies der einfachen Addition mod p.
        @param a  Erstes Element
        @param b  Zweites Element
        @return Summe in GF(p)
        @raises ValueError  Wenn a oder b nicht in {0,...,p-1}
        @lastModified 2026-03-10
        """
        if not (0 <= a < self.p and 0 <= b < self.p):
            raise ValueError(f"Elemente müssen in {{0,...,{self.p-1}}} liegen.")
        return (a + b) % self.p

    def multiply(self, a: int, b: int) -> int:
        """
        @brief Multiplikation in GF(p): (a * b) mod p.
        @description
            Für n=1: direkte Multiplikation mod p.
            Für n>1 wäre Polynommultiplikation mod irreduzibles Polynom nötig
            (hier: konzeptioneller Hinweis, Primkörper-Arithmetik).
        @param a  Erstes Element
        @param b  Zweites Element
        @return Produkt in GF(p)
        @raises ValueError  Wenn a oder b nicht in {0,...,p-1}
        @lastModified 2026-03-10
        """
        if not (0 <= a < self.p and 0 <= b < self.p):
            raise ValueError(f"Elemente müssen in {{0,...,{self.p-1}}} liegen.")
        return (a * b) % self.p

    def primitive_element(self) -> int:
        """
        @brief Findet ein primitives Element (Generator) der multiplikativen Gruppe GF(p)*.
        @description
            Ein Element g ∈ GF(p)* heißt primitiv, wenn:
                ord(g) = p - 1  (d.h. g^(p-1) = 1 und g^k ≠ 1 für 0 < k < p-1)
            Die multiplikative Gruppe GF(p)* ist zyklisch, ein primitives Element
            existiert immer.

            Algorithmus:
            1. Bestimme Primfaktoren von p-1
            2. Prüfe für jedes g ∈ {2,...,p-1}: g^((p-1)/q) ≠ 1 für alle Primfaktoren q
        @return Kleinstes primitives Element ≥ 2 (oder -1 falls GF(1))
        @lastModified 2026-03-10
        """
        if self.p == 2:
            # GF(2)* = {1}, einziges primitives Element ist 1
            return 1

        # Primfaktoren von p-1 bestimmen
        phi = self.p - 1
        factors = self._prime_factors(phi)

        # Jedes Element als Generator prüfen
        for g in range(2, self.p):
            is_primitive = True
            for q in factors:
                # Wenn g^(phi/q) == 1, ist die Ordnung von g ein Teiler von phi/q < phi
                if pow(g, phi // q, self.p) == 1:
                    is_primitive = False
                    break
            if is_primitive:
                return g
        return -1  # Sollte für Primzahlen p nie eintreten

    @staticmethod
    def _prime_factors(n: int) -> list[int]:
        """
        @brief Berechnet die Primfaktoren (ohne Vielfachheit) von n.
        @param n  Zu faktorisierende Zahl
        @return Sortierte Liste eindeutiger Primfaktoren
        @lastModified 2026-03-10
        """
        factors = set()
        d = 2
        while d * d <= n:
            while n % d == 0:
                factors.add(d)
                n //= d
            d += 1
        if n > 1:
            factors.add(n)
        return sorted(factors)

    def discrete_log_gf(self, base: int, element: int) -> int:
        """
        @brief Diskreter Logarithmus in GF(p): log_base(element) mod p.
        @description
            Findet x mit base^x ≡ element (mod p).
            Algorithmus: Baby-Step Giant-Step (BSGS), O(√p) Zeitkomplexität.

            Sei m = ⌈√(p-1)⌉:
            Schreibe x = a·m - b mit a ∈ {1,...,m}, b ∈ {0,...,m-1}
            baby steps: element · base^b für b = 0,...,m-1
            giant steps: (base^m)^a für a = 1,...,m
            Kollision: base^(am) = element · base^b → x = am - b
        @param base     Basis (primitives Element)
        @param element  Gesuchtes Element
        @return Diskreter Logarithmus x, oder -1 falls nicht existiert
        @lastModified 2026-03-10
        """
        if element == 0:
            return -1  # log(0) ist undefiniert
        if element == 1:
            return 0   # base^0 = 1

        p = self.p
        m = int((p - 1) ** 0.5) + 1

        # Baby steps: Tabelle {element * base^b : b = 0,...,m-1}
        baby_steps: dict[int, int] = {}
        val = element % p
        for b in range(m):
            baby_steps[val] = b
            val = (val * base) % p

        # Giant steps: base^(m·a) für a = 1,...,m
        base_m = pow(base, m, p)
        giant = base_m
        for a in range(1, m + 1):
            if giant in baby_steps:
                b = baby_steps[giant]
                result = (a * m - b) % (p - 1)
                return result
            giant = (giant * base_m) % p

        return -1  # Kein Logarithmus gefunden

    def is_primitive_poly(self, coeffs: list[int]) -> bool:
        """
        @brief Prüft, ob ein Polynom über GF(2) primitiv ist.
        @description
            Ein Polynom f(x) ∈ GF(2)[x] vom Grad r heißt primitiv, wenn:
            1. f(x) irreduzibel über GF(2)
            2. Die Ordnung von x in GF(2)[x]/(f(x)) gleich 2^r - 1 ist
               (d.h. x ist primitives Element von GF(2^r)*)

            Anwendung: BCH-Codes, CRC-Polynome.

            Algorithmus:
            1. Prüfe Irreduzibilität (kein Faktor in GF(2)[x])
            2. Prüfe x^(2^r-1) ≡ 1 mod f(x) und x^d ≢ 1 für Teiler d < 2^r-1

        @param coeffs  Koeffizienten [a0, a1, ..., ar] mit ar = Leitkoeffizient
                       Beispiel: x^4 + x + 1 → [1, 1, 0, 0, 1]
        @return True wenn primitives Polynom über GF(2), sonst False
        @lastModified 2026-03-10
        """
        # Nur für GF(2) implementiert
        if self.p != 2:
            raise NotImplementedError("is_primitive_poly nur für GF(2) implementiert.")

        # Grad des Polynoms
        r = len(coeffs) - 1
        if r < 1:
            return False

        # Leitkoeffizient muss 1 sein (normiertes Polynom)
        c = [int(x) % 2 for x in coeffs]
        if c[-1] != 1:
            return False

        # Hilfsfunktion: Polynom-Multiplikation mod f über GF(2)
        def poly_mod(p_poly: list[int], mod_poly: list[int]) -> list[int]:
            """Reduziert p_poly modulo mod_poly über GF(2)."""
            result = list(p_poly)
            deg_mod = len(mod_poly) - 1
            while len(result) > deg_mod:
                # Führenden Term eliminieren
                coef = result[-1]
                if coef:
                    for i, m in enumerate(mod_poly):
                        result[len(result) - deg_mod - 1 + i] ^= m
                result.pop()
            return result

        def poly_mul_mod(a: list[int], b: list[int], mod_poly: list[int]) -> list[int]:
            """Multipliziert a und b, dann mod mod_poly über GF(2)."""
            # Produktpolynom
            prod = [0] * (len(a) + len(b) - 1)
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    prod[i + j] ^= (ai * bj) % 2
            return poly_mod(prod, mod_poly)

        # Irreduzibilität prüfen: kein Teilerpolynom für Grad 1..r//2
        for deg in range(1, r // 2 + 1):
            # Alle Polynome vom Grad deg über GF(2)
            for bits in range(1 << deg, 1 << (deg + 1)):
                factor = [(bits >> i) & 1 for i in range(deg + 1)]
                # Division c durch factor prüfen
                remainder = poly_mod(list(c), factor)
                # rest ist 0 wenn factor teilt c
                if all(x == 0 for x in remainder) and sum(factor) > 0:
                    # factor(0) = 1 muss sein (kein Nullteiler)
                    if factor[0] == 1:
                        return False

        # Primitivitätsprüfung: x^(2^r - 1) ≡ 1 mod f
        order = (1 << r) - 1
        factors_ord = self._prime_factors(order)

        # Starte mit x (als Polynom [0, 1] = x^1)
        x_poly = [0] * (r + 1)
        x_poly[1] = 1

        # x^order mod f prüfen (schnelle Exponentiation)
        def poly_pow_mod(base_p: list[int], exp: int, mod_p: list[int]) -> list[int]:
            """Schnelle Polynom-Potenzierung: base_p^exp mod mod_p über GF(2)."""
            result = [0] * len(mod_p)
            result[0] = 1  # 1-Polynom
            b = list(base_p)
            while exp > 0:
                if exp & 1:
                    result = poly_mul_mod(result, b, mod_p)
                b = poly_mul_mod(b, b, mod_p)
                exp >>= 1
            return result

        # x^(2^r-1) ≡ 1 (mod f)?
        xn = poly_pow_mod(x_poly, order, c)
        if xn[0] != 1 or any(xn[i] != 0 for i in range(1, len(xn))):
            return False

        # Für jeden Primfaktor q von 2^r-1: x^(order/q) ≢ 1 (mod f)
        for q in factors_ord:
            xd = poly_pow_mod(x_poly, order // q, c)
            if xd[0] == 1 and all(xd[i] == 0 for i in range(1, len(xd))):
                return False  # Ordnung von x teilt order/q → nicht primitiv

        return True


# ============================================================
# 2. LinearCode – Linearer Code [n, k, d] über GF(2)
# ============================================================

class LinearCode:
    """
    @class LinearCode
    @brief Linearer [n,k,d]-Code über GF(2) via Generatormatrix.
    @description
        Ein linearer Code C ist ein k-dimensionaler Untervektorraum von GF(2)^n.
        Jedes Codewort entsteht als c = m·G (m ∈ GF(2)^k, G: k×n-Matrix).

        Formeln:
            Codierung: c = m·G (Matrixprodukt über GF(2))
            Prüfmatrix: H·cᵀ = 0 für alle c ∈ C
            Syndrom: s = H·rᵀ (r: empfangener Vektor)
            Minimaldistanz: d = min{ wt(c) : c ∈ C, c ≠ 0 }

        Fehlerkorrektur:
            Ein Code mit Minimaldistanz d kann ⌊(d-1)/2⌋ Fehler korrigieren.
            Syndrom-Dekodierung: Syndrom s → Fehlermuster e → korrigiertes Wort r-e

    @param generator_matrix  Generatormatrix G als 2D-Liste oder numpy-Array (k×n über GF(2))

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, generator_matrix: list[list[int]] | np.ndarray) -> None:
        """
        @brief Initialisiert den linearen Code mit der Generatormatrix G.
        @param generator_matrix  k×n-Matrix über GF(2)
        @raises ValueError  Wenn Matrix-Einträge nicht in {0,1}
        @lastModified 2026-03-10
        """
        # Generatormatrix als numpy-Array speichern (mod 2 zur Sicherheit)
        self.G = np.array(generator_matrix, dtype=int) % 2
        self.k, self.n = self.G.shape  # k Informationsbits, n Codewortlänge

        # Grundlegende Validierung
        if not np.all((self.G == 0) | (self.G == 1)):
            raise ValueError("Generatormatrix muss nur Einträge aus {0,1} enthalten.")

    def encode(self, message: list[int] | np.ndarray) -> np.ndarray:
        """
        @brief Kodiert eine Nachricht m zu einem Codewort c = m·G (mod 2).
        @description
            Vektormultiplikation: c = m·G über GF(2)
            m ∈ GF(2)^k, G ∈ GF(2)^(k×n) → c ∈ GF(2)^n
        @param message  Nachrichtenvektor der Länge k (Einträge in {0,1})
        @return Codewort der Länge n
        @raises ValueError  Wenn Länge der Nachricht ≠ k
        @lastModified 2026-03-10
        """
        m = np.array(message, dtype=int) % 2
        if len(m) != self.k:
            raise ValueError(
                f"Nachrichtenlänge {len(m)} stimmt nicht mit k={self.k} überein."
            )
        # Matrixprodukt mod 2
        codeword = (m @ self.G) % 2
        return codeword

    def parity_check_matrix(self) -> np.ndarray:
        """
        @brief Berechnet die Prüfmatrix H mit H·cᵀ = 0 für alle c ∈ C.
        @description
            Für einen systematischen Code G = [Ik | P]:
                H = [-Pᵀ | In-k]  (in GF(2): H = [Pᵀ | In-k])

            Algorithmus (allgemein):
            1. G in systematische Form [Ik | P] bringen (Gauss-Elimination mod 2)
            2. H = [Pᵀ | In-k] konstruieren

            Die Prüfmatrix hat Dimensionen (n-k) × n.
        @return Prüfmatrix H der Dimension (n-k) × n
        @lastModified 2026-03-10
        """
        # G in Zeilenstufenform bringen (systematische Form)
        G_sys, pivot_cols = self._row_reduce_gf2(self.G.copy())

        # Redundanzspalten (nicht-Pivot-Spalten) bestimmen
        all_cols = list(range(self.n))
        non_pivot_cols = [c for c in all_cols if c not in pivot_cols]

        r = len(pivot_cols)  # Rang von G (= k bei vollem Rang)
        m_dim = self.n - r   # Anzahl Prüfbits

        # Prüfmatrix H aufbauen
        H = np.zeros((m_dim, self.n), dtype=int)

        # Für jede Nicht-Pivot-Spalte eine Zeile in H
        for row_idx, nc in enumerate(non_pivot_cols):
            H[row_idx, nc] = 1  # Einheitsmatrix-Teil
            # P-Teil: G_sys[:r, nc] enthält die Abhängigkeiten
            for pi, pc in enumerate(pivot_cols):
                if pi < r and nc < self.n:
                    H[row_idx, pc] = G_sys[pi, nc]

        return H

    @staticmethod
    def _row_reduce_gf2(
        matrix: np.ndarray,
    ) -> tuple[np.ndarray, list[int]]:
        """
        @brief Gauss-Elimination über GF(2) (Zeilenstufenform).
        @param matrix  Zu reduzierende Matrix (wird modifiziert)
        @return (reduzierte Matrix, Liste der Pivot-Spalten)
        @lastModified 2026-03-10
        """
        m_rows, n_cols = matrix.shape
        pivot_cols: list[int] = []
        row = 0

        for col in range(n_cols):
            if row >= m_rows:
                break
            # Pivot-Zeile suchen
            pivot = -1
            for r in range(row, m_rows):
                if matrix[r, col] == 1:
                    pivot = r
                    break
            if pivot == -1:
                continue

            # Zeilen tauschen
            matrix[[row, pivot]] = matrix[[pivot, row]]
            pivot_cols.append(col)

            # Eliminierung: alle anderen Zeilen mit Einsen in dieser Spalte
            for r in range(m_rows):
                if r != row and matrix[r, col] == 1:
                    matrix[r] = (matrix[r] + matrix[row]) % 2
            row += 1

        return matrix, pivot_cols

    def syndrome(self, received: list[int] | np.ndarray) -> np.ndarray:
        """
        @brief Berechnet das Syndrom s = H·rᵀ eines empfangenen Vektors r.
        @description
            s = H·rᵀ (mod 2)
            Falls s = 0: kein (erkannter) Fehler
            Falls s ≠ 0: Fehler vorhanden, s zeigt auf Fehlerspalte

        @param received  Empfangener Vektor der Länge n
        @return Syndromvektor der Länge n-k
        @raises ValueError  Wenn Länge des Vektors ≠ n
        @lastModified 2026-03-10
        """
        r = np.array(received, dtype=int) % 2
        if len(r) != self.n:
            raise ValueError(
                f"Empfangene Länge {len(r)} stimmt nicht mit n={self.n} überein."
            )
        H = self.parity_check_matrix()
        return (H @ r) % 2

    def minimum_distance(self) -> int:
        """
        @brief Berechnet die Minimaldistanz d des Codes (Brute-Force).
        @description
            d = min{ wt(c) : c ∈ C, c ≠ 0 }
            wt(c) = Hamming-Gewicht (Anzahl Einsen)

            Für kleine k iteriert diese Methode über alle 2^k - 1 Codewörter.
            Achtung: Für k > 20 nicht praktikabel (exponentiell).

        @return Minimaldistanz d (oder ∞ falls Code nur aus 0 besteht)
        @lastModified 2026-03-10
        """
        min_dist = float('inf')

        # Alle 2^k Nachrichtenvektoren (außer dem Nullvektor)
        for i in range(1, 1 << self.k):
            # Bit-Darstellung des Nachrichtenvektors
            m = np.array([(i >> j) & 1 for j in range(self.k)], dtype=int)
            c = self.encode(m)
            weight = int(np.sum(c))
            if 0 < weight < min_dist:
                min_dist = weight

        return int(min_dist) if min_dist != float('inf') else 0

    def error_correct(self, received: list[int] | np.ndarray) -> np.ndarray:
        """
        @brief Fehlerkorrektur via Syndrom-Tabelle (für 1-Fehler-Korrektur).
        @description
            Algorithmus:
            1. Berechne Syndrom s = H·rᵀ
            2. Falls s = 0: kein Fehler, gib r zurück
            3. Suche Spalte j in H mit H[:,j] = s → Fehlerposition j
            4. Korrigiere: c = r ⊕ e_j (Einheitsvektor in Position j)

            Funktioniert zuverlässig für Codes mit d ≥ 3 (1 Fehler korrigierbar).

        @param received  Empfangener Vektor der Länge n
        @return Korrigiertes Codewort
        @lastModified 2026-03-10
        """
        r = np.array(received, dtype=int) % 2
        s = self.syndrome(r)

        # Kein Fehler (oder Fehler nicht erkennbar)
        if np.all(s == 0):
            return r.copy()

        # Syndrom mit Spalten der Prüfmatrix vergleichen
        H = self.parity_check_matrix()
        for j in range(self.n):
            if np.array_equal(H[:, j], s):
                # Fehler in Position j: Bit flippen
                corrected = r.copy()
                corrected[j] ^= 1
                return corrected

        # Fehler nicht korrigierbar (mehr als 1 Fehler bei d=3)
        return r.copy()

    def is_systematic(self) -> bool:
        """
        @brief Prüft, ob der Code systematisch ist (G enthält eine k×k-Einheitsmatrix).
        @description
            Ein systematischer Code enthält die Nachrichtenbits unverändert
            im Codewort: G = [Ik | P] oder G = [P | Ik].
            Prüft, ob Ik als zusammenhängende Teilmatrix in G vorkommt.
        @return True wenn systematisch, sonst False
        @lastModified 2026-03-10
        """
        Ik = np.eye(self.k, dtype=int)

        # Prüfe ob Einheitsmatrix als erste k Spalten vorkommt
        for start in range(self.n - self.k + 1):
            submatrix = self.G[:, start:start + self.k]
            if np.array_equal(submatrix, Ik):
                return True
        return False

    def dual_code(self) -> "LinearCode":
        """
        @brief Berechnet den dualen Code C⊥.
        @description
            Der duale Code C⊥ hat die Prüfmatrix H von C als Generatormatrix:
                C⊥ = { x ∈ GF(2)^n : x·cᵀ = 0 für alle c ∈ C }
            Dimensionsformel: dim(C⊥) = n - k
        @return LinearCode-Objekt, das C⊥ repräsentiert
        @lastModified 2026-03-10
        """
        H = self.parity_check_matrix()
        return LinearCode(H)


# ============================================================
# 3. HammingCode – Hamming-Codes [2^r-1, 2^r-r-1, 3]
# ============================================================

class HammingCode:
    """
    @class HammingCode
    @brief Hamming-Code der Ordnung r: [2^r-1, 2^r-r-1, 3]-Code.
    @description
        Hamming-Codes sind perfekte 1-Fehler-korrigierende lineare Codes.
        Sie erfüllen die Hamming-Schranke mit Gleichheit (perfekte Kugelfüllung).

        Parameter:
            n = 2^r - 1     (Codewortlänge)
            k = 2^r - r - 1 (Informationsbits)
            d = 3           (Minimaldistanz → 1 Fehler korrigierbar)

        Prüfmatrix H:
            Spalten = alle Bitvektoren von 1 bis 2^r-1 (Länge r)
            H hat Dimension r × (2^r-1)

        HINWEIS: Der (7,4)-Hamming-Code und grundlegende Hamming-Methoden sind
        in information_theory.py (ErrorCorrection) implementiert. Diese Klasse
        bietet die allgemeine algebraische Konstruktion für beliebiges r.

    @param r  Ordnung des Hamming-Codes (r ≥ 2)

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, r: int) -> None:
        """
        @brief Initialisiert den Hamming-Code der Ordnung r.
        @param r  Ordnung ≥ 2
        @raises ValueError  Wenn r < 2
        @lastModified 2026-03-10
        """
        if r < 2:
            raise ValueError(f"Ordnung r muss ≥ 2 sein, ist aber {r}.")
        self.r = r
        self.n = (1 << r) - 1       # 2^r - 1
        self.k = self.n - r          # 2^r - r - 1
        self.d = 3                   # Minimaldistanz immer 3

    def parameters(self) -> tuple[int, int, int]:
        """
        @brief Gibt die Code-Parameter (n, k, d) zurück.
        @return Tupel (n, k, d) mit n=2^r-1, k=2^r-r-1, d=3
        @lastModified 2026-03-10
        """
        return (self.n, self.k, self.d)

    def parity_check_matrix(self) -> np.ndarray:
        """
        @brief Konstruiert die Prüfmatrix H des Hamming-Codes.
        @description
            H ist eine r × n-Matrix, deren Spalten alle Bitvektoren
            {1, 2, ..., 2^r-1} (je r Bit lang) in aufsteigender Reihenfolge sind.

            Beispiel (r=2, [3,1,3]-Code):
                H = [[0, 1, 1],
                     [1, 0, 1]]
            (Spalten: 01, 10, 11)

        @return Prüfmatrix H der Dimension r × (2^r-1)
        @lastModified 2026-03-10
        """
        H = np.zeros((self.r, self.n), dtype=int)
        for col, val in enumerate(range(1, self.n + 1)):
            # Bit-Darstellung von val (von MSB nach LSB)
            for row in range(self.r):
                H[row, col] = (val >> (self.r - 1 - row)) & 1
        return H

    def generator_matrix(self) -> np.ndarray:
        """
        @brief Berechnet die Generatormatrix G des Hamming-Codes.
        @description
            G wird aus der Prüfmatrix H abgeleitet:
            G hat Dimension k × n, H·Gᵀ = 0 (mod 2).

            Algorithmus:
            1. Prüfmatrix H aufstellen
            2. Generatormatrix via LinearCode-Klasse bestimmen (dualer Ansatz)
               Alternativ: direkte Konstruktion über systematische Form.

        @return Generatormatrix G der Dimension k × n
        @lastModified 2026-03-10
        """
        H = self.parity_check_matrix()
        # Kern von H bestimmen (Zeilenstufenform über GF(2))
        # G = Menge aller Vektoren x mit H·xᵀ = 0

        # Erweiterte Matrix [H | I_n] zum Bestimmen des Kerns
        # Einfacherer Weg: alle 2^n Vektoren testen (nur für kleine n)
        # Für allgemeines r: systematische Konstruktion

        # Systematische Prüfmatrix: sortiere Spalten so, dass [A | I_r] entsteht
        # Standard-Hamming: Spalten 1..r sind Einheitsvektoren e_r...e_1
        # Wir nehmen direkt H und leiten G systematisch ab

        # Finde Pivot-Spalten (Spalten die Einheitsvektoren bilden)
        H_work = H.copy()
        pivot_cols: list[int] = []
        row = 0
        for col in range(self.n):
            if row >= self.r:
                break
            pivot = -1
            for r_idx in range(row, self.r):
                if H_work[r_idx, col] == 1:
                    pivot = r_idx
                    break
            if pivot == -1:
                continue
            H_work[[row, pivot]] = H_work[[pivot, row]]
            pivot_cols.append(col)
            for r_idx in range(self.r):
                if r_idx != row and H_work[r_idx, col] == 1:
                    H_work[r_idx] = (H_work[r_idx] + H_work[row]) % 2
            row += 1

        non_pivot_cols = [c for c in range(self.n) if c not in pivot_cols]

        # G: Zeilen aus Nicht-Pivot-Spalten aufbauen
        G = np.zeros((self.k, self.n), dtype=int)
        for i, nc in enumerate(non_pivot_cols):
            G[i, nc] = 1
            for j, pc in enumerate(pivot_cols):
                G[i, pc] = H_work[j, nc]

        return G

    def encode(self, message: list[int] | np.ndarray) -> np.ndarray:
        """
        @brief Kodiert eine k-Bit-Nachricht zum Hamming-Codewort.
        @param message  Nachrichtenvektor der Länge k
        @return Codewort der Länge n
        @lastModified 2026-03-10
        """
        G = self.generator_matrix()
        code = LinearCode(G)
        return code.encode(message)

    def decode(self, received: list[int] | np.ndarray) -> dict[str, Any]:
        """
        @brief Dekodiert ein empfangenes Wort und korrigiert bis zu 1 Fehler.
        @description
            Syndrom-Dekodierung:
            1. s = H·rᵀ (mod 2)
            2. Falls s = 0: kein Fehler
            3. s als Binärzahl interpretiert → Fehlerposition (1-indiziert)
            4. Bit in dieser Position flippen
            5. Nachricht aus den k Nicht-Prüfbit-Positionen extrahieren

        @param received  Empfangener Vektor der Länge n
        @return Dictionary mit:
                - 'corrected': korrigiertes Codewort
                - 'error_position': Fehlerposition (1-basiert, 0 = kein Fehler)
                - 'message': extrahierte Nachricht (k Bits)
        @raises ValueError  Wenn Länge des Vektors ≠ n
        @lastModified 2026-03-10
        """
        r_vec = np.array(received, dtype=int) % 2
        if len(r_vec) != self.n:
            raise ValueError(
                f"Empfangene Länge {len(r_vec)} ≠ n={self.n}."
            )

        H = self.parity_check_matrix()
        # Syndrom berechnen
        syndrome = (H @ r_vec) % 2

        # Fehlerposition aus Syndrom (Binärzahl, MSB zuerst)
        error_pos = 0
        for bit in syndrome:
            error_pos = (error_pos << 1) | int(bit)

        corrected = r_vec.copy()
        if error_pos != 0:
            # Fehlerposition ist 1-basiert
            corrected[error_pos - 1] ^= 1

        # Nachricht aus den k Nicht-Prüfbit-Positionen extrahieren
        G = self.generator_matrix()
        code = LinearCode(G)
        H2, pivot_cols = LinearCode._row_reduce_gf2(H.copy())
        non_pivot_cols = [c for c in range(self.n) if c not in pivot_cols]

        message = corrected[non_pivot_cols]

        return {
            'corrected': corrected,
            'error_position': error_pos,
            'message': message,
        }


# ============================================================
# 4. CyclicCode – Zyklische Codes
# ============================================================

class CyclicCode:
    """
    @class CyclicCode
    @brief Zyklischer Code der Länge n mit Generatorpolynom g(x).
    @description
        Ein zyklischer Code C der Länge n über GF(2) ist ein linearer Code,
        bei dem für jedes c = (c_0,...,c_{n-1}) ∈ C auch die zyklische
        Verschiebung (c_{n-1}, c_0,...,c_{n-2}) ∈ C gilt.

        Algebraisch: C ist ein Ideal in GF(2)[x]/(x^n - 1),
        erzeugt von g(x) mit g(x) | (x^n - 1).

        Codierung: c(x) = m(x) · g(x) (mod x^n - 1)
        Syndrom: s(x) = r(x) mod g(x)

    @param n              Codewortlänge
    @param generator_poly Generatorpolynom als Koeffizientenliste [a0, a1, ..., ak]
                          (a0 = Konstantterm, ak = Leitterm)

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, n: int, generator_poly: list[int]) -> None:
        """
        @brief Initialisiert den zyklischen Code.
        @param n              Codewortlänge
        @param generator_poly Generatorpolynom-Koeffizienten über GF(2)
        @raises ValueError  Wenn Generatorpolynom nicht (x^n-1) teilt
        @lastModified 2026-03-10
        """
        self.n = n
        # Koeffizienten mod 2, führende Nullen entfernen
        self.g = [int(c) % 2 for c in generator_poly]
        while len(self.g) > 1 and self.g[-1] == 0:
            self.g.pop()

        # Grad des Generatorpolynoms
        self.deg_g = len(self.g) - 1
        self.k = n - self.deg_g  # Informationsbits

    # ----------------------------------------------------------
    # Polynom-Hilfsmethoden (statisch, über GF(2))
    # ----------------------------------------------------------

    @staticmethod
    def _poly_add(a: list[int], b: list[int]) -> list[int]:
        """Addition zweier Polynome über GF(2)."""
        length = max(len(a), len(b))
        result = [0] * length
        for i, c in enumerate(a):
            result[i] ^= c
        for i, c in enumerate(b):
            result[i] ^= c
        # Führende Nullen entfernen
        while len(result) > 1 and result[-1] == 0:
            result.pop()
        return result

    @staticmethod
    def _poly_mul(a: list[int], b: list[int]) -> list[int]:
        """Multiplikation zweier Polynome über GF(2)."""
        if not a or not b:
            return [0]
        result = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                result[i + j] ^= (ai * bj) % 2
        return result

    @staticmethod
    def _poly_divmod(dividend: list[int], divisor: list[int]) -> tuple[list[int], list[int]]:
        """
        Division mit Rest für Polynome über GF(2).
        @return (Quotient, Rest)
        """
        if all(c == 0 for c in divisor):
            raise ZeroDivisionError("Division durch Nullpolynom.")

        # Kopie anlegen und führende Nullen entfernen
        rem = list(dividend)
        while len(rem) > 1 and rem[-1] == 0:
            rem.pop()
        div = list(divisor)
        while len(div) > 1 and div[-1] == 0:
            div.pop()

        deg_div = len(div) - 1
        quot: list[int] = []

        while len(rem) >= len(div):
            # Koeffizient des Führungsterms
            coef = rem[-1]  # in GF(2) immer 1 (wenn Polynom nicht 0)
            quot.insert(0, coef)
            if coef:
                # Subtrahiere (= Addition in GF(2)) divisor * x^(deg_rem - deg_div)
                shift = len(rem) - len(div)
                for i, d in enumerate(div):
                    rem[shift + i] ^= (coef * d) % 2
            rem.pop()

        if not rem:
            rem = [0]

        return quot, rem

    # ----------------------------------------------------------
    # Öffentliche Methoden
    # ----------------------------------------------------------

    def encode_poly(self, message_poly: list[int]) -> list[int]:
        """
        @brief Kodiert ein Nachrichtenpolynom m(x) → c(x) = m(x) · g(x).
        @description
            c(x) = m(x) · g(x)  (mod x^n - 1, Grad ≤ n-1)
            Das Ergebnis wird auf Länge n aufgefüllt/gekürzt.

        @param message_poly  Nachrichtenpolynom der Länge k (Koeffizientenliste)
        @return Codewort der Länge n als Liste
        @lastModified 2026-03-10
        """
        m = [int(c) % 2 for c in message_poly]
        # Codierung: c = m * g
        c = self._poly_mul(m, self.g)
        # Auf Länge n trimmen/auffüllen
        c = c[:self.n]
        c += [0] * (self.n - len(c))
        return c

    def syndrome_poly(self, received_poly: list[int]) -> list[int]:
        """
        @brief Berechnet das Syndrompolynom s(x) = r(x) mod g(x).
        @description
            Falls s(x) = 0: kein (erkannter) Fehler
            Falls s(x) ≠ 0: Fehler vorhanden
        @param received_poly  Empfangenes Polynom (Länge n)
        @return Syndrompolynom (Grad < deg(g))
        @lastModified 2026-03-10
        """
        r = [int(c) % 2 for c in received_poly]
        _, remainder = self._poly_divmod(r, self.g)
        return remainder

    def cyclic_shift(self, codeword: list[int]) -> list[int]:
        """
        @brief Zyklische Verschiebung eines Codewortes um eine Position.
        @description
            (c_0, c_1, ..., c_{n-1}) → (c_{n-1}, c_0, ..., c_{n-2})
            Entspricht Multiplikation mit x modulo x^n - 1.
        @param codeword  Codewort der Länge n
        @return Zyklisch verschobenes Codewort
        @lastModified 2026-03-10
        """
        c = list(codeword)
        if not c:
            return c
        return [c[-1]] + c[:-1]

    def is_cyclic(self, code_words: list[list[int]]) -> bool:
        """
        @brief Prüft, ob eine Menge von Codewörtern zyklisch abgeschlossen ist.
        @description
            Ein Code C heißt zyklisch, wenn für jedes c ∈ C auch
            alle zyklischen Verschiebungen von c in C liegen.
        @param code_words  Liste aller Codewörter
        @return True wenn zyklisch, sonst False
        @lastModified 2026-03-10
        """
        # Menge für schnellen Lookup
        code_set = {tuple(c) for c in code_words}

        for cw in code_words:
            shifted = self.cyclic_shift(list(cw))
            if tuple(shifted) not in code_set:
                return False
        return True


# ============================================================
# 5. BCHCode – BCH-Codes und Berlekamp-Massey
# ============================================================

class BCHCode:
    """
    @class BCHCode
    @brief BCH-Codes (Bose-Chaudhuri-Hocquenghem), Reed-Solomon-Demo,
           Berlekamp-Massey-Algorithmus.
    @description
        BCH-Codes sind eine verallgemeinerung der Hamming-Codes und erlauben
        die Korrektur von t Fehlern. Sie sind zyklische Codes.

        BCH-Konstruktion (binär):
        - Wähle primitive n-te Einheitswurzel α in GF(2^m)
        - Generatorpolynom g(x) = kgV der Minimalpolynome von α, α², ..., α^(2t)
        - Code hat Länge n = 2^m - 1, korrigiert t Fehler

        Reed-Solomon-Code:
        - Spezialfall über GF(q) mit n = q-1
        - Maximale Distanz separierbar (MDS-Code): d = n - k + 1 (Singleton-Bound)

        Berlekamp-Massey-Algorithmus:
        - Findet das kürzeste LFSR (Linear Feedback Shift Register), das eine
          gegebene Sequenz erzeugt
        - Grundlage der BCH/RS-Dekodierung

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def bch_generator_poly(n: int, t: int, q: int = 2) -> dict[str, Any]:
        """
        @brief Berechnet das BCH-Generatorpolynom (Demo für binäre BCH-Codes).
        @description
            Für einen binären BCH-Code der Länge n = 2^m - 1 mit Korrekturfähigkeit t:
            g(x) = kgV(m_1(x), m_2(x), ..., m_{2t}(x))
            wobei m_i(x) das Minimalpolynom der primitiven Einheitswurzel α^i ist.

            Diese Demo-Implementierung gibt die Codeparameter zurück.

            Wichtige bekannte BCH-Codes:
            - [7,4,3]: r=3, t=1 (= Hamming-Code)
            - [15,7,5]: r=4, t=2
            - [15,5,7]: r=4, t=3
            - [31,21,5]: r=5, t=2

        @param n  Codewortlänge (sollte 2^m - 1 für primitive BCH-Codes sein)
        @param t  Korrekturfähigkeit (Anzahl korrigierbarer Fehler)
        @param q  Grundkörpergröße (Standard: 2 für binäre Codes)
        @return Dictionary mit Code-Parametern und Beschreibung
        @lastModified 2026-03-10
        """
        # Mindestdistanz: d ≥ 2t + 1 (BCH-Schranke)
        d_min = 2 * t + 1

        # Anzahl Prüfbits: mindestens m*t (m = log2(n+1) für primitive Codes)
        import math
        m = math.ceil(math.log2(n + 1)) if n > 0 else 1
        r_bits = m * t  # Anzahl Redundanzbits

        # Informationsbits
        k = max(0, n - r_bits)

        # Bekannte primitive BCH-Generatorpolynome (GF(2), Demo)
        known_polys: dict[tuple[int, int], list[int]] = {
            (7, 1): [1, 0, 1, 1],           # x³+x+1 = [1,1,0,1]
            (7, 2): [1, 1, 0, 1, 1, 1, 1],   # [7,1,7]-BCH
            (15, 1): [1, 0, 0, 1, 1],         # x⁴+x+1
            (15, 2): [1, 0, 1, 0, 0, 1, 1, 0, 1],  # Grad 8
            (15, 3): [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1],  # Grad 10
        }

        gen_poly = known_polys.get((n, t), None)

        return {
            'n': n,
            'k': k,
            't': t,
            'd_min': d_min,
            'q': q,
            'm': m,
            'redundancy_bits': r_bits,
            'generator_poly': gen_poly,
            'description': (
                f"BCH-Code [{n},{k},{d_min}] über GF({q}), "
                f"korrigiert {t} Fehler, BCH-Schranke: d ≥ {d_min}"
            ),
        }

    @staticmethod
    def reed_solomon_demo(n: int, k: int, q: int) -> dict[str, Any]:
        """
        @brief Informations-Demo für Reed-Solomon-Codes.
        @description
            Reed-Solomon-Codes sind MDS-Codes (Maximum Distance Separable):
                d = n - k + 1  (Singleton-Bound wird exakt erfüllt)

            Eigenschaften:
            - Länge n ≤ q-1 (über GF(q))
            - Korrigiert t = ⌊(d-1)/2⌋ = ⌊(n-k)/2⌋ Fehler
            - Wichtig für: CD, DVD, QR-Code, RAID-6, Deep Space Communication

            Standard-RS: n = q-1 = 255 (GF(2^8) = GF(256))
            RS(255, 223): t = 16 Fehler korrigierbar (NASA Voyager, CD-ROM)

        @param n  Codewortlänge (n ≤ q-1)
        @param k  Informationssymbole
        @param q  Grundkörpergröße (q = p^m, Primzahlpotenz)
        @return Dictionary mit RS-Code-Informationen
        @lastModified 2026-03-10
        """
        d = n - k + 1       # Minimaldistanz (Singleton-Schranke)
        t = (d - 1) // 2   # Korrekturfähigkeit

        return {
            'n': n,
            'k': k,
            'd': d,
            't': t,
            'q': q,
            'rate': k / n if n > 0 else 0,
            'is_mds': True,  # RS-Codes sind immer MDS
            'description': (
                f"Reed-Solomon-Code RS({n},{k}) über GF({q}): "
                f"d={d}, korrigiert {t} Fehler/Löschungen"
            ),
            'applications': [
                'CD/DVD Fehlerkorrektur',
                'QR-Codes',
                'RAID-6 Speichersysteme',
                'NASA Deep Space Communication',
                'DSL/ADSL Modems',
            ],
        }

    @staticmethod
    def berlekamp_massey(sequence: list[int]) -> dict[str, Any]:
        """
        @brief Berlekamp-Massey-Algorithmus für binäre Sequenzen.
        @description
            Findet das kürzeste LFSR (Linear Feedback Shift Register), das die
            gegebene binäre Sequenz erzeugt.

            Algorithmus (Massey 1969):
            Initialisierung: C = [1], B = [1], L = 0, m = 1, b = 1
            Für jedes Symbol s[n]:
                d = s[n] ⊕ (C * s[n-1..n-L])  (Diskrepanz)
                Falls d = 0: m++
                Falls d = 1:
                    T = C, C = C ⊕ (d/b) * x^m * B, B = T, L = n+1-L, b = d, m = 1

            Ausgabe:
            - Verbindungspolynom C(x) (LFSR-Feedback-Polynom)
            - LFSR-Länge L (= Linearkomplexität der Sequenz)
            - Linearkomplexität λ(s) = L

        @param sequence  Binäre Sequenz (Liste aus 0 und 1)
        @return Dictionary mit:
                - 'connection_poly': Verbindungspolynom C(x) als Koeffizientenliste
                - 'lfsr_length': LFSR-Länge (Linearkomplexität)
                - 'sequence_length': Länge der Eingangssequenz
        @lastModified 2026-03-10
        """
        s = [int(x) % 2 for x in sequence]
        N = len(s)

        if N == 0:
            return {
                'connection_poly': [1],
                'lfsr_length': 0,
                'sequence_length': 0,
            }

        # Initialisierung
        C = [1]   # Verbindungspolynom (startet mit 1)
        B = [1]   # Hilfpolynom
        L = 0     # Aktuelle LFSR-Länge
        m = 1     # Verschiebungszähler

        for n in range(N):
            # Diskrepanz d berechnen: d = s[n] + Σ C[i]*s[n-i] für i=1..L
            d = s[n]
            for i in range(1, L + 1):
                if i < len(C) and n - i >= 0:
                    d ^= C[i] * s[n - i]
            d %= 2

            if d == 0:
                # Kein Fehler: Zähler erhöhen
                m += 1
            else:
                # Verbindungspolynom anpassen
                T = list(C)
                # C = C + B verschoben um m Positionen (x^m * B)
                B_shifted = [0] * m + B
                # C = C XOR B_shifted (über GF(2))
                new_len = max(len(C), len(B_shifted))
                C_new = [0] * new_len
                for i in range(len(C)):
                    C_new[i] ^= C[i]
                for i in range(len(B_shifted)):
                    C_new[i] ^= B_shifted[i]
                C = C_new

                if 2 * L <= n:
                    L = n + 1 - L
                    B = T
                    m = 1
                else:
                    m += 1

        return {
            'connection_poly': C,
            'lfsr_length': L,
            'sequence_length': N,
        }


# ============================================================
# 6. FormalLanguages – Formale Sprachen und Automaten
# ============================================================

class DFA:
    """
    @class DFA
    @brief Deterministischer Endlicher Automat (Deterministic Finite Automaton).
    @description
        Ein DFA ist ein 5-Tupel M = (Q, Σ, δ, q₀, F) mit:
        - Q: endliche Menge von Zuständen
        - Σ: Eingabealphabet
        - δ: Q × Σ → Q: Übergangsfunktion
        - q₀ ∈ Q: Startzustand
        - F ⊆ Q: Menge akzeptierender Zustände

        Akzeptanz: M akzeptiert w = a₁a₂...aₙ genau dann, wenn
        δ*(q₀, w) ∈ F (δ* = erweiterte Übergangsfunktion).

    @param states      Menge der Zustände
    @param alphabet    Eingabealphabet (Menge von Symbolen)
    @param transitions Übergangsfunktion als Dict {(state, symbol): next_state}
    @param start       Startzustand
    @param accepting   Menge akzeptierender Zustände

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        states: set,
        alphabet: set,
        transitions: dict[tuple, Any],
        start: Any,
        accepting: set,
    ) -> None:
        """
        @brief Initialisiert den DFA.
        @param states      Zustandsmenge Q
        @param alphabet    Eingabealphabet Σ
        @param transitions Übergangsfunktion δ: {(q, a): q'}
        @param start       Startzustand q₀
        @param accepting   Akzeptierende Zustände F
        @lastModified 2026-03-10
        """
        self.states = set(states)
        self.alphabet = set(alphabet)
        self.transitions = dict(transitions)
        self.start = start
        self.accepting = set(accepting)

    def accepts(self, word: str | list) -> bool:
        """
        @brief Prüft, ob der DFA das Wort w akzeptiert.
        @description
            Simulation: Startzustand q₀, lese Symbol für Symbol,
            wende Übergangsfunktion an.
            Akzeptanz: Endzustand ∈ F
            Leeres Wort ε: Akzeptiert genau dann, wenn q₀ ∈ F.

        @param word  Eingabewort (String oder Liste von Symbolen)
        @return True wenn akzeptiert, sonst False
        @lastModified 2026-03-10
        """
        current = self.start

        for symbol in word:
            key = (current, symbol)
            if key not in self.transitions:
                # Kein Übergang definiert → Senkzustand (Verwerfung)
                return False
            current = self.transitions[key]

        return current in self.accepting

    def minimize(self) -> dict[str, Any]:
        """
        @brief Minimierung des DFA via Myhill-Nerode-Äquivalenz (Demo).
        @description
            Algorithmus (Tabellenausfüll-Methode / Hopcroft):
            1. Initiale Partition: {F, Q\\F}
            2. Verfeinere: zwei Zustände p,q äquivalent falls
               ∀a ∈ Σ: δ(p,a) und δ(q,a) in selber Klasse
            3. Wiederhole bis stabil

            Diese Demo-Implementierung gibt die Äquivalenzklassen zurück.

        @return Dictionary mit:
                - 'equivalence_classes': Liste von Äquivalenzklassen
                - 'num_states_original': Anzahl ursprünglicher Zustände
                - 'num_states_minimal': Anzahl Zustände im minimalen DFA
        @lastModified 2026-03-10
        """
        # Initiale Partition: akzeptierende vs. nicht-akzeptierende Zustände
        non_accepting = self.states - self.accepting
        partition: list[frozenset] = []
        if self.accepting:
            partition.append(frozenset(self.accepting))
        if non_accepting:
            partition.append(frozenset(non_accepting))

        # Iterative Verfeinerung
        changed = True
        while changed:
            changed = False
            new_partition: list[frozenset] = []

            for group in partition:
                # Versuche Gruppe zu verfeinern
                # Zwei Zustände sind äquivalent, wenn ihre Nachfolger
                # in denselben Klassen liegen
                sub_groups: dict[tuple, set] = {}

                for state in group:
                    # Signatur: für jedes Symbol die Partitionsklasse des Nachfolgers
                    sig: list[int] = []
                    for symbol in sorted(self.alphabet, key=str):
                        key = (state, symbol)
                        if key in self.transitions:
                            next_state = self.transitions[key]
                            # Klassen-Index des Nachfolgers
                            cls_idx = -1
                            for idx, cls in enumerate(partition):
                                if next_state in cls:
                                    cls_idx = idx
                                    break
                        else:
                            cls_idx = -2  # Kein Übergang = Senkzustand
                        sig.append(cls_idx)

                    sig_tuple = tuple(sig)
                    if sig_tuple not in sub_groups:
                        sub_groups[sig_tuple] = set()
                    sub_groups[sig_tuple].add(state)

                if len(sub_groups) > 1:
                    changed = True
                for sub in sub_groups.values():
                    new_partition.append(frozenset(sub))

            partition = new_partition

        return {
            'equivalence_classes': [set(cls) for cls in partition],
            'num_states_original': len(self.states),
            'num_states_minimal': len(partition),
        }


class NFA:
    """
    @class NFA
    @brief Nichtdeterministischer Endlicher Automat (Nondeterministic Finite Automaton).
    @description
        Ein NFA ist ein 5-Tupel M = (Q, Σ, δ, q₀, F) mit:
        - δ: Q × Σ → 2^Q: Übergangsfunktion (liefert Menge von Zuständen)
        - Akzeptanz: Es existiert eine Berechnung, die in F endet.

        Simulation: Teilmengenkonstruktion (Powerset-Konstruktion).
        Aus NFA mit n Zuständen entsteht äquivalenter DFA mit bis zu 2^n Zuständen.

    @param states      Zustandsmenge Q
    @param alphabet    Alphabet Σ
    @param transitions Übergangsfunktion {state: {symbol: set_of_states}}
    @param start       Startzustand q₀
    @param accepting   Akzeptierende Zustände F

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        states: set,
        alphabet: set,
        transitions: dict,
        start: Any,
        accepting: set,
    ) -> None:
        """
        @brief Initialisiert den NFA.
        @param transitions  Dict {state: {symbol: set(next_states)}}
        @lastModified 2026-03-10
        """
        self.states = set(states)
        self.alphabet = set(alphabet)
        self.transitions = transitions  # {state: {symbol: set}}
        self.start = start
        self.accepting = set(accepting)

    def _step(self, current_states: set, symbol: Any) -> set:
        """
        @brief Berechnet die Menge der Folgezustände für ein Symbol.
        @param current_states  Aktuelle Zustandsmenge
        @param symbol          Eingabesymbol
        @return Menge der Folgezustände
        @lastModified 2026-03-10
        """
        next_states: set = set()
        for state in current_states:
            if state in self.transitions:
                state_trans = self.transitions[state]
                if symbol in state_trans:
                    targets = state_trans[symbol]
                    if isinstance(targets, set):
                        next_states |= targets
                    else:
                        next_states.add(targets)
        return next_states

    def accepts(self, word: str | list) -> bool:
        """
        @brief Prüft, ob der NFA das Wort w akzeptiert (Teilmengenkonstruktion).
        @description
            Simuliert alle möglichen Berechnungen parallel.
            Startet mit {q₀}, verarbeitet Symbole, akzeptiert wenn
            Endzustandsmenge ∩ F ≠ ∅.
        @param word  Eingabewort
        @return True wenn mindestens ein akzeptierender Pfad existiert
        @lastModified 2026-03-10
        """
        current = {self.start}

        for symbol in word:
            current = self._step(current, symbol)
            if not current:
                return False  # Alle Berechnungen sind abgestorben

        # Akzeptanz: mindestens ein akzeptierender Zustand erreichbar
        return bool(current & self.accepting)


class RegularExpression:
    """
    @class RegularExpression
    @brief Einfache reguläre Ausdrücke mit *, +, |, () und Literalen.
    @description
        Implementiert eine vereinfachte Teilmenge regulärer Ausdrücke:
        - Literale: a, b, c, ... (einzelne Zeichen)
        - Konkatenation: ab (a gefolgt von b)
        - Alternation: a|b (a oder b)
        - Kleene-Stern: a* (beliebig viele a's, auch 0)
        - Plus: a+ (mindestens ein a)
        - Gruppierung: (ab)* (Gruppe wiederholen)

        Implementierung: Python's re-Modul (vollständig kompatibel).

    @param pattern  Regulärer Ausdruck als String

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, pattern: str) -> None:
        """
        @brief Initialisiert den regulären Ausdruck.
        @param pattern  Muster (Python-RE-Syntax)
        @raises re.error  Wenn das Muster syntaktisch ungültig ist
        @lastModified 2026-03-10
        """
        self.pattern = pattern
        # Vollständiges Match erzwingen (kein Teilstring-Match)
        self._compiled = _re.compile(f'^(?:{pattern})$')

    def matches(self, word: str) -> bool:
        """
        @brief Prüft ob das Wort vollständig zum Muster passt.
        @param word  Zu prüfendes Wort
        @return True wenn das Wort vom Muster akzeptiert wird, sonst False
        @lastModified 2026-03-10
        """
        return bool(self._compiled.match(word))


class ContextFreeGrammar:
    """
    @class ContextFreeGrammar
    @brief Kontextfreie Grammatik G = (T, N, P, S).
    @description
        Eine KFG besteht aus:
        - T: terminale Symbole (Alphabet)
        - N: nichtterminale Symbole (Variablen)
        - P: Produktionsregeln {A → α | A ∈ N, α ∈ (T ∪ N)*}
        - S ∈ N: Startsymbol

        Erzeugungsrelation: S ⊢* w  (Ableitung in beliebig vielen Schritten)
        Eine KFG erzeugt die Sprache L(G) = {w ∈ T* : S ⊢* w}.

        Wichtige Normalformen:
        - Chomsky-Normalform (CNF): A → BC | a
        - Greibach-Normalform (GNF): A → aα

        Beispiel: Palindrome über {a, b}:
            S → aSa | bSb | a | b | ε

    @param terminals      Menge terminaler Symbole
    @param nonterminals   Menge nichtterminaler Symbole
    @param rules          Produktionsregeln {A: [['B', 'C'], ['a'], ...]}
    @param start          Startsymbol

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        terminals: set[str],
        nonterminals: set[str],
        rules: dict[str, list[list[str]]],
        start: str,
    ) -> None:
        """
        @brief Initialisiert die kontextfreie Grammatik.
        @param terminals     Menge terminaler Symbole (z.B. {'a', 'b'})
        @param nonterminals  Menge nichtterminaler Symbole (z.B. {'S', 'A'})
        @param rules         Produktionen {NT: [[Symbol, ...], ...]}
                             Leere Produktion: [[]] oder [['']]
        @param start         Startsymbol
        @lastModified 2026-03-10
        """
        self.terminals = set(terminals)
        self.nonterminals = set(nonterminals)
        self.rules = dict(rules)
        self.start = start

    def generates(self, word: str, max_depth: int = 10) -> bool:
        """
        @brief Prüft per Backtracking, ob die Grammatik das Wort erzeugt.
        @description
            Rekursives Backtracking über alle Ableitungsschritte.
            Startet mit S und versucht alle Produktionen anzuwenden.

            Abbruchbedingungen:
            - Aktuelle Ableitung = word (Erfolg)
            - Maximale Tiefe erreicht (Abbruch)
            - Aktuelle Ableitung länger als word (Abbruch bei nicht-löschenden Regeln)

            Leeres Wort ε: word = '' → prüft ob S →* ε

        @param word       Zu prüfendes Wort (Zeichenkette über T)
        @param max_depth  Maximale Ableitungstiefe
        @return True wenn w ∈ L(G), sonst False
        @lastModified 2026-03-10
        """
        def derive(sentential: list[str], depth: int) -> bool:
            """Rekursive Ableitung einer Satzform."""
            if depth > max_depth:
                return False

            # Satzform zu String konvertieren
            s = ''.join(sentential)
            if s == word:
                return True

            # Wenn alle Symbole terminal sind und ≠ word: Fehlschlag
            if all(sym in self.terminals or sym == '' for sym in sentential):
                return False

            # Erstes Nichtterminal in der Satzform finden und expandieren
            for i, sym in enumerate(sentential):
                if sym in self.nonterminals:
                    # Alle Produktionen für dieses Nichtterminal versuchen
                    for production in self.rules.get(sym, []):
                        # Neue Satzform: sym durch Produktion ersetzen
                        if production == [''] or production == []:
                            new_sent = sentential[:i] + sentential[i+1:]
                        else:
                            new_sent = sentential[:i] + production + sentential[i+1:]

                        # Längenbasierter Pruning (nur wenn Wort nicht leer)
                        joined = ''.join(new_sent)
                        if len(joined) > len(word) + max_depth:
                            continue

                        if derive(new_sent, depth + 1):
                            return True
                    # Nur erstes Nichtterminal verarbeiten (leftmost derivation)
                    return False

            return False

        return derive([self.start], 0)

    def is_ambiguous_demo(self) -> dict[str, Any]:
        """
        @brief Demo: Erläuterung von Ambiguität in KFGs.
        @description
            Eine Grammatik G heißt ambig, wenn es ein Wort w ∈ L(G) gibt,
            das mindestens zwei verschiedene Ableitungsbäume (Linkableitungen) hat.

            Klassisches Beispiel (Dangling-else-Problem):
                S → if B then S | if B then S else S | a
            Das Wort "if b then if b then a else a" hat zwei Ableitungsbäume.

            Ambiguität ist im Allgemeinen unentscheidbar (KFG-Äquivalenzproblem
            ist unentscheidbar, Ambiguitätsproblem ebenfalls).

        @return Dictionary mit Beispiel und Erläuterung
        @lastModified 2026-03-10
        """
        return {
            'definition': (
                'Eine KFG G heißt ambig, wenn ∃w ∈ L(G) mit zwei verschiedenen '
                'Ableitungsbäumen (Linksableitungen).'
            ),
            'classic_example': {
                'grammar': 'S → if B then S | if B then S else S | a',
                'ambiguous_word': 'if b then if b then a else a',
                'interpretation_1': 'if b then (if b then a else a)',
                'interpretation_2': '(if b then if b then a) else a',
                'context': 'Dangling-Else-Problem in Programmiersprachen',
            },
            'decidability': (
                'Ambiguität von KFGs ist im Allgemeinen UNENTSCHEIDBAR '
                '(Semi-Entscheidbarkeit: kann beweisen, aber nicht widerlegen).'
            ),
            'current_grammar_nonterminals': list(self.nonterminals),
            'current_grammar_rules': {k: v for k, v in self.rules.items()},
        }


class PushdownAutomaton:
    """
    @class PushdownAutomaton
    @brief Kellerautomat (Pushdown Automaton, PDA) – Demo-Implementierung.
    @description
        Ein PDA ist ein NFA mit zusätzlichem Kellerspeicher (Stack).
        Er akzeptiert genau die kontextfreien Sprachen (Chomsky-Typ 2).

        Formell: P = (Q, Σ, Γ, δ, q₀, Z₀, F) mit:
        - Γ: Kelleralphabet
        - Z₀ ∈ Γ: initialer Kellerinhalt
        - δ: Q × (Σ ∪ {ε}) × Γ → 2^(Q × Γ*): Übergangsfunktion

        Übergangsfunktion: (q, a, A) → {(q', γ) | ...}
        Bedeutung: Im Zustand q, bei Symbol a, oberstem Kellerzeichen A:
                   wechsle zu q', ersetze A durch γ

        Demo: Akzeptanz der Sprache {aⁿbⁿ | n ≥ 0} (klassisches PDA-Beispiel)

    @param states        Zustandsmenge Q
    @param alphabet      Eingabealphabet Σ
    @param stack_alphabet Kelleralphabet Γ
    @param transitions   Übergänge {(q, symbol, stack_top): [(q', new_stack_top_list), ...]}
    @param start         Startzustand q₀
    @param accepting     Akzeptierende Zustände F

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        states: set,
        alphabet: set,
        stack_alphabet: set,
        transitions: dict,
        start: Any,
        accepting: set,
    ) -> None:
        """
        @brief Initialisiert den Kellerautomaten.
        @lastModified 2026-03-10
        """
        self.states = set(states)
        self.alphabet = set(alphabet)
        self.stack_alphabet = set(stack_alphabet)
        self.transitions = transitions
        self.start = start
        self.accepting = set(accepting)

    def accepts_demo(self, word: str | list) -> dict[str, Any]:
        """
        @brief Demonstriert die Verarbeitung des Wortes durch den PDA.
        @description
            Simuliert den PDA mittels rekursiver Tiefensuche über alle
            möglichen Konfigurationen (q, restliches Eingabewort, Kellerinhalt).

            Konfiguration: (Zustand, restliche Eingabe, Keller)
            Anfangskonfiguration: (q₀, w, [Z₀])
            Akzeptanzkonfiguration: (q ∈ F, ε, beliebig) oder (q, ε, [])

            Achtung: Diese Demo-Implementierung hat keine Garantie der Vollständigkeit
            für alle PDA-Typen (Endlosrekursion bei ε-Übergängen möglich).

        @param word  Eingabewort
        @return Dictionary mit Akzeptanzentscheidung und Ablaufprotokoll
        @lastModified 2026-03-10
        """
        word_list = list(word)

        # Initialer Stack: ['Z0'] (Kellerbodenmarker)
        init_stack = ['Z0']

        trace: list[dict] = []
        visited: set[tuple] = set()

        def simulate(state: Any, pos: int, stack: list) -> bool:
            """Rekursive PDA-Simulation."""
            # Konfiguration als Tupel für Zykluserkennung
            config = (state, pos, tuple(stack))
            if config in visited:
                return False
            visited.add(config)

            remaining = word_list[pos:]
            step_info = {
                'state': state,
                'remaining_input': ''.join(remaining),
                'stack': list(stack),
            }
            trace.append(step_info)

            # Akzeptanzbedingung: Eingabe vollständig gelesen und Zustand akzeptierend
            if pos == len(word_list) and state in self.accepting:
                return True

            # Übergang mit aktuellem Symbol (oder ε)
            stack_top = stack[-1] if stack else None
            current_symbol = word_list[pos] if pos < len(word_list) else None

            # Übergänge mit aktuellem Symbol prüfen
            for sym in ([current_symbol] if current_symbol is not None else []) + [None]:
                key = (state, sym, stack_top)
                if key in self.transitions:
                    for (next_state, push_symbols) in self.transitions[key]:
                        new_stack = stack[:-1]  # Stack-Top entfernen
                        # Neue Symbole auf Stack legen (rechts = oben)
                        if push_symbols:
                            new_stack = new_stack + list(push_symbols)
                        new_pos = pos + (1 if sym is not None else 0)
                        if simulate(next_state, new_pos, new_stack):
                            return True

            return False

        accepted = simulate(self.start, 0, init_stack)

        return {
            'word': ''.join(word_list),
            'accepted': accepted,
            'trace_steps': len(trace),
            'description': (
                f"PDA {'akzeptiert' if accepted else 'verwirft'} "
                f"das Wort '{word}'"
            ),
        }


# Namespace-Klasse für Formale Sprachen (Bündelung)
class FormalLanguages:
    """
    @class FormalLanguages
    @brief Namespace-Klasse, die alle formalen Sprachen-Klassen bündelt.
    @description
        Bequeme Re-Exporte von DFA, NFA, RegularExpression,
        ContextFreeGrammar und PushdownAutomaton.
    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """
    DFA = DFA
    NFA = NFA
    RegularExpression = RegularExpression
    ContextFreeGrammar = ContextFreeGrammar
    PushdownAutomaton = PushdownAutomaton


# ============================================================
# 7. ComputabilityTheory – Berechenbarkeitstheorie
# ============================================================

class ComputabilityTheory:
    """
    @class ComputabilityTheory
    @brief Berechenbarkeitstheorie: primitiv-rekursiv, μ-rekursiv,
           Church-Turing-These, unentscheidbare Probleme.
    @description
        Ergänzung zu recursion_theory.py mit Schwerpunkt auf:
        - Primitiv-rekursive Funktionen (Grundfunktionen + Komposition + Rekursion)
        - μ-rekursive (allgemein-rekursive) Funktionen
        - Church-Turing-These und Äquivalenz verschiedener Berechnungsmodelle
        - Unentscheidbare Probleme (Halteproblem, Post, Rice, ...)

    @author Kurt Ingwer
    @date 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def primitive_recursive_demo() -> dict[str, Any]:
        """
        @brief Demo: Primitiv-rekursive Funktionen und ihre Grundoperationen.
        @description
            Primitiv-rekursive Funktionen entstehen durch:
            1. Grundfunktionen:
               - Nullfunktion: Z(n) = 0
               - Nachfolgerfunktion: S(n) = n+1
               - Projektion: P_i^k(x_1,...,x_k) = x_i
            2. Komposition: h(x⃗) = f(g_1(x⃗),...,g_m(x⃗))
            3. Primitive Rekursion:
               h(x⃗, 0) = f(x⃗)
               h(x⃗, n+1) = g(x⃗, n, h(x⃗, n))

            Beispiele primitiv-rekursiver Funktionen:
            - Addition, Multiplikation, Exponentiation
            - Vorgänger, Subtraktion (abgeschnitten)
            - Fibonacci, Fakultät

        @return Dictionary mit Beispielen und Implementierungen
        @lastModified 2026-03-10
        """
        # Grundfunktionen
        zero: Callable = lambda *args: 0
        successor: Callable = lambda n: n + 1
        projection: Callable = lambda i: (lambda *args: args[i])

        # Primitiv-rekursive Implementierungen
        def add(x: int, y: int) -> int:
            """Addition via primitiver Rekursion: add(x,0)=x, add(x,n+1)=S(add(x,n))"""
            result = x
            for _ in range(y):
                result = successor(result)
            return result

        def multiply(x: int, y: int) -> int:
            """Multiplikation: mul(x,0)=0, mul(x,n+1)=add(mul(x,n),x)"""
            result = 0
            for _ in range(y):
                result = add(result, x)
            return result

        def factorial(n: int) -> int:
            """Fakultät: fac(0)=1, fac(n+1)=(n+1)*fac(n)"""
            result = 1
            for k in range(1, n + 1):
                result = multiply(result, k)
            return result

        def fibonacci(n: int) -> int:
            """Fibonacci: fib(0)=0, fib(1)=1, fib(n+1)=fib(n)+fib(n-1)"""
            if n == 0:
                return 0
            a, b = 0, 1
            for _ in range(n - 1):
                a, b = b, add(a, b)
            return b

        return {
            'description': (
                'Primitiv-rekursive Funktionen: geschlossen unter Komposition '
                'und primitiver Rekursion. Alle primitiv-rekursiven Funktionen '
                'sind total (immer terminierend).'
            ),
            'grundfunktionen': {
                'zero': 'Z(n) = 0 für alle n',
                'successor': 'S(n) = n+1',
                'projection': 'P_i^k(x_1,...,x_k) = x_i',
            },
            'examples': {
                'add(3,4)': add(3, 4),
                'multiply(5,6)': multiply(5, 6),
                'factorial(6)': factorial(6),
                'fibonacci(10)': fibonacci(10),
            },
            'limitation': (
                'Nicht alle totalen berechenbaren Funktionen sind primitiv-rekursiv. '
                'Gegenbeispiel: Ackermannfunktion A(m,n) (wächst schneller als jede '
                'primitiv-rekursive Funktion).'
            ),
            'ackermann_values': {
                'A(0,0)': 1,
                'A(1,1)': 3,
                'A(2,2)': 7,
                'A(3,3)': 61,
                'note': 'A(4,4) = 2^(2^(2^(65536)))-3 (nicht praktisch berechenbar)',
            },
        }

    @staticmethod
    def mu_recursive_demo() -> dict[str, Any]:
        """
        @brief Demo: μ-rekursive (allgemein-rekursive) Funktionen.
        @description
            μ-rekursive Funktionen erweitern primitiv-rekursive Funktionen um
            den μ-Operator (Minimaloperator):

            μ-Operator: μy[f(x⃗,y) = 0] = kleinste y mit f(x⃗,y) = 0
                        (partiell falls kein solches y existiert)

            Satz: Die Klasse der μ-rekursiven Funktionen = die Klasse der
            Turing-berechenbaren Funktionen (Church-Turing-These).

            Partiell-rekursive Funktionen können undefiniert sein (divergieren).

        @return Dictionary mit Erläuterungen und Beispielen
        @lastModified 2026-03-10
        """
        def mu_operator(f: Callable, max_search: int = 1000) -> Callable:
            """
            μ-Operator: Findet kleinstes y mit f(y) = 0.
            max_search begrenzt die Suche (Approximation der partiellen Funktion).
            """
            def search(*args: int) -> int | None:
                for y in range(max_search):
                    try:
                        if f(*args, y) == 0:
                            return y
                    except Exception:
                        pass
                return None  # Divergenz (im echten Modell: undefiniert)
            return search

        # Beispiel: μy[y² - x = 0] ≈ floor(√x) (ganzzahlige Wurzel)
        def integer_sqrt(x: int) -> int | None:
            """Ganzzahlige Quadratwurzel via μ-Operator."""
            for y in range(x + 1):
                if y * y == x:
                    return y
                if y * y > x:
                    return None
            return None

        return {
            'description': (
                'μ-rekursive Funktionen = Turing-berechenbare Funktionen. '
                'Der μ-Operator kann Divergenz erzeugen (partielle Funktionen).'
            ),
            'mu_operator': {
                'definition': 'μy[f(x,y)=0] = kleinstes y ≥ 0 mit f(x,y)=0',
                'partiality': 'Undefiniert falls kein solches y existiert',
            },
            'examples': {
                'integer_sqrt(16)': integer_sqrt(16),
                'integer_sqrt(25)': integer_sqrt(25),
                'integer_sqrt(2)': integer_sqrt(2),  # None = nicht exakt
            },
            'equivalence': (
                'μ-rekursiv ≡ Turing-berechenbar ≡ λ-definierbar ≡ WHILE-berechenbar '
                '(Church-Turing-These, 1936)'
            ),
        }

    @staticmethod
    def church_turing_thesis_demo() -> dict[str, Any]:
        """
        @brief Demo: Church-Turing-These und Äquivalenz von Berechnungsmodellen.
        @description
            Church-Turing-These (1936):
            Jede intuitiv berechenbare Funktion ist Turing-berechenbar.

            Äquivalente Berechnungsmodelle (alle berechnen dieselbe Klasse):
            - Turing-Maschinen (Alan Turing, 1936)
            - λ-Kalkül (Alonzo Church, 1936)
            - μ-rekursive Funktionen (Kleene, 1936)
            - RAM-Maschinen
            - WHILE-Programme
            - Registermaschinen

            Diese Thesis ist NICHT beweisbar (sie ist eine Aussage über
            intuitive Berechenbarkeit, nicht über formale Systeme).

        @return Dictionary mit These und Äquivalenzen
        @lastModified 2026-03-10
        """
        return {
            'thesis': (
                'Church-Turing-These (1936): Eine Funktion f: ℕ → ℕ ist '
                'genau dann intuitiv berechenbar, wenn sie Turing-berechenbar ist.'
            ),
            'equivalent_models': {
                'Turing-Maschinen': 'Turing, 1936: Lese/Schreibkopf auf unendlichem Band',
                'Lambda-Kalkül': 'Church, 1936: Funktionale Berechnung via β-Reduktion',
                'μ-Rekursion': 'Kleene, 1936: Primitive Rekursion + μ-Operator',
                'RAM-Maschinen': 'Random Access Machine (abstraktes Computermodell)',
                'WHILE-Programme': 'Imperative Programme mit WHILE-Schleifen',
                'Registermaschinen': 'Maschinen mit unbegrenzt vielen Registern',
            },
            'status': (
                'Die These ist NICHT formal beweisbar – sie ist eine '
                'philosophische Aussage über intuitive Berechenbarkeit. '
                'Bisher wurde kein intuitiv berechenbares Problem gefunden, '
                'das nicht Turing-berechenbar ist.'
            ),
            'implications': {
                'universality': 'Jede ausreichend mächtige Programmiersprache ist Turing-vollständig',
                'limits': 'Gibt es absolut unberechenbare Probleme (Halteproblem)',
                'quantum': 'Quantencomputer berechnen dieselbe Klasse (nur schneller für manche Probleme)',
            },
        }

    @staticmethod
    def undecidable_problems() -> dict[str, Any]:
        """
        @brief Liste bekannter unentscheidbarer Probleme.
        @description
            Ein Problem ist unentscheidbar, wenn keine Turing-Maschine existiert,
            die für jede Eingabe in endlicher Zeit korrekt Ja/Nein antwortet.

            Unentscheidbarkeit via Reduktion:
            Problem A ≤ₘ Problem B: Wenn A auf B reduzierbar, dann ist
            A entscheidbar falls B entscheidbar.

            Beweis-Technik: Diagonalisierung (Cantor/Turing) oder Reduktion vom Halteproblem.

        @return Dictionary mit unentscheidbaren Problemen und Beweisideen
        @lastModified 2026-03-10
        """
        return {
            'halting_problem': {
                'name': 'Halteproblem (Entscheidungsproblem)',
                'statement': 'Gegeben Turing-Maschine M und Eingabe w: Hält M auf w?',
                'proof': 'Diagonalisierung (Turing, 1936): Annahme entscheidbar → Widerspruch',
                'undecidable': True,
                'semi_decidable': True,  # Hält es, kann man ja sagen; hält es nicht: unbekannt
            },
            'post_correspondence': {
                'name': "Post'sches Korrespondenzproblem (PCP)",
                'statement': 'Gegeben Wortpaare (u_i, v_i): Existiert Sequenz i_1,...,i_k mit u_{i1}...u_{ik} = v_{i1}...v_{ik}?',
                'proof': 'Reduktion vom Halteproblem (Post, 1946)',
                'undecidable': True,
            },
            'cfg_equivalence': {
                'name': 'KFG-Äquivalenz',
                'statement': 'Gegeben KFGs G1, G2: Gilt L(G1) = L(G2)?',
                'proof': 'Reduktion vom PCP',
                'undecidable': True,
            },
            'cfg_ambiguity': {
                'name': 'Ambiguität von KFGs',
                'statement': 'Gegeben KFG G: Ist G ambig?',
                'undecidable': True,
            },
            'tiling_problem': {
                'name': 'Kachelproblem (Wang-Tiles)',
                'statement': 'Gegeben Kachelset: Kann die Ebene damit gekachelt werden?',
                'undecidable': True,
            },
            'rices_theorem': {
                'name': 'Satz von Rice',
                'statement': (
                    'Jede nicht-triviale Eigenschaft der von TMs berechneten Funktionen '
                    'ist unentscheidbar.'
                ),
                'examples': [
                    'Hält M auf jeder Eingabe?',
                    'Berechnet M die leere Funktion?',
                    'Gilt L(M) = Σ*?',
                ],
                'trivial_exceptions': ['Eigenschaft gilt für ALLE oder KEINE TMs'],
            },
            'hilbert_10th': {
                'name': 'Hilberts 10. Problem (Diophantische Gleichungen)',
                'statement': 'Gegeben Diophantische Gleichung: Hat sie eine ganzzahlige Lösung?',
                'proof': 'Matiyasevich, 1970 (aufbauend auf Davis, Putnam, Robinson)',
                'undecidable': True,
            },
        }
