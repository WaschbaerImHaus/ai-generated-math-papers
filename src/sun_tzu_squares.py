"""
Sun-Tzu-Quadrate: Aufeinanderfolgende Quadrate mit vorgegebenen Resten.

Dieses Modul implementiert:
    - Chinesischer Restsatz (CRT) nach Sun-Tzu (~3. Jh. n.Chr.)
    - Legendre-Symbol und quadratische Reste mod p
    - Euler-Kriterium: a^{(p-1)/2} ≡ (a/p) (mod p)
    - Sun-Tzu-Quadrate-Vermutung #285: Aufeinanderfolgende Quadrate mit CRT

Mathematischer Hintergrund:
    Chinesischer Restsatz (THEOREM, Sun-Tzu ~3. Jh.):
        Sind m₁, ..., mₙ paarweise teilerfremd, so hat das System
            x ≡ r₁ (mod m₁), ..., x ≡ rₙ (mod mₙ)
        genau eine Lösung mod M = m₁·...·mₙ.

    Legendre-Symbol (a/p):
        (a/p) = 0  falls p | a
        (a/p) = 1  falls a ist QR mod p
        (a/p) = −1 falls a ist kein QR mod p

    Euler-Kriterium (THEOREM):
        (a/p) ≡ a^{(p-1)/2} (mod p)

    Sun-Tzu-Quadrate (CONJECTURE #285):
        Für paarweise teilerfremde Moduli m₁,...,mₙ und vorgegebene Reste
        r₁,...,rₙ mit rᵢ = quadratischer Rest mod mᵢ existieren
        n aufeinanderfolgende Quadrate a², (a+1)², ..., (a+n-1)²,
        sodass (a+i)² ≡ rᵢ (mod mᵢ) für alle i.

        Verbindung zu CRT: (a+i)² ≡ rᵢ (mod mᵢ) ⟺
            a+i ≡ ±√rᵢ (mod mᵢ) — das gibt ein System für a.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
from typing import List, Tuple, Optional, Dict
from functools import reduce


def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    Erweiterter Euklidischer Algorithmus: ax + by = gcd(a, b).

    @param a: Erste ganze Zahl
    @param b: Zweite ganze Zahl
    @return: Tripel (gcd, x, y) mit a·x + b·y = gcd
    @lastModified: 2026-03-12
    """
    if b == 0:
        return a, 1, 0
    g, x, y = _extended_gcd(b, a % b)
    return g, y, x - (a // b) * y


def _modinv(a: int, m: int) -> Optional[int]:
    """
    Berechnet das modulare Inverse von a mod m.

    Existiert nur falls gcd(a, m) = 1.

    @param a: Zahl, deren Inverses gesucht wird
    @param m: Modulus
    @return: a^{-1} mod m oder None falls nicht existent
    @lastModified: 2026-03-12
    """
    g, x, _ = _extended_gcd(a % m, m)
    if g != 1:
        return None
    return x % m


def _gcd(a: int, b: int) -> int:
    """Euklidischer Algorithmus."""
    while b:
        a, b = b, a % b
    return abs(a)


def _isqrt(n: int) -> int:
    """Ganzzahliger Quadratwurzel via Newton-Methode."""
    if n < 0:
        raise ValueError("Keine reelle Quadratwurzel für negative Zahlen")
    if n == 0:
        return 0
    x = int(math.isqrt(n))
    return x


def _is_perfect_square(n: int) -> Tuple[bool, int]:
    """
    Prüft ob n ein vollständiges Quadrat ist.

    @param n: Ganze Zahl ≥ 0
    @return: (True, sqrt(n)) falls n = k², sonst (False, 0)
    @lastModified: 2026-03-12
    """
    if n < 0:
        return False, 0
    s = _isqrt(n)
    return s * s == n, s


class ChineseRemainderTheorem:
    """
    Klassischer Chinesischer Restsatz (CRT).

    Löst simultane Kongruenzen:
        x ≡ r₀ (mod m₀)
        x ≡ r₁ (mod m₁)
        ...
        x ≡ rₙ₋₁ (mod mₙ₋₁)

    unter der Bedingung gcd(mᵢ, mⱼ) = 1 für i ≠ j.

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    @staticmethod
    def solve(remainders: List[int], moduli: List[int]) -> Optional[Tuple[int, int]]:
        """
        Löst das CRT-System und gibt die kleinste positive Lösung zurück.

        Algorithmus:
            M = Π mᵢ
            Mᵢ = M / mᵢ
            yᵢ = Mᵢ^{-1} mod mᵢ (modulares Inverses)
            x = Σ rᵢ · Mᵢ · yᵢ (mod M)

        @param remainders: Liste der Reste r₀, ..., rₙ₋₁
        @param moduli: Liste der paarweise teilerfremden Moduli m₀, ..., mₙ₋₁
        @return: (x mod M, M) oder None falls keine Lösung existiert
        @raises ValueError: Falls Moduli nicht paarweise teilerfremd
        @lastModified: 2026-03-12
        """
        if len(remainders) != len(moduli):
            raise ValueError("Längen von remainders und moduli müssen übereinstimmen")

        # Überprüfe paarweise Teilerfremdheit
        n = len(moduli)
        for i in range(n):
            for j in range(i + 1, n):
                if _gcd(moduli[i], moduli[j]) != 1:
                    raise ValueError(
                        f"Moduli {moduli[i]} und {moduli[j]} sind nicht teilerfremd "
                        f"(gcd = {_gcd(moduli[i], moduli[j])})"
                    )

        # Berechne Produkt M
        M = reduce(lambda a, b: a * b, moduli)

        x = 0
        for r, m in zip(remainders, moduli):
            Mi = M // m
            yi = _modinv(Mi, m)
            if yi is None:
                return None
            x += r * Mi * yi

        return x % M, M

    @staticmethod
    def are_coprime(moduli: List[int]) -> bool:
        """
        Prüft ob alle Moduli paarweise teilerfremd sind.

        @param moduli: Liste von ganzen Zahlen
        @return: True falls paarweise teilerfremd
        @lastModified: 2026-03-12
        """
        n = len(moduli)
        for i in range(n):
            for j in range(i + 1, n):
                if _gcd(moduli[i], moduli[j]) != 1:
                    return False
        return True


class SunTzuSquares:
    """
    Aufeinanderfolgende Quadrate mit vorgegebenen Resten (Sun-Tzu-Vermutung #285).

    Sucht ganzzahlige a, sodass:
        (a+0)² ≡ r₀ (mod m₀)
        (a+1)² ≡ r₁ (mod m₁)
        ...
        (a+n-1)² ≡ rₙ₋₁ (mod mₙ₋₁)

    Methode:
        Für jedes i: (a+i)² ≡ rᵢ (mod mᵢ)
        → a+i ≡ ±√rᵢ (mod mᵢ) (falls rᵢ quadratischer Rest mod mᵢ)
        → a ≡ (±√rᵢ - i) (mod mᵢ)
        → CRT löst simultanes System für a.

    CONJECTURE #285:
        Für paarweise teilerfremde Moduli m₁,...,mₙ und quadratische Reste
        r₁,...,rₙ existiert immer ein solches a.

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, moduli: List[int], remainders: List[int]):
        """
        Initialisiert das Sun-Tzu-Quadrate-Problem.

        @param moduli: Paarweise teilerfremde Moduli m₀, ..., mₙ₋₁
        @param remainders: Gewünschte Reste r₀, ..., rₙ₋₁
        @raises ValueError: Falls Längen nicht übereinstimmen
        @lastModified: 2026-03-12
        """
        if len(moduli) != len(remainders):
            raise ValueError("Moduli und Reste müssen gleich lang sein")
        self.moduli = moduli
        self.remainders = remainders
        self.n = len(moduli)
        self.crt = ChineseRemainderTheorem()

    def sqrt_mod_prime(self, a: int, p: int) -> List[int]:
        """
        Berechnet √a mod p für ungerade Primzahl p via Tonelli-Shanks.

        Falls a kein quadratischer Rest mod p ist, wird leere Liste zurückgegeben.

        Algorithmus (Tonelli-Shanks):
            1. Prüfe (a/p) = a^{(p-1)/2} mod p
            2. Schreibe p-1 = Q·2^S
            3. Finde Quadratnichttrest n, setze M=S, c=n^Q, t=a^Q, R=a^{(Q+1)/2}
            4. Iterate bis t ≡ 1

        @param a: Zahl unter der Wurzel
        @param p: Ungerade Primzahl
        @return: Liste der Quadratwurzeln [r, p-r] oder []
        @lastModified: 2026-03-12
        """
        a = a % p
        if a == 0:
            return [0]

        # Euler-Kriterium: a ist QR ↔ a^{(p-1)/2} ≡ 1 (mod p)
        if pow(a, (p - 1) // 2, p) != 1:
            return []  # a ist kein quadratischer Rest mod p

        # Spezialfall: p ≡ 3 (mod 4) → √a = a^{(p+1)/4} mod p
        if p % 4 == 3:
            r = pow(a, (p + 1) // 4, p)
            return sorted(set([r, p - r]))

        # Allgemeiner Fall: Tonelli-Shanks
        # Schreibe p-1 = Q·2^S
        Q = p - 1
        S = 0
        while Q % 2 == 0:
            Q //= 2
            S += 1

        # Finde Quadratnichttrest z
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1

        M = S
        c = pow(z, Q, p)
        t = pow(a, Q, p)
        R = pow(a, (Q + 1) // 2, p)

        while True:
            if t == 1:
                return sorted(set([R, p - R]))
            # Finde kleinste i > 0 mit t^{2^i} ≡ 1
            i = 1
            temp = (t * t) % p
            while temp != 1:
                temp = (temp * temp) % p
                i += 1
            b = pow(c, pow(2, M - i - 1, p - 1), p)
            M = i
            c = (b * b) % p
            t = (t * c) % p
            R = (R * b) % p

    def sqrt_mod_composite(self, a: int, m: int) -> List[int]:
        """
        Berechnet √a mod m für beliebiges m (faktorisiert m in Primpotenzen).

        Für m = p (prim): Tonelli-Shanks
        Für m = p^k: Hensels Lemma (hier: einfache Iterationsmethode)

        @param a: Zahl unter der Wurzel
        @param m: Modulus
        @return: Liste aller Quadratwurzeln mod m
        @lastModified: 2026-03-12
        """
        a = a % m
        if m == 1:
            return [0]
        if m == 2:
            return [a % 2]

        # Prüfe ob m prim
        def is_prime(n: int) -> bool:
            if n < 2:
                return False
            if n == 2:
                return True
            if n % 2 == 0:
                return False
            for d in range(3, int(n**0.5) + 1, 2):
                if n % d == 0:
                    return False
            return True

        if is_prime(m):
            return self.sqrt_mod_prime(a, m)

        # Brute-Force für kleine m
        results = []
        for x in range(m):
            if (x * x) % m == a:
                results.append(x)
        return results

    def legendre_symbol(self, a: int, p: int) -> int:
        """
        Berechnet das Legendre-Symbol (a/p) via Euler-Kriterium.

        (a/p) = 0 falls p | a
        (a/p) = 1 falls a QR mod p
        (a/p) = -1 falls a kein QR mod p

        THEOREM (Euler-Kriterium): (a/p) ≡ a^{(p-1)/2} (mod p)

        @param a: Zahl
        @param p: Ungerade Primzahl
        @return: Legendre-Symbol ∈ {-1, 0, 1}
        @lastModified: 2026-03-12
        """
        if p < 2:
            raise ValueError(f"p={p} ist keine gültige Primzahl")
        a = a % p
        if a == 0:
            return 0
        result = pow(a, (p - 1) // 2, p)
        # Ergebnis ist 1 oder p-1 ≡ -1 (mod p)
        return 1 if result == 1 else -1

    def find_a(self, search_bound: int = 10000) -> Optional[int]:
        """
        Sucht das kleinste positive a, das die aufeinanderfolgenden Quadrat-Kongruenzen erfüllt.

        Methode:
            1. Für jedes i: berechne √rᵢ mod mᵢ
            2. Für jede Kombination von Vorzeichen: löse CRT für a
            3. Prüfe, ob a alle Kongruenzen erfüllt

        @param search_bound: Obere Schranke für a (nur für Brute-Force-Fallback)
        @return: Kleinstes nicht-negatives a oder None falls nicht gefunden
        @lastModified: 2026-03-12
        """
        # Berechne mögliche Quadratwurzeln für jedes rᵢ mod mᵢ
        sqrt_options: List[List[int]] = []
        for i, (m, r) in enumerate(zip(self.moduli, self.remainders)):
            roots = self.sqrt_mod_composite(r, m)
            if not roots:
                return None  # rᵢ ist kein QR mod mᵢ → keine Lösung möglich
            # a ≡ (root - i) mod mᵢ
            a_options = [(root - i) % m for root in roots]
            sqrt_options.append(list(set(a_options)))

        # Erzeuge alle Kombinationen und löse CRT
        from itertools import product as iproduct
        best_a = None

        for combination in iproduct(*sqrt_options):
            crt_result = ChineseRemainderTheorem.solve(list(combination), self.moduli)
            if crt_result is None:
                continue
            a_val, M = crt_result

            # Prüfe Lösung
            valid = True
            for i in range(self.n):
                if ((a_val + i) ** 2) % self.moduli[i] != self.remainders[i] % self.moduli[i]:
                    valid = False
                    break

            if valid:
                if best_a is None or a_val < best_a:
                    best_a = a_val

        return best_a

    def verify_solution(self, a: int) -> bool:
        """
        Verifiziert, dass a alle Kongruenzen (a+i)² ≡ rᵢ (mod mᵢ) erfüllt.

        @param a: Zu prüfender Wert
        @return: True wenn alle Kongruenzen erfüllt
        @lastModified: 2026-03-12
        """
        for i in range(self.n):
            if ((a + i) ** 2) % self.moduli[i] != self.remainders[i] % self.moduli[i]:
                return False
        return True

    def brute_force_search(self, a_max: int = 100000) -> Optional[int]:
        """
        Brute-Force-Suche nach kleinstem a ≥ 0 mit (a+i)² ≡ rᵢ (mod mᵢ).

        @param a_max: Suchbereich für a
        @return: Kleinstes a oder None
        @lastModified: 2026-03-12
        """
        for a in range(a_max):
            if self.verify_solution(a):
                return a
        return None

    def quadratic_residues_mod_p(self, p: int) -> List[int]:
        """
        Berechnet alle quadratischen Reste mod p.

        QR(p) = {a² mod p : a = 0, 1, ..., p-1}

        @param p: Primzahl
        @return: Sortierte Liste der quadratischen Reste
        @lastModified: 2026-03-12
        """
        return sorted(set((a * a) % p for a in range(p)))

    def crt_connection_explanation(self) -> str:
        """
        Erklärt die Verbindung zwischen Sun-Tzu-Quadraten und dem CRT.

        @return: Erklärungstext
        @lastModified: 2026-03-12
        """
        lines = [
            "Sun-Tzu (Chinese Remainder Theorem) Verbindung:",
            "",
            "Gesucht: a mit (a+i)² ≡ rᵢ (mod mᵢ), i=0,...,n-1",
            "",
            "Schritt 1: Quadratwurzeln berechnen",
            "   (a+i)² ≡ rᵢ (mod mᵢ)",
            "   ↔ a+i ≡ ±√rᵢ (mod mᵢ)  (wenn rᵢ ist QR mod mᵢ)",
            "   ↔ a ≡ (±√rᵢ - i) (mod mᵢ)",
            "",
            "Schritt 2: CRT löst das simultane System",
            "   a ≡ c₀ (mod m₀), a ≡ c₁ (mod m₁), ...",
            "   (eindeutige Lösung mod M = m₀·m₁·...·mₙ₋₁)",
            "",
            "CONJECTURE #285: Für paarweise teilerfremde Moduli existiert",
            "immer ein solches a (wenn die rᵢ QR modulo mᵢ sind).",
            "",
            f"Aktuelles Problem: Moduli={self.moduli}, Reste={self.remainders}"
        ]
        return "\n".join(lines)


class QuadraticResidueAnalysis:
    """
    Analyse quadratischer Reste mod p und verwandte Strukturen.

    Quadratische Reste mod p (p ungerade Primzahl):
        QR_p = {1², 2², ..., ((p-1)/2)²} mod p
        |QR_p| = (p-1)/2

    Gauss-Lemma (THEOREM):
        (a/p) = (-1)^n, wobei n = |{k ∈ {1,...,(p-1)/2} : k·a mod p > p/2}|

    Quadratisches Reziprozitätsgesetz (THEOREM, Gauss 1796):
        (p/q)·(q/p) = (-1)^{(p-1)/2·(q-1)/2} für ungerade Primzahlen p ≠ q

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    @staticmethod
    def all_qr(p: int) -> List[int]:
        """
        Alle quadratischen Reste mod p (ohne 0).

        @param p: Ungerade Primzahl
        @return: Sortierte Liste der QR mod p
        @lastModified: 2026-03-12
        """
        return sorted(set((a * a) % p for a in range(1, p)))

    @staticmethod
    def euler_criterion(a: int, p: int) -> int:
        """
        Berechnet das Legendre-Symbol via Euler-Kriterium.

        THEOREM: a^{(p-1)/2} ≡ (a/p) (mod p)

        @param a: Ganze Zahl
        @param p: Ungerade Primzahl
        @return: 1, -1 oder 0
        @lastModified: 2026-03-12
        """
        a = a % p
        if a == 0:
            return 0
        result = pow(a, (p - 1) // 2, p)
        return 1 if result == 1 else -1

    @staticmethod
    def quadratic_reciprocity(p: int, q: int) -> Dict[str, object]:
        """
        Berechnet das quadratische Reziprozitätsgesetz für ungerade Primzahlen p, q.

        THEOREM (Gauss 1796):
            (p/q)·(q/p) = (-1)^{(p-1)/2·(q-1)/2}

        Speziell:
            Falls p ≡ 1 (mod 4) oder q ≡ 1 (mod 4): (p/q) = (q/p)
            Falls p ≡ q ≡ 3 (mod 4): (p/q) = −(q/p)

        @param p: Erste ungerade Primzahl
        @param q: Zweite ungerade Primzahl (p ≠ q)
        @return: Dictionary mit Symbolen und Reziprozitätsergebnis
        @lastModified: 2026-03-12
        """
        leg_pq = QuadraticResidueAnalysis.euler_criterion(p, q)
        leg_qp = QuadraticResidueAnalysis.euler_criterion(q, p)
        exponent = ((p - 1) // 2) * ((q - 1) // 2)
        expected_product = (-1) ** exponent

        return {
            "p": p, "q": q,
            "legendre_p_mod_q": leg_pq,
            "legendre_q_mod_p": leg_qp,
            "product": leg_pq * leg_qp,
            "expected_(-1)^exp": expected_product,
            "reciprocity_holds": leg_pq * leg_qp == expected_product,
            "theorem": "Quadratisches Reziprozitätsgesetz (Gauss 1796)"
        }
