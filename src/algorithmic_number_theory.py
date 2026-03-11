"""
Algorithmische Zahlentheorie - Siebe, Primzahltests, Faktorisierung,
diskreter Logarithmus, Kettenbrüche und kryptographische Algorithmen.

Dieses Modul implementiert klassische und moderne Algorithmen der algorithmischen
Zahlentheorie. Es dient als mathematische Grundlage für Kryptographie, Primzahl-
forschung und zahlentheoretische Berechnungen.

Wichtig: miller_rabin_primality_test, sieve_of_eratosthenes, legendre_symbol und
jacobi_symbol stammen aus proof_theory.py und werden dort importiert, um Duplikate
zu vermeiden.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import math
import random
from typing import Optional
from proof_theory import (
    miller_rabin_primality_test,
    sieve_of_eratosthenes as _sieve_eratosthenes,
    legendre_symbol as _legendre_symbol,
    jacobi_symbol,
)


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _gcd(a: int, b: int) -> int:
    """
    Berechnet den größten gemeinsamen Teiler via euklidischem Algorithmus.

    @param a: Erste ganze Zahl
    @param b: Zweite ganze Zahl
    @return: ggT(a, b)
    @lastModified: 2026-03-10
    """
    while b:
        a, b = b, a % b
    return abs(a)


def _mod_inverse(a: int, m: int) -> Optional[int]:
    """
    Berechnet das modulare Inverse von a modulo m via erweitertem euklidischen Algo.

    @param a: Zahl, deren Inverses gesucht wird
    @param m: Modulus
    @return: x mit a*x ≡ 1 (mod m), oder None wenn kein Inverses existiert
    @lastModified: 2026-03-10
    """
    if _gcd(a, m) != 1:
        return None
    # Erweiterter euklidischer Algorithmus
    old_r, r = a % m, m
    old_s, s = 1, 0
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    return old_s % m


def _is_perfect_square(n: int) -> tuple[bool, int]:
    """
    Prüft ob n eine perfekte Quadratzahl ist.

    @param n: Nicht-negative ganze Zahl
    @return: (True, Wurzel) wenn n = k², (False, 0) sonst
    @lastModified: 2026-03-10
    """
    if n < 0:
        return False, 0
    root = math.isqrt(n)
    if root * root == n:
        return True, root
    return False, 0


# ===========================================================================
# 1. PRIMZAHLSIEBE
# ===========================================================================

class PrimeSieve:
    """
    Klasse mit verschiedenen Siebverfahren zur Primzahlberechnung.

    Enthält das klassische Eratosthenes-Sieb, das moderne Atkin-Sieb,
    sowie Funktionen zur Analyse von Primzahllücken und Zwillingsprimzahlen.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def sieve_of_eratosthenes(n: int) -> list[int]:
        """
        Sieb des Eratosthenes - findet alle Primzahlen bis n.

        Delegiert an die Implementierung in proof_theory.py, um Duplikate
        zu vermeiden. Zeitkomplexität: O(n log log n).

        @param n: Obere Grenze (inklusiv)
        @return: Sortierte Liste aller Primzahlen ≤ n
        @lastModified: 2026-03-10
        """
        return _sieve_eratosthenes(n)

    @staticmethod
    def sieve_of_atkin(n: int) -> list[int]:
        """
        Sieb von Atkin - modernes Primzahlsieb, effizienter für große n.

        Der Atkin-Algorithmus (2003) hat Zeitkomplexität O(n / log log n)
        im Vergleich zu O(n log log n) beim Eratosthenes-Sieb.

        Grundprinzip: Basierend auf quadratischen Formen:
        - 4x² + y² ≡ 1 (mod 4)  → n ist möglicherweise prim
        - 3x² + y² ≡ 7 (mod 12) → n ist möglicherweise prim
        - 3x² - y² ≡ 11 (mod 12) → n ist möglicherweise prim
        Dann werden Quadratzahlvielfache ausgesiebt.

        @param n: Obere Grenze (inklusiv)
        @return: Sortierte Liste aller Primzahlen ≤ n
        @lastModified: 2026-03-10
        """
        if n < 2:
            return []
        if n == 2:
            return [2]
        if n == 3:
            return [2, 3]

        # Initialisierung: alle auf False setzen
        is_prime = [False] * (n + 1)

        # Schritt 1: Quadratische Formen auswerten
        sqrt_n = math.isqrt(n)
        for x in range(1, sqrt_n + 1):
            for y in range(1, sqrt_n + 1):
                # Form 1: 4x² + y²
                val = 4 * x * x + y * y
                if val <= n and val % 12 in (1, 5):
                    is_prime[val] = not is_prime[val]

                # Form 2: 3x² + y²
                val = 3 * x * x + y * y
                if val <= n and val % 12 == 7:
                    is_prime[val] = not is_prime[val]

                # Form 3: 3x² - y² (nur wenn x > y)
                if x > y:
                    val = 3 * x * x - y * y
                    if val <= n and val % 12 == 11:
                        is_prime[val] = not is_prime[val]

        # Schritt 2: Quadratzahlvielfache aussieben
        for r in range(5, sqrt_n + 1):
            if is_prime[r]:
                # Alle Vielfachen von r² als nicht-prim markieren
                r_sq = r * r
                for k in range(r_sq, n + 1, r_sq):
                    is_prime[k] = False

        # Schritt 3: 2 und 3 manuell hinzufügen, Rest sammeln
        primes = [2, 3]
        primes.extend(i for i in range(5, n + 1) if is_prime[i])
        return primes

    @staticmethod
    def prime_gaps(n: int) -> list[tuple[int, int, int]]:
        """
        Berechnet alle Primzahllücken bis n.

        Eine Primzahllücke ist die Differenz zwischen aufeinanderfolgenden Primzahlen.
        Die erste Primzahllücke ist g(2) = 3 - 2 = 1.

        @param n: Obere Grenze (inklusiv)
        @return: Liste von (p1, p2, gap) für aufeinanderfolgende Primzahlpaare
        @lastModified: 2026-03-10
        """
        primes = _sieve_eratosthenes(n)
        gaps = []
        # Lücken zwischen je zwei aufeinanderfolgenden Primzahlen berechnen
        for i in range(len(primes) - 1):
            p1 = primes[i]
            p2 = primes[i + 1]
            gaps.append((p1, p2, p2 - p1))
        return gaps

    @staticmethod
    def twin_primes(n: int) -> list[tuple[int, int]]:
        """
        Findet alle Zwillingsprimzahlpaare (p, p+2) mit p ≤ n.

        Zwillingsprimzahlen sind Primzahlpaare mit Abstand 2.
        Die Zwillingsprimzahlvermutung besagt, dass es unendlich viele gibt.

        Beispiele: (3,5), (5,7), (11,13), (17,19), (29,31), ...

        @param n: Obere Grenze für die erste Primzahl des Paares
        @return: Liste aller Zwillingsprimzahlpaare (p, p+2) mit p ≤ n
        @lastModified: 2026-03-10
        """
        if n < 3:
            return []
        # Sieb bis n+2, damit auch (n, n+2) erfasst wird
        primes_set = set(_sieve_eratosthenes(n + 2))
        primes = _sieve_eratosthenes(n)
        twins = []
        for p in primes:
            if (p + 2) in primes_set:
                twins.append((p, p + 2))
        return twins

    @staticmethod
    def prime_counting_exact(n: int) -> int:
        """
        Berechnet π(n) - die exakte Anzahl der Primzahlen ≤ n via Sieb.

        Die Primzahlzählfunktion π(x) ist ein zentrales Objekt der analytischen
        Zahlentheorie. Der Primzahlsatz besagt: π(x) ~ x / ln(x).

        @param n: Obere Grenze
        @return: π(n) = exakte Anzahl der Primzahlen ≤ n
        @lastModified: 2026-03-10
        """
        return len(_sieve_eratosthenes(n))


# ===========================================================================
# 2. PRIMZAHLTESTS
# ===========================================================================

class PrimalityTests:
    """
    Sammlung von Primzahltests: deterministisch und probabilistisch.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def aks_test_demo(n: int) -> dict:
        """
        AKS-Primzahltest - Demonstration für kleine n.

        Der AKS-Test (Agrawal, Kayal, Saxena, 2002) ist der erste deterministische
        Primzahltest in Polynomialzeit: O(log(n)^{12}) (ursprünglich), verbessert
        auf O(log(n)^{6}).

        Kernidee: n ist genau dann prim, wenn für ein geeignetes r gilt:
            (X + a)^n ≡ X^n + a (mod X^r - 1, n)  für alle a ≤ √φ(r) · log(n)

        Diese Demo-Implementierung prüft zunächst einfache Fälle und nutzt
        für kleine n (< 1000) eine vereinfachte Version des Perfekte-Potenz-Tests
        und anschließend den Miller-Rabin-Test deterministisch.

        @param n: Zu prüfende Zahl (n ≥ 2, für Demo: n ≤ 10000 empfohlen)
        @return: Dict mit 'is_prime', 'method', 'steps'
        @lastModified: 2026-03-10
        """
        steps = []

        if n < 2:
            return {'is_prime': False, 'method': 'AKS-Demo', 'steps': ['n < 2: nicht prim']}

        # Schritt 1: Ist n eine perfekte Potenz? (n = a^b mit b ≥ 2)
        steps.append(f"Schritt 1: Prüfe ob {n} eine perfekte Potenz ist")
        for b in range(2, math.floor(math.log2(n)) + 2 if n > 1 else 2):
            a = round(n ** (1.0 / b))
            for candidate in [a - 1, a, a + 1]:
                if candidate >= 2 and candidate ** b == n:
                    steps.append(f"  {n} = {candidate}^{b} → zusammengesetzt")
                    return {'is_prime': False, 'method': 'AKS-Demo', 'steps': steps}
        steps.append(f"  Kein perfekte Potenz gefunden")

        # Schritt 2: Für Demo: nutze deterministischen Miller-Rabin
        steps.append(f"Schritt 2: Miller-Rabin (deterministisch für n < 3.3×10²⁴)")
        result = miller_rabin_primality_test(n)
        steps.append(f"  Miller-Rabin: {'prim' if result else 'zusammengesetzt'}")

        return {
            'is_prime': result,
            'method': 'AKS-Demo (vereinfacht)',
            'steps': steps
        }

    @staticmethod
    def solovay_strassen(n: int, k: int = 10) -> bool:
        """
        Solovay-Strassen-Primzahltest (probabilistisch).

        Fehlerwahrscheinlichkeit: ≤ (1/2)^k pro Runde.
        Grundlage: Euler-Kriterium via Jacobi-Symbol:
            a^((n-1)/2) ≡ (a/n) (mod n)  für alle a coprim zu n

        Euler-Pseudoprimzahlen täuschen diesen Test, aber mit Wahrscheinlichkeit ≤ 1/2
        pro Basis. Im Gegensatz zu Fermat-Pseudoprimzahlen gibt es keine analogen
        universellen Gegenbeispiele (Carmichael-Zahlen schlagen HIER NICHT durch).

        @param n: Zu prüfende Zahl (n ≥ 2)
        @param k: Anzahl der Testrunden (Standard: 10)
        @return: True wenn wahrscheinlich prim, False wenn sicher zusammengesetzt
        @lastModified: 2026-03-10
        """
        if n < 2:
            return False
        if n == 2:
            return True
        if n % 2 == 0:
            return False

        for _ in range(k):
            # Zufällige Basis a im Bereich [2, n-1]
            a = random.randint(2, n - 1)

            # Wenn ggT(a, n) ≠ 1 → n ist zusammengesetzt
            if _gcd(a, n) != 1:
                return False

            # Jacobi-Symbol (a/n)
            j = jacobi_symbol(a, n)

            # Euler-Kriterium: a^((n-1)/2) mod n muss gleich j mod n sein
            euler = pow(a, (n - 1) // 2, n)
            j_mod_n = j % n  # Jacobi-Symbol kann -1 sein → mod n macht es positiv

            if euler != j_mod_n:
                return False  # Definitiv zusammengesetzt

        return True  # Wahrscheinlich prim

    @staticmethod
    def fermat_test(n: int, a: int) -> bool:
        """
        Fermatscher Primzahltest: prüft ob a^(n-1) ≡ 1 (mod n).

        Kleiner Satz von Fermat: Wenn p prim ist und p∤a, dann gilt a^(p-1) ≡ 1 (mod p).
        Der Umkehrschluss gilt NICHT allgemein (Carmichael-Zahlen als Gegenbeispiel).

        VORSICHT: Carmichael-Zahlen (z.B. 561, 1105) täuschen diesen Test für ALLE
        Basen a mit ggT(a, n) = 1!

        @param n: Zu testende Zahl (n ≥ 2)
        @param a: Fermat-Zeuge/Basis (2 ≤ a ≤ n-2 empfohlen)
        @return: True wenn n den Fermat-Test zur Basis a besteht (n "Fermat-prim zur Basis a")
        @lastModified: 2026-03-10
        """
        if n < 2:
            return False
        if n == 2:
            return True
        if n % 2 == 0:
            return False
        if _gcd(a, n) != 1:
            # Falls ggT(a,n) ≠ 1, ist n sicher zusammengesetzt
            return False
        # Fermatscher Test: a^(n-1) ≡ 1 (mod n)?
        return pow(a, n - 1, n) == 1

    @staticmethod
    def lucas_primality_test(n: int) -> dict:
        """
        Lucas-Primzahltest (vereinfacht).

        Ein Lucas-Primzahltest beweist die Primalität von n, wenn gilt:
        Es existiert eine Zahl a mit:
          1. a^(n-1) ≡ 1 (mod n)
          2. a^((n-1)/q) ≢ 1 (mod n) für alle Primteiler q von n-1

        Das heißt: a ist eine primitive Wurzel modulo n, was nur möglich ist,
        wenn n prim ist.

        @param n: Zu prüfende Zahl (n ≥ 2)
        @return: Dict mit 'is_prime', 'witness', 'factors_of_n_minus_1'
        @lastModified: 2026-03-10
        """
        if n < 2:
            return {'is_prime': False, 'witness': None, 'factors_of_n_minus_1': []}
        if n == 2:
            return {'is_prime': True, 'witness': 1, 'factors_of_n_minus_1': []}
        if n % 2 == 0:
            return {'is_prime': False, 'witness': None, 'factors_of_n_minus_1': []}

        # Primfaktorisierung von n-1
        n_minus_1 = n - 1
        factors = _prime_factors_list(n_minus_1)
        unique_factors = list(set(factors))

        # Suche Basis a, die als Lucas-Zeuge dient
        for a in range(2, min(n, 100)):
            if _gcd(a, n) != 1:
                continue

            # Bedingung 1: a^(n-1) ≡ 1 (mod n)
            if pow(a, n_minus_1, n) != 1:
                continue

            # Bedingung 2: a^((n-1)/q) ≢ 1 (mod n) für alle Primteiler q von n-1
            is_primitive = True
            for q in unique_factors:
                if pow(a, n_minus_1 // q, n) == 1:
                    is_primitive = False
                    break

            if is_primitive:
                return {
                    'is_prime': True,
                    'witness': a,
                    'factors_of_n_minus_1': unique_factors
                }

        # Kein Zeuge gefunden → wahrscheinlich nicht prim
        # (Für wirklich große n kann die Suche bis n-1 gehen, hier vereinfacht)
        return {
            'is_prime': miller_rabin_primality_test(n),
            'witness': None,
            'factors_of_n_minus_1': unique_factors
        }

    @staticmethod
    def is_carmichael_number(n: int) -> bool:
        """
        Prüft ob n eine Carmichael-Zahl ist.

        Carmichael-Zahlen (auch absolute Fermat-Pseudoprimzahlen genannt) sind
        zusammengesetzte Zahlen n, für die gilt:
            a^n ≡ a (mod n) für alle ganzen Zahlen a

        Äquivalent (Korselt-Kriterium):
          1. n ist quadratfrei (kein Primzahlquadrat teilt n)
          2. Für jeden Primteiler p von n gilt: (p-1) | (n-1)

        Kleinste Carmichael-Zahl: 561 = 3 × 11 × 17

        @param n: Zu prüfende Zahl
        @return: True wenn n eine Carmichael-Zahl ist
        @lastModified: 2026-03-10
        """
        # Carmichael-Zahlen sind zusammengesetzt
        if n < 2 or miller_rabin_primality_test(n):
            return False

        # Primfaktorisierung von n
        factors = _prime_factors_list(n)

        # Korselt-Kriterium: n muss quadratfrei sein
        if len(factors) != len(set(factors)):
            return False  # Quadratbehaftet → keine Carmichael-Zahl

        # Korselt-Kriterium: (p-1) | (n-1) für alle Primteiler p
        n_minus_1 = n - 1
        for p in set(factors):
            if n_minus_1 % (p - 1) != 0:
                return False

        return True


def _prime_factors_list(n: int) -> list[int]:
    """
    Gibt die Primfaktoren von n als Liste zurück (mit Wiederholung).

    Interne Hilfsfunktion für Primfaktorisierung.

    @param n: Positive ganze Zahl ≥ 2
    @return: Liste aller Primfaktoren (mit Vielfachheit)
    @lastModified: 2026-03-10
    """
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


# ===========================================================================
# 3. FAKTORISIERUNGSALGORITHMEN
# ===========================================================================

class IntegerFactorization:
    """
    Verschiedene Algorithmen zur Faktorisierung ganzer Zahlen.

    Von einfacher Probedivision bis zu modernen Algorithmen wie
    Pollard's ρ-Algorithmus und Fermat-Faktorisierung.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def trial_division(n: int) -> dict:
        """
        Faktorisierung durch Probedivision bis √n.

        Einfachster Faktorisierungsalgorithmus: teste alle möglichen Teiler
        von 2 bis √n. Zeitkomplexität: O(√n).

        Effizient nur für kleine n (n < 10^12). Für größere Zahlen sind
        spezialisiertere Algorithmen nötig.

        @param n: Zu faktorisierende Zahl (n ≥ 1)
        @return: Dict mit 'factors' (Primfaktoren mit Vielfachheit), 'is_prime'
        @lastModified: 2026-03-10
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")
        if n == 1:
            return {'factors': {}, 'is_prime': False}

        factors = {}
        remaining = n

        # Teste alle Teiler von 2 bis √remaining
        d = 2
        while d * d <= remaining:
            while remaining % d == 0:
                factors[d] = factors.get(d, 0) + 1
                remaining //= d
            d += 1

        # Falls remaining > 1, ist es ein Primteiler
        if remaining > 1:
            factors[remaining] = factors.get(remaining, 0) + 1

        is_prime = (len(factors) == 1 and list(factors.values())[0] == 1)
        return {'factors': factors, 'is_prime': is_prime}

    @staticmethod
    def pollard_rho(n: int, seed: int = 2) -> Optional[int]:
        """
        Pollard's ρ-Algorithmus (Floyd-Zykluserkennung) für Faktorisierung.

        Findet einen nicht-trivialen Teiler von n.
        Erwartet O(n^{1/4}) Schritte dank Birthday-Paradoxon.

        Idee: Iteriere f(x) = x² + c (mod n), erzeuge eine zufällig-wirkende
        Folge. Zykluserkennung nach Floyd: "Hase" und "Schildkröte" laufen
        mit doppelter Geschwindigkeit. Wenn ggT(|x_hase - x_schildkröte|, n) > 1
        → Teiler gefunden.

        @param n: Zu faktorisierende Zahl (n ≥ 4, n zusammengesetzt)
        @param seed: Startwert (Standard: 2)
        @return: Ein nicht-trivialer Teiler von n, oder None wenn fehlgeschlagen
        @lastModified: 2026-03-10
        """
        if n < 4:
            return None
        if n % 2 == 0:
            return 2

        # Versuche mehrere c-Werte falls nötig
        for c in range(1, 20):
            x = seed % n
            y = seed % n
            d = 1

            # Floyd's Zyklus-Erkennung
            while d == 1:
                # Schildkröte: ein Schritt
                x = (x * x + c) % n
                # Hase: zwei Schritte
                y = (y * y + c) % n
                y = (y * y + c) % n
                # ggT berechnen
                d = _gcd(abs(x - y), n)

            if d != n:
                return d  # Nicht-trivialer Teiler gefunden

        return None  # Fehlgeschlagen (sehr selten)

    @staticmethod
    def pollard_p_minus_1(n: int, B: int = 100) -> Optional[int]:
        """
        Pollard's p-1 Algorithmus zur Faktorisierung.

        Effektiv wenn n einen Primteiler p hat, sodass p-1 B-glatt ist,
        d.h. alle Primfaktoren von p-1 sind ≤ B.

        Idee (Fermatscher kleiner Satz): Wenn p prim und p | n,
        dann gilt a^(p-1) ≡ 1 (mod p). Setze M = LCM(1,...,B).
        Falls p-1 | M, dann p | (a^M - 1), also p | ggT(a^M - 1, n).

        @param n: Zu faktorisierende Zahl
        @param B: Glattheitsg-renze (B-smooth-Schranke)
        @return: Nicht-trivialer Teiler von n, oder None
        @lastModified: 2026-03-10
        """
        if n < 4:
            return None
        if n % 2 == 0:
            return 2

        # Berechne a = 2^(∏ Primpotenzen ≤ B) mod n
        a = 2
        primes = _sieve_eratosthenes(B)

        for p in primes:
            # Maximale Potenz von p ≤ B
            pk = p
            while pk * p <= B:
                pk *= p
            # a = a^pk mod n
            a = pow(a, pk, n)

        # ggT(a-1, n) berechnen
        d = _gcd(a - 1, n)
        if 1 < d < n:
            return d

        return None  # Kein Teiler gefunden

    @staticmethod
    def fermat_factoring(n: int) -> Optional[tuple[int, int]]:
        """
        Fermat-Faktorisierung über die Darstellung n = a² - b².

        Funktioniert gut wenn n zwei ähnlich große Primteiler hat.
        Idee: n = (a+b)(a-b) → suche a ≥ √n mit a² - n = b² perfekte Quadratzahl.

        Zeitkomplexität: O((p-q)²/8 · 1/√n) - effizient wenn p ≈ q.

        @param n: Zu faktorisierende ungerade Zahl (n ≥ 3)
        @return: Paar (p, q) mit n = p·q, oder None wenn n prim/keine Faktorisierung
        @lastModified: 2026-03-10
        """
        if n < 3 or n % 2 == 0:
            return None
        if miller_rabin_primality_test(n):
            return None  # n ist prim

        # Starte bei a = ceil(√n)
        a = math.isqrt(n)
        if a * a < n:
            a += 1

        # Maximale Anzahl Versuche (für praktische Nutzung begrenzt)
        max_attempts = min(100000, n)

        for _ in range(max_attempts):
            b_sq = a * a - n
            is_sq, b = _is_perfect_square(b_sq)
            if is_sq:
                p, q = a - b, a + b
                if p > 1 and q > 1:
                    return (p, q)
            a += 1

        return None  # Kein Teiler gefunden innerhalb Grenzen

    @staticmethod
    def factor_completely(n: int) -> dict:
        """
        Vollständige Primfaktorisierung durch Kombination aller Methoden.

        Strategie:
        1. Probedivision bis min(√n, 1000) für kleine Faktoren
        2. Miller-Rabin-Primalitätstest für übrig gebliebene Faktoren
        3. Pollard's ρ für mittelgroße zusammengesetzte Faktoren
        4. Fermat-Faktorisierung als Fallback

        @param n: Positive ganze Zahl ≥ 2
        @return: Dict {'factors': {p: exp}, 'factorization_str': '...'}
        @lastModified: 2026-03-10
        """
        if n < 2:
            raise ValueError(f"n muss ≥ 2 sein, erhalten: {n}")

        factors = {}
        remaining = n

        def add_factor(p: int):
            """Primfaktor p zur Faktorisierung hinzufügen."""
            factors[p] = factors.get(p, 0) + 1

        def factorize_recursive(m: int):
            """Rekursive vollständige Faktorisierung."""
            if m == 1:
                return
            if miller_rabin_primality_test(m):
                add_factor(m)
                return

            # Methode 1: Kleine Faktoren via Probedivision
            for p in range(2, min(1000, m)):
                if m % p == 0:
                    add_factor(p)
                    factorize_recursive(m // p)
                    return

            # Methode 2: Pollard's ρ
            d = IntegerFactorization.pollard_rho(m)
            if d is not None and 1 < d < m:
                factorize_recursive(d)
                factorize_recursive(m // d)
                return

            # Methode 3: Fermat-Faktorisierung
            result = IntegerFactorization.fermat_factoring(m)
            if result is not None:
                factorize_recursive(result[0])
                factorize_recursive(result[1])
                return

            # Fallback: m als "wahrscheinlich prim" behandeln
            add_factor(m)

        factorize_recursive(remaining)

        # Faktorisierungsstring erzeugen
        parts = [f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factors.items())]
        factorization_str = " × ".join(parts)

        return {'factors': factors, 'factorization_str': factorization_str}

    @staticmethod
    def smooth_number_check(n: int, B: int) -> bool:
        """
        Prüft ob n B-glatt ist (alle Primfaktoren ≤ B).

        Eine B-glatte Zahl hat ausschließlich Primfaktoren ≤ B.
        B-glatte Zahlen spielen eine wichtige Rolle im Quadratischen Sieb
        und im Indexkalkül-Algorithmus.

        Beispiel: 360 = 2³ × 3² × 5 ist 5-glatt.

        @param n: Zu prüfende Zahl (n ≥ 2)
        @param B: Glattheitsg-renze
        @return: True wenn n B-glatt ist
        @lastModified: 2026-03-10
        """
        if n < 2:
            return False
        remaining = n
        # Alle Primfaktoren ≤ B herausteilen
        for p in _sieve_eratosthenes(B):
            while remaining % p == 0:
                remaining //= p
        # Wenn remaining = 1: alle Faktoren waren ≤ B → B-glatt
        return remaining == 1


# ===========================================================================
# 4. DISKRETER LOGARITHMUS
# ===========================================================================

class DiscreteLogarithm:
    """
    Algorithmen zur Berechnung des diskreten Logarithmus g^x ≡ h (mod p).

    Der diskrete Logarithmus ist die Grundlage vieler kryptographischer Systeme
    (Diffie-Hellman, ElGamal, DSA). Die Schwierigkeit seiner Berechnung ist
    ein offenes Problem der Komplexitätstheorie.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def baby_step_giant_step(g: int, h: int, p: int) -> Optional[int]:
        """
        Baby-Step-Giant-Step (BSGS) Algorithmus für diskreten Logarithmus.

        Findet x mit g^x ≡ h (mod p) in O(√p) Zeit und Speicher.

        Idee: Schreibe x = i·m - j mit m = ⌈√p⌉:
        - Baby Steps: Berechne h·g^j für j = 0,...,m-1 (Tabelle)
        - Giant Steps: Berechne (g^m)^i für i = 1,...,m
        - Suche Kollision: g^(im) = h·g^j → x = im - j

        @param g: Basis (Generator)
        @param h: Zielwert
        @param p: Primzahlmodulus
        @return: x mit g^x ≡ h (mod p), oder None wenn nicht gefunden
        @lastModified: 2026-03-10
        """
        if p <= 1:
            return None

        # Trivialer Fall: h = 1 → x = 0
        if h % p == 1:
            return 0

        # Schrittweite m = ⌈√(p-1)⌉
        m = math.isqrt(p - 1) + 1

        # Baby Steps: Tabelle {h·g^j mod p: j} für j = 0,...,m-1
        baby_steps = {}
        g_j = 1  # g^0 = 1
        for j in range(m):
            # Speichere h·g^j mod p mit Index j
            key = (h * g_j) % p
            baby_steps[key] = j
            g_j = (g_j * g) % p  # g^(j+1)

        # Giant Steps: berechne g^(m·i) für i = 1,...,m
        g_m = pow(g, m, p)  # g^m mod p
        g_mi = g_m  # g^(m·1) = g^m

        for i in range(1, m + 1):
            # Suche ob g^(mi) in Baby-Steps-Tabelle steht
            if g_mi in baby_steps:
                j = baby_steps[g_mi]
                x = (i * m - j) % (p - 1)
                # Verifikation
                if pow(g, x, p) == h % p:
                    return x
            g_mi = (g_mi * g_m) % p  # g^(m·(i+1))

        return None  # Kein Logarithmus gefunden

    @staticmethod
    def pohlig_hellman(g: int, h: int, p: int) -> Optional[int]:
        """
        Pohlig-Hellman-Algorithmus für diskreten Logarithmus.

        Effizient wenn p-1 nur kleine Primfaktoren hat (B-glatt ist).
        Nutzt Chinesischen Restsatz (CRT) und Baby-Step-Giant-Step.

        Wenn ord(g) = q₁^e₁ · q₂^e₂ · ... (Primfaktorisierung von p-1),
        berechne x mod qᵢ^eᵢ für jeden Primfaktor separat,
        dann kombiniere via CRT.

        @param g: Generator/Basis
        @param h: Zielwert
        @param p: Primzahlmodulus
        @return: x mit g^x ≡ h (mod p), oder None
        @lastModified: 2026-03-10
        """
        if p <= 1:
            return None

        order = p - 1  # Ordnung der Gruppe Z_p*

        # Primfaktorisierung von order
        factors_dict = IntegerFactorization.trial_division(order)['factors']

        # Für jeden Primpotenzteil der Gruppenordnung
        remainders = []
        moduli = []

        for q, e in factors_dict.items():
            q_e = q ** e

            # Berechne x mod q^e via Baby-Step-Giant-Step in der q^e-Untergruppe
            # g_q = g^(order/q^e) ist ein Element der Ordnung q^e
            g_q = pow(g, order // q_e, p)
            h_q = pow(h, order // q_e, p)

            # Löse diskreten Logarithmus in der Untergruppe der Ordnung q^e
            x_part = DiscreteLogarithm.baby_step_giant_step(g_q, h_q, q_e + 1)

            if x_part is None:
                # Fallback: Brute Force für kleine q^e
                x_part = 0
                curr = 1
                for k in range(q_e):
                    if curr == h_q:
                        x_part = k
                        break
                    curr = (curr * g_q) % p

            remainders.append(x_part)
            moduli.append(q_e)

        # Chinesischer Restsatz: kombiniere alle x mod q^e
        x = _crt(remainders, moduli)
        if x is None:
            return None

        # Verifikation
        if pow(g, x, p) == h % p:
            return x

        return None

    @staticmethod
    def index_calculus_demo(g: int, h: int, p: int, B: int = 20) -> dict:
        """
        Indexkalkül-Algorithmus (vereinfachte Demo) für diskreten Logarithmus.

        Der Indexkalkül ist der beste bekannte Algorithmus für diskrete Logarithmen
        in Zp* (subexponentiell: L[1/2]).

        Idee:
        1. Faktorbasis: kleine Primzahlen P = {p₁,...,pₖ} mit pᵢ ≤ B
        2. Beziehungsphase: finde s mit g^s B-glatt, schreibe g^s = ∏ pᵢ^aᵢ
        3. Lineare Algebra mod (p-1): löse log-Gleichungen
        4. Individuelle Logarithmen: finde r mit h·g^r B-glatt

        Diese Demo-Version zeigt die Struktur des Algorithmus für kleine p.

        @param g: Generator/Basis
        @param h: Zielwert
        @param p: Primzahlmodulus (klein für Demo)
        @param B: Größe der Faktorbasis
        @return: Dict mit 'result', 'factor_base', 'relations', 'log_found'
        @lastModified: 2026-03-10
        """
        # Faktorbasis: Primzahlen ≤ B
        factor_base = _sieve_eratosthenes(B)

        def try_smooth(val: int, fb: list) -> Optional[list]:
            """Prüft ob val über Faktorbasis glatt ist, gibt Exponentvektor zurück."""
            exponents = [0] * len(fb)
            v = val % p
            if v == 0:
                return None
            for i, prime in enumerate(fb):
                while v % prime == 0:
                    exponents[i] += 1
                    v //= prime
            return exponents if v == 1 else None

        # Beziehungsphase: sammle glatte Potenzen von g
        relations = []
        for s in range(1, min(p - 1, 500)):
            val = pow(g, s, p)
            exp = try_smooth(val, factor_base)
            if exp is not None:
                relations.append({'s': s, 'value': val, 'exponents': exp})
            if len(relations) >= len(factor_base) + 5:
                break

        # Fallback: nutze BSGS für das eigentliche Ergebnis
        result = DiscreteLogarithm.baby_step_giant_step(g, h, p)

        return {
            'result': result,
            'factor_base': factor_base,
            'relations': relations[:5],  # Zeige nur ersten 5 Beziehungen
            'log_found': result is not None,
            'algorithm': 'Indexkalkül-Demo (vereinfacht)'
        }

    @staticmethod
    def discrete_log_brute_force(g: int, h: int, p: int) -> Optional[int]:
        """
        Brute-Force-Suche für diskreten Logarithmus (nur für kleine p).

        Testet alle x von 0 bis p-2 bis g^x ≡ h (mod p) gilt.
        Zeitkomplexität: O(p) - nur für sehr kleine p (< 10^6) praktikabel.

        @param g: Basis/Generator
        @param h: Zielwert
        @param p: Primzahlmodulus (muss klein sein, < 10^6 empfohlen)
        @return: x mit g^x ≡ h (mod p), oder None
        @lastModified: 2026-03-10
        """
        if p > 10**6:
            raise ValueError(f"p zu groß für Brute Force: {p} > 10^6")

        h_mod = h % p
        current = 1  # g^0 mod p = 1
        for x in range(p):
            if current == h_mod:
                return x
            current = (current * g) % p
        return None


def _crt(remainders: list[int], moduli: list[int]) -> Optional[int]:
    """
    Chinesischer Restsatz (CRT) - kombiniert kongruente Gleichungen.

    Löst: x ≡ rᵢ (mod mᵢ) für i = 1,...,k (paarweise teilerfremde Moduli).

    @param remainders: Liste der Reste [r₁, r₂, ..., rₖ]
    @param moduli: Liste der Moduli [m₁, m₂, ..., mₖ] (müssen teilerfremd sein)
    @return: Eindeutige Lösung x mod M (M = m₁·m₂·...·mₖ), oder None
    @lastModified: 2026-03-10
    """
    M = 1
    for m in moduli:
        M *= m

    x = 0
    for r, m in zip(remainders, moduli):
        Mi = M // m
        inv = _mod_inverse(Mi, m)
        if inv is None:
            return None
        x += r * Mi * inv

    return x % M


# ===========================================================================
# 5. KETTENBRÜCHE UND FAKTORISIERUNG
# ===========================================================================

class ContinuedFractionFactoring:
    """
    Kettenbruchentwicklung und kettenbruchbasierte Faktorisierung.

    Kettenbrüche sind eine wichtige Verbindung zwischen Analysis und Zahlentheorie.
    Die Kettenbruchentwicklung von √n liefert die beste rationale Approximation.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def continued_fraction_expansion(n: int, max_terms: int = 50) -> list[int]:
        """
        Kettenbruchentwicklung von √n.

        Jede irrationale Quadratwurzel √n (n nicht-quadratische Zahl) hat eine
        periodische Kettenbruchentwicklung:
            √n = [a₀; a₁, a₂, ..., aₖ, a₁, a₂, ..., aₖ, ...]

        Algorithmus:
            m₀ = 0, d₀ = 1, a₀ = ⌊√n⌋
            mᵢ₊₁ = dᵢ·aᵢ - mᵢ
            dᵢ₊₁ = (n - mᵢ₊₁²) / dᵢ
            aᵢ₊₁ = ⌊(a₀ + mᵢ₊₁) / dᵢ₊₁⌋

        @param n: Positive ganze Zahl (kein perfektes Quadrat für Irrationalität)
        @param max_terms: Maximale Anzahl Kettenbruchterme
        @return: Liste der Kettenbruchkoeffizienten [a₀; a₁, a₂, ...]
        @lastModified: 2026-03-10
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        a0 = math.isqrt(n)
        # Prüfe ob n ein perfektes Quadrat ist
        if a0 * a0 == n:
            return [a0]  # Rationale Zahl → endlicher Kettenbruch

        coefficients = [a0]
        m, d, a = 0, 1, a0

        for _ in range(max_terms):
            # Nächster Schritt der Kettenbruchentwicklung
            m = d * a - m
            d = (n - m * m) // d
            if d == 0:
                break
            a = (a0 + m) // d
            coefficients.append(a)

            # Periode erkannt wenn a = 2·a₀
            if a == 2 * a0:
                break

        return coefficients

    @staticmethod
    def convergents(cf_coefficients: list[int]) -> list[tuple[int, int]]:
        """
        Berechnet die Konvergenten pₖ/qₖ einer Kettenbruchentwicklung.

        Die k-te Konvergente ist die beste rationale Approximation mit Nenner ≤ qₖ.

        Rekursionsformel:
            p₋₁ = 1, p₀ = a₀
            q₋₁ = 0, q₀ = 1
            pₖ = aₖ · pₖ₋₁ + pₖ₋₂
            qₖ = aₖ · qₖ₋₁ + qₖ₋₂

        @param cf_coefficients: Liste der Kettenbruchkoeffizienten [a₀, a₁, ...]
        @return: Liste von Konvergenten als (Zähler, Nenner)-Paare
        @lastModified: 2026-03-10
        """
        if not cf_coefficients:
            return []

        convs = []
        # Initialisierung: p₋₁=1, p₀=a₀; q₋₁=0, q₀=1
        p_prev, p_curr = 1, cf_coefficients[0]
        q_prev, q_curr = 0, 1

        convs.append((p_curr, q_curr))

        for a in cf_coefficients[1:]:
            # Rekursion: pₖ = aₖ·pₖ₋₁ + pₖ₋₂
            p_next = a * p_curr + p_prev
            q_next = a * q_curr + q_prev
            convs.append((p_next, q_next))
            p_prev, p_curr = p_curr, p_next
            q_prev, q_curr = q_curr, q_next

        return convs

    @staticmethod
    def pell_equation_solve(D: int) -> tuple[int, int]:
        """
        Löst die Pell'sche Gleichung x² - D·y² = 1 via Kettenbruchentwicklung.

        Die fundamentale Lösung (x₁, y₁) ist die kleinste positive Lösung
        und erzeugt alle weiteren Lösungen via:
            xₙ + yₙ·√D = (x₁ + y₁·√D)ⁿ

        Algorithmus: Die fundamentale Lösung findet man unter den Konvergenten
        der Kettenbruchentwicklung von √D. Falls die Periode k gerade ist,
        nimm die (k-1)-te Konvergente. Falls k ungerade, die (2k-1)-te.

        @param D: Positive ganze Zahl (kein perfektes Quadrat)
        @return: (x₁, y₁) - fundamentale Lösung der Pell-Gleichung x²-Dy²=1
        @raises ValueError: Wenn D ein perfektes Quadrat ist
        @lastModified: 2026-03-10
        """
        # Perfekte Quadrate haben keine Pell-Gleichung
        if _is_perfect_square(D)[0]:
            raise ValueError(f"D = {D} ist ein perfektes Quadrat - keine Pell-Gleichung")

        # Kettenbruchentwicklung von √D berechnen
        cf = ContinuedFractionFactoring.continued_fraction_expansion(D, max_terms=200)
        period_len = len(cf) - 1  # Periodenlänge (ohne a₀)

        if period_len == 0:
            raise ValueError(f"Keine Periode gefunden für D = {D}")

        # Bestimme welche Konvergente die Lösung liefert
        # Bei gerader Periode: (periode-1)-te Konvergente
        # Bei ungerader Periode: (2·periode-1)-te Konvergente
        if period_len % 2 == 1:
            # Ungerade Periode: brauche volle zwei Perioden
            target_idx = 2 * period_len - 1
            # Kettenbruch mit zwei Perioden generieren
            cf_extended = [cf[0]]
            period = cf[1:]
            for _ in range(2 * period_len):
                cf_extended.append(period[(len(cf_extended) - 1) % period_len])
        else:
            # Gerade Periode: eine Periode reicht
            target_idx = period_len - 1
            cf_extended = cf

        # Konvergenten berechnen bis zum Zielindex
        convs = ContinuedFractionFactoring.convergents(cf_extended[:target_idx + 1])

        if not convs:
            raise ValueError(f"Keine Konvergente gefunden für D = {D}")

        x, y = convs[-1]

        # Verifikation
        if x * x - D * y * y == 1:
            return x, y

        # Fallback: brute-force Suche für kleine D
        y_val = 1
        while y_val < 10**6:
            val = D * y_val * y_val + 1
            is_sq, x_val = _is_perfect_square(val)
            if is_sq:
                return x_val, y_val
            y_val += 1

        raise ValueError(f"Keine fundamentale Lösung für D = {D} gefunden")

    @staticmethod
    def cfrac_factoring_demo(n: int) -> dict:
        """
        CFRAC (Continued Fraction Factoring) - Faktorisierungsmethode via Kettenbrüche.

        CFRAC (Morrison & Brillhart, 1975) war der erste Faktorisierungsalgorithmus
        mit subexponentieller Laufzeit (vor dem Quadratischen Sieb).

        Idee: Suche x, y mit x² ≡ y² (mod n) und x ≢ ±y (mod n).
        Dann ist ggT(x-y, n) ein nicht-trivialer Teiler von n.

        Die Kettenbruchentwicklung von √n liefert Konvergenten pₖ/qₖ mit
        pₖ² ≡ (-1)^(k+1) · Rₖ (mod n), wobei Rₖ oft klein und damit glatt ist.

        Diese Demo zeigt das Prinzip für kleine n.

        @param n: Zu faktorisierende Zahl (n ≥ 4, ungerade, zusammengesetzt)
        @return: Dict mit 'factor', 'method', 'convergents_checked', 'success'
        @lastModified: 2026-03-10
        """
        if n < 4 or n % 2 == 0:
            return {'factor': None, 'method': 'CFRAC-Demo', 'convergents_checked': 0, 'success': False}

        if miller_rabin_primality_test(n):
            return {'factor': None, 'method': 'CFRAC-Demo', 'convergents_checked': 0, 'success': False}

        # Kettenbruchentwicklung von √n mit mehr Termen
        cf = ContinuedFractionFactoring.continued_fraction_expansion(n, max_terms=200)

        # Periodische Erweiterung für mehr Konvergenten
        period = cf[1:]
        cf_long = [cf[0]]
        if period:
            for i in range(200):
                cf_long.append(period[i % len(period)])

        # Konvergenten berechnen
        convs = ContinuedFractionFactoring.convergents(cf_long)
        checked = 0

        for p_k, q_k in convs:
            checked += 1
            # Prüfe ob pₖ² - n·qₖ² ein brauchbares Residuum liefert
            residuum = (p_k * p_k) % n
            # Suche Paare (i, j) mit residuum_i · residuum_j = Quadratzahl
            is_sq, sqrt_res = _is_perfect_square(residuum)

            if is_sq and sqrt_res != 0:
                # x = pₖ mod n, y = sqrt(pₖ² mod n) mod n
                x = p_k % n
                y = sqrt_res % n
                if x != y and x != n - y:
                    d = _gcd(abs(x - y), n)
                    if 1 < d < n:
                        return {
                            'factor': d,
                            'cofactor': n // d,
                            'method': 'CFRAC-Demo',
                            'convergents_checked': checked,
                            'success': True
                        }

        # Fallback auf Pollard's ρ
        d = IntegerFactorization.pollard_rho(n)
        return {
            'factor': d,
            'cofactor': n // d if d else None,
            'method': 'CFRAC-Demo (Fallback: Pollard-ρ)',
            'convergents_checked': checked,
            'success': d is not None
        }


# ===========================================================================
# 6. KRYPTOGRAPHISCHE ALGORITHMEN (MATHEMATISCHE SEITE)
# ===========================================================================

class CryptographicAlgorithms:
    """
    Mathematische Grundlagen kryptographischer Protokolle.

    Implementiert die zahlentheoretischen Kernalgorithmen hinter
    Diffie-Hellman, ElGamal und anderen Kryptosystemen.

    HINWEIS: Diese Implementierung dient ausschließlich Lernzwecken.
    Für produktive Kryptographie sind etablierte Bibliotheken zu verwenden.

    @author: Michael Fuhrmann
    @since: 2026-03-10
    @lastModified: 2026-03-10
    """

    @staticmethod
    def diffie_hellman_demo(p: int = 23, g: int = 5) -> dict:
        """
        Demonstration des Diffie-Hellman-Schlüsselaustausches.

        DH ermöglicht zwei Parteien (Alice, Bob) einen geheimen Schlüssel
        über einen unsicheren Kanal auszutauschen, ohne diesen direkt zu übertragen.

        Protokoll:
            Alice wählt privates a → sendet A = g^a mod p
            Bob wählt privates b   → sendet B = g^b mod p
            Gemeinsamer Schlüssel: s = A^b = B^a = g^(ab) mod p

        Sicherheit basiert auf der Schwierigkeit des diskreten Logarithmus.

        @param p: Primzahlmodulus (Standard: 23)
        @param g: Generator/Primitivwurzel modulo p (Standard: 5)
        @return: Dict mit allen DH-Parametern und dem gemeinsamen Schlüssel
        @lastModified: 2026-03-10
        """
        # Private Schlüssel (in der Praxis: zufällig und groß)
        a = random.randint(2, p - 2)  # Alice's privater Schlüssel
        b = random.randint(2, p - 2)  # Bob's privater Schlüssel

        # Öffentliche Schlüssel berechnen
        A = pow(g, a, p)  # Alice's öffentlicher Schlüssel
        B = pow(g, b, p)  # Bob's öffentlicher Schlüssel

        # Gemeinsamer Schlüssel (beide Seiten berechnen dasselbe)
        shared_alice = pow(B, a, p)  # Alice: B^a mod p = g^(ab) mod p
        shared_bob = pow(A, b, p)    # Bob: A^b mod p = g^(ab) mod p

        assert shared_alice == shared_bob, "DH-Fehler: Schlüssel stimmen nicht überein!"

        return {
            'p': p,
            'g': g,
            'alice_private': a,
            'bob_private': b,
            'alice_public': A,
            'bob_public': B,
            'shared_key': shared_alice,
            'protocol': 'Diffie-Hellman'
        }

    @staticmethod
    def elgamal_keygen(p: int = 23, g: int = 5) -> dict:
        """
        ElGamal-Schlüsselgenerierung.

        ElGamal ist ein asymmetrisches Kryptosystem basierend auf dem
        diskreten Logarithmus. Entwickelt von Taher ElGamal (1985).

        Schlüsselerzeugung:
            Wähle private Zahl x mit 1 < x < p-1
            Berechne öffentlichen Schlüssel y = g^x mod p

        @param p: Primzahlmodulus
        @param g: Generator (primitiver Root modulo p)
        @return: Dict mit 'public_key' (p, g, y) und 'private_key' x
        @lastModified: 2026-03-10
        """
        # Privater Schlüssel: zufälliges x
        x = random.randint(2, p - 2)
        # Öffentlicher Schlüssel: y = g^x mod p
        y = pow(g, x, p)

        return {
            'public_key': {'p': p, 'g': g, 'y': y},
            'private_key': x,
            'algorithm': 'ElGamal'
        }

    @staticmethod
    def elgamal_encrypt(message: int, pub_key: dict, p: int, g: int) -> tuple[int, int]:
        """
        ElGamal-Verschlüsselung.

        Verschlüsselung einer Nachricht m:
            Wähle zufälliges k mit 1 < k < p-1
            c₁ = g^k mod p       (ephemerer Teil)
            c₂ = m · y^k mod p   (verschlüsselte Nachricht)

        @param message: Nachricht als Zahl (0 ≤ m < p)
        @param pub_key: Öffentlicher Schlüssel y
        @param p: Primzahlmodulus
        @param g: Generator
        @return: Ciphertext (c₁, c₂)
        @lastModified: 2026-03-10
        """
        y = pub_key  # Öffentlicher Schlüssel y
        # Zufälliges ephemeres k
        k = random.randint(2, p - 2)
        # Verschlüsselung
        c1 = pow(g, k, p)           # g^k mod p
        c2 = (message * pow(y, k, p)) % p   # m · y^k mod p
        return c1, c2

    @staticmethod
    def elgamal_decrypt(ciphertext: tuple[int, int], priv_key: int, p: int) -> int:
        """
        ElGamal-Entschlüsselung.

        Entschlüsselung von (c₁, c₂):
            s = c₁^x mod p                     (geteiltes Geheimnis)
            m = c₂ · s⁻¹ mod p                 (Nachricht)

        Korrektheit: c₂ · s⁻¹ = m · y^k · (g^(kx))⁻¹ = m · g^(kx) · g^(-kx) = m

        @param ciphertext: Tupel (c₁, c₂)
        @param priv_key: Privater Schlüssel x
        @param p: Primzahlmodulus
        @return: Entschlüsselte Nachricht als ganze Zahl
        @lastModified: 2026-03-10
        """
        c1, c2 = ciphertext
        x = priv_key
        # Geteiltes Geheimnis: s = c₁^x mod p
        s = pow(c1, x, p)
        # Inverses von s
        s_inv = _mod_inverse(s, p)
        if s_inv is None:
            raise ValueError("Entschlüsselung fehlgeschlagen: kein Inverses")
        # Nachricht: m = c₂ · s⁻¹ mod p
        return (c2 * s_inv) % p

    @staticmethod
    def rsa_attack_small_e(e: int, n: int, c: int) -> Optional[int]:
        """
        RSA-Angriff bei kleinem Exponent e (Low Public Exponent Attack).

        Falls e klein ist (z.B. e=3) und die Nachricht m klein ist,
        sodass m^e < n, dann gilt: c = m^e (ohne mod-Reduktion).
        Man kann direkt die e-te Wurzel berechnen: m = c^(1/e).

        Dieser Angriff zeigt die Schwäche von zu kleinen RSA-Exponenten
        ohne Padding (wie z.B. OAEP).

        @param e: Öffentlicher Exponent (z.B. 3)
        @param n: RSA-Modulus
        @param c: Ciphertext
        @return: Klartext m falls Angriff erfolgreich, sonst None
        @lastModified: 2026-03-10
        """
        # Versuche e-te Ganzzahlwurzel von c
        # Binärsuche für die e-te Wurzel
        low, high = 0, min(c, n)

        # e-te Wurzel von c berechnen
        m = round(c ** (1.0 / e))

        # Überprüfe Kandidaten in der Umgebung
        for candidate in range(max(0, m - 2), m + 3):
            if candidate ** e == c:
                return candidate  # Angriff erfolgreich

        # Falls c < n^e, direkte Kubikwurzel möglich
        # Binärsuche für große Zahlen
        lo, hi = 0, c
        while lo <= hi:
            mid = (lo + hi) // 2
            val = mid ** e
            if val == c:
                return mid
            elif val < c:
                lo = mid + 1
            else:
                hi = mid - 1

        return None  # Angriff fehlgeschlagen (m^e wurde mod n reduziert)

    @staticmethod
    def baby_step_giant_step_ecc(
        P: tuple[int, int],
        Q: tuple[int, int],
        n: int,
        a: int,
        p: int
    ) -> Optional[int]:
        """
        Baby-Step-Giant-Step auf elliptischer Kurve (vereinfachte Demo).

        Löst das diskrete Logarithmusproblem auf einer elliptischen Kurve:
        Finde k mit k·P = Q auf der Kurve y² ≡ x³ + ax (mod p).

        Idee (analog zu BSGS in Zp*):
        - Baby Steps: Berechne j·P für j = 0,...,m-1
        - Giant Steps: Berechne Q - i·(m·P) für i = 0,...,m
        - Kollision: j·P = Q - i·(m·P) → k = im + j

        @param P: Basispunkt (x, y) auf der Kurve
        @param Q: Zielpunkt (x, y) auf der Kurve
        @param n: Ordnung der Kurvengruppe (Anzahl Punkte)
        @param a: Kurvenparameter (y² = x³ + ax)
        @param p: Primzahlmodulus
        @return: k mit k·P = Q, oder None
        @lastModified: 2026-03-10
        """

        def point_add(P1: Optional[tuple], P2: Optional[tuple]) -> Optional[tuple]:
            """Addiert zwei Punkte auf der Kurve y² = x³ + ax (mod p)."""
            if P1 is None:
                return P2  # Neutrales Element (Punkt im Unendlichen)
            if P2 is None:
                return P1
            x1, y1 = P1
            x2, y2 = P2

            if x1 == x2:
                if y1 != y2:
                    return None  # P + (-P) = Neutralelement
                # Verdoppelung: Steigung = (3x²+a) / (2y)
                if y1 == 0:
                    return None
                lam_num = (3 * x1 * x1 + a) % p
                lam_den = _mod_inverse(2 * y1 % p, p)
                if lam_den is None:
                    return None
            else:
                # Addition: Steigung = (y2-y1) / (x2-x1)
                lam_num = (y2 - y1) % p
                lam_den = _mod_inverse((x2 - x1) % p, p)
                if lam_den is None:
                    return None

            lam = (lam_num * lam_den) % p
            x3 = (lam * lam - x1 - x2) % p
            y3 = (lam * (x1 - x3) - y1) % p
            return (x3, y3)

        def point_neg(P1: Optional[tuple]) -> Optional[tuple]:
            """Negiert einen Punkt: (x, y) → (x, -y mod p)."""
            if P1 is None:
                return None
            return (P1[0], (-P1[1]) % p)

        def scalar_mult(k: int, point: Optional[tuple]) -> Optional[tuple]:
            """Skalarmultiplikation k·P via Double-and-Add."""
            result = None
            addend = point
            while k > 0:
                if k % 2 == 1:
                    result = point_add(result, addend)
                addend = point_add(addend, addend)
                k //= 2
            return result

        # Baby Steps: j·P für j = 0,...,m-1
        m = math.isqrt(n) + 1
        baby_steps = {}
        current = None  # 0·P = Neutralelement
        for j in range(m):
            baby_steps[current] = j
            current = point_add(current, P)

        # Giant Steps: Q + i·(-m·P) = Q - i·(m·P)
        mP = scalar_mult(m, P)
        neg_mP = point_neg(mP)
        current = Q
        for i in range(m + 1):
            if current in baby_steps:
                k = (i * m + baby_steps[current]) % n
                # Verifikation
                if scalar_mult(k, P) == Q:
                    return k
            current = point_add(current, neg_mP)

        return None


# ===========================================================================
# 7. STANDALONE-FUNKTIONEN
# ===========================================================================

def euler_product_approximation(s: float, n_primes: int = 100) -> float:
    """
    Approximiert das Euler-Produkt ∏_{p prim} (1 - p^{-s})^{-1} = ζ(s).

    Das Euler-Produkt ist eine der fundamentalen Verbindungen zwischen
    Primzahlen und der Riemann-Zeta-Funktion:
        ζ(s) = ∑_{n=1}^∞ n^{-s} = ∏_{p prim} (1 - p^{-s})^{-1}

    @param s: Reeller Parameter (s > 1 für Konvergenz)
    @param n_primes: Anzahl der Primzahlen für die Approximation
    @return: Approximierter Wert des Euler-Produkts
    @raises ValueError: Wenn s ≤ 1 (Divergenz)
    @lastModified: 2026-03-10
    """
    if s <= 1.0:
        raise ValueError(f"s muss > 1 sein für Konvergenz, erhalten: {s}")

    # Schätze obere Grenze für n_primes Primzahlen (PNT: p_n ≈ n·ln(n))
    limit = max(100, int(n_primes * (math.log(n_primes) + math.log(math.log(n_primes + 2)) + 2)))
    primes = _sieve_eratosthenes(limit)[:n_primes]

    # Euler-Produkt: ∏ (1 - p^{-s})^{-1}
    product = 1.0
    for p in primes:
        product *= 1.0 / (1.0 - p ** (-s))

    return product


def mertens_first_theorem(n: int) -> dict:
    """
    Berechnet Mertens' 1. Satz: ∑_{p≤n} ln(p)/p ≈ ln(n).

    Mertens' erster Satz (1874) besagt:
        ∑_{p≤n, p prim} ln(p) / p = ln(n) + O(1)

    Das heißt die Summe der gewichteten reziproken Primzahlen wächst wie ln(n).

    @param n: Obere Grenze der Summation
    @return: Dict mit 'mertens_sum', 'log_n', 'difference', 'primes_used'
    @lastModified: 2026-03-10
    """
    primes = _sieve_eratosthenes(n)

    # ∑ ln(p)/p für alle Primzahlen p ≤ n
    mertens_sum = sum(math.log(p) / p for p in primes)
    log_n = math.log(n) if n >= 1 else 0.0

    return {
        'mertens_sum': mertens_sum,
        'log_n': log_n,
        'difference': mertens_sum - log_n,
        'primes_used': len(primes),
        'theorem': 'Mertens 1. Satz: Σ ln(p)/p ≈ ln(n)'
    }


def primitive_root(p: int) -> Optional[int]:
    """
    Findet die kleinste primitive Wurzel modulo p (p prim).

    Eine primitive Wurzel g modulo p ist eine Zahl, deren Potenzen
    g⁰, g¹, g², ..., g^{p-2} alle Reste 1,...,p-1 erzeugen.
    (g hat die maximale Ordnung p-1 in der Gruppe ℤ_p*)

    Algorithmus:
        Für jeden Kandidaten g: prüfe ob g^{(p-1)/q} ≢ 1 (mod p)
        für alle Primteiler q von p-1.

    @param p: Primzahl (p ≥ 3)
    @return: Kleinste primitive Wurzel modulo p, oder None wenn p nicht prim
    @lastModified: 2026-03-10
    """
    if p < 2:
        return None
    if p == 2:
        return 1  # Triviale Gruppe
    if not miller_rabin_primality_test(p):
        return None  # p muss prim sein

    # Primfaktorisierung von p-1
    phi = p - 1
    factors = set(_prime_factors_list(phi))

    # Teste alle Kandidaten von 2 bis p-1
    for g in range(2, p):
        # g ist primitiv wenn g^((p-1)/q) ≢ 1 (mod p) für alle Primteiler q
        is_primitive = all(pow(g, phi // q, p) != 1 for q in factors)
        if is_primitive:
            return g

    return None  # Sollte nie eintreten (primitive Wurzel existiert immer für Primes)


def legendre_symbol_jacobi(a: int, n: int) -> int:
    """
    Berechnet das Jacobi-Symbol (a/n) - Verallgemeinerung des Legendre-Symbols.

    Für ungerade Primzahlen p: Jacobi-Symbol = Legendre-Symbol.
    Für zusammengesetzte n: Jacobi-Symbol ist Produkt der Legendre-Symbole
    über alle Primfaktoren (mit Vielfachheit).

    WICHTIG: (a/n) = 1 bedeutet nicht zwingend, dass a ein quadratischer
    Rest mod n ist (nur wenn n prim)!

    Delegiert an jacobi_symbol aus proof_theory.py.

    @param a: Ganzzahl
    @param n: Ungerade positive ganze Zahl > 1
    @return: Jacobi-Symbol ∈ {-1, 0, 1}
    @lastModified: 2026-03-10
    """
    return jacobi_symbol(a, n)
