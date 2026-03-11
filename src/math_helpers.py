"""
@file math_helpers.py
@brief Zentrale Hilfsfunktionen für das gesamte specialist-maths-Projekt.
@description
    Dieses Modul bündelt mathematische Hilfsfunktionen, die in vielen
    anderen Modulen benötigt werden. Es verhindert Codeduplikation und
    stellt eine einheitliche, gut getestete Implementierung bereit.

    Bisher waren diese Funktionen in folgenden Modulen dupliziert:
    - algebra_numbertheory.py (is_prime, prime_factorization, euler_phi)
    - algebraic_number_theory.py (_is_prime, _gcd, _euler_totient, _prime_factorization)
    - galois_representations.py (_is_prime, _primes_up_to)
    - l_functions.py (_is_prime, _gcd, _euler_totient)
    - elliptic_curves.py (_is_prime)
    - und weitere...

    Neue Module sollen von hier importieren statt eigene Kopien zu pflegen.
    Bestehende Module werden schrittweise migriert (rückwärtskompatibel).

    **Refactoring-Stand**: Build 102 (2026-03-11)
    Schritt 1: Dieses Modul erstellt. Neue Module (langlands_program.py etc.)
               importieren direkt von hier.
    Schritt 2 (geplant): Bestehende Module auf math_helpers umstellen.

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
from functools import lru_cache
from typing import Iterator


# =============================================================================
# PRIMZAHLEN
# =============================================================================

@lru_cache(maxsize=10000)
def is_prime(n: int) -> bool:
    """
    @brief Prüft ob n eine Primzahl ist (6k±1-Methode, gecacht).
    @description
        Optimierter Primzahltest:
        1. Sonderfälle n < 2 (keine Primzahl), n=2 und n=3 (Primzahlen)
        2. Gerade Zahlen und Vielfache von 3 ausschließen
        3. Nur Kandidaten der Form 6k±1 bis √n prüfen

        Alle Primzahlen > 3 haben die Form 6k±1, was die Kandidatenmenge
        auf 1/3 der natürlichen Zahlen reduziert.

    @param n Die zu prüfende Zahl (n ≥ 0).
    @return True wenn n prim ist, sonst False.
    @lastModified 2026-03-11
    """
    # Sonderfälle
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    # Gerade Zahlen und Vielfache von 3 sofort ausschließen
    if n % 2 == 0 or n % 3 == 0:
        return False
    # Nur 6k±1-Kandidaten bis √n testen
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


@lru_cache(maxsize=1000)
def prime_factorization(n: int) -> dict:
    """
    @brief Berechnet die Primfaktorzerlegung n = p1^e1 * p2^e2 * ... (gecacht).
    @description
        Zerlegt n durch Probedivision:
        1. Teile durch 2 solange möglich
        2. Teile durch ungerade Zahlen 3, 5, 7, ... bis √n
        3. Falls Restzahl > 1: sie selbst ist ein Primfaktor

    @param n Ganze Zahl ≥ 2.
    @return dict {primzahl: exponent}, z.B. {2:3, 3:2, 5:1} für n=360.
    @lastModified 2026-03-11
    """
    if n < 2:
        return {}
    factors: dict[int, int] = {}
    # Faktor 2 separat herausziehen
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    # Ungerade Faktoren ab 3
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2
    # Verbleibender Primfaktor
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def primes_up_to(limit: int) -> list[int]:
    """
    @brief Gibt alle Primzahlen ≤ limit via Sieb des Eratosthenes zurück.
    @description
        Klassisches Siebverfahren in O(n·log·log·n) Zeit und O(n) Speicher.
        Für kleine limit (<10000) oft schneller als Prime-Test-Schleife.

    @param limit Obere Grenze (inklusive).
    @return Sortierte Liste aller Primzahlen ≤ limit.
    @lastModified 2026-03-11
    """
    if limit < 2:
        return []
    # Sieb initialisieren: True = noch nicht ausgesiebt
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    # Vielfache jeder Primzahl aussieben
    i = 2
    while i * i <= limit:
        if sieve[i]:
            # Alle Vielfachen ab i² streichen
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [i for i in range(2, limit + 1) if sieve[i]]


def primes_generator(start: int = 2) -> Iterator[int]:
    """
    @brief Generator für Primzahlen ab start (unbegrenzt, lazy).
    @description
        Liefert Primzahlen der Reihe nach via is_prime()-Test.
        Nützlich wenn die Grenze unbekannt ist oder nur wenige Primzahlen
        benötigt werden.

    @param start Startwert (Standard: 2).
    @return Iterator über Primzahlen ≥ start.
    @lastModified 2026-03-11
    """
    n = max(2, start)
    while True:
        if is_prime(n):
            yield n
        n += 1


def first_n_primes(n: int) -> list[int]:
    """
    @brief Gibt die ersten n Primzahlen zurück.
    @param n Anzahl der gewünschten Primzahlen.
    @return Liste der ersten n Primzahlen.
    @lastModified 2026-03-11
    """
    result: list[int] = []
    gen = primes_generator()
    for _ in range(n):
        result.append(next(gen))
    return result


# =============================================================================
# GRUNDLEGENDE ZAHLENTHEORIE
# =============================================================================

def gcd(a: int, b: int) -> int:
    """
    @brief Größter gemeinsamer Teiler via euklidischem Algorithmus.
    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return ggT(a, b) ≥ 0.
    @lastModified 2026-03-11
    """
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def lcm(a: int, b: int) -> int:
    """
    @brief Kleinstes gemeinsames Vielfaches.
    @param a Erste ganze Zahl (≠ 0).
    @param b Zweite ganze Zahl (≠ 0).
    @return kgV(a, b) > 0.
    @lastModified 2026-03-11
    """
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """
    @brief Erweiterter euklidischer Algorithmus: ggT und Bézout-Koeffizienten.
    @description
        Berechnet (g, x, y) mit g = ggT(a, b) und a·x + b·y = g.

    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return Tupel (g, x, y): Bezout-Koeffizienten a·x + b·y = g.
    @lastModified 2026-03-11
    """
    if b == 0:
        return abs(a), (1 if a >= 0 else -1), 0
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def mod_inverse(a: int, m: int) -> int:
    """
    @brief Modulares Inverses: a^{-1} mod m via erweitertem euklidischem Alg.
    @description
        Existiert nur wenn ggT(a, m) = 1.

    @param a Zahl, deren Inverses gesucht wird.
    @param m Modulus (m > 1).
    @return x mit a·x ≡ 1 (mod m), 0 ≤ x < m.
    @raises ValueError wenn ggT(a, m) ≠ 1.
    @lastModified 2026-03-11
    """
    g, x, _ = extended_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"Kein modulares Inverses: ggT({a}, {m}) = {g} ≠ 1")
    return x % m


@lru_cache(maxsize=500)
def euler_phi(n: int) -> int:
    """
    @brief Eulersche Phi-Funktion φ(n) — Anzahl der zu n teilerfremden Zahlen 1..n.
    @description
        Berechnung via Primfaktorzerlegung:
            φ(n) = n · Π_{p|n} (1 - 1/p)

        Eigenschaften:
        - φ(1) = 1
        - φ(p) = p−1 für Primzahlen p
        - φ(p^k) = p^{k-1}·(p-1)
        - φ ist multiplikativ: φ(mn) = φ(m)φ(n) wenn ggT(m,n)=1

    @param n Positive ganze Zahl.
    @return φ(n).
    @lastModified 2026-03-11
    """
    if n <= 0:
        raise ValueError(f"phi({n}): n muss positiv sein")
    if n == 1:
        return 1
    # Phi via Primfaktorzerlegung berechnen
    result = n
    for p in prime_factorization(n):
        result -= result // p
    return result


def legendre_symbol(a: int, p: int) -> int:
    """
    @brief Legendre-Symbol (a/p) für Primzahl p.
    @description
        (a/p) = 0  wenn p|a
        (a/p) = 1  wenn a quadratischer Rest mod p
        (a/p) = -1 wenn a quadratischer Nichtrest mod p

        Berechnung via Euler-Kriterium: (a/p) ≡ a^{(p-1)/2} (mod p).

    @param a Ganze Zahl.
    @param p Ungerade Primzahl.
    @return 0, 1 oder -1.
    @raises ValueError wenn p gerade oder < 2.
    @lastModified 2026-03-11
    """
    if p < 2 or p % 2 == 0:
        raise ValueError(f"Legendre-Symbol: p={p} muss ungerade Primzahl sein")
    a = a % p
    if a == 0:
        return 0
    val = pow(a, (p - 1) // 2, p)
    return 1 if val == 1 else -1


def jacobi_symbol(a: int, n: int) -> int:
    """
    @brief Jacobi-Symbol (a/n) — Verallgemeinerung des Legendre-Symbols.
    @description
        Für ungerades n > 0 definiert als Produkt der Legendre-Symbole
        über die Primfaktoren von n.

        Achtung: (a/n)=1 bedeutet NICHT zwingend, dass a ein QR mod n ist!

    @param a Ganze Zahl.
    @param n Positive ungerade ganze Zahl.
    @return 0, 1 oder -1.
    @lastModified 2026-03-11
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError(f"Jacobi-Symbol: n={n} muss positive ungerade Zahl sein")
    a = a % n
    result = 1
    while a != 0:
        # Faktor 2 herausziehen
        while a % 2 == 0:
            a //= 2
            # Quadratisches Reziprozitätsgesetz für den Faktor 2
            if n % 8 in (3, 5):
                result = -result
        # a und n tauschen (Reziprozitätsgesetz)
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a = a % n
    return result if n == 1 else 0


def chinese_remainder_theorem(remainders: list[int], moduli: list[int]) -> tuple[int, int]:
    """
    @brief Chinesischer Restsatz: löst x ≡ r_i (mod m_i) für paarweise teilerfremde m_i.
    @description
        Konstruktive Lösung des simultanen Kongruenzsystems via CRT:
            x = Σ_i (r_i · M_i · M_i^{-1}) mod M
        wobei M = Π m_i und M_i = M / m_i.

    @param remainders Liste der Reste [r_0, r_1, ..., r_{k-1}].
    @param moduli     Liste der paarweise teilerfremden Moduln [m_0, ..., m_{k-1}].
    @return Tupel (x, M): kleinste nichtnegative Lösung x und Gesamtmodulus M.
    @raises ValueError wenn Moduln nicht paarweise teilerfremd sind.
    @lastModified 2026-03-11
    """
    if len(remainders) != len(moduli):
        raise ValueError("Anzahl der Reste und Moduln muss übereinstimmen")
    # Gesamtmodulus M = Produkt aller Moduln
    M = 1
    for m in moduli:
        M *= m
    # CRT-Konstruktion
    x = 0
    for r, m in zip(remainders, moduli):
        # M_i = M / m_i
        Mi = M // m
        # M_i^{-1} mod m_i
        inv = mod_inverse(Mi, m)
        x += r * Mi * inv
    return x % M, M


# =============================================================================
# EXPORT-LISTE (für from math_helpers import *)
# =============================================================================

__all__ = [
    # Primzahlen
    "is_prime",
    "prime_factorization",
    "primes_up_to",
    "primes_generator",
    "first_n_primes",
    # Grundlegende Zahlentheorie
    "gcd",
    "lcm",
    "extended_gcd",
    "mod_inverse",
    "euler_phi",
    "legendre_symbol",
    "jacobi_symbol",
    "chinese_remainder_theorem",
]
