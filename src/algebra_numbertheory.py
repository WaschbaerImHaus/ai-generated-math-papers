"""
@file algebra_numbertheory.py
@brief Zahlentheorie: Primzahlen, Euler-Phi, RSA-Kryptosystem.
@description
    Enthält zahlentheoretische Funktionen und das RSA-Kryptosystem:

    - is_prime()            – Primzahltest via Probedivision (gecacht)
    - prime_factorization() – Primfaktorzerlegung n = p1^e1 * p2^e2 * ...
    - euler_phi()           – Eulersche Phi-Funktion φ(n)
    - rsa_keygen()          – RSA-Schlüsselerzeugung aus zwei Primzahlen
    - rsa_encrypt()         – RSA-Verschlüsselung: c = m^e mod n
    - rsa_decrypt()         – RSA-Entschlüsselung: m = c^d mod n

    Importiert extended_gcd und mod_inverse aus algebra_core.
    Ausgelagert aus algebra.py für bessere Modularität.

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import math
import functools
from typing import Callable

# Abhängigkeiten aus dem Kern-Modul importieren
from algebra_core import extended_gcd, mod_inverse


# =============================================================================
# PRIMZAHLTEST UND PRIMFAKTORZERLEGUNG
# =============================================================================

@functools.lru_cache(maxsize=10000)
def is_prime(n: int) -> bool:
    """
    @brief Prüft, ob eine natürliche Zahl eine Primzahl ist.
    @description
        Ein Primzahltest durch Probedivision mit Optimierungen:
        1. Sonderfälle: n < 2 ist keine Primzahl; 2 und 3 sind Primzahlen
        2. Zahlen der Form 6k ± 1 überprüfen (alle Primzahlen > 3 haben diese Form)
        3. Nur bis sqrt(n) testen (wenn n einen Teiler > sqrt(n) hat,
           hat es auch einen < sqrt(n))

        Laufzeit: O(sqrt(n)) – Ergebnis wird via lru_cache gecacht (bis 10000 Einträge).
        Deutlicher Geschwindigkeitsvorteil bei wiederholten Aufrufen (z.B. in proof_theory).

    @param n Die zu prüfende natürliche Zahl.
    @return True wenn n eine Primzahl ist, sonst False.
    @author Michael Fuhrmann
    @lastModified 2026-03-10

    Beispiele:
    >>> is_prime(7)
    True
    >>> is_prime(4)
    False
    >>> is_prime(2)
    True
    >>> is_prime(1)
    False
    """
    if n < 2:
        return False
    if n < 4:
        return True   # 2 und 3 sind Primzahlen
    if n % 2 == 0 or n % 3 == 0:
        return False  # Geradzahlige und Vielfache von 3

    # Alle möglichen Teiler der Form 6k ± 1 bis sqrt(n) testen
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6

    return True


@functools.lru_cache(maxsize=1000)
def prime_factorization(n: int) -> dict[int, int]:
    """
    @brief Berechnet die Primfaktorzerlegung einer natürlichen Zahl.
    @description
        Zerlegt n in seine Primfaktoren:
            n = p1^e1 * p2^e2 * ... * pk^ek

        Algorithmus:
        1. Teile durch 2 solange möglich
        2. Teile durch ungerade Zahlen bis sqrt(n)
        3. Falls Rest > 1: Rest selbst ist ein Primfaktor

        Beispiel: 360 = 2^3 * 3^2 * 5^1

        Laufzeit: O(sqrt(n)) – Ergebnis wird via lru_cache gecacht (bis 1000 Einträge).
        Vorteil bei Euler-Phi und anderen zahlentheoretischen Funktionen, die dieselbe
        Primfaktorzerlegung mehrfach benötigen.

    @param n Die zu zerlegende natürliche Zahl (n >= 2).
    @return Dictionary {primzahl: exponent}.
    @author Michael Fuhrmann
    @lastModified 2026-03-10

    Beispiele:
    >>> prime_factorization(12)
    {2: 2, 3: 1}
    >>> prime_factorization(360)
    {2: 3, 3: 2, 5: 1}
    >>> prime_factorization(7)
    {7: 1}
    """
    factors = {}

    # Faktor 2 separat behandeln (häufigster kleiner Primfaktor)
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2

    # Ungerade Faktoren ab 3 testen
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors[i] = factors.get(i, 0) + 1
            n //= i
        i += 2

    # Restlicher Faktor > 1 ist selbst eine Primzahl
    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


@functools.lru_cache(maxsize=1000)
def euler_phi(n: int) -> int:
    """
    @brief Berechnet Eulers Phi-Funktion phi(n).
    @description
        phi(n) zählt die Anzahl der ganzen Zahlen von 1 bis n,
        die zu n teilerfremd sind.

        Berechnung über die Primfaktorzerlegung:
            phi(n) = n * Produkt(1 - 1/p) für alle Primfaktoren p von n

        Eigenschaften:
        - phi(1) = 1
        - phi(p) = p-1 für Primzahlen p
        - phi(p^k) = p^(k-1) * (p-1)
        - phi ist multiplikativ: phi(m*n) = phi(m)*phi(n) wenn ggT(m,n)=1

        Anwendung: Eulers Satz (a^phi(n) ≡ 1 mod n), RSA

        Laufzeit: O(sqrt(n)) – Ergebnis wird via lru_cache gecacht (bis 1000 Einträge).
        Wichtig bei der Farey-Folge-Länge (Σ φ(k)) und RSA-Berechnungen.

    @param n Positive ganze Zahl.
    @return phi(n) - Anzahl der teilerfremden Zahlen.
    @author Michael Fuhrmann
    @lastModified 2026-03-10

    Beispiele:
    >>> euler_phi(12)
    4
    >>> euler_phi(7)
    6
    >>> euler_phi(1)
    1
    """
    if n == 1:
        return 1

    # Primfaktorzerlegung über prime_factorization() wiederverwenden (gecacht!)
    # Spart doppelte Traversierung und profitiert vom lru_cache
    factors = prime_factorization(n)

    # Euler-Phi-Formel: phi(n) = n * prod(1 - 1/p) für alle distinkten Primfaktoren p
    # Äquivalent: result = result * (p-1) // p für jeden Primfaktor p
    result = n
    for p in factors:
        # Multiplikation mit (p-1)/p: erst durch p dividieren (ganzzahlig), dann mal (p-1)
        # Reihenfolge wichtig: result ist immer durch p teilbar wegen Primfaktorzerlegung
        result = result // p * (p - 1)

    return result


# =============================================================================
# RSA-KRYPTOSYSTEM
# =============================================================================

def rsa_keygen(p: int, q: int) -> tuple[tuple[int, int], tuple[int, int]]:
    """
    @brief Erzeugt RSA-Schlüsselpaar aus zwei Primzahlen.
    @description
        RSA (Rivest–Shamir–Adleman, 1977) ist das bekannteste Public-Key-Kryptosystem.
        Sicherheit beruht auf der Schwierigkeit der Primfaktorzerlegung großer Zahlen.

        Schlüsselerzeugung:
            1. n = p·q (Modulus)
            2. λ(n) = lcm(p-1, q-1)  (Carmichael-Funktion, hier vereinfacht: φ(n))
            3. e wählen: 1 < e < φ(n), gcd(e, φ(n)) = 1  (oft e = 65537)
            4. d = e⁻¹ mod φ(n)  (privater Exponent via erweitertem eukl. Algorithmus)

        Korrektheit: m^{ed} ≡ m (mod n)  für alle m mit gcd(m,n)=1  (Euler-Satz)

    @param p Erste Primzahl.
    @param q Zweite Primzahl (p ≠ q).
    @return ((e, n), (d, n)) – (public_key, private_key).
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    n = p * q
    phi_n = (p - 1) * (q - 1)  # Eulersche Phi-Funktion für n = p*q

    # Öffentlichen Exponenten e wählen: gcd(e, phi_n) = 1
    # Standard: e = 65537; für kleine n kleinere Werte probieren
    e = 65537
    if math.gcd(e, phi_n) != 1:
        # Fallback: kleinstes e > 1 mit gcd(e, phi_n) = 1
        e = 2
        while e < phi_n and math.gcd(e, phi_n) != 1:
            e += 1

    # Privaten Exponenten d = e⁻¹ mod phi_n berechnen
    d = mod_inverse(e, phi_n)

    return (e, n), (d, n)


def rsa_encrypt(message: int, public_key: tuple[int, int]) -> int:
    """
    @brief RSA-Verschlüsselung: c = m^e mod n.
    @description
        Schnelle modulare Exponentiation via Square-and-Multiply (Python built-in pow).

    @param message Klartextnachricht als ganze Zahl (0 ≤ m < n).
    @param public_key (e, n) – öffentlicher Schlüssel.
    @return Chiffretext c.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    e, n = public_key
    # pow(m, e, n) nutzt schnelle modulare Exponentiation (Square-and-Multiply)
    return pow(message, e, n)


def rsa_decrypt(ciphertext: int, private_key: tuple[int, int]) -> int:
    """
    @brief RSA-Entschlüsselung: m = c^d mod n.
    @description
        Schnelle modulare Exponentiation via Square-and-Multiply (Python built-in pow).

    @param ciphertext Chiffretext c.
    @param private_key (d, n) – privater Schlüssel.
    @return Klartextnachricht m.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    d, n = private_key
    return pow(ciphertext, d, n)
