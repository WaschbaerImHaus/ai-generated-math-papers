"""
@file algebra.py
@brief Algebra-Modul: Polynome, Gleichungslöser, Zahlentheorie.
@description
    Implementiert grundlegende und fortgeschrittene algebraische Konzepte:

    - Polynomial: Klasse für Polynomoperationen (Addition, Multiplikation,
      Ableitung, Auswertung, Darstellung)
    - Gleichungslöser: linear, quadratisch
    - Zahlentheorie: ggT, kgV, erweiterter euklidischer Algorithmus,
      modulares Inverses, Primzahltest, Primfaktorzerlegung, Euler-Phi

    Alle Algorithmen sind mit mathematischen Erklärungen dokumentiert.

@author Reisen macht Spass... mit Pia und Dirk e.Kfm.
@date 2026-03-05
"""

import math
import cmath
from typing import Union, Optional


# =============================================================================
# KLASSE: Polynomial
# =============================================================================

class Polynomial:
    """
    @brief Repräsentation und Operationen mit Polynomen.
    @description
        Ein Polynom der Form:
            p(x) = a_n * x^n + a_{n-1} * x^{n-1} + ... + a_1 * x + a_0

        Die Koeffizienten werden als Liste gespeichert, wobei Index 0 dem
        höchsten Grad entspricht:
            coefficients[0] * x^n + coefficients[1] * x^{n-1} + ... + coefficients[n]

    @example
        p = Polynomial([1, -3, 2])  # x^2 - 3x + 2
        p.evaluate(1)  # ergibt 0 (Nullstelle)
    @date 2026-03-05
    """

    def __init__(self, coefficients: list):
        """
        @brief Erstellt ein Polynom aus einer Koeffizientenliste.
        @param coefficients Liste von Koeffizienten (höchster Grad zuerst).
                           Beispiel: [1, -3, 2] für x^2 - 3x + 2
        @date 2026-03-05
        """
        if not coefficients:
            # Leere Liste: Nullpolynom
            self.coefficients = [0]
        else:
            self.coefficients = list(coefficients)

        # Führende Nullen entfernen (Normalisierung)
        while len(self.coefficients) > 1 and self.coefficients[0] == 0:
            self.coefficients.pop(0)

    @property
    def degree(self) -> int:
        """
        @brief Gibt den Grad des Polynoms zurück.
        @description
            Der Grad ist der höchste Exponent mit einem Koeffizient != 0.
            Das Nullpolynom hat per Konvention Grad -1 (oder 0 je nach Kontext).
        @return Grad des Polynoms.
        @date 2026-03-05
        """
        if len(self.coefficients) == 1 and self.coefficients[0] == 0:
            return 0  # Nullpolynom
        return len(self.coefficients) - 1

    def evaluate(self, x: float) -> float:
        """
        @brief Wertet das Polynom an der Stelle x aus (Horner-Schema).
        @description
            Das Horner-Schema ist eine effiziente Methode zur Polynomauswertung.
            Statt x^n direkt zu berechnen, werden die Berechnungen verschachtelt:
                a*x^3 + b*x^2 + c*x + d = ((a*x + b)*x + c)*x + d

            Das reduziert die Anzahl der Multiplikationen von O(n^2) auf O(n).

        @param x Der Auswertungspunkt.
        @return Wert des Polynoms an der Stelle x.
        @date 2026-03-05
        """
        # Horner-Schema: Von links beginnend akkumulieren
        result = 0
        for coeff in self.coefficients:
            result = result * x + coeff
        return result

    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """
        @brief Addition zweier Polynome.
        @description
            Polynome werden komponentenweise addiert. Unterschiedliche Grade
            werden durch Auffüllen mit Nullen ausgeglichen.
        @param other Das zu addierende Polynom.
        @return Summe der beiden Polynome.
        @date 2026-03-05
        """
        # Längere Liste bestimmt den Grad des Ergebnisses
        n1, n2 = len(self.coefficients), len(other.coefficients)
        # Kürzere Liste vorne mit Nullen auffüllen
        a = [0] * max(0, n2 - n1) + self.coefficients
        b = [0] * max(0, n1 - n2) + other.coefficients
        result = [x + y for x, y in zip(a, b)]
        return Polynomial(result)

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """
        @brief Subtraktion zweier Polynome.
        @param other Das abzuziehende Polynom.
        @return Differenz der beiden Polynome.
        @date 2026-03-05
        """
        # Negation des zweiten Polynoms, dann addieren
        neg_other = Polynomial([-c for c in other.coefficients])
        return self + neg_other

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        @brief Multiplikation zweier Polynome (Faltungsoperation).
        @description
            Das Produkt p(x) * q(x) hat Grad deg(p) + deg(q).
            Jeder Koeffizient des Produkts ist:
                c_k = Summe von a_i * b_j für alle i+j=k

        @param other Das zu multiplizierende Polynom.
        @return Produkt der beiden Polynome.
        @date 2026-03-05
        """
        n1 = len(self.coefficients)
        n2 = len(other.coefficients)
        # Ergebnis hat n1 + n2 - 1 Koeffizienten
        result = [0] * (n1 + n2 - 1)
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                result[i + j] += a * b
        return Polynomial(result)

    def derivative(self) -> 'Polynomial':
        """
        @brief Berechnet die Ableitung des Polynoms.
        @description
            Die Ableitung nutzt die Potenzregel:
                d/dx (a * x^n) = a * n * x^(n-1)

            Beispiel: d/dx (x^3 - 2x + 1) = 3x^2 - 2

        @return Ableitungspolynom.
        @date 2026-03-05
        """
        if self.degree == 0:
            return Polynomial([0])

        # Für jeden Koeffizient: Multiplikation mit dem Exponenten
        new_coeffs = []
        n = self.degree
        for i, coeff in enumerate(self.coefficients[:-1]):  # letzten weglassen (Konstante wird 0)
            # Exponent des Terms: n - i
            exponent = n - i
            new_coeffs.append(coeff * exponent)
        return Polynomial(new_coeffs)

    def __str__(self) -> str:
        """
        @brief Gibt das Polynom als lesbaren String aus.
        @description
            Formatiert das Polynom in mathematischer Notation, z.B.:
            [1, -3, 2] -> "x^2 - 3x + 2"
        @return String-Darstellung des Polynoms.
        @date 2026-03-05
        """
        if not self.coefficients or all(c == 0 for c in self.coefficients):
            return "0"

        terms = []
        n = self.degree

        for i, coeff in enumerate(self.coefficients):
            if coeff == 0:
                continue
            exp = n - i  # Exponent dieses Terms

            # Vorzeichen bestimmen
            sign = "+" if coeff > 0 and terms else ""
            if coeff < 0:
                sign = "-"
                coeff = -coeff

            # Term aufbauen
            if exp == 0:
                terms.append(f"{sign}{coeff}")
            elif exp == 1:
                coeff_str = "" if coeff == 1 else str(coeff)
                terms.append(f"{sign}{coeff_str}x")
            else:
                coeff_str = "" if coeff == 1 else str(coeff)
                terms.append(f"{sign}{coeff_str}x^{exp}")

        return " ".join(terms) if terms else "0"

    def __repr__(self) -> str:
        return f"Polynomial({self.coefficients})"


# =============================================================================
# GLEICHUNGSLÖSER
# =============================================================================

def solve_linear(a: float, b: float) -> float:
    """
    @brief Löst eine lineare Gleichung ax + b = 0.
    @description
        Lösung: x = -b / a

        Sonderfälle:
        - a = 0, b != 0: keine Lösung (Widerspruch)
        - a = 0, b = 0: unendlich viele Lösungen

    @param a Koeffizient von x.
    @param b Konstanter Term.
    @return Lösung x.
    @raises ValueError Wenn keine oder unendlich viele Lösungen existieren.
    @date 2026-03-05
    """
    if a == 0:
        if b == 0:
            raise ValueError("unendlich viele Lösungen: 0 = 0 ist immer wahr")
        else:
            raise ValueError(f"keine Lösung: {b} = 0 ist ein Widerspruch")
    return -b / a


def solve_quadratic(a: float, b: float, c: float) -> list:
    """
    @brief Löst eine quadratische Gleichung ax^2 + bx + c = 0.
    @description
        Verwendet die Lösungsformel (pq-Formel / Mitternachtsformel):
            x_{1,2} = (-b ± sqrt(b^2 - 4ac)) / (2a)

        Die Diskriminante D = b^2 - 4ac entscheidet über die Art der Wurzeln:
        - D > 0: zwei verschiedene reelle Wurzeln
        - D = 0: eine doppelte reelle Wurzel
        - D < 0: zwei komplexkonjugierte Wurzeln

    @param a Koeffizient von x^2 (darf nicht 0 sein).
    @param b Koeffizient von x.
    @param c Konstanter Term.
    @return Liste mit zwei Wurzeln (reell oder komplex).
    @raises ValueError Wenn a = 0.
    @date 2026-03-05
    """
    if a == 0:
        raise ValueError("a darf nicht 0 sein (dann ist es keine quadratische Gleichung)")

    # Diskriminante berechnen
    discriminant = b**2 - 4 * a * c

    if discriminant >= 0:
        # Reelle Wurzeln (sqrt(D) ist reell)
        sqrt_d = math.sqrt(discriminant)
        x1 = (-b + sqrt_d) / (2 * a)
        x2 = (-b - sqrt_d) / (2 * a)
    else:
        # Komplexe Wurzeln (sqrt(D) ist imaginär)
        sqrt_d = cmath.sqrt(discriminant)
        x1 = (-b + sqrt_d) / (2 * a)
        x2 = (-b - sqrt_d) / (2 * a)

    return [x1, x2]


# =============================================================================
# ZAHLENTHEORIE
# =============================================================================

def gcd(a: int, b: int) -> int:
    """
    @brief Berechnet den größten gemeinsamen Teiler (ggT) mit dem Euklidischen Algorithmus.
    @description
        Der Euklidische Algorithmus basiert auf der Eigenschaft:
            ggT(a, b) = ggT(b, a mod b)

        Terminiert wegen abnehmender Reste (a mod b < b).
        Laufzeit: O(log(min(a,b)))

        Beispiel: ggT(48, 18)
            = ggT(18, 12)   [48 = 2*18 + 12]
            = ggT(12, 6)    [18 = 1*12 + 6]
            = ggT(6, 0)     [12 = 2*6 + 0]
            = 6

    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return ggT(a, b), immer >= 0.
    @date 2026-03-05
    """
    # Negative Zahlen behandeln
    a, b = abs(a), abs(b)
    # Euklidischer Algorithmus
    while b:
        a, b = b, a % b
    return a


def lcm(a: int, b: int) -> int:
    """
    @brief Berechnet das kleinste gemeinsame Vielfache (kgV).
    @description
        Verwendet die Beziehung:
            kgV(a, b) = |a * b| / ggT(a, b)

        Diese Formel verhindert Integer-Überlauf durch frühe Division:
            kgV(a, b) = a // ggT(a, b) * b

    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return kgV(a, b).
    @date 2026-03-05
    """
    if a == 0 or b == 0:
        return 0
    return abs(a) // gcd(a, b) * abs(b)


def extended_gcd(a: int, b: int) -> tuple:
    """
    @brief Erweiterter Euklidischer Algorithmus: findet x, y mit ax + by = ggT(a,b).
    @description
        Der erweiterte Euklidische Algorithmus löst die Bezout-Gleichung:
            a*x + b*y = ggT(a, b)

        Er ist die Grundlage für modulare Inverse und RSA-Kryptographie.

        Rekursive Formulierung:
            extended_gcd(a, 0) = (a, 1, 0)
            extended_gcd(a, b) = (g, y, x - (a//b)*y)
            wobei (g, x, y) = extended_gcd(b, a%b)

    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return Tupel (g, x, y) mit g=ggT(a,b) und a*x + b*y = g.
    @date 2026-03-05
    """
    if b == 0:
        # Basisfall: ggT(a,0) = a, mit a*1 + 0*0 = a
        return a, 1, 0

    # Rekursiver Schritt
    g, x, y = extended_gcd(b, a % b)
    # Rücksubstitution der Koeffizienten
    return g, y, x - (a // b) * y


def mod_inverse(a: int, m: int) -> int:
    """
    @brief Berechnet das modulare Inverse von a modulo m.
    @description
        Das modulare Inverse von a (mod m) ist die Zahl x mit:
            a * x ≡ 1 (mod m)

        Existiert genau dann, wenn ggT(a, m) = 1 (also a und m teilerfremd sind).
        Berechnet mit dem erweiterten Euklidischen Algorithmus.

        Anwendung: RSA-Entschlüsselung, Modulo-Division

    @param a Die Zahl, deren Inverses gesucht wird.
    @param m Der Modul.
    @return x mit (a * x) % m == 1.
    @raises ValueError Wenn ggT(a, m) != 1 (kein Inverses existiert).
    @date 2026-03-05
    """
    g, x, _ = extended_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"Kein modulares Inverses: ggT({a}, {m}) = {g} ≠ 1")
    # x könnte negativ sein, auf positiven Bereich normieren
    return x % m


def is_prime(n: int) -> bool:
    """
    @brief Prüft, ob eine natürliche Zahl eine Primzahl ist.
    @description
        Ein Primzahltest durch Probedivision mit Optimierungen:
        1. Sonderfälle: n < 2 ist keine Primzahl; 2 und 3 sind Primzahlen
        2. Zahlen der Form 6k ± 1 überprüfen (alle Primzahlen > 3 haben diese Form)
        3. Nur bis sqrt(n) testen (wenn n einen Teiler > sqrt(n) hat,
           hat es auch einen < sqrt(n))

        Laufzeit: O(sqrt(n))

    @param n Die zu prüfende natürliche Zahl.
    @return True wenn n eine Primzahl ist, sonst False.
    @date 2026-03-05
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


def prime_factorization(n: int) -> dict:
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

    @param n Die zu zerlegende natürliche Zahl (n >= 2).
    @return Dictionary {primzahl: exponent}.
    @date 2026-03-05
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

    @param n Positive ganze Zahl.
    @return phi(n) - Anzahl der teilerfremden Zahlen.
    @date 2026-03-05
    """
    if n == 1:
        return 1

    # Primfaktoren bestimmen (eindeutig, ohne Exponenten)
    result = n
    temp = n

    # Für jeden Primfaktor p: result = result * (1 - 1/p) = result * (p-1)/p
    p = 2
    while p * p <= temp:
        if temp % p == 0:
            # p ist ein Primfaktor
            while temp % p == 0:
                temp //= p
            # phi-Formel anwenden
            result -= result // p
        p += 1

    # Letzter verbleibender Primfaktor > sqrt(n)
    if temp > 1:
        result -= result // temp

    return result
