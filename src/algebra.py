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

@author Kurt Ingwer
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


# =============================================================================
# RSA-KRYPTOSYSTEM
# =============================================================================

def rsa_keygen(p: int, q: int) -> tuple:
    """
    @brief Erzeugt RSA-Schlüsselpaar aus zwei Primzahlen.

    RSA (Rivest–Shamir–Adleman, 1977) ist das bekannteste Public-Key-Kryptosystem.
    Sicherheit beruht auf der Schwierigkeit der Primfaktorzerlegung großer Zahlen.

    Schlüsselerzeugung:
        1. n = p·q (Modulus)
        2. λ(n) = lcm(p-1, q-1)  (Carmichael-Funktion, hier vereinfacht: φ(n))
        3. e wählen: 1 < e < φ(n), gcd(e, φ(n)) = 1  (oft e = 65537)
        4. d = e⁻¹ mod φ(n)  (privater Exponent via erweitertem eukl. Algorithmus)

    Korrektheit: m^{ed} ≡ m (mod n)  für alle m mit gcd(m,n)=1  (Euler-Satz)

    @param p: Erste Primzahl
    @param q: Zweite Primzahl (p ≠ q)
    @return: ((e, n), (d, n)) – (public_key, private_key)
    @lastModified: 2026-03-08
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


def rsa_encrypt(message: int, public_key: tuple) -> int:
    """
    @brief RSA-Verschlüsselung: c = m^e mod n.

    @param message: Klartextnachricht als ganze Zahl (0 ≤ m < n)
    @param public_key: (e, n) – öffentlicher Schlüssel
    @return: Chiffretext c
    @lastModified: 2026-03-08
    """
    e, n = public_key
    # pow(m, e, n) nutzt schnelle modulare Exponentiation (Square-and-Multiply)
    return pow(message, e, n)


def rsa_decrypt(ciphertext: int, private_key: tuple) -> int:
    """
    @brief RSA-Entschlüsselung: m = c^d mod n.

    @param ciphertext: Chiffretext c
    @param private_key: (d, n) – privater Schlüssel
    @return: Klartextnachricht m
    @lastModified: 2026-03-08
    """
    d, n = private_key
    return pow(ciphertext, d, n)


# ===========================================================================
# DIOPHANTISCHE GLEICHUNGEN
# ===========================================================================

def solve_linear_diophantine(a: int, b: int, c: int) -> Optional[tuple]:
    """
    @brief Löst die lineare Diophantische Gleichung ax + by = c in ganzen Zahlen.
    @description
        Eine lineare Diophantische Gleichung hat ganzzahlige Lösungen genau dann,
        wenn gcd(a, b) | c (d.h. gcd(a,b) teilt c).

        Algorithmus:
        1. g = gcd(a, b) berechnen via erweitertem euklidischen Algorithmus
        2. Prüfen ob g | c
        3. Partikuläre Lösung: x0 = (c/g)*x', y0 = (c/g)*y'
           wobei a*x' + b*y' = g (aus Bezout-Identität)
        4. Allgemeine Lösung: x = x0 + (b/g)*t, y = y0 - (a/g)*t  für alle t ∈ ℤ

        Beispiel: 3x + 5y = 1
        - gcd(3,5) = 1, 1 | 1 ✓
        - Bezout: 3*2 + 5*(-1) = 1 → x0=2, y0=-1
        - Allgemein: x=2+5t, y=-1-3t

    @param a Koeffizient von x.
    @param b Koeffizient von y.
    @param c Rechte Seite.
    @return Tupel (x0, y0, gcd_ab) wobei (x0 + b/g*t, y0 - a/g*t) alle Lösungen sind,
            oder None wenn keine Lösung existiert.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    # Erweiterter Euklidischer Algorithmus liefert g, x', y' mit a*x' + b*y' = g
    g, x_bezout, y_bezout = extended_gcd(a, b)

    # Notwendige und hinreichende Bedingung: g | c
    if c % g != 0:
        return None

    # Skalierungsfaktor für partikuläre Lösung
    factor = c // g

    # Partikuläre Lösung berechnen
    x0 = factor * x_bezout
    y0 = factor * y_bezout

    return (x0, y0, g)


def solve_quadratic_diophantine_pell(D: int, n_solutions: int = 5) -> list:
    """
    @brief Löst die Pell-Gleichung x² - D·y² = 1 (D kein perfektes Quadrat).
    @description
        Die Pell-Gleichung x² - D·y² = 1 hat unendlich viele ganzzahlige Lösungen,
        falls D kein perfektes Quadrat ist.

        Fundamentallösung (x₁, y₁) via Kettenbruchentwicklung von √D:
        - √D = a₀ + 1/(a₁ + 1/(a₂ + ...)) ist eine periodische Kettenbruchentwicklung
        - Die Fundamentallösung ergibt sich aus dem (r-1)-ten Konvergenten,
          wobei r die Periode des Kettenbruchs ist
        - Alle weiteren Lösungen via Rekursion:
          x_{k+1} = x₁·x_k + D·y₁·y_k
          y_{k+1} = x₁·y_k + y₁·x_k

        Historisch: Brahmagupta (7. Jh.), Bhaskara II (12. Jh.), Lagrange (1768)

    @param D Ganzzahl, kein perfektes Quadrat.
    @param n_solutions Anzahl der zu berechnenden Lösungen.
    @return Liste von (x, y) Tupeln, erste Lösungen der Pell-Gleichung.
    @raises ValueError Wenn D ein perfektes Quadrat ist.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    # D darf kein perfektes Quadrat sein
    sqrt_D = int(math.isqrt(D))
    if sqrt_D * sqrt_D == D:
        raise ValueError(f"D={D} ist ein perfektes Quadrat, keine Pell-Gleichung")

    # Kettenbruchentwicklung von √D finden
    # a0 = floor(√D), dann periodisch
    a0 = sqrt_D

    # Konvergenten berechnen bis Pell-Gleichung erfüllt ist
    # Algorithmus: m, d, a initialisieren, dann iterieren
    m = 0
    d = 1
    a = a0

    # Konvergenten-Numerator und -Denominator berechnen
    p_prev, p_curr = 1, a0  # h_{-1}=1, h_0=a0
    q_prev, q_curr = 0, 1   # k_{-1}=0, k_0=1

    # Fundamentallösung durch Iteration finden
    for _ in range(10000):
        # Nächstes Kettenbruch-Glied berechnen
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d

        # Konvergenten aktualisieren
        p_prev, p_curr = p_curr, a * p_curr + p_prev
        q_prev, q_curr = q_curr, a * q_curr + q_prev

        # Pell-Gleichung prüfen: x²-D·y²=1?
        if p_curr * p_curr - D * q_curr * q_curr == 1:
            x1, y1 = p_curr, q_curr
            break
    else:
        raise RuntimeError(f"Fundamentallösung für D={D} nicht gefunden")

    # Alle Lösungen via Rekursionsformel erzeugen
    solutions = [(x1, y1)]
    xk, yk = x1, y1

    for _ in range(n_solutions - 1):
        # Rekursion: x_{k+1} = x1·x_k + D·y1·y_k, y_{k+1} = x1·y_k + y1·x_k
        xk_new = x1 * xk + D * y1 * yk
        yk_new = x1 * yk + y1 * xk
        xk, yk = xk_new, yk_new
        solutions.append((xk, yk))

    return solutions


def solve_pythagorean_triples(n: int) -> list:
    """
    @brief Findet alle primitiven Pythagoräischen Tripel (a, b, c) mit c ≤ n.
    @description
        Ein Pythagoräisches Tripel (a, b, c) erfüllt a² + b² = c².
        Es ist primitiv wenn gcd(a, b, c) = 1.

        Parametrische Darstellung (Euklid's Formel):
        Für m > k > 0, gcd(m,k) = 1 und m - k ungerade:
            a = m² - k²
            b = 2mk
            c = m² + k²

        Dies erzeugt alle primitiven Tripel genau einmal.

        Beispiele: (3,4,5), (5,12,13), (8,15,17), (7,24,25)

    @param n Obere Schranke für die Hypotenuse c.
    @return Sortierte Liste von (a, b, c) Tupeln mit a < b < c.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    triples = []

    # m von 2 bis ca. sqrt(n) iterieren (c = m²+k² ≤ n → m ≤ sqrt(n))
    m_max = int(math.isqrt(n)) + 1

    for m in range(2, m_max + 1):
        for k in range(1, m):
            # Bedingungen für primitive Tripel:
            # 1. gcd(m, k) = 1 (teilerfremd)
            # 2. m - k ist ungerade (verschiedene Parität)
            if gcd(m, k) != 1:
                continue
            if (m - k) % 2 == 0:
                continue

            # Euklid's Formel anwenden
            a = m * m - k * k
            b = 2 * m * k
            c = m * m + k * k

            # Schranke prüfen
            if c > n:
                break

            # a < b sicherstellen (sortiert)
            if a > b:
                a, b = b, a

            triples.append((a, b, c))

    # Nach c, dann a sortieren
    triples.sort(key=lambda t: (t[2], t[0]))
    return triples


def solve_diophantine_two_squares(n: int) -> list:
    """
    @brief Schreibt n als Summe zweier Quadrate: n = a² + b².
    @description
        Fermats Theorem der zwei Quadrate:
        Eine natürliche Zahl n ist genau dann als Summe zweier Quadrate darstellbar,
        wenn in der Primfaktorzerlegung jede Primzahl der Form 4k+3 mit gerader
        Vielfachheit vorkommt.

        Primzahlen p = 4k+1 sind immer als Summe zweier Quadrate darstellbar (Fermat).
        Primzahlen p = 4k+3 können nur in gerader Potenz auftreten.

        Algorithmus: Direktes Durchsuchen aller a von 0 bis floor(√n).
        Für jeden a prüfen ob n - a² ein perfektes Quadrat ist.

        Beispiele:
        - 25 = 0² + 5² = 3² + 4²
        - 50 = 1² + 7² = 5² + 5²
        - 3 → nicht darstellbar

    @param n Die zu zerlegende natürliche Zahl.
    @return Liste von (a, b) Tupeln mit a ≤ b und a² + b² = n, leer wenn nicht möglich.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    representations = []

    # a von 0 bis floor(√n) durchsuchen
    a_max = int(math.isqrt(n))

    for a in range(a_max + 1):
        remainder = n - a * a
        if remainder < 0:
            break

        # Prüfen ob remainder ein perfektes Quadrat ist
        b = int(math.isqrt(remainder))
        if b * b == remainder and a <= b:
            representations.append((a, b))

    return representations


def markov_numbers(n_terms: int = 20) -> list:
    """
    @brief Berechnet Markov-Zahlen (Lösungen der Markov-Gleichung x² + y² + z² = 3xyz).
    @description
        Die Markov-Gleichung x² + y² + z² = 3xyz hat unendlich viele ganzzahlige
        Lösungen, die Markov-Tripel genannt werden. Die zugehörigen Zahlen heißen
        Markov-Zahlen.

        Baumstruktur der Lösungen:
        - Fundamentallösung: (1, 1, 1)
        - Aus jedem Tripel (x, y, z) erhält man durch "Spiegelung" ein neues:
          z' = 3xy - z  (wenn x² + y² + z² = 3xyz, dann auch mit z' statt z)
        - Aus (1,1,1): → (1,1,2) → (1,2,5) → (2,5,29) → ...

        Markov-Zahlen: 1, 2, 5, 13, 29, 34, 89, 169, 194, 233, ...

        Die Markov-Eindeutigkeitsvermutung (offen seit 1879) besagt:
        Jede Markov-Zahl bestimmt ihr Tripel eindeutig.

    @param n_terms Anzahl der zu berechnenden Markov-Zahlen.
    @return Sortierte Liste der ersten n_terms Markov-Zahlen.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    # BFS-Ansatz: Alle Tripel per Breitensuche erzeugen
    from collections import deque

    markov_set = set()
    # Starttripel
    queue = deque([(1, 1, 1)])
    visited = set()
    visited.add((1, 1, 1))

    while queue and len(markov_set) < n_terms * 3:
        x, y, z = queue.popleft()

        # Alle drei Zahlen des Tripels zur Menge hinzufügen
        markov_set.add(x)
        markov_set.add(y)
        markov_set.add(z)

        # Neue Tripel durch Spiegelung erzeugen
        # Spiegeln von z: z' = 3xy - z
        candidates = [
            tuple(sorted((3 * y * z - x, y, z))),
            tuple(sorted((x, 3 * x * z - y, z))),
            tuple(sorted((x, y, 3 * x * y - z))),
        ]

        for candidate in candidates:
            if candidate not in visited:
                visited.add(candidate)
                queue.append(candidate)

    # Sortiert zurückgeben
    result = sorted(markov_set)
    return result[:n_terms]


# ===========================================================================
# QUADRATISCHES REZIPROZITÄTSGESETZ
# ===========================================================================

def is_quadratic_residue(a: int, p: int) -> bool:
    """
    @brief Prüft ob a ein quadratischer Rest modulo der Primzahl p ist.
    @description
        a ist ein quadratischer Rest mod p wenn es ein x gibt mit x² ≡ a (mod p).

        Euler-Kriterium: a ist QR mod p ⟺ a^((p-1)/2) ≡ 1 (mod p)
        (für p ungerade Primzahl und p ∤ a)

        Das Euler-Kriterium folgt aus dem Kleinen Fermatschen Satz:
        a^(p-1) ≡ 1 (mod p) → (a^((p-1)/2))² ≡ 1 → a^((p-1)/2) ≡ ±1

        Falls a ≡ 0 (mod p): kein QR (Trivialfall)

    @param a Die zu prüfende Zahl.
    @param p Ungerade Primzahl.
    @return True wenn a ein quadratischer Rest mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    a = a % p
    # 0 ist kein quadratischer Rest im üblichen Sinne
    if a == 0:
        return False
    # Euler-Kriterium: a^((p-1)/2) ≡ 1 (mod p)
    return pow(a, (p - 1) // 2, p) == 1


def quadratic_residues(p: int) -> list:
    """
    @brief Gibt alle quadratischen Reste modulo p zurück.
    @description
        Die quadratischen Reste mod p sind die Elemente von {1,...,p-1}
        die ein Quadrat modulo p sind.

        Eigenschaft: Es gibt genau (p-1)/2 quadratische Reste mod p (für Primzahl p).
        Diese sind die Werte {1², 2², ..., ((p-1)/2)²} mod p (alle verschieden).

        Beispiel p=7: 1²=1, 2²=4, 3²=2 mod 7 → QR={1,2,4}

    @param p Ungerade Primzahl.
    @return Sortierte Liste der quadratischen Reste mod p.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    residues = set()
    # Nur bis (p-1)/2 iterieren nötig, da k² ≡ (p-k)² mod p
    for i in range(1, p):
        residues.add((i * i) % p)
    # 0 entfernen falls vorhanden
    residues.discard(0)
    return sorted(residues)


def quadratic_reciprocity(p: int, q: int) -> dict:
    """
    @brief Wendet das Quadratische Reziprozitätsgesetz an.
    @description
        Das Quadratische Reziprozitätsgesetz (Gauß, 1796) beschreibt die Beziehung
        zwischen dem Legendre-Symbol (p|q) und (q|p):

            (p/q) · (q/p) = (-1)^{((p-1)/2) · ((q-1)/2)}

        Das bedeutet:
        - Wenn p ≡ 1 (mod 4) oder q ≡ 1 (mod 4): (p/q) = (q/p)
        - Wenn p ≡ q ≡ 3 (mod 4): (p/q) = -(q/p)

        Das Legendre-Symbol (a/p) = 1 wenn a QR mod p, -1 wenn NQR mod p, 0 wenn p|a.

    @param p Ungerade Primzahl.
    @param q Ungerade Primzahl (p ≠ q).
    @return Dict mit Legendre-Symbolen, Produkt, erwartetem Vorzeichen und Verifikation.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    def legendre_symbol(a: int, prime: int) -> int:
        """Berechnet das Legendre-Symbol (a/prime) via Euler-Kriterium."""
        a = a % prime
        if a == 0:
            return 0
        val = pow(a, (prime - 1) // 2, prime)
        # val ist 1 oder prime-1 (≡ -1)
        return 1 if val == 1 else -1

    # Legendre-Symbole berechnen
    lpq = legendre_symbol(p, q)  # (p/q)
    lqp = legendre_symbol(q, p)  # (q/p)

    # Vorzeichen gemäß Reziprozitätsgesetz
    exponent = ((p - 1) // 2) * ((q - 1) // 2)
    expected_sign = (-1) ** exponent

    # Produkt der beiden Legendre-Symbole
    product = lpq * lqp

    # Gesetz verifizieren
    law_verified = (product == expected_sign)

    return {
        'legendre_p_q': lpq,
        'legendre_q_p': lqp,
        'product': product,
        'expected_sign': expected_sign,
        'law_verified': law_verified
    }


def tonelli_shanks(n: int, p: int) -> Optional[int]:
    """
    @brief Berechnet √n (mod p) via Tonelli-Shanks-Algorithmus.
    @description
        Der Tonelli-Shanks-Algorithmus findet x mit x² ≡ n (mod p).
        Er funktioniert für beliebige ungerade Primzahlen p.

        Spezialfall: p ≡ 3 (mod 4): x = n^((p+1)/4) mod p (direkte Formel)

        Allgemeiner Fall (Tonelli-Shanks):
        1. Schreibe p-1 = Q · 2^S mit Q ungerade
        2. Finde z mit z QNR mod p (quadratischer Nicht-Rest)
        3. Initialisiere M=S, c=z^Q, t=n^Q, R=n^((Q+1)/2) mod p
        4. Iteriere: Reduziere M bis Konvergenz

    @param n Zahl deren Quadratwurzel gesucht wird.
    @param p Ungerade Primzahl.
    @return x mit x² ≡ n (mod p), oder None wenn n kein QR mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    n = n % p

    # Triviale Fälle
    if n == 0:
        return 0
    if p == 2:
        return n

    # Kein quadratischer Rest?
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Spezialfall: p ≡ 3 (mod 4) → direkte Formel
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Allgemeiner Tonelli-Shanks
    # Schritt 1: p-1 = Q * 2^S zerlegen
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # Schritt 2: Quadratischen Nicht-Rest z finden
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    # Schritt 3: Initialisierung
    M = S
    c = pow(z, Q, p)       # c = z^Q mod p
    t = pow(n, Q, p)       # t = n^Q mod p
    R = pow(n, (Q + 1) // 2, p)  # R = n^((Q+1)/2) mod p

    # Schritt 4: Hauptschleife
    while True:
        if t == 0:
            return 0
        if t == 1:
            return R

        # i finden: kleinste i > 0 mit t^(2^i) ≡ 1 (mod p)
        i = 1
        temp = (t * t) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1

        # Update
        b = pow(c, pow(2, M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p


def cipolla_algorithm(n: int, p: int) -> Optional[int]:
    """
    @brief Berechnet √n (mod p) via Cipolla-Algorithmus.
    @description
        Der Cipolla-Algorithmus (Michele Cipolla, 1903) ist eine Alternative
        zu Tonelli-Shanks zur Berechnung von Quadratwurzeln modulo Primzahlen.

        Algorithmus:
        1. Finde a, so dass a² - n ein QNR (quadratischer Nicht-Rest) mod p ist
        2. Arbeite in GF(p²) = GF(p)[x]/(x² - (a²-n))
           (erweitertes Körper über GF(p) mit irreduziblem Element ω mit ω² = a²-n)
        3. Berechne (a + ω)^((p+1)/2) in GF(p²)
        4. Das Ergebnis hat imaginären Teil 0 und reellen Teil = √n mod p

    @param n Zahl deren Quadratwurzel gesucht wird.
    @param p Ungerade Primzahl.
    @return x mit x² ≡ n (mod p), oder None wenn n kein QR mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-08
    """
    n = n % p

    # Trivialfälle
    if n == 0:
        return 0
    if p == 2:
        return n

    # n muss quadratischer Rest sein
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Schritt 1: a finden so dass a²-n ein QNR ist
    # (a²-n ist QNR ⟺ (a²-n)^((p-1)/2) ≡ -1 ≡ p-1 mod p)
    a = 0
    omega_sq = 0  # ω² = a² - n
    for a in range(p):
        omega_sq = (a * a - n) % p
        if omega_sq == 0:
            # a² ≡ n mod p → a ist direkt die Wurzel
            return a
        if pow(omega_sq, (p - 1) // 2, p) == p - 1:
            # omega_sq ist QNR → gefunden
            break
    else:
        return None  # Sollte nicht vorkommen

    # Schritt 2+3: Potenzierung in GF(p²)
    # Elemente der Form (x, y) entsprechen x + y*ω in GF(p²)
    # Multiplikation: (x1+y1·ω)·(x2+y2·ω) = (x1x2 + y1y2·ω²) + (x1y2+x2y1)·ω

    def mul_gf(x1: int, y1: int, x2: int, y2: int) -> tuple:
        """Multiplikation in GF(p²)."""
        real = (x1 * x2 + y1 * y2 * omega_sq) % p
        imag = (x1 * y2 + x2 * y1) % p
        return real, imag

    # Schnelle Potenzierung: (a + ω)^((p+1)/2) mod p
    exp = (p + 1) // 2
    result_real, result_imag = 1, 0  # neutrale Element = 1
    base_real, base_imag = a % p, 1  # Basis = a + 1*ω

    while exp > 0:
        if exp % 2 == 1:
            result_real, result_imag = mul_gf(result_real, result_imag, base_real, base_imag)
        base_real, base_imag = mul_gf(base_real, base_imag, base_real, base_imag)
        exp //= 2

    # Imaginärteil muss 0 sein (Ergebnis liegt in GF(p))
    # result_real ist die gesuchte Quadratwurzel
    return result_real
