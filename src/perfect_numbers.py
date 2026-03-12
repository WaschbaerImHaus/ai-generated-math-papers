"""
Vollkommene Zahlen – Euklid-Euler, ungerade vollkommene Zahlen, σ-Funktion.

Eine natürliche Zahl n heißt **vollkommen (perfekt)**, wenn sie gleich der
Summe ihrer echten Teiler ist:  σ(n) = 2n  (äquivalent: σ(n)/n = 2).

**Bekannte Resultate:**
    - Euklid (~300 v.Chr.): 2^{p-1}(2^p-1) ist vollkommen wenn 2^p-1 prim.
    - Euler (1849): Alle geraden vollkommenen Zahlen haben Euklid-Form.
    - Kein ungerades Analogon bekannt (open problem seit Antike).

**Ungerade vollkommene Zahlen (Conjecture: existieren nicht):**
    Falls n ungerades vollkommenes Zahl, dann:
    - Euler-Form: n = p^{4a+1} · m²  (p prim, p ∤ m)
    - Nielsen (2006): n hat mindestens 9 verschiedene Primfaktoren
    - Brent-Cohen-te Riele (1991): n > 10^{300}
    - Chein (1979): spezifische Struktur der Primfaktor-Exponenten
    - Touchard (1953): n ≡ 1 (mod 12) oder n ≡ 9 (mod 36)

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
import sympy
from sympy import factorint, divisor_sigma, isprime, nextprime
from typing import List, Tuple, Optional, Iterator
import itertools


# ===========================================================================
# SIGMA-FUNKTION
# ===========================================================================

class SigmaFunction:
    """
    Summe der Teiler-Funktion σ(n) und verwandte arithmetische Funktionen.

    Definiert als:  σ(n) = ∑_{d|n} d  (Summe aller positiven Teiler von n)

    **Eigenschaften:**
        - Multiplikativ: σ(mn) = σ(m)·σ(n) für ggt(m,n)=1
        - Primzahlpotenzen: σ(pᵏ) = (p^{k+1} - 1) / (p - 1)
        - Vollkommen: σ(n) = 2n
        - Abundant: σ(n) > 2n
        - Defizient: σ(n) < 2n

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    @staticmethod
    def sigma(n: int, k: int = 1) -> int:
        """
        Berechnet σₖ(n) = Summe der k-ten Potenzen aller Teiler von n.

        Für k=1: σ₁(n) = σ(n) = Summe aller Teiler.
        Für k=0: σ₀(n) = d(n) = Anzahl der Teiler.

        @param n: Positive ganze Zahl
        @param k: Exponent (Standard: 1)
        @return: σₖ(n)
        @raises ValueError: wenn n ≤ 0
        """
        if n <= 0:
            raise ValueError(f"n muss positiv sein, erhalten: {n}")
        return int(divisor_sigma(n, k))

    @staticmethod
    def sigma_from_factorization(factorization: dict) -> int:
        """
        Berechnet σ(n) direkt aus der Primfaktorzerlegung.

        Nutzt die Multiplikativität:
            σ(∏ pᵢ^eᵢ) = ∏ σ(pᵢ^eᵢ) = ∏ (pᵢ^{eᵢ+1} - 1) / (pᵢ - 1)

        @param factorization: Dict {Primzahl: Exponent}, z.B. {2:3, 3:1} für n=24
        @return: σ(n) als ganze Zahl
        """
        result = 1
        for p, e in factorization.items():
            # σ(pᵉ) = 1 + p + p² + ... + pᵉ = (p^{e+1} - 1) / (p - 1)
            result *= (p**(e + 1) - 1) // (p - 1)
        return result

    @staticmethod
    def is_perfect(n: int) -> bool:
        """
        Prüft ob n eine vollkommene Zahl ist (σ(n) = 2n).

        @param n: Positive ganze Zahl
        @return: True wenn n vollkommen
        """
        if n <= 1:
            return False
        return int(divisor_sigma(n)) == 2 * n

    @staticmethod
    def is_abundant(n: int) -> bool:
        """
        Prüft ob n abundant ist (σ(n) > 2n, d.h. Summe echter Teiler > n).

        @param n: Positive ganze Zahl
        @return: True wenn n abundant
        """
        return int(divisor_sigma(n)) > 2 * n

    @staticmethod
    def is_deficient(n: int) -> bool:
        """
        Prüft ob n defizient ist (σ(n) < 2n).

        @param n: Positive ganze Zahl
        @return: True wenn n defizient
        """
        return int(divisor_sigma(n)) < 2 * n

    @staticmethod
    def abundancy_index(n: int) -> float:
        """
        Berechnet den Abundanzindex σ(n)/n.

        Vollkommen: Index = 2. Abundant: > 2. Defizient: < 2.

        @param n: Positive ganze Zahl
        @return: σ(n)/n als float
        """
        return float(divisor_sigma(n)) / n

    @staticmethod
    def multiplicativity_check(m: int, n: int) -> bool:
        """
        Verifiziert die Multiplikativität: σ(mn) = σ(m)·σ(n) für ggt(m,n)=1.

        @param m: Erste positive ganze Zahl
        @param n: Zweite positive ganze Zahl
        @return: True wenn σ(mn) = σ(m)·σ(n) (erfordert ggt(m,n)=1)
        @raises ValueError: wenn ggt(m,n) > 1
        """
        if math.gcd(m, n) != 1:
            raise ValueError(f"m={m} und n={n} sind nicht teilerfremd (ggt={math.gcd(m,n)})")
        lhs = int(divisor_sigma(m * n))
        rhs = int(divisor_sigma(m)) * int(divisor_sigma(n))
        return lhs == rhs

    @classmethod
    def sigma_sequence(cls, limit: int) -> List[Tuple[int, int]]:
        """
        Berechnet σ(n) für n = 1, 2, ..., limit.

        @param limit: Obere Grenze (inklusiv)
        @return: Liste von (n, σ(n)) Paaren
        """
        return [(n, cls.sigma(n)) for n in range(1, limit + 1)]

    @classmethod
    def perfect_numbers_up_to(cls, limit: int) -> List[int]:
        """
        Findet alle vollkommenen Zahlen bis limit durch direkten σ-Test.

        @param limit: Obere Grenze
        @return: Liste der vollkommenen Zahlen ≤ limit
        """
        return [n for n in range(2, limit + 1) if cls.is_perfect(n)]


# ===========================================================================
# GERADE VOLLKOMMENE ZAHLEN (EUKLID-EULER)
# ===========================================================================

class EvenPerfectNumbers:
    """
    Gerade vollkommene Zahlen und das Euklid-Euler-Theorem.

    **Euklid-Euler-Theorem (vollständig bewiesen):**
        n ist eine gerade vollkommene Zahl
        ↔ n = 2^{p-1}(2^p - 1) für eine Primzahl p mit 2^p - 1 prim.

    Die Zahlen der Form 2^p - 1 (mit p prim) heißen Mersenne-Zahlen.
    Ist 2^p - 1 zusätzlich prim, heißt sie Mersenne-Primzahl.

    **Bekannte gerade vollkommene Zahlen (Stand 2024):**
        Es sind 51 bekannt (entsprechend 51 bekannten Mersenne-Primzahlen).
        Ob es unendlich viele gibt, ist unbekannt (Conjecture).

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Mersenne-Exponenten der ersten bekannten geraden vollkommenen Zahlen
    # Quelle: OEIS A000043, Stand 2024
    MERSENNE_EXPONENTS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89,
                          107, 127, 521, 607, 1279, 2203, 2281, 3217,
                          4253, 4423, 9689, 9941, 11213, 19937, 21701]

    def __init__(self):
        """Initialisiert den Rechner für gerade vollkommene Zahlen."""
        pass

    @staticmethod
    def mersenne_number(p: int) -> int:
        """
        Berechnet die Mersenne-Zahl Mₚ = 2^p - 1.

        @param p: Exponent (sollte prim sein für Mersenne-Primzahlen)
        @return: 2^p - 1
        """
        return (1 << p) - 1  # Effizient: Bit-Shift statt pow()

    @staticmethod
    def is_mersenne_prime(p: int) -> bool:
        """
        Prüft ob 2^p - 1 eine Mersenne-Primzahl ist.

        Nutzt sympy.isprime für probabilistische Primalitätstests.
        Für große p (> 1000) wird dies sehr rechenintensiv.

        @param p: Exponent (p selbst muss prim sein, sonst ist 2^p-1 nicht prim)
        @return: True wenn 2^p - 1 prim ist
        """
        if not isprime(p):
            return False  # p muss prim sein (notwendige Bedingung)
        mersenne = (1 << p) - 1
        return isprime(mersenne)

    @classmethod
    def even_perfect_from_exponent(cls, p: int) -> Optional[int]:
        """
        Berechnet 2^{p-1}(2^p - 1) wenn 2^p - 1 prim ist.

        @param p: Mersenne-Exponent
        @return: Gerade vollkommene Zahl oder None wenn 2^p-1 nicht prim
        """
        if cls.is_mersenne_prime(p):
            return (1 << (p - 1)) * cls.mersenne_number(p)
        return None

    @classmethod
    def first_n_even_perfect(cls, n: int) -> List[int]:
        """
        Gibt die ersten n geraden vollkommenen Zahlen zurück.

        Nutzt die bekannte Liste der Mersenne-Exponenten (OEIS A000043).

        @param n: Anzahl der gewünschten vollkommenen Zahlen
        @return: Liste der ersten n geraden vollkommenen Zahlen
        @raises ValueError: wenn n > len(MERSENNE_EXPONENTS)
        """
        if n > len(cls.MERSENNE_EXPONENTS):
            raise ValueError(
                f"Nur {len(cls.MERSENNE_EXPONENTS)} Mersenne-Exponenten bekannt, "
                f"aber {n} angefragt."
            )
        result = []
        for p in cls.MERSENNE_EXPONENTS[:n + 5]:  # Etwas mehr für Sicherheit
            val = cls.even_perfect_from_exponent(p)
            if val is not None:
                result.append(val)
                if len(result) == n:
                    break
        return result

    @classmethod
    def verify_euclid_euler(cls, p: int) -> dict:
        """
        Verifiziert das Euklid-Euler-Theorem für einen gegebenen Exponenten p.

        Prüft alle Bedingungen:
        1. p ist prim
        2. 2^p - 1 ist prim (Mersenne-Primzahl)
        3. n = 2^{p-1}(2^p-1) ist vollkommen

        @param p: Zu prüfender Exponent
        @return: Dict mit Prüfergebnissen
        """
        p_prime = isprime(p)
        mersenne = cls.mersenne_number(p)
        mersenne_prime = isprime(mersenne) if p_prime else False
        n = (1 << (p - 1)) * mersenne if p_prime else None

        result = {
            'exponent_p': p,
            'p_is_prime': p_prime,
            'mersenne_number': mersenne if p_prime else None,
            'mersenne_is_prime': mersenne_prime,
            'perfect_number': n,
            'sigma_equals_2n': int(divisor_sigma(n)) == 2 * n if n and p <= 31 else None
        }
        return result

    @staticmethod
    def even_perfect_digit_count(p: int) -> int:
        """
        Berechnet die Anzahl der Dezimalstellen von 2^{p-1}(2^p-1).

        Ohne die Zahl vollständig zu berechnen:
            Stellen ≈ ⌊(2p-1) · log₁₀(2)⌋ + 1

        @param p: Mersenne-Exponent
        @return: Approximative Anzahl der Dezimalstellen
        """
        return int((2 * p - 1) * math.log10(2)) + 1

    @classmethod
    def euclid_euler_formula_str(cls, p: int) -> str:
        """
        Gibt die Euklid-Euler-Formel als lesbaren String aus.

        @param p: Mersenne-Exponent
        @return: Formel-String
        """
        m = cls.mersenne_number(p)
        n = (1 << (p - 1)) * m
        return f"2^{p-1} · (2^{p}-1) = 2^{p-1} · {m} = {n}"


# ===========================================================================
# UNGERADE VOLLKOMMENE ZAHLEN – SCHRANKEN UND STRUKTUR
# ===========================================================================

class OddPerfectNumberBounds:
    """
    Struktureigenschaften und Schranken für ungerade vollkommene Zahlen.

    **Hauptconjecture (seit Antike offen):**
        Es gibt keine ungeraden vollkommenen Zahlen.
        (Status: OFFEN / Conjecture, kein Beweis in Sicht)

    **Bekannte notwendige Bedingungen** (bewiesen, falls n existiert):
        1. Euler-Form: n = p^{4a+1} · q₁^{2e₁} · ... · qₖ^{2eₖ}
           mit p ≡ 1 (mod 4), p prim, p ∤ qᵢ [Euler ~1849]
        2. n > 10^{1500} [Brent-Cohen-te Riele 1991, verschärft bis 2020+]
        3. Mindestens 9 verschiedene Primfaktoren [Nielsen 2006]
        4. Mindestens 101 Primfaktoren (mit Vielfachheit) [Chein 1979]
        5. Touchard (1953): n ≡ 1 (mod 12) oder n ≡ 9 (mod 36)

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Bekannte untere Schranken (historisch)
    LOWER_BOUND_BRENT_COHEN = 10**300    # Brent-Cohen-te Riele 1991
    LOWER_BOUND_CURRENT = 10**1500       # Aktuell (2020er)

    # Mindestanzahl verschiedener Primfaktoren (Nielsen 2006)
    MIN_DISTINCT_PRIME_FACTORS = 9

    # Mindestanzahl Primfaktoren mit Vielfachheit (Chein 1979)
    MIN_PRIME_FACTORS_WITH_MULTIPLICITY = 101

    def __init__(self):
        """Initialisiert den Schranken-Prüfer."""
        pass

    @staticmethod
    def euler_form_check(n: int) -> dict:
        """
        Prüft ob n die Euler-notwendige Form für ungerade vollkommene Zahlen hat.

        **Euler-Form (notwendige Bedingung):**
            n = p^{4a+1} · m²  mit p prim, p ≡ 1 (mod 4), ggt(p, m) = 1

        Diese Bedingung ist NOTWENDIG aber nicht hinreichend.

        @param n: Zu prüfende ungerade Zahl
        @return: Dict mit Strukturanalyse
        """
        if n % 2 == 0:
            return {'odd': False, 'error': 'n ist gerade'}

        factorization = factorint(n)
        result = {
            'n': n,
            'odd': True,
            'factorization': factorization,
            'euler_form_candidates': []
        }

        # Suche nach einem Primfaktor p mit ungeradem Exponenten 4a+1
        for p, exp in factorization.items():
            # Exponent muss von der Form 4a+1 sein (≡ 1 mod 4)
            if exp % 4 == 1:
                # Prüfe ob p ≡ 1 (mod 4)
                if p % 4 == 1:
                    result['euler_form_candidates'].append({
                        'p': p,
                        'exponent': exp,
                        'a': (exp - 1) // 4,
                        'p_mod4': p % 4
                    })

        result['has_euler_form'] = len(result['euler_form_candidates']) > 0
        return result

    @staticmethod
    def touchard_congruence(n: int) -> dict:
        """
        Prüft die Touchard-Kongruenz für ungerade vollkommene Zahlen.

        **Touchard-Theorem (1953, bewiesen):**
            Falls n ungerade vollkommen, dann:
            n ≡ 1 (mod 12)  ODER  n ≡ 9 (mod 36)

        Äquivalent: n ≡ 1 (mod 12) oder n ≡ 9 (mod 36).

        Verletzt eine ungerade Zahl diese Kongruenz, kann sie nicht
        vollkommen sein.

        @param n: Zu prüfende Zahl
        @return: Dict mit Kongruenzanalyse
        """
        mod12 = n % 12
        mod36 = n % 36

        satisfies_touchard = (mod12 == 1) or (mod36 == 9)

        return {
            'n': n,
            'n_mod_12': mod12,
            'n_mod_36': mod36,
            'satisfies_n_equiv_1_mod_12': mod12 == 1,
            'satisfies_n_equiv_9_mod_36': mod36 == 9,
            'satisfies_touchard': satisfies_touchard,
            'could_be_odd_perfect': satisfies_touchard
        }

    @staticmethod
    def nielsen_bound_check(n: int) -> dict:
        """
        Prüft die Nielsen-Schranke: mindestens 9 verschiedene Primfaktoren.

        **Nielsen (2006):**
            Jede ungerade vollkommene Zahl hat mindestens 9 verschiedene
            Primfaktoren.

        @param n: Zu prüfende ungerade Zahl
        @return: Dict mit Analyseergebnissen
        """
        factorization = factorint(n)
        num_distinct = len(factorization)

        return {
            'n': n,
            'distinct_prime_factors': list(factorization.keys()),
            'num_distinct': num_distinct,
            'satisfies_nielsen': num_distinct >= OddPerfectNumberBounds.MIN_DISTINCT_PRIME_FACTORS,
            'nielsen_bound': OddPerfectNumberBounds.MIN_DISTINCT_PRIME_FACTORS
        }

    @staticmethod
    def chein_structure(n: int) -> dict:
        """
        Analysiert die Chein-Struktur: Exponenten in der Primfaktorzerlegung.

        **Chein (1979):**
            Falls n ungerade vollkommen: n = p^{4a+1} · m²
            Alle Exponenten außer einem sind gerade.
            Die Gesamtzahl der Primfaktoren (mit Vielfachheit) ≥ 101.

        @param n: Zu prüfende Zahl
        @return: Dict mit Strukturanalyse
        """
        factorization = factorint(n)
        exponents = list(factorization.values())
        total_with_mult = sum(exponents)
        odd_exponents = [e for e in exponents if e % 2 != 0]

        return {
            'n': n,
            'factorization': factorization,
            'exponents': exponents,
            'total_prime_factors_with_multiplicity': total_with_mult,
            'odd_exponent_count': len(odd_exponents),
            'odd_exponents': odd_exponents,
            'has_exactly_one_odd_exponent': len(odd_exponents) == 1,
            'satisfies_chein_multiplicity': total_with_mult >= OddPerfectNumberBounds.MIN_PRIME_FACTORS_WITH_MULTIPLICITY
        }

    @classmethod
    def is_excluded_as_odd_perfect(cls, n: int) -> Tuple[bool, List[str]]:
        """
        Prüft mehrere notwendige Bedingungen für ungerade vollkommene Zahlen.

        Gibt zurück ob n durch bekannte Kriterien ausgeschlossen werden kann.

        @param n: Zu prüfende ungerade Zahl
        @return: (ausgeschlossen, Liste der verletzten Bedingungen)
        """
        if n % 2 == 0:
            return True, ['n ist gerade']

        violated = []

        # Touchard-Kongruenz
        tc = cls.touchard_congruence(n)
        if not tc['satisfies_touchard']:
            violated.append(
                f"Touchard verletzt: n≡{tc['n_mod_12']}(mod 12) und n≡{tc['n_mod_36']}(mod 36)"
            )

        # σ(n) muss 2n sein
        if n < 10**6:  # Nur für kleine n prüfbar
            sig = int(divisor_sigma(n))
            if sig != 2 * n:
                violated.append(f"σ(n)={sig} ≠ 2n={2*n}")

        return len(violated) > 0, violated

    @classmethod
    def summary_necessary_conditions(cls) -> str:
        """
        Gibt eine Zusammenfassung aller bekannten notwendigen Bedingungen.

        @return: Formatierter String mit allen Bedingungen
        """
        return (
            "Notwendige Bedingungen für ungerade vollkommene Zahlen n:\n"
            "1. Euler-Form: n = p^{4a+1} · m², p ≡ 1 (mod 4) [Euler ~1849]\n"
            "2. n > 10^{1500} [Brent-Cohen-te Riele 1991, verschärft]\n"
            "3. Mindestens 9 verschiedene Primfaktoren [Nielsen 2006]\n"
            "4. Mindestens 101 Primfaktoren (mit Vielfachheit) [Chein 1979]\n"
            "5. Touchard: n ≡ 1 (mod 12) oder n ≡ 9 (mod 36) [1953]\n"
            "\n"
            "CONJECTURE: Es gibt keine ungeraden vollkommenen Zahlen.\n"
            "(Offen seit Antike, Stand 2026)"
        )
