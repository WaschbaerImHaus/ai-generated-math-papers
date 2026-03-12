"""
Waring-Goldbach-Problem – Darstellung von Zahlen als Summen von Primzahlpotenzen.

Das Waring-Goldbach-Problem fragt: Wie viele Potenzen k-ter Ordnung von
Primzahlen genügen, um jede hinreichend große natürliche Zahl darzustellen?

**Verallgemeinerte Frage:**
    G(k) = Minimale Anzahl k-ter Potenzen von Primzahlen, die jede
           hinreichend große ganze Zahl darstellen.
    g(k) = Minimale Anzahl k-ter Potenzen von Primzahlen für ALLE ganzen Zahlen.

**Bekannte Resultate:**
    k=1 (Vinogradov 1937): Jede hinreichend große ungerade Zahl = p₁+p₂+p₃.
        → Goldbach-Vermutung (Conjecture): Jede gerade Zahl ≥ 4 = p₁+p₂.
    k=2 (Hua 1938): Jede hinreichend große Zahl = p₁²+...+p₅².
        G(2) ≤ 5, und 5 ist optimal für bestimmte Restklassen.
    k=3: G(3) ≤ 9 (Hua 1938, verbessert durch Vaughan-Wooley).
    k=4: G(4) ≤ 15 (Hua).
    Allgemein: G(k) ≤ k(3 log k + 4.2 log log k + C) (Vinogradov-Wooley).

**Lokale Bedingungen (Singular Series):**
    Nicht alle Zahlen sind als Summen von k Primzahlquadraten darstellbar.
    Für k=2: Nur Zahlen ≠ 0, 4, 7 (mod 8) haben Darstellungen.
    Die singuläre Reihe S(n) beschreibt die Anzahl der Darstellungen asymptotisch.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
import numpy as np
import sympy
from sympy import isprime, primerange, nextprime
from typing import List, Optional, Tuple, Dict
import itertools


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _sieve_primes(limit: int) -> List[int]:
    """
    Sieb des Eratosthenes bis limit.

    @param limit: Obere Grenze (inklusiv)
    @return: Sortierte Liste aller Primzahlen ≤ limit
    """
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i, v in enumerate(sieve) if v]


# ===========================================================================
# WARING-GOLDBACH k=2 (PRIMZAHLQUADRATE)
# ===========================================================================

class WaringGoldbachK2:
    """
    Waring-Goldbach-Problem für k=2: Darstellung als Summen von Primzahlquadraten.

    **Hua Looeng Keng (1938):**
        Jede hinreichend große natürliche Zahl n ist Summe von ≤ 5 Primzahlquadraten,
        sofern n die lokalen Bedingungen erfüllt.

    **Lokale Bedingungen:**
        n muss ≡ 5 (mod 24) sein für die 5-Quadrate-Darstellung in voller Allgemeinheit.
        (Genauer: Jede n ≥ N₀ mit n ≠ 0, 4, 7 (mod 8) und gewissen mod 3 Bedingungen.)

    **G(2) = 5** (Hua, exakt unter den lokalen Bedingungen).

    **Hinweis:** Nicht jede Zahl ist als Summe von 4 Primzahlquadraten darstellbar;
    die Darstellung mit ≤ 5 Termen gilt für hinreichend große Zahlen.

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Hua's Theorem: Jede hinreichend große Zahl = Summe von ≤ G(2) Primzahlquadraten
    HUA_G2 = 5  # G(2) ≤ 5

    def __init__(self, max_prime: int = 400):
        """
        Initialisiert den Rechner mit Primzahlen bis max_prime.

        @param max_prime: Obergrenze für Primzahlen (Standard: 400, Quadrat bis 160000)
        """
        self.max_prime = max_prime
        self.primes = _sieve_primes(max_prime)
        # Primzahlquadrate vorberechnen
        self.prime_squares = [p * p for p in self.primes]

    def local_conditions(self, n: int) -> dict:
        """
        Analysiert die lokalen Bedingungen für die Darstellung von n.

        Lokale Bedingungen aus der Kreismethode:
        - n mod 8: Bestimmt welche Kongruenzklassen darstellbar sind
        - n mod 3: Zusätzliche Einschränkungen

        @param n: Zu analysierende Zahl
        @return: Dict mit lokalen Bedingungen
        """
        return {
            'n': n,
            'n_mod_8': n % 8,
            'n_mod_24': n % 24,
            'n_mod_3': n % 3,
            'local_obstruction_mod8': n % 8 in [0, 4, 7],
            'comment': (
                'Lokal darstellbar (keine offensichtliche Hürde)'
                if n % 8 not in [0, 4, 7]
                else 'Möglicherweise lokale Obstruktion (mod 8)'
            )
        }

    def represent_as_prime_squares(self, n: int, max_terms: int = 5) -> Optional[List[int]]:
        """
        Findet eine Darstellung n = p₁² + p₂² + ... + pₖ² mit k ≤ max_terms.

        Verwendet BFS/dynamische Suche über Primzahlquadrate.

        @param n: Darzustellende Zahl (n > 0)
        @param max_terms: Maximale Anzahl der Summanden
        @return: Liste [p₁, p₂, ..., pₖ] oder None wenn keine Darstellung gefunden
        """
        if n <= 0:
            return None

        # Filtere relevante Primzahlquadrate
        relevant_sq = [sq for sq in self.prime_squares if sq <= n]
        if not relevant_sq:
            return None

        def search(remaining: int, terms_left: int, min_idx: int) -> Optional[List[int]]:
            """Rekursive Suche mit Pruning."""
            if remaining == 0:
                return []
            if terms_left == 0:
                return None
            # Pruning: Mit terms_left Termen können wir höchstens terms_left * max(sq) erreichen
            for i in range(min_idx, len(relevant_sq)):
                sq = relevant_sq[i]
                if sq > remaining:
                    break  # Alle weiteren Quadrate sind größer
                result = search(remaining - sq, terms_left - 1, i)
                if result is not None:
                    # Primzahl aus Quadrat rückgewinnen
                    return [self.primes[i]] + result
            return None

        return search(n, max_terms, 0)

    def verify_hua_theorem(self, limit: int = 1000) -> dict:
        """
        Verifiziert Hua's Theorem numerisch für alle n ≤ limit.

        Überprüft ob jede Zahl n mit n ≡ 5 (mod 24) und n ≤ limit
        als Summe von ≤ 5 Primzahlquadraten darstellbar ist.

        @param limit: Obere Prüfgrenze
        @return: Dict mit Statistiken
        """
        successes = 0
        failures = []
        tested = 0

        for n in range(4, limit + 1):
            rep = self.represent_as_prime_squares(n, max_terms=5)
            tested += 1
            if rep is not None:
                successes += 1
            else:
                failures.append(n)

        return {
            'tested': tested,
            'successes': successes,
            'failures': failures,
            'failure_count': len(failures),
            'success_rate': successes / tested if tested > 0 else 0
        }

    def count_representations(self, n: int, num_terms: int = 5) -> int:
        """
        Zählt die Anzahl der Darstellungen von n als Summe von num_terms Primzahlquadraten.

        Darstellungen werden als GEORDNETE Mengen gezählt (p₁ ≤ p₂ ≤ ... ≤ pₖ).

        @param n: Darzustellende Zahl
        @param num_terms: Anzahl der Summanden
        @return: Anzahl der Darstellungen
        """
        if num_terms == 0:
            return 1 if n == 0 else 0

        relevant_sq = [sq for sq in self.prime_squares if sq <= n]
        count = 0

        def count_recursive(remaining: int, terms_left: int, min_idx: int) -> int:
            if remaining == 0 and terms_left == 0:
                return 1
            if terms_left == 0 or remaining == 0:
                return 0
            total = 0
            for i in range(min_idx, len(relevant_sq)):
                sq = relevant_sq[i]
                if sq > remaining:
                    break
                total += count_recursive(remaining - sq, terms_left - 1, i)
            return total

        return count_recursive(n, num_terms, 0)

    def minimal_representation(self, n: int) -> Tuple[int, Optional[List[int]]]:
        """
        Findet die minimale Anzahl von Primzahlquadraten für n.

        @param n: Darzustellende Zahl
        @return: (Minimale Anzahl, Liste der Primzahlen)
        """
        for k in range(1, self.HUA_G2 + 2):
            rep = self.represent_as_prime_squares(n, max_terms=k)
            if rep is not None:
                return k, rep
        return -1, None  # Nicht darstellbar mit ≤ HUA_G2+1 Termen

    def singular_series_heuristic(self, n: int) -> float:
        """
        Heuristische Schätzung der singulären Reihe S(n) für k=2.

        Die singuläre Reihe beschreibt die Anzahl der Darstellungen:
            R₅(n) ~ S(n) · n^{3/2} / (log n)^5  (für n → ∞)

        @param n: Zu analysierende Zahl
        @return: Heuristischer Wert von S(n)
        """
        if n <= 1:
            return 0.0

        # Vereinfachte Euler-Produkt-Schätzung
        # S(n) = ∏_p (1 + χ_p(n) / (p-1)^5)  (grobe Näherung)
        s = 1.0
        for p in _sieve_primes(min(50, n)):
            n_mod_p = n % p
            # Einfache lokale Dichte
            local_density = sum(1 for a in range(p) if isprime(a) or a == 0)
            s *= (1 + local_density / (p * 4))
        return max(0.0, s)


# ===========================================================================
# WARING-GOLDBACH ALLGEMEIN (k≥3)
# ===========================================================================

class WaringGoldbachGeneral:
    """
    Allgemeines Waring-Goldbach-Problem für beliebige Potenzen k ≥ 3.

    **Bekannte Schranken für G(k)** (Primzahlpotenzen, alle hinr. große Zahlen):
        k=2:  G(2) ≤ 5   (Hua 1938)
        k=3:  G(3) ≤ 9   (Hua 1938, Vaughan-Wooley verschärft)
        k=4:  G(4) ≤ 15  (Hua 1938)
        k=5:  G(5) ≤ 21  (Hua 1938)
        k=6:  G(6) ≤ 33  (Hua 1938)
        General: G(k) ≤ k·(3 log k + 4.2 log log k + C)  [Vinogradov-Wooley]

    **Wooley (1992-2019, Efficient Congruencing):**
        Signifikante Verbesserungen für große k durch Vaughan-Wooley-Methode.

    **Statt G(k) = Schranke für hinr. große Zahlen:**
        g(k) = Schranke für ALLE Zahlen (schwieriger, größere Werte).

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Bekannte obere Schranken G(k) für hinreichend große Zahlen (Hua 1938)
    G_BOUNDS_HUA: Dict[int, int] = {
        1: 3,   # Vinogradov 1937 (für ungerade n)
        2: 5,   # Hua 1938
        3: 9,   # Hua 1938
        4: 15,  # Hua 1938
        5: 21,  # Hua 1938
        6: 33,  # Hua 1938
        7: 45,  # Hua 1938
        8: 63,  # Hua 1938
    }

    # Verbesserte Schranken (Vaughan-Wooley, Efficient Congruencing)
    G_BOUNDS_IMPROVED: Dict[int, int] = {
        2: 5,
        3: 7,   # Vinogradov-Wooley (für k=3 gibt es bessere Resultate)
        4: 13,
        5: 15,
        6: 19,
    }

    def __init__(self, k: int, max_prime: int = 200):
        """
        Initialisiert den Rechner für Waring-Goldbach der Ordnung k.

        @param k: Potenzordnung (k ≥ 1)
        @param max_prime: Obergrenze für Primzahlen in Berechnungen
        @raises ValueError: wenn k < 1
        """
        if k < 1:
            raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")

        self.k = k
        self.max_prime = max_prime
        self.primes = _sieve_primes(max_prime)
        # k-te Potenzen der Primzahlen vorberechnen
        self.prime_powers = [p**k for p in self.primes]

    def g_bound_hua(self) -> Optional[int]:
        """
        Gibt die Hua (1938) Schranke G(k) zurück (falls bekannt).

        @return: Obere Schranke G(k) oder None wenn k > 8
        """
        return self.G_BOUNDS_HUA.get(self.k)

    def g_bound_vinogradov_wooley(self) -> int:
        """
        Gibt die Vinogradov-Wooley Schranke für G(k) zurück.

        Allgemeine Formel:  G(k) ≤ k(3 log k + 4.2 log log k + C)
        Für k=1,2,3 werden explizite bessere Schranken genutzt.

        @return: Obere Schranke für G(k) nach Vinogradov-Wooley
        """
        if self.k == 1:
            return 3
        if self.k == 2:
            return 5
        if self.k in self.G_BOUNDS_IMPROVED:
            return self.G_BOUNDS_IMPROVED[self.k]

        # Allgemeine Formel (Vinogradov-Wooley Heuristik)
        k = self.k
        if k >= 2:
            bound = int(k * (3 * math.log(k) + 4.2 * math.log(math.log(k) + 1) + 10))
            return bound
        return k * 4

    def represent_as_k_prime_powers(
        self, n: int, max_terms: Optional[int] = None
    ) -> Optional[List[int]]:
        """
        Findet Darstellung n = p₁ᵏ + p₂ᵏ + ... + pₜᵏ mit t ≤ max_terms.

        @param n: Darzustellende Zahl
        @param max_terms: Maximale Terme (Standard: G(k)-Schranke)
        @return: Liste der Primzahlen [p₁, ..., pₜ] oder None
        """
        if max_terms is None:
            max_terms = self.g_bound_hua() or self.g_bound_vinogradov_wooley()

        # Relevante Primzahlpotenzen ≤ n
        relevant = [(i, pw) for i, pw in enumerate(self.prime_powers) if pw <= n]
        if not relevant:
            return None

        def search(remaining: int, terms_left: int, min_idx: int) -> Optional[List[int]]:
            if remaining == 0:
                return []
            if terms_left == 0:
                return None
            for i in range(min_idx, len(relevant)):
                orig_i, pw = relevant[i]
                if pw > remaining:
                    break
                result = search(remaining - pw, terms_left - 1, i)
                if result is not None:
                    return [self.primes[orig_i]] + result
            return None

        return search(n, max_terms, 0)

    def vinogradov_three_primes_check(self, n: int) -> dict:
        """
        Prüft Vinogradov's Drei-Primzahl-Theorem für k=1 (n = p₁+p₂+p₃).

        **Vinogradov (1937):**
            Jede hinreichend große ODD Zahl ist Summe von 3 Primzahlen.
            (Beweis für n > e^{e^{11.503}} ≈ 10^{43000})

        **Helfgott (2013):**
            Für alle ungeraden n > 5 gilt die Darstellung (Vollbeweis).

        @param n: Zu prüfende ungerade Zahl
        @return: Dict mit Darstellungsergebnis
        """
        if self.k != 1:
            return {'error': 'Nur für k=1 anwendbar'}

        if n % 2 == 0:
            return {
                'n': n,
                'is_odd': False,
                'note': 'Vinogradov gilt nur für ungerade Zahlen'
            }

        primes = self.primes
        prime_set = set(primes)

        # Suche p₁ + p₂ + p₃ = n
        representation = None
        for i, p1 in enumerate(primes):
            if p1 >= n:
                break
            for p2 in primes[i:]:
                p3 = n - p1 - p2
                if p3 < p2:
                    break
                if p3 in prime_set:
                    representation = [p1, p2, p3]
                    break
            if representation:
                break

        return {
            'n': n,
            'is_odd': True,
            'representation': representation,
            'found': representation is not None,
            'sum_check': sum(representation) == n if representation else None
        }

    def goldbach_conjecture_check(self, n: int) -> dict:
        """
        Prüft die Goldbach-Vermutung für eine gerade Zahl n.

        **Goldbach-Conjecture (1742, UNBEWIESEN):**
            Jede gerade Zahl ≥ 4 ist Summe von zwei Primzahlen.
            n = p₁ + p₂

        Status: OFFEN / Conjecture (numerisch bis 4·10^{18} verifiziert).

        @param n: Gerade Zahl ≥ 4
        @return: Dict mit Darstellungsergebnis
        """
        if n % 2 != 0 or n < 4:
            return {'n': n, 'error': 'Nur für gerade Zahlen ≥ 4'}

        prime_set = set(self.primes)
        representation = None

        for p in self.primes:
            if p > n // 2:
                break
            q = n - p
            if q in prime_set:
                representation = [p, q]
                break

        return {
            'n': n,
            'representation': representation,
            'found': representation is not None,
            'goldbach_status': 'CONJECTURE (unbewiesen)',
            'verified_up_to': '4 * 10^18 (numerisch)'
        }

    def waring_goldbach_statistics(self, limit: int) -> dict:
        """
        Statistiken über Waring-Goldbach-Darstellungen bis limit.

        @param limit: Obere Grenze für die Untersuchung
        @return: Dict mit Statistiken
        """
        bound = self.g_bound_hua() or self.g_bound_vinogradov_wooley()
        results = {
            'k': self.k,
            'g_bound': bound,
            'representations': {},
            'min_terms_needed': {}
        }

        for n in range(2, limit + 1):
            for t in range(1, bound + 1):
                rep = self.represent_as_k_prime_powers(n, max_terms=t)
                if rep is not None:
                    results['min_terms_needed'][n] = t
                    break

        if results['min_terms_needed']:
            counts = list(results['min_terms_needed'].values())
            results['avg_terms'] = sum(counts) / len(counts)
            results['max_terms_used'] = max(counts)
            results['distribution'] = {
                t: counts.count(t) for t in range(1, bound + 2)
                if counts.count(t) > 0
            }

        return results

    def explicit_vinogradov_wooley_bound(self) -> str:
        """
        Gibt die explizite Vinogradov-Wooley Formel als String zurück.

        @return: Formel-String mit Erklärung
        """
        k = self.k
        bound = self.g_bound_vinogradov_wooley()
        return (
            f"Waring-Goldbach G({k}) nach Vinogradov-Wooley:\n"
            f"  G({k}) ≤ k·(3·log(k) + 4.2·log(log(k)) + C)\n"
            f"  Für k={k}: G({k}) ≤ {bound}\n"
            f"  Methode: Efficient Congruencing (Wooley 1992-2019)"
        )

    def hua_theorem_statement(self) -> str:
        """
        Gibt Hua's Theorem als formalen String zurück.

        @return: Theorem-Formulierung
        """
        g = self.g_bound_hua()
        k = self.k
        if g is None:
            return f"Hua's Schranke für k={k} nicht direkt tabelliert."

        return (
            f"Hua's Theorem (1938) für k={k}:\n"
            f"  Jede hinreichend große natürliche Zahl n kann als\n"
            f"  n = p₁^{k} + p₂^{k} + ... + p_{g}^{k}\n"
            f"  mit Primzahlen p₁,...,p_{g} dargestellt werden.\n"
            f"  G({k}) ≤ {g}"
        )
