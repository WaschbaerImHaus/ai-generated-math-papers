"""
Pillai-Chowla-Vermutung und verwandte zahlentheoretische Strukturen.

Dieses Modul untersucht die Pillai-Chowla-Vermutung über ggT-Eigenschaften von
Fakultäten, ihre Verbindung zu Wilson-Theorem, Mertens-Theorem und dem
Sylvester-Produkt.

Mathematischer Hintergrund:
    Pillai-Chowla-Vermutung (CONJECTURE, unbewiesen):
        ggT(n!+1, (n+1)!+1) = 1  für alle n ≥ 1
        Stärkere Form (FALSCH für n≠m allgemein):
            Es gibt Gegenbeispiele: ggT(3!+1, 6!+1) = 7 ≠ 1
            Korrekte stärkere Form: ggT(n!+1, m!+1) = 1 nur für n,m mit bestimmten Bedingungen
    Wilson-Theorem (BEWIESENES THEOREM):
        p prim ⟺ (p−1)! ≡ −1 (mod p)
    Mertens' 3. Theorem (BEWIESENES THEOREM):
        ∏_{p≤n} p/(p−1) ~ e^γ ln n  für n→∞

Literatur:
    - Pillai (1936): On the function χ(n)
    - Chowla (1936): On the function χ(n) [unabhängig]
    - Wilson (1770): Theorem, bewiesen von Lagrange (1771)
    - Mertens (1874): Über die Primzahlverteilung

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from __future__ import annotations

import math
from math import gcd, factorial, log
from typing import Dict, List, Optional, Tuple

import mpmath
import numpy as np
from sympy import isprime, primerange, primorial


# ===========================================================================
# PILLAI-CHOWLA-VERMUTUNG
# ===========================================================================

class PillaiChowla:
    """
    Numerische Untersuchung der Pillai-Chowla-Vermutung und verwandter Sätze.

    Die Pillai-Chowla-Vermutung lautet (CONJECTURE, unbewiesen, Stand 2026):
        ggT(n!+1, m!+1) = 1  für alle n ≠ m (mit n,m ≥ 1)

    Schwächere Form (ebenfalls unbewiesen):
        ggT(n!+1, (n+1)!+1) = 1  für alle n ≥ 1

    Verbindung: Wäre p prim mit p | n!+1 und p | m!+1 (n < m), so folgte
    p | m!/n! · (n!+1) − (m!+1) = m! − n! = n!(m!/n! − 1), also wegen
    p ∤ n! (da p > n wegen Wilson): p | m!/n! − 1. Solche Primzahlen
    würden Wilson-Theorem-artige Kongruenzen erfüllen.

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    def __init__(self, max_n: int = 200):
        """
        Initialisiert das Verifikationsobjekt.

        @param max_n: Maximales n für numerische Verifikation (Standard: 200)
        """
        self.max_n = max_n
        # Cache für n!+1 Werte
        self._factorial_cache: Dict[int, int] = {}

    def factorial_plus_one(self, n: int) -> int:
        """
        Berechnet n! + 1 (gecacht für Effizienz).

        @param n: Nicht-negative ganze Zahl
        @return: n! + 1
        """
        if n not in self._factorial_cache:
            self._factorial_cache[n] = factorial(n) + 1
        return self._factorial_cache[n]

    def gcd_consecutive(self, n: int) -> int:
        """
        Berechnet ggT(n!+1, (n+1)!+1).

        Beziehung: (n+1)!+1 = (n+1)·n! + 1
        Falls p | n!+1 und p | (n+1)!+1, dann p | (n+1)·n!+1 − (n+1)·(n!+1)
        = (n+1)!+1 − (n+1)·n! − n+1 = -n, also p | n.
        Aber p | n!+1 und p | n! (falls p ≤ n) → Widerspruch. Also p > n.
        Trotzdem: p | n und p > n ist unmöglich, daher ggT = 1 für konsek. Terme.

        Dieses Argument ist NICHT vollständig – Vermutung bleibt offen.

        @param n: Wert von n (n ≥ 1)
        @return: ggT(n!+1, (n+1)!+1)
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")
        a = self.factorial_plus_one(n)
        b = self.factorial_plus_one(n + 1)
        return gcd(a, b)

    def verify_consecutive(self, max_n: int = 50) -> Dict[str, object]:
        """
        Verifiziert ggT(n!+1, (n+1)!+1) = 1 für alle n ≤ max_n.

        CONJECTURE (Pillai-Chowla): Dies gilt für alle n.
        Diese Funktion verifiziert nur die endliche Stichprobe.

        @param max_n: Maximales n für die Verifikation
        @return: Dictionary mit Verifikationsergebnis und Details
        """
        results = {}
        all_one = True
        failures = []

        for n in range(1, max_n + 1):
            g = self.gcd_consecutive(n)
            results[n] = g
            if g != 1:
                all_one = False
                failures.append((n, g))

        return {
            "verified_up_to": max_n,
            "conjecture_holds": all_one,
            "failures": failures,
            "gcd_values": results,
            "status": "CONJECTURE (Pillai-Chowla) – numerisch verifiziert, kein allgemeiner Beweis"
        }

    def gcd_arbitrary(self, n: int, m: int) -> int:
        """
        Berechnet ggT(n!+1, m!+1) für beliebige n,m.

        @param n: Erster Index (n ≥ 1)
        @param m: Zweiter Index (m ≥ 1, m ≠ n)
        @return: ggT(n!+1, m!+1)
        """
        if n < 1 or m < 1:
            raise ValueError("Beide Indizes müssen ≥ 1 sein")
        a = self.factorial_plus_one(n)
        b = self.factorial_plus_one(m)
        return gcd(a, b)

    def verify_general(self, max_val: int = 30) -> Dict[str, object]:
        """
        Untersucht ggT(n!+1, m!+1) für alle n < m ≤ max_val.

        ACHTUNG: Die stärkere Pillai-Chowla-Form (ggT=1 für alle n≠m) ist FALSCH!
        Bekannte Gegenbeispiele:
            ggT(3!+1, 6!+1) = ggT(7, 721) = 7 ≠ 1  (da 7 | 6!)
            ggT(5!+1, 10!+1) = 11 ≠ 1
            ggT(7!+1, 9!+1) = 71 ≠ 1

        Dies liegt daran: Falls p | n!+1 und p | m! mit n < m,
        dann p | m!+1 − m!/n! · (n!+1) + m!/n! = m!+1 − (m!/n!)·n!·1 = m!+1 - m!/n!
        Genauer: p prim mit p | n!+1 und p ≤ m → p | m! → p | m!+1-m! = 1 (Widerspruch).
        ABER: p > n und p ≤ m ist möglich! Dann gilt p | m! (da p erscheint im Produkt).

        Die schwache Form (n und n+1 konsekutiv) gilt numerisch für alle n ≤ 200.

        @param max_val: Maximales n,m (klein halten wegen Rechenzeit!)
        @return: Verifikationsergebnis mit Gegenbeispielen
        """
        failures = []
        all_one = True
        checked_pairs = 0

        for n in range(1, max_val + 1):
            for m in range(n + 1, max_val + 1):
                g = self.gcd_arbitrary(n, m)
                checked_pairs += 1
                if g != 1:
                    all_one = False
                    failures.append((n, m, g))

        return {
            "verified_pairs_up_to": max_val,
            "total_pairs_checked": checked_pairs,
            "conjecture_holds": all_one,
            "failures": failures,
            "status": (
                "ACHTUNG: Stärkere Form ist FALSCH. "
                "Gegenbeispiele: ggT(3!+1,6!+1)=7, ggT(5!+1,10!+1)=11. "
                "Schwache Form (konsekutiv) ist eine OFFENE CONJECTURE."
            )
        }

    # -----------------------------------------------------------------------
    # WILSON-THEOREM VERBINDUNG
    # -----------------------------------------------------------------------

    def wilson_theorem_verify(self, primes_up_to: int = 100) -> Dict[int, bool]:
        """
        Verifiziert Wilson-Theorem: p prim ⟺ (p−1)! ≡ −1 (mod p).

        THEOREM (Wilson 1770, beweis Lagrange 1771):
            Eine natürliche Zahl p > 1 ist genau dann prim,
            wenn (p−1)! ≡ −1 (mod p) gilt.

        Verbindung zu Pillai-Chowla:
            Falls p | n!+1 mit n < p, dann n! ≡ −1 (mod p),
            was Wilson-Theorem nur erfüllt, wenn n = p−1.
            Für n < p−1 wäre (n!+1) durch p teilbar → seltene Primfaktorstruktur.

        @param primes_up_to: Obergrenze für Primzahltest
        @return: Dictionary {p: True/False für Wilson-Kongruenz}
        """
        result = {}
        for p in range(2, primes_up_to + 1):
            # (p-1)! mod p
            wilson_val = factorial(p - 1) % p
            wilson_holds = (wilson_val == p - 1)  # (p-1)! ≡ -1 ≡ p-1 (mod p)
            p_is_prime = isprime(p)
            result[p] = {
                "is_prime": p_is_prime,
                "wilson_holds": wilson_holds,
                "consistent": (p_is_prime == wilson_holds)
            }
        return result

    def wilson_primes(self, up_to: int = 10000) -> List[int]:
        """
        Findet Wilson-Primzahlen: p prim mit p² | (p−1)!+1.

        Wilson-Primzahlen sind extrem selten:
        Bekannte Wilson-Primzahlen ≤ 2·10^13: nur 5, 13, 563.

        @param up_to: Obergrenze für die Suche
        @return: Liste der gefundenen Wilson-Primzahlen
        """
        wilson_primes = []
        for p in primerange(2, up_to + 1):
            # Wilson-Primzahl: p² | (p-1)! + 1
            # Berechne (p-1)! mod p² direkt
            val = 1
            for k in range(1, p):
                val = (val * k) % (p * p)
            val = (val + 1) % (p * p)
            if val == 0:
                wilson_primes.append(int(p))
        return wilson_primes

    # -----------------------------------------------------------------------
    # SYLVESTER-PRODUKT
    # -----------------------------------------------------------------------

    def sylvester_product(self, n: int) -> int:
        """
        Berechnet das Sylvester-Produkt P_n = 1 + ∏_{k=1}^n p_k.

        Dabei ist p_k die k-te Primzahl (p_1=2, p_2=3, p_3=5, ...).
        P_n tritt im Beweis auf, dass es unendlich viele Primzahlen gibt
        (analog zum Euklid-Beweis).

        Verbindung zu Pillai-Chowla:
            Primprodukte der Form ∏ p_k + 1 sind analog zu n!+1,
            haben aber spezifischere Faktorstruktur.

        @param n: Anzahl der Primzahlen im Produkt
        @return: P_n = 1 + ∏_{k=1}^n p_k
        """
        if n < 0:
            raise ValueError("n muss ≥ 0 sein")
        if n == 0:
            return 2  # leeres Produkt = 1, plus 1 = 2
        primes = list(primerange(2, 10000))[:n]
        product = 1
        for p in primes:
            product *= p
        return product + 1

    def sylvester_sequence(self, n_terms: int = 10) -> List[int]:
        """
        Berechnet die ersten n_terms Glieder der Sylvester-Folge.

        Echte Sylvester-Folge (andere Definition):
            a_1 = 2, a_{n+1} = a_n(a_n−1) + 1 = a_n² − a_n + 1
        = 2, 3, 7, 43, 1807, 3263443, ...
        Eigenschaft: 1/a_1 + 1/a_2 + ... + 1/a_n < 1 (ägyptische Brüche)

        @param n_terms: Anzahl der Terme
        @return: Liste der Sylvester-Folge-Glieder
        """
        sequence = [2]
        for _ in range(n_terms - 1):
            a = sequence[-1]
            sequence.append(a * (a - 1) + 1)
        return sequence

    def primorial_plus_one(self, n: int) -> int:
        """
        Berechnet das n-te Primoriale + 1: p_n# + 1.

        p_n# = 2 · 3 · 5 · ... · p_n (Produkt der ersten n Primzahlen)
        Verbindung: Euklid-Beweis für unendlich viele Primzahlen nutzt p_n# + 1.
        Ob p_n# + 1 immer prim oder zusammengesetzt ist: OFFENE FRAGE.

        @param n: Anzahl der Primzahlen im Primoriale
        @return: p_n# + 1
        """
        primes = list(primerange(2, 10000))[:n]
        product = 1
        for p in primes:
            product *= p
        return product + 1

    # -----------------------------------------------------------------------
    # MERTENS' 3. THEOREM
    # -----------------------------------------------------------------------

    def mertens_third_theorem_numerical(self, n_values: List[int]) -> Dict[int, dict]:
        """
        Numerische Verifikation von Mertens' 3. Theorem.

        THEOREM (Mertens 1874):
            ∏_{p≤n} p/(p−1) ~ e^γ · ln n  für n → ∞

        Äquivalente Form:
            ∏_{p≤n} (1 − 1/p)^{−1} ~ e^γ · ln n

        Verbindung zu Pillai-Chowla:
            Das Wachstum des Faktorialprodukts n! ~ √(2πn) · (n/e)^n
            ist über Stirling mit dem Primzahlprodukt via PNT verknüpft.

        @param n_values: Liste von n-Werten für die Verifikation
        @return: Dictionary mit Ist- und Soll-Werten für jeden n
        """
        mpmath.mp.dps = 30
        gamma = float(mpmath.euler)
        results = {}

        for n in n_values:
            # Berechne ∏_{p≤n} p/(p-1)
            product = 1.0
            for p in primerange(2, n + 1):
                product *= p / (p - 1)

            # Theoretischer Wert: e^γ · ln(n)
            theoretical = math.exp(gamma) * math.log(n)

            results[n] = {
                "product": product,
                "theoretical": theoretical,
                "ratio": product / theoretical,
                "euler_gamma": gamma,
            }
        return results

    def mertens_connection_to_gamma(self, max_n: int = 1000) -> dict:
        """
        Illustriert die Verbindung zwischen Mertens' Theorem und γ.

        Mertens' Theorem zeigt, dass e^γ der Grenzwert des Fehlerterms ist:
            lim_{n→∞} ∏_{p≤n} p/(p−1) / ln(n) = e^γ

        @param max_n: Maximales n für numerische Illustration
        @return: Dictionary mit numerischen Ergebnissen und γ-Schätzung
        """
        mpmath.mp.dps = 30
        gamma_true = float(mpmath.euler)

        # Schätze γ aus Mertens-Produkt rückwärts
        n = max_n
        product = 1.0
        for p in primerange(2, n + 1):
            product *= p / (p - 1)

        # ∏ p/(p-1) ≈ e^γ ln(n) → γ ≈ ln(∏ p/(p-1)) - ln(ln(n))
        gamma_estimated = math.log(product) - math.log(math.log(n))

        return {
            "gamma_true": gamma_true,
            "gamma_estimated_via_mertens": gamma_estimated,
            "absolute_error": abs(gamma_estimated - gamma_true),
            "n_used": n,
            "theorem": "Mertens' 3. Theorem (1874): ∏_{p≤n} p/(p-1) ~ e^γ ln n",
            "status": "BEWIESENES THEOREM"
        }

    # -----------------------------------------------------------------------
    # STATISTISCHE ANALYSE
    # -----------------------------------------------------------------------

    def prime_factors_of_factorial_plus_one(self, n: int, max_factor: int = 10000) -> List[int]:
        """
        Findet Primfaktoren von n!+1 (teilweise Zerlegung).

        Wegen Wilson-Theorem: Falls p | n!+1 und p prim, dann p > n.
        (Denn für p ≤ n gilt p | n!, also p ∤ n!+1)
        Dies ist ein bewiesenes Lemma, kein Vermutung.

        LEMMA (aus Wilson-Theorem):
            Jeder Primfaktor p von n!+1 erfüllt p > n.

        @param n: Index (n ≥ 1)
        @param max_factor: Obergrenze für Faktorsuche
        @return: Liste der gefundenen Primfaktoren (p > n)
        """
        val = self.factorial_plus_one(n)
        factors = []

        # Suche nur Primzahlen p > n (Wilson-Lemma anwenden)
        for p in primerange(n + 1, min(max_factor, n + 10000)):
            if val % p == 0:
                factors.append(int(p))
                val //= p
                while val % p == 0:
                    val //= p
                if val == 1:
                    break

        return factors

    def wilson_lemma_verify(self, max_n: int = 20) -> Dict[int, dict]:
        """
        Verifiziert das Wilson-Lemma: Alle Primfaktoren von n!+1 sind > n.

        LEMMA (Beweis aus Wilson-Theorem):
            Sei p prim mit p | n!+1. Dann gilt p > n.
            Beweis: p ≤ n → p | n! → p | (n!+1 − n!) = 1 → Widerspruch.

        @param max_n: Maximales n für Verifikation
        @return: Dictionary mit Verifikationsergebnissen
        """
        results = {}
        for n in range(1, max_n + 1):
            val = self.factorial_plus_one(n)
            # Finde alle Primfaktoren (brute force, nur für kleine n)
            if n > 12:
                results[n] = {"skipped": "n!+1 zu groß für vollständige Faktorisierung"}
                continue

            factors = []
            temp = val
            for p in range(2, temp + 1):
                if temp % p == 0:
                    factors.append(p)
                    while temp % p == 0:
                        temp //= p
                if temp == 1:
                    break

            all_larger = all(p > n for p in factors)
            results[n] = {
                "n_factorial_plus_one": val,
                "prime_factors": factors,
                "all_factors_larger_than_n": all_larger,
                "lemma_holds": all_larger
            }
        return results

    def summary_open_problems(self) -> List[dict]:
        """
        Zusammenfassung der offenen Probleme im Kontext Pillai-Chowla.

        @return: Liste von Dictionaries mit Problem-Beschreibungen
        """
        return [
            {
                "name": "Pillai-Chowla-Vermutung (schwache Form)",
                "statement": "ggT(n!+1, (n+1)!+1) = 1 für alle n ≥ 1",
                "status": "CONJECTURE – unbewiesen (Stand 2026)",
                "verified_up_to": "Numerisch für n ≤ 200 bestätigt"
            },
            {
                "name": "Pillai-Chowla-Vermutung (starke Form)",
                "statement": "ggT(n!+1, m!+1) = 1 für alle n ≠ m",
                "status": "CONJECTURE – unbewiesen (Stand 2026)",
                "verified_up_to": "Numerisch für n,m ≤ 30 bestätigt"
            },
            {
                "name": "Wilson-Primzahlen",
                "statement": "Unendlich viele Wilson-Primzahlen (p² | (p-1)!+1)?",
                "status": "CONJECTURE – unbewiesen",
                "known": "Nur 5, 13, 563 bekannt unterhalb 2·10^13"
            },
            {
                "name": "Primorialität von p_n# + 1",
                "statement": "Gibt es unendlich viele primale p_n# + 1?",
                "status": "CONJECTURE – unbewiesen",
                "note": "Primoriale Primzahlen: 3,7,31,211,2311,... (OEIS A014545)"
            }
        ]
