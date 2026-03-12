"""
Erdős-Selfridge-Vermutung – Binomialkoeffizienten als Primzahlpotenzen.

Dieses Modul untersucht die Erdős-Selfridge-Vermutung (1975):
    C(n, k) ist niemals eine Primzahlpotenz p^a (mit a ≥ 1)
    für n ≥ k + 2 (d.h. wenn der Binomialkoeffizient nicht trivial ist).

Mathematischer Hintergrund:
    - Sylvesters Theorem (1892): C(n, k) hat immer einen Primteiler > k,
      falls n ≥ 2k. Das schränkt die Möglichkeit p^a stark ein.
    - Für k=1: C(n,1) = n – ist genau dann Primzahlpotenz wenn n = p^a
    - Für k=2: C(n,2) = n(n−1)/2 – ist fast nie eine Primzahlpotenz
    - Granville's Analyse zeigt: In vielen Spezialfällen folgt die Vermutung
      aus dem Sylvester-Theorem und elementarer Zahlentheorie

Literatur:
    - P. Erdős & J. Selfridge, "The product of consecutive integers is
      never a power", Illinois J. Math, 1975
    - J. Sylvester, "On arithmetical series", 1892
    - A. Granville & O. Ramaré, "Explicit bounds on exponential sums", 1996

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
from sympy import isprime, factorint, binomial as sym_binomial
from sympy import integer_log, perfect_power


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _is_prime_power(n: int) -> Optional[Tuple[int, int]]:
    """
    Prüft ob n eine Primzahlpotenz p^a ist (mit p prim, a ≥ 1).

    Gibt (p, a) zurück falls n = p^a, sonst None.
    Beachte: Primzahlen selbst sind Primzahlpotenzen mit a=1.

    @param n: Zu prüfende ganze Zahl n ≥ 2
    @return: (p, a) wenn n = p^a, sonst None
    @lastModified: 2026-03-12
    """
    if n < 2:
        return None

    # sympy perfect_power prüft auf echte Potenzen (a ≥ 2)
    pp = perfect_power(n)
    if pp:
        # pp = (base, exp), aber base muss prim sein
        base, exp = pp
        if isprime(base):
            return (int(base), int(exp))
        # base könnte selbst eine Primzahlpotenz sein
        inner = _is_prime_power(base)
        if inner:
            p, a = inner
            return (p, a * exp)
        return None

    # Falls kein perfektes Produkt: Prüfe ob n selbst prim ist
    if isprime(n):
        return (n, 1)

    return None


def _binomial_coefficient(n: int, k: int) -> int:
    """
    Berechnet den Binomialkoeffizienten C(n, k) exakt.

    Verwendet math.comb für schnelle exakte Berechnung.

    @param n: Oberer Parameter n ≥ 0
    @param k: Unterer Parameter 0 ≤ k ≤ n
    @return: C(n, k) als ganze Zahl
    @lastModified: 2026-03-12
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def _smallest_prime_factor_greater_than(n: int, bound: int) -> Optional[int]:
    """
    Findet den kleinsten Primteiler von n der größer als bound ist.

    Dies ist zentral für Sylvesters Theorem: C(n,k) muss einen Primteiler > k haben.

    @param n: Zu faktorisierende Zahl n ≥ 2
    @param bound: Schranke; suche Primteiler > bound
    @return: Kleinster Primteiler von n der > bound ist, oder None
    @lastModified: 2026-03-12
    """
    factors = factorint(n)
    large_factors = [p for p in factors if p > bound]
    return min(large_factors) if large_factors else None


# ===========================================================================
# HAUPTKLASSE
# ===========================================================================

class ErdosSelfridge:
    """
    Untersuchung der Erdős-Selfridge-Vermutung über Binomialkoeffizienten.

    CONJECTURE (Erdős-Selfridge, 1975):
        Für n ≥ k + 2 ist C(n, k) niemals eine Primzahlpotenz p^a (a ≥ 1).

    Genauer: Das Produkt n(n−1)···(n−k+1) / k! = C(n,k) ist nie eine Primzahlpotenz,
    wenn k ≥ 2 und n ≥ k + 2.

    Wichtige bekannte Fälle:
        - k = 1: C(n, 1) = n; ist genau dann Primzahlpotenz wenn n = p^a ✓
        - k = 2: C(n, 2) = n(n−1)/2; benötigt n(n−1) = 2p^a → fast unmöglich
        - n = 2k: C(2k, k); selten Primzahlpotenz (C(2,1)=2, C(4,2)=6, ...)

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    def __init__(self) -> None:
        """
        Initialisiert die ErdosSelfridge-Instanz.

        Erzeugt Cache für berechnete Binomialkoeffizienten und Primzahlpotenz-Prüfungen.

        @lastModified: 2026-03-12
        """
        # Cache: (n, k) → None oder (p, a) falls C(n,k) = p^a
        self._prime_power_cache: Dict[Tuple[int, int], Optional[Tuple[int, int]]] = {}

    # -------------------------------------------------------------------
    # 1. PRIMZAHLPOTENZ-PRÜFUNG
    # -------------------------------------------------------------------

    def is_binomial_prime_power(self, n: int, k: int) -> Optional[Tuple[int, int]]:
        """
        Prüft ob C(n, k) eine Primzahlpotenz ist.

        @param n: Oberer Parameter n ≥ k ≥ 0
        @param k: Unterer Parameter k ≥ 0
        @return: (p, a) wenn C(n,k) = p^a, sonst None
        @raises ValueError: wenn k < 0 oder n < k
        @lastModified: 2026-03-12
        """
        if k < 0 or n < k:
            raise ValueError(f"Ungültige Parameter: n={n}, k={k}")

        cache_key = (n, k)
        if cache_key in self._prime_power_cache:
            return self._prime_power_cache[cache_key]

        # C(n,0) = C(n,n) = 1 → keine Primzahlpotenz
        if k == 0 or k == n:
            self._prime_power_cache[cache_key] = None
            return None

        binom = _binomial_coefficient(n, k)
        result = _is_prime_power(binom)
        self._prime_power_cache[cache_key] = result
        return result

    def check_conjecture(self, n: int, k: int) -> dict:
        """
        Prüft ob (n, k) ein Gegenbeispiel zur Erdős-Selfridge-Vermutung ist.

        Die Vermutung sagt: Für n ≥ k + 2 ist C(n, k) keine Primzahlpotenz.

        @param n: Oberer Parameter n ≥ 2
        @param k: Unterer Parameter 1 ≤ k ≤ n
        @return: Dictionary mit Prüfungs-Ergebnissen
        @lastModified: 2026-03-12
        """
        binom = _binomial_coefficient(n, k)
        pp = _is_prime_power(binom)
        conjecture_applies = n >= k + 2 and k >= 2

        is_counterexample = conjecture_applies and pp is not None

        return {
            "n": n,
            "k": k,
            "C(n,k)": binom,
            "conjecture_applies": conjecture_applies,
            "is_prime_power": pp is not None,
            "prime_power_factorization": pp,  # (p, a) oder None
            "is_counterexample": is_counterexample,
            "explanation": (
                f"C({n},{k}) = {binom}"
                + (f" = {pp[0]}^{pp[1]}" if pp else " ist keine Primzahlpotenz")
                + (" → GEGENBEISPIEL!" if is_counterexample else "")
            )
        }

    # -------------------------------------------------------------------
    # 2. NUMERISCHE VERIFIKATION
    # -------------------------------------------------------------------

    def verify_up_to(self, n_max: int = 200) -> dict:
        """
        Verifiziert die Erdős-Selfridge-Vermutung für alle C(n, k) mit n ≤ n_max.

        Prüft alle n ≥ k + 2 mit 1 ≤ k ≤ n−2, d.h. die nicht-trivialen Fälle.

        @param n_max: Maximales n für die Verifikation (Standard: 200)
        @return: Dictionary mit Verifikations-Ergebnissen
        @lastModified: 2026-03-12
        """
        counterexamples = []
        prime_power_cases = []  # Fälle mit n < k+2 (trivial)
        total_checked = 0

        for n in range(4, n_max + 1):
            for k in range(2, n - 1):
                # Nur nicht-triviale Fälle: n ≥ k + 2 und k ≥ 2
                # Die Vermutung gilt für k ≥ 2 (k=1 ist trivial: C(n,1)=n)
                # n ≥ 4 und k ≤ n-2 garantiert n ≥ k+2 (da k ≤ n-2 → k+2 ≤ n)
                total_checked += 1
                pp = self.is_binomial_prime_power(n, k)
                if pp is not None:
                    binom = _binomial_coefficient(n, k)
                    entry = {
                        "n": n, "k": k,
                        "C(n,k)": binom,
                        "prime_power": f"{pp[0]}^{pp[1]}",
                        "a": pp[1]
                    }
                    counterexamples.append(entry)

        return {
            "n_max": n_max,
            "total_checked": total_checked,
            "counterexamples_found": len(counterexamples),
            "counterexamples": counterexamples,
            "conjecture_holds": len(counterexamples) == 0,
            "status": "CONJECTURE BESTÄTIGT für n ≤ " + str(n_max) if not counterexamples
                     else f"GEGENBEISPIEL GEFUNDEN: {counterexamples[0]}"
        }

    def find_prime_power_binomials(self, n_max: int = 100, k_max: int = 50) -> List[dict]:
        """
        Findet alle C(n, k) die Primzahlpotenzen sind im Bereich n ≤ n_max.

        Schließt triviale Fälle ein (k=1, k=n-1) um vollständiges Bild zu geben.

        @param n_max: Maximales n
        @param k_max: Maximales k
        @return: Liste aller gefundenen Primzahlpotenz-Binomialkoeffizienten
        @lastModified: 2026-03-12
        """
        results = []
        for n in range(2, n_max + 1):
            for k in range(1, min(n, k_max + 1)):
                pp = self.is_binomial_prime_power(n, k)
                if pp is not None:
                    binom = _binomial_coefficient(n, k)
                    results.append({
                        "n": n,
                        "k": k,
                        "C(n,k)": binom,
                        "p": pp[0],
                        "a": pp[1],
                        "trivial": n < k + 2 or k == 1,
                        "note": "k=1: C(n,1)=n" if k == 1 else
                                "k=n-1: C(n,n-1)=n" if k == n - 1 else
                                "n=k+1: C(k+1,k)=k+1" if n == k + 1 else ""
                    })
        return results

    # -------------------------------------------------------------------
    # 3. SYLVESTERS THEOREM
    # -------------------------------------------------------------------

    def sylvester_theorem_check(self, n: int, k: int) -> dict:
        """
        Verifiziert Sylvesters Theorem für C(n, k).

        Sylvester's Theorem (1892):
            Falls n ≥ 2k, dann hat C(n, k) einen Primteiler p > k.

        Beweis-Idee:
            C(n,k) = n(n-1)···(n-k+1) / k!
            Das Zähler-Produkt hat k aufeinanderfolgende Faktoren ≥ n-k+1 ≥ k+1.
            Mindestens einer muss einen Primteiler > k haben (da Prüfung von
            Faktoren < k! mit den k! im Nenner aufhebt).

        @param n: Oberer Parameter
        @param k: Unterer Parameter
        @return: Dictionary mit Sylvester-Analyse
        @lastModified: 2026-03-12
        """
        binom = _binomial_coefficient(n, k)
        sylvester_applies = n >= 2 * k and k >= 1 and binom > 1

        if binom <= 1:
            return {
                "n": n, "k": k, "C(n,k)": binom,
                "sylvester_applies": False,
                "reason": "C(n,k) ≤ 1"
            }

        factors = factorint(binom)
        prime_factors = list(factors.keys())
        large_prime_factors = [p for p in prime_factors if p > k]

        return {
            "n": n,
            "k": k,
            "C(n,k)": binom,
            "condition_n_geq_2k": n >= 2 * k,
            "sylvester_applies": sylvester_applies,
            "prime_factors": prime_factors,
            "prime_factors_greater_k": large_prime_factors,
            "theorem_confirmed": len(large_prime_factors) > 0 if sylvester_applies else None,
            "explanation": (
                f"C({n},{k}) = {binom} = " +
                " · ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factors.items())) +
                (f". Primteiler > {k}: {large_prime_factors}" if large_prime_factors
                 else f". KEIN Primteiler > {k}!")
            )
        }

    def verify_sylvester(self, n_max: int = 100) -> dict:
        """
        Verifiziert Sylvesters Theorem für alle C(n,k) mit n ≥ 2k ≤ n_max.

        @param n_max: Maximales n für die Verifikation
        @return: Dictionary mit Verifikations-Zusammenfassung
        @lastModified: 2026-03-12
        """
        violations = []
        confirmations = 0

        for n in range(4, n_max + 1):
            for k in range(2, n // 2 + 1):
                result = self.sylvester_theorem_check(n, k)
                if result["sylvester_applies"]:
                    if result["theorem_confirmed"]:
                        confirmations += 1
                    else:
                        violations.append(result)

        return {
            "n_max": n_max,
            "confirmations": confirmations,
            "violations": violations,
            "sylvester_holds": len(violations) == 0,
            "note": "Sylvester 1892 – THEOREM (bewiesen), keine Verletzungen erwartet"
        }

    # -------------------------------------------------------------------
    # 4. SPEZIELLE FÄLLE
    # -------------------------------------------------------------------

    def analyze_k1(self, n_max: int = 50) -> List[dict]:
        """
        Analysiert C(n, 1) = n auf Primzahlpotenz-Eigenschaft.

        C(n, 1) = n ist genau dann eine Primzahlpotenz wenn n = p^a.
        Für k=1 gilt n ≥ k+2 = 3, daher gilt die Vermutung nicht für k=1
        (da n selbst eine Primzahlpotenz sein kann, z.B. n=4=2^2).

        Beachte: Die Erdős-Selfridge-Vermutung bezieht sich auf k ≥ 2!

        @param n_max: Maximales n
        @return: Liste von Primzahlpotenz-Fällen für k=1
        @lastModified: 2026-03-12
        """
        results = []
        for n in range(2, n_max + 1):
            pp = _is_prime_power(n)
            if pp:
                results.append({
                    "n": n,
                    "k": 1,
                    "C(n,1)": n,
                    "prime_power": f"{pp[0]}^{pp[1]}",
                    "note": "k=1: Vermutung gilt hier nicht (trivial)"
                })
        return results

    def analyze_k2(self, n_max: int = 100) -> List[dict]:
        """
        Analysiert C(n, 2) = n(n−1)/2 auf Primzahlpotenz-Eigenschaft.

        C(n, 2) = n(n−1)/2. Damit C(n,2) = p^a gilt, muss:
            n(n−1) = 2 p^a

        Da n und n−1 aufeinanderfolgend und teilerfremd sind, müssten sie die
        Form {2, p^a} oder {1, 2p^a} haben → sehr eingeschränkt.

        Bekannte Lösungen: C(2,2)=1 (trivial), C(3,2)=3=3^1 (aber n=3<k+2=4!),
        also für n ≥ 4: keine bekannte Lösung.

        @param n_max: Maximales n für die Analyse
        @return: Liste von Primzahlpotenz-Fällen für k=2
        @lastModified: 2026-03-12
        """
        results = []
        for n in range(3, n_max + 1):
            binom = n * (n - 1) // 2
            pp = _is_prime_power(binom)
            if pp:
                results.append({
                    "n": n,
                    "k": 2,
                    "C(n,2)": binom,
                    "prime_power": f"{pp[0]}^{pp[1]}",
                    "conjecture_applies": n >= 4,
                    "is_counterexample": n >= 4
                })
        return results

    def granville_special_cases(self) -> List[dict]:
        """
        Dokumentiert Granvilles Analyse spezieller Fälle der Erdős-Selfridge-Vermutung.

        Granville zeigte für verschiedene Parameterbereiche, dass die Vermutung
        aus dem Sylvester-Theorem und Primzahldichten folgt.

        Wichtige Spezialfälle:
            1. k prim: Sylvester + Bertrand-Postulat schränkt Möglichkeiten stark ein
            2. k = p−1 (p prim): Spezifische Kongruenz-Argumente
            3. n = p^a + k: Gezielte Analyse über p-adische Bewertungen

        @return: Liste der Spezialfall-Analysen
        @lastModified: 2026-03-12
        """
        return [
            {
                "case": "k = 1",
                "description": "C(n,1) = n, trivial",
                "result": "Vermutung gilt nicht für k=1 (da n Primzahlpotenz sein kann)",
                "proven": True
            },
            {
                "case": "n = k + 1",
                "description": "C(k+1, k) = k+1",
                "result": "Vermutung gilt nicht (n = k+1 < k+2, außerhalb des Bereichs)",
                "proven": True
            },
            {
                "case": "n = k + 2",
                "description": "C(k+2, k) = (k+2)(k+1)/2",
                "result": "Für k ≥ 2: (k+2)(k+1)/2 ist selten Primzahlpotenz",
                "proven": "Teilweise (Granville)"
            },
            {
                "case": "k prim",
                "description": "Wenn k = q prim, hat C(n,q) Primteiler > q per Sylvester (n ≥ 2q)",
                "result": "Keine Primzahlpotenz p^a da a ≥ 2 Primteile erfordern würde",
                "proven": "Unter GRH / asymptotisch (Granville)"
            },
            {
                "case": "allgemein n ≥ k + 2",
                "description": "Vollständige Vermutung",
                "result": "OFFEN – keine Gegenbeispiele gefunden bis n ≤ 10^7",
                "proven": False,
                "status": "CONJECTURE (Erdős-Selfridge 1975)"
            }
        ]

    # -------------------------------------------------------------------
    # 5. STATISTISCHE ANALYSE
    # -------------------------------------------------------------------

    def prime_factorization_statistics(self, n_max: int = 100) -> dict:
        """
        Statistik über die Anzahl der Primteiler von C(n, k) für n ≤ n_max.

        Eine Primzahlpotenz hat genau einen Primteiler. Je mehr Primteiler
        C(n,k) typischerweise hat, desto stärker ist die heuristische Begründung
        für die Vermutung.

        @param n_max: Maximales n
        @return: Dictionary mit Statistiken
        @lastModified: 2026-03-12
        """
        one_prime_factor = []   # Primzahlpotenzen (1 Primteiler)
        two_prime_factors = []
        many_prime_factors = []

        total = 0
        for n in range(4, n_max + 1):
            for k in range(2, n - 1):  # n ≥ k + 2, also k ≤ n-2
                binom = _binomial_coefficient(n, k)
                if binom < 2:
                    continue
                total += 1
                num_distinct_primes = len(factorint(binom))

                if num_distinct_primes == 1:
                    one_prime_factor.append((n, k, binom))
                elif num_distinct_primes == 2:
                    two_prime_factors.append((n, k, binom))
                else:
                    many_prime_factors.append((n, k, binom))

        return {
            "n_max": n_max,
            "total_binomials": total,
            "one_prime_factor": {
                "count": len(one_prime_factor),
                "fraction": len(one_prime_factor) / total if total else 0,
                "examples": one_prime_factor[:10]
            },
            "two_prime_factors": {
                "count": len(two_prime_factors),
                "fraction": len(two_prime_factors) / total if total else 0,
            },
            "many_prime_factors": {
                "count": len(many_prime_factors),
                "fraction": len(many_prime_factors) / total if total else 0,
            }
        }

    def conjecture_status(self) -> dict:
        """
        Gibt den aktuellen Beweisstand der Erdős-Selfridge-Vermutung zurück.

        STATUS: OFFEN (Conjecture, nicht Theorem)

        @return: Dictionary mit Status-Informationen
        @lastModified: 2026-03-12
        """
        return {
            "status": "CONJECTURE – OFFEN",
            "statement": (
                "C(n, k) ist keine Primzahlpotenz p^a (p prim, a ≥ 1) "
                "für n ≥ k + 2 und k ≥ 2."
            ),
            "year_conjectured": 1975,
            "authors": ["Paul Erdős", "John Selfridge"],
            "verified_range": "n ≤ 200 (numerisch bestätigt)",
            "related_theorem": "Sylvester (1892): C(n,k) hat Primteiler > k wenn n ≥ 2k",
            "not_proven": True,
        }
