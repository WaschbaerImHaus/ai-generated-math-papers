"""
Artin-Vermutung über primitive Wurzeln – Dichte und Artin-Konstante.

Dieses Modul implementiert die Artin-Vermutung über primitive Wurzeln:
    Für jedes a ∈ ℤ mit a ≠ −1 und a kein vollständiges Quadrat
    gilt: Es gibt unendlich viele Primzahlen p, für die a eine
    primitive Wurzel modulo p ist.

Mathematischer Hintergrund:
    - a ist primitive Wurzel mod p, wenn ord_p(a) = p − 1, d.h.
      a erzeugt die multiplikative Gruppe (ℤ/pℤ)*.
    - Artin-Konstante: C_A = ∏_{p prim} (1 − 1/(p(p−1))) ≈ 0.3739558136...
    - Hooley (1967): Unter GRH gilt Dichte = C_A * Anpassungsfaktor.
    - Heath-Brown (1986): Unbedingt für mindestens eine von drei Basen {2,3,5}.

Literatur:
    - Emil Artin, "Über die Zetafunktionen gewisser algebraischer Zahlkörper",
      1924 (ursprüngliche Vermutung)
    - Christopher Hooley, "On Artin's conjecture", J. reine angew. Math, 1967
    - D. R. Heath-Brown, "Artin's conjecture for primitive roots", 1986

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
from sympy import isprime, factorint, primitive_root, n_order
from sympy.ntheory import totient


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _sieve_primes(limit: int) -> List[int]:
    """
    Erzeugt alle Primzahlen bis limit via Sieb des Eratosthenes.

    @param limit: Obere Schranke (inklusiv)
    @return: Sortierte Liste aller Primzahlen ≤ limit
    @lastModified: 2026-03-12
    """
    if limit < 2:
        return []
    # Boolean-Array: sieve[i] = True bedeutet i ist prim
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            # Alle Vielfachen von i als nicht-prim markieren
            sieve[i * i::i] = False
    return list(np.where(sieve)[0])


def _multiplicative_order(a: int, n: int) -> Optional[int]:
    """
    Berechnet die multiplikative Ordnung ord_n(a) = min{k ≥ 1 : a^k ≡ 1 (mod n)}.

    Falls gcd(a, n) ≠ 1, existiert keine multiplikative Ordnung.

    @param a: Basis (muss teilerfremd zu n sein)
    @param n: Modulus n ≥ 2
    @return: ord_n(a) oder None falls gcd(a,n) ≠ 1
    @lastModified: 2026-03-12
    """
    if math.gcd(a % n, n) != 1:
        return None
    # Nutze sympy für korrekte Berechnung (schneller als naive Iteration)
    try:
        return int(n_order(a, n))
    except Exception:
        return None


# ===========================================================================
# HAUPTKLASSE
# ===========================================================================

class ArtinPrimitiveRoots:
    """
    Untersuchung der Artin-Vermutung über primitive Wurzeln.

    CONJECTURE (Artin, 1927):
        Sei a ∈ ℤ mit a ≠ 0, ±1 und a kein vollständiges Quadrat.
        Dann ist a primitive Wurzel für unendlich viele Primzahlen p.
        Die Dichte dieser Primzahlen ist (unter GRH):
            δ(a) = C_A = ∏_{p prim} (1 − 1/(p(p−1))) ≈ 0.37396

    Besondere Fälle:
        - Für a = −1: a ≡ p−1 (mod p) ist IMMER eine primitive Wurzel wenn p ≡ 3 (mod 4)
          → Aber Artin-Vermutung gilt für a = −1 nicht in der Standardform
        - Für a = Quadratzahl (4, 9, 16, ...): Vermutung gilt möglicherweise nicht

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Bekannter Näherungswert der Artin-Konstante
    # C_A = ∏_{p prim} (1 - 1/(p(p-1)))
    ARTIN_CONSTANT_APPROX: float = 0.3739558136192022880547280543074

    def __init__(self) -> None:
        """
        Initialisiert die ArtinPrimitiveRoots-Instanz.

        Erzeugt interne Caches für Primzahlen und berechnete Ordnungen.

        @lastModified: 2026-03-12
        """
        # Cache: prime_list[limit] → Liste der Primzahlen bis limit
        self._prime_cache: Dict[int, List[int]] = {}

    def _get_primes(self, limit: int) -> List[int]:
        """
        Gibt gecachte Primzahlen bis limit zurück.

        @param limit: Obere Schranke
        @return: Liste der Primzahlen ≤ limit
        @lastModified: 2026-03-12
        """
        if limit not in self._prime_cache:
            self._prime_cache[limit] = _sieve_primes(limit)
        return self._prime_cache[limit]

    # -------------------------------------------------------------------
    # 1. PRIMITIVE WURZEL PRÜFEN
    # -------------------------------------------------------------------

    def is_primitive_root(self, a: int, p: int) -> bool:
        """
        Prüft ob a eine primitive Wurzel modulo der Primzahl p ist.

        a ist primitive Wurzel mod p genau dann wenn:
            ord_p(a) = p − 1

        Äquivalent: a^((p-1)/q) ≢ 1 (mod p) für jeden Primteiler q von p−1.

        @param a: Basis (beliebige ganze Zahl)
        @param p: Primzahl p ≥ 3
        @return: True wenn a primitive Wurzel mod p ist
        @raises ValueError: wenn p keine Primzahl ist
        @lastModified: 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"p={p} ist keine Primzahl")

        a = a % p
        if a == 0:
            return False

        # Ordnung berechnen und mit p-1 vergleichen
        order = _multiplicative_order(a, p)
        return order == p - 1

    def primitive_root_check_detailed(self, a: int, p: int) -> dict:
        """
        Führt eine detaillierte Analyse durch ob a primitive Wurzel mod p ist.

        Gibt die Ordnung, die Faktorisierung von p−1 und den Vergleich aus.

        @param a: Basis
        @param p: Primzahl
        @return: Dictionary mit Analyse-Ergebnissen
        @lastModified: 2026-03-12
        """
        if not isprime(p):
            raise ValueError(f"p={p} ist keine Primzahl")

        a_red = a % p
        order = _multiplicative_order(a_red, p)
        phi_p = p - 1  # Für Primzahl p gilt φ(p) = p-1
        factors_phi = factorint(phi_p)

        return {
            "a": a,
            "a_mod_p": a_red,
            "p": p,
            "phi(p)": phi_p,
            "factors_of_phi(p)": factors_phi,
            "ord_p(a)": order,
            "is_primitive_root": order == phi_p,
            "explanation": (
                f"ord_{p}({a}) = {order}, φ({p}) = {phi_p}. "
                f"{'✓ Primitive Wurzel' if order == phi_p else '✗ Keine primitive Wurzel'}"
            )
        }

    # -------------------------------------------------------------------
    # 2. EMPIRISCHE DICHTE
    # -------------------------------------------------------------------

    def empirical_density(self, a: int, limit: int = 10000) -> dict:
        """
        Berechnet die empirische Dichte der Primzahlen p ≤ limit mit a als primitive Wurzel.

        Die Dichte ist der Anteil der Primzahlen p (mit p ∤ a) für die a
        eine primitive Wurzel mod p ist.

        Erwarteter Wert (Artin-Vermutung unter GRH): C_A ≈ 0.3739...
        (mit möglichen Anpassungsfaktoren für spezifische a)

        @param a: Basis (a ≠ 0, ±1, kein Quadrat für Standardvermutung)
        @param limit: Obere Grenze für Primzahlsuche
        @return: Dictionary mit Dichte-Informationen
        @lastModified: 2026-03-12
        """
        primes = self._get_primes(limit)
        # Nur Primzahlen die a nicht teilen
        relevant_primes = [p for p in primes if p > 2 and p != abs(a)]

        count_prim_root = 0
        prim_root_primes = []

        for p in relevant_primes:
            if self.is_primitive_root(a, p):
                count_prim_root += 1
                if len(prim_root_primes) < 20:
                    prim_root_primes.append(p)

        total = len(relevant_primes)
        density = count_prim_root / total if total > 0 else 0.0

        return {
            "a": a,
            "limit": limit,
            "total_primes_checked": total,
            "count_primitive_root": count_prim_root,
            "empirical_density": density,
            "artin_constant": self.ARTIN_CONSTANT_APPROX,
            "deviation": abs(density - self.ARTIN_CONSTANT_APPROX),
            "first_20_primes": prim_root_primes,
        }

    def count_primitive_root_primes(self, a: int, limit: int) -> int:
        """
        Zählt die Primzahlen p ≤ limit für die a eine primitive Wurzel ist.

        @param a: Basis
        @param limit: Obere Grenze
        @return: Anzahl der Primzahlen p ≤ limit mit a als primitive Wurzel
        @lastModified: 2026-03-12
        """
        primes = self._get_primes(limit)
        return sum(
            1 for p in primes
            if p > 2 and p != abs(a) and self.is_primitive_root(a, p)
        )

    # -------------------------------------------------------------------
    # 3. ARTIN-KONSTANTE
    # -------------------------------------------------------------------

    def artin_constant(self, num_primes: int = 1000) -> float:
        """
        Berechnet die Artin-Konstante C_A = ∏_{p prim} (1 − 1/(p(p−1))).

        Das unendliche Produkt wird über die ersten num_primes Primzahlen
        approximiert. Konvergiert langsam aber sicher gegen den wahren Wert
        C_A ≈ 0.3739558136...

        Herleitung: Jeder Primfaktor q von p−1 "sabotiert" mit
        Wahrscheinlichkeit 1/q die primitive-Wurzel-Eigenschaft.
        Das Produkt über alle q ergibt genau den Mobius-Inklusionsexklusionsfaktor.

        @param num_primes: Anzahl der Primzahlen für die Produktapproximation
        @return: Approximierter Wert der Artin-Konstante
        @lastModified: 2026-03-12
        """
        primes = _sieve_primes(num_primes * 15)[:num_primes]  # Genug Primzahlen
        product = 1.0
        for p in primes:
            # Faktor: (1 - 1/(p*(p-1)))
            product *= (1.0 - 1.0 / (p * (p - 1)))
        return product

    def artin_constant_convergence(self, steps: List[int]) -> List[dict]:
        """
        Zeigt die Konvergenz der Artin-Konstante bei zunehmend mehr Primzahlen.

        @param steps: Liste von Primzahl-Anzahlen für die Konvergenz-Messung
        @return: Liste von Dictionaries mit Konvergenz-Daten
        @lastModified: 2026-03-12
        """
        results = []
        for n in steps:
            approx = self.artin_constant(n)
            results.append({
                "num_primes": n,
                "approximation": approx,
                "true_value": self.ARTIN_CONSTANT_APPROX,
                "relative_error": abs(approx - self.ARTIN_CONSTANT_APPROX) / self.ARTIN_CONSTANT_APPROX
            })
        return results

    # -------------------------------------------------------------------
    # 4. HOOLEY'S GRH-RESULTAT
    # -------------------------------------------------------------------

    def hooley_density(self, a: int) -> dict:
        """
        Berechnet Hooley's GRH-basierte Dichte für Basis a.

        Hooley (1967) bewies unter der verallgemeinerten Riemann-Vermutung (GRH):
            Falls a nicht ±1 und kein vollständiges Quadrat, dann gilt:
            #{p ≤ x : a ist primitive Wurzel mod p} ~ C_A · x / ln(x)

        Die Dichte C_A muss ggf. angepasst werden wenn a spezielle Eigenschaften hat:
            - Wenn a = b^k für k ≥ 2: Korrekturfaktor nötig
            - Wenn a negativ: Vorzeichen-Korrektur

        @param a: Basis für die Dichte-Berechnung
        @return: Dictionary mit Hooley-Dichte-Informationen
        @lastModified: 2026-03-12
        """
        # Prüfe ob a ein vollständiges Quadrat ist
        if a > 0:
            sqrt_a = int(math.isqrt(a))
            is_perfect_square = (sqrt_a * sqrt_a == a)
        else:
            is_perfect_square = False

        # Prüfe ob a = -1
        is_minus_one = (a == -1)

        # Standarddichte (Korrekturterme für spezielle Fälle werden vereinfacht)
        density = self.ARTIN_CONSTANT_APPROX

        return {
            "a": a,
            "is_perfect_square": is_perfect_square,
            "is_minus_one": is_minus_one,
            "artin_applies": not is_perfect_square and not is_minus_one,
            "hooley_density": density if not is_perfect_square else 0.0,
            "grh_assumption": "GRH (Generalized Riemann Hypothesis) – nicht bewiesen!",
            "formula": "#{p ≤ x : a primitiv root mod p} ~ C_A · x / ln(x)",
            "note": (
                "CONJECTURE unter GRH. Hooley 1967 bewies dies unter GRH. "
                "Unbedingt (ohne GRH) nur für einige Basen bekannt (Heath-Brown 1986)."
            )
        }

    def heath_brown_result(self) -> dict:
        """
        Dokumentiert Heath-Browns unbedingtes Resultat (1986).

        Heath-Brown bewies ohne GRH:
            Unter den drei Zahlen {2, 3, 5} ist mindestens eine für
            unendlich viele Primzahlen p eine primitive Wurzel.

        Dies ist das stärkste bekannte unbedingte Resultat zur Artin-Vermutung.

        @return: Dictionary mit Heath-Brown-Resultat
        @lastModified: 2026-03-12
        """
        return {
            "author": "D. R. Heath-Brown",
            "year": 1986,
            "result": (
                "Mindestens eine der Zahlen {2, 3, 5} ist primitive Wurzel "
                "für unendlich viele Primzahlen p."
            ),
            "unconditional": True,  # Kein GRH nötig
            "status": "THEOREM (bewiesen)",
            "implication": "Schwächere Form der Artin-Vermutung für spezifische Basen",
            "full_artin": "CONJECTURE (offen, unter GRH von Hooley 1967 bewiesen)"
        }

    # -------------------------------------------------------------------
    # 5. SPEZIFISCHE BASEN a = 2, 3, 5, 7
    # -------------------------------------------------------------------

    def analyze_base(self, a: int, limit: int = 5000) -> dict:
        """
        Vollständige Analyse der Basis a für die Artin-Vermutung.

        Berechnet:
            - Empirische Dichte der Primzahlen mit a als primitive Wurzel
            - Vergleich mit Artin-Konstante
            - Liste der ersten Primzahlen
            - Hooley-Dichte-Schätzung

        @param a: Basis (typisch 2, 3, 5, 7)
        @param limit: Obere Grenze für Primzahlsuche
        @return: Vollständige Analyse als Dictionary
        @lastModified: 2026-03-12
        """
        density_data = self.empirical_density(a, limit)
        hooley_data = self.hooley_density(a)

        return {
            "base": a,
            "analysis": density_data,
            "hooley": hooley_data,
            "is_perfect_square": hooley_data["is_perfect_square"],
            "conjecture_applies": hooley_data["artin_applies"],
            "summary": (
                f"Basis a={a}: {density_data['count_primitive_root']} von "
                f"{density_data['total_primes_checked']} Primzahlen ≤ {limit} "
                f"haben a als primitive Wurzel. "
                f"Empirische Dichte: {density_data['empirical_density']:.4f} "
                f"(Artin-Konstante: {self.ARTIN_CONSTANT_APPROX:.4f})"
            )
        }

    def compare_bases(self, bases: List[int] = None, limit: int = 2000) -> List[dict]:
        """
        Vergleicht mehrere Basen bezüglich ihrer empirischen Dichten.

        @param bases: Liste der Basen (Standard: [2, 3, 5, 7])
        @param limit: Obere Grenze für Primzahlsuche
        @return: Liste von Analyse-Dictionaries, sortiert nach Dichte
        @lastModified: 2026-03-12
        """
        if bases is None:
            bases = [2, 3, 5, 7]

        results = []
        for a in bases:
            density_info = self.empirical_density(a, limit)
            results.append({
                "base": a,
                "density": density_info["empirical_density"],
                "count": density_info["count_primitive_root"],
                "total": density_info["total_primes_checked"],
            })

        # Sortierung nach empirischer Dichte (absteigend)
        return sorted(results, key=lambda x: x["density"], reverse=True)

    def first_primitive_root_primes(self, a: int, count: int = 10) -> List[int]:
        """
        Gibt die ersten 'count' Primzahlen zurück, für die a eine primitive Wurzel ist.

        @param a: Basis
        @param count: Anzahl der gewünschten Primzahlen
        @return: Liste der ersten count Primzahlen mit a als primitive Wurzel
        @lastModified: 2026-03-12
        """
        result = []
        p = 3
        while len(result) < count:
            if isprime(p) and p != abs(a):
                if self.is_primitive_root(a, p):
                    result.append(p)
            p += 2 if p > 2 else 1
        return result

    def conjecture_status(self) -> dict:
        """
        Gibt den aktuellen Beweisstand der Artin-Vermutung zurück.

        STATUS: OFFEN (Conjecture, nicht Theorem – außer unter GRH)

        @return: Dictionary mit Status-Informationen
        @lastModified: 2026-03-12
        """
        return {
            "status": "CONJECTURE – teilweise bewiesen (unter GRH)",
            "statement": (
                "Für jedes a ≠ 0, ±1 das kein vollständiges Quadrat ist, "
                "gibt es unendlich viele Primzahlen p mit ord_p(a) = p−1. "
                "Die Dichte dieser Primzahlen ist C_A ≈ 0.3739558..."
            ),
            "artin_constant": self.ARTIN_CONSTANT_APPROX,
            "grh_result": "Hooley 1967: Unter GRH gilt die Vermutung für alle zulässigen a",
            "unconditional": "Heath-Brown 1986: Gilt für mindestens eine aus {2, 3, 5}",
            "not_proven_in_general": True,
            "year_conjectured": 1927,
        }
