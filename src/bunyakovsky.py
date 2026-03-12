"""
@file bunyakovsky.py
@brief Bunyakovsky-Vermutung: Irreduzible Polynome und unendlich viele Primwerte.
@description
    Dieses Modul implementiert die Bunyakovsky-Vermutung und verwandte
    quantitative Aussagen über primzahlerzeugende Polynome.

    **Bunyakovsky-Vermutung** (Viktor Bunyakovsky, 1857 — Conjecture, OFFEN):
        Sei f(x) ∈ ℤ[x] ein irreduzibles Polynom mit ganzzahligen Koeffizienten
        und positivem Leitkoeffizienten.
        Falls gcd({f(n) : n ∈ ℕ}) = 1 (keine festen Primteiler),
        dann nimmt f(n) unendlich viele Primwerte an.

    **Bewiesene Spezialfälle**:
    - deg(f) = 1: f(x) = ax + b mit gcd(a,b) = 1 → Dirichlet 1837 (BEWEIS)
    - deg(f) ≥ 2: Alle offen, auch f(x) = x² + 1

    **Offene Fälle** (alle Conjectures):
    - f(x) = x² + 1: Primwerte: 2, 5, 17, 37, 41, 53, 61, 73, 89, 97, ...
    - f(x) = x² + x + 1: Primwerte: 3, 7, 13, 31, 43, 73, 157, ...
    - f(x) = x⁴ + 1: Primwerte: 2 (da 1⁴+1=2), 5 (nicht!), ... selten

    **Hardy-Littlewood-Konstante C_f** (heuristische Dichte):
        π_f(x) ~ C_f · x / ln(x)
        C_f = ∏_{p prim} (1 − χ_f(p)/p) / (1 − 1/p)
        mit χ_f(p) = #{a mod p : f(a) ≡ 0 mod p} / p

    **Bateman-Horn-Vermutung** (1962, Conjecture):
        Quantitative Verallgemeinerung der Bunyakovsky-Vermutung auf
        Systeme von Polynomen: Π_f(x) ~ C · x / (ln x)^k

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from sympy import (
    Poly, Symbol, factorint, gcd, isprime, nextprime, primerange, symbols
)
from sympy.polys.polyclasses import DMF


# ===========================================================================
# HAUPTKLASSE: BunyakovskyConjecture
# ===========================================================================

class BunyakovskyConjecture:
    """
    Implementierung der Bunyakovsky-Vermutung und quantitativer Heuristiken.

    Untersucht irreduzible Polynome f(x) ∈ ℤ[x] auf ihre Eigenschaft,
    unendlich viele Primwerte anzunehmen.

    **Bedingungen der Bunyakovsky-Vermutung**:
    1. f(x) ist irreduzibel über ℤ
    2. Der Leitkoeffizient von f ist positiv
    3. gcd({f(n) : n ≥ 1}) = 1 (kein fester Primteiler)

    Wenn alle drei Bedingungen erfüllt sind: VERMUTLICH unendlich viele Primwerte.

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self, poly_func: Callable[[int], int], degree: int = 2,
                 description: str = "f(x)"):
        """
        Initialisiert die Bunyakovsky-Analyse für ein Polynom f.

        @param poly_func: Python-Funktion die f(n) berechnet
        @param degree: Grad des Polynoms
        @param description: Textuelle Beschreibung des Polynoms
        @lastModified: 2026-03-12
        """
        self.f = poly_func
        self.degree = degree
        self.description = description

    @classmethod
    def from_coefficients(cls, coeffs: List[int]) -> "BunyakovskyConjecture":
        """
        Erstellt eine BunyakovskyConjecture-Instanz aus Koeffizientenliste.

        Koeffizienten: [a_0, a_1, ..., a_n] für f(x) = a_0 + a_1·x + ... + a_n·x^n

        @param coeffs: Koeffizientenliste (von niedrigstem zu höchstem Grad)
        @return: BunyakovskyConjecture-Instanz
        @lastModified: 2026-03-12
        """
        def poly_func(n: int) -> int:
            result = 0
            for i, c in enumerate(coeffs):
                result += c * (n ** i)
            return result
        deg = len(coeffs) - 1
        # Beschreibungstext erstellen
        terms = []
        for i in range(len(coeffs) - 1, -1, -1):
            c = coeffs[i]
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                terms.append(f"{c}x" if c != 1 else "x")
            else:
                terms.append(f"{c}x^{i}" if c != 1 else f"x^{i}")
        desc = " + ".join(terms) if terms else "0"
        return cls(poly_func, degree=deg, description=f"f(x) = {desc}")

    def fixed_prime_divisor(self, check_limit: int = 1000) -> Optional[int]:
        """
        Sucht nach festen Primteilern von f (falls vorhanden).

        Ein fester Primteiler p teilt f(n) für ALLE n ∈ ℕ.
        → Dann nimmt f(n) nur endlich viele Primwerte an (nämlich p selbst).

        Methode: Prüfe für kleine Primzahlen p, ob p | f(n) für alle n mod p.

        @param check_limit: Bis zu welcher Primzahl wird geprüft
        @return: Erster fester Primteiler oder None
        @lastModified: 2026-03-12
        """
        for p in primerange(2, check_limit):
            # Prüfe ob p | f(n) für alle n mod p
            all_divisible = True
            for n in range(p):
                if self.f(n) % p != 0:
                    all_divisible = False
                    break
            if all_divisible:
                return int(p)
        return None

    def gcd_of_values(self, limit: int = 1000) -> int:
        """
        Berechnet gcd({f(n) : 1 ≤ n ≤ limit}).

        Falls gcd > 1: f hat einen festen Primteiler → Bunyakovsky-Bedingung verletzt.
        Falls gcd = 1: Bunyakovsky-Bedingung erfüllt (notwendige Bedingung).

        @param limit: Obere Grenze für die gcd-Berechnung
        @return: gcd aller f(n)-Werte für n = 1, ..., limit
        @lastModified: 2026-03-12
        """
        values = [abs(self.f(n)) for n in range(1, limit + 1)]
        result = values[0]
        for v in values[1:]:
            if v == 0:
                continue
            g = math.gcd(result, v)
            result = g
            if result == 1:
                break
        return result

    def find_prime_values(self, n_limit: int = 10_000) -> List[Tuple[int, int]]:
        """
        Sucht alle n ≤ n_limit, für die f(n) prim ist.

        @param n_limit: Obere Grenze für n
        @return: Liste von (n, f(n)) mit f(n) prim
        @lastModified: 2026-03-12
        """
        result = []
        for n in range(1, n_limit + 1):
            val = self.f(n)
            if val > 1 and isprime(val):
                result.append((n, int(val)))
        return result

    def prime_density(self, x: int) -> float:
        """
        Numerische Dichte der Primwerte von f: #{n ≤ x : f(n) prim} / x.

        @param x: Obere Schranke
        @return: Anteil der n ≤ x mit f(n) prim
        @lastModified: 2026-03-12
        """
        prime_values = self.find_prime_values(x)
        return len(prime_values) / x if x > 0 else 0.0

    def hardy_littlewood_constant(self, prime_limit: int = 1000) -> float:
        """
        Berechnet die Hardy-Littlewood-Konstante C_f für das Polynom f.

        Heuristische Dichtformel (Hardy-Littlewood, Conjecture):
            π_f(x) ~ C_f · x / ln(x)  (für deg(f) = 1 korrekt via Dirichlet)
            C_f = ∏_{p prim} (1 − ρ_f(p)/p) / (1 − 1/p)

        ρ_f(p) = #{a ∈ {0,...,p-1} : f(a) ≡ 0 (mod p)} (Anzahl der Nullstellen mod p)

        Für Gradpolynome deg > 1: Conjecture (unbewiesen).

        @param prime_limit: Obere Schranke für das Produkt
        @return: Näherung von C_f
        @lastModified: 2026-03-12
        """
        product = 1.0
        for p in primerange(2, prime_limit):
            p = int(p)
            # ρ_f(p): Anzahl der Nullstellen von f mod p
            rho = sum(1 for a in range(p) if self.f(a) % p == 0)
            # Faktor: (1 − ρ/p) / (1 − 1/p) = (p − ρ) / (p − 1)
            if p == 2:
                # Für p=2: Sonderbehandlung
                factor = (p - rho) / (p - 1) if p > 1 else 1.0
            else:
                factor = (p - rho) / (p - 1)
            if abs(factor) > 1e-15:
                product *= factor
        return product

    def bunyakovsky_criterion_check(self, check_limit: int = 500) -> Dict:
        """
        Prüft alle Bedingungen der Bunyakovsky-Vermutung für f.

        Prüft:
        1. gcd({f(n)}) = 1 (keine festen Primteiler)
        2. f hat keinen festen Primteiler (äquivalent zu 1)
        3. Ob f(n) > 0 für alle n ≥ 1 (Leitkoeffizient positiv)

        Note: Irreduzibilitäts-Prüfung ist für allgemeine Polynome schwieriger
        und wird hier nicht vollständig implementiert.

        @param check_limit: Prüfgrenze
        @return: Dictionary mit Prüfergebnissen
        @lastModified: 2026-03-12
        """
        gcd_val = self.gcd_of_values(check_limit)
        fixed_prime = self.fixed_prime_divisor(check_limit)
        # Überprüfe ob f(n) > 0 für n = 1..100
        positive = all(self.f(n) > 0 for n in range(1, 101))
        prime_vals = self.find_prime_values(min(check_limit, 200))

        return {
            "polynomial": self.description,
            "gcd_of_values": gcd_val,
            "has_fixed_prime_divisor": fixed_prime is not None,
            "fixed_prime_divisor": fixed_prime,
            "bunyakovsky_condition_1": gcd_val == 1,
            "positive_values": positive,
            "prime_count_up_to_200": len(prime_vals),
            "first_prime_values": prime_vals[:10],
            "conjecture_applicable": (gcd_val == 1 and positive),
            "note": (
                "Conjecture: Wenn gcd=1 und irreduzibel, dann ∞ viele Primwerte "
                "(Bunyakovsky 1857 — unbewiesen für deg ≥ 2). "
                "Für deg=1: Beweis via Dirichlet 1837."
            ),
        }

    def count_primes_up_to(self, x: int) -> int:
        """
        Zählt #{n ≤ x : f(n) prim}.

        @param x: Obere Schranke
        @return: Anzahl der Primwerte
        @lastModified: 2026-03-12
        """
        return len(self.find_prime_values(x))

    def compare_with_heuristic(self, x: int) -> Dict:
        """
        Vergleicht tatsächliche Anzahl von Primwerten mit Hardy-Littlewood-Heuristik.

        @param x: Obere Schranke
        @return: Vergleichs-Dictionary
        @lastModified: 2026-03-12
        """
        actual = self.count_primes_up_to(x)
        # Für Grad 1: π_f(x) ~ x/ln(x) · C_f
        # Für Grad d: π_f(x) ~ C_f · x^{1/d} / (ln x) nach Bateman-Horn (Conjecture)
        if x > 1:
            c_f = self.hardy_littlewood_constant(prime_limit=200)
            if self.degree == 1:
                predicted = c_f * x / math.log(x)
            else:
                # Bateman-Horn (Conjecture, unbewiesen):
                predicted = c_f * x ** (1.0 / self.degree) / math.log(x)
        else:
            predicted = 0.0

        return {
            "x": x,
            "actual_prime_count": actual,
            "hardy_littlewood_prediction": predicted,
            "ratio_actual_predicted": actual / predicted if predicted > 0 else None,
        }


# ===========================================================================
# VORDEFINIERTE POLYNOME
# ===========================================================================

def make_x_squared_plus_1() -> BunyakovskyConjecture:
    """
    Erstellt f(x) = x² + 1.

    Primwerte: f(1)=2, f(2)=5, f(4)=17, f(6)=37, f(10)=101, ...
    Vermutung: Unendlich viele (Bunyakovsky-Vermutung, deg=2 — OFFEN).

    @return: BunyakovskyConjecture für x²+1
    @lastModified: 2026-03-12
    """
    return BunyakovskyConjecture(
        lambda n: n * n + 1,
        degree=2,
        description="f(x) = x^2 + 1"
    )


def make_x_squared_plus_x_plus_1() -> BunyakovskyConjecture:
    """
    Erstellt f(x) = x² + x + 1.

    Primwerte: f(1)=3, f(2)=7, f(3)=13, f(5)=31, ...
    Vermutung: Unendlich viele — OFFEN.

    @return: BunyakovskyConjecture für x²+x+1
    @lastModified: 2026-03-12
    """
    return BunyakovskyConjecture(
        lambda n: n * n + n + 1,
        degree=2,
        description="f(x) = x^2 + x + 1"
    )


def make_linear(a: int, b: int) -> BunyakovskyConjecture:
    """
    Erstellt f(x) = ax + b (Dirichlet-Fall, BEWEIS 1837).

    Falls gcd(a,b) = 1: Unendlich viele Primwerte (Dirichlet 1837 — BEWIESEN).
    Falls gcd(a,b) > 1: Nur endlich viele (oder keine Primwerte außer p|gcd).

    @param a: Lineare Koeffizient
    @param b: Konstante
    @return: BunyakovskyConjecture für ax+b
    @lastModified: 2026-03-12
    """
    return BunyakovskyConjecture(
        lambda n: a * n + b,
        degree=1,
        description=f"f(x) = {a}x + {b}"
    )
