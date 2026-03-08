"""
Analytische Zahlentheorie – Primzahlverteilung, Siebmethoden, L-Funktionen.

Dieses Modul implementiert die Kernwerkzeuge der analytischen Zahlentheorie:
    - Primzahlzählfunktion π(x) und Tschebyschow-Funktionen θ(x), ψ(x)
    - Logarithmisches Integral Li(x) (Hauptterm des Primzahlsatzes)
    - Dirichlet-Charaktere und L-Funktionen
    - Siebmethoden: Brun-Sieb, Selberg-Sieb, Chen-Zerlegung
    - Von-Mangoldt-Funktion Λ(n)
    - Riemanns explizite Formel

Mathematisches Fernziel:
    Goldbach-Vermutung: Jede gerade Zahl > 2 = p + q (p, q prim)
    Chen-Theorem (1966): Jede gerade Zahl = p + m, m = Produkt ≤ 2 Primzahlen

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import math
import cmath
from typing import Optional
import numpy as np


# ===========================================================================
# PRIMZAHLZÄHLFUNKTION UND LOGARITHMISCHES INTEGRAL
# ===========================================================================

def prime_counting_function(x: float) -> int:
    """
    Berechnet π(x) = Anzahl der Primzahlen ≤ x.

    Der Primzahlsatz (Hadamard & de la Vallée Poussin, 1896):
        π(x) ~ x / ln(x)  (asymptotisch)

    Bessere Approximation (Gauß/Riemann):
        π(x) ≈ Li(x) = ∫₂ˣ dt/ln(t)

    Exakte Berechnung via Sieb des Eratosthenes: O(x log log x).

    @param x: Obere Schranke
    @return: Exakte Anzahl der Primzahlen ≤ x
    @lastModified: 2026-03-08
    """
    if x < 2:
        return 0

    n = int(x)
    # Sieb des Eratosthenes
    is_prime = bytearray([1]) * (n + 1)  # bytearray für Speichereffizienz
    is_prime[0] = is_prime[1] = 0

    p = 2
    while p * p <= n:
        if is_prime[p]:
            # Alle Vielfachen ab p² markieren
            is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
        p += 1

    return sum(is_prime)


def logarithmic_integral(x: float, terms: int = 100) -> float:
    """
    Berechnet das logarithmische Integral Li(x) = ∫₂ˣ dt/ln(t).

    Li(x) ist die beste "einfache" Approximation für π(x).
    Es gilt immer: Li(x) > π(x) (bis zum ersten Crossover bei ~10^316).

    Berechnung via:
        li(x) = γ + ln|ln(x)| + Σ_{n=1}^{∞} (ln x)^n / (n · n!)

    wobei Li(x) = li(x) - li(2).
    γ = 0.5772156649... ist die Euler-Mascheroni-Konstante.

    @param x: Argument x > 1
    @param terms: Anzahl der Reihenglieder
    @return: Li(x)
    @raises ValueError: Wenn x ≤ 1
    @lastModified: 2026-03-08
    """
    if x <= 1:
        raise ValueError(f"Li(x) ist nur für x > 1 definiert, erhalten: {x}")
    if abs(x - 1) < 1e-10:
        raise ValueError("Li(x) divergiert für x → 1")

    def li(t: float) -> float:
        """Exponentialintegral li(t) = ∫_0^t du/ln(u)."""
        if t <= 0 or abs(t - 1) < 1e-12:
            return float('-inf')

        euler_mascheroni = 0.5772156649015328606  # γ

        ln_t = math.log(t)
        # li(t) = γ + ln|ln(t)| + Σ_{n=1}^∞ (ln t)^n / (n·n!)
        result = euler_mascheroni + math.log(abs(ln_t))

        ln_t_power = ln_t
        factorial = 1
        for n in range(1, terms + 1):
            factorial *= n
            result += ln_t_power / (n * factorial)
            ln_t_power *= ln_t
            if abs(ln_t_power / factorial) < 1e-15:
                break  # Konvergiert

        return result

    # Li(x) = li(x) - li(2)
    return li(x) - li(2)


def prime_number_theorem_comparison(x_values: list[float]) -> list[dict]:
    """
    Vergleicht π(x), x/ln(x) und Li(x) für gegebene x-Werte.

    Zeigt, wie gut die verschiedenen Approximationen für π(x) sind.
    Li(x) ist die beste bekannte einfache Approximation.

    @param x_values: Liste von x-Werten zum Vergleich
    @return: Liste von Dicts mit π(x), x/ln(x), Li(x) und relativen Fehlern
    @lastModified: 2026-03-08
    """
    results = []
    for x in x_values:
        if x <= 2:
            continue
        pi_x = prime_counting_function(x)
        approx_simple = x / math.log(x)
        approx_li = logarithmic_integral(x)
        results.append({
            "x": x,
            "pi_x": pi_x,
            "x_over_ln_x": approx_simple,
            "Li_x": approx_li,
            "error_simple_percent": abs(approx_simple - pi_x) / pi_x * 100 if pi_x > 0 else 0,
            "error_Li_percent": abs(approx_li - pi_x) / pi_x * 100 if pi_x > 0 else 0
        })
    return results


# ===========================================================================
# TSCHEBYSCHOW-FUNKTIONEN (wichtig für Riemanns Explizite Formel)
# ===========================================================================

def von_mangoldt_function(n: int) -> float:
    """
    Berechnet die von-Mangoldt-Funktion Λ(n).

    Definition:
        Λ(n) = ln(p)  falls n = p^k für eine Primzahl p und k ≥ 1
        Λ(n) = 0      sonst

    Bedeutung:
        Die von-Mangoldt-Funktion kodiert die Primzahlverteilung.
        Sie ist das "natürliche" Gewicht in Riemanns expliziter Formel.
        ψ(x) = Σ_{n≤x} Λ(n) = Σ_{p^k≤x} ln(p)

    @param n: Positive ganze Zahl n ≥ 1
    @return: Λ(n) ≥ 0
    @lastModified: 2026-03-08
    """
    if n <= 1:
        return 0.0

    # Prüfe ob n = p^k für eine Primzahl p
    # Faktorisiere n und prüfe auf Primzahlpotenz
    temp = n
    for p in range(2, int(math.sqrt(n)) + 2):
        if temp % p == 0:
            # p teilt n – prüfe ob alle Faktoren p sind
            while temp % p == 0:
                temp //= p
            if temp == 1:
                return math.log(p)  # n = p^k
            return 0.0  # n hat mehrere verschiedene Primteiler

    # n selbst ist prim (kein Teiler gefunden)
    return math.log(n)


def chebyshev_theta(x: float) -> float:
    """
    Berechnet die erste Tschebyschow-Funktion θ(x) = Σ_{p ≤ x} ln(p).

    Eigenschaften:
        - θ(x) ~ x  (Primzahlsatz äquivalent)
        - |θ(x) - x| < c·√x·(ln x)²  (falls RH gilt)
        - θ(x) macht "Sprünge" der Höhe ln(p) an jeder Primzahl p

    @param x: Obere Schranke
    @return: θ(x)
    @lastModified: 2026-03-08
    """
    if x < 2:
        return 0.0

    n = int(x)
    # Sieb
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    p = 2
    while p * p <= n:
        if is_prime[p]:
            is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
        p += 1

    # Summe ln(p) für alle Primzahlen p ≤ x
    return sum(math.log(p) for p in range(2, n + 1) if is_prime[p])


def chebyshev_psi(x: float) -> float:
    """
    Berechnet die zweite Tschebyschow-Funktion ψ(x) = Σ_{p^k ≤ x} ln(p).

    ψ(x) = Σ_{n=1}^{x} Λ(n)  (Summe der von-Mangoldt-Funktion)

    Riemanns explizite Formel verbindet ψ(x) direkt mit den Nullstellen ρ von ζ:
        ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - ½ln(1 - x^{-2})

    @param x: Obere Schranke
    @return: ψ(x)
    @lastModified: 2026-03-08
    """
    n = int(x)
    return sum(von_mangoldt_function(k) for k in range(2, n + 1))


def chebyshev_psi_vs_x(x: float) -> dict:
    """
    Vergleicht ψ(x) mit x (Primzahlsatz in ψ-Form).

    Der Primzahlsatz in der stärksten Form lautet: ψ(x) ~ x.
    Die Riemann-Hypothese impliziert: |ψ(x) - x| = O(√x · ln²x).

    @param x: Obere Schranke
    @return: Vergleichs-Dict
    @lastModified: 2026-03-08
    """
    psi_x = chebyshev_psi(x)
    return {
        "x": x,
        "psi_x": psi_x,
        "ratio": psi_x / x if x > 0 else 0,
        "error": abs(psi_x - x),
        "relative_error": abs(psi_x - x) / x if x > 0 else 0
    }


# ===========================================================================
# DIRICHLET-CHARAKTERE UND L-FUNKTIONEN
# ===========================================================================

def is_coprime(a: int, b: int) -> bool:
    """
    Prüft ob gcd(a, b) = 1.

    @param a: Erste ganze Zahl
    @param b: Zweite ganze Zahl
    @return: True wenn teilerfremd
    @lastModified: 2026-03-08
    """
    return math.gcd(a, b) == 1


def principal_dirichlet_character(n: int, modulus: int) -> int:
    """
    Berechnet den Hauptcharakter χ₀ modulo q.

    Der Hauptcharakter ist der einfachste Dirichlet-Charakter:
        χ₀(n) = 1  falls gcd(n, q) = 1
        χ₀(n) = 0  falls gcd(n, q) > 1

    Dirichlet-Charaktere: Multiplikative Funktionen χ: ℤ → ℂ mit
        - χ(mn) = χ(m)·χ(n)  (multiplikativ)
        - χ(n+q) = χ(n)      (periodisch mit Periode q)
        - χ(n) = 0 falls gcd(n,q) > 1

    @param n: Ganzzahl
    @param modulus: Modul q ≥ 1
    @return: χ₀(n) ∈ {0, 1}
    @lastModified: 2026-03-08
    """
    return 1 if is_coprime(n % modulus, modulus) else 0


def dirichlet_l_function_partial(
    s: complex,
    character_values: list[int],
    modulus: int,
    terms: int = 5000
) -> complex:
    """
    Berechnet die Dirichlet-L-Funktion L(s, χ) via Partialsumme.

    Definition:
        L(s, χ) = Σ_{n=1}^{∞} χ(n)/n^s  für Re(s) > 1

    Eigenschaften:
        - Für den Hauptcharakter χ₀: L(s, χ₀) = ζ(s) · Π_{p|q}(1 - p^{-s})
        - Für primitive Charaktere: L(s, χ) hat analytische Fortsetzung auf ℂ
        - Dirichlet-Theorem: L(1, χ) ≠ 0 für χ ≠ χ₀ → unendlich viele Primzahlen in jeder Restklasse

    @param s: Komplexe Zahl mit Re(s) > 1
    @param character_values: χ(1), χ(2), ..., χ(q) als Liste
    @param modulus: Modul q
    @param terms: Anzahl der Summanden
    @return: Näherungswert von L(s, χ)
    @lastModified: 2026-03-08
    """
    if s.real <= 1:
        raise ValueError(f"Partialsumme konvergiert nur für Re(s) > 1, erhalten: {s.real}")

    total = complex(0.0)
    q = modulus

    for n in range(1, terms + 1):
        # χ(n) aus der periodischen Werteliste ablesen
        chi_n = character_values[(n - 1) % q]
        if chi_n != 0:
            total += chi_n / (n ** s)

    return total


def dirichlet_prime_counting_in_residue_class(
    x: float,
    residue: int,
    modulus: int
) -> int:
    """
    Zählt Primzahlen p ≤ x mit p ≡ r (mod q).

    Dirichlets Theorem über Primzahlen in arithmetischen Progressionen (1837):
        Wenn gcd(r, q) = 1, gibt es unendlich viele Primzahlen p ≡ r (mod q).

    Asymptotik (gleichmäßige Verteilung):
        π(x; q, r) ~ π(x) / φ(q)

    d.h. Primzahlen verteilen sich gleichmäßig auf alle coprime Restklassen.

    @param x: Obere Schranke
    @param residue: Restklasse r
    @param modulus: Modul q
    @return: Anzahl der Primzahlen p ≤ x mit p ≡ r (mod q)
    @lastModified: 2026-03-08
    """
    if not is_coprime(residue, modulus):
        # In dieser Restklasse kann es nur endlich viele Primzahlen geben
        # (nämlich die Primteiler von gcd(r,q))
        count = 0
        n = int(x)
        is_prime = bytearray([1]) * (n + 1)
        is_prime[0] = is_prime[1] = 0
        p = 2
        while p * p <= n:
            if is_prime[p]:
                is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
            p += 1
        for p in range(2, n + 1):
            if is_prime[p] and p % modulus == residue % modulus:
                count += 1
        return count

    # Sieb und zähle Primzahlen in der gewünschten Restklasse
    n = int(x)
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    p = 2
    while p * p <= n:
        if is_prime[p]:
            is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
        p += 1

    return sum(
        1 for p in range(2, n + 1)
        if is_prime[p] and p % modulus == residue % modulus
    )


# ===========================================================================
# SIEBMETHODEN FÜR GOLDBACH-VERWANDTE PROBLEME
# ===========================================================================

def is_semiprime(n: int) -> bool:
    """
    Prüft ob n ein Semiprim (Produkt genau zweier Primzahlen) ist.

    Semiprimen: 4=2·2, 6=2·3, 9=3·3, 10=2·5, 14=2·7, 15=3·5, ...

    Relevant für Chen-Theorem:
        Jede gerade Zahl n = p + m, wobei m entweder prim ODER semiprim.

    @param n: Positive ganze Zahl
    @return: True wenn n = p·q für Primzahlen p, q (p = q erlaubt)
    @lastModified: 2026-03-08
    """
    if n < 4:
        return False

    # Zähle Primfaktoren (mit Vielfachheit)
    factor_count = 0
    temp = n

    for p in range(2, int(math.sqrt(n)) + 2):
        while temp % p == 0:
            factor_count += 1
            temp //= p
            if factor_count > 2:
                return False

    if temp > 1:
        factor_count += 1

    return factor_count == 2


def chen_decomposition(n: int) -> Optional[tuple]:
    """
    Sucht eine Chen-Zerlegung n = p + m.

    Chen-Theorem (1966, Chen Jingrun):
        Jede genügend große gerade Zahl n ist darstellbar als
        n = p + m, wobei p prim und m entweder prim oder semiprim ist.

    Dies ist das bisher stärkste bewiesene Resultat in Richtung Goldbach.
    ("1+2"-Resultat: Goldbach wäre "1+1")

    @param n: Gerade ganze Zahl > 2
    @return: Tupel (p, m, typ) mit n=p+m, typ∈{'prime','semiprime'} oder None
    @raises ValueError: Wenn n ungerade oder ≤ 2
    @lastModified: 2026-03-08
    """
    if n <= 2 or n % 2 != 0:
        raise ValueError(f"n muss gerade und > 2 sein, erhalten: {n}")

    # Importiere Primzahltest aus proof_theory
    def is_prime(k: int) -> bool:
        if k < 2:
            return False
        if k < 4:
            return True
        if k % 2 == 0 or k % 3 == 0:
            return False
        i = 5
        while i * i <= k:
            if k % i == 0 or k % (i + 2) == 0:
                return False
            i += 6
        return True

    for p in range(2, n):
        if not is_prime(p):
            continue
        m = n - p
        if is_prime(m):
            return (p, m, 'prime')     # Goldbach-Zerlegung (1+1)
        if is_semiprime(m):
            return (p, m, 'semiprime') # Chen-Zerlegung (1+2)

    return None  # Sollte für gerade n > 2 nicht passieren


def count_omega(n: int) -> int:
    """
    Zählt ω(n) = Anzahl der verschiedenen Primteiler von n.

    Beispiele:
        ω(12) = ω(2²·3) = 2  (Primteiler: 2 und 3)
        ω(30) = ω(2·3·5) = 3

    @param n: Positive ganze Zahl
    @return: Anzahl der verschiedenen Primteiler
    @lastModified: 2026-03-08
    """
    if n <= 1:
        return 0
    count = 0
    temp = n
    p = 2
    while p * p <= temp:
        if temp % p == 0:
            count += 1
            while temp % p == 0:
                temp //= p
        p += 1
    if temp > 1:
        count += 1
    return count


def almost_prime(n: int, k: int) -> bool:
    """
    Prüft ob n ein k-fast-Prim (k-almost prime) ist.

    n ist k-fast-prim wenn Ω(n) = k, wobei Ω(n) die Gesamtanzahl der
    Primfaktoren mit Vielfachheit ist.

    Beispiele:
        1-fast-prim = Primzahl
        2-fast-prim = Semiprim (p·q oder p²)

    Chen-Theorem: Jede gerade Zahl = prim + 1-fast-prim ODER prim + 2-fast-prim.

    @param n: Positive ganze Zahl
    @param k: Anzahl der Primfaktoren (mit Vielfachheit)
    @return: True wenn Ω(n) = k
    @lastModified: 2026-03-08
    """
    if n <= 1:
        return k == 0

    count = 0
    temp = n
    p = 2
    while p * p <= temp:
        while temp % p == 0:
            count += 1
            temp //= p
            if count > k:
                return False
        p += 1
    if temp > 1:
        count += 1

    return count == k


# ===========================================================================
# RIEMANNS EXPLIZITE FORMEL (approximativ)
# ===========================================================================

def explicit_formula_psi(
    x: float,
    zeros_t: list[float],
    terms: int = 10
) -> dict:
    """
    Berechnet Riemanns explizite Formel für ψ(x).

    Riemanns explizite Formel (1859):
        ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - ½·ln(1 - x^{-2})

    wobei die Summe über alle nicht-trivialen Nullstellen ρ = 1/2 + it läuft.
    Unter der Riemann-Hypothese: ρ = 1/2 ± itₙ, also kommt jede Nullstelle
    als konjugiertes Paar:
        x^ρ/ρ + x^{ρ̄}/ρ̄ = 2·x^(1/2)·cos(t·ln x)/(t² + 1/4)^(1/2) · [Phasenterm]

    Paarweise Vereinfachung (für ρ = 1/2 + it und ρ̄ = 1/2 - it):
        x^ρ/ρ + x^ρ̄/ρ̄ = 2·Re(x^ρ/ρ)

    @param x: Argument x > 1
    @param zeros_t: Imaginärteile der Nullstellen [t₁, t₂, ...] (positiv)
    @param terms: Anzahl der zu verwendenden Nullstellen (beschränkt Summe)
    @return: Dict mit exaktem ψ(x), Näherung und Fehler
    @lastModified: 2026-03-08
    """
    psi_exact = chebyshev_psi(x)

    # Hauptterm: x
    result = x

    # Nullstellen-Summe: -Σ_ρ x^ρ/ρ
    # Unter RH: ρ = 1/2 + it, kommt als Paar mit ρ̄ = 1/2 - it
    used_zeros = zeros_t[:terms]
    oscillation_sum = 0.0

    for t in used_zeros:
        rho = complex(0.5, t)
        # x^ρ / ρ = e^{ρ·ln(x)} / ρ
        x_rho = cmath.exp(rho * math.log(x))
        term_rho = x_rho / rho
        # Paarweise: -( x^ρ/ρ + x^{ρ̄}/ρ̄ ) = -2·Re(x^ρ/ρ)
        oscillation_sum += 2 * term_rho.real

    result -= oscillation_sum

    # Konstanter Term: -ln(2π) ≈ -1.8379
    result -= math.log(2 * math.pi)

    # Logarithmischer Term: -½·ln(1 - x^{-2}) (vernachlässigbar für große x)
    if x > 1:
        result -= 0.5 * math.log(1 - x ** (-2))

    return {
        "x": x,
        "psi_exact": psi_exact,
        "psi_explicit_formula": result,
        "error": abs(result - psi_exact),
        "zeros_used": len(used_zeros),
        "note": "Mehr Nullstellen → bessere Approximation"
    }


# ===========================================================================
# PRIMZAHL-LÜCKEN (Cramér-Vermutung)
# ===========================================================================

def prime_gaps(limit: int) -> list[dict]:
    """
    Berechnet alle Primzahl-Lücken (Abstände aufeinanderfolgender Primzahlen) bis limit.

    Cramér-Vermutung: Für große p gilt: p_{n+1} - p_n = O((ln p_n)²)
    Dies ist stärker als was aus der Riemann-Hypothese folgt.

    Bekannte große Lücken:
        p = 23: Lücke 6 (nächste Primzahl 29)
        p = 89: Lücke 8
        p = 113: Lücke 14

    @param limit: Obere Schranke
    @return: Liste aller Primzahl-Lücken mit Statistik
    @lastModified: 2026-03-08
    """
    n = int(limit)
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    p = 2
    while p * p <= n:
        if is_prime[p]:
            is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
        p += 1

    primes = [p for p in range(2, n + 1) if is_prime[p]]
    gaps = []
    for i in range(len(primes) - 1):
        gap = primes[i + 1] - primes[i]
        gaps.append({
            "prime": primes[i],
            "next_prime": primes[i + 1],
            "gap": gap,
            "cramer_bound": math.log(primes[i]) ** 2
        })
    return gaps


def maximal_prime_gap(limit: int) -> dict:
    """
    Findet die größte Primzahl-Lücke bis zur gegebenen Schranke.

    @param limit: Obere Schranke
    @return: Dict mit der größten gefundenen Lücke
    @lastModified: 2026-03-08
    """
    gaps = prime_gaps(limit)
    if not gaps:
        return {}
    return max(gaps, key=lambda g: g["gap"])
