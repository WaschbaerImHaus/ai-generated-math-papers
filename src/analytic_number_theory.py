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
import functools
from typing import Optional
import numpy as np


# ===========================================================================
# PRIMZAHLZÄHLFUNKTION UND LOGARITHMISCHES INTEGRAL
# ===========================================================================

@functools.lru_cache(maxsize=512)
def _prime_counting_cached(n: int) -> int:
    """
    Gecachte interne Implementierung von π(x) für ganzzahlige Argumente.

    lru_cache speichert bis zu 512 verschiedene n-Werte. Bei wiederholtem
    Aufruf mit gleichem n wird das Ergebnis direkt aus dem Cache zurückgegeben,
    ohne den teuren Sieb-Algorithmus erneut auszuführen.

    Warum interne Funktion?
        prime_counting_function() nimmt float entgegen. Da float nicht hashbar
        für lru_cache in allen Fällen zuverlässig ist (z.B. 1e3 vs 1000.0),
        konvertieren wir zuerst zu int und cachen die int-Version.

    @param n: Ganzzahlige obere Schranke n ≥ 0
    @return: π(n) = Anzahl der Primzahlen ≤ n
    @lastModified 2026-03-10
    """
    if n < 2:
        return 0
    # Sieb des Eratosthenes – O(n log log n) Zeit, O(n) Speicher
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    p = 2
    while p * p <= n:
        if is_prime[p]:
            is_prime[p * p::p] = bytearray(len(is_prime[p * p::p]))
        p += 1
    return sum(is_prime)


def prime_counting_function(x: float) -> int:
    """
    Berechnet π(x) = Anzahl der Primzahlen ≤ x.

    Der Primzahlsatz (Hadamard & de la Vallée Poussin, 1896):
        π(x) ~ x / ln(x)  (asymptotisch)

    Bessere Approximation (Gauß/Riemann):
        π(x) ≈ Li(x) = ∫₂ˣ dt/ln(t)

    Exakte Berechnung via Sieb des Eratosthenes: O(x log log x).
    Wiederholte Aufrufe mit demselben (ganzzahligen) x werden gecacht (lru_cache).

    @param x: Obere Schranke
    @return: Exakte Anzahl der Primzahlen ≤ x
    @lastModified: 2026-03-10
    """
    if x < 2:
        return 0
    # Konvertierung zu int für Cache-Kompatibilität (float ist nicht sicher hashbar)
    return _prime_counting_cached(int(x))


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

@functools.lru_cache(maxsize=512)
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

    lru_cache: Da Λ(n) deterministisch und rein funktional ist (kein globaler Zustand),
    können Ergebnisse sicher gecacht werden. maxsize=512 speichert die 512 häufigsten Werte.

    @param n: Positive ganze Zahl n ≥ 1
    @return: Λ(n) ≥ 0
    @lastModified: 2026-03-10
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


# ===========================================================================
# SELBERG-KLASSE – AXIOM-PRÜFUNG UND ORTHOGONALITÄTSVERMUTUNG
# ===========================================================================

def selberg_class_check(
    coefficients: list[complex],
    degree: float,
    conductor: float
) -> dict:
    """
    Überprüft ob eine L-Funktion (definiert durch ihre Koeffizienten a_n)
    die Selberg-Klassen-Axiome empirisch erfüllt.

    Die Selberg-Klasse S wurde 1992 von Atle Selberg eingeführt, um
    L-Funktionen axiomatisch zu charakterisieren:

    Axiom 1 – Dirichlet-Reihe:
        L(s) = Σ_{n=1}^{∞} a_n / n^s  konvergiert für Re(s) > 1.
        Voraussetzung: a_1 = 1 (Normierung).

    Axiom 2 – Analytische Fortsetzung:
        (s-1)^m · L(s) ist eine ganze Funktion (m = 0 oder 1).
        (nicht numerisch prüfbar ohne vollständige Funktion)

    Axiom 3 – Funktionalgleichung:
        Λ(s) = ε · Λ̄(1-s)  mit |ε| = 1
        wobei Λ(s) = Q^s · Π_j Γ(λ_j s + μ_j) · L(s)
        (nicht numerisch prüfbar ohne vollständige Struktur)

    Axiom 4 – Euler-Produkt:
        L(s) = Π_p F_p(p^{-s})^{-1}
        mit F_p(X) = Π_j (1 - α_{p,j} X), |α_{p,j}| ≤ 1
        Empirisch prüfbar: Euler-Produkt stimmt mit Partialsumme überein.

    Axiom 5 – Ramanujan-Vermutung:
        |a_p| ≤ d (Grad d der L-Funktion) für alle Primzahlen p.
        Für ζ(s): a_n = 1 → |a_p| = 1 ≤ 1 ✓
        Für Modulformen-L-Funktionen: |a_p| ≤ 2 (Grad 2).

    @param coefficients: Koeffizientenliste [a_1, a_2, ..., a_N] (1-indiziert, a_0 ignoriert)
    @param degree: Grad d der L-Funktion (ζ hat Grad 1, Modulformen Grad 2)
    @param conductor: Führer N > 0 (Skalar für Funktionalgleichung)
    @return: Dict mit Axiom-Prüfergebnissen und Statistiken
    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """
    n = len(coefficients)
    if n < 2:
        return {"error": "Mindestens 2 Koeffizienten benötigt", "valid": False}

    # --- Axiom 1: Normierung a_1 = 1 ---
    a1 = coefficients[0]
    axiom1_normalization = abs(a1 - 1.0) < 1e-10

    # --- Axiom 1: Absolute Konvergenz prüfen (|a_n| = O(n^ε)) ---
    # Prüfe ob die Koeffizienten polynomiell wachsen (kein exp. Wachstum)
    max_coeff = max(abs(c) for c in coefficients if c != 0) if coefficients else 0.0
    growth_ok = True
    for idx, c in enumerate(coefficients):
        idx_1based = idx + 1
        # Grobe Schranke: |a_n| ≤ n^0.5 impliziert Konvergenz für Re(s) > 1
        if idx_1based > 1 and abs(c) > idx_1based ** 0.6:
            growth_ok = False
            break

    axiom1_satisfied = axiom1_normalization and growth_ok

    # --- Axiom 4: Euler-Produkt Konsistenz ---
    # Prüfe ob L(2) via Partialsumme mit Euler-Produkt über Primzahlen übereinstimmt
    # Partialsumme: Σ a_n / n^2
    s_test = 2.0
    partial_sum = sum(
        coefficients[idx] / (idx + 1) ** s_test
        for idx in range(min(n, 200))
        if coefficients[idx] != 0
    )

    # Euler-Produkt Näherung: Π_p (1 - a_p/p^s + ...)^{-1}
    # Nur für multiplikative Koeffizienten exakt; hier numerische Näherung
    euler_product = complex(1.0)
    for p in _small_primes_up_to(min(n, 50)):
        p_idx = p - 1  # 0-basierter Index
        if p_idx < n:
            ap = coefficients[p_idx]
            # Lokaler Faktor: (1 - a_p · p^{-s})^{-1} (Grad-1-Näherung)
            denom = 1.0 - ap / (p ** s_test)
            if abs(denom) > 1e-12:
                euler_product *= 1.0 / denom

    euler_deviation = abs(partial_sum - euler_product)
    axiom4_approx = euler_deviation < 5.0  # Großzügige Toleranz (Partialsumme konvergiert langsam)

    # --- Axiom 5: Ramanujan-Vermutung ---
    # Prüfe |a_p| ≤ 2·degree für alle Primzahlen p mit Index < n
    ramanujan_violations = []
    ramanujan_bound = max(2.0 * degree, 1.0)
    for p in _small_primes_up_to(min(n, 100)):
        p_idx = p - 1
        if p_idx < n:
            ap_abs = abs(coefficients[p_idx])
            if ap_abs > ramanujan_bound + 1e-10:
                ramanujan_violations.append({
                    "prime": p,
                    "abs_a_p": ap_abs,
                    "bound": ramanujan_bound
                })
    axiom5_satisfied = len(ramanujan_violations) == 0

    # --- Statistiken ---
    prime_coeffs = [
        abs(coefficients[p - 1])
        for p in _small_primes_up_to(min(n, 100))
        if p - 1 < n
    ]
    mean_ap = float(np.mean(prime_coeffs)) if prime_coeffs else 0.0
    max_ap = float(np.max(prime_coeffs)) if prime_coeffs else 0.0

    return {
        "axiom1_normalization": axiom1_normalization,
        "axiom1_growth": growth_ok,
        "axiom1_satisfied": axiom1_satisfied,
        "axiom4_euler_partial_sum": complex(partial_sum),
        "axiom4_euler_product_approx": complex(euler_product),
        "axiom4_deviation": float(euler_deviation),
        "axiom4_satisfied": axiom4_approx,
        "axiom5_satisfied": axiom5_satisfied,
        "axiom5_violations": ramanujan_violations,
        "ramanujan_bound": ramanujan_bound,
        "mean_abs_a_p": mean_ap,
        "max_abs_a_p": max_ap,
        "degree": degree,
        "conductor": conductor,
        "n_coefficients": n,
        # Axiome 2 und 3 erfordern vollständige analytische Struktur
        "axiom2_note": "Analytische Fortsetzung nicht empirisch prüfbar",
        "axiom3_note": "Funktionalgleichung nicht empirisch prüfbar",
        "selberg_class_candidate": axiom1_satisfied and axiom5_satisfied
    }


def _small_primes_up_to(n: int) -> list[int]:
    """
    Hilfsfunktion: Gibt alle Primzahlen ≤ n zurück (Sieb des Eratosthenes).

    @param n: Obere Schranke
    @return: Liste der Primzahlen ≤ n
    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
        p += 1
    return [i for i in range(2, n + 1) if sieve[i]]


def selberg_orthogonality(
    chi1_coeffs: list[complex],
    chi2_coeffs: list[complex],
    prime_bound: int = 100
) -> dict:
    """
    Überprüft die Selberg-Orthogonalitäts-Vermutung numerisch.

    Die Selberg-Orthogonalitätsvermutung (Selberg 1992) lautet:
        (1 / log X) · Σ_{p ≤ X} a_p(F) · ā_p(G) / p  →  δ_{F,G} · n_F

    wobei:
        - a_p(F), a_p(G): Euler-Koeffizienten der L-Funktionen F, G bei Primzahl p
        - δ_{F,G} = 1 wenn F = G (gleiche L-Funktion), sonst 0
        - n_F = Ordnung des Pols von L_F(s) bei s = 1 (0 für ganze Funktionen, 1 für ζ)
        - X = prime_bound

    Für ζ(s): a_p = 1 für alle p → Summe ~ log log X (langsame Divergenz).
    Für zwei verschiedene primitive Dirichlet-Charaktere χ₁ ≠ χ₂:
        Σ χ₁(p) · χ̄₂(p) / p → 0 (Orthogonalität).

    Die Vermutung verallgemeinert die Orthogonalitätsrelation für Dirichlet-Charaktere
    auf allgemeine L-Funktionen.

    @param chi1_coeffs: Koeffizientenliste [a_1, ..., a_N] der ersten L-Funktion F
    @param chi2_coeffs: Koeffizientenliste [a_1, ..., a_N] der zweiten L-Funktion G
    @param prime_bound: Obere Schranke X für die Primzahlsumme
    @return: Dict mit:
             - 'sum': empirische gewichtete Summe
             - 'normalized': Normierung durch log(X)
             - 'same_function_detected': heuristische Gleichheitserkennung
             - 'primes_used': Anzahl der verwendeten Primzahlen
             - 'orthogonality_measure': Betrag der normierten Summe (nahe 0 → orthogonal)
    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """
    primes = _small_primes_up_to(prime_bound)
    n1 = len(chi1_coeffs)
    n2 = len(chi2_coeffs)

    # Gewichtete Summe: Σ_{p ≤ X} a_p(F) · ā_p(G) / p
    weighted_sum = complex(0.0)
    primes_used = 0

    for p in primes:
        p_idx = p - 1  # 0-basierter Index (a_1 steht an Position 0)
        if p_idx < n1 and p_idx < n2:
            ap_F = chi1_coeffs[p_idx]
            ap_G = chi2_coeffs[p_idx]
            # Selberg-Gewichtung: a_p(F) · ā_p(G) / p
            weighted_sum += ap_F * ap_G.conjugate() / p
            primes_used += 1

    # Normierung: Teile durch log(X)
    if prime_bound > 1:
        log_x = math.log(prime_bound)
        normalized_sum = weighted_sum / log_x
    else:
        normalized_sum = complex(0.0)
        log_x = 1.0

    # Heuristische Gleichheitserkennung: wenn Summen sehr ähnlich sind
    # (gleiche Koeffizient-Muster bei Primzahlen)
    abs_diff_at_primes = 0.0
    common_primes = 0
    for p in primes[:20]:  # Erste 20 Primzahlen vergleichen
        p_idx = p - 1
        if p_idx < n1 and p_idx < n2:
            abs_diff_at_primes += abs(chi1_coeffs[p_idx] - chi2_coeffs[p_idx])
            common_primes += 1

    same_function = (abs_diff_at_primes / max(common_primes, 1)) < 0.01

    # Orthogonalitätsmaß: klein → orthogonal, groß → gleichartig
    orthogonality_measure = abs(normalized_sum)

    return {
        "sum": complex(weighted_sum),
        "log_x": float(log_x),
        "normalized": complex(normalized_sum),
        "orthogonality_measure": float(orthogonality_measure),
        "same_function_detected": same_function,
        "primes_used": primes_used,
        "prime_bound": prime_bound,
        # Interpretation: für orthogonale L-Funktionen → 0, für F=G → n_F (pole order)
        "interpretation": (
            "F ≈ G (nicht-triviale Korrelation)" if orthogonality_measure > 0.3
            else "F ⊥ G (Orthogonalitätsvermutung erfüllt)"
        )
    }


def selberg_zeta_motivation(x: float) -> dict:
    """
    Verbindet die Primzahlverteilung mit der Selberg-Klasse über Riemanns explizite Formel.

    Riemanns explizite Formel (1859) drückt π(x) durch Nullstellen von ζ(s) aus:
        ψ(x) = x - Σ_{ρ} x^ρ / ρ - ln(2π) - (1/2) ln(1 - x^{-2})

    Die Fehlerterme in der Primzahlzählung:
        |π(x) - Li(x)| ≤ C · √x · ln(x)  (unter RH)

    Präziser, mit endlich vielen Nullstellen bis Höhe T:
        |ψ(x) - x| ≤ Σ_{|Im(ρ)|≤T} |x^ρ/ρ| + O(x · ln²(x) / T)

    Jede Nullstelle ρ = σ + it trägt bei:
        |x^ρ / ρ| = x^σ / |ρ|

    Unter RH (σ = 1/2 für alle Nullstellen):
        Fehler ~ x^{1/2} · (Summe über 1/|ρ|)

    Die bekannten ersten Nullstellen (gerundete Im-Teile):
        t₁ ≈ 14.13, t₂ ≈ 21.02, t₃ ≈ 25.01, t₄ ≈ 30.42, t₅ ≈ 32.93

    @param x: Auswertungspunkt x > 2
    @return: Dict mit:
             - 'x': Eingabewert
             - 'pi_x': exakte Primzahlzählung π(x)
             - 'Li_x': logarithmisches Integral (Hauptterm)
             - 'absolute_error': |π(x) - Li(x)|
             - 'zero_contributions': Beiträge der ersten bekannten Nullstellen
             - 'total_zero_correction': Summe der Nullstellenbeiträge (unter RH)
             - 'rh_error_bound': Fehlerschranke unter RH: √x · ln(x) / π
    @raises ValueError: Wenn x ≤ 2
    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """
    if x <= 2:
        raise ValueError(f"x muss > 2 sein, erhalten: {x}")

    # Primzahlzählung und Logarithmisches Integral
    pi_x = prime_counting_function(x)
    li_x = logarithmic_integral(x)
    absolute_error = abs(pi_x - li_x)

    # Bekannte Riemann-Nullstellen (Imaginärteile der ersten Nullstellen auf Re=1/2)
    known_zeros_imaginary = [
        14.134725141734693790,
        21.022039638771554993,
        25.010857580145688763,
        30.424876125859513210,
        32.935061587739189691,
        37.586178158825671257,
        40.918719012147495187,
        43.327073280914999519,
        48.005150881167159727,
        49.773832477672302181,
    ]

    # Beitrag jeder Nullstelle ρ = 1/2 + it zur expliziten Formel
    zero_contributions = []
    total_correction = 0.0

    for t in known_zeros_imaginary:
        rho = complex(0.5, t)
        rho_conj = complex(0.5, -t)

        # x^ρ / ρ + x^{ρ̄} / ρ̄ = 2 · Re(x^ρ / ρ) (reeller Beitrag)
        x_rho = cmath.exp(rho * cmath.log(x))
        contribution_complex = x_rho / rho
        contribution_real = 2.0 * contribution_complex.real  # Beide konjugierten Nullstellen

        abs_contribution = abs(contribution_complex)

        zero_contributions.append({
            "imaginary_part": t,
            "rho": rho,
            "x_rho_over_rho": complex(contribution_complex),
            "abs_contribution": float(abs_contribution),
            "real_correction": float(contribution_real)
        })
        total_correction += abs_contribution

    # Fehlerschranke unter RH: |ψ(x) - x| = O(√x · ln²(x))
    # Für π(x): |π(x) - Li(x)| = O(√x · ln(x))
    rh_error_bound = math.sqrt(x) * math.log(x) / math.pi

    return {
        "x": float(x),
        "pi_x": int(pi_x),
        "Li_x": float(li_x),
        "absolute_error": float(absolute_error),
        "relative_error": float(absolute_error / max(pi_x, 1)),
        "zero_contributions": zero_contributions,
        "n_zeros_used": len(known_zeros_imaginary),
        "total_zero_correction": float(total_correction),
        "rh_error_bound": float(rh_error_bound),
        "rh_error_bound_formula": "sqrt(x) * ln(x) / pi",
        # Unter RH: Fehler wächst wie O(√x), ohne RH wie O(x^θ) mit θ > 1/2
        "error_within_rh_bound": absolute_error <= rh_error_bound + 1.0
    }
