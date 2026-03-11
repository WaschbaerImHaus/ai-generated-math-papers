"""
@file modular_forms_hecke.py
@brief Hecke-Algebra, Eigenformen und L-Funktionen für Modulformen.
@description
    Implementiert die Hecke-Algebra-Theorie für Modulformen:
    - Hecke-Operatoren T_p (gecacht via lru_cache)
    - Hecke-Algebra-Struktur und Kommutativität
    - Hecke-Eigenformen und Eigenwerte
    - Petersson-Skalarprodukt
    - L-Funktionen von Modulformen (Dirichlet-Reihe + Euler-Produkt)
    - Vollständige L-Funktion mit Funktionalgleichung

    Die Hecke-Algebra ist eine kommutative Algebra von Operatoren auf dem
    Raum M_k der Modulformen vom Gewicht k:
        T_p: M_k → M_k,  (T_p f)(τ) = Σ b(n) q^n

    Für eine Hecke-Eigenform f gilt: T_p f = λ_p · f mit Eigenwert λ_p = a_p/a_1.

    Diese Theorie ist zentral für:
    - den Beweis des Großen Fermatschen Satzes (via Shimura-Taniyama-Wiles)
    - die BSD-Vermutung (L-Funktion elliptischer Kurven)
    - die Ramanujan-Vermutung (|τ(p)| ≤ 2p^{11/2}, bewiesen von Deligne 1974)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
import functools
from typing import Callable


# ===========================================================================
# INTERNE HILFSFUNKTIONEN (Modul-privat)
# ===========================================================================

def _primes_up_to_hecke(n: int) -> list[int]:
    """
    Gibt alle Primzahlen bis n zurück (Sieb des Eratosthenes).

    Interne Hilfsfunktion für die Hecke-Algebra-Berechnung.

    @param n: Obere Grenze
    @return: Liste aller Primzahlen ≤ n
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 2:
        return []
    # Sieb initialisieren: True = Primzahlkandidat
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            # Alle Vielfachen von i als zusammengesetzt markieren
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return [i for i in range(2, n + 1) if sieve[i]]


# ===========================================================================
# HECKE-OPERATOREN
# ===========================================================================

@functools.lru_cache(maxsize=512)
def _hecke_operator_cached(coefficients_tuple: tuple, p: int, k: int) -> tuple:
    """
    Gecachte interne Implementierung des Hecke-Operators.

    Nimmt ein unveränderliches Tupel statt einer Liste an (hashbar, damit
    lru_cache den Aufruf als Schlüssel verwenden kann).

    lru_cache: Hecke-Operatoren sind deterministisch – gleiche Koeffizienten,
    gleiches p und k ergeben immer dasselbe Ergebnis. Der Cache vermeidet
    Neuberechnungen bei wiederholtem Aufruf mit identischen Parametern.

    Der Hecke-Operator T_p wirkt auf Fourier-Koeffizienten a(n) via:
        b(n) = a(pn) + p^{k-1} · a(n/p)   falls p | n
        b(n) = a(pn)                        falls p ∤ n

    @param coefficients_tuple: Fourier-Koeffizienten als hashbares Tupel
    @param p: Primzahl p ≥ 2
    @param k: Gewicht der Modulform
    @return: Transformierte Koeffizienten als Tupel
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    N = len(coefficients_tuple)
    M = max(1, N // p)
    b = []

    for n in range(M):
        # Index p*n (für a(pn))
        idx_pn = p * n
        a_pn = coefficients_tuple[idx_pn] if idx_pn < N else 0

        if n % p == 0:
            # p teilt n: zusätzlicher Term a(n/p) * p^{k-1}
            idx_np = n // p
            a_np = coefficients_tuple[idx_np] if idx_np < N else 0
            bn = a_pn + (p ** (k - 1)) * a_np
        else:
            # p teilt n nicht: nur a(pn)
            bn = a_pn

        b.append(bn)

    return tuple(b)


def hecke_operator(coefficients: list[float | int], p: int, k: int = 12) -> list[float]:
    """
    Wendet den Hecke-Operator T_p auf eine Modulform an.

    Die Modulform ist gegeben durch ihre Fourier-Koeffizienten:
        f(τ) = Σ_{n=0}^∞ a(n) q^n   mit q = e^{2πiτ}

    Der Hecke-Operator T_p wirkt auf Fourier-Koeffizienten via:
        (T_p f)(τ) = Σ_{n=0}^∞ b(n) q^n

    mit:
        b(n) = a(pn) + p^{k-1} · a(n/p)    falls p | n
        b(n) = a(pn)                         falls p ∤ n

    wobei a(n/p) = 0 falls p ∤ n (nicht ganzzahlig).

    Eigenschaften:
        - Hecke-Operatoren kommutieren: T_p ∘ T_q = T_q ∘ T_p für p ≠ q
        - Hecke-Eigenformen: T_p f = λ_p · f mit λ_p = a(p)/a(1)
        - Multiplizität-Eins-Satz: Eigenformen sind durch ihre Eigenwerte bestimmt
        - Verbindung zu L-Funktionen: L(f, s) = Π_p (1 - a(p)p^{-s} + p^{k-1-2s})^{-1}

    Caching: Da Hecke-Operatoren deterministisch sind, werden Ergebnisse via
    lru_cache (in _hecke_operator_cached) gespeichert. Die Liste wird für den
    Cache zu einem Tupel konvertiert.

    @param coefficients: Liste [a(0), a(1), ..., a(N)] der Fourier-Koeffizienten
    @param p: Primzahl p ≥ 2 (Hecke-Operator T_p)
    @param k: Gewicht der Modulform (Standard: 12 für Δ)
    @return: Liste [b(0), b(1), ..., b(M)] der transformierten Koeffizienten
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Liste → Tupel für lru_cache-Kompatibilität (Listen sind nicht hashbar)
    result_tuple = _hecke_operator_cached(tuple(coefficients), p, k)
    return list(result_tuple)


# ===========================================================================
# HECKE-ALGEBRA VERTIEFUNG
# ===========================================================================

def hecke_algebra_structure(k: int, level: int = 1) -> dict[str, object]:
    """
    Untersucht die Struktur der Hecke-Algebra für Modulformen vom Gewicht k.

    Die Hecke-Algebra ist die von den Hecke-Operatoren T_n erzeugte kommutative
    Algebra. Für die Gruppe SL(2,Z) gilt:

        - Kommutativität: T_p ∘ T_q = T_q ∘ T_p für alle Primzahlen p, q
        - Multiplikativität: T_mn = T_m · T_n für gcd(m, n) = 1
        - Rekursion für Primzahlpotenzen: T_{p^{r+1}} = T_p · T_{p^r} - p^{k-1} · T_{p^{r-1}}

    Die Hecke-Eigenwerte für die Ramanujan-Delta-Funktion Δ sind die
    Ramanujan-Tau-Zahlen: T_p Δ = τ(p) · Δ.

    Diese Funktion berechnet:
        1. Hecke-Eigenwerte (tau-Funktion) für kleine Primzahlen
        2. Kommutatorscheck: T_p T_q - T_q T_p = 0 (numerisch auf Koeffizienten)
        3. Multiplikativitätsscheck: τ(p·q) = τ(p)·τ(q) für verschiedene Primzahlen p≠q

    @param k: Gewicht der Modulform (gerade ≥ 4; Standard Delta hat k=12)
    @param level: Niveau (Standard: 1 = volle Modulgruppe SL(2,Z))
    @return: Dictionary mit Hecke-Struktur-Informationen:
             'hecke_eigenvalues': {p: tau(p)} für kleine Primzahlen
             'commutativity_check': {(p,q): bool} – T_p T_q = T_q T_p?
             'multiplicativity_check': {(p,q): bool} – tau(p*q) = tau(p)*tau(q)?
             'is_commutative': bool – Gesamturteil Kommutativität
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Importiere fourier_coefficients_delta aus dem Haupt-Modul (zirkulärer Import
    # wird vermieden, da modular_forms als Wrapper alles re-exportiert)
    from modular_forms import fourier_coefficients_delta

    # Ramanujan-Tau-Koeffizienten bis zu einem hinreichend hohen Index berechnen
    max_prime = 13  # Kleine Primzahlen: 2, 3, 5, 7, 11, 13
    small_primes = _primes_up_to_hecke(max_prime)

    # Tau-Koeffizienten (Hecke-Eigenwerte der Delta-Funktion) berechnen
    # Wir brauchen genug Koeffizienten: bis zu p*q = 13*11 = 143
    n_coeff = 200
    tau_list = fourier_coefficients_delta(n_coeff)

    # Hecke-Eigenwerte: τ(p) für jede Primzahl p ≤ max_prime
    hecke_eigenvalues = {}
    for p in small_primes:
        if p <= n_coeff:
            # tau(p) ist der p-te Fourier-Koeffizient der Delta-Funktion
            hecke_eigenvalues[p] = tau_list[p - 1]

    # Kommutativitätsprüfung: T_p(T_q(f)) == T_q(T_p(f))
    # Wir wenden beide Reihenfolgen auf die tau-Koeffizienten an
    commutativity_check = {}
    for i, p in enumerate(small_primes[:3]):       # Nur erste 3 Primes für Effizienz
        for q in small_primes[i+1:i+3]:            # Nächste zwei Primes
            # Wende T_p dann T_q an
            tpq = hecke_operator(hecke_operator(tau_list, p, k), q, k)
            # Wende T_q dann T_p an
            tqp = hecke_operator(hecke_operator(tau_list, q, k), p, k)
            # Vergleiche erste Koeffizienten
            min_len = min(len(tpq), len(tqp), 10)
            comm = all(tpq[i] == tqp[i] for i in range(min_len))
            commutativity_check[(p, q)] = comm

    # Multiplikativitätsprüfung: τ(p*q) = τ(p) * τ(q) für verschiedene Primzahlen
    multiplicativity_check = {}
    for i, p in enumerate(small_primes[:3]):
        for q in small_primes[i+1:i+3]:
            pq = p * q
            if pq <= n_coeff:
                tau_pq = tau_list[pq - 1]
                tau_p = tau_list[p - 1]
                tau_q = tau_list[q - 1]
                # Multiplikativität: τ(pq) = τ(p) * τ(q) da gcd(p,q)=1 für zwei verschiedene Primzahlen
                multiplicativity_check[(p, q)] = (tau_pq == tau_p * tau_q)

    # Gesamturteil Kommutativität
    is_commutative = all(commutativity_check.values()) if commutativity_check else True

    return {
        'hecke_eigenvalues': hecke_eigenvalues,
        'commutativity_check': commutativity_check,
        'multiplicativity_check': multiplicativity_check,
        'is_commutative': is_commutative,
        'weight': k,
        'level': level,
        'description': (
            'Die Hecke-Algebra ist kommutativ. Eigenwerte sind Ramanujan-Tau-Zahlen. '
            'Multiplikativität τ(mn)=τ(m)τ(n) für gcd(m,n)=1 bestätigt die '
            'Hecke-Eigenform-Eigenschaft der Delta-Funktion.'
        )
    }


def hecke_eigenform(coefficients: list[float], primes: list[int]) -> dict[str, object]:
    """
    Prüft, ob eine Modulform (gegeben durch Fourier-Koeffizienten) eine Hecke-Eigenform ist.

    Eine Hecke-Eigenform ist eine Modulform f mit:
        T_p f = λ_p · f   für alle Primzahlen p

    Für eine normierte Eigenform (a_1 = 1) gilt:
        λ_p = a_p   (p-ter Fourier-Koeffizient)

    Prüfkriterium: Wenn f eine Hecke-Eigenform ist, müssen die Koeffizienten
    der transformierten Form T_p(f) Vielfache der Originalkoeffizienten sein.
    Genauer: (T_p f)_n = a_p · a_n für alle n.

    Für die Delta-Funktion (Ramanujan-Tau): T_p(Δ) = τ(p) · Δ, also ist Δ eine
    normierte Hecke-Eigenform mit Eigenwerten λ_p = τ(p).

    Multiplikativität als notwendige Bedingung:
        a_{mn} = a_m · a_n   für gcd(m,n) = 1

    @param coefficients: Liste [a_1, a_2, ..., a_N] der Fourier-Koeffizienten (ab n=1)
    @param primes: Liste der zu prüfenden Primzahlen
    @return: Dictionary mit:
             'is_eigenform': bool – Gesamturteil (Eigenform?)
             'eigenvalues': {p: lambda_p} – geschätzte Eigenwerte
             'eigenvalue_test': {p: bool} – besteht jede T_p-Prüfung?
             'multiplicativity': {(m,n): bool} – a_{mn} = a_m * a_n für coprime m,n?
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Koeffizienten ab Index 0 = a_1 (Konvention: coefficients[i] = a_{i+1})
    N = len(coefficients)

    # Schätze Eigenwerte: Falls a_1 ≠ 0, ist λ_p = a_p / a_1
    a1 = coefficients[0] if N > 0 else 0

    eigenvalues = {}
    eigenvalue_test = {}

    for p in primes:
        if p > N:
            continue
        # Eigenwert-Schätzung: λ_p = a_p (für normierte Form mit a_1 = 1)
        a_p = coefficients[p - 1] if p - 1 < N else 0

        if a1 != 0:
            # Eigenwert λ_p = a_p / a_1
            eigenvalues[p] = a_p / a1
        else:
            eigenvalues[p] = 0

        # Eigenform-Test: Wende T_p auf die Koeffizienten an
        # Ergebnis sollte λ_p-fache der Originalkoeffizienten sein
        # Interne Koeffizientendarstellung: coefficients[i] = a_{i+1}
        # hecke_operator erwartet [a(0), a(1), ..., a(N)] → füge a(0)=0 vorne an
        full_coeffs = [0] + list(coefficients)
        t_p_result = hecke_operator(full_coeffs, p, k=12)

        # Prüfe ob T_p(f) = λ_p * f (auf verfügbaren Koeffizienten)
        lambda_p = eigenvalues[p]
        is_eigen = True
        for i in range(min(len(t_p_result), min(N, 10))):
            expected = lambda_p * (full_coeffs[i] if i < len(full_coeffs) else 0)
            actual = t_p_result[i]
            # Relative Toleranz für numerischen Vergleich
            if abs(expected) > 1e-10:
                if abs(actual - expected) / abs(expected) > 1e-6:
                    is_eigen = False
                    break
            else:
                if abs(actual - expected) > 1e-6:
                    is_eigen = False
                    break
        eigenvalue_test[p] = is_eigen

    # Multiplikativitätstest: a_{mn} = a_m * a_n für gcd(m,n) = 1
    multiplicativity = {}
    for i, p in enumerate(primes[:3]):
        for q in primes[i+1:i+3]:
            if p != q and p * q <= N:
                a_p = coefficients[p - 1] if p - 1 < N else 0
                a_q = coefficients[q - 1] if q - 1 < N else 0
                a_pq = coefficients[p * q - 1] if p * q - 1 < N else 0
                # Normierung: a_{pq} = a_p * a_q / a_1² (für normierte Form a_1=1)
                if a1 != 0:
                    expected_pq = (a_p * a_q) / (a1)
                else:
                    expected_pq = 0
                # Prüfe Multiplikativität mit Toleranz
                if abs(expected_pq) > 1e-10:
                    mult_ok = abs(a_pq - expected_pq) / abs(expected_pq) < 1e-6
                else:
                    mult_ok = abs(a_pq - expected_pq) < 1e-6
                multiplicativity[(p, q)] = mult_ok

    # Gesamturteil: Eigenform wenn alle T_p-Tests bestehen
    is_eigenform = all(eigenvalue_test.values()) if eigenvalue_test else False

    return {
        'is_eigenform': is_eigenform,
        'eigenvalues': eigenvalues,
        'eigenvalue_test': eigenvalue_test,
        'multiplicativity': multiplicativity,
        'a1_coefficient': a1
    }


def petersson_inner_product_estimate(f_coeffs: list[float], g_coeffs: list[float], k: int) -> complex:
    """
    Schätzt das Petersson-Skalarprodukt ⟨f, g⟩ numerisch.

    Das Petersson-Skalarprodukt ist definiert als:
        ⟨f, g⟩ = ∫∫_F f(τ) · conj(g(τ)) · (Im τ)^k · dx dy / y²

    wobei F der Fundamentalbereich von SL(2,Z) ist:
        F = {τ ∈ H : |Re(τ)| ≤ 1/2, |τ| ≥ 1, Im(τ) > 0}

    Numerische Approximation via Rankin-Selberg-Methode (für zwei Cusp-Formen):
        ⟨f, g⟩ ≈ Σ_n a_n · conj(b_n) · (4πn)^{-(k-1)} · Γ(k-1) / normierung

    Da diese Formel auf dem Rankkin-Selberg-Integral basiert und die genaue
    Normierung vom Niveau N abhängt, gibt diese Implementierung eine rohe Schätzung zurück.

    @param f_coeffs: Fourier-Koeffizienten [a_1, a_2, ..., a_N] von f
    @param g_coeffs: Fourier-Koeffizienten [b_1, b_2, ..., b_N] von g
    @param k: Gewicht der Modulformen (muss für beide gleich sein)
    @return: Schätzwert für ⟨f, g⟩ als komplexe Zahl
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Rankin-Selberg-Approximation:
    # ⟨f, g⟩ ~ Σ_{n=1}^{N} a_n * conj(b_n) * weight(n)
    # mit weight(n) = n^{-(k-1)} * Γ(k-1) (vereinfacht)

    N = min(len(f_coeffs), len(g_coeffs))
    if N == 0:
        return complex(0.0)

    # Gammafunktion Γ(k-1) (für k ≥ 2: Γ(k-1) = (k-2)!)
    if k >= 2:
        gamma_k_minus_1 = math.factorial(k - 2) if k >= 2 else 1.0
    else:
        gamma_k_minus_1 = 1.0

    # Normierungsfaktor (vereinfacht, ohne genauen Fundamentalbereich-Faktor)
    norm_factor = gamma_k_minus_1 / ((4 * math.pi) ** (k - 1))

    # Summiere Beiträge: Σ a_n * conj(b_n) * n^{-(k-1)}
    result = complex(0.0)
    for n in range(1, N + 1):
        a_n = complex(f_coeffs[n - 1])
        b_n = complex(g_coeffs[n - 1])
        # Gewichtsfaktor: n^{-(k-1)}
        weight = n ** (-(k - 1))
        result += a_n * b_n.conjugate() * weight

    return result * norm_factor


# ===========================================================================
# L-FUNKTIONEN VON MODULFORMEN
# ===========================================================================

def l_function_modular_form(coefficients: list[float], s: complex, k: int) -> complex:
    """
    Berechnet die L-Funktion einer Modulform f = Σ a_n q^n.

    Definition:
        L(f, s) = Σ_{n=1}^∞ a_n / n^s

    Die Reihe konvergiert absolut für Re(s) > (k+1)/2 (nach dem Phragmen-Lindelöf-Prinzip).

    Für eine Hecke-Eigenform hat L(f, s) ein Euler-Produkt:
        L(f, s) = Π_p (1 - a_p · p^{-s} + p^{k-1-2s})^{-1}

    Für die Ramanujan-Delta-Funktion (k=12):
        L(Δ, s) = Π_p (1 - τ(p) · p^{-s} + p^{11-2s})^{-1}

    Verbindung zur vollständigen L-Funktion:
        Λ(f, s) = (√N / 2π)^s · Γ(s) · L(f, s) = ε · Λ(f, k-s)

    @param coefficients: Liste [a_1, a_2, ..., a_N] der Fourier-Koeffizienten
    @param s: Komplexes Argument der L-Funktion (Re(s) > (k+1)/2)
    @param k: Gewicht der Modulform (für Konvergenzhinweis)
    @return: Wert L(f, s) als komplexe Zahl (Partialsumme)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    result = complex(0.0)
    N = len(coefficients)

    # Dirichlet-Reihe: L(f, s) = Σ_{n=1}^{N} a_n / n^s
    for n in range(1, N + 1):
        a_n = complex(coefficients[n - 1])
        if a_n != 0:
            # Summand: a_n / n^s
            result += a_n / (n ** s)

    return result


def functional_equation_l_function(coefficients: list[float], k: int, N: int, s: complex) -> dict[str, object]:
    """
    Berechnet die vollständige L-Funktion und verifiziert die Funktionalgleichung.

    Vollständige L-Funktion (completed L-function):
        Λ(f, s) = (√N / (2π))^s · Γ(s) · L(f, s)

    Funktionalgleichung:
        Λ(f, s) = ε · Λ(f̄, k - s)

    wobei:
        ε = ±1: Vorzeichen (Atkin-Lehner-Eigenwert)
        f̄: Komplex-konjugierte Modulform (für reelle Koeffizienten: f̄ = f)
        N: Niveau (conductor) der Modulform
        k: Gewicht

    Für normierte Hecke-Eigenformen mit reellen Koeffizienten (wie Δ):
        Λ(f, s) = ε · Λ(f, k - s)

    Das Vorzeichen ε hängt von der Atkin-Lehner-Involution ab:
        ε = (-1)^{k/2} für SL(2,Z) (Niveau 1)

    @param coefficients: Fourier-Koeffizienten [a_1, ..., a_M]
    @param k: Gewicht der Modulform
    @param N: Niveau (conductor), Standardfall: N = 1 für SL(2,Z)
    @param s: Komplexes Argument
    @return: Dictionary mit:
             'lambda_s': Λ(f, s) – vollständige L-Funktion bei s
             'lambda_k_minus_s': Λ(f, k-s) – vollständige L-Funktion bei k-s
             'epsilon': Vorzeichen ε = (-1)^{k/2}
             'functional_equation_ratio': Λ(f,s) / (ε · Λ(f, k-s))
             'functional_equation_verified': bool – |ratio - 1| < 0.1?
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Vorzeichen: ε = (-1)^{k/2} für die volle Modulgruppe SL(2,Z)
    epsilon = (-1) ** (k // 2)

    def gamma_approx(z: complex) -> complex:
        """
        Näherung der Gamma-Funktion via Stirling-Formel.

        Für Re(z) > 0.5: Stirling-Approximation.
        Für Re(z) ≤ 0.5: Reflexionsformel Γ(z)·Γ(1-z) = π/sin(πz).

        @param z: Komplexes Argument
        @return: Näherungswert von Γ(z)
        @lastModified: 2026-03-11
        """
        if z.real > 0.5:
            # Stirling: Γ(z) ≈ sqrt(2π/z) · (z/e)^z
            return cmath.sqrt(2 * math.pi / z) * (z / math.e) ** z
        else:
            # Reflexionsformel: Γ(z) · Γ(1-z) = π / sin(πz)
            return math.pi / (cmath.sin(math.pi * z) * gamma_approx(1 - z))

    # Vorfaktor: (√N / (2π))^s
    prefactor_s = (math.sqrt(N) / (2 * math.pi)) ** s
    prefactor_ks = (math.sqrt(N) / (2 * math.pi)) ** (k - s)

    # L-Funktion bei s und k-s
    l_s = l_function_modular_form(coefficients, s, k)
    l_ks = l_function_modular_form(coefficients, k - s, k)

    # Gamma-Faktoren
    gamma_s = gamma_approx(s)
    gamma_ks = gamma_approx(k - s)

    # Vollständige L-Funktionen: Λ(f, s) = Vorfaktor · Γ(s) · L(f, s)
    lambda_s = prefactor_s * gamma_s * l_s
    lambda_ks = prefactor_ks * gamma_ks * l_ks

    # Verhältnis: Λ(f,s) / (ε · Λ(f, k-s)) soll = 1 sein
    if abs(lambda_ks) > 1e-30:
        ratio = lambda_s / (epsilon * lambda_ks)
        verified = abs(ratio - 1.0) < 0.1  # Toleranz wegen Partialsummen-Approximation
    else:
        ratio = complex(float('nan'))
        verified = False

    return {
        'lambda_s': lambda_s,
        'lambda_k_minus_s': lambda_ks,
        'epsilon': epsilon,
        'functional_equation_ratio': ratio,
        'functional_equation_verified': verified,
        'l_s': l_s,
        'l_k_minus_s': l_ks,
        's': s,
        'k': k,
        'N': N
    }
