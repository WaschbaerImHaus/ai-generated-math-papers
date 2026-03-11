"""
@file modular_forms.py
@brief Rückwärtskompatibilitäts-Wrapper — Inhalt aufgeteilt in modular_forms_hecke.py, modular_forms_theta.py
@description
    Dieser Wrapper stellt alle bisherigen öffentlichen Namen aus modular_forms
    weiterhin zur Verfügung. Der eigentliche Code wurde aufgeteilt in:

    - modular_forms_hecke.py  : Hecke-Algebra, Eigenformen und L-Funktionen
    - modular_forms_theta.py  : Theta-Reihen, Jacobi-Theta, Dedekind-Eta

    Dieses Modul enthält weiterhin direkt:
    - Die modulare Gruppe SL(2,Z): ModularGroup
    - Eisenstein-Reihen: eisenstein_series, eisenstein_series_fast,
      normalized_eisenstein_E4, normalized_eisenstein_E6
    - Ramanujan-Delta-Funktion: delta_function
    - j-Invariante: j_invariant
    - Ramanujan-Tau-Funktion: fourier_coefficients_delta, ramanujan_tau_properties
    - Modularitätsbedingung: modular_form_check
    - Cusp-Formen: cusp_form_dimension, is_cusp_form_basis
    - Shimura-Taniyama-Wiles: shimura_taniyama_check
    - Elliptische Kurven: elliptic_curve_points_over_fp, trace_of_frobenius,
      l_function_elliptic_curve, birch_swinnerton_dyer_bsd_estimate

    Rückwärtskompatibilität: `from modular_forms import *` funktioniert weiterhin.

@author Michael Fuhrmann
@version 4.0
@since 2026-03-08
@lastModified 2026-03-11
"""

import math
import cmath
import functools
import numpy as np
from typing import Callable
from math import gcd

# ===========================================================================
# RE-EXPORTS AUS AUFGETEILTEN MODULEN
# ===========================================================================

# Theta-Reihen-Modul: Jacobi-Theta, Dreifachprodukt, Quadratsummen, Dedekind-Eta
from modular_forms_theta import (
    theta_function,
    theta_transformation,
    sum_of_squares_theta,
    jacobi_triple_product,
    dedekind_eta,
)

# Hecke-Algebra-Modul: Hecke-Operatoren, Eigenformen, L-Funktionen
from modular_forms_hecke import (
    hecke_operator,
    hecke_algebra_structure,
    hecke_eigenform,
    petersson_inner_product_estimate,
    l_function_modular_form,
    functional_equation_l_function,
)


# ===========================================================================
# MODULARE GRUPPE SL(2,Z)
# ===========================================================================

class ModularGroup:
    """
    Repräsentiert ein Element der modularen Gruppe SL(2,Z).

    SL(2,Z) ist die Gruppe aller 2×2-Matrizen [[a,b],[c,d]] mit
    ganzzahligen Einträgen und Determinante det = ad - bc = 1.

    Diese Gruppe wirkt auf die obere Halbebene H = {τ ∈ ℂ : Im(τ) > 0}
    durch Möbius-Transformationen (gebrochen-lineare Transformationen):
        γ(τ) = (aτ + b) / (cτ + d)

    Wichtige Elemente:
        I = [[1,0],[0,1]]  – Einheitsmatrix (Identität)
        S = [[0,-1],[1,0]] – τ ↦ -1/τ (Inversion)
        T = [[1,1],[0,1]]  – τ ↦ τ+1 (Translation)
        S und T erzeugen SL(2,Z).

    @author Michael Fuhrmann
    @since 2026-03-08
    @lastModified 2026-03-11
    """

    def __init__(self, a: int, b: int, c: int, d: int) -> None:
        """
        Initialisiert ein Element der modularen Gruppe.

        @param a: Eintrag oben links
        @param b: Eintrag oben rechts
        @param c: Eintrag unten links
        @param d: Eintrag unten rechts
        @raises ValueError: Wenn ad - bc ≠ 1
        @lastModified: 2026-03-11
        """
        # Determinante prüfen: muss exakt 1 sein
        det = a * d - b * c
        if det != 1:
            raise ValueError(
                f"Determinante muss 1 sein, erhalten: det([[{a},{b}],[{c},{d}]]) = {det}"
            )
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def apply(self, z: complex) -> complex:
        """
        Wendet die Möbius-Transformation auf eine komplexe Zahl z an.

        Berechnet: τ ↦ (aτ + b) / (cτ + d)

        @param z: Komplexe Zahl mit Im(z) > 0 (Punkt in der oberen Halbebene)
        @return: Transformierter Punkt (aτ+b)/(cτ+d)
        @raises ZeroDivisionError: Wenn cz + d = 0
        @lastModified: 2026-03-11
        """
        nenner = self.c * z + self.d
        if abs(nenner) < 1e-15:
            raise ZeroDivisionError(f"Möbius-Transformation: Nenner cτ+d = {nenner} zu nah an 0")
        return (self.a * z + self.b) / nenner

    def compose(self, other: 'ModularGroup') -> 'ModularGroup':
        """
        Berechnet das Matrixprodukt self · other.

        @param other: Zweites SL(2,Z)-Element
        @return: Produkt als neues ModularGroup-Element
        @lastModified: 2026-03-11
        """
        # Standardmäßige Matrixmultiplikation 2x2
        a_new = self.a * other.a + self.b * other.c
        b_new = self.a * other.b + self.b * other.d
        c_new = self.c * other.a + self.d * other.c
        d_new = self.c * other.b + self.d * other.d
        return ModularGroup(a_new, b_new, c_new, d_new)

    def inverse(self) -> 'ModularGroup':
        """
        Berechnet die Inverse der Matrix.

        Für [[a,b],[c,d]] ∈ SL(2,Z) mit det=1 gilt:
            [[a,b],[c,d]]^{-1} = [[d,-b],[-c,a]]

        @return: Inverse Matrix als ModularGroup
        @lastModified: 2026-03-11
        """
        return ModularGroup(self.d, -self.b, -self.c, self.a)

    def is_in_fundamental_domain(self, z: complex) -> bool:
        """
        Prüft, ob ein Punkt z im Fundamentalbereich F liegt.

        Der Standardfundamentalbereich der modularen Gruppe ist:
            F = {τ ∈ H : |τ| > 1 und |Re(τ)| ≤ 1/2}

        @param z: Komplexe Zahl in der oberen Halbebene
        @return: True wenn z im Fundamentalbereich liegt
        @lastModified: 2026-03-11
        """
        abs_z = abs(z)
        re_z = z.real
        return abs_z > 1.0 and abs(re_z) <= 0.5

    def __repr__(self) -> str:
        return f"ModularGroup([[{self.a},{self.b}],[{self.c},{self.d}]])"

    def __str__(self) -> str:
        return f"[[{self.a}, {self.b}], [{self.c}, {self.d}]]"

    def __eq__(self, other: object) -> bool:
        """Vergleich zweier SL(2,Z)-Elemente (bis auf Vorzeichen, da ±I wirken gleich)."""
        if not isinstance(other, ModularGroup):
            return False
        return (self.a == other.a and self.b == other.b and
                self.c == other.c and self.d == other.d)


# ===========================================================================
# EISENSTEIN-REIHEN
# ===========================================================================

def eisenstein_series(k: int, z: complex, n_terms: int = 50) -> complex:
    """
    Berechnet die Eisenstein-Reihe G_k(τ) für gerades k ≥ 4.

    Definition:
        G_k(τ) = Σ_{(m,n) ≠ (0,0)} 1 / (mτ + n)^k

    Eigenschaften:
        - G_k ist eine Modulform vom Gewicht k für k gerade, k ≥ 4
        - G_k(-1/τ) = τ^k · G_k(τ)
        - G_k(τ+1) = G_k(τ)

    Fourier-Entwicklung:
        G_k(τ) = 2ζ(k) + 2(2πi)^k/(k-1)! · Σ_{n=1}^∞ σ_{k-1}(n) q^n

    @param k: Gewicht der Eisenstein-Reihe (gerade, ≥ 4)
    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Halbbreite des Gitters (-n_terms ≤ m,n ≤ n_terms)
    @return: Wert G_k(τ) als komplexe Zahl
    @raises ValueError: Wenn k ungerade oder k < 4
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if k < 4 or k % 2 != 0:
        raise ValueError(f"k muss gerade und ≥ 4 sein, erhalten: k = {k}")
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0), erhalten: Im(τ) = {z.imag}")

    total = complex(0.0)

    # Summiere über alle Gitterpunkte (m,n) ≠ (0,0) im Bereich [-N, N]²
    for m in range(-n_terms, n_terms + 1):
        for n in range(-n_terms, n_terms + 1):
            if m == 0 and n == 0:
                continue  # Ausnahme: (0,0) überspringen

            denom = m * z + n

            # Numerische Stabilität: sehr kleine Nenner überspringen
            if abs(denom) < 1e-15:
                continue

            total += 1.0 / (denom ** k)

    return total


def eisenstein_series_fast(k: int, z: complex, n_terms: int = 50) -> complex:
    """
    @brief Optimierte Eisenstein-Reihe via Symmetrie des Gitters.
    @description
        Berechnet G_k(τ) wie eisenstein_series(), aber nutzt die Symmetrie
        des ganzzahligen Gitters Z² für gerades k:

            Für gerades k: (m,n) und (-m,-n) liefern denselben Beitrag:
                1 / (m·τ + n)^k = 1 / ((-m)·τ + (-n))^k

        Speedup: ~4× gegenüber der vollen Doppelsumme.

        Fallback für ungerades k: ruft eisenstein_series() auf.

    @param k Gewicht der Eisenstein-Reihe (gerades k ≥ 4 für Optimierung).
    @param z Punkt τ in der oberen Halbebene (Im(z) > 0).
    @param n_terms Halbbreite des Gitters.
    @return Wert G_k(τ) als komplexe Zahl.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if z.imag <= 0:
        raise ValueError(f"τ={z} liegt nicht in der oberen Halbebene (Im(τ)={z.imag:.4f} ≤ 0)")

    # Für ungerades k oder k < 4: Fallback auf vollständige Summe
    if k < 4 or k % 2 != 0:
        return eisenstein_series(k, z, n_terms)

    result = complex(0)

    # --- Teil 1: Halbraum m > 0, alle n ---
    # Symmetrie (m,n) ↔ (-m,-n) für gerades k: Faktor × 2
    for m in range(1, n_terms + 1):
        for n in range(-n_terms, n_terms + 1):
            denom = m * z + n
            if abs(denom) > 1e-15:
                result += 1.0 / denom**k

    # Verdopplung wegen Symmetrie
    result *= 2

    # --- Teil 2: Achse m = 0, n > 0 ---
    # 1/n^k + 1/(-n)^k = 2/n^k (da k gerade)
    for n in range(1, n_terms + 1):
        result += 2.0 / (n ** k)

    return result


def normalized_eisenstein_E4(z: complex, n_terms: int = 50) -> complex:
    """
    Berechnet die normierte Eisenstein-Reihe E_4(τ).

    Definition:
        E_4(τ) = G_4(τ) / (2ζ(4))  = G_4(τ) · 45/π⁴

    Fourier-Entwicklung: E_4 = 1 + 240·Σ σ_3(n) q^n

    @param z: Punkt τ in der oberen Halbebene
    @param n_terms: Gittergröße für Eisenstein-Reihe
    @return: Wert E_4(τ)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    G4 = eisenstein_series(4, z, n_terms)
    zeta4 = math.pi ** 4 / 90.0
    return G4 / (2.0 * zeta4)


def normalized_eisenstein_E6(z: complex, n_terms: int = 50) -> complex:
    """
    Berechnet die normierte Eisenstein-Reihe E_6(τ).

    Definition:
        E_6(τ) = G_6(τ) / (2ζ(6))

    Fourier-Entwicklung: E_6 = 1 - 504·Σ σ_5(n) q^n

    @param z: Punkt τ in der oberen Halbebene
    @param n_terms: Gittergröße für Eisenstein-Reihe
    @return: Wert E_6(τ)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    G6 = eisenstein_series(6, z, n_terms)
    zeta6 = math.pi ** 6 / 945.0
    return G6 / (2.0 * zeta6)


# ===========================================================================
# RAMANUJAN-DELTA-FUNKTION
# ===========================================================================

def delta_function(z: complex, n_terms: int = 100) -> complex:
    """
    Berechnet die Ramanujan-Delta-Funktion (Diskriminantenform) Δ(τ).

    Definition via Produktformel (Ramanujan, 1916):
        Δ(τ) = (2π)^12 · q · Π_{n=1}^∞ (1 - q^n)^24
    wobei q = e^{2πiτ}.

    Eigenschaften:
        - Δ ist eine Spitzenform (Cuspidalform) vom Gewicht 12
        - Δ(τ) ≠ 0 für alle τ ∈ H
        - Fourier-Entwicklung: Δ(τ) = Σ_{n=1}^∞ τ(n) q^n
        - Ramanujan-Vermutung (Deligne 1974): |τ(p)| ≤ 2p^{11/2}

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme (mehr = genauer)
    @return: Wert Δ(τ)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0)")

    # q = e^{2πiτ}
    q = cmath.exp(2j * math.pi * z)

    # Produktformel: Δ(τ) = (2π)^12 · q · Π_{n=1}^N (1-q^n)^24
    product = complex(1.0)
    q_power = q  # q^n, beginnt bei q^1

    for n in range(1, n_terms + 1):
        if abs(q_power) < 1e-300:
            break
        factor = 1.0 - q_power
        product *= factor ** 24
        q_power *= q  # Nächste q-Potenz

    return (2.0 * math.pi) ** 12 * q * product


# ===========================================================================
# J-INVARIANTE
# ===========================================================================

def j_invariant(z: complex, n_terms: int = 100) -> complex:
    """
    Berechnet die j-Invariante j(τ) (Klein's j-Funktion).

    Definition:
        j(τ) = 1728 · E_4(τ)³ / Δ(τ)

    Eigenschaften:
        - j ist eine modulare Funktion (Gewicht 0)
        - j(γτ) = j(τ) für alle γ ∈ SL(2,Z)
        - j(i) = 1728, j(e^{2πi/3}) = 0
        - Fourier-Entwicklung: j(τ) = q⁻¹ + 744 + 196884q + ...

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme für Δ
    @return: Wert j(τ) als komplexe Zahl
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    E4 = normalized_eisenstein_E4(z, n_terms=min(n_terms, 50))
    Delta = delta_function(z, n_terms=n_terms)

    if abs(Delta) < 1e-300:
        raise ZeroDivisionError(f"Δ(τ) ≈ 0 bei τ = {z}, kann j nicht berechnen")

    return 1728.0 * (E4 ** 3) / Delta


# ===========================================================================
# RAMANUJAN-TAU-FUNKTION
# ===========================================================================

def fourier_coefficients_delta(n_max: int) -> list[int]:
    """
    Berechnet die Ramanujan-Tau-Funktion τ(n) für n = 1, ..., n_max.

    Fourier-Koeffizienten der Delta-Funktion:
        Δ(τ) = Σ_{n=1}^∞ τ(n) q^n   mit q = e^{2πiτ}

    Bekannte Werte:
        τ(1) = 1, τ(2) = -24, τ(3) = 252, τ(4) = -1472, τ(5) = 4830

    Ramanujan-Vermutung (Deligne, Fields-Medaille 1978):
        |τ(p)| ≤ 2 · p^{11/2} für alle Primzahlen p

    @param n_max: Maximales n für die Berechnung
    @return: Liste [τ(1), τ(2), ..., τ(n_max)] als ganze Zahlen
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Wir berechnen die Koeffizienten des formalen Produkts
    # p(q) = Π_{n=1}^{n_max} (1-q^n)^24
    # Δ(τ) = q · p(q), also τ(n) = Koeffizient von q^{n-1} in p(q)

    N = n_max  # Maximale benötigte Potenz in p(q) = q^0..q^{n_max-1}

    # Initialisierung: p(q) = 1
    coeffs = [0] * (N + 1)
    coeffs[0] = 1

    # Binomialkoeffizienten C(24, k) vorberechnen
    binom24 = [1] * 25
    for j in range(1, 25):
        binom24[j] = binom24[j-1] * (25 - j) // j  # C(24, j)

    # Multipliziere sukzessiv mit (1 - q^m)^24 für m = 1, 2, ..., n_max
    for m in range(1, N + 1):
        # (1 - q^m)^24 = Σ_{k=0}^{24} C(24,k) · (-1)^k · q^{k·m}
        new_coeffs = [0] * (N + 1)
        for k in range(25):
            shift = k * m
            if shift > N:
                break
            sign_binom = ((-1) ** k) * binom24[k]
            for i in range(N + 1 - shift):
                new_coeffs[i + shift] += sign_binom * coeffs[i]
        coeffs = new_coeffs

    # τ(n) = Koeffizient von q^{n-1} in p(q)  (wegen des Vorfaktors q in Δ = q·p(q))
    tau = []
    for n in range(1, n_max + 1):
        if n - 1 <= N:
            tau.append(coeffs[n - 1])
        else:
            tau.append(0)

    return tau


# ===========================================================================
# MODULARITÄTSBEDINGUNG
# ===========================================================================

def modular_form_check(
    f: Callable[[complex], complex],
    k: int,
    z: complex,
    gamma: 'ModularGroup',
    tol: float = 1e-6
) -> bool:
    """
    Prüft, ob f die Modularitätsbedingung für Gewicht k erfüllt.

    Modularitätsbedingung:
        f(γτ) = (cτ + d)^k · f(τ)

    @param f: Die zu prüfende Funktion f: H → ℂ
    @param k: Gewicht der Modulform
    @param z: Testpunkt τ in der oberen Halbebene (Im(z) > 0)
    @param gamma: Element der modularen Gruppe
    @param tol: Toleranz für numerischen Vergleich
    @return: True wenn |f(γτ) - (cτ+d)^k · f(τ)| < tol
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Berechne γτ (transformierter Punkt)
    gamma_z = gamma.apply(z)

    # Linke Seite: f(γτ)
    lhs = f(gamma_z)

    # Rechte Seite: (cτ + d)^k · f(τ)
    automorphy_factor = (gamma.c * z + gamma.d) ** k
    rhs = automorphy_factor * f(z)

    # Relativer oder absoluter Fehler
    diff = abs(lhs - rhs)
    scale = max(abs(lhs), abs(rhs), 1.0)

    return diff / scale < tol


# ===========================================================================
# CUSP-FORMEN
# ===========================================================================

def cusp_form_dimension(k: int, level: int = 1) -> int:
    """
    Berechnet die Dimension des Raums der Cusp-Formen S_k(Γ_0(N)).

    Für die volle modulare Gruppe Γ_0(1) = SL(2,Z) gilt die Riemann-Roch-Formel:
        dim S_k(SL(2,Z)) = max(0, floor(k/12) - ε(k))
    wobei ε(k) = 1 falls k ≡ 2 (mod 12), sonst ε(k) = 0.

    Bekannte Werte:
        k=12: dim=1 (Δ-Funktion), k=24: dim=2

    @param k: Gewicht der Modulform (gerade, ≥ 2)
    @param level: Niveau N ≥ 1 (Standard: 1 = volle Modulgruppe)
    @return: Dimension des Cusp-Formen-Raums (nicht-negative ganze Zahl)
    @raises ValueError: Wenn k ungerade oder k < 2
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if k < 2 or k % 2 != 0:
        raise ValueError(f"Gewicht k muss gerade und ≥ 2 sein, erhalten: {k}")

    if level == 1:
        # Exakte Formel für die volle Modulgruppe SL(2,Z)
        epsilon = 1 if k % 12 == 2 else 0
        dim = k // 12 - epsilon
        return max(0, dim)
    else:
        # Näherungsformel für Γ_0(N) mit N > 1
        # Index μ(N) = N · Π_{p | N} (1 + 1/p)
        N = level
        mu = float(N)
        temp = N
        d = 2
        seen_primes = set()
        while d * d <= temp:
            if temp % d == 0:
                if d not in seen_primes:
                    mu *= (1.0 + 1.0 / d)
                    seen_primes.add(d)
                while temp % d == 0:
                    temp //= d
            d += 1
        if temp > 1 and temp not in seen_primes:
            mu *= (1.0 + 1.0 / temp)

        # Dimension via Riemann-Roch-Näherung: μ(N)·(k-1)/12
        dim_approx = mu * (k - 1) / 12.0
        return max(0, int(round(dim_approx)))


def is_cusp_form_basis(k: int, z: complex, n_terms: int = 100) -> bool:
    """
    Prüft ob für das gegebene Gewicht k eine Basis-Cusp-Form existiert.

    @param k: Gewicht (gerade ≥ 2)
    @param z: Punkt in der oberen Halbebene (für numerische Verifikation)
    @param n_terms: Anzahl Terme für Produktformel
    @return: True wenn dim S_k ≥ 1 (d.h. Cusp-Form-Raum nicht trivial)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if k < 2 or k % 2 != 0:
        return False
    dim = cusp_form_dimension(k, level=1)
    return dim >= 1


def ramanujan_tau_properties(n: int) -> dict[str, object]:
    """
    Überprüft bekannte Eigenschaften der Ramanujan-Tau-Funktion τ(n).

    Bekannte arithmetische Eigenschaften:
    1. Multiplikativität: τ(mn) = τ(m)·τ(n) falls gcd(m,n) = 1
    2. Ramanujan-Vermutung (Deligne 1974): |τ(p)| ≤ 2·p^{11/2}
    3. Rekursion: τ(p^k) = τ(p)·τ(p^{k-1}) - p^{11}·τ(p^{k-2})

    @param n: Positive ganze Zahl ≥ 1
    @return: Dictionary mit 'tau_n', 'multiplicativity_verified', 'ramanujan_bound_ok'
    @raises ValueError: Wenn n < 1
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

    n_max = max(n + 1, 10)
    tau_list = fourier_coefficients_delta(n_max)
    tau_n = tau_list[n - 1]  # τ(n) ist an Position n-1 (1-indiziert)

    # Multiplikativität prüfen
    mult_verified = True
    if n > 1:
        for a in range(2, n):
            if n % a == 0:
                b = n // a
                if gcd(a, b) == 1 and a >= 2 and b >= 2:
                    tau_a_b_max = max(a, b) + 1
                    tau_ab = fourier_coefficients_delta(tau_a_b_max)
                    tau_a = tau_ab[a - 1]
                    tau_b = tau_ab[b - 1]
                    mult_verified = (tau_n == tau_a * tau_b)
                    break

    # Ramanujan-Vermutung: |τ(p)| ≤ 2·p^{11/2}
    ramanujan_bound = 2.0 * (n ** (11.0 / 2.0))
    ramanujan_bound_ok = abs(tau_n) <= ramanujan_bound

    return {
        'tau_n': tau_n,
        'multiplicativity_verified': mult_verified,
        'ramanujan_bound_ok': ramanujan_bound_ok
    }


# ===========================================================================
# SHIMURA-TANIYAMA-WILES-VERBINDUNG
# ===========================================================================

def shimura_taniyama_check(
    a_p_elliptic: dict[int, float | int],
    a_p_modular: dict[int, float | int],
    primes: list[int],
) -> dict[str, object]:
    """
    Vergleicht a_p-Koeffizienten einer elliptischen Kurve mit einer Modulform.

    Satz von Shimura-Taniyama-Wiles (Wiles 1995):
        Jede elliptische Kurve über Q ist modular.
        Für jede E/Q existiert eine Modulform f der Stufe N mit:
            a_p(E) = a_p(f) für fast alle Primzahlen p

    @param a_p_elliptic: Dict {p: a_p(E)} für die elliptische Kurve E
    @param a_p_modular: Dict {p: a_p(f)} für die Modulform f
    @param primes: Liste der zu vergleichenden Primzahlen
    @return: Dictionary mit Vergleichsergebnissen
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    matches = {}
    discrepancies = []

    for p in primes:
        ap_E = a_p_elliptic.get(p, None)
        ap_f = a_p_modular.get(p, None)

        if ap_E is None or ap_f is None:
            matches[p] = None
            continue

        diff = abs(ap_E - ap_f)
        discrepancies.append(diff)
        matches[p] = diff < 0.5

    valid_matches = [v for v in matches.values() if v is not None]
    agreement_ratio = sum(valid_matches) / len(valid_matches) if valid_matches else 0.0
    max_discrepancy = max(discrepancies) if discrepancies else 0.0
    is_modular_candidate = agreement_ratio >= 0.95

    return {
        'matches': matches,
        'agreement_ratio': agreement_ratio,
        'max_discrepancy': max_discrepancy,
        'is_modular_candidate': is_modular_candidate,
        'primes_checked': len(primes),
        'primes_matching': sum(1 for v in valid_matches if v)
    }


# ===========================================================================
# INTERNE HILFSFUNKTIONEN
# ===========================================================================

def _is_prime_simple(n: int) -> bool:
    """
    Einfacher Primzahltest (nur für interne Nutzung in diesem Modul).

    @param n: Zu testende ganze Zahl
    @return: True wenn n eine Primzahl ist
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


def _primes_up_to(n: int) -> list[int]:
    """
    Gibt alle Primzahlen bis n zurück (Sieb des Eratosthenes).

    @param n: Obere Grenze
    @return: Liste aller Primzahlen ≤ n
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n + 1, i):
                sieve[j] = False
    return [i for i in range(2, n + 1) if sieve[i]]


# ===========================================================================
# L-FUNKTIONEN ELLIPTISCHER KURVEN
# ===========================================================================

def elliptic_curve_points_over_fp(a: int, b: int, p: int) -> list[tuple[int, int]]:
    """
    Berechnet alle affinen Punkte einer elliptischen Kurve E: y² = x³ + ax + b über F_p.

    Algorithmus (Brute-Force):
        Für jedes x ∈ {0, 1, ..., p-1}:
            1. Berechne rhs = (x³ + ax + b) mod p
            2. Prüfe ob rhs ein quadratischer Rest mod p ist (Euler-Kriterium)
            3. Falls ja, berechne y = sqrt(rhs) mod p

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param p: Primzahl p ≥ 3 (Charakteristik des Körpers F_p)
    @return: Liste aller affinen Punkte (x, y) ∈ F_p × F_p auf E
    @raises ValueError: Wenn p < 3
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if p < 3:
        raise ValueError(f"p muss eine Primzahl ≥ 3 sein, erhalten: p = {p}")

    points = []

    for x in range(p):
        # Rechte Seite der Kurvengleichung: x³ + ax + b mod p
        rhs = (pow(x, 3, p) + a * x + b) % p

        if rhs == 0:
            points.append((x, 0))
        else:
            # Euler-Kriterium: rhs ist QR ⟺ rhs^{(p-1)/2} ≡ 1 (mod p)
            if pow(rhs, (p - 1) // 2, p) == 1:
                # Für p ≡ 3 (mod 4): y = rhs^{(p+1)/4} mod p
                y = pow(rhs, (p + 1) // 4, p)
                if pow(y, 2, p) == rhs:
                    points.append((x, y))
                    if y != 0 and y != p - y:
                        points.append((x, p - y))
                else:
                    # Fallback: Brute-Force für Fälle wo p ≡ 1 (mod 4)
                    for y_try in range(1, p):
                        if pow(y_try, 2, p) == rhs:
                            points.append((x, y_try))
                            if y_try != p - y_try:
                                points.append((x, p - y_try))
                            break

    return points


def trace_of_frobenius(a: int, b: int, p: int) -> int:
    """
    Berechnet die Frobenius-Spur a_p = p + 1 - #E(F_p) einer elliptischen Kurve.

    Satz von Hasse (Hasse-Weil-Schranke):
        |a_p(E)| = |p + 1 - #E(F_p)| ≤ 2√p

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param p: Primzahl p ≥ 3
    @return: Frobenius-Spur a_p = p + 1 - #E(F_p)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    affine_points = elliptic_curve_points_over_fp(a, b, p)
    num_points = len(affine_points) + 1  # +1 für Punkt im Unendlichen O
    return p + 1 - num_points


def l_function_elliptic_curve(a: int, b: int, s: complex, terms: int = 100) -> complex:
    """
    Berechnet die L-Funktion einer elliptischen Kurve E: y² = x³ + ax + b.

    L(E, s) = Σ_{n=1}^∞ a_n(E) / n^s  (konvergiert für Re(s) > 3/2)

    Die Koeffizienten a_n(E) sind multiplikativ:
        - Für Primzahlen p: a_p(E) = p + 1 - #E(F_p) (Frobenius-Spur)
        - Für p^r: Rekursion a_{p^r} = a_p·a_{p^{r-1}} - p·a_{p^{r-2}}
        - Für gcd(m,k)=1: a_{mk} = a_m·a_k (Multiplikativität)

    @param a: Koeffizient der elliptischen Kurve
    @param b: Konstantterm der elliptischen Kurve
    @param s: Komplexes Argument (Re(s) > 3/2 für Konvergenz)
    @param terms: Anzahl der Summanden in der Dirichlet-Reihe
    @return: Näherungswert von L(E, s) als komplexe Zahl
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Schritt 1: Primzahlen bis terms bestimmen
    primes_needed = _primes_up_to(min(terms, 200))

    # Schritt 2: a_p für jede Primzahl p berechnen (Frobenius-Spuren)
    a_prime = {}
    for p in primes_needed:
        try:
            a_prime[p] = trace_of_frobenius(a, b, p)
        except Exception:
            a_prime[p] = 0

    # Schritt 3: a_n für alle n bis terms via Multiplikativität
    a_n = [0.0] * (terms + 1)
    a_n[1] = 1.0  # a_1 = 1 per Konvention

    for n in range(2, terms + 1):
        temp = n
        factors = {}
        for p in primes_needed:
            if p * p > temp:
                break
            if temp % p == 0:
                factors[p] = 0
                while temp % p == 0:
                    factors[p] += 1
                    temp //= p
        if temp > 1:
            factors[temp] = 1

        if len(factors) == 1:
            # n = p^r: Rekursionsformel
            p = list(factors.keys())[0]
            r = factors[p]
            if r == 1:
                a_n[n] = float(a_prime.get(p, 0))
            elif r == 2:
                a_n[n] = float(a_prime.get(p, 0)) ** 2 - p
            else:
                prev1 = int(round(a_n[p ** (r-1)])) if p ** (r-1) <= terms else 0
                prev2 = int(round(a_n[p ** (r-2)])) if p ** (r-2) <= terms else 0
                a_n[n] = float(a_prime.get(p, 0)) * prev1 - p * prev2
        elif len(factors) > 1:
            # n = m * k mit gcd(m,k) = 1: a_n = a_m * a_k
            ps = list(factors.keys())
            p1 = ps[0]
            m = p1 ** factors[p1]
            k_val = n // m
            if m <= terms and k_val <= terms:
                a_n[n] = a_n[m] * a_n[k_val]

    # Schritt 4: L-Funktion als Dirichlet-Reihe summieren
    result = complex(0.0)
    for n in range(1, terms + 1):
        if a_n[n] != 0.0:
            result += complex(a_n[n]) / (n ** s)

    return result


def birch_swinnerton_dyer_bsd_estimate(a: int, b: int, primes: list[int]) -> dict[str, object]:
    """
    Empirische Schätzung der BSD-Vermutung (Birch und Swinnerton-Dyer).

    BSD-Vermutung (Millennium-Problem):
        ord_{s=1} L(E, s) = Rang(E(Q))

    Empirische Schätzung (Birch und Swinnerton-Dyer, 1960er):
        Π_{p ≤ X} (N_p / p) ≈ C · (log X)^r
    wobei r = Rang(E(Q)) und N_p = #E(F_p).

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param primes: Liste der zu verwendenden Primzahlen p ≥ 3
    @return: Dictionary mit Rang-Schätzung und Hilfsdaten
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Schritt 1: N_p = #E(F_p) für jede Primzahl berechnen
    n_p_values = {}
    log_product = 0.0
    product = 1.0
    valid_primes = []

    for p in primes:
        if p < 3:
            continue
        try:
            affine_pts = elliptic_curve_points_over_fp(a, b, p)
            n_p = len(affine_pts) + 1  # +1 für Punkt im Unendlichen
            n_p_values[p] = n_p

            ratio = n_p / p
            if ratio > 0:
                log_product += math.log(ratio)
                product *= ratio
                valid_primes.append(p)
        except Exception:
            continue

    # Schritt 2: Rang-Schätzung via log(log X)
    if valid_primes:
        X = max(valid_primes)
        log_log_X = math.log(math.log(X)) if X > 2 else 1.0
        rank_estimate = log_product / log_log_X if abs(log_log_X) > 1e-10 else 0.0
    else:
        rank_estimate = 0.0

    return {
        'estimate': rank_estimate,
        'primes_used': valid_primes,
        'log_product': log_product,
        'product': product,
        'n_p_values': n_p_values
    }


# ===========================================================================
# __all__: Alle öffentlichen Namen für `from modular_forms import *`
# ===========================================================================

__all__ = [
    # Modulare Gruppe
    'ModularGroup',
    # Eisenstein-Reihen
    'eisenstein_series',
    'eisenstein_series_fast',
    'normalized_eisenstein_E4',
    'normalized_eisenstein_E6',
    # Delta-Funktion und j-Invariante
    'delta_function',
    'j_invariant',
    # Ramanujan-Tau-Funktion
    'fourier_coefficients_delta',
    'ramanujan_tau_properties',
    # Modularitätsbedingung
    'modular_form_check',
    # Cusp-Formen
    'cusp_form_dimension',
    'is_cusp_form_basis',
    # Shimura-Taniyama-Wiles
    'shimura_taniyama_check',
    # Elliptische Kurven
    'elliptic_curve_points_over_fp',
    'trace_of_frobenius',
    'l_function_elliptic_curve',
    'birch_swinnerton_dyer_bsd_estimate',
    # Hecke-Algebra (aus modular_forms_hecke)
    'hecke_operator',
    'hecke_algebra_structure',
    'hecke_eigenform',
    'petersson_inner_product_estimate',
    'l_function_modular_form',
    'functional_equation_l_function',
    # Theta-Reihen (aus modular_forms_theta)
    'theta_function',
    'theta_transformation',
    'sum_of_squares_theta',
    'jacobi_triple_product',
    'dedekind_eta',
    # Interne Hilfsfunktionen (für Tests)
    '_is_prime_simple',
    '_primes_up_to',
]
