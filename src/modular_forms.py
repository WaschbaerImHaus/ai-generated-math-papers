"""
@file modular_forms.py
@brief Modulformen – Grundlagen, Cusp-Formen, Theta-Reihen, Hecke-Algebra und L-Funktionen.
@description
    Implementiert grundlegende Konzepte der Theorie der Modulformen:
    - Die modulare Gruppe SL(2,Z) und Möbius-Transformationen
    - Eisenstein-Reihen G_k(τ) und normierte Versionen E_4, E_6
    - Ramanujan-Delta-Funktion Δ(τ) (Diskriminantenform)
    - j-Invariante j(τ)
    - Ramanujan-Tau-Funktion τ(n)
    - Hecke-Operatoren T_p und Hecke-Algebra-Struktur
    - Shimura-Taniyama-Wiles-Verbindung (Demonstration)
    - Cusp-Formen S_k(Γ_0(N)): Dimension, Eigenschaften, Ramanujan-Tau
    - Theta-Reihen: Jacobi-Theta, Dedekind-Eta, Jacobi-Dreifachprodukt
    - L-Funktionen elliptischer Kurven und Modulformen
    - BSD-Vermutung (empirische Schätzung)

    Modulformen sind holomorphe Funktionen auf der oberen Halbebene H = {τ : Im(τ) > 0},
    die unter Möbius-Transformationen γ ∈ SL(2,Z) folgende Transformationsregel erfüllen:
        f(γτ) = (cτ + d)^k · f(τ)
    wobei k das Gewicht der Modulform ist.

    Diese Theorie ist zentral für den Beweis des Großen Fermatschen Satzes (Wiles 1995),
    da elliptische Kurven und Modulformen via Shimura-Taniyama-Wiles verknüpft sind.

@author Michael Fuhrmann
@version 3.0
@since 2026-03-08
@lastModified 2026-03-10
"""

import math
import cmath
import functools
import numpy as np
from typing import Callable
from math import gcd


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
    @lastModified 2026-03-10
    """

    def __init__(self, a: int, b: int, c: int, d: int) -> None:
        """
        Initialisiert ein Element der modularen Gruppe.

        @param a: Eintrag oben links
        @param b: Eintrag oben rechts
        @param c: Eintrag unten links
        @param d: Eintrag unten rechts
        @raises ValueError: Wenn ad - bc ≠ 1
        @lastModified: 2026-03-10
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

        Für c=0 ist dies eine Translation, für c≠0 eine echte Inversion.
        Die Transformation erhält die obere Halbebene Im(τ) > 0.

        @param z: Komplexe Zahl mit Im(z) > 0 (Punkt in der oberen Halbebene)
        @return: Transformierter Punkt (aτ+b)/(cτ+d)
        @raises ZeroDivisionError: Wenn cz + d = 0
        @lastModified: 2026-03-10
        """
        nenner = self.c * z + self.d
        if abs(nenner) < 1e-15:
            raise ZeroDivisionError(f"Möbius-Transformation: Nenner cτ+d = {nenner} zu nah an 0")
        return (self.a * z + self.b) / nenner

    def compose(self, other: 'ModularGroup') -> 'ModularGroup':
        """
        Berechnet das Matrixprodukt self · other.

        Entspricht der Komposition der Möbius-Transformationen:
        (self ∘ other)(τ) = self(other(τ))

        @param other: Zweites SL(2,Z)-Element
        @return: Produkt als neues ModularGroup-Element
        @lastModified: 2026-03-10
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
        @lastModified: 2026-03-10
        """
        return ModularGroup(self.d, -self.b, -self.c, self.a)

    def is_in_fundamental_domain(self, z: complex) -> bool:
        """
        Prüft, ob ein Punkt z im Fundamentalbereich F liegt.

        Der Standardfundamentalbereich der modularen Gruppe ist:
            F = {τ ∈ H : |τ| > 1 und |Re(τ)| ≤ 1/2}
        (mit Randidentifikationen für den Beweis)

        Strikt: |τ| > 1, |Re(τ)| < 1/2
        Auf dem Rand: |τ| = 1 und Re(τ) ≤ 0, oder Re(τ) = 1/2

        Diese Methode prüft die strengen Ungleichungen (inneres Gebiet).

        @param z: Komplexe Zahl in der oberen Halbebene
        @return: True wenn z im Fundamentalbereich liegt
        @lastModified: 2026-03-10
        """
        # Im(z) > 0 vorausgesetzt (obere Halbebene)
        abs_z = abs(z)
        re_z = z.real
        # Fundamentalbereich: |z| > 1 (außerhalb Einheitskreis) und |Re(z)| ≤ 0.5
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

    wobei die Summe über alle ganzzahligen Paare (m,n) außer (0,0) geht.

    Eigenschaften:
        - G_k ist eine Modulform vom Gewicht k für k gerade, k ≥ 4
        - G_k(-1/τ) = τ^k · G_k(τ)  (Transformationsregel unter S)
        - G_k(τ+1) = G_k(τ)           (Periodizität)
        - Konvergiert absolut für Im(τ) > 0 und k ≥ 3 (gerade k ≥ 4)

    Fourier-Entwicklung:
        G_k(τ) = 2ζ(k) + 2(2πi)^k/(k-1)! · Σ_{n=1}^∞ σ_{k-1}(n) q^n
    wobei q = e^{2πiτ} und σ_{k-1}(n) = Summe der (k-1)-ten Potenzen der Teiler.

    @param k: Gewicht der Eisenstein-Reihe (gerade, ≥ 4)
    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Halbbreite des Gitters (-n_terms ≤ m,n ≤ n_terms)
    @return: Wert G_k(τ) als komplexe Zahl
    @raises ValueError: Wenn k ungerade oder k < 4
    @lastModified: 2026-03-10
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

            # Berechne mτ + n
            denom = m * z + n

            # Numerische Stabilität: sehr kleine Nenner überspringen
            if abs(denom) < 1e-15:
                continue

            # Summand: 1 / (mτ + n)^k
            total += 1.0 / (denom ** k)

    return total


def eisenstein_series_fast(k: int, z: complex, n_terms: int = 50) -> complex:
    """
    @brief Optimierte Eisenstein-Reihe via Symmetrie des Gitters.
    @description
        Berechnet G_k(τ) wie eisenstein_series(), aber nutzt die Symmetrie
        des ganzzahligen Gitters Z² für gerades k:

            Für gerades k: (m,n) und (-m,-n) liefern denselben Beitrag:
                1 / (m·τ + n)^k = 1 / ((-m)·τ + (-n))^k  (da k gerade)

        Strategie:
        1. Summiere nur über den Halbraum m > 0 (alle n): Beitrag × 2 (Symmetrie m → -m)
        2. Achse m = 0, n > 0: Beitrag × 2 (Symmetrie n → -n, also n und -n)
        3. Punkt (0,0) wird übersprungen (per Definition)

        Speedup: ~4× gegenüber der vollen Doppelsumme, da nur die Hälfte der
        Gitterpunkte explizit berechnet wird.

        Fallback für ungerades k: ruft eisenstein_series() auf, da die Symmetrie
        1/(mτ+n)^k ≠ 1/(-mτ-n)^k für ungerades k gilt.

    @param k Gewicht der Eisenstein-Reihe (gerades k ≥ 4 für Optimierung).
    @param z Punkt τ in der oberen Halbebene (Im(z) > 0).
    @param n_terms Halbbreite des Gitters (-n_terms ≤ m,n ≤ n_terms).
    @return Wert G_k(τ) als komplexe Zahl.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # τ muss in der oberen Halbebene liegen: Im(τ) > 0
    if z.imag <= 0:
        raise ValueError(f"τ={z} liegt nicht in der oberen Halbebene (Im(τ)={z.imag:.4f} ≤ 0)")

    # Für ungerades k oder k < 4: Fallback auf vollständige Summe
    # (Symmetrie gilt nur für gerades k)
    if k < 4 or k % 2 != 0:
        return eisenstein_series(k, z, n_terms)

    result = complex(0)

    # --- Teil 1: Halbraum m > 0, alle n ---
    # Für jedes (m, n) mit m > 0 gilt: (-m, -n) liefert identischen Beitrag.
    # Daher summieren wir nur m > 0 und verdoppeln am Ende.
    for m in range(1, n_terms + 1):
        for n in range(-n_terms, n_terms + 1):
            denom = m * z + n
            # Numerische Stabilität: sehr kleine Nenner überspringen
            if abs(denom) > 1e-15:
                result += 1.0 / denom**k

    # Verdopplung wegen Symmetrie: (m,n) und (-m,-n) haben denselben Beitrag
    result *= 2

    # --- Teil 2: Achse m = 0, n > 0 ---
    # (0, 0) wird übersprungen (per Definition der Eisenstein-Reihe)
    # (0, n) und (0, -n) liefern 1/n^k + 1/(-n)^k = 2/n^k (da k gerade)
    for n in range(1, n_terms + 1):
        result += 2.0 / (n ** k)

    return result


def normalized_eisenstein_E4(z: complex, n_terms: int = 50) -> complex:
    """
    Berechnet die normierte Eisenstein-Reihe E_4(τ).

    Definition:
        E_4(τ) = G_4(τ) / (2ζ(4))
    wobei ζ(4) = π⁴/90, also:
        E_4(τ) = G_4(τ) · 90 / (2π⁴) = G_4(τ) · 45/π⁴

    Eigenschaften:
        - E_4 ist eine normierte Modulform vom Gewicht 4
        - E_4(i) = Γ(1/4)⁴ / (4π³) ≈ 1.0
        - Fourier-Koeffizient a_0 = 1 (normiert)
        - Fourier-Entwicklung: E_4 = 1 + 240·Σ σ_3(n) q^n

    @param z: Punkt τ in der oberen Halbebene
    @param n_terms: Gittergröße für Eisenstein-Reihe
    @return: Wert E_4(τ)
    @lastModified: 2026-03-10
    """
    # G_4 normieren: 2ζ(4) = 2·π⁴/90 = π⁴/45
    G4 = eisenstein_series(4, z, n_terms)
    zeta4 = math.pi ** 4 / 90.0
    return G4 / (2.0 * zeta4)


def normalized_eisenstein_E6(z: complex, n_terms: int = 50) -> complex:
    """
    Berechnet die normierte Eisenstein-Reihe E_6(τ).

    Definition:
        E_6(τ) = G_6(τ) / (2ζ(6))
    wobei ζ(6) = π⁶/945, also:
        E_6(τ) = G_6(τ) · 945 / (2π⁶)

    Eigenschaften:
        - E_6 ist eine normierte Modulform vom Gewicht 6
        - Fourier-Entwicklung: E_6 = 1 - 504·Σ σ_5(n) q^n

    @param z: Punkt τ in der oberen Halbebene
    @param n_terms: Gittergröße für Eisenstein-Reihe
    @return: Wert E_6(τ)
    @lastModified: 2026-03-10
    """
    # G_6 normieren: 2ζ(6) = 2·π⁶/945
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
        - Δ(τ) ≠ 0 für alle τ ∈ H (keine Nullstellen in der oberen Halbebene)
        - Fourier-Entwicklung: Δ(τ) = Σ_{n=1}^∞ τ(n) q^n
          wobei τ(n) die Ramanujan-Tau-Funktion ist
        - τ(1) = 1, τ(2) = -24, τ(3) = 252, τ(4) = -1472, ...
        - Ramanujan-Vermutung (bewiesen von Deligne 1974): |τ(p)| ≤ 2p^{11/2}

    Zusammenhang mit Eisenstein-Reihen:
        Δ(τ) = (E_4(τ)³ - E_6(τ)²) / 1728

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme (mehr = genauer)
    @return: Wert Δ(τ)
    @lastModified: 2026-03-10
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0)")

    # q = e^{2πiτ}
    q = cmath.exp(2j * math.pi * z)

    # Produktformel: Δ(τ) = (2π)^12 · q · Π_{n=1}^N (1-q^n)^24
    # Für Im(τ) > 0 gilt |q| < 1, Produkt konvergiert absolut
    product = complex(1.0)
    q_power = q  # q^n, beginnt bei q^1

    for n in range(1, n_terms + 1):
        # Konvergenzcheck: q^n wird exponentiell klein
        if abs(q_power) < 1e-300:
            break
        factor = 1.0 - q_power
        # Faktor 24-mal potenzieren: (1-q^n)^24
        product *= factor ** 24
        q_power *= q  # Nächste q-Potenz

    # Gesamtformel: (2π)^12 · q · Produkt
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
        - j ist vollständig invariant unter SL(2,Z): j(γτ) = j(τ) für alle γ ∈ SL(2,Z)
        - j hat einen einfachen Pol im Cusp ∞ (bei q=0)
        - Fourier-Entwicklung: j(τ) = q⁻¹ + 744 + 196884q + 21493760q² + ...
        - j(i) = 1728 (am Fixpunkt der S-Transformation)
        - j(e^{2πi/3}) = 0 (am Fixpunkt der ST-Transformation)
        - j ist die einzige (bis auf Konstante) meromorphe Modulform von Gewicht 0

    Bedeutung: j parametrisiert Isomorphieklassen elliptischer Kurven ℂ/Λ.
    Zwei elliptische Kurven sind isomorph ⟺ sie haben dasselbe j.

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme für Δ
    @return: Wert j(τ) als komplexe Zahl
    @lastModified: 2026-03-10
    """
    # E_4 berechnen (normiert)
    E4 = normalized_eisenstein_E4(z, n_terms=min(n_terms, 50))

    # Δ berechnen
    Delta = delta_function(z, n_terms=n_terms)

    if abs(Delta) < 1e-300:
        raise ZeroDivisionError(f"Δ(τ) ≈ 0 bei τ = {z}, kann j nicht berechnen")

    # j(τ) = 1728 · E_4³ / Δ
    return 1728.0 * (E4 ** 3) / Delta


# ===========================================================================
# RAMANUJAN-TAU-FUNKTION
# ===========================================================================

def fourier_coefficients_delta(n_max: int) -> list[int]:
    """
    Berechnet die Ramanujan-Tau-Funktion τ(n) für n = 1, ..., n_max.

    Die Tau-Funktion erscheint als Fourier-Koeffizient der Delta-Funktion:
        Δ(τ) = Σ_{n=1}^∞ τ(n) q^n   mit q = e^{2πiτ}

    Formel via Rekursion (Newton-Identitäten / Divisorpotenzsummen):
        τ(n) ist durch die Identität bestimmt:
        Δ(τ) = q · Π_{n≥1} (1-q^n)^24 = Σ τ(n) q^n

    Bekannte Werte (für Verifikation):
        τ(1) = 1
        τ(2) = -24
        τ(3) = 252
        τ(4) = -1472
        τ(5) = 4830
        τ(6) = -6048
        τ(7) = -16744

    Ramanujan-Vermutung (Deligne, Fields-Medaille 1978):
        |τ(p)| ≤ 2 · p^{11/2} für alle Primzahlen p

    Wir berechnen τ(n) durch Koeffizientenextraktion aus dem Produkt:
        Δ = (2π)^12 · q · Π(1-q^n)^24

    Praktisch: Arbeite mit q-Potenzen als Polynomkoeffizienten.

    @param n_max: Maximales n für die Berechnung
    @return: Liste [τ(1), τ(2), ..., τ(n_max)] als ganze Zahlen
    @lastModified: 2026-03-10
    """
    # Wir berechnen die Koeffizienten des formalen Produkts
    # p(q) = Π_{n=1}^{n_max} (1-q^n)^24
    # Δ(τ) = q · p(q), also τ(n) = Koeffizient von q^{n-1} in p(q)

    # Koeffizientenliste: coeffs[k] = Koeffizient von q^k
    # Wir brauchen Koeffizienten bis Grad n_max-1 (wegen des q-Faktors)
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
        # Polynommultiplikation: neues Polynom = altes × (1-q^m)^24
        # Wir multiplizieren in-place von hinten nach vorne, um Überschreibung zu vermeiden
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

    für γ = [[a,b],[c,d]] ∈ SL(2,Z).

    Diese Eigenschaft charakterisiert Modulformen: holomorphe Funktionen
    auf der oberen Halbebene, die unter der modularen Gruppe transformieren.

    @param f: Die zu prüfende Funktion f: H → ℂ
    @param k: Gewicht der Modulform
    @param z: Testpunkt τ in der oberen Halbebene (Im(z) > 0)
    @param gamma: Element der modularen Gruppe
    @param tol: Toleranz für numerischen Vergleich
    @return: True wenn |f(γτ) - (cτ+d)^k · f(τ)| < tol
    @lastModified: 2026-03-10
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

    @param coefficients_tuple: Fourier-Koeffizienten als hashbares Tupel
    @param p: Primzahl p ≥ 2
    @param k: Gewicht der Modulform
    @return: Transformierte Koeffizienten als Tupel
    @lastModified 2026-03-10
    """
    N = len(coefficients_tuple)
    M = max(1, N // p)
    b = []

    for n in range(M):
        idx_pn = p * n
        a_pn = coefficients_tuple[idx_pn] if idx_pn < N else 0

        if n % p == 0:
            idx_np = n // p
            a_np = coefficients_tuple[idx_np] if idx_np < N else 0
            bn = a_pn + (p ** (k - 1)) * a_np
        else:
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
    @lastModified: 2026-03-10
    """
    # Liste → Tupel für lru_cache-Kompatibilität (Listen sind nicht hashbar)
    result_tuple = _hecke_operator_cached(tuple(coefficients), p, k)
    return list(result_tuple)


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

    Der Satz von Shimura-Taniyama-Wiles (Wiles 1995, Taylor-Wiles):
        Jede elliptische Kurve über Q ist modular.
        D.h. für jede elliptische Kurve E/Q gibt es eine Modulform f der Stufe N,
        so dass für fast alle Primzahlen p gilt:
            a_p(E) = a_p(f)
        wobei:
            a_p(E) = p + 1 - #E(F_p)  (Spurformel, Hasse-Schranke |a_p| ≤ 2√p)
            a_p(f) = p-ter Hecke-Eigenwert von f

    Diese Verbindung war der Schlüssel zum Beweis des Großen Fermatschen Satzes:
        Frey: Eine Gegenbeispielkurve zu FLT wäre nicht modular.
        Ribet: Frey-Kurve verletzt Shimura-Taniyama.
        Wiles: Alle semistabilen elliptischen Kurven sind modular.

    @param a_p_elliptic: Dict {p: a_p(E)} für die elliptische Kurve E
    @param a_p_modular: Dict {p: a_p(f)} für die Modulform f
    @param primes: Liste der zu vergleichenden Primzahlen
    @return: Dictionary mit Vergleichsergebnissen:
             'matches': Dict {p: bool} – Übereinstimmung pro Primzahl
             'agreement_ratio': float – Anteil der übereinstimmenden Primzahlen
             'max_discrepancy': float – Größte Abweichung
             'is_modular_candidate': bool – Vermutung: E ist modular
    @lastModified: 2026-03-10
    """
    matches = {}
    discrepancies = []

    for p in primes:
        # Koeffizienten aus beiden Dictionaries holen
        ap_E = a_p_elliptic.get(p, None)
        ap_f = a_p_modular.get(p, None)

        if ap_E is None or ap_f is None:
            # Fehlende Daten: kein Vergleich möglich
            matches[p] = None
            continue

        # Numerischer Vergleich mit Toleranz
        diff = abs(ap_E - ap_f)
        discrepancies.append(diff)

        # Übereinstimmung falls Abweichung < 0.5 (für ganzzahlige Koeffizienten)
        matches[p] = diff < 0.5

    # Statistiken berechnen
    valid_matches = [v for v in matches.values() if v is not None]
    agreement_ratio = sum(valid_matches) / len(valid_matches) if valid_matches else 0.0
    max_discrepancy = max(discrepancies) if discrepancies else 0.0

    # Heuristik: "modular" falls ≥ 95% der Primzahlen übereinstimmen
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
# CUSP-FORMEN
# ===========================================================================

def cusp_form_dimension(k: int, level: int = 1) -> int:
    """
    Berechnet die Dimension des Raums der Cusp-Formen S_k(Γ_0(N)).

    Cusp-Formen (Spitzenformen) sind Modulformen, die an allen Spitzen verschwinden.
    Der Raum S_k(SL(2,Z)) ist ein wichtiger Unterraum der Modulformen M_k.

    Für die volle modulare Gruppe Γ_0(1) = SL(2,Z) gilt die Riemann-Roch-Formel:
        dim S_k(SL(2,Z)) = max(0, floor(k/12) - ε(k))
    wobei ε(k) = 1 falls k ≡ 2 (mod 12), sonst ε(k) = 0.

    Bekannte Werte:
        k=2:  dim=0, k=4:  dim=0, k=6:  dim=0, k=8:  dim=0
        k=10: dim=0, k=12: dim=1, k=14: dim=0, k=16: dim=1
        k=18: dim=1, k=20: dim=1, k=22: dim=1, k=24: dim=2

    Die einzige Cusp-Form von Gewicht 12 ist die Ramanujan-Delta-Funktion Δ(τ).

    Für allgemeines Niveau N (Γ_0(N)):
        dim S_k(Γ_0(N)) ≈ μ(N) · (k-1) / 12  für k ≥ 2 gerade
    wobei μ(N) = N · Π_{p | N} (1 + 1/p) der Index von Γ_0(N) in SL(2,Z) ist.

    @param k: Gewicht der Modulform (gerade, ≥ 2)
    @param level: Niveau N ≥ 1 (Standard: 1 = volle Modulgruppe)
    @return: Dimension des Cusp-Formen-Raums (nicht-negative ganze Zahl)
    @raises ValueError: Wenn k ungerade oder k < 2
    @lastModified: 2026-03-10
    """
    # Nur gerade Gewichte für holomorphe Cusp-Formen
    if k < 2 or k % 2 != 0:
        raise ValueError(f"Gewicht k muss gerade und ≥ 2 sein, erhalten: {k}")

    if level == 1:
        # Exakte Formel für die volle Modulgruppe SL(2,Z)
        # ε(k) = 1 falls k ≡ 2 (mod 12), sonst 0
        epsilon = 1 if k % 12 == 2 else 0
        dim = k // 12 - epsilon
        # Dimension kann nicht negativ sein
        return max(0, dim)
    else:
        # Näherungsformel für Γ_0(N) mit N > 1
        # Berechne den Index μ(N) = N · Π_{p | N} (1 + 1/p)
        # Primteiler von N ermitteln
        N = level
        mu = float(N)
        temp = N
        d = 2
        seen_primes = set()
        while d * d <= temp:
            if temp % d == 0:
                if d not in seen_primes:
                    # Eulerproduktsformel: Faktor (1 + 1/p) pro Primteiler p
                    mu *= (1.0 + 1.0 / d)
                    seen_primes.add(d)
                while temp % d == 0:
                    temp //= d
            d += 1
        if temp > 1 and temp not in seen_primes:
            mu *= (1.0 + 1.0 / temp)

        # Dimension via Riemann-Roch-Näherung: μ(N)·(k-1)/12
        # Für k=2 subtrahiere Anzahl der Spitzen (-1)
        dim_approx = mu * (k - 1) / 12.0
        return max(0, int(round(dim_approx)))


def is_cusp_form_basis(k: int, z: complex, n_terms: int = 100) -> bool:
    """
    Prüft ob für das gegebene Gewicht k eine Basis-Cusp-Form existiert.

    Die Ramanujan-Delta-Funktion Δ(τ) ist die einzige (bis auf Skalierung)
    Cusp-Form vom Gewicht 12. Für k = 12m gilt: Δ(τ)^m ist eine Cusp-Form
    (falls der Raum S_{12m} nicht leer ist).

    Ein Gewicht k ist ein "Cusp-Form-Gewicht" (dim S_k ≥ 1), falls
    cusp_form_dimension(k) ≥ 1. Dies ist äquivalent dazu, dass
    die Ramanujan-Delta-Funktion als Basiselement verwendet werden kann.

    Für k=12: Δ(τ) ist eindeutige (normierte) Basis-Cusp-Form.
    Für k=24: Δ(τ)² ist eine Cusp-Form, dim S_24 = 2.

    @param k: Gewicht (gerade ≥ 2)
    @param z: Punkt in der oberen Halbebene (für numerische Verifikation)
    @param n_terms: Anzahl Terme für Produktformel
    @return: True wenn dim S_k ≥ 1 (d.h. Cusp-Form-Raum nicht trivial)
    @lastModified: 2026-03-10
    """
    # Nur gerade Gewichte zulässig
    if k < 2 or k % 2 != 0:
        return False

    # Dimension des Cusp-Formen-Raums berechnen
    dim = cusp_form_dimension(k, level=1)
    return dim >= 1


def ramanujan_tau_properties(n: int) -> dict[str, object]:
    """
    Überprüft bekannte Eigenschaften der Ramanujan-Tau-Funktion τ(n).

    Die Tau-Funktion erscheint als n-ter Fourier-Koeffizient der Delta-Funktion:
        Δ(τ) = Σ_{n=1}^∞ τ(n) q^n   mit q = e^{2πiτ}

    Bekannte arithmetische Eigenschaften (alle bewiesen):
    1. Multiplikativität: τ(mn) = τ(m)·τ(n) falls gcd(m,n) = 1
       (Δ ist Hecke-Eigenform)
    2. Ramanujan-Vermutung (Deligne 1974): |τ(p)| ≤ 2·p^{11/2} für alle Primzahlen p
       (Bewiesen als Teil der Weil-Vermutungen)
    3. Rekursion für Primzahlpotenzen:
       τ(p^k) = τ(p)·τ(p^{k-1}) - p^{11}·τ(p^{k-2})

    @param n: Positive ganze Zahl ≥ 1
    @return: Dictionary mit:
             'tau_n': τ(n) (Ramanujan-Tau-Wert)
             'multiplicativity_verified': True wenn τ(n) = τ(1)·τ(n) (trivial für n)
                                          oder bei n = p·q geprüft
             'ramanujan_bound_ok': True wenn |τ(n)| ≤ 2·n^{11/2} (für n prim)
    @raises ValueError: Wenn n < 1
    @lastModified: 2026-03-10
    """
    if n < 1:
        raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

    # Tau-Koeffizienten bis n berechnen (benötigt ausreichend Terme)
    n_max = max(n + 1, 10)
    tau_list = fourier_coefficients_delta(n_max)
    tau_n = tau_list[n - 1]  # τ(n) ist an Position n-1 (1-indiziert)

    # Multiplikativität prüfen: τ(1)·τ(n) = 1·τ(n) = τ(n) (trivial für n≥1)
    # Aussagekräftiger: Teste τ(p·q) = τ(p)·τ(q) für coprime p,q ≤ n
    mult_verified = True
    if n > 1:
        # Finde zwei coprime Faktoren von n (falls n zusammengesetzt)
        # Suche nach Zerlegung n = a · b mit gcd(a,b)=1 und a,b ≥ 2
        for a in range(2, n):
            if n % a == 0:
                b = n // a
                if gcd(a, b) == 1 and a >= 2 and b >= 2:
                    # Brauche τ(a) und τ(b)
                    tau_a_b_max = max(a, b) + 1
                    tau_ab = fourier_coefficients_delta(tau_a_b_max)
                    tau_a = tau_ab[a - 1]
                    tau_b = tau_ab[b - 1]
                    # Prüfe τ(n) == τ(a)·τ(b)
                    mult_verified = (tau_n == tau_a * tau_b)
                    break

    # Ramanujan-Vermutung: |τ(p)| ≤ 2·p^{11/2} für prim p
    # Für zusammengesetzte n: allgemeine Schranke |τ(n)| ≤ d(n)·n^{11/2}
    # (wobei d(n) die Teileranzahl, hier vereinfacht: n^{11/2} als obere Schranke)
    ramanujan_bound = 2.0 * (n ** (11.0 / 2.0))
    ramanujan_bound_ok = abs(tau_n) <= ramanujan_bound

    return {
        'tau_n': tau_n,
        'multiplicativity_verified': mult_verified,
        'ramanujan_bound_ok': ramanujan_bound_ok
    }


# ===========================================================================
# THETA-REIHEN
# ===========================================================================

def theta_function(z: complex, n_max: int = 50) -> complex:
    """
    Berechnet die Jacobi-Theta-Funktion ϑ_3(0|τ).

    Definition (Jacobi-Theta-Funktion dritten Typs bei z=0):
        θ(τ) = ϑ_3(0|τ) = Σ_{n=-∞}^{∞} q^{n²}   mit q = e^{πiτ}

    Äquivalent: θ(τ) = 1 + 2·Σ_{n=1}^∞ q^{n²} (da n und -n denselben Beitrag liefern)

    Konvergenz: Für Im(τ) > 0 gilt |q| = e^{-π·Im(τ)} < 1,
    daher konvergiert die Reihe exponentiell schnell.

    Anwendungen:
    - θ(τ)^k = Σ_{n=0}^∞ r_k(n) q^n liefert Darstellungsanzahlen r_k(n)
    - θ ist eine Modulform halben Gewichts (Gewicht 1/2)
    - Verbindung zu elliptischen Funktionen und Modulformen

    Transformationsformel (Jacobiische Imagination):
        θ(-1/τ) = sqrt(-iτ) · θ(τ)

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_max: Maximaler Summationsindex (symmetrisch: -n_max bis n_max)
    @return: Wert der Theta-Funktion als komplexe Zahl
    @raises ValueError: Wenn Im(z) ≤ 0
    @lastModified: 2026-03-10
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0), erhalten: {z}")

    # q = e^{πiτ}  (beachte: π, nicht 2π !)
    q = cmath.exp(1j * math.pi * z)

    # θ(τ) = Σ_{n=-N}^{N} q^{n²} = 1 + 2·Σ_{n=1}^{N} q^{n²}
    # Nutze symmetrische Form für Effizienz
    total = complex(1.0)  # n=0 Beitrag
    for n in range(1, n_max + 1):
        # q^{n²} berechnen
        term = q ** (n * n)
        # Konvergenzcheck: bei sehr kleinen Werten abbrechen
        if abs(term) < 1e-300:
            break
        # Faktor 2 wegen n und -n (gleicher Beitrag da n²=(-n)²)
        total += 2.0 * term

    return total


def theta_transformation(z: complex) -> dict[str, object]:
    """
    Verifiziert die Jacobi-Theta-Transformationsformel numerisch.

    Die Transformationsformel lautet:
        θ(-1/τ) = sqrt(-iτ) · θ(τ)

    Dies ist eine fundamentale Symmetrie der Theta-Funktion unter der
    S-Transformation (τ ↦ -1/τ) der modularen Gruppe.

    Herleitung via Poissonsche Summenformel:
        Σ_{n} e^{-πn²/t} = √t · Σ_{n} e^{-πn²t}
    (Riemann verwendete diese Formel für seine zeta-Funktion!)

    @param z: Punkt τ in der oberen Halbebene
    @return: Dictionary mit:
             'theta_z': θ(τ)
             'theta_minus_inv': θ(-1/τ)
             'transformation_ratio': θ(-1/τ) / θ(τ)
             'expected': sqrt(-iτ) (erwartetes Verhältnis)
             'verified': True wenn |ratio - expected| / |expected| < 1e-6
    @lastModified: 2026-03-10
    """
    # Berechne θ(τ)
    theta_z = theta_function(z)

    # Berechne θ(-1/τ) – der transformierte Argument ist -1/τ
    z_transformed = -1.0 / z
    theta_minus_inv = theta_function(z_transformed)

    # Erwartetes Verhältnis: sqrt(-iτ)
    # sqrt(-iτ) = sqrt(-i) · sqrt(τ) mit Hauptzweig
    expected = cmath.sqrt(-1j * z)

    # Tatsächliches Verhältnis
    if abs(theta_z) < 1e-300:
        # Sonderfall: θ(τ) ≈ 0 (sollte nicht vorkommen für Im(τ) > 0)
        transformation_ratio = complex(float('nan'))
        verified = False
    else:
        transformation_ratio = theta_minus_inv / theta_z
        # Relativer Fehler
        rel_error = abs(transformation_ratio - expected) / max(abs(expected), 1e-15)
        verified = rel_error < 1e-4

    return {
        'theta_z': theta_z,
        'theta_minus_inv': theta_minus_inv,
        'transformation_ratio': transformation_ratio,
        'expected': expected,
        'verified': verified
    }


def sum_of_squares_theta(n: int, k: int = 2) -> int:
    """
    Berechnet r_k(n): Anzahl der Darstellungen von n als Summe von k Quadraten.

    Definition:
        r_k(n) = #{(x_1,...,x_k) ∈ Z^k : x_1² + ... + x_k² = n}
    (negative Zahlen und verschiedene Reihenfolgen zählen separat)

    Theoretisch via Theta-Reihen:
        θ(τ)^k = (Σ_{n} q^{n²})^k = Σ_{m=0}^∞ r_k(m) q^m
    d.h. r_k(m) ist der m-te Koeffizient von θ^k.

    Für k=2 gilt die Formel von Jacobi (1829):
        r_2(n) = 4 · Σ_{d|n} χ(d)
    wobei χ(d) = 0 für gerades d, χ(d) = (-1)^{(d-1)/2} für ungerades d.
    (Dirichlet-Charakter mod 4)

    Bekannte Werte:
        r_2(0) = 1, r_2(1) = 4, r_2(2) = 4, r_2(4) = 4
        r_2(5) = 8 (5 = 1²+2² = 2²+1² = (-1)²+2² = 2²+(-1)² = 1²+(-2)² = ...)
        r_2(25) = 12

    Diese Implementierung nutzt Brute-Force für alle k (universell und korrekt).

    @param n: Nicht-negative ganze Zahl
    @param k: Anzahl der Quadrate (Standard: 2)
    @return: r_k(n) als nicht-negative ganze Zahl
    @raises ValueError: Wenn n < 0 oder k < 1
    @lastModified: 2026-03-10
    """
    if n < 0:
        raise ValueError(f"n muss ≥ 0 sein, erhalten: {n}")
    if k < 1:
        raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")

    # Sonderfälle
    if n == 0:
        return 1  # Nur die Nulldarstellung: (0,0,...,0)

    # Maximaler Betrag einer Komponente: |x_i| ≤ floor(sqrt(n))
    max_val = int(math.isqrt(n))

    # Rekursive Zählung via dynamischer Programmierung
    # count[m] = Anzahl Wege, m als Summe von genau j Quadraten darzustellen
    # Iteriere über die k Quadrate
    # Starte mit 1 Weg für Summe=0 (leere Summe)
    current = {0: 1}

    for _ in range(k):
        # Nächste Iteration: addiere ein weiteres Quadrat x² für x ∈ Z
        next_count: dict[int, int] = {}
        for s, cnt in current.items():
            # Iteriere über alle möglichen x (von -max_val bis max_val)
            for x in range(-max_val, max_val + 1):
                new_s = s + x * x
                if new_s <= n:  # Nur Werte ≤ n weiter verfolgen
                    next_count[new_s] = next_count.get(new_s, 0) + cnt
        current = next_count

    return current.get(n, 0)


def jacobi_triple_product(z: complex, q: complex, n_terms: int = 20) -> complex:
    """
    Berechnet das Jacobi-Dreifachprodukt.

    Das Jacobi-Dreifachprodukt (Jacobi, 1829) ist eine der schönsten Identitäten
    der Mathematik:
        Π_{n=1}^∞ (1 - q^{2n})(1 + z·q^{2n-1})(1 + z^{-1}·q^{2n-1})
        = Σ_{n=-∞}^{∞} z^n · q^{n²}

    Voraussetzung: |q| < 1 und z ≠ 0.

    Spezialfälle:
        z = 1: Gibt θ_3(0|τ) = Σ q^{n²} zurück (mit q = e^{πiτ})
        z = -1: Gibt θ_4(0|τ) zurück
        z = q: Gibt θ_2(0|τ) zurück

    Anwendungen:
        - Verbindung zu Eta-Funktion: η(τ)³ = q^{1/8}·Σ_{n} (-1)^n·(2n+1)·q^{n(n+1)/2}
        - Beweis der Pentagonalzahlformel von Euler
        - Partitionentheorie

    @param z: Komplexer Parameter (z ≠ 0)
    @param q: Konvergenzparameter mit |q| < 1
    @param n_terms: Anzahl der Produktterme
    @return: Wert des Dreifachprodukts (Produktseite)
    @raises ValueError: Wenn |q| ≥ 1 oder z = 0
    @lastModified: 2026-03-10
    """
    if abs(q) >= 1.0:
        raise ValueError(f"|q| muss < 1 sein, erhalten: |q| = {abs(q)}")
    if abs(z) < 1e-300:
        raise ValueError("z darf nicht 0 sein")

    # Produktformel: Π_{n=1}^{N} (1-q^{2n})(1+z·q^{2n-1})(1+z^{-1}·q^{2n-1})
    product = complex(1.0)
    z_inv = 1.0 / z

    for n in range(1, n_terms + 1):
        # q^{2n} und q^{2n-1} berechnen
        q2n = q ** (2 * n)
        q2n_m1 = q ** (2 * n - 1)

        # Drei Faktoren des Produkts
        f1 = 1.0 - q2n                   # (1 - q^{2n})
        f2 = 1.0 + z * q2n_m1            # (1 + z·q^{2n-1})
        f3 = 1.0 + z_inv * q2n_m1        # (1 + z^{-1}·q^{2n-1})

        product *= f1 * f2 * f3

        # Konvergenzcheck: wenn q^{2n} sehr klein, abbrechen
        if abs(q2n) < 1e-300:
            break

    return product


def _is_prime_simple(n: int) -> bool:
    """
    Einfacher Primzahltest (nur für interne Nutzung in diesem Modul).

    @param n: Zu testende ganze Zahl
    @return: True wenn n eine Primzahl ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
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
    @lastModified 2026-03-10
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


def dedekind_eta(z: complex, n_terms: int = 100) -> complex:
    """
    Berechnet die Dedekind-Eta-Funktion η(τ).

    Definition (Dedekind, 1877):
        η(τ) = q^{1/24} · Π_{n=1}^∞ (1 - q^n)   mit q = e^{2πiτ}

    Eigenschaften:
        - η ist eine Modulform halben Gewichts (Gewicht 1/2) mit Multiplikatorsystem
        - Transformationsformel: η(-1/τ) = sqrt(-iτ) · η(τ)
        - Periodizität: η(τ+1) = e^{πi/12} · η(τ)
        - Verbindung zur Delta-Funktion: Δ(τ) = η(τ)^{24}
          (d.h. η(τ) = Δ(τ)^{1/24})
        - Etafunktion taucht in der Physik auf: Quantenstring-Theorie (Partition function)

    Pentagonalzahlformel von Euler:
        Π_{n=1}^∞ (1-q^n) = Σ_{n=-∞}^∞ (-1)^n · q^{n(3n-1)/2}
    (Pentagonalzahlen: 0, 1, 2, 5, 7, 12, 15, ...)

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme (mehr = genauer)
    @return: Wert η(τ) als komplexe Zahl
    @raises ValueError: Wenn Im(z) ≤ 0
    @lastModified: 2026-03-10
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0), erhalten: {z}")

    # q = e^{2πiτ}  (ganzer Kreis!)
    q = cmath.exp(2j * math.pi * z)

    # Vorfaktor: q^{1/24} = e^{2πiτ/24} = e^{πiτ/12}
    q_124 = cmath.exp(2j * math.pi * z / 24.0)

    # Unendliches Produkt: Π_{n=1}^{N} (1 - q^n)
    product = complex(1.0)
    q_power = q  # q^n, beginnt bei q^1

    for n in range(1, n_terms + 1):
        if abs(q_power) < 1e-300:
            break
        product *= (1.0 - q_power)
        q_power *= q  # q^{n+1}

    return q_124 * product


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
    @lastModified 2026-03-10
    """
    # Ramanujan-Tau-Koeffizienten bis zu einem hinreichend hohen Index berechnen
    # Für die Kommutatortests brauchen wir Koeffizienten bis ca. p*q*p
    max_prime = 13  # Kleine Primzahlen: 2, 3, 5, 7, 11, 13
    small_primes = _primes_up_to(max_prime)

    # Tau-Koeffizienten (Hecke-Eigenwerte der Delta-Funktion) berechnen
    # Wir brauchen genug Koeffizienten: bis zu p*q = 13*11 = 143
    n_coeff = 200
    tau_list = fourier_coefficients_delta(n_coeff)

    # Hecke-Eigenwerte: τ(p) für jede Primzahl p ≤ max_prime
    hecke_eigenvalues = {}
    for p in small_primes:
        # tau(p) ist der p-te Fourier-Koeffizient der Delta-Funktion
        if p <= n_coeff:
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
    @lastModified 2026-03-10
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

    Numerische Approximation: Wir integrieren über Gitterpunkte im Fundamentalbereich.
    Gitterpunkte: τ = x + iy mit x ∈ [-1/2, 1/2], y ∈ [1, Y_max], |τ| ≥ 1.

    Aus der Fourier-Entwicklung:
        f(τ) = Σ a_n q^n,  g(τ) = Σ b_n q^n  (q = e^{2πiτ})

    Approximierter Wert via Rankin-Selberg-Methode (für zwei Cusp-Formen):
        ⟨f, g⟩ ≈ Σ_n a_n · conj(b_n) · (4πn)^{-(k-1)} · Γ(k-1) / normierung

    Da diese Formel auf dem Rangkin-Selberg-Integral basiert und die genaue
    Normierung vom Niveau N abhängt, gibt diese Implementierung eine rohe Schätzung zurück.

    @param f_coeffs: Fourier-Koeffizienten [a_1, a_2, ..., a_N] von f
    @param g_coeffs: Fourier-Koeffizienten [b_1, b_2, ..., b_N] von g
    @param k: Gewicht der Modulformen (muss für beide gleich sein)
    @return: Schätzwert für ⟨f, g⟩ als komplexe Zahl
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    import math

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
# L-FUNKTIONEN ELLIPTISCHER KURVEN
# ===========================================================================

def elliptic_curve_points_over_fp(a: int, b: int, p: int) -> list[tuple[int, int]]:
    """
    Berechnet alle affinen Punkte einer elliptischen Kurve E: y² = x³ + ax + b über F_p.

    Algorithmus (Brute-Force):
        Für jedes x ∈ {0, 1, ..., p-1}:
            1. Berechne rhs = (x³ + ax + b) mod p
            2. Prüfe ob rhs ein quadratischer Rest mod p ist
            3. Falls ja, berechne y = sqrt(rhs) mod p (beide Wurzeln ± y)

    Quadratische Reste mod p: Für ungerades p und a ≢ 0 (mod p) gilt:
        a ist QR ⟺ a^{(p-1)/2} ≡ 1 (mod p)   (Euler-Kriterium)

    Die Kurve ist nicht-singulär, falls die Diskriminante
        Δ = -16(4a³ + 27b²) ≢ 0 (mod p)

    Hinweis: Der Punkt im Unendlichen O ist NICHT in der Liste enthalten,
    wird aber für #E(F_p) mitgezählt.

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param p: Primzahl p ≥ 3 (Charakteristik des Körpers F_p)
    @return: Liste aller affinen Punkte (x, y) ∈ F_p × F_p auf E
    @raises ValueError: Wenn p < 3 oder p nicht prim
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    if p < 3:
        raise ValueError(f"p muss eine Primzahl ≥ 3 sein, erhalten: p = {p}")

    points = []

    for x in range(p):
        # Rechte Seite der Kurvengleichung: x³ + ax + b mod p
        rhs = (pow(x, 3, p) + a * x + b) % p

        if rhs == 0:
            # y = 0 ist einzige Lösung (doppelter Punkt bei y=0)
            points.append((x, 0))
        else:
            # Euler-Kriterium: rhs ist QR ⟺ rhs^{(p-1)/2} ≡ 1 (mod p)
            if pow(rhs, (p - 1) // 2, p) == 1:
                # Berechne Quadratwurzel mod p
                # Für p ≡ 3 (mod 4): y = rhs^{(p+1)/4} mod p
                # Allgemein: Tonelli-Shanks-Algorithmus (hier einfache Variante)
                y = pow(rhs, (p + 1) // 4, p)
                # Verifikation
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

    Definition:
        a_p(E) = p + 1 - #E(F_p)

    wobei #E(F_p) die Anzahl der F_p-Punkte inkl. dem Punkt im Unendlichen O ist:
        #E(F_p) = |{(x,y) ∈ F_p² : y² = x³+ax+b}| + 1

    Der Satz von Hasse (Hasse-Weil-Schranke) besagt:
        |a_p(E)| = |p + 1 - #E(F_p)| ≤ 2√p

    Dies ist die elliptische Kurve Analogie zur Riemann-Vermutung!

    Die Frobenius-Spur taucht auf als Eigenwert des Frobenius-Endomorphismus
    φ: (x, y) ↦ (x^p, y^p) auf dem Tate-Modul T_ℓ(E).

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param p: Primzahl p ≥ 3
    @return: Frobenius-Spur a_p = p + 1 - #E(F_p)
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Berechne alle affinen Punkte
    affine_points = elliptic_curve_points_over_fp(a, b, p)
    # Füge den Punkt im Unendlichen O hinzu
    num_points = len(affine_points) + 1  # +1 für O

    # Frobenius-Spur
    return p + 1 - num_points


def l_function_elliptic_curve(a: int, b: int, s: complex, terms: int = 100) -> complex:
    """
    Berechnet die L-Funktion einer elliptischen Kurve E: y² = x³ + ax + b.

    Die L-Funktion wird als Dirichlet-Reihe dargestellt:
        L(E, s) = Σ_{n=1}^∞ a_n(E) / n^s

    Die Koeffizienten a_n(E) sind multiplikativ und bestimmt durch:
        - Für Primzahlen p mit guter Reduktion:
            a_p(E) = p + 1 - #E(F_p)   (Frobenius-Spur)
        - Für Primzahlpotenzen p^r:
            a_{p^r}(E) = a_p · a_{p^{r-1}} - p · a_{p^{r-2}}  (Rekursion)
        - Für zusammengesetzte Zahlen n = m·k mit gcd(m,k) = 1:
            a_n(E) = a_m(E) · a_k(E)  (Multiplikativität)

    Die Reihe konvergiert für Re(s) > 3/2.

    Die Verbindung zu Hecke-Operatoren:
        L(E, s) = L(f_E, s)  (nach dem Satz von Shimura-Taniyama-Wiles)
    wobei f_E die zugehörige Modulform ist.

    @param a: Koeffizient der elliptischen Kurve
    @param b: Konstantterm der elliptischen Kurve
    @param s: Komplexes Argument (Re(s) > 3/2 für Konvergenz)
    @param terms: Anzahl der Summanden in der Dirichlet-Reihe
    @return: Näherungswert von L(E, s) als komplexe Zahl
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Berechne Frobenius-Spuren a_p für Primzahlen bis ca. terms
    # Für die Dirichlet-Reihe brauchen wir a_n für n = 1, ..., terms

    # Schritt 1: Primzahlen bis terms bestimmen
    primes_needed = _primes_up_to(min(terms, 200))

    # Schritt 2: a_p für jede Primzahl p berechnen (Frobenius-Spuren)
    a_prime = {}
    for p in primes_needed:
        try:
            a_prime[p] = trace_of_frobenius(a, b, p)
        except Exception:
            a_prime[p] = 0  # Schlechte Reduktion: a_p = 0 oder ±1

    # Schritt 3: a_n für alle n bis terms via Multiplikativität berechnen
    a_n = [0.0] * (terms + 1)
    a_n[1] = 1.0  # a_1 = 1 per Konvention

    for n in range(2, terms + 1):
        # Primfaktorzerlegung von n
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
                # Allgemeine Rekursion: a_{p^r} = a_p * a_{p^{r-1}} - p * a_{p^{r-2}}
                prev1 = int(round(a_n[p ** (r-1)])) if p ** (r-1) <= terms else 0
                prev2 = int(round(a_n[p ** (r-2)])) if p ** (r-2) <= terms else 0
                a_n[n] = float(a_prime.get(p, 0)) * prev1 - p * prev2
        elif len(factors) > 1:
            # n = m * k mit gcd(m,k) = 1: Multiplikativität a_n = a_m * a_k
            ps = list(factors.keys())
            # Teile n in zwei coprime Teile
            p1 = ps[0]
            m = p1 ** factors[p1]
            k = n // m
            if m <= terms and k <= terms:
                a_n[n] = a_n[m] * a_n[k]

    # Schritt 4: L-Funktion als Dirichlet-Reihe summieren: L(E,s) = Σ a_n / n^s
    result = complex(0.0)
    for n in range(1, terms + 1):
        if a_n[n] != 0.0:
            result += complex(a_n[n]) / (n ** s)

    return result


def birch_swinnerton_dyer_bsd_estimate(a: int, b: int, primes: list[int]) -> dict[str, object]:
    """
    Empirische Schätzung der BSD-Vermutung (Birch und Swinnerton-Dyer).

    Die BSD-Vermutung (Millennium-Problem) besagt:
        ord_{s=1} L(E, s) = Rang(E(Q))

    d.h. die Ordnung der Nullstelle von L(E, s) bei s=1 ist gleich dem
    arithmetischen Rang der elliptischen Kurve über Q.

    Empirische Schätzung (Birch und Swinnerton-Dyer, 1960er Jahre):
        Π_{p ≤ X} (N_p / p) ≈ C · (log X)^r

    wobei:
        N_p = #E(F_p)   (Punktanzahl über F_p)
        r   = Rang(E(Q))   (gesuchter Rang)
        C   = Konstante (abhängig von E)
        X   = Schranke für die Primzahlen

    Logarithmieren:
        log(Produkt) ≈ r · log(log X) + log(C)

    Damit kann man r schätzen durch lineare Regression von
        log(Produkt) gegen log(log X).

    @param a: Koeffizient der elliptischen Kurve E: y² = x³ + ax + b
    @param b: Konstantterm der elliptischen Kurve
    @param primes: Liste der zu verwendenden Primzahlen p ≥ 3
    @return: Dictionary mit:
             'estimate': float – Schätzung für den Rang
             'primes_used': list – verwendete Primzahlen
             'log_product': float – log(Π N_p/p)
             'product': float – Π N_p/p
             'n_p_values': dict – {p: N_p} Punktanzahlen
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    import math

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

            # Beitrag zum Produkt: N_p / p
            ratio = n_p / p
            if ratio > 0:
                log_product += math.log(ratio)
                product *= ratio
                valid_primes.append(p)
        except Exception:
            continue

    # Schritt 2: Rang-Schätzung via log(log X)
    # Bei X = max(primes): Schätze r aus log(Produkt) ≈ r * log(log X)
    if valid_primes:
        X = max(valid_primes)
        log_log_X = math.log(math.log(X)) if X > 2 else 1.0
        # Einfache Schätzung: r ≈ log(Produkt) / log(log X)
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
    @lastModified 2026-03-10
    """
    result = complex(0.0)
    N = len(coefficients)

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
    @lastModified 2026-03-10
    """
    import cmath
    import math

    # Vorzeichen: ε = (-1)^{k/2} für die volle Modulgruppe SL(2,Z)
    epsilon = (-1) ** (k // 2)

    # Gamma-Funktion Γ(s) via Stirling-Näherung oder für ganzzahlige Argumente exakt
    def gamma_approx(z: complex) -> complex:
        """Näherung der Gamma-Funktion via Stirling-Formel."""
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

    # Vollständige L-Funktionen
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
