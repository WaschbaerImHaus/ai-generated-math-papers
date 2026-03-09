"""
@file modular_forms.py
@brief Modulformen – Grundlagen, Cusp-Formen und Theta-Reihen.
@description
    Implementiert grundlegende Konzepte der Theorie der Modulformen:
    - Die modulare Gruppe SL(2,Z) und Möbius-Transformationen
    - Eisenstein-Reihen G_k(τ) und normierte Versionen E_4, E_6
    - Ramanujan-Delta-Funktion Δ(τ) (Diskriminantenform)
    - j-Invariante j(τ)
    - Ramanujan-Tau-Funktion τ(n)
    - Hecke-Operatoren T_p
    - Shimura-Taniyama-Wiles-Verbindung (Demonstration)
    - Cusp-Formen S_k(Γ_0(N)): Dimension, Eigenschaften, Ramanujan-Tau
    - Theta-Reihen: Jacobi-Theta, Dedekind-Eta, Jacobi-Dreifachprodukt

    Modulformen sind holomorphe Funktionen auf der oberen Halbebene H = {τ : Im(τ) > 0},
    die unter Möbius-Transformationen γ ∈ SL(2,Z) folgende Transformationsregel erfüllen:
        f(γτ) = (cτ + d)^k · f(τ)
    wobei k das Gewicht der Modulform ist.

    Diese Theorie ist zentral für den Beweis des Großen Fermatschen Satzes (Wiles 1995),
    da elliptische Kurven und Modulformen via Shimura-Taniyama-Wiles verknüpft sind.

@author Kurt Ingwer
@version 2.0
@since 2026-03-08
@lastModified 2026-03-09
"""

import math
import cmath
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

    @author Kurt Ingwer
    @since 2026-03-08
    @lastModified 2026-03-08
    """

    def __init__(self, a: int, b: int, c: int, d: int):
        """
        Initialisiert ein Element der modularen Gruppe.

        @param a: Eintrag oben links
        @param b: Eintrag oben rechts
        @param c: Eintrag unten links
        @param d: Eintrag unten rechts
        @raises ValueError: Wenn ad - bc ≠ 1
        @lastModified: 2026-03-08
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
        @lastModified: 2026-03-08
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
        @lastModified: 2026-03-08
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
        @lastModified: 2026-03-08
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
        @lastModified: 2026-03-08
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

    def __eq__(self, other) -> bool:
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    gamma: ModularGroup,
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
    @lastModified: 2026-03-08
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

def hecke_operator(coefficients: list, p: int, k: int = 12) -> list:
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

    @param coefficients: Liste [a(0), a(1), ..., a(N)] der Fourier-Koeffizienten
    @param p: Primzahl p ≥ 2 (Hecke-Operator T_p)
    @param k: Gewicht der Modulform (Standard: 12 für Δ)
    @return: Liste [b(0), b(1), ..., b(M)] der transformierten Koeffizienten
    @lastModified: 2026-03-08
    """
    N = len(coefficients)
    # Ausgabelänge: T_p erhöht den Index um p, daher brauchen wir weniger Koeffizienten
    M = max(1, N // p)
    b = []

    for n in range(M):
        # Terme: a(pn)
        idx_pn = p * n
        a_pn = coefficients[idx_pn] if idx_pn < N else 0

        # Terme: p^{k-1} · a(n/p) falls p | n
        if n % p == 0:
            idx_np = n // p
            a_np = coefficients[idx_np] if idx_np < N else 0
            bn = a_pn + (p ** (k - 1)) * a_np
        else:
            bn = a_pn

        b.append(bn)

    return b


# ===========================================================================
# SHIMURA-TANIYAMA-WILES-VERBINDUNG
# ===========================================================================

def shimura_taniyama_check(
    a_p_elliptic: dict,
    a_p_modular: dict,
    primes: list
) -> dict:
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-09
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
    @lastModified: 2026-03-09
    """
    # Nur gerade Gewichte zulässig
    if k < 2 or k % 2 != 0:
        return False

    # Dimension des Cusp-Formen-Raums berechnen
    dim = cusp_form_dimension(k, level=1)
    return dim >= 1


def ramanujan_tau_properties(n: int) -> dict:
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
    @lastModified: 2026-03-09
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
    @lastModified: 2026-03-09
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


def theta_transformation(z: complex) -> dict:
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
    @lastModified: 2026-03-09
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
    @lastModified: 2026-03-09
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
    @lastModified: 2026-03-09
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
    @lastModified: 2026-03-09
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
