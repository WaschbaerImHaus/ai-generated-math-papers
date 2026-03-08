"""
@file modular_forms.py
@brief Modulformen – Grundlagen und fundamentale Konzepte.
@description
    Implementiert grundlegende Konzepte der Theorie der Modulformen:
    - Die modulare Gruppe SL(2,Z) und Möbius-Transformationen
    - Eisenstein-Reihen G_k(τ) und normierte Versionen E_4, E_6
    - Ramanujan-Delta-Funktion Δ(τ) (Diskriminantenform)
    - j-Invariante j(τ)
    - Ramanujan-Tau-Funktion τ(n)
    - Hecke-Operatoren T_p
    - Shimura-Taniyama-Wiles-Verbindung (Demonstration)

    Modulformen sind holomorphe Funktionen auf der oberen Halbebene H = {τ : Im(τ) > 0},
    die unter Möbius-Transformationen γ ∈ SL(2,Z) folgende Transformationsregel erfüllen:
        f(γτ) = (cτ + d)^k · f(τ)
    wobei k das Gewicht der Modulform ist.

    Diese Theorie ist zentral für den Beweis des Großen Fermatschen Satzes (Wiles 1995),
    da elliptische Kurven und Modulformen via Shimura-Taniyama-Wiles verknüpft sind.

@author Kurt Ingwer
@version 1.0
@since 2026-03-08
@lastModified 2026-03-08
"""

import math
import cmath
import numpy as np
from typing import Callable


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
