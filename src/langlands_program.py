"""
@file langlands_program.py
@brief Umfassendes Python-Modul für das Langlands-Programm.
@description
    Implementiert zentrale Konzepte des Langlands-Programms:
    - Satake-Isomorphismus für GL_n
    - Lokale Langlands-Korrespondenz für GL_1 und GL_2
    - Globale Langlands-Korrespondenz (Modularitätssatz)
    - Langlands-Funktorialität (Basiswechsel, sym. Potenzen)
    - Artin-Reziprozität und Artin-L-Funktionen
    - Arthur-Selberg-Spurformel (numerische Demo)
    - Langlands-Dualgruppen

    Das Langlands-Programm verbindet:
    - Zahlentheorie (Galois-Darstellungen)
    - Automorphe Formen (Modulformen, π auf GL_n)
    - Darstellungstheorie (Weil-Gruppe, L-Gruppen)

    Grundlegende Korrespondenz:
        Galois-Darstellungen ρ: Gal(Q̄/Q) → GL_n(ℂ)
        ↔ Automorphe Darstellungen π auf GL_n(A_Q)

    Für GL_1: Klassenkörpertheorie (Artin-Reziprozität)
        Art_p: Q_p^× → Gal(Q_p^ab / Q_p)^∨ ≅ W_p^ab

    Für GL_2: Modulformen f ∈ S_k(Γ_0(N))
        ↔ cuspidal automorphe Darstellungen π auf GL_2(A_Q)

@author Michael Fuhrmann
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
import numpy as np
import sympy as sp
from sympy import (
    symbols, primerange, isprime, factorint, mobius,
    sqrt as sp_sqrt, Rational, pi as sp_pi, I as sp_I,
    cos, sin, exp, log, Abs, simplify, conjugate, zeta,
    Symbol, Poly, GF, cyclotomic_poly, totient, gcd as sp_gcd
)
from typing import Dict, List, Tuple, Optional, Union, Any, Callable
from functools import lru_cache
from dataclasses import dataclass

# Gemeinsame Hilfsfunktionen aus zentralem Modul importieren
from math_helpers import (
    is_prime as _is_prime,
    prime_factorization as _prime_factorization,
    primes_up_to as _primes_up_to,
    euler_phi as _euler_phi,
    gcd as _gcd,
    mod_inverse as _mod_inverse,
    legendre_symbol as _legendre_symbol
)


# ---------------------------------------------------------------------------
# Hilfsfunktionen (Modul-Ebene)
# ---------------------------------------------------------------------------

def local_langlands_parameter(
    p: int,
    a_p: float,
    b_p: float = 0.0
) -> Dict[str, Any]:
    """
    @brief Berechnet den Weil-Deligne-Parameter aus Hecke-Eigenwerten.
    @description
        Für eine automorphe Darstellung π_p auf GL_2(Q_p) mit
        Hecke-Eigenwerten a_p (T_p-Eigenwert) und b_p (Diagonalbeitrag)
        liefert der lokale Langlands-Parameter ein Element der
        Weil-Deligne-Gruppe W'_{Q_p} = W_{Q_p} × SL_2(ℂ).

        Satake-Parameter für unitäre Normierung (Gewicht k=2):
            α_p + β_p = a_p,  α_p · β_p = p  (Ramanujan: |α_p|=√p)
            Charakt. Polynom: X² - a_p·X + p = 0

        Discriminante: Δ = a_p² - 4p
        Falls Δ < 0: α_p, β_p ∈ ℂ (komplexe Satake-Parameter, typisch)

    @param p Primzahl (lokale Stelle)
    @param a_p Hecke-Eigenwert (Frobenius-Spur)
    @param b_p optionaler zweiter Hecke-Eigenwert (Standard 0)
    @return Dict mit Schlüsseln: p, a_p, alpha, beta, |alpha|, discriminant, type
    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """
    # Charakteristisches Polynom X² - a_p·X + p (unitäre Normierung Gewicht 2)
    disc = a_p ** 2 - 4 * p

    if disc < 0:
        # Komplexe Satake-Parameter (Ramanujan-Vermutung: |α_p| = √p)
        real_part = a_p / 2.0
        imag_part = math.sqrt(-disc) / 2.0
        alpha = complex(real_part, imag_part)
        beta = complex(real_part, -imag_part)
        param_type = "elliptic"
    elif disc == 0:
        # Doppelte reelle Wurzel
        alpha = complex(a_p / 2.0, 0.0)
        beta = alpha
        param_type = "parabolic"
    else:
        # Zwei verschiedene reelle Wurzeln
        sqrt_disc = math.sqrt(disc)
        alpha = complex((a_p + sqrt_disc) / 2.0, 0.0)
        beta = complex((a_p - sqrt_disc) / 2.0, 0.0)
        param_type = "hyperbolic"

    return {
        "p": p,
        "a_p": a_p,
        "alpha": alpha,
        "beta": beta,
        "|alpha|": abs(alpha),
        "discriminant": disc,
        "type": param_type,
        # Ramanujan-Vermutung: |α_p| = √p erwartet
        "ramanujan_satisfied": abs(abs(alpha) - math.sqrt(p)) < 1e-8,
    }


def epsilon_factor_normalization(
    conductor: int,
    s: complex
) -> complex:
    """
    @brief Normierung des ε-Faktors ε(s, π) in der Funktionalgleichung.
    @description
        Die vollständige L-Funktion Λ(π, s) erfüllt:
            Λ(π, s) = ε(π) · Λ(π̌, 1-s)

        Für ein automorphe Darstellung π mit Führer N:
            ε(π, s) = ε(π, 1/2) · N^{1/2 - s}

        Dabei ist ε(π, 1/2) ∈ {±1, ±i} das Wurzelzahl ("root number").
        Für Modulformen: ε(f) = (-1)^k · W_N, wo W_N das Atkin-Lehner-Vorzeichen.

        Normierung (Deligne): ε(s, π) = ε(1/2, π) · N^{1/2-s}

    @param conductor Führer N der automorphen Darstellung
    @param s komplexe Variable
    @return Normierter ε-Faktor als komplexe Zahl
    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """
    # Root number: für Demonstration ε(1/2, π) = +1 (typisch für CM-Kurven)
    root_number = 1.0 + 0j

    # N^{1/2 - s}: Exponentialfunktion in s
    exponent = 0.5 - s
    eps = root_number * (conductor ** exponent)
    return eps


def langlands_correspondence_table(
    E_ap_list: List[Tuple[int, int]]
) -> List[Dict[str, Any]]:
    """
    @brief Erstellt eine Tabelle: p → (a_p, α_p, β_p, |α_p|) für elliptische Kurven.
    @description
        Für eine elliptische Kurve E/Q und die zugehörige Modulform f:
        - a_p = p + 1 - #E(F_p) (Frobenius-Spur)
        - α_p, β_p: Satake-Parameter (Wurzeln von X² - a_p·X + p)
        - Ramanujan-Petersson-Vermutung: |α_p| = |β_p| = √p

        Beispiel: E: y² = x³ - x (Leiter N=32)
            p≡1(4): a_p = 2·Re(π_p), wo π_p∈Z[i] mit π_p·π̄_p=p, π_p≡1 mod (1+i)³
            p≡3(4): a_p = 0 (Frobenius wirkt als komplexe Konjugation → Spur=0)
            p=2:    a_p = 0 (schlechte Reduktion, additive Reduktion)

    @param E_ap_list Liste von (p, a_p)-Paaren
    @return Liste von Dicts mit vollständigen Langlands-Parametern
    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """
    table = []
    for p, a_p in E_ap_list:
        params = local_langlands_parameter(p, float(a_p))
        entry = {
            "p": p,
            "a_p": a_p,
            "alpha_p": params["alpha"],
            "beta_p": params["beta"],
            "|alpha_p|": params["|alpha|"],
            "sqrt_p": math.sqrt(p),
            "ramanujan": params["ramanujan_satisfied"],
            "type": params["type"],
        }
        table.append(entry)
    return table


# ---------------------------------------------------------------------------
# 1. Satake-Isomorphismus
# ---------------------------------------------------------------------------

class SatakeIsomorphism:
    """
    @brief Satake-Isomorphismus H(G//K) ≅ ℂ[X*(T̂)^W] für GL_n.
    @description
        Der Satake-Isomorphismus beschreibt die sphärische Hecke-Algebra
        H(G(Q_p)//G(Z_p)) als kommutative Algebra:

            S: H(G//K) → ℂ[X*(T̂)^W]

        Für G = GL_n(Q_p), K = GL_n(Z_p):
        - T = diagonale Tori, T̂ = ℂ^× × ... × ℂ^× (n-fach)
        - W = S_n (symmetrische Gruppe, Weyl-Gruppe von GL_n)
        - X*(T̂)^W = Polynome in e_1,...,e_n (elementare Symmetrische)

        Sphärische Funktion φ_z für z = (z_1,...,z_n) ∈ (ℂ^×)^n:
            φ_z(g) = ∫_K χ_z(kg) dk  (Harish-Chandra-Integral)

        Für GL_2 mit Satake-Parameter z = (α, β):
            φ_{α,β}(p^k) = Σ_{j=0}^{k} α^j · β^{k-j}  (k ≥ 0)

        Satake-Parameter aus Hecke-Eigenwerten (Gewicht k=2, Modulform):
            a_p = α_p + β_p,  p = α_p · β_p

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, n: int = 2):
        """
        @brief Initialisiert den Satake-Isomorphismus für GL_n.
        @param n Rang der Gruppe GL_n (Standard: 2)
        """
        self.n = n  # Dimension der Gruppe

    def satake_parameters_from_hecke(
        self,
        a_p: float,
        p: int,
        weight: int = 2
    ) -> Tuple[complex, complex]:
        """
        @brief Berechnet Satake-Parameter (α_p, β_p) aus Hecke-Eigenwert a_p.
        @description
            Für eine Modulform f ∈ S_k(Γ_0(N)) mit T_p-Eigenwert a_p:
                Satake-Parameter (normiert auf unitäre Darstellung):
                α_p = p^{(k-1)/2} · e^{iθ_p}
                β_p = p^{(k-1)/2} · e^{-iθ_p}
                wobei a_p = p^{(k-1)/2} · 2·cos(θ_p)

            Für k=2 (elliptische Kurven): α_p·β_p = p, α_p+β_p = a_p.
            Ramanujan-Petersson: |α_p| = p^{(k-1)/2}

        @param a_p Hecke-Eigenwert (Frobenius-Spur)
        @param p Primzahl
        @param weight Gewicht k der Modulform (Standard: 2)
        @return Tupel (α_p, β_p) der Satake-Parameter
        """
        # Normierungsfaktor p^{(k-1)/2}
        norm = p ** ((weight - 1) / 2.0)

        # a_p = 2·norm·cos(θ_p)
        cos_theta = a_p / (2.0 * norm)

        if abs(cos_theta) <= 1.0:
            # Elliptischer Fall (Ramanujan erfüllt): komplexe Satake-Parameter
            theta = math.acos(cos_theta)
            alpha = norm * cmath.exp(1j * theta)
            beta = norm * cmath.exp(-1j * theta)
        else:
            # Hyperbolischer Fall (Ramanujan verletzt, z.B. Eisenstein-Reihen)
            # α, β reell und verschieden
            disc = a_p ** 2 - 4 * (norm ** 2)
            alpha = complex((a_p + math.sqrt(abs(disc))) / 2.0, 0.0)
            beta = complex((a_p - math.sqrt(abs(disc))) / 2.0, 0.0)

        return alpha, beta

    def spherical_function(
        self,
        alpha: complex,
        beta: complex,
        p: int,
        k: int
    ) -> complex:
        """
        @brief Berechnet die sphärische Funktion φ_{α,β}(p^k) für GL_2.
        @description
            Die sphärische Funktion ist das Harish-Chandra-Integral:
                φ_z(p^k) = Σ_{j=0}^{k} α^j · β^{k-j}  (für k ≥ 0)

            Dies ist das Schur-Polynom s_{(k)}(α, β):
                s_k(α,β) = (α^{k+1} - β^{k+1}) / (α - β)  falls α ≠ β
                s_k(α,β) = (k+1) · α^k                    falls α = β

            Für die Hecke-Algebra: T_{p^k} hat Eigenwert φ_z(p^k).

        @param alpha erster Satake-Parameter α
        @param beta zweiter Satake-Parameter β
        @param p Primzahl
        @param k Exponent in p^k
        @return Wert der sphärischen Funktion φ_{α,β}(p^k)
        """
        if abs(alpha - beta) > 1e-12:
            # Schur-Polynom via Quotientenformel
            val = (alpha ** (k + 1) - beta ** (k + 1)) / (alpha - beta)
        else:
            # Entarteter Fall α = β: Ableitung nach α
            val = (k + 1) * (alpha ** k)
        return val

    def hecke_algebra_generators(
        self,
        p: int
    ) -> Dict[str, Any]:
        """
        @brief Gibt die Generatoren der Hecke-Algebra H(GL_n(Q_p)//GL_n(Z_p)).
        @description
            Für GL_2(Q_p): Die Hecke-Algebra wird erzeugt durch:
            - T_p  = char. Funktion von GL_2(Z_p)·diag(p,1)·GL_2(Z_p)
            - T_{p²} = char. Funktion von GL_2(Z_p)·diag(p,p)·GL_2(Z_p)
              (= "Diamant-Operator" p·<p>)

            Satake-Bild (Isomorphismus nach ℂ[z₁±1, z₂±1]^{S_2}):
                S(T_p) = z_1 + z_2  (elementar-symmetrisch e_1)
                S(T_p·<p>) = z_1·z_2  (elementar-symmetrisch e_2)

            Für GL_n(Q_p): n Generatoren T_{p,k} für k=1,...,n,
                die auf elementare Symmetrische Polynome in (z_1,...,z_n) abbilden.

        @param p Primzahl
        @return Dict mit Generator-Namen und ihren Satake-Bildern
        """
        # Symbolische Variablen für Satake-Parameter
        z_vars = [sp.Symbol(f"z_{i+1}") for i in range(self.n)]

        # Elementare symmetrische Polynome als Satake-Bilder
        from sympy.functions.combinatorial.numbers import nC
        gens = {}
        for k in range(1, self.n + 1):
            # e_k(z_1,...,z_n) = Σ_{i_1<...<i_k} z_{i_1}·...·z_{i_k}
            from itertools import combinations
            e_k = sum(
                sp.Mul(*[z_vars[i] for i in combo])
                for combo in combinations(range(self.n), k)
            )
            gens[f"T_{{p^{k}}}"] = {
                "generator": f"GL_{self.n}(Z_p)·diag({','.join(['p']*k + ['1']*(self.n-k))})·GL_{self.n}(Z_p)",
                "satake_image": e_k,
                "satake_image_str": str(e_k),
            }

        return {
            "group": f"GL_{self.n}",
            "prime": p,
            "weyl_group": f"S_{self.n}",
            "generators": gens,
            "satake_ring": f"C[z_1^±1,...,z_{self.n}^±1]^S_{self.n}",
        }

    def satake_transform(
        self,
        f_values: Dict[int, complex],
        alpha: complex,
        beta: complex
    ) -> complex:
        """
        @brief Berechnet den Satake-Transform S(f) für GL_2.
        @description
            S(f)(z) = Σ_k f(p^k) · φ_{α,β}(p^k) · (Maß-Normierung)

            Der Satake-Transform ist ein Algebra-Homomorphismus:
                S: H(GL_2(Q_p)//GL_2(Z_p)) → ℂ[α±1, β±1]^{S_2}

        @param f_values Dict {k: f(p^k)} der Hecke-Funktion-Werte
        @param alpha erster Satake-Parameter
        @param beta zweiter Satake-Parameter
        @return Wert des Satake-Transforms
        """
        result = 0.0 + 0j
        for k, fval in f_values.items():
            # φ_z(p^k): sphärische Funktion
            phi = self.spherical_function(alpha, beta, p=2, k=k)
            result += fval * phi
        return result


# ---------------------------------------------------------------------------
# 2. Lokale Langlands-Korrespondenz für GL_1
# ---------------------------------------------------------------------------

class LocalLanglandsGL1:
    """
    @brief Lokale Langlands-Korrespondenz für GL_1: Klassenkörpertheorie.
    @description
        Für GL_1 ist die lokale Langlands-Korrespondenz äquivalent zur
        lokalen Klassenkörpertheorie (LCFT):

            Art_p: Q_p^× → Gal(Q_p^{ab} / Q_p)

        Genauer: Die lokale Artin-Abbildung liefert einen Isomorphismus
            Art_p: Q_p^× / N_{L/Q_p}(L^×) → Gal(L/Q_p)

        für abelsche Erweiterungen L/Q_p.

        Dualer Standpunkt (Langlands GL_1):
        - Ein Charakter χ: Q_p^× → ℂ^× (glatter Charakter)
        - Entspricht einer eindimensionalen Darstellung ρ_χ: W_{Q_p} → GL_1(ℂ)
        - Durch Art_p^{-1}: Gal → Q_p^× → ℂ^×

        Normrestsymbol (Artin, für lokales Körper):
            (a, L/K)_p ∈ Gal(L/K)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, p: int):
        """
        @brief Initialisiert GL_1-Korrespondenz für die Stelle p.
        @param p Primzahl (lokale Stelle)
        """
        self.p = p

    def class_field_theory_map(
        self,
        chi_values: Dict[int, complex],
        mod: int
    ) -> Dict[str, Any]:
        """
        @brief Lokale Klassenkörpertheorie-Abbildung χ → ρ_χ.
        @description
            Ein Dirichlet-Charakter χ mod q mit χ(-1) = ±1 entspricht
            einem Grössencharakter der Weil-Gruppe:
                χ: A_Q^× / Q^× → ℂ^×

            Lokal bei p: χ_p: Q_p^× → ℂ^×
            - Wenn p ∤ q: χ_p unramifiziert, χ_p(p) = χ(p)
            - Wenn p | q: χ_p ramifiziert (Level = v_p(q))

            Die Korrespondenz Art_p schickt:
                Frobenius Frob_p ↔ p ∈ Q_p^×
                Trägheitsgenerator I_p ↔ (1+p)^{1/p} ∈ Q_p^×

        @param chi_values Dict {n: χ(n)} der Charakter-Werte
        @param mod Modulus des Charakters
        @return Dict mit Korrespondenz-Daten
        """
        p = self.p

        # Ramifizierung prüfen: p | mod?
        is_ramified = (mod % p == 0)

        # Charakter-Wert bei Frobenius
        if not is_ramified and p in chi_values:
            chi_p_frobenius = chi_values[p]
            level = 0
        else:
            # Ramifizierter Fall: χ_p(p) nicht wohldefiniert (p im Kern der Norm)
            chi_p_frobenius = None
            level = 1

        # v_p-Valuation des Modulus
        q_temp = mod
        v_p_mod = 0
        while q_temp % p == 0:
            q_temp //= p
            v_p_mod += 1

        # Weil-Gruppe-Parameter: unramifiziert → 1-dim. Darstellung via Frobenius
        rho_frobenius = chi_p_frobenius if chi_p_frobenius is not None else 0.0 + 0j

        return {
            "prime": p,
            "modulus": mod,
            "is_ramified": is_ramified,
            "conductor_exponent": v_p_mod,
            "chi_p_frobenius": chi_p_frobenius,
            "rho_frobenius_eigenvalue": rho_frobenius,
            "artin_conductor": p ** v_p_mod if is_ramified else 1,
            "description": (
                f"χ_p: Q_{p}^× → ℂ^×, {'ramifiziert' if is_ramified else 'unramifiziert'}, "
                f"Frobenius-Eigenwert: {rho_frobenius}"
            ),
        }

    def artin_reciprocity_local(
        self,
        a: int,
        n: int
    ) -> Dict[str, Any]:
        """
        @brief Lokale Artin-Reziprozität: (a, Q_p(ζ_n)/Q_p).
        @description
            Das lokale Normrestsymbol für die zyklotomische Erweiterung Q_p(ζ_n)/Q_p:
                Art_p(a) ∈ Gal(Q_p(ζ_n)/Q_p) ≅ (Z/nZ)^×

            Für unramifizierte Stellen (p ∤ n):
                Art_p(p^k · u) = Frob_p^k, u ∈ Z_p^×
                Frob_p(ζ_n) = ζ_n^p  (Frobenius-Automorphismus)

            Lokale Artin-Abbildung:
                Art_p: Q_p^× → (Z/nZ)^×
                a = p^{v_p(a)} · u  ↦  p^{v_p(a)} mod n (vereinfacht)

        @param a ganzzahliges Element von Q_p^×
        @param n Ordnung der zyklotomischen Erweiterung
        @return Dict mit Artin-Symbol und Frobenius-Daten
        """
        p = self.p

        # p-adische Valuation von a
        v = 0
        temp = abs(a)
        if temp == 0:
            raise ValueError("a muss ungleich 0 sein (Q_p^×)")
        while temp % p == 0:
            temp //= p
            v += 1

        # Einheit u = a / p^v (mod p-Anteil)
        unit_part = a // (p ** v) if v > 0 else a

        # Frobenius-Exponent: v_p(a) mod φ(n)
        phi_n = totient(n)
        frobenius_exp = v % int(phi_n)

        # Artin-Symbol: Frob_p^{v_p(a)} mod n
        artin_symbol = pow(p, v, n)

        # Überprüfe Reziprozitäts-Konsistenz
        is_norm = (artin_symbol == 1)  # a ist Norm ↔ Art_p(a) = id

        return {
            "a": a,
            "n": n,
            "p": p,
            "v_p(a)": v,
            "unit_part": unit_part,
            "frobenius_exponent": frobenius_exp,
            "artin_symbol_mod_n": artin_symbol,
            "is_local_norm": is_norm,
            "frobenius_on_zeta": f"ζ_n ↦ ζ_n^{pow(p, v, n)}",
        }

    def local_norm_residue_symbol(
        self,
        a: int,
        p: int,
        n: int
    ) -> int:
        """
        @brief Berechnet das lokale Normrestsymbol (a, Q_p / p^n).
        @description
            Das Normrestsymbol (a, L/K)_p misst, ob a eine Norm ist:
            - (a, L/K)_p = 1 ↔ a ∈ N_{L/K}(L^×)
            - Für K = Q_p, L = Q_p(ζ_n): (a, Q_p(ζ_n)/Q_p)

            Berechnung via p-adische Potenzreihe (vereinfacht):
                (a, Q_p/p^n) = Legendre-Symbol-Verallgemeinerung
            Für p ungerade und a = p^{2k}·u, u ∈ Z_p^×:
                (a, Q_p(√p)/Q_p) = Legendre(u, p)

        @param a Ganzzahl (Element von Q_p^×)
        @param p Primzahl
        @param n Potenz (Erweiterungsordnung)
        @return Normrestsymbol als Integer (±1 oder n-te Einheitswurzel als int)
        """
        # p-adische Valuation
        v = 0
        temp = abs(a)
        while temp % p == 0 and temp > 0:
            temp //= p
            v += 1

        unit = temp  # |a| = p^v · unit, gcd(unit, p) = 1

        if n == 2:
            # Quadratisches Normrestsymbol (Hilbert-Symbol)
            # (a, p)_2 = Legendre-Symbol (vereinfacht für p ungerade)
            if p == 2:
                # Hilbert-Symbol bei p=2: komplexer Fall
                u_mod8 = unit % 8
                hilbert = 1 if u_mod8 in {1, 7} else -1
                if v % 2 == 1:
                    hilbert = -hilbert
                return hilbert
            else:
                # Hilbert-Symbol bei ungeradem p
                if v % 2 == 0:
                    # a = u · (p-adische Einheit)²·p^{2k}: Legendre(unit, p)
                    return int(sp.legendre_symbol(unit % p, p))
                else:
                    # v ungerade: (p, u)_p · Legendre(u, p) · Legendre(p, p)
                    leg_u = int(sp.legendre_symbol(unit % p, p))
                    return leg_u
        else:
            # Allgemeines n: Potenzrestsymbol (vereinfachte Rückgabe)
            return pow(int(unit), int(totient(p)) // n, p)


# ---------------------------------------------------------------------------
# 3. Lokale Langlands-Korrespondenz für GL_2
# ---------------------------------------------------------------------------

class LocalLanglandsGL2:
    """
    @brief Lokale Langlands-Korrespondenz für GL_2.
    @description
        Die lokale Langlands-Korrespondenz für GL_2(Q_p) ist ein Bijection:

            Irreduzible glatte Darstellungen π von GL_2(Q_p)
            ↔ 2-dimensionale Weil-Deligne-Darstellungen (ρ, N) von W_{Q_p}

        Die Darstellungstypen von GL_2(Q_p):
        1. Hauptreihen π(χ_1, χ_2) = Ind_B^G (χ_1 ⊗ χ_2), χ_1 ≠ χ_2·|·|^±1
           Weil-Deligne: ρ = χ_1 ⊕ χ_2 (direkte Summe), N = 0
        2. Spezielle Darstellung Sp(χ) = χ·StGL_2
           Weil-Deligne: ρ = χ·|·|^{1/2} ⊕ χ·|·|^{-1/2}, N ≠ 0
        3. Supercuspidal (keine Einbettung in Hauptreihe)
           Weil-Deligne: ρ irreduzibel, kommt von Induktion W_{Q_{p²}} → W_{Q_p}

        Lokale L-Faktor:
            L(s, π) = L(s, ρ) = det(1 - Frob_p·p^{-s} | ρ^{I_p})^{-1}
            Für Hauptreihe: L(s, π(χ_1,χ_2)) = L(s, χ_1)·L(s, χ_2)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, p: int):
        """
        @brief Initialisiert die GL_2-Korrespondenz für Stelle p.
        @param p Primzahl (lokale Stelle)
        """
        self.p = p

    def weil_deligne_to_admissible(
        self,
        rep_type: str,
        alpha: complex,
        beta: complex
    ) -> Dict[str, Any]:
        """
        @brief Weil-Deligne-Darstellung → irreduzible admissible Darstellung von GL_2.
        @description
            Klassifikation nach Typ der Weil-Deligne-Darstellung (ρ, N):

            - "principal_series": ρ = χ_1 ⊕ χ_2 reduzibel, N=0
              → π(χ_1, χ_2) Hauptreihen-Darstellung
            - "special": ρ = χ·|·|^{1/2} ⊕ χ·|·|^{-1/2}, N≠0 (Jordan-Block)
              → Sp(χ) = χ ⊗ StGL_2 (Steinberg-Darstellung)
            - "supercuspidal": ρ irreduzibel (kein inv. Unterraum)
              → Supercuspidal π_σ ohne Vorkommen in Parabolischer Ind.

        @param rep_type Typ der WD-Darstellung: "principal_series", "special", "supercuspidal"
        @param alpha erster Parameter (Frobenius-Eigenwert auf ρ, Komp. 1)
        @param beta zweiter Parameter (Frobenius-Eigenwert auf ρ, Komp. 2)
        @return Dict mit admissible Darstellungs-Daten
        """
        p = self.p

        if rep_type == "principal_series":
            # Hauptreihe: Ind_B^{GL_2}(χ_1 ⊗ χ_2) für χ_1(p)=α, χ_2(p)=β
            return {
                "type": "principal_series",
                "weil_deligne": {"rho": [alpha, beta], "N": 0},
                "admissible_rep": f"π(χ_1, χ_2) mit χ_1(p)={alpha:.4f}, χ_2(p)={beta:.4f}",
                "l_factor": f"(1 - {alpha:.4f}·p^{{-s}})^{{-1}}·(1 - {beta:.4f}·p^{{-s}})^{{-1}}",
                "conductor_exponent": 0,  # unramifiziert
                "is_generic": True,
            }
        elif rep_type == "special":
            # Steinberg: χ ⊗ St, χ(p) = α, N = [[0,1],[0,0]]
            chi_val = alpha  # χ(p)
            return {
                "type": "special",
                "weil_deligne": {
                    "rho": [chi_val * math.sqrt(p), chi_val / math.sqrt(p)],
                    "N": [[0, 1], [0, 0]]  # Nilpotenter Operator
                },
                "admissible_rep": f"Sp(χ) = χ ⊗ StGL_2 mit χ(p)={chi_val:.4f}",
                "l_factor": f"(1 - {chi_val:.4f}·p^{{-s}})^{{-1}}",
                "conductor_exponent": 1,  # Führer = p
                "is_generic": True,
            }
        elif rep_type == "supercuspidal":
            # Supercuspidal: Weil-Deligne ρ = Ind_{W_{Q_{p²}}}^{W_{Q_p}} σ
            return {
                "type": "supercuspidal",
                "weil_deligne": {
                    "rho": "Ind_{W_{Q_{p^2}}}^{W_{Q_p}} σ (2-dim irreduzibel)",
                    "N": 0
                },
                "admissible_rep": "π supercuspidal (aus Kompakter Induktion)",
                "l_factor": "1",  # L(s,π_sc) = 1 für tiefe Supercuspidals
                "conductor_exponent": 2,  # mind. p²
                "is_generic": True,
            }
        else:
            raise ValueError(f"Unbekannter rep_type: {rep_type}")

    def principal_series(
        self,
        chi1_p: complex,
        chi2_p: complex
    ) -> Dict[str, Any]:
        """
        @brief Berechnet Hauptreihen-Darstellung π(χ_1, χ_2) für GL_2(Q_p).
        @description
            π(χ_1, χ_2) = Ind_B^{GL_2}(χ_1 ⊗ χ_2):
            - B = obere Borel-Gruppe
            - Irreduzibel wenn χ_1/χ_2 ≠ |·|^±1
            - Zentralcharakter ω = χ_1·χ_2
            - Satake-Parameter: z_1 = χ_1(p), z_2 = χ_2(p)
            - Hecke-Eigenwert: a_p = p^{(k-1)/2}·(z_1 + z_2)

        @param chi1_p Wert χ_1(p) (Charakter-Wert bei Frobenius)
        @param chi2_p Wert χ_2(p)
        @return Dict mit Darstellungs-Eigenschaften
        """
        p = self.p

        # Hecke-Eigenwert T_p: a_p = chi1(p) + chi2(p) (normiert Gewicht 2)
        a_p = chi1_p + chi2_p

        # Zentralcharakter: ω(p) = χ_1(p)·χ_2(p)
        omega_p = chi1_p * chi2_p

        # Irreduzibilitätsbedingung: χ_1/χ_2 ≠ |p|^±1
        ratio = chi1_p / chi2_p if abs(chi2_p) > 1e-12 else float('inf')
        is_irreducible = (
            abs(abs(ratio) - p) > 1e-8 and
            abs(abs(ratio) - 1.0 / p) > 1e-8
        )

        return {
            "representation": "principal_series",
            "chi1_p": chi1_p,
            "chi2_p": chi2_p,
            "hecke_eigenvalue_a_p": a_p,
            "central_character_omega_p": omega_p,
            "is_irreducible": is_irreducible,
            "satake_parameters": (chi1_p, chi2_p),
            "l_factor_inverse": f"(1-{chi1_p:.4f}·X)(1-{chi2_p:.4f}·X)",
        }

    def special_representation(
        self,
        chi_p: complex
    ) -> Dict[str, Any]:
        """
        @brief Steinberg/spezielle Darstellung Sp(χ) = χ ⊗ StGL_2.
        @description
            Die Steinberg-Darstellung StGL_2 ist die einzige irreduzible
            Quotient-Darstellung von π(|·|^{1/2}, |·|^{-1/2}).
            Sp(χ) = χ ⊗ StGL_2 für einen Charakter χ.

            Entspricht Weil-Deligne:
                ρ = χ·|·|^{1/2} ⊕ χ·|·|^{-1/2}, N = [[0,1],[0,0]]

            L-Faktor: L(s, Sp(χ)) = L(s + 1/2, χ)
            ε-Faktor: ε(s, Sp(χ)) = ε(s+1/2, χ)·ε(s-1/2, χ)^{-1} · χ(-1)

        @param chi_p Wert des Charakters χ bei p
        @return Dict mit Sp(χ)-Eigenschaften
        """
        p = self.p
        sqrt_p = math.sqrt(p)

        # Weil-Deligne: ρ-Eigenwerte
        rho1 = chi_p * sqrt_p      # χ(p)·|p|^{-1/2} = χ(p)·√p
        rho2 = chi_p / sqrt_p      # χ(p)·|p|^{1/2} = χ(p)/√p

        return {
            "representation": "special",
            "chi_p": chi_p,
            "weil_deligne_eigenvalues": (rho1, rho2),
            "monodromy_N": [[0, 1], [0, 0]],
            "l_factor": f"(1 - {chi_p:.4f}·p^{{-s-1/2}})^{{-1}}",
            "conductor_exponent": 1,
            "epsilon_factor": f"ε(s+1/2, χ) für χ(p)={chi_p:.4f}",
        }

    def supercuspidal_marker(self) -> Dict[str, Any]:
        """
        @brief Markiert eine supercuspidale Darstellung von GL_2(Q_p).
        @description
            Supercuspidale Darstellungen π sind definiert durch:
            - Kein nicht-triviales Jacquet-Modul: J_B(π) = 0
            - Keine Einbettung in Parabolisch-Induzierte Darstellungen
            - Kommen aus kompakter Induktion: π = c-Ind_{Z·GL_2(Z_p)}^{GL_2(Q_p)} τ

            Korrespondenz: Supercuspidal π_σ ↔ irreduzible WD-Darstellung ρ_σ
            (durch lokale Langlands-Korrespondenz eindeutig)

            L-Faktor: L(s, π_sc) = 1 (für "tiefe" Supercuspidals)
            ε-Faktor: ε(s, π_sc) = ε(s, ρ_σ) (nach WD-Daten)

        @return Dict mit Klassifikations-Markierungen
        """
        return {
            "representation": "supercuspidal",
            "jacquet_module": "trivial (J_B(π) = 0)",
            "compact_induction": "c-Ind_{Z·K}^G τ für irred. K-Darst. τ",
            "weil_deligne_type": "irreduzibel 2-dimensional",
            "l_factor": "1",
            "conductor_exponent": "≥ 2",
            "depth": "positiv (tiefe Darstellung)",
        }

    def local_factor_l(
        self,
        pi_type: str,
        s: complex,
        alpha: complex = 1.0 + 0j,
        beta: complex = 1.0 + 0j
    ) -> complex:
        """
        @brief Berechnet den lokalen L-Faktor L(s, π_p).
        @description
            L(s, π_p) hängt vom Typ der Darstellung ab:
            - Hauptreihe π(χ_1, χ_2): L = (1-α·p^{-s})^{-1}·(1-β·p^{-s})^{-1}
            - Spezielle Sp(χ):         L = (1-α·p^{-s})^{-1}  (α = χ(p))
            - Supercuspidal:           L = 1

            Für elliptische Kurven E/Q an guter Stelle p:
                L(s, E_p) = (1 - a_p·p^{-s} + p^{1-2s})^{-1}
                = (1 - α_p·p^{-s})^{-1}·(1 - β_p·p^{-s})^{-1}

        @param pi_type Darstellungstyp: "principal_series", "special", "supercuspidal"
        @param s komplexe Variable
        @param alpha erster Satake-Parameter
        @param beta zweiter Satake-Parameter
        @return Wert des lokalen L-Faktors
        """
        p = self.p
        p_s = p ** (-s)  # p^{-s}

        if pi_type == "principal_series":
            # L(s, π) = (1-α·p^{-s})^{-1}·(1-β·p^{-s})^{-1}
            l_val = 1.0 / ((1 - alpha * p_s) * (1 - beta * p_s))
        elif pi_type == "special":
            # L(s, Sp(χ)) = (1-α·p^{-s})^{-1}
            l_val = 1.0 / (1 - alpha * p_s)
        elif pi_type == "supercuspidal":
            # L(s, π_sc) = 1
            l_val = 1.0 + 0j
        else:
            raise ValueError(f"Unbekannter Typ: {pi_type}")
        return l_val

    def local_epsilon_factor(
        self,
        pi_type: str,
        psi_level: int = 0
    ) -> complex:
        """
        @brief Berechnet den lokalen ε-Faktor ε(1/2, π_p, ψ).
        @description
            Der ε-Faktor ist das "Root Number" der lokalen funktionalen Gleichung:
                L(s, π) = ε(s, π)·L(1-s, π̌)

            Bei s=1/2: ε(1/2, π_p) ∈ {±1, ±i} (bei unitären zentralen Charakter)

            Formeln:
            - Hauptreihe (unramifiziert): ε(1/2, π) = 1
            - Spezielle Darstellung:       ε(1/2, Sp(χ)) = χ(-1) = ±1
            - Supercuspidal Niveau n:       ε(1/2, π) = Gauß-Summe · p^{-n/2}
              (komplex, |ε|=1)
            - ψ-Abhängigkeit: ε(s, π, ψ·a) = χ(a)·|a|^{...}·ε(s, π, ψ)

        @param pi_type Darstellungstyp
        @param psi_level Niveau des additiven Charakters ψ
        @return Wert des ε-Faktors bei s=1/2
        """
        p = self.p

        if pi_type == "principal_series":
            # Unramifizierte Hauptreihe: ε = 1
            eps = 1.0 + 0j
        elif pi_type == "special":
            # Sp(χ): ε(1/2, Sp) = χ(-1) = ±1, hier χ_triv → +1
            eps = 1.0 + 0j
        elif pi_type == "supercuspidal":
            # Supercuspidal Niveau n ≥ 1: ε = Gauß-Summe / p^{n/2}
            # Beispiel: Niveau 1 Gauß-Summe ≈ √p · e^{iπ/4}
            n = max(1, psi_level)
            gauss_approx = math.sqrt(p) * cmath.exp(1j * math.pi / 4)
            eps = gauss_approx / (p ** (n / 2.0))
        else:
            raise ValueError(f"Unbekannter Typ: {pi_type}")
        return eps


# ---------------------------------------------------------------------------
# 4. Globale Langlands-Korrespondenz
# ---------------------------------------------------------------------------

class GlobalLanglands:
    """
    @brief Globale Langlands-Korrespondenz: Modularitätssatz und Galois-Darstellungen.
    @description
        Die globale Langlands-Korrespondenz für GL_2/Q verbindet:

        (A) Galois-Darstellungen ρ: Gal(Q̄/Q) → GL_2(ℂ)
            (irreduzibel, unverweigt außerhalb endl. vieler Primes)

        (B) Automorphe Darstellungen π = ⊗'_p π_p auf GL_2(A_Q)
            (cuspidal, unitär, π_∞ = holomorphe Darst. vom Gewicht k)

        Modularitätssatz (Wiles 1995, Taylor-Wiles):
            Jede elliptische Kurve E/Q ist modular:
            ∃ f ∈ S_2(Γ_0(N)) mit L(E, s) = L(f, s)
            Dabei N = Führer von E.

        Globale Funktionalgleichung:
            Λ(s, π) = ε(π) · Λ(1-s, π̌)
            Λ(s, π) = N^{s/2} · Γ_ℂ(s)^d · L(s, π)
            Γ_ℂ(s) = 2·(2π)^{-s}·Γ(s)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def modular_curve_wiles_check(
        self,
        a_coeff: Dict[int, int],
        N: int
    ) -> Dict[str, Any]:
        """
        @brief Überprüft Wiles' Modularitätsbedingung für eine elliptische Kurve.
        @description
            Wiles' Theorem (1995): Jede semistabile elliptische Kurve E/Q ist modular.
            Später: Breuil-Conrad-Diamond-Taylor (2001): Alle E/Q sind modular.

            Modularitätsbedingung:
            - a_p(E) = a_p(f) für fast alle Primzahlen p
            - f ∈ S_2(Γ_0(N)) ist eine Hecke-Eigenform
            - L(E, s) = Σ a_n(E) n^{-s} = L(f, s)

            Frobenius-Spur für E: y² = x³ - x (Führer N=32):
                a_p = 0           für p=2, oder p≡3(4)
                a_p = 2·Re(π_p) für p≡1(4) (mit π_p∈Z[i], π_p·π̄_p=p)

        @param a_coeff Dict {p: a_p} der Frobenius-Spuren
        @param N Führer der elliptischen Kurve
        @return Dict mit Modularitäts-Checks
        """
        checks = {}
        all_ok = True

        for p, a_p in a_coeff.items():
            if not isprime(p):
                continue

            # Ramanujan-Schranke: |a_p| ≤ 2√p
            ramanujan_bound = 2.0 * math.sqrt(p)
            satisfies_ramanujan = abs(a_p) <= ramanujan_bound + 1e-8

            # Für Gewicht 2: a_p·p^{-1/2} auf Einheitskreis? (Ramanujan)
            if p != 2 and N % p != 0:
                # Gute Reduktion: Char.-Polynom X²-a_p·X+p
                disc = a_p ** 2 - 4 * p
                is_ordinary = disc >= 0  # reelle Satake-Parameter → gewöhnlich
                is_supersingular = (a_p % p == 0)  # a_p ≡ 0 mod p
            else:
                is_ordinary = None
                is_supersingular = None

            checks[p] = {
                "a_p": a_p,
                "ramanujan_bound": ramanujan_bound,
                "satisfies_ramanujan": satisfies_ramanujan,
                "is_ordinary": is_ordinary,
                "is_supersingular": is_supersingular,
                "good_reduction": (N % p != 0),
            }
            if not satisfies_ramanujan:
                all_ok = False

        return {
            "conductor": N,
            "modular_form_weight": 2,
            "modularity_theorem": "Breuil-Conrad-Diamond-Taylor (2001)",
            "prime_checks": checks,
            "all_ramanujan_satisfied": all_ok,
            "conclusion": (
                f"E mit Führer N={N} ist modular (BCDT 2001). "
                f"L(E,s) = L(f,s) für f ∈ S_2(Γ_0({N}))."
            ),
        }

    def modularity_lifting_theorem_demo(
        self,
        E_conductor: int,
        E_ap: Dict[int, int]
    ) -> Dict[str, Any]:
        """
        @brief Demo des Modularitäts-Lifting-Theorems (Wiles-Taylor).
        @description
            Modularitäts-Lifting:
            1. Residuelle Galois-Darstellung: ρ̄_E: Gal(Q̄/Q) → GL_2(F_ℓ)
               (Tate-Modul mod ℓ)
            2. Wenn ρ̄_E modular (Modulär mod ℓ)
               UND ρ̄_E|_{G_{Q_ℓ}} geeignet (flach oder Steinberg)
            3. DANN: ρ_E (ℓ-adisch) ist modular
            4. Also: E modular.

            Für ℓ=3: ρ̄_E,3 irreduzibel → Serre's Modularity Conj. (Khare-Wintenberger)
            Für ℓ=5: Fallback wenn ρ̄_E,3 reduzibel (Wiles 1995 Trick: ρ̄_E,5 betrachten)

            Beweis-Strategie Wiles:
            - R = T (Ring-Isomorphismus Galois-Ring ≅ Hecke-Algebra)
            - Gelingt über Euler-System (Kolyvagin) und Selmer-Gruppen

        @param E_conductor Führer der elliptischen Kurve
        @param E_ap Dict {p: a_p} der Frobenius-Spuren
        @return Dict mit Lifting-Theorem-Daten
        """
        # Überprüfe residuelle Darstellung mod 3 und mod 5
        residues = {}
        for ell in [3, 5, 7]:
            residues[ell] = {}
            for p, ap in E_ap.items():
                # Spur mod ℓ
                trace_mod = ap % ell
                # Determinante = p (Zyklotomischer Charakter)
                det_mod = p % ell
                residues[ell][p] = {
                    "trace_mod_ell": trace_mod,
                    "det_mod_ell": det_mod,
                }

        # Hecke-Algebra T_N (symbolisch)
        T_primes = sorted([p for p in E_ap.keys() if isprime(p)])[:5]
        hecke_data = {p: E_ap[p] for p in T_primes}

        return {
            "conductor": E_conductor,
            "theorem": "Modularitäts-Lifting (Wiles 1995 + BCDT 2001)",
            "strategy": {
                "step1": "ρ̄_{E,3}: Gal(Q̄/Q) → GL_2(F_3) (Residualdarst.)",
                "step2": "ρ̄_{E,3} irreduzibel → modular (Langlands-Tunnell)",
                "step3": "Lifting: R=T Isomorphismus (Wiles Hauptsatz)",
                "step4": "∴ ρ_E,3-adisch modular → E modular",
                "fallback": "Falls ρ̄_{E,3} reduzibel: ρ̄_{E,5} stattdessen (Wiles' Trick)",
            },
            "residual_mod_ell": residues,
            "hecke_eigenvalues": hecke_data,
            "R_equals_T": "R ≅ T (Galois-Deformationsring ≅ Hecke-Algebra)",
            "conclusion": f"E/Q mit N={E_conductor} ist modular.",
        }

    def langlands_conjecture_evidence(
        self,
        rep_dim: int
    ) -> Dict[str, Any]:
        """
        @brief Sammelt Belege für die Langlands-Vermutung für GL_n.
        @description
            Stand der Langlands-Vermutung für verschiedene Dimensionen:

            GL_1: Bewiesen (Klassenkörpertheorie, Artin 1927)
            GL_2: Bewiesen für Q (Wiles 1995 für E, allgemein via Langlands, Tunnell, Jacquet-Langlands)
            GL_3: Partiell (Kim-Shahidi: sym²-Lift von GL_2)
            GL_4: sym³-Lift von GL_2 (Kim 2003)
            GL_n allgemein: Offen (V. Lafforgue für Funktionskörper: bewiesen 2018)

            Über Funktionskörpern F_q(t) (Lafforgue, Fields-Medal 2002):
                Vollständige Langlands-Korrespondenz GL_n bewiesen!

        @param rep_dim Dimension n der GL_n-Darstellung
        @return Dict mit Beleg-Zusammenfassung
        """
        evidence = {
            1: {
                "status": "Bewiesen",
                "theorem": "Lokale und globale Klassenkörpertheorie",
                "year": 1927,
                "author": "E. Artin",
            },
            2: {
                "status": "Bewiesen (für ℚ)",
                "theorem": "Modularitätssatz (Wiles) + Fontaine-Mazur + Langlands-Tunnell",
                "year": 2001,
                "author": "Breuil, Conrad, Diamond, Taylor (BCDT)",
                "partial": "sym²: Kim-Shahidi 2000, sym³: Kim 2003",
            },
            3: {
                "status": "Partiell",
                "theorem": "Induktive Lifts von GL_2 via Sym-Funktorialität",
                "year": 2003,
                "author": "Kim-Shahidi",
            },
            4: {
                "status": "Partiell",
                "theorem": "Sym³-Lift von GL_2 → GL_4 (Kim 2003)",
                "year": 2003,
                "author": "H. Kim",
            },
        }

        if rep_dim not in evidence:
            result = {
                "status": "Offen (für Zahlkörper)",
                "theorem": "V. Lafforgue (2018): Vollständig über Funktionskörpern",
                "year": 2018,
                "note": f"Für GL_{rep_dim}/Q: aktives Forschungsgebiet",
            }
        else:
            result = evidence[rep_dim]

        return {
            "dimension": rep_dim,
            "group": f"GL_{rep_dim}",
            "evidence": result,
            "function_field": "V. Lafforgue (2018): GL_n über F_q(t) vollständig bewiesen",
            "p_adic": "Colmez-Fontaine: p-adische Langlands für GL_2(Q_p) bewiesen",
        }


# ---------------------------------------------------------------------------
# 5. Langlands-Funktorialität
# ---------------------------------------------------------------------------

class LanglandsFunctoriality:
    """
    @brief Langlands-Funktorialität: Lifts zwischen automorphen Darstellungen.
    @description
        Langlands-Funktorialitätsvermutung:
        Für einen L-Gruppen-Homomorphismus r: ^H → ^G
        existiert ein Lift F_r: Π(H) → Π(G) der automorphen Darstellungen:
            π_H automorph auf H(A) → F_r(π_H) automorph auf G(A)
            mit L(s, F_r(π_H)) = L(s, π_H, r)

        Bewiesene Fälle:
        1. Basiswechsel (Base Change) für GL_n (Arthur-Clozel 1989)
        2. sym²-Lift GL_2 → GL_3 (Gelbart-Jacquet 1978)
        3. sym³-Lift GL_2 → GL_4 (Kim-Shahidi 2002)
        4. sym⁴-Lift GL_2 → GL_5 (Kim 2003)
        5. Tensor-Produkt GL_2 × GL_2 → GL_4 (Ramakrishnan 2000)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def base_change_GL2(
        self,
        a_p: Dict[int, int],
        p: int,
        n: int
    ) -> Dict[str, Any]:
        """
        @brief Basiswechsel (Base Change) für GL_2: f/Q → f_K für [K:Q]=n.
        @description
            Für eine Modulform f mit Hecke-Eigenwerten a_p (p unramifiziert in K/Q):
            Basiswechsel BC_{K/Q}(π_f) ist eine automorphe Darstellung auf GL_2(A_K).

            Lokale Formel (unramifizierte Stelle q über p mit Trägheitsgrad f_q):
                a_q(BC(π)) = Spur des Frobenius Frob_q auf (α_p^{f_q}, β_p^{f_q})
                = α_p^{f_q} + β_p^{f_q}

            Für zyklische Erweiterung K/Q vom Grad n:
            - Primzahl p zerlegt (f_q=1): a_q = a_p (Kopie)
            - Primzahl p träge (f_q=n): a_q = α_p^n + β_p^n (Potenzsumme)
            - Primzahl p voll verzweigt: spezielle Behandlung (Steinberg)

        @param a_p Dict {prime: a_p} Hecke-Eigenwerte der Ausgangsform
        @param p konkrete Primzahl (Demonstrationsfall)
        @param n Grad der Erweiterung [K:Q]
        @return Dict mit Basiswechsel-Daten
        """
        results = {}
        for prime, ap in a_p.items():
            # Satake-Parameter (unitäre Normierung Gewicht 2)
            disc = float(ap) ** 2 - 4.0 * prime
            if disc < 0:
                alpha = complex(ap / 2.0, math.sqrt(-disc) / 2.0)
                beta = complex(ap / 2.0, -math.sqrt(-disc) / 2.0)
            else:
                alpha = complex((ap + math.sqrt(abs(disc))) / 2.0, 0)
                beta = complex((ap - math.sqrt(abs(disc))) / 2.0, 0)

            # Potenzsummen: α^n + β^n
            # Newton-Identitäten: p_n = a_p·p_{n-1} - p·p_{n-2}
            power_sums = [0.0 + 0j] * (n + 1)
            power_sums[0] = 2.0 + 0j   # p_0 = Rang = 2
            power_sums[1] = ap + 0j    # p_1 = a_p

            for k in range(2, n + 1):
                power_sums[k] = ap * power_sums[k - 1] - prime * power_sums[k - 2]

            results[prime] = {
                "original_a_p": ap,
                "alpha": alpha,
                "beta": beta,
                "base_change_a_q": power_sums[n],  # α^n + β^n
                "power_sums": [ps for ps in power_sums],
                "split_case": complex(alpha + beta),  # f_q=1 → Kopie
                "inert_case": power_sums[n],           # f_q=n → α^n+β^n
            }

        return {
            "theorem": f"Basiswechsel BC_{{K/Q}} für [K:Q]={n} (Arthur-Clozel 1989)",
            "prime_p": p,
            "degree_n": n,
            "results": results,
        }

    def symmetric_power_lift(
        self,
        a_p: Dict[int, int],
        p: int,
        k: int
    ) -> Dict[str, Any]:
        """
        @brief Symmetrische Potenz sym^k-Lift GL_2 → GL_{k+1}.
        @description
            Sym^k-Funktorialität: π_f auf GL_2 → Sym^k(π_f) auf GL_{k+1}

            Lokale Formel an guter Stelle p mit Satake-Parametern (α, β):
                Sym^k-Satake-Parameter: {α^k, α^{k-1}β, ..., β^k}
                (k+1 Stück, symmetrische k-te Potenzen)

            Hecke-Eigenwert von Sym^k(π) bei p:
                a_p(Sym^k) = Σ_{j=0}^{k} α^j · β^{k-j} = φ_{α,β}(p^k)
                           = (α^{k+1} - β^{k+1}) / (α - β)  (Schur-Polynom)

            Bewiesene Lifts (Gelbart-Jacquet, Kim-Shahidi, Kim):
                k=1: Identität (GL_2 → GL_2)
                k=2: Sym² (GL_2 → GL_3), Gelbart-Jacquet 1978
                k=3: Sym³ (GL_2 → GL_4), Kim-Shahidi 2002
                k=4: Sym⁴ (GL_2 → GL_5), Kim 2003
                k≥5: Offen

        @param a_p Dict {prime: a_p} Hecke-Eigenwerte
        @param p Primzahl (Beispiel)
        @param k Potenz (sym^k)
        @return Dict mit sym^k-Lift-Daten
        """
        results = {}
        for prime, ap in a_p.items():
            # Satake-Parameter
            disc = float(ap) ** 2 - 4.0 * prime
            if disc < 0:
                alpha = complex(ap / 2.0, math.sqrt(-disc) / 2.0)
                beta = complex(ap / 2.0, -math.sqrt(-disc) / 2.0)
            else:
                alpha = complex((ap + math.sqrt(abs(disc))) / 2.0, 0)
                beta = complex((ap - math.sqrt(abs(disc))) / 2.0, 0)

            # Sym^k Satake-Parameter: {α^k, α^{k-1}·β, ..., β^k}
            sym_params = [alpha ** (k - j) * beta ** j for j in range(k + 1)]

            # Sym^k L-Faktor-Eigenwert (Schur-Polynom)
            if abs(alpha - beta) > 1e-12:
                sym_k_eigenvalue = (alpha ** (k + 1) - beta ** (k + 1)) / (alpha - beta)
            else:
                sym_k_eigenvalue = (k + 1) * (alpha ** k)

            results[prime] = {
                "original_a_p": ap,
                "sym_k_satake_params": sym_params,
                "sym_k_eigenvalue": sym_k_eigenvalue,
                "sym_k_group": f"GL_{k + 1}",
            }

        status_map = {
            1: "trivial",
            2: "bewiesen (Gelbart-Jacquet 1978)",
            3: "bewiesen (Kim-Shahidi 2002)",
            4: "bewiesen (Kim 2003)",
        }
        status = status_map.get(k, f"offen für k={k}")

        return {
            "lift": f"Sym^{k}: GL_2 → GL_{k + 1}",
            "status": status,
            "prime_p": p,
            "results": results,
        }

    def automorphic_induction(
        self,
        chi_values: Dict[int, complex],
        p: int,
        K_degree: int
    ) -> Dict[str, Any]:
        """
        @brief Automorphe Induktion AI_{K/Q}: Π(GL_1/K) → Π(GL_n/Q).
        @description
            Automorphe Induktion ist die Umkehrung des Basiswechsels:
                BC_{K/Q}: Π(GL_n/Q) → Π(GL_n/K)
                AI_{K/Q}: Π(GL_1/K) → Π(GL_n/Q)  (n = [K:Q])

            Für einen Hecke-Grössencharakter χ: A_K^× / K^× → ℂ^×
            liefert AI_{K/Q}(χ) eine cuspidal automorphe Darstellung auf GL_n(A_Q).

            L-Funktion: L(s, AI(χ)) = L(s, χ) (als GL_1/K-L-Funktion)
            Bewiesen: K/Q abelsch (Langlands-Jacquet), allgemein offen.

        @param chi_values Dict {p: χ(p)} der Grössencharakter-Werte bei p
        @param p Primzahl (Demonstrationsfall)
        @param K_degree Grad [K:Q] = n (Dimension des Lifts)
        @return Dict mit AI-Lift-Daten
        """
        n = K_degree

        # Induzierende Galois-Darstellung: Ind_{Gal(K̄/K)}^{Gal(Q̄/Q)} χ
        if n == 2:
            # Quadratische Erweiterung K/Q: χ = Grössencharakter von K
            # AI(χ) = π_K auf GL_2/Q
            local_param = chi_values.get(p, 1.0 + 0j)
            ind_eigenvalues = [local_param, local_param.conjugate()]
            a_p_induced = local_param + local_param.conjugate()  # Spur = 2·Re
        else:
            # Allgemeiner Grad: n Eigenwerte
            local_param = chi_values.get(p, 1.0 + 0j)
            # Galois-Konjugierte (n Stück für zyklische Erw.)
            ind_eigenvalues = [
                local_param * cmath.exp(2j * math.pi * k / n)
                for k in range(n)
            ]
            a_p_induced = sum(ind_eigenvalues)

        return {
            "lift": f"AI_{{K/Q}}: GL_1/K → GL_{n}/Q",
            "K_degree": n,
            "prime_p": p,
            "chi_p": chi_values.get(p, "nicht angegeben"),
            "induced_satake_params": ind_eigenvalues,
            "a_p_induced": a_p_induced,
            "l_function_identity": f"L(s, AI(χ)) = L(s, χ/K)",
            "status": "bewiesen für K/Q abelsch (Arthur-Clozel, Langlands-Jacquet)",
        }

    def langlands_dual_group(
        self,
        group_type: str
    ) -> Dict[str, Any]:
        """
        @brief Bestimmt die Langlands-Dualgruppe ^G für eine reduktive Gruppe G.
        @description
            Die Langlands-Dualgruppe ^G ist definiert durch Dualisieren der
            Wurzeldaten: (X*, Φ, X^*, Φ^∨) → (X^*, Φ^∨, X*, Φ)

            Tabelle der Dualgruppen:
            G = GL_n       → ^G = GL_n(ℂ)
            G = SL_n       → ^G = PGL_n(ℂ) = GL_n(ℂ)/Z
            G = Sp_{2n}    → ^G = SO_{2n+1}(ℂ) (Typ C_n ↔ B_n)
            G = SO_{2n+1}  → ^G = Sp_{2n}(ℂ)   (Typ B_n ↔ C_n)
            G = SO_{2n}    → ^G = SO_{2n}(ℂ)   (Typ D_n ↔ D_n)
            G = G_2        → ^G = G_2(ℂ)        (selbstdual)
            G = F_4        → ^G = F_4(ℂ)        (selbstdual)
            G = E_6        → ^G = E_6(ℂ) (oder ^E_6)
            G = E_7        → ^G = E_7(ℂ)
            G = E_8        → ^G = E_8(ℂ)        (selbstdual)

            L-Gruppe: ^L G = ^G ⋊ Gal(Q̄/Q) (semi-direktes Produkt)

        @param group_type Gruppentyp als String (z.B. "GL_n", "Sp_2n", "SO_2n")
        @return Dict mit Dualgruppe und L-Gruppe
        """
        dual_map = {
            "GL_n":    {"dual": "GL_n(C)",    "type": "A_{n-1}", "self_dual": True},
            "SL_n":    {"dual": "PGL_n(C)",   "type": "A_{n-1}", "self_dual": False},
            "PGL_n":   {"dual": "SL_n(C)",    "type": "A_{n-1}", "self_dual": False},
            "Sp_2n":   {"dual": "SO_{2n+1}(C)","type": "C_n→B_n","self_dual": False},
            "SO_2n+1": {"dual": "Sp_{2n}(C)", "type": "B_n→C_n", "self_dual": False},
            "SO_2n":   {"dual": "SO_{2n}(C)", "type": "D_n→D_n", "self_dual": True},
            "G_2":     {"dual": "G_2(C)",     "type": "G_2",     "self_dual": True},
            "F_4":     {"dual": "F_4(C)",     "type": "F_4",     "self_dual": True},
            "E_6":     {"dual": "E_6sc(C)",   "type": "E_6",     "self_dual": False},
            "E_7":     {"dual": "E_7(C)",     "type": "E_7",     "self_dual": True},
            "E_8":     {"dual": "E_8(C)",     "type": "E_8",     "self_dual": True},
        }

        info = dual_map.get(group_type, {
            "dual": f"^{group_type}",
            "type": "unbekannt",
            "self_dual": None,
        })

        return {
            "group": group_type,
            "dual_group": info["dual"],
            "root_system_duality": info["type"],
            "is_self_dual": info["self_dual"],
            "L_group": f"^L {group_type} = {info['dual']} ⋊ Gal(Q̄/Q)",
            "functoriality": (
                f"r: ^{group_type} → ^G induziert F_r: Π({group_type}) → Π(G)"
            ),
        }


# ---------------------------------------------------------------------------
# 6. Artin-Reziprozität
# ---------------------------------------------------------------------------

class ArtinReciprocity:
    """
    @brief Artin-Reziprozität und Artin-L-Funktionen.
    @description
        Die Artin-Reziprozität (globale Klassenkörpertheorie) besagt:

        Für eine Galois-Erweiterung K/Q mit Gruppe G=Gal(K/Q):
        ∃ kanonischer Isomorphismus (globale Artin-Abbildung):
            Art_{K/Q}: A_Q^× / Q^× → G^{ab}

        Artin-L-Funktion für Darstellung ρ: G → GL_n(ℂ):
            L(s, ρ) = Π_p det(1 - Frob_p · p^{-s} | V^{I_p})^{-1}
            = Π_p L_p(s, ρ)

        Artin-Conductor-Formel:
            f(ρ) = Π_p p^{f_p(ρ)}
            f_p(ρ) = Σ_{i=0}^∞ |G_i|/|G_0| · codim(V^{G_i})
            (G_i = i-te Trägheitsgruppe)

        Artin-Satz: L(s, ρ) = Π_χ L(s, χ)^{m_χ} (Artin-Induktionssatz)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def artin_l_function(
        self,
        chi_values: Dict[int, complex],
        s: complex,
        primes: List[int]
    ) -> complex:
        """
        @brief Berechnet die Artin-L-Funktion L(s, χ) als Euler-Produkt.
        @description
            L(s, χ) = Π_p (1 - χ(Frob_p) · p^{-s})^{-1}
            für einen Dirichlet-Charakter χ oder eine 1-dimensionale Galois-Darst.

            Für den trivialen Charakter: L(s, 1) = ζ(s) (Riemann-Zeta)
            Für χ Dirichlet-Charakter mod q: L(s, χ) = Dirichlet-L-Funktion

            Konvergenz: absolut für Re(s) > 1.

        @param chi_values Dict {p: χ(Frob_p)} Charakter-Werte an Primzahlen
        @param s komplexe Variable (Re(s) > 1 für Konvergenz)
        @param primes Liste der Primzahlen für Euler-Produkt-Approximation
        @return Approximation von L(s, χ) via endliches Euler-Produkt
        """
        l_value = 1.0 + 0j
        for p in primes:
            if p not in chi_values:
                # Unramifizierter Fall ohne expliziten Wert: überspringen
                continue
            chi_frob = chi_values[p]
            p_s = p ** (-s)
            # Euler-Faktor: (1 - χ(Frob_p) · p^{-s})^{-1}
            factor = 1.0 / (1.0 - chi_frob * p_s)
            l_value *= factor
        return l_value

    def artin_conductor_formula(
        self,
        ramified_primes: Dict[int, Dict[str, Any]],
        chi_order: int
    ) -> int:
        """
        @brief Berechnet den Artin-Führer f(ρ) via Trägheitsgruppen-Formel.
        @description
            Artin-Führer (lokale Beiträge):
                f_p(ρ) = dim(V/V^{I_p}) + Swan-Exponent
                Swan(ρ, p) = Σ_{i≥1} |G_i|/|G_0| · codim(V^{G_i})

            Für tamm-ramifizierte Stelle p (G_1 = 1):
                f_p(ρ) = codim(V^{I_p}) = n - dim(V^{I_p})

            Für wild-ramifizierte Stelle p (G_1 ≠ 1):
                f_p(ρ) ≥ 2 (Swan-Beitrag ≥ 1)

            Globaler Artin-Führer:
                f(ρ) = Π_p p^{f_p(ρ)}

        @param ramified_primes Dict {p: {inertia_fixed_dim, wild_part}} Trägheitsdaten
        @param chi_order Ordnung des Charakters (Grad der Erweiterung)
        @return Artin-Führer als Integer
        """
        conductor = 1
        for p, data in ramified_primes.items():
            inertia_fixed_dim = data.get("inertia_fixed_dim", 0)
            wild_part = data.get("wild_part", 0)  # Swan-Exponent

            # codim(V^{I_p}) = n - inertia_fixed_dim (n=chi_order für 1-dim)
            tame_exp = chi_order - inertia_fixed_dim
            total_exp = tame_exp + wild_part

            conductor *= p ** total_exp

        return conductor

    def global_artin_map_demo(
        self,
        n: int
    ) -> Dict[str, Any]:
        """
        @brief Demonstration der globalen Artin-Abbildung für Q(ζ_n)/Q.
        @description
            Die zyklotomische Erweiterung Q(ζ_n)/Q ist abelsch mit:
                Gal(Q(ζ_n)/Q) ≅ (Z/nZ)^×

            Globale Artin-Abbildung für Q(ζ_n)/Q:
                Art: A_Q^× / Q^× → (Z/nZ)^×
                (a_p)_p ↦ Π_p Frob_p^{v_p(a)}
                Frob_p(ζ_n) = ζ_n^p

            Kronecker-Weber: Jede abelsche Erweiterung K/Q liegt in Q(ζ_n) für n=f(K).

            Dirichlet-Charaktere = Charaktere von (Z/nZ)^× = Charaktere von Gal(Q(ζ_n)/Q)
            → Verbindung Langlands GL_1 ↔ Dirichlet-Charaktere

        @param n Ordnung der zyklotomischen Erweiterung
        @return Dict mit Artin-Abbildungs-Daten
        """
        from sympy import totient, mod_inverse

        # Galois-Gruppe (Z/nZ)^×
        galois_group = [k for k in range(1, n) if math.gcd(k, n) == 1]

        # Frobenius-Elemente: Frob_p = p mod n (für p ∤ n)
        small_primes = list(primerange(2, min(50, n * 3 + 2)))
        frobenius_elements = {}
        for p in small_primes:
            if n % p != 0:  # p unramifiziert
                frob = p % n
                frobenius_elements[p] = {
                    "Frob_p": frob,
                    "action": f"ζ_n ↦ ζ_n^{frob}",
                    "order": self._order_mod_n(frob, n),
                }

        # Artin-Reziprozitäts-Verbindung zu Dirichlet-L-Funktionen
        phi_n = int(totient(n))

        return {
            "extension": f"Q(ζ_{n})/Q",
            "galois_group": f"(Z/{n}Z)^× (Ordnung φ({n})={phi_n})",
            "galois_elements": galois_group,
            "frobenius_at_primes": frobenius_elements,
            "kronecker_weber": f"Jede abelsche K/Q ⊆ Q(ζ_m) für geeignetes m",
            "dirichlet_connection": (
                f"Dirichlet-Char. mod {n} ↔ Char. von Gal(Q(ζ_{n})/Q)"
            ),
            "artin_map": "Art: A_Q^× → (Z/nZ)^× (globale Klassenkörpertheorie)",
        }

    def _order_mod_n(self, a: int, n: int) -> int:
        """
        @brief Berechnet die multiplikative Ordnung von a modulo n.
        @param a Ganzzahl mit gcd(a,n)=1
        @param n Modulus
        @return Kleinste k>0 mit a^k ≡ 1 (mod n)
        """
        if math.gcd(a, n) != 1:
            return 0
        k = 1
        power = a % n
        while power != 1:
            power = (power * a) % n
            k += 1
            if k > n:
                break
        return k


# ---------------------------------------------------------------------------
# 7. Arthur-Selberg-Spurformel (numerische Demo)
# ---------------------------------------------------------------------------

class TraceFormula:
    """
    @brief Arthur-Selberg-Spurformel (numerische Demonstration).
    @description
        Die Selberg-Spurformel für Γ\\H (H = obere Halbebene):
            Σ_π m(π)·tr(π(f)) = Σ_{[γ]} O_γ(f)

        Linke Seite (spektral):
            Σ über automorphe Darstellungen π mit Multiplizität m(π)
        Rechte Seite (geometrisch):
            Σ über Konjugationsklassen [γ] in Γ, orbitale Integrale O_γ(f)

        Einfachste Form (Selberg, 1956): Für Γ = SL_2(Z), f = Wärmekern:
            Σ_n e^{-t·λ_n} = (vol/4π) · t^{-1} + ...  +
            Σ_{[γ] hyp.} (log N(P_γ)/(N(P_γ)^{1/2}-N(P_γ)^{-1/2})) · e^{-l(γ)²/4t} + ...

        Anwendungen:
        - Weyl-Gesetz: #{λ_n ≤ X} ~ vol(Γ\\H)·X/(4π)
        - Selberg-Zeta-Funktion: Z_Γ(s) = Π_{[P] prim hyp.} Π_{k=0}^∞ (1-N(P)^{-(s+k)})

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def geometric_side_demo(
        self,
        gamma_conj_classes: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        @brief Berechnet die geometrische Seite der Spurformel (Demo).
        @description
            Geometrische Seite = Σ_{[γ]} Vol(Γ_γ\\G_γ) · O_γ(f)

            Beitragstypen:
            1. Identität [1]: Vol(Γ\\G)·f(1) (Volume-Term)
            2. Hyperbolische [γ_hyp]: Geodätische im Quotienten
               Beitrag: log(N(P)) / (N(P)^{1/2} - N(P)^{-1/2})
               N(P) = Norm des primitiven Hyperbols P
            3. Elliptische [γ_ell]: Fixpunkte im Quotienten
               Beitrag: π / Vol(Γ_γ) · (Integral über f)
            4. Parabolische [I, γ_par]: Cusps
               Beitrag: -1/4π · ∫ φ'(1/2+it) / φ(1/2+it) · ĥ(t) dt

        @param gamma_conj_classes Liste von Konjugationsklassen-Dicts
        @return Dict mit geometrischen Beiträgen
        """
        geometric_total = 0.0 + 0j
        contributions = []

        for gamma in gamma_conj_classes:
            gtype = gamma.get("type", "unknown")
            N_P = gamma.get("norm", 1.0)    # Norm N(P)
            volume = gamma.get("volume", 1.0)

            if gtype == "identity":
                # Identitäts-Beitrag: Vol(Γ\\H) / 4π (für f = Wärmekern)
                vol_gamma_H = gamma.get("vol_quotient", math.pi / 3)  # SL_2(Z): π/3
                contrib = vol_gamma_H / (4 * math.pi)
                contributions.append({"class": "identity", "contribution": contrib})
                geometric_total += contrib

            elif gtype == "hyperbolic":
                # Hyperbolischer Beitrag: log(N(P)) / (√N(P) - 1/√N(P))
                log_NP = math.log(N_P)
                sqrt_NP = math.sqrt(N_P)
                contrib = log_NP / (sqrt_NP - 1.0 / sqrt_NP)
                contributions.append({
                    "class": "hyperbolic",
                    "N_P": N_P,
                    "contribution": contrib,
                })
                geometric_total += contrib

            elif gtype == "elliptic":
                # Elliptischer Beitrag: ~ π / Vol(Γ_γ)
                order = gamma.get("order", 2)
                contrib = math.pi / (order * volume) if volume > 0 else 0
                contributions.append({
                    "class": "elliptic",
                    "order": order,
                    "contribution": contrib,
                })
                geometric_total += contrib

        return {
            "geometric_side_total": geometric_total,
            "contributions": contributions,
            "num_classes": len(gamma_conj_classes),
        }

    def spectral_side_demo(
        self,
        eigenvalues: List[float]
    ) -> Dict[str, Any]:
        """
        @brief Berechnet die spektrale Seite der Spurformel (Demo).
        @description
            Spektrale Seite = Σ_π m(π) · λ_π (für Test-Funktion f = Wärmekern e^{-t·λ})

            Für Γ = SL_2(Z) auf H:
            - Laplace-Eigenwerte: λ_n = s_n(1-s_n), s_n = 1/2 + it_n
            - Weyl-Gesetz: N(T) = #{λ_n ≤ T} ~ Area(Γ\\H)·T / (4π) = T/12 + O(T^{1/2}·log T)
            - Eisenstein-Spektrum: kontinuierlich, Beitrag ~ -1/(4π)·∫ M'(1/2+it)/M(1/2+it) dt
            - Diskrete Eigenwerte: λ_1 ≈ 91.14 (erste Maass-Form), λ_2 ≈ 190.13

        @param eigenvalues Liste der Laplace-Eigenwerte λ_n
        @return Dict mit spektralen Daten
        """
        if not eigenvalues:
            return {"spectral_side": 0.0, "eigenvalues": []}

        # t_n aus λ_n = 1/4 + t_n²
        t_values = []
        for lam in eigenvalues:
            if lam >= 0.25:
                t_n = math.sqrt(lam - 0.25)
            else:
                # Ausnahme-Eigenwert (kleine Eigenform)
                t_n = complex(0, math.sqrt(0.25 - lam))
            t_values.append(t_n)

        # Weyl-Gesetz Schätzung für T = max(eigenvalues)
        T = max(eigenvalues) if eigenvalues else 0
        area_sl2z = math.pi / 3  # Fläche des Fundamentalbereichs von SL_2(Z)
        weyl_estimate = area_sl2z * T / (4 * math.pi)

        # Selberg-Spurformel Test: Σ f(t_n) für f(t) = e^{-t²·tau}
        tau = 0.01  # Demo-Parameter
        spectral_sum = sum(
            math.exp(-abs(t) ** 2 * tau) if not isinstance(t, complex)
            else abs(cmath.exp(-t ** 2 * tau))
            for t in t_values
        )

        return {
            "spectral_side_approx": spectral_sum,
            "eigenvalues": eigenvalues,
            "t_values": [str(t) for t in t_values],
            "weyl_estimate_N(T)": weyl_estimate,
            "counting": len(eigenvalues),
            "first_maass_eigenvalue": "λ_1 ≈ 91.14 (bekannt für SL_2(Z))",
        }

    def trace_formula_identity_check(
        self,
        N: int,
        weight: int
    ) -> Dict[str, Any]:
        """
        @brief Überprüft die Spurformel-Identität via Weyl-Gesetz für Γ_0(N).
        @description
            Für Γ_0(N) ⊂ SL_2(Z):
            - Index [SL_2(Z):Γ_0(N)] = N · Π_{p|N} (1 + 1/p)
            - Geschlecht g(X_0(N)) via Euler-Charakteristik
            - Dimension S_k(Γ_0(N)) via Riemann-Roch:
              dim S_k ≈ (k-1)/12 · [SL_2(Z):Γ_0(N)] + O(1)  für k ≥ 2 gerade

            Weyl-Gesetz für Maass-Formen auf Γ_0(N)\\H:
                N(T) ~ Area(Γ_0(N)\\H) · T / (4π)
                Area = π/3 · [SL_2(Z):Γ_0(N)]

        @param N Führer (Level) der Gruppe Γ_0(N)
        @param weight Gewicht der Modulformen (k ≥ 2 gerade)
        @return Dict mit Dimensions- und Spurformel-Daten
        """
        # Index [SL_2(Z):Γ_0(N)]
        index = N
        for p, e in factorint(N).items():
            index = index * (1 + 1 / p)
        index = int(index)

        # Fläche des Fundamentalbereichs
        area_sl2z = math.pi / 3
        area_gamma0N = area_sl2z * index

        # Dimension S_k(Γ_0(N)) (Näherung Riemann-Roch für gerades k ≥ 2)
        if weight >= 2 and weight % 2 == 0:
            dim_approx = max(0, int((weight - 1) * index / 12 - index / 6 + 1))
        else:
            dim_approx = 0

        # Genauere Dimension via bekannte Formel für k=2:
        # dim S_2(Γ_0(N)) = g(X_0(N)) (Geschlecht)
        dim_s2_formula = f"dim S_2(Γ_0({N})) = Geschlecht von X_0({N})"

        return {
            "N": N,
            "weight": weight,
            "index_in_SL2Z": index,
            "area_fundamental_domain": area_gamma0N,
            "dim_S_k_approx": dim_approx,
            "dim_S2_note": dim_s2_formula,
            "weyl_law_area": area_gamma0N,
            "weyl_estimate": f"N(T) ~ {area_gamma0N / (4 * math.pi):.4f} · T",
            "trace_formula": "Σ_π m(π)·Spur(π(f)) = Σ_[γ] O_γ(f) (geometrische = spektrale Seite)",
        }


# ---------------------------------------------------------------------------
# 8. Langlands-Dualgruppe
# ---------------------------------------------------------------------------

class LDualGroup:
    """
    @brief Langlands-Dualgruppe ^G und Funktorialitäts-Checks.
    @description
        Die Langlands-Dualgruppe ^G ist die split reductive Gruppe über ℂ,
        deren Wurzeldaten das Dual der Wurzeldaten von G sind.

        Konstruktion:
        1. Wurzeldaten von G: (X*, Φ, X^*, Φ^∨, Δ, Δ^∨)
        2. Duale Wurzeldaten: vertausche Koroots und Roots
        3. ^G = split red. Gruppe mit diesen dualen Wurzeldaten

        Bedeutung für Langlands:
        - Lokale Langlands-Parameter: W_F → ^G(ℂ) (Weil-Gruppen-Homomorphismus)
        - Funktorialität: r: ^H → ^G → Lift F_r: Π(H) → Π(G)
        - L-Gruppe: ^LG = ^G ⋊ Gal(Q̄/Q)

    @author Michael Fuhrmann
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    # Vollständige Dualgruppen-Tabelle
    DUAL_GROUPS: Dict[str, str] = {
        # GL_n: selbst-dual
        "GL_1": "GL_1(C)",
        "GL_2": "GL_2(C)",
        "GL_3": "GL_3(C)",
        "GL_n": "GL_n(C)",
        # SL_n ↔ PGL_n
        "SL_2": "PGL_2(C)",
        "SL_n": "PGL_n(C)",
        "PGL_2": "SL_2(C)",
        "PGL_n": "SL_n(C)",
        # Sp_{2n} ↔ SO_{2n+1} (Typ C_n ↔ B_n)
        "Sp_4": "SO_5(C)",
        "Sp_6": "SO_7(C)",
        "Sp_2n": "SO_{2n+1}(C)",
        "Sp_{2n}": "SO_{2n+1}(C)",  # Variante mit geschweiften Klammern
        # SO ↔ Sp
        "SO_3": "PGL_2(C)",
        "SO_4": "SO_4(C)/Z_2",
        "SO_5": "Sp_4(C)",
        "SO_2n+1": "Sp_{2n}(C)",
        "SO_{2n+1}": "Sp_{2n}(C)",  # Variante mit geschweiften Klammern
        "SO_2n": "SO_{2n}(C)",
        "SO_{2n}": "SO_{2n}(C)",
        # Ausnahme-Gruppen (selbst-dual)
        "G_2": "G_2(C)",
        "F_4": "F_4(C)",
        "E_6": "E_6sc(C)",
        "E_7": "E_7(C)",
        "E_8": "E_8(C)",
        # Unitäre Gruppen
        "U_n": "GL_n(C) (mit Galois-Twisting)",
    }

    def l_group(self, group_type: str) -> Dict[str, Any]:
        """
        @brief Berechnet die L-Gruppe ^L G = ^G ⋊ Gal(Q̄/Q).
        @description
            Die L-Gruppe kodiert die Galois-Wirkung auf ^G:
            - Wenn G über Q split: ^L G = ^G × Gal(Q̄/Q) (direkt)
            - Wenn G über Q nicht split (z.B. U_n):
              ^L G = ^G ⋊ Gal(F/Q), F = Zerfällungskörper

            Beispiele:
            GL_n: ^L GL_n = GL_n(ℂ) × Gal(Q̄/Q)  (split, triviale Wirkung)
            U_n:  ^L U_n = GL_n(ℂ) ⋊ Gal(Q(i)/Q) (komplex-unipotent)
            SO_2n: ^L SO_2n = SO_2n(ℂ) × Gal(Q̄/Q) (outer twisting möglich)

        @param group_type Gruppentyp-String
        @return Dict mit L-Gruppen-Daten
        """
        dual = self.DUAL_GROUPS.get(group_type, f"^{group_type}(C)")

        # Galois-Wirkung: split oder nicht?
        non_split_groups = {"U_n", "U_2", "U_3", "U(n)", "SU_n"}
        is_split = group_type not in non_split_groups

        galois_action = "trivial" if is_split else "nicht-trivial (ϕ-Twisting)"
        product_type = "×" if is_split else "⋊"

        return {
            "group": group_type,
            "connected_dual": dual,
            "L_group": f"{dual} {product_type} Gal(Q̄/Q)",
            "galois_action_on_dual": galois_action,
            "is_split_over_Q": is_split,
            "langlands_parameters": f"W_F → {dual} {product_type} Gal(Q̄/Q)",
        }

    def dual_group(self, group_type: str) -> str:
        """
        @brief Gibt die Langlands-Dualgruppe ^G für G zurück.
        @param group_type Gruppentyp-String
        @return String-Darstellung der Dualgruppe
        """
        return self.DUAL_GROUPS.get(group_type, f"^{group_type}(C)")

    def check_functoriality(
        self,
        source_group: str,
        target_group: str,
        lift_type: str
    ) -> Dict[str, Any]:
        """
        @brief Überprüft die Verfügbarkeit einer Funktorialitäts-Abbildung.
        @description
            Langlands-Funktorialitätsvermutung:
            Gegeben r: ^H → ^G, ∃ F_r: Π(H) → Π(G).

            Bekannte Fälle:
            - sym²: GL_2 → GL_3 (Gelbart-Jacquet 1978)
            - sym³: GL_2 → GL_4 (Kim-Shahidi 2002)
            - sym⁴: GL_2 → GL_5 (Kim 2003)
            - Tensor: GL_2 × GL_2 → GL_4 (Ramakrishnan 2000)
            - BC:    GL_n/F → GL_n/E für E/F abelsch (Arthur-Clozel 1989)
            - AI:    GL_1/K → GL_n/Q für [K:Q]=n (Langlands-Tunnell)
            - ηΛ:   G_2 → GL_7 (Kim 2000, durch Theta-Lifting)
            - Sp_4→GL_4: via Theta-Lifting (Asgari-Shahidi 2006)

        @param source_group Quell-Gruppe H
        @param target_group Ziel-Gruppe G
        @param lift_type Art des Lifts (z.B. "sym2", "base_change", "tensor")
        @return Dict mit Status und Referenz
        """
        known_lifts = {
            ("GL_2", "GL_3", "sym2"): {
                "status": "bewiesen",
                "method": "Rankin-Selberg + Langlands-Shahidi",
                "reference": "Gelbart-Jacquet 1978",
                "r_map": "sym²: GL_2(C) → GL_3(C)",
            },
            ("GL_2", "GL_4", "sym3"): {
                "status": "bewiesen",
                "method": "Langlands-Shahidi Methode",
                "reference": "Kim-Shahidi 2002",
                "r_map": "sym³: GL_2(C) → GL_4(C)",
            },
            ("GL_2", "GL_5", "sym4"): {
                "status": "bewiesen",
                "method": "Converse Theorem + Langlands-Shahidi",
                "reference": "Kim 2003",
                "r_map": "sym⁴: GL_2(C) → GL_5(C)",
            },
            ("GL_2", "GL_6", "sym5"): {
                "status": "offen",
                "method": "Unbekannt",
                "reference": None,
                "r_map": "sym⁵: GL_2(C) → GL_6(C)",
            },
            ("GL_2xGL_2", "GL_4", "tensor"): {
                "status": "bewiesen",
                "method": "Rankin-Selberg Faltung",
                "reference": "Ramakrishnan 2000",
                "r_map": "std⊗std: GL_2(C)×GL_2(C) → GL_4(C)",
            },
            ("GL_2", "GL_2", "base_change"): {
                "status": "bewiesen (abelsch)",
                "method": "Arthur-Clozel Basiswechsel",
                "reference": "Arthur-Clozel 1989",
                "r_map": "BC: ^GL_2 → ^GL_2 (Galois-Aktion)",
            },
        }

        key = (source_group, target_group, lift_type)
        result = known_lifts.get(key, {
            "status": "unbekannt / offen",
            "method": "Kein bekannter Beweis",
            "reference": None,
            "r_map": f"{lift_type}: ^{source_group} → ^{target_group}",
        })

        return {
            "source": source_group,
            "target": target_group,
            "lift_type": lift_type,
            "dual_map": f"r: ^{source_group} → ^{target_group}",
            **result,
        }
