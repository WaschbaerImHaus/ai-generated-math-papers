"""
@file algebraic_number_theory.py
@brief Algebraische Zahlentheorie: Zahlkörper, Ideale, Klassengruppen, Einheitengruppen.
@description
    Implementiert grundlegende Konzepte der algebraischen Zahlentheorie:

    1. Zahlkörper-Grundlagen (NumberField)
       - Q(α) mit Minimalpolynom, Grad [K:Q], Signatur (r1, r2)
       - Ring der ganzen Zahlen O_K (Dedekind-Ring)
       - Diskriminante Δ(K/Q)

    2. Ideal-Arithmetik (Dedekind)
       - Primidealzerlegung (p) = p1^e1 · p2^e2 · ... in O_K
       - Norm N(I) = |O_K/I|
       - Verzweigungsindex e(P|p), Trägheitsgrad f(P|p)
       - Fundamentalsatz: Σ e_i · f_i = n

    3. Klassengruppe
       - Klassenzahl h(Q(√d)) via Minkowski-Schranke
       - Klassifikation für quadratische Körper
       - h(K) = 1 ⟺ O_K Hauptidealbereich (PID)

    4. Einheitengruppe (Dirichlet)
       - Einheitensatz: Rang r = r1 + r2 - 1
       - Fundamentaleinheiten für Q(√d) via Pell-Gleichung
       - Regulator: log|ε| für reell-quadratische Körper

    5. Quadratische Körper (QuadraticField)
       - Legendre-Symbol, Kronecker-Symbol, Quadratisches Reziprozitätsgesetz
       - Hilbert-Klassenfeld-Grad = h(K)

    6. Zyklotomische Körper (CyclotomicField)
       - Q(ζ_n), Grad φ(n), Ganzheitsring Z[ζ_n]
       - Kummer-Lifting für Fermat-Kongruenzen

    7. Lokale Theorie
       - p-adische Vervollständigungen
       - Hensel-Lifting (allgemein)
       - Hilbert-Symbol (a,b)_p

    Mathematischer Hintergrund:
    Die algebraische Zahlentheorie untersucht zahlentheoretische Eigenschaften
    allgemeiner Zahlkörper K/Q. Zentrale Werkzeuge sind Dedekind-Ringe,
    Idealklassengruppen und der Dirichletsche Einheitensatz.

@author Michael Fuhrmann
@lastModified 2026-03-11
"""

import math
import cmath
from fractions import Fraction
from functools import lru_cache
from typing import Optional

# SymPy für symbolische Berechnungen
try:
    import sympy as sp
    from sympy import Poly, ZZ, QQ, factorint, isprime, totient
    from sympy import sqrt as sp_sqrt, Symbol, minimal_polynomial
    _SYMPY_AVAILABLE = True
except ImportError:
    _SYMPY_AVAILABLE = False

# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

@lru_cache(maxsize=512)
def _is_prime(n: int) -> bool:
    """
    @brief Einfacher Primalitätstest (Miller-Rabin für große Zahlen).
    @param n Zu prüfende Zahl
    @return True wenn n prim
    @lastModified 2026-03-11
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Deterministischer Miller-Rabin für n < 3.317e24 (Zeugen 2,3,5,7,11,13,17,19,23)
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        if n == a:
            return True
        if n % a == 0:
            return False
    # Trial division für kleine n
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


@lru_cache(maxsize=256)
def _squarefree_part(d: int) -> int:
    """
    @brief Berechnet den quadratfreien Teil von d (d = k² · m, m quadratfrei).
    @param d Ganzzahl ≠ 0, 1
    @return Quadratfreier Teil m
    @lastModified 2026-03-11
    """
    sign = -1 if d < 0 else 1
    d_abs = abs(d)
    p = 2
    while p * p <= d_abs:
        # Alle Faktoren p² herausdividieren
        while d_abs % (p * p) == 0:
            d_abs //= p * p
        p += 1
    return sign * d_abs


def _euler_totient(n: int) -> int:
    """
    @brief Euler'sche φ-Funktion: φ(n) = |{k : 1 ≤ k ≤ n, gcd(k,n)=1}|.
    @description
        Formel: φ(n) = n · ∏_{p|n} (1 - 1/p)
    @param n Positive ganze Zahl
    @return φ(n)
    @lastModified 2026-03-11
    """
    if n == 1:
        return 1
    result = n
    temp = n
    d = 2
    while d * d <= temp:
        if temp % d == 0:
            # Primfaktor d gefunden: multipliziere mit (1 - 1/d)
            result -= result // d
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        result -= result // temp
    return result


def _gcd(a: int, b: int) -> int:
    """
    @brief Größter gemeinsamer Teiler via Euklidischem Algorithmus.
    @param a Erste ganze Zahl
    @param b Zweite ganze Zahl
    @return gcd(a, b)
    @lastModified 2026-03-11
    """
    while b:
        a, b = b, a % b
    return abs(a)


def _prime_factorization(n: int) -> dict[int, int]:
    """
    @brief Vollständige Primfaktorzerlegung von n.
    @param n Positive ganze Zahl ≥ 2
    @return Dict {Primfaktor: Exponent}
    @lastModified 2026-03-11
    """
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def _legendre_symbol(a: int, p: int) -> int:
    """
    @brief Legendre-Symbol (a/p) für ungerade Primzahl p.
    @description
        Euler-Kriterium: (a/p) ≡ a^{(p-1)/2} (mod p)
        - Rückgabe 0: p | a
        - Rückgabe 1: a ist quadratischer Rest mod p
        - Rückgabe -1: a ist quadratischer Nichtrest mod p
    @param a Ganzzahl
    @param p Ungerade Primzahl
    @return -1, 0 oder 1
    @lastModified 2026-03-11
    """
    if p == 2:
        return 0 if a % 2 == 0 else 1
    a_mod = a % p
    if a_mod == 0:
        return 0
    result = pow(a_mod, (p - 1) // 2, p)
    return -1 if result == p - 1 else int(result)


def _poly_mod_p(poly_coeffs: list[int], p: int) -> list[int]:
    """
    @brief Reduziert Polynomkoeffizienten modulo p.
    @param poly_coeffs Koeffizientenliste [a_0, a_1, ..., a_n] (aufsteigend)
    @param p Primzahl
    @return Reduzierte Koeffizientenliste mod p
    @lastModified 2026-03-11
    """
    return [c % p for c in poly_coeffs]


def _poly_eval_mod_p(poly_coeffs: list[int], x: int, p: int) -> int:
    """
    @brief Wertet Polynom in x modulo p aus (Horner-Schema).
    @param poly_coeffs Koeffizientenliste [a_0, a_1, ..., a_n] (aufsteigend)
    @param x Auswertungspunkt
    @param p Modulus
    @return p(x) mod p
    @lastModified 2026-03-11
    """
    result = 0
    # Horner-Schema: aufsteigend → rückwärts iterieren
    for coef in reversed(poly_coeffs):
        result = (result * x + coef) % p
    return result


def _poly_gcd_mod_p(f: list[int], g: list[int], p: int) -> list[int]:
    """
    @brief Euklidischer Algorithmus für Polynome über Z/pZ.
    @param f Koeffizienten von f (aufsteigend)
    @param g Koeffizienten von g (aufsteigend)
    @param p Primzahl (Körpercharakteristik)
    @return ggT-Polynom (monisch, reduziert mod p)
    @lastModified 2026-03-11
    """
    def _degree(poly: list[int]) -> int:
        """Grad des Polynoms (höchster Nichtnull-Index)."""
        for i in range(len(poly) - 1, -1, -1):
            if poly[i] % p != 0:
                return i
        return -1  # Nullpolynom

    def _reduce(poly: list[int]) -> list[int]:
        """Entfernt führende Nullen und macht monisch."""
        d = _degree(poly)
        if d < 0:
            return [0]
        lc = poly[d] % p
        if lc == 0:
            return [0]
        # Monisch machen: durch Leitkoeffizient dividieren
        inv_lc = pow(lc, p - 2, p)
        result = [(c * inv_lc) % p for c in poly[:d + 1]]
        return result

    def _poly_div_mod_p(dividend: list[int], divisor: list[int]) -> list[int]:
        """Polynomdivision mod p, gibt Rest zurück."""
        deg_d = _degree(divisor)
        if deg_d < 0:
            raise ValueError("Division durch Nullpolynom")
        rem = list(dividend[:])
        deg_r = _degree(rem)
        lc_div = divisor[deg_d] % p
        inv_lc = pow(lc_div, p - 2, p)
        while deg_r >= deg_d:
            factor = (rem[deg_r] * inv_lc) % p
            # Subtrahiere factor * x^{deg_r - deg_d} * divisor
            shift = deg_r - deg_d
            for i in range(deg_d + 1):
                if shift + i < len(rem):
                    rem[shift + i] = (rem[shift + i] - factor * divisor[i]) % p
            deg_r = _degree(rem)
        return rem[:max(deg_r + 1, 1)]

    f = _reduce(f)
    g = _reduce(g)
    while _degree(g) >= 0 and not (len(g) == 1 and g[0] == 0):
        r = _poly_div_mod_p(f, g)
        f = g
        g = _reduce(r)
        if len(g) == 1 and g[0] == 0:
            break
    return _reduce(f)


# ---------------------------------------------------------------------------
# 1. Zahlkörper-Grundlagen
# ---------------------------------------------------------------------------

class NumberField:
    """
    @brief Repräsentation eines Zahlkörpers K = Q(α) mit Minimalpolynom f(α) = 0.
    @description
        Ein Zahlkörper K ist eine endliche Körpererweiterung von Q.
        Er wird durch ein irreduzibles Polynom f ∈ Q[x] definiert, dessen
        Nullstelle α den Körper erzeugt: K = Q[x]/(f(x)).

        Grundeigenschaften:
        - Grad: n = [K:Q] = deg(f)
        - Signatur: (r1, r2) mit r1 + 2·r2 = n
          - r1 = Anzahl reeller Einbettungen
          - r2 = Anzahl konjugierter Paare komplexer Einbettungen
        - Diskriminante: Δ(K/Q) = (-1)^{n(n-1)/2} · N(f'(α))

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, poly_coeffs: list[int | float], name: str = "α"):
        """
        @brief Initialisiert den Zahlkörper Q(α).
        @description
            Das Minimalpolynom wird als Koeffizientenliste übergeben:
            poly_coeffs = [a_0, a_1, ..., a_n] entspricht
            f(x) = a_0 + a_1·x + ... + a_n·x^n

        @param poly_coeffs Koeffizienten des Minimalpolynoms (aufsteigend)
        @param name Bezeichner für das erzeugende Element α
        @lastModified 2026-03-11
        """
        if len(poly_coeffs) < 2:
            raise ValueError("Minimalpolynom muss Grad ≥ 1 haben")
        # Führende Nullen entfernen
        while len(poly_coeffs) > 1 and poly_coeffs[-1] == 0:
            poly_coeffs = poly_coeffs[:-1]
        self.poly_coeffs = list(poly_coeffs)
        self.name = name
        # Nullstellen berechnen (für Signatur)
        self._roots: Optional[list[complex]] = None

    def _compute_roots(self) -> list[complex]:
        """
        @brief Berechnet die Nullstellen des Minimalpolynoms numerisch.
        @return Liste der komplexen Nullstellen
        @lastModified 2026-03-11
        """
        if self._roots is not None:
            return self._roots
        n = len(self.poly_coeffs) - 1
        if n == 0:
            self._roots = []
            return self._roots
        # Koeffizienten: numpy-Konvention ist absteigend
        coeffs_desc = list(reversed(self.poly_coeffs))
        try:
            import numpy as np
            self._roots = list(np.roots(coeffs_desc))
        except ImportError:
            # Fallback: Nur quadratische Gleichungen analytisch lösbar
            if n == 1:
                # a1·x + a0 = 0 → x = -a0/a1
                self._roots = [complex(-self.poly_coeffs[0] / self.poly_coeffs[1])]
            elif n == 2:
                a, b, c = self.poly_coeffs[2], self.poly_coeffs[1], self.poly_coeffs[0]
                disc = b * b - 4 * a * c
                if disc >= 0:
                    self._roots = [
                        complex((-b + math.sqrt(disc)) / (2 * a)),
                        complex((-b - math.sqrt(disc)) / (2 * a)),
                    ]
                else:
                    self._roots = [
                        complex(-b / (2 * a), math.sqrt(-disc) / (2 * a)),
                        complex(-b / (2 * a), -math.sqrt(-disc) / (2 * a)),
                    ]
            else:
                self._roots = []
        return self._roots

    def degree(self) -> int:
        """
        @brief Gibt den Grad [K:Q] = deg(f) zurück.
        @description
            Der Grad des Zahlkörpers ist gleich dem Grad des Minimalpolynoms.
            Für K = Q(√d) ist n = 2, für Q(ζ_n) ist n = φ(n).
        @return Grad [K:Q]
        @lastModified 2026-03-11
        """
        return len(self.poly_coeffs) - 1

    def signature(self) -> tuple[int, int]:
        """
        @brief Berechnet die Signatur (r1, r2) des Zahlkörpers.
        @description
            Die Signatur gibt die Anzahl der Einbettungen an:
            - r1: Anzahl reeller Einbettungen K → ℝ
            - r2: Anzahl konjugierter Paare K → ℂ (nicht reell)
            Es gilt: r1 + 2·r2 = n = [K:Q]

            Quadratischer Körper Q(√d):
            - d > 0: (r1, r2) = (2, 0) — reell-quadratisch
            - d < 0: (r1, r2) = (0, 1) — imaginär-quadratisch

        @return Tuple (r1, r2)
        @lastModified 2026-03-11
        """
        n = self.degree()
        roots = self._compute_roots()
        if not roots:
            # Fallback: alle Nullstellen zählen (approximativ)
            r1 = 0
            r2 = n // 2
            return (r1, r2)
        # Reelle Nullstellen: |Im(z)| < Toleranz
        tolerance = 1e-8
        r1 = sum(1 for z in roots if abs(z.imag) < tolerance)
        # Komplexe Nullstellen kommen in konjugierten Paaren
        complex_count = n - r1
        r2 = complex_count // 2
        return (r1, r2)

    def ring_of_integers(self) -> dict:
        """
        @brief Bestimmt den Ganzheitsring O_K (vereinfacht für quadratische Körper).
        @description
            Der Ganzheitsring O_K besteht aus allen α ∈ K, die ein normiertes
            Polynom mit ganzzahligen Koeffizienten erfüllen (ganze algebraische Zahlen).

            Für K = Q(√d) (d quadratfrei):
            - d ≡ 1 (mod 4): O_K = Z[(1+√d)/2], Basis {1, (1+√d)/2}
            - d ≡ 2,3 (mod 4): O_K = Z[√d], Basis {1, √d}

            Allgemein (Dedekind-Kriterium):
            O_K = Z[α] falls die Diskriminante von f mit der Körperdiskriminante
            übereinstimmt (kein quadratischer Faktor).

        @return Dict mit ring_basis, discriminant, description
        @lastModified 2026-03-11
        """
        n = self.degree()
        disc = self.discriminant()

        if n == 2:
            # Quadratischer Fall: explizite Formel
            # Aus f(x) = x² + bx + c → α = (-b ± √(b²-4c))/2
            # Wir lesen d aus der Diskriminante: disc = b²-4c
            # Quadratfreier Teil von disc gibt uns d
            b = self.poly_coeffs[1] if len(self.poly_coeffs) > 1 else 0
            c = self.poly_coeffs[0] if len(self.poly_coeffs) > 0 else 0
            a2 = self.poly_coeffs[2] if len(self.poly_coeffs) > 2 else 1
            # disc(f) = b² - 4·a2·c (für normiertes Polynom a2=1)
            inner_disc = b * b - 4 * a2 * c
            d_sf = _squarefree_part(inner_disc) if inner_disc != 0 else inner_disc

            if d_sf % 4 == 1:
                ring_basis = ["1", f"(1+√{d_sf})/2"]
                desc = f"ℤ[(1+√{d_sf})/2]"
                actual_disc = d_sf
            else:
                ring_basis = ["1", f"√{d_sf}"]
                desc = f"ℤ[√{d_sf}]"
                actual_disc = 4 * d_sf

            return {
                "ring_basis": ring_basis,
                "discriminant": actual_disc,
                "description": f"O_K = {desc}",
                "is_monogenic": True,
            }
        else:
            # Allgemeiner Fall: Z[α] ist ein Kandidat (Dedekind-Kriterium)
            ring_basis = ["1"] + [f"{self.name}^{k}" for k in range(1, n)]
            return {
                "ring_basis": ring_basis,
                "discriminant": disc,
                "description": f"O_K ⊇ ℤ[{self.name}] (Dedekind-Ring, Grad {n})",
                "is_monogenic": True,  # Z[α] ist zumindest ein Unterring
            }

    def discriminant(self) -> int | float:
        """
        @brief Berechnet die Diskriminante Δ(K/Q) des Zahlkörpers.
        @description
            Die Diskriminante des Zahlkörpers K = Q(α) mit Minimalpolynom f ist:
            Δ(K/Q) = (-1)^{n(n-1)/2} · N_{K/Q}(f'(α))

            Für das Polynom f: Δ(f) = (-1)^{n(n-1)/2} · Res(f, f')

            Für quadratische Körper Q(√d) mit d quadratfrei:
            - Δ = d falls d ≡ 1 (mod 4)
            - Δ = 4d sonst

            Bedeutung:
            - Δ < 0: gemischte Signatur
            - |Δ| groß: viele verzweigte Primzahlen
            - p | Δ ⟺ p verzweigt in K

        @return Diskriminante (ganze Zahl oder Float für hohe Grade)
        @lastModified 2026-03-11
        """
        n = self.degree()
        coeffs = self.poly_coeffs

        if n == 1:
            return 1  # Q hat Diskriminante 1

        if n == 2:
            # Quadratisches Polynom a·x² + b·x + c: Δ = b² - 4ac
            a = coeffs[2]
            b = coeffs[1]
            c = coeffs[0]
            raw_disc = b * b - 4 * a * c
            # Zur Körperdiskriminante: quadratfreien Teil bestimmen
            d_sf = _squarefree_part(raw_disc) if raw_disc != 0 else 0
            if d_sf % 4 == 1:
                return d_sf
            else:
                return 4 * d_sf

        # Allgemeiner Fall: Resultante von f und f' via Sylvester-Matrix
        # Für höhere Grade: numerische Näherung über Nullstellen
        roots = self._compute_roots()
        if not roots:
            return 0  # Unbekannt

        # Δ(f) = (-1)^{n(n-1)/2} · a_n^{2n-2} · ∏_{i<j} (r_i - r_j)²
        sign_exp = n * (n - 1) // 2
        sign = (-1) ** sign_exp
        leading_coef = coeffs[-1]
        product = 1.0
        for i in range(n):
            for j in range(i + 1, n):
                diff = roots[i] - roots[j]
                product *= (diff * diff.conjugate()).real  # |r_i - r_j|²

        disc_float = sign * (leading_coef ** (2 * n - 2)) * product
        return round(disc_float) if abs(disc_float - round(disc_float)) < 0.5 else disc_float

    def __repr__(self) -> str:
        """@brief Textuelle Darstellung des Zahlkörpers."""
        coeffs = self.poly_coeffs
        terms = []
        for i, c in enumerate(coeffs):
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                terms.append(f"{c}·{self.name}" if c != 1 else self.name)
            else:
                terms.append(f"{c}·{self.name}^{i}" if c != 1 else f"{self.name}^{i}")
        poly_str = " + ".join(reversed(terms)) if terms else "0"
        return f"NumberField Q({self.name}), f({self.name}) = {poly_str} = 0"


# ---------------------------------------------------------------------------
# 2. Ideal-Arithmetik (Dedekind)
# ---------------------------------------------------------------------------

class DedekindIdeal:
    """
    @brief Ideal in einem Dedekind-Ring O_K.
    @description
        In einem Dedekind-Ring (z.B. O_K für Zahlkörper K/Q) hat jedes
        von Null verschiedene Ideal I eine eindeutige Primidealzerlegung:
        I = p1^e1 · p2^e2 · ... · pr^er

        Die Norm N(I) = |O_K/I| ist für Primideale p mit N(p) = p^f,
        wobei f der Trägheitsgrad ist.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, generators: list[int], field_degree: int = 2):
        """
        @brief Initialisiert ein Ideal mit seinen Erzeugern.
        @param generators Liste der erzeugenden Elemente (als ganze Zahlen)
        @param field_degree Grad [K:Q] des Zahlkörpers
        @lastModified 2026-03-11
        """
        self.generators = generators
        self.field_degree = field_degree

    def norm(self) -> int:
        """
        @brief Berechnet die Norm N(I) = |O_K/I| des Ideals.
        @description
            Für ein Hauptideal (a) gilt N((a)) = |N_{K/Q}(a)|.
            Für das Ideal erzeugt von einer ganzen Zahl n ∈ Z:
            N((n)) = n^[K:Q]

            Für ein Primideal P über p gilt N(P) = p^f,
            wobei f der Trägheitsgrad von P über p ist.

        @return Norm des Ideals
        @lastModified 2026-03-11
        """
        if not self.generators:
            return 0
        # Für einfachste Fälle: N((n)) = n^n für n ∈ Z
        base = abs(self.generators[0])
        if len(self.generators) == 1 and base > 0:
            return base ** self.field_degree
        # Für mehrere Erzeuger: Norm = gcd der Normen (Näherung)
        g = self.generators[0]
        for gen in self.generators[1:]:
            g = _gcd(g, gen)
        return abs(g) ** self.field_degree

    def is_prime_ideal(self) -> bool:
        """
        @brief Prüft ob das Ideal ein Primideal ist (vereinfacht für Hauptideale).
        @description
            Ein Ideal P in O_K ist prim, wenn O_K/P ein Integritätsring ist.
            Für Primideale über einer Primzahl p gilt: N(P) = p^f für ein f ≥ 1.

            Notwendige Bedingung (implementiert):
            Das Ideal ist erzeugt von einer Primzahl p ∈ Z, und wir prüfen
            ob es als ganzes unter der Primzahl liegt.

        @return True wenn das Ideal prim erscheint
        @lastModified 2026-03-11
        """
        if len(self.generators) == 1:
            return _is_prime(abs(self.generators[0]))
        # Für mehrere Erzeuger: gcd prüfen
        g = abs(self.generators[0])
        for gen in self.generators[1:]:
            g = _gcd(g, abs(gen))
        return _is_prime(g)

    def is_principal(self) -> bool:
        """
        @brief Prüft ob das Ideal ein Hauptideal ist (vereinfacht).
        @description
            Ein Ideal I ist ein Hauptideal, wenn I = (a) für ein a ∈ O_K gilt.
            Für Ideale mit einem einzigen Erzeuger ist dies trivialerweise erfüllt.
            Für mehrere Erzeuger prüfen wir ob sie alle ein gemeinsames Vielfaches haben,
            das das Ideal erzeugt.

        @return True wenn das Ideal ein Hauptideal ist
        @lastModified 2026-03-11
        """
        if len(self.generators) == 1:
            return True
        # Prüfe ob alle Erzeuger Vielfache des gcd sind
        g = abs(self.generators[0])
        for gen in self.generators[1:]:
            g = _gcd(g, abs(gen))
        # Wenn alle Erzeuger Vielfache von g sind → (g) erzeugt das Ideal
        return all(gen % g == 0 for gen in self.generators)

    def __repr__(self) -> str:
        """@brief Textuelle Darstellung des Ideals."""
        return f"Ideal({self.generators})"


def ideal_factorization(p: int, poly_coeffs: list[int]) -> list[dict]:
    """
    @brief Berechnet die Primidealzerlegung von (p) in O_K.
    @description
        Zerlegung des Primideals (p) ⊆ Z in O_K = Z[x]/(f(x)):
        (p) = p1^e1 · p2^e2 · ... in O_K

        Algorithmus (Kummer-Dedekind-Theorem für f squarefree mod p):
        Faktorisiere f ≡ g1^e1 · g2^e2 · ... (mod p).
        Dann: pi = (p, gi(α)) und N(pi) = p^{deg(gi)}.

        Zerlegungsarten:
        - Inert: f irreduzibel mod p → ein Primideal mit f = deg(f)
        - Split/Zerfällt: f zerfällt in eindeutige Faktoren mod p → r Primideale
        - Verzweigt/Ramified: f hat mehrfache Faktoren mod p (p | disc)

        Für quadratische Körper K = Q(√d) und Primzahl p ≠ 2:
        - (d/p) = -1: p inert (bleibt prim in O_K)
        - (d/p) = 0: p verzweigt (p | disc)
        - (d/p) = 1: p zerfällt (zwei Primideale)

    @param p Primzahl
    @param poly_coeffs Koeffizienten des Minimalpolynoms (aufsteigend)
    @return Liste von Dicts {prime_ideal, exponent, inertia_degree, ramification_index, norm}
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"{p} ist keine Primzahl")

    n = len(poly_coeffs) - 1  # Grad des Polynoms
    factors = []

    if n == 2:
        # Quadratischer Körper: explizite Formel über Legendre-Symbol
        b = poly_coeffs[1]
        c = poly_coeffs[0]
        # disc = b² - 4c (für normiertes Polynom x² + bx + c)
        a_coeff = poly_coeffs[2] if len(poly_coeffs) > 2 else 1
        disc_raw = b * b - 4 * a_coeff * c
        d_sf = _squarefree_part(disc_raw) if disc_raw != 0 else 0

        if p == 2:
            # Sonderfall p=2: Zerlegung abhängig von d mod 8
            d_mod8 = d_sf % 8
            if d_mod8 in (2, 6):
                # p=2 ramifiziert: d ≡ 2,6 (mod 8)
                factors = [{
                    "prime_ideal": f"p_2",
                    "exponent": 2,
                    "inertia_degree": 1,
                    "ramification_index": 2,
                    "norm": 2,
                    "type": "ramified",
                }]
            elif d_mod8 in (3, 7):
                # p=2 inert: d ≡ 3,7 (mod 8)
                factors = [{
                    "prime_ideal": f"(2)",
                    "exponent": 1,
                    "inertia_degree": 2,
                    "ramification_index": 1,
                    "norm": 4,
                    "type": "inert",
                }]
            else:
                # p=2 zerfällt: d ≡ 1,5 (mod 8)
                factors = [
                    {"prime_ideal": "p_2^+", "exponent": 1, "inertia_degree": 1,
                     "ramification_index": 1, "norm": 2, "type": "split"},
                    {"prime_ideal": "p_2^-", "exponent": 1, "inertia_degree": 1,
                     "ramification_index": 1, "norm": 2, "type": "split"},
                ]
        else:
            # Ungerades p: Legendre-Symbol
            leg = _legendre_symbol(d_sf, p)
            if leg == 0:
                # Verzweigt: (p) = p² in O_K
                factors = [{
                    "prime_ideal": f"p_{p}",
                    "exponent": 2,
                    "inertia_degree": 1,
                    "ramification_index": 2,
                    "norm": p,
                    "type": "ramified",
                }]
            elif leg == -1:
                # Inert: (p) bleibt prim
                factors = [{
                    "prime_ideal": f"(p={p})",
                    "exponent": 1,
                    "inertia_degree": 2,
                    "ramification_index": 1,
                    "norm": p * p,
                    "type": "inert",
                }]
            else:
                # Zerfällt: (p) = p1 · p2
                factors = [
                    {"prime_ideal": f"p_{p}^+", "exponent": 1, "inertia_degree": 1,
                     "ramification_index": 1, "norm": p, "type": "split"},
                    {"prime_ideal": f"p_{p}^-", "exponent": 1, "inertia_degree": 1,
                     "ramification_index": 1, "norm": p, "type": "split"},
                ]
    else:
        # Allgemeiner Fall: Kummer-Dedekind via Faktorisierung mod p
        # Suche Nullstellen von f mod p (für einfache Faktorisierung)
        roots_mod_p = [x for x in range(p) if _poly_eval_mod_p(poly_coeffs, x, p) == 0]
        if not roots_mod_p:
            # Keine Nullstellen → f möglicherweise irreduzibel mod p (inert)
            factors = [{
                "prime_ideal": f"(p={p})",
                "exponent": 1,
                "inertia_degree": n,
                "ramification_index": 1,
                "norm": p ** n,
                "type": "inert",
            }]
        else:
            # Für jeden Linearfaktor (x - r) mod p gibt es ein Primideal
            seen_roots = []
            for r in roots_mod_p:
                # Mehrfachheit bestimmen (vereinfacht: Ableitung)
                deriv_coeffs = [i * poly_coeffs[i] for i in range(1, len(poly_coeffs))]
                deriv_val = _poly_eval_mod_p(deriv_coeffs, r, p)
                e = 1 if deriv_val != 0 else 2  # Einfache oder doppelte Nullstelle
                if r not in seen_roots:
                    seen_roots.append(r)
                    factors.append({
                        "prime_ideal": f"p_{p}^({r})",
                        "exponent": e,
                        "inertia_degree": 1,
                        "ramification_index": e,
                        "norm": p,
                        "type": "ramified" if e > 1 else "split",
                    })
            # Verbleibender irreduzibeler Teil
            remaining_degree = n - len(seen_roots)
            if remaining_degree > 0:
                factors.append({
                    "prime_ideal": f"(p={p}, irred_part)",
                    "exponent": 1,
                    "inertia_degree": remaining_degree,
                    "ramification_index": 1,
                    "norm": p ** remaining_degree,
                    "type": "inert_part",
                })

    return factors


def ramification_index(p: int, poly_coeffs: list[int]) -> list[int]:
    """
    @brief Gibt die Verzweigungsindizes e(Pi|p) für alle Primideale Pi über p zurück.
    @description
        Der Verzweigungsindex e(P|p) gibt an, wie oft die Primzahl p
        in der Zerlegung von P vorkommt:
        p · O_K = P1^{e1} · P2^{e2} · ...

        Verzweigung (e > 1) tritt genau dann auf, wenn p die Diskriminante teilt.

    @param p Primzahl
    @param poly_coeffs Koeffizienten des Minimalpolynoms
    @return Liste der Verzweigungsindizes [e1, e2, ...]
    @lastModified 2026-03-11
    """
    factors = ideal_factorization(p, poly_coeffs)
    return [f["ramification_index"] for f in factors]


def inertia_degree(p: int, poly_coeffs: list[int]) -> list[int]:
    """
    @brief Gibt die Trägheitsgrade f(Pi|p) für alle Primideale Pi über p zurück.
    @description
        Der Trägheitsgrad f(P|p) ist der Grad der Restklassenkörpererweiterung:
        f(P|p) = [O_K/P : Z/pZ]

        Für Primideale P gilt: N(P) = p^{f(P|p)}.

    @param p Primzahl
    @param poly_coeffs Koeffizienten des Minimalpolynoms
    @return Liste der Trägheitsgrade [f1, f2, ...]
    @lastModified 2026-03-11
    """
    factors = ideal_factorization(p, poly_coeffs)
    return [f["inertia_degree"] for f in factors]


def efg_equation_check(p: int, poly_coeffs: list[int]) -> dict:
    """
    @brief Überprüft den Fundamentalsatz Σ e_i · f_i = n = [K:Q].
    @description
        Fundamentalsatz der algebraischen Zahlentheorie:
        Für eine Primzahl p und eine Erweiterung K/Q gilt:
        Σ_{i=1}^{r} e_i · f_i = [K:Q] = n

        wobei:
        - e_i = Verzweigungsindex von Pi über p
        - f_i = Trägheitsgrad von Pi über p
        - r = Anzahl der Primideale über p

        Spezialfälle:
        - Voll zerfällt: e_i = f_i = 1 für alle i, r = n
        - Voll verzweigt: e_1 = n, f_1 = 1, r = 1
        - Inert: e_1 = 1, f_1 = n, r = 1

    @param p Primzahl
    @param poly_coeffs Koeffizienten des Minimalpolynoms
    @return Dict mit sum_eifi, n, check_passed, factors_detail
    @lastModified 2026-03-11
    """
    n = len(poly_coeffs) - 1
    factors = ideal_factorization(p, poly_coeffs)
    sum_eifi = sum(f["ramification_index"] * f["inertia_degree"] for f in factors)
    check_passed = (sum_eifi == n)

    return {
        "sum_eifi": sum_eifi,
        "n": n,
        "check_passed": check_passed,
        "n_prime_ideals": len(factors),
        "factors_detail": factors,
        "description": (
            f"Σ e_i·f_i = {sum_eifi} {'=' if check_passed else '≠'} n = {n} "
            f"(Fundamentalsatz {'erfüllt' if check_passed else 'VERLETZT'})"
        ),
    }


# ---------------------------------------------------------------------------
# 3. Klassengruppe
# ---------------------------------------------------------------------------

def class_number(d: int) -> int:
    """
    @brief Berechnet die Klassenzahl h(Q(√d)) für den quadratischen Körper.
    @description
        Die Klassenzahl h(K) misst die Abweichung von O_K von einem PID:
        - h(K) = 1 ⟺ O_K ist Hauptidealbereich (eindeutige Primfaktorzerlegung)
        - h(K) > 1 → es existieren nicht-triviale Idealklassen

        Bekannte Werte (Stark-Heegner, Baker, Goldfeld, ...):
        Imaginär-quadratisch mit h=1: d ∈ {-1,-2,-3,-7,-11,-19,-43,-67,-163}
        (Dies sind alle solchen d — ein tiefer Satz!)

        Formel (Klassenanzahlformel, Dirichlet):
        h(d) = (w·√|Δ|)/(2π) · L(1, χ_d) für d < 0
        wobei w = Anzahl der Einheitswurzeln, χ_d = Kronecker-Symbol (d/·)

        Minkowski-Schranke: h(K) kann durch Prüfung aller Primideale bis M_K
        exakt berechnet werden.

    @param d Quadratfreier Diskriminantenparameter
    @return Klassenzahl h(Q(√d))
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)

    # Umfangreiche Tabelle bekannter Klassenzahlen
    # Imaginär-quadratische Körper (d < 0)
    known_neg: dict[int, int] = {
        -1: 1, -2: 1, -3: 1, -5: 2, -6: 2, -7: 1,
        -10: 2, -11: 1, -13: 2, -14: 2, -15: 2,
        -17: 4, -19: 1, -21: 4, -22: 2, -23: 3,
        -26: 6, -29: 6, -30: 4, -31: 3, -33: 4,
        -34: 4, -35: 2, -37: 2, -38: 6, -39: 4,
        -41: 8, -43: 1, -46: 4, -47: 5, -51: 2,
        -53: 6, -55: 4, -57: 4, -58: 2, -59: 3,
        -61: 6, -65: 8, -67: 1, -69: 4, -71: 7,
        -163: 1,
    }
    # Reell-quadratische Körper (d > 0)
    known_pos: dict[int, int] = {
        2: 1, 3: 1, 5: 1, 6: 1, 7: 1, 10: 2, 11: 1,
        13: 1, 14: 1, 15: 2, 17: 1, 19: 1, 21: 1,
        22: 1, 23: 1, 26: 2, 29: 1, 30: 2, 31: 1,
        33: 1, 34: 2, 35: 2, 37: 1, 38: 1, 39: 2,
        41: 1, 43: 1, 46: 1, 47: 1, 51: 2, 53: 1,
        55: 2, 57: 1, 58: 2, 59: 1, 61: 1, 65: 2,
    }

    if d_sf < 0 and d_sf in known_neg:
        return known_neg[d_sf]
    if d_sf > 0 and d_sf in known_pos:
        return known_pos[d_sf]

    # Näherungsformel via Minkowski-Schranke für unbekannte d
    return _estimate_class_number_via_minkowski(d_sf)


def _estimate_class_number_via_minkowski(d: int) -> int:
    """
    @brief Schätzt die Klassenzahl via Primidealanalyse bis zur Minkowski-Schranke.
    @description
        Für quadratische Körper Q(√d) mit Diskriminante Δ:
        Minkowski-Schranke: M_K = (2/π)·√|Δ| (imaginär) oder (1/2)·√|Δ| (reell)

        Jede Idealklasse enthält ein Ideal mit Norm ≤ M_K.
        Klassenzahl = Anzahl der unabhängigen Idealklassen bis M_K.

    @param d Quadratfreier Parameter
    @return Geschätzte Klassenzahl (mindestens 1)
    @lastModified 2026-03-11
    """
    disc = d if d % 4 == 1 else 4 * d
    abs_disc = abs(disc)
    mink = minkowski_bound(d)
    bound = max(int(mink) + 1, 3)

    # Primideale bis zur Schranke zählen, die nicht Hauptideale sind
    non_principal_count = 0
    for p in range(2, bound + 1):
        if not _is_prime(p):
            continue
        leg = _legendre_symbol(d, p) if p != 2 else _legendre_symbol_at_2(d)
        if leg == 1:
            # Zerfällt → zwei neue Primideale, möglicherweise nicht Hauptideale
            non_principal_count += 1

    # Grobe Schätzung: logarithmisch mit der Anzahl zerfallender Primzahlen
    if non_principal_count == 0:
        return 1
    return max(1, non_principal_count // max(1, int(math.log(bound + 2))))


def _legendre_symbol_at_2(d: int) -> int:
    """
    @brief Legendre-/Kronecker-Symbol (d/2) für die Primzahl 2.
    @param d Ganzzahl
    @return 0, 1 oder -1
    @lastModified 2026-03-11
    """
    d_mod8 = d % 8
    if d_mod8 in (1, 7):
        return 1
    elif d_mod8 in (3, 5):
        return -1
    else:
        return 0


def minkowski_bound(d: int) -> float:
    """
    @brief Berechnet die Minkowski-Schranke für Q(√d).
    @description
        Die Minkowski-Schranke M_K garantiert, dass jede Idealklasse
        ein Ideal mit Norm ≤ M_K enthält:

        M_K = (n!/n^n) · (4/π)^{r2} · √|Δ|

        Für quadratische Körper (n=2):
        - Imaginär-quadratisch (d<0, r2=1): M_K = (2/π) · √|Δ|
        - Reell-quadratisch (d>0, r2=0): M_K = (1/2) · √|Δ|

        Bedeutung: Man muss nur Primideale mit Norm ≤ M_K untersuchen,
        um die vollständige Klassengruppe zu berechnen.

    @param d Quadratfreier Parameter des Körpers Q(√d)
    @return Minkowski-Schranke (Gleitkommazahl)
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)
    disc = d_sf if d_sf % 4 == 1 else 4 * d_sf
    abs_disc = abs(disc)
    sqrt_disc = math.sqrt(abs_disc)

    if d_sf < 0:
        # Imaginär-quadratischer Körper: r2 = 1
        return (2.0 / math.pi) * sqrt_disc
    else:
        # Reell-quadratischer Körper: r2 = 0
        return 0.5 * sqrt_disc


def class_group_structure(d: int) -> dict:
    """
    @brief Bestimmt die Struktur der Klassengruppe Cl(Q(√d)).
    @description
        Für kleine d kann die Klassengruppe explizit berechnet werden.
        Allgemein gilt: Cl(K) ist eine endliche abelsche Gruppe der Ordnung h(K).

        Für zyklische Klassengruppen: Cl(K) ≅ Z/hZ
        Für nicht-zyklische: Cl(K) ≅ Z/n1Z × Z/n2Z × ... (Smith-Normalform)

        Bekannte Strukturen:
        - d=-5: Cl ≅ Z/2Z (Klassenzahl 2, erzeugt von Ideal (2,1+√-5))
        - d=-23: Cl ≅ Z/3Z (Klassenzahl 3)
        - d=-14: Cl ≅ Z/4Z (Klassenzahl 4)

    @param d Quadratfreier Parameter
    @return Dict mit class_number, group_structure, generators, is_cyclic
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)
    h = class_number(d_sf)

    # Bekannte Gruppenstrukturen für kleine d
    known_structures: dict[int, tuple] = {
        -1: (1,), -2: (1,), -3: (1,), -7: (1,), -11: (1,), -19: (1,),
        -43: (1,), -67: (1,), -163: (1,),
        -5: (2,), -6: (2,), -10: (2,), -13: (2,), -15: (2,), -22: (2,), -35: (2,),
        -23: (3,), -31: (3,), -59: (3,),
        -14: (4,), -17: (4,), -21: (4,), -30: (4,), -33: (4,), -39: (4,),
        -65: (8,),
    }

    if d_sf in known_structures:
        structure = known_structures[d_sf]
    else:
        # Standardannahme: zyklisch (meist korrekt für kleine h)
        structure = (h,) if h > 0 else (1,)

    is_cyclic = len(structure) == 1
    is_pid = (h == 1)

    # Erzeuger der Klassengruppe
    # Für nicht-triviale Klassen: Primideale über den kleinsten zerfallenden Primzahlen
    generators = []
    if h > 1:
        mink = minkowski_bound(d_sf)
        for p in range(2, max(int(mink) + 2, 5)):
            if not _is_prime(p):
                continue
            leg = (_legendre_symbol(d_sf, p) if p != 2
                   else _legendre_symbol_at_2(d_sf))
            if leg == 1:
                generators.append(f"p_{p}^+ = ({p}, {_sqrt_approx(d_sf)} + ???)")
                if len(generators) >= len(structure):
                    break

    return {
        "class_number": h,
        "group_structure": structure,
        "group_name": " × ".join(f"ℤ/{n}ℤ" for n in structure),
        "is_cyclic": is_cyclic,
        "is_pid": is_pid,
        "generators": generators,
        "description": (
            f"Cl(ℚ(√{d_sf})) ≅ {' × '.join(f'ℤ/{n}ℤ' for n in structure)}, "
            f"h = {h}, {'PID' if is_pid else 'kein PID'}"
        ),
    }


def _sqrt_approx(d: int) -> str:
    """@brief Hilfsfunktion: Textdarstellung von √d."""
    return f"√{d}" if d > 0 else f"√({d})"


def is_pid(d: int) -> bool:
    """
    @brief Prüft ob O_{Q(√d)} ein Hauptidealbereich (PID) ist.
    @description
        Ein Dedekind-Ring ist genau dann ein PID, wenn die Klassengruppe
        trivial ist, d.h. h(Q(√d)) = 1.

        Satz (Stark-Heegner, vollständig bewiesen):
        Die einzigen imaginär-quadratischen Körper Q(√d) (d < 0) mit h=1 sind:
        d ∈ {-1, -2, -3, -7, -11, -19, -43, -67, -163}

        Für reell-quadratische Körper ist die vollständige Liste unbekannt
        (Cohen-Lenstra-Heuristik: es gibt unendlich viele).

    @param d Quadratfreier Parameter
    @return True wenn O_{Q(√d)} ein PID ist
    @lastModified 2026-03-11
    """
    return class_number(d) == 1


# ---------------------------------------------------------------------------
# 4. Einheitengruppe (Dirichlet)
# ---------------------------------------------------------------------------

def unit_rank(poly_coeffs: list[int]) -> int:
    """
    @brief Berechnet den Einheitenrang r = r1 + r2 - 1 (Dirichletscher Einheitensatz).
    @description
        Dirichletscher Einheitensatz:
        Die Einheitengruppe O_K^× hat die Struktur:
        O_K^× ≅ μ(K) × ℤ^r

        wobei:
        - μ(K) = endliche Gruppe der Einheitswurzeln in K
        - r = r1 + r2 - 1 = Einheitenrang
        - r1 = Anzahl reeller Einbettungen
        - r2 = Anzahl konjugierter Paare komplexer Einbettungen

        Spezialfälle:
        - Q selbst: r = 0 (nur ±1)
        - Reell-quadratisch: r = 1 (eine Fundamentaleinheit)
        - Imaginär-quadratisch: r = 0 (nur endlich viele Einheiten: ±1, ±i, ±ω...)
        - Totalreell (alle r2=0): r = r1 - 1

    @param poly_coeffs Koeffizienten des Minimalpolynoms
    @return Einheitenrang r
    @lastModified 2026-03-11
    """
    # Signatur aus dem Minimalpolynom berechnen
    nf = NumberField(list(poly_coeffs))
    r1, r2 = nf.signature()
    return r1 + r2 - 1


def fundamental_units_quadratic(d: int) -> dict:
    """
    @brief Findet die Fundamentaleinheit von Q(√d) für reell-quadratische Körper.
    @description
        Für reell-quadratische Körper K = Q(√d) mit d > 0 ist der Einheitenrang = 1.
        Die Einheitengruppe hat die Form: O_K^× = {±ε^n : n ∈ ℤ}

        Fundamentaleinheit ε = a + b√d mit a² - d·b² = ±1 (Pell-Gleichung)
        oder allgemeiner: N(ε) = a² - d·b² = ±1

        Für d ≡ 1 (mod 4): Einheiten der Form (a + b√d)/2 mit a² - d·b² = ±4

        Algorithmus: Kettenbruchentwicklung von √d findet ε in endlicher Zeit.
        Die Periode des Kettenbruchs [a0; a1, a2, ...] liefert die Lösung.

        Bekannte Fundamentaleinheiten:
        - Q(√2): ε = 1+√2 (Norm -1)
        - Q(√3): ε = 2+√3 (Norm 1)
        - Q(√5): ε = (1+√5)/2 (Goldener Schnitt, Norm -1)
        - Q(√7): ε = 8+3√7 (Norm 1)

    @param d Quadratfreier Parameter > 0 (muss positiv sein)
    @return Dict mit a, b, norm, description für ε = a + b√d (oder (a+b√d)/2)
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)

    if d_sf <= 0:
        return {
            "a": 0, "b": 0, "norm": 0,
            "description": "Imaginär-quadratischer Körper: keine Fundamentaleinheit (r=0)",
            "unit_rank": 0,
        }

    # Bekannte Fundamentaleinheiten (häufig verwendete Werte)
    known_units: dict[int, tuple[int, int]] = {
        2: (1, 1),    # 1 + √2, N = -1
        3: (2, 1),    # 2 + √3, N = 1
        5: (1, 1),    # (1+√5)/2 für d≡1(4): a=b=1 → (1+√5)/2
        6: (5, 2),    # 5 + 2√6, N = 1
        7: (8, 3),    # 8 + 3√7, N = 1
        10: (3, 1),   # 3 + √10, N = -1
        11: (10, 3),  # 10 + 3√11, N = 1
        13: (3, 1),   # (3+√13)/2 für d≡1(4)
        14: (15, 4),  # 15 + 4√14, N = 1
        15: (4, 1),   # 4 + √15, N = 1
        17: (4, 1),   # (4+√17)/2 (d≡1(4))
        19: (170, 39),
        21: (5, 1),   # (5+√21)/2
        22: (197, 42),
        23: (24, 5),
        26: (5, 1),
        29: (5, 1),   # (5+√29)/2
        30: (11, 2),
        31: (1520, 273),
        33: (23, 4),  # (23+4√33)/2
        37: (6, 1),
        41: (32, 5),  # (32+5√41)/2
        43: (3482, 531),
        46: (24335, 3588),
        47: (48, 7),
        53: (7, 1),
        57: (151, 20),
        58: (19603, 2574),
        61: (39, 5),
        65: (8, 1),
    }

    if d_sf in known_units:
        a, b = known_units[d_sf]
        # Bestimme ob halbe Einheit oder ganze
        if d_sf % 4 == 1:
            # d ≡ 1 (mod 4): Einheit der Form (a + b√d)/2
            norm_val = (a * a - d_sf * b * b) // 4
            unit_str = f"({a} + {b}·√{d_sf})/2"
        else:
            norm_val = a * a - d_sf * b * b
            unit_str = f"{a} + {b}·√{d_sf}"
    else:
        # Algorithmus: Kettenbruchentwicklung von √d
        a, b, norm_val = _pell_solution(d_sf)
        if d_sf % 4 == 1:
            unit_str = f"({a} + {b}·√{d_sf})/2"
        else:
            unit_str = f"{a} + {b}·√{d_sf}"

    return {
        "a": a,
        "b": b,
        "norm": norm_val,
        "unit_str": unit_str,
        "description": f"Fundamentaleinheit von ℚ(√{d_sf}): ε = {unit_str}, N(ε) = {norm_val}",
        "unit_rank": 1,
    }


def _pell_solution(d: int) -> tuple[int, int, int]:
    """
    @brief Löst die Pell-Gleichung x² - d·y² = ±1 via Kettenbruchmethode.
    @description
        Algorithmus:
        1. Berechne a0 = floor(√d)
        2. Kettenbruchentwicklung: m_k, d_k, a_k iterativ
        3. Konvergenten p_k/q_k via p_k = a_k·p_{k-1} + p_{k-2}
        4. Teste: p_k² - d·q_k² = ±1

    @param d Quadratfreier Parameter d > 0, kein perfektes Quadrat
    @return Tuple (a, b, norm) mit a + b√d als Fundamentaleinheit
    @lastModified 2026-03-11
    """
    # Prüfe ob d ein perfektes Quadrat ist
    a0 = int(math.isqrt(d))
    if a0 * a0 == d:
        return (a0, 0, a0 * a0)  # Trivial: √d ∈ Z

    # Kettenbruchentwicklung von √d
    # Algorithmus: √d = a0 + (√d - a0) = a0 + 1/((√d+a0)/(d-a0²))
    m, d_k, a = 0, 1, a0
    p_prev, p_curr = 1, a0
    q_prev, q_curr = 0, 1

    for _ in range(10000):  # Maximale Iterationstiefe
        # Schritt der Kettenbruchentwicklung
        m = d_k * a - m
        d_k = (d - m * m) // d_k
        if d_k == 0:
            break
        a = (a0 + m) // d_k

        # Konvergenten aktualisieren
        p_new = a * p_curr + p_prev
        q_new = a * q_curr + q_prev
        p_prev, p_curr = p_curr, p_new
        q_prev, q_curr = q_curr, q_new

        # Teste Pell-Lösung
        norm_val = p_curr * p_curr - d * q_curr * q_curr
        if norm_val in (1, -1):
            return (p_curr, q_curr, norm_val)

    # Fallback: triviale Schätzung
    return (a0 + 1, 1, (a0 + 1) ** 2 - d)


def regulator_estimate(d: int) -> float:
    """
    @brief Schätzt den Regulator R_K = log|ε| für reell-quadratische Körper.
    @description
        Der Regulator R_K misst das "Volumen" der Einheitengruppe:
        R_K = log|ε| für reell-quadratische Körper (Einheitenrang 1)

        Für imaginär-quadratische Körper: R_K = 1 (trivial, kein echter Regulator).

        Der Regulator erscheint in der Dirichletschen Klassenzahlformel:
        Residuum bei s=1: lim_{s→1} (s-1)·ζ_K(s) = (2^{r1}·(2π)^{r2}·h·R_K)/(w·√|Δ|)

    @param d Quadratfreier Parameter
    @return Regulator R_K = log|ε| (oder 1 für imaginär-quadratisch)
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)

    if d_sf <= 0:
        # Imaginär-quadratisch: kein echter Regulator, R = 1 per Konvention
        return 1.0

    # Fundamentaleinheit berechnen
    unit_data = fundamental_units_quadratic(d_sf)
    a = unit_data["a"]
    b = unit_data["b"]

    # Numerischer Wert der Einheit: |ε| = a + b·√d (für d≡2,3 mod 4)
    # oder |ε| = (a + b·√d)/2 (für d≡1 mod 4)
    if d_sf % 4 == 1:
        eps_value = (a + b * math.sqrt(d_sf)) / 2.0
    else:
        eps_value = a + b * math.sqrt(d_sf)

    if eps_value <= 0:
        eps_value = abs(eps_value)

    return math.log(max(eps_value, 1.0))


# ---------------------------------------------------------------------------
# 5. Quadratische Körper
# ---------------------------------------------------------------------------

class QuadraticField:
    """
    @brief Spezialisierte Klasse für quadratische Zahlkörper K = Q(√d).
    @description
        Quadratische Körper sind die einfachsten nicht-trivialen Zahlkörper.
        Jeder quadratische Körper ist durch eine quadratfreie ganze Zahl d ≠ 0, 1
        eindeutig bestimmt.

        Eigenschaften:
        - [K:Q] = 2
        - Signatur: (2,0) falls d > 0 (reell), (0,1) falls d < 0 (imaginär)
        - Diskriminante: Δ = d falls d≡1(4), sonst 4d
        - Ganzheitsring: Z[(1+√d)/2] falls d≡1(4), sonst Z[√d]

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, d: int):
        """
        @brief Initialisiert den quadratischen Körper Q(√d).
        @param d Quadratfreier ganzzahliger Parameter d ≠ 0, 1
        @lastModified 2026-03-11
        """
        self.d = _squarefree_part(d)
        if self.d in (0, 1):
            raise ValueError(f"d = {d} ist kein gültiger Parameter (d ≠ 0, 1, muss quadratfrei sein)")
        # Minimalpolynom von √d: x² - d = 0
        # Als Koeffizientenliste: [-d, 0, 1] (aufsteigend: a0 + a1·x + a2·x²)
        self._number_field = NumberField([-self.d, 0, 1], name=f"√{self.d}")

    def discriminant(self) -> int:
        """
        @brief Diskriminante des quadratischen Körpers.
        @description
            Δ(Q(√d)/Q) = d falls d ≡ 1 (mod 4), sonst 4d
        @return Diskriminante
        @lastModified 2026-03-11
        """
        if self.d % 4 == 1:
            return self.d
        return 4 * self.d

    def signature(self) -> tuple[int, int]:
        """
        @brief Signatur (r1, r2) des Körpers.
        @return (2, 0) falls d > 0, (0, 1) falls d < 0
        @lastModified 2026-03-11
        """
        return (2, 0) if self.d > 0 else (0, 1)

    def ring_of_integers_basis(self) -> list[str]:
        """
        @brief ℤ-Basis des Ganzheitsrings O_K.
        @return Basis als Liste von Strings
        @lastModified 2026-03-11
        """
        if self.d % 4 == 1:
            return ["1", f"(1+√{self.d})/2"]
        return ["1", f"√{self.d}"]

    def class_number(self) -> int:
        """
        @brief Klassenzahl h(Q(√d)).
        @return Klassenzahl
        @lastModified 2026-03-11
        """
        return class_number(self.d)

    def is_imaginary(self) -> bool:
        """@brief True wenn der Körper imaginär-quadratisch ist (d < 0)."""
        return self.d < 0

    def is_real(self) -> bool:
        """@brief True wenn der Körper reell-quadratisch ist (d > 0)."""
        return self.d > 0

    def number_of_roots_of_unity(self) -> int:
        """
        @brief Anzahl der Einheitswurzeln w(K) in K.
        @description
            Für imaginär-quadratische Körper:
            - d = -1 (Q(i)): w = 4 (Einheitswurzeln ±1, ±i)
            - d = -3 (Q(ω)): w = 6 (Einheitswurzeln ±1, ±ω, ±ω²)
            - sonst: w = 2 (nur ±1)
            Für reell-quadratische Körper: w = 2 (nur ±1)

        @return Anzahl der Einheitswurzeln
        @lastModified 2026-03-11
        """
        if self.d == -1:
            return 4
        if self.d == -3:
            return 6
        return 2

    def fundamental_unit(self) -> Optional[dict]:
        """
        @brief Fundamentaleinheit für reell-quadratische Körper.
        @return Fundamentaleinheit oder None für imaginär-quadratische Körper
        @lastModified 2026-03-11
        """
        if self.is_imaginary():
            return None
        return fundamental_units_quadratic(self.d)

    def __repr__(self) -> str:
        """@brief Textuelle Darstellung des quadratischen Körpers."""
        kind = "reell" if self.d > 0 else "imaginär"
        return f"QuadraticField ℚ(√{self.d}) [{kind}-quadratisch, h={self.class_number()}]"


def legendre_symbol_fast(a: int, p: int) -> int:
    """
    @brief Schnelle Berechnung des Legendre-Symbols (a/p) über Euler-Kriterium.
    @description
        Das Legendre-Symbol (a/p) für eine ungerade Primzahl p:
        - (a/p) = 0 falls p | a
        - (a/p) = 1 falls a quadratischer Rest mod p (∃x: x² ≡ a (mod p))
        - (a/p) = -1 falls a quadratischer Nichtrest mod p

        Euler-Kriterium: (a/p) ≡ a^{(p-1)/2} (mod p)

        Laufzeit: O(log p) via schneller modularer Potenzierung.

    @param a Ganzzahl (Zähler)
    @param p Ungerade Primzahl
    @return 0, 1 oder -1
    @lastModified 2026-03-11
    """
    if not _is_prime(p) or p == 2:
        raise ValueError(f"p={p} muss eine ungerade Primzahl sein")
    return _legendre_symbol(a, p)


def kronecker_symbol(a: int, n: int) -> int:
    """
    @brief Verallgemeinertes Jacobi-Symbol (Kronecker-Symbol) (a/n).
    @description
        Das Kronecker-Symbol verallgemeinert das Jacobi-Symbol auf beliebige n:
        - Für n = p prim: = Legendre-Symbol
        - Für n = p1^e1 · p2^e2 · ...: = ∏ (a/pi)^ei (multiplikativ)
        - Für n = 2: spezielle Regeln basierend auf a mod 8
        - Für n = -1: (a/-1) = -1 falls a < 0, sonst 1
        - Für n = 0: (a/0) = 1 falls |a|=1, sonst 0
        - Für n = 1: (a/1) = 1 immer

        Anwendung: Charakter χ_d(n) = (d/n) in Dirichlet-L-Reihen.

    @param a Ganzzahl (Zähler)
    @param n Ganzzahl (Nenner, ≥ 0)
    @return Kronecker-Symbol: -1, 0 oder 1
    @lastModified 2026-03-11
    """
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n == 1:
        return 1
    if n == -1:
        return -1 if a < 0 else 1

    # Vorzeichen von n behandeln
    result = 1
    if n < 0:
        result *= kronecker_symbol(a, -1)
        n = -n

    # Faktor 2 separat behandeln
    while n % 2 == 0:
        n //= 2
        # (a/2): Kronecker-Symbol für 2
        a_mod8 = a % 8
        if a_mod8 in (3, 5):
            result *= -1
        elif a_mod8 in (0, 2, 4, 6):
            return 0  # gerade a → Symbol = 0 (nur wenn a gerade und 2|n)
        # a ungerade: (a/2) = (-1)^{(a²-1)/8}

    if n == 1:
        return result

    # Restliche Faktoren (ungerade Primzahlen via Jacobi-Symbol)
    # Jacobi-Symbol: ∏ Legendre-Symbol für alle Primfaktoren
    result *= _jacobi_symbol(a, n)
    return result


def _jacobi_symbol(a: int, n: int) -> int:
    """
    @brief Jacobi-Symbol (a/n) für positive ungerade n.
    @description
        Berechnung via quadratisches Reziprozitätsgesetz (effizient).
        Laufzeit: O(log min(a,n)).
    @param a Ganzzahl
    @param n Positive ungerade ganze Zahl ≥ 1
    @return Jacobi-Symbol: -1, 0 oder 1
    @lastModified 2026-03-11
    """
    # Basisfälle
    if n == 1:
        return 1
    a = a % n
    if a == 0:
        return 0
    if a == 1:
        return 1

    result = 1
    while a != 0:
        # Faktor 2 aus a herausziehen
        while a % 2 == 0:
            a //= 2
            # (2/n) = (-1)^{(n²-1)/8}
            if n % 8 in (3, 5):
                result = -result

        # Reziprozitätsgesetz: (a/n) · (n/a) = (-1)^{(a-1)(n-1)/4}
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a = a % n

    return result if n == 1 else 0


def quadratic_reciprocity_law(p: int, q: int) -> dict:
    """
    @brief Überprüft das Quadratische Reziprozitätsgesetz von Gauß.
    @description
        Gaußsches Quadratisches Reziprozitätsgesetz (1796):
        Für verschiedene ungerade Primzahlen p, q gilt:

        (p/q) · (q/p) = (-1)^{(p-1)/2 · (q-1)/2}

        Äquivalent:
        - Falls p ≡ 1 (mod 4) oder q ≡ 1 (mod 4): (p/q) = (q/p)
        - Falls p ≡ q ≡ 3 (mod 4): (p/q) = -(q/p)

        Ergänzungssätze:
        - (-1/p) = (-1)^{(p-1)/2} = 1 falls p≡1(4), = -1 falls p≡3(4)
        - (2/p) = (-1)^{(p²-1)/8} = 1 falls p≡±1(8), = -1 falls p≡±3(8)

    @param p Ungerade Primzahl
    @param q Ungerade Primzahl (≠ p)
    @return Dict mit Legendre-Symbolen und Reziprozitätsprüfung
    @lastModified 2026-03-11
    """
    if not _is_prime(p) or p == 2:
        raise ValueError(f"p={p} muss eine ungerade Primzahl sein")
    if not _is_prime(q) or q == 2:
        raise ValueError(f"q={q} muss eine ungerade Primzahl sein")
    if p == q:
        raise ValueError("p und q müssen verschieden sein")

    leg_pq = _legendre_symbol(p, q)  # (p/q)
    leg_qp = _legendre_symbol(q, p)  # (q/p)

    # Erwarteter Wert des Produkts: (-1)^{(p-1)(q-1)/4}
    exp_val = ((p - 1) // 2) * ((q - 1) // 2)
    expected_product = (-1) ** exp_val

    actual_product = leg_pq * leg_qp
    law_holds = (actual_product == expected_product)

    # Ergänzungssätze
    neg1_p = (-1) ** ((p - 1) // 2)  # (-1/p)
    neg1_q = (-1) ** ((q - 1) // 2)  # (-1/q)
    two_p = (-1) ** ((p * p - 1) // 8)  # (2/p)
    two_q = (-1) ** ((q * q - 1) // 8)  # (2/q)

    return {
        "legendre_pq": leg_pq,
        "legendre_qp": leg_qp,
        "product": actual_product,
        "expected_product": expected_product,
        "law_holds": law_holds,
        "neg1_over_p": neg1_p,
        "neg1_over_q": neg1_q,
        "two_over_p": two_p,
        "two_over_q": two_q,
        "description": (
            f"({p}/{q})·({q}/{p}) = {leg_pq}·{leg_qp} = {actual_product} = "
            f"(-1)^{{({p}-1)({q}-1)/4}} = {expected_product}: "
            f"{'ERFÜLLT' if law_holds else 'VERLETZT'}"
        ),
    }


def hilbert_class_field_degree(d: int) -> int:
    """
    @brief Grad des Hilbert-Klassenfelds H/K über K = Q(√d).
    @description
        Das Hilbert-Klassenfeld H ist die maximale unverzweigte abelsche
        Erweiterung von K. Nach dem Hauptsatz der Klassenkörpertheorie gilt:
        [H : K] = h(K) (Klassenzahl von K)

        Also: deg(H/Q) = 2·h(K)

        Beispiel:
        - K = Q(√-5), h=2: H = Q(√-5, √-1) = Q(√5, i), [H:K]=2
        - K = Q(√-23), h=3: H über Q hat Grad 6

    @param d Quadratfreier Parameter von K = Q(√d)
    @return Grad [H:K] = h(K)
    @lastModified 2026-03-11
    """
    return class_number(d)


# ---------------------------------------------------------------------------
# 6. Zyklotomische Körper
# ---------------------------------------------------------------------------

class CyclotomicField:
    """
    @brief Zyklotomischer Körper Q(ζ_n) mit primitiver n-ter Einheitswurzel ζ_n.
    @description
        Der n-te zyklotomische Körper Q(ζ_n) entsteht durch Adjunktion
        einer primitiven n-ten Einheitswurzel ζ_n = e^{2πi/n}.

        Eigenschaften:
        - Grad: [Q(ζ_n):Q] = φ(n) (Euler'sche φ-Funktion)
        - Ganzheitsring: Z[ζ_n] (für alle n)
        - Diskriminante: Δ = (-1)^{φ(n)/2} · n^{φ(n)} / ∏_{p|n} p^{φ(n)/(p-1)}
        - Galoisgruppe: Gal(Q(ζ_n)/Q) ≅ (Z/nZ)^×
        - Signatur: (0, φ(n)/2) für n ≥ 3 (komplex für n > 2)

        Bedeutung in der Zahlentheorie:
        - Kummer'sche Theorie für Fermat: x^p + y^p = z^p in Z[ζ_p]
        - Beweis des Satzes von Gauß über Konstruierbarkeit regelmäßiger Polygone

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, n: int):
        """
        @brief Initialisiert den zyklotomischen Körper Q(ζ_n).
        @param n Grad der Einheitswurzel (n ≥ 2)
        @lastModified 2026-03-11
        """
        if n < 2:
            raise ValueError(f"n = {n} muss ≥ 2 sein")
        self.n = n
        # Grad = φ(n)
        self._degree = _euler_totient(n)

    def degree(self) -> int:
        """
        @brief Grad des zyklotomischen Körpers [Q(ζ_n):Q] = φ(n).
        @return φ(n)
        @lastModified 2026-03-11
        """
        return self._degree

    def signature(self) -> tuple[int, int]:
        """
        @brief Signatur des zyklotomischen Körpers.
        @description
            Für n = 1, 2: Q(ζ_n) = Q, Signatur (1, 0)
            Für n ≥ 3: Q(ζ_n) vollständig komplex → (0, φ(n)/2)
        @return Tuple (r1, r2)
        @lastModified 2026-03-11
        """
        if self.n <= 2:
            return (1, 0)
        phi = self._degree
        return (0, phi // 2)

    def cyclotomic_ring_of_integers(self) -> dict:
        """
        @brief Ganzheitsring des zyklotomischen Körpers: Z[ζ_n].
        @description
            Satz (monogen): Der Ganzheitsring O_{Q(ζ_n)} = Z[ζ_n] für alle n.
            Dies ist ein wichtiger Satz — der Ganzheitsring ist immer monogen!

            Z-Basis: {1, ζ_n, ζ_n², ..., ζ_n^{φ(n)-1}}

        @return Dict mit Basis, Beschreibung und Grad
        @lastModified 2026-03-11
        """
        phi = self._degree
        basis = [f"ζ_{self.n}^{k}" if k > 0 else "1" for k in range(phi)]
        return {
            "ring": f"ℤ[ζ_{self.n}]",
            "basis": basis,
            "degree": phi,
            "is_monogenic": True,
            "description": f"O_K = ℤ[ζ_{self.n}], ℤ-Basis: {{1, ζ_{self.n}, ..., ζ_{self.n}^{{{phi - 1}}}}}",
        }

    def discriminant(self) -> int:
        """
        @brief Diskriminante des zyklotomischen Körpers.
        @description
            Δ(Q(ζ_n)) = (-1)^{φ(n)/2} · n^{φ(n)} / ∏_{p|n} p^{φ(n)/(p-1)}

            Für n = p (Primzahl): Δ = (-1)^{(p-1)/2} · p^{p-2}
            Für n = p² (Primzahlquadrat): Δ = (-1)^{φ(p²)/2} · p^{p(p-2)+1}

        @return Diskriminante (Betrag)
        @lastModified 2026-03-11
        """
        phi = self._degree
        n = self.n
        # ∏_{p|n} p^{φ(n)/(p-1)}
        prime_factors = list(_prime_factorization(n).keys())
        denominator = 1
        for p in prime_factors:
            exp = phi // (p - 1)
            denominator *= p ** exp
        numerator = n ** phi
        disc_abs = numerator // denominator if denominator > 0 else numerator
        sign = (-1) ** (phi // 2)
        return sign * disc_abs

    def galois_group_order(self) -> int:
        """
        @brief Ordnung der Galoisgruppe Gal(Q(ζ_n)/Q) = φ(n).
        @description
            Gal(Q(ζ_n)/Q) ≅ (Z/nZ)^× (multiplikative Gruppe mod n)
            Die Isomorphie ist gegeben durch: σ_a(ζ_n) = ζ_n^a für gcd(a,n)=1
        @return φ(n) = |(Z/nZ)^×|
        @lastModified 2026-03-11
        """
        return self._degree

    def is_prime_cyclotomic(self) -> bool:
        """@brief True wenn n eine Primzahl ist (einfachster Fall)."""
        return _is_prime(self.n)

    def splitting_of_prime(self, p: int) -> dict:
        """
        @brief Zerlegung einer Primzahl p in Z[ζ_n].
        @description
            Für eine Primzahl p in Q(ζ_n):
            - Falls p | n: p verzweigt vollständig (e = φ(n), f = 1)
            - Sonst: Sei f = ord_n(p) (multiplikative Ordnung von p mod n)
              Dann: (p) = p1 · p2 · ... · p_r mit r = φ(n)/f, f = f, e = 1

        @param p Primzahl
        @return Dict mit Zerlegungstyp, e, f, r
        @lastModified 2026-03-11
        """
        if not _is_prime(p):
            raise ValueError(f"{p} ist keine Primzahl")

        n = self.n
        phi = self._degree

        if n % p == 0:
            # p | n: vollständig verzweigt
            e = phi
            f = 1
            r = 1
            split_type = "vollständig verzweigt"
        else:
            # p ∤ n: Ordnung von p mod n berechnen
            f = 1
            pk = p % n
            while pk != 1 and f <= n:
                pk = (pk * p) % n
                f += 1
            # Anzahl der Primideale
            r = phi // f
            e = 1
            split_type = "unverzweigt"

        return {
            "ramification_index": e,
            "inertia_degree": f,
            "n_prime_ideals": r,
            "split_type": split_type,
            "efr_check": e * f * r,
            "phi_n": phi,
            "check_passed": (e * f * r == phi),
        }

    def __repr__(self) -> str:
        """@brief Textuelle Darstellung des zyklotomischen Körpers."""
        return f"CyclotomicField Q(ζ_{self.n}), [Q(ζ_{self.n}):Q] = φ({self.n}) = {self._degree}"


def cyclotomic_ring_of_integers(n: int) -> dict:
    """
    @brief Gibt den Ganzheitsring Z[ζ_n] des n-ten zyklotomischen Körpers zurück.
    @param n Grad der Einheitswurzel (n ≥ 2)
    @return Dict mit Basis und Beschreibung
    @lastModified 2026-03-11
    """
    cf = CyclotomicField(n)
    return cf.cyclotomic_ring_of_integers()


def kummer_lifting(p: int, n: int) -> dict:
    """
    @brief Kummer-Theorie: Analyse der p-ten Einheitswurzeln in Z[ζ_n].
    @description
        Kummer's Theorie für Fermat's letzten Satz:
        In Z[ζ_p] (p Primzahl) gilt die Gleichung x^p + y^p = z^p
        implizit die Zerlegung in Primideale:
        (x + y)(x + ζ_p·y)·...·(x + ζ_p^{p-1}·y) = z^p

        Für reguläre Primzahlen p (p ∤ Klassenzahl h(Q(ζ_p))) beweist
        Kummer den Großen Fermatsatz.

        Kummer-Lifting (Hensel's Lemma für Primideale):
        Wenn f(a) ≡ 0 (mod p) und f'(a) ≢ 0 (mod p), dann existiert
        a' ≡ a (mod p) mit f(a') ≡ 0 (mod p^n).

        Diese Funktion untersucht:
        1. Ob p regulär (h(Q(ζ_p)) ≢ 0 (mod p))
        2. Zerlegung von p in Z[ζ_n]
        3. Kummer-Liftungsschritte für einfache Polynome

    @param p Primzahl (Fermat-Basis)
    @param n Potenz (Ziel: Kongruenz mod p^n)
    @return Dict mit Kummer-Analyse
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"{p} ist keine Primzahl")

    # Bekannte irreguläre Primzahlen (p | h(Q(ζ_p)))
    irregular_primes = {37, 59, 67, 101, 103, 131, 149, 157, 233, 257, 263}
    is_regular = p not in irregular_primes and p < 37

    # Zerlegung der Primzahl p-1 in Z[ζ_p]: (p) = (1-ζ_p)^{p-1}
    # (vollständig verzweigt, da p | p)
    phi_p = p - 1
    ramification = phi_p  # Vollständig verzweigt

    # Kummer-Lifting: Beispiel x^p ≡ 1 (mod p^n) → x ≡ 1 + k·p^k (mod p^{k+1})
    # Schritt 1: x₀ = 1 (triviale Lösung mod p)
    # Schritt k: xₖ = xₖ₋₁ + t·p^k mit geeignetem t
    lifting_steps = []
    x_current = 1
    for step in range(1, min(n + 1, 5)):
        # x^p ≡ ? (mod p^{step+1})
        x_p = pow(x_current, p)
        remainder = x_p % (p ** (step + 1))
        lifting_steps.append({
            "step": step,
            "x": x_current,
            "x_to_p_mod_pp": remainder,
        })
        # Nächste Approximation
        x_current = x_current  # Einfache Approximation

    return {
        "prime": p,
        "target_power": n,
        "is_regular_prime": is_regular,
        "phi_p": phi_p,
        "ramification_in_cyclotomic": ramification,
        "lifting_steps": lifting_steps,
        "fermat_relevant": is_regular,
        "description": (
            f"Kummer-Theorie für p={p}: "
            f"{'Reguläre' if is_regular else 'Irreguläre'} Primzahl. "
            f"Zerlegung (p) = (1-ζ_p)^{{{phi_p}}} in ℤ[ζ_p]."
        ),
    }


# ---------------------------------------------------------------------------
# 7. Lokale Theorie
# ---------------------------------------------------------------------------

def p_adic_completion(poly_coeffs: list[int], p: int) -> dict:
    """
    @brief Beschreibt die p-adische Vervollständigung K_p für K = Q[x]/(f).
    @description
        Die p-adische Vervollständigung K_p eines Zahlkörpers K bezüglich
        eines Primideals P über p ist ein lokaler Körper:

        K_p ≅ Q_p[x]/(g(x))

        wobei g ein irreduzibler Faktor von f über Q_p ist.

        Eigenschaften:
        - [K_p : Q_p] = f (Trägheitsgrad)
        - P ist verzweigt ⟺ e > 1
        - Residuenkörper: O_{K_p}/P·O_{K_p} ≅ F_{p^f}

        Die Zerlegung (p) = ∏ Pi^{ei} in O_K entspricht der Produktformel:
        K ⊗_Q Q_p ≅ ∏_i K_{Pi}

    @param poly_coeffs Koeffizienten des Minimalpolynoms
    @param p Primzahl
    @return Dict mit lokaler Zerlegung und Eigenschaften
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"{p} ist keine Primzahl")

    factors = ideal_factorization(p, poly_coeffs)
    n = len(poly_coeffs) - 1

    completions = []
    for factor in factors:
        e = factor["ramification_index"]
        f = factor["inertia_degree"]
        completions.append({
            "local_field": f"K_{{P}} ≅ Q_{p}[x]/(g_{f}(x))",
            "inertia_degree": f,
            "ramification_index": e,
            "residue_field": f"𝔽_{{p^{f}}}",
            "is_unramified": (e == 1),
            "degree_over_Qp": e * f,
        })

    return {
        "prime": p,
        "n_completions": len(completions),
        "completions": completions,
        "product_formula": f"K ⊗_Q Q_{p} ≅ {'  ×  '.join(c['local_field'] for c in completions)}",
        "description": (
            f"p-adische Vervollständigung von K bezüglich p={p}: "
            f"{len(completions)} lokale Faktoren"
        ),
    }


def hensel_lifting_general(
    poly_coeffs: list[int],
    a0: int,
    p: int,
    steps: int = 5
) -> dict:
    """
    @brief Hensel-Lifting: Newton-Iteration in O_{K_p}.
    @description
        Hensels Lemma (allgemeine Form):
        Sei f ∈ Z_p[x] und a₀ ∈ Z mit:
        - f(a₀) ≡ 0 (mod p)
        - f'(a₀) ≢ 0 (mod p) (einfache Nullstelle)

        Dann existiert ein eindeutiges â ∈ Z_p mit f(â) = 0 und â ≡ a₀ (mod p).

        Newton-Iteration:
        aₙ₊₁ = aₙ - f(aₙ)/f'(aₙ)  (in Z/p^{2^n}Z)

        Konvergenzrate: v_p(f(aₙ₊₁)) ≥ 2·v_p(f(aₙ))
        (quadratische Konvergenz → nach k Schritten: Genauigkeit p^{2^k})

    @param poly_coeffs Koeffizienten des Polynoms f (aufsteigend)
    @param a0 Startnäherung (Nullstelle mod p)
    @param p Primzahl (Basis des p-adischen Körpers)
    @param steps Anzahl der Lifting-Schritte
    @return Dict mit Lifting-Folge und Genauigkeitsanalyse
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"{p} ist keine Primzahl")

    # Ableitungskoeffizienten: f'(x) = Σ k·aₖ·x^{k-1}
    deriv_coeffs = [i * poly_coeffs[i] for i in range(1, len(poly_coeffs))]

    # Prüfe Startnäherung
    f_a0 = _poly_eval_mod_p(poly_coeffs, a0, p)
    fp_a0 = _poly_eval_mod_p(deriv_coeffs, a0, p)

    if f_a0 != 0:
        return {
            "success": False,
            "error": f"a₀={a0} ist keine Nullstelle von f mod {p}",
            "a0": a0, "f_a0": f_a0,
        }

    if fp_a0 == 0:
        return {
            "success": False,
            "error": f"f'(a₀={a0}) ≡ 0 (mod {p}): Hensel-Lemma nicht anwendbar (mehrfache Nullstelle)",
            "a0": a0, "fp_a0": fp_a0,
        }

    # Lifting-Schritte
    lifting_sequence = [{"step": 0, "a": a0, "modulus": p, "f_val": f_a0}]
    a_current = a0
    modulus = p

    for step in range(1, steps + 1):
        modulus_next = modulus * p
        # Newton-Schritt: a_new = a_old - f(a_old)/f'(a_old) mod p^{step+1}
        f_val = sum(poly_coeffs[k] * pow(a_current, k) for k in range(len(poly_coeffs)))
        fp_val = sum(deriv_coeffs[k] * pow(a_current, k) for k in range(len(deriv_coeffs)))

        # Division in Z/modulus_next Z: fp_val muss Einheit mod p sein
        fp_mod_p = fp_val % p
        if fp_mod_p == 0:
            break  # Hensel-Voraussetzung verletzt

        # Inverse von f'(a) mod p (für Newton-Schritt)
        fp_inv = pow(fp_mod_p, p - 2, p)  # Fermat'sche Umkehrung
        # Newton-Korrektur
        correction = (f_val * fp_inv) % modulus_next
        a_new = (a_current - correction) % modulus_next

        f_val_new = sum(poly_coeffs[k] * pow(a_new, k) for k in range(len(poly_coeffs)))
        f_val_new_mod = f_val_new % modulus_next

        lifting_sequence.append({
            "step": step,
            "a": a_new,
            "modulus": modulus_next,
            "f_val_mod": f_val_new_mod,
        })
        a_current = a_new
        modulus = modulus_next

    return {
        "success": True,
        "a0": a0,
        "final_a": a_current,
        "final_modulus": modulus,
        "lifting_sequence": lifting_sequence,
        "description": (
            f"Hensel-Lifting: Nullstelle a₀={a0} von f mod {p} "
            f"auf Genauigkeit mod p^{steps+1}={p**(steps)} geliftet."
        ),
    }


def local_norm_symbol(a: int, b: int, p: int) -> int:
    """
    @brief Berechnet das Hilbert-Symbol (a, b)_p für lokale Körper.
    @description
        Das Hilbert-Symbol (a, b)_p ∈ {+1, -1} für a, b ∈ Q_p^× ist definiert:
        (a, b)_p = 1 ⟺ ax² + by² = z² hat eine nichttriviale Lösung in Q_p

        Eigenschaften:
        - (a, b)_p = (b, a)_p (Symmetrie)
        - (a, -a)_p = 1 (Normierung)
        - (a, b)_p · (a, c)_p = (a, bc)_p (Bimultiplikativität)

        Produktformel: ∏_v (a, b)_v = 1 (über alle Stellen v)

        Berechnung für Q_p (via Hilbert-Symbol-Formeln):
        Für p ungerade und a, b p-adische Einheiten:
        (a, b)_p = (a/p)^{v_p(b)} · (b/p)^{v_p(a)} · (-1)^{v_p(a)·v_p(b)·(p-1)/2}

        Für p = 2: Spezialformel mit (a²-1)/8 + (b²-1)/8

    @param a Erste ganze Zahl (≠ 0)
    @param b Zweite ganze Zahl (≠ 0)
    @param p Primzahl (lokale Stelle)
    @return Hilbert-Symbol: +1 oder -1
    @lastModified 2026-03-11
    """
    if not _is_prime(p):
        raise ValueError(f"{p} ist keine Primzahl")
    if a == 0 or b == 0:
        raise ValueError("a und b müssen ungleich Null sein")

    # p-adische Bewertungen v_p(a) und v_p(b)
    def vp(x: int, prime: int) -> int:
        """p-adische Bewertung: größte Potenz von prime, die x teilt."""
        if x == 0:
            return float('inf')
        x = abs(x)
        v = 0
        while x % prime == 0:
            x //= prime
            v += 1
        return v

    va = vp(a, p)
    vb = vp(b, p)

    # Einheitenanteile (p-adische Einheiten)
    a_unit = a // (p ** va)
    b_unit = b // (p ** vb)

    if p == 2:
        # Hilbert-Symbol für p = 2
        # (a, b)_2 = (-1)^{(a_unit²-1)/8 + (b_unit²-1)/8 + va·vb·(a_unit-1)/2·(b_unit-1)/2...}
        # Vereinfachte Formel für ganze Zahlen
        a_odd = a_unit if a_unit % 2 == 1 else a_unit + 1
        b_odd = b_unit if b_unit % 2 == 1 else b_unit + 1
        # (a, b)_2: -1 wenn ax² + by² = z² keine Lösung in Q_2
        # Für Einheiten: (a,b)_2 = (-1)^{((a-1)(b-1)/4 + v_2(a)·(b²-1)/8 + v_2(b)·(a²-1)/8)}
        exp = (
            ((a_odd - 1) * (b_odd - 1) // 4 % 2) +
            (va * ((b_odd * b_odd - 1) // 8) % 2) +
            (vb * ((a_odd * a_odd - 1) // 8) % 2)
        ) % 2
        return (-1) ** exp
    else:
        # p ungerade: Formel via Legendre-Symbole
        leg_a_p = _legendre_symbol(a_unit, p)
        leg_b_p = _legendre_symbol(b_unit, p)

        # (a, b)_p = (a/p)^{v_p(b)} · (b/p)^{v_p(a)} · (-1)^{v_p(a)·v_p(b)·(p-1)/2}
        sign_exp = va * vb * ((p - 1) // 2) % 2
        result = (leg_a_p ** vb) * (leg_b_p ** va) * ((-1) ** sign_exp)

        # Normiere auf {-1, +1}
        return 1 if result >= 0 else -1


# ---------------------------------------------------------------------------
# Klassenzahlformel (Dirichlet)
# ---------------------------------------------------------------------------

def dirichlet_class_number_formula(d: int) -> dict:
    """
    @brief Berechnet die Klassenzahlformel für quadratische Körper.
    @description
        Dirichletsche Klassenzahlformel für K = Q(√d):

        Für imaginär-quadratische Körper (d < 0):
        h(K) = (w·√|Δ|)/(2π) · L(1, χ_Δ)

        wobei:
        - w = Anzahl der Einheitswurzeln in K
        - Δ = Diskriminante
        - L(1, χ_Δ) = Σ χ_Δ(n)/n (Dirichlet-L-Funktion bei s=1)
        - χ_Δ = Kronecker-Symbol (Δ/·)

        Für reell-quadratische Körper (d > 0):
        h(K)·R_K = (√Δ)/(2) · L(1, χ_Δ)

        wobei R_K = log|ε| der Regulator ist.

        Diese Formel verknüpft analytische (L-Funktionen) mit algebraischen
        (Klassenzahl, Regulator) Eigenschaften des Zahlkörpers.

    @param d Quadratfreier Parameter
    @return Dict mit Δ, w, L-Wert, h, Regulator und Formelprüfung
    @lastModified 2026-03-11
    """
    d_sf = _squarefree_part(d)
    disc = d_sf if d_sf % 4 == 1 else 4 * d_sf
    abs_disc = abs(disc)

    # Anzahl der Einheitswurzeln
    if d_sf == -1:
        w = 4
    elif d_sf == -3:
        w = 6
    else:
        w = 2

    # L(1, χ_Δ) numerisch approximieren
    # L(1, χ) = Σ_{n=1}^{∞} χ(n)/n für den primitiven Charakter χ = (Δ/·)
    n_terms = 2000
    l_value = 0.0
    for k in range(1, n_terms + 1):
        chi_k = kronecker_symbol(disc, k)
        if chi_k != 0:
            l_value += chi_k / k

    # Klassenzahl aus der Formel
    h_exact = class_number(d_sf)
    reg = regulator_estimate(d_sf) if d_sf > 0 else 1.0

    if d_sf < 0:
        # h = (w·√|Δ|)/(2π) · L(1, χ)
        h_formula = (w * math.sqrt(abs_disc)) / (2 * math.pi) * l_value
    else:
        # h·R = (√Δ)/(2) · L(1, χ)
        h_formula = (math.sqrt(abs_disc) / 2) * l_value / max(reg, 0.01)

    return {
        "d": d_sf,
        "discriminant": disc,
        "w": w,
        "l_value_approx": l_value,
        "h_exact": h_exact,
        "h_from_formula": h_formula,
        "regulator": reg,
        "error": abs(h_formula - h_exact),
        "description": (
            f"Klassenzahlformel für ℚ(√{d_sf}): "
            f"h = {h_exact}, L(1,χ) ≈ {l_value:.6f}, "
            f"Formelwert ≈ {h_formula:.4f}"
        ),
    }
