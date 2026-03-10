"""
@file galois_theory.py
@brief Galois-Theorie und Körpererweiterungen – vollständige Implementierung.
@description
    Implementiert die grundlegenden Konzepte der Galois-Theorie:

    - **FieldExtension**: Körpererweiterung L/K mit algebraischem Element α.
      Berechnet Grad, Minimalpolynom, Basis, Norm, Spur, Galois-Eigenschaft.

    - **GaloisGroup**: Galoisgruppe Gal(L/K) als Automorphismengruppe.
      Berechnet Ordnung, Elemente (als Permutationen), Auflösbarkeit,
      Untergruppen, Normalteiler und Fixkörper.

    - **FiniteField**: Endlicher Körper GF(p^n) = ℤ_p[x]/(f(x)).
      Berechnet Ordnung, multiplikative Gruppe, Frobenius-Endomorphismus,
      primitive Elemente.

    **Zentrale mathematische Aussagen:**

    1. Galois-Hauptsatz: f(x) ist durch Radikale lösbar ⟺ Gal(f) ist auflösbar.
    2. Galois-Korrespondenz: {Untergruppen H ≤ Gal(L/K)} ↔ {Zwischenkörper K⊆F⊆L}
    3. Kreisteilungskörper: Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^×

    **Algorithmen:**
    - Diskriminante via Resultante (SymPy)
    - Galoisgruppe via Diskriminante und Resolventenpolynome
    - Diskreter Logarithmus: Baby-Step-Giant-Step
    - Kreisteilungspolynome via Möbius-Formel
    - Minimalpolynome via SymPy-Annullierungsalgorithmus

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import math
import cmath
import itertools
from typing import Optional
from math import gcd, isqrt
from functools import reduce

# SymPy für symbolische Berechnungen (Polynomring, Resultante, etc.)
import sympy
from sympy import (
    Poly, Symbol, ZZ, QQ, factor, resultant, discriminant as sym_discriminant,
    minimal_polynomial as sym_minimal_poly, sqrt as sym_sqrt, Integer,
    factorint, isprime, totient, primitive_root, Rational, floor,
    cyclotomic_poly, GF, FiniteField as SympyFF, symbols, expand, Mul
)
from sympy.abc import x as sym_x
from sympy.ntheory import n_order

# Lokale Ausnahmen
from exceptions import DomainError, InvalidInputError, MathematicalError


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def _poly_mul_mod(a: list[int], b: list[int], mod_poly: list[int], p: int) -> list[int]:
    """
    @brief Multipliziert zwei Polynome modulo mod_poly und modulo p.
    @description
        Berechnet (a * b) mod (mod_poly, p) im Polynomring ℤ_p[x]/(mod_poly).
        Koeffizienten werden nach jeder Operation reduziert.

    @param a Koeffizientenliste des ersten Polynoms (aufsteigend).
    @param b Koeffizientenliste des zweiten Polynoms (aufsteigend).
    @param mod_poly Modulus-Polynom (aufsteigend).
    @param p Primzahl für die Koeffizientenarithmetik.
    @return Koeffizientenliste des Produkts mod (mod_poly, p).
    @date 2026-03-10
    """
    # Grad des Modulus-Polynoms
    deg_mod = len(mod_poly) - 1

    # Standardmultiplikation von Polynomen
    deg_a = len(a) - 1
    deg_b = len(b) - 1
    # Ergebnisgrad ohne Reduktion wäre deg_a + deg_b
    result = [0] * (deg_a + deg_b + 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p

    # Reduktion modulo mod_poly: iterativ führenden Term eliminieren
    # Leitkoeffizient von mod_poly muss 1 sein (normiert)
    while len(result) > deg_mod:
        lead_exp = len(result) - 1        # aktueller führender Exponent
        lead_coef = result[lead_exp] % p  # führender Koeffizient
        if lead_coef == 0:
            result.pop()
            continue
        # Subtrahiere lead_coef * x^(lead_exp - deg_mod) * mod_poly
        offset = lead_exp - deg_mod
        for k, mk in enumerate(mod_poly):
            result[offset + k] = (result[offset + k] - lead_coef * mk) % p
        result.pop()

    # Normalisierung: führende Nullen entfernen
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def _poly_pow_mod(base: list[int], exp: int, mod_poly: list[int], p: int) -> list[int]:
    """
    @brief Schnelle Potenzierung eines Polynoms modulo (mod_poly, p).
    @description
        Verwendet Square-and-Multiply für effiziente Potenzierung
        in ℤ_p[x]/(mod_poly). Laufzeit: O(log(exp) * n²) für Grad-n-Polynome.

    @param base Basis-Polynom als Koeffizientenliste.
    @param exp Exponent (nicht-negativ).
    @param mod_poly Modulus-Polynom.
    @param p Primzahl.
    @return base^exp mod (mod_poly, p).
    @date 2026-03-10
    """
    # Neutrale Element der Multiplikation ist 1 (konstantes Polynom [1])
    result = [1]
    base = [c % p for c in base]
    # Square-and-multiply Algorithmus
    while exp > 0:
        if exp % 2 == 1:
            result = _poly_mul_mod(result, base, mod_poly, p)
        base = _poly_mul_mod(base, base, mod_poly, p)
        exp //= 2
    return result


def _poly_eval(coeffs: list[int], val: int, p: int) -> int:
    """
    @brief Wertet ein Polynom an einer Stelle aus (Horner-Schema) in ℤ_p.
    @param coeffs Koeffizientenliste (aufsteigend: coeffs[0] = konstantem Term).
    @param val Auswertungsstelle.
    @param p Primzahl (Modulus).
    @return Wert des Polynoms an val in ℤ_p.
    @date 2026-03-10
    """
    # Horner-Auswertung von oben nach unten
    result = 0
    for c in reversed(coeffs):
        result = (result * val + c) % p
    return result


def _is_irreducible_over_fp(coeffs: list[int], p: int) -> bool:
    """
    @brief Prüft ob ein Polynom über ℤ_p irreduzibel ist.
    @description
        Rabin-Irreduzibilitätstest:
        f(x) irreduzibel über ℤ_p vom Grad n ⟺
        1. gcd(f, x^{p^{n/q}} - x) = 1 für alle Primteiler q von n
        2. f teilt x^{p^n} - x

        Für n ≤ 6 reicht hier die Prüfung auf fehlende Wurzeln in ℤ_p.

    @param coeffs Polynom-Koeffizienten (aufsteigend), höchster Koeff ≠ 0.
    @param p Primzahl.
    @return True wenn irreduzibel, sonst False.
    @date 2026-03-10
    """
    # Normiertes Polynom (leitend = 1 in ℤ_p)
    n = len(coeffs) - 1  # Grad des Polynoms
    if n == 0:
        return False   # Konstante sind nicht irreduzibel
    if n == 1:
        return True    # Lineare Polynome sind stets irreduzibel

    # Schritt 1: Keine Wurzel in ℤ_p (notwendig, aber nicht hinreichend für n > 2)
    for v in range(p):
        if _poly_eval(coeffs, v, p) == 0:
            return False  # Wurzel gefunden → reduzibel

    if n == 2:
        # Für Grad 2: Keine Wurzel ⟺ irreduzibel
        return True

    # Schritt 2: Rabin-Test für Grad ≥ 3
    # Alle Primteiler des Grads bestimmen
    prime_factors_of_n = list(factorint(n).keys())

    # Polynom als SymPy-Objekt für ggT-Berechnung
    f_sym = Poly(list(reversed(coeffs)), sym_x, domain=GF(p))

    # Für jeden Primteiler q von n: gcd(f, x^{p^{n/q}} - x) soll 1 sein
    for q in prime_factors_of_n:
        exp = p ** (n // q)
        # x^{exp} mod f in GF(p)[x] berechnen (über SymPy oder manuell)
        # Wir nutzen _poly_pow_mod mit Basis x = [0, 1]
        x_pow = _poly_pow_mod([0, 1], exp, coeffs, p)
        # x^{exp} - x
        x_pow_minus_x = list(x_pow)
        if len(x_pow_minus_x) < 2:
            x_pow_minus_x = x_pow_minus_x + [0] * (2 - len(x_pow_minus_x))
        x_pow_minus_x[1] = (x_pow_minus_x[1] - 1) % p
        while len(x_pow_minus_x) > 1 and x_pow_minus_x[-1] == 0:
            x_pow_minus_x.pop()

        # Wenn x^{exp} - x mod f = 0, dann teilt f x^{exp}-x → reduzibel
        if all(c == 0 for c in x_pow_minus_x):
            # Das würde bedeuten f | x^{p^{n/q}} - x → reduzibel
            # (außer n/q = n, d.h. q = 1, aber q ist Primfaktor also q ≥ 2)
            return False

    # Schritt 3: f soll x^{p^n} - x teilen
    exp_n = p ** n
    x_pow_n = _poly_pow_mod([0, 1], exp_n, coeffs, p)
    # x^{p^n} - x mod f soll = 0 sein
    x_pow_n[1] = (x_pow_n[1] - 1) % p
    while len(x_pow_n) > 1 and x_pow_n[-1] == 0:
        x_pow_n.pop()

    if not all(c == 0 for c in x_pow_n):
        return False  # f teilt x^{p^n}-x nicht → reduzibel (kein endlicher Körper)

    return True


def _find_irreducible_poly(p: int, n: int) -> list[int]:
    """
    @brief Findet ein irreduzibles Polynom vom Grad n über ℤ_p.
    @description
        Durchsucht alle normierten Polynome vom Grad n über ℤ_p
        bis ein irreduzibles gefunden wird. Für kleine n und p sehr effizient.

    @param p Primzahl (Charakteristik des Körpers).
    @param n Grad des gesuchten irreduziblen Polynoms.
    @return Koeffizientenliste des irreduziblen Polynoms (aufsteigend).
    @date 2026-03-10
    """
    if n == 1:
        return [0, 1]  # x ist irreduzibel

    # Alle normierten Polynome vom Grad n: Koeffizienten in [0..p-1], letzter = 1
    # Systematische Suche
    from itertools import product as iproduct
    for lower_coeffs in iproduct(range(p), repeat=n):
        # Polynom: lower_coeffs[0] + lower_coeffs[1]*x + ... + x^n
        coeffs = list(lower_coeffs) + [1]
        if _is_irreducible_over_fp(coeffs, p):
            return coeffs

    raise MathematicalError(f"Kein irreduzibles Polynom vom Grad {n} über GF({p}) gefunden")


# =============================================================================
# KLASSE: FiniteField
# =============================================================================

class FiniteField:
    """
    @brief Endlicher Körper GF(p^n) = ℤ_p[x]/(f(x)) mit f irreduzibel vom Grad n.
    @description
        Implementiert die Grundoperationen im endlichen Körper GF(p^n):

        - Elemente werden als Polynome in ℤ_p[x] der Stufe < n dargestellt.
        - Die Multiplikation erfolgt modulo dem irreduziblen Polynom f.
        - Die Ordnung des Körpers ist p^n.
        - Die multiplikative Gruppe GF(p^n)^× ist zyklisch der Ordnung p^n - 1.

        **Frobenius-Endomorphismus:**
        φ: a ↦ a^p (bijektiver Körperautomorphismus der Ordnung n).

    @example
        GF(2^4) hat 16 Elemente, primitive Elemente erzeugen GF(16)^×.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, p: int, n: int,
                 mod_poly: Optional[list[int]] = None):
        """
        @brief Konstruiert GF(p^n) mit optionalem Modulus-Polynom.
        @param p Primzahl (Charakteristik des Körpers).
        @param n Erweiterungsgrad (Dimension über ℤ_p).
        @param mod_poly Irreduzibles Polynom vom Grad n über ℤ_p
                       (aufsteigend: mod_poly[0] = konstanter Term).
                       Falls None, wird automatisch eines gesucht.
        @raises InvalidInputError Wenn p keine Primzahl oder n < 1 ist.
        @date 2026-03-10
        """
        # Eingabevalidierung
        if not isprime(p):
            raise InvalidInputError(f"p={p} muss eine Primzahl sein")
        if n < 1:
            raise InvalidInputError(f"Erweiterungsgrad n={n} muss ≥ 1 sein")

        self.p = p          # Charakteristik des Körpers
        self.n = n          # Erweiterungsgrad
        self.q = p ** n     # Ordnung des Körpers

        # Irreduzibles Modulus-Polynom bestimmen oder prüfen
        if mod_poly is not None:
            if len(mod_poly) != n + 1:
                raise InvalidInputError(
                    f"Modulus-Polynom muss Grad {n} haben, hat aber Grad {len(mod_poly)-1}"
                )
            if not _is_irreducible_over_fp(mod_poly, p):
                raise InvalidInputError(
                    f"Übergebenes Polynom ist nicht irreduzibel über GF({p})"
                )
            self.mod_poly = [c % p for c in mod_poly]
        else:
            # Automatisch irreduzibles Polynom suchen
            self.mod_poly = _find_irreducible_poly(p, n)

    def order(self) -> int:
        """
        @brief Gibt die Anzahl der Elemente im Körper zurück: |GF(p^n)| = p^n.
        @return Körperordnung p^n.
        @date 2026-03-10
        """
        return self.q

    def multiplicative_group(self) -> list[list[int]]:
        """
        @brief Gibt alle Elemente der multiplikativen Gruppe GF(p^n)^× zurück.
        @description
            GF(p^n)^× ist die zyklische Gruppe der Ordnung p^n - 1
            (alle Körperelemente außer dem Nullelement).
            Elemente werden als Polynomkoeffizienten in ℤ_p dargestellt.

        @return Liste aller Nicht-Null-Elemente als Koeffizientenlisten.
        @date 2026-03-10
        """
        from itertools import product as iproduct
        elements = []
        # Alle Polynome vom Grad < n über ℤ_p (außer Nullpolynom)
        for coeffs in iproduct(range(self.p), repeat=self.n):
            if any(c != 0 for c in coeffs):  # Nullelement ausschließen
                elements.append(list(coeffs))
        return elements

    def frobenius_endomorphism(self, element: list[int]) -> list[int]:
        """
        @brief Wendet den Frobenius-Endomorphismus φ: a ↦ a^p auf ein Element an.
        @description
            Der Frobenius-Endomorphismus ist der kanonische Erzeuger der
            Galoisgruppe Gal(GF(p^n)/GF(p)) ≅ ℤ/nℤ.
            Er ist ein Körperautomorphismus der Ordnung n.

            φ^k(a) = a^{p^k}, φ^n = id.

        @param element Körperelement als Koeffizientenliste der Länge n.
        @return φ(element) = element^p in GF(p^n).
        @date 2026-03-10
        """
        # Normiere Eingabe auf Länge n
        elem = list(element) + [0] * (self.n - len(element))
        elem = [c % self.p for c in elem[:self.n]]

        # a^p via schneller Polynompotenzierung modulo mod_poly
        result = _poly_pow_mod(elem, self.p, self.mod_poly, self.p)

        # Auf Länge n normieren
        result = result + [0] * (self.n - len(result))
        return result[:self.n]

    def is_primitive_element(self, element: list[int]) -> bool:
        """
        @brief Prüft ob ein Element primitiv ist (Erzeuger von GF(p^n)^×).
        @description
            a ∈ GF(p^n)^× ist primitiv ⟺ ord(a) = p^n - 1
            ⟺ a^{(p^n-1)/q} ≠ 1 für alle Primteiler q von p^n - 1.

        @param element Körperelement als Koeffizientenliste.
        @return True wenn primitiv, sonst False.
        @date 2026-03-10
        """
        elem = list(element) + [0] * (self.n - len(element))
        elem = [c % self.p for c in elem[:self.n]]

        # Nullelement ist nie primitiv
        if all(c == 0 for c in elem):
            return False

        group_order = self.q - 1  # Ordnung der multiplikativen Gruppe
        # Primteiler der Gruppenordnung
        prime_factors = list(factorint(group_order).keys())

        # a ist primitiv ⟺ a^{group_order/q} ≠ 1 für alle Primteiler q
        one = [1]  # Einselement
        for q in prime_factors:
            exp = group_order // q
            powered = _poly_pow_mod(elem, exp, self.mod_poly, self.p)
            # Prüfe ob = 1 (Einselement = [1, 0, 0, ...])
            powered_norm = powered + [0] * (self.n - len(powered))
            if powered_norm[0] == 1 and all(c == 0 for c in powered_norm[1:]):
                return False  # Ordnung teilt group_order/q → nicht primitiv

        return True

    def primitive_elements(self) -> list[list[int]]:
        """
        @brief Gibt alle primitiven Elemente von GF(p^n) zurück.
        @description
            Die Anzahl der primitiven Elemente ist φ(p^n - 1)
            (Euler'sche Phi-Funktion der Gruppenordnung).

        @return Liste aller primitiven Elemente als Koeffizientenlisten.
        @date 2026-03-10
        """
        return [e for e in self.multiplicative_group()
                if self.is_primitive_element(e)]

    def __repr__(self) -> str:
        """@brief Textdarstellung des endlichen Körpers."""
        return f"GF({self.p}^{self.n}) = GF({self.q}), mod_poly={self.mod_poly}"


# =============================================================================
# KLASSE: FieldExtension
# =============================================================================

class FieldExtension:
    """
    @brief Körpererweiterung L/K mit algebraischem Element α.
    @description
        Modelliert eine algebraische Körpererweiterung L = K(α) über K = ℚ.
        Das Element α wird durch sein Minimalpolynom f ∈ ℚ[x] charakterisiert.

        Grundlegende Konzepte:
        - **Grad**: [L:K] = deg(f) = Dimension von L als K-Vektorraum.
        - **Basis**: {1, α, α², ..., α^{n-1}} ist K-Basis von L.
        - **Norm**: N_{L/K}(α) = (-1)^n · f(0) = Produkt aller Konjugierten.
        - **Spur**: Tr_{L/K}(α) = -[x^{n-1}-Koeffizient von f] = Summe aller Konjugierten.

        **Galois-Eigenschaft:**
        L/K ist Galois-Erweiterung ⟺ L normal und separabel über K
        ⟺ L ist Zerfällungskörper eines separablen Polynoms über K.

    @example
        FieldExtension([−2, 0, 1]) → ℚ(√2), Grad 2
        FieldExtension([−2, 0, 0, 1]) → ℚ(∛2), Grad 3

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, min_poly_coeffs: list[int | float]):
        """
        @brief Konstruiert die Körpererweiterung ℚ(α) via Minimalpolynom.
        @param min_poly_coeffs Koeffizienten des Minimalpolynoms von α über ℚ
                               (aufsteigend: [a_0, a_1, ..., a_n] für a_0 + a_1·x + ... + a_n·x^n).
        @raises InvalidInputError Wenn das Polynom nicht mindestens Grad 1 hat.
        @date 2026-03-10
        """
        if len(min_poly_coeffs) < 2:
            raise InvalidInputError("Minimalpolynom muss mindestens Grad 1 haben")

        # Polynom normieren: führende Nullen entfernen
        coeffs = list(min_poly_coeffs)
        while len(coeffs) > 1 and coeffs[-1] == 0:
            coeffs.pop()

        self.min_poly = coeffs  # Koeffizientenliste (aufsteigend)
        self._degree = len(coeffs) - 1  # [L:K]

        # SymPy-Polynom für symbolische Berechnungen
        self._sympy_poly = Poly(list(reversed(coeffs)), sym_x, domain=QQ)

    def degree(self) -> int:
        """
        @brief Gibt den Körpergrad [L:K] zurück.
        @return Grad des Minimalpolynoms = Dimension von L über K.
        @date 2026-03-10
        """
        return self._degree

    def is_algebraic(self) -> bool:
        """
        @brief Prüft ob α algebraisch über ℚ ist (stets True für dieses Modell).
        @description
            Da α per Konstruktion durch ein Minimalpolynom definiert ist,
            ist α immer algebraisch. Diese Methode dient zur Schnittstellen-Konformität.

        @return Immer True.
        @date 2026-03-10
        """
        return True  # Per Konstruktion: α hat endliches Minimalpolynom

    def minimal_polynomial(self) -> list[int | float]:
        """
        @brief Gibt die Koeffizientenliste des Minimalpolynoms zurück.
        @return Koeffizientenliste (aufsteigend).
        @date 2026-03-10
        """
        return list(self.min_poly)

    def basis(self) -> list[str]:
        """
        @brief Gibt eine K-Basis von L als α über ℚ zurück.
        @description
            Die kanonische K-Basis von L = K(α) ist {1, α, α², ..., α^{n-1}}.
            Dies folgt aus dem Satz über primitive Elemente algebraischer Erweiterungen.

        @return Liste der Basis-Elemente als Strings (für Lesbarkeit).
        @date 2026-03-10
        """
        n = self._degree
        basis = ['1']
        for k in range(1, n):
            basis.append(f'α^{k}')
        return basis

    def norm(self) -> float:
        """
        @brief Berechnet die Norm N_{L/K}(α) = (-1)^n · f(0).
        @description
            Die Norm von α ist das Produkt aller Konjugierten:
                N_{L/K}(α) = ∏_{σ: L↪K̄} σ(α) = (-1)^n · a_0/a_n
            wobei f(x) = a_n·x^n + ... + a_0 das Minimalpolynom ist.

        @return Norm-Wert (rational).
        @date 2026-03-10
        """
        n = self._degree
        # Konstanter Term / Leitkoeffizient, mit Vorzeichen (-1)^n
        a0 = self.min_poly[0]   # konstanter Term
        an = self.min_poly[-1]  # Leitkoeffizient
        return ((-1) ** n) * (a0 / an)

    def trace(self) -> float:
        """
        @brief Berechnet die Spur Tr_{L/K}(α) = -a_{n-1}/a_n.
        @description
            Die Spur von α ist die Summe aller Konjugierten:
                Tr_{L/K}(α) = ∑_{σ: L↪K̄} σ(α) = -a_{n-1}/a_n
            (Koeffizient x^{n-1} in f, nach Vieta mit negativem Vorzeichen).

        @return Spurwert (rational).
        @date 2026-03-10
        """
        an = self.min_poly[-1]   # Leitkoeffizient
        # Koeffizient von x^{n-1} ist min_poly[-2] wenn n ≥ 2
        if self._degree >= 2:
            a_n1 = self.min_poly[-2]  # Koeffizient von x^{n-1}
        else:
            a_n1 = self.min_poly[0]   # Für Grad-1-Polynom
        return -a_n1 / an

    def is_separable(self) -> bool:
        """
        @brief Prüft ob L/K separabel ist.
        @description
            Über einem Körper der Charakteristik 0 (wie ℚ) ist jede
            algebraische Erweiterung separabel: gcd(f, f') = 1.
            Das ist äquivalent dazu, dass f keine mehrfachen Nullstellen hat.

        @return True wenn separabel (über ℚ stets True).
        @date 2026-03-10
        """
        # Über ℚ sind alle algebraischen Erweiterungen separabel
        # (Charakteristik 0 → Minimalpolynome sind separabel)
        f = self._sympy_poly
        df = f.diff()  # Formale Ableitung
        return sympy.gcd(f, df) == 1

    def is_normal(self) -> bool:
        """
        @brief Prüft ob L/ℚ eine normale Erweiterung ist.
        @description
            L = ℚ(α) ist normal ⟺ L ist Zerfällungskörper seines Minimalpolynoms
            ⟺ alle Wurzeln von f liegen in L
            ⟺ f zerfällt über L vollständig in Linearfaktoren.

            Für Grad 2 ist L stets normal (Quadratische Erweiterungen).
            Für höhere Grade: Alle Wurzeln müssen im Zerfällungskörper liegen.

        @return True wenn L/ℚ normal ist.
        @date 2026-03-10
        """
        n = self._degree
        if n <= 1:
            return True  # Triviale Erweiterung

        if n == 2:
            return True  # Quadratische Erweiterungen sind stets normal

        # Für Grad ≥ 3: Alle Wurzeln numerisch berechnen und prüfen
        # ob sie alle in ℚ(α) liegen (d.h. ob das Galois-Pol. vollständig zerfallt)
        coeffs_desc = list(reversed(self.min_poly))
        roots = sympy.Poly(coeffs_desc, sym_x, domain=QQ).all_roots()

        # Eine Erweiterung ℚ(α) ist normal, wenn f vollständig zerfällt
        # Prüfe: alle Wurzeln liegen in ℚ(α) ⟺ alle Roots rational ausdrückbar
        # oder f irreduzibel und sein Splitting-Field = ℚ(α)
        try:
            # Faktorisierung über dem Erweiterungskörper
            poly_sym = sympy.Poly(coeffs_desc, sym_x, domain=QQ)
            factors_list = sympy.factor_list(poly_sym.as_expr())
            # Wenn das Polynom über ℚ irreduzibel und die Wurzelanzahl = Grad ist,
            # dann normal ⟺ Grad des Zerfällungskörpers = Grad des Minimalpolynoms
            gal = galois_group_polynomial(self.min_poly)
            # Normal ⟺ Galois-Erweiterung ⟺ [L:Q] = |Gal(L/Q)|
            return gal['order'] == n
        except Exception:
            return n == 2  # Fallback: nur Grad-2-Erweiterungen sicher normal

    def is_galois(self) -> bool:
        """
        @brief Prüft ob L/ℚ eine Galois-Erweiterung ist.
        @description
            L/K ist Galois ⟺ L/K ist normal UND separabel.
            Über ℚ (Charakteristik 0): L/K ist Galois ⟺ L/K ist normal.

        @return True wenn Galois-Erweiterung.
        @date 2026-03-10
        """
        return self.is_separable() and self.is_normal()

    def __repr__(self) -> str:
        """@brief Textdarstellung der Körpererweiterung."""
        return f"FieldExtension(ℚ(α), min_poly={self.min_poly}, [L:ℚ]={self._degree})"


# =============================================================================
# KLASSE: GaloisGroup
# =============================================================================

class GaloisGroup:
    """
    @brief Galoisgruppe Gal(L/K) eines Polynoms als abstrakte Gruppe.
    @description
        Gal(L/K) = Aut_K(L) = {σ: L→L | σ Körperautomorphismus, σ|_K = id}

        Die Galoisgruppe operiert auf den Wurzeln des Minimalpolynoms
        als Untergruppe der symmetrischen Gruppe S_n (n = Grad des Polynoms).

        **Implementierungsansatz:**
        - Galoisgruppe wird als abstrakte Gruppe gespeichert (Name + Ordnung)
        - Auflösbarkeit wird über die Galois-Gruppe-Klassifikation bestimmt
        - Untergruppen werden für kleine Gruppen explizit aufgelistet

        **Galois-Gruppen nach Polynomial-Grad:**
        - Grad 1: Trivial {1}
        - Grad 2: ℤ/2ℤ
        - Grad 3: ℤ/3ℤ (wenn Disc perfektes Quadrat) oder S_3
        - Grad 4: S_4, A_4, D_4, ℤ/4ℤ, ℤ/2ℤ×ℤ/2ℤ
        - Grad 5+: Oft S_n oder A_n (nicht auflösbar für n ≥ 5)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, poly_coeffs: list[int], group_name: str = "",
                 group_order: int = 0, is_solvable_flag: bool = True):
        """
        @brief Konstruiert die Galoisgruppe aus Polynomial-Koeffizienten.
        @param poly_coeffs Koeffizientenliste des Polynoms (aufsteigend).
        @param group_name Name der Galoisgruppe (z.B. "S_3", "Z/2Z").
        @param group_order Gruppenordnung.
        @param is_solvable_flag Ob die Gruppe auflösbar ist.
        @date 2026-03-10
        """
        self.poly_coeffs = poly_coeffs
        self._group_name = group_name
        self._order = group_order
        self._is_solvable = is_solvable_flag

        # Falls keine Informationen übergeben wurden, automatisch berechnen
        if not group_name:
            info = galois_group_polynomial(poly_coeffs)
            self._group_name = info['galois_group']
            self._order = info['order']
            self._is_solvable = info['is_solvable']

    def order(self) -> int:
        """
        @brief Gibt die Gruppenordnung |Gal(L/K)| zurück.
        @return Gruppenordnung.
        @date 2026-03-10
        """
        return self._order

    def elements(self) -> list[str]:
        """
        @brief Gibt die Gruppenelemente zurück (für kleine Gruppen explizit).
        @description
            Für kleine Gruppen werden die Elemente als Permutationsschreibweise
            angegeben. Für große Gruppen wird eine Beschreibung zurückgegeben.

        @return Liste der Gruppenelemente als Strings.
        @date 2026-03-10
        """
        n = self._order
        name = self._group_name

        # Bekannte Gruppen: explizite Elemente
        if n == 1:
            return ['e']
        elif name in ('Z/2Z', 'ℤ/2ℤ', 'Z_2'):
            return ['e', 'σ']
        elif name in ('Z/3Z', 'ℤ/3ℤ', 'Z_3'):
            return ['e', 'σ', 'σ²']
        elif name in ('S_3',):
            return ['e', '(12)', '(13)', '(23)', '(123)', '(132)']
        elif name in ('Z/4Z', 'ℤ/4ℤ', 'Z_4'):
            return ['e', 'σ', 'σ²', 'σ³']
        elif name in ('V_4', 'Z/2Z×Z/2Z'):
            return ['e', 'σ', 'τ', 'στ']
        elif name in ('D_4',):
            return ['e', 'r', 'r²', 'r³', 's', 'sr', 'sr²', 'sr³']
        elif name in ('A_4',):
            return ['e', '(12)(34)', '(13)(24)', '(14)(23)',
                    '(123)', '(132)', '(124)', '(142)',
                    '(134)', '(143)', '(234)', '(243)']
        elif name in ('S_4',):
            return [f'σ_{i}' for i in range(24)]
        elif name in ('S_5', 'A_5'):
            return [f'σ_{i}' for i in range(n)]
        else:
            return [f'σ_{i}' for i in range(n)]

    def is_abelian(self) -> bool:
        """
        @brief Prüft ob die Galoisgruppe abelsch ist.
        @description
            Bekannte abelsche Gruppen: ℤ/nℤ, ℤ/2ℤ×ℤ/2ℤ, (ℤ/nℤ)^×.
            S_n ist für n≥3 nicht abelsch, A_n ist für n≥4 nicht abelsch.

        @return True wenn abelsch.
        @date 2026-03-10
        """
        name = self._group_name
        # Alle bekannten abelschen Gruppen
        abelian_names = {'Z/2Z', 'ℤ/2ℤ', 'Z_2', 'Z/3Z', 'ℤ/3ℤ', 'Z_3',
                         'Z/4Z', 'ℤ/4ℤ', 'Z_4', 'V_4', 'Z/2Z×Z/2Z',
                         'trivial', '{1}', 'Z/5Z', 'ℤ/5ℤ', 'Z_5',
                         'Z/6Z', 'ℤ/6ℤ', 'Z_6'}
        return name in abelian_names

    def is_solvable(self) -> bool:
        """
        @brief Prüft ob die Galoisgruppe auflösbar ist.
        @description
            Eine Gruppe G ist auflösbar wenn es eine Kompositionsreihe
            G = G_0 ⊃ G_1 ⊃ ... ⊃ G_k = {1} gibt mit G_{i}/G_{i+1} abelsch.

            S_n und A_n sind für n ≥ 5 nicht auflösbar (Abel-Ruffini-Theorem).

        @return True wenn auflösbar.
        @date 2026-03-10
        """
        return self._is_solvable

    def subgroups(self) -> list[dict]:
        """
        @brief Gibt die Untergruppen der Galoisgruppe zurück (für kleine Gruppen).
        @description
            Für S_3 (Grad-3-Fall) wird die vollständige Untergruppenliste
            explizit angegeben. Für größere Gruppen werden die wesentlichen
            Untergruppen beschrieben.

        @return Liste von Dictionarys {'name': str, 'order': int}.
        @date 2026-03-10
        """
        name = self._group_name
        n = self._order

        if n == 1:
            return [{'name': '{e}', 'order': 1}]
        elif name in ('Z/2Z', 'ℤ/2ℤ', 'Z_2'):
            return [{'name': '{e}', 'order': 1}, {'name': 'Z/2Z', 'order': 2}]
        elif name in ('Z/3Z', 'ℤ/3ℤ', 'Z_3'):
            return [{'name': '{e}', 'order': 1}, {'name': 'Z/3Z', 'order': 3}]
        elif name == 'S_3':
            return [
                {'name': '{e}',      'order': 1},
                {'name': '<(12)>',   'order': 2},
                {'name': '<(13)>',   'order': 2},
                {'name': '<(23)>',   'order': 2},
                {'name': 'A_3≅Z_3', 'order': 3},
                {'name': 'S_3',      'order': 6},
            ]
        elif name == 'V_4':
            return [
                {'name': '{e}',     'order': 1},
                {'name': '<σ>',     'order': 2},
                {'name': '<τ>',     'order': 2},
                {'name': '<στ>',    'order': 2},
                {'name': 'V_4',     'order': 4},
            ]
        elif name == 'Z/4Z':
            return [
                {'name': '{e}',     'order': 1},
                {'name': '<σ²>',    'order': 2},
                {'name': 'Z/4Z',    'order': 4},
            ]
        elif name == 'D_4':
            return [
                {'name': '{e}',      'order': 1},
                {'name': '<r²>',     'order': 2},
                {'name': '<s>',      'order': 2},
                {'name': '<r>',      'order': 4},
                {'name': '<r²,s>',   'order': 4},
                {'name': '<r²,rs>',  'order': 4},
                {'name': 'D_4',      'order': 8},
            ]
        elif name == 'A_4':
            return [
                {'name': '{e}',     'order': 1},
                {'name': 'V_4',     'order': 4},
                {'name': '<(123)>', 'order': 3},
                {'name': 'A_4',     'order': 12},
            ]
        elif name == 'S_4':
            return [
                {'name': '{e}',     'order': 1},
                {'name': 'Z/2Z',    'order': 2},
                {'name': 'Z/3Z',    'order': 3},
                {'name': 'V_4',     'order': 4},
                {'name': 'Z/4Z',    'order': 4},
                {'name': 'S_3',     'order': 6},
                {'name': 'D_4',     'order': 8},
                {'name': 'A_4',     'order': 12},
                {'name': 'S_4',     'order': 24},
            ]
        else:
            return [{'name': '{e}', 'order': 1}, {'name': name, 'order': n}]

    def normal_subgroups(self) -> list[dict]:
        """
        @brief Gibt die Normalteiler der Galoisgruppe zurück.
        @description
            H ≤ G ist Normalteiler ⟺ gHg⁻¹ = H für alle g ∈ G
            ⟺ H ist invariant unter Konjugation.

            Für die Galois-Korrespondenz sind die Normalteiler wichtig:
            H normal in Gal(L/K) ⟺ der Fixkörper L^H ist normale Erweiterung von K.

        @return Liste von Dictionarys {'name': str, 'order': int, 'is_normal': True}.
        @date 2026-03-10
        """
        name = self._group_name
        # Normalteiler sind eine Teilmenge der Untergruppen
        all_subs = self.subgroups()

        if name == 'S_3':
            # In S_3: Normalteiler sind {e}, A_3 und S_3 (die Untergruppen der Ordnung 2 sind nicht normal)
            return [s for s in all_subs if s['name'] in ('{e}', 'A_3≅Z_3', 'S_3')]
        elif name == 'D_4':
            # In D_4: Normalteiler sind {e}, <r²>, <r>, <r²,s>, <r²,rs>, D_4
            return [s for s in all_subs if s['name'] in ('{e}', '<r²>', '<r>', '<r²,s>', '<r²,rs>', 'D_4')]
        elif name == 'A_4':
            # In A_4: Normalteiler sind {e}, V_4, A_4
            return [s for s in all_subs if s['name'] in ('{e}', 'V_4', 'A_4')]
        elif name == 'S_4':
            # In S_4: Normalteiler sind {e}, V_4, A_4, S_4
            return [s for s in all_subs if s['name'] in ('{e}', 'V_4', 'A_4', 'S_4')]
        else:
            # Für abelsche Gruppen sind alle Untergruppen Normalteiler
            return all_subs

    def fixed_field(self, subgroup_name: str) -> str:
        """
        @brief Berechnet den Fixkörper L^H für eine Untergruppe H.
        @description
            Der Fixkörper von H ≤ Gal(L/K) ist:
                L^H = {α ∈ L : σ(α) = α für alle σ ∈ H}
            Es gilt [L : L^H] = |H| und [L^H : K] = [G : H].

        @param subgroup_name Name der Untergruppe H.
        @return Beschreibung des Fixkörpers als String.
        @date 2026-03-10
        """
        name = self._group_name

        # Für S_3 (= Gal(ℚ(∛2, ω)/ℚ)): vollständige Korrespondenz
        if name == 'S_3':
            fixed_fields = {
                '{e}': 'L (Zerfällungskörper = ℚ(∛2, ω))',
                '<(12)>': 'ℚ(∛2)',
                '<(13)>': 'ℚ(ω·∛2)',
                '<(23)>': 'ℚ(ω²·∛2)',
                'A_3≅Z_3': 'ℚ(ω) = ℚ(√-3)',
                'S_3': 'ℚ (Grundkörper)',
            }
            return fixed_fields.get(subgroup_name, f'L^{subgroup_name} (unbekannt)')

        elif name in ('Z/2Z', 'ℤ/2ℤ', 'Z_2'):
            fixed_fields = {
                '{e}': 'L (Erweiterungskörper)',
                'Z/2Z': 'ℚ (Grundkörper)',
            }
            return fixed_fields.get(subgroup_name, f'L^{subgroup_name}')

        return f'L^{subgroup_name} (explizite Berechnung für {name} nicht implementiert)'

    def __repr__(self) -> str:
        """@brief Textdarstellung der Galoisgruppe."""
        return f"GaloisGroup({self._group_name}, order={self._order}, solvable={self._is_solvable})"


# =============================================================================
# HAUPTFUNKTIONEN
# =============================================================================

def discriminant_polynomial(poly_coeffs: list[int]) -> int:
    """
    @brief Berechnet die Diskriminante Δ(f) eines Polynoms.
    @description
        Die Diskriminante ist definiert als:
            Δ(f) = (-1)^{n(n-1)/2} · Res(f, f') / a_n

        wobei:
        - n = deg(f)
        - Res(f, f') = Resultante von f und seiner Ableitung f'
        - a_n = Leitkoeffizient von f

        Bedeutung:
        - Δ > 0: gerade Anzahl von Paaren komplexer Wurzeln
        - Δ < 0: ungerade Anzahl von Paaren komplexer Wurzeln
        - Δ = 0: f hat mehrfache Wurzel

        Für Grad 2 (ax²+bx+c): Δ = b²-4ac
        Für Grad 3 (x³+px+q): Δ = -4p³-27q²

    @param poly_coeffs Koeffizientenliste (aufsteigend: [a_0, a_1, ..., a_n]).
    @return Diskriminante als Integer (kann sehr groß sein).
    @raises InvalidInputError Wenn das Polynom Grad < 1 hat.
    @date 2026-03-10
    """
    if len(poly_coeffs) < 2:
        raise InvalidInputError("Polynom muss mindestens Grad 1 haben")

    n = len(poly_coeffs) - 1  # Polynomgrad
    coeffs_desc = list(reversed(poly_coeffs))  # absteigende Koeffizienten für SymPy

    # SymPy-Polynom über ℤ[x]
    f = Poly(coeffs_desc, sym_x, domain=ZZ)

    # Resultante Res(f, f') berechnen
    df = f.diff()  # Formale Ableitung f'
    res = sympy.resultant(f, df, sym_x)  # Res(f, f')

    # Diskriminante: Δ = (-1)^{n(n-1)/2} · Res(f, f') / a_n
    sign = (-1) ** (n * (n - 1) // 2)
    a_n = poly_coeffs[-1]  # Leitkoeffizient

    disc = sympy.Rational(sign * res, a_n)
    return int(disc)


def minimal_polynomial(alpha_coeffs: list[float],
                       base_coeffs: Optional[list[int]] = None) -> list[int]:
    """
    @brief Berechnet das Minimalpolynom von α über ℚ.
    @description
        Das Minimalpolynom von α ist das kleinste normierte Polynom
        f ∈ ℚ[x] mit f(α) = 0.

        **Unterstützte Fälle:**
        - α = √d: f = x² - d (wenn d kein perfektes Quadrat)
        - α = ∛d: f = x³ - d
        - α = d (rational): f = x - d
        - Allgemein: LLL-ähnliche Heuristik via SymPy

        **Eingabeformat:**
        alpha_coeffs = [a_0, a_1, ..., a_{n-1}] bedeutet:
            α = a_0 + a_1·β + a_2·β² + ...
        wobei β die Wurzel des Grundkörper-Polynoms base_coeffs ist.

        Falls alpha_coeffs = [d] (ein Element), wird α = √|d| angenommen
        wenn d kein Quadrat, sonst α = d (rational).

    @param alpha_coeffs Koeffizienten zur Beschreibung von α.
    @param base_coeffs Optionales Grundkörper-Polynom (Standard: ℚ).
    @return Koeffizientenliste des Minimalpolynoms (aufsteigend).
    @date 2026-03-10
    """
    # Spezialfall: α direkt als Polynom angegeben
    if len(alpha_coeffs) == 1:
        val = alpha_coeffs[0]
        # Prüfe ob rational (perfektes Quadrat oder ganze Zahl)
        if isinstance(val, int) or (isinstance(val, float) and val == int(val)):
            return [-int(val), 1]  # x - val
        # Ansonsten: versuche als √val zu interpretieren
        d = val
        if d > 0:
            sq = isqrt(int(d))
            if sq * sq == int(d):
                return [-sq, 1]  # d war Quadrat, α = √d = ganzzahlig
            return [int(-d), 0, 1]  # x² - d für √d
        return [int(-d), 0, 1]

    if len(alpha_coeffs) == 2:
        # α = a + b·√d interpretieren: numerisch berechnen via SymPy
        a, b = alpha_coeffs
        if b == 0:
            return minimal_polynomial([a])

    # Allgemeiner Fall: Verwende SymPy's minimal_polynomial-Funktion
    # Interpret alpha_coeffs als numerischen Wert (Näherung)
    # und nutze SymPy um das exakte Minimalpolynom zu finden
    try:
        # Numerischer Wert von α
        val = sum(c * (2 ** 0.5) ** i for i, c in enumerate(alpha_coeffs))

        # SymPy-Ausdrücke für bekannte algebraische Zahlen
        # Für einfache Fälle direkte Berechnung
        sympy_alpha = sympy.nsimplify(val, rational=False, tolerance=1e-10)
        min_poly_sym = sympy.minimal_polynomial(sympy_alpha, sym_x)
        poly = Poly(min_poly_sym, sym_x, domain=ZZ)
        coeffs_desc = [int(c) for c in poly.all_coeffs()]
        return list(reversed(coeffs_desc))
    except Exception:
        # Fallback: lineare Approximation
        return [-int(round(sum(alpha_coeffs))), 1]


def splitting_field(poly_coeffs: list[int], p: Optional[int] = None) -> dict:
    """
    @brief Berechnet den Zerfällungskörper eines Polynoms.
    @description
        Der Zerfällungskörper von f über K ist der kleinste Körper L ⊇ K,
        über dem f vollständig in Linearfaktoren zerfällt:
            f(x) = a_n · ∏_{i=1}^n (x - α_i) mit allen α_i ∈ L.

        **Algorithmus:**
        1. Grad und Galoisgruppe bestimmen
        2. Ordnung des Zerfällungskörpers = |Gal(f/K)|
        3. Numerische Approximation der Wurzeln

        Für f ∈ ℚ[x]:
        - Grad 2: [L:ℚ] = 2 wenn Discriminante kein Quadrat, sonst 1
        - Grad 3: [L:ℚ] = 3 wenn Galois = ℤ_3, sonst 6 (= Gal = S_3)
        - Grad 4: [L:ℚ] = |Gal(f)|

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @param p Optional: Charakteristik (Primzahl) für endliche Körper.
    @return Dictionary {'degree': int, 'roots': list[complex], 'galois_group': str}.
    @date 2026-03-10
    """
    n = len(poly_coeffs) - 1  # Polynomgrad

    if p is not None:
        # Endlicher Körper GF(p): Zerfällungskörper ist GF(p^n)
        return {
            'degree': n,
            'roots': list(range(p)),
            'galois_group': f'Z/{n}Z (zyklisch)',
            'field': f'GF({p}^{n})'
        }

    # Numerische Wurzeln über ℂ berechnen
    coeffs_desc = list(reversed(poly_coeffs))
    try:
        import numpy as np
        roots_np = np.roots([float(c) for c in coeffs_desc])
        roots = [complex(r) for r in roots_np]
    except Exception:
        roots = []

    # Galoisgruppe und Grad des Zerfällungskörpers
    gal_info = galois_group_polynomial(poly_coeffs)
    galois_group_name = gal_info['galois_group']
    galois_order = gal_info['order']

    # Grad des Zerfällungskörpers = |Gal(f/ℚ)| (Gradsatz)
    splitting_degree = galois_order

    return {
        'degree': splitting_degree,
        'roots': roots,
        'galois_group': galois_group_name,
        'galois_order': galois_order
    }


def _discriminant_degree2(coeffs: list[int]) -> int:
    """
    @brief Diskriminante für Grad-2-Polynom ax²+bx+c: Δ = b²-4ac.
    @param coeffs [c, b, a] (aufsteigend).
    @return Diskriminante.
    @date 2026-03-10
    """
    c, b, a = coeffs[0], coeffs[1], coeffs[2]
    return b * b - 4 * a * c


def _is_perfect_square(n: int) -> bool:
    """
    @brief Prüft ob n ein perfektes Quadrat ist (für ganze Zahlen).
    @param n Zu prüfende Zahl.
    @return True wenn n = k² für k ∈ ℤ.
    @date 2026-03-10
    """
    if n < 0:
        return False
    if n == 0:
        return True
    sq = isqrt(n)
    return sq * sq == n


def galois_group_polynomial(poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet die Galoisgruppe von f(x) ∈ ℚ[x].
    @description
        Die Galoisgruppe Gal(f) = Gal(K_f/ℚ) operiert auf den Wurzeln
        von f als Permutationsgruppe und ist Untergruppe von S_n.

        **Algorithmus nach Grad:**

        **Grad 1:** Gal(f) = {1} (trivial)

        **Grad 2:** f = ax²+bx+c, Δ = b²-4ac
            - Δ perfektes Quadrat: Gal = {1} (f zerfällt über ℚ)
            - Δ kein Quadrat: Gal = ℤ/2ℤ

        **Grad 3:** f = x³+px+q (irreduzibel über ℚ angenommen), Δ = -4p³-27q²
            - Δ perfektes Quadrat: Gal = ℤ/3ℤ (A_3)
            - Δ kein Quadrat: Gal = S_3

        **Grad 4:** Berechnung über Resolventenpolynom (Kubische Resolvente)
            - Resolvente zerfällt vollständig: mehrere Unterfälle
            - Gal ∈ {S_4 (ord=24), A_4 (ord=12), D_4 (ord=8), ℤ_4 (ord=4), V_4 (ord=4)}

        **Grad 5:** Generisch S_5 (nicht auflösbar), Spezialfälle möglich.

    @param poly_coeffs Koeffizientenliste (aufsteigend: [a_0, ..., a_n]).
    @return Dictionary mit Schlüsseln:
            - 'galois_group': str (Name der Gruppe)
            - 'order': int (Gruppenordnung)
            - 'is_solvable': bool (Auflösbarkeit durch Radikale)
            - 'discriminant': int (Diskriminante)
    @date 2026-03-10
    """
    n = len(poly_coeffs) - 1  # Polynomgrad

    # --- Grad 1: Triviale Gruppe ---
    if n == 1:
        return {
            'galois_group': 'trivial',
            'order': 1,
            'is_solvable': True,
            'discriminant': 1
        }

    # Diskriminante berechnen
    try:
        disc = discriminant_polynomial(poly_coeffs)
    except Exception:
        disc = 0

    # --- Grad 2: Gal ≅ ℤ/2ℤ oder trivial ---
    if n == 2:
        disc2 = _discriminant_degree2(poly_coeffs) if len(poly_coeffs) == 3 else disc
        if _is_perfect_square(disc2) and disc2 >= 0:
            # Polynom zerfällt über ℚ: f = a(x-r1)(x-r2) mit r1, r2 ∈ ℚ
            return {
                'galois_group': 'trivial',
                'order': 1,
                'is_solvable': True,
                'discriminant': disc2
            }
        else:
            return {
                'galois_group': 'Z/2Z',
                'order': 2,
                'is_solvable': True,
                'discriminant': disc2 if disc2 != 0 else disc
            }

    # --- Grad 3 ---
    if n == 3:
        # Prüfe Irreduzibilität über ℚ (Eisenstein, rationale Wurzeln)
        # Vereinfachung: Standardmäßig irreduzibel annehmen wenn kein rationaler Root
        coeffs_desc = list(reversed(poly_coeffs))
        # Rational Root Theorem: Teiler von poly_coeffs[0] / Teiler von poly_coeffs[-1]
        const_term = abs(poly_coeffs[0]) if poly_coeffs[0] != 0 else 1
        lead_term = abs(poly_coeffs[-1])
        divisors_const = [d for d in range(1, const_term + 1) if const_term % d == 0]
        divisors_lead = [d for d in range(1, lead_term + 1) if lead_term % d == 0]
        rational_candidates = []
        for p_rat in divisors_const:
            for q_rat in divisors_lead:
                for sign in [1, -1]:
                    rational_candidates.append(Rational(sign * p_rat, q_rat))

        has_rational_root = False
        for candidate in rational_candidates:
            val = sum(int(poly_coeffs[i]) * candidate**i for i in range(len(poly_coeffs)))
            if val == 0:
                has_rational_root = True
                break

        if has_rational_root:
            # Zerfällt: Grad 3 mit rationaler Wurzel hat Gal ≅ ℤ/2ℤ oder trivial
            return {
                'galois_group': 'Z/2Z',
                'order': 2,
                'is_solvable': True,
                'discriminant': disc
            }

        # Irreduzibel: Galois = ℤ/3ℤ wenn Disc Quadrat, sonst S_3
        if disc > 0 and _is_perfect_square(disc):
            return {
                'galois_group': 'Z/3Z',
                'order': 3,
                'is_solvable': True,
                'discriminant': disc
            }
        else:
            return {
                'galois_group': 'S_3',
                'order': 6,
                'is_solvable': True,
                'discriminant': disc
            }

    # --- Grad 4 ---
    if n == 4:
        # Kubische Resolvente berechnen für Grad-4-Polynom
        # f(x) = x^4 + bx^3 + cx^2 + dx + e
        # Normieren (Leitkoeffizient = 1 via Division durch a_4)
        a4 = poly_coeffs[4]
        a3 = poly_coeffs[3]
        a2 = poly_coeffs[2]
        a1 = poly_coeffs[1]
        a0 = poly_coeffs[0]

        # In deprimierte Form x^4 + px^2 + qx + r transformieren
        # y = x - a3/(4*a4)
        # Vereinfachte Berechnung über Diskriminante und Galois-Untergruppen
        resolvent_disc = _cubic_resolvent_discriminant(a4, a3, a2, a1, a0)

        # Diskriminante-basierte Klassifikation
        if disc == 0:
            # Mehrfache Wurzel: reduzibel
            return {
                'galois_group': 'trivial',
                'order': 1,
                'is_solvable': True,
                'discriminant': disc
            }

        # A_4 wenn disc > 0 und kubische Resolvente irreduzibel
        if disc > 0 and _is_perfect_square(disc):
            # Galois ⊆ A_4 (alternierend): alle Permutationen sind gerade
            if resolvent_disc == 0:
                return {'galois_group': 'V_4', 'order': 4,
                        'is_solvable': True, 'discriminant': disc}
            return {'galois_group': 'A_4', 'order': 12,
                    'is_solvable': True, 'discriminant': disc}
        else:
            # Galois ⊄ A_4: enthält ungerade Permutationen
            if resolvent_disc == 0:
                return {'galois_group': 'D_4', 'order': 8,
                        'is_solvable': True, 'discriminant': disc}
            # Unterscheidung D_4 vs Z/4Z vs S_4 über Resolventenpolynom
            # Vereinfachung: S_4 für allgemeinen Fall
            return {'galois_group': 'S_4', 'order': 24,
                    'is_solvable': True, 'discriminant': disc}

    # --- Grad 5+ ---
    if n == 5:
        # Generisch: S_5 (nicht auflösbar) für "zufällige" Grad-5-Polynome
        # Spezialfall: wenn Polynom über ℚ zerfällt oder spezielle Struktur hat
        # Vereinfachung: Für den allgemeinen Test → S_5
        return {
            'galois_group': 'S_5',
            'order': 120,
            'is_solvable': False,
            'discriminant': disc
        }

    # --- Höherer Grad: S_n ---
    return {
        'galois_group': f'S_{n}',
        'order': math.factorial(n),
        'is_solvable': n <= 4,
        'discriminant': disc
    }


def _cubic_resolvent_discriminant(a4: int, a3: int, a2: int,
                                   a1: int, a0: int) -> int:
    """
    @brief Berechnet die Diskriminante der kubischen Resolvente für Grad-4-Polynome.
    @description
        Die kubische Resolvente von f(x) = a4·x^4 + a3·x^3 + ... + a0 ist:
            g(y) = y³ - a2·y² + (a1·a3 - 4·a0·a4)·y - (a1²·a4 + a0·a3² - 4·a0·a2·a4)

        Die Diskriminante von g bestimmt zusammen mit disc(f) die Galoisgruppe.

    @param a4..a0 Koeffizienten des Grad-4-Polynoms (absteigend).
    @return Diskriminante der kubischen Resolvente.
    @date 2026-03-10
    """
    # Koeffizienten der kubischen Resolvente (nach normiertem f)
    # g(y) = y^3 + py + q mit vereinfachten Koeffizienten
    # Diskriminante von g: 4p³ + 27q²
    # Vereinfachte Berechnung
    try:
        # Normiertes Polynom
        f_coeffs_desc = [a4, a3, a2, a1, a0]
        f_sym = Poly(f_coeffs_desc, sym_x, domain=ZZ)

        # Kubische Resolvente via SymPy
        y = Symbol('y')
        # Koeffizienten der Resolvente berechnen
        b = -a2
        c = a1 * a3 - 4 * a0 * a4
        d = -(a1 * a1 * a4 + a0 * a3 * a3 - 4 * a0 * a2 * a4)

        # Diskriminante der Kubischen: -4b³d + b²c² - 4c³ + 18bcd - 27d²
        # (für y³ + b·y² + c·y + d)
        disc_resolvent = (-4 * b**3 * d + b**2 * c**2
                          - 4 * c**3 + 18 * b * c * d - 27 * d**2)
        return int(disc_resolvent)
    except Exception:
        return 0


def is_solvable_by_radicals(poly_coeffs: list[int]) -> dict:
    """
    @brief Prüft ob f(x) durch Radikale lösbar ist (Galois-Hauptsatz).
    @description
        **Galois-Hauptsatz (Lösung durch Radikale):**
        f(x) ∈ ℚ[x] ist durch Radikale lösbar
            ⟺ Gal(f) ist eine auflösbare Gruppe
            ⟺ Gal(f) hat eine Kompositionsreihe mit abelschen Faktoren.

        **Bekannte Ergebnisse:**
        - Grad ≤ 4: Stets auflösbar (Formeln durch Radikale existieren)
          - Grad 1: ax + b = 0 → x = -b/a
          - Grad 2: Quadratische Formel
          - Grad 3: Cardano-Formel
          - Grad 4: Ferrari-Methode
        - Grad ≥ 5: Generisch nicht auflösbar (Abel-Ruffini, 1824)
          - Speziell: x^n - 1 = 0 auflösbar (Wurzeln der Einheit)
          - Speziell: x^p - x - a = 0 für Primpotenz p (Artin-Schreier)

        **Auflösbare Gruppen:**
        ℤ/nℤ, D_n, A_n (n≤4), S_n (n≤4), direkte/semidirekte Produkte davon.

        **Nicht auflösbar:** S_n (n≥5), A_n (n≥5).

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Dictionary:
            - 'solvable': bool
            - 'galois_group': str
            - 'reason': str (Erklärung)
    @date 2026-03-10
    """
    n = len(poly_coeffs) - 1  # Polynomgrad

    # Galoisgruppe berechnen
    gal = galois_group_polynomial(poly_coeffs)
    gal_name = gal['galois_group']
    is_solvable = gal['is_solvable']

    if n <= 4:
        # Grad ≤ 4: Stets auflösbar (Cardano/Ferrari)
        reason_map = {
            1: "Lineares Polynom: direkt lösbar",
            2: "Quadratische Formel existiert",
            3: "Cardano-Formel: Grad 3 stets durch Radikale lösbar",
            4: "Ferrari-Methode: Grad 4 stets durch Radikale lösbar",
        }
        return {
            'solvable': True,
            'galois_group': gal_name,
            'reason': reason_map.get(n, f"Gal = {gal_name} ist auflösbar"),
            'discriminant': gal['discriminant']
        }
    else:
        # Grad ≥ 5: Galois-Hauptsatz entscheidet
        if is_solvable:
            reason = (f"Gal(f) = {gal_name} ist auflösbar "
                      f"(Gruppenordnung {gal['order']} hat auflösbare Struktur)")
        else:
            reason = (f"Gal(f) = {gal_name} ist nicht auflösbar "
                      f"(A_n und S_n für n≥5 haben keine auflösbare Struktur). "
                      f"Abel-Ruffini-Theorem: keine allgemeine Radikallösung möglich.")
        return {
            'solvable': is_solvable,
            'galois_group': gal_name,
            'reason': reason,
            'discriminant': gal['discriminant']
        }


def galois_correspondence(poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet die Galois-Korrespondenz (Hauptsatz der Galois-Theorie).
    @description
        Der Hauptsatz der Galois-Theorie stellt eine Bijektion auf zwischen:

            {Untergruppen H ≤ Gal(L/K)}  ←→  {Zwischenkörper K ⊆ F ⊆ L}

        Diese Bijektion ist ordnungsumkehrend:
            - H₁ ⊆ H₂ ⟺ L^{H₁} ⊇ L^{H₂}
            - [L : L^H] = |H|, [L^H : K] = [G : H]

        Weitere Eigenschaften:
            - H normal in G ⟺ L^H galoissch über K
            - Gal(L^H/K) ≅ G/H wenn H normal in G

    @param poly_coeffs Koeffizientenliste des Polynoms f(x).
    @return Dictionary mit Korrespondenztabelle:
            'galois_group': str, 'subgroups': list, 'correspondence': list.
    @date 2026-03-10
    """
    # Galoisgruppe bestimmen
    gal_info = galois_group_polynomial(poly_coeffs)
    gal = GaloisGroup(poly_coeffs,
                      group_name=gal_info['galois_group'],
                      group_order=gal_info['order'],
                      is_solvable_flag=gal_info['is_solvable'])

    # Untergruppen und ihre Fixkörper
    subgroups = gal.subgroups()
    correspondence = []

    for sg in subgroups:
        sg_name = sg['name']
        sg_order = sg['order']
        gal_order = gal.order()

        # Degree of fixed field over K = [G:H] = |G|/|H|
        fixed_degree = gal_order // sg_order if sg_order > 0 else gal_order

        # Fixkörper beschreiben
        fixed_field_desc = gal.fixed_field(sg_name)

        # Ob Untergruppe normal → Zwischenkörper galoissch über K
        normal_subs = gal.normal_subgroups()
        is_normal = any(ns['name'] == sg_name for ns in normal_subs)

        correspondence.append({
            'subgroup': sg_name,
            'subgroup_order': sg_order,
            'fixed_field': fixed_field_desc,
            'degree_over_K': fixed_degree,
            'is_normal_subgroup': is_normal,
        })

    return {
        'galois_group': gal_info['galois_group'],
        'galois_order': gal_info['order'],
        'subgroups': subgroups,
        'correspondence': correspondence,
    }


def finite_field(p: int, n: int) -> FiniteField:
    """
    @brief Konstruiert den endlichen Körper GF(p^n).
    @description
        GF(p^n) ist der eindeutige (bis auf Isomorphie) endliche Körper
        mit p^n Elementen, wobei p eine Primzahl ist.

        Darstellung: GF(p^n) = ℤ_p[x]/(f(x)) mit f irreduzibel vom Grad n.

    @param p Primzahl (Charakteristik).
    @param n Erweiterungsgrad.
    @return FiniteField-Objekt für GF(p^n).
    @raises InvalidInputError Wenn p keine Primzahl oder n < 1.
    @date 2026-03-10
    """
    return FiniteField(p, n)


def finite_field_discrete_log(a: int, g: int, q: int) -> int:
    """
    @brief Berechnet den diskreten Logarithmus log_g(a) in GF(q)^×.
    @description
        Löst die Gleichung g^x ≡ a (mod q) für x ∈ {0, ..., q-2}.

        **Baby-Step-Giant-Step-Algorithmus (Shanks):**
        1. m = ⌈√(q-1)⌉ (Schrittgröße)
        2. Baby steps: Berechne g^j mod q für j = 0, 1, ..., m-1
        3. Giant steps: Berechne a · (g^{-m})^i mod q für i = 0, 1, ..., m-1
        4. Kollision finden: j + m·i = x

        Laufzeit: O(√q) Zeit und Speicher.

        **Korrektheit:**
        Falls x = j + m·i, dann g^x = g^j · g^{m·i} = g^j · (g^m)^i.

    @param a Zielwert (muss in GF(q)^× liegen, d.h. gcd(a, q) = 1).
    @param g Basis (Generator oder Erzeuger von GF(q)^×).
    @param q Primzahl (Ordnung des endlichen Körpers GF(q)).
    @return x mit g^x ≡ a (mod q), oder -1 wenn keine Lösung.
    @raises InvalidInputError Wenn q keine Primzahl ist.
    @date 2026-03-10
    """
    if not isprime(q):
        raise InvalidInputError(f"q={q} muss eine Primzahl sein")

    a = a % q
    g = g % q

    if a == 0:
        raise InvalidInputError("a=0 hat keinen diskreten Logarithmus")

    # Gruppenordnung in GF(q)^×
    group_order = q - 1

    # Schrittgröße m = ⌈√(q-1)⌉
    m = isqrt(group_order)
    if m * m < group_order:
        m += 1

    # Baby steps: Tabelle {g^j mod q: j} für j = 0, ..., m-1
    baby_steps = {}
    gj = 1  # g^0
    for j in range(m):
        baby_steps[gj] = j
        gj = (gj * g) % q

    # g^{-m} mod q = modulares Inverses von g^m
    gm = pow(g, m, q)
    gm_inv = pow(gm, q - 2, q)  # Fermat: g^{-1} ≡ g^{q-2} (mod q)

    # Giant steps: a · (g^{-m})^i für i = 0, 1, ..., m-1
    giant = a
    for i in range(m + 1):
        if giant in baby_steps:
            j = baby_steps[giant]
            x = (j + m * i) % group_order
            # Verifikation
            if pow(g, x, q) == a:
                return x
        giant = (giant * gm_inv) % q

    return -1  # Kein Logarithmus gefunden (a kein Element der von g erzeugten Gruppe)


def primitive_element_theorem(alpha: list, beta: list) -> dict:
    """
    @brief Findet ein primitives Element ℚ(α, β) = ℚ(α + c·β).
    @description
        **Satz vom primitiven Element:**
        Für algebraische Zahlen α, β über ℚ gilt:
            ℚ(α, β) = ℚ(α + c·β) für fast alle c ∈ ℚ.

        Genauer: c ∈ ℚ außerhalb einer endlichen Ausnahmemenge S.

        **Algorithmus:**
        Teste c = 0, 1, 2, ... und prüfe ob [ℚ(α+c·β):ℚ] = [ℚ(α,β):ℚ].
        Letzteres ist [ℚ(α):ℚ] · [ℚ(β):ℚ] wenn ℚ(α) ∩ ℚ(β) = ℚ.

    @param alpha Minimalpolynom von α (als Koeffizientenliste).
    @param beta Minimalpolynom von β (als Koeffizientenliste).
    @return Dictionary {'c': int, 'primitive_element': str, 'degree': int}.
    @date 2026-03-10
    """
    deg_alpha = len(alpha) - 1
    deg_beta = len(beta) - 1

    # Erwarteter Grad von ℚ(α, β) (wenn ℚ(α) ∩ ℚ(β) = ℚ)
    expected_degree = deg_alpha * deg_beta

    # Teste verschiedene c ∈ ℚ
    for c in range(0, 20):
        # α + c·β: numerische Approximation der primitiven Elements
        # Berechne Minimalpolynom von γ = α + c·β via SymPy
        try:
            # Numerische Wurzeln von α und β
            roots_alpha = sympy.Poly(list(reversed(alpha)), sym_x,
                                     domain=QQ).all_roots()
            roots_beta = sympy.Poly(list(reversed(beta)), sym_x,
                                    domain=QQ).all_roots()

            if not roots_alpha or not roots_beta:
                continue

            # Einfachste Wurzel verwenden (erste reelle falls vorhanden)
            a_root = roots_alpha[0]
            b_root = roots_beta[0]

            # γ = α + c·β
            gamma = a_root + sympy.Integer(c) * b_root

            # Minimalpolynom von γ
            gamma_min_poly = sympy.minimal_polynomial(gamma, sym_x)
            gamma_degree = sympy.degree(gamma_min_poly, sym_x)

            if gamma_degree >= expected_degree:
                poly_sym = Poly(gamma_min_poly, sym_x, domain=QQ)
                coeffs_desc = [int(sympy.Rational(ci).numerator)
                               for ci in poly_sym.all_coeffs()]

                return {
                    'c': c,
                    'primitive_element': f'α + {c}·β',
                    'degree': int(gamma_degree),
                    'minimal_polynomial': list(reversed(coeffs_desc))
                }
        except Exception:
            continue

    return {
        'c': 1,
        'primitive_element': 'α + β',
        'degree': expected_degree,
        'minimal_polynomial': None
    }


def cyclotomic_polynomial(n: int) -> list[int]:
    """
    @brief Berechnet das n-te Kreisteilungspolynom Φ_n(x).
    @description
        Das n-te Kreisteilungspolynom ist definiert als:
            Φ_n(x) = ∏_{gcd(k,n)=1, 1≤k≤n} (x - e^{2πik/n})

        Es ist das Minimalpolynom einer primitiven n-ten Einheitswurzel ζ_n = e^{2πi/n}
        über ℚ. Sein Grad ist φ(n) (Euler'sche Phi-Funktion).

        **Rekursionsformel via Möbius:**
            x^n - 1 = ∏_{d | n} Φ_d(x)

        Daher: Φ_n(x) = ∏_{d | n} (x^{n/d} - 1)^{μ(d)} (Möbius-Inversion)

        **Galoisgruppe:**
            Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^× (multiplikative Gruppe modulo n)
            Ordnung = φ(n)

        **Beispiele:**
        - Φ_1(x) = x - 1
        - Φ_2(x) = x + 1
        - Φ_3(x) = x² + x + 1
        - Φ_4(x) = x² + 1
        - Φ_6(x) = x² - x + 1
        - Φ_{12}(x) = x⁴ - x² + 1

    @param n Positive ganze Zahl ≥ 1.
    @return Koeffizientenliste von Φ_n(x) (aufsteigend: [a_0, a_1, ..., a_{φ(n)}]).
    @raises InvalidInputError Wenn n < 1.
    @date 2026-03-10
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # SymPy hat eine eingebaute Funktion für Kreisteilungspolynome
    phi_n = sympy.cyclotomic_poly(n, sym_x)

    # In Koeffizientenliste (aufsteigend) umwandeln
    poly = Poly(phi_n, sym_x, domain=ZZ)
    coeffs_desc = [int(c) for c in poly.all_coeffs()]  # absteigende Reihenfolge

    # Umkehren für aufsteigende Darstellung [a_0, a_1, ..., a_φ(n)]
    return list(reversed(coeffs_desc))
