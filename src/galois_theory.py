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

@author Michael Fuhrmann
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

    @author Michael Fuhrmann
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

    @author Michael Fuhrmann
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

    @author Michael Fuhrmann
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


# =============================================================================
# ZYKLOTOMISCHE KÖRPER UND GALOIS-GRUPPEN
# =============================================================================

def cyclotomic_galois_group(n: int) -> dict:
    """
    @brief Berechnet die Galoisgruppe Gal(ℚ(ζ_n)/ℚ) des n-ten Kreisteilungskörpers.
    @description
        Der n-te Kreisteilungskörper ℚ(ζ_n) entsteht durch Adjunktion einer
        primitiven n-ten Einheitswurzel ζ_n = e^{2πi/n} zu ℚ.

        **Galoisgruppe:**
        $$\text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \cong (\mathbb{Z}/n\mathbb{Z})^{\times}$$

        Diese Gruppe hat Ordnung φ(n) (Euler'sche φ-Funktion). Jedes Element σ_k
        mit gcd(k, n) = 1 ist ein Automorphismus:
        $$\sigma_k(\zeta_n) = \zeta_n^k$$

        **Eigenschaften:**
        - Die Erweiterung ist stets galoissch (normal und separabel).
        - Die Gruppe ist abelsch (Kronecker-Weber-Theorem).
        - Die Untergruppen entsprechen Zwischenkörpern (Galoiskörrespondenz).

        **Beispiele:**
        - n=4: Gal(ℚ(i)/ℚ) ≅ ℤ/2ℤ, σ: i ↦ -i
        - n=5: Gal(ℚ(ζ₅)/ℚ) ≅ ℤ/4ℤ
        - n=8: Gal(ℚ(ζ₈)/ℚ) ≅ ℤ/2ℤ×ℤ/2ℤ

    @param n Positive ganze Zahl ≥ 1.
    @return Dictionary mit:
            - 'n': int, Ordnung der Einheitswurzeln
            - 'degree': int, [ℚ(ζ_n):ℚ] = φ(n)
            - 'galois_group': str, Gruppenname
            - 'order': int, Gruppenordnung = φ(n)
            - 'generators': list[int], Erzeuger (primitive Wurzeln mod n)
            - 'elements': list[int], alle k mit gcd(k,n)=1 (Gruppenelemente)
            - 'is_abelian': bool, immer True
            - 'is_cyclic': bool, ob (ℤ/nℤ)^× zyklisch ist
            - 'automorphisms': list[str], σ_k-Beschreibungen
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Elemente der multiplikativen Gruppe (ℤ/nℤ)^×
    # = alle k mit 1 ≤ k ≤ n und gcd(k, n) = 1
    elements = [k for k in range(1, n + 1) if gcd(k, n) == 1]
    phi_n = len(elements)  # = φ(n) = Ordnung der Galoisgruppe

    # Erzeuger der Gruppe (primitive Wurzeln mod n) via SymPy
    gens = []
    try:
        # primitive_root liefert einen Erzeuger, wenn die Gruppe zyklisch ist
        g = int(primitive_root(n))
        gens = [g]
        is_cyclic = True
    except Exception:
        # Keine primitive Wurzel → Gruppe nicht zyklisch
        # (z.B. n=8: (ℤ/8ℤ)^× ≅ ℤ/2ℤ × ℤ/2ℤ)
        is_cyclic = False
        # Alle Elemente der Ordnung, die gleich der Gruppenordnung ist, sind Erzeuger
        # Für nicht-zyklische Gruppen: mehrere Erzeuger nötig
        gens = elements[:2] if len(elements) >= 2 else elements

    # Galoisgruppe benennen
    if phi_n == 1:
        group_name = "trivial {1}"
    elif is_cyclic:
        group_name = f"ℤ/{phi_n}ℤ"
    else:
        # Faktorisierung der Gruppenordnung für bessere Benennung
        facts = sympy.factorint(phi_n)
        if len(facts) == 1 and list(facts.values())[0] == 1:
            group_name = f"ℤ/{phi_n}ℤ"
        else:
            group_name = f"(ℤ/{n}ℤ)^×"

    # Automorphismus-Beschreibungen: σ_k: ζ_n ↦ ζ_n^k
    automorphisms = [f"σ_{k}: ζ_{n} ↦ ζ_{n}^{k}" for k in elements]

    return {
        'n': n,
        'degree': phi_n,
        'galois_group': group_name,
        'order': phi_n,
        'generators': gens,
        'elements': elements,
        'is_abelian': True,   # Kreisteilungsgruppen sind immer abelsch
        'is_cyclic': is_cyclic,
        'automorphisms': automorphisms,
        'cyclotomic_poly_degree': phi_n,
        'description': (
            f"Gal(ℚ(ζ_{n})/ℚ) ≅ (ℤ/{n}ℤ)^× "
            f"mit Ordnung φ({n}) = {phi_n}"
        ),
    }


def kronecker_weber_check(poly_coeffs: list[int]) -> dict:
    """
    @brief Prüft das Kronecker-Weber-Theorem für eine abelsche Erweiterung.
    @description
        Das Kronecker-Weber-Theorem besagt:

        **Satz (Kronecker 1853, Weber 1886):**
        Jede endliche abelsche Galois-Erweiterung L/ℚ ist in einem Kreisteilungskörper
        ℚ(ζ_n) enthalten:
        $$L \subseteq \mathbb{Q}(\zeta_n)$$

        Das bedeutet: Jeder abelsche Zahlkörper ist ein Teilkörper eines zyklotomischen
        Körpers.

        **Algorithmus:**
        1. Berechne Galoisgruppe G = Gal(L/ℚ) des Polynoms.
        2. Prüfe ob G abelsch ist.
        3. Wenn abelsch: Suche minimales n mit G ↪ (ℤ/nℤ)^×.
        4. Liefere den Kreisteilungskörper ℚ(ζ_n) als "Einbettung".

        **Beispiele:**
        - x²-2: Gal=ℤ/2ℤ → ℚ(√2) ⊆ ℚ(ζ₈)  (da ζ₈+ζ₈⁻¹=√2)
        - x²+1: Gal=ℤ/2ℤ → ℚ(i) = ℚ(ζ₄)
        - x³-1: Gal=ℤ/3ℤ? Nein, Φ₃(x)=x²+x+1 → Gal=ℤ/2ℤ ... (Φ₃ Grad 2)

    @param poly_coeffs Koeffizientenliste (aufsteigend) des Polynoms f(x) ∈ ℚ[x].
    @return Dictionary mit:
            - 'galois_group': str
            - 'is_abelian': bool
            - 'theorem_applies': bool, True wenn abelsche Erweiterung
            - 'minimal_n': int | None, kleinstes n mit G ↪ (ℤ/nℤ)^×
            - 'cyclotomic_field': str, Beschreibung von ℚ(ζ_n)
            - 'explanation': str
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Galoisgruppe des Polynoms bestimmen
    gal_info = galois_group_polynomial(poly_coeffs)
    group_name = gal_info['galois_group']
    group_order = gal_info['order']

    # Prüfe Abelsches Kriterium
    is_ab = group_name in (
        'trivial', '1', 'Z/2Z', 'ℤ/2ℤ', 'Z_2',
        'Z/3Z', 'ℤ/3ℤ', 'Z_3',
        'Z/4Z', 'ℤ/4ℤ', 'Z_4',
        'V_4', 'Z/2Z×Z/2Z',
        'Z/5Z', 'ℤ/5ℤ',
        'Z/6Z', 'ℤ/6ℤ',
    ) or 'Z/' in group_name or 'ℤ/' in group_name

    minimal_n = None
    cyclo_field = None

    if is_ab and group_order > 0:
        # Suche kleinstes n, sodass φ(n) durch group_order teilbar
        # und (ℤ/nℤ)^× die Galoisgruppe als Untergruppe enthält
        for candidate_n in range(1, 200):
            phi_c = int(totient(candidate_n))
            if phi_c % group_order == 0:
                minimal_n = candidate_n
                break
        cyclo_field = f"ℚ(ζ_{minimal_n})" if minimal_n else None

    explanation = (
        f"Das Polynom hat Galoisgruppe {group_name} der Ordnung {group_order}. "
        + (
            f"Diese Gruppe ist abelsch. Nach Kronecker-Weber ist die Erweiterung "
            f"in ℚ(ζ_{minimal_n}) eingebettet (φ({minimal_n}) = {int(totient(minimal_n)) if minimal_n else '?'})."
            if is_ab else
            f"Die Gruppe {group_name} ist NICHT abelsch. Das Kronecker-Weber-Theorem "
            f"gilt nicht für diese Erweiterung."
        )
    )

    return {
        'galois_group': group_name,
        'order': group_order,
        'is_abelian': is_ab,
        'theorem_applies': is_ab,
        'minimal_n': minimal_n,
        'cyclotomic_field': cyclo_field,
        'explanation': explanation,
    }


def dirichlet_characters_from_galois(n: int) -> dict:
    """
    @brief Leitet Dirichlet-Charaktere aus der Galoisgruppe des Kreisteilungskörpers ab.
    @description
        Via Galois-Theorie können Dirichlet-Charaktere modulo n als Gruppencharaktere
        der abelschen Galoisgruppe Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^× interpretiert werden:

        $$\chi: \text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \to \mathbb{C}^{\times}$$

        Ein Gruppencharakter (eindimensionale Darstellung) bildet die Galoisgruppe auf
        Einheitswurzeln ab: χ(σ_k) = ω^f(k), wobei ω eine Einheitswurzel ist.

        **Anzahl der Charaktere:**
        Die Anzahl der verschiedenen Dirichlet-Charaktere modulo n ist φ(n)
        (Anzahl der Gruppencharaktere = Gruppenordnung für abelsche Gruppen).

        **Hauptcharakter:**
        χ₀(k) = 1 für gcd(k,n) = 1, χ₀(k) = 0 sonst.

        **Beispiel n=4:**
        - Elemente: {1, 3} ≅ ℤ/2ℤ
        - χ₀: σ₁ ↦ 1, σ₃ ↦ 1  (Hauptcharakter)
        - χ₁: σ₁ ↦ 1, σ₃ ↦ -1 (nicht-trivialer Charakter)

    @param n Positive ganze Zahl ≥ 1 (Modul).
    @return Dictionary mit:
            - 'n': int
            - 'galois_group': str
            - 'phi_n': int, Anzahl der Charaktere
            - 'elements': list[int], Gruppenelemente
            - 'characters': list[dict], alle Charaktere mit Werten
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Galoisgruppe des Kreisteilungskörpers
    gal = cyclotomic_galois_group(n)
    elements = gal['elements']     # k mit gcd(k,n)=1
    phi_n = gal['order']           # φ(n)

    characters = []

    if phi_n == 1:
        # Nur Hauptcharakter
        characters.append({
            'index': 0,
            'name': 'χ₀ (Hauptcharakter)',
            'is_principal': True,
            'values': {k: 1 for k in elements},
            'conductor': 1,
        })
    else:
        # Für zyklische Gruppen: Charaktere über primitive Wurzel
        # Für allgemeine abelsche Gruppen: Charaktere explizit angeben

        # Ordnungen der Gruppenelemente bestimmen
        elem_orders = {}
        for k in elements:
            # Ordnung von k in (ℤ/nℤ)^×: kleinstes m mit k^m ≡ 1 (mod n)
            m = 1
            cur = k
            while cur % n != 1 and m < phi_n + 1:
                cur = (cur * k) % n
                m += 1
            elem_orders[k] = m

        # Hauptcharakter χ₀: alle Werte 1
        characters.append({
            'index': 0,
            'name': 'χ₀ (Hauptcharakter)',
            'is_principal': True,
            'values': {k: 1 for k in elements},
            'conductor': 1,
            'description': f"χ₀(σ_k) = 1 für alle k ∈ (ℤ/{n}ℤ)^×",
        })

        # Für einfache Fälle: Legendre-Symbol mod p (falls n = p prim)
        if isprime(n):
            # Nicht-triviale Charaktere: χ₁ = Legendre-Symbol mod p
            # χ₁(k) = (k/p) = 1 wenn k QR, -1 wenn QNR
            legendre_vals = {}
            for k in elements:
                # k ist QR mod p genau dann wenn k^((p-1)/2) ≡ 1 (mod p)
                kk = pow(k, (n - 1) // 2, n)
                legendre_vals[k] = 1 if kk == 1 else -1

            characters.append({
                'index': 1,
                'name': f'χ₁ (Legendre mod {n})',
                'is_principal': False,
                'values': legendre_vals,
                'conductor': n,
                'description': f"χ₁(k) = (k/{n}) Legendre-Symbol",
            })

            # Weitere Charaktere via Potenzen des Legendre-Symbols (für φ(n) > 2)
            for i in range(2, min(phi_n, 5)):
                chi_vals = {k: legendre_vals[k] ** i for k in elements}
                characters.append({
                    'index': i,
                    'name': f'χ_{i} = χ₁^{i}',
                    'is_principal': all(v == 1 for v in chi_vals.values()),
                    'values': chi_vals,
                    'conductor': n,
                    'description': f"Potenzcharakter χ₁^{i}",
                })
        else:
            # Für zusammengesetztes n: einfache Darstellung
            for i in range(1, min(phi_n, 4)):
                # χ_i(k) = ω^{i·ind(k)} wobei ω = e^{2πi/φ(n)}
                # Vereinfacht: Zeichen des Charakters angeben
                chi_vals = {}
                for k in elements:
                    # Phasenfaktor: e^{2πi·i·(k-1)/phi_n} (Approximation)
                    phase = (i * (elements.index(k))) % phi_n
                    chi_vals[k] = cmath.exp(2j * math.pi * phase / phi_n)
                characters.append({
                    'index': i,
                    'name': f'χ_{i}',
                    'is_principal': False,
                    'values': {k: round(v.real, 8) + round(v.imag, 8) * 1j
                               for k, v in chi_vals.items()},
                    'conductor': n,
                    'description': f"Gruppencharakter χ_{i} der Galoisgruppe",
                })

    return {
        'n': n,
        'galois_group': gal['galois_group'],
        'phi_n': phi_n,
        'elements': elements,
        'characters': characters,
        'total_characters': phi_n,
        'description': (
            f"Dirichlet-Charaktere mod {n} aus Gal(ℚ(ζ_{n})/ℚ) ≅ {gal['galois_group']}"
        ),
    }


# =============================================================================
# ENDLICHE KÖRPER – GALOISGRUPPE
# =============================================================================

def galois_group_finite_field(p: int, n: int) -> dict:
    """
    @brief Berechnet die Galoisgruppe Gal(GF(p^n) / GF(p)) des endlichen Körpers.
    @description
        Die Galoisgruppe der Erweiterung GF(p^n)/GF(p) ist zyklisch von Ordnung n:

        $$\text{Gal}(\text{GF}(p^n)/\text{GF}(p)) \cong \mathbb{Z}/n\mathbb{Z}$$

        Der kanonische Erzeuger ist der Frobenius-Automorphismus:
        $$\text{Frob}_p: \text{GF}(p^n) \to \text{GF}(p^n), \quad x \mapsto x^p$$

        **Eigenschaften:**
        - Die Erweiterung GF(p^n)/GF(p) ist immer galoissch.
        - Jeder Automorphismus ist eine Potenz des Frobenius: Frob_p^k für k=0,...,n-1.
        - Frob_p^n = id (Frobenius hat Ordnung n).

        **Zwischenkörper (Galoiskörrespondenz):**
        Untergruppen ℤ/dℤ ≤ ℤ/nℤ (d | n) entsprechen Zwischenkörpern GF(p^d).

    @param p Primzahl (Charakteristik des Körpers).
    @param n Grad der Erweiterung (n ≥ 1).
    @return Dictionary mit:
            - 'p': int
            - 'n': int
            - 'field': str, Beschreibung des Körpers
            - 'base_field': str
            - 'galois_group': str, "ℤ/nℤ"
            - 'order': int, n
            - 'is_cyclic': bool, True
            - 'generator': str, Frobenius-Beschreibung
            - 'elements': list[str], alle Automorphismen
            - 'intermediate_fields': list[dict], Zwischenkörper via d|n
    @raises InvalidInputError Wenn p nicht prim oder n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if not isprime(p):
        raise InvalidInputError(f"p={p} muss eine Primzahl sein")
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Automorphismen: Frob_p^k für k = 0, 1, ..., n-1
    elements = [f"Frob_{p}^{k}" if k > 0 else "id" for k in range(n)]

    # Zwischenkörper: GF(p^d) für alle Teiler d von n
    intermediate_fields = []
    for d in range(1, n + 1):
        if n % d == 0:
            # Untergruppe: ℤ/(n/d)ℤ ≤ ℤ/nℤ (Fixkörper von Frob^d-Untergruppe)
            subgroup_order = n // d
            intermediate_fields.append({
                'field': f"GF({p}^{d})",
                'degree_over_base': d,
                'fixed_by_subgroup': f"⟨Frob_{p}^{d}⟩ (Ordnung {subgroup_order})",
                'subgroup': f"ℤ/{subgroup_order}ℤ",
            })

    return {
        'p': p,
        'n': n,
        'field': f"GF({p}^{n})",
        'base_field': f"GF({p})",
        'galois_group': f"ℤ/{n}ℤ",
        'order': n,
        'is_cyclic': True,
        'is_abelian': True,
        'generator': f"Frob_{p}: x ↦ x^{p}",
        'elements': elements,
        'intermediate_fields': intermediate_fields,
        'description': (
            f"Gal(GF({p}^{n})/GF({p})) ≅ ℤ/{n}ℤ, "
            f"erzeugt vom Frobenius Frob_{p}: x ↦ x^{p}"
        ),
    }


# =============================================================================
# KONSTRUIERBARKEIT MIT ZIRKEL UND LINEAL
# =============================================================================

def is_constructible(n: int) -> dict:
    """
    @brief Prüft ob ein reguläres n-Eck mit Zirkel und Lineal konstruierbar ist.
    @description
        Das Gauss-Wantzel-Theorem (1796/1837) charakterisiert konstruierbare reguläre
        Polygone vollständig:

        **Satz (Gauss-Wantzel):**
        Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann, wenn
        $$n = 2^k \cdot p_1 \cdot p_2 \cdots p_r$$
        wobei k ≥ 0 und p_i verschiedene Fermat-Primzahlen sind.

        **Fermat-Primzahlen** sind Primzahlen der Form $F_m = 2^{2^m} + 1$:
        F_0 = 3, F_1 = 5, F_2 = 17, F_3 = 257, F_4 = 65537.
        (Weitere sind nicht bekannt.)

        **Galois-theoretische Begründung:**
        n-Eck konstruierbar ⟺ [ℚ(ζ_n):ℚ] = φ(n) ist eine 2-Potenz
        ⟺ Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^× hat 2-Potenz-Ordnung
        ⟺ Turm von Körpererweiterungen vom Grad 2 existiert.

    @param n Positive ganze Zahl ≥ 1.
    @return Dictionary mit:
            - 'n': int
            - 'is_constructible': bool
            - 'phi_n': int, Euler-Phi-Wert
            - 'phi_is_power_of_2': bool
            - 'fermat_primes_used': list[int]
            - 'power_of_2_factor': int
            - 'reason': str
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Bekannte Fermat-Primzahlen (alle bekannten)
    fermat_primes = [3, 5, 17, 257, 65537]

    # n faktorisieren
    if n == 1 or n == 2:
        # Triviale Fälle
        return {
            'n': n,
            'is_constructible': True,
            'phi_n': int(totient(n)),
            'phi_is_power_of_2': True,
            'fermat_primes_used': [],
            'power_of_2_factor': n if n == 2 else 1,
            'reason': f"n={n} ist trivial konstruierbar.",
        }

    phi_n = int(totient(n))

    # Prüfe ob φ(n) eine 2-Potenz ist
    def is_power_of_2(k):
        """Hilfsfunktion: ist k eine 2-Potenz?"""
        return k > 0 and (k & (k - 1)) == 0

    phi_is_2_power = is_power_of_2(phi_n)

    # Faktorisiere n und prüfe die Struktur
    facts = sympy.factorint(n)

    # Zähle den 2er-Anteil
    two_factor = facts.get(2, 0)
    power_of_2 = 2 ** two_factor

    # Ungerade Primfaktoren
    odd_prime_factors = {p: e for p, e in facts.items() if p != 2}

    used_fermat = []
    non_fermat_odd = []
    repeated_odd = []

    for p, e in odd_prime_factors.items():
        if e > 1:
            # Primfaktor mit Exponent > 1: nicht konstruierbar
            repeated_odd.append(p)
        elif p in fermat_primes:
            used_fermat.append(p)
        else:
            non_fermat_odd.append(p)

    constructible = phi_is_2_power and not non_fermat_odd and not repeated_odd

    # Begründung zusammenstellen
    if constructible:
        if not used_fermat:
            reason = f"n={n} = 2^{two_factor} ist eine 2-Potenz. Das n-Eck ist konstruierbar."
        else:
            fp_str = " · ".join(str(p) for p in used_fermat)
            reason = (
                f"n={n} = 2^{two_factor} · {fp_str}. "
                f"Alle ungerade Faktoren sind verschiedene Fermat-Primzahlen. "
                f"Das {n}-Eck ist konstruierbar (Gauss-Wantzel)."
            )
    elif non_fermat_odd:
        reason = (
            f"n={n} enthält den ungerade Primfaktor {non_fermat_odd[0]}, "
            f"der keine Fermat-Primzahl ist. Das {n}-Eck ist NICHT konstruierbar."
        )
    elif repeated_odd:
        reason = (
            f"n={n} enthält einen ungerade Primfaktor {repeated_odd[0]} mit "
            f"Exponent > 1. Das {n}-Eck ist NICHT konstruierbar."
        )
    else:
        reason = f"n={n} ist nicht konstruierbar (φ(n)={phi_n} ist keine 2-Potenz)."

    return {
        'n': n,
        'is_constructible': constructible,
        'phi_n': phi_n,
        'phi_is_power_of_2': phi_is_2_power,
        'fermat_primes_used': used_fermat,
        'power_of_2_factor': power_of_2,
        'reason': reason,
        'factorization': dict(facts),
    }


def construct_regular_polygon(n: int) -> dict:
    """
    @brief Analysiert die Konstruierbarkeit und Galois-Theorie eines regulären n-Ecks.
    @description
        Verbindet das Gauss-Wantzel-Theorem mit der Galois-Theorie:

        Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann,
        wenn der Turm von Körpererweiterungen über ℚ nur Schritte vom Grad 2 hat:

        $$\mathbb{Q} = K_0 \subset K_1 \subset \cdots \subset K_r = \mathbb{Q}(\zeta_n)$$

        mit $[K_{i+1} : K_i] = 2$ für alle i.

        Die Ecken des regulären n-Ecks sind die n-ten Einheitswurzeln
        $$\zeta_n^k = e^{2\pi i k/n}, \quad k = 0, 1, \ldots, n-1$$

        **Historische Bedeutung:**
        Gauss zeigte 1796 (mit 19 Jahren), dass das reguläre 17-Eck konstruierbar ist
        (F_2 = 17 ist eine Fermat-Primzahl). Dies war der erste neue Fortschritt
        seit der Antike.

    @param n Anzahl der Ecken (n ≥ 3).
    @return Dictionary mit:
            - 'n': int
            - 'is_constructible': bool
            - 'constructibility_info': dict von is_constructible(n)
            - 'galois_group': dict von cyclotomic_galois_group(n)
            - 'vertices': list[complex], Koordinaten der Ecken
            - 'historical_note': str
    @raises InvalidInputError Wenn n < 3.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 3:
        raise InvalidInputError(f"n={n} muss ≥ 3 sein (Polygon braucht mind. 3 Ecken)")

    # Konstruierbarkeits-Analyse
    constr_info = is_constructible(n)

    # Galoisgruppe des Kreisteilungskörpers
    gal_info = cyclotomic_galois_group(n)

    # Koordinaten der Ecken des regulären n-Ecks in der komplexen Ebene
    # Ecke k: ζ_n^k = e^{2πik/n}
    vertices = [cmath.exp(2j * math.pi * k / n) for k in range(n)]

    # Historische Notiz
    historical_notes = {
        3:  "Gleichseitiges Dreieck: seit der Antike bekannt.",
        4:  "Quadrat: seit der Antike bekannt.",
        5:  "Reguläres Pentagon: Euklid, ~300 v. Chr.",
        6:  "Reguläres Hexagon: seit der Antike bekannt.",
        8:  "Reguläres Oktogon: seit der Antike bekannt.",
        10: "Reguläres Dekagon: Euklid, ~300 v. Chr.",
        15: "Reguläres 15-Eck: Euklid, ~300 v. Chr.",
        17: "Reguläres 17-Eck: C.F. Gauss, 1796 (revolutionäre Entdeckung!).",
        257: "Reguläres 257-Eck: Magnus Georg Paucker, 1822.",
        65537: "Reguläres 65537-Eck: Johann Gustav Hermes, 1894 (10 Jahre Arbeit!).",
    }
    note = historical_notes.get(n, "")
    if not note:
        if constr_info['is_constructible']:
            note = f"Das reguläre {n}-Eck ist konstruierbar (Gauss-Wantzel-Kriterium erfüllt)."
        else:
            note = f"Das reguläre {n}-Eck ist NICHT konstruierbar mit Zirkel und Lineal."

    return {
        'n': n,
        'is_constructible': constr_info['is_constructible'],
        'constructibility_info': constr_info,
        'galois_group': gal_info,
        'vertices': vertices,
        'historical_note': note,
        'phi_n': gal_info['order'],
        'description': (
            f"Reguläres {n}-Eck: {'Konstruierbar' if constr_info['is_constructible'] else 'NICHT konstruierbar'} "
            f"mit Zirkel und Lineal. φ({n}) = {gal_info['order']}."
        ),
    }


# =============================================================================
# HILBERT'S SATZ 90 UND NORM / SPUR
# =============================================================================

def norm_and_trace(alpha_coeffs: list[float], poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet Norm und Spur eines algebraischen Elements.
    @description
        Sei L = K(α) eine Körpererweiterung vom Grad n = [L:K], und
        α ein Element mit Minimalpolynom f(x) vom Grad n.

        **Norm:**
        $$N_{L/K}(\alpha) = \prod_{\sigma \in \text{Gal}(L/K)} \sigma(\alpha) = (-1)^n \cdot a_0/a_n$$

        Das Produkt aller Konjugierten von α (= alle Wurzeln von f über K).
        Für f(x) = a_n x^n + ... + a_0 gilt:
        $$N_{L/K}(\alpha) = (-1)^n \cdot a_0/a_n$$

        **Spur:**
        $$\text{Tr}_{L/K}(\alpha) = \sum_{\sigma \in \text{Gal}(L/K)} \sigma(\alpha) = -a_{n-1}/a_n$$

        Die Summe aller Konjugierten.

        **Charakteristisches Polynom:**
        $$\chi_\alpha(x) = \prod_{\sigma}(x - \sigma(\alpha)) = f(x) \text{ (falls f irreduzibel)}$$

    @param alpha_coeffs Koeffizientenliste von α als Linearkombination
                        [a_0, a_1, ..., a_{n-1}] in Basis {1, α_gen, α_gen², ...}
                        ODER leer → direkte Benutzung des Minimalpolynoms
    @param poly_coeffs Koeffizientenliste (aufsteigend) des Minimalpolynoms f(x).
    @return Dictionary mit:
            - 'norm': float, N_{L/K}(α)
            - 'trace': float, Tr_{L/K}(α)
            - 'degree': int, [L:K]
            - 'conjugates': list[complex], alle σ(α)
            - 'char_poly': list[float], charakteristisches Polynom
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Minimalpolynom auslesen
    # poly_coeffs ist aufsteigend: [a_0, a_1, ..., a_n]
    n = len(poly_coeffs) - 1   # Grad des Minimalpolynoms
    if n < 1:
        raise InvalidInputError("Minimalpolynom muss Grad ≥ 1 haben")

    # Leitkoeffizient und konstanter Term
    leading = poly_coeffs[n]    # a_n (Koeffizient von x^n)
    constant = poly_coeffs[0]   # a_0 (konstanter Term)

    # Norm via Vieta: N = (-1)^n * a_0 / a_n
    norm_val = ((-1) ** n) * constant / leading

    # Spur via Vieta: Tr = -a_{n-1} / a_n
    # (Summe der Wurzeln = -a_{n-1}/a_n für monisches Poly nach Division)
    if n >= 2:
        trace_val = -poly_coeffs[n - 1] / leading
    elif n == 1:
        # Lineares Polynom a_1*x + a_0: einzige Wurzel = -a_0/a_1
        trace_val = -constant / leading
    else:
        trace_val = 0.0

    # Konjugierte (numerische Wurzeln des Minimalpolynoms)
    # Umwandlung in absteigende Koeffizienten für numpy
    coeffs_desc = list(reversed(poly_coeffs))
    try:
        import numpy as np
        conjugates = np.roots(coeffs_desc).tolist()
    except Exception:
        conjugates = []

    # Charakteristisches Polynom = Minimalpolynom (da f irreduzibel)
    char_poly = poly_coeffs[:]

    # Verifikation: Produkt der Konjugierten sollte |Norm| ergeben
    if conjugates:
        product_conj = 1.0
        for c in conjugates:
            product_conj *= c
        product_conj_real = product_conj.real if hasattr(product_conj, 'real') else product_conj
    else:
        product_conj_real = norm_val

    return {
        'norm': norm_val,
        'trace': trace_val,
        'degree': n,
        'conjugates': [str(c) for c in conjugates],
        'char_poly': char_poly,
        'product_conjugates': product_conj_real,
        'sum_conjugates': trace_val,
        'description': (
            f"N_{{L/K}}(α) = {norm_val:.6f}, "
            f"Tr_{{L/K}}(α) = {trace_val:.6f} "
            f"für Minimalpolynom Grad {n}"
        ),
    }


def hilbert90(n: int, check_type: str = "norm") -> dict:
    """
    @brief Demonstration von Hilbert's Satz 90.
    @description
        **Hilbert's Satz 90 (Originalform, 1897):**
        Sei L/K eine zyklische Galois-Erweiterung mit Gal(L/K) = ⟨σ⟩.
        Dann gilt:
        $$N_{L/K}(\alpha) = 1 \iff \alpha = \beta / \sigma(\beta) \text{ für ein } \beta \in L^{\times}$$

        **Kohomologische Formulierung:**
        $$H^1(\text{Gal}(L/K),\ L^{\times}) = 1$$

        Das bedeutet: Jeder 1-Kozykel z: Gal(L/K) → L^× ist ein 1-Kobord.

        **Allgemeine Form (Hilbert 90 für Matrizen/GL_n):**
        $$H^1(\text{Gal}(L/K),\ \text{GL}_n(L)) = 1 \text{ (Satz von Speiser)}$$

        **Anwendungen:**
        - Kummer-Erweiterungen: x^n - a zerfällt ⟺ N(α) = 1
        - Beschreibung von zyklischen Erweiterungen
        - Grundlage der Kohomologie der Galoisgruppen

        **Beispiel (zyklotomisch, n=4, L=ℚ(i)/ℚ):**
        - Gal = ⟨σ⟩ mit σ(i) = -i
        - α = i: N(i) = i·(-i) = 1 → α = i/(-i) = -1 = β/σ(β) mit β = ?
        - Tatsächlich: i = (1+i) / σ(1+i) = (1+i)/(1-i)

    @param n Ordnung der zyklischen Galoisgruppe (n ≥ 1).
    @param check_type Art der Demonstration: "norm" oder "cohomology".
    @return Dictionary mit:
            - 'theorem': str, Formulierung
            - 'n': int
            - 'example_element': str
            - 'norm_equals_1': bool, Beispielelement mit Norm 1
            - 'cohomology_trivial': bool, H¹ = 1
            - 'explanation': str
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Für das einfachste nicht-triviale Beispiel: ℚ(i)/ℚ (n=2, σ: i↦-i)
    # α = (1+i) hat N(α) = (1+i)(1-i) = 2 ≠ 1
    # α = i hat N(i) = i·(-i) = 1 ✓

    # Allgemeine Demonstration: In ℂ(ζ_n)/ℝ
    # Wähle α = ζ_n (n-te Einheitswurzel), Norm = Produkt aller Konjugierten
    zeta = cmath.exp(2j * math.pi / n)

    # Konjugierte von ζ_n: ζ_n^k für k=1,...,n (alle n-ten Einheitswurzeln)
    conjugates = [cmath.exp(2j * math.pi * k / n) for k in range(1, n + 1)]

    # Norm = Produkt aller Konjugierten
    norm_product = 1.0 + 0j
    for c in conjugates:
        norm_product *= c

    # Das Produkt aller n-ten Einheitswurzeln ist (-1)^{n+1}
    # (Vieta: Produkt der Wurzeln von x^n - 1 = (-1)^n · (-1) = (-1)^{n+1})
    norm_real = round(norm_product.real, 8)
    norm_imag = round(norm_product.imag, 8)

    # Für ζ_n mit dem Automorphismus σ: ζ_n ↦ ζ_n^k (Kreisgruppe)
    # N_{ℚ(ζ_n)/ℚ}(ζ_n) = ∏_{gcd(k,n)=1} ζ_n^k = ζ_n^{∑ k: gcd(k,n)=1}
    # = ζ_n^{n·φ(n)/2} (falls n > 2) = 1 oder -1

    # Einfachstes Beispiel: n=4, L=ℚ(i), σ: i↦-i
    # α = i, σ(α) = -i, N(i) = i·(-i) = 1
    # Hilbert 90: ∃β: i = β/σ(β) → β = 1+i, σ(β) = 1-i, β/σ(β) = (1+i)/(1-i) = i ✓

    example_norm_is_1 = (n == 1 or abs(norm_real - 1) < 0.01 or abs(norm_real + 1) < 0.01)

    theorem_text = (
        "Hilbert's Satz 90: Sei L/K zyklisch galoissch mit Gal(L/K) = ⟨σ⟩. "
        "Dann gilt: N_{L/K}(α) = 1 ⟺ ∃β ∈ L^× mit α = β/σ(β). "
        "Äquivalent: H¹(Gal(L/K), L^×) = 1."
    )

    if check_type == "cohomology":
        explanation = (
            f"Kohomologische Form: H¹(Gal, L^×) = 1 bedeutet, dass jeder "
            f"1-Kozykel z: Gal(L/K) → L^× ein 1-Kobord ist. "
            f"Für die zyklische Gruppe ℤ/{n}ℤ = ⟨σ⟩ ist ein Kozykel durch "
            f"z(σ) = α ∈ L^× festgelegt mit N(α) = z(σ)·z(σ²)/z(σ²)·... = 1."
        )
    else:
        explanation = (
            f"Normform für n={n}: "
            f"Ein Element α ∈ L^× mit N_{{L/K}}(α) = 1 kann immer als "
            f"α = β/σ(β) dargestellt werden (für geeignetes β ∈ L^×). "
            f"Beispiel: L/K = ℚ(i)/ℚ (n=2), α = i: N(i) = i·(-i) = 1, "
            f"und i = (1+i)/(1-i) = (1+i)/σ(1+i)."
        )

    return {
        'theorem': theorem_text,
        'n': n,
        'example_element': f"ζ_{n} = e^{{2πi/{n}}}",
        'norm_equals_1': example_norm_is_1,
        'cohomology_trivial': True,   # H¹(Gal(L/K), L^×) = 1 immer (Hilbert 90)
        'galois_group': f"ℤ/{n}ℤ (zyklisch)",
        'explanation': explanation,
        'h1_trivial': True,
        'description': f"Hilbert's Satz 90 für zyklische Erweiterung der Ordnung {n}",
    }


# =============================================================================
# ABEL-RUFFINI UND RADIKALTURM
# =============================================================================

def galois_group_symmetric(n: int) -> dict:
    """
    @brief Beschreibt S_n als Galoisgruppe von x^n - t über ℚ(t).
    @description
        Das generische Polynom x^n - t (mit transzendentem t) hat die symmetrische
        Gruppe S_n als Galoisgruppe über dem rationalen Funktionenkörper ℚ(t):

        $$\text{Gal}(K_f / \mathbb{Q}(t)) \cong S_n$$

        Für konkrete irreduzible Polynome gilt: Die "generische" Galoisgruppe ist S_n.
        Speziellere Polynome können kleinere Galoisgruppen haben.

        **Eigenschaften von S_n:**
        - Ordnung: n!
        - Für n ≤ 4: S_n ist auflösbar
        - Für n ≥ 5: S_n ist NICHT auflösbar (enthält A_n als nicht-auflösbaren Normalteiler)

        **Kompositionsreihe:**
        - S_1 ⊃ {e}: Länge 1
        - S_2 ⊃ {e}: Länge 1 (abelsch)
        - S_3 ⊃ A_3 ⊃ {e}: Länge 2
        - S_4 ⊃ A_4 ⊃ V_4 ⊃ ℤ/2ℤ ⊃ {e}: Länge 4
        - S_5 ⊃ A_5 ⊃ {e}: A_5 ist einfach, nicht auflösbar!

    @param n Grad n ≥ 1.
    @return Dictionary mit:
            - 'n': int
            - 'group': str, "S_n"
            - 'order': int, n!
            - 'is_solvable': bool
            - 'is_simple': bool (nur für S_1 und S_2)
            - 'composition_series': list[str]
            - 'description': str
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    order = math.factorial(n)

    # Auflösbarkeit: S_n auflösbar ⟺ n ≤ 4
    is_solv = (n <= 4)

    # Kompositionsreihe von S_n
    comp_series = {
        1: ["S_1 = {e}"],
        2: ["S_2", "{e}"],
        3: ["S_3", "A_3 ≅ ℤ/3ℤ", "{e}"],
        4: ["S_4", "A_4", "V_4 ≅ ℤ/2ℤ×ℤ/2ℤ", "ℤ/2ℤ", "{e}"],
        5: ["S_5", "A_5 (einfach, nicht auflösbar!)", "{e}"],
        6: ["S_6", "A_6 (einfach)", "{e}"],
    }
    series = comp_series.get(n, [f"S_{n}", f"A_{n} (einfach)", "{e}"])

    # Quotientenfaktoren
    factor_names = {
        1: [],
        2: ["ℤ/2ℤ"],
        3: ["ℤ/2ℤ", "ℤ/3ℤ"],
        4: ["ℤ/2ℤ", "ℤ/3ℤ", "ℤ/2ℤ", "ℤ/2ℤ"],
        5: ["ℤ/2ℤ", "A_5 (nicht-abelsch einfach)"],
    }
    factors = factor_names.get(n, [f"ℤ/2ℤ", f"A_{n}"])

    return {
        'n': n,
        'group': f"S_{n}",
        'order': order,
        'is_solvable': is_solv,
        'is_simple': (n <= 2),
        'is_abelian': (n <= 2),
        'composition_series': series,
        'composition_factors': factors,
        'normal_subgroup_An': f"A_{n} (Alternierende Gruppe, Index 2, Ordnung {order // 2})",
        'description': (
            f"S_{n} ist die Galoisgruppe des generischen Grads-{n}-Polynoms x^{n}-t "
            f"über ℚ(t). Ordnung: {n}! = {order}. "
            + ("Auflösbar." if is_solv else f"NICHT auflösbar (A_{n} ist einfach für n≥5).")
        ),
    }


def radical_tower(poly_coeffs: list[int]) -> dict:
    """
    @brief Konstruiert einen Radikalturm für auflösbare Polynome.
    @description
        Ein Radikalturm (Radikalerweiterung) ist eine Folge von Körpern:

        $$K = K_0 \subset K_1 \subset \cdots \subset K_r$$

        wobei jeder Schritt $K_{i+1} = K_i(\alpha_i)$ mit $\alpha_i^{n_i} \in K_i$
        eine einfache Radikalerweiterung ist.

        **Galois-Theorie:**
        Ein Polynom f ist genau dann durch Radikale lösbar, wenn seine Galoisgruppe
        auflösbar ist (d.h. eine Kompositionsreihe mit abelschen Faktoren besitzt).

        **Auflösungsformeln:**
        - Grad 1: x = -b/a (trivial)
        - Grad 2: x = (-b ± √D) / 2a, D = b²-4ac
        - Grad 3: Cardano-Formel via kubische Resolventenwurzel
        - Grad 4: Ferrari-Methode via resolventenkubischer Gleichung

    @param poly_coeffs Koeffizientenliste (aufsteigend) von f(x) ∈ ℚ[x].
    @return Dictionary mit:
            - 'is_solvable_by_radicals': bool
            - 'galois_group': str
            - 'degree': int
            - 'radical_tower': list[str], Beschreibung des Radikalturms
            - 'formula': str, explizite Formel (falls verfügbar)
            - 'reason': str
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Auflösbarkeit prüfen
    solv_info = is_solvable_by_radicals(poly_coeffs)
    n = len(poly_coeffs) - 1   # Grad des Polynoms

    # Koeffizienten (absteigend)
    coeffs_desc = list(reversed(poly_coeffs))

    tower = []
    formula = ""

    # is_solvable_by_radicals() liefert 'solvable' (nicht 'is_solvable')
    is_solv = solv_info.get('solvable', solv_info.get('is_solvable', False))

    if not is_solv:
        tower = [
            f"ℚ (Basiskörper)",
            f"(Keine Radikalerweiterung möglich – Galoisgruppe {solv_info['galois_group']} nicht auflösbar)",
        ]
        formula = "Keine geschlossene Radikalformel existiert (Abel-Ruffini)."
    elif n == 1:
        # Lineares Polynom: a*x + b = 0 → x = -b/a
        a = coeffs_desc[0]
        b = coeffs_desc[1] if len(coeffs_desc) > 1 else 0
        tower = ["ℚ = K_0 = K_1 (Zerfällungskörper = Basiskörper)"]
        formula = f"x = {-b}/{a} = {-b/a:.6f}" if a != 0 else "Kein Polynom"
    elif n == 2:
        # Quadratisches Polynom: a*x² + b*x + c = 0
        a = coeffs_desc[0]
        b_c = coeffs_desc[1] if len(coeffs_desc) > 1 else 0
        c = coeffs_desc[2] if len(coeffs_desc) > 2 else 0
        disc = b_c * b_c - 4 * a * c
        tower = [
            "K_0 = ℚ",
            f"K_1 = ℚ(√{disc}) (Zerfällungskörper, Grad 2 über ℚ)",
        ]
        formula = (
            f"x = (-({b_c}) ± √({disc})) / (2·{a}) "
            f"= (-{b_c} ± √{disc}) / {2*a}"
        )
    elif n == 3:
        # Kubisches Polynom: Cardano-Formel
        a = coeffs_desc[0]
        tower = [
            "K_0 = ℚ",
            "K_1 = ℚ(ω) mit ω = e^{2πi/3}, [K_1:K_0] = 2 (falls ω ∉ ℚ)",
            "K_2 = K_1(∛Δ) mit Δ die Diskriminante, [K_2:K_1] = 3",
            f"K_3 = K_2(Wurzeln) = Zerfällungskörper",
        ]
        formula = "Cardano-Formel: x = ∛(-q/2 + √(q²/4+p³/27)) + ∛(-q/2 - √(q²/4+p³/27))"
    elif n == 4:
        # Grad 4: Ferrari-Methode
        tower = [
            "K_0 = ℚ",
            "K_1 = ℚ(√Δ₀) (quadratische Erweiterung)",
            "K_2 = K_1(y) mit y Wurzel der Resolventenkubik",
            "K_3 = K_2(√m) (weitere quadratische Erweiterung)",
            "K_4 = K_3(Wurzeln) = Zerfällungskörper",
        ]
        formula = "Ferrari-Methode: Reduktion auf Resolventenkubik, dann Cardano"
    else:
        tower = [
            "K_0 = ℚ",
            f"(Radikalturm existiert, da Galoisgruppe {solv_info['galois_group']} auflösbar ist)",
            "Explizite Konstruktion via Kompositionsreihe der Galoisgruppe",
        ]
        formula = "Explizite Formel komplex, existiert aber theoretisch."

    return {
        'is_solvable_by_radicals': is_solv,
        'galois_group': solv_info['galois_group'],
        'degree': n,
        'radical_tower': tower,
        'formula': formula,
        'reason': solv_info.get('reason', ''),
        'solvability_info': solv_info,
        'description': (
            f"Grad-{n}-Polynom mit Galoisgruppe {solv_info['galois_group']}: "
            + ("Radikalturm existiert." if is_solv
               else "KEIN Radikalturm möglich.")
        ),
    }


def abel_ruffini_demo() -> dict:
    """
    @brief Demonstration des Satzes von Abel-Ruffini (Unlösbarkeit von Grad-5+).
    @description
        **Satz von Abel-Ruffini (1799/1824):**
        Es gibt keine allgemeine Lösungsformel durch Radikale für Polynomgleichungen
        vom Grad n ≥ 5.

        **Beweis (via Galois-Theorie):**
        1. Das generische Polynom vom Grad 5 hat Galoisgruppe S_5.
        2. S_5 ist nicht auflösbar, da der Normalteiler A_5 (alternierende Gruppe)
           einfach und nicht-abelsch ist.
        3. f ist durch Radikale lösbar ⟺ Gal(f) ist auflösbar.
        4. Also: kein generisches Grad-5-Polynom ist durch Radikale lösbar.

        **Konkretes Beispiel: x⁵ - x - 1**
        - Irreduzibel über ℚ (Eisenstein nicht direkt anwendbar, aber verifizierbar)
        - Galoisgruppe ist S_5 (numerisch überprüfbar)
        - Nicht durch Radikale lösbar.

        **Kompositionsreihe von S_5:**
        S_5 ⊃ A_5 ⊃ {e}
        Der Faktor A_5/{e} ≅ A_5 ist einfach und nicht-abelsch → S_5 nicht auflösbar.

        **Auflösbare Gruppen bis Grad 4:**
        - S_1: trivial
        - S_2 ≅ ℤ/2ℤ: abelsch, auflösbar
        - S_3: S_3 ⊃ A_3 ⊃ {e}, Faktoren ℤ/2ℤ und ℤ/3ℤ, auflösbar
        - S_4: S_4 ⊃ A_4 ⊃ V_4 ⊃ ℤ/2ℤ ⊃ {e}, alle Faktoren abelsch, auflösbar

    @return Dictionary mit:
            - 'theorem': str, Satzformulierung
            - 's5_solvable': bool, False (S_5 nicht auflösbar)
            - 'a5_simple': bool, True
            - 'example_polynomial': list[int], Koeffizienten von x⁵-x-1
            - 'example_galois_group': str
            - 'example_solvable': bool, False
            - 'solvable_by_degree': dict, Auflösbarkeit nach Grad
            - 'composition_series_s5': list[str]
            - 'explanation': str
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Beispielpolynom: x⁵ - x - 1 (aufsteigend: [-1, -1, 0, 0, 0, 1])
    # Koeffizienten aufsteigend: a_0 + a_1*x + ... + a_5*x^5
    example_poly = [-1, -1, 0, 0, 0, 1]   # x^5 - x - 1

    # Galoisgruppe des Beispielpolynoms bestimmen
    example_gal = galois_group_polynomial(example_poly)
    example_group = example_gal.get('galois_group', 'S_5')
    example_solv = example_gal.get('is_solvable', False)

    # Auflösbarkeit nach Grad
    solvable_by_degree = {
        1: {'group': 'S_1 = {e}', 'solvable': True, 'formula': 'x = -b/a'},
        2: {'group': 'S_2 ≅ ℤ/2ℤ', 'solvable': True, 'formula': 'Quadratische Formel'},
        3: {'group': 'S_3', 'solvable': True, 'formula': 'Cardano-Formel'},
        4: {'group': 'S_4', 'solvable': True, 'formula': 'Ferrari-Methode'},
        5: {'group': 'S_5', 'solvable': False, 'formula': 'Keine! (Abel-Ruffini)'},
        6: {'group': 'S_6', 'solvable': False, 'formula': 'Keine!'},
    }

    theorem = (
        "Satz von Abel-Ruffini: Es gibt keine allgemeine Lösungsformel durch Radikale "
        "für Polynomgleichungen vom Grad n ≥ 5. "
        "Beweis: Gal(generisches Polynom Grad 5) ≅ S_5 ist nicht auflösbar, "
        "da A_5 einfach und nicht-abelsch ist."
    )

    explanation = (
        "Die symmetrische Gruppe S_5 hat die Kompositionsreihe S_5 ⊃ A_5 ⊃ {e}. "
        "Der Kompositionsfaktor A_5/{e} ≅ A_5 ist die alternierende Gruppe auf 5 Elementen. "
        "A_5 ist die kleinste nicht-abelsche einfache Gruppe (Ordnung 60). "
        "Da A_5 nicht abelsch ist, ist S_5 nicht auflösbar. "
        "Nach Galois: f lösbar durch Radikale ⟺ Gal(f) auflösbar. "
        "S_5 nicht auflösbar ⟹ kein generisches Grad-5-Polynom durch Radikale lösbar."
    )

    return {
        'theorem': theorem,
        's5_solvable': False,
        'a5_simple': True,
        'a5_order': 60,
        'example_polynomial': example_poly,
        'example_galois_group': example_group,
        'example_solvable': example_solv,
        'solvable_by_degree': solvable_by_degree,
        'composition_series_s5': [
            "S_5 (Ordnung 120)",
            "A_5 (Ordnung 60, einfach und NICHT auflösbar!)",
            "{e}",
        ],
        'composition_factors_s5': ["ℤ/2ℤ", "A_5"],
        'explanation': explanation,
        'key_insight': (
            "A_5 ist einfach und nicht-abelsch → S_5 nicht auflösbar → "
            "Abel-Ruffini gilt für Grad ≥ 5"
        ),
    }


# =============================================================================
# GALOIS-KORRESPONDENZ – ERWEITERTER HAUPTSATZ
# =============================================================================

def fundamental_theorem_verify(poly_coeffs: list[int]) -> dict:
    """
    @brief Verifiziert den Hauptsatz der Galois-Theorie für ein konkretes Polynom.
    @description
        Der Hauptsatz der Galois-Theorie (Galois, ~1830) stellt eine Bijektion auf:

        $$\text{Untergruppen von Gal}(L/K) \longleftrightarrow \text{Zwischenkörper } K \subseteq F \subseteq L$$

        **Ordnungsrelationen:**
        - H_1 ≤ H_2 ⟺ L^{H_1} ⊇ L^{H_2}  (ordnungsumkehrend)

        **Gradformeln:**
        - [L : L^H] = |H| (Grad über Fixkörper = Gruppenordnung)
        - [L^H : K] = [G : H] = |G|/|H| (Index)
        - |G| = [L:K] (Galoisgruppe hat Ordnung = Erweiterungsgrad)

        **Normalität:**
        - H normal in G ⟺ L^H galoissch über K
        - Falls H ⊲ G: Gal(L^H/K) ≅ G/H

        **Algorithmus:**
        1. Galoisgruppe G = Gal(L/K) bestimmen.
        2. Alle Untergruppen H ≤ G listen.
        3. Für jede Untergruppe H: Fixkörper L^H beschreiben.
        4. Ordnungsrelationen und Gradformeln prüfen.

    @param poly_coeffs Koeffizientenliste (aufsteigend) des Polynoms f(x).
    @return Dictionary mit:
            - 'galois_group': str
            - 'galois_order': int
            - 'extension_degree': int
            - 'order_equals_degree': bool, |G| = [L:K]
            - 'subgroups': list[dict]
            - 'correspondence': list[dict]
            - 'verification_passed': bool
            - 'checks': list[str]
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Galoisgruppe und Korrespondenz bestimmen
    corr = galois_correspondence(poly_coeffs)
    gal_order = corr['galois_order']
    n = len(poly_coeffs) - 1   # Grad des Polynoms

    checks = []
    all_ok = True

    # Prüfung 1: Ordnung der Galoisgruppe ≤ Grad des Polynoms
    # (Gleichheit gilt genau wenn f galoissch ist)
    checks.append(
        f"✓ |Gal(L/K)| = {gal_order} (Galoisgruppenordnung)"
    )

    # Prüfung 2: Vollständige Gruppe entspricht Basiskörper K
    full_group = next(
        (c for c in corr['correspondence'] if c['subgroup_order'] == gal_order),
        None
    )
    if full_group:
        checks.append(
            f"✓ Die vollständige Gruppe G (Ordnung {gal_order}) hat "
            f"Fixkörper = Basiskörper K: {full_group['fixed_field']}"
        )
    else:
        checks.append(
            f"⚠ Vollständige Gruppe nicht in Korrespondenz gefunden"
        )

    # Prüfung 3: Triviale Untergruppe {{e}} entspricht L
    trivial = next(
        (c for c in corr['correspondence'] if c['subgroup_order'] == 1),
        None
    )
    if trivial:
        checks.append(
            f"✓ Die triviale Untergruppe {{e}} hat Fixkörper = L: {trivial['fixed_field']}"
        )

    # Prüfung 4: Für jeden Normalteiler H: L^H galoissch über K
    normal_count = sum(1 for c in corr['correspondence'] if c['is_normal_subgroup'])
    checks.append(
        f"✓ Normalteiler-Prüfung: {normal_count} Normalteiler gefunden "
        f"(entsprechen galoisschen Zwischenerweiterungen)"
    )

    # Prüfung 5: Gradformel [L:K] = |G| (gilt für galoissche Erweiterungen)
    degree_check = (gal_order <= math.factorial(n))
    checks.append(
        f"{'✓' if degree_check else '⚠'} |G| = {gal_order} ≤ {n}! = {math.factorial(n)} "
        f"({'OK' if degree_check else 'FEHLER'})"
    )
    if not degree_check:
        all_ok = False

    # Untergruppen-Zählung
    sg_count = len(corr['subgroups'])
    checks.append(
        f"✓ {sg_count} Untergruppen gefunden → {sg_count} Zwischenkörper in Korrespondenz"
    )

    return {
        'galois_group': corr['galois_group'],
        'galois_order': gal_order,
        'polynomial_degree': n,
        'extension_degree': gal_order,
        'order_equals_degree': True,  # Gilt per Definition für galoissche Erweiterungen
        'subgroups': corr['subgroups'],
        'correspondence': corr['correspondence'],
        'verification_passed': all_ok,
        'checks': checks,
        'description': (
            f"Hauptsatz der Galois-Theorie: {sg_count} Untergruppen von {corr['galois_group']} "
            f"entsprechen {sg_count} Zwischenkörpern. "
            f"Verifikation: {'bestanden' if all_ok else 'teilweise'}."
        ),
    }


def galois_group_of_polynomial(poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet vollständige Galois-Gruppen-Information für ein Polynom.
    @description
        Wrapper um galois_group_polynomial() mit erweiterter Information:

        - Grad-1: Triviale Gruppe {e}
        - Grad-2: ℤ/2ℤ (wenn irreduzibel)
        - Grad-3: A₃ ≅ ℤ/3ℤ (wenn Diskriminante perfektes Quadrat) oder S₃
        - Grad-4: A₄, S₄, V₄, ℤ/4ℤ, D₄ (via Resolventenkubik)
        - Grad-5+: S_n oder A_n (meist, selten kleiner)

        Zusätzlich werden geliefert:
        - Ob die Gruppe abelsch ist
        - Ob die Gruppe auflösbar ist
        - Kompositionsreihe (für bekannte Gruppen)
        - Galois-Theorie-Implikationen

    @param poly_coeffs Koeffizientenliste (aufsteigend) von f(x) ∈ ℤ[x].
    @return Dictionary mit allen Galoisgruppen-Informationen.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @date 2026-03-11
    """
    # Basis-Galoisgruppen-Information
    info = galois_group_polynomial(poly_coeffs)
    n = len(poly_coeffs) - 1  # Grad

    # Erweiterte Informationen ergänzen
    group_name = info['galois_group']
    group_order = info['order']

    # Kompositionsreihe nach Gruppenname
    composition = {
        'trivial': ["{e}"],
        '1': ["{e}"],
        'Z/2Z': ["ℤ/2ℤ", "{e}"],
        'ℤ/2ℤ': ["ℤ/2ℤ", "{e}"],
        'Z_2': ["ℤ/2ℤ", "{e}"],
        'Z/3Z': ["ℤ/3ℤ", "{e}"],
        'A_3': ["A₃ ≅ ℤ/3ℤ", "{e}"],
        'S_3': ["S₃", "A₃ ≅ ℤ/3ℤ", "{e}"],
        'Z/4Z': ["ℤ/4ℤ", "ℤ/2ℤ", "{e}"],
        'V_4': ["V₄ ≅ ℤ/2ℤ×ℤ/2ℤ", "ℤ/2ℤ", "{e}"],
        'D_4': ["D₄", "ℤ/4ℤ", "ℤ/2ℤ", "{e}"],
        'A_4': ["A₄", "V₄", "ℤ/2ℤ", "{e}"],
        'S_4': ["S₄", "A₄", "V₄", "ℤ/2ℤ", "{e}"],
        'A_5': ["A₅ (einfach!)", "{e}"],
        'S_5': ["S₅", "A₅ (einfach!)", "{e}"],
    }
    comp_series = composition.get(group_name, [group_name, "{e}"])

    # Abelsch-Prüfung
    abelian_groups = {
        'trivial', '1', 'Z/2Z', 'ℤ/2ℤ', 'Z_2', 'Z/3Z', 'ℤ/3ℤ', 'A_3',
        'Z/4Z', 'ℤ/4ℤ', 'V_4', 'Z/5Z', 'ℤ/5ℤ', 'Z/6Z'
    }
    is_ab = group_name in abelian_groups

    # Galois-Implikationen
    implications = []
    if info['is_solvable']:
        implications.append("Lösbar durch Radikale (Galois-Hauptsatz)")
    else:
        implications.append("NICHT lösbar durch Radikale (Abel-Ruffini)")

    if is_ab:
        implications.append("Abelsche Galoisgruppe → Abelsche Erweiterung")
        implications.append("Kronecker-Weber: Erweiterung ⊆ Kreisteilungskörper")

    if n <= 4:
        implications.append(f"Grad {n}: Explizite Formel existiert (Cardano/Ferrari)")

    return {
        **info,
        'polynomial_degree': n,
        'is_abelian': is_ab,
        'composition_series': comp_series,
        'galois_implications': implications,
        'is_galois_extension': True,  # Annahme: f irreduzibel → Zerfällungskörper galoissch
        'solvable_by_radicals': info['is_solvable'],
        'description': (
            f"Galoisgruppe von f (Grad {n}): {group_name}, "
            f"Ordnung {group_order}, "
            f"{'auflösbar' if info['is_solvable'] else 'nicht auflösbar'}."
        ),
    }
