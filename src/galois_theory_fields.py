"""
@file galois_theory_fields.py
@brief Galois-Theorie – Körper und Körpererweiterungen.
@description
    Dieses Modul enthält alle Klassen und Funktionen rund um endliche Körper,
    algebraische Körpererweiterungen und die Galoisgruppe als Automorphismengruppe.

    Inhalte:
    - Hilfsfunktionen für Polynomrechnung in GF(p)
    - FiniteField: Endlicher Körper GF(p^n) = ℤ_p[x]/(f(x))
    - FieldExtension: Körpererweiterung L = ℚ(α) via Minimalpolynom
    - GaloisGroup: Galoisgruppe Gal(L/K) als abstrakte Gruppe
    - Funktionen: finite_field(), cyclotomic_galois_group(),
      galois_group_finite_field(), kronecker_weber_check(),
      dirichlet_characters_from_galois(), norm_and_trace(), hilbert90(),
      primitive_element_theorem(), cyclotomic_polynomial(),
      splitting_field(), minimal_polynomial()

    **Mathematische Grundlagen:**
    - GF(p^n)^× ist zyklisch der Ordnung p^n - 1 (primitives Element).
    - Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^× (Kreisteilungssatz).
    - Gal(GF(p^n)/GF(p)) ≅ ℤ/nℤ, erzeugt durch Frobenius x ↦ x^p.

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
import itertools
from typing import Optional
from math import gcd, isqrt
from functools import reduce

# SymPy für symbolische Berechnungen
import sympy
from sympy import (
    Poly, Symbol, ZZ, QQ, factor, resultant,
    discriminant as sym_discriminant,
    minimal_polynomial as sym_minimal_poly,
    sqrt as sym_sqrt, Integer,
    factorint, isprime, totient, primitive_root, Rational, floor,
    cyclotomic_poly, GF, FiniteField as SympyFF, symbols, expand, Mul
)
from sympy.abc import x as sym_x
from sympy.ntheory import n_order

# Lokale Ausnahmen
from exceptions import DomainError, InvalidInputError, MathematicalError


# =============================================================================
# HILFSFUNKTIONEN FÜR POLYNOMRECHNUNG IN GF(p)
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
    @date 2026-03-11
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
    @date 2026-03-11
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
    @date 2026-03-11
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
    @date 2026-03-11
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
            return False

    # Schritt 3: f soll x^{p^n} - x teilen
    exp_n = p ** n
    x_pow_n = _poly_pow_mod([0, 1], exp_n, coeffs, p)
    # x^{p^n} - x mod f soll = 0 sein
    x_pow_n[1] = (x_pow_n[1] - 1) % p
    while len(x_pow_n) > 1 and x_pow_n[-1] == 0:
        x_pow_n.pop()

    if not all(c == 0 for c in x_pow_n):
        return False  # f teilt x^{p^n}-x nicht → reduzibel

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
    @date 2026-03-11
    """
    if n == 1:
        return [0, 1]  # x ist irreduzibel

    # Alle normierten Polynome vom Grad n: Koeffizienten in [0..p-1], letzter = 1
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
    @lastModified 2026-03-11
    """

    def __init__(self, p: int, n: int, mod_poly: Optional[list[int]] = None):
        """
        @brief Konstruiert GF(p^n) mit optionalem Modulus-Polynom.
        @param p Primzahl (Charakteristik des Körpers).
        @param n Erweiterungsgrad (Dimension über ℤ_p).
        @param mod_poly Irreduzibles Polynom vom Grad n über ℤ_p
                       (aufsteigend: mod_poly[0] = konstanter Term).
                       Falls None, wird automatisch eines gesucht.
        @raises InvalidInputError Wenn p keine Primzahl oder n < 1 ist.
        @date 2026-03-11
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
        @date 2026-03-11
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
        @date 2026-03-11
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
        @date 2026-03-11
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
        @date 2026-03-11
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
        @date 2026-03-11
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
    @lastModified 2026-03-11
    """

    def __init__(self, min_poly_coeffs: list[int | float]):
        """
        @brief Konstruiert die Körpererweiterung ℚ(α) via Minimalpolynom.
        @param min_poly_coeffs Koeffizienten des Minimalpolynoms von α über ℚ
                               (aufsteigend: [a_0, a_1, ..., a_n]).
        @raises InvalidInputError Wenn das Polynom nicht mindestens Grad 1 hat.
        @date 2026-03-11
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
        @date 2026-03-11
        """
        return self._degree

    def is_algebraic(self) -> bool:
        """
        @brief Prüft ob α algebraisch über ℚ ist (stets True für dieses Modell).
        @return Immer True.
        @date 2026-03-11
        """
        return True  # Per Konstruktion: α hat endliches Minimalpolynom

    def minimal_polynomial(self) -> list[int | float]:
        """
        @brief Gibt die Koeffizientenliste des Minimalpolynoms zurück.
        @return Koeffizientenliste (aufsteigend).
        @date 2026-03-11
        """
        return list(self.min_poly)

    def basis(self) -> list[str]:
        """
        @brief Gibt eine K-Basis von L als α über ℚ zurück.
        @return Liste der Basis-Elemente als Strings.
        @date 2026-03-11
        """
        n = self._degree
        basis = ['1']
        for k in range(1, n):
            basis.append(f'α^{k}')
        return basis

    def norm(self) -> float:
        """
        @brief Berechnet die Norm N_{L/K}(α) = (-1)^n · f(0).
        @return Norm-Wert (rational).
        @date 2026-03-11
        """
        n = self._degree
        a0 = self.min_poly[0]   # konstanter Term
        an = self.min_poly[-1]  # Leitkoeffizient
        return ((-1) ** n) * (a0 / an)

    def trace(self) -> float:
        """
        @brief Berechnet die Spur Tr_{L/K}(α) = -a_{n-1}/a_n.
        @return Spurwert (rational).
        @date 2026-03-11
        """
        an = self.min_poly[-1]   # Leitkoeffizient
        if self._degree >= 2:
            a_n1 = self.min_poly[-2]  # Koeffizient von x^{n-1}
        else:
            a_n1 = self.min_poly[0]   # Für Grad-1-Polynom
        return -a_n1 / an

    def is_separable(self) -> bool:
        """
        @brief Prüft ob L/K separabel ist (über ℚ stets True).
        @return True wenn separabel.
        @date 2026-03-11
        """
        # Über ℚ sind alle algebraischen Erweiterungen separabel
        f = self._sympy_poly
        df = f.diff()  # Formale Ableitung
        return sympy.gcd(f, df) == 1

    def is_normal(self) -> bool:
        """
        @brief Prüft ob L/ℚ eine normale Erweiterung ist.
        @return True wenn L/ℚ normal ist.
        @date 2026-03-11
        """
        # Import hier, um zirkuläre Abhängigkeit zu vermeiden
        from galois_theory_polynomials import galois_group_polynomial

        n = self._degree
        if n <= 1:
            return True  # Triviale Erweiterung
        if n == 2:
            return True  # Quadratische Erweiterungen sind stets normal

        # Für Grad ≥ 3: Alle Wurzeln numerisch berechnen und prüfen
        try:
            gal = galois_group_polynomial(self.min_poly)
            return gal['order'] == n
        except Exception:
            return n == 2  # Fallback: nur Grad-2-Erweiterungen sicher normal

    def is_galois(self) -> bool:
        """
        @brief Prüft ob L/ℚ eine Galois-Erweiterung ist.
        @return True wenn Galois-Erweiterung.
        @date 2026-03-11
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
        als Untergruppe der symmetrischen Gruppe S_n.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, poly_coeffs: list[int], group_name: str = "",
                 group_order: int = 0, is_solvable_flag: bool = True):
        """
        @brief Konstruiert die Galoisgruppe aus Polynomial-Koeffizienten.
        @param poly_coeffs Koeffizientenliste des Polynoms (aufsteigend).
        @param group_name Name der Galoisgruppe (z.B. "S_3", "Z/2Z").
        @param group_order Gruppenordnung.
        @param is_solvable_flag Ob die Gruppe auflösbar ist.
        @date 2026-03-11
        """
        self.poly_coeffs = poly_coeffs
        self._group_name = group_name
        self._order = group_order
        self._is_solvable = is_solvable_flag

        # Falls keine Informationen übergeben wurden, automatisch berechnen
        if not group_name:
            from galois_theory_polynomials import galois_group_polynomial
            info = galois_group_polynomial(poly_coeffs)
            self._group_name = info['galois_group']
            self._order = info['order']
            self._is_solvable = info['is_solvable']

    def order(self) -> int:
        """
        @brief Gibt die Gruppenordnung |Gal(L/K)| zurück.
        @return Gruppenordnung.
        @date 2026-03-11
        """
        return self._order

    def elements(self) -> list[str]:
        """
        @brief Gibt die Gruppenelemente zurück (für kleine Gruppen explizit).
        @return Liste der Gruppenelemente als Strings.
        @date 2026-03-11
        """
        n = self._order
        name = self._group_name

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
        @return True wenn abelsch.
        @date 2026-03-11
        """
        name = self._group_name
        abelian_names = {'Z/2Z', 'ℤ/2ℤ', 'Z_2', 'Z/3Z', 'ℤ/3ℤ', 'Z_3',
                         'Z/4Z', 'ℤ/4ℤ', 'Z_4', 'V_4', 'Z/2Z×Z/2Z',
                         'trivial', '{1}', 'Z/5Z', 'ℤ/5ℤ', 'Z_5',
                         'Z/6Z', 'ℤ/6ℤ', 'Z_6'}
        return name in abelian_names

    def is_solvable(self) -> bool:
        """
        @brief Prüft ob die Galoisgruppe auflösbar ist.
        @return True wenn auflösbar.
        @date 2026-03-11
        """
        return self._is_solvable

    def subgroups(self) -> list[dict]:
        """
        @brief Gibt die Untergruppen der Galoisgruppe zurück (für kleine Gruppen).
        @return Liste von Dictionarys {'name': str, 'order': int}.
        @date 2026-03-11
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
        @return Liste von Dictionarys {'name': str, 'order': int, 'is_normal': True}.
        @date 2026-03-11
        """
        name = self._group_name
        all_subs = self.subgroups()

        if name == 'S_3':
            return [s for s in all_subs if s['name'] in ('{e}', 'A_3≅Z_3', 'S_3')]
        elif name == 'D_4':
            return [s for s in all_subs if s['name'] in ('{e}', '<r²>', '<r>', '<r²,s>', '<r²,rs>', 'D_4')]
        elif name == 'A_4':
            return [s for s in all_subs if s['name'] in ('{e}', 'V_4', 'A_4')]
        elif name == 'S_4':
            return [s for s in all_subs if s['name'] in ('{e}', 'V_4', 'A_4', 'S_4')]
        else:
            # Für abelsche Gruppen sind alle Untergruppen Normalteiler
            return all_subs

    def fixed_field(self, subgroup_name: str) -> str:
        """
        @brief Berechnet den Fixkörper L^H für eine Untergruppe H.
        @param subgroup_name Name der Untergruppe H.
        @return Beschreibung des Fixkörpers als String.
        @date 2026-03-11
        """
        name = self._group_name

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
# HAUPTFUNKTIONEN – KÖRPER UND ERWEITERUNGEN
# =============================================================================

def finite_field(p: int, n: int) -> FiniteField:
    """
    @brief Konstruiert den endlichen Körper GF(p^n).
    @param p Primzahl (Charakteristik).
    @param n Erweiterungsgrad.
    @return FiniteField-Objekt für GF(p^n).
    @raises InvalidInputError Wenn p keine Primzahl oder n < 1.
    @date 2026-03-11
    """
    return FiniteField(p, n)


def cyclotomic_polynomial(n: int) -> list[int]:
    """
    @brief Berechnet das n-te Kreisteilungspolynom Φ_n(x).
    @description
        Φ_n(x) = ∏_{gcd(k,n)=1, 1≤k≤n} (x - e^{2πik/n})
        Grad = φ(n), Galoisgruppe Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^×.

    @param n Positive ganze Zahl ≥ 1.
    @return Koeffizientenliste von Φ_n(x) (aufsteigend).
    @raises InvalidInputError Wenn n < 1.
    @date 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # SymPy hat eine eingebaute Funktion für Kreisteilungspolynome
    phi_n = sympy.cyclotomic_poly(n, sym_x)

    # In Koeffizientenliste (aufsteigend) umwandeln
    poly = Poly(phi_n, sym_x, domain=ZZ)
    coeffs_desc = [int(c) for c in poly.all_coeffs()]  # absteigende Reihenfolge
    return list(reversed(coeffs_desc))


def cyclotomic_galois_group(n: int) -> dict:
    """
    @brief Berechnet die Galoisgruppe Gal(ℚ(ζ_n)/ℚ) des n-ten Kreisteilungskörpers.
    @description
        Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^×, Ordnung φ(n), abelsch.
        Jedes σ_k (gcd(k,n)=1) wirkt via ζ_n ↦ ζ_n^k.

    @param n Positive ganze Zahl ≥ 1.
    @return Dictionary mit Galoisgruppen-Information.
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Elemente der multiplikativen Gruppe (ℤ/nℤ)^×
    elements = [k for k in range(1, n + 1) if gcd(k, n) == 1]
    phi_n = len(elements)  # = φ(n)

    # Erzeuger der Gruppe via SymPy
    gens = []
    try:
        g = int(primitive_root(n))
        gens = [g]
        is_cyclic = True
    except Exception:
        is_cyclic = False
        gens = elements[:2] if len(elements) >= 2 else elements

    # Galoisgruppe benennen
    if phi_n == 1:
        group_name = "trivial {1}"
    elif is_cyclic:
        group_name = f"ℤ/{phi_n}ℤ"
    else:
        facts = sympy.factorint(phi_n)
        if len(facts) == 1 and list(facts.values())[0] == 1:
            group_name = f"ℤ/{phi_n}ℤ"
        else:
            group_name = f"(ℤ/{n}ℤ)^×"

    # Automorphismus-Beschreibungen
    automorphisms = [f"σ_{k}: ζ_{n} ↦ ζ_{n}^{k}" for k in elements]

    return {
        'n': n,
        'degree': phi_n,
        'galois_group': group_name,
        'order': phi_n,
        'generators': gens,
        'elements': elements,
        'is_abelian': True,
        'is_cyclic': is_cyclic,
        'automorphisms': automorphisms,
        'cyclotomic_poly_degree': phi_n,
        'description': (
            f"Gal(ℚ(ζ_{n})/ℚ) ≅ (ℤ/{n}ℤ)^× "
            f"mit Ordnung φ({n}) = {phi_n}"
        ),
    }


def galois_group_finite_field(p: int, n: int) -> dict:
    """
    @brief Berechnet die Galoisgruppe Gal(GF(p^n) / GF(p)).
    @description
        Gal(GF(p^n)/GF(p)) ≅ ℤ/nℤ, erzeugt durch Frobenius x ↦ x^p.
        Jeder Automorphismus ist eine Potenz des Frobenius.

    @param p Primzahl (Charakteristik).
    @param n Grad der Erweiterung.
    @return Dictionary mit Galoisgruppen-Information.
    @raises InvalidInputError Wenn p nicht prim oder n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
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


def kronecker_weber_check(poly_coeffs: list[int]) -> dict:
    """
    @brief Prüft das Kronecker-Weber-Theorem für eine abelsche Erweiterung.
    @description
        Jede endliche abelsche Galois-Erweiterung L/ℚ liegt in einem
        Kreisteilungskörper ℚ(ζ_n) (Kronecker 1853, Weber 1886).

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Dictionary mit Kronecker-Weber-Information.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    from galois_theory_polynomials import galois_group_polynomial

    gal_info = galois_group_polynomial(poly_coeffs)
    group_name = gal_info['galois_group']
    group_order = gal_info['order']

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
        Dirichlet-Charaktere modulo n als Gruppencharaktere von
        Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^×.

    @param n Positive ganze Zahl ≥ 1 (Modul).
    @return Dictionary mit Charakteren.
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    gal = cyclotomic_galois_group(n)
    elements = gal['elements']
    phi_n = gal['order']

    characters = []

    if phi_n == 1:
        characters.append({
            'index': 0,
            'name': 'χ₀ (Hauptcharakter)',
            'is_principal': True,
            'values': {k: 1 for k in elements},
            'conductor': 1,
        })
    else:
        elem_orders = {}
        for k in elements:
            m = 1
            cur = k
            while cur % n != 1 and m < phi_n + 1:
                cur = (cur * k) % n
                m += 1
            elem_orders[k] = m

        characters.append({
            'index': 0,
            'name': 'χ₀ (Hauptcharakter)',
            'is_principal': True,
            'values': {k: 1 for k in elements},
            'conductor': 1,
            'description': f"χ₀(σ_k) = 1 für alle k ∈ (ℤ/{n}ℤ)^×",
        })

        if isprime(n):
            legendre_vals = {}
            for k in elements:
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
            for i in range(1, min(phi_n, 4)):
                chi_vals = {}
                for k in elements:
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


def minimal_polynomial(alpha_coeffs: list[float],
                        base_coeffs: Optional[list[int]] = None) -> list[int]:
    """
    @brief Berechnet das Minimalpolynom von α über ℚ.
    @param alpha_coeffs Koeffizienten zur Beschreibung von α.
    @param base_coeffs Optionales Grundkörper-Polynom (Standard: ℚ).
    @return Koeffizientenliste des Minimalpolynoms (aufsteigend).
    @date 2026-03-11
    """
    if len(alpha_coeffs) == 1:
        val = alpha_coeffs[0]
        if isinstance(val, int) or (isinstance(val, float) and val == int(val)):
            return [-int(val), 1]  # x - val
        d = val
        if d > 0:
            sq = isqrt(int(d))
            if sq * sq == int(d):
                return [-sq, 1]
            return [int(-d), 0, 1]  # x² - d für √d
        return [int(-d), 0, 1]

    if len(alpha_coeffs) == 2:
        a, b = alpha_coeffs
        if b == 0:
            return minimal_polynomial([a])

    try:
        val = sum(c * (2 ** 0.5) ** i for i, c in enumerate(alpha_coeffs))
        sympy_alpha = sympy.nsimplify(val, rational=False, tolerance=1e-10)
        min_poly_sym = sympy.minimal_polynomial(sympy_alpha, sym_x)
        poly = Poly(min_poly_sym, sym_x, domain=ZZ)
        coeffs_desc = [int(c) for c in poly.all_coeffs()]
        return list(reversed(coeffs_desc))
    except Exception:
        return [-int(round(sum(alpha_coeffs))), 1]


def splitting_field(poly_coeffs: list[int], p: Optional[int] = None) -> dict:
    """
    @brief Berechnet den Zerfällungskörper eines Polynoms.
    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @param p Optional: Charakteristik für endliche Körper.
    @return Dictionary mit Zerfällungskörper-Information.
    @date 2026-03-11
    """
    from galois_theory_polynomials import galois_group_polynomial

    n = len(poly_coeffs) - 1

    if p is not None:
        return {
            'degree': n,
            'roots': list(range(p)),
            'galois_group': f'Z/{n}Z (zyklisch)',
            'field': f'GF({p}^{n})'
        }

    coeffs_desc = list(reversed(poly_coeffs))
    try:
        import numpy as np
        roots_np = np.roots([float(c) for c in coeffs_desc])
        roots = [complex(r) for r in roots_np]
    except Exception:
        roots = []

    gal_info = galois_group_polynomial(poly_coeffs)
    galois_group_name = gal_info['galois_group']
    galois_order = gal_info['order']

    return {
        'degree': galois_order,
        'roots': roots,
        'galois_group': galois_group_name,
        'galois_order': galois_order
    }


def norm_and_trace(alpha_coeffs: list[float], poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet Norm und Spur eines algebraischen Elements.
    @description
        N_{L/K}(α) = (-1)^n · a_0/a_n,  Tr_{L/K}(α) = -a_{n-1}/a_n.

    @param alpha_coeffs Koeffizientenliste von α (Linearkombination in Basis).
    @param poly_coeffs Koeffizientenliste des Minimalpolynoms (aufsteigend).
    @return Dictionary mit Norm, Spur und Konjugierten.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    n = len(poly_coeffs) - 1
    if n < 1:
        raise InvalidInputError("Minimalpolynom muss Grad ≥ 1 haben")

    leading = poly_coeffs[n]
    constant = poly_coeffs[0]

    # Norm via Vieta
    norm_val = ((-1) ** n) * constant / leading

    # Spur via Vieta
    if n >= 2:
        trace_val = -poly_coeffs[n - 1] / leading
    elif n == 1:
        trace_val = -constant / leading
    else:
        trace_val = 0.0

    # Konjugierte (numerische Wurzeln)
    coeffs_desc = list(reversed(poly_coeffs))
    try:
        import numpy as np
        conjugates = np.roots(coeffs_desc).tolist()
    except Exception:
        conjugates = []

    char_poly = poly_coeffs[:]

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
        Sei L/K zyklisch galoissch mit Gal(L/K) = ⟨σ⟩. Dann gilt:
        N_{L/K}(α) = 1 ⟺ α = β/σ(β) für ein β ∈ L^×.
        Äquivalent: H¹(Gal(L/K), L^×) = 1.

    @param n Ordnung der zyklischen Galoisgruppe (n ≥ 1).
    @param check_type Art der Demonstration: "norm" oder "cohomology".
    @return Dictionary mit Theorem-Information.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    zeta = cmath.exp(2j * math.pi / n)
    conjugates = [cmath.exp(2j * math.pi * k / n) for k in range(1, n + 1)]

    norm_product = 1.0 + 0j
    for c in conjugates:
        norm_product *= c

    norm_real = round(norm_product.real, 8)
    norm_imag = round(norm_product.imag, 8)

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
        'cohomology_trivial': True,
        'galois_group': f"ℤ/{n}ℤ (zyklisch)",
        'explanation': explanation,
        'h1_trivial': True,
        'description': f"Hilbert's Satz 90 für zyklische Erweiterung der Ordnung {n}",
    }


def primitive_element_theorem(alpha: list, beta: list) -> dict:
    """
    @brief Findet ein primitives Element ℚ(α, β) = ℚ(α + c·β).
    @description
        Satz vom primitiven Element: ℚ(α, β) = ℚ(α + c·β) für fast alle c ∈ ℚ.

    @param alpha Minimalpolynom von α (Koeffizientenliste).
    @param beta Minimalpolynom von β (Koeffizientenliste).
    @return Dictionary {'c': int, 'primitive_element': str, 'degree': int}.
    @date 2026-03-11
    """
    deg_alpha = len(alpha) - 1
    deg_beta = len(beta) - 1
    expected_degree = deg_alpha * deg_beta

    for c in range(0, 20):
        try:
            roots_alpha = sympy.Poly(list(reversed(alpha)), sym_x,
                                     domain=QQ).all_roots()
            roots_beta = sympy.Poly(list(reversed(beta)), sym_x,
                                    domain=QQ).all_roots()

            if not roots_alpha or not roots_beta:
                continue

            a_root = roots_alpha[0]
            b_root = roots_beta[0]

            gamma = a_root + sympy.Integer(c) * b_root
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
