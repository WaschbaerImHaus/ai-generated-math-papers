"""
Mahler-Maß von Polynomen – Jensen-Formel, Lehmer's Problem, Smyth-Schranke.

Das Mahler-Maß M(P) eines Polynoms P misst die „arithmetische Größe" eines
Polynoms über die geometrische Mittelung seiner Werte auf dem Einheitskreis.

**Definition (Jensen-Formel):**
    M(P) = exp(∫₀¹ log|P(e^{2πit})| dt) = |aₙ| · ∏_{|αᵢ|>1} |αᵢ|

wobei aₙ der Leitkoeffizient und αᵢ die Wurzeln von P sind.

**Lehrers Problem (1933, ungelöst):**
    Gibt es ein Polynom P ∈ ℤ[x], P nicht Kreisteilungspolynom,
    mit 1 < M(P) < M(x¹⁰+x⁹-x⁷-x⁶-x⁵-x⁴-x³+x+1) ≈ 1.17628?

**Smyth (1971):**
    Für nicht-reziproke Polynome (P(x) ≠ ±xⁿP(1/x)) gilt M(P) ≥ 1.3247...
    (kleinste Pisot-Zahl, Wurzel von x³-x-1)

**Kronecker-Theorem:**
    M(P) = 1 ↔ alle Wurzeln von P liegen auf dem Einheitskreis oder sind 0.
    Äquivalent: P ist ein Produkt von Kreisteilungspolynomen und xᵏ.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import numpy as np
import sympy as sp
from sympy import Poly, Symbol, roots, Abs, log, exp, pi, I, N
from sympy import cyclotomic_poly, factorint, ZZ
import mpmath
from typing import Union, List, Optional, Tuple
from fractions import Fraction


# ===========================================================================
# MAHLER-MAß KLASSE
# ===========================================================================

class MahlerMeasure:
    """
    Berechnung des Mahler-Maßes für Polynome über ℤ oder ℚ.

    Das Mahler-Maß ist ein fundamentales Maß für die arithmetische Komplexität
    von Polynomen und verknüpft Polynom-Arithmetik mit Transzendenztheorie.

    **Wichtige Formeln:**
        - Produktformel:  M(P) = |aₙ| · ∏_{|αᵢ|>1} |αᵢ|
        - Jensen-Formel:  log M(P) = ∫₀¹ log|P(e^{2πit})| dt
        - Logarithm. Mahler-Maß: m(P) = log M(P)

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Lehmer's Dezim-Polynom: das Polynom mit kleinstem bekanntem Mahler-Maß > 1
    LEHMER_COEFFS = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
    # M(Lehmer) = 1.17628081825991...
    LEHMER_MAHLER = 1.1762808182599175

    # Smyth-Konstante: kleinste Pisot-Zahl (Wurzel von x³-x-1)
    SMYTH_CONSTANT = 1.3247179572447460259609088

    def __init__(self, coeffs: List[Union[int, float, Fraction]], precision: int = 50):
        """
        Initialisiert MahlerMeasure mit Polynomkoeffizienten.

        Die Koeffizienten werden von niedrigstem zum höchsten Grad angegeben:
            coeffs = [a₀, a₁, ..., aₙ]  →  P(x) = a₀ + a₁x + ... + aₙxⁿ

        @param coeffs: Liste der Koeffizienten [a₀, a₁, ..., aₙ]
        @param precision: Dezimalstellen für mpmath-Berechnungen (Standard: 50)
        @raises ValueError: wenn coeffs leer oder nur Nullen enthält
        """
        if not coeffs or all(c == 0 for c in coeffs):
            raise ValueError("Polynomial darf nicht das Nullpolynom sein.")

        # Entferne führende Nullen (höchste Koeffizienten)
        while len(coeffs) > 1 and coeffs[-1] == 0:
            coeffs = coeffs[:-1]

        self.coeffs = list(coeffs)
        self.degree = len(coeffs) - 1
        self.precision = precision
        mpmath.mp.dps = precision

        # sympy-Polynom (Koeffizienten von höchstem zu niedrigstem Grad)
        x = Symbol('x')
        self._poly_sympy = Poly(list(reversed(coeffs)), x, domain='QQ')

    @property
    def leading_coeff(self) -> float:
        """Leitkoeffizient aₙ des Polynoms."""
        return float(self.coeffs[-1])

    def get_roots(self) -> List[complex]:
        """
        Berechnet alle Wurzeln des Polynoms numerisch.

        Nutzt numpy für numerische Stabilität bei höheren Graden.

        @return: Liste aller Wurzeln als komplexe Zahlen
        """
        # numpy polyroots erwartet Koeffizienten von höchstem zu niedrigstem Grad
        coeffs_high_to_low = list(reversed(self.coeffs))
        return np.roots(coeffs_high_to_low).tolist()

    def compute_product_formula(self) -> float:
        """
        Berechnet M(P) = |aₙ| · ∏_{|αᵢ|>1} |αᵢ| (Produktformel).

        Dies ist die direkteste Berechnungsmethode: Multipliziere den
        Absolutwert des Leitkoeffizienten mit allen Wurzeln außerhalb
        des Einheitskreises.

        @return: Mahler-Maß M(P) als float
        """
        roots_list = self.get_roots()
        result = abs(self.leading_coeff)

        for r in roots_list:
            modulus = abs(r)
            if modulus > 1.0 + 1e-10:  # Toleranz für numerische Fehler
                result *= modulus

        return result

    def compute_jensen_integral(self, num_points: int = 10000) -> float:
        """
        Berechnet M(P) via Jensen-Formel: exp(∫₀¹ log|P(e^{2πit})| dt).

        Numerische Integration mittels Quadratur auf dem Einheitskreis.
        Diese Methode ist nützlich zur Verifikation und für Polynome,
        deren Wurzeln schwer zu bestimmen sind.

        @param num_points: Anzahl der Stützstellen für die numerische Integration
        @return: Mahler-Maß M(P) als float
        """
        # Stützstellen t ∈ [0, 1)
        t_values = np.linspace(0, 1, num_points, endpoint=False)
        # Punkte auf dem Einheitskreis: e^{2πit}
        z_values = np.exp(2j * np.pi * t_values)

        # Werte P(e^{2πit}) berechnen (Horner-Schema über numpy)
        coeffs_high = list(reversed(self.coeffs))
        poly_values = np.polyval(coeffs_high, z_values)

        # log|P(e^{2πit})| mit Behandlung von Nullstellen auf dem Kreis
        log_abs_values = np.log(np.abs(poly_values) + 1e-300)

        # Numerische Integration (Rechteckregel = Riemann-Summe)
        integral = np.mean(log_abs_values)

        return float(np.exp(integral))

    def logarithmic_mahler_measure(self) -> float:
        """
        Berechnet das logarithmische Mahler-Maß m(P) = log M(P).

        Das logarithmische Mahler-Maß addiert sich unter Produkten:
            m(P·Q) = m(P) + m(Q)

        @return: Logarithmisches Mahler-Maß m(P) = log M(P)
        """
        return np.log(self.compute_product_formula())

    def is_kronecker(self, tol: float = 1e-6) -> bool:
        """
        Prüft ob M(P) = 1 gilt (Kronecker-Theorem).

        Nach Kronecker gilt M(P) = 1 genau dann, wenn alle Wurzeln von P
        entweder auf dem Einheitskreis liegen oder 0 sind.
        Äquivalent: P ist ein Produkt aus Kreisteilungspolynomen und xᵏ.

        @param tol: Toleranz für den Vergleich mit 1
        @return: True wenn M(P) = 1 (im Rahmen der Toleranz)
        """
        mahler = self.compute_product_formula()
        return abs(mahler - 1.0) < tol

    def is_reciprocal(self) -> bool:
        """
        Prüft ob P reziprokes Polynom ist: P(x) = ±xⁿ · P(1/x).

        Reziproke Polynome erfüllen: a_k = ±a_{n-k} für alle k.
        Das Lehmer-Polynom ist reziprok. Smyth's Schranke gilt für
        NICHT-reziproke Polynome.

        @return: True wenn P reziprok ist
        """
        n = self.degree
        coeffs = self.coeffs

        # Prüfe a_k = a_{n-k} (palindromisch) oder a_k = -a_{n-k} (antipalindromisch)
        palindrome = all(
            abs(coeffs[k] - coeffs[n - k]) < 1e-10
            for k in range(n + 1)
        )
        anti_palindrome = all(
            abs(coeffs[k] + coeffs[n - k]) < 1e-10
            for k in range(n + 1)
        )
        return palindrome or anti_palindrome

    def smyth_bound_applies(self) -> bool:
        """
        Prüft ob Smyth's Schranke M(P) ≥ 1.3247... anwendbar ist.

        Smyth (1971) zeigte: Für ganzzahlige, nicht-reziproke Polynome
        gilt M(P) ≥ 1.3247... (die kleinste Pisot-Zahl).

        @return: True wenn P nicht reziprok ist (Smyth-Schranke anwendbar)
        """
        return not self.is_reciprocal()

    @classmethod
    def lehmer_polynomial(cls) -> 'MahlerMeasure':
        """
        Erstellt das Lehmer-Dezim-Polynom.

        Lehmer's Polynom (1933):
            L(x) = x¹⁰ + x⁹ - x⁷ - x⁶ - x⁵ - x⁴ - x³ + x + 1

        Dies ist das Polynom mit dem kleinsten bekannten Mahler-Maß > 1
        unter ganzzahligen Polynomen. Es ist reziprok.

        **Lehrers Problem (Conjecture):**
            Es gibt kein ganzzahliges Polynom P mit 1 < M(P) < M(L) ≈ 1.17628.
            (Status: OFFEN / Conjecture, 1933)

        @return: MahlerMeasure-Instanz für Lehmer's Polynom
        """
        # Koeffizienten von niedrigstem zu höchstem Grad:
        # x¹⁰+x⁹+0·x⁸-x⁷-x⁶-x⁵-x⁴-x³+0·x²+x+1
        return cls(cls.LEHMER_COEFFS)

    @classmethod
    def smyth_polynomial(cls) -> 'MahlerMeasure':
        """
        Erstellt das Smyth-Polynom x³ - x - 1.

        Dies ist das nicht-reziproke Polynom mit kleinstem Mahler-Maß.
        Sein Mahler-Maß ist die kleinste Pisot-Zahl: 1.3247...

        @return: MahlerMeasure-Instanz für Smyth's Polynom
        """
        # x³ - x - 1: Koeffizienten [a₀, a₁, a₂, a₃] = [-1, -1, 0, 1]
        return cls([-1, -1, 0, 1])

    @classmethod
    def cyclotomic(cls, n: int) -> 'MahlerMeasure':
        """
        Erstellt das n-te Kreisteilungspolynom Φₙ(x).

        Kreisteilungspolynome haben M = 1 (alle Wurzeln auf dem Einheitskreis).

        @param n: Index des Kreisteilungspolynoms (n ≥ 1)
        @return: MahlerMeasure-Instanz für Φₙ(x)
        """
        x = Symbol('x')
        cyc_poly = cyclotomic_poly(n, x)
        p = Poly(cyc_poly, x, domain='ZZ')
        # Koeffizienten von höchstem zu niedrigstem Grad → umkehren für interne Darstellung
        coeffs_high = [int(c) for c in p.all_coeffs()]
        return cls(list(reversed(coeffs_high)))

    def salem_number_search(self, degree: int = 10, bound: float = 1.3) -> List[float]:
        """
        Sucht Salem-Zahlen (Mahler-Maß nahe 1) unter reziproken Polynomen.

        Salem-Zahlen sind algebraische Ganzzahlen α > 1, bei denen alle
        konjugierten außer α und 1/α auf dem Einheitskreis liegen.
        Sie sind kandidaten für Lehrers Problem.

        **Note:** Dies ist eine Heuristik, keine vollständige Suche.

        @param degree: Grad der zu untersuchenden Polynome
        @param bound: Obere Schranke für M(P) zur Klassifikation als Salem
        @return: Liste von Mahler-Maßen nahe 1 (mögliche Salem-Zahlen)
        """
        salem_candidates = []

        # Untersuche einige palindromische Polynome mit kleinen Koeffizienten
        import itertools
        half = degree // 2

        for mid_coeffs in itertools.product([-1, 0, 1], repeat=half):
            # Baue palindromisches Polynom: [c₀, c₁, ..., c_h, ..., c₁, c₀]
            full = list(mid_coeffs) + list(reversed(mid_coeffs))
            if degree % 2 == 0:
                full = list(mid_coeffs) + list(reversed(mid_coeffs))
            if full[-1] == 0:
                continue

            try:
                mm = cls(full)
                m = mm.compute_product_formula()
                if 1.0 < m < bound:
                    salem_candidates.append(m)
            except Exception:
                continue

        return sorted(set(round(s, 8) for s in salem_candidates))

    def __repr__(self) -> str:
        """String-Darstellung des Polynoms mit Mahler-Maß."""
        m = self.compute_product_formula()
        return f"MahlerMeasure(deg={self.degree}, M={m:.6f})"


# ===========================================================================
# SCHUR-SIEGEL-SMYTH SPURPFAD (HILFSFUNKTIONEN)
# ===========================================================================

class SchurSiegelSmythTrace:
    """
    Schur-Siegel-Smyth-Pfad: Analyse von Spuren total positiver algebraischer Zahlen.

    Der Spurpfad untersucht die Minimalspuren total positiver algebraischer
    Ganzzahlen und ist eng mit Lehmer's Problem verknüpft.

    **Hauptresultat (Schur 1918, Siegel 1945):**
        Für eine total positive algebraische Ganzzahl α vom Grad d gilt:
            Spur(α)/d ≥ 1.7719...  (Schur-Siegel-Smyth Konstante)

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Schur-Siegel-Smyth Konstante (untere Schranke für Spur/Grad)
    SSS_CONSTANT = 1.7719  # Approximation

    def __init__(self):
        """Initialisiert den Spurpfad-Analysator."""
        pass

    def trace_ratio(self, poly_coeffs: List[int]) -> Optional[float]:
        """
        Berechnet das Verhältnis Spur(α)/Grad für ein gegebenes Polynom.

        @param poly_coeffs: Koeffizienten [a₀, ..., aₙ] des Minimalpolynoms
        @return: Spur(α)/d oder None bei Berechnungsfehler
        """
        if not poly_coeffs or poly_coeffs[-1] == 0:
            return None

        coeffs_high = list(reversed(poly_coeffs))
        try:
            roots_list = np.roots(coeffs_high)
            degree = len(poly_coeffs) - 1
            trace = sum(r.real for r in roots_list)  # Spur = Summe der Wurzeln
            return trace / degree
        except Exception:
            return None
