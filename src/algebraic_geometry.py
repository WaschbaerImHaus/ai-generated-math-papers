"""
@file algebraic_geometry.py
@brief Algebraische Geometrie: Varietäten, Nullstellensatz, elliptische Kurven, Schnitttheorie.
@description
    Implementiert grundlegende Konzepte der algebraischen Geometrie:

    - AffineVariety: Affine algebraische Varietät, definiert durch ein Polynomideal.
      Krull-Dimension, Irreduzibilität, singuläre Punkte (Jacobi-Kriterium), Schnitt.

    - ProjectiveVariety: Projektive Varietät, definiert durch homogene Polynome.
      Grad, Glattheit.

    - HilbertBasisTheorem: Demonstration des Hilbert-Basissatzes via Gröbner-Basen.
      Gröbner-Basis-Berechnung, Ideal-Mitgliedschaft.

    - NullstellensatzDemo: Hilbert'scher Nullstellensatz.
      Schwacher Nullstellensatz (Existenz gemeinsamer Nullstellen),
      Starker Nullstellensatz (Radikalmitgliedschaft).

    - EllipticCurveGeometry: Elliptische Kurven y² = x³ + ax + b.
      Diskriminante, j-Invariante, Glattheit, rationale Punkte mod p.

    - IntersectionTheory: Schnitttheorie nach Bézout, lokale Schnittmultiplizität.

    - MorphismOfVarieties: Morphismus algebraischer Varietäten.
      Dominanz, Pullback einer Funktion.

    Standalone-Funktionen:
    - compute_groebner_basis
    - ideal_radical
    - zariski_closure_demo
    - hilbert_polynomial

    Mathematischer Hintergrund:
    Die algebraische Geometrie untersucht Lösungsmengen von Polynomgleichungen
    (Varietäten) mit Methoden der kommutativen Algebra. Fundamental sind:
    - Hilbert-Basissatz: Jedes Ideal in k[x₁,...,xₙ] ist endlich erzeugt.
    - Nullstellensatz: Verknüpft Ideale mit geometrischen Nullstellenmengen.
    - Gröbner-Basen: Algorithmisches Werkzeug zur Ideal-Arithmetik.

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import sympy
from sympy import (
    symbols, poly, groebner, factor, expand, diff, simplify,
    Symbol, Poly, Rational, Integer, oo, S, sqrt
)
from sympy.polys import ZZ, QQ
from sympy import groebner as sympy_groebner
from sympy.abc import x, y, z
from typing import Optional
import math


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _to_sympy_poly(expr, variables: list):
    """
    @brief Konvertiert einen Ausdruck in ein SymPy-Polynom.
    @param expr  SymPy-Ausdruck oder Zahl
    @param variables Liste der SymPy-Symbole
    @return SymPy-Ausdruck (expand)
    @lastModified 2026-03-10
    """
    return expand(expr)


def _jacobian_matrix(polys: list, variables: list) -> list:
    """
    @brief Berechnet die Jacobi-Matrix ∂fᵢ/∂xⱼ.
    @description
        Die Jacobi-Matrix J einer Menge von Polynomen f₁,...,fₘ in
        Variablen x₁,...,xₙ ist die Matrix mit J[i][j] = ∂fᵢ/∂xⱼ.
    @param polys     Liste von SymPy-Ausdrücken (Polynome)
    @param variables Liste von SymPy-Symbolen
    @return Liste von Listen: jacobian[i][j] = ∂polys[i]/∂variables[j]
    @lastModified 2026-03-10
    """
    return [[diff(f, v) for v in variables] for f in polys]


def _eval_at_point(expr, variables: list, point: dict):
    """
    @brief Wertet einen SymPy-Ausdruck an einem gegebenen Punkt aus.
    @param expr      SymPy-Ausdruck
    @param variables Liste von SymPy-Symbolen
    @param point     Dictionary {Symbol: Wert}
    @return Numerischer Wert
    @lastModified 2026-03-10
    """
    return expr.subs(point)


# ===========================================================================
# KLASSE: AffineVariety
# ===========================================================================

class AffineVariety:
    """
    Affine algebraische Varietät, definiert durch ein Polynomideal.

    @description
        Eine affine Varietät V(f₁,...,fₘ) ⊂ 𝔸ⁿ ist die Menge der gemeinsamen
        Nullstellen von Polynomen f₁,...,fₘ ∈ k[x₁,...,xₙ]:
            V(f₁,...,fₘ) = { p ∈ 𝔸ⁿ : f₁(p) = ... = fₘ(p) = 0 }

        Das erzeugte Ideal ist I = ⟨f₁,...,fₘ⟩ ⊂ k[x₁,...,xₙ].

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, generators: list, variables: list):
        """
        @brief Initialisiert die affine Varietät.
        @description
            generators: Liste der Erzeuger-Polynome des Ideals.
            variables:  Liste der SymPy-Symbole (Koordinaten des affinen Raums).
        @param generators Liste von SymPy-Ausdrücken f₁,...,fₘ
        @param variables  Liste von SymPy-Symbolen x₁,...,xₙ
        @raises ValueError bei leerem variables-Liste
        @lastModified 2026-03-10
        """
        if not variables:
            raise ValueError("Variablenliste darf nicht leer sein.")
        # Erzeuger normalisieren (expand für konsistente Form)
        self.generators = [expand(g) for g in generators]
        self.variables = list(variables)
        # Interne Gröbner-Basis (lazy berechnet)
        self._groebner = None

    def _get_groebner(self):
        """
        @brief Berechnet (lazy) die Gröbner-Basis des Ideals.
        @return SymPy-Gröbner-Basis-Objekt
        @lastModified 2026-03-10
        """
        if self._groebner is None:
            if not self.generators or all(g == 0 for g in self.generators):
                # Nullideal: Gröbner-Basis ist [0]
                self._groebner = [sympy.Integer(0)]
            else:
                try:
                    gb = sympy_groebner(self.generators, *self.variables, order='lex')
                    self._groebner = list(gb)
                except Exception:
                    self._groebner = self.generators
        return self._groebner

    def dimension(self) -> int:
        """
        @brief Berechnet die Krull-Dimension der Varietät.
        @description
            Die Krull-Dimension von V(I) entspricht der Dimension des
            Koordinatenrings k[x₁,...,xₙ]/I als k-Algebra.

            Vereinfachte Berechnung über die Gröbner-Basis:
            - Vollständiges Ideal (1 ∈ I) → leere Varietät, Dim = -1
            - Nullideal → Dim = n (ganzer affiner Raum)
            - Ansonsten: n minus maximale Anzahl unabhängiger führender Monome

            Hier wird eine heuristische Näherung verwendet, die für
            typische algebraische Varietäten korrekte Werte liefert.
        @return Krull-Dimension (ganzzahlig, -1 für leere Varietät)
        @lastModified 2026-03-10
        """
        n = len(self.variables)
        gb = self._get_groebner()

        # Prüfe ob Ideal das gesamte Ideal ist (1 im Ideal)
        for g in gb:
            if g == sympy.Integer(1) or g == 1:
                return -1  # Leere Varietät

        # Nullideal → ganzer Raum
        if not self.generators or all(
            expand(g) == sympy.Integer(0) for g in self.generators
        ):
            return n

        # Heuristik: Zähle Variablen, die nicht als führende Variable in
        # irgendeinem Gröbner-Basis-Element erscheinen (freie Variablen).
        # Dies ist eine Näherung der echten Krull-Dimension.
        try:
            gb_obj = sympy_groebner(self.generators, *self.variables, order='lex')
            # Variablen, die als führende Terme (nur in einer Variablen) vorkommen
            constrained = set()
            for g in gb_obj:
                g_poly = Poly(g, *self.variables)
                lt = g_poly.LT()  # führendes Term als (Monom, Koeff)
                monom = lt[0]  # Tuple der Exponenten
                # Wenn das führende Monom nur eine Variable hat, ist diese gebunden
                nonzero = [i for i, e in enumerate(monom) if e > 0]
                if len(nonzero) == 1:
                    constrained.add(nonzero[0])
            return n - len(constrained)
        except Exception:
            # Fallback: Anzahl Generatoren als groben Schätzer
            return max(0, n - len(self.generators))

    def is_irreducible(self) -> bool:
        """
        @brief Prüft, ob die Varietät irreduzibel ist.
        @description
            Eine affine Varietät V ist irreduzibel, wenn sie nicht als Vereinigung
            zweier echter abgeschlossener Teilmengen geschrieben werden kann.

            Äquivalent: V(I) ist irreduzibel ⟺ √I ist ein Primideal.

            Hier wird geprüft, ob jedes Erzeuger-Polynom irreduzibel ist und
            ob die Menge prim ist (vereinfachte Heuristik für typische Fälle).
        @return True wenn irreduzibel, False sonst
        @lastModified 2026-03-10
        """
        # Trivialfall: Nullideal → ganzer affiner Raum ist irreduzibel
        if not self.generators or all(
            expand(g) == sympy.Integer(0) for g in self.generators
        ):
            return True

        # Einzelne Hyperfläche: irreduzibel gdw. Polynom irreduzibel (über ℚ)
        if len(self.generators) == 1:
            f = self.generators[0]
            if f == 0:
                return True
            # Faktorisiere über ℚ
            f_factored = factor(f)
            # Zähle verschiedene irreduzible Faktoren (ohne Vielfachheiten)
            from sympy import Mul, Pow
            if isinstance(f_factored, Mul):
                # Mehrere Faktoren → zerlegbar → reduzibel
                non_const = [
                    fac for fac in f_factored.args
                    if not fac.is_number and not (isinstance(fac, Pow) and fac.args[0].is_number)
                ]
                if len(non_const) > 1:
                    return False
            return True

        # Für mehrere Erzeuger: Prüfe ob Ideal ein Primideal ist
        # (vereinfachte Heuristik)
        try:
            # Versuche alle Erzeuger gemeinsam zu faktorisieren
            from sympy import groebner as gb_func
            gb = gb_func(self.generators, *self.variables, order='lex')
            if len(list(gb)) == 1:
                g = list(gb)[0]
                factored = factor(g)
                from sympy import Mul
                if isinstance(factored, Mul):
                    non_const = [
                        fac for fac in factored.args
                        if not fac.is_number
                    ]
                    return len(non_const) <= 1
        except Exception:
            pass
        return True

    def singular_points(self) -> list:
        """
        @brief Berechnet singuläre Punkte der Varietät (Jacobi-Kriterium).
        @description
            Ein Punkt p ∈ V(I) ist singulär, wenn der Rang der Jacobi-Matrix
            J(p) kleiner ist als die Kodimension der Varietät:
                rk(J(p)) < codim(V) = n - dim(V)

            Für eine Hyperfläche V(f) ist p singulär gdw.
                f(p) = 0  UND  ∂f/∂xᵢ(p) = 0 für alle i.

            Die Berechnung erfolgt symbolisch via SymPy:
            Das Gleichungssystem {f₁=0,...,fₘ=0, alle ∂fᵢ/∂xⱼ=0} wird gelöst.
        @return Liste von Punkten (als SymPy-Dictionaries {var: wert})
        @lastModified 2026-03-10
        """
        if not self.generators or all(
            expand(g) == sympy.Integer(0) for g in self.generators
        ):
            return []

        # Jacobi-Matrix aufbauen
        jac = _jacobian_matrix(self.generators, self.variables)

        # Gleichungssystem: alle f_i = 0 und alle ∂f_i/∂x_j = 0
        equations = list(self.generators)
        for row in jac:
            equations.extend(row)

        # Doppelte Nullen entfernen
        equations = [eq for eq in equations if expand(eq) != 0]

        try:
            solutions = sympy.solve(equations, self.variables, dict=True)
            return solutions
        except Exception:
            return []

    def intersect(self, other: 'AffineVariety') -> 'AffineVariety':
        """
        @brief Berechnet den Schnitt zweier Varietäten.
        @description
            Für zwei Varietäten V(I) und V(J) gilt:
                V(I) ∩ V(J) = V(I + J)
            wobei I + J das von den Erzeugern von I und J gemeinsam
            erzeugte Ideal ist.
        @param other Zweite affine Varietät (muss gleiche Variablen haben)
        @return Neue AffineVariety, die den Schnitt repräsentiert
        @raises ValueError wenn Variablen nicht kompatibel sind
        @lastModified 2026-03-10
        """
        if set(str(v) for v in self.variables) != set(str(v) for v in other.variables):
            raise ValueError(
                "Varietäten haben unterschiedliche Variablen und können nicht geschnitten werden."
            )
        # Schnitt entspricht der Summe der Ideale
        combined_gens = self.generators + other.generators
        return AffineVariety(combined_gens, self.variables)


# ===========================================================================
# KLASSE: ProjectiveVariety
# ===========================================================================

class ProjectiveVariety:
    """
    Projektive Varietät im projektiven Raum ℙⁿ.

    @description
        Eine projektive Varietät V(F₁,...,Fₘ) ⊂ ℙⁿ ist die Menge der
        gemeinsamen Nullstellen homogener Polynome F₁,...,Fₘ ∈ k[x₀,...,xₙ]:
            V(F₁,...,Fₘ) = { [p₀:...:pₙ] ∈ ℙⁿ : F₁(p) = ... = Fₘ(p) = 0 }

        Homogenität ist notwendig, damit die Nullstellenbedingung
        wohlbefiniert auf Äquivalenzklassen ist: Fᵢ(λp) = λ^d Fᵢ(p) ≠ 0
        falls Fᵢ(p) ≠ 0.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, homogeneous_poly: list, n_vars: int):
        """
        @brief Initialisiert die projektive Varietät.
        @param homogeneous_poly Liste homogener SymPy-Polynome
        @param n_vars           Anzahl der Variablen (Dimension = n_vars - 1)
        @raises ValueError wenn n_vars < 2 oder ein Polynom nicht homogen ist
        @lastModified 2026-03-10
        """
        if n_vars < 2:
            raise ValueError("n_vars muss mindestens 2 sein (projektiver Raum ≥ ℙ¹).")
        self.n_vars = n_vars
        # Erzeuger-Polynome normalisieren
        self.generators = [expand(f) for f in homogeneous_poly]
        # Standard-Variablen x₀,...,x_{n-1}
        self.variables = list(symbols(f'x0:{n_vars}'))

    def _poly_degree(self, f) -> int:
        """
        @brief Berechnet den Gesamtgrad eines Polynoms.
        @param f SymPy-Ausdruck
        @return Totaler Grad (max. Summe der Exponenten in einem Monom)
        @lastModified 2026-03-10
        """
        if f == 0:
            return 0
        try:
            p = Poly(f, *self.variables)
            return p.total_degree()
        except Exception:
            return 0

    def degree(self) -> int:
        """
        @brief Gibt den maximalen Grad der definierenden Polynome zurück.
        @description
            Der Grad einer projektiven Varietät (als Hyperfläche) ist der Grad
            des definierenden Polynoms. Bei mehreren Polynomen wird das Maximum
            zurückgegeben.

            Für eine ebene Kurve F(x,y,z) = 0 vom Grad d gilt nach Bézout:
            Zwei solche Kurven schneiden sich in d₁·d₂ Punkten (mit Vielfachheit).
        @return Maximaler Grad der Erzeuger-Polynome
        @lastModified 2026-03-10
        """
        if not self.generators:
            return 0
        return max(self._poly_degree(f) for f in self.generators)

    def is_projective_smooth(self) -> bool:
        """
        @brief Prüft, ob die projektive Varietät glatt ist.
        @description
            Eine projektive Hyperfläche V(F) ⊂ ℙⁿ ist glatt, wenn
            die Jacobi-Matrix an keinem Punkt der Varietät verschwindet:
                ∀p ∈ V(F): (∂F/∂x₀(p), ..., ∂F/∂xₙ(p)) ≠ 0

            Äquivalent: F und alle partiellen Ableitungen haben keine
            gemeinsame Nullstelle in ℙⁿ.

            Für die Hyperfläche-Bedingung wird geprüft, ob das System
            {F=0, ∂F/∂x₀=0, ..., ∂F/∂xₙ=0} keine nicht-triviale Lösung hat.
        @return True wenn glatt, False wenn singulär
        @lastModified 2026-03-10
        """
        if not self.generators:
            return True

        # Für jedes definierende Polynom Singularitätensystem lösen
        for f in self.generators:
            if f == 0:
                continue
            equations = [f] + [diff(f, v) for v in self.variables]
            equations = [expand(eq) for eq in equations]

            try:
                sols = sympy.solve(equations, self.variables, dict=True)
                # Nur nicht-triviale Lösungen zählen (nicht alle Variablen = 0)
                for sol in sols:
                    vals = [sol.get(v, 0) for v in self.variables]
                    if any(val != 0 for val in vals):
                        return False  # Singulärer Punkt gefunden
            except Exception:
                pass  # Bei Lösungsfehler: konservativ → als glatt annehmen

        return True


# ===========================================================================
# KLASSE: HilbertBasisTheorem
# ===========================================================================

class HilbertBasisTheorem:
    """
    Demonstration des Hilbert-Basissatzes.

    @description
        Hilbert-Basissatz: Jedes Ideal I ⊂ k[x₁,...,xₙ] ist endlich erzeugt.
        Das bedeutet, es gibt f₁,...,fₘ ∈ k[x₁,...,xₙ] mit I = ⟨f₁,...,fₘ⟩.

        Eine Gröbner-Basis ist eine spezielle endliche Erzeugermenge eines Ideals,
        die viele algorithmische Fragen lösbar macht:
        - Idealmitgliedschaft
        - Gleichungslösung
        - Eliminationstheorie

        Gröbner-Basen werden bezüglich einer Monomordnung berechnet
        (hier: lexikographische Ordnung 'lex').

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def groebner_basis(generators: list, variables: list) -> list:
        """
        @brief Berechnet die Gröbner-Basis des gegebenen Ideals.
        @description
            Verwendet den Buchberger-Algorithmus (via SymPy) zur Berechnung
            einer reduzierten Gröbner-Basis bezüglich der lexikographischen
            Monomordnung.

            Eigenschaften einer Gröbner-Basis G:
            - G erzeugt dasselbe Ideal wie die Eingabe
            - Für jeden Erzeuger g ∈ G ist kein anderes Leitmonom LM(gᵢ)
              ein Teiler von LM(g)
            - Reduktionseigenschaft: jedes f ∈ I reduziert auf 0 modulo G
        @param generators Liste von SymPy-Ausdrücken (Ideal-Erzeuger)
        @param variables  Liste von SymPy-Symbolen
        @return Liste der Gröbner-Basis-Elemente
        @lastModified 2026-03-10
        """
        if not generators or all(expand(g) == 0 for g in generators):
            return [sympy.Integer(0)]
        if not variables:
            return [sympy.Integer(1)]

        try:
            gb = sympy_groebner(generators, *variables, order='lex')
            return list(gb)
        except Exception as e:
            # Fallback: ursprüngliche Erzeuger zurückgeben
            return [expand(g) for g in generators if expand(g) != 0]

    @staticmethod
    def ideal_membership(poly, generators: list, variables: list) -> bool:
        """
        @brief Prüft, ob ein Polynom im Ideal der Erzeuger liegt.
        @description
            Verwendet die Gröbner-Basis-Reduktion:
            f ∈ ⟨f₁,...,fₘ⟩ ⟺ NF(f, G) = 0
            wobei G die Gröbner-Basis des Ideals und NF(f,G) die
            Normalform von f modulo G ist.
        @param poly       Zu prüfendes Polynom (SymPy-Ausdruck)
        @param generators Ideal-Erzeuger (Liste von SymPy-Ausdrücken)
        @param variables  Liste von SymPy-Symbolen
        @return True wenn poly im Ideal liegt, False sonst
        @lastModified 2026-03-10
        """
        f = expand(poly)
        if f == 0:
            return True
        if not generators or all(expand(g) == 0 for g in generators):
            return f == 0

        try:
            gb = sympy_groebner(generators, *variables, order='lex')
            # Normalform berechnen via gb.reduce(f) → (quotient_list, remainder)
            _, remainder = gb.reduce(f)
            return expand(remainder) == 0
        except Exception:
            # Fallback: direkter Vergleich mit Erzeugern
            try:
                gb_list = sympy_groebner(generators, *variables, order='lex')
                for g in gb_list:
                    if g != 0 and expand(f) == expand(g):
                        return True
                return False
            except Exception:
                return False


# ===========================================================================
# KLASSE: NullstellensatzDemo
# ===========================================================================

class NullstellensatzDemo:
    """
    Demonstration des Hilbert'schen Nullstellensatzes.

    @description
        Der Hilbert'sche Nullstellensatz verbindet algebraische Ideale mit
        geometrischen Varietäten über algebraisch abgeschlossenen Körpern.

        Schwacher Nullstellensatz:
            I ⊂ k[x₁,...,xₙ] hat keine gemeinsame Nullstelle (in k̄ⁿ)
            ⟺ 1 ∈ I ⟺ I = k[x₁,...,xₙ] (das Einheitsideal)

        Starker Nullstellensatz:
            Für ein Ideal I und ein Polynom f gilt:
            f|_{V(I)} = 0 ⟺ f^m ∈ I für ein m ≥ 1 ⟺ f ∈ √I

        Das Radikal √I = { f ∈ k[x] : ∃m≥1: f^m ∈ I }.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def weak_nullstellensatz(generators: list, variables: list) -> dict:
        """
        @brief Prüft den schwachen Nullstellensatz.
        @description
            Prüft, ob das Ideal I = ⟨f₁,...,fₘ⟩ das Einheitsideal ist
            (d.h. 1 ∈ I ⟺ V(I) = ∅).

            Methode: Berechne Gröbner-Basis G. Dann gilt 1 ∈ I ⟺ G = {1}.
        @param generators Liste von Ideal-Erzeugern
        @param variables  Liste von SymPy-Symbolen
        @return Dictionary mit:
                  'has_no_solution': bool (True wenn V(I) = ∅)
                  'groebner_basis':  Liste der Gröbner-Basis
                  'explanation':     Text-Erklärung
        @lastModified 2026-03-10
        """
        if not generators or not variables:
            return {
                'has_no_solution': False,
                'groebner_basis': [sympy.Integer(0)],
                'explanation': 'Leeres Ideal oder keine Variablen: V(I) = ganzer Raum.'
            }

        try:
            gb = sympy_groebner(generators, *variables, order='lex')
            gb_list = list(gb)
            # 1 ∈ I gdw. Gröbner-Basis = {1}
            is_unit = (len(gb_list) == 1 and gb_list[0] == sympy.Integer(1))
            return {
                'has_no_solution': is_unit,
                'groebner_basis': gb_list,
                'explanation': (
                    'Das Ideal ist das Einheitsideal (1 ∈ I). '
                    'Nach dem schwachen Nullstellensatz hat V(I) keine Lösung.'
                    if is_unit else
                    'Das Ideal ist kein Einheitsideal. V(I) hat mindestens eine Lösung.'
                )
            }
        except Exception as e:
            return {
                'has_no_solution': False,
                'groebner_basis': generators,
                'explanation': f'Fehler bei der Berechnung: {e}'
            }

    @staticmethod
    def strong_nullstellensatz(f, generators: list, variables: list) -> dict:
        """
        @brief Prüft den starken Nullstellensatz (Radikalmitgliedschaft).
        @description
            Prüft, ob f ∈ √I, d.h. ob f auf der Varietät V(I) verschwindet.

            Nach dem starken Nullstellensatz gilt:
            f ∈ √I ⟺ ∃m≥1: f^m ∈ I

            Methode: Rabinowitsch-Trick:
            Führe neue Variable t ein. Dann gilt:
            f ∈ √I ⟺ 1 ∈ ⟨f₁,...,fₘ, 1 - t·f⟩ in k[x₁,...,xₙ,t]
        @param f          Polynom, dessen Radikalmitgliedschaft geprüft wird
        @param generators Liste von Ideal-Erzeugern
        @param variables  Liste von SymPy-Symbolen
        @return Dictionary mit:
                  'in_radical': bool
                  'explanation': Text-Erklärung
        @lastModified 2026-03-10
        """
        f_exp = expand(f)

        # Trivialfall: f = 0 ist immer im Radikal
        if f_exp == 0:
            return {'in_radical': True, 'explanation': '0 liegt in jedem Radikal.'}

        if not generators or all(expand(g) == 0 for g in generators):
            return {
                'in_radical': False,
                'explanation': 'Nullideal: Radikal = {0}. f≠0 liegt nicht darin.'
            }

        # Rabinowitsch-Trick: neue Variable t einführen
        t = symbols('_rabinowitsch_t')
        extended_vars = list(variables) + [t]
        extended_gens = list(generators) + [1 - t * f_exp]

        try:
            gb = sympy_groebner(extended_gens, *extended_vars, order='lex')
            gb_list = list(gb)
            # f ∈ √I gdw. 1 ∈ erweitertes Ideal
            in_radical = (len(gb_list) == 1 and gb_list[0] == sympy.Integer(1))
            return {
                'in_radical': in_radical,
                'explanation': (
                    f'f liegt im Radikal √I (Rabinowitsch-Trick: 1 ∈ erweitertes Ideal).'
                    if in_radical else
                    f'f liegt NICHT im Radikal √I.'
                )
            }
        except Exception as e:
            return {
                'in_radical': False,
                'explanation': f'Fehler bei der Berechnung: {e}'
            }


# ===========================================================================
# KLASSE: EllipticCurveGeometry
# ===========================================================================

class EllipticCurveGeometry:
    """
    Elliptische Kurve in Weierstrass-Normalform y² = x³ + ax + b.

    @description
        Eine elliptische Kurve über einem Körper k ist eine nicht-singuläre
        projektive Kurve vom Geschlecht 1 mit einem ausgezeichneten Punkt.
        In Weierstrass-Normalform:
            E: y² = x³ + ax + b  (a, b ∈ k)

        Diskriminante: Δ = -16(4a³ + 27b²)
        - E ist glatt ⟺ Δ ≠ 0
        - E hat eine Nodalsingularität wenn 4a³ + 27b² = 0

        j-Invariante: j(E) = -1728 · (4a)³ / Δ = 1728 · 4a³·(-16) / (Δ)
        Klassifiziert elliptische Kurven bis auf Isomorphie.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, a, b):
        """
        @brief Initialisiert die elliptische Kurve y² = x³ + ax + b.
        @param a Koeffizient von x (SymPy-Ausdruck oder Zahl)
        @param b Konstantterm (SymPy-Ausdruck oder Zahl)
        @lastModified 2026-03-10
        """
        self.a = sympy.sympify(a)
        self.b = sympy.sympify(b)

    def discriminant(self):
        """
        @brief Berechnet die Diskriminante Δ = -16(4a³ + 27b²).
        @description
            Die Diskriminante entscheidet über die Singularität der Kurve:
            - Δ ≠ 0: glatte elliptische Kurve
            - Δ = 0: singuläre Kurve (Knoten oder Spitze)
        @return Diskriminante als SymPy-Ausdruck
        @lastModified 2026-03-10
        """
        return simplify(-16 * (4 * self.a**3 + 27 * self.b**2))

    def j_invariant(self):
        """
        @brief Berechnet die j-Invariante j(E) = -1728·(4a)³/Δ.
        @description
            Die j-Invariante klassifiziert elliptische Kurven über ℂ
            bis auf Isomorphie: E₁ ≅ E₂ ⟺ j(E₁) = j(E₂).

            Spezialfälle:
            - j = 0:    a = 0 (Kurve hat CM durch ℤ[e^{2πi/3}])
            - j = 1728: b = 0 (Kurve hat CM durch ℤ[i])
        @return j-Invariante als SymPy-Ausdruck (oder oo wenn Δ=0)
        @lastModified 2026-03-10
        """
        delta = self.discriminant()
        if delta == 0:
            return sympy.oo  # Singuläre Kurve: j-Invariante nicht definiert
        numerator = -1728 * (4 * self.a)**3
        return simplify(numerator / delta)

    def is_smooth(self) -> bool:
        """
        @brief Prüft, ob die elliptische Kurve glatt ist (Δ ≠ 0).
        @description
            Eine Weierstrass-Kurve ist genau dann eine elliptische Kurve
            (nicht-singulär), wenn ihre Diskriminante von Null verschieden ist.
        @return True wenn glatt (Δ ≠ 0), False wenn singulär
        @lastModified 2026-03-10
        """
        delta = self.discriminant()
        return delta != 0 and simplify(delta) != 0

    def rational_points_mod_p(self, p: int) -> list:
        """
        @brief Berechnet rationale Punkte auf E mod p (über 𝔽_p).
        @description
            Bestimmt alle affinen Punkte (x, y) ∈ 𝔽_p² mit y² ≡ x³+ax+b (mod p)
            sowie den Punkt im Unendlichen O = ∞.

            Algorithmus:
            Für jedes x ∈ {0,...,p-1}:
                Berechne rhs = x³ + ax + b (mod p)
                Für jedes y ∈ {0,...,p-1}: prüfe y² ≡ rhs (mod p)
        @param p Primzahl (Charakteristik des endlichen Körpers)
        @return Liste von Tupeln (x, y) ∈ 𝔽_p² plus 'O' für den Punkt ∞
        @raises ValueError wenn p keine positive Ganzzahl ist
        @lastModified 2026-03-10
        """
        if not isinstance(p, int) or p < 2:
            raise ValueError(f"p muss eine Primzahl ≥ 2 sein, erhalten: {p}")

        # Koeffizienten mod p auswerten
        try:
            a_val = int(self.a) % p
            b_val = int(self.b) % p
        except Exception:
            # Symbolische Koeffizienten: als ganze Zahlen auswerten
            a_val = int(float(self.a)) % p
            b_val = int(float(self.b)) % p

        points = ['O']  # Punkt im Unendlichen

        for xv in range(p):
            # Rechte Seite: x³ + ax + b mod p
            rhs = (xv**3 + a_val * xv + b_val) % p
            for yv in range(p):
                if (yv * yv) % p == rhs:
                    points.append((xv, yv))

        return points


# ===========================================================================
# KLASSE: IntersectionTheory
# ===========================================================================

class IntersectionTheory:
    """
    Schnitttheorie für algebraische Kurven.

    @description
        Untersucht Schnittpunkte algebraischer Kurven im projektiven Raum
        mit Methoden der algebraischen Geometrie.

        Bézout-Theorem: Zwei projektive Kurven vom Grad d₁ und d₂ über
        einem algebraisch abgeschlossenen Körper schneiden sich in genau
        d₁·d₂ Punkten (mit Multiplizitäten gezählt).

        Lokale Schnittmultiplizität: Misst, wie oft sich zwei Kurven
        in einem bestimmten Punkt berühren.

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    @staticmethod
    def bezout_theorem(d1: int, d2: int) -> int:
        """
        @brief Gibt die Bézout-Zahl zweier Kurven zurück.
        @description
            Bézout-Theorem (1779): Zwei algebraische Kurven vom Grad d₁ und d₂
            im projektiven Raum ℙ² über einem algebraisch abgeschlossenen Körper
            schneiden sich in genau d₁·d₂ Punkten, wenn man:
            1. Schnittpunkte mit Multiplizitäten zählt
            2. Punkte im Unendlichen berücksichtigt
            3. Komplexe Schnittpunkte einschließt

            Beispiele:
            - Zwei Geraden (d₁=d₂=1): 1 Schnittpunkt
            - Gerade und Kegelschnitt (d₁=1, d₂=2): 2 Schnittpunkte
            - Zwei Kegelschnitte (d₁=d₂=2): 4 Schnittpunkte
        @param d1 Grad der ersten Kurve (≥ 0)
        @param d2 Grad der zweiten Kurve (≥ 0)
        @return Anzahl der Schnittpunkte (mit Vielfachheiten) = d₁·d₂
        @lastModified 2026-03-10
        """
        return d1 * d2

    @staticmethod
    def intersection_multiplicity(f, g, point: dict, variables: list) -> int:
        """
        @brief Berechnet die lokale Schnittmultiplizität zweier Kurven in einem Punkt.
        @description
            Die Schnittmultiplizität I_p(f, g) zweier Kurven f=0 und g=0
            in einem Punkt p misst die Ordnung der Berührung.

            Definition via lokale Algebra:
            I_p(f,g) = dim_k( 𝒪_{ℙ²,p} / (f,g) )

            Praktische Berechnung: Zähle die Nullstellen des resultierenden
            Einvariablen-Polynoms, das durch Eliminierung entsteht.

            Für einfache Schnittpunkte gilt I_p(f,g) = 1.
            Tangentiale Berührungen haben I_p(f,g) ≥ 2.

            Methode: Substituiere x = point[x] + t in f und g, entwickle
            nach t und bestimme die Ordnung (Valuation) des GCD.
        @param f         Erstes Polynom (SymPy-Ausdruck)
        @param g         Zweites Polynom (SymPy-Ausdruck)
        @param point     Punkt als Dictionary {Symbol: Wert}
        @param variables Liste der SymPy-Symbole
        @return Schnittmultiplizität (natürliche Zahl ≥ 0)
        @lastModified 2026-03-10
        """
        # Verschiebe Koordinaten so dass der Punkt im Ursprung liegt
        subs_dict = {v: v + point.get(v, 0) for v in variables}
        f_shifted = expand(f.subs({v: v + point.get(v, 0) for v in variables}))
        g_shifted = expand(g.subs({v: v + point.get(v, 0) for v in variables}))

        # Multiplizität über das lokale Ideal berechnen
        # Methode: Dimension des Vektorraums k[x,y]/⟨f,g⟩_{lokalisiert an Ursprung}
        # Näherungsweise: Ordnung des GCD der Leitpolynome
        try:
            if len(variables) == 2:
                x_var, y_var = variables[0], variables[1]
                # Resultante bezüglich einer Variablen
                res = sympy.resultant(f_shifted, g_shifted, x_var)
                if res == 0:
                    return sympy.oo  # Unendliche Multiplizität (gemeinsame Komponente)
                # Ordnung der Resultante in y_var am Ursprung
                res_poly = Poly(res, y_var)
                coeffs = res_poly.all_coeffs()
                # Zähle Nullen am Anfang (Ordnung bei y=0)
                order = 0
                for c in reversed(coeffs):
                    if expand(c) == 0:
                        order += 1
                    else:
                        break
                return order
            else:
                # Allgemeine Heuristik: beide Polynome am Punkt auswerten
                f_val = expand(f_shifted.subs({v: 0 for v in variables}))
                g_val = expand(g_shifted.subs({v: 0 for v in variables}))
                if f_val == 0 and g_val == 0:
                    return 1  # Mindestens einfacher Schnittpunkt
                return 0
        except Exception:
            return 1  # Fallback: einfacher Schnittpunkt


# ===========================================================================
# KLASSE: MorphismOfVarieties
# ===========================================================================

class MorphismOfVarieties:
    """
    Morphismus zwischen algebraischen Varietäten.

    @description
        Ein Morphismus φ: V → W zwischen affinen Varietäten ist eine
        Abbildung, die durch Polynome gegeben ist:
            φ(p) = (φ₁(p), ..., φₘ(p))
        wobei φᵢ ∈ k[x₁,...,xₙ] Polynome in den Koordinaten von V sind.

        Morphismen entsprechen (kontravarianten) Ringhomomorphismen
        zwischen den Koordinatenringen:
            φ*: k[W] → k[V],  f ↦ f ∘ φ  (Pullback)

        Dominanz: φ ist dominant, wenn der Pullback φ* injektiv ist,
        äquivalent: φ(V) ist dicht in W (Zariski-Topologie).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, source: AffineVariety, target: AffineVariety,
                 map_polynomials: list):
        """
        @brief Initialisiert den Varietäten-Morphismus.
        @param source          Quell-Varietät
        @param target          Ziel-Varietät
        @param map_polynomials Liste von Polynomen [φ₁,...,φₘ] in den
                               Koordinaten der Quell-Varietät
        @lastModified 2026-03-10
        """
        self.source = source
        self.target = target
        # Abbildungspolynome normalisieren
        self.map_polynomials = [expand(p) for p in map_polynomials]

    def is_dominant(self) -> bool:
        """
        @brief Prüft, ob der Morphismus dominant ist.
        @description
            φ: V → W ist dominant, wenn das Bild φ(V) Zariski-dicht in W liegt.

            Äquivalente Bedingung: Der Pullback φ*: k[W] → k[V] ist injektiv.

            Praktische Heuristik: φ ist dominant, wenn die Jacobi-Matrix dφ
            an einem generischen Punkt vollen Rang hat (Dimensionsgleichheit).

            Vereinfacht: Prüfe ob die Anzahl der Abbildungspolynome mit der
            Zieldimension übereinstimmt und die Ableitung nicht identisch Null ist.
        @return True wenn dominant (heuristische Schätzung), False sonst
        @lastModified 2026-03-10
        """
        if not self.map_polynomials:
            return False

        # Prüfe ob Abbildungspolynome alle Variablen der Quelle effektiv nutzen
        src_vars = set(str(v) for v in self.source.variables)
        used_vars = set()
        for phi in self.map_polynomials:
            used_vars |= set(str(v) for v in phi.free_symbols)

        # Dominant wenn alle Quell-Variablen verwendet werden
        # und die Jacobi-Matrix nicht identisch Null ist
        jac_nonzero = False
        for phi in self.map_polynomials:
            for v in self.source.variables:
                if expand(diff(phi, v)) != 0:
                    jac_nonzero = True
                    break
            if jac_nonzero:
                break

        return jac_nonzero and (used_vars & src_vars) != set()

    def pullback(self, poly) -> object:
        """
        @brief Berechnet den Pullback einer Funktion auf der Ziel-Varietät.
        @description
            Der Pullback φ*(f) einer Funktion f auf W ist definiert durch:
                φ*(f) = f ∘ φ

            Konkret: Ersetze in f die Koordinaten y₁,...,yₘ von W durch
            die Abbildungspolynome φ₁,...,φₘ.
        @param poly Polynom in den Koordinaten der Ziel-Varietät
        @return φ*(poly) als Polynom in den Koordinaten der Quell-Varietät
        @lastModified 2026-03-10
        """
        f = expand(poly)
        # Ersetze Ziel-Variablen durch Abbildungspolynome
        target_vars = self.target.variables
        subs_dict = {}
        for i, var in enumerate(target_vars):
            if i < len(self.map_polynomials):
                subs_dict[var] = self.map_polynomials[i]
        result = f.subs(subs_dict)
        return expand(result)


# ===========================================================================
# STANDALONE-FUNKTIONEN
# ===========================================================================

def compute_groebner_basis(gens: list, vars: list) -> list:
    """
    @brief Berechnet die Gröbner-Basis eines Ideals (lex-Ordnung).
    @description
        Wrapper um SymPy's groebner-Funktion. Berechnet eine reduzierte
        Gröbner-Basis bezüglich der lexikographischen Monomordnung.
    @param gens Liste der Ideal-Erzeuger (SymPy-Ausdrücke)
    @param vars Liste der SymPy-Symbole
    @return Liste der Gröbner-Basis-Elemente
    @lastModified 2026-03-10
    """
    return HilbertBasisTheorem.groebner_basis(gens, vars)


def ideal_radical(gens: list, vars: list) -> list:
    """
    @brief Berechnet das Radikal √I eines Ideals.
    @description
        Das Radikal √I = { f ∈ k[x] : ∃m≥1: f^m ∈ I } ist das kleinste
        Radikal-Ideal, das I enthält.

        Für ein Primideal gilt √I = I.
        Nach Hilberts Nullstellensatz ist I(V(I)) = √I.

        Näherungsweise Berechnung: Für einige Spezialfälle exakt, sonst
        Rückgabe der Gröbner-Basis (konservative Näherung).
    @param gens Liste der Ideal-Erzeuger
    @param vars Liste der SymPy-Symbole
    @return Liste der Erzeuger des Radikals
    @lastModified 2026-03-10
    """
    if not gens or all(expand(g) == 0 for g in gens):
        return [sympy.Integer(0)]

    try:
        # Für ein durch ein Polynom erzeugtes Ideal: Radikal = quadratfreier Teil
        if len(gens) == 1:
            f = expand(gens[0])
            if f == 0:
                return [sympy.Integer(0)]
            # Quadratfreier Teil = Produkt aller irreduziblen Faktoren
            f_factored = factor(f)
            from sympy import Mul, Pow, Number
            radical_factors = []
            if isinstance(f_factored, Mul):
                for fac in f_factored.args:
                    if isinstance(fac, Pow):
                        base, exp = fac.args
                        if not base.is_number:
                            radical_factors.append(base)
                    elif not fac.is_number:
                        radical_factors.append(fac)
            elif isinstance(f_factored, Pow):
                base, exp = f_factored.args
                if not base.is_number:
                    radical_factors.append(base)
            else:
                radical_factors.append(f_factored)

            if radical_factors:
                result = sympy.Integer(1)
                for fac in radical_factors:
                    result = expand(result * fac)
                return [result]
            return [f_factored]
        else:
            # Allgemeiner Fall: Gröbner-Basis als Näherung
            gb = sympy_groebner(gens, *vars, order='lex')
            return list(gb)
    except Exception:
        return [expand(g) for g in gens]


def zariski_closure_demo(points: list, vars: list) -> dict:
    """
    @brief Demonstriert die Zariski-Abschluss-Konstruktion.
    @description
        Der Zariski-Abschluss einer Menge S ⊂ 𝔸ⁿ ist die kleinste
        affine Varietät, die S enthält:
            Z̄(S) = V(I(S))
        wobei I(S) = { f ∈ k[x] : f(p) = 0 ∀p ∈ S }.

        Für endlich viele Punkte p₁,...,pₖ ∈ 𝔸ⁿ ist der Zariski-Abschluss
        eine Varietät, die durch Lagrange-Interpolation beschrieben werden kann.

        Für einen einzigen Punkt p = (a₁,...,aₙ) ist:
        Z̄({p}) = V(x₁-a₁, ..., xₙ-aₙ) (maximales Ideal von p).
    @param points Liste von Punkten als Dictionaries {Symbol: Wert}
    @param vars   Liste der SymPy-Symbole
    @return Dictionary mit 'ideal_generators' und 'description'
    @lastModified 2026-03-10
    """
    if not points:
        return {
            'ideal_generators': [sympy.Integer(1)],
            'description': 'Zariski-Abschluss der leeren Menge ist ∅ = V(1).'
        }

    if len(points) == 1:
        # Einziger Punkt: erzeuge maximales Ideal
        pt = points[0]
        generators = [v - pt.get(v, 0) for v in vars]
        return {
            'ideal_generators': generators,
            'description': f'Zariski-Abschluss eines Punktes: V({", ".join(str(g) for g in generators)})'
        }

    # Mehrere Punkte: Produkt der Linearformen
    # I(p₁,...,pₖ) erzeugt durch Produkt der Punkt-Ideale
    # Für 1D: Polynom mit Nullstellen an den Punkten
    if len(vars) == 1:
        v = vars[0]
        poly_expr = sympy.Integer(1)
        for pt in points:
            poly_expr = expand(poly_expr * (v - pt.get(v, 0)))
        return {
            'ideal_generators': [poly_expr],
            'description': f'Zariski-Abschluss von {len(points)} Punkten: V({poly_expr})'
        }

    # Mehrdimensionaler Fall: Schnittpunkt der Punkt-Ideale
    all_gens = []
    for pt in points:
        for v in vars:
            all_gens.append(v - pt.get(v, 0))

    return {
        'ideal_generators': all_gens,
        'description': f'Zariski-Abschluss von {len(points)} Punkten in {len(vars)} Variablen.'
    }


def hilbert_polynomial(gens: list, vars: list):
    """
    @brief Berechnet das Hilbert-Polynom des Quotienten-Rings k[x]/I.
    @description
        Das Hilbert-Polynom P(t) des graduierten Rings R = k[x₁,...,xₙ]/I
        gibt für große t die Dimension des t-ten graduierten Teils an:
            dim_k(R_t) = P(t)  für t >> 0

        Für ein Ideal I ⊂ k[x,y] einer ebenen Kurve vom Grad d:
            P(t) = d·t - d(d-3)/2 - 1  (Riemann-Roch)

        Vereinfachte Berechnung via Hilbert-Funktion:
        P(t) = Σ_{k=0}^{t} dim(R_k) für kleine t.
    @param gens Liste der Ideal-Erzeuger
    @param vars Liste der SymPy-Symbole
    @return SymPy-Ausdruck (Polynom in t) oder Näherungswert
    @lastModified 2026-03-10
    """
    t = symbols('t')
    n = len(vars)

    if not gens or all(expand(g) == 0 for g in gens):
        # Nullideal: Hilbert-Polynom = C(t+n, n) (Binomialkoeffizient)
        # Für polynomialen Ring k[x₁,...,xₙ]: P(t) = C(t+n,n)
        from sympy import binomial
        return binomial(t + n, n)

    # Für Hauptideal ⟨f⟩: Grad d → P(t) = C(t+n,n) - C(t+n-d,n)
    if len(gens) == 1:
        f = expand(gens[0])
        if f == 0:
            from sympy import binomial
            return binomial(t + n, n)
        try:
            d = Poly(f, *vars).total_degree()
            from sympy import binomial
            p = binomial(t + n, n) - binomial(t + n - d, n)
            return sympy.expand(p)
        except Exception:
            pass

    # Allgemeiner Fall: Näherung über Gröbner-Basis
    try:
        gb = sympy_groebner(gens, *vars, order='grlex')
        # Höchster Grad in der Gröbner-Basis
        max_deg = max(Poly(g, *vars).total_degree() for g in gb if g != 0)
        from sympy import binomial
        p = binomial(t + n, n) - binomial(t + n - max_deg, n)
        return sympy.expand(p)
    except Exception:
        # Fallback: lineares Polynom
        return t + 1
