"""
@file differential_topology.py
@brief Differentialtopologie – Mannigfaltigkeiten, Differentialformen, Morse-Theorie und mehr.
@description
    Implementiert grundlegende Konzepte der Differentialtopologie:

    1. SmoothManifold – Glatte Mannigfaltigkeit (abstrakt)
       Repräsentiert eine n-dimensionale glatte Mannigfaltigkeit M.
       Tangential- und Kotangentialräume haben jeweils Dimension n.

    2. DifferentialForm – k-Differentialform auf einer offenen Menge des ℝⁿ
       Äußere Ableitung dω via SymPy.
       Geschlossene Formen (dω = 0) und exakte Formen (ω = dη).

    3. DeRhamCohomology – de Rham-Kohomologie H^k_dR(M)
       Betti-Zahlen für Sphären und Tori.
       Euler-Charakteristik χ = Σ (-1)^k b_k.

    4. MorseTheory – Morse-Theorie
       Kritische Punkte (∇f = 0), Morse-Index (Anzahl negativer Eigenwerte der Hesse-Matrix),
       Morse-Ungleichungen und Henkelzerlegung.

    5. TransversalityTheory – Transversalitätsbedingungen
       Sard'sches Theorem (kritische Werte bilden eine Menge vom Maß 0).

    6. VectorBundle – Vektorbündel E → B
       Totalraum, Rang, Trivialität.

    7. CharacteristicClasses – Charakteristische Klassen
       Chern-, Pontryagin-, Euler- und Stiefel-Whitney-Klassen (algebraisch/demo).

    8. SphereTopology – Topologie der n-Sphäre S^n
       Homotopiegruppen πₖ(Sⁿ) für k ≤ n+3, Homologiegruppen, Orientierbarkeit.

    9. Standalone-Funktionen:
       degree_of_map, hairy_ball_theorem_demo, poincare_hopf_theorem_demo,
       whitney_embedding_dimension, compute_jacobian, implicit_function_theorem_check.

    Mathematischer Hintergrund:
    Die Differentialtopologie verbindet Analysis und Topologie. Zentrale Sätze sind
    Sard (generische Transversalität), Morse (Topologie aus kritischen Punkten) und
    de Rham (Kohomologie via Differentialformen).

@author Michael Fuhrmann
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import numpy as np
import sympy as sp
from sympy import symbols, diff, Matrix, hessian, solve, simplify, Rational
from sympy import zeros as sp_zeros
from typing import List, Dict, Tuple, Optional, Any, Callable
import itertools


# ===========================================================================
# 1. GLATTE MANNIGFALTIGKEIT
# ===========================================================================

class SmoothManifold:
    """
    Repräsentiert eine abstrakte glatte (C∞-)Mannigfaltigkeit der Dimension n.

    Eine n-dimensionale glatte Mannigfaltigkeit M ist ein topologischer Raum,
    der lokal wie der ℝⁿ aussieht (d.h. jeden Punkt umgibt eine offene Menge,
    die homöomorph zu einer offenen Teilmenge des ℝⁿ ist). Alle Kartenwechsel
    sind glatte Abbildungen (C∞).

    Tangentialraum T_pM und Kotangentialraum T*_pM haben beide Dimension n.

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(self, dimension: int, name: str = 'M') -> None:
        """
        Initialisiert eine glatte Mannigfaltigkeit.

        @param dimension: Dimension n der Mannigfaltigkeit (n ≥ 0)
        @param name: Bezeichner der Mannigfaltigkeit (Standard: 'M')
        @raises ValueError: Wenn dimension < 0
        @lastModified 2026-03-10
        """
        if dimension < 0:
            raise ValueError(f"Dimension muss ≥ 0 sein, erhalten: {dimension}")
        # Dimension der Mannigfaltigkeit speichern
        self.dimension = dimension
        # Name zur lesbaren Ausgabe
        self.name = name

    def tangent_space_dim(self) -> int:
        """
        Gibt die Dimension des Tangentialraums T_pM zurück.

        Der Tangentialraum T_pM an einen Punkt p ∈ M hat dieselbe Dimension
        wie die Mannigfaltigkeit selbst: dim(T_pM) = dim(M) = n.

        @return: Dimension des Tangentialraums
        @lastModified 2026-03-10
        """
        # Tangentialraum hat gleiche Dimension wie die Mannigfaltigkeit
        return self.dimension

    def cotangent_space_dim(self) -> int:
        """
        Gibt die Dimension des Kotangentialraums T*_pM zurück.

        Der Kotangentialraum T*_pM = (T_pM)* ist der Dualraum des Tangentialraums.
        Da Dualräume endlichdimensionaler Vektorräume die gleiche Dimension haben,
        gilt: dim(T*_pM) = dim(T_pM) = n.

        @return: Dimension des Kotangentialraums
        @lastModified 2026-03-10
        """
        # Kotangentialraum hat als Dualraum gleiche Dimension
        return self.dimension

    def __repr__(self) -> str:
        """
        Lesbare Darstellung der Mannigfaltigkeit.

        @return: Beschreibungsstring
        @lastModified 2026-03-10
        """
        return f"SmoothManifold(name='{self.name}', dim={self.dimension})"


# ===========================================================================
# 2. DIFFERENTIALFORM
# ===========================================================================

class DifferentialForm:
    """
    Repräsentiert eine k-Differentialform ω auf einer offenen Teilmenge des ℝⁿ.

    Eine k-Form ω ist ein schiefsymmetrischer (0,k)-Tensor, also eine Abbildung,
    die k Tangentialvektoren antisymmetrisch in eine reelle Zahl überführt.

    In lokalen Koordinaten (x₁, …, xₙ) wird eine k-Form geschrieben als:
        ω = Σ_{i₁ < … < iₖ} f_{i₁…iₖ}(x) dx_{i₁} ∧ … ∧ dx_{iₖ}

    Die coefficients werden als Dict { (i₁,…,iₖ) : f_expr } übergeben,
    wobei die Indizes als aufsteigendes Tupel ganzer Zahlen (0-basiert) dienen.

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        degree: int,
        coefficients: Dict[Tuple, Any],
        variables: List
    ) -> None:
        """
        Initialisiert eine k-Differentialform.

        @param degree: Grad k der Form (k = 0: Funktion, k = 1: 1-Form, …)
        @param coefficients: Dict { multi_index_tuple : sympy_expression }
                             Multi-Indizes müssen aufsteigende Tupel sein.
                             Für k=0 wird der Schlüssel () (leeres Tupel) verwendet.
        @param variables: Liste der SymPy-Symbole [x₁, …, xₙ]
        @raises ValueError: Wenn k < 0 oder Koeffiziententupel falsche Länge haben
        @lastModified 2026-03-10
        """
        if degree < 0:
            raise ValueError(f"Grad der Differentialform muss ≥ 0 sein, erhalten: {degree}")
        # Grad der Form speichern
        self.degree = degree
        # Koeffizientendict: multi_index → SymPy-Ausdruck
        self.coefficients = coefficients
        # Koordinatenvariablen
        self.variables = list(variables)
        # Plausibilitätscheck: alle Schlüssel müssen Länge k haben
        for idx in coefficients:
            if len(idx) != degree:
                raise ValueError(
                    f"Multi-Index {idx} hat Länge {len(idx)}, erwartet {degree}"
                )

    def wedge(self, other: 'DifferentialForm') -> 'DifferentialForm':
        """
        Berechnet das äußere Produkt (Wedge-Produkt) ω ∧ η.

        Das Wedge-Produkt einer k-Form ω mit einer l-Form η ergibt eine (k+l)-Form.
        Es ist bilinear und schiefsymmetrisch: ω ∧ η = (-1)^{kl} η ∧ ω.

        Für Koeffizienten gilt:
            (ω ∧ η)_{i₁…iₖ j₁…jₗ} = ω_{i₁…iₖ} · η_{j₁…jₗ}
        unter Berücksichtigung der Schiefsymmetrie (Vorzeichen durch Sortierung).

        @param other: l-Form η
        @return: (k+l)-Form ω ∧ η
        @raises ValueError: Wenn die Variablenlisten nicht übereinstimmen
        @lastModified 2026-03-10
        """
        # Prüfen: gleiche Variablen
        if self.variables != other.variables:
            raise ValueError("Wedge-Produkt: Variablenlisten müssen übereinstimmen")

        new_degree = self.degree + other.degree
        new_coeffs: Dict[Tuple, Any] = {}

        # Alle Paare von Multi-Indizes kombinieren
        for idx1, f1 in self.coefficients.items():
            for idx2, f2 in other.coefficients.items():
                # Zusammenführen und sortieren
                combined = list(idx1) + list(idx2)
                # Vorzeichen durch Bubble-Sort bestimmen (Anzahl der Transpositionen)
                sign, sorted_idx = _sign_and_sort(combined)
                new_idx = tuple(sorted_idx)
                # Produkt mit Vorzeichen addieren (schiefsymmetrisch)
                contrib = sign * f1 * f2
                if new_idx in new_coeffs:
                    new_coeffs[new_idx] = simplify(new_coeffs[new_idx] + contrib)
                else:
                    new_coeffs[new_idx] = simplify(contrib)

        # Koeffizienten, die 0 sind, entfernen
        new_coeffs = {k: v for k, v in new_coeffs.items() if simplify(v) != 0}

        return DifferentialForm(new_degree, new_coeffs, self.variables)

    def exterior_derivative(self) -> 'DifferentialForm':
        """
        Berechnet die äußere Ableitung dω der k-Form ω.

        Die äußere Ableitung ist der fundamentale Differentialoperator auf Formen.
        Für eine k-Form ω = Σ f_I dx_I gilt:
            dω = Σ_I Σ_j (∂f_I/∂x_j) dx_j ∧ dx_I

        Eigenschaften:
            - d ∘ d = 0 (d² = 0)
            - d ist eine (k+1)-Form
            - Leibniz-Regel: d(ω ∧ η) = dω ∧ η + (-1)^k ω ∧ dη

        @return: (k+1)-Form dω
        @lastModified 2026-03-10
        """
        new_degree = self.degree + 1
        new_coeffs: Dict[Tuple, Any] = {}

        # Über alle Koeffizienten iterieren
        for idx, f_expr in self.coefficients.items():
            # Über alle Variablen partiell ableiten
            for j, var in enumerate(self.variables):
                # Partielle Ableitung via SymPy
                df_dxj = diff(f_expr, var)
                df_dxj_simplified = simplify(df_dxj)
                if df_dxj_simplified == 0:
                    # Null-Beitrag überspringen
                    continue
                # Neuer Multi-Index: j voranstellen, dann ursprünglicher Index
                new_idx_list = [j] + list(idx)
                sign, sorted_idx = _sign_and_sort(new_idx_list)
                new_idx = tuple(sorted_idx)
                contrib = sign * df_dxj_simplified
                if new_idx in new_coeffs:
                    new_coeffs[new_idx] = simplify(new_coeffs[new_idx] + contrib)
                else:
                    new_coeffs[new_idx] = contrib

        # Null-Koeffizienten entfernen
        new_coeffs = {k: v for k, v in new_coeffs.items() if simplify(v) != 0}

        return DifferentialForm(new_degree, new_coeffs, self.variables)

    def is_closed(self) -> bool:
        """
        Prüft, ob die Form geschlossen ist: dω = 0.

        Eine k-Form ω heißt geschlossen, wenn ihre äußere Ableitung verschwindet:
            dω = 0
        Exakte Formen sind stets geschlossen (d² = 0), aber nicht umgekehrt.
        Der Unterschied wird durch die de Rham-Kohomologie gemessen.

        @return: True wenn dω = 0, False sonst
        @lastModified 2026-03-10
        """
        # Äußere Ableitung berechnen
        d_omega = self.exterior_derivative()
        # Geschlossen ↔ alle Koeffizienten sind Null
        return len(d_omega.coefficients) == 0

    def is_exact(self) -> bool:
        """
        Heuristischer Test, ob die Form exakt ist: ω = dη für eine (k-1)-Form η.

        Eine k-Form ω heißt exakt, wenn es eine (k-1)-Form η gibt mit ω = dη.
        Für k = 0: Nur die Nullform ist exakt (es gibt keine (-1)-Formen).
        Für k = 1 auf ℝⁿ (sternförmig): ω = Σ fᵢ dxᵢ ist exakt genau dann, wenn
            ∂fᵢ/∂xⱼ = ∂fⱼ/∂xᵢ für alle i, j (Integrabilitätsbedingung).

        Diese Implementierung prüft die Integrabilitätsbedingung für 1-Formen
        und gibt für andere Grade None zurück (unbekannt).

        @return: True wenn exakt (oder heuristisch exakt), False wenn nicht exakt,
                 None wenn unbekannt (höhere Grade)
        @lastModified 2026-03-10
        """
        # 0-Formen: nur die Nullform ist exakt
        if self.degree == 0:
            # Prüfen ob alle Koeffizienten Null sind
            for f in self.coefficients.values():
                if simplify(f) != 0:
                    return False
            return True

        # 1-Formen auf ℝⁿ: Integrabilitätsbedingung (Poincaré-Lemma)
        if self.degree == 1:
            n = len(self.variables)
            # Koeffizientenvektor aufbauen: fᵢ = Koeffizient von dxᵢ
            f = {}
            for i in range(n):
                idx = (i,)
                f[i] = self.coefficients.get(idx, sp.Integer(0))

            # Symmetriebedingung: ∂fᵢ/∂xⱼ = ∂fⱼ/∂xᵢ für alle i < j
            for i in range(n):
                for j in range(i + 1, n):
                    dfi_dxj = simplify(diff(f[i], self.variables[j]))
                    dfj_dxi = simplify(diff(f[j], self.variables[i]))
                    if simplify(dfi_dxj - dfj_dxi) != 0:
                        return False
            return True

        # Für höhere Grade: geschlossen ist notwendige Bedingung, als Heuristik verwenden
        # (gilt exakt auf sternförmigen Gebieten via Poincaré-Lemma)
        return self.is_closed()

    def __repr__(self) -> str:
        """
        Lesbare Darstellung der Differentialform.

        @return: Beschreibungsstring
        @lastModified 2026-03-10
        """
        return (
            f"DifferentialForm(degree={self.degree}, "
            f"num_terms={len(self.coefficients)}, "
            f"vars={[str(v) for v in self.variables]})"
        )


def _sign_and_sort(indices: List[int]) -> Tuple[int, List[int]]:
    """
    Berechnet das Vorzeichen der Permutation, die eine Liste von Indizes sortiert,
    und gibt die sortierte Liste zurück.

    Das Vorzeichen ist +1 (gerade Permutation) oder -1 (ungerade Permutation).
    Bei doppelten Indizes wird das Vorzeichen auf 0 gesetzt (Schiefsymmetrie → 0).

    Algorithmus: Bubble-Sort, Vorzeichen wechselt bei jeder Transposition.

    @param indices: Liste ganzer Zahlen (Indizes des Multi-Index)
    @return: (sign, sorted_indices) – Vorzeichen (+1, -1 oder 0) und sortierte Liste
    @lastModified 2026-03-10
    """
    lst = list(indices)
    n = len(lst)
    sign = 1

    # Bubble-Sort mit Vorzeichenverfolgung
    for i in range(n):
        for j in range(n - 1 - i):
            if lst[j] > lst[j + 1]:
                # Transposition: Vorzeichen umkehren
                lst[j], lst[j + 1] = lst[j + 1], lst[j]
                sign *= -1
            elif lst[j] == lst[j + 1]:
                # Doppelter Index: Form ist Null (Schiefsymmetrie)
                return 0, lst

    return sign, lst


# ===========================================================================
# 3. DE RHAM-KOHOMOLOGIE
# ===========================================================================

class DeRhamCohomology:
    """
    Berechnet de Rham-Kohomologie-Invarianten.

    Die k-te de Rham-Kohomologiegruppe ist definiert als:
        H^k_dR(M) = ker(d : Ω^k → Ω^{k+1}) / im(d : Ω^{k-1} → Ω^k)
                  = {geschlossene k-Formen} / {exakte k-Formen}

    de Rham's Theorem besagt, dass H^k_dR(M) ≅ H^k_sing(M; ℝ) (singuläre Kohomologie).
    Die Dimensionen b_k = dim(H^k_dR(M)) heißen Betti-Zahlen.

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def closed_forms(degree: int, forms: List['DifferentialForm']) -> List['DifferentialForm']:
        """
        Filtert geschlossene k-Formen aus einer Liste von Formen.

        Eine Form ω ist geschlossen, wenn dω = 0.

        @param degree: Grad k der zu filternden Formen
        @param forms: Liste von DifferentialForm-Objekten
        @return: Liste der geschlossenen k-Formen
        @lastModified 2026-03-10
        """
        # Nur Formen des richtigen Grades und mit dω = 0 behalten
        return [omega for omega in forms
                if omega.degree == degree and omega.is_closed()]

    @staticmethod
    def exact_forms(degree: int, forms: List['DifferentialForm']) -> List['DifferentialForm']:
        """
        Filtert exakte k-Formen aus einer Liste von Formen.

        Eine Form ω ist exakt, wenn es η gibt mit ω = dη.

        @param degree: Grad k der zu filternden Formen
        @param forms: Liste von DifferentialForm-Objekten
        @return: Liste der exakten k-Formen
        @lastModified 2026-03-10
        """
        # Nur Formen des richtigen Grades und mit ω = dη behalten
        return [omega for omega in forms
                if omega.degree == degree and omega.is_exact()]

    @staticmethod
    def betti_numbers_sphere(n: int) -> List[int]:
        """
        Gibt die Betti-Zahlen der n-Sphäre S^n zurück.

        De Rham-Kohomologie der n-Sphäre:
            H^k(S^n) = ℝ  für k = 0 und k = n
            H^k(S^n) = 0  sonst

        Also: b_0 = 1, b_n = 1, alle anderen b_k = 0.
        Für n = 0 (zwei Punkte): b_0 = 2.

        @param n: Dimension der Sphäre (n ≥ 0)
        @return: Liste [b_0, b_1, …, b_n] der Betti-Zahlen
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"Dimension der Sphäre muss ≥ 0 sein, erhalten: {n}")

        # S^0 = zwei Punkte: b_0 = 2
        if n == 0:
            return [2]

        # S^n für n ≥ 1: b_0 = b_n = 1, Rest = 0
        betti = [0] * (n + 1)
        betti[0] = 1   # Zusammenhangkomponente
        betti[n] = 1   # Fundamentalklasse
        return betti

    @staticmethod
    def betti_numbers_torus(n: int) -> List[int]:
        """
        Gibt die Betti-Zahlen des n-Torus T^n = S^1 × … × S^1 zurück.

        De Rham-Kohomologie des n-Torus:
            H^k(T^n) = ℝ^{C(n,k)}
        wobei C(n,k) = n!/(k!(n-k)!) der Binomialkoeffizient ist.

        Also: b_k = C(n,k).

        Beispiele:
            T^1 = S^1: (1, 1)
            T^2: (1, 2, 1)
            T^3: (1, 3, 3, 1)

        @param n: Dimension des Torus (n ≥ 0)
        @return: Liste [b_0, b_1, …, b_n] der Betti-Zahlen
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"Dimension des Torus muss ≥ 0 sein, erhalten: {n}")

        import math
        # b_k = C(n,k) (Binomialkoeffizient)
        betti = [math.comb(n, k) for k in range(n + 1)]
        return betti

    @staticmethod
    def euler_characteristic(betti_numbers: List[int]) -> int:
        """
        Berechnet die Euler-Charakteristik aus den Betti-Zahlen.

        Die Euler-Charakteristik ist definiert als:
            χ(M) = Σ_{k=0}^{n} (-1)^k · b_k

        Beispiele:
            χ(S^2) = b_0 - b_1 + b_2 = 1 - 0 + 1 = 2
            χ(T^2) = b_0 - b_1 + b_2 = 1 - 2 + 1 = 0

        @param betti_numbers: Liste [b_0, b_1, …, b_n] der Betti-Zahlen
        @return: Euler-Charakteristik χ(M)
        @lastModified 2026-03-10
        """
        # Alternierende Summe der Betti-Zahlen
        chi = sum(
            (-1)**k * bk
            for k, bk in enumerate(betti_numbers)
        )
        return chi


# ===========================================================================
# 4. MORSE-THEORIE
# ===========================================================================

class MorseTheory:
    """
    Implementiert Konzepte der Morse-Theorie.

    Die Morse-Theorie verbindet die Analysis einer glatten Funktion f : M → ℝ
    mit der Topologie der Mannigfaltigkeit M.

    Ein Punkt p ∈ M heißt kritisch, wenn ∇f(p) = 0.
    Der Morse-Index μ(p) ist die Anzahl der negativen Eigenwerte der Hesse-Matrix.

    f heißt Morse-Funktion, wenn alle kritischen Punkte nicht-degeneriert sind
    (d.h. die Hesse-Matrix ist an jedem kritischen Punkt invertierbar).

    Morse-Ungleichungen:
        c_k ≥ b_k    (schwach)
        Σ_k (-1)^k c_k = χ(M)   (Euler-Charakteristik)

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def morse_function_critical_points(
        f_expr,
        variables: List
    ) -> List[Dict]:
        """
        Findet alle kritischen Punkte der Funktion f.

        Ein Punkt p ist kritisch, wenn der Gradient ∇f(p) = 0, d.h.
        alle partiellen Ableitungen gleichzeitig verschwinden:
            ∂f/∂x₁(p) = … = ∂f/∂xₙ(p) = 0

        @param f_expr: SymPy-Ausdruck für f
        @param variables: Liste der SymPy-Symbole
        @return: Liste von Dicts {var: Wert, …} (kritische Punkte)
        @lastModified 2026-03-10
        """
        # Gradient berechnen: alle partiellen Ableitungen
        gradient = [diff(f_expr, var) for var in variables]

        # Gleichungssystem ∇f = 0 lösen
        try:
            solutions = solve(gradient, variables, dict=True)
        except Exception:
            # Bei unlösbaren Systemen: leere Liste zurückgeben
            solutions = []

        return solutions

    @staticmethod
    def morse_index(
        f_expr,
        critical_point: Dict,
        variables: List
    ) -> int:
        """
        Berechnet den Morse-Index eines kritischen Punktes.

        Der Morse-Index μ(p) ist die Anzahl der negativen Eigenwerte
        der Hesse-Matrix H_f(p), wobei:
            H_f(p)_{ij} = ∂²f/(∂xᵢ∂xⱼ)(p)

        Er gibt die Anzahl der "absteigenden Richtungen" am kritischen Punkt an:
            - Index 0: Lokales Minimum
            - Index n: Lokales Maximum
            - 0 < Index < n: Sattelpunkt

        @param f_expr: SymPy-Ausdruck für f
        @param critical_point: Dict {var: Wert, …} des kritischen Punktes
        @param variables: Liste der SymPy-Symbole
        @return: Morse-Index (Anzahl negativer Eigenwerte der Hesse-Matrix)
        @raises ValueError: Wenn der Punkt degeneriert ist (det(H) = 0)
        @lastModified 2026-03-10
        """
        # Hesse-Matrix symbolisch berechnen
        H = hessian(f_expr, variables)

        # Hesse-Matrix am kritischen Punkt auswerten
        H_at_p = H.subs(critical_point)
        H_at_p_simplified = H_at_p.applyfunc(simplify)

        # Determinante prüfen (Nicht-Degeneriertheit)
        det_H = H_at_p_simplified.det()
        if simplify(det_H) == 0:
            raise ValueError(
                f"Kritischer Punkt {critical_point} ist degeneriert (det(H) = 0). "
                "Kein Morse-Index definiert."
            )

        # Eigenwerte berechnen
        eigenvals = H_at_p_simplified.eigenvals()

        # Anzahl negativer Eigenwerte zählen (mit Vielfachheit)
        morse_index = 0
        for eigenval, multiplicity in eigenvals.items():
            # Eigenwert numerisch auswerten
            ev_num = complex(eigenval.evalf())
            if ev_num.real < 0:
                morse_index += multiplicity

        return morse_index

    @staticmethod
    def is_morse_function(f_expr, variables: List) -> bool:
        """
        Prüft, ob f eine Morse-Funktion ist.

        f ist eine Morse-Funktion, wenn alle kritischen Punkte nicht-degeneriert sind,
        d.h. die Hesse-Matrix H_f(p) an jedem kritischen Punkt p invertierbar ist.

        @param f_expr: SymPy-Ausdruck für f
        @param variables: Liste der SymPy-Symbole
        @return: True wenn f eine Morse-Funktion ist, False wenn degenerierte kritische Punkte existieren
        @lastModified 2026-03-10
        """
        # Kritische Punkte finden
        critical_pts = MorseTheory.morse_function_critical_points(f_expr, variables)

        if not critical_pts:
            # Keine kritischen Punkte: trivialerweise Morse-Funktion
            return True

        # Hesse-Matrix symbolisch aufbauen
        H = hessian(f_expr, variables)

        # Jeden kritischen Punkt auf Nicht-Degeneriertheit prüfen
        for pt in critical_pts:
            H_at_p = H.subs(pt).applyfunc(simplify)
            det_H = simplify(H_at_p.det())
            if det_H == 0:
                # Degenerierter kritischer Punkt gefunden
                return False

        return True

    @staticmethod
    def morse_inequalities_check(
        betti_numbers: List[int],
        morse_data: Dict[int, int]
    ) -> Dict[str, Any]:
        """
        Überprüft die Morse-Ungleichungen.

        Morse-Ungleichungen besagen:
            c_k ≥ b_k  für alle k  (schwache Ungleichung)
            Σ_k (-1)^k c_k = χ(M)  (Euler-Charakteristik)

        wobei c_k die Anzahl der kritischen Punkte mit Morse-Index k ist
        und b_k die k-te Betti-Zahl.

        @param betti_numbers: Liste [b_0, b_1, …, b_n] der Betti-Zahlen
        @param morse_data: Dict {index: count} – Anzahl krit. Punkte pro Morse-Index
        @return: Dict mit 'weak_satisfied' (bool), 'euler_ok' (bool),
                 'euler_char' (int), 'alternating_sum_c' (int)
        @lastModified 2026-03-10
        """
        n = len(betti_numbers) - 1

        # Schwache Morse-Ungleichungen prüfen: c_k ≥ b_k
        weak_satisfied = True
        for k, bk in enumerate(betti_numbers):
            ck = morse_data.get(k, 0)
            if ck < bk:
                weak_satisfied = False
                break

        # Euler-Charakteristik über Betti-Zahlen
        chi = sum((-1)**k * bk for k, bk in enumerate(betti_numbers))

        # Alternierende Summe der kritischen Punkte
        alt_sum_c = sum(
            (-1)**k * morse_data.get(k, 0)
            for k in range(n + 1)
        )

        # Euler-Bedingung: Σ (-1)^k c_k = χ(M)
        euler_ok = (alt_sum_c == chi)

        return {
            'weak_satisfied': weak_satisfied,
            'euler_ok': euler_ok,
            'euler_char': chi,
            'alternating_sum_c': alt_sum_c
        }

    @staticmethod
    def handle_decomposition_demo(f_expr, variables: List) -> Dict[str, Any]:
        """
        Demonstriert die Henkelzerlegung einer Mannigfaltigkeit via Morse-Theorie.

        Nach Morse wird beim Überqueren eines kritischen Wertes f(p) mit Morse-Index k
        ein k-dimensionaler Henkel (Handle) angeklebt. Die Topologie entsteht durch
        sukzessives Ankleben von Henkeln, beginnend beim globalen Minimum.

        @param f_expr: SymPy-Ausdruck der Morse-Funktion
        @param variables: Liste der SymPy-Symbole
        @return: Dict mit 'critical_points', 'indices', 'handles', 'is_morse'
        @lastModified 2026-03-10
        """
        # Kritische Punkte bestimmen
        critical_pts = MorseTheory.morse_function_critical_points(f_expr, variables)

        # Für jeden kritischen Punkt Morse-Index berechnen
        handles = []
        indices = []

        for pt in critical_pts:
            try:
                idx = MorseTheory.morse_index(f_expr, pt, variables)
                indices.append(idx)
                # Kritischen Wert berechnen
                f_val = f_expr.subs(pt)
                handles.append({
                    'point': pt,
                    'critical_value': float(f_val.evalf()),
                    'morse_index': idx,
                    'handle_type': f"{idx}-Henkel"
                })
            except ValueError:
                # Degenerierter Punkt
                handles.append({
                    'point': pt,
                    'critical_value': None,
                    'morse_index': None,
                    'handle_type': 'degeneriert'
                })

        # Henkel nach kritischem Wert sortieren (Reihenfolge des Anklebens)
        handles_sorted = sorted(
            [h for h in handles if h['critical_value'] is not None],
            key=lambda h: h['critical_value']
        )

        return {
            'critical_points': critical_pts,
            'indices': indices,
            'handles': handles_sorted,
            'is_morse': MorseTheory.is_morse_function(f_expr, variables)
        }


# ===========================================================================
# 5. TRANSVERSALITÄTSTHEORIE
# ===========================================================================

class TransversalityTheory:
    """
    Implementiert Konzepte der Transversalitätstheorie.

    Zwei glatte Untermannigfaltigkeiten M, N ⊂ A heißen transversal (M ⋔ N),
    wenn an jedem Schnittpunkt p ∈ M ∩ N gilt:
        T_pM + T_pN = T_pA
    (die Tangentialräume spannen zusammen den Gesamtraum auf).

    Dies ist äquivalent zu:
        codim(M) + codim(N) ≤ dim(A)
    oder:    dim(M) + dim(N) ≥ dim(A)

    Sard's Theorem: Die Menge der kritischen Werte einer glatten Abbildung
    f : M → N hat Maß 0 in N.

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def are_transverse(
        M_codim: int,
        N_codim: int,
        ambient_dim: int
    ) -> bool:
        """
        Prüft die Transversalitätsbedingung für zwei Untermannigfaltigkeiten.

        M und N sind transversal in A (dim A = ambient_dim), wenn:
            codim(M) + codim(N) ≤ dim(A)
        Äquivalent: dim(M) + dim(N) ≥ dim(A)

        @param M_codim: Kodimension von M in A (= dim(A) - dim(M))
        @param N_codim: Kodimension von N in A (= dim(A) - dim(N))
        @param ambient_dim: Dimension des umgebenden Raums A
        @return: True wenn Transversalitätsbedingung erfüllt ist
        @lastModified 2026-03-10
        """
        # Transversalität: Summe der Kodimensionen ≤ ambient_dim
        return M_codim + N_codim <= ambient_dim

    @staticmethod
    def transverse_intersection_dimension(
        m_dim: int,
        n_dim: int,
        ambient_dim: int
    ) -> int:
        """
        Berechnet die erwartete Dimension des transversalen Schnitts M ∩ N.

        Bei transversalem Schnitt gilt:
            dim(M ∩ N) = dim(M) + dim(N) - dim(A)

        Falls dim(M) + dim(N) < dim(A), ist der Schnitt generisch leer (leere Menge
        hat formal Dimension -∞, hier -1 als Sentinel).

        @param m_dim: Dimension von M
        @param n_dim: Dimension von N
        @param ambient_dim: Dimension des umgebenden Raums A
        @return: Dimension des Schnitts, oder -1 wenn generisch leer
        @lastModified 2026-03-10
        """
        # Erwartete Schnittdimension
        expected_dim = m_dim + n_dim - ambient_dim
        # Negativ → generisch leer
        return expected_dim

    @staticmethod
    def sard_theorem_demo(
        f: Callable,
        df: Callable,
        domain_points: np.ndarray
    ) -> Dict[str, Any]:
        """
        Demonstriert das Sard'sche Theorem numerisch.

        Sard's Theorem: Sei f : M → N eine glatte Abbildung. Dann hat die Menge
            {f(p) : p kritischer Punkt von f} ⊂ N
        Maß 0 in N (die kritischen Werte sind "selten").

        Ein Punkt p heißt regulär, wenn df(p) : T_pM → T_{f(p)}N surjektiv ist.
        Sonst heißt p kritisch.

        @param f: Abbildung f : ℝⁿ → ℝᵐ als Python-Funktion
        @param df: Jacobi-Matrix df(p) als Funktion ℝⁿ → ℝ^{m×n}
        @param domain_points: Numpy-Array der Testpunkte (n_points × n)
        @return: Dict mit 'critical_points', 'regular_points',
                 'critical_values', 'ratio_critical' (Anteil kritischer Punkte)
        @lastModified 2026-03-10
        """
        critical_pts = []
        regular_pts = []
        critical_values = []

        for pt in domain_points:
            # Jacobi-Matrix an diesem Punkt auswerten
            J = np.array(df(pt), dtype=float)

            # Rang der Jacobi-Matrix bestimmen
            rank = np.linalg.matrix_rank(J)

            # Dimension des Zielraums (Anzahl Zeilen von J)
            m = J.shape[0] if J.ndim > 1 else 1

            # Kritisch ↔ Rang < m (nicht surjektiv)
            if rank < m:
                critical_pts.append(pt)
                critical_values.append(f(pt))
            else:
                regular_pts.append(pt)

        total = len(domain_points)
        ratio = len(critical_pts) / total if total > 0 else 0.0

        return {
            'critical_points': critical_pts,
            'regular_points': regular_pts,
            'critical_values': critical_values,
            'ratio_critical': ratio,
            'sard_verified': ratio < 0.5  # Heuristik: kritische Punkte sind selten
        }


# ===========================================================================
# 6. VEKTORBÜNDEL
# ===========================================================================

class VectorBundle:
    """
    Repräsentiert ein Vektorbündel E → B.

    Ein Vektorbündel über einer Basis B mit typischer Faser F ≅ ℝ^r besteht aus:
    - Totalraum E
    - Basisraum B (Dimension base_dim)
    - Projektion π : E → B
    - Typische Faser ℝ^r (Dimension fiber_dim = Rang r)

    Der Totalraum hat Dimension dim(E) = dim(B) + r.

    Beispiele:
    - Tangentbündel TM → M: fiber = T_pM ≅ ℝⁿ, rang = n
    - Triviales Bündel: E = B × ℝ^r (globales Produkt)
    - Möbiusband: nicht-triviales ℝ¹-Bündel über S¹

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    def __init__(
        self,
        base_dim: int,
        fiber_dim: int,
        name: str = 'E',
        trivial: bool = True
    ) -> None:
        """
        Initialisiert ein Vektorbündel.

        @param base_dim: Dimension des Basisraums B
        @param fiber_dim: Dimension der Faser (= Rang des Bündels)
        @param name: Bezeichner des Bündels (Standard: 'E')
        @param trivial: True wenn triviales Bündel E = B × ℝ^r (Standard: True)
        @raises ValueError: Wenn base_dim oder fiber_dim negativ
        @lastModified 2026-03-10
        """
        if base_dim < 0:
            raise ValueError(f"Basisdimension muss ≥ 0 sein, erhalten: {base_dim}")
        if fiber_dim < 0:
            raise ValueError(f"Faserdimension muss ≥ 0 sein, erhalten: {fiber_dim}")

        # Dimensionen der Basis und Faser speichern
        self.base_dim = base_dim
        self.fiber_dim = fiber_dim
        # Name des Bündels
        self.name = name
        # Trivialitätsstatus
        self._trivial = trivial

    def total_space_dim(self) -> int:
        """
        Berechnet die Dimension des Totalraums E.

        dim(E) = dim(B) + dim(Faser) = base_dim + fiber_dim

        @return: Dimension des Totalraums
        @lastModified 2026-03-10
        """
        return self.base_dim + self.fiber_dim

    def rank(self) -> int:
        """
        Gibt den Rang des Vektorbündels zurück.

        Der Rang eines Vektorbündels ist die Dimension der typischen Faser:
        rank(E) = dim(Faser) = fiber_dim.

        @return: Rang des Bündels
        @lastModified 2026-03-10
        """
        return self.fiber_dim

    def is_trivial(self) -> bool:
        """
        Gibt zurück, ob das Bündel trivial ist.

        Ein Vektorbündel E → B heißt trivial, wenn es einen globalen
        Bündelisomorphismus φ : E → B × ℝ^r gibt (globale Trivialisierung).

        @return: True wenn triviales Bündel
        @lastModified 2026-03-10
        """
        return self._trivial

    def __repr__(self) -> str:
        """
        Lesbare Darstellung des Vektorbündels.

        @return: Beschreibungsstring
        @lastModified 2026-03-10
        """
        trivial_str = "trivial" if self._trivial else "nicht-trivial"
        return (
            f"VectorBundle(name='{self.name}', "
            f"base_dim={self.base_dim}, fiber_dim={self.fiber_dim}, "
            f"rank={self.fiber_dim}, {trivial_str})"
        )


# ===========================================================================
# 7. CHARAKTERISTISCHE KLASSEN
# ===========================================================================

class CharacteristicClasses:
    """
    Berechnet charakteristische Klassen von Vektorbündeln (algebraisch/demo).

    Charakteristische Klassen sind kanonische Kohomologieklassen, die einem
    Vektorbündel zugeordnet werden und topologische Invarianten darstellen.

    - Chern-Klassen c_k(E) ∈ H^{2k}(M; ℤ) für komplexe Bündel
    - Pontryagin-Klassen p_k(E) ∈ H^{4k}(M; ℤ) für reelle Bündel
    - Euler-Klasse e(E) ∈ H^n(M; ℤ) für orientierte Bündel vom Rang n
    - Stiefel-Whitney-Klassen w_k(E) ∈ H^k(M; ℤ/2ℤ) für reelle Bündel

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def chern_class_demo(n: int) -> Dict[str, Any]:
        """
        Demonstriert die Chern-Klassen des Tangentbündels TCP^n des komplexen
        projektiven Raums ℂP^n.

        Für TCP^n gilt die Formel:
            c(TCP^n) = (1 + H)^{n+1}  in H*(ℂP^n; ℤ) = ℤ[H]/(H^{n+1})
        wobei H der erzeugende Hyperebenenklasse entspricht.

        Also: c_k(TCP^n) = C(n+1, k) · H^k

        @param n: Dimension von ℂP^n (n ≥ 0)
        @return: Dict mit 'chern_numbers' (Liste der Chern-Zahlen C(n+1,k)),
                 'total_chern_class' (symbolisch), 'top_chern_number' (= euler_characteristic)
        @lastModified 2026-03-10
        """
        import math
        # c_k = C(n+1, k) für k = 0, 1, …, n
        chern_numbers = [math.comb(n + 1, k) for k in range(n + 1)]

        # Euler-Charakteristik von ℂP^n ist n+1
        euler_char = n + 1

        return {
            'chern_numbers': chern_numbers,
            'manifold': f'CP^{n}',
            'top_chern_number': euler_char,
            'euler_characteristic': euler_char,
            'formula': 'c_k(TCP^n) = C(n+1, k) * H^k'
        }

    @staticmethod
    def pontryagin_class_demo(n: int) -> Dict[str, Any]:
        """
        Demonstriert die Pontryagin-Klassen des Tangentbündels TRP^n von RP^n.

        Für das Tangentbündel der n-Sphäre S^n gilt (da TS^n ⊕ trivial = trivial^{n+1}):
            p(TS^n) = 1  (alle Pontryagin-Klassen verschwinden)

        Für allgemeine Mannigfaltigkeiten: p_k ∈ H^{4k}(M; ℤ).
        Hier wird die formale Berechnung via Chern-Klassen demonstriert:
            p_k(E_ℝ) = (-1)^k c_{2k}(E_ℝ ⊗ ℂ)

        @param n: Dimension der Mannigfaltigkeit
        @return: Dict mit 'pontryagin_classes', 'description'
        @lastModified 2026-03-10
        """
        # Pontryagin-Klassen leben in H^{4k}, also nur für 4k ≤ n relevant
        max_k = n // 4

        # Für S^n: alle Pontryagin-Klassen verschwinden
        pontryagin = {k: 0 for k in range(1, max_k + 1)}
        # p_0 = 1 immer
        pontryagin[0] = 1

        return {
            'pontryagin_classes': pontryagin,
            'manifold': f'S^{n}',
            'description': (
                f'Für S^{n}: p_k(TS^{n}) = 0 für k ≥ 1, '
                'da TS^n ⊕ trivial = trivial^{n+1} (Stably trivial)'
            ),
            'max_non_trivial_index': max_k
        }

    @staticmethod
    def euler_class_demo(dimension: int) -> Dict[str, Any]:
        """
        Demonstriert die Euler-Klasse für das Tangentbündel einer orientieren Mannigfaltigkeit.

        Die Euler-Klasse e(E) ∈ H^r(M; ℤ) eines orientierten reellen Vektorbündels
        vom Rang r über M ist definiert als der primäre Obstruktionsklasse für
        nirgends-verschwindende Schnitte.

        Für das Tangentbündel TM: e(TM) integriert über M ergibt χ(M).

        Für S^n:
            e(TS^n) = χ(S^n) · [S^n]* = (1 + (-1)^n) · [S^n]*

        @param dimension: Dimension n der Mannigfaltigkeit
        @return: Dict mit 'euler_class_value', 'euler_characteristic',
                 'is_zero_section_avoidable'
        @lastModified 2026-03-10
        """
        # Euler-Charakteristik von S^n
        if dimension == 0:
            chi = 2
        else:
            # χ(S^n) = 1 + (-1)^n
            chi = 1 + (-1) ** dimension

        # Nirgends verschwindender Schnitt existiert ↔ Euler-Klasse = 0
        # Für S^n: möglich genau für n ungerade (hairy ball theorem für S^2, S^4, …)
        zero_section_avoidable = (chi == 0)

        return {
            'manifold': f'S^{dimension}',
            'euler_class_value': chi,
            'euler_characteristic': chi,
            'is_zero_section_avoidable': zero_section_avoidable,
            'description': (
                f'e(TS^{dimension}) = χ(S^{dimension}) = {chi}. '
                + ('Nirgends verschwindender Schnitt existiert (Euler-Klasse = 0).'
                   if zero_section_avoidable
                   else 'Kein nirgends verschwindender Schnitt (Satz vom Igel für gerades n).')
            )
        }

    @staticmethod
    def stiefel_whitney_demo(n: int) -> Dict[str, Any]:
        """
        Demonstriert die Stiefel-Whitney-Klassen w_k ∈ H^k(M; ℤ/2ℤ).

        Für das Tangentbündel der reellen projektiven Ebene RP^n gilt:
            w(TRP^n) = (1 + a)^{n+1}  in H*(RP^n; ℤ/2ℤ) = ℤ/2ℤ[a]/(a^{n+1})
        wobei a der erzeugende Klasse in H^1(RP^n; ℤ/2ℤ) ist.

        Also: w_k(TRP^n) = C(n+1, k) mod 2

        w_1 = 0 ↔ M ist orientierbar.
        w₁ = 0 und w₂ = 0 ↔ M trägt Spin-Struktur.

        @param n: Dimension von RP^n
        @return: Dict mit 'stiefel_whitney_classes', 'orientable', 'spin_structure'
        @lastModified 2026-03-10
        """
        import math
        # w_k = C(n+1, k) mod 2
        sw_classes = [math.comb(n + 1, k) % 2 for k in range(n + 1)]

        # Orientierbarkeit: w_1 = 0
        orientable = (sw_classes[1] == 0) if n >= 1 else True

        # Spin-Struktur: w_1 = 0 und w_2 = 0
        spin = (
            (sw_classes[1] == 0 and sw_classes[2] == 0)
            if n >= 2
            else (sw_classes[1] == 0 if n >= 1 else True)
        )

        return {
            'stiefel_whitney_classes': sw_classes,
            'manifold': f'RP^{n}',
            'orientable': orientable,
            'spin_structure': spin,
            'formula': 'w_k(TRP^n) = C(n+1, k) mod 2'
        }


# ===========================================================================
# 8. SPHÄRENTOPOLOGIE
# ===========================================================================

class SphereTopology:
    """
    Repräsentiert die n-dimensionale Sphäre S^n und ihre topologischen Invarianten.

    Die n-Sphäre ist definiert als:
        S^n = {x ∈ ℝ^{n+1} : |x| = 1}

    Sie ist eine kompakte, verbundene, orientierbare n-dimensionale Mannigfaltigkeit.

    @author Michael Fuhrmann
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    # Bekannte Homotopiegruppen πₖ(Sⁿ) für k ≤ n+3
    # Einträge: (n, k) → Gruppenname als String
    # Quellen: Hatcher "Algebraic Topology", Toda "Composition Methods in Homotopy Groups of Spheres"
    _HOMOTOPY_GROUPS = {
        # π_k(S^1)
        (1, 1): 'ℤ', (1, 2): '0', (1, 3): '0', (1, 4): '0',
        # π_k(S^2)
        (2, 2): 'ℤ', (2, 3): 'ℤ', (2, 4): 'ℤ/2ℤ', (2, 5): 'ℤ/2ℤ',
        # π_k(S^3)
        (3, 3): 'ℤ', (3, 4): 'ℤ/2ℤ', (3, 5): 'ℤ/2ℤ', (3, 6): 'ℤ/12ℤ',
        # π_k(S^4)
        (4, 4): 'ℤ', (4, 5): 'ℤ/2ℤ', (4, 6): 'ℤ/2ℤ', (4, 7): 'ℤ ⊕ ℤ/12ℤ',
        # π_k(S^5)
        (5, 5): 'ℤ', (5, 6): 'ℤ/2ℤ', (5, 7): 'ℤ/2ℤ', (5, 8): 'ℤ/24ℤ',
    }

    def __init__(self, n: int) -> None:
        """
        Initialisiert die n-Sphäre S^n.

        @param n: Dimension der Sphäre (n ≥ 0)
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"Dimension der Sphäre muss ≥ 0 sein, erhalten: {n}")
        # Dimension der Sphäre
        self.n = n

    def homotopy_groups_low(self, n: Optional[int] = None) -> Dict[int, str]:
        """
        Gibt bekannte Homotopiegruppen πₖ(Sⁿ) für k ≤ n+3 zurück.

        Homotopiegruppen von Sphären sind notorisch schwer zu berechnen.
        Bekannte Werte (Freudenthal, Toda, u.a.):
            πₙ(Sⁿ) = ℤ  (fundamentale Klasse)
            π_{n+1}(Sⁿ) = ℤ/2ℤ  für n ≥ 3  (Hopf-Abbildung)
            π_k(S^1) = 0  für k ≥ 2  (S^1 ist K(ℤ,1))

        @param n: Dimension der Sphäre (Standard: self.n)
        @return: Dict {k: Gruppenname} für verfügbare k-Werte
        @lastModified 2026-03-10
        """
        sphere_n = n if n is not None else self.n

        result = {}
        # πₖ(S^n) = 0 für k < n (Hurewicz-Theorem)
        for k in range(0, sphere_n):
            result[k] = '0'

        # Bekannte Werte aus Tabelle
        max_k = sphere_n + 3
        for k in range(sphere_n, max_k + 1):
            key = (sphere_n, k)
            if key in self._HOMOTOPY_GROUPS:
                result[k] = self._HOMOTOPY_GROUPS[key]
            elif k == sphere_n:
                # πₙ(Sⁿ) = ℤ immer
                result[k] = 'ℤ'
            else:
                result[k] = 'unbekannt'

        return result

    def is_orientable(self) -> bool:
        """
        Gibt zurück, ob die Sphäre orientierbar ist.

        Sphären S^n sind für alle n ≥ 0 orientierbar, da sie als Hyperflächen
        im ℝ^{n+1} mit dem nach außen zeigenden Normalenvektor orientiert werden können.

        @return: True (Sphären sind immer orientierbar)
        @lastModified 2026-03-10
        """
        # Alle Sphären sind orientierbar
        return True

    def homology_groups(self, n: Optional[int] = None) -> Dict[int, str]:
        """
        Gibt die singulären Homologiegruppen H_k(S^n; ℤ) zurück.

        Bekannte Werte (aus universellen Koeffizientensatz und CW-Struktur):
            H_0(S^n) = ℤ  (verbunden)
            H_n(S^n) = ℤ  (Fundamentalklasse)
            H_k(S^n) = 0  sonst

        Für S^0 = zwei Punkte: H_0 = ℤ ⊕ ℤ = ℤ².

        @param n: Dimension der Sphäre (Standard: self.n)
        @return: Dict {k: Gruppenname} für k = 0, 1, …, n
        @lastModified 2026-03-10
        """
        sphere_n = n if n is not None else self.n

        homology = {}

        if sphere_n == 0:
            # S^0 = {-1, +1}: zwei Zusammenhangskomponenten
            homology[0] = 'ℤ²'
            return homology

        # S^n für n ≥ 1
        for k in range(sphere_n + 1):
            if k == 0:
                homology[k] = 'ℤ'    # verbunden
            elif k == sphere_n:
                homology[k] = 'ℤ'    # Fundamentalklasse
            else:
                homology[k] = '0'    # alle anderen verschwinden

        return homology

    def __repr__(self) -> str:
        """
        Lesbare Darstellung der Sphäre.

        @return: Beschreibungsstring
        @lastModified 2026-03-10
        """
        return f"SphereTopology(n={self.n}, S^{self.n})"


# ===========================================================================
# 9. STANDALONE-FUNKTIONEN
# ===========================================================================

def degree_of_map(f_values: np.ndarray, g_values: np.ndarray) -> int:
    """
    Schätzt den topologischen Abbildungsgrad einer Abbildung numerisch.

    Der Abbildungsgrad deg(f) einer glatten Abbildung f : S^n → S^n ist
    die algebraische Anzahl der Urbilder eines regulären Wertes.

    Numerisch: Für Abbildungen S^1 → S^1 wird der Grad über das Umlaufintegral
    der Winkeländerung berechnet:
        deg(f) = (1/2π) · ∮ dθ

    @param f_values: Werte der ersten Abbildung (Array von Punkten auf S^1, als Winkel in [0, 2π))
    @param g_values: Werte der zweiten Abbildung (Referenz, gleichmäßig verteilt)
    @return: Abbildungsgrad (ganzzahlig)
    @lastModified 2026-03-10
    """
    # Winkelfolge der Abbildungswerte
    angles = np.array(f_values, dtype=float)

    # Gesamte Winkeländerung (Umlaufzahl)
    d_angles = np.diff(np.unwrap(angles))
    total_winding = np.sum(d_angles)

    # Grad = Gesamtwindung / 2π (gerundet)
    degree = int(np.round(total_winding / (2 * np.pi)))
    return degree


def hairy_ball_theorem_demo() -> Dict[str, Any]:
    """
    Demonstriert den Satz vom Igel (Hairy Ball Theorem).

    Satz vom Igel (Poincaré, 1885; Brouwer, 1912):
        Auf S^{2n} (geraddimensionale Sphären) gibt es kein nirgends
        verschwindendes stetiges Vektorfeld.

    Intuition: Man kann einen Igel nicht so kämmen, dass kein Haar steht.

    Äquivalenz: Ein nirgends verschwindendes Vektorfeld existiert genau dann,
    wenn die Euler-Charakteristik χ(M) = 0 ist.
        χ(S^2) = 2 ≠ 0 → kein nirgends verschwindendes Vektorfeld auf S^2
        χ(S^1) = 0 → Vektorfeld existiert (z.B. konstante Tangentialrichtung)

    @return: Dict mit Theorem-Informationen und Beispielen
    @lastModified 2026-03-10
    """
    return {
        'theorem': 'Satz vom Igel (Hairy Ball Theorem)',
        'statement': (
            'Auf S^{2n} (geraddimensionale Sphären) gibt es kein '
            'nirgends verschwindendes stetiges Vektorfeld.'
        ),
        'examples': {
            'S^1': {'euler_char': 0, 'has_nowhere_vanishing_field': True,
                    'example': 'Tangentialvektor (y, -x) auf S^1 verschwindet nie'},
            'S^2': {'euler_char': 2, 'has_nowhere_vanishing_field': False,
                    'reason': 'χ(S^2) = 2 ≠ 0'},
            'S^3': {'euler_char': 0, 'has_nowhere_vanishing_field': True,
                    'example': 'Hopf-Faserung liefert globalen Rahmen'},
            'S^4': {'euler_char': 2, 'has_nowhere_vanishing_field': False,
                    'reason': 'χ(S^4) = 2 ≠ 0'},
        },
        'criterion': 'Nirgends verschwindendes Vektorfeld existiert ↔ χ(M) = 0',
        'poincare_hopf': 'Summe der Indizes aller Singularitäten = χ(M)'
    }


def poincare_hopf_theorem_demo(
    euler_char: int,
    vector_field_index_sum: int
) -> Dict[str, Any]:
    """
    Demonstriert den Satz von Poincaré-Hopf.

    Satz von Poincaré-Hopf:
        Sei M eine kompakte orientierbare Mannigfaltigkeit ohne Rand,
        V ein Vektorfeld auf M mit isolierten Nullstellen p₁, …, pₖ.
        Dann gilt:
            Σᵢ ind(V, pᵢ) = χ(M)
        wobei ind(V, p) der Index des Vektorfeldes V an der Nullstelle p ist.

    @param euler_char: Euler-Charakteristik χ(M) der Mannigfaltigkeit
    @param vector_field_index_sum: Summe der Vektorfeldindizes Σᵢ ind(V, pᵢ)
    @return: Dict mit Verifikation und Interpretation
    @lastModified 2026-03-10
    """
    # Prüfen: Summe der Indizes = Euler-Charakteristik
    theorem_satisfied = (vector_field_index_sum == euler_char)

    return {
        'theorem': 'Satz von Poincaré-Hopf',
        'euler_characteristic': euler_char,
        'index_sum': vector_field_index_sum,
        'theorem_satisfied': theorem_satisfied,
        'statement': f'Σ ind(V, pᵢ) = χ(M) = {euler_char}',
        'verification': (
            f'Indexsumme = {vector_field_index_sum}, '
            f'Euler-Charakteristik = {euler_char}: '
            + ('Gleichung erfüllt ✓' if theorem_satisfied else 'Gleichung NICHT erfüllt ✗')
        ),
        'corollary': (
            'Nirgends verschwindendes Vektorfeld existiert'
            if euler_char == 0
            else f'Jedes Vektorfeld hat mindestens eine Nullstelle (χ ≠ 0)'
        )
    }


def whitney_embedding_dimension(manifold_dim: int) -> int:
    """
    Berechnet die von Whitney garantierte Einbettungsdimension.

    Whitney-Einbettungssatz (1944):
        Jede glatte n-dimensionale Mannigfaltigkeit M kann glatt eingebettet
        werden in den ℝ^{2n} (starke Form).

        Die schwache Form garantiert eine Einbettung in ℝ^{2n+1}.

    @param manifold_dim: Dimension n der Mannigfaltigkeit
    @return: Dimension 2n des kleinstmöglichen Einbettungsraums (Whitney stark)
    @raises ValueError: Wenn manifold_dim < 0
    @lastModified 2026-03-10
    """
    if manifold_dim < 0:
        raise ValueError(f"Mannigfaltigkeitsdimension muss ≥ 0 sein, erhalten: {manifold_dim}")

    # Whitney (stark): Einbettung in ℝ^{2n}
    return 2 * manifold_dim


def compute_jacobian(f_exprs: List, variables: List) -> Matrix:
    """
    Berechnet die Jacobi-Matrix einer vektorwertigen Abbildung f : ℝⁿ → ℝᵐ.

    Die Jacobi-Matrix (auch Funktionalmatrix) ist:
        J_f(x)_{ij} = ∂fᵢ/∂xⱼ

    Sie ist eine m×n-Matrix und stellt die beste lineare Approximation
    von f an einem Punkt dar: f(x + h) ≈ f(x) + J_f(x) · h.

    @param f_exprs: Liste von m SymPy-Ausdrücken [f₁, …, fₘ]
    @param variables: Liste von n SymPy-Symbolen [x₁, …, xₙ]
    @return: m×n SymPy-Matrix (Jacobi-Matrix)
    @lastModified 2026-03-10
    """
    m = len(f_exprs)
    n = len(variables)

    # Jacobi-Matrix Einträge: J[i,j] = ∂fᵢ/∂xⱼ
    J = Matrix(m, n, lambda i, j: diff(f_exprs[i], variables[j]))

    return J


def implicit_function_theorem_check(
    F_exprs: List,
    x_vars: List,
    y_vars: List,
    point: Dict
) -> Dict[str, Any]:
    """
    Überprüft die Voraussetzungen des Satzes über implizite Funktionen.

    Satz über implizite Funktionen (Dini):
        Sei F : ℝⁿ × ℝᵐ → ℝᵐ eine C¹-Abbildung mit F(x₀, y₀) = 0.
        Wenn die partielle Jacobi-Matrix J_y F(x₀, y₀) (m×m-Matrix)
        invertierbar ist, dann gibt es eine eindeutige C¹-Funktion
        g : U → V mit F(x, g(x)) = 0 in einer Umgebung U von x₀.

    @param F_exprs: Liste von m SymPy-Ausdrücken für F = (F₁, …, Fₘ)
    @param x_vars: Liste der "unabhängigen" Variablen (x-Variablen)
    @param y_vars: Liste der "abhängigen" Variablen (y-Variablen)
    @param point: Dict {var: Wert} für den Punkt (x₀, y₀)
    @return: Dict mit 'condition_satisfied' (bool),
             'Jy_det' (Determinante der partiellen Jacobi-Matrix),
             'implicit_function_exists' (bool)
    @lastModified 2026-03-10
    """
    m = len(F_exprs)
    if len(y_vars) != m:
        raise ValueError(
            f"Anzahl Gleichungen ({m}) muss mit Anzahl y-Variablen ({len(y_vars)}) übereinstimmen"
        )

    # Partielle Jacobi-Matrix J_y F: m×m-Matrix der Ableitungen nach y
    J_y = Matrix(m, m, lambda i, j: diff(F_exprs[i], y_vars[j]))

    # J_y an gegebenem Punkt auswerten
    J_y_at_point = J_y.subs(point).applyfunc(simplify)

    # Determinante berechnen
    det_Jy = simplify(J_y_at_point.det())

    # Bedingung: det(J_y) ≠ 0
    condition_satisfied = (det_Jy != 0)

    # F(x₀, y₀) auswerten
    F_at_point = [simplify(f.subs(point)) for f in F_exprs]
    f_zero = all(v == 0 for v in F_at_point)

    return {
        'condition_satisfied': bool(condition_satisfied),
        'Jy_det': det_Jy,
        'Jy_matrix': J_y_at_point,
        'F_at_point': F_at_point,
        'F_is_zero': f_zero,
        'implicit_function_exists': bool(condition_satisfied and f_zero),
        'theorem': (
            'Implizite Funktion g mit F(x,g(x))=0 existiert '
            + ('in Umgebung des Punktes ✓' if (condition_satisfied and f_zero)
               else '– Bedingungen nicht erfüllt ✗')
        )
    }
