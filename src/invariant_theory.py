"""
@file invariant_theory.py
@brief Invariantentheorie: Polynominvarianten unter Gruppenwirkungen.
       Reynolds-Operator, Molien-Reihe, elementarsymmetrische Polynome,
       Newton-Identitäten, Hilbert-Basissatz, Diskriminante.
@description
    Dieses Modul implementiert die klassische Invariantentheorie:

    Freie Funktionen:
    - reynolds_operator()                   – Mittelung über Gruppenorbit
    - is_invariant()                        – Prüft G-Invarianz eines Polynoms
    - polynomial_invariants_sn()            – Invarianten unter Sₙ
    - elementary_symmetric_polynomials()    – e₁,...,eₙ in x₁,...,xₙ
    - power_sum_symmetric_polynomials()     – pₖ = Σ xᵢᵏ
    - newton_identities()                   – Rekursion eₖ ↔ pₖ
    - molien_series()                       – Poincaré-Reihe des Invariantenrings
    - hilbert_basis_theorem_demo()          – Erzeuger für S₂-Invariantenring
    - discriminant_invariant()              – Diskriminante als Invariante
    - fundamental_theorem_symmetric_poly()  – Darstellung als Polynom in eₖ

    Mathematische Grundlagen:
    - Reynolds-Operator: R(f) = (1/|G|) Σ_{g∈G} g·f projiziert auf G-invariante Polynome
    - Molien: M(t) = (1/|G|) Σ_{g∈G} 1/det(I - t·g), Koeff. von tᵏ = dim(Inv_k)
    - Hilbert-Basissatz: Jeder Invariantenring k[x₁,...,xₙ]^G ist endlich erzeugt

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sympy as sp
from sympy import symbols, expand, factor, simplify, Rational, Matrix
from itertools import permutations
from typing import Callable, List, Optional
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))


def _apply_permutation(poly_expr: sp.Expr, variable_names: list,
                        perm: tuple) -> sp.Expr:
    """
    Wendet eine Permutation auf die Variablen eines Polynoms an.

    Ersetzt xᵢ durch x_{perm[i]} in poly_expr — simultan, um
    Kettensubstitutionsfehler zu vermeiden (z.B. x1→x2→... wird falsch bei
    sequenzieller Auswertung).

    @param poly_expr      SymPy-Ausdruck in Variablen x₁,...,xₙ
    @param variable_names Liste der Variablennamen ['x1','x2',...]
    @param perm           Tupel der Permutation (0-indiziert)
    @return               Permutiertes Polynom
    """
    syms = sp.symbols(variable_names)
    # Simultane Substitutionsliste: [(xᵢ, x_{perm[i]}), ...]
    substitution = [(syms[i], syms[perm[i]]) for i in range(len(syms))]
    # simultaneous=True verhindert Kettensubstitution (xᵢ→xⱼ→... Kette)
    return poly_expr.subs(substitution, simultaneous=True)


def reynolds_operator(poly_expr: sp.Expr, group_elements: list,
                      action: Callable) -> sp.Expr:
    """
    Reynolds-Operator: Projiziert ein Polynom auf den G-invarianten Unterraum.

    Der Reynolds-Operator ist definiert als:
    R(f) = (1/|G|) · Σ_{g∈G} g·f

    Er ist ein Projektor: R∘R = R und Bild(R) = k[x]^G (G-invariante Polynome).
    Für endliche Gruppen über Körpern mit char(k) ∤ |G| ist R wohldefiniert.

    @param poly_expr      SymPy-Ausdruck (das zu projizierende Polynom)
    @param group_elements Liste der Gruppenelemente g
    @param action         Callable: (poly, g) → g·poly (Gruppenoperation auf Polynomen)
    @return               G-invariantes Polynom R(f)
    """
    n = len(group_elements)
    if n == 0:
        raise ValueError("Leere Gruppe — mindestens ein Element erforderlich")

    # Summe aller Gruppenoperationen auf f
    total = sp.Integer(0)
    for g in group_elements:
        total = total + action(poly_expr, g)

    # Mittelwert = Projektion auf invarianten Unterraum
    result = sp.expand(total / n)
    return result


def is_invariant(poly_expr: sp.Expr, variable_names: list,
                 group_elements: list, action: Callable) -> bool:
    """
    Prüft ob ein Polynom G-invariant ist: g·f = f für alle g∈G.

    @param poly_expr      SymPy-Ausdruck
    @param variable_names Liste der Variablennamen
    @param group_elements Liste der Gruppenelemente
    @param action         Callable: (poly, g) → g·poly
    @return               True wenn f G-invariant ist
    """
    for g in group_elements:
        # Gruppenanwendung und Vereinfachung
        g_f = sp.expand(action(poly_expr, g))
        f_expanded = sp.expand(poly_expr)
        # Prüfe ob Differenz = 0
        if sp.expand(g_f - f_expanded) != 0:
            return False
    return True


def elementary_symmetric_polynomials(n: int) -> list:
    """
    Berechnet die elementarsymmetrischen Polynome e₁,...,eₙ in x₁,...,xₙ.

    Definition:
    e₁ = x₁ + x₂ + ... + xₙ
    e₂ = Σᵢ<ⱼ xᵢ·xⱼ
    e₃ = Σᵢ<ⱼ<ₖ xᵢ·xⱼ·xₖ
    ⋮
    eₙ = x₁·x₂·...·xₙ

    Die eₖ bilden eine Basis des Invariantenrings unter Sₙ (Hauptsatz).
    Zusammenhang mit charakteristischem Polynom: det(tI - A) = tⁿ - e₁tⁿ⁻¹ + e₂tⁿ⁻² - ...

    @param n Anzahl der Variablen
    @return  Liste [e₁, e₂, ..., eₙ] als SymPy-Ausdrücke
    """
    # Variablen x₁,...,xₙ als SymPy-Symbole erzeugen
    syms = sp.symbols([f'x{i}' for i in range(1, n + 1)])
    result = []

    for k in range(1, n + 1):
        # eₖ = Summe aller k-elementigen Teilmengenprodukte
        ek = sp.Integer(0)
        for subset in itertools.combinations(syms, k):
            # Produkt der Elemente der Teilmenge
            product = sp.Integer(1)
            for var in subset:
                product *= var
            ek += product
        result.append(sp.expand(ek))

    return result


# itertools wird unten importiert
import itertools


def power_sum_symmetric_polynomials(n: int, max_degree: int) -> list:
    """
    Berechnet die Potenzsummen pₖ = x₁ᵏ + x₂ᵏ + ... + xₙᵏ für k=1,...,max_degree.

    Die Potenzsummen sind ebenfalls symmetrische Polynome und durch die
    Newton-Identitäten mit den elementarsymmetrischen Polynomen verbunden.

    @param n          Anzahl der Variablen
    @param max_degree Maximaler Grad (berechnet p₁, ..., p_{max_degree})
    @return           Liste [p₁, p₂, ..., p_{max_degree}] als SymPy-Ausdrücke
    """
    syms = sp.symbols([f'x{i}' for i in range(1, n + 1)])
    result = []

    for k in range(1, max_degree + 1):
        # pₖ = Summe der k-ten Potenzen aller Variablen
        pk = sum(x**k for x in syms)
        result.append(sp.expand(pk))

    return result


def newton_identities(n: int) -> dict:
    """
    Berechnet die Newton-Identitäten: Rekursion zwischen eₖ und pₖ.

    Newton-Identitäten verbinden Potenzsummen pₖ mit elementarsymmetrischen
    Polynomen eₖ:

    Für k ≤ n:   pₖ - e₁·pₖ₋₁ + e₂·pₖ₋₂ - ... + (-1)^{k-1}·eₖ₋₁·p₁ + (-1)^k·k·eₖ = 0
    Für k > n:   pₖ - e₁·pₖ₋₁ + e₂·pₖ₋₂ - ... + (-1)^n·eₙ·pₖ₋ₙ = 0

    @param n Anzahl der Variablen
    @return  Dict mit 'e' (elementarsymm.), 'p' (Potenzsummen), 'identities' (Gleichungen)
    """
    e_polys = elementary_symmetric_polynomials(n)
    p_polys = power_sum_symmetric_polynomials(n, n)

    identities = []
    syms = sp.symbols([f'x{i}' for i in range(1, n + 1)])

    for k in range(1, n + 1):
        # Newton-Identität für k: pₖ - e₁·pₖ₋₁ + ... + (-1)^k·k·eₖ = 0
        lhs = p_polys[k - 1]  # pₖ
        for i in range(1, k):
            # eᵢ·pₖ₋ᵢ mit Vorzeichen (-1)^i
            lhs = lhs + ((-1)**i) * e_polys[i - 1] * p_polys[k - i - 1]
        lhs = lhs + ((-1)**k) * k * e_polys[k - 1]
        identities.append({
            'k': k,
            'identity': f'p_{k} - e_1·p_{k-1} + ... + (-1)^{k}·{k}·e_{k} = 0',
            'verification': sp.expand(lhs),  # Sollte 0 sein
        })

    return {
        'n': n,
        'elementary_symmetric': e_polys,
        'power_sums': p_polys,
        'identities': identities,
        'note': 'verification=0 bestätigt die Newton-Identität',
    }


def polynomial_invariants_sn(n: int, degree: int) -> list:
    """
    Findet invariante Monome unter Sₙ (Permutationen von x₁,...,xₙ) bis Grad 'degree'.

    Ein Polynom f ist Sₙ-invariant, wenn für alle σ∈Sₙ gilt:
    f(x_{σ(1)}, ..., x_{σ(n)}) = f(x₁, ..., xₙ)

    Die elementarsymmetrischen Polynome e₁,...,eₙ erzeugen den gesamten Invariantenring.

    @param n      Anzahl der Variablen
    @param degree Maximaler Grad der gesuchten Invarianten
    @return       Liste von SymPy-Ausdrücken (invariante Polynome)
    """
    syms = sp.symbols([f'x{i}' for i in range(1, n + 1)])
    all_perms = list(permutations(range(n)))

    def sn_action(poly: sp.Expr, perm: tuple) -> sp.Expr:
        """
        Permutationsoperation auf Polynomen — simultane Substitution.
        simultaneous=True verhindert Ketteneffekte beim Variablentausch.
        """
        sub = [(syms[i], syms[perm[i]]) for i in range(n)]
        return sp.expand(poly.subs(sub, simultaneous=True))

    # Generiere alle Monome bis Grad 'degree' und bilde Orbitssummen
    invariants = []
    seen_orbits = set()

    # Alle Monome bis Grad 'degree' aufzählen
    for total_deg in range(0, degree + 1):
        for exponents in itertools.product(range(total_deg + 1), repeat=n):
            if sum(exponents) != total_deg:
                continue
            # Monom x₁^e₁ · x₂^e₂ · ... · xₙ^eₙ
            monomial = sp.Integer(1)
            for i, e in enumerate(exponents):
                monomial *= syms[i]**e

            # Kanonische Form des Exponentenvektors (Orbit-ID = sortierte Exponenten)
            orbit_id = tuple(sorted(exponents, reverse=True))
            if orbit_id in seen_orbits:
                continue
            seen_orbits.add(orbit_id)

            # Orbit-Summe = Reynolds-Operator angewendet auf Monom (ohne Division)
            orbit_sum = sp.Integer(0)
            orbit_set = set()
            for perm in all_perms:
                perm_monomial = sn_action(monomial, perm)
                orbit_set.add(perm_monomial)
            for m in orbit_set:
                orbit_sum += m

            orbit_sum_expanded = sp.expand(orbit_sum)
            if orbit_sum_expanded != 0:
                invariants.append(orbit_sum_expanded)

    return invariants


def molien_series(group_elements: list, matrices: list, max_degree: int = 10) -> list:
    """
    Berechnet die Molien-Reihe (Poincaré-Reihe) eines Invariantenrings.

    Die Molien-Reihe ist definiert als:
    M(t) = (1/|G|) · Σ_{g∈G} 1/det(I - t·g)

    Der Koeffizient von tᵏ gibt die Dimension des Raums der G-invarianten
    Polynome vom Grad k an.

    @param group_elements Liste der Gruppenelemente (nur für Länge verwendet)
    @param matrices       Liste der Matrizen der Gruppenelemente (als numpy-Arrays oder Listen)
    @param max_degree     Maximaler Grad für Koeffizientenberechnung
    @return               Liste [c₀, c₁, ..., c_{max_degree}] mit cₖ = dim(Inv_k)
    """
    import numpy as np

    G = len(matrices)
    if G == 0:
        raise ValueError("Leere Matrizenliste")

    # Formale Potenzreihen-Koeffizienten durch Taylor-Entwicklung berechnen
    coefficients = [0.0] * (max_degree + 1)
    t_sym = sp.Symbol('t')

    for mat in matrices:
        mat_sp = sp.Matrix(mat)
        n = mat_sp.shape[0]
        # I - t·g als symbolische Matrix
        identity = sp.eye(n)
        I_minus_tg = identity - t_sym * mat_sp
        # Determinante
        det_val = I_minus_tg.det()
        # 1/det als Potenzreihe bis max_degree
        series = sp.series(1 / det_val, t_sym, 0, max_degree + 1)
        # Koeffizienten extrahieren
        for k in range(max_degree + 1):
            coeff = series.coeff(t_sym, k)
            coefficients[k] += float(coeff)

    # Durch Gruppenordnung teilen
    coefficients = [c / G for c in coefficients]
    # Runden auf ganze Zahlen (theoretisch ganzzahlig)
    coefficients = [int(round(c)) for c in coefficients]

    return coefficients


def hilbert_basis_theorem_demo() -> dict:
    """
    Demonstriert den Hilbert-Basissatz für S₂-Invariantenring ℝ[x,y]^{S₂}.

    Hilbert-Basissatz: Jeder Invariantenring k[x₁,...,xₙ]^G ist endlich erzeugt
    (als k-Algebra).

    Für S₂ = {id, (12)} wirkend auf ℝ[x,y]:
    - g·f(x,y) = f(y,x) (Vertauschung)
    - Invariantenring ℝ[x,y]^{S₂} wird erzeugt von:
      e₁ = x+y  (elementarsymmetrisch, Grad 1)
      e₂ = x·y  (elementarsymmetrisch, Grad 2)

    @return Dict mit Erzeugern, Beispielen und Beweisskizze
    """
    x, y = sp.symbols('x y')

    # Erzeuger des Invariantenrings
    e1 = x + y       # elementarsymmetrisch Grad 1
    e2 = x * y       # elementarsymmetrisch Grad 2

    # Beispiel: x² + y² ist invariant und durch e₁, e₂ darstellbar
    # x² + y² = (x+y)² - 2xy = e₁² - 2e₂
    poly1 = x**2 + y**2
    repr1 = e1**2 - 2*e2
    check1 = sp.expand(poly1 - repr1) == 0

    # x³ + y³ = (x+y)³ - 3xy(x+y) = e₁³ - 3e₂e₁
    poly2 = x**3 + y**3
    repr2 = e1**3 - 3*e2*e1
    check2 = sp.expand(poly2 - repr2) == 0

    # Diskriminante (x-y)² = e₁² - 4e₂
    discriminant = (x - y)**2
    repr_disc = e1**2 - 4*e2
    check_disc = sp.expand(discriminant - repr_disc) == 0

    return {
        'ring': 'ℝ[x,y]^{S₂}',
        'group': 'S₂ = {id, (12)}: (x,y) ↦ (y,x)',
        'generators': {'e1': str(e1), 'e2': str(e2)},
        'theorem': (
            'Hilbert-Basissatz: ℝ[x,y]^{S₂} = ℝ[e₁, e₂] '
            '(endlich erzeugt von e₁=x+y, e₂=xy)'
        ),
        'examples': [
            {
                'poly': str(poly1),
                'representation': 'e₁² - 2e₂',
                'verified': check1,
            },
            {
                'poly': str(poly2),
                'representation': 'e₁³ - 3e₂e₁',
                'verified': check2,
            },
            {
                'poly': str(discriminant),
                'representation': 'e₁² - 4e₂',
                'verified': check_disc,
            },
        ],
        'note': (
            'Jedes S₂-invariante Polynom lässt sich als Polynom in e₁ und e₂ schreiben '
            '(Hauptsatz der Theorie symmetrischer Polynome)'
        ),
    }


def discriminant_invariant(coefficients: list) -> sp.Expr:
    """
    Berechnet die Diskriminante eines Polynoms als Invariante.

    Die Diskriminante ist definiert als:
    Δ(f) = Πᵢ<ⱼ (xᵢ - xⱼ)²

    wobei x₁,...,xₙ die Wurzeln von f sind. Δ ist invariant unter Sₙ
    (da das Quadrat das Vorzeichen absorbiert) und Polynomausdruck
    in den Koeffizienten von f.

    Für quadratisches f = ax² + bx + c: Δ = b² - 4ac

    @param coefficients Koeffizientenliste [aₙ, aₙ₋₁, ..., a₁, a₀] (höchste Potenz zuerst)
    @return             Diskriminante als SymPy-Ausdruck
    """
    n = len(coefficients) - 1  # Grad des Polynoms
    if n < 1:
        raise ValueError("Mindestens lineares Polynom erforderlich")

    # SymPy-Symbole für formale Wurzeln x₁,...,xₙ
    roots = sp.symbols([f'x{i}' for i in range(1, n + 1)])

    # Δ = Π_{i<j} (xᵢ - xⱼ)²
    disc = sp.Integer(1)
    for i in range(n):
        for j in range(i + 1, n):
            disc *= (roots[i] - roots[j])**2

    disc_expanded = sp.expand(disc)

    # Für niedrige Grade: Diskriminante als Funktion der Koeffizienten
    if n == 2:
        # Quadratisch: ax² + bx + c → Δ = b² - 4ac
        a, b, c = sp.symbols('a b c')
        if len(coefficients) == 3:
            a_val, b_val, c_val = coefficients
            return sp.sympify(b_val**2 - 4*a_val*c_val)
        return b**2 - 4*a*c
    elif n == 3:
        # Kubisch: ax³ + bx² + cx + d → Δ = 18abcd - 4b³d + b²c² - 4ac³ - 27a²d²
        if len(coefficients) == 4:
            a_val, b_val, c_val, d_val = coefficients
            delta = (18*a_val*b_val*c_val*d_val
                     - 4*b_val**3*d_val
                     + b_val**2*c_val**2
                     - 4*a_val*c_val**3
                     - 27*a_val**2*d_val**2)
            return sp.sympify(delta)
        a, b, c, d = sp.symbols('a b c d')
        return 18*a*b*c*d - 4*b**3*d + b**2*c**2 - 4*a*c**3 - 27*a**2*d**2

    # Allgemein: gib das formale Produkt zurück
    return disc_expanded


def fundamental_theorem_symmetric_poly(poly_str: str, n: int) -> dict:
    """
    Demonstriert den Hauptsatz der Theorie symmetrischer Polynome.

    Hauptsatz: Jedes symmetrische Polynom in x₁,...,xₙ ist eindeutig als
    Polynom in den elementarsymmetrischen Polynomen e₁,...,eₙ darstellbar.

    Algorithmus: Sukzessive Elimination des führenden Monoms durch Subtraktion
    geeigneter Produkte von eₖ.

    @param poly_str Symmetrisches Polynom als String (z.B. 'x1**2 + x2**2')
    @param n        Anzahl der Variablen
    @return         Dict mit Original, Darstellung in eₖ, Verifikation
    """
    syms = sp.symbols([f'x{i}' for i in range(1, n + 1)])
    e_polys = elementary_symmetric_polynomials(n)
    e_syms = sp.symbols([f'e{i}' for i in range(1, n + 1)])

    poly = sp.expand(sp.sympify(poly_str, locals={f'x{i}': syms[i-1] for i in range(1, n+1)}))

    # Prüfe Symmetrie (simultane Substitution vermeidet Ketteneffekte)
    all_perms = list(permutations(range(n)))
    is_sym = True
    for perm in all_perms:
        # Simultane Substitutionsliste: alle Variablen werden gleichzeitig getauscht
        sub = [(syms[i], syms[perm[i]]) for i in range(n)]
        if sp.expand(poly.subs(sub, simultaneous=True) - poly) != 0:
            is_sym = False
            break

    if not is_sym:
        return {
            'input': poly_str,
            'is_symmetric': False,
            'error': 'Polynom ist nicht symmetrisch — Hauptsatz nicht anwendbar',
        }

    # Darstellung durch elementarsymmetrische Polynome (für einfache Fälle)
    # Bekannte Identitäten direkt anwenden
    representation = "Keine direkte Darstellung berechnet"
    verified = False

    # Spezifische Fälle implementieren
    if n == 2:
        x1, x2 = syms
        e1, e2 = e_polys
        # Teste ob poly = c₀ + c₁·e₁ + c₂·e₁² + c₃·e₂ + c₄·e₁·e₂ + c₅·e₂² + ...
        # Durch Koeffizientenvergleich mit Grad ≤ 4
        a, b, c, d, e, f = sp.symbols('a b c d e f')
        ansatz = a + b*e1 + c*e1**2 + d*e2 + e*e1*e2 + f*e2**2
        diff = sp.expand(poly - ansatz)
        # Löse nach Koeffizienten
        eq_system = []
        for monomial in sp.Poly(diff, [x1, x2]).as_dict():
            eq_system.append(diff.coeff(x1, monomial[0]).coeff(x2, monomial[1]))

        try:
            solution = sp.solve(eq_system, [a, b, c, d, e, f])
            if solution:
                rep_poly = ansatz.subs(solution)
                rep_poly_expanded = sp.expand(rep_poly)
                verified = sp.expand(rep_poly_expanded - poly) == 0
                representation = str(rep_poly_expanded)
        except Exception:
            pass

    return {
        'input': poly_str,
        'n': n,
        'is_symmetric': is_sym,
        'elementary_symmetric': [str(e) for e in e_polys],
        'representation_in_e': representation,
        'verified': verified,
        'theorem': (
            'Hauptsatz: Jedes symmetrische Polynom in x₁,...,xₙ ist eindeutig '
            'als Polynom in e₁,...,eₙ darstellbar.'
        ),
    }
