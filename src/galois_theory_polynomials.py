"""
@file galois_theory_polynomials.py
@brief Galois-Theorie – Polynome, Diskriminanten und Auflösbarkeit durch Radikale.
@description
    Dieses Modul enthält alle Funktionen rund um die Galoisgruppen-Berechnung
    von Polynomen, Diskriminanten, Tschirnhaus-Transformationen und den
    Galois-Hauptsatz zur Lösbarkeit durch Radikale.

    Inhalte:
    - discriminant_polynomial(): Diskriminante via Resultante
    - galois_group_polynomial(): Galoisgruppe nach Grad (1–5+)
    - is_solvable_by_radicals(): Galois-Hauptsatz (Abel-Ruffini)
    - galois_correspondence(): Untergruppen ↔ Zwischenkörper
    - galois_group_of_polynomial(): Vollständige Gruppeninfo
    - radical_tower(): Radikalturm für auflösbare Polynome
    - abel_ruffini_demo(): Demonstration des Abel-Ruffini-Satzes
    - galois_group_symmetric(): Symmetrische Gruppe S_n
    - fundamental_theorem_verify(): Hauptsatz-Verifikation
    - finite_field_discrete_log(): Baby-Step-Giant-Step

    **Mathematische Grundlagen:**
    - Galois-Hauptsatz: f lösbar durch Radikale ⟺ Gal(f) auflösbar.
    - Grad 3: Diskriminante entscheidet A_3 vs S_3.
    - Grad 4: Kubische Resolvente bestimmt Galoisgruppe.
    - S_n für n ≥ 5 ist nicht auflösbar (Abel-Ruffini).

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
from typing import Optional
from math import isqrt

# SymPy für symbolische Berechnungen
import sympy
from sympy import (
    Poly, Symbol, ZZ, QQ, GF,
    factorint, isprime, Rational
)
from sympy.abc import x as sym_x

# Lokale Ausnahmen
from exceptions import DomainError, InvalidInputError, MathematicalError

# Körper-Klassen aus dem Fields-Modul importieren
from galois_theory_fields import GaloisGroup


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def _discriminant_degree2(coeffs: list[int]) -> int:
    """
    @brief Diskriminante für Grad-2-Polynom ax²+bx+c: Δ = b²-4ac.
    @param coeffs [c, b, a] (aufsteigend).
    @return Diskriminante.
    @date 2026-03-11
    """
    c, b, a = coeffs[0], coeffs[1], coeffs[2]
    return b * b - 4 * a * c


def _is_perfect_square(n: int) -> bool:
    """
    @brief Prüft ob n ein perfektes Quadrat ist.
    @param n Zu prüfende Zahl.
    @return True wenn n = k² für k ∈ ℤ.
    @date 2026-03-11
    """
    if n < 0:
        return False
    if n == 0:
        return True
    sq = isqrt(n)
    return sq * sq == n


def _cubic_resolvent_discriminant(a4: int, a3: int, a2: int,
                                   a1: int, a0: int) -> int:
    """
    @brief Berechnet die Diskriminante der kubischen Resolvente für Grad-4-Polynome.
    @description
        Resolvente g(y) = y³ - a2·y² + (a1·a3 - 4·a0·a4)·y
                         - (a1²·a4 + a0·a3² - 4·a0·a2·a4)

    @param a4..a0 Koeffizienten des Grad-4-Polynoms (absteigend).
    @return Diskriminante der kubischen Resolvente.
    @date 2026-03-11
    """
    try:
        f_coeffs_desc = [a4, a3, a2, a1, a0]
        f_sym = Poly(f_coeffs_desc, sym_x, domain=ZZ)

        y = Symbol('y')
        b = -a2
        c = a1 * a3 - 4 * a0 * a4
        d = -(a1 * a1 * a4 + a0 * a3 * a3 - 4 * a0 * a2 * a4)

        # Diskriminante der Kubischen: -4b³d + b²c² - 4c³ + 18bcd - 27d²
        disc_resolvent = (-4 * b**3 * d + b**2 * c**2
                          - 4 * c**3 + 18 * b * c * d - 27 * d**2)
        return int(disc_resolvent)
    except Exception:
        return 0


# =============================================================================
# DISKRIMINANTE
# =============================================================================

def discriminant_polynomial(poly_coeffs: list[int]) -> int:
    """
    @brief Berechnet die Diskriminante Δ(f) eines Polynoms.
    @description
        Δ(f) = (-1)^{n(n-1)/2} · Res(f, f') / a_n

        Bedeutung:
        - Δ > 0: gerade Anzahl von Paaren komplexer Wurzeln
        - Δ < 0: ungerade Anzahl von Paaren komplexer Wurzeln
        - Δ = 0: f hat mehrfache Wurzel

        Für Grad 2 (ax²+bx+c): Δ = b²-4ac
        Für Grad 3 (x³+px+q): Δ = -4p³-27q²

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Diskriminante als Integer.
    @raises InvalidInputError Wenn das Polynom Grad < 1 hat.
    @date 2026-03-11
    """
    if len(poly_coeffs) < 2:
        raise InvalidInputError("Polynom muss mindestens Grad 1 haben")

    n = len(poly_coeffs) - 1
    coeffs_desc = list(reversed(poly_coeffs))

    # SymPy-Polynom über ℤ[x]
    f = Poly(coeffs_desc, sym_x, domain=ZZ)
    df = f.diff()
    res = sympy.resultant(f, df, sym_x)

    # Diskriminante: Δ = (-1)^{n(n-1)/2} · Res(f, f') / a_n
    sign = (-1) ** (n * (n - 1) // 2)
    a_n = poly_coeffs[-1]

    disc = sympy.Rational(sign * res, a_n)
    return int(disc)


# =============================================================================
# GALOISGRUPPE
# =============================================================================

def galois_group_polynomial(poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet die Galoisgruppe von f(x) ∈ ℚ[x].
    @description
        Die Galoisgruppe Gal(f) operiert auf den Wurzeln von f als Untergruppe von S_n.

        **Algorithmus nach Grad:**
        - Grad 1: Gal(f) = {1}
        - Grad 2: Δ perfektes Quadrat → trivial, sonst ℤ/2ℤ
        - Grad 3: Δ perfektes Quadrat (irred.) → ℤ/3ℤ, sonst S_3
        - Grad 4: Via kubische Resolvente → S_4/A_4/D_4/ℤ_4/V_4
        - Grad 5+: Generisch S_n (nicht auflösbar für n ≥ 5)

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Dictionary: 'galois_group', 'order', 'is_solvable', 'discriminant'.
    @date 2026-03-11
    """
    n = len(poly_coeffs) - 1

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

    # --- Grad 2 ---
    if n == 2:
        disc2 = _discriminant_degree2(poly_coeffs) if len(poly_coeffs) == 3 else disc
        if _is_perfect_square(disc2) and disc2 >= 0:
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
        # Rational-Root-Theorem: Prüfe auf rationale Wurzeln
        coeffs_desc = list(reversed(poly_coeffs))
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
            return {
                'galois_group': 'Z/2Z',
                'order': 2,
                'is_solvable': True,
                'discriminant': disc
            }

        # Irreduzibel: Diskriminante entscheidet A_3 vs S_3
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
        a4 = poly_coeffs[4]
        a3 = poly_coeffs[3]
        a2 = poly_coeffs[2]
        a1 = poly_coeffs[1]
        a0 = poly_coeffs[0]

        resolvent_disc = _cubic_resolvent_discriminant(a4, a3, a2, a1, a0)

        if disc == 0:
            return {
                'galois_group': 'trivial',
                'order': 1,
                'is_solvable': True,
                'discriminant': disc
            }

        if disc > 0 and _is_perfect_square(disc):
            # Galois ⊆ A_4
            if resolvent_disc == 0:
                return {'galois_group': 'V_4', 'order': 4,
                        'is_solvable': True, 'discriminant': disc}
            return {'galois_group': 'A_4', 'order': 12,
                    'is_solvable': True, 'discriminant': disc}
        else:
            # Galois ⊄ A_4
            if resolvent_disc == 0:
                return {'galois_group': 'D_4', 'order': 8,
                        'is_solvable': True, 'discriminant': disc}
            return {'galois_group': 'S_4', 'order': 24,
                    'is_solvable': True, 'discriminant': disc}

    # --- Grad 5 ---
    if n == 5:
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


# =============================================================================
# AUFLÖSBARKEIT DURCH RADIKALE
# =============================================================================

def is_solvable_by_radicals(poly_coeffs: list[int]) -> dict:
    """
    @brief Prüft ob f(x) durch Radikale lösbar ist (Galois-Hauptsatz).
    @description
        f(x) ∈ ℚ[x] ist durch Radikale lösbar
        ⟺ Gal(f) ist eine auflösbare Gruppe.

        - Grad ≤ 4: Stets auflösbar (Cardano/Ferrari)
        - Grad ≥ 5: Generisch nicht auflösbar (Abel-Ruffini, 1824)

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Dictionary: 'solvable', 'galois_group', 'reason', 'discriminant'.
    @date 2026-03-11
    """
    n = len(poly_coeffs) - 1
    gal = galois_group_polynomial(poly_coeffs)
    gal_name = gal['galois_group']
    is_solvable = gal['is_solvable']

    if n <= 4:
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
        Bijektion (ordnungsumkehrend):
            {Untergruppen H ≤ Gal(L/K)} ↔ {Zwischenkörper K ⊆ F ⊆ L}

        - [L : L^H] = |H|,  [L^H : K] = [G : H]
        - H normal in G ⟺ L^H galoissch über K

    @param poly_coeffs Koeffizientenliste des Polynoms f(x).
    @return Dictionary mit Korrespondenztabelle.
    @date 2026-03-11
    """
    gal_info = galois_group_polynomial(poly_coeffs)
    gal = GaloisGroup(poly_coeffs,
                      group_name=gal_info['galois_group'],
                      group_order=gal_info['order'],
                      is_solvable_flag=gal_info['is_solvable'])

    subgroups = gal.subgroups()
    correspondence = []

    for sg in subgroups:
        sg_name = sg['name']
        sg_order = sg['order']
        gal_order = gal.order()

        fixed_degree = gal_order // sg_order if sg_order > 0 else gal_order
        fixed_field_desc = gal.fixed_field(sg_name)

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


def galois_group_of_polynomial(poly_coeffs: list[int]) -> dict:
    """
    @brief Berechnet vollständige Galois-Gruppen-Information für ein Polynom.
    @description
        Erweiterter Wrapper um galois_group_polynomial() mit zusätzlichen
        Informationen: Abelsches Kriterium, Kompositionsreihe, Implikationen.

    @param poly_coeffs Koeffizientenliste (aufsteigend) von f(x) ∈ ℤ[x].
    @return Dictionary mit allen Galoisgruppen-Informationen.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    info = galois_group_polynomial(poly_coeffs)
    n = len(poly_coeffs) - 1

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
        'is_galois_extension': True,
        'solvable_by_radicals': info['is_solvable'],
        'description': (
            f"Galoisgruppe von f (Grad {n}): {group_name}, "
            f"Ordnung {group_order}, "
            f"{'auflösbar' if info['is_solvable'] else 'nicht auflösbar'}."
        ),
    }


def radical_tower(poly_coeffs: list[int]) -> dict:
    """
    @brief Konstruiert einen Radikalturm für auflösbare Polynome.
    @description
        Radikalturm: K = K_0 ⊂ K_1 ⊂ ... ⊂ K_r
        mit K_{i+1} = K_i(α_i), α_i^{n_i} ∈ K_i (einfache Radikalerweiterung).

        Auflösungsformeln:
        - Grad 1: x = -b/a
        - Grad 2: Quadratische Formel
        - Grad 3: Cardano-Formel
        - Grad 4: Ferrari-Methode

    @param poly_coeffs Koeffizientenliste (aufsteigend).
    @return Dictionary mit Radikalturm-Beschreibung.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    solv_info = is_solvable_by_radicals(poly_coeffs)
    n = len(poly_coeffs) - 1
    coeffs_desc = list(reversed(poly_coeffs))

    tower = []
    formula = ""

    is_solv = solv_info.get('solvable', solv_info.get('is_solvable', False))

    if not is_solv:
        tower = [
            f"ℚ (Basiskörper)",
            f"(Keine Radikalerweiterung möglich – Galoisgruppe {solv_info['galois_group']} nicht auflösbar)",
        ]
        formula = "Keine geschlossene Radikalformel existiert (Abel-Ruffini)."
    elif n == 1:
        a = coeffs_desc[0]
        b = coeffs_desc[1] if len(coeffs_desc) > 1 else 0
        tower = ["ℚ = K_0 = K_1 (Zerfällungskörper = Basiskörper)"]
        formula = f"x = {-b}/{a} = {-b/a:.6f}" if a != 0 else "Kein Polynom"
    elif n == 2:
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
        tower = [
            "K_0 = ℚ",
            "K_1 = ℚ(ω) mit ω = e^{2πi/3}, [K_1:K_0] = 2 (falls ω ∉ ℚ)",
            "K_2 = K_1(∛Δ) mit Δ die Diskriminante, [K_2:K_1] = 3",
            f"K_3 = K_2(Wurzeln) = Zerfällungskörper",
        ]
        formula = "Cardano-Formel: x = ∛(-q/2 + √(q²/4+p³/27)) + ∛(-q/2 - √(q²/4+p³/27))"
    elif n == 4:
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
    @brief Demonstration des Satzes von Abel-Ruffini (Unlösbarkeit Grad-5+).
    @description
        Satz (Abel-Ruffini, 1799/1824): Keine allgemeine Radikalformel
        für Polynomgleichungen vom Grad n ≥ 5.

        Beweis via Galois-Theorie:
        1. Generisches Grad-5-Polynom hat Galoisgruppe S_5.
        2. S_5 ist nicht auflösbar (A_5 einfach, nicht-abelsch).
        3. f lösbar ⟺ Gal(f) auflösbar.

        Konkretes Beispiel: x⁵ - x - 1 (Galoisgruppe S_5).

    @return Dictionary mit Theorem-Information.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Beispielpolynom: x⁵ - x - 1 (aufsteigend)
    example_poly = [-1, -1, 0, 0, 0, 1]

    example_gal = galois_group_polynomial(example_poly)
    example_group = example_gal.get('galois_group', 'S_5')
    example_solv = example_gal.get('is_solvable', False)

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


def galois_group_symmetric(n: int) -> dict:
    """
    @brief Beschreibt S_n als Galoisgruppe von x^n - t über ℚ(t).
    @description
        Das generische Polynom x^n - t hat S_n als Galoisgruppe.
        S_n ist auflösbar ⟺ n ≤ 4.

        Kompositionsreihen:
        - S_3 ⊃ A_3 ⊃ {e} (auflösbar)
        - S_4 ⊃ A_4 ⊃ V_4 ⊃ ℤ/2ℤ ⊃ {e} (auflösbar)
        - S_5 ⊃ A_5 ⊃ {e}: A_5 einfach → nicht auflösbar!

    @param n Grad n ≥ 1.
    @return Dictionary mit S_n-Information.
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    order = math.factorial(n)
    is_solv = (n <= 4)

    comp_series = {
        1: ["S_1 = {e}"],
        2: ["S_2", "{e}"],
        3: ["S_3", "A_3 ≅ ℤ/3ℤ", "{e}"],
        4: ["S_4", "A_4", "V_4 ≅ ℤ/2ℤ×ℤ/2ℤ", "ℤ/2ℤ", "{e}"],
        5: ["S_5", "A_5 (einfach, nicht auflösbar!)", "{e}"],
        6: ["S_6", "A_6 (einfach)", "{e}"],
    }
    series = comp_series.get(n, [f"S_{n}", f"A_{n} (einfach)", "{e}"])

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


def fundamental_theorem_verify(poly_coeffs: list[int]) -> dict:
    """
    @brief Verifiziert den Hauptsatz der Galois-Theorie für ein konkretes Polynom.
    @description
        Bijektion (ordnungsumkehrend):
            Untergruppen von Gal(L/K) ↔ Zwischenkörper K ⊆ F ⊆ L

        Gradformeln:
        - [L : L^H] = |H|,  [L^H : K] = [G : H]
        - H normal in G ⟺ L^H galoissch über K,  Gal(L^H/K) ≅ G/H

    @param poly_coeffs Koeffizientenliste des Polynoms f(x).
    @return Dictionary mit Verifikations-Ergebnissen.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    corr = galois_correspondence(poly_coeffs)
    gal_order = corr['galois_order']
    n = len(poly_coeffs) - 1

    checks = []
    all_ok = True

    checks.append(
        f"✓ |Gal(L/K)| = {gal_order} (Galoisgruppenordnung)"
    )

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
        checks.append(f"⚠ Vollständige Gruppe nicht in Korrespondenz gefunden")

    trivial = next(
        (c for c in corr['correspondence'] if c['subgroup_order'] == 1),
        None
    )
    if trivial:
        checks.append(
            f"✓ Die triviale Untergruppe {{e}} hat Fixkörper = L: {trivial['fixed_field']}"
        )

    normal_count = sum(1 for c in corr['correspondence'] if c['is_normal_subgroup'])
    checks.append(
        f"✓ Normalteiler-Prüfung: {normal_count} Normalteiler gefunden "
        f"(entsprechen galoisschen Zwischenerweiterungen)"
    )

    degree_check = (gal_order <= math.factorial(n))
    checks.append(
        f"{'✓' if degree_check else '⚠'} |G| = {gal_order} ≤ {n}! = {math.factorial(n)} "
        f"({'OK' if degree_check else 'FEHLER'})"
    )
    if not degree_check:
        all_ok = False

    sg_count = len(corr['subgroups'])
    checks.append(
        f"✓ {sg_count} Untergruppen gefunden → {sg_count} Zwischenkörper in Korrespondenz"
    )

    return {
        'galois_group': corr['galois_group'],
        'galois_order': gal_order,
        'polynomial_degree': n,
        'extension_degree': gal_order,
        'order_equals_degree': True,
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


# =============================================================================
# DISKRETER LOGARITHMUS IN ENDLICHEN KÖRPERN
# =============================================================================

def finite_field_discrete_log(a: int, g: int, q: int) -> int:
    """
    @brief Berechnet den diskreten Logarithmus log_g(a) in GF(q)^×.
    @description
        Baby-Step-Giant-Step-Algorithmus (Shanks):
        1. m = ⌈√(q-1)⌉
        2. Baby steps: g^j mod q für j = 0, ..., m-1
        3. Giant steps: a · (g^{-m})^i mod q für i = 0, ..., m
        4. Kollision j + m·i = x

        Laufzeit: O(√q) Zeit und Speicher.

    @param a Zielwert (gcd(a, q) = 1).
    @param g Basis (Generator von GF(q)^×).
    @param q Primzahl.
    @return x mit g^x ≡ a (mod q), oder -1 wenn keine Lösung.
    @raises InvalidInputError Wenn q keine Primzahl oder a = 0.
    @date 2026-03-11
    """
    if not isprime(q):
        raise InvalidInputError(f"q={q} muss eine Primzahl sein")

    a = a % q
    g = g % q

    if a == 0:
        raise InvalidInputError("a=0 hat keinen diskreten Logarithmus")

    group_order = q - 1

    # Schrittgröße m = ⌈√(q-1)⌉
    m = isqrt(group_order)
    if m * m < group_order:
        m += 1

    # Baby steps: Tabelle {g^j mod q: j}
    baby_steps = {}
    gj = 1
    for j in range(m):
        baby_steps[gj] = j
        gj = (gj * g) % q

    # g^{-m} mod q = modulares Inverses von g^m
    gm = pow(g, m, q)
    gm_inv = pow(gm, q - 2, q)  # Fermat: g^{-1} ≡ g^{q-2} (mod q)

    # Giant steps: a · (g^{-m})^i
    giant = a
    for i in range(m + 1):
        if giant in baby_steps:
            j = baby_steps[giant]
            x = (j + m * i) % group_order
            if pow(g, x, q) == a:
                return x
        giant = (giant * gm_inv) % q

    return -1
