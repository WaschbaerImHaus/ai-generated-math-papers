"""
@file galois_theory_constructibility.py
@brief Galois-Theorie – Konstruierbarkeit mit Zirkel und Lineal.
@description
    Dieses Modul implementiert das Gauss-Wantzel-Theorem und alle zugehörigen
    Konstruierbarkeits-Beweise für reguläre Polygone.

    **Gauss-Wantzel-Theorem (1796/1837):**
    Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann, wenn
    $$n = 2^k \\cdot p_1 \\cdot p_2 \\cdots p_r$$
    wobei k ≥ 0 und p_i verschiedene Fermat-Primzahlen sind.

    **Fermat-Primzahlen** sind Primzahlen der Form $F_m = 2^{2^m} + 1$:
    - F_0 = 3, F_1 = 5, F_2 = 17, F_3 = 257, F_4 = 65537.
    Weitere bekannte Fermat-Primzahlen existieren nicht (F_5 = 4294967297 = 641 × 6700417).

    **Galois-theoretische Begründung:**
    n-Eck konstruierbar ⟺ [ℚ(ζ_n):ℚ] = φ(n) ist eine 2-Potenz
    ⟺ Gal(ℚ(ζ_n)/ℚ) ≅ (ℤ/nℤ)^× hat 2-Potenz-Ordnung
    ⟺ Turm von Körpererweiterungen vom Grad 2 existiert.

    **Unmöglichkeitsbeweise:**
    - Würfelverdopplung: ∛2 nicht konstruierbar (φ(irred. Poly.) = 3, keine 2-Potenz)
    - Winkeldreiteilung: cos(20°) erfüllt 8x³-6x-1=0, irreduzibel über ℚ, Grad 3
    - Kreisquadratur: π transzendent (Lindemann 1882), nicht algebraisch

    Enthaltene Funktionen:
    - is_constructible(): Prüft Konstruierbarkeit nach Gauss-Wantzel
    - construct_regular_polygon(): Vollständige Analyse mit Galois-Gruppen-Info
    - gauss_wantzel_theorem(): Beschreibung des Satzes mit Beispielen
    - impossible_constructions(): Unmöglichkeitsbeweise für klassische Probleme
    - fermat_prime_check(): Prüft ob eine Zahl eine Fermat-Primzahl ist
    - constructible_numbers(): Listet konstruierbare n-Ecke bis zu einem Limit auf
    - degree_tower_decomposition(): Zerlegt den Körpererweiterungs-Turm

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath
from typing import Optional

# SymPy für symbolische Berechnungen (φ(n), Faktorisierung)
import sympy
from sympy import (
    totient, factorint, isprime, cyclotomic_poly
)

# Lokale Ausnahmen
from exceptions import InvalidInputError, MathematicalError


# =============================================================================
# BEKANNTE FERMAT-PRIMZAHLEN
# =============================================================================

# Alle bekannten Fermat-Primzahlen F_m = 2^(2^m) + 1 für m = 0,...,4
FERMAT_PRIMES = [3, 5, 17, 257, 65537]


def fermat_prime_check(p: int) -> dict:
    """
    @brief Prüft ob eine Zahl eine Fermat-Primzahl ist und berechnet den Index.
    @description
        Eine Fermat-Primzahl hat die Form $F_m = 2^{2^m} + 1$.
        Bekannte Fermat-Primzahlen: F_0=3, F_1=5, F_2=17, F_3=257, F_4=65537.

        Wegen Peymanns Computerberechnungen ist bekannt:
        - F_5 = 4294967297 = 641 × 6700417 (keine Primzahl)
        - F_6 bis F_32: zusammengesetzt
        - Status von F_33 und größeren: unbekannt

    @param p Zu prüfende Zahl (p ≥ 2).
    @return Dictionary mit:
            - 'p': int
            - 'is_prime': bool
            - 'is_fermat_prime': bool
            - 'fermat_index': int oder None (m, falls F_m = p)
            - 'form': str, z.B. "2^(2^2)+1 = 17"
    @raises InvalidInputError Wenn p < 2.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if p < 2:
        raise InvalidInputError(f"p={p} muss ≥ 2 sein")

    # Primtest
    is_p = bool(isprime(p))

    # Prüfe ob p = 2^(2^m) + 1 für ein m ≥ 0
    fermat_index = None
    form_str = ""

    if is_p:
        # Suche m mit 2^(2^m) + 1 = p
        # Nur für bekannte Fermat-Primzahlen prüfen (bis m=4 realistisch)
        for m in range(40):
            val = 2 ** (2 ** m) + 1
            if val == p:
                fermat_index = m
                form_str = f"2^(2^{m})+1 = {p}"
                break
            if val > p:
                # val ist bereits größer als p → kein Fermat-Primzahl
                break

    is_fermat = (fermat_index is not None)

    return {
        'p': p,
        'is_prime': is_p,
        'is_fermat_prime': is_fermat,
        'fermat_index': fermat_index,
        'form': form_str if is_fermat else f"{p} ist keine Fermat-Zahl",
        'known_fermat_primes': FERMAT_PRIMES,
        'description': (
            f"{p} ist {'eine' if is_fermat else 'keine'} Fermat-Primzahl."
            + (f" F_{fermat_index} = 2^(2^{fermat_index})+1" if is_fermat else "")
        ),
    }


def is_constructible(n: int) -> dict:
    """
    @brief Prüft ob ein reguläres n-Eck mit Zirkel und Lineal konstruierbar ist.
    @description
        Das Gauss-Wantzel-Theorem (1796/1837) charakterisiert konstruierbare reguläre
        Polygone vollständig:

        **Satz (Gauss-Wantzel):**
        Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann, wenn
        $$n = 2^k \\cdot p_1 \\cdot p_2 \\cdots p_r$$
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
            - 'factorization': dict
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Bekannte Fermat-Primzahlen (alle bekannten)
    fermat_primes_set = set(FERMAT_PRIMES)

    # Triviale Fälle: n=1 oder n=2
    if n == 1 or n == 2:
        return {
            'n': n,
            'is_constructible': True,
            'phi_n': int(totient(n)),
            'phi_is_power_of_2': True,
            'fermat_primes_used': [],
            'power_of_2_factor': n if n == 2 else 1,
            'reason': f"n={n} ist trivial konstruierbar.",
            'factorization': {n: 1} if n > 1 else {1: 1},
        }

    # Euler-Phi-Wert berechnen
    phi_n = int(totient(n))

    # Hilfsfunktion: ist k eine 2-Potenz?
    def is_power_of_2(k: int) -> bool:
        """Gibt True zurück, wenn k eine 2-Potenz ist (k > 0)."""
        return k > 0 and (k & (k - 1)) == 0

    phi_is_2_power = is_power_of_2(phi_n)

    # Faktorisierung von n berechnen
    facts = sympy.factorint(n)

    # 2er-Anteil bestimmen
    two_exp = facts.get(2, 0)
    power_of_2 = 2 ** two_exp

    # Ungerade Primfaktoren klassifizieren
    odd_prime_factors = {p: e for p, e in facts.items() if p != 2}

    used_fermat = []       # Fermat-Primzahlen, die vorkommen (Exponent 1)
    non_fermat_odd = []    # Ungerade Primfaktoren, die keine Fermat-Primzahlen sind
    repeated_odd = []      # Ungerade Primfaktoren mit Exponent > 1

    for p, e in odd_prime_factors.items():
        if e > 1:
            # Primfaktor mit Exponent > 1 → nicht konstruierbar
            repeated_odd.append(p)
        elif p in fermat_primes_set:
            used_fermat.append(p)
        else:
            non_fermat_odd.append(p)

    # Konstruierbar genau dann, wenn keine Nicht-Fermat-Faktoren und keine
    # wiederholten Faktoren
    constructible = (phi_is_2_power and not non_fermat_odd and not repeated_odd)

    # Begründung für das Ergebnis zusammenstellen
    if constructible:
        if not used_fermat:
            reason = (
                f"n={n} = 2^{two_exp} ist eine 2-Potenz. "
                f"Das {n}-Eck ist konstruierbar."
            )
        else:
            fp_str = " · ".join(str(p) for p in sorted(used_fermat))
            reason = (
                f"n={n} = 2^{two_exp} · {fp_str}. "
                f"Alle ungeraden Faktoren sind verschiedene Fermat-Primzahlen. "
                f"Das {n}-Eck ist konstruierbar (Gauss-Wantzel)."
            )
    elif non_fermat_odd:
        reason = (
            f"n={n} enthält den ungeraden Primfaktor {non_fermat_odd[0]}, "
            f"der keine Fermat-Primzahl ist. Das {n}-Eck ist NICHT konstruierbar."
        )
    elif repeated_odd:
        reason = (
            f"n={n} enthält einen ungeraden Primfaktor {repeated_odd[0]} mit "
            f"Exponent > 1. Das {n}-Eck ist NICHT konstruierbar."
        )
    else:
        reason = (
            f"n={n} ist nicht konstruierbar (φ(n)={phi_n} ist keine 2-Potenz)."
        )

    return {
        'n': n,
        'is_constructible': constructible,
        'phi_n': phi_n,
        'phi_is_power_of_2': phi_is_2_power,
        'fermat_primes_used': sorted(used_fermat),
        'power_of_2_factor': power_of_2,
        'reason': reason,
        'factorization': dict(facts),
        'non_fermat_primes': non_fermat_odd,
        'repeated_primes': repeated_odd,
    }


def construct_regular_polygon(n: int) -> dict:
    """
    @brief Analysiert die Konstruierbarkeit und Galois-Theorie eines regulären n-Ecks.
    @description
        Verbindet das Gauss-Wantzel-Theorem mit der Galois-Theorie:

        Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann,
        wenn der Turm von Körpererweiterungen über ℚ nur Schritte vom Grad 2 hat:

        $$\\mathbb{Q} = K_0 \\subset K_1 \\subset \\cdots \\subset K_r = \\mathbb{Q}(\\zeta_n)$$

        mit $[K_{i+1} : K_i] = 2$ für alle i.

        Die Ecken des regulären n-Ecks sind die n-ten Einheitswurzeln:
        $$\\zeta_n^k = e^{2\\pi i k/n}, \\quad k = 0, 1, \\ldots, n-1$$

        **Historische Bedeutung:**
        Gauss zeigte 1796 (mit 19 Jahren), dass das reguläre 17-Eck konstruierbar ist
        (F_2 = 17 ist eine Fermat-Primzahl). Dies war der erste neue Fortschritt
        seit der Antike.

    @param n Anzahl der Ecken (n ≥ 3).
    @return Dictionary mit:
            - 'n': int
            - 'is_constructible': bool
            - 'constructibility_info': dict von is_constructible(n)
            - 'vertices': list[complex], Koordinaten der Ecken
            - 'historical_note': str
            - 'phi_n': int
            - 'description': str
    @raises InvalidInputError Wenn n < 3.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 3:
        raise InvalidInputError(
            f"n={n} muss ≥ 3 sein (ein Polygon braucht mind. 3 Ecken)"
        )

    # Konstruierbarkeits-Analyse nach Gauss-Wantzel
    constr_info = is_constructible(n)

    # Koordinaten der Ecken des regulären n-Ecks in der komplexen Ebene
    # k-te Ecke: ζ_n^k = e^{2πik/n}
    vertices = [cmath.exp(2j * math.pi * k / n) for k in range(n)]

    # Historische Notizen zu besonderen Polygonen
    historical_notes = {
        3:     "Gleichseitiges Dreieck: seit der Antike bekannt.",
        4:     "Quadrat: seit der Antike bekannt.",
        5:     "Reguläres Pentagon: Euklid, ~300 v. Chr.",
        6:     "Reguläres Hexagon: seit der Antike bekannt.",
        8:     "Reguläres Oktogon: seit der Antike bekannt.",
        10:    "Reguläres Dekagon: Euklid, ~300 v. Chr.",
        15:    "Reguläres 15-Eck: Euklid, ~300 v. Chr.",
        17:    "Reguläres 17-Eck: C.F. Gauss, 1796 (revolutionäre Entdeckung!).",
        257:   "Reguläres 257-Eck: Magnus Georg Paucker, 1822.",
        65537: "Reguläres 65537-Eck: Johann Gustav Hermes, 1894 (10 Jahre Arbeit!).",
    }

    # Historische Notiz suchen oder generieren
    note = historical_notes.get(n, "")
    if not note:
        if constr_info['is_constructible']:
            note = (
                f"Das reguläre {n}-Eck ist konstruierbar "
                f"(Gauss-Wantzel-Kriterium erfüllt)."
            )
        else:
            note = (
                f"Das reguläre {n}-Eck ist NICHT konstruierbar mit Zirkel und Lineal."
            )

    # Euler-Phi-Wert (= Grad des Kreisteilungskörpers über ℚ)
    phi_n = constr_info['phi_n']

    return {
        'n': n,
        'is_constructible': constr_info['is_constructible'],
        'constructibility_info': constr_info,
        'vertices': vertices,
        'historical_note': note,
        'phi_n': phi_n,
        'description': (
            f"Reguläres {n}-Eck: "
            f"{'Konstruierbar' if constr_info['is_constructible'] else 'NICHT konstruierbar'} "
            f"mit Zirkel und Lineal. φ({n}) = {phi_n}."
        ),
        'galois_degree': phi_n,
        'field_extension': f"[ℚ(ζ_{n}):ℚ] = φ({n}) = {phi_n}",
    }


def gauss_wantzel_theorem() -> dict:
    """
    @brief Gibt eine vollständige Beschreibung des Gauss-Wantzel-Theorems zurück.
    @description
        Das Gauss-Wantzel-Theorem ist eines der schönsten Ergebnisse der
        klassischen und modernen Mathematik, da es antike Geometrie mit
        moderner Algebra (Galois-Theorie) verbindet.

        **Geschichte:**
        - 1796: Carl Friedrich Gauss zeigt, dass 17-Eck konstruierbar (mit 19 Jahren).
        - 1801: Gauss beweist die Notwendigkeit in seinen "Disquisitiones Arithmeticae".
        - 1837: Pierre Wantzel beweist die Unmöglichkeit für nicht-Fermat-Primzahlen
                und liefert den vollständigen Beweis.

        **Satz:**
        Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar ⟺
        n = 2^k · p_1 · p_2 · ... · p_r, wo p_i verschiedene Fermat-Primzahlen.

        **Beweisidee:**
        1. Werkzeuge: Zirkel+Lineal erlauben Addition, Subtraktion, Multiplikation,
           Division und Quadratwurzeln.
        2. Algebraisch: konstruierbar = in einem Turm von Grad-2-Erweiterungen über ℚ.
        3. Galois-Theorie: Konstruierbarkeit ⟺ φ(n) ist eine 2-Potenz.

    @return Dictionary mit vollständiger Beschreibung des Theorems.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Beispiele für konstruierbare und nicht-konstruierbare Polygone
    constructible_examples = [3, 4, 5, 6, 8, 10, 12, 15, 16, 17, 20, 24]
    non_constructible_examples = [7, 9, 11, 13, 14, 18, 19, 21, 22, 23, 25]

    # Alle bekannten Fermat-Primzahlen mit ihren Formeln
    fermat_prime_list = [
        {'m': 0, 'F_m': 3,     'formula': '2^(2^0)+1 = 2^1+1 = 3'},
        {'m': 1, 'F_m': 5,     'formula': '2^(2^1)+1 = 2^2+1 = 5'},
        {'m': 2, 'F_m': 17,    'formula': '2^(2^2)+1 = 2^4+1 = 17'},
        {'m': 3, 'F_m': 257,   'formula': '2^(2^3)+1 = 2^8+1 = 257'},
        {'m': 4, 'F_m': 65537, 'formula': '2^(2^4)+1 = 2^16+1 = 65537'},
    ]

    return {
        'theorem_name': "Gauss-Wantzel-Theorem",
        'year': "1796 (Gauss) / 1837 (vollständiger Beweis durch Wantzel)",
        'statement': (
            "Ein reguläres n-Eck ist mit Zirkel und Lineal konstruierbar genau dann, "
            "wenn n = 2^k · p_1 · p_2 · ... · p_r, wobei k ≥ 0 und p_1,...,p_r "
            "paarweise verschiedene Fermat-Primzahlen sind."
        ),
        'fermat_prime_condition': (
            "Eine ungerade Primzahl p ist genau dann als Faktor erlaubt, wenn "
            "p eine Fermat-Primzahl der Form 2^(2^m)+1 ist und nur mit Exponent 1 vorkommt."
        ),
        'galois_reformulation': (
            "φ(n) = [ℚ(ζ_n):ℚ] ist eine 2-Potenz ⟺ Gal(ℚ(ζ_n)/ℚ) hat 2-Potenz-Ordnung "
            "⟺ Turm von Grad-2-Erweiterungen ℚ = K_0 ⊂ K_1 ⊂ ... ⊂ K_r = ℚ(ζ_n) existiert."
        ),
        'known_fermat_primes': fermat_prime_list,
        'open_question': (
            "Gibt es unendlich viele Fermat-Primzahlen? "
            "Bisher nur F_0 bis F_4 bekannt (Stand 2026). "
            "F_5 = 4294967297 = 641 × 6700417 ist zusammengesetzt."
        ),
        'constructible_examples': constructible_examples,
        'non_constructible_examples': non_constructible_examples,
        'historical_note': (
            "Gauss' Entdeckung des konstruierbaren 17-Ecks (1796) war der erste Fortschritt "
            "in der Konstruierbarkeits-Theorie seit über 2000 Jahren (seit Euklid). "
            "Der vollständige Beweis des Nicht-Konstruierbarkeits-Teils gelang erst Wantzel (1837)."
        ),
    }


def impossible_constructions() -> dict:
    """
    @brief Beschreibt die drei klassischen unmöglichen Konstruktionen der Antike.
    @description
        Die drei klassischen geometrischen Probleme der Antike waren jahrhundertelang
        ungelöst. Mit Galois-Theorie und algebraischen Methoden konnten ihre
        Unmöglichkeiten bewiesen werden.

        **1. Würfelverdopplung (Delisches Problem):**
        Gegeben ein Würfel mit Seite 1, konstruiere einen Würfel mit doppeltem Volumen.
        → Erfordert ∛2, aber [ℚ(∛2):ℚ] = 3 (keine 2-Potenz). Unmöglich.

        **2. Winkeldreiteilung:**
        Teile einen beliebigen Winkel in drei gleiche Teile.
        → cos(20°) = cos(60°/3) erfüllt 8x³-6x-1=0, irreduzibel über ℚ, Grad 3.
        Unmöglich für allgemeine Winkel (z.B. 60°).

        **3. Kreisquadratur:**
        Konstruiere ein Quadrat mit gleichem Flächeninhalt wie ein gegebener Kreis.
        → Erfordert π, aber π ist transzendent (Lindemann, 1882).
        Transzendente Zahlen sind nicht algebraisch → nicht konstruierbar.

    @return Dictionary mit Beschreibungen aller drei unmöglichen Konstruktionen.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    return {
        'cube_duplication': {
            'name': "Würfelverdopplung (Delisches Problem)",
            'problem': "Konstruiere ∛2 mit Zirkel und Lineal.",
            'obstacle': "∛2 ist Wurzel von x³-2, welches irreduzibel über ℚ ist.",
            'degree': "[ℚ(∛2):ℚ] = 3",
            'reason': (
                "3 ist keine 2-Potenz. Der Erweiterungsgrad muss eine 2-Potenz sein "
                "für Konstruierbarkeit. Da [ℚ(∛2):ℚ] = 3, ist ∛2 nicht konstruierbar."
            ),
            'algebraic_proof': (
                "Angenommen ∛2 sei konstruierbar. Dann liegt ∛2 in einem Turm "
                "ℚ = K_0 ⊂ K_1 ⊂ ... ⊂ K_r mit je [K_{i+1}:K_i] = 2. "
                "Also 2^r = [K_r:ℚ] ≥ [ℚ(∛2):ℚ] = 3. Aber 3 teilt 2^r nicht. "
                "Widerspruch."
            ),
            'proved_by': "Pierre Wantzel, 1837",
            'impossible': True,
        },
        'angle_trisection': {
            'name': "Winkeldreiteilung",
            'problem': "Dreiteilung des 60°-Winkels (konstruiere 20°).",
            'obstacle': "cos(20°) ist Wurzel von 8x³-6x-1, irreduzibel über ℚ.",
            'degree': "[ℚ(cos(20°)):ℚ] = 3",
            'reason': (
                "Wenn man 60° dritteln könnte, könnte man cos(20°) konstruieren. "
                "Das Minimalpolynom von cos(20°) hat Grad 3 (keine 2-Potenz). "
                "Also ist die Dreiteilung des 60°-Winkels unmöglich."
            ),
            'algebraic_proof': (
                "cos(60°) = 1/2. Mittels Dreifachwinkel-Formel: "
                "cos(3θ) = 4cos³(θ) - 3cos(θ). "
                "Für θ=20°: 1/2 = 4cos³(20°) - 3cos(20°). "
                "Also: 8x³ - 6x - 1 = 0 (x = cos(20°)). "
                "Dieses Polynom ist irreduzibel über ℚ (Rational Root Test: "
                "mögliche Wurzeln ±1, ±1/2, ±1/4, ±1/8, keine ist Wurzel). "
                "Grad 3 ist keine 2-Potenz → nicht konstruierbar."
            ),
            'note': "Spezielle Winkel (z.B. 90°) können dennoch drittiert werden!",
            'proved_by': "Pierre Wantzel, 1837",
            'impossible': True,  # Für allgemeine Winkel
        },
        'circle_squaring': {
            'name': "Kreisquadratur",
            'problem': "Konstruiere ein Quadrat mit Fläche π (= Fläche des Einheitskreises).",
            'obstacle': "π ist transzendent.",
            'degree': "π ∉ algebraische Zahlen → kein endlicher Grad über ℚ",
            'reason': (
                "Konstruierbare Zahlen sind algebraisch über ℚ (Nullstellen von "
                "Polynomen mit rationalen Koeffizienten). Da π transzendent ist "
                "(d.h. kein Polynom über ℚ annulliert π), ist π nicht konstruierbar. "
                "Daher ist √π (für ein Quadrat der Seite √π) ebenfalls transzendent."
            ),
            'algebraic_proof': (
                "Lindemann (1882) bewies: π ist transzendent. "
                "Konstruierbare Zahlen ⊂ algebraische Zahlen. "
                "π ∉ algebraische Zahlen ⇒ π ∉ konstruierbare Zahlen. "
                "Damit ist Kreisquadratur unmöglich."
            ),
            'proved_by': "Ferdinand von Lindemann, 1882",
            'impossible': True,
        },
    }


def constructible_numbers(limit: int = 100) -> dict:
    """
    @brief Listet alle konstruierbaren n-Ecke für n von 3 bis limit auf.
    @description
        Berechnet für jedes n ∈ {3,...,limit}, ob das reguläre n-Eck konstruierbar ist.
        Grundlage ist das Gauss-Wantzel-Theorem.

    @param limit Obere Grenze der Suche (limit ≥ 3, Standard: 100).
    @return Dictionary mit:
            - 'constructible': list[int], alle konstruierbaren n
            - 'non_constructible': list[int], alle nicht-konstruierbaren n
            - 'count_constructible': int
            - 'count_non_constructible': int
            - 'details': dict[int, dict], Details für jedes n
    @raises InvalidInputError Wenn limit < 3.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if limit < 3:
        raise InvalidInputError(f"limit={limit} muss ≥ 3 sein")

    constructible_list = []
    non_constructible_list = []
    details = {}

    # Für jedes n von 3 bis limit prüfen
    for n in range(3, limit + 1):
        result = is_constructible(n)
        details[n] = result
        if result['is_constructible']:
            constructible_list.append(n)
        else:
            non_constructible_list.append(n)

    return {
        'limit': limit,
        'constructible': constructible_list,
        'non_constructible': non_constructible_list,
        'count_constructible': len(constructible_list),
        'count_non_constructible': len(non_constructible_list),
        'details': details,
        'description': (
            f"Von n=3 bis n={limit}: "
            f"{len(constructible_list)} konstruierbare, "
            f"{len(non_constructible_list)} nicht-konstruierbare reguläre Polygone."
        ),
    }


def degree_tower_decomposition(n: int) -> dict:
    """
    @brief Zerlegt den Körpererweiterungs-Turm ℚ ⊂ ... ⊂ ℚ(ζ_n) in Grad-2-Schritte.
    @description
        Wenn n-Eck konstruierbar ist, existiert ein Turm von Körpererweiterungen:
        $$\\mathbb{Q} = K_0 \\subset K_1 \\subset \\cdots \\subset K_r = \\mathbb{Q}(\\zeta_n)$$
        mit $[K_{i+1}:K_i] = 2$ für alle i, also $r = \\log_2(\\varphi(n))$.

        Dieser Turm entspricht einer auflösbaren Untergruppen-Kette in der Galoisgruppe
        $\\text{Gal}(\\mathbb{Q}(\\zeta_n)/\\mathbb{Q}) \\cong (\\mathbb{Z}/n\\mathbb{Z})^\\times$.

    @param n Positive ganze Zahl ≥ 1.
    @return Dictionary mit:
            - 'n': int
            - 'is_constructible': bool
            - 'phi_n': int
            - 'tower_height': int oder None
            - 'tower_description': str
            - 'galois_group_order': int (= phi_n)
    @raises InvalidInputError Wenn n < 1.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 1:
        raise InvalidInputError(f"n={n} muss ≥ 1 sein")

    # Konstruierbarkeits-Info holen
    constr = is_constructible(n)
    phi_n = constr['phi_n']

    # Turmhöhe berechnen (nur wenn konstruierbar)
    tower_height = None
    tower_desc = ""

    if constr['is_constructible'] and phi_n > 0:
        # Turmhöhe = log_2(φ(n)), da jeder Schritt Grad 2 hat
        tower_height = int(math.log2(phi_n)) if phi_n > 0 else 0
        tower_desc = (
            f"ℚ = K_0 ⊂ K_1 ⊂ ... ⊂ K_{tower_height} = ℚ(ζ_{n}), "
            f"je [K_{{i+1}}:K_i] = 2, Gesamtgrad φ({n}) = {phi_n} = 2^{tower_height}."
        )
    else:
        tower_desc = (
            f"ℚ(ζ_{n})/ℚ hat Grad φ({n}) = {phi_n}, "
            f"welcher keine 2-Potenz ist → Turm mit nur Grad-2-Schritten nicht möglich."
        )

    return {
        'n': n,
        'is_constructible': constr['is_constructible'],
        'phi_n': phi_n,
        'tower_height': tower_height,
        'tower_description': tower_desc,
        'galois_group_order': phi_n,
        'galois_group': f"Gal(ℚ(ζ_{n})/ℚ) ≅ (ℤ/{n}ℤ)^×, Ordnung φ({n}) = {phi_n}",
        'description': (
            f"Körpererweiterungs-Turm für n={n}: "
            f"{'konstruierbar' if constr['is_constructible'] else 'nicht konstruierbar'}."
        ),
    }
