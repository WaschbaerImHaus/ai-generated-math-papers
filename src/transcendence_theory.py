"""
transcendence_theory.py – Transzendenztheorie

Dieses Modul implementiert Algorithmen und Demonstrationen zur Transzendenztheorie:
- Klassifikation algebraischer und transzendenter Zahlen
- Liouville-Zahlen und Irrationali­tätsmaße
- Lindemann-Weierstrass-Satz (Demo)
- Gelfond-Schneider, Baker-Theorem (Demo)
- Algebraische Zahlen als Klasse

@author    Kurt Ingwer
@version   1.0.0
@timestamp 2026-03-10
"""

import math
import cmath
import fractions
import itertools
from typing import Optional

# SymPy für symbolische Berechnungen
try:
    import sympy
    from sympy import (
        Symbol, Rational, sqrt, Poly, minimal_polynomial as sympy_minimal_poly,
        cyclotomic_poly, factorint, isprime, cos, pi, exp, log, I,
        nsimplify, nroots, ZZ, QQ
    )
    from sympy.abc import x as sym_x
    SYMPY_AVAILABLE = True
except ImportError:
    SYMPY_AVAILABLE = False

# NumPy für numerische Berechnungen
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# 1. Klassifikation von Zahlen
# ---------------------------------------------------------------------------

def is_algebraic_integer_check(x: float, max_degree: int = 6) -> bool:
    """
    Prüft, ob x eine algebraische ganze Zahl ist, d.h. ob x Nullstelle eines
    normierten ganzzahligen Polynoms p(t) = tⁿ + aₙ₋₁tⁿ⁻¹ + … + a₀ mit
    aᵢ ∈ ℤ ist.

    Algorithmus:
    Für Grad d von 1 bis max_degree werden alle Koeffizientenkombinationen
    aus [-max_coeff, max_coeff]ⁿ  (max_coeff intern = 5) ausprobiert.
    Das führende Koeffizient ist 1 (normiert).

    @param x          Die zu prüfende Zahl
    @param max_degree Maximaler Polynomgrad der Suche (Standard 6)
    @return           True wenn x wahrscheinlich algebraische ganze Zahl

    @timestamp 2026-03-10
    """
    # Rationale Zahlen sind immer algebraische ganze Zahlen, falls ihr
    # Nenner = 1, d.h. sie sind selbst ganze Zahlen
    if isinstance(x, int):
        return True

    max_coeff = 5  # Koeffizientenbereich beschränken
    tol = 1e-9

    for degree in range(1, max_degree + 1):
        # Iteriere über alle Koeffizientenkombinationen
        coeff_range = range(-max_coeff, max_coeff + 1)
        for coeffs in itertools.product(coeff_range, repeat=degree):
            # Normiertes Polynom: t^degree + c_{d-1}*t^{d-1} + ... + c_0
            # poly_coeffs: [a0, a1, ..., a_{d-1}, 1] (aufsteigend)
            poly = list(coeffs) + [1]  # [a0, ..., a_{d-1}, 1]
            # Auswerten: p(x) = 1*x^d + a_{d-1}*x^{d-1} + ... + a_0
            val = sum(poly[i] * (x ** i) for i in range(degree + 1))
            if abs(val) < tol:
                return True
    return False


def minimal_polynomial(x: float, max_degree: int = 6, max_coeff: int = 10) -> Optional[list]:
    """
    Sucht das Minimalpolynom von x: das normierte ganzzahlige Polynom
    kleinsten Grades mit x als Nullstelle.

    Verwendet LLL-artige Idee: für jeden Grad d wird geprüft, ob eine
    ganzzahlige Linearkombination der Potenzen 1, x, x², ..., x^d ≈ 0.

    @param x          Die algebraische Zahl
    @param max_degree Maximaler Grad der Suche
    @param max_coeff  Maximaler Betrag der Koeffizienten
    @return           Koeffizientenliste [a0, a1, ..., 1] oder None

    @timestamp 2026-03-10
    """
    tol = 1e-8
    for degree in range(1, max_degree + 1):
        # Normiertes Polynom: führender Koeffizient = 1
        coeff_range = range(-max_coeff, max_coeff + 1)
        best = None
        best_val = float('inf')
        for coeffs in itertools.product(coeff_range, repeat=degree):
            poly = list(coeffs) + [1]
            val = abs(sum(poly[i] * (x ** i) for i in range(degree + 1)))
            if val < best_val:
                best_val = val
                best = poly
            if val < tol:
                return poly
        # Falls kein exakter Treffer, gehe zum nächsten Grad
    return None  # kein Minimalpolynom gefunden → wahrscheinlich transzendent


def algebraic_degree(x: float, max_degree: int = 8) -> int:
    """
    Bestimmt den algebraischen Grad von x: den kleinsten Grad d, sodass x
    Nullstelle eines normierten ganzzahligen Polynoms vom Grad d ist.

    @param x          Die zu untersuchende Zahl
    @param max_degree Maximaler Suchgrad
    @return           Algebraischer Grad (1 = rational, >1 = irrational algebraisch,
                      max_degree+1 bedeutet wahrscheinlich transzendent)

    @timestamp 2026-03-10
    """
    tol = 1e-10  # Strenger: nur echte Rationale erkennen
    # Grad 1: x = p/q → rational (mit strenger Toleranz)
    frac = fractions.Fraction(x).limit_denominator(100)
    if abs(float(frac) - x) < tol:
        return 1

    for degree in range(2, max_degree + 1):
        max_coeff = 10
        for coeffs in itertools.product(range(-max_coeff, max_coeff + 1), repeat=degree):
            poly = list(coeffs) + [1]
            val = abs(sum(poly[i] * (x ** i) for i in range(degree + 1)))
            if val < tol:
                return degree
    return max_degree + 1  # wahrscheinlich transzendent


def classify_number(x: float) -> str:
    """
    Klassifiziert eine Zahl x als:
    - 'rational'               : x ∈ ℚ
    - 'algebraic'              : x ist algebraisch (Nullstelle eines rationalen Polynoms)
    - 'probably_transcendental': kein ganzzahliges Polynom kleinen Grades gefunden

    @param x Die zu klassifizierende Zahl
    @return  Klassifikationsstring

    @timestamp 2026-03-10
    """
    tol = 1e-9
    # Rationale Prüfung: Näherungsbruch
    frac = fractions.Fraction(x).limit_denominator(10000)
    if abs(float(frac) - x) < tol:
        return 'rational'

    # Algebraische Prüfung: Minimalpolynom suchen
    poly = minimal_polynomial(x, max_degree=6, max_coeff=8)
    if poly is not None:
        return 'algebraic'

    return 'probably_transcendental'


# ---------------------------------------------------------------------------
# 2. Liouville-Zahlen
# ---------------------------------------------------------------------------

def liouville_number(n_terms: int = 10) -> float:
    """
    Berechnet die Liouville-Zahl L = Σ_{k=1}^{n_terms} 10^{-k!}

    Diese Zahl ist transzendent: Sie wird durch rationale Zahlen
    „zu gut" approximiert (Irrationali­tätsexponent = ∞).

    @param n_terms Anzahl der Summanden (Standard 10)
    @return        Näherungswert der Liouville-Zahl

    @timestamp 2026-03-10
    """
    result = 0.0
    for k in range(1, n_terms + 1):
        # 10^{-k!} – für k ≥ 16 ist k! > 308 → float-Unterlauf, OK in Python
        factorial_k = math.factorial(k)
        if factorial_k > 300:
            break  # float-Unterlauf: Beitrag vernachlässigbar
        result += 10.0 ** (-factorial_k)
    return result


def is_liouville_like(x: float, precision: int = 50) -> bool:
    """
    Einfache Heuristik: Prüft, ob x Liouville-artig ist, d.h. ob der
    Irrationali­tätsexponent sehr groß erscheint (> precision).

    @param x         Die zu prüfende Zahl
    @param precision Schwellenwert für den Exponenten
    @return          True wenn x wahrscheinlich Liouville-artig

    @timestamp 2026-03-10
    """
    exp_val = liouville_approximation_exponent(x, max_q=500)
    return exp_val >= precision


def liouville_approximation_exponent(x: float, max_q: int = 1000) -> float:
    """
    Schätzt den Irrationali­tätsexponenten μ(x) numerisch:
    μ(x) = sup{μ : |x - p/q| < 1/q^μ hat unendlich viele Lösungen p/q ∈ ℚ}

    Methode: Für alle q bis max_q wird der beste Approximationsbruch p/q
    gesucht und der Exponent μ = log(1/|x-p/q|) / log(q) berechnet.

    @param x      Die zu untersuchende Zahl
    @param max_q  Maximaler Nenner (Standard 1000)
    @return       Geschätzter Irrationali­tätsexponent

    @timestamp 2026-03-10
    """
    max_exponent = 2.0  # Ausgangswert (Hurwitz: mindestens 2)
    for q in range(2, max_q + 1):
        # Bester ganzzahliger Zähler p
        p = round(x * q)
        diff = abs(x - p / q)
        if diff < 1e-15:
            continue  # Vermeidung von log(0)
        # μ aus |x - p/q| ≈ 1/q^μ → μ = -log(diff)/log(q)
        mu = -math.log(diff) / math.log(q)
        if mu > max_exponent:
            max_exponent = mu
    return max_exponent


# ---------------------------------------------------------------------------
# 3. Lindemann-Weierstrass und verwandte Sätze (Demonstrationen)
# ---------------------------------------------------------------------------

def lindemann_weierstrass_demo() -> dict:
    """
    Demonstration des Lindemann-Weierstrass-Satzes:
    Sind α₁, …, αₙ verschiedene algebraische Zahlen, so sind
    e^α₁, …, e^αₙ linear unabhängig über ℚ.

    Konkret: Für jedes algebraische α ≠ 0 ist e^α transzendent.

    @return Wörterbuch mit Beispielen und Erklärungen

    @timestamp 2026-03-10
    """
    demo = {
        'satz': (
            "Lindemann-Weierstrass (1882/1885): Sind α₁,…,αₙ ∈ ℚ̄ paarweise verschieden, "
            "so sind e^α₁,…,e^αₙ linear unabhängig über ℚ̄. "
            "Insbesondere: e^α ∈ ℚ̄ ⟹ α = 0."
        ),
        'beispiele': {
            'e^1 = e': {
                'alpha': 1,
                'alpha_algebraisch': True,
                'wert': math.e,
                'transzendent': True,
                'erklaerung': 'α=1 algebraisch, also e^1=e transzendent (Hermite 1873).'
            },
            'e^(i*pi) = -1': {
                'alpha': '±iπ',
                'erklaerung': (
                    "e^(iπ) = -1 ist algebraisch. Nach L-W muss iπ nicht-algebraisch sein "
                    "oder α=0 – da iπ ≠ 0, folgt: π transzendent."
                ),
                'pi_transzendent': True
            },
            'e^sqrt(2)': {
                'alpha': math.sqrt(2),
                'alpha_algebraisch': True,
                'wert': math.exp(math.sqrt(2)),
                'transzendent': True,
                'erklaerung': '√2 ist algebraisch (Minpoly: x²-2), also e^√2 transzendent.'
            },
        },
        'folgerungen': [
            'e ist transzendent (Hermite, 1873)',
            'π ist transzendent (Lindemann, 1882)',
            'sin(1), cos(1), tan(1) sind transzendent',
            'log(2), log(3), ln(k) für k∈ℤ, k>0, k≠1 sind transzendent',
        ]
    }
    return demo


def e_is_transcendental_demo() -> dict:
    """
    Demonstration, dass e transzendent ist (Hermite 1873).

    Idee: Annahme e erfüllt aₙeⁿ + … + a₀ = 0 (aᵢ ∈ ℤ).
    Durch Konstruktion eines Hilfsintegrals I(t) = ∫₀ᵗ e^{t-u} f(u) du
    führt dies zu einem Widerspruch (ganzzahliger Wert ≠ 0, aber beliebig klein).

    @return Wörterbuch mit Erklärungen und numerischen Verifikationen

    @timestamp 2026-03-10
    """
    # Numerische Verifikation: kein Polynom kleinen Grades hat e als Nullstelle
    e_val = math.e
    poly_evals = {}
    for degree in range(1, 7):
        # Bestes normiertes ganzzahliges Polynom vom Grad `degree`
        best_val = float('inf')
        best_poly = None
        for coeffs in itertools.product(range(-5, 6), repeat=degree):
            poly = list(coeffs) + [1]
            val = abs(sum(poly[i] * (e_val ** i) for i in range(degree + 1)))
            if val < best_val:
                best_val = val
                best_poly = poly
        poly_evals[degree] = {'bester_wert': best_val, 'polynom': best_poly}

    return {
        'satz': 'e ist transzendent (Hermite, 1873).',
        'beweis_skizze': (
            "Annahme: e ist algebraisch, d.h. Σ aₖ eᵏ = 0. "
            "Konstruiere f(x) = xᵖ⁻¹(x-1)ᵖ…(x-n)ᵖ / (p-1)! für große Primzahl p. "
            "Das Integral J = Σ aₖ Σ_{j=0}^{k} e^{k-j} f^{(j)}(0) ist ganz "
            "und gleichzeitig |J| < 1 für großes p → Widerspruch."
        ),
        'numerische_verifikation': poly_evals,
        'schlussfolgerung': (
            "Kein normiertes ganzzahliges Polynom mit Koeffizienten in [-5,5] "
            "und Grad ≤ 6 hat e≈2.71828... als Nullstelle."
        )
    }


def pi_is_transcendental_demo() -> dict:
    """
    Demonstration, dass π transzendent ist (Lindemann 1882).

    Beweis-Idee: e^{iπ} = -1 ∈ ℚ̄. Nach Lindemann-Weierstrass gilt:
    Falls α algebraisch und α ≠ 0, dann ist e^α transzendent.
    Da e^{iπ} = -1 algebraisch ist, muss iπ = 0 sein – Widerspruch,
    also ist π (und damit iπ) nicht algebraisch → π transzendent.

    @return Wörterbuch mit Erklärungen

    @timestamp 2026-03-10
    """
    pi_val = math.pi
    # Beste Polynomannäherung suchen (numerisch)
    best_poly_vals = {}
    for degree in range(1, 7):
        best_val = float('inf')
        best_poly = None
        for coeffs in itertools.product(range(-5, 6), repeat=degree):
            poly = list(coeffs) + [1]
            val = abs(sum(poly[i] * (pi_val ** i) for i in range(degree + 1)))
            if val < best_val:
                best_val = val
                best_poly = poly
        best_poly_vals[degree] = {'bester_wert': best_val, 'polynom': best_poly}

    return {
        'satz': 'π ist transzendent (Lindemann, 1882).',
        'beweis_idee': (
            "e^{iπ} + 1 = 0 (Euler). Da -1 algebraisch ist, wäre iπ nach L-W "
            "nicht-algebraisch (sonst e^{iπ} transzendent). "
            "Da i algebraisch (Minpoly: x²+1), muss π transzendent sein."
        ),
        'folgerungen': [
            'Quadratur des Kreises unmöglich (mit Zirkel und Lineal)',
            'cos(r), sin(r) für r∈ℚ, r≠0 sind transzendent',
        ],
        'numerische_verifikation': best_poly_vals,
        'pi_wert': pi_val
    }


def gelfond_schneider_demo() -> dict:
    """
    Demonstration des Gelfond-Schneider-Satzes (1934):
    Sind a, b algebraisch, a ∉ {0,1}, b ∉ ℚ, dann ist a^b transzendent.

    Löst Hilberts 7. Problem (1900).

    @return Wörterbuch mit Beispielen

    @timestamp 2026-03-10
    """
    gelfond_const = 2 ** math.sqrt(2)  # Gelfond-Schneider-Konstante ≈ 2.6651...
    return {
        'satz': (
            "Gelfond-Schneider (1934): Sind a,b ∈ ℚ̄, a ∉ {0,1}, b ∉ ℚ, "
            "dann ist a^b transzendent."
        ),
        'hilbert_problem': 'Löst Hilberts 7. Problem aus dem Jahr 1900.',
        'beispiele': {
            '2^sqrt(2)': {
                'a': 2, 'b': math.sqrt(2),
                'wert': gelfond_const,
                'transzendent': True,
                'erklaerung': 'a=2 algebraisch, b=√2 irrational algebraisch → 2^√2 transzendent'
            },
            'e^pi (Gelfond-Konstante)': {
                'erklaerung': (
                    "e^π = (e^{iπ})^{-i} = (-1)^{-i}. "
                    "Mit a=-1, b=-i (irrational algebraisch) → e^π transzendent."
                ),
                'wert': math.exp(math.pi),
                'transzendent': True
            },
            'i^i': {
                'erklaerung': (
                    "i^i = e^{i·log(i)} = e^{i·(iπ/2)} = e^{-π/2}. "
                    "Da e^{-π/2} = e^{π/2}⁻¹ und e^π transzendent, ist i^i transzendent."
                ),
                'wert': (1j) ** (1j),
                'transzendent': True
            }
        }
    }


def baker_theorem_demo() -> dict:
    """
    Demonstration von Bakers Theorem (1966):
    Logarithmische Linearkombinationen algebraischer Zahlen.

    Baker: Sind α₁,…,αₙ algebraische Zahlen ≠ 0,1, und
    β₁,…,βₙ algebraisch mit β₁ log α₁ + … + βₙ log αₙ ≠ 0,
    dann ist β₁ log α₁ + … + βₙ log αₙ transzendent.

    @return Wörterbuch mit Erklärungen und Beispielen

    @timestamp 2026-03-10
    """
    return {
        'satz': (
            "Baker (1966): Sind α₁,…,αₙ algebraisch ≠ 0,1, und β₀,β₁,…,βₙ algebraisch, "
            "nicht alle 0, dann ist β₀ + β₁·log(α₁) + … + βₙ·log(αₙ) ≠ 0 "
            "(falls diese Summe nicht verschwindet – und sie ist transzendent)."
        ),
        'bedeutung': [
            'Effektive Lösung von Thue-Mahler-Gleichungen',
            'Schranken für elliptische Logarithmen',
            'Anwendung in der Kryptographie (Diskrete-Logarithmus-Probleme)',
            'Beweis der Endlichkeit ganzzahliger Punkte auf elliptischen Kurven',
        ],
        'beispiele': {
            'log(2) + log(3) = log(6)': {
                'beta': [1, 1], 'alpha': [2, 3],
                'wert': math.log(2) + math.log(3),
                'gleich_log6': abs(math.log(2) + math.log(3) - math.log(6)) < 1e-12,
                'transzendent': True,
                'erklaerung': 'log(2), log(3), log(6) sind alle transzendent (L-W)'
            },
            'pi = -i*log(-1)': {
                'erklaerung': 'π = Im(log(-1)) ist transzendent nach L-W',
                'wert': math.pi
            }
        },
        'fields_medal': 'Alan Baker erhielt 1970 die Fields-Medaille für diesen Satz.'
    }


# ---------------------------------------------------------------------------
# 4. Maßzahlen der Transzendenz
# ---------------------------------------------------------------------------

def irrationality_measure(x: float, max_q: int = 10000) -> float:
    """
    Schätzt das Irrationali­tätsmaß μ(x) numerisch:
    μ(x) = sup{μ : |x - p/q| < q^{-μ} für unendlich viele p/q ∈ ℚ}

    Jede irrationale Zahl hat μ(x) ≥ 2 (Dirichlet).
    Rationale Zahlen: μ(x) = 1.
    Liouville-Zahlen: μ(x) = ∞.

    @param x      Zu untersuchende Zahl
    @param max_q  Maximaler Nenner (Standard 10000)
    @return       Geschätztes Irrationali­tätsmaß

    @timestamp 2026-03-10
    """
    return liouville_approximation_exponent(x, max_q=max_q)


def known_irrationality_measures() -> dict:
    """
    Gibt bekannte Irrationali­tätsmaße wichtiger Konstanten zurück.

    @return Wörterbuch mit Konstante → (untere Schranke, Quelle)

    @timestamp 2026-03-10
    """
    return {
        'pi': {
            'untere_schranke': 7.103,
            'quelle': 'Salikhov (2008)',
            'erklaerung': 'μ(π) ≥ 7.103205...'
        },
        'e': {
            'untere_schranke': 2.0,
            'exakt': True,
            'quelle': 'Klassisch',
            'erklaerung': 'μ(e) = 2, e ist kein Liouville-Zahl'
        },
        'ln(2)': {
            'untere_schranke': 3.57,
            'quelle': 'Rukhadze (1987)',
            'erklaerung': 'μ(ln 2) ≥ 3.5702...'
        },
        'sqrt(2)': {
            'untere_schranke': 2.0,
            'exakt': True,
            'quelle': 'Klassisch (Lagrange)',
            'erklaerung': 'μ(√2) = 2, da Kettenbruch [1;2,2,2,...] periodisch'
        },
        'zeta(3)': {
            'untere_schranke': 5.51,
            'quelle': 'Rhin-Viola (2001)',
            'erklaerung': "μ(ζ(3)) ≥ 5.513891... (Apéry's Konstante)"
        },
        'liouville': {
            'untere_schranke': float('inf'),
            'quelle': 'Definition',
            'erklaerung': 'μ(L) = ∞ für Liouville-Zahlen'
        }
    }


def diophantine_approximation_quality(x: float, max_q: int = 1000) -> dict:
    """
    Analysiert die diophantische Approximationsqualität von x nach Hurwitz:
    |x - p/q| < 1/(√5 · q²) hat unendlich viele Lösungen (optimal).

    @param x      Zu untersuchende Zahl
    @param max_q  Maximaler Nenner
    @return       Wörterbuch mit Analyse-Ergebnissen

    @timestamp 2026-03-10
    """
    hurwitz_const = math.sqrt(5)
    hurwitz_violations = 0  # Brüche, die Hurwitz-Schranke unterschreiten
    best_approx = None
    best_error = float('inf')
    results = []

    for q in range(1, max_q + 1):
        p = round(x * q)
        error = abs(x - p / q)
        hurwitz_bound = 1.0 / (hurwitz_const * q * q)

        if error < hurwitz_bound:
            hurwitz_violations += 1

        if error < best_error and error > 0:
            best_error = error
            best_approx = (p, q, error)

        if error < 1e-10:
            results.append({'p': p, 'q': q, 'fehler': error, 'hurwitz': error < hurwitz_bound})

    return {
        'x': x,
        'beste_annaeherung': best_approx,
        'hurwitz_treffer': hurwitz_violations,
        'hurwitz_konstante': hurwitz_const,
        'gute_annaeherungen': results[:10]  # Top 10
    }


def continued_fraction_irrationality(x: float, n_terms: int = 20) -> dict:
    """
    Berechnet die Kettenbruchentwicklung [a₀; a₁, a₂, …] von x und schätzt
    das Irrationali­tätsmaß anhand der Größe der Kettenbruchkoeffizienten.

    Großen aₖ entsprechen gute Näherungen und größeren Irrationali­tätsmaß-Beiträgen.

    @param x       Zu untersuchende Zahl
    @param n_terms Anzahl der Kettenbruchglieder
    @return        Wörterbuch mit Kettenbruch und Analyse

    @timestamp 2026-03-10
    """
    coeffs = []
    convergents = []
    val = x

    for _ in range(n_terms):
        a = int(val)  # Ganzteilanteil
        coeffs.append(a)
        frac_part = val - a

        # Konvergente berechnen
        if len(coeffs) == 1:
            p, q = a, 1
        elif len(coeffs) == 2:
            p = coeffs[0] * coeffs[1] + 1
            q = coeffs[1]
        else:
            # Rekursion: pₙ = aₙ·pₙ₋₁ + pₙ₋₂
            p = coeffs[-1] * convergents[-1][0] + convergents[-2][0]
            q = coeffs[-1] * convergents[-1][1] + convergents[-2][1]

        error = abs(x - p / q) if q != 0 else float('inf')
        convergents.append((p, q, error))

        if frac_part < 1e-15:
            break  # x ist rational oder Präzision erschöpft
        val = 1.0 / frac_part

    # Irrationali­tätsmaß-Schätzung über große Koeffizienten
    # μ ≈ 2 + lim sup log(aₙ₊₁) / log(qₙ)
    mu_estimates = []
    for i in range(len(coeffs) - 1):
        if coeffs[i + 1] > 1 and convergents[i][1] > 1:
            mu_i = 2 + math.log(coeffs[i + 1]) / math.log(convergents[i][1])
            mu_estimates.append(mu_i)

    irrationality_est = max(mu_estimates) if mu_estimates else 2.0

    return {
        'x': x,
        'kettenbruch': coeffs,
        'konvergente': [(p, q) for p, q, _ in convergents],
        'fehler': [e for _, _, e in convergents],
        'groesster_koeff': max(coeffs[1:]) if len(coeffs) > 1 else 0,
        'irrationality_schaetzung': irrationality_est
    }


# ---------------------------------------------------------------------------
# 5. Algebraische Zahlen
# ---------------------------------------------------------------------------

class AlgebraicNumber:
    """
    Klasse zur Darstellung einer algebraischen Zahl α, definiert durch ihr
    Minimalpolynom p(t) = aₙtⁿ + … + a₁t + a₀ mit aᵢ ∈ ℤ, aₙ ≠ 0.

    @author    Kurt Ingwer
    @version   1.0.0
    @timestamp 2026-03-10
    """

    def __init__(self, poly_coeffs: list):
        """
        Initialisiert eine algebraische Zahl durch ihr Minimalpolynom.

        @param poly_coeffs Koeffizientenliste [a0, a1, ..., an] (aufsteigend in t),
                          z.B. [−2, 0, 1] für t² − 2 (√2)

        @timestamp 2026-03-10
        """
        if len(poly_coeffs) < 2:
            raise ValueError("Polynom muss mindestens Grad 1 haben.")
        self.poly_coeffs = list(poly_coeffs)
        self._roots = None  # Lazy evaluation

    def degree(self) -> int:
        """
        Gibt den algebraischen Grad zurück (= Grad des Minimalpolynoms).

        @return Algebraischer Grad

        @timestamp 2026-03-10
        """
        return len(self.poly_coeffs) - 1

    def _compute_roots(self) -> list:
        """
        Berechnet alle Nullstellen (Konjugierte) des Minimalpolynoms numerisch.

        @return Liste der komplexen Nullstellen

        @timestamp 2026-03-10
        """
        if self._roots is not None:
            return self._roots
        # NumPy: numpy.roots erwartet Koeffizienten absteigend
        coeffs_desc = list(reversed(self.poly_coeffs))
        if NUMPY_AVAILABLE:
            self._roots = list(np.roots(coeffs_desc))
        else:
            # Einfache Bisektion für Grad 1 und 2 (Fallback)
            if self.degree() == 1:
                # a1*t + a0 = 0 → t = -a0/a1
                self._roots = [-self.poly_coeffs[0] / self.poly_coeffs[1]]
            elif self.degree() == 2:
                a, b, c = self.poly_coeffs[2], self.poly_coeffs[1], self.poly_coeffs[0]
                disc = b * b - 4 * a * c
                self._roots = [
                    (-b + cmath.sqrt(disc)) / (2 * a),
                    (-b - cmath.sqrt(disc)) / (2 * a)
                ]
            else:
                self._roots = []
        return self._roots

    def conjugates(self) -> list:
        """
        Gibt alle Konjugierten (Nullstellen des Minimalpolynoms) zurück.

        @return Liste komplexer Konjugierter

        @timestamp 2026-03-10
        """
        return self._compute_roots()

    def norm(self) -> complex:
        """
        Berechnet die Norm N(α) = ∏ σ(α), Produkt über alle Konjugierte.

        Für normiertes Polynom gilt: N(α) = (-1)^d · a₀ (bis auf Vorzeichen).

        @return Norm der algebraischen Zahl

        @timestamp 2026-03-10
        """
        roots = self._compute_roots()
        if not roots:
            return complex(self.poly_coeffs[0])
        result = complex(1.0)
        for r in roots:
            result *= r
        return result

    def is_unit_root(self, tol: float = 1e-8) -> bool:
        """
        Prüft, ob α eine Einheitswurzel ist, d.h. ob |α| = 1 für alle Konjugierten
        und α^n = 1 für ein n ∈ ℕ.

        @param tol Toleranz für numerischen Vergleich
        @return    True wenn α Einheitswurzel

        @timestamp 2026-03-10
        """
        roots = self._compute_roots()
        if not roots:
            return False
        # Alle Konjugierte müssen Betrag 1 haben
        for r in roots:
            if abs(abs(r) - 1.0) > tol:
                return False
        # Prüfe ob eine Potenz ≈ 1
        root = roots[0]
        for n in range(1, 101):
            if abs(root ** n - 1.0) < tol:
                return True
        return False

    def __repr__(self) -> str:
        """Textdarstellung der algebraischen Zahl."""
        return f"AlgebraicNumber(poly={self.poly_coeffs})"


def cyclotomic_polynomial(n: int) -> list:
    """
    Berechnet das n-te Kreisteilungspolynom Φₙ(x).

    Φₙ(x) = ∏_{1≤k≤n, gcd(k,n)=1} (x - e^{2πik/n})

    Für kleine n direkt über die Rekursion:
    x^n - 1 = ∏_{d|n} Φ_d(x)

    @param n Ganzzahl ≥ 1
    @return  Koeffizientenliste [a0, a1, ..., aφ(n)] des Kreisteilungspolynoms

    @timestamp 2026-03-10
    """
    if SYMPY_AVAILABLE:
        # SymPy berechnet Kreisteilungspolynome exakt
        sym_x = Symbol('x')
        poly = cyclotomic_poly(n, sym_x)
        # Konvertiere zu Koeffizientenliste (aufsteigend)
        poly_obj = Poly(poly, sym_x)
        coeffs_desc = [int(c) for c in poly_obj.all_coeffs()]
        return list(reversed(coeffs_desc))
    else:
        # Fallback: direkte Berechnung über Primitivwurzeln
        euler_phi = sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)
        # Für n=1: Φ₁(x) = x - 1
        if n == 1:
            return [-1, 1]
        # Für Primzahlen p: Φ_p(x) = 1 + x + x² + ... + x^{p-1}
        if all(n % i != 0 for i in range(2, int(math.sqrt(n)) + 1)) and n > 1:
            return [1] * n  # [1, 1, 1, ..., 1] (n Terme)
        return [1]  # Fallback


def number_field_degree(poly_coeffs: list) -> int:
    """
    Bestimmt den Grad des Zahlkörpers ℚ(α), wobei α Nullstelle des Polynoms
    mit Koeffizienten poly_coeffs ist.

    Der Grad [ℚ(α):ℚ] = Grad des Minimalpolynoms über ℚ.

    @param poly_coeffs Koeffizientenliste [a0, a1, ..., an]
    @return            Grad des Zahlkörpers

    @timestamp 2026-03-10
    """
    # Der Grad des Zahlkörpers ist der Grad des irreduziblen Minimalpolynoms
    return len(poly_coeffs) - 1
