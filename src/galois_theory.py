"""
@file galois_theory.py
@brief Galois-Theorie – Rückwärtskompatibilitäts-Wrapper.
@description
    Dieses Modul ist ein Kompatibilitäts-Wrapper, der alle Symbole aus den drei
    spezialisierten galois_theory_*-Modulen re-exportiert.

    Die ursprüngliche galois_theory.py wurde aufgeteilt in:

    1. **galois_theory_fields.py** – Körper und Körpererweiterungen:
       - _poly_mul_mod, _poly_pow_mod, _poly_eval, _is_irreducible_over_fp, _find_irreducible_poly
       - FiniteField, FieldExtension, GaloisGroup
       - finite_field(), cyclotomic_polynomial(), cyclotomic_galois_group()
       - galois_group_finite_field(), kronecker_weber_check()
       - dirichlet_characters_from_galois(), minimal_polynomial(), splitting_field()
       - norm_and_trace(), hilbert90(), primitive_element_theorem()

    2. **galois_theory_polynomials.py** – Polynome und Auflösbarkeit:
       - _discriminant_degree2, _is_perfect_square, _cubic_resolvent_discriminant
       - discriminant_polynomial(), galois_group_polynomial()
       - is_solvable_by_radicals(), galois_correspondence()
       - galois_group_of_polynomial(), radical_tower(), abel_ruffini_demo()
       - galois_group_symmetric(), fundamental_theorem_verify()
       - finite_field_discrete_log()

    3. **galois_theory_constructibility.py** – Konstruierbarkeit mit Zirkel und Lineal:
       - FERMAT_PRIMES
       - fermat_prime_check(), is_constructible(), construct_regular_polygon()
       - gauss_wantzel_theorem(), impossible_constructions()
       - constructible_numbers(), degree_tower_decomposition()

    **Verwendung:**
    Bestehender Code, der `from galois_theory import ...` verwendet, funktioniert
    ohne Änderungen weiter, da alle Symbole hier re-exportiert werden.

@author Michael Fuhrmann
@version 2.0
@since 2026-03-11
@lastModified 2026-03-11
"""

# =============================================================================
# RE-EXPORT: galois_theory_fields
# =============================================================================

from galois_theory_fields import (
    # Hilfsfunktionen für Polynomrechnung in GF(p)
    _poly_mul_mod,
    _poly_pow_mod,
    _poly_eval,
    _is_irreducible_over_fp,
    _find_irreducible_poly,

    # Klassen
    FiniteField,
    FieldExtension,
    GaloisGroup,

    # Fabriken und Berechnungsfunktionen
    finite_field,
    cyclotomic_polynomial,
    cyclotomic_galois_group,
    galois_group_finite_field,
    kronecker_weber_check,
    dirichlet_characters_from_galois,
    minimal_polynomial,
    splitting_field,
    norm_and_trace,
    hilbert90,
    primitive_element_theorem,
)

# =============================================================================
# RE-EXPORT: galois_theory_polynomials
# =============================================================================

from galois_theory_polynomials import (
    # Hilfsfunktionen
    _discriminant_degree2,
    _is_perfect_square,
    _cubic_resolvent_discriminant,

    # Hauptfunktionen
    discriminant_polynomial,
    galois_group_polynomial,
    is_solvable_by_radicals,
    galois_correspondence,
    galois_group_of_polynomial,
    radical_tower,
    abel_ruffini_demo,
    galois_group_symmetric,
    fundamental_theorem_verify,
    finite_field_discrete_log,
)

# =============================================================================
# RE-EXPORT: galois_theory_constructibility
# =============================================================================

from galois_theory_constructibility import (
    # Konstanten
    FERMAT_PRIMES,

    # Funktionen
    fermat_prime_check,
    is_constructible,
    construct_regular_polygon,
    gauss_wantzel_theorem,
    impossible_constructions,
    constructible_numbers,
    degree_tower_decomposition,
)

# =============================================================================
# ÖFFENTLICHE API (für `from galois_theory import *`)
# =============================================================================

__all__ = [
    # Hilfsfunktionen Felder
    '_poly_mul_mod', '_poly_pow_mod', '_poly_eval',
    '_is_irreducible_over_fp', '_find_irreducible_poly',

    # Klassen
    'FiniteField', 'FieldExtension', 'GaloisGroup',

    # Körper-Funktionen
    'finite_field', 'cyclotomic_polynomial', 'cyclotomic_galois_group',
    'galois_group_finite_field', 'kronecker_weber_check',
    'dirichlet_characters_from_galois', 'minimal_polynomial',
    'splitting_field', 'norm_and_trace', 'hilbert90',
    'primitive_element_theorem',

    # Polynom-Hilfsfunktionen
    '_discriminant_degree2', '_is_perfect_square', '_cubic_resolvent_discriminant',

    # Polynom-Funktionen
    'discriminant_polynomial', 'galois_group_polynomial',
    'is_solvable_by_radicals', 'galois_correspondence',
    'galois_group_of_polynomial', 'radical_tower', 'abel_ruffini_demo',
    'galois_group_symmetric', 'fundamental_theorem_verify',
    'finite_field_discrete_log',

    # Konstruierbarkeits-Konstanten und -Funktionen
    'FERMAT_PRIMES',
    'fermat_prime_check', 'is_constructible', 'construct_regular_polygon',
    'gauss_wantzel_theorem', 'impossible_constructions',
    'constructible_numbers', 'degree_tower_decomposition',
]
