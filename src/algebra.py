"""
@file algebra.py
@brief Thin-Wrapper: Importiert alle Klassen für Rückwärtskompatibilität.
@description
    Dieser Wrapper stellt alle Klassen und Funktionen bereit, die vorher direkt in
    algebra.py definiert waren. Code der `from algebra import Polynomial`
    verwendet, muss nicht geändert werden.

    Sub-Module:
    - algebra_core.py         → Polynomial, solve_linear, solve_quadratic,
                                gcd, lcm, extended_gcd, mod_inverse
    - algebra_numbertheory.py → is_prime, prime_factorization, euler_phi,
                                rsa_keygen, rsa_encrypt, rsa_decrypt
    - algebra_diophantine.py  → solve_linear_diophantine, solve_quadratic_diophantine_pell,
                                solve_pythagorean_triples, solve_diophantine_two_squares,
                                markov_numbers, is_quadratic_residue, quadratic_residues,
                                quadratic_reciprocity, tonelli_shanks, cipolla_algorithm

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

# Rückwärtskompatibilität: Re-exportiere alle Klassen und Funktionen
from algebra_core import Polynomial, solve_linear, solve_quadratic, gcd, lcm, extended_gcd, mod_inverse
from algebra_numbertheory import is_prime, prime_factorization, euler_phi, rsa_keygen, rsa_encrypt, rsa_decrypt
from algebra_diophantine import (
    solve_linear_diophantine, solve_quadratic_diophantine_pell,
    solve_pythagorean_triples, solve_diophantine_two_squares, markov_numbers,
    is_quadratic_residue, quadratic_residues, quadratic_reciprocity,
    tonelli_shanks, cipolla_algorithm
)

# Für `from algebra import *`
__all__ = [
    'Polynomial', 'solve_linear', 'solve_quadratic', 'gcd', 'lcm',
    'extended_gcd', 'mod_inverse', 'is_prime', 'prime_factorization',
    'euler_phi', 'rsa_keygen', 'rsa_encrypt', 'rsa_decrypt',
    'solve_linear_diophantine', 'solve_quadratic_diophantine_pell',
    'solve_pythagorean_triples', 'solve_diophantine_two_squares', 'markov_numbers',
    'is_quadratic_residue', 'quadratic_residues', 'quadratic_reciprocity',
    'tonelli_shanks', 'cipolla_algorithm'
]
