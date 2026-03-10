"""
@file linear_algebra.py
@brief Thin-Wrapper: Importiert alle Klassen aus den Sub-Modulen für Rückwärtskompatibilität.
@description
    Dieser Wrapper stellt alle Klassen und Funktionen bereit, die vorher direkt in
    linear_algebra.py definiert waren. Code der `from linear_algebra import Vector`
    verwendet, muss nicht geändert werden.

    Sub-Module:
    - vectors.py      → Vector, gram_schmidt
    - matrix_ops.py   → Matrix
    - matrix_decomp.py → lu_decomposition, qr_decomposition, svd, matrix_rank,
                         condition_number, givens_rotation_matrix,
                         givens_qr_decomposition, givens_solve_least_squares

@author Kurt Ingwer
@lastModified 2026-03-10
"""

# Rückwärtskompatibilität: Re-exportiere alle Klassen und Funktionen
from vectors import Vector, gram_schmidt
from matrix_ops import Matrix
from matrix_decomp import (
    lu_decomposition, qr_decomposition, svd, matrix_rank, condition_number,
    givens_rotation_matrix, givens_qr_decomposition, givens_solve_least_squares
)

# Für `from linear_algebra import *`
__all__ = [
    'Vector', 'Matrix', 'gram_schmidt',
    'lu_decomposition', 'qr_decomposition', 'svd', 'matrix_rank', 'condition_number',
    'givens_rotation_matrix', 'givens_qr_decomposition', 'givens_solve_least_squares'
]
