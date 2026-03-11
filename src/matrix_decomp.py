"""
@file matrix_decomp.py
@brief Matrix-Zerlegungen: LU, QR (Householder & Givens), SVD, Rang, Konditionszahl.
@description
    Enthält Algorithmen für Matrixzerlegungen und verwandte Operationen:

    - lu_decomposition()         – LU-Zerlegung (Doolittle mit Teilpivotisierung)
    - qr_decomposition()         – QR-Zerlegung via Householder-Reflexionen
    - svd()                      – Singulärwertzerlegung (via numpy/LAPACK)
    - matrix_rank()              – Matrixrang via SVD
    - condition_number()         – Konditionszahl κ = σ_max / σ_min
    - givens_rotation_matrix()   – Givens-Rotationsmatrix G(i,j,θ)
    - givens_qr_decomposition()  – QR via Givens-Rotationen
    - givens_solve_least_squares() – Kleinste-Quadrate via Givens-QR

    Importiert Matrix aus matrix_ops.py und Vector aus vectors.py.
    Ausgelagert aus linear_algebra.py für bessere Modularität.

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import math
import numpy as np

# Sub-Module importieren
from vectors import Vector
from matrix_ops import Matrix

# Spezifische mathematische Ausnahmen importieren
from exceptions import SingularMatrixError


# =============================================================================
# LU-ZERLEGUNG (Doolittle mit Teilpivotisierung)
# =============================================================================

def lu_decomposition(matrix: 'Matrix') -> tuple['Matrix', 'Matrix', 'Matrix', int]:
    """
    LU-Zerlegung mit Teilpivotisierung (Doolittle-Algorithmus).

    Jede reguläre Matrix A lässt sich zerlegen in:
        P·A = L·U

    Wobei:
        P – Permutationsmatrix (Zeilenvertauschungen)
        L – untere Dreiecksmatrix (lower triangular), Diagonale = 1
        U – obere Dreiecksmatrix (upper triangular)

    Teilpivotisierung: Pro Spalte wird das betragsmäßig größte Element
    als Pivot gewählt → numerisch stabiler als ohne Pivotisierung.

    Anwendung:
        LGS lösen: Ax = b → LUx = Pb → Lz = Pb (Vorwärts), Ux = z (Rückwärts)
        Determinante: det(A) = det(U) · sign(P) = Π U[i,i] · (-1)^swaps

    @param matrix: Quadratische Matrix A (n×n)
    @return: (L, U, P, n_swaps) – untere/obere Dreiecksmatrix, Permutation, Anzahl Tausche
    @raises ValueError: Wenn Matrix nicht quadratisch oder singulär
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    n = matrix.rows
    if n != matrix.cols:
        raise ValueError(f"LU-Zerlegung erfordert quadratische Matrix, erhalten: {n}×{matrix.cols}")

    # Arbeitskopieen der Matrix-Daten
    A = [list(matrix._data[i]) for i in range(n)]
    # Permutationsvektor: P[i] = j bedeutet Zeile i des Originals ist jetzt Zeile j
    perm = list(range(n))
    n_swaps = 0  # Anzahl Zeilenvertauschungen (für Vorzeichen der Determinante)

    for k in range(n):
        # Teilpivotisierung: größtes Element in Spalte k (ab Zeile k) finden
        max_val = abs(A[k][k])
        max_row = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > max_val:
                max_val = abs(A[i][k])
                max_row = i

        # Zeilen tauschen (Permutation)
        if max_row != k:
            A[k], A[max_row] = A[max_row], A[k]
            perm[k], perm[max_row] = perm[max_row], perm[k]
            n_swaps += 1

        if abs(A[k][k]) < 1e-14:
            # Pivot nahe 0 bedeutet singuläre Matrix → LU-Zerlegung nicht möglich
            raise SingularMatrixError("LU-Zerlegung (Doolittle)")

        # Elimination: L[i,k] = A[i,k] / A[k,k], dann A[i,:] -= L[i,k] * A[k,:]
        for i in range(k + 1, n):
            factor = A[i][k] / A[k][k]
            A[i][k] = factor  # L-Eintrag im unteren Teil speichern
            for j in range(k + 1, n):
                A[i][j] -= factor * A[k][j]

    # L und U aus kombinierter Matrix extrahieren
    L_data = [[0.0] * n for _ in range(n)]
    U_data = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                L_data[i][j] = A[i][j]  # Unterhalb der Diagonale
            elif i == j:
                L_data[i][j] = 1.0       # L-Diagonale ist 1 (Doolittle)
                U_data[i][j] = A[i][j]  # U-Diagonale
            else:
                U_data[i][j] = A[i][j]  # Oberhalb der Diagonale

    # Permutationsmatrix aus perm-Vektor aufbauen
    P_data = [[0.0] * n for _ in range(n)]
    for i in range(n):
        P_data[i][perm[i]] = 1.0

    return Matrix(L_data), Matrix(U_data), Matrix(P_data), n_swaps


# =============================================================================
# QR-ZERLEGUNG (Householder-Reflexionen)
# =============================================================================

def qr_decomposition(matrix: 'Matrix') -> tuple['Matrix', 'Matrix']:
    """
    QR-Zerlegung via Householder-Reflexionen.

    Jede Matrix A (m×n, m ≥ n) lässt sich zerlegen in:
        A = Q·R

    Wobei:
        Q – unitäre/orthogonale Matrix (Q^T·Q = I), m×m
        R – obere Dreiecksmatrix, m×n

    Householder-Reflexion:
        H(v) = I - 2·v·vᵀ / (vᵀ·v)

    Jede Reflexion macht eine Spalte unterhalb der Diagonale null.
    Numerisch stabiler als Gram-Schmidt-QR (kein Auslöschungsproblem).

    Anwendung:
        LGS (überbestimmt, Methode der kleinsten Quadrate)
        Eigenwertberechnung (QR-Iteration)
        Singulärwertzerlegung

    @param matrix: Matrix A (m×n, m ≥ n)
    @return: (Q, R) – orthogonale Matrix und obere Dreiecksmatrix
    @raises ValueError: Wenn m < n (zu wenig Zeilen)
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    m = matrix.rows
    n = matrix.cols

    if m < n:
        raise ValueError(f"QR erfordert m ≥ n, erhalten {m}×{n}")

    # Numpy für numerische Effizienz verwenden
    A = np.array([[matrix._data[i][j] for j in range(n)] for i in range(m)], dtype=float)
    Q = np.eye(m)  # Startet als Einheitsmatrix, wird durch Reflexionen aufgebaut

    for k in range(n):
        # Spaltenvektor unterhalb (und auf) der Diagonale
        x = A[k:, k].copy()
        e1 = np.zeros(len(x))
        e1[0] = 1.0

        # Householder-Vektor: v = x ± ||x||·e₁ (+ für Stabilität wenn x[0] ≥ 0)
        sign = 1.0 if x[0] >= 0 else -1.0
        v = x + sign * np.linalg.norm(x) * e1

        if np.linalg.norm(v) < 1e-14:
            continue

        # Householder-Matrix H = I - 2·v·vᵀ / (vᵀ·v) (nur implizit verwendet)
        v_norm_sq = np.dot(v, v)

        # A aktualisieren: A[k:, k:] = A[k:, k:] - 2·v·(vᵀ·A[k:, k:]) / (vᵀ·v)
        A[k:, k:] -= 2 * np.outer(v, v @ A[k:, k:]) / v_norm_sq

        # Q aktualisieren: Q[:, k:] = Q[:, k:] - 2·(Q[:, k:]·v)·vᵀ / (vᵀ·v)
        Q[:, k:] -= 2 * np.outer(Q[:, k:] @ v, v) / v_norm_sq

    # In Matrix-Objekte konvertieren
    Q_matrix = Matrix(Q.tolist())
    R_matrix = Matrix(A.tolist())

    return Q_matrix, R_matrix


# =============================================================================
# SINGULÄRWERTZERLEGUNG (SVD) via numpy
# =============================================================================

def svd(matrix: 'Matrix') -> tuple['Matrix', list[float], 'Matrix']:
    """
    Singulärwertzerlegung (Singular Value Decomposition, SVD).

    Jede Matrix A (m×n) lässt sich zerlegen in:
        A = U·Σ·Vᵀ

    Wobei:
        U   – unitäre Matrix (m×m), Spalten sind linke Singulärvektoren
        Σ   – Diagonalmatrix (m×n) mit Singulärwerten σ₁ ≥ σ₂ ≥ ... ≥ 0
        Vᵀ  – unitäre Matrix (n×n), Zeilen sind rechte Singulärvektoren

    Eigenschaften:
        - Singulärwerte σᵢ sind immer ≥ 0
        - Für symmetrische positive Matrizen: σᵢ = |λᵢ|
        - Matrixrang = Anzahl σᵢ > 0
        - Pseudoinverse: A⁺ = V·Σ⁺·Uᵀ (Least-Squares-Lösung)
        - Niedrigrangapproximation (PCA, Komprimierung): Aₖ = Σᵢ₌₁ᵏ σᵢ·uᵢ·vᵢᵀ

    Anwendung:
        - Hauptkomponentenanalyse (PCA)
        - Dimensionsreduktion
        - Matrixnorm: ||A||₂ = σ₁ (größter Singulärwert)
        - Konditionszahl: κ = σ₁/σₙ

    @param matrix: Beliebige Matrix A (m×n)
    @return: (U, sigma, Vt) – linke Vektoren, Singulärwerte, rechte Vektoren transponiert
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Numpy-Array aus Matrix-Daten
    A = np.array([[matrix._data[i][j] for j in range(matrix.cols)]
                   for i in range(matrix.rows)], dtype=float)

    # SVD via numpy (verwendet LAPACK-Routinen)
    U_np, s_np, Vt_np = np.linalg.svd(A, full_matrices=True)

    # Konvertierung zurück in Matrix-Objekte und Liste
    U = Matrix(U_np.tolist())
    Vt = Matrix(Vt_np.tolist())
    sigma = list(s_np)  # Singulärwerte als Liste

    return U, sigma, Vt


def matrix_rank(matrix: 'Matrix', tol: float = 1e-10) -> int:
    """
    Berechnet den Rang einer Matrix via SVD.

    Rang = Anzahl der Singulärwerte > tol.
    Numerisch robuster als Gauss-Elimination bei fast-singulären Matrizen.

    @param matrix: Eingangsmatrix A
    @param tol: Toleranz für "Null" (Standard: 1e-10)
    @return: Rang der Matrix
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _, sigma, _ = svd(matrix)
    return sum(1 for s in sigma if s > tol)


def condition_number(matrix: 'Matrix') -> float:
    """
    Berechnet die Konditionszahl κ(A) = σ_max / σ_min.

    Die Konditionszahl misst die Empfindlichkeit des LGS Ax=b
    gegenüber Störungen in b oder A:
        ||Δx|| / ||x|| ≤ κ(A) · ||Δb|| / ||b||

    Kleine κ ≈ 1: gut konditioniert
    Große κ: schlecht konditioniert (numerische Probleme)

    @param matrix: Quadratische Matrix
    @return: Konditionszahl κ ≥ 1
    @raises ValueError: Wenn Matrix singulär (σ_min = 0)
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _, sigma, _ = svd(matrix)
    if len(sigma) == 0 or sigma[-1] < 1e-15:
        # σ_min = 0 bedeutet singuläre Matrix → Konditionszahl unendlich
        raise SingularMatrixError("Konditionszahl-Berechnung (Matrix singulär, σ_min = 0)")
    return sigma[0] / sigma[-1]


# ===========================================================================
# GIVENS-ROTATIONEN (QR-ZERLEGUNG ALTERNATIV)
# ===========================================================================

def givens_rotation_matrix(n: int, i: int, j: int, theta: float) -> 'Matrix':
    """
    @brief Erstellt eine n×n Givens-Rotationsmatrix G(i, j, θ).
    @description
        Eine Givens-Rotation rotiert einen Vektor in der (i,j)-Ebene
        um den Winkel θ. Die Matrix ist eine modifizierte Einheitsmatrix:

            G[i,i] = G[j,j] = cos(θ)
            G[i,j] = -sin(θ)
            G[j,i] =  sin(θ)
            alle anderen Einträge: Einheitsmatrix-Werte

        Anwendung: Um in einem Vektor das j-te Element auf 0 zu setzen,
        wählt man θ so, dass die Rotation den (i,j)-Block anpasst.

        Eigenschaften:
        - G^T = G^{-1} (orthogonale Matrix)
        - det(G) = 1 (eigentliche Rotation, keine Spiegelung)
        - G^T @ G = I

    @param n: Dimension der Rotationsmatrix (n×n)
    @param i: Index der ersten Achse (0-basiert)
    @param j: Index der zweiten Achse (0-basiert, j > i)
    @param theta: Rotationswinkel in Radiant
    @return: n×n Givens-Rotationsmatrix als Matrix-Objekt
    @raises ValueError: Wenn i, j ungültig oder i == j
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    if i < 0 or j < 0 or i >= n or j >= n:
        raise ValueError(f"Ungültige Indizes i={i}, j={j} für n={n}")
    if i == j:
        raise ValueError("i und j müssen verschieden sein")

    # Startwert: Einheitsmatrix
    data = [[1.0 if r == c else 0.0 for c in range(n)] for r in range(n)]

    c = math.cos(theta)  # Cosinus des Winkels
    s = math.sin(theta)  # Sinus des Winkels

    # Rotationseinträge setzen (modifiziert nur den (i,j)-Block)
    data[i][i] = c    # G[i,i] = cos(θ)
    data[j][j] = c    # G[j,j] = cos(θ)
    data[i][j] = -s   # G[i,j] = -sin(θ)
    data[j][i] = s    # G[j,i] =  sin(θ)

    return Matrix(data)


def givens_qr_decomposition(A: 'Matrix') -> tuple['Matrix', 'Matrix']:
    """
    @brief QR-Zerlegung via Givens-Rotationen (Alternative zu Householder).
    @description
        Zerlegt eine Matrix A (m×n) in:
            A = Q · R

        Wobei:
            Q – orthogonale Matrix (Q^T @ Q = I), m×m
            R – obere Dreiecksmatrix, m×n

        Algorithmus:
        - Iteriere über alle Spalten j = 0, ..., n-1
          und alle Zeilen i = m-1, ..., j+1 (von unten nach oben)
        - Für jeden Eintrag R[i,j] ungleich 0:
          * a = R[j,j], b = R[i,j]
          * r = sqrt(a² + b²)
          * cos(θ) = a/r,  sin(θ) = -b/r
          → Givens-Rotation G(j,i,θ) setzt R[i,j] auf 0
        - Q = G_1^T · G_2^T · ... (Produkt der transponierten Rotationen)

        Vorteil gegenüber Householder:
        - Gut für dünn besetzte Matrizen (spärlich)
        - Einfach parallelisierbar
        - Jede Rotation löscht genau ein Element

    @param A: Matrix (m×n, m ≥ n)
    @return: Tupel (Q, R) als Matrix-Objekte
    @raises ValueError: Wenn m < n
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    m = A.rows
    n = A.cols

    if m < n:
        raise ValueError(f"Givens-QR erfordert m ≥ n, erhalten {m}×{n}")

    # Arbeitskopie als numpy-Array für effiziente Berechnungen
    R = np.array([[A._data[i][j] for j in range(n)] for i in range(m)], dtype=float)
    # Q startet als Einheitsmatrix, wird durch transponierte Rotationen aufgebaut
    Q = np.eye(m)

    # Spalte für Spalte bearbeiten
    for j in range(n):
        # Von der untersten Zeile bis zur Zeile j+1 (von unten nach oben)
        for i in range(m - 1, j, -1):
            # Elemente auf die wir die Rotation anwenden
            a = R[j, j]  # Diagonalelement (soll erhalten bleiben)
            b = R[i, j]  # Element das auf 0 gesetzt werden soll

            # Wenn b bereits (nahe) 0 ist, keine Rotation nötig
            if abs(b) < 1e-15:
                continue

            # Rotationswinkel berechnen
            r = math.sqrt(a * a + b * b)  # Hypotenuse
            c = a / r    # cos(θ)
            s = -b / r   # sin(θ) — negativ, damit R[i,j] zu 0 wird

            # Givens-Rotation G auf die relevanten Zeilen von R anwenden
            # Nur Zeilen j und i betroffen: [j-te Zeile, i-te Zeile] * G^T
            row_j = R[j, :].copy()
            row_i = R[i, :].copy()
            R[j, :] = c * row_j - s * row_i   # Neue j-te Zeile
            R[i, :] = s * row_j + c * row_i   # Neue i-te Zeile (wird zu 0 an Pos j)

            # Q akkumulieren: Q = Q · G^T (transponierte Rotation)
            # G^T hat c an (j,j),(i,i) und s an (j,i), -s an (i,j)
            col_j = Q[:, j].copy()
            col_i = Q[:, i].copy()
            Q[:, j] = c * col_j - s * col_i
            Q[:, i] = s * col_j + c * col_i

    # Sehr kleine Werte unterhalb der Diagonalen auf 0 setzen (numerisches Rauschen)
    for i in range(m):
        for j in range(min(i, n)):
            if abs(R[i, j]) < 1e-12:
                R[i, j] = 0.0

    # In Matrix-Objekte konvertieren
    Q_matrix = Matrix(Q.tolist())
    R_matrix = Matrix(R.tolist())

    return Q_matrix, R_matrix


def givens_solve_least_squares(A: 'Matrix', b: list[float]) -> list[float]:
    """
    @brief Löst das Kleinste-Quadrate-Problem min||Ax - b||₂ via Givens-QR.
    @description
        Für überbestimmte Systeme (m > n) gibt es i.A. keine exakte Lösung.
        Gesucht ist x* = argmin ||Ax - b||₂.

        Methode: QR-Zerlegung und Rücksubstitution
        1. A = Q · R  via Givens-QR
        2. ||Ax - b||₂ = ||Q·R·x - b||₂ = ||R·x - Q^T·b||₂  (Q orthogonal)
        3. Rücksubstitution: Löse R₁·x = (Q^T·b)₁  (nur obere n×n-Teilmatrix)

        Das System R₁·x = c₁ ist eindeutig lösbar wenn rang(A) = n (voller Spaltenrang).

        Für quadratische Systeme (m = n) ist dies äquivalent zur exakten Lösung.

    @param A: Koeffizientenmatrix (m×n, m ≥ n)
    @param b: Rechte Seite als Liste der Länge m
    @return: Lösungsvektor x als Liste der Länge n
    @raises ValueError: Wenn Dimensionen nicht passen oder System singulär
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    m = A.rows
    n = A.cols

    if len(b) != m:
        raise ValueError(f"Dimension von b ({len(b)}) passt nicht zu A ({m}×{n})")
    if m < n:
        raise ValueError(f"Unterbestimmtes System: m={m} < n={n}")

    # QR-Zerlegung via Givens-Rotationen
    Q, R = givens_qr_decomposition(A)

    # c = Q^T · b berechnen
    b_np = np.array(b, dtype=float)
    Q_np = np.array([[Q._data[i][j] for j in range(Q.cols)] for i in range(Q.rows)], dtype=float)
    c = Q_np.T @ b_np  # Q^T @ b, Länge m

    # Rücksubstitution: Löse R₁ · x = c₁
    # R₁ ist das n×n Oberdreiecks-Teilblock von R
    R_np = np.array([[R._data[i][j] for j in range(R.cols)] for i in range(R.rows)], dtype=float)
    R1 = R_np[:n, :n]  # Nur die obere quadratische n×n-Teilmatrix
    c1 = c[:n]          # Entsprechende Einträge von c

    # Rücksubstitution (von unten nach oben)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if abs(R1[i, i]) < 1e-14:
            raise ValueError(f"Matrix A hat keinen vollen Spaltenrang (R[{i},{i}] ≈ 0)")
        # x[i] = (c1[i] - Σ_{j>i} R1[i,j]·x[j]) / R1[i,i]
        x[i] = (c1[i] - np.dot(R1[i, i+1:], x[i+1:])) / R1[i, i]

    return x.tolist()
