"""
@file linear_algebra.py
@brief Lineare-Algebra-Modul: Vektoren, Matrizen, Eigenwerte, Gram-Schmidt.
@description
    Implementiert fundamentale Konzepte der Linearen Algebra:

    Vector: n-dimensionaler Vektor mit Skalarprodukt, Kreuzprodukt, Norm
    Matrix: m×n-Matrix mit Multiplikation, Transponierung, Determinante,
            Inverser, Gauss-Elimination, Eigenwertberechnung

    Algorithmen:
    - Gauss-Jordan-Elimination für LGS und Inverse
    - Laplace-Entwicklung für Determinante (rekursiv)
    - QR-Iteration für Eigenwerte
    - Gram-Schmidt für Orthogonalisierung

@author Kurt Ingwer
@date 2026-03-05
"""

import math
import copy
from typing import List, Optional, Union
import numpy as np


# =============================================================================
# KLASSE: Vector
# =============================================================================

class Vector:
    """
    @brief n-dimensionaler reeller Vektor.
    @description
        Unterstützt: Skalarprodukt, Kreuzprodukt (3D), Norm,
        Normierung, Addition, Skalarmultiplikation.
    @date 2026-03-05
    """

    def __init__(self, components: list):
        """
        @brief Erstellt einen Vektor aus einer Komponentenliste.
        @param components Liste der Vektorkomponenten.
        @date 2026-03-05
        """
        self.components = list(components)
        self.dim = len(components)

    def dot(self, other: 'Vector') -> float:
        """
        @brief Berechnet das Skalarprodukt (inneres Produkt).
        @description
            Das Skalarprodukt von v und w:
                v · w = v_1*w_1 + v_2*w_2 + ... + v_n*w_n

            Geometrische Bedeutung:
                v · w = |v| * |w| * cos(θ)
            Wobei θ der Winkel zwischen den Vektoren ist.

            Eigenschaften:
            - v · w = 0 genau dann wenn v ⊥ w (senkrecht)
            - v · v = |v|^2

        @param other Der zweite Vektor.
        @return Skalarprodukt v · w.
        @raises ValueError Bei unterschiedlichen Dimensionen.
        @date 2026-03-05
        """
        if self.dim != other.dim:
            raise ValueError(f"Dimensionen passen nicht: {self.dim} vs {other.dim}")
        return sum(a * b for a, b in zip(self.components, other.components))

    def norm(self) -> float:
        """
        @brief Berechnet die euklidische Norm (Länge) des Vektors.
        @description
            Die euklidische Norm (L2-Norm):
                ||v|| = sqrt(v_1^2 + v_2^2 + ... + v_n^2) = sqrt(v · v)

        @return Euklidische Länge des Vektors.
        @date 2026-03-05
        """
        return math.sqrt(self.dot(self))

    def normalize(self) -> 'Vector':
        """
        @brief Gibt den normierten Einheitsvektor zurück.
        @description
            Einheitsvektor: v_hat = v / ||v||
            Eigenschaften: ||v_hat|| = 1, gleiche Richtung wie v.

        @return Normierter Vektor mit Länge 1.
        @raises ValueError Wenn v der Nullvektor ist.
        @date 2026-03-05
        """
        n = self.norm()
        if n < 1e-15:
            raise ValueError("Nullvektor kann nicht normiert werden")
        return Vector([c / n for c in self.components])

    def cross(self, other: 'Vector') -> 'Vector':
        """
        @brief Berechnet das Kreuzprodukt (nur für 3D-Vektoren).
        @description
            Das Kreuzprodukt v × w ist nur in 3 Dimensionen definiert:
                v × w = |i  j  k |
                        |v1 v2 v3| = (v2*w3 - v3*w2, v3*w1 - v1*w3, v1*w2 - v2*w1)
                        |w1 w2 w3|

            Eigenschaften:
            - v × w steht senkrecht auf v und w
            - |v × w| = |v| * |w| * sin(θ) = Flächeninhalt des Parallelogramms
            - v × w = -(w × v) (antikommutativ)
            - v × v = 0

        @param other Der zweite Vektor (muss ebenfalls 3D sein).
        @return Kreuzprodukt-Vektor.
        @raises ValueError Wenn Dimensionen nicht 3 sind.
        @date 2026-03-05
        """
        if self.dim != 3 or other.dim != 3:
            raise ValueError(f"Kreuzprodukt nur für 3D-Vektoren (nicht {self.dim}D)")
        a = self.components
        b = other.components
        return Vector([
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ])

    def __add__(self, other: 'Vector') -> 'Vector':
        """
        @brief Vektoraddition komponentenweise.
        @date 2026-03-05
        """
        if self.dim != other.dim:
            raise ValueError(f"Dimensionen passen nicht: {self.dim} vs {other.dim}")
        return Vector([a + b for a, b in zip(self.components, other.components)])

    def __sub__(self, other: 'Vector') -> 'Vector':
        """
        @brief Vektorsubtraktion komponentenweise.
        @date 2026-03-05
        """
        if self.dim != other.dim:
            raise ValueError(f"Dimensionen passen nicht: {self.dim} vs {other.dim}")
        return Vector([a - b for a, b in zip(self.components, other.components)])

    def __mul__(self, scalar: float) -> 'Vector':
        """
        @brief Skalarmultiplikation v * s.
        @date 2026-03-05
        """
        return Vector([c * scalar for c in self.components])

    def __rmul__(self, scalar: float) -> 'Vector':
        """@brief Skalarmultiplikation s * v (rechtsseitig). @date 2026-03-05"""
        return self.__mul__(scalar)

    def __repr__(self) -> str:
        return f"Vector({self.components})"


# =============================================================================
# KLASSE: Matrix
# =============================================================================

class Matrix:
    """
    @brief m×n-Matrix mit allen fundamentalen Operationen der Linearen Algebra.
    @description
        Intern gespeichert als Liste von Zeilen (List[List[float]]).
        Bietet: Multiplikation, Transponierung, Determinante, Inverse,
                Gauss-Elimination, LGS-Lösung, Eigenwerte.
    @date 2026-03-05
    """

    def __init__(self, data: List[List[float]]):
        """
        @brief Erstellt eine Matrix aus einer Liste von Zeilen.
        @param data Liste von Zeilen, jede Zeile eine Liste von Zahlen.
        @date 2026-03-05
        """
        if not data or not data[0]:
            raise ValueError("Matrix darf nicht leer sein")
        self.rows = len(data)
        self.cols = len(data[0])
        # Tiefe Kopie um ungewollte Seiteneffekte zu vermeiden
        self._data = [list(row) for row in data]

    def get(self, row: int, col: int) -> float:
        """@brief Gibt den Wert an Position (row, col) zurück. @date 2026-03-05"""
        return self._data[row][col]

    def set(self, row: int, col: int, value: float):
        """@brief Setzt den Wert an Position (row, col). @date 2026-03-05"""
        self._data[row][col] = value

    @classmethod
    def identity(cls, n: int) -> 'Matrix':
        """
        @brief Erstellt eine n×n-Einheitsmatrix.
        @description
            Die Einheitsmatrix I_n hat auf der Diagonale 1, sonst 0.
            Eigenschaft: A @ I = I @ A = A
        @param n Dimension.
        @return n×n-Einheitsmatrix.
        @date 2026-03-05
        """
        data = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        return cls(data)

    @classmethod
    def zeros(cls, rows: int, cols: int) -> 'Matrix':
        """@brief Erstellt eine Nullmatrix der Größe rows×cols. @date 2026-03-05"""
        return cls([[0.0] * cols for _ in range(rows)])

    def __matmul__(self, other: 'Matrix') -> 'Matrix':
        """
        @brief Matrixmultiplikation A @ B.
        @description
            Das Element (i,j) des Produkts:
                (A @ B)[i,j] = Summe_k A[i,k] * B[k,j]

            = Skalarprodukt der i-ten Zeile von A mit der j-ten Spalte von B.

            Voraussetzung: A ist m×n, B muss n×p sein -> Ergebnis ist m×p.

        @param other Die zweite Matrix.
        @return Produktmatrix.
        @raises ValueError Wenn Dimensionen nicht kompatibel.
        @date 2026-03-05
        """
        if self.cols != other.rows:
            raise ValueError(
                f"Matrixdimensionen inkompatibel: {self.rows}×{self.cols} @ {other.rows}×{other.cols}"
            )
        result = Matrix.zeros(self.rows, other.cols)
        for i in range(self.rows):
            for j in range(other.cols):
                val = sum(self.get(i, k) * other.get(k, j) for k in range(self.cols))
                result.set(i, j, val)
        return result

    def transpose(self) -> 'Matrix':
        """
        @brief Transponiert die Matrix (Zeilen und Spalten tauschen).
        @description
            (A^T)[i,j] = A[j,i]

            Eigenschaften:
            - (A^T)^T = A
            - (A@B)^T = B^T @ A^T
            - Für symmetrische Matrizen: A^T = A

        @return Transponierte Matrix.
        @date 2026-03-05
        """
        data = [[self.get(i, j) for i in range(self.rows)] for j in range(self.cols)]
        return Matrix(data)

    def trace(self) -> float:
        """
        @brief Berechnet die Spur (Summe der Hauptdiagonalelemente).
        @description
            Spur(A) = A[0,0] + A[1,1] + ... + A[n,n]

            Eigenschaften:
            - Spur(A+B) = Spur(A) + Spur(B)
            - Spur(AB) = Spur(BA)
            - Spur = Summe der Eigenwerte

        @return Spur der Matrix.
        @raises ValueError Wenn Matrix nicht quadratisch.
        @date 2026-03-05
        """
        if self.rows != self.cols:
            raise ValueError("Spur nur für quadratische Matrizen definiert")
        return sum(self.get(i, i) for i in range(self.rows))

    def determinant(self) -> float:
        """
        @brief Berechnet die Determinante via Gauss-Elimination.
        @description
            Die Determinante misst den "Volumenskalierungsfaktor" einer Transformation.

            Eigenschaften:
            - det(A) = 0 genau dann wenn A singulär (nicht invertierbar)
            - det(A @ B) = det(A) * det(B)
            - det(A^T) = det(A)
            - Zeilenvertauschung: Vorzeichen wechselt

            Berechnung: Gauss-Elimination zur oberen Dreiecksform,
            dann ist det = Produkt der Diagonalelemente (mit Vorzeichen).

        @return Determinante der Matrix.
        @raises ValueError Wenn Matrix nicht quadratisch.
        @date 2026-03-05
        """
        if self.rows != self.cols:
            raise ValueError("Determinante nur für quadratische Matrizen definiert")

        n = self.rows
        # Arbeitskopie erstellen
        m = [list(row) for row in self._data]
        sign = 1  # Vorzeichen (wechselt bei Zeilenvertauschung)

        for col in range(n):
            # Pivotelement suchen (größter Betrag in Spalte -> numerische Stabilität)
            max_row = max(range(col, n), key=lambda r: abs(m[r][col]))

            if abs(m[max_row][col]) < 1e-15:
                return 0.0  # Singuläre Matrix

            # Zeilen tauschen (falls nötig)
            if max_row != col:
                m[col], m[max_row] = m[max_row], m[col]
                sign *= -1  # Vorzeichen wechseln

            # Elimination: alle Zeilen unter dem Pivot nullen
            pivot = m[col][col]
            for row in range(col + 1, n):
                factor = m[row][col] / pivot
                for k in range(col, n):
                    m[row][k] -= factor * m[col][k]

        # Determinante = Produkt der Pivots (Diagonale der Dreiecksmatrix)
        det = sign
        for i in range(n):
            det *= m[i][i]

        return det

    def inverse(self) -> 'Matrix':
        """
        @brief Berechnet die Inverse via Gauss-Jordan-Elimination.
        @description
            Die inverse Matrix A^{-1} erfüllt: A @ A^{-1} = A^{-1} @ A = I

            Algorithmus (Gauss-Jordan):
            Erweiterte Matrix [A | I] durch Zeilenoperationen in [I | A^{-1}] umformen.

            Existiert genau dann wenn det(A) ≠ 0.

        @return Inverse Matrix A^{-1}.
        @raises ValueError Wenn Matrix singulär oder nicht quadratisch.
        @date 2026-03-05
        """
        if self.rows != self.cols:
            raise ValueError("Inverse nur für quadratische Matrizen definiert")

        n = self.rows
        # Erweiterte Matrix [A | I] erstellen
        aug = [list(self._data[i]) + [1.0 if i == j else 0.0 for j in range(n)]
               for i in range(n)]

        # Vorwärts-Elimination mit Pivotsuche
        for col in range(n):
            # Pivotzeile suchen
            max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
            if abs(aug[max_row][col]) < 1e-15:
                raise ValueError("Matrix ist singulär (nicht invertierbar)")
            aug[col], aug[max_row] = aug[max_row], aug[col]

            # Pivot normieren (Pivot-Element auf 1 setzen)
            pivot = aug[col][col]
            aug[col] = [v / pivot for v in aug[col]]

            # Alle anderen Zeilen eliminieren (Gauss-Jordan: auch oberhalb!)
            for row in range(n):
                if row != col:
                    factor = aug[row][col]
                    aug[row] = [aug[row][k] - factor * aug[col][k] for k in range(2 * n)]

        # Rechte Seite [I | A^{-1}] extrahieren
        inv_data = [row[n:] for row in aug]
        return Matrix(inv_data)

    def solve(self, b: 'Vector') -> 'Vector':
        """
        @brief Löst das lineare Gleichungssystem Ax = b.
        @description
            Lösung via Gauss-Elimination mit Rücksubstitution.
            Äquivalent zu x = A^{-1} @ b, aber numerisch stabiler.

        @param b Rechte Seite (Vektor).
        @return Lösungsvektor x.
        @raises ValueError Wenn System keine eindeutige Lösung hat.
        @date 2026-03-05
        """
        if self.rows != self.cols:
            raise ValueError("Gauss-Elimination nur für quadratische Matrizen")
        if self.rows != b.dim:
            raise ValueError(f"Dimension passt nicht: Matrix {self.rows}×{self.cols}, Vektor {b.dim}")

        n = self.rows
        # Erweiterte Matrix [A | b]
        aug = [list(self._data[i]) + [b.components[i]] for i in range(n)]

        # Vorwärts-Elimination mit Pivotsuche
        for col in range(n):
            max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
            if abs(aug[max_row][col]) < 1e-15:
                raise ValueError("Matrix ist singulär, LGS hat keine eindeutige Lösung")
            aug[col], aug[max_row] = aug[max_row], aug[col]

            pivot = aug[col][col]
            for row in range(col + 1, n):
                factor = aug[row][col] / pivot
                for k in range(col, n + 1):
                    aug[row][k] -= factor * aug[col][k]

        # Rücksubstitution
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            x[i] = aug[i][n]
            for j in range(i + 1, n):
                x[i] -= aug[i][j] * x[j]
            x[i] /= aug[i][i]

        return Vector(x)

    def eigenvalues(self) -> list:
        """
        @brief Berechnet die Eigenwerte via QR-Iteration (Power-Methode für 2x2).
        @description
            Ein Eigenwert λ und zugehöriger Eigenvektor v erfüllen:
                A @ v = λ * v

            Für 2x2-Matrizen: Charakteristisches Polynom lösen:
                det(A - λI) = 0
                λ^2 - Spur(A)*λ + det(A) = 0

            Die Eigenwerte sind die Wurzeln dieses Polynoms.
            Für größere Matrizen: iterative QR-Algorithmus.

        @return Liste der Eigenwerte (komplex möglich).
        @raises ValueError Wenn Matrix nicht quadratisch.
        @date 2026-03-05
        """
        if self.rows != self.cols:
            raise ValueError("Eigenwerte nur für quadratische Matrizen definiert")

        n = self.rows

        if n == 2:
            # Für 2x2: analytische Formel über das charakteristische Polynom
            # λ^2 - Spur*λ + det = 0
            trace = self.trace()
            det = self.determinant()
            # pq-Formel: λ = (Spur ± sqrt(Spur^2 - 4*det)) / 2
            discriminant = trace**2 - 4 * det
            import cmath
            sqrt_d = cmath.sqrt(discriminant)
            return [(trace + sqrt_d) / 2, (trace - sqrt_d) / 2]

        # Für größere Matrizen: QR-Iteration
        return self._qr_eigenvalues()

    def _qr_eigenvalues(self, max_iter: int = 1000, tol: float = 1e-10) -> list:
        """
        @brief QR-Iteration zur Eigenwertberechnung für n×n-Matrizen.
        @description
            Der QR-Algorithmus konvergiert die Matrix A iterativ zur
            oberen Dreiecksform (Schur-Form), deren Diagonale die Eigenwerte enthält:

            Schritt k:
                A_k = Q_k @ R_k  (QR-Zerlegung)
                A_{k+1} = R_k @ Q_k

            Konvergenz: Die Matrix A_k konvergiert gegen obere Dreiecksform.

        @return Näherungswerte der Eigenwerte.
        @date 2026-03-05
        """
        import numpy as np
        # numpy für numerische Stabilität verwenden
        A = np.array(self._data, dtype=float)
        for _ in range(max_iter):
            Q, R = np.linalg.qr(A)
            A_new = R @ Q
            # Konvergenztest: Unterdiagonale nahe 0?
            off_diag = np.sum(np.abs(np.tril(A_new, -1)))
            A = A_new
            if off_diag < tol:
                break
        return list(np.diag(A))

    def eigenvectors(self) -> list:
        """
        @brief Berechnet Eigenwert-Eigenvektor-Paare der Matrix.
        @description
            Gibt für jeden Eigenwert λ den zugehörigen normierten Eigenvektor v zurück.
            Es gilt: A @ v = λ * v

            Methode:
            - Eigenwerte via eigenvalues() berechnen
            - Für jeden Eigenwert: (A - λI)v = 0 lösen
              → Kern von (A - λI) via SVD bestimmen (letzter rechter Singulärvektor
                zur kleinsten Singulärwert = Kern-Basis)

            Vorteil SVD: Numerisch stabil auch bei fast-singulären Matrizen
            (A - λI ist per Definition singulär → Gauss-Elimination instabil).

        @return Liste von Tupeln (eigenwert, eigenvektor_als_Vector)
                Eigenvektoren sind normiert (||v|| = 1).
        @raises ValueError Wenn Matrix nicht quadratisch.
        @lastModified: 2026-03-08
        """
        if self.rows != self.cols:
            raise ValueError("Eigenvektoren nur für quadratische Matrizen definiert")

        eigenvals = self.eigenvalues()
        result = []
        A_np = np.array(self._data, dtype=float)

        for lam in eigenvals:
            # (A - λI) bilden: B = A - λ·I
            lam_real = lam.real if hasattr(lam, 'real') else float(lam)
            B = A_np - lam_real * np.eye(self.rows)

            # SVD von B: B = U·Σ·Vᵀ
            # Eigenvektor = letzte Spalte von V (kleinster Singulärwert → Kern)
            _, _, Vt = np.linalg.svd(B)
            eigenvec_arr = Vt[-1]  # Letzte Zeile von Vᵀ = letzter rechter Singulärvektor

            # Normieren (SVD liefert bereits normierte Vektoren, aber sicherheitshalber)
            norm = np.linalg.norm(eigenvec_arr)
            if norm > 1e-14:
                eigenvec_arr = eigenvec_arr / norm

            # In internen Vector-Typ konvertieren
            eigenvec = Vector(eigenvec_arr.tolist())
            result.append((complex(lam_real), eigenvec))

        return result

    def __repr__(self) -> str:
        rows_str = '\n  '.join(str(row) for row in self._data)
        return f"Matrix([\n  {rows_str}\n])"


# =============================================================================
# GRAM-SCHMIDT-ORTHOGONALISIERUNG
# =============================================================================

def gram_schmidt(vectors: List[Vector], normalize: bool = False) -> List[Vector]:
    """
    @brief Gram-Schmidt-Orthogonalisierung einer Menge von Vektoren.
    @description
        Das Gram-Schmidt-Verfahren erzeugt aus einer Basis {v_1, ..., v_k}
        eine orthogonale (oder orthonormale) Basis {u_1, ..., u_k}:

            u_1 = v_1
            u_k = v_k - Summe_{j=1}^{k-1} proj_{u_j}(v_k)

        Dabei ist die orthogonale Projektion:
            proj_u(v) = (v · u / u · u) * u

        Geometrische Idee: Jeder neue Vektor wird so modifiziert,
        dass er senkrecht auf allen vorherigen steht (Subtraktion der Projektionen).

        Bei normalize=True: Zusätzlich jeden Vektor auf Länge 1 normieren
        -> Orthonormalbasis (ONB).

    @param vectors Liste der zu orthogonalisierenden Vektoren.
    @param normalize True für ONB (Orthonormalbasis), False für OB.
    @return Liste orthogonaler (oder orthonormaler) Vektoren.
    @date 2026-03-05
    """
    orthogonal = []

    for v in vectors:
        # Startvektor für diesen Schritt: v_k
        u = Vector(list(v.components))

        # Alle bereits berechneten orthogonalen Vektoren abziehen
        for prev in orthogonal:
            # Projektion von u auf prev berechnen: proj = (u·prev / prev·prev) * prev
            proj_scalar = u.dot(prev) / prev.dot(prev)
            projection = prev * proj_scalar
            # Projektion abziehen -> u wird senkrecht zu prev
            u = u - projection

        # Null-Vektor überspringen (linear abhängige Vektoren)
        if u.norm() > 1e-10:
            if normalize:
                u = u.normalize()
            orthogonal.append(u)

    return orthogonal


# =============================================================================
# LU-ZERLEGUNG (Doolittle mit Teilpivotisierung)
# =============================================================================

def lu_decomposition(matrix: 'Matrix') -> tuple:
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
    @lastModified: 2026-03-08
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
            raise ValueError("Matrix ist singulär (LU-Zerlegung nicht möglich)")

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

def qr_decomposition(matrix: 'Matrix') -> tuple:
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
    @lastModified: 2026-03-08
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

def svd(matrix: 'Matrix') -> tuple:
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
    """
    _, sigma, _ = svd(matrix)
    if len(sigma) == 0 or sigma[-1] < 1e-15:
        raise ValueError("Matrix ist singulär (Konditionszahl = ∞)")
    return sigma[0] / sigma[-1]
