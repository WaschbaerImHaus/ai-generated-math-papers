"""
@file matrix_ops.py
@brief Matrix-Modul: m×n-Matrizen mit allen fundamentalen Operationen.
@description
    Enthält die Matrix-Klasse mit:
    - Grundoperationen: Multiplikation, Transponierung, Spur
    - Gauss-Jordan-Elimination für Inverse und LGS
    - Determinante via Gauss-Elimination
    - Eigenwerte (analytisch für 2×2, QR-Iteration für größere Matrizen)
    - Eigenvektoren via SVD-Kern

    Importiert Vector aus vectors.py.
    Ausgelagert aus linear_algebra.py für bessere Modularität.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import copy
from typing import List, Optional, Union
import numpy as np

# Vector-Klasse aus dem vectors-Submodul importieren
from vectors import Vector

# Spezifische mathematische Ausnahmen importieren
from exceptions import SingularMatrixError

# Zentrales Logging-System importieren
from math_logger import get_logger

# Modul-Logger für matrix_ops.py erstellen
_logger = get_logger("matrix_ops")

# Schwellenwert für Warnung bei schlecht konditionierten Matrizen.
# Konditionszahl > 1e10 bedeutet: Ergebnisfehler können um Faktor 10^10 verstärkt werden.
# Bei IEEE 754 double precision (eps ≈ 2.2e-16) verliert man dann ~10 signifikante Stellen.
_CONDITION_WARNING_THRESHOLD: float = 1e10


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
    @author Kurt Ingwer
    @lastModified 2026-03-10
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

    def _check_condition(self, operation: str = "Operation") -> float:
        """
        @brief Berechnet und loggt die Konditionszahl der Matrix.
        @description
            Die Konditionszahl κ(A) = ||A|| · ||A⁻¹|| misst, wie sensitiv das
            Ergebnis von Ax=b gegenüber kleinen Perturbationen in A oder b ist.

            Interpretation:
            - κ ≈ 1: gut konditioniert (numerisch stabil)
            - κ ≈ 1e8: 8 Stellen Genauigkeit verloren
            - κ > 1e10: schlecht konditioniert → Ergebnis unzuverlässig

            Berechnung über numpy.linalg.cond() (2-Norm, d.h. σ_max / σ_min der SVD).

        @param operation: Name der aufrufenden Operation (für Log-Ausgabe).
        @return Konditionszahl (float).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        # Nur für quadratische Matrizen sinnvoll
        if self.rows != self.cols:
            return 1.0

        # numpy-Array aus internen Daten aufbauen
        arr = np.array(self._data, dtype=float)
        try:
            cond = float(np.linalg.cond(arr))
        except Exception:
            # Bei numerischen Fehlern (z.B. NaN/Inf) keine Warnung ausgeben
            return 1.0

        # Warnung bei schlecht konditionierter Matrix ausgeben
        if cond > _CONDITION_WARNING_THRESHOLD:
            _logger._logger.warning(
                f"{operation}: Schlecht konditionierte Matrix! "
                f"Konditionszahl κ = {cond:.3e} > {_CONDITION_WARNING_THRESHOLD:.0e}. "
                f"Numerische Ergebnisse können stark fehlerbehaftet sein."
            )
        else:
            # Auf DEBUG-Level nur informationshalber ausgeben
            _logger._logger.debug(
                f"{operation}: Konditionszahl κ = {cond:.3e}"
            )
        return cond

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

        # Konditionszahl prüfen und bei schlecht konditionierter Matrix warnen
        self._check_condition("determinant()")

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

        # Konditionszahl prüfen: schlecht konditionierte Matrix → ungenaue Inverse
        self._check_condition("inverse()")

        n = self.rows
        # Erweiterte Matrix [A | I] erstellen
        aug = [list(self._data[i]) + [1.0 if i == j else 0.0 for j in range(n)]
               for i in range(n)]

        # Vorwärts-Elimination mit Pivotsuche
        for col in range(n):
            # Pivotzeile suchen
            max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
            if abs(aug[max_row][col]) < 1e-15:
                # Pivot nahe 0 bedeutet singuläre Matrix → keine Inverse berechenbar
                raise SingularMatrixError("Gauss-Jordan Inversion")
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

        # Konditionszahl prüfen: schlecht konditioniertes System → ungenaue Lösung
        self._check_condition("solve()")

        n = self.rows
        # Erweiterte Matrix [A | b]
        aug = [list(self._data[i]) + [b.components[i]] for i in range(n)]

        # Vorwärts-Elimination mit Pivotsuche
        for col in range(n):
            max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
            if abs(aug[max_row][col]) < 1e-15:
                # Singuläre Matrix → LGS hat keine eindeutige Lösung
                raise SingularMatrixError("LGS-Lösung (Gauss-Elimination)")
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
        @author Kurt Ingwer
        @lastModified 2026-03-08
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
