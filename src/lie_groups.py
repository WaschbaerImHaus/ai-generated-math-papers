"""
@file lie_groups.py
@brief Lie-Gruppen und Lie-Algebren: Matrixgruppen SO(n), SU(n), GL(n),
       Lie-Klammer, Exponentialabbildung, BCH-Formel, Darstellungstheorie,
       Klassifikation einfacher Lie-Algebren (Cartan-Klassifikation).
@description
    Dieses Modul implementiert die grundlegenden Strukturen der Theorie
    der Lie-Gruppen und Lie-Algebren ohne externe Lie-Bibliotheken:

    Klassen:
    - LieGroup               – Abstrakte Basisklasse für Lie-Gruppen
    - MatrixLieGroup         – Allgemeine Matrixgruppe (Untergruppe von GL(n))
    - SO                     – Spezielle orthogonale Gruppe SO(n)
    - SU                     – Spezielle unitäre Gruppe SU(n)
    - GL                     – Allgemeine lineare Gruppe GL(n, R/C)
    - LieAlgebra             – Lie-Algebra als Tangentialraum am Einselement
    - ExponentialMap         – Exponentialabbildung zwischen Gruppe und Algebra

    Freie Funktionen:
    - one_parameter_subgroup()          – Einparametrige Untergruppe exp(tX)
    - classify_simple_lie_algebra()     – Klassifikation nach Cartan (A_n..D_n)
    - exceptional_lie_algebras_info()   – Info zu G2, F4, E6, E7, E8
    - cartan_matrix()                   – Cartan-Matrix für klassische Typen
    - dynkin_diagram_info()             – Dynkin-Diagramm Beschreibung
    - fundamental_representation_su2()  – Spin-j Darstellung von SU(2)
    - casimir_element()                 – Casimir-Operator (quadratisch)

    Mathematische Grundlagen:
    - Lie-Gruppe G: differenzierbare Mannigfaltigkeit + Gruppenstruktur
    - Lie-Algebra g = T_e G: Tangentialraum am Einselement
    - Exponentialabbildung: exp: g → G, X ↦ e^X = Σ X^k/k!
    - Lie-Klammer: [X, Y] = XY - YX (Kommutator)
    - BCH-Formel: log(e^X e^Y) = X + Y + ½[X,Y] + 1/12([X,[X,Y]]+[Y,[Y,X]]) + …
    - Cartan-Klassifikation: A_n = su(n+1), B_n = so(2n+1),
                              C_n = sp(2n),   D_n = so(2n)

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import numpy as np
from scipy.linalg import expm, logm
from abc import ABC, abstractmethod
from typing import Optional

import sys
import os

# Pfad für eigene Module setzen
sys.path.insert(0, os.path.dirname(__file__))


# =============================================================================
# ABSTRAKTE BASISKLASSE: LieGroup
# =============================================================================

class LieGroup(ABC):
    """
    Abstrakte Basisklasse für Lie-Gruppen.

    Eine Lie-Gruppe ist eine Gruppe G, die gleichzeitig eine
    differenzierbare Mannigfaltigkeit ist, so dass die Gruppenoperationen
    (Multiplikation und Inversenbildung) glatte Abbildungen sind.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @abstractmethod
    def dimension(self) -> int:
        """
        Gibt die Dimension der Lie-Gruppe als Mannigfaltigkeit zurück.

        Die Dimension entspricht der Dimension der zugehörigen Lie-Algebra
        (Tangentialraum am Einselement).

        @return Dimension der Lie-Gruppe (int)
        @lastModified 2026-03-10
        """
        pass

    @abstractmethod
    def is_compact(self) -> bool:
        """
        Gibt an, ob die Lie-Gruppe kompakt ist.

        Eine kompakte Lie-Gruppe ist als topologischer Raum kompakt
        (beschränkt und abgeschlossen im übertragenen Sinn).
        Beispiele: SO(n), SU(n) sind kompakt; GL(n,R) ist nicht kompakt.

        @return True wenn kompakt, False sonst
        @lastModified 2026-03-10
        """
        pass

    @abstractmethod
    def is_connected(self) -> bool:
        """
        Gibt an, ob die Lie-Gruppe zusammenhängend ist.

        Eine zusammenhängende Lie-Gruppe lässt sich nicht in zwei
        disjunkte offene Mengen zerlegen.

        @return True wenn zusammenhängend, False sonst
        @lastModified 2026-03-10
        """
        pass


# =============================================================================
# KLASSE: MatrixLieGroup
# =============================================================================

class MatrixLieGroup(LieGroup):
    """
    Allgemeine Matrixgruppe: Untergruppe von GL(n, F).

    Eine Matrixlie-Gruppe ist eine abgeschlossene Untergruppe der
    allgemeinen linearen Gruppe GL(n, F), wobei F = ℝ oder ℂ.
    Die Gruppenoperation ist die Matrixmultiplikation.

    @param n      Größe der n×n-Matrizen
    @param field  Körper: 'real' (ℝ) oder 'complex' (ℂ)

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, n: int, field: str = 'real'):
        """
        Initialisiert eine n×n Matrixlie-Gruppe über dem angegebenen Körper.

        @param n      Größe der quadratischen Matrizen
        @param field  'real' für ℝ oder 'complex' für ℂ
        @lastModified 2026-03-10
        """
        if n < 1:
            raise ValueError("Matrixgröße n muss mindestens 1 sein.")
        if field not in ('real', 'complex'):
            raise ValueError("field muss 'real' oder 'complex' sein.")

        self.n = n          # Matrixgröße
        self.field = field  # Körper (real/complex)

    def dimension(self) -> int:
        """
        Standardimplementierung: Dimension = n² (für GL(n)).
        Unterklassen überschreiben diese Methode.

        @return n² (Dimension von GL(n) als Mannigfaltigkeit)
        @lastModified 2026-03-10
        """
        return self.n * self.n

    def is_compact(self) -> bool:
        """
        GL(n) ist nicht kompakt (Determinante kann beliebig groß werden).

        @return False
        @lastModified 2026-03-10
        """
        return False

    def is_connected(self) -> bool:
        """
        GL(n, ℝ) hat zwei Zusammenhangskomponenten (det > 0 und det < 0).
        GL(n, ℂ) ist zusammenhängend.

        @return True für ℂ, False für ℝ
        @lastModified 2026-03-10
        """
        return self.field == 'complex'

    def identity(self) -> np.ndarray:
        """
        Gibt die n×n Einheitsmatrix zurück (neutrales Element der Gruppe).

        Die Einheitsmatrix I_n ist das neutrale Element der Matrixmultiplikation:
            I_n · A = A · I_n = A für alle A ∈ G.

        @return n×n Einheitsmatrix als numpy-Array
        @lastModified 2026-03-10
        """
        if self.field == 'complex':
            return np.eye(self.n, dtype=complex)
        return np.eye(self.n, dtype=float)

    def compose(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """
        Gruppenoperation: Matrixmultiplikation A · B.

        @param A  erste Matrix (n×n)
        @param B  zweite Matrix (n×n)
        @return   Produkt A·B als numpy-Array
        @lastModified 2026-03-10
        """
        return A @ B

    def inverse(self, A: np.ndarray) -> np.ndarray:
        """
        Berechnet das Gruppeninverse A⁻¹ via numpy.linalg.inv.

        Falls A singulär ist, wird eine LinAlgError-Exception geworfen.

        @param A  invertierbare Matrix (n×n)
        @return   A⁻¹ als numpy-Array
        @raises   numpy.linalg.LinAlgError bei singulärer Matrix
        @lastModified 2026-03-10
        """
        return np.linalg.inv(A)

    def is_element(self, A: np.ndarray, tol: float = 1e-10) -> bool:
        """
        Prüft, ob die Matrix A in der Gruppe liegt.

        Für die Basisklasse MatrixLieGroup: A ∈ GL(n) ⟺ det(A) ≠ 0.
        Unterklassen überschreiben diese Methode mit spezielleren Bedingungen.

        @param A    zu prüfende Matrix (n×n)
        @param tol  Toleranz für numerische Vergleiche
        @return     True wenn A in der Gruppe liegt
        @lastModified 2026-03-10
        """
        if A.shape != (self.n, self.n):
            return False
        return bool(abs(np.linalg.det(A)) > tol)


# =============================================================================
# KLASSE: SO - Spezielle orthogonale Gruppe
# =============================================================================

class SO(MatrixLieGroup):
    """
    Spezielle orthogonale Gruppe SO(n).

    SO(n) = { A ∈ M(n×n, ℝ) | AᵀA = I_n und det(A) = 1 }

    SO(n) beschreibt alle Rotationen im ℝⁿ ohne Spiegelungen.
    Eigenschaften:
    - dim(SO(n)) = n(n-1)/2
    - SO(n) ist kompakt und zusammenhängend
    - SO(2) ≅ S¹ (Einheitskreis)
    - SO(3) beschreibt räumliche Rotationen (Euler-Winkel, Quaternionen)
    - Lie-Algebra: so(n) = schiefsymmetrische Matrizen

    @param n  Dimension des Raums

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, n: int):
        """
        Initialisiert SO(n).

        @param n  Größe der Rotationsmatrizen (n ≥ 2)
        @lastModified 2026-03-10
        """
        if n < 2:
            raise ValueError("SO(n) ist nur für n ≥ 2 definiert.")
        super().__init__(n, field='real')

    def dimension(self) -> int:
        """
        Dimension von SO(n) = n(n-1)/2.

        Dies entspricht der Anzahl der unabhängigen Rotationsparameter:
        - SO(2): dim = 1 (ein Winkel)
        - SO(3): dim = 3 (Euler-Winkel oder Achse+Winkel)
        - SO(4): dim = 6

        @return n(n-1)/2
        @lastModified 2026-03-10
        """
        return self.n * (self.n - 1) // 2

    def is_compact(self) -> bool:
        """
        SO(n) ist kompakt (als abgeschlossene beschränkte Untergruppe von GL(n)).

        @return True
        @lastModified 2026-03-10
        """
        return True

    def is_connected(self) -> bool:
        """
        SO(n) ist zusammenhängend für alle n ≥ 1.

        (O(n) hat zwei Komponenten: det = +1 und det = -1;
         SO(n) ist die Komponente mit det = +1.)

        @return True
        @lastModified 2026-03-10
        """
        return True

    def is_element(self, A: np.ndarray, tol: float = 1e-10) -> bool:
        """
        Prüft, ob A ∈ SO(n).

        Bedingungen:
        1. AᵀA = I_n  (orthogonal)
        2. det(A) ≈ +1  (spezielle orthogonale, keine Spiegelung)

        @param A    zu prüfende Matrix (n×n, reell)
        @param tol  Toleranz für numerische Abweichungen
        @return     True wenn A ∈ SO(n)
        @lastModified 2026-03-10
        """
        if A.shape != (self.n, self.n):
            return False
        # Prüfe Orthogonalität: AᵀA ≈ I
        ata = A.T @ A
        if not np.allclose(ata, np.eye(self.n), atol=tol):
            return False
        # Prüfe det(A) ≈ +1
        return bool(abs(np.linalg.det(A) - 1.0) < tol)

    def rotation_2d(self, theta: float) -> np.ndarray:
        """
        Erzeugt eine 2D-Rotationsmatrix für SO(2).

        R(θ) = [[cos θ, -sin θ],
                [sin θ,  cos θ]]

        @param theta  Rotationswinkel in Radiant
        @return       2×2 Rotationsmatrix
        @raises       ValueError wenn n ≠ 2
        @lastModified 2026-03-10
        """
        if self.n != 2:
            raise ValueError("rotation_2d() ist nur für SO(2) (n=2) definiert.")
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s],
                         [s,  c]])

    def rotation_3d_x(self, theta: float) -> np.ndarray:
        """
        3D-Rotation um die X-Achse für SO(3).

        R_x(θ) = [[1,    0,     0  ],
                  [0, cos θ, -sin θ],
                  [0, sin θ,  cos θ]]

        @param theta  Rotationswinkel in Radiant
        @return       3×3 Rotationsmatrix um X-Achse
        @raises       ValueError wenn n ≠ 3
        @lastModified 2026-03-10
        """
        if self.n != 3:
            raise ValueError("rotation_3d_x() ist nur für SO(3) (n=3) definiert.")
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[1, 0,  0],
                         [0, c, -s],
                         [0, s,  c]])

    def rotation_3d_y(self, theta: float) -> np.ndarray:
        """
        3D-Rotation um die Y-Achse für SO(3).

        R_y(θ) = [[ cos θ, 0, sin θ],
                  [   0,   1,    0  ],
                  [-sin θ, 0, cos θ]]

        @param theta  Rotationswinkel in Radiant
        @return       3×3 Rotationsmatrix um Y-Achse
        @raises       ValueError wenn n ≠ 3
        @lastModified 2026-03-10
        """
        if self.n != 3:
            raise ValueError("rotation_3d_y() ist nur für SO(3) (n=3) definiert.")
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[ c, 0, s],
                         [ 0, 1, 0],
                         [-s, 0, c]])

    def rotation_3d_z(self, theta: float) -> np.ndarray:
        """
        3D-Rotation um die Z-Achse für SO(3).

        R_z(θ) = [[cos θ, -sin θ, 0],
                  [sin θ,  cos θ, 0],
                  [  0,      0,   1]]

        @param theta  Rotationswinkel in Radiant
        @return       3×3 Rotationsmatrix um Z-Achse
        @raises       ValueError wenn n ≠ 3
        @lastModified 2026-03-10
        """
        if self.n != 3:
            raise ValueError("rotation_3d_z() ist nur für SO(3) (n=3) definiert.")
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0],
                         [s,  c, 0],
                         [0,  0, 1]])


# =============================================================================
# KLASSE: SU - Spezielle unitäre Gruppe
# =============================================================================

class SU(MatrixLieGroup):
    """
    Spezielle unitäre Gruppe SU(n).

    SU(n) = { A ∈ M(n×n, ℂ) | A†A = I_n und det(A) = 1 }

    wobei A† = ĀᵀP (konjugiert transponiert = hermitesch adjungiert).

    Eigenschaften:
    - dim(SU(n)) = n²-1
    - SU(n) ist kompakt und zusammenhängend
    - SU(2) ≅ S³ (3-Sphäre), ist die universelle Überlagerung von SO(3)
    - SU(3) ist die Eichgruppe der starken Kernkraft (QCD)
    - Lie-Algebra su(n): schief-hermitesche spurlose Matrizen

    @param n  Dimension der unitären Gruppe

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, n: int):
        """
        Initialisiert SU(n).

        @param n  Größe der unitären Matrizen (n ≥ 1)
        @lastModified 2026-03-10
        """
        if n < 1:
            raise ValueError("SU(n) ist nur für n ≥ 1 definiert.")
        super().__init__(n, field='complex')

    def dimension(self) -> int:
        """
        Dimension von SU(n) = n²-1.

        - SU(2): dim = 3 (Pauli-Matrizen σ₁, σ₂, σ₃)
        - SU(3): dim = 8 (Gell-Mann-Matrizen λ₁..λ₈)

        @return n²-1
        @lastModified 2026-03-10
        """
        return self.n * self.n - 1

    def is_compact(self) -> bool:
        """
        SU(n) ist kompakt.

        @return True
        @lastModified 2026-03-10
        """
        return True

    def is_connected(self) -> bool:
        """
        SU(n) ist zusammenhängend für alle n ≥ 1.

        @return True
        @lastModified 2026-03-10
        """
        return True

    def is_element(self, A: np.ndarray, tol: float = 1e-10) -> bool:
        """
        Prüft, ob A ∈ SU(n).

        Bedingungen:
        1. A†A = I_n  (unitär): konjugiert transponiert × A = Einheitsmatrix
        2. det(A) ≈ +1  (speziell unitär, Determinante = 1)

        @param A    zu prüfende komplexe Matrix (n×n)
        @param tol  Toleranz für numerische Abweichungen
        @return     True wenn A ∈ SU(n)
        @lastModified 2026-03-10
        """
        if A.shape != (self.n, self.n):
            return False
        # Hermitesch adjungiert: A† = Ā ᵀ
        A_dag = A.conj().T
        # Prüfe Unitarität: A†A ≈ I
        if not np.allclose(A_dag @ A, np.eye(self.n, dtype=complex), atol=tol):
            return False
        # Prüfe det(A) ≈ +1
        return bool(abs(np.linalg.det(A) - 1.0) < tol)

    def su2_from_angles(self, theta: float, phi: float, psi: float) -> np.ndarray:
        """
        Erzeugt ein SU(2)-Element aus Euler-Winkeln (θ, φ, ψ).

        Parametrisierung von SU(2) via Euler-Winkel:
            U(θ, φ, ψ) = [[cos(θ/2)·e^{i(φ+ψ)/2},  i·sin(θ/2)·e^{i(φ-ψ)/2}],
                           [i·sin(θ/2)·e^{-i(φ-ψ)/2}, cos(θ/2)·e^{-i(φ+ψ)/2}]]

        Dabei gilt: θ ∈ [0, π], φ ∈ [0, 2π), ψ ∈ [0, 4π).

        @param theta  Polar-Winkel θ ∈ [0, π]
        @param phi    Azimut-Winkel φ ∈ [0, 2π)
        @param psi    Dreh-Winkel ψ ∈ [0, 4π)
        @return       2×2 Matrix in SU(2)
        @raises       ValueError wenn n ≠ 2
        @lastModified 2026-03-10
        """
        if self.n != 2:
            raise ValueError("su2_from_angles() ist nur für SU(2) (n=2) definiert.")
        # Berechne die Matrixelemente
        half_theta = theta / 2.0
        cos_t = np.cos(half_theta)
        sin_t = np.sin(half_theta)
        # Euler-Winkel-Phasen
        phase_pp = np.exp(1j * (phi + psi) / 2.0)   # e^{i(φ+ψ)/2}
        phase_pm = np.exp(1j * (phi - psi) / 2.0)   # e^{i(φ-ψ)/2}
        # SU(2)-Matrix
        U = np.array([
            [ cos_t * phase_pp,       1j * sin_t * phase_pm ],
            [ 1j * sin_t * np.conj(phase_pm), cos_t * np.conj(phase_pp) ]
        ], dtype=complex)
        return U


# =============================================================================
# KLASSE: GL - Allgemeine lineare Gruppe
# =============================================================================

class GL(MatrixLieGroup):
    """
    Allgemeine lineare Gruppe GL(n, F).

    GL(n, F) = { A ∈ M(n×n, F) | det(A) ≠ 0 }

    GL(n, F) ist die Gruppe aller invertierbaren n×n-Matrizen über dem
    Körper F (ℝ oder ℂ) unter Matrixmultiplikation.

    Eigenschaften:
    - dim(GL(n, ℝ)) = n²
    - GL(n, ℝ) ist nicht kompakt (Einträge können beliebig groß werden)
    - GL(n, ℝ) hat zwei Zusammenhangskomponenten (det > 0 und det < 0)
    - GL(n, ℂ) ist zusammenhängend

    @param n      Matrixgröße
    @param field  'real' oder 'complex'

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, n: int, field: str = 'real'):
        """
        Initialisiert GL(n, F).

        @param n      Matrixgröße (n ≥ 1)
        @param field  Körper: 'real' (ℝ) oder 'complex' (ℂ)
        @lastModified 2026-03-10
        """
        super().__init__(n, field=field)

    def is_element(self, A: np.ndarray, tol: float = 1e-10) -> bool:
        """
        Prüft, ob A ∈ GL(n, F).

        Bedingung: det(A) ≠ 0 (d.h. A ist invertierbar).

        @param A    zu prüfende Matrix (n×n)
        @param tol  Schwellwert: |det(A)| > tol gilt als invertierbar
        @return     True wenn A invertierbar ist
        @lastModified 2026-03-10
        """
        if A.shape != (self.n, self.n):
            return False
        return bool(abs(np.linalg.det(A)) > tol)


# =============================================================================
# KLASSE: LieAlgebra
# =============================================================================

class LieAlgebra:
    """
    Lie-Algebra g als Tangentialraum am Einselement einer Lie-Gruppe.

    Eine Lie-Algebra ist ein Vektorraum g mit einer bilinearen Operation
    [·, ·]: g × g → g (Lie-Klammer / Kommutator), die folgende Axiome erfüllt:
    1. Antisymmetrie:     [X, Y] = -[Y, X]
    2. Jacobi-Identität:  [X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0

    Für Matrixgruppen gilt [X, Y] = XY - YX (Matrixkommutator).

    @param basis_matrices  Liste von n×n-Matrizen, die eine Basis der Lie-Algebra bilden

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    def __init__(self, basis_matrices: list):
        """
        Initialisiert die Lie-Algebra mit einer Basis.

        @param basis_matrices  Liste von numpy-Arrays (alle gleicher Form n×n)
        @raises ValueError     wenn die Liste leer ist oder Matrizen verschiedener Form
        @lastModified 2026-03-10
        """
        if not basis_matrices:
            raise ValueError("basis_matrices darf nicht leer sein.")
        # Prüfe, dass alle Matrizen die gleiche Form haben
        shape = basis_matrices[0].shape
        if len(shape) != 2 or shape[0] != shape[1]:
            raise ValueError("Alle Basismatrizen müssen quadratisch sein.")
        for M in basis_matrices[1:]:
            if M.shape != shape:
                raise ValueError("Alle Basismatrizen müssen die gleiche Form haben.")
        self.basis = [np.array(M, dtype=complex) for M in basis_matrices]
        self.n = shape[0]        # Matrixgröße
        self.dim = len(self.basis)  # Dimension der Lie-Algebra

    def bracket(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        """
        Berechnet die Lie-Klammer [X, Y] = XY - YX.

        Die Lie-Klammer (Kommutator) ist die grundlegende Operation
        in einer Lie-Algebra. Sie misst, wie sehr X und Y nicht kommutieren.

        Eigenschaften:
        - [X, X] = 0
        - [X, Y] = -[Y, X]  (Antisymmetrie)
        - [[X,Y],Z] + [[Y,Z],X] + [[Z,X],Y] = 0  (Jacobi)

        @param X  erste Matrix aus der Lie-Algebra
        @param Y  zweite Matrix aus der Lie-Algebra
        @return   [X, Y] = XY - YX
        @lastModified 2026-03-10
        """
        X = np.array(X, dtype=complex)
        Y = np.array(Y, dtype=complex)
        return X @ Y - Y @ X

    def adjoint_representation(self, X: np.ndarray) -> np.ndarray:
        """
        Berechnet die adjungierte Darstellung ad_X: g → g.

        Die adjungierte Darstellung ist definiert durch:
            ad_X(Y) = [X, Y]

        In Matrixform (bezüglich der gewählten Basis {e_i}):
            (ad_X)_{ij} = Koeffizienten von [X, e_j] in der Basis {e_i}

        Algorithmus:
        1. Berechne [X, e_j] für jede Basismatrix e_j
        2. Drücke das Ergebnis als Linearkombination der Basis aus
           (via Frobenius-Skalarprodukt ⟨A, B⟩ = Tr(A†B))

        @param X  Element der Lie-Algebra (n×n-Matrix)
        @return   dim×dim-Matrix der adjungierten Darstellung ad_X
        @lastModified 2026-03-10
        """
        X = np.array(X, dtype=complex)
        d = self.dim
        # Adjungierte Darstellungsmatrix (d×d)
        ad_X = np.zeros((d, d), dtype=complex)
        for j, ej in enumerate(self.basis):
            # Berechne [X, e_j]
            bracket_result = self.bracket(X, ej)
            # Projiziere auf jeden Basisvektor via Frobenius-Skalarprodukt
            for i, ei in enumerate(self.basis):
                # ⟨e_i, [X, e_j]⟩ = Tr(e_i† · [X, e_j])
                inner = np.trace(ei.conj().T @ bracket_result)
                ad_X[i, j] = inner
        return ad_X

    def killing_form(self, X: np.ndarray, Y: np.ndarray) -> complex:
        """
        Berechnet die Killing-Form B(X, Y) = Tr(ad_X ∘ ad_Y).

        Die Killing-Form ist eine symmetrische Bilinearform auf der Lie-Algebra.
        Sie spielt eine zentrale Rolle in der Klassifikation einfacher Lie-Algebren:
        - B nicht-degeneriert ⟺ g ist halbeinfach (Cartan-Kriterium)
        - B negativ definit ⟺ g ist kompakt (für reelle Formen)

        @param X  erstes Element der Lie-Algebra
        @param Y  zweites Element der Lie-Algebra
        @return   skalarer Wert der Killing-Form B(X,Y) ∈ ℂ
        @lastModified 2026-03-10
        """
        ad_X = self.adjoint_representation(X)
        ad_Y = self.adjoint_representation(Y)
        # Killing-Form = Spur der Komposition
        return np.trace(ad_X @ ad_Y)

    def is_semisimple(self, tol: float = 1e-10) -> bool:
        """
        Prüft, ob die Lie-Algebra halbeinfach (semisimple) ist.

        Nach dem Cartan-Kriterium:
            g ist halbeinfach ⟺ die Killing-Form B ist nicht-degeneriert.

        Die Killing-Form ist nicht-degeneriert, wenn die Gramsche Matrix
        G_{ij} = B(e_i, e_j) eine von Null verschiedene Determinante hat.

        @param tol  Toleranz: |det(G)| > tol gilt als nicht-degeneriert
        @return     True wenn die Lie-Algebra halbeinfach ist
        @lastModified 2026-03-10
        """
        d = self.dim
        # Berechne Gramsche Matrix der Killing-Form
        gram = np.zeros((d, d), dtype=complex)
        for i, ei in enumerate(self.basis):
            for j, ej in enumerate(self.basis):
                gram[i, j] = self.killing_form(ei, ej)
        # Nicht-degeneriert ⟺ det ≠ 0
        return bool(abs(np.linalg.det(gram)) > tol)


# =============================================================================
# KLASSE: ExponentialMap
# =============================================================================

class ExponentialMap:
    """
    Exponentialabbildung exp: g → G und Logarithmus log: G → g.

    Die Exponentialabbildung ist die zentrale Verbindung zwischen
    einer Lie-Gruppe G und ihrer Lie-Algebra g:

        exp: g → G,  X ↦ e^X = Σ_{k=0}^∞ X^k / k!
        log: G → g,  A ↦ log(A)  (Umkehrung, lokal definiert)

    Für Matrixgruppen wird scipy.linalg.expm/logm verwendet (Padé-Approximation).

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    @staticmethod
    def exp_map(X: np.ndarray) -> np.ndarray:
        """
        Berechnet die Matrixexponentialfunktion e^X.

        Verwendet scipy.linalg.expm (Padé-Approximation mit Skalierung):
            e^X = Σ_{k=0}^∞ X^k / k!

        @param X  quadratische Matrix (Element der Lie-Algebra)
        @return   e^X (Element der Lie-Gruppe)
        @lastModified 2026-03-10
        """
        return expm(np.array(X, dtype=complex))

    @staticmethod
    def log_map(A: np.ndarray) -> np.ndarray:
        """
        Berechnet den Matrixlogarithmus log(A).

        Verwendet scipy.linalg.logm. Achtung: Der Logarithmus ist nur
        lokal (in einer Umgebung der Einheitsmatrix) eindeutig definiert.
        Für Matrizen mit negativen Eigenwerten kann das Ergebnis komplex sein.

        @param A  invertierbare quadratische Matrix (Element der Lie-Gruppe)
        @return   log(A) (Element der Lie-Algebra, näherungsweise)
        @lastModified 2026-03-10
        """
        return logm(np.array(A, dtype=complex))

    @staticmethod
    def baker_campbell_hausdorff(X: np.ndarray, Y: np.ndarray,
                                 order: int = 3) -> np.ndarray:
        """
        Berechnet die Baker-Campbell-Hausdorff (BCH) Formel bis zur angegebenen Ordnung.

        Die BCH-Formel gibt log(e^X · e^Y) als Lie-Reihe in X und Y:
            Z = X + Y
              + (1/2)[X, Y]
              + (1/12)([X, [X, Y]] - [Y, [X, Y]])
              + (1/24)([X, [Y, [X, Y]]])      (Ordnung 4, optional)
              + …

        BCH ist fundamental für:
        - Numerische Lie-Gruppen-Integration
        - Quantenmechanik (Trotter-Produkt-Formel)
        - Kontrolltheorie auf Lie-Gruppen

        @param X      erstes Lie-Algebra-Element
        @param Y      zweites Lie-Algebra-Element
        @param order  maximale Ordnung (1, 2 oder 3)
        @return       Approximation von log(e^X · e^Y)
        @lastModified 2026-03-10
        """
        X = np.array(X, dtype=complex)
        Y = np.array(Y, dtype=complex)

        # Hilfsfunktion für Kommutator
        def comm(A, B):
            """Kommutator [A, B] = AB - BA."""
            return A @ B - B @ A

        # Ordnung 1: Z ≈ X + Y
        Z = X + Y

        if order >= 2:
            # Ordnung 2: + (1/2)[X, Y]
            Z = Z + 0.5 * comm(X, Y)

        if order >= 3:
            # Ordnung 3: + (1/12)([X,[X,Y]] + [Y,[Y,X]])
            # = (1/12)([X,[X,Y]] - [Y,[X,Y]])
            cXY = comm(X, Y)
            Z = Z + (1.0 / 12.0) * (comm(X, cXY) + comm(Y, comm(Y, X)))

        if order >= 4:
            # Ordnung 4: - (1/24)[Y,[X,[X,Y]]]
            cXY = comm(X, Y)
            cXcXY = comm(X, cXY)
            Z = Z - (1.0 / 24.0) * comm(Y, cXcXY)

        return Z


# =============================================================================
# FREIE FUNKTIONEN: Einparametrige Untergruppe
# =============================================================================

def one_parameter_subgroup(X: np.ndarray, t: float) -> np.ndarray:
    """
    Berechnet die einparametrige Untergruppe exp(tX).

    Eine einparametrige Untergruppe ist ein Lie-Gruppen-Homomorphismus
        γ: ℝ → G,  t ↦ exp(tX)
    für ein festes X ∈ g (Lie-Algebra-Element).

    Eigenschaften:
    - γ(0) = I (Einheitsmatrix)
    - γ(s+t) = γ(s) · γ(t)  (Homomorphismus-Eigenschaft)
    - γ'(0) = X  (X ist der Tangentialvektor bei t=0)

    @param X  Lie-Algebra-Element (quadratische Matrix)
    @param t  Parameter (reell)
    @return   exp(t·X) als numpy-Array
    @lastModified 2026-03-10
    """
    return expm(t * np.array(X, dtype=complex))


# =============================================================================
# KLASSIFIKATION: Einfache Lie-Algebren (Cartan)
# =============================================================================

def classify_simple_lie_algebra(n: int, type_char: str) -> dict:
    """
    Klassifiziert eine einfache Lie-Algebra nach Cartan (klassische Typen).

    Die Cartan-Klassifikation teilt einfache komplexe Lie-Algebren ein in:
    - A_n (n≥1): sl(n+1, ℂ),  Lie-Algebra von SL(n+1)
    - B_n (n≥2): so(2n+1, ℂ), Lie-Algebra von SO(2n+1)
    - C_n (n≥3): sp(2n, ℂ),   Lie-Algebra von Sp(2n)
    - D_n (n≥4): so(2n, ℂ),   Lie-Algebra von SO(2n)
    Plus 5 Ausnahmen: G₂, F₄, E₆, E₇, E₈

    @param n          Rang der Lie-Algebra (≥ 1)
    @param type_char  'A', 'B', 'C' oder 'D'
    @return           dict mit Klassifikationsinformationen
    @raises ValueError bei ungültigem Typ oder Rang
    @lastModified 2026-03-10
    """
    type_char = type_char.upper()
    if type_char not in ('A', 'B', 'C', 'D'):
        raise ValueError(f"Ungültiger Typ '{type_char}'. Erlaubt: A, B, C, D.")

    info = {
        'type': type_char,
        'rank': n,
        'name': f"{type_char}_{n}",
    }

    if type_char == 'A':
        # A_n = sl(n+1): spurlose (n+1)×(n+1)-Matrizen
        if n < 1:
            raise ValueError("A_n erfordert n ≥ 1.")
        info['lie_group'] = f"SL({n+1})"
        info['dimension'] = (n + 1) ** 2 - 1   # = n² + 2n
        info['lie_algebra'] = f"sl({n+1}, ℂ)"
        info['description'] = (
            f"Spurlose ({n+1})×({n+1})-Matrizen; "
            f"Dynkin-Diagramm: lineare Kette von {n} Knoten (A-Typ)"
        )
        info['roots_positive'] = n * (n + 1) // 2

    elif type_char == 'B':
        # B_n = so(2n+1): (2n+1)×(2n+1) schiefsymmetrische Matrizen
        if n < 2:
            raise ValueError("B_n erfordert n ≥ 2.")
        info['lie_group'] = f"SO({2*n+1})"
        info['dimension'] = n * (2 * n + 1)
        info['lie_algebra'] = f"so({2*n+1}, ℂ)"
        info['description'] = (
            f"Schiefsymmetrische ({2*n+1})×({2*n+1})-Matrizen; "
            f"Dynkin-Diagramm: {n-1} einfache Kanten + eine Doppelkante rechts"
        )
        info['roots_positive'] = n * n

    elif type_char == 'C':
        # C_n = sp(2n): symplektische Lie-Algebra
        if n < 3:
            raise ValueError("C_n erfordert n ≥ 3.")
        info['lie_group'] = f"Sp({2*n})"
        info['dimension'] = n * (2 * n + 1)
        info['lie_algebra'] = f"sp({2*n}, ℂ)"
        info['description'] = (
            f"Symplektische Lie-Algebra; Matrizen die die symplektische Form erhalten; "
            f"Dynkin-Diagramm: {n-1} einfache Kanten + eine Doppelkante links"
        )
        info['roots_positive'] = n * n

    elif type_char == 'D':
        # D_n = so(2n): gerades orthogonales
        if n < 4:
            raise ValueError("D_n erfordert n ≥ 4.")
        info['lie_group'] = f"SO({2*n})"
        info['dimension'] = n * (2 * n - 1)
        info['lie_algebra'] = f"so({2*n}, ℂ)"
        info['description'] = (
            f"Schiefsymmetrische ({2*n})×({2*n})-Matrizen; "
            f"Dynkin-Diagramm: {n-2} einfache Kanten + Gabelung am Ende (D-Typ)"
        )
        info['roots_positive'] = n * (n - 1)

    return info


def exceptional_lie_algebras_info() -> dict:
    """
    Gibt Informationen über die 5 außergewöhnlichen (exzeptionellen) Lie-Algebren zurück.

    Die exzeptionellen einfachen Lie-Algebren G₂, F₄, E₆, E₇, E₈ passen
    nicht in die Cartan-Klassifikation A, B, C, D. Sie wurden von
    Wilhelm Killing (1888) und Élie Cartan (1894) entdeckt.

    @return dict mit Informationen zu G2, F4, E6, E7, E8
    @lastModified 2026-03-10
    """
    return {
        'G2': {
            'rank': 2,
            'dimension': 14,
            'roots_positive': 6,
            'description': 'Automorphismen des Oktonionenkörpers ℚ; kleinste exzeptionelle Lie-Algebra',
            'dynkin': '⬤═══⬤ (Doppelkante mit Pfeil)',
            'compact_form': 'G₂ (kompakte reelle Form)',
        },
        'F4': {
            'rank': 4,
            'dimension': 52,
            'roots_positive': 24,
            'description': 'Isometriegruppe des Oktonionischen projektiven Raums OP²',
            'dynkin': '⬤─⬤═⬤─⬤',
            'compact_form': 'F₄',
        },
        'E6': {
            'rank': 6,
            'dimension': 78,
            'roots_positive': 36,
            'description': 'Erscheint in F-Theorie und Heterotischer Stringtheorie',
            'dynkin': 'Gabeldiagramm mit 6 Knoten (5 in Linie + 1 oben an Position 3)',
            'compact_form': 'E₆',
        },
        'E7': {
            'rank': 7,
            'dimension': 133,
            'roots_positive': 63,
            'description': '56-dimensionale Minimaldarstellung; wichtig in M-Theorie',
            'dynkin': 'Gabeldiagramm mit 7 Knoten',
            'compact_form': 'E₇',
        },
        'E8': {
            'rank': 8,
            'dimension': 248,
            'roots_positive': 120,
            'description': 'Größte exzeptionelle Lie-Algebra; Eichgruppe der Heterotischen Stringtheorie',
            'dynkin': 'Gabeldiagramm mit 8 Knoten',
            'compact_form': 'E₈',
            'note': 'E₈ × E₈ ist eine der Eichgruppen der Heterotischen String-Theorie',
        },
    }


def cartan_matrix(type_char: str, n: int) -> np.ndarray:
    """
    Berechnet die Cartan-Matrix für klassische Lie-Algebra-Typen A_n, B_n, C_n, D_n.

    Die Cartan-Matrix A_{ij} = 2⟨α_i, α_j⟩/⟨α_j, α_j⟩ kodiert die
    Winkelbeziehungen zwischen den einfachen Wurzeln α_1, ..., α_n.

    Diagonale: immer 2 (da ⟨α_i, α_i⟩/⟨α_i, α_i⟩ = 1, mal 2 = 2)
    Außerdiagonale: 0, -1, -2 oder -3

    @param type_char  'A', 'B', 'C' oder 'D'
    @param n          Rang (Anzahl einfacher Wurzeln)
    @return           n×n Cartan-Matrix
    @raises ValueError bei ungültigem Typ oder Rang
    @lastModified 2026-03-10
    """
    type_char = type_char.upper()
    if type_char not in ('A', 'B', 'C', 'D'):
        raise ValueError(f"Ungültiger Typ '{type_char}'.")

    # Initialisiere Cartan-Matrix mit 2 auf der Diagonalen
    A = np.zeros((n, n), dtype=int)
    np.fill_diagonal(A, 2)

    if type_char == 'A':
        # A_n: Dynkin-Diagramm ist eine lineare Kette
        # Alle Bindungen sind einfach (Winkel = 120°)
        for i in range(n - 1):
            A[i, i+1] = -1
            A[i+1, i] = -1

    elif type_char == 'B':
        # B_n: lineare Kette mit Doppelkante am Ende (rechts)
        # B_n: α_1─α_2─…─α_{n-1}⟹α_n (lange→kurze Wurzel)
        for i in range(n - 2):
            A[i, i+1] = -1
            A[i+1, i] = -1
        if n >= 2:
            # Letzte Bindung ist Doppelkante: -1 (lang→kurz) und -2 (kurz→lang)
            A[n-2, n-1] = -1
            A[n-1, n-2] = -2

    elif type_char == 'C':
        # C_n: lineare Kette mit Doppelkante am Anfang (links)
        # C_n: α_1⟸α_2─…─α_n (kurze→lange Wurzel)
        for i in range(n - 2):
            A[i, i+1] = -1
            A[i+1, i] = -1
        if n >= 2:
            # Erste Bindung Doppelkante: -2 (kurz→lang) und -1 (lang→kurz)
            A[n-2, n-1] = -2
            A[n-1, n-2] = -1

    elif type_char == 'D':
        # D_n: Gabelung am Ende
        # D_n: α_1─α_2─…─α_{n-2}─α_{n-1}
        #                           └─α_n
        if n < 2:
            raise ValueError("D_n erfordert n ≥ 2.")
        for i in range(n - 3):
            A[i, i+1] = -1
            A[i+1, i] = -1
        if n >= 3:
            # Gabelungsknoten α_{n-2} verbindet mit α_{n-1} und α_n
            A[n-3, n-2] = -1
            A[n-2, n-3] = -1
            A[n-3, n-1] = -1
            A[n-1, n-3] = -1

    return A


def dynkin_diagram_info(type_char: str, n: int) -> dict:
    """
    Gibt eine textuelle Beschreibung des Dynkin-Diagramms zurück.

    Dynkin-Diagramme kodieren die Geometrie der Wurzelsysteme einfacher
    Lie-Algebren:
    - Knoten = einfache Wurzeln
    - Einfache Kante (─): Winkel 120° zwischen Wurzeln
    - Doppelkante (═): Winkel 135° (ein Wurzelpaar)
    - Dreifachkante (≡): Winkel 150° (G₂, nicht hier)
    - Pfeil zeigt von längerer zu kürzerer Wurzel

    @param type_char  'A', 'B', 'C' oder 'D'
    @param n          Rang
    @return           dict mit Diagramm-Beschreibung und Informationen
    @lastModified 2026-03-10
    """
    type_char = type_char.upper()
    info = {
        'type': type_char,
        'rank': n,
        'nodes': n,
    }

    if type_char == 'A':
        # A_n: ○─○─○─…─○ (n Knoten in Reihe)
        info['diagram'] = '─'.join(['○'] * n)
        info['edges'] = [(i, i+1, 'single') for i in range(1, n)]
        info['description'] = f"Lineare Kette von {n} Knoten (alle Wurzeln gleich lang)"
        info['simply_laced'] = True

    elif type_char == 'B':
        # B_n: ○─○─…─○⟹○ (n-1 einfache Kanten + 1 Doppelkante rechts)
        chain = '─'.join(['○'] * (n - 1))
        info['diagram'] = chain + '⟹○' if n > 1 else '○'
        info['edges'] = (
            [(i, i+1, 'single') for i in range(1, n - 1)] +
            [(n-1, n, 'double')]
        )
        info['description'] = (
            f"Kette von {n} Knoten; letzte Kante doppelt mit Pfeil (lang→kurz); "
            f"kurze Endwurzel"
        )
        info['simply_laced'] = False

    elif type_char == 'C':
        # C_n: ○⟸○─…─○─○ (Doppelkante links, dann einfache Kanten)
        chain = '─'.join(['○'] * (n - 1))
        info['diagram'] = '○⟸' + chain if n > 1 else '○'
        info['edges'] = (
            [(1, 2, 'double')] +
            [(i, i+1, 'single') for i in range(2, n)]
        )
        info['description'] = (
            f"Kette von {n} Knoten; erste Kante doppelt mit Pfeil (kurz→lang); "
            f"lange Endwurzel"
        )
        info['simply_laced'] = False

    elif type_char == 'D':
        # D_n: Gabelung am vorletzten Knoten
        if n >= 3:
            chain = '─'.join(['○'] * (n - 2))
            info['diagram'] = chain + '─<○ (Gabelung)'
        else:
            info['diagram'] = '○─○' if n == 2 else '○'
        info['edges'] = (
            [(i, i+1, 'single') for i in range(1, n - 2)] +
            [(n-2, n-1, 'single'), (n-2, n, 'single')]
        )
        info['description'] = (
            f"Kette von {n-2} Knoten mit Gabelung am ({n-2})-ten Knoten; "
            f"alle Wurzeln gleich lang (simply laced)"
        )
        info['simply_laced'] = True

    else:
        raise ValueError(f"Ungültiger Typ '{type_char}'.")

    return info


# =============================================================================
# DARSTELLUNGSTHEORIE: SU(2) Spin-Darstellungen
# =============================================================================

def fundamental_representation_su2(j: float) -> dict:
    """
    Berechnet die Spin-j-Darstellung (fundamentale Darstellung) von SU(2).

    SU(2) hat für jedes j ∈ {0, 1/2, 1, 3/2, …} genau eine irreduzible
    Darstellung der Dimension 2j+1. Diese werden durch die Drehimpuls-Operatoren
    J_x, J_y, J_z (oder J_+, J_-, J_z) beschrieben.

    Basis: |j, m⟩ für m = -j, -j+1, …, j-1, j

    Matrixelemente:
        ⟨j, m'|J_z|j, m⟩ = m · δ_{m',m}
        ⟨j, m'|J_+|j, m⟩ = √(j(j+1) - m(m+1)) · δ_{m',m+1}
        ⟨j, m'|J_-|j, m⟩ = √(j(j+1) - m(m-1)) · δ_{m',m-1}

    @param j  Spin-Quantenzahl: 0, 1/2, 1, 3/2, 2, ...
    @return   dict mit Operatoren Jx, Jy, Jz, J+, J- und Metadaten
    @raises   ValueError wenn j < 0 oder nicht ganzzahlig/halbganzzahlig
    @lastModified 2026-03-10
    """
    # Validierung: j muss nicht-negativ und ganzzahlig oder halbganzzahlig sein
    if j < 0:
        raise ValueError("Spin j muss ≥ 0 sein.")
    two_j = int(round(2 * j))
    if abs(two_j - 2 * j) > 1e-10:
        raise ValueError("j muss ganzzahlig oder halbganzzahlig sein (j = 0, 1/2, 1, 3/2, ...).")

    # Dimension der Darstellung: d = 2j + 1
    d = two_j + 1
    # Magnetische Quantenzahlen m: von -j bis +j in Einerschritten
    m_values = [j - k for k in range(d)]  # j, j-1, ..., -j

    # Initialisiere Matrizen
    Jz = np.zeros((d, d), dtype=complex)
    Jplus = np.zeros((d, d), dtype=complex)
    Jminus = np.zeros((d, d), dtype=complex)

    for idx, m in enumerate(m_values):
        # J_z-Diagonalelement: ⟨m|J_z|m⟩ = m (in Einheiten von ℏ)
        Jz[idx, idx] = m

        # J_+ hebt m um 1 an: |m⟩ → |m+1⟩
        # ⟨m+1|J_+|m⟩ = √(j(j+1) - m(m+1))
        if idx > 0:  # m+1 existiert (idx-1 entspricht m+1)
            m_up = m_values[idx - 1]  # = m + 1
            coeff = np.sqrt(j * (j + 1) - m * m_up)
            Jplus[idx - 1, idx] = coeff

        # J_- senkt m um 1 ab: |m⟩ → |m-1⟩
        # ⟨m-1|J_-|m⟩ = √(j(j+1) - m(m-1))
        if idx < d - 1:  # m-1 existiert (idx+1 entspricht m-1)
            m_down = m_values[idx + 1]  # = m - 1
            coeff = np.sqrt(j * (j + 1) - m * m_down)
            Jminus[idx + 1, idx] = coeff

    # J_x = (J_+ + J_-) / 2
    Jx = (Jplus + Jminus) / 2.0
    # J_y = (J_+ - J_-) / (2i)
    Jy = (Jplus - Jminus) / (2.0j)

    return {
        'j': j,
        'dimension': d,
        'm_values': m_values,
        'Jz': Jz,
        'Jx': Jx,
        'Jy': Jy,
        'Jplus': Jplus,
        'Jminus': Jminus,
        'casimir_eigenvalue': j * (j + 1),  # J² = j(j+1)·I
    }


# =============================================================================
# CASIMIR-OPERATOR
# =============================================================================

def casimir_element(basis_matrices: list) -> np.ndarray:
    """
    Berechnet den quadratischen Casimir-Operator einer Lie-Algebra-Darstellung.

    Der Casimir-Operator ist definiert als:
        C₂ = Σ_i e_i · e_i  (Summe über Basiselemente, Matrizenprodukt)

    wobei die Basis {e_i} bezüglich der Killing-Form dual gewählt wird.
    Für kompakte halbeinfache Lie-Algebren (mit Standard-Normierung) gilt:
        C₂ = Σ_i (e_i)²

    Der Casimir-Operator pendelt mit allen Elementen der Algebra und ist
    daher nach dem Schur-Lemma in irreduziblen Darstellungen ein Vielfaches
    der Einheitsmatrix.

    Für SU(2): C₂ = J_x² + J_y² + J_z² = j(j+1)·I

    @param basis_matrices  Liste der Lie-Algebra-Basismatrizen
    @return                Casimir-Matrix als numpy-Array
    @lastModified 2026-03-10
    """
    if not basis_matrices:
        raise ValueError("basis_matrices darf nicht leer sein.")

    # Initialisiere mit Nullmatrix gleicher Form
    n = basis_matrices[0].shape[0]
    C = np.zeros((n, n), dtype=complex)

    # Summiere quadratische Beiträge: C₂ = Σ_i e_i²
    for e in basis_matrices:
        e = np.array(e, dtype=complex)
        C = C + e @ e

    return C
