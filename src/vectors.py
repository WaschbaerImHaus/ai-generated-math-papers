"""
@file vectors.py
@brief Vektor-Modul: n-dimensionale Vektoren und Gram-Schmidt-Orthogonalisierung.
@description
    Enthält die Vector-Klasse für n-dimensionale reelle Vektoren sowie
    die gram_schmidt()-Funktion zur Orthogonalisierung einer Vektormenge.

    Operationen:
    - Skalarprodukt (dot), Kreuzprodukt (cross, nur 3D), Norm, Normierung
    - Addition, Subtraktion, Skalarmultiplikation
    - Gram-Schmidt-Orthogonalisierung (klassisches Verfahren)

    Ausgelagert aus linear_algebra.py für bessere Modularität.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math


# =============================================================================
# KLASSE: Vector
# =============================================================================

class Vector:
    """
    @brief n-dimensionaler reeller Vektor.
    @description
        Unterstützt: Skalarprodukt, Kreuzprodukt (3D), Norm,
        Normierung, Addition, Skalarmultiplikation.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, components: list[float]) -> None:
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

        Beispiele:
        >>> Vector([1, 0, 0]).dot(Vector([0, 1, 0]))
        0
        >>> Vector([1, 2, 3]).dot(Vector([4, 5, 6]))
        32
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

        Beispiele:
        >>> Vector([3, 4]).norm()
        5.0
        >>> Vector([1, 0, 0]).norm()
        1.0
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

    def __len__(self) -> int:
        """
        @brief Gibt die Dimension des Vektors zurück.
        @return Anzahl der Komponenten.
        @date 2026-03-10
        """
        return self.dim

    def __getitem__(self, idx: int) -> float:
        """
        @brief Gibt die Komponente am Index idx zurück.
        @param idx 0-basierter Index der Komponente.
        @return Komponentenwert.
        @date 2026-03-10
        """
        return self.components[idx]

    def __repr__(self) -> str:
        return f"Vector({self.components})"


# =============================================================================
# GRAM-SCHMIDT-ORTHOGONALISIERUNG
# =============================================================================

def gram_schmidt(vectors: list['Vector'], normalize: bool = False) -> list['Vector']:
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
    @author Kurt Ingwer
    @lastModified 2026-03-10
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
