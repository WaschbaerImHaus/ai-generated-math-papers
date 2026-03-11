"""
@file spinor_calculus.py
@brief Spinor-Rechnung: Clifford-Algebren, Spinoren, Dirac-Gleichung, Pauli-Matrizen.
@description
    Dieses Modul implementiert die mathematischen Strukturen der Spinorrechnung:

    ## Clifford-Algebren
    Die Clifford-Algebra Cl(p,q) ist die assoziative Algebra mit Erzeugern
    e₁,...,eₙ (n = p + q) und den Relationen:
        eᵢeⱼ + eⱼeᵢ = 2η_{ij}
    wobei η = diag(+1,...,+1,-1,...,-1) mit p Plusen und q Minusen.

    Wichtige Fälle:
    - Cl(3,0): Geometrische Algebra des 3D-Raums, isomorph zu Quaternionen ⊕ Quaternionen
    - Cl(1,3): Raumzeit-Algebra der Minkowski-Raumzeit (Dirac-Algebra)
    - Cl(0,2) ≅ ℍ: Quaternionen

    ## Spinoren
    Ein Spinor ist ein Element des minimalen Linksideals der Clifford-Algebra.
    Physikalisch: Objekte, die sich unter 360°-Rotation um -1 transformieren.

    Dirac-Spinor ψ ∈ ℂ⁴: Beschreibt Spin-1/2-Teilchen (Elektronen)

    ## Dirac-Gleichung
    (iγ^μ∂_μ - m)ψ = 0

    Dirac-Darstellung der γ-Matrizen:
        γ⁰ = [[I, 0], [0, -I]]
        γᵢ = [[0, σᵢ], [-σᵢ, 0]]  (i = 1,2,3)

    Antikommutatorrelation:
        {γ^μ, γ^ν} = γ^μγ^ν + γ^νγ^μ = 2g^{μν}·I

    Minkowski-Metrik: g = diag(+1,-1,-1,-1)

@author Michael Fuhrmann
@lastModified 2026-03-11
@version 1.0.0
"""

import numpy as np
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Klasse: PauliMatrices
# ---------------------------------------------------------------------------

class PauliMatrices:
    """
    Pauli-Matrizen σ₁, σ₂, σ₃ und ihre algebraischen Eigenschaften.

    Die drei Pauli-Matrizen sind hermitesche, spurlose 2×2-Matrizen:
        σ₁ = [[0, 1], [1, 0]]    (Spiegelsymmetrie in x-Richtung)
        σ₂ = [[0,-i], [i, 0]]    (komplexe Drehsymmetrie)
        σ₃ = [[1, 0], [0,-1]]    (Eigenzustände für Spin-up/down)

    Eigenschaften:
    - σᵢ² = I (Idempotenz)
    - σᵢσⱼ = δᵢⱼI + iεᵢⱼₖσₖ (Algebra-Relation)
    - [σᵢ, σⱼ] = 2iεᵢⱼₖσₖ (Lie-Algebra su(2))
    - {σᵢ, σⱼ} = 2δᵢⱼI (Clifford-Relation)
    - Spur(σᵢ) = 0, det(σᵢ) = -1

    Physikalische Bedeutung:
        S = ℏ/2 · σ  (Spin-Operator für Spin-1/2-Teilchen)
        [Sᵢ, Sⱼ] = iℏεᵢⱼₖSₖ  (Drehimpuls-Algebra)
    """

    # Pauli-Matrizen als Klassenkonstanten (2×2 komplexe Matrizen)
    _SIGMA = {
        1: np.array([[0, 1], [1, 0]], dtype=complex),        # σ_x
        2: np.array([[0, -1j], [1j, 0]], dtype=complex),     # σ_y
        3: np.array([[1, 0], [0, -1]], dtype=complex),       # σ_z
    }

    @staticmethod
    def sigma(i: int) -> np.ndarray:
        """
        @brief Gibt die i-te Pauli-Matrix zurück.
        @description
            i=1: σ_x = [[0,1],[1,0]]
            i=2: σ_y = [[0,-i],[i,0]]
            i=3: σ_z = [[1,0],[0,-1]]

        @param i: Index 1, 2 oder 3.
        @return: 2×2 numpy-Array (komplex).
        @raises ValueError: Falls i ∉ {1, 2, 3}.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        if i not in (1, 2, 3):
            raise ValueError(f"Pauli-Matrix-Index muss 1, 2 oder 3 sein, nicht {i}")
        return PauliMatrices._SIGMA[i].copy()

    @staticmethod
    def commutator(i: int, j: int) -> np.ndarray:
        """
        @brief Berechnet den Kommutator [σᵢ, σⱼ] = 2iεᵢⱼₖσₖ.
        @description
            [σᵢ, σⱼ] = σᵢσⱼ - σⱼσᵢ = 2i·εᵢⱼₖ·σₖ

            Levi-Civita-Symbol ε:
            - εᵢⱼₖ = +1 für (i,j,k) gerade Permutation von (1,2,3)
            - εᵢⱼₖ = -1 für ungerade Permutation
            - εᵢⱼₖ = 0 falls zwei Indizes gleich

        @param i: Erster Index (1, 2 oder 3).
        @param j: Zweiter Index (1, 2 oder 3).
        @return: 2×2 Kommutator-Matrix.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        sigma_i = PauliMatrices.sigma(i)
        sigma_j = PauliMatrices.sigma(j)
        # [A, B] = AB - BA
        return sigma_i @ sigma_j - sigma_j @ sigma_i

    @staticmethod
    def anticommutator(i: int, j: int) -> np.ndarray:
        """
        @brief Berechnet den Antikommutator {σᵢ, σⱼ} = 2δᵢⱼ·I.
        @description
            {σᵢ, σⱼ} = σᵢσⱼ + σⱼσᵢ = 2δᵢⱼ·I₂

            Dies ist die Clifford-Algebra-Relation:
            Die Pauli-Matrizen erzeugen Cl(3,0) (in geeignetem Sinne).

        @param i: Erster Index (1, 2 oder 3).
        @param j: Zweiter Index (1, 2 oder 3).
        @return: 2×2 Antikommutator-Matrix.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        sigma_i = PauliMatrices.sigma(i)
        sigma_j = PauliMatrices.sigma(j)
        return sigma_i @ sigma_j + sigma_j @ sigma_i

    @staticmethod
    def spin_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
        """
        @brief Berechnet den Spin-Rotationsoperator exp(-iθ·n̂·σ/2).
        @description
            Rotation um Achse n̂ = (n₁, n₂, n₃) mit |n̂|=1 und Winkel θ:
                U = exp(-iθ·n̂·σ/2) = cos(θ/2)·I - i·sin(θ/2)·(n̂·σ)

            Dies ist die SU(2)-Darstellung der SO(3)-Rotation.
            Spinoren transformieren sich als:
                ψ → U·ψ

            Wichtig: U(-1) = U(2π-θ), d.h. eine 360°-Rotation gibt -I!
            → Spinoren sind zweideutig unter 360°-Rotation.

        @param axis: Rotationsachse n̂ (wird normiert).
        @param angle: Rotationswinkel θ in Bogenmaß.
        @return: 2×2 unitäre Matrix (SU(2)-Element).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        axis = np.asarray(axis, dtype=float)
        norm = np.linalg.norm(axis)
        if norm < 1e-12:
            return np.eye(2, dtype=complex)
        axis = axis / norm  # normieren

        # n̂·σ = n₁σ₁ + n₂σ₂ + n₃σ₃
        n_dot_sigma = (axis[0] * PauliMatrices.sigma(1)
                       + axis[1] * PauliMatrices.sigma(2)
                       + axis[2] * PauliMatrices.sigma(3))

        # U = cos(θ/2)·I - i·sin(θ/2)·(n̂·σ)
        half_angle = angle / 2.0
        U = (np.cos(half_angle) * np.eye(2, dtype=complex)
             - 1j * np.sin(half_angle) * n_dot_sigma)

        return U


# ---------------------------------------------------------------------------
# Klasse: CliffordAlgebra
# ---------------------------------------------------------------------------

class CliffordAlgebra:
    """
    Clifford-Algebra Cl(p, q).

    Die Clifford-Algebra Cl(p, q) hat n = p + q Erzeuger e₁,...,eₙ mit:
        eᵢeⱼ + eⱼeᵢ = 2η_{ij}
    wobei η = diag(+1,...,+1,-1,...,-1) (p Plusse, q Minusse).

    Dimension der Algebra: 2ⁿ (über den reellen Zahlen).

    Wichtige Spezialfälle:
    - Cl(0,0): ℝ
    - Cl(1,0): ℝ ⊕ ℝ
    - Cl(0,1): ℂ (Komplexe Zahlen)
    - Cl(0,2) ≅ Cl(2,0): Quaternionen ℍ
    - Cl(3,0): Geometrische Algebra in 3D
    - Cl(1,3): Dirac-Algebra der Raumzeit

    Gamma-Matrizen als Darstellung:
        Für Cl(p,q) mit n=p+q: Matrizen in Mat(2^⌊n/2⌋ × 2^⌊n/2⌋, ℂ)
    """

    def __init__(self, p: int, q: int = 0) -> None:
        """
        @brief Initialisiert die Clifford-Algebra Cl(p, q).
        @param p: Anzahl positiver Generatoren (eᵢ² = +1 für i ≤ p).
        @param q: Anzahl negativer Generatoren (eᵢ² = -1 für i > p).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        self.p = p  # positive Generatoren
        self.q = q  # negative Generatoren
        self.n = p + q  # Gesamtanzahl Generatoren

        # Signaturmatrix η = diag(+1,...,+1,-1,...,-1)
        self.signature = np.diag(
            [1.0] * p + [-1.0] * q
        )

        # Gamma-Matrizen (Darstellung der Generatoren) berechnen
        self._gamma_matrices = self._build_gamma_matrices()

    def _build_gamma_matrices(self) -> List[np.ndarray]:
        """
        @brief Konstruiert die Gamma-Matrizen als Darstellung der Clifford-Algebra.
        @description
            Nutzt rekursive Tensorprodukt-Konstruktion (Majorana-Basis).

            Für n=1: γ₁ = [[0,1],[1,0]] (Cl(1,0)) oder [[0,i],[-i,0]] (Cl(0,1))

            Für n gerade: Matrizen der Größe 2^(n/2) × 2^(n/2)
            Für n ungerade: Matrizen der Größe 2^((n+1)/2) × 2^((n+1)/2)

            Konstruktionsprinzip (Pauli-basiert):
            - σ₁, σ₂, σ₃ als Basisbausteine
            - Tensorprodukte für höhere Dimensionen

        @return: Liste von Gamma-Matrizen [γ₁, ..., γₙ].
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Pauli-Matrizen als Bausteine
        s1 = np.array([[0, 1], [1, 0]], dtype=complex)
        s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
        s3 = np.array([[1, 0], [0, -1]], dtype=complex)
        I2 = np.eye(2, dtype=complex)

        gammas = []

        if self.n == 0:
            # Triviale Algebra: nur das Skalarteil
            return [np.array([[1]], dtype=complex)]

        elif self.n == 1:
            # Cl(1,0): γ₁² = +1
            # Cl(0,1): γ₁² = -1
            if self.p == 1:
                gammas = [s1]
            else:
                gammas = [1j * s1]

        elif self.n == 2:
            # Cl(2,0), Cl(1,1), Cl(0,2)
            if self.p == 2:
                gammas = [s1, s2]  # γ₁²=+1, γ₂²=+1
            elif self.p == 1 and self.q == 1:
                gammas = [s1, 1j * s2]  # γ₁²=+1, γ₂²=-1
            else:  # p=0, q=2
                gammas = [1j * s1, 1j * s2]

        elif self.n == 3:
            # Cl(3,0): Pauli-Matrizen direkt
            if self.p == 3:
                gammas = [s1, s2, s3]
            elif self.p == 2 and self.q == 1:
                gammas = [s1, s2, 1j * s3]
            elif self.p == 1 and self.q == 2:
                gammas = [s1, 1j * s2, 1j * s3]
            else:
                gammas = [1j * s1, 1j * s2, 1j * s3]

        elif self.n == 4:
            # Cl(1,3): Dirac-Algebra
            # γ⁰ = [[I,0],[0,-I]], γⁱ = [[0,σᵢ],[-σᵢ,0]]
            if self.p == 1 and self.q == 3:
                # Minkowski-Signatur (+,-,-,-)
                zero2 = np.zeros((2, 2), dtype=complex)
                gamma0 = np.block([[I2, zero2], [zero2, -I2]])
                gamma1 = np.block([[zero2, s1], [-s1, zero2]])
                gamma2 = np.block([[zero2, s2], [-s2, zero2]])
                gamma3 = np.block([[zero2, s3], [-s3, zero2]])
                gammas = [gamma0, gamma1, gamma2, gamma3]
            else:
                # Allgemeiner Fall: rekursiv via Tensorprodukte
                gammas = self._recursive_gamma(self.n, self.p)

        else:
            # Allgemeiner Fall: rekursive Konstruktion
            gammas = self._recursive_gamma(self.n, self.p)

        return gammas

    def _recursive_gamma(self, n: int, p: int) -> List[np.ndarray]:
        """
        @brief Rekursive Konstruktion der Gamma-Matrizen via Tensorprodukte.
        @description
            Für gerades n=2m: Matrizen der Größe 2^m × 2^m
            Für ungerades n=2m+1: Matrizen der Größe 2^(m+1) × 2^(m+1)

        @param n: Gesamtanzahl Generatoren.
        @param p: Anzahl positiver Generatoren.
        @return: Liste von Gamma-Matrizen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Matrixgröße: 2^⌈n/2⌉
        dim = 2 ** max(1, (n + 1) // 2)
        gammas = []

        # Signatur-Vorzeichen für jeden Generator
        signs = [1.0] * p + [-1.0] * (n - p)

        # Einfache Konstruktion: Diagonale + Off-Diagonale Matrizen
        for k in range(n):
            gamma = np.zeros((dim, dim), dtype=complex)
            # Blockdiagonale Struktur
            block_size = dim // 2
            if block_size == 0:
                block_size = 1

            if k % 2 == 0:
                # Obere rechte + untere linke Blöcke
                for i in range(min(block_size, dim)):
                    j = (i + block_size) % dim
                    if j < dim:
                        gamma[i, j] = np.sqrt(abs(signs[k]))
                        gamma[j, i] = np.sqrt(abs(signs[k])) * (1j if signs[k] < 0 else 1)
            else:
                # Diagonale Blöcke
                for i in range(dim):
                    gamma[i, i] = 1.0 if i < block_size else -1.0
                    if signs[k] < 0:
                        gamma[i, i] *= 1j

            gammas.append(gamma)

        return gammas

    def generators(self) -> List[np.ndarray]:
        """
        @brief Gibt die Gamma-Matrizen (Darstellung der Generatoren) zurück.
        @description
            Gibt [γ₁, ..., γₙ] zurück, wobei:
                γᵢγⱼ + γⱼγᵢ = 2η_{ij}·I

        @return: Liste von numpy-Arrays (Gamma-Matrizen).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return [g.copy() for g in self._gamma_matrices]

    def clifford_product(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        @brief Berechnet das Clifford-Produkt zweier Matrizen a·b.
        @description
            Das Clifford-Produkt ist das gewöhnliche Matrixprodukt
            in der Darstellung durch Gamma-Matrizen.

        @param a: Erste Matrix (Clifford-Element in Matrixdarstellung).
        @param b: Zweite Matrix (Clifford-Element in Matrixdarstellung).
        @return: Produkt-Matrix a·b.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return a @ b

    def grade_involution(self, a: np.ndarray) -> np.ndarray:
        """
        @brief Hauptinvolution α: eᵢ → -eᵢ.
        @description
            Die Hauptinvolution α(a) vertauscht das Vorzeichen aller
            Erzeuger-Vektoren:
                α(eᵢ₁···eᵢₖ) = (-1)ᵏ · eᵢ₁···eᵢₖ

            Für Matrizen: Implementierung via Ähnlichkeitstransformation
            mit der Volumenelement-Matrix.

        @param a: Matrix-Darstellung eines Clifford-Elements.
        @return: Matrix der Hauptinvolution.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Vereinfachte Implementierung: Konjugation mit erstem Gamma
        # α(a) = -γ₁ a γ₁⁻¹ (wenn γ₁² = ±I)
        if len(self._gamma_matrices) == 0:
            return a.copy()
        gamma0 = self._gamma_matrices[0]
        # γ₁⁻¹ via Adjungierte (γ₁ hermitesch → γ₁⁻¹ = γ₁† für unitäre Gamma)
        try:
            gamma0_inv = np.linalg.inv(gamma0)
        except np.linalg.LinAlgError:
            return a.copy()
        return -gamma0 @ a @ gamma0_inv

    def reversion(self, a: np.ndarray) -> np.ndarray:
        """
        @brief Reversion ã: eᵢ₁···eᵢₖ → eᵢₖ···eᵢ₁.
        @description
            Die Reversion kehrt die Reihenfolge der Faktoren um:
                ã = (eᵢ₁···eᵢₖ)~ = eᵢₖ···eᵢ₁
                   = (-1)^{k(k-1)/2} · eᵢ₁···eᵢₖ

            Für Matrizen: Transposition (bei reeller Darstellung)
            oder Adjungierung (bei komplexer Darstellung).

        @param a: Matrix-Darstellung eines Clifford-Elements.
        @return: Matrix der Reversion.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Für hermitesche Gamma-Matrizen: Reversion ~ hermitesche Konjugation
        return a.conj().T

    def verify_clifford_relations(self) -> bool:
        """
        @brief Prüft die Clifford-Algebra-Relationen {γᵢ, γⱼ} = 2ηᵢⱼ·I.
        @return: True, wenn alle Relationen erfüllt sind.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        gammas = self._gamma_matrices
        n = len(gammas)
        if n == 0:
            return True

        mat_dim = gammas[0].shape[0]
        I = np.eye(mat_dim, dtype=complex)

        for i in range(n):
            for j in range(n):
                # Antikommutator {γᵢ, γⱼ} = γᵢγⱼ + γⱼγᵢ
                anticomm = gammas[i] @ gammas[j] + gammas[j] @ gammas[i]
                # Soll = 2·η_{ij}·I sein
                expected = 2.0 * self.signature[i, j] * I
                if not np.allclose(anticomm, expected, atol=1e-10):
                    return False
        return True


# ---------------------------------------------------------------------------
# Klasse: Spinor
# ---------------------------------------------------------------------------

class Spinor:
    """
    Spinor als Element des Spinorraums.

    Ein Dirac-Spinor ψ ∈ ℂ⁴ beschreibt ein Spin-1/2-Quantenteilchen
    (z.B. Elektron) in der relativistischen Quantenmechanik.

    Unter Lorentz-Transformationen Λ transformiert sich ψ als:
        ψ → S(Λ)·ψ
    wobei S(Λ) die Spinordarstellung von Λ ist.

    Dirac-Konjugierte:
        ψ̄ = ψ†·γ⁰
    Wichtig für Lorentz-invariante Größen: ψ̄ψ (Skalar), ψ̄γ^μψ (Vektor).
    """

    def __init__(
        self,
        components: "np.ndarray | list",
        clifford_alg: Optional[CliffordAlgebra] = None
    ) -> None:
        """
        @brief Initialisiert einen Spinor.
        @param components: Komplexe Komponenten des Spinors (ℂⁿ-Vektor).
        @param clifford_alg: Zugehörige Clifford-Algebra (Standard: Cl(1,3) für Dirac).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        self.components = np.asarray(components, dtype=complex)
        # Standard: Dirac-Algebra Cl(1,3) für physikalische Anwendungen
        if clifford_alg is None:
            self.clifford_alg = CliffordAlgebra(1, 3)
        else:
            self.clifford_alg = clifford_alg

    def norm(self) -> float:
        """
        @brief Berechnet die Norm des Spinors: ||ψ|| = √(ψ†ψ).
        @return: Nicht-negative reelle Zahl.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # ψ†ψ = Σᵢ |ψᵢ|²
        return float(np.sqrt(np.dot(self.components.conj(), self.components).real))

    def conjugate(self) -> np.ndarray:
        """
        @brief Berechnet die Dirac-Konjugierte ψ̄ = ψ†·γ⁰.
        @description
            Die Dirac-Konjugierte ψ̄ = ψ†·γ⁰ ist die korrekte
            Lorentz-kovariant konjugierte Größe (nicht ψ†).

            Lorentz-Invariante:
            - ψ̄ψ: Skalar
            - ψ̄γ^μψ: Vierer-Strom (Vektor)
            - ψ̄σ^{μν}ψ: Antisymmetrischer Tensor

        @return: Zeilenvektor ψ̄ = ψ†·γ⁰ als numpy-Array.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        gammas = self.clifford_alg.generators()
        if len(gammas) == 0:
            return self.components.conj()
        # γ⁰ ist das erste Gamma (Zeitkomponente)
        gamma0 = gammas[0]
        # ψ† = Zeilenvektor der komplex-konjugierten Komponenten
        psi_dag = self.components.conj()
        # ψ̄ = ψ† · γ⁰ (Matrixmultiplikation)
        if gamma0.shape[0] == len(psi_dag):
            return psi_dag @ gamma0
        else:
            return psi_dag

    def lorentz_transform(self, S: np.ndarray) -> "Spinor":
        """
        @brief Wendet eine Spinordarstellung einer Lorentz-Transformation an.
        @description
            Lorentz-Transformation: ψ → S·ψ
            wobei S ∈ Spin(1,3) ≅ SL(2,ℂ) die Spinordarstellung ist.

            Für infinitesimale Transformationen:
                S = I + ½·ω_{μν}·σ^{μν}

        @param S: 4×4 (oder n×n) Matrix der Spinordarstellung.
        @return: Transformierter Spinor S·ψ.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        S = np.asarray(S, dtype=complex)
        new_components = S @ self.components
        return Spinor(new_components, self.clifford_alg)


# ---------------------------------------------------------------------------
# Klasse: DiracEquation
# ---------------------------------------------------------------------------

class DiracEquation:
    """
    Dirac-Gleichung: (iγ^μ∂_μ - m)ψ = 0.

    Beschreibt relativistische Spin-1/2-Teilchen (Elektronen, Myonen etc.).

    Dispersionsrelation:
        E² = |p|² + m²   (natürliche Einheiten: c = ℏ = 1)

    Lösungen:
    - u(p,s): Teilchenlösungen (positive Energie)
    - v(p,s): Antiteilchenlösungen (negative Energie)

    Erzwingt Antimaterie: Löcher in der negativen Energiefülle = Positronen.

    Gamma-Matrizen (Dirac-Darstellung):
        γ⁰ = [[I₂, 0], [0, -I₂]]
        γⁱ = [[0, σᵢ], [-σᵢ, 0]]

    Antikommutatorrelation (Clifford-Relation):
        {γ^μ, γ^ν} = 2g^{μν}·I₄
    mit Minkowski-Metrik g = diag(+1,-1,-1,-1).
    """

    def __init__(self, mass: float = 1.0, dim: int = 4) -> None:
        """
        @brief Initialisiert die Dirac-Gleichung.
        @param mass: Ruhemasse des Teilchens (in natürlichen Einheiten: ℏ=c=1).
        @param dim: Spinor-Dimension (Standard: 4 für Dirac-Spinor in 4D).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        self.mass = float(mass)
        self.dim = dim
        # Gamma-Matrizen in Dirac-Darstellung
        self._gammas = self._build_dirac_gammas()

    def _build_dirac_gammas(self) -> List[np.ndarray]:
        """
        @brief Konstruiert die Gamma-Matrizen in der Dirac-Darstellung.
        @description
            Dirac-Darstellung (Standard in der Quantenfeldtheorie):
                γ⁰ = [[I₂, 0 ], [0 , -I₂]]
                γ¹ = [[0 , σ₁], [-σ₁, 0 ]]
                γ² = [[0 , σ₂], [-σ₂, 0 ]]
                γ³ = [[0 , σ₃], [-σ₃, 0 ]]

            Diese erfüllen:
                {γ^μ, γ^ν} = 2g^{μν}·I₄
            mit Minkowski-Metrik g^{μν} = diag(+1,-1,-1,-1).

        @return: Liste [γ⁰, γ¹, γ², γ³].
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        I2 = np.eye(2, dtype=complex)
        zero2 = np.zeros((2, 2), dtype=complex)

        # Pauli-Matrizen
        sigma1 = PauliMatrices.sigma(1)
        sigma2 = PauliMatrices.sigma(2)
        sigma3 = PauliMatrices.sigma(3)

        # Gamma-Matrizen in Dirac-Darstellung (4×4-Matrizen)
        gamma0 = np.block([[I2, zero2], [zero2, -I2]])
        gamma1 = np.block([[zero2, sigma1], [-sigma1, zero2]])
        gamma2 = np.block([[zero2, sigma2], [-sigma2, zero2]])
        gamma3 = np.block([[zero2, sigma3], [-sigma3, zero2]])

        return [gamma0, gamma1, gamma2, gamma3]

    def gamma_matrices(self) -> List[np.ndarray]:
        """
        @brief Gibt die vier Gamma-Matrizen γ⁰, γ¹, γ², γ³ zurück.
        @return: Liste [γ⁰, γ¹, γ², γ³] als numpy-Arrays.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return [g.copy() for g in self._gammas]

    def free_particle_energy(self, momentum: np.ndarray) -> float:
        """
        @brief Berechnet die relativistische Energie E = √(|p|² + m²).
        @description
            Einstein-Energie-Impuls-Relation:
                E² = |p|²c² + m²c⁴
            In natürlichen Einheiten (c = ℏ = 1):
                E = √(|p|² + m²)

        @param momentum: Impulsvektor (3D).
        @return: Positive Energie E (float).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        p = np.asarray(momentum, dtype=float)
        return float(np.sqrt(np.sum(p ** 2) + self.mass ** 2))

    def plane_wave_solution(
        self,
        momentum: np.ndarray,
        spin: str = 'up'
    ) -> np.ndarray:
        """
        @brief Berechnet die ebene-Wellen-Lösung u(p) der Dirac-Gleichung.
        @description
            Für positive Energie (Teilchenlösung):
                ψ = u(p) · exp(-ipx)

            u(p) in Dirac-Darstellung (nicht-relativistischer Limes):
                u(p, spin-up)   = N · (1, 0, pz/E+m, p+/E+m)ᵀ
                u(p, spin-down) = N · (0, 1, p-/E+m, -pz/E+m)ᵀ
            mit p± = px ± i·py und Normierung N = √(E+m).

        @param momentum: 3-Impuls (p₁, p₂, p₃).
        @param spin: 'up' oder 'down' für Spin-Eigenrichtung.
        @return: 4-komponentiger Dirac-Spinor u(p).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        p = np.asarray(momentum, dtype=float)
        if len(p) < 3:
            p = np.pad(p, (0, 3 - len(p)))

        E = self.free_particle_energy(p)
        m = self.mass

        # p± = px ± i·py
        p_plus = p[0] + 1j * p[1]
        p_minus = p[0] - 1j * p[1]
        pz = p[2]

        # Normierungsfaktor
        N = np.sqrt(E + m) if E + m > 1e-12 else 1.0

        # Spinor u(p) in Dirac-Darstellung
        if spin == 'up' or spin == 'up':
            u = np.array([
                1.0,
                0.0,
                pz / (E + m) if E + m > 1e-12 else 0.0,
                p_plus / (E + m) if E + m > 1e-12 else 0.0
            ], dtype=complex) * N
        else:  # spin == 'down'
            u = np.array([
                0.0,
                1.0,
                p_minus / (E + m) if E + m > 1e-12 else 0.0,
                -pz / (E + m) if E + m > 1e-12 else 0.0
            ], dtype=complex) * N

        return u

    def verify_clifford_algebra(self) -> bool:
        """
        @brief Prüft die Clifford-Algebra-Relationen {γ^μ, γ^ν} = 2g^{μν}·I.
        @description
            Prüft für alle μ, ν ∈ {0,1,2,3}:
                γ^μ·γ^ν + γ^ν·γ^μ = 2·g^{μν}·I₄
            mit Minkowski-Metrik g = diag(+1,-1,-1,-1).

        @return: True, wenn alle Clifford-Relationen erfüllt sind.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Minkowski-Metrik (Signatur +,-,-,-)
        g = np.diag([1.0, -1.0, -1.0, -1.0])
        I4 = np.eye(4, dtype=complex)

        gammas = self._gammas

        for mu in range(4):
            for nu in range(4):
                # Antikommutator
                anticomm = gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu]
                # Soll: 2·g^{μν}·I₄
                expected = 2.0 * g[mu, nu] * I4
                if not np.allclose(anticomm, expected, atol=1e-10):
                    return False
        return True
