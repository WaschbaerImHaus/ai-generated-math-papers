"""
@file spinors.py
@brief Spinoren, Clifford-Algebren und Dirac-Gleichung.
@description
    Dieses Modul implementiert die grundlegenden mathematischen Strukturen
    der Spinorgeometrie und der relativistischen Quantenmechanik:

    - Clifford-Algebra Cl(n,m): Erzeugt von {e_i} mit der anti-kommutativen
      Relation {e_i, e_j} = 2η_ij · I (η = Metrik-Tensor)
    - Gamma-Matrizen in der Dirac-Darstellung (4×4 für 4D-Raumzeit)
    - Spinoren: Elemente des Darstellungsraums der Spin-Gruppe Spin(n) ⊂ Cl(n)
    - Dirac-Gleichung im Impulsraum: (γ^μ p_μ - m·I)u(p) = 0
    - Weyl-Spinoren (chirale Projektion)
    - Majorana-Bedingung (Teilchen = eigenes Antiteilchen)
    - Wilson-Fermion-Hamiltonian auf dem 1D-Gitter

    ## Mathematischer Hintergrund

    Die Clifford-Algebra Cl(1,3) über Minkowski-Raumzeit mit Metrik
    η = diag(+1, -1, -1, -1) wird durch 4×4-Gamma-Matrizen realisiert.

    Die Dirac-Gleichung lautet in natürlichen Einheiten (ħ = c = 1):
        (iγ^μ ∂_μ - m)ψ = 0

    Im Impulsraum wird daraus:
        (γ^μ p_μ - m·I)u(p) = 0

    mit der Dispersionsrelation E² = |p|² + m² (on-shell-Bedingung).

@author Michael Fuhrmann
@lastModified 2026-03-10
@version 1.0.0
"""

import numpy as np
from typing import Optional


# ---------------------------------------------------------------------------
# Hilfsfunktionen: Pauli-Matrizen (Bausteine der Gamma-Matrizen)
# ---------------------------------------------------------------------------

def _pauli_matrices() -> list[np.ndarray]:
    """
    @brief Gibt die drei Pauli-Matrizen σ_1, σ_2, σ_3 zurück.
    @description
        Die Pauli-Matrizen sind die Generatoren der SU(2)-Gruppe und dienen
        als Bausteine für die 4×4-Gamma-Matrizen in der Dirac-Darstellung.

        σ_1 = [[0, 1], [1, 0]]
        σ_2 = [[0, -i], [i, 0]]
        σ_3 = [[1, 0], [0, -1]]

    @return list[np.ndarray] Liste mit drei 2×2-Pauli-Matrizen (komplex).
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    # σ_1: Spin-x-Komponente
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    # σ_2: Spin-y-Komponente (enthält imaginäre Einheit)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    # σ_3: Spin-z-Komponente
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [sigma1, sigma2, sigma3]


# ---------------------------------------------------------------------------
# Gamma-Matrizen der Clifford-Algebra
# ---------------------------------------------------------------------------

def gamma_matrices(dim: int = 4) -> list[np.ndarray]:
    """
    @brief Gibt die Gamma-Matrizen der Clifford-Algebra Cl(1,dim-1) zurück.
    @description
        Die Gamma-Matrizen erfüllen die fundamentale Clifford-Relation:
            {γ^μ, γ^ν} = γ^μ γ^ν + γ^ν γ^μ = 2η^{μν} · I

        mit der Minkowski-Metrik η = diag(+1, -1, ..., -1).

        Für dim=2 (2D-Raumzeit):
            γ^0 = σ_1 = [[0, 1], [1, 0]]
            γ^1 = i·σ_2 = [[0, 1], [-1, 0]]

        Für dim=4 (Dirac-Darstellung, 4D-Minkowski-Raumzeit):
            γ^0 = [[I₂, 0], [0, -I₂]]  (zeitartige Matrix, η^{00} = +1)
            γ^k = [[0, σ_k], [-σ_k, 0]] für k=1,2,3 (raumartige, η^{kk} = -1)

        Die γ^5-Matrix (Chiralitäts-Matrix) ist NICHT in dieser Liste enthalten
        — dafür siehe gamma5().

    @param dim int Raumzeit-Dimension. Unterstützt: 2 (Pauli-artig) und 4 (Dirac).
    @return list[np.ndarray] Liste mit 'dim' Gamma-Matrizen.
    @raises ValueError Wenn dim weder 2 noch 4 ist.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    if dim == 2:
        # 2D-Darstellung: Pauli-basierte 2×2-Gamma-Matrizen
        sigma = _pauli_matrices()
        # γ^0 = σ_1
        g0 = sigma[0].copy()
        # γ^1 = i·σ_2  →  {γ^0,γ^1} = 0, {γ^1,γ^1} = -2I
        g1 = 1j * sigma[1]
        return [g0, g1]

    elif dim == 4:
        # Standard Dirac-Darstellung: 4×4-Gamma-Matrizen
        sigma = _pauli_matrices()

        # 2×2-Einheitsmatrix
        I2 = np.eye(2, dtype=complex)
        # 2×2-Nullmatrix
        Z2 = np.zeros((2, 2), dtype=complex)

        # γ^0 = [[I₂, 0], [0, -I₂]]
        # Erzeugt: {γ^0,γ^0} = 2·(γ^0)² = 2·I₄  (positive Metrik-Komponente)
        g0 = np.block([[I2, Z2], [Z2, -I2]])

        # γ^k = [[0, σ_k], [-σ_k, 0]] für k=1,2,3
        # Erzeugt: {γ^k,γ^k} = -2·I₄  (negative Metrik-Komponenten)
        gammas = [g0]
        for sigma_k in sigma:
            g_k = np.block([[Z2, sigma_k], [-sigma_k, Z2]])
            gammas.append(g_k)

        return gammas

    else:
        raise ValueError(
            f"Dimension {dim} wird nicht unterstützt. "
            "Unterstützt: dim=2 (Pauli) oder dim=4 (Dirac)."
        )


def clifford_algebra_check(
    gammas: list[np.ndarray],
    metric: Optional[np.ndarray] = None
) -> dict:
    """
    @brief Überprüft die Clifford-Algebra-Relation {γ^μ, γ^ν} = 2η^{μν}·I.
    @description
        Die Clifford-Algebra-Relation (Anti-Kommutatorrelation) lautet:
            {γ^μ, γ^ν} = γ^μγ^ν + γ^νγ^μ = 2η^{μν}·I

        Für die Minkowski-Metrik η = diag(+1,-1,-1,-1) bedeutet das:
            {γ^0, γ^0} = +2·I
            {γ^k, γ^k} = -2·I für k=1,2,3
            {γ^μ, γ^ν} =  0   für μ≠ν

        Diese Methode berechnet für jedes Paar (μ,ν) den Frobenius-Fehler:
            error_{μν} = ||{γ^μ, γ^ν} - 2η^{μν}·I||_F

    @param gammas list[np.ndarray] Liste der Gamma-Matrizen.
    @param metric Optional[np.ndarray] Metrik-Tensor η (quadratisch).
                  Standard: Minkowski diag(+1,-1,...,-1).
    @return dict Mit Schlüsseln:
            - 'max_error' (float): Maximaler Frobenius-Fehler über alle (μ,ν).
            - 'is_clifford_algebra' (bool): True wenn max_error < 1e-10.
            - 'anticommutators' (dict): {(μ,ν): Frobenius-Fehler}.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    n = len(gammas)
    dim_matrix = gammas[0].shape[0]

    # Standard-Metrik: Minkowski diag(+1,-1,...,-1)
    if metric is None:
        metric = np.diag([1.0] + [-1.0] * (n - 1))

    # Einheitsmatrix mit passender Dimension
    I = np.eye(dim_matrix, dtype=complex)

    # Berechne Anti-Kommutator-Fehler für alle (μ,ν)-Paare
    anticommutators = {}
    max_error = 0.0

    for mu in range(n):
        for nu in range(n):
            # Anti-Kommutator: {γ^μ, γ^ν} = γ^μγ^ν + γ^νγ^μ
            anticomm = gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu]

            # Erwarteter Wert: 2η^{μν}·I
            expected = 2.0 * metric[mu, nu] * I

            # Frobenius-Norm des Fehlers
            error = np.linalg.norm(anticomm - expected, ord='fro')
            anticommutators[(mu, nu)] = float(error)
            max_error = max(max_error, error)

    return {
        'max_error': max_error,
        'is_clifford_algebra': max_error < 1e-10,
        'anticommutators': anticommutators
    }


def gamma5(dim: int = 4) -> np.ndarray:
    """
    @brief Berechnet die γ^5-Matrix (Chiralitäts-Matrix).
    @description
        Die γ^5-Matrix ist definiert als:
            γ^5 = i · γ^0 · γ^1 · γ^2 · γ^3

        Sie besitzt folgende wichtige Eigenschaften:
            (γ^5)² = I              (idempotent bis auf Sign)
            {γ^5, γ^μ} = 0         (anti-kommutiert mit allen Gamma-Matrizen)
            Tr(γ^5) = 0             (spurlos)
            (γ^5)† = γ^5            (hermitesch)

        In der Weyl/Chiral-Basis gilt γ^5 = diag(-I₂, +I₂).
        Sie teilt den Spinorraum in linkschirale (P_L) und rechtschirale (P_R)
        Komponenten:
            P_L = (I - γ^5)/2       (linkshändig)
            P_R = (I + γ^5)/2       (rechtshändig)

    @param dim int Raumzeit-Dimension (nur dim=4 unterstützt).
    @return np.ndarray Die 4×4-Matrix γ^5 (komplex).
    @raises ValueError Wenn dim ≠ 4.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    if dim != 4:
        raise ValueError("gamma5 ist nur für dim=4 definiert.")

    gammas = gamma_matrices(4)
    # γ^5 = i · γ^0 · γ^1 · γ^2 · γ^3
    g5 = 1j * gammas[0] @ gammas[1] @ gammas[2] @ gammas[3]
    return g5


# ---------------------------------------------------------------------------
# Dirac-Spinoren
# ---------------------------------------------------------------------------

def dirac_spinor(
    mass: float,
    momentum: np.ndarray,
    spin: int = 0
) -> np.ndarray:
    """
    @brief Berechnet einen ebene-Wellen-Dirac-Spinor u(p, s).
    @description
        Die positiv-Energie-Lösungen der Dirac-Gleichung lauten:
            u(p, s) = N · [[χ_s], [σ·p/(E+m) · χ_s]]

        mit:
            E = √(m² + |p|²)    (relativistische Energie, on-shell)
            N = √(E + m)         (Normierung, sodass ū·u = 2m)
            χ_0 = [[1], [0]]     (Spin-up, Spinor)
            χ_1 = [[0], [1]]     (Spin-down, Spinor)
            σ·p = p_x·σ₁ + p_y·σ₂ + p_z·σ₃  (Pauli-Skalarprodukt)

        Der Spinor erfüllt die Dirac-Gleichung im Impulsraum:
            (γ^μ p_μ - m·I) u(p, s) = 0
        mit p_μ = (E, -p_x, -p_y, -p_z) (kovarianter 4-Impuls).

    @param mass float Ruhemasse m ≥ 0 (in natürlichen Einheiten).
    @param momentum np.ndarray 3-Impuls [p_x, p_y, p_z].
    @param spin int Spin-Quantenzahl: 0 = Spin-up (↑), 1 = Spin-down (↓).
    @return np.ndarray 4-komponentiger Dirac-Spinor u(p,s) (komplex).
    @raises ValueError Wenn spin nicht 0 oder 1 ist.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    if spin not in (0, 1):
        raise ValueError("spin muss 0 (up) oder 1 (down) sein.")

    momentum = np.asarray(momentum, dtype=complex)
    px, py, pz = momentum[0], momentum[1], momentum[2]

    # Relativistische Energie aus Dispersionsrelation E² = m² + |p|²
    p_sq = float(np.real(px**2 + py**2 + pz**2))
    E = np.sqrt(mass**2 + p_sq)

    # Pauli-Matrizen für σ·p
    sigma = _pauli_matrices()
    # σ·p = p_x·σ₁ + p_y·σ₂ + p_z·σ₃
    sigma_dot_p = px * sigma[0] + py * sigma[1] + pz * sigma[2]

    # Basis-Spinor χ_s (2-komponentig)
    if spin == 0:
        chi = np.array([1.0, 0.0], dtype=complex)  # Spin-up
    else:
        chi = np.array([0.0, 1.0], dtype=complex)  # Spin-down

    # Kleine Komponente: σ·p/(E+m) · χ
    # Verhindere Division durch Null für masselose Grenzfälle mit E≈0
    denom = E + mass
    if abs(denom) < 1e-15:
        small_component = np.zeros(2, dtype=complex)
    else:
        small_component = (sigma_dot_p @ chi) / denom

    # Normierungsfaktor N = √(E+m)
    N = np.sqrt(E + mass)

    # Vollständiger 4-Spinor: u = N · [[χ], [σ·p/(E+m)·χ]]
    spinor = N * np.concatenate([chi, small_component])
    return spinor


def dirac_equation_check(
    psi: np.ndarray,
    mass: float,
    momentum: np.ndarray
) -> dict:
    """
    @brief Überprüft ob ein Spinor ψ die Dirac-Gleichung erfüllt.
    @description
        Im Impulsraum lautet die Dirac-Gleichung:
            (γ^μ p_μ - m·I) u(p) = 0

        mit dem kovarianten 4-Impuls p_μ = (E, -p_x, -p_y, -p_z)
        und der Slash-Notation: p̸ = γ^μ p_μ = γ^0·E - γ^1·p_x - γ^2·p_y - γ^3·p_z.

        Das Residuum ||( p̸ - m·I) ψ|| misst die Verletzung der Dirac-Gleichung.

    @param psi np.ndarray 4-komponentiger Dirac-Spinor ψ.
    @param mass float Ruhemasse m.
    @param momentum np.ndarray 3-Impuls [p_x, p_y, p_z].
    @return dict Mit Schlüsseln:
            - 'residual' (float): ||( p̸ - m·I) ψ||₂ (Euklidische Norm).
            - 'satisfies_dirac' (bool): True wenn residual < 1e-10.
            - 'energy' (float): Relativistische Energie E = √(m²+|p|²).
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    momentum = np.asarray(momentum, dtype=complex)
    px, py, pz = float(np.real(momentum[0])), float(np.real(momentum[1])), float(np.real(momentum[2]))

    # Relativistische Energie
    E = float(np.sqrt(mass**2 + px**2 + py**2 + pz**2))

    # Gamma-Matrizen in Dirac-Darstellung
    gammas = gamma_matrices(4)

    # Feynman-Slash: p̸ = γ^0·E - γ^1·p_x - γ^2·p_y - γ^3·p_z
    # (kovarianter Impuls p_μ = (E, -p_x, -p_y, -p_z))
    p_slash = (E * gammas[0]
               - px * gammas[1]
               - py * gammas[2]
               - pz * gammas[3])

    # Dirac-Operator: (p̸ - m·I)
    I4 = np.eye(4, dtype=complex)
    dirac_op = p_slash - mass * I4

    # Residuum: ||(p̸ - m·I) ψ||₂
    residual = float(np.linalg.norm(dirac_op @ psi))

    return {
        'residual': residual,
        'satisfies_dirac': residual < 1e-10,
        'energy': E
    }


def spin_group_element(theta: float, axis: np.ndarray) -> np.ndarray:
    """
    @brief Berechnet ein Element der Spin(3)-Gruppe SU(2).
    @description
        Die Spin(3)-Gruppe ist die doppelte Überlagerung von SO(3).
        Ein Element entspricht einer Rotation um den Einheitsvektor n̂ um den
        Winkel θ:
            U(θ, n̂) = cos(θ/2)·I + i·sin(θ/2)·(σ·n̂)
                     = cos(θ/2)·I + i·sin(θ/2)·(n_x·σ₁ + n_y·σ₂ + n_z·σ₃)

        Wichtige Eigenschaft (Spinor-Statistik):
            U(2π, n̂) = -I  (Spinor kehrt nach 2π-Rotation sein Vorzeichen um!)
            U(4π, n̂) = +I  (erst nach 4π-Rotation kehrt Spinor in Ausgangszustand)

        Dies unterscheidet Spinoren (halb-ganzzahliger Spin) fundamental von
        gewöhnlichen Vektoren, die sich unter 2π-Rotationen nicht ändern.

    @param theta float Rotationswinkel in Radiant.
    @param axis np.ndarray Rotationsachse [n_x, n_y, n_z] (wird normiert).
    @return np.ndarray 2×2-unitäre Matrix U ∈ SU(2) (komplex).
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    axis = np.asarray(axis, dtype=float)

    # Normiere die Rotationsachse auf Einheitslänge
    norm = np.linalg.norm(axis)
    if norm < 1e-15:
        # Triviale Rotation: identische Matrix
        return np.eye(2, dtype=complex)
    n_hat = axis / norm

    # Pauli-Matrizen
    sigma = _pauli_matrices()

    # σ·n̂ = n_x·σ₁ + n_y·σ₂ + n_z·σ₃
    sigma_dot_n = n_hat[0] * sigma[0] + n_hat[1] * sigma[1] + n_hat[2] * sigma[2]

    # U = cos(θ/2)·I + i·sin(θ/2)·(σ·n̂)
    I2 = np.eye(2, dtype=complex)
    U = np.cos(theta / 2) * I2 + 1j * np.sin(theta / 2) * sigma_dot_n

    return U


def weyl_spinors(psi: np.ndarray) -> dict:
    """
    @brief Zerlegt einen Dirac-Spinor in Weyl-Spinoren (chirale Projektion).
    @description
        Die chiralen Projektionsoperatoren P_L und P_R sind definiert als:
            P_L = (I - γ^5) / 2    (linkshändiger Projektor)
            P_R = (I + γ^5) / 2    (rechtshändiger Projektor)

        Sie erfüllen:
            P_L + P_R = I           (Vollständigkeit: ψ_L + ψ_R = ψ)
            P_L² = P_L              (Projektor-Eigenschaft)
            P_R² = P_R              (Projektor-Eigenschaft)
            P_L · P_R = 0           (Orthogonalität)

        Die Weyl-Spinoren lauten:
            ψ_L = P_L · ψ = (I - γ^5)/2 · ψ   (linkshändig)
            ψ_R = P_R · ψ = (I + γ^5)/2 · ψ   (rechtshändig)

        Im Standardmodell der Teilchenphysik transformieren linkschirale und
        rechtschirale Fermionen unterschiedlich unter schwachen Wechselwirkungen.

    @param psi np.ndarray 4-komponentiger Dirac-Spinor.
    @return dict Mit Schlüsseln:
            - 'left' (np.ndarray): Linkshändiger Weyl-Spinor ψ_L.
            - 'right' (np.ndarray): Rechtshändiger Weyl-Spinor ψ_R.
            - 'is_purely_left' (bool): True wenn ||ψ_R|| < 1e-10.
            - 'is_purely_right' (bool): True wenn ||ψ_L|| < 1e-10.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    psi = np.asarray(psi, dtype=complex)

    # γ^5-Matrix für chirale Projektion
    g5 = gamma5(4)
    I4 = np.eye(4, dtype=complex)

    # Linkshändiger Projektor: P_L = (I - γ^5) / 2
    P_L = (I4 - g5) / 2.0

    # Rechtshändiger Projektor: P_R = (I + γ^5) / 2
    P_R = (I4 + g5) / 2.0

    # Chirale Zerlegung
    psi_L = P_L @ psi
    psi_R = P_R @ psi

    return {
        'left': psi_L,
        'right': psi_R,
        'is_purely_left': float(np.linalg.norm(psi_R)) < 1e-10,
        'is_purely_right': float(np.linalg.norm(psi_L)) < 1e-10
    }


def majorana_condition_check(psi: np.ndarray) -> dict:
    """
    @brief Prüft ob ein Spinor die Majorana-Bedingung erfüllt.
    @description
        Ein Majorana-Fermion ist sein eigenes Antiteilchen. Die Majorana-Bedingung
        lautet:
            ψ = ψ^c = C · ψ̄^T

        wobei:
            ψ̄ = ψ† · γ^0           (Dirac-adjungierter Spinor)
            C = i · γ^2 · γ^0      (Ladungskonjugations-Matrix in Dirac-Darstellung)

        In der Dirac-Darstellung gilt explizit:
            C = i · γ^2 · γ^0

        Das Residuum ||ψ - C·(γ^0·ψ*)^T||₂ = 0 für Majorana-Fermionen.

        Bekannte Majorana-Kandidaten in der Physik:
        - Möglicherweise Neutrinos (experimentell noch ungeklärt)
        - Majorana-Fermionen in topologischen Supraleitern

    @param psi np.ndarray 4-komponentiger Dirac-Spinor.
    @return dict Mit Schlüsseln:
            - 'residual' (float): ||ψ - ψ^c||₂.
            - 'is_majorana' (bool): True wenn residual < 1e-10.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    psi = np.asarray(psi, dtype=complex)

    # Gamma-Matrizen in Dirac-Darstellung
    gammas = gamma_matrices(4)

    # Ladungskonjugations-Matrix: C = i · γ^2 · γ^0
    C = 1j * gammas[2] @ gammas[0]

    # Dirac-adjungierter Spinor: ψ̄ = ψ† · γ^0
    # Transponiert: ψ̄^T = (ψ† · γ^0)^T = (γ^0)^T · (ψ†)^T = (γ^0)^T · ψ*
    # In Dirac-Darstellung: (γ^0)^T = γ^0 (symmetrisch)
    psi_bar_T = gammas[0].T @ np.conj(psi)

    # Ladungskonjugierter Spinor: ψ^c = C · ψ̄^T
    psi_c = C @ psi_bar_T

    # Majorana-Residuum: ||ψ - ψ^c||₂
    residual = float(np.linalg.norm(psi - psi_c))

    return {
        'residual': residual,
        'is_majorana': residual < 1e-10
    }


# ---------------------------------------------------------------------------
# Gitter-Dirac-Operator (Wilson-Fermionen)
# ---------------------------------------------------------------------------

def dirac_hamiltonian_1d(
    n_sites: int = 10,
    mass: float = 1.0,
    lattice_spacing: float = 1.0
) -> np.ndarray:
    """
    @brief Erstellt den Wilson-Fermion-Hamiltonian auf einem 1D-Gitter.
    @description
        Der Wilson-Fermion-Hamiltonian auf einem 1D-Gitter mit periodischen
        Randbedingungen lautet:

            H = Σ_n [m · ψ̄_n ψ_n
                     - (1/2a) · (ψ̄_n γ^1 ψ_{n+1} - ψ̄_{n+1} γ^1 ψ_n)
                     + (r/2a) · (ψ̄_n ψ_{n+1} + ψ̄_{n+1} ψ_n - 2ψ̄_n ψ_n)]

        Der Wilson-Term r (hier r=1) verhindert das "Fermion-Doubling":
        Ohne ihn würden in d Dimensionen 2^d entartete Fermion-Moden auftreten
        (Nielsen-Ninomiya-Theorem), was physikalisch unerwünscht ist.

        In Matrixform hat H die Blockstruktur:
            H = (m + r/a)·I ⊗ I_sites  (Masse + Wilson-Term auf Diagonale)
            - (1/2a)·(γ^1 + r·I) ⊗ T  (Hopping-Term)
            - (1/2a)·(-γ^1 + r·I) ⊗ T^† (konjugierter Hopping-Term)

        wobei T die Translationsmatrix (zyklisch verschoben) ist.

    @param n_sites int Anzahl der Gitterplätze N (Gesamtsystem-Größe).
    @param mass float Ruhemasse m der Gitter-Fermionen.
    @param lattice_spacing float Gitterabstand a.
    @return np.ndarray (2·n_sites × 2·n_sites)-Hermitescher Hamiltonian (komplex).
                       Faktor 2 wegen 2-komponentiger Weyl-Spinoren im 1D-System.
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    # Wir verwenden 2×2-Gamma-Matrizen für das 1D-System (Weyl-Spinoren)
    gammas_2d = gamma_matrices(2)
    g1 = gammas_2d[1]   # γ^1 für räumliche Dimension

    # 2×2-Einheitsmatrix (Spinor-Raum)
    I2 = np.eye(2, dtype=complex)

    # Wilson-Parameter r=1 verhindert Fermion-Doubling
    r = 1.0
    a = lattice_spacing
    t = 1.0 / (2.0 * a)   # Hopping-Parameter (1/2a)

    # Gesamtgröße des Hilbert-Raums: 2 Spinor-Komponenten × n_sites Gitterplätze
    total_dim = 2 * n_sites
    H = np.zeros((total_dim, total_dim), dtype=complex)

    # Masse-Term + Wilson-Term auf der Diagonalen: (m + r/a)·I pro Gitterplatz
    mass_wilson = mass + r / a
    for n in range(n_sites):
        # Füge 2×2-Block auf Diagonale ein
        start = 2 * n
        H[start:start+2, start:start+2] += mass_wilson * I2

    # Hopping-Term mit periodischen Randbedingungen
    for n in range(n_sites):
        # Nächster Gitterplatz (periodisch): (n+1) mod n_sites
        n_plus = (n + 1) % n_sites

        start_n = 2 * n
        start_np = 2 * n_plus

        # Vorwärts-Hopping: -t·(γ^1 - r·I) von n nach n+1
        hop_forward = -t * (g1 - r * I2)
        H[start_n:start_n+2, start_np:start_np+2] += hop_forward

        # Rückwärts-Hopping (hermitesch konjugiert): von n+1 nach n
        H[start_np:start_np+2, start_n:start_n+2] += hop_forward.conj().T

    return H


def dirac_spectrum_1d(
    mass: float = 1.0,
    n_sites: int = 20
) -> dict:
    """
    @brief Berechnet das Spektrum des 1D-Gitter-Dirac-Operators.
    @description
        Das Spektrum des Dirac-Operators auf dem 1D-Gitter wird durch
        Diagonalisierung des Wilson-Fermion-Hamiltonians berechnet.

        Physikalische Interpretationen:
        - Die Massenlücke Δ = min(|Eigenwert|) gibt die minimale Anregungsenergie an.
        - Für freie Fermionen: Δ ≈ m (Ruhemasse dominiert)
        - Das Spektrum ist symmetrisch um Null (particle-hole-Symmetrie)

        Das Wilson-Fermion-Verfahren auf dem Gitter beseitigt spurious Moden
        (Fermion-Doubling), die im naiven Gitter-Ansatz auftreten würden.

    @param mass float Ruhemasse m der Gitter-Fermionen.
    @param n_sites int Anzahl der Gitterplätze.
    @return dict Mit Schlüsseln:
            - 'eigenvalues' (np.ndarray): Alle Eigenwerte (reell für herm. H).
            - 'mass_gap' (float): min(|Eigenwerte|) = Massenlücke.
            - 'is_gapped' (bool): True wenn mass_gap > 1e-10.
            - 'n_zero_modes' (int): Anzahl der Null-Moden (|λ| < 1e-8).
    @lastModified 2026-03-10
    @author Michael Fuhrmann
    """
    # Baue den Hamiltonian
    H = dirac_hamiltonian_1d(n_sites=n_sites, mass=mass)

    # Berechne Eigenwerte (Hermitescher Operator → reelle Eigenwerte)
    eigenvalues = np.linalg.eigvalsh(H)

    # Massenlücke: minimaler Betrag der Eigenwerte
    abs_eigenvalues = np.abs(eigenvalues)
    mass_gap = float(np.min(abs_eigenvalues))

    # Null-Moden: Eigenwerte mit |λ| < 1e-8
    n_zero_modes = int(np.sum(abs_eigenvalues < 1e-8))

    return {
        'eigenvalues': eigenvalues,
        'mass_gap': mass_gap,
        'is_gapped': mass_gap > 1e-10,
        'n_zero_modes': n_zero_modes
    }
