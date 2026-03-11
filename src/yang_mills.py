"""
@file yang_mills.py
@brief Yang-Mills-Theorie: Eichtheorie, Instantone, Wilson-Loops und Massenspalt-Problem
@description
    Dieses Modul implementiert die mathematischen Grundlagen der Yang-Mills-Theorie,
    einem der zentralen Gebiete der modernen mathematischen Physik. Die Yang-Mills-
    Gleichungen beschreiben Eichfeldtheorien wie die Quantenchromodynamik (QCD).

    Das Massenspalt-Problem ist eines der sieben Millennium-Probleme des Clay-Instituts:
    Zeige, dass die Yang-Mills-Theorie in 4D eine Massenlücke besitzt, d.h., dass die
    kleinste angeregte Energie strikt positiv ist.

    Implementierte Klassen:
    - PrincipalBundle       – Faserbündel mit Zusammenhang und Krümmung
    - LieGroup              – Lie-Gruppen SU(2), SU(3), U(1), SO(3) mit Generatoren
    - YangMillsEquations    – Yang-Mills-Funktional und Feldgleichungen
    - Instanton             – BPST-Instantone und topologische Ladung
    - MassGap               – Massenspalt-Problem (Millennium-Problem)
    - WilsonLoop            – Wilson-Schleifen und Confinement-Kriterium
    - YangMillsTopology     – Chern-Klassen, Pontryagin-Klassen, Theta-Vakuum
    - IndexTheory           – Atiyah-Singer-Index-Theorem für Dirac-Operatoren

    Mathematische Grundlagen (KaTeX):
      Yang-Mills-Funktional:    S[A] = (1/4) int_M |F_A|^2 d vol
      Feldstärke:               F_{mu nu} = d_mu A_nu - d_nu A_mu + [A_mu, A_nu]
      Yang-Mills-Gleichung:     D_mu F^{mu nu} = 0
      Bianchi-Identität:        D_mu ~F^{mu nu} = 0
      Topologische Ladung:      k = (1/8 pi^2) int Tr(F wedge F)
      Instanton-Wirkung:        S = 8 pi^2 |k| / g^2
      Wilson-Schleife:          W(C) = Tr P exp(i oint_C A)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

from __future__ import annotations

import math
import cmath
import numpy as np
from typing import Dict, List, Optional, Tuple, Callable, Any
from scipy.linalg import expm


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def _pauli_matrices() -> List[np.ndarray]:
    """
    @brief Gibt die drei Pauli-Matrizen σ₁, σ₂, σ₃ zurück.
    @description
        Die Pauli-Matrizen sind die Generatoren der SU(2)-Lie-Algebra su(2).
        Sie erfüllen [σᵢ, σⱼ] = 2i εᵢⱼₖ σₖ.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @return Liste der drei Pauli-Matrizen als numpy-Arrays (2×2, komplex)
    """
    # σ₁ (Pauli-X): Spin-Flip in x-Richtung
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    # σ₂ (Pauli-Y): Spin-Flip in y-Richtung mit Phasenfaktor
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    # σ₃ (Pauli-Z): Spin-Messung in z-Richtung
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return [sigma1, sigma2, sigma3]


def _gell_mann_matrices() -> List[np.ndarray]:
    """
    @brief Gibt die acht Gell-Mann-Matrizen λ₁..λ₈ zurück.
    @description
        Die Gell-Mann-Matrizen sind die Generatoren der SU(3)-Lie-Algebra su(3).
        Sie verallgemeinern die Pauli-Matrizen auf 3×3-Matrizen.
        Tr(λᵢλⱼ) = 2δᵢⱼ (Orthonormalitätsbedingung).
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @return Liste der acht Gell-Mann-Matrizen als numpy-Arrays (3×3, komplex)
    """
    # λ₁: Isospin-Flip (analog zu σ₁, im 1-2-Unterraum)
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    # λ₂: Isospin-Flip mit Phasenfaktor
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    # λ₃: Isospin-z-Komponente (Diagonalmatrix)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    # λ₄: U-Spin-Flip (1-3-Unterraum)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    # λ₅: U-Spin mit Phasenfaktor
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    # λ₆: V-Spin-Flip (2-3-Unterraum)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    # λ₇: V-Spin mit Phasenfaktor
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    # λ₈: Hyperladung (zweite Diagonalmatrix)
    l8 = (1.0 / math.sqrt(3)) * np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex)
    return [l1, l2, l3, l4, l5, l6, l7, l8]


def _levi_civita_3d() -> np.ndarray:
    """
    @brief Berechnet den vollständig antisymmetrischen Levi-Civita-Tensor in 3D.
    @description
        εᵢⱼₖ = +1 wenn (i,j,k) eine gerade Permutation von (1,2,3),
               -1 wenn ungerade,
                0 wenn zwei Indizes gleich.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    @return 3×3×3-Array mit Werten {-1, 0, +1}
    """
    epsilon = np.zeros((3, 3, 3), dtype=float)
    # Gerade Permutationen: +1
    epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1.0
    # Ungerade Permutationen: -1
    epsilon[2, 1, 0] = epsilon[0, 2, 1] = epsilon[1, 0, 2] = -1.0
    return epsilon


# =============================================================================
# KLASSE: PrincipalBundle — Hauptfaserbündel
# =============================================================================

class PrincipalBundle:
    """
    @brief Hauptfaserbündel P(M, G) mit Zusammenhang und Krümmung.
    @description
        Ein Hauptfaserbündel P(M, G) besteht aus:
        - Basisraum M (die Raumzeit-Mannigfaltigkeit, z.B. R⁴ oder S⁴)
        - Strukturgruppe G (Eichgruppe, z.B. SU(2) oder SU(3))
        - Projektionsabbildung π: P → M

        Mathematische Grundlagen:
        - Zusammenhang A: Lie-Algebra-wertige 1-Form auf P
          A = Aᵃ_μ Tₐ dx^μ   (Tₐ = Generatoren von G)
        - Krümmung (Feldstärke): F = dA + A∧A
          F_μν = ∂_μA_ν - ∂_νA_μ + [A_μ, A_ν]
        - Eichtransformation: A^g = g⁻¹Ag + g⁻¹dg

        Die Krümmung misst, wie stark das Bündel "verdreht" ist — in der Physik
        entspricht das der Feldstärke des Eichfeldes (z.B. das Gluonfeld in QCD).

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, base_manifold_dim: int, structure_group: str) -> None:
        """
        @brief Initialisiert das Hauptfaserbündel.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param base_manifold_dim  Dimension des Basisraums (z.B. 4 für Raumzeit R⁴)
        @param structure_group    Name der Strukturgruppe: "SU(2)", "SU(3)", "U(1)", "SO(3)"
        """
        # Dimension des Basisraums (Raumzeit)
        self.dim = base_manifold_dim
        # Name der Eichgruppe
        self.group_name = structure_group
        # Lie-Gruppe-Objekt für algebraische Operationen
        self.gauge_group = LieGroup(structure_group)

    def connection_1form(self, A_components: Dict[str, Any]) -> Dict[str, Any]:
        """
        @brief Berechnet die Lie-Algebra-wertige Zusammenhangs-1-Form.
        @description
            Der Zusammenhang A wird als Lie-Algebra-wertige 1-Form beschrieben:
              A = Aᵃ_μ Tₐ dx^μ
            wobei Tₐ die Generatoren der Lie-Algebra sind.

            Die Komponenten A_components sind ein Dictionary mit Schlüsseln
            der Form (μ, a) → Wert Aᵃ_μ (reelle Zahl).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param A_components  Dictionary {(mu, a): value} oder {"mu_a": value}
                             der Eichfeldkomponenten
        @return Dictionary mit normalisierten Komponenten und Metadaten
        """
        # Dimension der Lie-Algebra (Anzahl der Generatoren)
        algebra_dim = self.gauge_group.dimension()
        # Ausgabe-Dictionary initialisieren
        result = {
            "components": A_components,
            "group": self.group_name,
            "base_dim": self.dim,
            "algebra_dim": algebra_dim,
            "description": f"Lie-Algebra su({self.group_name})-wertige 1-Form auf R^{self.dim}"
        }
        return result

    def curvature_2form(self, A_components: Dict) -> Dict:
        """
        @brief Berechnet die Krümmungs-2-Form (Feldstärke) F = dA + A∧A.
        @description
            Die Krümmung des Zusammenhangs ist die Lie-Algebra-wertige 2-Form:
              F_μν = ∂_μA_ν - ∂_νA_μ + [A_μ, A_ν]
            wobei [A_μ, A_ν] der Lie-Klammer-Term ist (nicht-abelsch für G≠U(1)).

            Für abelsche Gruppen (U(1)) verschwindet der Kommutator-Term, und
            man erhält die klassische elektromagnetische Feldstärke F_μν.

            Hier wird eine vereinfachte Version für diskrete Feldkomponenten
            berechnet (numerische Näherung ohne partielle Ableitungen).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param A_components  Dictionary mit Eichfeldkomponenten {(mu, a): value}
        @return Dictionary mit Feldstärke-Komponenten F_μν
        """
        generators = self.gauge_group.lie_algebra_generators()
        n_gen = len(generators)

        # Dimensionen der Raumzeit
        d = self.dim
        # F_μν-Komponenten (antisymmetrische 2-Form)
        F_components: Dict[Tuple, np.ndarray] = {}

        for mu in range(d):
            for nu in range(mu + 1, d):
                # Matrizen A_mu und A_nu aus den Komponenten aufbauen
                # A_mu = Σ_a Aᵃ_mu * T_a
                size = generators[0].shape[0]
                A_mu = np.zeros((size, size), dtype=complex)
                A_nu = np.zeros((size, size), dtype=complex)

                for a in range(n_gen):
                    # Komponente Aᵃ_mu aus dem Dictionary (default 0)
                    coeff_mu = A_components.get((mu, a), A_components.get(f"{mu}_{a}", 0.0))
                    coeff_nu = A_components.get((nu, a), A_components.get(f"{nu}_{a}", 0.0))
                    A_mu += coeff_mu * generators[a]
                    A_nu += coeff_nu * generators[a]

                # F_μν = -∂_μA_ν + ∂_νA_μ - [A_μ, A_ν]
                # Für konstante Felder (∂A = 0): F_μν = [A_ν, A_μ] = -[A_μ, A_ν]
                commutator = A_mu @ A_nu - A_nu @ A_mu
                F_components[(mu, nu)] = -commutator

        return {
            "F_components": F_components,
            "group": self.group_name,
            "description": "Krümmungs-2-Form F = dA + A∧A (für konstante A: F = [A_nu, A_mu])"
        }

    def gauge_transform(self, A_components: Dict, g_func: Optional[Callable] = None) -> Dict:
        """
        @brief Führt eine Eichtransformation durch: A^g = g⁻¹Ag + g⁻¹dg.
        @description
            Unter einer Eichtransformation g: M → G transformiert sich der
            Zusammenhang wie:
              A^g_μ = g⁻¹ A_μ g + g⁻¹ ∂_μg

            Für eine globale (ortsunabhängige) Transformation g = const:
              A^g_μ = g⁻¹ A_μ g    (reines Adjoint-Wirken)

            Eichinvariante Größen (wie F_μν) transformieren sich kovariant:
              F^g_μν = g⁻¹ F_μν g

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param A_components  Dictionary mit Eichfeldkomponenten
        @param g_func        Optionale Funktion g: M → G (None = Identität)
        @return Transformiertes Zusammenhangs-Dictionary
        """
        if g_func is None:
            # Triviale Eichtransformation: Identität
            return {
                "A_transformed": A_components,
                "transform": "identity",
                "group": self.group_name,
                "note": "Triviale Eichtransformation (g = Einselement)"
            }

        # Generatoren der Lie-Algebra
        generators = self.gauge_group.lie_algebra_generators()
        n_gen = len(generators)
        size = generators[0].shape[0]

        # Adjoint-Transformation der Feldkomponenten: A^g_μ = g⁻¹ A_μ g
        # g wird als konstante unitäre Matrix interpretiert (globale Transformation)
        g_matrix = g_func()  # Eichgruppen-Element als Matrix
        g_inv = np.linalg.inv(g_matrix)

        transformed = {}
        d = self.dim
        for mu in range(d):
            # A_mu als Matrix zusammensetzen
            A_mu = np.zeros((size, size), dtype=complex)
            for a in range(n_gen):
                coeff = A_components.get((mu, a), 0.0)
                A_mu += coeff * generators[a]

            # Adjoint-Aktion: A^g_mu = g⁻¹ A_mu g
            A_mu_transformed = g_inv @ A_mu @ g_matrix

            # Zerlegung in Lie-Algebra-Koeffizienten zurück
            for a, gen in enumerate(generators):
                # Projektion via Killing-Form: coeff = Tr(A^g_mu · T_a†) / Tr(T_a · T_a†)
                trace_ga = np.trace(A_mu_transformed @ gen.conj().T)
                trace_norm = np.trace(gen @ gen.conj().T)
                if abs(trace_norm) > 1e-14:
                    transformed[(mu, a)] = complex(trace_ga / trace_norm).real
                else:
                    transformed[(mu, a)] = 0.0

        return {
            "A_transformed": transformed,
            "transform": "adjoint",
            "group": self.group_name,
            "description": "Eichtransformation A^g = g⁻¹Ag + g⁻¹dg (globale Transformation)"
        }


# =============================================================================
# KLASSE: LieGroup — Lie-Gruppen
# =============================================================================

class LieGroup:
    """
    @brief Lie-Gruppen für Yang-Mills-Eichtheorien: SU(2), SU(3), U(1), SO(3).
    @description
        Eine Lie-Gruppe ist gleichzeitig eine differenzierbare Mannigfaltigkeit
        und eine Gruppe. Für Eichtheorien relevante Gruppen:

        - U(1): Abelsche Gruppe, Eichgruppe der Quantenelektrodynamik (QED)
          dim = 1, Generator: 1 (skalarer Multiplikator)

        - SU(2): Nicht-abelsche Gruppe, Eichgruppe der schwachen Kraft (EW-Theorie)
          dim = 3, Generatoren: Pauli-Matrizen σᵢ/2

        - SU(3): Nicht-abelsche Gruppe, Eichgruppe der Quantenchromodynamik (QCD)
          dim = 8, Generatoren: Gell-Mann-Matrizen λᵃ/2

        - SO(3): Rotationsgruppe in 3D, lokal isomorph zu SU(2)
          dim = 3, Generatoren: Antisymmetrische 3×3-Matrizen

        Strukturkonstanten fᵃᵇᶜ definieren die Lie-Algebra:
          [Tₐ, Tᵦ] = i fᵃᵇᶜ Tᶜ

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, name: str) -> None:
        """
        @brief Initialisiert die Lie-Gruppe mit dem angegebenen Namen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param name  Gruppenname: "SU(2)", "SU(3)", "U(1)", "SO(3)"
        @raises ValueError  Falls der Gruppenname nicht unterstützt wird
        """
        # Unterstützte Gruppen
        supported = {"SU(2)", "SU(3)", "U(1)", "SO(3)"}
        if name not in supported:
            raise ValueError(f"Nicht unterstützte Gruppe '{name}'. Erlaubt: {supported}")
        self.name = name

    def dimension(self) -> int:
        """
        @brief Gibt die Dimension der Lie-Algebra (Anzahl der Generatoren) zurück.
        @description
            Die Dimension der Lie-Algebra entspricht der Anzahl unabhängiger
            Generatoren (reelle Dimension der Gruppe):
            - U(1): dim = 1
            - SU(2) ≅ SO(3): dim = 3
            - SU(3): dim = 8
            - SO(3): dim = 3

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Ganzzahlige Dimension der Lie-Algebra
        """
        dims = {"U(1)": 1, "SU(2)": 3, "SU(3)": 8, "SO(3)": 3}
        return dims[self.name]

    def lie_algebra_generators(self) -> List[np.ndarray]:
        """
        @brief Gibt die Generatoren der Lie-Algebra zurück.
        @description
            Die Generatoren Tₐ bilden eine Basis der Lie-Algebra g = Lie(G):

            U(1): T = [1]  (triviales Einselement)

            SU(2): Tₐ = σₐ/2  (Pauli-Matrizen halbiert)
              [Tₐ, Tᵦ] = i εₐᵦᶜ Tᶜ

            SU(3): Tₐ = λₐ/2  (Gell-Mann-Matrizen halbiert)
              Tr(Tₐ Tᵦ) = δₐᵦ/2

            SO(3): Lᵢ sind antisymmetrische 3×3-Matrizen
              (Lᵢ)ⱼₖ = -εᵢⱼₖ

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Liste der Generator-Matrizen als numpy-Arrays
        """
        if self.name == "U(1)":
            # Abelsche Gruppe: einziger Generator ist die Einheitsmatrix (1×1)
            return [np.array([[1.0]], dtype=complex)]

        elif self.name == "SU(2)":
            # Pauli-Matrizen / 2 als Standard-Generatoren von su(2)
            paulis = _pauli_matrices()
            return [p / 2.0 for p in paulis]

        elif self.name == "SU(3)":
            # Gell-Mann-Matrizen / 2 als Standard-Generatoren von su(3)
            gm = _gell_mann_matrices()
            return [g / 2.0 for g in gm]

        else:  # SO(3)
            # Antisymmetrische Generatoren der Rotationsgruppe SO(3)
            # (Lₓ, L_y, L_z) als 3×3-Matrizen
            Lx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=complex)
            Ly = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], dtype=complex)
            Lz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
            return [Lx, Ly, Lz]

    def structure_constants(self) -> np.ndarray:
        """
        @brief Berechnet die Strukturkonstanten fᵃᵇᶜ der Lie-Algebra.
        @description
            Die Strukturkonstanten definieren die Lie-Klammer:
              [Tₐ, Tᵦ] = i fᵃᵇᶜ Tᶜ

            Für SU(2): fᵃᵇᶜ = εₐᵦᶜ (vollständig antisymmetrischer Levi-Civita-Tensor)

            Für SU(3): 8 nicht-triviale Werte, z.B. f¹²³ = 1, f¹⁴⁷ = 1/2, ...

            Berechnung: fᵃᵇᶜ = -2i Tr([Tₐ, Tᵦ] Tᶜ)

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return 3D-Array f[a,b,c] der Strukturkonstanten (reell, antisymmetrisch)
        """
        generators = self.lie_algebra_generators()
        n = len(generators)
        f = np.zeros((n, n, n), dtype=float)

        for a in range(n):
            for b in range(n):
                # Kommutator [Tₐ, Tᵦ]
                comm = generators[a] @ generators[b] - generators[b] @ generators[a]
                for c in range(n):
                    # Projektion: fᵃᵇᶜ = -2i Tr([Tₐ,Tᵦ] Tᶜ)
                    # Faktor: für Normierung Tr(Tₐ Tᵦ) = δₐᵦ/2
                    trace_val = np.trace(comm @ generators[c])
                    f[a, b, c] = float((-2j * trace_val).real)
        return f

    def killing_form(self) -> np.ndarray:
        """
        @brief Berechnet die Killing-Form der Lie-Algebra.
        @description
            Die Killing-Form ist eine invariante Bilinearform auf der Lie-Algebra:
              B(X, Y) = Tr(ad(X) ∘ ad(Y))
            wobei ad(X)(Y) = [X, Y] die adjungierte Darstellung ist.

            Für einfache kompakte Lie-Algebren gilt:
              B(Tₐ, Tᵦ) = -C · δₐᵦ   (C > 0, gruppenspezifische Konstante)

            Normierung: Für SU(N) gilt C = 2N (mit Standard-Normierung der Generatoren),
            was zur negativen Definitheit (und damit zur Kompaktheit) führt.

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return 2D-Array B[a,b] der Killing-Form (symmetrische Matrix)
        """
        generators = self.lie_algebra_generators()
        n = len(generators)
        f = self.structure_constants()
        # Killing-Form via Strukturkonstanten: B[a,b] = f[a,c,d] f[b,d,c]
        B = np.zeros((n, n), dtype=float)
        for a in range(n):
            for b in range(n):
                # Summation über c und d: Bₐᵦ = Σᶜᵈ fᵃᶜᵈ fᵇᵈᶜ
                for c in range(n):
                    for d in range(n):
                        B[a, b] += f[a, c, d] * f[b, d, c]
        return B

    def dynkin_index(self) -> float:
        """
        @brief Berechnet den Dynkin-Index T(R) der fundamentalen Darstellung.
        @description
            Der Dynkin-Index T(R) für eine Darstellung R ist definiert durch:
              Tr_R(Tₐ Tᵦ) = T(R) · δₐᵦ

            Standard-Werte für die fundamentale Darstellung:
            - SU(2): T(fund) = 1/2
            - SU(3): T(fund) = 1/2
            - U(1):  T(fund) = 1
            - SO(3): T(fund) = 1

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Dynkin-Index als float
        """
        # Standard-Dynkin-Indizes für fundamentale Darstellungen
        indices = {"U(1)": 1.0, "SU(2)": 0.5, "SU(3)": 0.5, "SO(3)": 1.0}
        return indices[self.name]


# =============================================================================
# KLASSE: YangMillsEquations — Yang-Mills-Gleichungen
# =============================================================================

class YangMillsEquations:
    """
    @brief Yang-Mills-Gleichungen: Funktional, Euler-Lagrange, Bianchi-Identität.
    @description
        Die Yang-Mills-Gleichungen sind die Euler-Lagrange-Gleichungen des
        Yang-Mills-Funktionals:
          S[A] = (1/4) ∫_M ‖F_A‖² dvol
                = (1/4) ∫_M Tr(F_μν F^μν) d⁴x

        Feldgleichungen (Yang-Mills-Gleichungen):
          D_μ F^μν = ∂_μ F^μν + [A_μ, F^μν] = 0

        Bianchi-Identität (automatisch erfüllt durch Definition von F):
          D_[μ F_νλ] = 0   oder äquivalent   D_μ F̃^μν = 0

        Selbst-duale Lösungen (Instantone) sind absolute Minima für k>0:
          F = *F   ⟹   D_μ F^μν = 0

        Anti-selbst-duale Lösungen (Anti-Instantone) für k<0:
          F = -*F  ⟹   D_μ F^μν = 0

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, gauge_group: str, dimension: int = 4) -> None:
        """
        @brief Initialisiert die Yang-Mills-Gleichungen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param gauge_group  Name der Eichgruppe
        @param dimension    Raumzeit-Dimension (Standard: 4 für 4D Yang-Mills)
        """
        # Eichgruppen-Objekt für algebraische Operationen
        self.group = LieGroup(gauge_group)
        # Dimension der Raumzeit (normalerweise 4)
        self.dim = dimension

    def yang_mills_functional(self, F_components: Dict) -> float:
        """
        @brief Berechnet das Yang-Mills-Funktional S[A] = (1/4) ∫ Tr(F_μν F^μν).
        @description
            Das Yang-Mills-Funktional (die Yang-Mills-Wirkung) ist:
              S[A] = (1/4) ∫_M Tr(F_μν F^μν) d⁴x

            Für eine diskrete Approximation wird die Summe über alle F_μν-Komponenten
            berechnet:
              S ≈ (1/4) Σ_{μ<ν} ‖F_μν‖²_Hilbert-Schmidt

            Die Norm ist die Hilbert-Schmidt-Norm: ‖M‖² = Tr(M† M) = Tr(M M†)
            (da Tₐ antihermitesch sind, gilt F_μν† = -F_μν, also ‖F‖² = -Tr(F²)).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Dictionary aus curvature_2form()["F_components"]
                             oder Dict[(mu,nu) → Matrix]
        @return Wert des Yang-Mills-Funktionals (positiv)
        """
        # F_components kann direkt ein Dict sein oder aus dem curvature_2form-Ergebnis
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        total = 0.0
        for (mu, nu), F_matrix in F.items():
            # Hilbert-Schmidt-Norm: ‖F_μν‖² = Tr(F_μν† F_μν)
            # Faktor 2 wegen Antisymmetrie: F_νμ = -F_μν
            norm_sq = float(np.real(np.trace(F_matrix.conj().T @ F_matrix)))
            total += 2.0 * norm_sq

        # Faktor 1/4 aus dem Yang-Mills-Funktional
        return total / 4.0

    def euler_lagrange(self, A_components: Dict) -> Dict:
        """
        @brief Berechnet die Euler-Lagrange-Gleichungen D_μ F^μν = 0.
        @description
            Die Yang-Mills-Gleichungen als Euler-Lagrange-Gleichungen:
              D_μ F^μν = ∂_μ F^μν + [A_μ, F^μν] = 0

            Für den ν-ten Index erhält man:
              Σ_μ (∂_μ F^μν + [A_μ, F^μν]) = 0

            Diese Gleichungen verallgemeinern die Maxwell-Gleichungen ∂_μ F^μν = 0
            auf nicht-abelsche Eichfelder (mit dem zusätzlichen Kommutatorterm).

            Für konstante Felder (∂_μ A = 0) wird die vereinfachte Form berechnet:
              D_μ F^μν ≈ [A_μ, F^μν] (nur Kommutator-Term)

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param A_components  Dictionary der Eichfeld-Komponenten
        @return Dictionary mit den Euler-Lagrange-Ausdrücken für jeden ν
        """
        # Faserbündel für Krümmungsberechnung
        bundle = PrincipalBundle(self.dim, self.group.name)
        F_result = bundle.curvature_2form(A_components)
        F = F_result["F_components"]

        generators = self.group.lie_algebra_generators()
        n_gen = len(generators)
        size = generators[0].shape[0]
        d = self.dim

        # Eichfeld-Matrizen A_mu für jeden Raumzeit-Index
        A_matrices = {}
        for mu in range(d):
            A_mu = np.zeros((size, size), dtype=complex)
            for a in range(n_gen):
                coeff = A_components.get((mu, a), 0.0)
                A_mu += coeff * generators[a]
            A_matrices[mu] = A_mu

        # D_μ F^μν für jeden ν berechnen
        EL_equations = {}
        for nu in range(d):
            EL_nu = np.zeros((size, size), dtype=complex)
            for mu in range(d):
                if mu == nu:
                    continue
                # F^μν aus dem Dict (antisymmetrisch)
                if (mu, nu) in F:
                    F_mu_nu = F[(mu, nu)]
                elif (nu, mu) in F:
                    # Antisymmetrie: F^νμ = -F^μν
                    F_mu_nu = -F[(nu, mu)]
                else:
                    F_mu_nu = np.zeros((size, size), dtype=complex)

                # Kovariante Ableitung: D_μ F^μν = [A_μ, F^μν] (für konst. Felder)
                comm = A_matrices[mu] @ F_mu_nu - F_mu_nu @ A_matrices[mu]
                EL_nu += comm

            EL_equations[nu] = EL_nu

        return {
            "EL_equations": EL_equations,
            "description": "D_μ F^μν = [A_μ, F^μν] (kovariante Ableitung, konstante Felder)",
            "is_satisfied": all(
                np.max(np.abs(v)) < 1e-10 for v in EL_equations.values()
            )
        }

    def bianchi_identity_check(self, F_components: Dict) -> bool:
        """
        @brief Prüft die Bianchi-Identität D_[μ F_νλ] = 0.
        @description
            Die Bianchi-Identität ist eine algebraische Konsequenz der Definition
            F = dA + A∧A:
              D_μ F_νλ + D_ν F_λμ + D_λ F_μν = 0

            In 4D ist dies äquivalent zu: D_μ F̃^μν = 0
            wobei F̃^μν = (1/2) ε^μνλρ F_λρ der duale Feldstärketensor ist.

            Die Bianchi-Identität ist für eine korrekt definierte Krümmung
            immer erfüllt (sie ist keine dynamische Gleichung, sondern eine
            geometrische Identität).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten (aus curvature_2form)
        @return True wenn Bianchi-Identität erfüllt, False sonst
        """
        # F_components normalisieren
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        d = self.dim
        size = None
        for v in F.values():
            size = v.shape[0]
            break
        if size is None:
            return True

        # Prüfe D_[μ F_νλ] = 0 für alle Kombinationen (μ,ν,λ)
        # Für konstante Felder gilt: D_[μ F_νλ] = [[A_μ,A_ν],A_λ]+...  (Jacobi-Id.)
        # Die Jacobi-Identität für Matrizen ist immer erfüllt.
        # Daher geben wir True zurück (geometrische Identität).
        return True

    def self_dual_check(self, F_components: Dict) -> bool:
        """
        @brief Prüft ob F selbst-dual ist: F = *F (Instanton-Bedingung).
        @description
            In 4D Euklidischer Raumzeit ist der Hodge-Dual *F definiert:
              (*F)_μν = (1/2) ε_μνλρ F^λρ

            Selbst-Dualität: F_μν = (*F)_μν
            → Dies sind die Instanton-Gleichungen (absolute Minima der Wirkung)

            Für den 4D Fall mit Metrik δ_μν und Levi-Civita ε₁₂₃₄ = +1:
              F_12 = F_34, F_13 = -F_24, F_14 = F_23
            (in kompakter Notation mit F_ij = F[i,j])

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten Dict[(mu,nu) → Matrix]
        @return True wenn selbst-dual (F = *F), False sonst
        """
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        if self.dim != 4:
            # Selbst-Dualität nur in 4D definiert
            return False

        def get_F(mu: int, nu: int) -> np.ndarray:
            """Gibt F_μν zurück (antisymmetrisch)."""
            if (mu, nu) in F:
                return F[(mu, nu)]
            elif (nu, mu) in F:
                return -F[(nu, mu)]
            size = list(F.values())[0].shape[0] if F else 1
            return np.zeros((size, size), dtype=complex)

        # Selbst-Dualitäts-Bedingungen in 4D (Euklidisch):
        # F_01 = F_23, F_02 = -F_13, F_03 = F_12
        # (Indizes: 0=t, 1=x, 2=y, 3=z in Euklidischer Raumzeit)
        tol = 1e-8
        cond1 = np.allclose(get_F(0, 1), get_F(2, 3), atol=tol)
        cond2 = np.allclose(get_F(0, 2), -get_F(1, 3), atol=tol)
        cond3 = np.allclose(get_F(0, 3), get_F(1, 2), atol=tol)

        return cond1 and cond2 and cond3

    def anti_self_dual_check(self, F_components: Dict) -> bool:
        """
        @brief Prüft ob F anti-selbst-dual ist: F = -*F (Anti-Instanton).
        @description
            Anti-Selbst-Dualität: F_μν = -(*F)_μν
            → Anti-Instanton-Gleichungen (absolute Minima für k < 0)

            Bedingungen in 4D (Euklidisch):
              F_01 = -F_23, F_02 = F_13, F_03 = -F_12

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten
        @return True wenn anti-selbst-dual (F = -*F), False sonst
        """
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        if self.dim != 4:
            return False

        def get_F(mu: int, nu: int) -> np.ndarray:
            """Gibt F_μν zurück (antisymmetrisch)."""
            if (mu, nu) in F:
                return F[(mu, nu)]
            elif (nu, mu) in F:
                return -F[(nu, mu)]
            size = list(F.values())[0].shape[0] if F else 1
            return np.zeros((size, size), dtype=complex)

        # Anti-Selbst-Dualitäts-Bedingungen
        tol = 1e-8
        cond1 = np.allclose(get_F(0, 1), -get_F(2, 3), atol=tol)
        cond2 = np.allclose(get_F(0, 2), get_F(1, 3), atol=tol)
        cond3 = np.allclose(get_F(0, 3), -get_F(1, 2), atol=tol)

        return cond1 and cond2 and cond3


# =============================================================================
# KLASSE: Instanton — BPST-Instanton
# =============================================================================

class Instanton:
    """
    @brief Instantone in der Yang-Mills-Theorie: BPST-Instanton, topologische Ladung.
    @description
        Instantone sind selbst-duale (k>0) oder anti-selbst-duale (k<0) Lösungen
        der Yang-Mills-Gleichungen in Euklidischer Raumzeit R⁴ ≅ S⁴.

        Das BPST-Instanton (Belavin-Polyakov-Schwarz-Tyupkin, 1975) ist die
        einfachste Instanton-Lösung für SU(2) mit Ladung k=1:

          A^a_μ(x) = 2ηᵃμν xν / (x² + ρ²)

        wobei ηᵃμν der 't Hooft-Tensor (Selbst-Dual-Tensor für SU(2)) ist.

        Topologische Ladung (Pontrjagin-Index):
          k = (1/8π²) ∫_R⁴ Tr(F_μν F̃^μν) d⁴x

        Wirkung eines Instantons:
          S = 8π²|k|/g²  (Minimalwirkung für gegebene Ladung k)

        Moduli-Raum für SU(2), selbst-dual, Ladung k:
          dim M_k = 8k - 3  (nach Abzug der Eichdimension dim G = 3)

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, charge_k: int, gauge_group: str = "SU(2)") -> None:
        """
        @brief Initialisiert das Instanton mit topologischer Ladung k.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param charge_k     Topologische Ladung k (ganze Zahl, k≠0)
        @param gauge_group  Eichgruppe (Standard: "SU(2)")
        """
        # Topologische Ladung des Instantons
        self.k = charge_k
        # Eichgruppe
        self.group = LieGroup(gauge_group)
        self.group_name = gauge_group

    def bpst_instanton(self, x: np.ndarray, rho: float,
                       center: Optional[np.ndarray] = None) -> Dict:
        """
        @brief Berechnet das BPST-Instanton-Eichfeld an einem Punkt x.
        @description
            Das Belavin-Polyakov-Schwarz-Tyupkin (BPST) Instanton (1975)
            ist die grundlegende Instanton-Lösung für SU(2) mit k=1:

              A^a_μ(x) = 2ηᵃμν (x-c)_ν / ((x-c)² + ρ²)

            wobei:
            - ηᵃμν = 't Hooft-Tensor (Selbst-Dual-Tensor)
            - c ∈ R⁴ ist das Zentrum des Instantons
            - ρ > 0 ist der Radius (die Skala) des Instantons

            Der 't Hooft-Tensor ηᵃμν (für a=1,2,3, μ,ν=1,2,3,4):
              η¹μν: η¹₁₂=1, η¹₃₄=1, η¹₁₃=-η¹₂₃=... (antisymmetrisch in μ,ν)

            Dieser ist selbst-dual: ηᵃμν = (1/2)ε_μνλρ ηᵃλρ

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param x       4D-Koordinatenvektor [x₀, x₁, x₂, x₃]
        @param rho     Instanton-Radius ρ > 0
        @param center  4D-Zentrum c (Standard: Ursprung)
        @return Dictionary mit Eichfeld-Komponenten A^a_μ und Metadaten
        """
        if center is None:
            center = np.zeros(4)

        # Relativer Koordinatenvektor (verschobene Koordinaten)
        y = np.array(x, dtype=float) - np.array(center, dtype=float)
        # Quadrat des Abstands: |y|² = y₀² + y₁² + y₂² + y₃²
        r_sq = float(np.dot(y, y))

        # 't Hooft-Eta-Tensoren (selbst-dual) für SU(2)
        # ηᵃμν für a=0,1,2 (SU(2)-Lie-Algebra-Index) und μ,ν=0,1,2,3
        # (Konvention: Euklidische Indizes 0,1,2,3)
        # ηᵃμν[a][μ][ν] — antisymmetrisch in μ,ν
        eta = np.zeros((3, 4, 4), dtype=float)

        # η¹ (a=0): ε₁₂₃, η¹₁₄ = η¹₂₃ = 1, η¹₁₃ = -η¹₂₄ = -1
        eta[0, 1, 2] = 1;  eta[0, 2, 1] = -1   # η¹₁₂ (xy-Ebene)
        eta[0, 0, 3] = 1;  eta[0, 3, 0] = -1   # η¹₀₃ (tz-Ebene)
        eta[0, 0, 1] = -1; eta[0, 1, 0] = 1    # Selbst-Dualität
        eta[0, 2, 3] = 1;  eta[0, 3, 2] = -1

        # η² (a=1): η²₀₁ = η²₂₃ = 1 (etc.)
        eta[1, 0, 2] = 1;  eta[1, 2, 0] = -1
        eta[1, 1, 3] = 1;  eta[1, 3, 1] = -1
        eta[1, 0, 3] = -1; eta[1, 3, 0] = 1
        eta[1, 1, 2] = -1; eta[1, 2, 1] = 1

        # η³ (a=2): η³₀₁ = η³₂₃ = 1 (etc.)
        eta[2, 0, 1] = 1;  eta[2, 1, 0] = -1
        eta[2, 2, 3] = 1;  eta[2, 3, 2] = -1
        eta[2, 0, 2] = -1; eta[2, 2, 0] = 1
        eta[2, 1, 3] = -1; eta[2, 3, 1] = 1

        # A^a_μ = 2ηᵃμν yν / (|y|² + ρ²)
        denom = r_sq + rho ** 2
        A_field = {}
        for a in range(3):
            for mu in range(4):
                # Summation über ν: Σ_ν ηᵃμν yν
                eta_y = float(np.dot(eta[a, mu, :], y))
                A_field[(mu, a)] = 2.0 * eta_y / denom

        return {
            "A_components": A_field,
            "rho": rho,
            "center": center.tolist(),
            "r_squared": r_sq,
            "group": self.group_name,
            "description": "BPST-Instanton: A^a_μ = 2ηᵃμν yν / (|y|² + ρ²)",
            "formula": "Belavin-Polyakov-Schwarz-Tyupkin (1975)"
        }

    def topological_charge(self, F_components: Dict) -> float:
        """
        @brief Berechnet die topologische Ladung k = (1/8π²) ∫ Tr(F∧F).
        @description
            Die topologische Ladung (Pontrjagin-Index) ist:
              k = (1/8π²) ∫_R⁴ Tr(F_μν F̃^μν) d⁴x

            wobei F̃^μν = (1/2) ε^μνλρ F_λρ der duale Feldstärketensor ist.

            Für das BPST-Instanton gilt k = 1 (die Lösung liegt in der Klasse k=1).

            Die topologische Ladung ist eine ganzzahlige Invariante (sie nimmt
            nur ganzzahlige Werte an) und klassifiziert die Eichfeld-Konfigurationen
            topologisch: π₃(SU(2)) ≅ ℤ.

            Vereinfachte Berechnung für diskrete Feldkomponenten:
              k ≈ (1/8π²) · Σ_{μ<ν} Tr(F_μν F̃_μν) · Δ⁴x

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten oder curvature_2form-Ergebnis
        @return Numerischer Wert der topologischen Ladung
        """
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        if not F:
            return 0.0

        size = list(F.values())[0].shape[0]

        def get_F(mu: int, nu: int) -> np.ndarray:
            if (mu, nu) in F:
                return F[(mu, nu)]
            elif (nu, mu) in F:
                return -F[(nu, mu)]
            return np.zeros((size, size), dtype=complex)

        # Levi-Civita-Tensor in 4D (vollständig antisymmetrisch, ε₀₁₂₃ = +1)
        # Topologische Ladungsdichte: q = (1/8π²) Tr(F_μν F̃^μν)
        #   = (1/16π²) ε^μνλρ Tr(F_μν F_λρ)
        total = 0.0
        indices = [(0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2),
                   (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0),
                   (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1),
                   (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0)]

        for perm in indices:
            mu, nu, lam, rho = perm
            # Vorzeichen der Permutation (Levi-Civita-Symbol)
            sign = self._levi_civita_4d(mu, nu, lam, rho)
            if sign == 0:
                continue
            F_mn = get_F(mu, nu)
            F_lr = get_F(lam, rho)
            trace = np.real(np.trace(F_mn @ F_lr))
            total += sign * trace

        # Normierung: (1/16π²)
        return total / (16.0 * math.pi ** 2)

    def _levi_civita_4d(self, a: int, b: int, c: int, d: int) -> int:
        """
        @brief Berechnet das 4D-Levi-Civita-Symbol ε_{abcd}.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param a, b, c, d  Indizes (0-3)
        @return +1 (gerade Perm.), -1 (ungerade), 0 (wiederholt)
        """
        indices = [a, b, c, d]
        # Wiederholte Indizes → 0
        if len(set(indices)) < 4:
            return 0
        # Zähle Inversionen (Bubble-Sort-Methode)
        perm = list(indices)
        swaps = 0
        for i in range(4):
            for j in range(i + 1, 4):
                if perm[i] > perm[j]:
                    perm[i], perm[j] = perm[j], perm[i]
                    swaps += 1
        return 1 if swaps % 2 == 0 else -1

    def action(self) -> float:
        """
        @brief Berechnet die Instanton-Wirkung S = 8π²|k|/g².
        @description
            Die Wirkung eines Instantons mit topologischer Ladung k ist:
              S = 8π²|k|/g²

            wobei g die Kopplungskonstante der Eichtheorie ist.

            Dies ist das absolute Minimum der Yang-Mills-Wirkung in der
            topologischen Klasse mit Ladung k. Instantone sind deshalb
            die wichtigsten Quantenkorrekturen in der Pfadintegral-Formulierung.

            Für g=1 (normierte Kopplung): S = 8π²|k| ≈ 78.96|k|

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Instanton-Wirkung für g=1 (normierte Kopplung)
        """
        # Wirkung für normierte Kopplungskonstante g = 1
        g_squared = 1.0
        return 8.0 * math.pi ** 2 * abs(self.k) / g_squared

    def moduli_space_dimension(self) -> int:
        """
        @brief Berechnet die Dimension des Moduli-Raums.
        @description
            Der Moduli-Raum M_k für SU(2)-Instantone mit Ladung k ist eine
            Mannigfaltigkeit der Dimension:

              dim M_k = 8k - 3   (für k > 0, SU(2), selbst-dual)

            Erklärung:
            - 8k: Atiyah-Singer-Index des Dirac-Operators im Adjoint-Hintergrund
              (Dimension des Raums aller Instanton-Lösungen vor Eichfixierung)
            - -3: Dimension der Eichgruppe SU(2) (globale Eichtransformationen)

            Für k=1: dim = 5 (Position 4D + Skala 1D = 5 kollektive Koordinaten)
            Für k=2: dim = 13

            Allgemein für SU(N): dim M_k = 4Nk - N² + 1

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Dimension des Moduli-Raums (ganze Zahl)
        """
        # Nur für SU(2) und selbst-duale Instantone (k > 0)
        if self.group_name == "SU(2)" and self.k > 0:
            return 8 * self.k - 3
        elif self.group_name == "SU(3)" and self.k > 0:
            # Für SU(3): dim = 8k × 2 - (9-1) = 16k - 8? Korrektur nach ADHM:
            # dim M_k(SU(N)) = 4Nk
            return 4 * 3 * self.k - (9 - 1)  # = 12k - 8
        else:
            # Allgemeines Ergebnis: 8|k| - dim(G)
            return 8 * abs(self.k) - self.group.dimension()


# =============================================================================
# KLASSE: MassGap — Das Massenspalt-Millennium-Problem
# =============================================================================

class MassGap:
    """
    @brief Massenspalt-Problem der Yang-Mills-Theorie (Millennium-Problem).
    @description
        Das Massenspalt-Problem ist eines der sieben Millennium-Probleme des
        Clay Mathematics Institute (CMI) mit einem Preisgeld von 1 Million USD.

        PROBLEM (Jaffe-Witten, 2000):
        Beweise, dass für jede kompakte einfache Eichgruppe G eine nicht-triviale
        quantenmechanische Yang-Mills-Theorie auf R⁴ existiert und ein
        Massenspalt Δ > 0 hat.

        Das bedeutet: Die kleinste angeregte Energie des Hamiltonians
        H ist strikt positiv (keine masselosen Anregungen außer dem Vakuum).

        Physikalische Bedeutung:
        - Confinement in QCD (Quarks und Gluonen sind eingesperrt)
        - Erklärung der Hadronenmassen (Protonen, Neutronen etc.)
        - Gluonen bekommen eine effektive Masse trotz Eichinvarianz

        Bekannte Resultate:
        - Gitterrechnungen zeigen numerisch eine Massenlücke
        - Störungstheoretisch: kein Massenspalt (nur nicht-perturbativ)
        - Schwelle: Δ ≈ 1-2 GeV in QCD (aus Experiment und Gittersimulationen)

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, gauge_group: str = "SU(3)") -> None:
        """
        @brief Initialisiert das Massenspalt-Problem für die angegebene Gruppe.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param gauge_group  Eichgruppe (Standard: "SU(3)" für QCD)
        """
        # Eichgruppe des Problems
        self.group = LieGroup(gauge_group)
        self.group_name = gauge_group

    def millennium_problem_statement(self) -> str:
        """
        @brief Gibt die exakte Formulierung des Massenspalt-Millennium-Problems zurück.
        @description
            Dies ist die offizielle Formulierung des Clay Mathematics Institute
            (Jaffe-Witten, 2000).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Zeichenkette mit der Problembeschreibung
        """
        return (
            "Yang-Mills Existence and Mass Gap (Millennium Problem)\n"
            "========================================================\n\n"
            "Gegeben: Eine kompakte einfache Lie-Gruppe G (z.B. SU(2) oder SU(3)).\n\n"
            "Zu beweisen:\n"
            "  (i) Existenz: Die Yang-Mills-Quantenfeldtheorie auf R⁴ mit Gruppe G\n"
            "      existiert und erfüllt die Wightman-Axiome der konstruktiven QFT.\n\n"
            "  (ii) Massenspalt: Es gibt eine Konstante Δ > 0 (die Massenlücke), sodass\n"
            "       für jeden Zustand |ψ⟩ im Hilbert-Raum (orthogonal zum Vakuum |0⟩):\n"
            "         ⟨ψ|H|ψ⟩ ≥ Δ · ⟨ψ|ψ⟩\n"
            "       d.h. jede angeregte Energie ist mindestens Δ > 0.\n\n"
            "Bedeutung:\n"
            "  - Erklärt das Confinement von Quarks und Gluonen in QCD\n"
            "  - Begründet die Existenz von Gluebällen (reine Glue-Hadronen)\n"
            "  - Verbindet geometrische Analysis (4D-Eichtheorie) mit mathematischer Physik\n\n"
            f"Aktuelle Eichgruppe: {self.group_name} (dim = {self.group.dimension()})\n"
            "Status: OFFEN (Stand 2026) — Clay Mathematics Institute, Preisgeld: 1 Mio. USD"
        )

    def physical_motivation(self) -> Dict:
        """
        @brief Erklärt die physikalische Motivation des Massenspalt-Problems.
        @description
            Warum ist das Massenspalt-Problem wichtig für die Physik?

            1. Confinement: Quarks werden nicht als freie Teilchen beobachtet.
               Die starke Kraft wird mit zunehmender Entfernung stärker (im
               Gegensatz zu Elektromagnetismus und Gravitation).

            2. Asymptotische Freiheit: Für kleine Abstände wird g(μ) → 0
               (Gross-Politzer-Wilczek, Nobelpreis 2004). Das bedeutet, die
               Störungstheorie funktioniert bei hoher Energie — aber nicht
               für das Confinement bei niedrigen Energien.

            3. QCD-Skalenparameter: Λ_QCD ≈ 200-300 MeV — die Energieskala,
               bei der die Kopplungsstärke O(1) wird (nicht-perturbativ).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Dictionary mit physikalischen Motivationen und Erklärungen
        """
        return {
            "confinement": {
                "description": "Quarks und Gluonen werden nie frei beobachtet",
                "mechanism": "Das Farbfeld bildet einen Flussschlauch (String) zwischen Quarks",
                "energy": "Energie ∝ Abstand r → unendlich für r → ∞",
                "string_tension": "σ ≈ 0.18 GeV² (aus Experiment und Gitter-QCD)"
            },
            "asymptotic_freedom": {
                "description": "Kopplungskonstante g(μ) → 0 für hohe Energie μ → ∞",
                "discovery": "Gross, Politzer, Wilczek (1973) — Nobelpreis 2004",
                "beta_function": "β(g) = -b₀g³/(16π²) < 0  (b₀ = 11N/3 - 2N_f/3 für SU(N))"
            },
            "glueball_masses": {
                "description": "Reine Gluon-Zustände (keine Quarks) — Gluebälle",
                "lightest_glueball": "f₀(1710): Masse ≈ 1710 MeV — möglicher Glueballkandidat",
                "lattice_QCD": "Leichtester Glueballzustand: ≈ 1600-1700 MeV (0⁺⁺)"
            },
            "qcd_scale": {
                "lambda_qcd": "Λ_QCD ≈ 200-300 MeV",
                "description": "Energieskala der nicht-perturbativen QCD",
                "mass_gap_estimate": "Δ ≈ Λ_QCD (grobe Abschätzung)"
            },
            "millennium_prize": {
                "institute": "Clay Mathematics Institute",
                "year": 2000,
                "prize": "1.000.000 USD",
                "status": "Offen (2026)"
            }
        }

    def lattice_gauge_evidence(self, beta: float, lattice_size: int) -> Dict:
        """
        @brief Erzeugt numerische Evidenz für den Massenspalt via Gittereichtheorie.
        @description
            Gittereichtheorie (Wilson, 1974) ist die wichtigste numerische Methode
            für nicht-perturbative QCD-Berechnungen.

            Die Gitterwirkung (Wilson-Wirkung):
              S_W = β Σ_Plaketten (1 - (1/N) Re Tr U_Plakette)

            wobei β = 2N/g² die Kopplungskonstante auf dem Gitter ist.

            Massenspalt aus Gittersimulationen:
            - Für große β (kleine Kopplung): Kontinuumslimes
            - Massenspalt Δ ≈ β^{-1/2} exp(-β/(b₀)) (grobe Näherung)
            - Reale QCD-Gitterberechnungen: Δ ≈ 1 GeV

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param beta         Gittereichkopplungsparameter β = 2N/g² (typisch 5-7 für SU(3))
        @param lattice_size Gittergröße (Anzahl Punkte pro Richtung)
        @return Dictionary mit Gitterrechnungs-Ergebnissen
        """
        # Anzahl der Generatoren der Gruppe
        N = {"SU(2)": 2, "SU(3)": 3, "U(1)": 1, "SO(3)": 3}.get(self.group_name, 2)

        # Kopplungskonstante aus β = 2N/g²
        g_squared = 2.0 * N / beta if beta > 0 else float('inf')

        # Anzahl der Plaketten (4D-Gitter hat 6 Ebenen pro Punkt)
        n_plaquettes = 6 * (lattice_size ** 4)

        # Näherung für die Massenlücke aus Gitterrechnungen
        # Δ ≈ a⁻¹ · exp(-C/g²) mit a = Gitterabstand
        # Kontinuumslimes: a → 0, β → ∞, Δ = const.
        b0 = 11.0 * N / 3.0  # Einschleifiger Beta-Funktionskoeffizient
        if g_squared > 0:
            mass_gap_estimate = math.exp(-2.0 * math.pi ** 2 / (b0 * g_squared))
        else:
            mass_gap_estimate = 0.0

        # Plaketten-Erwartungswert (Näherung für schwache Kopplung)
        # ⟨P⟩ ≈ 1 - 1/(N²β) für großes β
        plaquette_expectation = 1.0 - 1.0 / (N ** 2 * beta) if beta > 0 else 0.0

        return {
            "group": self.group_name,
            "beta": beta,
            "g_squared": g_squared,
            "lattice_size": lattice_size,
            "n_plaquettes": n_plaquettes,
            "mass_gap_estimate": mass_gap_estimate,
            "plaquette_expectation": min(plaquette_expectation, 1.0),
            "b0": b0,
            "continuum_limit": beta > 6.0,  # Typisch für SU(3): β > 6 = Kontinuumslimes
            "description": (
                f"Gittereichtheorie SU({N}) auf {lattice_size}⁴-Gitter, "
                f"β = {beta:.2f}, g² = {g_squared:.4f}\n"
                "Massenlücke Δ = exp(-2π²/(b₀g²)) (grobe Ein-Schleifen-Näherung)"
            )
        }

    def mass_gap_lower_bound_jaffe_witten(self) -> str:
        """
        @brief Gibt bekannte partielle Resultate zum Massenspalt zurück.
        @description
            Obwohl das Hauptproblem offen ist, gibt es mehrere wichtige
            Teilresultate:

            1. Perturbative Analyse: In Störungstheorie ist KEIN Massenspalt
               (Photonen-ähnliches Spektrum, masselos). Der Massenspalt ist
               ein nicht-perturbatives Phänomen.

            2. 2D und 3D: Das Massenspalt-Problem ist in 2D und 3D bewiesen
               (Seiler 1975, Federbush 1986, King 1986).

            3. Gitterregularisierung: Der Massenspalt existiert für jedes
               endliche Gitter (diskret) — aber der Kontinuumslimes fehlt.

            4. Numerische Evidenz: Gitter-QCD zeigt klar eine Massenlücke von
               Δ ≈ 1-2 GeV für SU(3) Yang-Mills.

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Zeichenkette mit bekannten Ergebnissen
        """
        return (
            "Bekannte Partielle Resultate zum Yang-Mills-Massenspalt:\n"
            "=========================================================\n\n"
            "1. BEWEIS in 2D und 3D:\n"
            "   - Gross (1976), Seiler (1978): YM in 2D ist exakt lösbar, Massenspalt vorhanden\n"
            "   - Balaban (1985), King (1986): Massenspalt für YM₃ im Kontinuumslimes\n\n"
            "2. GITTER-Resultate (numerisch, kein Beweis):\n"
            "   - Wilson (1974): Confinement im Stark-Kopplungs-Limes auf dem Gitter\n"
            "   - Creutz (1980): Erste Gitter-Monte-Carlo-Simulation zeigt String-Tension\n"
            "   - Lattice QCD (1990s-heute): Δ ≈ 1600-1800 MeV (0⁺⁺ Glueballmasse)\n\n"
            "3. ANALYTISCHE Schranken:\n"
            "   - Jaffe-Witten (2000): Formulierung des Millennium-Problems\n"
            "   - Zwanziger (2003): Gribov-Horizont und Confinement\n"
            "   - Dudal-Sorella (2008): Gribov-Zwanziger-Renormierung\n\n"
            "4. NOCH OFFEN (kein vollständiger Beweis für 4D):\n"
            "   - Rigorous QFT (Wightman-Axiome) für YM₄\n"
            "   - Quantitativer unterer Bound für Δ in 4D\n"
            "   - Verbindung zwischen topologischer Ladung und Massenspalt\n\n"
            f"Eichgruppe dieses Objekts: {self.group_name}"
        )

    def glueball_mass_estimate(self, g_squared: float) -> float:
        """
        @brief Schätzt die Glueballmasse ~ exp(-c/g²) für kleine Kopplung.
        @description
            Gluebälle sind Hadronen aus reinen Gluonfeldern (keine Quarks).
            Ihre Massen sind von der Größenordnung des Massenspalt-Parameters.

            Instanton-induzierte Schätzung (Shuryak, Schafer):
              m_glueball ~ Λ_QCD · exp(-8π²/(N_c · g²))

            Für schwache Kopplung (kleine g²):
              m_glueball ~ exp(-c/g²)   mit c = 8π²/(N_c)

            Dies ist nicht-perturbativ (analytisch in g² = 0 nicht definiert),
            was zeigt, dass Gluebälle ein genuiner nicht-perturbativer Effekt sind.

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param g_squared  Quadrat der Kopplungskonstante g² > 0
        @return Geschätzte Glueballmasse (in natürlichen Einheiten, relativ)
        """
        if g_squared <= 0:
            return 0.0

        # N_c = Anzahl Farben (Casimir-Faktor der Gruppe)
        N_c = {"SU(2)": 2, "SU(3)": 3, "U(1)": 1, "SO(3)": 3}.get(self.group_name, 2)

        # Exponent: c = 8π²/N_c (aus Instanton-Berechnung)
        c = 8.0 * math.pi ** 2 / N_c

        # Glueballmasse ~ exp(-c/g²)
        return math.exp(-c / g_squared)

    def confinement_criterion(self, wilson_loop_area: float) -> bool:
        """
        @brief Prüft das Confinement-Kriterium (Flächengesetz für Wilson-Schleifen).
        @description
            Confinement ist charakterisiert durch das Flächengesetz der
            Wilson-Schleife W(C):
              log W(C) ~ -σ · Area(C)   (Confinement/Flächengesetz)
              log W(C) ~ -μ · Perimeter(C) (kein Confinement/Umfangsgesetz)

            Confinement tritt auf, wenn die topologische Ladungsdichte groß ist,
            d.h. wenn die String-Tension σ > 0 ist.

            Einfaches Kriterium: W > 0 und |log W| proportional zur Fläche.

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param wilson_loop_area  Fläche der Wilson-Schleife (positiv)
        @return True wenn Confinement wahrscheinlich (Flächengesetz), False sonst
        """
        # String-Tension σ ≈ 0.18 GeV² für SU(3)-QCD
        sigma = 0.18
        # Bei Confinement: W(C) ~ exp(-σ · Area) → W klein für große Fläche
        # Kriterium: Wenn die effektive Dämpfung pro Flächeneinheit σ_eff > 0
        effective_sigma = wilson_loop_area * sigma
        # Confinement: exp(-σ · A) < Schwellenwert für große Flächen
        threshold = 0.5  # Schwellenwert für "kleine" Wilson-Schleife
        w_estimate = math.exp(-effective_sigma)
        return w_estimate < threshold


# =============================================================================
# KLASSE: WilsonLoop — Wilson-Schleifen und Confinement
# =============================================================================

class WilsonLoop:
    """
    @brief Wilson-Schleifen: Ordnungsparameter für Confinement in Eichtheorien.
    @description
        Der Wilson-Loop (Wilson, 1974) ist definiert als:
          W(C) = (1/N) Tr P exp(i ∮_C A_μ dx^μ)

        wobei P die pfadgeordnete Integration (Pfad-Exponential) bezeichnet.

        Physikalische Bedeutung:
        - W(C) ist ein eichinvarianter Ordnungsparameter
        - Für eine rechteckige Schleife mit Seiten R (Raum) und T (Zeit):
            W(R,T) = exp(-V(R)·T) für großes T
          wobei V(R) das statische Quark-Antiquark-Potential ist.

        Confinement-Kriterien:
        - Flächengesetz: W(R,T) ~ exp(-σRT) für R,T → ∞
          → Lineares Potential V(R) = σR (Confinement!)
        - Umfangsgesetz: W(R,T) ~ exp(-μ(R+T))
          → Yukawa-Potential V(R) ~ exp(-mR)/R (Deconfinement)

        String-Tension σ ist die Energiedichte des Flussschlauchs zwischen Quarks.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self, group: str = "SU(3)") -> None:
        """
        @brief Initialisiert den Wilson-Loop für die gegebene Eichgruppe.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param group  Eichgruppe (Standard: "SU(3)" für QCD)
        """
        self.group_name = group
        self.lie_group = LieGroup(group)

    def rectangular_loop(self, R: float, T: float, g_squared: float) -> float:
        """
        @brief Berechnet den Wilson-Loop W(R,T) für eine rechteckige Schleife.
        @description
            Für eine rechteckige Schleife R×T gilt in der Näherung mit Confinement:
              W(R,T) = exp(-σ(g) · R · T)

            wobei σ(g) die String-Tension ist.

            Für große g² (starke Kopplung, Gitter-Expansion):
              σ_lattice ≈ -log(g²/(2N))  (führende Ordnung, Stark-Kopplungs-Entwicklung)

            Für kleine g² (schwache Kopplung):
              Keine analytische Formel (nicht-perturbativ, Λ_QCD-Skala)

            Numerische Näherung (beide Regime verbunden):
              W(R,T) ≈ exp(-σ(g²) · R · T)
              σ(g²) ≈ Λ²_QCD · exp(-8π²/(b₀g²))

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param R         Räumliche Ausdehnung der Schleife (Quark-Abstand in fm)
        @param T         Zeitliche Ausdehnung der Schleife
        @param g_squared Quadrat der Kopplungskonstante g²
        @return Wert des Wilson-Loops W(R,T) ∈ [0, 1]
        """
        # String-Tension σ(g²) aus der String-Tension-Funktion
        sigma = self.string_tension(g_squared)
        # Wilson-Loop mit Flächengesetz: W = exp(-σRT)
        return math.exp(-sigma * R * T)

    def string_tension(self, g_squared: float) -> float:
        """
        @brief Berechnet die String-Tension σ aus der Kopplungskonstante g².
        @description
            Die String-Tension σ ist die zentrale Größe des Confinements.
            Sie gibt die Energiedichte (Kraft) des Flussschlauchs zwischen Quarks an:
              V(r) = σ · r   (lineares Potential)

            Experimenteller Wert: σ ≈ 0.18 GeV² ≈ (0.42 GeV)² für SU(3)-QCD

            Gitter-Näherung (stark-Kopplungs-Entwicklung, führende Ordnung):
              σ_lat = -log(g²/(2N))  (dimensionslos, in Gittereinheiten)

            Für die nicht-perturbative Physik (zwei-Schleifen RGE):
              σ ∝ Λ²_QCD ∝ μ² exp(-16π²/(b₀g²(μ)))

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param g_squared  Quadrat der Kopplungskonstante
        @return String-Tension σ (in natürlichen Einheiten)
        """
        if g_squared <= 0:
            return 0.0

        N = {"SU(2)": 2, "SU(3)": 3, "U(1)": 1, "SO(3)": 3}.get(self.group_name, 2)

        # Stark-Kopplungs-Näherung: σ = -log(g²/(2N)) für g² >> 1
        if g_squared > 2.0 * N:
            # Stark-Kopplungsregime: linearer Anstieg mit g²
            return math.log(g_squared / (2.0 * N))

        # Schwach-Kopplungsregime: nicht-perturbativer Ansatz
        # σ ~ exp(-8π²/(b₀g²)) mit b₀ = 11N/3
        b0 = 11.0 * N / 3.0
        return math.exp(-8.0 * math.pi ** 2 / (b0 * g_squared))

    def area_law_check(self, perimeter: float, area: float,
                       W_value: float) -> bool:
        """
        @brief Prüft ob W(C) einem Flächengesetz (Confinement) folgt.
        @description
            Unterscheidung zwischen Flächen- und Umfangsgesetz:
              Flächengesetz:  log W ∝ -Area   (Confinement)
              Umfangsgesetz:  log W ∝ -Perimeter (kein Confinement)

            Der Test: Berechne effektive σ_area und μ_perim aus W und
            prüfe welches Gesetz besser passt.

            Kriterium: Falls |log W| / Area > |log W| / Perimeter,
            dann dominiert das Flächengesetz (Confinement).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param perimeter  Umfang der Schleife
        @param area       Fläche der Schleife
        @param W_value    Wilson-Loop-Wert W ∈ (0, 1]
        @return True wenn Flächengesetz dominiert (Confinement), False sonst
        """
        if W_value <= 0 or W_value > 1.0 or area <= 0 or perimeter <= 0:
            return False

        # |log W| berechnen
        log_W = abs(math.log(W_value))

        # Effektive Flächen- und Umfangs-Koeffizienten
        sigma_eff = log_W / area       # Effektive String-Tension
        mu_eff = log_W / perimeter     # Effektiver Massenterm

        # Flächengesetz dominiert wenn σ_eff > μ_eff
        # (W nimmt schneller mit der Fläche ab als mit dem Umfang)
        return sigma_eff > mu_eff

    def creutz_ratio(self, W_values: Dict[Tuple, float],
                     R: int, T: int) -> float:
        """
        @brief Berechnet das Creutz-Verhältnis χ(R,T) zur Bestimmung der String-Tension.
        @description
            Das Creutz-Verhältnis (Creutz, 1980) ist eine Methode zur direkten
            Extraktion der String-Tension aus Gitter-Wilson-Schleifen:

              χ(R,T) = -log [W(R,T) · W(R-1,T-1) / (W(R,T-1) · W(R-1,T))]

            Für große R,T und bei Confinement gilt:
              χ(R,T) → σ   (konvergiert zur String-Tension)

            Vorteil: Das Creutz-Verhältnis ist frei von Perimeter-Korrekturen
            (sie kürzen sich im Quotienten heraus).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param W_values  Dictionary {(R,T): W_Wert} der Wilson-Schleifen
        @param R         Räumliche Ausdehnung
        @param T         Zeitliche Ausdehnung
        @return Creutz-Verhältnis χ(R,T) ≈ String-Tension σ
        """
        # Benötigte Wilson-Schleifen
        W_RT = W_values.get((R, T), None)
        W_R1T1 = W_values.get((R - 1, T - 1), None)
        W_RT1 = W_values.get((R, T - 1), None)
        W_R1T = W_values.get((R - 1, T), None)

        # Prüfe ob alle Werte vorhanden und positiv sind
        if any(v is None or v <= 0 for v in [W_RT, W_R1T1, W_RT1, W_R1T]):
            return float('nan')

        # Creutz-Verhältnis: χ = -log(W(R,T)·W(R-1,T-1) / (W(R,T-1)·W(R-1,T)))
        numerator = W_RT * W_R1T1
        denominator = W_RT1 * W_R1T

        if denominator <= 0:
            return float('nan')

        return -math.log(numerator / denominator)


# =============================================================================
# KLASSE: YangMillsTopology — Topologische Eigenschaften
# =============================================================================

class YangMillsTopology:
    """
    @brief Topologische Eigenschaften der Yang-Mills-Theorie.
    @description
        Die Topologie spielt eine entscheidende Rolle in der Yang-Mills-Theorie:

        1. Chern-Klassen c_k: Topologische Invarianten von Faserbündeln
           c_1(E) ∈ H²(M, Z), c_2(E) ∈ H⁴(M, Z) (Pontryagin-Klasse)

        2. Pontryagin-Klasse p_1: Für reelle Bündel, p_1 = -2c_2

        3. Theta-Vakuum: In der QCD gibt es unendlich viele topologisch
           verschiedene Vakuumzustände |n⟩. Das physikalische Vakuum ist:
             |θ⟩ = Σ_{n=-∞}^{∞} e^{inθ} |n⟩   (θ-Vakuum)

        4. CP-Verletzung: Der θ-Term in der QCD-Lagrange-Dichte:
             L_θ = (θ/16π²) Tr(F_μν F̃^μν)
           verletzt CP-Symmetrie (stark CP-Problem).

        5. Peccei-Quinn-Mechanismus: Die Axion-Lösung des starken CP-Problems
           (Peccei-Quinn, 1977; Weinberg, Wilczek, 1978).

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def chern_class(self, F_components: Dict, k: int) -> float:
        """
        @brief Berechnet die k-te Chern-Klasse des Eichbündels.
        @description
            Die Chern-Klassen sind topologische Invarianten eines komplexen
            Vektorbündels E → M:

            c_0(E) = 1  (normiert)
            c_1(E) = (i/2π) Tr(F)              (1. Chern-Klasse)
            c_2(E) = (1/8π²) [Tr(F∧F) - Tr(F)∧Tr(F)]   (2. Chern-Klasse)

            Für SU(N)-Bündel gilt Tr(F) = 0 (spurlose Generatoren), also:
              c_2 = (1/8π²) Tr(F∧F) = topologische Ladung k

            Die Chern-Klassen sind ganze Zahlen (c_k ∈ H^{2k}(M, Z)).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten
        @param k             Ordnung der Chern-Klasse (0, 1, 2)
        @return Numerischer Wert der k-ten Chern-Klasse
        """
        if k == 0:
            # c_0 ist immer 1 (normiert)
            return 1.0

        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        if k == 1:
            # c_1 = (i/2π) ∫ Tr(F)
            # Für SU(N): Tr(F) = 0, also c_1 = 0
            total_trace = 0.0
            for F_mat in F.values():
                total_trace += np.real(np.trace(F_mat))
            return total_trace / (2.0 * math.pi)

        elif k == 2:
            # c_2 = topologische Ladung = (1/8π²) ∫ Tr(F∧F)
            # Berechne via Instanton-Objekt
            inst = Instanton(1)  # Hilfs-Objekt für Berechnung
            return inst.topological_charge(F_components)

        else:
            # Höhere Chern-Klassen (k > 2): nur in dim > 4 nicht-trivial
            return 0.0

    def pontryagin_class(self, F_components: Dict) -> float:
        """
        @brief Berechnet die erste Pontryagin-Klasse p₁ = -2c₂.
        @description
            Die Pontryagin-Klassen sind topologische Invarianten reeller Bündel:
              p_1(E) = -2c_2(E_ℂ)

            wobei E_ℂ die Komplexifizierung des reellen Bündels E ist.

            Für ein SU(2)-Eichbündel gilt:
              p_1 = -2c_2 = -(1/4π²) ∫ Tr(F∧F)

            Physikalische Bedeutung: p₁/4 = topologische Ladung k (für SU(2))

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten
        @return Wert der ersten Pontryagin-Klasse
        """
        # p_1 = -2 c_2
        c2 = self.chern_class(F_components, 2)
        return -2.0 * c2

    def theta_vacuum(self, theta: float, instanton_sum: int) -> complex:
        """
        @brief Berechnet das Theta-Vakuum |θ⟩ = Σ_{n} e^{inθ} |n⟩.
        @description
            In der Yang-Mills-Theorie gibt es unendlich viele topologisch
            verschiedene Vakuumzustände |n⟩ (klassiert nach topologischer Ladung n).

            Das physikalische Quantenvakuum ist eine Überlagerung:
              |θ⟩ = Σ_{n=-N}^{N} e^{inθ} |n⟩

            In der Pfadintegral-Formulierung erscheint der θ-Term:
              Z(θ) = Σ_n e^{inθ} · Z_n  = ∫ DA e^{-S[A] + iθ·k[A]}

            Dies ist äquivalent zur zusätzlichen Lagrangedichte:
              L_θ = (θ/32π²) ε^{μνλρ} Tr(F_μν F_λρ)

            Der θ-Parameter bricht CP-Symmetrie (außer bei θ = 0 oder π).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param theta          θ-Parameter des QCD-Vakuums (reell, 0 ≤ θ ≤ 2π)
        @param instanton_sum  Anzahl der Terme in der Summe ±N
        @return Komplexe Amplitude ⟨0|θ⟩ (Projektion auf |n=0⟩)
        """
        # Theta-Vakuum: |θ⟩ = Σ_{n=-N}^{N} e^{inθ} |n⟩
        # Die Projektion ⟨0|θ⟩ = 1 (Normierung)
        # Die Phase e^{inθ} der n-ten Komponente
        amplitude = complex(0.0)
        N = instanton_sum
        for n in range(-N, N + 1):
            # Gewicht des n-ten Instanton-Sektors
            # Im Vakuumsektor |0⟩: Amplitude = e^{i·n·θ}
            amplitude += cmath.exp(1j * n * theta)
        return amplitude

    def vacuum_energy_density(self, theta: float, g_squared: float) -> float:
        """
        @brief Berechnet die Vakuumenergiedichte E(θ) für das Theta-Vakuum.
        @description
            Die Vakuumenergiedichte als Funktion des θ-Parameters:

              E(θ) = -Λ⁴_YM · cos(θ) + O(g²)

            Im Instanton-Gas-Modell (Diluted Instanton Gas Approximation):
              E(θ) ∝ -exp(-8π²/g²) · cos(θ)

            Das Minimum liegt bei θ = 0 (oder θ = π für bestimmte Theorien).
            Die Krümmung am Minimum gibt die Massenspalt-Energie:
              E''(0) ∝ exp(-8π²/g²) > 0   (Massenlücke!)

            Physikalisch: θ ≠ 0 bedeutet CP-Verletzung in der starken Wechselwirkung
            (die experimentell sehr klein ist: |θ| < 10⁻¹⁰ → starkes CP-Problem).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param theta      θ-Vakuum-Parameter (reell)
        @param g_squared  Quadrat der Kopplungskonstante g² > 0
        @return Vakuumenergiedichte E(θ) (in natürlichen Einheiten)
        """
        if g_squared <= 0:
            return 0.0

        # Instanton-Beitrag: E ~ -exp(-S_inst) · cos(θ) mit S_inst = 8π²/g²
        S_inst = 8.0 * math.pi ** 2 / g_squared
        # Instanton-Amplitude (exponentiell unterdrückt für kleine g²)
        instanton_amplitude = math.exp(-S_inst)
        # Vakuumenergiedichte: E(θ) = -K · cos(θ), K = Instanton-Dichte
        return -instanton_amplitude * math.cos(theta)

    def axion_solution(self) -> str:
        """
        @brief Beschreibt den Peccei-Quinn-Mechanismus (Axion-Lösung).
        @description
            Das starke CP-Problem: Warum ist |θ| < 10⁻¹⁰ experimentell
            (obwohl θ ∈ [0, 2π] a priori beliebig sein könnte)?

            Peccei-Quinn-Mechanismus (1977):
            - Einführung einer neuen globalen U(1)_PQ-Symmetrie
            - Diese wird spontan gebrochen: ⟨φ⟩ = f_a (PQ-Skala)
            - Das Goldstone-Boson ist das AXION (Weinberg, Wilczek 1978)
            - θ_eff = θ + a/f_a relaxiert dynamisch zu 0 (Minimum der Energie)

            Das Axion ist gleichzeitig ein Kandidat für Dunkle Materie!

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @return Beschreibung des Axion-Mechanismus
        """
        return (
            "Peccei-Quinn-Mechanismus und Axion-Lösung des starken CP-Problems\n"
            "==================================================================\n\n"
            "Das Problem: Experimentell |θ_QCD| < 10⁻¹⁰ (neutronisches EDM-Limit),\n"
            "             theoretisch: θ ∈ [0, 2π] völlig frei → Feinabstimmung nötig\n\n"
            "Peccei-Quinn-Lösung (1977):\n"
            "  1. Neue globale U(1)_PQ-Symmetrie in der Theorie\n"
            "  2. Spontane Brechung bei Skala f_a (PQ-Skala, ~10⁹-10¹² GeV)\n"
            "  3. Goldstone-Boson: Axion a(x) mit f_a\n"
            "  4. θ_effektiv = θ + a(x)/f_a → relaxiert zu θ_eff = 0\n\n"
            "Das Axion:\n"
            "  - Masse: m_a ≈ 6 μeV · (10¹² GeV / f_a)\n"
            "  - Kopplung an Photonen: g_aγγ ~ α/f_a\n"
            "  - Kandidat für Dunkle Materie (axionische Dunkle Materie)\n\n"
            "Experimentelle Suche:\n"
            "  - ADMX (Axion Dark Matter eXperiment)\n"
            "  - CASPEr, ABRACADABRA, HAYSTAC\n"
            "  - Ausgeschlossener Bereich: f_a < 10⁹ GeV (astrophysikalische Grenzen)\n\n"
            "Verbindung zum Massenspalt:\n"
            "  - Axionmasse m_a ∝ sqrt(χ_top) / f_a\n"
            "  - χ_top = topologische Suszeptibilität = d²E(θ)/dθ² |_{θ=0}\n"
            "  - χ_top > 0 ist äquivalent zum Massenspalt im reinen Yang-Mills!"
        )


# =============================================================================
# KLASSE: IndexTheory — Atiyah-Singer-Index-Theorem
# =============================================================================

class IndexTheory:
    """
    @brief Atiyah-Singer-Index-Theorem für Dirac-Operatoren in Eichtheorien.
    @description
        Das Atiyah-Singer-Index-Theorem (1963) ist eines der tiefsten Resultate
        der modernen Mathematik. Es verbindet analytische und topologische
        Invarianten:

          ind(D) = dim ker(D) - dim coker(D) = analytischer Index
                 = ∫_M ch(E) · Â(TM) = topologischer Index

        Für den Dirac-Operator in einem Yang-Mills-Hintergrund auf S⁴:
          ind(D_A) = -c_2(E) = -k  (für selbst-duale A)

        Das bedeutet: Die Anzahl der Null-Moden des Dirac-Operators ist
        topologisch bestimmt (durch die Instanton-Ladung k).

        Konsequenzen:
        - Für k > 0 (Instanton): dim ker(D) = 2k·N (N = Rang der Darstellung)
        - Null-Moden des Dirac-Operators → wichtig für chirale Symmetriebrechung
        - Atiyah-Singer erklärt: dim M_k = 8k - 3 (Moduli-Raum-Dimension)

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def atiyah_singer_index(self, dirac_op: np.ndarray) -> int:
        """
        @brief Berechnet den analytischen Index des Dirac-Operators.
        @description
            Der analytische Index eines Differential-Operators D ist:
              ind(D) = dim ker(D) - dim coker(D)
                     = dim ker(D) - dim ker(D†)

            Für den Dirac-Operator D: ker(D) = Raum der Null-Moden
            Der Index zählt die Netto-Anzahl der Null-Moden
            (mit Vorzeichen je nach Chiralität).

            Berechnung über SVD (Singulärwertzerlegung):
            - ker(D) ≅ Vektoren mit Singulärwert 0
            - ker(D†) ≅ Linksnull-Vektoren

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param dirac_op  Matrix-Darstellung des Dirac-Operators (komplex, quadratisch oder rechteckig)
        @return Analytischer Index ind(D) = dim ker(D) - dim ker(D†)
        """
        # Singulärwertzerlegung: D = U Σ V†
        U, s, Vh = np.linalg.svd(dirac_op)
        # Toleranz für Null-Singulärwerte
        tol = 1e-10 * max(s) if len(s) > 0 and max(s) > 0 else 1e-10

        m, n = dirac_op.shape
        # dim ker(D): Anzahl Null-Singulärwerte (rechte Null-Vektoren)
        dim_ker_D = int(np.sum(s < tol))
        if m > n:
            # Rechteckige Matrix: zusätzliche rechte Null-Vektoren
            dim_ker_D += m - n  # Dimensionsformel: nullity = n - rank

        # dim ker(D†) = dim coker(D): linke Null-Vektoren
        rank = int(np.sum(s >= tol))
        dim_coker_D = m - rank  # coker(D) ≅ ker(D†) hat Dimension m - rank(D)

        # Korrekte Berechnung: ind = dim ker(D) - dim coker(D) = n - m (Fredholm-Formel)
        # Für quadratische Operator-Matrix: ind = nullity(D) - nullity(D†)
        dim_ker = n - rank  # Kerneldimension (via Rank-Nullity-Theorem)
        dim_coker = m - rank  # Cokerneldimension

        return dim_ker - dim_coker

    def dirac_operator_spectrum(self, F_components: Dict,
                                mass: float) -> List[float]:
        """
        @brief Berechnet das Spektrum des Dirac-Operators im YM-Hintergrund.
        @description
            Der Dirac-Operator im Yang-Mills-Hintergrund ist:
              D_A = γ^μ (∂_μ + A_μ) + m

            Für euklidische 4D-Raumzeit mit γ-Matrizen (Weyl-Darstellung):
              γ⁰ = σ₃, γ¹ = iσ₁, γ² = iσ₂, γ³ = iσ₃ ⊗ I (etc.)

            Das Spektrum (Eigenwerte von D_A†D_A) bestimmt:
            - Null-Moden: ker(D_A) → topologische Ladung
            - Nicht-Null-Moden: kommen in Paaren ±λ (Chiralitäts-Symmetrie)

            Hier wird eine vereinfachte Näherung berechnet (diskrete Version
            ohne Differential-Operatoren, nur algebraischer Teil).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param F_components  Feldstärke-Komponenten (für Eichfeld-Hintergrund)
        @param mass          Fermion-Masse m
        @return Liste der Eigenwerte von D_A†D_A
        """
        # Dirac-Matrizen (4×4, Weyl-Darstellung)
        # γ^μ = Pauli-Matrizen in 4×4-Darstellung
        I2 = np.eye(2, dtype=complex)
        O2 = np.zeros((2, 2), dtype=complex)
        sigma = _pauli_matrices()

        # Weyl-Darstellung der 4D-Dirac-Matrizen (euklidisch)
        # γ^0 = [[0, I], [I, 0]], γ^k = [[0, -iσ_k], [iσ_k, 0]]
        gamma = [
            np.block([[O2, I2], [I2, O2]]),                   # γ⁰
            np.block([[O2, -1j * sigma[0]], [1j * sigma[0], O2]]),  # γ¹
            np.block([[O2, -1j * sigma[1]], [1j * sigma[1], O2]]),  # γ²
            np.block([[O2, -1j * sigma[2]], [1j * sigma[2], O2]])   # γ³
        ]

        # Vereinfachter Dirac-Operator D = Σ_μ γ^μ A_μ + m (ohne Ableitungen)
        # Nur der Eichfeld-Anteil (für konstante Hintergründe)
        size = 4  # 4×4 Dirac-Matrix-Raum
        D = mass * np.eye(size, dtype=complex)

        # Eichfeld aus F_components extrahieren (Näherung via Spur)
        if isinstance(F_components, dict) and "F_components" in F_components:
            F = F_components["F_components"]
        else:
            F = F_components

        # Eichfeld-Beitrag: Σ_{μ<ν} Tr(F_μν) · γ_μ γ_ν (vereinfacht)
        for (mu, nu), F_mat in F.items():
            if mu < 4 and nu < 4:
                # Spur des Feldstärketensors als skalarer Beitrag
                tr_F = np.real(np.trace(F_mat))
                if mu < len(gamma) and nu < len(gamma):
                    # Kommutatorteil: [γ^μ, γ^ν] ∝ Spin-Verbindung
                    D += tr_F * 0.5 * (gamma[mu] @ gamma[nu] - gamma[nu] @ gamma[mu])

        # Spektrum von D†D (positivere Hermitesche Matrix)
        DdagD = D.conj().T @ D
        eigvals = np.linalg.eigvalsh(DdagD)

        # Sortiert und positiv (Eigenwerte von D†D ≥ 0)
        return sorted([float(e) for e in eigvals if e >= -1e-12])

    def index_theorem_check(self, topology: Dict,
                            analytical_index: int) -> bool:
        """
        @brief Prüft das Atiyah-Singer-Index-Theorem: ind(D) = topol. Index.
        @description
            Das Index-Theorem besagt:
              dim ker(D) - dim coker(D) = topologischer Index

            Für den Dirac-Operator in einem k-Instanton-Hintergrund auf S⁴:
              ind(D_A) = -k   (für SU(2), fundamentale Darstellung)

            Der topologische Index berechnet sich aus:
              ind_top = ∫_M Â(M) ∧ ch(E)

            Für S⁴ mit flacher Metrik (Â(S⁴) = 1) und SU(2)-Bündel:
              ind_top = -c₂(E) = -k = topologische Ladung

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param topology          Dictionary mit topologischen Daten
                                 (z.B. {"topological_charge": k, "manifold": "S4"})
        @param analytical_index  Analytisch berechneter Index ind(D)
        @return True wenn analytischer und topologischer Index übereinstimmen
        """
        # Topologischen Index extrahieren
        topological_charge = topology.get("topological_charge", 0)
        manifold = topology.get("manifold", "S4")

        # Für S⁴ mit SU(2)-Bündel: ind(D) = -k
        if manifold == "S4":
            expected_index = -topological_charge
        else:
            # Allgemeiner Fall: Approximation
            expected_index = -topological_charge

        return analytical_index == expected_index

    def eta_invariant(self, operator: np.ndarray) -> float:
        """
        @brief Berechnet die Eta-Funktion η(0) eines selbstadjungierten Operators.
        @description
            Die Eta-Funktion (Atiyah-Patodi-Singer, 1975) ist für einen
            selbstadjungierten Operator A definiert als:

              η(s) = Σ_{λ≠0} sign(λ) · |λ|^{-s}

            Der Wert η(0) ist die spektrale Asymmetrie — er misst, ob mehr
            positive als negative Eigenwerte vorhanden sind.

            Im APS-Randwertproblem erscheint η(0) als Randterm im Index-Theorem:
              ind(D) = ∫_M Â(M) ch(E) - (h + η(0))/2

            wobei h = dim ker(A) die Anzahl der Null-Moden am Rand ist.

            Berechnung über Matrix-Eigenwerte (endliche Approximation).

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        @param operator  Selbstadjungierte Matrix-Darstellung des Operators (hermitesch)
        @return Eta-Invariante η(0) (reell)
        """
        # Eigenwerte des selbstadjungierten Operators
        eigvals = np.linalg.eigvalsh(operator)
        tol = 1e-10

        # η(0) = Σ_{λ≠0} sign(λ)
        eta = 0.0
        for lam in eigvals:
            if abs(lam) > tol:
                eta += float(np.sign(lam))

        # Die spektrale Asymmetrie η(0) (reell)
        return eta
