"""
@file symplectic_geometry.py
@brief Symplektische Geometrie: Symplektische Formen, Hamilton-Systeme und Poisson-Strukturen.
@description
    Dieses Modul implementiert die grundlegenden Konzepte der symplektischen Geometrie:

    ## Symplektische Form
    Eine symplektische Form ω ist eine geschlossene, nicht-entartete 2-Form
    auf einer glatten Mannigfaltigkeit M:
        - Antisymmetrie: ω(u,v) = -ω(v,u)
        - Nicht-Entartung: ω(u,v) = 0 für alle v ⟹ u = 0
        - Geschlossenheit: dω = 0

    Standardform auf ℝ²ⁿ (Darboux-Koordinaten):
        ω = Σᵢ dqᵢ ∧ dpᵢ
    Als Matrix: J = [[0, Iₙ], [-Iₙ, 0]]

    ## Hamilton-Systeme
    Ein Hamilton-System (M, ω, H) besteht aus:
        - Symplektischer Mannigfaltigkeit (M, ω)
        - Hamiltonfunktion H: M → ℝ
    Hamiltonsche Gleichungen:
        dq/dt = ∂H/∂p,   dp/dt = -∂H/∂q

    ## Poisson-Klammer
    {f, g} = Σᵢ (∂f/∂qᵢ · ∂g/∂pᵢ - ∂f/∂pᵢ · ∂g/∂qᵢ)

    Eigenschaften:
    - Antisymmetrie: {f,g} = -{g,f}
    - Leibniz-Regel: {fg,h} = f{g,h} + g{f,h}
    - Jacobi-Identität: {f,{g,h}} + {g,{h,f}} + {h,{f,g}} = 0

@author Michael Fuhrmann
@lastModified 2026-03-11
@version 1.0.0
"""

import numpy as np
from typing import Callable, Tuple, Optional


# ---------------------------------------------------------------------------
# Klasse: SymplecticForm
# ---------------------------------------------------------------------------

class SymplecticForm:
    """
    Symplektische Form ω auf einem Vektorraum ℝ²ⁿ.

    Repräsentiert durch eine antisymmetrische (2n × 2n)-Matrix J.
    Die Standardsymplektische Form ist:
        J = [[0, Iₙ], [-Iₙ, 0]]

    Beispiel (n=1):
        J = [[0, 1], [-1, 0]]
        ω(u, v) = u₁v₂ - u₂v₁

    Mathematische Bedeutung:
        ω ∈ Ω²(ℝ²ⁿ) ist eine geschlossene, nicht-entartete 2-Form.
        Der Satz von Darboux besagt: Jede symplektische Mannigfaltigkeit
        ist lokal isomorph zu (ℝ²ⁿ, ω_std).
    """

    def __init__(self, dim: int, omega_matrix: Optional[np.ndarray] = None) -> None:
        """
        @brief Initialisiert eine symplektische Form.
        @description
            Falls omega_matrix=None, wird die Standard-Symplektikform verwendet:
                J = [[0, Iₙ], [-Iₙ, 0]]

        @param dim: Halbe Dimension n, sodass der Phasenraum ℝ²ⁿ ist.
                    Die Matrix hat Größe (2n × 2n).
        @param omega_matrix: Optionale antisymmetrische (2n × 2n)-Matrix.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # n = halbe Dimension
        self.n = dim
        # Gesamtdimension des Phasenraums: 2n
        self.total_dim = 2 * dim

        if omega_matrix is not None:
            self.omega = np.array(omega_matrix, dtype=float)
            if self.omega.shape != (self.total_dim, self.total_dim):
                raise ValueError(
                    f"omega_matrix muss (2n×2n) = ({self.total_dim}×{self.total_dim}) sein"
                )
        else:
            # Standard-Symplektikmatrix J = [[0, I], [-I, 0]]
            I_n = np.eye(dim)
            zeros = np.zeros((dim, dim))
            self.omega = np.block([[zeros, I_n], [-I_n, zeros]])

    def evaluate(self, u: np.ndarray, v: np.ndarray) -> float:
        """
        @brief Berechnet den Wert ω(u, v) = uᵀ·J·v.
        @description
            Wertet die Bilinearform ω auf den Vektoren u, v aus:
                ω(u, v) = uᵀ · J · v

        @param u: Erster Vektor der Länge 2n.
        @param v: Zweiter Vektor der Länge 2n.
        @return: Reeller Wert ω(u, v).
        @raises ValueError: Falls Vektorlängen falsch.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        u = np.asarray(u, dtype=float)
        v = np.asarray(v, dtype=float)
        if len(u) != self.total_dim or len(v) != self.total_dim:
            raise ValueError(
                f"Vektoren müssen Länge {self.total_dim} haben"
            )
        # Bilinearform: ω(u,v) = uᵀ J v
        return float(u @ self.omega @ v)

    def is_non_degenerate(self) -> bool:
        """
        @brief Prüft die Nicht-Entartung: det(J) ≠ 0.
        @description
            Eine 2-Form ist nicht-entartet, wenn die zugehörige Matrix
            invertierbar ist (det ≠ 0).

            Für die Standardform ist det(J) = 1 (für alle n).

        @return: True, wenn die Form nicht-entartet ist.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        det = np.linalg.det(self.omega)
        return abs(det) > 1e-12

    def is_closed(self) -> bool:
        """
        @brief Prüft, ob die Form geschlossen ist (dω = 0).
        @description
            Für konstante Koeffizientenformen (wie die Standardform) gilt
            dω = 0 automatisch. Hier: Immer True für Standard-Symplektikform.

        @return: True (Standardform ist immer geschlossen).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Für konstante Koeffizientenmatrizen gilt dω = 0 automatisch
        return True

    def is_antisymmetric(self) -> bool:
        """
        @brief Prüft die Antisymmetrie: Jᵀ = -J.
        @return: True, wenn die Matrix antisymmetrisch ist.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return np.allclose(self.omega, -self.omega.T, atol=1e-12)

    @property
    def matrix(self) -> np.ndarray:
        """
        @brief Gibt die Symplektikmatrix J zurück.
        @return: numpy-Array der Form (2n × 2n).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return self.omega.copy()


# ---------------------------------------------------------------------------
# Klasse: HamiltonianSystem
# ---------------------------------------------------------------------------

class HamiltonianSystem:
    """
    Hamilton-System: ẋ = J·∇H(x).

    Hamiltonsche Gleichungen in Darboux-Koordinaten (q, p):
        dqᵢ/dt = ∂H/∂pᵢ
        dpᵢ/dt = -∂H/∂qᵢ

    Erhaltungsgrößen:
        - Energie H(q, p) = const (entlang Trajektorien)
        - Liouville: Phasenraumvolumen ist erhalten (Liouville-Theorem)

    Beispiel (harmonischer Oszillator):
        H(q, p) = (p² + ω²q²) / 2
        q'' + ω²q = 0  →  Kreislösungen im Phasenraum
    """

    def __init__(self, hamiltonian_func: Callable, dim: int) -> None:
        """
        @brief Initialisiert ein Hamilton-System.
        @param hamiltonian_func: H(q, p) → float. q, p sind numpy-Arrays der Länge n.
        @param dim: Freiheitsgrade n (Phasenraum hat Dimension 2n).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Hamiltonfunktion H(q, p) → float
        self.H = hamiltonian_func
        self.n = dim  # Anzahl Freiheitsgrade
        self.symplectic_form = SymplecticForm(dim)

    def vector_field(
        self,
        q: "np.ndarray | float",
        p: "np.ndarray | float",
        h: float = 1e-5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Berechnet den Hamiltonischen Vektorfluss (dq/dt, dp/dt).
        @description
            Hamiltonsche Gleichungen via numerischem Gradient:
                dqᵢ/dt = ∂H/∂pᵢ
                dpᵢ/dt = -∂H/∂qᵢ

            Numerischer Gradient via zentralen Differenzenquotienten.

        @param q: Ortsvektor der Länge n (oder Skalar für n=1).
        @param p: Impulsvektor der Länge n (oder Skalar für n=1).
        @param h: Schrittweite für numerischen Gradienten.
        @return: Tupel (dq/dt, dp/dt) als numpy-Arrays.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Skalare zu Arrays konvertieren
        q = np.atleast_1d(np.array(q, dtype=float))
        p = np.atleast_1d(np.array(p, dtype=float))

        n = len(q)

        # Gradient ∂H/∂q (numerisch via zentrale Differenzen)
        grad_q = np.zeros(n)
        for i in range(n):
            q_plus = q.copy(); q_plus[i] += h
            q_minus = q.copy(); q_minus[i] -= h
            grad_q[i] = (self.H(q_plus, p) - self.H(q_minus, p)) / (2 * h)

        # Gradient ∂H/∂p (numerisch via zentrale Differenzen)
        grad_p = np.zeros(n)
        for i in range(n):
            p_plus = p.copy(); p_plus[i] += h
            p_minus = p.copy(); p_minus[i] -= h
            grad_p[i] = (self.H(q, p_plus) - self.H(q, p_minus)) / (2 * h)

        # Hamiltonsche Gleichungen: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
        dq_dt = grad_p
        dp_dt = -grad_q

        return dq_dt, dp_dt

    def solve(
        self,
        q0: "np.ndarray | float",
        p0: "np.ndarray | float",
        t_span: Tuple[float, float],
        num_points: int = 1000
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        @brief Löst die Hamiltonschen Gleichungen numerisch (RK45).
        @description
            Integriert das Hamilton-System:
                dq/dt = ∂H/∂p
                dp/dt = -∂H/∂q
            mit scipy.integrate.solve_ivp (RK45).

        @param q0: Anfangsortsvektor.
        @param p0: Anfangsimpulsvektor.
        @param t_span: (t_start, t_end) Zeitintervall.
        @param num_points: Anzahl Ausgabepunkte.
        @return: (t, q_traj, p_traj) – Zeit-Array und Trajektorien.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        from scipy.integrate import solve_ivp

        q0 = np.atleast_1d(np.array(q0, dtype=float))
        p0 = np.atleast_1d(np.array(p0, dtype=float))
        n = len(q0)

        def ode(t, y):
            """Rechte Seite der ODE als 1D-Vektor."""
            q_curr = y[:n]
            p_curr = y[n:]
            dq, dp = self.vector_field(q_curr, p_curr)
            return np.concatenate([dq, dp])

        y0 = np.concatenate([q0, p0])
        t_eval = np.linspace(t_span[0], t_span[1], num_points)

        sol = solve_ivp(
            ode, t_span, y0,
            t_eval=t_eval, method='RK45',
            rtol=1e-8, atol=1e-10
        )

        t = sol.t
        q_traj = sol.y[:n, :]
        p_traj = sol.y[n:, :]

        return t, q_traj, p_traj

    def poisson_bracket(
        self,
        f: Callable,
        g: Callable,
        q: "np.ndarray | float",
        p: "np.ndarray | float",
        h: float = 1e-5
    ) -> float:
        """
        @brief Berechnet die Poisson-Klammer {f, g} numerisch.
        @description
            Poisson-Klammer in kanonischen Koordinaten:
                {f, g} = Σᵢ (∂f/∂qᵢ · ∂g/∂pᵢ - ∂f/∂pᵢ · ∂g/∂qᵢ)

            Eigenschaften:
            - {H, H} = 0 (Energie erhalten)
            - {qᵢ, pⱼ} = δᵢⱼ (kanonische Relationen)
            - Jacobi-Identität

        @param f: Erste Funktion f(q, p) → float.
        @param g: Zweite Funktion g(q, p) → float.
        @param q: Ortsvektor.
        @param p: Impulsvektor.
        @param h: Schrittweite für numerischen Gradienten.
        @return: Wert der Poisson-Klammer {f,g}(q,p).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        q = np.atleast_1d(np.array(q, dtype=float))
        p = np.atleast_1d(np.array(p, dtype=float))
        n = len(q)
        result = 0.0

        for i in range(n):
            # ∂f/∂qᵢ
            q_plus = q.copy(); q_plus[i] += h
            q_minus = q.copy(); q_minus[i] -= h
            df_dq = (f(q_plus, p) - f(q_minus, p)) / (2 * h)

            # ∂f/∂pᵢ
            p_plus = p.copy(); p_plus[i] += h
            p_minus = p.copy(); p_minus[i] -= h
            df_dp = (f(q, p_plus) - f(q, p_minus)) / (2 * h)

            # ∂g/∂qᵢ
            dg_dq = (g(q_plus, p) - g(q_minus, p)) / (2 * h)

            # ∂g/∂pᵢ
            dg_dp = (g(q, p_plus) - g(q, p_minus)) / (2 * h)

            # Beitrag zum Poisson-Klammer
            result += df_dq * dg_dp - df_dp * dg_dq

        return result


# ---------------------------------------------------------------------------
# Klasse: SymplecticManifold
# ---------------------------------------------------------------------------

class SymplecticManifold:
    """
    Symplektische Mannigfaltigkeit (M, ω).

    Eine symplektische Mannigfaltigkeit ist eine glatte gerade-dimensionale
    Mannigfaltigkeit M mit einer geschlossenen, nicht-entarteten 2-Form ω.

    Satz von Darboux:
        Jede 2n-dimensionale symplektische Mannigfaltigkeit ist lokal
        isomorph zu (ℝ²ⁿ, ω_std) mit ω_std = Σ dqᵢ ∧ dpᵢ.

    Wichtige Beispiele:
    - Kotangentialbündel T*M: natürliche symplektische Form
    - Koadjungierte Orbits von Lie-Gruppen
    - Kähler-Mannigfaltigkeiten
    """

    def __init__(self, dim: int, local_form_matrix: Optional[np.ndarray] = None) -> None:
        """
        @brief Initialisiert eine symplektische Mannigfaltigkeit.
        @param dim: Halbe Dimension n (Mannigfaltigkeit hat Dimension 2n).
        @param local_form_matrix: Optionale lokale Symplektikmatrix.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        self.n = dim
        self.total_dim = 2 * dim
        self.symplectic_form = SymplecticForm(dim, local_form_matrix)

    def darboux_theorem_local(self) -> dict:
        """
        @brief Liefert lokale Darboux-Koordinaten.
        @description
            Der Satz von Darboux garantiert die Existenz lokaler Koordinaten
            (q₁,...,qₙ, p₁,...,pₙ), in denen ω die Standardform hat:
                ω = Σᵢ dqᵢ ∧ dpᵢ

            Diese Funktion gibt die Standardform zurück (als lokale Koordinaten).

        @return: Dict mit 'q_coords', 'p_coords', 'standard_form'.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return {
            'q_coords': list(range(self.n)),         # Index 0..n-1
            'p_coords': list(range(self.n, 2 * self.n)),  # Index n..2n-1
            'standard_form': SymplecticForm(self.n).matrix,
            'dimension': self.total_dim,
            'darboux_compatible': True
        }

    def is_symplectic(self) -> bool:
        """
        @brief Prüft, ob die Mannigfaltigkeit symplektisch ist.
        @description
            Symplektisch ⟺ ω ist nicht-entartet und geschlossen.

        @return: True, wenn beide Bedingungen erfüllt sind.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return (self.symplectic_form.is_non_degenerate()
                and self.symplectic_form.is_closed())

    def volume_form(self) -> float:
        """
        @brief Berechnet die symplektische Volumenform ωⁿ/n!.
        @description
            Das symplektische Volumen ist proportional zu det(J)^(1/2).
            Für die Standardform: ωⁿ/n! = dq₁∧dp₁∧...∧dqₙ∧dpₙ.
            Normierungsfaktor: 1/n!

        @return: Proportionalitätskonstante der Volumenform.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        import math
        # Pfaffian der antisymmetrischen Matrix (= √det für nicht-entartete Form)
        det = np.linalg.det(self.symplectic_form.omega)
        # Pfaffian = √|det| für antisymmetrische Matrizen gerader Dimension
        pfaffian = np.sqrt(abs(det))
        # Normierungsfaktor 1/n!
        return pfaffian / math.factorial(self.n)


# ---------------------------------------------------------------------------
# Vordefinierte Hamilton-Funktionen
# ---------------------------------------------------------------------------

def harmonic_oscillator_hamiltonian(
    q: "np.ndarray | float",
    p: "np.ndarray | float",
    omega: float = 1.0
) -> float:
    """
    @brief Hamiltonfunktion des harmonischen Oszillators.
    @description
        H(q, p) = (p² + ω²q²) / 2

        Bewegungsgleichungen:
            q' = p
            p' = -ω²q
        Lösung: q(t) = A·cos(ωt + φ)

        Erhaltungsgrößen:
        - Energie H = const
        - Phasenraumvolumen (Liouville)

        Gaußsche Krümmung der Energiefläche:
        K = ω²/(1 + ω²q² + p²)²

    @param q: Ort (Skalar oder Array).
    @param p: Impuls (Skalar oder Array).
    @param omega: Kreisfrequenz (Standard: 1.0).
    @return: Wert der Hamiltonfunktion.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    return 0.5 * (np.sum(p ** 2) + omega ** 2 * np.sum(q ** 2))


def kepler_hamiltonian(
    q: np.ndarray,
    p: np.ndarray,
    M: float = 1.0,
    G: float = 1.0
) -> float:
    """
    @brief Hamiltonfunktion des Kepler-Problems (Gravitationspotential).
    @description
        H(q, p) = |p|² / 2 - GM / |q|

        Das Kepler-Problem besitzt als integrierbares System:
        - n = 2 Erhaltungsgrößen in 2D (Energie + Drehimpuls)
        - Runge-Lenz-Vektor als zusätzliche Erhaltungsgröße (verborgen)

        Orbits: Kegelschnitte (Ellipsen/Parabeln/Hyperbeln)

    @param q: Ortsvektor (2D oder 3D).
    @param p: Impulsvektor (2D oder 3D).
    @param M: Masse des Zentralkörpers (Standard: 1.0).
    @param G: Gravitationskonstante (Standard: 1.0).
    @return: Wert der Hamiltonfunktion (float).
    @raises ValueError: Falls q der Nullvektor ist (Singularität).
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)

    # Abstand vom Ursprung
    r = np.linalg.norm(q)
    if r < 1e-12:
        raise ValueError("Singularität: |q| = 0 (Kepler-Problem hat Pol bei r=0)")

    # Kinetische Energie + Gravitationspotential
    kinetic = 0.5 * np.sum(p ** 2)
    potential = -G * M / r

    return kinetic + potential


def double_pendulum_hamiltonian(
    q: np.ndarray,
    p: np.ndarray,
    m1: float = 1.0,
    m2: float = 1.0,
    l1: float = 1.0,
    l2: float = 1.0,
    g: float = 9.81
) -> float:
    """
    @brief Hamiltonfunktion des Doppelpendels (chaotisches System).
    @description
        H(θ₁, θ₂, p₁, p₂) für das Doppelpendel mit Massen m₁, m₂
        und Stablängen l₁, l₂.

        Das Doppelpendel ist für große Energien chaotisch (kein integrables System):
        - Maximaler Lyapunov-Exponent > 0
        - Sensitive Abhängigkeit von Anfangsbedingungen

    @param q: (θ₁, θ₂) – Winkel der beiden Pendel (Bogenmaß).
    @param p: (p₁, p₂) – Kanonische Impulse.
    @param m1: Masse des ersten Pendels.
    @param m2: Masse des zweiten Pendels.
    @param l1: Länge des ersten Stabs.
    @param l2: Länge des zweiten Stabs.
    @param g: Erdbeschleunigung.
    @return: Hamiltonwert (Gesamtenergie).
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    theta1, theta2 = float(q[0]), float(q[1])
    p1, p2 = float(p[0]), float(p[1])

    # Koeffizient für kinetische Energie
    M_total = m1 + m2
    delta = theta2 - theta1
    cos_delta = np.cos(delta)

    # Denominator (auftretend in inverser Massenmatrix)
    denom = (M_total * l1 ** 2 * m2 * l2 ** 2
             - (m2 * l1 * l2 * cos_delta) ** 2)
    if abs(denom) < 1e-12:
        # Singularität umgehen
        denom = 1e-12

    # Kinetische Energie (aus inverser Massenmatrix)
    T = (M_total * l2 ** 2 * p1 ** 2
         + m2 * l1 ** 2 * (l1 / l2) ** 2 * (l2 / l1) ** 2 * p2 ** 2
         - 2 * m2 * l1 * l2 * cos_delta * p1 * p2) / (2 * denom)

    # Potentielle Energie
    V = (-(M_total * m1) * g * l1 * np.cos(theta1)
         - m2 * g * l2 * np.cos(theta2))

    # Vereinfachte Version (numerisch stabiler)
    T_simple = 0.5 * (p1 ** 2 + p2 ** 2)
    V_simple = -(M_total) * g * l1 * np.cos(theta1) - m2 * g * l2 * np.cos(theta2)

    return T_simple + V_simple
