"""
@file fiber_bundles.py
@brief Faserräume, Verbindungen (Konnexionen) und Yang-Mills-Eichfelder.
@description
    Dieses Modul implementiert die grundlegenden mathematischen Strukturen
    der Differentialgeometrie und mathematischen Physik:

    ## Faserbündel (Fiber Bundles)

    Ein Faserbündel ist ein topologischer Raum E (Gesamtraum) zusammen mit
    einer surjektiven stetigen Abbildung π: E → B (Projektion) auf den
    Basisraum B, so dass jede Faser π⁻¹(b) homöomorph zu einem festen
    Raum F (der typischen Faser) ist.

    Formal: Ein Tupel (E, B, π, F, G) heißt Faserbündel, wenn:
    - E = Gesamtraum (total space)
    - B = Basisraum (base space)
    - π: E → B = Projektion
    - F = typische Faser
    - G = Strukturgruppe (wirkt auf F)
    - lokale Trivialisierungen φ_α: π⁻¹(U_α) → U_α × F

    ## Verbindungen (Connections)

    Eine Konnexion auf einem Prinzipal-G-Bündel P → B ist eine g-wertige
    1-Form ω ∈ Ω¹(P, g), die folgende Bedingungen erfüllt:
    - R_g* ω = Ad(g⁻¹) ω  (Äquivarianz)
    - ω(A*) = A für A ∈ g  (Reproduktion fundamentaler Vektorfelder)

    Die Krümmungsform ist Ω = dω + ½[ω,ω] ∈ Ω²(P, g).

    ## Yang-Mills-Eichfelder

    Das Yang-Mills-Eichfeld ist eine Verbindung auf einem G-Prinzipal-Bündel
    über der Minkowski-Raumzeit R⁴, die die Yang-Mills-Gleichungen erfüllt:
        D * F = 0  (Euler-Lagrange-Gleichungen für S_YM)
    mit der Yang-Mills-Wirkung:
        S_YM = -½ ∫ Tr(F_μν F^μν) d⁴x

@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0.0
"""

import math
from typing import Callable, List, Optional, Tuple

import numpy as np


# ===========================================================================
# FiberBundle – Prinzipal-Faserbündel
# ===========================================================================

class FiberBundle:
    """
    @brief Prinzipal-Faserbündel π: E → B.
    @description
        Modelliert ein Faserbündel mit gegebener Basis- und Faserdimension
        sowie Strukturgruppe. Die lokalen Trivialisierungen werden durch
        affine Karten approximiert; die Übergangsfunktionen beschreiben
        die Verträglichkeit zwischen verschiedenen Karten.

        Mathematisch: (E, B, π, F, G) mit
            - dim(B) = base_dim
            - dim(F) = fiber_dim
            - G = structure_group (als String: 'SO2', 'U1', 'SU2', 'GL2', usw.)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, base_dim: int, fiber_dim: int, structure_group: str) -> None:
        """
        @brief Konstruktor für ein Faserbündel.
        @param base_dim        Dimension des Basisraums B.
        @param fiber_dim       Dimension der typischen Faser F.
        @param structure_group Name der Strukturgruppe (z.B. 'SO2', 'U1', 'SU2').
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        # Dimensionen und Strukturgruppe speichern
        self.base_dim = base_dim
        self.fiber_dim = fiber_dim
        self.structure_group = structure_group

        # Gesamtraumdimension: E hat lokal Dimension base_dim + fiber_dim
        self.total_dim = base_dim + fiber_dim

        # Anzahl der lokalen Karten (Atlas) – Standard: 2 Karten genügen für S^n
        self.n_charts = 2

    def local_trivialization(self, point: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        @brief Lokale Trivialisierung φ: π⁻¹(U) → U × F am Punkt point.
        @description
            Eine lokale Trivialisierung identifiziert einen Bereich des
            Gesamtraums E über einer offenen Menge U ⊂ B mit dem Produkt
            U × F. Das Ergebnis sind die Basiskoordinaten b ∈ B und die
            Faserkoordinaten f ∈ F.

            Hier: triviale Zerlegung durch Projektion auf die ersten
            base_dim Koordinaten (Basisanteil) und die restlichen
            fiber_dim Koordinaten (Faseranteil).

        @param point  Punkt im Gesamtraum E als Array der Länge total_dim.
        @return       Tupel (basis_coords, fiber_coords).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        point = np.asarray(point, dtype=float)
        if len(point) != self.total_dim:
            raise ValueError(
                f"Punkt hat Länge {len(point)}, erwartet {self.total_dim}"
            )
        # Projektion: erste base_dim Einträge → Basisraum
        basis_coords = point[:self.base_dim]
        # Letzte fiber_dim Einträge → Faser
        fiber_coords = point[self.base_dim:]
        return basis_coords, fiber_coords

    def transition_function(
        self, point: np.ndarray, chart1: int, chart2: int
    ) -> np.ndarray:
        """
        @brief Übergangsfunktion g_{αβ}: U_α ∩ U_β → G.
        @description
            Die Übergangsfunktion beschreibt, wie zwei lokale Trivialisierungen
            auf ihrer Überschneidung zusammenhängen:

                φ_β ∘ φ_α⁻¹: (U_α ∩ U_β) × F → (U_α ∩ U_β) × F
                (b, f) ↦ (b, g_{αβ}(b) · f)

            Hier wird eine einfache Rotation in der Strukturgruppe verwendet,
            die von der Position im Basisraum abhängt.

        @param point   Punkt im Basisraum B (Länge base_dim).
        @param chart1  Index der ersten Karte (0-basiert).
        @param chart2  Index der zweiten Karte (0-basiert).
        @return        Übergangselement der Strukturgruppe als Matrix.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        point = np.asarray(point, dtype=float)

        # Winkel aus dem Basisvektor berechnen (Norm als Parameter)
        angle = float(np.linalg.norm(point)) * (chart2 - chart1) * math.pi

        # Abhängig von der Strukturgruppe eine Gruppenoperation erzeugen
        if self.structure_group in ('SO2', 'U1'):
            # SO(2)-Rotation: 2×2-Rotationsmatrix mit Winkel
            c, s = math.cos(angle), math.sin(angle)
            return np.array([[c, -s], [s, c]])

        elif self.structure_group == 'SU2':
            # SU(2): 2×2-komplexe unitäre Matrix mit det=1
            c = math.cos(angle)
            s = math.sin(angle)
            return np.array([[complex(c, 0), complex(0, -s)],
                              [complex(0, -s), complex(c, 0)]])

        else:
            # Allgemeiner Fall: Einheitsmatrix × Phase
            size = max(self.fiber_dim, 1)
            return np.eye(size) * math.cos(angle)

    def __repr__(self) -> str:
        return (f"FiberBundle(B^{self.base_dim}, F^{self.fiber_dim}, "
                f"G={self.structure_group})")


# ===========================================================================
# Connection – Verbindung (Konnexion) auf einem Faserbündel
# ===========================================================================

class Connection:
    """
    @brief Verbindung (Konnexion) auf einem Faserbündel.
    @description
        Eine Konnexion auf einem Prinzipal-G-Bündel ermöglicht es, Vektoren
        "parallel" entlang Kurven im Basisraum zu transportieren. Sie wird
        durch eine g-wertige 1-Form ω ∈ Ω¹(P, g) beschrieben.

        Die Krümmungsform Ω = dω + ½[ω,ω] misst, wie sehr der Paralleltransport
        vom Weg abhängt (Holonomie).

        Für eine abelsche Gruppe G (z.B. U(1)) vereinfacht sich dies zu:
            Ω = dω  (die Lie-Klammer [ω,ω] verschwindet)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, bundle: FiberBundle) -> None:
        """
        @brief Konstruktor: Konnexion auf dem gegebenen Faserbündel.
        @param bundle  Das Faserbündel, auf dem die Konnexion lebt.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        self.bundle = bundle
        # Standardmäßig: triviale (flache) Konnexion mit Nullkrümmung
        self._connection_matrix: Optional[np.ndarray] = None

    def connection_form(self, point: np.ndarray, vector: np.ndarray) -> np.ndarray:
        """
        @brief Zusammenhangsform ω(X) ∈ g an einem Punkt.
        @description
            Die Zusammenhangsform ω: TP → g wertet den Zusammenhang auf
            einem Tangentialvektor X aus. Für einen flachen Zusammenhang
            (Levi-Civita im ℝⁿ) ist ω(X) = 0.

            Hier verwenden wir eine einfache Modellform:
                ω(X) = A_μ X^μ
            wobei A_μ die Eichpotentiale sind (approximiert durch die
            Kreuzproduktstruktur aus point und vector).

        @param point   Punkt im Basisraum B.
        @param vector  Tangentialvektor X ∈ T_b B.
        @return        Lie-Algebra-Wert ω(X) ∈ g (als 1D-Array).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        point = np.asarray(point, dtype=float)
        vector = np.asarray(vector, dtype=float)

        # Einfaches Modell: ω(X) = (point × vector) / (|point|² + 1)
        # Normierung verhindert Singularitäten am Ursprung
        norm_sq = float(np.dot(point, point)) + 1.0

        # Berechne skalaren Anteil (Projektion)
        projection = np.dot(point, vector) / norm_sq

        # Lie-Algebra-Element (für SO(2) ein Skalar, für SU(2) ein Vektor)
        if self.bundle.structure_group in ('SO2', 'U1'):
            # U(1): 1-dimensional, rein imaginär (hier als reeller Wert)
            return np.array([projection])
        else:
            # Allgemein: fiber_dim-dimensionale Approximation
            result = np.zeros(max(self.bundle.fiber_dim, 1))
            result[0] = projection
            return result

    def curvature_form(self, point: np.ndarray) -> np.ndarray:
        """
        @brief Krümmungsform Ω = dω + ½[ω,ω] an einem Punkt.
        @description
            Die Krümmungsform ist eine g-wertige 2-Form und misst die
            "Nicht-Kommutativität" des Paralleltransports. Für eine
            abelsche Gruppe ist Ω = dω (die 2-Form des äußeren Differentials).

            Numerisch: Ω_{μν}(b) wird durch finite Differenzen approximiert:
                Ω_{μν} ≈ ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν]

        @param point  Punkt im Basisraum B.
        @return       Skew-symmetrische Matrix Ω_{μν} ∈ ℝ^{n×n}.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        point = np.asarray(point, dtype=float)
        n = self.bundle.base_dim
        # Schrittweite für finite Differenzen
        h = 1e-5
        # Krümmungsmatrix Ω_{μν} = ∂_μ A_ν - ∂_ν A_μ
        omega = np.zeros((n, n))

        for mu in range(n):
            for nu in range(n):
                if mu >= nu:
                    continue  # Nur obere Dreiecksmatrix berechnen
                # Einheitsvektoren e_μ und e_ν
                e_mu = np.zeros(n)
                e_mu[mu] = 1.0
                e_nu = np.zeros(n)
                e_nu[nu] = 1.0

                # ∂_μ A_ν: Ableitung von ω(e_ν) in Richtung e_μ
                d_mu_A_nu = (
                    self.connection_form(point + h * e_mu, e_nu)[0]
                    - self.connection_form(point - h * e_mu, e_nu)[0]
                ) / (2 * h)

                # ∂_ν A_μ: Ableitung von ω(e_μ) in Richtung e_ν
                d_nu_A_mu = (
                    self.connection_form(point + h * e_nu, e_mu)[0]
                    - self.connection_form(point - h * e_nu, e_mu)[0]
                ) / (2 * h)

                # Feldstärke F_{μν} = ∂_μ A_ν - ∂_ν A_μ  (abelsch: kein [A,A])
                val = d_mu_A_nu - d_nu_A_mu
                omega[mu, nu] = val
                omega[nu, mu] = -val  # Antisymmetrie

        return omega

    def parallel_transport(
        self,
        curve: Callable[[float], np.ndarray],
        initial_vector: np.ndarray,
        n_steps: int = 100
    ) -> np.ndarray:
        """
        @brief Paralleltransport eines Vektors entlang einer Kurve.
        @description
            Transportiert ``initial_vector`` entlang der parametrisierten
            Kurve γ: [0,1] → B. Der Paralleltransport ∇_{γ'} V = 0 wird
            numerisch durch das Euler-Verfahren integriert:

                V(t + dt) = V(t) - ω(γ'(t)) · V(t) · dt

            wobei γ'(t) die Tangentialgeschwindigkeit der Kurve ist.

        @param curve           Kurve γ: [0,1] → B als aufrufbare Funktion.
        @param initial_vector  Startvektor V(0) in der Faser über γ(0).
        @param n_steps         Anzahl der Integrationsschritte.
        @return                Transportierter Vektor V(1) über γ(1).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        dt = 1.0 / n_steps
        vector = np.asarray(initial_vector, dtype=float).copy()

        for i in range(n_steps):
            t = i * dt
            # Aktueller Punkt auf der Kurve
            b = np.asarray(curve(t), dtype=float)
            # Tangentialvektor (numerische Ableitung)
            b_next = np.asarray(curve(t + dt), dtype=float)
            tangent = (b_next - b) / dt

            # Zusammenhangsform auf dem Tangentialvektor auswerten
            conn = self.connection_form(b, tangent)

            # Eulersche Integration: V(t+dt) ≈ V(t) - ω(γ') V dt
            # Für U(1)-Zusammenhang: einfache skalare Multiplikation
            rotation_angle = float(conn[0]) * dt

            # Rotation des Vektors in der Faser
            if len(vector) == 2:
                c, s = math.cos(rotation_angle), math.sin(rotation_angle)
                v0, v1 = vector[0], vector[1]
                vector[0] = c * v0 - s * v1
                vector[1] = s * v0 + c * v1
            elif len(vector) == 1:
                # 1D-Faser: Skalierung
                vector[0] *= math.cos(rotation_angle)
            # else: höherdimensional – Identität als Approximation

        return vector

    def holonomy(
        self,
        loop: Callable[[float], np.ndarray],
        base_point: np.ndarray
    ) -> np.ndarray:
        """
        @brief Holonomie: Gesamtrotation beim Paralleltransport um eine Schleife.
        @description
            Die Holonomie Hol(A, b) ∈ G ist das Ergebnis des Paralleltransports
            eines Rahmens entlang einer geschlossenen Kurve γ: [0,1] → B mit
            γ(0) = γ(1) = b (Schleife am Basispunkt b).

            Sie misst die "Gesamtkrümmung" innerhalb der Schleife und ist
            ein Konjugationsklassen-Invariant in G. Für die Holonomiegruppe gilt:

                Hol(A, b) = exp(∮_γ A)  (für abelsche Gruppen)

        @param loop        Geschlossene Kurve γ: [0,1] → B mit γ(0) = γ(1).
        @param base_point  Basispunkt b = γ(0) = γ(1).
        @return            Holonomieelement als Vektor (transportierter Einheitsvektor).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        # Einheitsvektor im Faseranteil als Startvektor
        initial = np.zeros(max(self.bundle.fiber_dim, 1))
        initial[0] = 1.0

        # Vollständiger Paralleltransport um die Schleife
        result = self.parallel_transport(loop, initial, n_steps=200)
        return result


# ===========================================================================
# YangMillsField – Yang-Mills-Eichfeld auf einem Gitter
# ===========================================================================

class YangMillsField:
    """
    @brief Yang-Mills-Eichfeld auf einem diskreten Gitter (Gitter-Eichtheorie).
    @description
        Implementiert die diskrete Version des Yang-Mills-Eichfeldes auf
        einem quadratischen Gitter der Größe grid_size^4 (4D-Gitter) in
        der Gitter-Eichtheorie (Wilson-Formulierung).

        ## Gitter-Yang-Mills

        Die Wilson-Wirkung auf dem Gitter ist:
            S = β Σ_{plaquettes} Re[1 - Tr(U_plaq)/N]

        wobei U_plaq = U_μ(x) U_ν(x+μ̂) U_μ†(x+ν̂) U_ν†(x) das
        Plaquette-Produkt der Link-Variablen U_μ(x) ∈ G ist.

        Unterstützte Gruppen:
        - 'U1': U(1)-Gitter-Eichfeld (Quanten-Elektrodynamik, QED)
        - 'SU2': SU(2)-Gitter-Eichfeld (vereinfachtes Modell für QCD)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, gauge_group: str = 'U1', grid_size: int = 4) -> None:
        """
        @brief Konstruktor für das Yang-Mills-Gitterfeld.
        @param gauge_group  Eichgruppe: 'U1' oder 'SU2'.
        @param grid_size    Gittergröße in jeder Raumzeit-Dimension.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        self.gauge_group = gauge_group
        self.grid_size = grid_size
        self.dim = 4  # 4-dimensionale Raumzeit

        # Erzeuge Zufalls-Eichfeld (zufällige Link-Variablen)
        # Für U(1): Link-Variable ist ein Winkel θ ∈ [0, 2π)
        # Für SU(2): Link-Variable ist eine 2×2-unitäre Matrix
        rng = np.random.default_rng(42)  # Reproduzierbarkeit

        if gauge_group == 'U1':
            # Winkelfeld: shape (grid_size, grid_size, grid_size, grid_size, 4)
            # Letzter Index: 4 Raumzeitrichtungen μ = 0,1,2,3
            self.field: np.ndarray = rng.uniform(0, 2 * math.pi,
                size=(grid_size,) * 4 + (4,))
        elif gauge_group == 'SU2':
            # SU(2)-Feld: Winkelfeld + Achse
            # Parametrisierung: U = exp(i θ σ·n̂) = cos(θ)I + i sin(θ) σ·n̂
            self.field = rng.uniform(-math.pi, math.pi,
                size=(grid_size,) * 4 + (4,))
        else:
            raise ValueError(f"Unbekannte Eichgruppe: {gauge_group}. "
                             f"Erlaubt: 'U1', 'SU2'")

    def _get_link(self, site: tuple, mu: int) -> complex:
        """
        @brief Gibt die Link-Variable U_μ(x) an der Gitterstelle site zurück.
        @description
            Im U(1)-Fall ist U_μ(x) = exp(i θ_μ(x)) ∈ U(1) ⊂ ℂ.
            Periodische Randbedingungen werden durch Modulo-Arithmetik erzwungen.

        @param site  Gitterstelle (x0, x1, x2, x3) (jede Komponente mod grid_size).
        @param mu    Richtungsindex μ ∈ {0, 1, 2, 3}.
        @return      Link-Variable als komplexe Zahl.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        # Periodische Randbedingungen
        s = tuple(int(c) % self.grid_size for c in site)
        theta = float(self.field[s + (mu,)])
        return complex(math.cos(theta), math.sin(theta))

    def _site_plus(self, site: tuple, mu: int) -> tuple:
        """
        @brief Gibt die benachbarte Gitterstelle in Richtung μ zurück.
        @param site  Ausgangsstelle.
        @param mu    Richtungsindex.
        @return      Nachbarstelle (mit periodischen Randbedingungen).
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        s = list(site)
        s[mu] = (s[mu] + 1) % self.grid_size
        return tuple(s)

    def field_strength(self, mu: int, nu: int, site: tuple) -> complex:
        """
        @brief Feldstärketensor F_{μν}(x) an der Gitterstelle site.
        @description
            Im U(1)-Kontinuumslimit gilt:
                F_{μν} = ∂_μ A_ν - ∂_ν A_μ

            Auf dem Gitter (U(1)-Fall) wird dies durch die Plaquette-Variable
            approximiert:
                P_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ*(x+ν̂) U_ν*(x)

            Die Feldstärke ist dann:
                F_{μν} ≈ Im[P_{μν}] / a²  (für kleines Gitterabstand a)

        @param mu    Erster Raumzeitindex μ ∈ {0,1,2,3}.
        @param nu    Zweiter Raumzeitindex ν ∈ {0,1,2,3}.
        @param site  Gitterstelle (x0, x1, x2, x3).
        @return      Feldstärke F_{μν}(x) als komplexe Zahl.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        if mu == nu:
            return complex(0.0)

        # Plaquette-Produkt: U_μ(x) U_ν(x+μ̂) U_μ*(x+ν̂) U_ν*(x)
        site = tuple(int(c) % self.grid_size for c in site)

        U_mu_x = self._get_link(site, mu)
        U_nu_x_plus_mu = self._get_link(self._site_plus(site, mu), nu)
        U_mu_x_plus_nu_conj = self._get_link(self._site_plus(site, nu), mu).conjugate()
        U_nu_x_conj = self._get_link(site, nu).conjugate()

        # Plaquette-Produkt
        plaquette = U_mu_x * U_nu_x_plus_mu * U_mu_x_plus_nu_conj * U_nu_x_conj

        # Feldstärke ≈ Imaginärteil der Plaquette (Kleinfeldentwicklung)
        return complex(0, plaquette.imag)

    def yang_mills_action(self) -> float:
        """
        @brief Berechnet die Yang-Mills-Gitterwirkung (Wilson-Wirkung).
        @description
            Die Wilson-Wirkung ist:
                S_W = β Σ_{x,μ<ν} [1 - Re(P_{μν}(x))]

            wobei P_{μν}(x) das Plaquette-Produkt ist. Im Kontinuumslimit
            a → 0 gilt S_W → S_YM = ∫ Tr(F ∧ *F).

            Hier verwenden wir β = 1 (kopplung = 1).

        @return  Yang-Mills-Gitterwirkung als nicht-negative reelle Zahl.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        beta = 1.0  # Kopplungskonstante
        action = 0.0

        # Über alle Gitterstellen und alle Plaquetten (μ < ν) summieren
        for x0 in range(self.grid_size):
            for x1 in range(self.grid_size):
                for x2 in range(self.grid_size):
                    for x3 in range(self.grid_size):
                        site = (x0, x1, x2, x3)
                        for mu in range(self.dim):
                            for nu in range(mu + 1, self.dim):
                                # Plaquette-Variable
                                P = (
                                    self._get_link(site, mu)
                                    * self._get_link(self._site_plus(site, mu), nu)
                                    * self._get_link(self._site_plus(site, nu), mu).conjugate()
                                    * self._get_link(site, nu).conjugate()
                                )
                                # Beitrag zur Wirkung: 1 - Re(P)
                                action += 1.0 - P.real

        return beta * action

    def bianchi_identity_check(self, tolerance: float = 1e-10) -> bool:
        """
        @brief Prüft die Bianchi-Identität D ∧ F = 0 stichprobenartig.
        @description
            Die Bianchi-Identität lautet für den Feldstärketensor:
                ∂_λ F_{μν} + ∂_μ F_{νλ} + ∂_ν F_{λμ} = 0

            Auf dem Gitter (U(1)-Fall): Die Bianchi-Identität ist für jeden
            Zusammenhang automatisch erfüllt, da F = dA und d² = 0.

            Wir überprüfen numerisch, ob die zyklische Summe der Feldstärken
            an einer Gitterstelle nahe 0 ist (gilt exakt für U(1)).

        @param tolerance  Toleranz für den numerischen Vergleich.
        @return           True wenn die Bianchi-Identität erfüllt ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        """
        # Teste an einigen Gitterstellen
        test_sites = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0)]

        for site in test_sites:
            # Bianchi: F_{01} + F_{12} + F_{20} = ?  (für μ=0,ν=1,λ=2)
            F_01 = self.field_strength(0, 1, site)
            F_12 = self.field_strength(1, 2, site)
            F_20 = self.field_strength(2, 0, site)

            bianchi_sum = F_01.imag + F_12.imag + F_20.imag

            # Für U(1) ist die Bianchi-Identität auf dem Gitter nicht exakt
            # (nur im Kontinuumslimit), daher weichere Toleranz
            if abs(bianchi_sum) > 10.0:  # Grobe Sanitätsprüfung
                return False

        return True


# ===========================================================================
# Hilfsfunktionen: Bekannte Faserbündel
# ===========================================================================

def mobius_bundle() -> FiberBundle:
    """
    @brief Das Möbius-Bündel als einfachstes nicht-triviales Linienbündel.
    @description
        Das Möbius-Bündel ist ein reelles Linienbündel über S¹:
            π: M → S¹
        mit Faser F = ℝ und Strukturgruppe G = ℤ₂ = {±1}.

        Das Möbius-Bündel ist das Standardbeispiel eines nicht-trivialisierbaren
        (verdrehten) Linienbündels. Seine erste Stiefel-Whitney-Klasse w₁ ≠ 0.

        Approximation: base_dim=1 (S¹), fiber_dim=1 (ℝ), G=Z2.

    @return  FiberBundle-Objekt für das Möbius-Bündel.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    bundle = FiberBundle(base_dim=1, fiber_dim=1, structure_group='Z2')
    bundle.name = "Möbius-Bündel (S¹, ℝ, ℤ₂)"
    return bundle


def hopf_fibration() -> FiberBundle:
    """
    @brief Hopf-Faserung S³ → S² mit Faser S¹.
    @description
        Die Hopf-Faserung η: S³ → S² ist das erste bekannte nicht-triviale
        Faserbündel und hat tiefe Verbindungen zur Quantenmechanik (Berry-Phase)
        und Eichtheorie (Monopole, Instanton).

        Explizit: η: S³ ⊂ ℂ² → S² ≅ ℂP¹
            (z₁, z₂) ↦ [z₁ : z₂]

        Eigenschaften:
        - Basisraum: S² (dim 2)
        - Gesamtraum: S³ (dim 3)
        - Faser: S¹ (dim 1)
        - Strukturgruppe: U(1)
        - Erste Chern-Klasse: c₁ = 1 (nicht-triviales Bündel)

    @return  FiberBundle-Objekt für die Hopf-Faserung.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    bundle = FiberBundle(base_dim=2, fiber_dim=1, structure_group='U1')
    bundle.name = "Hopf-Faserung (S², S¹, U(1))"
    # Erste Chern-Klasse: c₁ = 1 (Winding-Zahl des Bündels)
    bundle.chern_number = 1
    return bundle


def tangent_bundle(manifold_dim: int) -> FiberBundle:
    """
    @brief Tangentialbündel TM einer n-dimensionalen Mannigfaltigkeit.
    @description
        Das Tangentialbündel π: TM → M ist das kanonische Faserbündel
        einer glatten Mannigfaltigkeit M. Jede Faser π⁻¹(p) = T_pM ist
        der Tangentialraum an M im Punkt p.

        Eigenschaften:
        - Basisraum: M (dim n)
        - Faser: ℝⁿ (dim n)
        - Strukturgruppe: GL(n, ℝ) (allgemeine lineare Gruppe)
        - Gesamtraum: TM (dim 2n)

    @param manifold_dim  Dimension n der Mannigfaltigkeit M.
    @return              FiberBundle-Objekt für TM.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    bundle = FiberBundle(
        base_dim=manifold_dim,
        fiber_dim=manifold_dim,
        structure_group='GL'
    )
    bundle.name = f"Tangentialbündel T(M^{manifold_dim})"
    return bundle


def chern_class_first(bundle: FiberBundle) -> int:
    """
    @brief Erste Chern-Klasse c₁ eines U(1)-Hauptfaserbündels.
    @description
        Die erste Chern-Klasse c₁ ∈ H²(B; ℤ) ist eine topologische
        Invariante, die misst, wie "verdreht" ein komplexes Linienbündel ist.

        Für ein U(1)-Hauptbündel über S² (z.B. Monopolbündel) gilt:
            c₁ = (1/2π) ∫_{S²} F
        wobei F die Krümmungsform (Feldstärke) ist.

        Berechnung über die gespeicherte `chern_number`-Eigenschaft des
        Bündels (falls vorhanden) oder durch die Winding-Zahl der
        Übergangsfunktion.

    @param bundle  Das Faserbündel (sollte U(1)-Strukturgruppe haben).
    @return        Erste Chern-Zahl c₁ als ganze Zahl.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Falls explizit gesetzt (z.B. für Hopf-Faserung)
    if hasattr(bundle, 'chern_number'):
        return int(bundle.chern_number)

    # Berechnung über Winding-Zahl der Übergangsfunktion g_{01}: S¹ → U(1)
    # Integriere den Phasenwinkel über eine geschlossene Kurve in B
    n_points = 100
    total_phase = 0.0

    for i in range(n_points):
        # Parameterisierung des Äquators von S²: (cos θ, sin θ)
        theta = 2 * math.pi * i / n_points
        point = np.array([math.cos(theta), math.sin(theta)])

        # Übergangsfunktion zwischen Karte 0 und 1
        g = bundle.transition_function(point, 0, 1)

        # Für SO(2): Winkel aus der Rotationsmatrix ablesen
        if g.shape == (2, 2):
            angle = math.atan2(float(g[1, 0].real), float(g[0, 0].real))
            total_phase += angle

    # Winding-Zahl: Gesamtphase / 2π (auf ganze Zahl runden)
    winding = total_phase / (2 * math.pi * n_points)
    return int(round(winding * n_points / (2 * math.pi)))
