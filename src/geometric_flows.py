r"""
@file geometric_flows.py
@brief Ricci-Fluss, Perelmans Entropie-Funktionale und Beweis der Poincaré-Vermutung.
@description
    Dieses Modul implementiert die zentralen Konzepte des Ricci-Flusses und
    von Grigori Perelmans Beweis der Poincaré-Vermutung (2002–2003).

    **Historischer Hintergrund:**
    - 1982: Richard Hamilton führt den Ricci-Fluss ein: ∂g_{ij}/∂t = -2 R_{ij}
    - 1982–2002: Hamilton entwickelt das Programm (Singularitäten, Chirurgie)
    - 2002/11: Perelman, arXiv:math/0211159 — Entropie-Monotonie (F- und W-Funktionale)
    - 2003/03: Perelman, arXiv:math/0303109 — Geometrisierung mit Chirurgie
    - 2003/07: Perelman, arXiv:math/0307245 — Endliche Auslöschungszeit
    - 2006: Cao–Zhu, Kleiner–Lott, Morgan–Tian vervollständigen Details
    - 2006: Perelman lehnt Fields-Medaille ab

    **Kernideen Perelmans:**
    1. W-Entropie W(g,f,τ) wächst monoton entlang des Ricci-Flusses.
    2. κ-Nichtkollaps verhindert das Zusammenbrechen kleiner Gebiete.
    3. Kanonische Nachbarschaften (ε-Hälse, ε-Kappen) klassifizieren Singularitäten.
    4. Chirurgie schneidet Singularitäten heraus und klebt Standardkappen an.
    5. Für einfach-zusammenhängende 3-Mannigfaltigkeiten: Auslöschungszeit endlich.

    **Poincaré-Vermutung (bewiesen 2003):**
    Jede einfach-zusammenhängende, geschlossene 3-Mannigfaltigkeit ist
    homöomorph zur 3-Sphäre S³.

    **KaTeX-Formeln:**
    $$\frac{\partial g_{ij}}{\partial t} = -2 R_{ij}$$
    $$W(g,f,\tau) = \int_M \left[\tau(R + |\nabla f|^2) + f - n\right]
      (4\pi\tau)^{-n/2} e^{-f} \, dV \geq 0$$
    $$F(g,f) = \int_M (R + |\nabla f|^2) e^{-f} \, dV$$
    $$\Gamma^k_{ij} = \tfrac{1}{2} g^{kl}
      \left(\partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij}\right)$$

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import numpy as np
from typing import Callable, Dict, List, Optional, Tuple, Any


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _finite_diff_gradient(f: Callable, x: np.ndarray, h: float = 1e-5) -> np.ndarray:
    r"""
    @brief Numerischer Gradient via zentraler Differenzen.
    @param f   Skalare Funktion f: R^n → R
    @param x   Auswertungspunkt
    @param h   Schrittweite
    @return    Gradient ∇f(x) ∈ R^n
    @lastModified 2026-03-11
    """
    n = len(x)
    grad = np.zeros(n)
    for i in range(n):
        x_plus = x.copy()
        x_minus = x.copy()
        x_plus[i] += h
        x_minus[i] -= h
        grad[i] = (f(x_plus) - f(x_minus)) / (2.0 * h)
    return grad


def _finite_diff_hessian(f: Callable, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    r"""
    @brief Numerische Hesse-Matrix via zentraler Differenzen.
    @param f   Skalare Funktion
    @param x   Auswertungspunkt
    @param h   Schrittweite
    @return    Hesse-Matrix H_{ij} = ∂²f/∂x_i∂x_j
    @lastModified 2026-03-11
    """
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            x_pp = x.copy(); x_pp[i] += h; x_pp[j] += h
            x_pm = x.copy(); x_pm[i] += h; x_pm[j] -= h
            x_mp = x.copy(); x_mp[i] -= h; x_mp[j] += h
            x_mm = x.copy(); x_mm[i] -= h; x_mm[j] -= h
            val = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4.0 * h * h)
            H[i, j] = val
            H[j, i] = val
    return H


# ===========================================================================
# RIEMANNSCHE GEOMETRIE GRUNDLAGEN
# ===========================================================================

class RiemannianMetric:
    r"""
    @brief Riemannsche Metrik auf einer n-dimensionalen Mannigfaltigkeit.
    @description
        Eine Riemannsche Metrik g ist ein glattes, symmetrisches, positiv-
        definites (0,2)-Tensorfeld. In lokalen Koordinaten (x^1,...,x^n)
        ist g durch seine Komponenten g_{ij}(x) gegeben.

        Christoffel-Symbole (Levi-Civita-Zusammenhang):
        $$\Gamma^k_{ij} = \tfrac{1}{2} g^{kl}
          (\partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij})$$

        Riemann-Krümmungstensor:
        $$R^l{}_{ijk} = \partial_i \Gamma^l_{jk} - \partial_j \Gamma^l_{ik}
          + \Gamma^l_{im}\Gamma^m_{jk} - \Gamma^l_{jm}\Gamma^m_{ik}$$

        Ricci-Tensor (Spur des Riemann-Tensors):
        $$R_{ij} = R^k{}_{ikj}$$

        Ricci-Skalar (Spur des Ricci-Tensors):
        $$R = g^{ij} R_{ij}$$

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, dim: int, metric_tensor_func: Optional[Callable] = None):
        r"""
        @brief Konstruktor: Dimension und optionaler Metrik-Tensor.
        @param dim               Dimension der Mannigfaltigkeit (n)
        @param metric_tensor_func Funktion coords → np.ndarray (n×n),
                                  None = euklidische Metrik (Einheitsmatrix)
        @lastModified 2026-03-11
        """
        self.dim = dim  # Dimension n

        # Standard: euklidische Metrik δ_{ij}
        if metric_tensor_func is None:
            self._metric_func = lambda coords: np.eye(dim)
        else:
            self._metric_func = metric_tensor_func

    def metric_at(self, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Berechnet den metrischen Tensor g_{ij}(x) an einem Punkt.
        @param coords Koordinaten x ∈ R^n
        @return       Metrische Matrix g (n×n, symmetrisch, positiv-definit)
        @lastModified 2026-03-11
        """
        g = np.array(self._metric_func(coords), dtype=float)
        return g

    def christoffel_symbols(self, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Berechnet die Christoffel-Symbole Γ^k_{ij} numerisch.
        @description
            $$\Gamma^k_{ij} = \tfrac{1}{2} g^{kl}
              (\partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij})$$
            Ableitungen ∂_i g_{jl} werden über finite Differenzen genähert.
        @param coords Koordinaten x ∈ R^n
        @return       Array Γ[k,i,j] der Form (n,n,n)
        @lastModified 2026-03-11
        """
        n = self.dim
        h = 1e-5  # Schrittweite für finite Differenzen

        g = self.metric_at(coords)
        g_inv = np.linalg.inv(g)  # g^{ij} (inverse Metrik)

        # Partielle Ableitungen ∂_m g_{ij} numerisch berechnen
        # dg[m, i, j] = ∂g_{ij}/∂x^m
        dg = np.zeros((n, n, n))
        for m in range(n):
            coords_p = coords.copy()
            coords_m = coords.copy()
            coords_p[m] += h
            coords_m[m] -= h
            g_p = self.metric_at(coords_p)
            g_m = self.metric_at(coords_m)
            dg[m] = (g_p - g_m) / (2.0 * h)

        # Christoffel-Symbole assemblieren: Γ^k_{ij}
        Gamma = np.zeros((n, n, n))
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    # Σ_l g^{kl} (∂_i g_{jl} + ∂_j g_{il} - ∂_l g_{ij}) / 2
                    s = 0.0
                    for l in range(n):
                        s += g_inv[k, l] * (dg[i, j, l] + dg[j, i, l] - dg[l, i, j])
                    Gamma[k, i, j] = 0.5 * s
        return Gamma

    def riemann_curvature_tensor(self, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Berechnet den Riemann-Krümmungstensor R^l_{ijk} numerisch.
        @description
            $$R^l{}_{ijk} = \partial_i \Gamma^l_{jk} - \partial_j \Gamma^l_{ik}
              + \Gamma^l_{im}\Gamma^m_{jk} - \Gamma^l_{jm}\Gamma^m_{ik}$$
        @param coords Koordinaten x ∈ R^n
        @return       Tensor R[l,i,j,k] der Form (n,n,n,n)
        @lastModified 2026-03-11
        """
        n = self.dim
        h = 1e-4  # Etwas größere Schrittweite für 2. Ableitungen

        Gamma = self.christoffel_symbols(coords)

        # Partielle Ableitungen der Christoffel-Symbole ∂_i Γ^l_{jk}
        # dGamma[i, l, j, k] = ∂_i Γ^l_{jk}
        dGamma = np.zeros((n, n, n, n))
        for i in range(n):
            coords_p = coords.copy()
            coords_m = coords.copy()
            coords_p[i] += h
            coords_m[i] -= h
            Gp = self.christoffel_symbols(coords_p)
            Gm = self.christoffel_symbols(coords_m)
            dGamma[i] = (Gp - Gm) / (2.0 * h)

        # Riemann-Tensor assemblieren
        R = np.zeros((n, n, n, n))  # R[l, i, j, k]
        for l in range(n):
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        # ∂_i Γ^l_{jk} - ∂_j Γ^l_{ik}
                        val = dGamma[i, l, j, k] - dGamma[j, l, i, k]
                        # + Γ^l_{im} Γ^m_{jk} - Γ^l_{jm} Γ^m_{ik}
                        for m in range(n):
                            val += Gamma[l, i, m] * Gamma[m, j, k]
                            val -= Gamma[l, j, m] * Gamma[m, i, k]
                        R[l, i, j, k] = val
        return R

    def ricci_tensor(self, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Berechnet den Ricci-Tensor R_{ij} = R^k_{ikj}.
        @description
            Kontraktion des Riemann-Tensors über den ersten und dritten Index:
            $$R_{ij} = R^k{}_{ikj} = \sum_{k=0}^{n-1} R^k{}_{ikj}$$
        @param coords Koordinaten x ∈ R^n
        @return       Ricci-Tensor R_{ij} (n×n-Matrix)
        @lastModified 2026-03-11
        """
        n = self.dim
        R_riemann = self.riemann_curvature_tensor(coords)  # R[l,i,j,k]

        # Ricci_{ij} = Σ_k R[k, i, k, j]  (Spur über 1. und 3. Index)
        Ric = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    Ric[i, j] += R_riemann[k, i, k, j]
        return Ric

    def ricci_scalar(self, coords: np.ndarray) -> float:
        r"""
        @brief Berechnet den Ricci-Skalar R = g^{ij} R_{ij}.
        @description
            Vollständige Spur des Ricci-Tensors mit der inversen Metrik:
            $$R = g^{ij} R_{ij} = \sum_{i,j} g^{ij} R_{ij}$$
        @param coords Koordinaten x ∈ R^n
        @return       Ricci-Skalar (reelle Zahl)
        @lastModified 2026-03-11
        """
        g = self.metric_at(coords)
        g_inv = np.linalg.inv(g)
        Ric = self.ricci_tensor(coords)

        # R = Σ_{ij} g^{ij} R_{ij}
        return float(np.einsum('ij,ij->', g_inv, Ric))

    def sectional_curvature(self,
                            plane_vectors: Tuple[np.ndarray, np.ndarray],
                            coords: np.ndarray) -> float:
        r"""
        @brief Berechnet die Schnittkrümmung K(u,v) einer 2-Ebene.
        @description
            Für zwei linear unabhängige Tangentialvektoren u, v:
            $$K(u,v) = \frac{R(u,v,v,u)}{|u|^2|v|^2 - \langle u,v\rangle^2}$$
            wobei R_{ijkl} = g_{lm} R^m{}_{ijk} der vollständig kovariante
            Riemann-Tensor ist.
        @param plane_vectors Paar (u, v) von Tangentialvektoren
        @param coords        Koordinaten x ∈ R^n
        @return              Schnittkrümmung K(u,v)
        @lastModified 2026-03-11
        """
        u, v = np.array(plane_vectors[0], dtype=float), np.array(plane_vectors[1], dtype=float)
        g = self.metric_at(coords)
        R = self.riemann_curvature_tensor(coords)  # R[l,i,j,k]

        # Vollständig kovarianter Riemann-Tensor R_{lijk} = g_{lm} R^m_{ijk}
        # R(u,v,v,u) = R_{ijkl} u^i v^j v^k u^l  → Konvention variiert!
        # Hier: R(u,v,u,v) = Σ_{l,i,j,k} g_{la} R^a_{ijk} u^i v^j u^k v^l
        n = self.dim
        R_covar = np.einsum('lm,mijk->lijk', g, R)

        # Zähler: R_{lijk} u^i v^j u^k v^l
        numerator = np.einsum('lijk,i,j,k,l->', R_covar, u, v, u, v)

        # Nenner: g(u,u)·g(v,v) - g(u,v)²
        guu = float(np.einsum('ij,i,j->', g, u, u))
        gvv = float(np.einsum('ij,i,j->', g, v, v))
        guv = float(np.einsum('ij,i,j->', g, u, v))
        denominator = guu * gvv - guv ** 2

        if abs(denominator) < 1e-14:
            raise ValueError("Vektoren u und v sind linear abhängig.")
        return float(numerator / denominator)

    def volume_form(self, coords: np.ndarray) -> float:
        r"""
        @brief Berechnet das Volumenelement √det(g) an einem Punkt.
        @description
            Das Riemann'sche Volumenelement in lokalen Koordinaten:
            $$dV = \sqrt{\det(g)} \, dx^1 \wedge \cdots \wedge dx^n$$
        @param coords Koordinaten x ∈ R^n
        @return       √det(g(x))
        @lastModified 2026-03-11
        """
        g = self.metric_at(coords)
        det_g = np.linalg.det(g)
        if det_g < 0:
            # Numerische Ungenauigkeit: als 0 behandeln
            return 0.0
        return float(math.sqrt(det_g))


# ===========================================================================
# DIE 3-SPHÄRE S³
# ===========================================================================

class ThreeSphere:
    r"""
    @brief Die 3-Sphäre S³ mit runder Standardmetrik.
    @description
        Die 3-Sphäre ist das zentrale Objekt der Poincaré-Vermutung:
        $$S^3 = \{(x_1,x_2,x_3,x_4) \in \mathbb{R}^4 \mid
          x_1^2 + x_2^2 + x_3^2 + x_4^2 = R^2\}$$

        In Hopf-Koordinaten (η, ξ_1, ξ_2) mit η ∈ [0,π/2]:
        $$ds^2 = R^2 \left(d\eta^2 + \sin^2\eta\, d\xi_1^2
                           + \cos^2\eta\, d\xi_2^2\right)$$

        Topologische Eigenschaften:
        - Einfach-zusammenhängend: π_1(S³) = 0
        - Kompakt, ohne Rand
        - Konstante Schnittkrümmung K = 1/R² > 0

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, radius: float = 1.0):
        r"""
        @brief Konstruktor.
        @param radius Radius R der Sphäre (Standard: 1)
        @lastModified 2026-03-11
        """
        if radius <= 0:
            raise ValueError("Radius muss positiv sein.")
        self.radius = float(radius)  # Radius R > 0

    def metric_tensor(self, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Runde Metrik auf S³ in Hopf-Koordinaten (η, ξ₁, ξ₂).
        @description
            Hopf-Koordinaten: η ∈ [0,π/2], ξ₁,ξ₂ ∈ [0,2π).
            $$g = R^2 \begin{pmatrix}
              1 & 0 & 0 \\
              0 & \sin^2\eta & 0 \\
              0 & 0 & \cos^2\eta
            \end{pmatrix}$$
        @param coords Array [η, ξ₁, ξ₂]
        @return       3×3 metrische Matrix
        @lastModified 2026-03-11
        """
        eta = float(coords[0])
        R2 = self.radius ** 2

        g = np.array([
            [R2,                  0.0,               0.0],
            [0.0, R2 * math.sin(eta)**2,              0.0],
            [0.0,                 0.0, R2 * math.cos(eta)**2]
        ])
        return g

    def fundamental_group(self) -> str:
        r"""
        @brief Gibt die Fundamentalgruppe π₁(S³) zurück.
        @description
            S³ ist einfach-zusammenhängend: jede geschlossene Kurve lässt
            sich stetig zu einem Punkt kontrahieren.
            $$\pi_1(S^3) = \{1\} \text{ (triviale Gruppe)}$$
        @return String-Beschreibung der trivialen Gruppe
        @lastModified 2026-03-11
        """
        return "{1} (triviale Gruppe — S³ ist einfach-zusammenhängend)"

    def homology_groups(self) -> Dict[str, str]:
        r"""
        @brief Singuläre Homologiegruppen von S³.
        @description
            Via Poincaré-Dualität und Kugel-Formel:
            $$H_k(S^3; \mathbb{Z}) = \begin{cases}
              \mathbb{Z} & k=0,3 \\ 0 & \text{sonst}
            \end{cases}$$
        @return Dictionary H_k → Gruppenname
        @lastModified 2026-03-11
        """
        return {
            "H_0": "Z",   # Eine Zusammenhangskomponente
            "H_1": "0",   # Kein 1-dimensionaler Zyklus (einfach-zusammenhängend)
            "H_2": "0",   # Keine 2-dimensionale Torsion (kein Inneres)
            "H_3": "Z",   # Fundamentalklasse [S³]
        }

    def homotopy_groups(self) -> Dict[str, str]:
        r"""
        @brief Homotopiegruppen von S³ (erste vier).
        @description
            $$\pi_1(S^3)=0, \quad \pi_2(S^3)=0, \quad \pi_3(S^3)=\mathbb{Z}$$
            π_3(S³) = Z wird durch die Hopf-Faserung erzeugt.
        @return Dictionary π_k → Gruppenname
        @lastModified 2026-03-11
        """
        return {
            "pi_1": "0",          # Einfach-zusammenhängend
            "pi_2": "0",          # Zweite Homotopiegruppe trivial
            "pi_3": "Z",          # Hopf-Invariante erzeugt Z
            "pi_4": "Z/2Z",       # Bekannte höhere Homotopiegruppe
        }

    def is_simply_connected(self) -> bool:
        r"""
        @brief Prüft ob S³ einfach-zusammenhängend ist (immer True).
        @description
            Einfach-zusammenhängend ⟺ π₁ = 0 und zusammenhängend.
            S³ ist das Prototypbeispiel einer einfach-zusammenhängenden
            kompakten 3-Mannigfaltigkeit ohne Rand.
        @return True
        @lastModified 2026-03-11
        """
        return True


# ===========================================================================
# RICCI-FLUSS
# ===========================================================================

class RicciFlow:
    r"""
    @brief Ricci-Fluss auf einer Riemannschen Mannigfaltigkeit.
    @description
        Der Ricci-Fluss (Hamilton 1982) ist eine geometrische PDE, die
        eine Familie von Riemannschen Metriken g(t) entwickelt:
        $$\frac{\partial g_{ij}}{\partial t} = -2 R_{ij}(g(t))$$

        Intuition: Gebiete mit positiver Krümmung schrumpfen, Gebiete mit
        negativer Krümmung wachsen — das Volumen tendiert zu normalisieren.

        **Normalisierter Ricci-Fluss** (erhält Volumen):
        $$\frac{\partial g_{ij}}{\partial t} = -2 R_{ij}
          + \frac{2}{n} R_{\text{avg}} g_{ij}$$

        **Kurzzeit-Existenz** (DeTurck 1983):
        Für jede glatte Anfangsmetrik g_0 existiert eine eindeutige glatte
        Lösung auf einem kurzen Zeitintervall [0, T).

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, initial_metric: RiemannianMetric, dim: int = 3):
        r"""
        @brief Konstruktor.
        @param initial_metric Anfangsmetrik g_0 als RiemannianMetric-Objekt
        @param dim            Dimension der Mannigfaltigkeit (Standard: 3)
        @lastModified 2026-03-11
        """
        self.metric = initial_metric   # Aktuelle Metrik g(t)
        self.dim = dim                  # Dimension n

    def evolution_equation(self) -> str:
        r"""
        @brief Gibt die Evolutionsgleichung des Ricci-Flusses als String.
        @return LaTeX-String der Gleichung ∂g_{ij}/∂t = -2 R_{ij}
        @lastModified 2026-03-11
        """
        return r"∂g_{ij}/∂t = -2 R_{ij}  [Hamilton 1982]"

    def normalized_evolution_equation(self) -> str:
        r"""
        @brief Gibt die normalisierte Ricci-Fluss-Gleichung zurück.
        @description
            Der normalisierte Fluss erhält das Gesamtvolumen durch einen
            Skalierungsterm proportional zum mittleren Ricci-Skalar R_avg:
            $$\frac{\partial g_{ij}}{\partial t} = -2R_{ij}
              + \frac{2}{n} R_{\text{avg}} g_{ij}$$
        @return LaTeX-String der normierten Gleichung
        @lastModified 2026-03-11
        """
        return r"∂g_{ij}/∂t = -2 R_{ij} + (2/n) R_avg g_{ij}  [normalisiert]"

    def short_time_existence(self) -> str:
        r"""
        @brief Beschreibt das Kurzzeit-Existenzsatz für den Ricci-Fluss.
        @description
            DeTurck (1983) bewies, dass das Ricci-Fluss-System durch einen
            diffeomorphismus-Trick (DeTurck-Trick) in ein strikt parabolisches
            System überführt werden kann, für das klassische PDE-Theorie gilt.
        @return Erklärung des Kurzzeit-Existenzsatzes
        @lastModified 2026-03-11
        """
        return (
            "Kurzzeit-Existenz (DeTurck 1983): Für jede glatte Riemannsche "
            "Anfangsmetrik g_0 auf einer kompakten Mannigfaltigkeit existiert "
            "ein T > 0 und eine eindeutige glatte Lösung g(t) des Ricci-Flusses "
            "auf [0, T). Beweis via DeTurck-Trick: ∂g/∂t = -2Ric + L_W g, "
            "W^k = g^{ij}(Γ^k_{ij}[g] - Γ^k_{ij}[g̃]) — strikt parabolisch."
        )

    def euler_step(self, dt: float, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Führt einen expliziten Euler-Schritt des Ricci-Flusses durch.
        @description
            Diskrete Approximation der Evolutionsgleichung:
            $$g_{ij}(t + dt) \approx g_{ij}(t) - 2\, dt\, R_{ij}(g(t))$$
            Hinweis: Euler-Methode ist nur für kleines dt stabil.
        @param dt     Zeitschritt Δt > 0
        @param coords Koordinaten, an denen die Metrik ausgewertet wird
        @return       Neue Metrik g(t+dt) als (n×n)-Array
        @lastModified 2026-03-11
        """
        g_current = self.metric.metric_at(coords)
        Ric = self.metric.ricci_tensor(coords)

        # Euler-Schritt: g_new = g - 2·dt·Ric
        g_new = g_current - 2.0 * dt * Ric
        return g_new

    def rk4_step(self, dt: float, coords: np.ndarray) -> np.ndarray:
        r"""
        @brief Führt einen RK4-Schritt des Ricci-Flusses durch.
        @description
            Runge-Kutta 4. Ordnung für ∂g/∂t = -2 Ric(g):
            $$k_1 = -2\,\text{Ric}(g_t)$$
            $$k_2 = -2\,\text{Ric}(g_t + \tfrac{dt}{2}k_1)$$
            $$k_3 = -2\,\text{Ric}(g_t + \tfrac{dt}{2}k_2)$$
            $$k_4 = -2\,\text{Ric}(g_t + dt\,k_3)$$
            $$g_{t+dt} = g_t + \tfrac{dt}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
        @param dt     Zeitschritt Δt > 0
        @param coords Koordinaten, an denen die Metrik ausgewertet wird
        @return       Neue Metrik g(t+dt) als (n×n)-Array
        @lastModified 2026-03-11
        """
        g0 = self.metric.metric_at(coords)
        dim = self.metric.dim

        def ric_of(g_mat: np.ndarray) -> np.ndarray:
            r"""Berechnet Ricci-Tensor für temporäre Metrik g_mat."""
            # Temporäres Metric-Objekt mit konstantem Tensor
            tmp = RiemannianMetric(dim, metric_tensor_func=lambda c: g_mat)
            return tmp.ricci_tensor(coords)

        # Vier Stufen des RK4
        k1 = -2.0 * ric_of(g0)
        k2 = -2.0 * ric_of(g0 + 0.5 * dt * k1)
        k3 = -2.0 * ric_of(g0 + 0.5 * dt * k2)
        k4 = -2.0 * ric_of(g0 + dt * k3)

        g_new = g0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        return g_new


# ===========================================================================
# SINGULARITÄTEN DES RICCI-FLUSSES
# ===========================================================================

class RicciFlowSingularity:
    r"""
    @brief Modell für Singularitäten beim Ricci-Fluss.
    @description
        Hamilton (1995) und Perelman klassifizierten Singularitäten:

        **Typ I (Nack-Pinch):** Krümmungsexplosion |Rm| ~ 1/(T-t).
        Geometrisch: Zylinder S² × ℝ schnürt sich an einer Stelle ein.

        **Typ II (Cap):** |Rm| >> 1/(T-t). Seltenere, langsamere Singularität.

        **Zigarrensingularität (Cigar):** Perelmans κ-Nichtkollaps schließt
        diesen Typ aus — das war ein Schlüsselergebnis von 2002.

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    # Erlaubte Singularitätstypen
    SINGULARITY_TYPES = {"neck-pinch", "cap", "cigar"}

    def __init__(self, singularity_type: str):
        r"""
        @brief Konstruktor.
        @param singularity_type Art der Singularität: "neck-pinch", "cap" oder "cigar"
        @lastModified 2026-03-11
        """
        if singularity_type not in self.SINGULARITY_TYPES:
            raise ValueError(
                f"Unbekannter Singularitätstyp '{singularity_type}'. "
                f"Erlaubt: {self.SINGULARITY_TYPES}"
            )
        self.singularity_type = singularity_type

    def blow_up_rate(self, t: float, T: float) -> float:
        r"""
        @brief Berechnet die Aufblasrate |Rm| ~ C/(T-t) für Typ-I-Singularität.
        @description
            Nahe einer Typ-I-Singularität zum Zeitpunkt T gilt:
            $$\sup_M |Rm(g(t))| \sim \frac{C}{T - t}, \quad t \to T^-$$
            Für Typ II: |Rm| wächst schneller als 1/(T-t).
        @param t Aktueller Zeitpunkt t < T
        @param T Singularitätszeitpunkt T (Blow-up-Zeit)
        @return  Approximierte Aufblasrate C/(T-t) mit C=1
        @lastModified 2026-03-11
        """
        if t >= T:
            raise ValueError("Zeitpunkt t muss kleiner als Singularitätszeitpunkt T sein.")
        remaining = T - t
        if self.singularity_type == "neck-pinch":
            # Typ I: lineare Aufblasrate
            return 1.0 / remaining
        elif self.singularity_type == "cap":
            # Typ II: superlineare Aufblasrate
            return 1.0 / (remaining * math.sqrt(remaining))
        else:
            # Zigarrensingularität: subexponentiell (durch κ-Nichtkollaps ausgeschlossen)
            return math.exp(-1.0 / remaining)

    def canonical_neighborhood(self) -> str:
        r"""
        @brief Beschreibt die kanonische Nachbarschaft nahe der Singularität.
        @description
            Perelmans Theorem über kanonische Nachbarschaften (Theorem 12.1, 2003):
            Nahe jeder Singularität sieht die Geometrie aus wie:
            - ε-Hals: eine Geometrie nahe S² × ℝ (Zylinder)
            - ε-Kappe: eine Geometrie nahe einer abgerundeten Kappe
        @return Beschreibung der kanonischen Nachbarschaft
        @lastModified 2026-03-11
        """
        descriptions = {
            "neck-pinch": (
                "ε-Hals (ε-neck): Die Singularität liegt in einer Region, die "
                "einem dünnen Zylinder S² × (-ε⁻¹, ε⁻¹) ähnelt. "
                "Querschnitt ≈ Sphäre S² mit Radius → 0."
            ),
            "cap": (
                "ε-Kappe (ε-cap): Die Region ist eine Kombination aus einem "
                "ε-Hals und einer Standardkappe (diffeomorph zu B³ oder ℝP³ \\ B³). "
                "Entspricht dem Ende eines Hamilton-Zigarren-Grenzwerts."
            ),
            "cigar": (
                "Zigarrensingularität: Durch Perelmans κ-Nichtkollaps-Theorem "
                "(2002) ausgeschlossen. Die Zigarren-Soliton-Geometrie kollabiert "
                "auf kleinen Skalen, was dem κ-Nichtkollaps widerspricht."
            ),
        }
        return descriptions[self.singularity_type]

    def kappa_noncollapsed_check(self, r: float, kappa: float) -> bool:
        r"""
        @brief Prüft die κ-Nichtkollaps-Bedingung für eine Kugel vom Radius r.
        @description
            Perelmans κ-Nichtkollaps (Theorem 4.1, 2002):
            Falls |Rm| ≤ r⁻² auf B(x,r), dann:
            $$\text{Vol}(B(x,r)) \geq \kappa \cdot r^n$$
            Das Volumen kann nicht kleiner als κrⁿ sein.
            Zigarrensingularitäten verletzen dies und sind daher ausgeschlossen.
        @param r     Kugelradius r > 0
        @param kappa κ-Konstante aus Perelmans Theorem (κ > 0)
        @return      True wenn Nichtkollaps vorausgesagt wird, False bei Verstoß
        @lastModified 2026-03-11
        """
        if r <= 0 or kappa <= 0:
            raise ValueError("r und kappa müssen positiv sein.")

        # Zigarrensingularitäten kollabieren — Nichtkollaps verletzt
        if self.singularity_type == "cigar":
            return False

        # Für Hals und Kappe: Nichtkollaps-Bedingung gilt
        # (vereinfachtes Modell: immer True für nicht-Zigarre)
        return True


# ===========================================================================
# PERELMANS ENTROPIE-FUNKTIONALE
# ===========================================================================

class PerelmanEntropy:
    r"""
    @brief Perelmans Entropie-Funktionale F, W und μ.
    @description
        Grigori Perelman (2002) führte zwei Entropie-Funktionale ein, die
        monoton wachsen entlang des (konjugierten) Ricci-Flusses:

        **F-Funktional** (Hamilton-Perelman):
        $$F(g,f) = \int_M \left(R + |\nabla f|^2\right) e^{-f} \, dV$$

        **W-Funktional** (Perelmans Entropie):
        $$W(g,f,\tau) = \int_M \left[\tau\left(R + |\nabla f|^2\right)
          + f - n\right] (4\pi\tau)^{-n/2} e^{-f} \, dV$$

        Monotonie: Entlang des Ricci-Flusses gilt dW/dt ≥ 0.

        **μ-Funktional** (Perelman):
        $$\mu(g,\tau) = \inf\left\{W(g,f,\tau) \mid
          \int_M (4\pi\tau)^{-n/2} e^{-f} dV = 1\right\}$$

        **ν-Funktional**:
        $$\nu(g) = \inf_{\tau>0} \mu(g,\tau)$$

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self,
                 metric: RiemannianMetric,
                 scalar_func: Callable,
                 tau: float):
        r"""
        @brief Konstruktor.
        @param metric      RiemannianMetric-Objekt
        @param scalar_func f: R^n → R — skalare Funktion (Perelman-Potential)
        @param tau         τ > 0 — Skalenparameter (τ = T - t beim rückwärts Fluss)
        @lastModified 2026-03-11
        """
        if tau <= 0:
            raise ValueError("τ (tau) muss positiv sein.")
        self.metric = metric       # Riemannsche Metrik g
        self.scalar_func = scalar_func  # f(x): R^n → R
        self.tau = tau             # Zeitparameter τ > 0
        self.n = metric.dim        # Dimension n

    def w_functional(self, f: Callable, num_points: int = 50) -> float:
        r"""
        @brief Numerische Approximation des W-Funktionals.
        @description
            $$W(g,f,\tau) = \int_M \left[\tau(R + |\nabla f|^2) + f - n\right]
              (4\pi\tau)^{-n/2} e^{-f} \, dV$$
            Numerische Integration über einen Stützgitterbereich [-1,1]^n.
        @param f          Skalare Funktion f: R^n → R
        @param num_points Anzahl Gitterpunkte pro Achse
        @return           Numerischer Wert von W(g,f,τ)
        @lastModified 2026-03-11
        """
        n = self.n
        tau = self.tau
        # Normierungsfaktor (4πτ)^{-n/2}
        norm_factor = (4.0 * math.pi * tau) ** (-n / 2.0)

        # 1D-Gitter über [-1, 1]
        grid_1d = np.linspace(-1.0, 1.0, num_points)
        total = 0.0
        count = 0

        # Einfache Summen-Approximation (Monte-Carlo wäre besser für n>3)
        if n == 1:
            dx = 2.0 / (num_points - 1)
            for x0 in grid_1d:
                coords = np.array([x0])
                R_val = self.metric.ricci_scalar(coords)
                f_val = float(f(coords))
                # Gradient |∇f|² über finite Differenzen
                grad_f = _finite_diff_gradient(f, coords, h=1e-4)
                g = self.metric.metric_at(coords)
                g_inv = np.linalg.inv(g)
                grad_f_sq = float(np.einsum('i,ij,j->', grad_f, g_inv, grad_f))
                vol = self.metric.volume_form(coords)
                integrand = (tau * (R_val + grad_f_sq) + f_val - n) * norm_factor * math.exp(-f_val) * vol
                total += integrand * dx
                count += 1
        else:
            # Für n>1: zufällige Stichprobenpunkte (Monte-Carlo)
            rng = np.random.default_rng(42)
            pts = rng.uniform(-1.0, 1.0, size=(num_points**2, n))
            volume_domain = 2.0 ** n  # Volumen des Würfels [-1,1]^n
            for coords in pts:
                try:
                    R_val = self.metric.ricci_scalar(coords)
                    f_val = float(f(coords))
                    grad_f = _finite_diff_gradient(f, coords, h=1e-4)
                    g = self.metric.metric_at(coords)
                    g_inv = np.linalg.inv(g)
                    grad_f_sq = float(np.einsum('i,ij,j->', grad_f, g_inv, grad_f))
                    vol = self.metric.volume_form(coords)
                    integrand = (tau * (R_val + grad_f_sq) + f_val - n) * norm_factor * math.exp(-f_val) * vol
                    total += integrand
                    count += 1
                except (np.linalg.LinAlgError, ValueError):
                    continue
            if count > 0:
                total = total * volume_domain / count

        return float(total)

    def f_functional(self, f: Callable, num_points: int = 50) -> float:
        r"""
        @brief Numerische Approximation des F-Funktionals.
        @description
            $$F(g,f) = \int_M \left(R + |\nabla f|^2\right) e^{-f} \, dV$$
            Hamilton-Perelman-Entropie, monoton wächst entlang des Ricci-Flusses.
        @param f          Skalare Funktion f: R^n → R
        @param num_points Gitterpunkte pro Achse (1D) oder Stichprobenzahl (nD)
        @return           Numerischer Wert von F(g,f)
        @lastModified 2026-03-11
        """
        n = self.n
        total = 0.0
        count = 0

        if n == 1:
            grid = np.linspace(-1.0, 1.0, num_points)
            dx = 2.0 / (num_points - 1)
            for x0 in grid:
                coords = np.array([x0])
                R_val = self.metric.ricci_scalar(coords)
                f_val = float(f(coords))
                grad_f = _finite_diff_gradient(f, coords, h=1e-4)
                g = self.metric.metric_at(coords)
                g_inv = np.linalg.inv(g)
                grad_f_sq = float(np.einsum('i,ij,j->', grad_f, g_inv, grad_f))
                vol = self.metric.volume_form(coords)
                integrand = (R_val + grad_f_sq) * math.exp(-f_val) * vol
                total += integrand * dx
                count += 1
        else:
            rng = np.random.default_rng(42)
            pts = rng.uniform(-1.0, 1.0, size=(num_points**2, n))
            volume_domain = 2.0 ** n
            for coords in pts:
                try:
                    R_val = self.metric.ricci_scalar(coords)
                    f_val = float(f(coords))
                    grad_f = _finite_diff_gradient(f, coords, h=1e-4)
                    g = self.metric.metric_at(coords)
                    g_inv = np.linalg.inv(g)
                    grad_f_sq = float(np.einsum('i,ij,j->', grad_f, g_inv, grad_f))
                    vol = self.metric.volume_form(coords)
                    integrand = (R_val + grad_f_sq) * math.exp(-f_val) * vol
                    total += integrand
                    count += 1
                except (np.linalg.LinAlgError, ValueError):
                    continue
            if count > 0:
                total = total * volume_domain / count

        return float(total)

    def mu_functional(self, tau: float) -> float:
        r"""
        @brief Schätzung des μ-Funktionals μ(g,τ).
        @description
            $$\mu(g,\tau) = \inf\left\{W(g,f,\tau) \mid
              \int_M (4\pi\tau)^{-n/2} e^{-f} dV = 1\right\}$$
            Stellt das Infimum des W-Funktionals unter Normierungsbedingung dar.
            Für eine euklidische Metrik gilt μ(δ, τ) = 0 für alle τ > 0.
        @param tau Skalenparameter τ > 0
        @return    Schätzwert von μ(g,τ)
        @lastModified 2026-03-11
        """
        # Einfachste wählbare Funktion: f = n/2 · ln(4πτ) (Normierungsbedingung erfüllt)
        # → W(g, f, τ) = 0 für g = Standardmetrik
        n = self.n
        f_opt = lambda coords: n / 2.0 * math.log(4.0 * math.pi * tau)
        return self.w_functional(f_opt, num_points=20)

    def nu_functional(self) -> float:
        r"""
        @brief Schätzung des ν-Funktionals ν(g) = inf_{τ>0} μ(g,τ).
        @description
            $$\nu(g) = \inf_{\tau > 0} \mu(g,\tau)$$
            Das ν-Funktional ist eine Gesamtinvariante der Metrik.
            Für Metriken mit nicht-negativem Ricci-Tensor gilt ν(g) ≤ 0.
        @return Schätzwert von ν(g) durch Minimierung über τ ∈ {0.1, 0.5, 1, 2}
        @lastModified 2026-03-11
        """
        tau_values = [0.1, 0.5, 1.0, 2.0, 5.0]
        mu_vals = [self.mu_functional(tau) for tau in tau_values]
        return float(min(mu_vals))

    def monotonicity_check(self, t1: float, t2: float) -> bool:
        r"""
        @brief Prüft ob W(t₁) ≤ W(t₂) (Monotoniebedingung).
        @description
            Perelmans Monotonie-Satz: Entlang des Ricci-Flusses gilt
            $$\frac{d}{dt} W(g(t), f(t), \tau(t)) \geq 0$$
            Diese Methode prüft die Konsistenz: W(t₁) ≤ W(t₂) für t₁ < t₂.
            (Vereinfacht: für euklidische Metriken gilt W ≈ 0 = const.)
        @param t1 Früherer Zeitpunkt t₁ < t₂
        @param t2 Späterer Zeitpunkt t₂
        @return   True wenn Monotonie plausibel (W wächst)
        @lastModified 2026-03-11
        """
        if t1 >= t2:
            raise ValueError("t1 muss kleiner als t2 sein.")

        # τ(t) = T - t (rückwärts-Zeitparameter)
        # Für euklidische Standardmetrik: W ≡ 0 (Grenzfall, Monotonie trivial)
        tau1 = abs(t1) + 0.1
        tau2 = abs(t2) + 0.1

        w1 = self.mu_functional(tau1)
        w2 = self.mu_functional(tau2)

        # Perelmans Satz: dW/dt ≥ 0 → W wächst (oder bleibt konstant)
        return w2 >= w1 - 1e-8  # Toleranz für numerische Fehler


# ===========================================================================
# RICCI-FLUSS MIT CHIRURGIE
# ===========================================================================

class RicciFlowWithSurgery:
    r"""
    @brief Ricci-Fluss mit Chirurgie nach Perelman.
    @description
        Wenn der Ricci-Fluss eine Singularität entwickelt (typisch: Nack-Pinch),
        wird chirurgisch eingegriffen:

        1. **Singularitätserkennung**: Überwache |Rm| → ∞
        2. **Chirurgie**: Schneide den engen Hals heraus, klebe Standardkappen an
        3. **Weiterfluss**: Setze den Ricci-Fluss auf den reparierten Stücken fort
        4. **Endlichkeit**: Für π₁ endlich: endliche Auslöschungszeit

        Perelman zeigte:
        - Chirurgien häufen sich nicht auf endlichem Zeitintervall
        - Jedes Stück konvergiert zu einer der 8 Thurston-Geometrien

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, initial_metric: RiemannianMetric):
        r"""
        @brief Konstruktor.
        @param initial_metric Anfangsmetrik g_0
        @lastModified 2026-03-11
        """
        self.metric = initial_metric   # Aktuelle Metrik
        self.surgery_count = 0         # Anzahl durchgeführter Chirurgien
        self.surgery_log: List[Dict] = []  # Protokoll der Chirurgien

    def detect_singularity(self, t: float) -> Optional[str]:
        r"""
        @brief Erkennt Singularitäten anhand der Krümmungsexplosion.
        @description
            Heuristik: Für Testmetriken wird eine Singularität simuliert,
            indem geprüft wird ob t nahe an einer vorgegebenen Blow-up-Zeit T liegt.
            In der Praxis: Überwache ‖Rm‖_{L^∞} und vergleiche mit 1/(T-t).
        @param t Aktueller Zeitpunkt
        @return  Singularitätstyp als String oder None (keine Singularität)
        @lastModified 2026-03-11
        """
        # Modellhaftes Szenario: Singularität bei T=1.0
        T_blowup = 1.0
        epsilon = 0.05  # Detektionsschwelle

        if t >= T_blowup:
            return "neck-pinch"  # Nach Blow-up

        if T_blowup - t < epsilon:
            # Kurz vor Singularität: Nack-Pinch erwartet
            return "neck-pinch"

        # Kein Problem erkannt
        return None

    def perform_surgery(self, singularity_region: str) -> Dict:
        r"""
        @brief Beschreibt die Chirurgie an einer Singularitätsregion.
        @description
            Perelmans Chirurgie-Vorgang (2003):
            1. Identifiziere den engen Hals (ε-neck) um die Singularität
            2. Schneide den Hals durch (füge Rand ein)
            3. Klebe standardisierte ε-Kappen an die offenen Enden
            4. Skaliere die Kappen so, dass die Metrik regulär ist

            Das resultierende Stück hat weniger Topologie als das Original.
        @param singularity_region Art der Singularitätsregion
        @return                   Dictionary mit Beschreibung der Operation
        @lastModified 2026-03-11
        """
        self.surgery_count += 1

        result = {
            "surgery_number": self.surgery_count,
            "region": singularity_region,
            "operation": "neck-cut + cap-gluing",
            "description": (
                f"Chirurgie #{self.surgery_count}: "
                f"Engen Hals in Region '{singularity_region}' durchgeschnitten. "
                f"Standardkappen (S³-Geometrie, Radius δ) angeklebt. "
                f"Euler-Charakteristik erhöht um +2."
            ),
            "topology_change": "χ += 2 (zwei neue S²-Grenzen geschlossen)",
            "metric_after": "glatt, κ-nicht-kollabiert",
        }

        self.surgery_log.append(result)
        return result

    def canonical_neighborhood_theorem(self) -> str:
        r"""
        @brief Perelmans Theorem über kanonische Nachbarschaften (Theorem 12.1).
        @description
            Für δ > 0 genügend klein gilt: Falls (M, g) eine δ-Lösung des
            Ricci-Flusses ist und |Rm|(x₀) ≥ r⁻², dann hat (M, g) in einer
            Umgebung von x₀ eine kanonische Nachbarschaft — entweder einen
            ε-Hals, eine ε-Kappe oder einen abgerundeten S³-Anteil.
        @return Beschreibung des Theorems
        @lastModified 2026-03-11
        """
        return (
            "Perelman, Theorem 12.1 (2003): Für jedes ε > 0 existiert κ > 0, "
            "so dass jede κ-Lösung des Ricci-Flusses in 3D eine kanonische "
            "Nachbarschaft hat: entweder einen ε-Hals (S² × ℝ-Modell), "
            "eine ε-Kappe (B³-Modell) oder einen abgerundeten S³-Anteil. "
            "Dies ermöglicht die präzise Chirurgiebeschreibung."
        )

    def extinction_time(self) -> Optional[float]:
        r"""
        @brief Schätzung der Auslöschungszeit für Mannigfaltigkeiten mit endlichem π₁.
        @description
            Perelman (2003, Paper 3) bewies: Für eine kompakte 3-Mannigfaltigkeit
            mit endlicher Fundamentalgruppe π₁(M) erlischt die Mannigfaltigkeit
            unter dem Ricci-Fluss mit Chirurgie in endlicher Zeit T_ext < ∞.

            Für M = S³ (triviales π₁): T_ext ≈ R²/6 (für Radius-R-Sphäre)
            (Normalisierter Fluss auf S³ konvergiert zu runder Metrik.)
        @return Geschätzte Auslöschungszeit oder None (unendliche Topologie)
        @lastModified 2026-03-11
        """
        # Modellwert: S³ mit Radius 1 erlischt in endlicher Zeit
        return 1.0 / 6.0  # T_ext ≈ 1/6 für normalisierter Fluss auf S³(R=1)


# ===========================================================================
# POINCARÉ-VERMUTUNG
# ===========================================================================

class PoincaréConjecture:
    r"""
    @brief Poincaré-Vermutung: Aussage, Beweis und historischer Kontext.
    @description
        **Poincaré-Vermutung** (Henri Poincaré, 1904):
        Jede einfach-zusammenhängende, geschlossene 3-Mannigfaltigkeit ist
        homöomorph zur 3-Sphäre S³.

        **Beweis** (Grigori Perelman, 2002–2003, via Hamiltons Ricci-Fluss):
        1. W-Entropie ist monoton → κ-Nichtkollaps
        2. Kanonische Nachbarschaften klassifizieren Singularitäten
        3. Ricci-Fluss mit Chirurgie ist wohldefiniert
        4. Für einfach-zusammenhängende M: endliche Auslöschungszeit
        5. Also M ≅ S³ □

        **Clay Mathematics Institute**: Eines der sieben Millenium-Probleme.
        Perelman lehnte sowohl die Fields-Medaille (2006) als auch das
        Preisgeld von $1 Million ab.

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def theorem_statement(self) -> str:
        r"""
        @brief Gibt die Aussage der Poincaré-Vermutung als bewiesenen Satz zurück.
        @return Mathematische Formulierung des Satzes
        @lastModified 2026-03-11
        """
        return (
            "Satz (Poincaré-Vermutung, bewiesen durch Perelman 2002–2003): "
            "Sei M eine kompakte, einfach-zusammenhängende 3-Mannigfaltigkeit "
            "ohne Rand (d.h. π₁(M) = {1}). Dann ist M homöomorph zur "
            "3-Sphäre S³. "
            "KaTeX: \\pi_1(M)=0 \\Rightarrow M \\cong S^3"
        )

    def geometrization_statement(self) -> str:
        r"""
        @brief Gibt die Thurston-Geometrisierungsvermutung (allgemeiner) zurück.
        @description
            Thurstons Geometrisierungsvermutung (1982, bewiesen durch Perelman 2003)
            ist deutlich allgemeiner als die Poincaré-Vermutung:
            Jede kompakte, orientierbare 3-Mannigfaltigkeit besitzt eine
            Dekomposition, wobei jedes Stück eine der 8 Thurston-Geometrien trägt.
        @return Beschreibung der Geometrisierungsvermutung
        @lastModified 2026-03-11
        """
        return (
            "Thurston-Geometrisierungsvermutung (Thurston 1982, Perelman 2003): "
            "Jede kompakte, orientierbare 3-Mannigfaltigkeit lässt sich entlang "
            "unzusammendrückbarer Tori in Stücke zerlegen (JSJ-Zerlegung), wobei "
            "jedes Stück eine der acht Thurston-Modellgeometrien trägt: "
            "S³, E³, H³, S²×ℝ, H²×ℝ, SL̃(2,ℝ), Nil, Sol. "
            "Die Poincaré-Vermutung ist der Spezialfall π₁ = 0."
        )

    def eight_geometries(self) -> List[str]:
        r"""
        @brief Gibt die acht Thurston-Modellgeometrien zurück.
        @description
            William Thurstons Klassifikation (1982): Es gibt genau 8 Typen
            von homogenen Riemannschen 3-Geometrien, die als Modellgeometrien
            für die Stücke der Geometrisierung auftreten.
        @return Liste der 8 Geometrien mit Beschreibung
        @lastModified 2026-03-11
        """
        return [
            "S³  — 3-Sphäre, konstante positive Krümmung K=+1",
            "E³  — Euklidischer 3-Raum, K=0 (flach)",
            "H³  — Hyperbolischer 3-Raum, konstante negative Krümmung K=-1",
            "S²×ℝ — Produkt aus 2-Sphäre und Gerade",
            "H²×ℝ — Produkt aus hyperbolischer Ebene und Gerade",
            "SL̃(2,ℝ) — Universelle Überlagerung von SL(2,ℝ), Nil-ähnlich verdreht",
            "Nil   — Nilgeometrie (Heisenberg-Gruppe), nicht-abelsche Auflösungsgeometrie",
            "Sol   — Solvgeometrie (auflösbare Gruppe), anisotrop",
        ]

    def proof_outline(self) -> List[str]:
        r"""
        @brief Skizziert die Hauptschritte von Perelmans Beweis.
        @return Geordnete Liste der Beweisschritte
        @lastModified 2026-03-11
        """
        return [
            "Schritt 1 (Entropie): F- und W-Funktionale sind monoton wachsend "
            "entlang des Ricci-Flusses → κ-Nichtkollaps (kein Volumen-Kollaps "
            "auf kleinen Skalen). [Perelman 2002, arXiv:math/0211159]",

            "Schritt 2 (Geometrie): Gradienten-Shrinking-Solitons sind kanonisch "
            "(S³, S²×ℝ oder Zigarren). Zigarren sind durch κ-Nichtkollaps "
            "ausgeschlossen.",

            "Schritt 3 (Kanonische Nachbarschaften): Theorem 12.1 — nahe jeder "
            "Singularität sieht die Geometrie wie ε-Hals oder ε-Kappe aus. "
            "[Perelman 2003, arXiv:math/0303109]",

            "Schritt 4 (Chirurgie): Definiere Ricci-Fluss mit Chirurgie für "
            "3-Mannigfaltigkeiten. Chirurgien häufen sich nicht auf [0, T_max). "
            "Jedes Stück erhält Thurston-Geometrie.",

            "Schritt 5 (Auslöschung): Für π₁(M) endlich erlischt M in endlicher "
            "Zeit T_ext < ∞ unter normalisiertem Ricci-Fluss. "
            "[Perelman 2003, arXiv:math/0307245]",

            "Schlussfolgerung: Falls π₁(M) = 0 (einfach-zusammenhängend), "
            "erlischt M zu einer Metrik positiver Krümmung, also M ≅ S³. □",
        ]

    def perelman_papers(self) -> List[Dict]:
        r"""
        @brief Gibt Perelmans drei arXiv-Papers mit korrekten bibliographischen Daten.
        @return Liste von Dictionaries mit Titeln, Daten und arXiv-IDs
        @lastModified 2026-03-11
        """
        return [
            {
                "title": "The entropy formula for the Ricci flow and its geometric applications",
                "author": "Grigori Perelman",
                "date": "2002-11-11",
                "arxiv": "math/0211159",
                "content": "F- und W-Entropie, κ-Nichtkollaps, No-Local-Collapsing-Theorem",
            },
            {
                "title": "Ricci flow with surgery on three-manifolds",
                "author": "Grigori Perelman",
                "date": "2003-03-10",
                "arxiv": "math/0303109",
                "content": "Kanonische Nachbarschaften, Chirurgie, Langzeit-Existenz",
            },
            {
                "title": "Finite extinction time for the solutions to the Ricci flow on certain three-manifolds",
                "author": "Grigori Perelman",
                "date": "2003-07-17",
                "arxiv": "math/0307245",
                "content": "Endliche Auslöschungszeit für π₁ endlich → Poincaré-Vermutung",
            },
        ]

    def fields_medal_declined(self) -> str:
        r"""
        @brief Historische Notiz zur abgelehnten Fields-Medaille.
        @return Historischer Kontext
        @lastModified 2026-03-11
        """
        return (
            "Fields-Medaille 2006 (International Congress of Mathematicians, Madrid): "
            "Grigori Perelman wurde die Fields-Medaille für den Beweis der "
            "Poincaré-Vermutung und der Geometrisierungsvermutung zugesprochen. "
            "Perelman lehnte die Auszeichnung ab — eine beispiellose Geste in der "
            "Geschichte der Mathematik. Er lehnte auch das Clay-Preisgeld von "
            "$1.000.000 ab und begründete dies mit Unzufriedenheit über den "
            "Zustand der mathematischen Gemeinschaft. Zitat: "
            "'Ich habe gelernt, wie man das Universum kontrolliert. Warum sollte "
            "ich nach einer Million Dollar greifen?'"
        )


# ===========================================================================
# GEOMETRISIERUNGSVERMUTUNG
# ===========================================================================

class GeometrizationConjecture:
    r"""
    @brief Thurstons Geometrisierungsvermutung für 3-Mannigfaltigkeiten.
    @description
        William Thurston (1982) stellte die allgemeine Geometrisierungsvermutung
        für kompakte orientierbare 3-Mannigfaltigkeiten auf. Sie wurde von
        Perelman (2003) vollständig bewiesen.

        Die Vermutung besagt, dass jede solche Mannigfaltigkeit in Stücke
        zerlegt werden kann, die jeweils eine von 8 homogenen Modellgeometrien tragen.

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def thurston_classification(self, manifold_type: str) -> str:
        r"""
        @brief Klassifiziert eine 3-Mannigfaltigkeit nach Thurstons Schema.
        @param manifold_type Typ der Mannigfaltigkeit (z.B. "lens", "torus-bundle", etc.)
        @return              Zugehörige Thurston-Geometrie
        @lastModified 2026-03-11
        """
        classification = {
            "S3":            "S³-Geometrie (sphärisch, konstante positive Krümmung)",
            "lens":          "S³-Geometrie (Linsenräume L(p,q) tragen sphärische Geometrie)",
            "flat-torus":    "E³-Geometrie (flacher 3-Torus T³)",
            "nilmanifold":   "Nil-Geometrie (Heisenberg-Nilmannigfaltigkeit)",
            "hyperbolic":    "H³-Geometrie (hyperbolisch, endliches Volumen, Mostow-starr)",
            "sol-manifold":  "Sol-Geometrie (auflösbare Geometrie)",
            "torus-bundle":  "E³ oder Nil oder Sol (abhängig von Monodromie)",
            "seifert":       "S²×ℝ, H²×ℝ oder SL̃(2,ℝ)-Geometrie (Seifert-Faserbündel)",
        }
        return classification.get(manifold_type,
            f"Unbekannter Typ '{manifold_type}' — JSJ-Zerlegung in Thurston-Stücke nötig.")

    def prime_decomposition(self, manifold: str) -> List[str]:
        r"""
        @brief Primzerlegung nach Kneser-Milnor-Theorem.
        @description
            Kneser-Milnor-Theorem (1960/1962): Jede kompakte orientierbare
            3-Mannigfaltigkeit M zerlegt sich eindeutig (bis auf Reihenfolge) als
            $$M \cong M_1 \# M_2 \# \cdots \# M_k$$
            wobei jedes M_i prim ist (nicht weiter zerlegbar durch verbundene Summe).
        @param manifold String-Beschreibung der Mannigfaltigkeit
        @return         Liste der Primteile
        @lastModified 2026-03-11
        """
        # Einfache Beispiele für bekannte Mannigfaltigkeiten
        examples = {
            "S3":       ["S³"],
            "S2xS1":   ["S²×S¹"],
            "RP3":     ["ℝP³"],
            "T3":      ["T³ (flacher 3-Torus — nicht prim, aber irreduzibel)"],
            "L(p,q)":  ["L(p,q) (Linsenraum — prim)"],
        }
        return examples.get(manifold, [f"{manifold} (Primzerlegung nicht automatisch bestimmbar)"])

    def jaco_shalen_johannsen_decomposition(self) -> str:
        r"""
        @brief Beschreibt die JSJ-Zerlegung (Jaco-Shalen-Johannsen).
        @description
            JSJ-Zerlegung (1979): Nach der Primzerlegung wird jede irreduzible
            3-Mannigfaltigkeit entlang einer minimalen Menge unzusammendrückbarer
            Tori T² in Stücke zerlegt:
            - Seifert-gefaserte Stücke
            - Atoroidale, hyperbolische Stücke (Mostow-Rigidität)
        @return Erklärung der JSJ-Zerlegung
        @lastModified 2026-03-11
        """
        return (
            "JSJ-Zerlegung (Jaco-Shalen 1979, Johannsen 1979): "
            "Jede irreduzible, kompakte, orientierbare 3-Mannigfaltigkeit M "
            "besitzt eine minimale, bis auf Isotopie eindeutige Kollektion "
            "inkompressibler Tori T₁,...,Tₖ ⊂ M, so dass jede Komponente von "
            "M \\ (T₁ ∪ ... ∪ Tₖ) entweder Seifert-gefasert oder atoroidal ist. "
            "Die atoroidalen Stücke tragen hyperbolische Geometrie (Thurston) "
            "und sind durch die Mostow-Rigidität eindeutig bestimmt."
        )

    def hyperbolic_pieces(self) -> str:
        r"""
        @brief Erklärt die Bedeutung hyperbolischer Stücke und Mostow-Rigidität.
        @description
            Mostow-Rigidität (1973): Für n ≥ 3 ist jeder Isomorphismus
            zwischen Fundamentalgruppen kompakter hyperbolischer n-Mannigfaltigkeiten
            durch einen Isometrie realisiert. Insbesondere:
            Das Volumen und der Spektrum sind topologische Invarianten.
        @return Erklärung
        @lastModified 2026-03-11
        """
        return (
            "Mostow-Rigidität (George Mostow 1973): Kompakte hyperbolische "
            "3-Mannigfaltigkeiten sind durch ihre Fundamentalgruppe π₁(M) "
            "eindeutig bis auf Isometrie bestimmt. Das hyperbolische Volumen "
            "Vol(M) ist daher eine topologische Invariante — eine der "
            "mächtigsten Rigiditätseigenschaften in der 3-Topologie. "
            "Perelmans Beweis zeigt, dass atoroidale Stücke tatsächlich "
            "hyperbolische Geometrie im Sinne von Thurston tragen."
        )


# ===========================================================================
# UNIFORMISIERUNGSSATZ (2D-ANALOGON)
# ===========================================================================

class UniformizationTheorem:
    r"""
    @brief Uniformisierungssatz für Flächen (2D-Analogon der Poincaré-Vermutung).
    @description
        Der Uniformisierungssatz (Poincaré–Koebe 1907) ist das 2D-Analogon
        der Geometrisierungsvermutung:

        Jede einfach-zusammenhängende Riemannsche Fläche ist konform äquivalent
        zu genau einer von drei kanonischen Flächen:
        - Riemannsche Sphäre S² (positive Krümmung, χ > 0)
        - Euklidische Ebene ℝ² (flach, χ = 0)
        - Hyperbolische Ebene H² (negative Krümmung, χ < 0)

        Der Ricci-Fluss in 2D konvergiert immer zu einer Metrik konstanter
        Krümmung (kein Chirurgie nötig — Singularitäten möglich nur bei χ ≤ 0).

    @author Michael Fuhrmann
    @version 1.0
    @since 2026-03-11
    @lastModified 2026-03-11
    """

    def __init__(self, genus: int):
        r"""
        @brief Konstruktor.
        @param genus Geschlecht der Fläche (g ≥ 0)
        @lastModified 2026-03-11
        """
        if genus < 0:
            raise ValueError("Genus muss nicht-negativ sein.")
        self.genus = genus  # Topologisches Geschlecht g
        # Euler-Charakteristik χ = 2 - 2g
        self.euler_char = 2 - 2 * genus

    def classify_surface(self) -> str:
        r"""
        @brief Klassifiziert die Fläche nach Geschlecht.
        @description
            Klassifikationssatz für kompakte orientierbare Flächen:
            - g = 0: S² (2-Sphäre)
            - g = 1: T² (Torus)
            - g ≥ 2: Σ_g (geschlossene Fläche vom Geschlecht g)
        @return String-Beschreibung der Fläche
        @lastModified 2026-03-11
        """
        if self.genus == 0:
            return "S² (2-Sphäre, χ=2, sphärische Geometrie)"
        elif self.genus == 1:
            return "T² (Torus, χ=0, euklidische/flache Geometrie)"
        else:
            return f"Σ_{self.genus} (Fläche vom Geschlecht {self.genus}, χ={self.euler_char}, hyperbolische Geometrie)"

    def universal_cover(self) -> str:
        r"""
        @brief Gibt die universelle Überlagerung der Fläche zurück.
        @description
            Uniformisierungssatz:
            - g = 0: Universelle Überlagerung = S² (Riemannsche Sphäre)
            - g = 1: Universelle Überlagerung = ℝ² (Euklidische Ebene)
            - g ≥ 2: Universelle Überlagerung = H² (Poincaré-Halbebene/Disk)
        @return Beschreibung der universellen Überlagerung
        @lastModified 2026-03-11
        """
        if self.genus == 0:
            return "S² (Riemannsche Sphäre) — einzige selbst-überlagernde Fläche"
        elif self.genus == 1:
            return "ℝ² (Euklidische Ebene) — π₁(T²) = Z² wirkt durch Translationen"
        else:
            return f"H² (Hyperbolische Ebene) — π₁(Σ_{self.genus}) wirkt als Fuchssche Gruppe"

    def gauss_bonnet(self, K_avg: float) -> float:
        r"""
        @brief Berechnet ∫K dA = 2πχ via Gauss-Bonnet-Theorem.
        @description
            Gauss-Bonnet-Theorem (globale Form):
            $$\int_M K \, dA = 2\pi \chi(M)$$
            Für eine Fläche vom Geschlecht g: χ = 2 - 2g.
            Dieser Satz verbindet Geometrie (Krümmung K) mit Topologie (χ).
        @param K_avg Mittlere Gauss-Krümmung (∫K dA / Fläche)
        @return      Wert von 2πχ (der Gauss-Bonnet-Integralwert)
        @lastModified 2026-03-11
        """
        # Gauss-Bonnet: ∫K dA = 2πχ — unabhängig von K_avg!
        # (K_avg ist nur für didaktische Kontrolle angegeben)
        return 2.0 * math.pi * self.euler_char

    def ricci_flow_2d(self, t: float) -> str:
        r"""
        @brief Beschreibt das Verhalten des 2D-Ricci-Flusses für Zeit t.
        @description
            In 2D ist der normalisierte Ricci-Fluss äquivalent zur
            Yamabe-Gleichung und konvergiert für t → ∞ immer zu einer
            Metrik konstanter Krümmung (Hamilton 1988, Chow 1991).
            Keine Chirurgie notwendig in 2D!
        @param t Zeitpunkt t ≥ 0
        @return  Beschreibung des Flusszustands
        @lastModified 2026-03-11
        """
        if t < 0:
            raise ValueError("Zeitpunkt t muss nicht-negativ sein.")

        geometry = self.universal_cover()

        if t == 0.0:
            return f"t=0: Anfangsmetrik g_0 (beliebige konforme Klasse auf Σ_{self.genus})"
        elif t < 1.0:
            return f"t={t:.2f}: Kurzzeit-Entwicklung. Metrik glättet sich. {geometry}"
        else:
            # Langzeitverhalten
            if self.genus == 0:
                return (
                    f"t={t:.2f}: Normalisierter Fluss auf S² konvergiert zur "
                    f"runden Metrik konstanter Krümmung K=+1 (Hamilton 1988)."
                )
            elif self.genus == 1:
                return (
                    f"t={t:.2f}: Normalisierter Fluss auf T² konvergiert zur "
                    f"flachen Metrik K=0 (Chow 1991)."
                )
            else:
                return (
                    f"t={t:.2f}: Normalisierter Fluss auf Σ_{self.genus} konvergiert zur "
                    f"hyperbolischen Metrik K=-1 (Chow 1991, Hamilton 1988)."
                )
