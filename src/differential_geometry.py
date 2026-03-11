"""
@file differential_geometry.py
@brief Differentialgeometrie – Kurven, Flächen, Riemannsche Geometrie.
@description
    Implementiert grundlegende Konzepte der klassischen Differentialgeometrie:

    - ParametricCurve: Parametrisierte Kurven im ℝⁿ (Frenet-Serret-Theorie)
    - ParametricSurface: Parametrisierte Flächen im ℝ³ (1. und 2. Fundamentalform)
    - GeodesicComputation: Geodätische auf Flächen (Sphäre, allgemein)
    - RiemannianGeometry: Schnittkrümmung, Ricci-Tensor, Skalarkrümmung
    - ClassicalSurfaces: Sphäre, Torus, Zylinder, Helikoid, Katenoid, Sattelfläche
    - CurveTheory2D: Ebene Kurventheorie (Evolute, Evolvente, isoperimetrische Ungleichung)

    Wichtige Formeln:
    - Frenet-Serret: T' = κN, N' = -κT + τB, B' = -τN
    - Gaußsche Krümmung: K = (LN - M²) / (EG - F²)
    - Mittlere Krümmung: H = (EN + GL - 2FM) / (2(EG - F²))
    - Gauß-Bonnet: ∫∫_M K dA + ∮_∂M κ_g ds = 2πχ(M)

    Hinweis: MetricTensor, christoffel_symbols, riemann_tensor, gaussian_curvature
    (tensorielle Version), einstein_tensor und geodesic_equation sind bereits in
    tensor_geometry.py definiert und werden dort importiert.

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import math
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from typing import Callable, Tuple, List, Optional


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def _norm(v: np.ndarray) -> float:
    """
    Berechnet die euklidische Norm eines Vektors.

    @param v: Eingabevektor als numpy-Array.
    @return Euklidische Norm |v|.
    @lastModified 2026-03-10
    """
    return float(np.linalg.norm(v))


def _normalize(v: np.ndarray) -> np.ndarray:
    """
    Normiert einen Vektor auf Einheitslänge.

    @param v: Eingabevektor.
    @return Einheitsvektor v / |v|; Nullvektor bleibt Nullvektor.
    @lastModified 2026-03-10
    """
    n = _norm(v)
    # Nullvektor abfangen
    if n < 1e-15:
        return np.zeros_like(v, dtype=float)
    return v / n


def _cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Kreuzprodukt zweier 3D-Vektoren.

    @param a: Erster 3D-Vektor.
    @param b: Zweiter 3D-Vektor.
    @return Kreuzprodukt a × b.
    @lastModified 2026-03-10
    """
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ], dtype=float)


# ---------------------------------------------------------------------------
# Klasse: ParametricCurve
# ---------------------------------------------------------------------------

class ParametricCurve:
    """
    Parametrisierte Kurve im ℝⁿ, definiert durch SymPy-Ausdrücke.

    Eine parametrisierte Kurve ist eine glatte Abbildung
        r: I ⊂ ℝ → ℝⁿ, t ↦ r(t) = (r₁(t), ..., rₙ(t))

    Die Frenet-Serret-Formeln beschreiben die Geometrie der Kurve:
        T' = κ N
        N' = -κ T + τ B
        B' = -τ N

    wobei T = Tangenteneinheitsvektor, N = Hauptnormalenvektor,
    B = Binormalvektor, κ = Krümmung, τ = Torsion.

    Beispiel:
        >>> import sympy as sp
        >>> t = sp.Symbol('t')
        >>> # Helix: r(t) = (cos t, sin t, t)
        >>> curve = ParametricCurve([sp.cos(t), sp.sin(t), t], t, (0, 2*sp.pi))
        >>> abs(curve.curvature(1.0) - 0.5) < 1e-8
        True
    """

    def __init__(
        self,
        components: List[sp.Expr],
        parameter: sp.Symbol,
        t_range: Tuple[float, float] = (0, 1),
    ) -> None:
        """
        Initialisiert eine parametrisierte Kurve.

        @param components: Liste von SymPy-Ausdrücken [r₁(t), ..., rₙ(t)].
        @param parameter: SymPy-Symbol für den Parameter (z.B. t).
        @param t_range: Parameterintervall (t_min, t_max).
        @lastModified 2026-03-10
        """
        self.components = components          # Koordinatenfunktionen
        self.param = parameter                # Parametervariable
        self.t_range = t_range                # Definitionsbereich
        self.dim = len(components)            # Dimension des Raums

        # Erste und zweite Ableitung symbolisch berechnen
        self._velocity_sym = [sp.diff(c, parameter) for c in components]
        self._accel_sym = [sp.diff(c, parameter, 2) for c in components]
        # Dritte Ableitung (für Torsion in 3D)
        self._jerk_sym = [sp.diff(c, parameter, 3) for c in components]

        # Lambdify für numerische Auswertung
        self._r_fn = sp.lambdify(parameter, components, modules="numpy")
        self._v_fn = sp.lambdify(parameter, self._velocity_sym, modules="numpy")
        self._a_fn = sp.lambdify(parameter, self._accel_sym, modules="numpy")
        self._j_fn = sp.lambdify(parameter, self._jerk_sym, modules="numpy")

    def _eval(self, fn: Callable, t: float) -> np.ndarray:
        """
        Wertet eine lambdifizierte Funktion bei t aus und gibt numpy-Array zurück.

        @param fn: Lambdifizierte SymPy-Funktion.
        @param t: Parameterwert.
        @return Ergebnis als float64-Array.
        @lastModified 2026-03-10
        """
        result = fn(float(t))
        # Skalare in Listen einbetten
        if not hasattr(result, '__iter__') or isinstance(result, (int, float)):
            return np.array([float(result)], dtype=float)
        return np.array([float(x) for x in result], dtype=float)

    def velocity(self, t: float) -> np.ndarray:
        """
        Berechnet den Geschwindigkeitsvektor r'(t) an der Stelle t.

        @param t: Parameterwert.
        @return Geschwindigkeitsvektor als numpy-Array der Dimension n.
        @lastModified 2026-03-10
        """
        return self._eval(self._v_fn, t)

    def speed(self, t: float) -> float:
        """
        Berechnet die Bogenlängengeschwindigkeit |r'(t)| an der Stelle t.

        @param t: Parameterwert.
        @return Betrag der Geschwindigkeit (skalarer Wert ≥ 0).
        @lastModified 2026-03-10
        """
        return _norm(self.velocity(t))

    def arc_length(self, t0: float, t1: float) -> float:
        """
        Berechnet die Bogenlänge der Kurve von t0 bis t1.

        Formel: L = ∫_{t0}^{t1} |r'(t)| dt

        @param t0: Startwert des Parameters.
        @param t1: Endwert des Parameters.
        @return Bogenlänge (numerisch mit adaptiver Quadratur).
        @lastModified 2026-03-10
        """
        result, _ = quad(self.speed, float(t0), float(t1))
        return float(result)

    def curvature(self, t: float) -> float:
        """
        Berechnet die Krümmung κ(t) der Kurve an der Stelle t.

        Formel im ℝ³: κ = |r' × r''| / |r'|³
        Formel im ℝ²: κ = |x'y'' - y'x''| / (x'² + y'²)^{3/2}
        Allgemein: κ = √(|r''|² |r'|² - (r'·r'')²) / |r'|³

        @param t: Parameterwert.
        @return Krümmung κ ≥ 0.
        @lastModified 2026-03-10
        """
        v = self.velocity(t)        # r'(t)
        a = self._eval(self._a_fn, t)  # r''(t)
        speed_val = _norm(v)

        # Gerade (Nullgeschwindigkeit) oder reguläre Kurve
        if speed_val < 1e-15:
            return 0.0

        if self.dim == 3:
            # 3D: Kreuzprodukt r' × r''
            cross = _cross3(v, a)
            return _norm(cross) / (speed_val ** 3)
        else:
            # Allgemeine Dimension: Lagrange-Identität
            v2 = np.dot(v, v)
            a2 = np.dot(a, a)
            va = np.dot(v, a)
            # Betrag des verallgemeinerten Kreuzprodukts
            discriminant = a2 * v2 - va ** 2
            if discriminant < 0:
                discriminant = 0.0
            return math.sqrt(discriminant) / (speed_val ** 3)

    def torsion(self, t: float) -> float:
        """
        Berechnet die Torsion τ(t) der Kurve (nur in 3D definiert).

        Formel: τ = (r' × r'') · r''' / |r' × r''|²

        Die Torsion misst, wie stark die Kurve aus ihrer Schmiegebene herausragt.
        Für ebene Kurven ist τ = 0.

        @param t: Parameterwert.
        @return Torsion τ (kann negativ sein); 0.0 wenn nicht 3D.
        @lastModified 2026-03-10
        """
        if self.dim != 3:
            return 0.0

        v = self.velocity(t)
        a = self._eval(self._a_fn, t)
        j = self._eval(self._j_fn, t)

        # r' × r''
        cross = _cross3(v, a)
        cross_norm_sq = np.dot(cross, cross)

        # Entartung: gerade Kurve
        if cross_norm_sq < 1e-15:
            return 0.0

        # τ = (r' × r'') · r''' / |r' × r''|²
        return float(np.dot(cross, j)) / float(cross_norm_sq)

    def unit_tangent(self, t: float) -> np.ndarray:
        """
        Berechnet den Tangenteneinheitsvektor T(t) = r'(t) / |r'(t)|.

        @param t: Parameterwert.
        @return Normierter Tangentialvektor (Einheitsvektor).
        @lastModified 2026-03-10
        """
        return _normalize(self.velocity(t))

    def principal_normal(self, t: float) -> np.ndarray:
        """
        Berechnet den Hauptnormalenvektor N(t) = T'(t) / |T'(t)|.

        Der Hauptnormalenvektor zeigt in Richtung der Krümmung (zum Krümmungsmittelpunkt).

        @param t: Parameterwert.
        @return Normierter Hauptnormalenvektor; Nullvektor bei κ=0.
        @lastModified 2026-03-10
        """
        # Numerische Ableitung von T
        h = 1e-7
        T_plus = self.unit_tangent(t + h)
        T_minus = self.unit_tangent(t - h)
        T_prime = (T_plus - T_minus) / (2 * h)
        return _normalize(T_prime)

    def binormal(self, t: float) -> np.ndarray:
        """
        Berechnet den Binormalvektor B(t) = T(t) × N(t) (nur 3D).

        Der Binormalvektor steht senkrecht auf der Schmiegebene.

        @param t: Parameterwert.
        @return Binormaleinheitsvektor; Nullvektor wenn nicht 3D oder κ=0.
        @lastModified 2026-03-10
        """
        if self.dim != 3:
            return np.zeros(self.dim)
        T = self.unit_tangent(t)
        N = self.principal_normal(t)
        return _normalize(_cross3(T, N))

    def frenet_serret_frame(self, t: float) -> dict:
        """
        Berechnet den vollständigen Frenet-Serret-Rahmen {T, N, B} an der Stelle t.

        @param t: Parameterwert.
        @return Dictionary mit Schlüsseln 'T', 'N', 'B', 'kappa', 'tau'.
        @lastModified 2026-03-10
        """
        return {
            "T": self.unit_tangent(t),
            "N": self.principal_normal(t),
            "B": self.binormal(t),
            "kappa": self.curvature(t),
            "tau": self.torsion(t),
        }

    def frenet_serret_formulas(self) -> dict:
        """
        Gibt die symbolischen Frenet-Serret-Formeln zurück.

        Formeln:
            T' = κ N
            N' = -κ T + τ B
            B' = -τ N

        @return Dictionary mit Beschreibungen der Frenet-Serret-Gleichungen.
        @lastModified 2026-03-10
        """
        return {
            "T_prime": "κ · N",
            "N_prime": "-κ · T + τ · B",
            "B_prime": "-τ · N",
            "description": (
                "Frenet-Serret-Formeln: T'=κN, N'=-κT+τB, B'=-τN. "
                "κ = Krümmung ≥ 0, τ = Torsion ∈ ℝ."
            ),
        }


# ---------------------------------------------------------------------------
# Klasse: ParametricSurface
# ---------------------------------------------------------------------------

class ParametricSurface:
    """
    Parametrisierte Fläche im ℝ³, definiert durch eine callable r(u, v).

    Eine Fläche ist eine glatte Abbildung
        r: U ⊂ ℝ² → ℝ³, (u,v) ↦ r(u,v) = (x(u,v), y(u,v), z(u,v))

    Die intrinsische Geometrie wird durch die erste Fundamentalform beschrieben:
        ds² = E du² + 2F du dv + G dv²
    mit E = r_u · r_u, F = r_u · r_v, G = r_v · r_v.

    Die zweite Fundamentalform beschreibt die extrinsische Krümmung:
        II = L du² + 2M du dv + N dv²
    mit L = r_uu · n, M = r_uv · n, N = r_vv · n (n = Einheitsnormale).

    Beispiel:
        >>> import numpy as np
        >>> def sphere(u, v): return np.array([np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u)])
        >>> S = ParametricSurface(sphere, (0, np.pi), (0, 2*np.pi))
        >>> abs(S.gaussian_curvature(np.pi/2, 0) - 1.0) < 0.01
        True
    """

    def __init__(
        self,
        r_func: Callable[[float, float], np.ndarray],
        u_range: Tuple[float, float],
        v_range: Tuple[float, float],
    ) -> None:
        """
        Initialisiert eine parametrisierte Fläche.

        @param r_func: Callable r(u, v) → np.ndarray der Form (3,).
        @param u_range: Wertebereich für u als Tupel (u_min, u_max).
        @param v_range: Wertebereich für v als Tupel (v_min, v_max).
        @lastModified 2026-03-10
        """
        self.r = r_func
        self.u_range = u_range
        self.v_range = v_range

    def partial_derivatives(
        self, u: float, v: float, h: float = 1e-7
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Berechnet numerisch die partiellen Ableitungen r_u und r_v.

        Formel (zentraler Differenzenquotient):
            r_u ≈ (r(u+h, v) - r(u-h, v)) / (2h)
            r_v ≈ (r(u, v+h) - r(u, v-h)) / (2h)

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @param h: Schrittweite für numerische Ableitung.
        @return Tupel (r_u, r_v) als numpy-Arrays.
        @lastModified 2026-03-10
        """
        # Partielle Ableitung nach u
        r_u = (np.array(self.r(u + h, v), dtype=float)
               - np.array(self.r(u - h, v), dtype=float)) / (2 * h)
        # Partielle Ableitung nach v
        r_v = (np.array(self.r(u, v + h), dtype=float)
               - np.array(self.r(u, v - h), dtype=float)) / (2 * h)
        return r_u, r_v

    def _second_partials(
        self, u: float, v: float, h: float = 1e-5
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Berechnet numerisch die zweiten partiellen Ableitungen r_uu, r_uv, r_vv.

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @param h: Schrittweite.
        @return Tupel (r_uu, r_uv, r_vv).
        @lastModified 2026-03-10
        """
        r0 = np.array(self.r(u, v), dtype=float)
        # Zweite partielle Ableitung nach u
        r_uu = (np.array(self.r(u + h, v), dtype=float)
                - 2 * r0
                + np.array(self.r(u - h, v), dtype=float)) / (h ** 2)
        # Zweite partielle Ableitung nach v
        r_vv = (np.array(self.r(u, v + h), dtype=float)
                - 2 * r0
                + np.array(self.r(u, v - h), dtype=float)) / (h ** 2)
        # Gemischte partielle Ableitung
        r_uv = (np.array(self.r(u + h, v + h), dtype=float)
                - np.array(self.r(u + h, v - h), dtype=float)
                - np.array(self.r(u - h, v + h), dtype=float)
                + np.array(self.r(u - h, v - h), dtype=float)) / (4 * h ** 2)
        return r_uu, r_uv, r_vv

    def normal_vector(self, u: float, v: float) -> np.ndarray:
        """
        Berechnet den Normalenvektor n = r_u × r_v (nicht normiert).

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Normalenvektor (nicht einheitlich).
        @lastModified 2026-03-10
        """
        r_u, r_v = self.partial_derivatives(u, v)
        return _cross3(r_u, r_v)

    def unit_normal(self, u: float, v: float) -> np.ndarray:
        """
        Berechnet den normierten Einheitsnormalenvektor n / |n|.

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Einheitsnormalenvektor.
        @lastModified 2026-03-10
        """
        return _normalize(self.normal_vector(u, v))

    def first_fundamental_form(self, u: float, v: float) -> Tuple[float, float, float]:
        """
        Berechnet die Koeffizienten der ersten Fundamentalform E, F, G.

        Formeln:
            E = r_u · r_u (Metrikkoeffizient in u-Richtung)
            F = r_u · r_v (Kreuzterm; F=0 bei orthogonalen Koordinaten)
            G = r_v · r_v (Metrikkoeffizient in v-Richtung)

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Tupel (E, F, G).
        @lastModified 2026-03-10
        """
        r_u, r_v = self.partial_derivatives(u, v)
        E = float(np.dot(r_u, r_u))
        F = float(np.dot(r_u, r_v))
        G = float(np.dot(r_v, r_v))
        return E, F, G

    def second_fundamental_form(self, u: float, v: float) -> Tuple[float, float, float]:
        """
        Berechnet die Koeffizienten der zweiten Fundamentalform L, M, N.

        Formeln (Weingarten):
            L = r_uu · n̂  (Normalenkrümmung in u-Richtung)
            M = r_uv · n̂  (gemischter Term)
            N = r_vv · n̂  (Normalenkrümmung in v-Richtung)

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Tupel (L, M, N_coeff).
        @lastModified 2026-03-10
        """
        n_hat = self.unit_normal(u, v)
        r_uu, r_uv, r_vv = self._second_partials(u, v)
        L = float(np.dot(r_uu, n_hat))
        M = float(np.dot(r_uv, n_hat))
        N_coeff = float(np.dot(r_vv, n_hat))
        return L, M, N_coeff

    def gaussian_curvature(self, u: float, v: float) -> float:
        """
        Berechnet die Gaußsche Krümmung K der Fläche an (u, v).

        Formel: K = (LN - M²) / (EG - F²)

        Die Gaußsche Krümmung ist eine intrinsische Größe (Theorema Egregium).
        Sphäre: K > 0, Ebene/Zylinder: K = 0, Sattelfläche: K < 0.

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Gaußsche Krümmung K.
        @lastModified 2026-03-10
        """
        E, F, G = self.first_fundamental_form(u, v)
        L, M, N_coeff = self.second_fundamental_form(u, v)
        denom = E * G - F ** 2
        if abs(denom) < 1e-15:
            return 0.0
        return (L * N_coeff - M ** 2) / denom

    def mean_curvature(self, u: float, v: float) -> float:
        """
        Berechnet die mittlere Krümmung H der Fläche an (u, v).

        Formel: H = (EN + GL - 2FM) / (2(EG - F²))

        H = 0 kennzeichnet Minimalflächen (z.B. Seifenblasen-Formen).

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Mittlere Krümmung H.
        @lastModified 2026-03-10
        """
        E, F, G = self.first_fundamental_form(u, v)
        L, M, N_coeff = self.second_fundamental_form(u, v)
        denom = 2 * (E * G - F ** 2)
        if abs(denom) < 1e-15:
            return 0.0
        return (E * N_coeff + G * L - 2 * F * M) / denom

    def principal_curvatures(self, u: float, v: float) -> Tuple[float, float]:
        """
        Berechnet die Hauptkrümmungen κ₁ und κ₂ (Eigenwerte der Formmatrix).

        Die Hauptkrümmungen sind die Extremalwerte der Normalkrümmung.
        Es gilt: K = κ₁ · κ₂ und H = (κ₁ + κ₂) / 2.

        Eigenwertproblem der Weingarten-Abbildung:
            [L  M] [du]   [E  F] [du]
            [M  N] [dv] = κ [F  G] [dv]

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Tupel (κ₁, κ₂) mit κ₁ ≤ κ₂.
        @lastModified 2026-03-10
        """
        E, F, G = self.first_fundamental_form(u, v)
        L, M, N_coeff = self.second_fundamental_form(u, v)

        # Charakteristisches Polynom: (EG-F²)κ² - (EN+GL-2FM)κ + (LN-M²) = 0
        A = E * G - F ** 2
        B_coeff = -(E * N_coeff + G * L - 2 * F * M)
        C = L * N_coeff - M ** 2

        # Quadratische Formel
        discriminant = B_coeff ** 2 - 4 * A * C
        if discriminant < 0:
            discriminant = 0.0
        sqrt_disc = math.sqrt(discriminant)
        if abs(A) < 1e-15:
            return 0.0, 0.0
        kappa1 = (-B_coeff - sqrt_disc) / (2 * A)
        kappa2 = (-B_coeff + sqrt_disc) / (2 * A)
        return min(kappa1, kappa2), max(kappa1, kappa2)

    def area_element(self, u: float, v: float) -> float:
        """
        Berechnet das Flächenelement dA = |r_u × r_v| an der Stelle (u, v).

        @param u: Parameterwert u.
        @param v: Parameterwert v.
        @return Flächenelement (positiver Skalar).
        @lastModified 2026-03-10
        """
        return _norm(self.normal_vector(u, v))

    def total_area(
        self,
        u_range: Optional[Tuple[float, float]] = None,
        v_range: Optional[Tuple[float, float]] = None,
        n_pts: int = 20,
    ) -> float:
        """
        Berechnet die Gesamtfläche der Fläche numerisch (Gaußsche Quadratur).

        A = ∫∫ |r_u × r_v| du dv

        @param u_range: Integrationsbereich für u (Standard: self.u_range).
        @param v_range: Integrationsbereich für v (Standard: self.v_range).
        @param n_pts: Anzahl Stützstellen pro Dimension (Gaußsche Quadratur).
        @return Approximierte Gesamtfläche.
        @lastModified 2026-03-10
        """
        if u_range is None:
            u_range = self.u_range
        if v_range is None:
            v_range = self.v_range

        # Gaußsche Legendre-Quadratur in 2D
        u_nodes, u_weights = np.polynomial.legendre.leggauss(n_pts)
        v_nodes, v_weights = np.polynomial.legendre.leggauss(n_pts)

        # Transformation von [-1,1] → [u_min, u_max]
        u_mid = 0.5 * (u_range[0] + u_range[1])
        u_half = 0.5 * (u_range[1] - u_range[0])
        v_mid = 0.5 * (v_range[0] + v_range[1])
        v_half = 0.5 * (v_range[1] - v_range[0])

        total = 0.0
        for i, (xi, wi) in enumerate(zip(u_nodes, u_weights)):
            u_val = u_mid + u_half * xi
            for j, (xj, wj) in enumerate(zip(v_nodes, v_weights)):
                v_val = v_mid + v_half * xj
                # Flächenelement an diesem Punkt
                dA = self.area_element(u_val, v_val)
                total += wi * wj * dA

        return float(total * u_half * v_half)


# ---------------------------------------------------------------------------
# Klasse: GeodesicComputation
# ---------------------------------------------------------------------------

class GeodesicComputation:
    """
    Berechnung von Geodätischen auf Flächen und in Riemannschen Räumen.

    Eine Geodätische ist eine verallgemeinerte "gerade Linie" auf einer
    gekrümmten Fläche – die kürzeste Verbindung zweier Punkte.

    Geodätengleichung:
        d²xᵏ/dt² + Γᵏᵢⱼ (dxⁱ/dt)(dxʲ/dt) = 0

    Auf der Sphäre S² sind die Geodätischen Großkreisbögen.

    Beispiel:
        >>> gc = GeodesicComputation()
        >>> result = gc.geodesic_sphere(np.pi/2, 0.0, [1.0, 0.0], t_max=np.pi)
        >>> result['success']
        True
    """

    def geodesic_sphere(
        self,
        theta0: float,
        phi0: float,
        direction: List[float],
        t_max: float = 1.0,
    ) -> dict:
        """
        Berechnet eine Geodätische auf der Einheitssphäre S².

        Auf S² sind Geodätische Großkreise. Die Geodätengleichung lautet in
        Kugelkoordinaten (θ, φ):
            θ'' = sin(θ)cos(θ) φ'²
            φ'' = -2 cot(θ) θ' φ'

        @param theta0: Startbreite θ₀ ∈ [0, π].
        @param phi0: Startlänge φ₀ ∈ [0, 2π].
        @param direction: Anfangsgeschwindigkeit [θ', φ'] in Tangentialkoordinaten.
        @param t_max: Endzeitpunkt der Integration.
        @return Dictionary mit 'theta', 'phi', 't', 'success'.
        @lastModified 2026-03-10
        """
        # Geodätengleichung auf S² in (θ, φ)-Koordinaten
        def geodesic_ode(t, y):
            theta, phi, dtheta, dphi = y
            # Christoffel-Symbole der Einheitssphäre:
            # Γ^θ_{φφ} = -sin(θ)cos(θ), Γ^φ_{θφ} = Γ^φ_{φθ} = cot(θ)
            ddtheta = math.sin(theta) * math.cos(theta) * dphi ** 2
            # Schutz vor θ ≈ 0 oder π (Polsingularität)
            if abs(math.sin(theta)) < 1e-12:
                ddphi = 0.0
            else:
                ddphi = -2.0 * (math.cos(theta) / math.sin(theta)) * dtheta * dphi
            return [dtheta, dphi, ddtheta, ddphi]

        # Anfangsbedingungen: [θ₀, φ₀, θ'₀, φ'₀]
        y0 = [theta0, phi0, float(direction[0]), float(direction[1])]
        sol = solve_ivp(
            geodesic_ode, [0, t_max], y0,
            method="RK45", dense_output=True, rtol=1e-10, atol=1e-12
        )
        return {
            "theta": sol.y[0],
            "phi": sol.y[1],
            "t": sol.t,
            "success": sol.success,
        }

    def geodesic_equations_general(
        self,
        metric_fn: Callable[[List[float]], np.ndarray],
        initial_pos: List[float],
        initial_vel: List[float],
        t_max: float = 1.0,
    ) -> dict:
        """
        Löst die allgemeine Geodätengleichung numerisch für eine gegebene Metrik.

        Geodätengleichung:
            d²xᵏ/dt² + Γᵏᵢⱼ (dxⁱ/dt)(dxʲ/dt) = 0

        Die Christoffel-Symbole werden numerisch aus der Metrik berechnet.

        @param metric_fn: Funktion g(x) → n×n-Matrix (Metriktensor).
        @param initial_pos: Startposition im Koordinatenraum.
        @param initial_vel: Startgeschwindigkeit (Tangentialvektor).
        @param t_max: Endzeitpunkt der Integration.
        @return Dictionary mit 'positions', 'velocities', 't', 'success'.
        @lastModified 2026-03-10
        """
        n = len(initial_pos)
        h = 1e-5  # Schrittweite für numerische Christoffel-Symbole

        def christoffel_numerical(x: np.ndarray) -> np.ndarray:
            """Berechnet Christoffel-Symbole Γᵏᵢⱼ numerisch."""
            g = np.array(metric_fn(list(x)), dtype=float)
            g_inv = np.linalg.inv(g)
            Gamma = np.zeros((n, n, n))
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        # ∂g_{ij}/∂xᵏ, ∂g_{ik}/∂xʲ, ∂g_{jk}/∂xⁱ
                        xi_plus = x.copy(); xi_plus[i] += h
                        xi_minus = x.copy(); xi_minus[i] -= h
                        xj_plus = x.copy(); xj_plus[j] += h
                        xj_minus = x.copy(); xj_minus[j] -= h
                        xk_plus = x.copy(); xk_plus[k] += h
                        xk_minus = x.copy(); xk_minus[k] -= h

                        dg_kij = (np.array(metric_fn(list(xk_plus)))[i, j]
                                  - np.array(metric_fn(list(xk_minus)))[i, j]) / (2 * h)
                        dg_jik = (np.array(metric_fn(list(xj_plus)))[i, k]
                                  - np.array(metric_fn(list(xj_minus)))[i, k]) / (2 * h)
                        dg_ivjk = (np.array(metric_fn(list(xi_plus)))[j, k]
                                   - np.array(metric_fn(list(xi_minus)))[j, k]) / (2 * h)

                        # Γ_{ij,k} = ½(∂_i g_{jk} + ∂_j g_{ik} - ∂_k g_{ij})
                        Gamma_lower = 0.5 * (dg_jik + dg_ivjk - dg_kij)
                        # Hochziehen: Γᵏᵢⱼ = gᵏˡ Γ_{ij,l}
                        for l in range(n):
                            Gamma[k, i, j] += g_inv[k, l] * Gamma_lower

            return Gamma

        def ode_system(t, y):
            """Geodätengleichung als ODE-System erster Ordnung."""
            x = y[:n]
            v = y[n:]
            Gamma = christoffel_numerical(np.array(x))
            # d²xᵏ/dt² = -Γᵏᵢⱼ vⁱ vʲ
            accel = np.zeros(n)
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        accel[k] -= Gamma[k, i, j] * v[i] * v[j]
            return list(v) + list(accel)

        y0 = list(initial_pos) + list(initial_vel)
        sol = solve_ivp(
            ode_system, [0, t_max], y0,
            method="RK45", rtol=1e-8, atol=1e-10
        )
        return {
            "positions": sol.y[:n].T,
            "velocities": sol.y[n:].T,
            "t": sol.t,
            "success": sol.success,
        }

    def conjugate_points_demo(self, surface_fn: str = "sphere") -> dict:
        """
        Demonstration konjugierter Punkte auf einer Fläche.

        Konjugierte Punkte sind Punkte, bei denen benachbarte Geodätische
        sich wieder treffen. Auf der Sphäre sind Nord- und Südpol konjugiert.

        @param surface_fn: Art der Fläche ('sphere', 'torus', 'plane').
        @return Dictionary mit Beschreibung und Beispiel.
        @lastModified 2026-03-10
        """
        examples = {
            "sphere": {
                "description": (
                    "Auf der Einheitssphäre S² sind zwei Punkte P und Q "
                    "konjugiert, wenn Q diametral gegenüber von P liegt. "
                    "Alle Meridiane (Großkreise) von Nordpol nach Südpol "
                    "treffen sich wieder im Südpol."
                ),
                "example": {
                    "P": [0, 0, 1],   # Nordpol
                    "Q": [0, 0, -1],  # Südpol
                    "distance": math.pi,
                },
                "focal_point_distance": math.pi,
            },
            "torus": {
                "description": (
                    "Auf einem Torus hängen konjugierte Punkte von der "
                    "Geodätischen ab. Durch die variable Krümmung entstehen "
                    "konjugierte Punkte in endlichem Abstand."
                ),
                "example": None,
            },
            "plane": {
                "description": (
                    "Auf der Ebene (K=0) gibt es keine konjugierten Punkte. "
                    "Geodätische (gerade Linien) divergieren immer."
                ),
                "example": None,
            },
        }
        return examples.get(surface_fn, examples["sphere"])

    def exponential_map_surface(
        self, p: List[float], v: List[float], t: float
    ) -> np.ndarray:
        """
        Approximiert die Exponentialabbildung exp_p(t·v) auf einer Fläche.

        Die Exponentialabbildung sendet einen Tangentialvektor v ∈ T_p M
        auf den Punkt exp_p(v) = γ(1), wobei γ die Geodätische mit γ(0)=p,
        γ'(0)=v ist.

        Diese Implementierung arbeitet im euklidischen ℝ³ mit linearer Approximation
        (erste Ordnung) und normiert auf die Sphäre für eine geometrisch sinnvolle Demo.

        @param p: Startpunkt p ∈ ℝ³ (sollte auf der Sphäre liegen).
        @param v: Tangentialvektor v ∈ T_p S².
        @param t: Zeitparameter t.
        @return Punkt auf der Sphäre exp_p(t·v).
        @lastModified 2026-03-10
        """
        p = np.array(p, dtype=float)
        v = np.array(v, dtype=float)
        # Linearer Schritt, dann Projektion auf Sphäre (erste Approximation)
        q = p + t * v
        norm_q = _norm(q)
        if norm_q < 1e-15:
            return p.copy()
        # Für Einheitssphäre: normiere das Ergebnis
        return q / norm_q * _norm(p)


# ---------------------------------------------------------------------------
# Klasse: RiemannianGeometry
# ---------------------------------------------------------------------------

class RiemannianGeometry:
    """
    Werkzeuge für Riemannsche Geometrie.

    Riemannsche Geometrie verallgemeinert die euklidische Geometrie auf
    gekrümmte Räume. Zentrale Objekte:
    - Metriktensor g_{ij}: definiert Abstände und Winkel
    - Riemann-Krümmungstensor R^l_{ijk}: misst Nichtkommutativität des Paralleltransports
    - Ricci-Tensor Ric_{ij} = R^k_{ikj}: Spur des Riemann-Tensors
    - Skalarkrümmung R = g^{ij} Ric_{ij}: vollständige Spur

    Hinweis: riemann_tensor() aus tensor_geometry.py wird hier für konkrete
    Berechnungen importiert/verwendet. Diese Klasse bietet höherstufige Operationen.

    Beispiel:
        >>> rg = RiemannianGeometry()
        >>> R = np.zeros((3, 3, 3, 3))  # flacher Raum
        >>> metric = np.eye(3)
        >>> rg.scalar_curvature(rg.ricci_curvature(R, metric), metric)
        0.0
    """

    def sectional_curvature(
        self,
        R_tensor: np.ndarray,
        v1: np.ndarray,
        v2: np.ndarray,
    ) -> float:
        """
        Berechnet die Schnittkrümmung K(σ) für den durch v1, v2 aufgespannten Schnitt σ.

        Formel: K(σ) = R(v1, v2, v1, v2) / (|v1|²|v2|² - (v1·v2)²)

        Die Schnittkrümmung verallgemeinert die Gaußsche Krümmung auf höhere Dimensionen.

        @param R_tensor: Riemannscher Krümmungstensor R^l_{ijk} als (n,n,n,n)-Array.
        @param v1: Erster Tangentialvektor.
        @param v2: Zweiter Tangentialvektor.
        @return Schnittkrümmung K(σ) (skalarer Wert).
        @lastModified 2026-03-10
        """
        v1 = np.array(v1, dtype=float)
        v2 = np.array(v2, dtype=float)
        n = len(v1)

        # R(v1, v2, v1, v2) = R_{ijkl} v1^i v2^j v1^k v2^l
        # Für R^l_{ijk}: R_{ijkl} = g_{lm} R^m_{ijk} (hier ohne Metrik: Euklidisch)
        numerator = 0.0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        # Konvention: R_{l,ijk} = R_tensor[l,i,j,k]
                        numerator += R_tensor[l, i, j, k] * v1[i] * v2[j] * v1[k] * v2[l]

        # Gram-Determinante: |v1|²|v2|² - (v1·v2)²
        denominator = np.dot(v1, v1) * np.dot(v2, v2) - np.dot(v1, v2) ** 2
        if abs(denominator) < 1e-15:
            return 0.0
        return float(numerator / denominator)

    def ricci_curvature(
        self, R_tensor: np.ndarray, metric: np.ndarray
    ) -> np.ndarray:
        """
        Berechnet den Ricci-Tensor Ric_{ij} = R^k_{ikj}.

        Der Ricci-Tensor ist eine Spur des Riemann-Tensors und beschreibt
        die Volumenverzerrung durch Krümmung (relevant für Einsteingleichungen).

        @param R_tensor: Riemannscher Krümmungstensor als (n,n,n,n)-Array (R^l_{ijk}).
        @param metric: Metriktensor g_{ij} als (n,n)-Matrix.
        @return Ricci-Tensor als (n,n)-Matrix.
        @lastModified 2026-03-10
        """
        R = np.array(R_tensor, dtype=float)
        n = R.shape[0]
        ricci = np.zeros((n, n))

        # Ric_{ij} = R^k_{ikj} = Spur über ersten und dritten Index
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    # R_tensor[k, i, k, j] = R^k_{ikj}
                    ricci[i, j] += R[k, i, k, j]

        return ricci

    def scalar_curvature(
        self, ricci: np.ndarray, metric: np.ndarray
    ) -> float:
        """
        Berechnet die Skalarkrümmung R = g^{ij} Ric_{ij}.

        Die Skalarkrümmung ist eine einzige reelle Zahl, die die Gesamt-
        krümmung des Raums an einem Punkt zusammenfasst.

        @param ricci: Ricci-Tensor Ric_{ij} als (n,n)-Matrix.
        @param metric: Metriktensor g_{ij} als (n,n)-Matrix.
        @return Skalarkrümmung (reelle Zahl).
        @lastModified 2026-03-10
        """
        g = np.array(metric, dtype=float)
        Ric = np.array(ricci, dtype=float)
        # Inverse Metrik g^{ij}
        g_inv = np.linalg.inv(g)
        # R = g^{ij} Ric_{ij} = Spur(g⁻¹ · Ric)
        return float(np.trace(g_inv @ Ric))

    def covariant_derivative_demo(
        self,
        metric_fn: Callable[[List[float]], np.ndarray],
        vector_field: np.ndarray,
    ) -> dict:
        """
        Demonstriert die kovariante Ableitung eines Vektorfeldes.

        Die kovariante Ableitung ∇_i V^j = ∂_i V^j + Γ^j_{ik} V^k
        ersetzt die gewöhnliche Ableitung in der Riemannschen Geometrie.

        @param metric_fn: Metrikfunktion g(x) → (n,n)-Matrix.
        @param vector_field: Konstantes Vektorfeld V^j (Demo: an einem Punkt).
        @return Dictionary mit Beschreibung und formaler Darstellung.
        @lastModified 2026-03-10
        """
        n = len(vector_field)
        return {
            "description": (
                "Kovariante Ableitung ∇_i V^j = ∂_i V^j + Γ^j_{ik} V^k. "
                "Sie berücksichtigt die Krümmung des Raums beim Ableiten."
            ),
            "formula": "∇_i V^j = ∂_i V^j + Γ^j_{ik} V^k",
            "vector_field": vector_field.tolist() if hasattr(vector_field, 'tolist') else list(vector_field),
            "dimension": n,
        }

    def parallel_transport_demo(
        self,
        curve: np.ndarray,
        vector: np.ndarray,
        metric: np.ndarray,
    ) -> dict:
        """
        Demonstration des Paralleltransports entlang einer Kurve.

        Beim Paralleltransport wird ein Vektor V längs einer Kurve γ so
        transportiert, dass ∇_{γ'} V = 0. Auf gekrümmten Flächen dreht
        sich der Vektor nach einem vollständigen Umlauf um den Holonomiewinkel.

        @param curve: Kurve γ als (N,n)-Array von Punkten.
        @param vector: Anfangsvektor V₀.
        @param metric: Metriktensor (konstant, für Demo).
        @return Dictionary mit Beschreibung und transportiertem Vektor.
        @lastModified 2026-03-10
        """
        # Demo-Implementierung: Paralleltransport auf flachem Raum (trivial)
        curve = np.array(curve, dtype=float)
        vector = np.array(vector, dtype=float)
        transported = vector.copy()  # auf flachem Raum bleibt Vektor konstant

        return {
            "description": (
                "Paralleltransport: ∇_{γ'} V = 0 längs γ. "
                "Auf gekrümmten Flächen entsteht nach geschlossenem Umlauf "
                "der Holonomiewinkel θ = ∫∫ K dA (Gauß-Bonnet)."
            ),
            "initial_vector": vector.tolist(),
            "transported_vector": transported.tolist(),
            "holonomy_formula": "θ = ∬_D K dA",
            "curve_length": float(np.sum(np.linalg.norm(np.diff(curve, axis=0), axis=1))),
        }

    def gauss_bonnet_theorem_demo(self, surface_type: str = "sphere") -> dict:
        """
        Demonstration des Gauß-Bonnet-Satzes für klassische Flächen.

        Gauß-Bonnet-Satz (ohne Rand):
            ∫∫_M K dA = 2π χ(M)

        wobei χ(M) die Euler-Charakteristik der Fläche ist:
            χ = V - E + F (Euler-Formel)

        Klassische Werte:
        - Sphäre: χ = 2, ∫K dA = 4π (K = 1/r²)
        - Torus: χ = 0, ∫K dA = 0 (positive und negative K heben sich auf)
        - Ebene/Zylinder: χ = 0, K = 0

        @param surface_type: 'sphere', 'torus', 'plane', 'genus_g'.
        @return Dictionary mit χ, ∫K dA und Erläuterung.
        @lastModified 2026-03-10
        """
        data = {
            "sphere": {
                "euler_characteristic": 2,
                "integral_K": 4 * math.pi,
                "expected_2pi_chi": 2 * math.pi * 2,
                "K_constant": 1.0,
                "description": "S² hat χ=2, K=1 (Einheitssphäre), ∫K dA = 4π = 2πχ.",
                "satisfied": True,
            },
            "torus": {
                "euler_characteristic": 0,
                "integral_K": 0.0,
                "expected_2pi_chi": 0.0,
                "K_constant": None,
                "description": "Torus hat χ=0. K ist variabel (positiv außen, negativ innen), ∫K dA = 0.",
                "satisfied": True,
            },
            "plane": {
                "euler_characteristic": 0,
                "integral_K": 0.0,
                "expected_2pi_chi": 0.0,
                "K_constant": 0.0,
                "description": "Ebene hat χ=0, K=0 (flach).",
                "satisfied": True,
            },
            "genus_g": {
                "euler_characteristic": None,
                "integral_K": None,
                "expected_2pi_chi": None,
                "K_constant": None,
                "description": "Fläche vom Geschlecht g hat χ = 2-2g, ∫K dA = 2π(2-2g).",
                "formula": "χ(Σ_g) = 2 - 2g, ∫∫ K dA = 4π(1-g)",
                "satisfied": True,
            },
        }
        return data.get(surface_type, data["sphere"])


# ---------------------------------------------------------------------------
# Klasse: ClassicalSurfaces
# ---------------------------------------------------------------------------

class ClassicalSurfaces:
    """
    Erzeugt klassische Flächen als ParametricSurface-Objekte.

    Beinhaltet: Sphäre, Torus, Zylinder, Helikoid, Katenoid, Sattelfläche
    sowie eine Methode zur Prüfung, ob eine Fläche eine Minimalfläche ist (H ≈ 0).

    Beispiel:
        >>> cs = ClassicalSurfaces()
        >>> sphere = cs.sphere(radius=1.0)
        >>> abs(sphere.gaussian_curvature(np.pi/2, 0.0) - 1.0) < 0.05
        True
    """

    def sphere(self, radius: float = 1.0) -> ParametricSurface:
        """
        Erzeugt eine Sphäre mit gegebenem Radius als ParametricSurface.

        Parametrisierung:
            r(θ, φ) = (R sin θ cos φ, R sin θ sin φ, R cos θ)
            θ ∈ [0, π], φ ∈ [0, 2π]

        Gaußsche Krümmung: K = 1/R²
        Mittlere Krümmung: H = 1/R

        @param radius: Kugelradius R > 0.
        @return ParametricSurface-Objekt der Sphäre.
        @lastModified 2026-03-10
        """
        R = float(radius)

        def r_sphere(theta: float, phi: float) -> np.ndarray:
            """Kugelkoordinaten nach kartesisch."""
            return np.array([
                R * math.sin(theta) * math.cos(phi),
                R * math.sin(theta) * math.sin(phi),
                R * math.cos(theta),
            ])

        return ParametricSurface(r_sphere, (0, math.pi), (0, 2 * math.pi))

    def torus(self, R: float = 2.0, r: float = 1.0) -> ParametricSurface:
        """
        Erzeugt einen Torus mit großem Radius R und kleinem Radius r.

        Parametrisierung:
            x = (R + r cos v) cos u
            y = (R + r cos v) sin u
            z = r sin v
            u, v ∈ [0, 2π]

        Gaußsche Krümmung: K = cos v / (r(R + r cos v))

        @param R: Großer Radius (Abstand Achse → Rohrenmitte).
        @param r: Kleiner Radius (Rohrradius).
        @return ParametricSurface-Objekt des Torus.
        @lastModified 2026-03-10
        """
        def r_torus(u: float, v: float) -> np.ndarray:
            """Torusparametrisierung."""
            return np.array([
                (R + r * math.cos(v)) * math.cos(u),
                (R + r * math.cos(v)) * math.sin(u),
                r * math.sin(v),
            ])

        return ParametricSurface(r_torus, (0, 2 * math.pi), (0, 2 * math.pi))

    def cylinder(self, radius: float = 1.0, height: float = 2.0) -> ParametricSurface:
        """
        Erzeugt einen Zylinder (Gaußsche Krümmung K = 0, da lokal flach).

        Parametrisierung:
            r(u, v) = (R cos u, R sin u, v)
            u ∈ [0, 2π], v ∈ [0, H]

        Der Zylinder ist lokal isometrisch zur Ebene.

        @param radius: Zylinderradius R.
        @param height: Zylinderhöhe H.
        @return ParametricSurface-Objekt des Zylinders.
        @lastModified 2026-03-10
        """
        R = float(radius)

        def r_cylinder(u: float, v: float) -> np.ndarray:
            """Zylinderparametrisierung."""
            return np.array([R * math.cos(u), R * math.sin(u), v])

        return ParametricSurface(r_cylinder, (0, 2 * math.pi), (0, float(height)))

    def helicoid(self) -> ParametricSurface:
        """
        Erzeugt einen Helikoid (klassische Minimalfläche).

        Parametrisierung:
            r(u, v) = (v cos u, v sin u, u)
            u ∈ [0, 2π], v ∈ [-1, 1]

        Der Helikoid ist eine Minimalfläche: H = 0 überall.
        Gauß-Krümmung: K = -1 / (1 + v²)²

        @return ParametricSurface-Objekt des Helikoids.
        @lastModified 2026-03-10
        """
        def r_helicoid(u: float, v: float) -> np.ndarray:
            """Helikoidparametrisierung."""
            return np.array([v * math.cos(u), v * math.sin(u), u])

        return ParametricSurface(r_helicoid, (-math.pi, math.pi), (-1.0, 1.0))

    def catenoid(self) -> ParametricSurface:
        """
        Erzeugt ein Katenoid (klassische Minimalfläche, Rotationsfläche).

        Parametrisierung:
            r(u, v) = (cosh(v) cos u, cosh(v) sin u, v)
            u ∈ [0, 2π], v ∈ [-2, 2]

        Das Katenoid ist die einzige Rotationsminimalfläche (neben der Ebene).
        H = 0, K = -sech⁴(v).

        @return ParametricSurface-Objekt des Katenoids.
        @lastModified 2026-03-10
        """
        def r_catenoid(u: float, v: float) -> np.ndarray:
            """Katenoidparametrisierung."""
            return np.array([
                math.cosh(v) * math.cos(u),
                math.cosh(v) * math.sin(u),
                v,
            ])

        return ParametricSurface(r_catenoid, (0, 2 * math.pi), (-2.0, 2.0))

    def saddle_surface(self) -> ParametricSurface:
        """
        Erzeugt eine Sattelfläche z = x² - y² (hyperbolisches Paraboloid).

        Parametrisierung:
            r(u, v) = (u, v, u² - v²)

        Gaußsche Krümmung: K < 0 (negativ definit).
        Mittlere Krümmung: H ≠ 0 im Allgemeinen.

        @return ParametricSurface-Objekt der Sattelfläche.
        @lastModified 2026-03-10
        """
        def r_saddle(u: float, v: float) -> np.ndarray:
            """Sattelflächenparametrisierung."""
            return np.array([u, v, u ** 2 - v ** 2])

        return ParametricSurface(r_saddle, (-1.0, 1.0), (-1.0, 1.0))

    def minimal_surface_check(
        self,
        surface: ParametricSurface,
        u_range: Tuple[float, float],
        v_range: Tuple[float, float],
        n_samples: int = 5,
        tol: float = 1e-3,
    ) -> dict:
        """
        Prüft, ob eine Fläche eine Minimalfläche ist (mittlere Krümmung H ≈ 0).

        Eine Minimalfläche hat H = 0 an jedem regulären Punkt.
        Sie sind stationäre Punkte des Flächeninhalts-Funktionals.

        @param surface: ParametricSurface-Objekt.
        @param u_range: Testbereich für u.
        @param v_range: Testbereich für v.
        @param n_samples: Anzahl Stichproben pro Dimension.
        @param tol: Toleranz für |H| < tol.
        @return Dictionary mit 'is_minimal', 'max_abs_H', 'mean_abs_H'.
        @lastModified 2026-03-10
        """
        u_vals = np.linspace(u_range[0] + 0.1, u_range[1] - 0.1, n_samples)
        v_vals = np.linspace(v_range[0] + 0.1, v_range[1] - 0.1, n_samples)
        H_values = []

        for u in u_vals:
            for v in v_vals:
                try:
                    H = surface.mean_curvature(u, v)
                    H_values.append(abs(H))
                except Exception:
                    pass  # Singuläre Punkte überspringen

        if not H_values:
            return {"is_minimal": False, "max_abs_H": float("nan"), "mean_abs_H": float("nan")}

        max_H = max(H_values)
        mean_H = sum(H_values) / len(H_values)
        return {
            "is_minimal": max_H < tol,
            "max_abs_H": float(max_H),
            "mean_abs_H": float(mean_H),
        }


# ---------------------------------------------------------------------------
# Klasse: CurveTheory2D
# ---------------------------------------------------------------------------

class CurveTheory2D:
    """
    Differentialgeometrie ebener Kurven im ℝ².

    Ebene Kurven haben eine besonders reiche Theorie:
    - Krümmung als vorzeichenbehaftete Größe
    - Evolute (Krümmungsmittelpunktskurve)
    - Evolvente (Abwicklungskurve)
    - Isoperimetrische Ungleichung: 4πA ≤ L²
    - Vier-Scheitel-Satz: Jede einfach geschlossene glatte Kurve hat ≥ 4 Scheitel
    - Gesamtkrümmung: ∫κ ds = 2πn (Umlaufzahl)

    Beispiel:
        >>> ct = CurveTheory2D()
        >>> import sympy as sp
        >>> x = sp.Symbol('x')
        >>> # Parabel y = x²: κ an x=0 maximal
        >>> abs(ct.curvature_from_cartesian(x**2, x, 0.0) - 2.0) < 0.01
        True
    """

    def curvature_from_cartesian(
        self,
        f_expr: sp.Expr,
        x_sym: sp.Symbol,
        x_val: float,
    ) -> float:
        """
        Berechnet die Krümmung einer kartesisch gegebenen Kurve y = f(x) an x = x_val.

        Formel: κ = |f''(x)| / (1 + f'(x)²)^{3/2}

        @param f_expr: SymPy-Ausdruck für f(x).
        @param x_sym: SymPy-Symbol für x.
        @param x_val: Auswertungsstelle x₀.
        @return Krümmung κ ≥ 0.
        @lastModified 2026-03-10
        """
        # Erste und zweite Ableitung symbolisch
        f_prime = sp.diff(f_expr, x_sym)
        f_double_prime = sp.diff(f_expr, x_sym, 2)

        # Numerische Auswertung
        fp = float(f_prime.subs(x_sym, x_val))
        fpp = float(f_double_prime.subs(x_sym, x_val))

        # Krümmungsformel
        denom = (1 + fp ** 2) ** 1.5
        if abs(denom) < 1e-15:
            return 0.0
        return abs(fpp) / denom

    def evolute(
        self,
        curve_components: List[sp.Expr],
        param: sp.Symbol,
    ) -> List[sp.Expr]:
        """
        Berechnet die Evolute einer ebenen parametrischen Kurve symbolisch.

        Die Evolute ist der geometrische Ort der Krümmungsmittelpunkte.
        Für eine Kurve r(t) = (x(t), y(t)) gilt:
            Evolute: e(t) = r(t) + (1/κ) · N(t)
            E_x = x - y'(x'² + y'²) / (x'y'' - y'x'')
            E_y = y + x'(x'² + y'²) / (x'y'' - y'x'')

        @param curve_components: [x(t), y(t)] als SymPy-Ausdrücke.
        @param param: SymPy-Symbol für den Parameter.
        @return [E_x(t), E_y(t)] als SymPy-Ausdrücke (symbolisch).
        @lastModified 2026-03-10
        """
        x, y = curve_components[0], curve_components[1]
        t = param

        # Ableitungen
        xp = sp.diff(x, t)
        yp = sp.diff(y, t)
        xpp = sp.diff(x, t, 2)
        ypp = sp.diff(y, t, 2)

        # Krümmungsform-Nenner
        kappa_denom = xp * ypp - yp * xpp
        speed_sq = xp ** 2 + yp ** 2

        # Evolute-Koordinaten
        E_x = x - yp * speed_sq / kappa_denom
        E_y = y + xp * speed_sq / kappa_denom

        # Vereinfachen
        E_x = sp.simplify(E_x)
        E_y = sp.simplify(E_y)

        return [E_x, E_y]

    def involute(
        self,
        curve_components: List[sp.Expr],
        param: sp.Symbol,
        c: float = 0.0,
    ) -> List[sp.Expr]:
        """
        Berechnet die Evolvente (Involute) einer parametrischen Kurve.

        Die Evolvente einer Kurve γ ist die Kurve, die ein aufgespulter Faden
        beschreibt. Formel:
            I(t) = r(t) - s(t) · T(t)

        wobei s(t) = ∫_{c}^{t} |r'(u)| du die Bogenlänge ist (symbolisch approximiert).

        @param curve_components: [x(t), y(t)] als SymPy-Ausdrücke.
        @param param: SymPy-Symbol für den Parameter.
        @param c: Integrationskonstante (Startwert der Bogenlänge).
        @return [I_x(t), I_y(t)] als SymPy-Ausdrücke.
        @lastModified 2026-03-10
        """
        x, y = curve_components[0], curve_components[1]
        t = param

        # Ableitungen
        xp = sp.diff(x, t)
        yp = sp.diff(y, t)

        # Geschwindigkeit |r'(t)|
        speed = sp.sqrt(xp ** 2 + yp ** 2)

        # Bogenlänge symbolisch (als Integral)
        # Für allgemeine Kurven: symbolisches Integral, kann komplex sein
        s = sp.Integral(speed, (t, c, t))

        # Tangenteneinheitsvektor T = (x', y') / |r'|
        T_x = xp / speed
        T_y = yp / speed

        # Evolvente: I = r - s * T
        I_x = x - s * T_x
        I_y = y - s * T_y

        return [I_x, I_y]

    def isoperimetric_inequality_check(
        self, perimeter: float, area: float
    ) -> dict:
        """
        Prüft die isoperimetrische Ungleichung: 4πA ≤ L².

        Die isoperimetrische Ungleichung besagt, dass unter allen geschlossenen
        Kurven gleicher Länge der Kreis den maximalen Flächeninhalt einschließt.

        Gleichheit gilt genau dann, wenn die Kurve ein Kreis ist.
        Isoperimetrisches Verhältnis: Q = 4πA/L² ∈ (0, 1].

        @param perimeter: Umfang L der Kurve.
        @param area: Flächeninhalt A der eingeschlossenen Region.
        @return Dictionary mit 'satisfied', 'ratio', 'is_circle'.
        @lastModified 2026-03-10
        """
        L = float(perimeter)
        A = float(area)

        # Isoperimetrisches Verhältnis Q = 4πA / L²
        if L < 1e-15:
            ratio = 0.0
        else:
            ratio = 4 * math.pi * A / (L ** 2)

        return {
            "satisfied": ratio <= 1.0 + 1e-10,  # 4πA ≤ L²
            "ratio": float(ratio),
            "is_circle": abs(ratio - 1.0) < 1e-6,
            "inequality": "4πA ≤ L²",
            "lhs": 4 * math.pi * A,
            "rhs": L ** 2,
        }

    def four_vertex_theorem_demo(self) -> dict:
        """
        Demonstration des Vier-Scheitel-Satzes.

        Vier-Scheitel-Satz (Mukhopadhyaya 1909):
            Jede einfach geschlossene, glatte, konvexe Kurve in der Ebene
            besitzt mindestens vier Scheitelpunkte, d.h. vier Punkte, an denen
            die Krümmung κ(t) ein lokales Extremum annimmt.

        Beweis-Idee: Die Krümmungsfunktion κ(t) ist eine periodische Funktion.
        Eine stetige, periodische Funktion mit Mittelwert ≠ const hat ≥ 2 Extrema.
        Durch den Vier-Scheitel-Satz folgen ≥ 4 Extrema.

        Für den Ellipsen-Beispiel:
            r(t) = (a cos t, b sin t) mit a ≠ b:
            κ(t) = ab / (a² sin²t + b² cos²t)^{3/2}
            4 Scheitel bei t = 0, π/2, π, 3π/2

        @return Dictionary mit Satzformulierung, Beispiel und Scheitelwerten.
        @lastModified 2026-03-10
        """
        # Ellipse a=2, b=1 als Beispiel
        a, b = 2.0, 1.0
        # endpoint=False: t=2π ≡ t=0 (geschlossene Kurve), letzten Punkt ausschließen
        t_vals = np.linspace(0, 2 * math.pi, 1000, endpoint=False)

        # Krümmung der Ellipse
        kappas = [
            a * b / (a ** 2 * math.sin(t) ** 2 + b ** 2 * math.cos(t) ** 2) ** 1.5
            for t in t_vals
        ]

        # Scheitel finden (lokale Extrema) – zyklische Randbehandlung für geschlossene Kurve
        kappas_arr = np.array(kappas)
        n_k = len(kappas_arr)
        vertices_idx = []
        for i in range(n_k):
            # Zyklische Nachbarn: letzter Punkt grenzt an ersten (Ringpuffer)
            prev_val = kappas_arr[(i - 1) % n_k]
            next_val = kappas_arr[(i + 1) % n_k]
            if (kappas_arr[i] > prev_val and kappas_arr[i] > next_val) or \
               (kappas_arr[i] < prev_val and kappas_arr[i] < next_val):
                vertices_idx.append(i)

        return {
            "theorem": (
                "Vier-Scheitel-Satz: Jede einfach geschlossene, glatte, konvexe "
                "Kurve in der Ebene hat mindestens 4 Scheitelpunkte (lokale Extrema von κ)."
            ),
            "example": "Ellipse r(t) = (2 cos t, sin t)",
            "a": a, "b": b,
            "num_vertices_found": len(vertices_idx),
            "satisfies_theorem": len(vertices_idx) >= 4,
            "kappa_max": float(np.max(kappas_arr)),
            "kappa_min": float(np.min(kappas_arr)),
        }

    def total_curvature_closed_curve(
        self,
        curve_fn: Callable[[float], List[float]],
        t_range: Tuple[float, float],
    ) -> float:
        """
        Berechnet die Gesamtkrümmung einer geschlossenen ebenen Kurve: ∫κ ds.

        Nach dem Umlaufzahlsatz gilt für einfach geschlossene Kurven:
            ∫₀ᴸ κ(s) ds = 2πn

        wobei n die Umlaufzahl (Windungszahl) der Kurve ist.
        Für konvexe Kurven: n = ±1, also Gesamtkrümmung = 2π.

        @param curve_fn: Funktion t → [x(t), y(t)] (parametrische 2D-Kurve).
        @param t_range: Parameterintervall (t_min, t_max).
        @return Gesamtkrümmung ∫κ ds (numerisch).
        @lastModified 2026-03-10
        """
        h = 1e-6  # Schrittweite für numerische Ableitungen

        def integrand(t: float) -> float:
            """κ(t) · |r'(t)| als Integrand."""
            # Numerische Ableitungen
            r0 = np.array(curve_fn(t), dtype=float)
            r_plus = np.array(curve_fn(t + h), dtype=float)
            r_minus = np.array(curve_fn(t - h), dtype=float)

            # Erste Ableitung
            rp = (r_plus - r_minus) / (2 * h)
            # Zweite Ableitung
            rpp = (r_plus - 2 * r0 + r_minus) / (h ** 2)

            xp, yp = rp[0], rp[1]
            xpp, ypp = rpp[0], rpp[1]

            # Vorzeichenbehaftete Krümmung × Geschwindigkeit = κ · |r'| = x'y'' - y'x''
            speed = math.sqrt(xp ** 2 + yp ** 2)
            if speed < 1e-15:
                return 0.0
            # Integrand: κ ds/dt = (x'y'' - y'x'') / |r'|² · |r'| · ... Nein:
            # ∫κ ds = ∫(x'y''-y'x'')/(x'²+y'²) · |r'| dt (Vorzeichenbehaftet)
            # Da wir Gesamtkrümmung (absolut) wollen:
            signed_kappa_speed = (xp * ypp - yp * xpp) / (xp ** 2 + yp ** 2)
            return signed_kappa_speed  # = κ · (positiv/negativ je nach Orientierung)

        # limit=200: mehr Unterteilungsschritte für glatte periodische Funktionen
        result, _ = quad(integrand, float(t_range[0]), float(t_range[1]), limit=200)
        return float(result)
