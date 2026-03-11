"""
Ricci-Fluss und Perelmansche Techniken

Dieses Modul implementiert numerische Methoden für den Ricci-Fluss auf
Riemannschen Mannigfaltigkeiten sowie die wesentlichen Konzepte aus
Grisha Perelmans Beweis der Poincaré-Vermutung.

Der Ricci-Fluss (R. Hamilton, 1982) ist eine partielle Differentialgleichung
für Metriken auf Mannigfaltigkeiten:

    ∂g_{ij}/∂t = -2 R_{ij}

Er deformiert die Geometrie so, dass positive Krümmung zusammengedrückt
und negative Krümmung ausgedehnt wird.

Referenzen:
- Hamilton, R.S. (1982): "Three-manifolds with positive Ricci curvature"
- Perelman, G. (2002): "The entropy formula for the Ricci flow and its geometric applications"
- Perelman, G. (2003): "Ricci flow with surgery on three-manifolds"

Autor: Michael Fuhrmann
Letzte Änderung: 2026-03-11
"""

import numpy as np
from typing import Callable, Optional


# ---------------------------------------------------------------------------
# 1. RicciFlowMetric — Numerische Krümmungsberechnungen
# ---------------------------------------------------------------------------

class RicciFlowMetric:
    """
    Riemannsche Metrik auf einer 2D- oder 3D-Mannigfaltigkeit.

    Die Metrik wird als Python-Funktion übergeben, die Koordinaten auf eine
    symmetrische (dim×dim)-Matrix abbildet.

    Alle geometrischen Größen werden numerisch über finite Differenzen
    berechnet.

    @param metric_func  Funktion coords → ndarray(dim,dim), gibt g_{ij}(p)
    @param dim          Dimension der Mannigfaltigkeit (2 oder 3)
    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def __init__(self, metric_func: Callable, dim: int = 2):
        """
        Initialisiert die Metrik.

        @param metric_func  Funktion, die Koordinaten auf Metrik-Matrix abbildet
        @param dim          Dimension (Standard: 2)
        """
        self.metric_func = metric_func
        self.dim = dim
        # Schrittweite für numerische Ableitungen
        self.h = 1e-4

    def _metric(self, coords: np.ndarray) -> np.ndarray:
        """
        Wertet die Metrik g_{ij} an den gegebenen Koordinaten aus.

        @param coords  Koordinatenvektor der Länge dim
        @return        Symmetrische (dim×dim)-Matrix
        """
        return np.array(self.metric_func(coords), dtype=float)

    def _inverse_metric(self, coords: np.ndarray) -> np.ndarray:
        """
        Berechnet die inverse Metrik g^{ij} = (g_{ij})^{-1}.

        @param coords  Koordinatenvektor
        @return        Inverse Metrik-Matrix
        """
        g = self._metric(coords)
        return np.linalg.inv(g)

    def _dg(self, coords: np.ndarray) -> np.ndarray:
        """
        Numerische partielle Ableitungen der Metrik: ∂g_{ij}/∂x^k

        Verwendet zentralen Differenzenquotienten der Ordnung O(h²).

        @param coords  Koordinatenvektor
        @return        Array der Form (dim, dim, dim), wobei [k, i, j] = ∂_k g_{ij}
        """
        n = self.dim
        h = self.h
        dg = np.zeros((n, n, n))  # [k, i, j] = d_k g_{ij}
        for k in range(n):
            coords_p = coords.copy()
            coords_m = coords.copy()
            coords_p[k] += h
            coords_m[k] -= h
            dg[k] = (self._metric(coords_p) - self._metric(coords_m)) / (2 * h)
        return dg

    def christoffel_symbols(self, coords: np.ndarray) -> np.ndarray:
        """
        Berechnet die Christoffel-Symbole zweiter Art Γ^k_{ij} numerisch.

        Die Formel lautet:
            Γ^k_{ij} = (1/2) g^{kl} (∂_i g_{jl} + ∂_j g_{il} - ∂_l g_{ij})

        @param coords  Koordinatenvektor
        @return        Array der Form (dim, dim, dim): Gamma[k, i, j]
        """
        n = self.dim
        g_inv = self._inverse_metric(coords)  # g^{kl}
        dg = self._dg(coords)                 # ∂_k g_{ij}

        Gamma = np.zeros((n, n, n))  # Gamma[k, i, j]
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    s = 0.0
                    for l in range(n):
                        # ∂_i g_{jl} + ∂_j g_{il} - ∂_l g_{ij}
                        s += g_inv[k, l] * (
                            dg[i, j, l] + dg[j, i, l] - dg[l, i, j]
                        )
                    Gamma[k, i, j] = 0.5 * s
        return Gamma

    def _d_christoffel(self, coords: np.ndarray) -> np.ndarray:
        """
        Numerische Ableitungen der Christoffel-Symbole: ∂_l Γ^k_{ij}

        @param coords  Koordinatenvektor
        @return        Array der Form (dim, dim, dim, dim): dGamma[l, k, i, j]
        """
        n = self.dim
        h = self.h
        dGamma = np.zeros((n, n, n, n))  # [l, k, i, j]
        for l in range(n):
            coords_p = coords.copy()
            coords_m = coords.copy()
            coords_p[l] += h
            coords_m[l] -= h
            Gp = self.christoffel_symbols(coords_p)
            Gm = self.christoffel_symbols(coords_m)
            dGamma[l] = (Gp - Gm) / (2 * h)
        return dGamma

    def ricci_tensor(self, coords: np.ndarray) -> np.ndarray:
        """
        Berechnet den Ricci-Tensor R_{ij} numerisch.

        Der Riemann-Krümmungstensor ist:
            R^l_{kij} = ∂_i Γ^l_{jk} - ∂_j Γ^l_{ik}
                       + Γ^l_{im} Γ^m_{jk} - Γ^l_{jm} Γ^m_{ik}

        Der Ricci-Tensor ergibt sich durch Kontraktion:
            R_{ij} = R^k_{ikj}

        @param coords  Koordinatenvektor
        @return        Symmetrische (dim×dim)-Matrix R_{ij}
        """
        n = self.dim
        Gamma = self.christoffel_symbols(coords)  # Γ^k_{ij}
        dGamma = self._d_christoffel(coords)       # ∂_l Γ^k_{ij}

        Ric = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                r = 0.0
                for k in range(n):
                    # ∂_k Γ^k_{ij} - ∂_j Γ^k_{ik}
                    r += dGamma[k, k, i, j] - dGamma[j, k, i, k]
                    # + Γ^k_{km} Γ^m_{ij} - Γ^k_{jm} Γ^m_{ik}
                    for m in range(n):
                        r += Gamma[k, k, m] * Gamma[m, i, j]
                        r -= Gamma[k, j, m] * Gamma[m, i, k]
                Ric[i, j] = r
        return Ric

    def ricci_scalar(self, coords: np.ndarray) -> float:
        """
        Berechnet den Ricci-Skalare R = g^{ij} R_{ij}.

        @param coords  Koordinatenvektor
        @return        Skalare Krümmung (float)
        """
        g_inv = self._inverse_metric(coords)
        Ric = self.ricci_tensor(coords)
        return float(np.einsum('ij,ij->', g_inv, Ric))

    def sectional_curvature(self, coords: np.ndarray,
                             v1: np.ndarray, v2: np.ndarray) -> float:
        """
        Berechnet die Schnittkriimmung K(v1, v2) für zwei Tangentialvektoren.

        Formel:
            K(v1, v2) = R(v1, v2, v1, v2) / (g(v1,v1)g(v2,v2) - g(v1,v2)²)

        wobei R_{ijkl} = g_{im} R^m_{jkl} der ko-variante Riemann-Tensor ist.

        Hier wird die vereinfachte 2D-Formel K = R / (det g) verwendet.

        @param coords  Koordinatenvektor
        @param v1      Erster Tangentialvektor
        @param v2      Zweiter Tangentialvektor
        @return        Schnittkriimmung K
        """
        g = self._metric(coords)
        # Gram-Gramian für die Normalisierung
        g11 = float(v1 @ g @ v1)
        g22 = float(v2 @ g @ v2)
        g12 = float(v1 @ g @ v2)
        denom = g11 * g22 - g12 ** 2
        if abs(denom) < 1e-15:
            return 0.0

        # Verwende R_{1212} aus dem Riemann-Tensor
        n = self.dim
        Gamma = self.christoffel_symbols(coords)
        dGamma = self._d_christoffel(coords)

        # R^l_{kij} berechnen
        R_up = np.zeros((n, n, n, n))  # R^l_{k i j}
        for l in range(n):
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        r = dGamma[i, l, j, k] - dGamma[j, l, i, k]
                        for m in range(n):
                            r += Gamma[l, i, m] * Gamma[m, j, k]
                            r -= Gamma[l, j, m] * Gamma[m, i, k]
                        R_up[l, k, i, j] = r

        # Senke ersten Index: R_{abij} = g_{al} R^l_{bij}
        g_mat = self._metric(coords)
        R_down = np.einsum('al,lbij->abij', g_mat, R_up)

        # Schnittkriimmung: K = R(v1,v2,v1,v2) / Gram-Det
        num = np.einsum('a,b,c,d,abcd->', v1, v2, v1, v2, R_down)
        return float(num / denom)


# ---------------------------------------------------------------------------
# 2. RicciFlow — Ricci-Fluss-Evolution
# ---------------------------------------------------------------------------

class RicciFlow:
    """
    Numerische Integration des Ricci-Flusses ∂g/∂t = -2 Ric(g).

    Der Ricci-Fluss deformiert die Riemannsche Metrik entlang ihrer
    Ricci-Krümmung. Bei positiver Krümmung (z.B. Sphäre) wird die Metrik
    homogenisiert; bei negativer Krümmung expandiert sie.

    Implementierung: Runge-Kutta 4. Ordnung für die Metrik-ODE.

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def __init__(self, initial_metric: Callable, dim: int = 2):
        """
        Initialisiert den Ricci-Fluss.

        @param initial_metric  Anfangsmetrik g_0 als Funktion coords → Matrix
        @param dim             Dimension der Mannigfaltigkeit
        """
        self.initial_metric = initial_metric
        self.dim = dim
        # Aktuelle Metrik (als Skalierungsfaktor der Einheitsmetrik)
        self._scale = 1.0
        # Interne Metrik-Funktion (wird bei evolve aktualisiert)
        self._current_metric = initial_metric

    def _flow_derivative(self, coords: np.ndarray, metric_func: Callable) -> np.ndarray:
        """
        Berechnet -2 Ric(g) an den gegebenen Koordinaten.

        @param coords       Koordinatenvektor
        @param metric_func  Aktuelle Metrik-Funktion
        @return             (dim×dim)-Matrix: rechte Seite der Fluss-ODE
        """
        rfm = RicciFlowMetric(metric_func, self.dim)
        Ric = rfm.ricci_tensor(coords)
        return -2.0 * Ric

    def evolve_step(self, coords: np.ndarray, dt: float) -> np.ndarray:
        """
        Führt einen RK4-Schritt des Ricci-Flusses durch.

        Gibt die aktualisierte Metrik g(t + dt) zurück.

        @param coords  Koordinatenvektor (Auswertungspunkt)
        @param dt      Zeitschritt
        @return        Neue Metrik-Matrix g_{ij}(t + dt)
        """
        g0 = np.array(self._current_metric(coords), dtype=float)

        # RK4-Koeffizienten für die Metrik-Matrix
        k1 = self._flow_derivative(coords, self._current_metric)

        def mk2_func(c):
            return np.array(self._current_metric(c), dtype=float) + 0.5 * dt * k1

        k2 = -2.0 * RicciFlowMetric(mk2_func, self.dim).ricci_tensor(coords)

        def mk3_func(c):
            return np.array(self._current_metric(c), dtype=float) + 0.5 * dt * k2

        k3 = -2.0 * RicciFlowMetric(mk3_func, self.dim).ricci_tensor(coords)

        def mk4_func(c):
            return np.array(self._current_metric(c), dtype=float) + dt * k3

        k4 = -2.0 * RicciFlowMetric(mk4_func, self.dim).ricci_tensor(coords)

        g_new = g0 + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        return g_new

    def evolve(self, coords: np.ndarray, T: float,
               steps: int = 100) -> list:
        """
        Integriert den Ricci-Fluss von t=0 bis t=T.

        @param coords  Auswertungskoordinaten
        @param T       Endzeit
        @param steps   Anzahl der Zeitschritte
        @return        Liste von (t, g_matrix) Paaren
        """
        dt = T / steps
        history = []
        g_current = np.array(self.initial_metric(coords), dtype=float)
        t = 0.0

        history.append((t, g_current.copy()))

        for _ in range(steps):
            # Einfacher Euler-Schritt für Geschwindigkeit (RK4 zu teuer hier)
            rfm = RicciFlowMetric(
                lambda c, g=g_current: g, self.dim
            )
            Ric = rfm.ricci_tensor(coords)
            g_current = g_current - 2.0 * dt * Ric
            t += dt
            history.append((t, g_current.copy()))

        return history

    def detect_singularity(self, coords: np.ndarray, t: float,
                           threshold: float = 1e6) -> bool:
        """
        Erkennt Singularitäten durch Überwachung der Ricci-Krümmung.

        Eine Singularität liegt vor, wenn |R_{ij}| > threshold (Krümmungsexplosion).

        @param coords     Auswertungskoordinaten
        @param t          Aktuelle Zeit
        @param threshold  Schwellwert für Krümmungsexplosion
        @return           True wenn Singularität erkannt
        """
        rfm = RicciFlowMetric(self._current_metric, self.dim)
        Ric = rfm.ricci_tensor(coords)
        max_curv = float(np.max(np.abs(Ric)))
        return max_curv > threshold


# ---------------------------------------------------------------------------
# 3. NormalizedRicciFlow — Volumen-erhaltender Ricci-Fluss
# ---------------------------------------------------------------------------

class NormalizedRicciFlow:
    """
    Normalisierter Ricci-Fluss mit Volumenerhaltung.

    Gleichung:
        ∂g_{ij}/∂t = -2 R_{ij} + (2r/n) g_{ij}

    wobei r = ∫R dV / ∫dV die mittlere skalare Krümmung ist.

    Dieser Fluss konvergiert auf Flächen (dim=2) gegen eine Metrik
    konstanter Gauss-Krümmung (Satz von Hamilton, 1988).

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def __init__(self, initial_metric: Callable, dim: int = 2):
        """
        Initialisiert den normalisierten Ricci-Fluss.

        @param initial_metric  Anfangsmetrik
        @param dim             Dimension
        """
        self.initial_metric = initial_metric
        self.dim = dim
        self._current_metric = initial_metric

    def mean_curvature(self, coords: np.ndarray) -> float:
        """
        Schätzt die mittlere skalare Krümmung r = ∫R dV / ∫dV.

        Vereinfachung: Auswertung am gegebenen Punkt (keine Integration).

        @param coords  Auswertungskoordinaten
        @return        Mittlere skalare Krümmung r
        """
        rfm = RicciFlowMetric(self._current_metric, self.dim)
        return rfm.ricci_scalar(coords)

    def evolve_step(self, coords: np.ndarray, dt: float) -> np.ndarray:
        """
        Führt einen normalisierten Ricci-Fluss-Schritt durch.

        Formel: g_{ij}(t+dt) = g_{ij}(t) + dt·(-2 R_{ij} + (2r/n) g_{ij})

        @param coords  Auswertungskoordinaten
        @param dt      Zeitschritt
        @return        Neue Metrik-Matrix
        """
        rfm = RicciFlowMetric(self._current_metric, self.dim)
        g = rfm._metric(coords)
        Ric = rfm.ricci_tensor(coords)
        r = rfm.ricci_scalar(coords)
        n = self.dim

        # Normalisierungsterm
        dg = -2.0 * Ric + (2.0 * r / n) * g
        return g + dt * dg

    def long_time_behavior(self, coords: np.ndarray, T: float,
                           steps: int = 200) -> dict:
        """
        Untersucht das Langzeitverhalten des normalisierten Ricci-Flusses.

        @param coords  Auswertungskoordinaten
        @param T       Endzeit
        @param steps   Anzahl der Schritte
        @return        Dict mit 'curvature_history', 'final_metric', 'converged'
        """
        dt = T / steps
        g_current = np.array(self.initial_metric(coords), dtype=float)
        curvature_history = []

        for _ in range(steps):
            rfm = RicciFlowMetric(lambda c, g=g_current: g, self.dim)
            R = rfm.ricci_scalar(coords)
            Ric = rfm.ricci_tensor(coords)
            r = R
            n = self.dim
            dg = -2.0 * Ric + (2.0 * r / n) * g_current
            g_current = g_current + dt * dg
            curvature_history.append(float(R))

        # Konvergenz: Krümmungsvariation in letzten 20% der Schritte
        tail = curvature_history[int(0.8 * steps):]
        converged = (max(tail) - min(tail)) < 0.1 * abs(np.mean(tail) + 1e-15)

        return {
            'curvature_history': curvature_history,
            'final_metric': g_current,
            'converged': converged
        }


# ---------------------------------------------------------------------------
# 4. SphericalMetric — Sphärische Geometrie
# ---------------------------------------------------------------------------

class SphericalMetric:
    """
    Runde Metrik auf der n-Sphäre S^n.

    In lokalen Karten (stereographische Projektion) hat S^n die Metrik
    g_{ij} = (2/(1+|x|²))² δ_{ij}.

    In der einfachsten Form: Standard-Euklid (flache Karte).

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def standard_metric(self, dim: int = 2) -> Callable:
        """
        Gibt die Standard-Einheitsmetrik als Funktion zurück.

        Im euklidischen Koordinatensystem: g_{ij} = δ_{ij}.

        @param dim  Dimension
        @return     Metrik-Funktion coords → ndarray(dim, dim)
        """
        def metric_func(coords):
            return np.eye(dim)
        return metric_func

    def stereographic_metric(self, dim: int = 2) -> Callable:
        """
        Stereographische Projektion: g_{ij} = (2/(1+|x|²))² δ_{ij}

        @param dim  Dimension der Sphäre
        @return     Metrik-Funktion
        """
        def metric_func(coords):
            r2 = float(np.dot(coords[:dim], coords[:dim]))
            factor = (2.0 / (1.0 + r2)) ** 2
            return factor * np.eye(dim)
        return metric_func

    def sectional_curvature_sphere(self, dim: int) -> float:
        """
        Gibt die konstante Schnittkriimmung der Einheitssphäre zurück.

        Für S^n mit Radius 1: K ≡ +1.

        @param dim  Dimension (nicht benötigt für Einheitssphäre)
        @return     +1.0
        """
        return 1.0

    def ricci_tensor_sphere(self, dim: int) -> np.ndarray:
        """
        Berechnet den Ricci-Tensor der Einheitssphäre.

        Für S^n gilt: R_{ij} = (n-1) g_{ij} = (dim-1) δ_{ij}.

        @param dim  Dimension der Sphäre
        @return     Diagonalmatrix (dim×dim)
        """
        return (dim - 1) * np.eye(dim)


# ---------------------------------------------------------------------------
# 5. HyperbolicMetric — Hyperbolische Geometrie
# ---------------------------------------------------------------------------

class HyperbolicMetric:
    """
    Hyperbolische Metrik auf dem Poincaré-Modell.

    Das Poincaré-Scheibenmodell hat die Metrik:
        g_{ij}(x) = 4 / (1 - |x|²)² · δ_{ij}

    Dies ist eine Riemannsche Metrik auf dem offenen Einheitskreis.
    Die skalare Krümmung ist konstant -1.

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def poincare_disk_metric(self, coords: np.ndarray) -> np.ndarray:
        """
        Berechnet die Poincaré-Scheibenmetrik g_{ij}(x).

        g_{ij} = 4 / (1 - |x|²)² · δ_{ij}

        @param coords  2D-Koordinatenvektor (muss |coords| < 1 gelten)
        @return        (2×2)-Metrik-Matrix
        """
        x = coords[:2]
        r2 = float(np.dot(x, x))
        # Sicherstellung: innerhalb der Einheitsscheibe
        r2 = min(r2, 0.99)
        factor = 4.0 / (1.0 - r2) ** 2
        return factor * np.eye(2)

    def poincare_disk_metric_func(self) -> Callable:
        """
        Gibt die Poincaré-Metrik als Callable zurück.

        @return  Funktion coords → (2×2)-Matrix
        """
        return self.poincare_disk_metric

    def sectional_curvature_hyperbolic(self) -> float:
        """
        Gibt die konstante Schnittkriimmung des hyperbolischen Raums zurück.

        Für das Poincaré-Modell: K ≡ -1.

        @return  -1.0
        """
        return -1.0

    def hyperbolic_distance(self, z1: complex, z2: complex) -> float:
        """
        Berechnet den hyperbolischen Abstand zweier Punkte.

        Formel (Poincaré-Scheibe):
            d(z1, z2) = 2 arctanh(|z1 - z2| / |1 - conj(z1)·z2|)

        @param z1  Erster Punkt (komplexe Zahl, |z| < 1)
        @param z2  Zweiter Punkt (komplexe Zahl, |z| < 1)
        @return    Hyperbolischer Abstand
        """
        num = abs(z1 - z2)
        denom = abs(1 - z1.conjugate() * z2)
        if denom < 1e-15:
            return float('inf')
        ratio = num / denom
        ratio = min(ratio, 1.0 - 1e-15)
        return 2.0 * np.arctanh(ratio)


# ---------------------------------------------------------------------------
# 6. RicciFlowSurgery — Perelmansche Chirurgie
# ---------------------------------------------------------------------------

class RicciFlowSurgery:
    """
    Konzeptuelle Implementierung der Perelmannschen Chirurgie-Methode.

    Die Chirurgie (Surgery) ist eine topologische Operation, die
    Singularitäten des Ricci-Flusses behandelt, indem sie Engstellen
    (ε-Hälse) abschneidet und durch Kappen ersetzt.

    Dies ist konzeptuell — eine vollständige Implementierung würde
    3D-Triangulierungen erfordern.

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def neck_detection(self, metric_data: list,
                       threshold: float = 0.5) -> list:
        """
        Erkennt ε-Hals-Strukturen (ε-neck) in den Metrikdaten.

        Ein ε-Hals ist eine Region, die einem Zylinder S²×[-L, L] ähnelt,
        wobei die Skalierungsfunktion klein ist.

        Heuristik: Regionen mit kleinem Volumen-Element werden als Hälse markiert.

        @param metric_data  Liste von Metrik-Matrizen entlang einer Kurve
        @param threshold    Schwellwert für das Volumen-Element
        @return             Liste der Indizes mit Hals-Strukturen
        """
        necks = []
        for i, g in enumerate(metric_data):
            g_mat = np.array(g)
            vol = compute_volume_element(g_mat)
            if vol < threshold:
                necks.append(i)
        return necks

    def surgery_step(self, metric_data: list,
                     singularity_time: float) -> dict:
        """
        Führt einen konzeptuellen Chirurgie-Schritt durch.

        Der Algorithmus:
        1. Erkenne Hals-Strukturen (ε-necks)
        2. Schneide entlang der Hälse ab
        3. Verschließe die Enden mit Kappen (3-Bällen)
        4. Setze den Ricci-Fluss fort

        @param metric_data        Liste von Metrik-Matrizen
        @param singularity_time   Zeit der Singularität
        @return                   Dict mit Beschreibung der Chirurgie
        """
        necks = self.neck_detection(metric_data)
        n_components = len(necks) + 1 if necks else 1

        return {
            'singularity_time': singularity_time,
            'neck_count': len(necks),
            'neck_positions': necks,
            'resulting_components': n_components,
            'surgery_performed': len(necks) > 0,
            'description': self.describe_surgery()
        }

    def describe_surgery(self) -> str:
        """
        Beschreibt die Perelmannsche Chirurgie-Methode.

        @return  Textbeschreibung der Methode
        """
        return (
            "Perelmannsche Ricci-Fluss-Chirurgie:\n"
            "1. Erkenne ε-Hals-Strukturen (zylinderähnliche Regionen) "
            "kurz vor einer Singularität.\n"
            "2. Schneide entlang der Hälse ab und verschließe die Enden "
            "mit 3-Bällen (Kappen).\n"
            "3. Setze den normalisierten Ricci-Fluss auf den Komponenten fort.\n"
            "4. Wiederhole bis keine Singularitäten mehr auftreten oder "
            "die Mannigfaltigkeit S³ ist.\n"
            "Schlüsselergebnis: Jede geschlossene, einfach zusammenhängende "
            "3-Mannigfaltigkeit ist homöomorph zu S³ (Poincaré-Vermutung)."
        )


# ---------------------------------------------------------------------------
# 7. PoincareConjecture — Formale Klasse
# ---------------------------------------------------------------------------

class PoincareConjecture:
    """
    Formale Klasse zur Poincaré-Vermutung und Perelmans Beweis.

    Die Poincaré-Vermutung (1904, bewiesen 2003 von Perelman):
    Jede geschlossene, einfach zusammenhängende 3-Mannigfaltigkeit ist
    homöomorph zur 3-Sphäre S³.

    Dies ist eines der Millennium-Probleme (Clay Mathematics Institute).

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def statement(self) -> str:
        """
        Gibt die formale Aussage der Poincaré-Vermutung zurück.

        @return  Aussage als String
        """
        return (
            "Poincaré-Vermutung (Grisha Perelman, 2003):\n\n"
            "Sei M eine geschlossene (kompakt, ohne Rand), "
            "einfach zusammenhängende 3-dimensionale Mannigfaltigkeit.\n"
            "Dann ist M homöomorph zur 3-Sphäre S³.\n\n"
            "Äquivalente Formulierung: "
            "Jede geschlossene 3-Mannigfaltigkeit mit trivialem "
            "Fundamentalgruppe π₁(M) = {e} ist homöomorph zu S³.\n\n"
            "Beweismethode: Ricci-Fluss mit Chirurgie (Hamilton-Perelman)."
        )

    def perelman_sketch(self) -> dict:
        """
        Gibt eine Beweisskizze von Perelmans Beweis zurück.

        @return  Dict mit Schlüsselschritten
        """
        return {
            'schritt_1': (
                "W-Funktional: Einführung von Perelmanns W-Entropie "
                "W(g, f, τ) = ∫[τ(R + |∇f|²) + f - n] e^{-f}/(4πτ)^{n/2} dV. "
                "Monoton wachsend entlang des Ricci-Flusses."
            ),
            'schritt_2': (
                "κ-Nicht-Kollaps-Theorem: Der Ricci-Fluss kollabiert nicht lokal "
                "(keine zu kleinen geodätischen Bälle), was Singularitäten kontrolliert."
            ),
            'schritt_3': (
                "Klassifikation der Singularitäten: Alle Singularitäten sind "
                "entweder runde Sphären oder Zylinder (κ-Lösungen)."
            ),
            'schritt_4': (
                "Chirurgie: Entlang ε-Hälsen wird die Mannigfaltigkeit aufgeschnitten "
                "und mit Kappen verschlossen. Das Geschlecht sinkt bei jeder Chirurgie."
            ),
            'schritt_5': (
                "Endlichkeit der Chirurgien: Es gibt nur endlich viele Chirurgien "
                "in endlicher Zeit (geometrisierungstheoretische Schranken)."
            ),
            'schritt_6': (
                "Geometrisierung: Nach dem Fluss zerfällt M in Stücke mit einer "
                "der 8 Thurston-Geometrien. Für einfach zusammenhängende M: S³."
            ),
            'konklusion': (
                "Da M einfach zusammenhängend ist und nach Chirurgie in S³-Stücke "
                "zerfällt: M ≅ S³. □"
            )
        }

    def verify_s3_topology(self, betti_numbers: list) -> bool:
        """
        Überprüft ob die Betti-Zahlen zur 3-Sphäre passen.

        Für S³ gilt: β = (β₀, β₁, β₂, β₃) = (1, 0, 0, 1).
        Euler-Charakteristik: χ(S³) = 1 - 0 + 0 - 1 = 0.

        @param betti_numbers  Liste [β₀, β₁, β₂, β₃]
        @return               True wenn Betti-Zahlen zu S³ passen
        """
        if len(betti_numbers) < 4:
            return False
        return (
            betti_numbers[0] == 1 and
            betti_numbers[1] == 0 and
            betti_numbers[2] == 0 and
            betti_numbers[3] == 1
        )

    def hamilton_program(self) -> dict:
        """
        Beschreibt Hamiltons Programm zur Poincaré-Vermutung.

        @return  Dict mit Hamiltons Ansatz
        """
        return {
            'idee': (
                "Richard Hamilton (1982) schlug vor, den Ricci-Fluss zu verwenden, "
                "um 3-Mannigfaltigkeiten zu klassifizieren."
            ),
            'ergebnis_1982': (
                "3-Mannigfaltigkeiten mit positiver Ricci-Krümmung sind "
                "diffeomorph zu S³ (Hamiltons erster Satz)."
            ),
            'hauptproblem': (
                "Singularitäten des Ricci-Flusses verhinderten die Durchführung "
                "des Programms für allgemeine 3-Mannigfaltigkeiten."
            ),
            'perelman_lösung': (
                "Perelman (2002-2003) löste das Singularitätsproblem durch "
                "W-Funktional, κ-Nicht-Kollaps und Chirurgie."
            ),
            'geometrisierung': (
                "Das vollständige Programm beweist auch Thurstons "
                "Geometrisierungsvermutung als Spezialfall."
            )
        }

    def perelman_contributions(self) -> list:
        """
        Listet Perelmans wichtigste Beiträge auf.

        @return  Liste der Schlüsselbeiträge
        """
        return [
            {
                'name': 'W-Funktional',
                'beschreibung': (
                    "W(g, f, τ) = ∫[τ(R + |∇f|²) + f - n] · "
                    "e^{-f}/(4πτ)^{n/2} dV. "
                    "Monoton wachsend — ersetzt fehlende Kompaktheit."
                )
            },
            {
                'name': 'L-Funktional',
                'beschreibung': (
                    "L(γ) = ∫₀^τ √t (R(γ(t)) + |γ'(t)|²) dt. "
                    "Definiert L-Geodäten und reduziertes Volumen."
                )
            },
            {
                'name': 'κ-Nicht-Kollaps',
                'beschreibung': (
                    "Wenn |Rm| ≤ 1/r² auf B(x,r), dann Vol(B(x,r)) ≥ κ·rⁿ. "
                    "Verhindert geometrisches Kollabieren."
                )
            },
            {
                'name': 'κ-Lösungen',
                'beschreibung': (
                    "Klassifikation der Blow-up-Grenzwerte: "
                    "Sphären, Zylinder oder Bryant-Solitone."
                )
            },
            {
                'name': 'Chirurgie-Algorithmus',
                'beschreibung': (
                    "Präzise Definition der Chirurgie entlang ε-Hälsen "
                    "mit Kontrolle aller geometrischen Parameter."
                )
            }
        ]


# ---------------------------------------------------------------------------
# 8. EntropyFunctional — Perelmansche Entropie-Funktionale
# ---------------------------------------------------------------------------

class EntropyFunctional:
    """
    Perelmansche Entropie-Funktionale für den Ricci-Fluss.

    Die Monotonizität der Funktionale ist der Schlüssel zu Perelmans Beweis:
    - F(g, f): Monoton wachsend entlang (gekoppeltem) Ricci-Fluss
    - W(g, f, τ): Monoton wachsend → κ-Nicht-Kollaps

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def w_functional(self, g_func: Callable, f_func: Callable,
                     tau: float, coords_list: list,
                     weights: Optional[np.ndarray] = None) -> float:
        """
        Berechnet das W-Funktional von Perelman numerisch.

        W(g, f, τ) = ∫[τ(R + |∇f|²) + f - n] · e^{-f}/(4πτ)^{n/2} dV

        Approximation: Summe über gegebene Koordinatenpunkte.

        @param g_func       Metrik-Funktion
        @param f_func       Skalare Funktion f (potentialähnlich)
        @param tau          Zeitparameter τ > 0
        @param coords_list  Liste von Koordinatenvektoren
        @param weights      Gewichte für numerische Integration
        @return             Numerischer Wert des W-Funktionals
        """
        dim = len(coords_list[0])
        n = dim
        normalization = (4 * np.pi * tau) ** (-n / 2.0)

        if weights is None:
            weights = np.ones(len(coords_list)) / len(coords_list)

        W = 0.0
        for coords, w in zip(coords_list, weights):
            coords = np.array(coords, dtype=float)
            rfm = RicciFlowMetric(g_func, dim)

            # Skalare Krümmung R
            R = rfm.ricci_scalar(coords)

            # f und |∇f|² numerisch
            f_val = float(f_func(coords))

            # Gradient von f numerisch
            h = 1e-5
            grad_f = np.zeros(dim)
            for k in range(dim):
                cp = coords.copy(); cp[k] += h
                cm = coords.copy(); cm[k] -= h
                grad_f[k] = (f_func(cp) - f_func(cm)) / (2 * h)

            # |∇f|² = g^{ij} ∂_i f ∂_j f
            g_inv = rfm._inverse_metric(coords)
            grad_f_sq = float(grad_f @ g_inv @ grad_f)

            # Volumen-Element
            g_mat = rfm._metric(coords)
            dV = compute_volume_element(g_mat)

            # Integrand
            integrand = (tau * (R + grad_f_sq) + f_val - n)
            integrand *= np.exp(-f_val) * normalization * dV

            W += w * integrand

        return float(W)

    def f_functional(self, g_func: Callable, f_func: Callable,
                     coords_list: list,
                     weights: Optional[np.ndarray] = None) -> float:
        """
        Berechnet das F-Funktional von Perelman numerisch.

        F(g, f) = ∫(R + |∇f|²) · e^{-f} dV

        @param g_func       Metrik-Funktion
        @param f_func       Skalare Funktion f
        @param coords_list  Liste von Koordinatenvektoren
        @param weights      Gewichte für numerische Integration
        @return             Numerischer Wert des F-Funktionals
        """
        dim = len(coords_list[0])

        if weights is None:
            weights = np.ones(len(coords_list)) / len(coords_list)

        F = 0.0
        for coords, w in zip(coords_list, weights):
            coords = np.array(coords, dtype=float)
            rfm = RicciFlowMetric(g_func, dim)

            R = rfm.ricci_scalar(coords)
            f_val = float(f_func(coords))

            # Gradient numerisch
            h = 1e-5
            grad_f = np.zeros(dim)
            for k in range(dim):
                cp = coords.copy(); cp[k] += h
                cm = coords.copy(); cm[k] -= h
                grad_f[k] = (f_func(cp) - f_func(cm)) / (2 * h)

            g_inv = rfm._inverse_metric(coords)
            grad_f_sq = float(grad_f @ g_inv @ grad_f)

            g_mat = rfm._metric(coords)
            dV = compute_volume_element(g_mat)

            integrand = (R + grad_f_sq) * np.exp(-f_val) * dV
            F += w * integrand

        return float(F)

    def monotonicity_statement(self) -> str:
        """
        Gibt die Monotonieeigenschaften der Funktionale aus.

        @return  Beschreibung der Monotonizität
        """
        return (
            "Perelmansche Monotonieeigenschaften:\n\n"
            "1. F-Funktional: d/dt F(g(t), f(t)) = 2∫|Ric + ∇²f|² e^{-f} dV ≥ 0\n"
            "   (monoton wachsend entlang des Ricci-Flusses mit ∂_t f = -R - Δf)\n\n"
            "2. W-Funktional: d/dt W(g(t), f(t), τ(t)) ≥ 0\n"
            "   für τ(t) = T - t (rückwärts in der Zeit).\n\n"
            "3. Konsequenz: λ(g) = inf_f F(g,f) ist monoton wachsend.\n"
            "   Dies verhindert lokales Kollabieren (κ-Nicht-Kollaps-Theorem).\n\n"
            "4. Minimales W: μ(g, τ) = inf_{f} W(g, f, τ) ist monoton → "
            "Sobolev-Konstanten bleiben kontrolliert."
        )


# ---------------------------------------------------------------------------
# 9. ThreeDimensionalTopology — 3D Topologie
# ---------------------------------------------------------------------------

class ThreeDimensionalTopology:
    """
    Hilfsmethoden für dreidimensionale Topologie.

    Im Kontext der Poincaré-Vermutung und Thurstons Geometrisierungsprogramm.

    @author Michael Fuhrmann
    @since  2026-03-11
    """

    def euler_characteristic(self, betti_numbers: list) -> int:
        """
        Berechnet die Euler-Charakteristik χ = Σ(-1)^i β_i.

        @param betti_numbers  Liste der Betti-Zahlen [β₀, β₁, β₂, ...]
        @return               Euler-Charakteristik
        """
        chi = 0
        for i, beta in enumerate(betti_numbers):
            chi += ((-1) ** i) * beta
        return chi

    def is_simply_connected(self, fundamental_group_trivial: bool) -> bool:
        """
        Überprüft ob eine Mannigfaltigkeit einfach zusammenhängend ist.

        Eine Mannigfaltigkeit M ist einfach zusammenhängend gdw.
        π₁(M) = {e} (triviale Fundamentalgruppe).

        @param fundamental_group_trivial  True wenn π₁(M) trivial ist
        @return                          True wenn M einfach zusammenhängend
        """
        return fundamental_group_trivial

    def geometrization_conjecture_cases(self) -> dict:
        """
        Listet die 8 Thurston-Geometrien für die Geometrisierungsvermutung.

        Thurston (1982): Jede geschlossene 3-Mannigfaltigkeit zerfällt in
        Stücke, die je eine der 8 Modellgeometrien tragen.

        @return  Dict mit den 8 Geometrien und ihren Eigenschaften
        """
        return {
            'S³': {
                'name': '3-Sphäre',
                'krümmung': 'positiv konstant',
                'beispiel': 'S³, Linsenräume, Prismräume',
                'isometriegruppe': 'O(4)'
            },
            'E³': {
                'name': 'Euklidischer Raum',
                'krümmung': '0 (flach)',
                'beispiel': '3-Torus T³, 10 flache 3-Mannigfaltigkeiten',
                'isometriegruppe': 'ISO(3)'
            },
            'H³': {
                'name': 'Hyperbolischer Raum',
                'krümmung': 'negativ konstant -1',
                'beispiel': 'Hyperbolische 3-Mannigfaltigkeiten',
                'isometriegruppe': 'Isom(H³)'
            },
            'S²×ℝ': {
                'name': 'Produkt S² × ℝ',
                'krümmung': 'gemischt',
                'beispiel': 'S²×S¹',
                'isometriegruppe': 'Isom(S²) × Isom(ℝ)'
            },
            'H²×ℝ': {
                'name': 'Produkt H² × ℝ',
                'krümmung': 'gemischt',
                'beispiel': 'Σ_g × S¹ (Fläche Geschlecht g≥2)',
                'isometriegruppe': 'Isom(H²) × Isom(ℝ)'
            },
            'SL₂(ℝ)': {
                'name': 'Universelle Überlagerung von PSL₂(ℝ)',
                'krümmung': 'gemischt',
                'beispiel': 'Einheitstangentenbündel über H²',
                'isometriegruppe': 'SL₂(ℝ)'
            },
            'Nil': {
                'name': 'Nil-Geometrie (Heisenberg-Gruppe)',
                'krümmung': 'gemischt',
                'beispiel': 'Nil-Mannigfaltigkeiten',
                'isometriegruppe': 'Nil ⋊ U(1)'
            },
            'Sol': {
                'name': 'Sol-Geometrie',
                'krümmung': 'gemischt',
                'beispiel': 'Sol-Mannigfaltigkeiten',
                'isometriegruppe': 'Sol ⋊ Z/2'
            }
        }

    def classify_3manifold(self, properties: dict) -> str:
        """
        Klassifiziert eine 3-Mannigfaltigkeit anhand ihrer Eigenschaften.

        @param properties  Dict mit Eigenschaften (z.B. 'simply_connected',
                           'positive_curvature', 'betti_numbers', usw.)
        @return            Klassifikationsergebnis als String
        """
        simply_connected = properties.get('simply_connected', False)
        positive_curvature = properties.get('positive_curvature', False)
        betti = properties.get('betti_numbers', [])
        flat = properties.get('flat', False)
        hyperbolic = properties.get('hyperbolic', False)

        if simply_connected and positive_curvature:
            return "S³ (3-Sphäre) — Poincaré-Vermutung"

        if simply_connected and len(betti) >= 4 and betti == [1, 0, 0, 1]:
            return "S³ (3-Sphäre) — via Betti-Zahlen"

        if flat:
            return "Flache 3-Mannigfaltigkeit — eine der 10 Typen (E³-Geometrie)"

        if hyperbolic:
            return "Hyperbolische 3-Mannigfaltigkeit — H³-Geometrie"

        if not simply_connected and len(betti) >= 2 and betti[1] > 0:
            return "Nicht-einfach-zusammenhängende 3-Mannigfaltigkeit — komplex"

        return "Unbestimmt — weitere Eigenschaften benötigt"


# ---------------------------------------------------------------------------
# 10. Hilfsfunktionen auf Modul-Ebene
# ---------------------------------------------------------------------------

def ricci_flow_on_surface(metric_func: Callable, T: float,
                          steps: int = 50) -> list:
    """
    Führt den Ricci-Fluss auf einer 2D-Fläche durch.

    Dies ist eine Vereinfachung: Der Fluss wird an einem Referenzpunkt
    (Ursprung) ausgewertet.

    @param metric_func  Anfangsmetrik auf der Fläche
    @param T            Endzeit
    @param steps        Anzahl der Zeitschritte
    @return             Liste von (t, g_matrix) Paaren
    @author Michael Fuhrmann
    @since  2026-03-11
    """
    flow = RicciFlow(metric_func, dim=2)
    coords = np.array([0.1, 0.1])  # Referenzpunkt
    return flow.evolve(coords, T, steps)


def normalize_metric(metric_matrix: np.ndarray) -> np.ndarray:
    """
    Normalisiert eine Metrik-Matrix so, dass det(g) = 1.

    Normalisierungsformel:
        g_norm = g / det(g)^(1/n)

    @param metric_matrix  Quadratische Metrik-Matrix
    @return               Normalisierte Metrik-Matrix
    @author Michael Fuhrmann
    @since  2026-03-11
    """
    g = np.array(metric_matrix, dtype=float)
    n = g.shape[0]
    det = np.linalg.det(g)
    if abs(det) < 1e-15:
        return g
    return g / (abs(det) ** (1.0 / n))


def compute_volume_element(metric_matrix: np.ndarray) -> float:
    """
    Berechnet das Volumen-Element dV = √det(g).

    @param metric_matrix  Metrik-Matrix g_{ij}
    @return               Volumen-Element √det(g)
    @author Michael Fuhrmann
    @since  2026-03-11
    """
    g = np.array(metric_matrix, dtype=float)
    det = np.linalg.det(g)
    if det < 0:
        return 0.0
    return float(np.sqrt(det))


def ricci_flow_convergence_rate(curvature_history: list) -> float:
    """
    Schätzt die Konvergenzrate des Ricci-Flusses.

    Verwendet lineare Regression im Log-Raum auf die Krümmungsvariation:
        |R(t)| ≈ C · e^{-λt}

    @param curvature_history  Liste von Krümmungswerten über die Zeit
    @return                   Geschätzte Konvergenzrate λ (positiv = Konvergenz)
    @author Michael Fuhrmann
    @since  2026-03-11
    """
    if len(curvature_history) < 3:
        return 0.0

    values = np.array([abs(r) for r in curvature_history], dtype=float)

    # Verhindere log(0)
    values = np.maximum(values, 1e-15)

    # Lineare Regression: log(|R(t)|) = log(C) - λ·t
    t = np.arange(len(values), dtype=float)
    log_vals = np.log(values)

    # Kleinste-Quadrate-Schätzung
    t_mean = np.mean(t)
    log_mean = np.mean(log_vals)
    slope = np.sum((t - t_mean) * (log_vals - log_mean)) / (
        np.sum((t - t_mean) ** 2) + 1e-15
    )

    # Negative Steigung = Konvergenz; λ = -slope
    return float(-slope)
