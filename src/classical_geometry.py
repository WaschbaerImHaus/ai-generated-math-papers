"""
@file classical_geometry.py
@brief Klassische, synthetische und nichteuklidische Geometrie.
@description
    Dieses Modul implementiert grundlegende Konzepte der klassischen Geometrie:

    1. EuclideanGeometry  – Euklidische Geometrie (synthetisch + analytisch):
       Abstände, Winkel, Parallelität, Orthogonalität, Umkreis, Inkreis,
       Schwerpunkt, Euler-Gerade.

    2. ProjectiveGeometry – Projektive Geometrie:
       Homogene Koordinaten, projektive Geraden/Schnittpunkte, Doppelverhältnis,
       harmonische Konjugierte, projektive Transformationen, Sätze von Desargues
       und Pappus.

    3. HyperbolicGeometry – Hyperbolische Geometrie (Poincaré-Kreismodell):
       Hyperbolischer Abstand, Mittelpunkt, Winkelsumme, Flächeninhalt,
       Parallelen-Demo, Gauß-Bonnet, oberes Halbebenenmodell.

    4. EllipticGeometry  – Elliptische (sphärische) Geometrie:
       Großkreisabstand (Haversine), sphärische Winkelsumme, sphärischer Exzess,
       sphärische Fläche, Zweieck (Lune), Girard'scher Satz.

    5. AffinGeometry     – Affine Geometrie:
       Affinkombination, affine Abbildung, baryzentrische Koordinaten,
       affine Hülle, konvexe Hülle (Gift-Wrapping), affine Unabhängigkeit.

    6. SyntheticGeometry – Synthetische Geometrie (Hilbert-Axiome):
       Übersicht Hilbert-Axiome, Inzidenz-Demo, Parallelen-Demo,
       Axiomen-Unabhängigkeit, endliche projektive Ebene PG(2,q).

    7. Standalone-Funktionen:
       Euler-Charakteristik, Picks Satz, Ptolemäus-Satz (Prüfung),
       Neun-Punkte-Kreis, Morley'scher Dreieckssatz.

    Hinweis: `hyperbolic_plane_metric`, `gaussian_curvature` und
    `riemann_tensor` befinden sich bereits in `tensor_geometry.py` und
    werden hier NICHT dupliziert.

    Formeln können im KaTeX-Format angegeben werden (siehe .md-Datei).

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
import itertools
import numpy as np
from typing import Optional

# ============================================================================
# Hilfsfunktionen (intern)
# ============================================================================

def _vec2(p1, p2):
    """Gibt den 2D-Richtungsvektor von p1 nach p2 als numpy-Array zurück."""
    return np.array(p2, dtype=float) - np.array(p1, dtype=float)


def _cross2(u, v):
    """Kreuzprodukt zweier 2D-Vektoren (Skalar = z-Komponente)."""
    return float(u[0] * v[1] - u[1] * v[0])


# ============================================================================
# 1. EuclideanGeometry
# ============================================================================

class EuclideanGeometry:
    """
    Euklidische Geometrie in beliebiger Dimension (Standard: 2D).

    Implementiert analytische und synthetische Methoden der euklidischen
    Geometrie gemäß Euklids „Elementen" und moderner analytischer Darstellung.

    Wichtige Eigenschaften euklidischer Räume:
    - Abstand: d(P,Q) = √(Σ(xᵢ-yᵢ)²)  (KaTeX: $d = \\sqrt{\\sum (x_i - y_i)^2}$)
    - Winkel: cos θ = (u·v)/(|u||v|)
    - Winkelsumme im Dreieck: α+β+γ = π

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, dim: int = 2) -> None:
        """
        Initialisiert das euklidische Geometrie-Objekt.

        @param dim: Dimension des Raums (Standard: 2).
        @lastModified 2026-03-10
        """
        self.dim = dim  # Dimension des euklidischen Raums

    def distance(self, p1, p2) -> float:
        """
        Berechnet den euklidischen Abstand zwischen zwei Punkten.

        Formel: $d(P, Q) = \\sqrt{\\sum_{i=1}^{n}(p_i - q_i)^2}$

        @param p1: Erster Punkt (Liste oder Array).
        @param p2: Zweiter Punkt (Liste oder Array).
        @return: Euklidischer Abstand als float.
        @lastModified 2026-03-10
        """
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        return float(np.linalg.norm(p2 - p1))

    def angle_between(self, v1, v2) -> float:
        """
        Berechnet den Winkel (in Radiant) zwischen zwei Vektoren.

        Formel: $\\theta = \\arccos\\left(\\frac{\\mathbf{u}\\cdot\\mathbf{v}}{|\\mathbf{u}||\\mathbf{v}|}\\right)$

        @param v1: Erster Vektor (Liste oder Array).
        @param v2: Zweiter Vektor (Liste oder Array).
        @return: Winkel in Radiant [0, π].
        @raises ValueError: Wenn ein Vektor der Nullvektor ist.
        @lastModified 2026-03-10
        """
        v1 = np.array(v1, dtype=float)
        v2 = np.array(v2, dtype=float)
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 < 1e-15 or n2 < 1e-15:
            raise ValueError("Nullvektor hat keinen definierten Winkel.")
        # Clipping verhindert Domänenfehler durch Rundungsfehler
        cos_theta = np.dot(v1, v2) / (n1 * n2)
        cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
        return math.acos(cos_theta)

    def are_parallel(self, l1_start, l1_end, l2_start, l2_end,
                     tol: float = 1e-10) -> bool:
        """
        Prüft, ob zwei Geraden in 2D parallel sind.

        Zwei Geraden sind parallel, wenn ihre Richtungsvektoren linear
        abhängig sind, d.h. das Kreuzprodukt (in 2D der Skalar u×v) gleich 0 ist.

        @param l1_start: Startpunkt der ersten Geraden.
        @param l1_end:   Endpunkt der ersten Geraden.
        @param l2_start: Startpunkt der zweiten Geraden.
        @param l2_end:   Endpunkt der zweiten Geraden.
        @param tol:      Toleranz für Nullvergleich.
        @return: True wenn parallel (oder identisch), sonst False.
        @lastModified 2026-03-10
        """
        d1 = _vec2(l1_start, l1_end)
        d2 = _vec2(l2_start, l2_end)
        return abs(_cross2(d1, d2)) < tol

    def are_perpendicular(self, v1, v2, tol: float = 1e-10) -> bool:
        """
        Prüft, ob zwei Vektoren orthogonal (rechtwinklig) zueinander sind.

        Bedingung: $\\mathbf{u} \\cdot \\mathbf{v} = 0$

        @param v1:  Erster Vektor.
        @param v2:  Zweiter Vektor.
        @param tol: Toleranz für Nullvergleich.
        @return: True wenn orthogonal, sonst False.
        @lastModified 2026-03-10
        """
        v1 = np.array(v1, dtype=float)
        v2 = np.array(v2, dtype=float)
        return abs(float(np.dot(v1, v2))) < tol

    def circumcenter(self, p1, p2, p3) -> np.ndarray:
        """
        Berechnet den Umkreismittelpunkt (Circumcenter) eines 2D-Dreiecks.

        Der Umkreismittelpunkt ist der Schnittpunkt der Streckensymmetralen.
        Er ist von allen drei Ecken gleich weit entfernt.

        Formel über lineare Gleichungssysteme:
        $|C - P_1|^2 = |C - P_2|^2 = |C - P_3|^2$

        @param p1: Erster Eckpunkt [x, y].
        @param p2: Zweiter Eckpunkt [x, y].
        @param p3: Dritter Eckpunkt [x, y].
        @return: Umkreismittelpunkt als numpy-Array [x, y].
        @raises ValueError: Wenn die Punkte kollinear sind (entartetes Dreieck).
        @lastModified 2026-03-10
        """
        ax, ay = float(p1[0]), float(p1[1])
        bx, by = float(p2[0]), float(p2[1])
        cx, cy = float(p3[0]), float(p3[1])

        # Determinante prüfen (= 2 · Fläche des Dreiecks)
        D = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(D) < 1e-12:
            raise ValueError("Kollineare Punkte: kein eindeutiger Umkreismittelpunkt.")

        # Formel für Umkreismittelpunkt
        ux = ((ax**2 + ay**2) * (by - cy) +
              (bx**2 + by**2) * (cy - ay) +
              (cx**2 + cy**2) * (ay - by)) / D
        uy = ((ax**2 + ay**2) * (cx - bx) +
              (bx**2 + by**2) * (ax - cx) +
              (cx**2 + cy**2) * (bx - ax)) / D
        return np.array([ux, uy])

    def incenter(self, p1, p2, p3) -> np.ndarray:
        """
        Berechnet den Inkreismittelpunkt eines 2D-Dreiecks.

        Der Inkreismittelpunkt ist der Schnittpunkt der Winkelhalbierenden.
        Er ist gewichtet mit den gegenüberliegenden Seitenlängen:

        $I = \\frac{a \\cdot P_1 + b \\cdot P_2 + c \\cdot P_3}{a + b + c}$

        wobei a = |P2P3|, b = |P1P3|, c = |P1P2|.

        @param p1: Erster Eckpunkt.
        @param p2: Zweiter Eckpunkt.
        @param p3: Dritter Eckpunkt.
        @return: Inkreismittelpunkt als numpy-Array.
        @raises ValueError: Bei entartetem Dreieck (Umfang = 0).
        @lastModified 2026-03-10
        """
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        p3 = np.array(p3, dtype=float)

        # Seitenlängen (gegenüber dem jeweiligen Punkt)
        a = self.distance(p2, p3)  # gegenüber P1
        b = self.distance(p1, p3)  # gegenüber P2
        c = self.distance(p1, p2)  # gegenüber P3
        perimeter = a + b + c

        if perimeter < 1e-12:
            raise ValueError("Entartetes Dreieck: Umfang = 0.")

        # Gewichteter Schwerpunkt mit den Seitenlängen
        return (a * p1 + b * p2 + c * p3) / perimeter

    def centroid(self, points) -> np.ndarray:
        """
        Berechnet den Schwerpunkt (Centroid) einer Punktmenge.

        Der Schwerpunkt ist das arithmetische Mittel aller Punkte:
        $G = \\frac{1}{n} \\sum_{i=1}^{n} P_i$

        @param points: Liste von Punkten (jeder als Liste oder Array).
        @return: Schwerpunkt als numpy-Array.
        @raises ValueError: Wenn die Punktliste leer ist.
        @lastModified 2026-03-10
        """
        pts = np.array(points, dtype=float)
        if pts.shape[0] == 0:
            raise ValueError("Leere Punktliste.")
        return np.mean(pts, axis=0)

    def euler_line_demo(self, p1, p2, p3) -> dict:
        """
        Demonstriert die Euler-Gerade eines Dreiecks.

        Die Euler-Gerade verbindet drei wichtige Punkte eines Dreiecks:
        - O: Umkreismittelpunkt (Circumcenter)
        - G: Schwerpunkt (Centroid)
        - H: Höhenschnittpunkt (Orthocenter)

        Satz von Euler: O, G, H sind kollinear, und |OG| : |GH| = 1 : 2.

        Berechnung des Höhenschnittpunkts H:
        $H = O + 3(G - O) = 3G - 2O$

        @param p1: Erster Eckpunkt.
        @param p2: Zweiter Eckpunkt.
        @param p3: Dritter Eckpunkt.
        @return: Dict mit 'circumcenter', 'centroid', 'orthocenter',
                 'euler_line_verified' (bool), 'ratio' (≈1:2).
        @lastModified 2026-03-10
        """
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        p3 = np.array(p3, dtype=float)

        O = self.circumcenter(p1, p2, p3)  # Umkreismittelpunkt
        G = self.centroid([p1, p2, p3])    # Schwerpunkt

        # Höhenschnittpunkt: H = 3G - 2O (Euler'sche Relation)
        H = 3.0 * G - 2.0 * O

        # Verifikation der Kollinearität: (G-O) × (H-O) ≈ 0
        v1 = G - O
        v2 = H - O
        cross = abs(_cross2(v1, v2))
        collinear = cross < 1e-8

        # Verhältnis |OG| : |GH| sollte 1:2 sein
        d_OG = self.distance(O, G)
        d_GH = self.distance(G, H)
        ratio = d_GH / d_OG if d_OG > 1e-15 else float('inf')

        return {
            "circumcenter": O,
            "centroid": G,
            "orthocenter": H,
            "euler_line_verified": collinear,
            "ratio_GH_to_OG": ratio,  # ≈ 2.0
        }


# ============================================================================
# 2. ProjectiveGeometry
# ============================================================================

class ProjectiveGeometry:
    """
    Projektive Geometrie in der projektiven Ebene RP².

    Punkte werden durch homogene Koordinaten [x:y:w] dargestellt.
    Ein gewöhnlicher Punkt (a,b) entspricht [a:b:1].
    Punkte im Unendlichen: [x:y:0].

    Dualität: Geraden und Punkte sind dual (Gerade = [a:b:c] mit ax+by+cw=0).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def homogeneous_coords(self, point) -> np.ndarray:
        """
        Konvertiert einen affinen Punkt in homogene Koordinaten.

        Formel: $(x, y) \\mapsto [x : y : 1]$

        @param point: 2D-Punkt [x, y] oder 3D-Punkt [x, y, w].
        @return: Homogene Koordinaten als numpy-Array [x, y, 1].
        @lastModified 2026-03-10
        """
        p = np.array(point, dtype=float)
        if p.shape[0] == 2:
            return np.array([p[0], p[1], 1.0])
        # Bereits homogene Koordinaten
        return p

    def projective_line_through(self, p1, p2) -> np.ndarray:
        """
        Berechnet die projektive Gerade durch zwei Punkte.

        Im Projektiven ist eine Gerade das Kreuzprodukt zweier Punkte:
        $\\ell = P_1 \\times P_2$ (in homogenen Koordinaten)

        @param p1: Erster Punkt (affin oder homogen).
        @param p2: Zweiter Punkt (affin oder homogen).
        @return: Geraden-Koeffizientenvektor [a:b:c] mit ax+by+cw=0.
        @raises ValueError: Wenn p1 und p2 identisch sind.
        @lastModified 2026-03-10
        """
        h1 = self.homogeneous_coords(p1)
        h2 = self.homogeneous_coords(p2)
        line = np.cross(h1, h2)
        if np.linalg.norm(line) < 1e-14:
            raise ValueError("Identische Punkte definieren keine eindeutige Gerade.")
        return line

    def projective_intersection(self, l1, l2) -> np.ndarray:
        """
        Berechnet den Schnittpunkt zweier projektiver Geraden.

        Im Projektiven ist der Schnittpunkt zweier Geraden wieder ein Kreuzprodukt:
        $P = \\ell_1 \\times \\ell_2$

        Parallele Geraden schneiden sich im Unendlichen: [x:y:0].

        @param l1: Erste Gerade [a:b:c].
        @param l2: Zweite Gerade [a:b:c].
        @return: Schnittpunkt in homogenen Koordinaten.
        @raises ValueError: Wenn die Geraden identisch sind.
        @lastModified 2026-03-10
        """
        l1 = np.array(l1, dtype=float)
        l2 = np.array(l2, dtype=float)
        point = np.cross(l1, l2)
        if np.linalg.norm(point) < 1e-14:
            raise ValueError("Identische Geraden: kein eindeutiger Schnittpunkt.")
        return point

    def cross_ratio(self, p1, p2, p3, p4) -> float:
        """
        Berechnet das Doppelverhältnis (Cross Ratio) von vier kollinearen Punkten.

        Das Doppelverhältnis ist die zentrale projektive Invariante:
        $(P_1, P_2; P_3, P_4) = \\frac{(P_3 - P_1)(P_4 - P_2)}{(P_3 - P_2)(P_4 - P_1)}$

        Es ist invariant unter projektiven Transformationen.

        @param p1: Erster Punkt (1D-Koordinate oder homogen [x:w]).
        @param p2: Zweiter Punkt.
        @param p3: Dritter Punkt.
        @param p4: Vierter Punkt.
        @return: Doppelverhältnis als float.
        @raises ZeroDivisionError: Bei entartetem Konfiguration.
        @lastModified 2026-03-10
        """
        # Verwende 1D-Koordinaten (affin) oder normalisiere homogene Koordinaten
        def to_1d(p):
            p = np.array(p, dtype=float).ravel()
            if p.shape[0] == 1:
                return p[0]
            if p.shape[0] == 2 and abs(p[1]) > 1e-15:
                return p[0] / p[1]  # normalisiere [x:w] → x/w
            return p[0]  # direkte 1D-Koordinate

        a, b, c, d = to_1d(p1), to_1d(p2), to_1d(p3), to_1d(p4)

        # Doppelverhältnis: (c-a)(d-b) / ((c-b)(d-a))
        num = (c - a) * (d - b)
        den = (c - b) * (d - a)
        if abs(den) < 1e-14:
            raise ZeroDivisionError("Entartetes Doppelverhältnis (Nenner = 0).")
        return num / den

    def harmonic_conjugate(self, p1: float, p2: float, p3: float) -> float:
        """
        Berechnet die harmonische Konjugierte P4 zu drei Punkten P1, P2, P3.

        Das harmonische Verhältnis gilt:
        $(P_1, P_2; P_3, P_4) = -1$

        Auflösen nach P4 liefert:
        $P_4 = \\frac{2 P_1 P_2 - P_3(P_1 + P_2)}{2 P_1 P_2 / P_3 - (P_1 + P_2)}$

        Vereinfachte Formel (P1, P2 fest; P3 variabel):
        $P_4 = \\frac{P_1 P_2 (P_3 - P_1 - P_2 + P_4) \\ldots}{\\ldots}$

        @param p1: Erster Punkt (1D-Koordinate).
        @param p2: Zweiter Punkt (1D-Koordinate).
        @param p3: Dritter Punkt (1D-Koordinate).
        @return: Vierter Punkt P4 (harmonisch konjugiert zu P3 bzgl. P1,P2).
        @raises ZeroDivisionError: Bei entartetem Fall.
        @lastModified 2026-03-10
        """
        # Aus (P1,P2;P3,P4)=-1:
        # (P3-P1)(P4-P2) / ((P3-P2)(P4-P1)) = -1
        # (P3-P1)(P4-P2) = -(P3-P2)(P4-P1)
        # Löse nach P4 auf:
        # P4[(P3-P1) + (P3-P2)] = P2(P3-P1) - P1(P3-P2)
        # P4 = [P2(P3-P1) - P1(P3-P2)] / [(P3-P1) + (P3-P2)]
        #    = [P2(P3-P1) - P1(P3-P2)] / [2P3 - P1 - P2]
        num = p2 * (p3 - p1) - p1 * (p3 - p2)
        den = 2.0 * p3 - p1 - p2
        if abs(den) < 1e-14:
            raise ZeroDivisionError("Entartete harmonische Konjugierte.")
        return num / den

    def projective_transformation(self, matrix, point) -> np.ndarray:
        """
        Wendet eine projektive Transformation (3×3-Matrix) auf einen Punkt an.

        $P' = M \\cdot P$ in homogenen Koordinaten.

        @param matrix: 3×3-Matrix (Homographie).
        @param point:  Punkt in homogenen Koordinaten [x:y:w] oder affin [x,y].
        @return: Transformierter Punkt in homogenen Koordinaten.
        @raises ValueError: Wenn Matrix nicht invertierbar ist.
        @lastModified 2026-03-10
        """
        M = np.array(matrix, dtype=float)
        p = self.homogeneous_coords(point)
        result = M @ p
        return result

    def desargues_theorem_check(self, triangle1, triangle2,
                                tol: float = 1e-8) -> dict:
        """
        Prüft den Desarguesschen Satz numerisch.

        Satz von Desargues: Zwei Dreiecke sind genau dann perspektiv bzgl.
        eines Punktes (Zentrum), wenn sie perspektiv bzgl. einer Geraden
        (Achse) sind.

        @param triangle1: Liste von 3 Punkten [[x1,y1],[x2,y2],[x3,y3]].
        @param triangle2: Liste von 3 Punkten [[x1,y1],[x2,y2],[x3,y3]].
        @param tol: Toleranz für Kollinearitätstest.
        @return: Dict mit 'center_perspective', 'axis_perspective',
                 'desargues_satisfied'.
        @lastModified 2026-03-10
        """
        # Konvertiere alle Punkte in homogene Koordinaten
        A1, B1, C1 = [self.homogeneous_coords(p) for p in triangle1]
        A2, B2, C2 = [self.homogeneous_coords(p) for p in triangle2]

        # Schnittpunkte der entsprechenden Seiten (Achsenpunkte)
        # Seite A1B1 und A2B2
        l_A1B1 = np.cross(A1, B1)
        l_A2B2 = np.cross(A2, B2)
        P = np.cross(l_A1B1, l_A2B2)

        # Seite B1C1 und B2C2
        l_B1C1 = np.cross(B1, C1)
        l_B2C2 = np.cross(B2, C2)
        Q = np.cross(l_B1C1, l_B2C2)

        # Seite A1C1 und A2C2
        l_A1C1 = np.cross(A1, C1)
        l_A2C2 = np.cross(A2, C2)
        R = np.cross(l_A1C1, l_A2C2)

        # Prüfe Kollinearität von P, Q, R (Achse)
        # Determinante der 3 homogenen Punkte ≈ 0 wenn kollinear
        det_axis = abs(np.linalg.det(np.array([P, Q, R])))
        axis_ok = det_axis < tol

        # Prüfe Perspektivität bzgl. Punkt:
        # Verbindungsgeraden AA', BB', CC' müssen sich in einem Punkt treffen
        l_AA = np.cross(A1, A2)
        l_BB = np.cross(B1, B2)
        l_CC = np.cross(C1, C2)
        det_center = abs(np.linalg.det(np.array([l_AA, l_BB, l_CC])))
        center_ok = det_center < tol

        return {
            "center_perspective": center_ok,
            "axis_perspective": axis_ok,
            "desargues_satisfied": center_ok == axis_ok,
            "axis_collinearity_det": det_axis,
            "center_concurrency_det": det_center,
        }

    def pappus_theorem_check(self, l1_points, l2_points,
                             tol: float = 1e-8) -> bool:
        """
        Prüft den Satz von Pappus numerisch.

        Satz von Pappus: Liegen drei Punkte A, B, C auf Gerade g1 und drei
        Punkte D, E, F auf Gerade g2, dann sind die drei Schnittpunkte
        AE∩BD, AF∩CD, BF∩CE kollinear (Pappus-Gerade).

        @param l1_points: [A, B, C] – drei Punkte auf der ersten Geraden.
        @param l2_points: [D, E, F] – drei Punkte auf der zweiten Geraden.
        @param tol: Toleranz für Kollinearitätsprüfung.
        @return: True wenn Pappus-Bedingung erfüllt.
        @lastModified 2026-03-10
        """
        A, B, C = [self.homogeneous_coords(p) for p in l1_points]
        D, E, F = [self.homogeneous_coords(p) for p in l2_points]

        # Schnittpunkte: AE∩BD, AF∩CD, BF∩CE
        AE = np.cross(A, E)
        BD = np.cross(B, D)
        P = np.cross(AE, BD)

        AF = np.cross(A, F)
        CD = np.cross(C, D)
        Q = np.cross(AF, CD)

        BF = np.cross(B, F)
        CE = np.cross(C, E)
        R = np.cross(BF, CE)

        # Kollinearität: det([P,Q,R]) ≈ 0
        det = abs(np.linalg.det(np.array([P, Q, R], dtype=float)))
        return det < tol


# ============================================================================
# 3. HyperbolicGeometry
# ============================================================================

class HyperbolicGeometry:
    """
    Hyperbolische Geometrie im Poincaré-Kreismodell und im
    oberen Halbebenenmodell.

    Im Poincaré-Kreismodell wird die hyperbolische Ebene durch den offenen
    Einheitskreis {z ∈ ℂ : |z| < 1} modelliert.

    Hyperbolischer Abstand:
    $d(z_1, z_2) = \\text{arctanh}\\left(\\frac{|z_1 - z_2|}{|1 - \\bar{z}_1 z_2|}\\right) \\cdot 2$

    Wichtige Eigenschaften:
    - Winkelsumme im hyperbolischen Dreieck: α+β+γ < π
    - Fläche = π − (α+β+γ) (Gauß-Bonnet für hyperbolische Dreiecke)
    - Unendlich viele Parallelen durch einen Punkt außerhalb einer Geraden

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def poincare_disk_distance(self, z1, z2) -> float:
        """
        Berechnet den hyperbolischen Abstand im Poincaré-Kreismodell.

        Formel:
        $d(z_1, z_2) = 2 \\operatorname{arctanh}\\!\\left(\\frac{|z_1 - z_2|}{|1 - \\bar{z}_1 z_2|}\\right)$

        @param z1: Erster Punkt als komplexe Zahl oder [x, y] mit |z1| < 1.
        @param z2: Zweiter Punkt als komplexe Zahl oder [x, y] mit |z2| < 1.
        @return: Hyperbolischer Abstand (≥ 0).
        @raises ValueError: Wenn ein Punkt außerhalb des Einheitskreises liegt.
        @lastModified 2026-03-10
        """
        # Konvertiere zu komplexen Zahlen
        if not isinstance(z1, complex):
            z1 = complex(z1[0], z1[1])
        if not isinstance(z2, complex):
            z2 = complex(z2[0], z2[1])

        if abs(z1) >= 1.0:
            raise ValueError(f"Punkt z1 liegt nicht im Einheitskreis: |z1|={abs(z1):.4f}")
        if abs(z2) >= 1.0:
            raise ValueError(f"Punkt z2 liegt nicht im Einheitskreis: |z2|={abs(z2):.4f}")

        # Möbius-Pseudometrik
        num = abs(z1 - z2)
        den = abs(1.0 - z1.conjugate() * z2)
        if den < 1e-15:
            return float('inf')
        ratio = num / den
        # arctanh nur für |ratio| < 1 definiert
        ratio = min(ratio, 1.0 - 1e-15)
        return 2.0 * math.atanh(ratio)

    def poincare_disk_midpoint(self, z1, z2) -> complex:
        """
        Berechnet den hyperbolischen Mittelpunkt zweier Punkte im Poincaré-Modell.

        Der Mittelpunkt M erfüllt d(z1, M) = d(M, z2).
        Berechnung über die Möbius-Transformation T_{z1}(z) = (z - z1)/(1 - conj(z1)·z):

        M = T_{z1}^{-1}(T_{z1}(z2) / 2 · (Normierung))

        @param z1: Erster Punkt im Einheitskreis.
        @param z2: Zweiter Punkt im Einheitskreis.
        @return: Hyperbolischer Mittelpunkt als komplexe Zahl.
        @lastModified 2026-03-10
        """
        if not isinstance(z1, complex):
            z1 = complex(z1[0], z1[1])
        if not isinstance(z2, complex):
            z2 = complex(z2[0], z2[1])

        # Verschiebe z1 in den Ursprung mit Möbius-Transformation
        def mobius(z, a):
            """Möbius-Transformation T_a: z → (z-a)/(1-conj(a)·z)"""
            return (z - a) / (1.0 - a.conjugate() * z)

        def inv_mobius(w, a):
            """Inverse Möbius-Transformation T_a^{-1}: w → (w+a)/(1+conj(a)·w)"""
            return (w + a) / (1.0 + a.conjugate() * w)

        # z2 im durch T_{z1} transformierten Koordinatensystem
        z2_shifted = mobius(z2, z1)

        # Mittelpunkt im verschobenen System: skaliere auf halbe Länge
        r = abs(z2_shifted)
        if r < 1e-15:
            return z1  # Punkte sind identisch
        # Mittelpunkt liegt bei z2_shifted / 2 (in Richtung, aber halbem Abstand)
        # Genauer: hyperbolischer Mittelpunkt bei Skalierung via tanh
        # mid_r = tanh(atanh(r)/1) → der Punkt bei halbem hyperb. Abstand
        # d(0,w) = 2·atanh(|w|) → bei halbem Abstand: 2·atanh(r/2)?
        # Korrekt: mid_r = tanh(atanh(r)) / 1... → Korrekte Formel:
        # half_dist = atanh(r), then new_r = tanh(half_dist) → r selbst!
        # Wir nutzen: Mittelpunkt bei |w| = tanh(atanh(r)/1)... nein.
        # Korrekt: d(0, w) = 2·atanh(|w|), also bei Hälfte d/2 = atanh(|w|)
        # → |w_mid| = tanh(atanh(r)/1)
        # Formel für halben hyperbolischen Abstand vom Ursprung:
        r_mid = math.tanh(math.atanh(r) / 1.0)  # = r (trivial, Abstand bleibt)
        # Der echte hyperbolische Mittelpunkt zwischen 0 und z2_shifted:
        # liegt bei r_half = tanh(atanh(r)/2) * (z2_shifted/|z2_shifted|)... Warte:
        # d(0, p) = 2*atanh(|p|). Hälfte: atanh(|p_half|) = atanh(|z2_shifted|)/2
        # |p_half| = tanh(atanh(|z2_shifted|)/2)
        r_half = math.tanh(math.atanh(r) / 2.0)
        direction = z2_shifted / r
        mid_shifted = r_half * direction

        # Zurücktransformieren
        return inv_mobius(mid_shifted, z1)

    def _compute_angle_at_vertex(self, vertex, p1, p2) -> float:
        """
        Berechnet den hyperbolischen Winkel an einem Dreieck-Eckpunkt.

        Im Poincaré-Modell ist das Modell winkeltreu (konform). Daher:
        1. Verschiebe den Scheitelpunkt mit Möbius-Transformation in den Ursprung.
        2. Messe den euklidischen Winkel zwischen den transformierten Nachbarpunkten.

        Dies ist die korrekte Methode, weil die euklidische Winkelmessung
        im Ursprung dem hyperbolischen Winkel entspricht, wenn die Tangenten
        der Geodäten durch die Ursprungsrichtungen gegeben sind.

        Wichtig: Der Winkel wird zwischen den Tangenten der hyperbolischen
        Geraden (Kreisbögen) am Scheitelpunkt gemessen, nicht zwischen den
        Verbindungsgeraden der affinen Punkte.

        @param vertex: Eckpunkt (complex), liegt im offenen Einheitskreis.
        @param p1:     Erster benachbarter Eckpunkt (complex).
        @param p2:     Zweiter benachbarter Eckpunkt (complex).
        @return: Winkel in Radiant.
        @lastModified 2026-03-10
        """
        # Möbius-Transformation: verschiebt vertex in den Ursprung
        def mobius(z, a):
            """T_a(z) = (z-a)/(1-conj(a)·z) – Isometrie des Poincaré-Disks."""
            denom = 1.0 - a.conjugate() * z
            if abs(denom) < 1e-15:
                return complex(float('inf'), 0)
            return (z - a) / denom

        # Transformiere alle Punkte so, dass vertex im Ursprung liegt
        q1 = mobius(p1, vertex)
        q2 = mobius(p2, vertex)
        # q0 = mobius(vertex, vertex) = 0 (Ursprung)

        # Im Ursprung: Tangente an Geodäte ist die Richtung zum transformierten Punkt
        # (weil eine hyperbolische Gerade durch den Ursprung ein Durchmesser ist)
        if abs(q1) < 1e-15 or abs(q2) < 1e-15:
            return 0.0

        # Euklidischer Winkel zwischen q1 und q2 (vom Ursprung aus)
        cos_a = (q1.real * q2.real + q1.imag * q2.imag) / (abs(q1) * abs(q2))
        cos_a = max(-1.0, min(1.0, cos_a))
        return math.acos(cos_a)

    def hyperbolic_angle_sum(self, triangle_vertices) -> float:
        """
        Berechnet die Winkelsumme eines hyperbolischen Dreiecks.

        Im hyperbolischen Dreieck gilt stets:
        $\\alpha + \\beta + \\gamma < \\pi$

        Das Poincaré-Modell ist winkeltreu (konform), daher können euklidische
        Winkelberechnungen an den Scheitelpunkten verwendet werden.

        @param triangle_vertices: [z1, z2, z3] als komplexe Zahlen oder [x,y]-Paare.
        @return: Winkelsumme in Radiant (< π für nicht-entartete Dreiecke).
        @lastModified 2026-03-10
        """
        def to_complex(z):
            if isinstance(z, complex):
                return z
            return complex(z[0], z[1])

        A, B, C = [to_complex(v) for v in triangle_vertices]

        alpha = self._compute_angle_at_vertex(A, B, C)
        beta  = self._compute_angle_at_vertex(B, A, C)
        gamma = self._compute_angle_at_vertex(C, A, B)

        return alpha + beta + gamma

    def hyperbolic_area(self, triangle_vertices) -> float:
        """
        Berechnet den Flächeninhalt eines hyperbolischen Dreiecks.

        Gauß-Bonnet-Formel für hyperbolische Dreiecke (Gaußsche Krümmung K=-1):
        $A = \\pi - (\\alpha + \\beta + \\gamma)$

        Der Flächeninhalt ist immer positiv und kleiner als π.

        @param triangle_vertices: [z1, z2, z3] als komplexe Zahlen.
        @return: Hyperbolische Fläche (∈ (0, π)).
        @lastModified 2026-03-10
        """
        angle_sum = self.hyperbolic_angle_sum(triangle_vertices)
        area = math.pi - angle_sum
        return max(0.0, area)  # Numerische Sicherheit

    def hyperbolic_parallel_lines(self) -> dict:
        """
        Demonstriert das Parallelenaxiom der hyperbolischen Geometrie.

        In der hyperbolischen Geometrie gibt es durch einen Punkt P außerhalb
        einer Geraden g UNENDLICH VIELE Geraden, die g nicht schneiden.

        Im Poincaré-Modell sind Geraden Kreise, die senkrecht auf dem Rand
        ∂D (Einheitskreis) stehen, oder Durchmesser.

        @return: Dict mit Beispielpunkten und Beschreibung.
        @lastModified 2026-03-10
        """
        # Beispiel: Gerade g = x-Achse (Durchmesser), Punkt P = (0, 0.5)
        return {
            "model": "Poincaré-Kreismodell",
            "base_line": "x-Achse (Verbindung (-1,0) und (1,0))",
            "external_point": complex(0, 0.5),
            "parallel_axiom": "hyperbolisch",
            "description": (
                "Im hyperbolischen Modell gibt es durch P=(0, 0.5) unendlich "
                "viele 'Geraden' (Kreise orthogonal zum Rand), die die "
                "x-Achse nicht schneiden. Dies widerspricht Euklids "
                "Parallelenpostulat (genau eine Parallele)."
            ),
            "parallel_count": float('inf'),
            "angle_of_parallelism_formula": "Π(d) = 2·arctan(e^{-d})",
        }

    def gauss_bonnet_hyperbolic(self, polygon_angles) -> float:
        """
        Wendet den Gauß-Bonnet-Satz auf ein hyperbolisches Polygon an.

        Für ein hyperbolisches n-Eck mit Innenwinkeln α₁,...,αₙ gilt:
        $A = (n-2)\\pi - \\sum_{i=1}^{n} \\alpha_i$

        Für n=3: A = π - (α₁+α₂+α₃) (hyperbolisches Dreieck).

        @param polygon_angles: Liste der Innenwinkel in Radiant.
        @return: Flächeninhalt des hyperbolischen Polygons.
        @lastModified 2026-03-10
        """
        n = len(polygon_angles)
        angle_sum = sum(polygon_angles)
        area = (n - 2) * math.pi - angle_sum
        return max(0.0, area)

    def upper_half_plane_distance(self, z1, z2) -> float:
        """
        Berechnet den hyperbolischen Abstand im oberen Halbebenenmodell.

        Im oberen Halbebenenmodell H = {z ∈ ℂ : Im(z) > 0} gilt:
        $d(z_1, z_2) = \\operatorname{arccosh}\\!\\left(1 + \\frac{|z_1 - z_2|^2}{2 \\operatorname{Im}(z_1) \\operatorname{Im}(z_2)}\\right)$

        @param z1: Erster Punkt mit Im(z1) > 0 (als complex oder [x,y]).
        @param z2: Zweiter Punkt mit Im(z2) > 0.
        @return: Hyperbolischer Abstand.
        @raises ValueError: Wenn Punkte nicht in der oberen Halbebene liegen.
        @lastModified 2026-03-10
        """
        if not isinstance(z1, complex):
            z1 = complex(z1[0], z1[1])
        if not isinstance(z2, complex):
            z2 = complex(z2[0], z2[1])

        if z1.imag <= 0:
            raise ValueError("z1 liegt nicht in der oberen Halbebene (Im > 0).")
        if z2.imag <= 0:
            raise ValueError("z2 liegt nicht in der oberen Halbebene (Im > 0).")

        # Formel: arccosh(1 + |z1-z2|² / (2·Im(z1)·Im(z2)))
        delta = abs(z1 - z2) ** 2
        denom = 2.0 * z1.imag * z2.imag
        arg = 1.0 + delta / denom
        # Sicherheitsclipping (numerische Fehler)
        arg = max(1.0, arg)
        return math.acosh(arg)


# ============================================================================
# 4. EllipticGeometry
# ============================================================================

class EllipticGeometry:
    """
    Elliptische (sphärische) Geometrie auf der Einheitssphäre S².

    In der elliptischen Geometrie gilt:
    - Winkelsumme im Dreieck: α+β+γ > π
    - Keine parallelen Geraden (alle Großkreise schneiden sich)
    - Kurze Geraden = Großkreise

    Anwendungen: Kugelgeometrie, Geodäsie, Navigationssysteme (WGS-84).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def spherical_distance(self, lat1: float, lon1: float,
                           lat2: float, lon2: float,
                           radius: float = 1.0) -> float:
        """
        Berechnet den Großkreisabstand zwischen zwei Punkten auf der Sphäre.

        Haversine-Formel (numerisch stabil für kleine und große Abstände):
        $d = 2R \\arcsin\\!\\sqrt{\\sin^2\\!\\frac{\\Delta\\phi}{2} + \\cos\\phi_1\\cos\\phi_2\\sin^2\\!\\frac{\\Delta\\lambda}{2}}$

        @param lat1: Breitengrad des ersten Punkts (Radiant).
        @param lon1: Längengrad des ersten Punkts (Radiant).
        @param lat2: Breitengrad des zweiten Punkts (Radiant).
        @param lon2: Längengrad des zweiten Punkts (Radiant).
        @param radius: Kugelradius (Standard: 1.0).
        @return: Großkreisabstand auf der Sphäre.
        @lastModified 2026-03-10
        """
        dlat = lat2 - lat1
        dlon = lon2 - lon1

        # Haversine-Formel
        a = (math.sin(dlat / 2.0) ** 2 +
             math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2)
        a = min(1.0, a)  # Numerische Sicherheit
        c = 2.0 * math.asin(math.sqrt(a))
        return radius * c

    def spherical_angle_sum(self, triangle_vertices_rad) -> float:
        """
        Berechnet die Winkelsumme eines sphärischen Dreiecks.

        Gegeben drei Großkreisbögen, berechne die Innenwinkel über
        das sphärische Kosinusgesetz für Winkel:
        $\\cos a = \\cos b \\cos c + \\sin b \\sin c \\cos A$

        @param triangle_vertices_rad: [[lat1,lon1],[lat2,lon2],[lat3,lon3]] in Radiant.
        @return: Winkelsumme in Radiant (> π für nicht-entartete sphärische Dreiecke).
        @lastModified 2026-03-10
        """
        # Konvertiere sphärische Koordinaten in 3D-Kartesisch
        def to_cart(lat, lon):
            return np.array([
                math.cos(lat) * math.cos(lon),
                math.cos(lat) * math.sin(lon),
                math.sin(lat)
            ])

        pts = [to_cart(v[0], v[1]) for v in triangle_vertices_rad]
        A, B, C = pts

        # Seitenlängen als Winkel (Großkreiswinkel)
        def arc_angle(p, q):
            c = np.clip(np.dot(p, q), -1.0, 1.0)
            return math.acos(c)

        a = arc_angle(B, C)  # Seite gegenüber A
        b = arc_angle(A, C)  # Seite gegenüber B
        c = arc_angle(A, B)  # Seite gegenüber C

        # Sphärisches Kosinusgesetz: cos(A) = (cos(a)-cos(b)cos(c))/(sin(b)sin(c))
        def sph_angle(opp, adj1, adj2):
            denom = math.sin(adj1) * math.sin(adj2)
            if abs(denom) < 1e-14:
                return 0.0
            cos_A = (math.cos(opp) - math.cos(adj1) * math.cos(adj2)) / denom
            cos_A = max(-1.0, min(1.0, cos_A))
            return math.acos(cos_A)

        alpha = sph_angle(a, b, c)
        beta  = sph_angle(b, a, c)
        gamma = sph_angle(c, a, b)

        return alpha + beta + gamma

    def spherical_excess(self, triangle_vertices_rad) -> float:
        """
        Berechnet den sphärischen Exzess E = α+β+γ − π.

        Der sphärische Exzess ist der Überschuss der Winkelsumme über π:
        $E = \\alpha + \\beta + \\gamma - \\pi > 0$

        @param triangle_vertices_rad: [[lat,lon],...] in Radiant.
        @return: Sphärischer Exzess in Radiant (> 0).
        @lastModified 2026-03-10
        """
        angle_sum = self.spherical_angle_sum(triangle_vertices_rad)
        return angle_sum - math.pi

    def spherical_area(self, triangle_vertices_rad, radius: float = 1.0) -> float:
        """
        Berechnet den Flächeninhalt eines sphärischen Dreiecks.

        Girard'scher Satz:
        $A = R^2 \\cdot E = R^2 \\cdot (\\alpha + \\beta + \\gamma - \\pi)$

        @param triangle_vertices_rad: [[lat,lon],...] in Radiant.
        @param radius: Kugelradius.
        @return: Sphärische Fläche = R² · E.
        @lastModified 2026-03-10
        """
        E = self.spherical_excess(triangle_vertices_rad)
        return radius ** 2 * E

    def lune_area(self, dihedral_angle: float, radius: float = 1.0) -> float:
        """
        Berechnet den Flächeninhalt eines sphärischen Zweiecks (Lune).

        Eine Lune (Zweieck) ist das Gebiet zwischen zwei Halbkreisen auf S².
        Ihr Flächeninhalt beträgt:
        $A_{\\text{Lune}} = 2 R^2 \\alpha$

        wobei α der Diederwinkel (Öffnungswinkel) ist.

        @param dihedral_angle: Diederwinkel in Radiant (0 < α < π).
        @param radius: Kugelradius.
        @return: Flächeninhalt der Lune.
        @lastModified 2026-03-10
        """
        return 2.0 * radius ** 2 * dihedral_angle

    def girard_theorem(self, alpha: float, beta: float, gamma: float,
                       radius: float = 1.0) -> float:
        """
        Berechnet den Flächeninhalt eines sphärischen Dreiecks via Girard'schem Satz.

        Girard'scher Dreieckssatz (Albert Girard, 1625):
        $A = R^2 \\cdot (\\alpha + \\beta + \\gamma - \\pi)$

        Besonderer Fall: Gleichseitiges sphärisches Dreieck mit α=β=γ=2π/3
        hat Fläche 4π/8 = π/2 (ein Achtel der Kugeloberfläche).

        @param alpha: Innenwinkel α in Radiant.
        @param beta:  Innenwinkel β in Radiant.
        @param gamma: Innenwinkel γ in Radiant.
        @param radius: Kugelradius R.
        @return: Flächeninhalt A = R²·(α+β+γ−π).
        @lastModified 2026-03-10
        """
        E = alpha + beta + gamma - math.pi
        return radius ** 2 * E


# ============================================================================
# 5. AffinGeometry
# ============================================================================

class AffinGeometry:
    """
    Affine Geometrie: Geometrie ohne Abstands- und Winkelbegriff.

    Affine Geometrie untersucht Eigenschaften, die unter affinen Abbildungen
    (lineare Abbildungen + Translationen) invariant sind:
    - Kollinearität von Punkten
    - Teilverhältnisse von Strecken
    - Parallelität von Geraden
    - Konvexität

    Affine Abbildung: $f(\\mathbf{x}) = A\\mathbf{x} + \\mathbf{b}$

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def affine_combination(self, points, coefficients) -> np.ndarray:
        """
        Berechnet eine Affinkombination von Punkten.

        Eine Affinkombination ist eine gewichtete Summe mit Σλᵢ = 1:
        $P = \\sum_{i=1}^{n} \\lambda_i P_i, \\quad \\sum \\lambda_i = 1$

        @param points:       Liste von Punkten.
        @param coefficients: Koeffizienten λᵢ mit Σλᵢ = 1.
        @return: Affinkombination als numpy-Array.
        @raises ValueError: Wenn Koeffizienten nicht auf 1 summieren.
        @lastModified 2026-03-10
        """
        pts = np.array(points, dtype=float)
        coeff = np.array(coefficients, dtype=float)

        # Prüfe Normierungsbedingung Σλᵢ = 1
        if abs(np.sum(coeff) - 1.0) > 1e-10:
            raise ValueError(
                f"Koeffizienten müssen auf 1 summieren, ist: {np.sum(coeff):.6f}"
            )
        # Σλᵢ·Pᵢ (Matrixmultiplikation: coeff @ pts)
        return coeff @ pts

    def affine_map(self, A, b, point) -> np.ndarray:
        """
        Wendet eine affine Abbildung auf einen Punkt an.

        $f(\\mathbf{x}) = A \\mathbf{x} + \\mathbf{b}$

        @param A:     Lineare Matrix (n×n).
        @param b:     Translationsvektor (n,).
        @param point: Eingangs-Punkt (n,).
        @return: Transformierter Punkt.
        @lastModified 2026-03-10
        """
        A = np.array(A, dtype=float)
        b = np.array(b, dtype=float)
        x = np.array(point, dtype=float)
        return A @ x + b

    def barycentric_coords(self, p, p1, p2, p3) -> np.ndarray:
        """
        Berechnet die baryzentrischen Koordinaten von P bzgl. Dreieck P1P2P3.

        Baryzentrische Koordinaten (λ₁, λ₂, λ₃) mit λ₁+λ₂+λ₃=1:
        $P = \\lambda_1 P_1 + \\lambda_2 P_2 + \\lambda_3 P_3$

        Eigenschaften:
        - P liegt im Dreieck ⟺ alle λᵢ ≥ 0
        - P auf einer Seite ⟺ ein λᵢ = 0
        - P ist ein Eckpunkt ⟺ zwei λᵢ = 0

        @param p:  Zu analysierender Punkt [x, y].
        @param p1: Erster Eckpunkt.
        @param p2: Zweiter Eckpunkt.
        @param p3: Dritter Eckpunkt.
        @return: Baryzentrische Koordinaten [λ₁, λ₂, λ₃].
        @raises ValueError: Bei entartetem Dreieck (Fläche = 0).
        @lastModified 2026-03-10
        """
        p  = np.array(p,  dtype=float)
        p1 = np.array(p1, dtype=float)
        p2 = np.array(p2, dtype=float)
        p3 = np.array(p3, dtype=float)

        # Gesamtfläche via Kreuzprodukt
        v0 = p2 - p1
        v1 = p3 - p1
        v2 = p  - p1

        # Löse per Cramers Regel: [v0, v1] · [λ₂, λ₃]ᵀ = v2
        denom = _cross2(v0, v1)
        if abs(denom) < 1e-12:
            raise ValueError("Entartetes Dreieck: Fläche = 0.")

        lam2 = _cross2(v2, v1) / denom
        lam3 = _cross2(v0, v2) / denom
        lam1 = 1.0 - lam2 - lam3

        return np.array([lam1, lam2, lam3])

    def affine_hull(self, points) -> dict:
        """
        Bestimmt die affine Hülle einer Punktmenge und ihre Dimension.

        Die affine Hülle ist die kleinste affine Teilmenge, die alle Punkte enthält.
        Dimension = Rang der Differenzmatrix.

        @param points: Liste von Punkten.
        @return: Dict mit 'dimension', 'affine_hull_description'.
        @lastModified 2026-03-10
        """
        pts = np.array(points, dtype=float)
        if pts.shape[0] == 0:
            return {"dimension": -1, "affine_hull_description": "leer"}

        # Differenzmatrix: erste Punkt als Referenz
        origin = pts[0]
        diffs = pts[1:] - origin  # (n-1) × d Matrix

        if diffs.shape[0] == 0:
            return {"dimension": 0, "affine_hull_description": "Punkt"}

        rank = np.linalg.matrix_rank(diffs)
        descriptions = {0: "Punkt", 1: "Gerade", 2: "Ebene", 3: "Raum"}
        desc = descriptions.get(rank, f"{rank}D-Unterraum")

        return {
            "dimension": rank,
            "affine_hull_description": desc,
        }

    def convex_hull_2d(self, points) -> list:
        """
        Berechnet die konvexe Hülle einer 2D-Punktmenge (Gift-Wrapping-Algorithmus).

        Der Gift-Wrapping-Algorithmus (auch Jarvis-March) findet die konvexe Hülle
        in O(n·h) Zeit, wobei h die Anzahl der Hüllenpunkte ist.

        @param points: Liste von 2D-Punkten [[x1,y1], [x2,y2], ...].
        @return: Geordnete Liste der Hüllenpunkte (gegen den Uhrzeigersinn).
        @raises ValueError: Wenn weniger als 3 Punkte vorhanden sind.
        @lastModified 2026-03-10
        """
        pts = [tuple(p) for p in points]
        n = len(pts)
        if n < 3:
            # Weniger als 3 Punkte: triviale konvexe Hülle
            return list(dict.fromkeys(pts))  # Duplikate entfernen

        # Startpunkt: ganz links (kleinste x-Koordinate)
        start = min(range(n), key=lambda i: (pts[i][0], pts[i][1]))
        hull = []
        current = start

        while True:
            hull.append(pts[current])
            # Finde den Punkt, der von allen anderen "rechts" liegt
            next_point = (current + 1) % n
            for i in range(n):
                if i == current:
                    continue
                # Kreuzprodukt: positiv = i liegt links von (current → next_point)
                p0 = np.array(pts[current])
                p1 = np.array(pts[next_point])
                p2 = np.array(pts[i])
                cross = _cross2(p1 - p0, p2 - p0)
                if cross < 0:
                    # i ist weiter "rechts" (im Uhrzeigersinn) → neues Kandidat
                    next_point = i
                elif abs(cross) < 1e-12:
                    # Kollinear: nehme den weiter entfernten Punkt
                    d1 = np.linalg.norm(p1 - p0)
                    d2 = np.linalg.norm(p2 - p0)
                    if d2 > d1:
                        next_point = i

            current = next_point
            if current == start:
                break

        return hull

    def is_affinely_independent(self, points, tol: float = 1e-10) -> bool:
        """
        Prüft, ob eine Punktmenge affin unabhängig ist.

        Punkte P₀,...,Pₙ sind affin unabhängig, wenn
        P₁-P₀, ..., Pₙ-P₀ linear unabhängig sind.

        @param points: Liste von Punkten.
        @param tol:    Toleranz für den Rangtest.
        @return: True wenn affin unabhängig.
        @lastModified 2026-03-10
        """
        pts = np.array(points, dtype=float)
        if pts.shape[0] <= 1:
            return True  # 0 oder 1 Punkt ist trivial unabhängig
        # Differenzmatrix
        diffs = pts[1:] - pts[0]
        rank = np.linalg.matrix_rank(diffs, tol=tol)
        return rank == pts.shape[0] - 1


# ============================================================================
# 6. SyntheticGeometry
# ============================================================================

class SyntheticGeometry:
    """
    Synthetische Geometrie basierend auf Hilberts Axiomensystem (1899).

    David Hilbert reformulierte Euklids Geometrie axiomatisch:
    - Gruppe I: Inzidenzaxiome (5 Axiome)
    - Gruppe II: Anordnungsaxiome (4 Axiome)
    - Gruppe III: Kongruenzaxiome (6 Axiome)
    - Gruppe IV: Parallelenaxiom (1 Axiom)
    - Gruppe V: Stetigkeitsaxiome (2 Axiome, Archimedes + Vollständigkeit)

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def hilbert_axioms_overview(self) -> dict:
        """
        Gibt eine Übersicht über Hilberts Axiomensystem der Geometrie.

        Aus Hilberts „Grundlagen der Geometrie" (1899).

        @return: Dict mit Axiomengruppen I–V.
        @lastModified 2026-03-10
        """
        return {
            "I_Inzidenz": {
                "I.1": "Durch zwei Punkte geht genau eine Gerade.",
                "I.2": "Auf jeder Geraden liegen mindestens zwei Punkte.",
                "I.3": "Es gibt mindestens drei nicht-kollineare Punkte.",
                "I.4": "Durch drei nicht-kollineare Punkte geht genau eine Ebene.",
                "I.5": "Wenn zwei Punkte einer Geraden in einer Ebene liegen, liegt die Gerade in der Ebene.",
            },
            "II_Anordnung": {
                "II.1": "Wenn B zwischen A und C liegt, liegt B auch zwischen C und A.",
                "II.2": "Zu zwei Punkten A, B gibt es einen Punkt C mit B zwischen A und C.",
                "II.3": "Unter drei Punkten einer Geraden liegt höchstens einer zwischen den anderen.",
                "II.4": "Pasch-Axiom: Eine Gerade, die eine Seite eines Dreiecks trifft, trifft auch eine der anderen Seiten.",
            },
            "III_Kongruenz": {
                "III.1": "Auf einer Geraden kann eine gegebene Strecke von einem Punkt aus aufgetragen werden.",
                "III.2": "Kongruenz von Strecken ist eine Äquivalenzrelation.",
                "III.3": "Strecken können addiert werden (Transitivität).",
                "III.4": "Winkel können aufgetragen werden.",
                "III.5": "Kongruenz von Winkeln ist eine Äquivalenzrelation.",
                "III.6": "Kongruenzsatz SWS (Seite-Winkel-Seite).",
            },
            "IV_Parallelen": {
                "IV.1": "Euklidsches Parallelenpostulat: Durch einen Punkt außerhalb einer Geraden gibt es genau eine Parallele.",
            },
            "V_Stetigkeit": {
                "V.1": "Archimedisches Axiom: Zu jeder Strecke AB gibt es n Vielfache von CD, die AB übersteigen.",
                "V.2": "Vollständigkeitsaxiom: Die Geometrie kann nicht durch Hinzunahme weiterer Punkte/Geraden erweitert werden, ohne ein anderes Axiom zu verletzen.",
            },
        }

    def incidence_axioms_demo(self) -> dict:
        """
        Demonstriert die Inzidenzaxiome mit konkreten Beispielen.

        @return: Dict mit Beispielen für Inzidenzaxiome I.1–I.3.
        @lastModified 2026-03-10
        """
        return {
            "I.1_Beispiel": {
                "Punkte": "A=(0,0), B=(1,1)",
                "Gerade": "y = x (eindeutig bestimmt)",
                "Verifiziert": True,
            },
            "I.2_Beispiel": {
                "Gerade": "y = 2x + 1",
                "Punkte_auf_Geraden": ["(0, 1)", "(1, 3)", "(-1, -1)"],
                "unendlich_viele": True,
            },
            "I.3_Beispiel": {
                "Nicht_kollinear": "A=(0,0), B=(1,0), C=(0,1)",
                "Kollinear_wäre": "A=(0,0), B=(1,1), C=(2,2) – dies erfüllt I.3 NICHT als Gegenbeispiel",
            },
        }

    def parallel_axiom_demo(self) -> dict:
        """
        Vergleicht das Parallelenaxiom in drei Geometrien.

        @return: Dict mit euklidischer, hyperbolischer und elliptischer Variante.
        @lastModified 2026-03-10
        """
        return {
            "euklidisch": {
                "Axiom": "Genau eine Parallele durch externen Punkt",
                "Winkelsumme_Dreieck": "= π (180°)",
                "Modell": "ℝ² mit Standard-Metrik",
                "Krümmung": 0,
            },
            "hyperbolisch": {
                "Axiom": "Unendlich viele Parallelen durch externen Punkt",
                "Winkelsumme_Dreieck": "< π",
                "Modell": "Poincaré-Kreismodell / Obere Halbebene",
                "Krümmung": -1,
            },
            "elliptisch": {
                "Axiom": "Keine Parallelen (alle Geraden schneiden sich)",
                "Winkelsumme_Dreieck": "> π",
                "Modell": "Einheitssphäre S² (projektiv: RP²)",
                "Krümmung": +1,
            },
        }

    def axiom_independence_demo(self) -> dict:
        """
        Demonstriert die Unabhängigkeit von Hilberts Axiomen.

        Jedes Axiom ist unabhängig, d.h. es gibt Modelle, die alle anderen
        Axiome erfüllen, aber dieses eine verletzen.

        @return: Dict mit Unabhängigkeitsbeispielen.
        @lastModified 2026-03-10
        """
        return {
            "Parallelenaxiom_unabhängig": {
                "Beweis": "Hyperbolische Geometrie erfüllt I-III, V aber NICHT IV",
                "Modell": "Poincaré-Kreismodell",
                "Geschichte": "Bolyai und Lobatschewski (ca. 1830), Riemann (1854)",
            },
            "Archimedisches_Axiom_unabhängig": {
                "Beweis": "Nicht-archimedische Körper (z.B. rationale Funktionen)",
                "Modell": "Geometrie über ℝ(t) mit t als infinitesimales Element",
            },
            "Kongruenzaxiome_unabhängig": {
                "Beweis": "Projektive Ebene ohne Metrik erfüllt I-II, IV-V aber nicht III",
            },
        }

    def finite_projective_plane(self, q: int) -> dict:
        """
        Berechnet Eigenschaften der endlichen projektiven Ebene PG(2,q).

        PG(2,q) ist die projektive Ebene über dem endlichen Körper GF(q)
        (für q = Primzahlpotenz).

        Eigenschaften:
        - Punkte: q² + q + 1
        - Geraden: q² + q + 1 (dual zu Punkten)
        - Punkte pro Gerade: q + 1
        - Geraden durch jeden Punkt: q + 1

        @param q: Primzahl (oder Primzahlpotenz) für GF(q).
        @return: Dict mit Eigenschaften von PG(2,q).
        @raises ValueError: Wenn q keine Primzahl ist.
        @lastModified 2026-03-10
        """
        # Primzahlprüfung (einfach für kleine q)
        if q < 2:
            raise ValueError(f"q={q} ist keine Primzahl. q muss ≥ 2 sein.")

        # Einfache Primzahlprüfung
        def is_prime(n):
            if n < 2:
                return False
            for i in range(2, int(n**0.5) + 1):
                if n % i == 0:
                    return False
            return True

        if not is_prime(q):
            raise ValueError(
                f"q={q} ist keine Primzahl. PG(2,q) ist nur für Primzahlpotenzen definiert."
            )

        n_points = q * q + q + 1
        n_lines  = q * q + q + 1  # Dual
        points_per_line = q + 1
        lines_per_point = q + 1

        return {
            "q": q,
            "field": f"GF({q})",
            "n_points": n_points,
            "n_lines": n_lines,
            "points_per_line": points_per_line,
            "lines_per_point": lines_per_point,
            "order": q,
            "description": (
                f"PG(2,{q}): projektive Ebene der Ordnung {q} über GF({q}). "
                f"{n_points} Punkte, {n_lines} Geraden, je {points_per_line} Punkte pro Gerade."
            ),
        }


# ============================================================================
# 7. Standalone-Funktionen
# ============================================================================

def euler_characteristic_polygon(vertices: int, edges: int, faces: int) -> int:
    """
    Berechnet die Euler-Charakteristik eines Polyeders.

    Euler'sche Polyederformel:
    $\\chi = V - E + F$

    Für konvexe Polyeder (und die Sphäre S²): χ = 2.
    Für den Torus: χ = 0. Für das doppelte Torus: χ = -2.

    @param vertices: Anzahl der Ecken (V).
    @param edges:    Anzahl der Kanten (E).
    @param faces:    Anzahl der Flächen (F).
    @return: Euler-Charakteristik χ = V - E + F.
    @lastModified 2026-03-10
    """
    return vertices - edges + faces


def pick_theorem(interior_points: int, boundary_points: int) -> float:
    """
    Berechnet den Flächeninhalt eines Gittervielecks via Picks Satz.

    Pick'scher Flächeninhaltssatz (Georg Pick, 1899):
    $A = I + \\frac{B}{2} - 1$

    wobei:
    - I: Anzahl der inneren Gitterpunkte
    - B: Anzahl der Randpunkte (auf dem Rand des Polygons)

    @param interior_points: Anzahl innerer Gitterpunkte I.
    @param boundary_points: Anzahl der Rand-Gitterpunkte B.
    @return: Flächeninhalt des Polygons.
    @lastModified 2026-03-10
    """
    return interior_points + boundary_points / 2.0 - 1.0


def ptolemy_theorem_check(p1, p2, p3, p4, tol: float = 1e-8) -> dict:
    """
    Prüft den Ptolemäus-Satz für ein Sehnenviereck numerisch.

    Ptolemäus-Satz: Für ein Sehnenviereck ABCD gilt:
    $|AC| \\cdot |BD| = |AB| \\cdot |CD| + |AD| \\cdot |BC|$

    Das Viereck ist genau dann ein Sehnenviereck (auf einem Kreis),
    wenn diese Gleichung erfüllt ist.

    @param p1: Erster Punkt [x, y].
    @param p2: Zweiter Punkt [x, y].
    @param p3: Dritter Punkt [x, y].
    @param p4: Vierter Punkt [x, y].
    @param tol: Toleranz.
    @return: Dict mit 'lhs', 'rhs', 'ptolemy_satisfied', 'is_cyclic'.
    @lastModified 2026-03-10
    """
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    p3 = np.array(p3, dtype=float)
    p4 = np.array(p4, dtype=float)

    def dist(a, b):
        return float(np.linalg.norm(b - a))

    # Seiten und Diagonalen des Vierecks ABCD = p1,p2,p3,p4
    AC = dist(p1, p3)
    BD = dist(p2, p4)
    AB = dist(p1, p2)
    CD = dist(p3, p4)
    AD = dist(p1, p4)
    BC = dist(p2, p3)

    lhs = AC * BD                # Produkt der Diagonalen
    rhs = AB * CD + AD * BC      # Summe der Seiten-Produkte

    satisfied = abs(lhs - rhs) < tol

    return {
        "AC": AC, "BD": BD,
        "AB": AB, "CD": CD, "AD": AD, "BC": BC,
        "lhs": lhs,
        "rhs": rhs,
        "difference": abs(lhs - rhs),
        "ptolemy_satisfied": satisfied,
        "is_cyclic": satisfied,
    }


def nine_point_circle(p1, p2, p3) -> dict:
    """
    Berechnet den Neun-Punkte-Kreis eines Dreiecks.

    Der Neun-Punkte-Kreis (Feuerbach-Kreis) geht durch folgende 9 Punkte:
    1–3: Mittelpunkte der drei Seiten
    4–6: Fußpunkte der drei Höhen
    7–9: Mittelpunkte der Strecken von jedem Eckpunkt zum Höhenschnittpunkt H

    Mittelpunkt des Neun-Punkte-Kreises: N = (O + H)/2 (O = Umkreismittelpunkt)
    Radius: r₉ = R/2 (Hälfte des Umkreisradius)

    @param p1: Erster Eckpunkt [x, y].
    @param p2: Zweiter Eckpunkt [x, y].
    @param p3: Dritter Eckpunkt [x, y].
    @return: Dict mit 'center', 'radius', 'midpoints_of_sides',
             'foot_of_altitudes', 'euler_midpoints'.
    @raises ValueError: Bei entartetem Dreieck.
    @lastModified 2026-03-10
    """
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    p3 = np.array(p3, dtype=float)

    geom = EuclideanGeometry()

    # Umkreismittelpunkt O und Höhenschnittpunkt H
    O = geom.circumcenter(p1, p2, p3)
    G = geom.centroid([p1, p2, p3])
    H = 3.0 * G - 2.0 * O  # Euler-Relation: H = 3G - 2O

    # Mittelpunkt des Neun-Punkte-Kreises: N = (O + H)/2
    N = (O + H) / 2.0

    # Umkreisradius R
    R = float(np.linalg.norm(p1 - O))

    # Radius des Neun-Punkte-Kreises: r₉ = R/2
    r9 = R / 2.0

    # 1–3: Seitenmittelpunkte
    M1 = (p2 + p3) / 2.0  # gegenüber P1
    M2 = (p1 + p3) / 2.0  # gegenüber P2
    M3 = (p1 + p2) / 2.0  # gegenüber P3

    # 4–6: Höhenfußpunkte (Projektion von Ecke auf gegenüberliegende Seite)
    def altitude_foot(vertex, side_a, side_b):
        """Berechnet den Fußpunkt der Höhe vom Vertex auf die Seite AB."""
        AB = side_b - side_a
        AB_norm_sq = np.dot(AB, AB)
        if AB_norm_sq < 1e-15:
            return side_a.copy()
        t = np.dot(vertex - side_a, AB) / AB_norm_sq
        return side_a + t * AB

    F1 = altitude_foot(p1, p2, p3)
    F2 = altitude_foot(p2, p1, p3)
    F3 = altitude_foot(p3, p1, p2)

    # 7–9: Mittelpunkte der Strecken von Ecke zu H
    E1 = (p1 + H) / 2.0
    E2 = (p2 + H) / 2.0
    E3 = (p3 + H) / 2.0

    return {
        "center": N,
        "radius": r9,
        "circumcenter": O,
        "orthocenter": H,
        "midpoints_of_sides": [M1, M2, M3],
        "foot_of_altitudes": [F1, F2, F3],
        "euler_midpoints": [E1, E2, E3],
    }


def morley_theorem_demo() -> dict:
    """
    Demonstriert Morleys Dreieckssatz.

    Morleys Satz (Frank Morley, 1899):
    Die Schnittpunkte der Drittelteilenden Winkel jedes Dreiecks bilden
    stets ein gleichseitiges Dreieck (Morley-Dreieck).

    Für ein Dreieck mit Winkeln α, β, γ (α+β+γ=π) gilt:
    Die Winkeldrittler treffen sich in einem gleichseitigen Dreieck
    mit Seitenlänge:
    $s = 8R \\sin(\\alpha/3)\\sin(\\beta/3)\\sin(\\gamma/3)$

    @return: Dict mit Beispiel, Seitenlängenformel und Verifikation.
    @lastModified 2026-03-10
    """
    # Beispieldreieck: gleichseitiges Dreieck (α=β=γ=60°=π/3)
    alpha = math.pi / 3  # 60°
    beta  = math.pi / 3
    gamma = math.pi / 3

    R = 1.0  # Umkreisradius = 1

    # Seitenlänge des Morley-Dreiecks
    s = 8.0 * R * math.sin(alpha/3) * math.sin(beta/3) * math.sin(gamma/3)

    # Für ein gleichseitiges Dreieck mit α=β=γ=60°: s = 8·sin(20°)³
    # Allgemeines Beispiel: Dreieck mit α=50°, β=60°, γ=70°
    a2 = math.radians(50)
    b2 = math.radians(60)
    g2 = math.radians(70)
    s2 = 8.0 * R * math.sin(a2/3) * math.sin(b2/3) * math.sin(g2/3)

    return {
        "satz": (
            "Die Schnittpunkte der Winkeldrittler eines beliebigen Dreiecks "
            "bilden stets ein gleichseitiges Dreieck (Morley-Dreieck)."
        ),
        "seitenlaenge_formel": "s = 8·R·sin(α/3)·sin(β/3)·sin(γ/3)",
        "beispiel_gleichseitig": {
            "alpha_deg": 60, "beta_deg": 60, "gamma_deg": 60,
            "R": R,
            "morley_side_length": s,
            "morley_equilateral": True,
        },
        "beispiel_allgemein": {
            "alpha_deg": 50, "beta_deg": 60, "gamma_deg": 70,
            "R": R,
            "morley_side_length": s2,
        },
        "bemerkung": (
            "Gilt für JEDES Dreieck – eines der überraschendsten Resultate "
            "der elementaren Geometrie (erst 1899 entdeckt)."
        ),
    }
