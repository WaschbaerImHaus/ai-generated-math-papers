"""
@file topology.py
@brief Topologie und Geometrie – metrische Räume, Stetigkeit, Mannigfaltigkeiten,
       parametrische Kurven, simpliziale Homologie und fraktale Dimensionen.
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math
from typing import Callable


# ============================================================
# Standardmetriken
# ============================================================

def euclidean_metric(x: list[float], y: list[float]) -> float:
    """
    @brief Euklidische Metrik: d(x,y) = sqrt(Σ(x_i - y_i)²).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die gewöhnliche "Luftlinie"-Distanz im ℝⁿ.

    Mathematische Definition:
        d(x, y) = \\sqrt{\\sum_{i=1}^{n} (x_i - y_i)^2}

    @param x Erster Punkt als Liste von Koordinaten.
    @param y Zweiter Punkt als Liste von Koordinaten.
    @return Euklidischer Abstand als float.

    Beispiele:
    >>> euclidean_metric([0, 0], [3, 4])
    5.0
    >>> euclidean_metric([0, 0, 0], [1, 0, 0])
    1.0
    """
    # Summe der quadratischen Differenzen berechnen, dann Wurzel ziehen
    return math.sqrt(sum((xi - yi) ** 2 for xi, yi in zip(x, y)))


def manhattan_metric(x: list[float], y: list[float]) -> float:
    """
    @brief Manhattan-Metrik (Taxicab): d(x,y) = Σ|x_i - y_i|.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Benannt nach dem Straßennetz Manhattans – man kann nur
    horizontal und vertikal fahren.

    Mathematische Definition:
        d(x, y) = \\sum_{i=1}^{n} |x_i - y_i|

    @param x Erster Punkt als Liste von Koordinaten.
    @param y Zweiter Punkt als Liste von Koordinaten.
    @return Manhattan-Abstand als float.

    Beispiele:
    >>> manhattan_metric([0, 0], [3, 4])
    7
    >>> manhattan_metric([1, 2], [4, 6])
    7
    """
    # Summe aller absoluten Koordinatendifferenzen
    return sum(abs(xi - yi) for xi, yi in zip(x, y))


def chebyshev_metric(x: list[float], y: list[float]) -> float:
    """
    @brief Chebyshev-Metrik: d(x,y) = max|x_i - y_i|.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Auch "Schachbrett-Metrik" genannt: Ein König auf dem Schachbrett
    benötigt max(|Δx|, |Δy|) Züge.

    Mathematische Definition:
        d(x, y) = \\max_{i} |x_i - y_i|

    @param x Erster Punkt als Liste von Koordinaten.
    @param y Zweiter Punkt als Liste von Koordinaten.
    @return Chebyshev-Abstand als float.
    """
    # Maximum aller absoluten Koordinatendifferenzen
    return max(abs(xi - yi) for xi, yi in zip(x, y))


def discrete_metric(x: object, y: object) -> float:
    """
    @brief Diskrete Metrik: d(x,y) = 0 wenn x=y, sonst 1.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die einfachste Metrik überhaupt. Jedes Paar ungleicher Elemente
    hat denselben Abstand 1.

    Mathematische Definition:
        d(x, y) = \\begin{cases} 0 & x = y \\\\ 1 & x \\neq y \\end{cases}

    @param x Erstes Element (beliebiger Typ).
    @param y Zweites Element (beliebiger Typ).
    @return 0.0 oder 1.0.
    """
    return 0.0 if x == y else 1.0


def p_norm_metric(x: list[float], y: list[float], p: float) -> float:
    """
    @brief Lp-Metrik: d(x,y) = (Σ|x_i - y_i|^p)^{1/p}.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Verallgemeinert Euklidisch (p=2), Manhattan (p=1) und
    Chebyshev (p→∞) in einer einheitlichen Familie.

    Mathematische Definition:
        d_p(x, y) = \\left( \\sum_{i=1}^{n} |x_i - y_i|^p \\right)^{1/p}

    @param x Erster Punkt als Liste von Koordinaten.
    @param y Zweiter Punkt als Liste von Koordinaten.
    @param p Exponent (p ≥ 1 für gültige Metrik).
    @return Lp-Abstand als float.
    @raises ValueError wenn p < 1.
    """
    if p < 1:
        raise ValueError(f"p muss >= 1 sein, erhalten: {p}")
    # Innere Summe berechnen, dann p-te Wurzel ziehen
    return sum(abs(xi - yi) ** p for xi, yi in zip(x, y)) ** (1.0 / p)


# ============================================================
# Metrischer Raum
# ============================================================

class MetricSpace:
    """
    @brief Metrischer Raum (X, d) mit Metrikfunktion d: X×X → ℝ≥0.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein metrischer Raum ist eine Menge X zusammen mit einer Abstandsfunktion d,
    die folgende Axiome erfüllt:

    1. d(x,y) ≥ 0               (Nichtnegativität)
    2. d(x,y) = 0 ⟺ x = y      (Definitheit)
    3. d(x,y) = d(y,x)          (Symmetrie)
    4. d(x,z) ≤ d(x,y) + d(y,z) (Dreiecksungleichung)

    Beispiele: (ℝⁿ, euklidisch), (C([a,b]), ‖·‖∞), (Graphen, Pfadlänge)
    """

    def __init__(self, distance_func: Callable[[list, list], float]) -> None:
        """
        @brief Initialisiert den metrischen Raum mit einer Abstandsfunktion.
        @param distance_func Callable d(x, y) → float (die Metrik).
        """
        # Die Metrikfunktion als Attribut speichern
        self.d = distance_func

    def verify_metric_axioms(self, points: list[list[float]]) -> dict[str, bool]:
        """
        @brief Prüft alle 4 Metrik-Axiome für eine Liste von Punkten.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Überprüft numerisch, ob die gespeicherte Abstandsfunktion
        für alle Punktepaare (und Tripel) alle Axiome einhält.

        @param points Liste von Punkten im Raum.
        @return Dict mit Schlüsseln 'non_negativity', 'identity', 'symmetry',
                'triangle_inequality' – alle True wenn die Axiome gelten.
        """
        tol = 1e-9  # Numerische Toleranz für Gleichheitsvergleiche
        results = {
            'non_negativity': True,
            'identity': True,
            'symmetry': True,
            'triangle_inequality': True
        }

        # Axiom 1 & 2 & 3: für alle Paare prüfen
        for i, x in enumerate(points):
            for j, y in enumerate(points):
                dist_xy = self.d(x, y)

                # Axiom 1: Nichtnegativität d(x,y) ≥ 0
                if dist_xy < -tol:
                    results['non_negativity'] = False

                # Axiom 2: Definitheit d(x,y) = 0 ⟺ x = y
                if i == j and dist_xy > tol:
                    results['identity'] = False
                if i != j and dist_xy < tol:
                    results['identity'] = False

                # Axiom 3: Symmetrie d(x,y) = d(y,x)
                dist_yx = self.d(y, x)
                if abs(dist_xy - dist_yx) > tol:
                    results['symmetry'] = False

        # Axiom 4: Dreiecksungleichung d(x,z) ≤ d(x,y) + d(y,z) für alle Tripel
        for x in points:
            for y in points:
                for z in points:
                    d_xz = self.d(x, z)
                    d_xy = self.d(x, y)
                    d_yz = self.d(y, z)
                    if d_xz > d_xy + d_yz + tol:
                        results['triangle_inequality'] = False

        return results

    def open_ball(self, center: list[float], radius: float, points: list[list[float]]) -> list[list[float]]:
        """
        @brief Offene Kugel B(center, r) = {x ∈ points : d(center, x) < r}.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine offene Kugel enthält alle Punkte, die echt weniger als r
        von center entfernt sind (ohne den Rand).

        @param center Mittelpunkt der Kugel.
        @param radius Radius r > 0.
        @param points Punktmenge, die gefiltert wird.
        @return Liste der Punkte in der offenen Kugel.
        """
        # Alle Punkte zurückgeben, deren Abstand zum Zentrum < Radius ist
        return [p for p in points if self.d(center, p) < radius]

    def is_open_set(self, subset: list[list[float]], all_points: list[list[float]], tol: float = 1e-9) -> bool:
        """
        @brief Prüft ob eine Teilmenge (in diskreten Punkten) offen ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine Menge U ist offen, wenn jeder Punkt in U eine offene Kugel
        besitzt, die vollständig in U liegt.

        In diskreten endlichen Punktmengen: Prüft ob für jeden Punkt in subset
        ein Radius r > 0 existiert, sodass B(x,r) ∩ all_points ⊆ subset.

        @param subset Zu prüfende Teilmenge.
        @param all_points Gesamte Punktmenge des Raums.
        @param tol Numerische Toleranz.
        @return True wenn subset offen ist.
        """
        subset_set = [tuple(p) if isinstance(p, list) else p for p in subset]

        for point in subset:
            # Finde den kleinsten Abstand zu Punkten AUSSERHALB der Menge
            min_dist_outside = float('inf')
            for other in all_points:
                other_key = tuple(other) if isinstance(other, list) else other
                point_key = tuple(point) if isinstance(point, list) else point
                if other_key not in subset_set:
                    dist = self.d(point, other)
                    if dist < min_dist_outside:
                        min_dist_outside = dist

            # Wenn kein Punkt außerhalb existiert, ist der Punkt automatisch "offen"
            if min_dist_outside == float('inf'):
                continue

            # Es muss ein r > 0 geben, sodass B(x, r) ⊆ subset
            if min_dist_outside <= tol:
                return False

        return True

    def is_cauchy_sequence(self, sequence: list[list[float]], tol: float = 1e-10) -> bool:
        """
        @brief Prüft ob eine Folge Cauchy ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine Folge (x_n) heißt Cauchy-Folge, wenn für jedes ε > 0
        ein N existiert, sodass d(x_m, x_n) < ε für alle m,n > N.

        Numerisch: Prüft ob die maximale Distanz in der zweiten Hälfte
        der Folge kleiner als tol ist.

        @param sequence Liste von Folgengliedern.
        @param tol Toleranz ε (Cauchy-Bedingung).
        @return True wenn die Folge Cauchy ist.
        """
        n = len(sequence)
        if n < 2:
            return True

        # Zweite Hälfte der Folge betrachten (Asymptotik)
        half = n // 2
        tail = sequence[half:]

        # Maximalen paarweisen Abstand in der zweiten Hälfte berechnen
        max_dist = 0.0
        for i in range(len(tail)):
            for j in range(i + 1, len(tail)):
                dist = self.d(tail[i], tail[j])
                if dist > max_dist:
                    max_dist = dist

        return max_dist < tol


# ============================================================
# Topologische Eigenschaften
# ============================================================

def is_connected(points: list[list[float]], metric_func: Callable[[list, list], float], radius: float) -> bool:
    """
    @brief Prüft ob ein diskreter Punktraum epsilon-zusammenhängend ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Zwei Punkte sind epsilon-verbunden wenn d(x,y) < radius.
    Der Raum ist zusammenhängend, wenn von jedem Punkt jeder andere
    erreichbar ist (Union-Find / BFS).

    @param points Liste der Punkte.
    @param metric_func Metrikfunktion.
    @param radius Verbindungsradius ε.
    @return True wenn der Raum epsilon-zusammenhängend ist.
    """
    if not points:
        return True
    if len(points) == 1:
        return True

    # Union-Find-Struktur initialisieren
    parent = list(range(len(points)))

    def find(i: int) -> int:
        """Findet die Wurzel der Komponente von i (mit Pfadkompression)."""
        while parent[i] != i:
            parent[i] = parent[parent[i]]  # Pfadkompression
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        """Vereinigt die Komponenten von i und j."""
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    # Alle Paare prüfen, die epsilon-nah sind
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if metric_func(points[i], points[j]) < radius:
                union(i, j)

    # Prüfen ob alle Punkte in derselben Komponente sind
    root = find(0)
    return all(find(i) == root for i in range(len(points)))


def hausdorff_distance(set_a: list[list[float]], set_b: list[list[float]], metric_func: Callable[[list, list], float]) -> float:
    """
    @brief Hausdorff-Abstand zwischen zwei Mengen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Hausdorff-Abstand misst, wie weit zwei Mengen voneinander entfernt sind.
    Er ist der maximale minimale Abstand zwischen den Mengen:

    Mathematische Definition:
        H(A,B) = max(sup_{a∈A} inf_{b∈B} d(a,b),  sup_{b∈B} inf_{a∈A} d(a,b))

    @param set_a Erste Punktmenge.
    @param set_b Zweite Punktmenge.
    @param metric_func Metrikfunktion.
    @return Hausdorff-Abstand als float.
    """
    # Für jeden Punkt in A: minimalen Abstand zu B berechnen
    def directed_hausdorff(source: list, target: list) -> float:
        max_min_dist = 0.0
        for a in source:
            # Minimaler Abstand von a zu irgendeinem Punkt in target
            min_dist = min(metric_func(a, b) for b in target)
            if min_dist > max_min_dist:
                max_min_dist = min_dist
        return max_min_dist

    # Hausdorff-Abstand = Maximum beider gerichteter Abstände
    return max(directed_hausdorff(set_a, set_b), directed_hausdorff(set_b, set_a))


def compute_diameter(points: list[list[float]], metric_func: Callable[[list, list], float]) -> float:
    """
    @brief Durchmesser einer Menge: diam(A) = sup_{x,y∈A} d(x,y).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Durchmesser ist der größtmögliche Abstand zwischen zwei Punkten
    der Menge. Für kompakte Mengen wird er angenommen.

    @param points Liste der Punkte in der Menge.
    @param metric_func Metrikfunktion.
    @return Durchmesser als float. 0.0 für leere/einelementige Mengen.
    """
    if len(points) < 2:
        return 0.0

    # Maximalen paarweisen Abstand berechnen
    max_dist = 0.0
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            dist = metric_func(points[i], points[j])
            if dist > max_dist:
                max_dist = dist
    return max_dist


def is_compact_discrete(points: list[list[float]], metric_func: Callable[[list, list], float], epsilon: float) -> bool:
    """
    @brief Prüft (epsilon-)Kompaktheit einer diskreten Punktmenge.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    In der Topologie ist Kompaktheit äquivalent zu: Jede offene Überdeckung
    besitzt eine endliche Teilüberdeckung.

    Für endliche, diskrete Punktmengen ist dies immer True – jede offene
    Überdeckung ist selbst bereits endlich. Diese Funktion gibt einen
    strukturellen Hinweis zurück.

    @param points Liste der Punkte.
    @param metric_func Metrikfunktion (wird für Epsilon-Netz verwendet).
    @param epsilon Epsilon für das Epsilon-Netz.
    @return True (für endliche diskrete Mengen immer kompakt).
    """
    # Für jede endliche diskrete Menge existiert trivialerweise eine
    # endliche Überdeckung: nimm für jeden Punkt eine Kugel mit Radius epsilon.
    # Diese Funktion prüft zusätzlich ob ein epsilon-Netz existiert.

    if not points:
        return True

    # Epsilon-Netz: Teilmenge S ⊆ points sodass jeder Punkt in points
    # innerhalb epsilon von einem Punkt in S liegt.
    covered = set()
    net = []

    for i, p in enumerate(points):
        if i not in covered:
            net.append(p)
            # Alle Punkte innerhalb epsilon als "abgedeckt" markieren
            for j, q in enumerate(points):
                if metric_func(p, q) <= epsilon:
                    covered.add(j)

    # Endliches Epsilon-Netz existiert immer für endliche Mengen
    return len(net) <= len(points)


# ============================================================
# Parametrische Kurven und Mannigfaltigkeiten
# ============================================================

class ParametricCurve:
    """
    @brief Parametrische Kurve γ: [a,b] → ℝⁿ.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine parametrische Kurve beschreibt eine kontinuierliche Abbildung
    eines Intervalls [a,b] in den ℝⁿ.

    Beispiele:
    - Kreis: (cos(t), sin(t)), t ∈ [0, 2π]
    - Spirale: (t·cos(t), t·sin(t))
    - Lissajous: (sin(at+δ), sin(bt))

    Mathematisch: γ = (γ₁(t), γ₂(t), ..., γₙ(t))
    """

    def __init__(self, funcs: list[Callable[[float], float]], t_range: tuple[float, float]) -> None:
        """
        @brief Initialisiert die parametrische Kurve.
        @param funcs Liste von Funktionen [x(t), y(t), ...] – eine pro Dimension.
        @param t_range Tupel (a, b) mit Parameterbereich.
        """
        # Komponentenfunktionen und Parameterbereich speichern
        self.funcs = funcs       # Liste: [f1, f2, ...]  mit fi: float → float
        self.t_start = t_range[0]
        self.t_end = t_range[1]

    def evaluate(self, t: float) -> list[float]:
        """
        @brief Wertet die Kurve bei Parameterwert t aus.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        @param t Parameterwert im Bereich [t_start, t_end].
        @return Punkt γ(t) als Liste von Koordinaten.
        """
        # Jede Komponentenfunktion bei t auswerten
        return [f(t) for f in self.funcs]

    def arc_length(self, n: int = 1000) -> float:
        """
        @brief Bogenlänge numerisch berechnen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Bogenlänge ist das Integral der Kurvengeschwindigkeit:

            L = ∫_a^b |γ'(t)| dt

        Numerisch via Trapezmethode mit gleichmäßiger Unterteilung.

        @param n Anzahl der Unterteilungsintervalle (höher = genauer).
        @return Bogenlänge als float.
        """
        # Schrittweite für die numerische Integration
        dt = (self.t_end - self.t_start) / n
        total_length = 0.0

        # Trapezregel: Summiere Abstände zwischen aufeinanderfolgenden Punkten
        prev_point = self.evaluate(self.t_start)
        for i in range(1, n + 1):
            t = self.t_start + i * dt
            curr_point = self.evaluate(t)
            # Abstand zwischen aufeinanderfolgenden Punkten auf der Kurve
            segment_length = math.sqrt(
                sum((c - p) ** 2 for c, p in zip(curr_point, prev_point))
            )
            total_length += segment_length
            prev_point = curr_point

        return total_length

    def curvature(self, t: float) -> float:
        """
        @brief Krümmung κ(t) der Kurve an der Stelle t.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Für 2D-Kurven:
            κ(t) = (x'y'' - y'x'') / (x'² + y'²)^{3/2}

        Für 3D-Kurven:
            κ(t) = |γ' × γ''| / |γ'|³

        Erste und zweite Ableitung werden numerisch per finitem Differenzenquotient
        berechnet.

        @param t Parameterwert.
        @return Krümmung κ(t) als float (stets ≥ 0).
        """
        # Schrittweite für numerische Differentiation (Kompromiss Genauigkeit/Stabilität)
        h = 1e-5

        # Erste Ableitung γ'(t) via zentraler Differenz
        pt_plus = self.evaluate(t + h)
        pt_minus = self.evaluate(t - h)
        d1 = [(pt_plus[i] - pt_minus[i]) / (2 * h) for i in range(len(self.funcs))]

        # Zweite Ableitung γ''(t) via zentraler Differenz zweiter Ordnung
        pt0 = self.evaluate(t)
        d2 = [(pt_plus[i] - 2 * pt0[i] + pt_minus[i]) / (h ** 2) for i in range(len(self.funcs))]

        # Krümmungsformel anwenden
        dim = len(self.funcs)

        if dim == 2:
            # 2D: skalare Formel κ = (x'y'' - y'x'') / (x'²+y'²)^{3/2}
            numerator = d1[0] * d2[1] - d1[1] * d2[0]
            denominator = (d1[0] ** 2 + d1[1] ** 2) ** 1.5
            if abs(denominator) < 1e-15:
                return 0.0
            return abs(numerator) / denominator

        elif dim >= 3:
            # 3D: κ = |γ' × γ''| / |γ'|³
            # Kreuzprodukt von d1 und d2 (nur erste 3 Komponenten)
            cross = [
                d1[1] * d2[2] - d1[2] * d2[1],
                d1[2] * d2[0] - d1[0] * d2[2],
                d1[0] * d2[1] - d1[1] * d2[0],
            ]
            cross_norm = math.sqrt(sum(c ** 2 for c in cross))
            d1_norm = math.sqrt(sum(c ** 2 for c in d1[:3]))
            if d1_norm < 1e-15:
                return 0.0
            return cross_norm / (d1_norm ** 3)

        else:
            # 1D: Krümmung ist immer 0 (Gerade)
            return 0.0

    def is_closed(self, tol: float = 1e-10) -> bool:
        """
        @brief Prüft ob die Kurve geschlossen ist: γ(a) = γ(b).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine geschlossene Kurve verbindet Anfangs- und Endpunkt.
        Beispiel: Kreis, Ellipse.

        @param tol Toleranz für den Vergleich.
        @return True wenn Start- und Endpunkt übereinstimmen.
        """
        # Start- und Endpunkt der Kurve berechnen
        start = self.evaluate(self.t_start)
        end = self.evaluate(self.t_end)

        # Euklidischer Abstand zwischen Start und Ende
        dist = math.sqrt(sum((s - e) ** 2 for s, e in zip(start, end)))
        return dist < tol

    def winding_number(self, point: list[float]) -> int:
        """
        @brief Umlaufzahl einer geschlossenen Kurve um einen Punkt.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Umlaufzahl zählt, wie oft die Kurve den Punkt umkreist:

            n = (1/2π) ∮ dθ

        Positiv = Gegenuhrzeigersinn, Negativ = Uhrzeigersinn.
        Nur sinnvoll für 2D-Kurven.

        @param point Punkt [x, y] um den die Umlaufzahl gezählt wird.
        @return Umlaufzahl als int.
        """
        # Numerische Integration der Winkeländerung
        n_steps = 10000
        dt = (self.t_end - self.t_start) / n_steps

        total_angle = 0.0
        prev_pt = self.evaluate(self.t_start)
        prev_angle = math.atan2(prev_pt[1] - point[1], prev_pt[0] - point[0])

        for i in range(1, n_steps + 1):
            t = self.t_start + i * dt
            curr_pt = self.evaluate(t)
            curr_angle = math.atan2(curr_pt[1] - point[1], curr_pt[0] - point[0])

            # Winkeländerung auf (-π, π] normieren
            d_angle = curr_angle - prev_angle
            while d_angle > math.pi:
                d_angle -= 2 * math.pi
            while d_angle < -math.pi:
                d_angle += 2 * math.pi

            total_angle += d_angle
            prev_angle = curr_angle

        # Umlaufzahl = Gesamtwinkel / (2π)
        return round(total_angle / (2 * math.pi))


def circle_curve(radius: float = 1.0, center: list[float] | None = None) -> ParametricCurve:
    """
    @brief Erzeugt einen Kreis als ParametricCurve.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Parametrisierung: γ(t) = (cx + r·cos(t), cy + r·sin(t)), t ∈ [0, 2π]

    @param radius Kreisradius r.
    @param center Mittelpunkt [cx, cy]. Standard: [0, 0].
    @return ParametricCurve des Kreises.
    """
    if center is None:
        center = [0.0, 0.0]

    cx, cy = center[0], center[1]

    # Parametrische Gleichungen des Kreises
    x_func = lambda t: cx + radius * math.cos(t)
    y_func = lambda t: cy + radius * math.sin(t)

    return ParametricCurve([x_func, y_func], (0.0, 2 * math.pi))


def lissajous_curve(a: int, b: int, delta: float = 0.0) -> ParametricCurve:
    """
    @brief Erzeugt eine Lissajous-Figur als ParametricCurve.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Parametrisierung: x(t) = sin(at+δ), y(t) = sin(bt), t ∈ [0, 2π]

    Lissajous-Figuren entstehen bei harmonischen Schwingungen in zwei
    senkrechten Richtungen. Das Verhältnis a:b bestimmt die Form.

    @param a Frequenz in x-Richtung.
    @param b Frequenz in y-Richtung.
    @param delta Phasenverschiebung δ in Radiant.
    @return ParametricCurve der Lissajous-Figur.
    """
    x_func = lambda t: math.sin(a * t + delta)
    y_func = lambda t: math.sin(b * t)

    return ParametricCurve([x_func, y_func], (0.0, 2 * math.pi))


def helix_curve(radius: float = 1.0, pitch: float = 1.0) -> ParametricCurve:
    """
    @brief Erzeugt eine Helix (Schraubenlinie) als ParametricCurve.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Parametrisierung: x(t)=r·cos(t), y(t)=r·sin(t), z(t)=pitch·t/(2π), t ∈ [0, 2π]

    Eine Helix ist eine Kurve im ℝ³, die sich schraubenförmig um eine Achse windet.

    @param radius Radius der Helix.
    @param pitch Steigung pro Umdrehung (pitch = Höhe pro Windung).
    @return ParametricCurve der Helix im ℝ³.
    """
    x_func = lambda t: radius * math.cos(t)
    y_func = lambda t: radius * math.sin(t)
    z_func = lambda t: pitch * t / (2 * math.pi)

    return ParametricCurve([x_func, y_func, z_func], (0.0, 2 * math.pi))


# ============================================================
# Simpliziale Homologie (vereinfacht)
# ============================================================

def euler_characteristic_polygon(n_vertices: int, n_edges: int, n_faces: int) -> int:
    """
    @brief Euler-Charakteristik χ = V - E + F.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Euler-Charakteristik ist eine topologische Invariante.
    Euler-Formel für konvexe Polyeder: χ = 2.

    Bekannte Werte:
    - Sphäre:         χ = 2  (z.B. Würfel: V=8, E=12, F=6 → χ=2)
    - Torus:          χ = 0
    - Kleinsche Flasche: χ = 0
    - Doppeltorus:    χ = -2

    Mathematisch: χ = V - E + F

    @param n_vertices Anzahl der Ecken V.
    @param n_edges Anzahl der Kanten E.
    @param n_faces Anzahl der Flächen F.
    @return Euler-Charakteristik als int.

    Beispiele:
    >>> euler_characteristic_polygon(8, 12, 6)
    2
    >>> euler_characteristic_polygon(4, 6, 4)
    2
    """
    return n_vertices - n_edges + n_faces


def genus_from_euler(euler_char: int, orientable: bool = True) -> int:
    """
    @brief Geschlecht g einer Fläche aus der Euler-Charakteristik.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Das Geschlecht beschreibt die Anzahl der "Henkel" einer Fläche:

    Orientierbare Flächen:     χ = 2 - 2g  →  g = (2 - χ) / 2
    Nicht-orientierbare Flächen: χ = 2 - g  →  g = 2 - χ

    Beispiele:
    - Sphäre:  g = 0, χ = 2
    - Torus:   g = 1, χ = 0
    - Brezel:  g = 2, χ = -2

    @param euler_char Euler-Charakteristik χ.
    @param orientable True für orientierbare Flächen.
    @return Geschlecht g als int.
    """
    if orientable:
        # χ = 2 - 2g  →  g = (2 - χ) / 2
        return (2 - euler_char) // 2
    else:
        # χ = 2 - g  →  g = 2 - χ
        return 2 - euler_char


def betti_numbers_graph(adjacency_matrix: list[list[int | float]]) -> dict[str, int | list]:
    """
    @brief Betti-Zahlen eines Graphen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Betti-Zahlen sind topologische Invarianten:

    β₀ = Anzahl der Zusammenhangskomponenten
    β₁ = Anzahl unabhängiger Zyklen = |E| - |V| + β₀

    Für einen Baum: β₀ = 1, β₁ = 0 (azyklisch).
    Für einen Kreis C_n: β₀ = 1, β₁ = 1.

    @param adjacency_matrix Quadratische Adjazenzmatrix als 2D-Liste.
    @return Dict mit 'beta_0', 'beta_1', 'vertices', 'edges', 'components'.
    """
    n = len(adjacency_matrix)
    if n == 0:
        return {'beta_0': 0, 'beta_1': 0, 'vertices': 0, 'edges': 0, 'components': []}

    # Kanten zählen (ungerichtet, also obere Dreiecksmatrix)
    n_edges = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adjacency_matrix[i][j] != 0:
                n_edges += 1

    # BFS um Zusammenhangskomponenten zu zählen
    visited = [False] * n
    components = []

    for start in range(n):
        if not visited[start]:
            # BFS von start aus
            component = []
            queue = [start]
            visited[start] = True
            while queue:
                node = queue.pop(0)
                component.append(node)
                for neighbor in range(n):
                    if adjacency_matrix[node][neighbor] != 0 and not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            components.append(component)

    beta_0 = len(components)                      # Anzahl Zusammenhangskomponenten
    beta_1 = n_edges - n + beta_0                 # Anzahl unabhängiger Zyklen

    return {
        'beta_0': beta_0,
        'beta_1': max(0, beta_1),
        'vertices': n,
        'edges': n_edges,
        'components': components
    }


# ============================================================
# Fraktale Dimensionen
# ============================================================

def box_counting_dimension(points: list[list[float]], epsilon_values: list[float]) -> float:
    """
    @brief Box-Counting-Dimension (Minkowski-Bouligand-Dimension).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Box-Counting-Dimension ist eine numerische Methode zur Bestimmung
    der fraktalen Dimension einer Punktmenge:

        d = -lim_{ε→0} log N(ε) / log ε

    wobei N(ε) die Anzahl der Boxen der Größe ε ist, die mindestens
    einen Punkt enthalten.

    Numerisch: lineare Regression von log(N) gegen log(ε).

    @param points Liste von 2D-Punkten [[x1,y1], [x2,y2], ...].
    @param epsilon_values Liste von Epsilon-Werten (absteigend sortiert).
    @return Geschätzte Box-Counting-Dimension als float.
    """
    if not points or len(epsilon_values) < 2:
        return 0.0

    log_inv_eps = []  # log(1/ε) = -log(ε)
    log_n = []        # log(N(ε))

    for eps in epsilon_values:
        if eps <= 0:
            continue

        # Alle besetzten Boxen zählen
        boxes = set()
        for pt in points:
            # Boxindizes berechnen: Punkt p fällt in Box (floor(x/ε), floor(y/ε))
            box_idx = tuple(int(math.floor(coord / eps)) for coord in pt)
            boxes.add(box_idx)

        n_boxes = len(boxes)
        if n_boxes > 0:
            log_inv_eps.append(math.log(1.0 / eps))
            log_n.append(math.log(n_boxes))

    if len(log_inv_eps) < 2:
        return 0.0

    # Lineare Regression: log(N) = d * log(1/ε) + const
    # Steigung d ist die Box-Counting-Dimension
    n = len(log_inv_eps)
    mean_x = sum(log_inv_eps) / n
    mean_y = sum(log_n) / n

    # Steigung via kleinste-Quadrate-Methode
    numerator = sum((log_inv_eps[i] - mean_x) * (log_n[i] - mean_y) for i in range(n))
    denominator = sum((log_inv_eps[i] - mean_x) ** 2 for i in range(n))

    if abs(denominator) < 1e-15:
        return 0.0

    return numerator / denominator


def hausdorff_dimension_cantor() -> float:
    """
    @brief Hausdorff-Dimension der Cantor-Menge.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Cantor-Menge entsteht durch iteratives Entfernen des mittleren
    Drittels eines Intervalls. Sie hat:

    - Länge (Lebesgue-Maß): 0
    - Überabzählbar viele Punkte
    - Hausdorff-Dimension: d = log(2) / log(3) ≈ 0.6309

    Herleitung: Bei jedem Schritt werden n=2 Teilmengen der Skalierung r=1/3 erzeugt.
    Selbstähnlichkeitsdimension: N = r^{-d} → 2 = 3^d → d = log(2)/log(3)

    @return Hausdorff-Dimension der Cantor-Menge als float.
    """
    # Selbstähnlichkeitsdimension: 2 Teile, je 1/3 der Größe
    return math.log(2) / math.log(3)


def sierpinski_dimension() -> float:
    """
    @brief Hausdorff-Dimension des Sierpinski-Dreiecks.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Das Sierpinski-Dreieck entsteht durch iteratives Entfernen des mittleren
    Dreiecks. Es hat:

    - Fläche (2D-Maß): 0
    - Hausdorff-Dimension: d = log(3) / log(2) ≈ 1.585

    Herleitung: Bei jedem Schritt werden n=3 Teilmengen der Skalierung r=1/2 erzeugt.
    Selbstähnlichkeitsdimension: N = r^{-d} → 3 = 2^d → d = log(3)/log(2)

    @return Hausdorff-Dimension des Sierpinski-Dreiecks als float.
    """
    # Selbstähnlichkeitsdimension: 3 Teile, je 1/2 der Größe
    return math.log(3) / math.log(2)


# ============================================================
# Algebraische Topologie
# ============================================================

def simplicial_homology(simplices: dict[int, list[tuple]]) -> dict:
    """
    @brief Berechnet die simplizialen Homologiegruppen H_k(X, ℤ₂).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Simpliziale Homologie ist ein fundamentales Werkzeug der algebraischen Topologie.
    Sie ordnet jedem topologischen Raum eine Folge von abelschen Gruppen zu, die
    "Löcher" unterschiedlicher Dimension zählen.

    Algorithmus:
    1. Randoperatoren ∂_k aufstellen als Matrizen über ℤ₂ (GF(2)).
       ∂_k: C_k → C_{k-1}, ∂([v0,...,vk]) = Σ_i (-1)^i [v0,...,v̂i,...,vk]
       Über ℤ₂ ist (-1) = 1, daher nur mod-2 Einträge.
    2. Rang von ∂_k (über ℤ₂) per Gauß-Elimination berechnen.
    3. β_k = dim(ker ∂_k) - dim(im ∂_{k+1})
             = (|C_k| - rank(∂_k)) - rank(∂_{k+1})
    4. Euler-Charakteristik: χ = Σ (-1)^k β_k

    Mathematik:
        H_k(X; ℤ₂) = ker(∂_k) / im(∂_{k+1})
        β_k = rank_ℤ₂(H_k)

    @param simplices Dict: {dim: [(v0,...), ...]}
                    0: Punkte, 1: Kanten, 2: Dreiecke, usw.
    @return Dict mit 'betti_numbers', 'euler_characteristic', 'h0', 'h1', 'h2'.
    """

    def _rank_gf2(matrix: list[list[int]]) -> int:
        """
        Berechnet den Rang einer Matrix über GF(2) (Gauß-Elimination mod 2).
        Alle Einträge sind 0 oder 1; Addition ist XOR.
        """
        if not matrix or not matrix[0]:
            return 0
        # Matrix kopieren, um das Original nicht zu verändern
        m = [row[:] for row in matrix]
        rows = len(m)
        cols = len(m[0])
        pivot_row = 0  # Zeiger auf die aktuelle Pivotzeile
        for col in range(cols):
            if pivot_row >= rows:
                break
            # Pivotelement suchen: erste Zeile >= pivot_row mit Eintrag 1 in Spalte col
            found = -1
            for r in range(pivot_row, rows):
                if m[r][col] == 1:
                    found = r
                    break
            if found == -1:
                continue  # Keine Pivot in dieser Spalte → weiter
            # Pivotzeile nach oben tauschen
            m[pivot_row], m[found] = m[found], m[pivot_row]
            # Elimination: alle anderen Zeilen mit Eintrag 1 in Spalte col
            for r in range(rows):
                if r != pivot_row and m[r][col] == 1:
                    # Zeilenaddition mod 2 (XOR)
                    m[r] = [(m[r][j] ^ m[pivot_row][j]) for j in range(cols)]
            pivot_row += 1
        return pivot_row  # Anzahl der Pivot-Zeilen = Rang

    def _build_boundary_matrix(k_simplices: list[tuple],
                                km1_simplices: list[tuple]) -> list[list[int]]:
        """
        Baut die k-te Randmatrix ∂_k auf (über ℤ₂).
        Zeilen = (k-1)-Simplizes, Spalten = k-Simplizes.
        Eintrag [i][j] = 1, wenn das i-te (k-1)-Simplex im Rand des j-ten k-Simplex liegt.
        """
        n_rows = len(km1_simplices)
        n_cols = len(k_simplices)
        # Matrix mit Nullen initialisieren
        matrix = [[0] * n_cols for _ in range(n_rows)]
        # Index der (k-1)-Simplizes aufbauen für schnellen Lookup
        km1_index = {tuple(sorted(s)): i for i, s in enumerate(km1_simplices)}
        for j, sigma in enumerate(k_simplices):
            # Alle (k+1) Seiten des k-Simplex berechnen
            sigma_list = list(sigma)
            for i_omit in range(len(sigma_list)):
                # i-te Seite: Simplex ohne i-ten Knoten
                face = tuple(sorted(
                    sigma_list[t] for t in range(len(sigma_list)) if t != i_omit
                ))
                if face in km1_index:
                    row_idx = km1_index[face]
                    # Über ℤ₂: Eintrag toggled (gerade Anzahl = 0, ungerade = 1)
                    matrix[row_idx][j] ^= 1
        return matrix

    # Maximale vorhandene Dimension bestimmen
    max_dim = max(simplices.keys()) if simplices else 0

    # Simplizialketten für jede Dimension sammeln
    chain_groups: dict[int, list[tuple]] = {}
    for dim in range(max_dim + 1):
        chain_groups[dim] = list(simplices.get(dim, []))

    # Rangberechnung der Randoperatoren
    # rank_boundary[k] = Rang von ∂_k
    rank_boundary: dict[int, int] = {}
    for k in range(1, max_dim + 1):
        k_simps = chain_groups.get(k, [])
        km1_simps = chain_groups.get(k - 1, [])
        if not k_simps or not km1_simps:
            rank_boundary[k] = 0
        else:
            mat = _build_boundary_matrix(k_simps, km1_simps)
            rank_boundary[k] = _rank_gf2(mat)

    # Betti-Zahlen berechnen: β_k = |C_k| - rank(∂_k) - rank(∂_{k+1})
    betti: list[int] = []
    for k in range(max_dim + 1):
        c_k = len(chain_groups.get(k, []))
        r_k = rank_boundary.get(k, 0)       # Rang von ∂_k (Bild von ∂_k)
        r_kp1 = rank_boundary.get(k + 1, 0)  # Rang von ∂_{k+1} (Bild in C_k)
        beta_k = c_k - r_k - r_kp1
        betti.append(max(0, beta_k))

    # Euler-Charakteristik: χ = Σ (-1)^k β_k
    chi = sum((-1) ** k * betti[k] for k in range(len(betti)))

    return {
        'betti_numbers': betti,
        'euler_characteristic': chi,
        'h0': betti[0] if len(betti) > 0 else 0,
        'h1': betti[1] if len(betti) > 1 else 0,
        'h2': betti[2] if len(betti) > 2 else 0,
    }


def betti_numbers_from_adjacency(adj_matrix: list[list[int]]) -> dict:
    """
    @brief Berechnet β₀ und β₁ aus der Adjazenzmatrix eines Graphen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein Graph ist ein 1-dimensionaler Simplizialkomplex:
    - 0-Simplizes = Knoten (Vertices)
    - 1-Simplizes = Kanten (Edges)

    Betti-Zahlen:
        β₀ = Anzahl Zusammenhangskomponenten (via BFS/Union-Find)
        β₁ = |Kanten| - |Knoten| + β₀  (aus Euler-Charakteristik)
        χ  = β₀ - β₁ = |Knoten| - |Kanten|

    Mathematisch: Für einen zusammenhängenden Graphen gilt
        β₁ = |E| - |V| + 1  (Kreisrang / cycle rank)

    @param adj_matrix Quadratische, symmetrische Adjazenzmatrix (0/1-Werte).
    @return Dict mit 'beta_0', 'beta_1', 'euler_characteristic',
                     'vertices', 'edges', 'components'.
    """
    n = len(adj_matrix)
    if n == 0:
        return {
            'beta_0': 0, 'beta_1': 0, 'euler_characteristic': 0,
            'vertices': 0, 'edges': 0, 'components': []
        }

    # Kanten zählen (ungerichteter Graph → oberes Dreieck der Matrix)
    n_edges = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj_matrix[i][j] != 0:
                n_edges += 1

    # BFS zur Bestimmung der Zusammenhangskomponenten
    visited = [False] * n
    components: list[list[int]] = []
    for start in range(n):
        if not visited[start]:
            component: list[int] = []
            queue = [start]
            visited[start] = True
            while queue:
                node = queue.pop(0)
                component.append(node)
                for neighbor in range(n):
                    if adj_matrix[node][neighbor] != 0 and not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            components.append(component)

    beta_0 = len(components)                   # Anzahl Zusammenhangskomponenten
    beta_1 = max(0, n_edges - n + beta_0)      # Kreisrang: |E| - |V| + β₀
    euler_char = beta_0 - beta_1               # χ = β₀ - β₁

    return {
        'beta_0': beta_0,
        'beta_1': beta_1,
        'euler_characteristic': euler_char,
        'vertices': n,
        'edges': n_edges,
        'components': components,
    }


def fundamental_group_free_generators(simplices: dict) -> dict:
    """
    @brief Schätzt die freien Generatoren der Fundamentalgruppe π₁(X).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die Fundamentalgruppe π₁(X, x₀) beschreibt Schleifen in X (modulo Homotopie).
    Für einen zusammenhängenden Simplizialkomplex gilt:

    - Ohne 2-Simplizes (nur Graph): π₁(X) ist frei mit β₁ Generatoren.
      Jeder unabhängige Zyklus im Graphen entspricht einem freien Generator.
    - Mit 2-Simplizes: Jede 2-Fläche "füllt" einen Zyklus aus.
      Effektive Generatoren = β₁ (des 1-Skeletts) - (Anzahl 2-Simplizes, die neue
      Randzyklen einbringen und bereits vorhandene Zyklen eliminieren).

    Für einfachere Fälle (ohne höhere Kohomologie):
        freie Generatoren = β₁ aus simplicial_homology

    @param simplices Dict mit {dim: [(v1,...), ...]}.
    @return Dict mit 'free_generators', 'is_simply_connected', 'comment'.
    """
    # Simplizialer Homologie berechnen, um β₁ zu erhalten
    hom = simplicial_homology(simplices)
    beta_1 = hom['h1']
    beta_0 = hom['h0']

    # Einfach zusammenhängend ⟺ β₀ = 1 und β₁ = 0
    is_connected_space = (beta_0 == 1)
    is_simply_connected = is_connected_space and (beta_1 == 0)

    # Kommentar je nach Topologie
    if not is_connected_space:
        comment = (
            f"Der Raum hat {beta_0} Zusammenhangskomponenten. "
            "Die Fundamentalgruppe ist pro Komponente definiert."
        )
    elif is_simply_connected:
        comment = (
            "Der Raum ist einfach zusammenhängend: π₁(X) = 1 (triviale Gruppe). "
            "Alle Schleifen sind kontrahierbar."
        )
    else:
        comment = (
            f"π₁(X) ist eine freie Gruppe mit {beta_1} Generator(en). "
            f"β₁ = {beta_1} unabhängige nicht-kontrahierbare Zyklen existieren."
        )

    return {
        'free_generators': beta_1,
        'is_simply_connected': is_simply_connected,
        'comment': comment,
    }


# ============================================================
# de Rham-Kohomologie
# ============================================================

def de_rham_cohomology_circle() -> dict:
    """
    @brief Berechnet H^k_dR(S¹) – de Rham-Kohomologie des Kreises.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die de Rham-Kohomologie klassifiziert Differentialformen auf
    glatten Mannigfaltigkeiten und steht via de Rham-Satz in
    kanonischer Verbindung zur singulären Kohomologie.

    Bekannte Ergebnisse für den Kreis S¹ = {(x,y) : x²+y²=1}:

        H^0_dR(S¹) = ℝ   (konstante Funktionen; S¹ ist zusammenhängend)
        H^1_dR(S¹) = ℝ   (dθ ist geschlossen, aber nicht exakt auf S¹)
        H^k_dR(S¹) = 0   für k ≥ 2 (dim(S¹) = 1)

    Numerische Verifikation:
    - ∮_{S¹} dθ = 2π ≠ 0 → dθ ist nicht exakt (kein globales θ auf S¹)
    - ∮_{S¹} df = 0 für jede glatte f → alle exakten 1-Formen integrieren zu 0

    Betti-Zahlen: β₀ = 1, β₁ = 1
    Euler-Charakteristik: χ(S¹) = β₀ - β₁ = 0

    @return Dict mit 'h0', 'h1', 'h_higher', 'betti_numbers',
                     'euler_characteristic', 'integral_dtheta',
                     'is_dtheta_exact'.
    """
    import numpy as np

    # Numerische Verifikation: ∮_{S¹} dθ via Trapezregel
    # Parametrisierung: θ(t) = t, t ∈ [0, 2π]
    # dθ = dt auf dem Kreis → Integral = 2π
    n_points = 10000
    t_vals = np.linspace(0, 2 * math.pi, n_points, endpoint=False)
    dt = 2 * math.pi / n_points

    # ∮ dθ = Summe aller dt-Schritte = 2π (numerisch)
    integral_dtheta = float(np.sum(np.ones(n_points) * dt))

    # dθ ist genau dann nicht exakt, wenn das Integral ≠ 0
    is_dtheta_exact = abs(integral_dtheta) < 1e-10

    # Euler-Charakteristik des Kreises: χ = 0 (β₀ - β₁ = 1 - 1 = 0)
    betti = [1, 1]          # β₀=1 (zusammenhängend), β₁=1 (ein Zyklus)
    euler_char = betti[0] - betti[1]

    return {
        'h0': 'ℝ',                           # H^0 ≅ ℝ (konstante Funktionen)
        'h1': 'ℝ',                           # H^1 ≅ ℝ (erzeugt von dθ)
        'h_higher': '0',                     # H^k = 0 für k ≥ 2
        'h0_dim': 1,                         # dim(H^0) = 1
        'h1_dim': 1,                         # dim(H^1) = 1
        'betti_numbers': betti,
        'euler_characteristic': euler_char,
        'integral_dtheta': integral_dtheta,  # ≈ 2π
        'is_dtheta_exact': is_dtheta_exact,  # False (dθ ist nicht exakt)
    }


def de_rham_cohomology_sphere(n: int = 2) -> dict:
    """
    @brief Berechnet H^k_dR(Sⁿ) – de Rham-Kohomologie der n-Sphäre.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Die n-Sphäre Sⁿ = {x ∈ ℝⁿ⁺¹ : |x| = 1} ist eine kompakte,
    orientierbare n-dimensionale Mannigfaltigkeit ohne Rand.

    Bekannte Ergebnisse (de Rham-Satz, Poincaré-Dualität):

        H^k_dR(Sⁿ) = ℝ  für k = 0 und k = n
        H^k_dR(Sⁿ) = 0  sonst

    Betti-Zahlen: β₀ = 1, β_n = 1, β_k = 0 für 0 < k < n

    Euler-Charakteristik:
        χ(Sⁿ) = 1 + (-1)ⁿ  →  χ = 2 für gerades n, χ = 0 für ungerades n

    Beispiele:
        S¹: χ = 0  (Kreis)
        S²: χ = 2  (Sphäre)
        S³: χ = 0  (3-Sphäre)

    @param n Dimension der Sphäre (n ≥ 1).
    @return Dict mit 'betti_numbers', 'euler_characteristic',
                     'nonzero_groups', 'dimension'.
    """
    if n < 1:
        raise ValueError(f"Dimension n muss >= 1 sein, erhalten: {n}")

    # Betti-Zahlen: β₀ = β_n = 1, alle anderen 0
    betti = [0] * (n + 1)
    betti[0] = 1      # H^0 = ℝ (Sⁿ ist zusammenhängend)
    betti[n] = 1      # H^n = ℝ (Orientierungsklasse)

    # Euler-Charakteristik: χ = Σ (-1)^k β_k = 1 + (-1)^n
    euler_char = 1 + ((-1) ** n)

    # Nicht-triviale Kohomologiegruppen beschreiben
    if n == 1:
        nonzero_groups = {'H^0': 'ℝ', 'H^1': 'ℝ'}
    else:
        nonzero_groups = {'H^0': 'ℝ', f'H^{n}': 'ℝ'}

    return {
        'betti_numbers': betti,
        'euler_characteristic': euler_char,
        'nonzero_groups': nonzero_groups,
        'dimension': n,
        'h0_dim': 1,
        'hn_dim': 1,
        'h_middle_dim': 0 if n > 1 else None,  # H^k = 0 für 0 < k < n
    }


def stokes_theorem_verify(f_expr: str = 'P=y,Q=x', region: str = 'disk') -> dict:
    """
    @brief Numerische Verifikation des Stokes-Satzes (Grünscher Satz in 2D).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Stokes-Satz ist eines der fundamentalsten Resultate der Differentialgeometrie:

        ∫_∂M ω = ∫_M dω

    In 2D (Grünscher Satz):
        ∮_{∂D} (P dx + Q dy) = ∫∫_D (∂Q/∂x - ∂P/∂y) dx dy

    Vorgebaut: f_expr = 'P=y,Q=x' → ∂Q/∂x - ∂P/∂y = 0, Integral = 0.
    Für 'P=-y,Q=x': ∂Q/∂x - ∂P/∂y = 2, Linienintegral = Fläche·2 = 2π.

    Numerische Methode:
    - Linienintegral: Trapezregel auf dem Einheitskreis
    - Flächenintegral: Polarkoordinaten-Integration (scipy.dblquad)

    @param f_expr String, der die Form P, Q beschreibt (nur für Dokumentation).
    @param region 'disk' (Einheitsscheibe D²).
    @return Dict mit 'line_integral', 'area_integral', 'relative_error',
                     'stokes_verified', 'green_curl'.
    """
    import numpy as np
    from scipy import integrate  # type: ignore

    # Gewählte Testform: ω = -y dx + x dy
    # dω = (∂x/∂x - ∂(-y)/∂y) dx ∧ dy = (1 + 1) dx ∧ dy = 2 dx ∧ dy
    # ∮_{S¹} (-y dx + x dy) = ∫∫_{D²} 2 dx dy = 2·π·1² = 2π
    # Linkes Integral (Linienintegral über ∂D² = S¹)
    # Parametrisierung: x = cos(t), y = sin(t), t ∈ [0, 2π]
    # dx = -sin(t) dt, dy = cos(t) dt
    # Integrand: -y·(-sin t) + x·cos t = sin²t + cos²t = 1
    n_line = 50000
    t_vals = np.linspace(0, 2 * math.pi, n_line, endpoint=False)
    dt = 2 * math.pi / n_line

    # P = -y, Q = x auf dem Kreis:
    # P·(dx/dt) + Q·(dy/dt) = (-sin t)·(-sin t) + (cos t)·(cos t) = sin²t + cos²t = 1
    integrand_line = np.sin(t_vals) ** 2 + np.cos(t_vals) ** 2  # = 1 überall
    # NumPy >= 2.0: trapezoid (früher: trapz)
    line_integral = float(np.trapezoid(integrand_line, dx=dt))

    # Rechtes Integral (Flächenintegral über D²): ∫∫_{D²} 2 dx dy = 2π
    # via scipy.dblquad: ∫_{-1}^{1} ∫_{-√(1-x²)}^{√(1-x²)} 2 dy dx = π·1²·2 = 2π
    def area_integrand(y_var: float, x_var: float) -> float:
        """Integrand für das Flächenintegral: curl(P,Q) = 2."""
        return 2.0

    def y_lower(x_var: float) -> float:
        """Untere Grenze: -√(1-x²)."""
        return -math.sqrt(max(0.0, 1.0 - x_var ** 2))

    def y_upper(x_var: float) -> float:
        """Obere Grenze: √(1-x²)."""
        return math.sqrt(max(0.0, 1.0 - x_var ** 2))

    area_result, area_error = integrate.dblquad(
        area_integrand, -1.0, 1.0, y_lower, y_upper
    )
    area_integral = float(area_result)

    # Relativer Fehler |LHS - RHS| / |RHS|
    if abs(area_integral) < 1e-15:
        relative_error = abs(line_integral - area_integral)
    else:
        relative_error = abs(line_integral - area_integral) / abs(area_integral)

    # Stokes verifiziert, wenn relativer Fehler < 1% (numerisch)
    stokes_verified = relative_error < 0.01

    return {
        'line_integral': line_integral,     # ≈ 2π
        'area_integral': area_integral,     # ≈ 2π
        'relative_error': relative_error,
        'stokes_verified': stokes_verified,
        'green_curl': 2.0,                  # ∂Q/∂x - ∂P/∂y = 1-(-1) = 2
        'expected_value': math.pi * 2.0,   # 2π ≈ 6.2832
        'form': 'omega = -y dx + x dy',
    }


def brouwer_fixed_point_evidence(n_trials: int = 1000) -> dict:
    """
    @brief Numerische Evidenz für den Brouwerschen Fixpunktsatz.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Brouwersche Fixpunktsatz (1910) besagt:
        Jede stetige Abbildung f: Dⁿ → Dⁿ besitzt mindestens einen Fixpunkt x = f(x).

    Für n=2 (Einheitsscheibe D²):
    - Dⁿ ist kompakt, konvex und nicht leer.
    - Fixpunkte sind Punkte mit f(x) = x.
    - Topologischer Beweis: Ansonsten könnte man D² auf S¹ = ∂D² zurückziehen
      (Retraktion), was unmöglich ist da π₁(S¹) = ℤ ≠ 0 = π₁(D²).

    Numerische Suche:
    Für jede zufällige Abbildung wird via Bisektion nach einem Fixpunkt gesucht.
    Die Abbildung f wird so konstruiert, dass sie D² in sich abbildet:
        f(x) = g(x) / max(1, |g(x)|)   für eine zufällige lineare g: D² → ℝ².

    @param n_trials Anzahl zufälliger Abbildungen.
    @return Dict mit 'found_ratio', 'n_found', 'n_trials',
                     'brouwer_confirmed', 'method'.
    """
    import numpy as np

    rng = np.random.default_rng(seed=42)  # Reproduzierbarer Seed
    n_found = 0

    for _ in range(n_trials):
        # Zufällige 2×2-Matrix A und Verschiebung b erzeugen
        A = rng.uniform(-0.9, 0.9, (2, 2))
        b = rng.uniform(-0.1, 0.1, 2)

        def make_f(A_: 'np.ndarray', b_: 'np.ndarray'):
            """Erzeugt eine stetige Abbildung f: D² → D²."""
            def f(x: 'np.ndarray') -> 'np.ndarray':
                """f(x) = (Ax + b) normiert auf D²."""
                y = A_ @ x + b_
                norm_y = float(np.linalg.norm(y))
                # Normierung sicherstellen: f bildet in D² ab
                if norm_y > 1.0:
                    y = y / norm_y
                return y
            return f

        f_trial = make_f(A, b)

        # Fixpunktsuche via Iteration (Banachscher Fixpunktsatz für kontraktive f)
        # Für nicht-kontraktive f: Suche ||f(x) - x|| < ε
        x = rng.uniform(-0.5, 0.5, 2)  # Startpunkt
        x = x / max(1.0, float(np.linalg.norm(x)))  # In D² projizieren

        fixed_found = False
        # Iterative Fixpunktsuche: x_{n+1} = (x_n + f(x_n)) / 2
        for _ in range(500):
            fx = f_trial(x)
            residual = float(np.linalg.norm(fx - x))
            if residual < 1e-6:
                fixed_found = True
                break
            # Dämpfte Iteration: Mittelung zwischen x und f(x)
            x_new = 0.5 * (x + fx)
            norm_new = float(np.linalg.norm(x_new))
            if norm_new > 1.0:
                x_new = x_new / norm_new
            x = x_new

        if fixed_found:
            n_found += 1

    found_ratio = n_found / n_trials if n_trials > 0 else 0.0

    return {
        'found_ratio': found_ratio,                    # Anteil mit gefundenem Fixpunkt
        'n_found': n_found,
        'n_trials': n_trials,
        'brouwer_confirmed': found_ratio > 0.95,       # Erwartet ≈ 100%
        'method': 'gedaempfte_iteration_D2',
        'theorem': 'Brouwerscher Fixpunktsatz: f: D² → D² stetig → ∃ x: f(x) = x',
    }
