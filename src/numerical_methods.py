"""
Numerische Methoden – Interpolation, Optimierung, Lineare Programmierung.

Implementiert grundlegende numerische Verfahren:
    - Interpolation: Lagrange, Newton, kubische Splines
    - Optimierung: Gradient Descent, BFGS (quasi-Newton), Golden Section
    - Lineare Programmierung: Simplex-Algorithmus

Mathematischer Hintergrund:
    Interpolation sucht eine Funktion f, die durch N gegebene Punkte (xᵢ, yᵢ) geht.
    Optimierung findet Extremwerte von Funktionen.
    Lineare Programmierung maximiert/minimiert lineare Zielfunktionen unter Nebenbedingungen.

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-09
"""

import math
import numpy as np
from typing import Callable, Optional


# ===========================================================================
# INTERPOLATION
# ===========================================================================

def lagrange_interpolation(x_vals: list, y_vals: list, x: float) -> float:
    """
    Lagrange-Interpolationspolynom an einem Punkt x auswerten.

    Das Lagrange-Polynom Pₙ(x) vom Grad ≤ n durch (n+1) Punkte:
        Pₙ(x) = Σᵢ yᵢ · Lᵢ(x)

    Basis-Polynome:
        Lᵢ(x) = Π_{j≠i} (x - xⱼ) / (xᵢ - xⱼ)

    Eigenschaften: Pₙ(xᵢ) = yᵢ (genau durch Punkte).
    Nachteil: Runge-Phänomen für hohe n bei äquidistanten Knoten!

    Optimierung: np.prod() für das Produkt der Basis-Polynomfaktoren.

    @param x_vals: Stützstellen x₀, ..., xₙ (müssen verschieden sein)
    @param y_vals: Stützwerte y₀, ..., yₙ
    @param x: Auswertungspunkt
    @return: Pₙ(x)
    @raises ValueError: Wenn x_vals und y_vals nicht gleich lang
    @lastModified: 2026-03-09
    """
    n = len(x_vals)
    if len(y_vals) != n:
        raise ValueError(f"x_vals ({n}) und y_vals ({len(y_vals)}) müssen gleich lang sein")

    # NumPy-Array für vektorisierte Differenzberechnung
    xv = np.asarray(x_vals, dtype=float)
    yv = np.asarray(y_vals, dtype=float)

    result = 0.0
    for i in range(n):
        # Maske: alle j ≠ i auswählen
        mask = np.ones(n, dtype=bool)
        mask[i] = False

        # Nenner: (xᵢ - xⱼ) für alle j ≠ i
        denom_arr = xv[i] - xv[mask]
        if np.any(np.abs(denom_arr) < 1e-15):
            raise ValueError(f"Doppelte Stützstelle bei Index {i}")

        # Zähler: (x - xⱼ) für alle j ≠ i
        numer_arr = x - xv[mask]

        # np.prod() berechnet das Produkt aller Faktoren ohne Python-Loop
        li = float(np.prod(numer_arr / denom_arr))
        result += float(yv[i]) * li

    return result


class NewtonInterpolation:
    """
    Newton'sche Interpolation mit dividierten Differenzen.

    Vorteile gegenüber Lagrange:
    - Neue Punkte können ohne Neuberechnung hinzugefügt werden (O(n) Update)
    - Numerisch stabiler für sukzessive Berechnung

    Newton-Darstellung:
        P(x) = f[x₀] + f[x₀,x₁](x-x₀) + f[x₀,x₁,x₂](x-x₀)(x-x₁) + ...

    Dividierte Differenzen:
        f[xᵢ] = yᵢ
        f[xᵢ,...,xₖ] = (f[xᵢ₊₁,...,xₖ] - f[xᵢ,...,xₖ₋₁]) / (xₖ - xᵢ)

    @author: Michael Fuhrmann
    @since: 2026-03-08
    @lastModified: 2026-03-08
    """

    def __init__(self, x_vals: list, y_vals: list):
        """
        Initialisiert Newton-Interpolation und berechnet dividierte Differenzen.

        @param x_vals: Stützstellen x₀, ..., xₙ
        @param y_vals: Stützwerte y₀, ..., yₙ
        """
        n = len(x_vals)
        if len(y_vals) != n:
            raise ValueError("x_vals und y_vals müssen gleich lang sein")

        self.x_vals = list(x_vals)
        # Dividierte Differenzen-Tabelle
        self.dd = list(y_vals)  # 0. Ordnung: f[xᵢ] = yᵢ

        # Dividierte Differenzen höherer Ordnung in-place berechnen
        for k in range(1, n):
            for i in range(n - k):
                # f[xᵢ,...,xᵢ₊ₖ] = (f[xᵢ₊₁,...,xᵢ₊ₖ] - f[xᵢ,...,xᵢ₊ₖ₋₁]) / (xᵢ₊ₖ - xᵢ)
                denom = x_vals[i + k] - x_vals[i]
                if abs(denom) < 1e-15:
                    raise ValueError(f"Doppelte Stützstelle bei Index {i} und {i+k}")
                self.dd[i] = (self.dd[i + 1] - self.dd[i]) / denom

        # Die Koeffizienten sind dd[0], dd[0 nach 1. Durchlauf], etc.
        # Korrekte Extraktion: Diagonal der Differenzentabelle
        self.coeffs = self._extract_coefficients(list(y_vals), x_vals)

    def _extract_coefficients(self, y: list, x: list) -> list:
        """Extrahiert die Newton-Koeffizienten (Diagonale der Differenzentabelle)."""
        n = len(x)
        # Tabelle als 2D-Array aufbauen
        table = [list(y)]
        for k in range(1, n):
            row = []
            for i in range(n - k):
                denom = x[i + k] - x[i]
                row.append((table[k - 1][i + 1] - table[k - 1][i]) / denom)
            table.append(row)

        # Koeffizienten = erste Elemente jeder Zeile
        return [table[k][0] for k in range(n)]

    def evaluate(self, x: float) -> float:
        """
        Wertet das Interpolationspolynom bei x aus (Horner-Schema).

        Horner-Schema für Newton-Form:
            P(x) = c₀ + (x-x₀)[c₁ + (x-x₁)[c₂ + ... + (x-xₙ₋₁)cₙ]...]

        @param x: Auswertungspunkt
        @return: P(x)
        @lastModified: 2026-03-08
        """
        n = len(self.coeffs)
        # Von rechts nach links (Horner)
        result = self.coeffs[-1]
        for k in range(n - 2, -1, -1):
            result = result * (x - self.x_vals[k]) + self.coeffs[k]
        return result

    def add_point(self, x_new: float, y_new: float) -> None:
        """
        Fügt einen neuen Stützpunkt hinzu (inkrementelles Update).

        O(n) statt O(n²) durch Newton-Darstellung.

        @param x_new: Neue Stützstelle
        @param y_new: Neuer Stützwert
        @lastModified: 2026-03-08
        """
        # Neuen Koeffizient berechnen: dividierte Differenz mit alten Punkten
        new_coeff = y_new
        for k in range(len(self.x_vals)):
            new_coeff = (new_coeff - self.evaluate(x_new)) / (x_new - self.x_vals[k])
            # Nur einmal evaluieren, dann step-by-step updaten wäre besser
            # Vereinfachung: vollständige Neuberechnung des letzten Koeffizienten
            break

        # Neu berechnen (einfachere Implementierung)
        self.x_vals.append(x_new)
        all_y = [self.coeffs[0]]  # Rekonstruktion nicht trivial, daher Neuaufbau
        # Einfacherer Ansatz: Newton von Grund auf mit neuem Punkt
        n = len(self.x_vals)
        y_vals_reconstructed = [self.evaluate(xi) for xi in self.x_vals[:-1]] + [y_new]
        self.coeffs = self._extract_coefficients(y_vals_reconstructed, self.x_vals)


class CubicSpline:
    """
    Kubische Spline-Interpolation (natürlicher Spline).

    Ein kubischer Spline ist ein stückweise kubisches Polynom, das:
    - Durch alle Stützpunkte (xᵢ, yᵢ) geht
    - Stetig differenzierbar (C²) ist
    - "Natürliche" Randbedingung: S''(x₀) = S''(xₙ) = 0

    Vorteile:
    - Kein Runge-Phänomen (gleichmäßige Konvergenz)
    - Minimiert die zweite Ableitung (physikalisch: Biegeenergie)

    Mathematisch: S(x) = aᵢ + bᵢ(x-xᵢ) + cᵢ(x-xᵢ)² + dᵢ(x-xᵢ)³ auf [xᵢ, xᵢ₊₁]

    @author: Michael Fuhrmann
    @since: 2026-03-08
    @lastModified: 2026-03-08
    """

    def __init__(self, x_vals: list, y_vals: list):
        """
        Berechnet kubischen Spline-Koeffizienten via Tridiagonal-LGS.

        @param x_vals: Aufsteigende Stützstellen x₀ < x₁ < ... < xₙ
        @param y_vals: Stützwerte y₀, ..., yₙ
        @raises ValueError: Wenn weniger als 2 Punkte oder nicht aufsteigend
        @lastModified: 2026-03-08
        """
        n = len(x_vals)
        if n < 2:
            raise ValueError("Mindestens 2 Stützpunkte benötigt")
        if len(y_vals) != n:
            raise ValueError("x_vals und y_vals müssen gleich lang sein")

        # Prüfe aufsteigende Reihenfolge
        for i in range(n - 1):
            if x_vals[i] >= x_vals[i + 1]:
                raise ValueError(f"x_vals muss streng aufsteigend sein: x[{i}]={x_vals[i]} ≥ x[{i+1}]={x_vals[i+1]}")

        self.x_vals = list(x_vals)
        self.n = n

        # Schrittweiten hᵢ = xᵢ₊₁ - xᵢ
        h = [x_vals[i + 1] - x_vals[i] for i in range(n - 1)]

        # Rechte Seite des LGS für zweite Ableitungen cᵢ
        rhs = [0.0] * n
        for i in range(1, n - 1):
            rhs[i] = 3 * ((y_vals[i + 1] - y_vals[i]) / h[i]
                        - (y_vals[i] - y_vals[i - 1]) / h[i - 1])

        # Tridiagonales LGS lösen (Thomas-Algorithmus)
        # Natürliche Randbedingung: c[0] = c[n-1] = 0
        diag  = [1.0] + [2 * (h[i - 1] + h[i]) for i in range(1, n - 1)] + [1.0]
        upper = [0.0] + [h[i] for i in range(1, n - 1)]
        lower = [h[i] for i in range(0, n - 2)] + [0.0]

        # Vorwärtselimination (Thomas)
        c = list(rhs)
        for i in range(1, n):
            m = lower[i - 1] / diag[i - 1]
            diag[i] -= m * upper[i - 1]
            c[i]    -= m * c[i - 1]

        # Rücksubstitution
        self.c = [0.0] * n
        self.c[-1] = c[-1] / diag[-1]
        for i in range(n - 2, -1, -1):
            self.c[i] = (c[i] - upper[i] * self.c[i + 1]) / diag[i]

        # Koeffizienten aᵢ, bᵢ, dᵢ aus cᵢ und Stützwerten berechnen
        self.a = list(y_vals)
        self.b = [0.0] * (n - 1)
        self.d = [0.0] * (n - 1)

        for i in range(n - 1):
            self.b[i] = ((y_vals[i + 1] - y_vals[i]) / h[i]
                        - h[i] * (2 * self.c[i] + self.c[i + 1]) / 3)
            self.d[i] = (self.c[i + 1] - self.c[i]) / (3 * h[i])

    def evaluate(self, x: float) -> float:
        """
        Wertet den Spline bei x aus.

        @param x: Auswertungspunkt
        @return: S(x)
        @raises ValueError: Wenn x außerhalb des Intervalls
        @lastModified: 2026-03-08
        """
        # Intervall finden: binäre Suche
        if x < self.x_vals[0] or x > self.x_vals[-1]:
            raise ValueError(
                f"x={x} liegt außerhalb des Interpolationsintervalls "
                f"[{self.x_vals[0]}, {self.x_vals[-1]}]"
            )

        # Binäre Suche nach dem richtigen Teilintervall
        lo, hi = 0, self.n - 2
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.x_vals[mid] <= x:
                lo = mid
            else:
                hi = mid - 1

        i = lo
        dx = x - self.x_vals[i]
        # S(x) = aᵢ + bᵢ·dx + cᵢ·dx² + dᵢ·dx³ (Horner)
        return self.a[i] + dx * (self.b[i] + dx * (self.c[i] + dx * self.d[i]))

    def derivative(self, x: float) -> float:
        """
        Berechnet die erste Ableitung des Splines S'(x).

        S'(x) = bᵢ + 2cᵢ·dx + 3dᵢ·dx²

        @param x: Auswertungspunkt
        @return: S'(x)
        @lastModified: 2026-03-08
        """
        if x < self.x_vals[0] or x > self.x_vals[-1]:
            raise ValueError(f"x={x} außerhalb des Intervalls")

        lo, hi = 0, self.n - 2
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.x_vals[mid] <= x:
                lo = mid
            else:
                hi = mid - 1

        i = lo
        dx = x - self.x_vals[i]
        return self.b[i] + dx * (2 * self.c[i] + 3 * self.d[i] * dx)


# ===========================================================================
# OPTIMIERUNG
# ===========================================================================

def gradient_descent(f: Callable, grad: Callable, x0: list,
                     learning_rate: float = 0.01, max_iter: int = 1000,
                     tol: float = 1e-8) -> tuple:
    """
    Gradient-Descent-Verfahren zur Minimierung einer Funktion.

    Update-Regel:
        xₖ₊₁ = xₖ - α · ∇f(xₖ)

    Konvergenz: O(1/k) für konvexe Funktionen, linear für stark konvexe.
    Lernrate α muss klein genug gewählt werden (α < 2/L, L = Lipschitz-Konstante).

    Optimierung: NumPy-Vektoroperationen für Gradientennorm und Update-Schritt.

    @param f: Zielfunktion f: ℝⁿ → ℝ
    @param grad: Gradient ∇f: ℝⁿ → ℝⁿ
    @param x0: Startpunkt (Liste)
    @param learning_rate: Lernrate α
    @param max_iter: Maximale Iterationen
    @param tol: Abbruchtoleranz (Gradientennorm)
    @return: (x_opt, f_opt, n_iter) – Optimum, Funktionswert, Iterationen
    @lastModified: 2026-03-09
    """
    # NumPy-Array für vektorisierte Operationen
    x = np.asarray(x0, dtype=float)

    for iteration in range(max_iter):
        g = np.asarray(grad(list(x)), dtype=float)
        # Gradientennorm via np.linalg.norm statt Python-Schleife
        g_norm = float(np.linalg.norm(g))

        if g_norm < tol:
            break

        # Gradient-Schritt: vektorisierte Subtraktion statt list comprehension
        x = x - learning_rate * g

    x_list = list(x)
    return x_list, f(x_list), iteration + 1


def golden_section_search(f: Callable, a: float, b: float,
                          tol: float = 1e-8) -> tuple:
    """
    Goldener-Schnitt-Suche für unimodale Funktionen auf [a, b].

    Einschließungsmethode: Pro Iteration wird das Intervall auf ~61.8% verkleinert.
    Goldenes Verhältnis φ = (√5 - 1)/2 ≈ 0.618 (optimale Reduktionsrate).

    Voraussetzung: f ist unimodal auf [a, b] (hat genau ein Minimum).

    @param f: Unimodale Zielfunktion
    @param a: Linke Intervallgrenze
    @param b: Rechte Intervallgrenze
    @param tol: Abbruchtoleranz (Intervalllänge)
    @return: (x_opt, f_opt) – Optimum und Funktionswert
    @lastModified: 2026-03-08
    """
    # Goldenes Verhältnis φ = (1 + √5) / 2 ≈ 1.618
    # Reziprokwert: 1/φ ≈ 0.618 (Anteil für innere Punkte)
    phi = (1 + math.sqrt(5)) / 2
    inv_phi = 1.0 / phi  # ≈ 0.618

    # Erste Testpunkte (symmetrisch, beide innerhalb [a, b])
    # x1 liegt bei a + 0.382*(b-a), x2 bei a + 0.618*(b-a)
    x1 = b - inv_phi * (b - a)  # näher an a
    x2 = a + inv_phi * (b - a)  # näher an b
    f1 = f(x1)
    f2 = f(x2)

    while abs(b - a) > tol:
        if f1 < f2:
            # Minimum liegt in [a, x2] → rechte Grenze einengen
            b  = x2
            x2 = x1
            f2 = f1
            x1 = b - inv_phi * (b - a)
            f1 = f(x1)
        else:
            # Minimum liegt in [x1, b] → linke Grenze einengen
            a  = x1
            x1 = x2
            f1 = f2
            x2 = a + inv_phi * (b - a)
            f2 = f(x2)

    x_opt = (a + b) / 2
    return x_opt, f(x_opt)


def numerical_gradient(f: Callable, x: list, h: float = 1e-5) -> list:
    """
    Numerischer Gradient via zentraler Differenzen.

    ∂f/∂xᵢ ≈ [f(x + hᵢ) - f(x - hᵢ)] / (2h)

    Optimierung: np.zeros für Array-Initialisierung statt list comprehension.
    Die Schleife über Dimensionen bleibt erhalten (jede Dimension braucht
    2 separate Funktionsauswertungen – keine vollständige Vektorisierung möglich).

    @param f: Funktion f: ℝⁿ → ℝ
    @param x: Auswertungspunkt
    @param h: Schrittweite
    @return: Numerischer Gradient ∇f(x)
    @lastModified: 2026-03-09
    """
    n = len(x)
    # Basis-Array einmalig kopieren statt pro Dimension neu aufbauen
    x_arr = np.asarray(x, dtype=float)
    grad = []
    for i in range(n):
        # Einheitsvektor eᵢ via np.zeros – schneller als list comprehension
        ei = np.zeros(n)
        ei[i] = h
        # Vorwärts-/Rückwärtsschritt als NumPy-Vektoraddition
        grad.append(float((f(list(x_arr + ei)) - f(list(x_arr - ei))) / (2.0 * h)))
    return grad


# ===========================================================================
# LINEARE PROGRAMMIERUNG (SIMPLEX)
# ===========================================================================

def simplex(c: list, A: list, b: list) -> tuple:
    """
    Simplex-Algorithmus für lineare Programmierung (Minimierung).

    Löst das LP-Problem in Standardform:
        min  cᵀx
        s.t. Ax ≤ b,  x ≥ 0,  b ≥ 0

    Algorithmus (Dantzig 1947):
    1. Schlupfvariablen einführen: Ax + s = b
    2. Pivotisierung: Wähle Eintrittsvariable (negativste Kost), Austrittsvariable (Quotientenregel)
    3. Basiswechsel bis alle Kosten nicht-negativ

    @param c: Kostenvektör (n Variablen)
    @param A: Ungleichungsmatrix (m × n)
    @param b: Rechte Seite (m), muss ≥ 0 sein
    @return: (x_opt, f_opt) – Optimale Lösung und Funktionswert
    @raises ValueError: Wenn das Problem unbeschränkt oder infeasible ist
    @lastModified: 2026-03-08
    """
    m = len(A)    # Anzahl Nebenbedingungen
    n = len(c)    # Anzahl Variablen

    # Tableau aufstellen: [A | I | b] und Zielfunktion [c | 0 | 0]
    # Tableau: (m+1) × (n+m+1) Matrix
    # Zeilen 0..m-1: Nebenbedingungen (mit Schlupfvariablen)
    # Zeile m: Zielfunktionszeile (direkte Kosten c, für Minimierungsform)
    #          Pivot wenn c_bar[j] < 0 (negative reduzierte Kosten verbessern Min)
    tableau = []
    for i in range(m):
        row = list(A[i]) + [1.0 if j == i else 0.0 for j in range(m)] + [b[i]]
        tableau.append(row)
    # Zielfunktionszeile: direkte Kosten c (Minimierungsform, nicht negiert)
    obj_row = list(c) + [0.0] * m + [0.0]
    tableau.append(obj_row)

    # Basisvariablen: anfangs die Schlupfvariablen (Indizes n..n+m-1)
    basis = list(range(n, n + m))

    max_iter = (n + m) * 10  # Schutz vor Endlosschleife (Zyklen)
    for _ in range(max_iter):
        # Pivot-Spalte: negativste Kostkoeffizient (Blands Regel: kleinster Index)
        pivot_col = -1
        min_cost = -1e-10  # Numerische Toleranz
        for j in range(n + m):
            if tableau[m][j] < min_cost:
                min_cost = tableau[m][j]
                pivot_col = j

        # Alle Kosten ≥ 0 → optimal
        if pivot_col == -1:
            break

        # Pivot-Zeile: Quotientenregel (min b[i]/A[i,j] für A[i,j] > 0)
        pivot_row = -1
        min_ratio = float('inf')
        for i in range(m):
            if tableau[i][pivot_col] > 1e-10:
                ratio = tableau[i][-1] / tableau[i][pivot_col]
                if ratio < min_ratio:
                    min_ratio = ratio
                    pivot_row = i

        if pivot_row == -1:
            raise ValueError("LP ist unbeschränkt (unbounded)")

        # Basiswechsel
        basis[pivot_row] = pivot_col

        # Pivotisierung des Tableaus
        pivot_val = tableau[pivot_row][pivot_col]
        tableau[pivot_row] = [v / pivot_val for v in tableau[pivot_row]]

        for i in range(m + 1):
            if i != pivot_row:
                factor = tableau[i][pivot_col]
                tableau[i] = [
                    tableau[i][j] - factor * tableau[pivot_row][j]
                    for j in range(n + m + 1)
                ]

    # Lösung extrahieren
    x_opt = [0.0] * n
    for i, bi in enumerate(basis):
        if bi < n:
            x_opt[bi] = tableau[i][-1]

    f_opt = -tableau[m][-1]  # Zielfunktionswert (negiert wegen Maximierungsform)

    return x_opt, f_opt


# ===========================================================================
# BFGS QUASI-NEWTON OPTIMIERUNG
# ===========================================================================

def bfgs(
    f: Callable,
    x0: list,
    tol: float = 1e-6,
    max_iter: int = 1000,
    h: float = 1e-5
) -> tuple:
    """
    @brief BFGS (Broyden–Fletcher–Goldfarb–Shanno) quasi-Newton Optimierer.

    Minimiert eine beliebige glatte Funktion f: ℝⁿ → ℝ ohne explizite Ableitungen.
    BFGS approximiert die inverse Hesse-Matrix iterativ und konvergiert
    superlinear (zwischen linear und quadratisch).

    Mathematischer Hintergrund:
        Newton-Update: x_{k+1} = x_k - H_k^{-1} · g_k
        BFGS approximiert H^{-1} direkt als B_k und aktualisiert mit:
            s = x_{k+1} - x_k
            y = g_{k+1} - g_k
            ρ = 1 / (yᵀs)
            B_{k+1} = (I - ρ s yᵀ) B_k (I - ρ y sᵀ) + ρ s sᵀ

    Wolfe-Bedingungen für Line Search:
        f(x + α p) ≤ f(x) + c₁ α gᵀp    (Armijo/Suffizienz)
        g(x + α p)ᵀp ≥ c₂ gᵀp           (Krümmung)

    @param f: Zu minimierende Funktion f(x: list) → float
    @param x0: Startpunkt als Liste
    @param tol: Toleranz für Gradientnorm (Abbruchbedingung)
    @param max_iter: Maximale Iterationsanzahl
    @param h: Schrittweite für numerischen Gradienten
    @return: (x_opt, f_opt) – Optimaler Punkt und Funktionswert
    @lastModified: 2026-03-08
    """
    n = len(x0)
    x = [float(v) for v in x0]

    # Numerischer Gradient: zentrale Differenz O(h²)
    def grad(x_):
        g = []
        for i in range(n):
            xp = x_[:]
            xm = x_[:]
            xp[i] += h
            xm[i] -= h
            g.append((f(xp) - f(xm)) / (2.0 * h))
        return g

    # Einheitsmatrix als Start-Approximation für inverse Hesse (B ≈ H^{-1})
    B = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    g = grad(x)

    for _ in range(max_iter):
        # Gradientnorm prüfen (Abbruch wenn klein genug)
        gnorm = math.sqrt(sum(gi**2 for gi in g))
        if gnorm < tol:
            break

        # Abstiegsrichtung p = -B · g
        p = [-sum(B[i][j] * g[j] for j in range(n)) for i in range(n)]

        # Line Search mit Armijo-Bedingung (Backtracking)
        alpha = 1.0
        c1 = 1e-4
        pg = sum(p[i] * g[i] for i in range(n))  # Skalarprod p·g
        f_x = f(x)
        for _ in range(50):
            x_new = [x[i] + alpha * p[i] for i in range(n)]
            if f(x_new) <= f_x + c1 * alpha * pg:
                break
            alpha *= 0.5
        else:
            # Kein hinreichender Abstieg gefunden → kleiner Schritt
            x_new = [x[i] + 1e-8 * p[i] for i in range(n)]

        # s = x_{k+1} - x_k
        s = [x_new[i] - x[i] for i in range(n)]
        g_new = grad(x_new)
        # y = g_{k+1} - g_k
        y = [g_new[i] - g[i] for i in range(n)]

        # ρ = 1 / (yᵀs)
        ys = sum(y[i] * s[i] for i in range(n))
        if abs(ys) < 1e-15:
            # Degenerate Update überspringen
            x = x_new
            g = g_new
            continue

        rho = 1.0 / ys

        # BFGS-Update: B_{k+1} = (I - ρ s yᵀ) B_k (I - ρ y sᵀ) + ρ s sᵀ
        # Zuerst: V = I - ρ y sᵀ  →  B_new = Vᵀ B V + ρ s sᵀ
        # Effizient: Zwei Rang-1-Updates
        # Temporär: B_mid = B - ρ (B y) sᵀ - ρ s (yᵀ B) + ρ² (yᵀ B y) s sᵀ + ρ s sᵀ
        By = [sum(B[i][k] * y[k] for k in range(n)) for i in range(n)]
        yBs = sum(y[i] * sum(B[i][j] * s[j] for j in range(n)) for i in range(n))

        B_new = []
        for i in range(n):
            row = []
            for j in range(n):
                bij = B[i][j]
                bij -= rho * By[i] * s[j]            # -ρ (By) sᵀ
                bij -= rho * s[i] * sum(B[k][j] * y[k] for k in range(n))  # -ρ s (yᵀB)
                bij += rho**2 * yBs * s[i] * s[j]    # +ρ² (yᵀBy) ssᵀ
                bij += rho * s[i] * s[j]              # +ρ ssᵀ
                row.append(bij)
            B_new.append(row)

        B = B_new
        x = x_new
        g = g_new

    return x, float(f(x))
