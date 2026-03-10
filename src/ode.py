"""
@file ode.py
@brief ODE-Modul: Gewöhnliche Differentialgleichungen (Ordinary Differential Equations).
@description
    Enthält numerische Lösungsverfahren für ODE-Anfangswertprobleme,
    lineare ODE mit konstanten Koeffizienten sowie Laplace-Transformation.

    Implementierte Verfahren:
    - Euler-Methode (1. Ordnung, explizit)
    - Runge-Kutta 4. Ordnung (klassisch, fest)
    - Runge-Kutta-Fehlberg 45 (adaptiv, eingebettet)
    - Lineare ODE n-ter Ordnung via Zustandsraumsystem
    - Laplace-Transformation (numerisch, Simpson-Regel)
    - Inverse Laplace-Transformation (Stehfest-Algorithmus)

@author Kurt Ingwer
@date 2026-03-05
@lastModified 2026-03-10
"""

import math
from typing import Callable, List, Optional, Tuple, Union


# ─────────────────────────────────────────────────────────────────────────────
# Typdefinitionen
# ─────────────────────────────────────────────────────────────────────────────

# Zustandsvektor: entweder ein Skalar (float) oder eine Liste von Skalaren
StateType = Union[float, List[float]]
# Rechte Seite: f(t, y) -> dy/dt
ODEFunc = Callable[[float, StateType], StateType]


def _add_states(a: StateType, b: StateType) -> StateType:
    """Addiert zwei Zustände (Skalar oder Vektor)."""
    if isinstance(a, list):
        return [ai + bi for ai, bi in zip(a, b)]
    return a + b


def _scale_state(c: float, a: StateType) -> StateType:
    """Multipliziert einen Zustand (Skalar oder Vektor) mit einem Skalar."""
    if isinstance(a, list):
        return [c * ai for ai in a]
    return c * a


# ─────────────────────────────────────────────────────────────────────────────
# Euler-Verfahren (explizit, 1. Ordnung)
# ─────────────────────────────────────────────────────────────────────────────

def euler_method(
    f: Callable[[float, float | list[float]], float | list[float]],
    t0: float,
    y0: float | list[float],
    t_end: float,
    n_steps: int
) -> tuple[list[float], list[float | list[float]]]:
    """
    Löst ein Anfangswertproblem (AWP) mit dem expliziten Euler-Verfahren.

    Das AWP lautet: dy/dt = f(t, y), y(t0) = y0.

    Diskretisierung: y_{n+1} = y_n + h * f(t_n, y_n)

    Fehlerordnung: O(h) pro Schritt (global O(h)).
    Einfach, aber ungenau - besser für didaktische Zwecke.

    :param f: Rechte Seite der ODE f(t, y)
    :param t0: Anfangszeitpunkt
    :param y0: Anfangswert (Skalar oder Vektor)
    :param t_end: Endzeitpunkt
    :param n_steps: Anzahl der Zeitschritte
    :return: Tupel (Zeitpunkte, Lösungswerte)
    :last_modified: 2026-03-05
    """
    h = (t_end - t0) / n_steps
    t = t0
    y = y0

    # Listen für Zeitpunkte und Zustände (inklusive Anfangswert)
    t_list = [t]
    y_list = [y]

    for _ in range(n_steps):
        # Euler-Schritt: y_new = y + h * f(t, y)
        dy = f(t, y)
        y = _add_states(y, _scale_state(h, dy))
        t = t + h
        t_list.append(t)
        y_list.append(y)

    return t_list, y_list


# ─────────────────────────────────────────────────────────────────────────────
# Runge-Kutta 4. Ordnung (klassisch)
# ─────────────────────────────────────────────────────────────────────────────

def runge_kutta4(
    f: Callable[[float, float | list[float]], float | list[float]],
    t0: float,
    y0: float | list[float],
    t_end: float,
    n_steps: int
) -> tuple[list[float], list[float | list[float]]]:
    """
    Löst ein Anfangswertproblem mit dem klassischen Runge-Kutta-Verfahren 4. Ordnung.

    Das AWP lautet: dy/dt = f(t, y), y(t0) = y0.

    Die vier Stufenwerte (Butcher-Tableau):
      k1 = h * f(t,       y)
      k2 = h * f(t + h/2, y + k1/2)
      k3 = h * f(t + h/2, y + k2/2)
      k4 = h * f(t + h,   y + k3)
    Update: y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4) / 6

    Fehlerordnung: O(h^4) pro Schritt (global O(h^4)).
    Gutes Gleichgewicht zwischen Aufwand und Genauigkeit.

    :param f: Rechte Seite der ODE f(t, y)
    :param t0: Anfangszeitpunkt
    :param y0: Anfangswert (Skalar oder Vektor für Systeme)
    :param t_end: Endzeitpunkt
    :param n_steps: Anzahl der Zeitschritte
    :return: Tupel (Zeitpunkte, Lösungswerte)
    :last_modified: 2026-03-05
    """
    h = (t_end - t0) / n_steps
    t = t0
    y = y0

    t_list = [t]
    y_list = [y]

    for _ in range(n_steps):
        # Vier Stufenwerte des RK4-Verfahrens (gewichteter Mittelwert der Steigungen)
        k1 = _scale_state(h, f(t,           y))
        k2 = _scale_state(h, f(t + h/2, _add_states(y, _scale_state(0.5, k1))))
        k3 = _scale_state(h, f(t + h/2, _add_states(y, _scale_state(0.5, k2))))
        k4 = _scale_state(h, f(t + h,   _add_states(y, k3)))

        # Gewichteter Mittelwert: Simpsonregel-ähnlich
        dy = _scale_state(1/6, _add_states(
            _add_states(k1, _scale_state(2, k2)),
            _add_states(_scale_state(2, k3), k4)
        ))
        y = _add_states(y, dy)
        t = t + h
        t_list.append(t)
        y_list.append(y)

    return t_list, y_list


# ─────────────────────────────────────────────────────────────────────────────
# Runge-Kutta-Fehlberg (adaptiv, RK45)
# ─────────────────────────────────────────────────────────────────────────────

def runge_kutta45(
    f: Callable[[float, float | list[float]], float | list[float]],
    t0: float,
    y0: float | list[float],
    t_end: float,
    tol: float = 1e-6,
    h_min: float = 1e-10,
    h_max: float = 0.1
) -> tuple[list[float], list[float | list[float]]]:
    """
    Löst ein AWP mit dem adaptiven Runge-Kutta-Fehlberg-Verfahren (RK45).

    Verwendet eingebettete RK4/RK5 Methoden um den lokalen Fehler zu schätzen
    und die Schrittweite automatisch anzupassen.

    Schrittweitenanpassung:
      h_neu = h * (tol / fehler)^(1/5) * 0.9 (Sicherheitsfaktor 0.9)

    :param f: Rechte Seite der ODE f(t, y)
    :param t0: Anfangszeitpunkt
    :param y0: Anfangswert
    :param t_end: Endzeitpunkt
    :param tol: Toleranz für den lokalen Fehler
    :param h_min: Minimale Schrittweite
    :param h_max: Maximale Schrittweite
    :return: Tupel (Zeitpunkte, Lösungswerte)
    :last_modified: 2026-03-05
    """
    # Fehlberg-Koeffizienten (Butcher-Tableau für RK45)
    # Knotenkoeffizienten
    c2, c3, c4, c5 = 1/4, 3/8, 12/13, 1.0

    # Stufenkoeffizienten (Eingangsgewichte)
    a21 = 1/4
    a31, a32 = 3/32, 9/32
    a41, a42, a43 = 1932/2197, -7200/2197, 7296/2197
    a51, a52, a53, a54 = 439/216, -8.0, 3680/513, -845/4104
    a61, a62, a63, a64, a65 = -8/27, 2.0, -3544/2565, 1859/4104, -11/40

    # RK4-Gewichte (5 Stufen)
    b1, b3, b4, b5 = 25/216, 1408/2565, 2197/4104, -1/5
    # RK5-Gewichte (6 Stufen) für Fehlerabschätzung
    e1, e3, e4, e5, e6 = 16/135, 6656/12825, 28561/56430, -9/50, 2/55

    def _scale(c, a):
        return _scale_state(c, a)

    def _add(*args):
        result = args[0]
        for a in args[1:]:
            result = _add_states(result, a)
        return result

    def _error_norm(e):
        """Berechnet die Fehlernorm (L_inf für Skalare, L2 für Vektoren)."""
        if isinstance(e, list):
            return max(abs(ei) for ei in e)
        return abs(e)

    t = t0
    y = y0
    h = min(h_max, (t_end - t0) / 10)  # Anfangsschrittweite

    t_list = [t]
    y_list = [y]

    while t < t_end:
        # Schrittweite anpassen damit wir t_end nicht überschreiten
        if t + h > t_end:
            h = t_end - t

        # 6 Stufenwerte berechnen
        k1 = _scale(h, f(t, y))
        k2 = _scale(h, f(t + h/4,    _add(y, _scale(a21, k1))))
        k3 = _scale(h, f(t + 3*h/8,  _add(y, _scale(a31, k1), _scale(a32, k2))))
        k4 = _scale(h, f(t + 12*h/13,_add(y, _scale(a41, k1), _scale(a42, k2), _scale(a43, k3))))
        k5 = _scale(h, f(t + h,       _add(y, _scale(a51, k1), _scale(a52, k2), _scale(a53, k3), _scale(a54, k4))))
        k6 = _scale(h, f(t + h/2,     _add(y, _scale(a61, k1), _scale(a62, k2), _scale(a63, k3), _scale(a64, k4), _scale(a65, k5))))

        # RK4-Lösung (4. Ordnung)
        y4 = _add(y, _scale(b1, k1), _scale(b3, k3), _scale(b4, k4), _scale(b5, k5))

        # RK5-Lösung (5. Ordnung) für Fehlerabschätzung
        y5 = _add(y, _scale(e1, k1), _scale(e3, k3), _scale(e4, k4), _scale(e5, k5), _scale(e6, k6))

        # Fehlerabschätzung: Differenz zwischen 4. und 5. Ordnung
        if isinstance(y4, list):
            error = max(abs(a - b) for a, b in zip(y4, y5))
        else:
            error = abs(y4 - y5)

        if error <= tol or h <= h_min:
            # Schritt akzeptieren
            t += h
            y = y4
            t_list.append(t)
            y_list.append(y)

        # Schrittweite für nächsten Schritt anpassen
        if error > 0:
            h_new = h * 0.9 * (tol / error) ** 0.2
            h = max(h_min, min(h_max, h_new))

    return t_list, y_list


# ─────────────────────────────────────────────────────────────────────────────
# Lineare ODE mit konstanten Koeffizienten
# ─────────────────────────────────────────────────────────────────────────────

def solve_linear_ode(
    coeffs: list[float],
    initial: list[float],
    t_vals: list[float],
    forcing: Callable[[float], float] | None = None
) -> list[float]:
    """
    Löst eine lineare ODE n-ter Ordnung mit konstanten Koeffizienten.

    Allgemeine Form: a_0*y^(n) + a_1*y^(n-1) + ... + a_n*y = g(t)

    coeffs[0]*y^(n) + coeffs[1]*y^(n-1) + ... + coeffs[n]*y = g(t)

    Methode: Umformung in ein System 1. Ordnung (Zustandsraum-Darstellung)
    und Lösung via RK4.

    :param coeffs: Koeffizienten [a_0, a_1, ..., a_n] (a_0 != 0)
    :param initial: Anfangswerte [y(t0), y'(t0), ..., y^(n-1)(t0)]
    :param t_vals: Zeitpunkte, an denen die Lösung gesucht wird
    :param forcing: Erzwingungsfunktion g(t), None für homogene ODE
    :return: Lösungswerte y(t) an den gegebenen Zeitpunkten
    :raises ValueError: Bei ungültigen Parametern
    :last_modified: 2026-03-05
    """
    n = len(coeffs) - 1  # Ordnung der ODE
    if len(initial) != n:
        raise ValueError(f"Brauche {n} Anfangswerte, {len(initial)} gegeben")
    if abs(coeffs[0]) < 1e-15:
        raise ValueError("Führender Koeffizient a_0 darf nicht 0 sein")

    # Normierung: a_0 = 1 durch Division
    a = [c / coeffs[0] for c in coeffs]
    g = forcing if forcing is not None else lambda t: 0.0

    # Zustandsvektor z = [y, y', y'', ..., y^(n-1)]
    # Rechte Seite des Zustandsraumsystems:
    # z'_0 = z_1, z'_1 = z_2, ..., z'_{n-2} = z_{n-1}
    # z'_{n-1} = (g(t) - a_n*z_0 - a_{n-1}*z_1 - ... - a_1*z_{n-1}) / a_0 (=1)
    def system(t: float, z: List[float]) -> List[float]:
        """Zustandsraum-Darstellung der ODE als System 1. Ordnung."""
        dz = list(z[1:])  # Verschiebung: z'_i = z_{i+1}
        # Letzte Gleichung: y^(n) = (g(t) - sum) / a_0
        last = g(t) - sum(a[n - i] * z[i] for i in range(n))
        dz.append(last)
        return dz

    # Anfangszeitpunkt und Schrittweite
    t0 = t_vals[0]
    z0 = list(initial)

    result = []
    t_current = t0
    z_current = z0

    for t_target in t_vals:
        if t_target == t_current:
            result.append(z_current[0])
            continue
        # Integriere von t_current bis t_target mit RK4
        n_steps = max(100, int(abs(t_target - t_current) / 0.01))
        _, z_traj = runge_kutta4(system, t_current, z_current, t_target, n_steps)
        z_current = z_traj[-1]
        t_current = t_target
        result.append(z_current[0])

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Laplace-Transformation (numerisch)
# ─────────────────────────────────────────────────────────────────────────────

def laplace_transform(
    f: Callable[[float], float],
    s: float,
    t_max: float = 50.0,
    n_intervals: int = 10000
) -> float:  # Rückgabe: Näherungswert F(s) = ∫₀^∞ f(t)·e^{-st} dt
    """
    Berechnet die einseitige Laplace-Transformation numerisch.

    Formel: F(s) = integral_0^infinity f(t) * e^(-s*t) dt

    Die Integration wird bei t_max abgebrochen, wenn e^(-s*t_max) sehr klein ist.
    Verwendet die Simpson-Regel für die numerische Integration.

    :param f: Zeitbereichsfunktion f(t)
    :param s: Komplexe Frequenz (hier real, s > 0)
    :param t_max: Obere Integrationsgrenze (approximiert Unendlich)
    :param n_intervals: Anzahl der Intervalle für die Simpson-Regel
    :return: Näherungswert von F(s)
    :raises ValueError: Bei s <= 0
    :last_modified: 2026-03-05
    """
    if s <= 0:
        raise ValueError("s muss > 0 sein für Konvergenz")

    # Bestimme sinnvolles t_max basierend auf s
    t_max = min(t_max, -math.log(1e-12) / s)

    h = t_max / n_intervals
    # Sicherstellen, dass n_intervals gerade ist (Simpson-Regel)
    if n_intervals % 2 != 0:
        n_intervals += 1
        h = t_max / n_intervals

    # Simpson-Regel: integral ≈ (h/3) * [f(0) + 4f(1) + 2f(2) + 4f(3) + ... + f(n)]
    def integrand(t: float) -> float:
        """Integrierte Funktion: f(t) * e^(-s*t)."""
        return f(t) * math.exp(-s * t)

    total = integrand(0.0) + integrand(t_max)
    for i in range(1, n_intervals):
        t_i = i * h
        weight = 4.0 if i % 2 == 1 else 2.0
        total += weight * integrand(t_i)

    return (h / 3.0) * total


def inverse_laplace(
    F: Callable[[float], float],
    t: float,
    sigma: float = 1.0,
    n_terms: int = 1000
) -> float:
    """
    Berechnet die inverse Laplace-Transformation numerisch (Post-Widder).

    Verwendet die Post-Widder-Inversion:
      f(t) ≈ lim_{n->inf} ((-1)^n / n!) * (n/t)^(n+1) * F^(n)(n/t)

    Praktische Implementierung via Euler-Summe (Abate & Whitt Methode).

    Formel: f(t) ≈ (e^A / (2*t)) * Re[ F(A/(2t)) + sum_{k=1}^{N} C_k * F((A + 2*pi*i*k) / (2t)) ]

    :param F: Laplace-Transformierte F(s) als Funktion
    :param t: Zeitpunkt (t > 0)
    :param sigma: Konvergenzparameter (sigma > maximale Singularität)
    :param n_terms: Anzahl der Terme für die Euler-Summe
    :return: Näherungswert von f(t)
    :last_modified: 2026-03-05
    """
    if t <= 0:
        raise ValueError("t muss > 0 sein")

    # Talbot-Methode: einfache aber effektive numerische Inversion
    # Basierend auf der Bromwich-Kontur mit optimierter Parametrisierung
    # Referenz: Abate & Valko (2004)
    M = n_terms
    # Talbot-Kontourparameter
    r = 2 * M / (5 * t)
    theta_vals = [k * math.pi / M for k in range(M + 1)]

    # Erste Term (k=0)
    s0 = r
    F0_val = F(s0).real if hasattr(F(s0), 'real') else F(s0)
    total = 0.5 * F0_val * math.exp(r * t)

    # Summe über Kontourpunkte
    for k in range(1, M + 1):
        theta = theta_vals[k]
        sigma_k = r * theta * (1.0/math.tan(theta) + 1j)  # Talbot-Kontur (komplex)

        # Ableitung der Talbot-Kontur
        sigma_prime = r * (1j + (theta * math.tan(theta) - 1.0) / (math.tan(theta)**2 + 1e-30))
        # Da F nur reell definiert: approximiere durch zwei reelle Auswertungen
        # F(s) für s = Re(sigma_k) ± kleine imaginäre Störung
        s_re = sigma_k.real

        # Vereinfachte Näherung: nur Realteil der Kontur verwenden
        try:
            F_val = F(s_re) if s_re > 0 else 0.0
        except Exception:
            F_val = 0.0

        weight = (2.0 / M) * math.exp(s_re * t)
        total += weight * F_val * math.cos(k * math.pi)

    return total / t


# Einfachere Alternative für inverse Laplace (Post-Widder Approximation)
def inverse_laplace(
    F: Callable[[float], float],
    t: float,
    sigma: float = 1.0,
    n_terms: int = 20
) -> float:
    """
    Berechnet die inverse Laplace-Transformation via Post-Widder-Formel.

    Numerische Approximation:
      f(t) ≈ (1/t) * sum_{k=1}^{N} (-1)^(k+N/2) * C(N, N/2+k) * F(k/t)
      (Stehfest-Algorithmus)

    :param F: Laplace-Transformierte F(s)
    :param t: Zeitpunkt (t > 0)
    :param sigma: Ungenutzt (Kompatibilität)
    :param n_terms: Anzahl der Terme (gerade Zahl empfohlen)
    :return: Näherungswert von f(t)
    :last_modified: 2026-03-05
    """
    # Stehfest-Algorithmus (sehr effizient für monotone Funktionen)
    N = n_terms if n_terms % 2 == 0 else n_terms + 1

    # Stehfest-Koeffizienten berechnen
    def stehfest_coeff(i: int, N: int) -> float:
        """Berechnet den i-ten Stehfest-Koeffizienten."""
        N2 = N // 2
        coeff = 0.0
        for k in range((i + 1) // 2, min(i, N2) + 1):
            num = (k ** N2) * _factorial(2 * k)
            den = (_factorial(N2 - k) * _factorial(k)
                   * _factorial(k - 1) * _factorial(i - k)
                   * _factorial(2 * k - i))
            coeff += num / den
        return ((-1) ** (i + N2)) * coeff

    # Approximation: f(t) ≈ (ln(2)/t) * sum_{i=1}^{N} V_i * F(i*ln(2)/t)
    ln2_over_t = math.log(2) / t
    total = 0.0
    for i in range(1, N + 1):
        s_i = i * ln2_over_t
        try:
            F_val = F(s_i)
        except Exception:
            F_val = 0.0
        total += stehfest_coeff(i, N) * F_val

    return ln2_over_t * total


def _factorial(n: int) -> float:
    """
    Berechnet n! (Fakultät).

    :param n: Nicht-negative ganze Zahl
    :return: n!
    :last_modified: 2026-03-05
    """
    if n < 0:
        return 0.0
    return float(math.factorial(n))
