"""
@file visualization.py
@brief Visualisierungsmodul für mathematische Funktionen und Strukturen.
@description
    Implementiert 2D/3D-Funktionsplotter, Vektorfelddarstellung,
    Phasenraum-Visualisierung für ODEs und Fraktal-Generator.

    Nutzt matplotlib für alle Visualisierungen.

    WICHTIG: Alle Plot-Funktionen haben einen Parameter save_path (Optional).
    Falls save_path angegeben, wird das Bild gespeichert statt angezeigt.
    Falls save_path=None, wird plt.show() aufgerufen.

    Unterstützte Visualisierungen:
    - 2D-Funktionsplotter (einzel, mehrere Funktionen, parametrisch)
    - 3D-Funktionsplotter (Surface, Konturplot, parametrische 3D-Kurven)
    - Vektorfelddarstellung (Quiver, Stromlinien)
    - Phasenraum-Visualisierung (Phasenportrait, Bifurkationsdiagramm)
    - Fraktal-Generator (Mandelbrot, Julia, Sierpinski, Newton)

@author Kurt Ingwer
@date 2026-03-08
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless-Modus (kein Display nötig)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  – wird für 3D-Plots gebraucht
from typing import Callable, Optional, List, Tuple
import math


# ===========================================================================
# HILFSFUNKTION: Plot speichern oder anzeigen
# ===========================================================================

def _save_or_show(fig: plt.Figure, save_path: Optional[str]) -> None:
    """
    @brief Speichert die Figur oder zeigt sie an, je nach save_path.
    @param fig: matplotlib Figure-Objekt
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    if save_path is not None:
        # Bild in Datei speichern (hohe Auflösung)
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    else:
        # Interaktiv anzeigen (braucht Display)
        plt.show()
    plt.close(fig)  # Speicher freigeben


# ===========================================================================
# 2D-FUNKTIONSPLOTTER
# ===========================================================================

def plot_function_2d(
    f: Callable,
    x_min: float,
    x_max: float,
    n_points: int = 500,
    title: str = '',
    xlabel: str = 'x',
    ylabel: str = 'y',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet eine einzelne 2D-Funktion y = f(x).
    @description
        Erzeugt einen einfachen Liniengraph der Funktion f über dem Intervall
        [x_min, x_max] mit n_points äquidistanten Auswertungspunkten.

        Punkte mit numerisch ungültigen Werten (NaN, Inf) werden
        automatisch ausgeblendet (z.B. bei Polstellen).

    @param f: Funktion f(x) → float (auswertbar für numpy-Arrays)
    @param x_min: Linke Intervallgrenze
    @param x_max: Rechte Intervallgrenze
    @param n_points: Anzahl Auswertungspunkte (Standard: 500)
    @param title: Titel des Diagramms
    @param xlabel: Beschriftung der x-Achse
    @param ylabel: Beschriftung der y-Achse
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # x-Werte gleichmäßig verteilen
    x = np.linspace(x_min, x_max, n_points)

    # Funktion auswerten, ungültige Werte (Polstellen etc.) maskieren
    y = np.array([f(xi) for xi in x], dtype=float)
    # NaN und Inf ausblenden
    mask = np.isfinite(y)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x[mask], y[mask], 'b-', linewidth=1.5)
    ax.set_title(title or f'f(x)', fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)  # x-Achse einzeichnen
    ax.axvline(0, color='k', linewidth=0.5)  # y-Achse einzeichnen

    _save_or_show(fig, save_path)


def plot_functions_2d(
    functions: List[Tuple[Callable, str]],
    x_min: float,
    x_max: float,
    n_points: int = 500,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet mehrere 2D-Funktionen in einem Diagramm.
    @description
        Erzeugt einen Liniengraphen mit mehreren Funktionen übereinander.
        Jede Funktion bekommt eine eigene Farbe und einen Legendeneintrag.

        Beispiel:
            functions = [(math.sin, 'sin(x)'), (math.cos, 'cos(x)')]

    @param functions: Liste von (Funktion, Label)-Tupeln
    @param x_min: Linke Intervallgrenze
    @param x_max: Rechte Intervallgrenze
    @param n_points: Anzahl Auswertungspunkte
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    x = np.linspace(x_min, x_max, n_points)

    fig, ax = plt.subplots(figsize=(9, 5))

    # Jede Funktion einzeln plotten
    for f, label in functions:
        y = np.array([f(xi) for xi in x], dtype=float)
        mask = np.isfinite(y)  # Polstellen etc. ausblenden
        ax.plot(x[mask], y[mask], linewidth=1.5, label=label)

    ax.set_title(title or 'Funktionenvergleich', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    _save_or_show(fig, save_path)


def plot_parametric_2d(
    x_func: Callable,
    y_func: Callable,
    t_min: float,
    t_max: float,
    n_points: int = 500,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet eine parametrische 2D-Kurve (x(t), y(t)).
    @description
        Eine parametrische Kurve wird durch zwei Funktionen beschrieben:
            x = x(t),  y = y(t),   t ∈ [t_min, t_max]

        Beispiel (Einheitskreis):
            x_func = math.cos,  y_func = math.sin,  t_min=0, t_max=2π

        Der Parameter t wird gleichmäßig über [t_min, t_max] verteilt,
        die Kurve durch Verbinden der erzeugten Punkte dargestellt.

    @param x_func: Funktion x(t)
    @param y_func: Funktion y(t)
    @param t_min: Startparameter
    @param t_max: Endparameter
    @param n_points: Anzahl Punkte
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    t = np.linspace(t_min, t_max, n_points)
    x = np.array([x_func(ti) for ti in t], dtype=float)
    y = np.array([y_func(ti) for ti in t], dtype=float)

    # Nur gültige Punkte plotten
    mask = np.isfinite(x) & np.isfinite(y)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(x[mask], y[mask], 'b-', linewidth=1.5)
    ax.set_title(title or 'Parametrische Kurve', fontsize=14)
    ax.set_xlabel('x(t)')
    ax.set_ylabel('y(t)')
    ax.set_aspect('equal')  # Gleiches Seitenverhältnis für Kreise etc.
    ax.grid(True, alpha=0.3)

    _save_or_show(fig, save_path)


# ===========================================================================
# 3D-FUNKTIONSPLOTTER
# ===========================================================================

def plot_function_3d(
    f: Callable,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    n_points: int = 50,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet eine 3D-Funktion z = f(x, y) als Surface-Plot.
    @description
        Erstellt ein Gitter aus n_points × n_points Punkten im Bereich
        [x_min, x_max] × [y_min, y_max] und berechnet dort z = f(x, y).

        Darstellung als gefärbte 3D-Oberfläche (Surface-Plot).
        Farbgebung entspricht der Höhe z (via colormap 'viridis').

    @param f: Funktion f(x, y) → float
    @param x_min: Untere x-Grenze
    @param x_max: Obere x-Grenze
    @param y_min: Untere y-Grenze
    @param y_max: Obere y-Grenze
    @param n_points: Gitterpunkte pro Achse (Standard: 50)
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Gitter erzeugen
    x = np.linspace(x_min, x_max, n_points)
    y = np.linspace(y_min, y_max, n_points)
    X, Y = np.meshgrid(x, y)

    # Funktion auswerten (sicher mit Fehlerbehandlung)
    Z = np.zeros_like(X)
    for i in range(n_points):
        for j in range(n_points):
            try:
                val = f(X[i, j], Y[i, j])
                Z[i, j] = val if math.isfinite(val) else np.nan
            except Exception:
                Z[i, j] = np.nan

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.85)
    fig.colorbar(surf, ax=ax, shrink=0.5)  # Farbskala
    ax.set_title(title or 'z = f(x, y)', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    _save_or_show(fig, save_path)


def plot_contour(
    f: Callable,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    levels: int = 20,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Konturplot (Höhenlinien) von f(x, y).
    @description
        Ein Konturplot zeigt Linien konstanter Funktionswerte (Isolinien):
            {(x,y) : f(x,y) = c}  für verschiedene c-Werte

        Nützlich um:
        - Minima/Maxima zu erkennen (konzentrische Kreise)
        - Sattelpunkte zu finden (Sattelform der Isolinien)
        - Gradientenrichtung abzuschätzen (senkrecht zu Isolinien)

    @param f: Funktion f(x, y) → float
    @param x_min: Untere x-Grenze
    @param x_max: Obere x-Grenze
    @param y_min: Untere y-Grenze
    @param y_max: Obere y-Grenze
    @param levels: Anzahl Höhenlinien (Standard: 20)
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    n = 100  # Feineres Gitter für Konturplot
    x = np.linspace(x_min, x_max, n)
    y = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(x, y)

    # Funktion über das Gitter auswerten
    Z = np.zeros_like(X)
    for i in range(n):
        for j in range(n):
            try:
                val = f(X[i, j], Y[i, j])
                Z[i, j] = val if math.isfinite(val) else np.nan
            except Exception:
                Z[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(8, 6))
    # Gefüllter Konturplot (farbig) + Linien mit Beschriftung
    cp = ax.contourf(X, Y, Z, levels=levels, cmap='RdYlBu_r')
    fig.colorbar(cp, ax=ax)
    ax.contour(X, Y, Z, levels=levels, colors='k', alpha=0.3, linewidths=0.5)
    ax.set_title(title or 'Konturplot f(x, y)', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    _save_or_show(fig, save_path)


def plot_parametric_3d(
    x_func: Callable,
    y_func: Callable,
    z_func: Callable,
    t_min: float,
    t_max: float,
    n_points: int = 500,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet eine 3D-Parameterkurve (x(t), y(t), z(t)).
    @description
        Eine Raumkurve wird durch drei Parameterfunktionen beschrieben:
            x = x(t),  y = y(t),  z = z(t),   t ∈ [t_min, t_max]

        Beispiel (Helix/Schraubenlinie):
            x(t) = cos(t),  y(t) = sin(t),  z(t) = t/(2π)

    @param x_func: Funktion x(t)
    @param y_func: Funktion y(t)
    @param z_func: Funktion z(t)
    @param t_min: Startparameter
    @param t_max: Endparameter
    @param n_points: Anzahl Punkte
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    t = np.linspace(t_min, t_max, n_points)
    x = np.array([x_func(ti) for ti in t], dtype=float)
    y = np.array([y_func(ti) for ti in t], dtype=float)
    z = np.array([z_func(ti) for ti in t], dtype=float)

    # Nur gültige Punkte
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x[mask], y[mask], z[mask], 'b-', linewidth=1.5)
    ax.set_title(title or '3D-Parameterkurve', fontsize=14)
    ax.set_xlabel('x(t)')
    ax.set_ylabel('y(t)')
    ax.set_zlabel('z(t)')

    _save_or_show(fig, save_path)


# ===========================================================================
# VEKTORFELDDARSTELLUNG
# ===========================================================================

def plot_vector_field_2d(
    fx: Callable,
    fy: Callable,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    n_grid: int = 20,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet ein 2D-Vektorfeld via Quiver-Plot.
    @description
        Ein 2D-Vektorfeld F(x,y) = (fx(x,y), fy(x,y)) wird als
        Pfeildiagramm dargestellt. An jedem Gitterpunkt wird ein
        Pfeil in Richtung und Länge des Feldvektors gezeichnet.

        Die Pfeillängen werden normiert (unit vectors) damit das
        Diagramm übersichtlich bleibt. Die Farbe codiert die Stärke.

    @param fx: x-Komponente des Vektorfelds F_x(x, y)
    @param fy: y-Komponente des Vektorfelds F_y(x, y)
    @param x_min: Untere x-Grenze
    @param x_max: Obere x-Grenze
    @param y_min: Untere y-Grenze
    @param y_max: Obere y-Grenze
    @param n_grid: Anzahl Gitterpunkte pro Achse (Standard: 20)
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Gitter erstellen
    x = np.linspace(x_min, x_max, n_grid)
    y = np.linspace(y_min, y_max, n_grid)
    X, Y = np.meshgrid(x, y)

    # Vektorfeldkomponenten auswerten
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(n_grid):
        for j in range(n_grid):
            try:
                U[i, j] = fx(X[i, j], Y[i, j])
                V[i, j] = fy(X[i, j], Y[i, j])
            except Exception:
                U[i, j] = 0.0
                V[i, j] = 0.0

    # Betrag (Länge) des Vektors für Farbkodierung
    magnitude = np.sqrt(U**2 + V**2)
    # Normieren (Pfeile einheitlicher Länge) – Nullvektoren sicher behandeln
    with np.errstate(invalid='ignore', divide='ignore'):
        U_norm = np.where(magnitude > 0, U / magnitude, 0.0)
        V_norm = np.where(magnitude > 0, V / magnitude, 0.0)

    fig, ax = plt.subplots(figsize=(8, 6))
    qv = ax.quiver(X, Y, U_norm, V_norm, magnitude, cmap='plasma', scale=n_grid)
    fig.colorbar(qv, ax=ax, label='|F|')
    ax.set_title(title or '2D-Vektorfeld', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    _save_or_show(fig, save_path)


def plot_stream_lines(
    fx: Callable,
    fy: Callable,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    density: float = 1.0,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Plottet Stromlinien eines 2D-Vektorfelds via streamplot.
    @description
        Stromlinien sind Kurven, die überall tangential zum Vektorfeld
        F(x,y) = (fx(x,y), fy(x,y)) verlaufen.

        Im Gegensatz zu Quiver-Plots zeigen Stromlinien den "Fluss"
        des Feldes anschaulich als Kurven.

        Der density-Parameter steuert die Dichte der Stromlinien
        (1.0 = Standard, größer = dichter).

    @param fx: x-Komponente F_x(x, y)
    @param fy: y-Komponente F_y(x, y)
    @param x_min: Untere x-Grenze
    @param x_max: Obere x-Grenze
    @param y_min: Untere y-Grenze
    @param y_max: Obere y-Grenze
    @param density: Dichte der Stromlinien (Standard: 1.0)
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Feines Gitter für Stromlinien
    n = 30
    x = np.linspace(x_min, x_max, n)
    y = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(x, y)

    # Vektorfeldkomponenten auswerten
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(n):
        for j in range(n):
            try:
                U[i, j] = fx(X[i, j], Y[i, j])
                V[i, j] = fy(X[i, j], Y[i, j])
            except Exception:
                U[i, j] = 0.0
                V[i, j] = 0.0

    # Betrag für Farbkodierung
    speed = np.sqrt(U**2 + V**2)

    fig, ax = plt.subplots(figsize=(8, 6))
    # streamplot zeichnet die Trajektorien des Vektorfelds
    strm = ax.streamplot(x, y, U, V, color=speed, cmap='cool',
                         density=density, linewidth=1.5)
    fig.colorbar(strm.lines, ax=ax, label='|F|')
    ax.set_title(title or 'Stromlinien', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    _save_or_show(fig, save_path)


# ===========================================================================
# PHASENRAUM-VISUALISIERUNG (für ODEs)
# ===========================================================================

def plot_phase_portrait(
    dx_dt: Callable,
    dy_dt: Callable,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    initial_conditions: Optional[List[Tuple[float, float]]] = None,
    t_max: float = 10.0,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Phasenporträt eines 2D-autonomen Differentialgleichungssystems.
    @description
        Das System lautet:
            dx/dt = f(x, y)
            dy/dt = g(x, y)

        Das Phasenporträt zeigt:
        1. Vektorfeld (Richtungsfeld): Gibt die Richtung des Flusses an
        2. Trajektorien: Lösungskurven für gegebene Anfangsbedingungen

        Wichtige Strukturen im Phasenraum:
        - Stabile Fixpunkte: alle Trajektorien laufen darauf zu (Attraktor)
        - Instabile Fixpunkte: alle Trajektorien fliehen davon
        - Grenzzyklen: periodische Trajektorien
        - Sattelpunkte: stabil in einer, instabil in anderer Richtung

        Trajektorien werden mit einfachem Euler-Verfahren berechnet.

    @param dx_dt: Rechte Seite dx/dt = f(x, y)
    @param dy_dt: Rechte Seite dy/dt = g(x, y)
    @param x_min: Linke x-Grenze des Phasenraums
    @param x_max: Rechte x-Grenze
    @param y_min: Untere y-Grenze
    @param y_max: Obere y-Grenze
    @param initial_conditions: Liste von (x0, y0) Anfangsbedingungen
    @param t_max: Maximale Simulationszeit für Trajektorien
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Vektorfeld-Gitter
    n_grid = 20
    x = np.linspace(x_min, x_max, n_grid)
    y = np.linspace(y_min, y_max, n_grid)
    X, Y = np.meshgrid(x, y)

    # Vektorfeldkomponenten auswerten
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(n_grid):
        for j in range(n_grid):
            try:
                U[i, j] = dx_dt(X[i, j], Y[i, j])
                V[i, j] = dy_dt(X[i, j], Y[i, j])
            except Exception:
                U[i, j] = 0.0
                V[i, j] = 0.0

    # Normieren für gleichmäßige Pfeildarstellung
    magnitude = np.sqrt(U**2 + V**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        U_norm = np.where(magnitude > 0, U / magnitude, 0.0)
        V_norm = np.where(magnitude > 0, V / magnitude, 0.0)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Vektorfeld plotten
    ax.quiver(X, Y, U_norm, V_norm, magnitude, cmap='Blues', alpha=0.6,
              scale=n_grid * 1.5)

    # Trajektorien berechnen und plotten (Euler-Verfahren)
    if initial_conditions:
        n_steps = 2000          # Anzahl Euler-Schritte
        dt = t_max / n_steps    # Zeitschritt

        for x0, y0 in initial_conditions:
            # Euler-Integration vorwärts
            traj_x = [x0]
            traj_y = [y0]
            cx, cy = float(x0), float(y0)

            for _ in range(n_steps):
                try:
                    # Euler-Schritt: x_{n+1} = x_n + dt * f(x_n, y_n)
                    fx_val = dx_dt(cx, cy)
                    fy_val = dy_dt(cx, cy)
                    cx += dt * fx_val
                    cy += dt * fy_val
                    # Abbruch wenn Trajektorie den Bereich verlässt
                    if not (x_min <= cx <= x_max and y_min <= cy <= y_max):
                        break
                    traj_x.append(cx)
                    traj_y.append(cy)
                except Exception:
                    break

            ax.plot(traj_x, traj_y, 'r-', linewidth=1.5, alpha=0.8)
            # Anfangspunkt markieren
            ax.plot(x0, y0, 'go', markersize=6)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(title or 'Phasenporträt', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    _save_or_show(fig, save_path)


def plot_bifurcation_diagram(
    f: Callable,
    r_min: float,
    r_max: float,
    n_r: int = 500,
    n_iter: int = 1000,
    n_discard: int = 500,
    title: str = '',
    save_path: Optional[str] = None
) -> None:
    """
    @brief Bifurkationsdiagramm für eine parametrische Abbildung f(x, r).
    @description
        Das Bifurkationsdiagramm zeigt, wie das Langzeitverhalten einer
        diskreten dynamischen Systems x_{n+1} = f(x_n, r) vom Parameter r abhängt.

        Klassisches Beispiel: Logistische Abbildung f(x, r) = r * x * (1 - x)
        - Für r < 3.0: Konvergenz zu einem Fixpunkt
        - Für 3.0 < r < 3.449: 2er-Periode
        - Für 3.449 < r < 3.544: 4er-Periode
        - Für r > 3.569: Chaos (Feigenbaum-Szenario)

        Algorithmus:
        1. Für jeden r-Wert: Starte mit x₀ = 0.5
        2. Iteriere n_iter + n_discard Mal
        3. Verwerfe die ersten n_discard Iterationen (Einschwingvorgang)
        4. Plotte die verbleibenden x-Werte als Punkte über r

    @param f: Abbildung f(x, r) → float
    @param r_min: Untere Parameterschranke
    @param r_max: Obere Parameterschranke
    @param n_r: Anzahl r-Werte (horizontale Auflösung)
    @param n_iter: Gesamtanzahl Iterationen pro r-Wert
    @param n_discard: Anzahl zu verwerfender Anfangsiterationen
    @param title: Diagrammtitel
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    r_values = np.linspace(r_min, r_max, n_r)

    # Alle r-Werte und zugehörige x-Attraktoren sammeln
    r_plot = []  # r-Werte für den Plot
    x_plot = []  # Attraktoren-x-Werte

    for r in r_values:
        x = 0.5  # Startwert (muss im Definitionsbereich liegen)

        # Einschwingvorgang: Erste n_discard Iterationen verwerfen
        for _ in range(n_discard):
            try:
                x = f(x, r)
                if not math.isfinite(x):
                    break
            except Exception:
                break

        # Attraktoren aufzeichnen
        for _ in range(n_iter - n_discard):
            try:
                x = f(x, r)
                if not math.isfinite(x):
                    break
                r_plot.append(r)
                x_plot.append(x)
            except Exception:
                break

    fig, ax = plt.subplots(figsize=(12, 7))
    # Sehr kleine Punkte für dichte Darstellung
    ax.plot(r_plot, x_plot, ',k', alpha=0.3, markersize=0.5)
    ax.set_title(title or 'Bifurkationsdiagramm', fontsize=14)
    ax.set_xlabel('Parameter r')
    ax.set_ylabel('x (Attraktor)')
    ax.grid(True, alpha=0.2)

    _save_or_show(fig, save_path)


# ===========================================================================
# FRAKTAL-GENERATOR
# ===========================================================================

def mandelbrot_set(
    x_min: float = -2.5,
    x_max: float = 1.0,
    y_min: float = -1.25,
    y_max: float = 1.25,
    width: int = 800,
    height: int = 500,
    max_iter: int = 100,
    save_path: Optional[str] = None
) -> np.ndarray:
    """
    @brief Berechnet und plottet die Mandelbrot-Menge.
    @description
        Die Mandelbrot-Menge M ist definiert als die Menge aller c ∈ ℂ,
        für die die Folge z_{n+1} = z_n² + c (mit z_0 = 0) beschränkt bleibt.

        Jeder Pixel entspricht einem komplexen Wert c.
        Die Iteration läuft, bis entweder:
        - |z_n| > 2 (Folge divergiert → Pixel außerhalb M, gefärbt nach Iterationszahl)
        - max_iter Iterationen erreicht (Pixel in M, schwarz)

        Die Grenze von M ist ein Fraktal mit unendlicher Komplexität.

    @param x_min: Realteil-Minimum der komplexen Ebene
    @param x_max: Realteil-Maximum
    @param y_min: Imaginärteil-Minimum
    @param y_max: Imaginärteil-Maximum
    @param width: Pixelbreite des Bildes
    @param height: Pixelhöhe des Bildes
    @param max_iter: Maximale Iterationszahl
    @param save_path: Pfad zum Speichern (None → plt.show())
    @return: 2D numpy-Array der Iterationszahlen (Höhe × Breite)
    @date 2026-03-08
    """
    # Komplexe Ebene diskretisieren
    x = np.linspace(x_min, x_max, width)
    y = np.linspace(y_min, y_max, height)
    C = x[np.newaxis, :] + 1j * y[:, np.newaxis]  # (height × width) komplexe Werte

    # Iterationszähler-Array
    iterations = np.zeros(C.shape, dtype=int)

    # Vektorisierte Mandelbrot-Iteration
    Z = np.zeros_like(C)           # Startwert z_0 = 0
    not_diverged = np.ones(C.shape, dtype=bool)  # Maske: noch nicht divergiert

    for i in range(max_iter):
        # z_{n+1} = z_n^2 + c (nur für noch nicht divergierte Punkte)
        Z[not_diverged] = Z[not_diverged]**2 + C[not_diverged]
        # Divergenzcheck: |z| > 2
        newly_diverged = not_diverged & (np.abs(Z) > 2.0)
        iterations[newly_diverged] = i + 1  # Iterationszahl merken
        not_diverged[newly_diverged] = False

    # Visualisierung
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.imshow(
        iterations,
        extent=[x_min, x_max, y_min, y_max],
        origin='lower',
        cmap='inferno',
        interpolation='bilinear'
    )
    ax.set_title('Mandelbrot-Menge', fontsize=14)
    ax.set_xlabel('Re(c)')
    ax.set_ylabel('Im(c)')

    _save_or_show(fig, save_path)
    return iterations


def julia_set(
    c: complex = -0.7 + 0.27j,
    x_min: float = -1.5,
    x_max: float = 1.5,
    y_min: float = -1.5,
    y_max: float = 1.5,
    width: int = 600,
    height: int = 600,
    max_iter: int = 100,
    save_path: Optional[str] = None
) -> np.ndarray:
    """
    @brief Berechnet und plottet die Julia-Menge für ein festes c.
    @description
        Die Julia-Menge J_c ist definiert als die Menge aller Startpunkte z_0 ∈ ℂ,
        für die die Folge z_{n+1} = z_n² + c (mit festem c) beschränkt bleibt.

        Im Gegensatz zur Mandelbrot-Menge (Startpunkt z_0=0, c variabel)
        wird hier c fest gewählt und der Startpunkt z_0 variiert.

        Interessante Werte von c:
        - c = -0.7 + 0.27j → gebundenes, zusammenhängendes Fraktal
        - c = 0.355 + 0.355j → dendritische Struktur
        - c = -1.417 + 0.107j → Kochkurven-ähnlich

    @param c: Fester komplexer Parameter c
    @param x_min: Realteil-Minimum des Startpunkts z_0
    @param x_max: Realteil-Maximum
    @param y_min: Imaginärteil-Minimum
    @param y_max: Imaginärteil-Maximum
    @param width: Pixelbreite
    @param height: Pixelhöhe
    @param max_iter: Maximale Iterationszahl
    @param save_path: Pfad zum Speichern (None → plt.show())
    @return: 2D numpy-Array der Iterationszahlen
    @date 2026-03-08
    """
    # Gitter der Startpunkte z_0
    x = np.linspace(x_min, x_max, width)
    y = np.linspace(y_min, y_max, height)
    Z = x[np.newaxis, :] + 1j * y[:, np.newaxis]  # Startpunkte (height × width)

    # Iterationszähler
    iterations = np.zeros(Z.shape, dtype=int)
    not_diverged = np.ones(Z.shape, dtype=bool)

    for i in range(max_iter):
        # z_{n+1} = z_n^2 + c (nur nicht-divergierte Punkte)
        Z[not_diverged] = Z[not_diverged]**2 + c
        # Divergenzcheck
        newly_diverged = not_diverged & (np.abs(Z) > 2.0)
        iterations[newly_diverged] = i + 1
        not_diverged[newly_diverged] = False

    # Visualisierung
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(
        iterations,
        extent=[x_min, x_max, y_min, y_max],
        origin='lower',
        cmap='magma',
        interpolation='bilinear'
    )
    ax.set_title(f'Julia-Menge (c = {c:.3f})', fontsize=14)
    ax.set_xlabel('Re(z₀)')
    ax.set_ylabel('Im(z₀)')

    _save_or_show(fig, save_path)
    return iterations


def sierpinski_triangle(
    n_iterations: int = 6,
    save_path: Optional[str] = None
) -> None:
    """
    @brief Sierpinski-Dreieck via Chaos-Spiel-Methode.
    @description
        Das Sierpinski-Dreieck ist ein selbstähnliches Fraktal.

        Chaos-Spiel-Methode (Barnsley):
        1. Definiere 3 Eckpunkte des Dreiecks: A, B, C
        2. Starte mit beliebigem Punkt P
        3. Wiederhole n_iterations mal:
           - Wähle zufällig eine Ecke E ∈ {A, B, C}
           - Setze P = Mittelpunkt von (P, E)
           - Zeichne P
        → Ergebnis konvergiert zum Sierpinski-Dreieck

        Für n_iterations = 100.000 entsteht ein sehr detailliertes Bild.

    @param n_iterations: Anzahl Punkte (mehr = feineres Fraktal, Standard: 6 → 6·10000=60000)
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Anzahl Punkte skalieren (n_iterations als Exponent-ähnlicher Parameter)
    n_points = n_iterations * 10000

    # Die 3 Eckpunkte des gleichseitigen Dreiecks
    vertices = np.array([
        [0.0, 0.0],    # A (links unten)
        [1.0, 0.0],    # B (rechts unten)
        [0.5, math.sqrt(3) / 2]  # C (oben, Höhe des gleichseitigen Dreiecks)
    ])

    # Zufälligen Startpunkt wählen
    rng = np.random.default_rng(42)  # Reproduzierbarer Zufallsgenerator
    points = np.zeros((n_points, 2))
    p = np.array([0.0, 0.0])  # Startpunkt

    for i in range(n_points):
        # Zufällige Ecke wählen (0, 1 oder 2)
        idx = rng.integers(0, 3)
        # Mittelpunkt zwischen aktuellem Punkt und gewählter Ecke
        p = (p + vertices[idx]) / 2.0
        points[i] = p

    fig, ax = plt.subplots(figsize=(8, 7))
    # Sehr kleine Punkte für feines Fraktal
    ax.scatter(points[:, 0], points[:, 1], s=0.01, c='darkblue', alpha=0.5)
    ax.set_title(f'Sierpinski-Dreieck ({n_points:,} Punkte)', fontsize=14)
    ax.set_aspect('equal')
    ax.axis('off')  # Achsen ausblenden für ästhetische Darstellung

    _save_or_show(fig, save_path)


def newton_fractal(
    f_coeffs: List[float],
    width: int = 400,
    height: int = 400,
    x_min: float = -2.0,
    x_max: float = 2.0,
    y_min: float = -2.0,
    y_max: float = 2.0,
    max_iter: int = 50,
    save_path: Optional[str] = None
) -> None:
    """
    @brief Newton-Fraktal für ein Polynom f.
    @description
        Das Newton-Fraktal entsteht durch Anwendung des Newton-Verfahrens
        auf ein komplexes Polynom f(z):
            z_{n+1} = z_n - f(z_n) / f'(z_n)

        Jeder Startpunkt z_0 (Pixel) konvergiert zu einer Nullstelle.
        Die Pixel werden nach der Nullstelle eingefärbt, zu der sie konvergieren.
        Je nach Startpunkt können sehr unterschiedliche Nullstellen erreicht werden,
        was zu den charakteristischen fraktalen Grenzen führt.

        Beispiel für x³ - 1 = 0 (f_coeffs = [1, 0, 0, -1]):
        Die 3 Nullstellen (3. Einheitswurzeln) erzeugen ein symmetrisches Fraktal.

    @param f_coeffs: Polynomkoeffizienten (höchster Grad zuerst)
                     Beispiel: [1, 0, 0, -1] → z³ - 1
    @param width: Pixelbreite
    @param height: Pixelhöhe
    @param x_min: Realteil-Minimum
    @param x_max: Realteil-Maximum
    @param y_min: Imaginärteil-Minimum
    @param y_max: Imaginärteil-Maximum
    @param max_iter: Maximale Iterationszahl des Newton-Verfahrens
    @param save_path: Pfad zum Speichern (None → plt.show())
    @date 2026-03-08
    """
    # Polynomkoeffizienten → Numpy poly1d Objekt
    p = np.poly1d(f_coeffs)
    dp = p.deriv()  # Ableitung des Polynoms (analytisch)

    # Nullstellen des Polynoms berechnen (für Farbzuordnung)
    roots = np.roots(f_coeffs)
    n_roots = len(roots)

    # Gitter der Startpunkte
    x = np.linspace(x_min, x_max, width)
    y = np.linspace(y_min, y_max, height)
    Z = x[np.newaxis, :] + 1j * y[:, np.newaxis]  # (height × width)

    # Farb-Array: speichert Index der konvergierten Nullstelle
    root_map = np.full(Z.shape, -1, dtype=int)    # -1 = nicht konvergiert
    iter_map = np.zeros(Z.shape, dtype=float)     # Iterationszahl (für Schattierung)

    # Toleranz für Konvergenzdetektion
    tol = 1e-6

    for i in range(max_iter):
        # Newton-Schritt: z = z - f(z)/f'(z)
        denom = dp(Z)
        # Nullteiler vermeiden
        safe_denom = np.where(np.abs(denom) < 1e-15, 1e-15, denom)
        Z = Z - p(Z) / safe_denom

        # Prüfen ob Z nahe einer Nullstelle ist
        for r_idx, root in enumerate(roots):
            # Maske: nicht konvergiert UND nahe dieser Nullstelle
            converged = (root_map == -1) & (np.abs(Z - root) < tol)
            root_map[converged] = r_idx
            iter_map[converged] = i + 1

    # Nicht konvergierte Pixel als letzten root markieren (Fallback)
    root_map[root_map == -1] = n_roots - 1 if n_roots > 0 else 0

    # Visualisierung: Farbe nach Nullstelle, Helligkeit nach Iterationszahl
    # Farb-Array erzeugen: Mischung aus root_index und iter_map für ästhetischen Effekt
    color_data = root_map.astype(float) / max(n_roots - 1, 1)  # Normiert auf [0,1]
    # Helligkeit modulieren (weniger Iterationen → heller)
    brightness = 1.0 - iter_map / (max_iter + 1)
    combined = color_data * 0.7 + brightness * 0.3

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(
        combined,
        extent=[x_min, x_max, y_min, y_max],
        origin='lower',
        cmap='hsv',
        interpolation='bilinear'
    )
    degree = len(f_coeffs) - 1
    ax.set_title(f'Newton-Fraktal (Grad {degree})', fontsize=14)
    ax.set_xlabel('Re(z₀)')
    ax.set_ylabel('Im(z₀)')

    _save_or_show(fig, save_path)
