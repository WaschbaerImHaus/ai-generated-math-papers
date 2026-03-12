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
    - Adaptiver Gitterplot (feineres Gitter nahe Nullstellen/Singularitäten)
    - Interaktiver Plot mit matplotlib.widgets.Slider

@author Michael Fuhrmann
@date 2026-03-08
@lastModified 2026-03-11
"""

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless-Modus (kein Display nötig)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  – wird für 3D-Plots gebraucht
from typing import Callable, Optional, List, Tuple


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


# ===========================================================================
# ODE-TRAJEKTORIE ANIMATIONEN
# ===========================================================================

import matplotlib.animation as animation  # noqa: E402  – nach matplotlib.use('Agg')
from scipy.integrate import solve_ivp     # noqa: E402
import os                                  # noqa: E402


def animate_ode_trajectory(
    f: Callable,
    y0: list,
    t_span: tuple,
    n_frames: int = 100,
    output_path: str = None,
    interval: int = 50,
    title: str = "ODE-Trajektorie"
) -> 'animation.FuncAnimation':
    """
    @brief Erstellt eine matplotlib-Animation einer ODE-Trajektorie.
    @description
        Löst eine gewöhnliche Differentialgleichung y' = f(t, y) und
        animiert die Lösung Frame für Frame.

        Wenn output_path angegeben:
        - .gif  → Pillow-Writer (immer verfügbar)
        - sonst → ffmpeg (falls nicht installiert: Warnung, kein Fehler)

        Unterstützt Systeme beliebiger Dimension (len(y0) bestimmt die
        Anzahl der dargestellten Komponenten).

    @param f: Funktion f(t, y) → dy/dt (rechte Seite der ODE)
    @param y0: Anfangswert-Vektor (Liste von Anfangswerten)
    @param t_span: (t_start, t_end) Zeitintervall
    @param n_frames: Anzahl Animationsframes (Standard: 100)
    @param output_path: Optionaler Pfad zum Speichern (.gif empfohlen)
    @param interval: Millisekunden zwischen Frames (Standard: 50)
    @param title: Animationstitel
    @return: FuncAnimation-Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # ODE lösen: t_eval liefert n_frames äquidistante Zeitpunkte
    t_eval = np.linspace(t_span[0], t_span[1], n_frames)
    sol = solve_ivp(f, t_span, y0, t_eval=t_eval, method='RK45', dense_output=True)

    t_vals = sol.t
    y_vals = sol.y  # Shape: (n_komponenten, n_frames)

    n_comp = y_vals.shape[0]  # Anzahl Zustandsvariablen

    # Figure und Achse erstellen
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('t')
    ax.set_ylabel('y(t)')

    # Achsengrenzen aus der vollständigen Lösung bestimmen
    y_min_val = float(y_vals.min())
    y_max_val = float(y_vals.max())
    y_margin = max(0.1 * (y_max_val - y_min_val), 0.1)
    ax.set_xlim(t_vals[0], t_vals[-1])
    ax.set_ylim(y_min_val - y_margin, y_max_val + y_margin)

    # Linienserie für jede Komponente anlegen
    farben = ['b', 'r', 'g', 'c', 'm', 'y']
    lines = []
    for i in range(n_comp):
        farbe = farben[i % len(farben)]
        linie, = ax.plot([], [], farbe + '-', linewidth=1.5, label=f'y{i}(t)')
        lines.append(linie)

    # Senkrechte Linie als Zeitmarkierung
    time_line = ax.axvline(x=t_vals[0], color='gray', linestyle='--', alpha=0.5)

    if n_comp > 1:
        ax.legend()

    ax.grid(True, alpha=0.3)

    def init():
        """Initialisierung: leere Linien, Zeitmarkierung am Startpunkt."""
        for linie in lines:
            linie.set_data([], [])
        time_line.set_xdata([t_vals[0]])
        return lines + [time_line]

    def update(frame):
        """Aktualisierung pro Frame: Trajektorie bis zum aktuellen Zeitpunkt."""
        idx = min(frame, len(t_vals) - 1)
        for i, linie in enumerate(lines):
            linie.set_data(t_vals[:idx + 1], y_vals[i, :idx + 1])
        time_line.set_xdata([t_vals[idx]])
        return lines + [time_line]

    # FuncAnimation erstellen
    anim = animation.FuncAnimation(
        fig, update,
        frames=n_frames,
        init_func=init,
        interval=interval,
        blit=True
    )

    # Speichern falls Pfad angegeben
    if output_path is not None:
        if output_path.endswith('.gif'):
            writer = animation.PillowWriter(fps=max(1, 1000 // interval))
            anim.save(output_path, writer=writer)
        else:
            try:
                anim.save(output_path, writer='ffmpeg')
            except Exception as e:
                import warnings
                warnings.warn(f"ffmpeg nicht verfügbar oder Fehler: {e}. GIF-Export verwenden.")

    plt.close(fig)
    return anim


def animate_phase_portrait(
    f: Callable,
    x_range: tuple,
    y_range: tuple,
    n_trajectories: int = 5,
    output_path: str = None,
    title: str = "Phasenporträt-Animation"
) -> 'animation.FuncAnimation':
    """
    @brief Animiertes Phasenporträt: Zeigt wie Trajektorien sich über Zeit entwickeln.
    @description
        Für ein 2D-System dx/dt = f(t, [x, y])[0], dy/dt = f(t, [x, y])[1]
        werden mehrere Trajektorien von verschiedenen Startpunkten aus
        gleichzeitig animiert.

        Hintergrund: statisches Quiver-Vektorfeld.
        Vordergrund: animierte Trajektorien.

        Startpunkte werden auf einem gleichmäßigen Gitter über den
        inneren 70% des Phasenraums verteilt.

    @param f: Funktion f(t, [x, y]) → [dx/dt, dy/dt]
    @param x_range: (x_min, x_max) Bereich der x-Achse
    @param y_range: (y_min, y_max) Bereich der y-Achse
    @param n_trajectories: Anzahl der Starttrajektorien (gerundet auf nächste Quadratzahl)
    @param output_path: Optionaler Speicherpfad (.gif)
    @param title: Animationstitel
    @return: FuncAnimation-Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    x_min, x_max = x_range
    y_min, y_max = y_range

    # Startpunkte: Gitter über innere 70% des Phasenraums
    n_sqrt = max(2, int(n_trajectories ** 0.5))
    xs = np.linspace(x_min * 0.7, x_max * 0.7, n_sqrt)
    ys = np.linspace(y_min * 0.7, y_max * 0.7, n_sqrt)
    start_points = [(x, y) for x in xs for y in ys][:n_trajectories]

    # Alle Trajektorien vorberechnen
    t_max_sim = 10.0
    n_frames = 80
    t_span = (0.0, t_max_sim)
    t_eval = np.linspace(0.0, t_max_sim, n_frames)

    trajs = []
    for x0, y0 in start_points:
        try:
            sol = solve_ivp(f, t_span, [x0, y0], t_eval=t_eval,
                            method='RK45', max_step=0.1)
            # Nur Punkte innerhalb des Plotbereichs
            mask = (
                (sol.y[0] >= x_min) & (sol.y[0] <= x_max) &
                (sol.y[1] >= y_min) & (sol.y[1] <= y_max)
            )
            if mask.any():
                # Ersten zusammenhängenden gültigen Block behalten
                end_idx = int(np.argmax(~mask)) if not mask.all() else len(mask)
                trajs.append((sol.t[:end_idx], sol.y[0, :end_idx], sol.y[1, :end_idx]))
            else:
                trajs.append((sol.t[:1], sol.y[0, :1], sol.y[1, :1]))
        except Exception:
            trajs.append((np.array([0.0]), np.array([x0]), np.array([y0])))

    # Figure erstellen
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True, alpha=0.3)

    # Hintergrund: statisches Vektorfeld
    n_grid = 15
    xg = np.linspace(x_min, x_max, n_grid)
    yg = np.linspace(y_min, y_max, n_grid)
    XG, YG = np.meshgrid(xg, yg)
    U_bg = np.zeros_like(XG)
    V_bg = np.zeros_like(YG)
    for i in range(n_grid):
        for j in range(n_grid):
            try:
                res = f(0.0, [XG[i, j], YG[i, j]])
                U_bg[i, j] = float(res[0])
                V_bg[i, j] = float(res[1])
            except Exception:
                pass
    mag = np.sqrt(U_bg**2 + V_bg**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        U_n = np.where(mag > 0, U_bg / mag, 0.0)
        V_n = np.where(mag > 0, V_bg / mag, 0.0)
    ax.quiver(XG, YG, U_n, V_n, mag, cmap='Blues', alpha=0.4, scale=n_grid * 1.5)

    # Dynamische Linien und Punkte für jede Trajektorie
    farben = plt.cm.Set1(np.linspace(0, 1, max(len(trajs), 1)))
    traj_lines = []
    traj_dots = []
    for i in range(len(trajs)):
        linie, = ax.plot([], [], '-', color=farben[i], linewidth=1.5, alpha=0.8)
        punkt, = ax.plot([], [], 'o', color=farben[i], markersize=6)
        traj_lines.append(linie)
        traj_dots.append(punkt)

    def init_pp():
        """Initialisierung: alle Linien und Punkte leer."""
        for linie in traj_lines:
            linie.set_data([], [])
        for punkt in traj_dots:
            punkt.set_data([], [])
        return traj_lines + traj_dots

    def update_pp(frame):
        """Aktualisierung: Trajektorien bis zum aktuellen Frame."""
        for i, (t_t, x_t, y_t) in enumerate(trajs):
            idx = min(frame, len(x_t) - 1)
            traj_lines[i].set_data(x_t[:idx + 1], y_t[:idx + 1])
            traj_dots[i].set_data([x_t[idx]], [y_t[idx]])
        return traj_lines + traj_dots

    anim = animation.FuncAnimation(
        fig, update_pp,
        frames=n_frames,
        init_func=init_pp,
        interval=60,
        blit=True
    )

    if output_path is not None:
        if output_path.endswith('.gif'):
            writer = animation.PillowWriter(fps=15)
            anim.save(output_path, writer=writer)
        else:
            try:
                anim.save(output_path, writer='ffmpeg')
            except Exception as e:
                import warnings
                warnings.warn(f"ffmpeg nicht verfügbar: {e}")

    plt.close(fig)
    return anim


def animate_wave_equation(
    k: float = 1.0,
    omega: float = 1.0,
    x_range: tuple = (-5, 5),
    t_max: float = 2 * math.pi,
    n_frames: int = 60,
    output_path: str = None
) -> 'animation.FuncAnimation':
    """
    @brief Animation der Wellengleichung u(x,t) = sin(k·x - ω·t).
    @description
        Zeigt die Ausbreitung einer harmonischen ebenen Welle.

        Die 1D-Wellengleichung (partielle DGL):
            ∂²u/∂t² = c² · ∂²u/∂x²  mit  c = ω/k (Phasengeschwindigkeit)

        Hat als Lösung u(x,t) = A · sin(k·x - ω·t + φ).

        Physikalische Bedeutung der Parameter:
            k = 2π/λ   (Wellenzahl, λ = Wellenlänge)
            ω = 2π/T   (Kreisfrequenz, T = Periodendauer)
            c = ω/k    (Phasengeschwindigkeit)

    @param k: Wellenzahl (räumliche Frequenz, Standard: 1.0)
    @param omega: Kreisfrequenz (zeitliche Frequenz, Standard: 1.0)
    @param x_range: (x_min, x_max) Raumbereich
    @param t_max: Animationsdauer in Zeiteinheiten (Standard: 2π)
    @param n_frames: Anzahl der Frames (Standard: 60)
    @param output_path: Optionaler Speicherpfad (.gif)
    @return: FuncAnimation-Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    x_min, x_max = x_range
    n_points = 500  # Räumliche Auflösung
    x = np.linspace(x_min, x_max, n_points)

    # Zeitpunkte für Animation
    t_values = np.linspace(0.0, t_max, n_frames)

    # Phasengeschwindigkeit c = ω/k
    c_phase = omega / k if k != 0.0 else 1.0

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-1.5, 1.5)  # Amplitude der Welle = 1
    ax.set_title(
        f'Wellengleichung: u(x,t) = sin({k:.1f}·x − {omega:.1f}·t)  '
        f'[Phasengeschw. c = {c_phase:.2f}]',
        fontsize=11
    )
    ax.set_xlabel('x (Ort)')
    ax.set_ylabel('u(x, t) (Auslenkung)')
    ax.axhline(0, color='k', linewidth=0.5)
    ax.grid(True, alpha=0.3)

    # Wellenlinie und Zeitanzeige
    linie, = ax.plot([], [], 'b-', linewidth=2.0)
    zeittext = ax.text(0.02, 0.90, '', transform=ax.transAxes, fontsize=11)

    def init_wave():
        """Initialisierung: leere Linie."""
        linie.set_data([], [])
        zeittext.set_text('')
        return linie, zeittext

    def update_wave(frame):
        """Frame-Update: harmonische Welle zum Zeitpunkt t_values[frame]."""
        t = t_values[frame]
        # Harmonische Welle: u(x,t) = sin(k·x - ω·t)
        u = np.sin(k * x - omega * t)
        linie.set_data(x, u)
        zeittext.set_text(f't = {t:.3f}')
        return linie, zeittext

    anim = animation.FuncAnimation(
        fig, update_wave,
        frames=n_frames,
        init_func=init_wave,
        interval=50,
        blit=True
    )

    if output_path is not None:
        if output_path.endswith('.gif'):
            writer = animation.PillowWriter(fps=20)
            anim.save(output_path, writer=writer)
        else:
            try:
                anim.save(output_path, writer='ffmpeg')
            except Exception as e:
                import warnings
                warnings.warn(f"ffmpeg nicht verfügbar: {e}")

    plt.close(fig)
    return anim


# ===========================================================================
# SVG/PDF-EXPORT
# ===========================================================================

def save_figure_svg(fig: plt.Figure, filepath: str, dpi: int = 300) -> str:
    """
    @brief Speichert eine matplotlib-Figure als SVG-Datei.
    @description
        SVG (Scalable Vector Graphics) ist ein verlustfreies Vektorformat,
        das in jedem Maßstab gestochen scharf bleibt.

        Eigenschaften:
        - XML-basiertes Format (lesbar und bearbeitbar)
        - Beliebig skalierbar ohne Qualitätsverlust
        - Kann in Inkscape, Adobe Illustrator etc. nachbearbeitet werden
        - Ideal für Webseiten und Präsentationen

        Der Parameter dpi betrifft nur eingebettete Rasterelemente.

    @param fig: matplotlib Figure-Objekt
    @param filepath: Ausgabepfad (sollte auf .svg enden)
    @param dpi: Auflösung für Rasterelemente innerhalb der SVG (Standard: 300)
    @return: Absoluter Pfad der gespeicherten Datei
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Absoluten Pfad auflösen
    abs_path = os.path.abspath(filepath)
    # Ausgabeverzeichnis anlegen falls nötig
    dir_name = os.path.dirname(abs_path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    # SVG speichern (bbox_inches='tight' schneidet überschüssige Ränder ab)
    fig.savefig(abs_path, format='svg', dpi=dpi, bbox_inches='tight')
    return abs_path


def save_figure_pdf(fig: plt.Figure, filepath: str) -> str:
    """
    @brief Speichert eine matplotlib-Figure als PDF-Datei.
    @description
        PDF ist das Standard-Format für wissenschaftliche Publikationen.

        Eigenschaften:
        - Matplotlib erzeugt echte Vektorgrafiken (keine Rasterisierung)
        - Gestochen scharf in jeder Druckgröße
        - Kann direkt in LaTeX eingebunden werden: \\includegraphics{...}
        - Unterstützt eingebettete Schriften (keine Zeichensatzprobleme)

    @param fig: matplotlib Figure-Objekt
    @param filepath: Ausgabepfad (sollte auf .pdf enden)
    @return: Absoluter Pfad der gespeicherten Datei
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    abs_path = os.path.abspath(filepath)
    dir_name = os.path.dirname(abs_path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    # PDF speichern (format='pdf' für echte Vektorgrafik ohne Rasterisierung)
    fig.savefig(abs_path, format='pdf', bbox_inches='tight')
    return abs_path


def plot_and_export(
    func: Callable,
    x_range: tuple,
    filepath: str,
    title: str = "",
    format: str = "svg"
) -> str:
    """
    @brief Kombiniert plot_function_2d() und Export in einem Schritt.
    @description
        Plottet eine Funktion und speichert sie direkt als SVG oder PDF,
        ohne interaktive Anzeige.

        Typische Verwendung:
            import math
            path = plot_and_export(
                math.sin, (-3.14, 3.14),
                '/tmp/sinus.svg', 'Sinusfunktion', 'svg'
            )

    @param func: Funktion f(x) → float (skalare Auswertung)
    @param x_range: (x_min, x_max) Plotintervall
    @param filepath: Ausgabepfad (Dateiendung wird durch format bestimmt)
    @param title: Diagrammtitel
    @param format: Exportformat: 'svg' oder 'pdf'
    @return: Absoluter Pfad der gespeicherten Datei
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # x-Werte erzeugen, Funktion auswerten, ungültige Werte maskieren
    x = np.linspace(x_range[0], x_range[1], 500)
    y = np.array([func(xi) for xi in x], dtype=float)
    mask = np.isfinite(y)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x[mask], y[mask], 'b-', linewidth=1.5)
    ax.set_title(title or 'f(x)', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    # Export je nach gewähltem Format
    if format.lower() == 'pdf':
        result = save_figure_pdf(fig, filepath)
    else:
        # Standardmäßig SVG
        result = save_figure_svg(fig, filepath)

    plt.close(fig)
    return result


def export_all_formats(fig: plt.Figure, base_path: str) -> dict:
    """
    @brief Exportiert eine Figure in alle Standard-Formate: PNG, SVG, PDF.
    @description
        Speichert dieselbe Figure in drei verschiedenen Formaten.

        Verwendungsszenarien:
        - PNG: Schnelle Vorschau, Web-Thumbnails, Rasterformat (150 dpi)
        - SVG: Skalierbare Web-Grafiken, Vektorgrafik-Bearbeitung
        - PDF: LaTeX-Integration, Druckvorbereitung, Publikationen

        Dateinamen werden automatisch erzeugt:
            base_path = '/pfad/figur'
            → Dateien: /pfad/figur.png, /pfad/figur.svg, /pfad/figur.pdf

    @param fig: matplotlib Figure-Objekt
    @param base_path: Basispfad OHNE Dateiendung
    @return: Dictionary mit Format als Schlüssel und absolutem Pfad als Wert:
             {'png': '/pfad/figur.png', 'svg': '/pfad/figur.svg', 'pdf': '/pfad/figur.pdf'}
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    abs_base = os.path.abspath(base_path)
    dir_name = os.path.dirname(abs_base)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)

    # PNG (Rasterformat, 150 dpi für gute Bildqualität)
    png_path = abs_base + '.png'
    fig.savefig(png_path, format='png', dpi=150, bbox_inches='tight')

    # SVG (Vektorformat)
    svg_path = save_figure_svg(fig, abs_base + '.svg')

    # PDF (Vektorformat für Publikationen)
    pdf_path = save_figure_pdf(fig, abs_base + '.pdf')

    return {
        'png': png_path,
        'svg': svg_path,
        'pdf': pdf_path
    }


# ===========================================================================
# MANDELBROT SMOOTH (vollständig NumPy-vektorisiert, Smooth Iteration Count)
# ===========================================================================

def mandelbrot_smooth(
    x_min: float = -2.5,
    x_max: float = 1.0,
    y_min: float = -1.25,
    y_max: float = 1.25,
    width: int = 800,
    height: int = 600,
    max_iter: int = 100
) -> np.ndarray:
    """
    @brief Mandelbrot-Menge mit Smooth Iteration Count für weiche Farbverläufe.
    @description
        Berechnet für jeden Pixel c = x + iy die kontinuierliche Escape-Zeit
        der Mandelbrot-Iteration z_{n+1} = z_n² + c (z_0 = 0).

        Standard-Coloring hat Treppeneffekte (ganzzahlige Iteration).
        Smooth Coloring berechnet einen kontinuierlichen Wert:

            n_smooth = n + 1 - log₂(log(|z_n|))

        Herleitung (Linas Vepstas, 2004):
            Wenn |z| > 2: tatsächliche Escape-Zeit ≈ n + 1 - log(log|z|)/log(2)
            Dieser Wert interpoliert glatt zwischen ganzzahligen Iterationen.

        Vollständig NumPy-vektorisiert: Keine Python-Schleifen über Pixel,
        nur eine äußere Schleife über max_iter Iterationsschritte.
        Durch Maskenoperationen werden nur aktive (noch nicht divergierte)
        Pixel pro Schritt aktualisiert.

    @param x_min: Realteil-Minimum der komplexen Ebene (Standard: -2.5)
    @param x_max: Realteil-Maximum (Standard: 1.0)
    @param y_min: Imaginärteil-Minimum (Standard: -1.25)
    @param y_max: Imaginärteil-Maximum (Standard: 1.25)
    @param width: Pixelbreite des Ergebnisarrays (Standard: 800)
    @param height: Pixelhöhe des Ergebnisarrays (Standard: 600)
    @param max_iter: Maximale Iterationszahl (Standard: 100)
    @return: 2D float-Array (height × width) mit Smooth-Escape-Werten.
             Punkte in M (beschränkt) haben Wert max_iter.
             Außerhalb liegende Pixel haben Werte in [0, max_iter).
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Komplexe Ebene diskretisieren: C[i,j] = (x_j) + i·(y_i)
    x = np.linspace(x_min, x_max, width)
    y = np.linspace(y_min, y_max, height)
    C = x[np.newaxis, :] + 1j * y[:, np.newaxis]  # (height × width)

    # Z: Iterationszustand, startet bei 0 für alle Pixel
    Z = np.zeros_like(C, dtype=complex)

    # smooth_count: Ergebnis-Array, initialisiert mit max_iter (= "in M")
    smooth_count = np.full(C.shape, float(max_iter), dtype=float)

    # not_diverged: Maske der noch aktiven Pixel (True = noch nicht divergiert)
    not_diverged = np.ones(C.shape, dtype=bool)

    for i in range(max_iter):
        # Iteration: z → z² + c (nur für noch aktive Pixel)
        Z[not_diverged] = Z[not_diverged] ** 2 + C[not_diverged]

        # Divergenzcheck: |z| > 2 (Escape-Radius = 2)
        newly_diverged = not_diverged & (np.abs(Z) > 2.0)

        if newly_diverged.any():
            # Smooth Iteration Count für neu divergierte Pixel
            abs_z = np.abs(Z[newly_diverged])
            # Numerisch sicher: abs_z > 2 garantiert, also log(abs_z) > log(2) > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                log_abs = np.log(np.maximum(abs_z, 1e-10))
                log_log_abs = np.log(np.maximum(log_abs, 1e-10))
            # n_smooth = (i+1) - log₂(log(|z|)) = (i+1) - log(log(|z|))/log(2)
            smooth_count[newly_diverged] = (i + 1) - log_log_abs / math.log(2)
            # Werte in [0, max_iter] clippen (numerische Ausreißer abfangen)
            smooth_count[newly_diverged] = np.clip(
                smooth_count[newly_diverged], 0.0, float(max_iter)
            )

        # Divergierte Pixel aus aktiver Menge entfernen
        not_diverged[newly_diverged] = False

        # Frühzeitig abbrechen wenn alle Pixel divergiert sind
        if not not_diverged.any():
            break

    return smooth_count


# ===========================================================================
# HILFSFUNKTION: Universeller Figure-Export (PNG/SVG/PDF)
# ===========================================================================

def export_figure(fig: plt.Figure, path: str) -> str:
    """
    @brief Exportiert eine matplotlib-Figure in das durch die Dateiendung
           bestimmte Format (PNG, SVG oder PDF).
    @description
        Erkennt das Zielformat automatisch anhand der Dateiendung:
            .png  → Rastergrafik (150 dpi)
            .svg  → Vektorgrafik (XML)
            .pdf  → PDF-Vektorgrafik für Publikationen

        Das Ausgabeverzeichnis wird automatisch angelegt, falls es nicht
        existiert. Gibt den absoluten Pfad der gespeicherten Datei zurück.

        Beispiel:
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])
            export_figure(fig, '/tmp/ausgabe.svg')

    @param fig: matplotlib Figure-Objekt, das exportiert werden soll
    @param path: Ausgabepfad inkl. Dateiname und Endung (.png/.svg/.pdf)
    @return: Absoluter Pfad der gespeicherten Datei
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    import os  # noqa: F811 – lokaler Import für Klarheit
    abs_path = os.path.abspath(path)
    # Verzeichnis anlegen falls nötig
    dir_name = os.path.dirname(abs_path)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    # Format aus Dateiendung bestimmen
    ext = os.path.splitext(abs_path)[1].lower().lstrip('.')
    if ext == 'pdf':
        # PDF: verlustfreies Vektorformat
        fig.savefig(abs_path, format='pdf', bbox_inches='tight')
    elif ext == 'svg':
        # SVG: XML-basiertes Vektorformat
        fig.savefig(abs_path, format='svg', dpi=300, bbox_inches='tight')
    else:
        # Standard: PNG mit hoher Auflösung (auch für unbekannte Endungen)
        fig.savefig(abs_path, format='png', dpi=150, bbox_inches='tight')
    return abs_path


# ===========================================================================
# 3D-VISUALISIERUNG: GAUẞSCHE KRÜMMUNG
# ===========================================================================

def plot_gaussian_curvature_3d(
    surface_name: str,
    param_range: Tuple[float, float] = (-2.0, 2.0),
    resolution: int = 50
) -> plt.Figure:
    """
    @brief Visualisiert die Gaußsche Krümmung K einer parametrischen Fläche
           als farbkodierten 3D-Surface-Plot.
    @description
        Die Gaußsche Krümmung K = κ₁ · κ₂ ist das Produkt der beiden
        Hauptkrümmungen einer Fläche:
            K > 0: elliptischer Punkt (z.B. Sphäre → K = 1/R²)
            K = 0: parabolischer Punkt (z.B. Zylinder)
            K < 0: hyperbolischer Punkt (z.B. Sattel)

        Unterstützte Flächen:
            'sphere':               Einheitssphäre, K = 1 überall
            'torus':                Torus (R=2, r=1), K variiert zwischen -1 und positiv
            'saddle':               z = x² - y², K = -4/(1 + 4x² + 4y²)²
            'hyperbolic_paraboloid': z = x·y, K = -1/(1+x²+y²)²

        Die Krümmung wird numerisch über die zweiten partiellen Ableitungen
        (diskrete Differenzen) und die erste Fundamentalform berechnet.

    @param surface_name: Name der Fläche ('sphere','torus','saddle','hyperbolic_paraboloid')
    @param param_range: (min, max) Parameterbereich für u und v
    @param resolution: Anzahl Gitterpunkte pro Parameterachse (Standard: 50)
    @return: matplotlib Figure-Objekt mit 3D-Surface-Plot
    @raises ValueError: wenn surface_name unbekannt ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Parameterwerte gleichmäßig verteilen
    u = np.linspace(param_range[0], param_range[1], resolution)
    v = np.linspace(param_range[0], param_range[1], resolution)
    U, V = np.meshgrid(u, v)

    # --- Fläche in kartesischen Koordinaten berechnen ---
    if surface_name == 'sphere':
        # Einheitssphäre: (u,v) = (θ, φ) ∈ [0,π] × [0,2π]
        # Für allgemeinen param_range normieren wir auf [0,π] × [0,2π]
        t1 = np.pi * (U - param_range[0]) / (param_range[1] - param_range[0])
        t2 = 2 * np.pi * (V - param_range[0]) / (param_range[1] - param_range[0])
        X = np.sin(t1) * np.cos(t2)
        Y = np.sin(t1) * np.sin(t2)
        Z = np.cos(t1)
        # Gaußsche Krümmung der Einheitssphäre: K = 1 überall
        K = np.ones_like(X)

    elif surface_name == 'torus':
        # Torus mit großem Radius R=2 und kleinem Radius r=1
        # x(u,v) = (R + r·cos(u))·cos(v), u,v ∈ [0, 2π]
        R_torus, r_torus = 2.0, 1.0
        # Normieren auf [0, 2π]
        t1 = 2 * np.pi * (U - param_range[0]) / (param_range[1] - param_range[0])
        t2 = 2 * np.pi * (V - param_range[0]) / (param_range[1] - param_range[0])
        X = (R_torus + r_torus * np.cos(t1)) * np.cos(t2)
        Y = (R_torus + r_torus * np.cos(t1)) * np.sin(t2)
        Z = r_torus * np.sin(t1)
        # Gaußsche Krümmung des Torus:
        #   K(u,v) = cos(u) / (r · (R + r·cos(u)))
        K = np.cos(t1) / (r_torus * (R_torus + r_torus * np.cos(t1)))

    elif surface_name == 'saddle':
        # Affines Sattelflächenparaboloid: z = x² - y²
        X = U
        Y = V
        Z = U**2 - V**2
        # Gaußsche Krümmung des Sattelparaboloids z = f(x,y) mit f=x²-y²:
        #   K = (f_xx·f_yy - f_xy²) / (1 + f_x² + f_y²)²
        #   f_x = 2x, f_y = -2y, f_xx = 2, f_yy = -2, f_xy = 0
        f_x = 2 * U
        f_y = -2 * V
        f_xx = np.full_like(U, 2.0)
        f_yy = np.full_like(V, -2.0)
        f_xy = np.zeros_like(U)
        denom = (1 + f_x**2 + f_y**2)**2
        K = (f_xx * f_yy - f_xy**2) / denom

    elif surface_name == 'hyperbolic_paraboloid':
        # Hyperbolisches Paraboloid: z = x·y
        X = U
        Y = V
        Z = U * V
        # Gaußsche Krümmung: f_x = y, f_y = x, f_xx = 0, f_yy = 0, f_xy = 1
        #   K = (0·0 - 1²) / (1 + y² + x²)² = -1 / (1 + x² + y²)²
        f_x = V   # ∂(x·y)/∂x = y
        f_y = U   # ∂(x·y)/∂y = x
        denom = (1 + f_x**2 + f_y**2)**2
        K = -1.0 / denom

    else:
        raise ValueError(
            f"Unbekannte Fläche: '{surface_name}'. "
            "Unterstützt: 'sphere', 'torus', 'saddle', 'hyperbolic_paraboloid'"
        )

    # --- Plot erstellen ---
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Krümmungswerte für Farbkodierung normieren (NaN-sicher)
    K_clean = np.where(np.isfinite(K), K, np.nan)

    # Surface-Plot mit Gaußscher Krümmung als Farbkodierung
    surf = ax.plot_surface(
        X, Y, Z,
        facecolors=plt.cm.RdYlBu(
            (K_clean - np.nanmin(K_clean)) /
            max(np.nanmax(K_clean) - np.nanmin(K_clean), 1e-12)
        ),
        alpha=0.85,
        linewidth=0,
        antialiased=True
    )
    # Farbskala als separates Mappable
    sm = plt.cm.ScalarMappable(
        cmap='RdYlBu',
        norm=plt.Normalize(vmin=float(np.nanmin(K_clean)),
                           vmax=float(np.nanmax(K_clean)))
    )
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.5, label='Gaußsche Krümmung K')

    ax.set_title(f'Gaußsche Krümmung: {surface_name}', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    return fig


# ===========================================================================
# GEODÄTEN-VISUALISIERUNG
# ===========================================================================

def plot_geodesic_on_sphere(
    start_angle: Tuple[float, float],
    direction: Tuple[float, float],
    n_steps: int = 200
) -> plt.Figure:
    """
    @brief Visualisiert eine Geodäte (Großkreis) auf der Einheitssphäre.
    @description
        Eine Geodäte auf der Einheitssphäre ist ein Großkreis.
        Die Geodätengleichung auf S² lautet:
            θ'' = sin(θ)·cos(θ)·(φ')²
            φ'' = -2·cos(θ)/sin(θ) · θ'·φ'

        Diese wird mit Runge-Kutta 4. Ordnung numerisch integriert.
        Der Pfad wird als rote Kurve auf der halbtransparenten Sphäre gezeichnet.

        Koordinaten: sphärische Polarkoordinaten (θ, φ) mit
            x = sin(θ)·cos(φ)
            y = sin(θ)·sin(φ)
            z = cos(θ)

    @param start_angle: (θ₀, φ₀) Startwinkel in Bogenmaß (Polar-, Azimutwinkel)
    @param direction: (dθ/dt, dφ/dt) Anfangsgeschwindigkeit im Parameterraum
    @param n_steps: Anzahl Integrationsschritte (Standard: 200)
    @return: matplotlib Figure-Objekt mit 3D-Sphäre und Geodäten-Pfad
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # --- Geodätengleichung auf S² als System 1. Ordnung ---
    def geodaeten_rhs(state: np.ndarray) -> np.ndarray:
        """
        @brief Rechte Seite der Geodäten-DGL auf S².
        @param state: [θ, φ, θ', φ'] – aktueller Zustand
        @return: Zeitableitung [θ', φ', θ'', φ'']
        """
        theta, phi, dtheta, dphi = state
        # Regularisierung: sin(θ) != 0 vermeiden (Pol-Singularität)
        sin_theta = np.sin(theta)
        if abs(sin_theta) < 1e-8:
            sin_theta = 1e-8
        # θ'' = sin(θ)·cos(θ)·(φ')²
        ddtheta = np.sin(theta) * np.cos(theta) * dphi**2
        # φ'' = -2·cos(θ)/sin(θ) · θ'·φ'
        ddphi = -2.0 * (np.cos(theta) / sin_theta) * dtheta * dphi
        return np.array([dtheta, dphi, ddtheta, ddphi])

    # Gesamtlänge des Pfades: ein vollständiger Großkreis bei Einheitsgeschwindigkeit
    t_end = 2.0 * math.pi / max(
        math.sqrt(direction[0]**2 + direction[1]**2), 1e-8
    )
    dt = t_end / n_steps

    # Anfangszustand
    state = np.array([
        float(start_angle[0]),
        float(start_angle[1]),
        float(direction[0]),
        float(direction[1])
    ])

    # --- Runge-Kutta 4 Integration ---
    thetas = [state[0]]
    phis = [state[1]]
    for _ in range(n_steps):
        k1 = geodaeten_rhs(state)
        k2 = geodaeten_rhs(state + 0.5 * dt * k1)
        k3 = geodaeten_rhs(state + 0.5 * dt * k2)
        k4 = geodaeten_rhs(state + dt * k3)
        state = state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        thetas.append(state[0])
        phis.append(state[1])

    # Sphärische Koordinaten → kartesisch
    thetas = np.array(thetas)
    phis = np.array(phis)
    gx = np.sin(thetas) * np.cos(phis)
    gy = np.sin(thetas) * np.sin(phis)
    gz = np.cos(thetas)

    # --- Sphäre zeichnen ---
    u_sph = np.linspace(0, 2 * np.pi, 40)
    v_sph = np.linspace(0, np.pi, 20)
    Xs = np.outer(np.cos(u_sph), np.sin(v_sph))
    Ys = np.outer(np.sin(u_sph), np.sin(v_sph))
    Zs = np.outer(np.ones(40), np.cos(v_sph))

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    # Halbtransparente Sphäre
    ax.plot_surface(Xs, Ys, Zs, color='lightblue', alpha=0.3, linewidth=0)
    # Geodätenpfad als rote Kurve
    ax.plot(gx, gy, gz, 'r-', linewidth=2.5, label='Geodäte')
    # Startpunkt markieren
    ax.scatter([gx[0]], [gy[0]], [gz[0]], color='green', s=60, zorder=5, label='Start')

    ax.set_title('Geodäte auf der Einheitssphäre', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    return fig


def plot_geodesic_on_torus(
    R: float = 2.0,
    r: float = 1.0,
    start_params: Tuple[float, float] = (0.0, 0.0),
    direction: Tuple[float, float] = (1.0, 0.3),
    n_steps: int = 500
) -> plt.Figure:
    """
    @brief Visualisiert eine Geodäte auf einem Torus.
    @description
        Ein Torus mit großem Radius R und kleinem Radius r wird durch die
        Parametrisierung beschrieben:
            x(u,v) = (R + r·cos(u))·cos(v)
            y(u,v) = (R + r·cos(u))·sin(v)
            z(u,v) = r·sin(u)

        Die Geodätengleichungen auf dem Torus (Christoffel-Symbole):
            u'' = sin(u)·cos(u)·(v')²·(R + r·cos(u))/r
            v'' = -2·sin(u)·u'·v' / (R + r·cos(u)) · r

        Integration per Runge-Kutta 4. Ordnung.

    @param R: Großer Torusradius (Mittelpunkt des Rohres vom Zentrum)
    @param r: Kleiner Torusradius (Rohrquerschnitt)
    @param start_params: (u₀, v₀) Startparameter in Bogenmaß
    @param direction: (du/dt₀, dv/dt₀) Anfangsgeschwindigkeit im Parameterraum
    @param n_steps: Anzahl Integrationsschritte (Standard: 500)
    @return: matplotlib Figure-Objekt mit 3D-Torus und Geodäten-Pfad
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # --- Geodätengleichung auf dem Torus als System 1. Ordnung ---
    def torus_geodaete_rhs(state: np.ndarray) -> np.ndarray:
        """
        @brief Rechte Seite der Geodäten-DGL auf dem Torus.
        @param state: [u, v, u', v'] – aktueller Zustand
        @return: Zeitableitung [u', v', u'', v'']
        """
        u, v, du, dv = state
        cos_u = np.cos(u)
        sin_u = np.sin(u)
        denom = R + r * cos_u
        # Regularisierung: Nenner != 0
        if abs(denom) < 1e-8:
            denom = 1e-8
        # u'' = sin(u)·cos(u)·(v')² · (R + r·cos(u)) / r
        ddu = sin_u * cos_u * dv**2 * denom / r
        # v'' = -2·r·sin(u)·u'·v' / (R + r·cos(u))
        ddv = -2.0 * r * sin_u * du * dv / denom
        return np.array([du, dv, ddu, ddv])

    # Zeitschritt
    t_end = 4.0 * math.pi / max(
        math.sqrt(direction[0]**2 + direction[1]**2), 1e-8
    )
    dt = t_end / n_steps

    # Anfangszustand
    state = np.array([
        float(start_params[0]),
        float(start_params[1]),
        float(direction[0]),
        float(direction[1])
    ])

    # --- Runge-Kutta 4 ---
    us = [state[0]]
    vs = [state[1]]
    for _ in range(n_steps):
        k1 = torus_geodaete_rhs(state)
        k2 = torus_geodaete_rhs(state + 0.5 * dt * k1)
        k3 = torus_geodaete_rhs(state + 0.5 * dt * k2)
        k4 = torus_geodaete_rhs(state + dt * k3)
        state = state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        us.append(state[0])
        vs.append(state[1])

    us = np.array(us)
    vs = np.array(vs)

    # Torus-Parametrisierung → kartesisch
    gx = (R + r * np.cos(us)) * np.cos(vs)
    gy = (R + r * np.cos(us)) * np.sin(vs)
    gz = r * np.sin(us)

    # --- Torus-Oberfläche zeichnen ---
    u_surf = np.linspace(0, 2 * np.pi, 40)
    v_surf = np.linspace(0, 2 * np.pi, 40)
    U_s, V_s = np.meshgrid(u_surf, v_surf)
    Xt = (R + r * np.cos(U_s)) * np.cos(V_s)
    Yt = (R + r * np.cos(U_s)) * np.sin(V_s)
    Zt = r * np.sin(U_s)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(Xt, Yt, Zt, color='lightcyan', alpha=0.3, linewidth=0)
    ax.plot(gx, gy, gz, 'r-', linewidth=2.0, label='Geodäte')
    ax.scatter([gx[0]], [gy[0]], [gz[0]], color='green', s=60, zorder=5, label='Start')

    ax.set_title(f'Geodäte auf dem Torus (R={R}, r={r})', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    return fig


# ===========================================================================
# ADAPTIVES GITTER FÜR 2D-PLOTS
# ===========================================================================

def plot_function_2d_adaptive(
    func_str: str,
    x_range: Tuple[float, float] = (-5.0, 5.0),
    base_points: int = 200,
    refinement_levels: int = 3,
    return_points: bool = False
):
    """
    @brief Plottet eine Funktion mit adaptivem Gitter: In Regionen hoher
           Krümmung werden automatisch mehr Auswertungspunkte eingefügt.
    @description
        Algorithmus:
            1. base_points gleichmäßig verteilt → Basispunkte x₀
            2. Zweite numerische Ableitung |f''(x)| schätzen (diskrete Differenz)
            3. Punkte mit |f''| > threshold werden verfeinert
            4. Rekursiv bis refinement_levels Iterationen

        Dadurch wird der Plot in "glatten" Regionen mit weniger Punkten
        dargestellt und in stark gekrümmten Regionen (Sattelpunkte, Extrema,
        steile Anstiege) mit mehr Punkten, was Aliasing verhindert.

        Unterstützte Ausdrücke in func_str (per eval):
            'sin(x)', 'cos(x)', 'exp(x)', 'log(x)', '1/x', 'x**2', etc.
            (Alle numpy-Funktionen sind per 'np.' erreichbar)
            Kurzformen sin/cos/exp/log werden automatisch auf np.sin etc. gemappt.

    @param func_str: Funktionsausdruck als String, z.B. 'sin(x)' oder 'x**2'
    @param x_range: (x_min, x_max) Plot-Intervall
    @param base_points: Anzahl äquidistanter Startpunkte (Standard: 200)
    @param refinement_levels: Maximale Verfeinerungstiefe (Standard: 3)
    @param return_points: Falls True, werden (x_array, y_array) zurückgegeben
                          statt einer Figure (für Tests)
    @return: matplotlib Figure-Objekt ODER (x_array, y_array) wenn return_points=True
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Funktion sicher mit sympy.lambdify parsen (kein eval auf Benutzereingabe)
    try:
        import sympy as _sp
        from sympy.parsing.sympy_parser import parse_expr, standard_transformations, \
            implicit_multiplication_application
        _x_sym = _sp.Symbol('x')
        _transformations = standard_transformations + (implicit_multiplication_application,)
        _expr = parse_expr(func_str, local_dict={'x': _x_sym},
                           transformations=_transformations)
        _lambdified = _sp.lambdify(_x_sym, _expr, modules=['numpy'])
    except Exception as _parse_err:
        raise ValueError(f"Ungültiger Funktionsausdruck '{func_str}': {_parse_err}")

    def _eval_func(x_val: np.ndarray) -> np.ndarray:
        """Wertet func_str an einem numpy-Array x_val aus (via lambdify)."""
        try:
            result = _lambdified(x_val)
            return np.asarray(result, dtype=float)
        except Exception:
            return np.full_like(x_val, np.nan, dtype=float)

    # --- Schritt 1: Basispunkte erzeugen ---
    x_pts = np.linspace(x_range[0], x_range[1], base_points)
    y_pts = _eval_func(x_pts)

    # --- Schritte 2–4: Iterative Verfeinerung ---
    for _level in range(refinement_levels):
        if len(x_pts) < 3:
            break

        # Zweite Ableitung (diskrete Differenz zweiter Ordnung) schätzen
        # f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
        # Hier approximiert über benachbarte Punkte im Array
        y_clean = np.where(np.isfinite(y_pts), y_pts, 0.0)
        d2y = np.abs(np.diff(y_clean, n=2))  # Länge: n-2

        # Threshold: 95. Perzentil der Krümmungswerte → obere 5% verfeinern
        if d2y.max() < 1e-12:
            break  # Funktion ist perfekt linear, keine Verfeinerung nötig
        threshold = np.percentile(d2y, 95)

        # Indizes der hochgekrümmten Segmente (Index bezieht sich auf Mittelpunkte)
        hoch_indices = np.where(d2y > threshold)[0] + 1  # +1: Versatz auf x_pts

        # Neue Punkte zwischen hochgekrümmten Nachbarn einfügen
        neue_x = []
        for idx in hoch_indices:
            if idx < len(x_pts) - 1:
                # Mittelpunkt zwischen x[idx] und x[idx+1]
                neue_x.append(0.5 * (x_pts[idx] + x_pts[idx + 1]))
            if idx > 0:
                # Mittelpunkt zwischen x[idx-1] und x[idx]
                neue_x.append(0.5 * (x_pts[idx - 1] + x_pts[idx]))

        if not neue_x:
            break

        # Neue Punkte zum Array hinzufügen und sortieren
        alle_x = np.unique(np.concatenate([x_pts, np.array(neue_x)]))
        alle_y = _eval_func(alle_x)
        x_pts = alle_x
        y_pts = alle_y

    # Nur gültige Werte behalten
    mask = np.isfinite(y_pts)
    x_final = x_pts[mask]
    y_final = y_pts[mask]

    # Falls nur Rohdaten gewünscht (für Tests)
    if return_points:
        return x_final, y_final

    # --- Plot erzeugen ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_final, y_final, 'b-', linewidth=1.5)
    ax.set_title(f'f(x) = {func_str} (adaptiv, {len(x_final)} Punkte)', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)

    return fig


# ===========================================================================
# INTERAKTIVE PLOTS (matplotlib.widgets)
# ===========================================================================

def interactive_plot(
    func,
    x_range: tuple,
    param_name: str = 'a',
    param_range: tuple = (0, 10),
    save_path: str = None
):
    """
    @brief Erstellt einen interaktiven Plot mit Slider-Widget.
    @description
        Erstellt einen Matplotlib-Plot mit einem Slider für einen Parameter.
        Der Slider ermöglicht es, den Parameter interaktiv zu verändern
        und die Auswirkungen auf den Graphen live zu sehen.

        Nur verfügbar wenn matplotlib.widgets importierbar ist.

    @param func: Funktion f(x, param) mit Slider-Parameter.
    @param x_range: (x_min, x_max) – Wertebereich der x-Achse.
    @param param_name: Name des Parameters (erscheint auf Slider).
    @param param_range: (param_min, param_max) – Wertebereich des Parameters.
    @param save_path: Speicherpfad für PNG (None = interaktiv anzeigen).
    @return: matplotlib.figure.Figure-Objekt.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    try:
        from matplotlib.widgets import Slider
    except ImportError:
        raise ImportError("matplotlib.widgets ist nicht verfügbar.")

    # Gitter aufbauen
    x_arr = np.linspace(x_range[0], x_range[1], 500)
    param_init = (param_range[0] + param_range[1]) / 2.0

    # Figure und Plot-Achse erstellen (Platz für Slider unten einplanen)
    fig, ax = plt.subplots(figsize=(9, 6))
    plt.subplots_adjust(bottom=0.2)

    # Initialen Plot zeichnen
    y_init = np.array([func(x, param_init) for x in x_arr])
    line, = ax.plot(x_arr, y_init, 'b-', linewidth=2)
    ax.set_xlabel('x')
    ax.set_ylabel(f'f(x, {param_name})')
    ax.set_title(f'Interaktiver Plot: f(x, {param_name})')
    ax.grid(True, alpha=0.3)

    # Slider erstellen
    ax_slider = plt.axes([0.1, 0.05, 0.8, 0.04])
    slider = Slider(
        ax_slider,
        param_name,
        param_range[0],
        param_range[1],
        valinit=param_init
    )

    # Update-Funktion für Slider-Ereignis
    def update(val):
        """Aktualisiert den Plot wenn Slider bewegt wird."""
        param = slider.val
        y_new = np.array([func(x, param) for x in x_arr])
        line.set_ydata(y_new)
        fig.canvas.draw_idle()

    slider.on_changed(update)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def interactive_3d(
    func_z,
    x_range: tuple,
    y_range: tuple,
    n: int = 50,
    save_path: str = None
):
    """
    @brief Erstellt einen interaktiven 3D-Plot mit Elevation/Azimuth-Slidern.
    @description
        Erzeugt eine 3D-Oberfläche mit matplotlib und zwei Slider-Widgets
        zum interaktiven Rotieren der Ansicht (Elevation und Azimuth).

    @param func_z: Funktion f(x, y) → z-Wert.
    @param x_range: (x_min, x_max) – Wertebereich der x-Achse.
    @param y_range: (y_min, y_max) – Wertebereich der y-Achse.
    @param n: Gitterpunkte pro Achse (n×n-Gitter).
    @param save_path: Speicherpfad (None = interaktiv anzeigen).
    @return: matplotlib.figure.Figure-Objekt.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    try:
        from matplotlib.widgets import Slider
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        raise ImportError("matplotlib.widgets oder mpl_toolkits nicht verfügbar.")

    # Gitter aufbauen
    x_arr = np.linspace(x_range[0], x_range[1], n)
    y_arr = np.linspace(y_range[0], y_range[1], n)
    X, Y = np.meshgrid(x_arr, y_arr)

    # Z-Werte berechnen
    Z = np.vectorize(func_z)(X, Y)

    # Figure erstellen
    fig = plt.figure(figsize=(10, 8))
    plt.subplots_adjust(bottom=0.2)
    ax3d = fig.add_subplot(111, projection='3d')
    surf = [ax3d.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)]

    ax3d.set_xlabel('x')
    ax3d.set_ylabel('y')
    ax3d.set_zlabel('z')
    ax3d.set_title('3D-Plot (interaktiv)')

    # Slider für Elevation (Neigungswinkel)
    ax_elev = plt.axes([0.1, 0.08, 0.35, 0.03])
    slider_elev = Slider(ax_elev, 'Elevation', 0, 90, valinit=30)

    # Slider für Azimuth (Drehwinkel)
    ax_azim = plt.axes([0.55, 0.08, 0.35, 0.03])
    slider_azim = Slider(ax_azim, 'Azimuth', -180, 180, valinit=-60)

    def update_view(val):
        """Aktualisiert Blickwinkel beim Slider-Bewegen."""
        ax3d.view_init(elev=slider_elev.val, azim=slider_azim.val)
        fig.canvas.draw_idle()

    slider_elev.on_changed(update_view)
    slider_azim.on_changed(update_view)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


# ===========================================================================
# ADAPTIVES GITTER
# ===========================================================================

def adaptive_plot(
    func,
    x_range: tuple,
    n_initial: int = 100,
    max_refinements: int = 5,
    tol: float = 1e-3,
    save_path: str = None,
    return_points: bool = False
):
    """
    @brief Adaptiver Plot mit automatischer Gitterverfeinerung.
    @description
        Verfeinert das Gitter automatisch in Bereichen mit starker Krümmung
        (nahe Nullstellen, Singularitäten oder schnell wechselnden Funktionswerten).

        Algorithmus:
        1. Anfangsgitter mit n_initial Punkten
        2. Gradient berechnen mit numpy.gradient()
        3. Punkte mit hohem Gradienten verfeinern (Mittelpunkte einfügen)
        4. Bis zu max_refinements Verfeinerungsrunden

    @param func: Funktion f(x) → float.
    @param x_range: (x_min, x_max) – Wertebereich.
    @param n_initial: Anfangsgitter-Größe.
    @param max_refinements: Maximale Verfeinerungsrunden.
    @param tol: Toleranz für Gradienten-Schwellenwert.
    @param save_path: Speicherpfad (None = anzeigen).
    @param return_points: Wenn True, gibt (x, y) zurück statt Figure.
    @return: Figure oder (x_arr, y_arr) wenn return_points=True.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Anfangsgitter erstellen
    x_pts = np.linspace(x_range[0], x_range[1], n_initial)
    y_pts = np.array([func(xi) for xi in x_pts])

    # Adaptive Verfeinerung
    for _ in range(max_refinements):
        # Gradient berechnen (Änderungsrate)
        grad = np.abs(np.gradient(y_pts, x_pts))

        # Schwellenwert: tol-Quantil der Gradienten
        threshold = np.quantile(grad, 1.0 - tol) if len(grad) > 1 else 0

        # Indizes mit hohem Gradienten finden
        high_grad_idx = np.where(grad > threshold)[0]

        if len(high_grad_idx) == 0:
            break  # Keine Verfeinerung mehr nötig

        # Neue Punkte zwischen hochgradigen Nachbarn einfügen
        neue_x = []
        for idx in high_grad_idx:
            if idx < len(x_pts) - 1:
                neue_x.append(0.5 * (x_pts[idx] + x_pts[idx + 1]))
            if idx > 0:
                neue_x.append(0.5 * (x_pts[idx - 1] + x_pts[idx]))

        if not neue_x:
            break

        # Neue Punkte einfügen und sortieren
        alle_x = np.unique(np.concatenate([x_pts, np.array(neue_x)]))
        alle_y = np.array([func(xi) for xi in alle_x])
        x_pts = alle_x
        y_pts = alle_y

    # Nur gültige (finite) Werte behalten
    mask = np.isfinite(y_pts)
    x_final = x_pts[mask]
    y_final = y_pts[mask]

    if return_points:
        return x_final, y_final

    # Plot erstellen
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_final, y_final, 'b-', linewidth=1.5)
    ax.scatter(x_final[::max(1, len(x_final)//20)], y_final[::max(1, len(y_final)//20)],
               color='red', s=10, alpha=0.5, label='Gitterpunkte')
    ax.set_title(f'Adaptiver Plot ({len(x_final)} Punkte nach {max_refinements} Verfeinerungen)')
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.grid(True, alpha=0.3)
    ax.legend()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


# ===========================================================================
# KRÜMMUNGS-VISUALISIERUNG
# ===========================================================================

def plot_gaussian_curvature(
    surface_func,
    u_range: tuple,
    v_range: tuple,
    n: int = 50,
    save_path: str = None,
    return_data: bool = False
):
    """
    @brief 3D-Plot mit Farbkodierung der Gaußschen Krümmung.
    @description
        Visualisiert eine parametrische Fläche und kodiert die Gaußsche Krümmung K
        als Farbe. Die Gaußsche Krümmung ist das Produkt der Hauptkrümmungen:
            K = κ₁ · κ₂

        Farbkodierung:
        - Blau: K < 0 (Sattelform, hyperbolischer Punkt)
        - Grün: K ≈ 0 (zylindrisch, parabolischer Punkt)
        - Rot: K > 0 (kugelförmig, elliptischer Punkt)

    @param surface_func: Funktion (u, v) → (x, y, z).
    @param u_range: (u_min, u_max) – Parameterbereich für u.
    @param v_range: (v_min, v_max) – Parameterbereich für v.
    @param n: Gitterpunkte pro Parameter.
    @param save_path: Speicherpfad (None = anzeigen).
    @param return_data: Wenn True, gibt (fig, K_values) zurück.
    @return: Figure oder (Figure, K_array).
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    from mpl_toolkits.mplot3d import Axes3D

    # Parametergitter aufbauen
    u_arr = np.linspace(u_range[0], u_range[1], n)
    v_arr = np.linspace(v_range[0], v_range[1], n)
    U, V = np.meshgrid(u_arr, v_arr)

    # Koordinaten berechnen
    # surface_func gibt ein Tupel (x,y,z) zurück, daher manuell entpacken
    X = np.zeros((n, n))
    Y = np.zeros((n, n))
    Z = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            xyz = surface_func(U[i, j], V[i, j])
            X[i, j] = xyz[0]
            Y[i, j] = xyz[1]
            Z[i, j] = xyz[2]

    # Gaußsche Krümmung numerisch approximieren
    # Über die Fundamentalformen der Differentialgeometrie
    h = 1e-5  # Schrittweite für numerische Ableitungen

    # Numerische partiell Ableitungen
    def surf_u(u, v):
        """Partielle Ableitung nach u."""
        p1 = np.array(surface_func(u + h, v))
        p2 = np.array(surface_func(u - h, v))
        return (p1 - p2) / (2 * h)

    def surf_v(u, v):
        """Partielle Ableitung nach v."""
        p1 = np.array(surface_func(u, v + h))
        p2 = np.array(surface_func(u, v - h))
        return (p1 - p2) / (2 * h)

    def surf_uu(u, v):
        """Zweite Ableitung nach u."""
        p1 = np.array(surface_func(u + h, v))
        p0 = np.array(surface_func(u, v))
        p2 = np.array(surface_func(u - h, v))
        return (p1 - 2*p0 + p2) / (h**2)

    def surf_vv(u, v):
        """Zweite Ableitung nach v."""
        p1 = np.array(surface_func(u, v + h))
        p0 = np.array(surface_func(u, v))
        p2 = np.array(surface_func(u, v - h))
        return (p1 - 2*p0 + p2) / (h**2)

    def surf_uv(u, v):
        """Gemischte Ableitung."""
        p1 = np.array(surface_func(u + h, v + h))
        p2 = np.array(surface_func(u + h, v - h))
        p3 = np.array(surface_func(u - h, v + h))
        p4 = np.array(surface_func(u - h, v - h))
        return (p1 - p2 - p3 + p4) / (4 * h**2)

    # Krümmung für jeden Gitterpunkt berechnen
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            u, v = U[i, j], V[i, j]

            # Erste Fundamentalform
            r_u = surf_u(u, v)
            r_v = surf_v(u, v)
            E = np.dot(r_u, r_u)
            F = np.dot(r_u, r_v)
            G = np.dot(r_v, r_v)

            # Normalenvektor
            N = np.cross(r_u, r_v)
            N_norm = np.linalg.norm(N)
            if N_norm < 1e-12:
                continue
            N = N / N_norm

            # Zweite Fundamentalform
            r_uu = surf_uu(u, v)
            r_vv = surf_vv(u, v)
            r_uv = surf_uv(u, v)
            e = np.dot(r_uu, N)
            f = np.dot(r_uv, N)
            g_coeff = np.dot(r_vv, N)

            # Gaußsche Krümmung K = (eg - f²) / (EG - F²)
            denom = E * G - F**2
            if abs(denom) > 1e-12:
                K[i, j] = (e * g_coeff - f**2) / denom

    # Plot erstellen
    fig = plt.figure(figsize=(12, 5))

    # Linke Seite: Oberfläche
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot_surface(X, Y, Z, alpha=0.7, cmap='viridis')
    ax1.set_title('Oberfläche')

    # Rechte Seite: Krümmungskarte
    ax2 = fig.add_subplot(122)
    vmax = np.percentile(np.abs(K), 95)
    im = ax2.contourf(U, V, K, levels=50, cmap='RdBu_r',
                       vmin=-vmax, vmax=vmax)
    plt.colorbar(im, ax=ax2, label='Gaußsche Krümmung K')
    ax2.set_xlabel('u')
    ax2.set_ylabel('v')
    ax2.set_title('Gaußsche Krümmung K(u, v)')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    if return_data:
        return fig, K
    return fig


def plot_sphere_curvature(R: float = 1.0, n: int = 30, save_path: str = None):
    """
    @brief Visualisiert die Gaußsche Krümmung einer Kugel.
    @description
        Auf einer Kugel mit Radius R ist die Gaußsche Krümmung überall K = 1/R².
        Diese Funktion visualisiert die Kugel und bestätigt dies numerisch.

    @param R: Kugelradius.
    @param n: Gitterpunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    def sphere(u, v):
        """Kugelparametrisierung: (θ, φ) → (x, y, z)."""
        x = R * np.sin(u) * np.cos(v)
        y = R * np.sin(u) * np.sin(v)
        z = R * np.cos(u)
        return (x, y, z)

    return plot_gaussian_curvature(
        sphere,
        u_range=(0.1, np.pi - 0.1),
        v_range=(0, 2 * np.pi),
        n=n,
        save_path=save_path
    )


def plot_torus_curvature(R: float = 2.0, r: float = 1.0, n: int = 40, save_path: str = None):
    """
    @brief Visualisiert die Gaußsche Krümmung eines Torus.
    @description
        Auf einem Torus wechselt die Krümmung das Vorzeichen:
        - Außenseite (φ-Seite): K > 0 (wie Kugel)
        - Innenseite (θ-Seite): K < 0 (wie Sattel)

    @param R: Großer Radius (Mittelpunkt → Rohrmitten-Kreis).
    @param r: Kleiner Radius (Rohrquerschnitt).
    @param n: Gitterpunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    def torus(u, v):
        """Torusparametrisierung."""
        x = (R + r * np.cos(v)) * np.cos(u)
        y = (R + r * np.cos(v)) * np.sin(u)
        z = r * np.sin(v)
        return (x, y, z)

    return plot_gaussian_curvature(
        torus,
        u_range=(0, 2 * np.pi),
        v_range=(0, 2 * np.pi),
        n=n,
        save_path=save_path
    )


def plot_saddle_curvature(n: int = 40, save_path: str = None):
    """
    @brief Visualisiert die Gaußsche Krümmung einer Sattelfläche z = x² - y².
    @description
        Die Sattelfläche z = x² - y² hat negative Gaußsche Krümmung K < 0
        überall (außer am Sattelpunkt), da die Hauptkrümmungen entgegengesetzte
        Vorzeichen haben.

    @param n: Gitterpunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    def saddle(u, v):
        """Sattel-Parametrisierung: z = u² - v²."""
        return (u, v, u**2 - v**2)

    return plot_gaussian_curvature(
        saddle,
        u_range=(-2, 2),
        v_range=(-2, 2),
        n=n,
        save_path=save_path
    )


# ===========================================================================
# PDE-ANIMATION
# ===========================================================================

def animate_heat_equation(
    u0,
    alpha: float,
    x_range: tuple,
    t_range: tuple,
    nx: int = 100,
    nt: int = 200,
    save_path: str = None
):
    """
    @brief Animiert die Wärmeleitungsgleichung ∂u/∂t = α ∂²u/∂x².
    @description
        Löst die 1D-Wärmeleitungsgleichung numerisch via Crank-Nicolson-Methode
        (unbedingt stabil) und erzeugt eine GIF-Animation.

        Wärmeleitungsgleichung (parabolische PDE):
            ∂u/∂t = α · ∂²u/∂x²

        Crank-Nicolson-Methode (θ = 1/2):
            (u^{n+1} - u^n)/Δt = (α/2)(u_xx^n + u_xx^{n+1})

        Randbedingungen: Dirichlet (u = 0 an den Rändern).

    @param u0: Anfangsbedingung als Funktion u0(x) oder numpy-Array.
    @param alpha: Wärmediffusivitätskoeffizient α > 0.
    @param x_range: (x_min, x_max) – Raumintervall.
    @param t_range: (t_min, t_max) – Zeitintervall.
    @param nx: Anzahl der Raumpunkte.
    @param nt: Anzahl der Zeitschritte (und Animationsframes).
    @param save_path: Pfad zum Speichern der GIF-Animation.
    @return: matplotlib.animation.FuncAnimation-Objekt.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import matplotlib.animation as animation
    from scipy.linalg import solve_banded

    # Gitter aufbauen
    x = np.linspace(x_range[0], x_range[1], nx)
    dx = x[1] - x[0]
    dt = (t_range[1] - t_range[0]) / nt
    r = alpha * dt / (2 * dx**2)  # Crank-Nicolson-Parameter

    # Anfangsbedingung aufbauen
    if callable(u0):
        u = np.array([u0(xi) for xi in x])
    else:
        u = np.array(u0, dtype=float)
        if len(u) != nx:
            u = np.interp(x, np.linspace(x_range[0], x_range[1], len(u)), u)

    # Dirichlet-Randbedingungen: u[0] = u[-1] = 0
    u[0] = 0.0
    u[-1] = 0.0

    # Tridiagonale Matrizen für Crank-Nicolson
    # Linke Seite (implizit): A = I + r*T
    # Rechte Seite (explizit): B = I - r*T
    n_inner = nx - 2  # Innere Punkte

    # Bandmatrix-Format für scipy.linalg.solve_banded
    # ab[0] = obere Diagonale, ab[1] = Hauptdiagonale, ab[2] = untere Diagonale
    ab = np.zeros((3, n_inner))
    ab[0, 1:] = -r         # Obere Diagonale
    ab[1, :] = 1 + 2 * r   # Hauptdiagonale
    ab[2, :-1] = -r         # Untere Diagonale

    # Animation-Daten speichern (alle nt Frames)
    u_history = [u.copy()]

    # Zeitschritte berechnen
    u_inner = u[1:-1].copy()
    for step in range(nt):
        # Rechte Seite berechnen (expliziter Teil)
        rhs = np.zeros(n_inner)
        rhs[1:-1] = r * u_inner[:-2] + (1 - 2*r) * u_inner[1:-1] + r * u_inner[2:]
        rhs[0] = r * 0 + (1 - 2*r) * u_inner[0] + r * u_inner[1]  # linker Rand = 0
        rhs[-1] = r * u_inner[-2] + (1 - 2*r) * u_inner[-1] + r * 0  # rechter Rand = 0

        # Löse tridiagonales System
        u_inner = solve_banded((1, 1), ab, rhs)

        # Vollständiges u speichern
        u_full = np.concatenate([[0.0], u_inner, [0.0]])
        u_history.append(u_full.copy())

    # Animation erstellen
    fig, ax = plt.subplots(figsize=(9, 5))
    line, = ax.plot(x, u_history[0], 'b-', linewidth=2)
    ax.set_xlim(x_range)
    u_max = max(np.max(np.abs(u)) for u in u_history)
    ax.set_ylim(-u_max * 1.1, u_max * 1.1)
    ax.set_xlabel('x')
    ax.set_ylabel('u(x, t)')
    ax.grid(True, alpha=0.3)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    dt_total = t_range[1] - t_range[0]

    def update_frame(frame_idx):
        """Aktualisiert den Plot für jeden Animationsframe."""
        line.set_ydata(u_history[frame_idx])
        t_val = t_range[0] + frame_idx * dt_total / nt
        time_text.set_text(f't = {t_val:.4f}')
        ax.set_title(f'Wärmeleitungsgleichung: ∂u/∂t = {alpha}·∂²u/∂x²\n'
                     f't = {t_val:.4f}')
        return line, time_text

    # Nur jeden 5. Frame animieren für flüssigere Darstellung
    frames = list(range(0, nt + 1, max(1, nt // 50)))
    anim = animation.FuncAnimation(
        fig, update_frame, frames=frames,
        interval=50, blit=True
    )

    if save_path:
        try:
            anim.save(save_path, writer='pillow', fps=20)
        except Exception:
            fig.savefig(save_path.replace('.gif', '.png'), dpi=100)

    return anim


# ===========================================================================
# MAßTHEORIE-VISUALISIERUNG
# ===========================================================================

def plot_cantor_set(n_steps: int = 6, save_path: str = None, return_fig: bool = True):
    """
    @brief Visualisiert die Cantor-Mengen-Konstruktion.
    @description
        Die Cantor-Menge wird konstruiert, indem iterativ das mittlere Drittel
        jedes verbleibenden Intervals entfernt wird.

        Schritt 0: [0, 1]
        Schritt 1: [0, 1/3] ∪ [2/3, 1]
        Schritt 2: [0, 1/9] ∪ [2/9, 1/3] ∪ [2/3, 7/9] ∪ [8/9, 1]
        ...

        Die Cantor-Menge hat:
        - Maß (Länge) = 0
        - Überabzählbar viele Punkte (Mächtigkeit des Kontinuums)
        - Fraktale Dimension ≈ 0.6309 (log(2)/log(3))

    @param n_steps: Anzahl der Konstruktionsschritte.
    @param save_path: Speicherpfad.
    @param return_fig: Wenn True, Figure zurückgeben.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Cantor-Menge iterativ konstruieren
    # Starte mit Einheitsintervall [0, 1]
    intervals = [(0.0, 1.0)]

    # Jede Iteration: Mittleres Drittel jedes Intervals entfernen
    all_steps = [list(intervals)]
    for step in range(n_steps):
        new_intervals = []
        for a, b in intervals:
            third = (b - a) / 3.0
            # Linkes Drittel behalten: [a, a+1/3(b-a)]
            new_intervals.append((a, a + third))
            # Rechtes Drittel behalten: [b-1/3(b-a), b]
            new_intervals.append((b - third, b))
        intervals = new_intervals
        all_steps.append(list(intervals))

    # Plot erstellen
    fig, ax = plt.subplots(figsize=(10, max(4, n_steps + 1)))

    # Jeden Konstruktionsschritt als horizontale Balken zeigen
    for step_idx, step_intervals in enumerate(all_steps):
        y = n_steps - step_idx  # Obere Schritte zuerst
        for a, b in step_intervals:
            ax.barh(y, b - a, left=a, height=0.6, color='black', alpha=0.8)

    # Achsen konfigurieren
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_steps + 1)
    ax.set_xlabel('x')
    ax.set_ylabel('Konstruktionsschritt')
    ax.set_yticks(range(1, n_steps + 2))
    ax.set_yticklabels([f'Schritt {n_steps - i}' for i in range(n_steps + 1)])
    ax.set_title(f'Cantor-Menge ({n_steps} Schritte)\n'
                 f'Fraktale Dimension ≈ {np.log(2)/np.log(3):.4f}')
    ax.grid(True, alpha=0.2, axis='x')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_cantor_function(n: int = 100, save_path: str = None):
    """
    @brief Visualisiert die Cantor-Funktion (Teufels-Treppe).
    @description
        Die Cantor-Funktion f: [0,1] → [0,1] ist:
        - Monoton nicht-abnehmend
        - Stetig
        - Fast überall konstant (Ableitung = 0 f.ü.)
        - Aber NICHT absolut stetig

        Konstruktion via ternäre (Basis-3) Darstellung:
        - Wenn x in entferntem Drittel: f(x) = Mittelwert der Grenzwerte
        - Sonst: rekursiv aus ternärer Darstellung

    @param n: Anzahl der Auswertungspunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    def cantor_val(x: float, depth: int = 20) -> float:
        """
        @brief Wertet die Cantor-Funktion an x ∈ [0,1] aus.
        @description
            Algorithmus über ternäre Darstellung:
            1. Schreibe x in Basis 3
            2. Wenn eine 1 vorkommt: ersetze alle Stellen nach der ersten 1 durch 0
            3. Ersetze alle 2en durch 1en → das ist die Basis-2-Darstellung von f(x)
        """
        a, b = 0.0, 1.0
        result = 0.0
        step = 0.5

        for _ in range(depth):
            mid = (a + b) / 2.0
            third = (b - a) / 3.0

            if x < a + third:
                # Linkes Drittel: weiter mit [a, a+third]
                b = a + third
            elif x > b - third:
                # Rechtes Drittel: weiter mit [b-third, b]
                result += step
                a = b - third
            else:
                # Mittleres Drittel (entfernt): f = result + step
                return result + step

            step /= 2.0

        return result + step / 2.0

    # Auswertungspunkte
    x_vals = np.linspace(0, 1, n)
    y_vals = np.array([cantor_val(xi) for xi in x_vals])

    # Plot erstellen
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_vals, y_vals, 'b-', linewidth=2)
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.set_title("Cantor-Funktion (Teufels-Treppe)\n"
                 "f'(x) = 0 fast überall, aber nicht absolut stetig")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.4, label='Identität (zum Vergleich)')
    ax.legend()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


# ===========================================================================
# SPEZIELLE-FUNKTIONEN-GALERIE
# ===========================================================================

def plot_bessel_gallery(
    orders=None,
    x_range: tuple = (0, 20),
    n: int = 500,
    save_path: str = None
):
    """
    @brief Galerie der Bessel-Funktionen verschiedener Ordnungen.
    @description
        Stellt die Bessel-Funktionen J_n(x) erster Art für verschiedene
        Ordnungen n in einem Plot dar.

        Bessel-Differentialgleichung:
            x² y'' + x y' + (x² - n²) y = 0

        J_n(x) entsteht bei Problemen mit Zylindersymmetrie (z.B. Wellen
        in zylindrischen Gefäßen, Wärmeleitung in kreisförmigen Platten).

    @param orders: Liste der Ordnungen (Standard: [0, 1, 2, 3, 4, 5]).
    @param x_range: (x_min, x_max) – Wertebereich.
    @param n: Anzahl der Auswertungspunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if orders is None:
        orders = [0, 1, 2, 3, 4, 5]

    from scipy.special import jn  # Bessel-Funktionen

    x = np.linspace(max(0.001, x_range[0]), x_range[1], n)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Farben für verschiedene Ordnungen
    colors = plt.cm.tab10(np.linspace(0, 1, len(orders)))

    for order, color in zip(orders, colors):
        y = jn(order, x)
        ax.plot(x, y, color=color, linewidth=2, label=f'J_{order}(x)')

    # Nulllinie
    ax.axhline(0, color='black', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('x')
    ax.set_ylabel('J_n(x)')
    ax.set_title('Bessel-Funktionen erster Art J_n(x)')
    ax.set_ylim(-0.5, 1.1)
    ax.legend(loc='upper right', ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_legendre_gallery(
    orders=None,
    x_range: tuple = (-1, 1),
    n: int = 300,
    save_path: str = None
):
    """
    @brief Galerie der Legendre-Polynome verschiedener Ordnungen.
    @description
        Stellt die Legendre-Polynome P_n(x) für verschiedene Ordnungen n
        in einem Plot dar.

        Legendre-Differentialgleichung:
            (1-x²) y'' - 2x y' + n(n+1) y = 0

        P_n(x) entsteht bei Problemen mit Kugelkoordinaten (z.B. Elektrodynamik,
        Quantenmechanik des Wasserstoffatoms).

        Eigenschaften:
        - P_n(1) = 1, P_n(-1) = (-1)^n für alle n
        - Orthogonal: ∫_{-1}^{1} P_m P_n dx = 2/(2n+1) · δ_{mn}
        - Rekursion: (n+1)P_{n+1} = (2n+1)x P_n - n P_{n-1}

    @param orders: Liste der Ordnungen (Standard: [0, 1, 2, 3, 4, 5]).
    @param x_range: (x_min, x_max) – Wertebereich.
    @param n: Anzahl der Auswertungspunkte.
    @param save_path: Speicherpfad.
    @return: Figure.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if orders is None:
        orders = [0, 1, 2, 3, 4, 5]

    from scipy.special import legendre as legendre_poly

    x = np.linspace(x_range[0], x_range[1], n)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Farben für verschiedene Ordnungen
    colors = plt.cm.tab10(np.linspace(0, 1, len(orders)))

    for order, color in zip(orders, colors):
        # scipy.special.legendre gibt ein Polynom-Objekt zurück
        P = legendre_poly(order)
        y = P(x)
        ax.plot(x, y, color=color, linewidth=2, label=f'P_{order}(x)')

    # Nulllinie und Einheitslinie
    ax.axhline(0, color='black', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('x')
    ax.set_ylabel('P_n(x)')
    ax.set_title('Legendre-Polynome P_n(x)')
    ax.legend(loc='upper left', ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


# ===========================================================================
# GAUSS'SCHE KRÜMMUNG (Build 46)
# ===========================================================================

def sphere_surface(u: float, v: float, r: float = 1.0) -> tuple:
    """
    @brief Parametrisierung der Einheitssphäre.
    @description
        Gibt kartesische Koordinaten (x, y, z) eines Punktes auf der
        Sphäre mit Radius r zurück, parametrisiert durch Längen- und
        Breitenwinkel (u = φ, v = θ).

        Formel: x = r·cos(u)·sin(v), y = r·sin(u)·sin(v), z = r·cos(v)

    @param u: Azimutalwinkel φ ∈ [0, 2π]
    @param v: Polarwinkel θ ∈ [0, π]
    @param r: Sphärenradius (Standard: 1.0)
    @return: Tuple (x, y, z)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    return (r * np.cos(u) * np.sin(v), r * np.sin(u) * np.sin(v), r * np.cos(v))


def torus_surface(u: float, v: float, R: float = 2.0, r: float = 1.0) -> tuple:
    """
    @brief Parametrisierung des Torus.
    @description
        Gibt kartesische Koordinaten eines Punktes auf dem Torus zurück.
        R = Großradius (Abstand Zentrum-Röhre), r = Kleinradius (Röhrenradius).

        Formel:
            x = (R + r·cos(v))·cos(u)
            y = (R + r·cos(v))·sin(u)
            z = r·sin(v)

    @param u: Winkel um die Hauptachse ∈ [0, 2π]
    @param v: Winkel um die Röhre ∈ [0, 2π]
    @param R: Großradius (Standard: 2.0)
    @param r: Kleinradius (Standard: 1.0)
    @return: Tuple (x, y, z)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    return ((R + r * np.cos(v)) * np.cos(u),
            (R + r * np.cos(v)) * np.sin(u),
            r * np.sin(v))


def saddle_surface(u: float, v: float) -> tuple:
    """
    @brief Sattelfläche (hyperbolisches Paraboloid) z = u² - v².
    @description
        Die Sattelfläche hat überall negative Gaußsche Krümmung K < 0,
        da sie sich in zwei Richtungen entgegengesetzt krümmt.

        Gaußsche Krümmung: K = -4/(1 + 4u² + 4v²)² (für z = u² - v²)

    @param u: Erste Parameterkoordinate
    @param v: Zweite Parameterkoordinate
    @return: Tuple (u, v, u²-v²)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    return (u, v, u ** 2 - v ** 2)


def _compute_gaussian_curvature_numerical(
    surface_func,
    u_arr: "np.ndarray",
    v_arr: "np.ndarray",
    h: float = 1e-4
) -> "np.ndarray":
    """
    @brief Berechnet die Gaußsche Krümmung K numerisch via finite Differenzen.
    @description
        Algorithmus:
        1. Partielle Ableitungen r_u, r_v, r_uu, r_vv, r_uv via zentralem
           Differenzenquotienten berechnen.
        2. Erste Fundamentalform: E = r_u·r_u, F = r_u·r_v, G = r_v·r_v
        3. Normalenvektor: n = r_u × r_v / |r_u × r_v|
        4. Zweite Fundamentalform: L = r_uu·n, M = r_uv·n, N = r_vv·n
        5. Gaußsche Krümmung: K = (LN - M²) / (EG - F²)

    @param surface_func: Parametrische Fläche (u, v) → (x, y, z)
    @param u_arr: 1D-Array der u-Werte
    @param v_arr: 1D-Array der v-Werte
    @param h: Schrittweite für finite Differenzen (Standard: 1e-4)
    @return: 2D-Array der Krümmungswerte K(u, v)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    nu = len(u_arr)
    nv = len(v_arr)
    K = np.zeros((nv, nu))

    for i, v in enumerate(v_arr):
        for j, u in enumerate(u_arr):
            def r(uu, vv):
                """Hilfsfunktion: Flächenpunkt als numpy-Array."""
                res = surface_func(uu, vv)
                return np.array(res, dtype=float)

            # Erste Ableitungen (zentrale Differenzen)
            r_u = (r(u + h, v) - r(u - h, v)) / (2 * h)
            r_v = (r(u, v + h) - r(u, v - h)) / (2 * h)

            # Zweite Ableitungen (zentrale Differenzen 2. Ordnung)
            r0 = r(u, v)
            r_uu = (r(u + h, v) - 2 * r0 + r(u - h, v)) / (h ** 2)
            r_vv = (r(u, v + h) - 2 * r0 + r(u, v - h)) / (h ** 2)
            r_uv = (r(u + h, v + h) - r(u + h, v - h)
                    - r(u - h, v + h) + r(u - h, v - h)) / (4 * h ** 2)

            # Erste Fundamentalform E, F, G
            E = np.dot(r_u, r_u)
            F = np.dot(r_u, r_v)
            G = np.dot(r_v, r_v)
            EG_F2 = E * G - F ** 2

            # Normalenvektor (normiert)
            normal = np.cross(r_u, r_v)
            norm_len = np.linalg.norm(normal)
            if norm_len < 1e-12 or abs(EG_F2) < 1e-20:
                # Singulärer Punkt (Pol etc.) → K = 0 setzen
                K[i, j] = 0.0
                continue
            normal = normal / norm_len

            # Zweite Fundamentalform L, M, N
            L = np.dot(r_uu, normal)
            M = np.dot(r_uv, normal)
            N = np.dot(r_vv, normal)

            # Gaußsche Krümmung K = (LN - M²) / (EG - F²)
            K[i, j] = (L * N - M ** 2) / EG_F2

    return K


def plot_gaussian_curvature_3d(
    surface,
    u_range: tuple = (0, 2 * np.pi),
    v_range: tuple = (0, np.pi),
    param_range: tuple = None,
    title: str = "Gaußsche Krümmung",
    resolution: int = 50,
    save_path: str = None,
    **kwargs
):
    """
    @brief Plottet eine parametrische Fläche mit Gaußscher Krümmung als Farbkodierung.
    @description
        Erzeugt einen 3D-Surface-Plot, bei dem die Farbe eines Flächenpunkts
        der Gaußschen Krümmung K entspricht:

        K > 0 (rot/warm):   elliptische Geometrie (Sphäre)
        K = 0 (weiß/mittel): flache Geometrie (Zylinder, Ebene)
        K < 0 (blau/kalt):  hyperbolische Geometrie (Sattel)

        Gaußsche Krümmung wird numerisch via finiten Differenzen berechnet.

        Unterstützte Flächen als String:
            'sphere', 'torus', 'saddle', 'hyperbolic_paraboloid'

    @param surface: Flächenname (str) oder parametrische Funktion (u,v) → (x,y,z)
    @param u_range: (u_min, u_max) für Meshgrid
    @param v_range: (v_min, v_max) für Meshgrid
    @param param_range: Falls angegeben, wird für beide Parameter verwendet
    @param title: Diagrammtitel
    @param resolution: Anzahl Gitterpunkte je Parameter (Standard: 50)
    @param save_path: Pfad zum Speichern (None → kein Speichern)
    @return: (fig, ax) - Figure und Axes3D
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Vordefinierte Flächen
    _surfaces = {
        'sphere': (sphere_surface, (0, 2 * np.pi), (0.01, np.pi - 0.01)),
        'torus': (torus_surface, (0, 2 * np.pi), (0, 2 * np.pi)),
        'saddle': (saddle_surface, (-2.0, 2.0), (-2.0, 2.0)),
        'hyperbolic_paraboloid': (saddle_surface, (-2.0, 2.0), (-2.0, 2.0)),
    }

    # Flächenfunktion bestimmen
    if isinstance(surface, str):
        surface_lower = surface.lower()
        if surface_lower not in _surfaces:
            raise ValueError(
                f"Unbekannte Fläche '{surface}'. "
                f"Bekannte Flächen: {list(_surfaces.keys())}"
            )
        surface_func, default_u, default_v = _surfaces[surface_lower]
        if param_range is not None:
            u_range_use = param_range
            v_range_use = param_range
        else:
            u_range_use = default_u
            v_range_use = default_v
    else:
        surface_func = surface
        if param_range is not None:
            u_range_use = param_range
            v_range_use = param_range
        else:
            u_range_use = u_range
            v_range_use = v_range

    # Auflösung begrenzen für Berechnungszeit
    n = min(resolution, 50)

    # Parametergitter erstellen
    u_arr = np.linspace(u_range_use[0], u_range_use[1], n)
    v_arr = np.linspace(v_range_use[0], v_range_use[1], n)

    # Fläche berechnen
    X = np.zeros((n, n))
    Y = np.zeros((n, n))
    Z = np.zeros((n, n))
    for i, v in enumerate(v_arr):
        for j, u in enumerate(u_arr):
            x, y, z = surface_func(u, v)
            X[i, j] = x
            Y[i, j] = y
            Z[i, j] = z

    # Gaußsche Krümmung numerisch berechnen
    K = _compute_gaussian_curvature_numerical(surface_func, u_arr, v_arr)

    # Farbkodierung: Krümmung normieren
    K_finite = K[np.isfinite(K)]
    if len(K_finite) > 0:
        k_max = max(abs(K_finite.max()), abs(K_finite.min()), 1e-10)
        K_norm = np.clip(K / k_max, -1, 1)
        K_norm = (K_norm + 1) / 2.0  # [0, 1] für Colormap
    else:
        K_norm = np.full_like(K, 0.5)

    # NaN-Werte auf 0.5 (neutral) setzen
    K_norm = np.where(np.isfinite(K_norm), K_norm, 0.5)

    # 3D-Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Farbkarte für Krümmung ('RdBu_r': blau = negativ, rot = positiv)
    colormap = plt.cm.RdBu_r
    facecolors = colormap(K_norm)

    ax.plot_surface(
        X, Y, Z,
        facecolors=facecolors,
        shade=True,
        alpha=0.9
    )

    # Colorbar hinzufügen
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_array(K)
    sm.set_clim(K_finite.min() if len(K_finite) > 0 else -1,
                K_finite.max() if len(K_finite) > 0 else 1)
    plt.colorbar(sm, ax=ax, label='Gaußsche Krümmung K', shrink=0.5)

    ax.set_title(title, fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig, ax


# ===========================================================================
# GEODÄTEN AUF FLÄCHEN (Build 46)
# ===========================================================================

def plot_geodesic_on_sphere(
    start_point: tuple = None,
    direction: tuple = (1.0, 0.0),
    steps: int = 200,
    r: float = 1.0,
    start_angle: tuple = None,
    n_steps: int = None,
    save_path: str = None,
) -> plt.Figure:
    """
    @brief Visualisiert eine Geodäte (Großkreis) auf der Einheitssphäre.
    @description
        Eine Geodäte auf der Sphäre ist ein Großkreis.
        Numerische Integration der Geodätengleichung:

            θ'' = sin(θ)·cos(θ)·(φ')²
            φ'' = -2·cos(θ)/sin(θ)·θ'·φ'

        Anfangsbedingungen: (θ₀, φ₀) mit Tangentialvektor (dθ/dt, dφ/dt).

    @param start_point: (theta, phi) Startpunkt in Kugelkoordinaten (Bogenmaß)
    @param direction: (dtheta, dphi) Tangentialvektor am Startpunkt
    @param steps: Anzahl Integrationsschritte
    @param r: Sphärenradius
    @param start_angle: Alternative zu start_point (für Testkompatibilität)
    @param n_steps: Alternative zu steps (für Testkompatibilität)
    @param save_path: Pfad zum Speichern
    @return: matplotlib Figure
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    from scipy.integrate import solve_ivp

    # Parameter-Aliase
    if start_angle is not None:
        start_point = start_angle
    if n_steps is not None:
        steps = n_steps
    if start_point is None:
        start_point = (np.pi / 2, 0.0)

    theta0, phi0 = float(start_point[0]), float(start_point[1])
    dtheta0, dphi0 = float(direction[0]), float(direction[1])

    # Tangentialvektor normieren
    sin_t0 = max(abs(np.sin(theta0)), 1e-8)
    v_norm = np.sqrt(dtheta0 ** 2 + (sin_t0 * dphi0) ** 2)
    if v_norm > 1e-12:
        dtheta0 /= v_norm
        dphi0 /= v_norm

    def geodesic_ode(t, y):
        """Geodätengleichung auf der Sphäre (Christoffel-Symbole)."""
        theta, phi, dtheta, dphi = y
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        if abs(sin_theta) < 1e-10:
            sin_theta = 1e-10
        ddtheta = sin_theta * cos_theta * dphi ** 2
        ddphi = -2.0 * (cos_theta / sin_theta) * dtheta * dphi
        return [dtheta, dphi, ddtheta, ddphi]

    t_span = (0, 2 * np.pi)
    t_eval = np.linspace(0, 2 * np.pi, steps)
    y0 = [theta0, phi0, dtheta0, dphi0]

    sol = solve_ivp(
        geodesic_ode, t_span, y0,
        t_eval=t_eval, method='RK45',
        rtol=1e-8, atol=1e-10,
    )

    # Kugelkoordinaten → kartesische Koordinaten
    theta_sol = sol.y[0]
    phi_sol = sol.y[1]
    x_geo = r * np.sin(theta_sol) * np.cos(phi_sol)
    y_geo = r * np.sin(theta_sol) * np.sin(phi_sol)
    z_geo = r * np.cos(theta_sol)

    # Sphäre als Hintergrund
    u_bg = np.linspace(0, 2 * np.pi, 40)
    v_bg = np.linspace(0, np.pi, 30)
    U, V = np.meshgrid(u_bg, v_bg)
    X_bg = r * np.sin(V) * np.cos(U)
    Y_bg = r * np.sin(V) * np.sin(U)
    Z_bg = r * np.cos(V)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X_bg, Y_bg, Z_bg, alpha=0.2, color='lightblue', shade=True)
    ax.plot(x_geo, y_geo, z_geo, 'r-', linewidth=2.5, label='Geodäte (Großkreis)')
    ax.scatter([x_geo[0]], [y_geo[0]], [z_geo[0]],
               color='green', s=60, zorder=5, label='Startpunkt')

    ax.set_title('Geodäte auf der Sphäre', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_geodesic_on_torus(
    R: float = 2.0,
    r: float = 1.0,
    start_params: tuple = (0.0, 0.0),
    direction: tuple = (1.0, 0.3),
    n_steps: int = 300,
    t_max: float = 10.0,
    save_path: str = None,
) -> plt.Figure:
    """
    @brief Visualisiert eine Geodäte auf dem Torus.
    @description
        Geodätengleichungen auf dem Torus (Parametrisierung mit u, v):

            u'' = 2r·sin(v)/(R + r·cos(v)) · u'·v'
            v'' = -(sin(v)·(R + r·cos(v))/r) · (u')²

        Numerische Integration via RK45.

    @param R: Großradius des Torus
    @param r: Kleinradius des Torus
    @param start_params: (u₀, v₀) Startparameter
    @param direction: (du/dt, dv/dt) Anfangsgeschwindigkeit
    @param n_steps: Anzahl Zeitschritte
    @param t_max: Endzeit der Integration
    @param save_path: Speicherpfad
    @return: matplotlib Figure
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    from scipy.integrate import solve_ivp

    u0, v0 = float(start_params[0]), float(start_params[1])
    du0, dv0 = float(direction[0]), float(direction[1])

    def torus_geodesic_ode(t, y):
        """Geodätengleichungen auf dem Torus (Christoffel-Symbole)."""
        u, v, du, dv = y
        f = R + r * np.cos(v)
        sin_v = np.sin(v)
        if abs(f) < 1e-10:
            f = 1e-10
        ddu = 2.0 * r * sin_v / f * du * dv
        ddv = -sin_v * f / r * du ** 2
        return [du, dv, ddu, ddv]

    t_eval = np.linspace(0, t_max, n_steps)
    y0 = [u0, v0, du0, dv0]

    sol = solve_ivp(
        torus_geodesic_ode, (0, t_max), y0,
        t_eval=t_eval, method='RK45',
        rtol=1e-7, atol=1e-9
    )

    u_sol = sol.y[0]
    v_sol = sol.y[1]

    # Torus-Koordinaten
    x_geo = (R + r * np.cos(v_sol)) * np.cos(u_sol)
    y_geo = (R + r * np.cos(v_sol)) * np.sin(u_sol)
    z_geo = r * np.sin(v_sol)

    # Torus-Oberfläche als Hintergrund
    u_bg = np.linspace(0, 2 * np.pi, 50)
    v_bg = np.linspace(0, 2 * np.pi, 30)
    U_bg, V_bg = np.meshgrid(u_bg, v_bg)
    X_bg = (R + r * np.cos(V_bg)) * np.cos(U_bg)
    Y_bg = (R + r * np.cos(V_bg)) * np.sin(U_bg)
    Z_bg = r * np.sin(V_bg)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X_bg, Y_bg, Z_bg, alpha=0.2, color='lightblue', shade=True)
    ax.plot(x_geo, y_geo, z_geo, 'r-', linewidth=2.0, label='Geodäte')
    ax.scatter([x_geo[0]], [y_geo[0]], [z_geo[0]],
               color='green', s=60, zorder=5, label='Startpunkt')

    ax.set_title('Geodäte auf dem Torus', fontsize=13)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


# ===========================================================================
# PDE-ANIMATIONEN (Build 46)
# ===========================================================================

def animate_heat_equation(
    alpha: float = 0.01,
    nx: int = 50,
    nt: int = 200,
    dt: float = 0.001,
    save_path: str = None,
) -> "matplotlib.animation.FuncAnimation":
    """
    @brief Animiert die Wärmegleichung u_t = α·u_xx auf [0,1].
    @description
        Löst die 1D-Wärmegleichung mit:
        - Anfangsbedingung: u(x,0) = sin(π·x)
        - Dirichlet-Randbedingungen: u(0,t) = u(1,t) = 0
        - Explizites Euler-Schema (FTCS)

        Stabilitätsbedingung: r = α·dt/dx² ≤ 0.5

        Die analytische Lösung ist: u(x,t) = e^{-α·π²·t} · sin(π·x)

    @param alpha: Wärmediffusivität (Standard: 0.01)
    @param nx: Anzahl Gitterpunkte in x (Standard: 50)
    @param nt: Anzahl Zeitschritte (Standard: 200)
    @param dt: Zeitschrittgröße (Standard: 0.001)
    @param save_path: Falls angegeben, GIF-Datei speichern
    @return: matplotlib.animation.FuncAnimation Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import matplotlib.animation as mpl_animation

    # Gitter aufbauen
    dx = 1.0 / (nx - 1)
    x = np.linspace(0, 1, nx)

    # Stabilitätsprüfung: r = α·dt/dx² ≤ 0.5
    r = alpha * dt / dx ** 2
    if r > 0.5:
        dt = 0.4 * dx ** 2 / alpha
        r = alpha * dt / dx ** 2

    # Anfangsbedingung: u(x,0) = sin(π·x)
    u_curr = np.sin(np.pi * x).copy()
    u_curr[0] = 0.0
    u_curr[-1] = 0.0

    # Zeitschritte berechnen, jeden 10. als Frame speichern
    frames = []
    for step in range(nt):
        u_new = u_curr.copy()
        u_new[1:-1] = u_curr[1:-1] + r * (
            u_curr[2:] - 2 * u_curr[1:-1] + u_curr[:-2]
        )
        u_new[0] = 0.0
        u_new[-1] = 0.0
        u_curr = u_new
        if step % 10 == 0:
            frames.append((step * dt, u_curr.copy()))

    # Animation erstellen
    fig, ax = plt.subplots(figsize=(9, 5))
    line, = ax.plot(x, frames[0][1], 'b-', linewidth=2)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title(f'Wärmegleichung: u_t = {alpha}·u_xx')
    ax.grid(True, alpha=0.3)

    def update(frame_idx):
        """Aktualisiert den Plot für jeden Frame."""
        t_val, u_frame = frames[frame_idx]
        line.set_ydata(u_frame)
        time_text.set_text(f't = {t_val:.4f}')
        return line, time_text

    anim = mpl_animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=50,
        blit=True
    )

    if save_path:
        anim.save(save_path, writer='pillow', fps=20)
    else:
        plt.close(fig)

    return anim


def animate_wave_equation_pde(
    c: float = 1.0,
    nx: int = 100,
    nt: int = 500,
    save_path: str = None,
) -> "matplotlib.animation.FuncAnimation":
    """
    @brief Animiert die Wellengleichung u_tt = c²·u_xx auf [0,1].
    @description
        Löst die 1D-Wellengleichung:
        - Anfangsbedingung: u(x,0) = sin(π·x), u_t(x,0) = 0
        - Dirichlet-Randbedingungen: u(0,t) = u(L,t) = 0
        - Explizites Leap-Frog-Schema (2. Ordnung in Zeit und Raum)

        CFL-Bedingung: c·dt/dx ≤ 1 (für Stabilität)

    @param c: Wellengeschwindigkeit (Standard: 1.0)
    @param nx: Anzahl Gitterpunkte (Standard: 100)
    @param nt: Anzahl Zeitschritte (Standard: 500)
    @param save_path: GIF-Speicherpfad
    @return: FuncAnimation Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import matplotlib.animation as mpl_animation

    dx = 1.0 / (nx - 1)
    x = np.linspace(0, 1, nx)

    # CFL-Bedingung: c·dt/dx ≤ 1
    dt = 0.9 * dx / c
    r2 = (c * dt / dx) ** 2

    # Anfangsbedingungen
    u_prev = np.sin(np.pi * x).copy()
    u_curr = u_prev.copy()
    # Erster Schritt (u_t(x,0) = 0)
    u_curr[1:-1] = (u_prev[1:-1]
                    + 0.5 * r2 * (u_prev[2:] - 2 * u_prev[1:-1] + u_prev[:-2]))
    u_curr[0] = u_curr[-1] = 0.0

    # Frames sammeln
    frames = [(0.0, u_prev.copy()), (dt, u_curr.copy())]
    for step in range(2, nt):
        u_next = np.zeros_like(u_curr)
        u_next[1:-1] = (2 * u_curr[1:-1] - u_prev[1:-1]
                        + r2 * (u_curr[2:] - 2 * u_curr[1:-1] + u_curr[:-2]))
        u_next[0] = u_next[-1] = 0.0
        u_prev = u_curr
        u_curr = u_next
        if step % 10 == 0:
            frames.append((step * dt, u_curr.copy()))

    # Animation
    fig, ax = plt.subplots(figsize=(9, 5))
    line, = ax.plot(x, frames[0][1], 'b-', linewidth=2)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title(f'Wellengleichung: u_tt = {c}²·u_xx')
    ax.grid(True, alpha=0.3)

    def update(frame_idx):
        """Aktualisiert Frame."""
        t_val, u_frame = frames[frame_idx]
        line.set_ydata(u_frame)
        time_text.set_text(f't = {t_val:.4f}')
        return line, time_text

    anim = mpl_animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=40,
        blit=True
    )

    if save_path:
        anim.save(save_path, writer='pillow', fps=25)
    else:
        plt.close(fig)

    return anim


# ===========================================================================
# SPEZIELLE FUNKTIONEN GALERIE (Build 46)
# ===========================================================================

def plot_special_functions_gallery(save_path: str = None) -> plt.Figure:
    """
    @brief Galerie der wichtigsten speziellen Funktionen der Mathematik.
    @description
        Erstellt eine Galerie-Figur mit 3 Zeilen:

        Zeile 1: Bessel-Funktionen J_0 bis J_4 auf [0, 15]
            - Lösung der Bessel-DGL: x²y'' + xy' + (x² - n²)y = 0

        Zeile 2: Legendre-Polynome P_0 bis P_4 auf [-1, 1]
            - Orthogonale Polynome auf [-1, 1]

        Zeile 3: Hermite-Polynome H_0 bis H_4 auf [-3, 3]
            - Orthogonale Polynome bzgl. Gaußmaß e^{-x²}

    @param save_path: Pfad zum Speichern (None → kein Speichern)
    @return: matplotlib Figure (3×1 Subplot-Grid)
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import scipy.special as sp

    fig, axes = plt.subplots(3, 1, figsize=(12, 14))

    # --- Zeile 1: Bessel-Funktionen J_n(x) ---
    ax1 = axes[0]
    x_bessel = np.linspace(0, 15, 400)
    colors_b = plt.cm.tab10(np.linspace(0, 0.5, 5))
    for n, color in zip(range(5), colors_b):
        y = sp.jn(n, x_bessel)
        ax1.plot(x_bessel, y, color=color, linewidth=2, label=f'J_{n}(x)')
    ax1.axhline(0, color='k', linewidth=0.5)
    ax1.set_xlim(0, 15)
    ax1.set_ylim(-0.5, 1.1)
    ax1.set_title('Bessel-Funktionen J_n(x)', fontsize=12)
    ax1.set_xlabel('x')
    ax1.set_ylabel('J_n(x)')
    ax1.legend(loc='upper right', ncol=5, fontsize=9)
    ax1.grid(True, alpha=0.3)

    # --- Zeile 2: Legendre-Polynome P_n(x) ---
    ax2 = axes[1]
    x_legendre = np.linspace(-1, 1, 300)
    colors_l = plt.cm.tab10(np.linspace(0, 0.5, 5))
    for n, color in zip(range(5), colors_l):
        P = sp.legendre(n)
        y = P(x_legendre)
        ax2.plot(x_legendre, y, color=color, linewidth=2, label=f'P_{n}(x)')
    ax2.axhline(0, color='k', linewidth=0.5)
    ax2.set_xlim(-1, 1)
    ax2.set_title('Legendre-Polynome P_n(x)', fontsize=12)
    ax2.set_xlabel('x')
    ax2.set_ylabel('P_n(x)')
    ax2.legend(loc='upper left', ncol=5, fontsize=9)
    ax2.grid(True, alpha=0.3)

    # --- Zeile 3: Hermite-Polynome H_n(x) ---
    ax3 = axes[2]
    x_hermite = np.linspace(-3, 3, 300)
    colors_h = plt.cm.tab10(np.linspace(0, 0.5, 5))
    for n, color in zip(range(5), colors_h):
        H = sp.hermite(n)
        y = H(x_hermite)
        # Skalieren für bessere Darstellung
        scale = max(abs(y).max(), 1.0)
        ax3.plot(x_hermite, y / scale, color=color, linewidth=2,
                 label=f'H_{n}(x)/|max|')
    ax3.axhline(0, color='k', linewidth=0.5)
    ax3.set_xlim(-3, 3)
    ax3.set_ylim(-1.2, 1.2)
    ax3.set_title('Hermite-Polynome H_n(x) (normiert)', fontsize=12)
    ax3.set_xlabel('x')
    ax3.set_ylabel('H_n(x) / |max|')
    ax3.legend(loc='upper right', ncol=5, fontsize=9)
    ax3.grid(True, alpha=0.3)

    plt.suptitle('Galerie spezieller Funktionen', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


# ===========================================================================
# ADAPTIVER GITTERPLOT
# ===========================================================================

def plot_adaptive_grid(
    f: Callable,
    xmin: float,
    xmax: float,
    tol: float = 1e-3,
    max_depth: int = 12,
    initial_points: int = 50,
    save_path: str = None,
    title: str = None,
    return_fig: bool = True
) -> plt.Figure:
    """
    @brief Adaptiver Gitterplot: feineres Gitter nahe Nullstellen und Singularitäten.
    @description
        Klassische äquidistante Gitter können Nullstellen oder Singularitäten
        einer Funktion übersehen. Diese Funktion verwendet eine adaptive
        Bisektionsstrategie:

        Algorithmus:
        1. Starte mit initial_points äquidistanten Punkten auf [xmin, xmax]
        2. Für jedes Intervall [a, b]: Falls die lokale Variation > tol,
           halbiere das Intervall (rekursive Bisektion)
        3. Kriterien für Verfeinerung:
           - |f(a) - f(b)| > tol (hohe Variation)
           - Vorzeichenwechsel f(a)·f(b) < 0 (Nullstelle enthalten)
           - Extremwertschätzung |f(m) - (f(a)+f(b))/2| > tol (Extremum)
        4. Maximaltiefe max_depth verhindert unendliche Rekursion

        Die resultierenden Punkte sind äquidistant in "wichtigen" Bereichen.

    @param f              Zu plottende Funktion f(x)
    @param xmin           Linke Grenze des Plotbereichs
    @param xmax           Rechte Grenze des Plotbereichs
    @param tol            Toleranz für Adaptivität (Standard: 1e-3)
    @param max_depth      Maximale Rekursionstiefe (Standard: 12)
    @param initial_points Anzahl der Startpunkte (Standard: 50)
    @param save_path      Pfad zum Speichern (None = anzeigen)
    @param title          Titel des Plots (None = automatisch)
    @param return_fig     True = Figure zurückgeben, False = None
    @return               matplotlib Figure-Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Adaptive Punktesammlung via iterativem Stack-basiertem Algorithmus
    # (Stack statt Rekursion um RecursionError bei großem max_depth zu vermeiden)
    x_points = set()

    # Hilfsfunktion: Sicherer Funktionsaufruf (Singularitäten abfangen)
    def safe_f(x: float) -> float | None:
        """Ruft f(x) sicher auf; None bei Ausnahmen (Singularitäten)."""
        try:
            val = float(f(x))
            # Inf/NaN als Singularität behandeln
            if not math.isfinite(val):
                return None
            return val
        except Exception:
            return None

    # Initiale Punktemenge: äquidistante Startpunkte
    initial_x = np.linspace(xmin, xmax, initial_points)
    x_points.update(initial_x.tolist())

    # Stack-basierte Bisektion: [(a, b, tiefe), ...]
    stack = [(initial_x[i], initial_x[i + 1], 0)
             for i in range(len(initial_x) - 1)]

    # Adaptive Verfeinerung
    while stack:
        a, b, depth = stack.pop()

        # Maximale Tiefe erreicht → keine weitere Verfeinerung
        if depth >= max_depth:
            continue

        # Mittelpunkt berechnen
        mid = (a + b) / 2.0

        # Funktionswerte an a, b, mid berechnen
        fa = safe_f(a)
        fb = safe_f(b)
        fmid = safe_f(mid)

        # Falls Singularitäten → immer verfeinern (bis max_depth)
        if fa is None or fb is None or fmid is None:
            x_points.add(mid)
            stack.append((a, mid, depth + 1))
            stack.append((mid, b, depth + 1))
            continue

        # Kriterium 1: Vorzeichenwechsel → Nullstelle im Intervall
        has_sign_change = (fa * fb < 0)

        # Kriterium 2: Lineare Interpolation weicht von Mittelpunkt ab
        # Misst die Nicht-Linearität der Funktion im Intervall
        linear_interpolation = (fa + fb) / 2.0
        has_high_curvature = abs(fmid - linear_interpolation) > tol

        # Kriterium 3: Großer absoluter Sprung zwischen den Endpunkten
        has_large_jump = abs(fb - fa) > tol * max(1.0, abs(fa), abs(fb))

        # Verfeinern falls eines der Kriterien erfüllt ist
        if has_sign_change or has_high_curvature or has_large_jump:
            # Mittelpunkt zur Punktemenge hinzufügen
            x_points.add(mid)
            # Beide Hälften auf den Stack legen
            stack.append((a, mid, depth + 1))
            stack.append((mid, b, depth + 1))

    # Punkte sortieren und Funktionswerte berechnen
    x_sorted = sorted(x_points)
    x_plot = []
    y_plot = []

    for x in x_sorted:
        y = safe_f(x)
        if y is not None:
            x_plot.append(x)
            y_plot.append(y)

    # Plot erstellen
    fig, ax = plt.subplots(figsize=(10, 6))

    # Funktion mit adaptivem Gitter plotten
    ax.plot(x_plot, y_plot, 'b-', linewidth=1.5, label='f(x)')

    # Nulllinien markieren
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)

    # Gitterpunkte als kleine Markierungen anzeigen (optional)
    ax.plot(x_plot, y_plot, 'g.', markersize=2, alpha=0.3, label=f'{len(x_plot)} Gitterpunkte')

    # Nullstellen grob markieren (Vorzeichenwechsel)
    for i in range(len(y_plot) - 1):
        if y_plot[i] * y_plot[i + 1] < 0:
            # Lineare Interpolation der Nullstelle
            x_zero = x_plot[i] - y_plot[i] * (x_plot[i + 1] - x_plot[i]) / (y_plot[i + 1] - y_plot[i])
            ax.axvline(x_zero, color='r', linewidth=0.8, linestyle=':', alpha=0.6)

    # Achsenbeschriftungen und Titel
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('f(x)', fontsize=12)
    plot_title = title or f'Adaptiver Gitterplot von f(x) auf [{xmin}, {xmax}]'
    ax.set_title(plot_title, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    elif not return_fig:
        plt.show()
        plt.close(fig)
        return None

    return fig


# ===========================================================================
# INTERAKTIVER PLOT MIT SLIDER
# ===========================================================================

def create_interactive_plot(
    f: Callable,
    xmin: float,
    xmax: float,
    param_name: str = 'a',
    param_min: float = 0.0,
    param_max: float = 2.0,
    param_init: float = 1.0,
    n_points: int = 500,
    save_path: str = None,
    title: str = None
) -> plt.Figure:
    """
    @brief Interaktiver Plot mit matplotlib.widgets.Slider für Parameter.
    @description
        Erstellt einen 2D-Plot von f(x, param) mit einem Slider-Widget,
        der den Parameter interaktiv ändert. Die Funktion wird als
        f(x, param) interpretiert wobei param der Slider-Wert ist.

        Hinweis: Der interaktive Slider funktioniert nur in interaktiven
        Matplotlib-Backends (z.B. TkAgg, Qt5Agg). Im Headless/Agg-Modus
        wird ein statischer Plot mit Slider-Bild erzeugt.

        Die Funktion f muss zwei Argumente akzeptieren: f(x, param).
        Beispiel: f = lambda x, a: a * np.sin(x)

    @param f          Funktion f(x, param) mit zwei Argumenten
    @param xmin       Linke Grenze
    @param xmax       Rechte Grenze
    @param param_name Name des Parameters (für Slider-Label)
    @param param_min  Minimalwert des Parameters
    @param param_max  Maximalwert des Parameters
    @param param_init Startwert des Parameters
    @param n_points   Anzahl der Stützpunkte
    @param save_path  Pfad zum Speichern (None = plt.show())
    @param title      Titel des Plots
    @return           matplotlib Figure-Objekt
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    from matplotlib.widgets import Slider

    # Ausgangspunkte berechnen
    x_vals = np.linspace(xmin, xmax, n_points)

    # Hilfsfunktion: f(x, param) sicher auswerten
    def eval_f(param_val: float) -> np.ndarray:
        """Berechnet f(x, param) für alle x-Werte."""
        y = np.zeros(len(x_vals))
        for i, xi in enumerate(x_vals):
            try:
                val = float(f(xi, param_val))
                y[i] = val if math.isfinite(val) else np.nan
            except Exception:
                y[i] = np.nan
        return y

    # Figure mit Platz für den Slider erstellen
    # Unteres Drittel für den Slider reservieren
    fig, ax = plt.subplots(figsize=(10, 7))
    plt.subplots_adjust(bottom=0.25)  # Platz für Slider unten

    # Initialer Plot
    y_init = eval_f(param_init)
    line, = ax.plot(x_vals, y_init, 'b-', linewidth=2, label=f'f(x, {param_name}={param_init:.2f})')
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.3)

    # Achsenbeschriftungen
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('f(x)', fontsize=12)
    ax.set_title(title or f'Interaktiver Plot f(x, {param_name})', fontsize=13)
    ax.grid(True, alpha=0.3)

    # y-Achse automatisch skalieren (Anfangswert)
    y_finite = y_init[np.isfinite(y_init)]
    if len(y_finite) > 0:
        y_range = max(abs(y_finite.max()), abs(y_finite.min())) * 1.2
        if y_range > 0:
            ax.set_ylim(-y_range, y_range)

    # Slider-Achse erstellen (unterhalb des Hauptplots)
    slider_ax = fig.add_axes([0.15, 0.08, 0.70, 0.03])
    slider = Slider(
        ax=slider_ax,
        label=param_name,
        valmin=param_min,
        valmax=param_max,
        valinit=param_init,
        color='steelblue'
    )

    # Callback-Funktion für Slider-Änderungen
    def update_plot(val: float) -> None:
        """Aktualisiert den Plot wenn der Slider bewegt wird."""
        # Neuen y-Wert berechnen
        y_new = eval_f(slider.val)
        # Plot aktualisieren
        line.set_ydata(y_new)
        line.set_label(f'f(x, {param_name}={slider.val:.3f})')
        ax.legend(fontsize=10)
        # Automatische y-Skalierung
        y_fin = y_new[np.isfinite(y_new)]
        if len(y_fin) > 0:
            yr = max(abs(y_fin.max()), abs(y_fin.min())) * 1.2
            if yr > 0:
                ax.set_ylim(-yr, yr)
        # Zeichnen neu auslösen
        fig.canvas.draw_idle()

    # Slider an Callback anschließen
    slider.on_changed(update_plot)

    # Legende anzeigen
    ax.legend(fontsize=10)

    # Speichern oder anzeigen
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    else:
        # Interaktiv nur im richtigen Backend möglich
        # Im Headless-Modus (Agg) wird plt.show() nichts anzeigen
        try:
            plt.show()
        except Exception:
            pass  # Headless-Modus: Keine Anzeige möglich

    # Figure zurückgeben (damit Tests darauf zugreifen können)
    return fig


# ===========================================================================
# EXPORT: Figure in mehrere Formate speichern (PNG, SVG, PDF)
# ===========================================================================

def export_figure(fig, filepath: str, formats: list = None, dpi: int = 150) -> dict:
    """
    @file visualization.py
    @brief Exportiert eine matplotlib-Figure in mehrere Formate gleichzeitig.
    @description
        Speichert eine matplotlib-Figure in einem oder mehreren Formaten
        (PNG, SVG, PDF) mit einem einzigen Aufruf.

        LaTeX-Kompatibilität:
        - SVG und PDF eignen sich direkt für \\includegraphics in LaTeX-Dokumenten.
          Bei pdflatex: PDF-Format nutzen.
          Bei xelatex/lualatex: SVG über \\includesvg (svg-Paket) oder PDF.
        - PNG für schnelle Vorschau oder Web-Einbindung.

        Beispiel:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot([1, 2, 3], [4, 5, 6])
            ergebnisse = export_figure(fig, "/tmp/mein_plot", formats=["pdf", "svg"])
            # ergebnisse: {"pdf": "/tmp/mein_plot.pdf", "svg": "/tmp/mein_plot.svg"}

    @param fig:      matplotlib Figure-Objekt, das gespeichert werden soll.
    @param filepath: Basispfad ohne Dateiendung (z.B. "/tmp/plot").
                     Die Dateiendung wird automatisch je Format angehängt.
    @param formats:  Liste der gewünschten Formate (Kleinbuchstaben).
                     Erlaubte Werte: "png", "svg", "pdf".
                     Standard: ["png", "svg", "pdf"] (alle drei).
    @param dpi:      Auflösung in DPI für Rasterformate (PNG). Standard: 150.
                     Vektorformate (SVG, PDF) ignorieren diesen Wert.
    @return          Dictionary mit Format als Schlüssel und absolutem Dateipfad
                     als Wert, z.B. {"png": "/tmp/plot.png", "svg": "/tmp/plot.svg"}.
    @author Michael Fuhrmann
    @date 2026-03-11
    @lastModified 2026-03-11
    """
    import os

    # Standard: alle drei gängigen Formate exportieren
    if formats is None:
        formats = ["png", "svg", "pdf"]

    # Ergebnis-Dictionary vorbereiten (Format → gespeicherter Pfad)
    saved_files = {}

    for fmt in formats:
        # Vollständigen Dateipfad mit Endung zusammensetzen
        full_path = filepath + "." + fmt

        # Figure in das jeweilige Format exportieren
        # bbox_inches="tight": Beschneidet leere Ränder automatisch
        fig.savefig(full_path, format=fmt, dpi=dpi, bbox_inches="tight")

        # Absoluten Pfad im Ergebnis-Dict speichern (plattformübergreifend)
        saved_files[fmt] = os.path.abspath(full_path)

    return saved_files


# ===========================================================================
# MASSTHEORIE: Cantor-Menge visualisieren
# ===========================================================================

def plot_cantor_set(steps: int = 6, figsize=(12, 4)) -> plt.Figure:
    """
    @file visualization.py
    @brief Visualisiert die Cantor-Menge durch iteratives Entfernen des mittleren Drittels.
    @description
        Die Cantor-Menge C ist definiert als:

            C = [0,1] \\ ⋃_{n=0}^{∞} ⋃_{k=0}^{3^n - 1} ((3k+1)/3^{n+1}, (3k+2)/3^{n+1})

        d.h. beginnend mit [0,1] wird in jedem Schritt das offene mittlere Drittel
        aus jedem verbleibenden Intervall entfernt.

        Eigenschaften der Cantor-Menge:
        - Überabzählbar (Mächtigkeit des Kontinuums)
        - Lebesgue-Maß 0 (obwohl überabzählbar!)
        - Hausdorff-Dimension: log(2)/log(3) ≈ 0.631 (echt gebrochen → Fraktal)
        - Nirgends dicht, aber abgeschlossen und perfekt (jeder Punkt ist Häufungspunkt)

        Jede Iterationsstufe wird als farbige horizontale Balkenreihe dargestellt.
        Stufe 0 = voller Einheitsbalken, Stufe k = 2^k verbleibende Intervalle.

    @param steps:   Anzahl der Iterationsschritte (sinnvoll: 0–8). Standard: 6.
                    Bei steps > 8 werden die Intervalle zu schmal zum Erkennen.
    @param figsize: Größe der Figure als (Breite, Höhe) in Zoll. Standard: (12, 4).
    @return         matplotlib Figure-Objekt mit der Visualisierung.
    @author Michael Fuhrmann
    @date 2026-03-11
    @lastModified 2026-03-11
    """
    # Startintervall: das vollständige Einheitsintervall [0, 1]
    intervals = [(0.0, 1.0)]

    # Alle Iterationsstufen vorberechnen und speichern
    # stages[0] = [(0.0, 1.0)], stages[1] = [(0.0, 1/3), (2/3, 1.0)], ...
    stages = [list(intervals)]

    for _ in range(steps):
        # Nächste Stufe: mittleres Drittel aus jedem Intervall entfernen
        next_intervals = []
        for a, b in intervals:
            # Länge des aktuellen Intervalls
            length = b - a
            # Linkes Drittel behalten: [a, a + length/3]
            next_intervals.append((a, a + length / 3.0))
            # Rechtes Drittel behalten: [b - length/3, b]
            next_intervals.append((b - length / 3.0, b))
        # Aktualisierte Intervalle für nächste Iteration übernehmen
        intervals = next_intervals
        stages.append(list(intervals))

    # Figure und Achsen anlegen
    fig, ax = plt.subplots(figsize=figsize)

    # Farbpalette: dunkler werdende Blautöne je tiefer die Stufe
    cmap = plt.get_cmap("Blues")

    for step_idx, stage_intervals in enumerate(stages):
        # Farbe: gleichmäßig über die Stufen verteilt (0.3–0.9 für guten Kontrast)
        color = cmap(0.3 + 0.6 * step_idx / max(steps, 1))

        # Intervalle als horizontale Balken darstellen
        # broken_barh erwartet: [(x_start, x_width), ...] und (y_start, y_height)
        bar_data = [(start, end - start) for start, end in stage_intervals]
        ax.broken_barh(bar_data, (step_idx - 0.4, 0.8), facecolors=color, edgecolors='none')

    # Achsenbeschriftungen und Formatierung
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.6, steps + 0.6)
    ax.set_xlabel('x ∈ [0, 1]', fontsize=12)
    ax.set_ylabel('Iterationsstufe', fontsize=12)
    # Hausdorff-Dimension der Cantor-Menge: log(2)/log(3) ≈ 0.630930
    _hausdorff_dim = math.log(2) / math.log(3)
    ax.set_title(
        f'Cantor-Menge (Mitteldrittels-Konstruktion, {steps} Schritte)\n'
        f'Hausdorff-Dimension: log(2)/log(3) ≈ {_hausdorff_dim:.6f}',
        fontsize=13
    )

    # y-Achse: Stufennummern als ganzzahlige Beschriftungen
    ax.set_yticks(range(steps + 1))
    ax.set_yticklabels([f'Stufe {i}  (2^{i}={2**i} Intervalle)' for i in range(steps + 1)],
                       fontsize=8)
    ax.grid(True, axis='x', alpha=0.3)

    fig.tight_layout()
    return fig

    return fig
