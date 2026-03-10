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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Sichere Auswertungsumgebung mit numpy-Funktionen
    _eval_env = {
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'exp': np.exp, 'log': np.log, 'log2': np.log2, 'log10': np.log10,
        'sqrt': np.sqrt, 'abs': np.abs, 'pi': np.pi, 'e': np.e,
        'arcsin': np.arcsin, 'arccos': np.arccos, 'arctan': np.arctan,
        'sinh': np.sinh, 'cosh': np.cosh, 'tanh': np.tanh,
        'np': np,
    }

    def _eval_func(x_val: np.ndarray) -> np.ndarray:
        """Wertet func_str an einem numpy-Array x_val aus."""
        env = dict(_eval_env)
        env['x'] = x_val
        try:
            result = eval(func_str, {"__builtins__": {}}, env)  # noqa: S307
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
