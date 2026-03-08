# visualization.py – Mathematische Visualisierung

## Übersicht

`visualization.py` bietet umfassende Werkzeuge zur grafischen Darstellung
mathematischer Objekte und Strukturen mittels **matplotlib** (Headless-Modus via `Agg`).

Das Modul arbeitet ohne Display (kein X-Server nötig) und speichert alle Plots
als Bilddateien (PNG, PDF, SVG etc.) wenn `save_path` angegeben wird.

---

## Technischer Aufbau

| Bibliothek | Zweck |
|---|---|
| `matplotlib (Agg)` | Headless-Rendering aller Plots |
| `numpy` | Vektorisierte Berechnungen (Gitter, Fraktale) |
| `mpl_toolkits.mplot3d` | 3D-Surface- und Kurvenplots |

---

## Allgemeines Nutzungsmuster

```python
from visualization import plot_function_2d

# Mit save_path: Bild in Datei speichern
plot_function_2d(math.sin, -3.14, 3.14, save_path='sin.png')

# Ohne save_path: plt.show() aufrufen (braucht Display)
plot_function_2d(math.sin, -3.14, 3.14)
```

Der Parameter `save_path` ist bei allen Funktionen optional.

---

## Funktionsreferenz

### 2D-Funktionsplotter

#### `plot_function_2d`
Plottet eine einzelne Funktion `y = f(x)`.

```python
plot_function_2d(f, x_min, x_max, n_points=500, title='', xlabel='x', ylabel='y', save_path=None)
```

- Punkte mit NaN oder Inf werden automatisch ausgeblendet (Polstellen-Toleranz)
- `n_points`: Anzahl äquidistanter Auswertungspunkte (Standard: 500)

#### `plot_functions_2d`
Plottet mehrere Funktionen in einem Diagramm (mit Legende).

```python
functions = [(math.sin, 'sin(x)'), (math.cos, 'cos(x)')]
plot_functions_2d(functions, x_min, x_max, n_points=500, title='', save_path=None)
```

#### `plot_parametric_2d`
Parametrische 2D-Kurve `(x(t), y(t))`.

```python
# Einheitskreis
plot_parametric_2d(math.cos, math.sin, 0, 2*math.pi, save_path='kreis.png')
```

---

### 3D-Funktionsplotter

#### `plot_function_3d`
Surface-Plot einer 2D-Funktion `z = f(x, y)`.

```python
plot_function_3d(f, x_min, x_max, y_min, y_max, n_points=50, title='', save_path=None)
```

- Colormap: `viridis` (Höhe → Farbe)
- Fehlerhafte Funktionswerte (Ausnahmen, NaN) werden mit NaN ersetzt

#### `plot_contour`
Konturplot (Isolinien) von `f(x, y)`.

```python
plot_contour(f, x_min, x_max, y_min, y_max, levels=20, title='', save_path=None)
```

- Gefüllter Konturplot (`contourf`) + überlagerte schwarze Linien
- Nützlich für Erkennung von Minima, Maxima, Sattelpunkten

#### `plot_parametric_3d`
3D-Raumkurve `(x(t), y(t), z(t))`.

```python
# Helix
plot_parametric_3d(math.cos, math.sin, lambda t: t/(2*pi), 0, 4*pi, save_path='helix.png')
```

---

### Vektorfelddarstellung

#### `plot_vector_field_2d`
2D-Vektorfeld als Pfeildiagramm (Quiver).

```python
# Rotationsfeld F = (-y, x)
plot_vector_field_2d(
    fx=lambda x, y: -y,
    fy=lambda x, y: x,
    x_min=-2, x_max=2, y_min=-2, y_max=2,
    n_grid=20, save_path='rotation.png'
)
```

- Pfeile werden normiert (einheitliche Länge), Farbe codiert den Betrag `|F|`
- `n_grid`: Gitterpunkte pro Achse (Standard: 20)

#### `plot_stream_lines`
Stromlinien (Trajektorien des Vektorfelds) via `streamplot`.

```python
plot_stream_lines(fx, fy, x_min, x_max, y_min, y_max, density=1.0, save_path=None)
```

- `density`: Dichte der Stromlinien (größer = dichter)
- Farbe codiert Feldstärke

---

### Phasenraum-Visualisierung

#### `plot_phase_portrait`
Phasenporträt eines 2D-autonomen ODE-Systems:

```
dx/dt = f(x, y)
dy/dt = g(x, y)
```

```python
# Harmonischer Oszillator: dx/dt = y, dy/dt = -x
plot_phase_portrait(
    dx_dt=lambda x, y: y,
    dy_dt=lambda x, y: -x,
    x_min=-3, x_max=3, y_min=-3, y_max=3,
    initial_conditions=[(2.0, 0.0), (1.0, 1.0)],
    t_max=10.0,
    save_path='harmonic.png'
)
```

- Vektorfeld + Trajektorien für gegebene Anfangsbedingungen
- Trajektorien via **Euler-Verfahren** (einfach, schnell)
- Anfangspunkte als grüne Punkte markiert
- Trajektorien brechen ab wenn sie den dargestellten Bereich verlassen

#### `plot_bifurcation_diagram`
Bifurkationsdiagramm `x vs r` für `x_{n+1} = f(x_n, r)`.

```python
# Logistische Abbildung
plot_bifurcation_diagram(
    f=lambda x, r: r * x * (1 - x),
    r_min=2.5, r_max=4.0,
    n_r=500, n_iter=1000, n_discard=500,
    title='Logistische Abbildung',
    save_path='bifurcation.png'
)
```

Sichtbare Strukturen:
- Fixpunkt (r < 3)
- Periodenverdopplungen (3 < r < 3.57)
- Chaos (r > 3.57, Feigenbaum-Szenario)
- Periodische Fenster im Chaos

---

### Fraktal-Generator

#### `mandelbrot_set`
Mandelbrot-Menge: `z_{n+1} = z_n² + c`, Startpunkt `z_0 = 0`.

```python
iterations = mandelbrot_set(
    x_min=-2.5, x_max=1.0, y_min=-1.25, y_max=1.25,
    width=800, height=500, max_iter=100,
    save_path='mandelbrot.png'
)
```

- **Rückgabe**: `numpy.ndarray` (height × width) mit Iterationszahlen
- Punkte in der Menge: `iterations[i,j] = 0` (nicht divergiert)
- Punkte außerhalb: Wert = Iterationszahl bei Divergenz
- Vektorisierte Berechnung (schnell)

#### `julia_set`
Julia-Menge: `z_{n+1} = z_n² + c`, festes `c`, variabler Startpunkt `z_0`.

```python
iterations = julia_set(
    c=-0.7 + 0.27j,
    x_min=-1.5, x_max=1.5, y_min=-1.5, y_max=1.5,
    width=600, height=600, max_iter=100,
    save_path='julia.png'
)
```

Interessante `c`-Werte:
- `-0.7 + 0.27j` → zusammenhängendes gebundenes Fraktal
- `0.355 + 0.355j` → dendritische Struktur
- `-1.417 + 0.107j` → Kochkurven-ähnlich

#### `sierpinski_triangle`
Sierpinski-Dreieck via **Chaos-Spiel-Methode** (Barnsley).

```python
sierpinski_triangle(n_iterations=6, save_path='sierpinski.png')
```

- `n_iterations` × 10.000 Punkte werden gezeichnet
- Startet mit zufälligem Punkt, springt iterativ zur zufälligen Ecke

Mathematischer Hintergrund:
- Das Chaos-Spiel konvergiert fast sicher (für Lebesgue-fast-alle Startpunkte)
  gegen das Sierpinski-Dreieck (ein iteriertes Funktionensystem / IFS)

#### `newton_fractal`
Newton-Fraktal für ein Polynom via Newton-Verfahren in ℂ.

```python
# z³ - 1 = 0
newton_fractal(
    f_coeffs=[1, 0, 0, -1],
    width=400, height=400,
    max_iter=50,
    save_path='newton_cubic.png'
)
```

- `f_coeffs`: Koeffizienten in absteigender Grad-Reihenfolge
- Farbe kodiert, zu welcher Nullstelle Newton konvergiert
- Grenzen zwischen Einzugsbereichen bilden das Fraktal

---

## Interne Hilfsfunktion

### `_save_or_show(fig, save_path)`
Interne Funktion (nicht öffentlich):
- `save_path` angegeben → `fig.savefig(path, dpi=150)`
- `save_path = None` → `plt.show()`
- Schließt immer die Figure (`plt.close(fig)`) zum Speicherschutz

---

## Vollständiges Beispiel

```python
import math
from visualization import (
    plot_function_2d,
    plot_functions_2d,
    mandelbrot_set,
    julia_set,
    sierpinski_triangle,
    plot_phase_portrait,
    plot_bifurcation_diagram,
)

# 1. Sinusfunktion
plot_function_2d(math.sin, -math.pi, math.pi,
                 title='Sinus', save_path='out/sin.png')

# 2. Mandelbrot-Menge
arr = mandelbrot_set(width=1920, height=1080, max_iter=200,
                     save_path='out/mandelbrot.png')
print(f"Arrayform: {arr.shape}")

# 3. Phasenporträt gedämpfter Oszillator
plot_phase_portrait(
    dx_dt=lambda x, y: y,
    dy_dt=lambda x, y: -x - 0.3*y,
    x_min=-3, x_max=3, y_min=-3, y_max=3,
    initial_conditions=[(2.5, 0), (2, 1), (-2, 0)],
    t_max=30.0,
    title='Gedämpfter Oszillator',
    save_path='out/damped.png'
)

# 4. Bifurkationsdiagramm
plot_bifurcation_diagram(
    f=lambda x, r: r * x * (1 - x),
    r_min=2.5, r_max=4.0,
    save_path='out/bifurcation.png'
)
```

---

## Hinweise zur Leistung

| Funktion | Aufwand | Tipp |
|---|---|---|
| `mandelbrot_set` | O(width × height × max_iter) | Vektorisiert – schnell |
| `julia_set` | O(width × height × max_iter) | Vektorisiert – schnell |
| `plot_function_3d` | O(n_points²) | n_points=30 für Drafts |
| `sierpinski_triangle` | O(n_iterations × 10000) | Linear – immer schnell |
| `newton_fractal` | O(width × height × max_iter) | Nutzt numpy poly1d |
| `plot_bifurcation_diagram` | O(n_r × n_iter) | Für n_r=500, n_iter=1000 ca. 1s |

---

## Autor

Kurt Ingwer – Stand: 2026-03-08
