"""
@file test_visualization.py
@brief Tests für das Visualisierungsmodul (visualization.py)
@description
    Testet alle Plot-Funktionen mit save_path in einem temporären Verzeichnis.
    Prüft dass Dateien erzeugt werden und Rückgabewerte korrekt sind.

    Da Tests headless laufen (kein Display), wird matplotlib.use('Agg') genutzt.

@author Kurt Ingwer
@date 2026-03-08
"""

import pytest
import math
import os
import tempfile
import numpy as np
import sys

# Sicherstellen dass src/ im Suchpfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import visualization as vis


# ===========================================================================
# HILFSFUNKTIONEN UND FIXTURES
# ===========================================================================

@pytest.fixture
def tmp_dir():
    """Temporäres Verzeichnis für Plot-Dateien."""
    with tempfile.TemporaryDirectory() as d:
        yield d


def make_path(tmp_dir: str, filename: str) -> str:
    """Erstellt vollständigen Pfad zur Ausgabedatei."""
    return os.path.join(tmp_dir, filename)


# Einfache Testfunktionen
def f_sin(x):
    """sin(x)"""
    return math.sin(x)


def f_cos(x):
    """cos(x)"""
    return math.cos(x)


def f_quadratic(x, y):
    """x² + y²"""
    return x**2 + y**2


def f_saddle(x, y):
    """x² - y² (Sattelform)"""
    return x**2 - y**2


# ===========================================================================
# TESTS: 2D-FUNKTIONSPLOTTER
# ===========================================================================

class TestPlotFunction2D:
    """Tests für plot_function_2d."""

    def test_creates_file(self, tmp_dir):
        """Plot muss eine Datei erzeugen."""
        path = make_path(tmp_dir, 'sin_2d.png')
        vis.plot_function_2d(f_sin, -math.pi, math.pi, save_path=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_custom_points(self, tmp_dir):
        """n_points-Parameter muss akzeptiert werden."""
        path = make_path(tmp_dir, 'sin_100.png')
        vis.plot_function_2d(f_sin, 0, 2 * math.pi, n_points=100, save_path=path)
        assert os.path.exists(path)

    def test_with_title(self, tmp_dir):
        """Titel-Parameter muss akzeptiert werden."""
        path = make_path(tmp_dir, 'sin_title.png')
        vis.plot_function_2d(f_sin, 0, math.pi, title='Sinusfunktion', save_path=path)
        assert os.path.exists(path)

    def test_lambda_function(self, tmp_dir):
        """Lambda-Funktionen müssen funktionieren."""
        path = make_path(tmp_dir, 'lambda_2d.png')
        vis.plot_function_2d(lambda x: x**2, -2, 2, save_path=path)
        assert os.path.exists(path)

    def test_function_with_infinity(self, tmp_dir):
        """Funktionen mit Polstellen dürfen nicht crashen."""
        path = make_path(tmp_dir, 'tan_2d.png')

        def safe_tan(x):
            """tan(x) – kann ±∞ sein"""
            return math.tan(x)

        # Darf nicht werfen (Polstellen werden maskiert)
        vis.plot_function_2d(safe_tan, 0, math.pi, save_path=path)
        assert os.path.exists(path)


class TestPlotFunctions2D:
    """Tests für plot_functions_2d (mehrere Funktionen)."""

    def test_creates_file(self, tmp_dir):
        """Plot mit zwei Funktionen muss Datei erzeugen."""
        path = make_path(tmp_dir, 'multi_2d.png')
        functions = [(f_sin, 'sin(x)'), (f_cos, 'cos(x)')]
        vis.plot_functions_2d(functions, -math.pi, math.pi, save_path=path)
        assert os.path.exists(path)

    def test_single_function(self, tmp_dir):
        """Auch eine einzelne Funktion muss funktionieren."""
        path = make_path(tmp_dir, 'single_multi.png')
        vis.plot_functions_2d([(f_sin, 'sin')], 0, 1, save_path=path)
        assert os.path.exists(path)

    def test_three_functions(self, tmp_dir):
        """Drei Funktionen gleichzeitig plotten."""
        path = make_path(tmp_dir, 'triple_2d.png')
        fns = [
            (lambda x: x, 'linear'),
            (lambda x: x**2, 'quadratisch'),
            (lambda x: x**3, 'kubisch'),
        ]
        vis.plot_functions_2d(fns, -1, 1, save_path=path)
        assert os.path.exists(path)


class TestPlotParametric2D:
    """Tests für parametrische 2D-Kurven."""

    def test_circle(self, tmp_dir):
        """Einheitskreis: x(t)=cos(t), y(t)=sin(t)."""
        path = make_path(tmp_dir, 'circle_2d.png')
        vis.plot_parametric_2d(math.cos, math.sin, 0, 2 * math.pi, save_path=path)
        assert os.path.exists(path)

    def test_lissajous(self, tmp_dir):
        """Lissajous-Figur."""
        path = make_path(tmp_dir, 'lissajous.png')
        vis.plot_parametric_2d(
            lambda t: math.sin(2 * t),
            lambda t: math.sin(3 * t),
            0, 2 * math.pi,
            n_points=1000,
            title='Lissajous',
            save_path=path
        )
        assert os.path.exists(path)


# ===========================================================================
# TESTS: 3D-FUNKTIONSPLOTTER
# ===========================================================================

class TestPlotFunction3D:
    """Tests für 3D-Surface-Plots."""

    def test_creates_file(self, tmp_dir):
        """3D-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'quadratic_3d.png')
        vis.plot_function_3d(f_quadratic, -2, 2, -2, 2, n_points=20, save_path=path)
        assert os.path.exists(path)

    def test_saddle(self, tmp_dir):
        """Sattelpunkt-Funktion."""
        path = make_path(tmp_dir, 'saddle_3d.png')
        vis.plot_function_3d(f_saddle, -2, 2, -2, 2, n_points=15, save_path=path)
        assert os.path.exists(path)

    def test_with_title(self, tmp_dir):
        """Titel muss akzeptiert werden."""
        path = make_path(tmp_dir, '3d_titled.png')
        vis.plot_function_3d(
            lambda x, y: math.sin(x) * math.cos(y),
            -math.pi, math.pi, -math.pi, math.pi,
            n_points=20, title='sin(x)·cos(y)', save_path=path
        )
        assert os.path.exists(path)


class TestPlotContour:
    """Tests für Konturplots."""

    def test_creates_file(self, tmp_dir):
        """Konturplot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'contour.png')
        vis.plot_contour(f_quadratic, -2, 2, -2, 2, levels=10, save_path=path)
        assert os.path.exists(path)

    def test_saddle_contour(self, tmp_dir):
        """Konturplot des Sattels."""
        path = make_path(tmp_dir, 'saddle_contour.png')
        vis.plot_contour(f_saddle, -3, 3, -3, 3, levels=15, save_path=path)
        assert os.path.exists(path)


class TestPlotParametric3D:
    """Tests für parametrische 3D-Kurven."""

    def test_helix(self, tmp_dir):
        """Schraubenlinie (Helix)."""
        path = make_path(tmp_dir, 'helix_3d.png')
        vis.plot_parametric_3d(
            math.cos,
            math.sin,
            lambda t: t / (2 * math.pi),
            0, 4 * math.pi,
            n_points=300,
            title='Helix',
            save_path=path
        )
        assert os.path.exists(path)

    def test_circle_3d(self, tmp_dir):
        """3D-Kreis (z=0)."""
        path = make_path(tmp_dir, 'circle_3d.png')
        vis.plot_parametric_3d(
            math.cos, math.sin, lambda t: 0.0,
            0, 2 * math.pi, save_path=path
        )
        assert os.path.exists(path)


# ===========================================================================
# TESTS: VEKTORFELDDARSTELLUNG
# ===========================================================================

class TestPlotVectorField2D:
    """Tests für Quiver-Plots."""

    def test_creates_file(self, tmp_dir):
        """Vektorfeld-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'vector_field.png')
        vis.plot_vector_field_2d(
            fx=lambda x, y: -y,
            fy=lambda x, y: x,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            n_grid=10, save_path=path
        )
        assert os.path.exists(path)

    def test_gradient_field(self, tmp_dir):
        """Gradientenfeld von f(x,y) = x² + y²."""
        path = make_path(tmp_dir, 'gradient_field.png')
        vis.plot_vector_field_2d(
            fx=lambda x, y: 2 * x,
            fy=lambda x, y: 2 * y,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            n_grid=15, title='Gradientenfeld',
            save_path=path
        )
        assert os.path.exists(path)


class TestPlotStreamLines:
    """Tests für Stromlinien."""

    def test_creates_file(self, tmp_dir):
        """Stromlinien-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'streamlines.png')
        vis.plot_stream_lines(
            fx=lambda x, y: 1.0,
            fy=lambda x, y: 0.0,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            density=1.0, save_path=path
        )
        assert os.path.exists(path)

    def test_rotational_field(self, tmp_dir):
        """Rotationsfeld."""
        path = make_path(tmp_dir, 'rotation_stream.png')
        vis.plot_stream_lines(
            fx=lambda x, y: -y,
            fy=lambda x, y: x,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            density=1.5, title='Rotation',
            save_path=path
        )
        assert os.path.exists(path)


# ===========================================================================
# TESTS: PHASENRAUM-VISUALISIERUNG
# ===========================================================================

class TestPlotPhasePortrait:
    """Tests für Phasenporträts."""

    def test_creates_file_without_trajectories(self, tmp_dir):
        """Phasenporträt ohne Trajektorien muss Datei erzeugen."""
        path = make_path(tmp_dir, 'phase_no_traj.png')
        vis.plot_phase_portrait(
            dx_dt=lambda x, y: -y,
            dy_dt=lambda x, y: x,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            save_path=path
        )
        assert os.path.exists(path)

    def test_creates_file_with_trajectories(self, tmp_dir):
        """Phasenporträt mit Trajektorien muss Datei erzeugen."""
        path = make_path(tmp_dir, 'phase_with_traj.png')
        vis.plot_phase_portrait(
            dx_dt=lambda x, y: -y,
            dy_dt=lambda x, y: x,
            x_min=-2, x_max=2, y_min=-2, y_max=2,
            initial_conditions=[(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)],
            t_max=5.0,
            save_path=path
        )
        assert os.path.exists(path)

    def test_damped_oscillator(self, tmp_dir):
        """Gedämpfter Oszillator (stabiler Fixpunkt)."""
        path = make_path(tmp_dir, 'damped_osc.png')
        vis.plot_phase_portrait(
            dx_dt=lambda x, y: y,
            dy_dt=lambda x, y: -x - 0.5 * y,
            x_min=-3, x_max=3, y_min=-3, y_max=3,
            initial_conditions=[(2.0, 0.0), (1.5, 1.0)],
            t_max=10.0,
            title='Gedämpfter Oszillator',
            save_path=path
        )
        assert os.path.exists(path)


class TestPlotBifurcationDiagram:
    """Tests für Bifurkationsdiagramm."""

    def test_logistic_map(self, tmp_dir):
        """Bifurkationsdiagramm der logistischen Abbildung."""
        path = make_path(tmp_dir, 'bifurcation.png')
        # Logistische Abbildung: f(x, r) = r * x * (1 - x)
        vis.plot_bifurcation_diagram(
            f=lambda x, r: r * x * (1 - x),
            r_min=2.5, r_max=4.0,
            n_r=100,       # Wenige Werte für schnelle Tests
            n_iter=200,
            n_discard=100,
            title='Logistische Abbildung',
            save_path=path
        )
        assert os.path.exists(path)


# ===========================================================================
# TESTS: FRAKTAL-GENERATOR
# ===========================================================================

class TestMandelbrotSet:
    """Tests für die Mandelbrot-Menge."""

    def test_returns_numpy_array(self):
        """mandelbrot_set muss numpy.ndarray zurückgeben."""
        result = mandelbrot_set_small()
        assert isinstance(result, np.ndarray)

    def test_correct_shape(self):
        """Array muss korrekte Form haben (height × width)."""
        width, height = 100, 60
        result = vis.mandelbrot_set(width=width, height=height, max_iter=20)
        assert result.shape == (height, width)

    def test_values_in_range(self):
        """Iterationswerte müssen zwischen 0 und max_iter liegen."""
        max_iter = 30
        result = vis.mandelbrot_set(width=50, height=50, max_iter=max_iter)
        assert result.min() >= 0
        assert result.max() <= max_iter

    def test_creates_file(self, tmp_dir):
        """Mandelbrot-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'mandelbrot.png')
        vis.mandelbrot_set(width=100, height=80, max_iter=20, save_path=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_origin_in_mandelbrot_set(self):
        """Der Ursprung c=0 muss in der Mandelbrot-Menge sein (Iteration = 0)."""
        # Bei c=0: z → 0 → 0 → ... bleibt beschränkt → iterations[center] = 0
        result = vis.mandelbrot_set(
            x_min=-0.1, x_max=0.1, y_min=-0.1, y_max=0.1,
            width=11, height=11, max_iter=50
        )
        # Mittlerer Pixel (Index 5,5) sollte c≈0 entsprechen → nicht divergiert (=0)
        assert result[5, 5] == 0

    def test_far_outside_diverges_immediately(self):
        """Weit außerhalb der Menge liegende Punkte divergieren sofort."""
        result = vis.mandelbrot_set(
            x_min=10.0, x_max=11.0, y_min=10.0, y_max=11.0,
            width=10, height=10, max_iter=50
        )
        # Alle Punkte bei c=10+10j divergieren in der 1. Iteration
        assert result.max() <= 2  # Divergiert sehr früh


def mandelbrot_set_small(tmp_path=None):
    """Hilfsfunktion: Kleine Mandelbrot-Menge ohne sichtbaren Plot."""
    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
        path = f.name
    try:
        result = vis.mandelbrot_set(width=50, height=40, max_iter=10, save_path=path)
    finally:
        if os.path.exists(path):
            os.unlink(path)
    return result


class TestJuliaSet:
    """Tests für Julia-Mengen."""

    def test_returns_numpy_array(self, tmp_dir):
        """julia_set muss numpy.ndarray zurückgeben."""
        result = vis.julia_set(width=50, height=50, max_iter=20,
                               save_path=make_path(tmp_dir, 'julia_arr.png'))
        assert isinstance(result, np.ndarray)

    def test_correct_shape(self, tmp_dir):
        """Array muss korrekte Form haben."""
        width, height = 80, 60
        result = vis.julia_set(
            width=width, height=height, max_iter=20,
            save_path=make_path(tmp_dir, 'julia_shape.png')
        )
        assert result.shape == (height, width)

    def test_creates_file(self, tmp_dir):
        """Julia-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'julia.png')
        vis.julia_set(c=-0.7 + 0.27j, width=80, height=80, max_iter=30, save_path=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_different_c_values(self, tmp_dir):
        """Verschiedene c-Werte müssen ohne Fehler funktionieren."""
        c_values = [0 + 1j, -1.0 + 0j, 0.355 + 0.355j]
        for i, c in enumerate(c_values):
            path = make_path(tmp_dir, f'julia_c{i}.png')
            result = vis.julia_set(c=c, width=40, height=40, max_iter=15, save_path=path)
            assert isinstance(result, np.ndarray)
            assert os.path.exists(path)

    def test_values_in_range(self, tmp_dir):
        """Iterationswerte müssen zwischen 0 und max_iter liegen."""
        max_iter = 25
        result = vis.julia_set(
            width=50, height=50, max_iter=max_iter,
            save_path=make_path(tmp_dir, 'julia_range.png')
        )
        assert result.min() >= 0
        assert result.max() <= max_iter


class TestSierpinskiTriangle:
    """Tests für das Sierpinski-Dreieck."""

    def test_no_exception(self, tmp_dir):
        """Sierpinski-Plot darf keine Exception werfen."""
        path = make_path(tmp_dir, 'sierpinski.png')
        vis.sierpinski_triangle(n_iterations=1, save_path=path)  # Wenige Iterationen für Test

    def test_creates_file(self, tmp_dir):
        """Sierpinski-Plot muss Datei erzeugen."""
        path = make_path(tmp_dir, 'sierpinski_file.png')
        vis.sierpinski_triangle(n_iterations=1, save_path=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_more_iterations(self, tmp_dir):
        """Auch mit mehr Iterationen darf kein Fehler auftreten."""
        path = make_path(tmp_dir, 'sierpinski_more.png')
        vis.sierpinski_triangle(n_iterations=3, save_path=path)
        assert os.path.exists(path)


class TestNewtonFractal:
    """Tests für Newton-Fraktale."""

    def test_cubic_polynomial(self, tmp_dir):
        """Newton-Fraktal für z³ - 1."""
        path = make_path(tmp_dir, 'newton_cubic.png')
        # z³ - 1 → Koeffizienten: [1, 0, 0, -1]
        vis.newton_fractal(
            f_coeffs=[1, 0, 0, -1],
            width=80, height=80,
            max_iter=20,
            save_path=path
        )
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0

    def test_quadratic_polynomial(self, tmp_dir):
        """Newton-Fraktal für z² - 1."""
        path = make_path(tmp_dir, 'newton_quadratic.png')
        vis.newton_fractal(
            f_coeffs=[1, 0, -1],  # z² - 1
            width=60, height=60,
            max_iter=15,
            save_path=path
        )
        assert os.path.exists(path)

    def test_degree_4_polynomial(self, tmp_dir):
        """Newton-Fraktal für z⁴ - 1."""
        path = make_path(tmp_dir, 'newton_degree4.png')
        vis.newton_fractal(
            f_coeffs=[1, 0, 0, 0, -1],  # z⁴ - 1
            width=70, height=70,
            max_iter=25,
            save_path=path
        )
        assert os.path.exists(path)

    def test_custom_region(self, tmp_dir):
        """Benutzerdefinierter Ausschnitt der komplexen Ebene."""
        path = make_path(tmp_dir, 'newton_zoom.png')
        vis.newton_fractal(
            f_coeffs=[1, 0, 0, -1],
            width=60, height=60,
            x_min=-1.5, x_max=1.5, y_min=-1.5, y_max=1.5,
            max_iter=30,
            save_path=path
        )
        assert os.path.exists(path)


# ===========================================================================
# EDGE CASES
# ===========================================================================

class TestEdgeCases:
    """Tests für Randfälle und Sondersituationen."""

    def test_constant_function_2d(self, tmp_dir):
        """Konstante Funktion (horizontale Linie)."""
        path = make_path(tmp_dir, 'constant.png')
        vis.plot_function_2d(lambda x: 5.0, -1, 1, save_path=path)
        assert os.path.exists(path)

    def test_very_small_interval(self, tmp_dir):
        """Sehr kleines Intervall darf nicht crashen."""
        path = make_path(tmp_dir, 'tiny_interval.png')
        vis.plot_function_2d(f_sin, 0.0, 0.001, save_path=path)
        assert os.path.exists(path)

    def test_mandelbrot_minimal_resolution(self, tmp_dir):
        """Mandelbrot mit minimaler Auflösung."""
        path = make_path(tmp_dir, 'mandelbrot_tiny.png')
        result = vis.mandelbrot_set(
            width=5, height=5, max_iter=5,
            save_path=path
        )
        assert result.shape == (5, 5)

    def test_julia_default_params(self, tmp_dir):
        """Julia-Set mit Standard-Parametern."""
        path = make_path(tmp_dir, 'julia_default.png')
        result = vis.julia_set(save_path=path)
        assert isinstance(result, np.ndarray)
        assert result.shape == (600, 600)

    def test_vector_field_zero_field(self, tmp_dir):
        """Nullvektorfeld (alle Vektoren = 0)."""
        path = make_path(tmp_dir, 'zero_field.png')
        vis.plot_vector_field_2d(
            fx=lambda x, y: 0.0,
            fy=lambda x, y: 0.0,
            x_min=-1, x_max=1, y_min=-1, y_max=1,
            n_grid=5,
            save_path=path
        )
        assert os.path.exists(path)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
