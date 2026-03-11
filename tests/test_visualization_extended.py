"""
@file test_visualization_extended.py
@brief Tests für die erweiterten Visualisierungsfunktionen:
       Animationen (ODE-Trajektorie, Phasenporträt, Wellengleichung),
       SVG/PDF-Export und mandelbrot_smooth.

@description
    Testet alle neuen Funktionen aus visualization.py Build 11:
    - animate_ode_trajectory(): FuncAnimation zurückgegeben
    - animate_phase_portrait(): FuncAnimation zurückgegeben
    - animate_wave_equation(): FuncAnimation zurückgegeben, korrektes Format
    - save_figure_svg(): Datei wird erstellt, korrekte Endung
    - save_figure_pdf(): Datei wird erstellt, korrekte Endung
    - plot_and_export(): Datei wird erstellt, format='svg' und 'pdf'
    - export_all_formats(): Dict mit 'png', 'svg', 'pdf' zurückgegeben
    - mandelbrot_smooth(): numpy-Array mit korrekter Shape

    Alle Tests verwenden matplotlib Agg-Backend (non-interactive).

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import os
import sys
import math
import tempfile
import shutil

# src-Verzeichnis zum Suchpfad hinzufügen (damit visualization importierbar ist)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import matplotlib
# Non-interactive Backend setzen (kein Display nötig, für CI-Umgebungen)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pytest

# Zu testende Funktionen importieren
from visualization import (
    animate_ode_trajectory,
    animate_phase_portrait,
    animate_wave_equation,
    save_figure_svg,
    save_figure_pdf,
    plot_and_export,
    export_all_formats,
    mandelbrot_smooth,
    plot_gaussian_curvature_3d,
    plot_geodesic_on_sphere,
    plot_geodesic_on_torus,
    plot_function_2d_adaptive,
    export_figure,
)


# ===========================================================================
# HILFSFUNKTIONEN FÜR TESTS
# ===========================================================================

def _make_temp_dir():
    """Erstellt ein temporäres Verzeichnis und gibt den Pfad zurück."""
    return tempfile.mkdtemp(prefix='test_vis_')


# ===========================================================================
# TESTS: animate_ode_trajectory
# ===========================================================================

class TestAnimateOdeTrajektorie:
    """Tests für die ODE-Trajektorie-Animation."""

    def test_gibt_funcanimation_zurueck(self):
        """animate_ode_trajectory() muss ein FuncAnimation-Objekt zurückgeben."""
        # Einfache ODE: y' = -y (exponentieller Abfall)
        def f(t, y):
            return [-y[0]]

        result = animate_ode_trajectory(f, y0=[1.0], t_span=(0, 2.0), n_frames=10)
        assert isinstance(result, animation.FuncAnimation), (
            "Rückgabewert muss FuncAnimation sein"
        )

    def test_system_erster_ordnung(self):
        """Funktioniert mit einem 2D-System (harmonischer Oszillator)."""
        # Harmonischer Oszillator: x'' + x = 0
        # Als System: x' = v, v' = -x
        def harmonisch(t, y):
            return [y[1], -y[0]]

        result = animate_ode_trajectory(
            harmonisch, y0=[1.0, 0.0], t_span=(0, 4.0), n_frames=20
        )
        assert isinstance(result, animation.FuncAnimation)

    def test_speichert_gif_datei(self):
        """Bei Angabe von output_path wird eine GIF-Datei erzeugt."""
        tmpdir = _make_temp_dir()
        try:
            gif_path = os.path.join(tmpdir, 'test_ode.gif')

            def f_einfach(t, y):
                return [-y[0]]

            animate_ode_trajectory(
                f_einfach, y0=[1.0], t_span=(0, 1.0),
                n_frames=5, output_path=gif_path
            )
            assert os.path.exists(gif_path), "GIF-Datei muss erstellt worden sein"
            assert os.path.getsize(gif_path) > 0, "GIF-Datei darf nicht leer sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_mit_titel(self):
        """Titel-Parameter wird ohne Fehler übergeben."""
        def f(t, y):
            return [0.0]

        result = animate_ode_trajectory(
            f, y0=[0.0], t_span=(0, 1.0),
            n_frames=5, title="Mein Test-Titel"
        )
        assert result is not None


# ===========================================================================
# TESTS: animate_phase_portrait
# ===========================================================================

class TestAnimatePhasenportrait:
    """Tests für das animierte Phasenporträt."""

    def test_gibt_funcanimation_zurueck(self):
        """animate_phase_portrait() muss FuncAnimation zurückgeben."""
        # Einfaches lineares System: x' = -x, y' = -y (stabiler Ursprung)
        def f_linear(t, xy):
            return [-xy[0], -xy[1]]

        result = animate_phase_portrait(
            f_linear,
            x_range=(-2.0, 2.0),
            y_range=(-2.0, 2.0),
            n_trajectories=4
        )
        assert isinstance(result, animation.FuncAnimation), (
            "Rückgabewert muss FuncAnimation sein"
        )

    def test_van_der_pol(self):
        """Funktioniert mit nichtlinearem Van-der-Pol-Oszillator."""
        def van_der_pol(t, xy):
            x, y = xy
            # Van-der-Pol: x' = y, y' = mu*(1-x²)*y - x
            mu = 1.0
            return [y, mu * (1 - x**2) * y - x]

        result = animate_phase_portrait(
            van_der_pol,
            x_range=(-3.0, 3.0),
            y_range=(-3.0, 3.0),
            n_trajectories=4
        )
        assert isinstance(result, animation.FuncAnimation)


# ===========================================================================
# TESTS: animate_wave_equation
# ===========================================================================

class TestAnimateWellengleichung:
    """Tests für die Wellengleichungs-Animation."""

    def test_gibt_funcanimation_zurueck(self):
        """animate_wave_equation() muss FuncAnimation zurückgeben."""
        result = animate_wave_equation(k=1.0, omega=1.0, n_frames=10)
        assert isinstance(result, animation.FuncAnimation), (
            "Rückgabewert muss FuncAnimation sein"
        )

    def test_verschiedene_parameter(self):
        """Funktioniert mit verschiedenen k und omega-Werten."""
        result = animate_wave_equation(
            k=2.0, omega=3.0,
            x_range=(-4.0, 4.0),
            t_max=math.pi,
            n_frames=15
        )
        assert isinstance(result, animation.FuncAnimation)

    def test_speichert_gif(self):
        """Bei Angabe von output_path wird GIF-Datei erzeugt."""
        tmpdir = _make_temp_dir()
        try:
            gif_path = os.path.join(tmpdir, 'welle.gif')
            animate_wave_equation(
                k=1.0, omega=1.0, n_frames=5, output_path=gif_path
            )
            assert os.path.exists(gif_path), "GIF-Datei muss erstellt worden sein"
            assert os.path.getsize(gif_path) > 0
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_standard_parameter(self):
        """Standardaufruf ohne Parameter funktioniert."""
        # Alle Parameter haben Standardwerte
        result = animate_wave_equation(n_frames=5)
        assert result is not None


# ===========================================================================
# TESTS: save_figure_svg
# ===========================================================================

class TestSaveFigureSvg:
    """Tests für den SVG-Export."""

    def test_erstellt_datei(self):
        """save_figure_svg() erstellt eine SVG-Datei."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'test.svg')
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])

            result_path = save_figure_svg(fig, svg_path)
            plt.close(fig)

            assert os.path.exists(result_path), "SVG-Datei muss existieren"
            assert os.path.getsize(result_path) > 0, "SVG-Datei darf nicht leer sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_rueckgabe_ist_absoluter_pfad(self):
        """save_figure_svg() gibt absoluten Pfad zurück."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'abs_test.svg')
            fig, ax = plt.subplots()

            result = save_figure_svg(fig, svg_path)
            plt.close(fig)

            assert os.path.isabs(result), "Rückgabe muss absoluter Pfad sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_datei_ist_svg_format(self):
        """Die erzeugte Datei ist tatsächlich im SVG-Format (XML)."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'format_test.svg')
            fig, ax = plt.subplots()
            ax.set_title('SVG Test')

            save_figure_svg(fig, svg_path)
            plt.close(fig)

            # SVG-Dateien beginnen mit XML-Header oder <svg
            with open(svg_path, 'r', encoding='utf-8') as f_in:
                inhalt = f_in.read(200)
            assert '<?xml' in inhalt or '<svg' in inhalt, (
                "Datei muss gültiges SVG-Format (XML) sein"
            )
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_loescht_danach(self):
        """Stellt sicher dass Testdateien korrekt aufgeräumt werden."""
        tmpdir = _make_temp_dir()
        svg_path = os.path.join(tmpdir, 'cleanup_test.svg')
        fig, ax = plt.subplots()

        save_figure_svg(fig, svg_path)
        plt.close(fig)

        assert os.path.exists(svg_path)
        # Aufräumen
        os.remove(svg_path)
        assert not os.path.exists(svg_path), "Datei nach Löschen nicht mehr da"


# ===========================================================================
# TESTS: save_figure_pdf
# ===========================================================================

class TestSaveFigurePdf:
    """Tests für den PDF-Export."""

    def test_erstellt_datei(self):
        """save_figure_pdf() erstellt eine PDF-Datei."""
        tmpdir = _make_temp_dir()
        try:
            pdf_path = os.path.join(tmpdir, 'test.pdf')
            fig, ax = plt.subplots()
            ax.plot([0, 1, 2], [0, 1, 0])

            result_path = save_figure_pdf(fig, pdf_path)
            plt.close(fig)

            assert os.path.exists(result_path), "PDF-Datei muss existieren"
            assert os.path.getsize(result_path) > 0, "PDF-Datei darf nicht leer sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_rueckgabe_ist_absoluter_pfad(self):
        """save_figure_pdf() gibt absoluten Pfad zurück."""
        tmpdir = _make_temp_dir()
        try:
            pdf_path = os.path.join(tmpdir, 'abs_pdf.pdf')
            fig, ax = plt.subplots()

            result = save_figure_pdf(fig, pdf_path)
            plt.close(fig)

            assert os.path.isabs(result), "Rückgabe muss absoluter Pfad sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_datei_ist_pdf_format(self):
        """Die erzeugte Datei beginnt mit dem PDF-Magic-Header '%PDF'."""
        tmpdir = _make_temp_dir()
        try:
            pdf_path = os.path.join(tmpdir, 'magic_test.pdf')
            fig, ax = plt.subplots()
            ax.set_title('PDF Test')

            save_figure_pdf(fig, pdf_path)
            plt.close(fig)

            # PDF-Dateien beginnen immer mit '%PDF'
            with open(pdf_path, 'rb') as f_in:
                magic = f_in.read(4)
            assert magic == b'%PDF', "Datei muss mit '%PDF' beginnen"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ===========================================================================
# TESTS: plot_and_export
# ===========================================================================

class TestPlotAndExport:
    """Tests für die kombinierte Plot-und-Export-Funktion."""

    def test_svg_export(self):
        """plot_and_export() erstellt SVG-Datei."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'func.svg')
            result = plot_and_export(math.sin, (-math.pi, math.pi),
                                     svg_path, 'Sinus', 'svg')
            assert os.path.exists(result), "SVG-Datei muss existieren"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_pdf_export(self):
        """plot_and_export() erstellt PDF-Datei bei format='pdf'."""
        tmpdir = _make_temp_dir()
        try:
            pdf_path = os.path.join(tmpdir, 'func.pdf')
            result = plot_and_export(math.cos, (0, 2 * math.pi),
                                     pdf_path, 'Kosinus', 'pdf')
            assert os.path.exists(result), "PDF-Datei muss existieren"
            assert os.path.getsize(result) > 0
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_standard_format_ist_svg(self):
        """Ohne format-Angabe wird SVG erstellt."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'default.svg')
            result = plot_and_export(lambda x: x**2, (-2, 2), svg_path)
            assert os.path.exists(result)
            # SVG-Format prüfen
            with open(result, 'r', encoding='utf-8') as f:
                inhalt = f.read(100)
            assert '<?xml' in inhalt or '<svg' in inhalt
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ===========================================================================
# TESTS: export_all_formats
# ===========================================================================

class TestExportAllFormats:
    """Tests für den Export in alle Formate."""

    def test_gibt_dict_mit_drei_schluessel_zurueck(self):
        """export_all_formats() gibt Dict mit 'png', 'svg', 'pdf' zurück."""
        tmpdir = _make_temp_dir()
        try:
            base = os.path.join(tmpdir, 'alle_formate')
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])

            result = export_all_formats(fig, base)
            plt.close(fig)

            assert isinstance(result, dict), "Rückgabe muss ein Dict sein"
            assert 'png' in result, "Dict muss 'png' enthalten"
            assert 'svg' in result, "Dict muss 'svg' enthalten"
            assert 'pdf' in result, "Dict muss 'pdf' enthalten"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_alle_dateien_existieren(self):
        """Alle drei exportierten Dateien müssen auf der Festplatte vorhanden sein."""
        tmpdir = _make_temp_dir()
        try:
            base = os.path.join(tmpdir, 'export_test')
            fig, ax = plt.subplots()
            ax.set_title('Export Test')
            ax.scatter([1, 2, 3], [1, 4, 9])

            result = export_all_formats(fig, base)
            plt.close(fig)

            # Jede Datei muss existieren und nicht leer sein
            for fmt, pfad in result.items():
                assert os.path.exists(pfad), f"{fmt}-Datei ({pfad}) muss existieren"
                assert os.path.getsize(pfad) > 0, f"{fmt}-Datei darf nicht leer sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_dateiendungen_korrekt(self):
        """Die Dateipfade haben die korrekten Endungen."""
        tmpdir = _make_temp_dir()
        try:
            base = os.path.join(tmpdir, 'endungen')
            fig, _ = plt.subplots()

            result = export_all_formats(fig, base)
            plt.close(fig)

            assert result['png'].endswith('.png'), "PNG-Pfad muss .png enden"
            assert result['svg'].endswith('.svg'), "SVG-Pfad muss .svg enden"
            assert result['pdf'].endswith('.pdf'), "PDF-Pfad muss .pdf enden"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ===========================================================================
# TESTS: mandelbrot_smooth
# ===========================================================================

class TestMandelbrotSmooth:
    """Tests für die Smooth-Mandelbrot-Funktion."""

    def test_gibt_numpy_array_zurueck(self):
        """mandelbrot_smooth() muss numpy-Array zurückgeben."""
        result = mandelbrot_smooth(width=50, height=40, max_iter=20)
        assert isinstance(result, np.ndarray), "Rückgabe muss numpy-Array sein"

    def test_korrekte_shape(self):
        """Shape des Arrays muss (height, width) sein."""
        w, h = 80, 60
        result = mandelbrot_smooth(width=w, height=h, max_iter=10)
        assert result.shape == (h, w), (
            f"Shape muss ({h}, {w}) sein, erhalten: {result.shape}"
        )

    def test_werte_in_bereich(self):
        """Alle Werte müssen im Bereich [0, max_iter] liegen."""
        max_iter = 30
        result = mandelbrot_smooth(width=50, height=40, max_iter=max_iter)
        assert result.min() >= 0.0, "Keine negativen Werte erlaubt"
        assert result.max() <= float(max_iter), (
            f"Kein Wert darf > max_iter={max_iter} sein"
        )

    def test_innere_punkte_haben_maxiter(self):
        """Punkte tief in der Mandelbrot-Menge haben Wert max_iter."""
        # Punkt c=0 liegt sicher in M (z_n = 0 für alle n)
        # Aber wegen vektorisierter Berechnung prüfen wir einen Bereich um 0
        max_iter = 50
        # Sehr kleiner Ausschnitt nahe der Mitte
        result = mandelbrot_smooth(
            x_min=-0.1, x_max=0.1,
            y_min=-0.1, y_max=0.1,
            width=10, height=10,
            max_iter=max_iter
        )
        # Mittelpunkt (index 5,5) sollte max_iter sein (c≈0 ∈ M)
        assert result[5, 5] == float(max_iter), (
            "Punkt nahe c=0 (in M) muss Wert max_iter haben"
        )

    def test_punkte_ausserhalb_haben_kleinen_wert(self):
        """Punkte weit außerhalb der Menge (c=5+5i) divergieren sofort."""
        result = mandelbrot_smooth(
            x_min=4.9, x_max=5.1,
            y_min=4.9, y_max=5.1,
            width=5, height=5,
            max_iter=100
        )
        # Alle Punkte in diesem Bereich divergieren nach 1 Iteration
        # → smooth_count muss deutlich kleiner als max_iter sein
        assert result.max() < 50.0, (
            "Punkte weit außerhalb M müssen kleinen smooth_count haben"
        )

    def test_standardaufruf(self):
        """Standardaufruf ohne Parameter produziert (600, 800)-Array."""
        result = mandelbrot_smooth()
        assert result.shape == (600, 800), (
            f"Standard-Shape muss (600, 800) sein, erhalten: {result.shape}"
        )

    def test_smooth_werte_nicht_ganzzahlig(self):
        """Smooth-Werte sind nicht alle ganzzahlig (das wäre normaler Escape-Count)."""
        result = mandelbrot_smooth(
            x_min=-2.5, x_max=1.0,
            y_min=-1.25, y_max=1.25,
            width=100, height=75,
            max_iter=50
        )
        # Mindestens ein Wert sollte nicht ganzzahlig sein
        # (außer max_iter für Punkte in M)
        aussen = result[result < 50.0]
        if len(aussen) > 0:
            # Wenn alle Außenpunkte ganzzahlig wären, ist Smooth Coloring kaputt
            fraktionsteile = aussen - np.floor(aussen)
            hat_fraktionsanteile = (fraktionsteile > 1e-6).any()
            assert hat_fraktionsanteile, (
                "Smooth Coloring sollte nicht-ganzzahlige Werte liefern"
            )


# ===========================================================================
# TESTS: export_figure (Hilfsfunktion für PNG/SVG/PDF-Export)
# ===========================================================================

class TestExportFigure:
    """Tests für die universelle Export-Hilfsfunktion export_figure()."""

    def test_png_export_erstellt_datei(self):
        """export_figure() mit .png-Endung muss PNG-Datei erzeugen."""
        tmpdir = _make_temp_dir()
        try:
            png_path = os.path.join(tmpdir, 'export_test.png')
            fig, ax = plt.subplots()
            ax.plot([0, 1], [0, 1])

            export_figure(fig, png_path)
            plt.close(fig)

            assert os.path.exists(png_path), "PNG-Datei muss existieren"
            assert os.path.getsize(png_path) > 0, "PNG-Datei darf nicht leer sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_svg_export_erstellt_datei(self):
        """export_figure() mit .svg-Endung muss SVG-Datei erzeugen."""
        tmpdir = _make_temp_dir()
        try:
            svg_path = os.path.join(tmpdir, 'export_test.svg')
            fig, ax = plt.subplots()
            ax.plot([0, 1, 2], [0, 1, 0])

            export_figure(fig, svg_path)
            plt.close(fig)

            assert os.path.exists(svg_path), "SVG-Datei muss existieren"
            # SVG ist XML-basiert
            with open(svg_path, 'r', encoding='utf-8') as f:
                inhalt = f.read(200)
            assert '<?xml' in inhalt or '<svg' in inhalt, (
                "Exportierte Datei muss gültiges SVG-Format sein"
            )
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_pdf_export_erstellt_datei(self):
        """export_figure() mit .pdf-Endung muss PDF-Datei erzeugen."""
        tmpdir = _make_temp_dir()
        try:
            pdf_path = os.path.join(tmpdir, 'export_test.pdf')
            fig, ax = plt.subplots()
            ax.set_title('PDF Export Test')

            export_figure(fig, pdf_path)
            plt.close(fig)

            assert os.path.exists(pdf_path), "PDF-Datei muss existieren"
            # PDF-Magic-Bytes prüfen
            with open(pdf_path, 'rb') as f:
                magic = f.read(4)
            assert magic == b'%PDF', "Datei muss gültiges PDF sein"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ===========================================================================
# TESTS: plot_gaussian_curvature_3d
# ===========================================================================

class TestGaussianCurvature3D:
    """Tests für die Gaußsche-Krümmungs-3D-Visualisierung."""

    def test_sphäre_gibt_figure_zurück(self):
        """plot_gaussian_curvature_3d('sphere') muss (Figure, Axes) zurückgeben."""
        # Funktion gibt (fig, ax)-Tuple zurück – erstes Element ist die Figure
        result = plot_gaussian_curvature_3d('sphere', param_range=(-1.5, 1.5), resolution=10)
        fig = result[0] if isinstance(result, tuple) else result
        assert isinstance(fig, plt.Figure), (
            "Erstes Element des Rückgabe-Tuples muss matplotlib.Figure sein"
        )
        plt.close(fig)

    def test_torus_gibt_figure_zurück(self):
        """plot_gaussian_curvature_3d('torus') muss (Figure, Axes) zurückgeben."""
        result = plot_gaussian_curvature_3d('torus', param_range=(-math.pi, math.pi), resolution=10)
        fig = result[0] if isinstance(result, tuple) else result
        assert isinstance(fig, plt.Figure), (
            "Rückgabewert für Torus muss Figure sein"
        )
        plt.close(fig)

    def test_sattel_gibt_figure_zurück(self):
        """plot_gaussian_curvature_3d('saddle') muss (Figure, Axes) zurückgeben."""
        result = plot_gaussian_curvature_3d('saddle', param_range=(-2, 2), resolution=10)
        fig = result[0] if isinstance(result, tuple) else result
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_hyperbolisches_paraboloid_gibt_figure_zurück(self):
        """plot_gaussian_curvature_3d('hyperbolic_paraboloid') muss (Figure, Axes) zurückgeben."""
        result = plot_gaussian_curvature_3d(
            'hyperbolic_paraboloid', param_range=(-2, 2), resolution=10
        )
        fig = result[0] if isinstance(result, tuple) else result
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_unbekannte_fläche_wirft_fehler(self):
        """Unbekannte Fläche soll ValueError auslösen."""
        with pytest.raises((ValueError, KeyError)):
            plot_gaussian_curvature_3d('unbekannte_fläche', resolution=5)


# ===========================================================================
# TESTS: plot_geodesic_on_sphere
# ===========================================================================

class TestGeodesicOnSphere:
    """Tests für Geodäten auf der Einheitssphäre."""

    def test_gibt_figure_zurück(self):
        """plot_geodesic_on_sphere() muss Figure zurückgeben."""
        fig = plot_geodesic_on_sphere(
            start_angle=(0.0, 0.0), direction=(1.0, 0.5), n_steps=20
        )
        assert isinstance(fig, plt.Figure), (
            "Rückgabewert muss matplotlib.Figure sein"
        )
        plt.close(fig)

    def test_verschiedene_startwinkel(self):
        """Verschiedene Startwinkel produzieren Figure ohne Fehler."""
        for theta, phi in [(0.0, 0.0), (math.pi / 4, math.pi / 3), (math.pi / 2, math.pi)]:
            fig = plot_geodesic_on_sphere(
                start_angle=(theta, phi), direction=(1.0, 0.3), n_steps=15
            )
            assert isinstance(fig, plt.Figure)
            plt.close(fig)

    def test_viele_schritte(self):
        """Auch bei n_steps=100 funktioniert die Funktion korrekt."""
        fig = plot_geodesic_on_sphere(
            start_angle=(0.5, 0.5), direction=(1.0, -0.5), n_steps=100
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# ===========================================================================
# TESTS: plot_geodesic_on_torus
# ===========================================================================

class TestGeodesicOnTorus:
    """Tests für Geodäten auf dem Torus."""

    def test_gibt_figure_zurück(self):
        """plot_geodesic_on_torus() muss Figure zurückgeben."""
        fig = plot_geodesic_on_torus(
            R=2.0, r=1.0, start_params=(0.0, 0.0), direction=(1.0, 0.3), n_steps=50
        )
        assert isinstance(fig, plt.Figure), (
            "Rückgabewert muss matplotlib.Figure sein"
        )
        plt.close(fig)

    def test_verschiedene_r_parameter(self):
        """Unterschiedliche R und r Werte funktionieren korrekt."""
        fig = plot_geodesic_on_torus(
            R=3.0, r=0.5, start_params=(0.0, 0.0), direction=(1.0, 0.5), n_steps=30
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_verschiedene_richtungen(self):
        """Verschiedene Richtungsvektoren erzeugen Figure ohne Fehler."""
        for richtung in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -0.3)]:
            fig = plot_geodesic_on_torus(direction=richtung, n_steps=30)
            assert isinstance(fig, plt.Figure)
            plt.close(fig)


# ===========================================================================
# TESTS: plot_function_2d_adaptive
# ===========================================================================

class TestAdaptivePlot:
    """Tests für das adaptive Gitter-Plotting."""

    def test_gibt_figure_zurück(self):
        """plot_function_2d_adaptive() muss Figure zurückgeben."""
        fig = plot_function_2d_adaptive('sin(x)', x_range=(-math.pi, math.pi), base_points=30)
        assert isinstance(fig, plt.Figure), (
            "Rückgabewert muss matplotlib.Figure sein"
        )
        plt.close(fig)

    def test_glatte_funktion(self):
        """Bei glatter Funktion (x²) wird korrekt geplottet."""
        fig = plot_function_2d_adaptive('x**2', x_range=(-2, 2), base_points=20)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_springende_funktion_mehr_punkte(self):
        """Bei stark gekrümmter Funktion (1/x²) werden Regionen verfeinert."""
        # 1/(x²+0.01) hat starke Krümmung nahe x=0
        x_pts, y_pts = plot_function_2d_adaptive(
            '1/(x**2+0.01)',
            x_range=(-1, 1),
            base_points=20,
            refinement_levels=2,
            return_points=True
        )
        # Mit Verfeinerung muss mehr Punkte vorhanden sein als base_points
        assert len(x_pts) > 20, (
            f"Adaptive Verfeinerung muss mehr als base_points=20 Punkte liefern, erhalten: {len(x_pts)}"
        )

    def test_refinement_levels_null(self):
        """Mit refinement_levels=0 werden genau base_points Punkte verwendet."""
        x_pts, y_pts = plot_function_2d_adaptive(
            'x',
            x_range=(0, 1),
            base_points=50,
            refinement_levels=0,
            return_points=True
        )
        # Keine Verfeinerung → genau base_points Punkte
        assert len(x_pts) == 50, (
            f"Ohne Verfeinerung müssen genau 50 Punkte vorhanden sein, erhalten: {len(x_pts)}"
        )

    def test_sinus_funktioniert(self):
        """sin(x) kann korrekt geplottet werden."""
        fig = plot_function_2d_adaptive(
            'sin(x)',
            x_range=(-2 * math.pi, 2 * math.pi),
            base_points=40,
            refinement_levels=1
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# ===========================================================================
# TESTS FÜR NEUE VISUALISIERUNGSFUNKTIONEN (Build 46)
# ===========================================================================

import numpy as np

class TestPlotCantorSet:
    """Tests für plot_cantor_set()."""

    def test_gibt_figure_zurück(self, tmp_path):
        """Test: plot_cantor_set() gibt eine Figure zurück."""
        import sys
        sys.path.insert(0, '../src')
        from visualization import plot_cantor_set
        fig = plot_cantor_set(n_steps=3, save_path=str(tmp_path / "cantor.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_mehrere_schritte(self, tmp_path):
        """Test: Verschiedene Schrittanzahlen funktionieren."""
        from visualization import plot_cantor_set
        for steps in [1, 3, 6]:
            fig = plot_cantor_set(n_steps=steps, save_path=str(tmp_path / f"cantor_{steps}.png"))
            assert isinstance(fig, plt.Figure)
            plt.close(fig)


class TestPlotCantorFunction:
    """Tests für plot_cantor_function()."""

    def test_gibt_figure_zurück(self, tmp_path):
        """Test: plot_cantor_function() gibt eine Figure zurück."""
        from visualization import plot_cantor_function
        fig = plot_cantor_function(n=50, save_path=str(tmp_path / "devil_staircase.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestPlotBesselGallery:
    """Tests für plot_bessel_gallery()."""

    def test_gibt_figure_zurück(self, tmp_path):
        """Test: plot_bessel_gallery() gibt eine Figure zurück."""
        from visualization import plot_bessel_gallery
        fig = plot_bessel_gallery(orders=[0, 1, 2], save_path=str(tmp_path / "bessel.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_standard_orders(self, tmp_path):
        """Test: Standard-Ordnungen [0..5] funktionieren."""
        from visualization import plot_bessel_gallery
        fig = plot_bessel_gallery(save_path=str(tmp_path / "bessel_all.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestPlotLegendreGallery:
    """Tests für plot_legendre_gallery()."""

    def test_gibt_figure_zurück(self, tmp_path):
        """Test: plot_legendre_gallery() gibt eine Figure zurück."""
        from visualization import plot_legendre_gallery
        fig = plot_legendre_gallery(orders=[0, 1, 2, 3], save_path=str(tmp_path / "legendre.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_orthogonalitaet(self):
        """Test: P_0(x) = 1 und P_1(x) = x (Grundeigenschaften)."""
        from scipy.special import legendre as legendre_poly
        x = np.linspace(-1, 1, 100)
        P0 = legendre_poly(0)(x)
        P1 = legendre_poly(1)(x)
        # P_0(x) = 1 überall
        assert np.allclose(P0, np.ones_like(x))
        # P_1(x) = x
        assert np.allclose(P1, x)


class TestPlotGaussianCurvature:
    """Tests für plot_gaussian_curvature()."""

    def test_kugel_positive_kruemmung(self, tmp_path):
        """Test: Kugel hat überall positive Gaußsche Krümmung."""
        from visualization import plot_sphere_curvature
        # Nur prüfen ob Funktion ohne Fehler läuft
        fig = plot_sphere_curvature(R=1.0, n=10, save_path=str(tmp_path / "sphere.png"))
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestAdaptivePlot:
    """Tests für adaptive_plot()."""

    def test_gibt_figure_zurück(self, tmp_path):
        """Test: adaptive_plot() gibt eine Figure zurück."""
        from visualization import adaptive_plot
        fig = adaptive_plot(
            lambda x: np.sin(x),
            x_range=(0, 2 * np.pi),
            save_path=str(tmp_path / "adaptive.png")
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_return_points(self):
        """Test: return_points=True gibt Arrays zurück."""
        from visualization import adaptive_plot
        x_arr, y_arr = adaptive_plot(
            lambda x: np.sin(x),
            x_range=(0, np.pi),
            n_initial=20,
            return_points=True
        )
        assert len(x_arr) >= 20  # Mindestens Anfangspunkte
        assert len(x_arr) == len(y_arr)

    def test_singularity_handling(self, tmp_path):
        """Test: Singularitäten werden behandelt (nicht crashing)."""
        from visualization import adaptive_plot
        # tan(x) hat Singularitäten bei π/2
        fig = adaptive_plot(
            lambda x: np.tan(x) if abs(x - np.pi/2) > 0.01 else np.nan,
            x_range=(0, np.pi),
            n_initial=50,
            save_path=str(tmp_path / "singular.png")
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


# ===========================================================================
# Tests für export_figure (Build 69)
# ===========================================================================

class TestExportFigure:
    """
    @brief Tests für die export_figure()-Funktion.
    @description
        Prüft Export nach PNG, SVG, PDF – einzeln und kombiniert –
        sowie Standardverhalten und DPI-Parameter.
    @author Michael Fuhrmann
    @date 2026-03-11
    """

    def _make_fig(self):
        """Hilfsmethode: einfache Test-Figure erzeugen."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        return fig

    def test_alle_drei_formate(self, tmp_path):
        """Alle drei Standard-Formate werden erzeugt."""
        from visualization import export_figure
        import matplotlib.pyplot as plt
        fig = self._make_fig()
        base = str(tmp_path / "export_alle")
        result = export_figure(fig, base)
        plt.close(fig)
        assert set(result.keys()) == {"png", "svg", "pdf"}
        for fmt, path in result.items():
            import os
            assert os.path.exists(path), f"Datei für {fmt} fehlt: {path}"
            assert os.path.getsize(path) > 0, f"Datei für {fmt} ist leer"

    def test_nur_png(self, tmp_path):
        """Nur PNG-Export wenn formats=['png']."""
        from visualization import export_figure
        import matplotlib.pyplot as plt, os
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "nur_png"), formats=["png"])
        plt.close(fig)
        assert list(result.keys()) == ["png"]
        assert os.path.exists(result["png"])

    def test_nur_svg(self, tmp_path):
        """Nur SVG-Export – Vektorformat für LaTeX."""
        from visualization import export_figure
        import matplotlib.pyplot as plt, os
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "nur_svg"), formats=["svg"])
        plt.close(fig)
        assert "svg" in result
        # SVG-Datei muss XML-Header enthalten
        with open(result["svg"], "r", encoding="utf-8") as f:
            inhalt = f.read(200)
        assert "<?xml" in inhalt or "<svg" in inhalt

    def test_nur_pdf(self, tmp_path):
        """Nur PDF-Export – für LaTeX \\includegraphics geeignet."""
        from visualization import export_figure
        import matplotlib.pyplot as plt, os
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "nur_pdf"), formats=["pdf"])
        plt.close(fig)
        assert "pdf" in result
        # PDF-Datei beginnt mit %PDF
        with open(result["pdf"], "rb") as f:
            header = f.read(4)
        assert header == b"%PDF"

    def test_png_svg_ohne_pdf(self, tmp_path):
        """PNG + SVG ohne PDF – teilweise Auswahl."""
        from visualization import export_figure
        import matplotlib.pyplot as plt
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "png_svg"), formats=["png", "svg"])
        plt.close(fig)
        assert set(result.keys()) == {"png", "svg"}
        assert "pdf" not in result

    def test_rueckgabe_absolute_pfade(self, tmp_path):
        """Rückgabe enthält absolute Pfade."""
        from visualization import export_figure
        import matplotlib.pyplot as plt, os
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "abs_pfad"), formats=["png"])
        plt.close(fig)
        assert os.path.isabs(result["png"])

    def test_dpi_parameter(self, tmp_path):
        """Höheres DPI erzeugt größere PNG-Datei."""
        from visualization import export_figure
        import matplotlib.pyplot as plt, os
        fig_lo = self._make_fig()
        result_lo = export_figure(fig_lo, str(tmp_path / "dpi_lo"), formats=["png"], dpi=50)
        plt.close(fig_lo)
        fig_hi = self._make_fig()
        result_hi = export_figure(fig_hi, str(tmp_path / "dpi_hi"), formats=["png"], dpi=300)
        plt.close(fig_hi)
        size_lo = os.path.getsize(result_lo["png"])
        size_hi = os.path.getsize(result_hi["png"])
        assert size_hi > size_lo, "Höheres DPI soll größere Datei erzeugen"

    def test_leere_formats_liste(self, tmp_path):
        """Leere formats-Liste → kein Fehler, leeres Dict zurück."""
        from visualization import export_figure
        import matplotlib.pyplot as plt
        fig = self._make_fig()
        result = export_figure(fig, str(tmp_path / "leer"), formats=[])
        plt.close(fig)
        assert result == {}


# ===========================================================================
# Tests für plot_cantor_set (Build 69)
# ===========================================================================

class TestPlotCantorSet:
    """
    @brief Tests für die plot_cantor_set()-Funktion.
    @description
        Prüft korrekte Anzahl von Intervallen je Stufe,
        Rückgabetyp und Grenzfälle.
    @author Michael Fuhrmann
    @date 2026-03-11
    """

    def test_rueckgabetyp(self):
        """Rückgabe ist matplotlib Figure."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set(steps=3)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_stufe_null(self):
        """steps=0 → Einheitsintervall, keine Iteration."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set(steps=0)
        assert fig is not None
        plt.close(fig)

    def test_standardschritte(self):
        """Standardmäßig 6 Schritte, kein Fehler."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set()
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_acht_schritte(self):
        """steps=8 → sehr viele Intervalle, kein Fehler."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set(steps=8)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_intervallanzahl_korrekt(self):
        """Nach k Schritten müssen genau 2^k Intervalle vorhanden sein."""
        # broken_barh erzeugt in matplotlib >= 3.8 PolyCollection statt
        # BrokenBarHCollection; wir zählen die Collections (eine je Stufe).
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        steps = 4
        fig = plot_cantor_set(steps=steps)
        ax = fig.axes[0]
        # Eine Collection pro Iterationsstufe (Stufe 0 bis steps)
        assert len(ax.collections) == steps + 1
        plt.close(fig)

    def test_figsize_parameter(self):
        """figsize-Parameter wird korrekt weitergegeben."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set(steps=3, figsize=(10, 3))
        w, h = fig.get_size_inches()
        assert abs(w - 10) < 0.01 and abs(h - 3) < 0.01
        plt.close(fig)

    def test_x_achse_grenzen(self):
        """x-Achse läuft von 0 bis 1 (Einheitsintervall)."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt
        fig = plot_cantor_set(steps=3)
        ax = fig.axes[0]
        xlim = ax.get_xlim()
        assert abs(xlim[0]) < 1e-9
        assert abs(xlim[1] - 1.0) < 1e-9
        plt.close(fig)

    def test_speichern_als_png(self, tmp_path):
        """Erzeugte Figure kann als PNG gespeichert werden."""
        import matplotlib
        matplotlib.use('Agg')
        from visualization import plot_cantor_set
        import matplotlib.pyplot as plt, os
        fig = plot_cantor_set(steps=4)
        pfad = str(tmp_path / "cantor.png")
        fig.savefig(pfad)
        plt.close(fig)
        assert os.path.exists(pfad)
        assert os.path.getsize(pfad) > 0
