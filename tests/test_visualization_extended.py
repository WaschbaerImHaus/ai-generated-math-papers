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
import math
import tempfile
import shutil

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
