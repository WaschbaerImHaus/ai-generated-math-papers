"""
@file test_plugin_system.py
@brief Tests für das Plugin-System (src/plugin_system.py).
@description
    Testet alle Funktionen des Plugin-Systems:
    - Registrierung und Abruf von Plugins
    - Laden aus Dateien
    - Automatische Erkennung im Verzeichnis
    - Fehlerbehandlung bei ungültigen Plugins
    - Template-Erstellung
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import os
import sys
import tempfile
import types

import pytest

# Projekt-Pfad einfügen, damit Module gefunden werden
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

# Plugin-Registry vor jedem Test leeren (isoliert Tests voneinander)
import plugin_system


@pytest.fixture(autouse=True)
def clear_registry():
    """Leert die Plugin-Registry vor und nach jedem Test."""
    plugin_system.PLUGIN_REGISTRY.clear()
    yield
    plugin_system.PLUGIN_REGISTRY.clear()


# ---------------------------------------------------------------------------
# Tests: register_plugin und get_plugin
# ---------------------------------------------------------------------------

class TestRegisterAndGet:
    """Tests für register_plugin() und get_plugin()."""

    def test_register_simple_module(self):
        """Registrierung eines einfachen Modul-Objekts."""
        # Einfaches Modul-Objekt erstellen
        m = types.ModuleType("testmod")
        m.value = 42

        plugin_system.register_plugin("testmod", m)
        retrieved = plugin_system.get_plugin("testmod")

        assert retrieved is m
        assert retrieved.value == 42

    def test_register_overwrites_existing(self):
        """Überschreiben eines vorhandenen Plugins."""
        m1 = types.ModuleType("mod")
        m1.version = 1
        m2 = types.ModuleType("mod")
        m2.version = 2

        plugin_system.register_plugin("myplugin", m1)
        plugin_system.register_plugin("myplugin", m2)  # Überschreiben

        assert plugin_system.get_plugin("myplugin").version == 2

    def test_get_nonexistent_raises_key_error(self):
        """get_plugin() für nicht-existentes Plugin wirft KeyError."""
        with pytest.raises(KeyError, match="not_there"):
            plugin_system.get_plugin("not_there")

    def test_register_multiple_plugins(self):
        """Mehrere Plugins gleichzeitig registrieren."""
        for i in range(5):
            m = types.ModuleType(f"plugin_{i}")
            plugin_system.register_plugin(f"plugin_{i}", m)

        assert len(plugin_system.PLUGIN_REGISTRY) == 5


# ---------------------------------------------------------------------------
# Tests: list_plugins und unload_plugin
# ---------------------------------------------------------------------------

class TestListAndUnload:
    """Tests für list_plugins() und unload_plugin()."""

    def test_list_plugins_empty(self):
        """Leere Registry → leere Liste."""
        assert plugin_system.list_plugins() == []

    def test_list_plugins_sorted(self):
        """list_plugins() gibt sortierte Namen zurück."""
        for name in ["zebra", "alpha", "mango"]:
            plugin_system.register_plugin(name, types.ModuleType(name))

        assert plugin_system.list_plugins() == ["alpha", "mango", "zebra"]

    def test_unload_removes_plugin(self):
        """unload_plugin() entfernt Plugin aus Registry."""
        m = types.ModuleType("temp")
        plugin_system.register_plugin("temp", m)
        assert "temp" in plugin_system.list_plugins()

        plugin_system.unload_plugin("temp")
        assert "temp" not in plugin_system.list_plugins()

    def test_unload_nonexistent_raises_key_error(self):
        """unload_plugin() für nicht-existentes Plugin wirft KeyError."""
        with pytest.raises(KeyError):
            plugin_system.unload_plugin("does_not_exist")


# ---------------------------------------------------------------------------
# Tests: load_plugin
# ---------------------------------------------------------------------------

class TestLoadPlugin:
    """Tests für load_plugin()."""

    def test_load_valid_plugin_file(self, tmp_path):
        """Laden einer gültigen Plugin-Datei."""
        plugin_file = tmp_path / "my_plugin.py"
        plugin_file.write_text(
            'PLUGIN_NAME = "my_plugin"\n'
            'def add(a, b): return a + b\n',
            encoding="utf-8"
        )

        module = plugin_system.load_plugin(str(plugin_file))
        assert module.PLUGIN_NAME == "my_plugin"
        assert module.add(3, 4) == 7

    def test_load_nonexistent_file_raises(self):
        """Laden einer nicht-existenten Datei wirft FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            plugin_system.load_plugin("/tmp/nonexistent_plugin_xyz.py")

    def test_load_plugin_executes_code(self, tmp_path):
        """Beim Laden wird der Plugin-Code ausgeführt."""
        plugin_file = tmp_path / "math_plugin.py"
        plugin_file.write_text(
            'import math\n'
            'PI_APPROX = round(math.pi, 4)\n'
            'def circle_area(r): return math.pi * r * r\n',
            encoding="utf-8"
        )

        module = plugin_system.load_plugin(str(plugin_file))
        assert module.PI_APPROX == 3.1416
        assert abs(module.circle_area(1) - 3.14159) < 0.001


# ---------------------------------------------------------------------------
# Tests: discover_plugins
# ---------------------------------------------------------------------------

class TestDiscoverPlugins:
    """Tests für discover_plugins()."""

    def test_discover_empty_directory(self, tmp_path):
        """Leeres Verzeichnis → leere Ergebnisliste."""
        result = plugin_system.discover_plugins(str(tmp_path))
        assert result == []

    def test_discover_nonexistent_directory(self):
        """Nicht-existentes Verzeichnis → leere Ergebnisliste."""
        result = plugin_system.discover_plugins("/tmp/nonexistent_plugin_dir_xyz")
        assert result == []

    def test_discover_loads_all_py_files(self, tmp_path):
        """Alle .py-Dateien im Verzeichnis werden geladen."""
        for i in range(3):
            (tmp_path / f"plugin_{i}.py").write_text(
                f'PLUGIN_NAME = "plugin_{i}"\n',
                encoding="utf-8"
            )

        result = plugin_system.discover_plugins(str(tmp_path))
        assert len(result) == 3
        assert sorted(result) == ["plugin_0", "plugin_1", "plugin_2"]

    def test_discover_skips_init_file(self, tmp_path):
        """__init__.py wird nicht als Plugin geladen."""
        (tmp_path / "__init__.py").write_text("# init\n", encoding="utf-8")
        (tmp_path / "real_plugin.py").write_text(
            'PLUGIN_NAME = "real"\n',
            encoding="utf-8"
        )

        result = plugin_system.discover_plugins(str(tmp_path))
        assert result == ["real"]

    def test_discover_uses_plugin_name_attribute(self, tmp_path):
        """PLUGIN_NAME-Attribut wird als Registrierungsname verwendet."""
        (tmp_path / "myfile.py").write_text(
            'PLUGIN_NAME = "fancy_name"\n',
            encoding="utf-8"
        )

        result = plugin_system.discover_plugins(str(tmp_path))
        assert "fancy_name" in result
        assert plugin_system.get_plugin("fancy_name") is not None

    def test_discover_falls_back_to_filename(self, tmp_path):
        """Ohne PLUGIN_NAME wird der Dateiname als Plugin-Name verwendet."""
        (tmp_path / "unnamed_plugin.py").write_text(
            'VALUE = 99\n',
            encoding="utf-8"
        )

        result = plugin_system.discover_plugins(str(tmp_path))
        assert "unnamed_plugin" in result

    def test_discover_skips_broken_plugins(self, tmp_path):
        """Syntaktisch fehlerhafte Plugins werden übersprungen."""
        (tmp_path / "broken.py").write_text(
            'this is not valid python !!!\n',
            encoding="utf-8"
        )
        (tmp_path / "good.py").write_text(
            'PLUGIN_NAME = "good"\n',
            encoding="utf-8"
        )

        result = plugin_system.discover_plugins(str(tmp_path))
        assert "good" in result
        assert "broken" not in result


# ---------------------------------------------------------------------------
# Tests: create_plugin_template
# ---------------------------------------------------------------------------

class TestCreatePluginTemplate:
    """Tests für create_plugin_template()."""

    def test_creates_file(self, tmp_path):
        """Template-Datei wird erstellt."""
        path = plugin_system.create_plugin_template("my_math", str(tmp_path))
        assert os.path.isfile(path)

    def test_template_is_valid_python(self, tmp_path):
        """Erzeugte Template-Datei ist gültiges Python."""
        path = plugin_system.create_plugin_template("test_template", str(tmp_path))
        module = plugin_system.load_plugin(path)
        assert module is not None

    def test_template_has_plugin_name(self, tmp_path):
        """Template enthält PLUGIN_NAME-Variable."""
        path = plugin_system.create_plugin_template("algebra_ext", str(tmp_path))
        module = plugin_system.load_plugin(path)
        assert hasattr(module, "PLUGIN_NAME")
        assert module.PLUGIN_NAME == "algebra_ext"

    def test_template_has_setup_function(self, tmp_path):
        """Template enthält setup()-Funktion."""
        path = plugin_system.create_plugin_template("my_ext", str(tmp_path))
        module = plugin_system.load_plugin(path)
        assert hasattr(module, "setup")
        assert callable(module.setup)

    def test_creates_output_directory(self, tmp_path):
        """Ausgabeverzeichnis wird ggf. erstellt."""
        new_dir = str(tmp_path / "subdir" / "plugins")
        plugin_system.create_plugin_template("nested", new_dir)
        assert os.path.isdir(new_dir)


# ---------------------------------------------------------------------------
# Tests: Beispiel-Plugin (example_plugin.py)
# ---------------------------------------------------------------------------

class TestExamplePlugin:
    """Tests für das mitgelieferte Beispiel-Plugin."""

    @pytest.fixture
    def example_mod(self):
        """Lädt das Beispiel-Plugin direkt."""
        plugins_dir = os.path.join(
            os.path.dirname(__file__), '..', 'plugins', 'example_plugin.py'
        )
        return plugin_system.load_plugin(plugins_dir)

    def test_perfect_number_6(self, example_mod):
        """6 ist eine vollkommene Zahl."""
        assert example_mod.perfect_number_check(6) is True

    def test_perfect_number_28(self, example_mod):
        """28 ist eine vollkommene Zahl."""
        assert example_mod.perfect_number_check(28) is True

    def test_not_perfect_number(self, example_mod):
        """12 ist keine vollkommene Zahl."""
        assert example_mod.perfect_number_check(12) is False

    def test_digit_sum(self, example_mod):
        """Quersumme von 1234 = 10."""
        assert example_mod.digit_sum(1234) == 10

    def test_collatz_length_1(self, example_mod):
        """Collatz-Länge für 1 ist 0 Schritte."""
        assert example_mod.collatz_length(1) == 0

    def test_abundant_numbers(self, example_mod):
        """12 ist die erste abundante Zahl."""
        abundant = example_mod.abundant_numbers(20)
        assert 12 in abundant
        assert 18 in abundant
        assert 10 not in abundant
