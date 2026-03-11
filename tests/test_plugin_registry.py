"""
@file test_plugin_registry.py
@brief Tests für das Plugin-Registrierungs-System.
@description
    Testet alle Methoden der PluginRegistry-Klasse:
    - register(): Plugin registrieren
    - load(): Plugin laden
    - list_plugins(): Alle Plugins auflisten
    - is_loaded(): Prüfen ob Plugin geladen
    - unregister(): Plugin entfernen

@author Kurt Ingwer
@date 2026-03-11
@lastModified 2026-03-11
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
from plugin_registry import PluginRegistry, DEFAULT_REGISTRY


class TestPluginRegistryRegister:
    """Tests für register() und list_plugins()."""

    def setup_method(self):
        """Erstellt eine frische Registry vor jedem Test."""
        self.registry = PluginRegistry()

    def test_register_single_plugin(self):
        """Test: Ein Plugin kann registriert werden."""
        self.registry.register("test_plugin", "os.path")
        assert "test_plugin" in self.registry.list_plugins()

    def test_register_multiple_plugins(self):
        """Test: Mehrere Plugins können registriert werden."""
        self.registry.register("plugin_a", "os.path")
        self.registry.register("plugin_b", "sys")
        plugins = self.registry.list_plugins()
        assert "plugin_a" in plugins
        assert "plugin_b" in plugins

    def test_register_overwrites_existing(self):
        """Test: Erneutes Registrieren überschreibt den Modulpfad."""
        self.registry.register("my_plugin", "os.path")
        self.registry.register("my_plugin", "sys")
        # Plugin noch vorhanden (ein Eintrag)
        assert self.registry.list_plugins().count("my_plugin") == 1

    def test_register_empty_name_raises(self):
        """Test: Leerer Name wirft ValueError."""
        with pytest.raises(ValueError, match="Name"):
            self.registry.register("", "os.path")

    def test_register_empty_path_raises(self):
        """Test: Leerer Modulpfad wirft ValueError."""
        with pytest.raises(ValueError, match="Modulpfad"):
            self.registry.register("plugin", "")

    def test_list_plugins_alphabetical(self):
        """Test: list_plugins() gibt alphabetisch sortierte Liste zurück."""
        self.registry.register("z_plugin", "sys")
        self.registry.register("a_plugin", "os")
        self.registry.register("m_plugin", "math")
        plugins = self.registry.list_plugins()
        assert plugins == sorted(plugins)

    def test_list_plugins_empty(self):
        """Test: Leere Registry gibt leere Liste zurück."""
        assert self.registry.list_plugins() == []


class TestPluginRegistryLoad:
    """Tests für load() und is_loaded()."""

    def setup_method(self):
        """Erstellt eine frische Registry vor jedem Test."""
        self.registry = PluginRegistry()

    def test_load_standard_module(self):
        """Test: Standardbibliotheks-Modul kann geladen werden."""
        self.registry.register("math_module", "math")
        modul = self.registry.load("math_module")
        assert modul is not None
        # Prüfen ob es wirklich das math-Modul ist
        assert hasattr(modul, 'pi')
        assert hasattr(modul, 'sqrt')

    def test_load_caches_module(self):
        """Test: Zweimaliges Laden gibt dasselbe Objekt zurück (Cache)."""
        self.registry.register("math_module", "math")
        first_load = self.registry.load("math_module")
        second_load = self.registry.load("math_module")
        assert first_load is second_load

    def test_load_unregistered_raises_key_error(self):
        """Test: Laden eines nicht registrierten Plugins wirft KeyError."""
        with pytest.raises(KeyError, match="nicht registriert"):
            self.registry.load("nicht_vorhanden")

    def test_load_nonexistent_module_raises_import_error(self):
        """Test: Nicht-existierendes Modul wirft ImportError."""
        self.registry.register("fake_module", "dieser.modul.existiert.nicht")
        with pytest.raises(ImportError):
            self.registry.load("fake_module")

    def test_is_loaded_before_load(self):
        """Test: is_loaded() gibt False zurück vor dem ersten Laden."""
        self.registry.register("math_module", "math")
        assert self.registry.is_loaded("math_module") is False

    def test_is_loaded_after_load(self):
        """Test: is_loaded() gibt True zurück nach erfolgreichem Laden."""
        self.registry.register("math_module", "math")
        self.registry.load("math_module")
        assert self.registry.is_loaded("math_module") is True

    def test_is_loaded_unregistered_plugin(self):
        """Test: is_loaded() gibt False für unregistrierte Plugins zurück."""
        assert self.registry.is_loaded("nicht_registriert") is False


class TestPluginRegistryUnregister:
    """Tests für unregister()."""

    def setup_method(self):
        """Erstellt eine frische Registry vor jedem Test."""
        self.registry = PluginRegistry()

    def test_unregister_removes_plugin(self):
        """Test: unregister() entfernt Plugin aus der Liste."""
        self.registry.register("to_remove", "math")
        self.registry.unregister("to_remove")
        assert "to_remove" not in self.registry.list_plugins()

    def test_unregister_clears_cache(self):
        """Test: unregister() löscht auch den Modul-Cache."""
        self.registry.register("math_module", "math")
        self.registry.load("math_module")  # In Cache laden
        assert self.registry.is_loaded("math_module") is True

        self.registry.unregister("math_module")
        assert self.registry.is_loaded("math_module") is False

    def test_unregister_nonexistent_no_error(self):
        """Test: unregister() für nicht-registriertes Plugin wirft keinen Fehler."""
        # Soll idempotent sein
        self.registry.unregister("existiert_nicht")  # Kein Fehler erwartet

    def test_unregister_prevents_load(self):
        """Test: Nach unregister() kann das Plugin nicht mehr geladen werden."""
        self.registry.register("to_remove", "math")
        self.registry.unregister("to_remove")
        with pytest.raises(KeyError):
            self.registry.load("to_remove")


class TestDefaultRegistry:
    """Tests für die globale DEFAULT_REGISTRY."""

    def test_default_registry_exists(self):
        """Test: DEFAULT_REGISTRY ist vorhanden und ist PluginRegistry."""
        assert isinstance(DEFAULT_REGISTRY, PluginRegistry)

    def test_default_registry_has_algebra(self):
        """Test: Standard-Plugins (algebra_core) sind registriert."""
        assert "algebra_core" in DEFAULT_REGISTRY.list_plugins()

    def test_default_registry_has_fourier(self):
        """Test: fourier ist in der Standard-Registry."""
        assert "fourier" in DEFAULT_REGISTRY.list_plugins()

    def test_default_registry_has_many_plugins(self):
        """Test: Standard-Registry hat mindestens 10 vorregistrierte Plugins."""
        assert len(DEFAULT_REGISTRY.list_plugins()) >= 10

    def test_default_registry_load_algebra(self):
        """Test: algebra_core-Plugin kann aus der Standard-Registry geladen werden."""
        modul = DEFAULT_REGISTRY.load("algebra_core")
        assert modul is not None
        # algebra-Modul hat is_prime-Funktion
        assert hasattr(modul, 'is_prime')

    def test_default_registry_load_fourier(self):
        """Test: fourier-Plugin kann aus der Standard-Registry geladen werden."""
        modul = DEFAULT_REGISTRY.load("fourier")
        assert modul is not None
        # fourier-Modul hat fft-Funktion
        assert hasattr(modul, 'fft')


class TestPluginRegistryRepr:
    """Tests für die __repr__()-Methode."""

    def test_repr_shows_count(self):
        """Test: __repr__() zeigt Anzahl registrierter Plugins."""
        registry = PluginRegistry()
        registry.register("a", "math")
        registry.register("b", "sys")
        repr_str = repr(registry)
        assert "2" in repr_str  # 2 registrierte Plugins
        assert "PluginRegistry" in repr_str
