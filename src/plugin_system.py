"""
@file plugin_system.py
@brief Plugin-System für mathematische Erweiterungsmodule.
@description
    Dieses Modul ermöglicht das dynamische Laden, Registrieren und Verwalten
    von Plugin-Modulen. Plugins sind eigenständige Python-Dateien, die
    mathematische Funktionalität bereitstellen und zur Laufzeit geladen
    werden können, ohne den Kern-Code zu verändern.

    ## Konzept

    Jedes Plugin ist eine Python-Datei mit einer optionalen
    ``PLUGIN_NAME``-Variable und einer ``setup()``-Funktion. Das Plugin-System
    sucht automatisch im angegebenen Verzeichnis nach solchen Dateien und
    registriert sie in einer zentralen Registry.

    ## Verwendung

    ```python
    from plugin_system import discover_plugins, get_plugin

    # Alle Plugins im Verzeichnis "plugins/" laden
    discover_plugins("plugins/")

    # Einzelnes Plugin abrufen
    plugin = get_plugin("example_plugin")
    result = plugin.my_function(42)
    ```

@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0.0
"""

import importlib
import importlib.util
import os
from typing import Any, Dict, List, Optional


# Globale Plugin-Registry: Name → Modul-Objekt
PLUGIN_REGISTRY: Dict[str, Any] = {}


def register_plugin(name: str, module: Any) -> None:
    """
    @brief Registriert ein Modul-Objekt als Plugin unter einem Namen.
    @description
        Speichert das übergebene Modul-Objekt in der globalen PLUGIN_REGISTRY
        unter dem angegebenen Namen. Ein vorhandenes Plugin gleichen Namens
        wird überschrieben.

    @param name   Bezeichner, unter dem das Plugin registriert wird.
    @param module Das geladene Python-Modul-Objekt.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Plugin in der globalen Registry unter dem Namen speichern
    PLUGIN_REGISTRY[name] = module


def load_plugin(path: str) -> Any:
    """
    @brief Lädt ein Plugin aus einer Datei und gibt das Modul-Objekt zurück.
    @description
        Verwendet ``importlib.util.spec_from_file_location``, um eine beliebige
        Python-Datei als Modul zu laden, ohne sie im Python-Suchpfad haben zu
        müssen. Das Modul wird ausgeführt (``exec_module``) und zurückgegeben.

    @param path  Absoluter oder relativer Pfad zur Plugin-Datei (.py).
    @return      Das geladene Modul-Objekt.
    @raises FileNotFoundError  Wenn die Datei nicht existiert.
    @raises ImportError        Wenn das Modul nicht geladen werden kann.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Sicherstellen, dass die Datei existiert
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Plugin-Datei nicht gefunden: {path}")

    # Modulnamen aus dem Dateinamen ableiten (ohne .py-Endung)
    module_name = os.path.splitext(os.path.basename(path))[0]

    # Modul-Spezifikation aus dem Dateipfad erstellen
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Konnte Modul-Spec nicht erstellen für: {path}")

    # Leeres Modul-Objekt erzeugen
    module = importlib.util.module_from_spec(spec)

    # Modul-Code ausführen (Definitionen und globale Anweisungen laden)
    spec.loader.exec_module(module)  # type: ignore[union-attr]

    return module


def discover_plugins(plugin_dir: str = "plugins/") -> List[str]:
    """
    @brief Sucht im angegebenen Verzeichnis nach Python-Dateien und lädt sie.
    @description
        Durchsucht ``plugin_dir`` nach ``.py``-Dateien (außer ``__init__.py``),
        lädt jede gefundene Datei als Plugin und registriert sie in der
        PLUGIN_REGISTRY. Gibt eine Liste der Namen aller erfolgreich geladenen
        Plugins zurück.

    @param plugin_dir  Pfad zum Plugin-Verzeichnis (Standard: "plugins/").
    @return            Liste der Namen der geladenen Plugins.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    loaded: List[str] = []

    # Verzeichnis existenz prüfen
    if not os.path.isdir(plugin_dir):
        return loaded

    # Alle Dateien im Verzeichnis durchgehen
    for filename in sorted(os.listdir(plugin_dir)):
        # Nur Python-Dateien laden, __init__.py ausschließen
        if not filename.endswith(".py") or filename == "__init__.py":
            continue

        file_path = os.path.join(plugin_dir, filename)
        module_name = os.path.splitext(filename)[0]

        try:
            # Plugin laden
            module = load_plugin(file_path)

            # Registrierungsname: PLUGIN_NAME-Attribut nutzen oder Dateiname
            plugin_name = getattr(module, "PLUGIN_NAME", module_name)

            # In der Registry speichern
            register_plugin(plugin_name, module)
            loaded.append(plugin_name)
        except Exception as exc:
            # Fehlerhafte Plugins überspringen, aber melden
            print(f"Warnung: Plugin '{filename}' konnte nicht geladen werden: {exc}")

    return loaded


def get_plugin(name: str) -> Any:
    """
    @brief Gibt ein registriertes Plugin anhand seines Namens zurück.
    @description
        Sucht in der globalen PLUGIN_REGISTRY nach dem Plugin mit dem
        angegebenen Namen und gibt das Modul-Objekt zurück.

    @param name  Name des registrierten Plugins.
    @return      Das Modul-Objekt des Plugins.
    @raises KeyError  Wenn kein Plugin mit diesem Namen registriert ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    if name not in PLUGIN_REGISTRY:
        raise KeyError(f"Plugin '{name}' ist nicht registriert. "
                       f"Verfügbare Plugins: {list(PLUGIN_REGISTRY.keys())}")
    return PLUGIN_REGISTRY[name]


def list_plugins() -> List[str]:
    """
    @brief Gibt eine Liste der Namen aller registrierten Plugins zurück.
    @return  Sortierte Liste der Plugin-Namen.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    return sorted(PLUGIN_REGISTRY.keys())


def unload_plugin(name: str) -> None:
    """
    @brief Entfernt ein Plugin aus der Registry.
    @description
        Löscht den Eintrag des Plugins aus der PLUGIN_REGISTRY. Das Modul-Objekt
        wird nicht im Python-Interpreter-Cache (sys.modules) entfernt, da es
        möglicherweise noch von anderem Code referenziert wird.

    @param name  Name des zu entfernenden Plugins.
    @raises KeyError  Wenn kein Plugin mit diesem Namen registriert ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    if name not in PLUGIN_REGISTRY:
        raise KeyError(f"Plugin '{name}' ist nicht registriert.")
    # Plugin aus der Registry entfernen
    del PLUGIN_REGISTRY[name]


def create_plugin_template(name: str, output_dir: str = "plugins/") -> str:
    """
    @brief Erzeugt eine Plugin-Template-Datei mit der Mindest-Struktur.
    @description
        Erstellt eine neue Plugin-Datei im angegebenen Verzeichnis mit einem
        vollständigen Docblock, der ``PLUGIN_NAME``-Variable, einer
        ``setup()``-Funktion und einer Beispiel-Funktion. Das Verzeichnis
        wird ggf. automatisch erstellt.

    @param name        Name des neuen Plugins (wird als Dateiname verwendet).
    @param output_dir  Zielverzeichnis für die Template-Datei.
    @return            Absoluter Pfad zur erstellten Datei.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Verzeichnis erstellen, falls es noch nicht existiert
    os.makedirs(output_dir, exist_ok=True)

    # Dateiname aus Plugin-Namen ableiten (Leerzeichen → Unterstrich)
    safe_name = name.replace(" ", "_").lower()
    file_path = os.path.join(output_dir, f"{safe_name}.py")

    # Template-Inhalt mit vollständigem Docblock
    template = f'''"""
@file {safe_name}.py
@brief Plugin: {name}
@description
    Beschreibung des Plugins hier einfügen.
    Dieses Plugin wurde automatisch mit create_plugin_template() erstellt.
@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0.0
"""

# Plugin-Name: wird von discover_plugins() zur Registrierung verwendet
PLUGIN_NAME = "{safe_name}"

# Plugin-Metadaten
PLUGIN_VERSION = "1.0.0"
PLUGIN_DESCRIPTION = "Beschreibung des Plugins"
PLUGIN_AUTHOR = "Kurt Ingwer"


def setup() -> None:
    """
    @brief Initialisierungsfunktion des Plugins.
    @description
        Wird beim Laden des Plugins aufgerufen. Hier können Ressourcen
        vorbereitet, Konfigurationen geladen oder Abhängigkeiten geprüft werden.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Initialisierungslogik hier einfügen
    pass


def example_function(x: float) -> float:
    """
    @brief Beispiel-Funktion des Plugins.
    @description
        Diese Funktion dient als Vorlage. Ersetze sie durch die eigentliche
        mathematische Funktionalität des Plugins.

    @param x  Eingabewert.
    @return   Berechnetes Ergebnis.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Beispiel: Identitätsfunktion – durch echte Logik ersetzen
    return x
'''

    # Template in Datei schreiben
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(template)

    return os.path.abspath(file_path)
