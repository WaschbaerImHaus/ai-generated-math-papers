# Plugin-System – Anleitung

## Überblick

Das specialist-maths Plugin-System ermöglicht das **dynamische Laden mathematischer
Erweiterungsmodule**, ohne den Kern-Code zu verändern. Plugins sind eigenständige
Python-Dateien, die in das `plugins/`-Verzeichnis gelegt werden.

## Aufbau eines Plugins

Jede Plugin-Datei muss folgende Mindeststruktur aufweisen:

```python
# Plugin-Name für die Registry
PLUGIN_NAME = "mein_plugin"

def setup():
    """Optionale Initialisierung."""
    pass

def meine_funktion(x):
    """Eigentliche Funktionalität."""
    return x * 2
```

### Pflichtfelder

| Element | Typ | Beschreibung |
|---------|-----|--------------|
| `PLUGIN_NAME` | str | Eindeutiger Bezeichner für die Registry |

### Optionale Felder

| Element | Typ | Beschreibung |
|---------|-----|--------------|
| `PLUGIN_VERSION` | str | Versionsnummer |
| `PLUGIN_DESCRIPTION` | str | Kurzbeschreibung |
| `PLUGIN_AUTHOR` | str | Autorenname |
| `setup()` | Funktion | Wird beim Laden aufgerufen |

## Verwendung

### Alle Plugins laden

```python
from src.plugin_system import discover_plugins, get_plugin

# Alle Plugins im plugins/-Verzeichnis laden
loaded = discover_plugins("plugins/")
print("Geladene Plugins:", loaded)
```

### Einzelnes Plugin laden

```python
from src.plugin_system import load_plugin, register_plugin

module = load_plugin("plugins/mein_plugin.py")
register_plugin("mein_plugin", module)
```

### Plugin abrufen und nutzen

```python
plugin = get_plugin("example_plugin")
result = plugin.perfect_number_check(28)  # True
```

### Alle registrierten Plugins auflisten

```python
from src.plugin_system import list_plugins
print(list_plugins())
```

### Plugin-Template erstellen

```python
from src.plugin_system import create_plugin_template
path = create_plugin_template("mein_neues_plugin", "plugins/")
print(f"Template erstellt: {path}")
```

## Beispiel-Plugin

Das mitgelieferte `example_plugin.py` stellt folgende Funktionen bereit:

- `perfect_number_check(n)` – Vollkommene Zahlen (6, 28, 496, ...)
- `digit_sum(n)` – Quersumme einer Zahl
- `collatz_length(n)` – Länge der Collatz-Folge
- `abundant_numbers(limit)` – Alle abundanten Zahlen bis `limit`

## Eigene Plugins schreiben

1. Neue Datei in `plugins/` anlegen (oder `create_plugin_template()` nutzen)
2. `PLUGIN_NAME` setzen
3. Mathematische Funktionen implementieren
4. Mit `discover_plugins()` automatisch laden lassen

## Hinweise

- Fehlerhafte Plugins werden übersprungen (Warnung in der Konsole)
- `__init__.py` wird niemals als Plugin geladen
- `unload_plugin(name)` entfernt ein Plugin aus der Registry
- Die Registry wird nicht persistent gespeichert (flüchtig im Speicher)
