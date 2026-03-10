# plugin_system.py – Dokumentation

## Übersicht

Das Modul `plugin_system.py` implementiert ein dynamisches Plugin-System für
das specialist-maths Projekt. Es ermöglicht das Laden, Registrieren und Verwalten
von mathematischen Erweiterungsmodulen zur Laufzeit.

## Grundkonzept

Ein Plugin ist eine eigenständige Python-Datei, die in das `plugins/`-Verzeichnis
gelegt wird. Das Plugin-System erkennt diese Dateien automatisch und stellt sie
über eine zentrale Registry zur Verfügung.

```
plugins/
├── example_plugin.py    ← Beispiel-Plugin
├── README.md            ← Anleitung
└── mein_plugin.py       ← Eigenes Plugin
```

## API-Referenz

### `register_plugin(name: str, module: Any) → None`

Registriert ein Modul-Objekt unter einem Namen in der globalen Registry
`PLUGIN_REGISTRY`.

```python
import types
m = types.ModuleType("meinmod")
m.compute = lambda x: x**2
register_plugin("meinmod", m)
```

### `load_plugin(path: str) → Any`

Lädt eine Python-Datei als Modul-Objekt, ohne den Python-Suchpfad zu verändern.
Intern wird `importlib.util.spec_from_file_location` verwendet.

```python
module = load_plugin("plugins/example_plugin.py")
result = module.perfect_number_check(28)  # True
```

**Fehler:**
- `FileNotFoundError` – Datei nicht gefunden
- `ImportError` – Modul-Spec kann nicht erstellt werden

### `discover_plugins(plugin_dir: str = "plugins/") → list`

Sucht im angegebenen Verzeichnis nach `.py`-Dateien und lädt alle als Plugins.
`__init__.py` wird übersprungen. Fehlerhafte Plugins werden mit Warnung übersprungen.

```python
loaded = discover_plugins("plugins/")
# → ['example_plugin', ...]
```

### `get_plugin(name: str) → Any`

Gibt das Modul-Objekt eines registrierten Plugins zurück.

```python
plugin = get_plugin("example_plugin")
print(plugin.digit_sum(1234))  # 10
```

**Fehler:**
- `KeyError` – Plugin nicht registriert

### `list_plugins() → list`

Gibt eine sortierte Liste aller Registrierungsnamen zurück.

```python
print(list_plugins())
# → ['example_plugin']
```

### `unload_plugin(name: str) → None`

Entfernt ein Plugin aus der Registry. Das Python-Modul-Objekt bleibt im Speicher.

```python
unload_plugin("example_plugin")
```

**Fehler:**
- `KeyError` – Plugin nicht registriert

### `create_plugin_template(name: str, output_dir: str = "plugins/") → str`

Erstellt eine Template-Datei mit der Mindest-Struktur eines Plugins.
Gibt den absoluten Pfad der erstellten Datei zurück.

```python
path = create_plugin_template("mein_plugin")
# → '/absoluter/pfad/plugins/mein_plugin.py'
```

## Globale Variable

### `PLUGIN_REGISTRY: Dict[str, Any]`

Die zentrale Registry aller geladenen Plugins. Kann direkt gelesen werden:

```python
from plugin_system import PLUGIN_REGISTRY
print(list(PLUGIN_REGISTRY.keys()))
```

## Plugin-Schnittstelle

Ein Plugin-Modul sollte folgende Konventionen einhalten:

| Element | Pflicht | Beschreibung |
|---------|---------|--------------|
| `PLUGIN_NAME` | empfohlen | Registrierungsname (sonst Dateiname) |
| `PLUGIN_VERSION` | optional | Versionsnummer |
| `PLUGIN_DESCRIPTION` | optional | Kurzbeschreibung |
| `PLUGIN_AUTHOR` | optional | Autorenname |
| `setup()` | optional | Initialisierungsfunktion |

## Beispiel: Vollständiges Plugin

```python
# plugins/statistics_ext.py
"""Erweiterungs-Plugin: Statistische Funktionen."""

PLUGIN_NAME = "statistics_ext"
PLUGIN_VERSION = "1.0.0"
PLUGIN_DESCRIPTION = "Erweiterte statistische Funktionen"
PLUGIN_AUTHOR = "Kurt Ingwer"

def setup():
    pass

def trimmed_mean(data: list, percent: float = 0.1) -> float:
    """Getrimmtes Mittel (entfernt extremste Werte)."""
    n = len(data)
    k = int(n * percent)
    sorted_data = sorted(data)
    trimmed = sorted_data[k:n-k] if k > 0 else sorted_data
    return sum(trimmed) / len(trimmed)
```

## Implementierungsdetails

- **Isolation**: Jedes Plugin wird als eigenständiges Modul-Objekt geladen
- **Keine Seiteneffekte**: Die Registry wird nicht persistent gespeichert
- **Fehlertoleranz**: `discover_plugins()` überspringt fehlerhafte Plugins
- **Überschreiben**: Ein Aufruf von `register_plugin()` mit gleichem Namen
  überschreibt den vorhandenen Eintrag

## Tests

Tests befinden sich in `tests/test_plugin_system.py` (30+ Tests):

```bash
python3 -m pytest tests/test_plugin_system.py -v
```
