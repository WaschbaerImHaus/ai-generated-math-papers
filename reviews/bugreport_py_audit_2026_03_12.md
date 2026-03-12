# Bugreport: Python-Code-Audit src/py/

**Datum:** 2026-03-12
**Build:** 123
**Umfang:** 161 `.py`-Dateien, ~127.000 Zeilen Code
**Autor:** Michael Fuhrmann

---

## Zusammenfassung

| Schwere | Anzahl |
|---|---|
| KRITISCH | 2 |
| HOCH | 3 |
| MITTEL | 5 |
| GERING | 4 |

---

## KRITISCH

### BUG-PY-001: `eval()` ohne ausreichende Härtung — Code-Injection möglich
- **Datei:** `src/py/webapp/app.py`, Zeile 4444
- **Funktion:** Route `/api/logic/truth_table`
- **Beschreibung:** `eval(formula_str, {"__builtins__": {}}, allowed)` wertet Benutzereingaben aus. Der `{"__builtins__": {}}` Ansatz ist nicht sandbox-sicher — über Attribut-Zugriffsketten (z.B. `().__class__.__bases__[0].__subclasses__()`) kann ein Angreifer trotzdem auf gefährliche Funktionen zugreifen, wenn Objekte im `allowed`-Dict enthalten sind.
- **Empfehlung:** Ersetzen durch einen dedizierten Parser (z.B. `pyparsing`) statt `eval()`.
- **Status:** Teilweise mitigiert (Build 123: Längen-Limit 512 Zeichen, Zeichen-Whitelist, Unterstriche verboten).

### BUG-PY-002: `eval()` in visualization.py — Code-Injection möglich
- **Datei:** `src/py/visualization.py`, Zeile 2071
- **Funktion:** `plot_function_2d_adaptive()`
- **Beschreibung:** Gleiche Problematik wie BUG-PY-001. `eval(func_str, {"__builtins__": {}}, env)` auf numpy-Ausdrücken.
- **Empfehlung:** Ersetzen durch `numexpr.evaluate()` oder `sympy.lambdify()` mit expliziter Whitelist.
- **Status:** Teilweise mitigiert (Build 123: Längen-Limit + Zeichen-Whitelist).

---

## HOCH

### BUG-PY-003: `sympify()` auf Benutzereingaben — SymPy Code-Injection
- **Dateien:**
  - `src/py/webapp/app.py`, Zeilen 540, 5844, 7410
  - `src/py/parallel_compute.py`, Zeilen 54, 69, 86
  - `src/py/formula_simplifier.py`, Zeilen 64, 111, 161
- **Beschreibung:** `sympy.sympify()` ist nicht sandbox-sicher und kann Code-ähnliche Strukturen ausführen. Keine Längen- oder Zeichen-Validierung vor den Aufrufen.
- **Empfehlung:** Ersetzen durch `sympy.parsing.sympy_parser.parse_expr()` mit `transformations` und vorheriger Eingabe-Validierung (Länge, erlaubte Zeichen).

### BUG-PY-004: Plugin-System führt beliebigen Python-Code aus
- **Datei:** `src/py/plugin_system.py`, Zeile 93
- **Funktion:** `load_plugin()`
- **Beschreibung:** `spec.loader.exec_module(module)` führt den Code jeder übergebenen `.py`-Datei ungeprüft aus. Wenn `plugin_dir` aus einer Benutzereingabe stammt oder manipuliert wird, ist Remote Code Execution (RCE) möglich.
- **Empfehlung:** Plugin-Verzeichnis auf einen fest codierten absoluten Pfad (`PLUGIN_BASE_DIR`) beschränken; nur Dateien aus diesem Verzeichnis laden.
- **Status:** Teilweise mitigiert (Build 123: Pfad-Normalisierung via `realpath()`, Dateiname-Whitelist).

### BUG-PY-005: Fehlerbehandlung in `parallel_compute.py` — Fehler werden stillschweigend verschluckt
- **Datei:** `src/py/parallel_compute.py`, Zeilen 136–138, 219–220, 259–261
- **Beschreibung:** Exceptions werden als Werte ins Ergebnis-Dictionary geschrieben (`results[args] = exc`), statt sie zu werfen oder als strukturierten Fehler zurückzugeben. Aufrufer können Fehler nicht zuverlässig von erfolgreichen Berechnungen unterscheiden.
- **Empfehlung:** Rückgabe-Struktur auf `{"success": bool, "value"/"error": ...}` umstellen oder Exceptions explizit weiterwerfen.

---

## MITTEL

### BUG-PY-006: Fehlende Input-Validierung für alle API-Routen
- **Datei:** `src/py/webapp/app.py` (alle 228+ Routes)
- **Beschreibung:** JSON-Eingaben werden ohne Typ-Prüfung, Range-Prüfung oder Längen-Limits direkt verarbeitet (z.B. `data.get('x')`). Ermöglicht Denial-of-Service durch übermäßig große Eingaben sowie unerwartetes Verhalten bei falschen Typen.
- **Beispiel:** `generators = [sympy.sympify(g) for g in gen_strs]` — keine Längen-Begrenzung für `gen_strs`.
- **Empfehlung:** Pydantic-Modelle für alle Route-Eingaben einführen.

### BUG-PY-007: `sagemath_bridge.py` — unsicherer Fallback auf `sympify()`
- **Datei:** `src/py/sagemath_bridge.py`, Zeilen 290–298
- **Beschreibung:** Wenn `sage_eval()` fehlschlägt, wird stillschweigend auf `sp.sympify()` zurückgegriffen — ohne Validierung des Eingabe-Strings. Kombiniert die Risiken aus BUG-PY-003.
- **Empfehlung:** `sage_str` vor dem Fallback explizit validieren; alternativ Fallback entfernen.

### BUG-PY-008: `__init__.py` — direktes Schreiben in `globals()`
- **Datei:** `src/py/__init__.py`, Zeilen 131–134
- **Beschreibung:** `globals()[name] = module` modifiziert den globalen Namensraum direkt, was zu Namespace-Pollution und schwer nachvollziehbaren Überschreibungen führen kann.
- **Empfehlung:** Privates Cache-Dictionary verwenden statt `globals()`.

### BUG-PY-009: `ungeprüfte sympify()`-Konvertierung ohne Exception-Handling
- **Datei:** `src/py/parallel_compute.py`, Zeilen 54, 69, 71, 86
- **Beschreibung:** `sp.sympify()` wird auf Benutzereingaben ohne Try-Except aufgerufen. Eine ungültige Eingabe bricht die gesamte parallele Berechnung ab.
- **Empfehlung:** `sympify()`-Aufrufe in Try-Except wrappen oder Eingaben vorab validieren.

### BUG-PY-010: `__import__()` direkt im f-String — Anti-Pattern
- **Datei:** `src/py/visualization.py`, Zeile 4212
- **Beschreibung:** `f'{__import__("math").log(2) / __import__("math").log(3):.6f}'` lädt ein Modul dynamisch im f-String. Schlechte Lesbarkeit, verhindert statische Analyse.
- **Empfehlung:** Reguläres `import math` am Modulanfang verwenden.

---

## GERING

### BUG-PY-011: `import *` in mehreren Modulen
- **Dateien:** `math_helpers.py`, `linear_algebra.py`, `algebra.py`, `galois_theory.py`, `modular_forms.py`
- **Beschreibung:** `from module import *` verschmutzt den Namensraum und erschwert statische Analyse.
- **Empfehlung:** Explizite Imports oder `__all__`-Listen verwenden.

### BUG-PY-012: Keine Python-Test-Suite vorhanden
- **Beschreibung:** Für die 161 `.py`-Dateien existieren keine Tests unter `tests/`. Kritische Module wie `webapp/app.py`, `parallel_compute.py`, `plugin_system.py` sind vollständig ungetestet.
- **Empfehlung:** `pytest`-Test-Suite anlegen; Mindest-Coverage: `test_webapp_security.py`, `test_sympify_input.py`, `test_plugin_system.py`, `test_parallel_compute.py`.

### BUG-PY-013: Fehlende Docstrings in `formula_simplifier.py`
- **Datei:** `src/py/formula_simplifier.py`
- **Beschreibung:** Mehrere interne Methoden haben keine oder unvollständige DocBlock-Kommentare.
- **Empfehlung:** Deutsche DocBlocks mit `@lastModified`-Timestamp ergänzen (gemäß Projektstandard).

### BUG-PY-014: Build-Nummer in `src/py/build.txt` von Gesamt-Build entkoppelt
- **Beschreibung:** `src/py/build.txt` führt eine separate Build-Nummerierung (aktuell 123), die unabhängig von der Projekt-Build-Nummer der Paper-Reviews läuft. Kann zu Verwirrung führen.
- **Empfehlung:** Entweder zusammenführen oder Namenskonvention (z.B. `py-build.txt`) klarer trennen.
