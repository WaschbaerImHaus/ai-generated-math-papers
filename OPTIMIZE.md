# OPTIMIZE.md - specialist-maths

## Optimierungsvorschläge

### Architektur
- **Modul-Trennung**: Jeden mathematischen Bereich als eigenes Python-Paket
- **Plugin-System**: Neue Gebiete als Plugins hinzufügbar
- **Caching**: Berechnete Ergebnisse zwischenspeichern (functools.lru_cache)
- **Lazy-Loading**: Module erst bei Bedarf laden

### Geschwindigkeit
- **NumPy-Vektorisierung**: Schleifen durch Array-Operationen ersetzen
- **Cython/Numba**: Performance-kritische Teile kompilieren
- **Parallele Berechnungen**: multiprocessing für unabhängige Teilprobleme
- **Symbolische vs. Numerische**: Wahl der richtigen Methode je nach Problem

### Code-Qualität
- **Type Hints**: Vollständige Typ-Annotationen für bessere IDE-Unterstützung
- **Doctests**: Beispiele direkt in Docstrings als ausführbare Tests
- **Protokoll-Logging**: Alle Berechnungsschritte nachvollziehbar machen
- **Fehlerbehandlung**: Spezifische Ausnahmen für mathematische Fehler (Division durch 0, etc.)

### Benutzerfreundlichkeit
- **REPL-Modus**: Interaktive Kommandozeile für Berechnungen
- **LaTeX-Export**: Ergebnisse als LaTeX-Formeln ausgeben
- **Notebook-Integration**: Jupyter-Notebook-kompatibel
- **Schritt-für-Schritt-Modus**: Lösungswege detailliert anzeigen

### Tests
- **Property-Based Testing**: Hypotheses-Bibliothek für automatische Testfälle
- **Fuzzing**: Zufällige Eingaben testen
- **Performance-Benchmarks**: Geschwindigkeitsmessungen für alle Algorithmen

---

## Neue Optimierungsvorschläge (Build 9, 2026-03-08)

### Sicherheit
- **sympify()-Härtung**: In analysis.py alle `sp.sympify(expr_str)` durch
  `parse_expr(expr_str, transformations='all', local_dict={...})` mit Whitelist ersetzen
  → Verhindert Code-Injection bei Web-Nutzung (aktuell MEDIUM-Risiko)

### Geschwindigkeit (konkrete Hot-Spots)
- **Mandelbrot/Julia**: Aktuell NumPy-vektorisiert, aber Python-Schleife über Pixel.
  → numba.jit oder Cython für 50-100x Speedup (wichtig bei großen Auflösungen)
- **Eisenstein-Reihe**: Doppelsumme O(n²) über (m,n)-Gitterpunkte.
  → Nur einen Quadranten berechnen + Symmetrie nutzen → 4x schneller
- **Fourier-Koeffizienten delta**: Polynommultiplikation O(n²) pro Term.
  → FFT-basierte Faltung: O(n log n)
- **Primzahlsieb für Goldbach**: Bei wiederholten Aufrufen Cache der Primzahlen nutzen
  (lru_cache oder Modul-Level-Variable)

### Architektur (nach Build 9 mit 13 Modulen)
- **Paket-Struktur**: `src/` in Python-Package umwandeln (`__init__.py`)
  → Imports werden konsistenter: `from specialist_maths import algebra`
- **Konfigurationsdatei**: `config.py` für globale Konstanten (Epsilon, max_iter, etc.)
  → Verhindert Magic Numbers in verschiedenen Modulen
- **Modulaufteilung linear_algebra.py**: Datei hat >1000 Zeilen
  → Aufteilen in `vectors.py`, `matrix_decomp.py`, `eigensolver.py`
- **Modulaufteilung algebra.py**: Ebenfalls sehr groß
  → Aufteilen in `polynomials.py`, `number_theory.py`, `diophantine.py`

### Visualisierung
- **Interactive Mode**: matplotlib-Widget statt statische Bilder (ipywidgets für Jupyter)
- **Animation**: matplotlib.animation für ODE-Trajektorien und iterative Algorithmen
- **Export-Formate**: SVG/PDF-Export neben PNG für Vektorgrafiken (LaTeX-ready)

### Mathematische Vertiefung
- **Adaptives Gitter für Visualisierungen**: Feineres Gitter in interessanten Regionen
  (Nullstellen, Singularitäten) → Bessere Fraktal-Qualität
- **Arbitrary Precision**: mpmath für Rechnungen mit >64-Bit Genauigkeit
  (wichtig für Riemann-Nullstellen-Verifikation)
