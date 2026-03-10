# OPTIMIZE.md - specialist-maths

Optimierungsvorschläge und -status für alle Bereiche des Projekts.
Letzte Aktualisierung: 2026-03-10

---

## Bereits erledigt (Build 1–13)

### Architektur
- [x] **Modul-Trennung**: `algebra.py` → `algebra_core.py`, `algebra_numbertheory.py`, `algebra_diophantine.py`
- [x] **Modul-Trennung**: `linear_algebra.py` → `vectors.py`, `matrix_ops.py`, `matrix_decomp.py`
- [x] **Paket-Struktur**: `src/__init__.py` angelegt; konsistente Imports via `from specialist_maths import ...`
- [x] **Konfigurationsdatei**: `config.py` mit globalen Konstanten (H_DERIVATIVE_1/2, NEWTON_TOL, BISECTION_TOL, SIMPSON_N, MAX_ITERATIONS, etc.)
- [x] **Fehlerbehandlung**: `exceptions.py` mit `MathematicalError`-Hierarchie (`ConvergenceError`, `DomainError`, etc.)

### Geschwindigkeit
- [x] **NumPy-Vektorisierung**: Mandelbrot-Fraktal vollständig vektorisiert (numpy-Broadcasting statt Pixel-Schleife)
- [x] **Eisenstein-Reihe Optimierung**: `eisenstein_series_fast()` nutzt Symmetrie des (m,n)-Gitters → ~4x Speedup gegenüber naiver Doppelsumme
- [x] **Multiprocessing**: `goldbach_verification_range()` verwendet `multiprocessing.Pool.map` für parallele Primzahl-Checks
- [x] **Caching**: `lru_cache` für `is_prime`, `euler_phi`, `prime_factorization`, `p_adic_valuation`, `bernoulli_number`

### Code-Qualität
- [x] **Type Hints**: Vollständige Typ-Annotationen in `algebra_core.py`, `algebra_numbertheory.py`, `algebra_diophantine.py`, `analysis.py`, `linear_algebra.py`, `vectors.py`, `matrix_ops.py`, `matrix_decomp.py`, `topology.py`, `graph_theory.py`, `config.py`, `exceptions.py`
- [x] **Doctests**: Ausführbare `>>>` Beispiele in allen Hauptmodulen eingefügt (diese Aufgabe, Build 13)
- [x] **sympify()-Härtung**: `safe_parse_expr()` in `analysis.py` – verwendet `parse_expr()` mit Whitelist statt `sympify()` (verhindert Code-Injection)

### Benutzerfreundlichkeit
- [x] **REPL-Modus**: `repl.py` implementiert (interaktive Kommandozeile für Berechnungen)
- [x] **LaTeX-Export**: `latex_export.py` implementiert (Ergebnisse als LaTeX-Formeln ausgeben)

### Tests
- [x] **Property-Based Testing**: `tests/test_property_based.py` mit Hypotheses-Bibliothek
- [x] **Performance-Benchmarks**: `debugging/profile_performance.py` + `PERFORMANCE_REPORT.md`

---

## Offen / Noch nicht erledigt

### Code-Qualität (hohe Priorität)
- [ ] **Type Hints**: Fehlende Annotationen in `modular_forms.py`, `complex_analysis.py`, `ode.py`, `statistics_math.py`, `fourier.py`, `analytic_number_theory.py`, `proof_theory.py`
- [ ] **Logging-System**: `logging`-Modul statt `print()` für Debug-Ausgaben; konfigurierbarer Log-Level je Modul
- [ ] **Schritt-für-Schritt-Modus**: Zwischenergebnisse bei iterativen Algorithmen ausgeben (Newton, Bisektion, Runge-Kutta etc.)

### Tests
- [ ] **Fuzzing**: Zufällige Eingaben mit `atheris` oder `hypothesis.extra.numpy` testen; besonders für Randwerte (n=0, negative Zahlen, sehr große Zahlen)

### Geschwindigkeit (konkrete Hot-Spots)
- [ ] **Fourier-Koeffizienten delta**: Polynommultiplikation O(n²) pro Term → FFT-basierte Faltung O(n log n) via `numpy.fft.fft`
- [ ] **Cython/Numba**: Performance-kritische Schleifen in `proof_theory.py` (Sieb), `complex_analysis.py` (ζ-Iteration) → Numba-JIT für 10–50x Speedup; Numba derzeit nicht installiert
- [ ] **Symbolische vs. Numerische Wahl**: Automatische Entscheidung je nach Problem-Typ (kleine Polynome → SymPy symbolisch; große Matrizen → NumPy numerisch)

### Architektur
- [ ] **Lazy-Loading in `__init__.py`**: Module erst bei Bedarf laden statt alle beim Import (spart Startzeit bei teilweiser Nutzung)
- [ ] **Plugin-System**: Neue mathematische Gebiete als Python-Packages einbinden ohne Core zu ändern
- [ ] **Differential-Geometrie / Tensorrechnung**: Neues Modul `tensor_geometry.py` mit Christoffel-Symbolen, Riemannscher Krümmung, Lie-Ableitungen

### Visualisierung
- [ ] **Interactive Mode**: `matplotlib`-Widgets statt statischer Bilder (`ipywidgets` für Jupyter)
- [ ] **Animation**: `matplotlib.animation` für ODE-Trajektorien und iterative Algorithmen (Newton-Konvergenz, Fraktal-Zoom)
- [ ] **Export-Formate**: SVG/PDF-Export neben PNG für Vektorgrafiken (LaTeX-ready)
- [ ] **Adaptives Gitter**: Feineres Gitter in interessanten Regionen (Nullstellen, Singularitäten) → bessere Fraktal-Qualität

### Mathematische Vertiefung
- [ ] **Arbitrary Precision**: `mpmath` für Rechnungen mit >64-Bit-Genauigkeit (wichtig für Riemann-Nullstellen-Verifikation mit >100 Stellen)
- [ ] **Notebook-Integration**: Jupyter-Notebook-kompatible Ausgaben (HTML-Darstellung für Matrizen, Polynome, LaTeX in Notebooks)
- [ ] **Differential-Formen**: Äußere Ableitung, Stokes-Satz, de Rham-Kohomologie als neues Modul

---

## Neue Optimierungsideen (entdeckt beim Durchlesen, 2026-03-10)

### Sicherheit
- [ ] **analysis.py `_safe_parse()`**: Fallback auf `sp.sympify()` bei Parse-Fehler sollte geloggt werden (stiller Fallback ist schwer debuggbar)
- [ ] **repl.py Eingabe-Validierung**: Prüfen ob auch dort unsichere Eval-Aufrufe vorhanden sind

### Performance
- [ ] **`prime_factorization()` + `euler_phi()`**: Derzeit zwei separate Primfaktor-Traversierungen; `euler_phi` könnte `prime_factorization()` wiederverwenden statt selbst zu iterieren
- [ ] **`bisection()` + `newton_raphson()`**: Hybridmethode (Illinois/Brent) wäre konvergenzgarantiert UND schnell → Kombination aus beiden

### Wartbarkeit
- [ ] **Einheitlicher Docblock-Stil**: `statistics_math.py` nutzt `:param:` (Sphinx), alle anderen `@param` (Doxygen) → Vereinheitlichung auf einen Stil
- [ ] **Zentrales Test-Discovery-Skript**: Ein Skript das sowohl `pytest`, Doctests und Property-Tests in einem Lauf ausführt
