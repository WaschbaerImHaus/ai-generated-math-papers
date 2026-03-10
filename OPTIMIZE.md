# OPTIMIZE.md - specialist-maths

Optimierungsvorschläge und -status für alle Bereiche des Projekts.
Letzte Aktualisierung: 2026-03-10 (Build 15)

---

## Erledigt (Build 1–15)

### Architektur
- [x] **Modul-Trennung**: `algebra.py` → `algebra_core.py`, `algebra_numbertheory.py`, `algebra_diophantine.py`
- [x] **Modul-Trennung**: `linear_algebra.py` → `vectors.py`, `matrix_ops.py`, `matrix_decomp.py`
- [x] **Paket-Struktur**: `src/__init__.py` angelegt; konsistente Imports via `from specialist_maths import ...`
- [x] **Konfigurationsdatei**: `config.py` mit globalen Konstanten (H_DERIVATIVE_1/2, NEWTON_TOL, etc.)
- [x] **Fehlerbehandlung**: `exceptions.py` mit `MathematicalError`-Hierarchie (`ConvergenceError`, `DomainError`, etc.)
- [x] **Lazy-Loading in `__init__.py`**: Module erst bei Bedarf laden (spart Startzeit bei teilweiser Nutzung)
- [x] **Differential-Geometrie / Tensorrechnung**: `tensor_geometry.py` mit Christoffel-Symbolen, Riemannscher Krümmung, Geodäten, Schwarzschild-Metrik, Differentialformen (Build 15)

### Geschwindigkeit
- [x] **NumPy-Vektorisierung**: Mandelbrot-Fraktal vollständig vektorisiert (numpy-Broadcasting statt Pixel-Schleife)
- [x] **Eisenstein-Reihe Optimierung**: `eisenstein_series_fast()` nutzt Symmetrie des (m,n)-Gitters → ~4x Speedup
- [x] **Multiprocessing**: `goldbach_verification_range()` verwendet `multiprocessing.Pool.map`
- [x] **Caching**: `lru_cache` für `is_prime`, `euler_phi`, `prime_factorization`, `p_adic_valuation`, `bernoulli_number`

### Code-Qualität
- [x] **Type Hints**: Vollständige Typ-Annotationen in `algebra_core.py`, `algebra_numbertheory.py`, `algebra_diophantine.py`, `analysis.py`, `linear_algebra.py`, `vectors.py`, `matrix_ops.py`, `matrix_decomp.py`, `topology.py`, `graph_theory.py`, `config.py`, `exceptions.py`
- [x] **Doctests**: Ausführbare `>>>` Beispiele in allen Hauptmodulen eingefügt
- [x] **sympify()-Härtung**: `safe_parse_expr()` in `analysis.py` – verhindert Code-Injection via Whitelist
- [x] **Logging-System**: `math_logger.py` mit konfigurierbarem Log-Level je Modul
- [x] **Schritt-für-Schritt-Modus**: `step_by_step.py` – Zwischenergebnisse bei iterativen Algorithmen

### Benutzerfreundlichkeit
- [x] **REPL-Modus**: `repl.py` implementiert (interaktive Kommandozeile für Berechnungen)
- [x] **LaTeX-Export**: `latex_export.py` implementiert (Ergebnisse als LaTeX-Formeln ausgeben)

### Tests
- [x] **Property-Based Testing**: `tests/test_property_based.py` mit Hypotheses-Bibliothek
- [x] **Performance-Benchmarks**: `debugging/profile_performance.py` + `PERFORMANCE_REPORT.md`
- [x] **Fuzzing**: `tests/test_fuzzing.py` mit zufälligen Eingaben und Randwert-Tests

---

## Offen / Noch nicht erledigt

### Code-Qualität (hohe Priorität)
- [x] **Type Hints fehlend**: Alle Module vollständig annotiert (Build 20–21)
- [x] **Docblock-Stil vereinheitlichen**: Alles auf Doxygen `@param/@return` (Build 20)
- [x] **Zentrales Test-Discovery-Skript**: `Makefile` mit `make test/test-fast/test-coverage/lint` (Build 20)

### Sicherheit
- [ ] **`analysis.py` `_safe_parse()` Fallback**: Stiller Fallback auf `sp.sympify()` bei Parse-Fehler sollte geloggt werden (schwer debuggbar ohne Log-Ausgabe)
- [ ] **`repl.py` Eingabe-Validierung**: Prüfen ob unsichere Eval-Aufrufe vorhanden sind; ggf. Sandboxing ergänzen

### Geschwindigkeit (konkrete Hot-Spots)
- [ ] **Fourier-Koeffizienten delta-Funktion**: Polynommultiplikation O(n²) pro Term → FFT-basierte Faltung O(n log n) via `numpy.fft.fft`
- [ ] **Cython/Numba JIT**: Performance-kritische Schleifen in `proof_theory.py` (Sieb), `complex_analysis.py` (ζ-Iteration) → Numba-JIT für 10–50x Speedup; Numba derzeit nicht installiert
- [ ] **Symbolische vs. Numerische Wahl**: Automatische Entscheidung je nach Problem-Typ (kleine Polynome → SymPy symbolisch; große Matrizen → NumPy numerisch)
- [ ] **`prime_factorization()` + `euler_phi()`**: Zwei separate Primfaktor-Traversierungen; `euler_phi` könnte `prime_factorization()` direkt wiederverwenden
- [x] **`bisection()` + `newton_raphson()` Hybridmethode**: `brent_method()` implementiert (Build 19)
- [ ] **Christoffel-Symbole cachen**: `christoffel_symbols()` wird in `riemann_tensor()` mehrfach an benachbarten Punkten aufgerufen – Memoization könnte Rechenzeit halbieren

### Architektur
- [ ] **Plugin-System**: Neue mathematische Gebiete als Python-Packages einbinden ohne Core zu ändern (z.B. via `importlib` + Konfigurations-Registry)
- [x] **Differentialformen-Modul vertiefen**: de Rham-Kohomologie + Stokes-Satz in `topology.py` (Build 25)
- [ ] **Kategorientheorie-Modul**: Objekte, Morphismen, Funktoren, natürliche Transformationen – Verbindung zu abstrakter Algebra

### Visualisierung
- [ ] **Interactive Mode**: `matplotlib`-Widgets statt statischer Bilder (`ipywidgets` für Jupyter, `matplotlib.widgets` für CLI)
- [ ] **Animation**: `matplotlib.animation` für ODE-Trajektorien und iterative Algorithmen (Newton-Konvergenz, Geodäten auf Sphäre)
- [ ] **Export-Formate**: SVG/PDF-Export neben PNG für Vektorgrafiken (LaTeX-ready)
- [ ] **Adaptives Gitter**: Feineres Gitter in interessanten Regionen (Nullstellen, Singularitäten)
- [ ] **3D-Visualisierung Krümmung**: Gaußsche Krümmung als Farbkodierung auf Mannigfaltigkeiten (Sphäre, Torus, Sattelfläche) in `visualization.py` ergänzen

### Mathematische Vertiefung
- [ ] **Arbitrary Precision mit mpmath**: `mpmath` für Rechnungen mit >64-Bit-Genauigkeit (wichtig für Riemann-Nullstellen-Verifikation mit >100 Stellen)
- [ ] **Notebook-Integration**: Jupyter-kompatible Ausgaben (HTML-Darstellung für Matrizen, Polynome, LaTeX-Rendering in Notebooks)
- [ ] **Symplektische Geometrie**: Hamilton-Mechanik, symplektische Mannigfaltigkeiten, Poisson-Klammern (Erweiterung von `tensor_geometry.py`)
- [ ] **Faserräume und Verbindungen**: Prinzipal-Faserbündel, Gauge-Theorie (Yang-Mills) – Bezug zu Standardmodell der Teilchenphysik
- [ ] **Spinor-Rechnung**: Clifford-Algebren, Dirac-Gleichung – Erweiterung von `tensor_geometry.py` Richtung Quantenfeldtheorie
- [ ] **Algebraische Topologie**: Homologiegruppen, Betti-Zahlen, Euler-Charakteristik (Ergänzung zu `topology.py`)

### Neue Ideen (entdeckt 2026-03-10)
- [ ] **Paralleles symbolisches Rechnen**: SymPy-Berechnungen (Grenzwerte, Integrale) parallelisieren via `concurrent.futures.ProcessPoolExecutor`
- [ ] **Benchmark-Regressions-Tests**: Automatischer Vergleich von Laufzeiten zwischen Builds; Alarm bei >10% Verlangsamung
- [ ] **Interoperabilität SageMath**: Export/Import-Schnittstelle zu SageMath für Berechnungen die SymPy nicht beherrscht
- [ ] **Automatische Formel-Vereinfachung**: Ergebnisse standardmäßig durch `sp.simplify()` + `sp.nsimplify()` vereinfachen und schönste Darstellung wählen
- [ ] **Geodäten-Visualisierung**: Aus `tensor_geometry.geodesic_equation()` direkt 3D-Plots auf Sphäre/Torus erzeugen via `visualization.py`
- [ ] **Numerische Stabilitäts-Analyse**: Konditionszahlen aller Matrixoperationen automatisch protokollieren; Warnung bei Kondition > 1e10

---

## Prioritäten (empfohlene Reihenfolge)

1. Type Hints in den verbleibenden Modulen ergänzen (einfach, hoher Nutzen für IDE-Support)
2. Docblock-Stil vereinheitlichen (Wartbarkeit)
3. Zentrales Test-Discovery-Skript (Entwicklungskomfort)
4. Geodäten-Visualisierung (zeigt tensor_geometry.py in Aktion)
5. Bisection/Newton Hybridmethode (numerische Verbesserung)
6. mpmath Arbitrary Precision (für Millennium-Problem-Recherche wichtig)
