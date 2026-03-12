# CLAUDE.md - Projekt: specialist-maths

## Projektziel
Dieses Projekt dient dazu, Claude als Mathematik-Spezialisten weiterzuentwickeln.
Das Ziel ist das Erlernen, Üben und Anwenden von Mathematik auf allen Gebieten,
von Grundlagen bis hin zu fortgeschrittenen Themen.
**Fernziel**: Untersuchung der Millennium-Probleme (RH, Goldbach, BSD, P≠NP).

## Technologie-Stack (Stand Build 131)
- **Sprache**: Python 3.13.x
- **Bibliotheken**:
  - sympy 1.14 (symbolische Mathematik)
  - numpy 2.2 (numerische Berechnungen)
  - scipy 1.17 (wissenschaftliche Algorithmen)
  - matplotlib 3.10 (Visualisierung)
  - mpmath (hochpräzise Arithmetik)
  - numba (JIT-Beschleunigung)
- **Tests**: pytest 9.0 + pytest-xdist (-n auto) + pytest-timeout

## Themengebiete (Lernpfad)
1. **Algebra** - Gleichungen, Polynome, abstrakte Algebra ✓
2. **Analysis** - Differential- und Integralrechnung, Grenzwerte ✓
3. **Lineare Algebra** - Vektoren, Matrizen, Eigenwerte ✓
4. **Zahlentheorie** - Primzahlen, Kongruenzen, Diophantische Gleichungen ✓
5. **Statistik & Wahrscheinlichkeit** - Verteilungen, Hypothesentests ✓
6. **Differentialgleichungen** - ODE, PDE, Systeme ✓
7. **Diskrete Mathematik** - Graphen, Kombinatorik, Logik ✓
8. **Numerische Methoden** - Approximation, Interpolation, Iteration ✓
9. **Komplexe Analysis** - Komplexe Zahlen, Fourier-Transformation ✓
10. **Topologie & Geometrie** - Räume, Mannigfaltigkeiten ✓
11. **L-Funktionen & analytische Zahlentheorie** ✓
12. **Modulformen & elliptische Kurven** ✓
13. **p-adische Zahlen & Iwasawa-Theorie** ✓
14. **Algebraische Zahlentheorie** (Paper 37) ✓
15. **Darstellungstheorie & Langlands-Programm** (geplant)

## Dateistruktur
- `/src/` - Python-Module für mathematische Algorithmen (93+ Module)
- `/tests/` - Unit-Tests (7000+ Tests, TDD)
- `/research/` - Rechercheergebnisse, OPTIMIZED_WORKER.md, Notizen
- `/debugging/` - Debug-Skripte
- `/build/` - Kompilierte/ausführbare Dateien
- `/dev-log/` - Entwicklungs-Tagebuch
- `/claude-generated/` - Sonstige generierte Dateien
- `/papers/` - LaTeX-Papers (batch9 aktuell); reviewed/ (Batches 1-9 fertig)
- `/reviews/` - Externe und Selbst-Reviews

## Aktueller Fortschritt (Build 131)
- [x] Projektstruktur angelegt
- [x] Python-Bibliotheken installiert (aktuell: sympy 1.14, numpy 2.2, scipy 1.17)
- [x] Algebra-Modul implementiert (algebra.py, algebraic_structures.py, ...)
- [x] Analysis-Modul implementiert (analysis.py, complex_analysis.py, ...)
- [x] Lineare Algebra-Modul implementiert (linear_algebra.py)
- [x] Zahlentheorie-Module (proof_theory.py, analytic_number_theory.py, ...)
- [x] Statistik-Modul (statistics_math.py)
- [x] ODE/PDE-Module (ode.py, pde.py)
- [x] Fourier/FFT-Modul (fourier.py)
- [x] Numerik-Modul (numerical_methods.py)
- [x] Modulformen (modular_forms.py)
- [x] p-adische Zahlen (p_adic.py)
- [x] L-Funktionen (l_functions.py)
- [x] Beweisversuche (beweisversuche.py): Giuga+Lehmer 3-Prim-Korollar, Erdős-Straus, Kurepa
- [x] Papers: 38 Papers in 9 Batches (alle reviewed, in papers/reviewed/)
- [x] Test-Suite aufgebaut (7000+ Tests, pytest-xdist)
- [x] Webapp (Flask, 220+ Routen, 46+ Templates)
- [x] Visualisierung (visualization.py: 2D/3D, Fraktale, Geodäten, Animationen)
- [ ] Langlands-Programm (Galois-Darstellungen ρ: Gal(Q̄/Q)→GL_n(ℤ_p))
- [x] Batch 10: Paper 39 Langlands-Programm (EN + DE) ✓
- [x] Batch 11: Paper 40 Giuga 4-Prim-Vermutung (EN + DE) ✓

## Besondere Hinweise
- Alle mathematischen Beweise werden kommentiert und erklärt
- Code soll auch als Lernmaterial dienen
- Jede Methode enthält mathematische Dokumentation
- erzeuge für jede programmcode datei ein .md datei mit umfassenden erläuterungen
- Formeln können im KaTeX-Format angegeben werden
- **Review-Workflow**: Neue Papers → review_batchN_YYYY-MM-DD.md → Bugs beheben → nach reviewed/batchN/ verschieben
- **Autor**: Michael Fuhrmann in allen Dateien
- **Theorem vs Conjecture**: Unbewiesene Aussagen IMMER als Conjecture/Vermutung deklarieren!
