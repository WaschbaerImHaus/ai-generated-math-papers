# MEMORY.md - Projekterinnerungen

## Session 2026-03-05

### Projektstart
- Erstes Setup: Verzeichnisstruktur, Bibliotheken, alle Pflichtdateien
- Python 3.13.7 auf LXC (Proxmox), pip mit --break-system-packages

### Implementiert
- algebra.py: Polynome, Gleichungslöser, vollständige Zahlentheorie
- analysis.py: Numerische Diff/Integration, Newton-Raphson, Taylor (SymPy)
- linear_algebra.py: Vector, Matrix, Gram-Schmidt, QR-Iteration

### Lektionen
- h=1e-7 zu klein für 2. Ableitung -> Auslöschungsfehler -> h=1.5e-4
- Finite Differenzen k-ter Ordnung für k>8 numerisch instabil -> SymPy nutzen
- Tests immer auf mathematische Korrektheit prüfen (Kommentarfehler im Test!)

### Nächste Sitzung: Geplant
- Statistik-Modul (Normalverteilung, t-Test, Bayes)
- Differentialgleichungen (Runge-Kutta 4. Ordnung)
- Visualisierungsmodul (matplotlib Plotter)

## Session 2026-03-12 (Build 149) — ABSCHLUSS BATCHES 11–18

### Papers-Gesamtstand (Build 149)
- **72 Papers** insgesamt: Batches 1–10 (Papers 1–39) + Batches 11–18 (Papers 40–71)
- **Batch 18 (Papers 68–71):** Yang-Mills, Hodge, Ricci-Fluss, P-vs-NP — DRUCKREIF nach Audit
- Alle Papers reviewed, alle Batches abgeschlossen

### Batch 11–18 Themen (Papers 40–71)
| Batch | Papers | Thema | Status |
|-------|--------|-------|--------|
| 11 | 40–43 | Algebraische Topologie, Homologie, Kohomologie | DRUCKREIF |
| 12 | 44–47 | Zahlentheorie: Siebert-Primzahlen, Waring, Legendre | DRUCKREIF |
| 13 | 48–51 | Diophantische Approximation, Baker, Roth, Schmidt | DRUCKREIF |
| 14 | 52–55 | Algebraische Geometrie, Weil-Vermutungen | DRUCKREIF |
| 15 | 56–59 | Analytische Zahlentheorie, L-Funktionen, BSD | DRUCKREIF |
| 16 | 60–63 | Kombinatorik, Ramsey, Extremale Graphentheorie | DRUCKREIF |
| 17 | 64–67 | Dynamische Systeme, Ergodische Theorie, Chaos | DRUCKREIF |
| 18 | 68–71 | Yang-Mills, Hodge, Ricci-Fluss, P-vs-NP | DRUCKREIF ✓ |

### Batch-18-Audit-Korrekturen (Build 149)
- paper69: `\title{}` fehlte in EN + DE — BEHOBEN
- paper71: `\documentclass[reqno,11pt]` → `[12pt,a4paper]` in EN + DE — BEHOBEN
- paper71: `\tableofcontents` in EN + DE eingefügt — BEHOBEN
- paper70: `remark`/`bemerkung` unter korrektem `\theoremstyle{remark}` — BEHOBEN

### Kritische mathematische Distinktionen Batch 18
- Yang-Mills Massenlücke: OFFEN (Conjecture) ✓
- Hodge-Vermutung (rational): OFFEN ✓
- Ganzzahlige Hodge (Voisin 2002): WIDERLEGT ✓
- Poincaré dim 3 (Perelman 2003): BEWIESEN ✓
- Smooth 4D Poincaré: OFFEN ✓
- P≠NP: OFFEN (Conjecture) ✓

## Session 2026-03-12 (Build 122)

### Neue Module implementiert
- `erdos_moser.py`: Erdős-Moser-Vermutung, Bernoulli-Zahlen, Moser-Schranke
- `artin_primitive_roots.py`: Primitive Wurzeln, Artin-Konstante, Hooley GRH
- `erdos_selfridge.py`: Binomialkoeffizienten, Sylvester-Theorem, Primzahlpotenzen
- 108 Tests in `tests/test_moser_artin_selfridge.py`

### Wichtige Hinweise
- sympy B_1 = +1/2 (nicht -1/2 wie in manchen Lehrbüchern)
- Erdős-Selfridge gilt für k≥2 (k=1 trivial: C(n,1)=n)
- verify_up_to() muss k ab 2 beginnen, sonst falsche Gegenbeispiele

### Beweisstand-Korrektheit
- Alle drei Conjectures korrekt als OFFEN markiert
- Hooley: unter GRH bewiesen (nicht unbedingt!)
- Heath-Brown: unbedingt nur für mindestens eine aus {2,3,5}
