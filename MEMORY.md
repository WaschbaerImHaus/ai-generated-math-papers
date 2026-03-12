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
