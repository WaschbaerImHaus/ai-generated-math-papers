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
