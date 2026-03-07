# OPTIMIZED_WORKER.md - Mathematik-Spezialist
## Aktualisiert: 2026-03-05

## Was ich über dieses Projekt weiß

### Ziel
Ich werde zum Mathematik-Spezialisten ausgebildet. Das Projekt dient als
Wissensbasis, Übungsumgebung und Dokumentationssammlung für mathematische Konzepte.

### Implementierter Stack
- **Python 3.13** mit sympy 1.13, numpy 2.2, scipy 1.15
- **Algebra-Modul** (`src/algebra.py`): Polynome, Gleichungslöser, Zahlentheorie
- **Analysis-Modul** (`src/analysis.py`): Differentiation, Integration, Taylor, Nullstellen
- **Lineare Algebra** (`src/linear_algebra.py`): Vektoren, Matrizen, Eigenwerte, Gram-Schmidt
- **63 Tests** alle grün (Stand 2026-03-05)

### Mathematische Erkenntnisse dieser Session

#### Numerische Stabilität (wichtige Lektion)
- Für die 1. numerische Ableitung: optimales h ≈ ε^(1/3) ≈ 6×10⁻⁶
- Für die 2. Ableitung: h ≈ ε^(1/4) ≈ 1.2×10⁻⁴ (wegen h² im Nenner)
- Finite Differenzen k-ter Ordnung sind für k > 8 INSTABIL (h^k Nenner)
- Lösung: SymPy für symbolische Taylor-Ableitungen verwenden

#### Algorithmen-Qualität
- **Horner-Schema**: O(n) statt O(n²) für Polynomauswertung
- **Gauss mit Pivotsuche**: numerisch stabiler als ohne Pivotsuche
- **Simpson-Regel**: O(h⁴) Fehler, deutlich besser als Trapez O(h²)
- **Newton-Raphson**: quadratische Konvergenz vs. lineare bei Bisektion

### Geplante Erweiterungen (Priorität)
1. Statistik-Modul (Verteilungen, Hypothesentests)
2. Differentialgleichungen (Runge-Kutta, Laplace)
3. Komplexe Analysis (Fourier-Transformation)
4. Visualisierung (matplotlib-Funktionsplotter)
5. Kryptographie-Anwendungen (RSA via Zahlentheorie)

## Brainstorming: Verbindungen

### "Fraktal" → Projekt
Mandelbrot-Menge verwendet komplexe Zahlen z_{n+1} = z_n² + c.
Konvergenzanalyse mit Funktionenfolgen aus der Analysis.
Implementierung mit numpy für schnelle Berechnung.

### "Fourier" → Projekt
Fourier-Transformation zerlegte Signale in Frequenzanteile.
Basis: orthogonale Funktionen (wie Gram-Schmidt für Vektoren).
Implementierung: scipy.fft oder eigene DFT.

### "Primzahlen" → Kryptographie
RSA-Kryptographie basiert auf Primfaktorzerlegung (bereits implementiert).
Miller-Rabin-Primzahltest für große Zahlen effizienter als Probedivision.

### "Topologie" → Fixpunktsätze
Banachscher Fixpunktsatz begründet Newton-Raphson-Konvergenz.
Brouwerscher Fixpunktsatz: Existenz von Gleichgewichten in Spieltheorie.

## Quellen & Recherche
- Numerik der Differentialgleichungen: Hairer, Nørsett, Wanner
- Lineare Algebra: Gilbert Strang "Introduction to Linear Algebra"
- Zahlentheorie: Hardy & Wright "An Introduction to the Theory of Numbers"
- Numerische Analysis: Stoer & Bulirsch
