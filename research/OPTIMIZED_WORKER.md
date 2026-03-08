# OPTIMIZED_WORKER.md - Mathematik-Spezialist
## Aktualisiert: 2026-03-08

## Was ich über dieses Projekt weiß

### Ziel
Ich werde zum Mathematik-Spezialisten ausgebildet. Das Projekt dient als
Wissensbasis, Übungsumgebung und Dokumentationssammlung für mathematische Konzepte.
Jedes Modul enthält Algorithmen mit ausführlicher mathematischer Dokumentation,
sodass der Code auch als Lernmaterial dient.

### Implementierter Stack
- **Python 3.13** mit sympy 1.13, numpy 2.2, scipy 1.15, pytest 9.0.2
- **Algebra** (`src/algebra.py`): Polynome, Gleichungslöser, Zahlentheorie
- **Analysis** (`src/analysis.py`): Differentiation, Integration, Taylor, Nullstellen
- **Lineare Algebra** (`src/linear_algebra.py`): Vektoren, Matrizen, Eigenwerte, Gram-Schmidt
- **Statistik** (`src/statistics_math.py`): Deskriptiv, Verteilungen, Tests, Bayes, Monte-Carlo
- **ODE** (`src/ode.py`): Euler, RK4, RK45 adaptiv, Laplace, inverse Laplace
- **134 Tests** alle grün (Stand 2026-03-07)

### Mathematische Erkenntnisse

#### Numerische Stabilität (wichtige Lektion)
- Für die 1. numerische Ableitung: optimales h ≈ ε^(1/3) ≈ 6×10⁻⁶
- Für die 2. Ableitung: h ≈ ε^(1/4) ≈ 1.2×10⁻⁴ (wegen h² im Nenner)
- Finite Differenzen k-ter Ordnung sind für k > 8 INSTABIL (h^k Nenner)
- Lösung: SymPy für symbolische Taylor-Ableitungen verwenden

#### Algorithmen-Qualität
- **Horner-Schema**: O(n) statt O(n²) für Polynomauswertung
- **Gauss mit Pivotsuche**: numerisch stabiler als ohne Pivotsuche
- **Simpson-Regel**: O(h⁴) Fehler, deutlich besser als Trapez O(h²)
- **Newton-Raphson**: quadratische Konvergenz, Ableitung wird numerisch bestimmt
- **RK45 adaptiv**: Fehlerkontrolle mit Schrittweitensteuerung
- **Stehfest-Algorithmus**: Numerische inverse Laplace-Transformation

#### Statistische Methoden
- **t-Test**: Welch-Variante für ungleiche Varianzen verwendet (robuster)
- **Chi-Quadrat**: Anwendbar für Unabhängigkeits- und Goodness-of-Fit-Tests
- **Monte-Carlo**: π-Schätzung als Einstiegsbeispiel, skalierbar auf beliebige Probleme
- **Normalverteilung PPF**: Numerische Umkehrung der CDF via Bisektionsverfahren

### Build 4: Beweistheorie-Modul (2026-03-08)
- `src/proof_theory.py` implementiert: Collatz, Goldbach, Zwillingsprimzahlen,
  Riemann-Zeta, Sieb, Miller-Rabin, Legendre/Jacobi, CRT, ProofByInduction
- 42 neue Tests, alle grün (Gesamt: 176/176)
- `research/OPEN_CONJECTURES.md`: Alle Millennium-Probleme dokumentiert
- `research/PROOF_STRATEGIES.md`: Beweis-Strategien und Lernpfad

### FERNZIEL (ab Build 4)
Das Projekt richtet sich auf das **Beweisen oder Widerlegen offener Vermutungen**:
1. **Priorität 1:** Collatz (empirisch, Musteranalyse, Tao-Ansatz)
2. **Priorität 2:** Goldbach (Sieb-Methoden, Kreismethode)
3. **Priorität 3:** Riemann-Hypothese (analytische Zahlentheorie)
4. **Langfristig:** P vs NP, Navier-Stokes

### Geplante Erweiterungen (Priorität)
1. **Analytische Zahlentheorie** (L-Funktionen, Explizite Formeln, Siebmethoden)
2. **Komplexe Analysis** (vollständige ζ-Funktion, Residuensatz)
3. **Visualisierungs-Modul** (matplotlib 2D/3D, Riemann-Flächen)
4. **Optimierung** (Gradient Descent, BFGS, Simplex-Algorithmus)

## Brainstorming: Verbindungen

### "Fraktal" → Projekt
Mandelbrot-Menge: z_{n+1} = z_n² + c, Konvergenz-Check mit |z| < 2.
Benötigt: komplexe Zahlen (algebra.py), numpy-Vektorisierung für Geschwindigkeit.
Visualisierung mit matplotlib (geplant).

### "Fourier" → Projekt
Fourier-Transformation zerlegt Signale in Frequenzanteile.
Mathematisch: orthogonale Funktionen (analog zu Gram-Schmidt für Vektoren).
Implementierung: scipy.fft oder eigene DFT via Komplexe Analysis.
Anwendung: Signalverarbeitung, Differentialgleichungen lösen (Frequenzraum).

### "Primzahlen" → Kryptographie
RSA basiert auf: Schwierigkeit der Primfaktorzerlegung großer Zahlen.
Bereits implementiert: prime_factorization, euler_phi (Grundlagen).
Nächster Schritt: Miller-Rabin-Test für kryptographisch sichere Primzahlen.

### "Topologie" → Fixpunktsätze
Banachscher Fixpunktsatz: Kontraktionsabbildung hat genau einen Fixpunkt.
→ Begründet Newton-Raphson-Konvergenz (bereits implementiert).
Brouwerscher Fixpunktsatz: Jede stetige Selbstabbildung auf kompaktem Raum hat Fixpunkt.
→ Anwendung: Spieltheorie, Nash-Gleichgewichte.

### "Optimierung" → Maschinelles Lernen
Gradient Descent minimiert Verlustfunktionen in ML.
Mathematisch: Gradientenberechnung (Ableitungen, bereits implementiert).
Erweiterung: BFGS (quasi-Newton), L-BFGS-B (gedächtnissparend).

### "Laplace" → Regelungstechnik
Laplace-Transformation wandelt ODE in algebraische Gleichungen um.
Bereits implementiert: numerische Laplace + inverse Laplace (Stehfest).
Anwendung: Übertragungsfunktionen, Stabilitätsanalyse (Polstellen).

## Quellen & Recherche
- Numerik der DGL: Hairer, Nørsett, Wanner "Solving ODEs"
- Lineare Algebra: Gilbert Strang "Introduction to Linear Algebra"
- Zahlentheorie: Hardy & Wright "An Introduction to the Theory of Numbers"
- Numerische Analysis: Stoer & Bulirsch "Numerische Mathematik"
- Statistik: Casella & Berger "Statistical Inference"
- Fourier-Analysis: Körner "Fourier Analysis"
