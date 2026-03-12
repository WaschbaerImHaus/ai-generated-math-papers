# lonely_runner.py — Dokumentation

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `lonely_runner.py`

---

## Überblick

Das Modul implementiert die Verifikation und Analyse der **Lonely-Runner-Vermutung** (Wills 1967, unabhängig Cusick 1973).

> **Vermutung (CONJECTURE — unbewiesen im Allgemeinen)**:
> Gegeben $n$ Läufer auf einem Einheitskreis (Umfang 1), die alle bei 0 starten und verschiedene ganzzahlige Geschwindigkeiten haben. Dann gibt es für jeden Läufer $i$ einen Zeitpunkt $t > 0$, zu dem er von **allen** anderen Läufern mindestens $\frac{1}{n}$ entfernt ist.

### Formale Aussage

Für $n \geq 2$ und paarweise verschiedene ganze Zahlen $v_0, \ldots, v_{n-1}$ (mit $v_0 = 0$ o.B.d.A.) existiert $t > 0$ mit:

$$\| t \cdot (v_i - v_j) \| \geq \frac{1}{n} \quad \text{für alle } j \neq i$$

wobei $\|x\| = \min(x \bmod 1, \; 1 - x \bmod 1)$ der **Kreisabstand** ist.

---

## Beweisstatus

| $n$ | Status | Quelle |
|-----|--------|--------|
| $2$ | Bewiesen (trivial) | — |
| $3, 4$ | Bewiesen | Cusick & Pomerance (1984) |
| $5$ | Bewiesen | Goddyn et al. (2006) |
| $6$ | Bewiesen | Bohman, Fonoberova & Pikhurko (2011) |
| $7$ | Bewiesen | Tao (2018) |
| $\geq 8$ | **OFFEN** (Conjecture) | — |

**Wichtig**: Die Vermutung ist für $n \geq 8$ unbewiesen. Dieses Modul implementiert Verifikation für $n \leq 6$.

---

## Klassen

### `LonelyRunnerConjecture`

Hauptklasse für Verifikation und Simulation einer Lonely-Runner-Instanz.

#### Initialisierung

```python
lrc = LonelyRunnerConjecture(velocities=[0, 1, 2, 3])
```

- `velocities`: Liste paarweise verschiedener ganzzahliger Geschwindigkeiten
- Normalisierung: $v_0$ wird auf 0 gesetzt (alle anderen relativ dazu)
- `threshold`: $\frac{1}{n}$ (Mindestabstand)

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `circle_distance(x)` | Kreisabstand $\|x\| = \min(x \bmod 1, 1 - x \bmod 1)$ |
| `circle_distance_fraction(x)` | Exakte Fraction-Version |
| `is_lonely_at_time(i, t)` | Ist Läufer $i$ zum Zeitpunkt $t$ lonely? |
| `is_lonely_at_time_exact(i, t)` | Exakte Version mit Fraction-Arithmetik |
| `find_lonely_time_simulation(i, t_max, steps)` | Numerische Suche |
| `find_lonely_time_exact(i, max_denom)` | Rationale Suche |
| `verify_conjecture_for_all_runners(t_max, steps)` | Alle Läufer prüfen |
| `min_distance_at_time(i, t)` | Minimaler Abstand bei Zeit $t$ |
| `describe_diophantine_reformulation()` | Textbeschreibung |
| `blocki_sketch_n4()` | Blocki-Beweisskizze für $n=4$ |

---

### `LonelyRunnerVerifier`

Exakter Verifikator für kleine Fälle ($n \leq 6$) mittels rationaler Arithmetik.

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `verify_n2(v1)` | Trivialfall $n=2$: liefert immer $(True, \frac{1}{2|v_1|})$ |
| `verify_runner_exact(velocities, runner_idx, max_denom)` | Exakte Verifikation eines Läufers |
| `verify_all_runners(velocities, max_denom)` | Alle Läufer einer Instanz prüfen |
| `generate_test_cases(n, max_vel)` | Test-Instanzen generieren |

---

## Mathematische Details

### Kreisabstand

Zwei Läufer mit Geschwindigkeiten $v_i$ und $v_j$ haben zum Zeitpunkt $t$ den Kreisabstand:

$$d_{ij}(t) = \|t \cdot (v_i - v_j)\| = \min\bigl(t(v_i - v_j) \bmod 1, \; 1 - t(v_i - v_j) \bmod 1\bigr)$$

### Diophantische Reformulierung

Die Vermutung ist äquivalent zur Aussage:

Für alle paarweise verschiedenen ganzen Zahlen $w_1, \ldots, w_{n-1}$ existiert $t \in \mathbb{R}$ mit:

$$\|t \cdot w_k\| \geq \frac{1}{n} \quad \text{für } k = 1, \ldots, n-1$$

### Algorithmus für `verify_runner_exact`

1. Berechne alle relativen Differenzen $d_j = |v_i - v_j|$ für Läufer $i$
2. Sammle kritische Zeitpunkte: $t = \frac{k}{d_j}$ für $k = 0, \ldots, d_j$
3. Bilde Mittelpunkte zwischen je zwei aufeinanderfolgenden kritischen Punkten
4. Prüfe zusätzlich alle Brüche $\frac{a}{b}$ mit $b \leq$ `max_denom`
5. Teste jeden Kandidaten mit `is_lonely_at_time_exact`

**Korrektheit**: Die Loneliness-Funktion ist zwischen kritischen Punkten konstant (stückweise konstant), daher reicht ein Punkt pro Intervall.

### Blocki-Beweis-Sketch für $n=4$

Blocki (2010) bewies den Fall $n=4$ mit folgendem Ansatz:

1. **Periodizität**: Mit $v_0=0, v_1=a, v_2=b, v_3=c$ und $\gcd(a,b,c)=1$ sind alle Positionen periodisch mit Periode $T = \text{lcm}(a,b,c)^{-1}$
2. **Kritische Intervalle**: Definiere $B_{ij} = \{t : \|t(v_i - v_j)\| < 1/4\}$ — jedes hat Maß $1/2$
3. **Überlappungslemma**: Das Maß von $[0,1) \setminus (B_{01} \cup B_{02} \cup B_{03})$ ist positiv (Cusick & Pomerance 1984)
4. **Symmetrie**: Gilt analog für alle Läufer

---

## Verwendungsbeispiele

```python
from lonely_runner import LonelyRunnerConjecture, LonelyRunnerVerifier
from fractions import Fraction

# n=2: Trivialfall
lrc = LonelyRunnerConjecture([0, 1])
print(lrc.is_lonely_at_time(0, 0.5))  # True (Abstand = 1/2 ≥ 1/2)

# n=3: Exakte Verifikation
lrc = LonelyRunnerConjecture([0, 1, 2])
t = Fraction(1, 3)
print(lrc.is_lonely_at_time_exact(0, t))  # True

# n=4: Simulation
lrc = LonelyRunnerConjecture([0, 1, 2, 3])
t = lrc.find_lonely_time_simulation(0, t_max=10.0, steps=100000)
print(t)  # z.B. 0.25

# Alle Läufer prüfen
verifier = LonelyRunnerVerifier()
results = verifier.verify_all_runners([0, 1, 2, 3])
for runner, (ok, t) in results.items():
    print(f"Läufer {runner}: lonely={'Ja' if ok else 'Nein'} bei t={t}")

# Diophantische Reformulierung beschreiben
lrc = LonelyRunnerConjecture([0, 1, 3, 5])
print(lrc.describe_diophantine_reformulation())

# Blocki-Beweis-Sketch für n=4
print(lrc.blocki_sketch_n4())
```

---

## Abhängigkeiten

- `numpy` (numerische Berechnungen)
- `sympy` (`Rational`, `gcd`, `lcm`, `factorint`)
- `fractions.Fraction` (exakte rationale Arithmetik)

---

## Referenzen

- Wills, J.M. (1967). Zwei Sätze über inhomogene diophantische Approximation.
- Cusick, T.W. (1973). View-obstruction problems in n-dimensional geometry.
- Cusick, T.W. & Pomerance, C. (1984). View-obstruction problems, III.
- Blocki, J. (2010). The Lonely Runner Problem (Master's Thesis).
- Bohman, T., Fonoberova, M. & Pikhurko, O. (2011). The Lonely Runner Problem for $n \leq 6$ Runners.
- Tao, T. (2018). Some remarks on the lonely runner conjecture.
