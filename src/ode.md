# ode.py – ODE-Modul (Gewöhnliche Differentialgleichungen)

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/ode.py`

---

## Überblick

Dieses Modul implementiert numerische Lösungsverfahren für gewöhnliche Differentialgleichungen (ODEs) sowie die Laplace-Transformation.

Ein **Anfangswertproblem (AWP)** hat die Form:
```
dy/dt = f(t, y),   y(t₀) = y₀
```

### Inhaltsübersicht

| Funktion | Methode | Fehlerordnung |
|----------|---------|--------------|
| `euler_method` | Explizites Euler-Verfahren | O(h) |
| `runge_kutta4` | Klassisches RK4 | O(h⁴) |
| `runge_kutta45` | Adaptives RK45 (Fehlberg) | O(h⁴/h⁵) |
| `solve_linear_ode` | Zustandsraum + RK4 | O(h⁴) |
| `laplace_transform` | Numerische Laplace-Transformation | – |
| `inverse_laplace` | Stehfest-Algorithmus | – |

---

## Typdefinitionen

```python
StateType = Union[float, List[float]]
ODEFunc   = Callable[[float, StateType], StateType]
```

Alle Verfahren unterstützen sowohl **skalare** ODEs als auch **Systeme** (Vektoren).

---

## Hilfsfunktionen

### `_add_states(a, b)` und `_scale_state(c, a)`
Abstraktion für Zustandsvektoren – funktioniert für Skalare und Listen:
```python
_add_states([1,2], [3,4])  # [4, 6]
_scale_state(2, [1,2])     # [2, 4]
```

---

## Euler-Verfahren

```python
euler_method(f, t0, y0, t_end, n_steps) → (t_list, y_list)
```

**Explizites Euler-Verfahren** – einfachstes ODE-Lösungsverfahren.

### Formel

```
y_{n+1} = y_n + h · f(t_n, y_n)
```

### Eigenschaften

| Eigenschaft | Wert |
|-------------|------|
| Fehlerordnung (lokal) | O(h²) |
| Fehlerordnung (global) | O(h) |
| Stabilität | Bedingt stabil (h muss klein genug) |
| Aufwand | 1 Funktionsauswertung pro Schritt |

**Verwendung:** Hauptsächlich für **didaktische Zwecke**. Für numerische Anwendungen sind RK4 oder RK45 deutlich besser.

---

## Runge-Kutta 4 (RK4)

```python
runge_kutta4(f, t0, y0, t_end, n_steps) → (t_list, y_list)
```

**Klassisches Runge-Kutta-Verfahren 4. Ordnung** – das Standard-ODE-Verfahren.

### Formel (Butcher-Tableau)

```
k₁ = h · f(t,       y)
k₂ = h · f(t + h/2, y + k₁/2)
k₃ = h · f(t + h/2, y + k₂/2)
k₄ = h · f(t + h,   y + k₃)

y_{n+1} = y_n + (k₁ + 2k₂ + 2k₃ + k₄) / 6
```

### Warum gewichteter Mittelwert?

| Stufenwert | Gewicht | Beschreibung |
|-----------|---------|-------------|
| k₁ | 1/6 | Steigung am Anfang des Intervalls |
| k₂ | 2/6 | Steigung in der Mitte (via k₁) |
| k₃ | 2/6 | Steigung in der Mitte (via k₂) |
| k₄ | 1/6 | Steigung am Ende des Intervalls |

Das Schema erinnert an die **Simpson-Regel**: Randpunkte zählen einmal, Mittelpunkte doppelt.

### Eigenschaften

| Eigenschaft | Wert |
|-------------|------|
| Fehlerordnung (lokal) | O(h⁵) |
| Fehlerordnung (global) | O(h⁴) |
| Aufwand | 4 Funktionsauswertungen pro Schritt |

---

## Runge-Kutta-Fehlberg (RK45, adaptiv)

```python
runge_kutta45(f, t0, y0, t_end, tol=1e-6, h_min=1e-10, h_max=0.1) → (t_list, y_list)
```

**Adaptives RK-Verfahren** mit automatischer Schrittweitenanpassung.

### Funktionsprinzip

Verwendet **eingebettete RK4/RK5-Methoden** für denselben Schritt:
- RK4-Lösung (y4) mit 4. Ordnung
- RK5-Lösung (y5) mit 5. Ordnung
- **Fehlerabschätzung:** `error = |y4 - y5|`

### Schrittweitenanpassung

```
h_neu = h · 0.9 · (tol / error)^(1/5)
```

Der Faktor 0.9 ist ein **Sicherheitsfaktor**, der übermäßig aggressive Vergrößerungen verhindert.

### Butcher-Tableau (Fehlberg-Koeffizienten)

Das Verfahren verwendet 6 Stufenwerte k₁...k₆. Die Knotenkoeffizienten (c-Werte) sind:
`c₂=1/4, c₃=3/8, c₄=12/13, c₅=1`

### Vorteile gegenüber RK4

| Kriterium | RK4 | RK45 |
|-----------|-----|------|
| Schrittweite | Fest | Automatisch |
| Effizienz bei glattem f | Gut | Sehr gut |
| Effizienz bei steilem f | Schlechte (klein h nötig) | Sehr gut |

---

## Lineare ODE mit konstanten Koeffizienten

```python
solve_linear_ode(coeffs, initial, t_vals, forcing=None) → List[float]
```

Löst eine ODE n-ter Ordnung:
```
a₀·y^(n) + a₁·y^(n-1) + ... + aₙ·y = g(t)
```

### Methode: Zustandsraum-Darstellung

Die ODE n-ter Ordnung wird in ein **System 1. Ordnung** umgewandelt:
```
z = [y, y', y'', ..., y^(n-1)]
z' = [z₁, z₂, ..., zₙ₋₁, (g(t) - Σ aₖ·zₖ) / a₀]
```

Dieses System wird mit RK4 gelöst.

---

## Laplace-Transformation

```python
laplace_transform(f, s, t_max=50.0, n_intervals=10000) → float
```

**Numerische einseitige Laplace-Transformation:**
```
F(s) = ∫₀^∞ f(t) · e^{-s·t} dt
```

Die Integration wird bei `t_max` abgebrochen (wenn `e^{-s·t_max}` vernachlässigbar klein).
Verwendet **Simpson-Regel** für die numerische Integration.

**Voraussetzung:** `s > 0` für Konvergenz.

### Bekannte Transformationspaare

| f(t) | F(s) |
|------|------|
| 1 | 1/s |
| t | 1/s² |
| e^{at} | 1/(s-a) |
| sin(ωt) | ω/(s²+ω²) |
| cos(ωt) | s/(s²+ω²) |

---

## Inverse Laplace-Transformation (Stehfest)

```python
inverse_laplace(F, t, sigma=1.0, n_terms=20) → float
```

**Numerische Inversion via Stehfest-Algorithmus:**
```
f(t) ≈ (ln2/t) · Σᵢ₌₁ᴺ Vᵢ · F(i·ln2/t)
```

**Stehfest-Koeffizienten Vᵢ:** Basieren auf Binomialkoeffizienten und alternierenden Vorzeichen.

**Eigenschaften:**
- Sehr effizient für **monotone Funktionen**
- `n_terms` gerade empfohlen
- Fehler nimmt mit n_terms zunächst ab, dann zu (Rundungsfehler)

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from ode import euler_method, runge_kutta4, runge_kutta45, solve_linear_ode

import math

# Beispiel: y' = -y, y(0) = 1  → analytische Lösung: y(t) = e^{-t}
f = lambda t, y: -y

# Euler-Verfahren
t_euler, y_euler = euler_method(f, t0=0, y0=1.0, t_end=5, n_steps=100)

# RK4
t_rk4, y_rk4 = runge_kutta4(f, t0=0, y0=1.0, t_end=5, n_steps=50)

# Vergleich bei t=1: y(1) = e^{-1} ≈ 0.3679
print(y_euler[20])   # Euler (weniger genau)
print(y_rk4[10])     # RK4 (sehr genau)

# Adaptives RK45
t_rk45, y_rk45 = runge_kutta45(f, t0=0, y0=1.0, t_end=5, tol=1e-8)

# System-ODE: harmonischer Oszillator y'' + ω²y = 0
# Als System: z = [y, y'], z' = [y', -ω²·y]
omega = 2.0
f_sys = lambda t, z: [z[1], -omega**2 * z[0]]
t_s, y_s = runge_kutta4(f_sys, t0=0, y0=[1.0, 0.0], t_end=2*math.pi, n_steps=1000)

# Lineare ODE: y'' + y = 0, y(0)=1, y'(0)=0
solution = solve_linear_ode(
    coeffs=[1, 0, 1],    # 1·y'' + 0·y' + 1·y = 0
    initial=[1.0, 0.0],  # y(0)=1, y'(0)=0
    t_vals=[0, 1, 2, 3]  # Auswertungspunkte
)
# Lösung: cos(t) → 1, 0.540, -0.416, -0.990
print(solution)
```
