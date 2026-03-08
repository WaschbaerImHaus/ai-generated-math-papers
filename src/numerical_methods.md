# numerical_methods.py – Numerische Methoden

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/numerical_methods.py`

---

## Überblick

Dieses Modul implementiert grundlegende numerische Verfahren aus drei Bereichen:

| Bereich | Inhalt |
|---------|--------|
| **Interpolation** | Lagrange, Newton, Kubischer Spline |
| **Optimierung** | Gradient Descent, Goldener Schnitt, Numerischer Gradient |
| **Lineare Programmierung** | Simplex-Algorithmus (Dantzig) |

---

## Interpolation

**Ziel:** Gesucht ist eine Funktion f, die durch N gegebene Punkte (xᵢ, yᵢ) geht.

### `lagrange_interpolation(x_vals, y_vals, x) → float`

**Lagrange-Interpolationspolynom** vom Grad ≤ n durch (n+1) Punkte:

```
Pₙ(x) = Σᵢ yᵢ · Lᵢ(x)
```

**Basis-Polynome:**
```
Lᵢ(x) = Π_{j≠i} (x - xⱼ) / (xᵢ - xⱼ)
```

**Eigenschaften:**
- `Pₙ(xᵢ) = yᵢ` (geht exakt durch alle Stützpunkte)
- Eindeutig bestimmt (für n+1 Punkte gibt es genau ein Polynom vom Grad ≤ n)

**Nachteil:** **Runge-Phänomen** – bei äquidistanten Knoten und hohem Grad kann das Polynom stark zwischen den Stützstellen schwingen.

---

### Klasse: `NewtonInterpolation`

**Newton-Interpolation mit dividierten Differenzen** – eine numerisch stabilere Alternative zu Lagrange.

**Newton-Darstellung:**
```
P(x) = f[x₀] + f[x₀,x₁]·(x-x₀) + f[x₀,x₁,x₂]·(x-x₀)(x-x₁) + ...
```

**Dividierte Differenzen:**
```
f[xᵢ] = yᵢ
f[xᵢ,...,xₖ] = (f[xᵢ₊₁,...,xₖ] - f[xᵢ,...,xₖ₋₁]) / (xₖ - xᵢ)
```

#### Konstruktor: `__init__(x_vals, y_vals)`
Berechnet alle dividierten Differenzen und speichert die Newton-Koeffizienten.

#### `evaluate(x) → float`
Auswertung via **Horner-Schema**:
```
P(x) = c₀ + (x-x₀)[c₁ + (x-x₁)[c₂ + ... + (x-xₙ₋₁)cₙ]...]
```

#### `add_point(x_new, y_new)`
Fügt einen neuen Stützpunkt hinzu – O(n) Update statt O(n²) Neuberechnung.

**Vorteile gegenüber Lagrange:**
| Kriterium | Lagrange | Newton |
|-----------|---------|--------|
| Neue Punkte hinzufügen | O(n²) Neuberechnung | O(n) Update |
| Numerische Stabilität | Gut | Besser |

---

### Klasse: `CubicSpline`

**Kubische Spline-Interpolation (natürlicher Spline):**

Ein kubischer Spline ist ein stückweise kubisches Polynom:
```
S(x) = aᵢ + bᵢ(x-xᵢ) + cᵢ(x-xᵢ)² + dᵢ(x-xᵢ)³   auf [xᵢ, xᵢ₊₁]
```

**Eigenschaften:**
- Geht durch alle Stützpunkte (xᵢ, yᵢ)
- Stetig differenzierbar bis zur 2. Ordnung (C²)
- **Natürliche Randbedingung:** S''(x₀) = S''(xₙ) = 0

**Physikalische Bedeutung:** Minimiert die **Biegeenergie** ∫(S'')² dx → "natürlichste" Kurve durch die Punkte.

**Kein Runge-Phänomen!** Gleichmäßige Konvergenz garantiert.

#### Konstruktor: `__init__(x_vals, y_vals)`
Löst ein **tridiagonales LGS** (Thomas-Algorithmus) für die Koeffizienten cᵢ.

**Schrittweiten:** `hᵢ = xᵢ₊₁ - xᵢ`

**Koeffizientenbeziehungen:**
```
bᵢ = (yᵢ₊₁ - yᵢ)/hᵢ - hᵢ(2cᵢ + cᵢ₊₁)/3
dᵢ = (cᵢ₊₁ - cᵢ) / (3hᵢ)
```

#### `evaluate(x) → float`
Auswertung via **binäre Suche** für das richtige Teilintervall, dann Horner-Schema:
```
S(x) = aᵢ + dx·(bᵢ + dx·(cᵢ + dx·dᵢ))    mit dx = x - xᵢ
```

#### `derivative(x) → float`
Erste Ableitung des Splines:
```
S'(x) = bᵢ + 2cᵢ·dx + 3dᵢ·dx²
```

---

## Optimierung

### `gradient_descent(f, grad, x0, learning_rate=0.01, max_iter=1000, tol=1e-8) → (x_opt, f_opt, n_iter)`

**Gradient-Descent** zur Minimierung einer Funktion f: ℝⁿ → ℝ:

```
xₖ₊₁ = xₖ - α · ∇f(xₖ)
```

**Konvergenz:**
- Konvexe Funktionen: O(1/k)
- Stark konvexe Funktionen: Linear

**Wichtig:** Lernrate `α < 2/L` (L = Lipschitz-Konstante) für Konvergenz.

**Abbruchbedingung:** `||∇f(xₖ)|| < tol`

---

### `golden_section_search(f, a, b, tol=1e-8) → (x_opt, f_opt)`

**Goldener-Schnitt-Suche** für **unimodale** Funktionen auf [a, b].

**Voraussetzung:** f hat genau ein Minimum auf [a, b].

**Goldenes Verhältnis:**
```
φ = (1 + √5) / 2 ≈ 1.618
1/φ ≈ 0.618
```

**Algorithmus:**
```
x₁ = b - (1/φ)·(b-a)   (38.2% von links)
x₂ = a + (1/φ)·(b-a)   (61.8% von links)

Falls f(x₁) < f(x₂): Minimum in [a, x₂] → b = x₂, x₂ = x₁
Sonst:                Minimum in [x₁, b] → a = x₁, x₁ = x₂
```

**Pro Iteration:** Das Intervall wird auf ~61.8% verkleinert.
**Konvergenz:** Linear mit Rate 1/φ ≈ 0.618.

---

### `numerical_gradient(f, x, h=1e-5) → list`

**Numerischer Gradient via zentraler Differenzen:**
```
∂f/∂xᵢ ≈ [f(x + h·eᵢ) - f(x - h·eᵢ)] / (2h)
```

Nützlich wenn der analytische Gradient nicht bekannt ist.

---

## Lineare Programmierung (Simplex)

### `simplex(c, A, b) → (x_opt, f_opt)`

**Simplex-Algorithmus** (Dantzig 1947) für lineare Programmierung.

**Minimierungsproblem in Standardform:**
```
min  cᵀx
s.t. Ax ≤ b,   x ≥ 0,   b ≥ 0
```

### Algorithmus-Schritte

1. **Schlupfvariablen einführen:** `Ax + s = b`, Basis = Schlupfvariablen
2. **Tableau aufstellen:** `[A | I | b]` mit Zielfunktionszeile `[c | 0 | 0]`
3. **Pivotisierung:**
   - **Eintrittvariable:** Negativste Kosten in der Zielfunktionszeile (Blands Regel für endliche Terminierung)
   - **Austrittsvariable:** Quotientenregel (min bᵢ/Aᵢⱼ für Aᵢⱼ > 0)
4. **Basiswechsel:** Tableau durch Zeilenoperationen aktualisieren
5. Wiederholen bis alle Zielfunktionskoeffizienten ≥ 0

### Tableau-Struktur

```
          x₁  x₂  ...  xₙ  s₁  s₂  ...  sₘ  | b
Zeile 1:  [A₁₁ ...     Aₘₙ  1   0   ...  0   | b₁]
...
Zeile m:  [...                    0   ...  1   | bₘ]
Obj:      [c₁  ...     cₙ  0    0   ...  0   | 0 ]
```

### Konvergenz und Komplexität
- Im Worst-case exponentiell (aber in der Praxis sehr effizient)
- **Blands Regel** verhindert Zyklen (Wahl des kleinsten Indexes)

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from numerical_methods import (
    lagrange_interpolation, NewtonInterpolation, CubicSpline,
    gradient_descent, golden_section_search, simplex
)
import math

# Lagrange-Interpolation: sin-Stützpunkte
x_vals = [0, math.pi/6, math.pi/3, math.pi/2]
y_vals = [math.sin(x) for x in x_vals]
print(lagrange_interpolation(x_vals, y_vals, math.pi/4))  # ≈ sin(π/4) = 0.7071

# Newton-Interpolation
newton = NewtonInterpolation(x_vals, y_vals)
print(newton.evaluate(math.pi/4))   # Gleiche Ergebnis

# Kubischer Spline
spline = CubicSpline(x_vals, y_vals)
print(spline.evaluate(math.pi/4))   # Bessere Genauigkeit

# Gradient Descent: minimiere f(x,y) = x² + y²
f = lambda x: x[0]**2 + x[1]**2
grad_f = lambda x: [2*x[0], 2*x[1]]
x_opt, f_opt, iters = gradient_descent(f, grad_f, x0=[5.0, 3.0])
print(x_opt)   # ≈ [0, 0]

# Goldener Schnitt: minimiere (x-2)² auf [0, 5]
x_m, f_m = golden_section_search(lambda x: (x-2)**2, 0, 5)
print(x_m)    # ≈ 2.0

# Simplex: min 2x + y, s.t. x+y≤4, 2x+y≤6, x,y≥0
c = [2.0, 1.0]
A = [[1, 1], [2, 1]]
b = [4.0, 6.0]
x_lp, f_lp = simplex(c, A, b)
print(x_lp, f_lp)   # Optimale LP-Lösung
```
