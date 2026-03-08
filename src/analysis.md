# analysis.py – Analysis-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/analysis.py`

---

## Überblick

Dieses Modul implementiert numerische und symbolische Methoden der mathematischen Analysis:
- Numerische Differentiation (1. und 2. Ordnung)
- Numerische Integration (Simpson-Regel)
- Nullstellensuche (Newton-Raphson und Bisektionsverfahren)
- Taylor-Reihenentwicklung (symbolisch via SymPy)
- Stetigkeitstest (heuristisch)

---

## Funktion: `numerical_derivative`

```python
numerical_derivative(f, x, h=None, order=1) → float
```

Berechnet die numerische Ableitung einer Funktion im Punkt `x`.

### Methode: Zentraler Differenzenquotient

**1. Ableitung:**
```
f'(x) ≈ [f(x+h) - f(x-h)] / (2h)     Fehler: O(h²)
```

**2. Ableitung:**
```
f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h²     Fehler: O(h²)
```

### Herleitung (Taylor-Entwicklung)
```
f(x+h) = f(x) + h·f'(x) + (h²/2)·f''(x) + O(h³)
f(x-h) = f(x) - h·f'(x) + (h²/2)·f''(x) + O(h³)
Subtraktion: f(x+h) - f(x-h) = 2h·f'(x) + O(h³)
→ f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
```

Der zentrale Differenzenquotient ist **genauer als der einseitige** (O(h²) vs. O(h)).

### Optimale Schrittweiten

| Ordnung | Optimale h | Wert | Begründung |
|---------|-----------|------|------------|
| 1 | ε^(1/3) | ~6×10⁻⁶ | Gleichgewicht Auslöschung/Abbruch |
| 2 | ε^(1/4) | ~1.5×10⁻⁴ | h² im Nenner → mehr Auslöschung |

(ε = Maschinengenauigkeit ≈ 2.2×10⁻¹⁶)

---

## Funktion: `numerical_integral`

```python
numerical_integral(f, a, b, n=10000) → float
```

Berechnet `∫ₐᵇ f(x) dx` via **Simpson-Regel**.

### Formel

```
∫ₐᵇ f(x) dx ≈ (h/3) · [f(x₀) + 4f(x₁) + 2f(x₂) + 4f(x₃) + ... + 4f(x_{n-1}) + f(xₙ)]
```

Dabei ist `h = (b-a)/n` und `n` muss **gerade** sein.

### Gewichtungsschema

| Index | Gewicht |
|-------|---------|
| 0 (Startpunkt) | 1 |
| Ungerade Indizes | 4 |
| Gerade Indizes (außer Enden) | 2 |
| n (Endpunkt) | 1 |

### Genauigkeit

| Methode | Fehlerordnung |
|---------|--------------|
| Rechteckregel | O(h) |
| Trapezregel | O(h²) |
| **Simpson-Regel** | **O(h⁴)** |

Die Simpson-Regel approximiert das Integral durch **stückweise Parabeln** durch je 3 Punkte. Die exakte Integration einer Parabel liefert den Faktor h/3·(f(a) + 4f(m) + f(b)).

---

## Funktion: `newton_raphson`

```python
newton_raphson(f, x0, tol=1e-10, max_iter=1000, h=1e-7) → float
```

**Newton-Raphson-Verfahren** zur Nullstellensuche.

### Iterationsformel

```
x_{n+1} = x_n - f(x_n) / f'(x_n)
```

### Geometrische Interpretation
In jedem Schritt wird die **Tangente** an f im Punkt (x_n, f(x_n)) gezogen und deren Nullstelle als nächste Näherung verwendet.

### Konvergenzeigenschaften

| Eigenschaft | Beschreibung |
|-------------|-------------|
| Konvergenzrate | Quadratisch (Fehler wird quadriert) |
| Voraussetzung | f'(x) ≠ 0 nahe der Nullstelle |
| Garantie | Keine globale Konvergenz (abhängig von x₀) |

### Implementierungsdetail
Die Ableitung `f'(x)` wird **numerisch berechnet** (kein `df`-Parameter nötig).

### Abbruchbedingungen
- `|f(x)| < tol` → Nullstelle gefunden
- `|f'(x)| < 1e-15` → Divergenz (Tangente waagerecht)
- `max_iter` überschritten → `RuntimeError`

---

## Funktion: `bisection`

```python
bisection(f, a, b, tol=1e-10, max_iter=1000) → float
```

**Bisektionsverfahren** (Intervallhalbierungsverfahren) zur Nullstellensuche.

### Voraussetzung
`f(a)` und `f(b)` müssen **verschiedene Vorzeichen** haben (Zwischenwertsatz: dann existiert mindestens eine Nullstelle in [a,b]).

### Algorithmus
```
1. Berechne Mittelpunkt m = (a+b)/2
2. Wenn |f(m)| < tol: fertig
3. Wenn f(a)·f(m) < 0: Nullstelle in [a,m] → b = m
4. Sonst: Nullstelle in [m,b] → a = m
5. Wiederhole
```

### Fehlergarantie
Nach n Schritten: `|Fehler| ≤ (b-a) / 2ⁿ`

### Vergleich mit Newton-Raphson

| Kriterium | Newton-Raphson | Bisektionsverfahren |
|-----------|---------------|---------------------|
| Konvergenz | Quadratisch (schnell) | Linear (langsam) |
| Globalität | Nicht garantiert | Garantiert |
| Robustheit | Kann divergieren | Immer konvergent |

---

## Funktion: `taylor_series`

```python
taylor_series(f, center, degree, evaluate_at) → float
```

Berechnet die **Taylor-Reihe** einer Funktion symbolisch via SymPy.

### Taylor-Formel

```
f(x) = Σ_{k=0}^{n} [f^(k)(c) / k!] · (x-c)^k
```

### Wichtige Beispiele (Entwicklung um 0)

| Funktion | Taylor-Reihe |
|----------|-------------|
| eˣ | 1 + x + x²/2! + x³/3! + ... |
| sin(x) | x - x³/3! + x⁵/5! - ... |
| cos(x) | 1 - x²/2! + x⁴/4! - ... |
| ln(1+x) | x - x²/2 + x³/3 - ... (|x| < 1) |

### Warum symbolisch?
Numerische finite Differenzen für die k-te Ableitung haben h^k im Nenner – für **k > 8** dominieren Auslöschungsfehler exponentiell. SymPy berechnet dagegen exakte symbolische Ableitungen.

### Unterstützte Funktionen
`math.sin`, `math.cos`, `math.exp`, `math.log`, `math.tan`, `math.sqrt` sowie beliebige Lambdas mit SymPy-Kompatibilität.

---

## Funktion: `is_continuous`

```python
is_continuous(f, x, delta=1e-6) → bool
```

Heuristischer Stetigkeitstest: Überprüft ob `lim_{t→x} f(t) = f(x)`.

**Hinweis:** Dies ist ein **heuristischer Test** – er kann Unstetigkeitsstellen mit zu kleinem `delta` übersehen.

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `sympy` | Symbolische Differentiation für Taylor-Reihen |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from analysis import numerical_derivative, numerical_integral, newton_raphson, bisection

import math

# Ableitung von sin(x) bei x=π/4
f = math.sin
dfdx = numerical_derivative(f, math.pi/4)
print(dfdx)   # ≈ 0.7071 (= cos(π/4))

# Integral von x² von 0 bis 1
result = numerical_integral(lambda x: x**2, 0, 1)
print(result)   # ≈ 0.3333 (= 1/3)

# Nullstelle von x² - 2 = 0 (Wurzel aus 2)
root = newton_raphson(lambda x: x**2 - 2, x0=1.5)
print(root)   # ≈ 1.41421356 (= √2)

# Bisektionsverfahren
root2 = bisection(lambda x: x**2 - 2, a=1, b=2)
print(root2)  # ≈ 1.41421356

# Taylor-Reihe: sin(x) bis Grad 5 bei x=0.5
from analysis import taylor_series
approx = taylor_series(math.sin, center=0, degree=5, evaluate_at=0.5)
print(approx)  # ≈ 0.47943 (= sin(0.5))
```

---

## Numerische Stabilitätshinweise

- **Zentraler Differenzenquotient** ist dem einseitigen immer vorzuziehen (doppelte Genauigkeit bei gleichem Aufwand)
- **Schrittweite h zu klein:** Auslöschungsfehler in `f(x+h) - f(x-h)` dominieren
- **Schrittweite h zu groß:** Abbruchfehler der Differenzenapproximation dominieren
- **Simpson vs. Trapez:** Simpson ist mit O(h⁴) deutlich besser, aber erfordert gerades n
