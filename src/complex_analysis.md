# complex_analysis.py – Komplexe Analysis

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/complex_analysis.py`

---

## Überblick

Dieses Modul implementiert die Kernobjekte der komplexen Analysis mit Fokus auf die **Riemann-Zeta-Funktion** und verwandte Funktionen:

| Funktion | Beschreibung |
|----------|-------------|
| `gamma_lanczos` | Gamma-Funktion Γ(z) für komplexe z |
| `log_gamma` | ln(Γ(z)) via Stirling-Reihe |
| `riemann_zeta` | Vollständige ζ(s) für beliebige komplexe s ≠ 1 |
| `xi_function` | Riemanns ξ-Funktion (symmetrische Form) |
| `xi_symmetry_check` | Verifikation ξ(s) = ξ(1-s) |
| `zeta_on_critical_line` | ζ(1/2 + it) auf der kritischen Geraden |
| `riemann_siegel_z` | Riemann-Siegel-Z-Funktion Z(t) |
| `find_zeta_zeros` | Nullstellensuche via Vorzeichenwechsel |
| `N_count_formula` | N(T): Anzahl der Nullstellen bis T |
| `functional_equation_chi` | Faktor χ(s) der Funktionalgleichung |
| `verify_functional_equation` | Numerische Verifikation ζ(s) = χ(s)·ζ(1-s) |
| `residue_at_pole` | Residuum aus Laurent-Koeffizienten |
| `cauchy_integral_numerical` | Numerisches Cauchy-Integral |

---

## Gamma-Funktion

### `gamma_lanczos(z: complex) → complex`

Berechnet **Γ(z)** für komplexe Zahlen via Lanczos-Approximation.

**Die Gamma-Funktion** verallgemeinert die Fakultät auf komplexe Zahlen:
```
Γ(n) = (n-1)!    für positive ganze n
Γ(z+1) = z·Γ(z) (Funktionalgleichung)
Γ(1/2) = √π
```

**Pole:** Γ hat Pole bei z = 0, -1, -2, -3, ...

**Lanczos-Approximation (1964):**
Für `Re(z) > 0.5`:
```
Γ(z) ≈ √(2π) · t^{z+0.5} · e^{-t} · Σ cᵢ/(z+i)    mit t = z + g + 0.5
```
(g=7, 9 Koeffizienten nach Paul Godfrey)

**Reflexionsformel** für `Re(z) < 0.5`:
```
Γ(z)·Γ(1-z) = π / sin(πz)
```

**Genauigkeit:** ~15 signifikante Stellen.

### `log_gamma(z: complex) → complex`

Berechnet **ln(Γ(z))** für `Re(z) > 0` via Stirling-Reihe.

**Stirling-Näherung** (asymptotisch für |z| → ∞):
```
ln Γ(z) ≈ (z-½)·ln(z) - z + ½·ln(2π) + 1/(12z) - 1/(360z³) + ...
```

Numerisch stabiler als `gamma_lanczos` für sehr große |z| (vermeidet Overflow).

---

## Vollständige Riemann-Zeta-Funktion

### `riemann_zeta(s: complex, precision=200) → complex`

Berechnet **ζ(s) für beliebige komplexe s ≠ 1** via dreiteiliger Strategie:

**Definition:**
```
ζ(s) = Σ_{n=1}^∞ n^{-s}    (konvergiert nur für Re(s) > 1)
```

### Dreifallstrategie

| Region | Methode |
|--------|---------|
| Re(s) > 1 | Euler-Maclaurin (direkte Konvergenz) |
| 0 < Re(s) ≤ 1 | Dirichlet-Eta + Euler-Knopp-Beschleunigung |
| Re(s) ≤ 0 | Funktionalgleichung (Spiegelung nach Re ≥ 1) |

### Fall 2: Euler-Knopp-Beschleunigung

**Problem:** Die Eta-Reihe `η(s) = Σ (-1)^{n-1}/n^s` konvergiert für Re(s) > 0, aber sehr langsam auf der kritischen Geraden (O(1/√N) Konvergenz → ~10^6 Terme für 3 Stellen).

**Lösung (Euler-E-Transformation):**
```
η(s) ≈ Σ_{k=0}^{n-1} (-1)^k · (1/2^{k+1}) · Δ^k c₀
```
mit Vorwärtsdifferenzen `Δ^k c_j = c_{j+k} - k·c_{j+k-1} + ...` und **unsigned** Termen `cⱼ = 1/(j+1)^s`.

Mit n=60 Termen: Fehler ~2^{-60} ≈ 10^{-18} (Maschinengenauigkeit)!

### Fall 3: Funktionalgleichung

**Riemanns Funktionalgleichung (1859):**
```
ζ(s) = 2^s · π^{s-1} · sin(πs/2) · Γ(1-s) · ζ(1-s)
```

Für Re(s) ≤ 0 ist Re(1-s) ≥ 1 → Fall 1 anwendbar.

### Interne Hilfsfunktion: `_eta_euler_accelerated(s, n=60)`

Berechnet η(s) via Euler-Knopp. Die korrekte Implementierung verwendet **unsigned Terme** (ohne Vorzeichen):
```
c[j] = 1/(j+1)^s    (NICHT (-1)^j · c[j])
```
Das Vorzeichen wird erst bei der Gewichtung `(-1)^k / 2^{k+1}` eingeführt.

### Interne Hilfsfunktion: `_zeta_euler_maclaurin(s, terms=200)`

Direkte Summe + Euler-Maclaurin-Korrektur:
```
ζ(s) ≈ Σ_{n=1}^{N} n^{-s} + N^{1-s}/(s-1) + N^{-s}/2
        + (1/12)·(-s)·N^{-s-1} + ...
```

---

## Riemanns ξ-Funktion

### `xi_function(s: complex) → complex`

**Definition:**
```
ξ(s) = ½ · s · (s-1) · π^{-s/2} · Γ(s/2) · ζ(s)
```

**Eigenschaften:**
- `ξ(s) = ξ(1-s)` **(Symmetrie um Re(s) = 1/2)**
- ξ(s) ist ganz (keine Pole, obwohl ζ einen Pol bei s=1 hat)
- Nullstellen von ξ = nicht-triviale Nullstellen von ζ
- **Riemann-Hypothese:** Alle Nullstellen von ξ haben Re(s) = 1/2

Die ξ-Funktion macht die Symmetrie besonders deutlich: Da ξ(s) = ξ(1-s), kommen Nullstellen als **konjugierte Paare** bezüglich Re = 1/2.

### `xi_symmetry_check(s: complex) → dict`

Verifiziert die Symmetrie ξ(s) = ξ(1-s) numerisch. Wenn der Fehler < 1e-6 ist, ist die Implementierung korrekt.

---

## Analyse auf der kritischen Geraden

### `zeta_on_critical_line(t: float, precision=300) → complex`

Berechnet **ζ(1/2 + it)** auf der kritischen Geraden.

Die **Riemann-Hypothese** besagt: Alle nicht-trivialen Nullstellen haben genau diese Form (Re = 1/2).

### `riemann_siegel_z(t: float) → float`

**Riemann-Siegel-Z-Funktion** – reellwertig, Nullstellen = Nullstellen von ζ auf der kritischen Geraden:

```
Z(t) = e^{iθ(t)} · ζ(1/2 + it)   ∈ ℝ
```

mit dem **Riemann-Siegel-Winkel:**
```
θ(t) = Im(ln Γ(1/4 + it/2)) - (t/2)·ln(π)
```

**Warum Z(t)?** Vorzeichenwechsel von Z(t) zeigen Nullstellen an – einfacher zu detektieren als |ζ| nahe Null.

**Bekannte Nullstellen:** t ≈ 14.135, 21.022, 25.011, 30.425, 32.935, ...

### `find_zeta_zeros(t_min, t_max, steps=1000) → list[dict]`

**Nullstellensuche** via:
1. Berechne Z(t) auf feinem Gitter
2. Suche Vorzeichenwechsel
3. Verfeinere mit Bisektionsverfahren (50 Iterationen → ~2^{-50} Genauigkeit)

**Rückgabe** pro gefundener Nullstelle:
```python
{
    "t": 14.1347...,        # Imaginärteil
    "s": (0.5 + 14.1347j), # Vollständige komplexe Zahl
    "on_critical_line": True,
    "abs_zeta": 1.2e-8      # Wie nah ist |ζ| an 0?
}
```

### `N_count_formula(T: float) → float`

**Riemann-von-Mangoldt-Formel** für die Anzahl der Nullstellen:
```
N(T) = T/(2π) · ln(T/(2π)) - T/(2π) + 7/8 + O(ln T)
```

Konsistenzprüfung: Wenn wir numerisch N Nullstellen finden und die Formel ~N liefert, haben wir wahrscheinlich alle erfasst.

---

## Funktionalgleichung

### `functional_equation_chi(s: complex) → complex`

Berechnet den Faktor **χ(s)** der Funktionalgleichung `ζ(s) = χ(s)·ζ(1-s)`:
```
χ(s) = 2^s · π^{s-1} · sin(πs/2) · Γ(1-s)
```

**Eigenschaft auf der kritischen Geraden:** `|χ(1/2 + it)| = 1` (unitär).

### `verify_functional_equation(s: complex) → dict`

Numerische Verifikation der Funktionalgleichung. Fehler < 1e-4 bestätigt korrekte Implementierung.

---

## Residuensatz und Cauchy-Integral

### `residue_at_pole(f_coeffs, order) → complex`

Berechnet das **Residuum** einer meromorphen Funktion aus ihrer Laurent-Reihe.

**Laurent-Reihe** um z=0: `f(z) = Σ_{n=-m}^∞ aₙ·zⁿ`

Das Residuum ist der Koeffizient `a_{-1}`.

**Residuensatz:**
```
∮_C f(z) dz = 2πi · Σ Res(f, zₖ)
```

### `cauchy_integral_numerical(f, center, radius, n_points=1000) → complex`

**Numerisches Cauchy-Integral:**
```
f(z₀) = 1/(2πi) · ∮_C f(z)/(z-z₀) dz
```

**Implementierung:** Mittelwerteigenschaft holomorpher Funktionen:
```
f(z₀) = 1/(2π) · ∫_0^{2π} f(z₀ + r·e^{iθ}) dθ
       ≈ (1/n) · Σ_{k=0}^{n-1} f(z₀ + r·e^{i·2πk/n})
```

Die diskrete Trapezregel ist **exakt für trigonometrische Polynome** bis Grad n/2.

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `cmath` | Komplexe Exponential-/Winkelfunktionen |
| `math` | Reelle math. Funktionen |
| `numpy` | Linspace für Gitterpunkte |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from complex_analysis import (
    gamma_lanczos, riemann_zeta, xi_function,
    riemann_siegel_z, find_zeta_zeros, N_count_formula
)

# Gamma-Funktion
import math
print(gamma_lanczos(complex(0.5)))    # √π ≈ 1.7725
print(gamma_lanczos(complex(5)))      # 4! = 24.0

# Zeta-Funktion (drei Regionen)
print(riemann_zeta(complex(2)))       # π²/6 ≈ 1.6449
print(riemann_zeta(complex(0.5, 14))) # Nahe Nullstelle

# Xi-Symmetrie verifizieren
result = xi_function(complex(0.25, 5))
# xi_function(0.25+5j) ≈ xi_function(0.75-5j)

# Z-Funktion: Vorzeichenwechsel bei Nullstellen
z_vals = [(t, riemann_siegel_z(t)) for t in [13, 14, 15, 16]]
print(z_vals)  # Z wechselt Vorzeichen bei t ≈ 14.135

# Nullstellen suchen
zeros = find_zeta_zeros(t_min=10, t_max=35, steps=500)
print(f"Gefundene Nullstellen: {len(zeros)}")
for z in zeros:
    print(f"t = {z['t']:.6f}, |ζ| = {z['abs_zeta']:.2e}")

# Anzahl der Nullstellen schätzen
print(N_count_formula(100))  # ≈ 29.1 (tatsächlich: 29 Nullstellen bis T=100)
```
