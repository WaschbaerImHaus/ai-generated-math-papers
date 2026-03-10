# Modul: Differentialgeometrie (`differential_geometry.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `differential_geometry.py` implementiert die klassische Differentialgeometrie für parametrisierte Kurven und Flächen im $\mathbb{R}^3$. Es umfasst die Frenet-Serret-Theorie, die erste und zweite Fundamentalform von Flächen, geodätische Berechnungen, Riemannsche Krümmungsgrößen sowie klassische Minimalflächen und ebene Kurventheorie.

**Hinweis:** Die Klassen `MetricTensor`, `christoffel_symbols`, `riemann_tensor` (tensorielle Version), `einstein_tensor` und `geodesic_equation` befinden sich in `tensor_geometry.py` und werden dort importiert. Dieses Modul ergänzt mit geometrisch-analytischer Perspektive.

---

## Mathematischer Hintergrund

### Frenet-Serret-Formeln

Für eine reguläre Raumkurve $\gamma(t)$ mit Einheitstangentenvektor $\mathbf{T}$, Hauptnormalvektor $\mathbf{N}$ und Binormalvektor $\mathbf{B}$:

$$\mathbf{T}' = \kappa \mathbf{N}, \quad \mathbf{N}' = -\kappa \mathbf{T} + \tau \mathbf{B}, \quad \mathbf{B}' = -\tau \mathbf{N}$$

wobei $\kappa$ die **Krümmung** und $\tau$ die **Torsion** bezeichnen.

**Krümmung einer parametrisierten Kurve:**
$$\kappa = \frac{|\gamma' \times \gamma''|}{|\gamma'|^3}$$

**Torsion einer Raumkurve:**
$$\tau = \frac{(\gamma' \times \gamma'') \cdot \gamma'''}{|\gamma' \times \gamma''|^2}$$

**Bogenlänge:**
$$s(t_0, t_1) = \int_{t_0}^{t_1} |\gamma'(t)|\, dt$$

### Erste Fundamentalform

Für eine parametrisierte Fläche $\mathbf{r}(u, v)$ mit $\mathbf{r}_u = \partial \mathbf{r}/\partial u$, $\mathbf{r}_v = \partial \mathbf{r}/\partial v$:

$$I = E\, du^2 + 2F\, du\, dv + G\, dv^2$$

$$E = \mathbf{r}_u \cdot \mathbf{r}_u, \quad F = \mathbf{r}_u \cdot \mathbf{r}_v, \quad G = \mathbf{r}_v \cdot \mathbf{r}_v$$

### Zweite Fundamentalform

Mit dem Einheitsnormalenvektor $\hat{\mathbf{n}} = (\mathbf{r}_u \times \mathbf{r}_v)/|\mathbf{r}_u \times \mathbf{r}_v|$:

$$II = L\, du^2 + 2M\, du\, dv + N\, dv^2$$

$$L = \mathbf{r}_{uu} \cdot \hat{\mathbf{n}}, \quad M = \mathbf{r}_{uv} \cdot \hat{\mathbf{n}}, \quad N = \mathbf{r}_{vv} \cdot \hat{\mathbf{n}}$$

### Gaußsche und Mittlere Krümmung

**Gaußsche Krümmung (Theorema Egregium):**
$$K = \frac{LN - M^2}{EG - F^2}$$

**Mittlere Krümmung:**
$$H = \frac{EN + GL - 2FM}{2(EG - F^2)}$$

**Hauptkrümmungen** $\kappa_1, \kappa_2$:
$$K = \kappa_1 \kappa_2, \quad H = \frac{\kappa_1 + \kappa_2}{2}$$

### Gauß-Bonnet-Theorem

$$\iint_M K\, dA + \oint_{\partial M} \kappa_g\, ds = 2\pi \chi(M)$$

Für eine geschlossene Fläche ohne Rand ($\partial M = \emptyset$):
$$\iint_M K\, dA = 2\pi \chi(M)$$

### Isoperimetrische Ungleichung (ebene Kurven)

$$L^2 \geq 4\pi A$$

Gleichheit genau dann, wenn die Kurve ein Kreis ist.

---

## Klassen und Methoden

### `ParametricCurve`

Parametrisierte Kurven $\gamma: [a,b] \to \mathbb{R}^n$ mit vollständiger Frenet-Serret-Theorie.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `velocity` | `(t: float) → np.ndarray` | Tangentialvektor $\gamma'(t)$ |
| `speed` | `(t: float) → float` | Geschwindigkeit $|\gamma'(t)|$ |
| `arc_length` | `(t0, t1: float) → float` | Bogenlänge $\int_{t_0}^{t_1} |\gamma'|\, dt$ |
| `curvature` | `(t: float) → float` | Krümmung $\kappa(t)$ |
| `torsion` | `(t: float) → float` | Torsion $\tau(t)$ (nur $\mathbb{R}^3$) |
| `unit_tangent` | `(t: float) → np.ndarray` | Einheitstangente $\mathbf{T}(t)$ |
| `principal_normal` | `(t: float) → np.ndarray` | Hauptnormalvektor $\mathbf{N}(t)$ |
| `binormal` | `(t: float) → np.ndarray` | Binormalvektor $\mathbf{B}(t)$ |
| `frenet_serret_frame` | `(t: float) → dict` | Frenet-Rahmen $\{T, N, B\}$ als dict |
| `frenet_serret_formulas` | `() → dict` | Symbolische Frenet-Gleichungen |

### `ParametricSurface`

Parametrisierte Flächen $\mathbf{r}(u, v): U \to \mathbb{R}^3$ mit Krümmungstheorie.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `partial_derivatives` | `(u, v: float) → tuple` | $\mathbf{r}_u$, $\mathbf{r}_v$ |
| `normal_vector` | `(u, v: float) → np.ndarray` | Normalenvektor $\mathbf{r}_u \times \mathbf{r}_v$ |
| `unit_normal` | `(u, v: float) → np.ndarray` | Einheitsnormale $\hat{\mathbf{n}}$ |
| `first_fundamental_form` | `(u, v: float) → (E, F, G)` | Koeffizienten der 1. Fundamentalform |
| `second_fundamental_form` | `(u, v: float) → (L, M, N)` | Koeffizienten der 2. Fundamentalform |
| `gaussian_curvature` | `(u, v: float) → float` | Gaußsche Krümmung $K$ |
| `mean_curvature` | `(u, v: float) → float` | Mittlere Krümmung $H$ |
| `principal_curvatures` | `(u, v: float) → (κ₁, κ₂)` | Hauptkrümmungen |
| `area_element` | `(u, v: float) → float` | Flächenelement $|\mathbf{r}_u \times \mathbf{r}_v|$ |
| `total_area` | `(u0, u1, v0, v1, ...) → float` | Gesamtfläche per numerischer Integration |

### `GeodesicComputation`

Geodätische Kurven auf Flächen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `geodesic_sphere` | `(radius, p0, v0, t_end, ...) → dict` | Geodätische auf der Sphäre (Großkreis) |
| `geodesic_equations_general` | `(metric_fn, x0, v0, t_span) → dict` | Allgemeine Geodätengleichung per ODE |
| `conjugate_points_demo` | `(surface_fn) → dict` | Demo: Konjugierte Punkte |
| `exponential_map_surface` | `(surface_fn, p, tangent_vectors, ...) → dict` | Exponentialabbildung auf Fläche |

### `RiemannianGeometry`

Riemannsche Krümmungsgrößen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `sectional_curvature` | `(metric_fn, x, u, v, ...) → float` | Schnittkrümmung $K(u, v)$ |
| `ricci_curvature` | `(metric_fn, x, ...) → np.ndarray` | Ricci-Tensor $R_{ij}$ |
| `scalar_curvature` | `(metric_fn, x, ...) → float` | Skalarkrümmung $R$ |
| `covariant_derivative_demo` | `(surface_fn, ...) → dict` | Kovariante Ableitung (Demo) |
| `parallel_transport_demo` | `(surface_fn, ...) → dict` | Paralleltransport auf Fläche |
| `gauss_bonnet_theorem_demo` | `(surface_type) → dict` | Numerische Gauß-Bonnet-Verifikation |

### `ClassicalSurfaces`

Fabrik für klassische Flächen als `ParametricSurface`-Objekte.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `sphere` | `(radius=1.0) → ParametricSurface` | Sphäre $S^2$ |
| `torus` | `(R=2.0, r=1.0) → ParametricSurface` | Torus mit Haupt- und Nebenradius |
| `cylinder` | `(radius=1.0, height=2.0) → ParametricSurface` | Zylinder |
| `helicoid` | `() → ParametricSurface` | Helikoid (minimale Wendelfläche) |
| `catenoid` | `() → ParametricSurface` | Katenoid (Minimalfläche aus Kettenlinie) |
| `saddle_surface` | `() → ParametricSurface` | Sattelfläche $z = x^2 - y^2$ |
| `minimal_surface_check` | `(surface, ...) → dict` | Prüft ob $H = 0$ (Minimalfläche) |

### `CurveTheory2D`

Ebene Kurventheorie – Evoluten, Evolventen, isoperimetrische Ungleichung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `curvature_from_cartesian` | `(f, x, ...) → float` | Krümmung von $y = f(x)$: $\kappa = y''/(1+y'^2)^{3/2}$ |
| `evolute` | `(curve, t_vals) → np.ndarray` | Evolute (Mittelpunkt der Schmiegkreise) |
| `involute` | `(curve, t_vals, ...) → np.ndarray` | Evolvente einer ebenen Kurve |
| `isoperimetric_inequality_check` | `(curve, t0, t1) → dict` | Prüft $L^2 \geq 4\pi A$ |
| `four_vertex_theorem_demo` | `() → dict` | Demo: Vier-Scheitelpunkte-Satz |
| `total_curvature_closed_curve` | `(curve, t0, t1) → float` | Gesamtkrümmung $\oint \kappa\, ds$ |

---

## Klassische Flächen – Übersicht

| Fläche | Gaußsche Krümmung $K$ | Mittlere Krümmung $H$ | Typ |
|--------|----------------------|----------------------|-----|
| Sphäre $S^2(r)$ | $K = 1/r^2 > 0$ | $H = 1/r$ | elliptisch |
| Zylinder | $K = 0$ | $H = 1/(2r)$ | parabolisch |
| Helikoid | $K < 0$ | $H = 0$ | Minimalfläche |
| Katenoid | $K < 0$ | $H = 0$ | Minimalfläche |
| Sattelfläche | $K < 0$ | $H = 0$ (am Ursprung) | hyperbolisch |
| Torus | $K$ wechselt Vorzeichen | variiert | gemischt |

---

## Beispiele

```python
import numpy as np
from differential_geometry import ParametricCurve, ParametricSurface, ClassicalSurfaces, CurveTheory2D

# -- Raumkurve: Helix --
helix = ParametricCurve(
    fn=lambda t: np.array([np.cos(t), np.sin(t), t]),
    dim=3
)
print(helix.curvature(0.0))       # κ = 0.5
print(helix.torsion(0.0))         # τ = 0.5
frame = helix.frenet_serret_frame(0.0)
print(frame["T"], frame["N"], frame["B"])

# -- Parametrisierte Fläche: Sphäre --
factory = ClassicalSurfaces()
sphere = factory.sphere(radius=1.0)
K = sphere.gaussian_curvature(1.0, 1.0)  # K ≈ 1.0 überall
H = sphere.mean_curvature(1.0, 1.0)       # H ≈ 1.0

# -- Ebene Kurventheorie: Kreis --
ct = CurveTheory2D()
circle = ParametricCurve(fn=lambda t: np.array([np.cos(t), np.sin(t)]), dim=2)
result = ct.isoperimetric_inequality_check(circle, 0, 2 * np.pi)
print(result["L_squared"], result["4piA"])  # beide ≈ 4π²
```

---

## Tests

**Testdatei:** `tests/test_differential_geometry.py`
**Abdeckung:** Alle Klassen – Frenet-Rahmen, Fundamentalformen, Krümmungen, Geodätische, Gauß-Bonnet-Verifikation, Minimalflächen, isoperimetrische Ungleichung

---

## Implementierungshinweise

- **Numerische Ableitungen:** Kurven und Flächen werden mit zentralen Differenzenquotienten ($h = 10^{-5}$) differenziert.
- **Torsion:** Wird nur für Kurven im $\mathbb{R}^3$ sinnvoll berechnet. Für $n \neq 3$ wird `0.0` zurückgegeben.
- **Geodäten:** Allgemeine Geodäten werden als ODE-System $\ddot{\gamma}^k + \Gamma^k_{ij}\dot{\gamma}^i\dot{\gamma}^j = 0$ mit `scipy.integrate.solve_ivp` gelöst.
- **Gauß-Bonnet:** Die numerische Verifikation integriert $K\, dA$ über das Parametergebiet und vergleicht mit $2\pi\chi$.
- **Katenoid vs. Helikoid:** Beide sind isometrisch – sie haben dieselbe 1. Fundamentalform, aber unterschiedliche 2. Fundamentalformen.
