# Modul: Klassische Geometrie (`classical_geometry.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10 | **Python:** 3.13.7 | **Abhängigkeiten:** numpy, math, itertools

---

## Überblick

Das Modul `classical_geometry.py` implementiert klassische, synthetische und nichteuklidische Geometrie. Es deckt sechs geometrische Rahmenwerke sowie sieben nützliche Standalone-Funktionen ab.

**Wichtiger Hinweis:** Dieses Modul ergänzt `tensor_geometry.py`. Die dort bereits vorhandenen Funktionen `hyperbolic_plane_metric`, `gaussian_curvature` und `riemann_tensor` werden hier **nicht** dupliziert.

---

## Klassen und Funktionen

### 1. `EuclideanGeometry`

Euklidische Geometrie in beliebiger Dimension (Standard: 2D), sowohl analytisch als auch synthetisch.

| Methode | Beschreibung | Formel (KaTeX) |
|---------|-------------|----------------|
| `distance(p1, p2)` | Euklidischer Abstand | $d = \sqrt{\sum (p_i - q_i)^2}$ |
| `angle_between(v1, v2)` | Winkel zwischen Vektoren | $\theta = \arccos\frac{u \cdot v}{|u||v|}$ |
| `are_parallel(l1_start, l1_end, l2_start, l2_end)` | Parallelität zweier Geraden | Kreuzprodukt = 0 |
| `are_perpendicular(v1, v2)` | Orthogonalität | $u \cdot v = 0$ |
| `circumcenter(p1, p2, p3)` | Umkreismittelpunkt (2D) | Schnittpunkt der Streckensymmetralen |
| `incenter(p1, p2, p3)` | Inkreismittelpunkt | $I = (a P_1 + b P_2 + c P_3)/(a+b+c)$ |
| `centroid(points)` | Schwerpunkt | $G = \frac{1}{n}\sum P_i$ |
| `euler_line_demo(p1, p2, p3)` | Euler-Gerade | $H = 3G - 2O$ |

**Euler-Gerade:** Verbindet Umkreismittelpunkt O, Schwerpunkt G und Höhenschnittpunkt H. Es gilt $|OG| : |GH| = 1 : 2$.

---

### 2. `ProjectiveGeometry`

Projektive Ebene $\mathbb{RP}^2$ mit homogenen Koordinaten $[x:y:w]$.

| Methode | Beschreibung |
|---------|-------------|
| `homogeneous_coords(point)` | $(x,y) \mapsto [x:y:1]$ |
| `projective_line_through(p1, p2)` | Gerade $\ell = P_1 \times P_2$ (Kreuzprodukt) |
| `projective_intersection(l1, l2)` | Schnittpunkt $P = \ell_1 \times \ell_2$ |
| `cross_ratio(p1, p2, p3, p4)` | Doppelverhältnis $(P_1,P_2;P_3,P_4)$ |
| `harmonic_conjugate(p1, p2, p3)` | Harmonische Konjugierte (CR = −1) |
| `projective_transformation(matrix, point)` | $P' = M \cdot P$ (Homographie) |
| `desargues_theorem_check(triangle1, triangle2)` | Desarguesscher Satz (numerisch) |
| `pappus_theorem_check(l1_points, l2_points)` | Satz von Pappus (numerisch) |

**Doppelverhältnis:**
$$\!(P_1, P_2; P_3, P_4) = \frac{(P_3 - P_1)(P_4 - P_2)}{(P_3 - P_2)(P_4 - P_1)}$$

Das Doppelverhältnis ist die fundamentale projektive Invariante – es bleibt unter jeder projektiven Transformation erhalten.

**Parallelität im Projektiven:** Zwei parallele Geraden schneiden sich in einem „Ferndpunkt" $[x:y:0]$ (dritte Komponente = 0).

---

### 3. `HyperbolicGeometry`

Hyperbolische Geometrie im **Poincaré-Kreismodell** $\{z \in \mathbb{C} : |z| < 1\}$ und im **oberen Halbebenenmodell** $\{z \in \mathbb{C} : \text{Im}(z) > 0\}$.

| Methode | Beschreibung |
|---------|-------------|
| `poincare_disk_distance(z1, z2)` | Hyperbolischer Abstand im Disk |
| `poincare_disk_midpoint(z1, z2)` | Hyperbolischer Mittelpunkt |
| `hyperbolic_angle_sum(triangle_vertices)` | Winkelsumme (< π) |
| `hyperbolic_area(triangle_vertices)` | Fläche = π − (α+β+γ) |
| `hyperbolic_parallel_lines()` | Demo: ∞ viele Parallelen |
| `gauss_bonnet_hyperbolic(polygon_angles)` | Gauß-Bonnet-Satz |
| `upper_half_plane_distance(z1, z2)` | Abstand im oberen Halbebenenmodell |

**Poincaré-Abstand:**
$$d(z_1, z_2) = 2 \operatorname{arctanh}\!\left(\frac{|z_1 - z_2|}{|1 - \bar{z}_1 z_2|}\right)$$

**Gauß-Bonnet für hyperbolische Polygone:**
$$A = (n-2)\pi - \sum_{i=1}^{n} \alpha_i$$

**Oberes Halbebenenmodell:**
$$d(z_1, z_2) = \operatorname{arccosh}\!\left(1 + \frac{|z_1 - z_2|^2}{2\,\text{Im}(z_1)\,\text{Im}(z_2)}\right)$$

**Winkelberechnung:** Da das Poincaré-Modell konform ist, werden Winkel korrekt durch Möbius-Transformation des Scheitelpunkts in den Ursprung berechnet. Geodäten durch den Ursprung sind Durchmesser, daher ist der euklidische Winkel dort identisch mit dem hyperbolischen.

---

### 4. `EllipticGeometry`

Elliptische (sphärische) Geometrie auf der Einheitssphäre $S^2$.

| Methode | Beschreibung |
|---------|-------------|
| `spherical_distance(lat1, lon1, lat2, lon2)` | Großkreisabstand (Haversine) |
| `spherical_angle_sum(triangle_vertices_rad)` | Winkelsumme (> π) |
| `spherical_excess(triangle_vertices_rad)` | Sphärischer Exzess E = α+β+γ−π |
| `spherical_area(triangle_vertices_rad)` | A = R²·E |
| `lune_area(dihedral_angle)` | A = 2R²·α |
| `girard_theorem(alpha, beta, gamma)` | A = R²·(α+β+γ−π) |

**Haversine-Formel:**
$$d = 2R\arcsin\!\sqrt{\sin^2\!\tfrac{\Delta\phi}{2} + \cos\phi_1\cos\phi_2\sin^2\!\tfrac{\Delta\lambda}{2}}$$

**Girard'scher Dreieckssatz (Albert Girard, 1625):**
$$A = R^2 \cdot (\alpha + \beta + \gamma - \pi)$$

---

### 5. `AffinGeometry`

Affine Geometrie – Geometrie ohne Abstands- und Winkelbegriff.

| Methode | Beschreibung |
|---------|-------------|
| `affine_combination(points, coefficients)` | $P = \sum \lambda_i P_i$, $\sum \lambda_i = 1$ |
| `affine_map(A, b, point)` | $f(x) = Ax + b$ |
| `barycentric_coords(p, p1, p2, p3)` | Baryzentrische Koordinaten |
| `affine_hull(points)` | Dimension der affinen Hülle |
| `convex_hull_2d(points)` | Konvexe Hülle (Gift-Wrapping / Jarvis-March) |
| `is_affinely_independent(points)` | Affine Unabhängigkeit (via Rangtest) |

**Baryzentrische Koordinaten:** $P = \lambda_1 P_1 + \lambda_2 P_2 + \lambda_3 P_3$ mit $\lambda_1 + \lambda_2 + \lambda_3 = 1$
- $P$ im Inneren ⟺ alle $\lambda_i > 0$
- $P$ auf einer Seite ⟺ ein $\lambda_i = 0$
- $P$ ist Eckpunkt ⟺ zwei $\lambda_i = 0$

**Gift-Wrapping-Algorithmus:** $O(n \cdot h)$ Laufzeit, wobei $h$ die Hüllenpunktanzahl ist.

---

### 6. `SyntheticGeometry`

Synthetische Geometrie basierend auf Hilberts Axiomensystem (1899).

| Methode | Beschreibung |
|---------|-------------|
| `hilbert_axioms_overview()` | Alle 5 Axiomengruppen I–V |
| `incidence_axioms_demo()` | Beispiele für Inzidenzaxiome |
| `parallel_axiom_demo()` | Vergleich der drei Geometrien |
| `axiom_independence_demo()` | Unabhängigkeitsbeispiele |
| `finite_projective_plane(q)` | PG(2,q) für Primzahl q |

**Hilberts Axiomengruppen:**
| Gruppe | Inhalt | Anzahl |
|--------|--------|--------|
| I | Inzidenz | 5 |
| II | Anordnung | 4 |
| III | Kongruenz | 6 |
| IV | Parallelen (Euklid) | 1 |
| V | Stetigkeit | 2 |

**Endliche projektive Ebene PG(2,q):**
$$\text{Punkte} = \text{Geraden} = q^2 + q + 1, \quad \text{Punkte/Gerade} = q + 1$$

Die Fano-Ebene PG(2,2) ist die kleinste projektive Ebene mit 7 Punkten und 7 Geraden.

---

### 7. Standalone-Funktionen

| Funktion | Beschreibung | Formel |
|----------|-------------|--------|
| `euler_characteristic_polygon(V, E, F)` | Euler-Charakteristik | $\chi = V - E + F$ |
| `pick_theorem(I, B)` | Picks Satz | $A = I + B/2 - 1$ |
| `ptolemy_theorem_check(p1, p2, p3, p4)` | Ptolemäus-Satz | $|AC||BD| = |AB||CD| + |AD||BC|$ |
| `nine_point_circle(p1, p2, p3)` | Neun-Punkte-Kreis | $N = (O+H)/2$, $r_9 = R/2$ |
| `morley_theorem_demo()` | Morley'scher Dreieckssatz | $s = 8R\sin(\alpha/3)\sin(\beta/3)\sin(\gamma/3)$ |

**Euler'sche Polyederformel:** Für konvexe Polyeder (und topologische Sphären): $\chi = V - E + F = 2$. Für den Torus: $\chi = 0$.

**Pick'scher Satz (Georg Pick, 1899):** Berechnet exakt den Flächeninhalt von Gittervielecken ohne Koordinatenintegration.

**Neun-Punkte-Kreis:** Geht durch 9 ausgezeichnete Punkte des Dreiecks:
1–3: Seitenmittelpunkte, 4–6: Höhenfußpunkte, 7–9: Mittelpunkte der Ecke-zu-H-Strecken.

---

## Implementierungshinweise

- **Numerische Stabilität:** Alle Kreuzprodukte und Determinanten verwenden eine Toleranz von $10^{-8}$ bis $10^{-14}$ je nach Kontext.
- **Hyperbolische Winkel:** Die korrekte Berechnung erfordert Möbius-Transformation des Scheitelpunkts in den Ursprung (nicht bloße euklidische Winkelberechnung).
- **Projektives Unendlich:** Punkte mit $w = 0$ in homogenen Koordinaten repräsentieren Fernpunkte (parallele Geraden treffen sich dort).
- **Sphärisches Kosinusgesetz:** Wird für die Winkelberechnung in der sphärischen Geometrie genutzt.

---

## Tests

**Testdatei:** `tests/test_classical_geometry.py`
**Anzahl Tests:** 122
**Abdeckung:** Alle Klassen und Funktionen, Edge Cases (entartete Dreiecke, Fernpunkte, Randpunkte des Poincaré-Disks)
