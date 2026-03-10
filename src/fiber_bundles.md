# fiber_bundles.py – Dokumentation

## Übersicht

Das Modul `fiber_bundles.py` implementiert die grundlegenden mathematischen
Strukturen der Differentialgeometrie: Faserbündel, Verbindungen (Konnexionen),
Yang-Mills-Eichfelder und klassische Beispiele wie das Möbius-Bündel und die
Hopf-Faserung.

## Mathematischer Hintergrund

### Faserbündel

Ein Faserbündel ist ein Tupel $(E, B, \pi, F, G)$ mit:
- $E$ = Gesamtraum (total space)
- $B$ = Basisraum (base space)
- $\pi: E \to B$ = Projektion (surjektiv, stetig)
- $F$ = typische Faser (Homöomorphieklasse von $\pi^{-1}(b)$)
- $G$ = Strukturgruppe (wirkt auf $F$)

**Lokale Trivialisierung**: Für jede offene Menge $U \subset B$ existiert ein
Homöomorphismus:
$$\varphi_U: \pi^{-1}(U) \xrightarrow{\sim} U \times F$$

**Übergangsfunktionen**: Für zwei Karten $U_\alpha, U_\beta$ mit
$U_\alpha \cap U_\beta \neq \emptyset$:
$$g_{\alpha\beta}: U_\alpha \cap U_\beta \to G, \quad g_{\alpha\beta}(b) = \varphi_\alpha \circ \varphi_\beta^{-1}|_b$$

### Verbindungen (Konnexionen)

Eine Konnexion auf einem Prinzipal-$G$-Bündel $P \to B$ ist eine
$\mathfrak{g}$-wertige 1-Form:
$$\omega \in \Omega^1(P, \mathfrak{g})$$

Die **Krümmungsform** $\Omega \in \Omega^2(P, \mathfrak{g})$ ist:
$$\Omega = d\omega + \frac{1}{2}[\omega, \omega]$$

Für abelsche Gruppen (z.B. $U(1)$) verschwindet der Lie-Klammern-Term:
$$\Omega = d\omega = F \quad \text{(Feldstärke)}$$

### Yang-Mills-Theorie

Die Yang-Mills-Wirkung ist:
$$S_{\text{YM}} = -\frac{1}{2} \int_M \text{Tr}(F \wedge *F)$$

Auf dem Gitter (Wilson-Formulierung) mit Link-Variablen $U_\mu(x) \in G$:
$$S_W = \beta \sum_{x,\mu<\nu} \left[1 - \text{Re}\,\text{Tr}(P_{\mu\nu}(x))\right]$$

wobei das Plaquette-Produkt:
$$P_{\mu\nu}(x) = U_\mu(x)\, U_\nu(x+\hat{\mu})\, U_\mu^\dagger(x+\hat{\nu})\, U_\nu^\dagger(x)$$

## API-Referenz

### Klasse `FiberBundle`

```python
FiberBundle(base_dim: int, fiber_dim: int, structure_group: str)
```

**Parameter:**
| Parameter | Typ | Beschreibung |
|-----------|-----|--------------|
| `base_dim` | int | Dimension des Basisraums $B$ |
| `fiber_dim` | int | Dimension der typischen Faser $F$ |
| `structure_group` | str | Strukturgruppe: `'SO2'`, `'U1'`, `'SU2'`, `'GL'`, ... |

**Attribute:**
- `total_dim = base_dim + fiber_dim` – Dimension des Gesamtraums

**Methoden:**

#### `local_trivialization(point) → (basis_coords, fiber_coords)`

Zerlegt einen Punkt $p \in E$ in seine Basis- und Faserkoordinaten.

```python
b = FiberBundle(2, 1, 'U1')
point = np.array([1.0, 2.0, 3.0])
base, fiber = b.local_trivialization(point)
# base = [1.0, 2.0], fiber = [3.0]
```

#### `transition_function(point, chart1, chart2) → np.ndarray`

Berechnet die Übergangsfunktion $g_{\alpha\beta}$ zwischen zwei Karten.

```python
g = b.transition_function(np.array([1.0, 0.0]), 0, 1)
# → 2×2-Rotationsmatrix (für SO(2)/U(1))
```

---

### Klasse `Connection`

```python
Connection(bundle: FiberBundle)
```

**Methoden:**

#### `connection_form(point, vector) → np.ndarray`

Berechnet $\omega(X)$ am Punkt $b$ für den Tangentialvektor $X$.

Modell: $\omega(X) = \frac{b \cdot X}{|b|^2 + 1}$

#### `curvature_form(point) → np.ndarray`

Berechnet die Krümmungsmatrix $\Omega_{\mu\nu}$ durch finite Differenzen:
$$\Omega_{\mu\nu} \approx \partial_\mu A_\nu - \partial_\nu A_\mu$$

Rückgabe: antisymmetrische $n \times n$-Matrix.

#### `parallel_transport(curve, initial_vector, n_steps=100) → np.ndarray`

Transportiert `initial_vector` entlang der Kurve $\gamma: [0,1] \to B$
durch numerische Integration von $\nabla_{\gamma'} V = 0$.

```python
c = Connection(FiberBundle(2, 2, 'U1'))
circle = lambda t: np.array([np.cos(2*np.pi*t), np.sin(2*np.pi*t)])
transported = c.parallel_transport(circle, np.array([1.0, 0.0]))
```

#### `holonomy(loop, base_point) → np.ndarray`

Berechnet die Holonomie: Paralleltransport des Einheitsvektors um eine
geschlossene Kurve.

---

### Klasse `YangMillsField`

```python
YangMillsField(gauge_group: str = 'U1', grid_size: int = 4)
```

Implementiert das Yang-Mills-Gitterfeld in der Wilson-Formulierung.

**Unterstützte Gruppen:** `'U1'`, `'SU2'`

**Methoden:**

#### `field_strength(mu, nu, site) → complex`

Berechnet den Feldstärketensor $F_{\mu\nu}(x)$ via Plaquette-Produkt.
$$F_{\mu\nu} \approx \text{Im}[P_{\mu\nu}]$$

#### `yang_mills_action() → float`

Berechnet die Wilson-Gitterwirkung:
$$S = \beta \sum_{x, \mu<\nu} [1 - \text{Re}(P_{\mu\nu}(x))]$$

Rückgabe ist immer $\geq 0$.

#### `bianchi_identity_check(tolerance=1e-10) → bool`

Stichprobenartige Überprüfung der Bianchi-Identität $D \wedge F = 0$.

---

## Hilfsfunktionen

### `mobius_bundle() → FiberBundle`

Das **Möbius-Bündel**: einfachstes nicht-triviales reelles Linienbündel.
- Basisraum: $S^1$ (dim 1)
- Faser: $\mathbb{R}$ (dim 1)
- Strukturgruppe: $\mathbb{Z}_2 = \{±1\}$
- Nicht-trivial: $w_1 \neq 0$

### `hopf_fibration() → FiberBundle`

Die **Hopf-Faserung** $\eta: S^3 \to S^2$:
- Basisraum: $S^2$ (dim 2)
- Faser: $S^1$ (dim 1)
- Strukturgruppe: $U(1)$
- Erste Chern-Zahl: $c_1 = 1$

Explizit: $(z_1, z_2) \mapsto [z_1 : z_2] \in \mathbb{CP}^1 \cong S^2$

### `tangent_bundle(manifold_dim: int) → FiberBundle`

Das **Tangentialbündel** $TM$ einer $n$-dimensionalen Mannigfaltigkeit:
- Basisraum: $M$ (dim $n$)
- Faser: $T_pM \cong \mathbb{R}^n$ (dim $n$)
- Strukturgruppe: $GL(n, \mathbb{R})$
- Gesamtraum: $TM$ (dim $2n$)

### `chern_class_first(bundle: FiberBundle) → int`

Berechnet die **erste Chern-Klasse** $c_1 \in H^2(B; \mathbb{Z})$ als
topologische Invariante eines $U(1)$-Bündels.

Für die Hopf-Faserung: $c_1 = 1$.

## Bekannte Faserbündel (Übersicht)

| Bündel | Basisraum | Faser | Strukturgruppe | $c_1$ |
|--------|-----------|-------|----------------|-------|
| Möbius-Bündel | $S^1$ | $\mathbb{R}$ | $\mathbb{Z}_2$ | – |
| Hopf-Faserung | $S^2$ | $S^1$ | $U(1)$ | 1 |
| Tangentialbündel $TM^n$ | $M^n$ | $\mathbb{R}^n$ | $GL(n)$ | – |
| Kanonisches Linienbündel | $\mathbb{CP}^n$ | $\mathbb{C}$ | $U(1)$ | 1 |
| Instanton-Bündel | $S^4$ | $SU(2)$ | $SU(2)$ | – |

## Tests

Tests befinden sich in `tests/test_fiber_bundles.py` (28+ Tests):

```bash
python3 -m pytest tests/test_fiber_bundles.py -v
```

## Literaturhinweise

- Nakahara, M.: *Geometry, Topology and Physics* (2003)
- Bleecker, D.: *Gauge Theory and Variational Principles* (1981)
- Atiyah, M.F., Bott, R.: *The Yang-Mills equations over Riemann surfaces* (1983)
- Wilson, K.G.: *Confinement of quarks*, Phys. Rev. D (1974)
