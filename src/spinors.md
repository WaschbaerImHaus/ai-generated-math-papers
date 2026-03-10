# spinors.py — Dokumentation

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Version:** 1.0.0

---

## Überblick

Das Modul `spinors.py` implementiert die grundlegenden mathematischen Strukturen
der Spinorgeometrie und der relativistischen Quantenmechanik. Es umfasst:

- **Clifford-Algebren** $\mathrm{Cl}(1,3)$ mit Gamma-Matrizen in Dirac-Darstellung
- **Spinoren** als Elemente des Darstellungsraums der Spin-Gruppe $\mathrm{Spin}(n)$
- **Dirac-Gleichung** im Impulsraum
- **Weyl-Spinoren** (chirale Projektion)
- **Majorana-Fermionen** (Teilchen = eigenes Antiteilchen)
- **Gitter-Dirac-Operator** (Wilson-Fermionen, 1D)

---

## Mathematischer Hintergrund

### Clifford-Algebra und Gamma-Matrizen

Die Clifford-Algebra $\mathrm{Cl}(1,3)$ über der Minkowski-Raumzeit wird durch
vier $4 \times 4$-Matrizen $\gamma^\mu$ ($\mu = 0,1,2,3$) realisiert, die die
**fundamentale Anti-Kommutatorrelation** erfüllen:

$$\{\gamma^\mu, \gamma^\nu\} = \gamma^\mu \gamma^\nu + \gamma^\nu \gamma^\mu = 2\eta^{\mu\nu} \cdot I$$

mit der Minkowski-Metrik $\eta = \mathrm{diag}(+1, -1, -1, -1)$.

In der **Dirac-Darstellung** lauten die Gamma-Matrizen explizit:

$$\gamma^0 = \begin{pmatrix} I_2 & 0 \\ 0 & -I_2 \end{pmatrix}, \quad
\gamma^k = \begin{pmatrix} 0 & \sigma_k \\ -\sigma_k & 0 \end{pmatrix} \; (k=1,2,3)$$

wobei $\sigma_k$ die Pauli-Matrizen sind.

### Chiralitäts-Matrix $\gamma^5$

$$\gamma^5 = i \cdot \gamma^0 \gamma^1 \gamma^2 \gamma^3$$

**Eigenschaften:**
- $(\gamma^5)^2 = I$
- $\{\gamma^5, \gamma^\mu\} = 0$ für alle $\mu$
- $\mathrm{Tr}(\gamma^5) = 0$
- $(\gamma^5)^\dagger = \gamma^5$ (hermitesch)

### Dirac-Gleichung

Im Impulsraum (natürliche Einheiten $\hbar = c = 1$):

$$(\gamma^\mu p_\mu - m \cdot I) u(p) = 0$$

mit dem **Feynman-Slash** $\not{p} = \gamma^\mu p_\mu = \gamma^0 E - \gamma^1 p_x - \gamma^2 p_y - \gamma^3 p_z$
und der Dispersionsrelation $E = \sqrt{m^2 + |\mathbf{p}|^2}$.

### Ebene-Wellen-Lösungen (Dirac-Spinoren)

Die positiv-Energie-Lösungen lauten:

$$u(p, s) = N \begin{pmatrix} \chi_s \\ \frac{\boldsymbol{\sigma} \cdot \mathbf{p}}{E+m} \chi_s \end{pmatrix}, \quad N = \sqrt{E+m}$$

mit $\chi_0 = (1, 0)^T$ (Spin-up) und $\chi_1 = (0, 1)^T$ (Spin-down).

### Spin-Gruppe SU(2) ≅ Spin(3)

Ein Element der Spin-Gruppe entspricht einer Rotation um Achse $\hat{n}$ mit Winkel $\theta$:

$$U(\theta, \hat{n}) = \cos\frac{\theta}{2} \cdot I + i\sin\frac{\theta}{2} \cdot (\boldsymbol{\sigma} \cdot \hat{n})$$

**Spinor-Statistik (fundamentales Merkmal):**

$$U(2\pi, \hat{n}) = -I \quad \text{(Spinor ändert Vorzeichen nach } 2\pi\text{-Rotation!)}$$
$$U(4\pi, \hat{n}) = +I \quad \text{(erst nach } 4\pi \text{ kehrt Spinor zurück)}$$

### Weyl-Spinoren (chirale Projektion)

$$\psi_L = P_L \psi = \frac{I - \gamma^5}{2} \psi \quad \text{(linkshändig)}$$
$$\psi_R = P_R \psi = \frac{I + \gamma^5}{2} \psi \quad \text{(rechtshändig)}$$

Projektoreigenschaften: $P_L^2 = P_L$, $P_R^2 = P_R$, $P_L P_R = 0$, $P_L + P_R = I$.

### Majorana-Bedingung

Ein Majorana-Fermion ist sein eigenes Antiteilchen:

$$\psi = \psi^c = C \bar{\psi}^T, \quad C = i\gamma^2\gamma^0$$

wobei $\bar{\psi} = \psi^\dagger \gamma^0$ der Dirac-adjungierte Spinor ist.

### Wilson-Fermion-Hamiltonian (1D-Gitter)

Der Wilson-Fermion-Hamiltonian verhindert das **Fermion-Doubling**
(Nielsen-Ninomiya-Theorem) durch einen zusätzlichen Wilson-Term $r/a$:

$$H = \sum_n \left[ \left(m + \frac{r}{a}\right) \bar{\psi}_n \psi_n
- \frac{1}{2a}\left(\bar{\psi}_n(\gamma^1 - r\cdot I)\psi_{n+1} + \text{h.c.}\right) \right]$$

---

## API-Referenz

### `gamma_matrices(dim=4) -> list[np.ndarray]`

Gibt die Gamma-Matrizen der Clifford-Algebra zurück.

| Parameter | Typ | Beschreibung |
|-----------|-----|--------------|
| `dim` | `int` | Raumzeit-Dimension: 2 (2×2) oder 4 (4×4 Dirac) |

**Rückgabe:** Liste mit `dim` Gamma-Matrizen.

---

### `clifford_algebra_check(gammas, metric=None) -> dict`

Überprüft die Clifford-Algebra-Relation.

**Rückgabe:**
```python
{
    'max_error': float,          # Maximaler Frobenius-Fehler
    'is_clifford_algebra': bool, # True wenn max_error < 1e-10
    'anticommutators': dict      # {(μ,ν): Fehler}
}
```

---

### `gamma5(dim=4) -> np.ndarray`

Berechnet die Chiralitäts-Matrix $\gamma^5 = i\gamma^0\gamma^1\gamma^2\gamma^3$.

---

### `dirac_spinor(mass, momentum, spin=0) -> np.ndarray`

Berechnet einen Dirac-Spinor $u(p, s)$ als Lösung der Dirac-Gleichung.

| Parameter | Typ | Beschreibung |
|-----------|-----|--------------|
| `mass` | `float` | Ruhemasse $m \geq 0$ |
| `momentum` | `np.ndarray` | 3-Impuls $[p_x, p_y, p_z]$ |
| `spin` | `int` | 0 = Spin-up, 1 = Spin-down |

---

### `dirac_equation_check(psi, mass, momentum) -> dict`

Überprüft ob $(\not{p} - m\cdot I)\psi = 0$.

**Rückgabe:**
```python
{
    'residual': float,       # ||( p̸ - m·I) ψ||₂
    'satisfies_dirac': bool, # True wenn residual < 1e-10
    'energy': float          # E = √(m²+|p|²)
}
```

---

### `spin_group_element(theta, axis) -> np.ndarray`

Berechnet $U(\theta, \hat{n}) \in \mathrm{SU}(2)$.

---

### `weyl_spinors(psi) -> dict`

Zerlegt in chirale Komponenten $\psi_L, \psi_R$.

---

### `majorana_condition_check(psi) -> dict`

Prüft $\psi = C\bar{\psi}^T$ (Majorana-Bedingung).

---

### `dirac_hamiltonian_1d(n_sites, mass, lattice_spacing) -> np.ndarray`

Wilson-Fermion-Hamiltonian auf 1D-Gitter der Größe $(2 n_\text{sites} \times 2 n_\text{sites})$.

---

### `dirac_spectrum_1d(mass, n_sites) -> dict`

Berechnet das Spektrum des 1D-Gitter-Dirac-Operators.

**Rückgabe:**
```python
{
    'eigenvalues': np.ndarray, # Alle Eigenwerte (reell)
    'mass_gap': float,         # min(|Eigenwerte|)
    'is_gapped': bool,         # True wenn mass_gap > 1e-10
    'n_zero_modes': int        # Anzahl der Nullmoden
}
```

---

## Tests

50 Tests in `tests/test_spinors.py`, alle grün (Build 23):

| Testklasse | Anzahl | Was wird geprüft |
|------------|--------|------------------|
| `TestGammaMatrices4D` | 8 | Clifford-Relationen, Hermitizität |
| `TestCliffordAlgebraCheck` | 5 | Vollständige Algebra-Überprüfung |
| `TestGamma5` | 5 | $(γ^5)^2=I$, Anti-Kommutator, Spur |
| `TestDiracSpinor` | 6 | Dirac-Gleichung, Normierung |
| `TestSpinGroupElement` | 5 | Unitarität, Spinor-Statistik |
| `TestWeylSpinors` | 5 | Vollständigkeit, Projektoreigenschaft |
| `TestMajoranaCondition` | 4 | Majorana-Bedingung, Residuum |
| `TestLatticeDirac` | 7 | Hermitizität, Massenlücke, Spektrum |
| `TestGammaMatrices2D` | 3 | 2D-Darstellung |
| `TestGamma5EdgeCases` | 2 | Fehlerbehandlung |

---

## Physikalische Anwendungen

- **Quanten-Elektrodynamik (QED):** Elektronen als Dirac-Spinoren, Photonen als Eichbosonen
- **Schwache Wechselwirkung:** Nur linkshändige Quarks/Leptonen koppeln an W-Bosonen
- **Gitter-QCD:** Wilson-Fermionen auf dem Raumzeit-Gitter
- **Majorana-Neutrinos:** Experimentell noch ungeklärte Frage ob Neutrinos Majorana-Fermionen sind
- **Topologische Supraleiter:** Majorana-Zustände an Randflächen als Quanten-Computing-Ressource
