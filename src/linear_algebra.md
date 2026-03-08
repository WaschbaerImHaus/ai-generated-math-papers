# linear_algebra.py – Lineare-Algebra-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/linear_algebra.py`

---

## Überblick

Dieses Modul implementiert fundamentale Konzepte der Linearen Algebra:

| Klasse/Funktion | Beschreibung |
|----------------|-------------|
| `Vector` | n-dimensionaler Vektor mit allen Operationen |
| `Matrix` | m×n-Matrix mit Gauss-Elimination, Eigenwerten |
| `gram_schmidt` | Gram-Schmidt-Orthogonalisierung |
| `lu_decomposition` | LU-Zerlegung (Doolittle mit Teilpivotisierung) |
| `qr_decomposition` | QR-Zerlegung (Householder-Reflexionen) |
| `svd` | Singulärwertzerlegung (via NumPy LAPACK) |
| `matrix_rank` | Rang-Berechnung via SVD |
| `condition_number` | Konditionszahl κ(A) = σ_max/σ_min |

---

## Klasse: `Vector`

Repräsentiert einen n-dimensionalen reellen Vektor.

### Konstruktor

```python
v = Vector([1, 2, 3])   # 3D-Vektor
```

### Methoden

#### `dot(other) → float`
**Skalarprodukt (inneres Produkt):**
```
v · w = v₁·w₁ + v₂·w₂ + ... + vₙ·wₙ
```

Geometrisch: `v · w = |v|·|w|·cos(θ)`, wobei θ der Winkel zwischen den Vektoren ist.

**Wichtige Eigenschaften:**
- `v · w = 0` genau dann wenn v ⊥ w (senkrecht)
- `v · v = |v|²`

#### `norm() → float`
**Euklidische Norm (L2-Norm):**
```
||v|| = √(v₁² + v₂² + ... + vₙ²) = √(v · v)
```

#### `normalize() → Vector`
Gibt den **Einheitsvektor** zurück: `v̂ = v / ||v||`
Eigenschaften: `||v̂|| = 1`, gleiche Richtung wie v.

#### `cross(other) → Vector`
**Kreuzprodukt** (nur für 3D-Vektoren):
```
v × w = | i  j  k  |
        | v₁ v₂ v₃ |
        | w₁ w₂ w₃ |
```

**Eigenschaften:**
- `v × w ⊥ v` und `v × w ⊥ w`
- `|v × w| = |v|·|w|·sin(θ)` = Flächeninhalt des Parallelogramms
- Antikommutativ: `v × w = -(w × v)`

---

## Klasse: `Matrix`

Repräsentiert eine m×n-Matrix. Intern gespeichert als `List[List[float]]`.

### Konstruktor

```python
A = Matrix([[1, 2], [3, 4]])   # 2×2-Matrix
```

### Klassenmethoden

| Methode | Beschreibung |
|---------|-------------|
| `Matrix.identity(n)` | n×n-Einheitsmatrix (1 auf Diagonale) |
| `Matrix.zeros(rows, cols)` | Nullmatrix |

### Operationen

#### `__matmul__(other)` (A @ B)
**Matrixmultiplikation:**
```
(A @ B)[i,j] = Σₖ A[i,k] · B[k,j]
```
= Skalarprodukt der i-ten Zeile von A mit der j-ten Spalte von B.
**Voraussetzung:** A ist m×n, B muss n×p sein → Ergebnis ist m×p.

#### `transpose() → Matrix`
**Transposition:** `(A^T)[i,j] = A[j,i]`

**Eigenschaften:**
- `(A^T)^T = A`
- `(A @ B)^T = B^T @ A^T`
- Für symmetrische Matrizen: `A^T = A`

#### `trace() → float`
**Spur** = Summe der Diagonalelemente.
Eigenschaft: Spur = Summe aller Eigenwerte.

#### `determinant() → float`
**Determinante** via Gauss-Elimination mit Pivotsuche.

Die Determinante misst den "Volumenskalierungsfaktor" der linearen Abbildung:
- `det(A) = 0` ⟺ A ist singulär (nicht invertierbar)
- `det(A @ B) = det(A) · det(B)`
- Zeilenvertauschung → Vorzeichenwechsel

**Algorithmus:** Elimination zur oberen Dreiecksform, dann `det = ∏ Diagonalelemente · (-1)^Tausche`.

#### `inverse() → Matrix`
**Inverse** via Gauss-Jordan-Elimination.

Erweiterte Matrix `[A | I]` wird durch Zeilenoperationen in `[I | A⁻¹]` überführt.
Existiert genau dann wenn `det(A) ≠ 0`.

#### `solve(b) → Vector`
Löst das lineare Gleichungssystem `Ax = b` via Gauss-Elimination + Rücksubstitution.
Numerisch stabiler als `x = A⁻¹ @ b`.

#### `eigenvalues() → list`
Berechnet die Eigenwerte. Ein Eigenvektor `v` und Eigenwert `λ` erfüllen `A @ v = λ·v`.

**Für 2×2-Matrizen:** Analytische Formel über charakteristisches Polynom:
```
λ² - Spur(A)·λ + det(A) = 0
```

**Für größere Matrizen:** **QR-Iteration** – iteriert `A_{k+1} = R_k @ Q_k` bis zur Schur-Form.

---

## Funktion: `gram_schmidt`

```python
gram_schmidt(vectors, normalize=False) → List[Vector]
```

**Gram-Schmidt-Orthogonalisierung:** Erzeugt aus einer Basis eine orthogonale (oder orthonormale) Basis.

### Algorithmus

```
u₁ = v₁
uₖ = vₖ - Σⱼ₌₁^{k-1} proj_{uⱼ}(vₖ)
```

Projektion: `proj_u(v) = (v · u / u · u) · u`

**Geometrische Idee:** Jeder neue Vektor wird so modifiziert, dass er senkrecht auf allen vorherigen steht.

Mit `normalize=True` → **Orthonormalbasis (ONB)**: `||uᵢ|| = 1` und `uᵢ · uⱼ = 0` für `i ≠ j`.

---

## Funktion: `lu_decomposition`

```python
lu_decomposition(matrix) → (L, U, P, n_swaps)
```

**LU-Zerlegung (Doolittle-Algorithmus mit Teilpivotisierung):**

```
P·A = L·U
```

| Matrix | Beschreibung |
|--------|-------------|
| P | Permutationsmatrix (Zeilenvertauschungen) |
| L | Untere Dreiecksmatrix (Diagonale = 1) |
| U | Obere Dreiecksmatrix |

**Teilpivotisierung:** Pro Spalte wird das betragsmäßig größte Element als Pivot gewählt → numerisch stabiler.

### Anwendungen
- **LGS lösen:** `Ax = b` → `LUx = Pb` → `Lz = Pb` (Vorwärts) → `Ux = z` (Rückwärts)
- **Determinante:** `det(A) = det(U) · (-1)^n_swaps = ∏ U[i,i] · (-1)^n_swaps`

### Doolittle-Formel
```
L[i,k] = A[i,k] / A[k,k]    (Eliminationsfaktor)
A[i,j] -= L[i,k] · A[k,j]   (in-place Update)
```

---

## Funktion: `qr_decomposition`

```python
qr_decomposition(matrix) → (Q, R)
```

**QR-Zerlegung via Householder-Reflexionen:**

```
A = Q·R
```

| Matrix | Beschreibung |
|--------|-------------|
| Q | Orthogonale Matrix (Q^T·Q = I), m×m |
| R | Obere Dreiecksmatrix, m×n |

### Householder-Reflexion
```
H(v) = I - 2·v·vᵀ / (vᵀ·v)
```

Jede Reflexion macht eine **Spalte unterhalb der Diagonale null**.
**Numerisch stabiler als Gram-Schmidt-QR** (kein Auslöschungsproblem).

### Householder-Vektor
```
v = x + sign(x₀)·||x||·e₁
```
Das Vorzeichen wird so gewählt, dass Auslöschung vermieden wird.

---

## Funktion: `svd`

```python
svd(matrix) → (U, sigma, Vt)
```

**Singulärwertzerlegung (Singular Value Decomposition):**

```
A = U·Σ·Vᵀ
```

| Symbol | Beschreibung |
|--------|-------------|
| U | Unitäre Matrix (m×m), linke Singulärvektoren |
| Σ | Diagonalmatrix mit Singulärwerten σ₁ ≥ σ₂ ≥ ... ≥ 0 |
| Vᵀ | Unitäre Matrix (n×n), rechte Singulärvektoren transponiert |

**Wichtige Eigenschaften:**
- Singulärwerte σᵢ sind immer ≥ 0
- Matrixrang = Anzahl σᵢ > 0
- Pseudoinverse: `A⁺ = V·Σ⁺·Uᵀ`
- `||A||₂ = σ₁` (Spektralnorm)

### Anwendungen
- **PCA (Hauptkomponentenanalyse)**
- **Dimensionsreduktion / Datenkomprimierung**
- **Niedrigrangapproximation:** `Aₖ = Σᵢ₌₁ᵏ σᵢ·uᵢ·vᵢᵀ`

---

## Funktion: `matrix_rank`

```python
matrix_rank(matrix, tol=1e-10) → int
```

Rang = Anzahl der Singulärwerte > tol.
Numerisch robuster als Gauss-Elimination bei fast-singulären Matrizen.

---

## Funktion: `condition_number`

```python
condition_number(matrix) → float
```

**Konditionszahl:**
```
κ(A) = σ_max / σ_min
```

Misst die Empfindlichkeit des LGS gegenüber Störungen:
```
||Δx|| / ||x|| ≤ κ(A) · ||Δb|| / ||b||
```

| κ | Bedeutung |
|---|-----------|
| κ ≈ 1 | Gut konditioniert |
| κ >> 1 | Schlecht konditioniert (numerische Probleme) |

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `copy` | Tiefe Kopien für Sicherheit |
| `numpy` | QR-Iteration für Eigenwerte, LAPACK für SVD |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from linear_algebra import Vector, Matrix, gram_schmidt, lu_decomposition, svd

# Vektoren
v = Vector([1, 2, 3])
w = Vector([4, 5, 6])
print(v.dot(w))       # 32
print(v.cross(w))     # [-3, 6, -3]
print(v.norm())       # 3.7417...

# Matrix-Operationen
A = Matrix([[2, 1], [5, 3]])
b = Vector([4, 7])
x = A.solve(b)        # Lösung von Ax = b
print(x)              # [5, -6]

print(A.determinant())    # 1.0
print(A.eigenvalues())    # [4.303, 0.697] (näherungsweise)

# LU-Zerlegung
L, U, P, swaps = lu_decomposition(A)

# SVD
U_mat, sigma, Vt = svd(A)
print(sigma)   # Singulärwerte

# Gram-Schmidt
v1 = Vector([1, 1, 0])
v2 = Vector([1, 0, 1])
ortho = gram_schmidt([v1, v2], normalize=True)
```
