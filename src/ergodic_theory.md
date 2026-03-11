# ergodic_theory.py — Dokumentation

**Autor:** Michael Fuhrmann
**Letzte Änderung:** 2026-03-11
**Build:** 85

---

## Übersicht

Das Modul `ergodic_theory.py` implementiert die zentralen Konzepte der modernen
Ergodentheorie als eigenständiges, mathematisch vollständiges Python-Modul. Es deckt
folgende Themenbereiche ab:

1. **Grundlagen**: Maß-erhaltende dynamische Systeme
2. **Ergodensätze**: Birkhoff (1931), von Neumann (1932)
3. **Entropie**: KS-Entropie, topologische Entropie, Variationsprinzip
4. **Mixing-Eigenschaften**: Starkes/schwaches Mixing, Mischungsraten
5. **Collatz-Dynamik** (Tao 2019): Ergodentheorie des 3n+1-Problems
6. **Furstenberg-Theorie**: Zahlentheorie via Ergodentheorie
7. **Symbolische Dynamik**: Shift-Systeme, SFT, Zeta-Funktionen

---

## Mathematische Grundlagen

### 1. Maß-erhaltendes System

Ein **maß-erhaltendes dynamisches System** ist ein Tripel $(X, T, \mu)$ bestehend aus:
- **Zustandsraum** $X$ (Maßraum)
- **Transformation** $T: X \to X$ (messbar)
- **invariantes Maß** $\mu$ mit $\mu(T^{-1}A) = \mu(A)$ für alle messbaren $A \subseteq X$

**Beispiele:**

| System | $T(x)$ | Invariantes Maß |
|--------|---------|-----------------|
| Kreisrotation | $x + \alpha \pmod{1}$ | Lebesgue |
| Verdopplungsabbildung | $2x \pmod{1}$ | Lebesgue |
| Logistische Abbildung | $4x(1-x)$ | Arkussinus-Maß $\frac{dx}{\pi\sqrt{x(1-x)}}$ |
| Gauß-Abbildung | $\{1/x\}$ | Gauß-Maß $\frac{dx}{(1+x)\log 2}$ |
| Zeltabbildung | $1-|2x-1|$ | Lebesgue |

---

### 2. Ergodizität

Ein System $(X, T, \mu)$ heißt **ergodisch**, wenn alle $T$-invarianten messbaren
Mengen Maß 0 oder 1 haben:

$$\forall A: T^{-1}A = A \Rightarrow \mu(A) \in \{0, 1\}$$

**Equivalent:** Für alle $f, g \in L^2(\mu)$ gilt:
$$\frac{1}{N}\sum_{n=0}^{N-1} \int f \cdot g \circ T^n \, d\mu \to \int f \, d\mu \cdot \int g \, d\mu$$

---

### 3. Birkhoff'sches Ergodentheorem (1931)

Sei $(X, T, \mu)$ maß-erhaltend und $f \in L^1(\mu)$. Dann gilt für $\mu$-fast alle $x$:

$$\lim_{n \to \infty} \frac{1}{n} \sum_{k=0}^{n-1} f(T^k x) = \int_X f \, d\mu$$

Das **Zeitmittel** konvergiert fast überall gegen das **Raummittel**.

Für ergodische Systeme gilt die Konvergenz für **alle** Startpunkte aus einer
Menge vom vollen Maß.

---

### 4. Von-Neumann-Ergodentheorem (1932)

Für $f \in L^2(\mu)$ konvergiert die Folge der Zeitdurchschnitte im $L^2$-Sinne:

$$\left\| \frac{1}{n}\sum_{k=0}^{n-1} f \circ T^k - P_\mathcal{I} f \right\|_{L^2} \to 0$$

wobei $P_\mathcal{I}$ die orthogonale Projektion auf den Raum der $T$-invarianten
Funktionen ist. Für ergodische Systeme gilt $P_\mathcal{I} f = \int f \, d\mu$.

---

### 5. Partitionsentropie

Für eine endliche messbare Partition $\mathcal{P} = \{A_1, \ldots, A_k\}$ von $X$:

$$H_\mu(\mathcal{P}) = -\sum_{i=1}^k \mu(A_i) \log_2 \mu(A_i)$$

Maximiert bei Gleichverteilung: $H_{\max} = \log_2 k$.

---

### 6. Kolmogorov-Sinai-Entropie

Die **KS-Entropie** (Kolmogorov 1958, Sinai 1959) ist definiert als:

$$h_\mu(T) = \sup_{\mathcal{P}} h_\mu(T, \mathcal{P})$$

wobei die **Entropierate** einer Partition ist:

$$h_\mu(T, \mathcal{P}) = \lim_{n \to \infty} \frac{1}{n} H_\mu\!\left(\bigvee_{k=0}^{n-1} T^{-k}\mathcal{P}\right)$$

**Pesin-Formel** (für $C^1$-Diffeomorphismen):

$$h_\mu(T) = \sum_{\lambda_i > 0} \lambda_i$$

(Summe der positiven Lyapunov-Exponenten).

---

### 7. Topologische Entropie

Die **topologische Entropie** misst die Komplexität der Trajektorien:

$$h_{\text{top}}(T) = \lim_{\varepsilon \to 0} \limsup_{n \to \infty} \frac{1}{n} \log s_n(\varepsilon)$$

wobei $s_n(\varepsilon)$ die maximale Anzahl $(n, \varepsilon)$-separierter Punkte ist.

---

### 8. Variationsprinzip (Goodman/Dinaburg 1970/71)

$$h_{\text{top}}(T) = \sup_{\mu \text{ invariant}} h_\mu(T)$$

Das Supremum wird vom **Maß maximaler Entropie** (Parry-Maß) angenommen.

---

### 9. Mixing-Hierarchie

$$\text{Stark mischend} \Rightarrow \text{Schwach mischend} \Rightarrow \text{Ergodisch}$$

**Stark mischend:**
$$\lim_{n \to \infty} \mu(A \cap T^{-n}B) = \mu(A) \cdot \mu(B) \quad \forall A, B$$

**Schwach mischend** (Cesàro):
$$\frac{1}{N} \sum_{n=0}^{N-1} |\mu(A \cap T^{-n}B) - \mu(A)\mu(B)| \to 0$$

---

### 10. Collatz-Abbildung und Tao (2019)

Die **Collatz-Abbildung**:
$$T(n) = \begin{cases} n/2 & n \text{ gerade} \\ 3n+1 & n \text{ ungerade} \end{cases}$$

**Tao (2019):** Für fast alle $n \in \{1, \ldots, N\}$ gilt:
es existiert $k$ mit $T^k(n) \leq f(n)$ für jedes noch so langsam wachsende $f(n) \to \infty$.

**Lyapunov-Exponent des Zufallsmodells:**
$$\lambda = \frac{\log 3 - 2\log 2}{2} \approx -0.0235 < 0$$

Das System ist im Mittel kontraktiv — jeder Schritt reduziert den Wert im Erwartungswert.

**2-adisches Modell** (Lagarias 1985):
Die Collatz-Abbildung auf $\mathbb{Z}_2$ (2-adische Zahlen) ist maß-erhaltend
bezüglich des Haar-Maßes.

---

### 11. Furstenberg-Korrespondenzprinzip (1977)

Jeder Menge $A \subseteq \mathbb{N}$ mit positiver oberer Dichte
$$\bar{d}(A) = \limsup_{N \to \infty} \frac{|A \cap \{1,\ldots,N\}|}{N} > 0$$

entspricht ein maß-erhaltendes System $(X, T, \mu)$ mit $B \subseteq X$, $\mu(B) = \bar{d}(A)$, so dass:
$$n \in A - A \iff \mu(B \cap T^{-n}B) > 0$$

Dies ermöglicht den Beweis des **Szemerédi-Theorems** via Ergodentheorie.

---

### 12. Symbolische Dynamik: SFT

Ein **Subshift of Finite Type** $X_M$ via Übergangsmatrix $M \in \{0,1\}^{d \times d}$:

$$X_M = \{(x_n)_{n \in \mathbb{Z}} \in \{0,\ldots,d-1\}^\mathbb{Z} : M[x_n, x_{n+1}] = 1 \; \forall n\}$$

**Topologische Entropie:**
$$h_{\text{top}}(\sigma|_{X_M}) = \log \lambda_{\max}(M)$$

wobei $\lambda_{\max}$ der Perron-Frobenius-Eigenwert ist.

**Dynamische Zeta-Funktion** (Artin-Mazur):
$$\zeta_T(z) = \exp\!\left(\sum_{n=1}^\infty \frac{|\text{Fix}(T^n)|}{n} z^n\right) = \frac{1}{\det(I - zM)}$$

**Golden-Mean-Shift** ($M = \begin{pmatrix}1&1\\1&0\end{pmatrix}$):
$$h = \log\!\left(\frac{1+\sqrt{5}}{2}\right) = \log \varphi \approx 0.481$$

---

## Klassen und Funktionen

### `MeasurePreservingSystem`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `orbit` | `(x0, n) → np.ndarray` | Berechnet Orbit $x_0, T(x_0), \ldots, T^n(x_0)$ |
| `is_measure_preserving` | `(test_points, n_bins) → bool` | $\chi^2$-Test auf Maß-Erhaltung |
| `is_ergodic` | `(x0, n_iter, n_bins) → bool` | Histogramm-Test auf Ergodizität |
| `ergodic_decomposition` | `(initial_points, n_iter, n_clusters) → List` | Zerlegung in Komponenten |

### `BirkhoffErgodicTheorem`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `birkhoff_ergodic_theorem` | `(T, f, x0, n) → np.ndarray` | Zeitdurchschnitte $A_1, \ldots, A_n$ |
| `von_neumann_mean_ergodic` | `(f_values, space_mean) → Dict` | $L^2$-Konvergenz |
| `ergodic_average_convergence` | `(T, f, x0, space_mean, n_steps) → Dict` | Konvergenzverifikation |

### `ErgodicEntropy`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `partition_entropy` | `(probs) → float` | $H(\mathcal{P}) = -\sum p_i \log_2 p_i$ |
| `metric_entropy` | `(T, x0, bins, n_iter) → float` | KS-Entropie-Schätzung |
| `topological_entropy` | `(T, n_orbits, n_iter, eps) → float` | $h_{\text{top}}$ via separierte Mengen |
| `variational_principle_check` | `(T, x0_list, bins, n_iter) → Dict` | $h_{\text{top}} = \sup h_\mu$ |

### `MixingProperties`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `correlation_function` | `(trajectory, max_lag) → np.ndarray` | Normierte ACF $\rho(0), \ldots, \rho(n)$ |
| `is_mixing` | `(T, f, x0, n_iter, lag) → bool` | Starkes Mixing-Kriterium |
| `is_weak_mixing` | `(T, f, x0, n_iter, max_lag) → bool` | Cesàro-Konvergenz-Kriterium |
| `mixing_rate_estimate` | `(trajectory, max_lag) → Dict` | Exp. Fit: $C(n) \approx C_0 e^{-\gamma n}$ |

### `CollatzErgodicSystem`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `collatz_map` | `(n) → int` | $T(n) = n/2$ oder $3n+1$ |
| `collatz_iteration` | `(n, max_steps) → List[int]` | Vollständige Trajektorie bis 1 |
| `collatz_stopping_time` | `(n, max_steps) → int` | $\sigma(n) = \min\{k: T^k(n) < n\}$ |
| `collatz_total_stopping_time` | `(n, max_steps) → int` | Schritte bis $T^k(n) = 1$ |
| `collatz_density_tao` | `(N, eps, C) → Dict` | Taos Dichte-1-Resultat |
| `collatz_lyapunov_exponent` | `(n_samples) → float` | Empirischer Lyapunov-Exp. |
| `collatz_theoretical_lyapunov` | `() → float` | $\lambda = (\log 3 - 2\log 2)/2$ |
| `collatz_2adic_analysis` | `(n, depth) → Dict` | 2-adische Darstellung |
| `collatz_ergodic_model` | `(N) → Dict` | Gleichvert. mod $2^k$ |
| `collatz_stopping_time_distribution` | `(N) → Dict` | Normalverteilungs-Fit |

### `FurstenbergTheory`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `upper_density` | `(subset, N) → float` | $\bar{d}(A) = |A \cap \{1,\ldots,N\}|/N$ |
| `furstenberg_correspondence` | `(subset, N) → Dict` | Korrespondenzprinzip |
| `szemeredi_via_ergodic` | `(k, N) → Dict` | Szemerédi AP-Länge $k$ |
| `van_der_waerden_bound` | `(k, r) → Dict` | $W(k;r)$ Schranken |

### `ShiftSystem`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `shift` | `(sequence, steps) → List` | $\sigma^n(x)$ |
| `is_in_shift_space` | `(sequence, forbidden_words) → bool` | Zulässigkeitsprüfung |
| `count_words` | `(length, M) → int` | $|W_n(X_M)|$ |
| `entropy_shift_sft` | `(M) → float` | $h = \log \lambda_{\max}(M)$ |
| `zeta_function_shift` | `(M, max_n) → Dict` | $\zeta_T(z)$ |

### `SubshiftOfFiniteType`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `is_irreducible` | `() → bool` | Prüft Irreduzibilität von $M$ |
| `is_mixing` | `() → bool` | Topologisches Mixing |
| `entropy` | `() → float` | $h = \log \lambda_{\max}$ |
| `parry_measure` | `() → np.ndarray` | Parry-Maß (max. Entropie-Maß) |
| `generate_sequence` | `(length, start, seed) → List` | Zufälliger Pfad im SFT |

### Hilfsfunktionen

| Funktion | Formel | Eigenschaften |
|----------|--------|---------------|
| `doubling_map(x)` | $2x \bmod 1$ | Ergodisch, stark mischend, $h = \log 2$ |
| `circle_rotation(α)(x)` | $x + \alpha \bmod 1$ | Ergodisch gdw. $\alpha$ irrational, $h=0$ |
| `tent_map(x)` | $1-|2x-1|$ | Ergodisch, stark mischend, $h = \log 2$ |
| `logistic_map_r4(x)` | $4x(1-x)$ | Ergodisch, Arkussinus-Maß, $h=\log 2$ |
| `gauss_map(x)` | $\{1/x\}$ | Gauß-Maß, $h = \pi^2/(6\log 2) \approx 2.37$ |

---

## Beispiele

### Birkhoff-Ergodentheorem numerisch verifizieren

```python
from ergodic_theory import BirkhoffErgodicTheorem, logistic_map_r4

bw = BirkhoffErgodicTheorem()

# Zeitmittel → Raummittel (E[x] = 0.5 unter Arkussinus-Maß)
avgs = bw.birkhoff_ergodic_theorem(logistic_map_r4, lambda x: x, 0.3, 5000)
print(f"Zeitmittel nach 5000 Schritten: {avgs[-1]:.4f}")  # → ~0.5
```

### Collatz-Lyapunov-Exponent

```python
from ergodic_theory import CollatzErgodicSystem

cs = CollatzErgodicSystem()
lam_theor = cs.collatz_theoretical_lyapunov()
print(f"Theoretischer Lyapunov: {lam_theor:.6f}")  # ≈ -0.0235
```

### Golden-Mean-Shift Entropie

```python
import numpy as np
from ergodic_theory import SubshiftOfFiniteType

M = np.array([[1, 1], [1, 0]])  # Verbotenes Wort: "11"
sft = SubshiftOfFiniteType(M)
print(f"Entropie: {sft.entropy():.6f}")  # ≈ log(φ) ≈ 0.4812
print(f"Mischend: {sft.is_mixing()}")   # True
```

### Szemerédi-Theorem numerisch

```python
from ergodic_theory import FurstenbergTheory

ft = FurstenbergTheory()
results = ft.szemeredi_via_ergodic(k=3, N=500)
# Bei Dichte 0.3: AP der Länge 3 gefunden?
print(results)
```

---

## Literaturhinweise

| Referenz | Inhalt |
|----------|--------|
| Birkhoff (1931) | Ergodentheorem (fast überall Konvergenz) |
| von Neumann (1932) | Mittelwertergodentheorem ($L^2$-Konvergenz) |
| Kolmogorov (1958), Sinai (1959) | KS-Entropie |
| Parry (1964) | Maß maximaler Entropie für SFT |
| Furstenberg (1977) | Korrespondenzprinzip, Szemerédi via Ergodentheorie |
| Lagarias (1985) | Collatz auf 2-adischen Zahlen |
| Tao (2019) | Almost all orbits of Collatz attain almost bounded values |
| Artin-Mazur (1965) | Dynamische Zeta-Funktion |
| Goodman (1971), Dinaburg (1970) | Variationsprinzip |
