# Funktionalanalysis – Modul `functional_analysis.py`

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Sprache:** Python 3.13.7

---

## Überblick

Das Modul `functional_analysis.py` implementiert die zentralen Konzepte der
Funktionalanalysis – dem Teilgebiet der Mathematik, das Analysis und lineare
Algebra auf unendlich-dimensionale Räume verallgemeinert.

Die Funktionalanalysis bildet das theoretische Fundament für:
- **Partielle Differentialgleichungen** (schwache Formulierungen, Sobolev-Räume)
- **Quantenmechanik** (Hilbert-Raum-Formalismus, Observable als Operatoren)
- **Numerische Mathematik** (Approximationstheorie, Fehleranalyse)
- **Signalverarbeitung** (Fourier-Analyse in L²)

---

## Hierarchie der Räume

```
NormedSpace (V, ‖·‖)
    └── BanachSpace (vollständig)
            └── HilbertSpace (mit ⟨·,·⟩)
```

---

## 1. Normierte Räume (`NormedSpace`)

### Definition

Ein **normierter Raum** ist ein Vektorraum $V$ über $\mathbb{R}$ mit einer Norm
$\|\cdot\| : V \to \mathbb{R}_{\geq 0}$, die folgende Axiome erfüllt:

$$\|v\| \geq 0, \quad \|v\| = 0 \Leftrightarrow v = 0$$
$$\|\alpha v\| = |\alpha| \cdot \|v\| \quad \text{(Homogenität)}$$
$$\|u + v\| \leq \|u\| + \|v\| \quad \text{(Dreiecksungleichung)}$$

### Induzierte Metrik

Jede Norm induziert eine Metrik:
$$d(u, v) = \|u - v\|$$

### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `norm(v)` | Berechnet $\|v\|$ |
| `distance(u, v)` | Berechnet $\|u - v\|$ |
| `is_normed_space()` | Prüft Normaxiome |
| `unit_ball(n_points)` | Punkte auf dem Rand der Einheitskugel |

---

## 2. Banach-Räume (`BanachSpace`)

### Definition

Ein **Banach-Raum** ist ein vollständiger normierter Raum: jede Cauchy-Folge konvergiert.

**Cauchy-Folge:** $(x_n)$ mit $\|x_m - x_n\| \to 0$ für $m, n \to \infty$

### Die $\ell^p$-Normen

$$\|v\|_p = \left(\sum_i |v_i|^p\right)^{1/p}, \quad 1 \leq p < \infty$$
$$\|v\|_\infty = \max_i |v_i|$$

**Spezialfälle:**
- $p = 1$: Betragssumme (Manhattan-Norm)
- $p = 2$: Euklidische Norm
- $p = \infty$: Tschebyschew-Norm (Maximum)

### Klassische Banach-Räume

| Raum | Norm |
|------|------|
| $\mathbb{R}^n$ | $\|\cdot\|_p$ |
| $\ell^p = \{(x_n) : \sum|x_n|^p < \infty\}$ | $\|x\|_p$ |
| $C[a,b]$ (stetige Funktionen) | $\|f\|_\infty = \max|f(x)|$ |
| $L^p[a,b]$ | $\|f\|_p = (\int|f|^p)^{1/p}$ |

### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `cauchy_sequence_test(seq, tol)` | Prüft ob Folge Cauchy-Folge ist |
| `is_separable()` | Gibt an ob Raum separabel ist |
| `lp_norm(v, p)` | Berechnet $\ell^p$-Norm |

---

## 3. Hilbert-Räume (`HilbertSpace`)

### Definition

Ein **Hilbert-Raum** ist ein Banach-Raum mit einem Skalarprodukt
$\langle \cdot, \cdot \rangle : H \times H \to \mathbb{R}$, das die Norm induziert:

$$\|v\| = \sqrt{\langle v, v \rangle}$$

### Skalarprodukt-Axiome

$$\langle \alpha u + \beta v, w \rangle = \alpha\langle u,w\rangle + \beta\langle v,w\rangle$$
$$\langle u, v \rangle = \langle v, u \rangle \quad \text{(Symmetrie)}$$
$$\langle v, v \rangle \geq 0, \quad \langle v, v \rangle = 0 \Leftrightarrow v = 0$$

### Wichtige Ungleichungen

**Cauchy-Schwarz:**
$$|\langle u, v \rangle| \leq \|u\| \cdot \|v\|$$

**Parallelogramm-Gesetz:**
$$\|u + v\|^2 + \|u - v\|^2 = 2(\|u\|^2 + \|v\|^2)$$

### Orthogonalität und Gram-Schmidt

Vektoren $u, v$ heißen **orthogonal** ($u \perp v$), wenn $\langle u, v \rangle = 0$.

Das **Gram-Schmidt-Verfahren** erzeugt eine Orthonormalbasis (ONB):
$$u_k = v_k - \sum_{j<k} \langle v_k, e_j \rangle e_j, \quad e_k = \frac{u_k}{\|u_k\|}$$

### Projektionen

Die **Orthogonalprojektion** auf $U = \text{span}\{e_1, \ldots, e_n\}$:
$$P_U(v) = \sum_i \langle v, e_i \rangle \cdot e_i$$

**Orthogonalitätsprinzip:** $v - P_U(v) \perp U$

### Parseval-Gleichung

Für eine vollständige ONB $\{e_i\}$:
$$\|v\|^2 = \sum_i |\langle v, e_i \rangle|^2$$

### Bessel-Ungleichung

Für jede orthonormale Menge $\{e_1, \ldots, e_n\}$:
$$\sum_i |\langle v, e_i \rangle|^2 \leq \|v\|^2$$

### Riesz-Darstellungssatz

**Satz (Riesz 1907):** Zu jedem stetigen linearen Funktional $f \in H^*$ existiert
genau ein $y \in H$ mit:
$$f(x) = \langle x, y \rangle \quad \text{für alle } x \in H$$

Außerdem gilt $\|f\|_{H^*} = \|y\|_H$.

### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `inner_product(u, v)` | Berechnet $\langle u, v \rangle$ |
| `is_orthogonal(u, v)` | Prüft $\langle u, v \rangle = 0$ |
| `gram_schmidt_orthonormalize(basis)` | ONB-Erzeugung |
| `projection(v, onto)` | Orthogonalprojektion |
| `pythagoras(u, v)` | Prüft $\|u+v\|^2 = \|u\|^2 + \|v\|^2$ |
| `riesz_representation(f, basis)` | Riesz-Repräsentant |
| `parseval_identity(v, onb)` | Parseval-Prüfung |
| `bessel_inequality(v, ortho_set)` | Bessel-Prüfung |

---

## 4. Lineare Operatoren (`LinearOperator`)

### Definition

$T: X \to Y$ heißt **linear**, wenn:
$$T(\alpha u + \beta v) = \alpha T(u) + \beta T(v)$$

$T$ heißt **beschränkt**, wenn:
$$\|T\| = \sup\left\{\frac{\|Tv\|_Y}{\|v\|_X} : v \neq 0\right\} < \infty$$

### Operatornormen (Matrizen)

$$\|T\|_2 = \text{größter Singulärwert}$$
$$\|T\|_1 = \max_j \sum_i |t_{ij}| \quad \text{(max Spaltensumme)}$$
$$\|T\|_\infty = \max_i \sum_j |t_{ij}| \quad \text{(max Zeilensumme)}$$

### Adjungierter Operator

Der **adjungierte Operator** $T^*: Y \to X$ ist definiert durch:
$$\langle Tv, w \rangle_Y = \langle v, T^* w \rangle_X$$

Für Matrizen: $T^* = T^\top$ (reell).

### Klassen linearer Operatoren

| Eigenschaft | Definition |
|-------------|------------|
| Selbstadjungiert | $T = T^*$ |
| Unitär | $T^* T = T T^* = I$ |
| Normal | $T^* T = T T^*$ |
| Projektion | $T^2 = T$ |
| Kompakt | Einheitskugel → präkompakt |

### Spektrum eines Operators

$$\sigma(T) = \{\lambda \in \mathbb{C} : (T - \lambda I)^{-1} \text{ existiert nicht oder unbeschränkt}\}$$

Zerlegung:
- **Punktspektrum** $\sigma_p(T)$: Eigenwerte ($Tv = \lambda v$, $v \neq 0$)
- **Stetiges Spektrum** $\sigma_c(T)$: $(T-\lambda I)^{-1}$ unbeschränkt
- **Residualspektrum** $\sigma_r(T)$: $(T-\lambda I)^{-1}$ nicht dicht definiert

### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `apply(v)` | $y = Tv$ |
| `is_bounded(vecs)` | Beschränktheitsprüfung |
| `operator_norm()` | $\|T\|$ |
| `is_compact()` | Kompaktheitsprüfung |
| `adjoint()` | Gibt $T^*$ zurück |
| `is_self_adjoint()` | Prüft $T = T^*$ |
| `is_unitary()` | Prüft $T^*T = I$ |
| `is_normal()` | Prüft $T^*T = TT^*$ |
| `is_projection()` | Prüft $T^2 = T$ |
| `kernel()` | Basisvektoren von $\ker(T)$ |
| `range_space()` | Basisvektoren von $\text{range}(T)$ |

---

## 5. Spektralsatz (`spectral_theorem_demo`)

**Spektralsatz (für kompakte selbstadjungierte Operatoren):**
Sei $T: H \to H$ kompakt und selbstadjungiert. Dann gibt es eine ONB $\{e_n\}$
aus Eigenvektoren und reelle Eigenwerte $\lambda_n$:

$$T = \sum_n \lambda_n \langle \cdot, e_n \rangle e_n$$

**Rekonstruktion:** $A = Q \Lambda Q^\top$ (Eigenwertzerlegung)

---

## 6. Fredholm-Alternative (`fredholm_alternative`)

**Satz (Fredholm 1903):**
$Ax = b$ ist lösbar $\Longleftrightarrow$ $b \perp \ker(A^*)$

$$\text{Fredholm-Index: } \text{ind}(A) = \dim \ker(A) - \dim \ker(A^*)$$

---

## 7. Funktionenräume

### C[a,b] – Stetige Funktionen

$$\|f\|_\infty = \max_{x \in [a,b]} |f(x)|$$

**Stone-Weierstraß-Satz:** Polynome liegen dicht in $C[a,b]$ bzgl. $\|\cdot\|_\infty$.

**Bernstein-Approximation:**
$$B_n f(x) = \sum_{k=0}^n f\!\left(\frac{k}{n}\right) \binom{n}{k} x^k (1-x)^{n-k}$$

### Sobolev-Räume W^{k,p}

$$W^{k,p}(\Omega) = \{f \in L^p : D^\alpha f \in L^p \text{ für } |\alpha| \leq k\}$$

$$\|f\|_{W^{k,p}} = \left(\sum_{|\alpha| \leq k} \|D^\alpha f\|_p^p\right)^{1/p}$$

**Einbettungssatz (Sobolev 1938):**
$$W^{k,p}(\Omega) \hookrightarrow C^m(\bar{\Omega}) \quad \text{wenn } k - \frac{n}{p} > m$$

$H^k = W^{k,2}$: Hilbert-Raum (zentral für schwache PDE-Formulierungen).

### L²[a,b]

$$L^2[a,b] = \left\{f \text{ messbar} : \int_a^b |f|^2 \, dx < \infty\right\}$$

$$\langle f, g \rangle = \int_a^b f(x) g(x) \, dx$$

**Fourier-ONB:**
$$e_0(x) = \frac{1}{\sqrt{b-a}}, \quad e_n(x) = \sqrt{\frac{2}{b-a}} \cos\!\left(\frac{n\pi(x-a)}{b-a}\right)$$

---

## 8. Klassische Sätze

### Hahn-Banach-Satz

Jedes beschränkte lineare Funktional $f: U \to \mathbb{R}$ auf einem Unterraum $U$
lässt sich auf den Gesamtraum $X$ fortsetzen ohne die Norm zu vergrößern:
$$\|F\|_{X^*} = \|f\|_{U^*}$$

### Satz von der offenen Abbildung (Banach-Schauder)

Surjektive beschränkte lineare Abbildung zwischen Banach-Räumen ist offen.
**Korollar:** Bijektive beschränkte Abbildung $\Rightarrow$ $T^{-1}$ beschränkt.

### Satz vom abgeschlossenen Graphen

$$\text{Graph}(T) \text{ abgeschlossen} \Longleftrightarrow T \text{ beschränkt}$$

### Gleichmäßiges Beschränktheitsprinzip (Banach-Steinhaus)

$$\sup_\alpha \|T_\alpha x\| < \infty \; \forall x \in X \Longrightarrow \sup_\alpha \|T_\alpha\| < \infty$$

### Schwache Konvergenz

$x_n \rightharpoonup x$ (schwach): $f(x_n) \to f(x)$ für alle $f \in X^*$

**Gegenbeispiel in $\ell^2$:** Einheitsvektoren $e_n \rightharpoonup 0$ (schwach), aber $\|e_n\| = 1$ (keine starke Konvergenz).

---

## 9. Fixpunktsätze

### Banach-Fixpunktsatz (1922)

Sei $T: X \to X$ eine **Kontraktion** ($d(Tx, Ty) \leq q \cdot d(x,y)$, $q < 1$)
auf einem vollständigen metrischen Raum. Dann:

1. $T$ hat genau einen Fixpunkt $x^*$
2. $T^n(x_0) \to x^*$ für jeden Startpunkt

**A-priori-Schranke:**
$$d(x_n, x^*) \leq \frac{q^n}{1-q} \cdot d(x_1, x_0)$$

### Schauder-Fixpunktsatz (1930)

Sei $K \subseteq X$ nichtleer, konvex, kompakt. Sei $T: K \to K$ stetig.
Dann hat $T$ mindestens einen Fixpunkt.

---

## Tests

Die Datei `tests/test_functional_analysis.py` enthält **117 Tests** in 11 Klassen:

| Klasse | Anzahl | Bereich |
|--------|--------|---------|
| `TestNormedSpace` | 11 | Normaxiome, Einheitsball |
| `TestBanachSpace` | 9 | Cauchy-Folgen, ℓᵖ-Normen |
| `TestHilbertSpace` | 18 | Skalarprodukt, Gram-Schmidt, Parseval |
| `TestLinearOperator` | 16 | Anwendung, Norm, Adjungierte |
| `TestCompactOperator` | 3 | Spektralzerlegung |
| `TestSpectralTheory` | 6 | Spektrum, Spektralsatz |
| `TestFredholmAlternative` | 4 | Lösbarkeit, Index |
| `TestContinuousFunctions` | 4 | C[a,b], Stone-Weierstraß |
| `TestSobolevSpace` | 4 | Sobolev-Norm, Einbettung |
| `TestL2Space` | 8 | L²-Skalarprodukt, Fourier-ONB |
| `TestClassicalTheorems` | 8 | Hahn-Banach, Offene Abbildung |
| `TestAdvancedTheorems` | 7 | Riesz, Schwache Konvergenz |
| `TestFixedPointTheorems` | 9 | Banach, Schauder |
| `TestEdgeCases` | 8 | Randfälle, Edge Cases |

**Alle 117 Tests grün.**
