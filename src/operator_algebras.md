# operator_algebras.py — Operatoralgebren

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Abhängigkeiten:** `numpy`

---

## Überblick

Das Modul `operator_algebras.py` implementiert abstrakte Strukturen der Funktionalanalysis und nichtkommutativen Geometrie. Alle Operatoren werden durch endlichdimensionale komplexe Matrizen (NumPy-Arrays) dargestellt, was numerische Verifikation algebraischer Eigenschaften ermöglicht.

---

## 1. C*-Algebren

### 1.1 Axiome

Eine **C*-Algebra** ist ein Banach-Raum $A$ mit einer assoziativen Multiplikation und einer konjugierter Involution $a \mapsto a^*$, die folgende Axiome erfüllt:

| Axiom | Eigenschaft |
|-------|-------------|
| $(a^*)^* = a$ | Involution ist involutiv |
| $(ab)^* = b^*a^*$ | Anti-Multiplikativität |
| $\|a\| \geq 0$ | Norm positiv definit |
| $\|ab\| \leq \|a\|\cdot\|b\|$ | Submultiplikativität |
| $\|a^*a\| = \|a\|^2$ | **C*-Identität** (zentrales Axiom) |

### 1.2 Kanonisches Beispiel: $M_n(\mathbb{C})$

Die Algebra der $n\times n$-Matrizen mit:
- Multiplikation: Matrizenmultiplikation $a \cdot b = ab$
- Involution: $a^* = \bar{a}^T$ (konjugierte Transposition)
- Norm: $\|a\| = \sigma_{\max}(a)$ (größter Singulärwert)

**Verifikation der C*-Identität:**
$$\|a^*a\| = \sigma_{\max}(a^*a) = \sigma_{\max}(a)^2 = \|a\|^2 \checkmark$$

### 1.3 Spektrum

Das **Spektrum** eines Elements $a \in A$:
$$\sigma(a) = \{\lambda \in \mathbb{C} : (a - \lambda I) \text{ nicht invertierbar}\}$$

In endlichen Dimensionen: $\sigma(a) = $ Menge der Eigenwerte von $a$.

**Spektralradius:**
$$r(a) = \max\{|\lambda| : \lambda \in \sigma(a)\} = \lim_{n\to\infty}\|a^n\|^{1/n}$$

Für normale Elemente ($a^*a = aa^*$): $r(a) = \|a\|$.

---

## 2. Von-Neumann-Algebren

### 2.1 Definition

Eine **von-Neumann-Algebra** $M \subseteq B(H)$ ist eine schwach abgeschlossene $*$-Unteralgebra von $B(H)$, die die Einheit enthält.

### 2.2 Kommutant

Der **Kommutant** einer Menge $S \subseteq B(H)$:
$$S' = \{T \in B(H) : TS = ST \;\forall S \in \mathcal{S}\}$$

Eigenschaften:
- $S \subseteq S''$ (doppelter Kommutant enthält $S$)
- $S' = S'''$ (dreifacher = einfacher Kommutant)
- $S'$ ist immer eine von-Neumann-Algebra

### 2.3 Bikommutant-Satz (von Neumann, 1929)

$$M \text{ ist von-Neumann-Algebra} \iff M = M''$$

**Interpretation:** Eine $*$-Unteralgebra ist genau dann schwach abgeschlossen, wenn sie mit ihrem doppelten Kommutanten übereinstimmt.

**Numerische Umsetzung:** In endlichen Dimensionen wird das Kommutantensystem
$$[T, S] = TS - ST = 0$$
als lineares Gleichungssystem für $\text{vec}(T)$ gelöst:
$$\bigl[(S^T \otimes I) - (I \otimes S)\bigr]\,\text{vec}(T) = 0$$

### 2.4 Zentrum und Faktoren

Das **Zentrum** $Z(M) = M \cap M'$ besteht aus allen Elementen, die mit allen anderen kommutieren.

$M$ ist ein **Faktor** genau dann, wenn:
$$Z(M) = \{\lambda I : \lambda \in \mathbb{C}\}$$

---

## 3. Projektionen in $B(H)$

### 3.1 Definition

$P \in B(H)$ ist eine **orthogonale Projektion**, falls:
$$P = P^* = P^2$$

### 3.2 Ordnung auf Projektionen

$$P \leq Q \iff PQ = P \iff \text{Bild}(P) \subseteq \text{Bild}(Q)$$

### 3.3 Partiell-Isometrische Operatoren

$V \in B(H)$ ist eine **partielle Isometrie**, falls $V^*V$ eine Projektion ist (Initial-Projektion). Äquivalent: $VV^*V = V$.

Jede Isometrie ist partiell-isometrisch: Isometrie $\Rightarrow$ $V^*V = I$ (Projektion).

---

## 4. Gelfand-Isomorphismus

### 4.1 Satz (Gelfand, 1943)

Sei $A$ eine kommutative C*-Algebra mit Einheit. Dann ist:
$$A \cong C(X)$$
wobei $X = \hat{A}$ der **Gelfand-Raum** (Raum der multiplikativen Funktionale/Charaktere) ist.

**Gelfand-Transformation:**
$$\hat{a}: X \to \mathbb{C}, \quad \hat{a}(\chi) = \chi(a)$$

Isometrisch: $\|\hat{a}\|_{\sup} = \|a\|_A$.

### 4.2 Beispiel: Diagonalmatrizen

$$D_n = \{\text{diag}(d_1,\ldots,d_n)\} \subset M_n(\mathbb{C}) \cong C(\{1,\ldots,n\}) = \mathbb{C}^n$$

Charaktere: $\chi_k(\text{diag}(d_1,\ldots,d_n)) = d_k$ für $k = 1,\ldots,n$.

Gelfand-Transformation: $\widehat{\text{diag}(d)} = (d_1,\ldots,d_n)$.

Isometrie: $\sup_k |d_k| = \|\text{diag}(d)\|_{\rm op}$. ✓

---

## 5. GNS-Konstruktion (Gelfand-Naimark-Segal)

### 5.1 Satz

Sei $A$ eine C*-Algebra und $\phi: A \to \mathbb{C}$ ein **Zustand** (positives lineares Funktional mit $\phi(I) = 1$). Dann existiert eine Hilbert-Raum-Darstellung:
$$\pi: A \to B(H_\phi)$$

### 5.2 Konstruktion

1. Prä-Hilbert-Raum: $(A, \langle a, b \rangle_\phi) = \phi(a^* b)$
2. Quotienten: $H_\phi = A / N_\phi$ mit $N_\phi = \{a : \phi(a^*a) = 0\}$ (Vervollständigung)
3. Darstellung: $\pi(a)[b] = [ab]$ (Linksmultiplikation)
4. Zyklischer Vektor: $\Omega = [I] \in H_\phi$

### 5.3 Beispiel: $M_2(\mathbb{C})$ mit normierter Spur

Zustand: $\phi(a) = \tfrac{1}{2}\text{Tr}(a)$

Hilbert-Schmidt-Skalarprodukt:
$$\langle A, B \rangle = \tfrac{1}{2}\text{Tr}(A^*B)$$

Isometrie: $\|\pi(a)\|_{B(H)} = \|a\|_{M_2}$ (GNS-Darstellung ist isometrisch für injektiven Zustand).

---

## 6. Murray-von Neumann Klassifikation der Faktoren

| Typ | Charakterisierung | Beispiel |
|-----|-------------------|---------|
| $I_n$ | Endlichdimensional, $n < \infty$ | $M_n(\mathbb{C})$ |
| $I_\infty$ | $B(H)$, separabel | $B(\ell^2(\mathbb{N}))$ |
| $II_1$ | Endlicher Faktor, normierte Spur: $\text{Tr}(I) = 1$ | Hyperendlicher Faktor $R$ |
| $II_\infty$ | Halbbeschränkt: $II_1 \otimes I_\infty$ | $R \otimes B(\ell^2)$ |
| $III$ | Rein unendlich, keine halbendliche Spur | Powers-Faktoren $R_\lambda$ |

**Connes-Klassifikation** (1973–76) der Typ-III-Faktoren:
- $III_0$: keine weitere Invariante (Krieger-Faktoren)
- $III_\lambda$ ($0 < \lambda < 1$): Spektrum des Modularoperators $\{0\} \cup \lambda^{\mathbb{Z}}$
- $III_1$: Spektrum $= [0,\infty)$

---

## 7. K-Theorie für C*-Algebren

### 7.1 $K_0$-Gruppe

Die Grothendieck-Gruppe der stabilen Äquivalenzklassen von Projektionen:

- Zwei Projektionen $P, Q \in M_\infty(A)$ sind **stabil äquivalent**, falls $P \oplus 0 \sim Q \oplus 0$ (unitäre Äquivalenz).
- $K_0(A) = $ Grothendieck-Gruppe dieser Halbgruppe.

**Beispiele:**
$$K_0(M_n(\mathbb{C})) \cong \mathbb{Z}, \quad K_0(C(X)) \cong K^0(X)$$

### 7.2 $K_1$-Gruppe

Unitäre Äquivalenzklassen in $GL_\infty(A)^+ = \varinjlim GL_n(A)$:
$$K_1(A) = GL_\infty(A) / [GL_\infty(A), GL_\infty(A)]$$

**Beispiele:**
$$K_1(M_n(\mathbb{C})) = 0, \quad K_1(C(S^1)) \cong \mathbb{Z} \text{ (Windungszahl)}$$

### 7.3 Bott-Periodizität

$$K_{n+2}(A) \cong K_n(A) \quad \forall n \geq 0$$

**Sechsterm-Folge** (für kurze exakte Sequenz $0 \to I \to A \to A/I \to 0$):
$$K_0(I) \to K_0(A) \to K_0(A/I) \xrightarrow{\partial} K_1(I) \to K_1(A) \to K_1(A/I) \xrightarrow{\partial} K_0(I)$$

---

## 8. Cuntz-Algebren $\mathcal{O}_n$

### 8.1 Definition (Cuntz, 1977)

$\mathcal{O}_n$ ist die universelle C*-Algebra, erzeugt von $n$ Isometrien $S_1,\ldots,S_n$ mit:

$$S_i^* S_j = \delta_{ij} I, \quad \sum_{i=1}^n S_i S_i^* = I$$

### 8.2 Eigenschaften

- **Einfach**: kein nichttriviales abgeschlossenes zweiseitiges Ideal
- **Rein unendlich**: jede Projektion $\neq 0$ ist äquivalent zu $I$ (Murray-von Neumann)
- **Typ $III_1$-Faktor** in der GNS-Darstellung bzgl. der Spurzustands (für $n \geq 2$)

### 8.3 K-Theorie

$$K_0(\mathcal{O}_n) \cong \mathbb{Z}/(n-1)\mathbb{Z}, \quad K_1(\mathcal{O}_n) = 0$$

### 8.4 Toy-Modell in $M_4(\mathbb{C})$

Isometrien $S_1, S_2: \mathbb{C}^2 \to \mathbb{C}^4$ mit $S_1 S_1^* + S_2 S_2^* = I_4$ und $S_i^* S_j = \delta_{ij} I_2$.

---

## 9. Spurklasse-Operatoren $S_1(H)$

### 9.1 Definition

$T \in B(H)$ ist in der **Spurklasse** $S_1(H)$, falls:
$$\|T\|_1 = \text{Tr}(|T|) = \sum_{i} \sigma_i(T) < \infty$$

wobei $|T| = \sqrt{T^*T}$ und $\sigma_i$ die Singulärwerte sind.

### 9.2 Eigenschaften

| Eigenschaft | Formel |
|-------------|--------|
| Spurformel | $\text{Tr}(T) = \sum_n \langle Te_n, e_n\rangle$ (absolut konvergent) |
| Zyklizität | $\text{Tr}(AB) = \text{Tr}(BA)$ |
| Spurklassen-Norm | $\|T\|_1 = \sum \sigma_i$ |
| Dualität | $S_1(H)^* \cong B(H)$ |

### 9.3 Schatten-$p$-Klassen

$$S_p(H) = \Bigl\{T \in K(H) : \|T\|_p = \Bigl(\sum_i \sigma_i^p\Bigr)^{1/p} < \infty\Bigr\}$$

Inklusion: $S_1 \subset S_2 \subset \cdots \subset K(H) \subset B(H)$

Norm-Ungleichungen: $\|T\| \leq \|T\|_2 \leq \|T\|_1$

- $S_1$: Spurklasse-Operatoren
- $S_2$: Hilbert-Schmidt-Operatoren ($\|T\|_2^2 = \text{Tr}(T^*T) < \infty$)
- $S_\infty = K(H)$: kompakte Operatoren
