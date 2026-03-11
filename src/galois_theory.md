# galois_theory.py – Galois-Theorie Dokumentation

**Autor:** Michael Fuhrmann
**Letzte Änderung:** 2026-03-11
**Build:** 84

---

## Übersicht

Das Modul `galois_theory.py` implementiert die vollständige Galois-Theorie –
eines der tiefgründigsten Gebiete der modernen Algebra. Es verbindet Körpertheorie,
Gruppentheorie und Polynomtheorie zu einem mächtigen Werkzeug.

**Kernaussage:**
Ein Polynom $f(x) \in \mathbb{Q}[x]$ ist genau dann durch Radikale lösbar, wenn
seine Galoisgruppe $\text{Gal}(f)$ eine **auflösbare** Gruppe ist.

---

## Mathematische Grundlagen

### 1. Körpererweiterungen

Eine **Körpererweiterung** $L/K$ ist eine Inklusionsbeziehung von Körpern $K \subseteq L$.

**Algebraisches Element:** $\alpha \in L$ ist algebraisch über $K$, wenn ein
Polynom $0 \neq f \in K[x]$ existiert mit $f(\alpha) = 0$.

**Minimalpolynom:** Das normierte Polynom minimalen Grads in $K[x]$ mit Nullstelle $\alpha$.

**Erweiterungsgrad:**
$$[L:K] = \dim_K L$$

**Zerfällungskörper:** Der kleinste Körper, über dem $f$ vollständig in Linearfaktoren zerfällt.

---

### 2. Galois-Gruppe

Die **Galois-Gruppe** einer Erweiterung $L/K$ ist:
$$\text{Gal}(L/K) = \text{Aut}_K(L) = \{\sigma : L \to L \mid \sigma \text{ Körperaut.}, \sigma|_K = \text{id}\}$$

**Klassifikation nach Grad:**

| Grad | Galoisgruppe | Auflösbar |
|------|-------------|-----------|
| 1 | trivial $\{e\}$ | Ja |
| 2 | $\mathbb{Z}/2\mathbb{Z}$ | Ja |
| 3 | $\mathbb{Z}/3\mathbb{Z}$ (Disk. QR) oder $S_3$ | Ja |
| 4 | $A_4, S_4, V_4, \mathbb{Z}/4\mathbb{Z}, D_4$ | Ja |
| 5+ | meist $S_n$ oder $A_n$ | Nein (für $n \geq 5$) |

**Diskriminante:**
$$\Delta(f) = a_n^{2n-2} \prod_{i < j} (\alpha_i - \alpha_j)^2$$

Für Grad 3: $\text{Gal}(f) = A_3$ gdw. $\sqrt{\Delta} \in \mathbb{Q}$.

---

### 3. Galois-Korrespondenz (Hauptsatz)

**Satz (Galois, ~1830):** Sei $L/K$ galoissch. Dann gibt es eine ordnungsumkehrende Bijektion:

$$\{\text{Untergruppen } H \leq \text{Gal}(L/K)\} \xleftrightarrow{1:1} \{\text{Zwischenkörper } K \subseteq F \subseteq L\}$$

**Gradformeln:**
- $[L : L^H] = |H|$
- $[L^H : K] = [G : H] = |G|/|H|$
- $|G| = [L : K]$

**Normalität:**
$H \trianglelefteq G \iff L^H / K$ galoissch, und dann $\text{Gal}(L^H/K) \cong G/H$.

---

### 4. Auflösbarkeit durch Radikale

**Definition (Radikalturm):**
$$K = K_0 \subset K_1 \subset \cdots \subset K_r$$
wobei $K_{i+1} = K_i(\alpha_i)$ mit $\alpha_i^{n_i} \in K_i$.

**Galois-Kriterium:**
$$f \text{ durch Radikale lösbar} \iff \text{Gal}(f) \text{ ist auflösbar}$$

**Abel-Ruffini (1799/1824):**
Es gibt keine allgemeine Lösungsformel für Polynome $n \geq 5$, da $S_5$ nicht auflösbar ist.

**Auflösungsformeln nach Grad:**
- Grad 1: $x = -b/a$
- Grad 2: $x = \dfrac{-b \pm \sqrt{b^2-4ac}}{2a}$ (Quadratische Formel)
- Grad 3: Cardano-Formel via $\sqrt[3]{\cdot}$
- Grad 4: Ferrari-Methode via Resolventenkubik

---

### 5. Zyklotomische Körper

Der **n-te Kreisteilungskörper** $\mathbb{Q}(\zeta_n)$ mit $\zeta_n = e^{2\pi i/n}$:

$$\text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \cong (\mathbb{Z}/n\mathbb{Z})^\times$$

**Ordnung:** $[\mathbb{Q}(\zeta_n):\mathbb{Q}] = \varphi(n)$ (Euler-Phi-Funktion)

**Kreisteilungspolynom:**
$$\Phi_n(x) = \prod_{\substack{k=1 \\ \gcd(k,n)=1}}^{n} (x - \zeta_n^k)$$

**Beispiele:**
| $n$ | $\Phi_n(x)$ | Gal-Gruppe |
|-----|-------------|------------|
| 3 | $x^2+x+1$ | $\mathbb{Z}/2\mathbb{Z}$ |
| 4 | $x^2+1$ | $\mathbb{Z}/2\mathbb{Z}$ |
| 5 | $x^4+x^3+x^2+x+1$ | $\mathbb{Z}/4\mathbb{Z}$ |
| 12 | $x^4-x^2+1$ | $\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/2\mathbb{Z}$ |

---

### 6. Kronecker-Weber-Theorem

**Satz (Kronecker 1853, Weber 1886):**
Jede endliche abelsche Galois-Erweiterung $L/\mathbb{Q}$ ist in einem
Kreisteilungskörper enthalten:
$$L \subseteq \mathbb{Q}(\zeta_n) \text{ für ein } n \in \mathbb{N}$$

---

### 7. Endliche Körper

Der eindeutige (bis auf Isomorphie) endliche Körper mit $q = p^n$ Elementen:

$$\text{GF}(p^n) = \mathbb{F}_p[x]/(f(x))$$

wobei $f$ irreduzibel vom Grad $n$ über $\mathbb{F}_p$.

**Galoisgruppe:**
$$\text{Gal}(\text{GF}(p^n)/\text{GF}(p)) \cong \mathbb{Z}/n\mathbb{Z}$$

**Frobenius-Automorphismus:** $\text{Frob}_p : x \mapsto x^p$ (Erzeuger der Galoisgruppe)

---

### 8. Konstruierbarkeit (Gauss-Wantzel)

**Satz (Gauss 1796, Wantzel 1837):**
Ein reguläres $n$-Eck ist mit Zirkel und Lineal konstruierbar genau dann, wenn:
$$n = 2^k \cdot p_1 \cdot p_2 \cdots p_r$$
mit verschiedenen **Fermat-Primzahlen** $p_i = 2^{2^m}+1$.

**Bekannte Fermat-Primzahlen:** $F_0 = 3$, $F_1 = 5$, $F_2 = 17$, $F_3 = 257$, $F_4 = 65537$.

**Galois-Begründung:** $n$-Eck konstruierbar $\iff \varphi(n) = 2^k$ (2-Potenz).

| $n$ | Konstruierbar? | Warum |
|-----|---------------|-------|
| 3 | Ja | $3 = F_0$ |
| 5 | Ja | $5 = F_1$ |
| 7 | Nein | 7 keine Fermat-Prim |
| 9 | Nein | $9 = 3^2$ (Exponent > 1) |
| 17 | Ja | $17 = F_2$ (Gauss 1796) |
| 15 | Ja | $15 = 3 \cdot 5$ (beide Fermat-Prims) |

---

### 9. Hilbert's Satz 90

**Originalform (Hilbert, 1897):**
Sei $L/K$ zyklisch galoissch mit $\text{Gal}(L/K) = \langle\sigma\rangle$. Dann:
$$N_{L/K}(\alpha) = 1 \iff \exists \beta \in L^\times : \alpha = \frac{\beta}{\sigma(\beta)}$$

**Kohomologische Form:**
$$H^1(\text{Gal}(L/K),\ L^\times) = 1$$

**Anwendungen:**
- Charakterisierung zyklischer Erweiterungen (Kummer-Theorie)
- Grundlage der Galois-Kohomologie
- Beschreibung von Lösbarkeit von Normengleichungen

---

### 10. Norm und Spur

Für $\alpha \in L = K(\alpha)$ mit Minimalpolynom $f(x) = a_n x^n + \cdots + a_0$:

**Norm:**
$$N_{L/K}(\alpha) = \prod_{\sigma \in \text{Gal}(L/K)} \sigma(\alpha) = (-1)^n \cdot \frac{a_0}{a_n}$$

**Spur:**
$$\text{Tr}_{L/K}(\alpha) = \sum_{\sigma \in \text{Gal}(L/K)} \sigma(\alpha) = -\frac{a_{n-1}}{a_n}$$

---

## Klassen

### `FieldExtension`

Körpererweiterung $L = K(\alpha)$ mit Minimalpolynom $f$.

```python
ext = FieldExtension(base_poly=[0, 1], min_poly=[-2, 0, 1])
# Repräsentiert Q(sqrt(2)) mit Minimalpolynom x^2 - 2

ext.degree()           # [L:K] = 2
ext.is_algebraic()     # True
ext.minimal_polynomial()  # [-2, 0, 1]
ext.norm()             # N_{L/K}(alpha)
ext.trace()            # Tr_{L/K}(alpha)
ext.is_galois()        # True/False
```

### `GaloisGroup`

Galoisgruppe $\text{Gal}(L/K)$ als abstrakte Gruppe.

```python
gal = GaloisGroup(poly_coeffs=[-2, 0, 1])

gal.order()            # |Gal(L/K)| = 2
gal.elements()         # ['e', 'sigma']
gal.is_abelian()       # True
gal.is_solvable()      # True
gal.subgroups()        # [{'name': 'G', ...}, {'name': '{e}', ...}]
gal.fixed_field('G')   # 'K (Basiskörper)'
```

### `FiniteField`

Endlicher Körper $\text{GF}(p^n)$.

```python
gf4 = FiniteField(p=2, n=2)   # GF(4) = GF(2^2)

gf4.order()                    # 4
gf4.multiplicative_group()     # Elemente von GF(4)*
gf4.frobenius(element)         # x -> x^2
gf4.primitive_elements()       # Erzeuger von GF(4)*
```

---

## Funktionen (Build 84 – neu hinzugefügt)

### `cyclotomic_galois_group(n)`

Berechnet $\text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \cong (\mathbb{Z}/n\mathbb{Z})^\times$.

```python
result = cyclotomic_galois_group(5)
# {'n': 5, 'order': 4, 'galois_group': 'Z/4Z',
#  'elements': [1,2,3,4], 'is_abelian': True}
```

### `kronecker_weber_check(poly_coeffs)`

Prüft, ob Erweiterung abelsch ist und gibt minimalen Kreisteilungskörper an.

```python
result = kronecker_weber_check([-2, 0, 1])
# {'is_abelian': True, 'minimal_n': 8,
#  'cyclotomic_field': 'Q(zeta_8)'}
```

### `dirichlet_characters_from_galois(n)`

Leitet Dirichlet-Charaktere mod $n$ aus der Galoisgruppe ab.

### `galois_group_finite_field(p, n)`

Berechnet $\text{Gal}(\text{GF}(p^n)/\text{GF}(p)) \cong \mathbb{Z}/n\mathbb{Z}$.

```python
result = galois_group_finite_field(2, 4)
# {'galois_group': 'Z/4Z', 'order': 4,
#  'generator': 'Frob_2: x -> x^2'}
```

### `is_constructible(n)` / `construct_regular_polygon(n)`

```python
is_constructible(17)  # True  (Gauss 1796)
is_constructible(7)   # False

polygon = construct_regular_polygon(17)
# {'is_constructible': True, 'vertices': [...17 Punkte...]}
```

### `norm_and_trace(alpha_coeffs, poly_coeffs)`

Berechnet $N_{L/K}(\alpha)$ und $\text{Tr}_{L/K}(\alpha)$ via Vieta-Formeln.

```python
result = norm_and_trace([], [-2, 0, 1])
# {'norm': -2.0, 'trace': 0.0, 'degree': 2}
```

### `hilbert90(n)`

Demonstration von $H^1(\text{Gal}(L/K), L^\times) = 1$.

### `galois_group_symmetric(n)`

Beschreibt $S_n$ als Galoisgruppe des generischen Grad-$n$-Polynoms.

```python
result = galois_group_symmetric(5)
# {'group': 'S_5', 'order': 120, 'is_solvable': False}
```

### `radical_tower(poly_coeffs)`

Konstruiert expliziten Radikalturm für auflösbare Polynome.

### `abel_ruffini_demo()`

Demonstration des Satzes von Abel-Ruffini (Unlösbarkeit von Grad $\geq 5$).

```python
result = abel_ruffini_demo()
# {'s5_solvable': False, 'a5_simple': True, 'a5_order': 60,
#  'solvable_by_degree': {1: ..., 5: False}}
```

### `fundamental_theorem_verify(poly_coeffs)`

Verifiziert den Hauptsatz der Galois-Theorie konkret.

### `galois_group_of_polynomial(f)`

Vollständige Galoisgruppen-Information mit Kompositionsreihe und Implikationen.

---

## Tests

Die Test-Datei `tests/test_galois_theory.py` enthält **184 Tests** in 19 Klassen.

Alle Tests ausführen:
```bash
python3 -m pytest tests/test_galois_theory.py -v
```

| Testklasse | Anzahl | Thema |
|------------|--------|-------|
| `TestDiscriminantPolynomial` | 8 | Diskriminante |
| `TestGaloisGroupPolynomial` | 8 | Galoisgruppe nach Grad |
| `TestIsSolvableByRadicals` | 7 | Auflösbarkeit |
| `TestCyclotomicPolynomial` | 9 | Kreisteilungspolynome |
| `TestSplittingField` | 6 | Zerfällungskörper |
| `TestFiniteField` | 11 | Endliche Körper |
| `TestFiniteFieldDiscreteLog` | 6 | Diskreter Logarithmus |
| `TestMinimalPolynomial` | 3 | Minimalpolynome |
| `TestFieldExtension` | 9 | Körpererweiterungen |
| `TestGaloisGroup` | 8 | Galoisgruppe-Klasse |
| `TestGaloisCorrespondence` | 5 | Galois-Korrespondenz |
| `TestPrimitiveElementTheorem` | 3 | Primitives Element |
| `TestHelperFunctions` | 5 | Hilfsfunktionen |
| `TestCyclotomicGaloisGroup` | 8 | NEU: Kreisteilungs-Galoisgruppe |
| `TestKroneckerWeberCheck` | 6 | NEU: Kronecker-Weber |
| `TestDirichletCharactersFromGalois` | 6 | NEU: Dirichlet-Charaktere |
| `TestGaloisGroupFiniteField` | 8 | NEU: Galoisgruppe endl. Körper |
| `TestIsConstructible` | 10 | NEU: Konstruierbarkeit |
| `TestConstructRegularPolygon` | 6 | NEU: Reguläre Polygone |
| `TestNormAndTrace` | 8 | NEU: Norm und Spur |
| `TestHilbert90` | 6 | NEU: Hilbert 90 |
| `TestGaloisGroupSymmetric` | 8 | NEU: S_n |
| `TestRadicalTower` | 7 | NEU: Radikalturm |
| `TestAbelRuffiniDemo` | 9 | NEU: Abel-Ruffini |
| `TestFundamentalTheoremVerify` | 5 | NEU: Hauptsatz-Verifikation |
| `TestGaloisGroupOfPolynomial` | 7 | NEU: Erweiterte Galoisgruppen-Info |

---

## Referenzen

- Galois, É. (1846): *Oeuvres mathématiques d'Évariste Galois*
- Artin, E. (1944): *Galois Theory*, Notre Dame Math. Lectures
- Lang, S. (2002): *Algebra*, Springer GTM 211
- Stewart, I. (2015): *Galois Theory*, 4th ed., CRC Press
- Neukirch, J. (1999): *Algebraic Number Theory*, Springer
- Gauss, C.F. (1801): *Disquisitiones Arithmeticae*
