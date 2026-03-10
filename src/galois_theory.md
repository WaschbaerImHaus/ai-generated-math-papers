# galois_theory.py – Galois-Theorie und Körpererweiterungen

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Build:** 28

---

## Übersicht

Das Modul `galois_theory.py` implementiert die grundlegenden Konzepte der **Galois-Theorie** –
einem zentralen Zweig der modernen Algebra, der Körpererweiterungen und Polynomgleichungen
mittels Gruppentheorie untersucht.

Die Galois-Theorie erklärt, warum es für Polynome vom Grad ≥ 5 **keine allgemeine Radikallösung**
gibt (Abel-Ruffini-Theorem, 1824) und führt die Lösbarkeit von Polynomgleichungen auf
gruppentheoretische Eigenschaften zurück.

---

## Klassen

### `FieldExtension`

Modelliert eine algebraische Körpererweiterung $L = K(\alpha)$ über $K = \mathbb{Q}$.

Das Element $\alpha$ wird durch sein **Minimalpolynom** $f \in \mathbb{Q}[x]$ charakterisiert.

**Wichtige Konzepte:**
- **Grad** $[L:K] = \deg(f)$: Dimension von $L$ als $K$-Vektorraum
- **Basis**: $\{1, \alpha, \alpha^2, \ldots, \alpha^{n-1}\}$ ist $K$-Basis von $L$
- **Norm**: $N_{L/K}(\alpha) = (-1)^n \cdot \frac{a_0}{a_n}$ (Produkt aller Konjugierten)
- **Spur**: $\text{Tr}_{L/K}(\alpha) = -\frac{a_{n-1}}{a_n}$ (Summe aller Konjugierten)

**Methoden:**
| Methode | Beschreibung |
|---------|-------------|
| `degree()` | Körpergrad $[L:K]$ |
| `is_algebraic()` | Immer True (per Konstruktion) |
| `minimal_polynomial()` | Koeffizientenliste des Minimalpolynoms |
| `basis()` | $K$-Basis von $L$ |
| `norm()` | Norm $N_{L/K}(\alpha)$ |
| `trace()` | Spur $\text{Tr}_{L/K}(\alpha)$ |
| `is_separable()` | Separabilität (über $\mathbb{Q}$ stets True) |
| `is_normal()` | Normalität der Erweiterung |
| `is_galois()` | Galois-Eigenschaft (normal + separabel) |

**Beispiel:**
```python
ext = FieldExtension([-2, 0, 1])  # x² - 2 → ℚ(√2)
ext.degree()   # → 2
ext.norm()     # → -2  (√2 · (-√2) = -2)
ext.trace()    # → 0   (√2 + (-√2) = 0)
```

---

### `GaloisGroup`

Repräsentiert die Galoisgruppe $\text{Gal}(L/K) = \text{Aut}_K(L)$ eines Polynoms.

Die Galoisgruppe operiert auf den Wurzeln des Polynoms als **Untergruppe der symmetrischen Gruppe** $S_n$.

**Galoisgruppen nach Polynomgrad:**
| Grad | Mögliche Galoisgruppen |
|------|----------------------|
| 1 | Trivial $\{1\}$ |
| 2 | $\mathbb{Z}/2\mathbb{Z}$ oder trivial |
| 3 | $\mathbb{Z}/3\mathbb{Z}$ (wenn $\Delta$ Quadrat) oder $S_3$ |
| 4 | $S_4, A_4, D_4, \mathbb{Z}/4\mathbb{Z}, V_4 = \mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/2\mathbb{Z}$ |
| $\geq 5$ | Generisch $S_n$ (nicht auflösbar für $n \geq 5$) |

**Methoden:**
| Methode | Beschreibung |
|---------|-------------|
| `order()` | Gruppenordnung |
| `elements()` | Gruppenelemente (als Strings) |
| `is_abelian()` | Abelsch? |
| `is_solvable()` | Auflösbar? |
| `subgroups()` | Alle Untergruppen |
| `normal_subgroups()` | Alle Normalteiler |
| `fixed_field(subgroup_name)` | Fixkörper $L^H$ |

---

### `FiniteField`

Implementiert den endlichen Körper $\text{GF}(p^n) = \mathbb{Z}_p[x]/(f(x))$
mit $f$ irreduzibel vom Grad $n$.

**Mathematische Grundlage:**
- $|\text{GF}(p^n)| = p^n$ (Körperordnung)
- $\text{GF}(p^n)^\times$ ist **zyklisch** der Ordnung $p^n - 1$
- **Frobenius-Endomorphismus**: $\varphi: a \mapsto a^p$ (erzeugt $\text{Gal}(\text{GF}(p^n)/\text{GF}(p)) \cong \mathbb{Z}/n\mathbb{Z}$)

**Methoden:**
| Methode | Beschreibung |
|---------|-------------|
| `order()` | Körperordnung $p^n$ |
| `multiplicative_group()` | Alle Nicht-Null-Elemente |
| `frobenius_endomorphism(element)` | $\varphi(a) = a^p$ |
| `is_primitive_element(element)` | Primitives Element? |
| `primitive_elements()` | Alle primitiven Elemente |

---

## Hauptfunktionen

### `discriminant_polynomial(poly_coeffs)`

Berechnet die Diskriminante $\Delta(f)$ eines Polynoms:

$$\Delta(f) = (-1)^{n(n-1)/2} \cdot \frac{\text{Res}(f, f')}{a_n}$$

**Bedeutung:**
- $\Delta > 0$: gerade Anzahl von Paaren komplexer Wurzeln
- $\Delta < 0$: ungerade Anzahl von Paaren komplexer Wurzeln
- $\Delta = 0$: $f$ hat eine Mehrfachwurzel

**Beispiele:**
```python
discriminant_polynomial([-2, 0, 1])   # x²-2 → 8
discriminant_polynomial([-2, 0, 0, 1]) # x³-2 → -108
```

---

### `galois_group_polynomial(poly_coeffs)`

Berechnet die Galoisgruppe von $f(x) \in \mathbb{Q}[x]$.

**Algorithmus:**
1. **Grad 2**: Diskriminante $\Delta = b^2 - 4ac$. Quadrat → trivial, sonst $\mathbb{Z}/2\mathbb{Z}$.
2. **Grad 3**: $\Delta$ Quadrat → $\mathbb{Z}/3\mathbb{Z}$, sonst $S_3$.
3. **Grad 4**: Kubische Resolvente und Diskriminante.
4. **Grad $\geq 5$**: Generisch $S_n$ (nicht auflösbar).

**Rückgabe:**
```python
{
    'galois_group': 'S_3',     # Name der Gruppe
    'order': 6,                 # Gruppenordnung
    'is_solvable': True,        # Auflösbarkeit
    'discriminant': -108        # Diskriminante
}
```

---

### `is_solvable_by_radicals(poly_coeffs)`

Prüft ob $f(x)$ durch Radikale lösbar ist (**Galois-Hauptsatz**):

$$f(x) \text{ durch Radikale lösbar} \iff \text{Gal}(f) \text{ auflösbar}$$

**Bekannte Ergebnisse:**
- Grad $\leq 4$: **Stets auflösbar** (Cardano/Ferrari)
- Grad $\geq 5$: Generisch **nicht auflösbar** (Abel-Ruffini)

```python
is_solvable_by_radicals([-2, 0, 0, 1])  # x³-2
# → {'solvable': True, 'galois_group': 'S_3', 'reason': 'Cardano-Formel...'}

is_solvable_by_radicals([12, -5, 0, 0, 0, 1])  # x⁵-5x+12 (Gal=S_5)
# → {'solvable': False, 'galois_group': 'S_5', 'reason': 'Abel-Ruffini...'}
```

---

### `cyclotomic_polynomial(n)`

Berechnet das $n$-te **Kreisteilungspolynom**:

$$\Phi_n(x) = \prod_{\substack{1 \leq k \leq n \\ \gcd(k,n)=1}} \left(x - e^{2\pi i k/n}\right)$$

- Grad = $\varphi(n)$ (Euler-Phi)
- Irreduzibel über $\mathbb{Q}$
- $\text{Gal}(\mathbb{Q}(\zeta_n)/\mathbb{Q}) \cong (\mathbb{Z}/n\mathbb{Z})^\times$

```python
cyclotomic_polynomial(1)   # → [-1, 1]      (x - 1)
cyclotomic_polynomial(2)   # → [1, 1]       (x + 1)
cyclotomic_polynomial(3)   # → [1, 1, 1]    (x² + x + 1)
cyclotomic_polynomial(4)   # → [1, 0, 1]    (x² + 1)
cyclotomic_polynomial(6)   # → [1, -1, 1]   (x² - x + 1)
cyclotomic_polynomial(12)  # → [1, 0, -1, 0, 1]  (x⁴ - x² + 1)
```

---

### `splitting_field(poly_coeffs)`

Berechnet den Zerfällungskörper $K_f$ von $f$ über $\mathbb{Q}$:
der kleinste Körper, über dem $f$ vollständig in Linearfaktoren zerfällt.

$$[K_f : \mathbb{Q}] = |\text{Gal}(f)|$$

---

### `finite_field_discrete_log(a, g, q)`

Löst $g^x \equiv a \pmod{q}$ mittels **Baby-Step-Giant-Step** (Shanks):

1. $m = \lceil\sqrt{q-1}\rceil$
2. Baby steps: Tabelle $\{g^j \bmod q : j\}$ für $j = 0, \ldots, m-1$
3. Giant steps: $a \cdot (g^{-m})^i$ für $i = 0, \ldots, m$
4. Kollision: $x = j + m \cdot i$

**Laufzeit:** $O(\sqrt{q})$

---

### `galois_correspondence(poly_coeffs)`

Implementiert den **Hauptsatz der Galois-Theorie**:

$$\{\text{Untergruppen } H \leq \text{Gal}(L/K)\} \xleftrightarrow{1:1} \{\text{Zwischenkörper } K \subseteq F \subseteq L\}$$

Die Bijektion ist **ordnungsumkehrend**:
- $H_1 \subseteq H_2 \iff L^{H_1} \supseteq L^{H_2}$
- $[L : L^H] = |H|$, $[L^H : K] = [G:H]$

---

### `primitive_element_theorem(alpha, beta)`

**Satz vom primitiven Element:**
$$\mathbb{Q}(\alpha, \beta) = \mathbb{Q}(\alpha + c \cdot \beta) \text{ für fast alle } c \in \mathbb{Q}$$

---

## Algorithmische Details

### Irreduzibilitätstest (Rabin-Test)

Für $f \in \mathbb{F}_p[x]$ vom Grad $n$ gilt:

$f$ ist irreduzibel $\iff$
1. $\gcd(f, x^{p^{n/q}} - x) = 1$ für alle Primteiler $q$ von $n$
2. $f \mid x^{p^n} - x$

### Polynomoperationen in $\mathbb{F}_p[x]/(f)$

Alle Arithmetikoperationen (Multiplikation, Potenzierung) verwenden **Square-and-Multiply** mit Polynomreduktion modulo $f$ und Koeffizienten modulo $p$.

---

## Testabdeckung

Die Datei `tests/test_galois_theory.py` enthält **90 Tests** (alle grün):

| Testklasse | Anzahl | Beschreibung |
|-----------|--------|-------------|
| `TestDiscriminantPolynomial` | 8 | Diskriminanten für Grad 1-3 |
| `TestGaloisGroupPolynomial` | 8 | Galoisgruppen Grad 1-5 |
| `TestIsSolvableByRadicals` | 7 | Auflösbarkeit Grad 1-5 |
| `TestCyclotomicPolynomial` | 9 | Φ_n für n=1..12 |
| `TestSplittingField` | 6 | Zerfällungskörper |
| `TestFiniteField` | 11 | GF(p^n) Eigenschaften |
| `TestFiniteFieldDiscreteLog` | 6 | BSGS-Algorithmus |
| `TestMinimalPolynomial` | 3 | Minimalpolynome |
| `TestFieldExtension` | 10 | Körpererweiterungen |
| `TestGaloisGroup` | 9 | Galoisgruppen-Klasse |
| `TestGaloisCorrespondence` | 5 | Hauptsatz |
| `TestPrimitiveElementTheorem` | 3 | Primitives Element |
| `TestHelperFunctions` | 5 | Interne Hilfsfunktionen |
