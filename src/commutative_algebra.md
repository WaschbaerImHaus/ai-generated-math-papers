# commutative_algebra.py — Kommutative Algebra

## Überblick

Dieses Modul implementiert grundlegende Konzepte der **kommutativen Algebra**,
dem Fundament der algebraischen Geometrie und der algebraischen Zahlentheorie.

---

## Mathematischer Hintergrund

Die kommutative Algebra untersucht kommutative Ringe mit Eins und ihre Ideale.
Zentrale Sätze sind:

- **Hilbert-Basissatz**: Polynomringe über noetherschen Ringen sind noetherisch
- **Nakayama-Lemma**: Charakterisierung von Moduln über lokalen Ringen
- **Noether-Normalisierung**: Jede affine Algebra ist ganz über einem Polynomring

---

## Funktionen

### `localization(ring_elems, mult_set, mod)`

**Lokalisierung** $S^{-1}R$ eines Rings $R$ an einer multiplikativen Menge $S$:

Formale Brüche $\frac{a}{s}$ mit $a \in R$, $s \in S$, wobei
$\frac{a}{s} = \frac{a'}{s'}$ wenn $t(as' - a's) = 0$ für ein $t \in S$.

Für $R = \mathbb{Z}$, $S = \{p^n : n \geq 0\}$: Lokalisierung bei $p$
$$\mathbb{Z}_{(p)} = \left\{\frac{a}{b} \in \mathbb{Q} : p \nmid b\right\}$$

---

### `prime_spectrum(n)`

**Primspektrum** $\operatorname{Spec}(\mathbb{Z}/n\mathbb{Z})$:

Die Menge aller Primideale, versehen mit der Zariski-Topologie.

Für $R = \mathbb{Z}/n\mathbb{Z}$:
$$\operatorname{Spec}(\mathbb{Z}/n\mathbb{Z}) = \{(p) : p \text{ Primzahl}, p \mid n\}$$

Maximale Ideale $=$ Primideale (da $\mathbb{Z}/n\mathbb{Z}$ artinsch für $n > 0$).

**Beispiel:**
```python
prime_spectrum(12)
# → prime_ideals = [2, 3]
# → maximal_ideals = [2, 3]
```

---

### `nakayama_lemma_check(M_generators, I_generators, p)`

**Nakayama-Lemma:** Sei $R$ ein lokaler Ring mit maximalem Ideal $\mathfrak{m}$,
$M$ ein endlich erzeugter $R$-Modul. Dann gilt:
$$M = \mathfrak{m} \cdot M \implies M = 0$$

Äquivalent: Elemente $m_1, \ldots, m_n \in M$ erzeugen $M$, wenn ihre
Bilder in $M/\mathfrak{m}M$ den $R/\mathfrak{m}$-Vektorraum erzeugen.

---

### `integral_closure(d)`

**Ganzer Abschluss** von $\mathbb{Z}$ im quadratischen Zahlkörper $\mathbb{Q}(\sqrt{d})$:

Der **Ring der ganzen Zahlen** $\mathcal{O}_K$ hängt von $d \bmod 4$ ab:

$$\mathcal{O}_K = \begin{cases}
\mathbb{Z}\left[\frac{1+\sqrt{d}}{2}\right], & \text{Basis: } \left\{1, \frac{1+\sqrt{d}}{2}\right\}, \quad \Delta = d & \text{falls } d \equiv 1 \pmod{4} \\
\mathbb{Z}[\sqrt{d}], & \text{Basis: } \{1, \sqrt{d}\}, \quad \Delta = 4d & \text{sonst}
\end{cases}$$

**Bekannte Beispiele:**
- $d = -1$: $\mathcal{O}_K = \mathbb{Z}[i]$ (Gaußsche ganze Zahlen), $\Delta = -4$
- $d = -3$: $\mathcal{O}_K = \mathbb{Z}[\omega]$ (Eisensteinsche ganze Zahlen), $\Delta = -3$

```python
integral_closure(-1)
# → ring_basis = ['1', '√-1'], discriminant = -4, is_pid = True
integral_closure(-3)
# → discriminant = -3, is_pid = True
```

---

### `class_group_estimate(d)`

**Klassengruppe** $\operatorname{Cl}(\mathbb{Q}(\sqrt{d}))$ und Klassenzahl $h(d)$:

Die Klassengruppe misst die Abweichung von der eindeutigen Primfaktorzerlegung:
$$h(d) = 1 \iff \mathcal{O}_K \text{ ist Hauptidealring (PID)}$$

**Minkowski-Schranke** (obere Grenze für Normgen. prüfender Primideale):
$$M_K = \left(\frac{2}{\pi}\right)^{r_2} \cdot \frac{n!}{n^n} \cdot \sqrt{|\Delta|}$$

Für imaginär-quadratische Körper ($r_2 = 1$, $n = 2$):
$$M_K = \frac{2}{\pi} \sqrt{|\Delta|}$$

**Bekannte Klassenzahlen:**

| $d$ | $h(d)$ | Ring |
|-----|--------|------|
| $-1$ | 1 | $\mathbb{Z}[i]$ (PID) |
| $-3$ | 1 | $\mathbb{Z}[\omega]$ (PID) |
| $-5$ | 2 | $\mathbb{Z}[\sqrt{-5}]$ (kein PID) |
| $-23$ | 3 | — |
| $-163$ | 1 | (größtes imaginär-quad. PID) |

**Stark-Heegner-Theorem:** Die einzigen imaginär-quadratischen Körper mit $h=1$ sind
$d \in \{-1, -2, -3, -7, -11, -19, -43, -67, -163\}$.

```python
class_group_estimate(-1)   # → class_number = 1, is_pid = True
class_group_estimate(-5)   # → class_number = 2, is_pid = False
```

---

### `noether_normalization(poly_system, n_vars)`

**Noether-Normalisierungslemma:** Für jede endlich erzeugte $k$-Algebra
$A = k[x_1, \ldots, x_n]/I$ existieren algebraisch unabhängige Elemente
$y_1, \ldots, y_d \in A$, sodass $A$ über $k[y_1, \ldots, y_d]$ ganz ist.

Dabei ist $d = \operatorname{Krull-Dim}(A) = \operatorname{Transzendenzgrad}(\operatorname{Quot}(A)/k)$.

Geometrisch: $V(I) \subseteq \mathbb{A}^n$ hat eine endliche Abbildung auf $\mathbb{A}^d$.

**Beispiel:** Für eine Kurve $V(f) \subseteq \mathbb{A}^2$: $d = 1$.

---

### `hilbert_basis_theorem_verify(ideal_gens, p)`

**Hilbert-Basissatz:** Wenn $R$ noetherisch ist, dann ist $R[x]$ noetherisch.

**Folgerung:** $k[x_1, \ldots, x_n]$ ist noetherisch für jeden Körper $k$,
also ist jedes Ideal endlich erzeugt.

---

## Tests

Die Testsuite (`tests/test_commutative_algebra.py`) umfasst **32 Tests**:

| Testklasse | Tests | Abdeckung |
|-----------|-------|-----------|
| `TestLocalization` | 3 | Lokalisierung |
| `TestPrimeSpectrum` | 6 | Primspektrum |
| `TestNakayamaLemma` | 3 | Nakayama-Lemma |
| `TestIntegralClosure` | 8 | Ganzer Abschluss |
| `TestClassGroupEstimate` | 6 | Klassengruppe |
| `TestNoetherNormalization` | 3 | Noether-Normalisierung |
| `TestHilbertBasisTheorem` | 3 | Hilbert-Basissatz |

---

## Autor

**Kurt Ingwer** — @lastModified 2026-03-10
