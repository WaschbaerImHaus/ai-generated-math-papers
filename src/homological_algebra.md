# homological_algebra.py — Homologische Algebra

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Datei:** `src/homological_algebra.py`

---

## Übersicht

Dieses Modul implementiert die grundlegenden Strukturen der **homologischen Algebra**:

1. **Smith-Normalform** — Ganzzahlige Matrixzerlegung für Modulrechnung
2. **ChainComplex** — Kettenkomplexe mit ∂∘∂ = 0
3. **CochainComplex** — Duale Kokettenkomplexe
4. **Simpliziale Komplexe** — Konvertierung zu Kettenkomplexen
5. **Exakte Sequenzen** — Exaktheit, kurze exakte Sequenzen
6. **5-Lemma und Schlangen-Lemma** — Fundamentale Homomorphiesätze
7. **Freie Auflösungen** — Minimale projektive Auflösungen
8. **Ext und Tor** — Ableitete Funktoren für zyklische Moduln
9. **Universeller Koeffizienten-Satz** — UCT für Kohomologie

---

## 1. Das Fundamentale Axiom

Das zentrale Axiom der homologischen Algebra ist:

$$
\partial \circ \partial = 0
$$

Das bedeutet: Das Bild von $\partial_{n+1}$ liegt immer im Kern von $\partial_n$:

$$
\text{im}(\partial_{n+1}) \subseteq \ker(\partial_n)
$$

Dies ermöglicht die Definition der **Homologiegruppen** als Quotienten:

$$
H_n = \ker(\partial_n) / \text{im}(\partial_{n+1})
$$

---

## 2. Smith-Normalform

### 2.1 Definition

Jede ganzzahlige Matrix $M$ lässt sich schreiben als:

$$
M = P \cdot D \cdot Q
$$

wobei:
- $P, Q$ invertierbare ganzzahlige Matrizen sind ($|\det P| = |\det Q| = 1$)
- $D$ eine Diagonalmatrix mit $d_1 \mid d_2 \mid \cdots \mid d_r$ ist (**Teilbarkeitskette**)

Die $d_i$ heißen **Elementarteiler** von $M$.

### 2.2 Algorithmus

Der Algorithmus verwendet die euklidische Division zur schrittweisen Eliminierung:

1. Finde das kleinste Nicht-Null-Element als Pivotelement
2. Eliminiere alle anderen Einträge in der Pivot-Zeile und -Spalte via Divisionsalgorithmus
3. Prüfe Teilbarkeitsbedingung für die Smith-NF
4. Wiederhole für die Untermatrix

**Anwendung:** Die Homologiegruppen $H_n = \mathbb{Z}^\beta \oplus \bigoplus_i \mathbb{Z}/d_i\mathbb{Z}$ werden direkt aus den Elementarteilern berechnet.

---

## 3. Kettenkomplex

### 3.1 Definition

Ein **Kettenkomplex** $C_*$ ist eine Folge abelscher Gruppen mit Randoperatoren:

$$
\cdots \xrightarrow{\partial_{n+2}} C_{n+1} \xrightarrow{\partial_{n+1}} C_n \xrightarrow{\partial_n} C_{n-1} \xrightarrow{\partial_{n-1}} \cdots
$$

mit dem Axiom $\partial_n \circ \partial_{n+1} = 0$ für alle $n$.

### 3.2 Homologiegruppen

$$
H_n(C_*) = \ker(\partial_n) / \text{im}(\partial_{n+1})
$$

Für freie abelsche Gruppen gilt nach dem **Klassifikationssatz für endlich erzeugte abelsche Gruppen**:

$$
H_n \cong \mathbb{Z}^{\beta_n} \oplus \mathbb{Z}/d_1\mathbb{Z} \oplus \cdots \oplus \mathbb{Z}/d_k\mathbb{Z}
$$

- $\beta_n$: **Betti-Zahl** (freier Rang)
- $d_i$: **Torsionskoeffizienten** mit $d_1 \mid d_2 \mid \cdots \mid d_k$

### 3.3 Euler-Charakteristik

$$
\chi = \sum_n (-1)^n \cdot \text{rank}(C_n) = \sum_n (-1)^n \cdot \beta_n
$$

Beide Formeln liefern dasselbe Ergebnis — ein tiefes Theorem der homologischen Algebra.

**Beispiele:**

| Raum | $\chi$ |
|------|--------|
| Punkt $\{pt\}$ | 1 |
| Kreis $S^1$ | 0 |
| 2-Sphäre $S^2$ | 2 |
| Torus $T^2$ | 0 |
| Klein'sche Flasche | 0 |
| $\mathbb{RP}^2$ | 1 |

---

## 4. Kokettenkomplex

Der **duale** Komplex zu $C_*$ ist der Kokettenkomplex $C^*$:

$$
\cdots \xrightarrow{d^{n-2}} C^{n-1} \xrightarrow{d^{n-1}} C^n \xrightarrow{d^n} C^{n+1} \xrightarrow{d^{n+1}} \cdots
$$

Der Kobrandoperator ist die Transponierte: $d^n = (\partial_{n+1})^T$.

**Kohomologiegruppen:**

$$
H^n(C^*) = \ker(d^n) / \text{im}(d^{n-1})
$$

---

## 5. Simplizialer Komplex

### 5.1 Randoperator

Der **simpliziale Randoperator** $\partial_k$ wirkt auf $k$-Simplizes:

$$
\partial_k[v_0, v_1, \ldots, v_k] = \sum_{i=0}^{k} (-1)^i [v_0, \ldots, \hat{v}_i, \ldots, v_k]
$$

wobei $\hat{v}_i$ bedeutet, dass $v_i$ weggelassen wird.

**Beispiele:**
- $\partial_1([v_0, v_1]) = [v_1] - [v_0]$
- $\partial_2([v_0, v_1, v_2]) = [v_1, v_2] - [v_0, v_2] + [v_0, v_1]$

### 5.2 Randmatrix

Die Randmatrix $\partial_k$ hat:
- **Zeilen:** $(k-1)$-Simplizes
- **Spalten:** $k$-Simplizes
- **Einträge:** $\pm 1$ (Inzidenz mit Vorzeichen)

---

## 6. Exakte Sequenzen

### 6.1 Definition

Eine Folge abelscher Gruppen ist **exakt** bei $B$, wenn:

$$
\text{im}(f) = \ker(g) \quad \text{in} \quad A \xrightarrow{f} B \xrightarrow{g} C
$$

### 6.2 Kurze Exakte Sequenz

$$
0 \to A \xrightarrow{f} B \xrightarrow{g} C \to 0
$$

Dies bedeutet:
- $f$ ist **injektiv** (ker $f = 0$)
- $g$ ist **surjektiv** (im $g = C$)
- $\text{im}(f) = \ker(g)$

**Aufspaltungssatz:** Für freie abelsche Gruppen zerfällt jede kurze exakte Sequenz: $B \cong A \oplus C$.

---

## 7. Das Schlangen-Lemma

Das **Schlangen-Lemma** (Snake Lemma) ist eines der wichtigsten Werkzeuge der homologischen Algebra.

Gegeben ein kommutatives Diagramm mit exakten Zeilen:

$$
\begin{array}{ccccccc}
0 & \to & A & \xrightarrow{f} & B & \xrightarrow{g} & C & \to 0 \\
& & \downarrow\alpha & & \downarrow\beta & & \downarrow\gamma & \\
0 & \to & A' & \xrightarrow{f'} & B' & \xrightarrow{g'} & C' & \to 0
\end{array}
$$

Dann gibt es eine lange exakte Sequenz:

$$
0 \to \ker\alpha \to \ker\beta \to \ker\gamma \xrightarrow{\delta} \text{coker}\,\alpha \to \text{coker}\,\beta \to \text{coker}\,\gamma \to 0
$$

Der **Verbindungsmorphismus** $\delta$ wird durch Diagrammjagd konstruiert.

---

## 8. Das 5-Lemma

In einem kommutativen Diagramm mit exakten Zeilen:

$$
\begin{array}{ccccccccc}
A_1 & \to & A_2 & \to & A_3 & \to & A_4 & \to & A_5 \\
\downarrow f_1 & & \downarrow f_2 & & \downarrow f_3 & & \downarrow f_4 & & \downarrow f_5 \\
B_1 & \to & B_2 & \to & B_3 & \to & B_4 & \to & B_5
\end{array}
$$

**5-Lemma:** Wenn $f_1, f_2, f_4, f_5$ Isomorphismen sind, dann ist auch $f_3$ ein Isomorphismus.

---

## 9. Freie Auflösungen

Eine **freie Auflösung** eines Moduls $M$ ist ein exakter Komplex aus freien Moduln:

$$
\cdots \to F_2 \to F_1 \to F_0 \to M \to 0
$$

**Minimale Auflösung von $\mathbb{Z}/n\mathbb{Z}$:**

$$
0 \to \mathbb{Z} \xrightarrow{\cdot n} \mathbb{Z} \to \mathbb{Z}/n\mathbb{Z} \to 0
$$

**Projektive Dimension:**
- $\text{pd}(\mathbb{Z}) = 0$ (bereits frei/projektiv)
- $\text{pd}(\mathbb{Z}/n\mathbb{Z}) = 1$ (Auflösung der Länge 1)
- $\text{pd}(0) = -\infty$ (Konvention $-1$)

---

## 10. Ext-Gruppen

Die **Ext-Gruppen** klassifizieren Erweiterungen:

Aus der freien Auflösung $0 \to \mathbb{Z} \xrightarrow{n} \mathbb{Z} \to \mathbb{Z}/n\mathbb{Z} \to 0$ und Anwendung von $\text{Hom}(-, \mathbb{Z}/m\mathbb{Z})$ ergibt sich:

$$
\text{Ext}^k(\mathbb{Z}/n\mathbb{Z},\, \mathbb{Z}/m\mathbb{Z}) \cong \begin{cases}
\mathbb{Z}/\gcd(n,m)\mathbb{Z} & k = 0 \text{ (Hom)} \\
\mathbb{Z}/\gcd(n,m)\mathbb{Z} & k = 1 \\
0 & k \geq 2
\end{cases}
$$

---

## 11. Tor-Gruppen

Die **Tor-Gruppen** messen Torsion im Tensorprodukt:

$$
\text{Tor}_k(\mathbb{Z}/n\mathbb{Z},\, \mathbb{Z}/m\mathbb{Z}) \cong \begin{cases}
\mathbb{Z}/\gcd(n,m)\mathbb{Z} & k = 0 \text{ (Tensorprodukt)} \\
\mathbb{Z}/\gcd(n,m)\mathbb{Z} & k = 1 \\
0 & k \geq 2
\end{cases}
$$

Aus der freien Auflösung und Tensorierung mit $\mathbb{Z}/m\mathbb{Z}$:

$$
\mathbb{Z}/m\mathbb{Z} \xrightarrow{\cdot n} \mathbb{Z}/m\mathbb{Z}
$$

Der Kern ist $\text{Tor}_1 \cong \ker(\cdot n: \mathbb{Z}/m\mathbb{Z} \to \mathbb{Z}/m\mathbb{Z}) \cong \mathbb{Z}/\gcd(n,m)\mathbb{Z}$.

---

## 12. Universeller Koeffizienten-Satz (UCT)

Der **Universelle Koeffizienten-Satz** verknüpft Homologie und Kohomologie:

$$
H^n(X;\,\mathbb{Z}) \cong \text{Hom}(H_n(X),\,\mathbb{Z}) \oplus \text{Ext}^1(H_{n-1}(X),\,\mathbb{Z})
$$

**Interpretation:**
- Freier Teil: $\text{Hom}(\mathbb{Z}^\beta \oplus \text{Tor},\, \mathbb{Z}) \cong \mathbb{Z}^\beta$
- Torsion: $\text{Ext}^1(\mathbb{Z}/d\mathbb{Z},\, \mathbb{Z}) \cong \mathbb{Z}/d\mathbb{Z}$

**Merkregel:** Die Kohomologie "erbt" die Torsion um einen Grad höher als die Homologie.

**Beispiel für $S^2$ (Sphäre):**

| Grad $n$ | $H_n(S^2)$ | $H^n(S^2)$ |
|----------|-----------|-----------|
| 0 | $\mathbb{Z}$ | $\mathbb{Z}$ |
| 1 | $0$ | $0$ |
| 2 | $\mathbb{Z}$ | $\mathbb{Z}$ |

**Beispiel für $\mathbb{RP}^2$:**

| Grad $n$ | $H_n(\mathbb{RP}^2)$ | $H^n(\mathbb{RP}^2)$ |
|----------|-----------|-----------|
| 0 | $\mathbb{Z}$ | $\mathbb{Z}$ |
| 1 | $\mathbb{Z}/2\mathbb{Z}$ | $0$ |
| 2 | $0$ | $\mathbb{Z}/2\mathbb{Z}$ |

---

## Verwendete Bibliotheken

| Bibliothek | Verwendung |
|------------|-----------|
| `numpy` | Matrixoperationen, Randmatrizen |
| `math` | Größter gemeinsamer Teiler, Faktorielle |

---

## Tests

Alle Tests in `tests/test_homological_algebra.py`. Abdeckung: 60+ Tests über alle Klassen und Funktionen.
