# ramsey_numbers.py – Dokumentation

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `src/ramsey_numbers.py`

---

## Überblick

Das Modul implementiert die Theorie der **Ramsey-Zahlen** $R(s,t)$ – eine der zentralen Fragestellungen der extremalen Graphentheorie und kombinatorischen Mathematik.

---

## Mathematischer Hintergrund

### Definition: Ramsey-Zahl

Die **Ramsey-Zahl** $R(s,t)$ ist die kleinste natürliche Zahl $N$, sodass jeder Graph $G$ auf $N$ Knoten entweder

- eine **Clique der Größe $s$** (vollständiger Teilgraph $K_s$), oder
- eine **unabhängige Menge der Größe $t$** (keine zwei Knoten verbunden)

enthält.

**Äquivalente Formulierung (Kantenfarbenversion):**
Jede 2-Färbung der Kanten von $K_N$ enthält entweder eine rote $K_s$ oder eine blaue $K_t$.

### Fundamentalsatz (Ramsey 1930)

> **Satz (Ramsey 1930):** Für alle $s, t \geq 1$ existiert $R(s,t) < \infty$.

---

## Bekannte exakte Werte

| $R(s,t)$ | Wert | Quelle |
|----------|------|--------|
| $R(3,3)$ | 6 | Greenwood-Gleason 1955 |
| $R(3,4)$ | 9 | Greenwood-Gleason 1955 |
| $R(3,5)$ | 14 | Greenwood-Gleason 1955 |
| $R(3,6)$ | 18 | Kéry 1964 |
| $R(3,7)$ | 23 | Kalbfleisch 1966 |
| $R(3,8)$ | 28 | McKay-Min 1992 |
| $R(3,9)$ | 36 | Grinstead-Roberts 1982 |
| $R(4,4)$ | 18 | Greenwood-Gleason 1955 |
| $R(4,5)$ | 25 | McKay-Radziszowski 1995 |

**Aktuelle Schranken (2025):**

$$43 \leq R(5,5) \leq 48$$

---

## Schranken für R(k,k)

### Obere Schranke: Binomialkoeffizient (Erdős-Szekeres 1935)

$$R(s,t) \leq \binom{s+t-2}{s-1}$$

**Beweis:** Induktion über $s+t$: $R(s,t) \leq R(s-1,t) + R(s,t-1)$.

Beispiel: $R(5,5) \leq \binom{8}{4} = 70$.

### Untere Schranke: Probabilistische Methode (Erdős 1947)

$$R(k,k) > 2^{k/2}$$

**Beweis (probabilistisch):**
Sei $G$ ein zufälliger Graph auf $n$ Knoten (jede Kante unabhängig mit $p = 1/2$). Für eine $k$-Teilmenge $S$:

$$\Pr[S \text{ ist Clique oder I-Menge}] = 2 \cdot 2^{-\binom{k}{2}}$$

Für $n = 2^{k/2}$: Erwartete Anzahl monochromatischer $K_k$ $< 1$ → Zeuge existiert.

### Verbesserte obere Schranke (Sah 2023)

$$R(k,k) \leq (4 - \varepsilon)^k \quad \text{für ein explizites } \varepsilon > 0$$

**Bedeutung:** Erster Durchbruch seit 1935, der die $4^k$-Schranke übertrifft.

### Spencer-Untergrenze (1975)

$$R(k,k) > c \cdot k \cdot \frac{2^{k/2}}{\sqrt{\ln k}}$$

Nutzt das **Lovász Lokale Lemma** für eine explizite Konstruktion.

---

## Klasse: `RamseyNumbers`

### Methoden

#### `get(s, t) → Optional[int]`

Gibt bekannten exakten Wert $R(s,t)$ zurück oder `None`.

#### `get_bounds(s, t) → Tuple[int, int]`

Gibt Schranken `(lower, upper)` zurück. Exakte Werte → `(R, R)`.

#### `has_clique(graph, size) → bool`

Prüft Clique-Existenz via Bron-Kerbosch-Algorithmus (NetworkX).

#### `has_independent_set(graph, size) → bool`

Äquivalent zu `has_clique(complement(graph), size)`.

#### `is_ramsey_witness(graph, s, t) → bool`

Prüft ob $G$ ein Untergrenzen-Zeugnis für $R(s,t) > |V(G)|$ ist.

#### `build_paley_graph(q) → nx.Graph`

Konstruiert Paley-Graphen $P(q)$ für Primzahlen $q \equiv 1 \pmod{4}$:

$$\{i,j\} \in E(P(q)) \iff (i-j) \text{ ist quadratischer Rest mod } q$$

| $P(q)$ | $\omega$ | $\alpha$ | Zeuge für |
|--------|----------|----------|-----------|
| $P(5)$ | 2 | 2 | $R(3,3) > 5$ |
| $P(13)$ | 3 | 3 | $R(4,4) > 13$ |
| $P(17)$ | 3 | 3 | selbst-komplementär |
| $P(29)$ | 4 | 4 | $R(5,5) > 29$ |

#### `build_r55_lower_bound_witness() → nx.Graph`

Gibt $P(29)$ zurück: Zeugnis für $R(5,5) > 29$ (omega = alpha = 4).

#### `verify_r33_equals_6() → Dict`

Vollständige Verifikation: $C_5$ ist Zeuge $R(3,3) > 5$, $K_6$ erzwingt $K_3$.

#### `verify_via_sat(s, t, n) → Optional[bool]`

SAT-basierte Verifikation über PySAT (MiniSat22):

- **SAT** → Zeugengraph gefunden → $R(s,t) > n$
- **UNSAT** → kein Zeuge → $R(s,t) \leq n$

**Klauselkodierung:**
- Variable $x_{i,j}$ = 1 gdw. Kante $\{i,j\}$ vorhanden.
- Keine $s$-Clique: Für jede $s$-Teilmenge mind. eine fehlende Kante.
- Keine $t$-Unabhängige-Menge: Für jede $t$-Teilmenge mind. eine Kante.

---

## Klasse: `RamseyBounds`

### Methoden

| Methode | Formel | Typ |
|---------|--------|-----|
| `upper_binomial(s,t)` | $\binom{s+t-2}{s-1}$ | Obere Schranke |
| `lower_erdos(k)` | $\lfloor 2^{k/2} \rfloor$ | Untere Schranke |
| `upper_sah_2023(k)` | $(4-\varepsilon)^k$ | Obere Schranke |
| `upper_spencer_1975(k)` | $c \cdot k \cdot 2^{k/2} / \sqrt{\ln k}$ | Untere Schranke |
| `recursive_upper(s,t)` | $R(s-1,t)+R(s,t-1)$ | Obere Schranke |
| `diagonal_bounds(k_max)` | Tabelle für $k=3\ldots k_{\max}$ | Übersicht |
| `multiplicity(s,t,n)` | Goodman-Formel | Multiplizität |

### Goodman-Formel für $M(3,3,n)$

Die Anzahl monochromatischer Dreiecke in jeder 2-Färbung von $K_n$:

$$M(3,3,n) \geq \begin{cases} \frac{n(n-1)(n-5)}{24} & n \text{ gerade} \\ \frac{n(n-2)(n-4)}{24} & n \text{ ungerade} \end{cases}$$

---

## Beispiel-Verwendung

```python
from ramsey_numbers import RamseyNumbers, RamseyBounds

r = RamseyNumbers()
b = RamseyBounds()

# Exakter Wert
print(r.get(3, 3))   # 6
print(r.get(4, 4))   # 18

# Schranken
lower, upper = r.get_bounds(5, 5)
print(lower, upper)  # 43, 48

# Binomial-Schranke
print(b.upper_binomial(5, 5))   # 70
print(b.upper_binomial(3, 4))   # 10

# Paley-Zeugnis
p29 = r.build_r55_lower_bound_witness()
print(r.is_ramsey_witness(p29, 5, 5))  # True (R(5,5) > 29)

# SAT-Verifikation
print(r.verify_via_sat(3, 3, 5))  # True  (Zeuge existiert, R(3,3) > 5)
print(r.verify_via_sat(3, 3, 6))  # False (kein Zeuge, R(3,3) ≤ 6)
```

---

## Offene Probleme

1. **$R(5,5)$**: Genauer Wert unbekannt ($43 \leq R(5,5) \leq 48$). Seit 1995 keine Verbesserung der unteren Schranke.
2. **$R(6,6)$**: $102 \leq R(6,6) \leq 165$.
3. **Asymptotisch**: Ist $\lim_{k \to \infty} R(k,k)^{1/k}$ wohldefiniert? Vermutlich $\sqrt{2} \leq c \leq 4$.

---

## Literatur

- Ramsey, F.P. (1930): On a Problem of Formal Logic. *Proc. London Math. Soc.* 30, 264–286.
- Erdős, P. (1947): Some Remarks on the Theory of Graphs. *Bull. AMS* 53, 292–294.
- Erdős, P. und Szekeres, G. (1935): A Combinatorial Problem in Geometry. *Compositio Math.* 2, 463–470.
- Radziszowski, S.P. (2021): Small Ramsey Numbers. *Electronic J. Combinatorics DS6* (dynamisch aktualisiert).
- Sah, A. (2023): Diagonal Ramsey via effective quasirandomness. *Duke Math. J.* 172(3).
- Spencer, J. (1975): Ramsey's Theorem – A New Lower Bound. *J. Combinatorial Theory (A)* 18, 108–115.
