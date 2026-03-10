# recursion_theory.py — Rekursionstheorie / Berechenbarkeitstheorie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Datei:** `src/recursion_theory.py`

---

## Übersicht

Dieses Modul implementiert die wichtigsten Konzepte der **Rekursionstheorie** (auch: Berechenbarkeitstheorie), einem Kerngebiet der theoretischen Informatik und mathematischen Logik. Es behandelt die Grenzen algorithmischer Berechenbarkeit, die Klassifikation von Problemen nach ihrer Lösbarkeit und die Komplexität von Berechnungen.

---

## 1. Turing-Maschinen

### Formale Definition

Eine **deterministische Turing-Maschine** (DTM) ist ein 7-Tupel:

$$M = (Q, \Sigma, \Gamma, \delta, q_0, q_{acc}, q_{rej})$$

| Symbol | Bedeutung |
|--------|-----------|
| $Q$ | Endliche Zustandsmenge |
| $\Sigma$ | Eingabealphabet ($\Sigma \subseteq \Gamma \setminus \{\sqcup\}$) |
| $\Gamma$ | Bandalphabet (enthält Leerzeichen $\sqcup$) |
| $\delta$ | Übergangsfunktion: $Q \times \Gamma \to Q \times \Gamma \times \{L, R\}$ |
| $q_0$ | Startzustand |
| $q_{acc}$ | Akzeptierzustand |
| $q_{rej}$ | Verwerfzustand |

### Konfiguration und Berechnung

Eine **Konfiguration** der TM ist ein Tripel $(w_1, q, w_2)$, wobei:
- $w_1 \in \Gamma^*$ das Band links vom Kopf (inklusive aktuelle Zelle)
- $q \in Q$ der aktuelle Zustand
- $w_2 \in \Gamma^*$ der Rest des Bandes

Die **Berechnungsrelation** $\vdash$ definiert:

$$w_1 a q b w_2 \vdash w_1 q' a' b w_2 \quad \text{falls } \delta(q, b) = (q', b', L), a' = b'$$

### Implementierte Klassen

#### `TuringMachine`

```python
tm = TuringMachine(
    states={'q0', 'accept', 'reject'},
    input_alphabet={'a', 'b'},
    tape_alphabet={'a', 'b', '_'},
    transitions={(state, symbol): (new_state, write, direction)},
    initial_state='q0',
    accept_state='accept',
    reject_state='reject'
)
result = tm.run('abba')  # {'accepted': True, 'steps': 12, 'tape': '...', 'halted': True}
```

#### Methoden

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `run()` | `(input_str, max_steps=10000) → dict` | Vollständige Simulation |
| `accepts()` | `(input_str, max_steps=10000) → bool\|None` | Akzeptiert/verwirft/Timeout |
| `configuration_sequence()` | `(input_str, max_steps=20) → list[str]` | Berechnungsverlauf |

### Beispiel-TMs

#### Palindrome über $\{a, b\}$

```
Idee: Lösche äußerste Zeichen paarweise, prüfe Gleichheit.
Zeitkomplexität: O(n²) Schritte
```

#### Binärinkrement

$$\text{bin}(n) \to \text{bin}(n+1)$$

```
Algorithmus: Addiere 1 von rechts mit Übertragspropagation.
Zeitkomplexität: O(n) im Durchschnitt, O(n) im Worst Case
```

#### Unäre Addition

$$1^n \cdot 0 \cdot 1^m \to 1^{n+m}$$

```
Ersetze Trenner 0 durch 1, lösche letzte 1.
Zeitkomplexität: O(n+m)
```

---

## 2. μ-rekursive Funktionen

### Grundfunktionen

Die drei **Grundfunktionen** der μ-Rekursion:

$$Z(n) = 0 \quad \text{(Nullfunktion)}$$

$$S(n) = n + 1 \quad \text{(Nachfolgerfunktion)}$$

$$U_i^n(x_1, \ldots, x_n) = x_i \quad \text{(Projektionsfunktion)}$$

### Primitive Rekursion

Gegeben $g: \mathbb{N}^k \to \mathbb{N}$ und $h: \mathbb{N}^{k+2} \to \mathbb{N}$, definiert **primitive Rekursion**:

$$f(\bar{x}, 0) = g(\bar{x})$$
$$f(\bar{x}, n+1) = h(\bar{x}, n, f(\bar{x}, n))$$

**Beispiel — Addition:**
$$\text{add}(x, 0) = x, \quad \text{add}(x, n+1) = S(\text{add}(x, n))$$

**Beispiel — Multiplikation:**
$$\text{mul}(x, 0) = 0, \quad \text{mul}(x, n+1) = \text{add}(\text{mul}(x, n), x)$$

### μ-Operator (Minimierungsoperator)

Der **μ-Operator** erweitert primitive Rekursion zu partieller μ-Rekursion:

$$\mu y[P(\bar{x}, y)] = \min\{y \in \mathbb{N} : P(\bar{x}, y) = 0\}$$

Falls kein solches $y$ existiert, **divergiert** die Funktion (Partiellheit).

**Churchsche These:** Jede intuitiv berechenbare Funktion ist μ-rekursiv.

### Ackermann-Funktion

$$A(0, n) = n + 1$$
$$A(m+1, 0) = A(m, 1)$$
$$A(m+1, n+1) = A(m, A(m+1, n))$$

Wachstumsverhalten:

| $m$ | $A(m, n)$ | Wachstumsklasse |
|-----|-----------|-----------------|
| 0 | $n+1$ | Linear |
| 1 | $n+2$ | Linear |
| 2 | $2n+3$ | Linear |
| 3 | $2^{n+3} - 3$ | Exponentiell |
| 4 | $\underbrace{2^{2^{\cdot^{\cdot^{\cdot^2}}}}}_{n+3} - 3$ | Tetrational |

**Satz:** $A$ ist total rekursiv, aber **nicht** primitiv rekursiv.
*Beweis:* Für jede pr.-rek. Funktion $f$ gibt es $m_0$ mit $A(m_0, n) > f(n)$ für fast alle $n$.

---

## 3. Halteproblem und Entscheidbarkeit

### Halteproblem (HALT)

$$\text{HALT} = \{\langle M, w \rangle : M \text{ hält auf Eingabe } w\}$$

**Satz (Turing 1936):** HALT ist **unentscheidbar**.

**Beweis** (Cantor-Diagonalargument):

1. **Annahme:** Es gibt TM $H$ die HALT entscheidet.
2. **Konstruktion:** Definiere $D(\langle M \rangle)$:
   - Falls $H(\langle M, M \rangle) = \text{accept}$: Schleife (halte nicht)
   - Falls $H(\langle M, M \rangle) = \text{reject}$: Halte (akzeptiere)
3. **Widerspruch:** Was passiert bei $D(\langle D \rangle)$?
   - $D$ hält auf $\langle D \rangle$ $\Rightarrow$ $H$ akzeptiert $\Rightarrow$ $D$ schleift $\Rightarrow$ Widerspruch
   - $D$ schleift auf $\langle D \rangle$ $\Rightarrow$ $H$ verwirft $\Rightarrow$ $D$ hält $\Rightarrow$ Widerspruch
4. **Schluss:** $H$ kann nicht existieren. $\square$

### Satz von Rice

**Satz (Rice 1951):** Sei $\mathcal{C}$ eine nicht-triviale Klasse partiell berechenbarer Funktionen. Dann ist

$$L_\mathcal{C} = \{\langle M \rangle : \varphi_M \in \mathcal{C}\}$$

unentscheidbar.

*Nicht-trivial:* $\emptyset \subsetneq \mathcal{C} \subsetneq \{\text{alle partiell ber. Funktionen}\}$

---

## 4. Arithmetische Hierarchie

### Definition

Die **arithmetische Hierarchie** klassifiziert Teilmengen von $\mathbb{N}$ nach der logischen Komplexität ihrer Definitionen:

$$\Sigma_0 = \Pi_0 = \Delta_1: \text{Entscheidbare Mengen}$$
$$\Sigma_1: \exists y \, R(x, y) \text{ mit } R \in \Delta_1$$
$$\Pi_1: \forall y \, R(x, y) \text{ mit } R \in \Delta_1$$
$$\Sigma_{n+1}: \exists y \, \varphi(x, y) \text{ mit } \varphi \in \Pi_n$$
$$\Pi_{n+1}: \forall y \, \varphi(x, y) \text{ mit } \varphi \in \Sigma_n$$
$$\Delta_{n+1} = \Sigma_{n+1} \cap \Pi_{n+1}$$

### Struktur

```
         Σ₃
        /   \
      Δ₃
      / \
    Σ₂   Π₂
     \   /
      Δ₂
      / \
    Σ₁   Π₁
     \   /
      Δ₁ = Σ₀ = Π₀ (Entscheidbar)
```

### Vollständige Probleme

| Klasse | Vollständiges Problem |
|--------|----------------------|
| $\Sigma_1$ | ATM = $\{\langle M, w \rangle : M \text{ akzeptiert } w\}$ |
| $\Pi_1$ | co-HALT |
| $\Sigma_2$ | FIN = $\{\langle M \rangle : |L(M)| < \infty\}$ |
| $\Pi_2$ | TOTAL = $\{\langle M \rangle : M \text{ hält auf allen Eingaben}\}$ |
| $\Sigma_3$ | REC = $\{\langle M \rangle : L(M) \text{ entscheidbar}\}$ |

---

## 5. Post'sches Korrespondenzproblem (PCP)

**Problem:** Gegeben Dominos $\{(t_1, b_1), \ldots, (t_n, b_n)\}$ mit $t_i, b_i \in \Sigma^+$.
Existiert eine nicht-leere Sequenz $i_1, \ldots, i_k$ mit:

$$t_{i_1} t_{i_2} \cdots t_{i_k} = b_{i_1} b_{i_2} \cdots b_{i_k}$$

**Satz:** PCP ist **unentscheidbar** (Post 1946).

Die Implementierung verwendet BFS/Backtracking für kleine Instanzen (Sequenzlänge $\leq 8$).

---

## 6. Komplexitätstheorie

### Zeitkomplexitätsklassen

$$\mathbf{P} = \bigcup_{k \geq 1} \text{DTIME}(n^k)$$

$$\mathbf{NP} = \bigcup_{k \geq 1} \text{NTIME}(n^k) = \{L : \exists \text{ polynom. Verifizierer für } L\}$$

$$\mathbf{PSPACE} = \bigcup_{k \geq 1} \text{DSPACE}(n^k)$$

$$\mathbf{EXPTIME} = \bigcup_{k \geq 1} \text{DTIME}(2^{n^k})$$

**Bekannte Inklusionen:**
$$\mathbf{P} \subseteq \mathbf{NP} \cap \mathbf{co\text{-}NP} \subseteq \mathbf{NP} \subseteq \mathbf{PSPACE} \subseteq \mathbf{EXPTIME}$$

**Offen:** $\mathbf{P} \overset{?}{=} \mathbf{NP}$ (Millennium-Problem, Clay Institute, $1.000.000\$$)

### Cook-Levin-Satz

**Satz:** SAT ist **NP-vollständig** (Cook 1971, Levin 1973).

Beweis-Idee (Tableau-Konstruktion):
- Für NP-TM $M$ und Eingabe $w$ ($|w| = n$): Erstelle Tableau $T[0..p(n)][0..p(n)]$
- Variablen $C[t, j, a]$: "Zur Zeit $t$ steht an Position $j$ das Symbol $a$"
- Klauseln kodieren: Eindeutigkeit, Startzustand, Übergänge, Akzeptanz
- Formelgröße: $O(p(n)^2)$ — polynomiell

### Polynomielle Hierarchie

$$\Sigma^p_0 = \Pi^p_0 = \mathbf{P}$$
$$\Sigma^p_1 = \mathbf{NP}, \quad \Pi^p_1 = \mathbf{co\text{-}NP}$$
$$\Sigma^p_{n+1} = \mathbf{NP}^{\Sigma^p_n}$$

$$\mathbf{PH} = \bigcup_{n \geq 0} \Sigma^p_n \subseteq \mathbf{PSPACE}$$

---

## 7. Kolmogorov-Komplexität

### Definition

Die **Kolmogorov-Komplexität** eines Strings $s$ bezüglich universeller TM $U$:

$$K(s) = \min\{|p| : U(p) = s\}$$

### Eigenschaften

- $K(s) \leq |s| + c$ (triviale obere Schranke)
- $K(s)$ ist **nicht berechenbar** (Reduktion auf Halteproblem)
- Fast alle Strings sind **inkompressibel**: $|\{s \in \{0,1\}^n : K(s) < n - c\}| \leq 2^{n-c}$

### Zählargument (Inkompressibilität)

Es gibt $2^n$ Strings der Länge $n$, aber nur $\sum_{k=0}^{n-1} 2^k = 2^n - 1 < 2^n$ kürzere Programme.
Also hat mindestens ein String der Länge $n$ keine kürzere Beschreibung. $\square$

### Implementierung

Die Funktion `kolmogorov_complexity_approx()` gibt eine **obere Schranke** via `zlib` (DEFLATE) zurück:

```python
result = kolmogorov_complexity_approx('a' * 1000)
# {'compression_ratio': 0.012, 'complexity_class': 'sehr niedrig (hochgradig strukturiert)', ...}
```

---

## Schnellreferenz

| Funktion | Beschreibung |
|----------|--------------|
| `TuringMachine` | DTM-Simulator |
| `tm_recognizes_palindromes()` | Beispiel-TM für Palindrome |
| `tm_binary_increment()` | Binärzahl um 1 erhöhen |
| `tm_unary_addition()` | Unäre Addition |
| `zero_function(n)` | $Z(n) = 0$ |
| `successor_function(n)` | $S(n) = n+1$ |
| `projection(args, i)` | $U_i^n$ |
| `primitive_recursion(g, h)` | PR-Schema |
| `mu_operator(P)` | $\mu y[P(\cdot, y) = 0]$ |
| `ackermann_function(m, n)` | Ackermann-Funktion |
| `halting_problem_undecidability_proof()` | Beweis HALT unentscheidbar |
| `rice_theorem_demo()` | Satz von Rice |
| `arithmetical_hierarchy()` | $\Sigma_n/\Pi_n/\Delta_n$ |
| `post_correspondence_problem(dominoes)` | PCP (Backtracking) |
| `np_complete_problems()` | Liste NP-vollständiger Probleme |
| `kolmogorov_complexity_approx(s)` | Obere Schranke $K(s)$ via zlib |

---

## Literatur

- Turing, A.M. (1936): *On Computable Numbers, with an Application to the Entscheidungsproblem*
- Rice, H.G. (1951): *Classes of Recursively Enumerable Sets and Their Decision Problems*
- Cook, S. (1971): *The Complexity of Theorem-Proving Procedures*
- Rogers, H. (1987): *Theory of Recursive Functions and Effective Computability*
- Sipser, M. (2013): *Introduction to the Theory of Computation* (3. Aufl.)
