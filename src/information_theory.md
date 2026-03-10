# Modul: Informationstheorie (`information_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `information_theory.py` implementiert die fundamentalen Konzepte der Informationstheorie nach Claude Shannon (1948). Es umfasst Shannon-Entropie und verwandte Maße, Kullback-Leibler-Divergenz und Informationsabstandsmaße, Kanalkapazität (inkl. Blahut-Arimoto-Algorithmus), Quellencodierung (Huffman, Arithmetik, Lempel-Ziv), Fehlerkorrektur-Grundlagen sowie differentielle Entropie für kontinuierliche Zufallsvariablen.

---

## Mathematischer Hintergrund

### Shannon-Entropie

Die Shannon-Entropie quantifiziert den mittleren Informationsgehalt einer diskreten Zufallsvariablen $X$ mit Wahrscheinlichkeiten $\{p_i\}$:

$$H(X) = -\sum_{i} p_i \log_2 p_i \quad [\text{Bit}]$$

$$H(X) = -\sum_{i} p_i \ln p_i \quad [\text{Nat}]$$

**Maximale Entropie** bei $n$ gleichwahrscheinlichen Ereignissen: $H_\mathrm{max} = \log_2 n$

**Binäre Entropie:** $h(p) = -p\log_2 p - (1-p)\log_2(1-p)$

**Verbundentropie:** $H(X,Y) = -\sum_{i,j} p(x_i, y_j)\log_2 p(x_i, y_j)$

**Bedingte Entropie:** $H(Y|X) = H(X,Y) - H(X)$

**Gegenseitige Information:**
$$I(X;Y) = H(X) + H(Y) - H(X,Y) = H(X) - H(X|Y)$$

### Kullback-Leibler-Divergenz

Die KL-Divergenz (relative Entropie) misst den Informationsverlust beim Approximieren von $P$ durch $Q$:

$$D_\mathrm{KL}(P \| Q) = \sum_i p_i \log_2 \frac{p_i}{q_i} \geq 0$$

(Gleichheit genau dann, wenn $P = Q$; nicht symmetrisch!)

**Jensen-Shannon-Divergenz** (symmetrisch, beschränkt auf $[0,1]$):
$$\mathrm{JS}(P \| Q) = \frac{1}{2} D_\mathrm{KL}(P \| M) + \frac{1}{2} D_\mathrm{KL}(Q \| M), \quad M = \frac{P+Q}{2}$$

**Kreuzentropie:** $H(P, Q) = H(P) + D_\mathrm{KL}(P \| Q) = -\sum_i p_i \log_2 q_i$

**Rényi-Entropie der Ordnung $\alpha$:**
$$H_\alpha(X) = \frac{1}{1-\alpha} \log_2 \sum_i p_i^\alpha$$

**Tsallis-Entropie der Ordnung $q$:**
$$S_q(X) = \frac{1}{q-1}\left(1 - \sum_i p_i^q\right)$$

**Total-Variationsabstand:**
$$d_\mathrm{TV}(P, Q) = \frac{1}{2}\sum_i |p_i - q_i|$$

### Kanalkapazität

Die Kapazität $C$ eines Informationskanals ist das Maximum der gegenseitigen Information über alle Eingabeverteilungen:

$$C = \max_{p(x)} I(X;Y) \quad [\text{Bit/Kanalnutzung}]$$

**Binärer symmetrischer Kanal (BSC)** mit Fehlerwahrscheinlichkeit $p$:
$$C_\mathrm{BSC} = 1 - h(p) = 1 + p\log_2 p + (1-p)\log_2(1-p)$$

**Binärer Löschkanal (BEC)** mit Löschwahrscheinlichkeit $\varepsilon$:
$$C_\mathrm{BEC} = 1 - \varepsilon$$

**AWGN-Kanal** (Shannon-Hartley-Theorem) mit Signal-Rausch-Verhältnis $\text{SNR}$:
$$C_\mathrm{AWGN} = \log_2(1 + \text{SNR}) \quad [\text{Bit/s/Hz}]$$

**Blahut-Arimoto-Algorithmus:** Iterativer Algorithmus zur Berechnung von $C$ für beliebige diskrete Kanäle; konvergiert garantiert gegen das globale Maximum.

### Quellencodierungssatz (Shannon, 1948)

Für eine Quelle mit Entropie $H(X)$ und einem Huffman-Code mit mittlerer Codelänge $\bar{L}$ gilt:
$$H(X) \leq \bar{L} < H(X) + 1 \quad [\text{Bit/Symbol}]$$

**Huffman-Code:** Optimaler präfixfreier Code – konstruiert als Binärbaum durch schrittweises Zusammenfassen der beiden unwahrscheinlichsten Symbole.

**Arithmetische Codierung:** Codiert eine Folge als einzelne reelle Zahl im Intervall $[0,1)$ – erreicht $H(X)$ asymptotisch ohne die $+1$-Schranke des Huffman-Codes.

**Lempel-Ziv-Komplexität:** Misst die Beschreibungskomplexität einer Folge durch Zählen neuer Teilworte.

**Kolmogorov-Komplexität** $K(x)$: Länge des kürzesten Programms, das $x$ ausgibt. Nicht berechenbar, aber approximierbar durch Komprimierungsalgorithmen.

### Fehlerkorrektur (Grundlagen)

**Hamming-Distanz:** $d_H(c_1, c_2) = |\{i : c_1[i] \neq c_2[i]\}|$

**Hamming-Gewicht:** $w_H(c) = d_H(c, \mathbf{0})$

**Singleton-Schranke:** $d \leq n - k + 1$ (Codes, die sie erreichen: MDS-Codes)

**Hamming-Schranke (Kugelpackung):** $\sum_{i=0}^t \binom{n}{i} \leq 2^{n-k}$

### Differentielle Entropie

Für kontinuierliche Zufallsvariablen mit Dichte $f(x)$:

$$h(X) = -\int_{-\infty}^\infty f(x) \log_2 f(x)\, dx$$

**Normalverteilung** $\mathcal{N}(\mu, \sigma^2)$: $h(X) = \frac{1}{2}\log_2(2\pi e \sigma^2)$

**Exponentialverteilung** $\text{Exp}(\lambda)$: $h(X) = 1 - \log_2 \lambda$

**Gleichverteilung** $\mathcal{U}(a,b)$: $h(X) = \log_2(b-a)$

Die Normalverteilung maximiert die differentielle Entropie bei gegebener Varianz.

---

## Klassen und Methoden

### `ShannonEntropy`

Shannon-Entropie, Verbund- und bedingte Entropie, gegenseitige Information.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `entropy` | `(probabilities, base=2) → float` | Shannon-Entropie $H(X)$ |
| `joint_entropy` | `(joint_probs, base=2) → float` | Verbundentropie $H(X,Y)$ |
| `conditional_entropy` | `(joint_probs, base=2) → float` | Bedingte Entropie $H(Y|X)$ |
| `mutual_information` | `(joint_probs, base=2) → float` | Gegenseitige Information $I(X;Y)$ |
| `entropy_maximum` | `(n, base=2) → float` | Maximale Entropie $\log_2 n$ |
| `entropy_binary` | `(p, base=2) → float` | Binäre Entropie $h(p)$ |

### `KLDivergence`

Kullback-Leibler-Divergenz und verwandte Informationsabstandsmaße.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `kl_divergence` | `(p, q, base=2) → float` | $D_\mathrm{KL}(P \| Q)$ |
| `js_divergence` | `(p, q, base=2) → float` | Jensen-Shannon-Divergenz |
| `cross_entropy` | `(p, q, base=2) → float` | Kreuzentropie $H(P, Q)$ |
| `renyi_entropy` | `(probabilities, alpha, base=2) → float` | Rényi-Entropie $H_\alpha$ |
| `tsallis_entropy` | `(probabilities, q, base=2) → float` | Tsallis-Entropie $S_q$ |
| `total_variation_distance` | `(p, q) → float` | Total-Variationsabstand $d_\mathrm{TV}$ |

### `ChannelCapacity`

Kanalkapazität und Informationsübertragung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `binary_symmetric_channel_capacity` | `(p) → float` | $C_\mathrm{BSC} = 1 - h(p)$ |
| `binary_erasure_channel_capacity` | `(epsilon) → float` | $C_\mathrm{BEC} = 1 - \varepsilon$ |
| `awgn_channel_capacity` | `(snr_db) → float` | $C_\mathrm{AWGN} = \log_2(1 + \mathrm{SNR})$ |
| `channel_capacity_general` | `(transition_matrix, ...) → float` | Kapazität für allgemeinen diskreten Kanal |
| `blahut_arimoto` | `(transition_matrix, max_iter, tol) → dict` | Blahut-Arimoto-Iteration |
| `noisy_channel_coding_theorem_demo` | `(capacity, rate) → dict` | Shannon'sches Kanalcodierungstheorem (Demo) |

### `SourceCoding`

Quellencodierung und Datenkompression.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `huffman_code` | `(symbols, probabilities) → dict` | Huffman-Codebaum und -codes |
| `huffman_average_length` | `(huffman_dict, probabilities) → float` | Mittlere Codelänge $\bar{L}$ |
| `arithmetic_coding_demo` | `(symbols, probabilities, message) → dict` | Arithmetische Codierung (Demo) |
| `lempel_ziv_complexity` | `(sequence) → int` | LZ-Komplexität einer Folge |
| `kolmogorov_complexity_demo` | `(n) → dict` | Approximation der Kolmogorov-Komplexität |
| `source_coding_theorem_bound` | `(probabilities, base=2) → dict` | Shannon-Schranken für Quellencodierung |

### `ErrorCorrection`

Grundlagen der Fehlerkorrektur (einfache Codes, Schranken).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `hamming_distance` | `(codeword1, codeword2) → int` | Hamming-Distanz $d_H$ |
| `hamming_weight` | `(codeword) → int` | Hamming-Gewicht $w_H$ |
| `minimum_distance` | `(code) → int` | Minimaldistanz eines Codes |
| `singleton_bound` | `(n, k) → int` | Singleton-Schranke $d \leq n-k+1$ |
| `hamming_bound` | `(n, t) → int` | Hamming-Schranke (Kugelpackung) |
| `hamming_code_74` | `() → dict` | Vollständiger (7,4)-Hamming-Code |
| `parity_check_code` | `(n) → dict` | Paritätsprüfcode der Länge $n$ |
| `repetition_code` | `(message, k) → list[int]` | $k$-fache Wiederholungscodierung |

### `DifferentialEntropy`

Differentielle Entropie für kontinuierliche Verteilungen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `differential_entropy_normal` | `(sigma, base=2) → float` | $h(\mathcal{N}) = \frac{1}{2}\log_2(2\pi e\sigma^2)$ |
| `differential_entropy_exponential` | `(lambda_, base=2) → float` | $h(\mathrm{Exp}) = 1 - \log_2\lambda$ |
| `differential_entropy_uniform` | `(a, b, base=2) → float` | $h(\mathcal{U}) = \log_2(b-a)$ |
| `mutual_information_gaussian` | `(rho, base=2) → float` | $I = -\frac{1}{2}\log_2(1-\rho^2)$ für bivariates Gauss |

### Standalone-Funktionen

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `surprisal` | `(p, base=2) → float` | Informationsgehalt $I(x) = -\log_2 p$ |
| `entropy_rate_markov` | `(transition_matrix, stationary_dist) → float` | Entropierate $H = -\sum_i \pi_i \sum_j T_{ij}\log T_{ij}$ |
| `data_processing_inequality_demo` | `() → dict` | Demo: $I(X;Z) \leq I(X;Y)$ für Markov-Kette $X-Y-Z$ |

---

## Wichtige Theoreme

**Noisy Channel Coding Theorem (Shannon 1948):** Für jede Rate $R < C$ existieren Codes mit beliebig kleiner Fehlerwahrscheinlichkeit. Für $R > C$ ist fehlerfreie Übertragung unmöglich.

**Quellencodierungssatz:** Jede Quelle mit Entropie $H$ kann nicht verlustfrei unter $H$ Bit/Symbol komprimiert werden; Huffman erreicht $[H, H+1)$.

**Datenverarbeitungs-Ungleichung:** Für eine Markov-Kette $X \to Y \to Z$ gilt: $I(X;Z) \leq I(X;Y)$. Verarbeitung kann keine Information hinzufügen.

**Informations-Ungleichung (Gibbs):** $D_\mathrm{KL}(P \| Q) \geq 0$, Gleichheit genau bei $P = Q$.

---

## Beispiele

```python
from information_theory import (ShannonEntropy, KLDivergence, ChannelCapacity,
                                 SourceCoding, ErrorCorrection, surprisal)

# -- Entropie einer Münze --
se = ShannonEntropy()
print(se.entropy([0.5, 0.5]))          # 1.0 Bit (maximale Entropie)
print(se.entropy([0.9, 0.1]))          # ≈ 0.469 Bit
print(se.entropy_binary(0.1))          # h(0.1) ≈ 0.469 Bit

# -- KL-Divergenz --
kl = KLDivergence()
P = [0.4, 0.3, 0.2, 0.1]
Q = [0.25, 0.25, 0.25, 0.25]
print(kl.kl_divergence(P, Q))         # ≈ 0.124 Bit
print(kl.js_divergence(P, Q))         # ≤ 1 Bit, symmetrisch

# -- Kanalkapazität --
cc = ChannelCapacity()
print(cc.binary_symmetric_channel_capacity(0.1))  # ≈ 0.531 Bit
print(cc.awgn_channel_capacity(10))               # log₂(1+10) ≈ 3.459 Bit

# -- Huffman-Code --
sc = SourceCoding()
symbols = ["a", "b", "c", "d"]
probs   = [0.5, 0.25, 0.125, 0.125]
codes   = sc.huffman_code(symbols, probs)
L_bar   = sc.huffman_average_length(codes, probs)
print(codes)    # {'a': '0', 'b': '10', 'c': '110', 'd': '111'}
print(L_bar)    # 1.75 Bit (optimal, H = 1.75 Bit)

# -- Informationsgehalt --
print(surprisal(0.5))    # 1.0 Bit
print(surprisal(0.125))  # 3.0 Bit

# -- (7,4)-Hamming-Code --
ec = ErrorCorrection()
hc = ec.hamming_code_74()
print(hc["n"], hc["k"], hc["d"])  # 7, 4, 3
```

---

## Tests

**Testdatei:** `tests/test_information_theory.py`
**Abdeckung:** Entropie-Axiome, KL-Nicht-Negativität, Symmetrie von JS, Kanalkapazität-Schranken, Huffman-Optimalität, Hamming-Schranken, differentielle Entropie-Formeln

---

## Implementierungshinweise

- **Entropie:** Terme mit $p_i = 0$ werden als $0 \cdot \log 0 = 0$ behandelt (stetige Erweiterung).
- **Blahut-Arimoto:** Konvergenz in $O(1/\varepsilon)$ Iterationen; bei schlechter konditionierten Kanälen kann $10^{-6}$-Toleranz mehr als 1000 Iterationen erfordern.
- **Huffman:** Implementiert als Min-Heap über `heapq`; bei Gleichwahrscheinlichkeit nicht eindeutig.
- **LZ-Komplexität:** Martin-Algorithmus; misst strukturelle Komplexität, nicht direkt Kompressionsrate.
- **Differentielle Entropie:** Kann negativ sein (im Gegensatz zur diskreten Entropie) – kein Informationsgehalt-Paradox, da sie nicht unitär ist.
