# Modul: Kodierungstheorie (`coding_theory.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `coding_theory.py` implementiert die algebraische Kodierungstheorie sowie die Theorie formaler Sprachen und Automaten. Es behandelt Galois-Körper, lineare Codes, Hamming-Codes, zyklische Codes, BCH-Codes (mit Berlekamp-Massey-Algorithmus), deterministische und nichtdeterministische endliche Automaten, kontextfreie Grammatiken, Kellerautomaten und die Berechenbarkeitstheorie.

**Hinweis:** Hamming-Distanz, Singleton-Bound, der (7,4)-Hamming-Code und Paritätscodes befinden sich in `information_theory.py` (Klasse `ErrorCorrection`). Dieses Modul ergänzt mit algebraischer Körpertheorie und formalen Sprachen.

---

## Mathematischer Hintergrund

### Lineare Codes $[n, k, d]$

Ein **linearer Code** $C$ über $\mathrm{GF}(2)$ ist ein $k$-dimensionaler Unterraum von $\mathrm{GF}(2)^n$:

$$n = \text{Codewortlänge}, \quad k = \text{Informationsbits}, \quad d = \text{Minimaldistanz}$$

**Codierung:** $\mathbf{c} = \mathbf{m} \cdot G$ mit Generatormatrix $G \in \mathrm{GF}(2)^{k \times n}$

**Prüfmatrix:** $H \in \mathrm{GF}(2)^{(n-k) \times n}$, $H \cdot \mathbf{c}^T = \mathbf{0}$ für alle Codewörter

**Syndrom:** $\mathbf{s} = H \cdot \mathbf{r}^T$ für ein empfangenes Wort $\mathbf{r}$; $\mathbf{s} = \mathbf{0}$ bedeutet fehlerfrei

### Hamming-Code $[2^r - 1,\; 2^r - r - 1,\; 3]$

Perfekter 1-fehlerkorrigierender Code. Die Prüfmatrix $H$ hat als Spalten alle $2^r - 1$ Bitvektoren $\neq \mathbf{0}$. Das Syndrom identifiziert direkt die Fehlerposition.

**Hamming-Schranke (Kugelpackungs-Schranke):**
$$\sum_{i=0}^t \binom{n}{i} \leq 2^{n-k}$$

### Zyklische Codes

Ein Code $C$ heißt **zyklisch**, wenn mit jedem Codewort $(c_0, c_1, \ldots, c_{n-1})$ auch die zyklische Verschiebung $(c_{n-1}, c_0, \ldots, c_{n-2})$ in $C$ liegt.

**Darstellung:** Durch Generatorpolynom $g(x)$ mit $g(x) \mid (x^n - 1)$:
$$\deg g(x) = n - k, \quad \text{Codewort: } c(x) = m(x) \cdot g(x)$$

**Syndrompolynom:** $s(x) = r(x) \bmod g(x)$; $s(x) = 0$ bedeutet kein Fehler.

### BCH-Codes und Reed-Solomon

**BCH-Code** (Bose–Chaudhuri–Hocquenghem): Zyklischer Code mit $t$-Fehlerkorrektur. Das Generatorpolynom ist das kgV der Minimalpolynome der ersten $2t$ Potenzen eines primitiven Elements $\alpha$.

**Reed-Solomon-Code:** Spezieller BCH-Code über $\mathrm{GF}(q)$ mit $n = q-1$, $k$ Informationssymbolen, Minimaldistanz $d = n - k + 1$ (erreicht Singleton-Schranke).

**Berlekamp-Massey-Algorithmus:** Findet das kürzeste lineare Rückkopplungsschieberegister (LFSR), das eine gegebene Folge erzeugt. Laufzeit $O(n^2)$.

### Galois-Körper $\mathrm{GF}(p^n)$

**Körpergröße:** $|\mathrm{GF}(p^n)| = p^n$ ($p$ Primzahl)

**Multiplikative Gruppe:** $\mathrm{GF}(p^n)^* = \mathrm{GF}(p^n) \setminus \{0\}$ ist zyklisch der Ordnung $p^n - 1$.

**Primitives Element** $g$: erzeugt die gesamte multiplikative Gruppe, $\text{ord}(g) = p^n - 1$.

**Diskreter Logarithmus:** Finde $x$ mit $g^x \equiv a \pmod{p}$ – algorithmisch schweres Problem (Baby-Step-Giant-Step).

### Formale Sprachen und Chomsky-Hierarchie

| Typ | Grammatiktyp | Automat | Sprachen |
|-----|-------------|---------|---------|
| Typ 0 | Unrestriktiert | Turing-Maschine | Rekursiv aufzählbar |
| Typ 1 | Kontextsensitiv | Linear-bounded TM | Kontextsensitiv |
| Typ 2 | Kontextfrei | Kellerautomat | Kontextfreie Sprachen |
| Typ 3 | Regulär | Endlicher Automat | Reguläre Sprachen |

### Berechenbarkeitstheorie

**Primitiv-rekursive Funktionen:** Aufgebaut aus Nullfunktion, Nachfolgerfunktion, Projektion via Komposition und primitiver Rekursion.

**$\mu$-Operator:** Minimierungsoperator $\mu y[f(x, y) = 0]$ – kleinste $y$ mit $f(x,y)=0$.

**Church-Turing-These:** Jede intuitiv berechenbare Funktion ist Turing-berechenbar.

---

## Klassen und Methoden

### `GaloisField`

Galois-Körper $\mathrm{GF}(p^n)$ mit Körperarithmetik.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `elements` | `() → list[int]` | Alle $p^n$ Körperelemente |
| `add` | `(a, b: int) → int` | Addition $(a + b) \bmod p$ |
| `multiply` | `(a, b: int) → int` | Multiplikation $(a \cdot b) \bmod p$ |
| `primitive_element` | `() → int` | Kleinstes primitives Element $g$ |
| `discrete_log_gf` | `(base, element: int) → int` | Diskreter Logarithmus $\log_g(\text{element})$ |
| `is_primitive_poly` | `(coeffs: list[int]) → bool` | Prüft ob Polynom über $\mathrm{GF}(p)$ primitiv ist |

### `LinearCode`

Linearer $[n,k,d]$-Code über $\mathrm{GF}(2)$.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `encode` | `(message) → np.ndarray` | Codierung $\mathbf{c} = \mathbf{m} \cdot G$ |
| `parity_check_matrix` | `() → np.ndarray` | Berechnet $H$ aus $G$ (systematische Form) |
| `syndrome` | `(received) → np.ndarray` | Syndrom $\mathbf{s} = H \cdot \mathbf{r}^T$ |
| `minimum_distance` | `() → int` | Minimaldistanz (Brute-Force über alle Codeworte) |
| `error_correct` | `(received) → np.ndarray` | Fehlerkorrektur per Syndromdekodierung |
| `is_systematic` | `() → bool` | Prüft ob Code in systematischer Form vorliegt |
| `dual_code` | `() → LinearCode` | Dualer Code $C^\perp$ |

### `HammingCode`

Hamming-Code $[2^r-1, 2^r-r-1, 3]$ für gegebenes $r$.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `parameters` | `() → (n, k, d)` | Codeparameter $(n, k, d)$ |
| `parity_check_matrix` | `() → np.ndarray` | Prüfmatrix $H$ (Spalten = alle Bitvektoren ≠ 0) |
| `generator_matrix` | `() → np.ndarray` | Generatormatrix $G$ (aus $H$ berechnet) |
| `encode` | `(message) → np.ndarray` | Codierung einer $k$-Bit-Nachricht |
| `decode` | `(received) → dict` | Dekodierung mit Fehlerkorrektur und -position |

### `CyclicCode`

Zyklischer Code mit Generatorpolynom.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `encode_poly` | `(message_poly) → list[int]` | Codierung als Polynom $c(x) = m(x) \cdot g(x)$ |
| `syndrome_poly` | `(received_poly) → list[int]` | Syndrompolynom $s(x) = r(x) \bmod g(x)$ |
| `cyclic_shift` | `(codeword) → list[int]` | Zyklische Verschiebung eines Codewortes |
| `is_cyclic` | `(code_words) → bool` | Prüft ob eine Menge von Wörtern zyklisch abgeschlossen ist |

### `BCHCode`

BCH-Codes und Reed-Solomon (als Demo).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `bch_generator_poly` | `(n, t, q=2) → dict` | BCH-Generatorpolynom für $t$-Fehlerkorrektur |
| `reed_solomon_demo` | `(n, k, q) → dict` | Reed-Solomon-Parameter und Eigenschaften |
| `berlekamp_massey` | `(sequence) → dict` | Kürzestes LFSR für gegebene Folge |

### `DFA`

Deterministischer Endlicher Automat.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `accepts` | `(word) → bool` | Prüft ob Wort akzeptiert wird |
| `minimize` | `() → dict` | Minimierung per Äquivalenzklassen (Hopcroft/Moore) |

### `NFA`

Nichtdeterministischer Endlicher Automat.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `accepts` | `(word) → bool` | Prüft ob Wort akzeptiert wird (Teilmengenkonstruktion) |

### `RegularExpression`

Reguläre Ausdrücke über Python-`re`.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `matches` | `(word: str) → bool` | Prüft ob Wort auf regulären Ausdruck passt |

### `ContextFreeGrammar`

Kontextfreie Grammatik $G = (V, \Sigma, P, S)$.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `generates` | `(word, max_depth=10) → bool` | Prüft per CYK/Backtracking ob Wort ableitbar |
| `is_ambiguous_demo` | `() → dict` | Demo zur Mehrdeutigkeit einer CFG |

### `PushdownAutomaton`

Kellerautomat (nichtdeterministisch).

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `accepts_demo` | `(word) → dict` | Simuliert PDA und gibt Konfigurationsfolge zurück |

### `FormalLanguages`

Überblicksklasse für formale Sprachen.

Enthält statische Methoden zur Illustration der Chomsky-Hierarchie, Pumpinglemma-Demo, Abschlusseigenschaften regulärer Sprachen.

### `ComputabilityTheory`

Berechenbarkeitstheorie – Ergänzung zu `recursion_theory.py`.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `primitive_recursive_demo` | `() → dict` | Demo primitiv-rekursiver Funktionen (add, mul, fac, fib) |
| `mu_recursive_demo` | `() → dict` | Demo des $\mu$-Operators (integer square root) |
| `church_turing_thesis_demo` | `() → dict` | Illustration der Church-Turing-These |
| `undecidable_problems` | `() → dict` | Liste unentscheidbarer Probleme (Halteprobleme u.a.) |

---

## Beispiele

```python
from coding_theory import GaloisField, LinearCode, HammingCode, CyclicCode, DFA
import numpy as np

# -- Galois-Körper GF(7) --
gf = GaloisField(p=7)
print(gf.elements())          # [0, 1, 2, 3, 4, 5, 6]
print(gf.add(5, 4))           # (5+4) mod 7 = 2
print(gf.primitive_element()) # 3 (kleinste primitive Wurzel mod 7)

# -- Linearer Code: (6,3)-Code --
G = np.array([[1,0,0,1,1,0],
              [0,1,0,1,0,1],
              [0,0,1,0,1,1]], dtype=int)
code = LinearCode(G)
c = code.encode([1, 0, 1])    # [1, 0, 1, 1, 0, 1]
r = [1, 0, 1, 0, 0, 1]        # Fehler an Position 3
s = code.syndrome(r)
corrected = code.error_correct(r)

# -- Hamming-Code r=3: [7,4,3]-Code --
hc = HammingCode(r=3)
print(hc.parameters())        # (7, 4, 3)
cw = hc.encode([1, 0, 1, 1])
result = hc.decode(cw)        # fehlerfrei
print(result["corrected"])

# -- Zyklischer Code --
cc = CyclicCode(n=7, generator_poly=[1, 0, 1, 1])  # g(x) = x³+x+1
encoded = cc.encode_poly([1, 1, 0, 1])
print(cc.syndrome_poly(encoded))  # [0, 0, 0] = fehlerfrei

# -- DFA: Sprache {w | w enthält "ab"} --
dfa = DFA(
    states={0, 1, 2},
    alphabet={"a", "b"},
    transitions={(0,"a"):1, (0,"b"):0, (1,"a"):1, (1,"b"):2, (2,"a"):1, (2,"b"):0},
    initial=0,
    accepting={2}
)
print(dfa.accepts("aab"))  # True
print(dfa.accepts("ba"))   # False
```

---

## Tests

**Testdatei:** `tests/test_coding_theory.py`
**Abdeckung:** Galois-Körper-Arithmetik, Kodierung/Dekodierung, Syndromkorrektur, Hamming-Parameter, zyklische Eigenschaften, Automaten-Akzeptanz, Grammatikableitungen, Berechenbarkeitsdemos

---

## Implementierungshinweise

- **GF(2)-Arithmetik:** Addition und Multiplikation per XOR bzw. AND; `_row_reduce_gf2()` implementiert Gauss-Jordan-Elimination über GF(2).
- **Parity-Check-Matrix:** Wird aus der Generatormatrix in systematischer Form $(I_k | P)$ als $H = (-P^T | I_{n-k})$ abgeleitet; über GF(2) ist $-1 = 1$.
- **Berlekamp-Massey:** Iterativer Algorithmus; Ausgabe enthält LFSR-Polynom und -Länge.
- **DFA-Minimierung:** Myhill-Nerode-Äquivalenzklassen per iterativer Verfeinerung.
- **NFA:** Simulation per Potenzmengenkonstruktion (Teilmengenkonstruktion) on-the-fly.
