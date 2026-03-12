# schur_numbers.py – Schur-Zahlen und sum-freie Partitionen

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Letzte Änderung**: 2026-03-12

---

## Mathematischer Hintergrund

### Schur-Zahlen S(k)

**Definition** (Schur 1916): $S(k)$ ist die größte natürliche Zahl $n$, sodass die Menge $\{1, 2, \ldots, n\}$ in $k$ Klassen aufgeteilt werden kann, ohne dass eine Klasse ein *Schur-Tripel* $(a, b, c)$ mit $a + b = c$ enthält.

Eine solche Partition heißt **sumfreie Partition** (sum-free partition). Eine einzelne Menge $M$ heißt **sumfrei**, falls es kein $a, b, c \in M$ mit $a + b = c$ gibt (Beachte: $a = b$ ist erlaubt).

### Bekannte Werte

| $k$ | $S(k)$ | Beweis |
|-----|--------|--------|
| 1 | 1 | Trivial: $\{1\}$ sumfrei, $\{1,2\}$ nicht ($1+1=2$) |
| 2 | 4 | Partition: $\{1,4\}, \{2,3\}$ |
| 3 | 13 | Baumert (1965) |
| 4 | 44 | Baumert (1965) |
| 5 | 160 | Heule (2017, SAT-Beweis) |
| 6 | ? | $536 \leq S(6) \leq 1836$ |

### Ramsey-Theorie-Kontext

Schur bewies: Für jede $k$-Färbung der positiven ganzen Zahlen existiert ein einfarbiges Schur-Tripel, wenn $n > S(k)$. Dies ist ein **Ramsey-artiger Satz**: Ab einer gewissen Größe können Strukturen nicht vermieden werden.

### Rekursive Schranke

$$S(k) \geq 3 \cdot S(k-1) + 1$$

Diese Schranke liefert $S(6) \geq 3 \cdot 160 + 1 = 481$. Die bekannte untere Schranke $S(6) \geq 536$ ist besser und stammt aus einer expliziten Computer-Partition.

### SAT-Problem-Verbindung

Das Problem "$\{1,\ldots,n\}$ ist $k$-färbbar" ist ein SAT-Problem:
- **Variablen**: $x_{i,c} \in \{0,1\}$ für $i \in \{1,\ldots,n\}$, $c \in \{0,\ldots,k-1\}$
- **Klauseln**: At-least-one, At-most-one, Schur-Tripel-Verbote

Heule et al. (2017) bewiesen $S(5) = 160$ mit einem SAT-Solver (21 GB Zertifikat).

---

## Klassen-API

### `SchurZahlen`

#### `ist_schursch(partition: Partition) -> bool`

Prüft ob eine Partition sumfrei ist. Eine Partition ist ein `List[Set[int]]`.

```python
sz = SchurZahlen()
p = [{1, 4}, {2, 3}]
assert sz.ist_schursch(p) is True   # S(2)=4: gültig
```

#### `generiere_schursche_partition(n, k, methode, seed) -> Optional[Partition]`

Versucht $\{1,\ldots,n\}$ in $k$ sumfreie Klassen zu partitionieren.

- **greedy**: Nicht vollständig, aber schnell. Kann `None` zurückgeben obwohl Lösung existiert.
- **backtracking**: Vollständig. Findet immer eine Lösung wenn eine existiert.

```python
p = sz.generiere_schursche_partition(44, 4, methode="backtracking")
assert sz.ist_schursch(p)
```

#### `verifiziere_schur_zahl(k, kandidat_n) -> bool`

Prüft via Backtracking ob $S(k) \geq \text{kandidat\_n}$.

#### `bekannte_partitionen() -> Dict[int, Partition]`

Gibt verifizierte Partitionen für $S(1)$ bis $S(5)$ zurück.

#### `partition_s5() -> Partition`

Konstruiert eine sumfreie Partition von $\{1,\ldots,160\}$ in 5 Klassen (zeigt $S(5) \geq 160$). Methode: Rekursive Erweiterung der $S(4)=44$-Partition.

#### `untere_schranke_s6() -> Tuple[int, Optional[Partition]]`

Konstruiert eine Partition von $\{1,\ldots,536\}$ (oder möglichst viel) in 6 Klassen. Zeigt $S(6) \geq n$.

#### `schranken_analyse() -> Dict`

Gibt bekannte Schranken für $S(6)$ zurück: $536 \leq S(6) \leq 1836$.

#### `sat_encoding(n, k) -> Dict`

Beschreibt die SAT-Kodierung für "$\{1,\ldots,n\}$ ist $k$-färbbar?" mit Variablenanzahl und Klausel-Statistiken.

#### `monte_carlo_partition(n, k, versuche, seed) -> Dict`

Probabilistische Suche: Mehrere zufällige Greedy-Versuche. Nützlich für große $n$ wo Backtracking zu langsam ist.

---

## Interne Hilfsfunktionen

### `_verletzt_sumfreiheit(zahl, klasse) -> bool`

Prüft alle drei Tripel-Typen die entstehen wenn `zahl` zu `klasse` hinzugefügt wird:
- **Typ A**: $\text{zahl} + a = c$ mit $a, c \in \text{klasse}$
- **Typ B**: $a + b = \text{zahl}$ mit $a, b \in \text{klasse}$
- **Typ C**: $2 \cdot \text{zahl} \in \text{klasse}$ (doppelter Summand)

### `_ist_sumfrei(menge) -> bool`

Prüft ob eine Menge sumfrei ist (kein Tripel $a + b = c$).

### `_hat_schur_tripel(menge) -> Optional[Tuple]`

Gibt das erste Schur-Tripel zurück oder `None`.

---

## Bekannte Einschränkungen

- **Greedy nicht vollständig**: Greedy kann für lösbare Instanzen `None` zurückgeben. Verwende Backtracking für exakte Ergebnisse.
- **Backtracking exponentiell**: Für $n \geq 45, k = 4$ kann Backtracking sehr lange dauern.
- **partition_s5**: Die Methode kann mehr als 5 Klassen zurückgeben wenn Greedy für $\{134,\ldots,160\}$ eine neue Klasse öffnet.
- **SAT für S(6)**: Der SAT-Beweis würde $\sim 10^{12}$ Operationen erfordern.

---

## Quellen

- I. Schur (1916): *Über die Kongruenz $x^m + y^m \equiv z^m \pmod{p}$*, Jahresber. Dtsch. Math.-Ver.
- L.D. Baumert (1965): *A note on Schur's problem*, unpublished
- M. Heule (2017): *Schur Number Five*, arXiv:1711.08076
- S. Eliahou et al. (2011): *Schur numbers and the Ramsey multiplicity*
