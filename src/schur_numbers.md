# schur_numbers.py — Schur-Zahlen, sum-freie Partitionen und Ramsey-Theorie

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12
**Modul:** `src/schur_numbers.py`

---

## Inhaltsverzeichnis

1. [Mathematischer Hintergrund](#mathematischer-hintergrund)
2. [Klassen- und Methodenübersicht](#klassen--und-methodenübersicht)
3. [Wichtigste Algorithmen](#wichtigste-algorithmen)
4. [Beispielanwendungen](#beispielanwendungen)
5. [SAT-Kodierung](#sat-kodierung)
6. [Bekannte Schur-Zahlen](#bekannte-schur-zahlen)

---

## Mathematischer Hintergrund

### Definition der Schur-Zahlen

**Issai Schur (1916)** definierte:

> $S(k)$ ist die größte natürliche Zahl $n$, sodass die Menge $\{1, 2, \ldots, n\}$ in $k$ Klassen aufgeteilt werden kann, ohne dass eine Klasse drei Zahlen $a, b, c$ mit
> $$a + b = c$$
> enthält. Eine solche Klasse heißt **sum-free** (summenfrei).

Mit anderen Worten: $S(k)$ ist die größte Zahl $n$, für die eine **sum-freie $k$-Färbung** von $\{1, \ldots, n\}$ existiert.

### Schur-Theorem

**Schur's Theorem (1916):** Für jede $k$-Färbung der positiven ganzen Zahlen existiert eine einfarbige Lösung von $x + y = z$, sobald $n > S(k)$.

Äquivalent: Für jedes $k$ gibt es eine endliche Schranke $S(k)$.

**Ursprüngliche Motivation:** Schur wollte zeigen, dass der Satz von Fermat (also $x^n + y^n = z^n$ hat keine positiven ganzzahligen Lösungen) nicht durch Modularrechnung beweisbar ist: Für jede Primzahl $p > S(k)$ gibt es eine einfarbige Lösung von $x + y \equiv z \pmod{p}$.

### Sum-freie Mengen

Eine Menge $M \subseteq \mathbb{N}$ heißt **sum-frei**, falls:
$$\forall a, b, c \in M: a + b \neq c$$

Beachte: $a = b$ ist erlaubt ($a + a = 2a$, also muss $2a \notin M$ gelten).

**Beispiele:**
- $\{2, 3\}$ ist sum-frei: $2+2=4\notin M$, $2+3=5\notin M$, $3+3=6\notin M$ — ja.
- $\{1, 2, 3\}$ ist **nicht** sum-frei: $1+2=3\in M$.

### Bekannte Schur-Zahlen

| $k$ | $S(k)$ | Status |
|-----|--------|--------|
| 1 | 1 | exakt |
| 2 | 4 | exakt |
| 3 | 13 | exakt |
| 4 | 44 | exakt |
| 5 | 160 | exakt (Heule 2017, SAT-Beweis) |
| 6 | ? | $536 \leq S(6) \leq 1836$ (offen) |

### Schranken

**Rekursive untere Schranke:**
$$S(k) \geq 3 \cdot S(k-1) + 1$$

Denn: Ist $(K_1, \ldots, K_{k-1})$ eine $(k-1)$-Partition von $\{1,\ldots,n\}$, so ist die Menge $\{n+1,\ldots,2n+1\}$ sum-frei (für $a,b \geq n+1$ gilt $a+b \geq 2n+2 > 2n+1$). Damit erhält man eine $k$-Partition von $\{1,\ldots,3n+1\}$.

**Obere Schranke (Schur 1916):**
$$S(k) \leq \lfloor k! \cdot e \rfloor$$

### Verbindung zur Ramsey-Theorie

Schur-Zahlen sind ein Spezialfall der **Ramsey-Zahlen** $R(3,3,\ldots,3)$ (mit $k$ Argumenten). Es gilt:
$$S(k) \leq R_k(3) - 2$$

wobei $R_k(3)$ die $k$-dimensionale Ramsey-Zahl für Dreiecke ist.

---

## Klassen- und Methodenübersicht

### Typ-Alias

```python
Partition = List[Set[int]]  # Eine Liste von Mengen (eine Menge pro Farbe)
```

### Modul-interne Hilfsfunktionen

| Funktion | Beschreibung |
|----------|--------------|
| `_verletzt_sumfreiheit(zahl, klasse)` | Prüft ob Hinzufügen von `zahl` die Sumfreiheit verletzt |
| `_ist_sumfrei_partition(partition)` | Prüft ob alle Klassen sum-frei sind |
| `_ist_sumfrei(menge)` | Prüft ob eine einzelne Menge sum-frei ist |
| `_hat_schur_tripel(menge)` | Gibt erstes Schur-Tripel $(a,b,c)$ mit $a+b=c$ zurück |

### Klasse `SchurZahlen`

Klassenattribut:
- `BEKANNTE_SCHUR_ZAHLEN: Dict[int,int]` — $S(1)=1, S(2)=4, S(3)=13, S(4)=44, S(5)=160$

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `ist_schursch` | `(partition: Partition) → bool` | Prüft ob alle Klassen sum-frei sind |
| `generiere_schursche_partition` | `(n, k, methode?, seed?) → Optional[Partition]` | Hauptmethode für Partitionierung |
| `verifiziere_schur_zahl` | `(k, kandidat_n) → bool` | Prüft ob $S(k) \geq$ kandidat_n |
| `bekannte_partitionen` | `() → Dict[int, Partition]` | Gibt bekannte Partitionen für $S(1)$–$S(5)$ zurück |
| `partition_s5` | `() → Partition` | Explizite Partition von $\{1,\ldots,160\}$ in 5 Klassen |
| `untere_schranke_s6` | `() → Tuple[int, Optional[Partition]]` | Konstruiert sum-freie 6-Partition |
| `schranken_analyse` | `() → Dict` | Zusammenfassung der bekannten $S(6)$-Schranken |
| `sat_encoding` | `(n, k) → Dict` | Beschreibt SAT-Kodierung des Problems |
| `monte_carlo_partition` | `(n, k, versuche?, seed?) → Dict` | Probabilistische Partitionssuche |

---

## Wichtigste Algorithmen

### 1. Greedy-Algorithmus

**Idee:** Weise jede Zahl $1, 2, \ldots, n$ der ersten Klasse zu, die nach der Aufnahme noch sum-frei bleibt.

```
Für zahl = 1 bis n:
    Für jede Klasse K:
        Falls _verletzt_sumfreiheit(zahl, K) = False:
            K.add(zahl)
            break
    Falls keine Klasse passt: return None
```

**Vorteil:** Sehr schnell, $O(n \cdot k)$ Schrittzahl.
**Nachteil:** Nicht vollständig — kann auch bei lösbaren Instanzen None zurückgeben.

### 2. Backtracking-Algorithmus

**Idee:** Systematische Tiefensuche über alle Farbzuweisungen.

```
farbe[0..n-1] = [-1, ..., -1]

backtrack(pos):
    falls pos == n: return True
    für c = 0 bis k-1:
        falls verletzt_sumfreiheit(pos, c) = False:
            farbe[pos] = c
            falls backtrack(pos+1): return True
            farbe[pos] = -1
    return False
```

**Symmetriebrechung:** Zahl 1 wird immer Farbe 0 zugewiesen (reduziert Suchraum um Faktor $k$).

**Komplexität:** Im schlimmsten Fall $O(k^n)$, praktisch für $n \leq 50$ effizient.

### 3. Verletzungsprüfung `_verletzt_sumfreiheit`

Prüft drei Fälle wenn `zahl` zu `klasse` hinzugefügt würde:

- **Typ A:** $\text{zahl} + a = c$ mit $a, c \in \text{klasse}$ (zahl als erster Summand)
- **Typ B:** $a + b = \text{zahl}$ mit $a, b \in \text{klasse}$ (zahl als neue Summe)
- **Typ C:** $\text{zahl} + \text{zahl} = c$ mit $c \in \text{klasse}$ (beide Summanden gleich)

### 4. Rekursive Konstruktion für $S(5)$

Aus einer $4$-Partition $(K_1,K_2,K_3,K_4)$ von $\{1,\ldots,44\}$ wird eine $5$-Partition von $\{1,\ldots,133\}$ konstruiert:

1. **Mittlere sum-freie Menge:** $K_5 = \{45, \ldots, 89\}$ (sum-frei: $a+b \geq 90 > 89$ für $a,b \geq 45$)
2. **Verschobene Kopie:** $K_i' = \{x + 89 : x \in K_i\}$ für $i=1,2,3,4$ — deckt $\{90,\ldots,133\}$
3. **Vereinigung:** $K_i \cup K_i'$ für $i=1,2,3,4$, plus $K_5$

Für $\{134,\ldots,160\}$ wird Greedy auf die bestehenden 5 Klassen angewendet.

### 5. Monte-Carlo-Suche

Mehrfache zufällige Greedy-Läufe mit zufällig permutierter Reihenfolge der Zahlen. Für jede Zahl wird zufällig unter allen gültigen Klassen gewählt.

---

## Beispielanwendungen

### Partition verifizieren

```python
from schur_numbers import SchurZahlen

sz = SchurZahlen()

# S(2) = 4: Partition {{1,4}, {2,3}}
partition_s2 = [{1, 4}, {2, 3}]
print(sz.ist_schursch(partition_s2))  # True

# S(3) = 13: bekannte Partition
partitionen = sz.bekannte_partitionen()
k3 = partitionen[3]
print(sz.ist_schursch(k3))            # True
```

### Partition erzeugen (Greedy)

```python
# {1,...,13} in 3 sum-freie Klassen
p = sz.generiere_schursche_partition(13, 3, methode="greedy")
if p:
    print(sz.ist_schursch(p))         # True
    for i, k in enumerate(p):
        print(f"Klasse {i+1}: {sorted(k)}")
```

### Partition erzeugen (Backtracking)

```python
# Beweise S(4) ≥ 44
p4 = sz.generiere_schursche_partition(44, 4, methode="backtracking")
print(p4 is not None)                 # True
print(sz.ist_schursch(p4))            # True

# S(4) ≥ 45 ist falsch (Backtracking liefert None)
p45 = sz.generiere_schursche_partition(45, 4, methode="backtracking")
print(p45 is None)                    # True  →  S(4) < 45  →  S(4) = 44
```

### Schur-Zahl verifizieren

```python
print(sz.verifiziere_schur_zahl(3, 13))   # True  (S(3) ≥ 13)
print(sz.verifiziere_schur_zahl(3, 14))   # False (S(3) < 14)
```

### Schranken-Analyse für S(6)

```python
analyse = sz.schranken_analyse()
print(f"Untere Schranke: S(6) ≥ {analyse['untere_schranke']}")   # ≥ 536
print(f"Obere Schranke:  S(6) ≤ {analyse['obere_schranke']}")    # ≤ 1836
print(analyse['status_2026'])
```

### Monte-Carlo-Suche

```python
mc = sz.monte_carlo_partition(200, 6, versuche=500, seed=42)
print(f"Erfolgsrate: {mc['erfolgsrate']:.1%}")
print(f"Partition für {{1,...,200}} in 6 Klassen: {mc['partition_gefunden']}")
```

---

## SAT-Kodierung

Das Problem "Existiert eine sum-freie $k$-Partition von $\{1,\ldots,n\}$?" lässt sich als SAT-Problem kodieren:

**Variablen:** $x_{i,c} \in \{0,1\}$ für $i \in \{1,\ldots,n\}$, $c \in \{0,\ldots,k-1\}$

Bedeutung: $x_{i,c} = 1$ gdw. Zahl $i$ hat Farbe $c$.

**Klausel-Typen:**

| Typ | Anzahl | Länge | Formel |
|-----|--------|-------|--------|
| Mindestens eine Farbe | $n$ | $k$ | $\bigvee_{c=0}^{k-1} x_{i,c}$ für jedes $i$ |
| Höchstens eine Farbe | $n \cdot \binom{k}{2}$ | 2 | $\neg x_{i,c} \vee \neg x_{i,c'}$ für $c \neq c'$ |
| Schur-Bedingung | $\#\text{Tripel} \cdot k$ | 3 | $\neg x_{a,c} \vee \neg x_{b,c} \vee \neg x_{a+b,c}$ |

**Symmetriebrechung:** $x_{1,0} = 1$ (Zahl 1 hat immer Farbe 0).

```python
sat = sz.sat_encoding(536, 6)
print(f"Variablen:  {sat['variablen']['anzahl']}")        # 536*6 = 3216
print(f"Klauseln:   {sat['klauseln']['gesamt']}")
print(f"Schur-Tripel: {sat['klauseln']['typ3_schur_bedingung']['tripel_count']}")
```

---

## Bekannte Schur-Zahlen und Referenzen

| $k$ | $S(k)$ | Beweis / Quelle |
|-----|--------|-----------------|
| 1 | 1 | trivial |
| 2 | 4 | Schur 1916 |
| 3 | 13 | Baumert 1965 |
| 4 | 44 | Baumert 1965 |
| 5 | 160 | **Heule 2017** (SAT-Solver, 2PB-Zertifikat) |
| 6 | ? | $[536, 1836]$ — offen |

**Quellen:**
- L.D. Baumert (1965): *A note on Schur's problem*
- S. Eliahou, J.M. Marín, M.P. Revuelta, M.I. Sanz (2011)
- M. Heule, O. Kullmann, V. Manthey (2012): *Cube and Conquer*
- **M. Heule (2017):** *Schur Number Five* — SAT-Beweis mit gigantischem Zertifikat (2 Petabyte)

---

## Komplexität

- **Entscheidungsproblem** "Hat $\{1,\ldots,n\}$ eine sum-freie $k$-Partition?": NP-vollständig (SAT-Reduktion)
- **Greedy:** $O(n \cdot k \cdot |K_{\max}|)$ pro Zuweisung, praktisch sehr schnell
- **Backtracking:** Im schlimmsten Fall $O(k^n)$, für $n \leq 50$ und $k \leq 4$ nutzbar
- **Monte-Carlo:** $O(\text{versuche} \cdot n \cdot k \cdot |K_{\max}|)$, gut parallelisierbar

---

*Dokumentation generiert für das specialist-maths Projekt. Autor: Michael Fuhrmann.*
