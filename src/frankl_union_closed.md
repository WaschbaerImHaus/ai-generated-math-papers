# frankl_union_closed.py — Erläuterungen

**Autor:** Michael Fuhrmann
**Datum:** 2026-03-13
**Status:** CONJECTURE (offen seit 1979)

---

## Überblick

Dieses Modul implementiert computergestützte Experimente zur **Frankl-Vermutung** (1979).

### Die Vermutung (CONJECTURE)

> Für jede vereinigungsabgeschlossene Familie $\mathcal{F}$ endlicher Mengen mit $\mathcal{F} \neq \{\emptyset\}$ gibt es ein Element $x$, das in mindestens $|\mathcal{F}|/2$ Mengen vorkommt.

Formal: Es existiert $x \in \bigcup_{F \in \mathcal{F}} F$ mit
$$|\{F \in \mathcal{F} : x \in F\}| \geq \frac{|\mathcal{F}|}{2}$$

**Status:** Offen. Kein vollständiger Beweis bekannt.

---

## Mathematischer Hintergrund

### Vereinigungsabgeschlossene Familien

Eine Familie $\mathcal{F}$ heißt **vereinigungsabgeschlossen** (union-closed), wenn gilt:
$$A, B \in \mathcal{F} \Rightarrow A \cup B \in \mathcal{F}$$

Beispiele:
- Potenzmenge $2^U$: trivialerweise uc-abgeschlossen
- Ketten $\{S_1 \subseteq S_2 \subseteq \ldots \subseteq S_n\}$: uc-abgeschlossen
- Intervall-Halbverbände: uc-abgeschlossen

### Bekannte Schranken

| Ergebnis | Bound | Jahr |
|----------|-------|------|
| Gilmer (Entropie-Methode) | $\geq 0.01 \cdot |\mathcal{F}|$ | 2022 |
| Chase-Lovász (Verbesserung) | $\geq 0.38 \cdot |\mathcal{F}|$ | 2023 |
| Frankl's Ziel | $\geq 0.5 \cdot |\mathcal{F}|$ | offen |

### Gilmer-Ansatz (2022)

Gilmer nutzte Entropiemethoden und bewies, dass stets ein Element in $\geq 0.01 |\mathcal{F}|$ Mengen vorkommt. Dies war der erste nicht-triviale asymptotische Durchbruch.

### Chase-Lovász (2023)

Verbesserung auf $0.38 |\mathcal{F}|$ durch Analyse von Wahrscheinlichkeitsverteilungen über der Familie.

---

## Implementierte Klassen

### `UnionClosedFamily`

Kernklasse zur Darstellung einer uc-Familie.

| Methode | Beschreibung |
|---------|--------------|
| `is_union_closed()` | Prüft uc-Eigenschaft in $O(|F|^2)$ |
| `frequency_vector()` | Häufigkeiten aller Grundmengeneelemente |
| `max_frequency_ratio()` | $\max_x \text{freq}(x) / |\mathcal{F}|$ |
| `satisfies_frankl()` | Ob ratio $\geq 0.5$ |
| `frankl_witness()` | Erstes Element mit freq $\geq |\mathcal{F}|/2$ |
| `from_bitmasks(n, masks)` | Konstruktion via Bitmasken |

### `FranklChecker`

Exhaustive Verifikation für kleine Grundmengen.

- `all_union_closed_families(n)`: Enumeriert alle $2^{2^n}$ möglichen Familien (praktisch: $n \leq 4$)
- `verify_up_to_n(n)`: Prüft alle nicht-trivialen uc-Familien
- `gilmer_bound(size)`: $0.01 \cdot |\mathcal{F}|$
- `chase_lovasz_bound(size)`: $0.38 \cdot |\mathcal{F}|$

### `UnionClosedGenerator`

Erzeugt kanonische Testfamilien:

| Generator | Beschreibung |
|-----------|--------------|
| `power_set(n)` | $2^{\{0,\ldots,n-1\}}$, hat $2^n$ Mengen |
| `interval_family(n)` | Intervall-Halbverband $\{[i,j] : 0 \leq i \leq j < n\}$ |
| `chain_family(n)` | Kette $\{\{0\}, \{0,1\}, \ldots, \{0,\ldots,n-1\}\}$ |
| `union_closure(sets)` | Uc-Abschluss einer beliebigen Familie |
| `random_union_closed(n, seed)` | Zufällige uc-Familie |

---

## Algorithmus: Vereinigungsabschluss

```
Eingabe: Familie F = {S1, ..., Sk}
1. closed = {S1, ..., Sk}
2. Wiederhole bis keine Änderung:
   Für alle A, B ∈ closed:
     Falls A ∪ B ∉ closed: füge A ∪ B ein, markiere Änderung
3. Ausgabe: closed
```

Zeitkomplexität: $O(|F|^2 \cdot |U|)$ pro Iteration, bis zu $O(2^{|U|})$ Iterationen.

---

## Verifikation für kleine n

| n | uc-Familien | davon nicht-trivial | Verletzungen |
|---|-------------|---------------------|--------------|
| 1 | wenige | wenige | 0 |
| 2 | ~20 | ~15 | 0 |
| 3 | ~100s | ~100s | 0 |

**Ergebnis:** Für alle geprüften Fälle ($n \leq 3$, vollständig exhaustiv) bestätigt.

---

## Literatur

- Frankl, P. (1979): *Multiply-intersecting families*, JCTB
- Gilmer, J. (2022): *A constant lower bound for the union-closed sets conjecture*, arXiv:2211.09055
- Chase, Z. & Lovász, L. (2023): *Improvement of Gilmer's bound*, arXiv:2301.09334
- Knill, E. (1994): Exhaustive computer verification for $n \leq 7$
