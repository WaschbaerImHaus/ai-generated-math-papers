# brocard_extension.py — Erweiterte Analyse der Brocard-Ramanujan-Vermutung

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Status**: Vermutung offen; numerisch bestätigt bis $n = 10^9$ (Berndt-Galway 2000)

---

## Mathematischer Hintergrund

### Die Brocard-Ramanujan-Vermutung

**Vermutung** (Brocard 1876, Ramanujan 1913):

$$n! + 1 = m^2$$

hat nur die drei Lösungen:

$$(n, m) \in \{(4, 5), (5, 11), (7, 71)\}$$

**Verifikation**:
- $4! + 1 = 24 + 1 = 25 = 5^2$ ✓
- $5! + 1 = 120 + 1 = 121 = 11^2$ ✓
- $7! + 1 = 5040 + 1 = 5041 = 71^2$ ✓

### Bekannte Fakten

| Fakt | Quelle |
|------|--------|
| Verifiziert bis $n = 10^9$ | Berndt & Galway (2000) |
| Vollständige Restklassenabdeckung mod 840 | Bekannt |
| Für $n \geq p$ (Primzahl): $n! \equiv 0 \pmod{p}$, also $n! + 1 \equiv 1 \pmod{p}$ | Elementar |

### Methode: Modulare Ausschlussanalyse

**Idee**: $m^2 \bmod p$ kann nur Werte aus der Menge der **quadratischen Reste** mod $p$ annehmen:

$$QR(p) = \{k^2 \bmod p \mid k = 0, 1, \ldots, p-1\}$$

Für eine Primzahl $p$ gilt $|QR(p)| = \frac{p+1}{2}$ (inkl. 0).

**Ausschluss**: Falls $(n! + 1) \bmod p \notin QR(p)$, dann ist $n! + 1$ kein Quadrat.

**Einschränkung**: Für $n \geq p$ gilt $p \mid n!$, also $n! + 1 \equiv 1 \pmod{p}$.
Da $1 \in QR(p)$ immer, liefern Primzahlen $p \leq n$ **keinen** Ausschluss.
→ Nützliche Moduli sind Primzahlen $p > n$ oder Primzahlquadrate $p^2$ mit geeignetem $p$.

### Theoretische Schranken

Wachstum von $m$: Falls $n! + 1 = m^2$, dann gilt (Stirling-Näherung):

$$m = \sqrt{n! + 1} \approx \sqrt{n!} \approx \left(\frac{n}{e}\right)^{n/2} \cdot (2\pi n)^{1/4}$$

$$\log m \approx \frac{n}{2} \log n - \frac{n}{2}$$

Für $n = 100$ hat $m$ ca. 80 Stellen; für $n = 1000$ ca. 1285 Stellen.

### p-adische Bewertungsanalyse

Die **p-adische Bewertung** $v_p(k)$ ist der größte Exponent $j$ mit $p^j \mid k$.

Für $n! + 1 = m^2$ muss gelten: $v_p(n! + 1)$ ist gerade für alle Primzahlen $p$.

**Legendres Formel** für $v_p(n!)$:

$$v_p(n!) = \sum_{k=1}^{\infty} \left\lfloor \frac{n}{p^k} \right\rfloor$$

Für $n \geq p$: $p \mid n!$, also $n! + 1 \equiv 1 \pmod{p}$, damit $v_p(n!+1) = 0$ (gerade) → kein Ausschluss durch $p \leq n$.

---

## Klassen- und Methodenübersicht

### Hilfsfunktionen (modullevel)

| Funktion | Beschreibung |
|----------|--------------|
| `_quadratische_reste(m)` | Berechnet $QR(m) = \{k^2 \bmod m \mid k \in \mathbb{Z}_m\}$ |
| `_ist_quadrat(n)` | Exakte Quadratprüfung via `math.isqrt` (ganzzahlig) |
| `_fakultaet_mod(n, m)` | Berechnet $n! \bmod m$ effizient (Abbruch bei 0) |

### Klasse `BrocardExtension`

#### Klassenkonstante

```python
BEKANNTE_LOESUNGEN: Dict[int, int] = {4: 5, 5: 11, 7: 71}
```

#### Modulare Ausschlussanalyse

| Methode | Beschreibung |
|---------|--------------|
| `modular_ausschluss(n, moduli_list)` | Schließt $n$ durch modulare Argumente aus |
| `analysiere_restklassen(max_mod)` | Findet nützliche Moduli für $n \in [8, 50]$ |
| `finde_ausschlusskandidaten(n)` | Sucht möglichst viele Ausschluss-Moduli für ein $n$ |
| `vollstaendige_modular_analyse(n_min, n_max)` | Vollständige modulare Analyse über einen Bereich |

#### Numerische Suche

| Methode | Beschreibung |
|---------|--------------|
| `numerische_suche(n_max, verbose)` | Sucht $n! + 1 = m^2$ für alle $n \leq n_{\max}$ |

#### Theoretische Schranken

| Methode | Beschreibung |
|---------|--------------|
| `schranken_analyse()` | Stirling-Schranken und numerische Beispiele |
| `p_adische_analyse(p)` | p-adische Bewertung von $n! + 1$ für Primzahl $p$ |

---

## Wichtigste Algorithmen

### 1. Modulare Ausschlussanalyse (`modular_ausschluss`)

```
Eingabe: n (Kandidat), moduli_list (Liste von Moduli)
Für jeden Modulus m:
  1. r ← n! + 1 (mod m)  [exakt via _fakultaet_mod]
  2. QR ← quadratische Reste mod m
  3. Falls r ∉ QR: → Ausschluss! n! + 1 ist kein Quadrat.
```

**Wichtig**: Ein einziger Ausschluss genügt. Kein Ausschluss ist kein Beweis.

### 2. Effiziente Fakultät mod m (`_fakultaet_mod`)

```python
ergebnis = 1
for k in range(1, n + 1):
    ergebnis = (ergebnis * k) % m
    if ergebnis == 0:
        return 0  # Alle weiteren Faktoren sind nutzlos
```

Für $m$ zusammengesetzt und $m \leq n$: $m \mid n!$, daher früher Abbruch.

### 3. Ausschlusskandidaten für ein bestimmtes $n$ (`finde_ausschlusskandidaten`)

Strategie: Primzahlen $p > n$ sind besonders geeignet (da $p \nmid n!$).
Zusätzlich: Primzahlquadrate $p^2$ mit $p \leq n$ aber $p^2 > n$.

```
Kandidaten: [2..49] ∪ {Primzahlen in [n+1, n+200]} ∪ {p² | p prim, p ≤ n, p² > n}
```

### 4. Numerische Suche (`numerische_suche`)

```
fak = 1  (0! = 1, dann inkrementell)
Für n = 1, 2, ..., n_max:
  fak ← fak * n        (inkrementell: kein erneutes math.factorial nötig)
  val ← fak + 1
  r   ← isqrt(val)
  Falls r² == val: Lösung gefunden!
```

Inkrementelle Berechnung ist entscheidend für Performance bei großen $n$.

### 5. p-adische Analyse (`p_adische_analyse`)

Für jedes $n \in [1, 100]$:
- Berechne $v_p(n!)$ via Legendres Formel: $\sum_{k \geq 1} \lfloor n / p^k \rfloor$
- Berechne $v_p(n! + 1)$ durch direkte Division
- **Ausschluss**: Falls $v_p(n! + 1)$ ungerade → $n! + 1$ ist kein Quadrat

---

## Beispielanwendungen

```python
from brocard_extension import BrocardExtension
import time

ext = BrocardExtension()

# Numerische Suche bis n = 1000
res = ext.numerische_suche(1000, verbose=False)
print(f"Lösungen: {[(s['n'], s['m']) for s in res['loesungen']]}")
print(f"Status: {res['vermutung_status']}")
# Lösungen: [(4, 5), (5, 11), (7, 71)]
# Status: BESTÄTIGT (im geprüften Bereich)

# Modulare Ausschlussanalyse für n=10
kandidaten = ext.finde_ausschlusskandidaten(10)
print(f"n=10 ausgeschlossen: {kandidaten['ausgeschlossen']}")
stärkstes = kandidaten['stärkstes_argument']
if stärkstes:
    print(f"Durch Modulus {stärkstes['modulus']}: "
          f"10!+1 ≡ {stärkstes['rest']} (mod {stärkstes['modulus']})")

# p-adische Analyse für p=5
pad = ext.p_adische_analyse(5)
print(pad['erklaerung'])
# Für n ≥ 5: 5 | n!, also n!+1 ≡ 1 (mod 5) → v_5(n!+1) = 0 (gerade). Kein Ausschluss.
ausschluesse = [e['n'] for e in pad['ausschluss_faelle']]
print(f"Ausschlussfälle (p=5): {ausschluesse[:10]}")

# Vollständige modulare Analyse n ∈ [8, 50]
voll = ext.vollstaendige_modular_analyse(8, 50)
print(f"Ausgeschlossen: {len(voll['ausgeschlossen'])}/{50-8+1} "
      f"({voll['anteil_ausgeschlossen']}%)")
print(f"Nicht ausgeschlossen: {voll['nicht_ausgeschlossen']}")
```

### Konkretes Beispiel: Modularer Ausschluss von $n = 10$

$$10! = 3628800, \quad 10! + 1 = 3628801$$

Für Modulus $m = 121 = 11^2$:
- $10! \bmod 121$: Da $11 \leq 10$? Nein ($11 > 10$), also muss $10! \bmod 121$ direkt berechnet werden.
- $QR(121)$: Quadratische Reste mod 121.
- Falls $10! + 1 \bmod 121 \notin QR(121)$: Ausschluss!

### Schranken-Beispiele

| $n$ | Stellen von $n!$ | $m = \lfloor\sqrt{n!+1}\rfloor$ | Stellen von $m$ | Lösung? |
|-----|-----------------|--------------------------------|-----------------|---------|
| 4   | 2               | 5                              | 1               | **ja** |
| 5   | 3               | 11                             | 2               | **ja** |
| 7   | 4               | 71                             | 2               | **ja** |
| 10  | 7               | 1904                           | 4               | nein |
| 20  | 19              | $\approx 1.5 \times 10^9$      | 10              | nein |
| 50  | 65              | $\approx 3.0 \times 10^{32}$   | 33              | nein |

---

## Wichtige Hinweise

- **Nützliche Moduli**: Primzahlen $p > n$ (da $p \nmid n!$) und Primzahlquadrate sind effektiver als Primzahlen $\leq n$.
- **Bekannte Lösungen**: Das Modul prüft automatisch, dass die 3 bekannten Lösungen $(4, 5)$, $(5, 11)$, $(7, 71)$ niemals fälschlich ausgeschlossen werden (`assert` in `p_adische_analyse`).
- **Status**: Die Brocard-Ramanujan-Vermutung ist **offen**. Modulare Argumente können $n$ als Nicht-Lösung ausschließen, aber nicht beweisen, dass $n! + 1$ ein Quadrat ist.
- **Performance**: Die inkrementelle Fakultätsberechnung ermöglicht Suchen bis $n \approx 10000$ in vernünftiger Zeit.
