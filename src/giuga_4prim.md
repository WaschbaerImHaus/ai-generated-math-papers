# giuga_4prim.py — Formaler Beweisversuch: Kein 4-Primfaktor-Giuga-Pseudoprime

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Status**: Numerisch bestätigt (bis $n < 10^6$), algebraisch offen (Stand 2026)

---

## Mathematischer Hintergrund

### Giuga-Zahlen und Giuga-Pseudoprimes

Eine **Giuga-Zahl** ist eine quadratfreie zusammengesetzte Zahl $n = p_1 \cdot p_2 \cdots p_k$ mit paarweise verschiedenen Primfaktoren, die die **schwache Giuga-Bedingung** erfüllt:

$$\forall i: p_i \mid \left(\frac{n}{p_i} - 1\right)$$

Ein **Giuga-Pseudoprime** erfüllt zusätzlich die **starke Giuga-Bedingung**:

$$\forall i: (p_i - 1) \mid \left(\frac{n}{p_i} - 1\right)$$

Die Giuga-Vermutung (1950) besagt, dass es kein Giuga-Pseudoprime gibt — oder äquivalent:

$$n \text{ ist prim} \iff \sum_{k=1}^{n-1} k^{n-1} \equiv -1 \pmod{n}$$

### Bekannte Giuga-Zahlen (nur schwache Bedingung)

Die ersten bekannten Giuga-Zahlen sind: $\{30, 858, 1722, 66198, \ldots\}$

Keines dieser Beispiele erfüllt zusätzlich die starke Bedingung.

### Literaturstand

| Quelle | Ergebnis |
|--------|----------|
| Borwein et al. (1996) | Jedes Giuga-Pseudoprime hat $\geq 59$ Primfaktoren |
| Borwein (erweitert) | $\geq 13635$ Primfaktoren notwendig |
| Bednarek (2014) | Jedes Giuga-Pseudoprime hat $> 10^{19907}$ Dezimalstellen |

Diese Resultate schließen implizit den 4-Prim-Fall aus — aber mittels tiefer analytischer Methoden, nicht elementar.

### Aufbau auf beweisversuche.py

Dieses Modul baut auf den Sätzen 1–4 aus `beweisversuche.py` auf:

- **Satz 1**: $n$ ist quadratfrei
- **Satz 2**: Kein 2-Prim-Giuga-Pseudoprime ($n = p \cdot q$)
- **Satz 3**: Kein 3-Prim-Giuga-Pseudoprime mit $p_1 = 2$
- **Satz 4**: Kein 3-Prim-Giuga-Pseudoprime mit allen ungeraden Faktoren
- **Korollar**: Kein 3-Prim-Giuga-Pseudoprime existiert

Das Ziel ist die Erweiterung auf $k = 4$ Primfaktoren.

---

## Klassen- und Methodenübersicht

### Freie Funktionen

| Funktion | Beschreibung |
|----------|--------------|
| `berechne_giuga_bedingung(n)` | Prüft schwache und starke Giuga-Bedingung für $n$ |
| `numerische_suche_4prim(grenze)` | Sucht alle 4-Prim-Kandidaten bis zur Grenze |
| `schranke_fuer_4prim_fall_p_gleich_2(q_max)` | Schrankenanalyse für $n = 2 \cdot q \cdot r \cdot s$ |
| `schranke_fuer_4prim_fall_p_gleich_3(q_max)` | Schrankenanalyse für $n = 3 \cdot q \cdot r \cdot s$ |
| `schranke_fuer_4prim_allgemeiner_fall(p_max)` | Schrankenanalyse für $n = p \cdot q \cdot r \cdot s$, $p \geq 5$ |

### Klasse `Giuga4PrimBeweis`

| Methode | Beschreibung |
|---------|--------------|
| `fall1_p_gleich_2_analyse()` | Vollständige Analyse für $p_1 = 2$ |
| `fall2_p_gleich_3_analyse()` | Vollständige Analyse für $p_1 = 3$ |
| `fall3_p_ungerade_analyse()` | Allgemeiner Fall $p_1 \geq 5$ |
| `schranken_analyse()` | Zusammenfassung aller Schranken + Literaturvergleich |
| `numerische_verifikation(grenze)` | Numerische Suche bis `grenze` |
| `korollar_kein_4prim_pseudoprime()` | Status-Zusammenfassung des Beweisversuchs |
| `vollstaendige_analyse()` | Führt alle Analysen aus |

---

## Wichtigste Algorithmen

### 1. Schwache und starke Giuga-Bedingung (`berechne_giuga_bedingung`)

Für $n = p_1 \cdot p_2 \cdots p_k$ (quadratfrei, zusammengesetzt):

**Schwach**: $p_i \mid \left(\frac{n}{p_i} - 1\right)$ für alle $i$
**Stark**: $(p_i - 1) \mid \left(\frac{n}{p_i} - 1\right)$ für alle $i$ (für $p_i = 2$ trivial)

Die Funktion gibt ein Dictionary zurück mit `ist_giuga_zahl` (nur schwach) und `ist_giuga_pseudoprime` (beide Bedingungen).

### 2. Effiziente 4-Prim-Suche (`numerische_suche_4prim`)

Statt über alle $n \leq$ Grenze zu iterieren, werden direkt alle Produkte $p_1 \cdot p_2 \cdot p_3 \cdot p_4 \leq$ Grenze mit $p_1 < p_2 < p_3 < p_4$ erzeugt:

- Maximaler $p_1$: $\lfloor \text{Grenze}^{1/4} \rfloor$
- Vorzeitige Abbrüche durch Schranken (z. B. $p_1 \cdot p_2 \cdot (p_2+2) \cdot (p_2+4) > \text{Grenze}$)

Dies ist deutlich effizienter als eine naive Iteration über alle zusammengesetzten Zahlen.

### 3. Schrankenargument für Fall $p_1 = 2$

Sei $n = 2 \cdot q \cdot r \cdot s$ mit $2 < q < r < s$ (alle ungerade prim).

Aus der starken Bedingung für $s$ (Kombination von $s \mid (2qr - 1)$ und $(s-1) \mid (2qr-1)$) folgt wegen $\gcd(s, s-1) = 1$:

$$s \cdot (s-1) \mid (2qr - 1) \implies s^2 - s \leq 2qr - 1 < 2qr \implies s < \sqrt{2qr}$$

Analog ergeben sich Schranken für $q$ und $r$. Das System ist jedoch **konsistent** — kein Widerspruch entsteht.

**Warum der 3-Prim-Beweis versagt**: Im 3-Prim-Fall $2 \cdot q \cdot r$ gilt $r \mid (2q-1)$ mit $r > q$, also muss $r = 2q - 1$ sein (eindeutig!). Dann ist $r - 1 = 2(q-1)$ gerade, aber $2q - 1$ ungerade → Paritätswiderspruch. Im 4-Prim-Fall kann $2qr - 1$ ein Vielfaches $s \cdot k$ mit $k \geq 1$ sein — keine Eindeutigkeit.

### 4. Allgemeines Schrankenargument für $p_1 \geq 5$

Für den größten Faktor $s$ gilt:

$$s(s-1) \mid (pqr - 1) \implies s(s-1) \leq pqr - 1 < pqr$$

Da $s \geq r + 2$ (nächste ungerade Primzahl):

$$(r+2)(r+1) \leq s(s-1) < pqr \implies pq > r + 3 \quad (*)$$

Im Gegensatz zum 3-Prim-Fall $p \cdot q \cdot r$ (wo $p \geq q + 4$ folgt, aber $p < q$ → Widerspruch) ist $(*)$ mit $p < q < r$ vereinbar (z. B. $p = 5, q = 7 \Rightarrow pq = 35 > r + 3$ für $r \leq 31$).

---

## Fallunterscheidung: 3-Prim vs. 4-Prim

| Eigenschaft | 3-Prim-Fall | 4-Prim-Fall |
|-------------|-------------|-------------|
| Schlüsselschranke | $s(s-1) \mid pq - 1$ | $s(s-1) \mid pqr - 1$ |
| Folgerung | $p \geq q + 4$ → Widerspruch | $pq > r + 3$ → konsistent |
| Algebraischer Beweis | Vollständig (beweisversuche.py) | Offen (Stand 2026) |
| Numerisch bis | $10^6$ geprüft | $10^6$ geprüft |

---

## Beispielanwendungen

```python
from giuga_4prim import berechne_giuga_bedingung, Giuga4PrimBeweis

# Bekannte Giuga-Zahlen prüfen (nur schwache Bedingung)
for gz in [30, 858, 1722, 66198]:
    r = berechne_giuga_bedingung(gz)
    print(f"n={gz}: Giuga-Zahl={r['ist_giuga_zahl']}, "
          f"Pseudoprime={r['ist_giuga_pseudoprime']}")
# n=30:    Giuga-Zahl=True,  Pseudoprime=False
# n=858:   Giuga-Zahl=True,  Pseudoprime=False
# n=1722:  Giuga-Zahl=True,  Pseudoprime=False

# Hauptklasse: Vollständige Analyse
beweis = Giuga4PrimBeweis()

# Numerische Suche bis 10^6
ergebnis = beweis.numerische_verifikation(1_000_000)
print(ergebnis['fazit'])
# → "Kein 4-Prim-Giuga-Pseudoprime bis 1000000 gefunden."

# Status des Beweisversuchs
korollar = beweis.korollar_kein_4prim_pseudoprime()
print(korollar['status'])
# → "NUMERISCH BESTÄTIGT (bis 10^6), algebraisch OFFEN"

# Schrankenanalyse
schranken = beweis.schranken_analyse()
print(schranken['beispiel_schranke'])
# → {'p': 3, 'q': 5, 'r': 7, 'pqr_minus_1': 104, 's_kandidaten': [13]}
```

### Schranken-Beispiel: p=3, q=5, r=7

$$s \cdot (s-1) \mid (3 \cdot 5 \cdot 7 - 1) = 104$$

Primteiler von 104 die $> 7$ sind: $\{13\}$.
Kandidat $n = 3 \cdot 5 \cdot 7 \cdot 13 = 1365$.
Prüfung der Giuga-Bedingungen ergibt: **kein Giuga-Pseudoprime**.

---

## Bekannte Lücken und offene Fragen

Das Modul dokumentiert explizit, wo der Beweis noch unvollständig ist:

1. **Kein algebraischer Widerspruch für $k = 4$**: Die Schranken-Methode aus Satz 4 ist nicht direkt übertragbar.
2. **Notwendige nächste Schritte**:
   - Implementierung der Borwein-Methodik (vollständige Primorialschranken) für $k = 4$
   - Untersuchung zusätzlicher Kongruenzbedingungen (mod kleiner Primzahlen)
   - Prüfung ob Induktion von $k = 3$ auf $k = 4$ möglich ist

---

## Wichtige Hinweise

- **Theorem vs. Conjecture**: Die Giuga-Vermutung ist **offen**. Das Modul markiert alle unvollständigen Beweise explizit.
- Bekannte Literaturergebnisse (Borwein, Bednarek) schließen $k = 4$ **implizit** aus, aber ohne elementaren Direktbeweis.
- Das Modul dient als Lernmaterial zur Schrankenanalyse in der elementaren Zahlentheorie.
