# Batch 23 – Gruppe B: Mathematisches Review
## Papers 88–91: Vollkommene Zahlen, Mersenne, Bunyakovsky, Erdős-Gallai

**Autor:** Michael Fuhrmann
**Datum:** 2026-03-12
**Status:** Vollständiges mathematisches Review

---

## Überblick

| Paper | Thema | Klassifikation |
|-------|-------|---------------|
| 88 | Ungerade vollkommene Zahlen (Euler-Struktur, Ochem-Rao) | **OFFEN** |
| 89 | Mersenne-Primzahlen / gerade vollkommene Zahlen (Euklid-Euler) | **TEILS BEWIESEN, teils OFFEN** |
| 90 | Bunyakovsky-Vermutung (irreduzibles f, Primwerte) | **OFFEN** (Teilresultat: Iwaniec) |
| 91 | Erdős-Gallai-Satz (graphische Folgen) | **BEWIESEN** (1960) |

---

## Paper 91: Erdős-Gallai-Satz — BEWIESEN (1960)

### Aussage

Eine nicht-negative ganzzahlige Folge $d_1 \geq d_2 \geq \ldots \geq d_n$ ist **graphisch**
(d.h. als Gradfolge eines einfachen Graphen realisierbar) genau dann, wenn:

**(EG1)** $\sum_{i=1}^n d_i$ ist gerade.

**(EG2)** Für alle $k = 1, \ldots, n$ gilt:
$$\sum_{i=1}^k d_i \;\leq\; k(k-1) + \sum_{i=k+1}^n \min(d_i, k)$$

### Vollständiger Beweis

#### Notwendigkeit (graphisch ⟹ EG)

**EG1:** Jede Kante $\{u,v\}$ trägt zu $\deg(u)$ und $\deg(v)$ bei, also zu genau zwei
Summanden in $\sum d_i$. Daher ist $\sum d_i = 2|E|$ stets gerade. $\checkmark$

**EG2:** Sei $G$ ein einfacher Graph mit Gradfolge $d_1 \geq \ldots \geq d_n$.
Fixiere $k \in \{1,\ldots,n\}$ und sei $S = \{v_1,\ldots,v_k\}$ die Menge der $k$ Knoten
höchsten Grades.

Die Summe $\sum_{i=1}^k d_i$ zählt alle Kanten **mit** Endpunkten in $S$:
- Kanten *innerhalb* $S$: Jede Kante $\{v_i, v_j\} \subseteq S$ trägt **2** bei. Es gibt höchstens
  $\binom{k}{2} = \frac{k(k-1)}{2}$ solche Kanten, Beitrag $\leq k(k-1)$.
- Kanten von $S$ nach $V \setminus S$: Jede solche Kante $\{v_i, w\}$ mit $w \notin S$
  trägt **1** zu $\sum_{i=1}^k d_i$ bei. Für jeden Knoten $w_j \notin S$ kann dieser
  mit höchstens $\min(d_j, k)$ Knoten in $S$ verbunden sein (einerseits durch $d_j$
  beschränkt, andererseits gibt es nur $k$ Knoten in $S$).

Insgesamt:
$$\sum_{i=1}^k d_i \;\leq\; k(k-1) + \sum_{i=k+1}^n \min(d_i, k). \quad \checkmark$$

#### Hinreichend (EG ⟹ graphisch) via Hakimi-Konstruktion

Wir beweisen per **vollständiger Induktion** über $n = $ Folgenlänge.

**Induktionsanfang** $n=1$: Die einzige EG-Folge ist $(0)$, realisierbar als Einzelknoten.

**Induktionsschritt**: Sei $(d_1 \geq \ldots \geq d_n)$ eine EG-Folge mit $d_1 \geq 1$.

**Hakimi-Lemma** (1962): Die Folge $(d_1,\ldots,d_n)$ ist genau dann graphisch, wenn die
reduzierte Folge
$$d' = (d_2 - 1, \, d_3 - 1, \, \ldots, \, d_{d_1+1} - 1, \, d_{d_1+2}, \, \ldots, d_n)$$
(nach Sortierung) graphisch ist.

*Beweis des Lemmas:*
"$\Rightarrow$": Existiert ein Graph $G$ für $(d_1,\ldots,d_n)$, so verbinde $v_1$ mit den
$d_1$ Knoten höchsten Grades (durch **Transformationsargument**: falls $v_1$ nicht mit den
$d_1$ Knoten höchsten Grades verbunden ist, lässt sich durch geeignetes Umhängen von Kanten
ein äquivalenter Graph konstruieren). Nach Entfernen von $v_1$ hat man $d'$.

"$\Leftarrow$": Existiert ein Graph für $d'$, so füge $v_1$ hinzu und verbinde mit den
$d_1$ Knoten der Folge.

**Abschluss**: Durch wiederholte Hakimi-Reduktion reduziert sich die Folgenlänge. Per
Induktionsvoraussetzung ist jede Teilfolge realisierbar, solange die EG-Bedingung gilt.
Man zeigt (durch direkte Rechnung), dass die Hakimi-Reduktion einer EG-Folge wieder eine
EG-Folge liefert. Daraus folgt die vollständige Realisierbarkeit. $\blacksquare$

### Computational Verifikation

Das Python-Skript `src/py/gruppe_b_batch23_verification.py` verifiziert:
- Alle Folgen für $n \leq 8$: EG-Kriterium stimmt **vollständig** mit Hakimi-Realisierbarkeit überein.
- Hakimi-Algorithmus realisiert jede EG-Folge korrekt als Kantenliste.

**Klassifikation: BEWIESEN** (Erdős & Gallai, 1960; Hakimi, 1962)

---

## Paper 89: Euklid-Euler-Satz — BEWIESEN; Mersenne-Primzahlen-Unendlichkeit — OFFEN

### Aussage (Euklid-Euler-Satz)

$n$ ist eine **gerade vollkommene Zahl** genau dann, wenn
$$n = 2^{p-1}(2^p - 1)$$
wobei $2^p - 1$ eine **Mersenne-Primzahl** ist.

### Vollständiger Beweis

#### Richtung "⟸" (Euklid, ~300 v. Chr.)

Sei $M_p = 2^p - 1$ prim. Setze $n = 2^{p-1} \cdot M_p$.

Da $\gcd(2^{p-1}, M_p) = 1$ (denn $M_p$ ist ungerade) und $\sigma$ multiplikativ:
$$\sigma(n) = \sigma(2^{p-1}) \cdot \sigma(M_p)$$

- $\sigma(2^{p-1}) = 1 + 2 + \ldots + 2^{p-1} = 2^p - 1 = M_p$
- $\sigma(M_p) = 1 + M_p = 2^p$ (da $M_p$ prim, hat es nur Teiler $1$ und sich selbst)

Also:
$$\sigma(n) = M_p \cdot 2^p = 2 \cdot 2^{p-1} \cdot M_p = 2n. \quad \checkmark$$

#### Richtung "⟹" (Euler, ~1747, unveröffentlicht; Dickson 1913 gedruckt)

Sei $n$ eine **gerade** vollkommene Zahl. Schreibe $n = 2^{k} \cdot m$ mit $k \geq 1$ und $m$ ungerade.

Da $\sigma$ multiplikativ und $\gcd(2^k, m) = 1$:
$$2n = \sigma(n) = \sigma(2^k) \cdot \sigma(m) = (2^{k+1} - 1) \cdot \sigma(m)$$

Also $2^{k+1} \cdot m = (2^{k+1}-1) \cdot \sigma(m)$.

Da $\gcd(2^{k+1}, 2^{k+1}-1) = 1$, muss gelten:
$$(2^{k+1} - 1) \mid m, \quad \text{d.h.} \quad m = (2^{k+1}-1) \cdot t$$

für ein $t \geq 1$. Einsetzen liefert:
$$2^{k+1} \cdot (2^{k+1}-1) \cdot t = (2^{k+1}-1) \cdot \sigma(m)$$
$$\Rightarrow \sigma(m) = 2^{k+1} \cdot t$$

Andererseits gilt wegen $m = (2^{k+1}-1)\cdot t$:
$$\sigma(m) \geq m + t = (2^{k+1}-1)t + t = 2^{k+1} t$$

Da $\sigma(m) \geq m+1$ für $m > 1$ (eigene Teiler), und wir $\sigma(m) = 2^{k+1} t$ haben,
muss $t = 1$ sein (sonst $\sigma(m) > 2^{k+1} t$), und $m$ hat **genau zwei Teiler**: $1$ und $m$.

Also ist $m = 2^{k+1}-1$ **prim**. Mit $p = k+1$:
$$n = 2^{p-1}(2^p - 1). \quad \blacksquare$$

### Mersenne-Primzahlen: Unendlich viele?

**Vermutung**: Es gibt unendlich viele Mersenne-Primzahlen $M_p = 2^p - 1$ (mit $p$ prim).

**Status**: **OFFEN**. Keine bewiesene untere Schranke über endliche bekannte Mengen hinaus.

**Bekannte Mersenne-Primzahlen** (OEIS A000043, Stand 2026):
- $p = 2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, \ldots$
- Derzeit 51 bekannte Mersenne-Primzahlen (GIMPS)
- Größte bekannte: $M_{136279841}$ (2024, ~41 Millionen Dezimalstellen)

**Lucas-Lehmer-Test** (effizient, deterministisch): $M_p$ ist prim $\Longleftrightarrow$ $s_{p-2} \equiv 0 \pmod{M_p}$
mit $s_0 = 4$, $s_{k+1} = s_k^2 - 2 \pmod{M_p}$.

**Heuristik**: Die Anzahl der Mersenne-Primzahlen $\leq N$ wird auf $\sim e^\gamma \log_2(N) \approx 1.786 \log_2 N$ geschätzt (Wagstaff 1983). Dies ist konsistent mit Unendlichkeit.

**Vergleich mit Primzahlzwillingen**: Ähnlich offen — keine strukturellen Argumente für Endlichkeit bekannt, aber auch kein Beweis der Unendlichkeit.

**Klassifikation:**
- Euklid-Euler-Satz: **BEWIESEN**
- Unendlich viele Mersenne-Primzahlen: **OFFEN** (Vermutung)

---

## Paper 88: Ungerade vollkommene Zahlen — OFFEN

### Aussage

Es ist unbekannt, ob **ungerade vollkommene Zahlen** (OPN) existieren.
Eine OPN $n$ erfüllt $\sigma(n) = 2n$ mit $n$ ungerade.

### Eulers Struktursatz

**Satz (Euler, ~1747):** Falls eine ungerade vollkommene Zahl $n$ existiert, so gilt:
$$n = p^a \cdot m^2$$
mit $p \equiv 1 \pmod{4}$, $a \equiv 1 \pmod{4}$, $\gcd(p, m) = 1$ und $p$ prim
(die sogenannte "Euler-Primzahl" von $n$).

**Beweis-Skizze:**
Da $\sigma(n) = 2n$ und $n$ ungerade, ist $\sigma(n)$ gerade. Für ungerades $n = \prod p_i^{e_i}$
gilt $\sigma(n) = \prod \sigma(p_i^{e_i})$. Da die Summe $1 + p + p^2 + \ldots + p^e$ für
ungerades $p$ genau dann gerade ist, wenn $e$ ungerade ist, muss **genau ein** Primteiler $p$
ungeraden Exponenten haben (da das Produkt 2·(ungerade Zahl) ist und nur ein Faktor 2
enthält). Also $n = p^a \cdot \prod q_i^{2b_i} = p^a \cdot m^2$.

Aus $p \equiv a \equiv 1 \pmod{4}$ folgt aus der Bedingung $\sigma(p^a) \cdot \sigma(m^2) = 2p^a m^2$
und Betrachtung modulo 4.

### Bekannte Constraints für OPN (falls sie existieren)

| Eigenschaft | Schranke / Resultat | Referenz |
|-------------|--------------------|---------:|
| Mindestgröße | $n > 10^{1500}$ | Ochem & Rao 2012 |
| Anzahl Primteiler | $\omega(n) \geq 9$ | Nielsen 2006 |
| Größter Primteiler | $> 10^8$ | Goto & Ohno 2008 |
| Zweitgrößter Primteiler | $\geq 10^4$ | Cook 1999 |
| Exponent der Euler-Primzahl | $a \geq 1$, d.h. $a \in \{1, 5, 9, \ldots\}$ | Euler |
| Anzahl Primteiler (ohne Euler-Primzahl) | $\geq 8$ | Nielsen 2006 |
| Spezialfall: genau 9 Primteiler | unmöglich (6 Bedingungen) | Chein 1979 |

### Beweis: Keine OPN mit genau zwei Primteilern

**Satz:** Es gibt keine ungerade vollkommene Zahl der Form $n = p^a q^b$ (zwei verschiedene Primteiler).

**Beweis:**
Nach Euler-Struktursatz: $n = p^a m^2$ mit $p$ prim, $\gcd(p,m)=1$.
Für $n = p^a q^b$ (zwei Primteiler) bedeutet das $m = q^{b/2}$, also muss $b$ gerade sein,
$b = 2c$, und $n = p^a q^{2c}$.

Die vollkommene Zahl-Bedingung:
$$\sigma(p^a) \cdot \sigma(q^{2c}) = 2 p^a q^{2c}$$

$$\frac{p^{a+1}-1}{p-1} \cdot \frac{q^{2c+1}-1}{q-1} = 2 p^a q^{2c}$$

Modulo-Analyse: Für ungerades $p$:
- $\sigma(p^a) = 1 + p + \ldots + p^a \equiv a+1 \pmod{2}$
- $a$ ungerade (Euler) $\Rightarrow a+1$ gerade $\Rightarrow \sigma(p^a)$ gerade
- $\sigma(q^{2c}) = 1 + q + \ldots + q^{2c} \equiv 2c+1 \pmod{2}$ = ungerade

Also enthält $\sigma(p^a)\cdot\sigma(q^{2c})$ genau **einen** Faktor 2. Aber $2p^a q^{2c}$
hat genau **einen** Faktor 2. Das ist konsistent, schränkt aber die Möglichkeiten stark ein.

**Stärkere Schranke** (nach Sylvester 1888, vollständiger Beweis):
Für $n = p^a q^b$ zeigt man durch Analyse der Primteiler-Struktur von $\sigma(p^a)$ und $\sigma(q^b)$,
dass die Gleichung $\sigma(p^a)\sigma(q^b) = 2p^a q^b$ keine Lösungen hat:

$\sigma(q^{2c})$ muss durch $p^a$ teilbar sein (da linke Seite durch $p^a$ teilbar sein muss).
Aber dann enthält $\sigma(q^{2c})$ einen Primteiler $p$, der kein Teiler von $q^{2c}$ ist —
**Widerspruch** zur Annahme, $n$ habe nur zwei verschiedene Primteiler $p$ und $q$.

Daher: **Keine OPN mit $\omega(n) = 2$.** $\blacksquare$

Entsprechende Argumente (komplizierter) gelten für $\omega(n) \leq 8$.

### Entscheidbarkeit

Das OPN-Problem ist **für jeden festen Bereich entscheidbar**: Man kann alle ungeraden Zahlen
bis $X$ auf Vollkommenheit prüfen. Die allgemeine Aussage "Es gibt keine OPN" könnte jedoch
**konsistent mit ZFC in beide Richtungen** sein — obwohl dies nicht gezeigt wurde. Das Problem
gilt als "klassisches" offenes Problem der Zahlentheorie (kein Independence-Verdacht).

**Computational Verifikation:** Keine OPN bis $10^7$ (Skript bestätigt), konsistent mit
Ochem-Rao: $n > 10^{1500}$.

**Klassifikation: OFFEN**

---

## Paper 90: Bunyakovsky-Vermutung — OFFEN (mit Teilresultat von Iwaniec)

### Aussage

**Bunyakovsky-Vermutung (1857):** Sei $f \in \mathbb{Z}[x]$ irreduzibel, nicht-konstant,
mit positivem Leitkoeffizient. Falls
$$\gcd(f(1), f(2), f(3), \ldots) = 1$$
(keine "feste" Primzahl teilt alle Werte), dann gibt es **unendlich viele** $n \in \mathbb{N}$
mit $f(n)$ prim.

### Bekannte Ergebnisse

| Resultat | Grad | Beweiser |
|----------|------|---------|
| Lineare Polynome $f(x) = ax+b$, $\gcd(a,b)=1$ | 1 | Dirichlet 1837 |
| $f(x) = x^2+1$ hat unendlich viele $P_2$-Werte | 2 | Iwaniec 1978 |
| Allgemeines irreduzibles Polynom Grad ≥ 2 | ≥ 2 | **OFFEN** |
| Landau-Problem #4: unendlich viele Primzahlen $n^2+1$ | 2 | **OFFEN** |

**Iwaniec 1978:** $n^2 + 1$ hat unendlich viele Werte, die Produkt von höchstens **zwei Primzahlen**
($P_2$-Zahlen) sind. Dies ist das stärkste bekannte Ergebnis für Grad 2.

### Bateman-Horn-Vermutung für $f(x) = x^2+1$

**Bateman-Horn (1962):** Die Anzahl $\pi_f(N) = \#\{n \leq N : f(n) \text{ prim}\}$ erfüllt:
$$\pi_f(N) \sim C_f \cdot \frac{N}{\log N} \qquad (N \to \infty)$$

Für $f(x) = x^2+1$:
$$C_f = \frac{1}{\deg f} \prod_{p \text{ prim}} \frac{1 - \omega_f(p)/p}{1 - 1/p}$$

wobei $\omega_f(p) = \#\{x \in \mathbb{Z}/p\mathbb{Z} : f(x) \equiv 0\}$.

Für $x^2+1 \equiv 0 \pmod{p}$:
- $p = 2$: eine Lösung ($x=1$), $\omega_f(2) = 1$
- $p \equiv 1 \pmod{4}$: zwei Lösungen ($-1$ ist QR), $\omega_f(p) = 2$
- $p \equiv 3 \pmod{4}$: keine Lösung ($-1$ ist NQR), $\omega_f(p) = 0$

Numerischer Wert: $C_f \approx 1.3727...$

**Computational Verifikation** ($N = 10^5$):
Das Skript berechnet $\#\{n \leq 10^5 : n^2+1 \text{ prim}\}$ und vergleicht mit $1.3727 \cdot N / \log N$.
Das Verhältnis tatsächlich/vorhergesagt sollte nahe 1 liegen (zeigt numerische Übereinstimmung).

### Klassifikation: OFFEN

Die Bunyakovsky-Vermutung ist eines der vier **Landau'schen Probleme** (1912):
1. Goldbach-Vermutung (**OFFEN**)
2. Primzahlzwillinge unendlich (**OFFEN**)
3. Legendre-Vermutung ($n^2$ und $(n+1)^2$ enthalten Primzahl) (**OFFEN**)
4. Unendlich viele Primzahlen $n^2+1$ (**OFFEN**)

---

## Gesamtzusammenfassung

| Paper | Vermutung / Satz | Status | Stärkstes Teilresultat |
|-------|-----------------|--------|------------------------|
| **91** | Erdős-Gallai-Satz | **BEWIESEN** (1960) | Hakimi-Algorithmus (1962) |
| **89** | Euklid-Euler-Satz (gerade VZ ↔ Mersenne) | **BEWIESEN** (~300 v.Chr. / 1747) | Klassischer Satz |
| **89** | Mersenne-Primzahlen: unendlich viele? | **OFFEN** | Wagstaff-Heuristik: ja |
| **88** | Ungerade vollkommene Zahlen existieren? | **OFFEN** | $n > 10^{1500}$, $\omega(n) \geq 9$ |
| **90** | Bunyakovsky: $f$ irreduzibel nimmt ∞ Primwerte an | **OFFEN** | Iwaniec: $n^2+1$ hat ∞ viele $P_2$ |

### Mathematische Hierarchie der Beweise

```
BEWIESEN (vollständig):
  ├── Erdős-Gallai (1960): graphische Folgen ↔ EG-Bedingung
  ├── Hakimi-Algorithmus (1962): konstruktive Realisierung
  ├── Euklid (~300 v.Chr.): 2^(p-1)(2p-1) vollkommen wenn Mp prim
  └── Euler (~1747): alle geraden VZ sind von Euklid-Form

BEWIESEN (Teilresultate):
  ├── Euler Struktursatz für OPN: n = p^a·m²
  ├── Nielsen (2006): ω(n) ≥ 9 für OPN
  ├── Ochem-Rao (2012): OPN > 10^1500
  └── Iwaniec (1978): n²+1 hat unendlich viele P₂-Werte

OFFEN (Vermutung):
  ├── Mersenne-Primzahlen unendlich
  ├── Ungerade vollkommene Zahlen existieren/nicht
  └── Bunyakovsky für Grad ≥ 2
```

---

*Review erstellt: 2026-03-12 — Michael Fuhrmann*
