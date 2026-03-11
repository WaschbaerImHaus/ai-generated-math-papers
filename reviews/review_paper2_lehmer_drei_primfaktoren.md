# Gutachten: Paper 2
**Titel:** „Lehmers Totient-Problem für Produkte dreier Primzahlen"
**Dateien:** `paper2_lehmer_three_primes.tex` / `paper2_lehmer_drei_primfaktoren_de.tex`
**Gutachter:** Mathematik-Review-Spezialist (Zahlentheorie)
**Datum Erstgutachten:** 2026-03-11
**Datum Revisionsgutachten:** 2026-03-11
**Gesamturteil:** ✅ EN: Annahme ohne Auflagen | ⚠️ DE: Annahme mit einer kleinen Auflage (neuer redaktioneller Überrest)

---

## Revisionsgutachten (2026-03-11)

### BUG-004 (EN) — Status: ✅ BEHOBEN

Die Lücke in der Ungleichungskette für den Semiprim-Fall ($p \ge 7$) wurde vollständig behoben.
Der neue Beweis zeigt:
$$(p-1)(q-1) - (p+q-2) = (p-2)(q-2) - 1 \ge 1 \cdot 2 - 1 = 1 > 0$$
(gültig für alle $p \ge 3$, $q \ge p+2 \ge 5$). ✅ Lückenlos korrekt.

### BUG-002 (DE) — Status: ✅ BEHOBEN

Die fehlerhafte Formelkette „$r = 2q-1+1-1 = 2q$" wurde entfernt.
Der Fall $b=1$ lautet jetzt korrekt: „aus $r-1 = 2q-1$ folgt $r = 2q$." ✅

### BUG-003 (DE) — Status: ✅ BEHOBEN

Lemma `lem:grenze_q` hat jetzt einen vollständigen Beweis via Schlüsselidentität:
$$pqr - 1 = 8abc + 4(ab+bc+ca) + 2(a+b+c)$$
Daraus folgt korrekt $(p-1)(q-1) \mid (p+q-1)+k$, was zu $q \le \frac{3(p-1)}{p-2}$ führt. ✅

**Anmerkung:** ✅ Angepasst (2026-03-11) — Der Verweis auf die Schlüsselidentität
(„Reduktion aus Schritt 2") wurde explizit eingefügt. Lesbarkeit damit verbessert.

### NEU GEFUNDENES PROBLEM (DE): Redaktioneller Überrest in Section 3

**Schwere:** Gering (stilistisch)
**Datei:** `paper2_lehmer_drei_primfaktoren_de.tex`, Zeilen 145–147

Nach dem vollständigen Paritätsargument (Zeile 141–143) folgt:
```latex
Alternativ: Sei $b = \tfrac{2qr-1}{r-1}$. Dann $2q = 1 + \tfrac{b-1}{1}\cdot\ldots$
```
Die Formel bricht mit `\ldots` ab — ein offensichtlicher Entwurfsüberrest.
Da der Paritätsbeweis (Zeile 141–143) bereits vollständig und korrekt ist,
sollte der gesamte „Alternativ:"-Block (Zeilen 144–156) entweder vollständig
ausgeführt oder ersatzlos gestrichen werden.

**Empfehlung EN: Druckreif.**
**Empfehlung DE: Kleine Überarbeitung (Alternativ-Block entfernen oder vervollständigen).**

---

---

## 1. Einordnung

Das Paper beweist, dass kein Produkt genau dreier Primzahlen eine Lehmer-Zahl sein kann
(d.h. $\varphi(n) \mid (n-1)$ impliziert Primalität für solche $n$). Der englische Text
präsentiert einen sauberen, vollständigen Beweis mit einer eleganten „Key Identity"-Technik.
Die deutsche Version enthält darüber hinaus denselben Kern, leidet aber unter
Rechendarstellungsfehlern und einem unvollständigen Beweis-Sketch.

---

## 2. Lemma: Quadratfreiheit

**Behauptung:** Jede Lehmer-Zahl ist quadratfrei.

**Beweis:** Sei $p^2 \mid n$. Dann gilt $p(p-1) \mid \varphi(p^2)$ und damit $p \mid \varphi(n)$.
Aus $\varphi(n) \mid (n-1)$ folgt $p \mid (n-1)$. Da auch $p \mid n$: $p \mid 1$. Widerspruch. ✅

---

## 3. Proposition: Kein Semiprim

**Behauptung:** Es gibt keine Lehmer-Zahl der Form $n = pq$, $p < q$ prim.

**Kernidentität:** $pq - 1 = (p-1)(q-1) + (p-1) + (q-1)$

**Verifikation:**
$(p-1)(q-1) + (p-1) + (q-1) = pq - p - q + 1 + p - 1 + q - 1 = pq - 1$. ✅

Damit: $(p-1)(q-1) \mid (pq-1)$ genau dann wenn $(p-1)(q-1) \mid (p+q-2)$.

**Fall $p = 2$:** $(q-1) \mid q$. Da $\gcd(q-1, q) = 1$ (aufeinanderfolgend): $q-1 \mid 1$,
also $q = 2 = p$. Widerspruch zu $p < q$. ✅

**Fall $p \ge 3$:** Es muss $(p-1)(q-1) \le p+q-2$ gelten, damit Teilbarkeit möglich ist.
Das ist äquivalent zu $pq \le 2p + 2q - 3$.

- Für $p = 3, q \ge 5$: $3q \le 2q + 3 \implies q \le 3$. Widerspruch zu $q \ge 5$. ✅
- Für $p \ge 5, q > p$: $(p-1)(q-1) \ge 2(q-1) > q-1+(q-p) = q+q-p-1 > q+q-2(q-1)$...

**Klarerer Beweis:** $(p-1)(q-1) > p+q-2$ für $p < q$ (beide $\ge 2$), denn:
$$pq - p - q + 1 > p + q - 2 \iff pq > 2p + 2q - 3.$$
Für $p = 3, q \ge 5$: $15 > 13$. ✓ Für $p = 5, q \ge 7$: $35 > 17$. ✓ Allgemein
$pq \ge p(p+2) > 2p + 2(p+2) - 3 = 4p + 1$ für $p \ge 3$ (da $p(p+2) - 4p - 1 = p^2 - 2p - 1 > 0$ für $p \ge 3$).

Die englische Version liefert für $p \ge 5$ die Kette:
$$pq \ge 5q > 2q+7 \ge 2p+2q-3 \text{ (für } p \ge 5\text{)}.$$
Der letzte Schritt $2q+7 \ge 2p+2q-3$ bedeutet $7 \ge 2p-3$, d.h. $p \le 5$.
Damit ist die Kette nur für $p = 5$ gültig — **für $p \ge 7$ bricht die Kette** an der
letzten Ungleichung ab. Die Schlussfolgerung ist trotzdem korrekt (Proposition stimmt),
aber der Beweis hat für $p \ge 7$ eine Lücke in der Darstellung.

⚠️ **Kleinere Darstellungsschwäche** (englische Version): Die Ungleichungskette für
$p \ge 7$ ist so wie präsentiert unvollständig. Empfehlung: Entweder induktiv oder
durch $pq > 2p + 2q - 3 \iff (p-2)(q-2) > 1$ (für $p,q \ge 3$) zeigen.

---

## 4. Theorem: Fall $n = 2qr$ ($3 \le q < r$ prim)

### 4.1 Englische Version

$\varphi(2qr) = (q-1)(r-1)$. Bedingung: $(q-1)(r-1) \mid (2qr-1)$.

Insbesondere $(r-1) \mid (2qr-1)$. Da $r \equiv 1 \pmod{r-1}$: $(r-1) \mid (2q-1)$.
Setze $2q-1 = b(r-1)$, $b \ge 1$:

**Fall $b \ge 2$:** $2q = 1 + b(r-1) \ge 1 + 2(r-1) = 2r-1$, also $q \ge r - \frac{1}{2}$,
d.h. $q \ge r$. Widerspruch zu $q < r$. ✅

**Fall $b = 1$:** $r - 1 = 2q - 1$, also $r = 2q$. Da $q \ge 3$ ungerade Primzahl,
ist $r = 2q$ gerade und $> 2$, also zusammengesetzt. Widerspruch. ✅

**Bewertung Englisch:** ✅ Vollständig korrekt.

### 4.2 Deutsche Version

Zunächst ein eleganterer Einzeiler-Beweis:
> $(q-1)(r-1)$ ist **gerade** (da $q, r > 2$), aber $2qr - 1$ ist **ungerade**.
> Eine gerade Zahl kann keine ungerade teilen. Widerspruch.

✅ **Dieser Beweis ist korrekt und vollständig** — er übertrifft die Länge des
englischen Beweises deutlich.

Danach folgt „Alternativ:" mit einem längeren Beweis, der einen **Darstellungsfehler**
enthält:

> `Fall $b = 1$: $r = 2q - 1 + 1 - 1 = 2q$`

Die Formelkette $2q - 1 + 1 - 1 = 2q - 1 \ne 2q$ ist **arithmetisch falsch**.
Korrekt wäre: aus $r - 1 = 2q - 1$ folgt $r = 2q$.
Das Endergebnis $r = 2q$ stimmt, aber der Lösungsweg ist fehlerhaft notiert.

❌ **Arithmetikfehler in der Formeldarstellung** (deutsches Alternativargument).

---

## 5. Theorem: Fall $n = pqr$ (alle ungerade) — Englische Version

### Key Identity
$\alpha = p + q - 1$, $\beta = \frac{pq-1}{2}$, $a = \frac{p-1}{2}$, $b = \frac{q-1}{2}$, $c = \frac{r-1}{2}$.

**Behauptung:** $2(ab+bc+ca) + (a+b+c) = \alpha c + \beta$.

**Verifikation:**
$$\alpha c = (p+q-1) \cdot \frac{r-1}{2} = (2a+2b+1)c = 2ac + 2bc + c.$$
$$\beta = \frac{(2a+1)(2b+1)-1}{2} = \frac{4ab+2a+2b}{2} = 2ab + a + b.$$
$$\alpha c + \beta = 2ac + 2bc + c + 2ab + a + b = 2(ab+bc+ca) + (a+b+c). \quad ✅$$

**Expansion:** $pqr - 1 = 8abc + 4(ab+bc+ca) + 2(a+b+c)$.

**Verifikation:**
$(2a+1)(2b+1)(2c+1) = (4ab+2a+2b+1)(2c+1) = 8abc + 4ab + 4ac + 4bc + 2a + 2b + 2c + 1$.
Also $pqr - 1 = 8abc + 4(ab+bc+ca) + 2(a+b+c)$. ✅

**Hieraus:** $8abc \mid (pqr-1) \iff 4abc \mid 2(ab+bc+ca) + (a+b+c) = \alpha c + \beta$. ✅

### Schritt 1: Reduktion modulo $c$
Da $c \mid 4abc$: $c \mid (\alpha c + \beta)$ iff $c \mid \beta = \frac{pq-1}{2}$,
d.h. $(r-1) \mid (pq-1)$.
Setze $k = \frac{pq-1}{r-1}$ (ganzzahlig). Dann $\beta = kc$.
$4abc \mid c(\alpha + k)$, also $4ab \mid (\alpha + k)$, d.h. $(p-1)(q-1) \mid (p+q-1+k)$. ✅

### Schritt 2: Schranke für $k$
$r \ge q+2$ (ungerade Primzahlen), also $r-1 \ge q+1 > q$.
$$k = \frac{pq-1}{r-1} < \frac{pq}{q} = p \implies k \le p-1. \quad ✅$$

### Schritt 3: Größenschranke
$(p-1)(q-1) \le (p+q-1) + k \le (p+q-1) + (p-1) = 2p+q-2$.
$$pq - p - q + 1 \le 2p + q - 2 \implies q(p-2) \le 3(p-1) \implies q \le \frac{3(p-1)}{p-2}. \quad ✅$$

### Schritt 4: Fallanalyse
- $p = 3$: $q \le 6$, also $q = 5$ (einzige Primzahl mit $3 < q \le 6$). ✅
- $p = 5$: $q \le 4 < 5 < q$. Widerspruch. ✅
- $p \ge 7$: $\frac{3(p-1)}{p-2} = 3\left(1 + \frac{1}{p-2}\right) \le 3 \cdot \frac{6}{5} = 3{,}6 < 4$.
  Da $q > p \ge 7 > 4$: Widerspruch. ✅

### Schritt 5: Ausschluss $(p,q) = (3,5)$
$(r-1) \mid (pq-1) = 14$. Teiler von 14: $\{1, 2, 7, 14\}$, also $r \in \{2, 3, 8, 15\}$.
- $r = 2, 3$: nicht $> q = 5$.
- $r = 8 = 2^3$: zusammengesetzt.
- $r = 15 = 3 \cdot 5$: zusammengesetzt.
Kein gültiges Prim-$r$. Widerspruch. ✅

**Bewertung Englisch:** ✅ Vollständig korrekt und lückenlos.

---

## 6. Theorem: Fall $n = pqr$ (alle ungerade) — Deutsche Version

Die deutsche Version führt denselben Beweis durch, jedoch:

### Schritt 3 (Schritt 3 im deutschen Text)
Der Text beginnt: „Aus $k \le p-1$ und $r-1 = \frac{pq-1}{k} \ge \frac{pq-1}{p-1}$ folgt..."
Dies ist **nur für $k = p-1$ gültig**. Für $k < p-1$ ist $r-1$ größer (was unproblematisch
für den Widerspruch wäre, aber die Aussage ist nicht so allgemein formulierbar).
Anschließend folgt ein Sprung zu einer anderen Methode ($(q-1) \mid (p^2-1)$-Argument),
der nicht vollständig ausgeführt wird. ⚠️

### Lemma 4.1 (Schranke für $q$)
Als „Beweisskizze" markiert. Die Skizze verwendet:
$(q-1) \mid (pr - 1)$, $r \le pq$, also $(q-1) \mid (p^2 q - 1)$, damit
$(q-1) \mid (p^2 - 1)$.

Das Argument „$r \le pq$" ist zwar wahr ($r \mid (pq-1)/k$ ergibt keine scharfe Schranke
in diese Richtung), aber die Ableitung $(q-1) \mid (p^2q - 1)$ durch Einsetzen eines
beliebigen $r \le pq$ ist **logisch nicht korrekt** — man kann nicht einfach $r = pq$
einsetzen, da das kein zulässiger Wert ist.

Die **Schlussformel $q \le \frac{3(p-1)}{p-2}$** ist inhaltlich korrekt (folgt sauber aus
dem englischen Beweis), aber der deutsche Beweissketch liefert dafür **keine valide Herleitung**.

❌ **Beweis-Lücke:** Lemma 4.1 (deutsche Version) ist nur eine Skizze und der
Beweis-Weg ist nicht korrekt ausgeführt.

---

## 7. Zusammenfassung der Befunde

| Komponente | EN-Version | DE-Version | Anmerkung |
|---|---|---|---|
| Quadratfreiheit | ✅ | ✅ | Korrekt |
| Kein Semiprim | ✅ (⚠️) | ✅ | EN: Lücke für $p \ge 7$ |
| Fall n=2qr | ✅ | ✅ (❌) | DE: Formel-Darstellungsfehler |
| Key Identity | ✅ | – | Nicht in DE-Version |
| Schritt 1: Reduktion | ✅ | ✅ | Korrekt |
| Schritt 2: k-Schranke | ✅ | ✅ | Korrekt |
| Schritt 3: q-Schranke | ✅ | ⚠️ | DE-Beweis unvollständig |
| Schritt 4: Fallanalyse | ✅ | ✅ | Korrekt |
| Schritt 5: (3,5)-Ausschluss | ✅ | ✅ | Korrekt |
| Korollar 4 Primfaktoren | ✅ | ✅ | Korrekt |

---

## 8. Empfehlung

**Englische Version:** Annahme nach kleinerer Revision.
- Ungleichungskette für Semiprim-Fall ($p \ge 7$) explizit schließen.

**Deutsche Version:** Überarbeitung erforderlich.
1. Formeldarstellung „$r = 2q - 1 + 1 - 1 = 2q$" korrigieren zu „$(r-1) = 2q-1$, also $r = 2q$".
2. Lemma 4.1 (Schranke $q$): Entweder vollständigen Beweis nachliefern oder auf die
   englische Version verweisen. Der aktuelle „Beweissketch" ist mathematisch nicht valide.
3. Die Key Identity-Technik (aus dem englischen Paper) in die deutsche Version integrieren —
   sie macht den Beweis erheblich transparenter.
