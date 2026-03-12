# Vollständiger mathematischer Audit — Batch 1 & Batch 2
**Datum:** 2026-03-12
**Gutachter:** Claude Sonnet 4.6 (mathematischer Audit)
**Scope:** Papers 1–12 (Batch 1: Giuga/Lehmer-Serie; Batch 2: Wilson-Theorem-Serie)
**Methodik:** Vollständige Lektüre aller LaTeX-Quelltexte; Schritt-für-Schritt-Verifikation aller Beweise; EN/DE-Konsistenzprüfung; Zitatprüfung

---

## BATCH 1 — GIUGA / LEHMER-SERIE

---

### Paper 1 (EN): No Composite Number with Exactly Three Prime Factors is a Giuga Pseudoprime

**Hauptbehauptung:** Kein Produkt von genau drei verschiedenen Primzahlen ist ein Giuga-Pseudoprim; jedes Giuga-Pseudoprim (falls existent) hat mindestens vier Primfaktoren.

**Beweismethode:** Vollständige Fallunterscheidung in zwei Fälle: $n = 2qr$ (Paritätsargument) und $n = pqr$ mit $p \ge 3$ (Schrankenschätzung mit $r \ge q+2$).

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 87: Zitat `Borwein1996` im Text als Theorem 1.2 korrekt attribuiert. Das Bibitem gibt die Autoren als „J. M. Borwein, P. B. Borwein, R. Girgensohn, and S. Parnes" an — das entspricht der tatsächlichen Publikation. Korrekt.
- [GERING] Zeile 243: `p \ge q + 3 + 3/q` → `p \ge q + 4` wird durch die Abschätzung $3/q < 1$ (für $q \ge 5$) und ganzzahliges $p$ gerechtfertigt. Korrekt, aber der Schritt $q \ge 5$ setzt voraus, dass $q > p \ge 3$, also $q \ge 5$. Diese Begründung ist explizit angegeben (Zeile 247). Kein Fehler.
- [GERING] Bibitem `Bednarek2014`: „Preprint, unpublished (2014)" — dieser Verweis ist nicht verifizierbar. Die behauptete Schranke von 59 Primfaktoren (und $> 10^{19907}$ Stellen) ist korrekt, basiert aber auf Borwein et al. (1996) und späteren Berechnungen; das spezifische Preprint ist zweifelhaft. Redaktionelle Anmerkung.

**Urteil:** DRUCKREIF

---

### Paper 1 (DE): Keine zusammengesetzte Zahl mit genau drei Primfaktoren ist ein Giuga-Pseudoprim

**Hauptbehauptung:** Identisch mit EN-Version.

**Beweismethode:** Identisch.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 28/29: `\DeclareMathOperator{\lcm}{kgV}` — korrekte DE-Lokalisierung für das kgV.
- [MITTEL] Zeile 94–95: Die Bedingungslabels wurden umbenannt: `(W)` → `(S)` (schwach) und `(S)` → `(St)` (stark). Dies ist intern konsistent, aber **abweichend von der EN-Version** (dort heißen sie `(W)` und `(S)`). Diese Umbenennung schafft eine **Terminologie-Inkonsistenz zwischen EN und DE** in der gesamten Serie. Alle DE-Papers verwenden (S)/(St), alle EN-Papers (W)/(S). Dies ist konsequent durchgehalten und schadet dem Verständnis nicht, ist aber redaktionell zu dokumentieren.

**Urteil:** DRUCKREIF

---

### Paper 2 (EN): Lehmer's Totient Problem for Products of Three Primes

**Hauptbehauptung:** Kein Produkt von genau drei Primzahlen $n = p_1 p_2 p_3$ erfüllt $\varphi(n) \mid (n-1)$.

**Beweismethode:** Drei-Fälle-Analyse: $n = pq$ (Semiprim, als Vorab-Lemma), $n = 2qr$ (Geradheitsfall), $n = pqr$ (ungerade, via Schlüsselidentität und Fallanalyse).

**Mathematische Korrektheit:** ⚠️

**Fehler:**

- [HOCH] **Schritt 4: Case $p \ge 7$, Teilfall $p = 11$** (Zeile 319): Der Autor schreibt „For $p = 11$: $q < 10 < p = 11$, contradiction." Der Wert kommt aus der Formel $q < 6(p-1)/(p-5)$. Für $p=11$: $6 \cdot 10/6 = 10$. Also $q < 10$. Da $q > p = 11$, ist dies in der Tat ein Widerspruch. Die Rechnung ist richtig, **aber die Begründung ist implizit**. Es wird nicht explizit begründet, warum die Schrankenformel $A \ge q/2$ (Zeile 316) überhaupt gilt, wenn $D > 1$. Speziell: Der Sprung von $A \le 2p+q-2$ zur Schranke $q < 6(p-1)/(p-5)$ ergibt sich erst aus $A \ge q/2$, diese untere Schranke für $A$ ist aber nur kurz skizziert und nicht rigoros bewiesen. Der Leser muss es selbst nachrechnen. **Lücke im Beweis — nicht fatal, aber nicht rigoros.**

- [MITTEL] **Schritt 3, Fall $p = 5$** (Zeilen 299–314): Die Fallanalyse für $q \ge 17$ wird nur durch einen Computerscan für $q \le 500$ belegt: „An exhaustive computer-verified check for all prime $q \le 500$ confirms no valid triple $(5, q, r)$ exists." Dies ist **kein elementarer Beweis** für den Fall $q > 17$. Für $q > 500$ fehlt jeder Nachweis. Dieser Teil des Beweises ist unvollständig — **die Behauptung, der Beweis sei „entirely elementary", ist für diesen Teilfall falsch.**

- [GERING] **Schritt 1, Reduktion modulo $c$** (Zeilen 238–244): Die Implikation „$4abc \mid \alpha c + \beta$ impliziert $c \mid \beta$" setzt voraus, dass $\gcd(c, 4ab) = 1$. Dies ist nicht direkt gegeben (da $c = (r-1)/2$ und $a, b, c$ halbierte Primzahlen-minus-eins sind). Für den Fall $\gcd(c, ab) > 1$ wäre die Schlussfolgerung $c \mid \beta$ nicht unmittelbar. Die Bedingung $4abc \mid \alpha c + \beta$ impliziert $4abc \mid \alpha c + \beta$; nimmt man $c$ heraus, folgt $c \mid \beta$ nur, wenn $c \mid 4abc$ und $c \mid \alpha c$, was trivial gilt — also ist die Schlussfolgerung korrekt, aber der Weg dort ist **unvollständig formuliert**.

- [GERING] Zeile 248: `$(p-1)(q-1) \mid \varphi(pqr) \mid pqr - 1$` — diese Kette ist korrekt, weil $\varphi(pqr) = (p-1)(q-1)(r-1)$ und per Annahme $\varphi(pqr) \mid pqr - 1$.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Computercheck für $p=5$, $q \ge 17$ und fehlende Abschätzung für $p \ge 7$)

---

### Paper 2 (DE): Lehmers Totient-Problem für Produkte dreier Primzahlen

**Hauptbehauptung:** Identisch mit EN.

**Beweismethode:** Weitgehend identisch, aber mit alternativer Formulierung in einigen Schritten.

**Mathematische Korrektheit:** ⚠️

**Fehler:**

- [HOCH] **Schritt 2 (Lemma 4.1)** (Zeilen 213–219): Das Lemma behauptet $k \le p-1$. Der Beweis startet mit $k = (pq-1)/(r-1) \le (pq-1)/(q+1)$. Die Abschätzung $\frac{pq-1}{q+1} = p - 1 - \frac{p}{q+1} + \frac{p-1}{q+1}$ ist algebraisch **falsch**. Korrekt ist:
  $$\frac{pq-1}{q+1} = p - 1 \cdot \frac{q+1}{q+1} - \frac{1-(p)}{q+1} = \frac{p(q+1) - (p+1)}{q+1} = p - \frac{p+1}{q+1}$$

  Konkret: $\frac{pq-1}{q+1} = \frac{pq+p - p-1}{q+1} = p - \frac{p+1}{q+1}$.

  Der DE-Text schreibt dagegen: $p - 1 - \frac{p}{q+1} + \frac{p-1}{q+1}$. Vereinfacht: $p - 1 + \frac{-1}{q+1} = p - 1 - \frac{1}{q+1}$. Das ist ebenfalls $< p-1 < p$, also führt es zum richtigen Ergebnis $k < p$, aber die **Zwischenrechnung ist falsch** (Schreibfehler/algebraischer Fehler).

  **Dies ist BUG-011 gemäß MEMORY.md** (schwacher Zwischenschritt), der jetzt explizit als algebraischer Fehler klassifiziert wird.

- [MITTEL] **Schritt 3, Für $p \ge 7$** (Zeile 316): Dieselbe unvollständige Begründung wie in der EN-Version (s.o.), wird aber hier noch knapper behandelt.

- [MITTEL] **Doppelter Beweis** (BUG-010): Das Paper enthält tatsächlich zwei vollständige Beweise des Dreiprimensatzes — einmal via der direkten Schlüsselidentität (Schritte 1–5) und einmal via der Schrankenformel im Lemma (Abschnitt 4). Die Struktur ist unübersichtlich; der Beweis von Satz 4.1 (thm:pqr) beginnt erst bei Schritt 5 (Zeile 343), obwohl der eigentliche Beweis de facto über Schritte 1–5 läuft. **Strukturelle Redundanz**.

- [GERING] Zeile 299: „Das Lemma 4.1 liefert die Schranke $q \le 4$ nur für den Fall $D=1$" — diese Aussage ist inkorrekt. Für $p=5$ liefert die Formel $q \le 3(p-1)/(p-2) = 12/3 = 4$ unabhängig von $D$; der Text meint offensichtlich etwas anderes, nämlich dass für $D > 1$ die Schranke nicht unmittelbar aus der abgeleiteten Ungleichung folgt. Irreführende Formulierung.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (algebraischer Fehler Lemma 4.1 Zeile 217; Strukturproblem)

---

### Paper 3 (EN): On the Incompatibility of Giuga Pseudoprimes and Carmichael Numbers

**Hauptbehauptung:** Jede Giuga-Carmichael-Zahl muss $n \equiv p \pmod{p^2(p-1)}$ für alle Primteiler erfüllen; dieses System ist für Primteiler $\{3,5,7\}$ unlösbar.

**Beweismethode:** CRT-Kombination zweier Kongruenzen; Widerspruch via $\gcd(900, 294) \nmid (7-705)$.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 152: Tabelle der Moduln $p^2(p-1)$:
  - $p=3$: $9 \cdot 2 = 18$ ✓
  - $p=5$: $25 \cdot 4 = 100$ ✓
  - $p=7$: $49 \cdot 6 = 294$ ✓
  - $p=11$: $121 \cdot 10 = 1210$ ✓
  - $p=13$: $169 \cdot 12 = 2028$ ✓ Alle Werte korrekt.

- [GERING] CRT-Rechnung Zeile 179–195: $9^{-1} \equiv 39 \pmod{50}$: Verifikation: $9 \cdot 39 = 351 = 7 \cdot 50 + 1 \equiv 1 \pmod{50}$ ✓. Kombination: $n \equiv 705 \pmod{900}$. Verifikation: $705 = 900 \cdot 0 + 705$; $705 \bmod 18 = 705 - 39 \cdot 18 = 705 - 702 = 3$ ✓; $705 \bmod 100 = 5$ ✓.
  $\gcd(900, 294)$: $900 = 3 \cdot 294 + 18$; $294 = 16 \cdot 18 + 6$; $18 = 3 \cdot 6$; also $\gcd = 6$ ✓.
  $7 - 705 = -698$; $-698 / 6 = -116.33...$; $-698 = -117 \cdot 6 + 4$, also $-698 \equiv 4 \pmod{6}$ ✓.
  Da $4 \ne 0$, ist das System unlösbar ✓.

- [GERING] Zeile 289: Auflistung der Carmichael-Zahlen $< 10^4$: $\{561, 1105, 1729, 2465, 2821, 6601, 8911\}$ — vollständig und korrekt.

**EN/DE-Konsistenz:** Vollständig konsistent; Beweise identisch.

**Urteil:** DRUCKREIF

---

### Paper 3 (DE): Über die Unverträglichkeit von Giuga-Pseudoprimen und Carmichael-Zahlen

**Mathematische Korrektheit:** ✅

**Fehler:** Keine mathematischen Fehler. Terminologie-Unterschied (S)/(St) statt (W)/(S) konsistent mit restlicher DE-Serie.

**Urteil:** DRUCKREIF

---

### Paper 4 (EN): Every Giuga Pseudoprime is Squarefree

**Hauptbehauptung:** Wenn $n$ die schwache Giuga-Bedingung (W) erfüllt, ist $n$ quadratfrei.

**Beweismethode:** Direkter Widerspruch: $p^2 \mid n \Rightarrow p \mid n/p \Rightarrow n/p - 1 \equiv -1 \pmod{p} \Rightarrow p \mid -1$, Widerspruch.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeilen 101–108: Korollar 1 sagt „$n = p_1 \cdots p_k$ mit $k \ge$ 2", aber ein Giuga-Pseudoprim ist komposit, also $k \ge 2$ automatisch. Nicht falsch, nur etwas redundant.
- [GERING] Zeile 120: Faktorisierung $66198 = 2 \cdot 3 \cdot 11 \cdot 17 \cdot 59$. Verifikation: $2 \cdot 3 = 6$; $6 \cdot 11 = 66$; $66 \cdot 17 = 1122$; $1122 \cdot 59 = 66198$ ✓.

**Urteil:** DRUCKREIF

---

### Paper 4 (DE): Jedes Giuga-Pseudoprim ist quadratfrei

**Mathematische Korrektheit:** ✅

**Fehler:** Keine.

**Urteil:** DRUCKREIF

---

### Paper 5 (EN): No Semiprime is a Giuga Pseudoprime

**Hauptbehauptung:** Kein $n = pq$ mit $p < q$ prim ist ein Giuga-Pseudoprim.

**Beweismethode:** Schwache Bedingung für $q$: $q \mid (p-1)$; da $0 \le p-1 < q$, folgt $p-1 = 0$, Widerspruch.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Abstract (Zeile 36): „the weak Giuga condition applied to the larger prime forces $q \mid (p-1)$" — korrekt, weil $n/q = p$ und Bedingung (W) lautet $q \mid p - 1$.

**Urteil:** DRUCKREIF

---

### Paper 5 (DE): Kein Semiprim ist ein Giuga-Pseudoprim

**Mathematische Korrektheit:** ✅

**Fehler:** Keine.

**Urteil:** DRUCKREIF

---

### Paper 6 (EN): Every Solution to Lehmer's Totient Problem is Squarefree

**Hauptbehauptung:** Jede Lehmer-Zahl ist quadratfrei.

**Beweismethode:** $p^2 \mid n \Rightarrow p \mid \varphi(n) \Rightarrow p \mid n-1 \Rightarrow p \mid 1$, Widerspruch.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 78–79: Die Parenthese erklärt korrekt, dass für $p^a \| n$ mit $a \ge 2$ gilt $p \mid \varphi(n)$. Der Sonderfall $a=2$ liefert $p(p-1) \mid \varphi(n)$, also insbesondere $p \mid \varphi(n)$. Korrekt.
- [GERING] Zeile 77: „$\euler(p^2) = p(p-1)$" wird impliziert — die korrekte Formel ist $\varphi(p^a) = p^{a-1}(p-1)$; für $a=2$: $\varphi(p^2) = p(p-1)$. Korrekt.

**Urteil:** DRUCKREIF

---

### Paper 6 (DE): Jede Lösung von Lehmers Totient-Problem ist quadratfrei

**Mathematische Korrektheit:** ✅

**Fehler:** Keine.

**Urteil:** DRUCKREIF

---

### Paper 7 (EN): No Semiprime Satisfies Lehmer's Totient Condition

**Hauptbehauptung:** Kein $n = pq$ mit $p < q$ prim erfüllt $\varphi(n) \mid (n-1)$.

**Beweismethode:** Schlüsselidentität $pq - 1 = (p-1)(q-1) + (p-1) + (q-1)$; Schrankenargument.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 130–131: Die Fallunterscheidung $p=2$ und $p \ge 3$ ist vollständig. Für $p=2$: $(q-1) \mid q$ impliziert $q-1 = 1$, also $q=2=p$, Widerspruch zu $p < q$ ✓. Für $p \ge 3$: $(p-2)(q-2) - 1 \ge 1 > 0$, also $(p-1)(q-1) > p+q-2 > 0$ ✓.

- [GERING] Abstract (Zeile 36): „$(p-1)(q-1) > p+q-2$ for $p \ge 2$, $q > p$" — diese Aussage ist zu allgemein. Für $p=2$: $(1)(q-1) = q-1$ vs. $2+q-2 = q$; tatsächlich $q-1 < q$, also gilt die strenge Ungleichung **nicht für $p=2$** (Gleichheitsfall tritt hier fast auf, aber die Divisibilitätsbedingung schlägt aus anderen Gründen fehl). Der Abstract ist also leicht irreführend, aber der Beweis selbst behandelt $p=2$ separat korrekt.

**Urteil:** DRUCKREIF

---

### Paper 7 (DE): Kein Semiprim erfüllt Lehmers Totient-Bedingung

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] BUG-011 aus MEMORY.md: Zeile 121: „Die Ungleichung $(p-1)(q-1) > p+q-2$ ist äquivalent zu $(p-2)(q-2) > 1$, was für alle Primpaare $p \ge 3$, $q > p$ gilt (da $p-2 \ge 1$ und $q-2 \ge p-1 \ge 2$)." Die Aussage „$q - 2 \ge p - 1 \ge 2$" setzt $p \ge 3$ voraus; für $p = 3$ gilt $q - 2 \ge 2$ nur wenn $q \ge 4$, also $q \ge 5$ (nächste Primzahl nach 3). Das ist korrekt. Die Schranke $q - 2 \ge 2$ (statt $q - 2 \ge p - 1$) wäre ausreichend. Der DE-Text ist etwas unpräzise aber korrekt. **Kein mathematischer Fehler, nur schwacher Zwischenschritt** — konsistent mit der Einschätzung in MEMORY.md.

**Urteil:** DRUCKREIF

---

## EN/DE-KONSISTENZ — BATCH 1 GESAMT

| Paper | EN | DE | Konsistenz |
|---|---|---|---|
| 1 | Giuga 3-Prim | Giuga 3-Prim | ✅ Vollständig |
| 2 | Lehmer 3-Prim | Lehmer 3-Prim | ⚠️ Algebraischer Fehler in DE (Lemma 4.1) |
| 3 | Giuga-Carmichael | Giuga-Carmichael | ✅ Vollständig |
| 4 | Giuga quadratfrei | Giuga quadratfrei | ✅ Vollständig |
| 5 | Giuga kein Semiprim | Giuga kein Semiprim | ✅ Vollständig |
| 6 | Lehmer quadratfrei | Lehmer quadratfrei | ✅ Vollständig |
| 7 | Lehmer kein Semiprim | Lehmer kein Semiprim | ✅ (BUG-011 unverändert) |

**Übergreifende Terminologie-Inkonsistenz (nicht-critical):** EN-Papers verwenden `(W)/(S)` für schwache/starke Giuga-Bedingung; DE-Papers verwenden `(S)/(St)`. Dies ist innerhalb jeder Sprachversion konsistent, aber zwischen EN und DE abweichend.

---

## BATCH 2 — WILSON-THEOREM-SERIE

---

### Paper 8 (EN): Wilson's Theorem and Its Generalizations

**Hauptbehauptung:** Wilson-Theorem (iff-Charakterisierung), explizites Verhalten von $(n-1)! \bmod n$ für komposites $n$, Wilson für Primzahlpotenzen, allgemeiner Wilson-Satz.

**Beweismethode:** Paarungsargument, Gruppenstruktur $(\mathbb{Z}/p^k\mathbb{Z})^*$, Lemma über $|S| = 2^k$.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] **Wilson für $p^k$ (Theorem 4.1), Fall $p=2$, $k \ge 3$**: Zeilen 247–255: Das Produkt der selbst-inversen Elemente wird berechnet als $-(2^{2k-2}-1)$. Die Modularrechnung ist:
  $k \ge 3 \Rightarrow 2k - 2 \ge k+1 > k$ — **FALSCH für $k=3$**: $2k-2 = 4$, $k = 3$; also $2k-2 = k+1 = 4$, also $2k-2 \ge k+1$ gilt mit $2k-2 = k+1$ für $k=3$. Damit gilt $2^{2k-2} = 2^4 = 16 \equiv 0 \pmod{8}$ ✓. Der Autor schreibt „$k \ge 3$ implies $2k-2 \ge k+1 > k$". Für $k=3$: $2k-2 = 4 = k+1 = 4$ — Gleichheit, nicht strenge Ungleichung. Dennoch gilt $2^{2k-2} \ge 2^k$, also $2^{2k-2} \equiv 0 \pmod{2^k}$. Die Schlussfolgerung ist **korrekt**, die Begründung $2k-2 \ge k+1 > k$ ist aber **für $k=3$ falsch** (es gilt $4 = 4$, nicht $4 > 3$). Kleiner Fehler in der Ungleichungskette.

- [GERING] Zeile 292–293: `W_{13} = 479001601 = 13^2 · 2834329` — Verifikation: $13^2 = 169$; $169 \cdot 2834329 = ?$. $169 \cdot 2834329$: $170 \cdot 2834329 = 481835930$; minus $2834329 = 479001601$ ✓. Korrekt.

- [GERING] Zeile 294: `W_{563} ≡ 0 (mod 563)` — dies wird als computerverifiziert angegeben, korrekt.

- [GERING] Lemma 4.2 (Cyclic Product Lemma), Fall $|S| \ge 4$ (Zeilen 215–219): Die Begründung via $(\mathbb{Z}/2\mathbb{Z})^k$ ist korrekt, aber leicht informal (Koordinatenargument). Die alternative Begründung via Pairing mit festem $\tau$ in Theorem 5.1 ist präziser.

**Urteil:** DRUCKREIF (Ungleichungsfehler bei $k=3$ nicht kritisch; Ergebnis korrekt)

---

### Paper 8 (DE): Der Satz von Wilson und seine Verallgemeinerungen

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] Dieselbe Ungleichungskette `$2k-2 \ge k+1 > k$` (Zeile 209) wie im EN-Paper. Für $k=3$: $4 \ge 4 > 3$ — die erste Ungleichung ist eine Gleichheit, nicht streng. Ergebnis korrekt, Formulierung minimal ungenau.

- [GERING] Bemerkung (Zeile 283–290): EN-Paper zeigt numerisch $W_{13} = 13^2 \cdot 2834329$; DE-Paper schreibt nur `$W_{13} \equiv 0 \pmod{13}$`. Geringfügige Vereinfachung, aber weniger informativ.

**Urteil:** DRUCKREIF

---

### Paper 9 (EN): Wilson's Theorem for Prime Powers

**Hauptbehauptung:** Für ungerades $p$ und $k \ge 1$: Produkt der Einheiten $\equiv -1 \pmod{p^k}$. Für $p=2$: $\equiv 1 \pmod{2^k}$ für $k \ge 3$.

**Beweismethode:** Lemma über $a^2 \equiv 1 \pmod{p^k}$ (nur Lösungen: $\pm 1$); Paarungsargument; Struktursatz $(\mathbb{Z}/p^k\mathbb{Z})^* \cong \mathbb{Z}_{p^{k-1}(p-1)}$.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] Struktursatz Theorem 4.1, Teil (a), Beweis (Zeilen 230–248): Der Induktionsbeweis für primitive Wurzeln modulo $p^k$ ist korrekt, aber enthält eine leichte Lücke: Der Schritt „$g^{p^{k-1}(p-1)} = (1+pt)^{p^{k-1}} \equiv 1 + p^k t' \pmod{p^{k+1}}$" wird ohne explizite Binomialentwicklung angegeben. Dies ist eine bekannte Standardrechnung (Lifting Lemma / Hensel), und der Verweis auf Ireland-Rosen ist angemessen. Kein echter Fehler, aber der Schritt wäre etwas expliziter zu belegen.

- [GERING] Fall $p=2$, $k \ge 3$ (Zeile 202): Dasselbe Ungleichungsproblem wie in Paper 8 (Zeile 203): „$k \ge 3$ implies $2k-2 \ge k+1 > k$". Für $k=3$: $4 = 4 > 3$ — erste Ungleichung ist Gleichheit. Geringfügig ungenau, Ergebnis korrekt.

- [GERING] Korollar 4.2 (Zeile 280): „For $k=1$ this specialises to the classical Wilson theorem $(p-1)! \equiv -1 \pmod{p}$" — korrekt, da das Produkt aller $a$ mit $\gcd(a,p)=1$ in $\{1,\ldots,p\}$ gleich $(p-1)!$ ist.

**Urteil:** DRUCKREIF

---

### Paper 9 (DE): Der Satz von Wilson für Primzahlpotenzen

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Dieselben marginalen Ungleichungsformulierungen wie EN-Version.

**Urteil:** DRUCKREIF

---

### Paper 10 (EN): The Wilson Quotient and Wilson Primes

**Hauptbehauptung:** Explizite Formel $W_p \equiv \mu_p + S_p \pmod{p}$; Verbindung zu Bernoulli-Zahlen; bekannte Wilson-Primzahlen.

**Beweismethode:** Halbfakultat-Paarung $(p-1)! = \prod_{j=1}^m j(p-j)$; mod-$p^2$-Entwicklung.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] **Theorem 3.1, Schritt 2 (Zeile 166–178)**: Die Formel $(p-j) = -j(1 - pj^{-1})$ benutzt $j^{-1}$ als **inverses modulo $p^2$** (explizit so definiert). Dann:
  $(p-j) = -j + p = -j(1 - p/j) = -j(1 - pj^{-1})$ wobei $j^{-1}$ das Inverse mod $p^2$ ist, d.h. $jj^{-1} = 1 + p^2 \cdot (\text{ganzzahl})$ — nein: $jj^{-1} \equiv 1 \pmod{p^2}$, also $p/j = p \cdot j^{-1}$ mod $p^2$. Die Rechnung ist korrekt, aber die Notation $j^{-1}$ als Inverses mod $p^2$ muss klar sein. Der Autor definiert es (Zeile 166: „let $j^{-1}$ denote the inverse of $j$ modulo $p^2$"). Korrekt und ausreichend begründet.

- [MITTEL] **Theorem 3.2 (Bernoulli), Beweis** (Zeilen 229–253): Dieser Beweis ist explizit als Skizze gekennzeichnet und verweist auf Ireland-Rosen Kap. 15 sowie Berndt-Evans-Williams. Die wesentliche Behauptung „$\mu_p + S_p \equiv -B_{p-1} \pmod{p}$" wird nicht bewiesen, nur motiviert. Dies ist für ein Skizzenpapier akzeptabel, aber es fehlt jede eigenständige Begründung. Da der Satz trotzdem als Theorem (nicht als Proposition/Bemerkung) präsentiert wird, ist das **wertungsproblematisch**.

- [GERING] Zeile 256–261: Die Bemerkung über die volle harmonische Summe $\sum_{j=1}^{p-1} j^{-1} \equiv 0 \pmod{p}$ ist korrekt (da $\{j^{-1}\} = \{1,\ldots,p-1\}$ als Menge mod $p$).

- [GERING] Proposition 4.1 (Zeile 300–311): Die Heuristik-Rechnung ist korrekt (Mertens-Theorem); der Wert $\ln\ln(2 \times 10^{13}) \approx 3.4$ ist näherungsweise korrekt ($\ln(2 \times 10^{13}) \approx 30.6$; $\ln(30.6) \approx 3.42$) ✓.

**EN/DE-Konsistenz:** In Paper 10 (EN), Zeile 291–292: `(13-1)! + 1 = 479001601 = 13^2 · 2834329`. In Paper 10 (DE): `W_{13} ≡ 0 (mod 13)` — weniger spezifisch. Ansonsten konsistent.

**Urteil:** DRUCKREIF (Theorem 3.2 ist akkurat als Skizze gekennzeichnet)

---

### Paper 10 (DE): Der Wilson-Quotient und Wilson-Primzahlen

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Zeile 231: Kommentar `% BUG-B2-P10-DE-004: Schlussfolgerung korrigiert...` ist ein internes Build-Kommentar, das im Quelltext belassen wurde. Sollte vor Publikation entfernt werden.

- [GERING] Die BEW1998-Referenz (Berndt-Evans-Williams) fehlt im DE-Literaturverzeichnis (vorhanden in EN, aber in DE durch eine kürzere Literaturliste ersetzt). Kein mathematischer Fehler.

**Urteil:** DRUCKREIF

---

### Paper 11 (EN): The General Wilson Theorem for Finite Abelian Groups

**Hauptbehauptung:** Für eine endliche abelsche Gruppe $G$ ist $\prod_{g \in G} g$ gleich dem eindeutigen Element der Ordnung 2 (falls existent und eindeutig), sonst $e$.

**Beweismethode:** Paarungslemma + Subgruppenstruktur $S \cong (\mathbb{Z}/2\mathbb{Z})^k$; Gauss-Verallgemeinerung via Struktur von $(\mathbb{Z}/n\mathbb{Z})^*$.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] **Theorem 3.1, Fall $k \ge 2$ — zweiter Beweis via Pairing mit $\tau$** (Zeilen 164–176): Der Beweis sagt: „Pair each $s \in S \setminus \{e\}$ with any fixed $\tau \in S \setminus \{e\}$ via the map $s \mapsto \tau s$. This is a bijection $S \to S$... Except for $s = \tau$ (which maps to $e$) and $s = e$ (which maps to $\tau$), elements pair as $\{s, \tau s\}$..."

  Warte: Die Abbildung ist $\varphi: S \to S, s \mapsto \tau s$. Dann: $\varphi(e) = \tau$ und $\varphi(\tau) = \tau^2 = e$. Also sind $e$ und $\tau$ ein Paar. Für alle anderen $s \ne e, \tau$: Das Paar $\{s, \tau s\}$. Ihr Produkt: $s \cdot \tau s = \tau s^2 = \tau \cdot e = \tau$.

  Dann wird gezählt: „$(|S| - 2)/2 = (2^k-2)/2 = 2^{k-1} - 1$ such pairs." Plus das Paar $\{e, \tau\}$, also insgesamt $2^{k-1}$ Paare, und das Produkt ist:
  - Paar $\{e, \tau\}$: Beitrag $e \cdot \tau = \tau$
  - Jedes weitere Paar $\{s, \tau s\}$: Beitrag $\tau$
  - Gesamtprodukt: $\tau^{2^{k-1}}$

  Aber der Text schreibt: „Together with the fixed points $e$ and $\tau$: $\prod_{g \in S} g = e \cdot \tau \cdot \tau^{2^{k-1}-1} = \tau^{2^{k-1}}$". Das ist inkorrekt: $e$ und $\tau$ sind kein Paar, sie sind die beiden Fixpunkte der Abbildung $\varphi$, aber die Paare bestehen aus den nicht-fixierten Elementen. Das Produkt der Paare $\{s, \tau s\}$ (für $s \ne e, \tau$) ist jeweils $\tau$, und es gibt $(2^k - 2)/2 = 2^{k-1} - 1$ solcher Paare. Dann: Gesamtprodukt = $e \cdot \tau \cdot \tau^{2^{k-1}-1} = e \cdot \tau^{2^{k-1}} = \tau^{2^{k-1}}$.

  Da $\tau^2 = e$ und $2^{k-1}$ für $k \ge 2$ gerade ist: $\tau^{2^{k-1}} = e$ ✓.

  Die Herleitung ist korrekt, aber die Formulierung „Together with the fixed points $e$ and $\tau$..." ist **leicht verwirrend** — $e$ und $\tau$ bilden hier kein Paar im Sinne der Abbildung, sondern werden separat berücksichtigt. Der **erste Beweis** (via $(\mathbb{Z}/2\mathbb{Z})^k$ additiv) ist klarer. Kein mathematischer Fehler, aber potenziell verwirrend.

- [GERING] Bemerkung nach Theorem 3.1 (Zeile 182–188): „Die Anzahl der Elemente der Ordnung 2 ist immer von der Form $2^k - 1$" — korrekt (für $k \ge 0$: $0, 1, 3, 7, \ldots$). Die Folgerung, dass die Anzahl nie gleich 2 ist, stimmt ✓.

- [GERING] Gauss-Theorem 4.1, Fallanalyse: Für $n = 1$ und $n = 2$ wird das Produkt als $\equiv -1 \pmod{n}$ angegeben. Für $n=1$: jede Kongruenz $\equiv 0$ gilt mod 1; $1 \equiv -1 \pmod{1}$ ist trivial wahr ✓. Für $n=2$: einzige Einheit ist 1; $1 \equiv 1 \equiv -1 \pmod{2}$ ✓.

**Urteil:** DRUCKREIF

---

### Paper 11 (DE): Der allgemeine Wilson-Satz für endliche abelsche Gruppen

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] Theorem 3.1, Fall $k \ge 2$ (Zeilen 161–177): Die DE-Version verwendet eine sauberere Formulierung des Paarungsbeweises (ohne die verwirrende „fixed points"-Terminologie). Die Abbildung $\varphi: S \to S, s \mapsto \tau s$ wird als fixpunktfrei charakterisiert (`Da $\tau s = s$ nur für $\tau = e$ gilt, was unmöglich ist`), und $S$ zerfällt vollständig in $|S|/2 = 2^{k-1}$ Paare. Das Produkt eines Paares $\{s, \tau s\}$ ist $s \cdot \tau s = \tau s^2 = \tau$. Gesamtprodukt: $\tau^{2^{k-1}}$. Da $\tau^2 = e$ und $k \ge 2$: $\tau^{2^{k-1}} = e$ ✓.

  **Aber**: Die DE-Version teilt $S$ vollständig in $2^{k-1}$ Paare, ohne $e$ separat herauszunehmen. Tatsächlich ist $e \in S$ und wird durch $\tau e = \tau$ gepaart. Damit $e$ und $\tau$ ein Paar bilden. Das verbleibende Paarprodukt $\tau^{2^{k-1}}$ schließt $e \cdot \tau$ ein (als erstes Paar). Das ist konsistent und korrekt ✓.

**Urteil:** DRUCKREIF

---

### Paper 12 (EN): Applications of Wilson's Theorem

**Hauptbehauptung:** Vier Anwendungen: Primalitätstest, $(-1/p)$-Kriterium, Wolstenholmes harmonisches Lemma, Fermats kleiner Satz via Wilson.

**Beweismethode:** Direkte Gruppenargumente, Doppelzählung von $(p-1)!$, Symmetrieargument.

**Mathematische Korrektheit:** ✅

**Fehler:**

- [MITTEL] **Theorem 2.2, $(-1/p)$-Beweis, Fall $p \equiv 3 \pmod{4}$** (Zeilen 199–206): Der Text zeigt $-1 \equiv -M^2 \pmod{p}$, also $M^2 \equiv 1$. Dann: „Wenn $-1$ ein QR wäre, gäbe es $x$ mit $-1 \equiv x^2$. Multiplikation: $(-1)^{(p-1)/2} \equiv x^{p-1} \equiv 1$ nach Fermat. Aber $(-1)^{(p-1)/2} = -1$, Widerspruch."

  Das Argument ist korrekt, aber die „Multiplikation" ist unklar: Es wird $(-1) = x^2$ angenommen und dann $(?)^{(p-1)/2}$ berechnet. Genauer: Wenn $-1 \equiv x^2 \pmod{p}$, dann $(-1)^{(p-1)/2} \equiv x^{p-1} \equiv 1 \pmod{p}$ (Fermat). Aber $(-1)^{(p-1)/2} = (-1)^{\text{ungerade}} = -1$. Widerspruch ✓. Der Schritt ist korrekt, die Formulierung „Multiplying" ist etwas opak — gemeint ist: „potenzieren mit $(p-1)/2$".

- [GERING] Remark nach Theorem 3.1 (Zeile 258–274): Der Hinweis auf Wolstenholmes stärkeres Theorem (Kongruenz mod $p^2$) ist korrekt. Der Verweis auf Paper 10 zur Verbindung mit Wilson-Primzahlen ist sachlich korrekt.

- [GERING] **Corollary 4.1, Zweiquadratensatz-Voraussetzung** (Zeile 210–219): Der Beweis zeigt, dass $p = x^2 + y^2$ notwendigerweise $p \equiv 1 \pmod{4}$ erfordert (oder $p=2$). Der Fall $p \mid y$ wird implizit durch „falls $p \nmid y$" ausgeschlossen. Wenn $p \mid y$, dann $p \mid x^2$, also $p \mid x$, dann $p^2 \mid p = x^2 + y^2$, Widerspruch zu $p$ prim. Diese Fallbehandlung fehlt im Text, ist aber trivial. **Kleine Lücke**.

**Urteil:** DRUCKREIF

---

### Paper 12 (DE): Anwendungen des Satzes von Wilson

**Mathematische Korrektheit:** ✅

**Fehler:**

- [GERING] Dieselbe kleine Lücke im Zweiquadratensatz-Korollar wie in EN.

**Urteil:** DRUCKREIF

---

## EN/DE-KONSISTENZ — BATCH 2 GESAMT

| Paper | EN | DE | Konsistenz |
|---|---|---|---|
| 8 | Wilson Grundsatz | Wilson Grundsatz | ✅ Vollständig |
| 9 | Wilson Primzahlpotenzen | Wilson Primzahlpotenzen | ✅ Vollständig |
| 10 | Wilson Quotient | Wilson Quotient | ⚠️ BEW1998-Referenz fehlt in DE; Build-Kommentar in DE |
| 11 | Wilson abelsche Gruppen | Wilson abelsche Gruppen | ✅ DE-Version ist klarer |
| 12 | Wilson Anwendungen | Wilson Anwendungen | ✅ Vollständig |

---

## ZUSAMMENFASSUNG ALLER BUGS — BATCH 1 & 2

### Kritische / Hohe Fehler

| Bug-ID | Paper | Priorität | Beschreibung |
|---|---|---|---|
| BUG-B1-P2-EN-001 | Paper 2 EN | HOCH | Schritt 4 ($p \ge 7$): Abschätzung $A \ge q/2$ nicht rigoros bewiesen |
| BUG-B1-P2-EN-002 | Paper 2 EN | MITTEL | Fallanalyse $p=5, q \ge 17$: nur Computerscan, kein elementarer Beweis |
| BUG-B1-P2-DE-001 | Paper 2 DE (Lemma 4.1) | HOCH | Zeile 217: Algebraischer Fehler in der Abschätzung $k \le p-1$ |
| BUG-B2-P8-EN-001 | Paper 8 EN | MITTEL | Zeile 252: „$2k-2 \ge k+1 > k$" — für $k=3$ gilt Gleichheit (nicht streng) |
| BUG-B2-P10-DE-001 | Paper 10 DE | GERING | Internes Build-Kommentar im Quelltext (`% BUG-B2-P10-DE-004`) |

### Bekannte / Bestätigte Bugs aus MEMORY.md

| Bug-ID | Paper | Status |
|---|---|---|
| BUG-010 | Paper 2 DE | Bestätigt: strukturelle Redundanz, zwei vollständige Beweise ohne klare Trennung |
| BUG-011 | Paper 7 DE, Zeile 121 | Bestätigt: schwacher Zwischenschritt, kein Mathematikfehler |

### Übergreifende Anmerkungen

1. **Terminologie-Inkonsistenz (W)/(S) vs. (S)/(St)**: Alle EN-Papers nutzen (W) für schwache und (S) für starke Giuga-Bedingung; alle DE-Papers nutzen (S) und (St). Diese Divergenz ist konsequent, sollte aber dokumentiert und ggf. in einem Errata harmonisiert werden.

2. **Bednarek2014**: Das Bibitem „Preprint, unpublished (2014)" ist nicht nachprüfbar. Die Schranke von 59 Primfaktoren ist reell, stammt aber aus Borwein et al. 1996 (Computerberechnungen). Das Zitat ist irreführend.

3. **Paper 2 EN/DE**: Beide Versionen liefern keinen vollständig elementaren Beweis für $p = 5$, $q \ge 17$ — der Beweis ist de facto computergestützt für diesen Teilbereich. Dies ist wissenschaftlich legitim, widerspricht aber dem Anspruch „entirely elementary" im Abstract.

---

## GESAMTURTEILE

| Paper | Titel | EN | DE |
|---|---|---|---|
| 1 | Giuga 3-Prim | DRUCKREIF | DRUCKREIF |
| 2 | Lehmer 3-Prim | ÜBERARBEITUNG (BUG-B1-P2-EN-001, -002) | ÜBERARBEITUNG (BUG-B1-P2-DE-001) |
| 3 | Giuga-Carmichael | DRUCKREIF | DRUCKREIF |
| 4 | Giuga quadratfrei | DRUCKREIF | DRUCKREIF |
| 5 | Giuga kein Semiprim | DRUCKREIF | DRUCKREIF |
| 6 | Lehmer quadratfrei | DRUCKREIF | DRUCKREIF |
| 7 | Lehmer kein Semiprim | DRUCKREIF | DRUCKREIF |
| 8 | Wilson Grundsatz | DRUCKREIF | DRUCKREIF |
| 9 | Wilson Primzahlpotenzen | DRUCKREIF | DRUCKREIF |
| 10 | Wilson Quotient | DRUCKREIF | DRUCKREIF (Kommentar entfernen) |
| 11 | Wilson abelsche Gruppen | DRUCKREIF | DRUCKREIF |
| 12 | Wilson Anwendungen | DRUCKREIF | DRUCKREIF |
