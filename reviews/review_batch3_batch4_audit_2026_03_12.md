# Vollständiger mathematischer Audit — Batch 3 & 4 (Papers 13–20)

**Datum:** 2026-03-12
**Build:** 15
**Gutachter:** Claude Opus (Specialist Peer Review)
**Umfang:** 16 LaTeX-Dateien (Papers 13–20, jeweils EN + DE)

---

## Zusammenfassung

### Batch 3 (Papers 13–16): Sieblehre
- **Paper 13 (Brunscher Satz):** EN + DE mathematisch korrekt, druckreif.
- **Paper 14 (Selbergs Sieb):** EN + DE mathematisch korrekt, druckreif. Stilanmerkung: $|\lambda_d^*| \leq 1$ wird behauptet aber nicht bewiesen.
- **Paper 15 (Große Sieb-Ungleichung):** EN + DE mathematisch korrekt, druckreif. Build-Nummer veraltet (56 statt 80).
- **Paper 16 (Chens Satz):** EN + DE mathematisch korrekt, druckreif. Build-Nummer veraltet (56 statt 80).

### Batch 4 (Papers 17–20): Kreismethode
- **Paper 17 (Kreismethode):** EN + DE mathematisch korrekt, druckreif.
- **Paper 18 (Vinogradov Drei-Primzahlen):** EN enthält falschen Faktor 1/6 (Remark 2.1). DE enthält diesen Fehler nicht. Beide Versionen überarbeitungsbedürftig wegen dieses und weiterer bekannter Bugs.
- **Paper 19 (Waringsches Problem):** EN druckreif. DE hat fehlende Inhalte gegenüber EN.
- **Paper 20 (Goldbach Singuläre Reihe):** EN + DE druckreif. Geringe Mängel: Vorhersagewert n=100 in EN, Bibliografie-Inkonsistenz.

---

## Strukturierte Mängelliste

### BESTÄTIGTE BEKANNTE BUGS

| # | Paper-ID | Sprache | Schwere | Zeile | Bug-ID | Beschreibung | Fix-Empfehlung |
|---|----------|---------|---------|-------|--------|-------------|----------------|
| 1 | Paper 14 | EN | GERING | 341, 348 | BUG-B3-P14-EN-LAMBDA | Behauptung $\|\lambda_d^*\| \leq 1$ wird zweimal verwendet (Zeilen 341, 348), aber nie bewiesen. Die Aussage folgt aus der Cauchy-Schwarz-Optimierung der Selberg-Gewichte, wird aber nur verbal erwähnt ("which follows from the Cauchy-Schwarz bound on the optimal weights"). | Kurzen Beweis einfügen: Die optimalen Selberg-Gewichte sind $\lambda_d = \mu(d) \sum_{d\|l, l \leq z} \mu(l)/g(l) / \sum_{l \leq z} \mu^2(l)/g(l)$. Aus der Struktur folgt $\|\lambda_d\| \leq 1$. Alternativ: Verweis auf Halberstam-Richert, Sieve Methods, Chapter 3. |
| 2 | Paper 14 | DE | GERING | analog | BUG-B3-P14-DE-LAMBDA | Identischer Mangel in der deutschen Version. | Identisches Fix. |
| 3 | Paper 15 | EN | KOSMETISCH | 40 | BUG-B3-P15-BUILD | Build-Nummer zeigt "Build 56" statt aktueller Build 80. | `\address{Specialist Mathematics Project, Build 80}` |
| 4 | Paper 15 | DE | KOSMETISCH | analog | BUG-B3-P15-BUILD | Identisch in DE-Version. | Identisches Fix. |
| 5 | Paper 16 | EN | KOSMETISCH | 46 | BUG-B3-P16-BUILD | Build-Nummer zeigt "Build 56" statt aktueller Build 80. | Aktualisieren. |
| 6 | Paper 16 | DE | KOSMETISCH | analog | BUG-B3-P16-BUILD | Identisch in DE-Version. | Aktualisieren. |
| 7 | Paper 18 | EN | MITTEL | 136 | BUG-B4-P18-EN-002 | **Faktor 1/6 falsch.** In Remark nach Definition 2.1 steht: "$R_3(n) \sim (\log n)^3 r_3(n) / 6$". Hier ist $r_3(n) = \#\{(p_1,p_2,p_3): p_1+p_2+p_3 = n\}$ definiert als Zählung **geordneter** Tripel. Die Funktion $R_3(n) = \sum_{a+b+c=n} \Lambda(a)\Lambda(b)\Lambda(c)$ zählt ebenfalls **geordnete** Tripel. Daher ist der Faktor 1/6 (für Umrechnung geordnet→ungeordnet) hier deplatziert. Die korrekte Asymptotik ist $R_3(n) \sim (\log n)^3 \cdot r_3(n)$ (ohne Division durch 6), da beide Seiten geordnete Zählungen sind. | Zeile 136 ändern zu: "$R_3(n) \sim (\log n)^3 r_3(n)$" und den Klammerkommentar anpassen: "(each prime $p$ contributing roughly $\log p$ via the von Mangoldt weight)". |
| 8 | Paper 18 | DE | MITTEL | — | BUG-B4-P18-DE-003 | Die DE-Version enthält den 1/6-Fehler **nicht** explizit (keine analoge Remark-Passage). Allerdings ist die Konsistenz zum EN-Paper wichtig — nach Fix in EN sollte sichergestellt werden, dass auch DE keine implizite 1/6-Referenz enthält. | Nach Fix in EN: DE-Version auf Konsistenz prüfen. Derzeit kein aktiver Fehler in DE. **Status: Herabgestuft auf OFFEN/ABHÄNGIG.** |
| 9 | Paper 19 | DE | MITTEL | 222 | BUG-B4-P19-DE-001 | Wooleys "Effizientes Kongruenzrechnen" wird in Zeile 222 nur als `\begin{bemerkung}` eingeführt, während die EN-Version einen formalen `\begin{theorem}` mit vollständigem Beweissketch hat. | Bemerkung in Zeilen 222–235 zu einem Theorem hochstufen, analog zur EN-Version (dort: Theorem 5.3 mit Proof-Umgebung). |
| 10 | Paper 19 | DE | MITTEL | 245–254 | BUG-B4-P19-DE-002 | Es gibt zwar ein zweites `\begin{theorem}` (Zeile 245, "Effizientes Kongruenzrechnen, Wooley 2012–2016"), aber der Beweissketch fehlt im Vergleich zur EN-Version. | Beweissketch analog zur EN-Version einfügen (p-adische Induktionsstruktur, Hensel-Lift, Iteration). |
| 11 | Paper 19 | DE | HOCH | — | BUG-B4-P19-DE-003 | Die g(k)-Formel wird in der DE-Version weniger detailliert hergeleitet als in EN. Insbesondere fehlt die explizite Formel $g(k) = 2^k + \lfloor(3/2)^k\rfloor - 2$ mit vollständiger Herleitung (EN hat dies in Section 2 mit Dickson-Beweis). Die Tabelle für $G(k)$ hat zudem weniger Einträge als in EN. | Fehlende Herleitung der $g(k)$-Formel und zusätzliche Tabelleneinträge aus EN übernehmen. |
| 12 | Paper 20 | EN | GERING | 478 | BUG-B4-P20-EN-PREDICT | Tabelle zeigt für $n=100$: Pred. $r_2(n) \approx 8$, Actual = 6. Die Hardy-Littlewood-Vorhersage für $n=100$ mit $\mathfrak{S}(100) = 2C_2 \cdot 4/3 \approx 1.76$ ergibt: $r_2(100) \approx \mathfrak{S}(100) \cdot 100 / (2(\log 100)^2) \approx 1.76 \cdot 100 / (2 \cdot 21.2) \approx 4.15$. Weder 8 noch 6 stimmen mit der asymptotischen Formel exakt überein — für kleine $n$ ist die Abweichung erwartungsgemäß groß. Die DE-Version hat korrekt $\approx 6$. | EN-Wert von $\approx 8$ auf $\approx 6$ korrigieren (konsistent mit DE und dem tatsächlichen Wert). Alternativ: Fußnote ergänzen, dass die asymptotische Formel für kleine $n$ unzuverlässig ist. |
| 13 | Paper 20 | EN+DE | GERING | EN:556–559, DE:435–438 | BUG-B4-P20-BIBLIO | Vinogradov 1937 wird in Paper 18 als "Representation of an odd number as a sum of three primes" (Dokl. Akad. Nauk SSSR 15:169–172) zitiert, aber in Paper 20 DE als "Some theorems concerning the theory of primes" (Mat. Sb. (N.S.) 2(44):179–195). Dies sind **zwei verschiedene** Vinogradov-1937-Arbeiten. In Paper 20 DE ist die falsche zitiert: Für die Goldbach-Vermutung ist die Dokl.-Arbeit relevant, nicht die Mat. Sb.-Arbeit. Paper 20 EN zitiert Vinogradov 1937 gar nicht direkt. | In Paper 20 DE: Bibitem `Vinogradov1937` durch die Dokl.-Version ersetzen: "Darstellung einer ungeraden Zahl als Summe dreier Primzahlen. Dokl. Akad. Nauk SSSR, 15:169–172, 1937." (konsistent mit Paper 18 DE). |

### BESTÄTIGTE BEHOBENE BUGS (Verifiziert)

| Bug-ID | Paper | Status | Verifizierung |
|--------|-------|--------|---------------|
| BUG-B4-P18-EN-001 | Paper 18 EN | BEHOBEN | Vorzeichen in Remark 3.6 korrekt (Build 80) |
| BUG-B4-P18-DE-001 | Paper 18 DE | BEHOBEN | Identischer Fix in DE (Build 80) |
| BUG-B4-P18-DE-002 | Paper 18 DE | BEHOBEN | Vollständige Fallunterscheidung Simplex-Integral vorhanden |
| BUG-B4-P19-EN-001 | Paper 19 EN | BEHOBEN | Querverweise korrekt |
| BUG-B4-P20-EN-001 | Paper 20 EN | BEHOBEN | Kein Satzfragment mehr vorhanden |
| BUG-B4-P20-EN-002 | Paper 20 EN | BEHOBEN | Tabellenwerte (außer n=100) korrekt |
| BUG-B4-P20-DE-001 | Paper 20 DE | BEHOBEN | Tabellenwerte korrekt |

### NEUE BEFUNDE

| # | Paper-ID | Sprache | Schwere | Zeile | Bug-ID | Beschreibung | Fix-Empfehlung |
|---|----------|---------|---------|-------|--------|-------------|----------------|
| 14 | Paper 20 | DE | GERING | 436–438 | BUG-B4-P20-DE-AUTHOR | Autorname in Bibitem: "I. M. Vinogradov" statt der in Paper 18 DE verwendeten deutschen Transliteration "I. M. Winogradow". Inkonsistenz innerhalb der DE-Paperserie. | Entweder konsequent "Winogradow" (deutsche Konvention) oder "Vinogradov" (internationale Konvention) verwenden. Empfehlung: "Winogradow" in allen DE-Papers, "Vinogradov" in allen EN-Papers. |
| 15 | Paper 20 | EN | GERING | 478–479 | BUG-B4-P20-EN-TABLE-N1000 | Für $n=1000$: Pred. $r_2(n) \approx 18$, Actual = 28. Die Abweichung (Faktor ~1.5) ist für $n=1000$ noch signifikant und könnte den Leser irreführen. Die Formel $r_2(n) \sim \mathfrak{S}(n) \cdot n / (\log n)^2$ mit variierendem $\mathfrak{S}$ macht die "Pred."-Spalte ohne expliziten $\mathfrak{S}(1000)$-Wert schwer nachprüfbar. | $\mathfrak{S}(1000)$-Wert explizit angeben oder Fußnote ergänzen: "Asymptotic predictions become reliable only for $n \gg 10^4$." |

---

## Mathematische Verifikation — Detailergebnisse

### Paper 13: Brunscher Satz (EN + DE)
- **Brunscher Satz:** Die Konvergenz $\sum_{p, p+2 \text{ prim}} 1/p < \infty$ ist korrekt formuliert. Die Beweisstruktur über das kombinatorische Sieb ist standard.
- **Sieb-Schranken:** Die Darstellung des Brun-Siebs mit Bonferroni-Ungleichungen ist korrekt.
- **Numerischer Wert:** Die Brunsche Konstante $B_2 \approx 1.902$ ist korrekt.
- **Urteil:** DRUCKREIF (EN + DE).

### Paper 14: Selbergs Sieb (EN + DE)
- **Lambda-Quadrat-Methode:** Korrekt dargestellt. Die Optimierung der $\lambda_d$-Gewichte über quadratische Formen ist standard.
- **Hauptergebnis:** $S(\mathcal{A}, z) \leq X \cdot V(z) + R$ ist korrekt.
- **V(z)/V_h(z)-Notation:** Konsistent verwendet nach früheren Fixes.
- **Brun-Titchmarsh:** $\pi(x; q, a) \leq 2x / (\phi(q) \log(x/q))$ korrekt formuliert.
- **Urteil:** DRUCKREIF (mit Stilanmerkung zu $|\lambda_d^*| \leq 1$).

### Paper 15: Große Sieb-Ungleichung (EN + DE)
- **Montgomery-Vaughan-Form:** $\sum_{q \leq Q} \sum_{a \pmod{q}}^* |S(a/q)|^2 \leq (N + Q^2 - 1) \sum |a_n|^2$ korrekt.
- **Fejér-Kern:** Korrekte Darstellung als Approximation der Identität.
- **Bombieri-Vinogradov:** $\sum_{q \leq Q} \max_{(a,q)=1} |\pi(x;q,a) - \text{li}(x)/\phi(q)| \ll x/(\log x)^A$ für $Q \leq x^{1/2}/(\log x)^B$ korrekt.
- **Urteil:** DRUCKREIF (Build-Nummer kosmetisch veraltet).

### Paper 16: Chens Satz (EN + DE)
- **Chen'sches Resultat:** $p + 2 = P_2$ (Primzahl oder Produkt zweier Primzahlen) korrekt formuliert.
- **Buchstab-Identität:** Korrekt dargestellt.
- **Switching-Prinzip:** Die Umkehrtechnik Type I/Type II-Summen ist korrekt skizziert.
- **Paritätsproblem:** Korrekt als fundamentale Grenze linearer Siebe beschrieben.
- **Urteil:** DRUCKREIF (Build-Nummer kosmetisch veraltet).

### Paper 17: Kreismethode (EN + DE)
- **Major/Minor Arcs:** Farey-Zerlegung korrekt.
- **Singuläre Reihe:** $\mathfrak{S}(n) = \sum_q c_q(n)/\phi(q)^s$ korrekt.
- **Singuläres Integral:** $\mathfrak{J}(n) = \int_{[0,1]^{s-1}} \dots$ korrekt.
- **Urteil:** DRUCKREIF (EN + DE).

### Paper 18: Vinogradov Drei-Primzahlen (EN + DE)
- **Exponentialsumme über Primzahlen:** $S(\alpha) = \sum_{n \leq N} \Lambda(n) e(n\alpha)$ korrekt.
- **Major-Arc-Asymptotik:** $S(\alpha) = \mu(q)/\phi(q) \cdot I(\beta) + O(N e^{-c\sqrt{\log N}})$ korrekt (Siegel-Walfisz).
- **Minor-Arc-Schranke (Vinogradov):** $|S(\alpha)| \ll N/(\log N)^A$ korrekt.
- **Helfgott 2013:** Korrekt als Vervollständigung für alle ungeraden $n \geq 7$ beschrieben.
- **1/6-Faktor:** FALSCH in EN, Zeile 136 (siehe Bug #7).
- **Urteil:** EN: ÜBERARBEITUNG (Bug #7). DE: DRUCKREIF (kein 1/6-Fehler vorhanden).

### Paper 19: Waringsches Problem (EN + DE)
- **$g(k)$-Formel:** $g(k) = 2^k + \lfloor(3/2)^k\rfloor - 2$ korrekt (für $k \geq 6$ bewiesen, für $k \leq 5$ einzeln bekannt).
- **$G(k)$-Schranken:** $G(2)=4$, $G(4)=15$ korrekt. $G(k) \leq k(\log k + 4.20)$ (Wooley) korrekt.
- **Vinogradov Mean Value Theorem:** $J_{s,k}(N) \ll N^{s+\varepsilon}$ für $s \leq k(k+1)/2$ korrekt (bewiesen BDG 2015 / Wooley 2016).
- **Urteil:** EN: DRUCKREIF. DE: ÜBERARBEITUNG (fehlende Inhalte, Bugs #9–11).

### Paper 20: Goldbach Singuläre Reihe (EN + DE)
- **Singuläre Reihe für Goldbach:** $\mathfrak{S}(n) = 2C_2 \prod_{p|n, p>2} (p-1)/(p-2)$ für gerade $n$ korrekt.
- **Zwillingsprimzahlkonstante:** $C_2 = \prod_{p>2} (1 - 1/(p-1)^2) \approx 0.6602$ korrekt.
- **Hardy-Littlewood Hypothesis H:** Korrekt als Verallgemeinerung formuliert.
- **Urteil:** EN + DE: DRUCKREIF (mit geringen Mängeln in Tabelle/Bibliografie).

---

## Gesamturteil Batch 3 & 4

| Paper | EN | DE |
|-------|----|----|
| 13 Brunscher Satz | DRUCKREIF | DRUCKREIF |
| 14 Selbergs Sieb | DRUCKREIF (Stil) | DRUCKREIF (Stil) |
| 15 Große Sieb-Ungleichung | DRUCKREIF (Build) | DRUCKREIF (Build) |
| 16 Chens Satz | DRUCKREIF (Build) | DRUCKREIF (Build) |
| 17 Kreismethode | DRUCKREIF | DRUCKREIF |
| 18 Vinogradov 3-Prim | ÜBERARBEITUNG | DRUCKREIF |
| 19 Waringsches Problem | DRUCKREIF | ÜBERARBEITUNG |
| 20 Goldbach Sing. Reihe | DRUCKREIF (gering) | DRUCKREIF (gering) |

**Offene Fehler gesamt:** 15 (0 KRITISCH, 1 HOCH, 4 MITTEL, 6 GERING, 4 KOSMETISCH)
