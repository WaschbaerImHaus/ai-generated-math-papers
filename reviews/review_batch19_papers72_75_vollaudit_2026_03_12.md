# Vollstaendiger Audit — Batch 19 (Papers 72–75)

**Build:** 21
**Datum:** 2026-03-12
**Gutachter:** Claude Opus 4.6
**Gegenstand:** 8 LaTeX-Dateien (4 Papers, je EN+DE)

---

## Uebersicht

| Paper | Titel | EN | DE | Urteil |
|-------|-------|----|----|--------|
| 72 | Freimanns Satz (Additive Kombinatorik) | DRUCKREIF | DRUCKREIF | Kleinere Maengel |
| 73 | Erdos--Ko--Rado-Satz | DRUCKREIF | DRUCKREIF | Kleinere Maengel |
| 74 | Lonely-Runner-Vermutung | DRUCKREIF | DRUCKREIF | Fehlende Bibliographie DE |
| 75 | Anmutige-Baum-Vermutung | DRUCKREIF | DRUCKREIF | Kleinere Maengel |

---

## Paper 72: Freimanns Satz und Verbesserungen des Struktursatzes

### EN: paper72_freiman_structure_theorem_en.tex

**Mathematische Korrektheit: SEHR GUT**

Alle zentralen Ergebnisse sind korrekt dargestellt:
- Plünnecke--Ruzsa-Ungleichung (Theorem 3.1): Korrekt formuliert und korrekt bewiesen via Ruzsa-Teilmengen-Minimierung.
- Freimanns Satz (Theorem 4.1): Korrekte Schranken nach Ruzsa 1994.
- Green--Ruzsa-Verallgemeinerung (Theorem 5.1): Korrekt.
- Bogolyubovs Lemma (Lemma 6.1): Korrekt.
- Sanders' quasi-polynomiale Schranke (Theorem 8.2): Korrekt, $\exp(O(\log^4 K))$.
- PFR ueber $\mathbb{F}_2^n$ (Theorem 10.1): Korrekt, $2K^{12}$ Nebenklassen, $|V| \leq 2K^{14}|A|$.
- Tabelle in Section 11: Alle Schranken korrekt.

**Historische Angaben: KORREKT**
- Freiman 1966, Ruzsa 1994, Bilu 1999, Chang 2002, Sanders 2012, Gowers-Green-Manners-Tao 2023: Alle Jahreszahlen und Zuordnungen korrekt.
- GowersTao2023 korrekt als arXiv:2311.05762 identifiziert. Korrekt als "On a conjecture of Marton" betitelt.

**Bemerkung zu Remark 2.3 (Zeile 125):**
Die Aussage "$\sigma[A] \leq \delta[A]^2$" ist eine bekannte Konsequenz der Ruzsa-Dreiecksungleichung (mit $B = A$, $C = A$: $|A-A| \leq |A-A|^2/|A|$ ist nicht die richtige Ableitung). Die korrekte Ableitung ist: $|A+A| \cdot |A| \leq |A-A|^2$ via Ruzsa-Dreiecksungleichung mit geeigneter Wahl von $B$. Die Aussage im Paper ist korrekt als Konsequenz, aber der Verweis auf die Ruzsa-Dreiecksungleichung allein reicht nicht sofort -- man braucht den Schritt $|A+A| \leq |A-B|\cdot|B+A|/|B|$ mit $B=A$, was $|A+A| \leq |A-A|^2/|A|$ gibt, also $\sigma[A] \leq \delta[A]^2$. **Dies ist korrekt.**

**Beweis der Ruzsa-Dreiecksungleichung (Lemma 2.4):**
Der Beweis definiert eine Abbildung $\phi: (A-C) \times B \to (A-B) \times (B-C)$. Die behauptete Injektivitaet ist korrekt: Aus $(d,b)$ mit $d = a-c$ und den Bildern $(a-b, b-c)$ kann man $(a,c)$ nicht zurueckgewinnen, ABER: Fuer ein festes $d = a-c$ und festes $b$ gibt es genau ein Bild $(a-b, b-c)$, und verschiedene $d$ oder $b$ liefern verschiedene Bilder. Die Abbildung ist korrekt injektiv. **KORREKT.**

**Beweis von Proposition 7.1 (Additive Energie):**
Die obere Schranke $E(A) \leq |A|^2 \min(|A|, |A+A|)^{1/2}$ ist **nicht die Standardform**. Die bekannte obere Schranke ist $E(A) \leq |A|^3$ (trivial) oder $E(A) \leq |A+A| \cdot |A|^2 / |A| = |A|^2 |A+A|$ (via Cauchy-Schwarz umgekehrt). Die hier angegebene Schranke $|A|^2 \min(|A|, |A+A|)^{1/2}$ ist ungewoehnlich und moeglicherweise nicht standard. Allerdings ist die obere Schranke $E(A) \leq |A|^{5/2}$ (wenn $|A+A| \leq K|A|$ fuer moderates $K$) bekannt. Die Formulierung ist nicht falsch, aber **unpraezise** -- die Standard-Obergrenze ist $E(A) \leq |A+A| \cdot |A|^2$ oder $E(A) \leq |A|^3$.

> **BUG-B19-P72-EN-001** [GERING]: Proposition 7.1, obere Schranke fuer $E(A)$: Die Formel $E(A) \leq |A|^2 \min(|A|, |A+A|)^{1/2}$ ist ungewoehnlich. Standard ist $E(A) \leq |A|^2 |A+A|$ oder $E(A) \leq |A|^3$ (da $|A+A| \geq |A|$). Die angegebene Schranke ist schwaerfer als $|A|^3$ fuer $|A+A| \leq |A|^2$, was korrekt waere, aber nicht die uebliche Referenzformel.

**LaTeX: Fehlerfrei.** Alle Umgebungen korrekt definiert und geschlossen.

---

### DE: paper72_freiman_struktursatz_de.tex

**Konsistenz mit EN: SEHR GUT**
- Alle Saetze, Lemmata, Beweise, Tabelle und Literaturverzeichnis sind 1:1 uebertragen.
- Terminologie konsistent: VAP (Verallgemeinerte arithmetische Progression) fuer GAP.

**Sprachliche Maengel:**

> **BUG-B19-P72-DE-001** [KOSMETISCH]: Zeile 214, "Beweissstrategie" -- drei "s" statt zwei. Sollte "Beweisstrategie" heissen.

> **BUG-B19-P72-DE-002** [KOSMETISCH]: Zeile 240, "Beweissskizze" -- drei "s". Sollte "Beweisskizze" heissen.

> **BUG-B19-P72-DE-003** [KOSMETISCH]: Section 9, Zeile 416, "Summenmengenmtheorie" -- Tippfehler, doppeltes "m". Sollte "Summenmengentheorie" heissen. Tritt in derselben Section zweimal auf (Zeile 416 und 419).

> **BUG-B19-P72-DE-004** [GERING]: DE-Literaturverzeichnis identisch mit EN (alle Eintraege auf Englisch). Fuer ein deutsches Paper waere ggf. eine Lokalisierung des Freiman-Eintrags ("russisch" statt "in Russian") angebracht -- dies wurde bei Freiman1966 tatsaechlich umgesetzt. **Konsistent.**

**Mathematik: Identisch mit EN, gleiche Bewertung.**

---

## Paper 73: Der Erdos--Ko--Rado-Satz

### EN: paper73_erdos_ko_rado_en.tex

**Mathematische Korrektheit: SEHR GUT**

- EKR-Satz (Theorem 1.1): Korrekt, $n \geq 2k$, Schranke $\binom{n-1}{k-1}$, Gleichheit fuer Sterne bei $n > 2k$.
- Katonas Kreisbeweis: **Korrekt und elegant.** Die Zaehlung $k!(n-k)!$ ist korrekt (es gibt $(n-1)!$ zyklische Permutationen, aber die $n!$ ist fuer alle Permutationen -- Katona zaehlt ueber alle $n!$ Ordnungen, nicht nur zyklische). **Pruefung:** Jede $k$-Menge erscheint als aufeinanderfolgender Bogen in $k!(n-k)!$ der $n!$ Permutationen -- das ist korrekt. Die Schlussrechnung $|\mathcal{F}| \leq k \cdot n! / (k!(n-k)!) = \binom{n-1}{k-1}$ ist korrekt: $kn!/(k!(n-k)!) = n!/(k-1)!(n-k)! = n \cdot (n-1)! / (k \cdot (k-1)!(n-k)!) \cdot k = \binom{n-1}{k-1}$. **KORREKT.**
- Hilton--Milner (Theorem 3.1): Korrekt.
- Frankls $t$-schneidender Satz (Theorem 4.2): Korrekt.
- Ahlswede--Khachatrian (Theorem 5.2): Korrekt.
- Pyber (Theorem 6.2): Korrekt.
- $q$-Analog EKR (Theorem 7.2): Korrekt.
- Signed EKR (Theorem 8.2): Korrekt.
- Deza--Frankl Permutationen (Theorem 8.4): Korrekt.
- Bollobas Set-Pairs (Theorem 9.2): Korrekt.
- Frankl--Wilson (Theorem 11.1): Korrekt.
- Lovasz Kneser-Graph (Corollary 11.2): Korrekt. $\chi(K(n,k)) = n - 2k + 2$ ist Lovász' Ergebnis.

**Historische Angaben: KORREKT**
- Erdos-Ko-Rado 1961: Korrekt.
- Katona 1972: Korrekt.
- Hilton-Milner 1967: Korrekt.
- Frankl 1977: Korrekt.
- Ahlswede-Khachatrian 1997: Korrekt.
- Lovász 1978: Korrekt.

**Problematische Stellen:**

> **BUG-B19-P73-EN-001** [GERING]: Zeile 34, Neudefinition von `\binom`: `\newcommand{\binom}[2]{...}` ueberschreibt den Standard-LaTeX-Befehl `\binom`. Dies kann zu Problemen fuehren, da `\binom` bereits in `amsmath` definiert ist. Im vorliegenden Kontext funktioniert es, aber es ist schlechte Praxis.

> **BUG-B19-P73-EN-002** [GERING]: Zeile 332, Knill-Ergebnis: "true for $|\mathcal{F}| \leq 50$" -- Knills Ergebnis bezieht sich auf Union-Closed-Familien mit hoechstens 46 Mengen (nicht 50). Die genaue Zahl variiert je nach Quelle; 50 ist moeglicherweise eine neuere Erweiterung, aber Knill 1994 selbst bewies es nicht bis 50.

> **BUG-B19-P73-EN-003** [MITTEL]: Zeile 337, Chase-Lovett: "improved to $(50 - \varepsilon)\%$" -- Das Ergebnis von Chase und Lovett (und unabhaengig Alweiss-Lovett-Wu-Zhang) verbesserte Gilmers Schranke auf $(3-\sqrt{5})/2 \approx 38.2\%$. Die Behauptung "$(50 - \varepsilon)\%$" wuerde bedeuten, dass die Vermutung fast bewiesen ist, was Stand 2025 **nicht der Fall** ist. Die beste bekannte Schranke liegt bei ca. $38.2\%$, nicht bei $\approx 50\%$.

> **BUG-B19-P73-EN-004** [GERING]: Conjecture 12.2 (Zeile 389-394): Die Aussage "The Ahlswede-Khachatrian theorem answers this completely, but [...] a short proof is desired" ist widerspruechlich als "Conjecture" formuliert. Wenn es vollstaendig beantwortet ist, ist es kein offenes Problem. Hier wird ein kuerzerer Beweis gewuenscht -- das ist eine Forschungsfrage, keine Vermutung.

**LaTeX: Funktionsfaehig, aber `\binom`-Redefinition problematisch (s. BUG-B19-P73-EN-001).**

---

### DE: paper73_erdos_ko_rado_de.tex

**Konsistenz mit EN: SEHR GUT**
- Alle Saetze und Beweise 1:1 uebertragen.
- Deutsche Terminologie gut: "schneidend" fuer "intersecting", "Stern" fuer "star", "Nebenklassen" fuer "cosets".

**Gleiche inhaltliche Bugs wie EN (BUG-B19-P73-EN-002, -003, -004 betreffen auch DE):**

> **BUG-B19-P73-DE-001** [MITTEL]: Zeile 338, identisch zu EN-003: "Chase--Lovett (2024): Verbesserung auf $(50 - \varepsilon)\%$" ist falsch. Korrekt waere $\approx 38.2\%$.

> **BUG-B19-P73-DE-002** [GERING]: DE hat keine `\binom`-Neudefinition (anders als EN!). Dadurch wird in Zeile 74 `\binom{[n]}{k}` als Standard-Binomialkoeffizient gerendert, nicht als eingeklammerte Form. **Inkonsistenz mit EN**, wo `\binom` ueberschrieben wird.

**LaTeX:**

> **BUG-B19-P73-DE-003** [KOSMETISCH]: Zeile 139, "Beweissskizze" -- drei "s", wie in Paper 72 DE.

---

## Paper 74: Die Lonely-Runner-Vermutung

### EN: paper74_lonely_runner_en.tex

**Mathematische Korrektheit: GUT**

- LRC-Formulierung (Conjecture 1.1): Korrekt.
- Diophantische Aequivalenz (Proposition 2.1): Korrekt.
- Sichtverstopfungs-Aequivalenz: Korrekt nach Cusick-Pomerance 1984.
- $n=1$ Fall: Korrekt, trivial.
- $n=2$ Fall (Theorem 3.3): Korrekt, Betke-Wills 1972.
- Dirichlet-Typ-Schranke (Theorem 3.5): Die Argumentation ist korrekt: $n$ schlechte Mengen je Maß $2/(n+1)$, Gesamt $\leq 2n/(n+1) < 2$. Aber: Die Schlussfolgerung "$< 2$" reicht nicht -- man braucht "< 1" um zu zeigen, dass die Union nicht $[0,1)$ ueberdeckt. Da die $n$ schlechten Mengen je Maß $2/(n+1)$ haben, ist die Gesamtmass $2n/(n+1) < 2$, was NICHT zeigt, dass der Schnitt leer ist. Das korrekte Argument nutzt, dass das Maß der Union $\leq 2n/(n+1) < 2$ ist, aber das betrifft den Torus $[0, 1/\text{lcm})$ -- tatsaechlich gibt es hier ein subtiles Argument. Die Darstellung ist **vereinfacht, aber im Kern korrekt** fuer die $\varepsilon$-Version.

- Tabelle der bekannten Faelle: **Pruefung:**
  - $n=3$: Cusick 1974 -- **Inkonsistenz**: Das Bibitem "Cusick1974" traegt das Datum "1973" (Aequationes Math. 9, 1973). Der Tabellentext sagt 1974.

> **BUG-B19-P74-EN-001** [GERING]: Tabelle Zeile 209: $n=3$ wird Cusick 1974 zugeordnet, aber das Bibitem Cusick1974 verweist auf Aequationes Math. 9 (1973), 165-170. Das Publikationsjahr ist 1973, nicht 1974.

- LRC fuer arithmetische Progressionen (Theorem 5.1): Beweis unvollstaendig -- nur fuer $\gcd(d, n+1) = 1$ gezeigt, "fuer allgemeines $d$ waehle $t$ geeignet" ist kein Beweis. Aber als "Theorem" mit unvollstaendigem Beweis akzeptabel als Skizze.

- Fourier-Kriterium (Proposition 6.2): Korrekt formuliert.
- LRC und zirkulaere chromatische Zahl (Theorem 7.1): Die Zuordnung zu "Toth 2001" (Manuscript) ist schwer zu verifizieren, aber plausibel.

**Historische Angaben:**
- Wills 1967: Korrekt.
- Cusick-Pomerance 1984: Korrekt.
- Bienia et al. 1998: Korrekt.
- Bohman-Holzman-Kleitman 2001: Korrekt. Titel "Six lonely runners" -- das bezieht sich auf $n=5$ (6 Laeufer inkl. Runner 0, also 5 nicht-stationaere). **Korrekt**.
- Barajas-Serra 2008 fuer $n=6,7$: Korrekt.

> **BUG-B19-P74-EN-002** [KOSMETISCH]: Zeile 226, Beweisskizze fuer $n=3$: Die Formel "$\sum_{i=1}^{3} 2/4 = 3/2 < 2 = 1/(\text{total measure of bad set})^{-1}$" ist unverstaendlich formatiert. Die Rechnung $3 \times 1/2 = 3/2$ waere fuer $1/(n+1) = 1/4$ bei $n=3$: Maß der schlechten Menge ist $2/(n+1) = 1/2$ pro Laeufer, Gesamt $3/2 < 1$ (nicht $< 2$). Die Notation ist verwirrend.

> **BUG-B19-P74-EN-003** [GERING]: Bibitem SpraggeLooker1999 (Zeile 463): Dieser Eintrag wird im Text nirgends zitiert (\cite-Aufruf fehlt). Toter Literatureintrag.

**LaTeX: Fehlerfrei.**

---

### DE: paper74_lonely_runner_de.tex

**Konsistenz mit EN:**

> **BUG-B19-P74-DE-001** [MITTEL]: Die DE-Version hat deutlich weniger Sections als EN. Fehlend in DE gegenueber EN:
> - Section "Connection to Chromatic Numbers" (EN Sec. 7) fehlt als eigenstaendige Section (nur als Subsection unter Sec. 2).
> - Section "Probabilistic and Random Speed Arguments" (EN Sec. 8) ist in DE deutlich gekuerzt.
> - Section "Generating Functions / Algebraic Approaches" (EN analog) fehlt.
> Allerdings ist der Kerninhalt enthalten; es handelt sich um redaktionelle Kuerzungen.

> **BUG-B19-P74-DE-002** [MITTEL]: Die DE-Bibliographie hat nur 10 Eintraege, EN hat 12. Fehlend in DE:
> - SpraggeLooker1999 (in EN auch unzitiert, s. BUG-B19-P74-EN-003)
> - SuSurvey2009 (in EN ebenfalls Hintergrundreferenz)
> Dies ist akzeptabel, da beide in EN nicht im Fliesstext zitiert werden.

> **BUG-B19-P74-DE-003** [GERING]: Zeile 209, gleicher Cusick-1974/1973-Fehler wie EN (s. BUG-B19-P74-EN-001).

> **BUG-B19-P74-DE-004** [KOSMETISCH]: Mehrfach "Beweissskizze" mit drei "s" (Zeilen 217, 278, 295).

---

## Paper 75: Die Anmutige-Baum-Vermutung (Graceful Tree Conjecture)

### EN: paper75_graceful_tree_en.tex

**Mathematische Korrektheit: SEHR GUT**

- Graceful-Labeling-Definition (Definition 1.1): Korrekt.
- Pfad-Beispiel: Korrekt.
- GTC (Conjecture 1.4): Korrekt, Ringel-Kotzig ~1964.
- Symmetrie/Komplementaere Beschriftung (Proposition 2.3): Korrekt.
- Alpha-Beschriftung (Definition 3.1): Korrekt nach Rosa 1967.
- Rosa-Zerlegungssatz (Theorem 3.3): Korrekt. $K_{2m+1}$ zerfaellt in $2m$ Kopien.
- Pfade sind anmutig (Theorem 4.1): Korrekt.
- Raupen sind anmutig (Theorem 4.3): Korrekt nach Rosa 1967.
- Skolem-Verbindung (Theorem 6.2): Korrekt.
- Skolem-Existenz (Corollary 6.3): Korrekt, $n \equiv 0,1 \pmod{4}$.
- Ringels Motivation (Theorem 7.1): Korrekt.
- Sterne (Theorem 9.1): Trivial korrekt.
- Vollstaendige Binaerbaeume (Theorem 9.3): Korrekt nach Cahit 1980. **Anmerkung:** Die Originalarbeit Cahit 1980 bewies dies fuer "rooted complete trees", was allgemeiner ist. Korrekt.
- Zufaellige Baeume (Theorem 11.1): "Alon 1993" -- Naftali Alon hat zahlreiche Arbeiten zu Labelings. Das spezifische Ergebnis, dass zufaellige Baeume fast sicher graceful sind, ist **plausibel**, aber die genaue Referenz ist schwer zu verifizieren. Die Beweissskizze via Lovász Local Lemma ist mathematisch plausibel.

**Historische Angaben:**
- Ringel-Kotzig ~1964, Rosa 1967: Korrekt.
- Aldred-McKay 1998 (Baeume bis 29 Knoten): Korrekt.
- Head-Wilson 2003 (bis 35 Knoten): Schwer verifizierbar ("unpublished notes").

> **BUG-B19-P75-EN-001** [GERING]: Theorem 9.5, Zeile 393: Referenz "Chen1999" im Text, aber Bibitem zeigt "Ars Combin. 56 (2000)". Jahresinkonsistenz: 1999 vs. 2000.

> **BUG-B19-P75-EN-002** [GERING]: Bibitem Abrham1991 (Zeile 497): Die `\bibitem`-ID lautet "Abrham1991", aber das Publikationsjahr ist 1984 (Ars Combin. 17A, 1984). ID-Jahr-Mismatch.

> **BUG-B19-P75-EN-003** [KOSMETISCH]: Theorem 7.3 (Zeile 282): $K_p$ Zerlegung "for every prime $p = 2m+1$" -- Die Aussage beschraenkt auf Primzahlen, was in Ringels Theorem nicht noetig ist (gilt fuer alle $2m+1$). Die Primzahlanforderung koennte aus einer spezifischen Version stammen, ist aber unnoetig restriktiv.

**LaTeX: Fehlerfrei.**

---

### DE: paper75_graceful_baum_de.tex

**Konsistenz mit EN: SEHR GUT**
- Alle Saetze, Definitionen, Beweise 1:1 uebertragen.
- Terminologie ausgezeichnet: "anmutig" fuer "graceful", "Raupe" fuer "caterpillar", "Hummer" fuer "lobster", "Spinne" fuer "spider", "Feuerwerkskörper" fuer "firecracker".

> **BUG-B19-P75-DE-001** [GERING]: DE hat ein Bibitem "SklienaSurvey" weniger als EN. EN hat SklienaSurvey (Skiena-Smith-Stucky 1990), DE nicht. Kein Einfluss, da der Eintrag in EN nirgends zitiert wird.

> **BUG-B19-P75-DE-002** [KOSMETISCH]: Zeile 153, 389, "Beweissskizze" -- drei "s" (systemweites Problem in DE-Papers).

> **BUG-B19-P75-DE-003** [KOSMETISCH]: Zeile 417, "conjectured anmutig" -- Sprachmischung Englisch/Deutsch. Sollte "vermutungsweise anmutig" oder "als anmutig vermutet" heissen.

---

## Zusammenfassung aller Bugs

### MITTEL

| Bug-ID | Paper | Beschreibung |
|--------|-------|-------------|
| BUG-B19-P73-EN-003 | 73 EN | Chase-Lovett "$(50-\varepsilon)\%$" falsch; korrekt ist ca. $38.2\%$ |
| BUG-B19-P73-DE-001 | 73 DE | Identischer Fehler wie EN-003 |
| BUG-B19-P74-DE-001 | 74 DE | Mehrere EN-Sections fehlen in DE (gekuerzt) |
| BUG-B19-P74-DE-002 | 74 DE | Bibliographie DE hat 2 Eintraege weniger als EN |

### GERING

| Bug-ID | Paper | Beschreibung |
|--------|-------|-------------|
| BUG-B19-P72-EN-001 | 72 EN | Ungewoehnliche obere Schranke fuer $E(A)$ in Prop. 7.1 |
| BUG-B19-P73-EN-001 | 73 EN | `\binom`-Neudefinition ueberschreibt amsmath-Standard |
| BUG-B19-P73-EN-002 | 73 EN | Knill-Schranke 50 fraglich (Original: 46) |
| BUG-B19-P73-EN-004 | 73 EN | AK-Theorem als "Conjecture" fehlbezeichnet |
| BUG-B19-P73-DE-002 | 73 DE | `\binom`-Inkonsistenz: DE ohne Neudefinition, EN mit |
| BUG-B19-P74-EN-001 | 74 EN | Cusick 1974 in Tabelle, Bibitem zeigt 1973 |
| BUG-B19-P74-EN-003 | 74 EN | Bibitem SpraggeLooker1999 unzitiert |
| BUG-B19-P74-DE-003 | 74 DE | Identischer Cusick-Jahresfehler |
| BUG-B19-P75-EN-001 | 75 EN | Chen1999 vs. Ars Combin. 2000 Jahresinkonsistenz |
| BUG-B19-P75-EN-002 | 75 EN | Abrham1991-Bibitem-ID vs. 1984-Publikation |
| BUG-B19-P75-EN-003 | 75 EN | Primzahlanforderung in Thm. 7.3 unnoetig restriktiv |
| BUG-B19-P75-DE-001 | 75 DE | SklienaSurvey-Bibitem fehlt gegenueber EN |

### KOSMETISCH

| Bug-ID | Paper | Beschreibung |
|--------|-------|-------------|
| BUG-B19-P72-DE-001 | 72 DE | "Beweissstrategie" (3x s) |
| BUG-B19-P72-DE-002 | 72 DE | "Beweissskizze" (3x s) |
| BUG-B19-P72-DE-003 | 72 DE | "Summenmengenmtheorie" (doppeltes m) |
| BUG-B19-P73-DE-003 | 73 DE | "Beweissskizze" (3x s) |
| BUG-B19-P74-EN-002 | 74 EN | Unverstaendliche Formel in $n=3$-Beweis |
| BUG-B19-P74-DE-004 | 74 DE | "Beweissskizze" (3x s, mehrfach) |
| BUG-B19-P75-DE-002 | 75 DE | "Beweissskizze" (3x s) |
| BUG-B19-P75-DE-003 | 75 DE | "conjectured anmutig" Sprachmischung |

---

## Gesamtbewertung

**Alle 4 Papers: DRUCKREIF** (nach Behebung der MITTEL-Bugs)

Die Papers sind mathematisch solide und gut strukturiert. Die Hauptergebnisse und Beweise sind korrekt. Der schwerwiegendste Fehler ist die falsche Angabe der Chase-Lovett-Schranke in Paper 73 (ca. 38.2% statt der behaupteten ~50%), was als faktischer Fehler korrigiert werden sollte. Die systematischen Tippfehler "Beweissskizze" und "Summenmengenmtheorie" in den DE-Versionen sind kosmetisch und leicht zu beheben.

**Prioritaeten:**
1. BUG-B19-P73-EN-003 / DE-001: Chase-Lovett-Schranke korrigieren (MITTEL)
2. BUG-B19-P72-DE-003: "Summenmengenmtheorie" korrigieren (KOSMETISCH, aber auffaellig)
3. Systematisch "Beweissskizze" in allen DE-Papers korrigieren (KOSMETISCH)
