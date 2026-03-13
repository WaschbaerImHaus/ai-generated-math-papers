# Vollstaendiger Audit: Batch 21 (Papers 80--83)

**Datum:** 2026-03-12
**Build:** 21
**Gutachter:** Claude Opus 4.6 (Mathematik-Peer-Review)

---

## Uebersicht

| Paper | Titel | EN | DE | Urteil |
|-------|-------|----|----|--------|
| 80 | Jacobian-Vermutung | OK | Maengel | UEBERARBEITUNG (DE) |
| 81 | Hadwigers Vermutung | Maengel | OK | UEBERARBEITUNG (EN) |
| 82 | Andrews--Curtis-Vermutung | OK | OK | DRUCKREIF |
| 83 | Beal-Vermutung | Maengel | OK | UEBERARBEITUNG (EN) |

---

## Paper 80: The Jacobian Conjecture / Die Jacobian-Vermutung

### EN-Version: paper80_jacobian_conjecture_en.tex

**Urteil: DRUCKREIF (mit geringfuegigen Anmerkungen)**

**Mathematische Korrektheit:**
- Theorem 2.1 (n=1): Korrekt und vollstaendig bewiesen.
- Theorem 4.1 (BCW): Korrekt dargestellt. Die Reduktion auf kubisch-lineare Abbildungen ist standard (Bass--Connell--Wright 1982).
- Theorem 5.1 (Yagzhev): Korrekt.
- Theorem 6.1 (Pinchuk): Korrekt. Pinchuk 1994 ist das Standardgegenbeispiel.
- Theorem 7.1 (Wang): Korrekt fuer Grad <= 2.
- Theorem 8.1 (JC=DC): Korrekt. Tsuchimoto 2005 und BK 2007 korrekt zitiert.
- Die stabile Aequivalenz JC_n => DC_n und DC_{2n} => JC_n ist korrekt dargestellt.

**Anmerkungen:**
- Zeile 96--98: Die Aussage "by Hadamard's theorem on holomorphic functions" ist etwas ungenau. Gemeint ist, dass eine lokale Biholomorphie mit konstanter Jacobi-Determinante global injektiv ist -- das folgt nicht direkt aus Hadamards Satz, sondern aus der Tatsache, dass det(JF) != 0 ueberall die lokale Invertierbarkeit garantiert (Umkehrsatz). Der globale Aspekt (Surjektivitaet, Injektivitaet) ist genau das, was die Vermutung behauptet. Kein Fehler, aber unpraesize.
- Zeile 321: "degree approximately 10^3" fuer Pinchuks Gegenbeispiel. Die tatsaechliche Konstruktion hat Grad 25 (Pinchuk 1994, Math. Z. 217). Der Grad 10^3 ist falsch.

#### BUG-B21-P80-EN-001 [MITTEL]
**Zeile 321:** Pinchuks Gegenbeispiel hat **nicht** Grad "approximately 10^3". Der tatsaechliche Grad der Abbildung ist 25 (Pinchuk, Math. Z. 217, 1994, S. 1--4). Der Kommentar "the precise degree depends on the substitution" ist irrefuehrend -- der Grad ist in der Literatur klar angegeben.

### DE-Version: paper80_jacobian_vermutung_de.tex

**Urteil: UEBERARBEITUNG**

**Mathematische Korrektheit:** Identisch zur EN-Version, korrekt.

**Konsistenz EN vs. DE:**

#### BUG-B21-P80-DE-001 [MITTEL]
**Fehlende Sections:** Die DE-Version hat **3 Sections weniger** als die EN-Version:
- EN Section 9 "Further Algebraic Equivalences" (Poisson conjecture, stable tameness) fehlt komplett in DE.
- EN Section 10 "Open Problems and Current Research" fehlt komplett in DE (nur kurze Aufzaehlung in Section 9 vorhanden, nicht so ausfuehrlich).
- EN Conjecture 10.1 "Generalized BCW" fehlt in DE.
Die DE-Version hat zwar "Offene Probleme" (Section 9), aber ohne die Formulierung der generalisierten BCW-Vermutung und ohne das Theorem ueber die Aequivalenz mit der Poisson-Vermutung.

#### BUG-B21-P80-DE-002 [GERING]
**Zeile 94:** Tippfehler "Kelllers" statt "Kellers" (dreifaches l).

#### BUG-B21-P80-DE-003 [KOSMETISCH]
Die DE-Version verwendet `definition_de` als Umgebungsname, was ungewoehnlich, aber funktional ist. In der EN-Version wird `definition` verwendet. Nicht kritisch, aber inkonsistent.

**Bibliographie:** EN hat 12 bibitems, DE hat 12 bibitems -- identisch und korrekt.

---

## Paper 81: Hadwiger's Conjecture / Hadwigers Vermutung

### EN-Version: paper81_hadwiger_conjecture_en.tex

**Urteil: UEBERARBEITUNG**

**Mathematische Korrektheit:**
- Conjecture 1.1 (Hadwiger): Korrekt formuliert.
- Theorem 3.1 (k <= 3): Korrekt. Beweis fuer k=3 via ungerade Kreise ist standard.
- Theorem 3.3 (k=4, Dirac): Korrekt.
- Theorem 4.1 (Wagner): **PROBLEM** -- siehe unten.
- Theorem 4.3 (4CT): Korrekt.
- Theorem 5.1 (RST, k=6): Korrekt.
- Theorem 6.1 (Kostochka/Thomason): Korrekt.
- Theorem 7.1 (Kuehn--Osthus): Korrekt.

#### BUG-B21-P81-EN-001 [HOCH]
**Zeile 153--158, Beweis von Theorem 3.3 (Dirac, k=4):** Der Beweis behauptet: "4-chromatic graphs cannot be planar (by the planar 5-color theorem)". Das ist **falsch**. 4-chromatische Graphen **koennen** planar sein -- zum Beispiel ist jeder planare Graph mit chromatischer Zahl genau 4 ein Gegenbeispiel (und es gibt viele, z.B. Graphen, die den 4-Farben-Satz benoetigen). Der 5-Farben-Satz besagt, dass planare Graphen 5-faerbbar sind, was hier irrelevant ist. Diracs Beweis fuer k=4 nutzt eine strukturelle Analyse mit Kempe-Ketten und die Tatsache, dass Graphen mit Minimalgrad >= 3 eine Unterteilung von K_4 enthalten muessen, NICHT eine Planaritaet-Argumentation.

#### BUG-B21-P81-EN-002 [MITTEL]
**Zeile 169:** Wagners Satz (Theorem 4.1) wird als "A graph G is planar if and only if it contains neither K_5 nor K_{3,3} as a minor" formuliert. Das ist zwar korrekt, aber es ist im Text als "Wagner, 1937" zitiert. Wagners 1937-Paper bewies primaer den Struktursatz fuer K_5-Minor-freie Graphen (Aufbau aus planaren Graphen und V_8 ueber 3-Summen). Die aequivalente Formulierung als "Minor-Version von Kuratowski" wird oft Wagner zugeschrieben, aber die explizite Aequivalenz K_5/K_{3,3}-Minor <=> Planaritaet ist ein Korollar aus Wagners Struktursatz + Kuratowskis Theorem. Die Remark in Zeile 173--177 versucht das klarzustellen, aber sagt "Wagner's theorem is slightly stronger since minor is a weaker condition than topological minor" -- das ist logisch inkorrekt. Dass Minor eine schwaechere Bedingung als topologischer Minor ist, bedeutet, dass **Kuratowskis** Theorem staerker ist (da es eine staerkere Bedingung ausschliesst), nicht Wagners. Die Implikation laeuft andersherum: "kein Subgraph isomorph zu Subdivision von K_5/K_{3,3}" (Kuratowski) => "kein K_5/K_{3,3}-Minor" (Wagner), aber nicht umgekehrt. Da beide aequivalent zu Planaritaet sind, sind sie gleich stark, nur in der Formulierung unterschiedlich.

#### BUG-B21-P81-EN-003 [MITTEL]
**Zeile 193--200, Beweis von Theorem 4.4 (Hadwiger k=5 <=> 4CT):** Der Beweis enthaelt einen abgebrochenen Gedankengang: "If G ⊇ K_{3,3}, then since chi(K_{3,3}) = 2... but wait, G is (5)-chromatic. We need a more careful argument." Dies ist fuer ein peer-reviewed Paper inakzeptabel -- es liest sich wie eine Entwurfsnotiz. Das korrekte Argument nutzt Wagners Struktursatz (K_5-Minor-freie Graphen sind 3-Summen aus planaren Graphen und V_8), was dann im Text auch steht, aber der "but wait"-Abschnitt muss entfernt werden.

#### BUG-B21-P81-EN-004 [GERING]
**Zeile 221--222:** Der RST-Beweis behauptet, minimale Gegenbeispiele seien "vertex-transitive (by minimality arguments)". Das ist **falsch**. Minimale Gegenbeispiele zu Hadwiger fuer k=6 sind 6-vertex-critical (jeder echte Teilgraph hat chromatische Zahl <= 5), aber NICHT vertex-transitiv. Vertex-Transitivitaet ist eine viel staerkere Eigenschaft (jeder Knoten wird durch einen Automorphismus auf jeden anderen abgebildet).

#### BUG-B21-P81-EN-005 [GERING]
**Zeile 366:** "Mader's bound gives 5-colorability (since 2^{4-2}+1 = 5)" -- hier wird k=5 eingesetzt: 2^{5-2}+1 = 9, nicht 5. Fuer k=5: 2^{k-2} = 2^3 = 8, also 8-degenerate, also 9-faerbbar. Der Wert 5 stimmt nicht. Korrekterweise liefert Mader fuer k=5 die Schranke 9 (was schwaecher ist als der 5-Farben-Satz fuer planare Graphen).

Korrektur: Fuer k=4 liefert Mader 2^{4-2}+1 = 5, und Hadwiger fuer k=4 sagt 4-faerbbar. Fuer k=5 liefert Mader 2^{5-2}+1 = 9. Der Text verwechselt offenbar k=4 und k=5.

### DE-Version: paper81_hadwiger_vermutung_de.tex

**Urteil: DRUCKREIF (abhaengig von EN-Fixes)**

**Konsistenz EN vs. DE:**
- Die DE-Version enthaelt die problematischen Passagen aus EN NICHT: kein falscher Dirac-Beweis (nur Beweisidee ohne die falsche Planaritaets-Behauptung), kein "but wait"-Abschnitt, keine Vertex-Transitivitaets-Behauptung.
- Die DE-Version ist insgesamt etwas kuerzer und vermeidet dadurch mehrere Fehler der EN-Version.
- Die DE-Version hat KEINE Section "Odd Hadwiger Conjecture" (EN Section 9), was ein inhaltliches Defizit, aber kein Fehler ist.

#### BUG-B21-P81-DE-001 [GERING]
**Fehlende Section:** Die "Ungerade Hadwiger-Vermutung" (Gerards--Seymour) fehlt in DE, ist aber in EN (Section 9) vorhanden. Inhaltliches Defizit.

#### BUG-B21-P81-DE-002 [KOSMETISCH]
**Zeile 323:** "fundamentell" statt "fundamental" -- Tippfehler (korrekt: "fundamental" oder "grundlegend").

**Bibliographie:** Beide Versionen identisch, korrekt.

---

## Paper 82: Andrews--Curtis Conjecture / Andrews--Curtis-Vermutung

### EN-Version: paper82_andrews_curtis_en.tex

**Urteil: DRUCKREIF**

**Mathematische Korrektheit:**
- Conjecture 1.1 (AC): Korrekt formuliert mit den 5 AC-Zuegen.
- AK(n)-Praesentation und Beweis der Trivialitaet: Korrekt.
- Tietze-Transformationen: Korrekt dargestellt.
- h-Kobordismussatz (Smale 1960): Korrekt.
- s-Kobordismussatz (Barden--Mazur--Stallings 1964): Korrekt.
- Boone--Novikov-Satz: Korrekt.
- Markov-Satz: Korrekt.
- Deficiency-Ergebnisse: Korrekt.

**Historische Angaben:**
- Andrews--Curtis 1965: Korrekt (Proc. Amer. Math. Soc. 16).
- Akbulut--Kirby 1979: Korrekt.
- Smale "1960" vs. Bibitem "1961" -- der Satz heisst traditionell "Smale 1960" (Ankuendigung), Publikation 1961.

#### BUG-B21-P82-EN-001 [GERING]
**Zeile 433 (Bibitem MS99):** Die Bibliographie-Referenz "MS99" ist falsch datiert. Der Titel im Bibitem ist "The geometry of finitely presented infinite simple groups" in MSRI Publ. 23 (1992), S. 1--30. Aber der Schluessel heisst "MS99" (1999) und der Autor ist "C. F. Miller, P. Schupp". Tatsaechlich: Miller und Schupp publizierten zur AC-Vermutung in verschiedenen Arbeiten, aber die hier zitierte Arbeit aus MSRI 1992 ist von Miller allein oder anders betitelt. Die Referenz ist inkonsistent (Schluessel suggeriert 1999, Inhalt sagt 1992).

### DE-Version: paper82_andrews_curtis_de.tex

**Urteil: DRUCKREIF**

**Konsistenz EN vs. DE:**
- DE-Version ist vollstaendig, alle Sections vorhanden.
- EN hat 10 Bibitems, DE hat 8 Bibitems -- es fehlen in DE: `\bibitem{Berge}` (Berge 2000, unpublished) und `\bibitem{Stallings65}` (Stallings 1965). Kein inhaltlicher Verlust, da diese nicht zitiert werden in DE.
- Inhaltlich korrekt und konsistent.

#### BUG-B21-P82-DE-001 [GERING]
**Zeile 413 (Bibitem LS77):** Der Titel wird als "Kombinatorische Gruppentheorie" (Deutsch) angegeben, aber das Originalbuch von Lyndon/Schupp (Springer 1977) heisst "Combinatorial Group Theory" (Englisch). Die deutsche Ausgabe erschien erst spaeter. Geringfuegig, aber inkorrekt fuer die Erstausgabe 1977.

---

## Paper 83: The Beal Conjecture / Die Beal-Vermutung

### EN-Version: paper83_beal_conjecture_en.tex

**Urteil: UEBERARBEITUNG**

**Mathematische Korrektheit:**
- Conjecture 1.1 (Beal): Korrekt formuliert.
- Theorem 2.1 (Fermat--Wiles): Korrekt.
- Theorem 4.1 (Darmon--Granville): Korrekt und wichtig.
- abc => Beal (Theorem 6.1): Beweis korrekt.
- Theorem 7.1 (PSS): Korrekt.
- Theorem 9.2 (Beukers): Korrekt.
- Theorem 9.3 (Euler-Charakteristik-Klassifikation): Korrekt.

#### BUG-B21-P83-EN-001 [MITTEL]
**Zeile 133--153 (Example 3.1):** Das Beispiel ist chaotisch und enthaelt Arbeitsnotizen: "Wait, this needs checking" (Zeile 138) und inkonsequente Rechnungen mit Fragezeichen. Fuer ein peer-reviewed Paper voellig inakzeptabel. Konkrete Fehler:
- Zeile 136--137: Die Formel "(2^m)^3 + (2^n)^3 = (2^{m-1} + 2^{n-1})^3?" ist falsch fuer allgemeine m,n (es gilt (a+b)^3 != a^3 + b^3 im Allgemeinen).
- Zeile 147: "2^6 + 2^6 = 2^7" wird als "(2^2)^3 + (2^2)^3 = 2^7" interpretiert, also x=y=3, z=7, mit gcd=4. Das ist korrekt: 64+64=128=2^7, und gcd(4,4,2)=2.
- Zeile 151--152: "(3^n)^3(1+8) = 9 * 3^{3n}" -- falsch, es waere (3^n)^3 * 9 = 9 * 3^{3n} = 3^{3n+2}, was keine Potenz eines einzelnen Werts ergibt. Unklar, was gezeigt werden soll.

#### BUG-B21-P83-EN-002 [MITTEL]
**Zeile 196--200, Beweis von Theorem 4.1 (Darmon--Granville):** Die Geschlechtsformel ist falsch. Es wird angegeben:
```
g = 1 + (pqr/2) * (1/p + 1/q + 1/r - 1) * (n-1)/n
```
fuer "suitable n". Diese Formel ist nicht die Standard-Geschlechtsformel fuer die Fermat-Kurve x^p + y^q = z^r. Die korrekte Formel haengt von der Ramifizierung ab und laeuft ueber Riemann--Hurwitz. Die Variable "n" ist undefiniert und taucht in der Formel unvermittelt auf. Ausserdem: wenn 1/p + 1/q + 1/r - 1 < 0 (wie im hyperbolischen Fall), waere der Faktor in Klammern negativ, und mit pqr/2 > 0 ergibt sich g < 1, was dem Anspruch g >= 2 widerspricht. Der Beweis sagt sogar selbst: "When 1/p + 1/q + 1/r < 1, the expression inside is negative, giving g >= 2" -- das ist ein direkter Widerspruch in sich (negativ => g >= 2 folgt NICHT aus der angegebenen Formel).

#### BUG-B21-P83-EN-003 [GERING]
**Zeile 237--243 (Tabelle):** Die Tabelle enthaelt "14^2??" mit doppeltem Fragezeichen (Zeile 237). Arbeitsentwurfsnotiz, nicht publikationsreif. Ausserdem: 6^3 + 10^3 = 216 + 1000 = 1216, und 1216 ist NICHT 14^2 = 196. Der Eintrag ist schlicht falsch.

#### BUG-B21-P83-EN-004 [GERING]
**Zeile 303--305 (Theorem 7.1 PSS):** Der Satz besagt "27 primitive solutions (up to signs)". Poonen--Schaefer--Stoll finden tatsaechlich 16 Familien primitiver Loesungen von x^2 + y^3 = z^7 (PSS, Duke Math. J. 137, 2007). Die Zahl 27 bezieht sich moeglicherweise auf x^2 + y^3 + z^5 = 0 (Beukers, Theorem 9.2). Hier werden zwei verschiedene Gleichungen verwechselt -- PSS behandelt x^2 + y^3 = z^7, Beukers behandelt x^2 + y^3 + z^5 = 0.

**Historische Angaben:**
- Beal 1993: Korrekt.
- Wiles 1995: Korrekt.
- Darmon--Granville 1995: Korrekt.
- Ribet Bibitem: Schluessel "Ribet04" aber Jahreszahl 1997 im Eintrag -- inkonsistent (BUG unten).

#### BUG-B21-P83-EN-005 [KOSMETISCH]
**Bibitem Ribet04:** Schluessel suggeriert 2004, aber der Eintrag ist Acta Arith. 79 (1997). Korrekt waere Schluessel "Ribet97".

### DE-Version: paper83_beal_vermutung_de.tex

**Urteil: DRUCKREIF (abhaengig von EN-Fixes)**

**Konsistenz EN vs. DE:**
- Die DE-Version vermeidet die chaotischen Passagen der EN-Version: Kein "Wait, this needs checking", keine falschen Tabelleneintraege, keine falschen Geschlechtsformeln.
- Die DE-Geschlechtsformel in Theorem 4.1 (Darmon--Granville) ist korrekt vereinfacht ("Geschlecht g >= 2, wenn 1/p + 1/q + 1/r < 1"), ohne die falsche explizite Formel.
- DE hat weniger Bibitems als EN (9 vs. 10): Es fehlt `\bibitem{Ribet04}` (Ribet 1997). Kein Verlust, da der zugehoerige Satz (Ribet--Wiles) trotzdem in DE vorhanden ist, nur ohne Zitat.

#### BUG-B21-P83-DE-001 [GERING]
**Zeile 97:** "Computerssuche" -- Tippfehler, korrekt: "Computersuche" (doppeltes s).

**Bibliographie:** EN hat 10 bibitems, DE hat 9 -- fehlendes Ribet-Bibitem in DE. Geringfuegig.

---

## Zusammenfassung aller Bugs

| Bug-ID | Paper | Sprache | Schweregrad | Beschreibung |
|--------|-------|---------|-------------|--------------|
| BUG-B21-P80-EN-001 | 80 | EN | MITTEL | Pinchuks Grad falsch: 10^3 statt 25 |
| BUG-B21-P80-DE-001 | 80 | DE | MITTEL | 3 Sections fehlen gegenueber EN (Poisson, BCW-Gen., Open Problems) |
| BUG-B21-P80-DE-002 | 80 | DE | GERING | Tippfehler "Kelllers" statt "Kellers" |
| BUG-B21-P80-DE-003 | 80 | DE | KOSMETISCH | Inkonsistenter Umgebungsname definition_de |
| BUG-B21-P81-EN-001 | 81 | EN | HOCH | Falscher Beweis: 4-chromatische Graphen koennen planar sein |
| BUG-B21-P81-EN-002 | 81 | EN | MITTEL | Wagners Satz vs. Kuratowski: Logik der Staerke-Behauptung falsch |
| BUG-B21-P81-EN-003 | 81 | EN | MITTEL | "but wait" Entwurfsnotiz im Beweis von k=5 |
| BUG-B21-P81-EN-004 | 81 | EN | GERING | Minimale Gegenbeispiele sind nicht vertex-transitiv |
| BUG-B21-P81-EN-005 | 81 | EN | GERING | Maders Schranke fuer k=5 ist 9, nicht 5 (k=4/k=5 verwechselt) |
| BUG-B21-P81-DE-001 | 81 | DE | GERING | Fehlende Section "Ungerade Hadwiger-Vermutung" |
| BUG-B21-P81-DE-002 | 81 | DE | KOSMETISCH | Tippfehler "fundamentell" statt "fundamental" |
| BUG-B21-P82-EN-001 | 82 | EN | GERING | Bibitem MS99: Schluessel/Jahreszahl inkonsistent (1999 vs. 1992) |
| BUG-B21-P82-DE-001 | 82 | DE | GERING | Buchtitel "Kombinatorische Gruppentheorie" vs. engl. Original 1977 |
| BUG-B21-P83-EN-001 | 83 | EN | MITTEL | Chaotische Arbeitsnotizen ("Wait, this needs checking", Fragezeichen) |
| BUG-B21-P83-EN-002 | 83 | EN | MITTEL | Geschlechtsformel falsch und selbstwiderspruchlich |
| BUG-B21-P83-EN-003 | 83 | EN | GERING | Tabelle: 6^3+10^3 != 14^2, Fragezeichen im Text |
| BUG-B21-P83-EN-004 | 83 | EN | GERING | PSS hat nicht 27 Loesungen fuer x^2+y^3=z^7 (Verwechslung mit Beukers) |
| BUG-B21-P83-EN-005 | 83 | EN | KOSMETISCH | Bibitem-Schluessel Ribet04 vs. Jahreszahl 1997 |
| BUG-B21-P83-DE-001 | 83 | DE | GERING | Tippfehler "Computerssuche" |

---

## Statistik

- **KRITISCH:** 0
- **HOCH:** 1 (BUG-B21-P81-EN-001)
- **MITTEL:** 6
- **GERING:** 9
- **KOSMETISCH:** 3
- **Gesamt:** 19 Bugs

---

## Empfehlungen

1. **Paper 81 EN (Hadwiger) -- Prioritaet HOCH:** Der falsche Beweis fuer k=4 (4-chromatische Graphen seien nicht planar) muss dringend korrigiert werden. Diracs Originalargument nutzt Kempe-Ketten und die Existenz von K_4-Unterteilungen in Graphen mit Minimalgrad >= 3, NICHT ein Planaritaetsargument. Ausserdem die "but wait"-Entwurfsnotiz entfernen und die Vertex-Transitivitaets-Behauptung korrigieren.

2. **Paper 83 EN (Beal) -- Prioritaet MITTEL:** Die Arbeitsnotizen ("Wait, this needs checking", "??") muessen entfernt werden. Die falsche Geschlechtsformel korrigieren oder entfernen. Die PSS-Loesungszahl pruefen (16 Familien, nicht 27). Die Tabelle in Example 5.1 korrigieren (6^3+10^3 != 14^2).

3. **Paper 80 DE (Jacobian) -- Prioritaet MITTEL:** Die drei fehlenden Sections (Poisson-Aequivalenz, generalisierte BCW-Vermutung, detaillierte offene Probleme) sollten ergaenzt werden. Tippfehler "Kelllers" beheben.

4. **Paper 82 (Andrews--Curtis) EN+DE:** Geringfuegige Bibliographie-Korrekturen empfohlen, aber inhaltlich druckreif.
