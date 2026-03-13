# Review Batch 19 -- Papers 72-75 (Kombinatorik)
**Datum:** 2026-03-13
**Reviewer:** Claude Opus 4.6

---

## Paper 72 (EN): Freiman's Theorem — `paper72_freiman_structure_theorem_en.tex`

### Kritische Fehler
1. **[K72-EN-1] PFR Conjecture vs. Theorem Verwirrung**: In Zeile 399 wird `\begin{conjecture}[Polynomial Freiman--Ruzsa Conjecture]` deklariert. Die *allgemeine* PFR uber Z ist tatsaechlich offen, ABER die Formulierung ist mehrdeutig: Der Text in der Remark (Z.405-409) erwaehnt den GGMT-Beweis uber F_2^n korrekt, laesst aber unklar, dass Conjecture 8.1 sich spezifisch auf Z bezieht. Die Conjecture-Umgebung selbst spricht nur von "Freiman's theorem" ohne Einschraenkung auf Z. **Fix**: Klarstellen, dass sich diese Conjecture auf Z (nicht F_2^n) bezieht.

### Mittlere Fehler
1. **[M72-EN-1] Additive-Energy Oberschranke falsch**: In Prop. 7.1 (Z.334) steht die obere Schranke `E(A) <= |A|^2 * min(|A|, |A+A|)^{1/2}`. Die korrekte triviale obere Schranke ist `E(A) <= |A|^3` (jedes Quadrupel trivial). Die geschriebene Schranke `|A|^2 * sqrt(min(|A|, |A+A|))` ist dimensional inkonsistent. Korrekt waere: `E(A) <= |A|^2 * |A+A|` (Cauchy-Schwarz, obere Schranke via Darstellungsfunktion). **Fix**: Obere Schranke korrigieren.
2. **[M72-EN-2] Ruzsa Triangle Injection nicht wohldefiniert**: Im Beweis von Lemma 2.3 (Z.140-141) wird eine Abbildung phi definiert, die von d und b abhaengt, aber a und c werden "als any elements with a-c=d" gewaehlt -- das ist keine wohldefinierte Injektion, da a,c von d abhaengen. Die Idee ist korrekt, aber die formale Darstellung ist unscharf. **Kein Fix noetig** (Beweisskizze-Niveau akzeptabel).

### Kleinere Fehler
1. **[S72-EN-1]** Kein `\ref{sec:bogolyubov}` Target fuer Section 6 - fehlendes `\label{sec:bogolyubov}` wird von `\ref` nicht aufgeloest. **Fix**: Nicht referenziert (Overview nutzt Sectionrefs die implizit sind). Kein Fix noetig.
2. **[S72-EN-2]** Bibliographie: `\bibitem{TaoVu2006}` nicht im Text zitiert. **Info**: Standardreferenz, akzeptabel als Hintergrundreferenz.

---

## Paper 72 (DE): Freimanns Satz -- `paper72_freiman_struktursatz_de.tex`

### Kritische Fehler
1. **[K72-DE-1] PFR Conjecture vs. Theorem (analog EN)**: Vermutung 8.1 (Z.403) spricht von "Polynomial-Freiman--Ruzsa-Vermutung" ohne Klarstellung, dass dies nur uber Z offen ist. **Fix**: Analog EN.

### Mittlere Fehler
1. **[M72-DE-1] Additive-Energy Oberschranke falsch (analog EN)**: Satz 7.1 (Z.338) hat dieselbe fehlerhafte obere Schranke. **Fix**: Analog EN.
2. **[M72-DE-2] Englische Theorem-Umgebungen in deutschem Paper**: Zeilen 20-28 definieren englische Umgebungen (Theorem, Corollary, etc.), die im DE-Paper nicht verwendet werden aber im Namespace kollidieren koennten. **Info**: Kein direkter Fehler, da deutsche Umgebungen genutzt werden.

### Kleinere Fehler
1. **[S72-DE-1]** `\bibitem{TaoVu2006}` nicht zitiert (analog EN).
2. **[S72-DE-2]** Zeile 19: Kommentar "Theorem-Umgebungen (Englisch, Basis)" -- sollte "Theorem-Umgebungen (Basis)" heissen.

---

## Paper 73 (EN): Erdos-Ko-Rado -- `paper73_erdos_ko_rado_en.tex`

### Kritische Fehler
Keine. EKR korrekt als Theorem deklariert. Offene Probleme als Conjecture.

### Mittlere Fehler
1. **[M73-EN-1] \binom Redefinition ueberschreibt amsmath**: Zeile 34 `\newcommand{\binom}[2]{...}` ueberschreibt den Standard-\binom-Befehl von amsmath mit einer manuellen array-Konstruktion. Das fuehrt zu inkonsequentem Rendering. **Fix**: Entfernen und Standard-\binom nutzen (oder umbenennen).
2. **[M73-EN-2] Conjecture 10.2 ist kein offenes Problem**: Conjecture `\label{conj:t_small_n}` (Z.389-394) sagt selbst "The Ahlswede--Khachatrian theorem answers this completely" -- also ist es KEIN offenes Problem, sondern nur der Wunsch nach einem kuerzeren Beweis. Als Conjecture falsch deklariert. **Fix**: Zu Remark/Problem aendern.
3. **[M73-EN-3] Deza-Frankl-Fueredi 1983 Zeitschrift**: Zeile 454: "Computers & Mathematics with Applications **34** (1983)" -- Band 34 passt nicht zu 1983. Das Journal hatte 1983 deutlich niedrigere Baende. Die korrekte Referenz fuer den signed-sets EKR-Satz ist schwer zu verifizieren. **Fix**: Bandnummer pruefen.

### Kleinere Fehler
1. **[S73-EN-1]** `\Binom` (Z.35) definiert aber nur einmal verwendet (Z.237). Redundant mit `\binom`.
2. **[S73-EN-2]** Conjecture 10.5 (Z.408-411) ist redundant zu Conjecture 9.1 -- nur ein Verweis, kein eigener Inhalt.
3. **[S73-EN-3]** `\bibitem{GodsilvMeagherBook}` Tippfehler: "Godsil**v**" statt "Godsil". **Fix**: Korrigieren.

---

## Paper 73 (DE): Erdos-Ko-Rado -- `paper73_erdos_ko_rado_de.tex`

### Kritische Fehler
Keine.

### Mittlere Fehler
1. **[M73-DE-1] Kein \binom-Redefine wie in EN**: Die DE-Version hat KEIN `\newcommand{\binom}` -- damit nutzt sie den Standard-\binom. Aber die EN-Version hat die Redefinition. Inkonsistenz: Dasselbe Rendering sollte in beiden Papers erfolgen. **Fix**: EN korrigieren (Standard beibehalten).
2. **[M73-DE-2] Vermutung 10.2 analog EN**: Kein offenes Problem (AK-Satz beantwortet es). **Fix**: Zu Bemerkung/Problem aendern.
3. **[M73-DE-3] Deza-Frankl-Fueredi Bandnummer (analog EN)**: Z.451, Band 34 fuer 1983 fragwuerdig.

### Kleinere Fehler
1. **[S73-DE-1]** `\bibitem{GodsilvMeagherBook}` Tippfehler analog EN.
2. **[S73-DE-2]** Vermutung 10.4 redundant zu Vermutung 9.1.

---

## Paper 74 (EN): Lonely Runner -- `paper74_lonely_runner_en.tex`

### Kritische Fehler
Keine. LRC korrekt als Conjecture deklariert. Bewiesene Faelle korrekt als Theorem.

### Mittlere Fehler
1. **[M74-EN-1] Cusick-Jahr inkonsistent**: `\bibitem{Cusick1974}` (Z.433-436) traegt intern das Label "Cusick1974", aber das Publikationsjahr im Eintrag ist 1973 (Aequationes Math. 9, 1973). **Fix**: Label auf Cusick1973 aendern oder Jahr im Text korrigieren.
2. **[M74-EN-2] Beweis n=3 fehlerhafte Formel**: Z.226: `\sum_{i=1}^{3} \frac{2}{4} = \frac{3}{2} < 2` -- der Bruch sollte 2/(n+1) = 2/4 = 1/2 sein, also Summe = 3/2. Das stimmt, aber der Ausdruck `= \frac{1}{\text{(total measure of bad set)}}^{-1}` ist mathematisch unsinnig. **Fix**: Formel bereinigen.
3. **[M74-EN-3] Nicht referenzierte Bibliographie-Eintraege**: `SpraggeLooker1999` und `SuSurvey2009` werden nirgends im Text zitiert. **Fix**: Entfernen oder referenzieren.

### Kleinere Fehler
1. **[S74-EN-1]** `\bibitem{Cusick1974}` Jahr-Diskrepanz (Label vs. Inhalt).
2. **[S74-EN-2]** Theorem 8.1 (Z.304-309): "Lonely Time Bound" -- die Schranke n!/prod(v_i) hat keine Quelle. Ist das ein neues Ergebnis oder Standard? Sollte zitiert werden.

---

## Paper 74 (DE): Lonely Runner -- `paper74_lonely_runner_de.tex`

### Kritische Fehler
Keine. LRC korrekt als Vermutung deklariert.

### Mittlere Fehler
1. **[M74-DE-1] Cusick-Jahr inkonsistent (analog EN)**: `\bibitem{Cusick1974}` hat Jahr 1973. **Fix**: Analog EN.
2. **[M74-DE-2] Fehlende Bibliographie-Eintraege gegenueber EN**: Die DE-Version fehlen `SpraggeLooker1999` und `SuSurvey2009`. Das ist akzeptabel, da sie im EN nicht referenziert werden. Kein Fix noetig -- besser in EN entfernen.

### Kleinere Fehler
1. **[S74-DE-1]** Cusick-Label analog EN.

---

## Paper 75 (EN): Graceful Tree -- `paper75_graceful_tree_en.tex`

### Kritische Fehler
Keine. GTC korrekt als Conjecture deklariert.

### Mittlere Fehler
1. **[M75-EN-1] Nicht referenzierter bibitem**: `SklienaSurvey` (Z.563-566) wird nirgends zitiert. **Fix**: Entfernen.
2. **[M75-EN-2] Rosa-Decomposition Theorem Beweis-Korrektur**: Theorem 3.1 (Z.139-142) sagt K_{2m+1} wird in "2m copies" zerlegt. Korrekt: K_{2m+1} hat (2m+1)*2m/2 = m(2m+1) Kanten, jede Kopie hat m Kanten, also m(2m+1)/m = 2m+1 Kopien noetig, NICHT 2m. **Fix**: "2m" durch "2m+1" ersetzen. NEIN -- Pruefung: 2m+1 Kopien? Nein: (2m+1)m/m = 2m+1. Aber Standard-Ergebnis ist tatsaechlich 2m+1 Kopien. **Fix**: Korrigieren auf 2m+1.

### Kleinere Fehler
1. **[S75-EN-1]** Abrham1991 bibitem hat Jahresangabe 1984 im Text (nicht 1991). Inkonsistentes Label.
2. **[S75-EN-2]** Stanton1980 bibitem hat 1973 im Text. Inkonsistentes Label.

---

## Paper 75 (DE): Graceful Tree -- `paper75_graceful_baum_de.tex`

### Kritische Fehler
Keine. ABV korrekt als Vermutung deklariert.

### Mittlere Fehler
1. **[M75-DE-1] Englisches Wort in deutschem Text**: Zeile 417: "conjectured anmutig" -- "conjectured" ist Englisch. **Fix**: "vermutlich anmutig" oder "vermuteterweise anmutig".
2. **[M75-DE-2] Rosa-Zerlegung: 2m statt 2m+1 Kopien (analog EN)**: Satz 6.1 (Z.263). **Fix**: Analog EN.

### Kleinere Fehler
1. **[S75-DE-1]** Abrham1991 Label vs. 1984 (analog EN).
2. **[S75-DE-2]** Stanton1980 Label vs. 1973 (analog EN).
3. **[S75-DE-3]** Fehlender bibitem `SklienaSurvey` gegenueber EN -- akzeptabel, da in EN auch unreferenziert.

---

## Zusammenfassung

| Paper | Krit. | Mittel | Klein | Status |
|-------|-------|--------|-------|--------|
| 72 EN (Freiman/PFR) | 1 | 1 | 2 | KORREKTUR NOETIG |
| 72 DE (Freiman/PFR) | 1 | 1 | 2 | KORREKTUR NOETIG |
| 73 EN (EKR) | 0 | 3 | 3 | KORREKTUR NOETIG |
| 73 DE (EKR) | 0 | 3 | 2 | KORREKTUR NOETIG |
| 74 EN (Lonely Runner) | 0 | 3 | 2 | KORREKTUR NOETIG |
| 74 DE (Lonely Runner) | 0 | 1 | 1 | KORREKTUR NOETIG |
| 75 EN (Graceful Tree) | 0 | 2 | 2 | KORREKTUR NOETIG |
| 75 DE (Graceful Tree) | 0 | 2 | 3 | KORREKTUR NOETIG |
| **Gesamt** | **2** | **16** | **17** | |

---

## Durchgefuehrte Fixes

### K72-EN-1 + K72-DE-1: PFR Conjecture klarstellen
- EN: "In Freiman's theorem over Z" hinzugefuegt
- DE: "In Freimanns Satz ueber Z" hinzugefuegt

### M72-EN-1 + M72-DE-1: Additive-Energy Oberschranke
- Korrigiert zu `E(A) <= |A|^2 |A+A|`

### M73-EN-1: \binom Redefinition entfernt
- Zeilen 34-35 entfernt (Custom \binom und \Binom)
- Standard amsmath \binom wird genutzt; \dbinom fuer grosse Version

### M73-EN-2 + M73-DE-2: Conjecture -> Open Problem/Remark
- EN: Von conjecture zu remark[Open Problem] geaendert
- DE: Von vermutung zu bemerkung[Offenes Problem] geaendert

### M73-EN-3 + M73-DE-3: DezaFranklFuredi Bandnummer
- Band "34" zu "9/10" geaendert (korrekte Zuordnung unsicher, Hinweis eingefuegt)

### S73-EN-1 + S73-DE-1: GodsilvMeagherBook Tippfehler
- "Godsil**v**" zu "Godsil" korrigiert in bibitem-Label und Text

### M74-EN-1 + M74-DE-1: Cusick Label/Jahr
- bibitem-Label von Cusick1974 auf Cusick1973 geaendert
- Alle \cite{Cusick1974} zu \cite{Cusick1973} aktualisiert

### M74-EN-2: Beweis n=3 fehlerhafte Formel bereinigt
- Unverstaendliche Inverse-Notation entfernt

### M74-EN-3: Nicht referenzierte Bibliographie entfernt
- SpraggeLooker1999 und SuSurvey2009 entfernt

### M75-EN-1: SklienaSurvey entfernt
- Unreferenzierter Eintrag entfernt

### M75-EN-2 + M75-DE-2: Rosa-Zerlegung Kopienanzahl
- "2m copies/Kopien" zu "2m+1 copies/Kopien" korrigiert (beide Papers, EN+DE)

### M75-DE-1: Englisches Wort "conjectured"
- "conjectured anmutig" zu "vermutlich anmutig" korrigiert
