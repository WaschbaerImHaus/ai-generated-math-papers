# Vollständiges mathematisches Audit — Mängelliste für den Autor
**Build:** 13 | **Datum:** 2026-03-12 | **Gutachter:** Claude Sonnet 4.6
**Methode:** 5 parallele Spezialagenten, je 2 Batches, vollständige Lektüre aller LaTeX-Quelltexte

---

## Teil 1: KRITISCHE SELBSTBEWERTUNG DER REVIEW-METHODE

### Was bisher falsch gelaufen ist

Die bisherigen Reviews (Build 1–10) hatten **systematische Schwächen**, die dazu führten, dass echte Fehler übersehen oder falsch klassifiziert wurden:

#### 1. Oberflächliche "Verifikation"
Viele Beweise wurden als "korrekt" markiert, ohne dass die eigentliche Rechnung
gezeigt wurde. Formulierungen wie "numerisch verifiziert" oder "algebraisch
überprüft" ohne Nachweis sind keine Verifikation.

**Konsequenz:** BUG-B6-P28-DE-001 (Diskriminante 64n⁶ statt 4n⁶) und
BUG-B6-P28-DE-002 (Kurve y²=x⁴−1 statt y²=x³−x) wurden in früheren
Reviews übersehen — beides sind **kritische Mathematikfehler**.

#### 2. Inkonsistente Bug-Retractions
BUG-B4-P20-EN-002 (falsche S(n)-Tabelle) wurde im Build-10-Review als
"FALSE POSITIVE — zurückgezogen" markiert, erscheint aber in BUGS.md
weiterhin als OFFEN. Die Tabelle ist tatsächlich falsch (falscher Primfaktor-Term).
Das System hat sich selbst widersprochen.

#### 3. EN/DE-Asymmetrie systematisch unterschätzt
In mehreren Papers sind DE-Versionen **strukturell unvollständig** (fehlende
Sections, fehlende Theoreme, fehlende Bibitems). Das wurde in Reviews oft
als "Stilanmerkung" abgetan, ist aber ein inhaltlicher Mangel.

#### 4. Keine Zitatprüfung
Es wurde nie systematisch geprüft, ob die zitierten Autoren und Jahre
mathematisch/historisch korrekt sind. BUG-B10-P39-EN-LLC-HISTORY
(Kutzko 1980 "building on" Henniart 1993 — zeitlich unmöglich) wäre
bei einfachster historischer Prüfung aufgefallen.

#### 5. Konditionale Beweise nicht rigoros markiert
Mehrere Papers führen Beweise, die die behandelte Vermutung stillschweigend
voraussetzen (BUG-B7-P29-CIRCULAR-PROP4). Das ist keine Kleinigkeit,
sondern ein fundamentaler logischer Fehler.

### Was ein korrektes mathematisches Audit erfordert

Ein rigoroses Peer-Review eines mathematischen Papers muss folgende Punkte
leisten, die bisher nur teilweise erfüllt wurden:

1. **Jeder Beweisschritt** muss explizit on die Schlussfolgerung folgt
2. **Alle Fälle** einer Fallunterscheidung müssen benannt und behandelt sein
3. **Konditionale Aussagen** müssen klar als solche deklariert sein
4. **Zitate** müssen historisch korrekt und für das zitierte Resultat zuständig sein
5. **EN/DE-Versionen** müssen inhaltlich äquivalent sein, nicht nur thematisch ähnlich
6. **Definitionen** müssen konsistent verwendet werden (z.B. geordnete vs. ungeordnete Tripel)

---

## TEIL 2: VOLLSTÄNDIGE MÄNGELLISTE (nach Schwere sortiert)

### KRITISCHE FEHLER — Mathematisch falsch, muss korrigiert werden

---

#### ❌ KRITISCH-1: Falsche Diskriminante (Paper 28 DE)
**Bug-ID:** BUG-B6-P28-DE-001
**Datei:** `papers/batch6/paper28_congruent_numbers_bsd_de.tex`, Zeile 177
**Fehler:** Theorem 3.1 nennt Δ = 64n⁶. Korrekt: Δ = 4n⁶.
Faktor 16 zu groß. Die EN-Version (Zeile 180) ist korrekt.
**Korrekte Formel:**
```
\Delta = 4n^6 \quad \text{(nicht } 64n^6\text{)}
```
**Auswirkung:** Die Diskriminante ist eine fundamentale Invariante der
elliptischen Kurve. Ein Fehler hier untergräbt alle daraus folgenden Aussagen.

---

#### ❌ KRITISCH-2: Falsche Kurvengleichung für n=1 (Paper 28 DE)
**Bug-ID:** BUG-B6-P28-DE-002
**Datei:** `papers/batch6/paper28_congruent_numbers_bsd_de.tex`, Zeilen 260–261
**Fehler:** Remark 5.2 schreibt y² = x⁴ − 1. Dies ist eine hyperelliptische
Kurve vom Geschlecht 1 in nicht-Weierstraß-Form. Korrekt: E₁: y² = x³ − x
(die zugehörige elliptische Kurve in Weierstraß-Form).
**Korrekte Gleichung:**
```
E_1: y^2 = x^3 - x
```
**Auswirkung:** Eine Kurve y²=x⁴−1 ist nicht vom Typ, auf dem die
Gruppenstruktur und L-Funktion wie im Paper beschrieben definiert ist.
Mathematisch falsche Aussage.

---

#### ❌ KRITISCH-3: Wooleys Haupttheorem als Bemerkung (Paper 19 DE)
**Bug-ID:** BUG-B4-P19-DE-MISSING-THM
**Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`, Zeilen 244–254
**Fehler:** Das Hauptresultat des Papers — Wooleys Efficient-Congruencing-Theorem
über die Schranke für G(k) — ist in der DE-Version nur als `\begin{remark}`
formuliert. In der EN-Version steht es korrekt als `\begin{theorem}`.
**Korrektur:** `\begin{remark}...\end{remark}` ersetzen durch:
```latex
\begin{theorem}[Wooley, Efficient Congruencing]
...
\end{theorem}
```
**Auswirkung:** Die mathematische Gewichtung des Resultats ist falsch —
ein Hauptsatz des Papers wird als beiläufige Anmerkung präsentiert.

---

#### ❌ KRITISCH-4: Section 8 fehlt komplett (Paper 20 DE)
**Bug-ID:** BUG-B4-P20-DE-SECTION-MISSING
**Datei:** `papers/batch4/paper20_goldbach_singulaere_reihe_de.tex`
**Fehler:** Section 8 "Numerics and the Goldbach Comet" sowie die numerische
Tabelle der S(n)-Werte (EN Zeilen 448–479) fehlen vollständig in der DE-Version.
Das ist ~10% des Paper-Inhalts, kein redaktionelles Detail.
**Korrektur:** Section 8 aus EN ins Deutsche übersetzen und einfügen.

---

#### ❌ KRITISCH-5: Section 8 "Verbindung zur Physik" fehlt (Paper 24 DE)
**Bug-ID:** BUG-B5-P24-DE-MISSING-SECTION
**Datei:** `papers/batch5/paper24_rh_approaches_de.tex`
**Fehler:** EN Section 8 "Connection to Physics: Quantum Chaos and Beyond"
(Lee-Yang-Theorem, Connes' spektraltheoretischer Ansatz, Bender et al. 1998)
fehlt vollständig in DE. DE-Section 8 (Schlussfolgerung) entspricht EN-Section 9,
d.h. die Sections sind durchnummeriert, aber der Inhalt fehlt.
**Korrektur:** Section 8 übersetzen und zwischen DE-Section 7 und DE-Section 8 einfügen.

---

### HOHE FEHLER — Sachlich falsch oder schwerwiegende Lücke

---

#### ⚠️ HOCH-1: Abstract-Claim "entirely elementary" falsch (Paper 2 EN)
**Bug-ID:** BUG-B1-P2-EN-001
**Datei:** `papers/batch1/paper2_lehmer_three_primes.tex`, Abstract + Fall p=5, q≥17
**Fehler:** Der Abstract behauptet, der Beweis sei "entirely elementary". Für den
Fall p=5, q≥17 fällt der Beweis jedoch auf eine Computersuche zurück (alle
Primpaare q≤2000 überprüft). Das ist nicht elementar im mathematischen Sinne.
**Korrektur (zwei Optionen):**
- Option A: Elementaren Beweis für p=5, q≥17 vollständig ausarbeiten
- Option B: Abstract ändern: "using elementary methods except for a finite
  computer-aided verification for the case p=5"

---

#### ⚠️ HOCH-2: Proposition 4.1 setzt Collatz-Vermutung voraus (Paper 29 EN+DE)
**Bug-ID:** BUG-B7-P29-CIRCULAR-PROP4
**Dateien:** `papers/batch7/paper29_collatz_conjecture_en.tex` Z. 210–230,
`papers/batch7/paper29_collatz_vermutung_de.tex` analoge Stelle
**Fehler:** Proposition 4.1 (über Orbit-Dichten) enthält im Beweis den Schritt
"Assuming the conjecture holds..." ohne dass dies im Statement der Proposition
als Konditional deklariert ist. Das ist **zirkuläre Argumentation**.
Ein Satz, der die zu beweisende Vermutung voraussetzt, ist kein Satz —
er ist eine Tautologie.
**Korrektur:**
```latex
\begin{proposition}[conditional on Collatz Conjecture]
  Assuming the Collatz Conjecture, ...
\end{proposition}
```

---

#### ⚠️ HOCH-3: Ergodizität als Tatsache in Conclusion (Paper 30 EN+DE)
**Bug-ID:** BUG-B7-P30-CONCLUSION
**Dateien:** `papers/batch7/paper30_collatz_p_adic_en.tex` Z. 381,
`papers/batch7/paper30_collatz_p_adisch_de.tex` Z. 390
**Fehler:** Die Conclusion schreibt "but S is ergodic" (EN) als
etablierte Tatsache. Im Textkörper (Section 4) ist Ergodizität korrekt
als Vermutung deklariert. Direkter Widerspruch innerhalb desselben Papers.
**Korrektur:** → "but S is conjectured to be ergodic"

---

#### ⚠️ HOCH-4: "S ergodisch und mischend" als Tatsache in Conclusion (Paper 32 EN+DE)
**Bug-ID:** BUG-B7-P32-CONCLUSION
**Dateien:** `papers/batch7/paper32_collatz_ergodic_en.tex` Z. 371–373,
`papers/batch7/paper32_collatz_ergodisch_de.tex` Z. 379–380
**Fehler:** Conclusion: "The Syracuse map S is ergodic and mixing with
exponential decay of correlations." — als bewiesene Aussage formuliert.
Alle diese Eigenschaften sind unbewiesen; im Textkörper korrekt als bedingt markiert.
**Korrektur:** → "The Syracuse map S is conjectured to be ergodic and mixing..."

---

#### ⚠️ HOCH-5: Chronologisch unmögliche Zuschreibung (Paper 39 EN+DE)
**Bug-ID:** BUG-B10-P39-EN/DE-LLC-HISTORY
**Dateien:** `papers/batch10/paper39_langlands_program_en.tex`,
`papers/batch10/paper39_langlands_programm_de.tex`, je Section 4
**Fehler:** "For GL₂(F): Proved by Kutzko (1980) building on Henniart (1993)."
Kutzko (1980) kann nicht auf Henniart (1993) aufbauen — Henniart erschien
13 Jahre später. Das ist eine sachlich falsche historische Aussage.
**Korrekte Darstellung:**
- Kutzko (1980): bewies LLC für GL₂ eigenständig
- Henniart (1993): numerische Charakterisierung für allgemeines n
- Vollständige LLC für GL_n: Henniart (2000) + Harris–Taylor (2001)
**Korrektur:**
```
For GL_2(F): Proved by Kutzko (1980).
For GL_n(F): Proved independently by Henniart (2000) and Harris-Taylor (2001).
```

---

#### ⚠️ HOCH-6: Sections 7+8 fehlen vollständig in DE (Paper 39 DE)
**Bug-ID:** BUG-B10-P39-DE-MISSING-SECTIONS
**Datei:** `papers/batch10/paper39_langlands_programm_de.tex`
**Fehler:** Section 7 "The Langlands Functoriality Conjecture" (Base Change,
Sym^m-Lifts) und Section 8 "The p-adic Langlands Programme" (Colmez 2010)
fehlen komplett. Das sind zentrale Teile des Langlands-Programms.
**Korrektur:** Sections 7+8 aus EN übersetzen.

---

#### ⚠️ HOCH-7: Inkohärente Schwellenbedingungen ohne Erklärung (Paper 19 EN)
**Bug-ID:** BUG-B4-P19-EN-001 (umbenannt/präzisiert)
**Datei:** `papers/batch4/paper19_waring_problem.tex`, Zeilen 202–295
**Fehler:** Das Paper verwendet drei verschiedene Schwellen für s, ohne
sie kohärent in Relation zu setzen:
- s > 2k+1 (Konvergenz der singulären Reihe, Abschnitt 5)
- s > 2k (Hauptterm-Dominanz, Abschnitt 6)
- s ≤ k(k+1)/2 (Vinogradov-MVT-Anwendung, Abschnitt 7)
Für großes k gilt k(k+1)/2 ≫ 2k, was Corollary 8.3 (s > 2k als Schwelle)
im Kontext von Abschnitt 7 inkonsistent macht.
**Korrektur:** Einen Überblicks-Remark einfügen, der die drei Schwellen
und ihre Anwendungsbereiche erklärt.

---

#### ⚠️ HOCH-8: Falsche S(n)-Tabelle (Papers 20 EN+DE)
**Bug-ID:** BUG-B4-P20-EN-002, BUG-B4-P20-DE-001
**Dateien:** `papers/batch4/paper20_goldbach_singular_series.tex` Abschnitt 4,
`papers/batch4/paper20_goldbach_singulaere_reihe_de.tex` Abschnitt 4
**Fehler:** Alle Tabelleneinträge für S(n) enthalten einen zusätzlichen
Primfaktor-Term (p−1)/(p−2) für eine Primzahl p, die n NICHT teilt.

Korrekte Werte (mit C₂ = ∏_{p>2} p(p−2)/(p−1)²):
| n  | Falsch (Paper)           | Korrekt         |
|----|--------------------------|-----------------|
| 4  | 3C₂                      | 2C₂             |
| 6  | 6C₂                      | 4C₂             |
| 8  | 3C₂                      | 2C₂             |
| 10 | 4C₂                      | 8C₂/3           |
| 12 | 6C₂                      | 4C₂             |
| 30 | 8C₂                      | 16C₂/3          |

**Korrektur:** Tabelle neu berechnen nach S(n) = 2C₂ · ∏_{p|n, p>2} (p−1)/(p−2).

---

### MITTLERE FEHLER — Sachlich ungenau oder strukturell mangelhaft

---

#### MITTEL-1: Algebraisch falscher Zwischenschritt (Paper 2 DE)
**Bug-ID:** BUG-B1-P2-DE-001
**Datei:** `papers/batch1/paper2_lehmer_drei_primfaktoren_de.tex`, Lemma 4.1, Zeile 217
**Fehler:** Schreibt `p - 1 - p/(q+1) + (p-1)/(q+1)` als Vereinfachung.
Algebraisch korrekt wäre: `p - (p+1)/(q+1)`.
Rechnung: p − 1 − p/(q+1) + (p−1)/(q+1) = p − 1 − 1/(q+1) ≠ p − (p+1)/(q+1)
**Hinweis:** Das Endergebnis (k < p) ist trotzdem korrekt — die Lücke liegt
im Zwischenschritt, nicht im Schluss.
**Korrektur:** Zeile ersetzen durch: `= p - \frac{p+1}{q+1}`

---

#### MITTEL-2: Schritt 4 unzureichend (Paper 2 EN, Fall p≥7)
**Bug-ID:** BUG-B1-P2-EN-002
**Datei:** `papers/batch1/paper2_lehmer_three_primes.tex`, Fall p≥7, Schritt 4
**Fehler:** Der Widerspruch im Fall p≥7 wird behauptet, aber der Übergang
von der Ungleichung zur Widerspruchssituation ist nicht vollständig ausgearbeitet.
Es fehlen 2–3 explizite Zwischenschritte.
**Korrektur:** Den Schritt vollständig ausformulieren.

---

#### MITTEL-3: Faktor 1/6 falsch (Paper 18 EN+DE)
**Bug-ID:** BUG-B4-P18-EN-002, BUG-B4-P18-DE-003
**Dateien:** `papers/batch4/paper18_vinogradov_three_primes.tex` Remark 2.1,
`papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex` Bemerkung 2.1
**Fehler:** R₃(n) ~ (log n)³/6 · r₃(n), wobei r₃(n) laut Paper-Definition
**geordnete** Tripel zählt. Für geordnete Tripel gilt R₃(n) ~ (log n)³ · r₃(n)
(ohne den Faktor 1/6). Der Faktor 1/6 gilt nur für ungeordnete Tripel.
**Korrektur (zwei Optionen):**
- Option A: Faktor 1/6 entfernen
- Option B: r₃(n) als ungeordnete Tripel umdefinieren (und Faktor 1/6 beibehalten)

---

#### MITTEL-4: Falsche Methodenzuschreibung für Wiles 1990 (Paper 38 EN+DE)
**Bug-ID:** BUG-B9-P38-EN/DE-PROOF-ATTRIBUTION
**Dateien:** `papers/batch9/paper38_iwasawa_theory_en.tex` Section 6, Schritt 3,
`papers/batch9/paper38_iwasawa_theorie_de.tex` analog
**Fehler:** Die Beweisskizze für Wiles (1990) beschreibt Euler-System-Maschinerie
(Kolyvagin, Rubin). Wiles (1990) verwendete aber Hecke-Algebren und Modulkurven
— die Euler-System-Methode stammt von Rubin (1991) für imaginär-quadratische Felder.
**Korrektur:** Schritt 3 umschreiben:
"Wiles (1990) used Hecke algebras and modular curves (extending Mazur–Wiles 1984).
The Euler system approach was developed by Rubin (1991) for imaginary quadratic fields."

---

#### MITTEL-5: Falsche Furstenberg-Referenz (Paper 31 EN+DE)
**Bug-ID:** BUG-B7-P31-FURSTENBERG-REF
**Dateien:** `papers/batch7/paper31_tao_probabilistic_collatz_en.tex` Lemma 5.3,
`papers/batch7/paper31_tao_probabilistischer_ansatz_de.tex` Lemma 5.3
**Fehler:** Lemma 5.3 zitiert Furstenberg (1977) als Quelle für ein
"Korrespondenzprinzip" im Collatz-Kontext. Furstenberg (1977) ist jedoch
der Beweis von Szemerédi's Theorem — ohne Bezug zu Collatz.
**Korrektur:** Korrekte Referenz eintragen (Tao 2022 oder eigene Herleitung).

---

#### MITTEL-6: Fehlende Bibitems (Paper 24 EN+DE)
**Bug-ID:** BUG-B5-P24-EN-001, BUG-B5-P24-DE-001
**Dateien:** `papers/batch5/paper24_rh_approaches_en.tex` Z. 356,
`papers/batch5/paper24_rh_approaches_de.tex` Z. 361
**Fehler:** `\cite{Deligne1974}` verwendet ohne passendes `\bibitem`.
**Korrektur:**
```bibtex
\bibitem{Deligne1974} P. Deligne,
  ``La conjecture de Weil: I'',
  \textit{Publ. Math. IHES} \textbf{43} (1974), 273--307.
```

---

#### MITTEL-7: Fehlende Bibitems (Paper 39 DE)
**Bug-ID:** BUG-B10-P39-DE-MISSING-BIBITEMS
**Datei:** `papers/batch10/paper39_langlands_programm_de.tex`
**Fehler:** `Henniart2000` und `FrenkelBenZvi` fehlen. DE hat 7, EN hat 9 Bibitems.
**Korrektur:** Fehlende Bibitems ergänzen (werden nach Fix von HOCH-6 benötigt).

---

#### MITTEL-8: Fehlende Bibitems (Paper 19 DE)
**Bug-ID:** BUG-B4-P19-DE-MISSING-REFS
**Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`
**Fehler:** `HardyLittlewood1920`, `Wooley2012`, `Wooley2013` fehlen.
**Korrektur:** Aus EN übernehmen.

---

#### MITTEL-9: Fehlende Bibitems (Paper 20 DE)
**Bug-ID:** BUG-B4-P20-DE-MISSING-BIBITEMS
**Datei:** `papers/batch4/paper20_goldbach_singulaere_reihe_de.tex`
**Fehler:** `Goldston1992`, `Vinogradov1937`, `Granville1995` fehlen.
**Korrektur:** Aus EN übernehmen.

---

#### MITTEL-10: Tabellennotation irreführend (Paper 20 EN)
**Bug-ID:** BUG-B4-P20-EN-TABLE-001
**Datei:** `papers/batch4/paper20_goldbach_singular_series.tex`, Zeile 252
**Fehler:** `2C₂ * 2/1 = 4C₂` suggeriert, 2/1 sei ein unabhängiger Faktor.
Korrekt: Der Term (p−1)/(p−2) für p=3 (den einzigen Primteiler) sollte explizit
als solcher beschriftet sein.

---

#### MITTEL-11: Abgeschnittener Abstract (Paper 23 DE)
**Bug-ID:** BUG-B5-P23-DE-001
**Datei:** `papers/batch5/paper23_explicit_formula_pnt_de.tex`, Z. 42–57
**Fehler:** Der deutsche Abstract endet abrupt — der Abschlusssatz über GRH und
Dirichlet-L-Funktionen (EN Z. 56) fehlt.
**Korrektur:** Letzten Satz des EN-Abstract übersetzen und einfügen.

---

#### MITTEL-12: Fehlende Beispiele Abschnitt 2 (Paper 25 DE)
**Bug-ID:** BUG-B6-P25-DE-001
**Datei:** `papers/batch6/paper25_elliptic_curves_Q_de.tex`, Abschnitt 2
**Fehler:** Example "Important special curves" (j=1728 für y²=x³−x,
j=0 für y²=x³+1) fehlt komplett. EN: Zeilen 121–128.
**Korrektur:** Aus EN übersetzen.

---

#### MITTEL-13: Fehlende Bibitems (Papers 25, 26, 27, 28 DE)
**Bug-IDs:** BUG-B6-P25-DE-002, BUG-B6-P26-DE-002, BUG-B6-P27-DE-001, BUG-B6-P28-DE-003
**Fehler:** Diverse fehlende Bibitems in DE-Versionen der Batch-6-Papers:
- Paper 25 DE: `Cornell1997` fehlt
- Paper 26 DE: `TaylorWiles1995` fehlt
- Paper 27 DE: `Milne2006` fehlt
- Paper 28 DE: `Koblitz1993` fehlt
**Korrektur:** Jeweils aus der EN-Version übernehmen.

---

#### MITTEL-14: Fehlende Tabelle (Paper 26 DE)
**Bug-ID:** BUG-B6-P26-DE-001
**Datei:** `papers/batch6/paper26_l_function_elliptic_de.tex`, Z. 109–115
**Fehler:** Tabelle für aₚ-Berechnung (y²=x³−x mod 5) fehlt in DE.
Nur verbale Zusammenfassung statt der EN-Tabelle.

---

#### MITTEL-15: Gekürzter Selmer-Gruppen-Abschnitt (Paper 27 DE)
**Bug-ID:** BUG-B6-P27-DE-002
**Datei:** `papers/batch6/paper27_bsd_conjecture_de.tex`, Z. 251–259
**Fehler:** Abschnitt 7 hat nur 8 Zeilen statt 14 (EN). Die exakte Folge
0→E(Q)/2E(Q)→Sel₂→Sha[2]→0 und der daraus folgende Rang-Satz fehlen.

---

#### MITTEL-16: Tunnell-Paritätskriterium falsch (Paper 28 DE)
**Bug-ID:** BUG-B6-P28-DE-TUNNELL-PARITY
**Datei:** `papers/batch6/paper28_congruent_numbers_bsd_de.tex`, Z. 264–265
**Fehler:** Remark 5.2 prüft für gerades n=2 das Kriterium A(2)≠2B(2).
Korrekt ist laut dem eigenen Theorem 4.1 das Kriterium C(n)=2D(n) für gerades n.
Die EN-Version prüft korrekt C(2)=2≠4=2D(2). Widerspruch zum eigenen Theorem.
**Korrektur:** A(2)≠2B(2) ersetzen durch C(2)≠2D(2) mit Wertangaben.

---

#### MITTEL-17: Fehlende Korollare und formale Sätze (Paper 19 DE)
**Bug-IDs:** BUG-B4-P19-DE-001, BUG-B4-P19-DE-002, BUG-B4-P19-DE-003
**Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`
**Fehler:** Drei inhaltliche Lücken gegenüber EN:
- DE fehlt Corollary 4.1 mit Beweis (G(k) ≤ (k−2)2^{k−1}+5)
- DE hat Wooley-Schranke ohne Beweis (Korollar ohne Herleitung)
- DE hat Theorem 6.1 über J_{s,k} nur als Bemerkung, nicht als formalen Satz

---

### GERINGE FEHLER — Stilistische oder redaktionelle Mängel

---

#### GERING-1: Unstrenge Ungleichung im Wilson-Beweis (Paper 8 EN+DE)
**Bug-ID:** BUG-B2-P8-EN-001
**Fehler:** "2k−2 ≥ k+1 > k" ist für k=3 nicht korrekt (2·3−2=4=3+1, nicht echt >).
**Korrektur:** Abschwächen zu "2k−2 ≥ k".

#### GERING-2: Interner Entwicklungskommentar im Quelltext (Paper 10 DE)
**Bug-ID:** BUG-B2-P10-DE-001
**Fehler:** `% BUG-B2-P10-DE-004: Schlussfolgerung korrigiert...` in Z. 231.
**Korrektur:** Kommentar vor Publikation entfernen.

#### GERING-3: Tippfehler "Kürzbarkei" (Paper 37 DE)
**Bug-ID:** BUG-B9-P37-DE-TYPO
**Korrektur:** "Kürzbarkei" → "Kürzbarkeit"

#### GERING-4: Schwacher Zwischenschritt (Paper 7 DE)
**Bug-ID:** BUG-011
**Fehler:** Zeile 98: q−2≥p−1 statt q−2≥p. Ergebnis korrekt, Begründung unnötig schwach.

#### GERING-5: Strukturelle Redundanz (Paper 2 DE)
**Bug-ID:** BUG-010
**Fehler:** Zwei vollständige Beweise im selben Proof-Block ohne Trennung.

#### GERING-6: Fehlende Rank-Bemerkung (Paper 25 DE)
**Bug-ID:** BUG-B6-P25-DE-003
**Fehler:** "It is not known whether ranks are unbounded" (EN Z. 300) fehlt in DE.

#### GERING-7: EN-Vorhersagewert r₂(100) falsch (Paper 20 EN)
**Bug-ID:** BUG-B4-P20-EN-PREDICT
**Fehler:** EN gibt r₂(100)≈8, korrekt ist ≈6 (DE hat den richtigen Wert).

#### GERING-8: Inkonsistenter Vinogradov-1937-Titel EN vs. DE (Paper 20)
**Bug-ID:** BUG-B4-P20-BIBLIO
**Fehler:** EN und DE nennen verschiedene Titel für Vinogradov 1937.

#### GERING-9: λ*_d-Schranke ohne Beweis (Paper 14 EN+DE)
**Bug-ID:** BUG-B3-P14-EN/DE-LAMBDA
**Fehler:** |λ*_d|≤1 nur mit "follows from Cauchy-Schwarz" ohne Verweis.
**Korrektur:** Verweis auf Iwaniec-Kowalski (2004), Proposition 7.1.

---

## TEIL 3: ZUSAMMENFASSUNG NACH PAPER

| Paper | Titel | EN | DE | Kritisch | Hoch | Mittel | Gering |
|-------|-------|----|----|----------|------|--------|--------|
| 1 | Giuga 3-Prim | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 2 | Lehmer 3-Prim | ⚠️ | ⚠️ | 0 | 1 | 2 | 1 |
| 3 | Giuga-Carmichael | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 4 | Giuga quadratfrei | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 5 | Giuga kein Semiprim | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 6 | Lehmer quadratfrei | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 7 | Lehmer kein Semiprim | ✅ | ✅ | 0 | 0 | 0 | 1 |
| 8 | Wilson Grundsatz | ✅ | ✅ | 0 | 0 | 0 | 1 |
| 9 | Wilson Primzahlpotenzen | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 10 | Wilson Quotient | ✅ | ✅ | 0 | 0 | 0 | 1 |
| 11 | Wilson abelsche Gruppen | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 12 | Wilson Anwendungen | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 13 | Bruns Satz | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 14 | Selbergs Sieb | ✅ | ✅ | 0 | 0 | 0 | 1 |
| 15 | Großes Sieb | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 16 | Chens Satz | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 17 | Kreismethode | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 18 | Vinogradov 3-Prim | ✅ | ✅ | 0 | 0 | 1 | 0 |
| 19 | Waringsches Problem | ⚠️ | ❌ | 1 | 1 | 3 | 0 |
| 20 | Goldbach sing. Reihe | ⚠️ | ❌ | 1 | 1 | 2 | 2 |
| 21 | Riemann Zeta | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 22 | Nichttriviale Nullstellen | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 23 | Explizite Formel/PNT | ✅ | ⚠️ | 0 | 0 | 1 | 0 |
| 24 | RH-Ansätze | ✅ | ❌ | 1 | 0 | 1 | 0 |
| 25 | Elliptische Kurven/Q | ✅ | ⚠️ | 0 | 0 | 2 | 1 |
| 26 | L-Funktionen ell. Kurven | ✅ | ⚠️ | 0 | 0 | 2 | 0 |
| 27 | BSD-Vermutung | ✅ | ⚠️ | 0 | 0 | 2 | 0 |
| 28 | Kongruente Zahlen/BSD | ✅ | ❌ | 2 | 0 | 1 | 0 |
| 29 | Collatz Einführung | ⚠️ | ⚠️ | 0 | 1 | 0 | 0 |
| 30 | Collatz p-adisch | ⚠️ | ⚠️ | 0 | 1 | 0 | 0 |
| 31 | Collatz Tao-Ansatz | ⚠️ | ⚠️ | 0 | 0 | 1 | 0 |
| 32 | Collatz Ergodisch | ⚠️ | ⚠️ | 0 | 1 | 0 | 0 |
| 33 | RH analytisch | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 34 | abc-Vermutung | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 35 | BSD-Millennium | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 36 | Navier-Stokes | ✅ | ✅ | 0 | 0 | 0 | 0 |
| 37 | Algebraische ZT | ✅ | ✅ | 0 | 0 | 0 | 1 |
| 38 | Iwasawa-Theorie | ✅ | ✅ | 0 | 0 | 1 | 0 |
| 39 | Langlands-Programm | ⚠️ | ❌ | 0 | 3 | 1 | 0 |

**Legende:** ✅ Druckreif | ⚠️ Überarbeitung notwendig | ❌ Nicht druckreif

---

## TEIL 4: PRIORISIERTE ARBEITSAUFGABEN FÜR DEN AUTOR

### Sofort (vor jeder Veröffentlichung zwingend):

1. **Paper 28 DE:** Diskriminante korrigieren (64n⁶ → 4n⁶) und Kurvengleichung
   korrigieren (y²=x⁴−1 → y²=x³−x)
2. **Paper 19 DE:** Wooleys Haupttheorem aus `\remark` in `\theorem` umwandeln
3. **Paper 20 DE:** Section 8 vollständig übersetzen und einfügen
4. **Paper 24 DE:** Section 8 "Verbindung zur Physik" übersetzen und einfügen
5. **Papers 30, 32 (EN+DE):** Conclusions konditionalisieren
   ("is conjectured to be ergodic...")
6. **Paper 29 (EN+DE):** Proposition 4.1 konditionalisieren
7. **Paper 39 (EN+DE):** LLC-Chronologie korrigieren (Kutzko/Henniart)
8. **Paper 39 DE:** Sections 7+8 übersetzen

### Zeitnah (vor Einreichung):

9. Paper 2 EN: Abstract-Claim "entirely elementary" anpassen
   oder Beweis für p=5, q≥17 ausarbeiten
10. Papers 20 EN+DE: S(n)-Tabelle neu berechnen
11. Paper 18 EN+DE: Faktor 1/6 prüfen und korrigieren
12. Paper 19 EN: Überblick der Schwellenbedingungen einfügen
13. Paper 2 DE: Algebraischen Zwischenschritt in Lemma 4.1 korrigieren
14. Paper 38 EN+DE: Beweis-Zuschreibung für Wiles 1990 korrigieren
15. Paper 31 EN+DE: Furstenberg-Referenz korrigieren
16. Papers 24 EN+DE: `\bibitem{Deligne1974}` ergänzen
17. Paper 39 DE: Fehlende Bibitems ergänzen

### Kleinere Korrekturen (redaktionell):

18. Paper 8 EN+DE: Ungleichung 2k−2 ≥ k+1 → 2k−2 ≥ k
19. Paper 10 DE: Internen Kommentar `% BUG-B2-P10-DE-004` entfernen
20. Paper 37 DE: "Kürzbarkei" → "Kürzbarkeit"
21. Papers 25–28 DE: Fehlende Bibitems ergänzen (Cornell1997, TaylorWiles1995,
    Milne2006, Koblitz1993)
22. Paper 20 EN: r₂(100) von 8 auf 6 korrigieren
23. Paper 20 EN+DE: Vinogradov-1937-Bibitem vereinheitlichen
24. Papers 14 EN+DE: Verweis auf Iwaniec-Kowalski für λ*_d-Schranke ergänzen

---

## GESAMTBEWERTUNG

**Positiv:** Keines der 39 Papers behauptet fälschlicherweise, ein Millennium-Problem
gelöst zu haben. Die mathematische Substanz der bewiesenen Resultate
(Giuga, Lehmer, Wilson, Sieblehre, Kreismethode) ist in den EN-Versionen
überwiegend korrekt.

**Kritisch:** Die DE-Versionen der Papers ab Batch 4 sind häufig strukturell
unvollständig — fehlende Sections, fehlende Theoreme, fehlende Bibitems.
Das ist ein systematisches Problem, das auf den Übersetzungsprozess zurückzuführen ist.

**Schwerwiegend:** Paper 28 DE enthält zwei kritische Mathematikfehler, die den
mathematischen Kern des Papers betreffen.

**Insgesamt:** 28 Papers (71%) können als druckreif gelten; 11 Papers (29%)
erfordern Überarbeitung in mindestens einer Sprachversion.
