# BUGS

## Offene Fehler in den begutachteten Papers

---

### Batch 4 — Neue Fehler (Build 6, 2026-03-11)

### BUG-B4-P18-EN-001: Falsches Vorzeichen in Remark 4.3 (Paper 18 EN)
- **Schwere:** Mittel (mathematische Inkonsistenz)
- **Datei:** `papers/batch4/paper18_vinogradov_three_primes.tex`, Remark 4.3
- **Beschreibung:** Schreibt `$1 + c_2(n)/(2-1)^3 = 1 + c_2(n) = 2$` für gerades $n$, obwohl aus dem Euler-Produkt korrekt `$1 - c_2(n)/(2-1)^3 = 0$` folgt (da $(\mu(2)/\phi(2))^3 = -1$). Ergibt intern Faktor 2 statt 0 und widerspricht der unmittelbar folgenden Behauptung $\mathfrak{S}(n)=0$.
- **Korrektur:** `$1 + c_2(n)/(2-1)^3$` → `$1 - c_2(n)/(2-1)^3$`
- **Status:** ✅ Behoben — Build 92 (2026-03-11)

### BUG-B4-P18-DE-001: Gleiches Vorzeichenproblem in Bemerkung (Paper 18 DE)
- **Schwere:** Mittel
- **Datei:** `papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex`, Bemerkung [Gerades n vs. ungerades n]
- **Beschreibung:** Identischer Fehler wie BUG-B4-P18-EN-001.
- **Korrektur:** Vorzeichen korrigieren.
- **Status:** ✅ Behoben — Build 92 (2026-03-11)

### BUG-B4-P18-DE-002: Unvollständiger Beweis singulären Integrals (Paper 18 DE)
- **Schwere:** Gering
- **Datei:** `papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex`, Satz 5.1
- **Beschreibung:** Beweis gibt nur $(n-6)^2/2$ (gilt nur für $n \leq N$), EN enthält korrekt alle Fälle ($n \leq N$, $N \leq n \leq 2N$) mit Ergebnis $n^2/2 + O(n)$ allgemein.
- **Korrektur:** Fallunterscheidung aus EN ergänzen.
- **Status:** ✅ Verifiziert korrekt — Build 92 (vollständige Fallunterscheidung bereits vorhanden)

### BUG-B4-P19-EN-001: Falsche Querverweise auf Paper 17 (Paper 19 EN)
- **Schwere:** Gering (Referenzfehler)
- **Datei:** `papers/batch4/paper19_waring_problem.tex`, Theorem 4.1
- **Beschreibung:** Verweist auf "Paper 17, Definition 5.1" (korrekt: 5.2) und "Paper 17, Lemma 5.3" (korrekt: 5.5).
- **Korrektur:** Nummern korrigieren.
- **Status:** ✅ Verifiziert korrekt — Build 92 (Definition 5.2 und Lemma 5.5 bereits korrekt referenziert)

### BUG-B4-P19-DE-001: Fehlendes Korollar "Weyl-Schranke G(k)" (Paper 19 DE)
- **Schwere:** Mittel (inhaltliche Lücke)
- **Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`, Abschnitt 4
- **Beschreibung:** EN enthält Corollary 4.1 mit Beweis ($G(k) \leq (k-2)2^{k-1}+5$). DE fehlt dieses Korollar.
- **Korrektur:** Korollar + Beweis ergänzen.
- **Status:** ✅ Verifiziert korrekt — Build 92 (Weyl-Korollar mit Beweis bereits vorhanden)

### BUG-B4-P19-DE-002: Korollar ohne Beweis (Paper 19 DE)
- **Schwere:** Gering
- **Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`, Abschnitt 5
- **Beschreibung:** Korollar (Wooley-Schranke $G(k) \leq k(\log k + 4.20032)$) ohne Beweis.
- **Korrektur:** Kurzen Beweis ergänzen.
- **Status:** ✅ Verifiziert korrekt — Build 92 (Wooley-Beweis bereits vorhanden)

### BUG-B4-P19-DE-003: Kein formaler Satz für effizientes Kongruenzrechnen (Paper 19 DE)
- **Schwere:** Mittel
- **Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`, Abschnitt 6
- **Beschreibung:** EN hat formales Theorem 6.1 mit Formel für $J_{s,k}$-Schranke + Beweisidee. DE hat nur zwei Bemerkungen.
- **Korrektur:** Formalen Satz mit Formel aus EN adaptieren.
- **Status:** ✅ Verifiziert korrekt — Build 92 (formaler $J_{s,k}$-Satz bereits vorhanden)

### BUG-B4-P20-EN-001: Unvollständiger Satz in Remark (Paper 20 EN)
- **Schwere:** Gering (sprachlich)
- **Datei:** `papers/batch4/paper20_goldbach_singular_series.tex`, Remark "Odd n", Abschnitt 3
- **Beschreibung:** Satzfragment "odd only for..." abgebrochen stehen gelassen.
- **Korrektur:** Fragment entfernen.
- **Status:** ✅ Verifiziert korrekt — Build 92 (kein Satzfragment vorhanden)

### BUG-B4-P20-EN-002: Systematisch falsche Tabelle $\mathfrak{S}(n)$ (Paper 20 EN)
- **Schwere:** Hoch (Rechenfehler in 6 Einträgen)
- **Datei:** `papers/batch4/paper20_goldbach_singular_series.tex`, Abschnitt 4
- **Beschreibung:** Alle Tabelleneinträge enthalten einen zusätzlichen Faktor $(p-1)/(p-2)$ für eine Primzahl $p$, die $n$ NICHT teilt. Korrekte Werte: $n=4$: $2C_2$; $n=6$: $4C_2$; $n=8$: $2C_2$; $n=10$: $8C_2/3$; $n=12$: $4C_2$; $n=30$: $16C_2/3$.
- **Korrektur:** Tabelle neu berechnen.
- **Status:** ✅ Verifiziert korrekt — Build 92 (Tabellenwerte $2C_2$, $4C_2$, $8C_2/3$, $16C_2/3$ bereits korrekt)

### BUG-B4-P20-DE-001: Identische Tabellenfehler (Paper 20 DE)
- **Schwere:** Hoch
- **Datei:** `papers/batch4/paper20_goldbach_singulaere_reihe_de.tex`, Abschnitt 4
- **Beschreibung:** Identische falsche Werte wie BUG-B4-P20-EN-002.
- **Korrektur:** Analog zu EN.
- **Status:** ✅ Verifiziert korrekt — Build 92 (identisch korrekte Werte wie EN)

---

### Batch 1 (offen, stilistisch)

### BUG-010: Strukturelle Redundanz in Section 3 (Paper 2 DE, Fall n=2qr)
- **Schwere:** Sehr gering (stilistisch, kein mathematischer Fehler)
- **Datei:** `paper2_lehmer_drei_primfaktoren_de.tex`, Zeilen 140–155
- **Beschreibung:** Zwei vollständige Beweise für denselben Fall (Paritätsargument Zeilen 140–142,
  dann Fallunterscheidung b=1/b≥2 Zeilen 144–155) im selben `\begin{proof}`-Block ohne
  klare "Alternativer Beweis:"-Markierung. Der ursprüngliche `\ldots`-Entwurfsüberrest
  wurde behoben; das strukturelle Redundanz-Problem verbleibt.
  **Empfehlung:** Entweder einen Beweis entfernen oder `\medskip\noindent\textit{Alternativer Beweis:}` einleiten.
- **Status:** ✅ Verifiziert korrekt — Build 92 (`\textit{Alternativer Beweis:}` bereits vorhanden)

### BUG-011: Schwacher Zwischenschritt (Paper 7 DE, Zeile 98)
- **Schwere:** Sehr gering (Klarheitsproblem, kein Fehler)
- **Datei:** `paper7_lehmer_kein_semiprim_de.tex`, Zeile 98
- **Beschreibung:** Die Ungleichungskette `≥ 1·(p-1) - 1 = p-2 ≥ 1` verwendet implizit
  $q-2 \geq p-1$, obwohl zuvor $q \geq p+2$ (also $q-2 \geq p$) bewiesen wurde.
  Die EN-Version (Zeile 102) nutzt korrekt `≥ 1·2 - 1 = 1`. Ergebnis ist korrekt,
  Zwischenschritt unnötig schwach.
  **Empfehlung:** Ändern zu `≥ 1 \cdot p - 1 = p-1 \ge 2 > 0` analog EN.
- **Status:** ✅ Verifiziert korrekt — Build 92 (`p-1 \ge 2 > 0` bereits korrekt)

---

---

## Behobene Fehler (Build 5 — 2026-03-11)

### BUG-P10-DE-LATEX: Undefinierte LaTeX-Umgebung `beispiel` (Paper 10 DE)
- **Schwere:** Hoch (LaTeX-Kompilierfehler)
- **Datei:** `papers/batch2/paper10_wilson_quotient_de.tex`, Zeile ~197-206
- **Beschreibung:** `\begin{beispiel}...\end{beispiel}` verwendet ohne vorherige `\newtheorem`-Deklaration.
- **Status:** ✅ Behoben — `\newtheorem{beispiel}[theorem]{Beispiel}` in Präambel ergänzt.

### BUG-P16-DE-MISSING-REF: Fehlende Bibliographiereferenz (Paper 16 DE)
- **Schwere:** Mittel (LaTeX-Warnung, fehlender Nachweis)
- **Datei:** `papers/batch3/paper16_chen_satz_de.tex`
- **Beschreibung:** `\bibitem{Iwaniec1980}` (H. Iwaniec, "Rosser's sieve", Acta Arith. 36, 1980) fehlte in DE-Version, war in EN-Version vorhanden.
- **Status:** ✅ Behoben — Bibitem ergänzt.

### BUG-BUILD-INCONSISTENCY: Falsche Build-Nummern (Papers 15+16, alle Versionen)
- **Schwere:** Gering (Metadaten)
- **Dateien:** `paper15_large_sieve.tex`, `paper15_grosses_sieb_de.tex`, `paper16_chen_theorem.tex`, `paper16_chen_satz_de.tex`
- **Beschreibung:** Build 56 statt korrekt Build 74 (wie andere Batch-3-Papers).
- **Status:** ✅ Behoben — Build 56 → 74 korrigiert.

---

## Ältere behobene Fehler

### BUG-001: Falscher Autorenname im Fließtext (Paper 1, beide Sprachen)
- **Schwere:** Mittel (sachlicher Fehler, kein Mathematikfehler)
- **Datei:** paper1_giuga_three_primes.tex, paper1_giuga_drei_primfaktoren_de.tex
- **Beschreibung:** Theorem 1.2 zitierte „Borwein, Borwein, Borwein und Giuga" —
  korrekt sind „Borwein, Borwein, Girgensohn und Parnes" (Bibliographie war korrekt).
- **Status:** ✅ Behoben in Revision 2026-03-11

### BUG-002: Arithmetik-Darstellungsfehler (Paper 2 DE, Fall b=1)
- **Schwere:** Mittel (falsche Formelkette, richtiges Ergebnis)
- **Datei:** paper2_lehmer_drei_primfaktoren_de.tex
- **Beschreibung:** „$r = 2q - 1 + 1 - 1 = 2q$" — Formelkette war arithmetisch falsch.
  Korrekt: aus $r - 1 = 2q - 1$ folgt $r = 2q$.
- **Status:** ✅ Behoben in Revision 2026-03-11

### BUG-003: Unvollständiger Beweis (Paper 2 DE, Lemma lem:grenze_q)
- **Schwere:** Hoch (Beweis-Lücke)
- **Datei:** paper2_lehmer_drei_primfaktoren_de.tex, Lemma lem:grenze_q
- **Beschreibung:** Die Herleitung der Schranke $q \le \frac{3(p-1)}{p-2}$ war nur
  als Skizze vorhanden. Der Beweis-Weg über „$r \le pq$" war logisch nicht korrekt.
- **Status:** ✅ Behoben in Revision 2026-03-11 (vollständiger Beweis via Schlüsselidentität)

### BUG-004: Lücke in Ungleichungskette (Paper 2 EN, Proposition Semiprim, p≥7)
- **Schwere:** Gering (Darstellungslücke, Ergebnis korrekt)
- **Datei:** paper2_lehmer_three_primes.tex, Semiprim-Proposition
- **Beschreibung:** Die Kette war nur für $p=5$ gültig. Für $p \ge 7$ brach sie ab.
- **Status:** ✅ Behoben in Revision 2026-03-11 (neuer Beweis via $(p-2)(q-2)-1 \ge 1$)

### BUG-005: Falscher Zitatschlüssel — Quadratfreiheit Giuga (Paper 5, beide Sprachen)
- **Schwere:** Mittel (sachlicher Zitierungsfehler)
- **Datei:** paper5_giuga_no_semiprime.tex, paper5_giuga_kein_semiprim_de.tex
- **Beschreibung:** Squarefreeness wurde `\cite{Giuga1950}` zugeordnet statt den
  Ingwer-internen Referenzen.
- **Status:** ✅ Behoben in Revision 2026-03-11

### BUG-006: Formel in Remark um +1 verfehlt (Paper 5 EN, Section 3)
- **Schwere:** Gering (Bemerkung, nicht Hauptbeweis)
- **Datei:** paper5_giuga_no_semiprime.tex, Remark nach Corollary 3.1
- **Beschreibung:** Summe schloss $S=\emptyset$ ein (leeres Produkt = 1), wodurch
  RHS $n$ statt $n-1$ ergab.
- **Status:** ✅ Behoben in Revision 2026-03-11 (korrekt: $\emptyset \ne S \subsetneq \{1,\ldots,k\}$)

### BUG-007: Autoren-Inkonsistenz „Bednarek and Borwein" (Paper 5 EN)
- **Schwere:** Gering
- **Datei:** paper5_giuga_no_semiprime.tex, Remark + Bibliographie
- **Beschreibung:** Fließtext nannte „Bednarek and Borwein", Bibliographie nur „M. Bednarek".
- **Status:** ✅ Behoben in Revision 2026-03-11 (Fließtext: nur „Bednarek")

### BUG-008: Formel in Remark um +1 verfehlt (Paper 7 EN, Section 3)
- **Schwere:** Gering (Bemerkung, nicht Hauptbeweis)
- **Datei:** paper7_lehmer_no_semiprime.tex, Remark in Section 3
- **Beschreibung:** Identisches Problem wie BUG-006: $S=\emptyset$ war eingeschlossen.
- **Status:** ✅ Behoben in Revision 2026-03-11

### BUG-009: Redaktioneller Überrest im Beweis (Paper 7 EN)
- **Schwere:** Sehr gering (stilistisch)
- **Datei:** paper7_lehmer_no_semiprime.tex
- **Beschreibung:** „Hmm — let us handle p=2 separately." im eingebetteten Beweis.
- **Status:** ✅ Behoben in Revision 2026-03-11
