# BUGS.md - specialist-maths

## Offene Bugs
_(keine bekannten Bugs)_

## Behobene Bugs

### BUG-B8-P33-EN-001: paper33 EN — Fehlende fontenc/inputenc-Pakete
- **Datei**: papers/batch8/paper33_riemann_hypothesis_analytic_en.tex
- **Problem**: `\usepackage[T1]{fontenc}` und `\usepackage[utf8]{inputenc}` fehlten im EN-Paper (DE hatte sie)
- **Behoben**: Beide Pakete nach `\usepackage{enumitem,booktabs,array}` eingefügt

### BUG-B8-P34-EN-001: paper34 EN — Fehlende fontenc/inputenc-Pakete
- **Datei**: papers/batch8/paper34_abc_conjecture_en.tex
- **Problem**: `\usepackage[T1]{fontenc}` und `\usepackage[utf8]{inputenc}` fehlten im EN-Paper (DE hatte sie — Inkonsistenz zwischen EN und DE)
- **Behoben**: Beide Pakete nach `\usepackage{enumitem,booktabs,array}` eingefügt

### BUG-B8-P34-DE-001: paper34 DE — Grammatikfehler „unbewiesener"
- **Datei**: papers/batch8/paper34_abc_vermutung_de.tex, Abschnitt 5
- **Problem**: „gilt im Mainstream-Konsens als unbewiesener" — falsches Adjektiv im Prädikativum (kein Komparativ gemeint)
- **Behoben**: „als unbewiesener" → „als nicht bewiesen"

### BUG-B8-P34-DE-002: paper34 DE — Tippfehler „Beweissstrategie"
- **Datei**: papers/batch8/paper34_abc_vermutung_de.tex, Abschnitt 4.2
- **Problem**: Dreifaches „s" in „Beweissstrategie"
- **Behoben**: „Beweissstrategie" → „Beweisstrategie"

### BUG-B8-P35-EN-001: paper35 EN — Falsches Euler-Produkt mit doppeltem Inversen
- **Datei**: papers/batch8/paper35_bsd_conjecture_en.tex, Definition in Abschnitt 2
- **Problem**: `L(E,s) = ∏_p L_p(E,s)^{-1}` ist falsch. Die lokalen Faktoren `L_p(E,s)` sind bereits als `(1-a_p p^{-s}+...)^{-1}` (inverse) definiert. Das Produkt `∏_p L_p(E,s)` (ohne weiteres Inverses) ergibt korrekt die globale L-Funktion.
- **Behoben**: `\prod_p L_p(E,s)^{-1}` → `\prod_p L_p(E,s)`

### BUG-B8-P35-DE-001: paper35 DE — Falsches Euler-Produkt (analog EN)
- **Datei**: papers/batch8/paper35_bsd_vermutung_de.tex, Definition in Abschnitt 2
- **Problem**: Identischer Fehler wie BUG-B8-P35-EN-001
- **Behoben**: `\prod_p L_p(E,s)^{-1}` → `\prod_p L_p(E,s)`

### BUG-B8-P36-EN-001: paper36 EN — Falsche Attribuierung in Ergebnis-Tabelle (Leray vs. Ladyzhenskaya)
- **Datei**: papers/batch8/paper36_navier_stokes_en.tex, Tabelle in Abschnitt 8
- **Problem**: „2D Global, smooth (Leray 1934)" ist falsch. Leray 1934 bewies globale *schwache* Lösungen auch in 2D; die Existenz glatter globaler 2D-Lösungen bewies Ladyzhenskaya 1969.
- **Behoben**: Tabelleneintrag korrigiert auf „Global weak (Leray 1934), smooth (Ladyzhenskaya 1969)"

### BUG-B8-P36-DE-001: paper36 DE — Falsche Attribuierung in Ergebnis-Tabelle (analog EN)
- **Datei**: papers/batch8/paper36_navier_stokes_de.tex, Tabelle in Abschnitt 8
- **Problem**: Identisch zu BUG-B8-P36-EN-001
- **Behoben**: „Global, glatt (Leray 1934)" → „Global schwach (Leray 1934), glatt (Ladyzhenskaya 1969)"

### BUG-B8-P36-EN-002: paper36 EN — bibkey-Jahreszahl-Diskrepanz ConstantinDoering
- **Datei**: papers/batch8/paper36_navier_stokes_en.tex, Bibliographie
- **Problem**: `\bibitem{ConstantinDoering1994}` verwies auf Phys. Rev. E 51 (1995) — Jahreszahl und bibkey stimmten nicht überein
- **Behoben**: Quellenangabe korrigiert auf das korrekte Paper von 1994 (Phys. Rev. E 49)

### BUG-B8-P36-DE-002: paper36 DE — bibkey-Jahreszahl-Diskrepanz ConstantinDoering (analog EN)
- **Datei**: papers/batch8/paper36_navier_stokes_de.tex, Bibliographie
- **Problem**: Identisch zu BUG-B8-P36-EN-002
- **Behoben**: Quellenangabe analog korrigiert



### BUG-B6-P26-EN-NEW-001 (Review Batch6 2026-03-11): paper26 EN — Abstract Euler-Produkt falsch
- **Datei**: papers/batch6/paper26_l_function_elliptic_en.tex, Abstract Zeile 51
- **Problem**: `$L(E,s) = \prod_p L_p(E,s)^{-1}$` im Abstract war inkonsistent mit Definition 3.1, wo `L_p(E,s)^{-1}` als das Polynom `1 - a_p p^{-s} + p^{1-2s}` definiert ist. Mit dieser Notation ist `L_p(E,s) = 1/(1-...)` der konvergente Faktor, also ist `\prod_p L_p(E,s)` (ohne `^{-1}`) das korrekte Euler-Produkt.
- **Behoben**: `\prod_p L_p(E,s)^{-1}` → `\prod_p L_p(E,s)` im Abstract

### BUG-B6-P26-EN-NEW-002 (Review Batch6 2026-03-11): paper26 EN — Fehlende example-Umgebungsdefinition
- **Datei**: papers/batch6/paper26_l_function_elliptic_en.tex, Präambel
- **Problem**: `\begin{example}` wurde 3-mal verwendet (Zeilen 114, 295, 304), aber `\newtheorem{example}[theorem]{Example}` fehlte in der Präambel. DE-Version hatte es korrekt definiert.
- **Behoben**: `\newtheorem{example}[theorem]{Example}` nach `\newtheorem{remark}` eingefügt

### BUG-B6-P26-DE-NEW-001 (Review Batch6 2026-03-11): paper26 DE — Sektionsreihenfolge inkonsistent mit EN
- **Datei**: papers/batch6/paper26_l_function_elliptic_de.tex
- **Problem**: Sektion 6 "Numerische Beispiele" stand VOR Sektion 7 "Die L-Funktion bei s=1". In EN ist die Reihenfolge umgekehrt (erst Definition der Verschwindungsordnung, dann Beispiele), was didaktisch und strukturell korrekt ist.
- **Behoben**: Sektionen getauscht — "Die L-Funktion bei s=1" ist nun Sektion 6, "Numerische Beispiele" Sektion 7

### BUG-B6-P26-DE-NEW-002 (Review Batch6 2026-03-11): paper26 DE — Kolyvagin-Jahr 1989 statt 1990
- **Datei**: papers/batch6/paper26_l_function_elliptic_de.tex, Numerische Beispiele (Rang-0)
- **Problem**: "von Kolyvagin (1989) bedingungslos bewiesen" — bibitem heißt `Kolyvagin1990`, alle anderen Stellen im Paper sagen "1990". Die korrekte Veröffentlichung ist 1990.
- **Behoben**: `(1989)` → `(1990)`

### BUG-B6-P28-EN-004 (Build 102): paper28 EN — Falscher bibitem-Schlüssel Koblitz
- **Datei**: papers/batch6/paper28_congruent_numbers_bsd_en.tex
- **Problem**: bibitem-Schlüssel `Koblitz1993` widersprach dem Erscheinungsjahr der 2. Auflage (1994)
- **Behoben**: Schlüssel auf `Koblitz1994` korrigiert

### BUG-B6-P28-DE-001 (Build 102): paper28 DE — Gleicher falscher bibitem-Schlüssel
- **Datei**: papers/batch6/paper28_congruent_numbers_bsd_de.tex
- **Problem**: Analog zu BUG-B6-P28-EN-004
- **Behoben**: Schlüssel auf `Koblitz1994` korrigiert

### BUG-B6-P28-EN-003 (Build 102): paper28 EN — Falsche Substitutionsvariable im Abstract
- **Datei**: papers/batch6/paper28_congruent_numbers_bsd_en.tex, Zeile 48
- **Problem**: Abstract nannte `x = (b/2)^2` als Substitution, wobei `b` eine Kathete ist. Korrekt ist die Standardparametrisierung `x = (c/2)^2` (Hypotenuse `c`), wie im Beweis von Theorem 2.1 explizit ausgeführt (`x_0 = c^2/4`)
- **Behoben**: `(b/2)^2` → `(c/2)^2`

### BUG-B6-P27-DE-002 (Build 102): paper27 DE — Grammatikfehler "eine endliche Vektorraum"
- **Datei**: papers/batch6/paper27_bsd_conjecture_de.tex, Zeile 263
- **Problem**: "Der Selmer-Raum ist eine endliche $\mathbb{F}_2$-Vektorraum" — maskulines Nomen mit falschem Artikel und Adjektivendung
- **Behoben**: "eine endliche" → "ein endlicher"

### BUG-B6-P27-DE-001 (Build 102): paper27 DE — Falsche Genitivform "Kolyvagings"
- **Datei**: papers/batch6/paper27_bsd_conjecture_de.tex, Zeile 238
- **Problem**: "Kolyvagings Beweis" ist eine falsche Genitivbildung des Eigennamens
- **Behoben**: "Kolyvagings" → "Kolyvagins"

### BUG-B6-P26-DE-001 (Build 102): paper26 DE — Falsches Euler-Produkt
- **Datei**: papers/batch6/paper26_l_function_elliptic_de.tex, Definition 3.2
- **Problem**: `L(E,s) = ∏_p L_p(E,s)^{-1}` ist mathematisch falsch. Da `L_p(E,s)^{-1}` als Polynom `1 - a_p p^{-s} + p^{1-2s}` definiert ist, bedeutet `L_p(E,s) = 1/(1-...)` (der konvergente Faktor). Das Produkt `∏_p L_p(E,s)^{-1} = ∏_p Polynom` ergibt kein konvergentes L-Funktion-Produkt. Korrekt ist `L(E,s) = ∏_p L_p(E,s)`.
- **Behoben**: `\prod_p L_p(E,s)^{-1}` → `\prod_p L_p(E,s)`

### BUG-B6-P26-EN-001 (Build 102): paper26 EN — Falsches Euler-Produkt
- **Datei**: papers/batch6/paper26_l_function_elliptic_en.tex, Definition 3.2
- **Problem**: Analog zu BUG-B6-P26-DE-001: `L(E,s) = ∏_p L_p(E,s)^{-1}` falsch
- **Behoben**: `\prod_p L_p(E,s)^{-1}` → `\prod_p L_p(E,s)`

### BUG-B6-P25-DE-001 (Build 102): paper25 DE — Grammatikfehler Genus "das Invariante"
- **Datei**: papers/batch6/paper25_elliptic_curves_Q_de.tex, Zeile 287
- **Problem**: "ist das geheimnisvollste Invariante" — "Invariante" ist feminin, erfordert Artikel "die"
- **Behoben**: "das geheimnisvollste" → "die geheimnisvollste"

### BUG-B4-P18-EN-001 (Build 92): paper18 EN — Vereinfachung Remark 4.3
- **Datei**: papers/batch4/paper18_vinogradov_three_primes.tex
- **Problem**: In Remark 4.3 stand `$1 + \mu(2)^3 c_2(n)/\phi(2)^3 = 1 - c_2(n)$`
  — das `\mu`-Vorzeichen war implizit, was zu Verwirrung führen konnte
- **Behoben**: Auf `$1 - c_2(n)/(2-1)^3 = 1 - c_2(n)$` vereinfacht (explizites Vorzeichen)

### BUG-B4-P18-DE-001 (Build 92): paper18 DE — Gleiche Korrektur Bemerkung 4.3
- **Datei**: papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: Analog zu BUG-B4-P18-EN-001
- **Behoben**: Analog korrigiert

### BUG-B4-P18-DE-002 (Build 91): paper18 DE — Fallunterscheidung Satz 5.1
- **Datei**: papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: Beweis enthielt ursprünglich nur den Fall n≤N
- **Behoben**: Alle drei Fälle (n≤N, N<n≤2N, 2N<n≤3N) ergänzt

### BUG-B4-P19-EN-001 (Build 91): paper19 EN — Querverweise auf Paper 17
- **Datei**: papers/batch4/paper19_waring_problem.tex
- **Problem**: Falsche Referenzen Definition 5.1 und Lemma 5.3
- **Behoben**: Auf Definition 5.2 und Lemma 5.5 korrigiert

### BUG-B4-P19-DE-001 (Build 91): paper19 DE — Fehlendes Korollar Weyl-Schranke G(k)
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: Korollar G(k) ≤ (k-2)2^{k-1}+5 fehlte
- **Behoben**: Korollar mit Beweis in Abschnitt 4 ergänzt

### BUG-B4-P19-DE-002 (Build 91): paper19 DE — Wooley-Korollar ohne Beweis
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: G(k) ≤ k(log k + 4.20032) hatte keinen Beweis
- **Behoben**: Kurzen Beweis ergänzt

### BUG-B4-P19-DE-003 (Build 91): paper19 DE — Fehlender formaler Satz J_{s,k}
- **Datei**: papers/batch4/paper19_waringsches_problem_de.tex
- **Problem**: Nur Bemerkungen statt formalem Satz
- **Behoben**: Formaler Satz für J_{s,k}-Schranke aus EN adaptiert

### BUG-B4-P20-EN-001 (Build 91): paper20 EN — Satzfragment in Remark "Odd n"
- **Datei**: papers/batch4/paper20_goldbach_singular_series.tex
- **Problem**: "odd only for..." Satzfragment
- **Behoben**: Fragment entfernt, Remark vollständig formuliert

### BUG-B4-P20-EN-002 (Build 91): paper20 EN — Falsche Tabellenwerte S(n)
- **Datei**: papers/batch4/paper20_goldbach_singular_series.tex
- **Problem**: Tabellenwerte enthielten extra Faktor (p-1)/(p-2) für p∤n
- **Behoben**: Korrekte Werte n=4: 2C_2; n=6: 4C_2; n=8: 2C_2; n=10: 8C_2/3; n=12: 4C_2; n=30: 16C_2/3

### BUG-B4-P20-DE-001 (Build 91): paper20 DE — Identische Tabellenfehler
- **Datei**: papers/batch4/paper20_goldbach_singulaere_reihe_de.tex
- **Problem**: Analog zu BUG-B4-P20-EN-002
- **Behoben**: Korrekte Werte wie EN eingetragen

### BUG-010 (Build 91): paper2 DE — Redundanz in Abschnitt 3
- **Datei**: papers/reviewed/batch1/paper2_lehmer_drei_primfaktoren_de.tex
- **Problem**: Zwei vollständige Beweise ohne Kennzeichnung
- **Behoben**: `\medskip\noindent\textit{Alternativer Beweis:}` vor zweiten Beweis eingefügt

### BUG-011-Papers (Build 91): paper7 DE — Schwacher Zwischenschritt
- **Datei**: papers/reviewed/batch1/paper7_lehmer_kein_semiprim_de.tex
- **Problem**: `≥ 1·(p-1) - 1 = p-2 ≥ 1` war zu schwach
- **Behoben**: Auf `≥ 1 \cdot p - 1 = p-1 \ge 2 > 0` verschärft

### BUG-B1-P2-EN-001 (Build 116, Review Batch 1): paper2 EN — Acknowledgements fehlte
- **Datei**: papers/reviewed/batch1/paper2_lehmer_three_primes.tex
- **Problem**: DE-Version hatte \section*{Danksagung}, EN-Version fehlte entsprechende \section*{Acknowledgements}
- **Behoben**: Acknowledgements-Section vor Bibliographie eingefügt

### BUG-B1-P2-DE-002 (Build 116, Review Batch 1): paper2 DE — Falscher Divisor in Lemma lem:grenze_q
- **Datei**: papers/reviewed/batch1/paper2_lehmer_drei_primfaktoren_de.tex
- **Problem**: Schritt 3 benutzte "(p-1)(q-1) | (p+q-1)+k" statt des korrekten "A=(p-1)(q-1)/D | (p+q-1)+k" mit D=gcd(r-1,(p-1)(q-1))
- **Behoben**: D und A=(p-1)(q-1)/D eingeführt; Spezialfall D=1 für Schranke explizit benannt

### BUG-B1-P5-DE-001 (Build 116, Review Batch 1): paper5 DE — Grammatikfehler im Titel
- **Datei**: papers/reviewed/batch1/paper5_giuga_kein_semiprim_de.tex
- **Problem**: Titel lautete "Kein Semiprimes ist" (falsche Pluralform)
- **Behoben**: Titel korrigiert zu "Kein Semiprim ist"

### BUG-011 (Build 62): test_config.py – Falscher Autor-Wert
- **Datei**: tests/test_config.py
- **Problem**: TestMetadata::test_author_value prüfte auf 'Kurt Ingwer' statt 'Michael Fuhrmann'
- **Behoben**: Test korrigiert auf korrekten Autor 'Michael Fuhrmann'
