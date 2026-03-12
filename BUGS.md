# BUGS.md - specialist-maths

## Behobene Bugs (Build 153)

### BUG-B11-P40-EN/DE-001 (KRITISCH — BEHOBEN Build 153): Falsche "Giuga-Identität"
- **Dateien**: papers/batch11/paper40_giuga_four_primes_en.tex,
  papers/batch11/paper40_giuga_vier_primfaktoren_de.tex
- **Problem**: EN-Theorem 4.1 behauptete e₄ - e₃ + e₂ - e₁ + 1 = 0 als
  Äquivalenzbedingung. Das ist algebraisch unmöglich, da die linke Seite
  gleich φ(n) = ∏(pᵢ-1) ≥ 16 > 0 ist.
- **Behebung**: Theorem 4.1 in "Totient Identity" umstrukturiert.
  Korrekte Aussage: φ(n) = n - e₃ + e₂ - e₁ + 1 > 0.
  Das korrekte Kriterium ist das System (G1)-(G4). DE-Version war bereits
  korrekt formuliert; Abstract, Einleitung und Zusammenfassung beider
  Versionen bereinigt.

### BUG-B11-P40-EN/DE-002 (KRITISCH — BEHOBEN Build 153): Falsche Schranke p₄² ≤ p₁p₂p₃
- **Dateien**: papers/batch11/paper40_giuga_four_primes_en.tex (Theorem 5.3,
  Corollary 5.4), papers/batch11/paper40_giuga_vier_primfaktoren_de.tex
  (Satz 5.2, Korollar 5.3)
- **Problem**: Beide Dateien behaupteten p₄² ≤ p₁p₂p₃. Gegenbeispiel:
  n = 1722 = 2·3·7·41 ist eine Giuga-Zahl mit p₄ = 41, p₁p₂p₃ = 42,
  also p₄² = 1681 >> 42. Außerdem war der Beweis falsch:
  aus pᵢ | n/pᵢ - 1 folgt nur pᵢ ≤ n/pᵢ - 1, NICHT pᵢ(pᵢ+1) ≤ n/pᵢ.
- **Behebung**: Falsches Theorem + Korollar in beiden Dateien ersetzt.
  Korrekte lineare Schranken aus (G2)-(G4): p₄ ≤ p₁p₂p₃-1 etc.
  Neues Korollar: n < (p₁p₂p₃)². Erklärende Bemerkung über die
  Nichtexistenz der Quadratschranke hinzugefügt.

## Offene Bugs

### BUG-B1-P2-EN-HOCH-001: paper2 EN — Schritt 4 (p≥7) Abschätzung nicht rigoros
- **Datei**: papers/reviewed/batch1/paper2_lehmer_three_primes.tex
- **Problem**: Übergang von A≤2p+q-2 zur Schranke q<6(p-1)/(p-5) nur skizziert; A≥q/2 nicht rigoros bewiesen
- **Priorität**: HOCH — Beweislücke, kein fataler Fehler

### BUG-B1-P2-EN-MITTEL-002: paper2 EN — Fall p=5, q≥17 nur durch Computerscan
- **Datei**: papers/reviewed/batch1/paper2_lehmer_three_primes.tex
- **Problem**: Abstract behauptet "entirely elementary" — aber p=5, q≥17 wird nur durch Computerscan für q≤500 belegt; für q>500 kein Nachweis
- **Priorität**: MITTEL — Abstract-Claim falsch oder Beweis unvollständig

### BUG-B4-P18-EN/DE-003: paper18 EN+DE — Faktor 1/6 bei R₃(n) fragwürdig
- **Dateien**: papers/reviewed/batch4/paper18_vinogradov_three_primes.tex + paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: R₃(n) ~ (log n)³/6 · r₃(n) — falls r₃(n) geordnete Tripel zählt, gilt Faktor 1/6 NICHT (nur für ungeordnete Tripel)
- **Priorität**: MITTEL

### BUG-B4-P20-EN-TABLE: paper20 EN — r₂(100)-Vorhersagewert ≈8 falsch
- **Datei**: papers/reviewed/batch4/paper20_goldbach_singular_series.tex
- **Problem**: EN gibt r₂(100)≈8, korrekt ist ≈6 (DE hat den richtigen Wert)
- **Priorität**: GERING

### BUG-B5-P23-DE-001: paper23 DE — Abstract abgeschnitten
- **Datei**: papers/reviewed/batch5/paper23_explicit_formula_pnt_de.tex
- **Problem**: Letzter Satz über GRH und Dirichlet-L-Funktionen fehlt im DE-Abstract
- **Priorität**: MITTEL

### BUG-B6-P25-DE-003: paper25 DE — Fehlende Rank-Bemerkung
- **Datei**: papers/reviewed/batch6/paper25_elliptic_curves_Q_de.tex
- **Problem**: "It is not known whether ranks are unbounded" (EN Z. 300) fehlt in DE
- **Priorität**: GERING

### BUG-B6-P26-DE-TABELLE: paper26 DE — Fehlende aₚ-Tabelle
- **Datei**: papers/reviewed/batch6/paper26_l_function_elliptic_de.tex
- **Problem**: Tabelle für aₚ-Berechnung (y²=x³−x mod 5) fehlt in DE (Z. 109–115)
- **Priorität**: MITTEL

### BUG-B6-P27-DE-003: paper27 DE — Gekürzter Selmer-Gruppen-Abschnitt
- **Datei**: papers/reviewed/batch6/paper27_bsd_conjecture_de.tex
- **Problem**: Abschnitt 7 hat nur 8 Zeilen statt 14 (EN); exakte Folge 0→E(Q)/2E(Q)→Sel₂→Sha[2]→0 und Rang-Satz fehlen
- **Priorität**: MITTEL

## Behobene Bugs

### Build 124 (Audit 2026-03-12: Batch 1/2/4/6/9 Fixes)

### BUG-B1-P2-DE-MITTEL-001 (Build 124): paper2 DE — Algebraischer Fehler Lemma 4.1
- **Datei**: papers/reviewed/batch1/paper2_lehmer_drei_primfaktoren_de.tex
- **Problem**: Zeile 217: `p - 1 - p/(q+1) + (p-1)/(q+1)` algebraisch falsch; korrekt: `p - (p+1)/(q+1)`
- **Behoben**: Zwischenschritt korrigiert; Endergebnis k<p bleibt korrekt

### BUG-B2-P10-DE-GERING-001 (Build 124): paper10 DE — Interner Build-Kommentar
- **Datei**: papers/reviewed/batch2/paper10_wilson_quotient_de.tex
- **Problem**: `% BUG-B2-P10-DE-004: Schlussfolgerung korrigiert...` im Quelltext belassen
- **Behoben**: Kommentar entfernt

### BUG-B4-P19-DE-THM (Build 124): paper19 DE — Wooleys Haupttheorem als Bemerkung
- **Datei**: papers/reviewed/batch4/paper19_waringsches_problem_de.tex
- **Problem**: Wooleys Efficient-Congruencing-Theorem über G(k)-Schranke stand als `\begin{bemerkung}` statt `\begin{theorem}`
- **Behoben**: `\begin{bemerkung}` → `\begin{theorem}[Wooley, Efficient Congruencing]`

### BUG-B4-P20-DE-SECTION (Build 124): paper20 DE — Section 8 Einleitung ergänzt
- **Datei**: papers/reviewed/batch4/paper20_goldbach_singulaere_reihe_de.tex
- **Problem**: Einleitungsabsatz zu Section 8 "Der Goldbach-Komet" fehlte
- **Behoben**: Goldbach-Komet-Beschreibung aus EN übersetzt und eingefügt; Section-Titel präzisiert

### BUG-B5-P24-DE-ABSTRACT (Build 124): paper24 DE — Abstract unvollständig
- **Datei**: papers/reviewed/batch5/paper24_rh_approaches_de.tex
- **Problem**: Letzter Satz über GRH und Dirichlet-L-Funktionen fehlte
- **Behoben**: Abschließenden Satz eingefügt

### BUG-B5-P24-EN/DE-DELIGNE (Build 124): paper24 EN+DE — Fehlender Deligne1974-Bibitem
- **Dateien**: papers/reviewed/batch5/paper24_rh_approaches_en.tex + _de.tex
- **Problem**: `\cite{Deligne1974}` verwendet ohne passendes `\bibitem`
- **Behoben**: `\bibitem{Deligne1974}` in beide Versionen eingefügt

### BUG-B9-P37-DE-TYPO (Build 124): paper37 DE — Tippfehler "Kürzbarkei"
- **Datei**: papers/reviewed/batch9/paper37_algebraische_zahlentheorie_de.tex
- **Problem**: "Kürzbarkei" statt "Kürzbarkeit"
- **Behoben**: Tippfehler korrigiert

### Build 123 (Batch 10 Paper39 + Batch 9 Paper38 Bugfixes)

### BUG-B10-P39-EN/DE-LLC-HISTORY (Build 123): HOCH — Chronologisch unmögliche Aussage
- **Dateien**: papers/batch10/paper39_langlands_program_en.tex + paper39_langlands_programm_de.tex
- **Problem**: "Kutzko (1980) building on Henniart (1993)" — Kutzko kann nicht auf einer 13 Jahre späteren Arbeit aufgebaut haben.
- **Behoben**: EN+DE korrekt: GL₂(F) → Kutzko (1980); GL_n(F) → Henniart (2000) und Harris--Taylor (2001) unabhängig.

### BUG-B10-P39-DE-MISSING-SECTIONS (Build 123): HOCH — Fehlende Sektionen 7 und 8
- **Datei**: papers/batch10/paper39_langlands_programm_de.tex
- **Problem**: Sektion 7 (Langlands-Funktorialitätsvermutung) und Sektion 8 (p-adisches Langlands-Programm, Colmez 2010) fehlten vollständig.
- **Behoben**: Beide Sektionen aus EN übersetzt und eingefügt; Folgesektionen umnummeriert.

### BUG-B10-P39-DE-MISSING-BIBITEMS (Build 123): MITTEL — Fehlende Bibliographie-Einträge
- **Datei**: papers/batch10/paper39_langlands_programm_de.tex
- **Problem**: \bibitem{Henniart2000} und \bibitem{FrenkelBenZvi} fehlten in der deutschen Version.
- **Behoben**: Beide Einträge aus der EN-Version übernommen.

### BUG-B9-P38-EN/DE-WILES1990-STEP3 (Build 123): MITTEL — Falsche Beweismethode für Wiles 1990
- **Dateien**: papers/reviewed/batch9/paper38_iwasawa_theory_en.tex + paper38_iwasawa_theorie_de.tex
- **Problem**: Schritt 3 der Beweisskizze beschrieb Wiles 1990 mit "Euler system machinery (Kolyvagin, Rubin)"; Wiles 1990 verwendete aber Hecke-Algebren und Modulkurven. Der Euler-System-Ansatz stammt von Rubin (1991).
- **Behoben**: Schritt 3 in EN+DE korrigiert: Wiles 1990 → Hecke-Algebren + Modulkurven (Erweiterung Mazur–Wiles 1984); Euler-System-Ansatz → Rubin (1991) für imaginär-quadratische Körper.

### Build 122 (Batch 7 Collatz Bugfixes — Konditionale Formulierungen)

### BUG-B7-P29-CIRCULAR-PROP4 (Build 122): paper29 EN+DE — Proposition 4.1 zirkulär
- **Dateien**: papers/reviewed/batch7/paper29_collatz_conjecture_en.tex + paper29_collatz_vermutung_de.tex
- **Problem**: Beweis von Proposition 4.1 (Orbit-Dichten) setzte Collatz-Vermutung voraus ohne Konditional-Deklaration im Satzstatement
- **Behoben**: Proposition umbenannt zu "conditional on the Collatz Conjecture"; Beweis beginnt mit explizitem Vermutungshinweis

### BUG-B7-P30-CONCLUSION (Build 122): paper30 EN+DE — "S ist ergodisch" als Tatsache
- **Dateien**: papers/reviewed/batch7/paper30_collatz_p_adic_en.tex + paper30_collatz_p_adisch_de.tex
- **Problem**: Conclusion: "but S is ergodic" (EN) / "aber S ist ergodisch" (DE) als bewiesene Aussage; im Textkörper korrekt als Vermutung
- **Behoben**: → "is conjectured to be ergodic" (EN) / "wird vermutet, ergodisch zu sein" (DE)

### BUG-B7-P32-CONCLUSION (Build 122/150): paper32 EN+DE — "S ergodisch und mischend" als Tatsache
- **Dateien**: papers/reviewed/batch7/paper32_collatz_ergodic_en.tex + paper32_collatz_ergodisch_de.tex
- **Problem**: Conclusion: "S is ergodic and mixing..." als bewiesene Aussage; alle diese Eigenschaften sind unbewiesen
- **Behoben (Build 150)**: Konjunktur-Formulierung + "numerical evidence strongly supports this" hinzugefügt (EN+DE)

### BUG-B7-P31-FURSTENBERG-REF (Build 122): paper31 EN+DE — Falsche Furstenberg-Referenz
- **Dateien**: papers/reviewed/batch7/paper31_tao_probabilistic_collatz_en.tex + paper31_tao_probabilistischer_ansatz_de.tex
- **Problem**: Lemma 5.3 zitiert Furstenberg (1977) für Korrespondenzprinzip im Collatz-Kontext; Furstenberg 1977 ist der Beweis von Szemerédi's Theorem, ohne Collatz-Bezug
- **Behoben**: Lemma-Titel zu "Korrespondenzprinzip im Collatz-Kontext" umbenannt; `\cite{Furstenberg1977}` → `\cite{Tao2022}`

### Build 120 (Batch 9 Review)

### BUG-B9-P37-EN-001 (Build 120): paper37 EN — Theorem 6.1 widersprüchliche Hypothese
- **Datei**: papers/reviewed/batch9/paper37_algebraic_number_theory_en.tex
- **Problem**: "An odd prime $p$ not dividing $D$" — aber Fall 2 deckt "$p \mid D$ (ramified)" ab. Direkter Widerspruch.
- **Behoben**: "not dividing $D$" aus der Hypothese entfernt → gilt jetzt für alle ungeraden Primzahlen $p$.

### BUG-B9-P37-DE-001 (Build 120): paper37 DE — gleicher Fehler
- **Datei**: papers/reviewed/batch9/paper37_algebraische_zahlentheorie_de.tex
- **Problem**: "Eine ungerade Primzahl $p \nmid D$" aber Fall "verzweigt" setzt $p \mid D$ voraus.
- **Behoben**: "$p \nmid D$" aus der Hypothese entfernt.

### BUG-B9-P38-EN-001 (Build 120): paper38 EN — Kubota-Leopoldt Interpolation zu restriktiv
- **Datei**: papers/reviewed/batch9/paper38_iwasawa_theory_en.tex
- **Problem**: Interpolationsformel nur für "$n \equiv 1 \pmod{p-1}$" behauptet; Standardsatz (Washington Theorem 5.11) gilt für alle $n \geq 1$.
- **Behoben**: "$n \equiv 1 \pmod{p-1}$" → "$n \geq 1$".

### BUG-B9-P38-DE-001 (Build 120): paper38 DE — gleicher Fehler
- **Datei**: papers/reviewed/batch9/paper38_iwasawa_theorie_de.tex
- **Problem**: "für alle positiven ganzen Zahlen $n \equiv 1 \pmod{p-1}$" — zu restriktiv.
- **Behoben**: Kongruenzbedingung entfernt → "für alle positiven ganzen Zahlen $n \geq 1$".

## Behobene Bugs

### Build 118 (Mittel-Bugs)

### BUG-B4-P18-DE-002 (Build 118): paper18 DE — Simplex-Fallunterscheidung unvollständig (NEU)
- **Datei**: papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex
- **Problem**: Simplex-Fallunterscheidung war unvollständig; Substitution u_i=t_i-2 fehlte
- **Behoben**: Substitution u_i=t_i-2 eingefügt, zeigt n≤N-Formel exakt

### BUG-B4-P19-EN-001 (Build 118): paper19 EN — Querverweis "Paper 17, Theorem 6.1"
- **Datei**: papers/batch4/paper19_waring_problem.tex
- **Problem**: Querverweis "Paper 17, Theorem 6.1" — FALSE POSITIVE: Theorem 6.1 ist korrekt
- **Behoben**: Fußnote zur Klarstellung eingefügt (kein inhaltlicher Fehler)

### BUG-B4-P20-EN-001 (Build 118): paper20 EN — GRH-Beweisskizze ohne q-Abhängigkeit (NEU)
- **Datei**: papers/batch4/paper20_goldbach_singular_series.tex
- **Problem**: GRH-Beweisskizze: Fehlerterm ohne q-Abhängigkeit angegeben
- **Behoben**: O(√N(log N)²/φ(q)) pro Hauptbogen explizit ergänzt

### BUG-B5-P23-EN/DE-001 (Build 118): paper23 EN+DE — Beitrag bei s=0 falsch als "Pol von Γ"
- **Datei**: papers/batch5/paper23_explicit_formula_pnt_en.tex + _de.tex
- **Problem**: Beitrag bei s=0 als "Pol von Γ" beschrieben — falsch, da ζ(0)=-1/2≠0
- **Behoben**: "Residuum des 1/s-Faktors bei s=0" (korrekter Ursprung)

### BUG-B5-P23-EN/DE-002 (Build 118): paper23 EN — PNT-Beweis "fourth-order pole of ζ^4" falsch
- **Datei**: papers/batch5/paper23_explicit_formula_pnt_en.tex
- **Problem**: "fourth-order pole of ζ^4" falsch — ζ hat 4-fache Nullstelle, keinen Pol
- **Behoben**: "zero of order 4m, contributing (3-4m)·∞ → -∞"

### BUG-B7-P30-EN/DE-001 (Build 118): paper30 EN+DE — -1 als Fixpunkt vs. 2-Zyklus
- **Datei**: papers/batch7/paper30_collatz_p_adic_en.tex + paper30_collatz_p_adisch_de.tex
- **Problem**: -1 als Fixpunkt von S vs. 2-Zyklus unter T unklar dargestellt
- **Behoben**: Erklärende Remark-Box eingefügt

### BUG-B7-P30-EN/DE-003 (Build 118): paper30 EN+DE — "Orbits aller Perioden" unbewiesen
- **Datei**: papers/batch7/paper30_collatz_p_adic_en.tex + paper30_collatz_p_adisch_de.tex
- **Problem**: Aussage "Orbits aller Perioden" ohne Beweis
- **Behoben**: Beweisskizze (CRT+Hensel-Lift) + Lagarias-Zitat ergänzt

### BUG-B7-P30-EN/DE-004 (Build 118): paper30 EN+DE — Ergodizität von S unbewiesen
- **Datei**: papers/batch7/paper30_collatz_p_adic_en.tex + paper30_collatz_p_adisch_de.tex
- **Problem**: Ergodizität von S als Theorem behauptet ohne Beweis
- **Behoben**: Theorem → Conjecture umklassifiziert + Status-Remark

### BUG-B7-P32-EN/DE-001 (Build 118): paper32 EN+DE — ν-Eindeutigkeit unbegründet
- **Datei**: papers/batch7/paper32_collatz_ergodic_en.tex + paper32_collatz_ergodisch_de.tex
- **Problem**: Eindeutigkeit des invarianten Maßes ν unbegründet
- **Behoben**: Krylov-Bogolyubov statt "unique ergodicity of shift" als Begründung

### BUG-B7-P32-EN/DE-002 (Build 118): paper32 EN+DE — Mixing-Proof zirkulär
- **Datei**: papers/batch7/paper32_collatz_ergodic_en.tex + paper32_collatz_ergodisch_de.tex
- **Problem**: Mixing-Beweis enthielt zirkuläre Argumentation
- **Behoben**: Bedingter Beweis, Abhängigkeitskette explizit dargestellt

### BUG-B7-P32-EN/DE-003 (Build 118): paper32 EN+DE — Spektralradius Transfer-Operator nicht rigoros
- **Datei**: papers/batch7/paper32_collatz_ergodic_en.tex + paper32_collatz_ergodisch_de.tex
- **Problem**: Argument zum Spektralradius des Transfer-Operators nicht rigoros
- **Behoben**: Als informal markiert + erklärende Remark

### Build 117 (Kritische/Hohe Bugs)

### BUG-B5-P24-EN/DE-002 (Build 117): paper24 EN+DE — Weil "bewies" Weil-Vermutungen falsch
- **Datei**: papers/batch5/paper24_rh_approaches_en.tex + paper24_rh_approaches_de.tex
- **Problem**: "Weil (1940) proved Weil conjectures" — falsch. Weil VERMUTETE (1949), bewies nur Kurvenfall (1940); Deligne BEWIES allgemeinen Fall (1974)
- **Behoben**: Historische Attribuierung korrigiert: Vermutung 1949, Kurvenfall-Beweis 1940 (Weil), allgemeiner Beweis 1974 (Deligne)

### BUG-B7-P30-EN/DE-002 (Build 117): paper30 EN+DE — Syracuse Repr. Theorem c_{k+1} ∈ Z₂ unbewiesen
- **Datei**: papers/batch7/paper30_collatz_p_adic_en.tex + paper30_collatz_p_adisch_de.tex
- **Problem**: Im Induktionsbeweis des Syracuse Repr. Theorem: c_{k+1} ∈ Z₂ nicht explizit nachgewiesen
- **Behoben**: 2-adisches Bewertungsargument ergänzt

### BUG-B7-P31-EN/DE-001 (Build 117): paper31 EN+DE — Key Estimate zu informell
- **Datei**: papers/batch7/paper31_tao_probabilistic_collatz_en.tex + paper31_tao_probabilistischer_ansatz_de.tex
- **Problem**: Proof Sketch zu informell; ξ ∉ Z/log₂3 irreführend (falsche Formulierung der irrationalen Bedingung)
- **Behoben**: Als "[informal — kein rigoroser Beweis]" deklariert; Bedingung auf ‖ξ‖_{R/Z}>0 korrigiert

### BUG-B7-P32-EN/DE-004 (Build 117): paper32 EN+DE — Hausdorff-Dim Theorem falsch klassifiziert
- **Datei**: papers/batch7/paper32_collatz_ergodic_en.tex + paper32_collatz_ergodisch_de.tex
- **Problem**: Hausdorff-Dim Theorem: logdens 0 ⟹ dim<1 gilt NICHT allgemein — Satz zu stark
- **Behoben**: Theorem → Conjecture umklassifiziert + erklärende Remark

### BUG-B8-P34-STATUS (Build 117): paper34 EN+DE — abc-Vermutung Status unklar
- **Datei**: papers/batch8/paper34_abc_conjecture_en.tex + paper34_abc_vermutung_de.tex
- **Problem**: abc-Vermutung als "disputed" statt klar als "OFFEN" dargestellt; Scholze-Stix-Einwand unzureichend erläutert
- **Behoben**: Explizite Remark-Box eingefügt: abc OFFEN, Scholze-Stix von Mehrheit der führenden Zahlentheoretiker anerkannt

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
