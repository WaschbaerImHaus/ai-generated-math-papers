# Gutachten: Papers 18–28 (Build 100)
**Datum:** 2026-03-11
**Gutachter:** Claude (Sonnet 4.6)
**Dateipfade:** `papers/batch4/` (18–20), `papers/batch5/` (21–24), `papers/batch6/` (25–28)

---

## Zusammenfassungstabelle

| Paper | Titel (EN) | Formal EN | Formal DE | Math. korrekt | Vollständig | EN/DE-Parität |
|-------|-----------|-----------|-----------|---------------|-------------|---------------|
| 18 | Vinogradov Three-Primes | OK | OK | OK | OK | DE hat Extra-§8 "Paritätsproblem"; `\newtheorem*{theorem*}` fehlt in DE |
| 19 | Waring's Problem | BUG: Wooley2013/2016 doppelter Subtitel | OK | OK | OK | DE fehlt Remark 4.1 (3 Schwellenbedingungen); G(k)-Tabelle in DE kürzer |
| 20 | Goldbach Singular Series | BUG: Granville1995 unused | OK | OK | OK | Granville-Bibitem-Titel unterschiedlich EN vs DE |
| 21 | Riemann Zeta Function | OK | OK | OK | OK | Sehr gute Parität, nahezu identisch |
| 22 | Nontrivial Zeros | OK | OK | OK | OK | DE Abstract: `7/8` fehlt (nur `O(log T)`); Haupttext korrekt |
| 23 | Explicit Formula PNT | OK | OK | OK | OK | DE Abstract hat Extra-Satz über GRH; kein Apostol1976 in DE (konsistent) |
| 24 | RH Approaches | OK | OK | OK | OK | DE §9 Physikverbindungen erheblich ausführlicher als EN |
| 25 | Elliptic Curves over Q | BUG: `\DeclareMathOperator{\gcd}{gcd}` Konflikt | OK | OK | OK | DE erklärt Neron–Ogg–Shafarevich expliziter |
| 26 | L-Function Elliptic | OK | BUG: TaylorWiles1995 Bd. 141 statt 142 | OK | OK | DE-Abschnittsreihenfolge verschieden; DE §6 Numerik detaillierter |
| 27 | BSD Conjecture | OK | OK | OK | OK | DE §7 Selmer-Gruppe ausführlicher erklärt; `\DeclareMathOperator{\Sha}{\Sha}` — potenz. Konflikt |
| 28 | Congruent Numbers BSD | OK | OK | OK | OK | DE fehlt `\Fp`-Makro (aber nicht benutzt); Beweis EN ausführlicher |

**Legende:** OK = kein Befund; BUG = Fehler gefunden

---

## Paper 18: Vinogradov Three-Primes Theorem

### EN: `papers/batch4/paper18_vinogradov_three_primes.tex`

**Formale Korrektheit:**
- Kein `\usepackage[utf8]{inputenc}` / `\usepackage[ngerman]{babel}` (korrekt für EN)
- `\newtheorem*{theorem*}{Vinogradov's Theorem}` korrekt definiert
- `\newcommand{\Q}{\mathbb{Q}}` vorhanden
- `proposition`, `corollary` Umgebungen definiert
- Alle 7 `\bibitem`-Schlüssel: `{Vinogradov1937}`, `{Hardy1923}`, `{IwaniecKowalski2004}`, `{Montgomery1994}`, `{Vaughan1997}`, `{Helfgott2013}`, `{Helfgott2014}` — alle im Text zitiert
- Keine undefinierten Referenzen

**Mathematische Schlüssigkeit:**
- Singuläre Reihe $\mathfrak{S}(N) = \prod_p (1 - (p-1)^{-2}) \cdot \ldots$ korrekt beschrieben
- Helfgott-Referenzen (2013 Preprint, 2014) für $N > 5$ vollständig abgedeckt
- Exponentielle Summe $S(\alpha) = \sum_{p \le N} e^{2\pi i p \alpha}$ korrekt

**Vollständigkeit:** 9 Abschnitte vorhanden (Intro, Exponential Sum, Major Arc, Singular Series, Singular Integral, Minor Arcs, Proof of Main Theorem, Helfgott, Conclusion). Keine abgebrochenen Sätze.

**EN/DE-Parität:** Siehe DE unten.

---

### DE: `papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex`

**Formale Korrektheit:**
- `\usepackage[utf8]{inputenc}` + `\usepackage[ngerman]{babel}` vorhanden
- `\newtheorem*{theorem*}` **fehlt** in DE — stattdessen wird `\textit{}` für Vinogradovs Satz verwendet. Dies ist kein Kompilierfehler, aber inkonsistent mit EN-Stil
- `\newcommand{\Q}{\mathbb{Q}}` **fehlt** in DE-Präambel — jedoch wird `\Q` im Text vermutlich nicht verwendet (Überprüfung: kein `\Q` in DE-Text festgestellt, daher kein Kompilierfehler)
- 7 `\bibitem`-Einträge, alle im Text zitiert

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** DE hat **10 Abschnitte**, darunter §8 "Das Paritätsproblem" (erklärt, warum der Kreis-Methode für das binäre Goldbach-Problem ein fundamentales Paritätsproblem im Weg steht). Dieser Abschnitt ist in EN **nicht vorhanden**.

**EN/DE-Parität:**
- DE §8 "Das Paritätsproblem" ist ein **Plus** gegenüber EN (kein Fehler, aber EN ist unvollständiger)
- DE verwendet `\textit{Vinogradovs Satz}` statt `\begin{theorem*}` — strukturelle Inkonsistenz
- Alle anderen Abschnitte inhaltlich äquivalent

---

## Paper 19: Waring's Problem

### EN: `papers/batch4/paper19_waring_problem.tex`

**Formale Korrektheit:**
- **BUG-19-EN-01:** `\bibitem{Wooley2013}` hat Titel: *"Vinogradov's mean value theorem via efficient congruencing, II"*. `\bibitem{Wooley2016}` hat ebenfalls Subtitel *"...via efficient congruencing, II"*. Dies ist ein bibliografischer Fehler: der 2016er Eintrag müsste einen anderen Untertitel tragen (vermutlich *"Efficient congruencing and the Vinogradov mean value"* oder ähnlich; die tatsächliche 2016-Publikation ist Wooley, *"The cubic case of the main conjecture in Vinogradov's mean value theorem"*, Adv. Math. 294 (2016)). Der Subtitel *"II"* gehört zum 2012/2013-Paper, nicht zum 2016-Paper.
- Alle anderen Bibitems korrekt
- Alle `\cite{...}` haben entsprechende `\bibitem`-Einträge

**Mathematische Schlüssigkeit:** Korrekt. G(k)-Tabelle mit k=2,3,4,5,6,7,8 korrekt.

**Vollständigkeit:**
- **Remark 4.1** "Three threshold conditions" erklärt die drei wichtigen Schwellenwerte für $s$: (a) $s \ge k^2$, (b) $s \ge k(k+1)/2$, (c) $s \ge 2k$. Vorhanden nur in EN, nicht in DE.

---

### DE: `papers/batch4/paper19_waringsches_problem_de.tex`

**Formale Korrektheit:** Korrekt. `inputenc` + `babel` vorhanden. Alle Bibitems zitiert.

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:**
- **BUG-19-DE-01 (Inhalt):** Remark 4.1 aus EN ("Drei Schwellenbedingungen") **fehlt** in DE. Die drei Schwellenwerte für $s$ werden in DE nicht gesondert als Bemerkung hervorgehoben.
- **BUG-19-DE-02 (Tabelle):** G(k)-Tabelle in DE hat nur 5 Zeilen (k=2,3,4,5, "große k"), EN hat 7 Zeilen (k=2,3,4,5,6,7,8). Fehlende Zeilen für k=6 ($G(6) \le 24$) und k=7 ($G(7) \le 33$) könnten für Vollständigkeit ergänzt werden.

**EN/DE-Parität:** Strukturell ähnlich, aber DE fehlt Remark 4.1 und hat kürzere G(k)-Tabelle.

---

## Paper 20: Goldbach's Conjecture — Singular Series

### EN: `papers/batch4/paper20_goldbach_singular_series.tex`

**Formale Korrektheit:**
- **BUG-20-EN-01:** `\bibitem{Granville1995}` ist definiert, wird aber **nirgendwo** im EN-Text via `\cite{Granville1995}` zitiert. Verwaistes Bibitem — LaTeX erzeugt keine Fehlermeldung, aber das Literaturverzeichnis enthält einen Eintrag ohne Textverweis.
- Alle anderen Bibitems zitiert

**Mathematische Schlüssigkeit:** Singuläre Reihe $\mathfrak{S}(N) = 2\prod_{p>2}(1 - (p-1)^{-2}) \prod_{p|N, p>2}((p-1)/(p-2))$ korrekt.

**Vollständigkeit:** Vollständig.

---

### DE: `papers/batch4/paper20_goldbach_singulaere_reihe_de.tex`

**Formale Korrektheit:**
- `\bibitem{Granville1995}` in DE vorhanden und im Text zitiert — **konsistent** innerhalb DE.
- Der Granville1995-Bibitem-Titel unterscheidet sich von EN:
  - EN: *"Some conjectures in additive number theory"*
  - DE: *"Refinements of Goldbach's conjecture, and the generalized Riemann hypothesis"*
  Diese zwei Einträge verweisen auf unterschiedliche Quellen. Es könnte sich um denselben Autor, aber unterschiedliche Papiere handeln (Granville hat mehrere relevante Papers). Inhaltlich muss geprüft werden, welche Quelle korrekt ist.

**BUG-20-DE-01:** Granville1995-Bibitem enthält einen anderen Papiertitel als EN — mögliche Quelleninkonsistenz zwischen EN und DE.

---

## Paper 21: Riemann Zeta Function

### EN: `papers/batch5/paper21_riemann_zeta_function_en.tex`

**Formale Korrektheit:** Vollständig korrekte Präambel. 6 Bibitems (`{Riemann1859}`, `{Titchmarsh1986}`, `{Edwards1974}`, `{Davenport1980}`, `{IwaniecKowalski2004}`, `{Montgomery1994}`), alle im Text zitiert.

**Mathematische Schlüssigkeit:**
- Euler-Produkt $\zeta(s) = \prod_p (1 - p^{-s})^{-1}$ für $\Re(s) > 1$ korrekt
- Funktionalgleichung $\xi(s) = \xi(1-s)$ korrekt formuliert
- Analytische Fortsetzung via Laurent-Entwicklung bei $s=1$ mit Residuum 1 korrekt
- Kritischer Streifen $0 < \Re(s) < 1$ korrekt beschrieben

**Vollständigkeit:** Alle Abschnitte vollständig. Keine abgebrochenen Sätze.

**EN/DE-Parität:** Nahezu identische Struktur und Inhalte — ausgezeichnete Parität.

---

### DE: `papers/batch5/paper21_riemann_zeta_function_de.tex`

**Formale Korrektheit:** `inputenc` + `babel` vorhanden. Alle 6 Bibitems zitiert.

**Mathematische Schlüssigkeit:** Korrekt, identisch zu EN.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:** Sehr gut. Inhaltlich und strukturell äquivalent zu EN. Alle Theoreme, Lemmata und Beweise in beiden Versionen vorhanden.

---

## Paper 22: Nontrivial Zeros of the Riemann Zeta Function

### EN: `papers/batch5/paper22_nontrivial_zeros_en.tex`

**Formale Korrektheit:** Korrekt. Alle Bibitems vorhanden und zitiert.

**Mathematische Schlüssigkeit:**
- $N(T) = \frac{T}{2\pi}\log\frac{T}{2\pi e} + \frac{7}{8} + O(\log T)$ korrekt im Abstract und im Haupttext
- GUE-Verbindung (Montgomery–Odlyzko) korrekt beschrieben

**Vollständigkeit:** Vollständig.

---

### DE: `papers/batch5/paper22_nontrivial_zeros_de.tex`

**Formale Korrektheit:** `inputenc` + `babel` vorhanden. Alle Bibitems korrekt.

**Mathematische Schlüssigkeit:**
- **BUG-22-DE-01 (Minor):** Im **Abstract** lautet die N(T)-Formel nur $\frac{T}{2\pi}\log\frac{T}{2\pi e} + O(\log T)$ — die explizite Konstante $\frac{7}{8}$ fehlt im Abstract. Im Haupttext (Satz 2.1) ist $\frac{7}{8}$ korrekt angegeben. Der Abstract ist daher geringfügig ungenau/unvollständig.

**Vollständigkeit:** Vollständig. Haupttext korrekt.

**EN/DE-Parität:** Gut. Einziger Unterschied: Abstract-Formel in DE lässt $\frac{7}{8}$ weg.

---

## Paper 23: Explicit Formula and Prime Number Theorem

### EN: `papers/batch5/paper23_explicit_formula_pnt_en.tex`

**Formale Korrektheit:**
- `\bibitem{Apostol1976}` vorhanden, im Beweis von Theorem 8.1 zitiert
- Alle anderen Bibitems korrekt

**Mathematische Schlüssigkeit:** Korrekt. Explizite Formel $\psi(x) = x - \sum_\rho x^\rho/\rho - \log 2\pi - \frac{1}{2}\log(1-x^{-2})$ korrekt formuliert.

**Vollständigkeit:** Vollständig.

---

### DE: `papers/batch5/paper23_explicit_formula_pnt_de.tex`

**Formale Korrektheit:**
- `\bibitem{Apostol1976}` **fehlt** in DE — aber `\cite{Apostol1976}` wird im DE-Text ebenfalls **nicht** verwendet. Daher kein Kompilierfehler: DE ist in sich konsistent.
- Alle im DE-Text verwendeten `\cite{...}` haben entsprechende `\bibitem`-Einträge.

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- DE Abstract enthält einen **Extra-Satz** über GRH-Implikationen (schärfere Fehlerterme unter GRH), der im EN-Abstract nicht vorhanden ist. Dies ist ein inhaltlicher Zusatz in DE, kein Fehler.
- EN Theorem 8.1 enthält eine Apostol1976-Zitation, die DE weggelassen hat. Inhaltlich äquivalent.

---

## Paper 24: Approaches to the Riemann Hypothesis

### EN: `papers/batch5/paper24_rh_approaches_en.tex`

**Formale Korrektheit:** Korrekt. Alle Bibitems vorhanden und zitiert.

**Mathematische Schlüssigkeit:** Korrekt. Verschiedene Ansätze (Spektraltheorie, Analytische Zahlentheorie, Algebraische Geometrie) korrekt beschrieben.

**Vollständigkeit:** §9 "Physics Connections" enthält 3 Stichpunkte (Zufallsmatrizen, Quantenchaos, GUE). Kompakt, aber vollständig.

---

### DE: `papers/batch5/paper24_rh_approaches_de.tex`

**Formale Korrektheit:** `inputenc` + `babel` vorhanden. Alle Bibitems korrekt.

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- DE §9 "Physikalische Verbindungen" ist **erheblich ausführlicher** als EN §9: Jeder Stichpunkt (Zufallsmatrizen GUE, Quantenchaos, Hamiltonoperator) wird mit mehreren Sätzen erläutert. Dies ist ein inhaltlicher Vorteil von DE gegenüber EN — kein Fehler, aber EN könnte entsprechend ergänzt werden.
- Alle anderen Abschnitte inhaltlich äquivalent.

---

## Paper 25: Elliptic Curves over Q

### EN: `papers/batch6/paper25_elliptic_curves_Q_en.tex`

**Formale Korrektheit:**
- **BUG-25-EN-01:** `\DeclareMathOperator{\gcd}{gcd}` ist definiert. Dies **überschreibt** den Standard-LaTeX-Operator `\gcd` (der bereits in `amsmath` definiert ist). Das führt zu einem "already defined" Warning oder stilistischen Inkonsistenzen. Besser wäre: Operator-Definition weglassen und Standard-`\gcd` verwenden.
- Alle Bibitems vorhanden und zitiert

**Mathematische Schlüssigkeit:** Korrekt. Nagell–Lutz, Mazur-Torsionssatz, Mordell–Weil korrekt formuliert.

**Vollständigkeit:** Vollständig.

---

### DE: `papers/batch6/paper25_elliptic_curves_Q_de.tex`

**Formale Korrektheit:**
- `\DeclareMathOperator{\gcd}{gcd}` **fehlt** in DE-Präambel — aber `\gcd` wird im DE-Text nicht benutzt (DE verwendet deutschen Fließtext statt `\gcd`-Makro). Kein Kompilierfehler.
- Alle Bibitems korrekt

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- DE beschreibt den Neron–Ogg–Shafarevich-Satz expliziter mit Galois-Darstellungssprache
- Inhaltlich sehr gute Parität; strukturell äquivalent

---

## Paper 26: L-Function of an Elliptic Curve

### EN: `papers/batch6/paper26_l_function_elliptic_en.tex`

**Formale Korrektheit:** Korrekt. Alle Bibitems vorhanden.
- `\bibitem{TaylorWiles1995}`: Band **142**, S. 553–572 — korrekt (Ann. of Math. 142 (1995))

**Mathematische Schlüssigkeit:** Korrekt. Frobenius-Spur $a_p = p+1-\#E(\mathbb{F}_p)$, Hasse-Schranke $|a_p| \le 2\sqrt{p}$, Euler-Produkt für L(E,s) korrekt.

**Vollständigkeit:** Vollständig. Abschnitte: Intro → Euler-Produkt → Hassesche Schranke → Analytische Fortsetzung → L-Funktion bei s=1 → Numerische Evidenz → Fazit.

---

### DE: `papers/batch6/paper26_l_function_elliptic_de.tex`

**Formale Korrektheit:**
- **BUG-26-DE-01:** `\bibitem{TaylorWiles1995}` in DE gibt Band **141** an. Korrekt ist Band **142** (Taylor–Wiles, Ann. of Math. (2) **142** (1995), 553–572). Dies ist ein bibliografischer Fehler.

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- DE-Abschnittsreihenfolge weicht von EN ab: DE hat §5 Modularität → §6 Numerik → §7 L-Funktion bei s=1 → §8 Fazit; EN hat §5 L-Funktion bei s=1 → §6 Numerik → §7 Fazit.
- DE §6 "Numerische Evidenz" ist inhaltlich detaillierter als EN (enthält Kolyvagin 1989-Zitation für Rang-0 und Rang-1-Fälle).
- Kein struktureller Fehler, aber Reihenfolge und Detailtiefe unterscheiden sich.

---

## Paper 27: The Birch–Swinnerton-Dyer Conjecture

### EN: `papers/batch6/paper27_bsd_conjecture_en.tex`

**Formale Korrektheit:**
- `\DeclareMathOperator{\Sha}{\Sha}` — Achtung: `\Sha` wird sowohl als Kommando (mit `\newcommand`) als auch als Operator-Name verwendet. Hier: `\DeclareMathOperator{\Sha}{\Sha}` definiert `\Sha` als Operator mit dem Symbol `\Sha` — das ist eine Selbstreferenz. Da `\Sha` in amsmath nicht vordefiniert ist, funktioniert dies technisch (der Operator-Name wird als `\Sha` gesetzt, was wiederum als Buchstaben-Sha-Symbol ausgegeben wird). **Potenzielle Warnung**: die Selbstreferenz ist zirkulär; besser wäre `\DeclareMathOperator{\Sha}{\text{\textesh}}` oder ein explizites Symbol. In der Praxis kompiliert es ohne Fehler.
- Alle 7 `\bibitem`-Schlüssel zitiert: `{BSD1965}`, `{CoatesWiles1977}`, `{Kolyvagin1990}`, `{GrossZagier1986}`, `{Silverman1986}`, `{Wiles1995}`, `{Milne2006}`

**Mathematische Schlüssigkeit:**
- Schwache BSD: $\text{ord}_{s=1} L(E,s) = r$ korrekt
- Starke BSD-Formel: $\lim_{s\to1} L(E,s)/(s-1)^r = \Omega_E \cdot |\Sha(E/\Q)| \cdot R_E \cdot \prod_p c_p / |E(\Q)_{\text{tors}}|^2$ korrekt
- BSD-Heuristikprodukt $B_E(X) \sim C_E(\log X)^r$ korrekt
- Kolyvagin-Beweisskizze inhaltlich korrekt

**Vollständigkeit:** 7 Abschnitte. Alle Beweise skizziert. Vollständig.

---

### DE: `papers/batch6/paper27_bsd_conjecture_de.tex`

**Formale Korrektheit:**
- `inputenc` + `babel` vorhanden
- `\DeclareMathOperator{\Sha}{\Sha}` — gleiche Anmerkung wie EN (technisch funktional)
- Zusätzlich: `\DeclareMathOperator{\rang}{rang}` für deutschen Text — korrekt und konsistent
- Alle 7 Bibitems vorhanden und zitiert
- Milne2006 URL in DE: `\texttt{https://...}` statt `\url{...}` — kein Fehler, aber stilistisch weniger robust (kein automatischer Zeilenumbruch bei URL)

**Mathematische Schlüssigkeit:** Korrekt.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- DE §7 "Abstieg und die Selmer-Gruppe" ist **wesentlich ausführlicher** als EN §7: DE erklärt die Kummer-Einbettung explizit, die exakte Sequenz wird vollständig ausgeführt, der Struktursatz $E(\Q)/2E(\Q) \cong (\Z/2\Z)^r \times E(\Q)_{\text{tors}}[2]$ wird hergeleitet. EN hat nur eine kurze Bemerkung.
- Dies ist ein inhaltlicher Vorteil von DE — kein Fehler, aber EN könnte von der ausführlicheren DE-Version profitieren.

---

## Paper 28: Congruent Numbers, Elliptic Curves, and BSD

### EN: `papers/batch6/paper28_congruent_numbers_bsd_en.tex`

**Formale Korrektheit:**
- `\Fp`-Makro (`\newcommand{\Fp}{\mathbb{F}_p}`) definiert und im BSD-Heuristik-Produktabsatz verwendet
- Alle 6 `\bibitem`-Schlüssel: `{Tunnell1983}`, `{Koblitz1984}`, `{Silverman1992}`, `{Silverman1986}`, `{Koblitz1993}`, `{Cassels1991}` — alle im Text zitiert
- `example`-Umgebung definiert und benutzt

**Mathematische Schlüssigkeit:**
- Äquivalenz kongruente Zahl ↔ Rang-≥-1-Punkt auf $E_n: y^2 = x^3 - n^2x$ vollständig bewiesen
- $x_0 = c^2/4$, $y_0 = c(a^2-b^2)/8$ Parametrisierung mit vollständiger Verifikation korrekt
- Torsion $E_n(\Q)_{\text{tors}} \cong \Z/2\Z \times \Z/2\Z$ korrekt bewiesen
- Tunnell-Satz (1983) mit allen 4 Zählfunktionen A(n), B(n), C(n), D(n) korrekt
- Verifikation für n=5,6,7 vollständig mit expliziten Rechnungen

**Vollständigkeit:** 7 Abschnitte. Alle Beweise ausgeführt. Sehr vollständig.

---

### DE: `papers/batch6/paper28_congruent_numbers_bsd_de.tex`

**Formale Korrektheit:**
- `inputenc` + `babel` vorhanden
- `\Fp`-Makro (`\newcommand{\Fp}{\mathbb{F}_p}`) **fehlt** in DE-Präambel — aber `\Fp` wird im DE-Text **nicht** verwendet (DE benutzt kein BSD-Heuristikprodukt mit `\Fp`). Kein Kompilierfehler.
- Alle 6 Bibitems vorhanden und zitiert (Reihenfolge leicht verschieden: `{Koblitz1993}` kommt in DE am Ende, in EN als vorletztes)

**Mathematische Schlüssigkeit:**
- Beweis der Äquivalenz korrekt, aber **kürzer** als EN: DE verweist auf Standardparametrisierung via `\cite{Koblitz1984}`, EN führt die Rechnung $c^4 - 16n^2 = (a^2-b^2)^2$ explizit aus.
- Torsion korrekt, Tunnell korrekt.
- **BUG-28-DE-01 (Minor):** In Bemerkung zu n=1,2,3 gibt DE an: "Tunnells Satz sagt unter BSD vorher, dass $C(2) \ne 2D(2)$... und $A(3) \ne 2B(3)$". Die Negation ist hier etwas irreführend formuliert: Tunnell sagt, dass wenn n **nicht** kongruent ist, dann gilt $A(n) \ne 2B(n)$ (Kontraposition). Streng genommen sagt die Negation des Kriteriums, dass n nicht kongruent ist. Die Formulierung "Tunnells Satz sagt... $C(2)\ne2D(2)$ vorher" ist logisch rückwärts — tatsächlich ist $C(2)=2\ne4=2D(2)$ für n=2 nach tatsächlicher Berechnung. Die DE-Formulierung könnte präziser sein, ist aber inhaltlich korrekt gemeint.

**Vollständigkeit:** Vollständig.

**EN/DE-Parität:**
- EN Beweis der Äquivalenz (Theorem 2.1) ist **ausführlicher** mit expliziter Rechnung
- DE verweist auf Literatur statt Rechnung ausführen — inhaltlich konsistent, aber weniger selbstständig
- Alle Hauptergebnisse in beiden Versionen vorhanden

---

## Gesamtbefunde und Priorisierung

### Kritische Fehler (Kompilierfehler oder schwerwiegende mathematische Fehler)
*Keine gefunden.*

### Signifikante Fehler (Bibliografische Fehler, die zur Verwirrung führen)

| ID | Paper | Datei | Beschreibung |
|----|-------|-------|--------------|
| BUG-19-EN-01 | Paper 19 | `paper19_waring_problem.tex` | `\bibitem{Wooley2016}` trägt denselben Subtitel wie `\bibitem{Wooley2013}` ("via efficient congruencing, II"). Der 2016er Titel muss korrigiert werden. |
| BUG-20-EN-01 | Paper 20 | `paper20_goldbach_singular_series.tex` | `\bibitem{Granville1995}` nie im EN-Text zitiert — verwaistes Bibitem. |
| BUG-20-DE-01 | Paper 20 | `paper20_goldbach_singulaere_reihe_de.tex` | Granville1995-Bibitem-Titel unterscheidet sich von EN (möglicherweise zwei verschiedene Quellen — Konsistenz prüfen). |
| BUG-26-DE-01 | Paper 26 | `paper26_l_function_elliptic_de.tex` | TaylorWiles1995: Band **141** in DE, korrekt ist **142**. |

### Kleinere Befunde (Minor)

| ID | Paper | Datei | Beschreibung |
|----|-------|-------|--------------|
| BUG-22-DE-01 | Paper 22 | `paper22_nontrivial_zeros_de.tex` | Abstract-Formel $N(T)$ lässt Konstante $7/8$ weg (Haupttext korrekt). |
| BUG-25-EN-01 | Paper 25 | `paper25_elliptic_curves_Q_en.tex` | `\DeclareMathOperator{\gcd}{gcd}` überschreibt Standard-amsmath-`\gcd`. |
| BUG-27-EN/DE-01 | Paper 27 | beide | `\DeclareMathOperator{\Sha}{\Sha}` ist Selbstreferenz — kompiliert, aber stilistisch suboptimal. |
| BUG-28-DE-01 | Paper 28 | `paper28_congruent_numbers_bsd_de.tex` | Formulierung "Tunnell sagt $C(2)\ne2D(2)$ vorher" logisch etwas unklar (Kontrapositionsrichtung). |

### Inhaltliche EN/DE-Unterschiede (kein Fehler, aber Hinweis)

| Paper | Richtung | Beschreibung |
|-------|----------|--------------|
| Paper 18 | DE > EN | DE §8 "Paritätsproblem" in EN fehlend |
| Paper 19 | EN > DE | EN Remark 4.1 ("Three threshold conditions") in DE fehlend; EN G(k)-Tabelle vollständiger |
| Paper 23 | DE > EN | DE Abstract hat Extra-Satz über GRH-Implikationen |
| Paper 24 | DE > EN | DE §9 Physik-Verbindungen erheblich ausführlicher |
| Paper 26 | DE ≠ EN | Abschnittsreihenfolge verschieden; DE Numerik-Abschnitt detaillierter |
| Paper 27 | DE > EN | DE §7 Selmer-Gruppe ausführlicher (Kummer-Einbettung, Struktursatz explizit) |
| Paper 28 | EN > DE | EN-Beweis von Theorem 2.1 mit expliziter Rechnung; DE kürzer und via Literatur |

---

## Empfehlungen

1. **Sofortiger Fix (BUG-26-DE-01):** In `paper26_l_function_elliptic_de.tex` den Band von TaylorWiles1995 von 141 auf **142** korrigieren.
2. **Sofortiger Fix (BUG-19-EN-01):** In `paper19_waring_problem.tex` den Titel von `\bibitem{Wooley2016}` auf den tatsächlichen 2016-Paper-Titel korrigieren (z.B. *"The cubic case of the main conjecture in Vinogradov's mean value theorem"*, Adv. Math. 294 (2016), 532–561).
3. **Empfohlen:** In `paper20_goldbach_singular_series.tex` (EN) entweder einen `\cite{Granville1995}`-Verweis im Text hinzufügen oder den verwaisten `\bibitem{Granville1995}` entfernen.
4. **Empfohlen:** Granville1995-Bibitem-Titel zwischen EN und DE harmonisieren (BUG-20-DE-01).
5. **Optional:** Abstract in `paper22_nontrivial_zeros_de.tex` um die $7/8$-Konstante ergänzen.
6. **Optional:** `\DeclareMathOperator{\gcd}{gcd}` in Paper 25 EN entfernen (Standard-`\gcd` verwenden).
7. **Informativ:** EN-Versionen von Papers 18, 23, 24, 27 könnten von den ausführlicheren DE-Abschnitten profitieren (Paritätsproblem, GRH-Hinweis, Physik-Verbindungen, Selmer-Gruppe).
