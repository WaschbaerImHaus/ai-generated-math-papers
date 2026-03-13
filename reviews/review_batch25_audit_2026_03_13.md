# Review Batch 25 -- Audit 2026-03-13

## Papers geprüft
| Paper | Thema | EN | DE |
|-------|-------|----|----|
| 96 | Grothendieck Standard-Vermutungen | paper96_grothendieck_standard_conjectures_en.tex | paper96_grothendieck_standard_vermutungen_de.tex |
| 97 | Kontsevich-Integral | paper97_kontsevich_integral_en.tex | paper97_kontsevich_integral_de.tex |

---

## Prüfkriterium 1: Theorem vs. Conjecture

### Paper 96 (Grothendieck)
- **Hard Lefschetz über C (Hodge 1950)**: Korrekt als Theorem. BEWIESEN. OK.
- **Conjecture A (Lefschetz)**: Korrekt als Conjecture. OFFEN (allgemein). OK.
- **Conjecture B (Hodge Standard)**: Korrekt als Conjecture. OFFEN (allgemein). OK.
- **Conjecture C (Künneth)**: Korrekt als Conjecture. OFFEN (allgemein). OK.
- **Conjecture D (Num=Hom)**: Korrekt als Conjecture. OFFEN (allgemein). OK.
- **Kleiman D => A**: Korrekt als Theorem. BEWIESEN. OK.
- **A => C**: Korrekt als Theorem (Grothendieck). BEWIESEN. OK.
- **B => A**: Korrekt als Proposition. BEWIESEN. OK.
- **B+A => D**: Korrekt als Theorem (Grothendieck). BEWIESEN. OK.
- **Jannsen 1992 (Halbeinfachheit)**: Korrekt als Theorem. BEWIESEN (ohne Standardvermutungen!). OK.
- **Andre 1996 (motivierte Zyklen)**: Korrekt als Theorem. BEWIESEN. OK.
- **Weil-Vermutungen (Deligne 1974)**: Korrekt als bewiesenes Theorem referenziert. OK.
- **Hodge-Vermutung**: Korrekt als Conjecture. OFFEN. OK.
- **Tate-Vermutung**: Korrekt als Conjecture. OFFEN. OK.
- **Murre-Vermutungen**: Korrekt als Conjecture. OFFEN. OK.

**Spezialfall-Status korrekt:**
- Conjecture A für Kurven: trivial. OK.
- Conjecture A für abelsche Varietäten: Lieberman. OK.
- Conjecture B über C: Hodge-Theorie. OK.
- Conjecture D für Divisoren: Neron-Severi. OK.
- Künneth-Standard-Vermutung: Im Text wird C nicht als "bewiesen" deklariert, sondern korrekt als Conjecture behandelt. Die Aussage im Aufgaben-Status war irreführend -- tatsächlich ist C allgemein OFFEN, und die Hodge-Zerlegung über C gibt Künneth-Komponenten, die aber nicht als algebraisch bekannt sind. Paper 96 handhabt das korrekt.

### Paper 97 (Kontsevich)
- **Kontsevich-Integral (Konstruktion)**: Korrekt als Definition + Invarianz-Theorem. OK.
- **Universalität des Kontsevich-Integrals**: Korrekt als Theorem. BEWIESEN. OK.
- **Birman-Lin 1993**: Korrekt als Theorem. BEWIESEN. OK.
- **Vassiliev-Kontsevich (weight systems)**: Korrekt als Theorem. BEWIESEN. OK.
- **Lie-algebraische Gewichtssysteme (Bar-Natan 1995)**: Korrekt als Theorem. BEWIESEN. OK.
- **Drinfeld Assoziator (1990)**: Korrekt als Theorem. BEWIESEN (Konstruktion). OK.
- **Le-Murakami Konvergenz (1996)**: Korrekt als Theorem. BEWIESEN. OK.
- **Wheeling Theorem (Bar-Natan-Le-Thurston 2003)**: Korrekt als Theorem. BEWIESEN. OK.
- **LMO-Invariante (1998)**: Korrekt als Theorem. BEWIESEN. OK.
- **Reshetikhin-Turaev (1991)**: Korrekt als Theorem. BEWIESEN. OK.
- **Completeness Conjecture (Kontsevich trennt Knoten)**: Korrekt als Conjecture. OFFEN. OK.
- **Vassiliev erkennt Unknoten**: Korrekt als Conjecture. OFFEN. OK.
- **Alle Gewichtssysteme Lie-algebraisch**: Korrekt als Conjecture. OFFEN. OK.
- **Volume Conjecture (Kashaev-Murakami)**: Korrekt als Conjecture/Remark. OFFEN. OK.
- **Jones Unknot Conjecture**: Nicht explizit als separate Conjecture, aber die Completeness Conjecture subsumiert sie. OK.

---

## Prüfkriterium 2: Korrekte Formeln

### Paper 96
- Weil-Kohomologie-Axiome: Poincare-Dualität, Künneth, Zyklenklasse. Alle korrekt. OK.
- Lefschetz-Operator: $\mathbf{L}^{d-i}: H^i(X) \to H^{2d-i}(X)$. OK.
- Hodge-Riemann-Bilinearform: $Q(\alpha,\beta) = (-1)^{i(i-1)/2} \int_X \alpha \cup \beta \cup L^{d-i}$. OK.
- Implikationsdiagramm: $B \Rightarrow A \Rightarrow C$, $D \Rightarrow A$, $A+B \Rightarrow D$. OK.
- Zetafunktion: $Z(X,T) = \exp(\sum |X(\mathbb{F}_{p^n})| T^n/n)$. OK.
- Tate: Zyklenklasse nach Galois-Fixpunkte. OK.

### Paper 97
- Kontsevich-Integral: $Z(K) = \sum_{m=0}^\infty \frac{1}{(2\pi i)^m} \int \sum (-1)^{\downarrow P} D_P \bigwedge \frac{dz_j - dz_j'}{z_j - z_j'}$. OK.
- 4T-Relation: $D_1 - D_2 + D_3 - D_4 = 0$. OK.
- KZ-Gleichung: $\frac{dG}{dz} = (\frac{A}{z} + \frac{B}{z-1})G$. OK.
- Volume Conjecture: $\lim \frac{2\pi}{N}\log|J_N(K;e^{2\pi i/N})| = \mathrm{vol}(S^3\setminus K)$. OK.
- Wheeling: $\Omega = \sum \frac{b_{2n}}{n!}\omega_{2n}$. OK.
- LMO: Korrekt referenziert (komplexe Formel). OK.
- Chern-Simons: $Z_{CS}(M,K) = \int e^{ik/4\pi \int \tr(A\wedge dA + \frac{2}{3}A^3)} W_R(K) \mathcal{D}A$. OK.

---

## Prüfkriterium 3: LaTeX-Syntax

### Paper 96
- EN und DE: Keine Syntaxfehler. OK.

### Paper 97
- **KRITISCHER LaTeX-FEHLER**: Doppelte Definition von `\deg` (Zeilen 38-39 in EN): `\newcommand{\deg}` und `\DeclareMathOperator{\deg}` -- ergibt LaTeX-Kompilierungsfehler. **KORRIGIERT** (zweite Definition in `\degree` umbenannt, erste entfernt).
- DE: Keine Syntaxfehler. OK.

---

## Prüfkriterium 4: Bibliographie

### Paper 96
18 (EN) / 18 (DE) Referenzen. Grothendieck 1968, Kleiman 1968, Deligne 1974+1980, Jannsen 1992, Andre 1996, Tate 1994, Milne 2002, Lieberman 1968, Mumford 1970, Bloch-Srinivas 1983, Lefschetz 1924, Voevodsky 2000, Bondal-Orlov 2001, Beilinson 1985, Kleiman 1994, Murre 1993, Scholl 1994 -- alle vorhanden und korrekt. OK.

### Paper 97
15 (EN) / 15 (DE) Referenzen. Vassiliev 1990, Kontsevich 1993, Bar-Natan 1995+2003, Jones 1985, Witten 1989, Drinfeld 1990, Birman-Lin 1993, LMO 1998, Reshetikhin-Turaev 1991, Khovanov 2000, Le-Murakami 1996, Chmutov-Duzhin-Mostovoy 2012, Ohtsuki 2002, Murakami-Murakami 2001 -- alle vorhanden. OK.

---

## Gefundene und korrigierte Fehler

| # | Paper | Schwere | Beschreibung | Status |
|---|-------|---------|--------------|--------|
| 1 | 97 EN | KRITISCH | Doppelte `\DeclareMathOperator{\deg}` / `\newcommand{\deg}` -- LaTeX-Kompilierungsfehler | KORRIGIERT |
| 2 | 96 DE | MITTEL | Tippfehler "Halbeinachheitssatz" statt "Halbeinfachheitssatz" (2x) | KORRIGIERT |

---

## Gesamtbewertung

Alle 4 Papers in Batch 25 sind nach Korrektur **DRUCKREIF**.

- Grothendieck Standard-Vermutungen: Alle vier (A, B, C, D) korrekt als OFFEN markiert.
- Künneth-Standard-Vermutung (C): Korrekt als allgemein OFFEN behandelt. Die Hodge-Zerlegung gibt Künneth-Komponenten über C, aber deren Algebraizität ist nicht bekannt.
- Jannsen Halbeinfachheit: Korrekt als BEWIESEN (ohne Standardvermutungen).
- Kontsevich-Integral: Korrekt als Konstruktion/Definition, Universalität BEWIESEN.
- Completeness Conjecture: Korrekt als OFFEN.
- Jones Unknot Conjecture: Korrekt implizit durch Completeness Conjecture abgedeckt.
- Volume Conjecture: Korrekt als OFFEN.
