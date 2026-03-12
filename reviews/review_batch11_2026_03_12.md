# Batch 11 Audit — 2026-03-12

Geprüfte Dateien: 8 LaTeX-Papers in `/papers/batch11/`
Prüfer: Claude Sonnet 4.6 (Specialist Mathematics Project)
Build bei Audit: 133

---

## Paper 40 (Giuga 4-Prim-Fall)

### EN — `paper40_giuga_four_primes_en.tex`

**A. Mathematische Korrektheit**
- [x] Giugas Vermutung korrekt als `Conjecture` deklariert (Abschnitt 8, `\begin{conjecture}[Giuga's Conjecture]`)
- [x] Vier-Prim-Fall ebenfalls als `Conjecture` (nicht als Theorem)
- [x] Giuga-Identität für vier Faktoren korrekt: `n - e3 + e2 - e1 + 1 = φ(n)` — numerisch verifiziert mit n=1722
- [x] Beweis der Identität korrekt: `φ(n) = prod(p_i - 1) = f(1) = n - e3 + e2 - e1 + 1`
- [x] Keine-3-Prim-Pseudoprimzahl als Theorem mit Quellenangabe `[FuhrmannPaper1]`
- [x] Beispiel n=1722=2·3·7·41: alle vier Giuga-Bedingungen korrekt nachgerechnet
- [x] Schranke p4² ≤ p1·p2·p3 korrekt bewiesen
- [x] Computerresultat n ≤ 10^15 korrekt als „Computational Non-Existence Result" gekennzeichnet

**Gefundener Bug (KRITISCH — behoben):**
- **Build-Nummer in `\address`**: Stand `Build 131` statt `Build 133` → **korrigiert auf Build 133**

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}` ← nach Korrektur
- [x] Alle Theorem-Environments definiert (theorem, lemma, corollary, proposition, definition, example, remark, conjecture)
- [x] `\tableofcontents` nach `\maketitle`
- [x] 10 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Englisch
- [x] Mathematische Notation konsistent

**D. Vollständigkeit**
- [x] 10 Sections (Intro, Background, Equation System, Giuga Identity, Arithmetic Constraints, Parity, Bounds, Computation, Open Problems, Conclusion)
- [x] Konkrete Beispiele (n=1722, Näherungsbeispiele)
- [x] Conclusion-Section vorhanden

**Status: DRUCKREIF** (nach Build-Nummern-Fix)

---

### DE — `paper40_giuga_vier_primfaktoren_de.tex`

**A. Mathematische Korrektheit**
- [x] Giugas Vermutung korrekt als `Vermutung` (`\begin{conjecture}`)
- [x] Vier-Prim-Fall ebenfalls als `Vermutung`
- [x] Totient-Identität korrekt: `φ(n) = n - e3 + e2 - e1 + 1` — Beweis über f(1) korrekt
- [x] Schranken p4 ≤ p1p2p3-1 und p4² ≤ p1p2p3 korrekt bewiesen
- [x] Beispiel n=1722 korrekt verifiziert
- [x] Computerresultat als Satz (nicht als bewiesene Aussage), nur rechnerisch
- [x] Paritätsanalyse: korrekt — kein automatischer Widerspruch im 4-Prim-Fall

**Gefundener Bug (KRITISCH — behoben):**
- **Build-Nummer in `\address`**: Stand `Build 131` statt `Build 133` → **korrigiert auf Build 133**

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}` ← nach Korrektur
- [x] Theorem-Environments auf Deutsch: Satz, Lemma, Korollar, Proposition, Definition, Beispiel, Bemerkung, Vermutung
- [x] `\tableofcontents` nach `\maketitle`
- [x] 10 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Deutsch
- [x] Keine englischen Theorem-Überschriften (Satz statt Theorem, Korollar statt Corollary, etc.)
- [x] `\DeclareMathOperator{\ggT}{ggT}` statt `\lcm` auf Deutsch

**D. Vollständigkeit**
- [x] 9 Sections (Einleitung, Grundlagen, System, Identität, Schranken, Paritätsanalyse fehlt als eigener Abschnitt… liegt aber in Abschnitt 7 subsections)
- [x] Konkrete Beispiele vorhanden
- [x] Zusammenfassung-Section vorhanden

**Hinweis (geringfügig):** Paper EN hat 10 Sections, DE hat 9 — die Paritäts-Sektion (sec:parity EN) fehlt als eigener Abschnitt im DE; sie ist in die Schranken-Sektion integriert (sec:paritaet ist ein eigener Abschnitt, also doch 9 Abschnitte). Minimal-Abweichung, kein Bug.

**Status: DRUCKREIF** (nach Build-Nummern-Fix)

---

## Paper 41 (Erdős-Straus-Vermutung)

### EN — `paper41_erdos_straus_en.tex`

**A. Mathematische Korrektheit**
- [x] Erdős-Straus als `Conjecture` deklariert (zweifach: conj:ES und conj:ES2)
- [x] Konjektur-Formulierung korrekt: `4/n = 1/a + 1/b + 1/c` für `n ≥ 2` — die Bedingung lautet n ≥ 2 (nicht n ≥ 5, wie in Audit-Aufgabe als möglicher Fehler angedeutet). Das Paper verwendet `n ≥ 2` korrekt; die Primzahlreduktion auf `p ≥ 5` ist eine separate Proposition, nicht die Aussage der Vermutung selbst. **Kein Fehler.**
- [x] Diophantische Gleichungsform: `4abc = n(ab+ac+bc)` korrekt
- [x] Konstruktionen Typ I (n≡0 mod 4), II (n≡2 mod 3), III (n≡3 mod 4) korrekt bewiesen
- [x] Vaughan 1970: O(N exp(-c√log N)) — korrekt angegeben
- [x] Fibonacci-Sylvester-Algorithmus als Theorem (bewiesen)
- [x] Keine unbewiesen-als-bewiesen deklarierten Aussagen gefunden
- [x] Beispielrechnungen (4/5, 4/7, 4/11, 4/13, 4/17, 4/41) nachgeprüft: alle korrekt

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Alle Theorem-Environments definiert inkl. `conjecture`
- [x] `\tableofcontents` nach `\maketitle`
- [x] 10 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Englisch

**D. Vollständigkeit**
- [x] 10 Sections (Introduction, Egyptian Fractions, Formal Statement, Equivalent Formulations, Constructive Approaches, Modular Reductions, Density Results, Computational Verification, Connections, Open Problems) + starred Conclusion
- [x] Konkrete Beispiele mit Zahlen
- [x] Conclusion-Section vorhanden (starred)

**Gefundene Bugs:** Keine

**Status: DRUCKREIF**

---

### DE — `paper41_erdos_straus_de.tex`

**A. Mathematische Korrektheit**
- [x] Erdős-Straus als `Vermutung` (conjecture) deklariert
- [x] Vermutungs-Formulierung korrekt: n ≥ 2
- [x] Konstruktionen Typ I–III korrekt mit deutschen Beweisen
- [x] Beispielrechnungen korrekt (4/5=1/5+1/2+1/10, 4/7=1/2+1/14, etc.)
- [x] `[ngerman]{babel}` Paket vorhanden — einziges DE-Paper mit explizitem babel (kein Bug, aber Inkonsistenz zu anderen DE-Papers)

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Theorem-Environments auf Deutsch: Satz, Lemma, Korollar, Proposition, Definition, Beispiel, Bemerkung, Vermutung
- [x] `\tableofcontents` nach `\maketitle`
- [x] 10 `\bibitem`-Einträge (10 inkl. RheindPapyrus — EN hat 10 ohne diesen)

**C. Sprachliche Qualität**
- [x] Vollständig auf Deutsch
- [x] Alle mathematischen Strukturen konsistent mit EN-Version

**D. Vollständigkeit**
- [x] 10 Sections
- [x] Konkrete Beispiele vorhanden
- [x] Fazit-Section vorhanden

**Gefundene Bugs:** Keine

**Status: DRUCKREIF**

---

## Paper 42 (Kurepa-Vermutung)

### EN — `paper42_kurepa_conjecture_en.tex`

**A. Mathematische Korrektheit**
- [x] Kurepa-Vermutung als `Conjecture` deklariert
- [x] Formel korrekt: `!p = sum_{k=0}^{p-1} k!` — korrekt
- [x] Wilson-Verbindung korrekt: `!p ≡ !(p-1) - 1 (mod p)` — korrekt bewiesen
- [x] Wilson-Theorem als Theorem mit Quellenangabe [Wilson1770, Lagrange1771]
- [x] Anfangswerte korrekt: !0=0, !1=1, !2=2, !3=4, !4=10, !5=34, !6=154, !7=874
- [x] Numerische Verifikation: 664.578 Primzahlen bis 10^7, 0 Verletzungen, ~8396s — korrekt (π(10^7)=664579, minus p=2 ergibt 664578)
- [x] Integraldarstellung korrekt hergeleitet
- [x] Äquivalenz Kurepa ↔ !(p-1) ≢ 1 (mod p) korrekt bewiesen

**Hinweis (geringfügig):** In der Historiektabelle erscheint „Specialist Mathematics Project (Build 131)" mit Schranke `p ≤ 10^7`, obwohl frühere Arbeiten bereits bis `p ≤ 10^8` gingen (nach 2010). Unsere Verifikation liegt damit UNTER dem aktuellen Stand der Literatur — das Paper erklärt dies korrekt als eigenständige Verifikation, aber die Tabellen-Reihenfolge (erst 10^8, dann unser 10^7) könnte verwirrend wirken. Kein mathematischer Fehler.

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Alle Theorem-Environments definiert inkl. `conjecture`
- [x] `\tableofcontents` nach `\maketitle`
- [x] 8 `\bibitem`-Einträge (Mindest-Anforderung: 8 — genau erfüllt)

**C. Sprachliche Qualität**
- [x] Vollständig auf Englisch

**D. Vollständigkeit**
- [x] 9 Sections (Introduction, Left Factorial, Kurepa Conjecture, Wilson's Theorem, Heuristics, Equivalences, Computation, Related Problems, Conclusion)
- [x] Konkrete Beispiele mit Zahlen (p=7, p=11)
- [x] Conclusion-Section vorhanden

**Gefundene Bugs:** Keine

**Status: DRUCKREIF**

---

### DE — `paper42_kurepa_vermutung_de.tex`

**A. Mathematische Korrektheit**
- [x] Kurepa-Vermutung als `Vermutung` deklariert (via `\newtheorem{vermutung}`)
- [x] Formeln identisch mit EN
- [x] Wilson-Satz als Satz (theorem) korrekt
- [x] Alle Werte korrekt

**Hinweis (geringfügig, kein Fehler):** Das DE-Paper verwendet `\newtheorem{vermutung}[theorem]{Vermutung}` und `\newtheorem{beispiel}[theorem]{Beispiel}` und `\newtheorem{bemerkung}[theorem]{Bemerkung}` und `\newtheorem{korollar}[theorem]{Korollar}` — abweichend von anderen DE-Papers (40, 41, 43), die `conjecture/corollary/example/remark` als Environment-Namen mit deutschen Labels verwenden. Inkonsistenz, aber kein LaTeX-Fehler.

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Theorem-Environments vorhanden (eigene Namensgebung, s.o.)
- [x] `\tableofcontents` nach `\maketitle`
- [x] 8 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Deutsch

**D. Vollständigkeit**
- [x] 9 Sections
- [x] Konkrete Beispiele (p=7, p=11, p=13)
- [x] Schluss-Section vorhanden

**Gefundene Bugs:** Keine

**Status: DRUCKREIF**

---

## Paper 43 (Lehmer τ-Vermutung)

### EN — `paper43_lehmer_tau_en.tex`

**A. Mathematische Korrektheit**
- [x] Lehmer-Vermutung korrekt als `Conjecture` deklariert (conj:lehmer)
- [x] τ-Werte korrekt: τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830, τ(6)=-6048, τ(7)=-16744, τ(8)=84480 — alle numerisch verifiziert
- [x] τ(100) = 37.534.859.200 — korrekt (verifiziert: τ(4)·τ(25) = -1472·(-25499225) = 37534859200)
- [x] τ(1000) = -30.328.412.970.240.000 — korrekt (verifiziert: τ(8)·τ(125) = 84480·(-359001100500))
- [x] Deligne-Theorem als Theorem mit Beweisskizze und Quellenangabe
- [x] Multiplikativität als Theorem (Mordell 1917) mit Quellenangabe
- [x] Hecke-Rekurrenz τ(p^{k+1}) = τ(p)τ(p^k) - p^{11}τ(p^{k-1}) korrekt
- [x] Eulerprodukt korrekt: L(s,Δ) = Π(1 - τ(p)p^{-s} + p^{11-2s})^{-1}
- [x] Kongruenz τ(n) ≡ σ_{11}(n) (mod 691) korrekt als Theorem (Ramanujan 1916)
- [x] Verifikation n ≤ 50.000, τ(n) ≠ 0, 4330s — konsistent mit Angaben
- [x] Sato-Tate (Taylor et al.) korrekt als bewiesener Satz erwähnt

**Gefundene Bugs (KRITISCH — behoben):**
1. **Fehlendes `\newtheorem{conjecture}`**: Das Environment `\begin{conjecture}` (Zeile 253) wurde verwendet, aber nicht definiert — würde LaTeX-Fehler verursachen. **Behoben: `\newtheorem{conjecture}[theorem]{Conjecture}` hinzugefügt.**
2. **`\maketitle` vor `\begin{abstract}`**: In amsart sollte abstract VOR maketitle stehen. **Behoben: Reihenfolge korrigiert.**

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Alle Theorem-Environments definiert (nach Korrektur inkl. conjecture)
- [x] `\tableofcontents` vorhanden (nach \maketitle nach Korrektur)
- [x] 10 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Englisch

**D. Vollständigkeit**
- [x] 11 Sections
- [x] Tabelle mit τ-Werten n=1..12
- [x] Conclusion-Section vorhanden

**Status: DRUCKREIF** (nach den zwei Korrekturen)

---

### DE — `paper43_lehmer_tau_de.tex`

**A. Mathematische Korrektheit**
- [x] Lehmer-Vermutung als `Vermutung` deklariert (conjecture environment, korrekt)
- [x] τ-Werte identisch mit EN-Version — korrekt
- [x] τ(100) = 37.534.859.200 und τ(1000) = -30.328.412.970.240.000 korrekt
- [x] Alle Theoreme korrekt übertragen
- [x] Konsequenzen einer hypothetischen Nullstelle korrekt: τ(p₀²) = -p₀^{11}, τ(p₀³) = 0, τ(p₀⁴) = p₀^{22} — korrekt

**Gefundener Bug (KRITISCH — behoben):**
- **`\maketitle` vor `\begin{abstract}`**: Gleiche Inkonsistenz wie EN-Version. **Behoben: Reihenfolge korrigiert.**

**B. LaTeX-Format**
- [x] `\documentclass[12pt,a4paper]{amsart}`
- [x] `\author{Michael Fuhrmann}`
- [x] `\address{Specialist Mathematics Project, Build 133}`
- [x] Theorem-Environments auf Deutsch (Satz, Lemma, Korollar, Proposition, Definition, Beispiel, Bemerkung, Vermutung)
- [x] `\tableofcontents` vorhanden (nach Korrektur der Reihenfolge)
- [x] 10 `\bibitem`-Einträge

**C. Sprachliche Qualität**
- [x] Vollständig auf Deutsch
- [x] Konsistent mit EN-Version

**D. Vollständigkeit**
- [x] 11 Sections
- [x] Tabelle mit τ-Werten auf Deutsch
- [x] Schlussfolgerung-Section vorhanden

**Status: DRUCKREIF** (nach Korrektur der abstract/maketitle-Reihenfolge)

---

## Zusammenfassung aller Korrekturen

### Durchgeführte Korrekturen

| Paper | Typ | Beschreibung |
|-------|-----|--------------|
| paper40 EN | KRITISCH | `\address`: Build 131 → Build 133 |
| paper40 DE | KRITISCH | `\address`: Build 131 → Build 133 |
| paper43 EN | KRITISCH | Fehlendes `\newtheorem{conjecture}` hinzugefügt |
| paper43 EN | MITTEL | abstract/maketitle-Reihenfolge korrigiert (amsart-Konvention) |
| paper43 DE | MITTEL | abstract/maketitle-Reihenfolge korrigiert (amsart-Konvention) |

### Statistik

- **Kritische Bugs**: 3 (Build-Nummern ×2, fehlendes conjecture-Environment)
- **Mittlere Bugs**: 2 (abstract/maketitle-Reihenfolge ×2)
- **Geringfügige Hinweise**: 3 (Kurepa-Tabellen-Reihenfolge, inkonsistente DE-Environment-Namen in Paper 42, babel-Paket nur in Paper 41 DE)
- **Mathematische Fehler**: 0

### Mathematische Verifikation

Alle geprüften numerischen Werte sind korrekt:
- Giuga-Identität für n=1722: φ(1722) = 480 ✓
- τ(100) = 37.534.859.200 ✓
- τ(1000) = -30.328.412.970.240.000 ✓
- Kurepa: π(10^7) - 1 = 664.578 ungerade Primzahlen ✓
- Alle τ(1)..τ(12) korrekt ✓

### Gesamtstatus: **DRUCKREIF** (nach durchgeführten Korrekturen)
