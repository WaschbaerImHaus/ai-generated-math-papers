# Review Batch 21 -- Audit 2026-03-13

**Reviewer**: Claude (Opus 4.6), mathematischer Gutachter
**Papers**: 80--83 (EN + DE), 8 Dateien
**Datum**: 2026-03-13

---

## Zusammenfassung

| Paper | Thema | Theorem/Conjecture | Gesamtnote |
|-------|-------|-------------------|------------|
| 80 EN | Jacobian Conjecture | Korrekt als Conjecture | A- |
| 80 DE | Jacobian-Vermutung | Korrekt als Vermutung | A- (nach Fix) |
| 81 EN | Hadwiger Conjecture | Korrekt differenziert (k<=6 bewiesen, k>=7 offen) | A- (nach Fix) |
| 81 DE | Hadwiger-Vermutung | Korrekt differenziert | A- (nach Fix) |
| 82 EN | Andrews-Curtis Conjecture | Korrekt als Conjecture | A (nach Fix) |
| 82 DE | Andrews-Curtis-Vermutung | Korrekt als Vermutung | A (nach Fix) |
| 83 EN | Beal Conjecture | Korrekt als Conjecture | A |
| 83 DE | Beal-Vermutung | Korrekt als Vermutung | A (nach Fix) |

---

## Gefundene und behobene Fehler

### KRITISCH (Theorem vs. Conjecture / Mathematische Korrektheit)

#### BUG K1: Paper 80 DE -- Poisson-Satz statt Vermutung + Duplikat
- **Datei**: `paper80_jacobian_vermutung_de.tex`
- **Problem**: In Abschnitt "Weitere algebraische Aequivalenzen" waren Poisson-Vermutung und stabile Zahmheit als `\begin{satz}` deklariert. Diese sind aber OFFEN (aequivalent zur Jacobian-Vermutung, die selbst offen ist). Ausserdem war die Poisson-Vermutung doppelt vorhanden (als Satz in Sec. 8 und als Vermutung in Sec. 9).
- **Fix**: `\begin{satz}` -> `\begin{vermutung}` fuer beide, Duplikat in Sec. 9 entfernt.
- **Schwere**: KRITISCH

#### BUG K2: Paper 81 DE -- Mader-Schranke fuer k=5 falsch
- **Datei**: `paper81_hadwiger_vermutung_de.tex`
- **Problem**: Bemerkung behauptete "Mader liefert 5-Faerbbarkeit ohne VFS" fuer k=5. Korrekt: Mader liefert $(2^{k-2}+1)$-Faerbbarkeit, also $2^3+1=9$ fuer k=5, nicht 5.
- **Fix**: Korrigiert zu "9-Faerbbarkeit (da $2^{5-2}+1=9$)" und fuer k=6 "17-Faerbbarkeit (da $2^{6-2}+1=17$)".
- **Schwere**: KRITISCH (mathematischer Fehler)

#### BUG K3: Paper 83 DE -- PSS-Loesungszahl falsch
- **Datei**: `paper83_beal_vermutung_de.tex`
- **Problem**: Satz~\ref{satz:PSS} behauptete "genau 27 primitive Loesungen" fuer $x^2+y^3=z^7$. Korrekt sind 16 Loesungsfamilien (Poonen-Schaefer-Stoll 2007). Die Zahl 27 gehoert zu Beukers' Resultat fuer $x^2+y^3+z^5=0$.
- **Fix**: "27" -> "16 Loesungsfamilien".
- **Schwere**: KRITISCH (Verwechslung zweier verschiedener Resultate)

### MITTEL

#### BUG M1: Paper 81 EN -- Reed-Seymour als Theorem fuer alle k deklariert
- **Datei**: `paper81_hadwiger_conjecture_en.tex`
- **Problem**: Theorem 6.4 (Reed-Seymour 1998) behauptete polynomial-time (k-1)-coloring fuer K_k-minor-freie Graphen, aber die Formulierung vermischte bewiesene Faelle (k<=6) mit der offenen Vermutung (k>=7).
- **Fix**: Umgewandelt zu `\begin{remark}[Algorithmic coloring]` mit klarer Differenzierung.
- **Schwere**: MITTEL

#### BUG M2: Paper 81 EN -- Proposition (iii) unklar formuliert
- **Datei**: `paper81_hadwiger_conjecture_en.tex`
- **Problem**: Aussage "$K_k \preceq G$ implies $\chi(G) \geq k$ is false in general" war syntaktisch mehrdeutig und mathematisch ungenau (fehlte Gegenbeispiel).
- **Fix**: Klare Umformulierung mit Gegenbeispiel ($K_4 \preceq K_{3,3}$, aber $\chi(K_{3,3})=2$).
- **Schwere**: MITTEL

#### BUG M3: Paper 81 EN + DE -- Appel-Haken Jahreszahl inkonsistent
- **Dateien**: `paper81_hadwiger_conjecture_en.tex`, `paper81_hadwiger_vermutung_de.tex`
- **Problem**: Text sagt "Appel-Haken 1977", Bibliographie hatte Jahr 1976 (Bulletin-Ankuendigung). Die Hauptarbeit erschien 1977 im Illinois J. Math.
- **Fix**: Bibliographie auf Illinois J. Math. 1977 aktualisiert, Bulletin-Ref als Kommentar.
- **Schwere**: MITTEL

#### BUG M4: Paper 82 EN + DE -- Miller-Schupp Referenz falsch
- **Dateien**: `paper82_andrews_curtis_en.tex`, `paper82_andrews_curtis_de.tex`
- **Problem**: Bibitem `MS99` verwies auf einen falschen Artikel ("geometry of finitely presented infinite simple groups", MSRI 1992). Der korrekte Miller-Schupp-Artikel zu AC-Praesentationen ist von 1999 in Contemp. Math.
- **Fix**: Korrekte Referenz: Miller III, Schupp, "Some presentations of the trivial group", Contemp. Math. 250 (1999).
- **Schwere**: MITTEL

### GERING (keine Korrektur noetig)

- Paper 80 EN: Gesamtqualitaet sehr hoch, keine Fehler gefunden.
- Paper 82 EN/DE: AC Move (AC4) "cyclic permutation" ist technisch ein Spezialfall von (AC3) Konjugation -- korrekt so dargestellt.
- Paper 83 EN: PSS mit "16 primitive solutions" ist korrekt.
- Paper 83 EN/DE: abc-Vermutung korrekt als offen markiert, Mochizuki-Beweis korrekt als umstritten beschrieben.

---

## Theorem vs. Conjecture Pruefung (KRITISCH)

| Aussage | Paper 80 | Paper 81 | Paper 82 | Paper 83 |
|---------|----------|----------|----------|----------|
| Hauptvermutung als Conjecture | OK | OK | OK | OK |
| Bewiesene Teilresultate als Theorem | OK | OK (k<=6) | OK (n=1) | OK (FLT) |
| Offene Faelle als Conjecture/Open | OK | OK (k>=7) | OK (n>=2) | OK |
| Verwandte offene Probleme | OK* | OK | OK | OK |

*Nach Fix von BUG K1.

---

## Mathematische Korrektheit

### Paper 80 (Jacobian)
- n=1 Beweis: Korrekt und vollstaendig.
- BCW-Reduktion: Korrekt dargestellt.
- Yagzhev: Korrekt.
- Pinchuk Gegenbeispiel: Korrekt (Grad 25, Math. Z. 217, 1994).
- Dixmier-Aequivalenz: Korrekt (Tsuchimoto 2005, BK 2007).
- Jung-van der Kulk: Korrekt (amalgamiertes Produkt fuer n=2).
- Nagata-Automorphismus: Korrekt als wilder Automorphismus (Shestakov-Umirbaev 2004).

### Paper 81 (Hadwiger)
- k<=4: Korrekte Beweise.
- k=5 <-> 4CT: Korrekt via Wagners Struktursatz.
- k=6 RST 1993: Korrekt dargestellt.
- k>=7: Korrekt als offen.
- Kostochka-Thomason: Korrekt ($c \cdot k\sqrt{\log k}$ Average-Degree-Schranke).
- Norin-Song 2023: Korrekt ($\eta(G) \geq k/(\log k)^{0.999}$).
- Mader-Schranke: Nach Fix korrekt.

### Paper 82 (Andrews-Curtis)
- AC-Moves: Korrekt und vollstaendig.
- AK(n)-Familie: Korrekt (triviale Gruppe fuer alle n).
- Beweis dass AK(n) triviale Gruppe praesentiert: Korrekt.
- h/s-Kobordismus: Korrekt.
- Boone-Novikov: Korrekt.
- Tietze-Transformationen: Korrekt.

### Paper 83 (Beal)
- FLT als Spezialfall: Korrekt.
- Darmon-Granville: Korrekt (Faltings + Riemann-Hurwitz).
- abc -> Beal: Beweis korrekt.
- Beukers $x^2+y^3+z^5=0$: 27 Loesungen, korrekt.
- PSS $x^2+y^3=z^7$: 16 Loesungsfamilien, nach Fix korrekt.
- Euler-Charakteristik-Klassifikation: Korrekt.

---

## Bibliographie-Pruefung

| Paper | Eintraege | Fehlende Refs | Fehlerhafte Refs |
|-------|-----------|---------------|------------------|
| 80 EN | 12 | -- | -- |
| 80 DE | 12 | -- | -- |
| 81 EN | 12 | -- | AppHak77 (behoben) |
| 81 DE | 12 | -- | AppHak77 (behoben) |
| 82 EN | 10 | -- | MS99 (behoben) |
| 82 DE | 8 | -- | MS99 (behoben) |
| 83 EN | 10 | -- | -- |
| 83 DE | 8 | -- | -- |

---

## Gesamtbewertung

Alle 8 Papers sind nach Korrektur der 4 kritischen und 4 mittleren Bugs **DRUCKREIF**.

- **Paper 80**: Hervorragende Darstellung der Jacobian-Vermutung mit vollstaendigem n=1-Beweis, BCW/Yagzhev-Reduktion und Dixmier-Aequivalenz.
- **Paper 81**: Sehr gute Differenzierung zwischen bewiesenen (k<=6) und offenen (k>=7) Faellen. Asymptotische Resultate (Kostochka-Thomason, Norin-Song) korrekt eingeordnet.
- **Paper 82**: Solide Darstellung der AC-Vermutung mit Verbindungen zur 4-Mannigfaltigkeitstopologie.
- **Paper 83**: Umfassende Behandlung der Beal-Vermutung mit korrekter Einordnung als FLT-Verallgemeinerung und abc-Verbindung.

**Status: DRUCKREIF (nach Korrekturen)**
**Verschoben nach: papers/reviewed/batch21/**
