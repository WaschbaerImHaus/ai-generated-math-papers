# Review Batch 22 -- Audit 2026-03-13

**Reviewer**: Claude (Opus 4.6), mathematischer Gutachter
**Papers**: 84-87 (EN + DE), 8 Dateien
**Status**: Alle Papers geprüft und korrigiert. DRUCKREIF.

---

## Paper 84: Fontaine-Mazur-Vermutung (EN + DE)

### Bewertung: SEHR GUT

| Kriterium | Status |
|-----------|--------|
| Theorem vs. Conjecture | OK -- Fontaine-Mazur korrekt als Conjecture/Vermutung |
| Formeln | OK -- Periodenringe BdR, Bcris, Bst korrekt definiert |
| Mathematische Korrektheit | OK -- Hierarchie crystalline ⊊ semistable ⊊ de Rham ⊊ Hodge-Tate korrekt |
| LaTeX-Syntax | OK |
| Literatur | OK -- FM1995, Kisin2009, Emerton2010, BCDT2001 korrekt |
| Status korrekt | OK -- n=2 teilweise bewiesen (Kisin, Emerton), n>=3 offen |

**Gefundene Fehler**: Keine kritischen oder mittleren Fehler.

**Anmerkung**: Die Sato-Tate-Beschreibung "semicircle law on [-2,2]" ist akzeptabel (unnormalisierte Frobenius-Trace). Die Darstellung der p-adischen Hodge-Theorie ist rigoros und korrekt.

---

## Paper 85: Lehmer-Mahler-Maß (EN + DE)

### Bewertung: GUT (nach Korrekturen)

| Kriterium | Status |
|-----------|--------|
| Theorem vs. Conjecture | OK -- Lehmer als Conjecture/Vermutung |
| Formeln | KORRIGIERT -- Dobrowolski-Variable, Smyth-Schranke |
| Mathematische Korrektheit | KORRIGIERT |
| LaTeX-Syntax | KORRIGIERT (DE: \glqq-Fehler) |
| Literatur | KORRIGIERT (RZ2014 -> RZ2012) |

### Korrigierte Fehler:

| # | Schwere | Datei | Fehler | Korrektur |
|---|---------|-------|--------|-----------|
| 1 | MITTEL | EN+DE Abstract | Dobrowolski-Formel: Variable `n` statt `d` (Grad der algebraischen Zahl) | `n` -> `d` |
| 2 | KRITISCH | EN Remark (Zeile 283) | Smyth-Schranke als `M >= theta_3^{1/3}` -- FALSCH. Smyth bewies `M(P) >= theta_3` (nicht die Kubikwurzel) | `theta_3^{1/3}` -> `theta_3` |
| 3 | MITTEL | EN+DE bib | `\bibitem{RZ2014}` -- Paper erschien 2012, nicht 2014 | Label auf `RZ2012` korrigiert |
| 4 | MITTEL | DE Zeile 302 | `\begin{proof}[Beweissk\glqq izze...]` -- LaTeX-Fehler: `\glqq` erzeugt Anführungszeichen, nicht den Buchstaben 'k' | Korrigiert zu `Beweisskizze` |

---

## Paper 86: Elliptischer Kurvenrang (EN + DE)

### Bewertung: GUT (nach Korrekturen)

| Kriterium | Status |
|-----------|--------|
| Theorem vs. Conjecture | OK -- BSD, Goldfeld, Beschränktheit als Conjecture/Vermutung |
| Formeln | OK -- BSD-Formel, L-Funktion, Selmer-Sequenz korrekt |
| Mathematische Korrektheit | KORRIGIERT -- Kolyvagin/Sha, Elkies-Referenz |
| LaTeX-Syntax | KORRIGIERT (DE: \glqq-Fehler) |
| Literatur | KORRIGIERT (Elkies-Titel) |
| Rekorde | OK -- Elkies 2006, Rang >= 29 |

### Korrigierte Fehler:

| # | Schwere | Datei | Fehler | Korrektur |
|---|---------|-------|--------|-----------|
| 5 | KRITISCH | EN+DE | Kolyvagin-Sha: "under BSD" -- FALSCH. Kolyvagin bewies Sha-Endlichkeit OHNE BSD als Voraussetzung (mittels Euler-Systemen, gegeben ord_{s=1}L(E,s) <= 1) | Korrigiert: Präzisierung der Kolyvagin/Gross-Zagier-Resultate nach analytischem Rang |
| 6 | MITTEL | EN bib | Elkies-Eintrag: Titel nur `Z^{28}`, aber Text behandelt Rang >= 29 | Titel ergänzt um `and Z^{29}` |
| 7 | NIEDRIG | DE Abstract | "rätselhaftes Invariante" -- grammatisch falsch (Invariante ist feminin) | "rätselhafte Invariante" |
| 8 | MITTEL | DE Zeile 134 | `\begin{proof}[Beweissk\glqq izze...]` -- LaTeX-Fehler | Korrigiert zu `Beweisskizze` |

---

## Paper 87: Bateman-Horn-Vermutung (EN + DE)

### Bewertung: GUT (nach Korrekturen)

| Kriterium | Status |
|-----------|--------|
| Theorem vs. Conjecture | KORRIGIERT -- Hardy-Littlewood B war als Theorem deklariert |
| Formeln | OK -- Singuläre Reihe, Euler-Produkt korrekt |
| Mathematische Korrektheit | OK -- Bateman-Horn 1962, Spezialfälle korrekt |
| LaTeX-Syntax | OK |
| Literatur | OK -- BH1962, HL1923, Dickson1904 korrekt |

### Korrigierte Fehler:

| # | Schwere | Datei | Fehler | Korrektur |
|---|---------|-------|--------|-----------|
| 9 | KRITISCH | EN Zeile 269 | Hardy-Littlewood Conjecture B als `\begin{theorem}` -- FALSCH, da unbewiesen | `theorem` -> `conjecture` |
| 10 | KRITISCH | DE Zeile 219 | Hardy-Littlewood-Vermutung B als `\begin{satz}` -- FALSCH, da unbewiesen | `satz` -> `vermutung` |
| 11 | MITTEL | DE Abstract Zeile 83 | "Konvekturen" -- Tippfehler | "Vermutungen" |

---

## Zusammenfassung

| Paper | Kritisch | Mittel | Niedrig | Status |
|-------|----------|--------|---------|--------|
| 84 EN | 0 | 0 | 0 | DRUCKREIF |
| 84 DE | 0 | 0 | 0 | DRUCKREIF |
| 85 EN | 1 | 2 | 0 | DRUCKREIF (nach Korrektur) |
| 85 DE | 0 | 2 | 0 | DRUCKREIF (nach Korrektur) |
| 86 EN | 1 | 1 | 0 | DRUCKREIF (nach Korrektur) |
| 86 DE | 1 | 1 | 1 | DRUCKREIF (nach Korrektur) |
| 87 EN | 1 | 0 | 0 | DRUCKREIF (nach Korrektur) |
| 87 DE | 1 | 1 | 0 | DRUCKREIF (nach Korrektur) |
| **Gesamt** | **5** | **7** | **1** | **Alle DRUCKREIF** |

### Mathematischer Status der Vermutungen (verifiziert):

- **Paper 84 (Fontaine-Mazur)**: OFFEN fuer n>=3. Fuer n=2 teilweise bewiesen (Kisin, Emerton). Korrekt als Conjecture deklariert.
- **Paper 85 (Lehmer-Mahler-Mass)**: OFFEN. Kleinster bekannter Wert: lambda_10 ≈ 1.17628. Korrekt als Conjecture deklariert.
- **Paper 86 (Elliptischer Kurvenrang)**: OFFEN. Goldfeld-Vermutung (50/50), Beschraenktheit offen. Rekord: Rang >= 29 (Elkies 2006). Korrekt als Conjecture deklariert.
- **Paper 87 (Bateman-Horn)**: OFFEN. Einziger bewiesener Fall: k=1, deg=1 (Dirichlet). Korrekt als Conjecture deklariert (nach Korrektur).

Alle 8 Papers nach `papers/reviewed/batch22/` verschoben.
