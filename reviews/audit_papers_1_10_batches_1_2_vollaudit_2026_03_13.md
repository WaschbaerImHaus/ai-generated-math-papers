# Vollstaendiger Audit: Papers 1-10 (Batches 1-2)
## Datum: 2026-03-13
## Auditor: Claude Opus 4.6

---

## PAPER 1 -- Giuga 3-Prim (EN + DE)
**Dateien:** `paper1_giuga_three_primes.tex`, `paper1_giuga_drei_primfaktoren_de.tex`

### BUG-SELFAUDIT-P1-EN-01: GERING
- Datei: `paper1_giuga_three_primes.tex`, Zeile 358
- Problem: Bibitem-Key `Agoh1990` suggeriert Jahr 1990, tatsaechlich steht im Eintrag Jahr 1995, Band 87.
- Korrektur: Key in `Agoh1995` aendern fuer Konsistenz.

### BUG-SELFAUDIT-P1-EN-02: MITTEL
- Datei: `paper1_giuga_three_primes.tex`, Zeilen 295-298
- Problem: Verweis auf "Bednarek (2014)" der zeige, dass Giuga-Pseudoprimen mind. 59 Primfaktoren und >10^19907 Dezimalstellen haben. Diese Referenz ist als "Preprint, unpublished" zitiert. Die Existenz und Korrektheit dieses Preprints ist nicht verifizierbar. In der gaengigen Literatur (Borwein et al. 1996) liegt die bekannte Schranke deutlich niedriger (mind. 12-13 Primfaktoren fuer Giuga-Zahlen, mind. 13660 fuer Giuga-Pseudoprimen nach Borwein et al.). Die Zahlen "59" und "10^19907" sind suspekt.
- Korrektur: Referenz entweder entfernen oder durch die belegte Borwein-et-al.-Schranke ersetzen. Falls Bednarek-Preprint existiert, DOI/arXiv-Link angeben.

### BUG-SELFAUDIT-P1-DE-01: GERING
- Datei: `paper1_giuga_drei_primfaktoren_de.tex`, Zeile 358
- Problem: Gleicher Agoh1990-Key wie EN.
- Korrektur: Analog zu EN.

### BUG-SELFAUDIT-P1-DE-02: MITTEL
- Datei: `paper1_giuga_drei_primfaktoren_de.tex`, Zeilen 296-298
- Problem: Gleiche Bednarek-Referenz wie EN.
- Korrektur: Analog zu EN.

**Mathematische Korrektheit:** KORREKT. Alle Beweisschritte nachgerechnet:
- Fall n=2qr: r=2q-1 erzwungen, dann 2(q-1)|(2q-1) Paritaetswiderspruch. CHECK.
- Fall n=pqr alle ungerade: r(r-1)|pq-1 => r(r-1) <= pq, aber r(r-1) >= (q+2)(q+1) => p >= q+4, Widerspruch zu p<q. CHECK.
- Giuga-Zahlen {30,858,1722,66198} alle verifiziert: erfuellen schwache Bedingung, NICHT starke. CHECK.

**Theorem/Conjecture:** KORREKT. Hauptsatz als `\begin{theorem}` (bewiesene Aussage). Giugas Vermutung korrekt als offene Vermutung.

**Konsistenz EN/DE:** KORREKT. Identische Struktur, gleiche Beispiele und Schranken.

**Bewertung: DRUCKREIF** (nach Bereinigung der Bednarek-Referenz)

---

## PAPER 2 -- Lehmer 3-Prim (EN + DE)
**Dateien:** `paper2_lehmer_three_primes.tex`, `paper2_lehmer_drei_primfaktoren_de.tex`

### BUG-SELFAUDIT-P2-EN-01: HOCH
- Datei: `paper2_lehmer_three_primes.tex`, Zeilen 78-79
- Problem: "Any Lehmer number has at least 15 distinct prime factors [Cohen1980]." -- FALSCHE ATTRIBUTION. Cohen & Hagis (1980) bewiesen >= 14 Primfaktoren. Die Schranke >= 15 stammt von Banks & Luca (2007/2009), nicht von Cohen (1980).
- Korrektur: Aufteilen: "Any Lehmer number has at least 14 distinct prime factors [Cohen1980]. Banks and Luca [Banks2009] improved this to 15."

### BUG-SELFAUDIT-P2-EN-02: MITTEL
- Datei: `paper2_lehmer_three_primes.tex`, Zeile 80
- Problem: "No Lehmer number exists below 10^13 (computational)." -- Ohne Quellenangabe.
- Korrektur: Quelle hinzufuegen.

### BUG-SELFAUDIT-P2-DE-01: HOCH
- Datei: `paper2_lehmer_drei_primfaktoren_de.tex`, Zeile 83-84
- Problem: "aktuelle Computerschranken liegen bei 10^30 mit mindestens 15 verschiedenen Primfaktoren [Banks2009]." -- Die 10^30-Computerschranke ist nicht durch Banks2009 belegt und INKONSISTENT mit EN (dort 10^13). Banks & Luca 2009 bewiesen eine Faktorenanzahl-Schranke, keine Computersuchschranke.
- Korrektur: Konsistente und belegbare Schranke verwenden.

### BUG-SELFAUDIT-P2-ENDE-01: MITTEL
- Problem: EN und DE verwenden UNTERSCHIEDLICHE Beweisstrategien fuer den Fall n=2qr:
  - EN (Theorem 3.1): Langer Beweis via (r-1)|(2q-1), dann b>=2 Widerspruch und b=1 gibt r=2q nicht prim.
  - DE (Satz 3.1): Kurzer Paritaetsbeweis -- phi(2qr)=(q-1)(r-1) gerade, 2qr-1 ungerade, gerade teilt nicht ungerade.
  - Beide Beweise sind mathematisch korrekt, aber die BEWEISSTRUKTUR differiert signifikant.
- Korrektur: Beweise in EN und DE vereinheitlichen.

### BUG-SELFAUDIT-P2-DE-02: MITTEL
- Datei: `paper2_lehmer_drei_primfaktoren_de.tex`, Zeilen 307-308
- Problem: Fuer p=5, q>=17: DE hat nur "vollstaendige rechnerische Ueberpruefung fuer q<=500", waehrend EN den vollstaendigen ANALYTISCHEN Beweis enthaelt (4(q-1)|k+5 unmoeglich fuer q>=17 weil k+5 <= 9 < 4(q-1) >= 64). Der EN-Beweis ist staerker, da er KEINEN Computer benoetigt.
- Korrektur: Den analytischen Beweis auch in DE aufnehmen.

**Mathematische Korrektheit:** KORREKT. Schluesselidentitaet pqr-1 = 8abc + 4(ab+bc+ca) + 2(a+b+c) nachgerechnet. Alle numerischen Fallanalysen verifiziert:
- (3,5): pq-1=14, Teiler: {1,2,7,14}, r in {2,3,8,15}, keiner prim > 5. CHECK.
- (5,11,19): phi=720, 1044/720 nicht ganzzahlig. CHECK.
- (5,13,17): phi=768, 1104/768 nicht ganzzahlig. CHECK.
- (7,13,19): phi=1296, 1728/1296 nicht ganzzahlig. CHECK.
- (7,13,31): phi=2160, 2820/2160 nicht ganzzahlig. CHECK.

**Theorem/Conjecture:** KORREKT. Lehmers Problem als offene Frage. Bewiesene Saetze als Theorem.

**Bewertung: UEBERARBEITUNGSBEDARF** (falsche Attribution Cohen1980, inkonsistente EN/DE-Beweise und Schranken)

---

## PAPER 3 -- Giuga-Carmichael CRT (EN + DE)
**Dateien:** `paper3_giuga_carmichael.tex`, `paper3_giuga_carmichael_de.tex`

### BUG-SELFAUDIT-P3-EN-01: GERING
- Datei: `paper3_giuga_carmichael.tex`, Zeile 219
- Problem: "{5, 7, 11} mit kleinstem Faktor 5 > 3: the system MAY be consistent" -- vage Formulierung ohne Berechnung.
- Korrektur: Rechnerisch pruefen und Ergebnis explizit angeben.

**Mathematische Korrektheit:** KORREKT. CRT-Rechnung vollstaendig nachgerechnet:
- n = 705 (mod 900): 705 mod 18 = 3 CHECK, 705 mod 100 = 5 CHECK
- gcd(900,294) = 6, (7-705) mod 6 = 4 != 0, also unloesbar. CHECK
- 9*39 = 351, 351 mod 50 = 1 CHECK (Inversberechnung)
- p^2(p-1) Werte: p=3: 18, p=5: 100, p=7: 294 CHECK

**Theorem/Conjecture:** KORREKT. Giuga-Carmichael-Existenz als offene Frage deklariert.

**Konsistenz EN/DE:** KORREKT. Identische Struktur und Rechnungen.

**Bewertung: DRUCKREIF**

---

## PAPER 4 -- Giuga Quadratfreiheit (EN + DE)
**Dateien:** `paper4_giuga_squarefree.tex`, `paper4_giuga_quadratfrei_de.tex`

Keine Bugs gefunden.

**Mathematische Korrektheit:** KORREKT. p^2|n => p|n/p => n/p-1 = -1 (mod p) => p|-1, Widerspruch. Elementar und lueckenlos.

**Theorem/Conjecture:** KORREKT.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF**

---

## PAPER 5 -- Giuga kein Semiprim (EN + DE)
**Dateien:** `paper5_giuga_no_semiprime.tex`, `paper5_giuga_kein_semiprim_de.tex`

### BUG-SELFAUDIT-P5-EN-01: MITTEL
- Datei: `paper5_giuga_no_semiprime.tex`, Zeilen 115-117
- Problem: Verweis auf Bednarek [Bednarek2014] mit "59 prime factors" -- gleiche problematische Referenz wie Paper 1.
- Korrektur: Analog zu Paper 1.

**Mathematische Korrektheit:** KORREKT. q|(p-1) mit p<q => p-1<q => p-1=0 => p=1, Widerspruch.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF** (nach Bereinigung der Bednarek-Referenz)

---

## PAPER 6 -- Lehmer Quadratfreiheit (EN + DE)
**Dateien:** `paper6_lehmer_squarefree.tex`, `paper6_lehmer_quadratfrei_de.tex`

### BUG-SELFAUDIT-P6-EN-01: GERING
- Datei: `paper6_lehmer_squarefree.tex`, Zeile 59
- Problem: "Cohen and Hagis proved at least 13 prime factors" -- Cohen & Hagis 1980 bewies >= 14 (je nach Interpretation der Originalarbeit). Die Zahl 13 ist moeglicherweise um 1 zu niedrig.
- Korrektur: Verifizieren und ggf. auf 14 korrigieren.

**Mathematische Korrektheit:** KORREKT. p^2|n => p|phi(n) => p|(n-1) => p|n-(n-1)=1, Widerspruch. Identisch korrekt wie Paper 4 fuer Giuga.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF**

---

## PAPER 7 -- Lehmer kein Semiprim (EN + DE)
**Dateien:** `paper7_lehmer_no_semiprime.tex`, `paper7_lehmer_kein_semiprim_de.tex`

Keine signifikanten Bugs.

**Mathematische Korrektheit:** KORREKT. pq-1 = (p-1)(q-1)+(p-1)+(q-1) nachgerechnet: pq-p-q+1+p-1+q-1 = pq-1. CHECK. Fall p=2: (q-1)|q => (q-1)|1 => q=2=p, Widerspruch. Fall p>=3: (p-2)(q-2)-1 >= 1*2-1 = 1 > 0.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF**

---

## PAPER 8 -- Wilson-Satz (EN + DE)
**Dateien:** `paper8_wilson_theorem.tex`, `paper8_wilson_satz_de.tex`

### BUG-SELFAUDIT-P8-ENDE-01: GERING
- Problem: Autor ist "Michael Fuhrmann" (Papers 1-7 verwenden "Kurt Ingwer"). Inkonsistenz im Autorennamen zwischen Batches.
- Korrektur: Autorennamen vereinheitlichen.

**Mathematische Korrektheit:** KORREKT.
- Wilson-Satz fuer alle Primzahlen p=2,3,5,7,11,13 verifiziert: (p-1)! mod p = p-1. CHECK.
- Zusammengesetzte Zahlen: (4-1)!=6 mod 4=2 (nicht -1), (6-1)!=120 mod 6=0, etc. CHECK.
- Primzahlpotenzen: Paarungsmethode und Selbstinverse korrekt.
- Allgemeiner Wilson-Satz (endliche abelsche Gruppen): |S|=1,2,>=4 Faelle korrekt.

**Theorem/Conjecture:** KORREKT. Wilson-Satz ist BEWIESEN und als Theorem deklariert. Wilson-Primzahlen-Unendlichkeit korrekt als Vermutung deklariert.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF**

---

## PAPER 9 -- Wilson fuer Primzahlpotenzen (EN + DE)
**Dateien:** `paper9_wilson_prime_powers.tex`, `paper9_wilson_primzahlpotenzen_de.tex`

Keine signifikanten Bugs.

**Mathematische Korrektheit:** KORREKT.
- Quadratwurzeln der Eins mod p^k (p ungerade): nur +-1. CHECK.
- p=2, k>=3: vier Loesungen {1, 2^{k-1}-1, 2^{k-1}+1, 2^k-1}. CHECK.
- Struktursatz (Z/p^kZ)*: zyklisch fuer ungerades p, Z/2 x Z/2^{k-2} fuer p=2, k>=3. CHECK.
- Primitivwurzel-Lifting via Hensel korrekt skizziert. CHECK.

**Konsistenz EN/DE:** KORREKT.

**Bewertung: DRUCKREIF**

---

## PAPER 10 -- Wilson-Quotient und Wilson-Primzahlen (EN + DE)
**Dateien:** `paper10_wilson_quotient.tex`, `paper10_wilson_quotient_de.tex`

### BUG-SELFAUDIT-P10-DE-01: GERING
- Datei: `paper10_wilson_quotient_de.tex`
- Problem: DE-Version fehlt die Referenz [BEW1998] (Berndt-Evans-Williams "Gauss and Jacobi Sums"), die in EN vorhanden ist. EN hat 6 Bibliographieeintraege, DE hat nur 5.
- Korrektur: [BEW1998] in DE ergaenzen.

**Mathematische Korrektheit:** KORREKT.
- Wilson-Quotienten: W_2=1, W_3=1, W_5=5, W_7=103, W_11=329891. Alle verifiziert. CHECK.
- Wilson-Primzahlen: 5^2|25 CHECK, 13^2|479001601 CHECK.
- Halbfakultat-Formel: W_p = mu_p + S_p (mod p) korrekt hergeleitet.
- p=5: mu_5=1, S_5=4, W_5=0 mod 5. CHECK.
- p=7: mu_7=-5=2 mod 7, S_7=10=3 mod 7, W_7=5 mod 7. CHECK (103 mod 7 = 5).
- Bernoulli-Zusammenhang W_p = -B_{p-1} (mod p): korrekt als Skizze dargestellt.
- Heuristik: ~ln(ln(x)) Wilson-Primzahlen bis x. CHECK (Mertens-Satz).

**Konsistenz EN/DE:** FAST KORREKT (fehlende BEW1998-Referenz in DE).

**Bewertung: DRUCKREIF**

---

## ZUSAMMENFASSUNG

| Paper | Thema | Bugs | Schwere | Bewertung |
|-------|-------|------|---------|-----------|
| 1 | Giuga 3-Prim | 2+2 | 2x MITTEL, 2x GERING | DRUCKREIF* |
| 2 | Lehmer 3-Prim | 2+2+1 | 2x HOCH, 3x MITTEL | UEBERARBEITUNGSBEDARF |
| 3 | Giuga-Carmichael CRT | 1 | 1x GERING | DRUCKREIF |
| 4 | Giuga Quadratfrei | 0 | - | DRUCKREIF |
| 5 | Giuga kein Semiprim | 1 | 1x MITTEL | DRUCKREIF* |
| 6 | Lehmer Quadratfrei | 1 | 1x GERING | DRUCKREIF |
| 7 | Lehmer kein Semiprim | 0 | - | DRUCKREIF |
| 8 | Wilson-Satz | 1 | 1x GERING | DRUCKREIF |
| 9 | Wilson Primzahlpotenzen | 0 | - | DRUCKREIF |
| 10 | Wilson-Quotient | 1 | 1x GERING | DRUCKREIF |

*DRUCKREIF nach Bereinigung der Bednarek-Referenz

### Uebergreifende Probleme

1. **Autorenwechsel:** Papers 1-7 (Batch 1) verwenden "Kurt Ingwer", Papers 8-10 (Batch 2) verwenden "Michael Fuhrmann". Sollte vereinheitlicht werden.

2. **Bednarek-Referenz:** In Papers 1 und 5 wird eine nicht verifizierbare Referenz "Bednarek 2014" (unpublished preprint) mit der Behauptung "mindestens 59 Primfaktoren, >10^19907 Dezimalstellen" zitiert. Diese Zahlen sind in der Standard-Literatur nicht bestaetigt.

3. **Cohen & Hagis Primfaktorzahl:** Inkonsistente Angaben:
   - Paper 2 EN: "at least 15" [Cohen1980] -- FALSCH (Cohen bewies 14, Banks & Luca 15)
   - Paper 6 EN: "at least 13" [Cohen1980] -- moeglicherweise 14
   - Paper 2 DE: "mindestens 15" [Banks2009] -- korrekt fuer Banks, aber Computerschranke 10^30 nicht belegt

### Gesamturteil:
- **9 von 10 Papers: DRUCKREIF** (teils mit geringfuegigen Korrekturen)
- **1 Paper (Nr. 2): UEBERARBEITUNGSBEDARF** (falsche Attribution, inkonsistente EN/DE-Beweisstrategien, inkonsistente Computerschranken)
- **Alle Beweise mathematisch korrekt** -- keine inhaltlichen Fehler in den Kernaussagen.
