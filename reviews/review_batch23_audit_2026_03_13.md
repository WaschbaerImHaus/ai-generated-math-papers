# Review Batch 23 -- Audit 2026-03-13

**Gutachter**: Claude Opus 4.6
**Papers**: 88--91 (EN + DE), 8 Dateien
**Ergebnis**: DRUCKREIF nach Korrekturen

---

## Paper 88: Odd Perfect Numbers / Ungerade vollkommene Zahlen

### Status: OFFEN (korrekt als Conjecture dargestellt)

| # | Schwere | Datei | Zeile | Befund | Korrektur |
|---|---------|-------|-------|--------|-----------|
| 1 | MITTEL | paper88_..._en.tex | 119 | "completely multiplicative" -- sigma ist NICHT vollstaendig multiplikativ, sondern nur multiplikativ (sigma(ab) = sigma(a)*sigma(b) nur fuer gcd(a,b)=1) | KORRIGIERT: "multiplicative on coprime integers and satisfies on prime powers" |
| 2 | MITTEL | paper88_..._de.tex | 124 | "vollstaendig multiplikativ" -- gleicher Fehler | KORRIGIERT: "multiplikativ auf teilerfremden ganzen Zahlen" |
| 3 | OK | Beide | -- | Ochem-Rao 2012: n > 10^1500 | Korrekt |
| 4 | OK | Beide | -- | Euler-Form n = p^a * m^2, p = a = 1 (mod 4) | Korrekt |
| 5 | OK | Beide | -- | omega(n) >= 9 (Chein/Hagis), Omega(n) >= 101 (Ochem-Rao 2014) | Korrekt |
| 6 | OK | Beide | -- | Nielsen 2015: omega(n) >= 10 unter Bedingung p^a > n^{1/3} | Korrekt als Teilresultat dargestellt |
| 7 | OK | Beide | -- | Conjecture "No OPN" korrekt als Conjecture markiert | Korrekt |
| 8 | OK | Beide | -- | Touchard 1953, Steuerwald 1937, Iannucci 1999 korrekt zitiert | Korrekt |

---

## Paper 89: Mersenne Primes / Mersenne-Primzahlen

### Status: OFFEN (korrekt als Conjecture fuer Unendlichkeit)

| # | Schwere | Datei | Zeile | Befund | Korrektur |
|---|---------|-------|-------|--------|-----------|
| 1 | OK | Beide | -- | 52. Mersenne-Primzahl M_{136,279,841}, Luke Durant, Oktober 2024 | Korrekt und aktuell |
| 2 | OK | Beide | -- | 41,024,320 Dezimalstellen | Korrekt |
| 3 | OK | Beide | -- | Vollstaendige Exponentenliste (52 Eintraege in EN, Auswahl in DE) | Korrekt |
| 4 | MITTEL | paper89_..._en.tex | 303-306 | EFF $150k Preis: Text suggerierte faelschlich Zusammenhang mit M_{82,589,933}. Tatsaechlich: Preis fuer 100-Mio-stellige Primzahl ist Stand 2024 noch offen. | KORRIGIERT: Klargestellt dass $150k-Preis noch offen, groesste Mersenne hat nur 41 Mio Stellen |
| 5 | MITTEL | paper89_..._de.tex | 275-276 | "zwei EFF-Preise qualifiziert" -- vage und ungenau | KORRIGIERT: Praezisiert mit $100k gewonnen, $150k noch offen |
| 6 | OK | Beide | -- | Euclid-Euler Theorem korrekt als Theorem (BEWIESEN) | Korrekt |
| 7 | OK | Beide | -- | Lucas-Lehmer Test korrekt dargestellt | Korrekt |
| 8 | OK | Beide | -- | Lenstra-Wagstaff Conjecture korrekt als Conjecture | Korrekt |

---

## Paper 90: Bunyakovsky Conjecture / Bunyakovsky-Vermutung

### Status: OFFEN (korrekt als Conjecture dargestellt)

| # | Schwere | Datei | Zeile | Befund | Korrektur |
|---|---------|-------|-------|--------|-----------|
| 1 | OK | Beide | -- | Bunyakovsky 1857 als Conjecture markiert | Korrekt |
| 2 | OK | Beide | -- | Dirichlet 1837 als Theorem (BEWIESEN) | Korrekt |
| 3 | OK | Beide | -- | Bateman-Horn 1962 als Conjecture | Korrekt |
| 4 | OK | Beide | -- | Schinzel Hypothesis H als Conjecture | Korrekt |
| 5 | OK | Beide | -- | Iwaniec 1978 (P_2-Ergebnis) als Theorem (BEWIESEN) | Korrekt |
| 6 | OK | Beide | -- | n^2 + 1 prime infinitely often als Conjecture | Korrekt |
| 7 | OK | EN | 283 | f(0) = 41, f(1) = 41 fuer n^2-n+41: mathematisch korrekt (0-0+41=41, 1-1+41=41) | Korrekt |
| 8 | OK | Beide | -- | Euler Lucky Numbers: 2,3,5,11,17,41 | Korrekt |
| 9 | OK | Beide | -- | Ono/Baker-Heegner-Stark korrekt dargestellt | Korrekt |
| 10 | OK | Beide | -- | Linnik/Xylouris L=5 (2011) | Korrekt |

**Keine Fehler gefunden. Paper ist inhaltlich einwandfrei.**

---

## Paper 91: Erdos-Gallai Theorem / Erdos-Gallai-Theorem

### Status: BEWIESEN (1960) -- korrekt als Theorem dargestellt

| # | Schwere | Datei | Zeile | Befund | Korrektur |
|---|---------|-------|-------|--------|-----------|
| 1 | OK | Beide | -- | Hauptsatz (EG 1960) als Theorem/Satz markiert | KORREKT -- KRITISCH GEPRUEFT |
| 2 | OK | Beide | -- | Hakimi 1962 als Theorem/Satz | Korrekt (bewiesen) |
| 3 | OK | Beide | -- | Gale-Ryser 1957 als Theorem/Satz | Korrekt (bewiesen) |
| 4 | OK | Beide | -- | Chvatal-Hammer 1977 als Theorem/Satz | Korrekt (bewiesen) |
| 5 | OK | Beide | -- | McKay-Wormald 1991 als Theorem/Satz | Korrekt (bewiesen) |
| 6 | OK | Beide | -- | Planare/Chordale/Ausserplanare Gradfolgen als Conjecture/Vermutung | Korrekt (offen) |
| 7 | OK | Beide | -- | Zaehlen von Realisierungen als Conjecture/Vermutung | Korrekt (offen) |
| 8 | NIEDRIG | paper91_..._en.tex | 508 | Bibitem-Key "Fulkerson1960" aber Artikel ist von 1965 | KORRIGIERT zu "Fulkerson1965" |
| 9 | OK | Beide | -- | EG-Ungleichungen mathematisch korrekt formuliert | Korrekt |
| 10 | OK | Beide | -- | Handshaking Lemma korrekt | Korrekt |
| 11 | OK | Beide | -- | Multigraph-Charakterisierung korrekt | Korrekt |

---

## Zusammenfassung

| Paper | Kritisch | Mittel | Niedrig | Status |
|-------|----------|--------|---------|--------|
| 88 EN | 0 | 1 | 0 | KORRIGIERT -> DRUCKREIF |
| 88 DE | 0 | 1 | 0 | KORRIGIERT -> DRUCKREIF |
| 89 EN | 0 | 1 | 0 | KORRIGIERT -> DRUCKREIF |
| 89 DE | 0 | 1 | 0 | KORRIGIERT -> DRUCKREIF |
| 90 EN | 0 | 0 | 0 | DRUCKREIF |
| 90 DE | 0 | 0 | 0 | DRUCKREIF |
| 91 EN | 0 | 0 | 1 | KORRIGIERT -> DRUCKREIF |
| 91 DE | 0 | 0 | 0 | DRUCKREIF |
| **Gesamt** | **0** | **4** | **1** | **ALLE DRUCKREIF** |

### Kritische Pruefpunkte (vom Auftraggeber angefordert):

1. **Theorem vs. Conjecture**: KORREKT in allen Papers. Paper 91 (Erdos-Gallai) verwendet Theorem/Satz. Paper 88 (OPN) und 90 (Bunyakovsky) verwenden Conjecture/Vermutung. Paper 89 verwendet Theorem fuer Euclid-Euler und Conjecture fuer Unendlichkeit.

2. **52. Mersenne-Primzahl**: KORREKT. M_{136,279,841}, Luke Durant, Oktober 2024, 41,024,320 Stellen -- in beiden Sprachversionen enthalten.

3. **Formeln und Beweise**: KORREKT. Euler-Struktursatz, EG-Ungleichungen, Lucas-Lehmer, Bateman-Horn-Formel alle mathematisch korrekt.

4. **LaTeX-Syntax**: Keine Fehler gefunden. Alle Umgebungen korrekt geschlossen.

### Korrigierte Fehler:
- 2x "completely multiplicative" -> "multiplicative" (Paper 88 EN+DE)
- 2x EFF-Preis-Information praezisiert (Paper 89 EN+DE)
- 1x Bibitem-Key korrigiert (Paper 91 EN)

**Alle 8 Papers nach reviewed/batch23/ verschoben.**
