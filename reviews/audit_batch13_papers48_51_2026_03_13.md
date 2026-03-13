# Audit Batch 13: Papers 48-51 (Build 193)
## Datum: 2026-03-13

## Gepruefte Papers
| Paper | Thema | EN | DE |
|-------|-------|----|----|
| 48 | Additive Zahlentheorie | paper48_additive_number_theory_en.tex | paper48_additive_zahlentheorie_de.tex |
| 49 | Mertens-Funktion | paper49_mertens_function_en.tex | paper49_mertens_funktion_de.tex |
| 50 | Algorithmische Zahlentheorie | paper50_algorithmic_number_theory_en.tex | paper50_algorithmische_zahlentheorie_de.tex |
| 51 | Primzahllücken/Zwillingsprimzahlen/Andrica-Cramér | paper51_twin_primes_andrica_en.tex | paper51_zwillingsprimzahlen_andrica_de.tex |

## Gefundene und behobene Fehler

### Paper 48 DE
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 1 | MITTEL | "das kleinste ganze Zahl" (falscher Artikel) | "die kleinste ganze Zahl" (feminin) |
| 2 | MITTEL | "$h$ Faktoren" bei Sumset-Definition | "$h$ Summanden" (EN: "$h$ terms") |

### Paper 49 EN
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 3 | GERING | $12=4\cdot 3$ (verschleiert Quadratfreiheits-Argument) | $12=2^2\cdot 3$ (konsistent mit DE) |
| 4 | KRITISCH | Pólya-Vermutung: "smallest counterexample is $n = 906316571$" | Falsches Gegenbeispiel. Haselgrove (1958) nicht-konstruktiv. Kleinste: $n = 906\,150\,257$ (Tanaka, 1980) |

### Paper 49 DE
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 5 | MITTEL | "Kancelierung" (Anglizismus) im Titel und Text (3x) | "Auslöschung" (korrektes Deutsch) |
| 6 | KRITISCH | Pólya-Vermutung: falsches Gegenbeispiel $906\,316\,571$ | Korrigiert auf $906\,150\,257$ (Tanaka, 1980) |
| 7 | GERING | "kein formalen Hinweise" | "keine formalen Hinweise" |

### Paper 50 EN
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 8 | KRITISCH | "The RSA problem is polynomial-time equivalent to factoring" | FALSCH. Nur eine Richtung bewiesen (Faktorisierung => RSA). Umkehrung offen. |

### Paper 50 DE
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 9 | KRITISCH | "Lehrers Algorithmus" (2x) | "Lehmers Algorithmus" (D.H. Lehmer) |
| 10 | KRITISCH | RSA-Äquivalenz fälschlich als bewiesen dargestellt | Korrigiert: nur eine Richtung bewiesen |

### Paper 51 DE
| # | Schwere | Fehler | Korrektur |
|---|---------|--------|-----------|
| 11 | MITTEL | Fehlender Tabelleneintrag in Maximal-Gap-Tabelle | Eintrag 1129/1151/22/0,485 eingefügt (konsistent mit EN) |

## Zusammenfassung

| Schwere | Anzahl |
|---------|--------|
| KRITISCH | 5 |
| MITTEL | 3 |
| GERING | 3 |
| **Gesamt** | **11** |

### Prüfkategorien
- **Mathematische Korrektheit**: 2 kritische Fehler (Pólya-Gegenbeispiel, RSA-Äquivalenz)
- **Attributionen**: 1 kritischer Fehler (Lehmer -> "Lehrer")
- **Theorem vs. Conjecture**: Korrekt in allen Papers. Alle offenen Vermutungen als Conjecture/Vermutung deklariert.
- **Tipp-/Grammatikfehler**: 4 (Artikel, Anglizismus, fehlende Tabellenzeile, Deklination)
- **EN/DE-Konsistenz**: 3 Inkonsistenzen behoben (12=4*3, Tabellenzeile, Summanden/Faktoren)

### Gesamtbewertung
Nach Korrekturen: **DRUCKREIF**
