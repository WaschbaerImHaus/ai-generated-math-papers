# Mathematisches Audit — Batch 7 & 8 (Papers 29–36)
## Collatz-Vermutung und Millennium-Probleme
### Build 11 — 2026-03-12

---

## BATCH 7 — Collatz-Vermutung (Papers 29–32)

---

### Paper 29 (EN): Collatz Conjecture — Introduction and Survey

**Hauptbehauptung:** Einführung in die Collatz-Vermutung T(n): Überblick über Iterationsverhalten, Heuristiken, Dichteresultate (Tao 2022) und offene Fragen. Die Vermutung selbst gilt explizit als UNBEWIESEN.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert, keine falschen Beweisansprüche.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [MITTEL] Proposition 4.1, Beweis (ca. Zeile 210–230): Der Beweis der Dichteaussage über Orbits enthält eine zirkuläre Begründung — „Assuming the conjecture holds, the density of…" wird als Zwischenschritt eines Proposition-Beweises verwendet. Die Aussage selbst ist korrekt formuliert als Proposition (nicht Theorem), aber der Beweis setzt die Hauptvermutung voraus, ohne dies explizit als Konditional-Satz in der Behauptung zu formulieren. Empfehlung: Proposition entweder zu „Proposition 4.1 (conditional on Collatz conjecture)" umbenennen oder bedingte Logik in den Beweis explizit einführen.
- [GERING] Remark 4.4 (heuristische Faktorberechnung): Erst wird Faktor 1 berechnet (für nicht-2er-Teilungen), dann korrigiert auf (3/4)^{1/2}. Diese Darstellung ist inkonsistent und suggeriert einen Rechenfehler, obwohl es eine schrittweise Verfeinerung sein soll. Empfehlung: Darstellung klarer als iterative Verbesserung kennzeichnen.
- [GERING] Tao 2022 korrekt formuliert als logarithmische Dichte 1 für „almost bounded" — keine Übertreibung.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Proposition 4.1 Beweis konditionalisieren)

---

### Paper 29 (DE): Collatz-Vermutung — Einführung und Überblick

**Hauptbehauptung:** Äquivalente deutsche Version von Paper 29 EN. Gleiche Themenabdeckung, gleiche Propositions und Theoreme.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [MITTEL] Satz 4.1 (äquivalent zu Proposition 4.1 EN): Identische zirkuläre Begründung wie in EN-Version. Beweis setzt Vermutung voraus ohne konditionalen Vorbehalt im Satzstatement.
- [GERING] Bemerkung 4.4: Identische inkonsistente Faktordarstellung wie EN.
- EN/DE-Parität: gut.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (identisch zu EN)

---

### Paper 30 (EN): Collatz Conjecture — p-adic Approach

**Hauptbehauptung:** Analyse der Collatz-Funktion T und der Syracuse-Funktion S im Rahmen der 2-adischen ganzen Zahlen Z₂. Fixpunkte, Zyklen, Dichte von Orbits, Ergodizitätsfrage für S auf Z₂ mit Haar-Maß.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert. S Ergodizität korrekt als Conjecture (im Haupttext) behandelt.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [HOCH] Conclusion (Zeile 381): „but S is ergodic" wird als etablierte Tatsache geschrieben, obwohl im Textkörper (Section 4) korrekt als unbewiesen/Conjecture deklariert. Dies widerspricht dem eigenen Textkörper und könnte als Beweisanspruch missgedeutet werden.
- [INFORMATION] Frühere kritische Fehler BUG-B7-P30-001 bis -004 (T vs. S Fixpunkt-Verwechslung, fehlender Beweis c_{k+1} ∈ Z₂, fehlende Beweisskizze für Perioden aller Ordnungen, Ergodizität als Theorem) wurden in früheren Review-Runden behoben.
- [GERING] Fixpunkt -1/2 ∉ Z₂ korrekt behandelt. Unterschied T vs. S bei x=-1 (Remark 3.3) korrekt.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Conclusion Zeile 381 korrigieren)

---

### Paper 30 (DE): Collatz-Vermutung — p-adischer Ansatz

**Hauptbehauptung:** Äquivalente deutsche Version von Paper 30 EN. Gleiche mathematische Inhalte.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert im Textkörper.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [HOCH] Zusammenfassung (Zeile 390): „aber S ist ergodisch" als Faktum geschrieben, obwohl im Textkörper korrekt als unbewiesene Vermutung behandelt. Identisches Problem wie EN Zeile 381.
- EN/DE-Parität: gut.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Zusammenfassung Zeile 390 korrigieren)

---

### Paper 31 (EN): Tao's Probabilistic Approach to the Collatz Conjecture

**Hauptbehauptung:** Darstellung von Taos 2022-Ergebnis: Die Menge der Startwerte, deren Collatz-Orbit „almost bounded" ist, hat logarithmische Dichte 1 in den positiven ganzen Zahlen. Explizite und korrekte Abgrenzung von Taos Resultat gegenüber einem vollständigen Beweis der Collatz-Vermutung.

**Status der Hauptvermutung:** OFFEN — korrekt und präzise deklariert. Section 5 expliziert klar die verbleibende Lücke zwischen Taos Ergebnis und vollem Beweis.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [MITTEL] Lemma 5.3 (Furstenberg-Referenz): Zitiert Furstenberg 1977 für ein „Korrespondenzprinzip" im Collatz-Kontext. Furstenberg 1977 ist jedoch der Ergodic-Theorem-Beweis von Szemerédi's Theorem über arithmetische Progressionen — nicht ein Collatz-Korrespondenzprinzip. Die Anwendung ist mathematisch ungenau. Korrekte Referenz für das verwendete Prinzip wäre entweder Tao 2022 direkt oder eine eigene Herleitung.
- [INFORMATION] BUG-B7-P31-EN-001 (informeller Beweis nicht markiert) wurde in früherer Review-Runde behoben.
- Tao-Resultat selbst korrekt: „logarithmische Dichte 1 für almost bounded orbits" — nicht als vollständiger Beweis dargestellt.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Lemma 5.3 Referenz korrigieren)

---

### Paper 31 (DE): Taos probabilistischer Ansatz zur Collatz-Vermutung

**Hauptbehauptung:** Äquivalente deutsche Version. Taos Ergebnis korrekt dargestellt, Lücke zum vollen Beweis klar.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [MITTEL] Lemma 5.3 (Furstenberg-Referenz): Identisches Zitierproblem wie EN. Furstenberg 1977 falsch kontextualisiert.
- EN/DE-Parität: ausgezeichnet.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Lemma 5.3 Referenz korrigieren, identisch zu EN)

---

### Paper 32 (EN): Collatz Conjecture — Ergodic Theory Approach

**Hauptbehauptung:** Analyse der Collatz-Syracuse-Dynamik mit Methoden der Ergodentheorie: invariante Maße auf Z₂ (Haar-Maß), Birkhoff-Ergodensatz-Anwendung, Transfer-Operator L_S, Mischungseigenschaft, Hausdorff-Dimension von Orbits. Ergodizität und Spektrallücke werden korrekt (im Textkörper) als bedingte/unbewiesen deklariert.

**Status der Hauptvermutung:** OFFEN — im Textkörper korrekt behandelt.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [HOCH] Conclusion (Zeilen 371–373): „The Syracuse map S is ergodic and mixing with exponential decay of correlations" als bewiesene Tatsache aufgelistet, obwohl alles im Textkörper korrekt als bedingt/unbewiesen markiert. Direkter Widerspruch zum eigenen Textkörper.
- [INFORMATION] Frühere Fehler BUG-B7-P32-001 bis -004 (Eindeutigkeit von ν, zirkuläre Begründung, Spektralradius nicht rigoros, Hausdorff-Dimension als Theorem) in früherer Runde behoben.
- Spektrallücke des Transfer-Operators L_S im Textkörper korrekt als unbewiesene Behauptung markiert.
- Krylov-Bogolyubov-Theorem (Existenz invarianter Maße) korrekt angewendet.
- Bedford-McMullen-Maße und symbolische Dynamik korrekt beschrieben.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Conclusion Zeilen 371–373 korrigieren)

---

### Paper 32 (DE): Collatz-Vermutung — Ergodischer Ansatz

**Hauptbehauptung:** Äquivalente deutsche Version. Identische mathematische Inhalte.

**Status der Hauptvermutung:** OFFEN — im Textkörper korrekt.

**Mathematische Korrektheit:** ⚠️

**Fehler:**
- [HOCH] Zusammenfassung (Zeilen 379–380): „Die Syracuse-Abbildung S ist ergodisch und mischend mit exponentiellem Korrelationsabfall" als Tatsache. Identisches Problem wie EN.
- EN/DE-Parität: sehr gut.

**Urteil:** ÜBERARBEITUNG NOTWENDIG (Zusammenfassung Zeilen 379–380 korrigieren)

---

## BATCH 8 — Millennium-Probleme (Papers 33–36)

---

### Paper 33 (EN): The Riemann Hypothesis — An Analytic Survey

**Hauptbehauptung:** Überblick über die Riemann-Hypothese (RH): ζ-Funktion, Nullstellen, nullstellenfreie Gebiete (de la Vallée Poussin 1899, Korobov-Vinogradov 1958), Hardy 1914 (unendlich viele Nullstellen auf kritischer Geraden), Levinson 1/3, Conrey 2/5, Platt-Trudgian 2021 (numerische Verifikation bis 3·10^{12}). RH als offenes Millennium-Problem.

**Status der Hauptvermutung:** OFFEN — korrekt und unmissverständlich deklariert. Kein Beweisanspruch.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Alle zitierten Resultate korrekt: Hardy (1914), Levinson (1/3), Conrey (2/5), Platt-Trudgian (2021, 3·10^12 Nullstellen), Korobov-Vinogradov-Gebiet korrekt.
- Explizite Formeln (Riemann-von Mangoldt) korrekt dargestellt.
- RH als Millennium-Problem des Clay-Instituts korrekt referenziert.

**Urteil:** DRUCKREIF

---

### Paper 33 (DE): Die Riemann-Hypothese — eine analytische Übersicht

**Hauptbehauptung:** Äquivalente deutsche Version. Identische mathematische Inhalte und Quellenangaben.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Ausgezeichnete EN/DE-Parität.
- Alle Fachbegriffe korrekt ins Deutsche übersetzt.

**Urteil:** DRUCKREIF

---

### Paper 34 (EN): The abc Conjecture

**Hauptbehauptung:** Einführung in die abc-Vermutung: Radikal rad(abc), Qualitäten bekannter Tripel, Mason-Stothers-Theorem (Funktionenkörperanalogon bewiesen), Konsequenzen (Fermat, Roth, Mordell). Mochizukis IUT-Beweis: veröffentlicht in PRIMS 2021, Konsens nicht erreicht, Scholze-Stix-Einwand 2018 nicht aufgelöst.

**Status der Hauptvermutung:** OFFEN — korrekt und präzise deklariert. abc-Vermutung ausdrücklich als unbewiesen. Mochizukis Anspruch korrekt als „claimed proof; published PRIMS 2021; consensus not achieved; Scholze-Stix objection unresolved" dargestellt.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Mason-Stothers korrekt als bewiesenes Funktionenkörper-Analogon.
- Scholze-Stix-Einwand korrekt als von der Mehrheit der Gemeinschaft als stichhaltig betrachtet beschrieben.
- Qualitätsbeispiele (Brocard, Beukers) korrekt.

**Urteil:** DRUCKREIF

---

### Paper 34 (DE): Die abc-Vermutung

**Hauptbehauptung:** Äquivalente deutsche Version. Mochizukis Situation und Scholze-Stix korrekt und vollständig wiedergegeben.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Ausgezeichnete EN/DE-Parität.

**Urteil:** DRUCKREIF

---

### Paper 35 (EN): The Birch and Swinnerton-Dyer Conjecture

**Hauptbehauptung:** BSD-Vermutung: Rang r der Mordell-Weil-Gruppe E(Q) ≅ Z^r ⊕ T entspricht der Ordnung der Nullstelle von L(E,s) bei s=1. Starke BSD-Formel mit allen 5 Invarianten. Teilresultate: Coates-Wiles (1977, CM, Rang 0), Kolyvagin (1989, Rang 0 allgemein), Gross-Zagier+Kolyvagin (1986–1989, Rang 1), Skinner-Urban (2014, erste Teilbarkeit), Wei Zhang (2014, p-Umkehrung).

**Status der Hauptvermutung:** OFFEN — korrekt deklariert. BSD als Millennium-Problem korrekt referenziert. Kein Beweisanspruch für den vollen Satz.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Mordell-Weil-Struktur korrekt. Sha (Shafarevich-Tate-Gruppe) endlich nur für Rang 0,1 bekannt — korrekt.
- Starke BSD-Formel mit Ω, c_v, |Sha|, |T|², R vollständig und korrekt.
- Alle Teilresultate (Kolyvagin, Gross-Zagier, Skinner-Urban, Wei Zhang) korrekt attributiert.

**Urteil:** DRUCKREIF

---

### Paper 35 (DE): Die Birch-und-Swinnerton-Dyer-Vermutung

**Hauptbehauptung:** Äquivalente deutsche Version. Identische mathematische Inhalte.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Ausgezeichnete EN/DE-Parität.

**Urteil:** DRUCKREIF

---

### Paper 36 (EN): The Navier-Stokes Existence and Smoothness Problem

**Hauptbehauptung:** Navier-Stokes-Problem: Existenz und Glattheit globaler Lösungen. Leray 1934 (globale schwache Lösungen), Ladyzhenskaya 1969 (2D vollständig gelöst), CKN 1982 (parabolische Hausdorff-Dimension ≤ 1 der singulären Menge), Tao 2016 (Blow-up für gemitteltes NS — explizit vom echten NS unterschieden), Buckmaster-Vicol 2019 (nicht-eindeutige schwache Lösungen bei niedrigerer Regularität), Onsager-Vermutung (De Lellis-Székelyhidi 2019).

**Status der Hauptvermutung:** OFFEN — korrekt deklariert. Existenz/Eindeutigkeit glatter Lösungen als Millennium-Problem explizit als offen behandelt. Kein Beweisanspruch.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Taos 2016-Ergebnis für gemitteltes NS korrekt und explizit vom echten NS-Problem abgegrenzt — dies ist ein kritischer Punkt, der hier korrekt behandelt wird.
- CKN-Theorem (Caffarelli-Kohn-Nirenberg) korrekt mit parabolischer Hausdorff-Dimension ≤ 1.
- 2D-Vollständigkeit (Ladyzhenskaya) korrekt.
- Buckmaster-Vicol korrekt als Nicht-Eindeutigkeit bei H^s, s<1/2 (nicht bei regulären Lösungen).

**Urteil:** DRUCKREIF

---

### Paper 36 (DE): Das Navier-Stokes-Existenz-und-Glattheitsproblem

**Hauptbehauptung:** Äquivalente deutsche Version. Identische mathematische Inhalte.

**Status der Hauptvermutung:** OFFEN — korrekt deklariert.

**Mathematische Korrektheit:** ✅

**Fehler:**
- Keine mathematischen Fehler gefunden.
- Ausgezeichnete EN/DE-Parität.

**Urteil:** DRUCKREIF

---

## ZUSAMMENFASSUNG

### Gesamtübersicht

| Paper | Titel (Kurzform) | EN | DE | Neue Bugs |
|---|---|---|---|---|
| 29 | Collatz Einführung | ⚠️ ÜBERARB. | ⚠️ ÜBERARB. | BUG-B7-P29-CIRCULAR-PROP4 |
| 30 | Collatz p-adisch | ⚠️ ÜBERARB. | ⚠️ ÜBERARB. | BUG-B7-P30-CONCLUSION |
| 31 | Collatz Tao | ⚠️ ÜBERARB. | ⚠️ ÜBERARB. | BUG-B7-P31-FURSTENBERG-REF |
| 32 | Collatz Ergodisch | ⚠️ ÜBERARB. | ⚠️ ÜBERARB. | BUG-B7-P32-CONCLUSION |
| 33 | Riemann-Hypothese | ✅ DRUCKREIF | ✅ DRUCKREIF | keine |
| 34 | abc-Vermutung | ✅ DRUCKREIF | ✅ DRUCKREIF | keine |
| 35 | BSD-Vermutung | ✅ DRUCKREIF | ✅ DRUCKREIF | keine |
| 36 | Navier-Stokes | ✅ DRUCKREIF | ✅ DRUCKREIF | keine |

### Kritischer Befund: Kein falscher Beweisanspruch

**KEIN Paper aus Batch 7 oder 8 behauptet fälschlicherweise, die behandelte Vermutung bewiesen zu haben.**

- Collatz-Vermutung (Papers 29–32): Durchgehend als offen deklariert.
- Riemann-Hypothese (Paper 33): Korrekt als Millennium-Problem offen.
- abc-Vermutung (Paper 34): Mochizukis Anspruch korrekt als kontrovers und nicht konsensfähig beschrieben.
- BSD-Vermutung (Paper 35): Korrekt als offen; Teilresultate korrekt abgegrenzt.
- Navier-Stokes (Paper 36): Korrekt als offen; Taos gemitteltes NS korrekt vom echten Problem unterschieden.

### Neue Bug-IDs

| Bug-ID | Schwere | Paper | Beschreibung |
|---|---|---|---|
| BUG-B7-P29-CIRCULAR-PROP4 | Mittel | 29 EN+DE | Prop. 4.1: zirkulärer Beweis, setzt Conjecture voraus |
| BUG-B7-P30-CONCLUSION | Hoch | 30 EN+DE | Conclusion: „S ist ergodisch" als Fakt statt Vermutung |
| BUG-B7-P31-FURSTENBERG-REF | Mittel | 31 EN+DE | Lemma 5.3: Furstenberg 1977 falsch kontextualisiert |
| BUG-B7-P32-CONCLUSION | Hoch | 32 EN+DE | Conclusion: „S ergodisch und mischend" als Fakt statt Vermutung |

---

*Gutachter: Mathematischer Peer-Review-Dienst | Build 11 | 2026-03-12*
