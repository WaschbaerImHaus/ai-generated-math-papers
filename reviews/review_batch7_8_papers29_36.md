# Vollständiger mathematischer Audit — Batch 7 & 8 (Papers 29–36)

**Datum:** 2026-03-12
**Build:** 15
**Gutachter:** Claude Opus (automatisierter mathematischer Peer-Review)
**Fokus:** Collatz-Vermutung (Papers 29–32), Riemannsche Hypothese (33), ABC-Vermutung (34), BSD-Vermutung (35), Navier-Stokes (36)

---

## Zusammenfassung

| Paper | Titel | Sprache | Urteil | Offene Fehler |
|-------|-------|---------|--------|----------------|
| 29 | Collatz-Grundlagen | EN | ⚠️ ÜBERARBEITUNG | BUG-B7-P29-CIRCULAR-PROP4, BUG-B7-P29-NOTATION |
| 29 | Collatz-Grundlagen | DE | ⚠️ ÜBERARBEITUNG | Gleiche Fehler wie EN |
| 30 | Collatz p-adisch | EN | ⚠️ ÜBERARBEITUNG | BUG-B7-P30-CONCLUSION |
| 30 | Collatz p-adisch | DE | ⚠️ ÜBERARBEITUNG | BUG-B7-P30-CONCLUSION |
| 31 | Tao probabilistisch | EN | ⚠️ ÜBERARBEITUNG | BUG-B7-P31-FURSTENBERG-REF, BUG-B7-P31-BIRKHOFF |
| 31 | Tao probabilistisch | DE | ⚠️ ÜBERARBEITUNG | Gleiche Fehler wie EN |
| 32 | Collatz ergodisch | EN | ⚠️ ÜBERARBEITUNG | BUG-B7-P32-CONCLUSION |
| 32 | Collatz ergodisch | DE | ⚠️ ÜBERARBEITUNG | BUG-B7-P32-CONCLUSION |
| 33 | Riemannsche Hypothese | EN | ⚠️ ÜBERARBEITUNG | BUG-B8-P33-EN-SELBERG |
| 33 | Riemannsche Hypothese | DE | ✅ DRUCKREIF | keine |
| 34 | ABC-Vermutung | EN | ✅ DRUCKREIF | keine |
| 34 | ABC-Vermutung | DE | ✅ DRUCKREIF | keine |
| 35 | BSD-Vermutung | EN | ✅ DRUCKREIF | keine |
| 35 | BSD-Vermutung | DE | ✅ DRUCKREIF | keine |
| 36 | Navier-Stokes | EN | ⚠️ ÜBERARBEITUNG | BUG-B8-P36-EN-ONSAGER |
| 36 | Navier-Stokes | DE | ⚠️ ÜBERARBEITUNG | BUG-B8-P36-DE-ONSAGER |

**Druckreif:** 4 von 16
**Überarbeitung nötig:** 12 von 16

---

## Strukturierte Mängelliste

---

### BUG-B7-P30-CONCLUSION — BESTÄTIGT (aus vorherigem Audit)

- **Paper-ID:** Paper 30
- **Sprache:** EN + DE
- **Schwere:** HOCH (inhaltlicher Widerspruch)
- **Zeile:** EN: 381 / DE: 390
- **Beschreibung:** Die Zusammenfassung schreibt:
  - EN: „$T$ is not measure-preserving with respect to Haar measure, but $S$ is ergodic."
  - DE: „$T$ ist nicht maßerhaltend bezüglich des Haar-Maßes, aber $S$ ist ergodisch."

  Dies stellt Ergodizität als **etablierte Tatsache** dar. Im Textkörper (Section 4, Zeile 347–363 EN, insb. BUG-B7-P30-EN/DE-004) ist Ergodizität korrekt als **Vermutung/unbewiesen** deklariert. Direkter Widerspruch zwischen Body und Conclusion.
- **Korrektur:**
  - EN: „but $S$ is ergodic" → „but $S$ is conjectured to be ergodic (see Conjecture 4.X)"
  - DE: „aber $S$ ist ergodisch" → „aber $S$ wird als ergodisch vermutet (siehe Vermutung 4.X)"

---

### BUG-B7-P32-CONCLUSION — BESTÄTIGT (aus vorherigem Audit)

- **Paper-ID:** Paper 32
- **Sprache:** EN + DE
- **Schwere:** HOCH (inhaltlicher Widerspruch)
- **Zeile:** EN: 371 / DE: 377–378
- **Beschreibung:** Die Zusammenfassung listet:
  - EN: „The Syracuse map $S$ is ergodic and mixing with exponential decay of correlations."
  - DE: „Die Syrakus-Abbildung $S$ ist ergodisch und mischend mit exponentiellem Abfall von Korrelationen."

  Im Textkörper sind alle drei Eigenschaften korrekt als **bedingt/unbewiesen** markiert:
  - Theorem 3.3 (Zeile 135–145): Ergodensatz explizit als „conditional" bezeichnet
  - Remark 5.4 (Zeile 202–207): Spektrale Lücke als „key unproved claim" markiert
  - Theorem 4.2 (Zeile 271): Als Conjecture herabgestuft (BUG-B7-P32-EN/DE-004 behoben)

  Die Conclusion ignoriert alle diese Einschränkungen und präsentiert die Ergebnisse als bewiesen.
- **Korrektur:**
  - EN: „$S$ is ergodic and mixing" → „$S$ is conjectured to be ergodic and mixing (conditional on the unproved spectral gap; see Remark 5.4)"
  - DE: analog

---

### BUG-B7-P29-CIRCULAR-PROP4 — BESTÄTIGT (aus vorherigem Audit)

- **Paper-ID:** Paper 29
- **Sprache:** EN + DE
- **Schwere:** MITTEL (logischer Fehler im Beweis)
- **Zeile:** EN: 184–199 / DE: 190–206
- **Beschreibung:** Proposition 4.4 (prop:pos_density) behauptet **unbedingt**: „The set $A_M$ has natural density 1." Der Beweis beginnt jedoch mit:
  - EN (Zeile 192): „Every positive integer eventually reaches a power of 2 under the conjecture (if true). Assuming the conjecture, all integers are in $A_M$."

  Der Beweis mischt bedingtes Reasoning (unter Annahme der Collatz-Vermutung) mit einer scheinbar unbedingten Schlussfolgerung über ergodische Abschätzungen (Zeile 197: „by ergodic-theoretic estimates (see Terras1976)"). Das Statement der Proposition ist unbedingt formuliert, aber der Beweis stützt sich teils auf die Vermutung.

  **Präzisierung:** Der zweite Teil des Beweises (Rückwärtsiterations-Baum, Terras-Abschätzung) ist tatsächlich unbedingt gültig. Das Problem ist die **verwirrende Beweisstruktur**: Zeilen 192–193 suggerieren, der Beweis benötige die Vermutung, obwohl der eigentliche Beweis (Zeilen 194–198) davon unabhängig ist.
- **Korrektur:** Den bedingten Teil (Zeilen 192–193) entfernen oder als separate Motivation kennzeichnen. Der Beweis via Rückwärtsbaum/Terras ist eigenständig und sollte allein stehen.

---

### BUG-B7-P31-FURSTENBERG-REF — BESTÄTIGT (aus vorherigem Audit)

- **Paper-ID:** Paper 31
- **Sprache:** EN + DE
- **Schwere:** MITTEL (mathematisch ungenaue Referenz)
- **Zeile:** EN: 269–276 / DE: 276–283
- **Beschreibung:** Lemma 5.6 (lem:furstenberg) ist als **Lemma** deklariert, enthält aber:
  - Keine präzise mathematische Aussage
  - Keinen Beweis
  - Nur eine vage Beschreibung: „provides a structural result connecting combinatorial properties of the Collatz orbit with the ergodic structure"

  Furstenberg 1977 beweist das Korrespondenzprinzip für Szemerédi-artige Aussagen über arithmetische Progressionen — es ist **kein** Collatz-spezifisches Werkzeug. Die Anwendung auf Collatz-Dichteaussagen ist nicht offensichtlich und wird nicht hergeleitet.

  Zum Vergleich: Paper 32 (Zeile 211–214) formuliert dasselbe Prinzip korrekt als allgemeines Theorem mit präziser mathematischer Aussage.
- **Korrektur:** Entweder (a) das Lemma durch eine präzise Formulierung mit Beweis/Beweisskizze ersetzen, oder (b) auf Paper 32, Theorem 5.1 (thm:furstenberg) verweisen, oder (c) als Remark mit korrekter Beschreibung umwandeln.

---

### BUG-B7-P29-NOTATION — NEU

- **Paper-ID:** Paper 29
- **Sprache:** EN + DE
- **Schwere:** MITTEL (Konventionsbruch, Verwechslungsgefahr)
- **Zeile:** EN: 88–91 / DE: 88–91
- **Beschreibung:** Definition 2.3 verwendet:
  - $\sigma(n)$ für die **totale** Stoppzeit (smallest $k$ with $T^k(n) = 1$)
  - $\sigma_\infty(n)$ für die (einfache) **Stoppzeit** (smallest $k$ with $T^k(n) < n$)

  Die Standardkonvention (Lagarias 1985, Terras 1976, Tao 2022) ist **umgekehrt**:
  - $\sigma(n)$ = Stoppzeit (first $k$ with $T^k(n) < n$)
  - $\sigma_\infty(n)$ = totale Stoppzeit (first $k$ with $T^k(n) = 1$)

  Der Index $\infty$ soll andeuten, dass man „bis zum Ende" (bis 1) iteriert — nicht „bis erstmals absteigend". Die Papers 30–32 verwenden implizit die gleiche (vertauschte) Notation.
- **Korrektur:** $\sigma(n)$ und $\sigma_\infty(n)$ tauschen, um mit der Standardliteratur konsistent zu sein.

---

### BUG-B7-P31-BIRKHOFF — NEU

- **Paper-ID:** Paper 31
- **Sprache:** EN + DE
- **Schwere:** MITTEL (mathematische Inkorrektheit)
- **Zeile:** EN: 252–267 (Definition 5.8 / thm:ergodic_conv) / DE: 258–274
- **Beschreibung:** Theorem 5.8 wendet Birkhoffs Ergodensatz auf das System $(\mathbb{Z}_2^{\text{odd}}, \mu, S)$ an, wobei $\mu$ das Haar-Maß ist. Birkhoffs Ergodensatz erfordert jedoch, dass die Transformation $S$ **maßerhaltend** bezüglich $\mu$ ist.

  Paper 30 (Zeile 381) und Paper 32 (Zeile 369) stellen explizit fest, dass $T$ (und damit auch $S$) **nicht** maßerhaltend bezüglich des Haar-Maßes ist. Die direkte Anwendung von Birkhoff auf $(Z_2^{\text{odd}}, \mu_{\text{Haar}}, S)$ ist daher **ungültig**.

  Paper 32 löst dieses Problem korrekt, indem es ein S-invariantes Maß $\nu$ postuliert (Theorem 3.3, Zeile 136) und Birkhoff auf $(\mathbb{Z}_2^{\text{odd}}, \nu, S)$ anwendet. Paper 31 versäumt diese Unterscheidung.
- **Korrektur:** Das Maß in der Anwendung von Birkhoff von $\mu$ (Haar) auf ein $S$-invariantes Maß $\nu$ ändern (wie in Paper 32), oder explizit erklären, welches Maß gemeint ist.

---

### BUG-B8-P33-EN-SELBERG — NEU

- **Paper-ID:** Paper 33
- **Sprache:** EN
- **Schwere:** MITTEL (historischer Fehler)
- **Zeile:** EN: 244–246
- **Beschreibung:** Der Text nach dem Selberg-Theorem (Theorem zur Nullstellenanzahl auf der kritischen Geraden) schreibt:

  > „this shows a positive proportion of zeros [...] lie on the critical line, but the proportion was only shown to be positive by Levinson in 1974."

  Dies ist **historisch falsch**. Selberg (1942) zeigte bereits einen **positiven Anteil** der Nullstellen auf der kritischen Geraden ($N_0(T) \geq A \cdot T \log T$ impliziert positiven Anteil, da $N(T) \sim \frac{T}{2\pi} \log \frac{T}{2\pi}$). Levinsons Beitrag (1974) war die **quantitative Verbesserung** auf $\geq 1/3$, nicht das erstmalige Zeigen eines positiven Anteils.

  Die DE-Version hat diesen Fehler **nicht** — dort fehlt der problematische Brückensatz, und der Text springt direkt von Selberg zu Levinson.
- **Korrektur:** EN-Satz ändern zu: „Combined with (eq), this shows a positive proportion of zeros lie on the critical line. Levinson (1974) improved this to the explicit bound of at least $1/3$."

---

### BUG-B8-P36-EN-ONSAGER — NEU

- **Paper-ID:** Paper 36
- **Sprache:** EN
- **Schwere:** MITTEL (Zuschreibungsfehler)
- **Zeile:** EN: 339
- **Beschreibung:** Der Text schreibt:

  > „De Lellis and Székelyhidi (2019) proved the Onsager conjecture"

  Die **scharfe Onsager-Vermutung** (Existenz energiedissipierender $C^{1/3-\varepsilon}$-Lösungen) wurde von **Philip Isett (2018)** bewiesen, aufbauend auf der Methode von De Lellis und Székelyhidi (konvexe Integration). De Lellis und Székelyhidi lieferten die grundlegende Methodik (2009–2013) und publizierten 2019 eine **Übersichtsarbeit**, nicht den endgültigen Beweis der scharfen Vermutung.

  Die vollständige Geschichte: Nash–Kuiper (1954) → De Lellis–Székelyhidi (2009, 2013, h-Prinzip) → Buckmaster–De Lellis–Székelyhidi–Vicol (2015, Dissipation) → Isett (2018, scharfe $C^{1/3-\varepsilon}$) → Isett korrigiert von Brue–De Lellis (2020).
- **Korrektur:**
  - EN: „De Lellis and Székelyhidi (2019)" → „Isett (2018), building on the convex integration framework of De Lellis and Székelyhidi"
  - Bibitem entsprechend anpassen.

---

### BUG-B8-P36-DE-ONSAGER — NEU

- **Paper-ID:** Paper 36
- **Sprache:** DE
- **Schwere:** MITTEL (Zuschreibungsfehler)
- **Zeile:** DE: 350
- **Beschreibung:** Identischer Fehler wie BUG-B8-P36-EN-ONSAGER.
- **Korrektur:** Analog zu EN.

---

### BUG-B7-P30-BUGCOMMENTS — NEU (alle Batch-7-Papers)

- **Paper-ID:** Papers 29–32
- **Sprache:** EN + DE
- **Schwere:** KOSMETISCH (redaktionell)
- **Zeile:** Diverse (siehe Grep-Ergebnisse oben)
- **Beschreibung:** Alle vier Collatz-Papers enthalten interne BUG-Kommentare wie `% BUG-B7-P30-EN/DE-001 behoben: ...` im LaTeX-Quelltext. Diese sind interne Entwicklungsnotizen und dürfen nicht in publizierbaren Dokumenten erscheinen. Betroffen:
  - Paper 30 EN: 4 BUG-Kommentare (Zeilen 203, 223, 288, 347)
  - Paper 30 DE: 4 BUG-Kommentare (Zeilen 209, 229, 298, 354)
  - Paper 31 EN: 1 BUG-Kommentar (Zeile 171)
  - Paper 31 DE: 1 BUG-Kommentar (Zeile 175)
  - Paper 32 EN: 4 BUG-Kommentare (Zeilen 133, 134, 181, 271)
  - Paper 32 DE: 4 BUG-Kommentare (Zeilen 136, 137, 184, 274)
  - Paper 34 EN: 1 BUG-Kommentar (Zeile 301)
  - Paper 34 DE: 1 BUG-Kommentar (Zeile 308)
- **Korrektur:** Alle `% BUG-...`-Kommentare vor Publikation entfernen.

---

## Detailprüfung nach Themengebiet

### Collatz-Papers (29–32)

**Korrekt umgesetzt:**
- Collatz-Funktion und Syracuse-Funktion korrekt definiert (Paper 29)
- 2-adische Erweiterung mathematisch korrekt (Paper 30)
- Fixpunkte und Zyklen korrekt behandelt (Papers 29, 30)
- Tao-Resultat (2022, logarithmische Dichte) korrekt dargestellt und korrekt als Teilergebnis eingeordnet (Paper 31)
- Birkhoffs Ergodensatz in Paper 32 korrekt als bedingt deklariert
- Exponentieller Korrelationsabfall als unbewiesen markiert (Paper 32, Remark 5.4)
- Bedford-McMullen-Maße und symbolische Dynamik korrekt formuliert (Paper 32)
- Markov-Ketten-Approximation mathematisch korrekt (Paper 32)
- Frühere Bugs (BUG-B7-P30-001 bis -004, BUG-B7-P32-001 bis -004, BUG-B7-P31-001) korrekt behoben

**Verbleibende Probleme:**
- Conclusions beider Papers 30 und 32 widersprechen dem Textkörper (HOCH)
- Proposition 4.4 in Paper 29 mit verwirrender Beweisstruktur (MITTEL)
- Furstenberg-Referenz in Paper 31 inhaltsleer (MITTEL)
- Birkhoff-Anwendung in Paper 31 auf falsches Maß (MITTEL)
- Notation $\sigma$/$\sigma_\infty$ vertauscht (MITTEL)

### Riemannsche Hypothese (Paper 33)

**Korrekt umgesetzt:**
- RH korrekt als offenes Problem markiert
- Funktionalgleichung, Euler-Produkt, Hadamard-Produktformel korrekt
- Nullstellenfreie Regionen (de la Vallée-Poussin, Vinogradov–Korobov) korrekt
- Explizite Formeln korrekt
- Levinson 1/3, Conrey 2/5 korrekt
- Numerische Verifikation (Odlyzko, Platt) korrekt dargestellt
- GRH korrekt als Verallgemeinerung eingeführt

**Verbleibendes Problem:**
- EN: Selberg-Zuschreibung falsch (nur EN, DE ist korrekt) (MITTEL)

### ABC-Vermutung (Paper 34)

**Korrekt umgesetzt:**
- ABC-Vermutung korrekt formuliert (Radikal, Qualität)
- Explizit als OFFEN klassifiziert (Remark, Zeile 302–316)
- Mochizuki/IUT-Status korrekt: publiziert in PRIMS 2021, aber nicht akzeptiert
- Scholze-Stix-Einwand korrekt und namentlich dokumentiert
- Mason-Stothers (Polynomfall) als bewiesenes Analogon korrekt
- Konsequenzen (Fermat, Mordell, Szpiro) korrekt als bedingt markiert
- ABC-Tripel-Tabelle korrekt

**Keine Fehler gefunden.** ✅

### BSD-Vermutung (Paper 35)

**Korrekt umgesetzt:**
- Mordell-Weil-Theorem korrekt formuliert ($E(\mathbb{Q}) \cong \mathbb{Z}^r \oplus T$)
- Schwache und starke BSD korrekt unterschieden
- Coates-Wiles (1977, Rang 0), Gross-Zagier (1986, Heegner-Punkte), Kolyvagin (1989, Euler-Systeme) korrekt
- Sha-Gruppe korrekt als endlich vermutet
- Iwasawa-theoretischer Zugang korrekt skizziert
- BSD als offenes Millennium-Problem markiert
- Numerische Daten (Cremona, LMFDB) korrekt referenziert

**Keine Fehler gefunden.** ✅

### Navier-Stokes (Paper 36)

**Korrekt umgesetzt:**
- Millennium-Problem-Formulierung korrekt (globale Existenz und Regularität starker Lösungen in 3D)
- Schwache Lösungen (Leray 1934) vs. starke Lösungen korrekt unterschieden
- Sobolev-Räume korrekt verwendet
- CKN-Partialregularität (1982) korrekt dargestellt
- Buckmaster-Vicol (2019) Nicht-Eindeutigkeit korrekt
- Stochastische Navier-Stokes korrekt erwähnt
- Energieungleichung korrekt

**Verbleibendes Problem:**
- Onsager-Zuschreibung falsch: Isett (2018), nicht De Lellis-Székelyhidi (2019) (MITTEL)

---

## Statistik der Befunde

| Schwere | Anzahl | Anteil |
|---------|--------|--------|
| KRITISCH | 0 | 0% |
| HOCH | 2 (bestätigt) | 18% |
| MITTEL | 6 (2 bestätigt, 4 neu) | 55% |
| GERING | 0 | 0% |
| KOSMETISCH | 1 (neu) | 9% |

**Bestätigte alte Bugs:** 4 von 4 (alle weiterhin offen)
**Neue Bugs:** 5 (BUG-B7-P29-NOTATION, BUG-B7-P31-BIRKHOFF, BUG-B8-P33-EN-SELBERG, BUG-B8-P36-EN-ONSAGER, BUG-B8-P36-DE-ONSAGER, BUG-B7-P30-BUGCOMMENTS)

---

*Audit durchgeführt am 2026-03-12 — Build 15*
