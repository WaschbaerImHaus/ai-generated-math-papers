# Vollständiger mathematischer Review — Papers 1–36 (Batches 1–8)
**Datum:** 2026-03-12
**Build:** 10
**Reviewer:** Claude Code (Anthropic Sonnet 4.6)
**Methode:** 8 parallele Spezialagenten, je ein Batch, peinlich genaue Beweisverifikation

---

## EXECUTIVE SUMMARY

36 Papers in 8 Batches geprüft. **Kein Paper behauptet fälschlicherweise, ein ungelöstes Millennium-Problem zu beweisen.** Die Mehrheit der Papers ist mathematisch korrekt. Kritische Fehler mit falschem Beweiswert existieren in **Paper 30, 31, 32 (Collatz)** und einem Darstellungsproblem in **Paper 34 (abc)**.

| Batch | Thema | Papers | Kritische Fehler | Lücken | Stil |
|-------|-------|--------|-----------------|--------|------|
| 1 | Giuga/Lehmer | 1–7 | 0 | 0 | 1 |
| 2 | Wilson-Theorem | 8–12 | 0 | 0 | 0 |
| 3 | Sieblehre | 13–16 | 0 | 2 (Skizzen) | 0 |
| 4 | Kreismethode/Vinogradov | 17–20 | 0 | 4 | 4 |
| 5 | Riemann Zeta | 21–24 | 1 | 3 | 1 |
| 6 | Elliptische Kurven/BSD | 25–28 | 0 | 4 | 1 |
| 7 | Collatz | 29–32 | 3 | 6 | 1 |
| 8 | Millennium-Probleme | 33–36 | 1 | 0 | 0 |
| **Gesamt** | | **36** | **5** | **19** | **7** |

---

## BATCH 1: Giuga-Pseudoprime & Lehmer-Totient-Problem (Papers 1–7)

### Gesamturteil: ✅ ALLE DRUCKREIF

#### Paper 1 EN/DE — Kein 3-Prim-Produkt ist Giuga-Pseudoprim
- Alle Beweise **korrekt und vollständig**: Quadratfreiheit, kein 2-Prim-Fall, Fall n=2qr (Paritätsargument), Fall n=pqr ungerade (Ordnungsabschätzung).
- CRT-Rechnungen numerisch verifiziert.
- **KEINE FEHLER.**

#### Paper 2 EN/DE — Lehmers Totient-Problem (3 Primen)
- Key Identity $(2a+1)(2b+1)(2c+1)-1 = 8abc + 4ab+4ac+4bc + 2a+2b+2c$ **algebraisch verifiziert**.
- Exhaustive Fallanalyse (p=3,5,7,≥11) — alle Fälle numerisch kontrolliert. Korrekt.
- **KEINE FEHLER.**

#### Paper 3 EN/DE — Giuga-Carmichael Unverträglichkeit
- CRT-Rechnung: gcd(900,294)=6, 6∤698 → System unlösbar. **Numerisch verifiziert.** ✅
- **KEINE FEHLER.**

#### Papers 4–6 EN/DE — Quadratfreiheit und Semiprim-Ausschlüsse
- Triviale, aber präzise Beweise via Ordnungsargument. **KORREKT.**

#### Paper 7 EN — Lehmer kein Semiprim
- **KORREKT.** Fall p≥3: $(p-2)(q-2)-1 \geq 1\cdot 2-1=1>0$.

#### Paper 7 DE — ⚠️ STILFEHLER
- **BUG-011 (bekannt):** Zeile 98: Zwischenschritt $(p-2)(q-2)-1 \geq 1\cdot p-1$ ist schwach begründet.
- Endergebnis ist korrekt, aber die Begründung ist ungenauer als die EN-Version.
- **Klassifikation: STILFEHLER** — kein Mathematikfehler, kein falsches Ergebnis.

---

## BATCH 2: Wilson-Theorem (Papers 8–12)

### Gesamturteil: ✅ ALLE DRUCKREIF — KEINE FEHLER

- **Paper 8:** Wilson-Theorem (p prim ↔ (p-1)!≡-1 mod p). Paarungsmethode und Fallanalyse vollständig korrekt. Behandlung p=2, p=4, p=p^k (ungerade), p=2^k (alle Fälle) präzise.
- **Paper 9:** Wilson für Primzahlpotenzen. Lift-Methode für Primitivwurzeln korrekt. 2^k-Struktur Z/2×Z/2^(k-2) korrekt.
- **Paper 10:** Wilson-Quotient. Harmonische Summenformel W_p ≡ μ_p + S_p (mod p) — numerisch für p=5,7 verifiziert. ✅ Bernoulli-Bezug als Skizze klar deklariert.
- **Paper 11:** Allgemeiner Wilson-Satz für abelsche Gruppen. Involutions-Untergruppe (Z/2Z)^k korrekt. Gauss-Verallgemeinerung: n∈{1,2,4,p^k,2p^k} korrekt hergeleitet.
- **Paper 12:** Anwendungen. Alle 4 Theoreme (Wilson-Primalitätstest, -1 als QR, Wolstenholme, Fermat via Wilson) mathematisch korrekt.

---

## BATCH 3: Sieblehre (Papers 13–16)

### Gesamturteil: ✅ KORREKT — Beweisskizzen klar deklariert

- **Papers 13–14 (Brun, Selberg):** Vollständige, rigorose Beweise. Mertens-Produkt, Cauchy-Schwarz in Lagrange-Minimierung, V_h(z)=1/V(z) algebraisch verifiziert. **KEINE FEHLER.**
- **Papers 15–16 (Großes Sieb, Chen):** Beweisskizzen — klar als solche deklariert. Mathematisch plausibel, alle Schritte referenzieren korrekte Literatur.
- Alle dokumentierten Bugs aus früheren Reviews sind behoben (BUG-P13 bis P16-Serien).

---

## BATCH 4: Kreismethode & Vinogradov (Papers 17–20)

### Gesamturteil: ⚠️ MEHRERE LÜCKEN — kein falscher Beweis

#### Paper 17 EN/DE — Kreismethode: ✅ KORREKT
Weylsche Ungleichung, Farey-Partition, Haupt/Nebenbogen-Methode — alle korrekt.

#### Paper 18 EN/DE — Vinogradov 3-Primzahlen-Satz: ⚠️ LÜCKEN
- **BUG-B4-P18-EN-001:** Remark 4.3 — Singulärreihenschranke S(n) ≫ (log log n)^{-1} ist korrekt, aber schwach formuliert.
- **BUG-B4-P18-DE-002:** Simplex-Fallunterscheidung im Beweis des singulären Integrals unvollständig — Formel für n≤N nur näherungsweise, nicht exakt begründet.
- **HINWEIS zu BUG-B4-P18-DE-003 (zurückgezogen):** Der Agent meldete zunächst c_2(n)=(-1)^n als falsch. Nach Überprüfung ist c_2(n) = e(n/2) = (-1)^n korrekt (für ungerades n: -1). **Kein Fehler in den Papers.**

#### Paper 19 EN/DE — Waringsches Problem: ⚠️ REFERENZFEHLER
- **BUG-B4-P19-EN-001:** Querverweis "Paper 17, Theorem 6.1" — Theorem-Nummer stimmt nicht mit Paper 17 überein.
- **BUG-B4-P19-DE-001,-002,-003:** Fehlende Inhalte aus früherer Review weiterhin offen.
- Weyls Schranke G(k) ≤ (k-2)2^{k-1}+5 ist korrekt für asymptotisches k, aber für kleine k nicht optimal (G(2)=4, Weyl ergibt 5).

#### Paper 20 EN/DE — Goldbach singuläre Reihe: ⚠️ LÜCKE
- **BUG-B4-P20-EN-002 (FALSE POSITIVE — zurückgezogen):** Tabelle der S(n)-Werte ist nach vollständiger Verifikation **korrekt**. Eintrag aus MEMORY.md ist veraltet.
- **BUG-B4-P20-EN-001:** GRH-Beweisskizze (Theorem R_2 unter GRH) enthält unvollständige Argumentation — Fehlerterm O(N^{1/2}(log N)^2) ohne explizite q-Abhängigkeit.

---

## BATCH 5: Riemann-Zeta-Funktion (Papers 21–24)

### Gesamturteil: ⚠️ Papers 21–22 OK; Papers 23–24 mit Fehlern

#### Papers 21–22 EN/DE — Zeta-Funktion, nichttriviale Nullstellen: ✅ DRUCKREIF
- Alle Definitionen, Funktionalgleichung, numerische Werte (γ₁=14.134725...) korrekt.
- Riemann-Siegel-Formel, Euler-Produkt, Gram-Punkte, GUE-Paarkorrelation — alle korrekt.

#### Paper 23 EN/DE — Explizite Formel & PNT: ⚠️ FEHLER

- **BUG-B5-P23-EN/DE-001 — LÜCKE:** Konstante -ζ'(0)/ζ(0) in der ψ-Formel konzeptuell unklar. Davenport gibt log(2π)≈1.8379. Die Schreibweise ist nicht falsch, aber verwirrend — wird nicht als Grenzterm vs. Pol-Beitrag unterschieden.

- **BUG-B5-P23-EN/DE-002 — LÜCKE:** PNT-Beweis (Mertens-Argument für ζ(s)≠0 auf Re(s)=1) zu schnell abgehandelt. Der Übergang "falls ζ(1+it₀)=0, ergibt der vierfache Pol -∞" ist nicht hinreichend präzise — einfache Nullstelle erzeugt nur einfachen Pol in log ζ, kein 4-facher.

- **BUG-B5-P23-EN/DE-003 — STILFEHLER:** Schoenfeld (1976)-Schranke |π(x)-Li(x)| < (1/8π)√x log x korrekt, aber schärfere Formeln aus Schoenfeld (1976) nicht zitiert.

#### Paper 24 EN/DE — Ansätze zur RH: ⚠️ FEHLER

- **BUG-B5-P24-EN/DE-002 — HISTORISCHER FEHLER:** "Weil (1940) ... proved the Weil conjectures" ist ungenau. **Weil stellte die Vermutungen 1949 auf; Deligne bewies sie 1974.** Weil (1940) bewies nur den Fall von Kurven über endlichen Körpern.

- **BUG-B5-P24-EN/DE-001 — LÜCKE:** "Barrier 1: Explicit formulas insufficient" zu vage — die tiefe Kausalitäts-Umkehr (Formel drückt π(x) durch Nullstellen aus, nicht umgekehrt) wird nicht erklärt.

- **BUG-B5-P24-EN-003 — LÜCKE:** PT-symmetrische Quantenmechanik (Bender et al.) in EN-Version nicht erwähnt, in DE-Version korrekt behandelt. Diskrepanz zwischen EN und DE.

---

## BATCH 6: Elliptische Kurven & BSD-Vermutung (Papers 25–28)

### Gesamturteil: ✅ KORREKT — Nur geringfügige Lücken

**Alle zentralen Sätze korrekt:**
- Gruppengesetz (Weierstraß, Chord-and-Tangent) ✅
- Mazur-Torsionssatz (alle 15 Typen) ✅
- Hasse-Schranke |a_p| ≤ 2√p ✅
- Mordell-Weil-Satz ✅
- Modularitätssatz (Wiles 1995, BCDT 2001) ✅
- BSD (schwach und stark): alle Formeln korrekt (Ω_E, R_E, c_p, Sha, Tors²) ✅
- Kolyvagin (Rang 0 und 1) ✅
- Kongruente Zahlen ↔ infinite-order-Punkte: Beweis numerisch verifiziert (n=5,6,7) ✅
- Tunnell-Kriterium (A(n)=2B(n)) korrekt formuliert und bewiesen ✅

**Geringfügige Lücken:**
- Paper 25: Nagell-Lutz — keine Diskussion der Ordnung bei y₀=0
- Paper 26: Konvergenz Euler-Produkt Re(s)>3/2 nicht rigorös motiviert
- Paper 27: Sel₂-Dimensionsformel angegeben, nicht bewiesen
- Paper 28: Tunnell-Beweisskizze — quantitative Beziehung L(E_n,1) ↔ A(n)-2B(n) fehlt

**Alle als GERING eingestuft. KEINE KRITISCHEN FEHLER.**

---

## BATCH 7: Collatz-Vermutung (Papers 29–32)

### Gesamturteil: ⚠️ LÜCKEN — kein falscher Vollbeweis, aber Fehler in Teilbeweisen

**Positiv:** Kein Paper behauptet, die Collatz-Vermutung vollständig zu beweisen. Taos Resultat ("almost all") wird korrekt als Teilresultat deklariert.

#### Paper 29 EN/DE — Grundlagen: ✅ DRUCKREIF
Existenz von Vorgängern, Zyklusstrukturgleichung, Tao-Ergebnis — alle korrekt.

#### Paper 30 EN/DE — p-adische Zahlen: ⚠️ FEHLER

- **BUG-B7-P30-EN-002 — HOCH (Beweislücke):** Syracuse Representation Theorem S^k(x) = (3^k x + c_k)/2^{s_k} — Induktionsbeweis zeigt NICHT, dass c_{k+1} ∈ Z₂ liegt. Lücke im Beweis.

- **BUG-B7-P30-EN-001 — MITTEL:** Verwechslung T vs. S: −1 ist Fixpunkt von S, aber Element eines 2-Zyklus unter T. Didaktisch verwirrend.

- **BUG-B7-P30-EN-003 — MITTEL:** "T hat periodische Orbits aller Perioden in Z₂" — unbewiesen behauptet.

- **BUG-B7-P30-EN-004 — MITTEL:** Ergodizität von S behauptet, nicht bewiesen.

- **BUG-B7-P30-DE-001 — MITTEL:** Direkte Verwirrung in Fixpunktherleitung (−1/2 vs. −1).

#### Paper 31 EN/DE — Taos probabilistischer Ansatz: ⚠️ FEHLER

- **BUG-B7-P31-EN/DE-001 — HOCH:** Key Estimate |F_k(ξ)| ≤ e^{-c(δ)k} — Proof Sketch zu informell. Der Übergang von der Summe zum Produkt ist nicht rigoros. Aussage über ξ ∉ Z/log₂3 ist für irrationales log₂3 trivial (immer erfüllt) und damit verwirrend.

#### Paper 32 EN/DE — Ergodentheorie: ⚠️ MEHRERE FEHLER

- **BUG-B7-P32-EN/DE-004 — HOCH:** Theorem über Hausdorff-Dimension der Ausnahmemenge — behauptet dim < 1, aber die Implikation "logdens 0 ⟹ Hausdorff dim < 1" ist NICHT allgemein gültig ohne zusätzliche Bedingungen. **Theorem ist unbewiesen/spekulativ.**

- **BUG-B7-P32-EN/DE-001 — MITTEL:** Existenz und Eindeutigkeit des S-invarianten Maßes ν wird vorausgesetzt, nicht begründet.

- **BUG-B7-P32-EN/DE-002 — MITTEL:** Mixing-Eigenschaft von S ist zirkulär begründet (verweist auf Lemma, das danach kommt).

- **BUG-B7-P32-EN/DE-003 — MITTEL:** Spektralradius des Transfer-Operators nicht rigoros berechnet.

---

## BATCH 8: Millennium-Probleme (Papers 33–36)

### Gesamturteil: ⚠️ Paper 34 mit Darstellungsfehler; Rest korrekt

#### Paper 33 EN/DE — Riemannsche Hypothese (analytischer Ansatz): ✅ DRUCKREIF
- Zero-free regions (de la Vallée Poussin, Korobov-Vinogradov) korrekt.
- Levinson (1/3), Conrey (2/5) korrekt zitiert.
- Platt-Trudgian 2021 (3×10¹² Nullstellen) korrekt.
- Kein falscher Beweisanspruch. **KEINE FEHLER.**

#### Paper 34 EN/DE — ABC-Vermutung: ⚠️ DARSTELLUNGSFEHLER

- **KRITISCH:** Das Paper behandelt Mochizukis IUT als "disputed claimed proof" — korrekt in der Sache, aber **nicht deutlich genug**, dass die abc-Vermutung Stand 2026 als **UNBEWIESEN** gilt.
- Scholze-Stix-Einwand (2018) wird erwähnt, aber nicht als von der Mehrheit der Gemeinschaft als valide eingestuft dargestellt.
- Alle mathematischen Details (Radikal, Höhenfunktion, Mason-Stothers, abc-Tripel-Verifikation) sind korrekt.
- **Empfehlung:** Klarstellenden Absatz hinzufügen: *"Stand März 2026: Die abc-Vermutung gilt als OFFEN. Mochizukis IUT-Beweis (PRIMS 2021) ist im Mainstream nicht akzeptiert."*

#### Paper 35 EN/DE — BSD-Vermutung: ✅ DRUCKREIF
- Mordell-Weil, Hasse, Modularitätssatz, BSD (schwach+stark), Kolyvagin, Gross-Zagier — alle korrekt.
- Bekannte Ergebnisse korrekt eingeschränkt (Rang 0 und 1 bewiesene Implikationen, Rang ≥2 offen).
- **KEINE FEHLER.**

#### Paper 36 EN/DE — Navier-Stokes: ✅ DRUCKREIF
- Leray (1934) schwache Lösungen, Ladyzhenskaya 2D-Eindeutigkeit, Beale-Kato-Majda, CKN (1982), Taos averaged NS — alle korrekt.
- Wichtig: Taos Ergebnis (2016) korrekt als **averaged** NS (nicht echte NS) deklariert.
- **KEINE FEHLER.**

---

## GESAMTLISTE ALLER OFFENEN BUGS (Stand Build 10)

### KRITISCH / HOCH (Beweis fehlerhaft oder Kernaussage falsch)

| BUG-ID | Paper | Typ | Beschreibung |
|--------|-------|-----|-------------|
| BUG-B5-P24-EN/DE-002 | Paper 24 | Historischer Fehler | "Weil (1940) proved Weil conjectures" → Weil KONJIZIERTE, Deligne BEWIES (1974) |
| BUG-B7-P30-EN/DE-002 | Paper 30 | Beweislücke | Syracuse Repr. Theorem: c_{k+1} ∈ Z₂ nicht nachgewiesen |
| BUG-B7-P31-EN/DE-001 | Paper 31 | Beweislücke | Key Estimate: Proof Sketch irreführend und nicht rigoros |
| BUG-B7-P32-EN/DE-004 | Paper 32 | Falsche Aussage | Hausdorff-dim-Theorem: logdens 0 ⟹ dim < 1 nicht allgemein gültig |
| BUG-B8-P34-STATUS | Paper 34 | Darstellung | abc-Vermutung als "disputed" statt "OFFEN" — Mochizuki-Status unklar |

### MITTEL (Lücke im Beweis, aber Ergebnis plausibel)

| BUG-ID | Paper | Beschreibung |
|--------|-------|-------------|
| BUG-B4-P18-DE-002 | Paper 18 DE | Simplex-Fallunterscheidung (sing. Integral) unvollständig |
| BUG-B4-P19-EN-001 | Paper 19 EN | Querverweis "Paper 17, Theorem 6.1" — Nummer falsch |
| BUG-B4-P20-EN-001 | Paper 20 EN | GRH-Beweisskizze lückenhaft (Fehlerterm ohne q-Abhängigkeit) |
| BUG-B5-P23-EN/DE-002 | Paper 23 | PNT-Beweis: Mertens-Argument zu schnell; einfache vs. mehrfache Pol-Verwechslung |
| BUG-B5-P23-EN/DE-001 | Paper 23 | Konstante -ζ'(0)/ζ(0) in ψ-Formel konzeptuell unklar |
| BUG-B7-P30-EN/DE-001 | Paper 30 | Verwechslung T vs. S bei Fixpunkten |
| BUG-B7-P30-EN/DE-003 | Paper 30 | "Orbits aller Perioden in Z₂" unbewiesen |
| BUG-B7-P30-EN/DE-004 | Paper 30 | Ergodizität von S behauptet, nicht bewiesen |
| BUG-B7-P32-EN/DE-001 | Paper 32 | Existenz S-invariantes Maß ν nicht begründet |
| BUG-B7-P32-EN/DE-002 | Paper 32 | Mixing-Eigenschaft zirkulär begründet |
| BUG-B7-P32-EN/DE-003 | Paper 32 | Spektralradius Transfer-Operator nicht rigoros |

### GERING / STILFEHLER

| BUG-ID | Paper | Beschreibung |
|--------|-------|-------------|
| BUG-011 | Paper 7 DE | Schwacher Zwischenschritt Zeile 98 (Endergebnis korrekt) |
| BUG-B4-P18-EN-001 | Paper 18 EN | S(n) ≫ (log log n)^{-1} — korrekt aber schwach formuliert |
| BUG-B5-P23-EN/DE-003 | Paper 23 | Schoenfeld-Schranke korrekt aber nicht schärfste Version |
| BUG-B5-P24-EN/DE-001 | Paper 24 | Barrier 1 zu vage erklärt |
| BUG-B5-P24-EN-003 | Paper 24 EN | PT-Symmetrie fehlt in EN-Version |
| BUG-B8-P34-WERT | Paper 34 | Rekord-Tripel korrekt, aber Schulze-Stix Einwand zu neutral dargestellt |

### ZURÜCKGEZOGENE BUGS (False Positives aus früheren Reviews)

| BUG-ID | Status | Begründung |
|--------|--------|------------|
| BUG-B4-P18-DE-003 | ❌ FALSE POSITIVE | c_2(n) = (-1)^n für ungerades n ist KORREKT (c₂(n) = e(n/2) = e^{πin}) |
| BUG-B4-P20-EN-002 | ❌ FALSE POSITIVE | Goldbach-Tabelle S(n) nach vollständiger Verifikation KORREKT |

---

## BEWERTUNGSMATRIX (alle 36 Papers)

| Paper | Titel | EN | DE | Urteil |
|-------|-------|----|----|--------|
| 1 | Giuga 3-Prim | ✅ | ✅ | DRUCKREIF |
| 2 | Lehmer 3-Prim | ✅ | ✅ | DRUCKREIF |
| 3 | Giuga-Carmichael | ✅ | ✅ | DRUCKREIF |
| 4 | Giuga quadratfrei | ✅ | ✅ | DRUCKREIF |
| 5 | Giuga kein Semiprim | ✅ | ✅ | DRUCKREIF |
| 6 | Lehmer quadratfrei | ✅ | ✅ | DRUCKREIF |
| 7 | Lehmer kein Semiprim | ✅ | ⚠️ Stilfehler Z.98 | DRUCKREIF (Anm.) |
| 8 | Wilson Grundsatz | ✅ | ✅ | DRUCKREIF |
| 9 | Wilson Primzahlpotenzen | ✅ | ✅ | DRUCKREIF |
| 10 | Wilson Quotient | ✅ | ✅ | DRUCKREIF |
| 11 | Wilson abelsche Gruppen | ✅ | ✅ | DRUCKREIF |
| 12 | Wilson Anwendungen | ✅ | ✅ | DRUCKREIF |
| 13 | Bruns Satz | ✅ | ✅ | DRUCKREIF |
| 14 | Selbergs Sieb | ✅ | ✅ | DRUCKREIF |
| 15 | Große Sieb-Ungleichung | ⚠️ Skizze | ⚠️ Skizze | DRUCKREIF (Skizze) |
| 16 | Chens Satz | ⚠️ Skizze | ⚠️ Skizze | DRUCKREIF (Skizze) |
| 17 | Kreismethode | ✅ | ✅ | DRUCKREIF |
| 18 | Vinogradov 3-Prim | ⚠️ Lücken | ⚠️ Lücken | ÜBERARBEITUNG |
| 19 | Waringsches Problem | ⚠️ Lücken | ⚠️ Lücken | ÜBERARBEITUNG |
| 20 | Goldbach sing. Reihe | ⚠️ Lücken | ⚠️ Lücken | ÜBERARBEITUNG |
| 21 | Riemann-Zeta-Funktion | ✅ | ✅ | DRUCKREIF |
| 22 | Nichttriviale Nullstellen | ✅ | ✅ | DRUCKREIF |
| 23 | Explizite Formel PNT | ⚠️ Fehler | ⚠️ Fehler | ÜBERARBEITUNG |
| 24 | RH-Ansätze | ⚠️ Fehler | ⚠️ Fehler | ÜBERARBEITUNG |
| 25 | Elliptische Kurven/Q | ✅ | ✅ | DRUCKREIF |
| 26 | L-Fkt. ell. Kurven | ✅ | ✅ | DRUCKREIF (Lücke) |
| 27 | BSD-Vermutung | ✅ | ✅ | DRUCKREIF (Lücke) |
| 28 | Kongruente Zahlen/BSD | ✅ | ✅ | DRUCKREIF (Lücke) |
| 29 | Collatz Grundlagen | ✅ | ✅ | DRUCKREIF |
| 30 | Collatz p-adisch | ⚠️ FEHLER | ⚠️ FEHLER | ÜBERARBEITUNG |
| 31 | Tao probabilistisch | ⚠️ FEHLER | ⚠️ FEHLER | ÜBERARBEITUNG |
| 32 | Collatz ergodisch | ⚠️ FEHLER | ⚠️ FEHLER | ÜBERARBEITUNG |
| 33 | RH analytisch | ✅ | ✅ | DRUCKREIF |
| 34 | ABC-Vermutung | ⚠️ Status | ⚠️ Status | KLARSTELLUNG NÖTIG |
| 35 | BSD Millennium | ✅ | ✅ | DRUCKREIF |
| 36 | Navier-Stokes | ✅ | ✅ | DRUCKREIF |

---

## EMPFEHLUNGEN

### Sofortiger Handlungsbedarf:

1. **Paper 24:** "Weil (1940) proved" → korrigieren zu "Weil (1940/49) conjectured; Deligne (1974) proved"
2. **Paper 30:** Syracuse Representation Theorem — Beweis von c_{k+1} ∈ Z₂ hinzufügen oder Satz als Behauptung deklarieren
3. **Paper 31:** Key Estimate — als expliziten "Proof Sketch (informal)" deklarieren; missverständliche ξ-Formulierung klären
4. **Paper 32:** Hausdorff-Dimension-Theorem — als Vermutung deklarieren oder rigoros beweisen
5. **Paper 34:** Statusabsatz zu abc-Vermutung einfügen: als OFFEN klassifizieren, Scholze-Stix als von Mehrheit anerkannt darstellen

### Mittelfristiger Handlungsbedarf:

6. **Paper 23:** PNT-Beweis präzisieren (Mertens-Pol-Argument)
7. **Papers 18–20:** Offene Querverweise und Lücken schließen
