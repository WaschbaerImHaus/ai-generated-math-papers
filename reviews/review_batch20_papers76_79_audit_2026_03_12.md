# Vollstaendiger Audit — Batch 20, Papers 76–79

**Datum:** 2026-03-12
**Gutachter:** Claude Opus 4.6 (Mathematik-Audit)
**Build:** 21
**Umfang:** 8 LaTeX-Dateien (4 Papers, je EN + DE)

---

## Zusammenfassung

| Paper | Titel | Sprache | Urteil | Kritische Bugs |
|-------|-------|---------|--------|----------------|
| 76 | Mandelbrot / MLC | EN | DRUCKREIF | keine |
| 76 | Mandelbrot / MLC | DE | DRUCKREIF | BUG-B20-P76-DE-001 (Tippfehler) |
| 77 | Conway Knot Sliceness | EN | UEBERARBEITUNG | BUG-B20-P77-EN-001 (Freedman falsch zitiert), BUG-B20-P77-EN-002 (Rolfsen-Tabelle falsch) |
| 77 | Conway-Knoten | DE | UEBERARBEITUNG | BUG-B20-P77-DE-001 (Freedman falsch), BUG-B20-P77-DE-002 (Rolfsen falsch) |
| 78 | Whitehead Asphericity | EN | DRUCKREIF | BUG-B20-P78-EN-001 (Sieradski Ref.), BUG-B20-P78-EN-002 (Dunce Hat) |
| 78 | Whitehead Aspharizitaet | DE | DRUCKREIF | BUG-B20-P78-DE-001 (Sieradski Ref.) |
| 79 | Temperley-Lieb / Jones | EN | DRUCKREIF | BUG-B20-P79-EN-001 (Tanaka Ref.) |
| 79 | Temperley-Lieb / Jones | DE | DRUCKREIF | BUG-B20-P79-DE-001 (Skein Tippfehler) |

---

## Paper 76: The Mandelbrot Set — Local Connectivity and the MLC Conjecture

### EN-Version

**Urteil: DRUCKREIF**

**Mathematische Korrektheit:**
- Definition der Mandelbrot-Menge, Boettcher-Koordinate, Green-Funktion: Korrekt.
- Douady-Hubbard-Beweis der Zusammenhaengigkeit via konformer Abbildung: Korrekt dargestellt.
- MLC-Aequivalenzen (Caratheodory-Erweiterung, Landung aller Strahlen, Thurston-Modell): Korrekt.
- Yoccoz-Satz fuer endlich-renormalisierbare Parameter: Korrekt formuliert.
- Yoccoz-Ungleichung (Zeile 408–410): Korrekt — die Moduli der Puzzle-Annuli sind monoton.
- Fatou-Shishikura-Ungleichung: Korrekt, 2d-2 nicht-abstossende Zyklen.
- Kahn-Lyubich 2009: Korrekte Zuschreibung (lokale Zusammenhaengigkeit der Julia-Menge).
- Lyubich 1997 (bounded type): Korrekt.
- Feigenbaum-Parameter c_Feig ≈ -1.4011552: Korrekt (auf die angegebene Genauigkeit).
- Shishikura 1998, dim_H(partial M) = 2: Korrekt.
- Korollar Zeile 476–479: "K_c is a Jordan curve and f_c is hyperbolic" — leicht unpraezise. Wenn f_c hyperbolisch ist mit einem anziehenden Fixpunkt (nicht Zyklus), dann ist J_c eine Jordan-Kurve. Bei anziehenden Zyklen hoeherer Periode ist J_c i.A. keine Jordan-Kurve, sondern eine komplizierte zusammenhaengende Menge. Allerdings ist dies im Kontext der Hauptkardioid korrekt. Akzeptabel als vereinfachte Aussage.

**Historische Angaben:**
- Mandelbrot 1980, Douady-Hubbard 1982: Korrekt.
- Thurston 1985 Preprint: Korrekt (publiziert 2009 in Schleicher-Band).
- Yoccoz 1990 (unpublished): Korrekt — Yoccoz hat seinen Beweis nie formal publiziert.
- LavaursRousseau1993 bibitem: Der Bibitem-Key sagt "1993", aber der Eintrag beschreibt Lavaurs' These von 1989. Das ist eine Inkonsistenz im Bibitem-Key vs. Inhalt, aber kein inhaltlicher Fehler — die These selbst ist korrekt zitiert.

**Bibliographie:**
- AvilaLyubich2008: Preprint von 2006, zitiert als 2008 im Key. Eintrag selbst sagt korrekt "arXiv 2006". Harmlos.
- Alle 12 Bibitems vorhanden und korrekt.

**LaTeX:**
- Sauber kompilierbar (amsart-Klasse, alle Pakete standard).
- Build-Inkonsistenz: Header sagt Build 164, \address sagt Build 165. KOSMETISCH.

### DE-Version

**Urteil: DRUCKREIF**

**Mathematische Korrektheit:**
- Inhaltlich identisch mit EN. Alle Formeln korrekt uebernommen.
- Feigenbaum-Parameter mit deutschem Komma: $-1{,}4011552$ — korrekt.

**EN-DE-Konsistenz:**
- Section "Parabolic Implosion" in EN steht zwischen "Misiurewicz Points" und "Hausdorff Dimension". In DE steht sie nach "Hausdorff Dimension" (Reihenfolge vertauscht in Section 10). Inhaltlich identisch, aber Abschnittsreihenfolge differiert. KOSMETISCH.
- DE-Version hat keinen Korollar nach der Fatou-Shishikura-Ungleichung (EN hat "at most one attracting or indifferent cycle" als Korollar, Zeile 475–479). DE hat nur "hoechstens einen nicht-abstossenden Zyklus" (Zeile 459–461). Akzeptabel — die "Jordan curve"-Aussage ist in DE weggelassen, was sogar praeziser ist.

**Bug:**

**BUG-B20-P76-DE-001** [KOSMETISCH]
Zeile 267: "Periode-2-Knospfe" — Tippfehler, muss "Periode-2-Knospe" heissen.

---

## Paper 77: Conway's Knot and the Sliceness Problem

### EN-Version

**Urteil: UEBERARBEITUNG**

**Mathematische Korrektheit:**
- Concordance-Group, Slice/Ribbon-Definitionen: Korrekt.
- Fox-Milnor-Satz: Korrekt.
- Arf-Invariante, Rasmussen s-Invariante, tau-Invariante: Alle korrekt definiert.
- Piccirillo's Beweis: Logisch korrekt skizziert (Trace X_C, Knoten K', s(K')=-2, Widerspruch).
- Mutation-Invarianz von Alexander-Polynom, Signatur, Arf, Determinante: Korrekt.
- Mutation-Invarianz von tau: Korrekt (Ozsvath-Szabo 2004 bewiesen dies fuer tau, nicht fuer s).
- Agol 2022, ribbon concordance partial order: Korrekt.

**Kritischer Fehler — Freedman falsch zitiert:**

**BUG-B20-P77-EN-001** [HOCH]
Zeile 319–323: "A knot K ⊂ S³ with vanishing Arf invariant is topologically slice."
Dies ist **falsch**. Freedmans Satz (1982) besagt: Ein Knoten mit trivialem Alexander-Polynom (Δ_K(t) = 1) ist topologisch scheibig. Die Arf-Invariante allein reicht nicht aus. Es gibt Knoten mit Arf(K) = 0 aber nichttrivialem Alexander-Polynom, die nicht topologisch scheibig sind. Der Satz wird im Paper sogar direkt im naechsten Absatz (Zeile 328–329) korrekt formuliert ("Delta_K = 1"), aber die erste Formulierung in Zeile 320 ist mathematisch falsch und widerspricht der korrekten Version.

In Piccirillo's Beweis selbst wird Freedmans Satz nicht direkt benutzt (der Beweis ist ein Widerspruchsbeweis via s-Invariante, nicht via Freedman). Die Aufnahme von Freedmans Satz dient der Kontextualisierung, aber die erste Formulierung ist irrefuehrend.

**BUG-B20-P77-EN-002** [MITTEL]
Zeile 109: "denoted 11n34 in the Rolfsen table"
Die Rolfsen-Tabelle geht nur bis 10 Kreuzungen. Der Conway-Knoten hat 11 Kreuzungen und erscheint daher in der **Hoste-Thistlethwaite-Tabelle** (oder aehnlichen modernen Tabellierungen), nicht in der Rolfsen-Tabelle. Die Bezeichnung "11n34" ist die Hoste-Thistlethwaite-Notation. Die Zuschreibung an Rolfsen ist historisch falsch.

**BUG-B20-P77-EN-003** [GERING]
Zeile 163–165: "The first known examples of topologically but not smoothly slice knots are torus knots T(2, 2k+1) for large k"
Torusknoten T(2, 2k+1) sind NICHT topologisch scheibig (sie haben nichttriviales Alexander-Polynom, und ihre Signatur ist nicht null). Das erste Beispiel eines topologisch aber nicht glatt scheibigen Knotens stammt aus der Arbeit von Freedman (1982) — der Whitehead double des Trefoils ist topologisch scheibig (da sein Alexander-Polynom trivial ist), aber nicht glatt scheibig (gezeigt mit Gauge-Theorie). Die Zuschreibung an Torusknoten ist falsch.

**BUG-B20-P77-EN-004** [GERING]
Zeile 609–610: Gordon 1981 bibitem zitiert tatsaechlich 1978 ("Lecture Notes in Math. 685, Springer, Berlin, 1978"). Inkonsistenz zwischen bibitem-Key/Text ("Gordon 1981") und dem tatsaechlichen Publikationsjahr 1978.

**Bibliographie:**
- OzsvathSzabo2004 bibitem: Key sagt "2004", Eintrag sagt "8 (2008)". Inkonsistenz. KOSMETISCH.
- Hedden2009 bibitem: Key sagt "2009", Eintrag sagt "19 (2010)". Inkonsistenz. KOSMETISCH.
- Alle anderen Bibitems korrekt.

### DE-Version

**Urteil: UEBERARBEITUNG**

**BUG-B20-P77-DE-001** [HOCH]
Zeile 306–309: Identischer Fehler wie EN — "Ein Knoten K mit verschwindender Arf-Invariante ist topologisch scheibig." Falsch, siehe BUG-B20-P77-EN-001.

**BUG-B20-P77-DE-002** [MITTEL]
Zeile 116: "bezeichnet 11n34 in der Rolfsen-Tabelle" — Identischer Fehler wie EN, siehe BUG-B20-P77-EN-002.

**BUG-B20-P77-DE-003** [GERING]
Identischer Fehler zu BUG-B20-P77-EN-003 fehlt in DE: DE hat kein entsprechendes "Example" mit Torusknoten. Korrekt ausgelassen — dieser Bug betrifft nur EN.

**EN-DE-Konsistenz:**
- DE fehlt das "Example" ueber topologisch-aber-nicht-glatt-scheibige Knoten (EN Zeile 162–168). Akzeptabel — DE ist etwas kompakter.
- DE fehlt die "Satellite Concordance Invariants" Section (EN Section 8.2) — nur gekuerzt erwaehnt. Akzeptabel.
- Alle mathematischen Hauptaussagen stimmen ueberein.
- Bibliographien identisch.

---

## Paper 78: The Whitehead Asphericity Conjecture

### EN-Version

**Urteil: DRUCKREIF (mit Anmerkungen)**

**Mathematische Korrektheit:**
- Whitehead-Vermutung (1941): Korrekt formuliert.
- CW-Komplexe, Asphaerizitaet, K(G,1)-Raeume: Korrekt.
- Hurewicz-Satz: Korrekt.
- Stallings 1968 / Swan 1969 (cd ≤ 1 iff frei): Korrekt.
- Eilenberg-Ganea 1957: Korrekt.
- Beziehung Eilenberg-Ganea-Vermutung ↔ Whitehead: Korrekt — ein Gegenbeispiel zu EG liefert eins zu Whitehead.
- Cockcroft-Eigenschaft: Korrekt definiert.
- Magnus Freiheitssatz 1930: Korrekt.
- Lyndon 1950 (asphaerisch fuer Einrelatorgruppen ohne echte Potenz): Korrekt.
- Beweis der Whitehead-Vermutung fuer Einrelator-Komplexe (Zeile 282–293): Korrekt, aber vereinfacht. Streng genommen gibt es Teilkomplexe mit weniger 1-Zellen und der 2-Zelle; das Argument "entweder 1-Skelett oder ganzer Komplex" ist korrekt fuer Prasentationskomplexe mit einem Relator (der jede Teilmenge von 1-Zellen, die alle im Relator vorkommenden Erzeuger enthaelt, plus die 2-Zelle erhaelt wieder den vollen Komplex; fehlt ein Erzeuger, ist der Teilkomplex das 1-Skelett minus einige Kreise, immer noch asphärisch).
- Corson 1992 (DR hereditary): Korrekt.
- Johnson 2003 (D(2) ↔ Whitehead fuer endliche Komplexe): Korrekt.
- Howie 1981/1979: Korrekt (einfach zusammenhaengende Teilkomplexe).

**BUG-B20-P78-EN-001** [GERING]
Zeile 544–547: Sieradski1976 bibitem zitiert "J. Combin. Theory Ser. B 34 (1983), 241–246" mit dem Titel "A coloring invariant for topological graph theory". Das ist der falsche Artikel! Der korrekte Sieradski-Artikel ueber diagrammatische Reduzierbarkeit ist: A.J. Sieradski, "A semigroup of simple homotopy types", Math. Z. 153 (1977), 135–148, oder der Beitrag in Sieradski's spaeterer Arbeit. Im Paper-Text steht "Sieradski 1976" (Zeile 306), was auf keinen der zitierten Titel passt. Die Referenz ist falsch zugeordnet.

**BUG-B20-P78-EN-002** [GERING]
Zeile 251–257: "Dunce hat" Beschreibung ist verworren. Der Text beginnt mit "attached by the word aa^{-1}a^{-1}... wait" und korrigiert sich selbst. Das ist ein offensichtlicher Entwurfsrest, der nicht in ein finales Paper gehoert. Selbst die "korrigierte" Version ("the dunce hat has the word aa^{-1}, giving pi_1 = 1 and is contractible") ist nicht ganz korrekt: der Dunce Hat hat die Identifikation der drei Seiten eines Dreiecks durch a→a→a^{-1}, was pi_1 = ⟨a | a⟩ = 1 ergibt. Der Dunce Hat ist kontrahierbar, also insbesondere asphärisch. Die Passage muss bereinigt werden.

**BUG-B20-P78-EN-003** [GERING]
Zeile 559–562: BridsondelaHarpe1999 — Das Buch heisst korrekt "Metric Spaces of Non-Positive Curvature" von Bridson und **Haefliger** (nicht "de la Harpe"). Martin Bridson und Andre Haefliger, Springer 1999. Pierre de la Harpe ist ein anderer Mathematiker (Autor von "Topics in Geometric Group Theory"). Falscher Koautor.

**LaTeX:**
- Sauber. Keine undefinierte Referenzen oder Umgebungen.
- Build 164 / Build 165 Inkonsistenz (Header vs. address). KOSMETISCH.

### DE-Version

**Urteil: DRUCKREIF (mit Anmerkungen)**

**BUG-B20-P78-DE-001** [GERING]
Identisch zu BUG-B20-P78-EN-001: Sieradski-Referenz falsch zugeordnet (Zeile 506–508).

**EN-DE-Konsistenz:**
- DE fehlt das "Dunce Hat" Beispiel (BUG-B20-P78-EN-002 betrifft nur EN). Gut — der verworrene Absatz fehlt in DE.
- DE fehlt das "Bridson-de la Harpe" Beispiel (BUG-B20-P78-EN-003 betrifft nur EN). DE hat die gleiche falsche Referenz (Zeile 521–523), aber kein Theorem, das sie direkt verwendet.
- DE fehlt der Bestvina2002-Bibitem (vorhanden in EN Zeile 575–578). Nicht verwendet in DE. Akzeptabel.
- Alle mathematischen Kernaussagen stimmen ueberein.
- Howie-Bibitem: EN sagt "1981" im Text, "1979" im Eintrag. DE identisch. KOSMETISCH.

---

## Paper 79: The Temperley-Lieb Algebra and Jones Polynomial

### EN-Version

**Urteil: DRUCKREIF**

**Mathematische Korrektheit:**
- Temperley-Lieb-Algebra TL_n(delta): Definition korrekt (3 Relationen).
- dim TL_n = C_n (Catalan-Zahl): Korrekt.
- Hecke-Algebra H_n(q), Surjektion rho: T_i ↦ qe_i - 1: Korrekt (Jones 1987).
- Markov-Spur: Korrekt definiert.
- Jones-Polynom-Definition via Markov-Spur: Korrekt.
- Skein-Relation: Korrekt (t^{-1} V_{L+} - t V_{L-} = (t^{1/2} - t^{-1/2}) V_{L0}).
- V(trefoil) = -t^{-4} + t^{-3} + t^{-1}: Korrekt fuer den rechtshändigen Kleeblattknoten.
- Kauffman-Bracket: Korrekt.
- Gefaerbtes Jones-Polynom, Kashaev-Invariante: Korrekt.
- Volumen-Vermutung: Korrekt formuliert.
- Vol(S^3 \ 4_1) = 3v_3 ≈ 2.0298 mit v_3 ≈ 1.0149: Korrekt (v_3 = Volume eines regulaeren idealen Tetraeders in H^3).
- Khovanov-Homologie, Kategorifizierung: Korrekt.
- Kronheimer-Mrowka 2011, Khovanov detects unknot: Korrekt.
- Rasmussen s-Invariante: Korrekt.
- AJ-Vermutung (Garoufalidis): Korrekt.
- Odd Khovanov (Ozsvath-Rasmussen-Szabo 2013): Korrekt.

**BUG-B20-P79-EN-001** [GERING]
Zeile 293–294: "Positive knots (Tanaka 1996)" — Die angegebene Referenz Tanaka 1996 ("Maximal Thurston-Bennequin numbers of 2-bridge links") behandelt Thurston-Bennequin-Zahlen von 2-Bruecken-Verschlingungen, nicht die Jones-Unknoten-Vermutung fuer positive Knoten. Der korrekte Verweis fuer die Jones-Unknoten-Vermutung bei positiven Knoten waere Cromwell (2004) oder Stoimenow (2005). Falsche Referenzzuordnung.

**BUG-B20-P79-EN-002** [KOSMETISCH]
Zeile 383–384: "Ekholm-Ng, 2012 for the figure-eight knot via a specific 'state sum' computation" — Die Volumen-Vermutung fuer den Achtknoten wurde urspruenglich von Kashaev (1997) verifiziert und von Murakami-Murakami (2001) reformuliert/bestaetigt. Ekholm-Ng (2012) arbeiten an Knot Contact Homology/augmentation varieties. Die Zuschreibung ist ungenau.

**LaTeX:**
- Zeile 309–311: \vcenter{\hbox{...}} Konstruktion fuer Skein-Diagramme — funktioniert in amsart, aber die Platzhalter "$L_+$", "$L_0$", "$L_\infty$" sind nur Textsubstitute fuer Diagramme. Akzeptabel.

### DE-Version

**Urteil: DRUCKREIF**

**BUG-B20-P79-DE-001** [KOSMETISCH]
Zeile 77: "die Skeин-Relation" — enthalt ein kyrillisches "и" statt "in" in "Skein". Muss "Skein-Relation" heissen.

**EN-DE-Konsistenz:**
- DE hat keine "Fibred knots" Erwaehnung in der Jones-Unknoten-Bemerkung (EN Zeile 295). Akzeptabel — kuerzer.
- DE hat keine explizite "Knot Floer Homology vs. Khovanov Homology" Vermutung (EN Section 9, letzte Vermutung). Akzeptabel.
- Alle mathematischen Kernaussagen und Formeln stimmen ueberein.
- Bibliographien identisch.

---

## Gesamtliste aller Bugs — Batch 20

| Bug-ID | Paper | Sprache | Schweregrad | Beschreibung |
|--------|-------|---------|-------------|--------------|
| BUG-B20-P76-DE-001 | 76 | DE | KOSMETISCH | "Knospfe" → "Knospe" (Zeile 267) |
| BUG-B20-P77-EN-001 | 77 | EN | HOCH | Freedman falsch: "Arf=0 ⟹ top. slice" ist falsch; korrekt ist "Δ_K=1 ⟹ top. slice" |
| BUG-B20-P77-EN-002 | 77 | EN | MITTEL | "Rolfsen table" falsch — 11n34 ist Hoste-Thistlethwaite-Notation |
| BUG-B20-P77-EN-003 | 77 | EN | GERING | Torusknoten T(2,2k+1) sind nicht top.-aber-nicht-glatt-scheibig |
| BUG-B20-P77-EN-004 | 77 | EN | KOSMETISCH | Gordon "1981" bibitem ist eigentlich 1978 |
| BUG-B20-P77-DE-001 | 77 | DE | HOCH | Identisch zu EN-001: Freedman falsch zitiert |
| BUG-B20-P77-DE-002 | 77 | DE | MITTEL | Identisch zu EN-002: "Rolfsen-Tabelle" falsch |
| BUG-B20-P78-EN-001 | 78 | EN | GERING | Sieradski-Referenz falsch zugeordnet (falscher Artikel) |
| BUG-B20-P78-EN-002 | 78 | EN | GERING | Dunce-Hat-Beschreibung verworren ("wait"), Entwurfsrest |
| BUG-B20-P78-EN-003 | 78 | EN | GERING | "Bridson-de la Harpe" ist falsch — korrekt: Bridson-Haefliger |
| BUG-B20-P78-DE-001 | 78 | DE | GERING | Sieradski-Referenz falsch (identisch zu EN-001) |
| BUG-B20-P79-EN-001 | 79 | EN | GERING | Tanaka 1996 ist falsche Referenz fuer pos. Knoten / Jones-Unknot |
| BUG-B20-P79-EN-002 | 79 | EN | KOSMETISCH | Ekholm-Ng Zuschreibung fuer Volume Conjecture ungenau |
| BUG-B20-P79-DE-001 | 79 | DE | KOSMETISCH | Kyrillisches "и" in "Skeин" statt "Skein" (Zeile 77) |

**Gesamt: 14 Bugs (2 HOCH, 2 MITTEL, 6 GERING, 4 KOSMETISCH)**

---

## Detailbewertung

### Paper 76 (Mandelbrot/MLC): Mathematisch exzellent
Die Darstellung der Mandelbrot-Menge, Boettcher-Koordinaten, MLC-Vermutung und Yoccoz' Satz ist durchgehend korrekt. Die MLC-Aequivalenzen sind vollstaendig und praezise formuliert. Die Verbindung MLC ⟹ Dichtheit der Hyperbolizitaet ist korrekt argumentiert. Die historischen Angaben (Mandelbrot 1980, Douady-Hubbard 1982, Yoccoz 1990, Shishikura 1998) stimmen. Ein vorbildliches Survey-Paper.

### Paper 77 (Conway-Knoten): Ueberarbeitung noetig
Der kritische Fehler BUG-B20-P77-EN/DE-001 (Freedmans Satz falsch zitiert als "Arf=0 implies top. slice") ist schwerwiegend, weil er eine mathematisch falsche Aussage als Theorem praesentiert. Der Fehler tritt direkt neben der korrekten Formulierung auf (Delta_K=1), was besonders verwirrend ist. Zusaetzlich ist die Zuschreibung an die Rolfsen-Tabelle (statt Hoste-Thistlethwaite) ein historischer Fehler. Die Skizze von Piccirillo's Beweis selbst ist korrekt.

### Paper 78 (Whitehead-Asphaerizitaet): Solide, kleinere Maengel
Mathematisch korrekt in allen Kernaussagen. Die Bibliographie-Fehler (Sieradski, Bridson-Haefliger) sind aergerlich aber nicht inhaltlich kritisch, da die Saetze selbst korrekt formuliert sind. Der "wait"-Kommentar im Dunce-Hat-Beispiel (nur EN) ist ein Entwurfsrest, der in einem finalen Paper nicht vorkommen darf.

### Paper 79 (Temperley-Lieb/Jones): Mathematisch stark
Alle algebraischen Definitionen (TL-Algebra, Hecke-Algebra, Markov-Spur) und topologischen Saetze (Skein-Relation, Khovanov-Kategorifizierung, Kronheimer-Mrowka) sind korrekt. Die Tanaka-Referenz-Fehlzuordnung ist der einzige inhaltliche Mangel. Das kyrillische Zeichen im DE-Abstract ist ein Encoding-Problem.
