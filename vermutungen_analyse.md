# Durchführbarkeitsanalyse offener mathematischer Vermutungen

> Stand: 2026-03-10
> Autor: Kurt Ingwer
> Bezug: [vermutungen.md](./vermutungen.md)

---

## Einleitung

Diese Analyse bewertet die offenen mathematischen Vermutungen aus `vermutungen.md` nach
ihrer Lösbarkeit im Rahmen des aktuellen mathematischen Wissens. Die Bewertung ist
naturgemäß subjektiv und spiegelt den Konsens der mathematischen Gemeinschaft wider.

Bereits bewiesene Vermutungen (z. B. Fermat, Catalan, Poincaré) werden hier **nicht**
erneut aufgeführt – sie sind in `vermutungen.md` mit ✓ gekennzeichnet.

### Bewertungsskala

Für jede Vermutung werden vier Dimensionen bewertet:

- **Werkzeuge**: Welche mathematischen Methoden und Theorien stehen bereits zur Verfügung?
- **Hindernisse**: Was fehlt noch, oder warum ist das Problem so schwer?
- **Wahrscheinlichkeit**: `Sehr hoch` / `Hoch` / `Mittel` / `Niedrig` / `Sehr niedrig` / `Definitiv unmöglich`
- **Zeithorizont**: `1–10 Jahre` / `10–50 Jahre` / `50–100 Jahre` / `>100 Jahre` / `Unbestimmbar`

---

## Kategorien

---

## A. Möglicherweise in den nächsten Jahrzehnten beweisbar (Hohe Wahrscheinlichkeit)

Vermutungen mit klaren Angriffspunkten, aktiver Forschung und substanziellen Teilresultaten.
Hier existieren konkrete Strategien, auch wenn der abschließende Beweis noch aussteht.

---

**Zwillingsprimzahl-Vermutung** (~1849, de Polignac)
- Werkzeuge: Siebmethoden (Brun, Selberg), Zhang-Maynard-Polymath-8-Resultat (Primzahllücken < 246), GPY-Siebtechnik, Goldbach-Strukturanalysen
- Hindernisse: Die Lücke zwischen "< 246" und "= 2" ist qualitativ – alle bekannten Siebmethoden scheitern am sogenannten *Parity Problem* von Selberg
- Wahrscheinlichkeit: Hoch (Maynard-Breakthrough 2014 hat das Feld fundamental verschoben; weitere Siebverbesserungen erscheinen realistisch)
- Zeithorizont: 10–50 Jahre (optimistisch 10–20 Jahre, falls Parity Problem überwunden wird)

---

**Goldbach-Vermutung** (1742)
- Werkzeuge: Kreismethode (Hardy-Littlewood-Vinogradov), Chen-Theorem (jede gerade Zahl = p + p·q), Siebtheorie, Dichteresultate, Numerische Verifikation bis 4×10¹⁸
- Hindernisse: Chen (1973) bewies p + (p₁·p₂), aber nicht p + p; das Parity Problem der Siebtheorie blockiert den letzten Schritt; Hardy-Littlewood-Singulärreihenmethode liefert Asymptotik, aber keinen exakten Beweis
- Wahrscheinlichkeit: Hoch (Teilresultate sehr stark; ternäre Goldbach durch Helfgott 2013 vollständig bewiesen)
- Zeithorizont: 10–50 Jahre

---

**Legendres Vermutung** (1798)
- Werkzeuge: Primzahlsatz im kurzen Intervall, Riemannscher Nullstellenfreiheitsbereich, Baker-Methode, Ingham-Abschätzungen
- Hindernisse: Benötigt effektive Schranken für Primzahlabstände der Form O(n^(1/2)); diese folgen bedingt aus der Riemann-Hypothese, aber unbedingt ist die Schranke noch bei O(n^(0.525))
- Wahrscheinlichkeit: Hoch (bedingt auf RH sofort lösbar; auch unbedingt substanzielle Fortschritte möglich)
- Zeithorizont: 10–50 Jahre (falls RH bewiesen wird: sofort; unbedingt: 20–50 Jahre)

---

**Brocard-Vermutung** (1904)
- Werkzeuge: Primzahlsatz, Bertrand-Postulat (bewiesen), Schranken für Primzahllücken, Computerverifikation
- Hindernisse: Benötigt schärfere Schranken für Primzahllücken als Brocard-Niveau (pₙ² bis pₙ₊₁²); die Primzahllücken im quadratischen Bereich sind numerisch unproblematisch, aber ein allgemeiner Beweis fehlt
- Wahrscheinlichkeit: Hoch (für praktische Zwecke sicher richtig; formaler Beweis folgt aus RH)
- Zeithorizont: 10–50 Jahre

---

**Erdős-Straus-Vermutung** (1948, Erdős)
- Werkzeuge: Modulare Arithmetik, Siebtheorie, Computerverifikation bis 10¹⁴, Dichte-Argumente (fast alle n haben Lösungen)
- Hindernisse: Ein allgemeiner Existenzbeweis für alle n fehlt; die Ausnahmemengen werden immer kleiner, aber ein vollständiges Dichteargument ist noch nicht rigoros
- Wahrscheinlichkeit: Hoch (numerisch überwältigend bestätigt; Dichte der Ausnahmen tendiert nachweislich zu 0)
- Zeithorizont: 10–50 Jahre

---

**Artin-Vermutung über primitive Wurzeln** (1927)
- Werkzeuge: GRH (Generalisierte Riemann-Hypothese) impliziert die Vermutung (Hooley 1967); Methoden der Siebtheorie, Dichte-Argumente
- Hindernisse: Vollständiger Beweis ohne GRH unbekannt; bedingt auf GRH vollständig bewiesen
- Wahrscheinlichkeit: Hoch (bedingt auf GRH vollständig; nach GRH-Beweis automatisch gelöst)
- Zeithorizont: 10–50 Jahre (abhängig von GRH)

---

**Novikov-Vermutung** (1965)
- Werkzeuge: K-Theorie, C*-Algebren, Kasparov-Theorie (KK-Theorie), Baum-Connes-Methoden (bewiesen für viele Gruppen), geometrische Gruppentheorie
- Hindernisse: Allgemeiner Beweis für alle diskreten Gruppen fehlt; einige exotische Gruppen (mit Kazhdan-Eigenschaft T und großen "Löchern") sind schwer zu behandeln
- Wahrscheinlichkeit: Hoch (Mehrheit der relevanten Gruppenklassen bewiesen; Gegenbeispiel gilt als sehr unwahrscheinlich)
- Zeithorizont: 10–50 Jahre

---

**Baum-Connes-Vermutung** (1982)
- Werkzeuge: K-Theorie von C*-Algebren, Kasparov-Produkt, equivariante KK-Theorie; bewiesen für a-T-menable Gruppen (Higson-Kasparov 2001)
- Hindernisse: Gegenbeispiele für die Version mit Koeffizienten wurden gefunden (Higson et al.); allgemeine Gültigkeit für alle lokalkompakten Gruppen fraglich
- Wahrscheinlichkeit: Mittel bis Hoch (ohne Koeffizienten vermutlich wahr; mit Koeffizienten möglicherweise falsch)
- Zeithorizont: 10–50 Jahre

---

**Langlands-Programm** (1967, Teilweise bewiesen)
- Werkzeuge: Automorphe Formen, Galois-Darstellungen, Shimura-Varietäten, geometrisches Langlands (Fargues-Scholze 2021); Wiles bewies Spezialfall, Lafforgue für Funktionenkörper
- Hindernisse: Globales Langlands für Zahlkörper und allgemeine reduktive Gruppen noch weit offen; Funktorialität-Vermutung fehlt
- Wahrscheinlichkeit: Hoch (fundamentale Korrektheit nicht bezweifelt; vollständiger Beweis erfordert Generation neuer Methoden)
- Zeithorizont: 50–100 Jahre (vollständige Formulierung und Beweis)

---

**Farrell-Jones-Vermutung** (1993)
- Werkzeuge: Controlled Topology, algebraische K- und L-Theorie, geometrische Gruppentheorie; bewiesen für hyperbische Gruppen, CAT(0)-Gruppen (Bartels-Lück)
- Hindernisse: Allgemeiner Beweis für alle Gruppen; Gruppen mit pathologischen geometrischen Eigenschaften bereiten Schwierigkeiten
- Wahrscheinlichkeit: Hoch (Konsens: wahrscheinlich wahr; aktive Forschungsfront)
- Zeithorizont: 10–50 Jahre

---

**Gyárfás-Vermutung** (1975)
- Werkzeuge: Ramsey-Theorie, Probabilistische Methode, χ-Begrenzungsmethoden, strukturelle Graphentheorie
- Hindernisse: Induktive Methoden scheitern; Ramsey-Argumente liefern nur exponentielle Schranken statt linearer
- Wahrscheinlichkeit: Hoch (numerisch und durch Spezialfälle stark unterstützt)
- Zeithorizont: 10–50 Jahre

---

**Seymour's 2. Nachbar-Vermutung** (1990)
- Werkzeuge: Turnier-Theorie, probabilistische Graphentheorie; Seymour bewies Existenz eines 1. Nachbarn; Dean-Conjecture für spezielle Fälle bewiesen
- Hindernisse: Globale Struktur von Turnieren ist schwer zu kontrollieren; kein allgemeines induktives Prinzip bekannt
- Wahrscheinlichkeit: Hoch (starke numerische Evidenz; keine Gegenbeispiele bekannt)
- Zeithorizont: 10–50 Jahre

---

## B. Möglicherweise lösbar, aber schwer (Mittlere Wahrscheinlichkeit)

Vermutungen, die grundlegend neue mathematische Ideen erfordern. Aktive Forschung
existiert, aber der Abstand zum Ziel ist noch sehr groß.

---

**Riemann-Hypothese** (1859, Millennium-Problem)
- Werkzeuge: Analytische Zahlentheorie, L-Funktionen, Random-Matrix-Theorie (Montgomerys Pair-Correlation), Spektraltheorie (Hilbert-Pólya-Ansatz), explizite Formeln, GUE-Statistik der Nullstellen
- Hindernisse: Kein bekannter Ansatz liefert auch nur eine einzige Aussage über *alle* nichttrivialen Nullstellen; die kritische Linie Re(s)=1/2 ist numerisch enorm gut bestätigt (>10¹³ Nullstellen), aber kein struktureller Beweis; der Nullstellenfreiheitsbereich ist viel zu schmal
- Wahrscheinlichkeit: Mittel (mathematischer Konsens: fast sicher wahr, aber Beweismethode fundamental unklar)
- Zeithorizont: 10–50 Jahre (optimistisch), >100 Jahre (realistisch)

---

**Birch & Swinnerton-Dyer-Vermutung** (1965, Millennium-Problem)
- Werkzeuge: L-Funktionen elliptischer Kurven, Kolyvagin-Methode (Rang ≤ 1 vollständig bewiesen), Iwasawa-Theorie, modulare Formen (nach Wiles/Taylor), Gross-Zagier-Formel
- Hindernisse: Rang ≥ 2 ist vollständig offen; die Verbindung zwischen analytischem Rang (Nullstellenordnung) und algebraischem Rang (rationale Punkte) ist nur für kleine Ränge rigoros
- Wahrscheinlichkeit: Mittel (bedingt auf RH starke Fortschritte möglich; Rang-2-Fall erfordert neue Kohomologiemethoden)
- Zeithorizont: 50–100 Jahre

---

**Hodge-Vermutung** (1950, Millennium-Problem)
- Werkzeuge: Algebraische de Rham-Kohomologie, Deligne-Hodge-Theorie, Motiventheorie (Grothendieck), Kategorientheorie
- Hindernisse: Bereits für glatte projektive Varietäten der Dimension 3 weitgehend offen; das Motiven-Framework ist der vermutlich richtige Rahmen, aber selbst unvollständig
- Wahrscheinlichkeit: Mittel (technisch extrem schwer; einige Experten zweifeln an der vollen Allgemeinheit)
- Zeithorizont: 50–100 Jahre

---

**Navier-Stokes Glattheit** (Millennium-Problem)
- Werkzeuge: Funktionalanalysis, Sobolev-Räume, Energie-Methoden, Mikrolokale Analysis, Besov-Räume; Leray (1934) bewies schwache Lösungen; Serrin-Kriterium für Regularität
- Hindernisse: Ob Blow-ups (Singularitäten) in 3D entstehen können, ist vollständig offen; die Energie-Kaskade in turbulenten Strömungen ist mathematisch unkontrolliert; numerische Simulation zeigt keine Singularitäten, beweist aber nichts
- Wahrscheinlichkeit: Mittel (viele Experten vermuten Glattheit; einige erwarten Gegenbeispiele mit Singularitäten)
- Zeithorizont: 50–100 Jahre

---

**Yang-Mills Massenücke** (Millennium-Problem)
- Werkzeuge: Eichfeldtheorie, Funktionalintegrale, nicht-abelsche Yang-Mills-Gleichungen, Confinement-Mechanismus, Gitter-QCD (numerische Evidenz für Massenücke)
- Hindernisse: Keine rigoros mathematische Formulierung von Yang-Mills als Quantenfeldtheorie auf $\mathbb{R}^4$ existiert; Konstruktive Quantenfeldtheorie scheitert in 4D; der Begriff "Mass Gap" selbst muss noch rigoros definiert werden
- Wahrscheinlichkeit: Mittel (physikalisch fast sicher korrekt; mathematisch fehlt der gesamte Rahmen)
- Zeithorizont: 50–100 Jahre (möglicherweise >100 Jahre)

---

**abc-Vermutung** (1985, Masser/Oesterlé)
- Werkzeuge: Analytische Zahlentheorie, elliptische Kurven, Mochizukis Inter-Universal Teichmüller Theory (IUT, 2012, ~500 Seiten, nicht allgemein akzeptiert), Mason-Stothers für Polynome (bewiesen)
- Hindernisse: Mochizukis Beweis wurde von Scholze-Stix (2018) mit einem konkreten Gegenargument zurückgewiesen; keine alternativen Beweisansätze in Sicht; die Vermutung impliziert dutzende andere Sätze (Fermat, Goldbach-Schwächungen)
- Wahrscheinlichkeit: Mittel (vermutlich wahr; aber Beweis fehlt de facto vollständig)
- Zeithorizont: 10–50 Jahre (falls neuer Ansatz entsteht; derzeit 50+ Jahre realistischer)

---

**Collatz-Vermutung (3n+1)** (1937, Collatz)
- Werkzeuge: Ergodentheorie, probabilistische Heuristiken, Tao-Dichte-Resultat 2019 (fast alle Zahlen erreichen eine Schranke nahe 1), Computerverifikation bis ~2,95×10²⁰
- Hindernisse: Die Iteration erzeugt ein chaotisches dynamisches System; keine algebraische Struktur greifbar; Tao bewies "fast alle" – aber "alle" ist fundamental schwerer; möglicherweise Gödel-relevant (siehe Abschnitt E)
- Wahrscheinlichkeit: Mittel (numerisch extrem gut bestätigt; aber formaler Beweis erscheint sehr weit entfernt)
- Zeithorizont: >100 Jahre (oder unbestimmbar)

---

**P ≠ NP** (1971, Millennium-Problem)
- Werkzeuge: Boolesche Schaltkreiskomplexität, Orakel-Methoden, Natürliche Beweise (Razborov-Rudich), Algebrisierungsmethode (Aaronson-Wigderson); Nicht-uniforme Schranken für bestimmte Schaltkreisklassen bewiesen
- Hindernisse: Razborov-Rudich zeigten, dass die meisten "natürlichen" Beweistechniken scheitern müssen; Algebrisierungsbarriere; es fehlt ein fundamentaler Ansatz zur Komplexitätsseparation; möglicherweise in ZFC unentscheidbar (Baker-Gill-Solovay 1975)
- Wahrscheinlichkeit: Mittel (Konsens: P≠NP, aber Beweis extrem schwer; P=NP fast ausgeschlossen)
- Zeithorizont: Unbestimmbar (möglicherweise in ZFC unbeweisbar)

---

**Schanuel-Vermutung** (1960er)
- Werkzeuge: Transzendenztheorie (Baker-Methode), algebraische Unabhängigkeit (Nesterenko), p-adische Analysis; impliziert Lindemann-Weierstrass und Baker als Spezialfälle
- Hindernisse: Selbst einfachste Spezialfälle (z.B. Transzendenz von e+π) sind offen; die Vermutung ist so stark, dass fast keine der bekannten Methoden reicht
- Wahrscheinlichkeit: Mittel (mathematisch fast sicher korrekt; Beweis scheint Jahrhunderte entfernt)
- Zeithorizont: >100 Jahre

---

**Jacobian-Vermutung** (dim ≥ 2, 1939/1970)
- Werkzeuge: Algebraische Geometrie, Jet-Räume, Motiventheorie, Nash-Blow-ups; für n=1 trivial; numerische Algorithmen können Gegenbeispiele suchen; eng verwandt mit Dixmier-Vermutung
- Hindernisse: Dutzende von Beweisversuchen waren fehlerhaft; die Struktur polynomialer Automorphismen in ≥2 Dimensionen ist fundamental unverstanden; das Problem ist äquivalent zu mehreren anderen offenen Vermutungen
- Wahrscheinlichkeit: Mittel (weitgehend als wahr angesehen; aber Beweis fehlt für alle Dimensionen)
- Zeithorizont: 10–50 Jahre

---

**Elliptische-Kurven-Rang-Vermutung** (~1970)
- Werkzeuge: BSD-Vermutung (falls bewiesen), Kolyvagin-Systeme, Iwasawa-Theorie, explizite algebraische Konstruktionen (Elkies hat Kurven mit Rang ≥ 28 gefunden)
- Hindernisse: Ob der Rang unbeschränkt ist, hängt von tiefen Eigenschaften der L-Funktionen ab; keine konstruktive Methode für arbiträr hohe Ränge bekannt
- Wahrscheinlichkeit: Mittel (unbeschränkte Ränge mathematisch plausibel; Beweis in beiden Richtungen schwer)
- Zeithorizont: 50–100 Jahre

---

**Fontaine-Mazur-Vermutung** (1995)
- Werkzeuge: p-adische Hodge-Theorie (Fontaine), Galois-Darstellungen, potenzielle Automorphie (Taylor), Langlands-Korrespondenz
- Hindernisse: Vollständiger Beweis erfordert vollständiges Langlands; nur Spezialfälle zugänglich; die geometrische Seite der Vermutung (Berger-Colmez) ist noch offen
- Wahrscheinlichkeit: Mittel (technisch eng verwandt mit bereits bewiesenen Resultaten)
- Zeithorizont: 10–50 Jahre

---

**Hadwiger-Vermutung** (1943, für n≥7)
- Werkzeuge: Graphentheorie, Robertson-Seymour-Theorie (Graph Minors), Hadwiger-Beweis für n≤6; Vierfarbensatz (n=5 folgt daraus), Strukturelle Graphendekomposition
- Hindernisse: Für n≥7 grundlegend offen; jede Verbesserung erfordert neue Struktur-Theorie für Minor-freie Graphen
- Wahrscheinlichkeit: Mittel (für kleine n gut bestätigt; allgemeiner Fall erfordert neue Ideen)
- Zeithorizont: 50–100 Jahre

---

**Hadwiger-Nelson-Problem** (1950)
- Werkzeuge: Graphfärbung, Maßtheorie, Algebraische Graphentheorie; de Grey bewies 2018 χ≥5; bekannt: 4 ≤ χ ≤ 7
- Hindernisse: Beweis von χ=5, 6 oder 7 erfordert entweder konstruktive Graphen oder Unmöglichkeitsbeweise für Färbungen; die Unabhängigkeitszahl messbarer Graphen ist noch nicht charakterisiert
- Wahrscheinlichkeit: Mittel (de Grey 2018 war ein überraschender Durchbruch; weitere Schritte schwerer)
- Zeithorizont: 10–50 Jahre

---

**Mandelbrot-Lokal-Zusammenhangsvermutung (MLC)** (1982)
- Werkzeuge: Komplexe Dynamik, Thurston-Rigidität, Parabolisch-Implosion, Yoccoz-Puzzle (bewiesen für Misiurewicz- und parabolische Punkte), Renormalisierungstheorie
- Hindernisse: Unendlich oft renormierbare Parameter sind fundamental schwer; die Renormalisierungsoperator-Dynamik ist chaotisch und schlecht kontrollierbar
- Wahrscheinlichkeit: Mittel (Konsens: wahrscheinlich wahr; aktive Forschung in komplexer Dynamik)
- Zeithorizont: 10–50 Jahre

---

**Grothendieck-Standard-Vermutungen** (1969)
- Werkzeuge: Algebraische K-Theorie, Motiventheorie, Weil-Kohomologie (Deligne bewies Weil-Vermutungen); Lefschetz-Standard-Vermutung für einige Klassen bewiesen
- Hindernisse: Die Hodge-Standard-Vermutung ist stärker als Hodge-Vermutung; Motive als abelsche Kategorie noch nicht konstruiert; grundlegendste Vermutung (Lefschetz) für allgemeine Varietäten offen
- Wahrscheinlichkeit: Niedrig bis Mittel (teilweise wahr, teilweise möglicherweise falsch in voller Allgemeinheit)
- Zeithorizont: >100 Jahre

---

**Smooth-Poincaré-Vermutung dim=4** (1982)
- Werkzeuge: Donaldson-Invarianten, Seiberg-Witten-Invarianten, Freedman (topologische Version bewiesen 1982); Kirby-Siebenmann-Kriterium
- Hindernisse: Dimension 4 ist wegen der Fehlen glatter Strukturklassifikation einzigartig schwer; "exotic R⁴" und andere exotische 4-Mannigfaltigkeiten existieren; Methoden aus dim 3 und dim ≥ 5 sind hier nicht anwendbar
- Wahrscheinlichkeit: Niedrig bis Mittel (sowohl wahre als auch falsche Versionen von Experten für möglich gehalten)
- Zeithorizont: Unbestimmbar

---

**Kontsevich-Integral als vollständige Invariante** (1993)
- Werkzeuge: Quantengruppen, Vassiliev-Invarianten, Bar-Natan-Algebra der Chord-Diagramme, Knotentheorie
- Hindernisse: Ob das Kontsevich-Integral Knoten vollständig klassifiziert, ist offen; die Berechnung für nicht-triviale Knoten ist extrem aufwendig; kein bekannter Zusammenhang mit geometrischen Knoteninvarianten
- Wahrscheinlichkeit: Niedrig bis Mittel (technisch sehr schwer; keine klare Strategie)
- Zeithorizont: 50–100 Jahre

---

**Cramér-Vermutung** (1936)
- Werkzeuge: Probabilistische Heuristik (Cramér-Modell), Primzahlsatz in kurzen Intervallen, Maier-Phänomen (zeigt, dass Cramér-Modell zu naiv ist), Granville-Korrekturen
- Hindernisse: Das Cramér-Modell ist nachweislich falsch (Maier 1985); die richtige Konstante bleibt unklar; ein Beweis erfordert tiefe Kontrolle über Primzahllücken
- Wahrscheinlichkeit: Mittel (die qualitative Aussage ist wahrscheinlich korrekt; die genaue Konstante fraglich)
- Zeithorizont: >100 Jahre

---

**Waring-Problem: G(k) für alle k** (Verfeinerungen)
- Werkzeuge: Hardy-Littlewood-Kreismethode, Waring-Waring-Vinogradov-Methode; G(2)=4 (Lagrange), G(3) zwischen 4 und 7; effektive Schranken durch Wooley (Efficient Congruencing)
- Hindernisse: Die exakten Werte von G(k) für k ≥ 3 erfordern Kombinationen aus Siebtheorie und Kreismethode; G(3)=7 ist vermutet aber nicht bewiesen
- Wahrscheinlichkeit: Mittel (schrittweise Fortschritte wahrscheinlich; vollständige Lösung für alle k schwer)
- Zeithorizont: 10–50 Jahre (für einzelne k); >100 Jahre (für alle k)

---

**Beal-Vermutung** (1997)
- Werkzeuge: Modulare Formen (nach Wiles-Methode für Fermat), Frey-Kurven, Galois-Darstellungen; impliziert verallgemeinerte Fermat-Vermutung
- Hindernisse: Die Methoden aus Wiles' Fermat-Beweis decken nicht alle Exponenten-Kombinationen ab; die Fälle (2,3,7), (2,3,8) usw. erfordern separate Behandlung; einige Spezialfälle gelöst (Darmon-Granville)
- Wahrscheinlichkeit: Mittel (vermutlich wahr; Beweis erfordert fundamentale Erweiterung der Fermat-Methoden)
- Zeithorizont: 50–100 Jahre

---

**Gromov-Vermutung** (1983)
- Werkzeuge: Riemannsche Geometrie, Füllungsradius (Gromov), systolische Geometrie, Waist-Inequalities
- Hindernisse: Die genaue Beziehung zwischen Füllungsradius und Volumen in allgemeinen Dimensionen ist offen; Gromov bewies die Ungleichung selbst, aber die scharfe Konstante fehlt
- Wahrscheinlichkeit: Mittel
- Zeithorizont: 10–50 Jahre

---

**Lehmer-Vermutung** (1933)
- Werkzeuge: Algebraische Zahlentheorie, Mahler-Maß, Smyth-Resultat (für nicht-reziproke Polynome bewiesen), Kronecker-Theorem, Schur-Siegel-Smyth-Pfad
- Hindernisse: Reziproke Polynome sind fundamental schwerer; kein Weg bekannt, das Minimum des Mahler-Maßes global zu charakterisieren
- Wahrscheinlichkeit: Mittel (Lehmer-Polynomwert 1.1762... ist fast sicher das Minimum)
- Zeithorizont: 50–100 Jahre

---

**Perfekte Zahlen sind gerade** (~350 v. Chr.)
- Werkzeuge: Euler bewies: jede gerade perfekte Zahl hat Form 2^(p-1)·(2^p - 1); ungerade perfekte Zahlen hätten mindestens 9 verschiedene Primfaktoren (Nielsen 2015), mindestens 1500 Stellen (Ochem-Rao 2012)
- Hindernisse: Keines der bekannten Ausschlusskritierien schließt alle Fälle aus; die Existenz ist in ZFC konsistent; kein allgemeines algebraisches Argument bekannt
- Wahrscheinlichkeit: Niedrig bis Mittel (Konsens: vermutlich wahr, aber Beweis extrem schwer)
- Zeithorizont: >100 Jahre

---

**Es gibt unendlich viele vollkommene Zahlen** (~350 v. Chr.)
- Werkzeuge: Äquivalent zu: Es gibt unendlich viele Mersenne-Primzahlen; GIMPS-Projekt (numerische Suche); Heuristiken sprechen für unendlich viele
- Hindernisse: Kein analytischer Beweis für unendlich viele Mersenne-Primzahlen; scheitert an denselben Barrieren wie Artin-Vermutung
- Wahrscheinlichkeit: Mittel (heuristische Argumente stark; Beweis fehlt)
- Zeithorizont: Unbestimmbar

---

**Andrews-Curtis-Vermutung** (1965)
- Werkzeuge: Kombinatorische Gruppentheorie, Präsentation von Gruppen, Computer-gestützte Suche nach AC-Äquivalenzen; engverbunden mit Whitehead-Vermutung
- Hindernisse: Gegenbeispiele werden vermutet aber nicht gefunden; kein Invariante bekannt, die AC-Äquivalenz vollständig charakterisiert
- Wahrscheinlichkeit: Niedrig bis Mittel (einige Experten vermuten, dass Gegenbeispiele existieren)
- Zeithorizont: 50–100 Jahre

---

**Whitehead-Vermutung (Aspährizität)** (1941)
- Werkzeuge: Algebraische Topologie, Überlagerungsräume; für endliche Komplexe offen; Zeeman-Vermutung (in dim ≥ 5 bewiesen) impliziert Whitehead-Vermutung
- Hindernisse: Die 2-dimensionale CW-Komplex-Theorie ist fundamental unvollständig; kein allgemeines Kriterium für Aspährizität von Unterkomplexen bekannt
- Wahrscheinlichkeit: Mittel
- Zeithorizont: 10–50 Jahre

---

**Conway-Vermutung (Knotentheorie)** (1970, teils offen)
- Werkzeuge: Knotentheorie, Skein-Relationen, Jones-Polynom, HOMFLY-Polynom, Kategorifizierung (Khovanov-Homologie)
- Hindernisse: Verschiedene offene Teile; die Frage, ob bestimmte Knoteninvarianten vollständig sind, ist offen
- Wahrscheinlichkeit: Mittel (je nach Teilproblem variiert die Aussicht)
- Zeithorizont: 10–50 Jahre

---

**Erdős-Szemerédi-Vermutung** (1983)
- Werkzeuge: Additive Kombinatorik (Fourier-Analyse, Graphentheorie, Freiman-Ruzsa-Theorie); Solymosi bewies partielle Resultate (f(A) ≥ |A|^(4/3)); Elekes-Rónyai für Polynome
- Hindernisse: Das volle Ergebnis |A+A| + |A·A| ≥ |A|^(2-ε) ist weit entfernt; die Interaction zwischen Addition und Multiplikation ist schlecht verstanden
- Wahrscheinlichkeit: Mittel (Teilresultate kommen schrittweise; vollständiger Beweis erfordert neue Ideen)
- Zeithorizont: 10–50 Jahre

---

**Bateman-Horn-Vermutung** (1962)
- Werkzeuge: Kreismethode, Siebtheorie; bedingt auf GRH Teilergebnisse; quantitativ für spezielle Polynome bekannt (Lineare: Dirichlet; Quadratisch: teilweise)
- Hindernisse: Gleichzeitige Primalität mehrerer Polynome ist durch Parity Problem der Siebtheorie blockiert
- Wahrscheinlichkeit: Mittel (heuristische Evidenz überwältigend; Beweis fundamental schwer)
- Zeithorizont: >100 Jahre

---

**Bunyakovsky-Vermutung** (1857)
- Werkzeuge: Siebtheorie, Dirichlet-Reihen; Dirichlets Satz (Spezialfall für lineare Polynome, bewiesen 1837); Bateman-Horn als quantitative Version
- Hindernisse: Parity Problem der Siebtheorie; für quadratische Polynome (z. B. n²+1) vollständig offen
- Wahrscheinlichkeit: Mittel (fast sicher korrekt; Beweis erfordert Revolution in der Siebtheorie)
- Zeithorizont: 50–100 Jahre

---

**Lonely Runner-Vermutung** (1967, Wills, für ≥7 Läufer)
- Werkzeuge: Diophantische Approximation, Maßtheorie, Fourier-Analyse; bewiesen für ≤6 Läufer; probabilistische Methoden geben richtiges asymptotisches Bild
- Hindernisse: Jeder neue Läufer erfordert fundamental andere Argumente; keine Induktionsmethode bekannt
- Wahrscheinlichkeit: Mittel (numerisch stark bestätigt; formaler Beweis schwer)
- Zeithorizont: 10–50 Jahre

---

**Graceful Tree-Vermutung** (1966, Kotzig/Rosa)
- Werkzeuge: Graphentheorie, Computational Search; alle Bäume bis ~36 Knoten computergestützt verifiziert; verschiedene Baum-Klassen (Pfade, Sterne, Caterpillar-Graphen) bewiesen
- Hindernisse: Kein allgemeines strukturelles Argument; für allgemeine Bäume kein Beweisansatz
- Wahrscheinlichkeit: Mittel (wahrscheinlich wahr; Beweis könnte durch enumerative Methoden möglich sein)
- Zeithorizont: 10–50 Jahre

---

**Erdős–Gallai-Vermutung (Primzahlen als Summen)** (1959)
- Werkzeuge: Additive Zahlentheorie, Siebmethoden, Cramér-Typ-Schranken
- Hindernisse: Die genaue Struktur additiver Zerlegungen für Primzahlen ist schlecht kontrollierbar
- Wahrscheinlichkeit: Mittel
- Zeithorizont: 10–50 Jahre

---

**Turán-Vermutung über Nullstellen der ζ-Funktion** (1948)
- Werkzeuge: Analytische Zahlentheorie, Nullstellendichte-Sätze, Montgomery-Paar-Korrelation; Turán bewies Teilresultate über Nullstellenfreiheitsbereiche
- Hindernisse: Eng verwandt mit RH; ohne RH sind die Turán-Schranken sehr schwer zu verbessern
- Wahrscheinlichkeit: Mittel (Lösung wahrscheinlich nach RH-Beweis)
- Zeithorizont: Abhängig von RH

---

**Uniformisierungsvermutung für höhere Dimensionen** (diverse)
- Werkzeuge: Differentialgeometrie, Ricci-Fluss (Perelman), Kähler-Geometrie; in 2D vollständig (Poincaré/Köbe); Hamilton-Ricci-Fluss in 3D (Perelman)
- Hindernisse: Dimension ≥ 4 hat keine Analogie zum 3D-Ricci-Fluss-Beweis; Komplexität der holomorphen Struktur in höheren Dimensionen
- Wahrscheinlichkeit: Niedrig bis Mittel (je nach Präzisierung der Vermutung)
- Zeithorizont: >100 Jahre

---

**Hartshorne-Vermutung** (1970)
- Werkzeuge: Algebraische Geometrie, Kodaira-Einbettung, Positivity-Theorie (Ampleness, Nefness), Mori-Programm
- Hindernisse: Charakterisierung von PⁿC unter algebraischen Mannigfaltigkeiten durch Nef-Tangentialbündel; komplizierte Gegenbeispiele in positiver Charakteristik
- Wahrscheinlichkeit: Niedrig bis Mittel
- Zeithorizont: 50–100 Jahre

---

**Gromov-Vermutung (Füllungsradius)** (1983)
- Werkzeuge: Systolische Geometrie, Waist-Inequalities (Gromov), Fillingsradius-Theorie
- Hindernisse: Scharfe Konstanten in Füllungsungleichungen für allgemeine Riemannsche Mannigfaltigkeiten fehlen
- Wahrscheinlichkeit: Mittel
- Zeithorizont: 10–50 Jahre

---

**Temperley-Lieb-Vermutungen** (1971, teils offen)
- Werkzeuge: Kombinatorik, Algebra, Quantengruppen, Darstellungstheorie
- Hindernisse: Je nach Teilproblem verschieden; Verbindung zur statistischen Mechanik erschwert rein algebraische Ansätze
- Wahrscheinlichkeit: Mittel
- Zeithorizont: 10–50 Jahre

---

**Donaldson-Geometrie-Vermutung** (1983, teils bewiesen)
- Werkzeuge: Seiberg-Witten-Theorie, Gromov-Witten-Invarianten, Symplektische Topologie, Floer-Homologie
- Hindernisse: Offene Teile betreffen die vollständige Klassifikation glatter 4-Mannigfaltigkeiten; exotische Strukturen in Dimension 4 machen dies fundamental anders als andere Dimensionen
- Wahrscheinlichkeit: Mittel (Teile werden schrittweise bewiesen)
- Zeithorizont: 50–100 Jahre

---

**Freiman-Struktursatz Verbesserungen** (1964, teils offen)
- Werkzeuge: Additive Kombinatorik, Ruzsa-Covering-Lemma, Sanders-Quantifizierung (2012), Polynomial-Schranken (Gowers)
- Hindernisse: Die optimale Konstante im Freiman-Ruzsa-Satz ist unbekannt; polynomiale vs. exponentielle Schranken sind offen
- Wahrscheinlichkeit: Hoch (schrittweise Verbesserungen wahrscheinlich)
- Zeithorizont: 10–50 Jahre

---

**Erdős–Ko–Rado Verallgemeinerungen** (1961, teils offen)
- Werkzeuge: Extremale Mengenlehre, Shifting-Methode, Algebraische Methoden; Grundversion bewiesen; Verallgemeinerungen auf Permutationen (Ellis-Friedgut-Pilpel 2011)
- Hindernisse: Allgemeine intersektive Familien in abstrakten Strukturen; Extremalfälle für nicht-uniforme Hypergraphen
- Wahrscheinlichkeit: Hoch (schrittweise Fortschritte gut etabliert)
- Zeithorizont: 10–50 Jahre

---

## C. Extrem schwer / möglicherweise unentscheidbar (Niedrige Wahrscheinlichkeit)

Vermutungen, die möglicherweise unbeweisbar sind oder fundamentale neue mathematische
Strukturen erfordern, die noch nicht existieren.

---

**P ≠ NP** (1971) — *Hier nochmals als Grenzfall*
- Werkzeuge: Schaltkreiskomplexität, Orakel-Trennungen, Algebrisierung; Baker-Gill-Solovay zeigten: Orakel-Methoden reichen nicht
- Hindernisse: Razborov-Rudich-Barriere schließt "natürliche Beweise" aus; Algebrisierungsbarriere schließt relativisierbare Methoden aus; GCT (Geometric Complexity Theory, Mulmuley) ist der vielversprechendste Ansatz, aber Jahrzehnte entfernt
- Wahrscheinlichkeit: Niedrig (Lösung in naher Zukunft); die Frage selbst könnte in ZFC unentscheidbar sein
- Zeithorizont: Unbestimmbar

---

**Collatz-Vermutung** (1937) — *Hier nochmals als Grenzfall*
- Werkzeuge: Ergodentheorie, Tao-Methoden (2019), Probabilistische Modelle
- Hindernisse: Conway bewies, dass verallgemeinerte Collatz-Probleme unentscheidbar sind; die spezifische 3n+1-Version könnte ZFC-unentscheidbar sein
- Wahrscheinlichkeit: Sehr niedrig (Beweis in absehbarer Zeit)
- Zeithorizont: Unbestimmbar

---

**Grothendieck-Standard-Vermutungen (vollständig)** (1969)
- Werkzeuge: Motiventheorie; aber diese ist selbst noch nicht vollständig konstruiert
- Hindernisse: Die Hodge-Standard-Vermutung ist stärker als die Hodge-Vermutung selbst; Gegenbeispiele in positiver Charakteristik möglicherweise vorhanden
- Wahrscheinlichkeit: Sehr niedrig (vollständige Lösung)
- Zeithorizont: >100 Jahre

---

**Schanuel-Vermutung** (1960er)
- Werkzeuge: Baker-Theorie, Nesterenko-Methoden; selbst e+π-Transzendenz ist unbekannt
- Hindernisse: Impliziert so viele ungelöste Probleme, dass ein Beweis eine Revolution darstellen würde
- Wahrscheinlichkeit: Sehr niedrig (vollständiger Beweis)
- Zeithorizont: >100 Jahre

---

**Jacobian-Vermutung (allgemein, dim ≥ 3)** (1939/1970)
- Werkzeuge: Algebraische Geometrie; aber keine der Beweisstrategien für dim=2 übertragen sich auf dim≥3
- Hindernisse: Die Vermutung könnte in Dimension ≥ 3 falsch sein; Gegenbeispiele werden nicht ausgeschlossen
- Wahrscheinlichkeit: Niedrig (vollständiger allgemeiner Beweis)
- Zeithorizont: 50–100 Jahre

---

**Mochizuki's IUT / abc-Vermutung via neuer Methoden** (kontrovers)
- Werkzeuge: IUT (Mochizuki), aber nicht allgemein akzeptiert; kein alternativer Ansatz
- Hindernisse: Keine zugängliche Beweisstrategie; Scholze-Stix-Widerspruch noch nicht aufgelöst
- Wahrscheinlichkeit: Niedrig (kurzfristiger anerkannter Beweis)
- Zeithorizont: Unbestimmbar

---

**Smooth-Poincaré-Vermutung dim=4** (1982)
- Werkzeuge: Donaldson, Seiberg-Witten; aber exotische 4-Mannigfaltigkeiten machen den Ausgang unklar
- Hindernisse: Kann in beide Richtungen gehen; keine klare Strategie für Beweis oder Widerlegung
- Wahrscheinlichkeit: Niedrig (in absehbarer Zeit)
- Zeithorizont: Unbestimmbar

---

## D. Formal unwiderlegbar/unentscheidbar (Definitiv nicht durch konventionelle Mittel lösbar)

Vermutungen, die in ZFC wahrscheinlich unentscheidbar sind oder deren Lösung
übergeordnete Axiomensysteme erfordern würde.

> **Wichtig**: Keine der folgenden Vermutungen ist mit Sicherheit als ZFC-unentscheidbar
> bewiesen – das wäre selbst ein großes mathematisches Ergebnis. Die Einordnung reflektiert
> den Verdacht und die Analogie zu bekannten unentscheidbaren Sätzen.

---

**Verallgemeinerte Collatz-Vermutungen** (Conway 1987)
- Conway bewies rigoros, dass verallgemeinerte lineare Kongruenz-Iterationen (*Fractran*) unentscheidbar sind. Die spezifische 3n+1-Version ist vermutlich lösbar, aber es gibt keine Grenze, die sie von unentscheidbaren Varianten trennt.
- Werkzeuge: Turing-Maschinen-Simulation via Fractran; Gödel-Unvollständigkeitssatz
- Hindernisse: Jeder Beweis müsste das unendliche dynamische Verhalten einer nicht-periodischen Iteration kontrollieren – ähnlich wie beim Halteproblem
- Wahrscheinlichkeit: Die *spezifische* 3n+1-Frage ist möglicherweise in ZFC unentscheidbar
- Zeithorizont: Unbestimmbar

---

**P ≠ NP** (als mögliche ZFC-Unentscheidbarkeit)
- Scott Aaronson und andere haben diskutiert, dass P≠NP in ZFC unentscheidbar sein könnte, analog zu Gödels Unvollständigkeitssatz: Ein System, das P≠NP nicht beweisen kann, könnte konsistent mit P=NP oder P≠NP sein.
- Werkzeuge: Logische Relativierbarkeit, Orakel-Trennungen (Baker-Gill-Solovay), Unprovability-Resultate à la Paris-Harrington
- Hindernisse: Kein direkter Beweis der Unentscheidbarkeit; aber alle bekannten Beweisbarrieren schließen Standardmethoden aus
- Wahrscheinlichkeit: Möglicherweise in ZFC unentscheidbar; falls entscheidbar, dann P≠NP fast sicher
- Zeithorizont: Unbestimmbar

---

---

## E. Vermutungen, die formal nicht lösbar sein könnten

### Gödels Unvollständigkeitssatz und seine Relevanz

Kurt Gödel bewies 1931, dass jedes hinreichend mächtige konsistente formale System
(wie ZFC = Zermelo-Fraenkel-Mengenlehre mit Auswahlaxiom) **wahre Aussagen enthält,
die im System nicht beweisbar sind**. Dies hat direkte Konsequenzen für mathematische
Vermutungen.

**Zwei Formen der Unbeweisbarkeit:**

1. **Syntaktische Unentscheidbarkeit**: Die Vermutung ist weder beweisbar noch widerlegbar in ZFC. Ein Modell von ZFC, in dem sie wahr ist, und eines, in dem sie falsch ist, können beide konsistent sein.

2. **Wahre aber unbeweisbare Sätze**: Die Vermutung ist tatsächlich wahr (in der "Standardinterpretation"), aber kein Beweis in ZFC existiert.

### Bekannte ZFC-unentscheidbare Sätze (als Analogie)

- **Kontinuumshypothese (CH)**: Von Cantor 1878 aufgestellt; Cohen und Gödel bewiesen: CH ist in ZFC unabhängig (weder beweisbar noch widerlegbar).
- **Paris-Harrington-Satz**: Ein Satz über Ramsey-Theorie, der wahr, aber in Peano-Arithmetik nicht beweisbar ist.
- **Goodstein-Satz**: Wahr, aber in PA nicht beweisbar.

### Welche offenen Vermutungen könnten betroffen sein?

| Vermutung | Begründung des Verdachts |
|-----------|--------------------------|
| **Collatz-Vermutung** | Conway bewies Unentscheidbarkeit verallgemeinerter Varianten; die 3n+1-Frage könnte ähnliche logische Komplexität haben |
| **P ≠ NP** | Alle bekannten Beweismethoden scheitern an prinzipiellen Barrieren; Baker-Gill-Solovay zeigten Orakel-Relativität; möglicherweise analog zu CH |
| **Jacobian-Vermutung (dim≥3)** | Äquivalent zu Problemen über polynomial Automorphismen; könnte unentscheidbare Instanzen enthalten |
| **Schanuel-Vermutung** | Impliziert so viele transzendente Aussagen, dass ein ZFC-Beweis möglicherweise Axiome erfordern würde, die über ZFC hinausgehen |

### Was Gödel NICHT bedeutet

- **Goldbach und Zwillingsprimzahlen sind NICHT Gödel-relevant** in dem Sinne, dass sie klare arithmetische Bedeutung in ZFC haben und prinzipiell in ZFC entscheidbar sein können. Dass wir keinen Beweis kennen, liegt an mathematischer Schwierigkeit, nicht an logischer Unentscheidbarkeit.
- **Die meisten Vermutungen in der Liste sind in ZFC entscheidbar** – wir wissen es nur noch nicht.
- **Unentscheidbarkeit ist keine Niederlage**: Ein Beweis, dass P≠NP in ZFC unentscheidbar ist, wäre selbst ein Millennium-wüdiges Ergebnis.

### Empfehlung für die Forschungsstrategie

Wenn eine Vermutung alle konventionellen Beweisbarrieren überwunden hat, lohnt es sich,
zu untersuchen, ob sie möglicherweise:
1. Neue Axiome erfordert (wie das *Large Cardinal*-Programm)
2. Mit bekannten unentscheidbaren Problemen äquivalent ist
3. Durch Meta-Mathematik (Beweistheorie) charakterisiert werden kann

---

*Letzte Änderung: 2026-03-10*
*Autor: Kurt Ingwer*
*Quellen: Clay Mathematics Institute, AMS Bulletin, MathOverflow-Konsens, Einzelarbeiten zu den jeweiligen Vermutungen*
