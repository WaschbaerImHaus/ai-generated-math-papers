# Offene Probleme der Geometrie und Algebraischen Geometrie — Gruppe B

**Autor:** Michael Fuhrmann
**Erstellt:** 2026-03-12
**Letzte Änderung:** 2026-03-12
**Zweck:** Umfassende Recherche zu 6 zentralen offenen Problemen in Differentialtopologie, Kähler-Geometrie, algebraischer Geometrie, Knotentheorie und Systolischer Geometrie.

---

## Problem 1: Glattes Poincaré-Problem in Dimension 4 — Hat S⁴ eine exotische differenzierbare Struktur?

### Hintergrund und Problemstellung

Die glatte Poincaré-Vermutung in Dimension 4 (SPC4) fragt: Ist jede geschlossene glatte 4-Mannigfaltigkeit, die homotopieäquivalent zu S⁴ ist, auch diffeomorph zu S⁴? Dies ist das einzige verbliebene offene Problem der Poincaré-Vermutung in allen Dimensionen — in Dimensionen ≥5 wurde die Aussage von Smale (1961) bewiesen, die topologische Version in Dimension 4 von Freedman (1982), aber die glatte Variante in Dimension 4 ist vollständig offen.

### Gescheiterte Ansätze: Gluck-Twists

Herman Gluck beschrieb 1962 eine Klasse von Konstruktionen: Man schneide aus S⁴ eine Tubenumgebung einer eingebetteten 2-Sphäre heraus (deren Randmannigfaltigkeit S² × S¹ ist) und klebe sie mit dem nicht-trivialen Diffeomorphismus von S² × S¹ (Generator von π₁(SO(3)) = ℤ/2ℤ) zurück. Das Resultat ist stets homöomorph zu S⁴, aber die Frage, ob es diffeomorph zu S⁴ ist, blieb ungeklärt.

Das zentrale Problem: Es gibt keinen bekannten topologischen Invarianten, der verschiedene glatte Strukturen auf S⁴ unterscheiden könnte. Seiberg-Witten-Invarianten — die für Mannigfaltigkeiten mit b₂⁺ ≥ 2 äußerst effektiv exotische Strukturen nachweisen — sind für S⁴ (wo b₂ = 0 gilt) identisch null und liefern keinerlei Information. Dies ist keine Schwäche des Instruments, sondern ein fundamentales Hindernis: Für Mannigfaltigkeiten mit b₂⁺ = 0 können Seiberg-Witten-Invarianten prinzipiell keine exotischen Strukturen detektieren.

### Cappell-Shaneson-Sphären: Kandidaten, die sich als Standard erwiesen

Cappell und Shaneson konstruierten 1976 eine unendliche Familie von Homotopie-4-Sphären aus Matrizengruppen in GL(3, ℤ), die lange Zeit als Kandidaten für exotische S⁴ galten. Akbulut und Kirby zeigten 1979 für den einfachsten Fall (ungekoppelte Rahmung), dass er standard ist. Akbulut bewies 2009–2010 mittels Kirby-Kalkül, dass eine unendliche Teilfamilie dieser Sphären diffeomorph zu S⁴ ist (Annals of Mathematics, 2010). Gompf zeigte kurz darauf, dass eine noch größere Familie standard ist, indem er verborgene Symmetrien der Cappell-Shaneson-Konstruktion ausnutzte (ohne schwere Kirby-Kalkül-Berechnungen). Ein 2024er Ergebnis (Infinite families of standard Cappell-Shaneson spheres, arXiv:2404.05096) dehnte dies auf weitere unendliche Familien aus. Trotzdem: Alle diese Resultate zeigen nur, dass konkrete Kandidaten keine exotischen Sphären sind — sie schließen die Existenz einer exotischen S⁴ nicht aus.

### Warum versagt der Whitney-Trick in Dimension 4?

Der Whitney-Trick reduziert in Dimensionen ≥5 algebraische Schnittpunktzahlen auf geometrische: Bei zwei sich in einem einzigen positiven und einem negativen Punkt schneidenden Untermannigfaltigkeiten der Mitteldimension kann man eine eingebettete Whitney-Scheibe einziehen und die Schnittpunkte eliminieren. Das Problem in Dimension 4: Die Whitney-Scheiben selbst sind 2-dimensional und damit ebenfalls mitteldimensional — sie haben generisch neue Schnitte, die wiederum Whitney-Scheiben benötigen usw. (ein infiniter Regress). Freedman umging dies 1982 topologisch durch Casson-Towers und eine transzendente Grenzwertkonstruktion, die keine Differenzierbarkeit bewahrt. Das glatte Problem bleibt offen.

### Aktueller Stand (2022–2024)

Ein arXiv-Preprint (2209.09968, September 2022, mehrfach revidiert bis Juli 2024) behauptet einen Beweis der SPC4, ist aber noch nicht peer-reviewed verifiziert. 2023 wurde gezeigt: Falls einer von 5 konkreten Knoten scheibenartig (slice) ist, existiert eine exotische S⁴. 2024 wurde der erste Nachweis exotischer eingebetteter S³ in einer glatten geschlossenen 4-Mannigfaltigkeit mit diffeomorphen Komplementen erbracht. Das fundamentale Problem bleibt ungelöst: Es fehlt sowohl ein Beweis der Existenz als auch der Nichtexistenz exotischer S⁴.

---

## Problem 2: Yau-Uniformisierungsvermutung für nicht-kompakte Kähler-Mannigfaltigkeiten

### Problemstellung

Die Yau-Uniformisierungsvermutung (ca. 1974) lautet: Jede vollständige, nicht-kompakte Kähler-Mannigfaltigkeit M mit überall positiver bisektionaler (holomorpher) Krümmung ist biholomorph zu ℂⁿ. Dies ist das Kähler-Analogon des Riemann'schen Uniformisierungssatzes für Riemannflächen und deutlich schwieriger, da die komplexe Geometrie keine einfache Monodromiethorie hat.

### Frühe Teilergebnisse: Shi (1989) und das Volumenwachstum

Wan-Xiong Shi bewies 1989 (Journal of Differential Geometry, Bd. 30): Der Ricci-Fluss (in der Kähler-Variante) existiert auf vollständigen Mannigfaltigkeiten beschränkter Krümmung für kurze Zeit. Dies erlaubt es, die bisektionale Krümmung unter dem Ricci-Fluss zu deformieren. Für maximales Volumenwachstum (Vol(B(x,r)) ~ r^{2n}) und beschränkte Krümmung wurde die Verbindung zu Uniformisierungsaussagen hergestellt.

Das entscheidende technische Hindernis bei allgemeiner positiver (aber möglicherweise unbeschränkter) bisektionaler Krümmung: Shi's Existenzsatz für den Kähler-Ricci-Fluss setzt beschränkte Krümmung voraus. Ohne diese Schranke ist die Kurzzeit-Existenz des Flusses nicht gesichert, und die gesamte Fluss-Methodik bricht zusammen. Ni und Tam (2003) nutzten plurisubharmonische Funktionen und den Wärmfluss auf dem Kotangentialbündel (Journal of Differential Geometry, Bd. 64), um Strukturresultate auch ohne Volumenbeschränkung zu erhalten.

### Chau-Tam und Gang Liu: Entscheidende Fortschritte

Chau und Tam bewiesen (2006, 2011): Eine vollständige nicht-kompakte Kähler-Mannigfaltigkeit mit nichtnegativ beschränkter bisektionaler Krümmung und maximalem Volumenwachstum ist biholomorph zu ℂⁿ. Zum Beweis nutzten sie Li-Yau-Hamilton-Schätzungen und verifizierten ein Lückenprinzip (gap phenomenon): Solche Mannigfaltigkeiten haben automatisch quadratischen Krümmungsabfall im Mittel.

Gang Liu (2016, Journal of Differential Geometry, Bd. 102) bewies: Vollständige Kähler-Mannigfaltigkeiten mit nichtnegativ bisektionaler Krümmung und maximalem Volumenwachstum sind biholomorph zu ℂⁿ. Der Schlüssel: Cheeger-Colding-Theorie der Gromov-Hausdorff-Konvergenz kombiniert mit dem Dreikeis-Theorem für holomorphe Funktionen. Dies bestätigt Yaus Vermutung für den maximal-Volumenwachstum-Fall vollständig.

2019 zeigte Liu (Yau's Uniformization Conjecture, arXiv:1606.08958): Auch im Fall nicht-maximalen Volumenwachstums gibt es Teilergebnisse — unter zusätzlicher Annahme einer unteren Schranke für das Volumenwachstum sind strukturelle Aussagen möglich.

### Soul-Analogon und offene Fronten (2021–2024)

Im Riemannschen Fall liefert der Seele-Satz (Cheeger-Gromoll 1972) für Mannigfaltigkeiten nichtnegativer Schnittkrümmung eine vollständige Strukturbeschreibung. Im Kähler-Fall spielt die Theorie der Lelong-Zahlen und die Busemann-Funktion eine analoge Rolle: Sie misst das „Wachstum nach unendlich" und erlaubt unter günstigen Bedingungen, M als holomorphes Faserbündel über einer kompakten Kähler-Mannigfaltigkeit zu beschreiben. Der April-2024-Preprint (arXiv:2404.08537) „Complete Kähler manifolds with nonnegative Ricci curvature" und der November-2024-Artikel (arXiv:2511.01263) zu minimalen Graden und Volumenwachstum zeigen, dass das Gebiet 2024 hochaktiv ist. Vollständig offen bleibt die Vermutung ohne jede Volumenbeschränkung bei rein positiver bisektionaler Krümmung.

---

## Problem 3: Hartshorne Codimension-2-Vermutung in ℙⁿ (n ≥ 7)

### Problemstellung und historischer Kontext

Robin Hartshorne formulierte 1974 (Bulletin of the AMS) die Vermutung: Jede glatte projektive Untervarietät X ⊂ ℙⁿ der Kodimension c < n/3 ist ein vollständiger Durchschnitt. Im Spezialfall c = 2, n ≥ 7: Glatte Kodimension-2-Untervarietäten von ℙⁿ (n ≥ 7) sind vollständige Durchschnitte. Die Vermutung ist durch den Rang-2-Vektor-Bündel-Ansatz motiviert: X ist vollständiger Durchschnitt ⟺ das zugehörige Bündel (via Serre-Korrespondenz) ist eine direkte Summe von Geradenbündeln.

### Barth-Lefschetz-Theorie und ihre Grenzen

Der Barth-Lefschetz-Satz besagt: Für n ≥ 2c + 1 (d. h. n ≥ 5 für c = 2) hat eine glatte Untervarietät X der Kodimension c dieselben Betti-Zahlen und denselben Fundamentalgruppen wie ℙⁿ bis zur halben Dimension. Für X in ℙ⁷ (Kodimension 2) folgt: H¹(X, ℤ) = 0, Pic(X) ≅ ℤ, erzeugt von 𝒪(1)|_X. Dies sind notwendige Bedingungen für vollständige Durchschnitte, aber bei weitem nicht hinreichend — topologische Kohomologiebedingungen erzwingen keinen algebraischen Erzeugungssatz für das Ideal.

### Das Horrocks-Mumford-Bündel: Warum es für n ≥ 5 nicht hilft

Horrocks und Mumford konstruierten 1973 ein indekomposables Rang-2-Vektorbündel auf ℙ⁴ mit Nullstellenmengen als abelsche Flächen vom Grad 10 (Horrocks-Mumford-Flächen). Dieser Fund zeigte: Auf ℙ⁴ gibt es nicht-spaltende Rang-2-Bündel, also auch Kodimension-2-Untervarietäten (diese abelschen Flächen), die kein vollständiger Durchschnitt sind. Für n ≥ 5 ist kein einziges nicht-spaltendes Rang-2-Bündel bekannt. Hartshorne vermutet, dass für n ≥ 7 jedes Rang-2-Bündel spaltet — was die Vermutung implizieren würde. Klyachko bestätigte die Spalt-Vermutung für torus-äquivariante Bündel auf ℙⁿ für alle n. Ein 2020er Preprint (arXiv:2001.11075) lieferte einen nicht-kombinatorischen Beweis für torische Rang-2-Bündel, aber nicht für allgemeine Bündel.

### Fortschritte und Methoden (2000–2024)

Ionescu und Russo (ca. 2013) bewiesen Hartshornes Vermutung für Varietäten, die durch quadratische Gleichungen definiert werden. Dies nutzt die Theorie der Tangentialkegel und Secantvarietäten. Der arXiv-Preprint 1804.09730 (Eisenbud-Erman, 2018) behandelt die Stärke (strength) von Polynomen und Hartshornes Vermutung in hohem Grad — mit Fortschritten für große n relativ zum Grad. Für n ≥ 7 und allgemeine Kodimension-2-Untervarietäten ist die Frage 2024 vollständig offen. Ein Mai-2024-Preprint (arXiv:2405.11933) untersucht glatte Untervarietäten kleinen Grades und kleiner Kodimension, bleibt aber unter dem allgemeinen Fall.

### Ran-Resultat und aktuelle Lage

Ran (1990) bewies Kohomologieeinschränkungen für deformationstheoretische Ansätze, die zeigen, dass für „allgemeine" Einbettungen vollständige Durchschnitte erwartet werden. Ellia, Fiorentini und Miró-Roig haben systematisch Familien von Beispielen in kleinen Dimensionen untersucht und Moduli-Räume von Kodimension-2-Untervarietäten katalogisiert. Das Grundproblem bleibt: Es fehlt ein Satz, der aus den topologischen und kohomologischen Eigenschaften algebraische Erzeugungsaussagen ableitet.

---

## Problem 4: Grothendieck Standard-Vermutungen über algebraische Zyklen

### Die vier Vermutungen

Grothendieck formulierte um 1968 vier miteinander verflochtene Vermutungen über algebraische Zyklen auf glatten projektiven Varietäten über einem algebraisch abgeschlossenen Körper:

- **Vermutung A (Lefschetz-Typ):** Die Lefschetz-Involution Λ (Operator der harten Lefschetz-Zergliederung) wird durch einen algebraischen Zyklus induziert.
- **Vermutung B (algebraische Äquivalenz):** Äquivalente Bedingung zu A: Der Lambdaoperator ist algebraisch.
- **Vermutung C (Künneth-Zerlegung):** Die Künneth-Projektatoren auf H^{2k}(X × X) sind algebraisch.
- **Vermutung D (numerische vs. homologische Äquivalenz):** Numerische und homologische Äquivalenz von Zyklen fallen zusammen.

Vermutung A und B sind äquivalent. Vermutung D impliziert A/B. Zusammen implizierten sie den härtesten Teil der Weil-Vermutungen (das Riemann-Hypothesen-Analogon für Varietäten über endlichen Körpern).

### Warum beweist Deligne 1974 nicht Vermutung B?

Deligne bewies 1974 die Weil-Vermutungen (insbesondere Weil II: die Reinheitsaussage der Eigenwerte des Frobeniusautomorphismus) durch eine geniale Methode mit Leray-Spektralsequenzen und einer „tordue" (verdrehten) Version der Argumentation — ohne auf Grothendiecks Standard-Vermutungen zurückzugreifen. Deligne umging das Problem also, statt es zu lösen. Die Standard-Vermutungen blieben unbewiesen. Insbesondere Vermutung B ist für allgemeine Varietäten offen, weil sie eine algebraische Konstruktion des Λ-Operators erfordert — also einen expliziten Zyklus in X × X, der als Korrespondenz die Lefschetz-Involution realisiert.

### Bekannte Fälle und Fortschritte

Lieberman (1968) und Kleiman (1968) bewiesen Vermutung A/B für abelsche Varietäten: Dort ist Λ durch die Poincaré-Linientorsion algebraisch definiert. Moonen und Zarhin (1999, 2000er) untersuchten, für welche Klassen abelscher Varietäten (nach ihrer Endomorphismenalgebra klassifiziert) alle Standard-Vermutungen gelten — ein reiches Programm, das auch mit dem Mumford-Tate-Gruppenansatz zusammenhängt. Katz und Messing (1974) bewiesen die Künneth-Vermutung für Varietäten über endlichen Körpern via Weil-Vermutungen.

### Voevodsky und motivische Kohomologie

Voevodsky (1990er–2000er) konstruierte eine triangulierte Kategorie der Motive DM(k) und definierte motivische Kohomologie als Ext-Gruppen darin. Dies liefert ein konzeptuell befriedigendes Rahmenwerk, lässt aber die Standard-Vermutungen offen: Die Frage, ob DM(k) eine abelsche Kategorie (die Kategorie der reinen Motive) als Herz enthält, ist äquivalent zu den Standard-Vermutungen. Der Beweis der Bloch-Kato-Vermutung (Voevodsky 2010) betrifft Torsionsphenomene und liefert keinen direkten Fortschritt für die Standard-Vermutungen.

### Scholze und neueste Entwicklungen (2020–2024)

Scholzes Programm der Geometrisierung des lokalen Langlands-Programms (mit Fargues, 2021, publiziert in Annals 2024) nutzt Diamanten und v-Stacks — neue geometrische Objekte jenseits klassischer Schemata. Ein 2025er Preprint Scholzes über „Motivische Geometrisierung des lokalen Langlands" (in Crelle erscheinend) und sein 2024er Preprint über Berkovich-Motive zeigen, dass motivische Methoden in der Langlands-Geometrie zunehmend relevant werden. Direkte Auswirkungen auf Grothendiecks Standard-Vermutungen sind bisher nicht erzielt, aber die motivische Homotopietheorie nähert sich dem Problem von einer neuen Seite. Stand 2026: Alle vier Standard-Vermutungen offen für allgemeine Varietäten.

---

## Problem 5: Ist das Kontsevich-Integral eine vollständige Knoteninvariante?

### Das Problem

Maxim Kontsevich definierte 1993 ein Knotenintegral Z(K) als formale Potenzreihe in Chord-Diagrammen (Elemente des Vektorraums $\mathcal{A}$ der abgeschlossenen Diagramme). Z(K) ist eine universelle Vassiliev-Invariante über ℚ: Jede rationale Vassiliev-Invariante endlichen Typs faktorisiert über Z(K) durch einen Gewichtssystemhomomorphismus. Die offene Frage: Falls Z(K₁) = Z(K₂), folgt daraus K₁ ≅ K₂ (isotop)?

### Vassiliev-Invarianten und Trennung

Es ist bekannt, dass die Gesamtheit aller rationalen Vassiliev-Invarianten (äquivalent: Z selbst) Knoten von Unknoten zu unterscheiden vermag (für fast alle bekannten Knoten), aber es ist nicht bewiesen, dass Vassiliev-Invarianten alle Knotenpaare trennen. Das Trennungsproblem ist vollständig offen. Bar-Natan, Garoufalidis, Rozansky und Thurston haben den Aarhus-Integralkalkül entwickelt (2002), der Z mit dem LMO-Invarianten für 3-Mannigfaltigkeiten verbindet und Z eine universelle Rolle gibt.

### Mutanten und das Kontsevich-Integral

Ein Mutant von K ist ein Knoten, der durch Verdrehung eines 2-Strangtangels um 180° entsteht. Jones-Polynome, Alexander-Polynome und HOMFLY-Polynome können Mutanten nicht unterscheiden. Rozansky (2010) und andere untersuchten das Verhalten von Z unter Mutation: Das Kontsevich-Integral unterscheidet nach aktuellem Stand Mutanten ebenfalls nicht (oder zumindest nicht durch die bekannten Projektionsmethoden), obwohl dies nicht rigoros bewiesen ist. Dies ist ein starkes Indiz, dass Z kein vollständiger Invariant ist.

### Habiros Zyklotomische Entwicklung

Habiro konstruierte eine zyklotomische Entwicklung des gefärbten Jones-Polynoms: Für jeden Knoten K und jede Farbe n existieren eindeutige Laurent-Polynome $H_{K,k}(q) \in \mathbb{Z}[q^{\pm 1}]$ derart, dass $J_n(K) = \sum_{k=0}^{n-1} H_{K,k}(q) \binom{n}{k+1}_{q^2}$. Diese Entwicklung konvergiert in der zyklotomischen Vervollständigung und erzeugt einen einheitlichen Invarianten für ganzzahlige Homologie-3-Sphären. Die Beliakova-Gorsky-Konstruktion (2024, arXiv) erweitert dies auf $\mathfrak{gl}_N$-Invarianten via Interpolations-Macdonald-Polynome.

### Le-Murakami-Ohtsuki-Invariante

Der LMO-Invariant (Le, Murakami, Ohtsuki, 1998) ist der universelle endliche Typ-Invariant für ganzzahlige Homologie-3-Sphären, konstruiert als Grenzwert des Kontsevich-Integrals von Knoten in der Poincaré-Chirurgie-Präsentation. Eine 2025er Arbeit (Journal of Topology, Audoux) konstruiert einen universellen endlichen Typ-Invarianten für Knoten in Homologie-3-Sphären, was einen neuen Kandidaten für vollständige Invarianten liefert. Ob Z oder LMO vollständige Invarianten sind, ist fundamental offen.

---

## Problem 6: Systolische Geometrie — Optimale Konstanten für T^n und ℝP^n (n ≥ 3)

### Hintergrund und Gromovs Fundament

Die Systole sys(M, g) einer Riemannschen Mannigfaltigkeit (M, g) ist die Länge der kürzesten nicht-zusammenziehbaren geschlossenen Kurve. Loewner bewies 1949 für den 2-Torus T²: sys(T²)² ≤ (2/√3) · Vol(T², g), mit Gleichheit beim regelmäßig-sechseckigen Gitterflächen. Pu (1952) bewies das Analogon für ℝP²: sys(ℝP²)² ≤ (π/2) · Vol(ℝP², g). Gromov (1983) verallgemeinerte diese Ungleichungen auf alle wesentlichen Mannigfaltigkeiten (essential manifolds): sys(M)^n ≤ C_n · Vol(M, g), wobei die Konstante C_n aus dem Füllungsradius-Konzept hervorgeht.

### Warum T³ so viel schwerer als T² ist

Für T² ist die scharfe Konstante bekannt und durch die hexagonale Einbettung realisiert. Für T³ und höhere Tori n ≥ 3 ist die scharfe Konstante in Gromovs Ungleichung unbekannt. Der Grund liegt in der drastisch wachsenden Komplexität der Optimierungsaufgabe: Für T² ist der Modulraum der Metriken (Gitter-Modulraum) 1-dimensional, das Optimierungsproblem explizit lösbar. Für T³ ist der Gitter-Modulraum GL(3, ℤ) \ GL(3, ℝ) / O(3) — ein hochdimensionaler Raum ohne kanonisches Optimierungsprinzip. Sabourau bewies (ca. 2010) systemische Ungleichungen mit verbesserten (aber nicht scharfen) Konstanten für T³ durch das Studium der Längenspektren und Morse-Theorie auf dem Schlingenraum.

### Guths Minimax-Methode (2011)

Larry Guth bewies 2011 (Geometric and Functional Analysis) die Gromovsche systolische Ungleichung für T^n über minimale Hyperflächen — basierend auf dem Schoen-Yau-Beweis, dass T^n keine Metrik positiver Skalarkrümmung trägt. Guths kurzer Beweis liefert eine explizite Konstante, aber nicht die scharfe. Entscheidend ist Guths Füllungsradius-Abschätzung: FillRad(M) ≥ c · sys(M), deren Schärfe ebenfalls offen ist. Sein Minimax-Ansatz verbindet Variationsrechnung auf Zyklenräumen mit der Geometrie von Sweepouts (Flächen-Überstrichen) und wurde von Nabutovsky, Rotman und Sabourau weiter ausgebaut.

### Nabutovsky-Rotman-Sabourau: Sweepouts und Längenspektren

Nabutovsky und Rotman bewiesen (Journal of Differential Geometry, 2011): Auf jedem Riemannschen 2-Sphären-Paar existieren k geodätische Schleifen der Länge ≤ 20k · diam(M). Nabutovsky, Rotman und Sabourau (Geometric and Functional Analysis, 2021) untersuchten Sweepouts abgeschlossener Riemannscher Mannigfaltigkeiten und bewiesen min-max-Abschätzungen der Form W_k(M) ≤ c_n · Vol(M)^{1/n} für Breiten (widths). Ein Juli-2024er Preprint (arXiv:2407.12673) von Nabutovsky-Rotman über kurze einfache geodätische Schleifen auf S² zeigt die anhaltende Aktivität.

### Optimale Konstanten für ℝP^n und der heutige Stand

Für ℝP² ist die Pu-Konstante scharf (realisiert durch die runde Metrik). Für ℝP^n (n ≥ 3) ist die scharfe systolische Konstante offen. Katz zeigte (Systolic Geometry and Topology, AMS Monographs, 2007): Die stabilen Systolen unterliegen schärferen Bounds, aber die 1-dimensionalen Systolen in höheren Dimensionen bleiben unkontrolliert ohne Volumenwachstums-Annahmen. Eine 2023er Arbeit (Journal of Geometry, Springer) „Extending Gromov's optimal systolic inequality" (arXiv:2407.03803) erweitert Gromovs Ungleichung auf neue Klassen von Mannigfaltigkeiten, gibt aber keine scharfen Konstanten für T^n (n ≥ 3) oder ℝP^n (n ≥ 3). Das Feld ist 2024 hochaktiv, die fundamentale Frage nach scharfen Konstanten vollständig offen.

---

## Zusammenfassung und Querverbindungen

| Problem | Status | Haupthindernis | Aktivste Methode (2024) |
|---------|--------|----------------|------------------------|
| Glattes Poincaré S⁴ | Offen | Keine Invariante für b₂=0 | Kirby-Kalkül, Slice-Knoten |
| Yau-Uniformisierung | Teilweise gelöst (max. vol.) | Ohne Volumenschranke offen | Kähler-Ricci-Fluss, Gromov-Hausdorff |
| Hartshorne Kodim.-2 | Offen (n≥7) | Topologie ≠ Algebra | Vektorbündelmethoden, Strength |
| Standard-Vermutungen | Alle offen (allg. Fall) | Algebraischer Λ-Operator | Motivische Homotopie, Perfectoids |
| Kontsevich vollständig? | Offen | Mutanten, Integer-Koeff. | Zyklotomische Entwicklung, LMO |
| Systolische Konstanten | Offen (n≥3) | Hochdim. Modulraum | Sweepouts, Minimax-Theorie |

**Gemeinsames Thema:** In allen Problemen treffen auf ein- und derselben Bühne die Grenzen von Kohomologie-/Invariantentheorie auf das Versagen direkter Konstruktionsmethoden in mittlerer Dimension oder höheren Dimensionen.

---

## Quellen und Referenzen

1. Akbulut, S. (2009/2010). *Cappell-Shaneson homotopy spheres are standard.* Annals of Mathematics 171(3).
2. Akbulut, S. & Greene, R. (2010). *More Cappell-Shaneson spheres are standard.* AGT 10(3).
3. arXiv:2404.05096 — Infinite families of standard Cappell-Shaneson spheres (2024).
4. arXiv:2209.09968 — On 4-dimensional smooth Poincaré conjecture (2022–2024).
5. Freedman, M. H. (1982). *The topology of four-dimensional manifolds.* J. Differential Geometry 17(3).
6. Chau, A. & Tam, L.-F. (2006, 2011). *On the complex structure of Kähler manifolds with nonneg. curvature.* J. Differential Geometry.
7. Liu, G. (2016). *On the volume growth of Kähler manifolds with nonneg. bisectional curvature.* JDG 102(3), 485–500.
8. Liu, G. (2016). *Gromov-Hausdorff limits of Kähler manifolds and the finite generation conjecture.* Annals of Mathematics 184(3).
9. Ni, L. & Tam, L.-F. (2003). *Plurisubharmonic functions and complete Kähler manifolds.* JDG 64.
10. Hartshorne, R. (1974). *Varieties of small codimension in projective space.* Bull. AMS.
11. Horrocks, G. & Mumford, D. (1973). *A rank 2 vector bundle on P4 with 15,000 symmetries.* Topology.
12. arXiv:1804.09730 — Eisenbud-Erman: Strength and Hartshorne's Conjecture (2018).
13. arXiv:2405.11933 — Smooth projective varieties of small degree and codimension (2024).
14. Grothendieck, A. (1968). *Standard conjectures on algebraic cycles.* Algebraic Geometry, Bombay Colloquium.
15. Deligne, P. (1974). *La conjecture de Weil I.* Publ. Math. IHES 43.
16. Milne, J. S. (2002). *Polarizations and Grothendieck's standard conjectures.* arXiv:math/0103175.
17. Scholze, P. & Fargues, L. (2021/2024). *Geometrization of the local Langlands correspondence.* Annals.
18. Kontsevich, M. (1993). *Vassiliev's knot invariants.* Adv. Soviet Mathematics 16(2).
19. Le, T. Q. T., Murakami, J. & Ohtsuki, T. (1998). *On a universal perturbative invariant.* Topology.
20. Habiro, K. (2000). *On the quantum sl₂ invariants and cyclotomic expansions.*
21. Beliakova, A. & Gorsky, E. (2024). *Cyclotomic expansions for gl_N knot invariants.* arXiv.
22. Gromov, M. (1983). *Filling Riemannian manifolds.* J. Differential Geometry 18.
23. Guth, L. (2010). *Systolic inequalities and minimal hypersurfaces.* GAFA 20(2).
24. Nabutovsky, Rotman & Sabourau (2021). *Sweepouts of closed Riemannian manifolds.* GAFA.
25. arXiv:2407.03803 — Extending Gromov's optimal systolic inequality (2024).
