# Gruppe B: Vier offene topologische Probleme — Umfassende Forschungsübersicht

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Letzte Aktualisierung**: 2026-03-12
**Quellen**: 16 Webrecherchen, Arxiv-Preprints, Annals of Mathematics, IHÉS, Quanta Magazine

---

## Inhaltsverzeichnis

1. [Mandelbrot MLC — Lokale Zusammenhängigkeit](#1-mandelbrot-mlc--lokale-zusammenhängigkeit)
2. [Conway-Knoten: Glatte Scheibigkeit und Konkordanz-Invarianten](#2-conway-knoten-glatte-scheibigkeit-und-konkordanz-invarianten)
3. [Whitehead Asphärizitätsvermutung](#3-whitehead-asphärizitätsvermutung)
4. [Jones-Unknoten-Vermutung und Volumen-Vermutung](#4-jones-unknoten-vermutung-und-volumen-vermutung)

---

## 1. Mandelbrot MLC — Lokale Zusammenhängigkeit

### 1.1 Hintergrund

Die **MLC-Vermutung** (Mandelbrot Locally Connected) besagt, dass die Mandelbrot-Menge
M ⊂ ℂ an jedem Punkt lokal zusammenhängend ist. Obwohl M selbst zusammenhängend ist
(Douady–Hubbard, 1982, über Koebe-Verzerrungssätze und den Beweis der
Böttcher-Konjugation), ist lokale Zusammenhängigkeit eine stärkere, topologische
Aussage: Für jeden Punkt c ∈ M und jede Umgebung U ∋ c soll es eine
zusammenhängende offene Menge V mit c ∈ V ⊂ U ∩ M geben.

**Warum ist MLC so bedeutsam?** Die Vermutung ist äquivalent zur
**Dichte-der-Hyperbolizität-Vermutung**: Die Menge der hyperbolischen Parameter
(d. h. Parameter c, für die f_c(z) = z² + c einen endlich-anziehenden Periodenpunkt
hat) liegt dicht in M. Dies ist wiederum äquivalent zur Aussage, dass M durch seine
externen Winkel vollständig durch einen wohldefinierten Parameterraum-Kodierungsraum
parametrisiert werden kann — was die Mandelbrot-Menge kombinatorisch verstehbar machen
würde. Der **Feigenbaum-Punkt** c∞ ≈ −1.40115 ist der klassischste infinit
renormalisierbare Parameter: Hier verdoppelt sich die Periode unbegrenzt (2, 4, 8,
16, ...), und die Mandelbrot-Menge zeigt in seiner Umgebung selbstähnliche
Schachtelung ins Unendliche.

### 1.2 Gescheiterte Ansätze

**Yoccoz-Puzzle (1990–1992)**: Jean-Christophe Yoccoz entwickelte ein ausgeklügeltes
kombinatorisches Werkzeug (das "Yoccoz-Puzzle"), das durch successive Verfeinerung
von Fundamentaldomänen arbeitet. Er bewies MLC für alle **endlich renormalisierbaren**
Parameter (Acta Mathematica 1994). Das Schlüsselargument: Bei endlicher
Renormalisierbarkeit konvergieren die Puzzle-Stücke exponentiell schnell gegen den
Punkt, was Rigidität erzwingt.

**Warum scheitert Yoccoz-Puzzle bei infinit renormalisierbaren Parametern?**
Der Beweis bricht zusammen, weil die Puzzle-Stücke bei infinit renormalisierbaren
f_c niemals auf eine endliche Renormalisierungstiefe beschränkt sind — es gibt
unendlich viele verschachtelte "Baby-Mandelbrot-Kopien" in jeder Umgebung von c.
Konkret: Die diameters der Puzzle-Stücke konvergieren nicht gegen Null, es fehlen
sogenannte *a priori bounds* — quantitative Geometrieschranken, die kontrollieren, wie
stark sich die Renormalisierungen geometrisch verzerren. Ohne diese Schranken kann
die Rigidität (Starrheit der Dynamik) nicht geschlossen werden.

**Siegel-Schreiben und kleine Nenner**: Für Parameter nahe Siegel-Scheiben
(irrationale Rotationszahlen) scheitern Yoccoz-Methoden an Kleinnennerproblemen
(Liouville-Approximation). Buff und Chéritat zeigten 2012 in den Annals of
Mathematics (Vol. 176, S. 673–746), dass quadratische Polynome mit Cremer-Fixpunkten,
Siegel-Scheiben oder unendlich vielen Satellitenrenormalisierungen **Julia-Mengen mit
positivem Lebesgue-Maß** besitzen können. Dies war zunächst als potentielles
Gegenbeispiel zu MLC untersucht worden, widersprach MLC aber nicht direkt.

### 1.3 Aktuelle Ergebnisse

**Kahn–Lyubich (2006–2009)**: Lokale Zusammenhängigkeit für bestimmte infinit
renormalisierbare Parameter der "bounded satellite type" — dies sind Parameter mit
beschränkter kombinatorischer Komplexität.

**Avila–Kahn–Lyubich–Shen**: Kombinatorische Rigidität für unikritische Polynome
(Annals of Mathematics, 2009), d. h. für f_c(z) = z^d + c mit d ≥ 2. Dies liefert
MLC für unendlich viele neue Parameter, aber nicht für die klassischen Feigenbaum-Fälle.

**Dudko–Lyubich (2023/2025)**: Das Paper **"MLC at Feigenbaum points"** (arXiv:2309.02107,
v1 September 2023, v3 Dezember 2025) ist der bislang bedeutendste Durchbruch. Die
Autoren beweisen **a priori bounds für Feigenbaum-quadratische Polynome** (infinit
renormalisierbare Polynome f_c von beschränktem Typ). Kernmethode: **parabolische
Renormalisierung**, eine neue analytische Technik, die parabole Implosion kontrolliert
und quantitative Schranken für die Renormalisierungsoperatoren liefert. Daraus folgen:
(i) lokale Zusammenhängigkeit der entsprechenden Julia-Mengen J(f_c),
(ii) MLC am entsprechenden Parameter c,
(iii) **Skalenuniversalität** (Goldberg-Khanin-Sinai-Vermutung für komplexe
Dreifachrenormalisierungen, offen seit ~1983).

Das Paper löst das MLC-Problem für den klassischen Feigenbaum-Parameter sowie für
komplexe Tripling-Renormalisierungen — zwei der hartnäckigsten Spezialfälle.

**Kapiamba (2023)**: Neue Ergebnisse zur lokalen Zusammenhängigkeit bei Satelliten-
Parametern von beschränktem Typ (GAFA 2023), ergänzend zu Dudko–Lyubich.

### 1.4 Offene Teilprobleme

- **Allgemeiner Fall infiniter Renormalisierbarkeit**: MLC für beliebige (unbeschränkte)
  unendlich renormalisierbare Parameter bleibt offen. Die parabolische Renormalisierung
  deckt nur *beschränkte* Typen ab.
- **Dichte der Hyperbolizität (allgemein)**: Obwohl für reelle quadratische Polynome
  von Lyubich sowie Graczyk–Świątek (1997) bewiesen, fehlt der Beweis für die
  komplexe Familie.
- **Positive-Maß-Fragen**: Ist das Lebesgue-Maß von M gleich Null? Bekannt nur für
  die reelle Achse (Lyubich 2002).
- **Kubische und höhergradige Polynome**: MLC-Analoga für mehrdimensionale
  Parameterräume sind weitgehend unerforscht.
- **Parabolische Implosion als generelles Werkzeug**: Ob Dudko–Lyubichs Technik auf
  nicht-Feigenbaum-Typen verallgemeinert werden kann, ist unklar.

---

## 2. Conway-Knoten: Glatte Scheibigkeit und Konkordanz-Invarianten

### 2.1 Hintergrund

Der **Conway-Knoten** (11n34) wurde 1970 von John Horton Conway vorgeschlagen.
Er ist der Mutant des Kinoshita-Terasaka-Knotens (KT): Man wählt ein 2-Tangle im
Diagramm, spiegelt es, und erhält den anderen Knoten. Beide Knoten haben dasselbe
Jones-Polynom, dasselbe Alexander-Polynom, und dasselbe Conway-Polynom (nämlich das
des Unknotens).

**Topologische Scheibigkeit**: Ein Knoten K ⊂ S³ heißt *topologisch scheibartig*,
wenn er der Rand einer lokal flachen Scheibe D² ↪ B⁴ ist. Freedman's Theorem (1982)
sagt: Jeder Knoten mit Alexanderpolynom 1 ist topologisch scheibartig. Da Conway-Knoten
und KT-Knoten beide Alexander-Polynom 1 haben, sind beide **topologisch scheibartig**.

**Glatte Scheibigkeit**: Ein Knoten heißt *glatt scheibartig*, wenn eine glatte
Einbettung D² ↪ B⁴ existiert. Der KT-Knoten ist glatt scheibartig (konstruiert über
eine explizite Bändersumme). Beim Conway-Knoten war die glatte Scheibigkeit 50 Jahre
lang offen — das letzte offene Klassifizierungsproblem für Primknoten mit ≤ 12 Kreuzungen.

### 2.2 Gescheiterte Ansätze

**Klassische Konkordanz-Invarianten versagen**:

- **Alexander-Polynom**: Δ_Conway(t) = 1 — scheibartig-neutral, keine Information.
- **Arf-Invariant**: Arf(Conway) = 0 — topologisch kompatibel, keine Obstruktion
  gegen glatte Scheibigkeit.
- **Levine-Tristram-Signaturen** σ_K(ω): Die Signatur-Funktion weist dem Knoten für
  jedes ω ∈ S¹ eine ganzzahlige Invariante zu. Für Mutanten gilt:
  σ_Conway(ω) = σ_KT(ω) für alle ω, da die Mutationsoperation die Levine-Tristram-
  Signaturenfunktion erhält. Damit sind Signaturen komplett blind für den Unterschied.
- **Rasmussen s-Invariant** (direkt): s ist Mutationsinvariant — s(Conway) = s(KT) = 0.
  Da KT scheibartig ist, folgt s = 0 für beide, was keine Obstruktion liefert.
- **τ-Invariante (Ozsváth-Szabó, Heegaard-Floer)**: τ(Conway) = 0, ebenfalls
  Mutationsinvariant.
- **Ozsváth–Stipsicz–Szabó Υ-Invariante**: Υ_K(t) ist eine stückweise lineare Funktion
  auf [0, 2] und Konkordanzinvariant. Auch Υ(Conway)(t) = 0 für alle t — keine Hilfe.

Alle bekannten klassischen Invarianten scheitern, weil Mutation (im Sinne von
Conway-Mutation) sie invariant lässt.

### 2.3 Aktuelle Ergebnisse — Piccirillo 2020

**Lisa Piccirillo** (Universität Texas Austin, Dissertation 2019; Annals of Mathematics
2020) löste das Problem durch eine völlig neue Strategie:

**Kern: Spur-Äquivalenz (Knot Trace Technique)**.
Die **Spur** eines Knotens K ist die 4-Mannigfaltigkeit X(K), die durch Anheften eines
0-gerahmten 2-Henkels an B⁴ entlang K entsteht. Schlüsselobservation:
*K ist glatt scheibartig genau dann, wenn X(K) glatt in S⁴ einbettbar ist.*

Piccirillo konstruierte einen anderen Knoten K' mit X(K') ≅ X(Conway) (diffeomorph!),
aber K' ≠ Conway. Auf K' konnte sie dann den **Rasmussen s-Invariant** anwenden — und
dieser lieferte s(K') ≠ 0, also K' nicht glatt scheibartig. Da X(Conway) = X(K'), folgt
auch der Conway-Knoten ist **nicht glatt scheibartig**.

Die Konstruktion nutzte **Cork-Twists** und Kirby-Kalkül: Durch gezielte Manipulationen
an 4-Mannigfaltigkeits-Diagrammen (Handles, Blowups) erzeugte sie K' mit gleicher Spur.

**Einschränkung**: Diese Methode beweist nur die negative Aussage (nicht scheibartig).
Sie sagt nichts darüber aus, *wie weit* der Conway-Knoten von Scheibigkeit entfernt ist.

**Kinoshita-Terasaka Vergleich**: KT ist der einzige Conway-Mutant mit bekannten
Scheibartigkeitseigenschaften. Seit Piccirillo sind alle Primknoten mit ≤ 12 Kreuzungen
klassifiziert.

### 2.4 Offene Teilprobleme

- **4-Genus des Conway-Knotens**: g₄(Conway) = ? Bekannt: g₄ ≥ 1 (nicht scheibartig),
  aber ist g₄ = 1 oder größer?
- **Invarianten, die Mutation unterscheiden**: Welche schwächeren Konkordanz-Invarianten
  (z.B. höhere ν_n-Invarianten, Z₄-equivariante Floer-K-Theorie) differenzieren
  Conway-Mutanten?
- **Conway-Knoten in der Konkordanzgruppe**: Welche Ordnung hat [Conway] in C (der
  topologischen Konkordanzgruppe)? Ist er endliche Ordnung oder unendliche Ordnung?
- **Gordon-Litherland-Formel**: Läßt sich die Signatur-Invariante durch Spanning-
  Surfaces-Argumente (wie bei Gordon–Litherland 1978) auf den conway-Knoten ausbauen?
- **Verallgemeinerung auf Satellitenknoten**: Sind Satelliten des Conway-Knotens scheibartig?
- **Ozsváth–Stipsicz–Szabó κ-Invariante**: Diese von 2017 eingeführte
  Invariante für Legendre-Knoten könnte bei Mutanten differenzieren — unerforscht für Conway.

---

## 3. Whitehead Asphärizitätsvermutung

### 3.1 Hintergrund

**J. H. C. Whitehead** formulierte 1941 folgende Frage: *"Ist jeder zusammenhängende
Teilkomplex eines asphärischen 2-dimensionalen CW-Komplexes selbst asphärisch?"*
Ein CW-Komplex X heißt asphärisch, wenn πₙ(X) = 0 für alle n ≥ 2 — d. h. alle
höheren Homotopiegruppen verschwinden, der Komplex wird vollständig durch π₁(X)
beschrieben.

**Äquivalente Formulierung**: Ist jeder 2-dimensionale CW-Komplex, der homotopieäquivalent
zu einem K(G,1)-Raum ist, ein K(G,1)-Raum, wenn man eine 2-Zelle entfernt?

**D(2)-Problem (Wall 1965, Johnson 2003)**: Das D(2)-Problem fragt, ob jeder endliche
CW-Komplex X mit H_n(X; A) = 0 für alle n ≥ 3 und alle Koeffizienten A homotopieäquivalent
zu einem 2-dimensionalen Komplex ist. Johnson zeigte 2003 in seinem Buch *"Stable Modules
and the D(2)-Problem"* (Cambridge LMS Lecture Notes 301), dass D(2) äquivalent zu
Whiteheads Frage ist — genauer: das D(2)-Problem für eine endliche Gruppe G ist äquivalent
zur Aussage, dass gewisse algebraische 2-Zyklen über ℤ[G] geometrisch realisierbar sind.

**Lyndon-Satz (1950)**: Einrelatorgruppen sind asphärisch. Das liefert die Vermutung
für alle 2-Komplexe mit Fundamentalgruppe = Einrelatorgruppe.

**Cockroft-Eigenschaft**: Ein 2-Komplex hat die Cockroft-Eigenschaft, wenn jede
2-Zelle in der Fundamentalklasse von H₂ trivialen Beitrag leistet. Cockroft-Komplexe
sind ein hinreichendes Kriterium für Asphärizität.

### 3.2 Gescheiterte Ansätze

**CAT(0)-Methoden**: Man hoffte, dass wenn π₁(X) eine CAT(0)-Gruppe ist (also auf
einem CAT(0)-Raum negativ-gekrümmter Geometrie wirkt), die Asphärizität folgt.
CAT(0)-Gruppen sind kontraktierbar, was scheinbar Asphärizität liefert. Das Problem:
CAT(0)-Räume sind im Allgemeinen nicht 2-dimensional. Der überlagernde Universalraum
Ẽ eines 2-Komplexes muss kein CAT(0)-Raum sein, selbst wenn π₁ eine CAT(0)-Gruppe ist.
Das geometrische Realisierungsproblem bleibt ungeklärt.

**Wise's cubulierte Gruppen (2012)**: Dani Wise bewies 2012, dass viele Gruppen aus
der geometrischen Gruppentheorie *cubulierbar* sind (d. h. auf CAT(0)-Kubelkomplexen
wirken). Jedoch liefert Cubulierbarkeit keine direkte Information über π₂ eines
2-Komplexes. Die Universaläberlagerung eines 2-Komplexes ist kein CAT(0)-Kubel-Komplex
im Allgemeinen. Agol–Wise (2013) lieferten Resultate über Untergruppen, aber nicht
über Asphärizität von 2-Komplexen.

**Magnus-Problem (One-Relator Groups)**: Das Magnus-Problem fragt, ob alle
Untergruppen von Einrelatorgruppen selbst Einrelatorgruppen sind (oder zumindest
"good" für Asphärizität). Gersten zeigte in den 1980ern Asphärizität für bestimmte
Klassen von Präsentationen durch diagrammatische Reduzierbarkeit. Collins und Miller
(1998) zeigten, dass das Wortproblem für Einrelatorgruppen mit Torsion in bestimmten
Fällen lösbar ist, ohne dass dies die Whitehead-Vermutung erhellt.

**LOT-Präsentationen (Labeled Oriented Trees)**: LOTs modellieren Bänderkomplement-
Gruppen und sind ein wichtiger Spezialfall. Howie und Gersten zeigten Asphärizität
für große Klassen von LOTs (1983–1991). Aber der allgemeine LOT-Fall ist offen, und
damit auch die Whitehead-Vermutung für Bänderkomplementgruppen.

**Andrews-Curtis-Vermutung**: Falls AC (keine Vereinfachung von Präsentationen nötig)
und die Bänderscheiben-Vermutung beide gelten, gibt es keinen Typ-(a)-Gegenbeispiel
zu Whitehead (Howie 1984). Aber AC selbst ist offen — und der AC-Gegenbeispielen-Zoo
ist zahlreich.

### 3.3 Aktuelle Ergebnisse

- **Whitehead-Vermutung für quasi-konstruierbare Komplexe** ist bewiesen (Howie 1983).
- **Für alle Komplexe mit Fundamentalgruppe einer torsionsfreien One-Relator-Gruppe**:
  asphärisch (aus dem Magnus–Freiheitssatz).
- **Periodische Gruppen**: Die Whitehead-Asphärizitätsvermutung gilt für alle
  endlichen Gruppen G, d.h. kein endlicher K(G,1)-2-Komplex hat nicht-asphärischen
  Teilkomplex (trivial, da endliche Gruppen kein K(G,1) mit endlichem 2-Komplex besitzen).
- **Arxiv 2021**: Ein Preprint (arXiv:2107.12293) behauptete eine Antwort auf die
  Whitehead-Frage — der Paper wurde jedoch nicht in begutachteten Journalen publiziert
  und gilt als nicht verifiziert.
- **Wall's D(2) für spezifische Gruppen**: Für dihedral groups Dₙ und zyklische Gruppen
  ℤₙ ist D(2) positiv beantwortet (Johnson 2003 und nachfolgende Arbeiten).

### 3.4 Offene Teilprobleme

- **Allgemeiner Fall**: Die Vermutung für beliebige endliche asphärische 2-Komplexe
  bleibt vollständig offen.
- **D(2)-Problem für nicht-abelsche Gruppen der Ordnung > 1**: Für viele einfache
  Gruppen ist D(2) ungelöst.
- **Geometrische Realisierung algebraischer 2-Komplexe**: Gibt es eine algebraische
  Obstruction (in K₀(ℤ[G]))die Walls Endlichkeitshindernis liefert?
- **Verhältnis zu Andrews-Curtis**: Falls AC positiv gelöst würde, würde dies
  gewisse Typ-(a)-Gegenbeispiele ausschließen — aber AC selbst ist offen.
- **Bridson–de la Harpe**: Die Frage, ob CAT(0)-Gruppen immer asphärische
  Klassifizierungsräume haben, ist ein verwandtes offenes Problem.
- **Magnus-Problem (allgemein)**: Sind alle Untergruppen von Einrelatorgruppen selbst
  homologisch trivial in Grad ≥ 2?

---

## 4. Jones-Unknoten-Vermutung und Volumen-Vermutung

### 4.1 Hintergrund

**Jones-Polynom**: Vaughan Jones entdeckte 1984 ein Laurent-Polynom V_L(t) ∈ ℤ[t^{1/2},
t^{-1/2}] als topologische Invariante orientierter Knoten und Links. Das Jones-Polynom
des Unknotens ist V_∅(t) = 1.

**Jones-Unknoten-Vermutung**: Falls V_K(t) = 1 für einen Knoten K, ist K der Unknoten?
Kein Gegenbeispiel ist bekannt. Computationell verifiziert für alle Primknoten bis
24 Kreuzungen (Tuzun–Sikora 2020/2021, arXiv:2003.06724: über 300 Millionen Knoten
überprüft). Für *Links* (Verschlingungen) gibt es nicht-triviale Beispiele mit Jones-
Polynom eines Unlinks (Thistlethwaite 2001) — aber kein Knotengegenbeispiel.

**Volumen-Vermutung**: Kashaev (1997) vermutete für hyperbolische Knoten K:
$$\lim_{N \to \infty} \frac{2\pi}{N} \log |J_N(K; e^{2\pi i/N})| = \mathrm{Vol}(S^3 \setminus K)$$
wobei J_N(K; q) das N-farbige Jones-Polynom ist. Murakami–Murakami (2001) zeigten,
dass Kashaevs Invariante gleich dem N-dimensionalen farbigen Jones-Polynom
ausgewertet an N-ter Einheitswurzel ist.

### 4.2 Gescheiterte Ansätze (Jones-Unknoten-Vermutung)

**Kronheimer–Mrowka (2011) — Warum nicht ausreichend**: Kronheimer und Mrowka bewiesen
in *"Khovanov homology is an unknot-detector"* (Publications Mathématiques de l'IHÉS,
2011), dass Khovanov-Homologie Kh̃(K) den Unknoten erkennt: Kh̃(K) hat Rang 1 genau
dann, wenn K = Unknoten. Der Beweis nutzt eine Spektralsequenz von Kh̃(K) zu einer
Knoten-Homologie über singuläre Instantone (Knot Floer Homologie der vernähten
Knotenkomplementierung). Diese Instanton-Homologie erkennt bekanntermaßen den Unknoten.

**Warum dies die Jones-Vermutung NICHT löst**: Das Jones-Polynom V_K(t) ist die
Euler-Charakteristik (gradierte) der Khovanov-Homologie:
$$V_K(t) = \sum_{i,j} (-1)^i t^j \dim \mathrm{Kh}^{i,j}(K)$$
Aber die Euler-Charakteristik von Kh̃ kennt nicht die volle Struktur. V_K(t) = 1
bedeutet, dass diese alternierende Summe 1 ergibt — aber Kh̃(K) könnte nichtrivialen
Rang haben mit gegenseitiger Auslöschung. Ein Knoten mit V_K(t) = 1 könnte also
Kh̃(K) ≠ ℤ haben und trotzdem kein Unknoten sein — das wäre ein Gegenbeispiel zu
Jones, aber kein Widerspruch zu Kronheimer–Mrowka.

**Whitehead-Doppel Wh(T(2,3)) als Kandidat**: Das Whitehead-Doppel des Trefoil-Knotens
T(2,3) hat ein Jones-Polynom das schwer zu berechnen ist. Jedoch ist bekannt, dass
V_{Wh(K)}(t) kompliziert ist und sich von 1 unterscheidet — kein echter Kandidat.

**Virtuelle Knoten**: Für *virtuelle* Knoten existieren tatsächlich nicht-triviale
Beispiele mit Jones-Polynom 1 (Kauffman 2004). Diese liefern aber keine
Gegenbeispiele für klassische Knoten.

**Ansatz über Kauffman-Klammer**: Versuche, Jones' Nullstellen-Struktur auszunutzen
(z. B. an t = e^{2πi/k}) blieben erfolglos.

### 4.3 Aktuelle Ergebnisse zur Volumen-Vermutung

**Bewiesene Spezialfälle**:
- **Torusknoten T(2,n)**: Kashaev–Tirkkonen 2003, vollständig bewiesen.
- **Achterknoten (4₁)**: Numerisch bestätigt, nicht vollständig bewiesen.
- **Viele weitere Knoten**: Numerische Evidenz für alle bekannten hyperbolischen Knoten
  bis ~12 Kreuzungen.

**Garoufalidis–Lê / Garoufalidis–Zagier (2024)**:
- Paper *"Knots, Perturbative Series and Quantum Modularity"* (SIGMA 2024) verbindet
  die Volumen-Vermutung mit **Quantenmodularität**: Das farbige Jones-Polynom nahe
  Einheitswurzeln verhält sich wie eine Quantumodularform.
- Die **Resurgenz** (Borel-Transformation) des perturbativen Chern-Simons-Integrals
  liefert die Volumen-Asymptotik, falls gewisse Stokes-Konstanten nicht verschwinden.

**Garoufalidis–Gu–Mariño–Wheeler (2024/2025)**: *"Resurgence of Chern-Simons Theory
at the Trivial Flat Connection"* (Communications in Mathematical Physics 2024) gibt
eine analytische Fortsetzung des Kashaev-Invariants und komplettiert matrixwertige
holomorphe Quanten-Modularformen. Dies liefert eine **exakte Version der verfeinerten
Quantenmodularitätsvermutung**.

**Gukov–Murakami–Okounkov-Verbindung (SL(2,ℂ)-Chern-Simons)**: Die physikalische
Interpretation der Volumen-Vermutung liegt im komplexen Chern-Simons-Funktional:
Das farbige Jones-Polynom ist das Partitionsfunktional der SL(2,ℂ)-Chern-Simons-
Theorie auf S³\K, und das hyperbolische Volumen ist der dominante Sattelpunkt.

**Categorification (Khovanov, Bar-Natan)**: Der Jones-Unknoten-Ansatz über
Kategorifizierung: Wenn Kh(K) ↔ V_K(t), kann eine vollständige Kategorifizierung der
Jones-Nullstellen Informationen über Unknotenerkennung liefern? Bislang ohne Erfolg
für die direkte Jones-Frage.

### 4.4 Offene Teilprobleme

- **Jones-Unknoten (klassisch)**: Kein Gegenbeispiel, kein Beweis — komplett offen.
- **Volumen-Vermutung (allgemein)**: Für beliebige hyperbolische Knoten offen.
  Die Asymptotik ist numerisch bestätigt, analytisch weitgehend unbewiesen.
- **Komplexe Volumen-Vermutung**: Erweiterung auf Vol(S³\K) + i·CS(S³\K) (Chern-Simons-
  Term) — Murakami et al. 2002, offen.
- **Kashaev-Invariante für nicht-hyperbolische Knoten**: Was ist die Asymptotik für
  Satelliten- und Torus-Knoten (bekannt für Torusknoten, aber nicht für alle)?
- **Jones und virtuelle Knoten**: Warum scheitert die Vermutung für virtuelle Knoten?
  Kann das für klassische Knoten ausgeschlossen werden?
- **Knoten mit unbekanntem Volumen**: Die Vermutung gilt nur für hyperbolische Knoten.
  Was passiert bei Seifert-fibrierten Knoten?
- **Quantenmodularität und Volumen**: Ist die Garoufalidis-Zagier-Quantenmodularform
  ein vollständiger Ersatz für die Volumen-Vermutung?

---

## Zusammenfassung: Methodische Vergleiche und Querverbindungen

| Problem | Kernhindernis | Vielversprechendste Methode |
|---------|--------------|----------------------------|
| MLC (Mandelbrot) | A priori bounds für unbeschränkte ∞-Renorm. | Parabolische Renormalisierung (Dudko–Lyubich) |
| Conway-Scheibigkeit | Mutations-Invarianz aller bekannten Invarianten | Spur-Äquivalenz (Piccirillo), höhere ν-Invarianten |
| Whitehead-Asphärizität | D(2)-Realisierbarkeit algebraischer 2-Zyklen | Diagrammatische Reduktion, Cubulation |
| Jones/Volumen | Keine strukturelle Verbindung Kh ↔ V_K = 1 | Quantenmodularität, Resurgenz, Instanton-Floer |

---

## Literatur und Quellen

1. Dudko, D.; Lyubich, M.: *MLC at Feigenbaum points*, arXiv:2309.02107 (2023/2025)
2. Buff, X.; Chéritat, A.: *Quadratic Julia sets with positive area*, Annals of Mathematics 176/2 (2012), 673–746
3. Piccirillo, L.: *The Conway knot is not slice*, Annals of Mathematics 191/2 (2020), 581–591
4. Kronheimer, P.; Mrowka, T.: *Khovanov homology is an unknot-detector*, Publ. Math. IHÉS 113 (2011), 97–208
5. Johnson, F. E. A.: *Stable Modules and the D(2)-Problem*, Cambridge LMS Lecture Notes 301 (2003)
6. Garoufalidis, S.; Zagier, D.: *Knots, Perturbative Series and Quantum Modularity*, SIGMA 20 (2024), 055
7. Yoccoz, J.-C.: *Polynômes quadratiques et attracteur de Hénon*, Séminaire Bourbaki 734 (1990–91)
8. Avila, A.; Kahn, J.; Lyubich, M.; Shen, W.: *Combinatorial rigidity for unicritical polynomials*, Annals of Mathematics 170/2 (2009), 783–797
9. Kashaev, R. M.: *The hyperbolic volume of knots from the quantum dilogarithm*, Letters in Mathematical Physics 39/3 (1997), 269–275
10. Murakami, H.; Murakami, J.: *The colored Jones polynomials and the simplicial volume of a knot*, Acta Mathematica 186/1 (2001), 85–104
11. Freedman, M. H.: *The topology of four-dimensional manifolds*, J. Differential Geometry 17/3 (1982), 357–453
12. Whitehead, J. H. C.: *On adding relations to homotopy groups*, Annals of Mathematics 42/2 (1941), 409–428
13. Tuzun, R. E.; Sikora, A. S.: *Verification of the Jones unknot conjecture up to 24 crossings*, J. Knot Theory Ramifications (2021)
14. Garoufalidis, S.; Gu, J.; Mariño, M.; Wheeler, C.: *Resurgence of Chern-Simons Theory at the Trivial Flat Connection*, Communications in Mathematical Physics (2024)
15. Howie, J.: *On the asphericity of ribbon disc complements*, Transactions of the AMS 289/1 (1985), 281–302
16. Kapiamba, H.: *Local connectivity of Mandelbrot set at satellite parameters of bounded type*, GAFA 33/4 (2023)
