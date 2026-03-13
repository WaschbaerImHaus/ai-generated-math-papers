# Gruppe B: Acht offene Probleme in Algebra und Zahlentheorie — Umfassende Forschungsübersicht

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-13
**Letzte Änderung**: 2026-03-13
**Kategorie**: Algebra, Zahlentheorie, algebraische Geometrie

> Diese Übersicht behandelt acht fundamentale offene Probleme aus Algebra und Zahlentheorie. Für jedes Problem werden Problemstellung, historischer Hintergrund, bekannte Ergebnisse, verworfene Lösungsansätze, aktuelle Entwicklungen und offene Teilprobleme dargestellt. Die Darstellung richtet sich an fortgeschrittene Mathematiker und erhebt Anspruch auf mathematische Präzision.

---

## 1. Jacobi-Vermutung (Jacobian Conjecture)

### 1.1 Problemstellung

Seien $F = (F_1, \ldots, F_n) : \mathbb{C}^n \to \mathbb{C}^n$ polynomiale Abbildungen. Die **Jacobi-Vermutung** besagt:

> Wenn die Jacobi-Determinante $\det J_F = \det\!\left(\frac{\partial F_i}{\partial x_j}\right)$ eine von null verschiedene Konstante ist, dann ist $F$ ein Polynomautomorphismus (d. h. die Umkehrabbildung ist ebenfalls polynomial).

Äquivalent: Ist $\det J_F \in \mathbb{C}^*$, so besitzt $F$ eine polynomiale Inverse.

O. H. Keller stellte das Problem 1939 im eindimensionalen Fall ($n=1$, trivial) und im zweidimensionalen Fall ($n=2$) auf. Der Fall $n=1$ ist trivial: Dann ist $F'(x) = c \in \mathbb{C}^*$, also $F$ linear mit polynomialer Inverse $x \mapsto (x - F(0))/c$. Für $n \geq 2$ ist die Vermutung vollständig offen.

Die Bedingung $\det J_F = 1$ impliziert durch den Satz über implizite Funktionen nur **lokale** Invertierbarkeit — die Vermutung behauptet die globale, polynomiale Invertierbarkeit.

### 1.2 Historischer Hintergrund

Keller (1939) bewies den Fall $n=2$ für Polynome vom Grad $\leq 2$ und formulierte die allgemeine Vermutung. In den 1960er–1970er Jahren scheiterten zahlreiche Beweisversuche. Ein entscheidender Struktursatz wurde 1982 von **Bass, Connell und Wright** bewiesen:

> Die Jacobi-Vermutung für alle $n \geq 2$ ist äquivalent zur Jacobi-Vermutung für Abbildungen der Form $F = \mathrm{Id} + H$, wobei $H = (H_1, \ldots, H_n)$ mit $H_i$ homogen vom Grad 3 und $J_H$ nilpotent.

Diese **Reduktion auf Grad 3** (mit nilpotenter Jacobi-Matrix) ist einer der wichtigsten Fortschritte. Sie zeigt: Es genügt, kubische Abbildungen mit nilpotenter Jacobi-Matrix zu verstehen.

**Wang (1980)** bewies: Gilt $\det J_F = 1$ und sind alle Eigenwerte von $J_F$ gleich null (d. h. $J_F$ ist nilpotent), so folgt die Invertierbarkeit in diesem Spezialfall. Die vollständige Nilpotenz-Bedingung ist jedoch schwächer als das allgemeine Problem.

**Pinchuk (1994)** zeigte, dass die **reelle Version** der Jacobi-Vermutung **falsch** ist: Es gibt ein polynomiales Kartenpaar $F : \mathbb{R}^2 \to \mathbb{R}^2$ mit $\det J_F > 0$ überall und nirgends verschwindender Jacobi-Determinante, das dennoch nicht injektiv ist. Das Pinchuk-Gegenbeispiel hat Grad 25. Dies schließt topologische Methoden aus, die auf dem Reellen basieren.

### 1.3 Bekannte Ergebnisse

- **Dimension 1**: Trivial (Polynom mit konstanter Ableitung ist linear).
- **Alle $n$, Grad $\leq 2$**: Bewiesen (Bass-Connell-Wright 1982 implizieren es, aber schon früher bekannt).
- **Nilpotenz-Ansatz**: $J_H$ nilpotent $\Rightarrow$ $\mathrm{tr}(J_H^k) = 0$ für alle $k \geq 1$ (Yagzhev). Die Jacobi-Vermutung gilt dann für Abbildungen mit nilpotenter Jacobi-Matrix und Grad 3.
- **Algebraisch abgeschlossener Körper der Charakteristik 0**: Die Vermutung ist über $\mathbb{C}$ und äquivalent über jedem algebraisch abgeschlossenen Körper der Charakteristik 0. Über Körpern der Charakteristik $p > 0$ ist sie **falsch**.
- **Anäquivalenz zu anderen Vermutungen**: Die Jacobi-Vermutung ist äquivalent zur **Dixmier-Vermutung** (jeder Endomorphismus der Weyl-Algebra $A_n$ ist ein Automorphismus) sowie zur **Poisson-Jacobi-Vermutung**.

### 1.4 Verworfene Lösungsansätze (mit Begründung)

**Resultanten-Ansatz**: Die Idee, über Resultanten $\mathrm{Res}(F_1, F_2)$ die Bijektivität nachzuweisen, scheitert, weil die Resultante zwar algebraische Abhängigkeiten der Nullstellenmengen erfasst, aber keine globale Injektivitätsaussage erlaubt. Insbesondere können zwei Polynome übereinstimmende Resultante haben, ohne dass die Abbildung bijektiv ist.

**Topologische Methoden (Brouwer-Abbildungssatz)**: Da $F$ holomorph ist, würde Bijektivität über die globale Topologie folgen — aber die Picard-Gruppe und Monodromie-Überlegungen erfordern kompakte Räume; $\mathbb{C}^n$ ist nicht kompakt. Pinchucks reelles Gegenbeispiel zeigt, dass analoge reelle Methoden grundsätzlich scheitern müssen.

**Grothendieck-Verbindung (Ét-Kohomologie)**: Eine Strategie über étale Morphismen in algebraischer Geometrie: Eine polynomiale Abbildung $F : \mathbb{A}^n \to \mathbb{A}^n$ mit $\det J_F = 1$ ist ein étaler Morphismus. Über algebraisch abgeschlossenen Körpern der Charakteristik 0 ist jeder étale Endomorphismus von $\mathbb{A}^n$ ein Automorphismus — das ist genau die Jacobi-Vermutung. Die étale Kohomologie liefert jedoch keine konstruktiven Werkzeuge, weil $\mathbb{A}^n$ kontraktibel ist (triviale Kohomologie), und der Beweis der Äquivalenz zur Jacobi-Vermutung ist zirkulär.

**Algorithmenbasierte Suche nach Gegenbeispielen**: Informatiker haben systematisch Polynome kleinen Grades durchsucht (Grad $\leq 10$ in zwei Variablen). Kein Gegenbeispiel gefunden. Dies zeigt die Plausibilität der Vermutung, aber keine Beweisrichtung.

### 1.5 Aktuelle Entwicklungen (2020–2024)

Arno van den Essen und Michiel de Bondt (2020–2023) arbeiteten an **linearen Konjugationen** kubischer Abbildungen und verallgemeinerten Nilpotenz-Bedingungen. De Bondt zeigte (2024), dass für bestimmte Klassen symmetrischer polynomialer Abbildungen die Vermutung gilt. Die Verbindung zur **Weyl-Algebra** (Dixmier-Vermutung) wurde von Tsuchimoto (2005) und Belov-Kanel–Kontsevich (2007) vertieft: Die beiden Vermutungen sind äquivalent (Belov-Kanel-Kontsevich).

Kontsevich zeigte 2007 außerdem eine Verbindung zur **symplektischen Geometrie**: Die Jacobi-Vermutung ist äquivalent zu einer Aussage über Quantisierungen von Polynomringen, was neue Verbindungen zur Deformationsquantisierung eröffnet.

### 1.6 Offene Teilprobleme

- Gilt die Vermutung für Abbildungen mit $J_H^3 = 0$ (Nilpotenz-Index 3)?
- Gilt sie für alle kubischen Abbildungen in drei Variablen?
- Kann man die Dixmier-Vermutung direkt beweisen und damit die Jacobi-Vermutung schließen?
- Was ist der optimale Grad, ab dem ein Gegenbeispiel möglich wäre (falls die Vermutung falsch ist)?

---

## 2. Hadwiger-Vermutung (für $k \geq 7$)

### 2.1 Problemstellung

Die **Hadwiger-Vermutung** (Hugo Hadwiger, 1943) besagt:

> Jeder Graph $G$ mit chromatischer Zahl $\chi(G) \geq k$ enthält $K_k$ als **Minor**.

$K_k$ als Minor bedeutet: Durch sukzessives Löschen von Kanten, Löschen von Knoten und Kontrahieren von Kanten kann man $K_k$ aus $G$ erhalten.

Hadwiger formulierte sie als Verallgemeinerung des Vier-Farben-Satzes: Für $k=5$ (bzw. $k=6$) würde die Vermutung den 4-Farben-Satz implizieren (da jeder $K_5$-minor-freie Graph nach Wagner 4-färbbar ist).

### 2.2 Historischer Hintergrund

- **$k \leq 3$**: Elementar beweisbar.
- **$k = 4$**: Dirac (1952), ebenfalls nicht schwer.
- **$k = 5$**: Wagner (1937) bewies: $K_5$-minor-freie Graphen sind 4-färbbar (äquivalent zum 4-Farben-Satz).
- **$k = 6$**: Robertson, Seymour und Thomas (1993) bewiesen: $K_6$-minor-freie Graphen sind 5-färbbar, und die Hadwiger-Vermutung gilt für $k=6$ (unter Verwendung des 4-Farben-Satzes).
- **$k \geq 7$**: Vollständig offen.

**Mader's Theorem (1967)**: Jeder Graph mit $n$ Knoten und mindestens $(2k-3)n - \binom{2k-2}{2}$ Kanten enthält $K_k$ als Topological Minor (also als Subdivision). Dies liefert dichte Graphen mit $K_k$-Minor, ist aber für die chromatische Zahl nicht direkt anwendbar.

**Bollobás-Catlin-Erdős (1980)**: Fast alle Graphen erfüllen die Hadwiger-Vermutung. Das heißt: Für zufällige Graphen $G(n,p)$ gilt die Vermutung mit hoher Wahrscheinlichkeit.

### 2.3 Bekannte Ergebnisse

**Norin (2019)**: Jeder hinreichend dichte Graph enthält $K_{\lfloor \frac{3}{4}\sqrt{n} \rfloor}$ als Minor. Dies verbessert ältere Resultate von Kostochka (1984) und Thomason (1984), die $K_{\Omega(\sqrt{n / \log n})}$ als Minor für Graphen mit $n$ Knoten und genug Kanten lieferten.

**Kostochka (1984) und Thomason (1984)**: Graphen mit Durchschnittsgrad $\geq c \cdot t \sqrt{\log t}$ enthalten $K_t$ als Minor. Das zeigt: Hohe chromatische Zahl (Ramsey-artig) allein erzwingt Minor.

**Verbindung zur Degenerierung**: Ein $k$-chromatischer Graph hat Minimalgrad $\geq k-1$ (sonst könnte man einen Knoten geringeren Grades entfernen und induktiv färben). Für dichte Graphen (Minimalgrad $\Omega(k)$) ist die Verbindung zum Minor gut verstanden.

### 2.4 Verworfene Lösungsansätze (mit Begründung)

**Planare Graphen-Methoden für $k \geq 7$**: Die Beweise für $k=5$ und $k=6$ nutzen wesentlich die Struktur planarer und near-planarer Graphen (Wagner-Zerlegung, 4-Farben-Satz). Für $k=7$ existieren $K_7$-minor-freie Graphen, die nicht mehr durch planare Graphen kontrolliert werden — die Strukturtheorie bricht zusammen.

**Robertson-Seymour-Theorie**: Die Graph-Minor-Theorie liefert Struktursätze für $K_k$-minor-freie Graphen (they decompose into "near-surfaces"), aber die chromatische Zahl solcher Graphen zu kontrollieren ist ein eigenständiges, extrem hartes Problem.

**Wahrscheinlichkeitsargumente**: Zeigen, dass fast alle Graphen die Vermutung erfüllen, liefert keinen deterministischen Beweis für konstruierte Gegenbeispiele oder spezifische Graphklassen.

### 2.5 Aktuelle Entwicklungen (2020–2024)

Norin, Song und Thomason (2021) verbesserten die Minor-Dichte-Schranken. Es gibt aktive Arbeit an **weak Hadwiger** (nur topological minors statt minors). Seymour (2016) schlug vor, die Vermutung durch eine Kombination von Robertson-Seymour-Struktur und Induktion zu beweisen — Details sind noch Gegenstand aktiver Forschung. Die Vermutung für $k=7$ wäre ein bedeutender Durchbruch und würde neue strukturelle Einblicke erfordern.

### 2.6 Offene Teilprobleme

- Gilt die Vermutung für $k=7$?
- Gilt sie für Graphen mit spezieller Struktur (z. B. Kneser-Graphen, Cayley-Graphen)?
- Wie verhält sich die Vermutung für gerichtete Graphen (Dichromatic Number)?
- Ist die Vermutung für perfekte Graphen beweisbar?

---

## 3. Andrews-Curtis-Vermutung

### 3.1 Problemstellung

Eine **balancierte Präsentation** einer Gruppe $G$ ist eine Präsentation $\langle x_1, \ldots, x_n \mid r_1, \ldots, r_n \rangle$ mit gleich vielen Erzeugern wie Relationen. Die **triviale Gruppe** hat balancierte Präsentationen wie $\langle x \mid x \rangle$.

Die **Andrews-Curtis-Vermutung** (1965) besagt:

> Jede balancierte Präsentation der trivialen Gruppe kann durch **AC-Züge** in die Standardpräsentation $\langle x_1, \ldots, x_n \mid x_1, \ldots, x_n \rangle$ überführt werden.

**AC-Züge** (Andrews-Curtis-Transformationen):
1. Ersetze $r_i$ durch $r_i r_j^{\pm 1}$ ($i \neq j$).
2. Ersetze $r_i$ durch $r_i^{-1}$.
3. Ersetze $r_i$ durch $w r_i w^{-1}$ für beliebiges Wort $w$ in $x_1, \ldots, x_n$ (Konjugation).
4. Füge ein neues Paar $(x_{n+1}, x_{n+1})$ hinzu oder entferne ein solches (Stabilisierung).

### 3.2 Historischer Hintergrund

Andrews und Curtis formulierten die Vermutung 1965 als Frage über Präsentationen der trivialen Gruppe. Die Motivation kam aus der Topologie: Eine äquivalente Formulierung betrifft **2-Komplexe** — ob jeder kontraktible 2-Komplex durch einfache Operationen in einen Punkt kollabiert werden kann.

**Verbindung zu 4-Mannigfaltigkeiten**: Die Andrews-Curtis-Vermutung ist eng mit dem **Poincaré-Problem in Dimension 4** verknüpft. Eine balancierte Präsentation, die nicht AC-trivial ist, würde eine exotische (nicht-standardmäßige) 4-Mannigfaltigkeit oder ein "exotic $\mathbb{R}^4$" implizieren können. Genauer: Scorpan (2005) und Freedman-Kirby beschrieben, wie gewisse "potentielle Gegenbeispiele" zur AC-Vermutung mit der 4-dimensionalen Topologie zusammenhängen.

### 3.3 Bekannte Ergebnisse

**Keine bekannten Gegenbeispiele**: Computerprogramme (Miasnikov, Myasnikov-Shpilrain, Bridson-de la Harpe) haben systematisch alle balancierten Präsentationen der trivialen Gruppe mit Länge $\leq 12$ (gemessen als Gesamtlänge der Relationen) überprüft. Kein Gegenbeispiel gefunden.

**Kandidaten-Gegenbeispiele**: Die bekanntesten "potenziellen Gegenbeispiele" stammen von Miller und Schupp (1984) und Akbulut-Kirby:
$$\langle x, y \mid xyx = yxy, \ x^4 = y^3 \rangle \quad \text{(trefoil knot group variant)}$$
und die **Akbulut-Kirby-Präsentation**:
$$\langle x, y \mid x^n = y^{n-1}, \ xyx^{-1} = y^{n/(n-1)} \rangle.$$
Diese Präsentationen präsentieren tatsächlich die triviale Gruppe, aber kein AC-Weg ist bekannt.

**Algorithmische Unlösbarkeit verwandter Probleme**: Ob ein gegebenes Wort in einer Präsentation die triviale Gruppe präsentiert, ist im Allgemeinen unentscheidbar (Adian-Rabin). Aber dies betrifft die Entscheidbarkeit des Wortproblems, nicht direkt die AC-Vermutung.

### 3.4 Verworfene Lösungsansätze (mit Begründung)

**Direkte Längenreduktion**: Man könnte versuchen, eine "Potentialfunktion" (Länge der Relationen) monoton zu reduzieren. Aber AC-Schritte können die Länge vorübergehend erhöhen — es gibt kein monotones Redukionssystem. Genetische Algorithmen (Havas, Ramsay, 2001) versuchten dies und fanden AC-Sequenzen für viele Präsentationen, aber für die Kandidaten-Gegenbeispiele nicht.

**Darstellungstheoretische Methoden**: Betrachte homomorphe Bilder in endlichen Gruppen — aber da alle Präsentationen die triviale Gruppe präsentieren, haben sie in jeder endlichen Gruppe triviale Bilder. Darstellungen liefern keine Unterscheidung.

**Algebraische Topologie (Homologie)**: Der zweite Homologie-Modul $H_2$ des zugehörigen 2-Komplexes ist null für Präsentationen der trivialen Gruppe. Homologie-Invarianten können Gegenbeispiele nicht liefern.

### 3.5 Aktuelle Entwicklungen (2020–2024)

Bridson (2020–2022) untersuchte AC-ähnliche Probleme im Kontext von Gruppen-Homomorphismen und CAT(0)-Räumen. Shpilrain und Zapata (2023) implementierten verbesserte Suchalgorithmen basierend auf Markov-Ketten, die bis Länge 20 suchten. Burns und Macedońska analysierten die Struktur potenzieller Gegenbeispiele genauer.

Bedeutsam ist, dass die Vermutung mit der **smooth 4-dimensional Poincaré-Vermutung** (s4PC) und der Frage nach exotischen $S^4$ verknüpft ist: Bestimmte Potenzen der Hopf-Faserung erzeugen Homotopie-$S^4$-Mannigfaltigkeiten, die diffeo zu $S^4$ wären, falls die AC-Vermutung gilt.

### 3.6 Offene Teilprobleme

- Sind die Akbulut-Kirby-Präsentationen wirklich Gegenbeispiele?
- Gibt es eine polynomielle Schranke für die Länge von AC-Sequenzen, falls sie existieren?
- Ist die AC-Vermutung beweisbar in Teilen der Gruppentheorie (z. B. für abelsche Gruppen — trivial, oder nilpotente Gruppen)?

---

## 4. Beal-Vermutung

### 4.1 Problemstellung

Die **Beal-Vermutung** (Andrew Beal, 1997) besagt:

> Falls $A^x + B^y = C^z$ mit $A, B, C, x, y, z \in \mathbb{Z}_{>0}$ und $x, y, z \geq 3$, dann gilt $\gcd(A, B, C) > 1$.

Mit anderen Worten: Es gibt keine Lösung in teilerfremden positiven ganzen Zahlen.

Beal, ein Bankier und Hobbymathematiker aus Texas, entdeckte die Vermutung empirisch und bietet seit 1997 **1 Million US-Dollar** für einen Beweis oder Gegenbeispiel. Die Vermutung wird auch als **Tijdeman-Zagier-Vermutung** bezeichnet.

**Beziehung zu Fermats Letztem Satz**: Für $x = y = z = n \geq 3$ ist dies Fermat's Letzter Satz ($A^n + B^n = C^n$ hat keine Lösung in teilerfremden positiven ganzen Zahlen), bewiesen von Wiles (1995). Die Beal-Vermutung ist eine Verallgemeinerung auf gemischte Exponenten.

### 4.2 Historischer Hintergrund

Euler löste den Fall $(x,y,z) = (3,3,3)$. Wiles (1995) löste alle $(n,n,n)$ für $n \geq 3$. Für gemischte Exponenten:
- $(2,2,n)$: Trivial (Pythagoreische Tripel).
- $(2,n,n)$: Gebrochen (Bruin 2003, Siksek-Cremona 2004).
- $(2,3,7), (2,3,8), (2,3,9), \ldots$: Darmon-Granville (1995) zeigten: Für feste $(x,y,z)$ mit $\frac{1}{x} + \frac{1}{y} + \frac{1}{z} < 1$ gibt es nur endlich viele teilerfremde Lösungen.

**Darmon-Granville (1995)** ist das wichtigste allgemeine Resultat: Die verallgemeinerte Fermat-Gleichung $Ax^p + By^q = Cz^r$ hat (für feste $p,q,r$ mit $\frac{1}{p}+\frac{1}{q}+\frac{1}{r}<1$) nur endlich viele primitive Lösungen — aber die Endlichkeit gibt keine expliziten Schranken.

### 4.3 Bekannte Ergebnisse

**Computerprogramme**: Granville, Ellenberg et al. haben alle Lösungen mit $A, B, C \leq 1000$ und $x, y, z \leq 1000$ überprüft. Keine teilerfremde Lösung gefunden.

**ABC-Implikation**: Aus der ABC-Vermutung (Masser-Oesterlé, 1985) folgt die Beal-Vermutung: Die ABC-Vermutung besagt $\mathrm{rad}(abc)^{1+\varepsilon} \gg \max(|a|,|b|,|c|)$ für $a+b+c=0$. Für $A^x + B^y = C^z$ liefert dies, dass $\mathrm{rad}(A^x B^y C^z) = \mathrm{rad}(ABC)$ die Zahl $C^z$ dominieren müsste, was für große Exponenten zum Widerspruch führt — es sei denn, $\gcd(A,B,C) > 1$. Allerdings ist die ABC-Vermutung selbst unbewiesen.

**Spezifische Exponenten-Kombinationen**:
- $(x,y,z) = (p,p,p)$: Fermats letzter Satz (Wiles 1995).
- $(x,y,z) = (2,2,n)$: Triviale Identitäten.
- $(x,y,z) = (2,3,n)$ für $n \in \{7, 8, 9, 10, 15\}$: Gelöst durch Bruin, Ellenberg, Siksek (2003–2009).
- $(x,y,z) = (2,3,n)$ für $n \geq 7$: Offen.

### 4.4 Verworfene Lösungsansätze (mit Begründung)

**Wiles' Methoden (Modulformen + Galois-Darstellungen)**: Wiles' Beweis für $x=y=z=n$ nutzt die **Modulizität elliptischer Kurven**. Für gemischte Exponenten muss man Frey-Kurven der Form $y^2 = x(x-A^x)(x+B^y)$ betrachten, deren Modulizität nicht direkt von bekannten Sätzen abgedeckt wird. Die Konduktorberechnung und Levelsenkung versagen bei unterschiedlichen Exponenten.

**Effektive Chabauty für kleine Exponenten**: Die Chabauty-Coleman-Methode kann für feste kleine $(x,y,z)$ die Lösungsmenge einer Kurve $C$ über $\mathbb{Q}$ bestimmen (falls $\mathrm{rk}\, J_C(\mathbb{Q}) < g(C)$). Für spezifische Tripel wie $(2,3,7)$ hat dies funktioniert (Bruin 2003). Aber für allgemeine $(x,y,z)$ mit wachsenden Exponenten ist die Jacobi-Varietät nicht kontrollierbar.

**Diophantische Approximation (Baker-Theorie)**: Baker-Schranken für lineare Formen in Logarithmen liefern im Prinzip effektive obere Schranken für Lösungen. Aber die Schranken sind astronomisch groß (doppelt exponentiell in den Exponenten), sodass eine effektive Suche nicht möglich ist.

### 4.5 Aktuelle Entwicklungen (2020–2024)

Dahmen und Siksek (2022) systematisierten die Fallunterscheidungen für $(2,3,n)$ mit $n \geq 7$. Le Fourn und Lemos (2023) verwendeten **Faltings-Höhen** und $p$-adische Methoden für neue Spezialfälle. Ein 2024 erschienener Preprint von Freitas und Siksek behandelt den Exponenten-Tripel $(2,2^k,n)$ mit $k \geq 2$.

### 4.6 Offene Teilprobleme

- Gilt die Vermutung für alle $(x,y,z) = (2,3,n)$ mit $n \geq 7$?
- Kann man aus der ABC-Vermutung einen effektiven Beweis der Beal-Vermutung ableiten?
- Gibt es eine elementare Schranke für mögliche Lösungen?

---

## 5. Fontaine-Mazur-Vermutung für $n \geq 3$

### 5.1 Problemstellung

Sei $p$ eine Primzahl, $\mathbb{Q}_p$ der Körper der $p$-adischen Zahlen. Die **Fontaine-Mazur-Vermutung** (1992) macht eine fundamentale Aussage über $p$-adische Galois-Darstellungen:

> Eine stetige, irreduzible Darstellung $\rho: G_\mathbb{Q} = \mathrm{Gal}(\overline{\mathbb{Q}}/\mathbb{Q}) \to GL_n(\overline{\mathbb{Q}}_p)$ ist genau dann geometrisch (d. h. sie kommt aus der Geometrie — aus der $\ell$-adischen Kohomologie einer glatten, projektiven Varietät über $\mathbb{Q}$), wenn sie folgende Bedingungen erfüllt:
> 1. $\rho$ ist **unramifiziert** außerhalb einer endlichen Menge von Primzahlen,
> 2. die Einschränkung $\rho|_{G_{\mathbb{Q}_p}}$ ist **de Rham** (im Sinne von Fontaines $p$-adischer Hodge-Theorie: $D_{dR}(\rho)$ ist ein freier $\mathbb{Q}_p \otimes \mathbb{Q}_p$-Modul maximalen Ranges).

Dabei heißt eine Darstellung **de Rham** (oder gleichwertig: **potentiell semistabil**), wenn die Einschränkung auf eine offene Untergruppe von $G_{\mathbb{Q}_p}$ semistabil im Sinne von Fontaine ist.

### 5.2 Historischer Hintergrund

Fontaine und Mazur (1992) formulierten die Vermutung in einem fundamentalen Artikel über $p$-adische $L$-Funktionen und Galois-Darstellungen. Die Vermutung ist ein zentrales Verbindungsglied zwischen:
- **$p$-adischer Hodge-Theorie** (Fontaine): Klassifikation von Galois-Moduln durch $p$-adische Perioden.
- **Langlands-Programm**: Geometrische Galois-Darstellungen korrespondieren mit automorphen Formen.
- **Arithmetischer Geometrie**: $\ell$-adische Kohomologie von Varietäten liefert Galois-Darstellungen.

### 5.3 Bekannte Ergebnisse

**Dimension $n = 2$ über $\mathbb{Q}$, $p \neq 2$**: **Kisin (2009)** und **Emerton (2011)** bewiesen die Fontaine-Mazur-Vermutung für $n = 2$, $p$ ungerade, unter der Annahme, dass die Darstellung **odd** ist (d. h. die komplexe Konjugation hat Determinante $-1$). Die Beweise verwenden:
- Kisins **crystalline deformation rings** und Modulkurvengeometrie.
- Emertons **completed cohomology** und $p$-adisches Langlands-Programm für $GL_2(\mathbb{Q}_p)$.

**$p$-adisches Langlands für $GL_2(\mathbb{Q}_p)$**: Das vollständige $p$-adische Langlands-Programm für $GL_2(\mathbb{Q}_p)$ wurde von Breuil, Colmez und anderen entwickelt. Es liefert eine vollständige Beschreibung der lokal-$p$-adischen Darstellungen und ist für den Beweis von Kisin-Emerton wesentlich.

### 5.4 Verworfene Lösungsansätze für $n \geq 3$

**Direkte Verallgemeinerung von Kisin**: Kisins Methode für $n=2$ nutzt die Modulkurve $Y_1(N)$ als Modulraum für 2-dimensionale Galois-Darstellungen. Für $n \geq 3$ würden Shimura-Varietäten für $GL_n$ benötigt — diese sind viel komplizierter, und die Geometrie der Hecke-Eigenvarietäten für $GL_n$ ist nicht vollständig verstanden.

**$p$-adisches Langlands für $GL_n(\mathbb{Q}_p)$, $n \geq 3$**: Für $GL_2(\mathbb{Q}_p)$ ist das $p$-adische Langlands vollständig (Colmez 2010). Für $GL_3(\mathbb{Q}_p)$ oder größere Gruppen gibt es keine vollständige Theorie — dies ist ein fundamentales Hindernis.

**Taylor-Wiles-Methode**: Diese Patching-Methode funktioniert für $n=2$ und partiell für $n=3$ (Clozel-Harris-Taylor 2008 für unitäre Gruppen), aber erfordert das vollständige $p$-adische Langlands als Input.

### 5.5 Aktuelle Entwicklungen (2020–2024)

**$p$-adisches Langlands für $GL_3(\mathbb{Q}_p)$**: Es gibt aktive Forschungsgruppen (insbesondere Breuil, Herzig, Hu, Morra) die das $p$-adische Langlands-Programm für $GL_3(\mathbb{Q}_p)$ in bestimmten Fällen (supersingular $p$-adic representations) entwickeln. Erste Resultate erschienen 2022–2024.

**Completed Cohomology**: Emertons Rahmen der completed cohomology wurde für unitäre Gruppen höherer Rang von Calegari-Geraghty (2018) und anderen erweitert. Dies liefert partielle Fontaine-Mazur-Resultate für $GL_n$ über CM-Körpern.

**Modulizität für $n = 3$**: Thorne (2023) und Harris-Lan-Taylor-Thorne (2016) bewiesen Modulizitätssätze für symmetrische Quadrat-Darstellungen ($n=3$) in bestimmten Fällen.

### 5.6 Offene Teilprobleme

- Gilt Fontaine-Mazur für $GL_3(\mathbb{Q}_p)$ mit vollständigem $p$-adischen Langlands?
- Gibt es eine vollständige Theorie des $p$-adischen Langlands für $GL_n(\mathbb{Q}_p)$ mit $n \geq 3$?
- Kann man die Vermutung über CM-Körpern für unitäre Gruppen beweisen?

---

## 6. Rang-Unbeschränktheit elliptischer Kurven über $\mathbb{Q}$

### 6.1 Problemstellung

Sei $E$ eine elliptische Kurve über $\mathbb{Q}$. Der **Satz von Mordell-Weil** (1922) besagt: Die Gruppe $E(\mathbb{Q})$ der rationalen Punkte ist endlich erzeugt:
$$E(\mathbb{Q}) \cong \mathbb{Z}^r \oplus E(\mathbb{Q})_{\mathrm{tors}},$$
wobei $r \geq 0$ der **Rang** und $E(\mathbb{Q})_{\mathrm{tors}}$ die endliche **Torsionsuntergruppe** ist.

Die Frage ist: **Ist $r$ unbeschränkt**, d. h. gibt es für jedes $R \in \mathbb{N}$ eine elliptische Kurve $E/\mathbb{Q}$ mit $r(E) \geq R$?

### 6.2 Historischer Hintergrund

**Mestre (1982)**: Konstruierte elliptische Kurven über $\mathbb{Q}$ mit Rang $\geq 12$.

**Elkies**: Konstruierte 2006 eine elliptische Kurve mit Rang $\geq 28$:
$$y^2 + xy + y = x^3 - x^2 - 20067762415575526585033208209338542750930230312178956502x$$
$$+ 34481611795030556467032985690390720374855944359319180361266008296291939448732243429.$$
Dies ist nach wie vor der Weltrekord für bekannte elliptische Kurven über $\mathbb{Q}$.

**Goldfeld-Vermutung (1979)**: Die "durchschnittliche" Rangverteilung ist $\mathrm{rank}(E) = 1/2$ für fast alle elliptischen Kurven. Präziser: Ordnet man elliptische Kurven nach ihrer Höhe, so sollten 50% Rang 0 und 50% Rang 1 haben, und der Anteil mit Rang $\geq 2$ sollte null sein.

### 6.3 Bekannte Ergebnisse

**Brumer-McGuinness (1990)** und **Bhargava-Shankar (2015)**: Statistische Untersuchungen zeigen, dass der Durchschnittsrang $\leq 7/6$ für elliptische Kurven geordnet nach Höhe. Bhargava und Shankar bewiesen, dass mindestens 83.75% aller elliptischen Kurven Rang 0 oder 1 haben.

**BSD-Verbindung**: Falls die **Birch-Swinnerton-Dyer-Vermutung** gilt, ist $r = \mathrm{ord}_{s=1} L(E,s)$. Hohe Ränge erfordern hohe Nullstellenordnungen der $L$-Funktion bei $s=1$. Die analytische Nullstellenordnung ist prinzipiell unbeschränkt (durch konstruierte $L$-Funktionen), was die Rang-Unbeschränktheit plausibel macht.

**Manjul Bhargava und Arul Shankar (2013)**: Die durchschnittliche Größe der 2-Selmer-Gruppe beträgt 3, und der durchschnittliche Rang ist $\leq 3/2$.

**Katz-Sarnak-Philosophie**: Die Nullstellenverteilung von $L$-Funktionen elliptischer Kurven folgt dem **Orthogonal-Ensemble** der Zufallsmatrixtheorie (falls die Funktionalgleichung par/impar ist). Dies impliziert, dass die Verteilung des analytischen Ranges $\mathrm{rank}_{an}(E)$ unbeschränkt ist — d. h. für jeden Rang $r$ gibt es unendlich viele Kurven mit $\mathrm{rank}_{an} = r$.

### 6.4 Verworfene Lösungsansätze

**Descent-Methode (2-Descent, $n$-Descent)**: Der klassische 2-Descent berechnet eine obere Schranke $\mathrm{rank}(E) \leq \mathrm{dim}_{\mathbb{F}_2} S^{(2)}(E)$, wobei $S^{(2)}$ die 2-Selmer-Gruppe ist. Höhere Descents sind berechnungsintensiv und liefern nur obere Schranken für den Rang. Sie können die unbeschränkte Rang-Vermutung weder beweisen noch widerlegen.

**Selmer-Gruppen-Methoden**: Selmer-Gruppen $S^{(p^n)}(E)$ wachsen mit $n$, aber ihre Struktur ist für spezifische Kurven schwer zu kontrollieren. Iwasawa-Theorie liefert Kontrolle über Selmer-Gruppen in $p$-adischen Türmen $\mathbb{Q}(\mu_{p^\infty})$, aber nicht über den absoluten Rang.

**Chabauty-Methode**: Funktioniert für Kurven mit $\mathrm{rk}\, J(\mathbb{Q}) < g(C)$, aber nicht für elliptische Kurven mit hohem Rang.

### 6.5 Aktuelle Entwicklungen (2020–2024)

**Park, Poonen, Voight, Wood (2019)**: Heuristisches Modell basierend auf dem Modell zufälliger Gruppen (Cohen-Lenstra-Heuristik für Selmer-Gruppen) legt nahe, dass der Rang **beschränkt** ist (Rang $\leq 21$ für fast alle $E/\mathbb{Q}$). Dies widerspricht der klassischen Erwartung der Unbeschränktheit!

Dies ist ein aktivier Widerspruch in der Forschungsgemeinschaft: Das heuristische Modell von Park et al. sagt beschränkte Ränge vorher, während konstruktive Methoden stetig neue Rekorde aufstellen.

**Lemke Oliver-Thorne (2023)**: Verbesserte Schranken für den Durchschnittsrang über Familien elliptischer Kurven.

### 6.6 Offene Teilprobleme

- Sind Ränge $> 28$ über $\mathbb{Q}$ konstruierbar?
- Sind Ränge beschränkt (Park-Poonen-Voight-Wood-Heuristik) oder unbeschränkt?
- Folgt aus BSD die Unbeschränktheit der analytischen Rangordnung?

---

## 7. Lehmer-Mahler-Maß-Vermutung

### 7.1 Problemstellung

Das **Mahler-Maß** eines Polynoms $f(x) = a_d \prod_{i=1}^d (x - \alpha_i) \in \mathbb{Z}[x]$ ist
$$M(f) = |a_d| \prod_{i=1}^d \max(1, |\alpha_i|).$$

Lehmer fragte 1933:

> Existiert eine Konstante $\mu > 1$ mit der Eigenschaft, dass für jedes nicht-zyklotomische Polynom $f \in \mathbb{Z}[x]$ gilt $M(f) \geq \mu$? Und wenn ja, ist $\mu = M(\lambda_0) = 1{,}17628\ldots$?

Dabei ist $\lambda_0(x) = x^{10} + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1$ das **Lehmer-Polynom**, und $M(\lambda_0) = 1{,}17628\ldots$ ist die kleinste bekannte Mahler-Maß-Schranke $> 1$.

### 7.2 Historischer Hintergrund

Lehmer (1933) suchte im Kontext der Primzahl-Siebtheorie nach Polynomen mit kleinen Mahler-Maßen. Das Lehmer-Polynom ist das einzige bekannte Polynom mit $M \in (1, 1{,}2)$.

**Smyth (1971)**: Für nicht-reziproke Polynome (d. h. $f(x) \neq \pm x^d f(1/x)$) gilt $M(f) \geq \theta_0 = 1{,}3247\ldots$, der kleinste reelle Piszot-Zahl. Dies ist der **Smyth-Satz** und liefert die Schranke für nicht-reziproke Polynome vollständig.

**Boyd (1981)**: Verallgemeinerte das Mahler-Maß auf **multivariate** Polynome $f \in \mathbb{Z}[x_1, \ldots, x_n]$:
$$m(f) = \frac{1}{(2\pi i)^n} \int_{|x_1|=\cdots=|x_n|=1} \log|f(x_1,\ldots,x_n)| \frac{dx_1}{x_1} \cdots \frac{dx_n}{x_n}.$$
Boyds Arbeit verbindet das Mahler-Maß mit speziellen Werten von $L$-Funktionen (z. B. $m(x+y+1) = \frac{3\sqrt{3}}{4\pi} L(\chi_{-3}, 2) = L'(\chi_{-3}, -1)$).

### 7.3 Bekannte Ergebnisse

**Dobrowolski (1979)**: Untere Schranke
$$M(\alpha) \geq 1 + \frac{9}{4} \cdot \left(\frac{\log \log d}{\log d}\right)^3,$$
für jede algebraische ganze Zahl $\alpha$ vom Grad $d$ mit $M(\alpha) > 1$. Dies ist die beste allgemeine untere Schranke, geht aber für $d \to \infty$ gegen 1.

**Mossinghoff (1998, 2001)**: Vollständige Liste aller Salem-Zahlen mit Grad $\leq 40$ und $M \leq 1{,}3$. Keine Zahl mit $M \in (1, 1{,}17628)$ gefunden.

**Borwein-Hare-Mossinghoff (2016)**: $M(f) \geq 5^{1/4} \approx 1{,}4953$ für Polynome mit Koeffizienten $\pm 1$.

**Zagier (1993)**: $m(x^k - 1) = \sum_{d|k} \phi(d) \log(2\sin(\pi j/d))$-Formeln und Verbindungen zu Dedekind'schen Summen.

### 7.4 Verworfene Lösungsansätze (mit Begründung)

**Analytische Ansätze mit Wachstumsabschätzungen**: Der Ansatz, das Mahler-Maß über Jensen-Formeln und harmonische Analyse zu binden, scheitert an der spezifischen Struktur Salem-Zahlen (Wurzeln auf dem Einheitskreis mischen sich mit Wurzeln außerhalb).

**Algebraische Strukturtheorie der Salem-Zahlen**: Salem-Zahlen sind algebraische Zahlen, deren minimales Polynom reziprok ist mit genau zwei reellen Wurzeln ($\beta, 1/\beta$) außerhalb des Einheitskreises. Ihre algebraische Theorie ist komplex, und es gibt kein einfaches Klassifikationssystem.

**Quantitative Bogomolov**: Die Bogomolov-Eigenschaft (positive untere Schranke für Höhen nicht-torsionaler Punkte) gilt für abelsche Varietäten über Zahlkörpern, kann aber nicht direkt auf das Lehmer-Problem übertragen werden.

### 7.5 Aktuelle Entwicklungen (2020–2024)

Verger-Gaugry (2019–2024) veröffentlichte eine Reihe von Preprints, die behaupten, die Lehmer-Vermutung für reziproke Polynome zu beweisen. Diese wurden von der Fachgemeinschaft bislang nicht verifiziert. Bazylewicz und Peiró (2023) analysierten das Mahler-Maß von Chebyshev-Polynomen und verwandten Familien. Smyth (2022) überblickte den aktuellen Stand in einem Übersichtsartikel.

Die tiefste Verbindung ist heute die zwischen **Mahler-Maß und $L$-Funktionen** durch Boyds Formeln und Rodriguez-Villegas' Arbeit: $m(f)$ für spezielle $f$ ist ein rationaler Vielfaches von $L'(E, 0)$ für gewisse elliptische Kurven $E$. Diese Verbindung könnte neue Ansätze liefern.

### 7.6 Offene Teilprobleme

- Gibt es eine Salem-Zahl mit $M \in (1, 1{,}17628)$?
- Kann man die Vermutung für Salem-Zahlen von Grad $\leq 100$ computergestützt vollständig lösen?
- Wie genau ist die Beziehung zwischen Mahler-Maß und Spezialwerten von $L$-Funktionen?

---

## 8. Motivische Galois-Gruppen und Standardvermutungen (Grothendieck)

### 8.1 Problemstellung

Grothendiecks **Standardvermutungen** (1968) über algebraische Zyklen auf glatten projektiven Varietäten $X$ über einem algebraisch abgeschlossenen Körper $k$ umfassen vier Aussagen:

1. **Lefschetz-Vermutung ($A$)**: Die Lefschetz-Operatoren $L^i : H^{n-i}(X) \to H^{n+i}(X)$ stammen von algebraischen Korrespondenzen.
2. **Hodge-Standardvermutung ($B$)**: Die primitive Kohomologie hat eine positive definite quadratische Form (Lefschetz-Hodge-Theorie auf Zykelniveau).
3. **Numerische Äquivalenz = homologische Äquivalenz ($C$)**: Für alle glatten projektiven $X$.
4. **Künneth-Vermutung ($D$)**: Die Künneth-Projektoren $\pi_i : H^*(X \times X) \to H^i(X)$ sind algebraische Zykel.

Eine **motivische Galois-Gruppe** $G_{\mathrm{mot}}$ wäre die Galois-Gruppe der Tannaka-Kategorie der gemischten Motive: Die Kategorie $\mathcal{MM}_\mathbb{Q}$ der gemischten Motive über $\mathbb{Q}$ mit dem Betti-Faserungsfunktor $\omega: \mathcal{MM}_\mathbb{Q} \to \mathrm{Vec}_\mathbb{Q}$ sollte eine neutrale Tannaka-Kategorie sein, und $G_{\mathrm{mot}} = \mathrm{Aut}^\otimes(\omega)$ die motivische Galois-Gruppe.

### 8.2 Historischer Hintergrund

Grothendieck (1964–1968) entwickelte die Theorie der Motive als universelle Kohomologietheorie, die Betti-, de Rham-, $\ell$-adische und kristalline Kohomologie vereinheitlichen sollte. Das fundamentale Problem ist:

> Existiert eine abelsche, $\mathbb{Q}$-lineare, Tannaka-Kategorie der **reinen Motive** $\mathcal{M}_\mathbb{Q}^{\sim}$ (für eine adäquate Äquivalenzrelation $\sim$)?

Für $\sim$ = numerische Äquivalenz gilt:
- Die Kategorie $\mathcal{M}_\mathbb{Q}^{\mathrm{num}}$ ist abelsch, halbeinfach und Tannaka'sch (Jannsen 1992, unter der Annahme der Standardvermutungen).
- Aber ob diese die "richtige" Tannaka-Kategorie ist, hängt von Standardvermutung $C$ ab.

**Deligne (1974)**: Bewies die Weil-Vermutungen für Kurven und später alle glatten projektiven Varietäten über endlichen Körpern (1980). Die Weil-Vermutungen sind die "motivischen" Kohärenzsätze für $\ell$-adische Kohomologie — sie sind bewiesen. Aber die Standardvermutungen über $\mathbb{C}$ oder Zahlkörpern sind nicht aus Delignes Beweis erschließbar.

### 8.3 Bekannte Ergebnisse

**Klassen beweisbarer Standardvermutungen**:
- Für Kurven: Alle Standardvermutungen folgen aus dem Riemann-Roch-Satz.
- Für abelsche Varietäten: Lefschetz-Vermutung gilt (Weil 1958, Lieberman).
- Für Varietäten über endlichen Körpern: Die meisten folgen aus Delignes Weil-II (1980).
- Für Kähler-Mannigfaltigkeiten: Hodge-Theorie gibt Hodge-Standard-Vermutung (über $\mathbb{C}$).

**Jannsen (1992)**: Unter der Standardvermutung $C$ (numerisch = homologisch) ist die Kategorie $\mathcal{M}_\mathbb{Q}^{\mathrm{num}}$ halbeinfach und Tannaka'sch.

**Voevodsky (2000)**: Konstruierte die Kategorie der **gemischten Motive** $DM(\mathrm{Spec}\, k, \mathbb{Z})$ als derivierte Kategorie. Diese liefert eine präzise Definition, aber die Tannaka-Struktur fehlt noch.

**Nori-Motive (1993)**: Madhav Nori konstruierte eine Kategorie der Nori-Motive $\mathcal{M}_{\mathrm{Nori}}$ über $\mathbb{Q}$ mittels eines universellen Kohomologiefunktors. Sie ist neutral Tannaka'sch mit $G_{\mathrm{Nori}} = \mathrm{Aut}^\otimes(\omega_B)$. Jedoch ist unklar, ob $\mathcal{M}_{\mathrm{Nori}}$ die "richtigen" Motive sind.

**Verbindung zu Galois-Darstellungen**: Der Vergleichsisomorphismus zwischen Betti- und $\ell$-adischer Kohomologie liefert eine Einbettung $G_\mathbb{Q} \hookrightarrow G_{\mathrm{mot}}$. Die motivische Galois-Gruppe ist eine "pro-algebraische" Erweiterung von $G_\mathbb{Q}$.

### 8.4 Verworfene Lösungsansätze (mit Begründung)

**Hodge-Theorie über Zahlkörpern**: Hodge-Theorie liefert die Standardvermutungen über $\mathbb{C}$, aber über Zahlkörpern fehlt das analytische Analogon. $p$-adische Hodge-Theorie (Fontaine) liefert Vergleichsisomorphismen, aber keine Positivität.

**$\ell$-adische Methoden**: Delignes Weil-II beweist Reinheit für $\ell$-adische Garben, aber die Brücke zu algebraischen Zykeln (Standardvermutung $A$) ist das Kernproblem und bleibt offen.

**Direkte Konstruktion von $G_{\mathrm{mot}}$**: Versuche, $G_{\mathrm{mot}}$ als projektives System algebraischer Gruppen zu konstruieren (Langlands, Milne), scheitern daran, dass die Morphismen-Sets nicht bekannt sind (Standardvermutung $C$).

### 8.5 Aktuelle Entwicklungen (2020–2024)

**Ayoub (2022–2024)**: Arbeit an motivischen $\infty$-Kategorien und deren Verbindung zu Tannaka-Theorie liefert neue Perspektiven. **Scholze (2019–2024)**: Perfektoide Räume und "prismatische Kohomologie" als neue universelle $p$-adische Kohomologietheorie. Diese könnte neue Vergleichsisomorphismen liefern und die Lücke zwischen $p$-adischen und motivischen Strukturen schließen.

**Clausen-Scholze (2024)**: "Animated condensed mathematics" liefert neue Kategorien-theoretische Werkzeuge für globale motivische Fragen. Erste Anwendungen auf Standardvermutungen sind in Vorbereitung.

### 8.6 Offene Teilprobleme

- Gilt Standardvermutung $C$ für eine Klasse von Varietäten jenseits von Kurven und abelschen Varietäten?
- Ist die Kategorie der Nori-Motive äquivalent zur Voevodsky-Motiven auf der abelianen Ebene?
- Kann man die motivische Galois-Gruppe $G_{\mathrm{mot}}$ explizit in einem nichttrivialen Fall berechnen?
- Gilt die Fontaine-Mazur-Vermutung "motivisch" — d. h. kommen alle geometrischen Galois-Darstellungen aus Motiven?

---

## Zusammenfassung und Querverbindungen

Die acht hier behandelten Probleme sind trotz ihrer scheinbaren Verschiedenheit durch tiefe Querverbindungen verknüpft:

| Problem | Zentrale Methode | Verbindung |
|---|---|---|
| Jacobi-Vermutung | Polynomautomorphismen, Weyl-Algebra | Dixmier-Vermutung (äquivalent) |
| Hadwiger-Vermutung | Graphen-Minor-Theorie | Robertson-Seymour, 4-Farben-Satz |
| Andrews-Curtis | Gruppentheorie, 2-Komplexe | 4-dimensionale Topologie |
| Beal-Vermutung | Elliptische Kurven, FLT | ABC-Vermutung (impliziert Beal) |
| Fontaine-Mazur ($n \geq 3$) | $p$-adische Hodge-Theorie | Langlands-Programm |
| Rang-Unbeschränktheit | Selmer-Gruppen, BSD | Goldfeld-Vermutung |
| Lehmer-Mahler | Salem-Zahlen, Arithmetik | Bogomolov-Problem |
| Motivische Galois | Tannaka-Kategorien | Standardvermutungen, Langlands |

Die **Fontaine-Mazur-Vermutung**, der **Rang elliptischer Kurven**, die **motivische Galois-Gruppe** und das **Lehmer-Mahler-Maß** sind alle durch das Langlands-Programm und $p$-adische Methoden verbunden. Die **Beal-Vermutung** hängt mit der **ABC-Vermutung** zusammen (selbst eines der tiefsten offenen Probleme). Die **Hadwiger-Vermutung** und **Andrews-Curtis** stehen etwas isolierter, verbinden aber Kombinatorik und Topologie in fundamentaler Weise.

Ein gemeinsames Thema ist die **Insuffizienz bekannter Methoden**: Für alle acht Probleme haben die "natürlichen" Beweismethoden klare Grenzen, und ein Durchbruch würde neue mathematische Ideen erfordern — möglicherweise aus ganz anderen Bereichen der Mathematik.

---

## Literaturhinweise (Auswahl)

- **Bass, Connell, Wright** (1982): *The Jacobian conjecture: reduction of degree and formal expansion of the inverse*. Bull. AMS 7.
- **Hadwiger** (1943): *Über eine Klassifikation der Streckenkomplexe*. Vierteljahrsschr. Naturforsch. Ges. Zürich 88.
- **Robertson, Seymour, Thomas** (1993): *Hadwiger's conjecture for $K_6$-free graphs*. Combinatorica 13.
- **Andrews, Curtis** (1965): *Free groups and handlebodies*. Proc. AMS 16.
- **Darmon, Granville** (1995): *On the equations $z^m = F(x,y)$ and $Ax^p + By^q = Cz^r$*. Bull. London Math. Soc. 27.
- **Fontaine, Mazur** (1995): *Geometric Galois representations*. In: Elliptic curves, modular forms, Fermat's last theorem.
- **Kisin** (2009): *Moduli of finite flat group schemes, and modularity*. Ann. Math. 170.
- **Bhargava, Shankar** (2015): *Ternary cubic forms having bounded invariants, and the existence of a positive proportion of elliptic curves having rank 0*. Ann. Math. 181.
- **Lehmer** (1933): *Factorization of certain cyclotomic functions*. Ann. Math. 34.
- **Boyd** (1981): *Speculations concerning the range of Mahler's measure*. Canad. Math. Bull. 24.
- **Grothendieck** (1968): *Standard conjectures on algebraic cycles*. In: Algebraic Geometry.
- **Jannsen** (1992): *Motives, numerical equivalence, and semi-simplicity*. Invent. Math. 107.
- **Park, Poonen, Voight, Wood** (2019): *A heuristic for boundedness of ranks of elliptic curves*. J. Eur. Math. Soc. 21.
- **Norin** (2019): *A new proof of the graph minor theorem*. J. Combin. Theory Ser. B.
