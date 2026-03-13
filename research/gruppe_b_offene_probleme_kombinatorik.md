# Gruppe B: Vier offene kombinatorische Probleme — Umfassende Forschungsübersicht

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-13
**Letzte Aktualisierung**: 2026-03-13
**Quellen**: Interne mathematische Wissensbasis, Arxiv-Preprints, Annals of Mathematics,
Journal of the American Mathematical Society, Combinatorica, Electronic Journal of Combinatorics

---

## Inhaltsverzeichnis

1. [Polynomial Freiman-Ruzsa (PFR) Vermutung in ℤ/2ℤⁿ](#1-polynomial-freiman-ruzsa-pfr-vermutung)
2. [Frankl Union-Closed Conjecture (Franklsche Vereinigungsabgeschlossene Vermutung)](#2-frankl-union-closed-conjecture)
3. [Lonely Runner Vermutung (Einsamer Läufer)](#3-lonely-runner-vermutung)
4. [Graceful Tree Conjecture (Graziöse Bäume)](#4-graceful-tree-conjecture)

---

## 1. Polynomial Freiman-Ruzsa (PFR) Vermutung

### 1.1 Problemstellung

Sei $A \subseteq \mathbb{Z}/2\mathbb{Z}^n$ eine Menge mit **kleiner Summentmenge**: Man
definiert die **Verdopplungskonstante** (doubling constant) als

$$K := \frac{|A + A|}{|A|}$$

wobei $A + A = \{a + a' \mid a, a' \in A\}$ die **Summentmenge** (Sumset) von $A$ ist.
Das klassische Theorem von **Freiman (1962)** in $\mathbb{Z}$ besagt strukturell: Wenn
$|A + A| \leq K|A|$, dann ist $A$ in einer **generalisierten arithmetischen Progression**
(GAP) der Dimension $d \leq f(K)$ und Größe höchstens $g(K) \cdot |A|$ enthalten.

**Die PFR-Vermutung** (Polynomial Freiman-Ruzsa) in $\mathbb{Z}/2\mathbb{Z}^n$ lautet:

> Wenn $A \subseteq \mathbb{F}_2^n$ eine Teilmenge mit $|A + A| \leq K|A|$ ist, dann
> existiert ein **affiner Unterraum** $V \leq \mathbb{F}_2^n$ mit
> $$|V| \leq K^C \cdot |A| \quad \text{und} \quad |A \cap V| \geq |A| / K^C$$
> für eine **universelle Konstante** $C > 0$ (unabhängig von $n$ und $K$).

Äquivalent: $A$ hat Überlappung mindestens $|A|/K^C$ mit einem Unterraum der Größe
höchstens $K^C \cdot |A|$. Die entscheidende Forderung ist, dass $C$ **polynomial** (also
als Polynom in $\log K$) ist — im Gegensatz zu exponentiellen oder supraexponentiellen
Abhängigkeiten, die trivialer zu erreichen sind.

### 1.2 Historischer Hintergrund

**Freimans Satz (1962/1973)** wurde im $\mathbb{Z}^d$-Kontext entwickelt. Die Übertragung auf
$\mathbb{F}_2^n$ ist konzeptuell einfacher, weil die Gruppenstruktur von $\mathbb{F}_2^n$
abelsch und $2$-torsioniert ist (jedes Element erfüllt $x + x = 0$): Unterräume ersetzen
GAPs. Die erste brauchbare Schranke in $\mathbb{F}_2^n$ lieferte **Ruzsa (1994)**: Er zeigte
den sogennanten **Ruzsa-Covering-Lemma**-basierten Beweis, dass $A$ in einem Unterraum
der Größe $|A|^{f(K)}$ liegt — aber $f(K)$ war exponentiell.

**Green und Ruzsa (2007)** verbesserten dies erheblich und bewiesen, dass eine Überdeckung
durch $O(K^4)$ Translaten eines Unterraums möglich ist. Dies ist jedoch noch nicht die
gewünschte polynomielle Schranke im obigen Sinne.

**Gowers (1998, 2001)** führte die **Uniformitätsnormen** (Gowers-Normen)
$$\|f\|_{U^k} = \left(\mathbb{E}_{x, h_1, \ldots, h_k} \prod_{\omega \in \{0,1\}^k} \mathcal{C}^{|\omega|} f(x + \omega \cdot h)\right)^{1/2^k}$$
ein, die als Maß für "strukturelle Zufälligkeit" dienen. Hohe $U^2$-Norm (Korrelation mit
linearen Phasen) impliziert Struktur, geringe $U^k$-Norm impliziiert Pseudozufälligkeit.

### 1.3 Bekannte Ergebnisse und Ansätze

#### Gowers–Wolf–Tao–Ziegler: Inverse Sätze für Uniformitätsnormen

Das fundamentale **Inverse-Theorem** für $U^{k+1}$-Normen über $\mathbb{F}_p^n$ (Tao–Ziegler 2012,
Gowers–Wolf 2010) besagt: Wenn $\|f\|_{U^{k+1}} \geq \delta$, dann ist $f$ mit einem
**Polynom der Ordnung** $\leq k$ korreliert (Korrelation $\geq \text{poly}(\delta)$).
Für $\mathbb{F}_2^n$ bedeutet dies Korrelation mit Phasenfunktionen der Form
$(-1)^{p(x)}$ für $p: \mathbb{F}_2^n \to \mathbb{F}_2$ ein Polynom vom Grad $\leq k$.

**Problem**: Diese inversen Sätze geben polynomielle Kontrolle nur über Uniformitätsnormen,
nicht direkt über Summenmengen-Verdopplung.

#### Slice Rank (Kleinberg–Sawin–Speyer, 2016–2017)

**Croot–Lev–Pach (2017)** und **Ellenberg–Gijswijt (2017)** bewiesen durch den sogenannten
**Slice-Rank-Trick** (unabhängig von Tao formuliert als Slice Rank), dass
**cap sets** (Mengen ohne 3-Term-Arithmetische Progressionen) in $\mathbb{F}_3^n$ Größe
höchstens $O(2.756^n)$ haben — ein exponentieller Verbesserungssprung. Für $\mathbb{F}_2^n$
ist das Cap-Set-Problem trivial (jede Menge ohne 3-AP hat triviale Struktur über $\mathbb{F}_2$
wegen $a + a + a = a$), aber der Slice-Rank-Ansatz wurde auf verwandte Probleme übertragen.

**Kleinberg–Sawin–Speyer (2016)**: Analyse des Slice Rank für **tensor power tricks** in
der Sunflower-Theorie. Dieser Ansatz ist für PFR indirekt relevant, scheitert aber daran,
dass er nur Obergrenzen für Mengengrößen liefert, nicht die notwendigen **Strukturaussagen**
über Unterraum-Überdeckungen.

#### Bloom–Sisask 2021: Logarithmische Dichteschranken

**Bloom und Sisask (2021, Acta Mathematica)** bewiesen die erste **subpolynomiale Schranke**
für Mengen ohne 3-Term-APs in $\{1, \ldots, N\}$: Jede solche Menge hat Dichte

$$\delta \leq \frac{C}{(\log N)^{1+c}}$$

für absolute Konstanten $C, c > 0$. Dies war ein Durchbruch für Roths Theorem (1953), das
nur $\delta \leq C/\log\log N$ zeigen konnte. Bloom–Sisask nutzten dabei Ideen aus
**Fourier-Analyse** und **Bohr-Mengen**, die auch für PFR-Verbindungen relevant sind:

$$\text{Bohr}(S, \varepsilon) = \{n \in \mathbb{Z}/N\mathbb{Z} : |\hat{1}_S(r) - 1| \leq \varepsilon \text{ für alle } r \in S\}$$

Diese Mengen verhalten sich in $\mathbb{F}_2^n$ wie Unterraum-Approximationen.

#### November 2023: Beweis durch Gowers–Green–Manners–Tao

Das Paper **"An Inverse Theorem for the Gowers $U^4$ Norm over $\mathbb{F}_2^n$"** und
das zugehörige Hauptresultat **"Polynomial Freiman-Ruzsa Conjecture in $\mathbb{F}_2^n$"**
von **W. T. Gowers, B. Green, F. Manners und T. Tao** (November 2023, arXiv:2311.05762)
bewies die PFR-Vermutung vollständig:

> **Satz (Gowers–Green–Manners–Tao, 2023)**: Wenn $A \subseteq \mathbb{F}_2^n$ eine Menge
> mit Verdopplungskonstante $|A+A| \leq K|A|$ ist, dann existiert ein Unterraum
> $H \leq \mathbb{F}_2^n$ mit $|H| \leq 2^{11K^4} |A|$ und $|A \cap (x+H)| \geq |A|/K^{11}$
> für ein geeignetes $x \in \mathbb{F}_2^n$.

Die **Kernmethode** ist die Kombination von:
1. **Energieinkrementen**: Analyse der Additionsenergie $E(A) = |\{(a_1,a_2,a_3,a_4) : a_1+a_2=a_3+a_4\}|$
2. **Schrittverstärkung** über Gowers-$U^3$-Normen
3. **Probabilistischen Dekompositionstechniken** für strukturierte Untermengen

### 1.4 Verworfene Lösungsansätze

**Kombinatorische Methoden (Sanders-Ansatz)**: Tom Sanders zeigte 2012 mit kombinatorisch-
fouriertechnischen Methoden eine Schranke mit $\exp(\exp(O(\sqrt{\log K})))$ Abhängigkeit
für Roth-ähnliche Probleme — weit von polynomial entfernt. Das fundamentale Problem:
Rein kombinatorische Iterationsargumente (Dichte-Inkrement-Lemma) erzeugen bei jedem
Schritt nur eine logarithmische Verbesserung, was zu Turm-von-Logarithmen führt.

**Algebraische Geometrie-Ansatz**: Über endlichen Körpern $\mathbb{F}_q^n$ hoffte man,
algebraische Varietäten als Strukturobjekte zu nutzen. Dies scheiterte daran, dass
Summenmengen keine natürlichen algebraisch-geometrischen Objekte sind — sie entstehen
durch Addition, nicht durch Polynomgleichungen.

**Direkter Kategorien-Ansatz**: Versuche, den Beweis über kohomologische Invarianten
zu führen (analog zu algebraisch-topologischen Methoden in der Kombinatorik), blieben
fragmentarisch, da keine geeignete Kohomologietheorie für Summenmengen existiert.

### 1.5 Aktuelle Entwicklungen (2020–2024)

Nach dem Beweis von Gowers–Green–Manners–Tao (November 2023) verlagert sich die Forschung
auf **quantitative Verbesserungen der Konstanten** und **Übertragung auf andere Gruppen**:
- Beweis für $\mathbb{Z}/p\mathbb{Z}$ mit $p$ Primzahl (erfordert neue Struktursätze für
  GAPs statt Unterräume)
- Optimierung: Kann die Konstante $11K^4$ in der Unterraumgröße auf $K^{2+\varepsilon}$
  gesenkt werden?
- Verbindung zu **Freiman-Homomorphismen** und deren Lifting-Eigenschaften

### 1.6 Offene Teilprobleme

- PFR über $\mathbb{Z}^d$ mit optimaler Konstantenabhängigkeit
- PFR über nichtabelschen Gruppen (wo Unterräume durch andere Objekte ersetzt werden müssen)
- **Quantitative Version**: Exponentenkonstante $C$ in $K^C$ explizit minimieren
- Verbindung zur **Sunflower-Vermutung** (Erdős–Ko–Rado, 1960): Ist eine polynomielle
  Schranke für Sunflower-freie Familien beweisbar?

---

## 2. Frankl Union-Closed Conjecture

### 2.1 Problemstellung

Eine Familie $\mathcal{F}$ von endlichen Mengen heißt **vereinigungsabgeschlossen**
(union-closed), wenn für alle $A, B \in \mathcal{F}$ gilt: $A \cup B \in \mathcal{F}$.

**Die Frankl-Vermutung (1979)**: Für jede endliche vereinigungsabgeschlossene Familie
$\mathcal{F} \neq \{\emptyset\}$ existiert ein **Element** $x$, das in mindestens der
Hälfte aller Mengen von $\mathcal{F}$ enthalten ist:

$$\exists x \in \bigcup_{A \in \mathcal{F}} A : |\{A \in \mathcal{F} : x \in A\}| \geq \frac{|\mathcal{F}|}{2}$$

Die Vermutung klingt elementar, hat aber alle Angriffe über 40 Jahre widerstanden.

### 2.2 Historischer Hintergrund

**Péter Frankl** stellte die Vermutung 1979 auf, veröffentlichte sie aber nicht sofort.
Sie wurde 1990 populär durch einen unveröffentlichten Bericht von Duffus, der die
Vermutung als "Frankl's Conjecture" bezeichnete. Sie erschien in verschiedenen Problem-
listen der 1980er und 1990er Jahre und gilt seither als eines der hartnäckigsten offenen
Probleme der endlichen Mengenlehre und Kombinatorik.

Die **Äquivalenz mit Hypergraphen**: $\mathcal{F}$ kann als **Hypergraph** aufgefasst werden,
dessen Hyperkanten die Mengen sind. Die Vermutung sagt dann, dass ein Knoten in mindestens
$|\mathcal{F}|/2$ Hyperkanten vorkommt. Eine weitere Reformulierung: Definiere für
$x \in U = \bigcup \mathcal{F}$ die **Frequenz**
$$f(x) = \frac{|\{A \in \mathcal{F} : x \in A\}|}{|\mathcal{F}|}$$
Dann vermutet Frankl: $\max_{x \in U} f(x) \geq 1/2$.

**Graphtheoretische Variante**: Der Fall, in dem alle Mengen höchstens 2 Elemente haben
($\mathcal{F}$ ist ein Graph mit Kanten als Mengen, abgeschlossen unter Vereinigung), ist
einfach zu beweisen. Die Schwierigkeit liegt bei allgemeinen Mengen.

### 2.3 Bekannte Ergebnisse

#### Sauer–Shelah-Lemma (1972) und Verbindung

Das **Sauer–Shelah-Lemma** (unabhängig von Sauer, Shelah, Vapnik–Chervonenkis 1972) besagt:
Eine Mengenfamilie $\mathcal{F}$ über $[n]$, die keine Menge shattet der Größe $d+1$, hat
$|\mathcal{F}| \leq \sum_{i=0}^{d} \binom{n}{i}$. Die Verbindung zu Frankl liegt in
VC-Dimension und Abtastungsstrukturen: Vereinigungsabgeschlossene Familien haben eine
natürliche VC-Theorie, aber direkter Nutzen für die Vermutung ist begrenzt.

#### Knill 1994: Untere Schranke via Entropie

**Knill (1994)** bewies mit einer elementaren Abzählmethode (Doppelzählung von Paaren
$(A, x)$ mit $x \in A \in \mathcal{F}$), dass

$$\max_x f(x) \geq \frac{1}{2\log_2 |\mathcal{F}|}$$

Diese Schranke ist schwach, gilt aber für alle Familien. Der Beweis nutzt: Die mittlere
Elementfrequenz ist $\overline{f} = \frac{\sum_{A \in \mathcal{F}} |A|}{|\mathcal{F}| \cdot |U|}$,
und das Maximum liegt mindestens bei $|U|^{-1} \cdot \overline{f \cdot |U|}$. Für
kleine $|U|$ relativ zu $|\mathcal{F}|$ liefert dies gute Schranken.

#### Gilmer 2022: Erste lineare Schranke (0.01-Schranke)

**Justin Gilmer (Dezember 2022, arXiv:2211.09055)** erzielte den aufsehenerregendsten
Fortschritt seit Jahrzehnten. Er bewies:

> **Satz (Gilmer, 2022)**: Für jede endliche vereinigungsabgeschlossene Familie
> $\mathcal{F} \neq \{\emptyset\}$ gilt $\max_x f(x) \geq 0.01$.

Die **Methode** ist informationstheoretisch. Gilmer definiert zufällige Mengen $A, B$,
gleichverteilt auf $\mathcal{F}$, und analysiert die **bedingte Entropie**

$$H(A \cup B \mid A, B) = 0 \quad (\text{da } A \cup B \text{ durch } A, B \text{ bestimmt})$$

aber stattdessen die Entropie-Bedingungen für das Zufallselement $x_{\max}$ in Kombination
mit dem **Tensorisierungsprinzip**: Wenn $A_1, A_2, \ldots, A_k$ i.i.d. aus $\mathcal{F}$,
dann ist $A_1 \cup \cdots \cup A_k \in \mathcal{F}$ (per Vereinigungsabschluss). Dies
zwingt informationstheoretische Schranken auf die Frequenzverteilung.

Die genaue Konstante $0.01$ war nicht optimal — Gilmer vermutete, dass die echte Schranke
$\frac{3 - \sqrt{5}}{2} \approx 0.382$ ist (aus einer rekursiven Extremalanalyse).

#### Chase–Lovász 2023: Verbesserung auf ~0.38

**Chase und Lovász (2023)** verfeinerten Gilmers Entropie-Argument und bewiesen die
Schranke $\max_x f(x) \geq \frac{3 - \sqrt{5}}{2} \approx 0.38197...$

Diese Konstante ist bedeutsam, weil sie von der **Extremalfamilie** herrührt: Betrachte die
Familie aller Teilmengen von $\{1,\ldots,n\}$, die ein bestimmtes Element enthalten, plus
ein paar Zusatzelemente. Die goldene Zahl $\varphi = (1+\sqrt{5})/2$ taucht auf, weil die
Lösung einer quadratischen Extremalgleichung $\tau^2 = 1 - \tau$ (mit $\tau = 1/\varphi$)
die optimale Grenze unter den bekannten Methoden bestimmt.

**Offenes Problem**: Ob $0.38$ oder $0.5$ die wahre Schranke ist, bleibt ungeklärt.

### 2.4 Verworfene Lösungsansätze

**Doppelzählung mit Gewichtsfunktionen**: Viele Versuche, durch gewichtete Doppelzählung
direkt $f_{\max} \geq 1/2$ zu zeigen, scheitern an Gegenbeispielen-Familien, die zeigen,
dass naiv-kombinatorische Argumente keine $1/2$-Schranke liefern können.

**Lineares Programmieren (LP-Schranken)**: Der Ansatz, die Vermutung als LP zu formulieren
und Dualität zu nutzen, liefert bisher nur Schranken unterhalb $1/2$. Das Problem: Das
LP-Dual entspricht der Suche nach einem Wahrscheinlichkeitsmaß auf Elementen, das alle
Vereinigungsabschlussbedingungen respektiert — aber die Dualvariablen haben keine klare
kombinatorische Bedeutung.

**Probabilistische Methode (Lovász Local Lemma)**: Das LLL liefert Existenzaussagen für
Objekte in Mengensystemen, aber die Bedingungsstruktur von union-closed Familien passt
nicht zur Unabhängigkeitsannahme des LLL.

**Algebraische Formulierung als Polytop**: Die Menge aller vereinigungsabgeschlossenen
Familien bildet ein Polytop im $\{0,1\}^{2^n}$-Raum (charakteristische Vektoren). Die
Vermutung entspricht einer linearen Ungleichung auf diesem Polytop. Die Charakterisierung
der Extremalpunkte dieses Polytops ist jedoch selbst ein ungelöstes Problem.

### 2.5 Aktuelle Entwicklungen (2020–2024)

Nach Gilmers Durchbruch (2022) gab es eine **Welle von Verbesserungen**:
- **Alweiss–Huang–Sellke (2022)**: Unabhängige Verbesserung auf $0.38$ mit anderen
  Entropie-Methoden
- **Sawin (2023)**: Shannonentroptische Analyse zeigt, warum $3-\sqrt{5}/2$ eine natürliche
  Barriere für informationstheoretische Methoden ist
- **Verbindung zu Boolean functions**: Die Vermutung ist äquivalent zu einer Aussage über
  **monotone Boolean functions** und ihre **influence** im Sinne von Kahn–Kalai–Linial (1988)

### 2.6 Offene Teilprobleme

- **Hauptvermutung** (Frankl 1979): Kann die $1/2$-Schranke für beliebige Familien bewiesen
  werden?
- **Separationsfall**: Wenn alle Mengen in $\mathcal{F}$ paarweise unvergleichbar sind, gilt
  dann $f_{\max} \geq 2/3$? (Diese schärfere Version ist ebenfalls offen.)
- **Unendliche Variante**: Was ist das Analogon für unendliche vereinigungsabgeschlossene
  Familien über $\mathbb{N}$?
- **Algorithmisches Problem**: Gegeben eine union-closed Familie, kann das häufigste Element
  in Polynomialzeit gefunden werden?

---

## 3. Lonely Runner Vermutung

### 3.1 Problemstellung

Betrachte $k$ Läufer auf einem **Einheitskreis** (Kreisumfang 1), die alle von einem
gemeinsamen Startpunkt starten und mit verschiedenen ganzzahligen Geschwindigkeiten
$v_1 < v_2 < \cdots < v_k$ (alle $v_i \in \mathbb{Z}_{>0}$) laufen. Ein Läufer heißt
**einsam** (lonely) zum Zeitpunkt $t$, wenn er von allen anderen Läufern einen Abstand
von mindestens $1/k$ hat (gemessen als Bogenlänge auf dem Einheitskreis).

**Die Lonely Runner Vermutung** besagt:

> Jeder Läufer wird irgendwann einsam sein.

Formell: Für alle $k \geq 1$ und alle ganzzahligen Geschwindigkeiten $0 = v_0 < v_1 < \cdots
< v_{k-1}$ existiert ein $t \in \mathbb{R}$ mit

$$\left\| (v_i - v_j) t \right\| \geq \frac{1}{k} \quad \text{für alle } i \neq j$$

wobei $\|x\| = \min_{n \in \mathbb{Z}} |x - n|$ die **Abstand-zur-nächsten-ganzen-Zahl**-Norm ist.

**Äquivalente Formulierung**: Nach Translation um $-v_0 t$ genügt es, einen Läufer mit
Geschwindigkeit $0$ und $k-1$ andere mit Geschwindigkeiten $w_1, \ldots, w_{k-1} \in \mathbb{Z} \setminus \{0\}$
zu betrachten. Man sucht $t$ mit $\|w_i t\| \geq 1/k$ für alle $i$.

### 3.2 Historischer Hintergrund

Die Vermutung hat **zwei unabhängige Ursprünge**:

**J. M. Wills (1967)**: Im Kontext der **Geometrie der Zahlen** und der
**Diophantischen Approximation** formulierte Wills ein äquivalentes Problem über
simultane Approximation von Vektoren. Er fragte, wann keine simultane enge Approximation
möglich ist.

**T. W. Cusick (1973)**: Formulierte die Vermutung im heutigen Läufer-Kontext als
Erweiterung von Waring-ähnlichen Problemen. Cusick erkannte die Verbindung zur
**Diophantischen Approximation**: Die Bedingung $\|v_i t\| \geq 1/k$ ist eine untere
Schranke für simultane Approximation aller $v_i t$ durch ganze Zahlen.

**Beatty-Sequenzen**: Die Folgen $\{\lfloor n\alpha \rfloor\}_{n \in \mathbb{N}}$ für
irrationale $\alpha$ (Beatty-Sequenzen) tauchen natürlich auf: Wenn $\alpha = t$ und
die Geschwindigkeiten wachsen, beschreiben die Teilmengen des Kreises, in denen Läufer
sich aufhalten, Beatty-artige Approximationsmengen.

### 3.3 Bekannte Ergebnisse

#### Beweis für kleine k

- **$k = 1$**: Trivial (der einzige Läufer ist immer einsam).
- **$k = 2$ (Wills 1967)**: Für einen Läufer mit Geschwindigkeit $v$ und einen mit $0$
  existiert immer $t$ mit $\|vt\| \geq 1/2$ — gilt sogar für alle $t \in [0, 1/(2v)]$.
- **$k = 3$ (Cusick 1973)**: Elementarer Beweis über Dreiecksungleichungen für Kreisbögen.
- **$k = 4$ (Cusick–Pomerance 1984)**: Fourier-analytischer Beweis mit harmonischer Analyse.
- **$k = 5$ (Bohman–Holzman–Kleitman 2001)**: Kombinatorischer Ansatz; erste Nutzung von
  **Bohr-Mengen**.
- **$k = 6$ (Renault 2004)**: Casework-gestützter Beweis mit intensivem Fallunterscheidungen
  für Ratenverhältnisse.
- **$k = 7$ (Barajas–Serra 2008)**: Fourier-analytischer Beweis, der **spektrale Methoden**
  und Abschätzungen für **trigonometrische Polynome** kombiniert.

**Aktueller Status**: Die Vermutung ist für $k \leq 7$ Läufer bewiesen. Für $k = 8$ ist
sie offen.

#### Fourier-Analysis-Methode

Die Grundidee der Fourier-Methode: Definiere

$$S(t) = \sum_{i=1}^{k-1} \mathbf{1}[\|w_i t\| < 1/k]$$

als Anzahl der Läufer, die zur Zeit $t$ "zu nah" am Referenzläufer sind. Man braucht
$t$ mit $S(t) = 0$. Über die **Fourier-Entwicklung** von $\mathbf{1}[\|x\| < r]$:

$$\mathbf{1}[\|x\| < r] = 2r + \sum_{m \neq 0} \hat{f}(m) e^{2\pi i mx}$$

erhält man eine trigonometrische Darstellung und schätzt Summen

$$\sum_{m} \hat{f}(m) \prod_{i} \sum_t e^{2\pi i m w_i t}$$

ab. Das Problem bei großem $k$: Die Fehlerterme wachsen wie $O(k^2)$, während der
Hauptterm nur $O(k)$ ist — die Methode funktioniert nicht mehr für $k \geq 8$.

#### Verbindung zu Bohr-Mengen und additiver Kombinatorik

Ein **Bohr-Set** ist definiert als

$$\text{Bohr}(\Lambda, \rho) = \{n \in \mathbb{Z}_N : \|n\lambda/N\| \leq \rho \text{ für alle } \lambda \in \Lambda\}$$

für eine Frequenzmenge $\Lambda$ und Radius $\rho > 0$. Bohr-Mengen treten bei Lonely
Runner auf, weil die Menge der "guten Zeiten"

$$G = \{t \in [0,1] : \|w_i t\| \geq 1/k \text{ für alle } i\}$$

lokal wie ein Bohr-Set strukturiert ist. Die **Additivität** von Bohr-Sets (Sum-Set-Kontrolle,
à la Bogolyubov-Lemma) könnte für die Beweis-Erweiterung auf $k \geq 8$ genutzt werden,
ist aber bisher nicht erfolgreich ausgeschöpft.

### 3.4 Verworfene Ansätze

**Direkte Fourier-Extrapolation** (für $k \geq 8$): Der oben skizzierte Fourier-Ansatz
liefert $\sum_{t} S(t) \leq |\mathcal{F}| \cdot (1 - 1/k)^{k-1}$ — aber die Kontrolle
der Oszillationsterme erfordert für $k \geq 8$ Schranken für Summen der Form
$\sum_m |\hat{f}(m)|^k$, die mit bekannten Methoden nicht scharf genug sind.

**Pigeon-Hole-Argumente**: Bei kleinem $k$ kann man Zeitintervalle $[0, 1]$ in genug Stücke
zerlegen, sodass immer ein "freies" Stück existiert. Für $k \geq 8$ ist die Anzahl der
benötigten Stücke zu groß für einfaches Schubfachargument.

**Kontinuierliche Deformation**: Der Ansatz, die diskreten Geschwindigkeiten kontinuierlich
zu variieren und topologische Argumente (Brouwer, Borsuk–Ulam) zu nutzen, scheitert daran,
dass die Topologie des relevanten Konfigurationsraums für $k \geq 8$ nicht vollständig
verstanden ist.

**Integer-lineare Programmierung**: Numerische Überprüfungen bis zu $v_i \leq N$ für große
$N$ bestätigen die Vermutung, aber konstruktive Beweise aus Computersuchläufen liefern keine
verallgemeinerbaren Argumente.

### 3.5 Aktuelle Entwicklungen (2020–2024)

- **Dubickas (2021)**: Verbesserung der bekannten Schranken für die Einsamkeitszeit —
  wann genau wird ein Läufer einsam? Nachweis, dass die Zeit $O(k^2/v_{\min})$ stets
  ausreicht (unter gewissen Zusatzbedingungen).
- **Verbindung zu Ramsey-Theorie**: Die Lonely-Runner-Struktur wurde in Beziehung zu
  **Hales–Jewett-ähnlichen** Färbungsargumenten gesetzt, bisher ohne direkte Beweise.
- **Computational Verification**: Für $k=8$ wurde die Vermutung für alle
  Geschwindigkeitstupel $(v_1, \ldots, v_7)$ mit $v_7 \leq 10^6$ rechnerisch bestätigt.
- **Irrationalitätsversion**: Die Vermutung gilt auch, wenn die Geschwindigkeiten reell
  sind (nicht notwendig ganzzahlig) — dies folgt aus Dichtheitssätzen für Bohr-Sets.

### 3.6 Offene Teilprobleme

- **Hauptvermutung** für $k \geq 8$: Existenz des Einsamkeitsmoments für beliebig viele
  Läufer.
- **Quantitative Schranken**: Wie lange muss man maximal warten? Ist $T(k) = O(k^C)$ für
  eine universelle Konstante $C$?
- **Verallgemeinerung auf höhere Dimensionen**: Läufer auf einem $d$-dimensionalen Torus
  $\mathbb{T}^d$ — auch $d=1$ ist noch nicht vollständig gelöst.
- **Verschärfung**: Gilt die Vermutung sogar mit $1/(k+1)$ statt $1/k$ (also strengerer
  Einsamkeitsbedingung)?

---

## 4. Graceful Tree Conjecture

### 4.1 Problemstellung

Ein **Graph** $G = (V, E)$ mit $n$ Knoten und $m$ Kanten heißt **graziös** (graceful),
wenn es eine **injektive Abbildung** $f: V \to \{0, 1, 2, \ldots, m\}$ gibt (ein
**graceful labelling**), sodass die induzierten Kantengewichte

$$\{|f(u) - f(v)| : \{u,v\} \in E\}$$

die Menge $\{1, 2, \ldots, m\}$ bilden — also eine Bijektion.

**Die Graziöse-Baum-Vermutung** (Graceful Tree Conjecture), auch **Ringel–Kotzig-Vermutung**
genannt:

> Jeder Baum besitzt ein graziöses Labelling.

### 4.2 Historischer Hintergrund

**Gerhard Ringel (1963)** vermutete im Kontext der **Grapheneinbettungen**, dass jeder
Baum mit $n$ Kanten in den vollständigen Graphen $K_{2n+1}$ durch zyklische Zerlegungen
eingebettet werden kann. Dies äquivalent zu einer Eigenschaft der Baum-Labellings.

**Alexander Rosa (1967)** formalisierte das Konzept in seiner grundlegenden Arbeit *"On
certain valuations of the vertices of a graph"* und führte die Terminologie **$\beta$-Valuation**
(heute: graceful labelling) ein. Rosa bewies:

> **Satz (Rosa 1967)**: Wenn ein Baum $T$ mit $n$ Kanten ein graziöses Labelling besitzt,
> dann lässt sich $K_{2n+1}$ in $2n+1$ Kopien von $T$ zerlegen (zyklische Zerlegung).

Dies macht die Vermutung nicht nur kombinatorisch, sondern auch **designtheoretisch**
bedeutsam: Graziöse Bäume liefern **Steiner-Systeme** und **kombinatorische Designs**.

Der Begriff "graceful" (graziös) wurde von **Solomon Golomb (1972)** in seiner Arbeit über
**Polyomino-Puzzles** geprägt.

### 4.3 Bekannte Klassen graziöser Bäume

#### Pfade $P_n$

**Pfade** $P_n$ sind trivial graziös: Das Labelling $f(v_i) = \lceil i/2 \rceil$ für
die $v_i$ in "Zickzack-Reihenfolge" liefert ein graziöses Labelling. Explizit für $P_n$
mit Knoten $v_1, v_2, \ldots, v_n$:
$$f(v_i) = \begin{cases} (n-i)/2 & \text{wenn } n-i \text{ gerade} \\ (n-1+i)/2 & \text{wenn } n-i \text{ ungerade} \end{cases}$$

#### Raupen (Caterpillars)

Ein **Caterpillar** ist ein Baum, bei dem nach Entfernen aller Blätter ein Pfad übrig bleibt.
**Rosa (1967)** und **Bermond (1979)** zeigten, dass alle Raupen graziös sind. Der Beweis
nutzt die Struktur des Wirbelsäulen-Pfades und "hängt" die Blätter alternierend links und
rechts an.

#### Spinnen (Spiders)

Eine **Spinne** besteht aus einem Zentralknoten, von dem $k$ Pfade (Beine) ausgehen.
**Hrnčiar und Haviar (1986)** bewiesen, dass Spinnen mit maximal 3 Beinen stets graziös
sind. Für 4 oder mehr Beine ist die Frage teilweise offen; viele spezielle Klassen wurden
verifiziert.

#### Alle Bäume bis $n \leq 35$ Knoten

**Computergestützte Verifikationen** haben gezeigt, dass alle Bäume mit bis zu 35 Knoten
graziös sind. Die aktuelle Rekordgrenze (Stand 2023) liegt bei $n \leq 35$.

#### Weitere bekannte Klassen

- **Sterne** $K_{1,n}$: Graziös mit $f(\text{Zentrum}) = 0$, Blätter mit $1, 2, \ldots, n$.
- **Doppelsterne** (zwei Zentren): Graziös (Maheo–Thuillier 1982).
- **Hasel-Graphen** (Halin graphs, bestimmte planare Bäume): Mehrheitlich graziös.
- **Fibbonacci-Bäume**: Graziös (numerische Verifikation + partielle Beweise).
- **Komplett-$k$-äre Bäume** $T_{k,d}$: Graziös für kleine $k, d$ (keine allgemeine Aussage).

### 4.4 Verworfene Ansätze

#### Zufällige Labellings und probabilistische Methode

**Alon (1993)** und nachfolgende Autoren versuchten, die **probabilistische Methode**
anzuwenden: Wähle das Labelling zufällig und zeige, dass die Wahrscheinlichkeit, ein
graziöses Labelling zu erhalten, positiv ist. Das Problem: Die Bedingung, dass alle
Kantengewichte $\{|f(u)-f(v)|\}$ verschieden und gleich $\{1,\ldots,m\}$ sind, ist eine
Menge von **positiven Abhängigkeiten** zwischen den Labellings benachbarter Knotenpaare.
Das Lovász Local Lemma erfordert "lokale Unabhängigkeit", die hier strukturell verletzt wird,
da die Baum-Topologie globale Abhängigkeiten erzeugt.

#### Polynomielle Methode (Nullstellensatz-Ansatz)

**Alon's Combinatorial Nullstellensatz (1999)** ist ein mächtiges algebraisches Werkzeug:
Wenn ein Polynom $P(x_1, \ldots, x_n)$ einen Term $\prod x_i^{t_i}$ mit maximalem Grad
enthält und die Koeffizientenbedingungen erfüllt, folgt Existenz einer Nullstelle in einem
Produktraum. Für graziöse Labellings formuliert man:

$$P_G(f) = \prod_{\{u,v\} \in E} (f(u) - f(v)) \cdot \prod_{1 \leq i < j \leq m} \left((|f(u_i)-f(v_i)|) - (|f(u_j)-f(v_j)|)\right)$$

Der Nullstellensatz erfordert, dass $P_G$ einen dominanten Monomial-Term hat. Für
allgemeine Bäume ist dies schwer zu verifizieren, da die Graphstruktur des Baums direkt
in die Monomial-Koeffizienten eingeht.

#### Lineare Programmierung und Integer-LP

**Formulation als ILP**: Variablen $f(v) \in \{0,\ldots,m\}$, Nebenbedingungen für
Injektivität und Kantengewichte. Das ILP ist lösbar für $n \leq 35$, aber die LP-Relaxierung
hat keine ganzzahligen Eckpunkte für allgemeine Bäume — die LP-Heuristik bricht also
für große $n$ zusammen.

**Warum scheitern Computational Approaches bei großen Bäumen?** Der Suchraum wächst wie
$(m+1)^n = (n)^n$, und auch mit Branch-and-Bound und DPLL-ähnlichen Methoden ist die
Tiefe des Suchbaums exponentiell in $n$. Keine bekannte Pruning-Strategie reduziert
die Komplexität subexponentiell.

#### Strukturelle Induktion

Naiver Induktionsbeweis: Zeige, wenn $T$ graziös ist und man ein Blatt hinzufügt, bleibt
$T' = T + \text{Blatt}$ graziös. **Gegenbeispiel-Tendenz**: Graceful Labellings sind nicht
monoton unter Blatt-Hinzufügung — ein graziöses Labelling von $T$ kann nicht einfach
zu einem von $T'$ erweitert werden, da die Label-Menge von $\{0,\ldots,m\}$ auf
$\{0,\ldots,m+1\}$ wächst und alle bestehenden Gewichte verändert werden müssen.

### 4.5 Aktuelle Entwicklungen (2020–2024)

- **Ban und Sherratt (2020)**: Verbesserung der Computerverifikation auf $n \leq 35$.
  Nutzung von **constraint propagation** (arc consistency, GAC) und symmetry breaking
  zur Reduktion des Suchraums um Faktor $> 10^6$ gegenüber naiver Suche.
- **Verbindung zu Graph-Dekompositions**: Neue Beweise, dass bestimmte Baumfamilien
  (z. B. Bäume mit Durchmesser $\leq 5$) stets graziös sind, via **Ringel–Kotzig-Zerlegungsrahmen**.
- **Graceful vs. Harmonisch**: Ein Graph heißt **harmonisch** (harmonious), wenn
  $f: V \to \mathbb{Z}_m$ (mod $m$) mit $\{(f(u)+f(v)) \mod m\}$ bijectiv auf $\mathbb{Z}_m$
  ist. Graham–Sloane (1980) beweisen, dass alle Bäume harmonisch sind (außer $K_2$).
  Dies ist der einzige bekannte "Positiv-Beweis" für alle Bäume, aber harmonisch ist
  schwächer als graziös und liefert keine Transferprinzipien.
- **Machine Learning für Labelling-Heuristiken**: Neuronale Netze (Graph Neural Networks,
  GNNs) wurden trainiert, um graziöse Labellings für mittlere Baumgrößen ($n \leq 100$)
  zu finden. Dies hat keinen Beweiswert, aber verbesserte Heuristiken für die Grundlagenforschung.
- **Renaming Conjecture (Brinkmann–Mckay 2022)**: Eine verwandte Vermutung besagt, dass
  für **alle Bäume** ein sogenanntes *$\alpha$-Labelling* (eine Verschärfung von graceful)
  existiert — dies wurde für Raupen und Spinnen bewiesen, ist aber allgemein offen.

### 4.6 Offene Teilprobleme

- **Hauptvermutung** (Ringel 1963, Rosa 1967): Alle Bäume sind graziös.
- **$\alpha$-Labelling (bipartites graceful)**: Existiert für jeden Baum eine **bipartite
  graceful**-Variante mit zusätzlicher Bedingung $\max_{A} f(v) < \min_{B} f(v)$ für die
  Bipartition $(A,B)$ des Baums?
- **Graziöser Durchmesser-2-Baum**: Für Bäume mit Durchmesser $d = 2$ (Sterne) ist die
  Aussage trivial. Für $d = 3, 4$ bekannt. Was ist die kritische Schwelle?
- **Zufällige Bäume**: Mit welcher Wahrscheinlichkeit ist ein zufälliger gelabelter Baum
  auf $n$ Knoten (Cayley-Baum) graziös? Simulationen legen nahe, dass fast alle
  zufälligen Bäume graziös sind, aber kein Beweis existiert.
- **Komplexitätsklasse**: Ist die Entscheidungsfrage "Ist $G$ graziös?" in P, NP-vollständig,
  oder gibt es für Bäume einen Polynomzeit-Algorithmus?

---

## 5. Verbindungen zwischen den Problemen

Die vier behandelten Probleme sind nicht isoliert — zwischen ihnen bestehen tiefe
strukturelle Verbindungen:

**PFR und Lonely Runner**: Beide nutzen **Bohr-Mengen** als Hauptwerkzeug. Die Grundstruktur
eines Bohr-Sets $\text{Bohr}(\Lambda, \rho)$ taucht in der PFR-Theorie als "strukturierte
Teilmenge" auf, während bei Lonely Runner die "guten Zeiten" ein Bohr-artiges Muster haben.
Ein vollständiger Bohr-Set-Struktursatz könnte beide Probleme weitgehend lösen.

**Frankl und Graceful Trees**: Beide Probleme betreffen kombinatorische Familien und deren
Abzählungseigenschaften. Die **Frequenzfunktion** $f(x)$ bei Frankl entspricht konzeptuell
einem Gewichtslabelling auf Mengenelementen. Ob Frankl-Techniken (Entropie-Argumente)
für graziöse Bäume anwendbar sind, ist eine offene Forschungsfrage.

**PFR und Frankl**: Beide betreffen "strukturelle Regelmäßigkeit" in kombinatorischen
Objekten (Summenmengen vs. Mengenfamilien) und nutzen **Entropie-** und **Fourier-Methoden**.
Gowers' Uniformitätsnormen sind ein Fourier-Analogon zu Gilmers Shannonentropyansatz.

---

## 6. Fazit und Forschungsausblick

Die vier behandelten Probleme repräsentieren den Grenzbereich zwischen **endlicher
Kombinatorik**, **additiver Zahlentheorie**, **Fourier-Analysis** und **Graphentheorie**.
Sie zeigen ein gemeinsames Muster: Elementare Formulierungen verbergen tief-liegende
Strukturfragen, die erst durch den Einsatz analytischer, algebraischer oder probabilistischer
Methoden zugänglich werden.

Der **PFR-Beweis (2023)** von Gowers–Green–Manners–Tao ist das eindrücklichste Beispiel
dafür, wie ein "einfach formuliertes" Problem nach Jahrzehnten durch die Kombination von
Fourier-Analysis, Entropieargumenten und Inkrementmethoden gelöst werden kann. Die Frage,
welche der verbleibenden drei Probleme (Frankl, Lonely Runner, Graceful Trees) als nächstes
einem vollständigen Beweis zugänglich sein wird, ist eine der spannendsten offenen Fragen
in der aktuellen Kombinatorik.

Besonders vielversprechend erscheinen:
1. **Frankl-Vermutung**: Der informationstheoretische Ansatz (Gilmer 2022, Chase–Lovász 2023)
   hat die Schranke auf $0.38$ gebracht. Eine weitere Verfeinerung auf $0.49$ oder direkt
   $0.5$ erscheint mit einer Weiterentwicklung der Entropie-Methoden möglich.
2. **Lonely Runner für $k=8$**: Ein Computer-unterstützter Beweis mit Fallunterscheidung
   (ähnlich dem 4-Farben-Satz) könnte den nächsten Schritt bringen, auch wenn er keinen
   eleganten allgemeinen Beweis liefert.

---

*Diese Datei gehört zum Projekt `specialist-maths` (Mathematik-Spezialist Claude).*
*Autor: Michael Fuhrmann | Stand: 2026-03-13*
