# Gruppe B: Offene Probleme der analytischen Zahlentheorie
**Autor**: Michael Fuhrmann
**Erstellt**: 2026-03-12
**Letzte Änderung**: 2026-03-12

> Umfassende Forschungsübersicht zu fünf fundamentalen offenen Problemen der analytischen Zahlentheorie. Für jedes Problem werden Geschichte, gescheiterte Ansätze, aktuelle Ergebnisse und neue Forschungsrichtungen dargestellt.

---

## Problem 1: Das Lehmer-Mahler-Maß-Problem

### Geschichte und Problemstellung

Derrick H. Lehmer stellte 1933 die Frage, ob es algebraische ganze Zahlen $\alpha$ mit $1 < M(\alpha) < M(\lambda_0)$ gibt, wobei $M(\lambda_0) = 1{,}17628...$ das Mahler-Maß des sogenannten **Lehmer-Polynoms** ist:

$$\lambda_0(x) = x^{10} + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1.$$

Das Mahler-Maß eines Polynoms $f(x) = a_d \prod_{i=1}^d (x - \alpha_i)$ ist definiert als
$$M(f) = |a_d| \prod_{i=1}^d \max(1, |\alpha_i|).$$

Das Lehmer-Polynom ist das kleinste bekannte Salem-Polynom. **Salem-Zahlen** sind positive reelle algebraische ganze Zahlen $> 1$, deren sämtliche Konjugierte entweder auf dem Einheitskreis liegen oder mit ihrem reziproken Wert zusammenfallen — sie liefern die kleinsten bekannten Mahler-Maße $> 1$.

### Gescheiterte Ansätze

**Dobrowolski-Schranke (1979):** E. Dobrowolski bewies die bis heute beste asymptotische untere Schranke:
$$M(\alpha) \geq 1 + c \cdot \left(\frac{\log \log d}{\log d}\right)^3,$$
zunächst mit $c = 1/1200$, später verbessert auf $c = 9/4 - \varepsilon$. Diese Schranke wächst jedoch nur logarithmisch; sie beweist nicht die Lücke zwischen 1 und $1{,}17628$. Die Schranke ist "zu weich", um den konkreten Sprung zu 1,17628 zu erfassen — sie geht für $d \to \infty$ gegen null.

**Potts-Verfahren / numerische Suche:** Computerprogramme haben alle Polynome bis Grad 180 mit $M < 1{,}3$ klassifiziert (Mossinghoff 2001, mit Beiträgen von Boyd, Flammang, Smyth u. a.). Für Grad $\leq 40$ ist die Liste vollständig. Kein Polynom mit $M \in (1, 1{,}17628)$ wurde gefunden. Rein rechnerische Ansätze skalieren jedoch exponentiell und werden mit wachsendem Grad unbeherrschbar.

**Algebraische Potenz-Methoden:** Versuche, über den Ring $\mathbb{Z}[x]/(f)$ direkte untere Schranken zu erzwingen, scheitern an der Struktur reziprok-symmetrischer Polynome.

### Das Bogomolov-Problem und arakelov-geometrische Verbindung

Bogomolov (1981) formulierte ein verwandtes Problem in der arithmetischen Geometrie: Auf einer abelschen Varietät $A$ über einem Zahlkörper besitzt jede abgeschlossene Untervarietät, die keine Torsionsvarietät ist, eine untere Schranke für die **Weil-Höhe** $h$ ihrer Punkte. Ullmo (1998) und Zhang (1998) bewiesen die **Bogomolov-Vermutung** für abelsche Varietäten über Zahlkörpern mittels Arakelov-Geometrie und Äquidistributionssätzen. Das logarithmische Weil-Maß ist dabei eng mit dem Mahler-Maß verwandt: $h(\alpha) = \frac{1}{\deg \alpha} \log M(\alpha)$. Die Bogomolov-Eigenschaft einer Körpererweiterung — positive untere Schranke für alle nicht-Torsionselemente — ist eine direkte Verallgemeinerung von Lehmers Problem.

### Smyths Arbeit und Salem-Zahlen

Smyth klassifizierte nichtreziprokes Minimum: Für **nicht-reziproke** Polynome gilt $M \geq \theta_0 \approx 1{,}3247$ (Lehmer-Siegel-Satz für Pisot-Zahlen). Das eigentliche Problem liegt also bei den **reziproken** Polynomen, wo Salem-Zahlen auftreten. Smyth zeigte 1971, dass nicht-reziproke algebraische ganze Zahlen immer $M \geq \theta_0$ erfüllen. Die 47 kleinsten bekannten Salem-Zahlen haben alle Grad $\leq 44$ (Mossinghoff, 2001).

**Borwein-Hare-Mossinghoff (2016):** Für Polynome mit ausschließlich ungeraden Koeffizienten gilt $M \geq 5^{1/4} \approx 1{,}4953$, also deutlich über der Lehmer-Schranke. Dies zeigt, dass Koeffizientenbeschränkungen die Lücke schließen können — aber für allgemeine Polynome bleibt die Vermutung offen.

### Aktueller Stand (2022–2024)

Ein 2024 auf arXiv erschienener Artikel liefert neue Lehmер-artige untere Schranken für spezifische Polynomklassen (arXiv:2402.14771). Verger-Gaugry (2024) bewies eine nicht-triviale Minoration für die Menge der Salem-Zahlen (HAL:hal-04386425). Behauptete vollständige Beweise der Lehmer-Vermutung (Verger-Gaugry, HAL:hal-02322497, 2019) wurden von der Fachgemeinschaft bislang nicht als abschließend anerkannt. Die Vermutung gilt weiterhin als eines der tiefsten offenen Probleme der algebraischen Zahlentheorie.

---

## Problem 2: Ungerade vollkommene Zahlen (OPN)

### Geschichte und Problemstellung

Eine **vollkommene Zahl** $n$ erfüllt $\sigma(n) = 2n$, wobei $\sigma$ die Teilersummenfunktion ist. Alle bekannten vollkommenen Zahlen sind gerade (Euklid-Euler-Form: $2^{p-1}(2^p - 1)$). Die Frage nach ungeraden vollkommenen Zahlen ist eines der ältesten ungelösten Probleme der Mathematik — mindestens 2000 Jahre alt.

**Eulers Struktursatz:** Euler bewies, dass jede ungerade vollkommene Zahl die Form hat:
$$n = p^a \cdot m^2, \quad p \equiv a \equiv 1 \pmod{4},$$
wobei $p$ der sogenannte **Euler-Primteiler** ist ($\gcd(p, m) = 1$, $p$ prim).

### Stärkste bekannte Schranken

**Ochem-Rao (2012):** Jede ungerade vollkommene Zahl erfüllt $n > 10^{1500}$ und hat mindestens 101 (nicht notwendig verschiedene) Primteiler sowie mindestens 9 verschiedene Primteiler. Die größte Komponente (d. h. $p^k$ mit $p$ prim) ist $> 10^{62}$.

**Nielsen (2015):** Jede OPN hat mindestens 10 verschiedene Primteiler; falls sie nicht durch 3 teilbar ist, mindestens 12.

**Luca-Pomerance (2010):** Der zweitgrößte Primteiler einer OPN überschreitet $10^4$.

**Acquaah-Konyagin (2012):** Wenn $N = p_1^{a_1} \cdots p_k^{a_k}$ eine OPN mit $p_1 < \ldots < p_k$ ist, gilt $p_k < (3N)^{1/3}$. Dies liefert eine obere Schranke für den größten Primteiler in Abhängigkeit von $N$.

### Warum scheitern L-Funktionen-Ansätze?

L-Funktionen-Methoden würden eine analytische Kontrolle über $\sigma(n)/n$ erfordern. Das Problem ist multiplikativ, nicht additiv: $\sigma$ ist vollständig multiplikativ auf Primzahlpotenzen, aber die Gleichung $\sigma(n) = 2n$ verbindet viele Primfaktoren simultan. Es gibt keine sinnvolle L-Funktion, die die globale Gleichung $\prod_{p^a \| n} \frac{\sigma(p^a)}{p^a} = 2$ in eine analytische Aussage überführt. Zudem ist der Lösungsraum diskret, nicht kontinuierlich.

### Spiro-Korollar und verwandte Ergebnisse

Spiro zeigte, dass eine OPN mindestens $\Omega(n) \geq 75$ Primteiler (mit Vielfachheit) besitzen muss. Nielsen (2015) stärkte dies auf $\Omega(n) \geq 91$. Die Methode basiert auf "Effizienzargumenten": Für jeden Primteiler $p^a \| n$ misst man $\mathrm{eff}(p,a) = \frac{\ln(\sigma(p^a)/p^a)}{\ln(p^a)}$. Da $\sum \mathrm{eff}(p,a) = \ln 2$, können kleine Primteiler effizienter beitragen, was untere Schranken an die Anzahl der Faktoren erzwingt.

### Aktuelle Entwicklungen (2022–2024)

Eine 2024-Arbeit behauptet einen Widerspruchsbeweis, dass OPN nicht durch 3 teilbar sein können (Preprints.org, Oktober 2024). Ochers Website (LIRMM, Datenbank 12/2022) dokumentiert fortlaufende Faktorisierungsarbeiten mit 660 MB komprimierten Primteiler-Daten. Eine 2024-Arbeit in *INTEGERS* präsentiert verbesserte obere Schranken. Der Nachweis einer OPN oder ein vollständiger Unmöglichkeitsbeweis scheinen mit heutigen Methoden unerreichbar. Das Problem gilt als möglicherweise "unbeweisbar" mit klassischer algebraischer Zahlentheorie allein.

---

## Problem 3: Mersenne-Primzahlen — gibt es unendlich viele?

### Geschichte und Problemstellung

**Mersenne-Primzahlen** haben die Form $M_p = 2^p - 1$ für prim $p$. Die Geschichte reicht bis in die Antike; Marin Mersenne untersuchte sie systematisch im 17. Jahrhundert. Stand März 2026: Es gibt **52 bekannte** Mersenne-Primzahlen. Die bisher größte, $2^{136{.}279{.}841} - 1$ (mit 41.024.320 Stellen), wurde am 12. Oktober 2024 von Luke Durant via **GIMPS** (Great Internet Mersenne Prime Search) entdeckt — die erste mit einem Exponenten $> 10^8$.

### Die Lenstra-Pomerance-Wagstaff-Vermutung

Pomerance und Lenstra zeigten 1980 unabhängig voneinander, dass die erwartete Anzahl der Mersenne-Primzahlen $\leq x$ asymptotisch
$$\frac{e^\gamma}{\log 2} \cdot \log \log x$$
beträgt ($\gamma \approx 0{,}5772$ ist die Euler-Mascheroni-Konstante). Wagstaff formulierte dies 1983 präzise: Der Erwartungswert für den Quotienten aufeinanderfolgender Mersenne-Exponenten ist $e^{1/e^\gamma} \approx 1{,}4758$ (nicht $3/2$ wie Erhardt vermutete). Die Heuristik entsteht durch Kombination des **Primzahlsatzes** für $p$ und Mertens' Theorem für die Primteilerprodukte von $2^p - 1$.

### Warum scheitern Siebmethoden?

Das fundamentale Hindernis liegt in der **exponentiellen Wachstumsrate**: $M_p = 2^p - 1$ wächst exponentiell in $p$. Klassische Siebmethoden (Eratosthenes, Brun, Selberg) benötigen eine **multiplikative Dichtefunktion** für die untersuchte Folge. Bei $M_p$ ist diese Dichte nicht multiplikativ — das Prüfen auf $p$-Teilbarkeit von $2^p - 1$ involviert diskrete Logarithmen modulo variierender Primzahlen. Granville und andere zeigten, dass **stark divisible Folgen** wie Lucas-Folgen exponentiell wachsen und die Sieb-Axiome (nicht-negative multiplikative Gewichtsfunktion) systematisch verletzen. Die große Siebungleichung mit Exponentialfunktionen liefert zwar Verteilungsaussagen für Mersenne-Zahlen modulo Primzahlen (IEEE 2017), aber keinen Primzahlnachweis.

### Verbindung zum Zsygmondy-Satz

Der **Zsygmondy-Satz** (1892) besagt, dass für $a^n - b^n$ (mit $\gcd(a,b)=1$, $n \geq 3$) fast immer ein **primitiver Primteiler** existiert, d. h. einer, der $a^k - b^k$ für kein $k < n$ teilt. Für $2^p - 1$ liefert dies, dass der kleinste Primteiler $q$ von $M_p$ die Ordnung $p$ modulo $q$ hat, also $q \equiv 1 \pmod{p}$. Dies erklärt, warum primitive Primteiler von $M_p$ groß sein müssen — liefert aber keine Aussage über die Primalität von $M_p$ selbst.

### Ist das Problem "leichter" als Goldbach?

Goldbach betrifft additive Zerlegungen — ein additives Problem in einem multiplikativ strukturierten Raum. Mersenne-Primzahlen sind ein **multiplikatives Problem** in einem exponentiell strukturierten Raum. Beide gelten als vergleichbar schwer. Jedoch erlaubt die besondere Struktur von $2^p - 1$ (Lucas-Lehmer-Test!) eine deterministisch-effiziente Primalitätsprüfung, was GIMPS erst ermöglicht. Für Goldbach existiert kein analoges Testwerkzeug. In diesem Sinne ist das Mersenne-Problem experimentell zugänglicher, theoretisch aber nicht leichter.

### Aktueller Stand (2023–2026)

GIMPS hat seit 1996 alle 17 größten bekannten Mersenne-Primzahlen entdeckt. Der Sprung von $M_{82.589.933}$ (2018) zu $M_{136.279.841}$ (2024) war der erste via GPU-Beschleunigung. Lücken in der GIMPS-Suche (nicht alle Exponenten unter $10^8$ wurden geprüft) bedeuten, dass möglicherweise weitere Mersenne-Primzahlen mit Exponenten zwischen $82.589.933$ und $136.279.841$ noch unentdeckt sind.

---

## Problem 4: Bunyakovsky-Vermutung und das Selberg-Paritätsproblem

### Geschichte und Problemstellung

Viktor Bunyakovsky vermutete 1857: Sei $f \in \mathbb{Z}[x]$ irreduzibel mit positivem führenden Koeffizient und $\gcd(\{f(n) : n \in \mathbb{Z}\}) = 1$. Dann nimmt $f(n)$ unendlich viele Primwerte an. Der wichtigste Spezialfall ist $f(n) = n^2 + 1$: Diese Frage wurde von **Euler** gestellt und ist das **4. Landau-Problem** (1912).

### Iwaniec 1978: Fast-Primzahlen (P₂)

Iwaniec bewies 1978, dass $n^2 + 1$ unendlich oft Werte mit **höchstens zwei Primfaktoren** annimmt (d. h. $P_2$-Zahlen). Der Beweis verwendet den **Rosser-Sieb** mit einer Zwillingsprimzahl-artigen Struktur und exploitiert, dass $n^2 + 1$ im Körper $\mathbb{Q}(i)$ als Norm $N(n + i) = n^2 + 1$ geschrieben werden kann. Damit lassen sich Kreiszahlkörper-Methoden einsetzen.

### Warum reicht Iwaniec nicht für echte Primzahlen?

Das **Selberg-Paritätsproblem** (1949) ist das fundamentale Hindernis: Siebtechniken können die Parität der Anzahl der Primfaktoren eines Elements einer Folge nicht auflösen. Formal: Wenn alle Elemente einer Folge $A$ entweder ausschließlich Produkte einer **geraden** oder ausschließlich einer **ungeraden** Anzahl von Primfaktoren sind, kann ein Sieb keine nicht-triviale untere Schranke für $|A|$ liefern; obere Schranken sind um Faktor $\geq 2$ von der Wahrheit entfernt.

Für $n^2 + 1$: Der Übergang von $P_2$ (zwei Primfaktoren) zu $P_1$ (ein Primfaktor = echte Primzahl) erfordert genau die Auflösung der Parität. Das Sieb kann nicht unterscheiden, ob $n^2+1$ eher $p \cdot q$ oder $p$ ist — beide Typen erscheinen mit ähnlicher Häufigkeit, und das Sieb kann sie nicht trennen.

### Friedlander-Iwaniec (1997) und Heath-Brown (2001): Warum funktionieren diese?

**Friedlander-Iwaniec** bewiesen 1998, dass $x^2 + y^4$ unendlich viele Primwerte liefert. **Heath-Brown** bewies 2001 dasselbe für $x^3 + 2y^3$.

Der entscheidende Unterschied zu $n^2 + 1$: Beide Ausdrücke sind Normen in Zahlkörpern:
- $x^2 + y^4 = x^2 + (y^2)^2 = N_{\mathbb{Q}(i)/\mathbb{Q}}(x + iy^2)$: Die Variable $y^2$ ist durch eine vollständige Quadratvariable ersetzt — das Problem hat eine **bilineare Struktur**, die eine asymptotische Siebanalyse mit $x \asymp N^{1/2}$, $y \asymp N^{1/4}$ ermöglicht.
- $x^3 + 2y^3 = N_{\mathbb{Q}(\sqrt[3]{2})/\mathbb{Q}}(x + y\sqrt[3]{2})$: Hier wird das **Kubik-Zahlkörper-Sieb** eingesetzt.

Für $n^2 + 1$: Es gibt nur **eine** Variable. Es existiert keine bilineare oder multilineare Struktur, die das Paritätsproblem umgehen könnte. Der Raum $\{(x, y) : x^2 + y^4 \leq N\}$ hat Maß $\sim N^{3/4}$, deutlich dicker als die Menge $\{n : n^2 + 1 \leq N\} \sim N^{1/2}$. Diese "Dicke" erlaubt eine Kontrolle der Fehlterme, die im univariaten Fall fehlt.

### Granvilles Reformulierung

Tao (2007) und Granville haben präzisiert, dass das Paritätsproblem äquivalent zu einer **Selberg-Bedingung** ist: Siebe können keine Aussagen über Summen $\sum_{n \leq x} \lambda(f(n))$ machen, wobei $\lambda$ die Liouville-Funktion ist. Für $f(n) = n^2 + 1$ ist diese Summe nicht kontrollierbar.

### Aktueller Stand

Alle vier Landau-Probleme (Goldbach, Zwillingsprimzahlen, Legendres Vermutung, $n^2+1$) sind ungeklärt. Neueste Arbeiten (2022–2024) im Bereich des "parity-sensitive sieve" von Friedlander-Iwaniec erweitern die Liste der Formen, für die Primwerte bewiesen werden können, aber $n^2 + 1$ als univariates Polynom bleibt außerhalb der Reichweite.

---

## Problem 5: Bateman-Horn-Vermutung

### Geschichte und Problemstellung

Paul T. Bateman und Roger A. Horn postulierten 1962 die quantitative Version der Bunyakovsky-Vermutung für Systeme von Polynomen. Seien $f_1, \ldots, f_k \in \mathbb{Z}[x]$ paarweise irreduzibel, nicht-konstant, mit positiven Führungskoeffizienten. Definiere
$$\omega(p) = \#\{n \bmod p : p \mid f_1(n) \cdots f_k(n)\}.$$
Die **Bateman-Horn-Vermutung** besagt:
$$\#\{n \leq N : f_i(n) \text{ prim } \forall i\} \sim \frac{N}{(\log N)^k} \cdot \prod_p \frac{(1-\omega(p)/p)}{(1-1/p)^k} = \frac{N \cdot C(f_1,\ldots,f_k)}{(\log N)^k}.$$
Das Produkt $C(f_1, \ldots, f_k)$ konvergiert für zulässige Tupel.

### Spezialfälle und was bewiesen ist

| Spezialfall | Status |
|---|---|
| $k=1$, $f_1(n) = n$ (PNT) | Bewiesen (Hadamard/de la Vallée-Poussin, 1896) |
| $k=1$, $f_1(n) = an+b$, $\gcd(a,b)=1$ | Bewiesen (Dirichlet, 1837) |
| Funktion-Feld-Analogon (endliche Körper) | Bewiesen |
| $k \geq 2$ über $\mathbb{Z}$ | **Offen** |

### Warum ist Green-Tao kein Beweis für Bateman-Horn?

Der **Green-Tao-Satz** (2004): Die Primzahlen enthalten arithmetische Progressionen beliebiger Länge. Dies ist eine Aussage über **Existenz** (qualitativ), nicht über **Dichte** (quantitativ). Bateman-Horn macht präzise quantitative Vorhersagen über die Anzahl der $n \leq N$ mit allen $f_i(n)$ prim.

Green-Tao verwendet **Fourier-Analysis auf $\mathbb{Z}/N\mathbb{Z}$** und Szemerédi-Regularitäts-Lemma-Methoden. Diese erfassen globale additive Strukturen der Primzahlen, sagen aber nichts über die lokale Verteilung für spezifische Polynomfamilien.

Tao-Ziegler (2008) verallgemeinerten Green-Tao auf Polynomprogressionen (z.B. $n, n+d, n+d^2, \ldots$) — aber dies beweist Existenz, keine Dichte.

### Maynard-Tao bounded gaps (2013) und mehrdimensionale Bateman-Horn

Maynard und Tao bewiesen 2013 unabhängig, dass für **jedes $m$** unendlich viele $n$ existieren, sodass unter $n+h_1, \ldots, n+h_k$ (für geeignete Tupel) mindestens $m$ Primzahlen sind. Das Polymath-Projekt reduzierte den Mindestabstand von Primzahlen auf 246.

Für Bateman-Horn mit $k$ linearen Polynomen $f_i(n) = n + h_i$: Der Maynard-Tao-Ansatz zeigt, dass **viele** der $f_i(n)$ simultan prim sein können, liefert aber keine quantitative Dichte für **alle** $f_i(n)$ gleichzeitig. Die Methode ist kombinatorisch (weighted sieve), während Bateman-Horn eine multiplikative Produktformel benötigt. Der Übergang von "viele simultan prim" zu "alle simultan prim" mit exakten Asymptotiken ist das entscheidende Hindernis.

### Hardy-Littlewood Conjecture B und Granvilles Reformulierung

Die **Hardy-Littlewood Vermutung B** (1923) ist der Spezialfall $k=2$, $f_1(n) = n$, $f_2(n) = n+2$ (Zwillingsprimzahlen): $\#\{p \leq N : p+2 \text{ prim}\} \sim 2C_2 \cdot N/(\log N)^2$, wobei $C_2 \approx 0{,}6601618...$ die **Zwillingsprimzahlkonstante** ist.

Granville reformulierte die Bateman-Horn-Vermutung im Rahmen von **Galois-Kategorien**: Die lokalen Faktoren $\omega(p)/p$ können als Dichten von Frobenius-Konjugationsklassen in den Zerlegungsgruppen der $f_i$ interpretiert werden. Dies verbindet Bateman-Horn mit dem **Tschebotareff-Dichtesatz** für nicht-abelsche Erweiterungen.

### Soundararajan zu oberen Schranken

Montgomery und Soundararajan zeigten, dass für singuläre Serien $R_k(h) = \sum_{\substack{h_1+\ldots+h_k = h \\ \text{zulässig}}} C(n, n+h_1, \ldots)$ gilt:
$$R_k(h) = \mu_k(-h \log h + Ah)^{k/2} + O_k(h^{k/2-1/(7k)+\varepsilon}),$$
wobei $\mu_k$ Gauß'sche Momente sind. Für ungerades $k$ vermuten Soundararajan et al., dass $R_k(h) \asymp h^{(k-1)/2}(\log h)^{(k+1)/2}$; obere Schranken mit dem richtigen $h$-Exponenten wurden für $k=3$ bewiesen (2022, Forum of Mathematics Sigma). Dies gibt Evidenz für die statistische Verteilung, ersetzt aber keinen direkten Beweis.

### Funktion-Feld-Analogon: Was dort funktioniert

Über Funktionenkörpern $\mathbb{F}_q(t)$ ist Bateman-Horn ein **Theorem** (folgt aus Weil-Deligne, Etale-Kohomologie). Die Beweise nutzen algebraische Geometrie über endlichen Körpern. Der Wechsel zu $\mathbb{Z}$ scheitert daran, dass $\mathbb{Z}$ keinen absoluten Grundkörper $\mathbb{F}_q$ besitzt und die entsprechende motivische Kohomologie-Theorie fehlt.

### Aktueller Stand (2022–2026)

Eine 2024-Arbeit in den *AMS Notices* (Garcia) gab einen zugänglichen Überblick. Aktuelle Forschung konzentriert sich auf partielle Bateman-Horn-Resultate für spezielle Polynomfamilien, Funktion-Feld-Analogien als Testfeld, und die Verbindung zu Rankin-Selberg-L-Funktionen. Ein vollständiger Beweis erscheint in absehbarer Zeit nicht erreichbar.

---

## Querverbindungen zwischen den Problemen

| Problem | Verbindung |
|---|---|
| Lehmer ↔ Bogomolov | Mahler-Maß ist logarithmische Weil-Höhe; Bogomolov-Eigenschaft |
| OPN ↔ Mersenne | Alle geraden vollkommenen Zahlen sind Mersenne-basiert; OPN-Beweis würde Charakterisierung perfekter Zahlen vollenden |
| Bunyakovsky ↔ Bateman-Horn | Bateman-Horn ist die quantitative Version; Bunyakovsky der qualitative Spezialfall $k=1$ |
| Selberg-Parität ↔ alle | Das Paritätsproblem blockiert Siebansätze bei Mersenne ($2^p-1$), OPN (Primfaktoranzahl) und $n^2+1$ |
| Green-Tao ↔ Bateman-Horn | Green-Tao beweist Existenz in linearen Formen (AP), Bateman-Horn verlangt Dichten in polynomiellen Formen |

---

## Literaturverweise

- Lehmer, D.H. (1933): Factorization of certain cyclotomic functions. *Ann. Math.*
- Dobrowolski, E. (1979): On a question of Lehmer. *Acta Arithmetica*, 34.
- Borwein, P., Dobrowolski, E., Mossinghoff, M. (2007): Lehmer's problem for polynomials. *Ann. Math.*, 166.
- Verger-Gaugry, J.-L. (2024): Non-trivial minoration for Salem numbers. HAL:04386425.
- Euler, L. (1747): Variae observationes circa series infinitas. (Vollkommene Zahlen)
- Ochem, P., Rao, M. (2012): Odd perfect numbers are greater than $10^{1500}$. *Math. Comp.*
- Nielsen, P.P. (2015): Odd perfect numbers, Diophantine equations. *Math. Comp.*
- Acquaah, P., Konyagin, S. (2012): On prime factors of odd perfect numbers. *Int. J. Number Theory*, 8(6).
- Wagstaff, S.S. (1983): Divisors of Mersenne numbers. *Math. Comp.*
- Browning, T.D., Granville, A. (2024): Strong divisibility sequences and sieve methods. *Mathematika*.
- Iwaniec, H. (1978): Almost-primes represented by quadratic polynomials. *Acta Arithmetica*, 24.
- Friedlander, J., Iwaniec, H. (1997/1998): Using a parity-sensitive sieve. *PNAS* / *Ann. Math.*
- Heath-Brown, D.R. (2001): Primes represented by $x^3 + 2y^3$. *Acta Math.*, 186.
- Bateman, P.T., Horn, R.A. (1962): A heuristic asymptotic formula. *Math. Comp.*, 16.
- Green, B., Tao, T. (2008): The primes contain arbitrarily long arithmetic progressions. *Ann. Math.*
- Maynard, J. (2015): Small gaps between primes. *Ann. Math.*, 181.
- Montgomery, H., Soundararajan, K. (2004): Primes in short intervals. *Commun. Math. Phys.*
- Garcia, S.R. (2024): What is... the Bateman-Horn Conjecture? *AMS Notices*, Oktober 2024.
- Smyth, C.J.: The Mahler measure of algebraic numbers: a survey. Edinburgh, 2007.
- Mossinghoff, M.J. et al. (2001): Minimal Mahler measures. *Experimental Mathematics*, 17(4).
- Bogomolov, F.A. (1981): Points of finite order on abelian varieties.
- Ullmo, E. (1998): Positivité et discrétion des points algébriques des courbes. *Ann. Math.*
- Zhang, S. (1998): Equidistribution of small points on abelian varieties. *Ann. Math.*
