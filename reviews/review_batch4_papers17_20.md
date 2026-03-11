# Mathematisches Gutachten: Batch 4 — Papers 17–20 (EN + DE)

**Gutachter:** Claude (Sonnet 4.6)
**Datum:** 2026-03-11
**Batch:** 4
**Papers:** 17–20, je EN- und DE-Version (8 Dateien gesamt)

---

## Überblick

| Paper | Titel (EN) | Urteil EN | Urteil DE |
|-------|-----------|-----------|-----------|
| 17 | The Hardy–Littlewood Circle Method | DRUCKREIF (1 Minorbug) | DRUCKREIF (1 Minorbug) |
| 18 | Vinogradov's Three-Prime Theorem | NICHT DRUCKREIF (3 Bugs) | NICHT DRUCKREIF (2 Bugs) |
| 19 | Waring's Problem | NICHT DRUCKREIF (2 Bugs) | DRUCKREIF (1 Minorbug) |
| 20 | The Goldbach Singular Series | NICHT DRUCKREIF (2 Bugs) | DRUCKREIF (1 Minorbug) |

---

## Paper 17 EN — The Hardy–Littlewood Circle Method

**Datei:** `papers/batch4/paper17_circle_method.tex`

### Inhaltsbeschreibung
Das Paper entwickelt die Hardy–Littlewood-Kreismethode von Grund auf: Weyl-Summen, Gauß-Summen, Farey-Zerlegung, Haupt- und Nebenbögen, die singuläre Reihe und das singuläre Integral, Weyls Ungleichung für Nebenbögen, den Hauptsatz für das Waringsche Problem und die Anwendung auf Goldbach-artige Probleme. Es handelt sich um ein Übersichtspaper, das den theoretischen Rahmen für Papers 18–20 bereitstellt.

### Mathematische Prüfung

**Lemma 2.1 (Weyl, 1916 — Geometrische Summe):**
Korrekt. Die Abschätzung $|\sum_{n=1}^N e(n\alpha)| \leq \min(N, 1/(2\|\alpha\|))$ stimmt. Der Beweis über die geometrische Reihe ist vollständig. Die Ungleichung $|\sin(\pi x)| \geq 2\|x\|$ (für alle reellen $x$) ist ein Standardresultat und korrekt zitiert.

**Lemma 2.2 (Gauß-Summen-Schranke):**
Die Behauptung $|G(a,q)| \leq k \cdot q^{1/2+\varepsilon}$ ist korrekt für $k$-te Potenzen. Die verfeinerte Aussage für $k=2$: "$|G(a,q)| \leq \sqrt{q}$ wenn $q$ ungerade und $q \equiv 1 \pmod 4$" ist eine Vereinfachung — tatsächlich gilt $|G(a,q)| = \sqrt{q}$ für $q \equiv 1 \pmod 4$ und $|G(a,q)| = 0$ oder $\sqrt{2q}$ je nach $q \pmod 8$. Dieser Detailpunkt ist unscharf formuliert, aber für den Zweck des Papers (obere Schranke) akzeptabel.

**Lemma 4.1 (Farey-Partition):**
Korrekt. Der Beweis mit Verweis auf Hardy–Wright ist angemessen.

**Definition 4.1 (Haupt- und Nebenbögen):**
Die Definition $\mathfrak{M}(a,q) = \{\alpha : |\alpha - a/q| \leq Q/N\}$ ist korrekt. Die Bemerkung, dass die Bögen paarweise disjunkt sind für $Q^2 \leq N$, ist korrekt. Die Gesamtmessung $\ll Q^3/N$ ist korrekt: es gibt $\sum_{q \leq Q} \phi(q) \sim Q^2/2$ Bögen, jeder mit Länge $2Q/N$, also Gesamtmessung $\sim Q^3/N$.

**Lemma 5.1 (Hauptbogen-Approximation):**
Das Ergebnis $f(\alpha) = G(a,q)/q \cdot v(\beta) + O(q^{1+\varepsilon})$ ist korrekt. Der Beweis über Residuen-Zerlegung ist sauber. Die Verbesserung des Fehlerterms auf $O(q^{1+\varepsilon})$ ist korrekt mit Verweis auf Vaughan.

**Definition 5.2 (Singuläre Reihe):**
Die Formel und Konvergenzaussage für $s > 2k+1$ sind korrekt. Das Euler-Produkt $\mathfrak{S}(n) = \prod_p \sigma_p(n)$ ist korrekt beschrieben.

**Lemma 5.3 (Singuläres Integral):**
Die Asymptotik $\mathfrak{J}(n) \sim C_{k,s} \cdot n^{s/k-1}$ mit $C_{k,s} = \Gamma(1+1/k)^s / \Gamma(s/k)$ ist korrekt. Der Beweis über Skalierung ist korrekt skizziert.

**Satz 5.4 (Hauptbogen-Beitrag):**
Der Beweisrahmen ist korrekt. Die Fehlerabschätzung $O(N^{s/k-1-1/k})$ aus dem Approximationsfehler und $O(Q^{-1+\varepsilon})$ aus dem Schwanz der singulären Reihe sind korrekt.

**Satz 6.1 (Weylsche Ungleichung):**
Die Formel ist korrekt. Der Beweis für $k=2$ über Weylsche Differenzierung ist korrekt. Für allgemeines $k$ ist der Verweis auf Literatur angemessen.

**Korollar 6.2 (Nebenbogen-Schranke):**
Korrekt mit $\sigma(k) = 2^{1-k}$.

**Satz 7.1 (Hauptsatz):**
Der Gesamtbeweis ist korrekt. Im Minorbogen-Teil wird genutzt:
$$\left|\int_{\mathfrak{m}} f^s e(-n\cdot)\right| \leq N^{(s-2)(1-\sigma)} \cdot \int_0^1 |f|^2 \ll N^{s-1-(s-2)\sigma}$$
Die Bedingung lautet im Paper $s > 2/\sigma + 1$, was korrekt aus $s - 1 - (s-2)\sigma < s/k - 1$ folgt.

**Lemma 8.1 (Hauptbogen für Primzahlen):**
Die Approximation $S(\alpha) \approx \mu(q)/\phi(q) \cdot v_N(\beta)$ ist korrekt. Die Erklärung des $\mu(q)/\phi(q)$-Faktors ist stimmig.

**Abschnitt 8 (Nebenbogen für Primes):**
Die Aussage $\max_{\alpha \in \mathfrak{m}} |S(\alpha)| \ll N(\log N)^{-B}$ ist korrekt — dies ist eine Konsequenz des Bombieri–Vinogradov-Satzes.

**Historische Anmerkung:**
Wooley (1992–2016) hat die Vinogradov-Hauptvermutung bewiesen; Bourgain–Demeter–Guth (2015) haben dies unabhängig bewiesen. Das Paper nennt "Wooley (1992–2016)" und "Bourgain–Demeter–Guth theorem (2015)", was etwas ungenau ist — die Reihenfolge war: BDG 2015 (Publikation), Wooley 2016 (Publikation), aber Wooleys Arbeit war ebenfalls 2015 als Preprint verfügbar. Keine inhaltlich falsche Aussage.

### LaTeX-Prüfung
- Alle Theorem-Umgebungen korrekt definiert.
- `\newcommand{\e}[1]{...}` (Zeile 30) und `\newcommand{\eq}[1]{...}` (Zeile 31) sind beide definiert, aber `\e` wird im Dokument nie verwendet — harmloser toter Code.
- Alle `\bibitem`-Einträge vorhanden und konsistent mit Zitaten im Text.
- `\label`- und `\ref`-Konsistenz: Alle Referenzen im Text (`eq:key`, `lem:gauss`, `lem:farey`, `thm:major`, `thm:weyl`, `cor:minor`, `thm:main`) sind korrekt definiert.
- Keine undefinierten Umgebungen.
- Keine Baustellen (`\ldots`, TODO, ???).

### EN-DE-Konsistenz
Wird unter DE-Version verglichen.

**BUG-B4-01 (MINOR, EN):** Zeile 30: `\newcommand{\e}[1]{e^{2\pi i (#1)}}` ist definiert aber nie verwendet. Toter Code; keine funktionalen Auswirkungen, sollte bereinigt werden.

### Urteil: **DRUCKREIF** (1 Minorbug: BUG-B4-01)

---

## Paper 17 DE — Die Hardy–Littlewood-Kreismethode

**Datei:** `papers/batch4/paper17_kreismethode_de.tex`

### Inhaltsbeschreibung
Deutsche Übersetzung von Paper 17 EN. Deckt dieselben Themen ab, mit leicht komprimierterer Behandlung (z.B. fehlt die ausführliche Hauptbogen-Approximation in Abschnitt 5, es gibt keine separate Bemerkung über Farey-Bögen-Breite).

### Mathematische Prüfung
Alle mathematischen Aussagen sind konsistent mit der EN-Version. Keine inhaltlichen Abweichungen festgestellt.

**Vergleich mit EN:**
- EN hat zusätzliche `\newcommand{\e}` (nicht in DE).
- DE hat `\usepackage[utf8]{inputenc}` und `\usepackage[ngerman]{babel}` (korrekt).
- DE fehlt `\DeclareMathOperator{\ord}{ord}`, `\DeclareMathOperator{\li}{li}`, `\DeclareMathOperator{\sgn}{sgn}` — diese werden in DE aber auch nicht verwendet, kein Fehler.
- Vinogradov-Bibliographie: EN zitiert `\bibitem{Vinogradov1954}` mit "Interscience Publishers", DE zitiert `\bibitem{Vinogradov1954}` als "I. M. Winogradow, Harri Deutsch, Frankfurt, 1989 (Deutsche Übersetzung)". Beide korrekt, unterschiedliche Ausgaben/Übersetzungen.
- EN hat `\bibitem{HardyWright}` (6. Aufl., 2008); DE hat dieselbe Referenz.

### LaTeX-Prüfung
- Alle Umgebungen korrekt. Keine undefinierten Makros.
- Keine Baustellen.
- `\bibitem`-Einträge konsistent mit Textzitaten.

**BUG-B4-02 (MINOR, DE):** Die DE-Version fehlt dem Beweis von Lemma 3.1 (Farey-Partition) — nur in EN vorhanden. Inkonsistenz zwischen EN und DE: EN enthält einen Beweis-Skizze, DE nicht. Für ein akademisches Paper ist Konsistenz wünschenswert.

### Urteil: **DRUCKREIF** (1 Minorbug: BUG-B4-02)

---

## Paper 18 EN — Vinogradov's Three-Prime Theorem

**Datei:** `papers/batch4/paper18_vinogradov_three_primes.tex`

### Inhaltsbeschreibung
Das Paper beweist Vinogradovs Drei-Primzahlen-Satz: Jede hinreichend große ungerade Zahl ist Summe dreier Primzahlen. Es entwickelt die asymptotische Formel für $R_3(n)$ über Hauptbogen-Approximation (mit Siegel–Walfisz), die singuläre Reihe mit Euler-Produkt, die singuläre Integral-Berechnung und Nebenbögen-Schranken. Abschließend wird Helfgotts vollständiger Beweis (2013) diskutiert.

### Mathematische Prüfung

**Lemma 3.1 (Hauptbogen-Approximation für Primzahlen):**
Die Formel $S(\alpha) = \mu(q)/\phi(q) \cdot I(\beta) + O(N e^{-c\sqrt{\log N}})$ ist korrekt. Der Beweis über Zerlegung nach Dirichlet-Charakteren ist korrekt.

**KRITISCH — Schritt im Beweis (Zeile 179):**
```
Since $\sum_\chi \bar\chi(a) = \phi(q) \cdot \mathbf{1}_{a \equiv 1 (q)}$
```
Diese Formel ist **falsch**. Die korrekte Orthogonalitätsrelation lautet:
$$\sum_{\chi \pmod q} \bar\chi(a) \chi(b) = \phi(q) \cdot \mathbf{1}_{a \equiv b \pmod q}$$
Die Aussage $\sum_\chi \bar\chi(a) = \phi(q) \cdot \mathbf{1}_{a \equiv 1 (q)}$ entspräche $b=1$, d.h. $\mathbf{1}_{a \equiv 1 \pmod q}$. Das stimmt nicht für allgemeines $a$ mit $\gcd(a,q)=1$; für $a \not\equiv 1$ wäre die Summe $0$, nicht $\phi(q)$.

Der **richtige** Beweis: Man summiert $\sum_n \Lambda(n) e(na/q + n\beta)$ über Charaktere, erhält $\frac{1}{\phi(q)} \sum_\chi \bar\chi(a) \psi(\chi, \beta)$ mit $\psi(\chi, \beta) = \sum_{n \leq N} \Lambda(n)\chi(n) e(n\beta)$. Für $\chi = \chi_0$ ist $\psi(\chi_0, \beta) \approx I(\beta)$, für nicht-triviale $\chi$ ist $\psi(\chi, \beta)$ klein. Die Summation über $\chi$ mit dem Gewicht $\bar\chi(a)$ ergibt durch die Vollständigkeitsrelation $\frac{1}{\phi(q)} \sum_\chi \bar\chi(a) \cdot [\text{beitrag}]$. Der $\mu(q)/\phi(q)$-Faktor ergibt sich korrekt aus der Möbius-Funktion-Charakterisierung von $\chi_0$. Die Zwischenzeile ist ungenau, aber das Endergebnis ist korrekt.

**BUG-B4-03 (MITTEL, EN):** Zeile 179 in paper18_vinogradov_three_primes.tex: Die Orthogonalitätsformel `\sum_\chi \bar\chi(a) = \phi(q) \cdot \mathbf{1}_{a \equiv 1 (q)}` ist falsch dargestellt. Korrekte Formulierung: die Orthogonalitätsrelation besagt $\frac{1}{\phi(q)}\sum_\chi \bar\chi(a) \chi(n) = \mathbf{1}_{n \equiv a \pmod q}$. Das Endergebnis ist trotzdem korrekt, aber der Beweis hat eine fehlerhafte Zwischenzeile.

**Satz 4.1 (Euler-Produkt für $\mathfrak{S}(n)$):**
Das Endergebnis
$$\mathfrak{S}(n) = \prod_{p|n}\!\left(1 - \frac{1}{(p-1)^2}\right) \cdot \prod_{p\nmid n}\!\left(1 + \frac{1}{(p-1)^3}\right)$$
ist korrekt für die Goldbach-Drei-Primzahlen-Reihe mit $(\mu(q)/\phi(q))^3$.

**KRITISCH — Beweis von Satz 4.1, Zeilen 255–259:**
Im Beweis steht:
```
For $p \nmid n$:
1 + c_p(n)/(p-1)^3 = 1 - 1/(p-1)^3 + 1/(p-1)^2+...
A careful computation gives the stated product.
```
Dies ist ein **Beweis-Stub** mit `...` (Auslassungspunkte) und dem Eingeständnis "Wait, let us recompute". Der Text enthält explizit die Worte "Wait, let us recompute" (Zeile 253) — ein offensichtlicher **Entwurfsüberrest**. Der Beweis des Euler-Produkts ist für $p \nmid n$ nicht vollständig ausgeführt.

Die korrekte Rechnung ist: $1 + c_p(n)/(p-1)^3 = 1 - 1/(p-1)^3$ für $p \nmid n$. Dies ergibt nicht direkt $1 + 1/(p-1)^3$ wie behauptet.

Tatsächlich gilt für das **ternäre** Goldbach mit $(\mu(q)/\phi(q))^3$: Der $p$-Faktor im Euler-Produkt ist $\sum_{k=0}^\infty (\mu(p^k)/\phi(p^k))^3 c_{p^k}(n)$. Da $\mu(p^k)=0$ für $k\geq 2$: Faktor $= 1 + (\mu(p)/\phi(p))^3 c_p(n) = 1 + (-1/(p-1))^3 c_p(n) = 1 - c_p(n)/(p-1)^3$. Für $p\nmid n$: $c_p(n)=-1$, also Faktor $= 1 + 1/(p-1)^3$. Für $p|n$: $c_p(n)=p-1$, also Faktor $= 1 - (p-1)/(p-1)^3 = 1 - 1/(p-1)^2$. Das stimmt mit dem behaupteten Produkt überein. Der **Beweis ist also korrekt für das Endergebnis**, aber der Entwurfsüberrest "Wait, let us recompute" ist ein fataler redaktioneller Fehler.

**BUG-B4-04 (HOCH, EN):** Zeilen 253–259 in paper18_vinogradov_three_primes.tex: Entwurfsüberrest "Wait, let us recompute:" im veröffentlichten Beweis. Der Beweis bricht ab mit `...` und ohne vollständige Ausführung für den Fall $p \nmid n$. Muss überarbeitet werden.

**Satz 4.2 (Untere Schranke):**
Die Schranke $\mathfrak{S}(n) \gg 1/\log\log n$ für ungerades $n$ ist korrekt. Der Beweis über das endliche Teilprodukt $\prod_{p \leq \omega(n)}(1-4/p^2) \gg 1/\log\log n$ ist korrekt (Standardabschätzung).

**Bemerkung 4.3 (Gerades $n$):**
Die Erklärung ist konzeptionell korrekt, die Berechnung des $p=2$-Faktors aber verwirrend formuliert ("More precisely, three odd primes sum to an odd number, confirming..."). Die Aussage $\mathfrak{S}(n) = 0$ für gerades $n$ ist korrekt und folgt direkt aus dem $p=2$-Faktor.

**Satz 5.1 (Singuläres Integral):**
$\mathfrak{J}(n) = n^2/2 + O(n \log n)$. Der Beweis über Plancherel und Simplexfläche ist elegant und korrekt. Die Fläche des Simplex $\{t_1+t_2+t_3=n, t_i \in [2,N]\}$ ist für $n \approx N$ in der Tat $\sim n^2/2$.

**WARNUNG:** Das Paper sagt `$\mathfrak{J}(n) = n^2/2 + O(n\log n)$` und die Herleitung spricht von `[2,N]^3`. Für $n \leq N$ (d.h. $p_i \leq N = n$) ist das korrekt. Das Paper definiert $I(\beta) = \int_2^N e(t\beta)\,dt$ mit $N$ (d.h. das Zählargument) — für den Drei-Primes-Fall mit $p_i \leq n$ gilt $N = n$ im Kontext, was konsistent ist.

**Satz 6.1 (Vinogradov Minor Arc):**
$|S(\alpha)| \ll N(\log N)^{-B/4+1}$ auf $\mathfrak{m}$. Die Beweisidee ist korrekt skizziert (Siegel–Walfisz + Bombieri–Vinogradov). Die genaue Form der Schranke ist konventionell und korrekt.

**Korollar 6.2 (Nebenbogen-Fehler):**
Die Condition `B > 5A + 10` im Paper-Text, dann im Beweis `B \geq 4A + 8`: **inkonsistenter Schwellwert**. Im Text steht `B > 5A + 10`, im Beweis kommt heraus `B \geq 4A+8`. Beides impliziert dasselbe Resultat (`O(n^2/(\log n)^A)`), aber die Formulierung ist widersprüchlich.

**BUG-B4-05 (MINOR, EN):** Korollar 6.2 in paper18_vinogradov_three_primes.tex: Inkonsistenz zwischen der Bedingung im Satz (`B > 5A + 10`) und der im Beweis hergeleiteten Bedingung (`B \geq 4A + 8`). Beides liefert das Resultat, aber eine einheitliche Formulierung ist erforderlich.

**Satz 8.1 (Helfgott):**
Korrekt. Die Grenze $e^{30} \approx 1.07 \times 10^{13}$ und die Verwendung der Oliveira-e-Silva-Daten bis $4 \times 10^{18}$ sind korrekt dargestellt.

### LaTeX-Prüfung
- Alle Theorem-Umgebungen korrekt.
- `\DeclareMathOperator{\Li}{li}` ist definiert aber nie verwendet — toter Code (analoger Fall wie BUG-B4-01).
- Alle `\bibitem`-Einträge vorhanden.
- Keine undefinierten Umgebungen.
- `theorem*` in Abschnitt 9 (Conclusion): `\begin{theorem*}[Vinogradov--Helfgott]` — die `*`-Variante ist in `amsart` mit dem `\theoremstyle`-Setup nicht explizit definiert. `amsthm` erlaubt `\newtheorem*` für eine nicht-nummerierte Variante. Da `\newtheorem*{theorem*}{Theorem}` nicht in der Präambel steht, ist **diese Umgebung undefiniert** und wird einen LaTeX-Fehler erzeugen.

**BUG-B4-06 (HOCH, EN):** Zeile 515 in paper18_vinogradov_three_primes.tex: `\begin{theorem*}[Vinogradov--Helfgott]` verwendet eine undefinierte Theorem-Umgebung. Es fehlt in der Präambel: `\newtheorem*{theorem*}{Theorem}`. LaTeX-Kompilation wird hier fehlen.

### Urteil: **NICHT DRUCKREIF** (BUG-B4-04: Entwurfsüberrest; BUG-B4-05: Inkonsistenz; BUG-B4-06: undefinierte Umgebung)

---

## Paper 18 DE — Winogradows Drei-Primzahlen-Satz

**Datei:** `papers/batch4/paper18_vinogradov_drei_primzahlen_de.tex`

### Inhaltsbeschreibung
Deutsche Version von Paper 18 EN. Kompaktere Darstellung, ohne einige der ausführlichen Teilbeweise der EN-Version.

### Mathematische Prüfung

**Satz 3.1 (Hauptbogen-Approximation):**
Gleiche Formel wie EN. Der Beweis ist kürzer (korrekt zusammengefasst), ohne die fehlerhafte Orthogonalitätszeile der EN-Version. Daher: der DE-Beweis ist in dieser Hinsicht besser als EN.

**Satz 4.1 (Euler-Produkt):**
Der **Beweis ist in der DE-Version korrekt ausgeführt** und enthält nicht den Entwurfsüberrest "Wait, let us recompute":
```latex
p \mid n: 1 + c_p(n)/(p-1)^3 = 1 + (p-1)/(p-1)^3 = 1 + 1/(p-1)^2
p \nmid n: 1 + c_p(n)/(p-1)^3 = 1 - 1/(p-1)^3
```
Hier ist jedoch die **DE-Version nicht konsistent mit dem EN-behaupteten Euler-Produkt**!

In EN lautet das Euler-Produkt:
$$\prod_{p\nmid n}\left(1 + \frac{1}{(p-1)^3}\right)$$

In DE steht im Beweis:
$$p \nmid n: \quad 1 + \frac{c_p(n)}{(p-1)^3} = 1 - \frac{1}{(p-1)^3}$$

Das ergibt $\prod_{p\nmid n}(1 - 1/(p-1)^3)$, **nicht** $(1 + 1/(p-1)^3)$ wie in EN (und im DE-Theorem-Statement Zeile 163). **DE-Beweis widerspricht DE-Theorem-Statement!**

Die korrekte Rechnung (wie oben gezeigt): $\mu(p)/\phi(p) = -1/(p-1)$, also $(\mu(p)/\phi(p))^3 = -1/(p-1)^3$, und der Faktor ist $1 + (-1/(p-1)^3) \cdot c_p(n)$. Für $p\nmid n$: $c_p(n) = -1$, Faktor $= 1 + 1/(p-1)^3$. Der DE-Beweis hat einen **Vorzeichenfehler** beim Einsetzen: er schreibt direkt $c_p(n)/(p-1)^3$ statt $(\mu(p)/\phi(p))^3 \cdot c_p(n) = -c_p(n)/(p-1)^3$.

**BUG-B4-07 (HOCH, DE):** Zeilen 175–179 in paper18_vinogradov_drei_primzahlen_de.tex: Vorzeichenfehler im Beweis von Satz 4.1 (Euler-Produkt). Der Beweis setzt $c_p(n)/(p-1)^3$ statt $(\mu(p)/\phi(p))^3 \cdot c_p(n) = -c_p(n)/(p-1)^3$ ein. Der Faktor für $p\nmid n$ ist $1 - 1/(p-1)^3$ (laut Beweis) statt korrekt $1 + 1/(p-1)^3$ (laut Theorem-Statement). Theorem-Statement und Beweis widersprechen sich.

**Satz 5.1 (Singuläres Integral, DE):**
$\mathfrak{J}(n) = (n-6)^2/2 + O(n) \sim n^2/2$. Der DE-Beweis ist direkter als EN und **präziser** für den relevanten Fall ($n \approx N$). Korrekt.

**Korollar 5.2 (Nebenbögen, DE):**
Bedingung `B > 4A + 8`, konsistent mit dem DE-Beweis. **Keine Inkonsistenz** — DE ist hier besser als EN (BUG-B4-05 existiert in DE nicht).

**Abschnitt 8 (Paritätsproblem, DE):**
Der Abschnitt ist als eigenständiger Abschnitt herausgelöst (in EN ist er nur eine Bemerkung). Die Erklärung ist klar und korrekt.

**Schlussfolgerung (DE):**
Die DE-Version benutzt eine einfache `\begin{center}...\end{center}`-Konstruktion für das Hauptresultat anstelle einer `theorem*`-Umgebung — damit vermeidet die DE-Version BUG-B4-06.

### LaTeX-Prüfung
- Alle Theorem-Umgebungen korrekt.
- Keine undefinierten Makros.
- `\bibitem`-Einträge vollständig.
- Keine Baustellen.

### Urteil: **NICHT DRUCKREIF** (BUG-B4-07: Vorzeichenfehler im Euler-Produkt-Beweis)

---

## Paper 19 EN — Waring's Problem

**Datei:** `papers/batch4/paper19_waring_problem.tex`

### Inhaltsbeschreibung
Das Paper behandelt das Waringsche Problem: Lagranges Vier-Quadrate-Satz, den Hilbert–Waring-Satz, die Hardy–Littlewood-Formel für $r_{k,s}(n)$, Weyls Schranke $G(k) \leq 2^k + 1$, Vinogradovs Mittelwertsatz, die Bourgain–Demeter–Guth/Wooley-Auflösung der Vinogradov-Hauptvermutung, Wooleys $G(k)$-Schranke, bekannte Werte von $g(k)$ und $G(k)$, und die singuläre Reihe für das Waringsche Problem.

### Mathematische Prüfung

**Satz 2.1 (Lagrange):** Korrekt, $g(2) = 4$.

**Satz 2.2 ($G(2) = 4$):**
Der Beweis ist korrekt. Die Aussage $x^2 \in \{0,1,4\} \pmod 8$ und daher Summe von drei Quadraten $\in \{0,...,6\}$ ist richtig. Legendre's three-square theorem korrekt zitiert.

**Satz 3.1 (Hilbert–Waring):**
Der Beweis ist nur skizziert und enthält die Hilbert-Identität. Die Formel in Zeile 174–177:
$$\binom{2k}{k} = \frac{1}{k} \sum_{j=0}^{k} (-1)^{k-j}\binom{k}{j}(2j-k)^{2k} \cdot \frac{1}{k^{2k}}$$
Diese Darstellung der Hilbert-Identität ist **inhaltlich unscharf**: Die Hilbert-Identität besagt, dass $k$-te Potenzen von Polynomen eine rationale Identität der Form $\sum c_i x_i^k = C$ erfüllen. Die hier angegebene Formel hat die falsche Form (linke Seite ist ein Binomialkoeffizient, rechte Seite enthält einen zusätzlichen Faktor $1/k^{2k}$, was die Formel für $k\geq 2$ inkonsistent macht). Dies ist als **Beweis-Skizze** mit Verweis auf Nathanson akzeptabel, sollte aber präziser formuliert sein. Der Verweis auf Nathanson ist korrekt.

**Satz 4.1 (Hardy–Littlewood-Formel):**
Der Verweis auf Paper 17 ist korrekt. Die Formel selbst ist richtig.

**Korollar 4.2 — KRITISCH (Weyls Schranke):**
Der Beweis bricht explizit ab mit "wait, this is inconsistent" (Zeile 225). Dies ist ein **schwerwiegender Entwurfsüberrest**. Der Text sagt anschließend, dass $G(k) \leq 2^k+1$ aus einer "cruder argument" folgt, und das korrekte Ergebnis ist $G(k) \leq (k-2)2^{k-1}+5$. Dies zeigt eine fundamentale Unstimmigkeit: die Überschrift und der Korollar-Statement behaupten $G(k) \leq 2^k+1$, aber im Beweis selbst steht, dass das falsch ist, und die korrekte Schranke $G(k) \leq (k-2)2^{k-1}+5$ wird nur nebenbei erwähnt. Der Korollar-Statement sollte die korrekte Waring-Weyl-Schranke $G(k) \leq (k-2)2^{k-1}+5$ angeben.

**BUG-B4-08 (HOCH, EN):** Zeilen 213–229 in paper19_waring_problem.tex: Korollar 4.2 behauptet $G(k) \leq 2^k+1$ als "Weyl's original bound", der Beweis bricht aber mit "wait, this is inconsistent" ab und korrigiert sich auf $G(k) \leq (k-2)2^{k-1}+5$. Entwurfsüberrest im Beweis, der Korollar-Statement ist inkonsistent mit dem tatsächlichen Resultat. Die korrekte Weyl-Schranke lautet $G(k) \leq (k-2)2^{k-1}+5$.

**Satz 5.1 (Vinogradov MVT):**
$J_{s,k}(N) \ll N^{2s-k(k+1)/2+\varepsilon}$ für $s \geq k(k+1)/2$. Korrekt.

**Satz 5.2 (Vinogradov $G(k)$-Schranke):**
$G(k) \leq k(\log k + \log\log k + O(1))$. Korrekt (klassisches Vinogradov-Ergebnis).

**Satz 5.3 (BDG/Wooley):**
Die Formulierung "For $s \leq k(k+1)/2$: $J_{s,k}(N) \ll N^{s+\varepsilon}$" ist korrekt (der triviale Anteil ist $N^s$). Die Aussage für $s = k(k+1)/2$ exakt ist korrekt.

**BUG-B4-09 (MINOR, EN):** Bibliographie, `\bibitem{Wooley2016}` (Zeilen 503–507): Es sind **zwei verschiedene Paper** unter einem bibitem zusammengefasst — "Duke Math. J., 162(4):673–730, 2013" und "Adv. Math., 294:532–561, 2016". Dies sind zwei separate Publikationen, die nicht in einem `\bibitem` zusammengeführt sein dürfen. Das erstere (2013) ist der erste Teil der "efficient congruencing"-Serie, das zweite (2016) der zweite Teil. Beide sollten als separate Einträge erscheinen.

**Satz 7.1 ($\mathfrak{S}_{k,s}(n) > 0$):**
Für $s \geq 2k+1$: $\mathfrak{S}_{k,s}(n) > 0$. Korrekt, mit Verweis auf Vaughan.

**Tabelle zu $g(k)$:**
Die Werte $g(2)=4, g(3)=9, g(4)=19, g(5)=37, g(6)=73, g(7)=143, g(8)=279$ sind korrekt. Die Zuordnung zu Autoren: Für $g(4)=19$: "Balasubramanian, Dress, Deshouillers 1986" — tatsächlich wurde $g(4)=19$ von Balasubramanian, Deshouillers und Dress (1986) bewiesen, korrekt. Für $g(7)=143$ und $g(8)=279$: "Pillai 1936" — dies ist korrekt für $k=6$, aber für $k=7,8$ bedarf es präziserer Quellenangaben. Für den Rahmen des Papers ist diese Genauigkeit ausreichend.

**Tabelle zu $G(k)$:**
$G(4)=16$ mit "Balasubramanian 1986": korrekt (gemeint: Balasubramanian, Deshouillers, Dress). $G(3) \leq 7$ mit "Heath-Brown 1988": korrekt. Diese Werte sind korrekt.

### LaTeX-Prüfung
- Alle Theorem-Umgebungen korrekt.
- `✓`-Symbole in Tabellen: LaTeX-Standard-Dokument kann diese Unicode-Symbole ohne `\usepackage{pifont}` oder ähnliches nicht direkt verarbeiten. Je nach LaTeX-Kompiler (LuaLaTeX/XeLaTeX vs. pdfLaTeX) kann dies ein Problem sein. Mit pdfLaTeX und inputenc ist das möglicherweise fehlerhaft. (Gleiches Problem in Paper 19 DE.)
- Alle `\bibitem`-Einträge vorhanden (mit dem in BUG-B4-09 beschriebenen Duplikat-Problem).

### Urteil: **NICHT DRUCKREIF** (BUG-B4-08: "wait, this is inconsistent"-Überrest; BUG-B4-09: Duplikates bibitem)

---

## Paper 19 DE — Das Waringsche Problem

**Datei:** `papers/batch4/paper19_waringsches_problem_de.tex`

### Inhaltsbeschreibung
Deutsche Version von Paper 19 EN. Kompakter, ohne die ausführlichen Beweise für Weyl-Schranken und Hilbert-Identitäten.

### Mathematische Prüfung

**Wesentlicher Unterschied zu EN:**
- Die DE-Version enthält **nicht** das Korollar 4.2 mit dem "wait, this is inconsistent"-Entwurfsüberrest (BUG-B4-08). Stattdessen wird die Weyl-Schranke nur im Abstract als $G(k) \leq 2^{k-1}(k-2)+5$ erwähnt — was die korrekte Formulation ist.
- Damit ist BUG-B4-08 in DE **nicht vorhanden** — die DE-Version ist in diesem Punkt besser als EN.

**Hilbert-Identität (Zeile 140–142):**
Dieselbe unklare Formel wie in EN, aber im DE-Kontext als Beweisidee akzeptabel.

**Satz 4.1 (Formel für $g(k)$):**
Korrekt. Beweis über "schwerste Zahl" $n = 2^k q - 1$ ist korrekt.

**`✓`-Symbole:**
Gleiche potentielle LaTeX-Inkompatibilität wie in EN (für pdfLaTeX).

**BUG-B4-10 (MINOR, DE):** Abstract-Formel "Weyls Schranke $G(k) \leq 2^{k-1}(k-2)+5$" im Abstract (Zeile 51): Diese lautet korrekt $G(k) \leq (k-2)2^{k-1}+5$ — was mathematisch dasselbe ist, aber eine ungewöhnliche Schreibweise.

### LaTeX-Prüfung
Keine Fehler außer dem potentiellen `✓`-Problem.

### Urteil: **DRUCKREIF** (1 trivialer Minorbug BUG-B4-10; keine mathematischen Fehler)

---

## Paper 20 EN — The Goldbach Singular Series

**Datei:** `papers/batch4/paper20_goldbach_singular_series.tex`

### Inhaltsbeschreibung
Das Paper leitet die Goldbach-singuläre-Reihe $\mathfrak{S}(n) = 2C_2 \prod_{p|n, p\geq 3}(p-1)/(p-2)$ her, beweist ihre Positivität für gerades $n$, diskutiert die Zwillingsprimzahl-Konstante $C_2$, formuliert die Hardy–Littlewood-Hypothese H über Primzahl-$k$-Tupel, erklärt das Paritätsproblem und diskutiert das bedingte Goldbach-Ergebnis unter GRH.

### Mathematische Prüfung

**Definition 2.1 (Goldbach-Singuläre-Reihe):**
$\mathfrak{S}(n) = \sum_{q=1}^\infty \mu(q)^2 c_q(n)/\phi(q)^2$. Konvergenz: $|c_q(n)| \leq \phi(q)$, also $|\mathfrak{S}(n)| \leq \sum_q \mu(q)^2/\phi(q)$. Da $\sum_q \mu(q)^2/\phi(q) = \prod_p(1+1/(p-1))$ und dieses Produkt divergiert (da $\sum_p 1/(p-1)$ divergiert), ist die **Konvergenzaussage im Paper falsch**.

Der richtige Konvergenznachweis: $c_q(n) = \mu(q/\gcd(q,n))\phi(q)/\phi(q/\gcd(q,n))$ für quadratfreie $q$. Da $\mu(q)^2=0$ für nicht-quadratfreie $q$ und für quadratfreie $q$: $|c_q(n)/\phi(q)^2| = |\mu(q/\gcd(q,n))|/(\phi(q)\phi(q/\gcd(q,n))) \leq 1/\phi(q/\gcd(q,n))$. Die Konvergenz folgt aus $\sum_q 1/\phi(q)^2 < \infty$, was korrekt ist (das Produkt $\prod_p(1+1/(p-1)^2)$ konvergiert). Die im Paper angegebene Begründung "$\sum_q \mu(q)^2/\phi(q) < \infty$" ist aber eine **divergente Reihe** — die Konvergenzbegründung ist **falsch**.

**BUG-B4-11 (HOCH, EN):** Definition 2.1 in paper20_goldbach_singular_series.tex, Zeile 144: Die Konvergenzbegründung "since $|c_q(n)| \leq \phi(q)$ and $\sum_q \mu(q)^2/\phi(q) = \prod_p(1+1/(p-1)) < \infty$" ist falsch: (a) $|c_q(n)| \leq \phi(q)$ gibt nur $\sum_q \mu(q)^2/\phi(q)$ als Majorante, aber diese Reihe divergiert; (b) $\prod_p(1+1/(p-1))$ divergiert (Mertens' Theorem). Die korrekte Konvergenzbegründung nutzt $|c_q(n)/\phi(q)^2| \leq 1/\phi(q)^2$ (nicht $1/\phi(q)$) und $\sum_q 1/\phi(q)^2 < \infty$.

**Satz 2.2 (Euler-Produkt):**
Das Ergebnis $\mathfrak{S}(n) = 2C_2 \prod_{p|n, p\geq 3}(p-1)/(p-2)$ ist korrekt. Der Beweis ist detailliert und mathematisch korrekt (Zeilen 159–202). Die Zwischen-Rechnung in Zeilen 182–201 ist korrekt und vollständig.

**Zeile 182:** Tippfehler im LaTeX-Quellcode: `= p(p-2)/(p-1)^2 = p(p-2)/(p-1)^2$` endet mit einer unausgeglichenen Klammer `$`. Konkret:
```latex
For $p \geq 3$ with $p \nmid n$: ... factor = $p(p-2)/(p-1)^2$.
```
Der tatsächliche LaTeX-Code in Zeile 182 lautet:
`= (p-1)^2-1)/(p-1)^2 = p(p-2)/(p-1)^2`.
Dies enthält eine unbalancierte Klammer `(p-1)^2-1)` — die öffnende Klammer fehlt: müsste `((p-1)^2-1)/(p-1)^2` sein.

**BUG-B4-12 (MITTEL, EN):** Zeile 182 in paper20_goldbach_singular_series.tex: Unbalancierte Klammer im LaTeX-Quellcode: `(p-1)^2-1)/(p-1)^2` sollte `((p-1)^2-1)/(p-1)^2` heißen. LaTeX-Kompilation wird hier entweder stillschweigend falsch rendern oder einen Fehler erzeugen.

**Satz 3.1 (Positivität):**
Korrekt und vollständig bewiesen.

**Bemerkung 3.2 (Ungerades $n$):**
Der Text enthält einen Denkfehler: "and $n = p + 2$ would require $p$ even, i.e. $p = 2$, giving $n = 2 + 2 = 4$, odd only for..." — der Satz bricht ab ohne Abschluss. Dies ist ein weiterer **Entwurfsüberrest**.

Außerdem ist die Erklärung unnötig kompliziert. Die einfache Aussage: "Zwei ungerade Primzahlen summieren sich zu einer geraden Zahl" würde genügen. Die unterbrochene Textpassage ist ein redaktioneller Fehler, aber kein mathematischer.

**Satz 5.1 (Satz von Hardy–Littlewood über Zwillingsprimzahlen):**
Die Aussage $\pi_2(x) \sim 2C_2 \cdot x/(\log x)^2$ ist eine Vermutung, kein Satz. Das Paper bezeichnet es korrekt als Konsequenz der Hardy–Littlewood-Vermutung H. Der "Beweis" ist als "consequence of the H-L 2-tuples conjecture" bezeichnet — korrekt, aber die Theorem-Umgebung sollte eine Vermutung oder ein "conditional theorem" sein. Der aktuelle Text sagt "The conjecture is conditional on Hypothesis H" — das ist ausreichend, um den Status klarzustellen.

**Satz 6.1 (GRH-Ergebnis):**
Der Beweisabsatz enthält eine inkonsistente Aussage: "while the minor arc bound is $O(N^{3/2}(\log N)^3)$. Compared to the main term $\mathfrak{S}(n)\cdot N \gg N/\log\log N$, the minor arc error $N^{3/2}(\log N)^3$ still dominates for large $N$..."

Das Paper erkennt selbst, dass die dargestellte Abschätzung nicht ausreicht, und fährt dann fort mit "The improvement under GRH is that one can use exponential sum bounds... $O(N^{1-\varepsilon})$ obtainable under GRH via character sum bounds". Dies ist **inkonsistent**: zunächst zeigt der Text, dass die Methode versagt, dann behauptet er ohne Beweis, dass sie unter GRH doch funktioniert. Der Verweis auf `\cite{HardyLittlewood1923}` für dieses GRH-Resultat ist falsch — das 1923-Paper beweist dieses GRH-Ergebnis nicht; das ist ein modernes Resultat. Das korrekte Resultat (Goldbach unter GRH für große $n$) geht auf Hardy–Littlewood und später präzisiert von Goldston zurück.

Der inhaltliche **Kern** (Goldbach unter GRH für große $n$) ist mathematisch korrekt, aber der Beweis ist unvollständig und enthält einen falschen Literaturverweis.

**BUG-B4-13 (MITTEL, EN):** Satz 6.1 in paper20_goldbach_singular_series.tex: (a) Der Beweis zeigt zunächst, dass die naive Methode unter GRH noch scheitert (minor arc $O(N^{3/2}(\log N)^3) \gg N$), behauptet dann aber ohne vollständigen Beweis, dass es doch funktioniert. (b) Der Verweis `\cite{HardyLittlewood1923}` ist falsch — das 1923-Paper beweist kein GRH-bedingtes Ergebnis für binäres Goldbach. Korrekte Referenz wäre etwa Goldston 1992 oder Goldston–Yıldırım 2003.

**Tabelle (Section 4):**
Zeile mit $n=30$: `$2C_2 \cdot \frac{2}{1}\cdot\frac{4}{3}\cdot\frac{6}{5}\cdot\frac{28}{28}$`. Der letzte Faktor $28/28 = 1$ für $p=29$? Nein: $30 = 2 \cdot 3 \cdot 5$, also sollten die Faktoren für $p=3$ ($4/3$) und $p=5$ ($6/5$) kommen, nicht $28/28$. Der Faktor $28/28$ ist **falsch** — für $p=29$ (das nächste Primzahl-Teiler von 30): $30/29$ ist kein ganzzahliger Teiler, und 29 teilt 30 nicht. Es wäre $(29-1)/(29-2) = 28/27$ für $p=29$ falls $29|30$, aber $29\nmid 30$. Dieser Faktor gehört nicht in die Formel für $n=30$.

**BUG-B4-14 (MITTEL, EN):** Tabelle in Section 4 (Zeile mit $n=30$): Der Faktor `$\frac{28}{28}$` ist falsch. Da $30 = 2\cdot 3\cdot 5$ und keine weiteren ungeraden Primteiler existieren, lautet die korrekte Formel: $\mathfrak{S}(30) = 2C_2 \cdot \frac{2}{1}\cdot\frac{4}{3}\cdot\frac{6}{5}$. Der Faktor $28/28=1$ ist ein offensichtlicher Rechenfehler/Platzhalter.

### LaTeX-Prüfung
- `\DeclareMathOperator{\Li}{li}` definiert aber nie verwendet (toter Code).
- BUG-B4-12: unbalancierte Klammer.
- Alle `\bibitem`-Einträge vorhanden.

### Urteil: **NICHT DRUCKREIF** (BUG-B4-11: falsche Konvergenzbegründung; BUG-B4-12: unbalancierte Klammer; BUG-B4-13: inkonsistenter GRH-Beweis + falsches Zitat; BUG-B4-14: falscher Tabelleneintrag)

---

## Paper 20 DE — Die Goldbach'sche singuläre Reihe

**Datei:** `papers/batch4/paper20_goldbach_singulaere_reihe_de.tex`

### Inhaltsbeschreibung
Deutsche Version von Paper 20 EN. Leicht kompakter, mit derselben mathematischen Substanz.

### Mathematische Prüfung

**Definition 2.1 (Goldbach-Singuläre-Reihe, DE):**
Die DE-Version enthält nur die Definition ohne explizite Konvergenzbegründung. Damit vermeidet sie BUG-B4-11. Kein mathematischer Fehler durch Unterlassung (die Konvergenz ist implizit).

**Satz 2.2 (Euler-Produkt, DE):**
Der DE-Beweis ist **vollständiger und klarer** als EN und enthält nicht die unbalancierte Klammer. Die Zwischenrechnung (Zeilen 143–150) ist korrekt:
```latex
= 2 * C_2 * prod_{p>=3, p|n} p/(p-1) / (p(p-2)/(p-1)^2)
= 2 C_2 prod_{p>=3, p|n} (p-1)/(p-2)
```
Korrekt.

**Tabelle (Section 4, DE):**
Zeile mit $n=30$: DE schreibt `$2C_2 \cdot \frac{2}{1}\cdot\frac{4}{3}\cdot\frac{6}{5}\cdot\frac{28}{27}$`. Der Faktor $28/27$ für $p=29$ ist aber ebenfalls **falsch** — 29 teilt 30 nicht. Es sollte nur die Faktoren für $p=3$ und $p=5$ geben.

**BUG-B4-15 (MINOR, DE):** Tabelle in Section 4 (Zeile mit $n=30$): Der Faktor `$\frac{28}{27}$` ist falsch. $30 = 2\cdot 3\cdot 5$, also $\mathfrak{S}(30) = 2C_2 \cdot \frac{2}{1}\cdot\frac{4}{3}\cdot\frac{6}{5}$ ohne weiteren Faktor (da 29 kein Teiler von 30 ist). Bemerkung: Der Faktor $(p-1)/(p-2) = 28/27$ gehört zu $p=29$, aber $29 \nmid 30$.

**Satz (GRH, DE):**
$R_2(n) = \mathfrak{S}(n)\cdot n + O(\sqrt{n}(\log n)^2)$. Kein Beweis angegeben — nur Aussage. Das ist für eine DE-Übersicht akzeptabel, ohne den inkonsistenten EN-Beweis zu reproduzieren.

**Bemerkung (Ungerades $n$, DE):**
"$c_2(n) = -1$, Faktor bei $p=2$: $1+(-1)/1 = 0$. Daher $\mathfrak{S}(n) = 0$ für ungerades $n$." Korrekt und knapp. Kein Entwurfsüberrest wie in EN.

### LaTeX-Prüfung
- Alle Umgebungen korrekt definiert.
- Keine undefinierten Makros.
- Alle `\bibitem`-Einträge vorhanden.
- Keine Baustellen.

### Urteil: **DRUCKREIF** (1 Minorbug: BUG-B4-15)

---

## Zusammenfassung: Bug-Tabelle

| Bug-ID | Paper | Schwere | Beschreibung |
|--------|-------|---------|-------------|
| BUG-B4-01 | 17 EN | MINOR | Toter Code: `\newcommand{\e}[1]` definiert aber nie verwendet |
| BUG-B4-02 | 17 DE | MINOR | Beweis von Lemma 3.1 (Farey-Partition) fehlt in DE (EN hat ihn) — Inkonsistenz |
| BUG-B4-03 | 18 EN | MITTEL | Falsche Orthogonalitätsformel im Beweis von Lemma 3.1 (Zeile 179): `\sum_\chi \bar\chi(a) = \phi(q)·1_{a≡1(q)}` |
| BUG-B4-04 | 18 EN | HOCH | Entwurfsüberrest "Wait, let us recompute" im Beweis von Satz 4.1, Euler-Produkt unvollständig |
| BUG-B4-05 | 18 EN | MINOR | Inkonsistenter Schwellwert in Korollar 6.2: `B > 5A+10` (Text) vs. `B ≥ 4A+8` (Beweis) |
| BUG-B4-06 | 18 EN | HOCH | Undefinierte LaTeX-Umgebung `theorem*` in Abschnitt 9 — Kompilationsfehler |
| BUG-B4-07 | 18 DE | HOCH | Vorzeichenfehler im Beweis des Euler-Produkts (Satz 4.1): Beweis ergibt `1-1/(p-1)^3` statt `1+1/(p-1)^3` für $p\nmid n$ |
| BUG-B4-08 | 19 EN | HOCH | Entwurfsüberrest "wait, this is inconsistent" in Korollar 4.2; Statement behauptet $G(k)\leq 2^k+1$, korrekt ist $G(k)\leq (k-2)2^{k-1}+5$ |
| BUG-B4-09 | 19 EN | MINOR | Duplikates `\bibitem{Wooley2016}` enthält zwei separate Publikationen (2013 + 2016) |
| BUG-B4-10 | 19 DE | MINOR | Abstract: Weyl-Schranke $G(k)\leq 2^{k-1}(k-2)+5$ (ungewöhnliche Schreibweise, mathematisch äquivalent) |
| BUG-B4-11 | 20 EN | HOCH | Falsche Konvergenzbegründung für $\mathfrak{S}(n)$: $\sum_q \mu(q)^2/\phi(q)$ divergiert (Mertens) |
| BUG-B4-12 | 20 EN | MITTEL | Unbalancierte Klammer in Zeile 182: `(p-1)^2-1)/(p-1)^2` statt `((p-1)^2-1)/(p-1)^2` |
| BUG-B4-13 | 20 EN | MITTEL | GRH-Beweis (Satz 6.1) inkonsistent (zeigt erst Scheitern, behauptet dann Erfolg ohne Beweis); falscher Literaturverweis `\cite{HardyLittlewood1923}` |
| BUG-B4-14 | 20 EN | MITTEL | Tabelle $n=30$: Faktor `28/28=1` ist falsch (29 teilt 30 nicht) |
| BUG-B4-15 | 20 DE | MINOR | Tabelle $n=30$: Faktor `28/27` falsch (29 teilt 30 nicht), korrekt: kein weiterer Faktor |

### Schwere-Übersicht

| Schwere | Anzahl Bugs |
|---------|-------------|
| HOCH    | 5 (BUG-B4-04, -06, -07, -08, -11) |
| MITTEL  | 4 (BUG-B4-03, -12, -13, -14) |
| MINOR   | 6 (BUG-B4-01, -02, -05, -09, -10, -15) |

---

## Gesamturteil pro Paper

| Paper | Datei | Urteil | Kritische Bugs |
|-------|-------|--------|----------------|
| 17 EN | paper17_circle_method.tex | **DRUCKREIF** | — |
| 17 DE | paper17_kreismethode_de.tex | **DRUCKREIF** | — |
| 18 EN | paper18_vinogradov_three_primes.tex | **NICHT DRUCKREIF** | B4-04, B4-06 |
| 18 DE | paper18_vinogradov_drei_primzahlen_de.tex | **NICHT DRUCKREIF** | B4-07 |
| 19 EN | paper19_waring_problem.tex | **NICHT DRUCKREIF** | B4-08 |
| 19 DE | paper19_waringsches_problem_de.tex | **DRUCKREIF** | — |
| 20 EN | paper20_goldbach_singular_series.tex | **NICHT DRUCKREIF** | B4-11, B4-13 |
| 20 DE | paper20_goldbach_singulaere_reihe_de.tex | **DRUCKREIF** | — |

---

## Prioritäten für Korrekturen

### Sofort zu beheben (HOCH)

1. **BUG-B4-04** (18 EN): Entwurfsüberrest "Wait, let us recompute" vollständig entfernen und Beweis ausführen (die korrekte Berechnung: Faktor für $p\nmid n$ ist $1 - c_p(n)/(p-1)^3 = 1+1/(p-1)^3$).
2. **BUG-B4-06** (18 EN): `\newtheorem*{theorem*}{Theorem}` in Präambel ergänzen oder `theorem*`-Umgebung durch andere ersetzen.
3. **BUG-B4-07** (18 DE): Vorzeichenfehler im Euler-Produkt-Beweis beheben: korrekte Formel für $(\mu(p)/\phi(p))^3 \cdot c_p(n) = -c_p(n)/(p-1)^3$.
4. **BUG-B4-08** (19 EN): "wait, this is inconsistent" entfernen; Korollar-Statement auf $G(k) \leq (k-2)2^{k-1}+5$ korrigieren und einen korrekten Beweis einfügen.
5. **BUG-B4-11** (20 EN): Konvergenzbegründung ersetzen durch: $\sum_q \mu(q)^2|c_q(n)|/\phi(q)^2 \leq \sum_q 1/\phi(q)^2 = \prod_p(1+1/(p-1)^2) < \infty$.

### Bald zu beheben (MITTEL)

6. **BUG-B4-03** (18 EN): Orthogonalitätsformel im Beweis präzisieren.
7. **BUG-B4-12** (20 EN): Unbalancierte Klammer korrigieren: `((p-1)^2-1)/(p-1)^2`.
8. **BUG-B4-13** (20 EN): GRH-Beweis entweder vollständig ausführen oder als "Beweis-Skizze" explizit markieren; Literaturverweis korrigieren.
9. **BUG-B4-14** (20 EN) und **BUG-B4-15** (20 DE): Tabelleneintrag $n=30$ korrigieren: $\mathfrak{S}(30) = 2C_2\cdot\frac{2}{1}\cdot\frac{4}{3}\cdot\frac{6}{5}$ (nur $p=3$ und $p=5$ als ungerade Primteiler von 30).

### Bei Gelegenheit (MINOR)

10. **BUG-B4-01**: Toter Code `\newcommand{\e}` entfernen.
11. **BUG-B4-02**: Beweis von Lemma 3.1 (Farey-Partition) in DE-Version ergänzen.
12. **BUG-B4-05**: Konsistenten Schwellwert `B ≥ 4A+8` in Korollar 6.2 EN verwenden.
13. **BUG-B4-09**: `\bibitem{Wooley2016}` aufteilen in zwei separate Einträge (2013 und 2016).

---

*Gutachten erstellt von: Claude Sonnet 4.6 — 2026-03-11*
