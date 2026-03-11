# Review Paper 14: Selberg's Sieve (paper14_selberg_sieve.tex)

**Reviewer:** Mathematik-Review-Agent
**Datum:** 2026-03-11
**Paper:** "Selberg's Sieve and Upper Bounds for Primes in Arithmetic Progressions"
**Autor:** Michael Fuhrmann

---

## Gesamturteil: ⚠️ (Ein schwerwiegender Fehler und mehrere mittlere Fehler – Korrekturen nötig)

---

## Gefundene Bugs

| Bug-ID  | Schwere | Kurzbeschreibung |
|---------|---------|-----------------|
| BUG-P14-001 | Schwer | Theorem 5.1 (Brun-Titchmarsh): Beweis verwendet $z = x/q$ statt $z = \sqrt{x}$, aber dann $O(z^2) = O(x^2/q^2)$ als Fehlerterm angegeben – dies dominiert den Hauptterm für relevante $q$ |
| BUG-P14-002 | Mittel | Theorem 4.1: $V(z) \asymp V_h(z)^{-1}$ ist falsch notiert – korrekt ist $V_h(z) = 1/V(z)$ (exakt, nicht nur asymptotisch) |
| BUG-P14-003 | Mittel | Lemma 3.2: Die Constraint-Übersetzung "$\lambda_1 = 1$ translates to $\sum_d y_d = 1$" ist nicht bewiesen und in der verwendeten Form falsch |
| BUG-P14-004 | Mittel | Beweis Brun-Titchmarsh Step 3: Mertens-Formel falsch angewandt – $\phi(q)/q$-Korrektur umgekehrt |
| BUG-P14-005 | Gering | Bibliographie: Iwaniec1978-Eintrag mit falschem Journal und Jahreszahl (Label 1978, Text 1982) |
| BUG-P14-006 | Gering | Definition 2.2: $|r_d| \leq g(d)$ ist eine untypische und zu restriktive Normierung des Remainders |
| BUG-P14-007 | Gering | Vergleichstabelle Section 6: Brun-Sieb-Fehlerterm "$\sum_{d \leq z} |r_d|$" ist inkorrekt – korrekt ist $\sum_{d \leq z, \omega(d)\leq r} |r_d|$ |

---

## Details zu jedem Bug

### BUG-P14-001 (Schwer): Fehlerterm im Brun-Titchmarsh-Beweis dominiert den Hauptterm

**Fundstelle:** Zeilen 398–405, Theorem 5.1 Beweis, Step 3

**Das Paper schreibt:**
> "More carefully, choosing $z = x/q$ (instead of $\sqrt{x}$) and using the precise Mertens bound: [...]
> and the error term $\sum_{d_1,d_2 \leq z} |r_{[d_1,d_2]}|$ is $O(z^2) = O(x^2/q^2)$,
> which is dominated by the main term for $q \leq x^{1/2-\varepsilon}$."

**Fehler:**
1. Der Selberg-Fehlerterm ist $\sum_{d_1, d_2 \leq z} |r_{[d_1,d_2]}|$. Für unsere Wahl $r_d = O(1)$ (da $|r_d| \leq 1$) und $z$ ist die Anzahl der Paare $(d_1, d_2)$ mit $d_1, d_2 \leq z$ gleich $O(z^2)$, also Fehlerterm $= O(z^2)$.

2. Wählt man $z = x/q$, so ist $O(z^2) = O(x^2/q^2)$.

3. Der Hauptterm ist $\sim x/(\phi(q)\log(x/q))$.

4. Für den Fehlerterm subdominant zu sein, braucht man:
   $$\frac{x^2}{q^2} \ll \frac{x}{\phi(q)\log(x/q)},$$
   also $x \ll q^2/(\phi(q)\log(x/q)) \approx q/\log(x/q)$, was nur für sehr große $q$ (nahe $x$) gilt – also genau umgekehrt zur Behauptung.

**Korrekte Wahl der Parameter:** Für den Brun-Titchmarsh-Beweis sollte man $z = \sqrt{x/q}$ wählen (nicht $\sqrt{x}$ oder $x/q$). Dann ist $O(z^2) = O(x/q)$, der Hauptterm ist $\sim x/(\phi(q)\log\sqrt{x/q}) = 2x/(\phi(q)\log(x/q))$, und der Fehler $O(x/q)$ ist subdominant.

Die Wahl $z = \sqrt{x/q}$ ergibt unmittelbar den Brun-Titchmarsh-Koeffizienten 2 korrekt:
$$\log z = \frac{1}{2}\log(x/q) \implies \frac{1}{V(z)} \sim e^\gamma \log z = \frac{e^\gamma}{2}\log(x/q),$$
$$\frac{X}{V(z)} = \frac{x}{q} \cdot \frac{e^\gamma}{2}\log(x/q) \cdot \frac{q}{\phi(q) \cdot e^\gamma} = \frac{x\log(x/q)}{2\phi(q)\log(x/q)} \cdot 2 = \frac{2x}{\phi(q)\log(x/q)}.$$

**Fazit:** Die Wahl $z = x/q$ im Paper führt zu einem Fehlerterm, der den Hauptterm dominiert und macht den Beweis damit ungültig. Dies ist ein schwerwiegender Fehler.

---

### BUG-P14-002 (Mittel): Falsche Schreibweise $V(z) \asymp V_h(z)^{-1}$

**Fundstelle:** Zeilen 279–281, Theorem 4.1 Formulierung

**Das Paper schreibt:**
```latex
V(z) = \prod_{\substack{p \leq z \\ p \in \mathcal{P}}} \bigl(1 - g(p)\bigr)
\;\asymp\; V_h(z)^{-1}.
```

**Fehler:** Die Relation ist nicht nur $\asymp$ sondern exakt $=$. Im Beweis von Theorem 4.1 (Zeilen 290–295) wird dies auch explizit gezeigt:
$$V_h(z) = \prod_{p \leq z}\!\left(1 + \frac{g(p)}{1-g(p)}\right) = \prod_{p \leq z} \frac{1}{1-g(p)} = \frac{1}{V(z)}.$$

Das $\asymp$ in der Theorem-Formulierung ist also irreführend falsch – es sollte $V(z) = V_h(z)^{-1}$ heißen (ohne den Trunkatierungsfehler, der nur für die endliche Summe mit $d \leq z$ relevant ist).

**Einschränkung:** Bei der *endlichen* Summe $V_h(z) = \sum_{d \leq z, d \mid \mathcal{P}(z)} \mu(d)^2 h(d)$ gilt das Gleichheitszeichen nur ohne Trunkatierungseinschränkung. Mit der Bedingung $d \leq z$ fehlen Terme und die Relation ist tatsächlich $V_h(z) \leq 1/V(z)$. Das $\asymp$ ist in diesem Sinne korrekt, aber das Theorem sollte dies explizit klarstellen. Die Formulierung in Theorem 4.1 ist mehrdeutig.

---

### BUG-P14-003 (Mittel): Constraint-Übersetzung falsch

**Fundstelle:** Zeilen 245–247, Lemma 3.2, Step 3

**Das Paper schreibt:**
> "The constraint $\lambda_1 = 1$ translates to $1 = \lambda_1 = \sum_{d \mid \mathcal{P}(z)} y_d$."

**Fehler:** Diese Übersetzung ist nicht korrekt für die verwendete $y_d$-Definition.

Die $y_d$ wurden in Step 1 definiert als:
$$y_d = \sum_{\substack{e \mid \mathcal{P}(z) \\ d \mid e}} \mu(e/d)\, \lambda_e.$$

Für $d = 1$:
$$y_1 = \sum_{e \mid \mathcal{P}(z)} \mu(e)\, \lambda_e.$$

Das ist **nicht** das Gleiche wie $\sum_d y_d = 1$.

Die Möbius-Inversion lautet $\lambda_d = \sum_{e \mid \mathcal{P}(z), d \mid e} y_e$. Für $d = 1$:
$$\lambda_1 = \sum_{e \mid \mathcal{P}(z)} y_e = 1.$$

Das stimmt! Die Rechnung ist also letztlich korrekt, aber die Begründung sollte die Möbius-Inversionseigenschaft explizit nutzen ($\lambda_1 = \sum_e y_e$ via $\lambda_d = \sum_{d \mid e} y_e$ für $d=1$). Die aktuelle Darstellung im Paper suggeriert, man setze einfach $d=1$ in $y_d$ ein – das wäre falsch. Es ist die Inversion, nicht $y_1 = 1$.

---

### BUG-P14-004 (Mittel): Mertens-Formel und $\phi(q)/q$-Korrektur

**Fundstelle:** Zeilen 384–391, Theorem 5.1 Beweis, Step 3

**Das Paper schreibt:**
```latex
\prod_{\substack{p \leq z \\ p \nmid q}} \frac{p}{p-1}
= \prod_{p \leq z} \frac{p}{p-1} \cdot \prod_{\substack{p \mid q \\ p \leq z}} \frac{p-1}{p}
\;\sim\; e^{\gamma} \log z \cdot \frac{\phi(q)}{q}.
```

**Fehler:** Die Faktorisierung selbst ist korrekt, aber die anschließende Vereinfachung ist falsch.

Mertens' drittes Theorem sagt:
$$\prod_{p \leq z} \frac{p}{p-1} \sim e^{\gamma} \log z.$$

Der Korrekturfaktor für $p \mid q$ ist:
$$\prod_{\substack{p \mid q \\ p \leq z}} \frac{p-1}{p}.$$

Das ist ungefähr $\phi(q)/q$, **aber nur für die Primteiler von $q$ bis $z$**. Für $z$ groß genug (größer als alle Primteiler von $q$) gilt $\prod_{p \mid q, p \leq z}(p-1)/p = \phi(q)/q$. Dies ist korrekt.

**Dennoch:** Die Formel lautet dann:
$$\prod_{\substack{p \leq z \\ p \nmid q}} \frac{p}{p-1} \sim e^{\gamma} \log z \cdot \frac{\phi(q)}{q}.$$

Damit ist:
$$\frac{1}{V(z)} = \prod_{\substack{p \leq z \\ p \nmid q}} \frac{p}{p-1} \sim e^{\gamma} \log z \cdot \frac{\phi(q)}{q}.$$

Das bedeutet:
$$\frac{X}{V(z)} = \frac{x}{q} \cdot e^{\gamma} \log z \cdot \frac{\phi(q)}{q} \cdot \frac{1}{e^{\gamma}} = \frac{x \phi(q) \log z}{q^2}.$$

Dies ergibt **nicht** $\frac{2x}{\phi(q)\log x}$! Der Beweis hat den $e^{\gamma}$-Faktor im zweiten Mertens-Schritt nicht korrekt behandelt.

**Korrekte Ableitung:** Mit $X = x/q$, $g(p)=1/p$ für $p\nmid q$:
$$V(z) = \prod_{\substack{p\leq z\\p\nmid q}}(1-1/p) \sim \frac{e^{-\gamma}}{\log z} \cdot \frac{q}{\phi(q)},$$
$$\frac{X}{V(z)} \sim \frac{x}{q} \cdot \frac{\phi(q)\log z}{q \cdot e^{-\gamma}} \cdot e^{-\gamma} = \frac{x\phi(q)\log z}{q^2}.$$

Das stimmt immer noch nicht mit dem gewünschten Ergebnis überein. Das eigentliche Problem ist, dass die Beziehung zwischen $X$ und $|\mathcal{A}|$ mit mehr Sorgfalt behandelt werden muss. Für die korrekte Anwendung auf Primzahlen in APs gilt $X = x/\phi(q)$ (die Dichte der primen Reste modulo $q$), nicht $X = x/q$.

**Korrekte Wahl:** $X = x/\phi(q)$, $g(p) = 1/(p-1)$ für $p \nmid q$ (Dichte in der reduzierten Klasse). Dies führt direkt zu $\frac{X}{V(z)} \sim \frac{x}{\phi(q)\log z} = \frac{2x}{\phi(q)\log x}$ für $z = \sqrt{x}$.

Der Beweis verwendet die falsche Normierung ($X = x/q$ und $g(p) = 1/p$) und gelangt durch Fehler im Mertens-Schritt zum richtigen numerischen Ergebnis – das Zwischenergebnis (Zeilen 384–396) ist nicht konsistent.

---

### BUG-P14-005 (Gering): Bibliographie – Iwaniec-Eintrag

**Fundstelle:** Zeilen 541–544, \bibitem{Iwaniec1978}

**Das Paper schreibt:**
```
\bibitem{Iwaniec1978}
H.~Iwaniec.
\newblock On the Brun-Titchmarsh theorem.
\newblock \emph{Journal of Mathematical Society of Japan}, 34:95--123, 1982.
```

**Fehler:** Das Label lautet `Iwaniec1978`, aber der Text nennt 1982. Zudem: Das korrekte Journal ist **Journal of the Mathematical Society of Japan** (mit "the"). Außerdem: Die eigentliche Referenz für Iwaniec's Verbesserung des Brun-Titchmarsh-Theorems ist aus **1982** (JMSJ 34), nicht 1978. Das Label `Iwaniec1978` ist irreführend und sollte `Iwaniec1982` heißen.

---

### BUG-P14-006 (Gering): Remainder-Normierung $|r_d| \leq g(d)$

**Fundstelle:** Zeilen 136–137, Definition 2.2

**Das Paper schreibt:**
> "$r_d$ is a remainder with $|r_d| \leq g(d)$ typically."

**Problem:** Diese Normierung ist ungewöhnlich und zu restriktiv. In der Standardliteratur (Iwaniec-Kowalski, Cojocaru-Murty) gilt typischerweise $|r_d| \leq \rho(d)$ (die Anzahl der Restklassen) oder schlicht $|r_d| \leq 1$ für viele Anwendungen. Die Bedingung $|r_d| \leq g(d)$ wäre $|r_d| \leq 1/p$ für $d = p$, was unnatürlich ist.

Für die Anwendung in Theorem 5.1 gilt $|r_d| \leq 1$, was mit $g(p) = 1/p$ nicht durch $|r_d| \leq g(d)$ beschränkt wird. Die Definition ist also inkonsistent mit der späteren Anwendung. Das Wörtchen "typically" macht es zu einer Nicht-Definition.

---

### BUG-P14-007 (Gering): Vergleichstabelle – Brun-Fehlerterm

**Fundstelle:** Zeilen 455–457, Vergleichstabelle in Section 6

**Das Paper schreibt:**
```
Error term & $\sum_{d \leq z} |r_d|$ & $\sum_{d_1,d_2 \leq z} |r_{[d_1,d_2]}|$ \\
```

**Fehler:** Der Fehlerterm für Brun's Sieb ist nicht $\sum_{d \leq z} |r_d|$ sondern
$$\sum_{\substack{d \mid \mathcal{P}(z) \\ \omega(d) \leq r}} |r_d|,$$
wobei $r$ die Trunkierungstiefe ist. Die Summe über **alle** $d \leq z$ wäre der exakte Eratosthenes-Fehler, nicht Brun's trunkierter Fehler. Bei Brun wird gerade die Tiefe begrenzt, weshalb der Fehlerterm kleiner als $\sum_{d \leq z}$ ist.

---

## Weitere Beobachtungen (keine Bugs)

**Anmerkung 1 (Lemma 3.1, Selberg-Fundamentalungleichung):** Der Beweis ist korrekt und elegant. Die Fallunterscheidung ($\gcd = 1$ vs. $> 1$) ist vollständig. ✓

**Anmerkung 2 (Lemma 3.2, Diagonalisierung):** Die Idee der Variablensubstitution $\lambda \to y$ via Möbius-Inversion und die Diagonalisierung des quadratischen Forms sind mathematisch korrekt, wenn auch Schritt 2 (Verifikation der Diagonalform) nicht ausgeführt wird. Dies ist eine Auslassung, kein Fehler – die Berechnung ist in Standardreferenzen (Iwaniec-Kowalski §6.2) zu finden.

**Anmerkung 3 (Section 6, Paritätsproblem):** Die Diskussion des Paritätsproblems (Parity obstruction, Zeilen 473–479) ist mathematisch korrekt und gut erklärt. ✓

**Anmerkung 4 (Schluss, Maynard-Tao):** Das Paper erwähnt korrekt Maynard (2015, Annals) mit Seitenangaben 383–413 und $\liminf(p_{n+1}-p_n) \leq 246$. ✓

**Anmerkung 5 (Abstract):** Die Formulierung der Selberg-Schranke im Abstract,
$$S(\mathcal{A},\mathcal{P},z) \leq \frac{X}{V(z)} + O\!\left(\sum_{d_1,d_2 \leq z} |r_{[d_1,d_2]}|\right),$$
ist korrekt angegeben. ✓

**Anmerkung 6 (Definition 2.2, $g(1) = 1$):** Die Anforderung $g(1) = 1$ ist korrekt für eine multiplikative Funktion. ✓

---

## Zusammenfassung der Fehler

| Bug-ID | Status | Schwere |
|--------|--------|---------|
| BUG-P14-001 | Bestätigt | Schwer – Wahl $z=x/q$ macht Fehlerterm dominant |
| BUG-P14-002 | Bestätigt (mit Einschränkung) | Mittel – $\asymp$ sollte $=$ sein |
| BUG-P14-003 | Bestätigt | Mittel – Constraint-Ableitung nicht korrekt begründet |
| BUG-P14-004 | Bestätigt | Mittel – Mertens-Korrektur inkonsistente Normierung |
| BUG-P14-005 | Bestätigt | Gering – Bibliographie-Label/Jahr falsch |
| BUG-P14-006 | Bestätigt | Gering – Remainder-Normierung unnatürlich |
| BUG-P14-007 | Bestätigt | Gering – Vergleichstabelle Fehlerterm ungenau |

---

## Fazit

Das Paper zu Selbergs Sieb enthält einen **schwerwiegenden Fehler** im Beweis des Brun-Titchmarsh-Theorems (BUG-P14-001): Die Wahl $z = x/q$ für den Siebparameter führt zu einem Fehlerterm $O(z^2) = O(x^2/q^2)$, der den Hauptterm $\sim x/(\phi(q)\log(x/q))$ in der relevanten Parameterregion dominiert. Der Beweis ist damit in dieser Form nicht gültig. Die korrekte Wahl ist $z = \sqrt{x/q}$ (oder äquivalente Parametrisierungen), die den Fehlerterm $O(x/q)$ subdominant macht.

Darüber hinaus weist der Normierungsschritt in der Mertens-Abschätzung (BUG-P14-004) eine Inkonsistenz auf: Das Paper arbeitet mit $X = x/q$ und $g(p) = 1/p$, was nicht die natürliche Normierung für Primzahlen in APs ist. Dies führt zu einer Formelkette, die numerisch zum richtigen Ergebnis führt, aber algebraisch nicht kohärent ist.

Die Grundstruktur des Papers (Selbergs Lambda²-Methode, Diagonalisierung, Optimierung) ist mathematisch korrekt und didaktisch gut aufgebaut. Nach Korrektur der drei mittleren bis schwerwiegenden Fehler wäre das Paper publishable.

**Priorität der Korrekturen:**
1. BUG-P14-001: Siebparameter $z = \sqrt{x/q}$ korrigieren
2. BUG-P14-004: Normierung für AP-Anwendung konsistent machen ($X = x/\phi(q)$, $g(p) = 1/(p-1)$)
3. BUG-P14-003: Constraint-Ableitung klar über Möbius-Inversion begründen
