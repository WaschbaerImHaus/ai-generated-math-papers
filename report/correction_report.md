# Korrektur-Report: Ausstehende Änderungen an den Papers

**Erstellt:** 2026-03-11
**Autor des Reports:** Claude Code (Gutachter-System)
**Status:** Vom Autor eigenständig umzusetzen

---

## Übersicht

Die folgenden Fehler wurden durch das Peer-Review-System identifiziert. Die Korrekturen sind **nicht** in die LaTeX-Dateien eingepflegt worden — der Autor setzt sie eigenständig um.

---

## paper1_giuga_three_primes.tex (EN)

### BUG-001: Falscher Autorenname im Text
**Zeile:** ~86
**Ist:** `Borwein, Borwein, Borwein, and Giuga \cite{Borwein1996} established:`
**Soll:** `Borwein, Borwein, Girgensohn, and Parnes \cite{Borwein1996} established:`
**Begründung:** Die Originalveröffentlichung (Amer. Math. Monthly 103, 1996) hat die Autoren Borwein, Borwein, Girgensohn und Parnes — nicht drei Borweins und Giuga.

### BUG-002: Unklare Begründung für q ≥ 5
**Zeile:** ~247
**Ist:** `Since $q \ge 5$ (the smallest odd prime greater than~$p \ge 3$), we have`
**Soll:** `Since $q > p \ge 3$ we have $q \ge 5$, and thus`
**Begründung:** Die ursprüngliche Formulierung ist zirkulär — q ≥ 5 sollte aus q > p ≥ 3 gefolgert werden, nicht als separate Behauptung stehen.

### BUG-003: Falsche Autorenzuschreibung (Bednarek)
**Zeile:** ~295
**Ist:** `the much stronger result of Bednarek and Borwein \cite{Bednarek2014}, who`
**Soll:** `the much stronger result of Bednarek \cite{Bednarek2014}, who`
**Begründung:** Das Preprint \cite{Bednarek2014} hat nur einen Autor (M. Bednarek). Borwein ist nicht Co-Autor dieses Werks.

### BUG-003b: Unvollständige Literaturangabe
**Zeile:** ~353
**Ist:** `Preprint (2014).`
**Soll:** `Preprint, unpublished (2014).`
**Begründung:** Unveröffentlichte Preprints sollten als solche gekennzeichnet werden.

---

## paper1_giuga_drei_primfaktoren_de.tex (DE)

### BUG-001 (DE): Falscher Autorenname im Text
**Zeile:** ~86
**Ist:** `Borwein, Borwein, Borwein und Giuga \cite{Borwein1996} bewiesen:`
**Soll:** `Borwein, Borwein, Girgensohn und Parnes \cite{Borwein1996} bewiesen:`
**Begründung:** Siehe BUG-001 (EN).

### BUG-002 (DE): Unklare Begründung für q ≥ 5
**Zeile:** ~249
**Ist:** `Da $q \ge 5$ (die kleinste ungerade Primzahl größer als $p \ge 3$),`
**Soll:** `Da $q > p \ge 3$ gilt $q \ge 5$, also`
**Begründung:** Siehe BUG-002 (EN).

### BUG-003 (DE): Falsche Autorenzuschreibung (Bednarek)
**Zeile:** ~296
**Ist:** `das viel stärkere Ergebnis von Bednarek und Borwein \cite{Bednarek2014}, die zeigten,`
**Soll:** `das viel stärkere Ergebnis von Bednarek \cite{Bednarek2014}, der zeigte,`
**Begründung:** Einzelautor — grammatikalisch und inhaltlich anzupassen.

### BUG-003b (DE): Unvollständige Literaturangabe
**Zeile:** ~354
**Ist:** `Preprint (2014).`
**Soll:** `Preprint, unpublished (2014).`

---

## paper2_lehmer_three_primes.tex (EN)

### BUG-004: Uneleganter Fallunterscheidungsbeweis
**Zeile:** ~128–134
**Ist:** Beweis für den Fall p ≥ 3 über mehrere Ungleichungen mit expliziten Zahlenwerten (p=3, q≥5 und p≥5, q>p≥5)
**Soll:** Kompakter Beweis über algebraische Umformung:
```
Note that $(p-1)(q-1) - (p+q-2) = (p-2)(q-2) - 1$.
Since $p \ge 3$ we have $p - 2 \ge 1$, and since $q \ge p + 2 \ge 5$ we have $q - 2 \ge 3 \ge 2$.
Therefore $(p-2)(q-2) - 1 \ge 1 \cdot 2 - 1 = 1 > 0$.
Hence $(p-1)(q-1) > p+q-2 > 0$, so no divisibility can hold.
```
**Begründung:** Der bestehende Beweis verwendet unnötige Fallunterscheidungen; die algebraische Identität (p-2)(q-2)-1 macht ihn wesentlich klarer und kürzer.

---

## paper2_lehmer_drei_primfaktoren_de.tex (DE)

### BUG-004 (DE): Rechenfehler in Fall b=1
**Zeile:** ~150
**Ist:** `\noindent\textbf{Fall $b = 1$:} $r = 2q - 1 + 1 - 1 = 2q$.`
**Soll:** `\noindent\textbf{Fall $b = 1$:} aus $r - 1 = 2q - 1$ folgt $r = 2q$.`
**Begründung:** Die Herleitung `2q - 1 + 1 - 1 = 2q` ist rechnerisch falsch (ergibt 2q-1). Korrekte Herleitung: aus b(r-1) = 2q-1 mit b=1 folgt r-1 = 2q-1, also r = 2q.

### BUG-010: Abgebrochene Alternativformel (Entwurfsüberrest)
**Zeile:** ~144–156 (Section 3)
**Ist:** Unvollständiger Ausdruck mit `\ldots` als Platzhalter; erkennbar unfertige Alternativformulierung
**Soll:** Entweder vollständig ausformulieren oder den Abschnitt entfernen
**Begründung:** Ein `\ldots` ohne Kontext ist kein gültiger Beweisschritt und wirkt wie ein nicht entfernter Entwurfsüberrest.

### BUG-004b (DE): Unvollständiger Beweis Lemma 4.1
**Zeile:** ~263–280
**Ist:** Beweisskizze (`\begin{proof}[Beweisskizze]`) mit unvollständiger Argumentation — die Schranke $q \le \frac{3(p-1)}{p-2}$ wird nicht vollständig hergeleitet
**Soll:** Vollständiger Beweis:
```latex
\begin{proof}
Sei $a = \tfrac{p-1}{2}$, $b = \tfrac{q-1}{2}$, $c = \tfrac{r-1}{2}$.
Es gilt die Schlüsselidentität:
\[
  pqr - 1 = 8abc + 4(ab + bc + ca) + 2(a + b + c).
\]
Da $(p-1)(q-1)(r-1) = 8abc \mid (pqr - 1)$ vorausgesetzt wird, impliziert dies insbesondere
$(r-1) \mid (pq - 1)$. Sei $k = \tfrac{pq-1}{r-1}$.

Da $r > q$ und beide ungerade prim sind, gilt $r \ge q + 2$, also $r - 1 \ge q + 1 > q$.
Daher:
\[
  k = \frac{pq - 1}{r - 1} < \frac{pq}{q} = p, \quad\text{also}\quad k \le p - 1.
\]

Die Bedingung $(p-1)(q-1)(r-1) \mid (pqr-1)$ impliziert nach der Reduktion aus Schritt~2
auch $(p-1)(q-1) \mid (p + q - 1) + k$. Mit $k \le p - 1$:
\[
  (p-1)(q-1) \le (p + q - 1) + (p - 1) = 2p + q - 2.
\]
Umformen: $pq - p - q + 1 \le 2p + q - 2$, also $q(p-2) \le 3(p-1)$, d.\,h.:
\[
  q \le \frac{3(p-1)}{p-2}. \qedhere
\]
\end{proof}
```
**Begründung:** Eine Beweisskizze ist für ein Submission-Ready Paper nicht ausreichend.

---

## paper5_giuga_no_semiprime.tex (EN)

### BUG-005: Falsche Zitatreferenz (squarefree)
**Zeile:** ~57
**Ist:** `\cite{Giuga1950}`
**Soll:** `\cite{Ingwer2026sqfree}`
**Begründung:** Die Quadratfreiheit von Giuga-Pseudoprimen wird hier auf das eigene Paper \cite{Ingwer2026sqfree} gestützt, nicht auf Giuga 1950.

### BUG-006: Falsche Autorenzuschreibung (Bednarek)
**Zeile:** ~116
**Ist:** `due to Bednarek and Borwein \cite{Bednarek2014}, is $59$ prime factors.`
**Soll:** `due to Bednarek \cite{Bednarek2014}, is $59$ prime factors.`
**Begründung:** Einzelautor — siehe BUG-003.

### BUG-007: Unvollständige Literaturangabe
**Zeile:** ~151
**Ist:** `Preprint (2014).`
**Soll:** `Preprint, unpublished (2014).`

---

## paper5_giuga_kein_semiprim_de.tex (DE)

### BUG-005 (DE): Falsche Zitatreferenz (quadratfrei)
**Zeile:** ~59
**Ist:** `\cite{Giuga1950}`
**Soll:** `\cite{Ingwer2026qfrei}`
**Begründung:** Siehe BUG-005 (EN).

### BUG-006 (DE): Falsche Autorenzuschreibung (Bednarek)
**Zeile:** ~118
**Ist:** `von Bednarek und Borwein \cite{Bednarek2014}, liegt bei $59$ Primfaktoren.`
**Soll:** `von Bednarek \cite{Bednarek2014}, liegt bei $59$ Primfaktoren.`

### BUG-007 (DE): Unvollständige Literaturangabe
**Zeile:** ~155
**Ist:** `Preprint (2014).`
**Soll:** `Preprint, unpublished (2014).`

---

## paper7_lehmer_no_semiprime.tex (EN)

### BUG-008: Unstrukturierter Beweis mit falscher Zwischenrechnung
**Zeile:** ~90–106
**Problem:** Der bestehende Beweis von Ungleichung (2) enthält:
- Eine falsche Zwischenrechnung: `pq - p - q + 1 - p - q + 2 = pq - 2p - 2q + 3` — das Zwischenergebnis ist zwar korrekt, aber der Weg führt über einen selbst-widerlegenden Versuch (`Hmm — let us handle p=2 separately.`)
- Kommentar `Hmm` ist kein valider Beweistext
- Unnötige Fallunterscheidung mit sich überschneidenden Teilen

**Soll:** Saubere Fallunterscheidung:
```latex
\textbf{Case $p = 2$:} The condition $(p-1)(q-1) \mid (p+q-2)$ becomes $(q-1) \mid q$.
Since $2q - 1 = 2(q-1) + 1$, we get $(q-1) \mid 1$, so $q - 1 = 1$, i.e.\ $q = 2 = p$.
This contradicts $q > p$.

\textbf{Case $p \ge 3$:} We have $p - 2 \ge 1$ and, since $q \ge p + 2 \ge 5$, also $q - 2 \ge 3 \ge 2$.
Note that $(p-1)(q-1) - (p+q-2) = (p-2)(q-2) - 1$.  Therefore
\[
  (p-1)(q-1) - (p+q-2) = (p-2)(q-2) - 1 \ge 1 \cdot 2 - 1 = 1 > 0.
\]
Hence $(p-1)(q-1) > p + q - 2 \ge 1$, and so $(p-1)(q-1) \nmid (p+q-2)$.
```

### BUG-009: Unklare Formel in Remark (euler vs. varphi)
**Zeile:** ~123–126
**Ist:**
```latex
n - 1 = \euler(n) + \sum_{\substack{S \subseteq \{1,\ldots,k\} \\ S \ne \{1,\ldots,k\}}}
          \prod_{i \in S}(p_i - 1) \cdot \prod_{i \notin S} \mathbf{1}.
```
**Soll:**
```latex
n - 1 = \varphi(n) + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1).
```
**Begründung:** `\euler` ist kein Standard-LaTeX-Makro; `\varphi` ist die übliche Notation. Der Term `\prod_{i \notin S} \mathbf{1}` ist stets 1 und damit redundant — er verschleiert die Formel.

---

## Zusammenfassung

| Paper | BUG-ID | Typ | Priorität |
|---|---|---|---|
| paper1 EN+DE | BUG-001 | Falscher Autorenname | Hoch |
| paper1 EN+DE | BUG-002 | Unklare Begründung | Mittel |
| paper1 EN+DE | BUG-003 | Falscher Autorenname | Hoch |
| paper1 EN+DE | BUG-003b | Literaturangabe | Niedrig |
| paper2 EN | BUG-004 | Uneleganter Beweis | Mittel |
| paper2 DE | BUG-004 | Rechenfehler | Hoch |
| paper2 DE | BUG-004b | Unvollständiger Beweis | Hoch |
| paper2 DE | BUG-010 | Entwurfsüberrest | Hoch |
| paper5 EN+DE | BUG-005 | Falsche Zitatreferenz | Hoch |
| paper5 EN+DE | BUG-006 | Falscher Autorenname | Hoch |
| paper5 EN+DE | BUG-007 | Literaturangabe | Niedrig |
| paper7 EN | BUG-008 | Beweisfehler / Kommentar | Hoch |
| paper7 EN | BUG-009 | Formel / Makro | Mittel |
